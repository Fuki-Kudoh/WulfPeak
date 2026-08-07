"""Artifact validators used both during execution and for manifest checks."""

from __future__ import annotations

import gzip
import json
import shutil
import subprocess
import zipfile
from pathlib import Path

from .atomic import atomic_write_json
from .models import ValidationResult


class OutputValidationError(ValueError):
    """Raised when validation metadata is missing or malformed."""


BIGWIG_MAGIC = {b"\x26\xfc\x8f\x88", b"\x88\x8f\xfc\x26"}


def _regular_file_checks(path: Path) -> tuple[dict[str, object], int | None]:
    checks: dict[str, object] = {
        "exists": path.exists(),
        "regular_file": path.is_file(),
    }
    if not path.is_file():
        return checks, None
    size = path.stat().st_size
    checks["size"] = size
    return checks, size


def _validate_gzip_fastq(path: Path, checks: dict[str, object]) -> bool:
    records = 0
    try:
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            while True:
                header = handle.readline()
                if header == "":
                    break
                sequence = handle.readline()
                plus = handle.readline()
                quality = handle.readline()
                if not sequence or not plus or not quality:
                    checks["error"] = "truncated FASTQ record"
                    return False
                if not header.startswith("@") or not plus.startswith("+"):
                    checks["error"] = "invalid FASTQ record markers"
                    return False
                if len(sequence.rstrip("\r\n")) != len(quality.rstrip("\r\n")):
                    checks["error"] = "sequence and quality lengths differ"
                    return False
                records += 1
    except (OSError, EOFError, UnicodeError) as exc:
        checks["error"] = str(exc)
        return False
    checks["records"] = records
    return records > 0


def _validate_fastqc_html(path: Path, checks: dict[str, object]) -> bool:
    try:
        with path.open("rb") as handle:
            prefix = handle.read(4096).lower()
    except OSError as exc:
        checks["error"] = str(exc)
        return False
    html = b"<html" in prefix or b"<!doctype html" in prefix
    checks["html_signature"] = html
    return html


def _validate_fastqc_zip(path: Path, checks: dict[str, object]) -> bool:
    try:
        if not zipfile.is_zipfile(path):
            checks["zip_file"] = False
            return False
        with zipfile.ZipFile(path) as archive:
            names = archive.namelist()
            bad_member = archive.testzip()
    except (OSError, zipfile.BadZipFile) as exc:
        checks["error"] = str(exc)
        return False
    checks.update(
        {
            "zip_file": True,
            "members": len(names),
            "bad_member": bad_member,
        }
    )
    return bool(names) and bad_member is None


def _samtools_check(
    path: Path,
    checks: dict[str, object],
    *,
    samtools: str | Path | None,
    index: bool,
) -> bool:
    executable = str(samtools) if samtools is not None else shutil.which("samtools")
    if executable is None:
        checks["samtools"] = None
        checks["error"] = "samtools is required for BAM validation"
        return False
    argv = (
        [executable, "idxstats", str(Path(str(path)[: -len(".bai")]))]
        if index
        else [executable, "quickcheck", str(path)]
    )
    try:
        completed = subprocess.run(argv, check=False, capture_output=True, text=True)
    except OSError as exc:
        checks["error"] = str(exc)
        return False
    checks.update(
        {
            "samtools": executable,
            "samtools_command": argv[1],
            "samtools_exit_code": completed.returncode,
        }
    )
    if completed.stderr.strip():
        checks["samtools_stderr"] = completed.stderr.strip()
    return completed.returncode == 0


def validate_path(
    path: Path,
    kind: str,
    *,
    validator: str | None = None,
    samtools: str | Path | None = None,
) -> ValidationResult:
    selected = validator or kind
    if selected == "multiqc_output":
        checks: dict[str, object] = {
            "exists": path.exists(),
            "directory": path.is_dir(),
        }
        if not path.is_dir():
            return ValidationResult(False, checks)
        report = path / "multiqc_report.html"
        data = path / "multiqc_data"
        sources = data / "multiqc_sources.txt"
        sources_nonempty = False
        if sources.is_file() and sources.stat().st_size > 0:
            try:
                sources_nonempty = bool(sources.read_text(encoding="utf-8").strip())
            except (OSError, UnicodeError) as exc:
                checks["sources_error"] = str(exc)
        checks.update(
            {
                "report_exists": report.is_file(),
                "report_size": report.stat().st_size if report.is_file() else 0,
                "data_exists": data.is_dir(),
                "sources_exists": sources.is_file(),
                "sources_size": sources.stat().st_size if sources.is_file() else 0,
                "sources_nonempty": sources_nonempty,
            }
        )
        return ValidationResult(
            bool(checks["report_size"])
            and bool(checks["data_exists"])
            and sources_nonempty,
            checks,
        )
    if selected == "trim_galore_reports":
        checks = {"exists": path.exists(), "directory": path.is_dir()}
        if not path.is_dir():
            return ValidationResult(False, checks)
        combined = path / "trimming_report.txt"
        text_reports = sorted(path.glob("*_trimming_report.txt"))
        json_reports = sorted(path.glob("*_trimming_report.json"))
        text_nonempty = bool(text_reports) and all(
            report.stat().st_size > 0 for report in text_reports
        )
        json_valid = True
        for report in json_reports:
            try:
                json.loads(report.read_text(encoding="utf-8"))
            except (json.JSONDecodeError, OSError, UnicodeError) as exc:
                checks["json_error"] = f"{report.name}: {exc}"
                json_valid = False
                break
        checks.update(
            {
                "combined_exists": combined.is_file(),
                "combined_size": combined.stat().st_size if combined.is_file() else 0,
                "text_report_count": len(text_reports),
                "text_reports_nonempty": text_nonempty,
                "json_report_count": len(json_reports),
                "json_reports_valid": json_valid,
            }
        )
        return ValidationResult(
            bool(checks["combined_size"]) and text_nonempty and json_valid,
            checks,
        )
    if kind == "directory":
        checks = {"exists": path.exists(), "directory": path.is_dir()}
        return ValidationResult(path.is_dir(), checks)

    checks, size = _regular_file_checks(path)
    if size is None:
        return ValidationResult(False, checks)

    warnings: list[str] = []
    if selected in {
        "narrowPeak_or_empty",
        "broadPeak_or_empty",
        "bed_or_empty",
        "bed3_or_empty",
        "gappedPeak_or_empty",
    } or kind in {"narrowPeak", "broadPeak", "bed3"}:
        if size == 0:
            warnings.append("ZERO_PEAKS")
        return ValidationResult(True, checks, tuple(warnings))
    if size == 0:
        return ValidationResult(False, checks)
    if selected in {"gzip_fastq", "fastq_nonempty", "fastq"}:
        return ValidationResult(_validate_gzip_fastq(path, checks), checks)
    if selected in {"fastqc_html", "html_nonempty", "html"}:
        return ValidationResult(_validate_fastqc_html(path, checks), checks)
    if selected in {"fastqc_zip", "zip_nonempty", "zip"}:
        return ValidationResult(_validate_fastqc_zip(path, checks), checks)
    if selected in {"final_bam", "pooled_bam", "bam"}:
        return ValidationResult(
            _samtools_check(path, checks, samtools=samtools, index=False), checks
        )
    if selected in {"bam_index", "bai"}:
        try:
            with path.open("rb") as handle:
                magic = handle.read(4)
        except OSError as exc:
            checks["error"] = str(exc)
            return ValidationResult(False, checks)
        checks["magic"] = magic.hex()
        if magic != b"BAI\x01":
            return ValidationResult(False, checks)
        return ValidationResult(
            _samtools_check(path, checks, samtools=samtools, index=True), checks
        )
    if selected in {"bigwig_magic", "bigwig"} or kind == "bigwig":
        with path.open("rb") as handle:
            magic = handle.read(4)
        checks["magic"] = magic.hex()
        return ValidationResult(magic in BIGWIG_MAGIC, checks)
    if selected in {"text_nonempty", "summary_tsv", "warnings_tsv", "text"}:
        try:
            non_whitespace = bool(path.read_text(encoding="utf-8").strip())
        except (OSError, UnicodeError) as exc:
            checks["error"] = str(exc)
            return ValidationResult(False, checks)
        checks["non_whitespace"] = non_whitespace
        return ValidationResult(non_whitespace, checks)
    return ValidationResult(True, checks)


def validate_artifact(
    artifact: dict[str, object],
    *,
    temporary: bool,
    samtools: str | Path | None = None,
) -> ValidationResult:
    path_key = "temporary_path" if temporary else "canonical_path"
    return validate_path(
        Path(str(artifact[path_key])),
        str(artifact.get("kind", "file")),
        validator=str(artifact.get("validator", "file")),
        samtools=samtools,
    )


def validate_output_manifest(output_dir: str | Path) -> tuple[Path, bool, dict[str, object]]:
    output = Path(output_dir).expanduser().resolve()
    config_manifest = output / "config" / "manifest.json"
    output_manifest = output / "output_manifest.json"
    if not config_manifest.is_file():
        raise OutputValidationError(f"Missing config manifest: {config_manifest}")
    if not output_manifest.is_file():
        raise OutputValidationError(f"Missing output manifest: {output_manifest}")
    try:
        json.loads(config_manifest.read_text(encoding="utf-8"))
        payload = json.loads(output_manifest.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as exc:
        raise OutputValidationError(f"Cannot read validation metadata: {exc}") from exc
    entries = payload.get("outputs")
    if not isinstance(entries, list):
        raise OutputValidationError("output_manifest.json must contain an outputs list")
    all_valid = True
    for entry in entries:
        if not isinstance(entry, dict) or "path" not in entry:
            raise OutputValidationError("output manifest contains a malformed entry")
        if not entry.get("expected", True):
            entry.update({"exists": False, "validated": True, "checks": {}, "warnings": []})
            continue
        result = validate_path(
            Path(str(entry["path"])),
            str(entry.get("kind", "file")),
            validator=str(entry.get("validator", entry.get("kind", "file"))),
        )
        entry.update(
            {
                "exists": bool(result.checks.get("exists")),
                "validated": result.valid,
                "checks": result.checks,
                "warnings": list(result.warnings),
            }
        )
        if entry.get("required", False) and not result.valid:
            all_valid = False
    payload["validated"] = all_valid
    atomic_write_json(output_manifest, payload)
    return output_manifest, all_valid, payload
