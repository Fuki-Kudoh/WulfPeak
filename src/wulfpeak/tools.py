"""External executable discovery and best-effort version collection."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from .models import ToolInfo


REQUIRED_TOOLS = (
    "fastqc",
    "trim_galore",
    "bowtie2",
    "samtools",
    "bamCoverage",
    "macs3",
    "bedtools",
    "multiqc",
)


class ToolDetectionError(ValueError):
    """Raised when a required executable is not available."""


def detect_tools(
    required: tuple[str, ...] = REQUIRED_TOOLS,
) -> tuple[dict[str, ToolInfo], list[dict[str, str]]]:
    missing: list[str] = []
    tools: dict[str, ToolInfo] = {}
    warnings: list[dict[str, str]] = []
    for name in required:
        resolved = shutil.which(name)
        if resolved is None:
            missing.append(name)
            continue
        path = Path(resolved).resolve()
        version_command = (str(path), "--version")
        version: str | None = None
        try:
            completed = subprocess.run(
                list(version_command),
                check=False,
                capture_output=True,
                text=True,
                timeout=10,
            )
            combined = "\n".join(
                part.strip() for part in (completed.stdout, completed.stderr) if part.strip()
            )
            if combined:
                version = combined.splitlines()[0]
        except (OSError, subprocess.SubprocessError):
            version = None
        if version is None:
            warnings.append(
                {
                    "code": "TOOL_VERSION_UNKNOWN",
                    "message": f"Could not determine version for {name}",
                }
            )
        tools[name] = ToolInfo(name, path, version, version_command)
    if missing:
        raise ToolDetectionError(
            "Required external tools were not found on PATH: " + ", ".join(missing)
        )
    return tools, warnings
