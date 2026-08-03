"""Shell-free command execution with explicit logging and return-code checks."""

from __future__ import annotations

import os
import shlex
import subprocess
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


class CommandExecutionError(RuntimeError):
    """Raised when an external command exits unsuccessfully."""

    def __init__(self, message: str, *, returncode: int | None = None):
        super().__init__(message)
        self.returncode = returncode


@dataclass(frozen=True)
class CommandResult:
    argv: tuple[str, ...]
    returncode: int
    started_at: str
    finished_at: str
    log_path: Path
    stdout_path: Path | None = None


def display_command(argv: list[str] | tuple[str, ...]) -> str:
    return shlex.join(str(argument) for argument in argv)


class CommandRunner:
    def run(
        self,
        argv: list[str],
        log_path: str | Path,
        *,
        stdout_path: str | Path | None = None,
    ) -> CommandResult:
        if not argv:
            raise ValueError("command argv must not be empty")
        log = Path(log_path)
        log.parent.mkdir(parents=True, exist_ok=True)
        started = datetime.now(timezone.utc).isoformat()
        destination = Path(stdout_path) if stdout_path is not None else None
        temporary: Path | None = None
        output_handle = None
        try:
            if destination is not None:
                destination.parent.mkdir(parents=True, exist_ok=True)
                descriptor, temporary_name = tempfile.mkstemp(
                    dir=destination.parent,
                    prefix=f".{destination.name}.",
                    suffix=".tmp",
                )
                temporary = Path(temporary_name)
                output_handle = os.fdopen(descriptor, "wb")
            with log.open("a", encoding="utf-8") as handle:
                handle.write(
                    f"started_at: {started}\ncommand: {display_command(argv)}\n"
                )
                if destination is not None:
                    handle.write(f"stdout_path: {destination}\n")
                handle.flush()
                completed = subprocess.run(
                    argv,
                    check=False,
                    stdout=output_handle if output_handle is not None else handle,
                    stderr=handle,
                )
                if output_handle is not None:
                    output_handle.flush()
                    os.fsync(output_handle.fileno())
                    output_handle.close()
                    output_handle = None
                finished = datetime.now(timezone.utc).isoformat()
                handle.write(
                    f"finished_at: {finished}\nexit_code: {completed.returncode}\n"
                )
            if completed.returncode == 0 and temporary is not None and destination is not None:
                os.replace(temporary, destination)
                temporary = None
        finally:
            if output_handle is not None:
                output_handle.close()
            if temporary is not None:
                temporary.unlink(missing_ok=True)
        result = CommandResult(
            tuple(argv), completed.returncode, started, finished, log, destination
        )
        if completed.returncode != 0:
            raise CommandExecutionError(
                f"Command failed with exit code {completed.returncode}; see {log}",
                returncode=completed.returncode,
            )
        return result

    def run_spec(
        self, command: dict[str, object], log_path: str | Path
    ) -> CommandResult:
        """Execute one structured command-plan action."""

        if command.get("type") != "command":
            raise ValueError("command spec must use type='command'")
        argv = command.get("argv")
        if not isinstance(argv, list) or not all(isinstance(item, str) for item in argv):
            raise ValueError("command spec argv must be a list of strings")
        stdout_path = command.get("stdout_path")
        if stdout_path is not None and not isinstance(stdout_path, str):
            raise ValueError("command spec stdout_path must be a string or null")
        return self.run(argv, log_path, stdout_path=stdout_path)

    def run_pipeline(
        self, commands: list[list[str]], log_path: str | Path
    ) -> list[int]:
        """Run a stdout pipeline and require every process to succeed."""

        if not commands or any(not command for command in commands):
            raise ValueError("pipeline commands must not be empty")
        log = Path(log_path)
        log.parent.mkdir(parents=True, exist_ok=True)
        processes: list[subprocess.Popen[bytes]] = []
        previous_stdout = None
        with log.open("ab") as handle:
            handle.write(
                (
                    "pipeline_started_at: "
                    + datetime.now(timezone.utc).isoformat()
                    + "\ncommands:\n"
                    + "\n".join(f"  - {display_command(command)}" for command in commands)
                    + "\n"
                ).encode("utf-8")
            )
            try:
                for index, command in enumerate(commands):
                    process = subprocess.Popen(
                        command,
                        stdin=previous_stdout,
                        stdout=subprocess.PIPE if index < len(commands) - 1 else handle,
                        stderr=handle,
                    )
                    if previous_stdout is not None:
                        previous_stdout.close()
                    previous_stdout = process.stdout
                    processes.append(process)
                returncodes = [process.wait() for process in processes]
            except BaseException:
                for process in processes:
                    if process.poll() is None:
                        process.terminate()
                for process in processes:
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        process.kill()
                raise
            handle.write(
                (
                    "pipeline_finished_at: "
                    + datetime.now(timezone.utc).isoformat()
                    + f"\nexit_codes: {returncodes}\n"
                ).encode("utf-8")
            )
        if any(code != 0 for code in returncodes):
            first_failure = next(code for code in returncodes if code != 0)
            raise CommandExecutionError(
                f"Pipeline failed with exit codes {returncodes}; see {log}",
                returncode=first_failure,
            )
        return returncodes
