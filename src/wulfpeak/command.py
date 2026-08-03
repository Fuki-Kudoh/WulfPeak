"""Shell-free command execution with explicit logging and return-code checks."""

from __future__ import annotations

import shlex
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import TextIO


class CommandExecutionError(RuntimeError):
    """Raised when an external command exits unsuccessfully."""


@dataclass(frozen=True)
class CommandResult:
    argv: tuple[str, ...]
    returncode: int
    started_at: str
    finished_at: str
    log_path: Path


def display_command(argv: list[str] | tuple[str, ...]) -> str:
    return shlex.join(str(argument) for argument in argv)


class CommandRunner:
    def run(
        self,
        argv: list[str],
        log_path: str | Path,
        *,
        stdout: TextIO | int | None = None,
    ) -> CommandResult:
        if not argv:
            raise ValueError("command argv must not be empty")
        log = Path(log_path)
        log.parent.mkdir(parents=True, exist_ok=True)
        started = datetime.now(timezone.utc).isoformat()
        with log.open("a", encoding="utf-8") as handle:
            handle.write(f"started_at: {started}\ncommand: {display_command(argv)}\n")
            handle.flush()
            completed = subprocess.run(
                argv,
                check=False,
                stdout=stdout if stdout is not None else handle,
                stderr=handle,
            )
            finished = datetime.now(timezone.utc).isoformat()
            handle.write(f"finished_at: {finished}\nexit_code: {completed.returncode}\n")
        result = CommandResult(tuple(argv), completed.returncode, started, finished, log)
        if completed.returncode != 0:
            raise CommandExecutionError(
                f"Command failed with exit code {completed.returncode}; see {log}"
            )
        return result

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
            raise CommandExecutionError(
                f"Pipeline failed with exit codes {returncodes}; see {log}"
            )
        return returncodes
