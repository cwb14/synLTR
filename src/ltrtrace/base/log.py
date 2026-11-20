"""Logging and the exceptions that map to exit codes."""

import logging
import sys

log = logging.getLogger("ltrtrace")


class LtrtraceError(Exception):
    """A run failed (exit code 1)."""

    exit_code = 1


class UsageError(LtrtraceError):
    """The command line or an input file is unusable (exit code 2)."""

    exit_code = 2


class AnalysisError(LtrtraceError):
    """The analyses failed after genotyping finished (exit code 3)."""

    exit_code = 3


class _Lines(logging.Formatter):
    """A log line: the message, after "warning: " for a warning (as the command line's error line follows
    "error: ")."""

    def format(self, record: logging.LogRecord) -> str:
        text = super().format(record)
        return f"warning: {text}" if record.levelno == logging.WARNING else text


def setup(verbose: bool) -> None:
    """Send log messages to stderr: milestones by default, progress and checks with `verbose`; a warning's line
    starts with "warning: "."""
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(_Lines("%(message)s"))
    log.handlers[:] = [handler]
    log.setLevel(logging.DEBUG if verbose else logging.INFO)
    log.propagate = False
