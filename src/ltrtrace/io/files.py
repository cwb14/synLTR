"""Input files read whatever their compression: plain, gzip or bgzip told apart by the file's first bytes, never by its
name, so a gzip file without .gz in its name reads as gzip and a plain file as text (indexed genome access, which
needs bgzip, checks the kind itself: io.seqio.Fasta)."""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import IO

GZIP_MAGIC = b"\x1f\x8b"  # the first two bytes of every gzip member (RFC 1952)
_BGZF_EXTRA = slice(12, 14)  # a BGZF block's first extra subfield identifier, "BC" (SAM/BAM specification)


def gzip_kind(path: str | Path) -> str | None:
    """None (not gzip), "bgzip" (BGZF's 'BC' extra field) or "gzip" (plain gzip), from the file's first bytes."""
    with open(path, "rb") as fh:
        head = fh.read(14)
    if head[:2] != GZIP_MAGIC:
        return None
    return "bgzip" if head[_BGZF_EXTRA] == b"BC" else "gzip"


def open_input(path: str | Path, mode: str = "rt", **kwargs) -> IO:
    """A handle reading the file at `path` ("rt" text or "rb" bytes; `kwargs` such as encoding and newline go to the
    opener), decompressed when its first bytes say gzip or bgzip."""
    if mode not in ("rt", "rb"):
        raise ValueError(f"open_input reads: mode 'rt' or 'rb', not {mode!r}")
    return (gzip.open if gzip_kind(path) else open)(path, mode, **kwargs)
