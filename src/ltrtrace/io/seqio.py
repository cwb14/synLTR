"""FASTA access that never loads a genome: random access to indexed (plain or bgzip) files, a one-record-at-a-time
stream of plain, gzip or bgzip files, and reverse complements of text (every IUPAC letter, case kept) and of 2-bit
code arrays. Some modules reverse-complement on their own, over other alphabets: io.annotation.tsd (a TSD's
canonical orientation, ACGT), analyses.source (ACGTN), seq.align (its code arrays), seq.kmers (a k-mer's 2-bit value,
for its canonical form), sim.pop.library (ACGTN in either case) and sim.pop.genome (bytes, with this module's IUPAC
table)."""

from __future__ import annotations

import contextlib
import hashlib
import os
import re
from collections.abc import Iterator
from pathlib import Path

import numpy as np
import pysam

from ltrtrace.base.checkpoint import atomic_path
from ltrtrace.base.constants import FASTA_BLOCK
from ltrtrace.base.log import UsageError
from ltrtrace.io.files import gzip_kind, open_input

COMP_FROM, COMP_TO = "ACGTRYKMBVDHNacgtrykmbvdhn", "TGCAYRMKVBHDNtgcayrmkvbhdn"  # IUPAC letters and complements
_COMP = str.maketrans(COMP_FROM, COMP_TO)
_N_RUN = re.compile(r"[Nn]+")
_SPACE = (b" ", b"\t", b"\r", b"\x0b", b"\x0c")  # whitespace bytes.strip removes, besides the line end


def revcomp(seq: str) -> str:
    """Reverse complement (IUPAC-aware, case preserved)."""
    return seq.translate(_COMP)[::-1]


def revcomp_codes(codes: np.ndarray) -> np.ndarray:
    """Reverse complement of a uint8 code array (ACGT = 0-3 complement to 3-c; 4 stays 4)."""
    return np.where(codes < 4, 3 - codes, codes).astype(np.uint8)[::-1].copy()


def n_runs(seq: str) -> list[tuple[int, int]]:
    """0-based half-open runs of N/n in `seq`."""
    return [m.span() for m in _N_RUN.finditer(seq)]


_CODES = np.full(256, 4, dtype=np.uint8)
_CODES[np.frombuffer(b"ACGTacgt", dtype=np.uint8)] = [0, 1, 2, 3, 0, 1, 2, 3]
_CODE_TABLE = _CODES.tobytes()  # byte -> code, for bytes.translate


def encode(seq: bytes) -> np.ndarray:
    """Sequence bytes as uint8 codes: ACGT/acgt = 0-3, anything else = 4 (a read-only array)."""
    return np.frombuffer(bytes(seq).translate(_CODE_TABLE), dtype=np.uint8)


def read_fai(path: str | Path) -> dict[str, int]:
    """Sequence name -> length (bp) from a FASTA index (.fai: name, length, then samtools' offsets, tab-separated), in
    file order. Blank lines are skipped; any other line without a whole-number length is a ValueError naming the file
    and the line."""
    out: dict[str, int] = {}
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            if not line.strip():
                continue
            f = line.rstrip("\r\n").split("\t")
            if len(f) < 2 or not (f[1].isascii() and f[1].isdigit()):
                raise ValueError(f"{path}: line {n} is not a FASTA index line (name, length, ...)")
            out[f[0]] = int(f[1])
    return out


def _build_index(path: Path, fai: Path, gzi: Path | None) -> None:
    """Index `path` into `fai` (and `gzi` for bgzip), under temporary names renamed into place, the .gzi first:
    a reader that finds the .fai finds a complete index."""
    with contextlib.ExitStack() as stack:   # leaves in reverse: the .gzi is renamed first
        tmp_fai = stack.enter_context(atomic_path(fai))
        tmp_gzi = None if gzi is None else stack.enter_context(atomic_path(gzi))
        pysam.faidx(str(path), "--fai-idx", str(tmp_fai), *(["--gzi-idx", str(tmp_gzi)] if tmp_gzi else []))


class Fasta:
    """Indexed FASTA. Coordinates are 0-based half-open; output is uppercase.

    The index (.fai, plus .gzi for bgzip) is used if present next to the FASTA, built there if the directory is
    writable, else built in `cache_dir` under a name keyed by the FASTA's path, size and mtime (a changed FASTA
    gets a new index, never a stale one). An index is built under temporary names and renamed, so concurrent
    runs never read a partial one. Plain gzip is rejected up front (htslib reads only bgzip); an index pysam
    cannot build is a ValueError naming the file.
    """

    def __init__(self, path: str | Path, cache_dir: Path | None = None):
        self.path = Path(path)
        self.cache_dir = None if cache_dir is None else Path(cache_dir)  # kept: workers reopen alike
        if not self.path.is_file():
            raise FileNotFoundError(self.path)
        self._fh = self._open(cache_dir)
        self._lengths = dict(zip(self._fh.references, self._fh.lengths))

    def _open(self, cache_dir: Path | None) -> pysam.FastaFile:
        path = self.path
        kind = gzip_kind(path)
        if kind == "gzip":
            raise ValueError(f"{path} is gzip-compressed but not bgzip; recompress with bgzip")
        gz = kind is not None
        fai, gzi = Path(f"{path}.fai"), Path(f"{path}.gzi")
        try:
            if fai.exists() and (gzi.exists() or not gz):
                return pysam.FastaFile(str(path))
            if os.access(path.parent, os.W_OK):
                _build_index(path, fai, gzi if gz else None)
                return pysam.FastaFile(str(path))
            if cache_dir is None:
                raise ValueError(f"no index for {path} and its directory is not writable; the index cannot be "
                                 "written, pass a cache directory")
            cache_dir = Path(cache_dir)
            try:
                cache_dir.mkdir(parents=True, exist_ok=True)
            except OSError as exc:  # the cache (an outdir's work/fai), not the FASTA, is at fault: say so
                raise OSError(f"cannot create the FASTA index cache {cache_dir}: {exc}") from exc
            st = path.stat()
            tag = hashlib.sha1(f"{path.resolve()}\0{st.st_size}\0{st.st_mtime_ns}".encode()).hexdigest()[:12]
            fai = cache_dir / f"{tag}.{path.name}.fai"
            gzi = cache_dir / f"{tag}.{path.name}.gzi"
            if not fai.exists() or (gz and not gzi.exists()):
                _build_index(path, fai, gzi if gz else None)
            return pysam.FastaFile(str(path), filepath_index=str(fai),
                                   filepath_index_compressed=str(gzi) if gz else None)
        except pysam.SamtoolsError as exc:  # not an OSError or ValueError: callers would not catch it
            raise ValueError(f"{path}: cannot index the FASTA: {exc}") from exc

    def lengths(self) -> dict[str, int]:
        """Sequence lengths in index (file) order."""
        return dict(self._lengths)

    def length(self, seqid: str) -> int:
        """The length of sequence `seqid` (KeyError if the FASTA lacks it)."""
        return self._lengths[seqid]

    def fetch(self, seqid: str, start: int, end: int) -> str:
        """Uppercase sequence [start, end), clipped to the sequence."""
        if seqid not in self._lengths:
            raise KeyError(seqid)
        start = max(start, 0)
        end = min(end, self._lengths[seqid])
        if start >= end:
            return ""
        return self._fh.fetch(seqid, start, end).upper()

    def close(self) -> None:
        self._fh.close()

    def __enter__(self) -> Fasta:  # noqa: PYI034 (typing.Self needs py3.11)
        return self

    def __exit__(self, *exc) -> None:
        self.close()


def open_fasta(path: str | Path, label: str, cache_dir: Path | None = None) -> Fasta:
    """Genome `label`'s indexed FASTA (`Fasta`, its index cached in `cache_dir` when its folder is read-only); a FASTA
    pysam cannot open or index is one UsageError line naming the genome and the file."""
    try:
        return Fasta(path, cache_dir=cache_dir)
    except Exception as exc:  # pysam raises several kinds on a malformed or unindexable FASTA
        raise UsageError(f"{label}: cannot open FASTA {path}: {exc}") from exc


def _record(text: bytes) -> tuple[str, bytes]:
    """(name, sequence) of one record's text, its header line (without '>') and its sequence lines: the name is the
    header's first word; each sequence line is stripped of surrounding whitespace and the lines joined."""
    head, _, body = text.partition(b"\n")
    fields = head.split()
    name = fields[0].decode() if fields else ""
    if any(c in body for c in _SPACE):
        return name, b"".join(line.strip() for line in body.split(b"\n"))
    return name, body.replace(b"\n", b"")


def stream_fasta(path: str | Path) -> Iterator[tuple[str, bytes]]:
    """Yield (name, sequence) of each record of a plain, gzip or bgzip FASTA, in order.

    The name is the first word of the header line; each sequence line is stripped of surrounding whitespace. The file
    is read in blocks of FASTA_BLOCK bytes and split into whole records, so memory holds a block and one unfinished
    record. The compression is sniffed from the file's magic bytes, so a library FASTA in plain gzip streams as well as
    a bgzip one (only indexed access needs bgzip). Non-blank text before the first header is a ValueError.
    """
    pending, started = b"", False

    def fail() -> ValueError:
        return ValueError(f"{path}: not a FASTA file (sequence or text before the first '>' header)")

    with open_input(path, "rb") as fh:
        while block := fh.read(FASTA_BLOCK):
            text = pending + block
            if not started:
                i = 0 if text.startswith(b">") else text.find(b"\n>") + 1
                if i == 0 and not text.startswith(b">"):
                    if text.strip():
                        raise fail()
                    pending = text
                    continue
                if text[:i].strip():
                    raise fail()
                text, started = text[i:], True
            last = text.rfind(b"\n>")
            if last < 0:
                pending = text
                continue
            yield from map(_record, text[1:last].split(b"\n>"))
            pending = text[last + 1:]
    if not started:
        if pending.strip():
            raise fail()
        return
    yield from map(_record, pending[1:].split(b"\n>"))
