"""B's minimap2 index for the flank-mm2 route: built once per B FASTA, then loaded by mappy.

<workdir>/anchors/<B>.asm20.mmi, built by the minimap2 binary when it is on PATH, else by mappy (both write the
same single-part index), keyed by the FASTA's checksum, the preset, the minimap2 and mappy versions and the code of
this module only, so an edit elsewhere in the route never rebuilds an index (minutes and gigabytes for a large
genome).
"""
from __future__ import annotations

import logging
import shutil
import subprocess
import time
from pathlib import Path
from typing import TYPE_CHECKING

import mappy

from ltrtrace.base.checkpoint import (
    atomic_path,
    clear_done,
    code_fingerprint,
    file_checksum,
    is_done,
    mark_done,
    marker,
    stage_key,
    tool_version,
)
from ltrtrace.base.constants import ANCHOR_BEST_N, ANCHOR_PRESET, CHECKSUM_CACHE, MM_F_EQX

if TYPE_CHECKING:
    from ltrtrace.io.seqio import Fasta

log = logging.getLogger(__name__)


def index_versions() -> dict[str, str | None]:
    """The versions of what builds and reads B's index: the minimap2 binary (when on PATH) and mappy. Part of
    the index's key and of the flank-mm2 route's params."""
    return {"minimap2": tool_version("minimap2"), "mappy": mappy.__version__}


def load_index(label: str, fa: Fasta, workdir: Path, threads: int) -> mappy.Aligner:
    """mappy aligner on genome `label`'s asm20 index (its FASTA `fa`), built once at
    <workdir>/anchors/<label>.asm20.mmi with `threads`."""
    fasta = fa.path
    names = list(fa.lengths())
    mmi = Path(workdir) / "anchors" / f"{label}.{ANCHOR_PRESET}.mmi"
    done = marker(mmi)
    sha = file_checksum(fasta, Path(workdir) / CHECKSUM_CACHE)
    key = stage_key("mm2_index", {"B": sha}, {"preset": ANCHOR_PRESET, **index_versions()},
                    code_fingerprint("routes.assembly.index"))
    if is_done(done, key):
        log.debug("%s: reusing index %s", label, mmi)
        aligner = _load_index(mmi)
    else:
        clear_done(done)
        mmi.parent.mkdir(parents=True, exist_ok=True)
        t0 = time.perf_counter()
        aligner = None
        with atomic_path(mmi) as tmp:
            if shutil.which("minimap2"):
                _minimap2_index(fasta, tmp, sum(fa.lengths().values()), max(1, threads))
            else:
                log.info("%s: indexing %s with mappy (no minimap2 on PATH)", label, fasta)
                # a failed build raises inside the block, so no partial index is left behind
                aligner = _usable(mappy.Aligner(fn_idx_in=str(fasta), preset=ANCHOR_PRESET, best_n=ANCHOR_BEST_N,
                                                extra_flags=MM_F_EQX, n_threads=max(1, threads),
                                                fn_idx_out=str(tmp)), fasta, mmi)
        if aligner is None:
            aligner = _load_index(mmi)
        if aligner:
            mark_done(done, key, [mmi])
            log.debug("%s: indexed in %.1f s", label, time.perf_counter() - t0)
    aligner = _usable(aligner, fasta, mmi)
    if list(aligner.seq_names) != names:
        clear_done(done)  # rebuilt on the next run
        raise RuntimeError(f"index {mmi} does not hold the sequences of {fasta} (in order)")
    return aligner


def _usable(aligner: mappy.Aligner | None, fasta: Path, mmi: Path) -> mappy.Aligner:
    """`aligner`, or the one error of an index that could not be built or loaded."""
    if not aligner:
        raise RuntimeError(f"could not build or load the minimap2 index of {fasta} ({mmi})")
    return aligner


def _load_index(mmi: Path) -> mappy.Aligner:
    """The mapping aligner on index `mmi`: asm20, ANCHOR_BEST_N secondaries, =/X CIGARs."""
    return mappy.Aligner(fn_idx_in=str(mmi), preset=ANCHOR_PRESET, best_n=ANCHOR_BEST_N, extra_flags=MM_F_EQX)


def _minimap2_index(fasta: Path, out: Path, total_bp: int, threads: int) -> None:
    """Index `fasta` with minimap2 -x asm20 -d.

    -I above the genome size keeps the index in one part: mappy loads only the first part of an index.
    """
    cmd = ["minimap2", "-x", ANCHOR_PRESET, "-t", str(threads), "-I", str(max(total_bp, 1)),
           "-d", str(out), str(fasta)]
    log.debug("indexing: %s", " ".join(cmd))
    res = subprocess.run(cmd, check=False, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"minimap2 could not index {fasta} (exit {res.returncode}): "
                           f"{res.stderr.strip()[-2000:]}")
    log.debug("minimap2: %s", res.stderr.strip())
