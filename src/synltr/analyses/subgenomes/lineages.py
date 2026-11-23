"""k-mer lineages of the run: the library's clades (seq.kmers: k-mers carried by 2..CAP copies on >= 2 sequences,
FracMinHash-sampled above KMER_BUDGET bases), a checkpoint stage in the output folder (work/subgenomes/clades.npz)."""
from __future__ import annotations

from pathlib import Path

from synltr.analyses.subgenomes.inputs import log
from synltr.base.constants import CAP, KMER_BUDGET
from synltr.seq.kmers import auto_scale, clades

CHECKPOINTS = Path("work") / "subgenomes"  # under the output folder: the clade checkpoint


def load_or_build_clades(a, lib):
    """The clades of the run's library (a.threads numba threads), loaded from the output folder's checkpoint when its
    key holds (same library, parameters and code), else built and checkpointed there."""
    scale = auto_scale(int(lib.codes.size))
    cl = clades(lib.store, lib.seq, scale=scale, threads=a.threads, cache=Path(a.outdir) / CHECKPOINTS)
    log(f"{cl.n_occ:,} k-mers" + (f" (1/{scale} sampled: library > {KMER_BUDGET / 1e6:.0f} Mb)" if scale > 1 else
                                   " (exact)") + f"; {cl.n:,} clades of 2-{CAP} copies on >= 2 sequences")
    return cl
