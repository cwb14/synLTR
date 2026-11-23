"""The library route: evidence for gene flow from genomes' LTR-RT libraries alone (the screen), pair by pair.

For each pair of genomes: candidate record pairs from shared low-copy k-mers (seq.kmers.shared_kmers: one scan over
every library gives every pair's, with every trio's single-copy tallies); reciprocal best
pairs as anchors when they pass the LTR clock or share their TSD (clock); chromosome pairs enriched for anchors and
their collinear blocks (correspondence); each record's syntenic window in the other genome (absence); 1:1 orthologs
inside the windows (orthology); then, per pair, the call files of both directions, the ortholog divergences and the
blocks (run).

Route properties (kept as they were validated; internal, never options):
- record sequence: N-stripped (every base but A, C, G and T removed, so LTRquest's two clean conventions, nests cut
  out and nests N-masked, give one sequence), with the record-end LTR rule on every record (io.library's screen
  rule) and the own-LTR check at a fixed z, SELF_Z (clock.self_ltr_ok).
- element-vs-element divergence: edlib unit cost (seq.align.unit_cost_path), mismatches / aligned columns over the
  two whole records (clock.pair_counts).
- absence: A, inferred: a record inside a collinear window of the other genome without an ortholog there (run).
- inversions, for the gene-flow engine: each collinear block carries its orientation within its chromosome pair
  (blocks.tsv, sign).

Coordinates: genome positions are 1-based closed (a record's start, end and midpoint, the windows' anchor midpoints);
positions inside records are 0-based, in the N-stripped record.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from synltr.io.annotation.tsd import canon_tsd
from synltr.io.library import Library


@dataclass(frozen=True, eq=False)
class Records:
    """One genome's library as the screen reads it, one row per record in library order. An unplaced record (its name
    holds no coordinates) sits on a chromosome named after itself, from 0 to its length: it takes part in the k-mers
    and the anchors but lies in no collinear block, and no output names it."""

    label: str
    lib: Library
    names: np.ndarray  # object: the record names
    chrom: np.ndarray  # object: the seqid; an unplaced record's own name
    start: np.ndarray  # int64, 1-based closed; an unplaced record 0
    end: np.ndarray  # int64; an unplaced record its length
    mid: np.ndarray  # int64: (start + end) // 2
    strand: np.ndarray  # object: '+', '-' or '.' (records of a known strand are in element orientation)
    k2p: np.ndarray  # float64: the annotated LTR K2P, NaN unknown
    tsd: np.ndarray  # object: the annotated TSD in canonical orientation (canon_tsd), "" unknown
    ltr_len: np.ndarray  # int64 (n, 2): the own lengths of the genome-left and genome-right LTRs (0 unknown)

    @classmethod
    def of(cls, label: str, lib: Library) -> Records:
        """The records of `lib`, a library read under io.library's screen rule (LTR bounds in the stripped record)."""
        if (lib.bounds == "mapped").any():
            raise ValueError(f"{label}: not a screen library (io.library rule 'screen'): LTRs mapped from the genome")
        e = lib.elements
        placed = lib.placed
        start = np.where(placed, e.column("start"), 0).astype(np.int64)
        end = np.where(placed, e.column("end"), np.diff(lib.offs)).astype(np.int64)
        chrom = np.empty(lib.n, object)
        chrom[:] = [s if p else nm for s, p, nm in zip(e.column("seqid"), placed.tolist(), lib.names)]
        tsd = np.empty(lib.n, object)
        tsd[:] = [canon_tsd(t) for t in e.column("tsd")]
        return cls(label, lib, lib.names, chrom, start, end, (start + end) // 2, e.column("strand"),
                   e.column("k").astype(np.float64),
                   tsd, np.nan_to_num(np.stack([e.column("ltr_left_len"), e.column("ltr_right_len")], 1),
                                      nan=0.0).astype(np.int64))

    @property
    def n(self) -> int:
        return self.lib.n

    @property
    def placed(self) -> np.ndarray:
        return self.lib.placed

    def seq(self, i: int) -> np.ndarray:
        """Record i's codes, N-stripped (A C G T = 0-3 only)."""
        return self.lib.seq(int(i), stripped=True)

    def ltrs(self, i: int) -> tuple[tuple[int, int], tuple[int, int]] | None:
        """Record i's two LTRs in its N-stripped record (0-based half-open, 5' first), or None."""
        l5, l3 = self.lib.l5[i], self.lib.l3[i]
        if l5[0] < 0:
            return None
        return (int(l5[0]), int(l5[1])), (int(l3[0]), int(l3[1]))
