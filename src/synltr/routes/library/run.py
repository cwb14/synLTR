"""The library route over a samples sheet (column library): each genome's library, then each pair's orthology written
as pair call files, the screen's evidence for the gene-flow engine.

Stages, each checkpointed (base.checkpoint; a rerun skips what is done):
  elements  each genome's elements, its placed records in library order: <outdir>/work/elements/<label>.tsv
  k-mers    one bounded scan over every library (seq.kmers.shared_kmers), one unit of work: every pair's candidate
            record pairs, <outdir>/work/library/kmers/<A>__<B>.pairs.npz (a, b, w), and every trio's single-copy
            tallies of each member, <outdir>/work/library/kmers/<R>__<S>__<D>.tallies.npz (labels in sheet order),
            under one done marker (kmers.done.json) keyed by every library's files, the labels in sheet order and the
            code. The stage first removes the folder's earlier files; `read_candidates` and `read_tallies` read only
            the files its done marker lists (a KmerStageError otherwise)
  pairs     each pair of genomes (sheet order, A before B), one unit of work: orthology.orthology on its candidates,
            then
              <outdir>/pairs/<A>__<B>.calls.tsv.gz and <B>__<A>.calls.tsv.gz  both directions' calls
              <outdir>/pairs/<A>__<B>.div.tsv.gz and <B>__<A>.div.tsv.gz      the ortholog divergences
              <outdir>/work/library/<A>__<B>.blocks.tsv                         the pair's collinear blocks
            each with a done marker keyed by both libraries' files (checksums, in order), the labels and the code of
            this package's import closure. A pair is skipped when all five hold the key, else computed again whole.
then <outdir>/blocks.tsv: every pair's blocks in pair order (columns genome_a chrom_a start_a end_a genome_b chrom_b
start_b end_b sign n_anchor; positions are anchor midpoints, 1-based). A library is read (io.library.read_library)
only when a stage of its genome needs it, and once per caller that passes its `records` dict.

Call rows (calls.records.CALL_COLUMNS, meta engine "library"; the columns' meaning on each engine is in calls.records),
every element of A once, in library order:
  eid, a_seqid, a_start, a_end, parent, level   A's element (1-based closed)
  state      P an ortholog in B; A inside a collinear window of B without one (an inferred absence, never an observed
             empty site); U in no window (reason no_block)
  b_seqid, b_start, b_end   P: the ortholog (1-based closed; "." when the ortholog is a record whose name holds no
             coordinates); A: the window on B searched for the element (a search window, not a site): the midpoints of
             the anchors flanking the element, one anchor of slack each side
  b_strand   P, A: the window block's orientation within its chromosome pair (+, -)
  wid        P, A: the window block, its row among the pair's rows of blocks.tsv (0-based)
  b_ins      P, A: room, B's span between the midpoints of the element's two immediately flanking anchors (bp; not
             inserted bases)
  tsd_a, tsd_b   the annotated TSDs of A's element and of its ortholog, in canonical orientation
  tier "."; the other columns (flags, coverages, ident, anchors, margin) are the assembly route's: none here.
Divergence rows (calls.records.DIV_COLUMNS, meta engine "library", d_kind "p_whole_unit"), one per P call, in call
order: eid, b_eid (the ortholog; "." when its name holds no coordinates), n_sites (aligned columns of the two whole
N-stripped records, edlib unit cost), n_diff (mismatches among them), d = n_diff / n_sites to 6 significant digits
(none without a column; the exact value is n_diff / n_sites).
A record whose name holds no coordinates is no element: it has no row of its own in any of these files.
"""
from __future__ import annotations

import hashlib
import logging
import time
import zipfile
from pathlib import Path

import numpy as np

from synltr.base.checkpoint import (
    atomic_path,
    atomic_write,
    clear_done,
    code_fingerprint,
    done_outputs,
    file_checksum,
    is_done,
    mark_done,
    marker,
    stage_key,
)
from synltr.base.constants import CHECKSUM_CACHE
from synltr.base.log import SynltrError, UsageError
from synltr.calls.records import CALL_COLUMNS, call_line, write_call_lines, write_divergences
from synltr.io.elements import write_elements
from synltr.io.library import library_files, read_library
from synltr.io.samples import PAIR_SEP, Sample
from synltr.routes.library import Records
from synltr.routes.library.orthology import PairOrtho, orthology
from synltr.seq.kmers import PairWeights, shared_kmers

log = logging.getLogger(__name__)

ENGINE = "library"  # the meta of every call and divergence file this route writes
DIV_META = {"engine": ENGINE, "d_kind": "p_whole_unit"}  # d: whole-record mismatches / aligned columns, unit cost
BLOCK_COLUMNS = ("genome_a", "chrom_a", "start_a", "end_a", "genome_b", "chrom_b", "start_b", "end_b", "sign",
                 "n_anchor")
_STATE = {2: "P", 1: "A", 0: "U"}
_PAIRS, _TALLIES = ".pairs.npz", ".tallies.npz"  # the k-mer stage's file suffixes


class KmerStageError(SynltrError):
    """The k-mer stage's files cannot be read: the stage is not done in the folder, or it holds no such pair or trio."""


def screen_pairs(samples: list[Sample], outdir: Path, threads: int, *,
                 records: dict[str, Records] | None = None) -> list[Path]:
    """The screen's evidence for every pair of `samples` (each with a library) under `outdir`, with `threads` in all
    (module docstring). `records` (label -> Records), when given, is filled with every library read here and lends
    every library it already holds. Returns the call files, both directions of each pair in pair order. A sample
    without a library, or fewer than two samples, is a UsageError (`check_screen`)."""
    t0 = time.perf_counter()
    outdir = Path(outdir)
    check_screen(samples)
    records = {} if records is None else records
    labels = [s.label for s in samples]
    code = code_fingerprint("routes.library")
    digest = {s.label: _digest(s.library, outdir / "work" / CHECKSUM_CACHE) for s in samples}
    elements = {s.label: (_elements_path(outdir, s.label),
                          stage_key("screen_elements", {s.label: digest[s.label]}, {"genome": s.label}, code))
                for s in samples}
    kmer_key = stage_key("screen_kmers", digest, {"labels": labels}, code)
    kmer_done = is_done(_kmer_marker(outdir), kmer_key)
    pairs = [(x, y) for i, x in enumerate(samples) for y in samples[i + 1:]]
    keys = {(x.label, y.label): stage_key("screen_pair", {x.label: digest[x.label], y.label: digest[y.label]},
                                          {"pair": [x.label, y.label]}, code) for x, y in pairs}
    todo = [p for p in pairs if not all(is_done(marker(f), keys[p[0].label, p[1].label])
                                        for f in pair_files(outdir, p[0].label, p[1].label))]
    need = ({s.label for p in todo for s in p} | {x for x, (path, key) in elements.items()
                                                   if not is_done(marker(path), key)}
            | (set() if kmer_done else set(labels)))
    log.info("screen: %d genome(s), %d pair(s) (%d done), k-mer stage %s, %d thread(s) -> %s", len(samples),
             len(pairs), len(pairs) - len(todo), "done" if kmer_done else "to run", threads, outdir)
    for s in samples:
        if s.label in need and s.label not in records:
            records[s.label] = _load(s, threads)
    for label, (path, key) in elements.items():
        if not is_done(marker(path), key):
            _write_elements(records[label], path, key)
    if not kmer_done:
        _kmer_stage(outdir, [records[x] for x in labels], kmer_key, threads)
    out = []
    for n, (x, y) in enumerate(pairs, 1):
        files = pair_files(outdir, x.label, y.label)
        if (x, y) in todo:
            log.info("pair %d/%d: %s vs %s", n, len(pairs), x.label, y.label)
            for f in files:
                clear_done(marker(f))
            po = orthology(records[x.label], records[y.label], threads=threads,
                           candidates=read_candidates(outdir, x.label, y.label))
            _write_pair(po, records[x.label], records[y.label], files, keys[x.label, y.label])
        else:
            log.info("pair %d/%d: %s vs %s done, skip", n, len(pairs), x.label, y.label)
        out += files[:2]
    _write_blocks(outdir, [(x.label, y.label) for x, y in pairs])
    log.info("screen: %d pair(s) in %.1f s", len(pairs), time.perf_counter() - t0)
    return out


def check_screen(samples: list[Sample]) -> None:
    """The screen's checks on its samples, made before any work (UsageError): every sample has a library, and there
    are at least two."""
    none = [s.label for s in samples if not s.library]
    if none:
        raise UsageError(f"the library screen reads LTR-RT libraries: {', '.join(none)} has no library in the "
                         "samples sheet (column library: an LTR_peel prefix, a FASTA, or FASTA,TABLE)")
    if len(samples) < 2:
        raise UsageError("the library screen needs at least two genomes")


def _kmer_dir(outdir: Path) -> Path:
    return Path(outdir) / "work" / "library" / "kmers"


def _kmer_marker(outdir: Path) -> Path:
    return marker(_kmer_dir(outdir) / "kmers")


def _kmer_stage(outdir: Path, recs: list[Records], key: str, threads: int) -> None:
    """The k-mer stage (module docstring): the folder's earlier files removed, every pair's candidates and every
    trio's tallies written, then the done marker."""
    t0 = time.perf_counter()
    clear_done(_kmer_marker(outdir))
    d, files = _kmer_dir(outdir), []
    d.mkdir(parents=True, exist_ok=True)
    for p in d.iterdir():                               # an earlier stage's files (another sheet order, a kill)
        if p.is_file():
            p.unlink()
    shared = shared_kmers({r.label: r.lib for r in recs}, threads=threads)
    for (a, b), pw in shared.pairs.items():
        files.append(d / f"{a}{PAIR_SEP}{b}{_PAIRS}")
        _save_arrays(files[-1], {"a": pw.a, "b": pw.b, "w": pw.w})
    for trio, by_focal in shared.tallies.items():
        files.append(d / f"{PAIR_SEP.join(trio)}{_TALLIES}")
        _save_arrays(files[-1], {f"{focal}:{'+'.join(pattern)}": counts for focal, tally in by_focal.items()
                                 for pattern, counts in tally.items()})
    mark_done(_kmer_marker(outdir), key, files, {"passes": shared.passes})
    log.info("k-mer stage: %d librar%s, %d pair(s), %d trio(s), %d pass(es), %.1f s", len(recs),
             "y" if len(recs) == 1 else "ies", len(shared.pairs), len(shared.tallies), shared.passes,
             time.perf_counter() - t0)


def _save_arrays(path: Path, arrays: dict[str, np.ndarray]) -> None:
    """An uncompressed .npz of `arrays` (np.load reads it) with fixed member times, so the same arrays give the same
    bytes; written to a temporary file renamed into place."""
    with atomic_path(path) as tmp, zipfile.ZipFile(tmp, "w", zipfile.ZIP_STORED, allowZip64=True) as zf:
        for name, arr in arrays.items():
            with zf.open(zipfile.ZipInfo(f"{name}.npy", date_time=(1980, 1, 1, 0, 0, 0)), "w",
                         force_zip64=True) as fh:
                np.lib.format.write_array(fh, np.ascontiguousarray(arr), allow_pickle=False)


def _kmer_files(outdir: Path) -> dict[str, Path]:
    """The files the finished k-mer stage under `outdir` lists, by name; a KmerStageError naming the folder when the
    stage is not done there (no done marker, or a file it lists is gone or changed)."""
    files = done_outputs(_kmer_marker(outdir))
    if files is None:
        raise KmerStageError(f"{_kmer_dir(outdir)}: the k-mer stage is not done (no done marker, or a file it lists is "
                             "gone or changed): run the screen there first")
    return {p.name: p for p in files}


def read_candidates(outdir: Path, a: str, b: str) -> PairWeights:
    """Pair (a, b)'s candidate record pairs from the finished k-mer stage (a before b in the sheet the stage ran on),
    read only from a file its done marker lists; a KmerStageError naming the folder and the pair otherwise."""
    files, name = _kmer_files(outdir), f"{a}{PAIR_SEP}{b}"
    if name + _PAIRS not in files:
        raise KmerStageError(f"{_kmer_dir(outdir)}: the k-mer stage holds no pair {name} (a pair is named in the order "
                             "of the sheet the stage ran on)")
    with np.load(files[name + _PAIRS]) as z:
        return PairWeights(z["a"], z["b"], z["w"])


def read_tallies(outdir: Path, trio: tuple[str, str, str], focal: str) -> dict[tuple[str, ...], np.ndarray]:
    """One trio's single-copy tallies of one member (`focal`), as seq.kmers.single_copy_tallies gives them on the
    trio's three libraries in the order of the sheet the stage ran on: {(), (o1,), (o2,), (o1, o2)} -> int64 counts per
    record of the focal genome, o1 and o2 the other two in that order, whatever order `trio` lists the three in. Read
    only from a file the finished stage's done marker lists; a KmerStageError naming the folder and the trio when the
    stage holds no such trio or `focal` is not in it."""
    files = _kmer_files(outdir)
    found = [n[:-len(_TALLIES)].split(PAIR_SEP) for n in files if n.endswith(_TALLIES)]
    order = next((t for t in found if sorted(t) == sorted(trio)), None)
    if order is None or len(set(trio)) != 3 or focal not in order:
        raise KmerStageError(f"{_kmer_dir(outdir)}: the k-mer stage holds no tallies of {focal} in the trio "
                             f"{', '.join(trio)}")
    others = [x for x in order if x != focal]
    patterns = [(), (others[0],), (others[1],), tuple(others)]
    with np.load(files[PAIR_SEP.join(order) + _TALLIES]) as z:
        return {p: z[f"{focal}:{'+'.join(p)}"] for p in patterns}


def pair_files(outdir: Path, a: str, b: str) -> list[Path]:
    """The five files of pair (a, b): both directions' calls, both directions' divergences, the pair's blocks."""
    p, ab, ba = Path(outdir) / "pairs", f"{a}{PAIR_SEP}{b}", f"{b}{PAIR_SEP}{a}"
    return [p / f"{ab}.calls.tsv.gz", p / f"{ba}.calls.tsv.gz", p / f"{ab}.div.tsv.gz", p / f"{ba}.div.tsv.gz",
            Path(outdir) / "work" / "library" / f"{ab}.blocks.tsv"]


def _elements_path(outdir: Path, label: str) -> Path:
    return Path(outdir) / "work" / "elements" / f"{label}.tsv"


def _digest(spec: str, sums: Path) -> str:
    """sha256 of a library's files' checksums, in order (its FASTAs, then its tables)."""
    fastas, tables = library_files(spec)
    return hashlib.sha256("\n".join(file_checksum(f, sums) for f in (*fastas, *tables)).encode()).hexdigest()


def library_digest(spec: str, outdir: Path) -> str:
    """The key a library (`spec`, the samples sheet's library column) gives a checkpoint of the screen in `outdir`:
    sha256 of its files' checksums, in order, the checksums cached in the screen's work folder."""
    return _digest(spec, Path(outdir) / "work" / CHECKSUM_CACHE)


def _load(s: Sample, threads: int) -> Records:
    """Sample s's library as the screen reads it, with one summary line."""
    t0 = time.perf_counter()
    lib = read_library(s.library, s.label, threads=threads)
    rec = Records.of(s.label, lib)
    kinds = {k: int((lib.clean == k).sum()) for k in ("cut", "masked")}
    log.info("%s: %d records (%d unplaced), %.0f Mb library, K2P known %.2f, TSD known %.2f; clean records: %d cut, "
             "%d N-masked (%.1f s)", s.label, rec.n, int((~lib.placed).sum()), lib.codes.size / 1e6,
             float(np.isfinite(rec.k2p).mean()), float((rec.tsd != "").mean()), kinds["cut"], kinds["masked"],
             time.perf_counter() - t0)
    return rec


def _write_elements(rec: Records, path: Path, key: str) -> None:
    """The genome's elements (its placed records, library order) as elements.tsv, marked done."""
    clear_done(marker(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    with atomic_path(path) as tmp:
        write_elements((e for e, placed in zip(rec.lib.elements, rec.placed.tolist()) if placed), tmp)
    mark_done(marker(path), key, [path])


def _calls(po: PairOrtho, lx: Records, ly: Records, side: str):
    """The call lines and divergence rows of lx's elements in ly (module docstring)."""
    orth, status, win = (po.orth_b, po.status_a, po.win_a) if side == "a" else (po.orth_a, po.status_b, po.win_b)
    counts = po.counts_a if side == "a" else np.zeros((lx.n, 2), np.int64)
    if side == "b":
        has = orth >= 0
        counts[has] = po.counts_a[orth[has]]
    ex, ey = lx.lib.elements, ly.lib.elements
    eid, seqid, start, end, parent, level = (ex.column(c) for c in ("eid", "seqid", "start", "end", "parent",
                                                                     "level"))
    b_eid = ey.column("eid")
    lines, divs = [], []
    for i in np.flatnonzero(lx.placed).tolist():
        state, k = _STATE[int(status[i])], int(win.bid[i])
        row = dict.fromkeys(CALL_COLUMNS)
        row.update(eid=eid[i], a_seqid=seqid[i], a_start=int(start[i]), a_end=int(end[i]), parent=parent[i],
                   level=int(level[i]), state=state, tsd_a=lx.tsd[i] or None)
        if state == "U":
            row["reason"] = "no_block"
        elif k >= 0:
            blk = po.blocks[k]
            row.update(b_strand="+" if blk.sign > 0 else "-", wid=k, b_ins=int(win.room[i]))
            if state == "A":
                row.update(b_seqid=blk.chrom_b if side == "a" else blk.chrom_a, b_start=int(win.lo[i]),
                           b_end=int(win.hi[i]))
        j = int(orth[i])
        if state == "P":
            if ly.placed[j]:
                row.update(b_seqid=ly.chrom[j], b_start=int(ly.start[j]), b_end=int(ly.end[j]))
            row["tsd_b"] = ly.tsd[j] or None
            sites, diff = (int(v) for v in counts[i])
            divs.append({"eid": eid[i], "b_eid": b_eid[j], "n_sites": sites, "n_diff": diff,
                         "d": diff / sites if sites else None})
        lines.append(call_line(row[c] for c in CALL_COLUMNS))
    return lines, divs


def _write_pair(po: PairOrtho, la: Records, lb: Records, files: list[Path], key: str) -> None:
    """The pair's five files (pair_files order), each marked done with `key`."""
    meta = {"engine": ENGINE, "key": key}
    for side, (lx, ly), calls, div in (("a", (la, lb), files[0], files[2]), ("b", (lb, la), files[1], files[3])):
        lines, divs = _calls(po, lx, ly, side)
        write_call_lines(calls, lines, meta)
        write_divergences(div, divs, {**DIV_META, "key": key})
        status = (po.status_a if side == "a" else po.status_b)[lx.placed]
        log.debug("%s in %s: %d P, %d A, %d U -> %s", lx.label, ly.label, *np.bincount(status, minlength=3)[::-1],
                  calls)
    blocks = files[4]
    blocks.parent.mkdir(parents=True, exist_ok=True)
    atomic_write(blocks, "\t".join(BLOCK_COLUMNS) + "\n" + "".join(
        f"{la.label}\t{b.chrom_a}\t{b.a_start}\t{b.a_end}\t{lb.label}\t{b.chrom_b}\t{b.b_start}\t{b.b_end}\t"
        f"{b.sign}\t{b.n_anchor}\n" for b in po.blocks))
    mark_done(marker(blocks), key, [blocks])


def _write_blocks(outdir: Path, pairs: list[tuple[str, str]]) -> None:
    """blocks.tsv: the rows of every pair's blocks file, in pair order, under one header."""
    rows = []
    for a, b in pairs:
        with open(pair_files(outdir, a, b)[4]) as fh:
            next(fh)
            rows.extend(fh)
    atomic_write(Path(outdir) / "blocks.tsv", "\t".join(BLOCK_COLUMNS) + "\n" + "".join(rows))
