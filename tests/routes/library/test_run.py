"""The screen's evidence over a samples sheet (routes.library.run.screen_pairs): for each pair, both directions' call
files (P ortholog, A absent inside a collinear window, U outside every window), the ortholog divergences and the
pair's collinear blocks, checkpointed per pair; blocks.tsv; records without coordinates never written with one."""
import gzip
import os
import re
from pathlib import Path

import builders
import numpy as np
import pytest

from ltrtrace.base.checkpoint import done_meta, marker
from ltrtrace.base.log import UsageError
from ltrtrace.calls.records import read_calls, read_divergences
from ltrtrace.io.elements import ELEMENT_COLUMNS, make_eid, read_elements
from ltrtrace.io.samples import read_samples
from ltrtrace.routes.library import run as screen
from ltrtrace.routes.library.run import screen_pairs

PAIRS = [("R", "S"), ("R", "D"), ("S", "D")]
BLOCK_COLUMNS = "genome_a\tchrom_a\tstart_a\tend_a\tgenome_b\tchrom_b\tstart_b\tend_b\tsign\tn_anchor\n"
OLD_NS = 1_000_000_000_000_000_000


def _outputs(run):
    """The screen's files under `run` (relative path -> bytes; gzip files decompressed), done markers aside, and not
    the gene-flow analyses' (work/geneflow: tests share the screens they run on)."""
    out = {}
    for p in sorted(run.rglob("*")):
        if p.is_file() and not p.name.endswith(".done.json") and p.name != "checksums.json" and \
                not p.relative_to(run).as_posix().startswith("work/geneflow/"):
            data = p.read_bytes()
            out[str(p.relative_to(run))] = gzip.decompress(data) if p.suffix == ".gz" else data
    return out


def test_screen_outputs(admixed):
    """Both directions of each pair in CALL_COLUMNS (every element of A once, library order; engine library), a
    divergence row per P call (n_sites aligned columns of the two records, n_diff mismatches; d_kind p_whole_unit),
    the elements tables, and blocks.tsv (the pairs' blocks in pair order)."""
    _, run = admixed
    blocks = (run / "blocks.tsv").read_text().splitlines(keepends=True)
    assert blocks[0] == BLOCK_COLUMNS and len(blocks) > 9
    unresolved = 0
    for a, b in PAIRS:
        for x, y in ((a, b), (b, a)):
            path = run / "pairs" / f"{x}__{y}.calls.tsv.gz"
            elements = list(read_elements(run / "work" / "elements" / f"{x}.tsv"))
            calls = list(read_calls(path))
            assert [c["eid"] for c in calls] == [e.eid for e in elements]
            assert done_meta(marker(path))["engine"] == "library"
            assert {"P", "A"} <= {c["state"] for c in calls} <= {"P", "A", "U"}
            unresolved += sum(c["state"] == "U" for c in calls)
            for c, e in zip(calls, elements):
                assert (c["a_seqid"], c["a_start"], c["a_end"]) == (e.seqid, e.start, e.end)
                assert (c["state"] == "U") == (c["reason"] == "no_block") == (c["wid"] is None)
                if c["state"] != "U":
                    assert c["b_seqid"].startswith(f"{y}_chr") and 1 <= c["b_start"] <= c["b_end"]
                    assert c["b_strand"] in "+-" and c["b_ins"] >= 0
            divs = list(read_divergences(run / "pairs" / f"{x}__{y}.div.tsv.gz"))
            p_calls = [c for c in calls if c["state"] == "P"]
            assert [d["eid"] for d in divs] == [c["eid"] for c in p_calls]
            assert [d["b_eid"] for d in divs] == [make_eid(y, c["b_seqid"], c["b_start"], c["b_end"]) for c in p_calls]
            assert all(d["d"] is None if not d["n_sites"] else abs(d["d"] - d["n_diff"] / d["n_sites"]) < 1e-6
                       for d in divs)
            meta = done_meta(marker(run / "pairs" / f"{x}__{y}.div.tsv.gz"))
            assert (meta["engine"], meta["d_kind"]) == ("library", "p_whole_unit")
        assert sum(line.startswith(f"{a}\t") and f"\t{b}\t" in line for line in blocks) >= 3
    assert unresolved > 0
    head = (run / "work" / "elements" / "R.tsv").read_text().splitlines()[0]
    assert head.split("\t") == list(ELEMENT_COLUMNS)


def test_screen_resume(admixed, tmp_path, monkeypatch):
    """The k-mer stage is one unit, then each pair. A screen killed after the k-mer stage and two of its three pairs
    finished: the rerun does no k-mer work and computes only the third pair (the finished pairs' files untouched),
    and every output is byte-identical to an uninterrupted run. A pair whose last file was never marked done is
    computed again, alone. A screen killed inside the k-mer stage (some of its files written, none marked done)
    computes the stage again, then every pair, with the same outputs. A rerun on the sheet reordered (S before R)
    runs the k-mer stage again (its key holds the labels in sheet order), which removes the earlier order's files,
    and computes only the reordered pair: its outputs equal a fresh run on that order."""
    sc, reference = admixed
    sheet, out = read_samples(sc.sheet), tmp_path / "run"
    real, seen, scans = screen.orthology, [], []
    real_scan = screen.shared_kmers
    monkeypatch.setattr(screen, "shared_kmers", lambda *a, **kw: scans.append(1) or real_scan(*a, **kw))

    def killed_at_third(la, lb, **kw):
        seen.append((la.label, lb.label))
        if len(seen) == 3:
            raise KeyboardInterrupt
        return real(la, lb, **kw)

    monkeypatch.setattr(screen, "orthology", killed_at_third)
    with pytest.raises(KeyboardInterrupt):
        screen_pairs(sheet, out, 2)
    assert seen == PAIRS and scans == [1] and not (out / "blocks.tsv").exists()
    assert not (out / "pairs" / "S__D.calls.tsv.gz").exists()
    finished = [p for p in (out / "pairs").iterdir() if p.name.startswith(("R__", "S__R", "D__R"))]
    assert len(finished) == 16                          # 2 pairs x 2 directions x (calls, div) x (file, marker)
    kmer_files = sorted((out / "work" / "library" / "kmers").iterdir())
    for p in finished + kmer_files:
        os.utime(p, ns=(OLD_NS, OLD_NS))

    seen.clear()
    monkeypatch.setattr(screen, "orthology",
                        lambda la, lb, **kw: seen.append((la.label, lb.label)) or real(la, lb, **kw))
    paths = screen_pairs(sheet, out, 2)
    assert seen == [("S", "D")] and scans == [1]        # no k-mer work: the stage was done
    assert all(p.stat().st_mtime_ns == OLD_NS for p in finished + kmer_files)
    assert paths == [out / "pairs" / f"{x}__{y}.calls.tsv.gz" for a, b in PAIRS for x, y in ((a, b), (b, a))]
    assert _outputs(out) == _outputs(reference)

    marker(out / "pairs" / "D__R.div.tsv.gz").unlink()  # the R-D pair stopped before its last file was marked done
    seen.clear()
    screen_pairs(sheet, out, 2)
    assert seen == [("R", "D")] and scans == [1] and _outputs(out) == _outputs(reference)
    seen.clear()
    screen_pairs(sheet, out, 2)                         # everything done: nothing is computed
    assert seen == [] and scans == [1]

    inside = tmp_path / "inside"                        # killed while the k-mer stage writes its files
    real_save, saved = screen._save_arrays, []

    def killed_at_second_file(path, arrays):
        saved.append(path)
        if len(saved) == 2:
            raise KeyboardInterrupt
        return real_save(path, arrays)

    monkeypatch.setattr(screen, "_save_arrays", killed_at_second_file)
    with pytest.raises(KeyboardInterrupt):
        screen_pairs(sheet, inside, 2)
    assert scans == [1, 1] and len(saved) == 2 and not (inside / "pairs").exists()
    monkeypatch.setattr(screen, "_save_arrays", real_save)
    seen.clear()
    screen_pairs(sheet, inside, 2)
    assert scans == [1, 1, 1] and seen == PAIRS and _outputs(inside) == _outputs(reference)

    swapped = [sheet[1], sheet[0], sheet[2]]           # the sheet reordered, S before R, after a finished run
    seen.clear()
    screen_pairs(swapped, out, 2)
    assert scans == [1, 1, 1, 1] and seen == [("S", "R")]  # the k-mer stage again, and only the reordered pair
    fresh = tmp_path / "fresh"
    screen_pairs(swapped, fresh, 2)
    got, want = _outputs(out), _outputs(fresh)
    assert {name: got.get(name) for name in want} == want   # every output of a fresh run on that order, the same
    kmer_files = sorted(p.name for p in (out / "work" / "library" / "kmers").iterdir())
    assert kmer_files == sorted(p.name for p in (fresh / "work" / "library" / "kmers").iterdir())
    assert "S__R.pairs.npz" in kmer_files and "R__S.pairs.npz" not in kmer_files  # no file of the earlier order
    with pytest.raises(screen.KmerStageError, match="R__S"):
        screen.read_candidates(out, "R", "S")


def test_screen_records_cache(admixed, tmp_path, monkeypatch):
    """screen_pairs fills the caller's records dict with every library it loads and reads none it finds there; the
    trio's single-copy tallies of each member read back (read_tallies) as single_copy_tallies gives them."""
    from ltrtrace.seq.kmers import single_copy_tallies

    sc, _ = admixed
    sheet, loads = read_samples(sc.sheet), []
    real = screen.read_library
    monkeypatch.setattr(screen, "read_library", lambda spec, g, **kw: loads.append(g) or real(spec, g, **kw))
    cache = {}
    screen_pairs(sheet, tmp_path / "a", 2, records=cache)
    assert loads == ["R", "S", "D"] and list(cache) == ["R", "S", "D"]
    screen_pairs(sheet, tmp_path / "b", 2, records=cache)
    assert loads == ["R", "S", "D"]                     # all three reused
    for focal in "RSD":
        want = single_copy_tallies({g: cache[g].lib for g in "RSD"}, focal=focal)
        got = screen.read_tallies(tmp_path / "b", ("R", "S", "D"), focal)
        assert list(got) == list(want) and all(np.array_equal(got[p], want[p]) for p in want)


def test_kmer_files_read_only_as_the_stage_lists_them(admixed, tmp_path):
    """read_candidates and read_tallies read only the files the k-mer stage's done marker lists: before the stage is
    done (no marker, or a listed file changed), or for a pair or trio the stage does not hold (a stray file under its
    name aside), one KmerStageError naming the folder and the pair or trio. read_tallies takes the trio in any order,
    its keys listing the other two in sheet order."""
    import shutil

    _, reference = admixed
    run = tmp_path / "run"
    folder = run / "work" / "library" / "kmers"
    for read in (lambda: screen.read_candidates(run, "R", "S"), lambda: screen.read_tallies(run, ("R", "S", "D"), "R")):
        with pytest.raises(screen.KmerStageError, match=f"{re.escape(str(folder))}.*not done"):
            read()
    shutil.copytree(reference / "work" / "library" / "kmers", folder)
    want = screen.read_tallies(run, ("R", "S", "D"), "S")
    assert list(want) == [(), ("R",), ("D",), ("R", "D")]
    for trio in (("D", "S", "R"), ("S", "R", "D")):
        got = screen.read_tallies(run, trio, "S")
        assert list(got) == list(want) and all(np.array_equal(got[p], want[p]) for p in want)
    assert screen.read_candidates(run, "R", "S").a.size > 0
    shutil.copy(folder / "R__S.pairs.npz", folder / "S__R.pairs.npz")  # a stray file the stage does not list
    for read, name in ((lambda: screen.read_candidates(run, "S", "R"), "S__R"),
                       (lambda: screen.read_candidates(run, "R", "X"), "R__X"),
                       (lambda: screen.read_tallies(run, ("R", "S", "X"), "R"), "R, S, X"),
                       (lambda: screen.read_tallies(run, ("R", "S", "D"), "X"), "R, S, D")):
        with pytest.raises(screen.KmerStageError, match=f"{re.escape(str(folder))}.*{name}"):
            read()
    with open(folder / "R__S.pairs.npz", "ab") as fh:  # a listed file changed after the stage
        fh.write(b"x")
    with pytest.raises(screen.KmerStageError, match="not done"):
        screen.read_candidates(run, "R", "D")


def test_unplaced_ortholog_written_without_coordinates(tmp_path):
    """A P call whose ortholog is a record whose name holds no coordinates writes "." for b_seqid, b_start and b_end,
    and its divergence row "." for b_eid (the counts and d as for any ortholog)."""
    from types import SimpleNamespace

    from ltrtrace.calls.records import CALL_COLUMNS
    from ltrtrace.io.library import read_library
    from ltrtrace.routes.library import Records
    from ltrtrace.routes.library.absence import Windows

    rng = np.random.default_rng(4)
    seq = ["".join("ACGT"[c] for c in rng.integers(0, 4, 60)) for _ in range(3)]
    fa = builders.write_fasta(tmp_path / "A.fa", {"c1:1-60": seq[0], "c1:101-160": seq[1], "c1:201-260": seq[2]})
    fb = builders.write_fasta(tmp_path / "B.fa", {"c1:1-60": seq[0], "contig_7": seq[1], "c1:201-260": seq[2]})
    la, lb = Records.of("A", read_library(str(fa), "A")), Records.of("B", read_library(str(fb), "B"))
    assert lb.placed.tolist() == [True, False, True]
    po = SimpleNamespace(orth_b=np.array([0, 1, -1]), status_a=np.array([2, 2, 1], np.int8),
                         win_a=Windows(*(np.zeros(3, np.int64) for _ in range(4))),
                         counts_a=np.array([[60, 0], [60, 3], [0, 0]]),
                         blocks=[SimpleNamespace(sign=1, chrom_a="c1", chrom_b="c1")])
    lines, divs = screen._calls(po, la, lb, "a")
    calls = [dict(zip(CALL_COLUMNS, line.rstrip("\n").split("\t"))) for line in lines]
    assert [(c["state"], c["b_seqid"], c["b_start"], c["b_end"]) for c in calls] == \
        [("P", "c1", "1", "60"), ("P", ".", ".", "."), ("A", "c1", "0", "0")]
    assert [(d["b_eid"], d["n_sites"], d["n_diff"]) for d in divs] == \
        [(make_eid("B", "c1", 1, 60), 60, 0), (None, 60, 3)] and divs[1]["d"] == 0.05


def test_unplaced_record_never_written_with_a_coordinate(tmp_path):
    """A record whose name holds no coordinates (kept for its k-mers) is never written with one: no call, divergence
    or element row and no block names it, and no output holds a start below 1."""
    sc = builders.trio(tmp_path / "in", seed=3, admix=True, n_chrom=1, chrom_len=1_500_000)
    fa = sc.libraries["R"].split(",")[0]
    first = Path(fa).read_text().split(">")[1]
    with open(fa, "a") as fh:
        fh.write(">contig_unplaced\n" + first.split("\n", 1)[1].replace("\n", "") + "\n")
    sheet = tmp_path / "pair.tsv"
    sheet.write_text(f"label\tlibrary\nR\t{sc.libraries['R']}\nD\t{sc.libraries['D']}\n")
    screen_pairs(read_samples(sheet), tmp_path / "run", 2)
    outputs = _outputs(tmp_path / "run")
    assert not any(b"contig_unplaced" in data for data in outputs.values())
    for x, y in (("R", "D"), ("D", "R")):
        for c in read_calls(tmp_path / "run" / "pairs" / f"{x}__{y}.calls.tsv.gz"):
            assert all(c[k] is None or c[k] >= 1 for k in ("a_start", "a_end", "b_start", "b_end"))
            assert c["eid"] is not None and c["a_seqid"] is not None
    for line in (tmp_path / "run" / "blocks.tsv").read_text().splitlines()[1:]:
        f = line.split("\t")
        assert min(int(f[2]), int(f[3]), int(f[6]), int(f[7])) >= 1
    assert len(list(read_elements(tmp_path / "run" / "work" / "elements" / "R.tsv"))) == \
        Path(fa).read_text().count(">") - 1


def test_screen_needs_a_library_per_genome(tmp_path, admixed):
    sc, _ = admixed
    sheet = tmp_path / "s.tsv"
    sheet.write_text(f"label\tlibrary\tfasta\nR\t{sc.libraries['R']}\t.\nX\t.\t{tmp_path / 'x.fa'}\n")
    (tmp_path / "x.fa").write_text(">c\nACGT\n")
    with pytest.raises(UsageError, match="X has no library"):
        screen_pairs(read_samples(sheet), tmp_path / "run", 1)


def test_written_fields_match_their_sources(admixed):
    """The fields only the call files carry, against PairOrtho and blocks.tsv, both directions of the R-D pair: b_ins is
    the window's room for P and A calls; wid is the window block's row among the pair's rows of blocks.tsv; b_strand
    that block's sign; an A call's b_seqid, b_start and b_end the window on the other genome; a P call's the
    ortholog's coordinates."""
    from ltrtrace.io.library import read_library
    from ltrtrace.routes.library import Records
    from ltrtrace.routes.library.orthology import orthology

    sc, run = admixed
    rec = {g: Records.of(g, read_library(sc.libraries[g], g)) for g in "RD"}
    po = orthology(rec["R"], rec["D"], threads=2)
    rows = [line.split("\t") for line in (run / "blocks.tsv").read_text().splitlines()[1:]]
    rows = [r for r in rows if (r[0], r[4]) == ("R", "D")]
    assert len(rows) == len(po.blocks) > 0
    for k, (r, b) in enumerate(zip(rows, po.blocks)):
        assert (r[1], int(r[2]), int(r[3]), r[5], int(r[6]), int(r[7]), int(r[8])) == \
            (b.chrom_a, b.a_start, b.a_end, b.chrom_b, b.b_start, b.b_end, b.sign), k
    checked = 0
    for x, y, win, orth in (("R", "D", po.win_a, po.orth_b), ("D", "R", po.win_b, po.orth_a)):
        index = {e: i for i, e in enumerate(rec[x].lib.elements.column("eid"))}
        for c in read_calls(run / "pairs" / f"{x}__{y}.calls.tsv.gz"):
            i = index[c["eid"]]
            if c["state"] == "U":
                assert win.bid[i] < 0 and c["wid"] is c["b_ins"] is c["b_strand"] is None
                continue
            k = int(win.bid[i])
            assert (c["wid"], c["b_ins"]) == (k, int(win.room[i]))
            assert c["b_strand"] == ("+" if po.blocks[k].sign > 0 else "-")
            if c["state"] == "A":
                other = po.blocks[k].chrom_b if x == "R" else po.blocks[k].chrom_a
                assert (c["b_seqid"], c["b_start"], c["b_end"]) == (other, int(win.lo[i]), int(win.hi[i]))
            else:
                j = int(orth[i])
                assert (c["b_seqid"], c["b_start"], c["b_end"]) == (rec[y].chrom[j], rec[y].start[j], rec[y].end[j])
            checked += 1
    assert checked > 1000
