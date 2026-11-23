"""Inputs of the subgenomes analysis: the LTR-RT library joined to its annotation (GFF3 or a K2P table), LTR bounds of
clean records, the input checks and the .fai."""
import sys

import numpy as np
import pytest
from analyses.subgenomes import sgtools as sg

from synltr.analyses.subgenomes import inputs
from synltr.base.log import UsageError
from synltr.io.seqio import encode


def test_natural_key_orders_numbers():
    assert sorted(["c10", "c2", "c1"], key=inputs.natural_key) == ["c1", "c2", "c10"]


def test_labels_many_groups():
    assert inputs.group_label(11) == "SG12"


def _write_lib(tmp_path, gff=True, identity=False, drop=0):
    fa = tmp_path / "l.fa"
    fa.write_text(">c1:1-50#LTR/Gypsy/Tekay\n" + "A" * 50 + "\n>c2:101-160#LTR/Copia/Ale\n" + "C" * 60 + "\n"
                  ">c2:201-260\n" + "GT" * 30 + "\n")
    rows = [("c1", 1, 50, 0.01), ("c2", 101, 160, 0.02), ("c2", 201, 260, 0.03)][drop:]
    if gff:
        g = tmp_path / "a.gff3"
        lines = ["##gff-version 3"]
        for i, (c, s, e, k) in enumerate(rows):
            kk = 1 - k if identity else k
            lines.append(f"{c}\tLTRquest\tLTR_retrotransposon\t{s}\t{e}\t.\t+\t.\tID=e{i};Name={c}:{s}-{e};"
                         f"family=f{i};clade=Ale;superfamily=Copia;depth=0;ltr5_len=10;ltr3_len=10;aln_len=10;"
                         f"subs=1;ti=1;tv=0;K2P_d={kk}")
            lines.append(f"{c}\tLTRquest\tlong_terminal_repeat\t{s}\t{s + 9}\t.\t+\t.\tID=e{i}.l;Parent=e{i}")
            lines.append(f"{c}\tLTRquest\tlong_terminal_repeat\t{e - 9}\t{e}\t.\t+\t.\tID=e{i}.r;Parent=e{i}")
        g.write_text("\n".join(lines) + "\n")
        return fa, g
    t = tmp_path / "k.tsv"
    t.write_text("#seq_id\tk2p\tk2p_se\tn_sites\tn_ts\tn_tv\tfamily\n" +
                 "".join(f"{c}:{s}-{e}#x\t{k}\t0.001\t100\t1\t0\tfam{i}\n" for i, (c, s, e, k) in enumerate(rows)))
    return fa, t


@pytest.mark.parametrize("gff", [True, False])
def test_load_library_gff_and_tsv(tmp_path, gff):
    fa, an = _write_lib(tmp_path, gff=gff)
    a = sg.parse_args(["--ltr_fasta", str(fa), "--annotation", str(an), "--outdir", str(tmp_path / "o")])
    lib = inputs.load_library(a)
    assert lib.n == 3 and lib.seqnames == ["c1", "c2"] and list(lib.seq) == [0, 1, 1]
    assert np.allclose(lib.k2p, [0.01, 0.02, 0.03]) and lib.offs[-1] == 170
    assert not lib.fai_given and list(lib.seq_len) == [50, 260]


def test_load_library_rejects_identity_and_mismatch(tmp_path):
    """Divergences that look like identities, and a FASTA most of whose copies the annotation lacks, are usage errors
    naming the file."""
    fa, an = _write_lib(tmp_path, identity=True)
    a = sg.parse_args(["--ltr_fasta", str(fa), "--annotation", str(an), "--outdir", str(tmp_path / "o")])
    with pytest.raises(UsageError, match=f"{an}: the divergence column looks like identity"):
        inputs.load_library(a)
    fa, an = _write_lib(tmp_path, drop=2)                      # only 1 of 3 copies annotated
    a = sg.parse_args(["--ltr_fasta", str(fa), "--annotation", str(an), "--outdir", str(tmp_path / "o")])
    with pytest.raises(UsageError, match="only 1 of 3 FASTA copies are annotated"):
        inputs.load_library(a)


def test_fai_lengths_and_missing_sequence(tmp_path):
    fa, an = _write_lib(tmp_path)
    fai = tmp_path / "g.fai"
    fai.write_text("c1\t1000\t0\t0\t0\nc2\t2000\t0\t0\t0\nc3\t50\t0\t0\t0\n")
    a = sg.parse_args(["--ltr_fasta", str(fa), "--annotation", str(an), "--fai", str(fai), "--outdir", str(tmp_path / "o")])
    lib = inputs.load_library(a)
    assert lib.fai_given and lib.seqnames == ["c1", "c2", "c3"] and list(lib.seq_len) == [1000, 2000, 50]
    fai.write_text("c1\t1000\t0\t0\t0\n")
    with pytest.raises(UsageError, match=f"{fai}: 1 sequence.s. carrying LTR-RTs are absent from the .fai, e.g. c2"):
        inputs.load_library(a)


def _ltr_lib(tmp_path, strand, table=False, dup=False):
    rng = np.random.default_rng(5)
    seq = "".join(rng.choice(list("ACGT"), 1000))
    fa = tmp_path / "l.fa"
    fa.write_text(f">c1:100-1099\n{seq}\n" + (f">c1:100-1099\n{seq}\n" if dup else "") + f">c2:5-1004\n{seq}\n")
    if table:
        an = tmp_path / "k.tsv"
        an.write_text("#seq_id\tk2p\tltr5_start\tltr5_end\tltr3_start\tltr3_end\tstrand\n"
                      f"c1:100-1099\t0.01\t1\t50\t901\t1000\t{strand}\nc2:5-1004\t0.01\t1\t50\t901\t1000\t{strand}\n")
    else:
        an = tmp_path / "a.gff3"
        rows = []
        for i, (c, s, e) in enumerate((("c1", 100, 1099), ("c2", 5, 1004))):
            rows.append(f"{c}\tq\tLTR_retrotransposon\t{s}\t{e}\t.\t{strand}\t.\tID=e{i};K2P_d=0.01")
            rows.append(f"{c}\tq\tlong_terminal_repeat\t{s}\t{s + 49}\t.\t{strand}\t.\tID=e{i}.l;Parent=e{i}")
            rows.append(f"{c}\tq\tlong_terminal_repeat\t{e - 99}\t{e}\t.\t{strand}\t.\tID=e{i}.r;Parent=e{i}")
        an.write_text("\n".join(rows) + "\n")
    return sg.parse_args(["--ltr_fasta", str(fa), "--annotation", str(an), "--outdir", str(tmp_path / "o")])


@pytest.mark.parametrize("strand,table,l5,l3", [("+", False, (0, 50), (900, 1000)), ("-", False, (0, 100), (950, 1000)),
                                                ("+", True, (0, 50), (900, 1000))])
def test_load_library_ltr_offsets(tmp_path, strand, table, l5, l3):
    lib = inputs.load_library(_ltr_lib(tmp_path, strand, table))
    assert tuple(lib.l5[0]) == l5 and tuple(lib.l3[0]) == l3
    assert list(lib.strand_code) == [{"+": 1, "-": -1}[strand]] * 2


def test_duplicate_coordinates_warn_and_keep_first(tmp_path):
    """A record repeating an earlier record's coordinates is a minor input fault (warn, keep the first); more than
    5% of such records is a broken input (fail)."""
    rng = np.random.default_rng(6)
    seqs = ["".join(rng.choice(list("ACGT"), 60)) for _ in range(30)]
    recs = [(f"c{i % 3}", 100 * i + 1, 100 * i + 60) for i in range(30)]
    (tmp_path / "k.tsv").write_text("#seq_id\tk2p\n" + "".join(f"{c}:{s}-{e}\t0.01\n" for c, s, e in recs))
    fa = "".join(f">{c}:{s}-{e}\n{q}\n" for (c, s, e), q in zip(recs, seqs))
    (tmp_path / "l.fa").write_text(fa + f">c0:1-60\n{seqs[5]}\n")                # repeats record 0
    a = sg.parse_args(["--ltr_fasta", str(tmp_path / "l.fa"), "--annotation", str(tmp_path / "k.tsv"), "--outdir",
                       str(tmp_path / "o")])
    lib = inputs.load_library(a)
    assert lib.n == 30 and lib.codes[lib.offs[0]:lib.offs[1]].tolist() == encode(seqs[0].encode()).tolist()
    (tmp_path / "l.fa").write_text(fa + "".join(f">{c}:{s}-{e}\n{q}\n" for (c, s, e), q in zip(recs[:5], seqs)))
    with pytest.raises(UsageError, match="repeat the coordinates of an earlier record"):
        inputs.load_library(a)


LTR = "ACGTTGCAAC"                                         # both LTRs identical: annotated K2P 0


def _rec(n_int, tail=""):
    return LTR + "T" * n_int + LTR + tail


_R = np.random.default_rng(7)


LTR6, INT6 = ("".join(_R.choice(list("ACGT"), n)) for n in (100, 200))   # realistic sizes for the agreement test


def _clean_library(tmp_path):
    """Six copies (genome coordinates, record, GFF3 LTR features, ltr5_len/ltr3_len = genome-left/right LTR):
    e0 '+' clean: span 1-100, record 60 bp (40-bp nest cut from the internal region), LTRs 10/10
    e1 '-' clean: span 201-300, genome-left LTR 12 bp, genome-right 10 bp; record in element orientation, so its
       first LTR is the genome-right one (10 bp) and its last the genome-left one (12 bp)
    e2 '.' clean: like e0 (unknown-strand records are in genome orientation)
    e3 '+' clean, nest inside the 3' LTR: feature 471-500 (30 bp), ltr3_len 10 (after removal)
    e4 '+' clean, LTR lengths sum beyond the record (15 bp): no LTR bounds (whole element)
    e5 '+' whole-span record (100 bp): genome offsets as before
    e6 '+' clean (100-bp LTRs, 200-bp internal region), annotated LTR lengths (150) still include a nest: the record-end
       windows reach 50 bp into the internal region, their K2P does not reproduce the annotated 0: no LTR bounds."""
    rows = [("e0", 1, 100, "+", (1, 10), (91, 100), 10, 10, _rec(40)),
            ("e1", 201, 300, "-", (201, 212), (291, 300), 12, 10, _rec(38, "GG")),
            ("e2", 301, 400, ".", (301, 310), (391, 400), 10, 10, _rec(40)),
            ("e3", 401, 500, "+", (401, 410), (471, 500), 10, 10, _rec(40)),
            ("e4", 601, 700, "+", (601, 610), (691, 700), 10, 10, "ACGTTGCAACTTTTT"),
            ("e5", 801, 900, "+", (801, 810), (891, 900), 10, 10, _rec(80)),
            ("e6", 1001, 1600, "+", (1001, 1150), (1451, 1600), 150, 150, LTR6 + INT6 + LTR6)]
    fa, g = tmp_path / "l.fa", tmp_path / "a.gff3"
    fa.write_text("".join(f">c1:{s}-{e}#LTR/Gypsy/Tekay\n{seq}\n" for _, s, e, _, _, _, _, _, seq in rows))
    lines = ["##gff-version 3"]
    for i, s, e, st, l5, l3, n5, n3, _ in rows:
        lines.append(f"c1\tx\tLTR_retrotransposon\t{s}\t{e}\t.\t{st}\t.\tID={i};K2P_d=0;ltr5_len={n5};ltr3_len={n3}")
        lines.append(f"c1\tx\tlong_terminal_repeat\t{l5[0]}\t{l5[1]}\t.\t{st}\t.\tID={i}.l;Parent={i}")
        lines.append(f"c1\tx\tlong_terminal_repeat\t{l3[0]}\t{l3[1]}\t.\t{st}\t.\tID={i}.r;Parent={i}")
    g.write_text("\n".join(lines) + "\n")
    a = sg.parse_args(["--ltr_fasta", str(fa), "--annotation", str(g), "--outdir", str(tmp_path / "o")])
    return inputs.load_library(a)


@pytest.fixture(scope="module")
def clean_lib(tmp_path_factory):
    return _clean_library(tmp_path_factory.mktemp("clean"))


@pytest.mark.parametrize("i,l5,l3", [(0, (0, 10), (50, 60)),       # '+': record ends, ltr5_len first
                                     (1, (0, 10), (48, 60)),       # '-': ltr3_len (genome-right LTR) first
                                     (2, (0, 10), (50, 60)),       # '.': genome orientation, as '+'
                                     (3, (0, 10), (50, 60)),       # nest inside an LTR: attribute length, not feature
                                     (5, (0, 10), (90, 100))])     # whole-span record: genome offsets (unchanged)
def test_clean_record_ltr_bounds(clean_lib, i, l5, l3):
    assert tuple(clean_lib.l5[i]) == l5 and tuple(clean_lib.l3[i]) == l3


@pytest.mark.parametrize("i", [4, 6])                               # too short; lengths that include a nest
def test_clean_record_without_matching_ltrs_keeps_whole_element(clean_lib, i):
    assert tuple(clean_lib.l5[i]) == (-1, -1) and tuple(clean_lib.l3[i]) == (-1, -1)


def test_clean_record_ltrs_hold_the_ltr_sequence(clean_lib):
    enc = lambda x: encode(x.encode())
    for i, last in ((0, LTR), (1, LTR + "GG"), (2, LTR), (3, LTR)):
        rec = clean_lib.codes[clean_lib.offs[i]:clean_lib.offs[i + 1]]
        a5, b5 = clean_lib.l5[i]
        a3, b3 = clean_lib.l3[i]
        assert (rec[a5:b5] == enc(LTR)).all() and b5 - a5 == len(LTR)
        assert (rec[a3:b3] == enc(last)).all() and b3 - a3 == len(last)


def test_clean_bounds_with_threads_in_a_fresh_process(tmp_path):
    """The clean-record check is the first numba use of a run; with -t > 1 it must not first compile inside forked
    workers (they scan package metadata through shared file handles: 'BadZipFile: Overlapped entries')."""
    import subprocess
    rng = np.random.default_rng(3)
    fa, g = tmp_path / "l.fa", tmp_path / "a.gff3"
    recs, lines = [], ["##gff-version 3"]
    for k in range(900):                                  # 450 mapped + 450 clean copies: several chunks each
        ltr, inner = ("".join(rng.choice(list("ACGT"), n)) for n in (120, 300))
        s0 = 1 + k * 2000
        span = 1000 if k % 2 else len(ltr) * 2 + len(inner)
        e0 = s0 + span - 1
        recs.append(f">c1:{s0}-{e0}#LTR/Gypsy/Tekay\n{ltr + inner + ltr}\n")
        lines += [f"c1\tx\tLTR_retrotransposon\t{s0}\t{e0}\t.\t+\t.\tID=e{k};K2P_d=0;ltr5_len=120;ltr3_len=120",
                  f"c1\tx\tlong_terminal_repeat\t{s0}\t{s0 + 119}\t.\t+\t.\tID=e{k}.l;Parent=e{k}",
                  f"c1\tx\tlong_terminal_repeat\t{e0 - 119}\t{e0}\t.\t+\t.\tID=e{k}.r;Parent=e{k}"]
    fa.write_text("".join(recs))
    g.write_text("\n".join(lines) + "\n")
    code = ("import numpy as np; from synltr.analyses.subgenomes import inputs as v, run; "
            "from synltr.base import runtime; from synltr.seq import align; "
            "ppe = runtime.ProcessPoolExecutor; seen = []; "
            "runtime.ProcessPoolExecutor = lambda *x, **kw: "
            "(seen.append(len(align._band.signatures)), ppe(*x, **kw))[1]; "
            f"a = run.options([{str(fa)!r}], [{str(g)!r}], {str(tmp_path / 'o')!r}, threads=4); "
            "run.setup(a); lib = v.load_library(a); print('bounds', int((lib.l5[:, 0] >= 0).sum()), "
            "'compiled before forking', all(seen) and len(seen) > 0)")
    import os
    env = dict(os.environ, NUMBA_CACHE_DIR=str(tmp_path / "numba"))   # nothing compiled yet, as for a new install
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=600, env=env)
    assert r.returncode == 0, r.stderr[-2000:]
    assert "bounds 900" in r.stdout and "compiled before forking True" in r.stdout, r.stdout


def test_log_lines_go_through_the_package_logger(caplog):
    """The run's log lines are records of the package's logger (synltr.analyses.subgenomes): milestones at INFO,
    warnings at WARNING, per-step progress at DEBUG (shown with -v only); inside `quiet` none is emitted."""
    import logging

    caplog.set_level(logging.DEBUG, logger="synltr")
    inputs.log("a milestone")
    inputs.warn("a warning")
    inputs.vlog("a step")
    with inputs.quiet():
        inputs.log("hidden")
        inputs.warn("hidden")
    got = [(r.name, r.levelno, r.getMessage()) for r in caplog.records]
    assert got == [("synltr.analyses.subgenomes", logging.INFO, "a milestone"),
                   ("synltr.analyses.subgenomes", logging.WARNING, "a warning"),
                   ("synltr.analyses.subgenomes", logging.DEBUG, "a step")]
    caplog.clear()
    logging.getLogger("synltr").setLevel(logging.INFO)   # the default verbosity: no per-step lines
    inputs.vlog("a step")
    assert not caplog.records and not inputs.verbose()
