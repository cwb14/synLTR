import gzip
import logging
from dataclasses import replace

import pytest

from synltr.base.constants import RECORD_WARNINGS
from synltr.io.annotation.ltrquest import load_ltrquest, load_table
from synltr.io.elements import ELEMENT_COLUMNS, NO_LTR, Element, make_eid, read_elements, write_elements

HEADER = '#seq_id\tseq_len\tstatus\tltr5_start\tltr5_end\tltr3_start\tltr3_end\tltr5_len\tltr3_len\tflank5_len\tflank3_len\taln_len\tn_sites\tn_ts\tn_tv\tn_gapcols\tidentity\tp_dist\tk2p\tk2p_se\tbitscore\tflank_margin_bits\tcigar\tmotif\tk2p_time\torientation\ttsd\ttsd_offset\ttsd_input\tstrand\tfamily\tconfidence\tdomains\tnest_status'
ROW_INNER = 'B73v5_chr1:60664-70130#LTR/Copia/SIRE\t9467\tpass\t1\t1322\t8163\t9467\t1322\t1305\t0\t0\t1324\t1303\t0\t0\t21\t1\t0\t0\t0\t1253\tNA\t23=2D861=19I419=\ttg...ca\t0\t-\tCACTC\t0,0\tCACTC\t-\tmerged_fam00001\thigh\tRT|SIRE@63557-64321;RH|SIRE@62912-63289;INT|SIRE@65048-65650;PROT|SIRE@65837-66046;GAG|SIRE@67917-68198\tnest-inner:B73v5_chr1:51712-70182'
ROW_FAIL = 'B73v5_chr1:100-9566#LTR/Copia/SIRE\t9467\tfail\t1\t1322\t8163\t9467\t1322\t1305\t0\t0\t1324\t1303\t0\t0\t21\t1\t0\t0\t0\t1253\tNA\t23=2D861=19I419=\ttg...ca\t0\t-\tCACTC\t0,0\tCACTC\t-\tmerged_fam00001\thigh\tRT|SIRE@63557-64321;RH|SIRE@62912-63289;INT|SIRE@65048-65650;PROT|SIRE@65837-66046;GAG|SIRE@67917-68198\tnest-inner:B73v5_chr1:51712-70182'
ROW_HOST = 'B73v5_chr1:51712-70182#LTR/Copia/SIRE\t18471\tpass\t1\t1294\t7751\t18471\t1294\t1254\t0\t0\t1263\t1180\t4\t0\t83\t0.99661\t0.00338983\t0.00340137\t0.00170359\t987\tNA\t21=1I442=1X120=1I82=19D95=40I89=1X1=1X126=22I17=1X183=\ttg...ca\t56690\t-\tCACCT\t0,0\tCACCT\t-\tmerged_fam00001\thigh\tRT|SIRE@54531-55295;RH|SIRE@53886-54263;INT|SIRE@55974-56576;PROT|SIRE@56769-56978;GAG|SIRE@58654-58935\tnest-outer:B73v5_chr1:60664-70130'


@pytest.fixture(autouse=True)
def _files(tmp_path):
    (tmp_path / "X_depth0_clean_ltr.tsv").write_text(HEADER + "\n" + ROW_INNER + "\n" + ROW_FAIL + "\n")
    (tmp_path / "X_depth1_clean_ltr.tsv").write_text(HEADER + "\n" + ROW_HOST + "\n")


def test_load_ltrquest(tmp_path):
    els, skipped = load_ltrquest(str(tmp_path / "X"), "B73")
    inner = next(e for e in els if e.start == 60664)
    assert (inner.seqid, inner.end, inner.strand) == ("B73v5_chr1", 70130, "-")
    assert inner.ltr_left == (60664, 61985) and inner.ltr_right == (68826, 70130)
    assert inner.tsd == "CACTC" and inner.k == 0.0 and inner.family == "merged_fam00001"
    assert (inner.superfamily, inner.clade, inner.confidence) == ("Copia", "SIRE", "high")
    assert inner.eid == "B73|B73v5_chr1:60664-70130" and inner.name == "B73v5_chr1:60664-70130#LTR/Copia/SIRE"
    host = next(e for e in els if e.start == 51712)
    assert host.ltr_left == (51712, 53005) and host.ltr_right == (59462, 70182)
    assert skipped["status_not_pass"] == 1
    assert len(els) == 2


def test_missing_files(tmp_path):
    with pytest.raises(FileNotFoundError, match="_depth"):
        load_ltrquest(str(tmp_path / "nope"), "B73")


def test_bad_seq_id_and_na(tmp_path):
    c = ROW_INNER.split("\t")
    c[0] = "garbage#LTR/unknown"
    bad = "\t".join(c)
    c = ROW_INNER.split("\t")
    c[0] = "chrZ:5-99#LTR"
    c[HEADER.lstrip("#").split("\t").index("tsd")] = "NA"
    c[HEADER.lstrip("#").split("\t").index("k2p")] = "NA"
    c[HEADER.lstrip("#").split("\t").index("strand")] = "?"
    c[HEADER.lstrip("#").split("\t").index("family")] = "NA"
    na = "\t".join(c)
    (tmp_path / "X_depth0_clean_ltr.tsv").write_text(HEADER + "\n" + bad + "\n" + na + "\n")
    (tmp_path / "X_depth1_clean_ltr.tsv").unlink()
    els, skipped = load_ltrquest(str(tmp_path / "X"), "G")
    assert skipped["bad_seq_id"] == 1 and len(els) == 1
    e = els[0]
    assert (e.tsd, e.k, e.strand, e.family) == (None, None, ".", None)
    assert (e.superfamily, e.clade) == (None, None)


def test_elements_tsv_round_trip(tmp_path):
    els, _ = load_ltrquest(str(tmp_path / "X"), "B73")
    extra = Element("g", make_eid("g", "c", 1, 10), "c:1-10#LTR", "c", 1, 10, ".", (1, 2), (9, 10),
                    None, 0.0123456789012345, None, None, None, None, parent="g|c:0-99", level=2)
    els = els + [extra]
    p = tmp_path / "e.tsv"
    write_elements(els, p)
    assert p.read_text().splitlines()[0] == "\t".join(ELEMENT_COLUMNS)
    # elements.tsv holds LTRpav's columns only: the fields beyond them are not written
    bare = [replace(e, ltr_left_len=None, ltr_right_len=None, k_se=None, n_sites=None, n_ts=None, n_tv=None,
                    k_kind=None, tsd_source=None)
            for e in els]
    assert list(read_elements(p)) == bare
    assert ELEMENT_COLUMNS[:4] == ("genome", "eid", "name", "seqid")


# ---------------------------------------------------------------- malformed records

COLS = HEADER.lstrip("#").split("\t")


def _row(**fields):
    c = ROW_INNER.split("\t")
    for k, v in fields.items():
        c[COLS.index(k)] = v
    return "\t".join(c)


def test_missing_required_column_fails_with_its_name(tmp_path):
    cut = COLS.index("ltr3_end")
    drop = lambda line: "\t".join(x for i, x in enumerate(line.split("\t")) if i != cut)
    (tmp_path / "X_depth0_clean_ltr.tsv").write_text(drop(HEADER) + "\n" + drop(ROW_INNER) + "\n")
    with pytest.raises(ValueError, match=r"X_depth0_clean_ltr\.tsv: .*lacks required column\(s\) ltr3_end"):
        load_ltrquest(str(tmp_path / "X"), "G")


def test_malformed_rows_warn_and_skip(tmp_path, caplog):
    """A short row, LTR bounds or k2p that are not numbers, a seq_id with start > end or of another form are
    skipped with a warning naming the file and line; NA and nan k2p are unknown; CRLF line ends are read."""
    rows = [ROW_INNER,  # line 2
            "B73v5_chr1:1-2#LTR\tpass",  # line 3: short
            _row(ltr5_end="1,322"),  # line 4
            _row(k2p="0.0O1"),  # line 5
            _row(seq_id="B73v5_chr1:9000-100#LTR/Copia/SIRE"),  # line 6: start > end
            _row(seq_id="B73v5_chr1#LTR"),  # line 7
            _row(seq_id="B73v5_chr2:60664-70130#LTR/Copia/SIRE", k2p="nan")]  # line 8: good
    p = tmp_path / "X_depth0_clean_ltr.tsv"
    p.write_text(HEADER + "\r\n" + "\r\n".join(rows) + "\r\n")
    (tmp_path / "X_depth1_clean_ltr.tsv").unlink()
    with caplog.at_level(logging.WARNING, logger="synltr.io.annotation"):
        els, skipped = load_ltrquest(str(tmp_path / "X"), "G")
    assert [(e.seqid, e.k, e.confidence) for e in els] == [("B73v5_chr1", 0.0, "high"), ("B73v5_chr2", None, "high")]
    assert skipped == {"malformed_columns": 1, "malformed_ltr_bounds": 1, "malformed_k2p": 1,
                       "malformed_coordinates": 1, "bad_seq_id": 1}
    warned = [r.getMessage() for r in caplog.records]
    assert [w.split(": ")[0] for w in warned][:RECORD_WARNINGS] == [f"{p}:{n}" for n in (3, 4, 5, 6, 7)]


def test_depths_sort_as_numbers_and_prefix_is_literal(tmp_path):
    """depth10 comes after depth2 (integer order), and a prefix with glob characters is matched literally."""
    d = tmp_path / "run[1]"
    d.mkdir()
    for depth, seqid in ((10, "B73v5_chr3:60664-70130#LTR/Copia/SIRE"), (2, "B73v5_chr2:60664-70130#LTR/Copia/SIRE")):
        (d / f"X_depth{depth}_clean_ltr.tsv").write_text(HEADER + "\n" + _row(seq_id=seqid) + "\n")
    els, _ = load_ltrquest(str(d / "X"), "G")
    assert [e.seqid for e in els] == ["B73v5_chr2", "B73v5_chr3"]


def test_extra_columns_and_ltr_lengths(tmp_path):
    """k2p_se, n_sites, n_ts and n_tv are read; the LTRs' own lengths come from ltr5_len and ltr3_len, not from
    their extents (a nested insertion lies inside the right LTR here)."""
    els, _ = load_ltrquest(str(tmp_path / "X"), "B73")
    host = next(e for e in els if e.start == 51712)
    assert (host.ltr_left_len, host.ltr_right_len) == (1294, 1254) and host.ltr_right == (59462, 70182)
    assert (host.k_se, host.n_sites, host.n_ts, host.n_tv, host.k_kind) == (0.00170359, 1180.0, 4.0, 0.0, "K2P")


def test_ltr_offsets_count_along_the_span_one_based(tmp_path):
    """ltr5_start..ltr3_end count from 1 at the span's left end (its seq_id start) in genome orientation, whatever
    the strand; the element's LTRs are genome coordinates, 1-based closed."""
    (tmp_path / "X_depth0_clean_ltr.tsv").write_text(
        HEADER + "\n" + _row(seq_id="c:100-299#LTR/Gypsy/Tekay", ltr5_start="1", ltr5_end="10", ltr3_start="191",
                             ltr3_end="200") + "\n")
    (tmp_path / "X_depth1_clean_ltr.tsv").unlink()
    (e,), _ = load_ltrquest(str(tmp_path / "X"), "G")
    assert (e.start, e.end, e.ltr_left, e.ltr_right) == (100, 299, (100, 109), (290, 299))


def test_class_rules(tmp_path):
    seq = "c:5-99#LTR/Gypsy"
    (tmp_path / "X_depth0_clean_ltr.tsv").write_text(HEADER + "\n" + _row(seq_id=seq) + "\n")
    (tmp_path / "X_depth1_clean_ltr.tsv").unlink()
    (e,), _ = load_ltrquest(str(tmp_path / "X"), "G")
    assert (e.superfamily, e.clade) == ("Gypsy", None)
    (e,), _ = load_ltrquest(str(tmp_path / "X"), "G", class_rule="last")
    assert (e.superfamily, e.clade) == (None, "Gypsy")


def test_tsd_rules(tmp_path):
    (tmp_path / "X_depth0_clean_ltr.tsv").write_text(HEADER + "\n" + _row(tsd="acgtacgtacg") + "\n")
    (tmp_path / "X_depth1_clean_ltr.tsv").unlink()
    (e,), _ = load_ltrquest(str(tmp_path / "X"), "G")
    assert e.tsd == "ACGTACGTACG"
    (e,), _ = load_ltrquest(str(tmp_path / "X"), "G", tsd_rule="clean")
    assert e.tsd is None


def test_table_class_names_keep_everything_but_unknown(tmp_path):
    """LTRpav's table loader maps only "unknown" to none, and keeps blank-padded or NA TSDs as written but for
    case."""
    (tmp_path / "X_depth0_clean_ltr.tsv").write_text(HEADER + "\n" + _row(seq_id="c:5-99#LTR/NA/unknown") + "\n")
    (tmp_path / "X_depth1_clean_ltr.tsv").unlink()
    (e,), _ = load_ltrquest(str(tmp_path / "X"), "G")
    assert (e.superfamily, e.clade) == ("NA", None)


def test_unknown_rules_raise(tmp_path):
    for kw in ({"class_rule": "x"}, {"tsd_rule": "x"}, {"table_rule": "x"}):
        with pytest.raises(ValueError, match="must be"):
            load_ltrquest(str(tmp_path / "X"), "G", **kw)


def test_phasing_table_rule(tmp_path):
    """Phasing's K2P table (table_rule "phasing"): column names in any case, every row whatever its status, k2p read
    as a number (else unknown), n_sites from ltr5_len when the table has no n_sites column, superfamily and clade
    left unknown, seq_id read as a record name; the LTR columns count along the span (as LTRquest's) and are
    optional: a row without them has no LTR extents."""
    t = tmp_path / "k.tsv"
    t.write_text("#Seq_ID\tStatus\tK2P\tltr5_start\tltr5_end\tltr3_start\tltr3_end\tltr5_len\tltr3_len\tstrand\t"
                 "family\ttsd\n"
                 "c1:100-1099#LTR/Gypsy/Tekay\tfail\t0.01\t1\t50\t901\t1000\t50\t100\t-\tf1\tacgtn\n"
                 "# a comment line\n"
                 "\n"
                 "c1:2001..2500\tpass\tn/a\t.\t.\t.\t.\tNA\t.\t+\tNA\t.\n"
                 "elt::c2:9-60\tpass\t0.02\n")
    els, skipped = load_table(t, "G", tsd_rule="acgtn", table_rule="phasing")
    assert not skipped
    a, b, c = els
    assert (a.seqid, a.start, a.end, a.strand, a.k, a.k_kind) == ("c1", 100, 1099, "-", 0.01, "K2P")
    assert (a.ltr_left, a.ltr_right) == ((100, 149), (1000, 1099))
    assert (a.superfamily, a.clade, a.family, a.tsd, a.tsd_source) == (None, None, "f1", "ACGTN", "annot")
    assert (a.n_sites, a.ltr_left_len, a.ltr_right_len, a.name) == (50.0, 50, 100, "c1:100-1099#LTR/Gypsy/Tekay")
    assert (b.start, b.end, b.k, b.ltr_left, b.ltr_right, b.family, b.n_sites) == (2001, 2500, None, NO_LTR, NO_LTR,
                                                                                   None, None)
    assert (c.seqid, c.start, c.end, c.k, c.strand) == ("c2", 10, 60, 0.02, ".")
    t.write_text("seq_id\tk2p\tn_sites\tltr5_len\nc1:1-50\t0.1\tNA\t40\nnonsense\t0.1\tNA\t40\n")
    (e,), skipped = load_table(t, "G", table_rule="phasing")
    assert e.n_sites is None and skipped == {"bad_seq_id": 1}         # an n_sites column, even NA, is used


def test_phasing_table_needs_k2p_and_reads_any_first_column(tmp_path):
    """Without a seq_id column the first column names the record; a table without a k2p column fails naming it."""
    t = tmp_path / "k.tsv"
    t.write_text("name\tk2p\nc1:1-50\t0.1\n")
    (e,), _ = load_table(t, "G", table_rule="phasing")
    assert (e.seqid, e.start, e.end, e.k) == ("c1", 1, 50, 0.1)
    t.write_text("#seq_id\tdivergence\nc1:1-50\t0.1\n")
    with pytest.raises(ValueError, match=r"k\.tsv: no 'k2p' column"):
        load_table(t, "G", table_rule="phasing")


def test_one_table_file_under_the_ltrquest_rule(tmp_path):
    """A single LTRquest table read as a file is read as each of a prefix's tables."""
    els, skipped = load_table(tmp_path / "X_depth0_clean_ltr.tsv", "B73")
    assert [e.start for e in els] == [60664] and skipped == {"status_not_pass": 1}


@pytest.mark.parametrize("rule", ["ltrquest", "phasing"])
def test_gzip_table_reads_as_the_plain_one(tmp_path, rule):
    """A gzip table reads under either table rule as the same table uncompressed, named .gz or not (its first bytes
    tell gzip)."""
    plain = tmp_path / "X_depth0_clean_ltr.tsv"
    want = load_table(plain, "B73", table_rule=rule)
    assert len(want[0]) == {"ltrquest": 1, "phasing": 2}[rule]
    for name in ("table.tsv.gz", "table.tsv"):
        (tmp_path / name).write_bytes(gzip.compress(plain.read_bytes()))
        assert load_table(tmp_path / name, "B73", table_rule=rule) == want, name


def test_table_with_required_columns_only_and_short_rows(tmp_path):
    """Optional columns the table lacks, or a row that stops before them, read as unknown; a repeated TSD and a
    repeated classification read the same each time."""
    d = tmp_path / "m"
    d.mkdir()
    cols = "seq_id\tstatus\tltr5_start\tltr5_end\tltr3_start\tltr3_end"
    rows = ["c1:100-900#LTR/Gypsy/Tekay\tpass\t1\t50\t751\t801\tacgta",
            "c1:1000-1800#LTR/Gypsy/Tekay\tpass\t1\t50\t751\t801\tacgta",
            "c1:2000-2800#LTR/Copia/Ale\tpass\t1\t50\t751\t801"]
    (d / "Y_depth0_clean_ltr.tsv").write_text(cols + "\ttsd\n" + "\n".join(rows) + "\n")
    (d / "Z_depth0_clean_ltr.tsv").write_text(cols + "\n" + "\n".join(r[:r.rindex("\t")] for r in rows[:2]) + "\n")
    els, skipped = load_ltrquest(str(d / "Y"), "G")
    assert not skipped
    assert [(e.tsd, e.superfamily, e.clade, e.k, e.strand, e.family) for e in els] == [
        ("ACGTA", "Gypsy", "Tekay", None, ".", None), ("ACGTA", "Gypsy", "Tekay", None, ".", None),
        (None, "Copia", "Ale", None, ".", None)]
    els, _ = load_ltrquest(str(d / "Z"), "G")
    assert [(e.tsd, e.tsd_source, e.clade) for e in els] == [(None, None, "Tekay")] * 2
