from dataclasses import replace

import numpy as np
import pytest

from synltr.io.elements import (
    ELEMENT_COLUMNS,
    Element,
    ElementTable,
    elements_digest,
    make_eid,
    read_elements,
    write_elements,
)


def _a():
    return Element("G", make_eid("G", "chr1", 10, 200), "n1", "chr1", 10, 200, "+", (10, 30), (180, 200), "ACGTA",
                   0.1 + 0.2, "f", "Gypsy", "Tekay", "high", parent="G|chr1:1-999", level=1)


def _b():
    return Element("G", make_eid("G", "chr2", 5, 99), "n2", "chr2", 5, 99, "-", (5, 14), (90, 99), None, None,
                   None, None, None, None)


def test_element_defaults_and_zstart():
    a = _a()
    assert a.zstart == 9 and a.start == 10  # 0-based start of the 1-based closed interval
    b = _b()
    assert (b.ltr_left_len, b.ltr_right_len, b.tsd_source) == (None, None, None)  # unknown, never assumed
    assert (b.k_kind, b.k_se, b.n_sites, b.n_ts, b.n_tv) == ("K2P", None, None, None, None)
    with pytest.raises(AttributeError):
        a.start = 3  # frozen


def test_table_round_trips_rows_and_columns():
    t = ElementTable.from_rows([_a(), _b()])
    assert len(t) == 2 and list(t) == [_a(), _b()] and t.row(1) == _b()
    assert t.column("start").tolist() == [10, 5]
    assert t.column("ltr_left").tolist() == [[10, 30], [5, 14]]
    assert np.isnan(t.column("k")[1]) and t.column("k")[0] == 0.1 + 0.2
    assert t.column("seqid").tolist() == ["chr1", "chr2"]
    assert len(ElementTable.from_rows([])) == 0


def test_concat_keeps_rows_in_order():
    t = ElementTable.concat([ElementTable.from_rows(rows) for rows in ([_a()], [], [_b()])])
    assert list(t) == [_a(), _b()] and t.column("ltr_left").shape == (2, 2) and t.column("start").dtype == np.int64


def test_table_select_by_mask_and_index():
    t = ElementTable.from_rows([_a(), _b()])
    assert list(t.select(t.column("start") > 8)) == [_a()]
    assert list(t.select(np.array([1, 0]))) == [_b(), _a()]


def test_table_keeps_the_extra_fields():
    e = Element("G", "G|c:1-9", "x", "c", 1, 9, "+", (1, 3), (7, 9), "ACGTA", 0.2, None, None, None, None,
                ltr_left_len=2, ltr_right_len=3, tsd_source="inferred", k_kind="JC69", k_se=0.01, n_sites=300.0,
                n_ts=4.0, n_tv=1.0)
    assert ElementTable.from_rows([e]).row(0) == e


def test_digest_matches_the_recorded_value():
    """The digest of the two elements above, as LTRpav computed it (checkpoint keys stay valid)."""
    assert elements_digest([_a(), _b()]) == "7f72c775f9146d466b74006bd30a319b75cda03a3b706b11c6fa1972358ceece"
    assert elements_digest(ElementTable.from_rows([_a(), _b()])) == elements_digest([_a(), _b()])
    assert elements_digest(None) == "none"
    assert elements_digest([]) == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"


def test_elements_tsv_is_ltrpavs_columns_exactly(tmp_path):
    p = tmp_path / "e.tsv"
    e = Element("g", "g|c:1-10", "c:1-10", "c", 1, 10, ".", (1, 2), (9, 10), None, 0.0123456789012345, None, None,
                None, None, k_se=0.5, n_sites=10.0, k_kind="JC69")
    write_elements(ElementTable.from_rows([_a(), _b(), e]), p)
    lines = p.read_text().splitlines()
    assert lines[0] == "\t".join(ELEMENT_COLUMNS) and len(ELEMENT_COLUMNS) == 17
    assert lines[1] == ("G\tG|chr1:10-200\tn1\tchr1\t10\t200\t+\t10-30\t180-200\tACGTA\t0.30000000000000004\tf\tGypsy\t"
                        "Tekay\thigh\tG|chr1:1-999\t1")
    assert lines[2] == "G\tG|chr2:5-99\tn2\tchr2\t5\t99\t-\t5-14\t90-99\t.\t.\t.\t.\t.\t.\t.\t0"
    got = read_elements(p)
    assert isinstance(got, ElementTable)
    assert list(got)[:2] == [replace(_a(), k_kind=None), replace(_b(), k_kind=None)]  # kind and source unknown
    assert got.row(0).tsd_source is None and got.row(0).k_kind is None
    assert got.row(2).k == 0.0123456789012345


def test_read_rejects_a_foreign_header(tmp_path):
    p = tmp_path / "e.tsv"
    p.write_text("a\tb\n")
    with pytest.raises(ValueError, match="unexpected header"):
        read_elements(p)
