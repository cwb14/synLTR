"""The figure data (which sequences report.pdf draws) and findings.txt."""
from types import SimpleNamespace

import numpy as np

from synltr import __version__
from synltr.analyses.subgenomes import report
from synltr.base import constants as cst
from synltr.stats.report import Findings


def test_shown_keeps_the_longest_grouped():
    """More phased sequences than FIG_MAX_SEQ: the longest are drawn, grouped by subgenome in natural order."""
    n = cst.FIG_MAX_SEQ + 10
    names = [f"s{i}" for i in range(n)]
    lib = SimpleNamespace(seqnames=names, seq_len=np.arange(1, n + 1) * 1000)
    ph = dict(names=names, z=np.array([i % 2 for i in range(n)]))
    shown, L = report._shown(lib, ph)
    assert len(shown) == cst.FIG_MAX_SEQ and min(L[shown]) == 11000
    zz = ph["z"][shown]
    assert list(zz) == sorted(zz) and [names[i] for i in shown][:2] == ["s10", "s12"]


def _homoeolog_stub():
    """Two subgenomes of two sequences each: three of the four sequences paired, nothing else to say."""
    calls = {(0, 1): dict(status="homoeolog"), (2, 1): dict(status="homoeolog"), (1, 0): dict(status="homoeolog"),
             (3, 0): dict(status="none")}
    ph = dict(z=np.array([0, 1, 0, 1]), names=["a1", "b1", "a2", "b2"], N=2, homoeologs=dict(calls=calls))
    hb = dict(mu=1e-8, reversed=[], segments=[], halflife=[], asym=[], transloc=[], many=[], density={0: 1.0, 1: 1.0})
    return ph, hb


def test_findings_go_through_the_shared_writer(tmp_path, caplog):
    """What the homoeolog pairs say: statements in a Findings, findings.txt in the subgenomes layout (the header
    line, naming synltr and its version, a blank line, one statement per line), every line echoed to the log."""
    ph, hb = _homoeolog_stub()
    f = report.homoeolog_findings(ph, hb)
    assert isinstance(f, Findings)
    said = [text for _, text in f.items()]
    assert said[0] == ("Homoeologs SG1-SG2: 2 of 2 SG1 and 1 of 2 SG2 sequences paired; the others keep too few "
                       "inherited insertions in collinear order (transposon turnover erases them as subgenomes age).")
    assert said[1:] == ["Orientation: every called homoeolog pair runs in the same direction in this assembly.",
                        ("Exchanges: no stretch of near-identical inherited insertions on the called homoeologs "
                         "(replacing exchanges >= 1 Mb are found 9 times in 10 in simulations; smaller or older ones "
                         "can be missed).")]
    with caplog.at_level("INFO", logger="synltr"):
        report.write_homoeolog_biology(SimpleNamespace(outdir=str(tmp_path)), ph, hb)
    want = [(f"What the homoeologous chromosome pairs say about this genome (synltr {__version__}). Each "
             "statement passed a benchmark bar set in advance."), "", *said]
    assert (tmp_path / "findings.txt").read_text() == "\n".join(want) + "\n"
    assert [r.getMessage() for r in caplog.records] == want
