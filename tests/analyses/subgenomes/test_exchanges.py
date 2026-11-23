"""Near-identical inherited pairs on homoeologs (exchanges) and their runs."""
import numpy as np

from synltr.analyses.subgenomes import exchanges


def test_flag_homogenized_and_runs_on_planted_values():
    """Near-identical pairs among inherited ones are flagged (empirical null), a run of them is significant and
    scattered flags are not."""
    rng = np.random.default_rng(5)
    k = np.r_[np.exp(rng.normal(np.log(0.05), 0.3, 300)), np.full(12, 0.0005)]
    f, thr = exchanges.flag_homogenized(k, np.full(k.size, 600.0))
    assert f[300:].all() and f[:300].sum() <= 2 and 0.0005 < thr < 0.05
    f, _ = exchanges.flag_homogenized(k, np.full(k.size, 600.0), limit=0.001)    # merger limit flags them too
    assert f[300:].all()
    chain = np.zeros(120, bool)
    chain[40:48] = True
    chain[[5, 90]] = True
    runs = exchanges.runs_of_flags(chain)
    assert [(a, b) for a, b, _ in runs] == [(40, 47)] and runs[0][2] < 1e-4
    scattered = np.zeros(120, bool)
    scattered[::12] = True
    assert exchanges.runs_of_flags(scattered) == []
