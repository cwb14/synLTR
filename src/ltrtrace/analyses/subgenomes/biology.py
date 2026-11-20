"""Derived biology of the phased library: the nesting check (an inner element is younger than its host), insertion
histories per subgenome and clade (observed and deconvolved by a g-model of LTR ages), the clade table and the copies
inserted after a dated merger."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ltrtrace.analyses.subgenomes.inputs import _f, group_label, log, warn
from ltrtrace.analyses.subgenomes.merger import AGE_GRID, _age_lik, age_prior, ltr_counts
from ltrtrace.base.constants import HIST_K2P_GRID, MIN_HIST, MU_DEFAULT, PRIOR_POOL_MAX, SE_MULT
from ltrtrace.seq.divergence import k2p_se_copy

HIST_EDGES = np.concatenate(([0.0], np.geomspace(*HIST_K2P_GRID)))   # LTR K2P bins of insertion_history.tsv: 40
#                      log-spaced over 1e-4..0.5 (~11 per decade) plus one for identical LTRs


def clade_names(lib):
    """Annotated LTR-RT clade of each copy (Tekay, CRM, Athila, ...): the GFF clade, else the third field of the FASTA
    header class (LTR/Gypsy/Tekay); 'unknown' when neither names one."""
    out = []
    for c, k in zip(lib.clade, getattr(lib, "cls", [""] * len(lib.clade))):
        if not c:
            f = k.split("/") if k else []
            c = f[2] if len(f) > 2 and f[2] else ""
        out.append(c if c else "unknown")
    return np.array(out)


def merger_cut(r):
    """LTR K2P separating copies inserted before and after the merger of a subgenome pair: the merger estimate; None
    when the merger is only bounded (a bound limits the onset of cross-subgenome transposition, not the merger) or not
    datable."""
    mg = r["merger"]
    if mg is None or not mg["datable"] or mg["bound"]:
        return None
    return float(mg["estimate"])


def after_merger(lib, loc, dt, rng):
    """Per subgenome with a dated merger (the earliest one when several pairs include it): each copy's probability of
    having been inserted after it, from its LTR counts and the subgenome's g-model age distribution (a hard K2P cut
    would count pre-merger copies whose short LTRs carry no substitution) -> {subgenome: array over copies, NaN
    elsewhere}."""
    cut = {}
    for r in (dt["pairs"] if dt else []):
        T = merger_cut(r)
        if T is not None:
            for g in r["groups"]:
                cut[g] = max(cut.get(g, 0.0), T)
    kc, nc = ltr_counts(lib)
    out = {}
    for g, T in sorted(cut.items()):
        idx = np.flatnonzero((loc == g) & np.isfinite(kc) & (nc > 0))
        if idx.size < MIN_HIST:
            continue
        pool = rng.choice(idx, min(idx.size, PRIOR_POOL_MAX), replace=False)
        post = _age_lik(kc[idx], nc[idx]) * age_prior(kc[pool], nc[pool])[None]
        p = np.full(lib.n, np.nan)
        p[idx] = post[:, AGE_GRID < T].sum(1) / np.maximum(post.sum(1), 1e-300)
        out[g] = p
    return out


def nesting_check(lib):
    """Elements inside another element (host = the smallest element strictly containing it on the same sequence): the
    inner one inserted later, so its LTR K2P should not exceed its host's. Violation: inner - host > 2 se (each from
    its own LTR sites; pairs without LTR site counts or lengths are not checked) -> list of dicts (inner, host, K2Ps,
    se, violation, inside a host LTR)."""
    se = k2p_se_copy(lib.k2p, lib.k2p_se, lib.n_sites, lib.ltr5_len)
    has_se = (np.isfinite(lib.k2p_se) & (lib.k2p_se > 0)) | (np.nan_to_num(lib.n_sites) > 0) | \
        (np.nan_to_num(lib.ltr5_len) > 0)                    # K2P alone carries no standard error
    out = []
    order = np.lexsort((-lib.end, lib.start, lib.seq))
    stack = []
    for i in order:
        while stack and (lib.seq[stack[-1]] != lib.seq[i] or lib.end[stack[-1]] < lib.end[i]):
            stack.pop()
        if stack and lib.start[stack[-1]] <= lib.start[i] and (lib.start[stack[-1]], lib.end[stack[-1]]) != (lib.start[i], lib.end[i]):
            h = stack[-1]
            if np.isfinite(lib.k2p[i]) and np.isfinite(lib.k2p[h]) and has_se[i] and has_se[h]:
                s = float(np.hypot(se[i], se[h]))
                l5 = lib.ltr5_len[h] if np.isfinite(lib.ltr5_len[h]) else 0
                l3 = lib.ltr3_len[h] if np.isfinite(lib.ltr3_len[h]) else 0
                out.append(dict(inner=int(i), host=int(h), k_inner=float(lib.k2p[i]), k_host=float(lib.k2p[h]), se=s,
                                violation=bool(lib.k2p[i] - lib.k2p[h] > SE_MULT * s),
                                in_host_ltr=bool(lib.start[i] < lib.start[h] + l5 or lib.end[i] > lib.end[h] - l3)))
        stack.append(i)
    return out


def hist_mass(g):
    """Mass of the age distribution g (point masses on AGE_GRID) in each HIST_EDGES bin: each positive grid point's
    mass spread over its log-age cell, read at the bin edges, so a bin's mass follows its width rather than how many
    grid points fall in it; the age-0 point goes to the first bin, ages above the last edge to the last."""
    a = AGE_GRID[1:]
    r = np.sqrt(a[1] / a[0])
    x = np.log(np.concatenate(([a[0] / r], a * r)))
    C = np.concatenate(([g[0]], g[0] + np.cumsum(g[1:])))
    Ce = np.interp(np.log(HIST_EDGES[1:]), x, C)
    Ce[-1] = C[-1]
    return np.diff(np.concatenate(([0.0], Ce)))


def insertion_history(lib, loc, N, rng):
    """Copies per LTR K2P bin (HIST_EDGES) per subgenome with >= MIN_HIST copies, all clades and each clade with
    >= MIN_HIST copies there (not one holding the whole subgenome): observed, and deconvolved (g-model age distribution
    of age_prior scaled to the copies; removes the binomial noise of short LTRs, not survival or detection biases; NaN
    without LTR substitution counts or lengths)."""
    kc, nc = ltr_counts(lib)
    fin = np.isfinite(lib.k2p)
    cnt = fin & np.isfinite(kc) & (nc > 0)
    clades = clade_names(lib)
    rows = []
    for g in range(N):
        mg = fin & (loc == g)
        cand = ["all"] + sorted(c for c in set(clades[mg]) if MIN_HIST <= (mg & (clades == c)).sum() < mg.sum())
        for c in cand:
            m = mg if c == "all" else mg & (clades == c)
            if m.sum() < MIN_HIST:
                continue
            idx = np.flatnonzero(m & cnt)
            dens = np.full(HIST_EDGES.size - 1, np.nan)
            if idx.size >= MIN_HIST:
                pool = rng.choice(idx, min(idx.size, PRIOR_POOL_MAX), replace=False)
                dens = hist_mass(age_prior(kc[pool], nc[pool])) * idx.size
            obs = np.bincount(np.clip(np.searchsorted(HIST_EDGES, lib.k2p[m], "right") - 1, 0, HIST_EDGES.size - 2),
                              minlength=HIST_EDGES.size - 1)
            for b in range(HIST_EDGES.size - 1):
                rows.append(dict(subgenome=g, clade=c, lo=HIST_EDGES[b], hi=HIST_EDGES[b + 1], observed=int(obs[b]),
                                 deconvolved=float(dens[b])))
    return rows


def clade_table(lib, loc, N, after):
    """Per clade and subgenome: copies, kb, LTR K2P quartiles, and the expected number of copies inserted after the
    subgenome's dated merger (sum of after_merger's probabilities; None without a dated merger)."""
    clades = clade_names(lib)
    L = (lib.end - lib.start + 1).astype(float)
    rows = []
    for g in range(N):
        for c in sorted(set(clades[loc == g])):
            m = (loc == g) & (clades == c)
            k = lib.k2p[m][np.isfinite(lib.k2p[m])]
            q = np.quantile(k, [0.25, 0.5, 0.75]) if k.size else [np.nan] * 3
            rows.append(dict(clade=c, subgenome=g, n=int(m.sum()), kb=L[m].sum() / 1e3, q25=q[0], q50=q[1], q75=q[2],
                             n_after=float(np.nansum(after[g][m])) if g in after else None))
    return rows


def run_biology(a, lib, ph, dt):
    """Derived biology: clade table, insertion history and nesting check (any N)."""
    rng = np.random.default_rng([a.seed, 11])
    N = max(ph["N"], 1)
    ek = ph["ek"]
    loc = np.where(ek >= 0, ph["z"][np.maximum(ek, 0)], -1) if ph["N"] > 1 else np.where(ek >= 0, 0, -1)
    kc, nc = ltr_counts(lib)
    if not (np.isfinite(kc) & (nc > 0)).any():
        warn("no LTR substitution counts or LTR lengths: no deconvolved insertion history, no nesting check, "
             "no copies-after-merger counts")
    nest = nesting_check(lib)
    nv = sum(r["violation"] for r in nest)
    if nest:
        log(f"nesting check: inner element older than its host in {nv} of {len(nest)} nested pairs "
            f"({100 * nv / len(nest):.1f}%; > 2 se)")
    history = insertion_history(lib, loc, N, rng)                  # drawn first (results as before elements.tsv)
    after = after_merger(lib, loc, dt, rng)
    return dict(nesting=nest, history=history, clades=clade_table(lib, loc, N, after), after=after,
                mu=dt["mu"] if dt else (a.mu or MU_DEFAULT), N=N)


def biology_summary(bio):
    """Compact biology record for summary.json: nesting check."""
    nv = sum(r["violation"] for r in bio["nesting"])
    return dict(nesting=dict(n_nested=len(bio["nesting"]), n_inner_older=nv,
                             rate=nv / len(bio["nesting"]) if bio["nesting"] else None))


def write_biology(a, lib, bio):
    od = Path(a.outdir)
    yr = lambda k: _f(k / (2 * bio["mu"]) if k is not None and np.isfinite(k) else None)
    with open(od / "clades.tsv", "w") as f:
        f.write("clade\tsubgenome\tn_copies\tkb\tk2p_q25\tk2p_median\tk2p_q75\tyears_median\tn_after_merger\n")
        for r in bio["clades"]:
            f.write(f"{r['clade']}\t{group_label(r['subgenome'])}\t{r['n']}\t{r['kb']:.1f}\t{_f(r['q25'])}\t"
                    f"{_f(r['q50'])}\t{_f(r['q75'])}\t{yr(r['q50'])}\t{_f(r['n_after'], '{:.1f}')}\n")
    with open(od / "insertion_history.tsv", "w") as f:
        f.write(f"# mu={bio['mu']:g}; observed = copies per LTR K2P bin; deconvolved = smooth g-model fit (LTR "
                f"counting noise removed; survival and detection biases remain)\n")
        f.write("subgenome\tclade\tk2p_low\tk2p_high\tyears_low\tyears_high\tobserved\tdeconvolved\n")
        for r in bio["history"]:
            f.write(f"{group_label(r['subgenome'])}\t{r['clade']}\t{r['lo']:.4g}\t{r['hi']:.4g}\t{yr(r['lo'])}\t"
                    f"{yr(r['hi'])}\t{r['observed']}\t{_f(r['deconvolved'], '{:.2f}')}\n")
    with open(od / "nesting.tsv", "w") as f:
        f.write("inner\thost\tinner_k2p\thost_k2p\tdifference\tse\tinner_older\tin_host_ltr\n")
        for r in bio["nesting"]:
            f.write(f"{lib.ids[r['inner']]}\t{lib.ids[r['host']]}\t{r['k_inner']:.5f}\t{r['k_host']:.5f}\t"
                    f"{r['k_inner'] - r['k_host']:+.5f}\t{r['se']:.5f}\t{'yes' if r['violation'] else 'no'}\t"
                    f"{'yes' if r['in_host_ltr'] else 'no'}\n")
