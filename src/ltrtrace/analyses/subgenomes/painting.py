"""Ancestry painting: each copy's lineage ancestry (abundance and breadth of its clades on the other sequences,
calibrated to the copies' locations), an HMM over each phased sequence's informative copies, exchanges called when
enough independent lineages favor another subgenome decisively, the smallest detectable exchange per window, and the
ancestry of sequences too poor in LTR-RTs to phase. bedGraph and BED tiles are 0-based half-open."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from numba import njit

from ltrtrace.analyses.subgenomes.inputs import group_label, log
from ltrtrace.base.constants import (
    AGE_BINS,
    AGE_CONFINED_MIN,
    AGE_WEIGHT_MIN,
    BF_DECISIVE,
    BF_SUPPORT,
    BREADTH_FOLD,
    BREADTH_MIN_SEQS,
    GFF3_SOURCE,
    HALF_COUNT,
    MIN_LINEAGES,
    NMAX_KNN,
    PAINT_FIT_MAX,
    START_OTHER,
    SWITCH_ITERS,
    SWITCH_RANGE,
    SWITCH_START,
    SWITCH_TOL,
)
from ltrtrace.base.constants import PAINT_RES_WIN as RES_WIN
from ltrtrace.base.unionfind import find
from ltrtrace.io.annotation.gff3 import escape_seqid


@njit(cache=True)
def _ancestry_scores(members, indptr, mult, seq, zseq, N, T, Cn, nseq_tot):
    """S_ab, n_ab (abundance), S_br, n_br (breadth) per copy. seq: sequence of each copy (-1 unphased: scored on
    every phased carrier, nothing left out); zseq: group of each sequence; T: phased copies per group; Cn: phased
    sequences per group."""
    n = seq.size
    Sa = np.zeros((n, N))
    na = np.zeros(n)
    Sb = np.zeros((n, N))
    nb = np.zeros(n)
    cnt_seq = np.zeros(nseq_tot, np.int64)
    M = np.zeros(N)
    Mc = np.zeros(N)
    m = np.zeros(N)
    sc = np.zeros(N)
    for c in range(indptr.size - 1):
        lo = indptr[c]
        hi = indptr[c + 1]
        M[:] = 0.0
        Mc[:] = 0.0
        for p in range(lo, hi):
            s = seq[members[p]]
            if s >= 0:
                if cnt_seq[s] == 0:
                    Mc[zseq[s]] += 1.0
                cnt_seq[s] += 1
                M[zseq[s]] += 1.0
        w = mult[c]
        for p in range(lo, hi):
            i = members[p]
            s = seq[i]
            g = zseq[s] if s >= 0 else -1                  # unphased copy: scored on all phased carriers
            # abundance: other members, own sequence left out
            tot = 0.0
            for h in range(N):
                m[h] = M[h]
            if g >= 0:
                m[g] -= cnt_seq[s]
            for h in range(N):
                tot += m[h]
            if tot > 0:
                mean = 0.0
                for h in range(N):
                    sc[h] = math.log((m[h] + HALF_COUNT) / T[h])
                    mean += sc[h]
                mean /= N
                for h in range(N):
                    Sa[i, h] += w * (sc[h] - mean)
                na[i] += w
            # breadth: other sequences per group carrying the clade, own sequence left out
            mx = 0.0
            for h in range(N):
                m[h] = Mc[h]
            if g >= 0:
                m[g] -= 1.0
            for h in range(N):
                if m[h] > mx:
                    mx = m[h]
            if mx >= BREADTH_MIN_SEQS:
                smax = -1e300
                smin = 1e300
                for h in range(N):
                    avail = Cn[h] - (1.0 if h == g else 0.0)
                    sc[h] = math.log((m[h] + HALF_COUNT) / (avail + 1.0))
                    if sc[h] > smax:
                        smax = sc[h]
                    if sc[h] < smin:
                        smin = sc[h]
                if smax - smin >= math.log(BREADTH_FOLD):
                    mean = 0.0
                    for h in range(N):
                        mean += sc[h]
                    mean /= N
                    for h in range(N):
                        Sb[i, h] += w * (sc[h] - mean)
                    nb[i] += w
        for p in range(lo, hi):
            s = seq[members[p]]
            if s >= 0:
                cnt_seq[s] = 0
    return Sa, na, Sb, nb


def _calibrate(S, n, y):
    """P(o | copy) = softmax(lam * S / n^beta + b) fitted to the copies' location labels (foreign copies are rare
    label noise that only makes the fit conservative)."""
    from scipy.optimize import minimize
    from scipy.special import expit, logsumexp
    N = S.shape[1]
    ln = np.log(np.maximum(n, 1.0))
    scale = float(np.median(np.abs(S).sum(1))) + 1e-12

    def unpack(par):
        return math.exp(min(par[0], 50.0)) / scale, float(expit(par[1])), np.concatenate(([0.0], par[2:]))

    def nll(par):
        lam, beta, b = unpack(par)
        Z = lam * S * np.exp(-beta * ln)[:, None] + b[None, :]
        return -(Z[np.arange(len(y)), y] - logsumexp(Z, 1)).sum()
    best = min((minimize(nll, np.concatenate(([0.0, b0], np.zeros(N - 1))), method="L-BFGS-B")
                for b0 in (-2.0, 0.0, 2.0)), key=lambda r: r.fun)
    return unpack(best.x)


def _apply_cal(S, n, cal):
    from scipy.special import logsumexp
    lam, beta, b = cal
    Z = lam * S * np.exp(-beta * np.log(np.maximum(n, 1.0)))[:, None] + b[None, :]
    return Z - logsumexp(Z, 1, keepdims=True)


def copy_ancestry(cl, ek, z, N, rng, fit_max=PAINT_FIT_MAX):
    """Per-copy lineage ancestry on the clades: two experts (lineage abundance and breadth on the other sequences,
    own sequence left out), each calibrated to the location labels, averaged, recalibrated. -> log P (n x N), number
    of informative k-mers per copy."""
    own = np.where(ek >= 0, z[np.maximum(ek, 0)], -1)
    T = np.bincount(own[own >= 0], minlength=N).astype(float)
    Cn = np.bincount(z, minlength=N).astype(float)
    Sa, na, Sb, nb = _ancestry_scores(cl.members.astype(np.int64), cl.indptr.astype(np.int64),
                                      cl.mult.astype(np.float64), np.asarray(ek, np.int64), np.asarray(z, np.int64),
                                      N, T, Cn, int(z.size))
    parts = []
    for S, nn in ((Sa, na), (Sb, nb)):
        fit = np.flatnonzero((own >= 0) & (nn > 0))
        if fit.size > fit_max:
            fit = rng.choice(fit, fit_max, replace=False)
        lp = _apply_cal(S, nn, _calibrate(S[fit], nn[fit], own[fit]))
        lp[nn == 0] = -math.log(N)
        parts.append(lp)
    Z = np.mean(parts, 0)
    Z -= Z.mean(1, keepdims=True)
    ninf = np.maximum(na, nb)
    fit = np.flatnonzero((own >= 0) & (ninf > 0))
    if fit.size > fit_max:
        fit = rng.choice(fit, fit_max, replace=False)
    return _apply_cal(Z, np.ones(ek.size), _calibrate(Z[fit], np.ones(fit.size), own[fit])), ninf


def age_weights(lib, own, lP, ninf, N):
    """Weight in [0, 1] of each copy as evidence about its neighborhood's origin: 1 at LTR ages whose copies are
    most confined to their own subgenome, 0 at ages whose copies carry other-subgenome ancestry as often as
    uninformative copies would (1 - prior of their own subgenome: post-merger insertions, lineages older than the
    split); copies without K2P get 1; all 1 when no age class is confined by >= 0.05 or all ages are equal."""
    use = (own >= 0) & (ninf > 0) & np.isfinite(lib.k2p)
    w = np.ones(lib.n)
    if use.sum() < AGE_WEIGHT_MIN:
        return w
    prior = np.bincount(own[own >= 0], minlength=N) / max((own >= 0).sum(), 1)
    share = 1 - np.exp(lP[np.arange(lib.n), np.maximum(own, 0)])
    s_null = float(np.mean(1 - prior[own[use]]))
    k = lib.k2p[use]
    edges = np.unique(np.quantile(k, np.linspace(0, 1, AGE_BINS + 1)))
    if edges.size < 2:                                          # every copy of one age: no age classes to compare
        return w
    bi = np.clip(np.searchsorted(edges, k, "right") - 1, 0, edges.size - 2)
    cnt = np.bincount(bi, minlength=edges.size - 1)
    ok = cnt > 0
    A = np.array([np.median(k[bi == j]) for j in np.flatnonzero(ok)])
    Y = np.bincount(bi, weights=share[use], minlength=edges.size - 1)[ok] / cnt[ok]
    if s_null - Y.min() < AGE_CONFINED_MIN:
        return w
    w = np.clip((s_null - np.interp(np.nan_to_num(lib.k2p, nan=-1.0), A, Y)) / (s_null - Y.min()), 0.0, 1.0)
    w[~np.isfinite(lib.k2p)] = 1.0
    return w


@njit(cache=True)
def _lse(a):
    m = a.max()
    return m + np.log(np.exp(a - m).sum())


@njit(cache=True)
def hmm_posterior(E, l0, sw):
    """Forward-backward over copies (log emissions E, n x N; log initial l0), switch probability sw per step split
    evenly over the other states -> (posteriors n x N, log-likelihood, expected number of switches)."""
    n, N = E.shape
    lsw = np.log(sw / max(N - 1, 1))
    lst = np.log1p(-sw)
    f = np.empty((n, N))
    b = np.zeros((n, N))
    f[0] = l0 + E[0]
    tmp = np.empty(N)
    for t in range(1, n):
        for j in range(N):
            for i in range(N):
                tmp[i] = f[t - 1, i] + (lst if i == j else lsw)
            f[t, j] = _lse(tmp) + E[t, j]
    for t in range(n - 2, -1, -1):
        for i in range(N):
            for j in range(N):
                tmp[j] = (lst if i == j else lsw) + E[t + 1, j] + b[t + 1, j]
            b[t, i] = _lse(tmp)
    ll = _lse(f[n - 1])
    post = np.exp(f + b - ll)
    esw = 0.0
    for t in range(n - 1):
        for i in range(N):
            for j in range(N):
                if i != j:
                    esw += np.exp(f[t, i] + lsw + E[t + 1, j] + b[t + 1, j] - ll)
    return post, ll, esw


@njit(cache=True)
def hmm_viterbi(E, l0, sw):
    """Most likely state path (same model as hmm_posterior)."""
    n, N = E.shape
    lsw = np.log(sw / max(N - 1, 1))
    lst = np.log1p(-sw)
    v = l0 + E[0]
    bp = np.zeros((n, N), np.int64)
    nv = np.empty(N)
    for t in range(1, n):
        for j in range(N):
            best = -1e300
            arg = 0
            for i in range(N):
                s = v[i] + (lst if i == j else lsw)
                if s > best:
                    best = s
                    arg = i
            nv[j] = best + E[t, j]
            bp[t, j] = arg
        v[:] = nv
    path = np.empty(n, np.int64)
    path[n - 1] = np.argmax(v)
    for t in range(n - 1, 0, -1):
        path[t - 1] = bp[t, path[t]]
    return path


def switch_rate(seqs, N):
    """One genome-wide switch probability per copy step by Baum-Welch (emissions fixed)."""
    sw = SWITCH_START
    for _ in range(SWITCH_ITERS):
        es = tot = 0.0
        for E, l0 in seqs:
            es += hmm_posterior(E, l0, sw)[2]
            tot += E.shape[0] - 1
        new = float(np.clip(es / max(tot, 1.0), *SWITCH_RANGE))
        done = abs(new - sw) <= SWITCH_TOL * sw
        sw = new
        if done:
            break
    return sw


class Lineages:
    """Copies sharing a k-mer clade belong to one lineage, one piece of evidence: through small clades
    (<= NMAX_KNN carriers) anywhere, and through any clade whose carriers lie mostly (>= half) inside the segment being
    judged (a lineage amplified in one place has no small clade: its k-mers are private or family-wide)."""

    def __init__(self, cl, n):
        cid = np.repeat(np.arange(cl.n), cl.size)
        mem = cl.members.astype(np.int64)
        o = np.argsort(mem, kind="stable")
        self.mem, self.cid, self.size = mem[o], cid[o], cl.size
        self.ptr = np.searchsorted(self.mem, np.arange(n + 1))

    def labels(self, members):
        """Lineage label of each copy in members: connected components through qualifying shared clades."""
        members = np.asarray(members)
        parent = np.arange(members.size)
        hits = {}
        for p, c in enumerate(members):
            for k in self.cid[self.ptr[c]:self.ptr[c + 1]]:
                hits.setdefault(int(k), []).append(p)
        for k, ps in hits.items():
            if len(ps) < 2 or (self.size[k] > NMAX_KNN and 2 * len(ps) < self.size[k]):
                continue
            r0 = find(parent, ps[0])
            for p in ps[1:]:
                r = find(parent, p)
                if r != r0:
                    parent[r] = r0
        return np.array([find(parent, i) for i in range(members.size)], np.int64)


def resolution_track(mid, n_needed, L):
    """Smallest exchange detectable around each RES_WIN window of a sequence: the shortest span of n_needed consecutive
    informative copies (midpoints, sorted here) that covers the window's center -> (starts, ends, bp); no rows for
    windows centred outside the copies, none at all when the sequence holds fewer copies than needed."""
    mid = np.sort(np.asarray(mid, float))
    starts = np.arange(0, L, RES_WIN, dtype=np.int64)
    ends = np.minimum(starts + RES_WIN, L)
    if mid.size < n_needed:
        return starts[:0], ends[:0], np.zeros(0)
    span = mid[n_needed - 1:] - mid[:mid.size - n_needed + 1]            # copies j .. j + n_needed - 1
    x = (starts + ends) / 2
    j0 = np.maximum(np.searchsorted(mid, x, "left") - n_needed + 1, 0)   # last copy at or after the center
    j1 = np.minimum(np.searchsorted(mid, x, "right") - 1, span.size - 1)  # first copy at or before it
    keep = (x >= mid[0]) & (x <= mid[-1]) & (j1 >= j0)
    out = np.array([span[a:b + 1].min() for a, b in zip(j0[keep], j1[keep])])
    return starts[keep], ends[keep], out


def posterior_tiles(mid, post, L, R):
    """Ancestry posterior tiles of one sequence for the bedGraphs: each informative copy (midpoint order) owns the
    stretch to the midpoints between it and its neighbours, but at most R/2 on either side; longer stretches without
    copies stay unknown (no tile) -> (starts, ends, posteriors)."""
    o = np.argsort(mid, kind="stable")
    m, p = np.asarray(mid, float)[o], post[o]
    half = R / 2
    lo = np.empty(m.size)
    hi = np.empty(m.size)
    if m.size:
        between = (m[1:] + m[:-1]) / 2
        lo[0], hi[-1] = max(0.0, m[0] - half), min(float(L), m[-1] + half)
        lo[1:] = np.maximum(between, m[1:] - half)
        hi[:-1] = np.minimum(between, m[:-1] + half)
    lo, hi = np.floor(lo).astype(np.int64), np.floor(hi).astype(np.int64)
    keep = hi > lo
    return lo[keep], hi[keep], p[keep]


def run_painting(a, lib, ph, dt):
    """LTR-only ancestry along each phased sequence: HMM over informative copies, emissions = calibrated per-copy
    ancestry (log P - log prior) x age weight, one genome-wide switch rate (Baum-Welch); Viterbi segments of another
    subgenome are called exchanges when their likelihood ratio, each lineage counted once, is >= BF_DECISIVE and >=
    MIN_LINEAGES lineages favour the origin >= BF_SUPPORT-fold. Resolution: smallest detectable exchange per window
    (the span of the informative copies whose typical evidence pays the two switches and the call threshold)."""
    if ph["N"] < 2:
        return None
    N, z, ek = ph["N"], ph["z"], ph["ek"]
    ancestry = dt.get("ancestry") if dt else None
    if ancestry is None:
        ancestry = copy_ancestry(ph["clades"], ek, z, N, np.random.default_rng([a.seed, 7]))
    lP, ninf = ancestry
    own = np.where(ek >= 0, z[np.maximum(ek, 0)], -1)
    prior = np.bincount(own[own >= 0], minlength=N) / max((own >= 0).sum(), 1)
    w = age_weights(lib, own, lP, ninf, N)
    inf = (ninf > 0) & (ek >= 0)
    si = {s: i for i, s in enumerate(lib.seqnames)}
    seqs, meta = [], []
    for i, s in enumerate(ph["names"]):
        idx = np.flatnonzero(inf & (ek == i))
        if idx.size < 2:
            continue
        idx = idx[np.argsort(lib.start[idx], kind="stable")]
        E = np.ascontiguousarray((lP[idx] - np.log(prior)[None]) * w[idx][:, None])
        l0 = np.log(np.where(np.arange(N) == z[i], 1 - START_OTHER, START_OTHER / (N - 1)))
        seqs.append((E, l0))
        meta.append((i, s, idx))
    sw = switch_rate(seqs, N) if seqs else SWITCH_START
    lin = Lineages(ph["clades"], lib.n)
    home = [E[np.arange(E.shape[0]), int(z[i])] - np.max(np.delete(E, int(z[i]), 1), 1)
            for (E, _), (i, _, _) in zip(seqs, meta)]                  # evidence of copies for their own subgenome
    pos = np.concatenate([h[h > 0] for h in home]) if home else np.zeros(0)
    # a segment enters the best path only when its evidence pays both switches (in and out) and reaches the call
    # threshold: copies needed = that log-evidence / the median evidence of a copy for its own subgenome
    need = max(math.log(BF_DECISIVE), 2 * (math.log1p(-sw) - math.log(sw / (N - 1))))
    n_needed = int(max(MIN_LINEAGES, math.ceil(need / np.median(pos)))) if pos.size else MIN_LINEAGES
    segs, post_tracks, res = [], [], []
    for (E, l0), (i, s, idx) in zip(seqs, meta):
        g0 = int(z[i])
        post, _, _ = hmm_posterior(E, l0, sw)
        path = hmm_viterbi(E, l0, sw)
        mid = (lib.start[idx] + lib.end[idx]) // 2
        L = int(lib.seq_len[si[s]])
        post_tracks.append((s, mid, post, L))
        res.append((s,) + resolution_track(mid.astype(float), n_needed, L))
        brk = np.flatnonzero(np.diff(path)) + 1
        for s0, e0 in zip(np.concatenate(([0], brk)), np.concatenate((brk, [idx.size]))):
            g = int(path[s0])
            if g == g0:
                continue
            ii = idx[s0:e0]
            d = E[s0:e0, g] - E[s0:e0, g0]
            lab = lin.labels(ii)
            _, inv = np.unique(lab, return_inverse=True)
            llr = float((np.bincount(inv, weights=d) / np.bincount(inv)).sum())
            sup = d >= math.log(BF_SUPPORT)
            nlin = int(np.unique(lab[sup]).size)
            segs.append(dict(seq=s, start=int(lib.start[ii].min()), end=int(lib.end[ii].max()), own=g0, origin=g,
                             left_flank=int(lib.end[idx[s0 - 1]]) if s0 > 0 else 1,
                             right_flank=int(lib.start[idx[e0]]) if e0 < idx.size else L,
                             n_ltr=int(ii.size), n_support=int(sup.sum()), n_against=int((d <= -math.log(BF_SUPPORT)).sum()),
                             n_lineages=nlin, llr=llr, post=float(post[s0:e0, g].mean()),
                             called=bool(llr >= math.log(BF_DECISIVE) and nlin >= MIN_LINEAGES)))
    nc = sum(x["called"] for x in segs)
    unph = unphased_ancestry(a, lib, ph, prior, N)
    rr = np.concatenate([r[3] for r in res]) if res else np.zeros(0)
    log(f"painting: switch {sw:.2e} per copy; {nc} exchange segments called (likelihood ratio >= {BF_DECISIVE:g} with "
        f"each lineage counted once, >= {MIN_LINEAGES} lineages); smallest detectable exchange: median "
        f"{np.median(rr) / 1e3:,.0f} kb ({n_needed} informative copies)" if rr.size else f"painting: {nc} exchange segments called")
    return dict(segs=segs, switch=sw, n_needed=n_needed, post=post_tracks, res=res, N=N, ancestry=ancestry,
                unphased=unph)


def unphased_ancestry(a, lib, ph, prior, N):
    """Ancestry of each sequence too short or too poor in LTR-RTs to phase, from its informative copies: the
    copies' evidence (log posterior over the prior) summed over lineages, each lineage counted once; called when the
    best subgenome beats the next >= BF_DECISIVE-fold with >= MIN_LINEAGES lineages favoring it >= BF_SUPPORT-fold
    ('decisive'), or >= BF_SUPPORT-fold ('substantial'); otherwise none. Phasing's clades hold phased copies only:
    the unphased copies' ancestry and lineages come from all clades (phased results are untouched)."""
    phased = set(ph["names"])
    ek = np.asarray(ph["ek"])
    if not (ek < 0).any() or "clades_all" not in ph:
        return []
    lP, ninf = copy_ancestry(ph["clades_all"], ek, np.asarray(ph["z"]), N, np.random.default_rng([a.seed, 13]))
    lin = Lineages(ph["clades_all"], lib.n)
    E = lP - np.log(prior)[None]
    o = np.argsort(lib.seq, kind="stable")
    cuts = np.flatnonzero(np.diff(lib.seq[o])) + 1
    out = []
    for grp in np.split(o, cuts):
        if not grp.size or lib.seqnames[lib.seq[grp[0]]] in phased:
            continue
        idx = grp[ninf[grp] > 0]
        row = dict(seq=lib.seqnames[lib.seq[grp[0]]], n_ltr=int(grp.size), n_inf=int(idx.size), n_lin=0, g=-1, fav=-1,
                   lo=0.0, conf="none")
        if idx.size:
            lab = lin.labels(idx)
            _, inv = np.unique(lab, return_inverse=True)
            cnt = np.bincount(inv)
            tot = np.array([(np.bincount(inv, weights=E[idx, g]) / cnt).sum() for g in range(N)])
            best = np.argsort(-tot)
            g, lo = int(best[0]), float(tot[best[0]] - tot[best[1]])
            d = E[idx, g] - np.max(np.delete(E[idx], g, 1), 1)
            nlin = int(np.unique(lab[d >= math.log(BF_SUPPORT)]).size)
            conf = ("decisive" if lo >= math.log(BF_DECISIVE) and nlin >= MIN_LINEAGES else
                    "substantial" if lo >= math.log(BF_SUPPORT) else "none")
            row.update(n_lin=nlin, g=g if conf != "none" else -1, fav=g, lo=lo, conf=conf)
        out.append(row)
    return out


def painting_summary(pa):
    """Compact painting record for summary.json: switch rate, calls, exchanged bp per origin -> recipient."""
    bp = {}
    for x in pa["segs"]:
        if x["called"]:
            k = f"{group_label(x['origin'])}>{group_label(x['own'])}"
            bp[k] = bp.get(k, 0) + x["end"] - x["start"] + 1
    rr = np.concatenate([r[3] for r in pa["res"]]) if pa["res"] else np.zeros(0)
    return dict(switch=pa["switch"], n_segments=len(pa["segs"]), n_called=sum(x["called"] for x in pa["segs"]),
                n_unphased=len(pa.get("unphased", [])),
                n_unphased_called=sum(u["conf"] != "none" for u in pa.get("unphased", [])),
                exchanged_bp=bp, copies_needed=pa["n_needed"],
                resolution_median_bp=float(np.median(rr)) if rr.size else None)


def write_unphased(od, rows):
    """unphased_sequences.tsv: ancestry (called subgenome or NA), favored (the subgenome the odds favor, NA without
    informative copies), log10 odds of the favored over the next subgenome (NA without informative copies)."""
    with open(Path(od) / "unphased_sequences.tsv", "w") as f:
        f.write("sequence\tn_ltr\tn_informative\tn_lineages\tancestry\tfavored\tlog10_odds_vs_next\tconfidence\n")
        for u in rows:
            ev = u["n_inf"] > 0
            f.write(f"{u['seq']}\t{u['n_ltr']}\t{u['n_inf']}\t{u['n_lin']}\t"
                    f"{group_label(u['g']) if u['g'] >= 0 else 'NA'}\t{group_label(u['fav']) if ev else 'NA'}\t"
                    f"{u['lo'] / math.log(10):.2f}\t{u['conf']}\n" if ev else
                    f"{u['seq']}\t{u['n_ltr']}\t0\t0\tNA\tNA\tNA\t{u['conf']}\n")


def write_painting(a, pa):
    od = Path(a.outdir)
    if pa.get("unphased"):
        write_unphased(od, pa["unphased"])
    with open(od / "segments.tsv", "w") as f:
        f.write("sequence\tstart\tend\tlength\tleft_flank\tright_flank\tsubgenome\torigin\tn_ltr\tn_support\t"
                "n_against\tn_lineages\tlog10_lr_lineages\tmean_posterior\tcalled\n")
        for x in pa["segs"]:
            f.write(f"{x['seq']}\t{x['start']}\t{x['end']}\t{x['end'] - x['start'] + 1}\t{x['left_flank']}\t"
                    f"{x['right_flank']}\t{group_label(x['own'])}\t"
                    f"{group_label(x['origin'])}\t{x['n_ltr']}\t{x['n_support']}\t{x['n_against']}\t{x['n_lineages']}\t"
                    f"{x['llr'] / math.log(10):.2f}\t{x['post']:.3f}\t{'yes' if x['called'] else 'no'}\n")
    called = [x for x in pa["segs"] if x["called"]]
    with open(od / "exchanges.bed", "w") as f:
        for k, x in enumerate(called):
            f.write(f"{x['seq']}\t{x['start'] - 1}\t{x['end']}\texchange_{k + 1}:{group_label(x['origin'])}>"
                    f"{group_label(x['own'])}\t{min(1000, int(100 * x['llr'] / math.log(10)))}\t.\n")
    with open(od / "exchanges.gff3", "w") as f:
        f.write("##gff-version 3\n")
        for k, x in enumerate(called):
            f.write(f"{escape_seqid(x['seq'])}\t{GFF3_SOURCE}\tregion\t{x['start']}\t{x['end']}\t"
                    f"{x['llr'] / math.log(10):.2f}\t.\t.\tID=exchange_{k + 1};"
                    f"origin={group_label(x['origin'])};recipient={group_label(x['own'])};"
                    f"lineages={x['n_lineages']};copies={x['n_ltr']}\n")
    rr = np.concatenate([r[3] for r in pa["res"]]) if pa["res"] else np.zeros(0)
    R = float(np.median(rr)) if rr.size else float(RES_WIN)             # longer copy-free stretches stay unknown
    for g in range(pa["N"]):
        with open(od / f"ancestry_{group_label(g)}.bedGraph", "w") as f:
            f.write(f"track type=bedGraph name=\"P({group_label(g)})\" description=\"posterior ancestry {group_label(g)}"
                    f" per informative LTR-RT (tiles between copy midpoints; copy-free stretches longer than the median"
                    f" resolution, {R / 1e3:,.0f} kb, unknown)\"\n")
            for s, mid, post, L in pa["post"]:
                for p_, q_, v_ in zip(*posterior_tiles(mid, post, L, R)):
                    f.write(f"{s}\t{p_}\t{q_}\t{v_[g]:.4f}\n")
    with open(od / "resolution.bedGraph", "w") as f:
        f.write(f"track type=bedGraph name=\"resolution\" description=\"smallest detectable exchange (bp): span of "
                f"{pa['n_needed']} informative LTR-RTs\"\n")
        for s, st, en, v in pa["res"]:
            for p, q, r in zip(st, en, v):
                f.write(f"{s}\t{p}\t{q}\t{int(r)}\n")


def _favors(lib, ph, anc):
    """Each copy's evidence as the painting reads it (log posterior over the subgenome's share of copies): informative
    copies; those whose lineages favor another subgenome >= BF_SUPPORT-fold over their sequence's own (foreign), and
    their own >= BF_SUPPORT-fold over every other (home); the most favored other subgenome."""
    lP, ninf = anc
    N, ek, z = ph["N"], np.asarray(ph["ek"]), np.asarray(ph["z"])
    own = np.where(ek >= 0, z[np.maximum(ek, 0)], 0)
    prior = np.bincount(own[ek >= 0], minlength=N) / max((ek >= 0).sum(), 1)
    E = lP - np.log(np.maximum(prior, 1e-300))[None]
    r = np.arange(lib.n)
    Eo = E[r, own]
    Ex = E.copy()
    Ex[r, own] = -np.inf
    gbest, d = Ex.argmax(1), Ex.max(1) - Eo
    inf = (ninf > 0) & (ek >= 0)
    return inf, inf & (d >= math.log(BF_SUPPORT)), inf & (-d >= math.log(BF_SUPPORT)), gbest
