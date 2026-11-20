"""What the homoeolog pairs say: near-identical stretches of inherited insertions on called homoeologs (exchanges after
the split), the half-life of intact LTR-RTs in each subgenome, translocated blocks, sequences with several homoeologs
and homoeologs collinear in reverse. Each output passed a bar fixed before it was scored (simulations with exact truth,
real polyploids with gene-synteny and exchange truth)."""
from __future__ import annotations

import math
from collections import defaultdict
from statistics import NormalDist

import numpy as np

from ltrtrace.analyses.subgenomes.inputs import log, natural_key
from ltrtrace.base.constants import (
    CI95,
    HALF_COUNT,
    HB_ALPHA,
    HB_MIN_RUN,
    HL_BOOT,
    HL_BOOT_MIN,
    HL_CHECK_FACTOR,
    HL_FACTOR,
    HL_MIN_R,
    HL_PLATEAU,
    HOMOG_ITERS,
    HOMOG_MIN,
    HOMOG_SCALE_FLOOR,
    MIN_ANCHORS,
    MU_DEFAULT,
    SG_LTR_MIN_SITES,
    SPLIT_MIN,
)
from ltrtrace.seq import align, divergence
from ltrtrace.seq.chain import anchor_ranks, chain_test, longest_chain


def flag_homogenized(k2p, sites, limit=0.0, alpha=HB_ALPHA, it=HOMOG_ITERS):
    """Anchors too similar to be inherited from the progenitors' split -> (flags, threshold K2P). Empirical null on
    log divergence (LTR rates vary, so a binomial is too narrow): centre = median, scale = spread of the upper half
    (homogenized pairs cannot reach it); flagged below centre - z(alpha / n) x scale, or at/below the merger's upper
    limit (K2P; 0 = unknown); iterated on the unflagged. NaN divergences are never flagged."""
    k = np.nan_to_num(np.asarray(k2p, float), nan=-1.0)
    x = np.log(np.maximum(np.where(k >= 0, k, np.nan), HALF_COUNT / np.maximum(np.asarray(sites, float), 1.0)))
    ok = np.isfinite(x)
    n = int(ok.sum())
    if n < HOMOG_MIN:
        return np.zeros(x.size, bool), float("nan")
    zq = NormalDist().inv_cdf(1 - alpha / n)
    keep = ok & (k > limit)
    if keep.sum() < HOMOG_MIN:
        keep = ok.copy()
    thr = float("nan")
    for _ in range(it):
        centre = float(np.median(x[keep]))
        up = x[keep][x[keep] >= centre]
        # spread: the upper half's median distance over the normal's (0.6745, the z of the 75th percentile)
        thr = centre - zq * max(float(np.median(up) - centre) / 0.6745, HOMOG_SCALE_FLOOR)
        new = ok & (x > thr) & (k > limit)
        if new.sum() < HOMOG_MIN or np.array_equal(new, keep):
            break
        keep = new
    return ok & ~((x > thr) & (k > limit)), float(math.exp(thr))


def runs_of_flags(flags, min_run=HB_MIN_RUN):
    """Runs of >= min_run consecutive flags in chain order -> [(first, last, p)], p = chance of a run at least this
    long anywhere in the chain given the number of flags (Poisson approximation for runs)."""
    f = np.asarray(flags, bool)
    n = f.size
    q = f.sum() / n if n else 0.0
    out, i = [], 0
    while i < n:
        if not f[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and f[j + 1]:
            j += 1
        L = j - i + 1
        if L >= min_run:
            lam = (n - L + 1) * (1 - q) * q ** L + q ** L
            out.append((i, j, float(1 - math.exp(-lam))))
        i = j + 1
    return out


def _merger_k2p(dt, g, h):
    """Upper limit of the merger of groups g and h (K2P; 0 when the dating could not date or bound it)."""
    for r in (dt or {}).get("pairs", []):
        if set(int(x) for x in r["groups"]) == {g, h}:
            mg = r.get("merger")
            return float(mg["hi"]) if mg and mg.get("datable") and mg.get("hi") is not None else 0.0
    return 0.0


def run_homoeolog_biology(a, lib, ph, dt, pa):
    """Divergence of every shared old insertion; near-identical stretches on called
    homoeologs (exchanges after the split); half-life of intact LTR-RTs in each direction and the two lineages
    compared; blocks collinear with a non-homoeolog (translocations) and sequences paired with several homoeologs;
    homoeologs collinear in reverse. None without homoeolog calls. Reads ph, dt and pa; changes nothing."""
    ho = ph.get("homoeologs")
    if not ho or not len(ho["P"]) or not any(c["partner"] >= 0 for c in ho["calls"].values()):
        return None
    N, C, names = ph["N"], ph["C"], ph["names"]
    z, ek = np.asarray(ph["z"]), np.asarray(ph["ek"])
    P, si, sj = ho["P"], np.asarray(ho["si"]), np.asarray(ho["sj"])
    rng = np.random.default_rng([a.seed, 91])                  # own stream: no earlier draw changes
    cnt = align.pairs(lib.store, P, mode="ltr", threads=a.threads)
    k2p, _ = divergence.k2p(cnt.sites, cnt.ts, cnt.tv)
    k2p = np.where(cnt.sites >= SG_LTR_MIN_SITES, k2p, np.nan)
    hot = ho["in_chain"] & ho["pair_called"]
    flags = np.zeros(len(P), bool)
    mu = dt["mu"] if dt else (a.mu or MU_DEFAULT)
    key = np.minimum(si, sj) * C + np.maximum(si, sj)
    gp = np.minimum(z[si], z[sj]) * N + np.maximum(z[si], z[sj])
    pairs = {}
    for g in range(N):
        for h in range(g + 1, N):
            idx = np.flatnonzero(hot & (gp == g * N + h))
            lim = _merger_k2p(dt, g, h)
            f, thr = flag_homogenized(k2p[idx], cnt.sites[idx], lim)
            flags[idx] = f
            use = idx[np.isfinite(k2p[idx]) & ~f]
            pairs[(g, h)] = dict(limit=lim, threshold=thr, use=use, n_flagged=int(f.sum()),
                                 d_s=float(np.median(k2p[use])) if use.size >= SPLIT_MIN else float("nan"))
    called = sorted({(min(i, c["partner"]), max(i, c["partner"])) for (i, h), c in ho["calls"].items()
                     if c["partner"] >= 0})
    v8 = [x for x in (pa or {}).get("segs", []) if x["called"]]
    segs = []
    for i, j in called:
        q = np.flatnonzero(hot & (key == i * C + j))
        ci, cj = np.where(si[q] == i, P[q, 0], P[q, 1]), np.where(si[q] == i, P[q, 1], P[q, 0])
        o = np.argsort(lib.start[ci], kind="stable")
        fin = np.isfinite(k2p[q[o]])
        q, ci, cj = q[o][fin], ci[o][fin], cj[o][fin]
        for s0, s1, p in runs_of_flags(flags[q]):
            if p > HB_ALPHA / len(called):
                continue
            r = slice(s0, s1 + 1)
            seg = dict(i=i, j=j, start_i=int(lib.start[ci[r]].min()), end_i=int(lib.end[ci[r]].max()),
                       start_j=int(lib.start[cj[r]].min()), end_j=int(lib.end[cj[r]].max()), n=s1 - s0 + 1,
                       n_chain=int(q.size), k2p=float(np.median(k2p[q[r]])), p=p)
            lim = pairs[(min(z[i], z[j]), max(z[i], z[j]))]["limit"]
            seg["timing"] = ("after_merger" if seg["k2p"] <= lim else "before_merger") if lim > 0 else "undetermined"
            seg["v8"] = [f"exchange_{k + 1}" for k, x in enumerate(v8)
                         if (x["seq"] == names[i] and x["start"] <= seg["end_i"] and seg["start_i"] <= x["end"]) or
                         (x["seq"] == names[j] and x["start"] <= seg["end_j"] and seg["start_j"] <= x["end"])]
            segs.append(seg)
    # half-life per direction g -> h (h's loss rate, given g's copy survived)
    partner = {(i, h): c["partner"] for (i, h), c in ho["calls"].items() if c["partner"] >= 0}
    sidx = {s_: k_ for k_, s_ in enumerate(lib.seqnames)}
    Lph = np.array([lib.seq_len[sidx[n_]] for n_ in names], float)
    ncop = np.bincount(ek[ek >= 0], minlength=C)
    density = {g: float(ncop[z == g].sum() / max(Lph[z == g].sum() / 1e6, 1e-9)) for g in range(N)}
    fin_k = np.isfinite(lib.k2p)
    by_seq = defaultdict(list)
    for c_ in np.flatnonzero((ek >= 0) & fin_k):
        by_seq[int(ek[c_])].append(int(c_))
    hl, edge = [], float(ho["old_edge"])
    for (g0, h0), pr in pairs.items():
        ds = pr["d_s"]
        if not np.isfinite(ds):
            continue
        thr = max(pr["limit"], pr["threshold"] if np.isfinite(pr["threshold"]) else 0.0)
        inh = np.flatnonzero((gp == g0 * N + h0) & np.isfinite(k2p) & (k2p > thr))
        links = defaultdict(lambda: defaultdict(int))
        for q in inh:
            links[int(P[q, 0])][int(sj[q])] += 1
            links[int(P[q, 1])][int(si[q])] += 1
        chain_k = k2p[pr["use"]]
        for g, h in ((g0, h0), (h0, g0)):
            chroms = [i for i in range(C) if z[i] == g and (i, h) in partner]
            nB = int((z == h).sum())
            if not chroms or nB < 2:
                continue
            per = {}
            for i in chroms:
                j = partner[(i, h)]
                cs = np.array(by_seq[i], np.int64)
                ks = lib.k2p[cs]
                pp = np.array([links[int(c_)].get(j, 0) for c_ in cs], float)
                qq = np.array([sum(n_ for t_, n_ in links[int(c_)].items() if z[t_] == h and t_ != j) for c_ in cs],
                              float)
                o = np.argsort(-ks)
                per[i] = (-ks[o], np.r_[0, np.cumsum(pp[o])], np.r_[0, np.cumsum(qq[o])])

            def ret(t, chs):
                obs = paired = other = 0.0
                for i in chs:
                    nk, cp, cq = per[i]
                    m = int(np.searchsorted(nk, -t, side="left"))
                    obs, paired, other = obs + m, paired + cp[m], other + cq[m]
                return (paired - other / max(nB - 1, 1)) / max(obs, 1.0), obs
            r15, n15 = ret(max(edge, HL_FACTOR * ds), chroms)
            r2, _ = ret(max(edge, HL_CHECK_FACTOR * ds), chroms)
            est = math.log(2) * ds / -math.log(r15) if 0 < r15 < 1 else float("nan")
            boot = np.full(HL_BOOT, np.nan)
            for b in range(HL_BOOT):
                chs = list(rng.choice(chroms, len(chroms), replace=True))
                dsb = float(np.median(rng.choice(chain_k, chain_k.size, replace=True)))
                rb = ret(max(edge, HL_FACTOR * dsb), chs)[0]
                boot[b] = math.log(2) * dsb / -math.log(rb) if 0 < rb < 1 else np.nan
            good = np.isfinite(boot).sum() > HL_BOOT_MIN
            lo, up = np.nanquantile(boot, CI95) if good else (np.nan, np.nan)
            qc = dict(plateau=r15 > 0 and abs(r2 - r15) <= HL_PLATEAU * r15, signal=bool(good), range=r15 >= HL_MIN_R)
            hl.append(dict(g=g, h=h, d_s=ds, r15=float(r15), r2=float(r2), n_old=int(n15), est=est, lo=float(lo),
                           hi=float(up), ok=all(qc.values()), qc=qc, boot=boot))
    asym = []
    for e1 in hl:
        e2 = next((e for e in hl if e["g"] == e1["h"] and e["h"] == e1["g"]), None)
        if e2 is None or e1["g"] > e1["h"] or not (e1["ok"] and e2["ok"]):
            continue
        m = np.isfinite(e1["boot"]) & np.isfinite(e2["boot"]) & (e1["boot"] > 0) & (e2["boot"] > 0)
        lr = np.log(e1["boot"][m] / e2["boot"][m])
        lo, up = np.quantile(lr, CI95) if lr.size else (np.nan, np.nan)
        asym.append(dict(a=e1["h"], b=e2["h"], ratio=e1["est"] / e2["est"], lo=float(np.exp(lo)), hi=float(np.exp(up)),
                         differ=bool(lr.size) and not (lo <= 0 <= up)))
    # translocated blocks: a collinear chain between i and a chromosome other than its homoeolog
    pset = defaultdict(set)
    for (i, h), j in partner.items():
        pset[i].add(j)
        pset[j].add(i)
    called_set = set(called)
    o = np.argsort(key, kind="stable")
    st = np.r_[0, np.flatnonzero(np.diff(key[o])) + 1] if o.size else np.zeros(0, np.int64)
    tests = []
    for s0, e0 in zip(st, np.r_[st[1:], o.size]):
        q = o[s0:e0]
        i, j = int(key[q[0]] // C), int(key[q[0]] % C)
        if q.size < MIN_ANCHORS or (i, j) in called_set or not (pset[i] or pset[j]):
            continue
        ci, cj = np.where(si[q] == i, P[q, 0], P[q, 1]), np.where(si[q] == i, P[q, 1], P[q, 0])
        oi, x, y = anchor_ranks(lib.start[ci], lib.start[cj])
        ci, cj = ci[oi], cj[oi]
        tests.append((i, j, ci, cj, x, y, longest_chain(x, y)))
    tl, alpha = [], HB_ALPHA / max(len(tests), 1)
    for i, j, ci, cj, x, y, ch in tests:
        if ch.size < MIN_ANCHORS:
            continue
        p = chain_test(x, y, int(ch.size), alpha, rng)
        if p <= alpha:
            tl.append(dict(i=i, j=j, start_i=int(lib.start[ci[ch]].min()), end_i=int(lib.end[ci[ch]].max()),
                           start_j=int(lib.start[cj[ch]].min()), end_j=int(lib.end[cj[ch]].max()), n=int(ch.size),
                           n_tested=int(x.size), p=float(p)))
    many = []
    for s_, ps in sorted(pset.items()):
        for h in sorted({int(z[p_]) for p_ in ps}):
            mem = sorted((p_ for p_ in ps if z[p_] == h), key=lambda p_: natural_key(names[p_]))
            if len(mem) >= 2:
                many.append(dict(i=s_, h=h, partners=mem))
    rev = [(i, j) for i, j in called if ho["pair"].get((i, j), {}).get("orientation") == "-"]
    log(f"homoeolog biology: {len(segs)} near-identical stretches, {sum(e['ok'] for e in hl)} of {len(hl)} half-lives "
        f"estimable, {len(tl)} translocated blocks, {len(many)} sequences with several homoeologs, {len(rev)} "
        f"homoeolog pairs collinear in reverse")
    return dict(k2p=k2p, sites=cnt.sites, flags=flags, pairs=pairs, segments=segs, halflife=hl, asym=asym,
                density=density, transloc=tl, many=many, reversed=rev, mu=mu, n_called=len(called))
