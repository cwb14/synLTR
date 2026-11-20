"""The run's tables (chromosomes, structure, pairing map, links, summary.json, elements, homoeologs and what they say,
findings.txt) and report_data.json: the figure data report.pdf is drawn from (plot.subgenomes). Element coordinates are
1-based closed, as in the library's record names."""
from __future__ import annotations

import json
import math
from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from ltrtrace import __version__ as VERSION
from ltrtrace.analyses.subgenomes.biology import biology_summary, clade_names
from ltrtrace.analyses.subgenomes.dating import PROF_EDGES, cross_profile, dating_summary
from ltrtrace.analyses.subgenomes.inputs import _f, group_label, log, natural_key
from ltrtrace.analyses.subgenomes.links import expected
from ltrtrace.analyses.subgenomes.painting import _favors, copy_ancestry, painting_summary, posterior_tiles
from ltrtrace.base.constants import (
    FIG_MAX_SEQ,
    MIN_ELEMENTS,
    MU_DEFAULT,
    PROF_MIN_PAIRS,
    PROGRAM,
    REPORT_DATA,
    SUPPORT_OK,
    Z95,
    K,
)
from ltrtrace.base.constants import PAINT_RES_WIN as RES_WIN
from ltrtrace.stats.report import Findings, write_findings

# the run objects' fields report.pdf reads (figure_data keeps these and drops the rest)
_PH_KEYS = ("names", "z", "ek", "W", "halves", "support", "tree", "source", "verdict", "strata_edges")
# findings.txt in the subgenomes layout (stats.report.write_findings without sections): this header line (it names the
# package and its version), a blank line, then one statement per line; every statement sits in the one section SECTION
FINDINGS_HEADER = (f"What the homoeologous chromosome pairs say about this genome ({PROGRAM}). Each statement passed a "
                   "benchmark bar set in advance.")
SECTION = "Homoeologs"


def _ma(t):
    """Readable age from years."""
    return f"{t / 1e6:.2f} Myr" if t >= 1e6 else f"{t / 1e3:,.0f} kyr" if t >= 1e3 else f"{t:,.0f} years"


def _homoeolog_numbers(ph):
    """Number of each sequence's homoeolog set (1..K, sets ordered by their first member in their lowest
    subgenome); empty without homoeolog calls."""
    ho = ph.get("homoeologs")
    if not ho or not ho["sets"]:
        return {}
    z = np.asarray(ph["z"])
    rank = {i: r for r, i in enumerate(sorted(range(len(ph["names"])), key=lambda i: natural_key(ph["names"][i])))}

    def first(s_):
        g0 = min(int(z[x]) for x in s_["members"])
        return g0, min(rank[x] for x in s_["members"] if z[x] == g0)
    return {x: k + 1 for k, s_ in enumerate(sorted(ho["sets"], key=first)) for x in s_["members"]}


def _shown(lib, ph):
    """Phased sequences the figures draw, grouped by subgenome: all of them, or the FIG_MAX_SEQ longest. Within a
    subgenome, members of homoeolog sets come first in set order (homoeologs share a position in every group), then
    the rest in natural order -> (indices into ph['names'], lengths in bp of every phased sequence)."""
    si = {s: i for i, s in enumerate(lib.seqnames)}
    L = np.array([lib.seq_len[si[s]] for s in ph["names"]], float)
    keep = np.argsort(-L, kind="stable")[:FIG_MAX_SEQ]
    z = np.asarray(ph["z"])
    num = _homoeolog_numbers(ph)
    return sorted(keep.tolist(), key=lambda i: (z[i], num.get(i, 1 << 30), natural_key(ph["names"][i]))), L


def _oe(W):
    """Observed over expected lineage sharing (degree-corrected null), NaN on the diagonal and without expectation."""
    W = np.asarray(W, float)
    E = expected(W)
    with np.errstate(divide="ignore", invalid="ignore"):
        OE = np.where(E > 0, W / E, np.nan)
    np.fill_diagonal(OE, np.nan)
    return OE


def _wilson(k, n, zq=Z95):
    """95% Wilson interval of a binomial share."""
    n = np.maximum(np.asarray(n, float), 1.0)
    p = np.asarray(k, float) / n
    den = 1 + zq ** 2 / n
    c = (p + zq ** 2 / (2 * n)) / den
    h = zq * np.sqrt(p * (1 - p) / n + zq ** 2 / (4 * n ** 2)) / den
    return c - h, c + h


def _fig_data(a, lib, ph, dt, bio, pa):
    """Everything the figure draws, from the run's objects: sequences shown, painting runs, lineage sharing,
    closest-relative profiles, insertion histories and the range of the time axis (years)."""
    N, names, z, ek = ph["N"], list(ph["names"]), np.asarray(ph["z"]), np.asarray(ph["ek"])
    shown, L = _shown(lib, ph)
    mu = dt["mu"] if dt else (a.mu or MU_DEFAULT)
    P = dict(N=N, names=names, z=z, shown=shown, L=L, mu=mu, mu_given=bool(dt["mu_given"] if dt else a.mu), lib=lib,
             dotplot=N > 1 and bool(ph.get("homoeologs")),             # panel B shows the homoeolog evidence
             ph=ph, dt=dt, bio=bio or {}, pa=pa, loc=np.where(ek >= 0, z[np.maximum(ek, 0)], -1),
             si={s: i for i, s in enumerate(lib.seqnames)}, label=a.label or Path(a.ltr_fasta[0]).name.split(".")[0],
             unit="chromosomes" if np.median(L) >= 5e6 else "sequences")
    # painting: posterior tiles merged into runs of one painted subgenome and confidence level; copy-free stretches
    # longer than the median resolution, or than ~2 pixel rows at this scale, stay unknown (bedGraph rule)
    tiles, R = {}, float(RES_WIN)
    if pa and pa.get("post"):
        rr = np.concatenate([r[3] for r in pa["res"]]) if pa.get("res") else np.zeros(0)
        R = float(np.median(rr)) if rr.size else float(RES_WIN)
        R_draw = max(R, float(L[shown].max()) / 400)
        for s, mid, post, Ls in pa["post"]:
            lo, hi, p = posterior_tiles(mid, post, Ls, R_draw)
            if not lo.size:
                continue
            c, lev = p.argmax(1), np.digitize(p.max(1), [0.7, 0.9, 0.99])
            key = c * 4 + lev
            brk = np.flatnonzero((np.diff(key) != 0) | (lo[1:] != hi[:-1])) + 1
            st, en = np.r_[0, brk], np.r_[brk, lo.size]
            tiles[s] = (lo[st], hi[en - 1], c[st], lev[st])
    P.update(tiles=tiles, R=R, calls=[c for c in (pa["segs"] if pa else []) if c["called"]])
    W_src = dict(ph.get("W_src") or {"all": ph["W_raw"]})
    P["oe"] = {s: _oe(W_)[np.ix_(shown, shown)] for s, W_ in W_src.items()}
    prof = []                                          # share of closest-relative pairs in two subgenomes, by divergence
    if dt and dt.get("pt") is not None:
        Pp, d = dt["pt"]["P"], dt["pt"]["d"]
        la, lb = P["loc"][Pp[:, 0]], P["loc"][Pp[:, 1]]
        cen = 10 ** ((PROF_EDGES[:-1] + PROF_EDGES[1:]) / 2)
        for r in dt["pairs"]:
            x, y = r["groups"]
            m0 = np.isin(la, [x, y]) & np.isin(lb, [x, y])
            if not m0.any():
                continue
            n, k = cross_profile(d[m0], (la[m0] != lb[m0]).astype(float))
            ok = n >= PROF_MIN_PAIRS
            if ok.sum() >= 2:
                lo, hi = _wilson(k[ok], n[ok])
                prof.append(dict(pair=(x, y), t=cen[ok] / (2 * mu), n=n[ok], f=k[ok] / n[ok], lo=lo, hi=hi))
    P["prof"] = prof
    rows = defaultdict(list)                           # insertion history per (subgenome, clade)
    for r in P["bio"].get("history", []):
        rows[(int(r["subgenome"]), str(r["clade"]))].append((r["lo"], r["hi"], r["observed"], r["deconvolved"]))
    P["hist"] = {}
    for key, rr_ in rows.items():
        lo, hi, ob, de = (np.array(v, float) for v in zip(*rr_))
        t = np.where(lo > 0, np.sqrt(np.maximum(lo, 1e-300) * hi), hi / 2) / (2 * mu)
        P["hist"][key] = dict(t=t, y=de if np.isfinite(de).all() else ob, n=float(ob.sum()))
    old = [p_["t"].max() for p_ in prof]
    for r in (dt["pairs"] if dt else []):
        old += [v["hi"] / (2 * mu) for v in r["split"].values() if v and np.isfinite(v.get("hi", np.nan))]
    fin = lib.k2p[np.isfinite(lib.k2p)]
    old.append(float(np.quantile(fin, 0.995)) / (2 * mu) if fin.size else 1e6)
    t0 = 0.4e-4 / (2 * mu)                             # below the identical-LTR bin of the insertion histories
    P["trange"] = (t0, max(1.25 * max(old), 100 * t0))
    return P


def write_elements(a, lib, ph, dt, bio, pa):
    """elements.tsv: every LTR-RT with its sequence's subgenome (NA when unphased), clade, family, LTR divergence and
    age, the calibrated ancestry posterior per subgenome for informative copies, a call by the painting's evidence rule
    (foreign: its lineages favor another subgenome >= BF_SUPPORT-fold over its sequence's own; own: the reverse; weak:
    neither; NA for one subgenome or uninformative copies) with the favored other subgenome, and the probability that
    it was inserted after its subgenome's dated merger (NA without a dated merger)."""
    N, ek, z = ph["N"], np.asarray(ph["ek"]), np.asarray(ph["z"])
    sg = np.where(ek >= 0, z[np.maximum(ek, 0)], -1)
    mu = dt["mu"] if dt else (a.mu or MU_DEFAULT)
    clades = clade_names(lib)
    after = bio.get("after", {})
    if N > 1:
        anc = pa["ancestry"] if pa else (dt["ancestry"] if dt else None)
        if anc is None:
            anc = copy_ancestry(ph["clades"], ek, z, N, np.random.default_rng([a.seed, 7]))
        inf, foreign, home, gbest = _favors(lib, ph, anc)
        P = np.exp(anc[0])
    lines = ["element\tsequence\tstart\tend\tstrand\tclade\tfamily\tk2p\tage_years\tsubgenome\tinformative\t"
             + "\t".join(f"P_{group_label(g)}" for g in range(max(N, 1))) + "\tcall\tfavored\tp_after_merger"]
    for i in range(lib.n):
        k = lib.k2p[i]
        if N > 1 and sg[i] >= 0 and inf[i]:
            probs = "\t".join(f"{v:.4f}" for v in P[i])
            call = "foreign" if foreign[i] else "own" if home[i] else "weak"
            fav = group_label(int(gbest[i])) if foreign[i] else "NA"
            info = "yes"
        else:
            probs, call, fav = "\t".join("NA" for _ in range(max(N, 1))), "NA", "NA"
            info = "no" if N > 1 and sg[i] >= 0 else "NA"
        pa_ = after[sg[i]][i] if sg[i] in after else float("nan")
        lines.append(f"{lib.ids[i]}\t{lib.seqnames[lib.seq[i]]}\t{lib.start[i]}\t{lib.end[i]}\t{lib.strand[i] or '.'}\t"
                     f"{clades[i]}\t{lib.family[i] or 'NA'}\t{_f(k, '{:.5f}')}\t"
                     f"{_f(k / (2 * mu) if np.isfinite(k) else None, '{:.0f}')}\t"
                     f"{group_label(int(sg[i])) if sg[i] >= 0 else 'NA'}\t{info}\t{probs}\t{call}\t{fav}\t"
                     f"{_f(pa_, '{:.3f}')}")
    (Path(a.outdir) / "elements.tsv").write_text("\n".join(lines) + "\n")


def _jsonable(x):
    if isinstance(x, dict):
        return {k: _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        return float(x) if math.isfinite(x) else None
    return x


def write_phasing(a, lib, cl, ph, timings, dt=None, bio=None, pa=None):
    od = Path(a.outdir)
    n_el_all = np.bincount(lib.seq, minlength=len(lib.seqnames))
    idx = {s: i for i, s in enumerate(ph["names"])}
    partner, psup, pex = ph["partners"]
    with open(od / "chromosomes.tsv", "w") as f:
        f.write("sequence\tlength\tn_ltr\tn_informative\tsubgenome\tstatus\tsupport\toe_own\toe_other\tpvalue\t"
                "partner\tpartner_support\tpartner_excess\tcore\n")
        for s, L, n in zip(lib.seqnames, lib.seq_len, n_el_all):
            i = idx.get(s)
            if i is None:
                f.write(f"{s}\t{L}\t{n}\tNA\tNA\tunphased" + "\tNA" * 8 + "\n")
                continue
            if ph["source"] == "user":
                st = "user"
            elif ph["N"] == 1:
                st = "no_structure" if ph["verdict"] == "no structure detected" else "cannot_tell"
            elif ph["n_inf"][i] < MIN_ELEMENTS:
                st = "insufficient"
            else:
                st = "assigned" if ph["support"][i] >= SUPPORT_OK else "ambiguous"
            pj = partner[i]
            f.write(f"{s}\t{L}\t{n}\t{ph['n_inf'][i]}\t{group_label(ph['z'][i]) if ph['N'] > 1 else 'SG1'}\t{st}\t"
                    f"{_f(ph['support'][i], '{:.3f}')}\t{_f(ph['own_oe'][i])}\t{_f(ph['oth_oe'][i])}\t"
                    f"{_f(ph['pval'][i], '{:.3g}')}\t{ph['names'][pj] if pj >= 0 else 'NA'}\t{_f(psup[i], '{:.2f}')}\t"
                    f"{_f(pex[i], '{:.2f}')}\t{'yes' if ph.get('core', np.ones(ph['C'], bool))[i] else 'no'}\n")
    with open(od / "structure.tsv", "w") as f:
        src = [k[len("replicability_"):] for k in (ph["tree"][0] if ph["tree"] else {})
               if k.startswith("replicability_") and k != "replicability_se"]
        f.write("level\tn\tn_a\tn_b\tpart_a\tpart_b\tcopies\treplicability\treplicability_se\tn_pairs\tborderline\t"
                "strength\tparent_strength\t" + "\t".join(f"replicability_{x}\tstrength_{x}\tagrees_{x}\t"
                                                            f"parent_strength_{x}" for x in src)
                + "\trule\taccepted\tforced\n")
        yn = lambda v: "NA" if v is None else ("yes" if v else "no")
        for t in ph["tree"]:
            f.write(f"{t['level']}\t{t['n']}\t{len(t['part_a'])}\t{len(t['part_b'])}\t{','.join(t['part_a'])}\t"
                    f"{','.join(t['part_b'])}\t{t['copies'] or 'NA'}\t{_f(t['replicability'], '{:.3f}')}\t"
                    f"{_f(t['replicability_se'], '{:.3f}')}\t{t['n_pairs']}\t{yn(t['borderline'])}\t"
                    f"{_f(t['strength'])}\t{_f(t['parent_strength'])}\t"
                    + "\t".join(f"{_f(t[f'replicability_{x}'], '{:.3f}')}\t{_f(t[f'strength_{x}'])}\t"
                                 f"{yn(t[f'agrees_{x}'])}\t{_f(t[f'parent_strength_{x}'])}"
                                 for x in src)
                    + f"\t{t['rule']}\t{yn(t['accepted'])}\t{yn(t['forced'])}\n")
    with open(od / "pairing_map.tsv", "w") as f:
        f.write("sequence_a\tsequence_b\tn_candidates\tchain\tnull_mean\tz\tsignificant\tspan_a\tspan_b\t"
                "shared_share\tanchor_k2p_median\tcopy_k2p_median\tpaired\n")
        for r in sorted(ph["pairing"], key=lambda r: -r["z"]):
            f.write(f"{ph['names'][r['i']]}\t{ph['names'][r['j']]}\t{r['n_candidates']}\t{r['chain']}\t"
                    f"{r['null_mean']:.2f}\t{r['z']:.2f}\t{'yes' if r['significant'] else 'no'}\t"
                    f"{r['span_a']:.3f}\t{r['span_b']:.3f}\t{_f(r['shared'], '{:.2f}')}\t"
                    f"{_f(r['d_median'], '{:.5f}')}\t{_f(r['k_median'], '{:.5f}')}\t{'yes' if r['paired'] else 'no'}\n")
    Wr = ph["W_raw"]
    E = expected(Wr)
    with open(od / "links_oe.tsv", "w") as f:                    # observed / expected links, NA on the diagonal
        f.write("sequence\t" + "\t".join(ph["names"]) + "\n")
        for i, s in enumerate(ph["names"]):
            f.write(s + "\t" + "\t".join("NA" if i == j or E[i, j] <= 0 else f"{Wr[i, j] / E[i, j]:.4g}"
                                         for j in range(len(ph["names"]))) + "\n")
    sup = ph["support"][np.isfinite(ph["support"])]
    summ = dict(version=VERSION, ltr_fasta=[str(p) for p in a.ltr_fasta], annotation=[str(p) for p in a.annotation],
                k=K, scale=cl.scale, n_ltr=lib.n, n_sequences=len(lib.seqnames), n_phased=ph["C"],
                n_core=int(np.asarray(ph.get("core", np.ones(ph["C"], bool))).sum()),
                n_kmers=cl.n_occ, n_multi_kmers=cl.n_multi, n_clades=cl.n, n_clades_phased=ph["n_clades"],
                n_subgenomes=ph["N"], n_subgenomes_requested=ph["N_requested"], n_subgenomes_source=ph["source"],
                user_labels=ph.get("user_labels"), verdict=ph["verdict"],
                min_support=float(sup.min()) if sup.size else None, n_half_pairs=ph["n_pairs"],
                subgenome_sizes={group_label(g): int((ph["z"] == g).sum()) for g in range(ph["N"])},
                pairing=dict(n_tested=len(ph["pairing"]), n_significant=sum(r["significant"] for r in ph["pairing"]),
                             n_paired=sum(r["paired"] for r in ph["pairing"]),
                             cannot_link_used=[[ph["names"][x], ph["names"][y]] for x, y in ph["cl_pairs"]],
                             paired_groups=[[ph["names"][x] for x in g] for g in ph["paired_groups"]]),
                structure_tree=ph["tree"], strata_k2p_edges=ph["strata_edges"], strata_copies=ph["strata_n"],
                homoeologs=homoeolog_summary(ph), timings=timings,
                dating=dating_summary(dt) if dt else None, biology=biology_summary(bio) if bio else None,
                painting=painting_summary(pa) if pa else None)
    (od / "summary.json").write_text(json.dumps(_jsonable(summ), indent=1) + "\n")


def homoeolog_summary(ph):
    """summary.json 'homoeologs': evidence used and, per subgenome pair, sequences paired in each direction."""
    ho = ph.get("homoeologs")
    if not ho:
        return None
    z, N, names = np.asarray(ph["z"]), ph["N"], ph["names"]
    calls, pairs = ho["calls"], []
    for g in range(N):
        for h in range(g + 1, N):
            ig, ih = np.flatnonzero(z == g), np.flatnonzero(z == h)
            stat = lambda rows, k: {s_: sum(calls[(i, k)]["status"] == s_ for i in rows)
                                    for s_ in ("homoeolog", "unresolved", "none")}
            pairs.append(dict(subgenomes=[group_label(g), group_label(h)], sequences=[int(ig.size), int(ih.size)],
                              status=[stat(ig, h), stat(ih, g)],
                              reciprocal=int(sum(bool(calls[(i, h)]["reciprocal"]) for i in ig))))
    hb, bio = ho.get("bio"), None
    if hb:
        yr = lambda k: k / (2 * hb["mu"]) if k is not None and np.isfinite(k) else None
        nf = lambda x: float(x) if x is not None and np.isfinite(x) else None
        bio = dict(exchanges=[dict(a=names[s_["i"]], b=names[s_["j"]], start_a=s_["start_i"], end_a=s_["end_i"],
                                   start_b=s_["start_j"], end_b=s_["end_j"], pairs=s_["n"], years=yr(s_["k2p"]),
                                   timing=s_["timing"], p=s_["p"]) for s_ in hb["segments"]],
                   half_life=[dict(lost_in=group_label(e["h"]), copies_from=group_label(e["g"]), k2p=nf(e["est"]),
                                   low=nf(e["lo"]), high=nf(e["hi"]), years=yr(e["est"]), qc_pass=e["ok"])
                              for e in hb["halflife"]],
                   lineages=[dict(subgenomes=[group_label(s_["a"]), group_label(s_["b"])], ratio=nf(s_["ratio"]),
                                  low=nf(s_["lo"]), high=nf(s_["hi"]), differ=s_["differ"]) for s_ in hb["asym"]],
                   translocated_blocks=len(hb["transloc"]),
                   several_homoeologs=[dict(sequence=names[m["i"]], partners=[names[p_] for p_ in m["partners"]])
                                       for m in hb["many"]],
                   collinear_in_reverse=[[names[i], names[j]] for i, j in hb["reversed"]])
    return dict(old_k2p_edge=ho["old_edge"], old_copies=ho["n_old"], anchors=int(len(ho["P"])), pairs=pairs,
                sets=[{group_label(g): [names[x] for x in s_["members"] if z[x] == g] for g in range(N)}
                      | dict(complete=s_["complete"]) for s_ in ho["sets"]], biology=bio)


def write_homoeologs(a, lib, ph):
    """homoeologs.tsv (each phased sequence x other subgenome: status, the called homoeolog, the best candidate and its
    evidence: anchors, observed/expected, longest collinear chain with its z against NULL_PERMS permutations and its
    p (chain_test: Monte Carlo, resolution alpha/2, or the exact bound for <= CHAIN_W anchors), orientation, spans), homoeolog_sets.tsv (connected sets of calls, one column per subgenome) and homoeolog_anchors.tsv (every
    shared old insertion: both copies, whether on its pair's longest collinear chain, whether the pair is called).
    Nothing is written for one group."""
    ho = ph.get("homoeologs")
    if not ho:
        return
    od, z, N, names, ek = Path(a.outdir), np.asarray(ph["z"]), ph["N"], ph["names"], np.asarray(ph["ek"])
    order = sorted(range(len(names)), key=lambda i: natural_key(names[i]))
    with open(od / "homoeologs.tsv", "w") as f:
        f.write("sequence\tsubgenome\tother_subgenome\tstatus\thomoeolog\treciprocal\tbest_candidate\tanchors\t"
                "anchors_oe\tchain\tchain_z\tchain_p\torientation\tspan\tcandidate_span\n")
        for i in order:
            for h in range(N):
                if h == z[i]:
                    continue
                c = ho["calls"][(i, h)]
                j = c["candidate"]
                r = ho["pair"].get((i, j)) if j >= 0 else None
                yn = "NA" if c["reciprocal"] is None else ("yes" if c["reciprocal"] else "no")
                ev = ("\t".join(["NA"] * 9) if r is None else
                      f"{names[j]}\t{r['n']}\t{ho['M'][i, j] / ho['E'][i, j]:.3g}\t{r['chain']}\t{r['z']:.3g}\t"
                      f"{c['p']:.3g}\t{r['orientation']}\t{r['span']:.3f}\t{r['partner_span']:.3f}")
                f.write(f"{names[i]}\t{group_label(int(z[i]))}\t{group_label(h)}\t{c['status']}\t"
                        f"{names[c['partner']] if c['partner'] >= 0 else 'NA'}\t{yn}\t{ev}\n")
    with open(od / "homoeolog_sets.tsv", "w") as f:
        f.write("set\t" + "\t".join(group_label(g) for g in range(N)) + "\tcomplete\n")
        for k, s_ in enumerate(ho["sets"]):
            cols = [",".join(sorted((names[x] for x in s_["members"] if z[x] == g), key=natural_key)) or "-"
                    for g in range(N)]
            f.write(f"H{k + 1}\t" + "\t".join(cols) + f"\t{'yes' if s_['complete'] else 'no'}\n")
    hb = ho.get("bio")
    with open(od / "homoeolog_anchors.tsv", "w") as f:
        f.write("sequence_a\tstart_a\tend_a\tsubgenome_a\tsequence_b\tstart_b\tend_b\tsubgenome_b\tk2p_a\tk2p_b\t"
                "in_chain\tpair_called" + ("\taligned_sites\tdivergence\tnear_identical" if hb else "") + "\n")
        P = ho["P"]
        for q in np.lexsort((lib.start[P[:, 0]], ek[P[:, 0]])):
            x, y = P[q]
            f.write(f"{names[ek[x]]}\t{lib.start[x]}\t{lib.end[x]}\t{group_label(int(z[ek[x]]))}\t{names[ek[y]]}\t"
                    f"{lib.start[y]}\t{lib.end[y]}\t{group_label(int(z[ek[y]]))}\t{_f(lib.k2p[x], '{:.5f}')}\t"
                    f"{_f(lib.k2p[y], '{:.5f}')}\t{'yes' if ho['in_chain'][q] else 'no'}\t"
                    f"{'yes' if ho['pair_called'][q] else 'no'}" +
                    (f"\t{hb['sites'][q]}\t{_f(hb['k2p'][q], '{:.5f}')}\t{'yes' if hb['flags'][q] else 'no'}"
                     if hb else "") + "\n")
    if hb:
        write_homoeolog_biology(a, ph, hb)


def write_homoeolog_biology(a, ph, hb):
    """homoeolog_exchanges.tsv (near-identical stretches on called homoeologs), homoeolog_halflife.tsv
    (half-life of intact LTR-RTs per direction), homoeolog_rearrangements.tsv (translocated blocks; sequences with
    several homoeologs; homoeologs collinear in reverse) and findings.txt (plain-language results, homoeolog_findings,
    in the subgenomes layout of FINDINGS_HEADER; each of its lines echoed to the log)."""
    od, z, names = Path(a.outdir), np.asarray(ph["z"]), ph["names"]
    yr = lambda k: k / (2 * hb["mu"]) if k is not None and np.isfinite(k) else None
    lab = lambda i: f"{names[i]}\t{group_label(int(z[i]))}"
    with open(od / "homoeolog_exchanges.tsv", "w") as f:
        f.write("stretch\tsequence_a\tsubgenome_a\tstart_a\tend_a\tsequence_b\tsubgenome_b\tstart_b\tend_b\t"
                "near_identical_pairs\tchain_pairs\tdivergence\tyears\ttiming\tp\tv8_exchanges\n")
        for k, s_ in enumerate(hb["segments"]):
            f.write(f"X{k + 1}\t{lab(s_['i'])}\t{s_['start_i']}\t{s_['end_i']}\t{lab(s_['j'])}\t{s_['start_j']}\t"
                    f"{s_['end_j']}\t{s_['n']}\t{s_['n_chain']}\t{s_['k2p']:.5f}\t{_f(yr(s_['k2p']))}\t{s_['timing']}\t"
                    f"{s_['p']:.3g}\t{','.join(s_['v8']) or '-'}\n")
    with open(od / "homoeolog_halflife.tsv", "w") as f:
        f.write("lost_in\tcopies_from\tsplit_divergence\tinherited_pairs\tretention_1.5x\tretention_2x\told_copies\t"
                "half_life_k2p\tlow\thigh\thalf_life_years\tyears_low\tyears_high\tcopies_per_mb_lost_in\tqc\n")
        for e in hb["halflife"]:
            bad = ",".join(k for k, v in e["qc"].items() if not v)
            n_use = hb["pairs"][(min(e["g"], e["h"]), max(e["g"], e["h"]))]["use"].size
            f.write(f"{group_label(e['h'])}\t{group_label(e['g'])}\t{e['d_s']:.5f}\t{n_use}\t{e['r15']:.4f}\t"
                    f"{e['r2']:.4f}\t{e['n_old']}\t{_f(e['est'])}\t{_f(e['lo'])}\t{_f(e['hi'])}\t{_f(yr(e['est']))}\t"
                    f"{_f(yr(e['lo']))}\t{_f(yr(e['hi']))}\t{hb['density'][e['h']]:.2f}\t{'pass' if e['ok'] else bad}\n")
    with open(od / "homoeolog_rearrangements.tsv", "w") as f:
        f.write("type\tsequence\tsubgenome\tstart\tend\tpartner\tpartner_subgenome\tpartner_start\tpartner_end\t"
                "collinear_pairs\tpairs_tested\tp\n")
        for t in hb["transloc"]:
            f.write(f"translocated_block\t{lab(t['i'])}\t{t['start_i']}\t{t['end_i']}\t{lab(t['j'])}\t{t['start_j']}\t"
                    f"{t['end_j']}\t{t['n']}\t{t['n_tested']}\t{t['p']:.3g}\n")
        for m in hb["many"]:
            f.write(f"several_homoeologs\t{lab(m['i'])}\tNA\tNA\t{','.join(names[p_] for p_ in m['partners'])}\t"
                    f"{group_label(m['h'])}\tNA\tNA\tNA\tNA\tNA\n")
        for i, j in hb["reversed"]:
            f.write(f"collinear_in_reverse\t{lab(i)}\tNA\tNA\t{lab(j)}\tNA\tNA\tNA\tNA\tNA\n")
    for line in write_findings(od / "findings.txt", homoeolog_findings(ph, hb), header=FINDINGS_HEADER):
        log(line)


def homoeolog_findings(ph, hb):
    """Plain-language statements of what the homoeolog pairs say, each with its numbers and how far the benchmark
    trusts it: a Findings, every statement under SECTION."""
    z, names, N = np.asarray(ph["z"]), ph["names"], ph["N"]
    ho = ph["homoeologs"]
    yrs = lambda k: _ma(k / (2 * hb["mu"]))
    nm = lambda i: f"{names[i]} ({group_label(int(z[i]))})"
    out = []
    for g in range(N):
        for h in range(g + 1, N):
            ig, ih = np.flatnonzero(z == g), np.flatnonzero(z == h)
            cg = sum(ho["calls"][(i, h)]["status"] == "homoeolog" for i in ig)
            ch = sum(ho["calls"][(j, g)]["status"] == "homoeolog" for j in ih)
            rest = cg < ig.size or ch < ih.size
            out.append(f"Homoeologs {group_label(g)}-{group_label(h)}: {cg} of {ig.size} {group_label(g)} and {ch} of "
                       f"{ih.size} {group_label(h)} sequences paired" +
                       ("; the others keep too few inherited insertions in collinear order (transposon turnover "
                        "erases them as subgenomes age)." if rest else "."))
    if hb["reversed"]:
        out.append(f"Orientation: {len(hb['reversed'])} homoeolog pair(s) run in opposite directions in this assembly "
                   f"({', '.join(f'{names[i]}-{names[j]}' for i, j in hb['reversed'][:8])}"
                   f"{', ...' if len(hb['reversed']) > 8 else ''}); one member of each is assembled reverse-complemented "
                   f"relative to the other. Benchmark: 252 of 256 called pairs agree with gene order.")
    else:
        out.append("Orientation: every called homoeolog pair runs in the same direction in this assembly.")
    sg_ = hb["segments"]
    if sg_:
        tim = defaultdict(int)
        for s_ in sg_:
            tim[s_["timing"]] += 1
        desc = "; ".join(f"{nm(s_['i'])} {s_['start_i'] / 1e6:.1f}-{s_['end_i'] / 1e6:.1f} Mb = {nm(s_['j'])} "
                         f"{s_['start_j'] / 1e6:.1f}-{s_['end_j'] / 1e6:.1f} Mb ({s_['n']} pairs, "
                         f"{'identical copies: very recent' if s_['k2p'] <= 0 else '~' + yrs(s_['k2p']) + ' ago'})"
                         for s_ in sg_[:6])
        out.append(f"Exchanges: {len(sg_)} stretch(es) where homoeologs carry near-identical copies of inherited "
                   f"insertions, so one homoeolog's sequence replaced the other's after the progenitors split "
                   f"(homoeologous exchange, or introgression before the merger): {desc}"
                   f"{'; ...' if len(sg_) > 6 else ''}. Timing vs the merger: "
                   + ", ".join(f"{v} {k.replace('_', ' ')}" for k, v in sorted(tim.items())) +
                   ". Benchmark: every call checked was real (432 of 432 in simulations, 5 of 5 in real genomes); "
                   "about 9 in 10 simulated replacements >= 1 Mb are found, fewer in real genomes.")
    else:
        out.append("Exchanges: no stretch of near-identical inherited insertions on the called homoeologs (replacing "
                   "exchanges >= 1 Mb are found 9 times in 10 in simulations; smaller or older ones can be missed).")
    for e in hb["halflife"]:
        if e["ok"]:
            out.append(f"Half-life of intact LTR retrotransposons in {group_label(e['h'])}: {e['est']:.4f} LTR "
                       f"divergence [{e['lo']:.4f}-{e['hi']:.4f}] = {yrs(e['est'])} [{yrs(e['lo'])}-{yrs(e['hi'])}] "
                       f"at mu {hb['mu']:.3g} (share of {group_label(e['g'])}'s pre-split copies whose "
                       f"{group_label(e['h'])} ortholog is still intact: {100 * e['r15']:.0f}%); intact copies per Mb "
                       f"in {group_label(e['h'])}: {hb['density'][e['h']]:.1f}.")
        else:
            out.append(f"Half-life in {group_label(e['h'])} (from {group_label(e['g'])}): not estimable "
                       f"({', '.join(k for k, v in e['qc'].items() if not v)} check failed).")
    for s_ in hb["asym"]:
        if s_["differ"]:
            out.append(f"Lineages differ: intact copies persist {s_['ratio']:.2f}x as long in {group_label(s_['a'])} "
                       f"as in {group_label(s_['b'])} [{s_['lo']:.2f}-{s_['hi']:.2f}]. Read this beside each "
                       f"subgenome's intact copies per Mb: a less complete assembly or annotation of one subgenome "
                       f"also shortens its half-life.")
        else:
            out.append(f"Lineages alike: half-lives in {group_label(s_['a'])} and {group_label(s_['b'])} do not differ "
                       f"detectably (ratio {s_['ratio']:.2f} [{s_['lo']:.2f}-{s_['hi']:.2f}]).")
    if hb["halflife"]:
        out.append("Half-life notes: copies the assembly or annotation missed count as lost, so these are lower "
                   "bounds on persistence; in simulations 94% of estimates that pass the checks fall within 25% of the "
                   "truth (splits of 1-5 Myr; very young splits read ~25% short).")
    if hb["transloc"]:
        out.append(f"Translocations: {len(hb['transloc'])} block(s) of a chromosome collinear with a chromosome other "
                   f"than its homoeolog: " + "; ".join(
                       f"{nm(t['i'])} {t['start_i'] / 1e6:.1f}-{t['end_i'] / 1e6:.1f} Mb with {nm(t['j'])} "
                       f"{t['start_j'] / 1e6:.1f}-{t['end_j'] / 1e6:.1f} Mb ({t['n']} pairs)"
                       for t in hb["transloc"][:6]) + ("; ..." if len(hb["transloc"]) > 6 else "") +
                   ". Benchmark: 16 of 19 such blocks in real genomes are supported by gene synteny.")
    if hb["many"]:
        out.append("Sequences with several homoeologs (fusion, fission or translocation): " + "; ".join(
            f"{nm(m['i'])} pairs with {' + '.join(names[p_] for p_ in m['partners'])} ({group_label(m['h'])})"
            for m in hb["many"][:8]) + ("; ..." if len(hb["many"]) > 8 else "") +
            ". Benchmark: 17 of 17 such cases in real genomes are supported by gene synteny.")
    f = Findings()
    for text in out:
        f.add(SECTION, text)
    return f


def figure_data(a, lib, ph, dt, bio, pa):
    """What report.pdf draws: _fig_data's output, its run objects cut to the fields the page reads (lib: n, seqnames,
    seq, start, end; ph: _PH_KEYS and the homoeolog calls with their exchanged stretches and translocated blocks; dt:
    the subgenome pairs; bio: the clade table), without the painting and the copy locations (read by _fig_data only)."""
    P = _fig_data(a, lib, ph, dt, bio, pa)
    ho = ph.get("homoeologs")
    hb = ho.get("bio") if ho else None
    P.update(lib=SimpleNamespace(n=lib.n, seqnames=list(lib.seqnames), seq=lib.seq, start=lib.start, end=lib.end),
             ph=dict({k: ph[k] for k in _PH_KEYS if k in ph},
                     homoeologs=None if not ho else dict(
                         P=ho["P"], in_chain=ho["in_chain"], pair_called=ho["pair_called"], calls=ho["calls"],
                         sets=ho["sets"],
                         bio=None if not hb else dict(flags=hb["flags"], transloc=hb["transloc"],
                                                      segments=hb["segments"]))),
             dt=dict(pairs=dt["pairs"]) if dt else None, bio=dict(clades=P["bio"].get("clades", [])), pa=None)
    P.pop("loc", None)
    return P


def write_report_data(a, lib, ph, dt, bio, pa):
    """report_data.json in the run folder: figure_data, for plot.subgenomes to draw report.pdf from."""
    out = Path(a.outdir) / REPORT_DATA
    dump_data(out, figure_data(a, lib, ph, dt, bio, pa))
    return out


def _enc(x):
    """A JSON value standing for x: arrays (dtype, shape, values as a list; float32 values in their shortest exact
    form), tuples, namespaces and dicts keyed by anything but plain strings as tagged objects, numpy scalars as numbers.
    NaN and infinities are written as JSON's NaN, Infinity and -Infinity."""
    if isinstance(x, np.ndarray):
        flat = x.ravel()
        data = [float(str(v)) for v in flat] if x.dtype == np.float32 else flat.tolist()
        return {"__array__": x.dtype.str, "shape": list(x.shape), "data": data}
    if isinstance(x, np.generic):
        return x.item()
    if isinstance(x, SimpleNamespace):
        return {"__namespace__": {k: _enc(v) for k, v in vars(x).items()}}
    if isinstance(x, tuple):
        return {"__tuple__": [_enc(v) for v in x]}
    if isinstance(x, list):
        return [_enc(v) for v in x]
    if isinstance(x, dict):
        if all(isinstance(k, str) and not k.startswith("__") for k in x):
            return {k: _enc(v) for k, v in x.items()}
        return {"__items__": [[_enc(k), _enc(v)] for k, v in x.items()]}
    return x


def _dec(x):
    """The value _enc wrote."""
    if isinstance(x, list):
        return [_dec(v) for v in x]
    if isinstance(x, dict):
        if "__array__" in x:
            return np.array(x["data"], dtype=np.dtype(x["__array__"])).reshape(x["shape"])
        if "__tuple__" in x:
            return tuple(_dec(v) for v in x["__tuple__"])
        if "__namespace__" in x:
            return SimpleNamespace(**{k: _dec(v) for k, v in x["__namespace__"].items()})
        if "__items__" in x:
            return {_dec(k): _dec(v) for k, v in x["__items__"]}
        return {k: _dec(v) for k, v in x.items()}
    return x


def dump_data(path, obj):
    """Write obj (figure data) to path as one JSON document (_enc)."""
    tmp = Path(path).with_name(Path(path).name + ".tmp")
    tmp.write_text(json.dumps(_enc(obj), separators=(",", ":")) + "\n")
    tmp.replace(path)


def load_data(path):
    """The object dump_data wrote."""
    return _dec(json.loads(Path(path).read_text()))
