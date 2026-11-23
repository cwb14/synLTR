"""The phasing pipeline: sequences with enough LTR-RTs are phased on the links of their clades (a fixed partition with
--assign), near-identical sequence pairs are kept apart, the structure is decided on the core sequences and the others
are placed after; the verdict and the tested splits are recorded."""
from __future__ import annotations

import numpy as np

from synltr.analyses.subgenomes.halves import (
    Sources,
    assign_satellites,
    chrom_tests,
    core_sequences,
    find_subgenomes,
    partner_pairs,
    support_from_halves,
)
from synltr.analyses.subgenomes.homoeologs import mask_pairs, pair_splits, paired_sets, pairing_map, run_homoeologs
from synltr.analyses.subgenomes.inputs import group_label, log, quiet, vlog, warn
from synltr.analyses.subgenomes.links import LinkModel
from synltr.analyses.subgenomes.partition import ari, canon_labels
from synltr.base.constants import CORE_FRAC, FIG_HALVES_MAX, MAX_N, MIN_ELEMENTS, MIN_REP, SUPPORT_OK
from synltr.base.log import UsageError


def verdict_text(N, inside, big, tests, source):
    """One-line verdict. inside: paired sequence pairs left in one final group; big: groups of >= 3 near-identical
    sequences (paired map components)."""
    if N >= 2:
        v = "subgenomes found" + (" (borderline)" if any(t["accepted"] and t["chosen"]["borderline"] for t in tests)
                                  else "")
        nf = sum(t["forced"] for t in tests)
        if nf:
            v += f"; {nf} of the requested splits fail the split rule (forced by -n)"
        if inside >= 2:
            v += f"; {inside} pairs of near-identical sequences share a group: more subgenomes than resolved"
        if big:
            v += f"; {len(big)} groups of 3 or more near-identical sequences"
        return v
    if inside >= 2 or big:
        return (f"cannot tell: near-identical copies of {inside + len(big)} chromosomes (homologues of a young "
                f"polyploid, or both haplotypes of a diploid assembly); LTR-RT lineages do not separate them")
    root = tests[0]["cands"].get("all") if tests else None
    if source == "given" or (root is not None and not root["borderline"]):
        return "no structure detected"
    return "cannot tell"


def read_assign(path, seqnames):
    """--assign file: sequence<whitespace>label per line ('#' comments; further columns ignored) -> {sequence: label}
    in file order. Conflicting labels for one sequence, a file naming no sequence of the library and more than MAX_N
    groups are UsageErrors naming the file; one-token lines and sequences the library does not hold are warned about
    and skipped; a single group is warned about (nothing to date or paint)."""
    out, bad, conflict = {}, 0, []
    with open(path) as f:
        for line in f:
            x = line.split("#", 1)[0].split()
            if len(x) == 1:
                bad += 1
            elif len(x) >= 2:
                if x[0] in out and out[x[0]] != x[1]:
                    conflict.append(f"{x[0]} ({out[x[0]]}, {x[1]})")
                out.setdefault(x[0], x[1])
    if conflict:
        raise UsageError(f"--assign {path}: {len(conflict)} sequence(s) given two labels: " + ", ".join(conflict[:5]))
    if not out:
        raise UsageError(f"--assign {path}: no 'sequence label' lines")
    if bad:
        warn(f"--assign: {bad} line(s) with a single field skipped")
    known = set(seqnames)
    miss = [s for s in out if s not in known]
    if len(miss) == len(out):
        raise UsageError(f"--assign {path}: none of its sequences is in the library (file: {', '.join(miss[:3])}; "
                         f"library: {', '.join(list(seqnames)[:3])})")
    if miss:
        warn(f"--assign: {len(miss)} sequence(s) not in the library (skipped): " + ", ".join(miss[:10]))
    out = {s: g for s, g in out.items() if s in known}
    n = len(set(out.values()))
    if n > MAX_N:
        raise UsageError(f"--assign {path}: {n} groups; at most {MAX_N} (is the second column a subgenome label?)")
    if n == 1:
        warn("--assign: one group only; nothing to date or paint")
    return out


def _shown_candidate(t):
    """The candidate a tested split is reported by: the chosen one, else the pooled one, else the most replicable.
    A split rejected by the homoeolog check is reported by the age-class candidate the check rejected."""
    if t["chosen"] is not None:
        return t["chosen"]
    if t.get("gated") is not None:
        return t["gated"]
    if "all" in t["cands"]:
        return t["cands"]["all"]
    return max(t["cands"].values(), key=lambda c: (round(c["rep"], 2), c["con"])) if t["cands"] else None


def tree_row(t, names, sources):
    """structure.tsv / summary.json record of one tested split: the reported candidate and every source's numbers."""
    c = _shown_candidate(t)
    g = t["members"]
    row = dict(level=t["path"], n=int(len(g)), copies=c["src"] if c else None, rule=t["how"],
               accepted=bool(t["accepted"]), forced=bool(t["forced"]),
               part_a=[names[i] for i in g[c["z"] == 0]] if c else [names[i] for i in g],
               part_b=[names[i] for i in g[c["z"] == 1]] if c else [],
               replicability=c["rep"] if c else None, replicability_se=c["se"] if c else None,
               n_pairs=c["n_pairs"] if c else 0, borderline=bool(c["borderline"]) if c else False,
               strength=c["con"] if c else None,
               parent_strength=t["parent"]["con"][c["src"]] if c and t["parent"] else None)
    for s_ in sources:
        x = t["cands"].get(s_)
        row[f"replicability_{s_}"] = x["rep"] if x else None
        row[f"strength_{s_}"] = x["con"] if x else None
        row[f"agrees_{s_}"] = (bool(ari(x["z"], c["z"]) >= MIN_REP) if x and c else None)
        row[f"parent_strength_{s_}"] = t["parent"]["con"][s_] if t["parent"] else None
    return row


def split_log(r):
    """One log line per tested split: replicability/strength on each source (parent strength in brackets)."""
    src = [k[len("replicability_"):] for k in r if k.startswith("replicability_") and k != "replicability_se"]
    parts = []
    for s_ in src:
        if r[f"replicability_{s_}"] is None:
            continue
        par = f" ({r[f'parent_strength_{s_}']:.2f})" if r[f"parent_strength_{s_}"] is not None else ""
        parts.append(f"{s_} {r[f'replicability_{s_}']:.2f}/{r[f'strength_{s_}']:.2f}{par}")
    verdict = (("forced by -n" if r["forced"] else "accepted" + (" (age-bounded, " + r["copies"] + " copies)"
                                                               if r["rule"] == "age" else ""))
               if r["accepted"] or r["forced"] else "rejected: " + r["rule"])
    return (f"split {len(r['part_a'])}|{len(r['part_b'])} [{r['level']}] replicability/strength: " + "; ".join(parts)
            + f" -> {verdict}")


def split_homoeologs(a, lib, cl, ek, edges, g, zg):
    """Reciprocal homoeolog calls between the two parts (zg 0/1) of a candidate split of the sequences g: the
    homoeolog step (run_homoeologs, unchanged) on the sub-problem of g's sequences with the split as two subgenomes.
    It draws from its own random stream, so phasing draws are unchanged; its log lines (and anything else it prints)
    are silenced."""
    import contextlib
    import io
    g, zg = np.asarray(g, np.int64), np.asarray(zg).astype(int)
    pos = np.full(int(max(int(ek.max()), int(g.max()))) + 1, -1, np.int64)
    pos[g] = np.arange(g.size)
    ek_sub = np.where(ek >= 0, pos[np.maximum(ek, 0)], -1)
    with quiet(), contextlib.redirect_stderr(io.StringIO()):
        hom = run_homoeologs(a, lib, dict(N=2, C=int(g.size), strata_edges=edges, z=zg, ek=ek_sub, clades_all=cl))
    if not hom:
        return 0
    return int(sum(1 for (i, h), c in hom["calls"].items()
                   if zg[i] == 0 and c["status"] == "homoeolog" and c["reciprocal"]))


def phasable(a, lib):
    """(the --assign partition or None, the indexes of the sequences to phase): --assign read and checked against the
    library (read_assign), and the sequences with >= MIN_ELEMENTS LTR-RTs (those --assign names, when given). Fewer
    than two, or than -n, is a UsageError naming --assign, -n or the library. It needs only the library, so a run
    calls it before the k-mer stage: a fault in these inputs fails in seconds."""
    n_el_all = np.bincount(lib.seq, minlength=len(lib.seqnames))
    user = read_assign(a.assign, lib.seqnames) if a.assign else None      # fixed partition: no structure search
    phased = [i for i in range(len(lib.seqnames)) if n_el_all[i] >= MIN_ELEMENTS and
              (user is None or lib.seqnames[i] in user)]
    if user and len(phased) < len(user):
        warn(f"--assign: {len(user) - len(phased)} sequence(s) carry < {MIN_ELEMENTS} LTR-RTs (unphased)")
    N = a.n_subgenomes
    if len(phased) < max(N or 2, 2):
        what = f"--assign {a.assign}" if user else f"-n {N}" if N else ", ".join(map(str, a.ltr_fasta))
        raise UsageError(f"{what}: only {len(phased)} sequences carry >= {MIN_ELEMENTS} LTR-RTs; cannot form "
                         f"{max(N or 2, 2)} groups")
    return user, phased


def run_phasing(a, lib, cl, sel=None):
    """Phase the library's sequences on the links of its clades `cl`; `sel` is phasable(a, lib), computed here when
    not given."""
    rng = np.random.default_rng(a.seed)
    user, phased = phasable(a, lib) if sel is None else sel
    N = a.n_subgenomes
    names = [lib.seqnames[i] for i in phased]
    N_requested = N
    C = len(names)
    pos = {s: k for k, s in enumerate(phased)}
    ek = np.array([pos.get(int(s), -1) for s in lib.seq], np.int64)
    log(f"phasing {C} sequences with >= {MIN_ELEMENTS} LTR-RTs" +
        (f" ({len(lib.seqnames) - C} others listed as unphased)" if len(lib.seqnames) > C else ""))
    lm = LinkModel(cl, ek, C, a.threads)
    if lm.cl.n == 0:
        raise UsageError(f"{', '.join(map(str, a.ltr_fasta))}: no shared k-mers link two sequences; is this an LTR-RT "
                         "library?")
    informative = np.zeros(lib.n, bool)
    informative[lm.cl.members] = True
    n_inf = np.bincount(ek[(ek >= 0) & informative], minlength=C)
    pairing = pairing_map(lib, lm, np.random.default_rng([a.seed, 3]), a.threads)
    cl_pairs, big = paired_sets(pairing, C)
    log(f"pairing map: {sum(r['significant'] for r in pairing)} sequence pairs with collinear shared insertions, "
        f"{sum(r['paired'] for r in pairing)} paired (near-identical sequences)")
    src = Sources(lm, lib.k2p, rng)
    if src.edges:
        log(f"link sources: all copies, and LTR-K2P terciles young <= {src.edges[0]:.4g} < middle <= "
            f"{src.edges[1]:.4g} < old ({', '.join(f'{s_} {int(m.sum()):,}' for s_, m in src.masks.items())} copies)")
    W = W_raw = src.W["all"]
    n_el = np.bincount(ek[ek >= 0], minlength=C)
    tree, tests, source, use_pairs = [], [], "given" if N else "inferred", []
    support = np.full(C, np.nan)
    core = np.ones(C, bool) if user else core_sequences(n_el)       # the structure is decided on the core
    if not core.all():
        log(f"core: {int(core.sum())} sequences with >= {CORE_FRAC:g} x the copy-N50 LTR-RTs decide the structure "
            f"({int(n_el[core].sum()):,} LTR-RTs); {int((~core).sum())} smaller sequences "
            f"({int(n_el[~core].sum()):,} LTR-RTs) are placed afterwards")
    if user:
        here = {user[s] for s in names}
        order = [g for g in dict.fromkeys(user.values()) if g in here]    # labels in the file's order
        z, N, source = np.array([order.index(user[s]) for s in names]), len(order), "user"
        log(f"partition given by the user: {N} group(s) ({', '.join(order)} -> "
            f"{', '.join(group_label(g) for g in range(N))})")
    else:
        sets_of = None
        if N and cl_pairs:                                        # forced N: paired sequences kept apart
            use_pairs = list(cl_pairs)

            def sets_of(g):
                p_ = {x: i for i, x in enumerate(g)}
                return [[p_[x], p_[y]] for x, y in use_pairs if x in p_ and y in p_]
        age_ok = None
        if src.edges:                                        # an age-class sub-split leaves homoeologs across

            def age_ok(g_, zg_):
                n_ = split_homoeologs(a, lib, cl, ek, src.edges, g_, zg_)
                vlog(f"age-class split of {len(g_)} sequences: {n_} reciprocal homoeolog pairs across it")
                return n_ > 0
        lab, tests, support = find_subgenomes(src, rng, n_force=N, sets_of=sets_of, root=np.flatnonzero(core),
                                              age_ok=age_ok)
        if N is None and cl_pairs:
            lab, ptests = pair_splits(src, lab, cl_pairs, rng, support)
            tests += ptests
            use_pairs = [p for t in ptests if t["accepted"] for p in t["inside"]]
        if use_pairs:
            log(f"paired sequences: {len(use_pairs)} pairs kept in different groups")
        for t in tests:
            if t["how"] == "too small":                            # groups of < 4 sequences are never split
                continue
            tree.append(tree_row(t, names, src.names))
            log(split_log(tree[-1]))
        if not core.all():                                        # satellites join the core's groups
            lab, sat_sup = assign_satellites(src, lab, np.flatnonzero(~core))
            support[~core] = sat_sup[~core]
            if int(lab.max()) >= 1:
                vlog(f"satellites: {int((sat_sup[~core] >= SUPPORT_OK).sum())} of {int((~core).sum())} placed with "
                     f"support >= {SUPPORT_OK:g}")
        z = canon_labels(lab, n_el)
        if N and int(z.max()) + 1 < N:
            warn(f"only {int(z.max()) + 1} of the {N} requested groups could be formed")
        N = int(z.max()) + 1
        log(f"{'inferred' if source == 'inferred' else 'requested'} number of subgenomes: {N}")
    if use_pairs:                                                  # paired links carry no subgenome signal
        W = mask_pairs(W_raw, use_pairs)
    halves = [mask_pairs(np.asarray(w, np.float64), use_pairs) for p in src.pairs["all"] for w in p]
    one = N == 1
    inside = int(sum(z[x] == z[y] for x, y in cl_pairs))            # paired sequences left in one final group
    verdict = "partition given by the user" if user else verdict_text(N, inside, big, tests, source)
    log(f"verdict: {verdict}")
    user_labels = {group_label(g): order[g] for g in range(len(order))} if user else None
    common = dict(names=names, C=C, ek=ek, tree=tree, W=W, W_raw=W_raw, W_src=dict(src.W), n_el=n_el, n_inf=n_inf,
                  halves=[np.asarray(h, np.float32) for h in halves] if C <= FIG_HALVES_MAX else None,   # report.pdf C
                  source=source, user_labels=user_labels, clades=lm.cl, clades_all=cl, N_requested=N_requested or N,
                  n_clades=lm.cl.n, pairing=pairing, cl_pairs=use_pairs, paired_groups=big, verdict=verdict,
                  n_pairs=len(src.pairs["all"]), strata_edges=src.edges,
                  strata_n={s_: int(m.sum()) for s_, m in src.masks.items()}, core=core)
    if one:
        partner, psup, pex = partner_pairs(W, halves, z)
        return dict(common, z=z, N=1, support=np.full(C, np.nan), partners=(partner, psup, pex),
                    own_oe=np.full(C, np.nan), oth_oe=np.full(C, np.nan), pval=np.full(C, np.nan))
    if user:
        support = support_from_halves(halves, z, N, [[i] for i in range(len(names))], rng)
    own_oe, oth_oe, pval = chrom_tests(W, z, N)
    partner, psup, pex = partner_pairs(W, halves, z)
    for g in range(N):
        log(f"{group_label(g)}: {int((z == g).sum())} sequences, {int(n_el[z == g].sum()):,} LTR-RTs")
    return dict(common, z=z, N=N, support=support, partners=(partner, psup, pex), own_oe=own_oe, oth_oe=oth_oe,
                pval=pval)
