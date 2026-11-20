"""Fixtures of the subgenomes tests: run settings from the flag names the tests use, a command-line run on a one-row
samples sheet (cached per inputs and options), small simulated polyploids (the builders' polyploid scenario, cached
per scenario), synthetic libraries, and stub run objects for the figure."""
import argparse
import hashlib
import shutil
import tempfile
from pathlib import Path
from types import SimpleNamespace

import builders
import numpy as np

from ltrtrace import cli
from ltrtrace.analyses.subgenomes import inputs, run
from ltrtrace.base import constants as cst
from ltrtrace.io.seqio import encode

SIM_TINY = {"root_age": 6e6, "chromosomes": {"n": 4, "length": 4e6}, "families": {"n": 60},
            "activity": {"rate": 0.001, "p_burst": 0.0}, "exchanges": {"nonreciprocal": 0, "reciprocal": 0}}

_SIMS: dict = {}  # one simulation per scenario and seed per session
_RUNS: dict = {}  # one run per inputs and options per session


def _flags() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(allow_abbrev=False)
    ap.add_argument("--ltr_fasta", nargs="+", required=True)
    ap.add_argument("--annotation", nargs="+", required=True)
    ap.add_argument("--fai", nargs="+", default=None)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("-n", "--n_subgenomes", type=int, default=None)
    ap.add_argument("--assign", default=None)
    ap.add_argument("--mu", type=float, default=None)
    ap.add_argument("--label", default=None)
    ap.add_argument("-t", "--threads", type=int, default=4)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("-v", "--verbose", action="store_true")
    return ap


def parse_args(argv):
    """The run settings (run.options) of `--ltr_fasta F.. --annotation A.. [--fai X..] --outdir D [-n N] [--assign F]
    [--mu M] [--label L] [-t T] [--seed S]` (-v, which acts through the log level, only cli_args passes on)."""
    x = _flags().parse_args(argv)
    return run.options(x.ltr_fasta, x.annotation, x.outdir, fai=x.fai, n=x.n_subgenomes, assign=x.assign, mu=x.mu,
                       label=x.label, threads=x.threads, seed=x.seed)


def sheet(path, ltr_fasta, annotation, fai=None, label=None) -> Path:
    """A one-row samples sheet: label, library (the LTR-RT FASTA(s)), annotation and the .fai in the fasta column."""
    path = Path(path)
    label = label or Path(ltr_fasta[0]).name.split(".")[0]
    row = [label, ",".join(str(p) for p in ltr_fasta), str(annotation[0]), str(fai[0]) if fai else "."]
    path.write_text("label\tlibrary\tannotation\tfasta\n" + "\t".join(row) + "\n")
    return path


def cli_args(argv):
    """`ltrtrace subgenomes` arguments for the flags of parse_args (the sheet written beside the output folder)."""
    x = _flags().parse_args(argv)
    out = Path(x.outdir)
    out.parent.mkdir(parents=True, exist_ok=True)
    sh = sheet(out.parent / f"{out.name}.sheet.tsv", x.ltr_fasta, x.annotation, x.fai, x.label)
    args = ["subgenomes", str(sh), "-o", str(out), "-t", str(x.threads), "--seed", str(x.seed)]
    args += ["-n", str(x.n_subgenomes)] if x.n_subgenomes is not None else []
    args += ["--assign", str(x.assign)] if x.assign else []
    args += ["--mu", repr(x.mu)] if x.mu is not None else []
    return args + (["-v"] if x.verbose else [])


def main(argv):
    """Run `ltrtrace subgenomes` in this process on the flags of parse_args; SystemExit with the exit code when it is
    not 0."""
    code = cli.main(cli_args(argv))
    if code:
        raise SystemExit(code)


def small_library(d, counts=(12, 12), shared=True, k2p=None, fai_seqs=None, seed=1):
    """A small library in d: counts[s] 300-bp copies on sequence c<s+1>, copies of six random families with 2%
    substitutions (k-mers link the sequences; unrelated random records without `shared`), named
    c<s>:<start>-<end>#LTR/Gypsy/Tekay; its K2P table k2p.tsv (copy j of sequence s: k2p(s, j), default
    0.01 + 0.001 j) and genome.fai listing fai_seqs (default the library's sequences). Returns (ltr.fa, k2p.tsv,
    genome.fai)."""
    import random

    rng = random.Random(seed)
    d = Path(d)
    d.mkdir(parents=True, exist_ok=True)
    fams = ["".join(rng.choice("ACGT") for _ in range(300)) for _ in range(6)]
    fa, tab = [], ["#seq_id\tk2p"]
    for s, n in enumerate(counts, 1):
        for j in range(n):
            if shared:
                seq = "".join(b if rng.random() >= 0.02 else rng.choice("ACGT") for b in fams[j % len(fams)])
            else:
                seq = "".join(rng.choice("ACGT") for _ in range(300))
            name = f"c{s}:{1000 * j + 1}-{1000 * j + 300}#LTR/Gypsy/Tekay"
            fa.append(f">{name}\n{seq}\n")
            tab.append(f"{name}\t{0.01 + 0.001 * j if k2p is None else k2p(s, j)}")
    (d / "ltr.fa").write_text("".join(fa))
    (d / "k2p.tsv").write_text("\n".join(tab) + "\n")
    seqs = fai_seqs or [f"c{s}" for s in range(1, len(counts) + 1)]
    (d / "genome.fai").write_text("".join(f"{c}\t{1000 * max(counts) + 1000}\t0\t60\t61\n" for c in seqs))
    return d / "ltr.fa", d / "k2p.tsv", d / "genome.fai"


def simulate(tmp_path, seed=1, scenario=None, **kw):
    """A fresh copy of a simulated polyploid (`sim nucleus`; scenario SIM_TINY with kw's keys replaced) in
    tmp_path/in and its truth in tmp_path/truth; each scenario is simulated once per session. Returns (in, {chromosome:
    subgenome})."""
    sc = dict(scenario or SIM_TINY, **kw)
    key = (seed, repr(sorted(sc.items())))
    if key not in _SIMS:
        src = Path(tempfile.mkdtemp(prefix="sgsim_", dir=tmp_path.parent))
        builders.simulate(sc, src / "in", src / "truth", seed)
        _SIMS[key] = src
    for sub in ("in", "truth"):  # a fresh copy: tests may edit their inputs
        shutil.copytree(_SIMS[key] / sub, tmp_path / sub, dirs_exist_ok=True)
    return tmp_path / "in", builders.truth_subgenomes(tmp_path / "truth")


def _digest(paths) -> str:
    h = hashlib.sha1()
    for p in paths:
        h.update(Path(p).read_bytes() if Path(p).exists() else b"-")
    return h.hexdigest()


def run_tool(d, out, *extra, annotation="annotation.gff3", cache=True):
    """`ltrtrace subgenomes` on the library in d (ltr.fa, annotation, genome.fai) into out, with extra flags; with
    cache, a run with the same inputs and flags is made once per session and copied (its summary.json then names the
    first run's input paths). Returns (chromosomes.tsv rows, summary)."""
    import json

    d, out = Path(d), Path(out)
    argv = ["--ltr_fasta", str(d / "ltr.fa"), "--annotation", str(d / annotation), "--fai", str(d / "genome.fai"),
            "--outdir", str(out), "-t", "4", *extra]
    key = (_digest([d / "ltr.fa", d / annotation, d / "genome.fai"]), annotation, tuple(extra))
    if not cache:
        main(argv)
    elif key not in _RUNS:
        main(argv)
        _RUNS[key] = out
    elif _RUNS[key] != out:
        shutil.copytree(_RUNS[key], out, dirs_exist_ok=True)
    rows = [line.split("\t") for line in (out / "chromosomes.tsv").read_text().splitlines()]
    return [dict(zip(rows[0], r)) for r in rows[1:]], json.loads((out / "summary.json").read_text())


def wrong(rows, truth):
    """Sequences whose subgenome label maps to another true subgenome (labels matched to truths one to one)."""
    from collections import Counter

    pairs = Counter((r["subgenome"], truth[r["sequence"]]) for r in rows if r["subgenome"] != "NA")
    best = {}
    for (g, t), n in pairs.most_common():
        if g not in best and t not in best.values():
            best[g] = t
    return sum(1 for r in rows if r["subgenome"] != "NA" and best.get(r["subgenome"]) != truth[r["sequence"]])


def dating_rows(out):
    """dating.tsv's rows as dicts (comment lines skipped)."""
    lines = [line.split("\t") for line in (Path(out) / "dating.tsv").read_text().splitlines()
             if not line.startswith("#")]
    return [dict(zip(lines[0], line)) for line in lines[1:]]


def lib_from_seqs(seqs, seqof):
    """A library of the given sequences (2-bit codes, record offsets) on the given sequence indices."""
    lib = inputs.Library()
    lib.n = len(seqs)
    chunks = [encode(s.encode()) for s in seqs]
    lib.codes = np.concatenate(chunks)
    lib.offs = np.concatenate(([0], np.cumsum([c.size for c in chunks]))).astype(np.int64)
    lib.seq = np.array(seqof, np.int32)
    return lib


def random_family(rng, n_copies=40, L=300, mut=0.02):
    """n_copies copies of one random L-bp sequence, each with mut random substitutions per base."""
    base = rng.choice(list("ACGT"), L)
    out = []
    for _ in range(n_copies):
        s = base.copy()
        hit = rng.random(L) < mut
        s[hit] = rng.choice(list("ACGT"), hit.sum())
        out.append("".join(s))
    return out


def mutate(rng, s, n_ts, n_tv):
    """s with n_ts transitions and n_tv transversions at distinct random sites."""
    ts = {"A": "G", "G": "A", "C": "T", "T": "C"}
    tv = {"A": "C", "C": "A", "G": "T", "T": "G"}
    s = list(s)
    sites = rng.choice(len(s), n_ts + n_tv, replace=False)
    for k, p in enumerate(sites):
        s[p] = ts[s[p]] if k < n_ts else tv[s[p]]
    return "".join(s)


def figure_stub(N, S, prefix="c", per=40, dating=False, source="inferred"):
    """Synthetic run objects for the figure alone: S sequences in N groups (lengths rising with the index), per copies
    each with LTR ages, links stronger within groups (every source and 6 half-libraries), one called exchange, painting
    posteriors and resolution tracks, insertion histories and clade rows; optionally a dated subgenome pair."""
    rng = np.random.default_rng(5)
    names = [f"{prefix}{i}" for i in range(S)]
    z = np.arange(S) % N
    ek = np.repeat(np.arange(S), per)
    n = ek.size
    start = np.tile(np.arange(per) * 25000 + 1000, S)
    lib = SimpleNamespace(n=n, seq=ek.astype(np.int32), start=start, end=start + 8000, seqnames=names,
                          seq_len=per * 25000 + 5000 + np.arange(S) * 1000, k2p=rng.uniform(0, 0.05, n))
    W = rng.uniform(1, 2, (S, S)) + 5 * (z[:, None] == z[None, :])
    W = (W + W.T) / 2
    np.fill_diagonal(W, 0)
    src = ("all",) + cst.AGE_STRATA
    halves = [np.asarray(W * rng.uniform(0.4, 0.6, W.shape), np.float32) for _ in range(6)]
    halves = [np.asarray((h + h.T) / 2, np.float32) for h in halves]
    root = dict(level="root", n=S, copies="all", rule="pooled", accepted=N > 1, forced=False,
                part_a=[s for s, g in zip(names, z) if g < max(N // 2, 1)],
                part_b=[s for s, g in zip(names, z) if g >= max(N // 2, 1)] if N > 1 else names[S // 2:],
                **{f"{k}_{s}": v for s in src for k, v in (("replicability", 0.95 if N > 1 else 0.3),
                                                         ("strength", 1.0), ("agrees", True),
                                                         ("parent_strength", None))})
    if N == 1:
        root["part_a"] = names[:S // 2]
    ph = dict(N=N, names=names, z=z, ek=ek, W=W, W_raw=W, W_src={s: W for s in src}, halves=halves,
              strata_edges=[0.01, 0.03], verdict="subgenomes found" if N > 1 else "no structure detected",
              source=source, support=np.ones(S) if N > 1 else np.full(S, np.nan),
              tree=[root] if source != "user" else [])
    pa = None
    if N > 1:
        post = []
        for i, s in enumerate(names):
            p_ = np.full((per, N), 0.02 / max(N - 1, 1))
            p_[:, z[i]] = 0.98
            post.append((s, (start[:per] + 4000).astype(np.int64), p_, int(lib.seq_len[i])))
        segs = [dict(seq=names[0], start=101000, end=1209000, own=0, origin=1 % N, left_flank=90000,
                     right_flank=1220000, n_ltr=5, n_support=5, n_against=0, n_lineages=3, llr=20.0, post=0.9,
                     called=True)]
        pa = dict(segs=segs, post=post, res=[(s, np.array([0]), np.array([1]), np.array([200000.0])) for s in names])
    from ltrtrace.analyses.subgenomes.biology import HIST_EDGES

    hist, clades = [], []
    for g in range(N):
        for c in ("all", "Tekay", "Athila"):
            for b in range(HIST_EDGES.size - 1):
                hist.append(dict(subgenome=g, clade=c, lo=HIST_EDGES[b], hi=HIST_EDGES[b + 1],
                                 observed=int(50 * np.exp(-((b - 20 - g) / 5) ** 2)),
                                 deconvolved=float(50 * np.exp(-((b - 20 - g) / 5) ** 2))))
        clades += [dict(clade=c, subgenome=g, n=300 + 10 * g) for c in ("Tekay", "Athila", "unknown")]
    bio = dict(history=hist, clades=clades)
    dt = None
    if dating and N > 1:
        P = rng.integers(0, n, (4000, 2))
        P = P[ek[P[:, 0]] != ek[P[:, 1]]]
        d = rng.uniform(1e-4, 0.3, len(P))

        def e(est, lo, hi):
            return dict(estimate=est, lo=lo, hi=hi, n=100)

        pairs = [dict(groups=(x, y), split=dict(orthologs=None, confinement=e(0.04, 0.03, 0.05),
                                                profile=e(0.05, 0.03, 0.1)),
                      split_reported="profile", merger=dict(e(0.001, 0.0, 0.002), datable=True, bound=True,
                                                            source="genealogy_clock"),
                      merger_ltr=dict(e(0.002, 0.0, 0.003), datable=True, bound=True), merger_genealogy=None)
                 for x in range(N) for y in range(x + 1, N)]
        dt = dict(mu=1.3e-8, mu_given=True, pairs=pairs, pt=dict(P=P, d=d))
    a = SimpleNamespace(seed=1, label="Testus", mu=None, outdir=None, ltr_fasta=["x.fa"])
    return a, lib, ph, dt, bio, pa
