"""Small inputs for the `ltrtrace sim pop` tests.

small_demog: the default demography (Demog) with every effective size cut tenfold and the outgroup split at
90,000 generations, so that one ancestry of the seven default genomes on L = 200 kb takes well under a
second. The split times stay the defaults (t_ab 9,000, t_bc 60,000), so a test can place pulses and lineage
splits against them.
toy_ts: an exact one-tree genealogy from a nested tuple. toy_backbone: random bases with planted elements.
toy_element and toy_donors: elements and donors of given shapes. tskit is imported only inside toy_ts.

Hand-built events for the builder and truth tests (test_build, test_truth): TREE7 and toy7 (an exact
7-sample genealogy and its internal nodes), DONORS (three fixed donors), ins_at and loss_at (one event at a point
of the tree), events (a ChromEvents of them, with SNPs and small indels), build (build_chrom on them, EXACT rates
by default) and other (a base other than a given one)."""
from __future__ import annotations

import random
from collections import Counter
from dataclasses import replace

import numpy as np

from ltrtrace.cli import main as _cli_main
from ltrtrace.sim.pop.build import build_chrom
from ltrtrace.sim.pop.donors import Donor
from ltrtrace.sim.pop.events import ChromEvents, Ins, Loss
from ltrtrace.sim.pop.genome import SimElement
from ltrtrace.sim.pop.model import Demog, Rates

L = 200_000  # bp, the toy chromosome


def sim_main(argv: list[str]) -> int:
    """`ltrtrace sim ARGV...` through the command line; the exit code (usage errors and --help included)."""
    return _cli_main(["sim", *argv])

_SMALL = Demog(ne_a=3_000, ne_b=10_000, ne_c=7_000, ne_bc=15_000, ne_o=15_000, ne_root=15_000, t_o=90_000)


def small_demog(**changes) -> Demog:
    """The toy demography with `changes` applied (Demog checks them)."""
    return replace(_SMALL, **changes)


def _leaves(tree) -> list[int]:
    return [tree] if isinstance(tree, int) else [x for k in tree[1:] for x in _leaves(k)]


def toy_ts(tree, length: int = L):
    """A tree sequence of one tree on [0, length) (time units: generations). `tree` is a nested tuple: a leaf is
    a sample index (node of that id, time 0; the leaves are 0..n-1), an internal node is (time, child, child,
    ...). Internal nodes are numbered from n in postorder."""
    import tskit

    leaves = sorted(_leaves(tree))
    if leaves != list(range(len(leaves))):
        raise ValueError(f"the leaves must be 0..n-1 once each: {leaves}")
    tables = tskit.TableCollection(sequence_length=length)
    tables.time_units = "generations"
    for _ in leaves:
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)

    def build(t) -> tuple[int, float]:
        if isinstance(t, int):
            return t, 0.0
        time, *kids = t
        built = [build(k) for k in kids]
        if not all(time > kt for _, kt in built):
            raise ValueError(f"node time {time} not above its children's")
        u = tables.nodes.add_row(time=time)
        for c, _ in built:
            tables.edges.add_row(0, length, u, c)
        return u, time

    build(tree)
    tables.sort()
    return tables.tree_sequence()


def tree_length(tree) -> float:
    """Total branch length of a toy_ts nested tuple."""
    if isinstance(tree, int):
        return 0.0
    time, *kids = tree
    return sum(time - (0.0 if isinstance(k, int) else k[0]) + tree_length(k) for k in kids)


def rand_bases(n: int, rng: random.Random) -> bytes:
    return "".join(rng.choices("ACGT", k=n)).encode()


def toy_element(eid: str, start: int, ltr: int, internal: int) -> SimElement:
    """A backbone SimElement [start, start + 2 ltr + internal) with two `ltr`-bp LTRs (coordinates only)."""
    end = start + 2 * ltr + internal
    return SimElement(eid, start, end, (start, start + ltr), (end - ltr, end), "+", "", None, None, None, "backbone")


def toy_backbone(length: int, spans=(), seed: int = 0) -> tuple[bytes, list[SimElement]]:
    """Random ACGT bases of `length` with an element planted at each (start, ltr, internal) of `spans`: an LTR,
    an internal region and the same LTR, with a 5-bp TSD (the bases before it, copied after it). Elements are
    bb0001, bb0002, ... in (start, -end) order, origin "backbone"; a span inside another's internal region is
    nested in it (parent: the smallest element holding it)."""
    rng = random.Random(seed)
    seq = bytearray(rand_bases(length, rng))
    els = sorted((toy_element("", s, ltr, n) for s, ltr, n in spans), key=lambda e: (e.start, -e.end))
    out = []
    for i, e in enumerate(els, 1):
        if not (5 <= e.start and e.end + 5 <= length):
            raise ValueError(f"span at {e.start} does not fit the backbone with its TSD")
        ltr = rand_bases(e.ltr_left[1] - e.ltr_left[0], rng)
        seq[e.start:e.end] = ltr + rand_bases(e.ltr_right[0] - e.ltr_left[1], rng) + ltr
        seq[e.end:e.end + 5] = seq[e.start - 5:e.start]
        hosts = [h for h in out if h.ltr_left[1] <= e.start and e.end <= h.ltr_right[0]]
        parent = min(hosts, key=lambda h: h.end - h.start).eid if hosts else None
        out.append(replace(e, eid=f"bb{i:04d}", tsd=bytes(seq[e.start - 5:e.start]).decode(), parent=parent))
    return bytes(seq), out


def toy_donors(n: int = 4, seed: int = 0, ltr: tuple[int, int] = (200, 600),
               internal: tuple[int, int] = (2_000, 9_000)) -> list[Donor]:
    """n random donors in their own orientation: each LTR (left and right drawn apart) and the internal region
    uniform in the given bp ranges."""
    rng = random.Random(seed)
    out = []
    for i in range(n):
        left, right, internal_len = rng.randint(*ltr), rng.randint(*ltr), rng.randint(*internal)
        out.append(Donor(rand_bases(left + internal_len + right, rng).decode(), left, right, f"fam{i}", "Gypsy"))
    return out


# exact genealogy: 7 samples; internal nodes 01 (3,000), 23 (5,000), 0123 (12,000), 45 (8,000), 456 (20,000), root
TREE7 = (40_000.0, (12_000.0, (3_000.0, 0, 1), (5_000.0, 2, 3)), (20_000.0, (8_000.0, 4, 5), 6))
EXACT = Rates(mu=0.0, copy_noise=0.0)  # no substitutions anywhere: payloads are donor bodies


def fixed_donors() -> list[Donor]:
    """Bodies of 2,600, 2,000 (left LTR 250 twice; the donor's right LTR is 260) and 3,300 bp."""
    r = random.Random(5)
    shapes = [(300, 2_000, 300), (250, 1_500, 260), (400, 2_500, 400)]
    return [Donor(rand_bases(a + b + c, r).decode(), a, c, f"fam{i}", "Gypsy")
            for i, (a, b, c) in enumerate(shapes)]


DONORS = fixed_donors()


def toy7(length: int = L):
    ts = toy_ts(TREE7, length)
    tree = ts.first()
    n = {"01": tree.parent(0), "23": tree.parent(2), "45": tree.parent(4)}
    n["0123"], n["456"] = tree.parent(n["01"]), tree.parent(n["45"])
    return ts, tree, n


def below(tree, node: int) -> frozenset[int]:
    index = {int(s): i for i, s in enumerate(tree.tree_sequence.samples())}
    return frozenset(index[int(s)] for s in tree.samples(node))


def ins_at(tree, k: int, host, pos: int, node: int, time: float, depth: int = 0, strand: str = "+",
           donor: int = 0) -> Ins:
    return Ins(f"i00{k:06d}", host, pos, node, float(time), donor, strand, depth, below(tree, node))


def loss_at(tree, element: str, node: int, time: float, kind: str, a: int, b: int = 0) -> Loss:
    return Loss(element, node, float(time), kind, a, b, below(tree, node))


def events(bb: bytes, n: int = 7, ins=(), losses=(), snps=None, indels=()) -> ChromEvents:
    """Hand-built events. snps: position -> {sample: base} (the backbone base elsewhere); indels: (pos, length,
    is_insertion, carriers)."""
    snps = snps or {}
    pos = np.array(sorted(snps), dtype=np.int64)
    base = np.empty((len(pos), n), dtype=np.uint8)
    for j, p in enumerate(pos.tolist()):
        base[j] = bb[p]
        for s, b in snps[p].items():
            base[j, s] = ord(b)
    return ChromEvents(list(ins), list(losses), pos, base, np.array([i[0] for i in indels], dtype=np.int64),
                       np.array([i[1] for i in indels], dtype=np.int64), np.array([i[2] for i in indels], dtype=bool),
                       [frozenset(i[3]) for i in indels], Counter())


def build(ts, bb, els, ev, rates=EXACT, seed=0, donors=DONORS):
    return build_chrom(0, bb, list(els), ts, ev, donors, rates, np.random.default_rng(seed))


def other(b: int) -> str:
    """A base other than ASCII code b."""
    return "C" if chr(b) == "A" else "A"
