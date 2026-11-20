"""Population model of `ltrtrace sim pop`: the genealogy of haploid genomes sampled from three taxa and an outgroup,
simulated with msprime, and the true introgressed tracts of each genome.

Demography (diploid effective sizes, times in generations, backward in time):
  extant    A, B, C, O, and the optional lineage populations A1p (sub_a) and C1p (sub_c)
  splits    A from B at t_ab; B and C from BC at t_bc; BC and O from ROOT at t_o; A1p from A and C1p from C at
            their own split times. A1 is sampled from A1p when it exists, C1 from C1p.
  pulse     at ta, a fraction f of the recipient's lineages moves to the donor (forward: donor -> recipient):
            "C>A", "A>C", or a pulse into a lineage population, "C>A1", "C1>A1", "A1>C1"
  bottleneck  optional (Ne, generations): A at Ne for that many generations right after t_ab
Recombination: r, and r_low on [low[0] L, low[1] L) of each simulated sequence (by default its middle 40%): a
per-chromosome rule applied to each region, a heuristic.
Tracts: a sample's introgressed tracts are the intervals the pulse moved from the recipient to the donor in the
lineage of the sample or of one of its ancestors (msprime's migration records; record_migrations=True).

msprime and tskit (the optional extra `sim`) are imported inside the functions that need them, so that
`import ltrtrace.sim` works without them; those functions raise RuntimeError(MSPRIME_HELP).
"""
from __future__ import annotations

import json
import math
import numbers
from dataclasses import dataclass, fields
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from ltrtrace.base.constants import SIM_DEMOG, SIM_PULSES, SIM_RATES
from ltrtrace.seq.intervals import union

if TYPE_CHECKING:
    import msprime
    import tskit

MSPRIME_HELP = ("msprime is required for `ltrtrace sim pop`: install it with pip install -e '.[sim]' from the "
                "repository, or mamba install msprime, or use the conda environment in environment.yml")

PULSES = SIM_PULSES
_LINEAGE = {"A1p": "sub_a", "C1p": "sub_c"}  # lineage population -> the Demog field that creates it


def _number(name: str, v, positive: bool = True) -> float:
    """v, if it is a finite number > 0 (>= 0 when not `positive`); else a ValueError naming the field."""
    if isinstance(v, bool) or not isinstance(v, numbers.Real) or not math.isfinite(v) or v < 0 \
            or (positive and v == 0):
        raise ValueError(f"{name} must be a finite number {'> 0' if positive else '>= 0'}: {v!r}")
    return v


def _numbers(name: str, v, n: int, positive: bool = True) -> tuple[float, ...]:
    if isinstance(v, (str, bytes)) or not hasattr(v, "__len__") or len(v) != n:
        raise ValueError(f"{name} must be {n} numbers: {v!r}")
    return tuple(_number(name, x, positive) for x in v)


@dataclass(frozen=True)
class Demog:
    """Three taxa A, B and C (diploid Ne, generations) plus an outgroup O and optional lineage populations.
    Checked on construction: a bad value raises ValueError."""

    ne_a: float = SIM_DEMOG["ne_a"]
    ne_b: float = SIM_DEMOG["ne_b"]
    ne_c: float = SIM_DEMOG["ne_c"]
    ne_bc: float = SIM_DEMOG["ne_bc"]
    ne_o: float = SIM_DEMOG["ne_o"]
    ne_root: float = SIM_DEMOG["ne_root"]
    t_ab: float = SIM_DEMOG["t_ab"]
    t_bc: float = SIM_DEMOG["t_bc"]
    t_o: float = SIM_DEMOG["t_o"]  # calibrated on real genomes, from their outgroup conflict
    f: float = SIM_DEMOG["f"]
    ta: float = SIM_DEMOG["ta"]
    pulse: str = SIM_DEMOG["pulse"]  # "C>A" | "A>C" | "C>A1" | "C1>A1" | "A1>C1"
    sub_a: tuple[float, float] | None = None  # (split generations, Ne) of A1's own lineage population "A1p"
    sub_c: tuple[float, float] | None = None  # (split generations, Ne) of "C1p"
    bottleneck: tuple[float, float] | None = None  # (Ne, generations) of A right after t_ab
    per_taxon: int = SIM_DEMOG["per_taxon"]

    def __post_init__(self):
        for name in ("ne_a", "ne_b", "ne_c", "ne_bc", "ne_o", "ne_root", "t_ab", "t_bc", "t_o"):
            _number(name, getattr(self, name))
        _number("ta", self.ta, positive=False)
        if isinstance(self.f, bool) or not isinstance(self.f, numbers.Real) or not 0 <= self.f <= 1:
            raise ValueError(f"f must be in [0, 1]: {self.f!r}")
        if not self.t_ab < self.t_bc < self.t_o:
            raise ValueError(f"split times must increase, t_ab < t_bc < t_o: {self.t_ab}, {self.t_bc}, {self.t_o}")
        if isinstance(self.per_taxon, bool) or not isinstance(self.per_taxon, numbers.Integral) \
                or self.per_taxon < 1:
            raise ValueError(f"per_taxon must be an integer >= 1: {self.per_taxon!r}")
        for name, older in (("sub_a", "t_ab"), ("sub_c", "t_bc")):  # a lineage splits from its taxon before
            v = getattr(self, name)                                   # the taxon itself merges
            if v is not None:
                v = _numbers(name, v, 2)
                object.__setattr__(self, name, v)
                if v[0] >= getattr(self, older):
                    raise ValueError(f"{name}: the lineage split ({v[0]:g}) must be below {older} "
                                     f"({getattr(self, older):g})")
        if self.bottleneck is not None:
            v = _numbers("bottleneck", self.bottleneck, 2)
            object.__setattr__(self, "bottleneck", v)
            if v[1] > self.t_ab:
                raise ValueError(f"bottleneck of {v[1]:g} generations is longer than t_ab ({self.t_ab:g})")
        if not isinstance(self.pulse, str) or self.pulse not in PULSES:
            raise ValueError(f"unknown pulse {self.pulse!r}: expected one of {', '.join(PULSES)}")
        if self.f > 0 and not 0 < self.ta < self.t_ab:
            raise ValueError(f"the pulse time ta ({self.ta:g}) must be in (0, t_ab = {self.t_ab:g}) when f > 0")
        for pop in PULSES[self.pulse]:  # strict at f = 0 too: a lineage code always needs its lineage; ta is
            if pop in _LINEAGE:        # read only when a pulse happens (f > 0), as for ta < t_ab above
                sub = getattr(self, _LINEAGE[pop])
                if sub is None:
                    raise ValueError(f"pulse {self.pulse} needs {_LINEAGE[pop]} (the lineage population {pop})")
                if self.f > 0 and not self.ta < sub[0]:
                    raise ValueError(f"pulse {self.pulse}: ta ({self.ta:g}) must be below the split of {pop} "
                                     f"({sub[0]:g})")


@dataclass(frozen=True)
class Rates:
    """Rates per bp (or per element) per generation, and the loss mechanisms' mix. Checked on construction."""

    mu: float = SIM_RATES["mu"]
    r: float = SIM_RATES["r"]
    r_low: float = SIM_RATES["r_low"]
    low: tuple[float, float] = SIM_RATES["low"]  # the low-recombination interval, as fractions of each sequence
    lam: float = SIM_RATES["lam"]
    delta: float = SIM_RATES["delta"]
    indel: float = SIM_RATES["indel"]
    copy_noise: float = SIM_RATES["copy_noise"]
    # loss_split (solo, del_flank, del_partial, replace) and replace_mean (mean replacement length, bp):
    # calibrated on the element losses of real genomes
    loss_split: tuple[float, float, float, float] = SIM_RATES["loss_split"]
    replace_mean: float = SIM_RATES["replace_mean"]

    def __post_init__(self):
        for name in ("mu", "r", "r_low", "lam", "delta", "indel", "copy_noise"):
            _number(name, getattr(self, name), positive=False)
        _number("replace_mean", self.replace_mean)
        low = _numbers("low", self.low, 2, positive=False)
        if not 0 <= low[0] <= low[1] <= 1:
            raise ValueError(f"low must be fractions 0 <= start <= end <= 1: {self.low!r}")
        split = _numbers("loss_split", self.loss_split, 4, positive=False)
        if abs(sum(split) - 1) > 1e-6:
            raise ValueError(f"loss_split must sum to 1: {self.loss_split!r}")
        object.__setattr__(self, "low", low)
        object.__setattr__(self, "loss_split", split)


def _msprime():
    try:
        import msprime
    except ImportError:
        raise RuntimeError(MSPRIME_HELP) from None
    return msprime


def _sample_pops(d: Demog) -> list[tuple[str, str, str]]:
    """(label, taxon, population) of every sample, in sample-node order."""
    out = []
    for taxon, n, lineage in (("A", d.per_taxon, d.sub_a), ("B", d.per_taxon, None), ("C", d.per_taxon, d.sub_c),
                              ("O", 1, None)):
        for i in range(1, n + 1):
            out.append((f"{taxon}{i}", taxon, f"{taxon}1p" if i == 1 and lineage is not None else taxon))
    return out


def sample_labels(d: Demog) -> list[tuple[str, str]]:
    """(label, group) of every sample in sample-node order: A1..An, B1..Bn, C1..Cn (n = per_taxon), then O1;
    the group is the taxon."""
    return [(label, taxon) for label, taxon, _ in _sample_pops(d)]


def demography(d: Demog) -> msprime.Demography:
    """The msprime Demography of `d`. Every extant population is active at time zero (B and the taxa that lineage
    populations split from are ancestral in a split, which msprime would otherwise start inactive)."""
    msprime = _msprime()
    dem = msprime.Demography()
    for name, ne in (("A", d.ne_a), ("B", d.ne_b), ("C", d.ne_c), ("O", d.ne_o)):
        dem.add_population(name=name, initial_size=ne, initially_active=True)
    dem.add_population(name="BC", initial_size=d.ne_bc)
    dem.add_population(name="ROOT", initial_size=d.ne_root)
    for pop, sub, taxon in (("A1p", d.sub_a, "A"), ("C1p", d.sub_c, "C")):
        if sub is not None:
            dem.add_population(name=pop, initial_size=sub[1], initially_active=True)
            dem.add_population_split(time=sub[0], derived=[pop], ancestral=taxon)
    if d.f > 0:
        donor, recipient = PULSES[d.pulse]
        dem.add_mass_migration(time=d.ta, source=recipient, dest=donor, proportion=d.f)
    if d.bottleneck is not None:
        ne, gens = d.bottleneck
        dem.add_population_parameters_change(time=d.t_ab - gens, population="A", initial_size=ne)
    dem.add_population_split(time=d.t_ab, derived=["A"], ancestral="B")
    dem.add_population_split(time=d.t_bc, derived=["B", "C"], ancestral="BC")
    dem.add_population_split(time=d.t_o, derived=["BC", "O"], ancestral="ROOT")
    dem.sort_events()
    return dem


def recombination_map(rates: Rates, length: int) -> msprime.RateMap:
    """r on [0, length), with r_low on [low[0] * length, low[1] * length) (boundaries rounded to whole bp)."""
    msprime = _msprime()
    if length < 1:
        raise ValueError(f"sequence length must be >= 1 bp: {length}")
    a, b = round(rates.low[0] * length), round(rates.low[1] * length)
    position, rate = [0], []
    for end, x in ((a, rates.r), (b, rates.r_low), (length, rates.r)):
        if end > position[-1]:
            position.append(end)
            rate.append(x)
    return msprime.RateMap(position=position, rate=rate)


def ancestry(d: Demog, rates: Rates, length: int, seed: int) -> tskit.TreeSequence:
    """The genealogy of the samples of `d` (sample_labels order) on one sequence of `length` bp: haploid samples
    of diploid populations, the recombination map of `rates`, migrations recorded (for `tracts`). `seed` is
    msprime's random_seed, 1 to 2**32 - 1 (msprime raises ValueError otherwise)."""
    msprime = _msprime()
    counts: dict[str, int] = {}
    for _, _, pop in _sample_pops(d):
        counts[pop] = counts.get(pop, 0) + 1
    samples = [msprime.SampleSet(n, population=pop, ploidy=1) for pop, n in counts.items()]
    return msprime.sim_ancestry(samples=samples, demography=demography(d), ploidy=2,
                                recombination_rate=recombination_map(rates, length), record_migrations=True,
                                random_seed=seed)


def tracts(ts: tskit.TreeSequence, d: Demog) -> dict[int, list[tuple[int, int]]]:
    """Sample index -> the merged introgressed intervals [l, r) of that sample on this sequence (one call per
    simulated sequence); samples without a tract are left out, so no pulse gives {}. A tract is where a
    migration record of the pulse (time ta, recipient to donor) holds the sample or one of its ancestors."""
    if d.f == 0:
        return {}
    donor, recipient = PULSES[d.pulse]
    meta = [ts.population(i).metadata for i in range(ts.num_populations)]
    pid = {m.get("name"): i for i, m in enumerate(meta) if isinstance(m, dict)}
    if donor not in pid or recipient not in pid:
        raise ValueError(f"the tree sequence has no population {donor!r} or {recipient!r}: not made from this "
                         "demography")
    sel = np.flatnonzero((ts.migrations_source == pid[recipient]) & (ts.migrations_dest == pid[donor])
                         & (ts.migrations_time == d.ta))
    if not len(sel):
        return {}
    sel = sel[np.argsort(ts.migrations_left[sel], kind="stable")]
    left, right, node = (ts.migrations_left[sel].tolist(), ts.migrations_right[sel].tolist(),
                         ts.migrations_node[sel].tolist())
    index = {int(u): i for i, u in enumerate(ts.samples())}
    pieces: dict[int, list[tuple[int, int]]] = {}
    nxt, active = 0, []
    for tree in ts.trees():
        tl, tr = tree.interval
        while nxt < len(left) and left[nxt] < tr:
            active.append(nxt)
            nxt += 1
        active = [k for k in active if right[k] > tl]
        for k in active:
            lo, hi = max(left[k], tl), min(right[k], tr)
            if hi > lo:
                for s in tree.samples(node[k]):
                    pieces.setdefault(index[int(s)], []).append((int(lo), int(hi)))
        if nxt == len(left) and not active:
            break
    out = {}
    for s in sorted(pieces):
        out[s] = union(pieces[s])
    return out


def flat_fields(*classes) -> list[set[str]]:
    """The field names of each dataclass, for one flat JSON object holding them side by side; a name two of them
    share would feed both, so it raises ValueError naming it."""
    names = [{f.name for f in fields(c)} for c in classes]
    for i, a in enumerate(names):
        for j in range(i + 1, len(names)):
            if a & names[j]:
                raise ValueError(f"{classes[i].__name__} and {classes[j].__name__} share field(s) "
                                 f"{', '.join(sorted(a & names[j]))}: a flat JSON object cannot hold both")
    return names


_JSON_FIELDS = flat_fields(Demog, Rates)  # checked on import: a shared field name fails fast


def demog_from_json(path: Path) -> tuple[Demog, Rates]:
    """Demog and Rates from one JSON object whose keys are their field names (a field left out keeps its
    default); an unknown key, or a value Demog or Rates rejects, raises ValueError."""
    try:
        data = json.loads(Path(path).read_text())
    except json.JSONDecodeError as e:
        raise ValueError(f"{path}: not valid JSON ({e})") from None
    if not isinstance(data, dict):
        raise ValueError(f"{path}: expected a JSON object of Demog and Rates fields")  # noqa: TRY004 (bad input)
    dn, rn = _JSON_FIELDS
    unknown = sorted(set(data) - dn - rn)
    if unknown:
        raise ValueError(f"{path}: unknown field(s) {', '.join(unknown)}; the fields are Demog's "
                         f"({', '.join(f.name for f in fields(Demog))}) and Rates' "
                         f"({', '.join(f.name for f in fields(Rates))})")
    try:
        return (Demog(**{k: v for k, v in data.items() if k in dn}),
                Rates(**{k: v for k, v in data.items() if k in rn}))
    except ValueError as e:
        raise ValueError(f"{path}: {e}") from None
