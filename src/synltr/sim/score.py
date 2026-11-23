"""Score an ordered pair's calls (`synltr genotype`: <X>__<Y>.calls.tsv.gz) against the simulator's truth.tsv.

Calls join truth through call.eid -> Element.name (the simulator's eid, carried by the GFF3 Name) -> truth row
(genome X, name). Elements the annotation dropped (dropout) have a truth row but no element: they are excluded
from every count (`n_unannotated`). An annotated element without a call row is scored as call U
(`n_missing_call`). Truth never holds U. Two tier sets: "all" calls, and "high" (only high-tier calls count as
calls; a truth state they miss is a false negative).
"""
from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

from synltr.calls.records import stream_calls
from synltr.calls.states import ALLELES, OBSERVED_ABSENT, PRESENT
from synltr.io.elements import read_elements
from synltr.io.files import open_input

# the columns of a truth table (poptruth writes them, score reads them)
TRUTH_COLUMNS = ("genome", "eid", "name", "seqid", "start", "end", "origin", "state", "flags", "events_x",
                 "events_y")
STATES = tuple(a for a in ALLELES if a not in "AUN")  # what a truth table holds: no inferred absence, reads state or U
ORDER = (*STATES, "U", "not_high")
CLASSES = {"present": tuple(a for a in STATES if a in PRESENT),
           "absent": tuple(a for a in STATES if a in OBSERVED_ABSENT)}
RESOLVED = (*CLASSES["present"], *CLASSES["absent"])  # states with a defined answer (the rest: missing)
FLAGS = ("inv", "host_empty")
TIER_SETS = ("all", "high")
SCORE_COLUMNS = ("state", "tier_set", "tp", "fp", "fn", "precision", "recall", "n_scored",
                 "correct_resolved_frac", "unresolved_frac", "n_unannotated", "n_missing_call")
CONFUSION_COLUMNS = ("tier_set", "truth", "call", "n")
_CALL_FIELDS = ("eid", "state", "tier", "flags")  # what scoring reads of a pair call file


def _read(path: Path, needed: tuple[str, ...]) -> list[dict[str, str]]:
    """The rows of a TSV, plain or gzip (told by its first bytes, not its name)."""
    with open_input(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in needed if c not in (reader.fieldnames or ())]
        if missing:
            raise ValueError(f"{path}: missing column(s) {', '.join(missing)}")
        return list(reader)


def _unique(rows: list[dict[str, str]], key: str, path: Path, allowed: tuple[str, ...]) -> dict[str, dict]:
    """{row[key]: row}; fails on a key seen twice or a state outside `allowed` (no row silently lost)."""
    out: dict[str, dict] = {}
    for r in rows:
        if r[key] in out:
            raise ValueError(f"{path}: {key} {r[key]!r} appears more than once")
        if r["state"] not in allowed:
            raise ValueError(f"{path}: {key} {r[key]!r} has state {r['state']!r}, not one of {', '.join(allowed)}")
        out[r[key]] = r
    return out


def _flags(cell: str | None) -> set[str]:
    return {f for f in (cell or "").split(",") if f and f != "."}


def _pairs(calls_path: Path, elements_path: Path, truth_path: Path):
    """(n_unannotated, [(truth row, call row or None)]) per annotated element, after the join checks."""
    elements = read_elements(elements_path)
    genomes = {e.genome for e in elements}
    if len(genomes) != 1:
        raise ValueError(f"{elements_path}: expected elements of exactly one genome, found "
                         f"{sorted(genomes) or 'none'}")
    x = genomes.pop()
    truth = _unique([r for r in _read(truth_path, TRUTH_COLUMNS) if r["genome"] == x], "name", truth_path, STATES)
    name_of = {e.eid: e.name for e in elements}
    calls = {}
    rows = [dict(zip(_CALL_FIELDS, f, strict=True)) for _, f in stream_calls(calls_path, _CALL_FIELDS)]
    for eid, r in _unique(rows, "eid", calls_path, (*STATES, "U")).items():
        if eid not in name_of:
            raise ValueError(f"{calls_path}: call eid {eid!r} is not in {elements_path}")
        calls[name_of[eid]] = r
    names = {e.name for e in elements}
    absent = sorted(names - truth.keys())
    if absent:
        raise ValueError(f"{truth_path}: {len(absent)} annotated element(s) of genome {x} have no truth row "
                         f"(first: {', '.join(absent[:5])}); do the inputs belong to one simulation?")
    return len(truth) - len(names), [(truth[n], calls.get(n)) for n in sorted(names)]


def _call(row: dict[str, str] | None, tier_set: str) -> str:
    """The call state of `row` under the tier set: U for a missing row, not_high for a non-high call in
    "high"."""
    if row is None:  # scored as U, tier "."
        return "not_high" if tier_set == "high" else "U"
    return "not_high" if tier_set == "high" and row["tier"] != "high" else row["state"]


def _ratio(num: int, den: int) -> float | None:
    return num / den if den else None


def _prf(tp: int, fp: int, fn: int) -> tuple[float | None, float | None]:
    return _ratio(tp, tp + fp), _ratio(tp, tp + fn)


def score(calls_path: Path, elements_path: Path, truth_path: Path) -> list[dict]:
    """Rows (SCORE_COLUMNS; None where a cell does not apply): per tier set, per state P S T E D O G and per
    class (present = P/S/T, absent = E) tp/fp/fn/precision/recall, flag rows (inv, host_empty) over the
    elements called correctly, and one summary row."""
    n_unannotated, pairs = _pairs(calls_path, elements_path, truth_path)
    n_missing = sum(c is None for _, c in pairs)
    rows = []
    for ts in TIER_SETS:
        obs = [(t, c, t["state"], _call(c, ts)) for t, c in pairs]
        groups = [(s, (s,)) for s in STATES] + list(CLASSES.items())
        for label, members in groups:
            tp = sum(tr in members and ca in members for _, _, tr, ca in obs)
            fp = sum(tr not in members and ca in members for _, _, tr, ca in obs)
            fn = sum(tr in members and ca not in members for _, _, tr, ca in obs)
            rows.append(_row(label, ts, tp=tp, fp=fp, fn=fn))
        right = [(t, c) for t, c, tr, ca in obs if ca == tr]
        for flag in FLAGS:
            both = sum(flag in _flags(t["flags"]) and flag in _flags(c["flags"]) for t, c in right)
            call_only = sum(flag not in _flags(t["flags"]) and flag in _flags(c["flags"]) for t, c in right)
            truth_only = sum(flag in _flags(t["flags"]) and flag not in _flags(c["flags"]) for t, c in right)
            rows.append(_row(f"flag:{flag}", ts, tp=both, fp=call_only, fn=truth_only))
        resolved = [o for o in obs if o[2] in RESOLVED]
        rows.append(_row("summary", ts, n_scored=len(pairs),
                         correct_resolved_frac=_ratio(sum(tr == ca for _, _, tr, ca in resolved),
                                                      len(resolved)),
                         unresolved_frac=_ratio(sum(c is None or c["state"] == "U" for _, c, _, _ in obs),
                                                len(pairs)),
                         n_unannotated=n_unannotated, n_missing_call=n_missing))
    return rows


def _row(state: str, tier_set: str, **cells) -> dict:
    row = dict.fromkeys(SCORE_COLUMNS)
    row.update(state=state, tier_set=tier_set, **cells)
    if "tp" in cells:
        row["precision"], row["recall"] = _prf(cells["tp"], cells["fp"], cells["fn"])
    return row


def confusion(calls_path: Path, elements_path: Path, truth_path: Path) -> list[dict]:
    """Rows (CONFUSION_COLUMNS) for every observed (truth, call) pair per tier set; in "high" a call that is
    not high tier counts as "not_high". Sorted by tier set, truth, call (state order P S T E D O G U)."""
    _, pairs = _pairs(calls_path, elements_path, truth_path)
    out = []
    for ts in TIER_SETS:
        n = Counter((t["state"], _call(c, ts)) for t, c in pairs)
        out += [{"tier_set": ts, "truth": tr, "call": ca, "n": k}
                for (tr, ca), k in sorted(n.items(), key=lambda kv: (ORDER.index(kv[0][0]),
                                                                    ORDER.index(kv[0][1])))]
    return out
