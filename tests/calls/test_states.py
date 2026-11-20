import ast
from pathlib import Path

import pytest

from ltrtrace.calls.states import (
    ALLELES,
    FLAGS,
    INFERRED_ABSENT,
    OBSERVED_ABSENT,
    PAIR_STATES,
    PRESENT,
    call_class,
)

SRC = Path(__file__).resolve().parents[2] / "src" / "ltrtrace"


def test_alleles():
    assert ALLELES == "PSTEADOGNU" and PRESENT == {"P", "S", "T"}
    assert OBSERVED_ABSENT == {"E"} and INFERRED_ABSENT == {"A"}
    assert "".join(PAIR_STATES) == "PSTEADOGU"  # N is a reads-route state, never in a pair call


@pytest.mark.parametrize("alleles, want", [
    ("P", "present"), ("S", "present"), ("T", "present"), ("PE", "present"), ("TEA", "present"),
    ("E", "absent"), ("EE", "absent"), ("ED", "absent"), ("EG", "absent"),
    ("A", "missing"), ("AD", "missing"),
    ("D", "missing"), ("O", "missing"), ("G", "missing"), ("N", "missing"), ("U", "missing"), ("", "missing"),
])
def test_call_class(alleles, want):
    assert call_class(alleles) == want


@pytest.mark.parametrize("alleles, want", [
    ("A", "absent"), ("E", "absent"), ("AE", "absent"), ("AD", "absent"), ("PA", "present"),
    ("D", "missing"), ("O", "missing"), ("G", "missing"), ("N", "missing"), ("U", "missing"),
])
def test_call_class_with_inferred_absence(alleles, want):
    assert call_class(alleles, absent=OBSERVED_ABSENT | INFERRED_ABSENT) == want


def test_call_class_takes_any_iterable_and_rejects_unknown_letters():
    assert call_class(iter(["E", "D"])) == "absent"
    with pytest.raises(ValueError, match="unknown allele"):
        call_class("Q")


def test_flags_registry():
    assert {"tsd_match", "conflict", "no_call", "multi_member", "host_empty"} <= set(FLAGS)
    assert all(FLAGS.values())


def _allele_tables(tree: ast.AST):
    """Literal tuples, lists, sets and dicts of three or more allele letters, and the string "PST"."""
    letters = set(ALLELES)
    for node in ast.walk(tree):
        if isinstance(node, (ast.Tuple, ast.List, ast.Set)):
            items = node.elts
        elif isinstance(node, ast.Dict):
            items = node.keys
        elif isinstance(node, ast.Constant) and node.value == "PST":
            yield node
            continue
        else:
            continue
        if len(items) >= 3 and all(isinstance(x, ast.Constant) and x.value in letters for x in items):
            yield node


def test_one_state_table():
    """The mapping from states to classes is written once, in calls/states.py."""
    found = []
    for path in sorted(SRC.rglob("*.py")):
        if path == SRC / "calls" / "states.py":
            continue
        found += [f"{path.relative_to(SRC)}:{n.lineno}" for n in _allele_tables(ast.parse(path.read_text()))]
    assert not found, f"state tables outside calls/states.py: {found}"
