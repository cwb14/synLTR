"""synltr.seq.intervals.union: the one union of half-open intervals."""
from synltr.seq.intervals import union


def test_union_merges_overlapping_and_touching_and_drops_empty():
    assert union([]) == []
    assert union([(5, 9), (1, 3), (3, 4), (8, 12), (20, 20), (30, 25), (15, 18)]) == [(1, 4), (5, 12), (15, 18)]
    assert union(iter([(2, 6), (0, 1)])) == [(0, 1), (2, 6)]  # any iterable, sorted on the way
    assert union([(0, 10), (2, 3)]) == [(0, 10)]
