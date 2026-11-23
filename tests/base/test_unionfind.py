"""synltr.base.unionfind.find: one union-find root for the package's connected components."""
import numpy as np

from synltr.base.unionfind import find


def test_find_returns_the_root_and_halves_the_path():
    parent = [0, 0, 1, 2, 3, 5]          # a chain 4 -> 3 -> 2 -> 1 -> 0, and 5 its own root
    assert find(parent, 4) == 0 and find(parent, 5) == 5
    assert parent == [0, 0, 0, 2, 2, 5]  # 4 and 2 now point to their grandparents
    arr = np.array([1, 2, 2])
    assert find(arr, 0) == 2 and arr.tolist() == [2, 2, 2]
