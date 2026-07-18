import pytest

from gl3_reps import gt_patterns

REPS = [(1, 0, 0), (2, 0, 0), (1, 1, 0), (2, 1, 0)]


def weyl_dim(lam):
    l1, l2, l3 = lam
    return (l1 - l2 + 1) * (l2 - l3 + 1) * (l1 - l3 + 2) // 2


@pytest.mark.parametrize("lam", REPS)
def test_dimension_matches_weyl_formula(lam):
    patterns = gt_patterns(lam)
    assert len(patterns) == weyl_dim(lam)
    assert len(set(patterns)) == len(patterns)  # no duplicate patterns
