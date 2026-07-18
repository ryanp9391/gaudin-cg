import pytest
import sympy

from gl3_reps import GL3Rep, gt_patterns

REPS = [(1, 0, 0), (2, 0, 0), (1, 1, 0), (2, 1, 0)]


def weyl_dim(lam):
    l1, l2, l3 = lam
    return (l1 - l2 + 1) * (l2 - l3 + 1) * (l1 - l3 + 2) // 2


@pytest.mark.parametrize("lam", REPS)
def test_dimension_matches_weyl_formula(lam):
    patterns = gt_patterns(lam)
    assert len(patterns) == weyl_dim(lam)
    assert len(set(patterns)) == len(patterns)  # no duplicate patterns


def highest_weight_pattern(lam):
    return tuple(tuple(lam[:k]) for k in range(1, len(lam) + 1))


@pytest.mark.parametrize("lam", REPS)
def test_diagonal_generator_highest_weight_eigenvalue(lam):
    rep = GL3Rep(lam)
    hw = highest_weight_pattern(lam)
    hw_index = rep.patterns.index(hw)
    for k in range(1, 4):
        E_kk = rep.generator(k, k)
        assert E_kk[hw_index, hw_index] == lam[k - 1]


@pytest.mark.parametrize("lam", REPS)
def test_adjacent_commutators(lam):
    rep = GL3Rep(lam)
    comm12 = rep.generator(1, 2) * rep.generator(2, 1) - rep.generator(2, 1) * rep.generator(1, 2)
    assert comm12 == rep.generator(1, 1) - rep.generator(2, 2)
    comm23 = rep.generator(2, 3) * rep.generator(3, 2) - rep.generator(3, 2) * rep.generator(2, 3)
    assert comm23 == rep.generator(2, 2) - rep.generator(3, 3)
