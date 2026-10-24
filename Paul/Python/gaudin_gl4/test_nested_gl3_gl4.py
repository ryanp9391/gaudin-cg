"""test_nested_gl3_gl4.py -- exact regression tests for nested_gl3_gl4.py (WP-G8b); pytest-style functions AND runnable as `python -X utf8 test_nested_gl3_gl4.py` (< 60 s)."""
import os, sys, time
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mpmath as mp
import glN_reps as gr
import nested_gl3_gl4 as ng


def test_b0_identification():
    r = ng.ident_check((3, 2, 1, 0), verbose=False)
    assert r['bad'] == 0 and r['leave'] == 0 and r['entries'] > 500 and r['kappas'] == 8, r
    return 'B0 on [3,2,1,0]: kappa-block of E12,E21,E23,E32,E11,E22,E33,E13,E31 = gl(3) GLNRep(kappa) exactly (%d entries, 8 kappas)' % r['entries']


def test_within_path_tridiagonal():
    lam = (3, 2, 1, 0)
    n = 0
    for mu, m, w in ng.sectors(lam, lam):
        if w < 2:
            continue
        for kappa, d in ng.within_data(lam, lam, mu).items():
            assert all(v == 0 for v in d['tri'].values())
            n += 1
    return 'all %d path blocks of [3,2,1,0]^2 are tridiagonal for V3, V4a, V4b' % n


def test_b1_verdict():
    # V3 and V4a: J = a V^(3) + b exactly (dominant sector of [6,4,2,0]^2, 3-point path); V4b: not, but J = a V + e [xi,V] + b + c xi + g xi^2
    lam = (6, 4, 2, 0)
    d = ng.within_data(lam, lam, (10, 8, 5, 1))[(6, 4, 1)]
    assert d['n'] == 3
    f3 = ng.fit_ansatz(d['J']['V3'], d['V3'], d['D'])
    f4 = ng.fit_ansatz(d['J']['V4a'], d['V3'], d['D'])
    fb = ng.fit_ansatz(d['J']['V4b'], d['V3'], d['D'])
    assert f3['ok1'] and f3['a'] == Fr(10, 9) and f4['ok1'] and f4['a'] == Fr(20, 9) and not fb['ok2']
    x = ng.within_V4b_fit(d)
    assert x is not None and x[1] == Fr(-1, 72) and x[4] == Fr(1, 288) and x[4] == -x[1] / 4
    # smallest failing example of V4b
    lam2 = (2, 1, 1, 0)
    d2 = ng.within_data(lam2, lam2, (3, 2, 2, 1))[(2, 1, 0)]
    assert ng.fit_ansatz(d2['J']['V3'], d2['V3'], d2['D'])['ok1'] and not ng.fit_ansatz(d2['J']['V4b'], d2['V3'], d2['D'])['ok1']
    assert d2['J']['V4b'] == [[Fr(55, 6), Fr(15, 2)], [Fr(17, 2), Fr(-21, 2)]]
    return 'B1: V3, V4a = a V^(3) + b (a = 10/9, 20/9 at kappa=(6,4,1)); V4b not, fits a V + e[xi,V] + ... with e = -1/72, g = -e/4; failing example [2,1,1,0]^2 printed'


def test_mtilde_direct_sum():
    n = 0
    for lam, lamp, mu in (((4, 2, 1, 0), (4, 2, 1, 0), (6, 4, 2, 2)), ((6, 4, 2, 0), (6, 4, 2, 0), (10, 7, 5, 2)), ((5, 3, 2, 1), (5, 3, 2, 1), (9, 6, 5, 2))):
        b = ng.b2_exact(lam, lamp, mu)
        assert b['direct_sum'] and b['dim_sum'] == b['rs_sum'] == b['dimM'], b
        n += 1
    return 'M~ = (+)_kappa M^(3)_kappa (dim = sum of Racah-Speiser m3) exactly on %d sectors' % n


def test_triple_preserves_mtilde():
    for lam, mu in (((4, 2, 1, 0), (6, 4, 2, 2)), ((6, 4, 2, 0), (10, 7, 4, 3))):
        b = ng.b2_exact(lam, lam, mu)
        assert b['preserved'] and b['Mmu_in_Mt'] and b['dimM'] > b['m'], b       # dimM > m: the E34 string is non-void here
    return 'the triple preserves M~ and M_mu in M~ (exact), with dim M~ > m on n >= k3 sectors'


def test_mt_equals_m_iff_void():
    lam = (4, 2, 1, 0)
    for mu, m, w in ng.sectors(lam, lam):
        if m == 0:
            continue
        mt = ng.mtilde_space(lam, lam, mu).m
        assert (mt == m) == ng.void3_of(lam, lam, mu), (mu, m, mt)
    return 'dim M~ = m exactly when E34^k3 is void (all sectors of [4,2,1,0]^2)'


def test_depth1_hopping_universal():
    lam = (3, 2, 1, 0)
    mu = (5, 4, 2, 1)
    c = ng.lattice_hop_check(lam, lam, mu)
    assert c['tests'] >= 18 and c['bad'] == 0 and c['nonzero'] >= 16, dict(c)
    return 'depth-1 sector (5,4,2,1) of [3,2,1,0]^2: %d universal-formula tests (%d nonzero), 0 mismatches' % (c['tests'], c['nonzero'])


def test_depth1_nested_solution():
    n = 0
    for lam, mu in (((3, 2, 1, 0), (5, 4, 2, 1)), ((4, 2, 1, 0), (6, 4, 3, 1)), ((6, 4, 2, 0), (10, 8, 5, 1))):
        v = ng.validate_nested_copies(lam, lam, mu)
        assert v['ok'], v
        n += v['n']
    return 'depth-1 nested solution (nested basis + hopping) = so.copies on 3 sectors (%d copies, 1e-30)' % n


def run_all():
    tests = [test_b0_identification, test_within_path_tridiagonal, test_b1_verdict, test_mtilde_direct_sum, test_triple_preserves_mtilde, test_mt_equals_m_iff_void, test_depth1_hopping_universal,
             test_depth1_nested_solution]
    nfail = 0
    t00 = time.perf_counter()
    mp.mp.dps = 50
    for t in tests:
        t0 = time.perf_counter()
        try:
            msg = t()
            print('PASS  %-36s %.1fs  %s' % (t.__name__, time.perf_counter() - t0, msg), flush=True)
        except Exception as e:
            nfail += 1
            import traceback
            print('FAIL  %-36s %s' % (t.__name__, repr(e)), flush=True)
            traceback.print_exc()
    print('%d tests, %d failed (%.1f s)' % (len(tests), nfail, time.perf_counter() - t00))
    return nfail


if __name__ == '__main__':
    sys.exit(1 if run_all() else 0)
