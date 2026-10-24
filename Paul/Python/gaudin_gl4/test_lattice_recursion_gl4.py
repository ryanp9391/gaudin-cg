"""exact regression tests for lattice_recursion_gl4 (G8a).  Run:  python -X utf8 test_lattice_recursion_gl4.py   (or pytest), < 60 s."""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
import mpmath as mp
import lattice_recursion_gl4 as L
import slice_ops_gl4 as so

P1 = ((4, 2, 1, 0), (4, 2, 1, 0))
P2 = ((6, 4, 2, 0), (6, 4, 2, 0))


def sec_of(prod, mu):
    return L.sector(prod[0], prod[1], mu)


def test_cut_is_kernel_of_E12_string():
    # two sectors with a non-void cut (W1 != W), one with the cut void
    checked = 0
    for prod in (P1, P2):
        for s in L.product_sectors(*prod):
            assert s.cut_is_kernel(), (prod, s.mu)
            checked += 1
    assert any(len(s.in1) < s.w for s in L.product_sectors(*P2))
    assert checked > 50


def test_w1_invariant_and_contains_prv():
    for prod in (P1, P2):
        for s in L.product_sectors(*prod):
            assert s.w1_invariant()
            if s.m:
                in1 = set(s.in1)
                # every PRV basis vector is supported on W1
                for b in s.prv.basis:
                    assert all(b[a] == 0 for a in range(s.w) if a not in in1)


def test_vandermonde_determinant():
    ratios, n, ncand = L.vandermonde_check(nsamples=20, seed=3)
    assert n == 20
    for (orient, r), cnt in ratios.items():
        assert r == (1 if orient == (0, -1, 1) else -1)


def test_forward_determinacy_best_grading():
    s = sec_of(P1, (6, 4, 3, 1))                                   # w = 6, m = 5, E12 string active in the cut
    bg = L.best_gradings(s, 'W1')
    assert bg and bg[0][0] == 1 and bg[0][3] == (-1, 0, 0)
    a = L.analyze(s.space('W1'), (-1, 0, 0))
    assert a['ok'] and a['band'] == 1 and a['n0'] == 1 and sum(a['sizes']) == len(s.in1)
    s2 = [x for x in L.product_sectors(*P2) if x.mu == (9, 7, 6, 2)][0]          # a structural failure: first step has rank deficit 1
    a2 = L.analyze(s2.space('W1'), (1, 0, 0))
    assert not a2['ok'] and a2['fail'] == 0 and a2['deficit'] == 1 and a2['sizes'][:2] == [2, 6]


def _check_copies(s, grading, ncheck=2):
    cps = so.copies(s.lam, s.lamp, s.mu, dps=50)
    n = 0
    for cp in cps[:ncheck]:
        sf = L.slice_from_triple(s.lam, s.lamp, s.mu, cp['ev'], grading=grading, space='W1', dps=50)
        assert L.slice_err(sf['slice'], cp['slice']) < mp.mpf(10) ** -40
        assert sf['resid'] < mp.mpf(10) ** -40
        n += 1
    return n


def test_slice_from_triple_reproduces_copies():
    mp.mp.dps = 60
    s = sec_of(P1, (6, 4, 3, 1))
    assert _check_copies(s, 'best', 5) == 5
    s2 = [x for x in L.product_sectors(*P2) if x.m >= 7 and L.best_gradings(x, 'W1')[0][0] == 1][0]
    assert _check_copies(s2, 'best', 3) == 3


def test_band2_grading_recursion():
    mp.mp.dps = 60
    for s in L.product_sectors(*P1):
        a = L.analyze(s.space('W1'), (1, 1, 0))
        if a['ok'] and a['band'] == 2 and a['n0'] == 1 and s.m >= 2:
            assert _check_copies(s, (1, 1, 0), 3) >= 2
            return
    raise AssertionError('no band-2 example found')


def test_characteristic_equations():
    mp.mp.dps = 60
    s = sec_of(P1, (6, 4, 3, 1))
    g = (-1, 0, 0)
    sy = L.sym_recursion(s, 'W1', g)
    assert max(L.tdeg(v[0]) for v in sy['psi'].values()) == 2
    polys = [L.prim_int(r[0]) for r in sy['unused'] if r[0] != 0]
    assert len(polys) == 13
    cps = so.copies(s.lam, s.lamp, s.mu, dps=50)
    for cp in cps:
        for p in polys:
            assert abs(L.peval(p, cp['ev'])) < mp.mpf(10) ** -35
    # and at a perturbed triple at least one of them is far from zero
    pert = tuple(x + 1 for x in cps[0]['ev'])
    assert max(abs(L.peval(p, pert)) for p in polys) > 1
    # the ideal has dim W1 = 6
    assert L.groebner_dim(polys, timeout=50)['dim'] == len(s.in1)


def test_n0_gt_1_kernel():
    mp.mp.dps = 60
    s = [x for x in L.product_sectors(*P2) if L.best_gradings(x, 'W1')[0][0] == 2 and x.w <= 7 and x.m >= 2][0]
    r = L.check_n0_sector(s)
    assert r['rank_ok'] == r['ncopies'] and r['ok'] + r['skipped'] == r['ncopies'] and r['maxker'] < mp.mpf(10) ** -40 and r['maxfull'] < mp.mpf(10) ** -40


if __name__ == '__main__':
    t0 = time.perf_counter()
    nfail = 0
    for name, f in list(globals().items()):
        if name.startswith('test_') and callable(f):
            t = time.perf_counter()
            try:
                f()
                print('PASS %-45s %.2fs' % (name, time.perf_counter() - t))
            except Exception as e:
                nfail += 1
                import traceback
                traceback.print_exc()
                print('FAIL', name, e)
    print('total %.1f s, %d failures' % (time.perf_counter() - t0, nfail))
    sys.exit(1 if nfail else 0)
