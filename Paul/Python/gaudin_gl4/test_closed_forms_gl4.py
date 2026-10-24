"""test_closed_forms_gl4.py -- exact regression tests for closed_forms_gl4.py (WP-G8d); pytest-style functions AND runnable as `python -X utf8 test_closed_forms_gl4.py` (< 60 s once
closed_forms_gl4_qcache.pkl exists; the first run builds that cache from G4's symbolic term sums, about 3 minutes)."""
import os, sys, time, random
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mpmath as mp
import glN_reps as gr
import slice_ops_gl4 as so
import closed_forms_gl4 as cf


def test_a_formula():
    # held-out kappa of lam not used in the fits: a = 1 - sum R_i, two-point-path data of G3's operator
    lam = (10, 7, 4, 1)
    ks = [(8, 5, 3), (9, 6, 2), (7, 5, 3), (10, 6, 2), (8, 7, 2)]
    for k in ks:
        res, nu = cf.a_data(lam, lam, k)
        assert res['V3'][0] == cf.a_of(k, lam), (k, res['V3'][0], cf.a_of(k, lam))
    # explicit values of G8b's table for lam = (6,4,2,0)
    assert cf.a_of((5, 3, 1), (6, 4, 2, 0)) == Fr(32, 35) and cf.a_of((5, 2, 1), (6, 4, 2, 0)) == Fr(464, 315) and cf.a_of((6, 4, 2), (6, 4, 2, 0)) == 1
    return 'a(kappa; lam) = 1 - sum_i R_i on 5 held-out kappa of [10,7,4,1] + the G8b table values (32/35, 464/315, 1)'


def test_a4a_formula():
    lam = (9, 7, 4, 2)
    for k in [(8, 5, 3), (9, 6, 2), (7, 5, 4)]:
        res, nu = cf.a_data(lam, lam, k)
        assert res['V4a'][0] == cf.a4a_of(k, lam), (k, res['V4a'][0], cf.a4a_of(k, lam))
    assert cf.a4a_of((5, 3, 1), (6, 4, 2, 0)) == Fr(128, 35) and cf.a4a_of((6, 4, 2), (6, 4, 2, 0)) == 0
    return 'a_4a(kappa; lam) = d + sum_i R_i (2 y_i - |kappa| + 4) on 3 held-out kappa of [9,7,4,2] (lam4 > 0) + G8b values (128/35, 0)'


def test_b_shifts():
    rng = random.Random(8)
    lam = (11, 7, 4, 2)
    n = 0
    for k, a, a4, b3, b4, nu, lp in cf.a_samples(lam, rng, nrand=1)[:6]:
        assert cf.b_of('V3', k, lam, nu, lp) == b3 and cf.b_of('V4a', k, lam, nu, lp) == b4, (k, nu, lp)
        n += 1
    # the example of G8b section 2.3: lam = (6,4,2,0), kappa = (5,3,0): b4a = 80
    res, nu = cf.a_data((6, 4, 2, 0), (6, 4, 2, 0), (5, 3, 0))
    assert res['V4a'][1] == 80 == cf.b_of('V4a', (5, 3, 0), (6, 4, 2, 0), nu, (6, 4, 2, 0))
    return 'b3, b4a closed forms: %d random samples of [11,7,4,2] + G8b example b4a = 80' % n


def test_gl3_chain_general_lamp():
    rng = random.Random(3)
    k = (6, 4, 1)
    R3 = gr.GLNRep.of(k)
    n = 0
    for nu in [(3, 4, 4), (4, 3, 4)]:
        lamp = (Fr(1, 3), Fr(-5, 7), Fr(9, 2))
        W = R3.wsp[nu]
        M = so.weight_space_mats(R3, {'V': so.reduced_op('V', lamp)}, W)['V']
        pos = {int(c): i for i, c in enumerate(W)}
        P, al, up, dn = cf.gl3_chain(k, lamp, nu)
        idx = [pos[R3.index[((nu[0],), (m, nu[0] + nu[1] - m), k)]] for m in P]
        for i in range(len(P)):
            assert M[idx[i]][idx[i]] == al[i]
            if i + 1 < len(P):
                assert M[idx[i + 1]][idx[i]] == up[i] and M[idx[i]][idx[i + 1]] == dn[i + 1]
            n += 1
    return 'general-lam\' gl(3) chain (u, d, alpha) == exact reduced operator on 2 weight spaces (%d path points, generic rational lam\')' % n


def test_closed_mats_vs_g3():
    lam = (6, 4, 2, 0)
    mus = [(10, 8, 5, 1), (9, 7, 6, 2)]
    for mu in mus:
        ww, eq, t1, t2 = cf.compare_mats(lam, lam, mu)
        assert ww >= 4 and all(eq.values()), (mu, eq)
    # generic non-dominant lam'
    lamp = (Fr(1, 3), Fr(5, 7), Fr(-2), Fr(4, 5))
    nu = (4, 3, 3, 2)
    mu = tuple(nu[i] + lamp[i] for i in range(4))
    ww, eq, t1, t2 = cf.compare_mats((6, 4, 2, 0), lamp, mu)
    assert all(eq.values())
    return 'closed_form_mats (V3, V4a) == G3 exactly on 2 weight spaces of [6,4,2,0]^2 and one with generic lam\' (w = %d)' % ww


def test_slices_closed():
    mp.mp.dps = 50
    lam = (6, 4, 2, 0)
    mu = (10, 7, 5, 2)
    cps = so.copies(lam, lam, mu, dps=50)
    res = cf.slices_closed(lam, lam, mu, 50, trips=None)
    R = gr.GLNRep.of(lam)
    assert len(res) == len(cps) > 1
    worst = mp.mpf(0)
    for cp, (trip, nz, r) in zip(cps, res):
        ref = {R.pats[i]: v for i, v in cp['slice'].items()}
        dv = max([abs(ref.get(k, 0) - nz.get(k, 0)) for k in set(ref) | set(nz)])
        worst = max(worst, dv)
    assert worst < mp.mpf(10) ** -30, worst
    return 'slices_closed == G3 copies on the sector (10,7,5,2) of [6,4,2,0]^2 (%d copies, max deviation %s)' % (len(cps), mp.nstr(worst, 3))


def test_word_groups():
    lam = (4, 2, 1, 0)
    cnt, seen = cf.words_product(lam, lam, names=('V3', 'V4a'))
    for (n, g, kappa), vals in seen.items():
        d = Fr(sum(lam) - sum(kappa))
        ref = {('V3', 'alpha'): Fr(1), ('V3', 'beta'): cf.a_of(kappa, lam) - 1, ('V4a', 'alpha'): d, ('V4a', 'beta'): cf.a4a_of(kappa, lam) - d}[(n, g)]
        assert vals == {ref}, (n, g, kappa, vals, ref)
    assert cnt[('V3', 'alpha', 'fit')] == cnt[('V3', 'alpha', 'pairs')] > 10
    return 'word groups of [4,2,1,0]^2: alpha = (1, d) V^(3), beta = (a - 1, a4a - d) V^(3) + scalar on all %d path blocks' % cnt[('V3', 'beta', 'pairs')]


def test_degenerate_and_e():
    tot = cf.degenerate_check([(8, 5, 3, 1)], nper=3)
    assert tot['bad'] == 0 and tot['tests'] > 20, dict(tot)
    rng = random.Random(1)
    lam, k = (9, 6, 3, 0), (8, 5, 2)
    s = cf.v4b_sample(lam, k, 5, 3, 3, (Fr(1, 2), Fr(3), Fr(-1, 3), Fr(2)))
    assert s['e'] == (1 - cf.a_of(k, lam)) / 8 == Fr(1, 105)
    return 'degenerate kappa (1x1 paths, R_i dropped): %d tests exact; V4b e = (1 - a)/8 = 1/105 at kappa = (8,5,2) of [9,6,3,0]' % tot['tests']


def run_all():
    tests = [test_a_formula, test_a4a_formula, test_b_shifts, test_gl3_chain_general_lamp, test_closed_mats_vs_g3, test_slices_closed, test_word_groups, test_degenerate_and_e]
    nfail = 0
    t00 = time.perf_counter()
    mp.mp.dps = 50
    for t in tests:
        t0 = time.perf_counter()
        try:
            msg = t()
            print('PASS  %-34s %.1fs  %s' % (t.__name__, time.perf_counter() - t0, msg), flush=True)
        except Exception as e:
            nfail += 1
            import traceback
            print('FAIL  %-34s %s' % (t.__name__, repr(e)), flush=True)
            traceback.print_exc()
    print('%d tests, %d failed (%.1f s)' % (len(tests), nfail, time.perf_counter() - t00))
    return nfail


if __name__ == '__main__':
    sys.exit(1 if run_all() else 0)
