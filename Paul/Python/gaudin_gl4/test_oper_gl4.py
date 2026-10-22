"""test_oper_gl4.py -- few-second regression test of oper_gl4.py (WP-G5):  python -X utf8 test_oper_gl4.py"""
import sys, time, traceback
from fractions import Fraction as Fr
import mpmath as mp
import sympy as sp
import oper_gl4 as O

mp.mp.dps = 60
TESTS = []


def test(f):
    TESTS.append(f)
    return f


@test
def test_central_characters():
    # T_k(nu) = (-1)^k a_{4-k}: S_{2,2} = e2(lam) + sum (k-1) lam_k of G2; trace-form Casimir fits
    assert O.Tchar((2, 1, 0, 0), 2) == 3 and O.Tchar((1, 1, 0, 0), 2) == 2 and O.Tchar((4, 3, 2, 1), 2) == 45
    assert O.Tchar((3, 2, 1, 0), 3) == 15 and O.Tchar((2, 2, 1, 1), 4) == 40            # G2: S_{3,3} = 15 on [3,2,1,0], S_{4,4} = 40 on [2,2,1,1]
    assert O.casimirs((2, 1, 0, 0)) == (3, 12, 12, 48)
    # closed forms in the shifted weights y_c = nu_c + 4 - c with elementary symmetric functions e_k(y): T_2 = e_2 - 3 e_1 + 7, T_3 = e_3 - e_2 + e_1 - 1, T_4 = e_4  (all 70 partitions with nu_1 <= 4)
    import itertools as it
    import math
    for nu in O.partitions4(4):
        y = O.shifted(nu)
        e = [Fr(sum(math.prod(cc) for cc in it.combinations(y, k))) for k in range(5)]
        assert O.Tchar(nu, 1) == e[1] - 6 == sum(nu)
        assert O.Tchar(nu, 2) == e[2] - 3 * e[1] + 7 and O.Tchar(nu, 3) == e[3] - e[2] + e[1] - 1 and O.Tchar(nu, 4) == e[4]
    for nu in [(3, 1, 1, 0), (4, 2, 1, 0), (5, 2, 0, 0)]:
        n, c2, c3, c4 = O.casimirs(nu)
        assert O.Tchar(nu, 2) == -c2 / 2 + Fr(n * n, 2) + Fr(3 * n, 2)
        assert O.Tchar(nu, 3) == c3 / 3 - Fr(n) * c2 / 2 + Fr(n ** 3, 6) - c2 / 3 + Fr(5 * n * n, 6) + 2 * n


@test
def test_dictionary_sector():
    prod = O.Product.of((2, 1, 1, 0), (2, 1, 1, 0))
    nb = 0
    for t in (Fr(1), Fr(31, 35)):
        for mu, m in prod.irreps():
            sec = prod.sector(mu)
            bad, found = O.check_dictionary_sector(sec, t)
            nb += len(bad)
    assert nb == 0, nb
    # the m = 2 sector: the remainders are scalar matrices (not only on the diagonal)
    sec = prod.sector((3, 2, 2, 1))
    assert sec.m == 2 and O.check_dictionary_sector(sec)[0] == []


@test
def test_oper_structure():
    prod = O.Product.of((2, 1, 0, 0), (2, 1, 0, 0))
    for mu in [(4, 2, 0, 0), (3, 2, 1, 0), (2, 2, 2, 0)]:
        op = O.Oper((2, 1, 0, 0), (2, 1, 0, 0), mu, Fr(31, 35))
        assert all(a == O.AFF0 for a in op.residue_defects().values())
        a = O.falling_coeffs(mu)
        assert all(op.I0_poly_coeffs()[i] == O.aff_const(a[4 - i]) for i in range(5))            # indicial polynomial at infinity = prod (n - (mu_c + 4 - c)), v-independent
        for s, nu in ((1, (2, 1, 0, 0)), (2, (2, 1, 0, 0))):
            a2 = O.falling_coeffs(nu)
            assert all(op.local_indicial(s)[i] == O.aff_const(a2[4 - i]) for i in range(5))      # local exponents {lam_c + 4 - c}
        sec = prod.sector(mu)
        for v in sec.exact_triples() or []:
            nd, sols = O.poly_solutions_exact(op, v)
            assert nd == [4, 3, 2, 1], nd                                                         # null dimensions #{k : d_k <= d_j}
            W = O.p_trim(O.wronskian_poly(sols, Fr(0)), Fr(0))
            T = O.wronskian_target(op, Fr(1))
            assert len(W) == len(T) and all(W[i] == (W[-1] / T[-1]) * T[i] for i in range(len(W)))


@test
def test_copy_mp():
    prod = O.Product.of((2, 1, 1, 0), (2, 1, 1, 0))
    sec = prod.sector((3, 2, 2, 1))                    # m = 2, irrational triples
    op = O.Oper((2, 1, 1, 0), (2, 1, 1, 0), (3, 2, 2, 1), Fr(1))
    js = sec.joint(60)
    assert len(js['trip']) == 2
    for trip in js['trip']:
        r = O.copy_check(op, list(trip), None, 60)
        assert r['ok_null'] and r['wr_ok'] and r['gap'] > mp.mpf(10) ** 40
        assert r['orders1'] == sorted(op.y1) and r['orders2'] == sorted(op.y2)
        ng = O.negative_tests(op, trip, 60)
        assert ng['true'][0] < mp.mpf(10) ** (-50) and all(x[0] > mp.mpf(10) ** (-20) for x in ng['perturbed']) and ng['random'][0] > mp.mpf(10) ** (-20)


@test
def test_exact_solve():
    for (l1, l2, mu, m) in [((2, 1, 0, 0), (2, 1, 0, 0), (3, 2, 1, 0), 2), ((2, 1, 0, 0), (2, 1, 0, 0), (2, 2, 2, 0), 1)]:
        r = O.solve_sector(l1, l2, mu, 60, subsets=True)
        assert r['all_vanish_on_G2'] and r['full']['dim'] == 0 and r['full']['cnt'] == m and r['distinct'] == m and r['jac_ranks'] == [3] * m
        assert r['dev'] < mp.mpf(10) ** (-50)
        assert r['subsets'][(4, 3, 2, 1)]['equalsK']
    r = O.solve_sector((2, 1, 0, 0), (2, 1, 0, 0), (2, 2, 2, 0), 60, subsets=True)
    s = r['subsets']
    assert s[(4,)]['dim'] == 1 and s[(2,)]['cnt'] == 5 and not s[(2,)]['equalsK'] and s[(3, 2, 1)]['equalsK']         # pencil-bug analogue: j = 2 alone has 5 solutions, {4,3} has 3
    assert s[(4, 3)]['cnt'] == 3 and s[(4, 2)]['equalsK']


@test
def test_mu4_zero():
    # mu_4 = 0 (lam_4 = lam'_4 = 0): the constant polynomial is the degree-d_4 = 0 solution, D.1 = s_4 = 0, which forces v4a = v4b = 0 exactly (G2's triples are (v3, 0, 0) there)
    n = 0
    for l1, l2 in O.REQUIRED_PRODUCTS[:5]:
        prod = O.Product.of(l1, l2)
        for mu, m in prod.irreps():
            if mu[3] == 0:
                for row in prod.sector(mu).joint(60)['trip']:
                    n += 1
                    assert abs(row[1]) < mp.mpf(10) ** (-40) and abs(row[2]) < mp.mpf(10) ** (-40), (l1, l2, mu, row)
    assert n >= 15, n


@test
def test_wb():
    coeffs, body = O.transcribe_gaudinODE()
    z, S, P, w1, w2, C2, C4 = (O.WB_SYMS[n] for n in ('z', 'S', 'P', 'w1', 'w2', 'C2', 'C4'))
    assert sp.expand(coeffs[4] - 8 * (z - w1) ** 2 * (z - w2) ** 2) == 0
    assert O.wb_Cas2((2, 2, 0, 0)) == 16 and O.wb_Cas4((2, 2, 0, 0)) == 256 and O.wb_degQ(1, 1, (2, 1, 1, 0)) == [0, 2, 3, 5]
    # (a) degrees for (S,P) = (1,1)
    for mu in O.wb_irreps(1, 1):
        subs = {S: 1, P: 1, w1: 0, w2: 1, C2: O.wb_Cas2(mu), C4: O.wb_Cas4(mu)}
        jumps, tot = O.wb_poly_solution_degrees(coeffs, subs, 12)
        assert jumps == sorted(O.wb_degQ(1, 1, mu)) and tot == 4
    # (b) formal adjoint + gauge p = (S, P): unique solution = G2's triple, and the closed-form dictionary
    W1, W2 = sp.Rational(2, 7), sp.Rational(-3, 5)
    mu = (2, 1, 1, 0)
    C2v, C4v = O.wb_Cas2(mu), O.wb_Cas4(mu)
    unk, sol = O.wb_match(1, 1, 1, 1, W1, W2, (C2v, C4v), coeffs=coeffs)
    sol = list(sol)
    assert len(sol) == 1
    g2 = O.Product.of((1, 1, 0, 0), (1, 1, 0, 0)).sector(mu).exact_triples()[0]
    assert [Fr(int(x.p), int(x.q)) for x in sol[0][:3]] == list(g2) and [Fr(int(x.p), int(x.q)) for x in sol[0][3:]] == [O.Tchar(mu, k) for k in (2, 3, 4)]
    cf = O.wb_dictionary(1, 1, C2v, C4v)
    assert Fr(int(cf['v3'].p), int(cf['v3'].q)) == g2[0] and cf['v4a'] == g2[1] == 0 and Fr(int(cf['v4b'].p), int(cf['v4b'].q)) == g2[2]
    unk, sol = O.wb_match(1, 1, 0, 1, W1, W2, (C2v, C4v), coeffs=coeffs)          # wrong gauge: inconsistent
    assert sol == sp.EmptySet


def main():
    t0 = time.perf_counter()
    failed = 0
    for f in TESTS:
        t1 = time.perf_counter()
        try:
            f()
            print('  ok    %-28s %.2fs' % (f.__name__, time.perf_counter() - t1))
        except Exception:
            failed += 1
            print('  FAIL  %-28s' % f.__name__)
            traceback.print_exc()
    print('%d tests, %d failed  (%.1fs)' % (len(TESTS), failed, time.perf_counter() - t0))
    return failed


if __name__ == '__main__':
    sys.exit(1 if main() else 0)
