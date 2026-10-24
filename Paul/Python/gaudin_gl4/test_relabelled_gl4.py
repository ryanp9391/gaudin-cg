"""test_relabelled_gl4.py -- exact regression tests for relabelled_gl4.py (WP-G8c); pytest-style functions AND runnable as `python -X utf8 test_relabelled_gl4.py` (< 60 s).
Everything is exact (Fractions / exact charpolys / sympy)."""
import os, sys, time
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import glN_reps as gr
import slice_engine_gl4 as se
import slice_ops_gl4 as so
import structure_gl4 as sg
import relabelled_gl4 as rg

L3 = (3, 2, 1, 0)
K4 = (5, 4, 2, 1)


def test_identity_equals_g3():
    """sigma = id: the matrices are G3's weight-space matrices, the restriction to M_mu is G3's V_red in the PRV basis, and the strings are the standard ones"""
    R = gr.GLNRep.of(L3)
    for mu in ((5, 4, 2, 1), (6, 3, 2, 1), (5, 3, 3, 1)):
        nu = tuple(mu[a] - L3[a] for a in range(4))
        W = R.wsp[nu]
        ref = so.weight_space_mats(R, {n: so.reduced_op(n, L3) for n in so.TRIPLE}, W)
        sm = rg.sigma_mats(L3, L3, mu, rg.ID)
        assert sm['mats'] == ref and sm['W'] == [int(c) for c in W]
        prv = se.prv_space(L3, L3, mu)
        basis, free = rg.kernel_basis(sm)
        assert basis == prv.basis and free == prv.free
        for n in so.TRIPLE:
            assert rg.restrict_to(sm['mats'][n], basis, free) == so.V_red(L3, L3, mu, n)
    return 'sigma = id reproduces G3 (weight-space matrices, PRV basis, V_red on M_mu) on 3 sectors of [3,2,1,0]^2'


def test_convention_check():
    """the joint spectrum (exact charpolys of V3, V4a, V4b, generic combination + commutation) of the relabelled triple on M^sigma equals that on M_mu for ALL 24 sigma; the wrong weight
    permutation fails for non-involutions"""
    res, m = rg.check_sector(L3, L3, K4)
    assert m == 4 and set(res.values()) == {'ok'}
    res2, _ = rg.check_sector(L3, L3, K4, conv='badw')
    for s, st in res2.items():
        if s == rg.inv(s):
            assert st == 'ok'
    assert any(st != 'ok' for s, st in res2.items() if s != rg.inv(s))
    res3, _ = rg.check_sector(L3, L3, (5, 3, 2, 2), sigmas=[(2, 3, 4, 1), (4, 3, 2, 1)])
    assert set(res3.values()) == {'ok'}
    return 'convention tau = sigma^-1, weight nu_sigma(a): 24/24 sigma ok on (5,4,2,1); the wrong weight permutation fails on %d non-involutions' % sum(st != 'ok' for s, st in res2.items() if s != rg.inv(s))


def test_w0_cut():
    """sigma = w0: E_34^{k3} -> E_21^{k3}; kernel = {m22 > nu'_1 - k3} = {m12 <= mu_3 - lam'_4} exactly (coordinate cut) on every sector of [3,2,1,0]^2 and [4,2,1,0]^2"""
    n = 0
    for lam in ((3, 2, 1, 0), (4, 2, 1, 0), (6, 4, 2, 0)):
        chk, bad = rg.w0_cut_check(lam, lam)
        assert chk >= 2 and bad == 0
        n += chk
        # same in the m12 form
        for mu, m in sg.sector_list(lam, lam, 'rs'):
            sm = rg.sigma_mats(lam, lam, mu, rg.W0)
            if sm is None:
                continue
            s = sm['strings'][2]
            if s['void']:
                continue
            assert s['zero'] == {j for j, c in enumerate(sm['coords']) if c[0] <= mu[2] - lam[3]}
    return 'w0: E_34 string is the coordinate cut m22 > nu_4 - k_3 (= m12 <= mu_3 - lam_4) on %d non-void sectors' % n


def test_cut_classification():
    """[4,2,1,0]^2 over all sectors: every predicted cut ({sigma1,sigma2} = {i,i+1}) is a coordinate subspace (exact); two active strings are never both non-trivially coordinate in one sigma"""
    lam = (4, 2, 1, 0)
    npred = 0
    two = 0
    for mu, m in sg.sector_list(lam, lam, 'rs'):
        R = gr.GLNRep.of(lam)
        nu = tuple(mu[a] - lam[a] for a in range(4))
        W = R.wsp.get(nu)
        if W is None or len(W) < 2:
            continue
        A = rg.active_strings(R, lam, nu)
        good_all = []
        for s in rg.SIGMAS:
            sm = rg.sigma_mats(lam, lam, mu, s)
            for i in A:
                st = sm['strings'][i - 1]
                if {s[0], s[1]} == {i, i + 1}:
                    assert st['coord'] and not st['void']
                    npred += 1
            if len(A) >= 2 and all(sm['strings'][i - 1]['coord'] and sm['strings'][i - 1]['nullity'] not in (0, sm['w']) for i in A):
                good_all.append(s)
        if len(A) >= 2:
            two += 1
            assert good_all == []
    return 'predicted cuts coordinate in %d (sector, string, sigma) cases; %d sectors with >= 2 active strings, never all coordinate in one sigma' % (npred, two)


def test_k4_example():
    """(5,4,2,1) of [3,2,1,0]^2: exact lattice, matrices, graph K4, charpolys, the 2+2 rational block structure"""
    sm = rg.sigma_mats(L3, L3, K4, rg.ID)
    assert sm['coords'] == [(2, 2, 2), (3, 3, 1), (2, 3, 2), (3, 3, 2)]
    V3 = sm['mats']['V3']
    assert V3[0] == [Fr(381, 8), Fr(-9, 8), Fr(-11, 8), Fr(3, 8)] and V3[3][3] == Fr(261, 4)
    g = sg.graph_stats(V3)
    assert g['edges'] == 6 and g['maxdeg'] == 3 and not g['paths'] and g['bw'] == 3
    import sympy as sp
    x = sp.symbols('x')
    cp = sp.factor(rg.spmat(V3).charpoly(x).as_expr())
    assert cp == sp.factor((x - 57) ** 2 * (x ** 2 - 114 * x + 3141))
    blocks = rg.galois_blocks([sm['mats'][n] for n in rg.TRIPLE])
    assert sorted(b['deg'] for b in blocks) == [2, 2]
    mm = [sm['mats'][n] for n in rg.TRIPLE]
    for b in blocks:
        Rs = rg.block_matrices(b['basis'], mm)
        if Rs[0] == sp.eye(2) * 57:
            assert Rs[1] + Rs[2] == 72 * sp.eye(2)
    # no sigma makes V3 (or the triple) a path
    for s in rg.SIGMAS:
        s2 = rg.sigma_mats(L3, L3, K4, s)
        assert not rg.graph(s2['mats']['V3'])['paths'] and not rg.graph(rg.union_matrix(s2['mats']))['paths']
    return 'K4 example: lattice (2,2,2),(3,3,1),(2,3,2),(3,3,2); V3 charpoly (x-57)^2(x^2-114x+3141); rational blocks 2+2; no sigma gives a path'


def test_bandwidth_statistic():
    """[3,2,1,0]^2: 8 weight spaces of dimension >= 2; sigma and sigma o (12) have identical graph statistics; the number of V3-paths at id and its sigma-invariance"""
    recs = rg.compute_records([(L3, L3)], cache=False, verbose=False)
    assert len(recs) == 8
    for r in recs:
        for s in rg.SIGMAS:
            a, b = r['per'][s], r['per'][rg.swap12(s)]
            assert (a['g3'], a['gu']) == (b['g3'], b['gu'])
    npath = sum(r['per'][rg.ID]['g3']['paths'] for r in recs)
    some = sum(any(r['per'][s]['g3']['paths'] for s in rg.SIGMAS) for r in recs)
    assert npath == some
    # the stencil is the 20 off-diagonal points in every sigma
    prod = {(a, b, c) for a in (-1, 0, 1) for (b, c) in [(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]} - {(0, 0, 0)}
    for s in rg.SIGMAS:
        U = set()
        for r in recs:
            U |= set(r['per'][s]['stencil'])
        assert U <= prod
    return '[3,2,1,0]^2: %d weight spaces, V3 path at id in %d, at some sigma in %d; stencil inside the 20 off-diagonal points for all 24 sigma' % (len(recs), npath, some)


def run_all():
    t00 = time.perf_counter()
    for name, f in sorted((k, v) for k, v in globals().items() if k.startswith('test_') and callable(v)):
        t0 = time.perf_counter()
        msg = f()
        print('PASS %-30s %5.1fs  %s' % (name, time.perf_counter() - t0, msg), flush=True)
    print('all passed in %.1f s' % (time.perf_counter() - t00))


if __name__ == '__main__':
    run_all()
