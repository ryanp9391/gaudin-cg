"""test_gaudin_gl4_ops.py -- quick regression tests for gaudin_gl4_ops.py (pytest-style functions AND runnable as `python -X utf8 test_gaudin_gl4_ops.py`).
All checks are exact (Fractions / exact polynomial identities); the heavy runs (full commutativity tables, all sectors of all products) are the CLI stages of the module."""
import os, sys, time, itertools
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gaudin_gl4_ops as g


def test_convention_gl3_eq4():
    S = g.cdet_S(3)
    ref = g.eq4_gl3()
    for m in (1, 2, 3):
        assert S[m] == ref[m], m
    return 'S_1, S_2, S_3 of the cdet engine = eq. (4) of the gl(3) notes as noncommutative polynomials'


def test_homogeneity_and_counts():
    for N in (3, 4):
        L = g.laurent(N)
        for k in range(1, N + 1):
            for j in range(0, k + 1):
                pows = {p for pw in L[k].get(j, {}).values() for p in pw}
                assert pows <= {j - k}
    L = g.laurent(4)
    cnt = {(k, j): len(L[k][j]) for k in (2, 3, 4) for j in range(1, k + 1)}
    assert cnt[(3, 2)] == 88 and cnt[(4, 3)] == 152 and cnt[(4, 2)] == 324 and cnt[(4, 1)] == 412
    return 'every S_{k,j} homogeneous of degree j-k in t (N = 3, 4); term counts V3 88, V4a 152, V4b 324, S_{4,1} 412'


def test_classical_laurent():
    for N in (3, 4):
        bad, n = g.classical_laurent_check(N, seed=3)
        assert bad == 0
    return 'Laurent expansion (stage B) = sympy series with commuting random integer matrices, N = 3, 4'


def test_residue_relations():
    from math import comb
    for N in (3, 4):
        LU = g.laurent_up(N)
        for k in range(2, N + 1):
            for j in range(1, k):
                tot = {}
                for i in range(1, j + 1):
                    g.up_add(tot, LU[(k, i)][0], comb(j - 1, i - 1))
                g.up_add(tot, g.up_swap(LU[(k, j)][0]), (-1) ** (j - k))
                assert not g.up_pbw(g.up_clean(tot))
    return 'S_k = O(u^-k) at infinity: all residue relations hold exactly in U(gl_N)(x)U(gl_N), N = 3, 4 (checks the theta_2 swap rule)'


def test_gl3_facts():
    se_path = os.path.join(HERE, '..', 'gaudin_sl3')
    sys.path.insert(0, se_path)
    import slice_engine as se
    lam = (2, 1, 0)
    R = g.Rep.of(lam)
    t = Fr(31, 35)
    ops = g.build_ops(R, R, 3, t)
    mats = {k: o.matrix_exact() for k, o in ops.items()}
    keys = list(mats)
    for a, b in itertools.combinations(keys, 2):
        assert g.nnz_exact(g.comm_exact(mats[a], mats[b])) == 0
    LU = g.laurent_up(3)
    assert g.site_scalar(R, LU[(2, 2)][0]) == lam[1] * (lam[0] + 1)
    V = {}
    for coef, ij1, ij2, ij3 in se.VTERMS:
        V[((ij1, ij2), (ij3,))] = V.get(((ij1, ij2), (ij3,)), 0) + coef
    n1, n2, Om = g.inv_n(1, 3), g.inv_n(2, 3), g.inv_Omega(3)
    full = g.up_lin((1, V), (Fr(1, 2), g.up_lin((1, g.up_mul(n1, n2)), (-1, Om))))
    assert g.nnz_exact(g._sub(g.Op2(LU[(3, 2)][0], R, R).matrix_exact(), g.Op2(full, R, R).matrix_exact())) == 0
    return 'gl(3) [2,1,0]^2: all S_{k,j} commute (exact), S_{2,2} = lam2 (lam1+1), t S_{3,2} = V + (n1 n2 - Omega)/2 on the whole space'


def test_gl4_commutativity_small():
    for l1, l2 in [((1, 0, 0, 0), (1, 1, 0, 0)), ((1, 1, 0, 0), (2, 1, 0, 0))]:
        R1, R2 = g.Rep.of(l1), g.Rep.of(l2)
        ops = g.build_ops(R1, R2, 4, Fr(-3, 7))
        mats = {k: o.matrix_exact() for k, o in ops.items()}
        nz = [k for k in mats if g.nnz_exact(mats[k])]
        tot = g.tot_generators_exact(R1, R2)
        for a, b in itertools.combinations(nz, 2):
            assert g.nnz_exact(g.comm_exact(mats[a], mats[b])) == 0
        for k in nz:
            for ij, E in tot.items():
                assert g.nnz_exact(g.comm_exact(mats[k], E)) == 0
    return 'N = 4: all 18 Laurent coefficients commute pairwise and with E^tot on [1,0,0,0]x[1,1,0,0] (dim 24) and [1,1,0,0]x[2,1,0,0] (dim 120), t = -3/7'


def test_double_poles():
    LU = g.laurent_up(4)
    for lam in [(1, 0, 0, 0), (2, 1, 0, 0), (2, 1, 1, 0), (3, 2, 1, 0), (2, 2, 1, 1), (4, 2, 1, 0)]:
        R = g.Rep.of(lam)
        form = lam[1] * (lam[0] + 1) + lam[2] * (lam[0] + lam[1] + 2) + lam[3] * (lam[0] + lam[1] + lam[2] + 3)
        assert g.site_scalar(R, LU[(2, 2)][0]) == form == g.newton_Smm(lam, 2)
        assert g.site_scalar(R, LU[(3, 3)][0]) == g.newton_Smm(lam, 3)
        assert g.site_scalar(R, LU[(4, 4)][0]) == g.newton_Smm(lam, 4)
    return 'S_{2,2} = lam2(lam1+1) + lam3(lam1+lam2+2) + lam4(lam1+lam2+lam3+3); S_{3,3}, S_{4,4} = Capelli/falling-factorial closed forms (6 reps)'


def test_pbw_identities():
    N = 4
    LU = g.laurent_up(4)
    V3, V4a = LU[(3, 2)][0], LU[(4, 3)][0]
    n1, n2, Om = g.inv_n(1, N), g.inv_n(2, N), g.inv_Omega(N)
    eta = g.up_lin((1, g.up_mul(n1, n2)), (-1, Om))
    D21, D31 = g.antisym_sum_subsets(N, 2, 1), g.antisym_form(N, 3, 1)
    assert not g.up_pbw(g.up_lin((1, V3), (-1, D21), (-1, eta)))
    assert not g.up_pbw(g.up_lin((1, V4a), (-1, D31), (-1, V3), (Fr(1, 3), eta)))
    T112 = g.inv_T((1, 1, 2), N)
    C21 = g.inv_C2(1, N)
    rem = g.up_lin((Fr(1, 2), g.up_mul(g.up_mul(n1, n1), n2)), (Fr(-1, 2), g.up_mul(C21, n2)), (-1, g.up_mul(n1, Om)), (Fr(1, 2), g.up_mul(n1, n2)), (1, Om))
    assert not g.up_pbw(g.up_lin((1, V3), (-1, T112), (-1, rem)))
    # the delta form of the brief: 1/2 sum_{sigma,tau in S4} sgn sgn e1 e1 e2 delta_{sigma4 tau4} = sum_I D(2,1)_I
    delta = {}
    for s in itertools.permutations(range(1, 5)):
        for t in itertools.permutations(range(1, 5)):
            if s[3] != t[3]:
                continue
            k = (((s[0], t[0]), (s[1], t[1])), ((s[2], t[2]),))
            delta[k] = delta.get(k, 0) + Fr(g.perm_sign(s) * g.perm_sign(t), 2)
    assert g.up_clean(delta) == D21
    out = g.cubic_general_N(Ns=(3, 4, 5), verbose=False)
    assert [out[N][0] for N in (3, 4, 5)] == [Fr(1, 2), Fr(1), Fr(3, 2)]
    return 'V3 = sum_I D(2,1)_I + (N-2)/2 (n1 n2 - Omega) (N = 3, 4, 5) = T(1,1,2) + [1/2 n1^2 n2 + 1/2 n1 n2 - n1 Omega - 1/2 n2 C2_1 + Omega] (N = 4);  V4a = D(3,1) + V3 - (n1 n2 - Omega)/3  (exact, PBW)'


def test_triple_sectors():
    out = []
    for (l1, l2, mu, v3_simple) in [((2, 1, 1, 0), (2, 1, 1, 0), (3, 2, 2, 1), True), ((3, 2, 1, 0), (3, 2, 1, 0), (5, 4, 2, 1), False), ((3, 2, 1, 0), (3, 2, 1, 0), (4, 4, 2, 2), False),
                                    ((3, 2, 1, 0), (3, 2, 1, 0), (5, 3, 3, 1), True)]:
        R1, R2 = g.Rep.of(l1), g.Rep.of(l2)
        m = g.racah_speiser(l1, l2)[mu]
        res = g.analyse_sector(R1, R2, mu, m, g.triple_ops(R1, R2), others=g.build_ops(R1, R2, 4, Fr(31, 35)))
        assert res['faithful_viol'] == 0 and res['comm_bad'] == 0 and res['selfadj_bad'] == 0
        assert res['joint_real_simple'] and res['V3_real_simple'] == v3_simple
        assert res['others_noncomm'] == 0 and res['others_offdiag'] < 1e-30
        out.append('%s m=%d' % (mu, m))
    return 'triple commutes exactly, is self-adjoint, joint spectrum simple (exact certificates); V3 alone degenerate at (5,4,2,1) [m=4] and (4,4,2,2) [m=2] of [3,2,1,0]^2: ' + ', '.join(out)


def test_terms_roundtrip():
    import tempfile
    p = os.path.join(tempfile.mkdtemp(), 'ops.txt')
    g.write_terms(p)
    loaded = g.load_terms(p)
    gen = g.terms_for_file()
    assert set(loaded) == set(gen)
    for name in gen:
        assert loaded[name] == {ww: {pp: Fr(c) for pp, c in pw.items()} for ww, pw in gen[name].items()}
    return 'term-list file: write/load round trip identical (%d operators)' % len(loaded)


def test_scaffold_vs_glN_reps():
    try:
        import glN_reps as G
    except Exception as e:                                                  # G1's module may not be importable at this moment: not a failure of this package
        return 'glN_reps not importable (%s): skipped' % type(e).__name__
    for lam in [(1, 0, 0, 0), (1, 1, 0, 0), (2, 1, 0, 0), (2, 1, 1, 0), (3, 2, 1, 0)]:
        R, S = g.Rep.of(lam), G.GLNRep.of(lam)
        assert R.pats == S.pats and R.wt == S.wt and R.top == S.top
        for ij in itertools.product(range(1, 5), repeat=2):
            assert R.E[ij] == S.E[ij], (lam, ij)
        assert list(R.g) == list(S.g)
    return 'private scaffold == glN_reps.GLNRep (patterns, weights, all 16 generators, metric) on 5 gl(4) irreps'


TESTS = [test_convention_gl3_eq4, test_homogeneity_and_counts, test_classical_laurent, test_residue_relations, test_gl3_facts, test_gl4_commutativity_small, test_double_poles,
         test_pbw_identities, test_triple_sectors, test_terms_roundtrip, test_scaffold_vs_glN_reps]


def main():
    nfail = 0
    for f in TESTS:
        t0 = time.perf_counter()
        try:
            msg = f()
            print('  PASS %-32s (%5.2fs)  %s' % (f.__name__, time.perf_counter() - t0, msg), flush=True)
        except Exception as e:
            nfail += 1
            import traceback
            print('  FAIL %-32s (%5.2fs)  %r' % (f.__name__, time.perf_counter() - t0, e), flush=True)
            traceback.print_exc()
    print('%d tests, %d failed' % (len(TESTS), nfail))
    return nfail


if __name__ == '__main__':
    sys.exit(1 if main() else 0)
