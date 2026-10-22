"""test_glN_reps.py -- tests of glN_reps.py (WP-G1, task G1.1).  pytest-style functions AND runnable:  python -X utf8 test_glN_reps.py

(i)   N = 3 regression: all nine E_ij equal slice_engine.Rep.of(lam).E entry for entry (diagonal: Rep.wt) and gl3_reps.GL3Rep(lam).generator(i,j) (adjacent
      ones and the diagonal); pats / wt / wsp / top / metric / Racah-Speiser identical to slice_engine
(ii)  N = 4: all 256 commutation relations [E_ij,E_kl] = d_jk E_il - d_li E_kj exact on nine irreps; N = 5: all 625 on two; path independence of E_14, E_41
(iii) dim = Weyl formula, HW state last with weight lam, weight multiset invariant under S_N, sum_w m(w) = dim
(iv)  metric: E_ij^T G = G E_ji exactly on every matrix element of every off-diagonal pair (0 violations)
(v)   Casimir C_2 = Sum_ij E_ij E_ji is the scalar Sum_i lam_i (lam_i + N + 1 - 2i) on every basis vector
(vi)  Racah-Speiser: Sum m dim = dim dim' on every product of G1.2/G1.3, [S,S,0,0]^2 (S = 1..10), and the Munich SU(4) table (hw_timings.txt) reproduced
      for S = 1..8 (irreps, dimensions, multiplicities; Munich's labels are ours minus mu_4)
(vii) the GT basis is Munich's: pattern order = Munich's Q order, orthonormal lowering matrix elements = clebsch::pattern::lowering_coeff exactly, trivial sign gauge
"""
import os, sys, re, time, itertools, traceback
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.normpath(os.path.join(HERE, '..', 'gaudin_sl3')))
import numpy as np
import glN_reps as gr
from glN_reps import GLNRep

INFO = {}


def _info(name, msg):
    INFO[name] = msg


MUNICH = os.path.normpath(os.path.join(HERE, '..', '..', 'Cpp', 'AlexCGC', 'results', 'hw_timings.txt'))


# ---------------------------------------------------------------------------------------- (i) N = 3 regression
def test_n3_regression():
    import slice_engine as se
    from gl3_reps import GL3Rep
    nent = 0
    for lam in [(2, 1, 0), (3, 1, 0), (3, 2, 0), (4, 2, 0)]:
        R = GLNRep.of(lam)
        S = se.Rep.of(lam)
        assert R.pats == S.pats and R.dim == S.dim and R.wt == S.wt and R.top == S.top, lam
        assert {w: list(v) for w, v in R.wsp.items()} == {w: list(v) for w, v in S.wsp.items()}
        assert R.g == S.g, lam
        for ij, cols in S.E.items():                      # the six off-diagonal matrices, entry for entry (Fractions)
            assert R.E[ij] == cols, (lam, ij)
            nent += sum(len(c) for c in cols)
        for k in (1, 2, 3):                               # the diagonal, from Rep.wt
            assert R.E[(k, k)] == [({c: Fr(w[k - 1])} if w[k - 1] != 0 else {}) for c, w in enumerate(S.wt)], (lam, k)
        G = GL3Rep(lam)
        for i in (1, 2, 3):                               # GL3Rep.generator: adjacent + diagonal
            for j in (1, 2, 3):
                if abs(i - j) > 1:
                    continue
                M = G.generator(i, j)
                cols = R.E[(i, j)]
                for c in range(R.dim):
                    for r in range(R.dim):
                        v = M[r, c]
                        ref = Fr(int(v.p), int(v.q)) if v != 0 else 0
                        assert cols[c].get(r, 0) == ref, (lam, i, j, r, c)
        assert R.check_metric() == 0
        assert R.check_metric() == S.check_metric() == 0
    # Racah-Speiser: our convention (sum over the weights of the second factor) vs slice_engine's (weights of the first factor), incl. lam != lamp
    for lam, lamp in [((2, 1, 0), (2, 1, 0)), ((3, 2, 0), (3, 1, 0)), ((3, 1, 0), (3, 2, 0)), ((4, 2, 0), (4, 2, 0)), ((4, 1, 1), (2, 2, 0)), ((6, 3, 0), (6, 3, 0))]:
        assert gr.racah_speiser(lam, lamp) == se.racah_speiser(lam, lamp), (lam, lamp)
        assert gr.racah_speiser(lam, lamp) == gr.racah_speiser(lamp, lam), (lam, lamp)
    _info('n3_regression', '4 irreps: all nine E_ij identical to slice_engine.Rep (%d off-diagonal entries) and to GL3Rep (adjacent+diagonal), metric identical, '
          'Racah-Speiser identical to slice_engine on 6 products incl. lam != lam\'' % nent)


# ---------------------------------------------------------------------------------------- (ii) commutation relations
def _commutators_ok(R):
    """all N^4 relations [E_ij, E_kl] = d_jk E_il - d_li E_kj, exactly, on every basis vector; returns the number of relations checked"""
    N, n = R.N, R.dim
    cnt = 0
    for i, j, k, l in itertools.product(range(1, N + 1), repeat=4):
        for c in range(n):
            e = {c: Fr(1)}
            lhs = {}
            for r, v in R.apply((i, j), R.apply((k, l), e)).items():
                lhs[r] = lhs.get(r, 0) + v
            for r, v in R.apply((k, l), R.apply((i, j), e)).items():
                lhs[r] = lhs.get(r, 0) - v
            rhs = {}
            if j == k:
                for r, v in R.apply((i, l), e).items():
                    rhs[r] = rhs.get(r, 0) + v
            if l == i:
                for r, v in R.apply((k, j), e).items():
                    rhs[r] = rhs.get(r, 0) - v
            lhs = {r: v for r, v in lhs.items() if v != 0}
            rhs = {r: v for r, v in rhs.items() if v != 0}
            assert lhs == rhs, ('commutation relation violated', R.lam, (i, j, k, l), c)
        cnt += 1
    return cnt


def test_commutation_n4():
    lams = [(1, 0, 0, 0), (1, 1, 0, 0), (2, 0, 0, 0), (2, 1, 0, 0), (2, 1, 1, 0), (2, 2, 0, 0), (3, 2, 1, 0), (1, 1, 1, 0), (3, 3, 0, 0)]
    tot = 0
    for lam in lams:
        R = GLNRep.of(lam)
        tot += _commutators_ok(R)
    _info('commutation_n4', 'N = 4: all 256 relations exact on %d irreps (%d relations x all basis vectors, 0 violations): %s' % (len(lams), tot, [list(l) for l in lams]))


def test_commutation_n5():
    lams = [(1, 0, 0, 0, 0), (2, 1, 0, 0, 0)]
    tot = 0
    for lam in lams:
        tot += _commutators_ok(GLNRep.of(lam))
    _info('commutation_n5', 'N = 5: all 625 relations exact on %s (%d relations, 0 violations)' % ([list(l) for l in lams], tot))


def test_path_independence():
    """E_14 = [E_13,E_34] (built) equals [E_12,E_24] and E_41 = [E_43,E_31] (built) equals [E_42,E_21]; same for N = 5 (E_15 via three paths)"""
    for lam in [(2, 1, 0, 0), (2, 1, 1, 0), (3, 2, 1, 0)]:
        R = GLNRep.of(lam)
        E = R.E
        pairs = [((1, 4), [((1, 3), (3, 4)), ((1, 2), (2, 4))]), ((4, 1), [((4, 3), (3, 1)), ((4, 2), (2, 1))]),
                 ((1, 3), [((1, 2), (2, 3))]), ((3, 1), [((3, 2), (2, 1))]), ((2, 4), [((2, 3), (3, 4))]), ((4, 2), [((4, 3), (3, 2))])]
        for tgt, paths in pairs:
            ref = E[tgt]
            for (a, b) in paths:
                alt = gr._sub(gr._compose(E[a], E[b]), gr._compose(E[b], E[a]))
                assert alt == ref, (lam, tgt, a, b)
    R = GLNRep.of((2, 1, 0, 0, 0))
    E = R.E
    for a, b in [((1, 4), (4, 5)), ((1, 2), (2, 5)), ((1, 3), (3, 5))]:
        assert gr._sub(gr._compose(E[a], E[b]), gr._compose(E[b], E[a])) == E[(1, 5)]
    _info('path_independence', 'E_14, E_41, E_13, E_31, E_24, E_42 (N = 4, three irreps) and E_15 (N = 5) are the same matrix along every commutator path')


# ---------------------------------------------------------------------------------------- (iii) dims, weights
def test_dims_weights_hw():
    lams = [(1, 0, 0, 0), (1, 1, 0, 0), (2, 1, 0, 0), (2, 1, 1, 0), (2, 2, 0, 0), (3, 2, 1, 0), (3, 3, 0, 0), (4, 4, 0, 0), (5, 3, 1, 0),
            (2, 1, 0), (4, 2, 0), (6, 3, 0), (2, 1, 0, 0, 0), (2, 2, 1, 0, 0), (3, 2, 1)]
    for lam in lams:
        R = GLNRep.of(lam)
        N = len(lam)
        assert R.dim == gr.weyl_dim(lam), lam
        assert R.top == R.dim - 1 and R.wt[R.top] == tuple(lam), lam
        assert R.pats[R.top] == tuple(tuple(lam[:k]) for k in range(1, N + 1))
        wm = gr.weight_mults(lam)
        assert sum(wm.values()) == R.dim
        assert dict(wm) == {w: len(v) for w, v in R.wsp.items()}
        for w, m in wm.items():                                         # S_N invariance of the weight multiset
            for p in itertools.permutations(range(N)):
                assert wm.get(tuple(w[i] for i in p), 0) == m, (lam, w, p)
        assert all(sum(w) == sum(lam) for w in wm)
    _info('dims_weights_hw', '%d irreps (N = 3, 4, 5): dim = Weyl, HW pattern last with weight lam, weight multiset S_N-symmetric, sum_w m(w) = dim; '
          'dims e.g. %s' % (len(lams), {str(list(l)): gr.weyl_dim(l) for l in [(2, 1, 0, 0), (3, 2, 1, 0), (4, 4, 0, 0), (10, 10, 0, 0)]}))


# ---------------------------------------------------------------------------------------- (iv) metric
def test_metric():
    lams = [(1, 0, 0, 0), (1, 1, 0, 0), (2, 1, 0, 0), (2, 1, 1, 0), (2, 2, 0, 0), (3, 2, 1, 0), (3, 3, 0, 0), (4, 4, 0, 0), (2, 1, 0, 0, 0), (2, 1, 0), (4, 2, 0)]
    for lam in lams:
        R = GLNRep.of(lam)
        assert R.check_metric() == 0, lam
        g = R.g
        assert g[R.top] == 1 and all(x > 0 for x in g)
        Eh = R.float_ops()                                               # orthonormal basis: Eh_ji = Eh_ij^T
        for i in range(1, R.N + 1):
            for j in range(i + 1, R.N + 1):
                d = abs(Eh[(i, j)] - Eh[(j, i)].T)
                assert (d.max() if d.nnz else 0.0) < 1e-12 * max(1.0, abs(Eh[(i, j)]).max()), (lam, i, j)
    _info('metric', '%d irreps (N = 3, 4, 5): E_ij^T G = G E_ji exact on every element of all off-diagonal pairs (0 violations), g[HW] = 1, g > 0, orthonormal Eh_ji = Eh_ij^T' % len(lams))


# ---------------------------------------------------------------------------------------- (v) Casimir
def test_casimir():
    lams = [(1, 0, 0, 0), (1, 1, 0, 0), (2, 1, 0, 0), (2, 1, 1, 0), (2, 2, 0, 0), (3, 2, 1, 0), (3, 3, 0, 0), (2, 1, 0, 0, 0), (2, 1, 0), (4, 2, 0), (3, 2, 1)]
    vals = {}
    for lam in lams:
        R = GLNRep.of(lam)
        N = R.N
        c2 = gr.casimir2(lam)
        for c in range(R.dim):
            acc = {}
            e = {c: Fr(1)}
            for i in range(1, N + 1):
                for j in range(1, N + 1):
                    for r, v in R.apply((i, j), R.apply((j, i), e)).items():
                        acc[r] = acc.get(r, 0) + v
            acc = {r: v for r, v in acc.items() if v != 0}
            assert acc == {c: Fr(c2)}, (lam, c, c2, acc)
        vals[str(list(lam))] = c2
    # the formula for N = 3 is W4's cas2: sum mu_i (mu_i + 4 - 2i)
    assert gr.casimir2((2, 1, 0)) == 2 * 4 + 1 + 0 == 9
    _info('casimir', 'C_2 = Sum_ij E_ij E_ji acts as the scalar Sum_i lam_i(lam_i + N + 1 - 2i) on every basis vector of %d irreps: %s' % (len(lams), vals))


# ---------------------------------------------------------------------------------------- (vi) Racah-Speiser
PRODUCTS = [((1, 0, 0, 0),) * 2, ((1, 1, 0, 0),) * 2, ((2, 1, 0, 0),) * 2, ((2, 1, 1, 0),) * 2, ((2, 2, 0, 0),) * 2, ((3, 2, 1, 0),) * 2, ((2, 2, 1, 0),) * 2,
            ((3, 3, 0, 0),) * 2,
            ((2, 1, 0, 0), (1, 1, 0, 0)), ((1, 1, 0, 0), (2, 1, 0, 0)), ((3, 1, 0, 0), (2, 2, 0, 0)), ((2, 2, 0, 0), (3, 1, 0, 0)), ((2, 2, 1, 0), (2, 1, 1, 0)),
            ((2, 1, 1, 0), (2, 2, 1, 0)), ((3, 2, 1, 0), (2, 1, 0, 0)), ((2, 1, 0, 0), (3, 2, 1, 0))]


def test_racah_speiser_dimsums():
    nprod = 0
    for lam, lamp in PRODUCTS + [((S, S, 0, 0),) * 2 for S in range(1, 11)]:
        rs = gr.racah_speiser(lam, lamp)
        assert all(m > 0 for m in rs.values())
        assert sum(m * gr.weyl_dim(mu) for mu, m in rs.items()) == gr.weyl_dim(lam) * gr.weyl_dim(lamp), (lam, lamp)
        assert all(sum(mu) == sum(lam) + sum(lamp) for mu in rs)
        # symmetric in the two factors
        assert rs == gr.racah_speiser(lamp, lam), (lam, lamp)
        # Casimir/Weyl cross-check on a few: every mu dominant
        assert all(all(mu[i] >= mu[i + 1] for i in range(3)) for mu in rs)
        nprod += 1
    # known decompositions
    assert gr.racah_speiser((1, 0, 0, 0), (1, 0, 0, 0)) == {(2, 0, 0, 0): 1, (1, 1, 0, 0): 1}
    assert gr.racah_speiser((1, 1, 0, 0), (1, 1, 0, 0)) == {(2, 2, 0, 0): 1, (2, 1, 1, 0): 1, (1, 1, 1, 1): 1}
    _info('racah_speiser_dimsums', 'sum_mu m dim V_mu = dim V_lam dim V_lam\' on %d products (G1.2 list in both orders + [S,S,0,0]^2, S = 1..10), RS symmetric in the two factors, all m > 0' % nprod)


def parse_munich(path=MUNICH):
    """{S: {(m1,m2,m3,m4): (dim, mult)}} from results/hw_timings.txt (SU(4) labels, mu_4 = 0)"""
    out = {}
    cur = None
    for line in open(path, encoding='utf8'):
        m = re.match(r'SU\(4\):\s+irreps mu in\s+\[(\d+),(\d+),0,0\] \(x\) \[(\d+),(\d+),0,0\]', line)
        if m:
            cur = out[int(m.group(1))] = {}
            continue
        m = re.match(r'\((\d+) (\d+) (\d+) (\d+)\)\s+(\d+)\s+(\d+)\s+([\d.]+) s', line)
        if m and cur is not None:
            cur[tuple(int(m.group(k)) for k in (1, 2, 3, 4))] = (int(m.group(5)), int(m.group(6)))
    return out


def test_munich_table():
    mun = parse_munich()
    assert sorted(mun) == list(range(1, 9)), sorted(mun)
    nirr = []
    for S in range(1, 9):
        lam = (S, S, 0, 0)
        rs = gr.racah_speiser(lam, lam)
        ours = {tuple(x - mu[3] for x in mu): (gr.weyl_dim(mu), m) for mu, m in rs.items()}
        assert ours == mun[S], (S, set(ours.items()) ^ set(mun[S].items()))
        assert len(ours) == (S + 1) * (S + 2) // 2
        assert all(m == 1 for (_, m) in ours.values())               # multiplicity free
        nirr.append(len(ours))
    _info('munich_table', 'Munich hw_timings.txt reproduced for S = 1..8: same irreps (our mu minus mu_4), same dimensions, all multiplicities 1; #irreps = %s = (S+1)(S+2)/2' % nirr)


def test_munich_gt_basis():
    """the GT basis of glN_reps is Munich's (Paul/Cpp/AlexCGC/ClebschGordan.cpp): Q(M) - 1 = our index (pattern order of clebsch::pattern::operator++ re-implemented),
    clebsch::pattern::lowering_coeff^2 = E_{l+1,l}[q',q]^2 g_q'/g_q exactly for every matrix element, no extra nonzero entries, and the diagonal sign gauge between the two
    orthonormal bases is trivial (all plain-basis simple raising/lowering coefficients are >= 0, as Munich's sqrt(...) >= 0)"""
    lams = [(1, 0, 0, 0), (1, 1, 0, 0), (2, 0, 0, 0), (2, 1, 0, 0), (2, 1, 1, 0), (3, 2, 1, 0), (3, 1, 1, 0), (2, 1, 0), (4, 2, 0)] + [(S, S, 0, 0) for S in range(1, 9)]
    nel = 0
    for lam in lams:
        r = gr.munich_check(lam)
        assert r['order_ok'], lam
        assert r['bad_sq'] == 0 and r['extra'] == 0 and r['gauge_bad'] == 0 and r['trivial_gauge'], (lam, r['bad_sq'], r['extra'], r['gauge_bad'])
        nel += r['nelem']
    _info('munich_gt_basis', '%d irreps (N = 3, 4): pattern order = Munich Q order, %d orthonormal lowering matrix elements equal Munich lowering_coeff^2 exactly, '
          'sign gauge trivial' % (len(lams), nel))


# ---------------------------------------------------------------------------------------- runner
def run_all(verbose=True):
    tests = [(k, v) for k, v in sorted(globals().items(), key=lambda kv: kv[1].__code__.co_firstlineno if callable(kv[1]) and hasattr(kv[1], '__code__') else 0)
             if k.startswith('test_') and callable(v)]
    nfail = 0
    for name, fn in tests:
        t0 = time.perf_counter()
        try:
            fn()
            status = 'PASS'
        except Exception:
            status = 'FAIL'
            nfail += 1
            traceback.print_exc()
        key = name[len('test_'):]
        if verbose:
            print('  %-4s %-26s (%6.2fs)  %s' % (status, name, time.perf_counter() - t0, INFO.get(key, '')), flush=True)
    print('  %d tests, %d failed' % (len(tests), nfail), flush=True)
    return nfail


if __name__ == '__main__':
    sys.exit(1 if run_all() else 0)
