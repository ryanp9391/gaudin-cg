"""cg_fast_gl4.py -- WP-G9: the float64 gl(4) Clebsch-Gordan engine (slice operators on the E12-cut space W1, w1 x w1 eigenproblem, physical subspace from the string powers, level recursion).

Brief: SPEC-G9-fast-cg-ubuntu.md.  Results: RESULTS-G9-fast-cg.md.  Everything here is float64 (no Fractions in the production path); the exact twin is G3's `copies`.

Pipeline for V_lam (x) V_lamp (site 1 = V_lam, site 2 = V_lamp), per irrep mu (Racah-Speiser multiplicity m >= 1):
  S0 setup      GLNRep of both factors (exact, G1), float orthonormal generators (adjacent: G1's float_ops; non-adjacent: FLOAT commutators), G1's left inverses (se.Recursion);
                the three universal reduced words (G3's ReducedOp tries) are evaluated at lamp on first use (counted in t_ops, like G3)
  S1 operators  V3, V4a, V4b on W1 = {m12 <= mu_1 - lamp_2} (coordinate cut of W = V_lam[mu - lamp] = ker E12^{k1}) in the ORTHONORMAL basis:
                  route 'B' = apply the reduced words to the w1 unit vectors (float, dense weight-block arrays);  route 'A' = float re-evaluation of the closed forms of G8d (V3, V4a only)
  S2 slices     (i) the physical subspace Q = null space of the stacked ROW-NORMALISED string-power matrices of E23^{k2}, E34^{k3} on W1 (SVD);
                (ii) the brief's eigenproblem: eig of G = V3 + 7/19 V4a + 11/13 V4b on W1, every eigenvector classified by its distance to Q and given the string statistic
                     s(u) = max_i ||A_i u|| / ||u|| (statistics + cross-check, count must be m);
                (iii) Rayleigh-Ritz inside Q: T_k = Q^T A_k Q, the combination of the triple with the best eigen-gap is diagonalised (accuracy step), joint diagonality is the acceptance test;
                canonical order (v3, v4a, v4b) ascending, plain slice u_c = u~_c / sqrt(g_c) normalised to first nonzero component 1 (G3's normalisation)
  S3 vectors    se.hw_from_slice (G1's float level recursion) -- output-bound
  S4 descendants  full CG tables from the HW vectors by the triangular sweep of the total lowering operators (descendants, cg_table)
Fallback: a sector whose physical count differs from m, whose null space has the wrong dimension, or whose Ritz step fails the joint-diagonality test is computed by G3's exact `copies`
(logged in sector['fallback']).  Lean mode (lean=True): multiplicity-one sectors skip S1 (no operator applied; the triple of such a sector is only a label).

Public API
  hw_vectors(lam, lamp, mus=None, orthonormal=False, route='B', refine=None, do_hw=True, strict=False, timing=None, keep_hw=True, lean=False)   dict mu -> Sector
  physical_slices(lam, lamp, mu, route='B', ...)      Sector without the HW vectors (slices, canonical triples, diagnostics)
  operators_W1(lam, lamp, mu, route='B')              (W, positions of W1 in W, {name: w1 x w1 matrix, orthonormal basis})
  descendants(lam, lamp, mu, sector, k), cg_table(lam, lamp)    S4
CLI: python -X utf8 cg_fast_gl4.py run l1 l2 l3 l4 [l1p l2p l3p l4p] [--mus a,b,c,d;...] [--out file] [--orthonormal] [--descendants]
     python -X utf8 cg_fast_gl4.py validate | t1 | t2 | t4s | t5 | t6 | bench [l1,l2,l3,l4 ...] | benchlean [...]      (validate = T1 + T4 + T5; outputs cg_fast_gl4_output.txt / _timings.txt)
"""
import os, sys, re, time, math, itertools, platform
from fractions import Fraction as Fr
from functools import cmp_to_key
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import scipy.sparse as sps
import glN_reps as gr
from glN_reps import GLNRep, dij, wadd
import slice_engine_gl4 as se
import slice_ops_gl4 as so
import gaudin_gl4_ops as go

TRIPLE = ('V3', 'V4a', 'V4b')
COMBO = np.array([float(c) for c in go.GENERIC_COMBO])      # (1, 7/19, 11/13)
TAU = 1e-6                                                  # string-test threshold on the ROW-NORMALISED statistic (brief: 1e-8 on the raw one); used for the statistics / cross-check only
ZTOL = 1e-9                                                 # "first nonzero component" threshold (relative to the max of the slice)
TIE = 1e-8                                                  # tolerance of the lexicographic triple order


# ======================================================================================== S0: setup (cached per product)
class Setup:
    """per (lam, lamp): GLNRep of both factors, the extended float generators, tries of the reduced words, weight-block cache, G1's recursion (lazy)"""
    _cache = {}

    @classmethod
    def of(cls, lam, lamp):
        key = (tuple(lam), tuple(lamp))
        if key not in cls._cache:
            cls._cache[key] = cls(*key)
        return cls._cache[key]

    @classmethod
    def clear(cls):
        cls._cache.clear()
        so.clear_caches()          # G1's caches, normal-ordering memo, universal reductions

    def __init__(self, lam, lamp):
        t0 = time.perf_counter()
        self.lam, self.lamp = lam, lamp
        self.R1 = GLNRep.of(lam)
        self.R2 = GLNRep.of(lamp) if lamp != lam else self.R1
        self.N = self.R1.N
        self.F = FloatGens(self.R1)                              # extended float generators of site 1 (the reduced words act on site 1)
        self._tries = None
        self.blocks = {}
        self._done = set()
        self._posw = None
        self.t_ops_setup = time.perf_counter() - t0
        self._rec = None
        self.t_rec = 0.0

    @property
    def tries(self):
        """the three reduced operators (universal reduction, polynomial in lamp, evaluated at lamp) -- built on first use, so that the reduction is counted in t_ops like G3 counts it in its t_ops"""
        if self._tries is None:
            self._tries = {n: so.reduced_op(n, self.lamp) for n in TRIPLE}
        return self._tries

    @property
    def rec(self):
        if self._rec is None:
            t = time.perf_counter()
            self._rec = se.Recursion.of(self.lam, self.lamp)
            self.t_rec = time.perf_counter() - t
        return self._rec

    def blk(self, ij, w):
        """dense block of the orthonormal-basis generator E_ij from the weight space w to w + d_ij: (target weight, array rows = target states, cols = source states) or None.
        All blocks of one generator are cut out of its COO form in one pass (first use)."""
        if ij not in self._done:
            self._build(ij)
        return self.blocks.get((ij, w))

    def _build(self, ij):
        R = self.R1
        if self._posw is None:
            wl = list(R.wsp)
            self._wl = wl
            wid = np.zeros(R.dim, dtype=np.int64)
            posw = np.zeros(R.dim, dtype=np.int64)
            for k, w_ in enumerate(wl):
                wid[R.wsp[w_]] = k
                posw[R.wsp[w_]] = np.arange(len(R.wsp[w_]))
            self._wid, self._posw = wid, posw
        coo = self.F[ij].tocoo()
        keep = coo.data != 0
        r, c, v = coo.row[keep], coo.col[keep], coo.data[keep]
        wc = self._wid[c]
        order = np.argsort(wc, kind='stable')
        r, c, v, wc = r[order], c[order], v[order], wc[order]
        cuts = np.nonzero(np.diff(wc))[0] + 1
        starts = np.concatenate([[0], cuts])
        ends = np.concatenate([cuts, [len(wc)]])
        for a_, b_ in zip(starts, ends):
            w_ = self._wl[wc[a_]]
            rr, cc, vv = r[a_:b_], c[a_:b_], v[a_:b_]
            tw = self._wl[self._wid[rr[0]]]
            B = np.zeros((len(R.wsp[tw]), len(R.wsp[w_])))
            B[self._posw[rr], self._posw[cc]] = vv
            self.blocks[(ij, w_)] = (tw, B)
        self._done.add(ij)


class FloatGens(dict):
    """float orthonormal-basis generators of a GLNRep: adjacent ones and the diagonal from R.float_ops(); non-adjacent by FLOAT commutators (Eh_ij = [Eh_i,j-1, Eh_j-1,j], Eh_ji = Eh_ij^T),
    so the Fractions of the exact non-adjacent generators (G1's lazy commutators) are never built"""

    def __init__(self, R):
        super().__init__()
        self.R = R
        self.base = R.float_ops()

    def __missing__(self, key):
        i, j = key
        if i == j or abs(i - j) == 1:
            val = self.base[key]
        elif j > i:
            A, B = self[(i, j - 1)], self[(j - 1, j)]
            val = (A @ B - B @ A).tocsr()
        else:
            val = self[(j, i)].T.tocsr()
        self[key] = val
        return val


# ======================================================================================== S1: operators on W1
def cut_space(S, mu):
    """W (global GT indices of V_lam[mu - lamp], ascending), the positions of W1 = {m12 <= mu_1 - lamp_2} inside W"""
    R = S.R1
    w0 = tuple(mu[k] - S.lamp[k] for k in range(S.N))
    W = R.wsp.get(w0)
    if W is None:
        return w0, np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
    cut = mu[0] - S.lamp[1]
    pos = np.array([a for a, c in enumerate(W) if R.pats[int(c)][1][0] <= cut], dtype=np.int64)
    return w0, W, pos


def _vadd(out, w, arr, c=1.0):
    cur = out.get(w)
    if cur is None:
        out[w] = c * arr if c != 1.0 else arr.copy()
    else:
        cur += c * arr


def _gen(S, g, vec):
    out = {}
    if g[0] == g[1]:
        k = g[0] - 1
        for w, arr in vec.items():
            if w[k] != 0:
                out[w] = arr * w[k]
        return out
    for w, arr in vec.items():
        b = S.blk(g, w)
        if b is not None:
            _vadd(out, b[0], b[1] @ arr)
    return out


def apply_trie(S, op, X, cache):
    """float action of a ReducedOp's trie (G3's ReducedOp.apply, same recursion, weight-block dense arrays) on X = {weight: array (states x columns)}; `cache` shares the leaf images"""
    def ev(node):
        c, ch = node
        out = {}
        if c != 0:
            fc = float(c)
            for w, a in X.items():
                _vadd(out, w, a, fc)
        for g, child in ch.items():
            cc, cch = child
            if not cch:
                e = cache.get(g)
                if e is None:
                    e = cache[g] = _gen(S, g, X)
                fc = float(cc)
                for w, a in e.items():
                    _vadd(out, w, a, fc)
            else:
                sub = ev(child)
                if sub:
                    for w, a in _gen(S, g, sub).items():
                        _vadd(out, w, a)
        return out
    return ev(op.root)


def operators_B(S, mu, check=True):
    """route B: V3, V4a, V4b on W1 (orthonormal basis, w1 x w1 float) by applying the reduced words to the w1 unit vectors.  Returns (W, pos, {name: matrix}, leak) with
    leak = max over the three operators of ||image outside W1||_F / max(||image||_F, 1) (must be ~1e-13: the operators preserve W1)"""
    w0, W, pos = cut_space(S, mu)
    w, w1 = len(W), len(pos)
    X0 = np.zeros((w, w1))
    X0[pos, np.arange(w1)] = 1.0
    X = {w0: X0}
    cache = {}
    mats = {}
    leak = 0.0
    for n in TRIPLE:
        out = apply_trie(S, S.tries[n], X, cache)
        img = out.get(w0)
        if img is None:
            img = np.zeros((w, w1))
        mats[n] = img[pos, :]
        if check:
            outside = np.ones(len(W), bool)
            outside[pos] = False
            o2 = sum(float(np.sum(a * a)) for w, a in out.items() if w != w0) + float(np.sum(img[outside] ** 2))
            leak = max(leak, math.sqrt(o2) / max(float(np.linalg.norm(mats[n])), 1.0))
    return W, pos, mats, leak


def operators_W1(lam, lamp, mu, route='B'):
    """(W, pos, {name: w1 x w1 matrix}) in the ORTHONORMAL basis; route 'B' = reduced words (V3, V4a, V4b), route 'A' = closed forms (V3, V4a)"""
    lam, lamp, mu = tuple(lam), tuple(lamp), tuple(mu)
    S = Setup.of(lam, lamp)
    if route == 'B':
        W, pos, mats, _ = operators_B(S, mu)
        return W, pos, mats
    if route == 'A':
        return operators_A(S, mu)
    raise ValueError(route)


# ======================================================================================== route A: float re-evaluation of the closed forms (G8d)
def _fl_R(kappa, lam):
    y = [kappa[i] - (i + 1) for i in range(3)]
    Y = [lam[s] - (s + 1) for s in range(4)]
    out_ = []
    for i in range(3):
        den = 1.0
        for j in range(3):
            if j != i:
                den *= (y[i] - y[j]) * (y[i] - y[j] + 1)
        if den == 0:
            out_.append(0.0)
            continue
        num = 1.0
        for s in range(4):
            num *= (y[i] - Y[s])
        out_.append(num / den)
    return out_, y


def _fl_a(kappa, lam):
    R, y = _fl_R(kappa, lam)
    a3 = 1.0 - sum(R)
    d = sum(lam) - sum(kappa)
    a4 = d + sum(R[i] * (2 * y[i] - sum(kappa) + 4) for i in range(3))
    return a3, a4


def _fl_gl3_entries(kappa, q, t, s, m1):
    """float twin of closed_forms_gl4.gl3_entries"""
    import closed_forms_gl4 as cf
    k1, k2, k3 = kappa
    q1, q2, q3 = q
    m2 = s - m1
    x = m1 - m2
    up = down = None
    if (x + 1) * (x + 2) != 0:
        up = -float((k1 - m1) * (k2 - m1 - 1) * (k3 - m1 - 2) * (t + 1 - m2)) * (t - m1 + q1 - q2) / ((x + 1) * (x + 2))
    if x * (x + 1) != 0:
        down = -float((k2 - m2) * (t - m1) * (k1 - m2 + 1) * (k3 - m2 - 1)) * (t + 1 - m2 + q1 - q2) / (x * (x + 1))
    xi = x * (x + 2)
    A0 = cf_A0(kappa, q, t, s)
    g = 3.0 * s / 16 - (k1 + k2 + k3) / 8 + (q1 + q2) / 8 - q3 / 4 + 3 / 8
    alpha = A0 + g * xi
    if xi != 0:
        Kp = float((s - 2 * k2) * (s - 2 * t) * (s - 2 * k1 - 2) * (s - 2 * k3 + 2)) * (2 - s + 2 * t + 2 * q1 - 2 * q2) / 16
        alpha += Kp / xi
    return alpha, up, down


_A0 = {}


def cf_A0(kappa, q, t, s):
    import closed_forms_gl4 as cf
    if not _A0:
        _A0['terms'] = [(e, float(co)) for e, co in cf.A0_TERMS]
    vals = [float(s), float(t), float(kappa[0]), float(kappa[1]), float(kappa[2]), float(q[0]), float(q[1]), float(q[2])]
    tot = 0.0
    for e, co in _A0['terms']:
        t_ = co
        for v, k in zip(vals, e):
            if k:
                t_ *= v ** k
        tot += t_
    return tot


_BC = {}


def _fl_b(name, kappa, lam, nu, lamp):
    """float b3 / b4a: coefficient tables of G8d evaluated in float (cached per (name, kappa, lam))"""
    import closed_forms_gl4 as cf
    key = (name, tuple(kappa), tuple(lam))
    c = _BC.get(key)
    if c is None:
        c = _BC[key] = [(m, float(v)) for m, v in cf.b_coefficients(name, kappa, lam).items()]
    v = [float(nu[0]), float(nu[1])] + [float(x) for x in lamp]
    tot = 0.0
    for m, co in c:
        if co == 0.0:
            continue
        t_ = co
        for x, k in zip(v, m):
            if k:
                t_ *= x ** k
        tot += t_
    return tot


def operators_A(S, mu, stats=None):
    """route A: V3 and V4a on W1 from labels alone, float, converted to the orthonormal basis (M~[r,c] = M[r,c] sqrt(g_r/g_c)).  The 12 factorising hopping classes and the diagonal /
    within-path entries are evaluated in float; the 6 Dm12 = 0 classes (cubic Q) use G4's exact evaluator (Fractions) entry by entry -- those entries are counted in stats['q_exact']"""
    import closed_forms_gl4 as cf
    import structure_gl4 as sg
    lam, lamp = S.lam, S.lamp
    R = S.R1
    w0, W, pos = cut_space(S, mu)
    nu = w0
    w, w1 = len(W), len(pos)
    posW1 = {int(W[a]): j for j, a in enumerate(pos)}
    pats = [R.pats[int(c)] for c in W]
    idx = {int(c): a for a, c in enumerate(W)}
    q3 = tuple(lamp[:3])
    hops = sg.f_classes()
    ops = ('V3', 'V4a')
    M = {n: np.zeros((w1, w1)) for n in ops}
    ab = {}
    qrec = None
    for jc, a in enumerate(pos):
        p = pats[a]
        e = cf.entries6(p)
        kappa = p[2]
        t = e[0]
        s = e[1] + e[2]
        if kappa not in ab:
            ab[kappa] = _fl_a(kappa, lam)
        a3, a4 = ab[kappa]
        al, up, dn = _fl_gl3_entries(kappa, q3, t, s, e[1])
        asc = {'V3': a3, 'V4a': a4}
        for n in ops:
            M[n][jc, jc] = asc[n] * al + _fl_b(n, kappa, lam, nu, lamp)
        for dm, val in ((1, up), (-1, dn)):
            if val is None:
                continue
            tg = (e[0], e[1] + dm, e[2] - dm, e[3], e[4], e[5])
            r = R.index.get(cf.pat_of(tg, lam))
            if r is not None and r in posW1:
                for n in ops:
                    M[n][posW1[r], jc] = asc[n] * val
        # hopping classes (12): universal factorised formula, float
        pt = [float(x) for x in cf.pt_of(p, lam, lamp)]
        for d in hops:
            tg = tuple(x + y_ for x, y_ in zip(e, d))
            r = R.index.get(cf.pat_of(tg, lam))
            if r is None or r not in posW1:
                continue
            m3 = _universal_float(pt, d)
            M['V3'][posW1[r], jc] = m3
            M['V4a'][posW1[r], jc] = (pt[{1: 3, 2: 4, 3: 5}[_spectator(d)]] + 3) * m3
        for d in cf.QCLASSES:
            tg = tuple(x + y_ for x, y_ in zip(e, d))
            r = R.index.get(cf.pat_of(tg, lam))
            if r is None or r not in posW1:
                continue
            if qrec is None:
                qrec = {(n, dd): cf.q_record(n, dd) for n in ops for dd in cf.QCLASSES}
            ptx = cf.pt_of(p, lam, lamp)
            for n in ops:
                M[n][posW1[r], jc] = float(cf.eval_qrec(qrec[(n, d)], ptx, cf._RNG))
            if stats is not None:
                stats['q_exact'] = stats.get('q_exact', 0) + 1
    # plain -> orthonormal basis: M~[r,c] = M[r,c] sqrt(g_r / g_c)
    sg1 = S.R1.sqrt_metric()[W[pos]]
    D = sg1[:, None] / sg1[None, :]
    return W, pos, {n: M[n] * D for n in ops}


def _spectator(d):
    d3 = {1: d[3], 2: d[4], 3: d[5]}
    return next(k for k, v in d3.items() if v == 0)


def _universal_float(pt, d):
    row3 = {1: 3, 2: 4, 3: 5}
    a, b = (1, 2) if d[1] == 1 else (2, 1)
    d3 = {1: d[3], 2: d[4], 3: d[5]}
    i = next(k for k, v in d3.items() if v == 1)
    j = next(k for k, v in d3.items() if v == -1)
    k = next(k for k, v in d3.items() if v == 0)
    ya, yb, yi, yj, yk = pt[a], pt[b], pt[row3[i]], pt[row3[j]], pt[row3[k]]
    num = (pt[0] - yb) * (pt[0] - ya + pt[10] - pt[11]) * (ya - yj) * (ya - yj + 1) * (ya - yk)
    for s_ in range(4):
        num *= (yi - pt[6 + s_])
    den = (ya - yb) * (ya - yb + 1) * (yi - yj) * (yi - yj + 1) * (yj - yk) * (yi - yk)
    return num / den


# ======================================================================================== S2: eigenproblem, string test, canonical copies
def strings_of(S, mu):
    """[(i, k_i)] for the non-void strings i = 2, ..., N-1 (the string i = 1 is the coordinate cut)"""
    w0 = tuple(mu[k] - S.lamp[k] for k in range(S.N))
    out_ = []
    for i in range(2, S.N):
        k = S.lamp[i - 1] - S.lamp[i] + 1
        if wadd(w0, dij(S.N, i, i + 1), k) in S.R1.wsp:
            out_.append((i, k))
    return out_


def string_images(S, mu, W, pos):
    """[P_i]: the matrices (rows = states of the image weight space, cols = the w1 basis vectors of W1) of E_{i,i+1}^{k_i} for the non-void strings i = 2, ..., N-1, orthonormal basis
    (the images leave W1 and W: nothing is restricted).  A column vector u of W1 is physical iff P_i u = 0 for all i;  s(u) = max_i ||P_i u|| / ||u||."""
    w0 = tuple(mu[k] - S.lamp[k] for k in range(S.N))
    w1 = len(pos)
    out_ = []
    for i, k in strings_of(S, mu):
        full = np.zeros((len(W), w1))
        full[pos, np.arange(w1)] = 1.0
        vec = {w0: full}
        for _ in range(k):
            vec = _gen(S, (i, i + 1), vec)
            if not vec:
                break
        if vec:
            out_.append(np.vstack([vec[w] for w in sorted(vec)]))
    return out_


def _row_normalised(P):
    """stack the string-power matrices and scale every row to unit 2-norm (the null space is unchanged; the dynamic range of E^k's entries, up to 1e7 at k ~ 9, otherwise puts the first
    non-null singular value at 1e-7 and the subspace error at 1e-16 / 1e-7)"""
    A = np.vstack(P)
    nr = np.linalg.norm(A, axis=1)
    if nr.size == 0 or nr.max() == 0.0:
        return np.zeros((0, A.shape[1]))
    keep = nr > 1e-14 * nr.max()
    return A[keep] / nr[keep, None]


def _lex_key(tol=TIE):
    def cmp(a, b):
        for x, y in zip(a, b):
            if abs(x - y) > tol * max(1.0, abs(x), abs(y)):
                return -1 if x < y else 1
        return 0
    return cmp_to_key(cmp)


class Sector(dict):
    """result for one irrep: keys mu, m, w, w1, trip (m x 3, canonical order), slices (w x m plain, first nonzero component 1), uhat (w x m orthonormal, unit 2-norm), W, hw (se.HW or None),
    s_phys (max s over physical), s_spur (min s over spurious), fallback (None or reason), t_ops, t_eig, t_hw, ncomplex, gapG (min distance of the physical G-eigenvalues)"""
    __getattr__ = dict.__getitem__


def _finish_sector(S, mu, W, pos, trip, Uhat, m, info):
    """canonical order, normalisation; Uhat (w1 x m) orthonormal-basis eigenvectors (real), trip (m x 3).  Returns the Sector fields"""
    order = sorted(range(m), key=lambda k: _lex_key()(tuple(trip[k])))
    trip = trip[order]
    Uhat = Uhat[:, order]
    w = len(W)
    sg1 = S.R1.sqrt_metric()[W]
    full_hat = np.zeros((w, m))
    full_hat[pos, :] = Uhat
    plain = full_hat / sg1[:, None]
    struct0 = np.abs(plain).max(axis=1) <= ZTOL * np.abs(plain).max()
    sl = np.zeros((w, m))
    smallest = np.inf
    for k in range(m):
        col = plain[:, k]
        mx = np.abs(col).max()
        cand = np.nonzero((~struct0) & (np.abs(col) > ZTOL * mx))[0]
        first = cand[0]
        sl[:, k] = col / col[first]
        full_hat[:, k] = full_hat[:, k] / col[first]
        nz = np.abs(sl[:, k])
        nz = nz[nz > ZTOL * nz.max()]
        smallest = min(smallest, nz.min() / np.abs(sl[:, k]).max())
    return trip, sl, full_hat, smallest


def _slices_float(S, mu, W, pos, mats, m, strict=False, adapt=True):
    """S2 on the operator matrices of W1.
    (i) the physical SUBSPACE Q (w1 x m) = the null space of the stacked, row-normalised string-power matrices A_i (SVD; the count of singular values below 1e-10 sigma_max must be m);
    (ii) the brief's eigenproblem: G = V3 + 7/19 V4a + 11/13 V4b, numpy eig, every eigenvector u classified by its distance to Q (physical <=> u in Q) and given the string statistic s(u) = max_i ||A_i u|| / ||u||
        (statistics s_phys = max over physical, s_spur = min over spurious; count of physical eigenvectors must be m);
    (iii) Rayleigh-Ritz inside Q: T_k = Q^T A_k Q; the combination of the triple with the best eigen-gap is diagonalised, joint diagonality is the acceptance test.
    Returns (ok, info); ok False -> the caller falls back to the exact path."""
    w1 = len(pos)
    G = COMBO[0] * mats['V3'] + COMBO[1] * mats['V4a'] + COMBO[2] * mats['V4b']
    if w1 == 1:
        E = np.array([G[0, 0]])
        U = np.ones((1, 1))
    else:
        E, U = np.linalg.eig(G)
    ncomplex = int(np.sum(np.abs(np.imag(E)) > 1e-9 * max(1.0, np.abs(E).max())))
    U = U / np.linalg.norm(U, axis=0)
    trivial = (w1 == m and not strict)
    P = [] if trivial else [_row_normalised([Pi]) for Pi in string_images(S, mu, W, pos)]
    P = [Pi for Pi in P if Pi.shape[0] > 0]            # strings whose image is identically zero on W1 impose nothing
    info = dict(w1=w1, ncomplex=ncomplex)
    if not P and w1 != m:
        info.update(nphys=w1, ritz='no string condition left but w1 = %d != m = %d' % (w1, m))
        return False, info
    if P:
        _, sv, Vt = np.linalg.svd(np.vstack(P), full_matrices=True)
        sv = np.concatenate([sv, np.zeros(w1 - len(sv))])
        nnull = int(np.sum(sv < 1e-10 * sv[0]))
        info['sv_null'] = float(sv[w1 - m:].max()) / float(sv[0])
        info['sv_gap'] = float(sv[w1 - m - 1]) / float(sv[0]) if w1 > m else np.inf
        info['nnull'] = nnull
        if nnull != m:
            info['ritz'] = 'null space of dimension %d != m = %d' % (nnull, m)
            info['nphys'] = nnull
            return False, info
        Q = Vt[w1 - m:].T
        s = np.max([np.linalg.norm(Pi @ U, axis=0) for Pi in P], axis=0) / np.linalg.norm(U, axis=0)
        dist = np.linalg.norm(U - Q @ (Q.T @ U), axis=0) / np.linalg.norm(U, axis=0)      # distance of the eigenvector to the physical subspace (physical: ~eps cond, spurious: O(1))
        phys = np.nonzero(dist < 1e-3)[0]
        info['nphys_tau'] = int(np.sum(s < TAU))          # the brief's test (s < tau) on the row-normalised statistic: how often it alone would give m
        info['mixed'] = int(np.sum((dist > 1e-8) & (dist < 0.1)))
    else:
        Q = np.eye(w1)
        s = np.zeros(w1)
        phys = np.arange(w1)
        info['mixed'] = 0
        info['nphys_tau'] = w1
    spur = np.setdiff1d(np.arange(w1), phys)
    info['s_phys'] = float(s[phys].max()) if len(phys) else 0.0
    info['s_spur'] = float(s[spur].min()) if len(spur) else np.inf
    info['nphys'] = len(phys)
    if len(phys) != m:
        info['ritz'] = 'eigenvector classification found %d physical, m = %d' % (len(phys), m)
        return False, info
    Up = U[:, phys]
    if np.iscomplexobj(Up):
        big = np.abs(Up).argmax(axis=0)
        ph = Up[big, np.arange(m)]
        Up = Up / (ph / np.abs(ph))
        if np.abs(Up.imag).max() > 1e-8 * np.abs(Up).max():
            info['complex_phys'] = float(np.abs(Up.imag).max())
            return False, info
        Up = Up.real
    Up = Up / np.linalg.norm(Up, axis=0)
    trip = np.zeros((m, 3))
    for t, n in enumerate(TRIPLE):
        AU = mats[n] @ Up
        trip[:, t] = np.sum(Up * AU, axis=0) / np.sum(Up * Up, axis=0)
    ev = np.real(E[phys])
    info['gapG'] = float(min([abs(a - b) for a, b in itertools.combinations(ev, 2)] + [np.inf]))
    info['trip_raw'] = trip
    if not adapt:
        return True, dict(info, Up=Up, trip=trip)
    info['subspace_dev'] = float(np.linalg.norm(Q - Up @ np.linalg.lstsq(Up, Q, rcond=None)[0]) / math.sqrt(m))
    ok2, Up2, trip2, inf2 = _ritz(mats, Q, trip, m)
    info.update(inf2)
    if not ok2:
        return False, info
    return True, dict(info, Up=Up2, trip=trip2)


COMBOS = [tuple(COMBO)] + [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (1.0, -1.0, 0.0), (1.0, 0.0, -1.0), (0.0, 1.0, -1.0), (1.0, 1.0, 1.0), (1.0, -0.5, 0.3), (0.3, 1.0, -0.7), (-0.6, 0.4, 1.0)]
_ROT = {}


def _rot(m):
    if m not in _ROT:
        _ROT[m], _ = np.linalg.qr(np.random.default_rng(12345).standard_normal((m, m)))
    return _ROT[m]


def _ritz(mats, Q, trip, m):
    """Rayleigh-Ritz on the physical subspace Q (w1 x m, orthonormal): T_k = Q^T A_k Q; the combination c of the triple with the best  min pairwise gap / ||T_c||  (gaps from the
    first-pass triples) is diagonalised (eigenvector conditioning ~ ||T_c|| / gap_c; the spurious spectrum plays no role); joint diagonality of the three T_k in that eigenbasis is the
    acceptance test.  A fixed generic rotation of Q avoids triangular T (LAPACK's balancing then loses 1e-8)."""
    Q = Q @ _rot(m)
    T = [Q.T @ mats[n] @ Q for n in TRIPLE]
    invres = max(float(np.linalg.norm(mats[n] @ Q - Q @ T[k]) / max(np.linalg.norm(mats[n] @ Q), 1e-300)) for k, n in enumerate(TRIPLE))
    if m == 1:
        return True, Q, np.array([[T[0][0, 0], T[1][0, 0], T[2][0, 0]]]), dict(ritz_invres=invres)
    floor = 1e-6 * max(np.linalg.norm(Tk, 2) for Tk in T)        # a combination whose T_c is numerical noise (V4a = 0 on some sectors) is never chosen by its noise 'gap'
    best = None
    for c in COMBOS:
        d = trip @ np.array(c)
        gap = min(abs(a - b) for a, b in itertools.combinations(d, 2))
        Tc = sum(ck * Tk for ck, Tk in zip(c, T))
        sc = gap / max(np.linalg.norm(Tc, 2), floor)
        if best is None or sc > best[0]:
            best = (sc, c, Tc)
    sc, c, Tc = best
    E, V = np.linalg.eig(Tc)
    if np.iscomplexobj(V):
        if np.abs(E.imag).max() > 1e-9 * max(1.0, np.abs(E).max()):
            return False, None, None, dict(ritz='complex Ritz eigenvalues', ritz_combo=c)
        V = V.real
    V = V / np.linalg.norm(V, axis=0)
    try:
        Vi = np.linalg.inv(V)
    except np.linalg.LinAlgError:
        return False, None, None, dict(ritz='singular Ritz eigenvectors', ritz_combo=c)
    trip2 = np.zeros((m, 3))
    off = 0.0
    for k in range(3):
        D = Vi @ T[k] @ V
        trip2[:, k] = np.diag(D)
        off = max(off, float(np.abs(D - np.diag(np.diag(D))).max()) / max(1.0, float(np.abs(T[k]).max())))
    info = dict(ritz_combo=c, ritz_score=sc, ritz_invres=invres, ritz_off=off, ritz_cond=float(np.linalg.cond(V)))
    if off > 1e-7:
        info['ritz'] = 'joint diagonality %.1e' % off
        return False, None, None, info
    return True, Q @ V, trip2, info


def _fallback_sector(lam, lamp, mu, dps=50):
    """G3's exact copies, converted to float (the logged exact path)"""
    cps = so.copies(lam, lamp, mu, dps=dps)
    prv = cps.prv
    F = so.copy_float_slices(prv, cps)
    trip = np.array([[float(x) for x in cp['ev']] for cp in cps])
    return prv, F, trip, cps.info


def physical_slices(lam, lamp, mu, route='B', refine=None, strict=False, adapt=True, lean=False):
    """S1 + S2 for one irrep: a Sector without HW vectors.  route 'B' or 'A' (A: V4b is taken from route B because it has no closed form; the G of route A is V3 + 7/19 V4a + 11/13 V4b)"""
    if refine:
        raise NotImplementedError('S2-prime (high-precision refinement of the slices) is not implemented in G9; see RESULTS-G9 section 7')
    lam, lamp, mu = tuple(lam), tuple(lamp), tuple(mu)
    S = Setup.of(lam, lamp)
    m = gr.racah_speiser(lam, lamp)[mu] if not hasattr(S, '_rs') else S._rs[mu]
    if lean and m == 1:
        return _lean_sector(S, mu)
    t0 = time.perf_counter()
    W, pos, mats, leak = operators_B(S, mu)
    if route == 'A':
        _, _, matsA = operators_A(S, mu)
        mats = dict(mats, **matsA)
    t1 = time.perf_counter()
    ok, info = _slices_float(S, mu, W, pos, mats, m, strict=strict, adapt=adapt)
    sec = Sector(mu=mu, m=m, w=len(W), w1=len(pos), W=W, pos=pos, leak=leak, mats=mats, hw=None, fallback=None, t_ops=t1 - t0, info=info)
    if ok:
        trip, sl, full_hat, smallest = _finish_sector(S, mu, W, pos, info['trip'], info['Up'], m, info)
        sec.update(trip=trip, slices=sl, uhat=full_hat, smallest=smallest, s_phys=info['s_phys'], s_spur=info['s_spur'], ncomplex=info['ncomplex'],
                   gapG=info.get('gapG', np.inf))
    else:
        prv, F, trip, jinfo = _fallback_sector(lam, lamp, mu)
        sg1 = S.R1.sqrt_metric()[W]
        sec.update(trip=trip, slices=F, uhat=F * sg1[:, None], smallest=np.nan, s_phys=info.get('s_phys', np.nan), s_spur=info.get('s_spur', np.nan), ncomplex=info['ncomplex'],
                   gapG=np.inf, fallback=info.get('ritz') or ('nphys=%d != m=%d%s' % (info['nphys'], m, ' (complex physical eigenvector)' if 'complex_phys' in info else '')))
    sec['t_eig'] = time.perf_counter() - t1
    return sec


def _lean_sector(S, mu):
    """multiplicity one: M_mu is a line, the slice is the null vector of the stacked string powers on W1 (or the unit vector if w1 = 1); NO operator is applied (the triple of a
    multiplicity-one irrep is a label, not needed for the HW vector; trip = nan)"""
    t0 = time.perf_counter()
    w0, W, pos = cut_space(S, mu)
    w1 = len(pos)
    s_phys, s_spur = 0.0, np.inf
    if w1 == 1:
        u = np.ones((1, 1))
    else:
        P = string_images(S, mu, W, pos)
        _, sv, Vt = np.linalg.svd(_row_normalised(P), full_matrices=True)
        sv = np.concatenate([sv, np.zeros(w1 - len(sv))])
        u = Vt[-1:].T
        s_phys = float(sv[-1] / sv[0])
        s_spur = float(sv[-2] / sv[0]) if w1 > 1 else np.inf
    trip = np.full((1, 3), np.nan)
    trip, sl, full_hat, smallest = _finish_sector(S, mu, W, pos, trip, u, 1, {})
    return Sector(mu=mu, m=1, w=len(W), w1=w1, W=W, pos=pos, leak=0.0, mats=None, hw=None, fallback=None, t_ops=0.0, t_eig=time.perf_counter() - t0, info=dict(w1=w1, nphys=1, ncomplex=0, lean=True),
                  trip=trip, slices=sl, uhat=full_hat, smallest=smallest, s_phys=s_phys, s_spur=s_spur, ncomplex=0, gapG=np.inf)


# ======================================================================================== S3: the whole product
def hw_vectors(lam, lamp, mus=None, orthonormal=False, route='B', refine=None, do_hw=True, strict=False, timing=None, keep_hw=True, lean=False):
    """all irreps (or the list `mus`) of V_lam (x) V_lamp: dict mu -> Sector, with sec.hw = se.HW (COO rows, cols in the plain GT basis, vals[:, k] = the k-th copy in G3's canonical
    normalisation; orthonormal=True: vals = orthonormal-basis components, each copy G(x)G-normalised to 1) and the diagnostics.  timing (dict) receives t_setup, t_ops, t_eig, t_hw"""
    lam, lamp = tuple(lam), tuple(lamp)
    t0 = time.perf_counter()
    S = Setup.of(lam, lamp)
    rs = gr.racah_speiser(lam, lamp)
    S._rs = rs
    if do_hw:
        S.rec          # G1's left inverses: part of the setup
    t_setup = time.perf_counter() - t0
    mus = sorted(rs, reverse=True) if mus is None else [tuple(x) for x in mus]
    for mu in mus:
        if mu not in rs:
            raise ValueError('mu = %s is not an irrep of %s x %s (our gl(4) labels, |mu| = |lam| + |lamp|)' % (mu, lam, lamp))
    out = {}
    t_ops = t_eig = t_hw = 0.0
    for mu in mus:
        sec = physical_slices(lam, lamp, mu, route=route, refine=refine, strict=strict, lean=lean)
        t_ops += sec['t_ops']
        t_eig += sec['t_eig']
        if do_hw:
            t = time.perf_counter()
            hw = se.hw_from_slice(lam, lamp, mu, sec['slices'])
            if orthonormal:
                hw.vals_hat = hw.vals_hat / np.linalg.norm(hw.vals_hat, axis=0)
            sec['nent'] = hw.nW * sec['m']
            sec['hw'] = hw if keep_hw else None             # keep_hw=False: stream (the HW vectors of [15,10,5,0]^2 are 2.5e9 coefficients)
            sec['t_hw'] = time.perf_counter() - t
            t_hw += sec['t_hw']
        out[mu] = sec
    if timing is not None:
        timing.update(t_setup=t_setup + S.t_ops_setup * 0, t_ops=t_ops, t_eig=t_eig, t_slices=t_ops + t_eig, t_hw=t_hw)
    return out


# ======================================================================================== output
def write_vectors(path, lam, lamp, secs, orthonormal=False):
    """text: for every irrep and copy the canonical triple, the slice (plain) and, when present, the HW coefficients in Munich's rows (se.to_munich)"""
    with open(path, 'w', encoding='utf-8') as f:
        f.write('# cg_fast_gl4 HW vectors of %s x %s (float64)\n' % (lam, lamp))
        for mu, sec in secs.items():
            for k in range(sec['m']):
                f.write('mu %s copy %d/%d triple %s\n' % (mu, k + 1, sec['m'], ' '.join('%.12g' % x for x in sec['trip'][k])))
                f.write('slice ' + ' '.join('%d:%.15g' % (int(c), x) for c, x in zip(sec['W'], sec['slices'][:, k]) if x != 0) + '\n')
                if sec['hw'] is not None:
                    M = se.to_munich(lam, lamp, sec['hw'], k)
                    for (a, b), v in sorted(M.items()):
                        f.write('%d\t%d\t%d\t%.15f\n' % (a, b, gr.weyl_dim(mu), v))


# ======================================================================================== validation helpers (T1 twin)
def exact_twin(lam, lamp, mu, dps=50, cache=None):
    """G3's exact path: (prv, copies, float slices (w x m), float triples (m x 3))"""
    key = (lam, lamp, mu)
    if cache is not None and key in cache:
        return cache[key]
    cps = so.copies(lam, lamp, mu, dps=dps)
    F = so.copy_float_slices(cps.prv, cps)
    trip = np.array([[float(x) for x in cp['ev']] for cp in cps])
    res = (cps, F, trip)
    if cache is not None:
        cache[key] = res
    return res


def compare_sector_exact(lam, lamp, sec, ref):
    """T1: max relative deviation of every slice component (relative to the max of the copy), of the triples (relative to max(1,|v|)), for the copies in canonical order"""
    cps, F, trip = ref
    dev = np.abs(sec['slices'] - F).max(axis=0) / np.abs(F).max(axis=0)
    tdev = np.abs(sec['trip'] - trip) / np.maximum(1.0, np.abs(trip))
    return float(dev.max()), float(tdev.max())


def hw_checks(lam, lamp, sec):
    """HW residual (se.hw_residual, relative) and the G(x)G-orthogonality of the copies of the sector (orthonormal-basis Gram matrix of the normalised copies, off-diagonal max)"""
    hw = sec['hw']
    res, wd = se.hw_residual(hw)
    V = hw.vals_hat
    V = V / np.linalg.norm(V, axis=0)
    Gm = V.T @ V
    off = float(np.abs(Gm - np.eye(Gm.shape[0])).max()) if Gm.shape[0] > 1 else 0.0
    return res, wd, off


def mats_vs_exact(lam, lamp, mu, route='B'):
    """T1 twin: float operators on W1 against G3's exact weight_space_mats conjugated by sqrt(g): max |M~ - M~_ref| / max |M~_ref| per operator"""
    S = Setup.of(tuple(lam), tuple(lamp))
    W, pos, mats = operators_W1(lam, lamp, mu, route)
    R = S.R1
    names = list(mats)
    ex = so.weight_space_mats(R, {n: so.reduced_op(n, tuple(lamp)) for n in names}, W)
    sg1 = R.sqrt_metric()[W]
    out_ = {}
    for n in names:
        ref = np.array([[float(x) for x in row] for row in ex[n]]).reshape(len(W), len(W)) * (sg1[:, None] / sg1[None, :])
        sub = ref[np.ix_(pos, pos)]
        # exactness of the invariance of W1: the columns of W1 have no component outside W1
        outside = np.setdiff1d(np.arange(len(W)), pos)
        leak = float(np.abs(ref[np.ix_(outside, pos)]).max()) if len(outside) and len(pos) else 0.0
        scale = max(np.abs(sub).max(), 1.0)       # relative to max |M~| (to 1 if the operator is identically zero on the sector: V4a = 0 happens)
        out_[n] = (float(np.abs(mats[n] - sub).max() / scale), leak / scale)
    return out_


def check_cut(lam, lamp, mu):
    """W1 is exactly the kernel of E12^{k1} on W (float): returns (ok, max kernel residual, min image norm outside W1)"""
    S = Setup.of(tuple(lam), tuple(lamp))
    w0, W, pos = cut_space(S, mu)
    k1 = S.lamp[0] - S.lamp[1] + 1
    X = np.eye(len(W))
    vec = {w0: X}
    for _ in range(k1):
        vec = _gen(S, (1, 2), vec)
        if not vec:
            break
    if not vec:
        return True, 0.0, np.inf
    img = sum(np.sum(a ** 2, axis=0) for a in vec.values())
    img = np.sqrt(img)
    inside = np.zeros(len(W), bool)
    inside[pos] = True
    mx_in = float(img[inside].max()) if inside.any() else 0.0
    mn_out = float(img[~inside].min()) if (~inside).any() else np.inf
    return (mx_in < 1e-12 and mn_out > 1e-6), mx_in, mn_out


# ======================================================================================== S4: descendants
def descendants(lam, lamp, mu, sec, k=0):
    """S4: the full CG table of copy k of the irrep mu.  Psi_q (q a GT state of V_mu) = the product-space vector (orthonormal basis of both factors, unit norm) with
    E_{l+1,l}^{tot} Psi_q = sum_t L_l[t, q] Psi_t,  L_l = the orthonormal lowering matrix of V_mu (GLNRep float_ops), E^{tot} = E(1) (x) 1 + 1 (x) E(2); starting from the HW vector Psi_top the vectors
    are obtained by the triangular sweep: for a known Psi_q and level l with exactly ONE unknown target t,  Psi_t = (E^{tot} Psi_q - sum_{t' known} L_l[t', q] Psi_t') / L_l[t, q].
    Returns (Rm, {q: sparse (n1*n2 x 1) column of orthonormal components, flat index nu * dim2 + q2}); the sweep is output-bound (a few sparse matvecs per state of V_mu)."""
    S = Setup.of(tuple(lam), tuple(lamp))
    R1, R2 = S.R1, S.R2
    Rm = GLNRep.of(tuple(mu))
    N = R1.N
    Eh1, Eh2 = R1.float_ops(), R2.float_ops()
    n1, n2 = R1.dim, R2.dim
    hw = sec['hw']
    v = hw.vals_hat[:, k]
    top = sps.csr_matrix((v / np.linalg.norm(v), (hw.rows, hw.cols)), shape=(n1, n2))
    Lm = {l: Rm.float_ops()[(l + 1, l)].tocsc() for l in range(1, N)}
    tg = {}
    for l, M in Lm.items():
        for q in range(Rm.dim):
            col = M[:, q]
            tg[(l, q)] = (col.indices.copy(), col.data.copy())
    known = {Rm.top: top}
    stack = [Rm.top]
    pend = {}
    while stack:
        q = stack.pop()
        for l in range(1, N):
            idx, dat = tg[(l, q)]
            unk = [i for i in range(len(idx)) if int(idx[i]) not in known]
            if len(unk) == 0:
                continue
            if len(unk) > 1:
                pend.setdefault(q, set()).add(l)
                continue
            ij = (l + 1, l)
            Xq = known[q]
            Y = (Eh1[ij] @ Xq + Xq @ Eh2[ij].T).tocsr()
            for i in range(len(idx)):
                t = int(idx[i])
                if t in known:
                    Y = Y - dat[i] * known[t]
            tnew = int(idx[unk[0]])
            known[tnew] = (Y / dat[unk[0]]).tocsr()
            stack.append(tnew)
            # a state that was waiting for this target may be solvable now
            for q2 in list(pend):
                stack.append(q2)
            pend.clear()
    if len(known) != Rm.dim:
        raise ArithmeticError('descendant sweep reached %d of %d states of V_mu' % (len(known), Rm.dim))
    return Rm, {q: X.reshape((n1 * n2, 1)).tocsc() for q, X in known.items()}


def cg_table(lam, lamp, secs=None):
    """the whole CG matrix of the product (orthonormal bases): rows (irrep mu, copy k, GT state M of V_mu), columns the product basis nu * dim2 + q; sparse.  Unitarity: C C^T = 1.
    Returns the matrix (rows ordered by sorted irreps descending, copies in canonical order, states in the GT order of V_mu)"""
    lam, lamp = tuple(lam), tuple(lamp)
    secs = secs or hw_vectors(lam, lamp)
    blocks = []
    for mu in sorted(secs, reverse=True):
        sec = secs[mu]
        for k in range(sec['m']):
            Rm, cols = descendants(lam, lamp, mu, sec, k)
            blocks.append(sps.hstack([cols[q] for q in range(Rm.dim)]).T.tocsr())
    return sps.vstack(blocks).tocsr()


def munich_full_table(lam, lamp, mu, alpha=1, timeout=600):
    """all rows (Q(M), Q(M'), Q(M''), CGC) of the irrep mu (OUR label) of V_lam (x) V_lamp from `clebsch` menu 5 (not only the HW rows, unlike hwcgc.sh); alpha = outer multiplicity index.
    Returns {(Q(M''), Q(M'), Q(M)) ...} as {Q(M''): {(Q(M), Q(M')): value}}"""
    import subprocess
    import munich_bridge as mb
    mm = mb.munich_label(mu)
    inp = '5\n4\n%s\n%s\n%s\n%d\n\n0\n' % (mb._s(lam), mb._s(lamp), mb._s(mm), alpha)
    r = subprocess.run(['./clebsch'], cwd=mb.CGC_DIR, input=inp, capture_output=True, text=True, timeout=timeout, env=dict(os.environ, LC_ALL='C'))
    out = {}
    for line in r.stdout.splitlines():
        f = line.replace(',', ' ').split()
        if len(f) == 4 and all(re.fullmatch(r'\d+', x) for x in f[:3]):
            try:
                out.setdefault(int(f[2]), {})[(int(f[0]), int(f[1]))] = float(f[3])
            except ValueError:
                pass
    return out


# ======================================================================================== stages: T1 .. T5 (the acceptance runs of SPEC-G9)
OUTPUT = os.path.join(HERE, 'cg_fast_gl4_output.txt')
TIMINGS = os.path.join(HERE, 'cg_fast_gl4_timings.txt')
RAWDIR = os.path.normpath(os.path.join(HERE, '..', '..', 'Cpp', 'AlexCGC', 'results', 'su4-same-machine', 'g9'))
_LOG = {'out': None}


def say(s='', to='out'):
    print(s, flush=True)
    path = OUTPUT if to == 'out' else TIMINGS   # to = 'out' | 't'
    with open(path, 'a', encoding='utf-8') as f:
        f.write(s + '\n')


def hdr(s, to='out'):
    say('', to)
    say('=' * 140, to)
    say(s, to)
    say('=' * 140, to)


def lam_str(lam):
    return '[' + ','.join(str(x) for x in lam) + ']'


def prod_str(lam, lamp):
    return lam_str(lam) + ('^2' if tuple(lam) == tuple(lamp) else 'x' + lam_str(lamp))


P_SMALL = [((2, 1, 0, 0),) * 2, ((2, 1, 1, 0),) * 2, ((3, 2, 1, 0),) * 2, ((4, 2, 0, 0),) * 2, ((4, 2, 1, 0),) * 2, ((3, 2, 1, 0), (4, 2, 1, 0)), ((4, 2, 1, 0), (3, 2, 1, 0)),
           ((4, 2, 2, 0),) * 2, ((4, 3, 1, 0),) * 2]
P_BIG = [((6, 4, 2, 0),) * 2, ((5, 3, 1, 0),) * 2, ((5, 3, 2, 1),) * 2, ((9, 6, 3, 0),) * 2]
P_ROBUST = [((6, 4, 2, 0), (5, 3, 1, 0)), ((7, 4, 1, 0), (4, 2, 2, 0)), ((4, 3, 2, 1), (3, 2, 1, 1)), ((16, 8, 0, 0),) * 2]


def ref_path(lam, lamp):
    return os.path.join(HERE, 'cg_fast_gl4_refcache_%s__%s.pkl' % ('_'.join(map(str, lam)), '_'.join(map(str, lamp))))


def load_refs(lam, lamp):
    import pickle
    p = ref_path(lam, lamp)
    if os.path.exists(p):
        with open(p, 'rb') as f:
            return pickle.load(f)
    return {}


def save_refs(lam, lamp, d):
    import pickle
    with open(ref_path(lam, lamp), 'wb') as f:
        pickle.dump(d, f)


def exact_ref(lam, lamp, mu, refs):
    """G3's exact 50-digit copies reduced to float64: (slices (w x m), triples (m x 3), min joint gap, offmax) -- cached per product in a gitignored pickle"""
    if mu not in refs:
        cps = so.copies(lam, lamp, mu, dps=50)
        F = so.copy_float_slices(cps.prv, cps)
        trip = np.array([[float(x) for x in cp['ev']] for cp in cps])
        refs[mu] = (F, trip, float(cps.info['jgap']), float(cps.info['offmax']))
    return refs[mu]


def t1_product(lam, lamp, verbose=True, twin=False, routeA=False, do_ref=True, hw_sample=None):
    """T1 for one product.  Returns dict(rows=[per sector dict], summary).  The HW vectors are streamed (checked and dropped sector by sector); hw_sample = k checks the HW residual
    of every k-th sector and of the 5 sectors with the largest n_W only (memory / time for [12,8,4,0]^2)."""
    lam, lamp = tuple(lam), tuple(lamp)
    refs = load_refs(lam, lamp) if do_ref else {}
    n0 = len(refs)
    S = Setup.of(lam, lamp)
    rs = gr.racah_speiser(lam, lamp)
    S._rs = rs
    S.rec
    rows = []
    mus = sorted(rs, reverse=True)
    big = set(sorted(mus, key=lambda k: -gr.nW_of(lam, lamp, k))[:5])
    t0 = time.perf_counter()
    for i, mu in enumerate(mus):
        sec = physical_slices(lam, lamp, mu, strict=True)
        r = dict(mu=mu, m=sec['m'], w=sec['w'], w1=sec['w1'], fb=sec['fallback'], s_phys=sec['s_phys'], s_spur=sec['s_spur'], leak=sec['leak'], smallest=sec['smallest'],
                 gapG=sec['gapG'], ncomplex=sec['ncomplex'], nphys=sec['info']['nphys'], nphys_tau=sec['info'].get('nphys_tau', -1), mixed=sec['info'].get('mixed', 0),
                 nstr=len(strings_of(S, mu)), cutactive=sec['w1'] < sec['w'])
        if hw_sample is None or i % hw_sample == 0 or mu in big:
            sec['hw'] = se.hw_from_slice(lam, lamp, mu, sec['slices'])
            res, wd, off = hw_checks(lam, lamp, sec)
            r.update(resid=res, wdev=wd, orth=off)
            sec['hw'] = None
        if do_ref:
            F, trip, jgap, offm = exact_ref(lam, lamp, mu, refs)
            dev = float((np.abs(sec['slices'] - F).max(axis=0) / np.abs(F).max(axis=0)).max())
            tdev = float((np.abs(sec['trip'] - trip) / np.maximum(1.0, np.abs(trip))).max())
            r.update(dev=dev, tdev=tdev, jgap=jgap)
        rows.append(r)
    T = dict(t_total_check=time.perf_counter() - t0)
    if do_ref and len(refs) != n0:
        save_refs(lam, lamp, refs)
    if twin:
        for r in rows:
            d = mats_vs_exact(lam, lamp, r['mu'], 'B')
            r['twin'] = max(v[0] for v in d.values())
            r['twin_leak'] = max(v[1] for v in d.values())
    if routeA:
        ab = {}
        tA = tB = 0.0
        for r in rows:
            t = time.perf_counter()
            W, pos, B, _ = operators_B(S, r['mu'], check=False)
            tB += time.perf_counter() - t
            t = time.perf_counter()
            _, _, A = operators_A(S, r['mu'], ab)
            tA += time.perf_counter() - t
            r['AB'] = max(float(np.abs(A[n] - B[n]).max() / max(np.abs(B[n]).max(), 1.0)) for n in A)
        T['tA'], T['tB'] = tA, tB
        T['q_exact'] = ab.get('q_exact', 0)
    return dict(rows=rows, T=T, ncopies=sum(r['m'] for r in rows), nfb=sum(1 for r in rows if r['fb']))


def t1_report(lam, lamp, res):
    rows = res['rows']
    mx = lambda k: max((r[k] for r in rows if k in r), default=float('nan'))
    spur = [r['s_spur'] for r in rows if np.isfinite(r['s_spur'])]
    line = ('%-22s irreps %4d copies %5d maxw1 %3d | fallbacks %d | #phys=m in all: %s | dev slices %.1e  triples %.1e | HW resid %.1e orth %.1e wdev %d | s_phys max %.1e  s_spur min %s | W1 leak %.1e'
            % (prod_str(lam, lamp), len(rows), res['ncopies'], max(r['w1'] for r in rows), res['nfb'], all(r['nphys'] == r['m'] for r in rows), mx('dev'), mx('tdev'), mx('resid'), mx('orth'),
               max([r['wdev'] for r in rows if 'wdev' in r] + [0]), mx('s_phys'), ('%.2e' % min(spur)) if spur else 'n/a (no spurious)', mx('leak')))
    if any('twin' in r for r in rows):
        line += ' | twin(B vs exact) %.1e (W1-invariance leak of the exact matrices %.1e)' % (mx('twin'), mx('twin_leak'))
    if any('AB' in r for r in rows):
        line += ' | route A vs B %.1e' % mx('AB')
    say(line)
    return line


def stage_t1(extra_big=True):
    hdr('T1  agreement with the existing engine (float64 vs G3 copies at 50 digits; HW residual; G(x)G orthogonality; Route A vs B; exact twin of the operators)')
    all_rows = {}
    for lam, lamp in P_SMALL + P_BIG:
        res = t1_product(lam, lamp, twin=True, routeA=True)
        t1_report(lam, lamp, res)
        all_rows[(lam, lamp)] = res
    return all_rows


def stage_t4(all_rows):
    hdr('T4  precision: deviation of the float64 slices from the 50-digit copies as a function of w1 and of the min joint gap; string-test statistics')
    pts = []
    for (lam, lamp), res in all_rows.items():
        for r in res['rows']:
            if 'dev' in r:
                pts.append((r['w1'], r['m'], r['jgap'], r['dev'], r['smallest'], prod_str(lam, lamp), r['mu']))
    say('  by w1 (max over sectors): ' + '   '.join('w1<=%d: %.1e' % (b, max([p[3] for p in pts if lo < p[0] <= b] + [0.0])) for lo, b in zip((0, 4, 8, 16, 24), (4, 8, 16, 24, 40))))
    say('  by min joint gap of the sector (max dev): ' + '   '.join('gap in (%g,%g]: %.1e (n=%d)' % (lo, hi, max([p[3] for p in pts if lo < p[2] <= hi] + [0.0]), len([p for p in pts if lo < p[2] <= hi]))
                                                           for lo, hi in ((0, 0.1), (0.1, 1), (1, 10), (10, 100), (100, 1e9))))
    say('  worst 8 sectors (dev, w1, m, min joint gap, smallest relative slice component, product, mu):')
    for p in sorted(pts, key=lambda p: -p[3])[:8]:
        say('    %.2e  w1=%d m=%d gap=%.3g smallest=%.1e %s mu=%s' % (p[3], p[0], p[1], p[2], p[4], p[5], p[6]))
    say('  string test per product: min over spurious / max over physical of s(u) (orthonormal basis, 2-norm, tau = %g):' % TAU)
    for (lam, lamp), res in all_rows.items():
        sp = [r['s_spur'] for r in res['rows'] if np.isfinite(r['s_spur'])]
        ph = [r['s_phys'] for r in res['rows']]
        say('    %-20s min spurious %.3e   max physical %.3e   gap factor %.1e' % (prod_str(lam, lamp), min(sp) if sp else float('nan'), max(ph), (min(sp) / max(max(ph), 1e-300)) if sp else float('nan')))


def _t5_one(lam, mu):
    t = time.perf_counter()
    sec = physical_slices(lam, lam, mu, strict=True)
    t_f = time.perf_counter() - t
    t = time.perf_counter()
    refs = load_refs(lam, lam)
    F, trip, jgap, offm = exact_ref(lam, lam, mu, refs)
    t_e = time.perf_counter() - t
    save_refs(lam, lam, refs)
    dev = float((np.abs(sec['slices'] - F).max(axis=0) / np.abs(F).max(axis=0)).max())
    d34 = min(max(abs(a[0] - b[0]), abs(a[1] - b[1])) for a, b in itertools.combinations(sec['trip'], 2))
    say('  [12,8,4,0]^2, mu = %s: m = %d , w = %d, w1 = %d; min joint gap of the triple %.3g; min (V3,V4a) gap %.3g; min G gap %.3g; ritz combination %s; fallback %s; '
        'dev slices %.1e; dev triples %.1e   [float %.2f s, exact 50-digit copies of this sector %s]'
        % (mu, sec['m'], sec['w'], sec['w1'], jgap, d34, sec['gapG'], sec['info'].get('ritz_combo'), sec['fallback'], dev, float((np.abs(sec['trip'] - trip) / np.maximum(1, np.abs(trip))).max()), t_f, t_e))


def stage_t5(all_rows=None):
    hdr('T5  robustness: degenerate (V3, V4a) sectors, lambda != lambda\', lambda_4 > 0, void / full string sectors, complex-pair guard')
    all_rows = dict(all_rows or {})
    # (a) the two (V3, V4a)-degenerate sectors
    lam = (6, 4, 2, 0)
    mu = (9, 7, 5, 3)
    refs = load_refs(lam, lam)
    sec = physical_slices(lam, lam, mu, strict=True)
    F, trip, jgap, offm = exact_ref(lam, lam, mu, refs)
    save_refs(lam, lam, refs)
    dev = float((np.abs(sec['slices'] - F).max(axis=0) / np.abs(F).max(axis=0)).max())
    d34 = min(max(abs(a[0] - b[0]), abs(a[1] - b[1])) for a, b in itertools.combinations(sec['trip'], 2))
    d3 = min(abs(a[0] - b[0]) for a, b in itertools.combinations(sec['trip'], 2))
    say('  [6,4,2,0]^2, mu = (9,7,5,3): m = %d, w = %d, w1 = %d; min joint gap of the triple %.3g; min (V3,V4a) gap %.3g (V3 alone: %.3g); min G gap %.3g; ritz combination %s; fallback %s; dev slices %.1e; dev triples %.1e'
        % (sec['m'], sec['w'], sec['w1'], jgap, d34, d3, sec['gapG'], sec['info'].get('ritz_combo'), sec['fallback'], dev, float((np.abs(sec['trip'] - trip) / np.maximum(1, np.abs(trip))).max())))
    lam = (12, 8, 4, 0)
    for mu in ((18, 14, 10, 6), (19, 14, 10, 5)):
        _t5_one(lam, mu)
    # (b) lambda != lambda', lambda_4 > 0, lambda_3 = lambda_4 = 0 at large degree; and the whole of [12,8,4,0]^2 without an exact reference (slices vs exact only for the two sectors above)
    lam = (12, 8, 4, 0)
    res = t1_product(lam, lam, do_ref=False, hw_sample=15)
    all_rows[(lam, lam)] = res
    t1_report(lam, lam, res)
    say('  ([12,8,4,0]^2: HW residual / orthogonality checked on every 15th sector and the 5 largest n_W only; no exact reference except the two sectors above)')
    for lam, lamp in P_ROBUST:
        res = t1_product(lam, lamp, twin=False)
        all_rows[(lam, lamp)] = res
        t1_report(lam, lamp, res)
    # (c) classes of sectors over every product run so far
    rows = [r for res in all_rows.values() for r in res['rows']]
    def cat(f):
        sub = [r for r in rows if f(r)]
        return '%d sectors (%d copies), max dev %.1e, max dev triples %.1e, fallbacks %d' % (len(sub), sum(r['m'] for r in sub), max([r.get('dev', 0.0) for r in sub] + [0.0]),
                                                                                           max([r.get('tdev', 0.0) for r in sub] + [0.0]), sum(1 for r in sub if r['fb']))
    say('  sectors with a void E12 cut (W1 = W):            ' + cat(lambda r: not r['cutactive']))
    say('  sectors with an active E12 cut (W1 < W):         ' + cat(lambda r: r['cutactive']))
    for k in (0, 1, 2):
        say('  sectors with %d non-void strings i = 2, 3:         %s' % (k, cat(lambda r, k=k: r['nstr'] == k)))
    say('  every string active (E12 cut active and both i = 2, 3 non-void), |A| = 3: ' + cat(lambda r: r['cutactive'] and r['nstr'] == 2))
    nc = [r for r in rows if r['ncomplex'] > 0]
    say('  complex-pair guard (eig of G on W1 returned complex eigenvalues; physical eigenvectors must be real, else fallback): complex eigenvalues in %d of %d sectors (%d eigenvalues in total); fallbacks among them %d'
        % (len(nc), len(rows), sum(r['ncomplex'] for r in nc), sum(1 for r in nc if r['fb'])))
    say('  fallbacks in all runs of T1 + T5: %d of %d sectors' % (sum(1 for r in rows if r['fb']), len(rows)))
    say('  the brief\'s criterion alone (number of eigenvectors with row-normalised s < tau = %g equals m): holds in %d of %d sectors; eigenvectors of G at distance 1e-8..0.1 from the physical subspace (mixed): %d'
        % (TAU, sum(1 for r in rows if r['nphys_tau'] == r['m']), len(rows), sum(r['mixed'] for r in rows)))
    return all_rows


def stage_t2():
    """T2: Munich literal comparison with the new engine as source.  Cached raw Munich outputs of G7 (munich_bridge.run_hwcgc reuse=True); nothing of Munich is re-run unless a file is missing."""
    import munich_bridge as mb
    hdr('T2  Munich literal comparison (cached G7 outputs; new engine as the source; MATCH <= 1e-9 after one global sign, sectors: subspace max|1 - sv| <= 1e-9)')
    cache = {}

    def ours(lam, lamp, mu):
        key = (lam, lamp)
        if key not in cache:
            cache[key] = hw_vectors(lam, lamp)
        return cache[key][mu]

    missing = []
    # 2a: multiplicity-free, vector by vector
    ndev, nok = [], 0
    for lam, lamp, mu in mb.mf_list():
        mm = mb.munich_label(mu)
        fname = 'hw_%d_%s.txt' % (lam[0], '_'.join(map(str, mm)))
        path = os.path.join(mb.OUT_DIR, fname)
        sec = ours(lam, lamp, mu) if True else None
        assert sec['m'] == 1
        if not os.path.exists(path):          # G7: hwcgc.sh timed out on the S = 8 top irrep; the repaired profiler's output stands in
            btag = '%s__%s__%s' % ('_'.join(map(str, lam)), '_'.join(map(str, lamp)), '_'.join(map(str, mm)))
            bpath = os.path.join(mb.OUT_DIR, 'bench_eig_%s.txt' % btag)
            if not os.path.exists(bpath):
                missing.append(fname)
                continue
            r = sector_vs_texts(lam, lamp, mu, sec, parse_bench(bpath, gr.weyl_dim(mu)))
            nok += r['status'] == 'OK'
            ndev.append(r['maxdev'])
            say('   2a (%s from the repaired profiler output %s): %s, max |1 - sv| = %.2e' % (prod_str(lam, lamp), os.path.basename(bpath), r['status'], r['maxdev']))
            continue
        ok, dev, nsup = compare_mf_file(lam, lamp, mu, sec, path)
        nok += ok
        ndev.append(dev)
        if not ok:
            say('   2a FAIL %s mu=%s support %s dev %.2e' % (prod_str(lam, lamp), mu, nsup, dev))
    say('  2a (multiplicity-free, vector by vector): %d MATCH of %d  (missing cached files: %d); max deviation after one global sign %.2e (tol %g)'
        % (nok, len(ndev), len(missing), max(ndev + [0.0]), mb.TOL_O))
    # 2b: sectors by subspace
    res = []
    for lam, lamp in mb.PRODUCTS:
        for mu, m in sorted(gr.racah_speiser(lam, lamp).items(), reverse=True):
            mm = mb.munich_label(mu)
            tag = '%s__%s__%s' % ('_'.join(map(str, lam)), '_'.join(map(str, lamp)), '_'.join(map(str, mm)))
            texts = []
            for a in range(1, m + 1):
                path = os.path.join(mb.OUT_DIR, 'hw_%s_a%d.txt' % (tag, a))
                if not os.path.exists(path):
                    texts = None
                    break
                texts.append(_read(path))
            if texts is None:
                missing.append('hw_%s_a*.txt' % tag)
                continue
            res.append((lam, lamp, mu, m, sector_vs_texts(lam, lamp, mu, ours(lam, lamp, mu), texts)))
    bad = [r for r in res if r[4]['status'] != 'OK']
    say('  2b (multiplicity sectors of the 8 G7 products, subspace by subspace): %d sectors, %d OK, sum m = %d, max m = %d; worst max|1 - sv| = %.2e (tol %g); worst Munich orthonormality %.1e; worst residual %.1e'
        % (len(res), len(res) - len(bad), sum(r[3] for r in res), max([r[3] for r in res] + [0]), max([r[4]['maxdev'] for r in res] + [0.0]), 1e-9, max([r[4]["orth_munich"] for r in res] + [0.0]),
           max([r[4]['resid'] for r in res] + [0.0])))
    # [6,4,2,0]^2 big sectors from the repaired profiler outputs
    lam = (6, 4, 2, 0)
    big = []
    for mu in mb.BIG:
        mm = mb.munich_label(mu)
        tag = '%s__%s__%s' % ('_'.join(map(str, lam)), '_'.join(map(str, lam)), '_'.join(map(str, mm)))
        path = os.path.join(mb.OUT_DIR, 'bench_eig_%s.txt' % tag)
        if not os.path.exists(path):
            missing.append(os.path.basename(path))
            continue
        texts = parse_bench(path, gr.weyl_dim(mu))
        r = sector_vs_texts(lam, lam, mu, ours(lam, lam, mu), texts)
        big.append((mu, gr.racah_speiser(lam, lam)[mu], r))
    say('  [6,4,2,0]^2, the %d sectors G7 ran through the repaired profiler: %s' % (len(big), '; '.join('mu=%s m=%d maxdev %.1e %s' % (mu, m, r['maxdev'], r['status']) for mu, m, r in big)))
    say('  missing cached Munich files (not re-run): %s' % (missing if missing else 'none'))
    return dict(nok2a=nok, n2a=len(ndev), res2b=res, big=big, missing=missing)


def compare_mf_file(lam, lamp, mu, sec, path, tol=1e-9):
    """multiplicity-free irrep: our HW vector (Munich's GT gauge, se.to_munich) against a cached hwcgc.sh output, same support (|C| > tol) and equal after ONE global sign to tol.  Returns (ok, dev, (n_ours, n_munich))"""
    import munich_bridge as mb
    rows = mb.hw_rows(_read(path), gr.weyl_dim(mu))
    o = se.to_munich(lam, lamp, sec['hw'], 0)
    keys_o = {k for k, v in o.items() if abs(v) > tol}
    keys_m = {k for k, v in rows.items() if abs(v) > tol}
    kbig = max(o, key=lambda k: abs(o[k]))
    sigma = 1.0 if rows.get(kbig, 0.0) * o[kbig] >= 0 else -1.0
    dev = max([abs(sigma * o.get(k, 0.0) - rows.get(k, 0.0)) for k in keys_o | keys_m] + [0.0])
    return (keys_o == keys_m and dev <= tol), dev, (len(keys_o), len(keys_m))


def parse_bench(path, dimmu):
    """the per-alpha texts of a repaired-profiler output (cgc_hw_prof_ob, HW_MODE=eig; munich_bridge.bench_texts' parser)"""
    out_ = _read(path)
    texts, cur = [], None
    for line in out_.splitlines():
        if re.match(r'\s*alpha = \d+:', line):
            cur = []
            texts.append(cur)
        elif cur is not None:
            f = line.split()
            if len(f) == 3 and all(re.fullmatch(r'-?\d+', x) for x in f[:2]):
                cur.append('%s\t%s\t%d\t%s' % (f[0], f[1], dimmu, f[2]))
    return ['\n'.join(t) + '\n' for t in texts]


def sector_vs_texts(lam, lamp, mu, sec, texts, tol_sv=1e-9):
    """subspace comparison of our sector (the Sector of this module) with m Munich outputs (G7's compare_sector, second half)"""
    import munich_bridge as mb
    m = sec['m']
    dimmu = gr.weyl_dim(mu)
    Mrows = [mb.hw_rows(t, dimmu) for t in texts]
    O = [se.to_munich(lam, lamp, sec['hw'], k) for k in range(m)]
    keys = {}
    for d in Mrows + O:
        for k in d:
            keys.setdefault(k, len(keys))
    Mmat = np.zeros((len(keys), m))
    Omat = np.zeros((len(keys), m))
    for j in range(m):
        Mmat[:, j] = mb.sparse_from_rows(Mrows[j], keys)[0]
        Omat[:, j] = mb.sparse_from_rows(O[j], keys)[0]
    orth_m = float(np.abs(Mmat.T @ Mmat - np.eye(m)).max())
    Onorm = Omat / np.linalg.norm(Omat, axis=0)
    Q = Onorm.T @ Mmat
    sv = np.linalg.svd(Q, compute_uv=False)
    maxdev = float(np.abs(1 - sv).max())
    resid = float(np.linalg.norm(Mmat - Onorm @ Q))
    ok = maxdev <= tol_sv and orth_m <= mb.TOL_O and resid <= 1e-7
    return dict(status='OK' if ok else 'FAIL', maxdev=maxdev, orth_munich=orth_m, resid=resid, sv=sv, Q=Q)


# ======================================================================================== CLI
def _fmt_lam(lam):
    return '[' + ','.join(str(x) for x in lam) + ']'


def run_cli(args):
    lam = tuple(int(x) for x in args[:4])
    rest = args[4:]
    lamp = lam
    if len(rest) >= 4 and all(r.lstrip('-').isdigit() for r in rest[:4]):
        lamp = tuple(int(x) for x in rest[:4])
        rest = rest[4:]
    mus = None
    outf = None
    ortho = '--orthonormal' in rest
    desc = '--descendants' in rest
    for i, a in enumerate(rest):
        if a == '--mus':
            mus = [tuple(int(x) for x in s.split(',')) for s in rest[i + 1].split(';')]
        if a == '--out':
            outf = rest[i + 1]
    T = {}
    t0 = time.perf_counter()
    secs = hw_vectors(lam, lamp, mus, orthonormal=ortho, timing=T)
    dt = time.perf_counter() - t0
    ncop = sum(s['m'] for s in secs.values())
    nfb = sum(1 for s in secs.values() if s['fallback'])
    print('%s x %s: %d irreps, %d copies, %d fallbacks; t_setup %.3f t_ops %.3f t_eig %.3f t_hw %.3f total %.3f s' %
          (_fmt_lam(lam), _fmt_lam(lamp), len(secs), ncop, nfb, T['t_setup'], T['t_ops'], T['t_eig'], T['t_hw'], dt))
    if outf:
        write_vectors(outf, lam, lamp, secs, orthonormal=ortho)
        print('written', outf)
    if desc:
        sec = next(iter(secs.values()))
        Rm, cols = descendants(lam, lamp, sec['mu'], sec, 0)
        print('descendants of the first sector: %d of %d states of V_mu' % (len(cols), Rm.dim))
    return secs


# ======================================================================================== T3: same-machine benchmark
def loadavg():
    try:
        return open('/proc/loadavg').read().split()[:3]
    except OSError:
        return ['n/a']


def cold():
    Setup.clear()
    import gc
    gc.collect()


def bench_ours(lam, lamp, reps=5, route='B', lean=False):
    """cold runs of the whole float pipeline; columns of the minimum-total run.  t_setup = GLNRep x2 (+ float generators of the adjacent / diagonal operators, reductions, Racah-Speiser, G1's
    left inverses); t_ops = S1; t_eig = S2 (eig + string images + string test + null space + Ritz + canonical normalisation); t_hw = S3"""
    best = None
    done = 0
    while done < reps:
        done += 1
        cold()
        T = {}
        t0 = time.perf_counter()
        secs = hw_vectors(lam, lamp, timing=T, route=route, keep_hw=False, lean=lean)
        T['total'] = time.perf_counter() - t0
        T['nfb'] = sum(1 for x in secs.values() if x['fallback'])
        if T['total'] > 20.0:
            reps = min(reps, 2)           # long runs: min of 2
        T['reps'] = reps
        if best is None or T['total'] < best['total']:
            best = T
            best['nirrep'] = len(secs)
            best['ncopies'] = sum(x['m'] for x in secs.values())
            best['maxw1'] = max(x['w1'] for x in secs.values())
            best['maxw'] = max(x['w'] for x in secs.values())
            best['nent'] = sum(x['nent'] for x in secs.values())
    return best


def bench_routeA(lam, lamp, budget=60.0):
    """warm time of route A (V3, V4a closed forms, float, entry-by-entry Python; the Q classes through G4's exact evaluator) over all sectors, with a time budget"""
    S = Setup.of(tuple(lam), tuple(lamp))
    st = {}
    t0 = time.perf_counter()
    n = 0
    mus = sorted(gr.racah_speiser(lam, lamp), reverse=True)
    for mu in mus:
        operators_A(S, mu, st)
        n += 1
        if time.perf_counter() - t0 > budget:
            break
    return dict(t=time.perf_counter() - t0, n=n, of=len(mus), q_exact=st.get('q_exact', 0))


class _Timeout(Exception):
    pass


def bench_g3(lam, lamp, reps=1, exact=True, flt=True, limit=900):
    """the current engine (G1 + G3), same venv, cold: setup, PRV spaces (exact), the three operators on the PRV basis (EXACT, Fractions), joint diagonalisation (mpmath 50 digits and
    float64, from the SAME operator matrices of the run), G1's float level recursion.  n_W-free part = t_prv + t_ops + t_diag.  limit = seconds per run (SIGALRM); returns {'exact': T, 'float': T}
    with T = None when the limit was hit"""
    import signal

    def handler(signum, frame):
        raise _Timeout()
    signal.signal(signal.SIGALRM, handler)
    best = {'exact': None, 'float': None}
    for _ in range(reps):
        cold()
        T = {}
        signal.alarm(limit)
        try:
            t0 = time.perf_counter()
            R1 = GLNRep.of(lam)
            R2 = GLNRep.of(lamp)
            rs = gr.racah_speiser(lam, lamp)
            mus = sorted(rs, reverse=True)
            rec = se.Recursion.of(lam, lamp)
            T['t_setup'] = time.perf_counter() - t0
            t = time.perf_counter()
            prvs = {mu: se.prv_space(lam, lamp, mu) for mu in mus}
            T['t_prv'] = time.perf_counter() - t
            t = time.perf_counter()
            mats = {mu: so.triple_mats(lam, lamp, mu, prvs[mu])[1] for mu in mus}
            T['t_ops'] = time.perf_counter() - t
            Ff = {}
            t = time.perf_counter()
            for mu in mus:
                jd = so.joint_diag_float([mats[mu][n] for n in so.TRIPLE])
                Ff[mu] = prvs[mu].float_basis() @ jd['U']
            T['t_diag_float'] = time.perf_counter() - t
            t = time.perf_counter()
            for mu in mus:
                se.hw_from_slice(lam, lamp, mu, Ff[mu])
            T['t_hw'] = time.perf_counter() - t
            Fe = {}
            t = time.perf_counter()
            if exact:
                for mu in mus:
                    cps = so.copies(lam, lamp, mu, dps=50, prv=prvs[mu], mats=mats[mu])
                    Fe[mu] = so.copy_float_slices(prvs[mu], cps)
            T['t_diag_exact'] = time.perf_counter() - t
        except _Timeout:
            T['timeout'] = True
        finally:
            signal.alarm(0)
        for v, key in (('exact', 't_diag_exact'), ('float', 't_diag_float')):
            if (v == 'exact' and not exact) or (v == 'float' and not flt) or key not in T:
                continue
            R = dict(t_setup=T['t_setup'], t_prv=T['t_prv'], t_ops=T['t_ops'], t_diag=T[key], t_hw=T['t_hw'])
            R['t_slices'] = R['t_prv'] + R['t_ops'] + R['t_diag']
            R['total'] = R['t_setup'] + R['t_slices'] + R['t_hw']
            if best[v] is None or R['total'] < best[v]['total']:
                best[v] = R
        if T.get('timeout'):
            break
    return best


def _read(path):
    with open(path) as f:
        return f.read()


def _lines(path):
    with open(path) as f:
        return f.read().splitlines()


def _num(x):
    try:
        return float(x)
    except ValueError:
        return None


def munich_numbers():
    """Munich same-machine numbers of G7 (orig = cgc_hw_bench, eig = repaired profiler, 1 thread): rectangular [S,S,0,0]^2 from rect_table.md, products from products_summary.txt"""
    d = os.path.normpath(os.path.join(HERE, '..', '..', 'Cpp', 'AlexCGC', 'results', 'su4-same-machine'))
    out_ = {}
    for line in _lines(os.path.join(d, 'rect_table.md')):
        f = [x.strip() for x in line.strip().strip('|').split('|')]
        if len(f) >= 7 and f[0].isdigit():
            S = int(f[0])
            num = _num
            origs = [num(f[2]), num(f[3])]
            origs = [x for x in origs if x is not None]
            eigs = [num(f[4]), num(f[5])]
            eigs = [x for x in eigs if x is not None]
            out_[((S, S, 0, 0), (S, S, 0, 0))] = dict(orig=min(origs) if origs else None, eig=min(eigs) if eigs else None)
    for line in _lines(os.path.join(d, 'products_summary.txt')):
        if line.startswith('('):
            f = [x.strip() for x in line.split('|')]
            lam = tuple(int(x) for x in f[0].strip('()').split(','))
            num = _num
            out_[(lam, lam)] = dict(orig=num(f[6]), eig=num(f[7]), partial=(f[5], f[8]), nirrep=int(f[1]))
    return out_


BENCH_LIST = ([((S, S, 0, 0),) * 2 for S in (4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 20)] +
              [((4, 2, 1, 0),) * 2, ((6, 3, 0, 0),) * 2, ((6, 4, 2, 0),) * 2, ((8, 4, 0, 0),) * 2, ((9, 6, 3, 0),) * 2, ((12, 8, 4, 0),) * 2,
               ((15, 10, 5, 0),) * 2, ((16, 8, 0, 0),) * 2, ((12, 12, 0, 0),) * 2])


def stage_benchlean(products=None, reps=5):
    hdr('T3b  lean mode (multiplicity-one sectors: no operator applied, slice = null vector of the string powers; triples of those sectors not computed), cold, min of %d  [%s]' % (reps, time.strftime('%Y-%m-%d %H:%M')), to='t')
    say('LC_ALL=%s OPENBLAS_NUM_THREADS=%s OMP_NUM_THREADS=%s' % (os.environ.get('LC_ALL'), os.environ.get('OPENBLAS_NUM_THREADS'), os.environ.get('OMP_NUM_THREADS')), to='t')
    bench_ours((3, 2, 1, 0), (3, 2, 1, 0), reps=1, lean=True)
    for lam, lamp in (products or BENCH_LIST):
        say('--- %s   loadavg before: %s' % (prod_str(lam, lamp), ' '.join(loadavg())), to='t')
        o = bench_ours(lam, lamp, reps=reps, lean=True)
        say('  ours lean (min of %d)  t_setup %.4f  t_ops %.4f  t_eig %.4f  | t_slices %.4f | t_hw %.4f | t_total %.4f   [irreps %d copies %d  fallbacks %d]' %
            (o['reps'], o['t_setup'], o['t_ops'], o['t_eig'], o['t_slices'], o['t_hw'], o['total'], o['nirrep'], o['ncopies'], o['nfb']), to='t')


def stage_bench(products=None, reps=5, g3_reps=None):
    hdr('T3  same-machine benchmark: ours (float64, cold, min of %d) vs G3/G1 (exact 50-digit and float64) vs Munich (G7 numbers)   [%s]' % (reps, time.strftime('%Y-%m-%d %H:%M')), to='t')
    say('LC_ALL=%s OPENBLAS_NUM_THREADS=%s OMP_NUM_THREADS=%s; python %s, numpy %s; machine %s' % (os.environ.get('LC_ALL'), os.environ.get('OPENBLAS_NUM_THREADS'), os.environ.get('OMP_NUM_THREADS'),
                                                                                              platform.python_version(), np.__version__, platform.platform()), to='t')
    mun = munich_numbers()
    bench_ours((3, 2, 1, 0), (3, 2, 1, 0), reps=1)          # warm-up of imports / BLAS
    G3_LIMIT = 600
    rows = []
    for lam, lamp in (products or BENCH_LIST):
        say('', to='t')
        say('--- %s   loadavg before: %s' % (prod_str(lam, lamp), ' '.join(loadavg())), to='t')
        o = bench_ours(lam, lamp, reps=reps)
        say('  ours  (min of %d)  t_setup %.4f  t_ops %.4f  t_eig %.4f  | t_slices %.4f | t_hw %.4f | t_total %.4f   [irreps %d copies %d  max w1 %d  entries %d  fallbacks %d]' %
            (o['reps'], o['t_setup'], o['t_ops'], o['t_eig'], o['t_slices'], o['t_hw'], o['total'], o['nirrep'], o['ncopies'], o['maxw1'], o['nent'], o['nfb']), to='t')
        A = bench_routeA(lam, lamp, budget=45.0)
        say('  route A (V3, V4a closed forms, warm): %.3f s for %d of %d sectors (%d Q-class entries through the exact evaluator)%s' %
            (A['t'], A['n'], A['of'], A['q_exact'], '' if A['n'] == A['of'] else '  [time budget hit]'), to='t')
        g_reps = g3_reps if g3_reps is not None else (3 if o['total'] < 0.5 else 1)
        exact_ok = True
        g = bench_g3(lam, lamp, reps=g_reps, exact=exact_ok, flt=True, limit=G3_LIMIT if lam[0] >= 15 or lam == (12, 12, 0, 0) else 900)
        for v in ('exact', 'float'):
            T = g[v]
            if T is None:
                say('  G3/G1 %-5s: not finished within the limit' % v, to='t')
                continue
            say('  G3/G1 %-5s (min of %d)  t_setup %.4f  t_prv %.4f  t_ops(exact) %.4f  t_diag %.4f  | n_W-free %.4f | t_hw %.4f | total %.4f' %
                (v, g_reps, T['t_setup'], T['t_prv'], T['t_ops'], T['t_diag'], T['t_slices'], T['t_hw'], T['total']), to='t')
        m = mun.get((lam, lamp), {})
        ms = []
        for k in ('orig', 'eig'):
            if m.get(k) is not None:
                ms.append('Munich %s %.4g s (x%.0f of our t_total)' % (k, m[k], m[k] / o['total']))
        sp = []
        for v in ('exact', 'float'):
            if g[v] is not None:
                sp.append('vs G3 %s: t_slices x%.1f (n_W-free), t_total x%.1f' % (v, g[v]['t_slices'] / o['t_slices'], g[v]['total'] / o['total']))
        say('  speed-ups of ours: %s; t_hw/t_total = %.2f;  %s' % ('; '.join(sp), o['t_hw'] / o['total'], '; '.join(ms) if ms else '(no Munich number: not run)'), to='t')
        rows.append(dict(lam=lam, lamp=lamp, o=o, g=g, A=A, m=m))
    say('', to='t')
    say('loadavg after: %s' % ' '.join(loadavg()), to='t')
    return rows


def stage_t4s(Slist=(4, 8, 12, 14, 16, 20), nworst=2, rounded_from=14):
    """T4 large-S: worst HW residual of the pipeline (float slices + G1's float recursion) over [S,S,0,0]^2, and for S >= rounded_from the nworst worst irreps against G1's exact recursion
    rounded to double (hw_from_slice_rounded): residuals, time, max relative difference of the two vectors"""
    hdr('T4  large S: HW residual of the float64 pipeline on [S,S,0,0]^2 and against the exact recursion rounded to double')
    for S in Slist:
        lam = (S, S, 0, 0)
        res = t1_product(lam, lam, do_ref=False)
        rows = res['rows']
        worst = max(rows, key=lambda r: r['resid'])
        say('  S = %2d: %3d irreps, worst HW residual %.1e (mu = %s), worst orthogonality %.1e, weight deviation %d, string statistic s_phys max %.1e / s_spur min %.2e, fallbacks %d'
            % (S, len(rows), worst['resid'], worst['mu'], max(r['orth'] for r in rows), max(r['wdev'] for r in rows), max(r['s_phys'] for r in rows), min([r['s_spur'] for r in rows if np.isfinite(r['s_spur'])] + [np.inf]), res['nfb']))
        if S >= rounded_from:
            for r in sorted(rows, key=lambda r: -r['resid'])[:nworst]:
                mu = r['mu']
                sec = physical_slices(lam, lam, mu)
                t = time.perf_counter()
                hw = se.hw_from_slice(lam, lam, mu, sec['slices'])
                tf = time.perf_counter() - t
                rf, _ = se.hw_residual(hw)
                prv = se.prv_space(lam, lam, mu)
                t = time.perf_counter()
                hr = se.hw_from_slice_rounded(lam, lam, mu, prv.basis[0])
                te = time.perf_counter() - t
                rr, _ = se.hw_residual(hr)
                a_, b_ = hw.vals_hat[:, 0], hr.vals_hat[:, 0]
                if np.dot(a_, b_) < 0:
                    b_ = -b_
                scale = max(float(np.abs(b_).max()), 1e-300)
                dev = float(np.abs(a_ / np.linalg.norm(a_) - b_ / np.linalg.norm(b_)).max())
                say('       mu = %-16s n_W %7d: float recursion %.1e (%.3f s)   exact recursion rounded to double %.1e (%.2f s)   max |difference| of the two unit vectors %.1e'
                    % (mu, hw.nW, rf, tf, rr, te, dev))


def stage_t6(products=None):
    """T6: the descendants (full CG tables) against Munich's menu-5 tables (all rows, not only the HW rows) for the small products; unitarity of the whole CG matrix"""
    import munich_bridge as mb
    hdr('T6  descendants: full CG tables vs Munich (clebsch menu 5, every row) and unitarity of the CG matrix')
    os.makedirs(RAWDIR, exist_ok=True)
    for lam, lamp in (products or [((2, 1, 0, 0),) * 2, ((2, 1, 1, 0),) * 2, ((3, 2, 1, 0),) * 2]):
        t0 = time.perf_counter()
        secs = hw_vectors(lam, lamp)
        t1 = time.perf_counter()
        C = cg_table(lam, lamp, secs)
        t2 = time.perf_counter()
        n = C.shape[0]
        uni = max(abs((C @ C.T) - sps.identity(n)).max(), abs((C.T @ C) - sps.identity(n)).max())
        n1, n2 = gr.weyl_dim(lam), gr.weyl_dim(lamp)
        worst = 0.0
        nstates = 0
        nsec = 0
        for mu in sorted(secs, reverse=True):
            sec = secs[mu]
            m = sec['m']
            dimmu = gr.weyl_dim(mu)
            mm = mb.munich_label(mu)
            tabs = []
            for alpha in range(1, m + 1):
                path = os.path.join(RAWDIR, 'full_%s__%s__%s_a%d.txt' % ('_'.join(map(str, lam)), '_'.join(map(str, lamp)), '_'.join(map(str, mm)), alpha))
                if os.path.exists(path):
                    # cached raw output
                    txt = _read(path)
                    tab = {}
                    for line in txt.splitlines():
                        f = line.replace(',', ' ').split()
                        if len(f) == 4 and all(re.fullmatch(r'\d+', x) for x in f[:3]):
                            tab.setdefault(int(f[2]), {})[(int(f[0]), int(f[1]))] = float(f[3])
                else:
                    tab = munich_full_table(lam, lamp, mu, alpha)
                    with open(path, 'w') as f:
                        for q, d in sorted(tab.items()):
                            for (a_, b_), v in sorted(d.items()):
                                f.write('%d\t%d\t%d\t%.15f\n' % (a_, b_, q, v))
                tabs.append(tab)
            # rotation from the HW rows (Q = Munich dimmu rows)
            O = [descendants(lam, lamp, mu, sec, k) for k in range(m)]
            keys_all = {}
            for tab in tabs:
                for d in tab.values():
                    for k_ in d:
                        keys_all.setdefault(k_, len(keys_all))
            # our columns as dicts
            def ours_q(k, q):
                col = O[k][1][q]
                coo = col.tocoo()
                return {(int(r) // n2 + 1, int(r) % n2 + 1): float(v) for r, v in zip(coo.row, coo.data) if abs(v) > 1e-14}
            for k in range(m):
                for q in range(dimmu):
                    for key in ours_q(k, q):
                        keys_all.setdefault(key, len(keys_all))
            nk = len(keys_all)

            def mat(getter, q):
                A = np.zeros((nk, m))
                for j in range(m):
                    for key, v in getter(j, q).items():
                        A[keys_all[key], j] = v
                return A
            mun = lambda j, q: tabs[j].get(q + 1, {})
            ou = lambda j, q: ours_q(j, q)
            Mtop, Otop = mat(mun, dimmu - 1), mat(ou, dimmu - 1)
            Qr = Otop.T @ Mtop
            if m > 1:
                U_, _, Vt_ = np.linalg.svd(Qr)
                Qr = U_ @ Vt_
            else:
                Qr = np.array([[1.0 if Qr[0, 0] >= 0 else -1.0]])
            for q in range(dimmu):
                dev = np.abs(mat(mun, q) - mat(ou, q) @ Qr).max()
                worst = max(worst, float(dev))
                nstates += m
            nsec += 1
        say('  %-14s CG matrix %d x %d (dims %d x %d), built in %.2f s (HW part %.2f s): unitarity max|C C^T - 1|, |C^T C - 1| = %.1e; vs Munich menu-5 tables: %d irreps, %d (copy, state) vectors, max deviation after the HW-fixed rotation %.2e'
            % (prod_str(lam, lamp), n, n, n1, n2, t2 - t1, t1 - t0, uni, nsec, nstates, worst))


def main(a):
    if a[0] == 'run':
        return run_cli(a[1:])
    if a[0] == 'validate':
        rows = stage_t1()
        stage_t4(rows)
        stage_t5(rows)
        return
    if a[0] == 't1':
        rows = stage_t1()
        stage_t4(rows)
        return
    if a[0] == 't2':
        return stage_t2()
    if a[0] == 't5':
        return stage_t5()
    if a[0] == 't6':
        return stage_t6()
    if a[0] == 't4s':
        return stage_t4s()
    if a[0] == 'benchlean':
        return stage_benchlean([(tuple(int(x) for x in t.split(',')),) * 2 for t in a[1:]] or None)
    if a[0] == 'bench':
        sel = None
        if len(a) > 1:
            sel = [(tuple(int(x) for x in t.split(',')),) * 2 for t in a[1:]]
        return stage_bench(sel)
    print(__doc__)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
    else:
        main(sys.argv[1:])
