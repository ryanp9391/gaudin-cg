"""glN_reps.py -- exact sparse Gelfand-Tsetlin (GT) representations of gl(N), any N (WP-G1 of gl4-slice-engine-plan.md).

The general-N sibling of  Paul/Python/gl3_reps.GL3Rep  and of  Paul/Python/gaudin_sl3/slice_engine.Rep.  The GT formulas are NOT re-derived: the four
functions gt_patterns / _raise / _lower / _weight of gl3_reps are general-N (they only use the row lengths) and are imported read-only.

Conventions (binding; identical to gl3_reps / slice_engine wherever those had one)
  * a GT pattern of gl(N) = tuple of N rows, rows[0] = the singleton, rows[N-1] = lam  (gl3_reps convention);
  * the basis order is the order of gl3_reps.gt_patterns(lam) (itertools.product order); the HIGHEST-WEIGHT state is the LAST index (asserted);
  * generators: E_{k,k+1} from _raise(k, .), E_{k+1,k} from _lower(k, .) (targets outside the pattern set dropped, never clipped), E_kk = diag(_weight(k, .)),
    E_ij = [E_{i,j-1}, E_{j-1,j}] and E_ji = [E_{j,j-1}, E_{j-1,i}] for j > i+1  (N = 3: E13 = [E12,E23], E31 = [E32,E21]);  the commutators ARE the
    generators (the GT matrices form a representation), so any path gives the same matrix -- tested;
  * plain rational non-unitary GT normalisation (Fractions); matrices act on kets as columns, E e_col = sum_row M[row,col] e_row; E[(i,j)] = list over the
    columns of dicts {row: Fraction}  (the slice_engine.Rep layout; zero entries are not stored; the diagonal E[(k,k)] is stored the same way);
  * weights = (E_11,...,E_NN) eigenvalues, simple roots alpha_i = e_i - e_{i+1};
  * invariant metric: the diagonal G with E_ij^T G = G E_ji, G[HW] = 1 (walk from the HW state through the simple raising operators); orthonormal basis
    e~_q = e_q / sqrt(g_q), i.e. Eh[r,q] = E[r,q] sqrt(g_r/g_q) and Eh_ji = Eh_ij^T.

Public API
  GLNRep.of(lam) / GLNRep(lam)  : .N .lam .pats .dim .index .wt .wsp .top .E (lazy: non-adjacent generators are built on first use) .g .check_metric()
                                  .apply(ij, vec) .sqrt_metric() .float_ops() (lazy dict of CSR, incl. diagonal) .block_exact(ij, target_w) .dense(ij)
  weyl_dim(lam)  weight_mults(lam)  racah_speiser(lam, lamp)  candidates(lam, lamp)  casimir2(lam)  nW_of(lam, lamp, mu)  dij(N, i, j)  wadd(w, d, c)
  munich_patterns(lam)  munich_lowering_sq(M, k, l)  munich_check(lam)  munich_label_to_ours(lam, lamp, mu_munich)   (bridge to Paul/Cpp/AlexCGC: same GT basis, see below)
Racah-Speiser convention used here: mult(mu) = sum over the weights w of the SECOND factor V_lamp of  m_lamp(w) sgn(s)  where  lam + w + rho  has distinct entries,
s the permutation sorting it decreasingly, mu = s(lam + w + rho) - rho  (so the dominant candidates are mu = lam + w, w in wt(V_lamp)).
"""
import os, sys, math, itertools
from fractions import Fraction as Fr
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, '..')))
import numpy as np
import scipy.sparse as sps
from gl3_reps import gt_patterns, _raise, _lower, _weight   # read-only reuse of the (general-N) GT formulas


# ======================================================================================== small helpers
def dij(N, i, j):
    """e_i - e_j as an N-tuple (1-based i, j): the weight shift of E_ij"""
    return tuple((1 if k == i - 1 else 0) - (1 if k == j - 1 else 0) for k in range(N))


def wadd(w, d, c=1):
    return tuple(x + c * y for x, y in zip(w, d))


def perm_sign(p):
    """sign of the permutation given as a list of distinct comparable entries (parity of the number of inversions)"""
    s = 1
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


def _frac(c):
    return Fr(int(c.p), int(c.q))


def _t(x):
    return tuple(int(v) for v in x)


def _compose(A, B):
    """(A o B)[c] = A applied to (B e_c);  A, B = lists of dicts row -> Fraction (columns)"""
    out = []
    for col in B:
        acc = {}
        for r, v in col.items():
            for r2, v2 in A[r].items():
                acc[r2] = acc.get(r2, 0) + v * v2
        out.append({r: v for r, v in acc.items() if v != 0})
    return out


def _sub(A, B):
    out = []
    for ca, cb in zip(A, B):
        acc = dict(ca)
        for r, v in cb.items():
            acc[r] = acc.get(r, 0) - v
        out.append({r: v for r, v in acc.items() if v != 0})
    return out


def _check_lam(lam):
    lam = tuple(int(x) for x in lam)
    if len(lam) < 2 or any(lam[i] < lam[i + 1] for i in range(len(lam) - 1)) or lam[-1] < 0:
        raise ValueError('lam must be a weakly decreasing tuple of non-negative integers, got %s' % (lam,))
    return lam


# ======================================================================================== module-level representation theory
def weyl_dim(lam):
    """Weyl dimension formula for gl(N): prod_{i<j} (lam_i - lam_j + j - i)/(j - i)"""
    lam = tuple(lam)
    N = len(lam)
    num = den = 1
    for i in range(N):
        for j in range(i + 1, N):
            num *= lam[i] - lam[j] + j - i
            den *= j - i
    assert num % den == 0
    return num // den


_WM_CACHE = {}


def weight_mults(lam):
    """Counter: weight (N-tuple) -> multiplicity in V_lam (counted on the GT patterns; cached)"""
    lam = tuple(int(x) for x in lam)
    c = _WM_CACHE.get(lam)
    if c is None:
        N = len(lam)
        c = _WM_CACHE[lam] = Counter(tuple(_weight(k, p) for k in range(1, N + 1)) for p in gt_patterns(lam))
    return c


def racah_speiser(lam, lamp=None):
    """{mu: multiplicity} of V_lam (x) V_lamp (default lamp = lam), mu dominant.  Sum over the weights of V_lamp, shifted lam (see module docstring).
    Raises on a negative multiplicity; Sum m dim V_mu = dim V_lam dim V_lamp is checked by the tests."""
    lam = tuple(int(x) for x in lam)
    lamp = lam if lamp is None else tuple(int(x) for x in lamp)
    N = len(lam)
    assert len(lamp) == N
    rho = tuple(range(N - 1, -1, -1))
    acc = defaultdict(int)
    for w, m in weight_mults(lamp).items():
        x = tuple(lam[i] + w[i] + rho[i] for i in range(N))
        if len(set(x)) < N:
            continue
        order = sorted(range(N), key=lambda i: -x[i])
        sgn = perm_sign(order)
        mu = tuple(x[order[i]] - rho[i] for i in range(N))
        acc[mu] += sgn * m
    out = {}
    for mu, m in acc.items():
        if m < 0:
            raise AssertionError(('negative multiplicity', mu, m))
        if m > 0:
            out[mu] = m
    return out


def candidates(lam, lamp=None, both=False):
    """all dominant mu = lam + w, w a weight of V_lamp (a superset of the irreps of V_lam (x) V_lamp, includes mu of multiplicity 0), sorted descending.
    both=True also adds the dominant mu = lamp + w1, w1 a weight of V_lam (the mu for which the weight space V_lam[mu - lamp] of the slice is non-empty)."""
    lam = tuple(lam)
    lamp = lam if lamp is None else tuple(lamp)
    out = set()
    for w in weight_mults(lamp):
        mu = tuple(a + b for a, b in zip(lam, w))
        if all(mu[i] >= mu[i + 1] for i in range(len(mu) - 1)):
            out.add(mu)
    if both:
        for w in weight_mults(lam):
            mu = tuple(a + b for a, b in zip(lamp, w))
            if all(mu[i] >= mu[i + 1] for i in range(len(mu) - 1)):
                out.add(mu)
    return sorted(out, reverse=True)


def casimir2(lam):
    """scalar value of C_2 = Sum_ij E_ij E_ji on V_lam:  Sum_i lam_i (lam_i + N + 1 - 2 i)"""
    N = len(lam)
    return sum(lam[i - 1] * (lam[i - 1] + N + 1 - 2 * i) for i in range(1, N + 1))


def nW_of(lam, lamp, mu):
    """dimension of the weight-mu subspace of V_lam (x) V_lamp:  Sum_w m_lamp(w) m_lam(mu - w)"""
    wl, wp = weight_mults(lam), weight_mults(lamp)
    return sum(m * wl.get(tuple(a - b for a, b in zip(mu, w)), 0) for w, m in wp.items())


# ======================================================================================== the representation
class _Gens(dict):
    """E[(i,j)] for 1 <= i,j <= N: diagonal and adjacent matrices are built with the representation, every other one lazily by commutators"""

    def __init__(self, rep):
        super().__init__()
        self.rep = rep

    def __missing__(self, key):
        i, j = key
        if j > i + 1:
            A, B = self[(i, j - 1)], self[(j - 1, j)]      # E_ij = [E_{i,j-1}, E_{j-1,j}]
        elif i > j + 1:
            A, B = self[(i, i - 1)], self[(i - 1, j)]      # E_ij (lower) = [E_{i,i-1}, E_{i-1,j}]
        else:
            raise KeyError(key)
        val = _sub(_compose(A, B), _compose(B, A))
        self[key] = val
        return val


class _FloatOps(dict):
    """lazy dict ij -> scipy CSR matrix of E_ij in the ORTHONORMAL GT basis (Eh[r,q] = E[r,q] sqrt(g_r/g_q); Eh_ji = Eh_ij^T), incl. the diagonal"""

    def __init__(self, rep):
        super().__init__()
        self.rep = rep

    def __missing__(self, key):
        R = self.rep
        n = R.dim
        i, j = key
        if i == j:
            val = sps.diags([float(w[i - 1]) for w in R.wt]).tocsr()
        else:
            g = R.g
            cols = R.E[key]
            rows_, cols_, vals_ = [], [], []
            for q in range(n):
                for r, v in cols[q].items():
                    rows_.append(r)
                    cols_.append(q)
                    vals_.append(float(v) * math.sqrt(float(g[r] / g[q])))
            val = sps.csr_matrix((vals_, (rows_, cols_)), shape=(n, n))
        self[key] = val
        return val


class GLNRep:
    """exact sparse GT irrep lam of gl(N), N = len(lam).  See the module docstring for the conventions."""
    _cache = {}

    @classmethod
    def of(cls, lam):
        lam = tuple(int(x) for x in lam)
        if lam not in cls._cache:
            cls._cache[lam] = cls(lam)
        return cls._cache[lam]

    def __init__(self, lam):
        lam = _check_lam(lam)
        N = self.N = len(lam)
        self.lam = lam
        self.pats = gt_patterns(lam)
        self.dim = n = len(self.pats)
        self.index = {p: i for i, p in enumerate(self.pats)}
        self.wt = [tuple(_weight(k, p) for k in range(1, N + 1)) for p in self.pats]
        d = defaultdict(list)
        for i, w in enumerate(self.wt):
            d[w].append(i)
        self.wsp = {w: np.array(v, dtype=np.int64) for w, v in d.items()}
        self.top = self.index[tuple(tuple(lam[:k]) for k in range(1, N + 1))]
        assert self.top == n - 1, 'the highest-weight pattern must be the LAST index'
        E = _Gens(self)
        for k in range(1, N + 1):
            E[(k, k)] = [({c: Fr(w[k - 1])} if w[k - 1] != 0 else {}) for c, w in enumerate(self.wt)]
        for k in range(1, N):
            U, D = [], []
            for p in self.pats:
                u, dd = {}, {}
                for tg, c in _raise(k, p).items():
                    j = self.index.get(tg)          # the raw GT formula can propose a target outside the interlacing range: skip it (never clip)
                    if j is not None:
                        u[j] = _frac(c)
                for tg, c in _lower(k, p).items():
                    j = self.index.get(tg)
                    if j is not None:
                        dd[j] = _frac(c)
                U.append(u)
                D.append(dd)
            E[(k, k + 1)] = U
            E[(k + 1, k)] = D
        self.E = E
        self._g = None
        self._Eh = None

    # ---- invariant metric: diagonal G with E_ij^T G = G E_ji, G[top] = 1, exact
    @property
    def g(self):
        if self._g is None:
            n = self.dim
            g = [None] * n
            g[self.top] = Fr(1)
            rev = defaultdict(list)        # r -> [(q, E_up[r,q], E_dn[q,r])] for the simple raising operators  (g_r E_up[r,q] = g_q E_dn[q,r])
            for k in range(1, self.N):
                U, D = self.E[(k, k + 1)], self.E[(k + 1, k)]
                for q in range(n):
                    for r, v in U[q].items():
                        rev[r].append((q, v, D[r].get(q, 0)))
            stack = [self.top]
            while stack:
                r = stack.pop()
                for q, v, w in rev[r]:
                    if g[q] is None and w != 0:
                        g[q] = g[r] * v / w
                        stack.append(q)
            assert all(x is not None and x > 0 for x in g)
            self._g = g
        return self._g

    def check_metric(self, pairs=None):
        """exact: g_r E_ij[r,q] = g_q E_ji[q,r] for ALL matrix elements of E_ij and E_ji, all i < j (or the given list of pairs); returns the number of violations"""
        g = self.g
        bad = 0
        pairs = pairs or [(i, j) for i in range(1, self.N + 1) for j in range(i + 1, self.N + 1)]
        for (i, j) in pairs:
            U, D = self.E[(i, j)], self.E[(j, i)]
            for q in range(self.dim):
                for r, v in U[q].items():              # g_r E_ij[r,q] = g_q E_ji[q,r]
                    if g[r] * v != g[q] * D[r].get(q, 0):
                        bad += 1
                for r, w in D[q].items():              # entries of E_ji without partner in E_ij: g_q E_ij[q,r] = g_r E_ji[r,q]
                    if g[q] * U[r].get(q, 0) != g[r] * w:
                        bad += 1
        return bad

    # ---- exact vector operations on dicts {index: Fraction}
    def apply(self, ij, vec):
        out = {}
        cols = self.E[ij]
        for c, v in vec.items():
            for r, e in cols[c].items():
                out[r] = out.get(r, 0) + v * e
        return {r: v for r, v in out.items() if v != 0}

    def sqrt_metric(self):
        """float sqrt(g) per state (orthonormal GT basis e~_q = e_q / sqrt(g_q))"""
        if not hasattr(self, '_sg'):
            self._sg = np.array([math.sqrt(float(x)) for x in self.g])
        return self._sg

    def float_ops(self):
        """lazy dict ij -> scipy.sparse CSR of E_ij (incl. the diagonal) in the ORTHONORMAL GT basis"""
        if self._Eh is None:
            self._Eh = _FloatOps(self)
        return self._Eh

    def block_exact(self, ij, target_w):
        """exact dense block of E_ij: rows wsp[target_w], cols wsp[target_w - d_ij]  (list of lists of Fractions), plus (rows, cols) index arrays"""
        rows = self.wsp.get(tuple(target_w))
        cols = self.wsp.get(wadd(target_w, dij(self.N, *ij), -1))
        if rows is None or cols is None:
            return [], rows, cols
        pos = {int(r): a for a, r in enumerate(rows)}
        M = [[Fr(0)] * len(cols) for _ in rows]
        E = self.E[ij]
        for b, c in enumerate(cols):
            for r, v in E[int(c)].items():
                M[pos[r]][b] = v
        return M, rows, cols

    def dense(self, ij):
        """numpy float64 dense matrix of E_ij in the PLAIN GT basis (tests / independent re-verification only)"""
        n = self.dim
        M = np.zeros((n, n))
        for c, col in enumerate(self.E[ij]):
            for r, v in col.items():
                M[r, c] = float(v)
        return M


# ======================================================================================== bridge to Munich's code (Paul/Cpp/AlexCGC/ClebschGordan.cpp)
# The Munich binary cannot be built on the Windows machine; what can be settled offline is that the two codes use the SAME GT basis:
#  (i) Munich's pattern index Q(M) (clebsch::pattern::operator++, Eq. C7/C9) = our GT index + 1;  (ii) Munich's orthonormal GT matrix elements
#  (clebsch::pattern::lowering_coeff, Eq. 28) squared equal ours exactly (E_{l+1,l}[q',q]^2 g_q'/g_q), and our plain-basis simple raising/lowering coefficients all have the
#  sign of Munich's (>= 0), i.e. the diagonal sign gauge eps_q relating the two bases is trivial.  See munich_check().

def munich_patterns(lam):
    """re-implementation of clebsch::pattern(irrep, index) iterated with operator++ (ClebschGordan.cpp 341-395): the patterns in Munich's Q order, as our row tuples
    (Munich's M(k,l) = our rows[l-1][k-1])"""
    lam = _t(lam)
    N = len(lam)
    M = {}
    for i in range(1, N + 1):
        M[(i, N)] = lam[i - 1]
    for l in range(N - 1, 0, -1):
        for k in range(1, l + 1):
            M[(k, l)] = M[(k + 1, l + 1)]

    def rows():
        return tuple(tuple(M[(k, l)] for k in range(1, l + 1)) for l in range(1, N + 1))
    out = [rows()]
    while True:
        k, l = 1, 1
        while l < N and M[(k, l)] == M[(k, l + 1)]:
            k -= 1
            if k == 0:
                l += 1
                k = l
        if l == N:
            break
        M[(k, l)] += 1
        while k != 1 or l != 1:
            k += 1
            if k > l:
                k = 1
                l -= 1
            M[(k, l)] = M[(k + 1, l + 1)]
        out.append(rows())
    return out


def munich_lowering_sq(M, k, l):
    """(clebsch::pattern::lowering_coeff(k, l))^2 = -result, exact; M = our rows tuple (the state being lowered)"""
    def m(i, ll):
        return M[ll - 1][i - 1]
    res = Fr(1)
    for i in range(1, l + 2):
        res *= m(i, l + 1) - m(k, l) + k - i + 1
    for i in range(1, l):
        res *= m(i, l - 1) - m(k, l) + k - i
    for i in range(1, l + 1):
        if i == k:
            continue
        res /= (m(i, l) - m(k, l) + k - i + 1)
        res /= (m(i, l) - m(k, l) + k - i)
    return -res


def munich_check(lam):
    """(i) Munich's pattern order = gl3_reps.gt_patterns order; (ii) for EVERY state q, level l and entry k whose lowered pattern exists: (our orthonormal lowering matrix
    element)^2 = E_{l+1,l}[q',q]^2 g_q'/g_q  equals Munich's lowering_coeff^2 exactly, and E_{l+1,l} has no other nonzero entries; (iii) the diagonal sign gauge eps with
    eps_q' eps_q sign(E_{l+1,l}[q',q]) = +1 (Munich's coefficients are >= 0) exists (BFS from the HW state) and is consistent on all matrix elements.
    Returns dict(order_ok, nelem, bad_sq, extra, gauge_bad, eps)."""
    lam = _t(lam)
    R = GLNRep.of(lam)
    N = R.N
    g = R.g
    order_ok = (munich_patterns(lam) == R.pats)
    nelem = bad_sq = extra = 0
    edges = defaultdict(list)                         # q -> [(q', sign)]
    for q, M in enumerate(R.pats):
        for l in range(1, N):
            col = R.E[(l + 1, l)][q]
            matched = set()
            for k in range(1, l + 1):
                rows = list(M)
                row = list(rows[l - 1])
                row[k - 1] -= 1
                rows[l - 1] = tuple(row)
                qp = R.index.get(tuple(rows))
                if qp is None:
                    continue
                e = col.get(qp, 0)
                nelem += 1
                if e * e * g[qp] / g[q] != munich_lowering_sq(M, k, l):
                    bad_sq += 1
                if e != 0:
                    matched.add(qp)
                    edges[q].append((qp, 1 if e > 0 else -1))
            extra += len(set(col) - matched)
    eps = [0] * R.dim
    eps[R.top] = 1
    stack = [R.top]
    while stack:
        q = stack.pop()
        for qp, s in edges[q]:
            if eps[qp] == 0:
                eps[qp] = eps[q] * s
                stack.append(qp)
    gauge_bad = sum(1 for q in range(R.dim) if eps[q] == 0)
    for q in range(R.dim):
        for qp, s in edges[q]:
            if eps[qp] != 0 and eps[q] != 0 and eps[qp] * eps[q] * s != 1:
                gauge_bad += 1
    return dict(order_ok=order_ok, nelem=nelem, bad_sq=bad_sq, extra=extra, gauge_bad=gauge_bad, eps=eps, trivial_gauge=all(e == 1 for e in eps))


def munich_label_to_ours(lam, lamp, mum):
    """Munich's SU(N) label (mu_N normalised to 0) -> our gl(N) label (|mu| = |lam| + |lamp|): add c (1,...,1) with c = (|lam| + |lamp| - |mu_M|)/N"""
    N = len(lam)
    s = sum(lam) + sum(lamp) - sum(mum)
    assert s % N == 0, (lam, lamp, mum)
    return tuple(x + s // N for x in mum)
