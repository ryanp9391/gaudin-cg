"""slice_engine.py -- round 6, WP-R: the PRV/Zhelobenko slice engine for gl(3) highest-weight (HW) vectors of V_Lam (x) V_Lam.

Maths (SPEC-slice-engine.md, plan section 2.1):  write Psi0 = sum_q psi_q (x) e_q (q = GT states of site 2).  The slice psi_top = psi_{e_top}
(site 2 = highest weight state) determines Psi0 completely (level recursion over the weights of site 2), and the multiplicity space of
mu is  M_mu = { v in V_Lam[mu - Lam] : E12^(l1-l2+1) v = E23^(l2-l3+1) v = 0 }  (Zhelobenko/PRV).  The Gaudin vertex operator V acts on
M_mu as an m x m matrix V_red that needs the slice only (not n_W, not even the level-1 psi_q: the raising generators annihilate Psi0).

Public API (all take lam = (l1,l2,l3), mu = (m1,m2,m3) as tuples of ints; everything is exact unless it says float)
  Rep.of(lam)                         exact sparse GT representation (gl3_reps order), weight spaces, invariant metric g
  racah_speiser(lam)                  {mu: multiplicity} of V_lam (x) V_lam
  prv_space(lam, mu) -> PRV           exact basis of M_mu (.basis[k] aligned with .W = global indices of V_lam[mu-lam]; .m = multiplicity)
  V_red(lam, mu, prv) -> m x m        exact Fractions, column j = V_red b_j in PRV coordinates
  copies(prv, Vm, dps) -> list        mpmath eigen-decomposition: ascending eigenvalue, slice with first nonzero component 1
  hw_from_slice(lam, mu, v) -> HW     float level recursion (orthonormal GT basis), v of shape (w,) or (w,k); COO vectors in the plain GT basis
  hw_from_slice_exact(lam, mu, v)     exact twin (dict {(nu,q): Fraction}) + number of exact consistency conditions verified
  hw_from_slice_rounded(lam, mu, v)   exact recursion rounded to double (slow, HW residual ~1e-15 at any a; the float64 recursion degrades like 2.2^a)
  reference_irrep(lam, mu)            all of the above for one irrep
  hw_residual / check_V / exact_hw_check / validate_w4 / bench / write_ref / verify_ref / bruteforce_check / tridiag_check   (tests, benchmark, reference files)
CLI (python -X utf8 slice_engine.py <cmd>): table1 | hw | neg | w4 | bench [amax [amin]] | ref [l1 l2 | all] | verify | brute | tridiag | closed | rounded | exact <a> | eigcheck | all
(the full `all` takes ~8 min, dominated by `verify` (3.5 min) and `brute` (1.5 min); `exact 14` (11 min) is not part of `all`)

Layout of this file
  1. Rep            exact (Fraction) sparse GT representation in the gl3_reps order, weight spaces, invariant metric
  2. Racah-Speiser  multiplicities from GT weight multiplicities, dominant candidates
  3. prv_space      exact basis of M_mu  (nullspace over Q of the two chain constraints on the weight space V[mu-Lam])
  4. V_red          exact m x m matrix of the vertex operator on M_mu (slice only);  copies(): eigen-decomposition (mpmath)
  5. Recursion      float level recursion in the orthonormal GT basis (numpy, left inverses of the site-2 raising blocks);
                    hw_from_slice();  exact twin ExactRecursion / hw_from_slice_exact() (Fractions)
  6. checks         hw_residual(), apply_V_full() (the 36-term V on a full vector)
  7. W4             validation against ref/untwisted-exact-*.txt
  8. tables/tests/benchmark/structure checks   9. reference files   10. independent re-verification of the written files, brute-force null space   CLI at the bottom

Conventions (binding, copied from gaudin_chain.py / gl3_reps.py): GT states in the order of GL3Rep(lam).patterns (ascending lexicographic in
(m1,m2,k); the highest-weight state is the LAST index); matrices act on kets as columns, E e_col = sum_row M[row,col] e_row;
E13 = [E12,E23], E31 = [E32,E21]; weights are (E11,E22,E33) eigenvalues; two-site flat index = nu*dim + q, nu = site 1 (slow), q = site 2 (fast).
"""
import os, sys, time, itertools, math
from fractions import Fraction as Fr
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla
import mpmath as mp
from gl3_reps import gt_patterns, _raise, _lower, _weight   # read-only reuse of the GT formulas

A1 = (1, -1, 0)      # simple roots in (E11,E22,E33) weight coordinates
A2 = (0, 1, -1)
DIJ = {(i, j): tuple((1 if k == i - 1 else 0) - (1 if k == j - 1 else 0) for k in range(3)) for i in (1, 2, 3) for j in (1, 2, 3)}


def wadd(w, d, c=1):
    return (w[0] + c * d[0], w[1] + c * d[1], w[2] + c * d[2])


def perm_sign(p):
    s = 1
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


def dim_gl3(mu):
    return (mu[0] - mu[1] + 1) * (mu[1] - mu[2] + 1) * (mu[0] - mu[2] + 2) // 2


def _frac(c):
    return Fr(int(c.p), int(c.q))


# ======================================================================================== 1. exact GT representation
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


class Rep:
    """exact sparse GT irrep lam in the gl3_reps order.  E[(i,j)] (i != j) = list over columns of dicts {row: Fraction}."""
    _cache = {}

    @classmethod
    def of(cls, lam):
        lam = tuple(int(x) for x in lam)
        if lam not in cls._cache:
            cls._cache[lam] = cls(lam)
        return cls._cache[lam]

    def __init__(self, lam):
        lam = tuple(int(x) for x in lam)
        if len(lam) != 3 or not (lam[0] >= lam[1] >= lam[2] >= 0):
            raise ValueError(lam)
        self.lam = lam
        self.pats = gt_patterns(lam)
        self.dim = n = len(self.pats)
        self.index = {p: i for i, p in enumerate(self.pats)}
        self.wt = [tuple(_weight(k, p) for k in (1, 2, 3)) for p in self.pats]
        d = defaultdict(list)
        for i, w in enumerate(self.wt):
            d[w].append(i)
        self.wsp = {w: np.array(v, dtype=np.int64) for w, v in d.items()}
        self.top = self.index[tuple(tuple(lam[:k]) for k in (1, 2, 3))]
        E = {}
        for k in (1, 2):
            U, D = [], []
            for p in self.pats:
                u, dd = {}, {}
                for tg, c in _raise(k, p).items():
                    j = self.index.get(tg)          # raw GT formula can propose a target outside the interlacing range: skip (as GL3Rep does)
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
        E[(1, 3)] = _sub(_compose(E[(1, 2)], E[(2, 3)]), _compose(E[(2, 3)], E[(1, 2)]))      # [E12,E23]
        E[(3, 1)] = _sub(_compose(E[(3, 2)], E[(2, 1)]), _compose(E[(2, 1)], E[(3, 2)]))      # [E32,E21]
        self.E = E
        self._g = None
        self._blk = {}

    # ---- metric: diagonal G with E_ij^T G = G E_ji, G[top] = 1 (W4's GMET1); exact, verified on every matrix element
    @property
    def g(self):
        if self._g is None:
            n = self.dim
            g = [None] * n
            g[self.top] = Fr(1)
            rev = defaultdict(list)                      # r -> [(q, E_up[r,q], E_dn[q,r])]  for the simple raising operators
            for k in (1, 2):
                U, D = self.E[(k, k + 1)], self.E[(k + 1, k)]
                for q in range(n):
                    for r, v in U[q].items():
                        rev[r].append((q, v, D[r].get(q, 0)))
            stack = [self.top]
            while stack:
                r = stack.pop()
                for q, v, w in rev[r]:
                    if g[q] is None and w != 0:
                        g[q] = g[r] * v / w              # g_r E_up[r,q] = g_q E_dn[q,r]
                        stack.append(q)
            assert all(x is not None and x > 0 for x in g)
            self._g = g
        return self._g

    def check_metric(self):
        """exact: g_r E_ij[r,q] = g_q E_ji[q,r] for ALL six off-diagonal matrices; returns number of violations"""
        g = self.g
        bad = 0
        for (i, j) in ((1, 2), (2, 3), (1, 3)):
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
        i, j = ij
        out = {}
        if i == j:
            k = i - 1
            for c, v in vec.items():
                w = self.wt[c][k]
                if w != 0:
                    out[c] = v * w
            return out
        cols = self.E[ij]
        for c, v in vec.items():
            for r, e in cols[c].items():
                out[r] = out.get(r, 0) + v * e
        return {r: v for r, v in out.items() if v != 0}

    def sqrt_metric(self):
        """float sqrt(g) per state (orthonormal GT basis e~_q = e_q / sqrt(g_q)); overflow-safe via ratios elsewhere"""
        if not hasattr(self, '_sg'):
            self._sg = np.array([math.sqrt(float(x)) for x in self.g])
        return self._sg

    def float_ops(self):
        """scipy.sparse CSR matrices of all E_ij (incl. diagonal) in the ORTHONORMAL basis: Eh[r,q] = E[r,q] sqrt(g_r/g_q); Eh_ji = Eh_ij^T"""
        if not hasattr(self, '_Eh'):
            n = self.dim
            g = self.g
            Eh = {}
            for ij, cols in self.E.items():
                rows_, cols_, vals_ = [], [], []
                for q in range(n):
                    for r, v in cols[q].items():
                        rows_.append(r)
                        cols_.append(q)
                        vals_.append(float(v) * math.sqrt(float(g[r] / g[q])))
                Eh[ij] = sps.csr_matrix((vals_, (rows_, cols_)), shape=(n, n))
            for k in (1, 2, 3):
                Eh[(k, k)] = sps.diags([float(w[k - 1]) for w in self.wt]).tocsr()
            self._Eh = Eh
        return self._Eh

    def block_exact(self, ij, target_w):
        """exact dense block of E_ij: rows wsp[target_w], cols wsp[target_w - d_ij]  (list of lists of Fractions)"""
        rows = self.wsp.get(target_w)
        cols = self.wsp.get(wadd(target_w, DIJ[ij], -1))
        if rows is None or cols is None:
            return [], rows, cols
        pos = {int(r): a for a, r in enumerate(rows)}
        M = [[Fr(0)] * len(cols) for _ in rows]
        E = self.E[ij]
        for b, c in enumerate(cols):
            for r, v in E[int(c)].items():
                M[pos[r]][b] = v
        return M, rows, cols


# ======================================================================================== 2. Racah-Speiser
def weight_mults(lam):
    """Counter: weight -> multiplicity in V_lam (from the GT patterns)"""
    return Counter(tuple(_weight(k, p) for k in (1, 2, 3)) for p in gt_patterns(tuple(lam)))


def racah_speiser(lam, lam2=None):
    """decomposition of V_lam (x) V_lam2 (default lam2 = lam): {mu: multiplicity} by the Racah-Speiser/Klimyk rule from the GT weight
    multiplicities of V_lam;  also returns nothing else.  Checked: sum m dim V_mu = dim V_lam dim V_lam2."""
    lam = tuple(lam)
    lam2 = lam if lam2 is None else tuple(lam2)
    wm = weight_mults(lam)
    rho = (2, 1, 0)
    acc = defaultdict(int)
    for w, m in wm.items():
        x = tuple(lam2[i] + w[i] + rho[i] for i in range(3))
        if len(set(x)) < 3:
            continue
        order = sorted(range(3), key=lambda i: -x[i])
        sgn = perm_sign(order)
        xs = [x[i] for i in order]
        mu = tuple(xs[i] - rho[i] for i in range(3))
        acc[mu] += sgn * m
    out = {}
    for mu, m in acc.items():
        if m < 0:
            raise AssertionError(('negative multiplicity', mu, m))
        if m > 0:
            out[mu] = m
    return out


def nW_of(lam, mu, wm=None):
    """dimension of the weight-mu subspace of V_lam (x) V_lam"""
    wm = wm or weight_mults(lam)
    return sum(m * wm.get((mu[0] - w[0], mu[1] - w[1], mu[2] - w[2]), 0) for w, m in wm.items())


def candidates(lam):
    """all dominant mu = Lam + w, w a weight of V_lam (a superset of the irreps of V_lam (x) V_lam; includes mu with multiplicity 0)"""
    lam = tuple(lam)
    out = []
    for w in weight_mults(lam):
        mu = (lam[0] + w[0], lam[1] + w[1], lam[2] + w[2])
        if mu[0] >= mu[1] >= mu[2]:
            out.append(mu)
    return sorted(set(out), reverse=True)


# ======================================================================================== 3. PRV space
def _rref_nullspace(rows, ncols):
    """exact nullspace of the matrix 'rows' (list of lists of Fractions, ncols columns): basis from the RREF, one vector per free column f
    (entry 1 at f, zero at the other free columns).  Returns (basis, free) with basis[k] a list of ncols Fractions."""
    M = [list(r) for r in rows]
    piv = []
    r = 0
    for c in range(ncols):
        p = None
        for i in range(r, len(M)):
            if M[i][c] != 0:
                p = i
                break
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        inv = 1 / M[r][c]
        M[r] = [x * inv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        piv.append(c)
        r += 1
        if r == len(M):
            break
    free = [c for c in range(ncols) if c not in set(piv)]
    basis = []
    for f in free:
        v = [Fr(0)] * ncols
        v[f] = Fr(1)
        for k, pc in enumerate(piv):
            v[pc] = -M[k][f]
        basis.append(v)
    return basis, free


class PRV:
    """exact basis of M_mu inside the weight space W = V_lam[mu-lam]:  basis[k] = list of Fractions aligned with W (global GT indices), the
    k-th vector is 1 at the k-th free position and 0 at the other free positions (so the coordinates of any v in M_mu are its free entries)"""
    def __init__(self, lam, mu, W, basis, free):
        self.lam, self.mu, self.W, self.basis, self.free = lam, mu, W, basis, free
        self.m = len(basis)
        self.w = len(W)

    def vec(self, k):
        return {int(c): x for c, x in zip(self.W, self.basis[k]) if x != 0}

    def coords(self, v):
        """coordinates of a vector (list aligned with W) in the basis; raises if it is not in M_mu"""
        co = [v[f] for f in self.free]
        rec = [Fr(0)] * self.w
        for k, ck in enumerate(co):
            if ck != 0:
                for i, x in enumerate(self.basis[k]):
                    rec[i] += ck * x
        if rec != list(v):
            raise ValueError('vector not in M_mu')
        return co


def power_images(rep, ij, k, W):
    """exact images E_ij^k e_c for c in W (global indices): list of dicts"""
    out = []
    for c in W:
        vec = {int(c): Fr(1)}
        for _ in range(k):
            vec = rep.apply(ij, vec)
            if not vec:
                break
        out.append(vec)
    return out


def prv_space(lam, mu):
    """exact basis of M_mu = {v in V_lam[mu-lam] : E12^(l1-l2+1) v = 0, E23^(l2-l3+1) v = 0}  (lam = Lam, mu = a dominant weight)"""
    lam = tuple(int(x) for x in lam)
    mu = tuple(int(x) for x in mu)
    rep = Rep.of(lam)
    w0 = (mu[0] - lam[0], mu[1] - lam[1], mu[2] - lam[2])
    W = rep.wsp.get(w0)
    if W is None:
        return PRV(lam, mu, np.zeros(0, dtype=np.int64), [], [])
    k1, k2 = lam[0] - lam[1] + 1, lam[1] - lam[2] + 1
    rows = []
    for ij, k in (((1, 2), k1), ((2, 3), k2)):
        if wadd(w0, DIJ[ij], k) not in rep.wsp:
            continue                                           # E_ij^k maps into a non-weight: constraint void
        imgs = power_images(rep, ij, k, W)
        byrow = defaultdict(lambda: [Fr(0)] * len(W))
        for a, vec in enumerate(imgs):
            for r, x in vec.items():
                byrow[r][a] = x
        rows.extend(byrow.values())
    basis, free = _rref_nullspace(rows, len(W))
    return PRV(lam, mu, W, basis, free)


# ======================================================================================== 4. V_red and the copies
# the 36 terms {coef, (s1,t1), (s2,t2), (s3,t3)} of V = 1/2 sum_{s,t in S3} sgn(s) sgn(t) e1_{s1 t1} e1_{s2 t2} e2_{s3 t3}  (the notebook's vTriples)
VTERMS = [(Fr(perm_sign(s) * perm_sign(t), 2), (s[0], t[0]), (s[1], t[1]), (s[2], t[2]))
          for s in itertools.permutations((1, 2, 3)) for t in itertools.permutations((1, 2, 3))]


def _vadd(acc, vec, c=1):
    for r, x in vec.items():
        acc[r] = acc.get(r, 0) + c * x


def vred_apply(rep, v):
    """V_red on a slice v (dict {global GT index: Fraction}, weight mu - Lam; v in M_mu):
         V_red = 1/2 sum_{s,t} sgn(s) sgn(t) E_{s1 t1} E_{s2 t2} Phi_{s3 t3},   Phi_cc = lam_c,  Phi_ab = -E_ab (a < b),  Phi_ab = 0 (a > b),
       because sum_q (e^2_ab)_{top,q} psi_q = -E_ab psi_top for a < b (the raising generators annihilate Psi0 and project on q = top),
       = lam_c psi_top for a = b, and 0 for a > b (a lowering operator never reaches the top state).  No site-2 information is needed."""
    lam = rep.lam
    phi = {}
    out = {}
    for coef, ij1, ij2, ij3 in VTERMS:
        a, b = ij3
        if a > b:
            continue
        if ij3 not in phi:
            if a == b:
                phi[ij3] = {c: x * lam[a - 1] for c, x in v.items() if lam[a - 1] != 0}
            else:
                e = rep.apply(ij3, v)
                phi[ij3] = {c: -x for c, x in e.items()}
        ph = phi[ij3]
        if not ph:
            continue
        _vadd(out, rep.apply(ij1, rep.apply(ij2, ph)), coef)
    return {r: x for r, x in out.items() if x != 0}


def V_red(lam, mu, prv=None):
    """exact m x m matrix of V_red in the PRV basis (column j = V_red b_j in PRV coordinates);  list of lists of Fractions.
       Raises if some V_red b_j leaves M_mu (it must not)."""
    prv = prv or prv_space(lam, mu)
    rep = Rep.of(lam)
    m = prv.m
    if m == 0:
        return []
    pos = {int(c): a for a, c in enumerate(prv.W)}
    cols = []
    for k in range(m):
        img = vred_apply(rep, prv.vec(k))
        vec = [Fr(0)] * prv.w
        for r, x in img.items():
            vec[pos[r]] = x
        cols.append(prv.coords(vec))
    return [[cols[j][i] for j in range(m)] for i in range(m)]


def mpf_of(x, dps=None):
    return mp.mpf(x.numerator) / x.denominator


def copies(prv, Vm, dps=50):
    """eigen-decomposition of the exact m x m matrix Vm (mpmath, dps digits): list of dicts, ascending eigenvalue:
       {'ev': mpf, 'coords': [mpf]*m, 'slice': {global idx: mpf}} with the slice normalised to first nonzero component (ascending GT index) = 1.
       Raises if an eigenvalue is not real to 10^-(dps/2)  (it is: V is self-adjoint for the invariant form)."""
    old = mp.mp.dps
    mp.mp.dps = dps
    m = prv.m
    try:
        if m == 1:
            evs = [mpf_of(Vm[0][0])]
            vecs = [[mp.mpf(1)]]
        else:
            A = mp.matrix(m, m)
            for i in range(m):
                for j in range(m):
                    A[i, j] = mpf_of(Vm[i][j])
            E, ER = mp.eig(A)
            evs, vecs = [], []
            for k in range(m):
                if abs(mp.im(E[k])) > mp.mpf(10) ** (-dps // 2):
                    raise ArithmeticError('complex eigenvalue %s' % E[k])
                evs.append(mp.re(E[k]))
                col = [ER[i, k] for i in range(m)]
                vecs.append(col)
        order = sorted(range(m), key=lambda k: evs[k])
        out = []
        Bm = [[mpf_of(x) for x in row] for row in prv.basis]
        for k in order:
            col = vecs[k]
            # make the eigenvector real: divide by its largest component (complex phase)
            big = max(range(m), key=lambda i: abs(col[i]))
            col = [mp.re(c / col[big]) for c in col]
            sl = [mp.fsum(col[kk] * Bm[kk][a] for kk in range(m)) for a in range(prv.w)]
            tol = mp.mpf(10) ** (-(dps - 20)) * max(abs(x) for x in sl)          # eigenvector accuracy ~ 10^-(dps-few): smaller components are exact zeros
            first = next(a for a in range(prv.w) if abs(sl[a]) > tol)
            sc = sl[first]
            sl = [x / sc for x in sl]
            col = [x / sc for x in col]
            out.append(dict(ev=mp.re(evs[k]), coords=col, slice={int(prv.W[a]): sl[a] for a in range(prv.w) if abs(sl[a]) > tol}))
        return out
    finally:
        mp.mp.dps = old


# ======================================================================================== 5. the level recursion (float, orthonormal GT basis)
class HW:
    """full HW vector(s) of weight mu: COO data in the PLAIN GT basis.  rows = site-1 state nu, cols = site-2 state q (both global GT indices),
    vals[:, k] = k-th vector.  Flat two-site index = nu*dim + q.  vals_hat = the same in the orthonormal basis (what the recursion computed)."""
    def __init__(self, lam, mu, rows, cols, vals_hat, sg):
        self.lam, self.mu = lam, mu
        self.rows, self.cols = rows, cols
        self.vals_hat = vals_hat
        self.sg = sg
        self.dim = len(sg)
        self.m = vals_hat.shape[1]

    @property
    def vals(self):
        return self.vals_hat / (self.sg[self.rows] * self.sg[self.cols])[:, None]

    @property
    def nW(self):
        return len(self.rows)

    def flat(self):
        return self.rows * self.dim + self.cols

    def hat_matrix(self, k):
        return sps.csr_matrix((self.vals_hat[:, k], (self.rows, self.cols)), shape=(self.dim, self.dim))

    def dense(self, k=0):
        X = np.zeros((self.dim, self.dim))
        X[self.rows, self.cols] = self.vals[:, k]
        return X


class Recursion:
    """level recursion for HW vectors of V_lam (x) V_lam given the slice psi_top.

    Site-2 weights are labelled (a,b): wt(q) = Lam - a alpha1 - b alpha2.  With psi_q in V[mu - wt(q)], the HW conditions read, for the unknown
    block X_ab (rows q in Q(a,b)):   [E12(site 2): Q(a,b) -> Q(a-1,b) ; E23(site 2): Q(a,b) -> Q(a,b-1)] X_ab = -[E12(site 1) X_{a-1,b} ; E23(site 1) X_{a,b-1}]
    i.e. A_ab X_ab = RHS with A_ab = stacked site-2 raising blocks (full column rank: no vector of V_Lam except the top is killed by both),
    independent of mu.  X_ab = L_ab RHS with L_ab a left inverse of A_ab (computed once per Lam, SVD pseudo-inverse in the orthonormal basis)."""
    _cache = {}

    @classmethod
    def of(cls, lam):
        lam = tuple(lam)
        if lam not in cls._cache:
            cls._cache[lam] = cls(lam)
        return cls._cache[lam]

    def __init__(self, lam):
        t0 = time.perf_counter()
        self.rep = R = Rep.of(lam)
        self.lam = R.lam
        n = R.dim
        self.sg = R.sqrt_metric()
        Eh = R.float_ops()
        posw = np.zeros(n, dtype=np.int64)
        for w, idxs in R.wsp.items():
            posw[idxs] = np.arange(len(idxs))
        self.posw = posw
        # dense blocks of the orthonormal raising matrices keyed by ROW weight: B[ij][Wrow] : wsp[Wrow - d_ij] -> wsp[Wrow]
        B = {}
        for ij in ((1, 2), (2, 3)):
            coo = Eh[ij].tocoo()
            d = {}
            wt = R.wt
            for r, c, v in zip(coo.row, coo.col, coo.data):
                Wr = wt[r]
                blk = d.get(Wr)
                if blk is None:
                    blk = d[Wr] = np.zeros((len(R.wsp[Wr]), len(R.wsp[wadd(Wr, DIJ[ij], -1)])))
                blk[posw[r], posw[c]] = v
            B[ij] = d
        self.B = B
        lam_ = R.lam
        self.ab = {}
        for w, idxs in R.wsp.items():
            self.ab[(lam_[0] - w[0], w[2] - lam_[2])] = idxs
        self.order = sorted(self.ab)
        self.AL = {}
        self.smin = 1e300
        for (a, b) in self.order:
            if (a, b) == (0, 0):
                continue
            Q = self.ab[(a, b)]
            parts = []
            if (a - 1, b) in self.ab:
                parts.append(B[(1, 2)][self.wq(a - 1, b)])
            if (a, b - 1) in self.ab:
                parts.append(B[(2, 3)][self.wq(a, b - 1)])
            A = np.vstack(parts)
            U, s, Vt = np.linalg.svd(A, full_matrices=False)
            if len(s) < len(Q) or s[len(Q) - 1] < 1e-10 * s[0]:
                raise ArithmeticError('site-2 raising block (%d,%d) is not of full column rank (s=%s)' % (a, b, s))
            self.smin = min(self.smin, s[len(Q) - 1] / s[0])
            AL = (Vt.T / s) @ U.T                                     # |Q| x r : left inverse (pseudo-inverse)
            self.AL[(a, b)] = np.ascontiguousarray(AL.T)              # stored transposed, r x |Q|, so that  Rhs @ AL^T
        self.t_setup = time.perf_counter() - t0

    def wq(self, a, b):
        l = self.lam
        return (l[0] - a, l[1] + a - b, l[2] + b)

    def reconstruct(self, mu, vhat):
        """vhat: (w, m) orthonormal-basis components of m slices on W = V[mu-lam] (ascending global index).  Returns (rows, cols, vals_hat (nW, m))."""
        R = self.rep
        lam = self.lam
        w0 = (mu[0] - lam[0], mu[1] - lam[1], mu[2] - lam[2])
        m = vhat.shape[1]
        B12, B23 = self.B[(1, 2)], self.B[(2, 3)]
        wsp = R.wsp
        Y = {}
        out_rows, out_cols, out_vals = [], [], []
        for (a, b) in self.order:
            W = (w0[0] + a, w0[1] - a + b, w0[2] - b)
            S = wsp.get(W)
            if S is None:
                continue
            Q = self.ab[(a, b)]
            if (a, b) == (0, 0):
                y = np.ascontiguousarray(vhat.T[:, :, None])               # (m, |S|, 1)
            else:
                parts = []
                k1 = (a - 1, b)
                if k1 in self.ab:
                    y1 = Y.get(k1)
                    T1 = B12.get(W)
                    if y1 is None or T1 is None:
                        parts.append(np.zeros((m, len(S), len(self.ab[k1]))))
                    else:
                        parts.append(-(T1 @ y1))
                k2 = (a, b - 1)
                if k2 in self.ab:
                    y2 = Y.get(k2)
                    T2 = B23.get(W)
                    if y2 is None or T2 is None:
                        parts.append(np.zeros((m, len(S), len(self.ab[k2]))))
                    else:
                        parts.append(-(T2 @ y2))
                Rr = parts[0] if len(parts) == 1 else np.concatenate(parts, axis=2)
                y = Rr @ self.AL[(a, b)]                                    # (m, |S|, |Q|)
            Y[(a, b)] = y
            out_rows.append(np.repeat(S, len(Q)))
            out_cols.append(np.tile(Q, len(S)))
            out_vals.append(y.reshape(m, -1).T)
        return np.concatenate(out_rows), np.concatenate(out_cols), np.concatenate(out_vals, axis=0)


def hw_from_slice(lam, mu, v, check=None):
    """full HW vector(s) of weight mu in V_lam (x) V_lam from the slice(s) v (floats, plain GT basis components on the weight space V[mu-lam] in
       ascending global GT order; shape (w,) or (w,k)).  v must lie in M_mu (prv_space); consistency is verified by hw_residual (check=tol raises
       ArithmeticError if the relative HW residual exceeds tol; the float recursion itself cannot see an inconsistent slice, the exact twin can)."""
    lam = tuple(int(x) for x in lam)
    mu = tuple(int(x) for x in mu)
    rec = Recursion.of(lam)
    R = rec.rep
    w0 = (mu[0] - lam[0], mu[1] - lam[1], mu[2] - lam[2])
    W = R.wsp[w0]
    v = np.asarray(v, dtype=float)
    if v.ndim == 1:
        v = v[:, None]
    vhat = v * rec.sg[W][:, None]
    rows, cols, vals_hat = rec.reconstruct(mu, vhat)
    hw = HW(lam, mu, rows, cols, vals_hat, rec.sg)
    if check is not None:
        r, wd = hw_residual(hw, extra_E13=False)
        if r > check or wd != 0:
            raise ArithmeticError('HW residual %.2e (weight deviation %d) exceeds %.1e: the slice is not in M_mu' % (r, wd, check))
    return hw


# ---- exact twin of the recursion (Fractions, plain GT basis): small reps, exact HW vectors and the exact consistency conditions
def exact_left_inverse(A, nq):
    """A: list of r rows of nq Fractions with rank nq.  Returns (L, T): L (nq x r) with L A = 1 and T ((r-nq) x r) with T A = 0
       (the consistency conditions on a right-hand side)."""
    r = len(A)
    M = [list(A[i]) + [Fr(int(i == j)) for j in range(r)] for i in range(r)]
    row = 0
    for c in range(nq):
        p = next((i for i in range(row, r) if M[i][c] != 0), None)
        if p is None:
            raise ValueError('rank deficient exact block')
        M[row], M[p] = M[p], M[row]
        inv = 1 / M[row][c]
        M[row] = [x * inv for x in M[row]]
        for i in range(r):
            if i != row and M[i][c] != 0:
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[row])]
        row += 1
    return [M[i][nq:] for i in range(nq)], [M[i][nq:] for i in range(nq, r)]


def _obj(M, shape=None):
    a = np.empty((len(M), len(M[0]) if M else (shape[1] if shape else 0)), dtype=object)
    for i, row in enumerate(M):
        for j, x in enumerate(row):
            a[i, j] = x
    return a


class ExactRecursion:
    """the same level recursion in exact rational arithmetic (plain GT basis); also verifies, block by block, the exact consistency conditions
       T . RHS = 0 of the overdetermined system (these are the HW conditions beyond the ones used to solve)"""
    _cache = {}

    @classmethod
    def of(cls, lam):
        lam = tuple(lam)
        if lam not in cls._cache:
            cls._cache[lam] = cls(lam)
        return cls._cache[lam]

    def __init__(self, lam):
        self.rep = R = Rep.of(lam)
        self.lam = R.lam
        self.ab = {(self.lam[0] - w[0], w[2] - self.lam[2]): idxs for w, idxs in R.wsp.items()}
        self.order = sorted(self.ab)
        self.L = {}
        self.T = {}
        for (a, b) in self.order:
            if (a, b) == (0, 0):
                continue
            rows = []
            if (a - 1, b) in self.ab:
                rows += R.block_exact((1, 2), self.wq(a - 1, b))[0]
            if (a, b - 1) in self.ab:
                rows += R.block_exact((2, 3), self.wq(a, b - 1))[0]
            L, T = exact_left_inverse(rows, len(self.ab[(a, b)]))
            self.L[(a, b)] = _obj(L)
            self.T[(a, b)] = _obj(T, (0, len(rows))) if T else None

    def wq(self, a, b):
        l = self.lam
        return (l[0] - a, l[1] + a - b, l[2] + b)

    def reconstruct(self, mu, v):
        """v: list of Fractions aligned with W = V[mu-lam].  Returns dict {(nu, q): Fraction} (all entries of the weight-mu subspace; zeros kept)
           and the number of exact consistency conditions checked (raises if one fails)"""
        R = self.rep
        lam = self.lam
        w0 = (mu[0] - lam[0], mu[1] - lam[1], mu[2] - lam[2])
        Y = {}
        out = {}
        nchk = 0
        for (a, b) in self.order:
            W = (w0[0] + a, w0[1] - a + b, w0[2] - b)
            S = R.wsp.get(W)
            if S is None:
                continue
            Q = self.ab[(a, b)]
            if (a, b) == (0, 0):
                y = np.empty((len(S), 1), dtype=object)
                for i in range(len(S)):
                    y[i, 0] = Fr(v[i])
            else:
                parts = []
                k1 = (a - 1, b)
                if k1 in self.ab:
                    y1 = Y.get(k1)
                    M, rr, cc = R.block_exact((1, 2), W)
                    if y1 is None or cc is None or not M:
                        parts.append(np.full((len(S), len(self.ab[k1])), Fr(0), dtype=object))
                    else:
                        parts.append(-(_obj(M) @ y1))
                k2 = (a, b - 1)
                if k2 in self.ab:
                    y2 = Y.get(k2)
                    M, rr, cc = R.block_exact((2, 3), W)
                    if y2 is None or cc is None or not M:
                        parts.append(np.full((len(S), len(self.ab[k2])), Fr(0), dtype=object))
                    else:
                        parts.append(-(_obj(M) @ y2))
                Rr = parts[0] if len(parts) == 1 else np.concatenate(parts, axis=1)       # (|S|, r)
                y = Rr @ self.L[(a, b)].T                                                   # (|S|, |Q|)
                T = self.T[(a, b)]
                if T is not None:
                    cond = Rr @ T.T
                    nchk += cond.size
                    if any(x != 0 for x in cond.ravel()):
                        raise ValueError('exact HW consistency condition violated at block (%d,%d): the slice is not in M_mu' % (a, b))
            Y[(a, b)] = y
            for i, nu in enumerate(S):
                for j, q in enumerate(Q):
                    out[(int(nu), int(q))] = y[i, j]
        return out, nchk


def hw_from_slice_exact(lam, mu, v):
    """exact HW vector (dict {(nu,q): Fraction}) from an exact slice v (list of Fractions aligned with the weight space V[mu-lam]);
       also returns the number of exact consistency conditions verified."""
    return ExactRecursion.of(tuple(lam)).reconstruct(tuple(mu), v)


def hw_from_slice_rounded(lam, mu, v):
    """exact recursion followed by rounding to double: an HW object with ONE vector, for an exact slice v (list of Fractions aligned with the weight space).
       Slow (0.7-1.3 s per vector at [28,14,0]) but the HW residual stays ~1e-15 for every a, whereas the float64 recursion `hw_from_slice` loses a factor ~2.2
       per unit of a (1.6e-14 at a = 7, 5.3e-12 at a = 14)."""
    lam = tuple(int(x) for x in lam)
    mu = tuple(int(x) for x in mu)
    R = Rep.of(lam)
    psi, _ = hw_from_slice_exact(lam, mu, v)
    sg = R.sqrt_metric()
    rows = np.array([k[0] for k in psi], dtype=np.int64)
    cols = np.array([k[1] for k in psi], dtype=np.int64)
    vals = np.array([float(x) for x in psi.values()])
    return HW(lam, mu, rows, cols, (vals * sg[rows] * sg[cols])[:, None], sg)


def exact_hw_check(lam, psi):
    """apply the total E12, E23 (and E13) exactly to a two-site vector {(nu,q): Fraction}; returns the number of nonzero components left (0 = HW)"""
    R = Rep.of(lam)
    bad = 0
    for ij in ((1, 2), (2, 3), (1, 3)):
        E = R.E[ij]
        res = {}
        for (nu, q), x in psi.items():
            if x == 0:
                continue
            for r, e in E[nu].items():
                res[(r, q)] = res.get((r, q), 0) + x * e
            for r, e in E[q].items():
                res[(nu, r)] = res.get((nu, r), 0) + x * e
        bad += sum(1 for x in res.values() if x != 0)
    return bad


# ======================================================================================== 6. checks
def hw_residual(hw, extra_E13=True, per_vector=False):
    """max over vectors of ||(E_i (x) 1 + 1 (x) E_i) Psi|| / ||Psi|| (invariant norm = Frobenius norm in the orthonormal basis), i = 12, 23 (and 13),
       and the max deviation of the total weight from mu over all entries.  per_vector=True returns the list of residuals instead of their max."""
    R = Rep.of(hw.lam)
    Eh = R.float_ops()
    wt = np.array(R.wt)
    wdev = int(np.abs(wt[hw.rows] + wt[hw.cols] - np.array(hw.mu)).max())
    res = []
    ijs = [(1, 2), (2, 3)] + ([(1, 3)] if extra_E13 else [])
    for k in range(hw.m):
        X = hw.hat_matrix(k)
        nx = spla.norm(X)
        worst = 0.0
        for ij in ijs:
            Rm = Eh[ij] @ X + X @ Eh[ij].T
            worst = max(worst, (spla.norm(Rm) / nx) if Rm.nnz else 0.0)
        res.append(worst)
    return (res if per_vector else max(res)), wdev


def apply_V_full(lam, X):
    """the full 36-term V on a two-site vector X (scipy sparse (dim x dim), ORTHONORMAL basis; site 1 = rows): V X = 1/2 sum sgn sgn E1 E2 X E3^T"""
    R = Rep.of(lam)
    Eh = R.float_ops()
    cache = {}
    out = None
    for coef, ij1, ij2, ij3 in VTERMS:
        if ij3 not in cache:
            cache[ij3] = X @ Eh[ij3].T
        term = Eh[ij1] @ (Eh[ij2] @ cache[ij3])
        out = float(coef) * term if out is None else out + float(coef) * term
    return out


def check_V(hw, Vm):
    """compare V Psi_j (full 36-term V on the reconstructed vector) with sum_i Vm[i][j] Psi_i, all entries; returns the max relative deviation
       (relative to ||V Psi_j||) over all entries and over the slice (site-2 = top state) only"""
    R = Rep.of(hw.lam)
    top = R.top
    m = hw.m
    Xs = [hw.hat_matrix(k) for k in range(m)]
    full = 0.0
    sl = 0.0
    for j in range(m):
        VX = apply_V_full(hw.lam, Xs[j])
        pred = sum(float(Vm[i][j]) * Xs[i] for i in range(m))
        nrm = spla.norm(VX)
        d = VX - pred
        full = max(full, (spla.norm(d) / nrm) if nrm > 0 else spla.norm(d))
        dc = d[:, top]
        sl = max(sl, (np.linalg.norm(dc.toarray()) / nrm) if nrm > 0 else 0.0)
    return full, sl


# ======================================================================================== 7. validation against W4's exact files
W4_REPS = [(1, 1, 0), (2, 1, 0), (2, 2, 0), (3, 1, 0), (3, 2, 0), (4, 1, 0)]


def w4_index_map(lam, w4):
    """W4 one-site index (0-based, notebook order = Reverse@Sort of the GT patterns) -> index in the GL3Rep order"""
    R = Rep.of(lam)
    return [R.index[tuple(tuple(int(x) for x in row) for row in p)] for p in w4['gtbasis']]


def _pair(p, D0, dps):
    x0 = mp.mpf(p[0].numerator) / p[0].denominator if isinstance(p[0], Fr) else mp.mpf(p[0])
    x1 = mp.mpf(p[1].numerator) / p[1].denominator if isinstance(p[1], Fr) else mp.mpf(p[1])
    return x0 + x1 * mp.sqrt(D0) if D0 else x0


def _fr_of(x):
    return x if isinstance(x, Fr) else Fr(int(x))


def w4_sector_cas2(lam):
    """{mu: cas2(mu)} from the SECTOR records of W4's file (nu_general.load_w4 does not keep it)"""
    import nu_general as ng
    fn = os.path.join(ng.W4DIR, 'untwisted-exact-%d%d%d.txt' % tuple(lam))
    out = {}
    for line in open(fn, encoding='utf8'):
        f = line.rstrip('\n').split('\t')
        if f[0] == 'SECTOR':
            out[tuple(int(x) for x in f[1].strip('{}').split(','))] = int(f[5])
    return out


def validate_w4(lam, dps=50, verbose=True):
    """all branches of W4's untwisted-exact-<lam>.txt:
         (i)   dim M_mu = W4 multiplicity;
         (ii)  W4's exact HWBASIS kets: their slices lie in M_mu and the exact reconstruction from the slice reproduces the whole ket EXACTLY;
         (iii) V_red eigenvalues vs W4's alpha_1:  alpha_1 = ev + eta_1/2  (max deviation, dps digits) -- and the exact matrix identity
               A1 = V_red + eta_1/2 in W4's HWBASIS slice basis (LAURMAT);
         (iv)  copy slices vs W4's cg (up to scale, dps digits);
         (v)   full copy vectors (float recursion) vs W4's vec up to scale, and at dps digits via exact HW vectors.
       returns a dict of maxima and a list of per-branch records."""
    import nu_general as ng
    lam = tuple(lam)
    R = Rep.of(lam)
    d = R.dim
    w4 = ng.load_w4(lam)
    w2m = w4_index_map(lam, w4)
    assert w2m == [int(x) for x in ng.w4_to_mine(lam, w4, type('C', (), {'rep': type('P', (), {'patterns': R.pats})()})())], 'index map differs from nu_general.w4_to_mine'
    old = mp.mp.dps
    mp.mp.dps = dps
    recs = []
    out = dict(lam=lam, nbranch=0, mult_mismatch=0, hwket_mismatch=0, hwket_not_in_prv=0, ev_dev=mp.mpf(0), laurmat_mismatch=0, slice_dev=mp.mpf(0),
               vec_dev_float=0.0, vec_dev_mp=mp.mpf(0), hw_exact_viol=0, n_hwket=0, cas2_mismatch=0, eta_mismatch=0, order_mismatch=0, gnorm_dev=mp.mpf(0))
    cas_w4 = w4_sector_cas2(lam)
    top = R.top
    try:
        for mu in w4['order']:
            sec = w4['sectors'][mu]
            p = prv_space(lam, mu)
            if p.m != sec['mult']:
                out['mult_mismatch'] += 1
                continue
            Vm = V_red(lam, mu, p)
            pos = {int(c): a for a, c in enumerate(p.W)}
            if cas_w4[mu] != cas2(mu):
                out['cas2_mismatch'] += 1
            # (ii) exact HW kets of W4
            kets = []
            for hb in sec['hwbasis']:
                ket = {}
                for idx, val in hb.items():
                    i, j = divmod(idx - 1, d)
                    ket[(w2m[i], w2m[j])] = _fr_of(val)
                kets.append(ket)
                out['n_hwket'] += 1
                sl = [Fr(0)] * p.w
                for (nu, q), val in ket.items():
                    if q == top:
                        sl[pos[nu]] = val
                try:
                    p.coords(sl)
                except ValueError:
                    out['hwket_not_in_prv'] += 1
                    continue
                psi, nc = hw_from_slice_exact(lam, mu, sl)
                psi = {k: x for k, x in psi.items() if x != 0}
                if psi != ket:
                    out['hwket_mismatch'] += 1
            # exact matrix identity A1 = V_red + eta1/2 in the W4 HWBASIS slice basis
            slc = []
            for ket in kets:
                sl = [Fr(0)] * p.w
                for (nu, q), val in ket.items():
                    if q == top:
                        sl[pos[nu]] = val
                slc.append(p.coords(sl))
            m = p.m
            Bm = [[slc[j][i] for j in range(m)] for i in range(m)]            # columns = W4 kets' coordinates in the PRV basis
            lm = sec['laurmat']
            A1m = [[_fr_of(x) for x in row] for row in lm[0]]
            eta1_lm = _fr_of(lm[4][0][0])
            # V_red in W4 basis: Bm^-1 Vm Bm  (exact, via sympy Rational matrices)
            import sympy
            sB = sympy.Matrix([[sympy.Rational(x.numerator, x.denominator) for x in row] for row in Bm])
            sV = sympy.Matrix([[sympy.Rational(x.numerator, x.denominator) for x in row] for row in Vm])
            Vw = sB.inv() * sV * sB
            sA = sympy.Matrix([[sympy.Rational(x.numerator, x.denominator) for x in row] for row in A1m])
            sA1 = sA - sympy.Rational(eta1_lm.numerator, eta1_lm.denominator) / 2 * sympy.eye(m)
            ok_direct = (Vw - sA1) == sympy.zeros(m, m)
            ok_trans = (Vw - sA1.T) == sympy.zeros(m, m)
            if not ok_direct:
                out['laurmat_mismatch'] += 1
            # eigen-decomposition
            cps = copies(p, Vm, dps)
            exact_basis = None
            for br in sec['branches']:
                out['nbranch'] += 1
                D0 = br['D0']
                a1 = _pair(br['alpha'][0], D0, dps)
                e1 = _pair(br['eta'][0], D0, dps)
                if not (br['eta'][0][1] == 0 and _fr_of(br['eta'][0][0]) == eta1(lam, mu) == eta1_lm):
                    out['eta_mismatch'] += 1
                pred = a1 - e1 / 2
                k = min(range(m), key=lambda kk: abs(cps[kk]['ev'] - pred))
                if k != br['copy'] - 1:
                    out['order_mismatch'] += 1
                cp = cps[k]
                evdev = abs(cp['ev'] - pred)
                out['ev_dev'] = max(out['ev_dev'], evdev)
                # slice vs cg
                cgw = {w2m[nu - 1]: _pair((ca, cb), D0, dps) for (nu, ca, cb) in br['cg']}
                sd = mp.mpf(0)
                keys = set(cgw) | set(cp['slice'])
                i0 = min(keys, key=lambda c: -abs(cgw.get(c, 0)))
                sc = cgw.get(i0, 0) / cp['slice'].get(i0, mp.mpf('1e-300'))
                mx = max(abs(x) for x in cgw.values())
                for c in keys:
                    sd = max(sd, abs(cgw.get(c, 0) - sc * cp['slice'].get(c, 0)) / mx)
                out['slice_dev'] = max(out['slice_dev'], sd)
                # full vector, float
                vflt = np.array([float(cp['slice'].get(int(c), 0)) for c in p.W])
                hw = hw_from_slice(lam, mu, vflt)
                dense = hw.dense(0)
                wv = {}
                for (idx, va, vb) in br['vec']:
                    i, j = divmod(idx - 1, d)
                    wv[(w2m[i], w2m[j])] = _pair((va, vb), D0, dps)
                keys2 = set(wv) | set(zip(hw.rows.tolist(), hw.cols.tolist()))
                big = max(wv, key=lambda c: abs(wv[c]))
                scf = float(wv[big]) / dense[big]
                vm = max(abs(float(x)) for x in wv.values())
                vd = max(abs(float(wv.get(c, 0)) - scf * dense[c]) for c in keys2) / vm
                out['vec_dev_float'] = max(out['vec_dev_float'], vd)
                # full vector at dps digits: sum_j coords_j Psi_exact(b_j)
                if exact_basis is None:
                    exact_basis = [hw_from_slice_exact(lam, mu, p.basis[kk])[0] for kk in range(m)]
                    for kk in range(m):
                        out['hw_exact_viol'] += exact_hw_check(lam, exact_basis[kk])
                full = {}
                for kk in range(m):
                    ck = cp['coords'][kk]
                    for key, x in exact_basis[kk].items():
                        if x != 0:
                            full[key] = full.get(key, 0) + ck * mpf_of(x)
                big = max(wv, key=lambda c: abs(wv[c]))
                scm = wv[big] / full[big]
                vdm = max(abs(wv.get(c, 0) - scm * full.get(c, 0)) for c in set(wv) | set(full)) / max(abs(x) for x in wv.values())
                out['vec_dev_mp'] = max(out['vec_dev_mp'], vdm)
                # W4's exact G-norm^2 of its (unnormalised) vector = (scale W4/mine)^2 * G-norm^2 of ours (same scale as the slice)
                gn2 = mp.fsum(mpf_of(R.g[nu] * R.g[q]) * x * x for (nu, q), x in full.items())
                gpred = sc * sc * gn2
                gdev = abs(gpred - _pair(br['gnorm2'], D0, dps)) / abs(gpred)
                out['gnorm_dev'] = max(out['gnorm_dev'], gdev)
                recs.append(dict(mu=mu, copy=br['copy'], D0=D0, alpha1=a1, eta1=e1, ev=cp['ev'], ev_dev=evdev, slice_dev=sd, vec_dev_float=vd, vec_dev_mp=vdm,
                                 ratio_slice=sc))
            if verbose:
                print('  ', mu, 'm', m, 'LAURMAT A1 - eta1/2 == V_red (W4 basis): %s%s' % (ok_direct, '' if ok_direct else ' (transpose: %s)' % ok_trans), flush=True)
    finally:
        mp.mp.dps = old
    out['records'] = recs
    return out


# ======================================================================================== 8. tables, tests, benchmark
def cas2(mu):
    """W4's cas2: sum_i mu_i (mu_i + 4 - 2 i)"""
    return mu[0] * (mu[0] + 2) + mu[1] * mu[1] + mu[2] * (mu[2] - 2)


def eta1(lam, mu):
    """W4's eta_1 (theta12-scaled simple-pole datum of S2 at site 1) = |Lam|^2 - Omega = |Lam|^2 + cas2(Lam) - cas2(mu)/2, Omega = sum_ij E1_ij E2_ji
       = (cas2(mu) - 2 cas2(Lam))/2, |Lam| = l1+l2+l3;   alpha_1 = eigenvalue(V_red) + eta_1/2"""
    return Fr(sum(lam) ** 2 + cas2(lam)) - Fr(cas2(mu), 2)


def table1(lams, verbose=True):
    """R1: dim M_mu vs Racah-Speiser for every irrep (and every dominant candidate mu = Lam + w); returns rows"""
    rows = []
    for lam in lams:
        lam = tuple(lam)
        t0 = time.perf_counter()
        R = Rep.of(lam)
        rs = racah_speiser(lam)
        wm = weight_mults(lam)
        sumdim = sum(m * dim_gl3(mu) for mu, m in rs.items())
        cand = candidates(lam)
        mism = 0
        maxw = 0
        for mu in cand:
            p = prv_space(lam, mu)
            maxw = max(maxw, p.w)
            if p.m != rs.get(mu, 0):
                mism += 1
        extra = [mu for mu in rs if mu not in set(cand)]
        row = dict(lam=lam, dim=R.dim, nirrep=len(rs), summ=sum(rs.values()), maxm=max(rs.values()), maxw=maxw, maxnW=max(nW_of(lam, mu, wm) for mu in rs),
                   dimsum_ok=(sumdim == R.dim ** 2), ncand=len(cand), mism=mism + len(extra), t=time.perf_counter() - t0)
        rows.append(row)
        if verbose:
            print('  %-10s dim %5d  irreps %4d  sum m %5d  max m %2d  max w %2d  max nW %6d  sum m dimV = dim^2: %s  dominant candidates tested %4d  mismatches %d  (%.2fs)' %
                  (str(list(lam)), row['dim'], row['nirrep'], row['summ'], row['maxm'], row['maxw'], row['maxnW'], row['dimsum_ok'], row['ncand'], row['mism'], row['t']), flush=True)
    return rows


def float_slices(prv, cps):
    """float slices (w, m) of the copies (columns, ascending eigenvalue) on the weight space"""
    return np.array([[float(cp['slice'].get(int(c), 0)) for cp in cps] for c in prv.W]).reshape(prv.w, len(cps))


def test_hw(lams, big=None, verbose=True, check_v=True):
    """R2/R3 on all irreps: reconstruct the copies' vectors (float), HW residual (12,23,13), weight, V (36 terms) vs V_red; exact twin for the small reps"""
    out = []
    for lam in list(lams) + ([big] if big else []):
        lam = tuple(lam)
        t0 = time.perf_counter()
        rs = racah_speiser(lam)
        worst = 0.0
        wdev = 0
        worstV = 0.0
        worstVs = 0.0
        nvec = 0
        nex = nexchk = exviol = 0
        for mu in sorted(rs, reverse=True):
            p = prv_space(lam, mu)
            Vm = V_red(lam, mu, p)
            cps = copies(p, Vm, 40)
            hw = hw_from_slice(lam, mu, float_slices(p, cps))
            r, wd = hw_residual(hw)
            worst = max(worst, r)
            wdev = max(wdev, wd)
            nvec += hw.m
            if check_v:
                # V_red matrix in the basis of the copies is diagonal; compare the full V on the reconstructed copies with ev * copy
                fv, fs = check_V(hw, [[float(cps[j]['ev']) if i == j else 0.0 for j in range(p.m)] for i in range(p.m)])
                worstV = max(worstV, fv)
                worstVs = max(worstVs, fs)
            if lam[0] <= 8:
                for k in range(p.m):
                    psi, nc = hw_from_slice_exact(lam, mu, p.basis[k])
                    nex += 1
                    nexchk += nc
                    exviol += exact_hw_check(lam, psi)
        row = dict(lam=lam, nirrep=len(rs), nvec=nvec, resid=worst, wdev=wdev, vdev=worstV, vsl=worstVs, nex=nex, nexchk=nexchk, exviol=exviol, t=time.perf_counter() - t0)
        out.append(row)
        if verbose:
            print('  %-10s irreps %4d  copies %4d  max HW residual %.2e  max weight deviation %d  full 36-term V vs ev*copy (all entries) %.2e (slice only %.2e)  exact twin: %d vectors, %d exact consistency conditions, %d HW violations  (%.1fs)' %
                  (str(list(lam)), len(rs), nvec, worst, wdev, worstV, worstVs, nex, nexchk, exviol, row['t']), flush=True)
    return out


def neg_tests(lams, nrand=20, seed=1, verbose=True):
    """slices outside M_mu must be rejected by the exact recursion (exact consistency conditions); slices inside must pass"""
    import random
    rng = random.Random(seed)
    res = []
    for lam in lams:
        lam = tuple(lam)
        rs = racah_speiser(lam)
        ndet = ntot = nin = ninok = 0
        for mu in sorted(rs, reverse=True):
            p = prv_space(lam, mu)
            if p.w > p.m:
                for _ in range(nrand):
                    v = [Fr(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(p.w)]
                    try:
                        p.coords(v)
                        continue                                 # accidentally inside M_mu
                    except ValueError:
                        pass
                    ntot += 1
                    try:
                        hw_from_slice_exact(lam, mu, v)
                    except ValueError:
                        ndet += 1
            for _ in range(3):
                c = [Fr(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(p.m)]
                v = [sum(c[k] * p.basis[k][i] for k in range(p.m)) for i in range(p.w)]
                nin += 1
                try:
                    hw_from_slice_exact(lam, mu, v)
                    ninok += 1
                except ValueError:
                    pass
        res.append((lam, ntot, ndet, nin, ninok))
        if verbose:
            print('  %-10s random slices outside M_mu: %d, rejected by the exact consistency conditions: %d;  random slices inside M_mu: %d, accepted: %d' % (str(list(lam)), ntot, ndet, nin, ninok), flush=True)
    return res


def _clear_caches():
    Rep._cache.clear()
    Recursion._cache.clear()
    ExactRecursion._cache.clear()


def bench_one(a, reps=5):
    """R5 for Lam = [2a,a,0]: (a) = exact generators + RS + PRV + V_red + float eig for all irreps (no n_W object touched);
       (b) = (a) + metric/float matrices/left inverses + reconstruction of all copies' full vectors (batched per irrep).  min over `reps` cold runs."""
    lam = (2 * a, a, 0)
    best = None
    for _ in range(reps):
        _clear_caches()
        T = {}
        t0 = time.perf_counter()
        R = Rep.of(lam)
        T['gen'] = time.perf_counter() - t0
        t = time.perf_counter()
        rs = racah_speiser(lam)
        wm = weight_mults(lam)
        mus = sorted(rs, reverse=True)
        T['rs'] = time.perf_counter() - t
        t = time.perf_counter()
        prvs = {mu: prv_space(lam, mu) for mu in mus}
        T['prv'] = time.perf_counter() - t
        t = time.perf_counter()
        Vs = {mu: V_red(lam, mu, prvs[mu]) for mu in mus}
        T['vred'] = time.perf_counter() - t
        t = time.perf_counter()
        eig = {}
        for mu in mus:
            Vf = np.array([[float(x) for x in row] for row in Vs[mu]])
            if prvs[mu].m == 1:
                eig[mu] = (Vf[0], np.eye(1), 0.0)
            else:
                w_, U = np.linalg.eig(Vf)
                eig[mu] = (w_.real, U.real, np.abs(w_.imag).max() / max(1.0, np.abs(w_).max()))
        T['eig'] = time.perf_counter() - t
        T['a'] = time.perf_counter() - t0
        t = time.perf_counter()
        rec = Recursion.of(lam)                  # metric, orthonormal float matrices, left inverses of all A_ab
        T['setup_b'] = time.perf_counter() - t
        t = time.perf_counter()
        nent = 0
        for mu in mus:
            p = prvs[mu]
            B = np.array([[float(x) for x in b] for b in p.basis]).T
            sl = B @ eig[mu][1]
            hw = hw_from_slice(lam, mu, sl)
            nent += hw.nW * hw.m
        T['recon'] = time.perf_counter() - t
        T['b'] = T['a'] + T['setup_b'] + T['recon']
        T['nent'] = nent
        if best is None or T['b'] < best['b']:
            best = T
    # accuracy pass (not timed): worst HW residual over all irreps (E12, E23 only); float-eigen quality
    worst = 0.0
    mingap = np.inf
    maximag = 0.0
    eigres = 0.0
    for mu in mus:
        p = prvs[mu]
        B = np.array([[float(x) for x in b] for b in p.basis]).T
        hw = hw_from_slice(lam, mu, B @ eig[mu][1])
        worst = max(worst, hw_residual(hw, extra_E13=False)[0])
        if p.m > 1:
            ev = np.sort(eig[mu][0])
            mingap = min(mingap, np.min(np.diff(ev) / np.maximum(1.0, np.abs(ev[:-1]))))
            maximag = max(maximag, eig[mu][2])
            Vf = np.array([[float(x) for x in row] for row in Vs[mu]])
            for k in range(p.m):
                u = eig[mu][1][:, k]
                eigres = max(eigres, np.linalg.norm(Vf @ u - eig[mu][0][k] * u) / np.linalg.norm(u))
    best.update(a=a, dim=R.dim, nirrep=len(rs), ncopies=sum(rs.values()), maxw=max(p.w for p in prvs.values()), maxnW=max(nW_of(lam, mu, wm) for mu in mus),
                resid=worst, smin=rec.smin, ta=best['a'], tb=best['b'], mingap=mingap, maximag=maximag, eigres=eigres)
    return best


def bench(amax=14, reps=5, verbose=True, amin=1):
    rows = []
    if verbose:
        print('  %2s %5s %6s %6s %5s %7s | %7s %7s %7s %7s %7s | %7s | %7s %7s | %9s | %8s' %
              ('a', 'dimV', 'irreps', 'copies', 'max w', 'max nW', 't_gen', 't_rs', 't_prv', 't_Vred', 't_eig', 't_a', 't_setup', 't_recon', 't_b', 'resid'), flush=True)
    for a in range(amin, amax + 1):
        T = bench_one(a, reps)
        rows.append(T)
        if verbose:
            print('  %2d %5d %6d %6d %5d %7d | %7.3f %7.3f %7.3f %7.3f %7.3f | %7.3f | %7.3f %7.3f | %9.3f | %8.1e   entries %d  [float eig: min rel. gap %.3f, max imag %.1e, max eigen-residual %.1e; min A_ab sv ratio %.3f]' %
                  (a, T['dim'], T['nirrep'], T['ncopies'], T['maxw'], T['maxnW'], T['gen'], T['rs'], T['prv'], T['vred'], T['eig'], T['ta'], T['setup_b'], T['recon'], T['tb'], T['resid'], T['nent'],
                   T['mingap'], T['maximag'], T['eigres'], T['smin']), flush=True)
    return rows


def show_irreps(lam, mus):
    """R2 table for named irreps: m, w, nW, HW residual (E12,E23,E13) of the reconstructed copies, weight deviation, 36-term V vs ev*copy"""
    lam = tuple(lam)
    rs = racah_speiser(lam)
    wm = weight_mults(lam)
    for mu in mus:
        mu = tuple(mu)
        if mu not in rs:
            print('  %s: not an irrep of V (x) V' % (mu,))
            continue
        p = prv_space(lam, mu)
        Vm = V_red(lam, mu, p)
        cps = copies(p, Vm, 40)
        hw = hw_from_slice(lam, mu, float_slices(p, cps))
        r, wd = hw_residual(hw)
        fv, fs = check_V(hw, [[float(cps[j]['ev']) if i == j else 0.0 for j in range(p.m)] for i in range(p.m)])
        print('  %-10s m %d  w %d  nW %5d  HW residual %.1e  weight deviation %d  V(36 terms) vs ev*copy %.1e  eigenvalues %s' %
              (str(list(mu)), p.m, p.w, nW_of(lam, mu, wm), r, wd, fv, [mp.nstr(c['ev'], 8) for c in cps]), flush=True)


def reference_irrep(lam, mu, dps=30):
    """everything for one irrep: PRV basis, exact V_red, copies (mpmath, dps digits), and the float full vectors of the copies (HW object, plain GT basis)"""
    p = prv_space(lam, mu)
    Vm = V_red(lam, mu, p)
    cps = copies(p, Vm, dps)
    return dict(prv=p, Vm=Vm, copies=cps, hw=hw_from_slice(lam, mu, float_slices(p, cps)), eta1=eta1(tuple(lam), tuple(mu)))


def tridiag_check(lams, verbose=True):
    """structure of the V_red formula as an operator on the FULL weight space V[mu-Lam] (GT order of the weight space): bandwidth, and the signs of the
       off-diagonal products M[i,i+1] M[i+1,i].  (The formula is only valid on M_mu, but as a definition it makes sense on all of V[mu-Lam]; for w = m the two
       coincide, and a tridiagonal matrix with strictly positive off-diagonal products has a real SIMPLE spectrum.)"""
    out = []
    for lam in lams:
        lam = tuple(lam)
        R = Rep.of(lam)
        rs = racah_speiser(lam)
        maxband = 0
        ntri = 0
        cnt = {'eq': [0, 0, 0, 0], 'gt': [0, 0, 0, 0]}
        for mu, m in rs.items():
            W = R.wsp[(mu[0] - lam[0], mu[1] - lam[1], mu[2] - lam[2])]
            pos = {int(c): a for a, c in enumerate(W)}
            w = len(W)
            M = [[Fr(0)] * w for _ in range(w)]
            for j, c in enumerate(W):
                for r, x in vred_apply(R, {int(c): Fr(1)}).items():
                    M[pos[r]][j] = x
            band = max([abs(i - j) for i in range(w) for j in range(w) if M[i][j] != 0] + [0])
            maxband = max(maxband, band)
            ntri += (band <= 1)
            key = 'eq' if w == m else 'gt'
            cnt[key][3] += 1
            for i in range(w - 1):
                pr = M[i][i + 1] * M[i + 1][i]
                cnt[key][0 if pr > 0 else (1 if pr == 0 else 2)] += 1
        out.append((lam, len(rs), maxband, ntri, cnt))
        if verbose:
            print('  %-10s sectors %3d: max bandwidth of V_red on the full weight space V[mu-Lam] = %d (0 = diagonal, 1 = tridiagonal), tridiagonal in %d of %d sectors; '
                  'off-diagonal products (positive, zero, negative): sectors with w = m (%d): %s, sectors with w > m (%d): %s' %
                  (str(list(lam)), len(rs), maxband, ntri, len(rs), cnt['eq'][3], tuple(cnt['eq'][:3]), cnt['gt'][3], tuple(cnt['gt'][:3])), flush=True)
    return out


def closed_forms(amax=10, verbose=True):
    """closed forms seen in the exact V_red: for Lam = [2a,a,0] the singlet (2a,2a,2a) has V = 2a^2 (2a+1) = l1 l2 (l1+1), the top irrep (4a,2a,0) has
       V = -l1 l2 (alpha_1 = 0, eta_1 = 2 l1 l2).  Returns the number of mismatches."""
    bad = 0
    for a in range(1, amax + 1):
        lam = (2 * a, a, 0)
        sing = V_red(lam, (2 * a, 2 * a, 2 * a))
        top = V_red(lam, (4 * a, 2 * a, 0))
        ok = (sing == [[Fr(2 * a * a * (2 * a + 1))]]) and (top == [[Fr(-2 * a * a)]]) and eta1(lam, (4 * a, 2 * a, 0)) == 2 * lam[0] * lam[1]
        bad += (not ok)
        if verbose:
            print('  a = %2d  Lam = %-10s singlet V = %s (2a^2(2a+1) = %d)  top V = %s (-l1 l2 = %d)  eta_1(top) = %s (2 l1 l2 = %d)  %s' %
                  (a, str(list(lam)), sing[0][0], 2 * a * a * (2 * a + 1), top[0][0], -2 * a * a, eta1(lam, (4 * a, 2 * a, 0)), 4 * a * a, 'ok' if ok else 'MISMATCH'), flush=True)
    return bad


def rounded_vs_float(a=14, mus=((37, 24, 23), (39, 23, 22), (40, 22, 22), (28, 28, 28)), verbose=True):
    """HW residual of the float64 recursion vs the exact recursion rounded to double (first PRV basis vector of each irrep)"""
    lam = (2 * a, a, 0)
    out = []
    for mu in mus:
        p = prv_space(lam, mu)
        t = time.perf_counter()
        hw_r = hw_from_slice_rounded(lam, mu, p.basis[0])
        te = time.perf_counter() - t
        t = time.perf_counter()
        hw_f = hw_from_slice(lam, mu, np.array([float(x) for x in p.basis[0]]))
        tf = time.perf_counter() - t
        rr = hw_residual(hw_r, extra_E13=False)[0]
        rf = hw_residual(hw_f, extra_E13=False)[0]
        out.append((mu, hw_r.nW, rf, rr, tf, te))
        if verbose:
            print('  a = %d  mu = %-12s m %d  nW %6d  HW residual: float64 recursion %.1e (%.3fs), exact recursion rounded to double %.1e (%.2fs)' % (a, str(list(mu)), p.m, hw_r.nW, rf, tf, rr, te), flush=True)
    return out


def exact_all(a):
    """exact twin on ALL PRV basis vectors of Lam = [2a,a,0]: exact HW vectors, exact consistency conditions, exact HW check (E12,E23,E13)"""
    lam = (2 * a, a, 0)
    rs = racah_speiser(lam)
    t0 = time.perf_counter()
    nvec = nchk = bad = ent = 0
    for mu in sorted(rs, reverse=True):
        p = prv_space(lam, mu)
        for k in range(p.m):
            psi, nc = hw_from_slice_exact(lam, mu, p.basis[k])
            nchk += nc
            nvec += 1
            ent += len(psi)
            bad += exact_hw_check(lam, psi)
    print('  a = %d, Lam = %s: exact twin on ALL %d PRV basis vectors (%d entries): exact consistency conditions verified %d, exact HW violations (E12,E23,E13) %d, %.1fs' %
          (a, list(lam), nvec, ent, nchk, bad, time.perf_counter() - t0), flush=True)
    return nvec, ent, nchk, bad


def eig_check(lam, dps=50):
    """float64 eigen-decomposition of the exact V_red (used in the benchmark) vs mpmath: max eigenvalue deviation, and ordering of the copies"""
    lam = tuple(lam)
    rs = racah_speiser(lam)
    worst = 0.0
    nc = 0
    for mu in sorted(rs, reverse=True):
        p = prv_space(lam, mu)
        if p.m < 2:
            continue
        Vm = V_red(lam, mu, p)
        cps = copies(p, Vm, dps)
        Vf = np.array([[float(x) for x in row] for row in Vm])
        ev = np.sort(np.linalg.eigvals(Vf).real)
        for k in range(p.m):
            worst = max(worst, abs(ev[k] - float(cps[k]['ev'])) / max(1.0, abs(float(cps[k]['ev']))))
        nc += p.m
    print('  %s: float64 eigenvalues of V_red vs mpmath (%d digits), %d copies in sectors with m >= 2: max relative deviation %.2e' % (list(lam), dps, nc, worst), flush=True)
    return worst


# ======================================================================================== 9. reference files
REFDIR = os.path.normpath(os.path.join(HERE, '..', '..', 'Mathematica', 'Gaudin', 'ref'))
import sympy
_X = sympy.Symbol('x')


def fq(x):
    x = Fr(x)
    return str(x.numerator) if x.denominator == 1 else '%d/%d' % (x.numerator, x.denominator)


def f30(x):
    return mp.nstr(x, 30, min_fixed=-mp.inf, max_fixed=mp.inf)


def fmu(mu):
    return '{%d,%d,%d}' % tuple(mu)


def _hdr(lam, vec, date):
    l1, l2, l3 = lam
    name = 'cg-slices-%d%d%d.txt' % lam
    L = []
    A = L.append
    A('# %s   (WP-R of round 6; generator Paul/Python/gaudin_sl3/slice_engine.py `python -X utf8 slice_engine.py ref %d %d`; %s)' % (name, l1, l2, date))
    A('# Highest-weight (HW) Clebsch-Gordan data of gl(3), V_Lambda (x) V_Lambda with Lambda = {%d,%d,%d} (dim V = DIM): for EVERY irrep mu of the decomposition the PRV/Zhelobenko slice space M_mu, the exact matrix of the Gaudin vertex operator V on it, and its eigen-copies (the SoV copies).' % lam)
    A('# FORMAT: one record per line, TAB-separated fields, first field = tag, other fields nested braces = lists (Wolfram InputForm style). p/q = exact rationals. Reals marked "30 digits" are plain decimals (no exponent) with up to 30 significant digits, trailing zeros stripped, so exact values like 2.5 appear short (Mathematica: SetPrecision[ToExpression[s],30]; python: mpmath.mpf(s)); residuals are Python floats (e.g. 1.2e-15). Lines starting with # are comments.')
    A('# ORDER OF THE ONE-SITE GT BASIS (binding for every index below): the order of Paul/Python/gl3_reps.GL3Rep(Lambda).patterns = ascending lexicographic in (m1,m2,k) of the pattern {{k},{m1,m2},{l1,l2,l3}}; the HIGHEST-WEIGHT state is the LAST index (= DIM). All indices nu, idx are 1-BASED in this order. GTBASIS lists the patterns in exactly this order (same record format as W4).')
    A('#   The notebook/W4 order (ref/untwisted-exact-*.txt, Reverse@Sort of the patterns = descending lexicographic in (k,m1,m2), index 1 = highest weight) is DIFFERENT: map by pattern (python: nu_general.w4_to_mine). Matrices act on kets as columns, E e_col = sum_row M[row,col] e_row (plain rational non-unitary GT matrices of gl3_reps), E13 = [E12,E23], E31 = [E32,E21]. WEIGHTS1 = (E11,E22,E33) eigenvalues per index; GMET1 = invariant diagonal metric G (E_ij^T G = G E_ji, G[HW] = 1; equals W4 GMET1 up to the order).')
    A('# TWO SITES: Psi = sum Psi[nu,q] e_nu (x) e_q, site 1 = first factor (nu, slow), site 2 = second factor (q, fast); 1-based flat index idx = (nu-1)*DIM + q (nu,q 1-based, GL3Rep order). The SLICE of a vector Psi of weight mu is psi_top[nu] = Psi[nu, q_top] (site 2 = HW state q_top = index DIM), a vector in the weight space V_Lambda[mu-Lambda] (weight space = all nu whose GT weight equals mu-Lambda; W4 record `cg` = the same slice).')
    A('# THEOREM USED (Zhelobenko/PRV): Psi -> slice is an isomorphism from the HW vectors of weight mu onto M_mu = { v in V_Lambda[mu-Lambda] : E12^(l1-l2+1) v = 0, E23^(l2-l3+1) v = 0 }; dim M_mu = m = multiplicity of mu (Racah-Speiser, checked for every irrep). The full vector is recovered from the slice by a triangular recursion over the weights of site 2 (slice_engine.py, exact twin available).')
    A('# VERTEX OPERATOR: V = 1/2 Sum_{s,t in S3} sgn(s) sgn(t) e1_{s1 t1} e1_{s2 t2} e2_{s3 t3} (36 terms, e^k_ij = E_ij on site k, product in the written order), the notebook V (Minimal_complete_working-fast.wl vTriples). On a slice v in M_mu: V_red v = 1/2 Sum sgn sgn E_{s1 t1} E_{s2 t2} Phi_{s3 t3} v, Phi_cc = lambda_c, Phi_ab = -E_ab (a<b), Phi_ab = 0 (a>b) (no n_W-sized object needed; verified against the full 36-term V on the reconstructed vectors).')
    A('# RELATION TO W4/Laurent DATA (exact on all 45 branches of W4 files 110,210,220,310,320,410; see RESULTS-slice-engine.md): alpha_1(W4) = eigenvalue(V_red) + eta1/2, eta1 = |Lambda|^2 + cas2(Lambda) - cas2(mu)/2 = |Lambda|^2 - Omega with Omega = Sum_ij E1_ij E2_ji = (cas2(mu) - 2 cas2(Lambda))/2, |Lambda| = l1+l2+l3, cas2(mu) = Sum_i mu_i (mu_i + 4 - 2 i) (= W4 cas2); i.e. eigenvalue = alpha_1 - eta1/2 (alpha_1 = theta12 * a_1, a_1 = double-pole coefficient of S3 at theta_1). copy 1 = smaller eigenvalue = smaller alpha_1 (W4 convention).')
    A('# RECORDS: REP {l1,l2,l3} | DIM | GTBASIS | WEIGHTS1 | GMET1 | NSECTORS | NCOPIES | then per sector (descending lexicographic mu): SECTOR, WSPACE, PRVBASIS x m, VRED, CHARPOLY, then per copy COPY, HWRESID, GNORM2, then VEC x m (small reps) |')
    A('#   SECTOR mu m w nW dimVmu cas2mu eta1 : m = multiplicity = dim M_mu; w = dim V_Lambda[mu-Lambda] = length of a slice (NOT W4s third field); nW = dim of the weight-mu subspace of V (x) V (size of a full vector); dimVmu = dim of the irrep mu; cas2mu, eta1 as above |')
    A('#   PRVBASIS mu {{nu,val},...} : one line per basis vector of M_mu, k = 1..m in line order; EXACT rationals; nonzero entries of the slice on V_Lambda[mu-Lambda]; canonical basis = RREF nullspace with columns in ascending GT index, basis vector k = 1 at the k-th free column and 0 at the other free columns |')
    A('#   WSPACE mu {nu,...} : the 1-based GT indices (ascending) of the weight space V_Lambda[mu-Lambda]; the slice/PRVBASIS positions are aligned with this list |')
    A('#   VRED mu {{...},...} : exact m x m matrix of V_red in the PRVBASIS (list of ROWS; column j = V_red b_j in PRV coordinates, i.e. V_red b_j = Sum_i VRED[i][j] b_i) |')
    A('#   CHARPOLY mu {c0,c1,...,cm} : exact characteristic polynomial det(x 1 - V_red) = Sum_k c_k x^k (c_m = 1; ascending powers, exact rationals); its roots are the eigenvalues of the COPY records (= the V-spectrum, the pencil roots of the Q-system) |')
    A('#   COPY mu k eigenvalue {{nu,val},...} : the SoV copies = eigenvectors of V_red; k = 1..m by INCREASING eigenvalue (30 digits); the slice (30 digits) is scaled so that its FIRST NONZERO component (ascending nu) is exactly 1; components below 1e-40 relative to the largest (exact zeros of the algebraic eigenvector at 60 working digits) are dropped |')
    A('#   HWRESID mu k residual : max over i in {12,23,13} of ||(E_i (x) 1 + 1 (x) E_i) Psi|| / ||Psi|| for the full vector reconstructed (floating point, orthonormal GT basis, ||.|| = invariant norm) from COPY k |')
    A('#   GNORM2 mu k value : G-norm^2 (invariant norm, sum_{nu,q} g_nu g_q |Psi[nu,q]|^2, double precision, e-notation) of the FULL vector whose slice is COPY k (slice scaled with first nonzero component 1); multiply by the squared scale to compare with W4 gnorm2 |')
    if vec:
        A('#   VEC mu k {{idx,val},...} : the full HW vector of copy k (30 digits; exact HW vectors of the PRV basis combined with the 30-digit eigenvector; same scale as COPY k, i.e. its slice component at q_top is the COPY slice), nonzero entries only, idx 1-based flat index (nu-1)*DIM + q |')
    else:
        A('#   (no VEC records in this file: full vectors have nW up to the table value; regenerate with slice_engine.hw_from_slice / the exact twin)')
    return L


def write_ref(lam, vec=False, dps=60, date='2026-10-21', verbose=True):
    lam = tuple(lam)
    t0 = time.perf_counter()
    R = Rep.of(lam)
    d = R.dim
    rs = racah_speiser(lam)
    wm = weight_mults(lam)
    mus = sorted(rs, reverse=True)
    body = []
    stats = dict(maxres=0.0, maxeig=mp.mpf(0), mingap=None, ncopies=0, vecdev=0.0, nrep=0)
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        for mu in mus:
            p = prv_space(lam, mu)
            assert p.m == rs[mu]
            Vm = V_red(lam, mu, p)
            cps = copies(p, Vm, dps)
            m = p.m
            body.append('SECTOR\t%s\t%d\t%d\t%d\t%d\t%d\t%s' % (fmu(mu), m, p.w, nW_of(lam, mu, wm), dim_gl3(mu), cas2(mu), fq(eta1(lam, mu))))
            body.append('WSPACE\t%s\t{%s}' % (fmu(mu), ','.join(str(int(c) + 1) for c in p.W)))
            for k in range(m):
                v = p.vec(k)
                body.append('PRVBASIS\t%s\t{%s}' % (fmu(mu), ','.join('{%d,%s}' % (c + 1, fq(x)) for c, x in sorted(v.items()))))
            body.append('VRED\t%s\t{%s}' % (fmu(mu), ','.join('{%s}' % ','.join(fq(x) for x in row) for row in Vm)))
            body.append('CHARPOLY\t%s\t{%s}' % (fmu(mu), ','.join(fq(Fr(int(c.p), int(c.q))) for c in reversed(sympy.Matrix([[sympy.Rational(x.numerator, x.denominator) for x in row] for row in Vm]).charpoly(_X).all_coeffs()))))
            # float reconstruction of the copies for the residuals
            hw = hw_from_slice(lam, mu, float_slices(p, cps))
            res, wd = hw_residual(hw, per_vector=True)
            assert wd == 0
            # eigen-residuals at dps digits: Vm c = ev c (coordinates)
            Vmp = mp.matrix([[mpf_of(x) for x in row] for row in Vm])
            for k, cp in enumerate(cps):
                c = mp.matrix(cp['coords'])
                r = Vmp * c - cp['ev'] * c
                stats['maxeig'] = max(stats['maxeig'], max(abs(x) for x in r) / max(abs(x) for x in c))
                stats['maxres'] = max(stats['maxres'], res[k])
                body.append('COPY\t%s\t%d\t%s\t{%s}' % (fmu(mu), k + 1, f30(cp['ev']), ','.join('{%d,%s}' % (c_ + 1, f30(x)) for c_, x in sorted(cp['slice'].items()))))
                body.append('HWRESID\t%s\t%d\t%.2e' % (fmu(mu), k + 1, res[k]))
                body.append('GNORM2\t%s\t%d\t%.15e' % (fmu(mu), k + 1, float(np.sum(hw.vals_hat[:, k] ** 2))))
                stats['ncopies'] += 1
            for k in range(m - 1):
                gap = abs(cps[k + 1]['ev'] - cps[k]['ev'])
                scale = max(1, abs(cps[k]['ev']))
                g = float(gap / scale)
                stats['mingap'] = g if stats['mingap'] is None else min(stats['mingap'], g)
            if vec:
                exb = [hw_from_slice_exact(lam, mu, p.basis[kk])[0] for kk in range(m)]
                for k, cp in enumerate(cps):
                    full = {}
                    for kk in range(m):
                        ck = cp['coords'][kk]
                        for key, x in exb[kk].items():
                            if x != 0:
                                full[key] = full.get(key, 0) + ck * mpf_of(x)
                    mx = max(abs(x) for x in full.values())
                    ent = sorted(((nu * d + q + 1, x) for (nu, q), x in full.items() if abs(x) > mx * mp.mpf(10) ** (-(dps - 15))))
                    # compare with the float recursion
                    dense = hw.dense(k)
                    sc = float(max(abs(x) for x in full.values()))
                    dev = max(abs(float(full.get((i, j), 0)) - dense[i, j]) for (i, j) in zip(hw.rows.tolist(), hw.cols.tolist())) / sc
                    stats['vecdev'] = max(stats['vecdev'], dev)
                    body.append('VEC\t%s\t%d\t{%s}' % (fmu(mu), k + 1, ','.join('{%d,%s}' % (i, f30(x)) for i, x in ent)))
    finally:
        mp.mp.dps = old
    head = _hdr(lam, vec, date)
    head.append('# SUMMARY (self-check of this file): %d sectors, %d copies; max HWRESID %.2e; max eigen-residual |V_red c - ev c|/|c| %s (%d digits); min relative gap between consecutive eigenvalues of a sector %s%s' %
                (len(mus), stats['ncopies'], stats['maxres'], mp.nstr(stats['maxeig'], 3), dps, ('%.3e' % stats['mingap']) if stats['mingap'] is not None else 'n/a (all m = 1)',
                 ('; max |VEC - float recursion|/|VEC| %.2e' % stats['vecdev']) if vec else ''))
    head.append('REP\t{%d,%d,%d}' % lam)
    head.append('DIM\t%d' % d)
    head.append('GTBASIS\t{%s}' % ','.join('{%s}' % ','.join('{%s}' % ','.join(str(x) for x in row) for row in pat) for pat in R.pats))
    head.append('WEIGHTS1\t{%s}' % ','.join('{%d,%d,%d}' % w for w in R.wt))
    head.append('GMET1\t{%s}' % ','.join(fq(x) for x in R.g))
    head.append('NSECTORS\t%d' % len(mus))
    head.append('NCOPIES\t%d' % stats['ncopies'])
    path = os.path.join(REFDIR, 'cg-slices-%d%d%d.txt' % lam)
    with open(path, 'w', encoding='utf8', newline='\n') as f:
        f.write('\n'.join(head + body) + '\n')
    stats['t'] = time.perf_counter() - t0
    stats['path'] = path
    stats['bytes'] = os.path.getsize(path)
    if verbose:
        print('  wrote %s: %d sectors, %d copies, max HWRESID %.2e, max eigen-residual %s, min gap %s, VEC dev %s, %d bytes, %.1fs' %
              (path, len(mus), stats['ncopies'], stats['maxres'], mp.nstr(stats['maxeig'], 3), stats['mingap'], ('%.2e' % stats['vecdev']) if vec else '-', stats['bytes'], stats['t']), flush=True)
    return stats


# ======================================================================================== 10. independent re-verification of the reference files
import re
_TOK = re.compile(r'\{|\}|,|-?\d+/\d+|-?\d+\.\d+(?:[eE][-+]?\d+)?|-?\d+')


def parse_nested(s):
    """parse '{{1,2/3},{0.5}}' -> nested lists of int / Fraction / mpf (no eval)"""
    toks = _TOK.findall(s)
    pos = [0]

    def val():
        t = toks[pos[0]]
        pos[0] += 1
        if t == '{':
            out = []
            while toks[pos[0]] != '}':
                if toks[pos[0]] == ',':
                    pos[0] += 1
                    continue
                out.append(val())
            pos[0] += 1
            return out
        if '/' in t:
            a, b = t.split('/')
            return Fr(int(a), int(b))
        if '.' in t or 'e' in t.lower():
            return mp.mpf(t)
        return int(t)
    return val()


def parse_ref(path, dps=50):
    """parse a cg-slices-*.txt: dict(lam, dim, gtbasis, weights, g, sectors={mu: dict(m,w,nW,dimV,cas2,eta1,prv=[{nu:val}],vred,copies=[(ev,{nu:val})],hwres=[],vec=[{idx:val}])})"""
    old = mp.mp.dps
    mp.mp.dps = dps
    out = dict(sectors={})
    cur = None
    try:
        for line in open(path, encoding='utf8'):
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            f = line.split('\t')
            tag = f[0]
            if tag == 'REP':
                out['lam'] = tuple(parse_nested(f[1]))
            elif tag == 'DIM':
                out['dim'] = int(f[1])
            elif tag == 'GTBASIS':
                out['gtbasis'] = parse_nested(f[1])
            elif tag == 'WEIGHTS1':
                out['weights'] = [tuple(w) for w in parse_nested(f[1])]
            elif tag == 'GMET1':
                out['g'] = [Fr(x) for x in parse_nested(f[1])]
            elif tag == 'SECTOR':
                mu = tuple(parse_nested(f[1]))
                cur = out['sectors'][mu] = dict(m=int(f[2]), w=int(f[3]), nW=int(f[4]), dimV=int(f[5]), cas2=int(f[6]), eta1=Fr(f[7]), prv=[], copies=[], hwres=[], vec=[])
            elif tag == 'WSPACE':
                cur['wspace'] = parse_nested(f[2])
            elif tag == 'CHARPOLY':
                cur['charpoly'] = [Fr(x) for x in parse_nested(f[2])]
            elif tag == 'PRVBASIS':
                cur['prv'].append({nu: Fr(x) for nu, x in parse_nested(f[2])})
            elif tag == 'VRED':
                cur['vred'] = [[Fr(x) for x in row] for row in parse_nested(f[2])]
            elif tag == 'COPY':
                cur['copies'].append((mp.mpf(f[3]), {nu: x for nu, x in parse_nested(f[4])}))
            elif tag == 'HWRESID':
                cur['hwres'].append(float(f[3]))
            elif tag == 'GNORM2':
                cur.setdefault('gnorm2', []).append(float(f[3]))
            elif tag == 'VEC':
                cur['vec'].append({i: x for i, x in parse_nested(f[3])})
    finally:
        mp.mp.dps = old
    return out


def _gl3_float_ops(lam):
    """independent float matrices of E_ij (all 9) from gl3_reps.GL3Rep (sympy -> float), E13 = [E12,E23], E31 = [E32,E21]"""
    from gl3_reps import GL3Rep
    G = GL3Rep(lam)
    E = {}
    for ij in ((1, 2), (2, 3), (2, 1), (3, 2)):
        E[ij] = np.array(G.generator(*ij).tolist(), dtype=float)
    E[(1, 3)] = E[(1, 2)] @ E[(2, 3)] - E[(2, 3)] @ E[(1, 2)]
    E[(3, 1)] = E[(3, 2)] @ E[(2, 1)] - E[(2, 1)] @ E[(3, 2)]
    for k in (1, 2, 3):
        E[(k, k)] = np.diag([float(_weight(k, p)) for p in G.patterns])
    return G, E


def verify_ref(path, verbose=True, ncheck_big=None):
    """independent checks of a written reference file (parsed back from the TEXT): (1) GTBASIS = GL3Rep order, WEIGHTS1, GMET1 (E_ij^T G = G E_ji with GL3Rep
       float matrices); (2) PRVBASIS vectors satisfy the two chain conditions (GL3Rep matrices) and VRED is consistent: V_red from the generic formula agrees;
       (3) for every COPY: reconstruct the full vector by the float recursion from the FILE's slice, then in the plain GT basis with GL3Rep matrices: HW residual and
       the full 36-term V eigen-equation V Psi = ev Psi (G-norm), slice = COPY; (4) VEC records (if present) satisfy HW and V Psi = ev Psi with the independent
       gaudin_chain.Chain two-site operators (exact Fractions, kron) and have the COPY slice at site 2 = top."""
    D = parse_ref(path)
    lam = D['lam']
    G, E = _gl3_float_ops(lam)
    n = G.dim
    res = dict(path=os.path.basename(path), lam=lam, nsec=len(D['sectors']), ncopy=0)
    res['gtbasis_ok'] = [tuple(tuple(r) for r in p) for p in D['gtbasis']] == [tuple(tuple(r) for r in p) for p in G.patterns]
    res['weights_ok'] = D['weights'] == [tuple(_weight(k, p) for k in (1, 2, 3)) for p in G.patterns]
    g = np.array([float(x) for x in D['g']])
    sg = np.sqrt(g)
    # metric: E_ij^T G = G E_ji  (float, relative)
    dev = 0.0
    for (i, j) in ((1, 2), (2, 3), (1, 3)):
        L = E[(i, j)].T * g[None, :]          # (E^T G)[q,r] = E[r,q] g_r
        Rr = g[:, None] * E[(j, i)]           # (G E')[q,r] = g_q E'[q,r]
        dev = max(dev, np.abs(L - Rr).max() / np.abs(Rr).max())
    res['metric_dev'] = dev
    top = n - 1
    k1, k2 = lam[0] - lam[1] + 1, lam[1] - lam[2] + 1
    # (2) PRV conditions on the file's PRVBASIS, float, relative to the product of norms
    prv_dev = 0.0
    vred_dev = 0.0
    for mu, S in D['sectors'].items():
        for vd in S['prv']:
            v = np.zeros(n)
            for nu, x in vd.items():
                v[nu - 1] = float(x)
            a = np.linalg.matrix_power(E[(1, 2)], k1) @ v
            b = np.linalg.matrix_power(E[(2, 3)], k2) @ v
            sc = np.linalg.norm(v) * max(np.linalg.norm(np.linalg.matrix_power(E[(1, 2)], k1)), np.linalg.norm(np.linalg.matrix_power(E[(2, 3)], k2)))
            prv_dev = max(prv_dev, max(np.abs(a).max(), np.abs(b).max()) / sc)
    res['prv_dev'] = prv_dev
    # WSPACE = weight space of mu - Lambda; CHARPOLY(ev) = 0 at 50 digits for every COPY eigenvalue; CHARPOLY = charpoly(VRED) (exact)
    old = mp.mp.dps
    mp.mp.dps = 50
    ws_bad = 0
    cp_dev = mp.mpf(0)
    cp_exact_bad = 0
    for mu, S in D['sectors'].items():
        wexp = [c + 1 for c in range(n) if tuple(D['weights'][c]) == (mu[0] - lam[0], mu[1] - lam[1], mu[2] - lam[2])]
        if S['wspace'] != wexp:
            ws_bad += 1
        c = S['charpoly']
        for ev, sd in S['copies']:
            val = mp.fsum(mpf_of(ck) * ev ** k for k, ck in enumerate(c))
            scale = mp.fsum(abs(mpf_of(ck)) * abs(ev) ** k for k, ck in enumerate(c))
            cp_dev = max(cp_dev, abs(val) / scale)
        cexp = [Fr(int(q.p), int(q.q)) for q in reversed(sympy.Matrix([[sympy.Rational(x.numerator, x.denominator) for x in row] for row in S['vred']]).charpoly(_X).all_coeffs())]
        if cexp != c:
            cp_exact_bad += 1
    mp.mp.dps = old
    res.update(ws_bad=ws_bad, cp_dev=cp_dev, cp_exact_bad=cp_exact_bad)
    # V_red from the generic formula (GL3Rep float matrices) on the PRV basis vs VRED
    for mu, S in D['sectors'].items():
        m = S['m']
        vs = []
        for vd in S['prv']:
            v = np.zeros(n)
            for nu, x in vd.items():
                v[nu - 1] = float(x)
            vs.append(v)
        lam_ = lam
        for j in range(m):
            v = vs[j]
            out = np.zeros(n)
            for coef, ij1, ij2, ij3 in VTERMS:
                a_, b_ = ij3
                if a_ > b_:
                    continue
                phi = (lam_[a_ - 1] * v) if a_ == b_ else -(E[ij3] @ v)
                out += float(coef) * (E[ij1] @ (E[ij2] @ phi))
            pred = sum(float(S['vred'][i][j]) * vs[i] for i in range(m))
            vred_dev = max(vred_dev, np.abs(out - pred).max() / max(np.abs(pred).max(), 1e-300))
    res['vred_dev'] = vred_dev
    # (3) copies: HW + V eigen-equation in the plain basis, with GL3Rep float matrices
    worst_hw = worst_v = worst_sl = 0.0
    ncopy = 0
    gg = sg[:, None] * sg[None, :]
    for mu, S in D['sectors'].items():
        W = [c for c in range(n) if tuple(D['weights'][c]) == (mu[0] - lam[0], mu[1] - lam[1], mu[2] - lam[2])]
        sl_all = np.zeros((len(W), len(S['copies'])))
        posw = {c: a for a, c in enumerate(W)}
        for k, (ev, sd) in enumerate(S['copies']):
            for nu, x in sd.items():
                sl_all[posw[nu - 1], k] = float(x)
        hw = hw_from_slice(lam, mu, sl_all)
        for k, (ev, sd) in enumerate(S['copies']):
            X = np.zeros((n, n))
            X[hw.rows, hw.cols] = hw.vals[:, k]
            nx = np.linalg.norm(gg * X)
            for ij in ((1, 2), (2, 3), (1, 3)):
                Rm = E[ij] @ X + X @ E[ij].T
                worst_hw = max(worst_hw, np.linalg.norm(gg * Rm) / nx)
            VX = np.zeros((n, n))
            cache = {}
            for coef, ij1, ij2, ij3 in VTERMS:
                if ij3 not in cache:
                    cache[ij3] = X @ E[ij3].T
                VX += float(coef) * (E[ij1] @ (E[ij2] @ cache[ij3]))
            worst_v = max(worst_v, np.linalg.norm(gg * (VX - float(ev) * X)) / nx)
            col = X[:, top]
            # slice equals the COPY slice up to scale
            s_file = np.zeros(n)
            for nu, x in sd.items():
                s_file[nu - 1] = float(x)
            j0 = int(np.argmax(np.abs(s_file)))
            sc = col[j0] / s_file[j0]
            worst_sl = max(worst_sl, np.abs(col - sc * s_file).max() / np.abs(col).max())
            ncopy += 1
    res.update(ncopy=ncopy, hw_dev=worst_hw, v_dev=worst_v, slice_dev=worst_sl)
    # (4) VEC with the independent two-site Chain (exact kron operators)
    if any(S['vec'] for S in D['sectors'].values()):
        from gaudin_chain import Chain
        ch = Chain(lam=lam)
        assert ch.rep.patterns == G.patterns
        N2 = ch.dim
        vec_hw = vec_v = vec_sl = 0.0
        nv = 0
        for mu, S in D['sectors'].items():
            for k, vd in enumerate(S['vec']):
                x = [0.0] * N2
                for idx, val in vd.items():
                    x[idx - 1] = float(val)
                nx = max(abs(t) for t in x)
                hwr = 0.0
                for ij in ((1, 2), (2, 3), (1, 3)):
                    Etot = ch.Es[(1,) + ij] + ch.Es[(2,) + ij]
                    hwr = max(hwr, max(abs(t) for t in Etot.apply(x)))
                vec_hw = max(vec_hw, hwr / nx)
                # V from the chain operators
                Vx = [0.0] * N2
                for coef, ij1, ij2, ij3 in VTERMS:
                    y = ch.Es[(2,) + ij3].apply(x)
                    y = ch.Es[(1,) + ij2].apply(y)
                    y = ch.Es[(1,) + ij1].apply(y)
                    c = float(coef)
                    Vx = [a + c * b for a, b in zip(Vx, y)]
                ev = float(S['copies'][k][0])
                vec_v = max(vec_v, max(abs(a - ev * b) for a, b in zip(Vx, x)) / nx)
                # slice: site 2 = top -> entries nu*n + top
                sd = S['copies'][k][1]
                dev_s = 0.0
                for nu in range(n):
                    a = x[nu * n + top]
                    b = float(sd.get(nu + 1, 0))
                    dev_s = max(dev_s, abs(a - b))
                vec_sl = max(vec_sl, dev_s)
                nv += 1
        res.update(nvec=nv, vec_hw=vec_hw, vec_v=vec_v, vec_slice=vec_sl)
    if verbose:
        print('  %s: sectors %d copies %d | GTBASIS=GL3Rep order %s, WEIGHTS1 ok %s, metric E^T G = G E\' dev %.1e | PRVBASIS chain conditions rel. dev %.1e | WSPACE wrong %d | VRED vs generic formula (GL3Rep floats) dev %.1e | '
              'CHARPOLY(ev)/scale %s, CHARPOLY != charpoly(VRED) %d | copies: HW resid %.1e, 36-term V Psi = ev Psi dev %.1e, slice dev %.1e%s' %
              (res['path'], res['nsec'], res['ncopy'], res['gtbasis_ok'], res['weights_ok'], res['metric_dev'], res['prv_dev'], res['ws_bad'], res['vred_dev'], mp.nstr(res['cp_dev'], 2), res['cp_exact_bad'], res['hw_dev'], res['v_dev'], res['slice_dev'],
               (' | VEC (%d vectors, gaudin_chain two-site operators): HW %.1e, V Psi = ev Psi %.1e, slice entries dev %.1e' % (res['nvec'], res['vec_hw'], res['vec_v'], res['vec_slice'])) if 'nvec' in res else ''), flush=True)
    return res


def bruteforce_check(lam, mus=None, verbose=True, tol_null=1e-9):
    """INDEPENDENT of the PRV theorem and of the recursion: for each irrep mu build the constraint matrix of the total E12, E23 on the weight-mu subspace of
       V (x) V (n_W columns; orthonormal GT basis), find its null space by the eigen-decomposition of A^T A (what Munich/fastEig do), and compare:
       (i) dim null = Racah-Speiser m (and the next eigenvalue is >= 1e-3 away), (ii) our reconstructed copies lie in it (relative distance), span = m,
       (iii) the eigenvalues of the full 36-term V restricted to the null space equal ours.  Returns the worst numbers."""
    lam = tuple(lam)
    R = Rep.of(lam)
    n = R.dim
    Eh = R.float_ops()
    rs = racah_speiser(lam)
    mus = sorted(rs, reverse=True) if mus is None else [tuple(m) for m in mus]
    wt = R.wt
    cscs = {ij: Eh[ij].tocsc() for ij in ((1, 2), (2, 3))}
    out = dict(lam=lam, n=0, dimnull_bad=0, gap=np.inf, nullnorm=0.0, span_min=1.0, evdev=0.0, maxnW=0)
    t0 = time.perf_counter()
    for mu in mus:
        # pairs (nu,q) of weight mu
        pairs = []
        for nu in range(n):
            qw = (mu[0] - wt[nu][0], mu[1] - wt[nu][1], mu[2] - wt[nu][2])
            qs = R.wsp.get(qw)
            if qs is not None:
                pairs.extend((nu, int(q)) for q in qs)
        nW = len(pairs)
        out['maxnW'] = max(out['maxnW'], nW)
        pos = {pq: i for i, pq in enumerate(pairs)}
        rows_, cols_, vals_ = [], [], []
        for ci, (nu, q) in enumerate(pairs):
            for t, ij in enumerate(((1, 2), (2, 3))):
                C = cscs[ij]
                for k in range(C.indptr[nu], C.indptr[nu + 1]):
                    rows_.append(t * n * n + int(C.indices[k]) * n + q)
                    cols_.append(ci)
                    vals_.append(C.data[k])
                for k in range(C.indptr[q], C.indptr[q + 1]):
                    rows_.append(t * n * n + nu * n + int(C.indices[k]))
                    cols_.append(ci)
                    vals_.append(C.data[k])
        ur, inv = np.unique(np.array(rows_), return_inverse=True)
        A = sps.csr_matrix((vals_, (inv, cols_)), shape=(len(ur), nW))
        Gm = (A.T @ A).toarray()
        w_, U = np.linalg.eigh(Gm)
        m = rs[mu]
        nnull = int(np.sum(w_ < tol_null))
        if nnull != m:
            out['dimnull_bad'] += 1
        out['gap'] = min(out['gap'], w_[m] if m < nW else np.inf)
        N = U[:, :m]                                   # orthonormal null vectors (the m lowest)
        # our copies in the same pair indexing
        p = prv_space(lam, mu)
        Vm = V_red(lam, mu, p)
        cps = copies(p, Vm, 30)
        hw = hw_from_slice(lam, mu, float_slices(p, cps))
        idx = np.array([pos[(int(r), int(c))] for r, c in zip(hw.rows, hw.cols)])
        Psi = np.zeros((nW, m))
        Psi[idx, :] = hw.vals_hat
        Psi /= np.linalg.norm(Psi, axis=0)
        P = N.T @ Psi                                  # coordinates in the null space
        out['nullnorm'] = max(out['nullnorm'], max(np.linalg.norm(Psi[:, k] - N @ P[:, k]) for k in range(m)))
        out['span_min'] = min(out['span_min'], np.linalg.svd(P, compute_uv=False)[-1])
        # V restricted to the null space: N^T V N from the full 36-term V on each null vector
        VN = np.zeros((nW, m))
        for k in range(m):
            rr = np.array([pq[0] for pq in pairs])
            cc = np.array([pq[1] for pq in pairs])
            Xk = sps.csr_matrix((N[:, k], (rr, cc)), shape=(n, n))
            VX = apply_V_full(lam, Xk).tocsr()
            VN[:, k] = np.asarray(VX[rr, cc]).ravel()
        Vb = N.T @ VN
        evb = np.sort(np.linalg.eigvals(Vb).real)
        evo = np.sort([float(cp['ev']) for cp in cps])
        out['evdev'] = max(out['evdev'], np.max(np.abs(evb - evo) / np.maximum(1.0, np.abs(evo))))
        out['n'] += 1
    out['t'] = time.perf_counter() - t0
    if verbose:
        print('  %-10s brute-force null space (A^T A eigen-decomposition, no PRV, no recursion) for %d irreps (max nW %d): dim null != RS multiplicity for %d irreps; smallest non-null eigenvalue %.2e;'
              ' max distance of our copies from the null space %.1e; smallest singular value of the coordinates (span = m) %.6f; max |ev(V on null space) - ev(ours)| %.1e  (%.1fs)' %
              (str(list(lam)), out['n'], out['maxnW'], out['dimnull_bad'], out['gap'], out['nullnorm'], out['span_min'], out['evdev'], out['t']), flush=True)
    return out


# ======================================================================================== CLI
def main(argv):
    cmd = argv[1] if len(argv) > 1 else 'all'
    small = [(2, 1, 0), (3, 1, 0), (3, 2, 0), (4, 2, 0), (6, 3, 0), (8, 4, 0)]
    big = (14, 7, 0)
    if cmd in ('table1', 'all'):
        print('== R1: dim M_mu (PRV) vs Racah-Speiser, every irrep (and every dominant mu = Lam + w)')
        table1(small + [big])
        print('   extra reps with zero-multiplicity candidates (dominant Lam + w that are NOT irreps of V (x) V; there dim M_mu must be 0):')
        table1([(3, 0, 0), (4, 1, 0), (5, 2, 0), (5, 5, 0), (6, 2, 0), (7, 3, 0), (9, 1, 0), (10, 5, 0)])
    if cmd in ('hw', 'all'):
        print('== R2/R3: reconstructed copies: HW residual, weights, 36-term V vs V_red; exact twin (small reps)')
        test_hw(small, big)
        print('   named irreps of [14,7,0] x [14,7,0]:')
        show_irreps(big, [(14, 14, 14), (15, 14, 13), (28, 14, 0), (21, 14, 7), (17, 13, 12), (16, 14, 12), (17, 14, 11), (20, 12, 10), (18, 12, 12), (24, 12, 6), (27, 14, 1), (22, 13, 7)])
    if cmd in ('neg', 'all'):
        print('== negative tests of the exact consistency conditions')
        neg_tests([(2, 1, 0), (3, 1, 0), (3, 2, 0), (4, 2, 0)])
    if cmd in ('w4', 'all'):
        print('== R4: validation against W4 exact files (all branches)')
        tot = 0
        for lam in W4_REPS:
            r = validate_w4(lam, verbose=False)
            tot += r['nbranch']
            print('  %-10s branches %2d  mult mismatches %d  W4 HW kets %2d (slice not in M_mu: %d, exact reconstruction != W4 ket: %d)  LAURMAT A1 != V_red + eta1/2: %d  cas2 formula mismatches %d  eta1 formula mismatches %d  '
                  'copy-order mismatches %d  max |ev - (alpha1 - eta1/2)| %s  max slice dev %s  max vec dev (float) %.2e  (50 digits) %s  max G-norm^2 dev %s  exact HW violations %d' %
                  (str(list(lam)), r['nbranch'], r['mult_mismatch'], r['n_hwket'], r['hwket_not_in_prv'], r['hwket_mismatch'], r['laurmat_mismatch'], r['cas2_mismatch'], r['eta_mismatch'], r['order_mismatch'],
                   mp.nstr(r['ev_dev'], 3), mp.nstr(r['slice_dev'], 3), r['vec_dev_float'], mp.nstr(r['vec_dev_mp'], 3), mp.nstr(r['gnorm_dev'], 3), r['hw_exact_viol']), flush=True)
            if cmd == 'w4':
                for rec in r['records']:
                    print('      mu %s copy %d  D0 %3d  alpha_1 %s  eta_1 %s  ev %s  |ev-(alpha_1-eta_1/2)| %s  slice scale (W4/mine) %s' %
                          (rec['mu'], rec['copy'], rec['D0'], mp.nstr(rec['alpha1'], 12), mp.nstr(rec['eta1'], 6), mp.nstr(rec['ev'], 12), mp.nstr(rec['ev_dev'], 3), mp.nstr(rec['ratio_slice'], 8)))
        print('  total branches validated: %d' % tot)
    if cmd in ('bench', 'all'):
        amax = int(argv[2]) if len(argv) > 2 and cmd == 'bench' else 14
        amin = int(argv[3]) if len(argv) > 3 and cmd == 'bench' else 1
        print('== R5: benchmark Lam = [2a,a,0] (times in s; min of 5 cold runs; t_a = gen+rs+prv+Vred+eig, t_b = t_a + setup + recon)')
        bench(amax, amin=amin)
    if cmd in ('ref', 'all'):
        args = argv[2:] if cmd == 'ref' else ['all']
        todo = [(4, 2), (6, 3), (8, 4), (14, 7)] if args[0] == 'all' else [(int(args[0]), int(args[1]))]
        print('== R6: reference files')
        for (l1, l2) in todo:
            write_ref((l1, l2, 0), vec=(l1 <= 8) or (len(args) > 2 and args[2] == 'vec'))
    if cmd in ('verify', 'all'):
        print('== independent re-verification of the written reference files (parsed back from the text)')
        for name in ('420', '630', '840', '1470'):
            verify_ref(os.path.join(REFDIR, 'cg-slices-%s.txt' % name))
    if cmd in ('tridiag', 'all'):
        print('== structure: bandwidth of V_red in the GT basis of the weight space')
        tridiag_check([(4, 2, 0), (8, 4, 0), (14, 7, 0), (28, 14, 0), (10, 3, 0), (7, 7, 0), (6, 1, 0)])
    if cmd in ('brute', 'all'):
        print('== independent brute-force null space check (Munich-like)')
        for lam in ([(4, 2, 0), (6, 3, 0), (8, 4, 0)] + ([(14, 7, 0)] if True else [])):
            bruteforce_check(lam)
    if cmd in ('exact', 'all'):
        exact_all(int(argv[2]) if (cmd == 'exact' and len(argv) > 2) else 7)
    if cmd in ('closed', 'all'):
        print('== closed forms: singlet and top-irrep eigenvalues of V_red for [2a,a,0]')
        closed_forms(10)
    if cmd in ('rounded', 'all'):
        print('== float64 recursion vs exact recursion rounded to double (a = 14)')
        rounded_vs_float(14)
    if cmd in ('eigcheck', 'all'):
        eig_check((14, 7, 0))


if __name__ == '__main__':
    main(sys.argv)
