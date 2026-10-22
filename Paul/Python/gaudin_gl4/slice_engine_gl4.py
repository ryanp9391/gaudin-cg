"""slice_engine_gl4.py -- WP-G1: the PRV/Zhelobenko slice engine for gl(N) (tested at scale for N = 3, 4): highest-weight (HW) vectors of V_lam (x) V_lamp.

Maths (SPEC-G1-slice-engine-gl4.md; the gl(3) engine is Paul/Python/gaudin_sl3/slice_engine.py):
  Psi = sum Psi[nu,q] e_nu (x) e_q,  site 1 = V_lam (nu, slow), site 2 = V_lamp (q, fast); the SLICE of a vector of weight mu is psi_top = Psi[., q_top]
  (q_top = HW state of V_lamp = last GT index), a vector of the weight space V_lam[mu - lamp].
  Theorem (Kostant 1959 / PRV 1967):  Psi -> psi_top is an isomorphism of the HW vectors of weight mu onto
        M_mu = { v in V_lam[mu - lamp] : E_{i,i+1}^(lamp_i - lamp_{i+1} + 1) v = 0 , i = 1..N-1 },    dim M_mu = multiplicity of V_mu in V_lam (x) V_lamp.
  Level recursion (inverse map): label the weights of site 2 by a in Z>=0^(N-1), wt(q) = lamp - sum_i a_i alpha_i, Q(a) = GT states of that weight, X_a = (psi_q)_{q in Q(a)}
  with psi_q in V_lam[mu - wt(q)].  The HW conditions (E_i (x) 1 + 1 (x) E_i) Psi = 0 read at a site-2 component r in Q(a - e_i):
        E_i^(1) psi_r + sum_{q in Q(a)} (E_i^(2))_{rq} psi_q = 0       =>   A_a X_a = - [ E_i^(1) X_{a - e_i} ]_i   (stacked over the i with a - e_i a weight),
  A_a = [ (E_i^(2))_{Q(a-e_i), Q(a)} ]_i has full column rank, is independent of mu and of the slice: ONE left inverse per (lamp, a) (SVD pseudo-inverse in the
  orthonormal basis).  The remaining rows are consistency conditions T.RHS = 0 (exact twin), and the "dangling" rows E_i psi_r = 0 for r in Q(b) with b + e_i not a
  weight of site 2 (the PRV string conditions) -- the exact twin checks BOTH kinds (the T-conditions alone do not imply the dangling ones: N = 2 has no T-conditions).

Public API (lam, lamp, mu = tuples of ints; everything exact unless it says float)
  prv_space(lam, lamp, mu) -> PRV       exact basis of M_mu (.W .basis .free .m .w .vec(k) .coords(v))
  Recursion.of(lam, lamp)               float level recursion: left inverses once per (lam, lamp);  hw_from_slice(lam, lamp, mu, v) -> HW (COO, plain GT basis)
  ExactRecursion / hw_from_slice_exact(lam, lamp, mu, v)   exact twin (Fractions): (dict {(nu,q): Fraction}, #conditions verified [, violations])
  hw_from_slice_rounded(lam, lamp, mu, v)   exact recursion rounded to double (the accurate path)
  hw_residual(hw, extra=...)  exact_hw_check(lam, lamp, psi)  bruteforce_check(lam, lamp)  (Munich-style null space)
  table_prv / test_hw / neg_tests / regression_n3 / bench / families / write_ref / verify_ref
CLI:  python -X utf8 slice_engine_gl4.py <stage>   stages: tests | prv | hw | exact | neg | brute | bench [Smax [ext]] | families | ref [all | l1 l2 l3 l4 [l1p l2p l3p l4p]] | verify | munich | munichcmp S m1 m2 m3 m4 file | n5 | extras | all
"""
import os, sys, time, itertools, math, re, random
from fractions import Fraction as Fr
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla
import glN_reps as gr
from glN_reps import GLNRep, dij, wadd, weight_mults, racah_speiser, candidates, weyl_dim, casimir2, nW_of
from glN_reps import munich_patterns, munich_lowering_sq, munich_check, munich_label_to_ours      # the GT-basis side of the Munich bridge lives in glN_reps

REFDIR = os.path.normpath(os.path.join(HERE, '..', '..', 'Mathematica', 'Gaudin', 'ref'))
MUNICH = os.path.normpath(os.path.join(HERE, '..', '..', 'Cpp', 'AlexCGC', 'results', 'hw_timings.txt'))


def _t(x):
    return tuple(int(v) for v in x)


# ======================================================================================== 0. site-2 weight labels
def a_shift(a):
    """sum_i a_i alpha_i as an N-tuple (N = len(a) + 1), alpha_i = e_i - e_{i+1}"""
    N = len(a) + 1
    return tuple((a[k] if k < N - 1 else 0) - (a[k - 1] if k > 0 else 0) for k in range(N))


def site2_labels(R2):
    """dict a -> array of site-2 GT indices Q(a):  wt(q) = lamp - sum_i a_i alpha_i,  a_i = sum_{k<=i} (lamp_k - wt_k)"""
    lamp, N = R2.lam, R2.N
    ab = {}
    for w, idxs in R2.wsp.items():
        a, s = [], 0
        for k in range(N - 1):
            s += lamp[k] - w[k]
            a.append(s)
        ab[tuple(a)] = idxs
    return ab


def _prev(a, i):
    """a - e_i (i is 1-based) or None if a_i = 0"""
    if a[i - 1] < 1:
        return None
    return a[:i - 1] + (a[i - 1] - 1,) + a[i:]


# ======================================================================================== 1. PRV space (exact)
def _rref_nullspace(rows, ncols):
    """exact nullspace of the matrix 'rows' (list of lists of Fractions, ncols columns): basis from the RREF, one vector per free column f (entry 1 at f, zero at the
    other free columns).  Returns (basis, free), basis[k] a list of ncols Fractions."""
    M = [list(r) for r in rows if any(x != 0 for x in r)]
    piv = []
    r = 0
    for c in range(ncols):
        if r == len(M):
            break
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
    ps = set(piv)
    free = [c for c in range(ncols) if c not in ps]
    basis = []
    for f in free:
        v = [Fr(0)] * ncols
        v[f] = Fr(1)
        for k, pc in enumerate(piv):
            v[pc] = -M[k][f]
        basis.append(v)
    return basis, free


class PRV:
    """exact basis of M_mu inside W = V_lam[mu - lamp]:  basis[k] = list of Fractions aligned with W (global GT indices of site 1), the k-th vector is 1 at the k-th
    free position and 0 at the other free positions (the coordinates of any v in M_mu are its free entries)"""

    def __init__(self, lam, lamp, mu, W, basis, free):
        self.lam, self.lamp, self.mu, self.W, self.basis, self.free = lam, lamp, mu, W, basis, free
        self.m = len(basis)
        self.w = len(W)

    def vec(self, k):
        return {int(c): x for c, x in zip(self.W, self.basis[k]) if x != 0}

    def coords(self, v):
        """coordinates of a vector (list aligned with W) in the basis; raises ValueError if it is not in M_mu"""
        co = [v[f] for f in self.free]
        rec = [Fr(0)] * self.w
        for k, ck in enumerate(co):
            if ck != 0:
                for i, x in enumerate(self.basis[k]):
                    rec[i] += ck * x
        if rec != list(v):
            raise ValueError('vector not in M_mu')
        return co

    def float_basis(self):
        """(w, m) float array of the PRV basis vectors (plain GT basis)"""
        return np.array([[float(x) for x in b] for b in self.basis], dtype=float).reshape(self.m, self.w).T


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


def prv_space(lam, lamp, mu):
    """exact basis of M_mu = { v in V_lam[mu - lamp] : E_{i,i+1}^(lamp_i - lamp_{i+1} + 1) v = 0, i = 1..N-1 }  (dominant mu).  A constraint whose image weight is not
    a weight of V_lam is void."""
    lam, lamp, mu = _t(lam), _t(lamp), _t(mu)
    R = GLNRep.of(lam)
    N = R.N
    w0 = tuple(mu[k] - lamp[k] for k in range(N))
    W = R.wsp.get(w0)
    if W is None:
        return PRV(lam, lamp, mu, np.zeros(0, dtype=np.int64), [], [])
    rows = []
    for i in range(1, N):
        k = lamp[i - 1] - lamp[i] + 1
        if wadd(w0, dij(N, i, i + 1), k) not in R.wsp:
            continue                                         # E_i^k maps into a non-weight: constraint void
        imgs = power_images(R, (i, i + 1), k, W)
        byrow = defaultdict(lambda: [Fr(0)] * len(W))
        for a, vec in enumerate(imgs):
            for r, x in vec.items():
                byrow[r][a] = x
        rows.extend(byrow.values())
    basis, free = _rref_nullspace(rows, len(W))
    return PRV(lam, lamp, mu, W, basis, free)


# ======================================================================================== 2. the float level recursion (orthonormal GT basis)
class HW:
    """full HW vector(s) of weight mu: COO data in the PLAIN GT basis.  rows = site-1 state nu, cols = site-2 state q (global GT indices), vals[:, k] = k-th vector.
    Flat two-site index = nu * dim2 + q.  vals_hat = the same in the orthonormal basis (what the recursion computed)."""

    def __init__(self, lam, lamp, mu, rows, cols, vals_hat, sg1, sg2):
        self.lam, self.lamp, self.mu = lam, lamp, mu
        self.rows, self.cols, self.vals_hat = rows, cols, vals_hat
        self.sg1, self.sg2 = sg1, sg2
        self.dim1, self.dim2 = len(sg1), len(sg2)
        self.m = vals_hat.shape[1]

    @property
    def vals(self):
        return self.vals_hat / (self.sg1[self.rows] * self.sg2[self.cols])[:, None]

    @property
    def nW(self):
        return len(self.rows)

    def flat(self):
        return self.rows * self.dim2 + self.cols

    def hat_matrix(self, k):
        return sps.csr_matrix((self.vals_hat[:, k], (self.rows, self.cols)), shape=(self.dim1, self.dim2))

    def dense(self, k=0):
        X = np.zeros((self.dim1, self.dim2))
        X[self.rows, self.cols] = self.vals[:, k]
        return X


def _raising_blocks(R):
    """dense blocks of the SIMPLE raising operators in the orthonormal basis, keyed (i, row weight): E_i : wsp[Wr - alpha_i] -> wsp[Wr]   (cached on R)"""
    if getattr(R, '_B', None) is None:
        Eh = R.float_ops()
        N = R.N
        posw = np.zeros(R.dim, dtype=np.int64)
        for w, idxs in R.wsp.items():
            posw[idxs] = np.arange(len(idxs))
        B = {}
        for i in range(1, N):
            coo = Eh[(i, i + 1)].tocoo()
            d = dij(N, i, i + 1)
            for r, c, v in zip(coo.row, coo.col, coo.data):
                Wr = R.wt[r]
                blk = B.get((i, Wr))
                if blk is None:
                    blk = B[(i, Wr)] = np.zeros((len(R.wsp[Wr]), len(R.wsp[wadd(Wr, d, -1)])))
                blk[posw[r], posw[c]] = v
        R._B = B
    return R._B


class Recursion:
    """float level recursion for HW vectors of V_lam (x) V_lamp given the slice.  One left inverse of the stacked site-2 raising blocks A_a per site-2 weight a != 0
    (SVD pseudo-inverse, orthonormal basis), independent of mu and of the slice."""
    _cache = {}

    @classmethod
    def of(cls, lam, lamp=None):
        lam = _t(lam)
        lamp = lam if lamp is None else _t(lamp)
        key = (lam, lamp)
        if key not in cls._cache:
            cls._cache[key] = cls(lam, lamp)
        return cls._cache[key]

    def __init__(self, lam, lamp):
        t0 = time.perf_counter()
        self.R1, self.R2 = R1, R2 = GLNRep.of(lam), GLNRep.of(lamp)
        self.lam, self.lamp, self.N = lam, lamp, R1.N
        N = self.N
        self.sg1, self.sg2 = R1.sqrt_metric(), R2.sqrt_metric()
        self.B1 = _raising_blocks(R1)
        self.B2 = _raising_blocks(R2)
        self.ab2 = site2_labels(R2)
        self.order = sorted(self.ab2)
        self.shift = {a: a_shift(a) for a in self.order}
        self.wt2 = {a: wadd(lamp, self.shift[a], -1) for a in self.order}
        self.zero = (0,) * (N - 1)
        self.plan = {}
        self.smin = float('inf')
        for a in self.order:
            if a == self.zero:
                continue
            parts = []
            for i in range(1, N):
                pv = _prev(a, i)
                if pv is not None and pv in self.ab2:
                    blk = self.B2.get((i, self.wt2[pv]))
                    if blk is not None:
                        parts.append((i, pv, blk))
            nq = len(self.ab2[a])
            A = np.vstack([p[2] for p in parts])
            U, s, Vt = np.linalg.svd(A, full_matrices=False)
            if len(s) < nq or s[nq - 1] < 1e-10 * s[0]:
                raise ArithmeticError('site-2 raising block %s is not of full column rank (s = %s)' % (a, s))
            self.smin = min(self.smin, s[nq - 1] / s[0])
            AL = (Vt.T / s) @ U.T                                      # |Q| x r: left inverse (pseudo-inverse)
            off = 0
            lst = []
            for (i, pv, blk) in parts:
                r_i = blk.shape[0]
                lst.append((i, pv, np.ascontiguousarray(AL[:, off:off + r_i].T)))    # (|Q(prev)|, |Q(a)|)
                off += r_i
            self.plan[a] = lst
        self.t_setup = time.perf_counter() - t0

    def reconstruct(self, mu, vhat):
        """vhat: (w, m) orthonormal-basis components of m slices on W = V_lam[mu - lamp] (ascending global index).  Returns (rows, cols, vals_hat (nW, m))."""
        R1 = self.R1
        N = self.N
        w0 = tuple(mu[k] - self.lamp[k] for k in range(N))
        m = vhat.shape[1]
        wsp1, B1 = R1.wsp, self.B1
        Y = {}
        out_rows, out_cols, out_vals = [], [], []
        for a in self.order:
            sh = self.shift[a]
            W = tuple(w0[k] + sh[k] for k in range(N))
            S = wsp1.get(W)
            if S is None:
                continue
            Q = self.ab2[a]
            if a == self.zero:
                y = np.ascontiguousarray(vhat.T[:, :, None])               # (m, |S|, 1)
            else:
                y = None
                for (i, pv, ALi) in self.plan[a]:
                    yp = Y.get(pv)
                    if yp is None:
                        continue
                    T = B1.get((i, W))
                    if T is None:
                        continue
                    term = (T @ yp) @ ALi
                    y = -term if y is None else y - term
                if y is None:
                    y = np.zeros((m, len(S), len(Q)))
            Y[a] = y
            out_rows.append(np.repeat(S, len(Q)))
            out_cols.append(np.tile(Q, len(S)))
            out_vals.append(y.reshape(m, -1).T)
        return np.concatenate(out_rows), np.concatenate(out_cols), np.concatenate(out_vals, axis=0)


def hw_from_slice(lam, lamp, mu, v, check=None):
    """full HW vector(s) of weight mu in V_lam (x) V_lamp from the slice(s) v (floats, plain GT basis components on the weight space V_lam[mu - lamp] in ascending
    global GT order; shape (w,) or (w,k)).  v must lie in M_mu (the float recursion cannot see an inconsistent slice; the exact twin can).  check=tol raises
    ArithmeticError if the relative HW residual exceeds tol."""
    lam, lamp, mu = _t(lam), _t(lamp), _t(mu)
    rec = Recursion.of(lam, lamp)
    N = rec.N
    w0 = tuple(mu[k] - lamp[k] for k in range(N))
    W = rec.R1.wsp[w0]
    v = np.asarray(v, dtype=float)
    if v.ndim == 1:
        v = v[:, None]
    vhat = v * rec.sg1[W][:, None]
    rows, cols, vals_hat = rec.reconstruct(mu, vhat)
    hw = HW(lam, lamp, mu, rows, cols, vals_hat, rec.sg1, rec.sg2)
    if check is not None:
        r, wd = hw_residual(hw)
        if r > check or wd != 0:
            raise ArithmeticError('HW residual %.2e (weight deviation %d) exceeds %.1e' % (r, wd, check))
    return hw


# ======================================================================================== 3. the exact twin
def exact_left_inverse(A, nq):
    """A: list of r rows of nq Fractions with rank nq.  Returns (L, T): L (nq x r) with L A = 1 and T ((r-nq) x r) with T A = 0 (the consistency conditions)."""
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
    """the same level recursion in exact rational arithmetic (plain GT basis).  Verifies, block by block, (T) the exact consistency conditions T.RHS = 0 of the
    overdetermined stacked system and (D) the dangling conditions E_i psi_r = 0 at the ends of the site-2 strings; together they ARE the full HW conditions
    (every row (i, r) of E_i^tot Psi = 0 is either used to solve, or a T-row, or a D-row)."""
    _cache = {}

    @classmethod
    def of(cls, lam, lamp=None):
        lam = _t(lam)
        lamp = lam if lamp is None else _t(lamp)
        key = (lam, lamp)
        if key not in cls._cache:
            cls._cache[key] = cls(lam, lamp)
        return cls._cache[key]

    def __init__(self, lam, lamp):
        self.R1, self.R2 = R1, R2 = GLNRep.of(lam), GLNRep.of(lamp)
        self.lam, self.lamp, self.N = lam, lamp, R1.N
        N = self.N
        self.ab2 = site2_labels(R2)
        self.order = sorted(self.ab2)
        self.shift = {a: a_shift(a) for a in self.order}
        self.wt2 = {a: wadd(lamp, self.shift[a], -1) for a in self.order}
        self.zero = (0,) * (N - 1)
        self.L, self.T, self.prevs = {}, {}, {}
        self._blk = {}
        for a in self.order:
            if a == self.zero:
                continue
            rows, prevs = [], []
            for i in range(1, N):
                pv = _prev(a, i)
                if pv is not None and pv in self.ab2:
                    M, rr, cc = R2.block_exact((i, i + 1), self.wt2[pv])
                    rows += M
                    prevs.append((i, pv, len(M)))
            L, T = exact_left_inverse(rows, len(self.ab2[a]))
            self.L[a] = _obj(L)
            self.T[a] = _obj(T, (0, len(rows))) if T else None
            self.prevs[a] = prevs

    def _block1(self, i, W):
        key = (i, W)
        b = self._blk.get(key)
        if b is None:
            M, rr, cc = self.R1.block_exact((i, i + 1), W)
            b = self._blk[key] = (_obj(M) if M else None)
        return b

    def reconstruct(self, mu, v, report=False):
        """v: list of Fractions aligned with W = V_lam[mu - lamp].  Returns (dict {(nu, q): Fraction} (all entries of the weight-mu subspace, zeros kept),
        number of exact conditions verified, [#violated T-conditions, #violated dangling conditions]); raises ValueError at the first violation unless report=True."""
        R1 = self.R1
        N = self.N
        w0 = tuple(mu[k] - self.lamp[k] for k in range(N))
        Y, out = {}, {}
        nchk = 0
        viol = [0, 0]
        for a in self.order:
            W = wadd(w0, self.shift[a])
            S = R1.wsp.get(W)
            if S is None:
                continue
            Q = self.ab2[a]
            if a == self.zero:
                assert len(v) == len(S)
                y = np.empty((len(S), 1), dtype=object)
                for i in range(len(S)):
                    y[i, 0] = Fr(v[i])
            else:
                parts = []
                for (i, pv, r_i) in self.prevs[a]:
                    yp = Y.get(pv)
                    M = self._block1(i, W)
                    if yp is None or M is None:
                        parts.append(np.full((len(S), r_i), Fr(0), dtype=object))
                    else:
                        parts.append(-(M @ yp))
                Rr = parts[0] if len(parts) == 1 else np.concatenate(parts, axis=1)       # (|S|, r)
                y = Rr @ self.L[a].T                                                       # (|S|, |Q|)
                T = self.T[a]
                if T is not None:
                    cond = Rr @ T.T
                    nchk += cond.size
                    nz = sum(1 for x in cond.ravel() if x != 0)
                    if nz:
                        viol[0] += nz
                        if not report:
                            raise ValueError('exact consistency condition (T) violated at site-2 weight %s: the slice is not in M_mu' % (a,))
            Y[a] = y
            for i in range(1, N):                                      # dangling conditions: b + e_i is not a weight of site 2  =>  E_i psi_r = 0
                nxt = a[:i - 1] + (a[i - 1] + 1,) + a[i:]
                if nxt not in self.ab2:
                    M = self._block1(i, wadd(W, dij(N, i, i + 1)))
                    if M is not None:
                        cond = M @ y
                        nchk += cond.size
                        nz = sum(1 for x in cond.ravel() if x != 0)
                        if nz:
                            viol[1] += nz
                            if not report:
                                raise ValueError('exact dangling (PRV string) condition violated at site-2 weight %s, i = %d: the slice is not in M_mu' % (a, i))
            for ii, nu in enumerate(S):
                for jj, q in enumerate(Q):
                    out[(int(nu), int(q))] = y[ii, jj]
        return out, nchk, viol


def hw_from_slice_exact(lam, lamp, mu, v, report=False):
    """exact HW vector (dict {(nu,q): Fraction}) from an exact slice v (list of Fractions aligned with V_lam[mu - lamp]); also returns the number of exact conditions
    verified (T + dangling) and the violation counts [T, dangling]."""
    return ExactRecursion.of(lam, lamp).reconstruct(_t(mu), v, report=report)


def hw_from_slice_rounded(lam, lamp, mu, v):
    """exact recursion followed by rounding to double: an HW object with ONE vector, for an exact slice v (list of Fractions aligned with the weight space).
    The float64 recursion loses accuracy with the size of the representation; this path does not."""
    lam, lamp, mu = _t(lam), _t(lamp), _t(mu)
    R1, R2 = GLNRep.of(lam), GLNRep.of(lamp)
    psi, _, _ = hw_from_slice_exact(lam, lamp, mu, v)
    sg1, sg2 = R1.sqrt_metric(), R2.sqrt_metric()
    rows = np.array([k[0] for k in psi], dtype=np.int64)
    cols = np.array([k[1] for k in psi], dtype=np.int64)
    vals = np.array([float(x) for x in psi.values()])
    return HW(lam, lamp, mu, rows, cols, (vals * sg1[rows] * sg2[cols])[:, None], sg1, sg2)


def exact_hw_check(lam, lamp, psi, pairs=None):
    """apply the total E_ij (all i < j, or the given pairs) exactly to a two-site vector {(nu,q): Fraction}; returns the number of nonzero components left (0 = HW)"""
    R1, R2 = GLNRep.of(lam), GLNRep.of(lamp)
    N = R1.N
    pairs = pairs or [(i, j) for i in range(1, N + 1) for j in range(i + 1, N + 1)]
    bad = 0
    for ij in pairs:
        E1, E2 = R1.E[ij], R2.E[ij]
        res = {}
        for (nu, q), x in psi.items():
            if x == 0:
                continue
            for r, e in E1[nu].items():
                res[(r, q)] = res.get((r, q), 0) + x * e
            for r, e in E2[q].items():
                res[(nu, r)] = res.get((nu, r), 0) + x * e
        bad += sum(1 for x in res.values() if x != 0)
    return bad


def exact_gnorm2(lam, lamp, psi):
    """exact G-norm^2 = sum g1[nu] g2[q] psi[nu,q]^2  (G = GMET1 (x) GMET2)"""
    g1, g2 = GLNRep.of(lam).g, GLNRep.of(lamp).g
    return sum(g1[nu] * g2[q] * x * x for (nu, q), x in psi.items())


# ======================================================================================== 4. residuals
def hw_residual(hw, extra=False, per_vector=False):
    """max over vectors of ||E_{i,i+1}^tot Psi|| / ||Psi|| (invariant norm = Frobenius norm in the orthonormal basis) for the simple raising operators (extra=True:
    for ALL E_ij with i < j), and the max deviation of the total weight from mu over all entries.  per_vector=True returns the list of residuals instead of the max."""
    R1, R2 = GLNRep.of(hw.lam), GLNRep.of(hw.lamp)
    N = R1.N
    Eh1, Eh2 = R1.float_ops(), R2.float_ops()
    wt1, wt2 = np.array(R1.wt), np.array(R2.wt)
    wdev = int(np.abs(wt1[hw.rows] + wt2[hw.cols] - np.array(hw.mu)).max())
    ijs = [(i, i + 1) for i in range(1, N)]
    if extra:
        ijs = [(i, j) for i in range(1, N + 1) for j in range(i + 1, N + 1)]
    res = []
    for k in range(hw.m):
        X = hw.hat_matrix(k)
        nx = spla.norm(X)
        worst = 0.0
        for ij in ijs:
            Rm = Eh1[ij] @ X + X @ Eh2[ij].T
            worst = max(worst, (spla.norm(Rm) / nx) if Rm.nnz else 0.0)
        res.append(worst)
    return (res if per_vector else max(res)), wdev


# ======================================================================================== 5. tables, tests
def _clear_caches():
    GLNRep._cache.clear()
    gr._WM_CACHE.clear()
    Recursion._cache.clear()
    ExactRecursion._cache.clear()


def fmtl(lam):
    return '[' + ','.join(str(x) for x in lam) + ']'


def prod_name(lam, lamp):
    return fmtl(lam) + ('^2' if tuple(lam) == tuple(lamp) else 'x' + fmtl(lamp))


PRODUCTS = [((1, 0, 0, 0),) * 2, ((1, 1, 0, 0),) * 2, ((2, 1, 0, 0),) * 2, ((2, 1, 1, 0),) * 2, ((2, 2, 0, 0),) * 2, ((3, 2, 1, 0),) * 2, ((2, 2, 1, 0),) * 2,
            ((3, 3, 0, 0),) * 2,
            ((2, 1, 0, 0), (1, 1, 0, 0)), ((1, 1, 0, 0), (2, 1, 0, 0)), ((3, 1, 0, 0), (2, 2, 0, 0)), ((2, 2, 0, 0), (3, 1, 0, 0)), ((2, 2, 1, 0), (2, 1, 1, 0)),
            ((2, 1, 1, 0), (2, 2, 1, 0)), ((3, 2, 1, 0), (2, 1, 0, 0)), ((2, 1, 0, 0), (3, 2, 1, 0))]
RECT = [((S, S, 0, 0),) * 2 for S in range(1, 11)]
EXTRA = [((4, 2, 0, 0),) * 2, ((6, 4, 2, 0),) * 2]                     # the multiplicity families of G3: max m = 3 and 12
BIG = [((9, 6, 3, 0),) * 2]                                         # dim 4096, 1592 vectors, m up to 27: PRV and float recursion only
ALLPRODUCTS = list(dict.fromkeys(PRODUCTS + RECT + EXTRA))      # [S,S,0,0]^2 for S = 1,2,3 are also in PRODUCTS: no duplicates


def table_prv(products, verbose=True):
    """G1.2: dim M_mu vs Racah-Speiser for EVERY dominant candidate mu (mu = lam + w, w in wt(V_lamp), and mu = lamp + w1, w1 in wt(V_lam): the union), including those
    of multiplicity 0 (M_mu must be {0}); returns rows"""
    rows = []
    for lam, lamp in products:
        t0 = time.perf_counter()
        R1, R2 = GLNRep.of(lam), GLNRep.of(lamp)
        rs = racah_speiser(lam, lamp)
        sumdim = sum(m * weyl_dim(mu) for mu, m in rs.items())
        cand = candidates(lam, lamp, both=True)
        mism = 0
        for mu in cand:
            p = prv_space(lam, lamp, mu)
            if p.m != rs.get(mu, 0):
                mism += 1
        extra = [mu for mu in rs if mu not in set(cand)]
        maxw = max(len(R1.wsp.get(tuple(mu[k] - lamp[k] for k in range(R1.N)), ())) for mu in rs)
        row = dict(lam=lam, lamp=lamp, ncand=len(cand), nirrep=len(rs), summ=sum(rs.values()), maxm=max(rs.values()), maxw=maxw,
                   maxnW=max(nW_of(lam, lamp, mu) for mu in rs), dimsum_ok=(sumdim == R1.dim * R2.dim), mism=mism + len(extra), nzero=len(cand) - len(rs),
                   t=time.perf_counter() - t0)
        rows.append(row)
        if verbose:
            print('  %-26s dim %5d x %-5d candidates %4d  irreps %4d  sum m %5d  max m %3d  max w %3d  max nW %6d  sum m dimV = dim dim\': %s  mismatches %d  (%.2fs)' %
                  (prod_name(lam, lamp), R1.dim, R2.dim, row['ncand'], row['nirrep'], row['summ'], row['maxm'], row['maxw'], row['maxnW'], row['dimsum_ok'], row['mism'],
                   row['t']), flush=True)
    return rows


def regression_n3(verbose=True):
    """N = 3 regression: the PRV bases coincide with slice_engine.prv_space exactly (same W, basis, free) on [4,2,0]^2, [6,3,0]^2; [3,2,0]x[3,1,0]: multiplicities;
    the float recursion reproduces slice_engine.hw_from_slice entry for entry"""
    sys.path.insert(0, os.path.normpath(os.path.join(HERE, '..', 'gaudin_sl3')))
    import slice_engine as se
    out = []
    for lam in [(4, 2, 0), (6, 3, 0), (14, 7, 0)]:
        ncand = nbasis_bad = nvec = 0
        worst = 0.0
        for mu in candidates(lam, lam, both=True):
            p = prv_space(lam, lam, mu)
            q = se.prv_space(lam, mu)
            ncand += 1
            if not (list(p.W) == list(q.W) and p.basis == q.basis and p.free == q.free):
                nbasis_bad += 1
            if p.m > 0 and p.m == q.m:
                h1 = hw_from_slice(lam, lam, mu, p.float_basis())
                h2 = se.hw_from_slice(lam, mu, np.array([[float(x) for x in b] for b in q.basis]).T)
                assert np.array_equal(h1.rows, h2.rows) and np.array_equal(h1.cols, h2.cols)
                worst = max(worst, float(np.abs(h1.vals - h2.vals).max() / np.abs(h2.vals).max()))
                nvec += p.m
        out.append((lam, ncand, nbasis_bad, nvec, worst))
        if verbose:
            print('  %-14s N = 3 regression: %d dominant candidates, PRV (W, basis, free) different from slice_engine.prv_space: %d; float recursion vs slice_engine.hw_from_slice on %d vectors: '
                  'max relative deviation %.2e' % (prod_name(lam, lam), ncand, nbasis_bad, nvec, worst), flush=True)
    lam, lamp = (3, 2, 0), (3, 1, 0)
    rs = racah_speiser(lam, lamp)
    assert rs == se.racah_speiser(lam, lamp)
    bad = sum(1 for mu in candidates(lam, lamp, both=True) if prv_space(lam, lamp, mu).m != rs.get(mu, 0))
    bad += sum(1 for mu in candidates(lamp, lam, both=True) if prv_space(lamp, lam, mu).m != rs.get(mu, 0))
    out.append((lam, lamp, len(rs), bad))
    if verbose:
        print('  %-14s N = 3 lam != lam\': %d irreps (Racah-Speiser identical to slice_engine), dim M_mu mismatches in both orders: %d' % (prod_name(lam, lamp), len(rs), bad), flush=True)
    return out


def test_hw(products, verbose=True, nround=0):
    """G1.3 float path: reconstruct ALL m PRV basis vectors of every irrep, HW residual for ALL E_ij (i < j) and the total weight; sectors above 1e-12 are re-done with
    the exact recursion rounded to double"""
    out = []
    for lam, lamp in products:
        t0 = time.perf_counter()
        R1 = GLNRep.of(lam)
        rs = racah_speiser(lam, lamp)
        worst, wdev, nvec, nbig, worst_rounded, nround_done, nw_bad = -1.0, 0, 0, 0, 0.0, 0, 0
        worst_mu = None
        for mu in sorted(rs, reverse=True):
            p = prv_space(lam, lamp, mu)
            hw = hw_from_slice(lam, lamp, mu, p.float_basis())
            if hw.nW != nW_of(lam, lamp, mu):
                nw_bad += 1
            r, wd = hw_residual(hw, extra=True)
            if r > worst:
                worst, worst_mu = r, mu
            wdev = max(wdev, wd)
            nvec += hw.m
            if r > 1e-12:
                nbig += 1
                for k in range(p.m):
                    hr = hw_from_slice_rounded(lam, lamp, mu, p.basis[k])
                    rr, wr = hw_residual(hr, extra=True)
                    worst_rounded = max(worst_rounded, rr)
                    nround_done += 1
        row = dict(lam=lam, lamp=lamp, nirrep=len(rs), nvec=nvec, resid=worst, worst_mu=worst_mu, wdev=wdev, nbig=nbig, rounded=worst_rounded, nw_bad=nw_bad,
                   smin=Recursion.of(lam, lamp).smin, t=time.perf_counter() - t0)
        out.append(row)
        if verbose:
            print('  %-26s irreps %3d  vectors %4d  max HW residual (all E_ij, i<j) %.2e at mu = %s  weight deviation %d  nW mismatches %d  sectors > 1e-12: %d%s  min sv ratio of A_a %.3f  (%.1fs)' %
                  (prod_name(lam, lamp), len(rs), nvec, worst, worst_mu, wdev, nw_bad, nbig,
                   (' (rounded-exact path on %d vectors: %.1e)' % (nround_done, worst_rounded)) if nbig else '', row['smin'], row['t']), flush=True)
    return out


def hw_table_rect(S, verbose=True):
    """per-irrep table of [S,S,0,0]^2 (Munich's SU(4) labels): dim V_mu, nW, w, HW residual (float, all E_ij), HW residual of the rounded-exact path if > 1e-12"""
    lam = (S, S, 0, 0)
    rs = racah_speiser(lam, lam)
    rows = []
    for mu in sorted(rs, reverse=True):
        p = prv_space(lam, lam, mu)
        hw = hw_from_slice(lam, lam, mu, p.float_basis())
        r, wd = hw_residual(hw, extra=True)
        rows.append((mu, tuple(x - mu[3] for x in mu), weyl_dim(mu), hw.nW, p.w, r, wd))
    if verbose:
        print('  [%d,%d,0,0]^2: %d irreps (Munich SU(4) labels = ours minus mu_4)' % (S, S, len(rows)))
        print('   %-14s %-14s %7s %7s %3s %10s %3s' % ('mu (ours)', 'mu (Munich)', 'dim', 'nW', 'w', 'HW resid', 'wdev'))
        for (mu, mm, d, nW, w, r, wd) in rows:
            print('   %-14s %-14s %7d %7d %3d %10.2e %3d' % (mu, mm, d, nW, w, r, wd))
        print('   max HW residual %.2e, max weight deviation %d' % (max(r[5] for r in rows), max(r[6] for r in rows)), flush=True)
    return rows


def exact_tests(products, max_vec=None, verbose=True, label=None, sample=None, seed=3):
    """G1.3 exact twin: for the irreps of each product (all, or a sample of 'sample' irreps incl. the largest n_W) and every PRV basis vector: exact recursion (T +
    dangling conditions verified), exact HW check E_ij^tot Psi = 0 for ALL i < j.  Returns rows."""
    rng = random.Random(seed)
    out = []
    for lam, lamp in products:
        t0 = time.perf_counter()
        rs = racah_speiser(lam, lamp)
        mus = sorted(rs, reverse=True)
        if sample is not None and sample < len(mus):
            big = sorted(mus, key=lambda mu: -nW_of(lam, lamp, mu))[:max(1, sample // 2)]
            rest = [mu for mu in mus if mu not in big]
            mus = sorted(set(big + rng.sample(rest, sample - len(big))), reverse=True)
        nvec = nchk = viol = ent = 0
        for mu in mus:
            p = prv_space(lam, lamp, mu)
            for k in range(p.m):
                psi, nc, vi = hw_from_slice_exact(lam, lamp, mu, p.basis[k])
                nchk += nc
                nvec += 1
                ent += len(psi)
                viol += exact_hw_check(lam, lamp, psi)
        row = dict(lam=lam, lamp=lamp, nirrep=len(mus), nirrep_all=len(rs), nvec=nvec, nchk=nchk, viol=viol, ent=ent, t=time.perf_counter() - t0)
        out.append(row)
        if verbose:
            print('  %-26s irreps %3d of %3d  PRV basis vectors %4d  entries %7d  exact conditions verified (T + dangling) %8d  exact HW violations (all E_ij, i<j) %d  (%.1fs)' %
                  (prod_name(lam, lamp), row['nirrep'], row['nirrep_all'], nvec, ent, nchk, viol, row['t']), flush=True)
    return out


def neg_tests(products, nrand=20, seed=1, verbose=True):
    """G1.3 negative tests: random exact vectors of V_lam[mu - lamp] NOT in M_mu (m < w sectors) must be rejected by the exact conditions; random vectors inside M_mu accepted.
    Reports how many were rejected by the T-conditions and how many ONLY by the dangling (PRV string) conditions."""
    rng = random.Random(seed)
    res = []
    for lam, lamp in products:
        rs = racah_speiser(lam, lamp)
        ntot = nrej = nT = nDonly = nin = ninok = nsect = incons = 0
        minres = float('inf')                                # smallest float HW residual of an out-of-M_mu slice (the float recursion cannot see inconsistency, the residual check can)
        for mu in sorted(rs, reverse=True):
            p = prv_space(lam, lamp, mu)
            if p.w > p.m:
                nsect += 1
                for _ in range(nrand):
                    v = [Fr(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(p.w)]
                    try:
                        p.coords(v)
                        continue                                 # accidentally inside M_mu
                    except ValueError:
                        pass
                    ntot += 1
                    psi, nc, vi = hw_from_slice_exact(lam, lamp, mu, v, report=True)
                    minres = min(minres, hw_residual(hw_from_slice(lam, lamp, mu, np.array([float(x) for x in v])))[0])
                    rejected = bool(vi[0] or vi[1])
                    # cross-validation of the condition logic: (some exact condition violated)  <=>  (the reconstructed vector is NOT exactly HW)
                    if rejected != (exact_hw_check(lam, lamp, psi) > 0):
                        incons += 1
                    if rejected:
                        nrej += 1
                        if vi[0]:
                            nT += 1
                        else:
                            nDonly += 1
            for _ in range(3):
                c = [Fr(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(p.m)]
                v = [sum(c[k] * p.basis[k][i] for k in range(p.m)) for i in range(p.w)]
                nin += 1
                psi, nc, vi = hw_from_slice_exact(lam, lamp, mu, v, report=True)
                if not (vi[0] or vi[1]):
                    ninok += 1
                if (vi[0] or vi[1]) != (exact_hw_check(lam, lamp, psi) > 0):
                    incons += 1
        res.append((lam, lamp, nsect, ntot, nrej, nT, nDonly, nin, ninok, incons, minres))
        if verbose:
            print('  %-26s sectors with w > m: %3d  random slices outside M_mu: %4d, rejected: %4d (by a T-condition: %4d, only by a dangling PRV condition: %4d);  random slices inside M_mu: %3d, accepted: %3d;  '
                  '(conditions violated) != (vector not exactly HW): %d;  min float HW residual of a rejected slice %.1e' % (prod_name(lam, lamp), nsect, ntot, nrej, nT, nDonly, nin, ninok, incons, minres), flush=True)
    return res


# ======================================================================================== 6. Munich-style brute force null space
def bruteforce_check(lam, lamp, mus=None, verbose=True, tol_null=1e-9, maxnW=None):
    """INDEPENDENT of the PRV theorem and of the recursion: for each irrep mu build the constraint matrix of the stacked total simple raising operators
    E_{i,i+1}^tot = E (x) 1 + 1 (x) E on the weight-mu subspace of V_lam (x) V_lamp (n_W columns; orthonormal GT basis, scipy kron), find its null space by the
    eigen-decomposition of A^T A (what Munich/fastEig do), and compare: (i) dim null = Racah-Speiser m (and the first non-null eigenvalue), (ii) our reconstructed vectors
    (orthonormalised) lie in it and span it (largest principal angle = || (1 - N N^T) Q ||_2).  Returns the worst numbers."""
    lam, lamp = _t(lam), _t(lamp)
    R1, R2 = GLNRep.of(lam), GLNRep.of(lamp)
    N = R1.N
    n1, n2 = R1.dim, R2.dim
    Eh1, Eh2 = R1.float_ops(), R2.float_ops()
    rs = racah_speiser(lam, lamp)
    mus = sorted(rs, reverse=True) if mus is None else [tuple(m) for m in mus]
    I1, I2 = sps.identity(n1, format='csr'), sps.identity(n2, format='csr')
    Etot = [(sps.kron(Eh1[(i, i + 1)], I2) + sps.kron(I1, Eh2[(i, i + 1)])).tocsc() for i in range(1, N)]
    out = dict(lam=lam, lamp=lamp, n=0, dimnull_bad=0, gap=np.inf, nullnorm=0.0, maxangle=0.0, maxnW=0, t=0.0)
    t0 = time.perf_counter()
    for mu in mus:
        flat = []
        for w2, Q in R2.wsp.items():
            S = R1.wsp.get(tuple(mu[k] - w2[k] for k in range(N)))
            if S is None:
                continue
            for nu in S:
                flat.extend(int(nu) * n2 + int(q) for q in Q)
        flat = np.array(sorted(flat), dtype=np.int64)
        nW = len(flat)
        if maxnW is not None and nW > maxnW:
            continue
        out['maxnW'] = max(out['maxnW'], nW)
        Gm = None
        for E in Etot:
            A = E[:, flat]
            G1 = (A.T @ A)
            Gm = G1 if Gm is None else Gm + G1
        w_, U = np.linalg.eigh(Gm.toarray())
        m = rs.get(mu, 0)
        nnull = int(np.sum(w_ < tol_null))
        if nnull != m:
            out['dimnull_bad'] += 1
        out['gap'] = min(out['gap'], w_[nnull] if nnull < nW else np.inf)
        if m == 0:
            out['n'] += 1
            continue
        Nn = U[:, :m]
        p = prv_space(lam, lamp, mu)
        hw = hw_from_slice(lam, lamp, mu, p.float_basis())
        fl = hw.rows * n2 + hw.cols
        pos = np.searchsorted(flat, fl)
        assert np.array_equal(flat[pos], fl) and len(fl) == nW
        Psi = np.zeros((nW, m))
        Psi[pos, :] = hw.vals_hat
        Qm, _ = np.linalg.qr(Psi)
        resid = Qm - Nn @ (Nn.T @ Qm)
        out['nullnorm'] = max(out['nullnorm'], float(np.linalg.norm(resid, axis=0).max()))
        out['maxangle'] = max(out['maxangle'], float(np.linalg.norm(resid, 2)))
        out['n'] += 1
    out['t'] = time.perf_counter() - t0
    if verbose:
        print('  %-26s brute-force null space (A^T A eigen-decomposition, no PRV, no recursion) for %d irreps (max nW %d): dim null != RS multiplicity for %d irreps; '
              'smallest non-null eigenvalue %.3e; max distance of our vectors from the null space %.1e; largest principal angle (our span vs null space) %.1e  (%.1fs)' %
              (prod_name(lam, lamp), out['n'], out['maxnW'], out['dimnull_bad'], out['gap'], out['nullnorm'], out['maxangle'], out['t']), flush=True)
    return out


# ======================================================================================== 7. benchmark
def parse_munich(path=MUNICH):
    """{S: (n_irreps, total hw-solve s, singlet hw-solve s, {mu_munich: (dim, mult, t)})} from results/hw_timings.txt"""
    out = {}
    cur = None
    for line in open(path, encoding='utf8'):
        m = re.match(r'SU\(4\):\s+irreps mu in\s+\[(\d+),(\d+),0,0\] \(x\) \[(\d+),(\d+),0,0\]', line)
        if m:
            cur = out[int(m.group(1))] = dict(rows={}, n=None, total=None)
            continue
        m = re.match(r'\((\d+) (\d+) (\d+) (\d+)\)\s+(\d+)\s+(\d+)\s+([\d.]+) s', line)
        if m and cur is not None:
            cur['rows'][tuple(int(m.group(k)) for k in (1, 2, 3, 4))] = (int(m.group(5)), int(m.group(6)), float(m.group(7)))
        m = re.match(r'irreps mu:\s+(\d+)\s+total highest-weight-solve time:\s+([\d.]+) s', line)
        if m and cur is not None:
            cur['n'], cur['total'] = int(m.group(1)), float(m.group(2))
    return out


def bench_one(lam, lamp=None, reps=5):
    """G1.5 for the product lam (x) lamp (default lamp = lam):  (a) = exact GT generators (simple operators; the non-adjacent ones are lazy) + Racah-Speiser (incl. the GT
    weight multiplicities) + exact PRV spaces of all irreps  (no n_W-sized object is touched);  (b) = (a) + metric/orthonormal float matrices/left inverses + the full
    reconstruction of EVERY PRV basis vector of every irrep (batched per irrep).  min over `reps` cold runs (caches cleared, single process)."""
    lam = _t(lam)
    lamp = lam if lamp is None else _t(lamp)
    best = None
    for _ in range(reps):
        _clear_caches()
        T = {}
        t0 = time.perf_counter()
        R1 = GLNRep.of(lam)
        R2 = GLNRep.of(lamp)
        T['gen'] = time.perf_counter() - t0
        t = time.perf_counter()
        rs = racah_speiser(lam, lamp)
        mus = sorted(rs, reverse=True)
        T['rs'] = time.perf_counter() - t
        t = time.perf_counter()
        prvs = {mu: prv_space(lam, lamp, mu) for mu in mus}
        T['prv'] = time.perf_counter() - t
        T['a'] = time.perf_counter() - t0
        t = time.perf_counter()
        rec = Recursion.of(lam, lamp)
        T['setup_b'] = time.perf_counter() - t
        t = time.perf_counter()
        nent = 0
        for mu in mus:
            hw = hw_from_slice(lam, lamp, mu, prvs[mu].float_basis())
            nent += hw.nW * hw.m
        T['recon'] = time.perf_counter() - t
        T['b'] = T['a'] + T['setup_b'] + T['recon']
        T['nent'] = nent
        if best is None or T['b'] < best['b']:
            best = T
    # accuracy pass (not timed): worst HW residual over all irreps (simple raising operators + weights), and for small cases all E_ij
    worst, wdev = 0.0, 0
    worst_mu = None
    for mu in mus:
        hw = hw_from_slice(lam, lamp, mu, prvs[mu].float_basis())
        r, wd = hw_residual(hw, extra=(R1.dim <= 400))
        if r > worst:
            worst, worst_mu = r, mu
        wdev = max(wdev, wd)
    best.update(lam=lam, lamp=lamp, dim=R1.dim, nirrep=len(rs), ncopies=sum(rs.values()), maxm=max(rs.values()), maxw=max(p.w for p in prvs.values()),
                maxnW=max(nW_of(lam, lamp, mu) for mu in mus), resid=worst, worst_mu=worst_mu, wdev=wdev, smin=rec.smin, ta=best['a'], tb=best['b'])
    return best


def bench(Smin=1, Smax=10, reps=5, verbose=True, header=True):
    mun = parse_munich()
    rows = []
    if verbose and header:
        print('  %2s %5s %6s %6s %5s %7s | %7s %7s %7s | %8s | %7s %7s | %8s | %10s | %8s | %s' %
              ('S', 'dimV', 'irreps', 'copies', 'max w', 'max nW', 't_gen', 't_rs', 't_prv', 't_a', 't_setup', 't_recon', 't_b', 'entries', 'resid', 'Munich total hw-solve (i5-1345U) (x = Munich / t_b)'), flush=True)
    for S in range(Smin, Smax + 1):
        T = bench_one((S, S, 0, 0), reps=reps)
        T['munich'] = mun[S]['total'] if S in mun else None
        rows.append(T)
        if verbose:
            print('  %2d %5d %6d %6d %5d %7d | %7.4f %7.4f %7.4f | %8.4f | %7.4f %7.4f | %8.4f | %10d | %8.1e | %s   [min A_a sv ratio %.3f; worst resid at mu=%s]' %
                  (S, T['dim'], T['nirrep'], T['ncopies'], T['maxw'], T['maxnW'], T['gen'], T['rs'], T['prv'], T['ta'], T['setup_b'], T['recon'], T['tb'], T['nent'], T['resid'],
                   (('%.4g s  (x%.0f)' % (T['munich'], T['munich'] / T['tb'])) if T['munich'] is not None else 'n/a (not in hw_timings.txt)'), T['smin'], T['worst_mu']), flush=True)
    return rows


def rounded_vs_float(S, nworst=3, verbose=True):
    """the float64 recursion's HW residual grows with the size of the representation; for the nworst worst irreps of [S,S,0,0]^2 compare it with the exact recursion rounded to
    double (hw_from_slice_rounded): residual, time, and the relative difference of the two vectors"""
    lam = (S, S, 0, 0)
    res = []
    for mu in sorted(racah_speiser(lam, lam), reverse=True):
        p = prv_space(lam, lam, mu)
        hw = hw_from_slice(lam, lam, mu, p.float_basis())
        r, wd = hw_residual(hw)
        res.append((r, mu))
    res.sort(reverse=True)
    out = []
    for r, mu in res[:nworst]:
        p = prv_space(lam, lam, mu)
        t0 = time.perf_counter()
        t = time.perf_counter()
        hw = hw_from_slice(lam, lam, mu, p.float_basis())
        tf = time.perf_counter() - t
        t = time.perf_counter()
        hr = hw_from_slice_rounded(lam, lam, mu, p.basis[0])
        te = time.perf_counter() - t
        rr, _ = hw_residual(hr)
        dev = float(np.abs(hw.vals_hat[:, 0] - hr.vals_hat[:, 0]).max() / np.abs(hr.vals_hat[:, 0]).max())
        out.append((mu, hw.nW, r, rr, tf, te, dev))
        if verbose:
            print('  S = %2d  mu = %-16s nW %6d  HW residual: float64 recursion %.1e (%.3f s), exact recursion rounded to double %.1e (%.2f s); max rel. difference of the two vectors %.1e' %
                  (S, mu, hw.nW, r, tf, rr, te, dev), flush=True)
    return out


def time_singlet(S, reps=3):
    """warm timing of ONE irrep (the singlet, Munich's slowest row) after the per-lam setup: prv_space + hw_from_slice; and the cold single-irrep time (everything)"""
    lam = (S, S, 0, 0)
    mu = (S, S, S, S)
    cold = None
    for _ in range(reps):
        _clear_caches()
        t0 = time.perf_counter()
        p = prv_space(lam, lam, mu)
        hw = hw_from_slice(lam, lam, mu, p.float_basis())
        t = time.perf_counter() - t0
        cold = t if cold is None else min(cold, t)
    warm = None
    for _ in range(reps):
        t0 = time.perf_counter()
        p = prv_space(lam, lam, mu)
        hw = hw_from_slice(lam, lam, mu, p.float_basis())
        t = time.perf_counter() - t0
        warm = t if warm is None else min(warm, t)
    return cold, warm, hw.nW


# ======================================================================================== 8. multiplicity families
FAMILIES = [('[2a,a,0,0]^2', lambda a: ((2 * a, a, 0, 0),) * 2, range(1, 9)),
            ('[3a,2a,a,0]^2', lambda a: ((3 * a, 2 * a, a, 0),) * 2, range(1, 5)),
            ('[2a,a,a,0]^2', lambda a: ((2 * a, a, a, 0),) * 2, range(1, 9)),
            ('[a,a,a,0]x[a,a,0,0]', lambda a: ((a, a, a, 0), (a, a, 0, 0)), range(1, 5)),
            ('[a,a,0,0]x[a,0,0,0]', lambda a: ((a, a, 0, 0), (a, 0, 0, 0)), range(1, 7))]
BENCH_FAMILIES = [('[2a,a,0,0]^2', lambda a: ((2 * a, a, 0, 0),) * 2, range(1, 7)),
                  ('[3a,2a,a,0]^2', lambda a: ((3 * a, 2 * a, a, 0),) * 2, range(1, 4)),
                  ('[2a,a,a,0]^2', lambda a: ((2 * a, a, a, 0),) * 2, range(1, 6))]


def families(verbose=True):
    """multiplicity families (Racah-Speiser only, no vectors): (#irreps, sum m, max m, max w for either choice of site 1, max nW) vs a"""
    out = []
    for name, f, rng_ in FAMILIES:
        if verbose:
            print('  family %s' % name)
            print('   %2s %-22s %6s %6s %7s %7s %6s %9s %9s %8s %5s  %s' % ('a', 'lam (x) lam\'', 'dim', 'dim\'', 'irreps', 'sum m', 'max m', 'max w(1)', 'max w(2)', 'max nW', '#mmax',
                                                                    'first irrep with max m (descending lex.)'), flush=True)
        for a in rng_:
            lam, lamp = f(a)
            rs = racah_speiser(lam, lamp)
            wl, wp = weight_mults(lam), weight_mults(lamp)
            N = len(lam)
            w1 = max(wl.get(tuple(mu[k] - lamp[k] for k in range(N)), 0) for mu in rs)       # site 1 = lam (slice on V_lam[mu - lamp])
            w2 = max(wp.get(tuple(mu[k] - lam[k] for k in range(N)), 0) for mu in rs)         # site 1 = lamp (swapped)
            mm = max(rs.values())
            mumax = [mu for mu in sorted(rs, reverse=True) if rs[mu] == mm]
            row = dict(name=name, a=a, lam=lam, lamp=lamp, dim=weyl_dim(lam), dimp=weyl_dim(lamp), nirrep=len(rs), summ=sum(rs.values()), maxm=mm, w1=w1, w2=w2,
                       maxnW=max(nW_of(lam, lamp, mu) for mu in rs), nmax=len(mumax), mumax=mumax[0])
            out.append(row)
            if verbose:
                print('   %2d %-22s %6d %6d %7d %7d %6d %9d %9d %8d %5d  %s' % (a, prod_name(lam, lamp), row['dim'], row['dimp'], row['nirrep'], row['summ'], mm, w1, w2, row['maxnW'],
                                                                           row['nmax'], fmtl(mumax[0])), flush=True)
    return out


def bench_families(reps=5, verbose=True, fams=None):
    """t_a / t_b (bench_one: every irrep, every PRV basis vector, i.e. m vectors per irrep of multiplicity m) along the multiplicity families"""
    rows = []
    for name, f, rng_ in (fams or BENCH_FAMILIES):
        if verbose:
            print('  family %s' % name)
            print('   %2s %-14s %6s %6s %6s %5s %5s %8s | %8s | %8s | %11s | %8s' % ('a', 'lam', 'dimV', 'irreps', 'copies', 'max m', 'max w', 'max nW', 't_a', 't_b', 'entries', 'resid'), flush=True)
        for a in rng_:
            lam, lamp = f(a)
            T = bench_one(lam, lamp, reps=reps)
            T['family'], T['a'] = name, a
            rows.append(T)
            if verbose:
                print('   %2d %-14s %6d %6d %6d %5d %5d %8d | %8.4f | %8.4f | %11d | %8.1e   [t_gen %.4f t_rs %.4f t_prv %.4f t_setup %.4f t_recon %.4f]' %
                      (a, fmtl(T['lam']), T['dim'], T['nirrep'], T['ncopies'], T['maxm'], T['maxw'], T['maxnW'], T['ta'], T['tb'], T['nent'], T['resid'], T['gen'], T['rs'], T['prv'],
                       T['setup_b'], T['recon']), flush=True)
    return rows


# ======================================================================================== 9. reference files
def fq(x):
    x = Fr(x)
    return str(x.numerator) if x.denominator == 1 else '%d/%d' % (x.numerator, x.denominator)


def fmu(mu):
    return '{%s}' % ','.join(str(int(x)) for x in mu)


def ref_name(lam, lamp):
    s = ''.join(str(x) for x in lam)
    if _t(lam) != _t(lamp):
        s += '-' + ''.join(str(x) for x in lamp)
    return 'gl4-cg-slices-%s.txt' % s


def _hdr(lam, lamp, has_vec, exact_gn, date, name):
    N = len(lam)
    same = _t(lam) == _t(lamp)
    L = []
    A = L.append
    A('# %s   (WP-G1 of the gl(4) slice-engine programme; generator Paul/Python/gaudin_gl4/slice_engine_gl4.py `python -X utf8 slice_engine_gl4.py ref %s`; %s)' %
      (name, ' '.join(str(x) for x in lam) + ('' if same else ' ' + ' '.join(str(x) for x in lamp)), date))
    A('# Highest-weight (HW) Clebsch-Gordan data of gl(%d), V_Lambda (x) V_Lambda\' with Lambda = {%s} (REP, dim DIM) and Lambda\' = {%s} (REP2, dim DIM2): for EVERY irrep mu of the decomposition the PRV/Zhelobenko slice space M_mu (its exact canonical basis), and for the multiplicity-free irreps (m = 1) the normalised slice, the exact full HW vector (small products), its exact G-norm^2 and the float HW residual of the engine.' %
      (N, ','.join(str(x) for x in lam), ','.join(str(x) for x in lamp)))
    A('# FORMAT: one record per line, TAB-separated fields, first field = tag, other fields nested braces = lists (Wolfram InputForm style). p/q = exact rationals (integers without denominator); residuals are Python floats (e.g. 1.2e-15); NO other decimals are written (every slice/vector entry below is an exact rational). Lines starting with # are comments. Same record style as ref/cg-slices-420.txt (gl(3)).')
    A('# ORDER OF THE GT BASES (binding for every index below): the order of gl3_reps.gt_patterns(Lambda) = itertools.product order = ascending lexicographic in (row N-1 = second-to-top row, row N-2, ..., row 1 = singleton) of the pattern {{k},{m1,m2},{..},{l1,...,l%d}} (rows 1..%d); the HIGHEST-WEIGHT state is the LAST index (= DIM, resp. DIM2). All indices nu, q, idx are 1-BASED in this order. GTBASIS / GTBASIS2 list the patterns in exactly this order (nested braces). For N = 3 the order is that of GL3Rep(Lambda).patterns.' % (N, N))
    A('#   Matrices act on kets as columns, E e_col = sum_row M[row,col] e_row; plain rational non-unitary GT matrices of glN_reps.GLNRep (E_{k,k+1} = gl3_reps._raise(k,.), E_{k+1,k} = gl3_reps._lower(k,.), targets outside the pattern set dropped; E_ij = [E_{i,j-1},E_{j-1,j}], E_ji = [E_{j,j-1},E_{j-1,i}] for j > i+1; all 256 relations [E_ij,E_kl] = d_jk E_il - d_li E_kj verified exactly). WEIGHTS1/WEIGHTS2 = (E11,E22,E33,E44) eigenvalues per index. GMET1/GMET2 = invariant diagonal metric G (E_ij^T G = G E_ji for all i,j, G[HW] = 1), exact; orthonormal basis e~_q = e_q/sqrt(g_q).')
    A('# TWO SITES: Psi = sum Psi[nu,q] e_nu (x) e_q, site 1 = first factor V_Lambda (nu, slow), site 2 = second factor V_Lambda\' (q, fast); 1-based flat index idx = (nu-1)*DIM2 + q. The SLICE of a vector Psi of weight mu is psi_top[nu] = Psi[nu, q_top] (site 2 = HW state q_top = index DIM2), a vector in the weight space V_Lambda[mu-Lambda\'] of SITE 1 (weight space = all nu whose GT weight equals mu-Lambda\').')
    A('# THEOREM USED (Kostant 1959 / PRV 1967): Psi -> slice is an isomorphism from the HW vectors of weight mu onto M_mu = { v in V_Lambda[mu-Lambda\'] : E_{i,i+1}^(Lambda\'_i - Lambda\'_{i+1} + 1) v = 0, i = 1..%d }; dim M_mu = m = multiplicity of mu in V_Lambda (x) V_Lambda\' (Racah-Speiser, checked for every irrep and every dominant candidate). The string exponents come from Lambda\' (site 2), the weight space lives in V_Lambda (site 1). The full vector is recovered from the slice by the level recursion over the weights of site 2 (slice_engine_gl4.py; exact twin verifies all consistency and dangling conditions).' % (N - 1))
    A('# RECORDS: REP | REP2 | DIM | DIM2 | GTBASIS | GTBASIS2 | WEIGHTS1 | WEIGHTS2 | GMET1 | GMET2 | NSECTORS | NCOPIES | then per sector (descending lexicographic mu): SECTOR, WSPACE, PRVBASIS x m, then for m = 1: COPY, HWRESID, GNORM2, VEC (small products), and for m >= 2 (small products) PRVVEC x m |')
    A('#   SECTOR mu m w nW dimVmu cas2mu : m = multiplicity = dim M_mu; w = dim V_Lambda[mu-Lambda\'] = length of a slice; nW = dim of the weight-mu subspace of V_Lambda (x) V_Lambda\' (number of entries of a full vector); dimVmu = Weyl dimension of the irrep mu; cas2mu = Sum_i mu_i (mu_i + N + 1 - 2i) (N = %d) |' % N)
    A('#   WSPACE mu {nu,...} : the 1-based GT indices (ascending) of the weight space V_Lambda[mu-Lambda\'] (site 1); the slice / PRVBASIS positions are aligned with this list |')
    A('#   PRVBASIS mu {{nu,val},...} : one line per basis vector of M_mu, k = 1..m in line order; EXACT rationals; nonzero entries of the slice; canonical basis = RREF nullspace with columns in ascending GT index, basis vector k = 1 at the k-th free column and 0 at the other free columns (so the coordinates of a slice are its free entries). For m >= 2 these are NOT yet the canonical copies (WP-G3 will add the joint eigenvectors of the Gaudin missing-label operators) |')
    A('#   COPY mu 1 0 {{nu,val},...} : (m = 1 only) the normalised slice = the PRVBASIS vector divided by its first nonzero component (ascending nu), so the first nonzero value is exactly 1; EXACT rationals; the 0 is a PLACEHOLDER for the eigenvalue of the missing-label operators (filled by WP-G3; irrelevant for m = 1 up to a constant) |')
    A('#   HWRESID mu 1 residual : max over i < j of ||(E_ij (x) 1 + 1 (x) E_ij) Psi|| / ||Psi|| (invariant norm = Frobenius norm in the orthonormal basis) for the full vector reconstructed by the FLOAT recursion (double precision) from COPY mu 1, plus the total-weight check (deviation 0 asserted) |')
    A('#   GNORM2 mu 1 value : G-norm^2 = Sum_{nu,q} g_nu g_q Psi[nu,q]^2 (G = GMET1 (x) GMET2, invariant norm) of the FULL vector whose slice is COPY mu 1, %s |' %
      ('EXACT rational (from the exact recursion)' if exact_gn else 'double precision from the float recursion (e-notation); the exact recursion was not run for this file'))
    if has_vec:
        A('#   VEC mu 1 {{idx,val},...} : (m = 1, DIM*DIM2 <= 2000) the EXACT full HW vector whose slice at q_top is COPY mu 1 (same scale), nonzero entries only, EXACT rationals, idx = (nu-1)*DIM2 + q, 1-based |')
        A('#   PRVVEC mu k {{idx,val},...} : (m >= 2, DIM*DIM2 <= 2000; extension to the brief) the EXACT full HW vector of PRVBASIS vector k (slice at q_top = PRVBASIS k), same format as VEC |')
    else:
        A('#   (no VEC / PRVVEC records in this file: DIM*DIM2 > 2000; regenerate with slice_engine_gl4.hw_from_slice_exact / hw_from_slice)')
    return L


def write_ref(lam, lamp=None, vec=None, exact_gnorm=True, date='2026-10-22', verbose=True):
    """write ref/gl4-cg-slices-<lam>[-<lamp>].txt.  vec = None: write VEC/PRVVEC iff DIM*DIM2 <= 2000.  exact_gnorm: GNORM2 from the exact recursion (else float)."""
    lam = _t(lam)
    lamp = lam if lamp is None else _t(lamp)
    t0 = time.perf_counter()
    R1, R2 = GLNRep.of(lam), GLNRep.of(lamp)
    n1, n2 = R1.dim, R2.dim
    if vec is None:
        vec = (n1 * n2 <= 2000)
    rs = racah_speiser(lam, lamp)
    mus = sorted(rs, reverse=True)
    body = []
    stats = dict(maxres=0.0, nvec=0, ncopies=0, nexact=0, exact_viol=0, nchk=0, wdev=0, mism=0)
    for mu in mus:
        p = prv_space(lam, lamp, mu)
        m = p.m
        assert m == rs[mu], (mu, m, rs[mu])
        body.append('SECTOR\t%s\t%d\t%d\t%d\t%d\t%d' % (fmu(mu), m, p.w, nW_of(lam, lamp, mu), weyl_dim(mu), casimir2(mu)))
        body.append('WSPACE\t%s\t{%s}' % (fmu(mu), ','.join(str(int(c) + 1) for c in p.W)))
        for k in range(m):
            body.append('PRVBASIS\t%s\t{%s}' % (fmu(mu), ','.join('{%d,%s}' % (c + 1, fq(x)) for c, x in sorted(p.vec(k).items()))))
        if m == 1:
            v = p.basis[0]
            first = next(x for x in v if x != 0)
            vn = [x / first for x in v]
            sl = {int(c): x for c, x in zip(p.W, vn) if x != 0}
            body.append('COPY\t%s\t1\t0\t{%s}' % (fmu(mu), ','.join('{%d,%s}' % (c + 1, fq(x)) for c, x in sorted(sl.items()))))
            hw = hw_from_slice(lam, lamp, mu, np.array([float(x) for x in vn]))
            res, wd = hw_residual(hw, extra=True)
            assert wd == 0
            stats['maxres'] = max(stats['maxres'], res)
            stats['wdev'] = max(stats['wdev'], wd)
            stats['ncopies'] += 1
            body.append('HWRESID\t%s\t1\t%.2e' % (fmu(mu), res))
            if vec or exact_gnorm:
                psi, nc, vi = hw_from_slice_exact(lam, lamp, mu, vn)
                stats['nexact'] += 1
                stats['nchk'] += nc
                stats['exact_viol'] += exact_hw_check(lam, lamp, psi)
                gn = exact_gnorm2(lam, lamp, psi)
                body.append('GNORM2\t%s\t1\t%s' % (fmu(mu), fq(gn)))
                dn = float(np.sum(hw.vals_hat[:, 0] ** 2))
                stats['mism'] += (abs(dn - float(gn)) > 1e-9 * float(gn))
                if vec:
                    body.append('VEC\t%s\t1\t{%s}' % (fmu(mu), ','.join('{%d,%s}' % (nu * n2 + q + 1, fq(x)) for (nu, q), x in sorted(psi.items(), key=lambda kv: kv[0][0] * n2 + kv[0][1]) if x != 0)))
                    stats['nvec'] += 1
            else:
                dn = float(np.sum(hw.vals_hat[:, 0] ** 2))
                body.append('GNORM2\t%s\t1\t%.15e' % (fmu(mu), dn))
        elif vec:
            for k in range(m):
                psi, nc, vi = hw_from_slice_exact(lam, lamp, mu, p.basis[k])
                stats['nexact'] += 1
                stats['nchk'] += nc
                stats['exact_viol'] += exact_hw_check(lam, lamp, psi)
                body.append('PRVVEC\t%s\t%d\t{%s}' % (fmu(mu), k + 1, ','.join('{%d,%s}' % (nu * n2 + q + 1, fq(x)) for (nu, q), x in sorted(psi.items(), key=lambda kv: kv[0][0] * n2 + kv[0][1]) if x != 0)))
                stats['nvec'] += 1
    name = ref_name(lam, lamp)
    head = _hdr(lam, lamp, vec, bool(vec or exact_gnorm), date, name)
    head.append('# SUMMARY (self-check of this file): %d sectors, %d multiplicity-free irreps (COPY), %d exact full vectors written (VEC/PRVVEC), exact reconstructions run %d (%d exact conditions verified, %d exact HW violations over all E_ij with i<j); max HWRESID (float recursion, all E_ij) %.2e; weight deviation %d; GNORM2 vs float %s' %
                (len(mus), stats['ncopies'], stats['nvec'], stats['nexact'], stats['nchk'], stats['exact_viol'], stats['maxres'], stats['wdev'],
                 ('mismatches %d' % stats['mism']) if (vec or exact_gnorm) else 'n/a (float GNORM2 written)'))
    head.append('REP\t{%s}' % ','.join(str(x) for x in lam))
    head.append('REP2\t{%s}' % ','.join(str(x) for x in lamp))
    head.append('DIM\t%d' % n1)
    head.append('DIM2\t%d' % n2)
    gtb = lambda R: '{%s}' % ','.join('{%s}' % ','.join('{%s}' % ','.join(str(x) for x in row) for row in pat) for pat in R.pats)
    head.append('GTBASIS\t%s' % gtb(R1))
    head.append('GTBASIS2\t%s' % gtb(R2))
    head.append('WEIGHTS1\t{%s}' % ','.join('{%s}' % ','.join(str(x) for x in w) for w in R1.wt))
    head.append('WEIGHTS2\t{%s}' % ','.join('{%s}' % ','.join(str(x) for x in w) for w in R2.wt))
    head.append('GMET1\t{%s}' % ','.join(fq(x) for x in R1.g))
    head.append('GMET2\t{%s}' % ','.join(fq(x) for x in R2.g))
    head.append('NSECTORS\t%d' % len(mus))
    head.append('NCOPIES\t%d' % sum(rs.values()))
    path = os.path.join(REFDIR, name)
    with open(path, 'w', encoding='utf8', newline='\n') as f:
        f.write('\n'.join(head + body) + '\n')
    stats['t'] = time.perf_counter() - t0
    stats['path'] = path
    stats['bytes'] = os.path.getsize(path)
    stats['nsec'] = len(mus)
    stats['ncopies_all'] = sum(rs.values())
    if verbose:
        print('  wrote %s: %d sectors, %d copies (%d with m = 1), VEC/PRVVEC %d, exact conditions verified %d, exact HW violations %d, max HWRESID %.2e, GNORM2 float-vs-exact mismatches %d, %d bytes, %.1fs' %
              (path, len(mus), sum(rs.values()), stats['ncopies'], stats['nvec'], stats['nchk'], stats['exact_viol'], stats['maxres'], stats['mism'], stats['bytes'], stats['t']), flush=True)
    return stats


# ======================================================================================== 10. independent re-verification of the reference files
_TOK = re.compile(r'\{|\}|,|-?\d+/\d+|-?\d+\.\d+(?:[eE][-+]?\d+)?|-?\d+')


def parse_nested(s):
    """parse '{{1,2/3},{0.5}}' -> nested lists of int / Fraction / float (no eval)"""
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
            return float(t)
        return int(t)
    return val()


def parse_ref(path):
    """parse a gl4-cg-slices-*.txt: dict(lam, lamp, dim, dim2, gtbasis, gtbasis2, weights, weights2, g, g2, sectors={mu: dict(m, w, nW, dimV, cas2, wspace, prv=[{nu:val}], copy={nu:val}|None,
    hwres, gnorm2, vec={idx:val}|None, prvvec=[{idx:val}])})"""
    out = dict(sectors={})
    cur = None
    for line in open(path, encoding='utf8'):
        line = line.rstrip('\n')
        if not line or line.startswith('#'):
            continue
        f = line.split('\t')
        tag = f[0]
        if tag == 'REP':
            out['lam'] = tuple(parse_nested(f[1]))
        elif tag == 'REP2':
            out['lamp'] = tuple(parse_nested(f[1]))
        elif tag == 'DIM':
            out['dim'] = int(f[1])
        elif tag == 'DIM2':
            out['dim2'] = int(f[1])
        elif tag == 'GTBASIS':
            out['gtbasis'] = parse_nested(f[1])
        elif tag == 'GTBASIS2':
            out['gtbasis2'] = parse_nested(f[1])
        elif tag == 'WEIGHTS1':
            out['weights'] = [tuple(w) for w in parse_nested(f[1])]
        elif tag == 'WEIGHTS2':
            out['weights2'] = [tuple(w) for w in parse_nested(f[1])]
        elif tag == 'GMET1':
            out['g'] = [Fr(x) for x in parse_nested(f[1])]
        elif tag == 'GMET2':
            out['g2'] = [Fr(x) for x in parse_nested(f[1])]
        elif tag == 'NSECTORS':
            out['nsectors'] = int(f[1])
        elif tag == 'NCOPIES':
            out['ncopies'] = int(f[1])
        elif tag == 'SECTOR':
            mu = tuple(parse_nested(f[1]))
            cur = out['sectors'][mu] = dict(m=int(f[2]), w=int(f[3]), nW=int(f[4]), dimV=int(f[5]), cas2=int(f[6]), prv=[], copy=None, hwres=None, gnorm2=None, vec=None, prvvec=[])
        elif tag == 'WSPACE':
            cur['wspace'] = parse_nested(f[2])
        elif tag == 'PRVBASIS':
            cur['prv'].append({nu: Fr(x) for nu, x in parse_nested(f[2])})
        elif tag == 'COPY':
            cur['copy'] = {nu: Fr(x) for nu, x in parse_nested(f[4])}
            cur['copy_ev'] = f[3]
        elif tag == 'HWRESID':
            cur['hwres'] = float(f[3])
        elif tag == 'GNORM2':
            cur['gnorm2'] = Fr(f[3]) if '/' in f[3] or re.fullmatch(r'-?\d+', f[3]) else float(f[3])
        elif tag == 'VEC':
            cur['vec'] = {i: Fr(x) for i, x in parse_nested(f[3])}
        elif tag == 'PRVVEC':
            cur['prvvec'].append({i: Fr(x) for i, x in parse_nested(f[3])})
    return out


def verify_ref(path, verbose=True):
    """independent checks of a written reference file, parsed back from the TEXT (nothing of the engine's recursion/PRV code is used):
    (1) GTBASIS/GTBASIS2/WEIGHTS = glN_reps order, GMET exact and E_ij^T G = G E_ji with the DENSE float matrices GLNRep.dense (all i<j);
    (2) SECTOR: m = Racah-Speiser, dimVmu = Weyl, cas2 = formula, w = |WSPACE|, nW, WSPACE = the weight space computed from the parsed WEIGHTS;
    (3) PRVBASIS: dense-matrix string conditions E_i^k v = 0 (float, relative), rank m, and an independent float SVD of the string-condition matrix on the weight space
        gives null dimension m;
    (4) COPY = PRVBASIS / first nonzero component (exact); VEC / PRVVEC (exact rational): slice at q_top, total weight mu, HW conditions E_ij^tot Psi = 0 for all i<j
        in exact arithmetic AND with dense float matrices, GNORM2 exact;
    (5) records without VEC: the vector reconstructed by the float recursion from the FILE's slice satisfies the HW conditions with the dense float matrices (all i<j), GNORM2 vs float."""
    D = parse_ref(path)
    lam, lamp = D['lam'], D['lamp']
    R1, R2 = GLNRep.of(lam), GLNRep.of(lamp)
    N = R1.N
    n1, n2 = R1.dim, R2.dim
    res = dict(path=os.path.basename(path), lam=lam, lamp=lamp, nsec=len(D['sectors']), nprv=0, nvec=0, ncopy=0)
    bad = defaultdict(int)
    if [tuple(tuple(r) for r in p) for p in D['gtbasis']] != R1.pats or [tuple(tuple(r) for r in p) for p in D['gtbasis2']] != R2.pats:
        bad['gtbasis'] += 1
    if D['weights'] != R1.wt or D['weights2'] != R2.wt:
        bad['weights'] += 1
    if D['g'] != R1.g or D['g2'] != R2.g or D['dim'] != n1 or D['dim2'] != n2:
        bad['gmet_dim'] += 1
    if bad:                                                      # the basic data are wrong: nothing below is meaningful
        res.update(total_bad=sum(bad.values()), bad=dict(bad))
        if verbose:
            print('  %s: HEADER MISMATCH %s -> TOTAL MISMATCHES %d (remaining checks skipped)' % (res['path'], dict(bad), res['total_bad']), flush=True)
        return res
    E1 = {(i, j): R1.dense((i, j)) for i in range(1, N + 1) for j in range(1, N + 1)}
    E2 = E1 if lam == lamp else {(i, j): R2.dense((i, j)) for i in range(1, N + 1) for j in range(1, N + 1)}
    g1 = np.array([float(x) for x in D['g']])
    g2 = np.array([float(x) for x in D['g2']])
    sg1, sg2 = np.sqrt(g1), np.sqrt(g2)
    mdev = 0.0
    for (E, g) in ((E1, g1), (E2, g2)):
        for i in range(1, N + 1):
            for j in range(i + 1, N + 1):
                L = E[(i, j)].T * g[None, :]
                Rr = g[:, None] * E[(j, i)]
                mdev = max(mdev, np.abs(L - Rr).max() / np.abs(Rr).max())
    res['metric_dev'] = mdev
    rs = racah_speiser(lam, lamp)
    if set(rs) != set(D['sectors']) or D['nsectors'] != len(rs) or D['ncopies'] != sum(rs.values()):
        bad['sector_set'] += 1
    w1 = np.array(D['weights'])
    w2 = np.array(D['weights2'])
    top1, top2 = n1 - 1, n2 - 1
    # orthonormal-basis simple raising operators for the independent null-space dimension (from the PARSED metric)
    Eh1 = {ij: E1[ij] * sg1[:, None] / sg1[None, :] for ij in E1}
    prv_dev = 0.0
    nullbad = 0
    for mu, S in D['sectors'].items():
        m = S['m']
        if m != rs.get(mu, -1) or S['dimV'] != weyl_dim(mu) or S['cas2'] != casimir2(mu) or S['w'] != len(S['wspace']) or S['nW'] != nW_of(lam, lamp, mu):
            bad['sector_fields'] += 1
        wexp = [c + 1 for c in range(n1) if tuple(w1[c]) == tuple(mu[k] - lamp[k] for k in range(N))]
        if S['wspace'] != wexp:
            bad['wspace'] += 1
        if len(S['prv']) != m:
            bad['prv_count'] += 1
            continue
        vs = []
        for vd in S['prv']:
            v = np.zeros(n1)
            for nu, x in vd.items():
                v[nu - 1] = float(x)
                if nu not in wexp:
                    bad['prv_support'] += 1
            vs.append(v)
            for i in range(1, N):
                k = lamp[i - 1] - lamp[i] + 1
                Mk = np.linalg.matrix_power(E1[(i, i + 1)], k)
                sc = np.linalg.norm(v) * np.abs(Mk).max()
                if sc > 0:                                           # E_i^k = 0 identically on V_lam: the condition is void
                    prv_dev = max(prv_dev, np.abs(Mk @ v).max() / sc)
                ev = {nu - 1: x for nu, x in vd.items()}             # the same condition in exact arithmetic (GLNRep.apply, Fractions)
                for _ in range(k):
                    ev = R1.apply((i, i + 1), ev)
                    if not ev:
                        break
                if ev:
                    bad['prv_exact'] += 1
            res['nprv'] += 1
        if m:
            if np.linalg.svd(np.array(vs), compute_uv=False)[-1] < 1e-8 * np.linalg.norm(vs[0]):
                bad['prv_rank'] += 1
        # independent float null-space dimension of the string-condition matrix on the weight space (orthonormal basis)
        if wexp:
            cols = [c - 1 for c in wexp]
            rows = []
            for i in range(1, N):
                k = lamp[i - 1] - lamp[i] + 1
                Mk = np.linalg.matrix_power(Eh1[(i, i + 1)], k)
                rows.append(Mk[:, cols])
            Ast = np.vstack(rows)
            sv = np.linalg.svd(Ast, compute_uv=False)
            thr = 1e-9 * max(sv[0], 1.0) if len(sv) else 0
            nullity = len(cols) - int(np.sum(sv > thr))
            if nullity != m:
                nullbad += 1
        # COPY = normalised PRVBASIS (m = 1)
        if m == 1:
            v = S['prv'][0]
            first = v[min(v)]
            cp = {nu: x / first for nu, x in v.items()}
            if S['copy'] != cp:
                bad['copy'] += 1
            res['ncopy'] += 1
    res['prv_dev'] = prv_dev
    res['null_dim_bad'] = nullbad
    # full vectors
    hw_exact_viol = 0
    worst_float = 0.0
    worst_gn = 0.0
    worst_slice = 0
    gg = sg1[:, None] * sg2[None, :]
    for mu, S in D['sectors'].items():
        items = []
        if len(S['prv']) != S['m'] or (S['m'] == 1 and (S['copy'] is None or S['gnorm2'] is None)) or len(S['prvvec']) > len(S['prv']):
            bad['records_missing'] += 1
            continue
        if S['m'] == 1:
            items.append((S['vec'], S['copy'], S['gnorm2']))
        for k, vd in enumerate(S['prvvec']):
            items.append((vd, S['prv'][k], None))
        for (vd, sl, gn) in items:
            if vd is not None:
                res['nvec'] += 1
                psi = {((i - 1) // n2, (i - 1) % n2): x for i, x in vd.items()}
                # slice at q_top equals the slice record; total weight = mu
                got = {nu + 1: x for (nu, q), x in psi.items() if q == top2}
                if got != sl:
                    worst_slice += 1
                if any(tuple(w1[nu] + w2[q]) != tuple(mu) for (nu, q) in psi):
                    bad['vec_weight'] += 1
                hw_exact_viol += exact_hw_check(lam, lamp, psi)
                X = np.zeros((n1, n2))
                for (nu, q), x in psi.items():
                    X[nu, q] = float(x)
                nx = np.linalg.norm(gg * X)
                for i in range(1, N + 1):
                    for j in range(i + 1, N + 1):
                        Rm = E1[(i, j)] @ X + X @ E2[(i, j)].T
                        worst_float = max(worst_float, np.linalg.norm(gg * Rm) / nx)
                if gn is not None and isinstance(gn, Fr):
                    if gn != exact_gnorm2(lam, lamp, psi):
                        bad['gnorm2_exact'] += 1
            elif gn is not None:                      # no VEC: reconstruct from the file's slice (float recursion), check with the dense independent matrices
                wsp = [c for c in range(n1) if tuple(w1[c]) == tuple(mu[k] - lamp[k] for k in range(N))]
                v = np.array([float(sl.get(c + 1, 0)) for c in wsp])
                hw = hw_from_slice(lam, lamp, mu, v)
                X = hw.dense(0)
                nx = np.linalg.norm(gg * X)
                for i in range(1, N + 1):
                    for j in range(i + 1, N + 1):
                        Rm = E1[(i, j)] @ X + X @ E2[(i, j)].T
                        worst_float = max(worst_float, np.linalg.norm(gg * Rm) / nx)
                if np.abs(X[:, top2][wsp] - v).max() > 1e-12 * np.abs(v).max():
                    worst_slice += 1
                dn = float(np.sum((gg * X) ** 2))
                worst_gn = max(worst_gn, abs(dn - float(gn)) / float(gn))
    res.update(hw_exact_viol=hw_exact_viol, hw_float=worst_float, gnorm_dev=worst_gn, slice_bad=worst_slice, bad=dict(bad))
    res['total_bad'] = sum(bad.values()) + hw_exact_viol + worst_slice + nullbad + (1 if (mdev > 1e-12 or prv_dev > 1e-12 or worst_float > 1e-12 or worst_gn > 1e-9) else 0)
    if verbose:
        print('  %s: sectors %d, PRVBASIS vectors %d (m = 1 copies %d), full vectors checked from the text %d | GTBASIS/WEIGHTS/GMET/sector set mismatches %d | metric E^T G = G E\' (dense floats) dev %.1e |'
              ' PRVBASIS string conditions rel. dev %.1e, sectors whose independent float null dimension != m: %d | other field mismatches %s | exact full vectors: exact HW violations (all E_ij, i<j) %d, '
              'slice/weight mismatches %d | float HW residual (dense independent matrices, all E_ij) %.1e | GNORM2 vs float %.1e | TOTAL MISMATCHES %d' %
              (res['path'], res['nsec'], res['nprv'], res['ncopy'], res['nvec'], bad['gtbasis'] + bad['weights'] + bad['gmet_dim'] + bad['sector_set'], mdev, prv_dev, nullbad,
               {k: v for k, v in bad.items() if k not in ('gtbasis', 'weights', 'gmet_dim', 'sector_set')}, hw_exact_viol, worst_slice, worst_float, worst_gn, res['total_bad']), flush=True)
    return res


REF_LIST = [((2, 1, 0, 0), None), ((2, 1, 1, 0), None), ((3, 2, 1, 0), None), ((2, 1, 0, 0), (1, 1, 0, 0)), ((2, 2, 1, 0), (2, 1, 1, 0)), ((2, 2, 0, 0), None),
            ((3, 3, 0, 0), None), ((4, 4, 0, 0), None), ((6, 6, 0, 0), None), ((8, 8, 0, 0), None)]


# ======================================================================================== 11. bridge to Munich's code (Paul/Cpp/AlexCGC/ClebschGordan.cpp)
# The binary cannot be built on this machine, so the literal comparison of HW vectors is deferred to the Ubuntu laptop.  What CAN be settled offline is that the two
# codes use the same basis: (i) Munich's pattern index Q(M) (clebsch::pattern::operator++, Eq. C7/C9) = our GT index + 1; (ii) Munich's orthonormal GT matrix elements
# (clebsch::pattern::lowering_coeff, Eq. 28) equal ours in absolute value; they differ at most by a diagonal sign gauge eps_q = +-1 (Munich's lowering coefficients are
# sqrt(-...) >= 0).  Then   Munich_CGC[Q(M), Q(M'), Q(M'')=dim] = sigma * eps_nu eps_q Psi[nu,q] sqrt(g_nu g_q / GNORM2),  sigma = +-1 a global sign (LAPACK null vector).
def to_munich(lam, lamp, hw, k=0):
    """our HW vector k (HW object) -> {(Q(M), Q(M')): C} in Munich's orthonormal GT gauge, normalised sum C^2 = 1, global sign = that of our plain-basis vector
    (Munich's own sign is arbitrary).  C[nu,q] = eps_nu eps_q Psi_hat[nu,q] / ||Psi_hat||,  Psi_hat = orthonormal-basis components."""
    e1 = np.array(munich_check(lam)['eps'])
    e2 = np.array(munich_check(lamp)['eps']) if _t(lam) != _t(lamp) else e1
    v = hw.vals_hat[:, k] * e1[hw.rows] * e2[hw.cols]
    v = v / np.linalg.norm(v)
    return {(int(r) + 1, int(c) + 1): float(x) for r, c, x in zip(hw.rows, hw.cols, v)}


def parse_munich_rows(text):
    """rows 'Q(M)<TAB>Q(M')<TAB>Q(M'')<TAB>CGC' (the hwcgc.sh / clebsch output) -> {(Q(M), Q(M'), Q(M'')): value}"""
    out = {}
    for line in text.splitlines():
        f = line.replace(',', ' ').split()
        if len(f) == 4 and all(re.fullmatch(r'\d+', x) for x in f[:3]):
            try:
                out[(int(f[0]), int(f[1]), int(f[2]))] = float(f[3])
            except ValueError:
                pass
    return out


def munich_compare(lam, lamp, mu, text, tol=1e-9, alpha=1):
    """compare the HW rows printed by `bash hwcgc.sh 4 "lam" "lamp" "mu_munich"` (text = its stdout) with our vector: same support (|C| > tol), values equal up to ONE global
    sign, max deviation <= tol.  mu = OUR label (Munich's = mu minus mu_4); only valid for multiplicity-free mu.  Returns (n_ours, n_munich, support_mismatch, max_dev, sigma)."""
    lam, lamp, mu = _t(lam), _t(lamp), _t(mu)
    p = prv_space(lam, lamp, mu)
    assert p.m == 1, 'munich_compare needs a multiplicity-free irrep'
    hw = hw_from_slice(lam, lamp, mu, p.float_basis())
    ours = to_munich(lam, lamp, hw, 0)
    rows = parse_munich_rows(text)
    dimmu = weyl_dim(mu)
    mun = {(a, b): v for (a, b, c), v in rows.items() if c == dimmu}
    keys_o = {k for k, v in ours.items() if abs(v) > tol}
    keys_m = {k for k, v in mun.items() if abs(v) > tol}
    kbig = max(ours, key=lambda k: abs(ours[k]))
    sigma = 1.0 if mun.get(kbig, 0.0) * ours[kbig] >= 0 else -1.0
    dev = max([abs(sigma * ours.get(k, 0.0) - mun.get(k, 0.0)) for k in keys_o | keys_m] + [0.0])
    return len(keys_o), len(keys_m), len(keys_o ^ keys_m), dev, sigma


def munich_commands(Slist=(1, 2, 3), verbose=True):
    """the exact commands for the Ubuntu session: for every irrep mu of [S,S,0,0]^2 (Munich label) the hwcgc.sh call and the comparison with our vector"""
    out = []
    for S in Slist:
        lam = (S, S, 0, 0)
        for mu in sorted(racah_speiser(lam, lam), reverse=True):
            mm = tuple(x - mu[3] for x in mu)
            assert munich_label_to_ours(lam, lam, mm) == mu
            f = 'hw_%d_%s.txt' % (S, '_'.join(str(x) for x in mm))
            c1 = 'bash hwcgc.sh 4 "%d %d 0 0" "%d %d 0 0" "%s" > %s' % (S, S, S, S, ' '.join(str(x) for x in mm), f)
            c2 = 'python -X utf8 slice_engine_gl4.py munichcmp %d %s %s' % (S, ' '.join(str(x) for x in mm), f)
            out.append((S, mm, mu, c1, c2))
            if verbose:
                print('  %s    && %s' % (c1, c2), flush=True)
    return out


def time_all_generators(S, reps=3):
    """cost of what t_a does NOT include: the 6 non-adjacent exact generators E_13, E_24, E_14, E_31, E_42, E_41 (commutators) and then all 12 float CSR matrices (orthonormal)"""
    best = None
    for _ in range(reps):
        _clear_caches()
        R = GLNRep.of((S, S, 0, 0))
        t0 = time.perf_counter()
        for i in range(1, 5):
            for j in range(1, 5):
                if i != j:
                    R.E[(i, j)]
        t1 = time.perf_counter()
        g = R.g
        Eh = R.float_ops()
        for i in range(1, 5):
            for j in range(1, 5):
                Eh[(i, j)]
        t2 = time.perf_counter()
        tt = (t1 - t0, t2 - t1)
        if best is None or sum(tt) < sum(best):
            best = tt
    return best


def munich_selftest(verbose=True):
    """the Munich bridge without the binary: (a) pattern order and orthonormal lowering coefficients identical to ClebschGordan.cpp's formulas for all irreps used;
    (b) the comparer accepts synthetic hwcgc.sh output generated from our vectors (global sign flipped, 15 decimals) and rejects corrupted output (changed value,
    changed index, dropped row)"""
    out = []
    lams = [(1, 0, 0, 0), (1, 1, 0, 0), (2, 1, 0, 0), (2, 1, 1, 0), (3, 2, 1, 0), (2, 0, 0, 0)] + [(S, S, 0, 0) for S in range(1, 9)] + [(10, 10, 0, 0)]
    for lam in lams:
        r = munich_check(lam)
        out.append((lam, r))
        if verbose:
            print('  %-14s Munich pattern order = ours: %s | %6d lowering matrix elements: Munich lowering_coeff^2 != ours: %d, extra nonzero entries of ours: %d | sign gauge eps: inconsistent/undefined %d, trivial (all +1): %s' %
                  (fmtl(lam), r['order_ok'], r['nelem'], r['bad_sq'], r['extra'], r['gauge_bad'], r['trivial_gauge']), flush=True)
    nacc = nrej = nbad = ncorr = 0
    for S in (1, 2, 3):
        lam = (S, S, 0, 0)
        for mu in sorted(racah_speiser(lam, lam), reverse=True):
            p = prv_space(lam, lam, mu)
            hw = hw_from_slice(lam, lam, mu, p.float_basis())
            C = to_munich(lam, lam, hw, 0)
            dimmu = weyl_dim(mu)
            txt = "Q(M)\tQ(M')\tQ(M'')\tCGC\n" + '\n'.join('%d\t%d\t%d\t%.15f' % (a, b, dimmu, -v) for (a, b), v in sorted(C.items()) if abs(v) > 1e-12)       # Munich prints nonzero CGCs only
            r = munich_compare(lam, lam, mu, txt)                       # global sign flipped on purpose
            nacc += 1
            nbad += (r[2] != 0 or r[3] > 1e-12)
            lines = txt.split('\n')
            if len(lines) > 3:                                           # corrupted variants: change a value, change an index, drop a row
                f = lines[2].split('\t')
                mut1 = lines[:2] + ['\t'.join(f[:3] + ['%.15f' % (float(f[3]) + 1e-3)])] + lines[3:]
                mut2 = lines[:2] + ['\t'.join([str(int(f[0]) + 1)] + f[1:])] + lines[3:]
                mut3 = lines[:2] + lines[3:]
                for mut in (mut1, mut2, mut3):
                    ncorr += 1
                    r2 = munich_compare(lam, lam, mu, '\n'.join(mut))
                    nrej += (r2[2] != 0 or r2[3] > 1e-6)
    if verbose:
        print('  comparer self-test: %d synthetic hwcgc.sh outputs accepted with %d failures; corrupted outputs rejected: %d of %d' % (nacc, nbad, nrej, ncorr), flush=True)
    return out, (nacc, nbad, nrej, ncorr)


# ======================================================================================== 12. extras: per-irrep timings vs Munich, N = 5 smoke test, brute-force extension
def per_irrep_times(S, reps=3, verbose=True):
    """S = 8 style table: for every irrep of [S,S,0,0]^2 Munich's hw-solve time (i5-1345U) vs ours (warm: the per-lam setup already done; prv_space + hw_from_slice)"""
    mun = parse_munich()[S]['rows']
    lam = (S, S, 0, 0)
    Recursion.of(lam, lam)
    rows = []
    for mu in sorted(racah_speiser(lam, lam), reverse=True):
        mm = tuple(x - mu[3] for x in mu)
        best = None
        for _ in range(reps):
            t0 = time.perf_counter()
            p = prv_space(lam, lam, mu)
            t1 = time.perf_counter()
            hw = hw_from_slice(lam, lam, mu, p.float_basis())
            t2 = time.perf_counter()
            tt = (t1 - t0, t2 - t1)
            if best is None or sum(tt) < sum(best):
                best = tt
        rows.append((mu, mm, mun[mm][0], hw.nW, mun[mm][2], best[0], best[1]))
    if verbose:
        print('  [%d,%d,0,0]^2 per irrep (ours: warm, after the per-lam setup; Munich: i5-1345U, its hw-solve column)' % (S, S))
        print('   %-14s %-14s %7s %7s %12s %10s %10s %10s %9s' % ('mu (ours)', 'mu (Munich)', 'dim', 'nW', 'Munich [s]', 't_prv [s]', 't_recon[s]', 'ours [s]', 'speedup'))
        for (mu, mm, d, nW, tm, tp, tr) in rows:
            print('   %-14s %-14s %7d %7d %12.6f %10.6f %10.6f %10.6f %9.0f' % (mu, mm, d, nW, tm, tp, tr, tp + tr, tm / (tp + tr)))
        print('   sum Munich %.3f s, sum ours (warm, excl. setup) %.4f s' % (sum(r[4] for r in rows), sum(r[5] + r[6] for r in rows)), flush=True)
    return rows


N5_PRODUCTS = [((1, 0, 0, 0, 0),) * 2, ((1, 1, 0, 0, 0),) * 2, ((2, 1, 0, 0, 0), (1, 1, 0, 0, 0)), ((1, 1, 1, 0, 0), (1, 1, 0, 0, 0)), ((2, 1, 1, 0, 0),) * 2,
               ((2, 2, 0, 0, 0),) * 2, ((2, 1, 0, 0, 0),) * 2, ((2, 2, 1, 0, 0), (2, 1, 1, 1, 0)), ((3, 2, 0, 0, 0), (2, 1, 0, 0, 0))]


def smoke_n(products, label, brute=True, nrand=10):
    """the code is general in N: PRV = Racah-Speiser for every dominant candidate, float and exact recursion + exact HW check on every irrep, negative tests, brute force"""
    print('  %s: PRV vs Racah-Speiser' % label)
    table_prv(products)
    print('  %s: float recursion (all E_ij, i < j)' % label)
    test_hw(products)
    print('  %s: exact twin' % label)
    exact_tests(products)
    print('  %s: negative tests' % label)
    neg_tests(products, nrand=nrand)
    if brute:
        print('  %s: brute-force null space' % label)
        for p in products:
            bruteforce_check(*p)


# ======================================================================================== 13. mutation test of the independent verifier
def mutation_test_verify(name='gl4-cg-slices-2100.txt', verbose=True):
    """verify_ref must DETECT corruptions: copy a reference file, corrupt one record at a time (a PRVBASIS value, a COPY value, a VEC value, a GNORM2 value, a GMET1 value,
    a SECTOR multiplicity, a WSPACE index, a WEIGHTS1 entry, an out-of-range VEC index/weight, a dropped PRVBASIS), and check that verify_ref reports a nonzero number of mismatches.
    The pristine copy must report 0."""
    import tempfile
    src = os.path.join(REFDIR, name)
    lines = open(src, encoding='utf8').read().split('\n')

    def find(tag, nth=0, pred=None):
        c = 0
        for i, ln in enumerate(lines):
            if ln.startswith(tag + '\t') and (pred is None or pred(ln)):
                if c == nth:
                    return i
                c += 1
        raise KeyError(tag)

    def bump(field_idx, newval, tag, nth=0, pred=None):
        out = list(lines)
        i = find(tag, nth, pred)
        f = out[i].split('\t')
        f[field_idx] = newval(f[field_idx])
        out[i] = '\t'.join(f)
        return out

    def first_val(s):                                   # '{{4,-1/2},{8,1}}' -> change the first value
        m = re.search(r',(-?\d+(?:/\d+)?)\}', s)
        v = Fr(m.group(1)) + Fr(1, 7)
        return s[:m.start(1)] + fq(v) + s[m.end(1):]

    muts = {
        'pristine copy': lines,
        'PRVBASIS value changed (vector with >= 2 entries)': bump(2, first_val, 'PRVBASIS', 0, lambda ln: ln.split('\t')[2].count('{') > 2),
        'PRVBASIS scaled (m = 2 sector, pinned by PRVVEC)': bump(2, first_val, 'PRVBASIS', 0, lambda ln: ln.split('\t')[1] == '{3,2,1,0}'),
        'COPY value changed': bump(4, first_val, 'COPY', 1, lambda ln: ln.count('{') > 2),
        'VEC value changed': bump(3, first_val, 'VEC', 1, lambda ln: ln.count('{') > 3),
        'GNORM2 changed': bump(3, lambda s: '%s' % (Fr(s) + 1), 'GNORM2', 1),
        'GMET1 entry changed': bump(1, lambda s: s.replace('{', '{7777,', 1), 'GMET1'),
        'SECTOR multiplicity changed': bump(2, lambda s: str(int(s) + 1), 'SECTOR', 2),
        'WSPACE index changed': bump(2, lambda s: s.replace('{', '{1,', 1), 'WSPACE', 3),
        'WEIGHTS1 entry changed': bump(1, lambda s: s.replace('{{', '{{9,', 1), 'WEIGHTS1'),
        'VEC index moved': bump(3, lambda s: re.sub(r'\{(\d+),', lambda m: '{%d,' % (int(m.group(1)) + 1), s, count=1), 'VEC', 1, lambda ln: ln.count('{') > 3),
        'PRVBASIS record dropped': [ln for k, ln in enumerate(lines) if k != find('PRVBASIS', 6)],
    }
    res = []
    with tempfile.TemporaryDirectory() as td:
        for label, ls in muts.items():
            pth = os.path.join(td, name)
            with open(pth, 'w', encoding='utf8', newline='\n') as f:
                f.write('\n'.join(ls))
            try:
                r = verify_ref(pth, verbose=False)
                tot = r['total_bad']
            except Exception as ex:                                 # a crash on corrupted input counts as detection
                tot = -1
            res.append((label, tot))
            if verbose:
                print('  %-30s verify_ref total mismatches: %s -> %s' % (label, tot if tot >= 0 else 'exception', 'OK (0 expected)' if (label == 'pristine copy' and tot == 0) else ('DETECTED' if tot != 0 else 'NOT DETECTED')), flush=True)
    return res


# ======================================================================================== CLI
def run_tests():
    import test_glN_reps as T
    print('== G1.1 glN_reps.py tests (test_glN_reps.py)')
    return T.run_all()


BRUTE_EXT = [(S, S, 0, 0) for S in (6, 7, 8, 9, 10)]
NEG_PRODUCTS = [((2, 1, 0, 0),) * 2, ((2, 1, 1, 0),) * 2, ((2, 2, 0, 0),) * 2, ((3, 2, 1, 0),) * 2, ((2, 1, 0, 0), (1, 1, 0, 0)), ((3, 1, 0, 0), (2, 2, 0, 0)),
                ((2, 2, 1, 0), (2, 1, 1, 0)), ((3, 3, 0, 0),) * 2, ((3, 2, 1, 0), (2, 1, 0, 0)), ((4, 4, 0, 0),) * 2, ((6, 6, 0, 0),) * 2, ((8, 8, 0, 0),) * 2]


def main(argv):
    cmd = argv[1] if len(argv) > 1 else 'all'
    ALL = cmd == 'all'
    if cmd in ('tests', 'all'):
        run_tests()
    if cmd in ('prv', 'all'):
        print('== G1.2: dim M_mu (PRV, three string conditions) vs Racah-Speiser, every dominant candidate mu (union of lam + wt(lamp) and lamp + wt(lam)), incl. multiplicity 0')
        table_prv(ALLPRODUCTS + BIG)
        print('   N = 3 regression against slice_engine.prv_space / hw_from_slice / racah_speiser')
        regression_n3()
    if cmd in ('hw', 'all'):
        print('== G1.3: float level recursion on EVERY irrep (all m PRV basis vectors): HW residual of all E_ij (i<j), total weight')
        test_hw(ALLPRODUCTS + BIG)
        print('   [8,8,0,0]^2, all 45 irreps:')
        hw_table_rect(8)
    if cmd in ('exact', 'all'):
        print('== G1.3: exact twin on EVERY irrep, EVERY PRV basis vector of all products of G1.2 and [S,S,0,0]^2 (S <= 10) (T + dangling conditions verified; exact HW check of all E_ij, i<j)')
        exact_tests(ALLPRODUCTS)
    if cmd in ('neg', 'all'):
        print('== G1.3: negative tests (random exact slices outside M_mu must be rejected by the exact conditions; inside accepted; conditions <=> exact HW cross-checked)')
        neg_tests(NEG_PRODUCTS, nrand=20)
        neg_tests([((6, 4, 2, 0),) * 2], nrand=5)
    if cmd in ('brute', 'all'):
        print('== G1.4: Munich-style brute-force null space (no PRV, no recursion): all products of G1.2 with dim*dim\' <= 10^4, [S,S,0,0]^2 for S <= 5, and the extension S = 6..10 (Munich\'s whole range S <= 8, and beyond)')
        for p in PRODUCTS:
            if GLNRep.of(p[0]).dim * GLNRep.of(p[1]).dim <= 10 ** 4:
                bruteforce_check(*p)
        for S in range(1, 6):
            bruteforce_check((S, S, 0, 0), (S, S, 0, 0))
        for lam in BRUTE_EXT:
            bruteforce_check(lam, lam)
        print('   multiplicity families of G3 ([4,2,0,0]^2: m <= 3; [6,4,2,0]^2: m <= 12, the irreps with n_W <= 2600):')
        bruteforce_check((4, 2, 0, 0), (4, 2, 0, 0))
        bruteforce_check((6, 4, 2, 0), (6, 4, 2, 0), maxnW=2600)
    if cmd in ('bench', 'all'):
        Smax = int(argv[2]) if len(argv) > 2 and cmd == 'bench' else 12
        print('== G1.5: benchmark [S,S,0,0]^2 (times in s; min of 5 cold runs, caches cleared, single process; t_a = gen + RS + PRV (exact); t_b = t_a + setup + recon of every vector (float))')
        print('   Munich column: i5-1345U (13th gen, 12 threads), g++ 11.4 + LAPACK 3.10, reps = 1  -- a DIFFERENT machine from ours (i7-12700KF)')
        bench(1, Smax)
        if cmd == 'all' or (len(argv) > 3 and argv[3] == 'ext'):
            print('   extension beyond the requested range (no Munich data):')
            for S in (14, 16, 20):
                bench(S, S, header=False)
        print('== float64 recursion vs exact recursion rounded to double (the accurate path) on the worst irreps')
        for S in (16, 20, 24):
            rounded_vs_float(S)
        print('== G1.5: Munich\'s slowest rows: per irrep at S = 8 (ours warm after the per-lam setup; plus cold single-irrep times)')
        per_irrep_times(8)
        for S in (6, 8, 10):
            cold, warm, nW = time_singlet(S)
            print('  singlet (S,S,S,S) of [%d,%d,0,0]^2: nW = %d, cold (everything incl. generators, setup) %.4f s, warm (prv + recon only) %.4f s' % (S, S, nW, cold, warm), flush=True)
    if cmd in ('families', 'all'):
        print('== G1.5: multiplicity families (Racah-Speiser only)')
        families()
        print('== G1.5: t_a / t_b along the families (all irreps, all m PRV basis vectors per irrep; min of 5 cold runs)')
        bench_families()
    if cmd in ('ref', 'all'):
        args = (argv[2:] if cmd == 'ref' else ['all']) or ['all']
        if args[0] == 'all':
            todo = REF_LIST
        else:
            nums = [int(x) for x in args]
            lam = tuple(nums[:4])
            lamp = tuple(nums[4:8]) if len(nums) >= 8 else None
            todo = [(lam, lamp)]
        print('== G1.6: reference files')
        for lam, lamp in todo:
            write_ref(lam, lamp)
    if cmd in ('verify', 'all'):
        print('== G1.6: independent re-verification of the written reference files (parsed back from the text)')
        for lam, lamp in REF_LIST:
            verify_ref(os.path.join(REFDIR, ref_name(lam, lamp or lam)))
        print('== G1.6: mutation test of verify_ref (it must detect corrupted records)')
        mutation_test_verify()
    if cmd in ('munich', 'all'):
        print('== Munich bridge (binary not available here): index order, orthonormal lowering coefficients, sign gauge, comparer self-test')
        munich_selftest()
        print('   commands for the Ubuntu session (run from Paul/Cpp/AlexCGC after `bash build.sh`; the comparison script is run from Paul/Python/gaudin_gl4 with the output file copied there):')
        munich_commands((1, 2, 3))
    if cmd == 'munichcmp':
        S = int(argv[2])
        mm = tuple(int(x) for x in argv[3:7])
        text = open(argv[7], encoding='utf8').read()
        lam = (S, S, 0, 0)
        mu = munich_label_to_ours(lam, lam, mm)
        n_o, n_m, nsup, dev, sigma = munich_compare(lam, lam, mu, text)
        print('  Munich (%s) vs ours (mu = %s) for [%d,%d,0,0]^2: nonzero rows ours %d, Munich %d, support mismatches %d, max |sigma*ours - Munich| %.2e (global sign sigma = %+d) -> %s' %
              (' '.join(str(x) for x in mm), mu, S, S, n_o, n_m, nsup, dev, int(sigma), 'MATCH' if (nsup == 0 and dev < 1e-9) else 'MISMATCH'))
    if cmd in ('extras', 'all'):
        print('== extras: cost of the lazily built non-adjacent generators and of all 12 float matrices (not part of t_a / t_b)')
        for S in (4, 8, 10):
            te, tf = time_all_generators(S)
            print('  [%d,%d,0,0]: exact E_13,E_24,E_14,E_31,E_42,E_41 %.4f s;  metric + all 12 float CSR matrices %.4f s' % (S, S, te, tf), flush=True)
    if cmd in ('n5', 'all'):
        print('== generality: N = 5 smoke test (PRV = RS, float + exact recursion, negative tests, brute force)')
        smoke_n(N5_PRODUCTS, 'N = 5')


if __name__ == '__main__':
    main(sys.argv)
