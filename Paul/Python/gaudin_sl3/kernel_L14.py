"""kernel_L14.py -- agent K, round 6: the one-row master kernel at lambda_1 = 8, 14 by the operator route scaled up, EXACT.

Spec: SPEC-kernel-L14.md; results: RESULTS-kernel-L14.md.  Python only (numpy/scipy; sympy/mpmath only in the one-site exact part and the checks).
Run from Paul/Python/gaudin_sl3 with `python -X utf8 kernel_L14.py <command> ...`:
    M <L> [batch]                       exact integer matrix M^(L)(1) -> kernel_M_L<L>.npz
    N <L> [nprimes|- [nofloat] [nosave]] exact N = M^-1 (multi-modular) + float64 comparison -> kernel_N_L<L>.npz
    Theta <L> <l2> [nprimes|- [nosave]]  exact Theta^{[L,l2,0]} (block theorem) + exchange / top-state sum rule in F_p -> kernel_Theta_<L><l2>0.npz
    univ <L...>                          universality M^(L)|labels<=L'+1 == M^(L')
    accept [Lbig]                        the stored exact kernels (acceptance (i)) from M^(L') and from M^(Lbig) restricted

Pipeline (one-row Lambda = [L,0,0], two sites; sections of this file):
  A. operators in the MONOMIAL basis (E_ij = x_i d_j, integer matrices) at zeta = 0 (g = E21 + E32), th1 - th2 = 1 (t = 1, delta = +1 / -1):
     H_s, b_s = Laurent coefficients of the quantum characteristic determinant  (verified against jet_basis_sym.companion_ops, exactly)
  B. Frobenius coefficients applied to covectors as a vector recursion (never as a polynomial; verified against jet_basis_sym.nu_formulas)
  C. y (one-site jets, exact Fractions), Y(1) rows, M(1) = Y(1) (y x y)^-1: an INTEGER unit block-lower-triangular matrix (level rho), computed in exact modular
     arithmetic (primes < 2^28, 3-prime CRT + verification prime); plain float64 fails from L = 8
  D. N = M^-1 exactly: block forward substitution modulo primes p < 790000 with float64 BLAS (exact integer sums < 2^53), Garner lift (symmetric digits)
  E. Theta^{[L,l2,0]} = block theorem applied to N[S x S, S x S]: exact residues -> Garner -> float64 + exact digits; T-format / npz
  F. checks: exchange rule and top-state sum rule in F_p; numeric checks of a stored file (mpmath), first-order closed form (C4), corner closed form (C3);
     exact-integer certificate of N columns
"""
import sys, os, time, pickle, itertools, math
from fractions import Fraction as Fr
from math import comb, factorial
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
os.environ.setdefault('OPENBLAS_NUM_THREADS', '8')      # 8 threads are as fast as 20 here, and small-m dgemm with 20 threads is pathological (1 GF/s)
import numpy as np
import scipy.sparse as sps

PMAX = 3          # highest Laurent power kept in the operator expansion (2 is enough; 3 for safety)


# =====================================================================================================  A. operators
def monos(L):
    """monomial exponents (a1,a2,a3), |a| = L, in the order of jet_basis_sym.monos (x1^L first)"""
    return [(a, b, L - a - b) for a in range(L, -1, -1) for b in range(L - a, -1, -1)]


def pairs_L(L):
    """jet pairs (k,l), 0 <= k < l <= L+1, lexicographic (== general_kernel.pairs_of((L,0,0)))"""
    return [(k, l) for k in range(L + 2) for l in range(k + 1, L + 2)]


def single_E(L, dtype=np.int64):
    """E_ij = x_i d_j on C[x]_L, column convention (E x^alpha = alpha_j x^{alpha+e_i-e_j}); dict (i,j) -> csr (d x d)"""
    ms = monos(L)
    ix = {m: k for k, m in enumerate(ms)}
    d = len(ms)
    E = {}
    for i in (1, 2, 3):
        for j in (1, 2, 3):
            rows, cols, vals = [], [], []
            for s, m in enumerate(ms):
                if m[j - 1] == 0:
                    continue
                t = list(m)
                c = t[j - 1]
                t[j - 1] -= 1
                t[i - 1] += 1
                rows.append(ix[tuple(t)]); cols.append(s); vals.append(c)
            E[(i, j)] = sps.csr_matrix((np.array(vals, dtype=dtype), (rows, cols)), shape=(d, d))
    return E, ms


class LS:
    """Laurent series with sparse-matrix coefficients: dict power -> csr (product = ordered matrix product, truncated above PMAX)"""
    def __init__(self, c=None):
        self.c = c if c is not None else {}

    def __add__(self, o):
        out = dict(self.c)
        for p, m in o.c.items():
            out[p] = out[p] + m if p in out else m
        return LS(out)

    def scale(self, s):
        return LS({p: m * s for p, m in self.c.items()})

    def __mul__(self, o):
        out = {}
        for p, a in self.c.items():
            for q, b in o.c.items():
                if p + q > PMAX:
                    continue
                m = a @ b
                out[p + q] = out[p + q] + m if (p + q) in out else m
        return LS(out)

    def deriv(self):
        return LS({p - 1: m * p for p, m in self.c.items() if p != 0})


def ell_series(E_s, E_t, g, delta, n, dtype):
    """dict (i,j) -> LS of l_ij(u) = g_ij + E^(s)_ij/x + E^(t)_ij/(x+delta), x = u - theta_s"""
    I = sps.identity(n, dtype=dtype, format='csr')
    out = {}
    for i in (1, 2, 3):
        for j in (1, 2, 3):
            c = {-1: E_s[(i, j)]}
            for p in range(0, PMAX + 1):
                coef = (-1) ** p * delta ** (-(p + 1))
                c[p] = E_t[(i, j)] * coef
            if g[i - 1][j - 1] != 0:
                c[0] = c[0] + I * g[i - 1][j - 1]
            out[(i, j)] = LS(c)
    return out


def cdet_series(S, ell):
    """Laurent series of S = sum coef * prod l^{(k)}_{ij} (ordered products; S from gaudin_chain.symbolic_cdet)"""
    tot = None
    for mono, coef in S.items():
        el = None
        for (i, j, k) in mono:
            s = ell[(i, j)]
            for _ in range(k):
                s = s.deriv()
            el = s if el is None else el * s
        if el is None:
            continue
        el = el.scale(coef)
        tot = el if tot is None else tot + el
    return tot


G_NILP = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]     # companion twist at zeta = 0: g = E21 + E32


def build_ops(L, delta=1, g=G_NILP, twosite=True, dtype=np.int64):
    """operators H_s, b_s (and the vanishing a_s, [x^-2]S2, [x^-3]S3) at both poles; two-site (d^2 x d^2) or one-site (d x d, key 'H1','b1' only).
    delta = theta_1 - theta_2 (site 1), -delta at site 2."""
    import gaudin_chain as gch
    E, ms = single_E(L, dtype)
    d = len(ms)
    I = sps.identity(d, dtype=dtype, format='csr')
    if twosite:
        n = d * d
        E1 = {k: sps.kron(v, I, format='csr') for k, v in E.items()}
        E2 = {k: sps.kron(I, v, format='csr') for k, v in E.items()}
    else:
        n = d
        E1 = E
        E2 = {k: sps.csr_matrix((d, d), dtype=dtype) for k in E}
    S1, S2, S3 = gch.symbolic_cdet()
    res = {}
    for site, (Es, Et, dl) in enumerate([(E1, E2, delta), (E2, E1, -delta)], start=1):
        ell = ell_series(Es, Et, g, dl, n, dtype)
        s2 = cdet_series(S2, ell)
        s3 = cdet_series(S3, ell)
        res['H%d' % site] = s2.c.get(-1)
        res['c2%d' % site] = s2.c.get(-2)
        res['b%d' % site] = s3.c.get(-1)
        res['a%d' % site] = s3.c.get(-2)
        res['t%d' % site] = s3.c.get(-3)
        if not twosite:
            break
    for k in ('c21', 'c22', 'a1', 'a2', 't1', 't2'):
        v = res.get(k)
        assert v is None or v.nnz == 0 or abs(v).max() <= 1e-9, ('operator %s does not vanish' % k)      # exact zero for integer dtype, rounding for floats
    return res, ms


# =====================================================================================================  B. Frobenius recursion on covectors
def frob_batch(X0, kind, HsT, BsT, HtT, BtT, delta, L, mod=None, t0=False):
    """Frobenius coefficients c_m applied (from the right) to the block of covectors X0 (columns), as a list c_0..c_mmax of arrays like X0.
    kind 'f0': rho = 0, c_1 = 0 (gamma_k = c_k, mmax = L+1); 'f1': rho = 1 (delta_l = c_{l-1}, mmax = L).  Operators are given TRANSPOSED
    (covector * A = A^T * covector).  delta integer +-1 (mod != None: exact arithmetic mod the prime `mod`, int64) or any number (float64).
    t0=True: decoupled limit delta -> infinity (the partner terms vanish; Ht, Bt ignored)."""
    rho = 0 if kind == 'f0' else 1
    mmax = L + 1 if kind == 'f0' else L
    nk = lambda k: rho + k

    def inv1(p):
        return (-1) ** p * delta ** (-(p + 1)) if not t0 else 0

    if mod is None:
        mul = lambda coef, a: float(coef) * a
        red = lambda a: a
        invI = lambda I: 1.0 / I
    else:
        mul = lambda coef, a: (int(coef) % mod) * a
        red = lambda a: a % mod
        invI = lambda I: pow(int(I) % mod, -1, mod)

    if mod is not None:
        assert t0 or abs(delta) == 1, 'the exact modular path needs delta = +-1 (th12 = +-1); other delta only in float64'
    zeros = np.zeros_like(X0)
    c = [X0]
    Hs_c = [red(HsT @ X0)]
    Bs_c = [red(BsT @ X0)]
    Ht_c = [red(HtT @ X0)] if not t0 else None
    Bt_c = [red(BtT @ X0)] if not t0 else None
    for m in range(1, mmax + 1):
        if kind == 'f0' and m == 1:
            cm = zeros
        else:
            I = nk(m) * (nk(m) - 1) * (nk(m) - 2 - L)
            rest = zeros.copy()
            rest = red(rest + mul(nk(m - 1), Hs_c[m - 1]))
            if m >= 2:
                rest = red(rest - Bs_c[m - 2])
            if not t0:
                for k in range(0, m):
                    cf = -L * inv1(m - k - 1) * nk(k) * (nk(k) - 1)
                    if cf != 0:
                        rest = red(rest + mul(cf, c[k]))
                for k in range(0, m - 1):
                    cf = inv1(m - k - 2) * nk(k)
                    if cf != 0:
                        rest = red(rest + mul(cf, Ht_c[k]))
                for k in range(0, m - 2):
                    cf = -inv1(m - k - 3)
                    rest = red(rest + mul(cf, Bt_c[k]))
            # c_m I = -rest   (the operator terms were collected as 'rest' with the sign of the equation  c_m I + rest = 0)
            cm = red(mul(-1, rest) if mod is not None else -rest)
            cm = red(mul(invI(I), cm)) if mod is not None else cm * invI(I)
        c.append(cm)
        if m < mmax:
            Hs_c.append(red(HsT @ cm))
            Bs_c.append(red(BsT @ cm))
            if not t0:
                Ht_c.append(red(HtT @ cm))
                Bt_c.append(red(BtT @ cm))
    return c


def jets_from_runs(Xs, kind):
    """list gamma_k (kind f0, k = 0..L+1) or delta_l (f1, l = 0..L+1) of arrays from the Frobenius coefficient list"""
    if kind == 'f0':
        return [Xs[0], np.zeros_like(Xs[0])] + Xs[2:]
    return [np.zeros_like(Xs[0])] + Xs


# =====================================================================================================  C. y, Z = y^-1, Y(1), M(1)
def _dense_obj(A):
    """sparse integer matrix -> dense object array of Python ints (exact arithmetic for the one-site part)"""
    D = np.zeros(A.shape, dtype=object)
    D[:] = 0
    Ac = A.tocoo()
    for r, c, v in zip(Ac.row, Ac.col, Ac.data):
        D[r, c] = int(v)
    return D


def frob_exact(X0, kind, h, b, L):
    """one-site (t = 0) Frobenius recursion in exact Fractions: h, b dense object arrays (NOT transposed: applied as covector * A = (A^T x)^T);
    X0 object array (d, B)"""
    rho = 0 if kind == 'f0' else 1
    mmax = L + 1 if kind == 'f0' else L
    nk = lambda k: rho + k
    hT, bT = h.T.copy(), b.T.copy()
    zeros = np.zeros_like(X0)
    zeros[:] = Fr(0)
    c = [X0]
    Hs_c = [hT.dot(X0)]
    Bs_c = [bT.dot(X0)]
    for m in range(1, mmax + 1):
        if kind == 'f0' and m == 1:
            cm = zeros
        else:
            I = nk(m) * (nk(m) - 1) * (nk(m) - 2 - L)
            rest = Hs_c[m - 1] * nk(m - 1)
            if m >= 2:
                rest = rest - Bs_c[m - 2]
            cm = (-rest) * Fr(1, I)
        c.append(cm)
        if m < mmax:
            Hs_c.append(hT.dot(cm))
            Bs_c.append(bT.dot(cm))
    return c


def y_exact(L):
    """one-site jet matrix y (rows = jet pairs (k,l) lexicographic, columns = monomials in `monos(L)` order), exact Fractions (object array)"""
    res, ms = build_ops(L, twosite=False)
    h, b = _dense_obj(res['H1']), _dense_obj(res['b1'])
    d = len(ms)
    x0 = np.zeros((d, 1), dtype=object)
    x0[:] = Fr(0)
    x0[0, 0] = Fr(1)                                  # coefficient of x1^L
    D = jets_from_runs(frob_exact(x0, 'f1', h, b, L), 'f1')     # D[l] (d,1), l = 0..L+1
    V = np.zeros((d, L + 2), dtype=object)
    for l in range(L + 2):
        V[:, l] = D[l][:, 0]
    G = jets_from_runs(frob_exact(V, 'f0', h, b, L), 'f0')      # G[k] (d, L+2): column l = x0 delta_l gamma_k
    prs = pairs_L(L)
    y = np.zeros((len(prs), d), dtype=object)
    for i, (k, l) in enumerate(prs):
        y[i, :] = (G[k][:, l] - G[l][:, k]) * (factorial(k) * factorial(l))
    return y, prs, ms


def weight_of_mono(m):
    return m[1] + 2 * m[2]


def y_inverse_exact(y, prs, ms):
    """Z = y^-1 (exact Fractions, dense object array); y is block diagonal w.r.t. the principal weight (asserted)"""
    d = len(ms)
    wr = [k + l - 1 for (k, l) in prs]
    wc = [weight_of_mono(m) for m in ms]
    Z = np.zeros((d, d), dtype=object)
    Z[:] = Fr(0)
    for i in range(d):
        for j in range(d):
            if wr[i] != wc[j] and y[i, j] != 0:
                raise AssertionError('y is not block diagonal in the principal weight')
    for w in sorted(set(wr)):
        ri = [i for i in range(d) if wr[i] == w]
        ci = [j for j in range(d) if wc[j] == w]
        assert len(ri) == len(ci), (w, len(ri), len(ci))
        n = len(ri)
        A = [[Fr(y[i, j]) for j in ci] for i in ri]
        Mx = [A[r] + [Fr(int(r == c)) for c in range(n)] for r in range(n)]
        for col in range(n):
            piv = next(r for r in range(col, n) if Mx[r][col] != 0)
            Mx[col], Mx[piv] = Mx[piv], Mx[col]
            pv = Mx[col][col]
            Mx[col] = [v / pv for v in Mx[col]]
            for r in range(n):
                if r != col and Mx[r][col] != 0:
                    f = Mx[r][col]
                    Mx[r] = [a - f * bb for a, bb in zip(Mx[r], Mx[col])]
        Ainv = [row[n:] for row in Mx]      # A^-1: Ainv[a][b], a = column-monomial index, b = row-pair index
        for a_, j in enumerate(ci):
            for b_, i in enumerate(ri):
                Z[j, i] = Ainv[a_][b_]
    return Z


def mod_matrix(Zexact, p):
    """exact Fraction matrix -> csr int64 residues mod p (entries in [0, p))"""
    n, m = Zexact.shape
    rows, cols, vals = [], [], []
    for i in range(n):
        for j in range(m):
            v = Zexact[i, j]
            if v != 0:
                v = Fr(v)
                den = v.denominator % p
                assert den != 0, 'denominator divisible by the prime'
                rows.append(i); cols.append(j); vals.append((v.numerator % p) * pow(den, -1, p) % p)
    return sps.csr_matrix((np.array(vals, dtype=np.int64), (rows, cols)), shape=(n, m))


def Zfloat(Zexact):
    return np.array([[float(Fr(v)) for v in row] for row in Zexact], dtype=np.float64)


def ops_T(res, mod=None):
    """transposed operators needed by the recursion, csr, float64 or int64"""
    out = {}
    for k in ('H1', 'H2', 'b1', 'b2'):
        A = res[k].T.tocsr()
        out[k] = A.astype(np.float64) if mod is None else A.astype(np.int64)
    return out


def build_Y_rows(L, opsT, mod=None, verbose=False, cols_sel=None, delta=1):
    """Y(1) rows as an array Y[i, j, :] (i, j = pair indices of site 1, site 2; last axis = ket components x^alpha y^beta, index alpha*d + beta),
    int64 residues mod `mod` or float64.  cols_sel: optional list of site-1 pair indices to compute (default all)."""
    prs = pairs_L(L)
    npr = len(prs)
    d = npr
    n = d * d
    dt = np.float64 if mod is None else np.int64
    H1, H2, b1, b2 = opsT['H1'], opsT['H2'], opsT['b1'], opsT['b2']
    fact = [factorial(k) for k in range(L + 2)]
    # ---- site 1 (delta = +1): x0 -> D runs, then gamma runs on all D_l
    x0 = np.zeros((n, 1), dtype=dt)
    x0[0, 0] = 1
    D = jets_from_runs(frob_batch(x0, 'f1', H1, b1, H2, b2, delta, L, mod), 'f1')
    V = np.concatenate([D[l] for l in range(L + 2)], axis=1)           # (n, L+2)
    G = jets_from_runs(frob_batch(V, 'f0', H1, b1, H2, b2, delta, L, mod), 'f0')   # G[k][:, l] = x0 D_l G_k
    W1 = np.zeros((n, npr), dtype=dt)
    for i, (k, l) in enumerate(prs):
        w = (G[k][:, l] - G[l][:, k])
        f = fact[k] * fact[l]
        W1[:, i] = (w * (f % mod)) % mod if mod else w * f
    sel = list(range(npr)) if cols_sel is None else list(cols_sel)
    W1 = np.ascontiguousarray(W1[:, sel])
    B = len(sel)
    # ---- site 2 (delta = -1): D runs on the batch, then gamma runs on each D_l'
    D2 = jets_from_runs(frob_batch(W1, 'f1', H2, b2, H1, b1, -delta, L, mod), 'f1')
    Y = np.zeros((B, npr, n), dtype=dt)
    pindex = {p: j for j, p in enumerate(prs)}
    for lp in range(L + 2):
        if verbose:
            print('   site-2 run', lp, flush=True)
        G2 = jets_from_runs(frob_batch(D2[lp], 'f0', H2, b2, H1, b1, -delta, L, mod), 'f0')   # G2[k'] = (.. D2_l') gamma_k'
        for kp in range(L + 2):
            if kp == lp:
                continue
            a_, b_ = (kp, lp) if kp < lp else (lp, kp)
            sg = 1 if kp < lp else -1
            j = pindex[(a_, b_)]
            f = fact[a_] * fact[b_]
            if mod is None:
                Y[:, j, :] += sg * f * G2[kp].T
            else:
                Y[:, j, :] = (Y[:, j, :] + ((sg * f) % mod) * G2[kp].T) % mod
    return Y


def M_from_Y(Y, Zmat, mod=None):
    """M rows: M_r = Z^T Y_r Z (Y_r reshaped d x d): float64 (Zmat dense float) or exact mod p (Zmat csr residues, via sparse kron(Z,Z))"""
    B, npr, n = Y.shape
    d = npr
    R = B * npr
    Y2 = Y.reshape(R, n)
    if mod is None:
        Y3 = Y2.reshape(R, d, d)
        out = np.matmul(np.matmul(Zmat.T[None, :, :], Y3), Zmat[None, :, :])
        return out.reshape(B, npr, n)
    K = sps.kron(Zmat, Zmat, format='csr')
    K.data %= mod                      # kron multiplies residues (< p^2 < 2^63): reduce BEFORE the matrix product
    KT = K.T.tocsr()
    out = np.zeros((R, n), dtype=np.int64)
    ch = 2048
    for s in range(0, R, ch):
        blk = Y2[s:s + ch]
        out[s:s + ch] = (KT @ blk.T).T % mod
    return out.reshape(B, npr, n)


# =====================================================================================================  C2. exact M(1) by several primes + CRT
def primes_below(start, count):
    import sympy
    out = []
    p = start
    while len(out) < count:
        p = sympy.prevprime(p)
        out.append(int(p))
    return out


PRIMES = primes_below(2 ** 28, 5)     # p < 2^28: products < 2^56, sums of <= 64 products < 2^62 (int64 safe)


def crt3_signed(r1, r2, r3, p1, p2, p3):
    """CRT of three residue arrays (primes < 2^28) -> signed int64 array, valid when |x| < 2^62 (the true range is p1 p2 p3 / 2 ~ 2^83;
    the low 64 bits are formed with uint64 wrap-around and the sign is decided from a float64 estimate, so no entry may be within 1e-6 relative of p1p2p3/2)"""
    t2 = ((r2 - r1) % p2) * pow(p1, -1, p2) % p2
    x12 = r1 + p1 * t2                                            # in [0, p1 p2), < 2^56
    t3 = ((r3 - x12) % p3) * pow((p1 * p2) % p3, -1, p3) % p3
    m = p1 * p2 * p3
    vf = x12.astype(np.float64) + float(p1 * p2) * t3.astype(np.float64)
    neg = vf > float(m) / 2.0
    P12 = np.uint64(p1 * p2)
    mm = np.uint64(m % (1 << 64))
    with np.errstate(over='ignore'):
        val = x12.astype(np.uint64) + P12 * t3.astype(np.uint64)
        val = np.where(neg, val - mm, val)
    return val.view(np.int64)


def build_M_exact(L, primes=None, batch=None, verbose=True, check_prime=True):
    """exact integer M^(L)(1): returns scipy csr int64 (R x R), R = (C(L+2,2))^2, rows r = i*d + j.  The first three primes lift (|M| < 2^62 required),
    the fourth verifies every entry (an entry beyond the range would fail the check)."""
    primes = primes or PRIMES[:4]
    assert len(primes) >= 4 or not check_prime
    res, ms = build_ops(L)
    y, prs, ms2 = y_exact(L)
    Zex = y_inverse_exact(y, prs, ms2)
    d = len(prs)
    R = d * d
    n = R
    batch = batch or d
    pm = {}
    for p in primes:
        pm[p] = (ops_T(res, p), mod_matrix(Zex, p))
    chunks = []
    t0 = time.time()
    for s in range(0, d, batch):
        sel = list(range(s, min(d, s + batch)))
        resid = []
        for p in primes:
            opsT, Zm = pm[p]
            Y = build_Y_rows(L, opsT, p, cols_sel=sel)
            Mp = M_from_Y(Y, Zm, p)
            del Y
            resid.append(Mp)
        x = crt3_signed(resid[0], resid[1], resid[2], primes[0], primes[1], primes[2])
        if check_prime and len(primes) >= 4:
            ok = np.array_equal(x % primes[3], resid[3])
            if not ok:
                raise RuntimeError('CRT verification failed in batch starting at %d (|M| beyond the 3-prime range 2^62?)' % s)
        chunks.append(sps.csr_matrix(x.reshape(len(sel) * d, n)))
        if verbose:
            print('  L=%d batch %d..%d done, %.1fs' % (L, sel[0], sel[-1], time.time() - t0), flush=True)
    M = sps.vstack(chunks, format='csr')
    M.eliminate_zeros()
    return M, prs


# =====================================================================================================  utilities: memory, timing
def mem_info():
    """(current working set, peak working set) in bytes of this process (Windows ctypes; falls back to resource on POSIX)"""
    try:
        import ctypes
        from ctypes import wintypes

        class PMC(ctypes.Structure):
            _fields_ = [('cb', wintypes.DWORD), ('PageFaultCount', wintypes.DWORD), ('PeakWorkingSetSize', ctypes.c_size_t),
                        ('WorkingSetSize', ctypes.c_size_t), ('QuotaPeakPagedPoolUsage', ctypes.c_size_t),
                        ('QuotaPagedPoolUsage', ctypes.c_size_t), ('QuotaPeakNonPagedPoolUsage', ctypes.c_size_t),
                        ('QuotaNonPagedPoolUsage', ctypes.c_size_t), ('PagefileUsage', ctypes.c_size_t),
                        ('PeakPagefileUsage', ctypes.c_size_t)]
        pmc = PMC()
        pmc.cb = ctypes.sizeof(PMC)
        k32 = ctypes.WinDLL('kernel32')
        k32.GetCurrentProcess.restype = wintypes.HANDLE
        k32.K32GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(PMC), wintypes.DWORD]
        k32.K32GetProcessMemoryInfo.restype = wintypes.BOOL
        h = k32.GetCurrentProcess()
        if not k32.K32GetProcessMemoryInfo(h, ctypes.byref(pmc), pmc.cb):
            raise OSError('GetProcessMemoryInfo failed')
        return pmc.WorkingSetSize, pmc.PeakWorkingSetSize
    except Exception:
        try:
            import resource
            r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
            return r, r
        except Exception:
            return 0, 0


# =====================================================================================================  D. N = M^-1 (float64 forward substitution; exact multi-modular)
def rho_vector(L):
    """jet level rho_r = (k+l-1) + (k'+l'-1) of the rows r = i*d + j (two sites), int64"""
    prs = pairs_L(L)
    lev = np.array([k + l - 1 for (k, l) in prs], dtype=np.int64)
    return (lev[:, None] + lev[None, :]).reshape(-1)


def level_perm(rho):
    perm = np.argsort(rho, kind='stable')
    rs = rho[perm]
    levels = np.unique(rs)
    starts = np.searchsorted(rs, levels, side='left')
    ends = np.searchsorted(rs, levels, side='right')
    return perm, levels, starts, ends


def check_unit_triangular(M, rho):
    """M must be unit block-lower-triangular in the level: M_{rc} != 0 only if rho_r > rho_c, or r == c"""
    Mc = M.tocoo()
    dr = rho[Mc.row] - rho[Mc.col]
    assert np.all((dr > 0) | ((Mc.row == Mc.col) & (Mc.data == 1))), 'M is not unit block-lower-triangular in the level'
    nd = np.count_nonzero(Mc.row == Mc.col)
    assert nd == M.shape[0], 'diagonal of M is not the identity'
    return int(dr.max())


def N_float(M, rho, maxdiff=None):
    """N = M^-1 in float64 (dense, ROW/COLUMN ORDER SORTED BY LEVEL, returned with perm: N_sorted[a,b] = N[perm[a], perm[b]]) by block forward substitution
    N_{rho,:} = I_rho - M[rho, window] N[window, :]  (M entries are exact integers < 2^53 here; the arithmetic is plain float64, NOT exact for large L)"""
    n = M.shape[0]
    md = check_unit_triangular(M, rho)
    maxdiff = md if maxdiff is None else maxdiff
    perm, levels, starts, ends = level_perm(rho)
    Ms = M[perm][:, perm].tocsr()
    N = np.zeros((n, n))
    for li, (a, b) in enumerate(zip(starts, ends)):
        N[a:b, a:b] = np.eye(b - a)
        if a == 0:
            continue
        lv = levels[li]                                # window: levels >= lv - maxdiff
        j0 = int(np.searchsorted(levels, lv - maxdiff, side='left'))
        a0 = starts[j0]
        Mb = Ms[a:b, a0:a].toarray()
        N[a:b, :a] = -(Mb @ N[a0:a, :a])
    return N, perm


class ModSolver:
    """exact N[:, cols] mod p by block forward substitution with float64 BLAS (all intermediate sums are integers < 2^53 when p < ~7.9e5):
    M must be unit block-lower-triangular in the level.  The dense level blocks M[rho, window] mod p are prepared once per prime (stored int32)."""
    def __init__(self, M, rho):
        self.n = M.shape[0]
        self.maxdiff = check_unit_triangular(M, rho)
        self.perm, self.levels, self.starts, self.ends = level_perm(rho)
        self.inv = np.empty(self.n, dtype=np.int64)
        self.inv[self.perm] = np.arange(self.n)
        self.Ms = M[self.perm][:, self.perm].tocsr()
        self.prep = {}
        self.wmax = 0
        for li, (a, b) in enumerate(zip(self.starts, self.ends)):
            if a > 0:
                j0 = int(np.searchsorted(self.levels, self.levels[li] - self.maxdiff, side='left'))
                self.wmax = max(self.wmax, a - self.starts[j0])

    def prepare(self, p):
        if p in self.prep:
            return
        assert (self.wmax + 1) * p * p < 2 ** 53, 'prime too large for exact float64 GEMM (window %d)' % self.wmax
        blocks = []
        for li, (a, b) in enumerate(zip(self.starts, self.ends)):
            if a == 0:
                blocks.append(None)
                continue
            j0 = int(np.searchsorted(self.levels, self.levels[li] - self.maxdiff, side='left'))
            a0 = self.starts[j0]
            Mb = (self.Ms[a:b][:, a0:a].toarray() % p).astype(np.int32)
            blocks.append((a0, np.ascontiguousarray(Mb.T)))          # transposed: (window x n_rho), contiguous
        self.prep[p] = blocks

    def release(self, p):
        self.prep.pop(p, None)

    def solve(self, cols_orig, p):
        """residues of N[:, cols_orig] (original column indices) as a float64 array of shape (len(cols), n): row c = column cols[c] of N, entries in
        ORIGINAL row order.  (Transposed layout: the GEMM is (nc x window) @ (window x n_rho); OpenBLAS is pathologically slow for small-m products.)"""
        blocks = self.prep[p]
        cols_s = self.inv[np.asarray(cols_orig)]          # positions in the sorted order
        nc = len(cols_s)
        Xt = np.zeros((nc, self.n))
        order = np.argsort(cols_s, kind='stable')
        cs_sorted = cols_s[order]
        for li, (a, b) in enumerate(zip(self.starts, self.ends)):
            lo = np.searchsorted(cs_sorted, a, side='left')
            hi = np.searchsorted(cs_sorted, b, side='left')
            rhs = np.zeros((nc, b - a))
            if hi > lo:
                rhs[order[lo:hi], cs_sorted[lo:hi] - a] = 1.0
            if blocks[li] is not None:
                a0, MbT = blocks[li]
                rhs -= Xt[:, a0:a] @ MbT.astype(np.float64)
            Xt[:, a:b] = np.mod(rhs, p)
        out = np.empty_like(Xt)
        out[:, self.perm] = Xt                               # sorted -> original row order
        return out


def garner_signed(res_list, primes):
    """mixed-radix (Garner) reconstruction with SYMMETRIC digits of the signed integer x (|x| < prod(primes)/2) from its residues (int64/float arrays);
    returns (digits list, float64 value of x).  x = sum_k digits[k] * prod_{j<k} primes[j]  exactly."""
    K = len(primes)
    digits = []
    P = [1]
    for k in range(K):
        pk = primes[k]
        r = np.asarray(res_list[k]).astype(np.int64)
        if k == 0:
            v = r % pk
        else:
            t = np.zeros_like(r)
            for j in range(k):
                t += digits[j] * (P[j] % pk)          # |digit| < p/2 < 4e5, P mod pk < 8e5: sums of <= 20 terms stay far below 2^63
            v = ((r - t) % pk) * pow(P[k] % pk, -1, pk) % pk
        v = np.where(v > pk // 2, v - pk, v)
        digits.append(v)
        P.append(P[-1] * pk)
    xf = digits[K - 1].astype(np.float64)
    for k in range(K - 2, -1, -1):
        xf = xf * float(primes[k]) + digits[k]
    return digits, xf


def primes_for_blas(count, upper=790000):
    """primes just below `upper` (p^2 * 14400 < 2^53 requires p < 790,900)"""
    import sympy
    out = []
    p = upper
    while len(out) < count:
        p = int(sympy.prevprime(p))
        out.append(p)
    return out


def N_exact_columns(solver, cols_orig, primes, chunk=1000, verbose=False):
    """exact N[:, cols] as float64 of shape (n, len(cols)) in original row order (Garner lift on the union support of the residues, |N| < prod(primes)/2 required:
    'depth' = number of primes actually needed, must be <= len(primes) - 2 for a safe margin).  Blocks are prepared once per prime (kept int32), then column
    chunks are solved for all primes."""
    cols_orig = np.asarray(cols_orig)
    nc = len(cols_orig)
    n = solver.n
    out = np.zeros((n, nc))
    depth = 0
    for p in primes:
        solver.prepare(p)
    for s in range(0, nc, chunk):
        e = min(nc, s + chunk)
        res = [solver.solve(cols_orig[s:e], p) for p in primes]          # each (e-s, n)
        mask = res[0] != 0
        for r in res[1:]:
            mask |= (r != 0)
        vals = [r[mask] for r in res]
        del res
        digits, xf = garner_signed(vals, primes)
        for k in range(len(primes) - 1, -1, -1):
            if np.any(digits[k] != 0):
                depth = max(depth, k + 1)
                break
        blk = np.zeros(mask.shape)
        blk[mask] = xf
        out[:, s:e] = blk.T
        del digits, vals, blk, mask
        if verbose:
            print('   N chunk %d..%d done' % (s, e), flush=True)
    return out, depth


# =====================================================================================================  E. Theta from N (exact mod p), block theorem, T-format, npz
S0 = 360360                       # lcm(1..15): c_kl * S0 is an integer for l <= 15 (all labels of lambda_1 <= 14)


def ckl_int(L, k, l):
    """S0 * c_kl(L), c_kl = (-1)^{k+l-1} (l-k) C(L,l-1) C(L+1,k) / l  (zero for l > L+1); needs l | S0"""
    if l - 1 > L:
        return 0
    assert S0 % l == 0, 'S0 does not clear the denominator l = %d' % l
    return (-1) ** (k + l - 1) * (l - k) * comb(L, l - 1) * comb(L + 1, k) * (S0 // l)


def block_factor(L, l2):
    """[l2! (L-l2)! / L!]^2 as a Fraction"""
    return Fr(factorial(l2) * factorial(L - l2), factorial(L)) ** 2


def theta_exact(L, l2, M, primes, chunk=1000, verbose=True, rho=None):
    """Exact Theta^{[L,l2,0]} (block theorem) from the exact integer M^(L): for the non-forced pairs S = {(k,l): l >= l2+1},
    N[S x S, S x S] = (M[S x S, S x S])^-1 (label triangularity), C~[x,y,(p,q)] = sum_{rho_r=x, rho_s=y} S0^2-scaled c_r c_s N_{rs;pq}, exact mod each prime
    -> Garner.  Returns dict: S (list of pairs), Cmod (list over primes, float64 arrays (nx, nx, nc) of residues), C (float64 exact-lifted C~/S0^2, shape (nx,nx,nc)),
    depth (primes needed), digits (list over digit index, int64 arrays) ; columns c = iS*|S| + jS."""
    prs = pairs_L(L)
    d = len(prs)
    S = [i for i, (k, l) in enumerate(prs) if l >= l2 + 1]
    nS = len(S)
    U = np.array([i * d + j for i in S for j in S], dtype=np.int64)
    rho_full = rho_vector(L) if rho is None else rho
    MU = M[U][:, U].tocsr()
    rhoU = rho_full[U]
    solver = ModSolver(MU, rhoU)
    nx = 2 * L + 1
    lev = np.array([prs[i][0] + prs[i][1] - 1 for i in S])
    cols = np.arange(nS * nS)
    Cmod = []
    t0 = time.time()
    for p in primes:
        solver.prepare(p)
        A = np.zeros((nx, nS))
        for iS, i in enumerate(S):
            k, l = prs[i]
            A[lev[iS], iS] = ckl_int(L, k, l) % p
        Cp = np.zeros((nx, nx, nS * nS))
        for s in range(0, nS * nS, chunk):
            e = min(nS * nS, s + chunk)
            Xt = solver.solve(cols[s:e], p)                  # (e-s, nS*nS) residues, row index in U-order (= (iS, jS))
            X3 = Xt.reshape(e - s, nS, nS)                   # [c, r, s]
            T1 = np.mod(np.matmul(A, X3), p)                 # [c, x, s] = sum_r A[x, r] X[c, r, s]
            T2 = np.mod(np.matmul(T1, A.T), p)               # [c, x, y] = sum_s T1[c, x, s] A[y, s]
            Cp[:, :, s:e] = T2.transpose(1, 2, 0)
        Cmod.append(Cp)
        solver.release(p)
        if verbose:
            print('   prime %d done (%.1fs)' % (p, time.time() - t0), flush=True)
    digits, xf = garner_signed([c.reshape(-1) for c in Cmod], primes)
    depth = 0
    for k in range(len(primes) - 1, -1, -1):
        if np.any(digits[k] != 0):
            depth = k + 1
            break
    C = (xf / float(S0) ** 2).reshape(nx, nx, nS * nS)
    return dict(S=S, nS=nS, lev=lev, Cmod=Cmod, C=C, depth=depth, digits=digits, L=L, l2=l2, primes=list(primes))


def theta_to_arrays(th, rho_pairs=None):
    """Flatten the exact Theta into monomial arrays for the T-format / npz.  Entry column c = iS*nS + jS; monomial (n,p,q) with
    n = x+y-rho_i-rho_j, p = x - l2, q = y - l2 (block theorem shift), coefficient = block_factor * C[x,y,c] (float64)."""
    L, l2, nS, lev, C = th['L'], th['l2'], th['nS'], th['lev'], th['C']
    f = float(block_factor(L, l2))
    nx = C.shape[0]
    xs, ys, cs = np.nonzero(C)
    order = np.argsort(cs, kind='stable')
    xs, ys, cs = xs[order], ys[order], cs[order]
    iS, jS = cs // nS, cs % nS
    n = xs + ys - lev[iS] - lev[jS]
    p = xs - l2
    q = ys - l2
    assert n.min() >= 0, 'negative t-power'
    assert p.min() >= 0 and q.min() >= 0, 'block theorem shift produced negative a-power (min %d, %d)' % (p.min(), q.min())
    coef = C[xs, ys, cs] * f
    ptr = np.zeros(nS * nS + 1, dtype=np.int64)
    np.add.at(ptr, cs + 1, 1)
    ptr = np.cumsum(ptr)
    out = dict(n=n.astype(np.int16), p=p.astype(np.int16), q=q.astype(np.int16), coef=coef, ptr=ptr)
    if 'digits' in th and th['depth'] > 0:
        flat = (xs * nx + ys) * (nS * nS) + cs
        out['exact_digits'] = np.stack([th['digits'][k][flat] for k in range(th['depth'])]).astype(np.int32)      # (depth, nmon) Garner digits of C~
    return out


def theta_T_exact(th):
    """exact T-format dict {(i,j): {(n,p,q): Fraction}} (small L only: Python big ints from the Garner digits), (i,j) over the non-forced pairs"""
    L, l2, nS, lev = th['L'], th['l2'], th['nS'], th['lev']
    primes, digits = th['primes'], th['digits']
    f = block_factor(L, l2)
    nx = 2 * L + 1
    K = len(digits)
    shape = th['Cmod'][0].shape
    flat = [d_.reshape(shape) for d_ in digits]
    nz = np.nonzero(np.any([d_ != 0 for d_ in flat], axis=0))
    T = {}
    pp = [1]
    for p in primes:
        pp.append(pp[-1] * p)
    for x, y, c in zip(*nz):
        val = 0
        for k in range(K):
            dk = int(flat[k][x, y, c])
            if dk:
                val += dk * pp[k]
        iS, jS = divmod(int(c), nS)
        n = int(x + y - lev[iS] - lev[jS])
        T.setdefault((iS, jS), {})[(n, int(x) - l2, int(y) - l2)] = Fr(val, S0 * S0) * f
    for iS in range(nS):
        for jS in range(nS):
            T.setdefault((iS, jS), {})
    return T


# =====================================================================================================  F. checks in residue arithmetic: exchange rule, top-state sum rule
def exchange_rule_residues(th):
    """Theta_ij(a,a',t) = Theta_ji(a',a,-t):  C[x,y,(i,j)] = (-1)^(x+y-rho_i-rho_j) C[y,x,(j,i)] mod every prime (exact).  Returns (n_violations, n_checked)."""
    nS, lev = th['nS'], th['lev']
    nx = th['Cmod'][0].shape[0]
    idx = np.arange(nS * nS)
    iS, jS = idx // nS, idx % nS
    colperm = jS * nS + iS
    x = np.arange(nx)[:, None, None]
    y = np.arange(nx)[None, :, None]
    sign = np.where(((x + y - lev[iS][None, None, :] - lev[jS][None, None, :]) % 2) == 0, 1, -1)
    bad = 0
    tot = 0
    for p, C in zip(th['primes'], th['Cmod']):
        Ct = C.transpose(1, 0, 2)[:, :, colperm]
        diff = np.mod(C - sign * Ct, p)
        bad += int(np.count_nonzero(diff))
        tot += int(C.size)
    return bad, tot


def _polmul(a, b, p):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] = (out[i + j] + x * y) % p
    return out


def _poleval(a, u, p):
    r = 0
    for c in reversed(a):
        r = (r * u + c) % p
    return r


def _polder(a, p):
    return [(i * a[i]) % p for i in range(1, len(a))] or [0]


def top_state_site(l1, l2, z2, z3, th_both, th_site, p, pairs):
    """m_{kl}(th_site) mod p for pairs; p2 built from BOTH poles th_both = (th1, th2)"""
    d = (z3 - z2) % p
    dinv = pow(d, -1, p)
    Pi = [1]
    for th in th_both:
        for _ in range(l2):
            Pi = _polmul(Pi, [(-th) % p, 1], p)
    p2 = [0] * len(Pi)
    cur = Pi
    k = 0
    while True:
        coef = pow(dinv, k, p)
        for i, c in enumerate(cur):
            p2[i] = (p2[i] + coef * c) % p
        cur = _polder(cur, p)
        k += 1
        if all(c == 0 for c in cur) or k > len(Pi) + 2:
            break
    r = []
    cur = p2
    for k in range(l1 + 3):
        r.append(_poleval(cur, th_site, p))
        der = _polder(cur, p)
        n = max(len(der), len(cur))
        der = der + [0] * (n - len(der))
        cc = cur + [0] * (n - len(cur))
        cur = [(der[i] + z2 * cc[i]) % p for i in range(n)]
    out = {}
    for (k, l) in pairs:
        out[(k, l)] = (r[k] * pow(z3, l, p) - pow(z3, k, p) * r[l]) % p
    return out


def sum_rule_residues(th, seed=12345, trials=1):
    """top-state sum rule  sum_ij Theta_ij m_i^top(th1) m_j^top(th2) = M1 M2 R(a) R(a'),  R(a) = (1 - z2 a)^(l1-l2) (1 - z3 a)^l1,
    evaluated in F_p at random residues for every prime used in the CRT (exact identity of rational functions => holds mod p).  Returns list of (p, lhs-rhs mod p) (all 0)."""
    import random
    L, l2, nS, lev = th['L'], th['l2'], th['nS'], th['lev']
    prs = pairs_L(L)
    S = [prs[i] for i in th['S']]
    rnd = random.Random(seed)
    out = []
    f = block_factor(L, l2)
    for p, C in zip(th['primes'], th['Cmod']):
        for _ in range(trials):
            while True:
                z2, z3, t1, t2, a, ap = [rnd.randrange(2, p - 1) for _ in range(6)]
                if len({z2, z3}) == 2 and (t1 - t2) % p != 0 and (z3 - z2) % p != 0:
                    break
            t = pow((t1 - t2) % p, -1, p)
            m1 = top_state_site(L, l2, z2, z3, (t1, t2), t1, p, S)
            m2 = top_state_site(L, l2, z2, z3, (t1, t2), t2, p, S)
            M1, M2 = m1[(0, l2 + 1)], m2[(0, l2 + 1)]
            nx = C.shape[0]
            ainv = pow(a, -1, p)
            apinv = pow(ap, -1, p)
            # Wt[x,y] = a^(x-l2) a'^(y-l2) t^(x+y)
            Wt = np.zeros((nx, nx), dtype=object)
            for x in range(nx):
                for y in range(nx):
                    ex = pow(a, x - l2, p) if x >= l2 else pow(ainv, l2 - x, p)
                    ey = pow(ap, y - l2, p) if y >= l2 else pow(apinv, l2 - y, p)
                    Wt[x, y] = ex * ey % p * pow(t, x + y, p) % p
            Wf = np.array([[int(v) for v in row] for row in Wt], dtype=np.float64).reshape(-1)
            # per column: sum_{xy} C[x,y,c] Wt[x,y]  (exact in float64: nx^2 p^2 < 2^53)
            Cf = C.reshape(nx * nx, nS * nS)
            col = np.mod(Cf.T @ Wf, p)
            tinv = pow(t, -1, p)
            total = 0
            for c in range(nS * nS):
                iS, jS = divmod(c, nS)
                v = int(col[c])
                if v == 0:
                    continue
                v = v * pow(tinv, int(lev[iS] + lev[jS]), p) % p
                total = (total + v * m1[S[iS]] % p * m2[S[jS]]) % p
            # prefactor f / S0^2 mod p
            pref = (f.numerator % p) * pow(f.denominator % p, -1, p) % p * pow(S0 * S0 % p, -1, p) % p
            lhs = total * pref % p
            R = lambda x_: pow((1 - z2 * x_) % p, L - l2, p) * pow((1 - z3 * x_) % p, L, p) % p
            rhs = M1 * M2 % p * R(a) % p * R(ap) % p
            out.append((p, (lhs - rhs) % p))
    return out


# =====================================================================================================  G. files (npz) and loaders
def _csr_arrays(prefix, A):
    A = A.tocsr()
    return {prefix + '_indptr': A.indptr.astype(np.int64), prefix + '_indices': A.indices.astype(np.int32), prefix + '_data': A.data,
            prefix + '_shape': np.array(A.shape, dtype=np.int64)}


def _csr_from(d, prefix):
    return sps.csr_matrix((d[prefix + '_data'], d[prefix + '_indices'], d[prefix + '_indptr']), shape=tuple(int(x) for x in d[prefix + '_shape']))


def M_path(L):
    return os.path.join(HERE, 'kernel_M_L%d.npz' % L)


def save_M(L, M):
    np.savez(M_path(L), L=np.array(L), pairs=np.array(pairs_L(L), dtype=np.int16), rho=rho_vector(L).astype(np.int16), **_csr_arrays('M', M))


def load_M(L):
    d = np.load(M_path(L))
    return _csr_from(d, 'M')


def get_M(L, verbose=True, batch=None):
    """exact M^(L): from the cache kernel_M_L<L>.npz, else from the M arrays stored inside kernel_N_L<L>.npz, else built (and cached)"""
    if os.path.exists(M_path(L)):
        return load_M(L)
    npath = os.path.join(HERE, 'kernel_N_L%d.npz' % L)
    if os.path.exists(npath):
        return _csr_from(np.load(npath), 'M')
    batch = batch or max(1, 2400 // (comb(L + 2, 2)) * 1)
    M, prs = build_M_exact(L, verbose=verbose, batch=batch)
    save_M(L, M)
    return M


def save_N(path, L, N, M, meta):
    """N: dense float64 (n, n) in the ORIGINAL row order (r = i*d + j over pairs_L(L)); stored as CSR together with the exact integer M (CSR int64)"""
    arr = dict(L=np.array(L), pairs=np.array(pairs_L(L), dtype=np.int16), rho=rho_vector(L).astype(np.int16))
    arr.update(_csr_arrays('N', sps.csr_matrix(N)))
    arr.update(_csr_arrays('M', M))
    arr.update({k: np.array(v) for k, v in meta.items()})
    np.savez_compressed(path, **arr)


def load_N(path):
    """returns dict(L, pairs, rho, N (csr float64, lower block-triangular in the level), M (csr int64, exact), meta...)"""
    d = np.load(path)
    out = {k: d[k] for k in d.files if not (k.startswith('N_') or k.startswith('M_'))}
    out['N'] = _csr_from(d, 'N')
    out['M'] = _csr_from(d, 'M')
    return out


def save_Theta(path, th, extra=None):
    ar = theta_to_arrays(th)
    prs = pairs_L(th['L'])
    pairs = np.array([prs[i] for i in th['S']], dtype=np.int16)
    f = block_factor(th['L'], th['l2'])
    meta = dict(lam=np.array([th['L'], th['l2'], 0]), pairs=pairs, entry_ptr=ar['ptr'], n=ar['n'], p=ar['p'], q=ar['q'], coef=ar['coef'],
                prefactor=np.array(float(f)), primes=np.array(th['primes'][:th['depth']]), depth=np.array(th['depth']),
                exact_digits=ar['exact_digits'],
                # exact coefficient of monomial m = Ctilde_m * exact_num / exact_den, Ctilde_m = sum_k exact_digits[k,m] * prod_{j<k} primes[j]  (Python ints)
                exact_num=np.array(str(f.numerator)), exact_den=np.array(str(S0 * S0 * f.denominator)))
    if extra:
        meta.update({k: np.array(v) for k, v in extra.items()})
    np.savez_compressed(path, **meta)


def load_Theta(path, as_T=True):
    """Theta file -> dict(lam, pairs, arrays...) and, if as_T, T[(i,j)] = {(n,p,q): float coefficient} in the format of general_kernel (pairs = non-forced pairs of lam)"""
    d = np.load(path)
    out = {k: d[k] for k in d.files}
    if as_T:
        nS = len(out['pairs'])
        ptr, n, p, q, coef = out['entry_ptr'], out['n'], out['p'], out['q'], out['coef']
        T = {}
        for e in range(nS * nS):
            sl = slice(ptr[e], ptr[e + 1])
            T[divmod(e, nS)] = {(int(a), int(b), int(c)): float(v) for a, b, c, v in zip(n[sl], p[sl], q[sl], coef[sl])}
        out['T'] = T
    return out


def theta_exact_numerators(arrs, sl=None):
    """exact integers Ctilde (Python ints, object array) of the monomials (all, or slice `sl`) of a loaded kernel_Theta file: Ctilde = sum_k digit_k prod_{j<k} p_j;
    the exact coefficient is Ctilde * int(exact_num) / int(exact_den)"""
    dg = arrs['exact_digits']
    primes = [int(x) for x in arrs['primes']]
    sl = slice(None) if sl is None else sl
    K = dg.shape[0]
    val = dg[K - 1][sl].astype(object)
    for k in range(K - 2, -1, -1):
        val = val * primes[k] + dg[k][sl].astype(object)
    return val


def theta_coefficients_exact(arrs, sl=None):
    """exact Fractions of the coefficients (all monomials or slice)"""
    num, den = int(str(arrs['exact_num'])), int(str(arrs['exact_den']))
    return [Fr(int(c) * num, den) for c in theta_exact_numerators(arrs, sl)]


def eval_Theta_entry_mp(arrs, e, a, ap, t, dps=60):
    """entry e = i*nS + j of a loaded Theta file at (a, a', t) in mpmath arithmetic with the EXACT coefficients"""
    import mpmath as mp
    mp.mp.dps = dps
    sl = slice(int(arrs['entry_ptr'][e]), int(arrs['entry_ptr'][e + 1]))
    num, den = int(str(arrs['exact_num'])), int(str(arrs['exact_den']))
    ints = theta_exact_numerators(arrs, sl)
    A, AP, T = mp.mpf(a), mp.mpf(ap), mp.mpf(t)
    tot = mp.mpf(0)
    for c, n_, p_, q_ in zip(ints, arrs['n'][sl], arrs['p'][sl], arrs['q'][sl]):
        tot += (mp.mpf(int(c)) * num / den) * T ** int(n_) * A ** int(p_) * AP ** int(q_)
    return tot


def eval_Theta_entry(arrs, e, a, ap, t):
    """float64 value of entry e = i*nS + j of a loaded Theta file at (a, a', t) (plain double-precision power sum)"""
    sl = slice(arrs['entry_ptr'][e], arrs['entry_ptr'][e + 1])
    return float(np.sum(arrs['coef'][sl] * t ** arrs['n'][sl].astype(float) * a ** arrs['p'][sl].astype(float) * ap ** arrs['q'][sl].astype(float)))


# =====================================================================================================  H. restriction (universality), acceptance tests, CLI
def restrict_M(M, Lbig, Lsmall):
    """M^(Lbig) restricted to the jet labels <= Lsmall+1 (rows/cols i*d+j over pairs_L(Lbig) -> pairs_L(Lsmall))"""
    prs = pairs_L(Lbig)
    d = len(prs)
    ix = {p: i for i, p in enumerate(prs)}
    sel = [ix[p] for p in pairs_L(Lsmall)]
    idx = np.array([i * d + j for i in sel for j in sel])
    return M[idx][:, idx].tocsr()


def restriction_index(Lbig, Lsmall):
    prs = pairs_L(Lbig)
    d = len(prs)
    ix = {p: i for i, p in enumerate(prs)}
    sel = [ix[p] for p in pairs_L(Lsmall)]
    return np.array([i * d + j for i in sel for j in sel])


def N6_reference():
    """exact N^(6) (master_kernel2_N6.pkl) as a dense float64 array over pairs_L(6)^2 (entries are integers; those above 2^53 are exactly representable)"""
    N6 = pickle.load(open(os.path.join(HERE, 'master_kernel2_N6.pkl'), 'rb'))
    d = len(pairs_L(6))
    A = np.zeros((d * d, d * d))
    for (r, c), v in N6.items():
        A[r, c] = float(int(v))
    return A


def fmt_gb(b):
    return '%.2f GB' % (b / 1e9)


def cmd_M(L, batch=None):
    t0 = time.time()
    d = comb(L + 2, 2)
    batch = batch or max(1, min(d, 2400 // d))
    M, prs = build_M_exact(L, verbose=True, batch=batch)
    dt = time.time() - t0
    save_M(L, M)
    cur, peak = mem_info()
    print('M^(%d): dim %d, nnz %d (fill %.4f), max|M| %d (%.3e), time %.1fs, peak working set %s, batch %d' % (
        L, d * d, M.nnz, M.nnz / (d * d) ** 2, abs(M).max(), abs(M).max(), dt, fmt_gb(peak), batch))
    return M


def cmd_N(L, nprimes=None, chunk=1000, float_compare=True, save=True):
    M = get_M(L)
    n = M.shape[0]
    d = comb(L + 2, 2)
    rho = rho_vector(L)
    nprimes = nprimes or (L + 2)
    primes = primes_for_blas(nprimes)
    t0 = time.time()
    solver = ModSolver(M, rho)
    Nex, depth = N_exact_columns(solver, np.arange(n), primes, chunk=chunk, verbose=True)
    t_exact = time.time() - t0
    cur, peak = mem_info()
    assert depth <= nprimes - 2, 'not enough primes: depth %d of %d' % (depth, nprimes)
    nnz = int(np.count_nonzero(Nex))
    print('N^(%d) exact (multi-modular, %d primes, depth %d): time %.1fs, peak %s, nnz %d (fill %.4f), max|N| %.4e' % (
        L, nprimes, depth, t_exact, fmt_gb(peak), nnz, nnz / n ** 2, np.abs(Nex).max()), flush=True)
    meta = dict(depth=depth, nprimes=nprimes, maxabsN=np.abs(Nex).max(), build_seconds=t_exact)
    if float_compare:
        t0 = time.time()
        Nf, perm = N_float(M, rho)
        t_float = time.time() - t0
        worst = 0.0
        worst_abs = 0.0
        for s in range(0, n, 600):
            e = min(n, s + 600)
            ref = Nex[perm[s:e]][:, perm]
            msk = ref != 0
            with np.errstate(all='ignore'):
                rel = np.abs(Nf[s:e] - ref)[msk] / np.abs(ref)[msk]
            if rel.size:
                worst = max(worst, float(np.nanmax(rel)) if np.any(np.isfinite(rel)) else float('inf'))
        print('plain float64 forward substitution from the exact M: %.1fs, max relative error over the nonzero entries of N: %.3e' % (t_float, worst), flush=True)
        meta['float_solve_max_rel_err'] = worst
        meta['float_solve_seconds'] = t_float
        del Nf
    # acceptance (ii): restriction to labels <= 7 vs the exact N^(6)
    if L >= 6:
        idx = restriction_index(L, 6)
        ref6 = N6_reference()
        sub = Nex[np.ix_(idx, idx)]
        msk = ref6 != 0
        dev = float(np.max(np.abs(sub - ref6)[msk] / np.abs(ref6)[msk]))
        spurious = int(np.count_nonzero((~msk) & (sub != 0)))
        print('N^(%d)|labels<=7 vs exact N^(6): max relative deviation %.3e, spurious nonzeros %d' % (L, dev, spurious), flush=True)
        meta['restriction6_max_rel_dev'] = dev
    if save:
        path = os.path.join(HERE, 'kernel_N_L%d.npz' % L)
        t0 = time.time()
        save_N(path, L, Nex, M, meta)
        print('saved %s (%.1fs, %s)' % (path, time.time() - t0, fmt_gb(os.path.getsize(path))), flush=True)
    return Nex


def theta_report(th, name):
    ar = theta_to_arrays(th)
    c = np.abs(ar['coef'])
    nS = th['nS']
    print('%s: %d pairs, %d entries, %d monomials, max|coef| %.4e, min|coef| %.4e, depth %d of %d primes, a/a\'-degree max %d, t-degree max %d' % (
        name, nS, nS * nS, len(ar['coef']), c.max(), c[c > 0].min(), th['depth'], len(th['primes']), int(max(ar['p'].max(), ar['q'].max())), int(ar['n'].max())), flush=True)
    return ar


def cmd_Theta(L, l2, nprimes=None, chunk=1000, trials=2, save=True):
    M = get_M(L)
    nprimes = nprimes or (L + 8)
    primes = primes_for_blas(nprimes)
    t0 = time.time()
    th = theta_exact(L, l2, M, primes, chunk=chunk, verbose=True)
    dt = time.time() - t0
    cur, peak = mem_info()
    assert th['depth'] <= nprimes - 2, 'not enough primes: depth %d of %d' % (th['depth'], nprimes)
    print('Theta^[%d,%d,0]: time %.1fs, peak %s' % (L, l2, dt, fmt_gb(peak)), flush=True)
    theta_report(th, 'Theta^[%d,%d,0]' % (L, l2))
    bad, tot = exchange_rule_residues(th)
    print('exchange rule Theta_ij(a,a\',t) = Theta_ji(a\',a,-t): %d violations among %d residues (all %d primes, exact)' % (bad, tot, len(primes)), flush=True)
    sr = sum_rule_residues(th, trials=trials)
    nb = sum(1 for (p, v) in sr if v != 0)
    print('top-state sum rule (random residues, %d primes x %d trials): %d failures' % (len(primes), trials, nb), flush=True)
    if save:
        path = os.path.join(HERE, 'kernel_Theta_%d%d0.npz' % (L, l2))
        save_Theta(path, th, extra=dict(exchange_violations=bad, sum_rule_failures=nb, build_seconds=dt))
        print('saved %s (%s)' % (path, fmt_gb(os.path.getsize(path))), flush=True)
    return th


def cmd_univ(Ls):
    Ms = {L: get_M(L, verbose=False) for L in Ls}
    ok = True
    for L in Ls:
        for Lp in Ls:
            if Lp < L:
                R = restrict_M(Ms[L], L, Lp)
                diff = (R - Ms[Lp])
                diff.eliminate_zeros()
                good = diff.nnz == 0
                ok &= good
                print('M^(%d)|labels<=%d == M^(%d): %s' % (L, Lp + 1, Lp, good))
    print('universality of M (exact integer equality) for L in %s: %s' % (list(Ls), ok))


def compare_T_exact(Tnew, Tref):
    bad = 0
    n = 0
    for k in set(Tnew) | set(Tref):
        a = {kk: v for kk, v in Tnew.get(k, {}).items() if v != 0}
        b = {kk: Fr(v) for kk, v in Tref.get(k, {}).items() if v != 0}
        n += 1
        if a != b:
            bad += 1
    return bad, n


def cmd_accept(Lbig=None):
    """acceptance (i): kernels from the pipeline vs the stored exact kernels (exact Fraction equality of every monomial coefficient), built (a) from the
    independently computed M^(L') and (b) from the big M^(Lbig) restricted to labels <= L'+1."""
    import general_kernel as G
    import master_kernel2 as M2
    primes = primes_for_blas(10)
    cases = [(1, 0), (2, 0), (3, 0), (4, 0), (5, 0), (6, 0), (1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3), (4, 4)]
    for src in ('M(L)', 'M(%s)|' % Lbig if Lbig else None):
        if src is None:
            continue
        for (l1, l2) in cases:
            if src == 'M(L)':
                Ml = get_M(l1, verbose=False)
            else:
                Ml = restrict_M(get_M(Lbig, verbose=False), Lbig, l1)
            th = theta_exact(l1, l2, Ml, primes, verbose=False)
            T = theta_T_exact(th)
            Tref = M2.get_T(l1) if (l2 == 0 and l1 >= 5) else G.load_rep((l1, l2, 0))[0]
            bad, n = compare_T_exact(T, Tref)
            # float relative deviation of every monomial coefficient
            ar = theta_to_arrays(th)
            worst = 0.0
            nS = th['nS']
            for e in range(nS * nS):
                i, j = divmod(e, nS)
                ref = {kk: float(v) for kk, v in Tref[(i, j)].items() if v != 0}
                got = {(int(a), int(b), int(c)): v for a, b, c, v in zip(ar['n'][ar['ptr'][e]:ar['ptr'][e + 1]], ar['p'][ar['ptr'][e]:ar['ptr'][e + 1]],
                                                                         ar['q'][ar['ptr'][e]:ar['ptr'][e + 1]], ar['coef'][ar['ptr'][e]:ar['ptr'][e + 1]])}
                assert set(ref) == set(got), ('monomial support differs', l1, l2, e)
                for kk in ref:
                    worst = max(worst, abs(got[kk] - ref[kk]) / abs(ref[kk]))
            print('%-10s [%d,%d,0]: %d entries, exact mismatches %d, max rel dev of float64 coefficients %.2e' % (src, l1, l2, n, bad, worst), flush=True)


# =====================================================================================================  I. numerical checks of the stored (double-precision) kernels
def check_theta_file(path, dps=60, seed=3, verbose=True):
    """independent numerical test of a kernel_Theta_*.npz file (the stored float64 coefficients are treated as exact):
    (1) exchange rule  coef_ij(n,p,q) = (-1)^n coef_ji(n,q,p)  (relative, per monomial);
    (2) top-state sum rule at a random rational point in mpmath with `dps` digits: sum_ij Theta_ij m_i^top(th1) m_j^top(th2) vs M1 M2 R(a) R(a')
        (reference jets from general_kernel.top_monos_exact, exact rationals)."""
    import random
    import mpmath as mp
    import general_kernel as G
    d = np.load(path)
    lam = tuple(int(x) for x in d['lam'])
    l1, l2 = lam[0], lam[1]
    pairs = [tuple(int(v) for v in r) for r in d['pairs']]
    nS = len(pairs)
    ptr, nn, pp, qq, coef = d['entry_ptr'], d['n'].astype(np.int64), d['p'].astype(np.int64), d['q'].astype(np.int64), d['coef']
    # (1) exchange rule via a sorted key array
    ent = np.repeat(np.arange(nS * nS), np.diff(ptr))
    i_, j_ = ent // nS, ent % nS
    base = 64
    key = ((ent * base + nn) * base + pp) * base + qq
    key_ex = (((j_ * nS + i_) * base + nn) * base + qq) * base + pp
    order = np.argsort(key)
    skey = key[order]
    pos = np.searchsorted(skey, key_ex)
    ok_support = bool(np.all((pos < len(skey)) & (skey[np.minimum(pos, len(skey) - 1)] == key_ex)))
    partner = coef[order[np.minimum(pos, len(skey) - 1)]]
    sgn = np.where(nn % 2 == 0, 1.0, -1.0)
    rel = np.abs(coef - sgn * partner) / np.maximum(np.abs(coef), 1e-300)
    res = dict(exchange_support_ok=ok_support, exchange_max_rel=float(rel.max()))
    # (2) sum rule in mpmath
    mp.mp.dps = dps
    rnd = random.Random(seed)
    # moderate random rational point (|a|,|a'| <= 1, |z|,|theta| ~ 1): the double-rounded coefficients reproduce the identity up to ~1.1e-16 * kappa, kappa = sum|terms|/|rhs|
    z2, z3 = Fr(rnd.randint(-9, 9), rnd.randint(5, 9)), Fr(rnd.randint(-9, 9), rnd.randint(5, 9))
    while z2 == z3:
        z3 = Fr(rnd.randint(-9, 9), rnd.randint(5, 9))
    th = (Fr(rnd.randint(-9, 9), rnd.randint(5, 9)), Fr(rnd.randint(-9, 9), rnd.randint(5, 9)))
    while th[0] == th[1]:
        th = (th[0], Fr(rnd.randint(-9, 9), rnd.randint(5, 9)))
    a, ap = Fr(rnd.randint(-10, 10) or 1, rnd.randint(10, 20)), Fr(rnd.randint(-10, 10) or 1, rnd.randint(10, 20))
    t = 1 / (th[0] - th[1])
    (m1, M1), (m2, M2) = G.top_monos_exact(lam, z2, z3, th)
    mpf = lambda fr: mp.mpf(fr.numerator) / fr.denominator
    A, AP, T = mpf(a), mpf(ap), mpf(t)
    maxn = int(nn.max()); maxp = int(pp.max()); maxq = int(qq.max())
    Tpow = [T ** k for k in range(maxn + 1)]
    Apow = [A ** k for k in range(maxp + 1)]
    Qpow = [AP ** k for k in range(maxq + 1)]
    mm1 = [mp.mpf(sp_.p) / sp_.q for sp_ in (m1[pr] for pr in pairs)]
    mm2 = [mp.mpf(sp_.p) / sp_.q for sp_ in (m2[pr] for pr in pairs)]
    lhs = mp.mpf(0)
    scale = mp.mpf(0)
    for e in range(nS * nS):
        s, f_ = ptr[e], ptr[e + 1]
        if f_ == s:
            continue
        i, j = divmod(e, nS)
        acc = mp.mpf(0)
        accabs = mp.mpf(0)
        for k in range(s, f_):
            term = mp.mpf(float(coef[k])) * Tpow[nn[k]] * Apow[pp[k]] * Qpow[qq[k]]
            acc += term
            accabs += abs(term)
        lhs += acc * mm1[i] * mm2[j]
        scale += accabs * abs(mm1[i] * mm2[j])
    R = lambda x: (1 - mpf(z2) * x) ** (l1 - l2) * (1 - mpf(z3) * x) ** l1
    rhs = (mp.mpf(M1.p) / M1.q) * (mp.mpf(M2.p) / M2.q) * R(A) * R(AP)
    res['sum_rule_rel_dev'] = float(abs(lhs - rhs) / abs(rhs))
    res['sum_rule_condition'] = float(scale / abs(rhs))        # sum |terms| / |rhs|: the deviation is expected to be <~ 1.1e-16 * this
    res['sum_rule_dev_over_u_kappa'] = res['sum_rule_rel_dev'] / (1.11e-16 * res['sum_rule_condition'])
    res['sum_rule_lhs'] = mp.nstr(lhs, 12)
    res['sum_rule_rhs'] = mp.nstr(rhs, 12)
    if verbose:
        print(path, res, flush=True)
    return res


def eval_kappa(path, nent=60, points=None, dps=60, seed=11):
    """digits lost when a kernel entry is evaluated in plain double precision (power sum of the stored float64 coefficients, numpy) at the given points
    (a, a', t): for random entries compare with the mpmath value computed from the EXACT coefficients (the Garner digits stored in the file);
    returns list of (point, median/max relative error, median/max kappa = sum|float terms| / |exact value|)"""
    import mpmath as mp
    mp.mp.dps = dps
    d = np.load(path)
    arrs = {k: d[k] for k in d.files}
    ptr, nn, pp, qq, coef = arrs['entry_ptr'], arrs['n'].astype(np.int64), arrs['p'].astype(np.int64), arrs['q'].astype(np.int64), arrs['coef']
    rnd = np.random.RandomState(seed)
    nonempty = np.nonzero(np.diff(ptr) > 0)[0]
    ents = rnd.choice(nonempty, size=min(nent, len(nonempty)), replace=False)
    out = []
    for pt in points:
        a, ap, t = pt
        errs, kaps = [], []
        for e in ents:
            sl = slice(ptr[e], ptr[e + 1])
            terms = coef[sl] * np.power(t, nn[sl].astype(float)) * np.power(a, pp[sl].astype(float)) * np.power(ap, qq[sl].astype(float))
            fl = float(np.sum(terms))
            ex = eval_Theta_entry_mp(arrs, int(e), a, ap, t, dps)
            kap = float(np.sum(np.abs(terms)) / abs(float(ex))) if ex != 0 else float('inf')
            err = float(abs(mp.mpf(fl) - ex) / abs(ex)) if ex != 0 else float('inf')
            errs.append(err)
            kaps.append(kap)
        out.append((pt, float(np.median(errs)), float(np.max(errs)), float(np.median(kaps)), float(np.max(kaps))))
    return out


def check_first_order(path, nent=400, seed=5, dps=40):
    """first-order term (t^1 coefficient polynomial in a, a') of a kernel_Theta file vs the independent closed form `general_kernel.theta1_closed`
    (conjecture C4 of G1, fitted on l1 <= 4) at random rational (a, a') for random entries; returns the worst relative deviation"""
    import random
    import mpmath as mp
    import general_kernel as G
    mp.mp.dps = dps
    d = np.load(path)
    lam = tuple(int(x) for x in d['lam'])
    pairs = [tuple(int(v) for v in r) for r in d['pairs']]
    nS = len(pairs)
    ptr, nn, pp, qq, coef = d['entry_ptr'], d['n'], d['p'], d['q'], d['coef']
    rnd = random.Random(seed)
    worst = 0.0
    for _ in range(nent):
        i, j = rnd.randrange(nS), rnd.randrange(nS)
        a, ap = Fr(rnd.randint(-9, 9) or 1, rnd.randint(2, 9)), Fr(rnd.randint(-9, 9) or 1, rnd.randint(2, 9))
        e = i * nS + j
        sl = slice(ptr[e], ptr[e + 1])
        val = mp.mpf(0)
        for n_, p_, q_, c in zip(nn[sl], pp[sl], qq[sl], coef[sl]):
            if n_ == 1:
                val += mp.mpf(float(c)) * (mp.mpf(a.numerator) / a.denominator) ** int(p_) * (mp.mpf(ap.numerator) / ap.denominator) ** int(q_)
        ref = G.theta1_closed(lam, pairs[i], pairs[j], a, ap, 1)
        refm = mp.mpf(ref.numerator) / ref.denominator
        if refm == 0:
            dev = float(abs(val))
        else:
            dev = float(abs(val - refm) / abs(refm))
        worst = max(worst, dev)
    return worst


def digits_to_ints(digits, primes):
    """exact Python-int object array from Garner digits (list of int64 arrays)"""
    K = len(primes)
    val = digits[K - 1].astype(object)
    for k in range(K - 2, -1, -1):
        val = val * int(primes[k]) + digits[k].astype(object)
    return val


def verify_exact_columns(M, rho, cols, primes, chunk=64):
    """independent exactness certificate: solve N[:, cols] mod all primes, rebuild the EXACT integer columns (Python ints, no modular arithmetic left),
    and check M * N[:, c] == e_c with Python-int arithmetic.  Returns (number of columns checked, number of failures, max |N| entry over these columns)."""
    solver = ModSolver(M, rho)
    for p in primes:
        solver.prepare(p)
    Mc = M.tocsr()
    fails = 0
    mx = 0
    cols = np.asarray(cols)
    for s in range(0, len(cols), chunk):
        cc = cols[s:s + chunk]
        res = [solver.solve(cc, p) for p in primes]
        digits, xf = garner_signed([r.reshape(-1) for r in res], primes)
        ints = digits_to_ints(digits, primes).reshape(len(cc), -1)         # (ncols, n)
        for ci, c in enumerate(cc):
            v = ints[ci]
            mx = max(mx, int(max(abs(x) for x in v)))
            out = np.zeros(M.shape[0], dtype=object)
            for r in range(M.shape[0]):
                a, b = Mc.indptr[r], Mc.indptr[r + 1]
                if b > a:
                    out[r] = int(np.dot(Mc.data[a:b].astype(object), v[Mc.indices[a:b]]))
                else:
                    out[r] = 0
            out[c] -= 1
            if any(x != 0 for x in out):
                fails += 1
    return len(cols), fails, mx


def check_corner(path):
    """corner block (pairs (k, l1+1)) of a kernel_Theta file vs the closed form of the single-index family (conjecture C3 of G1, verified for l <= 5):
    Theta^{[l1,l2,0]}_{(k,l1+1),(k',l1+1)} = [l2!(l1-l2)!/l1!]^2 (a a')^(l1-l2) Theta^{[l1,l1,0]}_{kk'},
    Theta^{[l,l,0]}_{kk'} = (-1)^{k+k'} C(l,k) C(l,k') a^k a'^k' sum_{m,m'} C(l-k,m) C(l-k',m') F_{mm'}(k,k') (a t)^m (-a' t)^m'.
    Compared monomial by monomial (support and relative value).  Returns (entries compared, monomials compared, worst relative deviation, support mismatches)."""
    import general_kernel as G
    d = np.load(path)
    l1, l2 = int(d['lam'][0]), int(d['lam'][1])
    pairs = [tuple(int(v) for v in r) for r in d['pairs']]
    nS = len(pairs)
    ptr, nn, pp, qq, coef = d['entry_ptr'], d['n'], d['p'], d['q'], d['coef']
    corner = [i for i, (k, l) in enumerate(pairs) if l == l1 + 1]
    f = block_factor(l1, l2)
    worst = 0.0
    nmon = 0
    mism = 0
    for i in corner:
        for j in corner:
            k, kp = pairs[i][0], pairs[j][0]
            exp = {}
            for m in range(l1 - k + 1):
                for mp_ in range(l1 - kp + 1):
                    c = (-1) ** (k + kp + mp_) * comb(l1, k) * comb(l1, kp) * comb(l1 - k, m) * comb(l1 - kp, mp_) * G.F_mm(m, mp_, k, kp) * f
                    if c != 0:
                        exp[(m + mp_, k + m + (l1 - l2), kp + mp_ + (l1 - l2))] = exp.get((m + mp_, k + m + (l1 - l2), kp + mp_ + (l1 - l2)), 0) + c
            e = i * nS + j
            sl = slice(ptr[e], ptr[e + 1])
            got = {(int(a), int(b), int(c)): float(v) for a, b, c, v in zip(nn[sl], pp[sl], qq[sl], coef[sl])}
            exp = {kk: v for kk, v in exp.items() if v != 0}
            if set(exp) != set(got):
                mism += 1
            for kk, v in exp.items():
                if kk in got:
                    nmon += 1
                    worst = max(worst, abs(got[kk] - float(v)) / abs(float(v)))
    return len(corner) ** 2, nmon, worst, mism


def theta_T_from_file(path):
    """exact T-format {(i,j): {(n,p,q): Fraction}} of a kernel_Theta file (from the stored Garner digits); (i,j) index general_kernel.pairs_of(lam)"""
    arrs = {k: v for k, v in np.load(path).items()}
    nS = len(arrs['pairs'])
    fr = theta_coefficients_exact(arrs)
    ptr = arrs['entry_ptr']
    n, p, q = arrs['n'], arrs['p'], arrs['q']
    T = {}
    for e in range(nS * nS):
        sl = range(ptr[e], ptr[e + 1])
        T[divmod(e, nS)] = {(int(n[m]), int(p[m]), int(q[m])): fr[m] for m in sl}
    return T, arrs


def official_sum_rule(path, trials=2):
    """G1's exact top-state sum-rule checker (`general_kernel.top_sum_rule_check`, sympy rationals, random rational points) applied to the EXACT coefficients of a file;
    also asserts that the pair list equals general_kernel.pairs_of(lam).  Returns worst |lhs - rhs|/(|rhs|+1)."""
    import general_kernel as G
    T, arrs = theta_T_from_file(path)
    lam = tuple(int(x) for x in arrs['lam'])
    assert [tuple(int(v) for v in r) for r in arrs['pairs']] == G.pairs_of(lam), 'pair list differs from general_kernel.pairs_of'
    return G.top_sum_rule_check(lam, T, trials=trials)


def N_exact_ints(M, rho, cols, primes=None, chunk=64):
    """EXACT integer columns N[:, cols] of N = M^{-1} as a Python-int object array of shape (n, len(cols)) (original row order); the number of primes defaults
    to a safe value for the largest L here (16 primes cover |N| < 2^300 ~ 1e90)."""
    primes = primes or primes_for_blas(16)
    solver = ModSolver(M, rho)
    for p in primes:
        solver.prepare(p)
    cols = np.asarray(cols)
    out = np.zeros((M.shape[0], len(cols)), dtype=object)
    for s in range(0, len(cols), chunk):
        cc = cols[s:s + chunk]
        res = [solver.solve(cc, p) for p in primes]
        digits, xf = garner_signed([r.reshape(-1) for r in res], primes)
        assert not np.any(digits[-1] != 0) and not np.any(digits[-2] != 0), 'not enough primes'
        out[:, s:s + len(cc)] = digits_to_ints(digits, primes).reshape(len(cc), -1).T
    return out


def G_eps_coefficients(arrs, e, wh, whp, t, dps=60, exact=True):
    """Coefficients g_k (k = 0 .. 2(2 l1 - l2)) of the polynomial in eps
        G_ij(eps*wh, eps*wh'; th12) = (w w')^(2 l1 - l2) Theta_ij(1/w, 1/w', t),  w = eps*wh, w' = eps*wh',  t = 1/th12,
    = sum_m c_m t^n wh^(2l1-l2-p) wh'^(2l1-l2-q) eps^(2(2l1-l2)-p-q)  (entry e = i*nS + j of a loaded kernel_Theta file; exponents are >= 0 because p, q <= 2 l1 - l2).
    exact=True uses the exact coefficients (mpmath, `dps` digits); exact=False the float64 ones.  Returns a list of mpmath numbers (index = power of eps)."""
    import mpmath as mp
    mp.mp.dps = dps
    l1, l2 = int(arrs['lam'][0]), int(arrs['lam'][1])
    D = 2 * l1 - l2
    sl = slice(int(arrs['entry_ptr'][e]), int(arrs['entry_ptr'][e + 1]))
    WH, WHP, T = mp.mpf(wh), mp.mpf(whp), mp.mpf(t)
    if exact:
        num, den = int(str(arrs['exact_num'])), int(str(arrs['exact_den']))
        coefs = [mp.mpf(int(c)) * num / den for c in theta_exact_numerators(arrs, sl)]
    else:
        coefs = [mp.mpf(float(c)) for c in arrs['coef'][sl]]
    g = [mp.mpf(0)] * (2 * D + 1)
    for c, n_, p_, q_ in zip(coefs, arrs['n'][sl], arrs['p'][sl], arrs['q'][sl]):
        k = 2 * D - int(p_) - int(q_)
        g[k] += c * T ** int(n_) * WH ** (D - int(p_)) * WHP ** (D - int(q_))
    return g


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else ''
    args = sys.argv[2:]
    if cmd == 'M':
        cmd_M(int(args[0]), int(args[1]) if len(args) > 1 else None)
    elif cmd == 'N':
        cmd_N(int(args[0]), int(args[1]) if len(args) > 1 and args[1] != '-' else None, float_compare=('nofloat' not in args), save=('nosave' not in args))
    elif cmd == 'Theta':
        cmd_Theta(int(args[0]), int(args[1]), int(args[2]) if len(args) > 2 and args[2] != '-' else None, save=('nosave' not in args))
    elif cmd == 'univ':
        cmd_univ([int(a) for a in args])
    elif cmd == 'accept':
        cmd_accept(int(args[0]) if args else None)
    else:
        print(__doc__)
        print('commands: M <L> [batch] | N <L> [nprimes|- [nofloat] [nosave]] | Theta <L> <l2> [nprimes|- [nosave]] | univ <L...> | accept [Lbig]')
