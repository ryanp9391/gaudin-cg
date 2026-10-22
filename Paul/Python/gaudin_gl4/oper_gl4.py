"""oper_gl4.py -- WP-G5: the oper side of the gl(4) missing labels.

On a joint eigenvector of G2's commuting family the quantum characteristic polynomial cdet(d_u - L(u)) of the two-site gl(4) Gaudin model becomes a SCALAR
4th-order Fuchsian operator (the oper)
        D = d^4 - s_1 d^3 + s_2 d^2 - s_3 d + s_4,    s_k(u) = sum_{s=1,2} sum_{j=1}^{k} S^{(s)}_{k,j} / (u - theta_s)^j     (exact),
with three accessory parameters v = (v3, v4a, v4b), the eigenvalues of G2's triple V3 = t S_{3,2}, V4a = t S_{4,3}, V4b = t^2 S_{4,2}.

Conventions (binding, = gaudin_gl4_ops.py = notes section 1.1 generalised): [E_ij,E_kl] = d_jk E_il - d_li E_kj; two sites V_lam (site 1) (x) V_lam' (site 2);
L(u)_ij = e1_ij/(u-theta_1) + e2_ij/(u-theta_2); t = theta_12 = theta_1 - theta_2; S^{(2)} = site swap and t -> -t of S^{(1)}.  Here theta_2 = 0, theta_1 = t (any t; v is t-free).

Layout
  1. central characters  T_k(nu) = (-1)^k a_{4-k}(nu),  prod_c (th - y_c) = sum_k a_k th^{falling k},  y_c = nu_c + 4 - c;  tau_k = T_k(mu) - T_k(lam) - T_k(lam')
     THE DICTIONARY (proved from the residue relations S_k = O(u^-k), verified here on every sector): all 18 Laurent coefficients are
        S^{(1)}_{k,j} = t^(j-k) (a V3 + b V4a + c V4b + scalar),  S^{(2)}_{k,j} = (-t)^(j-k) (a' V3 + b' V4a + c' V4b + scalar')   with the table DICT_TABLE and the scalars
        T_k(lam), T_k(lam'), +-tau_k   (closed forms, see DICT_TABLE)
  2. Casimir (trace-form) eigenvalues C_2, C_3, C_4 via GLNRep, and the fits T_k(nu) = polynomial(n, C_2, C_3, C_4)
  3. sectors: exact HW spaces, restricted triples / Laurent coefficients (G2's machinery, read-only), joint eigenvalue triples (mpmath)
  4. the oper: affine-in-v Laurent data, expansion at infinity (indicial polynomial, Frobenius recursion), the pencil M_j(v) of D~ = (u-theta_1)^4 (u-theta_2)^4 D on polynomials
  5. G5.2: exponents, polynomial solutions of every copy, Wronskian, negative tests
  6. G5.3: the independent solve: exact (sympy rings / Groebner) and numerical certification
  7. G5.4: the July-2026 gaudinODE of GaudinSolver.wb
  8. stages (CLI: python -X utf8 oper_gl4.py dict | exps | solve | search | wb | big | all; `wbsym` = the 150 s symbolic derivation of the July-2026 Casimir dictionary, not in `all`; `all` = dict exps solve search wb big)
"""
import os, sys, time, itertools, math, json, re, random
from fractions import Fraction as Fr
from math import comb, factorial
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
import numpy as np
import mpmath as mp
import gaudin_gl4_ops as G2          # read-only: Rep, HWExact, build_ops, triple_ops, restrict_exact, joint_spectrum, racah_speiser, ...
import glN_reps as G1                # read-only: GLNRep (G1's representations)

REFDIR = G2.REFDIR
DICT_FILE = os.path.join(REFDIR, 'gl4-oper-dictionary.txt')
NN = 4


def lam_str(lam):
    return '[' + ','.join(str(x) for x in lam) + ']'


def hdr(s):
    print('\n' + '=' * 150 + '\n' + s + '\n' + '=' * 150, flush=True)


# ====================================================================================================== 1. central characters and the dictionary
def shifted(nu):
    """y_c = nu_c + 4 - c, c = 1..4"""
    return [nu[c] + NN - 1 - c for c in range(NN)]


def falling_coeffs(nu):
    """a_k(nu), k = 0..4:  prod_c (th - y_c) = sum_k a_k th^{falling k}   (a_k = Delta^k F(0) / k!)"""
    ys = shifted(nu)

    def F(th):
        r = Fr(1)
        for y in ys:
            r *= (th - y)
        return r
    d = [F(i) for i in range(NN + 1)]
    a = []
    for k in range(NN + 1):
        a.append(d[0] / factorial(k))
        d = [d[i + 1] - d[i] for i in range(len(d) - 1)]
    return a


def Tchar(nu, m):
    """T_m(nu) = (-1)^m a_{4-m}(nu): the central value of the single-site / total coefficient of u^-m of S_m (m = 1..4)"""
    return (-1) ** m * falling_coeffs(nu)[NN - m]


def tau_of(l1, l2, mu, k):
    return Tchar(mu, k) - Tchar(l1, k) - Tchar(l2, k)


# DICT_TABLE[(k, j, s)] = (integer coefficients (a, b, c) of (V3, V4a, V4b), scalar kind) with
#   S^{(1)}_{k,j} = t^(j-k) [a V3 + b V4a + c V4b + sc],   S^{(2)}_{k,j} = (-t)^(j-k) [a V3 + b V4a + c V4b + sc]   (the hat of site 2 = swapped hat)
# scalar kinds: ('T', site, k) = T_k(lam) / T_k(lam')  (the central top poles);  ('tau', k, sign) = sign * (T_k(mu) - T_k(lam) - T_k(lam'));  ('zero',)
DICT_TABLE = {
    (2, 2, 1): ((0, 0, 0), ('T', 1, 2)), (2, 2, 2): ((0, 0, 0), ('T', 2, 2)),
    (2, 1, 1): ((0, 0, 0), ('tau', 2, +1)), (2, 1, 2): ((0, 0, 0), ('tau', 2, +1)),
    (3, 3, 1): ((0, 0, 0), ('T', 1, 3)), (3, 3, 2): ((0, 0, 0), ('T', 2, 3)),
    (3, 2, 1): ((1, 0, 0), ('zero',)), (3, 2, 2): ((-1, 0, 0), ('tau', 3, +1)),
    (3, 1, 1): ((-2, 0, 0), ('tau', 3, +1)), (3, 1, 2): ((2, 0, 0), ('tau', 3, -1)),
    (4, 4, 1): ((0, 0, 0), ('T', 1, 4)), (4, 4, 2): ((0, 0, 0), ('T', 2, 4)),
    (4, 3, 1): ((0, 1, 0), ('zero',)), (4, 3, 2): ((0, -2, -1), ('tau', 4, +1)),
    (4, 2, 1): ((0, 0, 1), ('zero',)), (4, 2, 2): ((0, 3, 2), ('tau', 4, -1)),
    (4, 1, 1): ((0, -3, -3), ('tau', 4, +1)), (4, 1, 2): ((0, -3, -3), ('tau', 4, +1)),
}


def scalar_of(kind, l1, l2, mu):
    if kind[0] == 'zero':
        return Fr(0)
    if kind[0] == 'T':
        return Tchar(l1 if kind[1] == 1 else l2, kind[2])
    if kind[0] == 'tau':
        return kind[2] * tau_of(l1, l2, mu, kind[1])
    raise ValueError(kind)


def dictionary(l1, l2, mu):
    """{(k, j, s): ((a, b, c), scalar Fraction)} : the t-free hat of S^{(s)}_{k,j} = (a V3 + b V4a + c V4b + scalar)"""
    return {key: (coef, scalar_of(kind, l1, l2, mu)) for key, (coef, kind) in DICT_TABLE.items()}


# ---- affine forms in v = (x, y, z) = (v3, v4a, v4b):  tuple (cx, cy, cz, c0) of Fractions
AFF0 = (Fr(0), Fr(0), Fr(0), Fr(0))


def aff_const(c):
    return (Fr(0), Fr(0), Fr(0), Fr(c))


def aff_add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def aff_scale(a, c):
    c = Fr(c)
    return tuple(x * c for x in a)


def aff_eval(a, v):
    return a[0] * v[0] + a[1] * v[1] + a[2] * v[2] + a[3]


# ====================================================================================================== 2. Casimir (trace-form) eigenvalues and fits of T_k
def casimir_trace(nu, k):
    """eigenvalue of C_k = sum e_{i2 i1} e_{i3 i2} ... e_{i1 i_k} (G2's trace form, leftmost factor acts last) on V_nu, from the HW vector of G1's GLNRep"""
    R = G1.GLNRep.of(tuple(nu))
    hw = {R.top: Fr(1)}
    tot = {}
    for idx in itertools.product(range(1, NN + 1), repeat=k):
        # word = e_{i2 i1} e_{i3 i2} ... e_{i1 ik}  -> factor f: (i_{f+1}, i_f), last factor (i_1, i_k); rightmost acts first
        facs = [(idx[(f + 1) % k], idx[f]) for f in range(k)]
        vec = hw
        for ij in reversed(facs):
            vec = R.apply(ij, vec)
            if not vec:
                break
        for r, v in vec.items():
            tot[r] = tot.get(r, 0) + v
    tot = {r: v for r, v in tot.items() if v != 0}
    assert set(tot) <= {R.top}, ('trace-form Casimir is not central?', nu, k)
    return tot.get(R.top, Fr(0))


_CAS_CACHE = {}


def casimirs(nu):
    """(n, C2, C3, C4) of V_nu (exact)"""
    nu = tuple(nu)
    if nu not in _CAS_CACHE:
        _CAS_CACHE[nu] = (sum(nu),) + tuple(casimir_trace(nu, k) for k in (2, 3, 4))
    return _CAS_CACHE[nu]


def partitions4(maxpart):
    return [(a, b, c, d) for a in range(maxpart + 1) for b in range(a + 1) for c in range(b + 1) for d in range(c + 1)]


def cas_monomials(k):
    """monomials n^a C2^b C3^c C4^d of weight a + 2b + 3c + 4d <= k, as exponent tuples"""
    out = []
    for a in range(k + 1):
        for b in range(k // 2 + 1):
            for c in range(k // 3 + 1):
                for d in range(k // 4 + 1):
                    if a + 2 * b + 3 * c + 4 * d <= k:
                        out.append((a, b, c, d))
    return out


def solve_exact_ls(rows, rhs):
    """exact solution of the (possibly overdetermined, consistent) system rows x = rhs over Q by Gaussian elimination; returns (solution list or None if inconsistent, rank)"""
    n = len(rows[0])
    A = [list(r) + [b] for r, b in zip(rows, rhs)]
    piv = []
    r = 0
    for c in range(n):
        p = None
        for i in range(r, len(A)):
            if A[i][c] != 0:
                p = i
                break
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        inv = 1 / A[r][c]
        A[r] = [x * inv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c)
        r += 1
    for i in range(r, len(A)):
        if A[i][n] != 0:
            return None, r
    sol = [Fr(0)] * n
    for i, c in enumerate(piv):
        sol[c] = A[i][n]
    return sol, r


def fit_T_in_casimirs(ks=(2, 3, 4), maxpart=4, heldout=None):
    """exact fit T_k(nu) = sum c_mono n^a C2^b C3^c C4^d (weight <= k) on all partitions with nu_1 <= maxpart (fit set), verified on the held-out partitions.
       Returns {k: ({mono: coeff}, n_fit, n_heldout, rank, ok)}"""
    parts = partitions4(maxpart)
    out = {}
    for k in ks:
        monos = cas_monomials(k)

        def row(nu):
            n, c2, c3, c4 = casimirs(nu)
            return [Fr(n) ** a * c2 ** b * c3 ** c * c4 ** d for (a, b, c, d) in monos]
        fitset = [p for p in parts if sum(p) <= 12]
        sol, rank = solve_exact_ls([row(p) for p in fitset], [Tchar(p, k) for p in fitset])
        ok = sol is not None
        nh = 0
        if ok and heldout:
            for p in heldout:
                nh += 1
                v = sum((c * r for c, r in zip(sol, row(p))), Fr(0))
                ok = ok and (v == Tchar(p, k))
        out[k] = ({m: c for m, c in zip(monos, sol) if c != 0} if sol is not None else None, len(fitset), nh, rank, ok)
    return out


def fmt_cas_poly(coefs):
    names = ['n', 'C2', 'C3', 'C4']
    terms = []
    for mono in sorted(coefs, key=lambda m: (-(m[0] + 2 * m[1] + 3 * m[2] + 4 * m[3]), m)):
        c = coefs[mono]
        fac = '*'.join((names[i] + ('^%d' % e if e > 1 else '')) for i, e in enumerate(mono) if e > 0)
        s = str(c) if c.denominator == 1 else '(%s)' % c
        terms.append(s + ('*' + fac if fac else ''))
    return ' + '.join(terms).replace('+ -', '- ')


# ====================================================================================================== 3. sectors
class Product:
    """V_lam (x) V_lam' with G2's scaffold; exact operators built lazily at t = 1 (hats; theta_2 hats = swapped, (-1)^p scaling is inside G2.build_ops)"""
    _cache = {}

    @classmethod
    def of(cls, l1, l2):
        key = (tuple(l1), tuple(l2))
        if key not in cls._cache:
            cls._cache[key] = cls(*key)
        return cls._cache[key]

    def __init__(self, l1, l2):
        self.l1, self.l2 = tuple(l1), tuple(l2)
        self.R1, self.R2 = G2.Rep.of(self.l1), G2.Rep.of(self.l2)
        self.rs = G2.racah_speiser(self.l1, self.l2)
        assert G2.sum_dim_check(self.l1, self.l2, self.rs)
        self._ops = None
        self._trip = None
        self._sectors = {}

    def all_ops(self, t=Fr(1)):
        """the 18 Laurent coefficients S^{(s)}_{k,j} (Op2) at the numeric t (G2.build_ops), key (k, j, s)"""
        if self._ops is None or self._ops[0] != t:
            self._ops = (t, G2.build_ops(self.R1, self.R2, 4, t))
        return self._ops[1]

    def triple(self):
        if self._trip is None:
            self._trip = G2.triple_ops(self.R1, self.R2)
        return self._trip

    def sector(self, mu):
        mu = tuple(mu)
        if mu not in self._sectors:
            self._sectors[mu] = Sector(self, mu)
        return self._sectors[mu]

    def irreps(self):
        return sorted(self.rs.items(), reverse=True)


class Sector:
    """HW space of weight mu in V_lam (x) V_lam' (exact), the restricted triple (exact m x m rational matrices) and the joint eigenvalue triples"""

    def __init__(self, prod, mu):
        self.prod, self.mu = prod, tuple(mu)
        self.l1, self.l2 = prod.l1, prod.l2
        self.m = prod.rs[self.mu]
        self.hw = G2.HWExact(prod.R1, prod.R2, self.mu)
        assert self.hw.m == self.m, ('HW dim != multiplicity', mu, self.hw.m, self.m)
        self.R = {}
        tr = prod.triple()
        viol = 0
        for name in ('V3', 'V4a', 'V4b'):
            Rm, v = G2.restrict_exact(tr[name], self.hw)
            self.R[name] = Rm
            viol += v
        self.faithful_viol = viol
        self.comm_bad = sum(0 if G2.mat_fr_zero(G2.mat_fr_comm(self.R[a], self.R[b])) else 1 for a, b in itertools.combinations(('V3', 'V4a', 'V4b'), 2))
        self._js = {}

    def Rmats(self):
        return [self.R['V3'], self.R['V4a'], self.R['V4b']]

    def joint(self, dps=60):
        """joint eigenvalue triples (mp, dps digits), sorted by the generic combination; m = 1: exact rationals"""
        if dps not in self._js:
            if self.m == 1:
                trip = [[mp.mpf(self.R[n][0][0].numerator) / self.R[n][0][0].denominator for n in ('V3', 'V4a', 'V4b')]]
                self._js[dps] = dict(trip=trip, offmax=mp.mpf(0), exact=[tuple(self.R[n][0][0] for n in ('V3', 'V4a', 'V4b'))])
            else:
                js = G2.joint_spectrum(self.Rmats(), G2.GENERIC_COMBO, dps)
                js['exact'] = None
                self._js[dps] = js
        return self._js[dps]

    def exact_triples(self):
        """list of exact rational triples if all eigenvalues are rational (m = 1, or all roots rational), else None"""
        js = self.joint(60)
        return js['exact']

    def restrict(self, op):
        Rm, v = G2.restrict_exact(op, self.hw)
        assert v == 0
        return Rm


def mat_lin(R3, R4a, R4b, coef, sc):
    """a V3 + b V4a + c V4b + sc*1 as an exact matrix"""
    m = len(R3)
    M = [[coef[0] * R3[i][j] + coef[1] * R4a[i][j] + coef[2] * R4b[i][j] for j in range(m)] for i in range(m)]
    for i in range(m):
        M[i][i] += sc
    return M


def check_dictionary_sector(sec, t=Fr(1)):
    """G5.1 on one sector: restrict all 18 Laurent coefficients (G2.build_ops at the numeric t), divide by the power of t, subtract the integer combination of the restricted triple,
       verify the remainder is the scalar matrix of the closed form.  Returns (bad keys, {key: scalar actually found (exact) or None if the remainder is not scalar})"""
    prod = sec.prod
    ops = prod.all_ops(t)
    m = sec.m
    bad = []
    found = {}
    D = dictionary(sec.l1, sec.l2, sec.mu)
    for (k, j, s), (coef, sc) in D.items():
        Rm = sec.restrict(ops[(k, j, s)])
        p = j - k
        fac = (t ** p) if s == 1 else ((-t) ** p)
        hat = [[x / fac for x in row] for row in Rm]
        rem = [[hat[i][jj] - (coef[0] * sec.R['V3'][i][jj] + coef[1] * sec.R['V4a'][i][jj] + coef[2] * sec.R['V4b'][i][jj]) for jj in range(m)] for i in range(m)]
        off = any(rem[i][jj] != 0 for i in range(m) for jj in range(m) if i != jj)
        diag = {rem[i][i] for i in range(m)}
        if off or len(diag) != 1:
            found[(k, j, s)] = None
            bad.append((k, j, s))
        else:
            c = next(iter(diag))
            found[(k, j, s)] = c
            if c != sc:
                bad.append((k, j, s))
    return bad, found


# ====================================================================================================== 4. the oper: affine Laurent data, expansion at infinity, pencil
def falling(n, k):
    r = 1
    for i in range(k):
        r *= (n - i)
    return r


def poly_mul(a, b):
    out = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x != 0:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def poly_pow_lin(root, n):
    """(u - root)^n, ascending Fractions"""
    out = [Fr(1)]
    for _ in range(n):
        out = poly_mul(out, [Fr(-root), Fr(1)])
    return out


class Oper:
    """The oper D = d^4 - s1 d^3 + s2 d^2 - s3 d + s4 of the sector (l1, l2, mu), theta_2 = 0, theta_1 = t, with the accessory parameters v = (x, y, z) = (v3, v4a, v4b) kept as the
       unknowns: every Laurent coefficient, every expansion coefficient at infinity and every entry of the pencil is an AFFINE form in v (tuple (cx, cy, cz, c0) of Fractions)."""

    def __init__(self, l1, l2, mu, t=Fr(1)):
        self.l1, self.l2, self.mu = tuple(l1), tuple(l2), tuple(mu)
        self.t = Fr(t)
        self.theta = {1: self.t, 2: Fr(0)}
        self.n = {1: sum(self.l1), 2: sum(self.l2)}
        self.d = [self.mu[c] + NN - 1 - c for c in range(NN)]          # mu_c + 4 - c  (decreasing)
        self.y1 = shifted(self.l1)                                      # exponents at theta_1
        self.y2 = shifted(self.l2)                                      # exponents at theta_2
        D = dictionary(self.l1, self.l2, self.mu)
        S = {}
        for s in (1, 2):
            S[(s, 1, 1)] = aff_const(self.n[s])
        for (k, j, s), (coef, sc) in D.items():
            base = (Fr(coef[0]), Fr(coef[1]), Fr(coef[2]), Fr(sc))
            fac = (self.t if s == 1 else -self.t) ** (j - k)
            S[(s, k, j)] = aff_scale(base, fac)
        self.S = S
        self._sig = {}
        self._I = {}
        self._P = None

    # ---- expansion at infinity
    def sigma(self, k, r):
        """coefficient of u^{-r} in s_k(u): sum_s sum_{j<=min(r,k)} C(r-1, j-1) theta_s^{r-j} S^{(s)}_{k,j}   (affine form)"""
        key = (k, r)
        if key not in self._sig:
            tot = AFF0
            if r >= 1:
                for s in (1, 2):
                    th = self.theta[s]
                    for j in range(1, min(r, k) + 1):
                        c = comb(r - 1, j - 1) * (th ** (r - j))
                        if c != 0:
                            tot = aff_add(tot, aff_scale(self.S[(s, k, j)], c))
            self._sig[key] = tot
        return self._sig[key]

    def I(self, R, n):
        """I_R(n):  D u^n = sum_{R >= 0} I_R(n) u^{n-4-R}   (affine form; I_0 is the indicial polynomial at infinity)"""
        key = (R, n)
        if key not in self._I:
            tot = aff_const(falling(n, 4)) if R == 0 else AFF0
            for k in range(1, NN + 1):
                f = falling(n, NN - k)
                if f != 0:
                    tot = aff_add(tot, aff_scale(self.sigma(k, k + R), (-1) ** k * f))
            self._I[key] = tot
        return self._I[key]

    def I0_poly_coeffs(self):
        """the indicial polynomial at infinity as coefficients of the falling factorials n^{(4)}, ..., n^{(0)}: [1, -sigma_{1,1}, sigma_{2,2}, -sigma_{3,3}, sigma_{4,4}] (affine forms)"""
        return [aff_const(1)] + [aff_scale(self.sigma(k, k), (-1) ** k) for k in range(1, NN + 1)]

    def local_indicial(self, s):
        """indicial polynomial at theta_s from the most singular terms: [1, -S_{1,1}, S_{2,2}, -S_{3,3}, S_{4,4}] in falling factorials of the exponent e (affine forms)"""
        return [aff_const(1)] + [aff_scale(self.S[(s, k, k)], (-1) ** k) for k in range(1, NN + 1)]

    # ---- the polynomials P_k(u) = (u-theta_1)^4 (u-theta_2)^4 s_k(u)  (affine-form valued), and the pencil of D~ = (u-theta_1)^4 (u-theta_2)^4 D on monomials
    def Ppolys(self):
        """P_k, k = 1..4, as {degree: affine form}  (degree <= 8 - k)"""
        if self._P is None:
            t = self.t
            out = {}
            for k in range(1, NN + 1):
                P = defaultdict(lambda: AFF0)
                for j in range(1, k + 1):
                    # site 1 (u - t): (u-t)^{4-j} u^4 ;  site 2 (u - 0): u^{4-j} (u-t)^4
                    q1 = poly_mul(poly_pow_lin(t, 4 - j), poly_pow_lin(0, 4))
                    q2 = poly_mul(poly_pow_lin(0, 4 - j), poly_pow_lin(t, 4))
                    for deg, c in enumerate(q1):
                        if c != 0:
                            P[deg] = aff_add(P[deg], aff_scale(self.S[(1, k, j)], c))
                    for deg, c in enumerate(q2):
                        if c != 0:
                            P[deg] = aff_add(P[deg], aff_scale(self.S[(2, k, j)], c))
                out[k] = dict(P)
            self._P = out
        return self._P

    def pencil_column(self, i):
        """D~ u^i = sum_m M[m][i] u^m, m = 0..i+4, as a list of affine forms"""
        P = self.Ppolys()
        t = self.t
        col = [AFF0] * (i + 5)
        pi4 = poly_mul(poly_pow_lin(0, 4), poly_pow_lin(t, 4))             # (u(u-t))^4, degree 8
        f4 = falling(i, 4)
        if f4 != 0:
            for deg, c in enumerate(pi4):
                col[i - 4 + deg] = aff_add(col[i - 4 + deg], aff_const(f4 * c))
        for k in range(1, NN + 1):
            f = falling(i, NN - k)
            if f == 0:
                continue
            sh = i - (NN - k)
            sgn = (-1) ** k
            for deg, a in P[k].items():
                if a == AFF0:
                    continue                      # (the would-be top coefficients cancel by the residue relations S_k = O(u^-k))
                col[sh + deg] = aff_add(col[sh + deg], aff_scale(a, sgn * f))
        return col

    def pencil(self, dmax):
        """the (dmax+5) x (dmax+1) pencil as four exact Fraction matrices (Mx, My, Mz, M0): M(v) = M0 + x Mx + y My + z Mz"""
        cols = [self.pencil_column(i) for i in range(dmax + 1)]
        rows = dmax + 5
        mats = [[[Fr(0)] * (dmax + 1) for _ in range(rows)] for _ in range(4)]
        for i, col in enumerate(cols):
            for m, a in enumerate(col):
                for c in range(4):
                    mats[c][m][i] = a[c]
        return mats

    def matrix(self, v, dmax, mode='exact', mats=None):
        """M(v): exact (v Fractions) or mpmath (v mpf)"""
        Mx, My, Mz, M0 = mats if mats is not None else self.pencil(dmax)
        rows, cols = dmax + 5, dmax + 1
        if mode == 'exact':
            return [[M0[m][i] + v[0] * Mx[m][i] + v[1] * My[m][i] + v[2] * Mz[m][i] for i in range(cols)] for m in range(rows)]
        conv = lambda f: mp.mpf(f.numerator) / f.denominator
        out = mp.zeros(rows, cols)
        for m in range(rows):
            for i in range(cols):
                out[m, i] = conv(M0[m][i]) + v[0] * conv(Mx[m][i]) + v[1] * conv(My[m][i]) + v[2] * conv(Mz[m][i])
        return out

    def eval_laurent(self, v):
        """{(s,k,j): value} at the numeric v"""
        return {key: aff_eval(a, v) for key, a in self.S.items()}

    def residue_defects(self):
        """affine forms of the residue defects sigma_{k,r}, r = 1..k-1, k = 2..4 (must all be identically zero: S_k = O(u^-k))"""
        return {(k, r): self.sigma(k, r) for k in range(2, NN + 1) for r in range(1, k)}


# ---- exact linear algebra over Q on small dense matrices (lists of lists of Fractions)
def mat_rank_exact(M):
    A = [list(r) for r in M]
    rows, cols = len(A), len(A[0]) if A else 0
    r = 0
    for c in range(cols):
        p = None
        for i in range(r, rows):
            if A[i][c] != 0:
                p = i
                break
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        inv = 1 / A[r][c]
        A[r] = [x * inv for x in A[r]]
        for i in range(r + 1, rows):
            if A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        r += 1
        if r == rows:
            break
    return r


def nullspace_exact_dense(M):
    """exact nullspace basis (list of vectors) of the dense matrix M (rows x cols), RREF based"""
    A = [list(r) for r in M]
    rows, cols = len(A), len(A[0])
    piv = []
    r = 0
    for c in range(cols):
        p = None
        for i in range(r, rows):
            if A[i][c] != 0:
                p = i
                break
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        inv = 1 / A[r][c]
        A[r] = [x * inv for x in A[r]]
        for i in range(rows):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c)
        r += 1
        if r == rows:
            break
    free = [c for c in range(cols) if c not in piv]
    basis = []
    for f in free:
        vec = [Fr(0)] * cols
        vec[f] = Fr(1)
        for i, c in enumerate(piv):
            vec[c] = -A[i][f]
        basis.append(vec)
    return basis


# ====================================================================================================== 6a. exact elimination: Frobenius conditions at infinity, Groebner bases (sympy rings)
from sympy.polys.rings import ring as _sym_ring
from sympy.polys.domains import QQ as _QQ
from sympy.polys.orderings import grevlex as _grevlex, lex as _lex
from sympy.polys.groebnertools import groebner as _groebner


def qq(fr):
    fr = Fr(fr)
    return _QQ(fr.numerator, fr.denominator)


class RingCtx:
    """sympy polynomial ring Q[x, y, z (, f2, f3, f4)] with the accessory parameters v = (x, y, z) and (optionally) the free resonance parameters f_k"""

    def __init__(self, nfree=0, order=_grevlex):
        names = ['x', 'y', 'z'] + ['f%d' % k for k in (2, 3, 4)][:nfree]
        out = _sym_ring(','.join(names), _QQ, order)
        self.R = out[0]
        self.gens = out[1:]
        self.zero, self.one = self.R.zero, self.R.one
        self.v = self.gens[:3]
        self.f = {k: g for k, g in zip((2, 3, 4), self.gens[3:])}

    def from_aff(self, a):
        return qq(a[0]) * self.v[0] + qq(a[1]) * self.v[1] + qq(a[2]) * self.v[2] + qq(a[3]) * self.one

    def const(self, c):
        return qq(c) * self.one


def I0_value(op, i):
    """I_0(i) = prod_c (i - d_c)  (v-independent)"""
    r = Fr(1)
    for dc in op.d:
        r *= (i - dc)
    return r


def frob_conditions(op, j, ctx, free=()):
    """Frobenius recursion at infinity for a polynomial solution p = sum p_i u^i of degree d_j (p_{d_j} = 1):  c_m = sum_{q=0}^{m} p_{d-q} I_{m-q}(d-q) = 0 for m = 0..d+4
       (these d+5 equations are equivalent to D~ p = 0).  m <= d: if d - m is not an exponent, it determines p_{d-m}; if d - m = d_k (k > j) it is a RESONANCE: I_0 = 0, the equation is a
       CONDITION on v and p_{d_k} is free: gauge-fixed to 0 (valid when the solution of degree d_k exists) unless k is in `free`, then p_{d_k} = f_k (ring variable).  m > d: bottom conditions.
       Returns (conds, p): condition polynomials in ctx.R (polynomial in v, linear in the f_k) and the coefficient dict."""
    d = op.d[j - 1]
    exps = {op.d[k - 1]: k for k in range(j + 1, NN + 1)}
    zero, one = ctx.zero, ctx.one
    p = {d: one}
    conds = []
    Icache = {}

    def Ir(R, n):
        key = (R, n)
        if key not in Icache:
            Icache[key] = ctx.from_aff(op.I(R, n))
        return Icache[key]
    for m in range(1, d + 5):
        acc = zero
        for q in range(0, min(m, d) + 1):
            if q == m:
                continue
            i = d - q
            pi = p.get(i)
            if pi is None or pi == 0:
                continue
            acc = acc + pi * Ir(m - q, i)
        if m <= d:
            i = d - m
            if i in exps:
                if acc != 0:
                    conds.append(acc)
                k = exps[i]
                p[i] = ctx.f[k] if k in free else zero
            else:
                p[i] = acc * qq(Fr(-1) / I0_value(op, i))
        else:
            if acc != 0:
                conds.append(acc)
    return conds, p


def monomial_support(m):
    return frozenset(i for i, e in enumerate(m) if e > 0)


def ideal_dimension(lms, nvars=3):
    """Krull dimension from the leading monomials of a Groebner basis (max size of an independent set of variables: no leading monomial supported inside it)"""
    best = -1
    for r in range(nvars, -1, -1):
        for S in itertools.combinations(range(nvars), r):
            S = set(S)
            if not any(monomial_support(m) <= S for m in lms):
                best = max(best, r)
        if best >= 0:
            return best
    return -1


def standard_monomials_count(lms, nvars=3):
    """number of standard monomials (= dim_Q of the quotient algebra) of a zero-dimensional ideal, or None"""
    lms = [tuple(m) for m in lms]
    bounds = []
    for i in range(nvars):
        pures = [m[i] for m in lms if all((m[k] == 0) for k in range(nvars) if k != i) and m[i] > 0]
        if not pures:
            return None
        bounds.append(min(pures))
    cnt = 0
    for mon in itertools.product(*[range(b) for b in bounds]):
        if not any(all(mon[i] >= m[i] for i in range(nvars)) for m in lms):
            cnt += 1
    return cnt


def gb_info(gens, ctx):
    """Groebner basis (grevlex) of the ideal generated by gens; returns (basis list, dimension, #standard monomials or None)"""
    gens = [g for g in gens if g != 0]
    if not gens:
        return [], 3, None
    G = _groebner(gens, ctx.R)
    lms = [g.LM for g in G]
    dim = ideal_dimension(lms, 3)
    cnt = standard_monomials_count(lms, 3) if dim == 0 else None
    return G, dim, cnt


def eval_poly_on_matrices(poly, mats):
    """evaluate the ring polynomial (in x, y, z) at commuting exact m x m matrices (V3, V4a, V4b); returns the exact matrix"""
    m = len(mats[0])
    I = [[Fr(1) if i == j else Fr(0) for j in range(m)] for i in range(m)]
    zero = [[Fr(0)] * m for _ in range(m)]
    pow_cache = {}

    def mpow(a, e):
        key = (a, e)
        if key not in pow_cache:
            if e == 0:
                pow_cache[key] = I
            else:
                pow_cache[key] = G2.mat_fr_mul(mpow(a, e - 1), mats[a])
        return pow_cache[key]
    tot = zero
    for mon, c in poly.terms():
        c = Fr(int(c.numerator), int(c.denominator))
        P = I
        for a in range(3):
            if mon[a]:
                P = G2.mat_fr_mul(P, mpow(a, mon[a]))
        tot = [[tot[i][j] + c * P[i][j] for j in range(m)] for i in range(m)]
    return tot


def vanishes_on_matrices(poly, mats):
    M = eval_poly_on_matrices(poly, mats)
    return all(x == 0 for row in M for x in row)


# ====================================================================================================== 5. G5.2: exponents, polynomial solutions of every copy, Wronskian, negative tests
def mpf_of(fr):
    fr = Fr(fr)
    return mp.mpf(fr.numerator) / fr.denominator


def to_mp_matrix_list(M):
    return [[mpf_of(x) for x in row] for row in M]


def mp_singular_values(A):
    """singular values (descending list of mpf) of an mp matrix"""
    s = mp.svd_r(A, compute_uv=False)
    return sorted([s[i] for i in range(len(s))], reverse=True)


def abs_scale(op, v, dmax, mats):
    """entrywise scale of the pencil at v: Frobenius norm of |M0| + max(|x|,1)|Mx| + max(|y|,1)|My| + max(|z|,1)|Mz| (mp)"""
    Mx, My, Mz, M0 = mats
    w = [max(abs(x), mp.mpf(1)) for x in v]
    tot = mp.mpf(0)
    for m in range(dmax + 5):
        for i in range(dmax + 1):
            a = abs(mpf_of(M0[m][i])) + w[0] * abs(mpf_of(Mx[m][i])) + w[1] * abs(mpf_of(My[m][i])) + w[2] * abs(mpf_of(Mz[m][i]))
            tot += a * a
    return mp.sqrt(tot)


def poly_solutions_exact(op, v):
    """rational v: for j = 1..4 the exact null space dimension of the pencil on polynomials of degree <= d_j (list), and the gauge-fixed solutions P_j (coefficient vectors, ascending):
       p_{d_j} = 1, p_{d_k} = 0 for k > j.  Returns (nulldims, sols or None)"""
    nd = []
    sols = []
    for j in range(1, NN + 1):
        d = op.d[j - 1]
        M = op.matrix(v, d)
        ns = nullspace_exact_dense(M)
        nd.append(len(ns))
        # gauge fix: constraints p_{d_j} = 1 and p_{d_k} = 0 (k > j): solve for the combination of the null vectors
        if not ns:
            sols.append(None)
            continue
        idx = [d] + [op.d[k - 1] for k in range(j + 1, NN + 1)]
        # find c with sum_a c_a ns[a][idx[0]] = 1, sum_a c_a ns[a][idx[i]] = 0
        rows = [[ns[a][ii] for a in range(len(ns))] for ii in idx]
        rhs = [Fr(1)] + [Fr(0)] * (len(idx) - 1)
        c, rk = solve_exact_ls(rows, rhs)
        if c is None:
            sols.append(None)
            continue
        vec = [sum((c[a] * ns[a][i] for a in range(len(ns))), Fr(0)) for i in range(d + 1)]
        sols.append(vec)
    return nd, sols


def mp_lstsq(B, rhs):
    """least squares solution of B x = rhs (mp matrices, B tall, full column rank) via the SVD pseudo-inverse (mp.qr_solve divides by zero when a leading entry vanishes);
       returns (x, residual norm)"""
    U, S, V = mp.svd_r(B, full_matrices=False)
    k = len(S)
    smax = max(S[i] for i in range(k))
    tmp = U.T * rhs
    y = mp.matrix(k, 1)
    for i in range(k):
        y[i] = tmp[i] / S[i] if S[i] > smax * mp.mpf(10) ** (-(mp.mp.dps // 2)) else mp.mpf(0)
    x = V.T * y
    return x, mp.norm(B * x - rhs)


def poly_solutions_mp(op, v, dps=60):
    """irrational v (mp numbers): for each j the singular values of the pencil matrix on degree <= d_j; gauge-fixed solution via least squares.
       Returns list over j of dict(d, ncols, nnull_expected, svals (desc), gap, null_rel, kept_rel, sol (mp list), resid)"""
    out = []
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        mats_cache = {}
        for j in range(1, NN + 1):
            d = op.d[j - 1]
            if d not in mats_cache:
                mats_cache[d] = op.pencil(d)
            mats = mats_cache[d]
            A = op.matrix(v, d, mode='mp', mats=mats)
            sc = abs_scale(op, v, d, mats)
            sv = mp_singular_values(A)
            nnull = 5 - j
            ncols = d + 1
            r = ncols - nnull
            null_rel = sv[r] / sc if r < len(sv) else mp.mpf(0)
            gap = (sv[r - 1] / sv[r]) if (r > 0 and sv[r] != 0) else (mp.inf if r > 0 else None)
            kept_rel = (sv[r - 1] / sc) if r > 0 else None
            # gauge-fixed solution: unknown p_i, i not in {d_k, k > j}; p_{d} = 1
            fixed = {op.d[k - 1]: mp.mpf(0) for k in range(j + 1, NN + 1)}
            fixed[d] = mp.mpf(1)
            freecols = [i for i in range(ncols) if i not in fixed]
            rhs = mp.matrix([-sum((A[m, i] * fixed[i] for i in fixed), mp.mpf(0)) for m in range(A.rows)])
            sol = [None] * ncols
            resid = mp.mpf(0)
            if freecols:
                B = mp.matrix(A.rows, len(freecols))
                for m in range(A.rows):
                    for jj, i in enumerate(freecols):
                        B[m, jj] = A[m, i]
                x, res = mp_lstsq(B, rhs)
                resid = res
                for jj, i in enumerate(freecols):
                    sol[i] = x[jj]
            for i, val in fixed.items():
                sol[i] = val
            if not freecols:
                resid = mp.norm(rhs)
            out.append(dict(j=j, d=d, ncols=ncols, nnull=nnull, svals=sv, gap=gap, null_rel=null_rel, kept_rel=kept_rel, sol=sol, resid=resid, scale=sc))
    finally:
        mp.mp.dps = old
    return out


def p_trim(a, zero):
    a = list(a)
    while len(a) > 1 and a[-1] == zero:
        a.pop()
    return a


def p_mul(a, b, zero):
    out = [zero] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] = out[i + j] + x * y
    return out


def p_add(a, b, zero):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else zero) + (b[i] if i < len(b) else zero) for i in range(n)]


def p_der(a, zero):
    return [a[i] * i for i in range(1, len(a))] or [zero]


def p_neg(a):
    return [-x for x in a]


def wronskian_poly(sols, zero):
    """W = det[P_i^{(k)}]_{i = 1..4, k = 0..3} as a coefficient list (numbers of one kind)"""
    n = len(sols)
    ders = []
    for P in sols:
        row = [list(P)]
        for _ in range(n - 1):
            row.append(p_der(row[-1], zero))
        ders.append(row)
    tot = [zero]
    for perm in itertools.permutations(range(n)):
        sgn = G2.perm_sign(perm)
        prod = [zero + 1]
        for i in range(n):
            prod = p_mul(prod, ders[i][perm[i]], zero)
        tot = p_add(tot, prod if sgn > 0 else p_neg(prod), zero)
    return tot


def wronskian_target(op, one):
    """(u - theta_1)^{|lam|} (u - theta_2)^{|lam'|} with theta_1 = t, theta_2 = 0"""
    out = [one]
    for _ in range(op.n[1]):
        out = p_mul(out, [-op.t * one, one], one * 0)
    for _ in range(op.n[2]):
        out = p_mul(out, [one * 0, one], one * 0)
    return out


def taylor_shift(P, a, zero, one):
    """coefficients of P(a + x) in x (ascending)"""
    n = len(P)
    out = [zero] * n
    for i, c in enumerate(P):
        for k in range(i + 1):
            out[k] = out[k] + c * comb(i, k) * (a ** (i - k))
    return out


def vanishing_orders_exact(sols_vec, a):
    """exact: the set of vanishing orders at u = a of the span of the given polynomial vectors (Fractions): pivot columns of the row-reduced Taylor matrix"""
    rows = [taylor_shift(P, Fr(a), Fr(0), Fr(1)) for P in sols_vec]
    n = max(len(r) for r in rows)
    rows = [r + [Fr(0)] * (n - len(r)) for r in rows]
    piv = []
    A = [list(r) for r in rows]
    r = 0
    for c in range(n):
        p = None
        for i in range(r, len(A)):
            if A[i][c] != 0:
                p = i
                break
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        inv = 1 / A[r][c]
        A[r] = [x * inv for x in A[r]]
        for i in range(r + 1, len(A)):
            if A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c)
        r += 1
        if r == len(A):
            break
    return piv


def vanishing_orders_mp(sols_vec, a, tol=mp.mpf(10) ** (-25)):
    """mp: vanishing orders at u = a of the span: smallest r columns whose rank increases (SVD rank with a relative tolerance), rows normalised"""
    zero, one = mp.mpf(0), mp.mpf(1)
    rows = [taylor_shift(P, mpf_of(a) if not isinstance(a, mp.mpf) else a, zero, one) for P in sols_vec]
    n = max(len(r) for r in rows)
    rows = [r + [zero] * (n - len(r)) for r in rows]
    # normalise each row by its largest entry
    rows = [[x / max(abs(y) for y in r) for x in r] for r in rows]
    piv = []
    prev = 0
    for c in range(1, n + 1):
        A = mp.matrix([[rows[i][jj] for jj in range(c)] for i in range(len(rows))])
        sv = mp_singular_values(A)
        rk = sum(1 for s in sv if s > tol)
        if rk > prev:
            piv.append(c - 1)
            prev = rk
        if prev == len(rows):
            break
    return piv


def copy_check(op, v, exact=None, dps=60):
    """G5.2 (iii)-(iv) for one copy.  v = triple of mp numbers (exact: the same triple as Fractions or None).  Returns a result dict."""
    res = {}
    if exact is not None:
        nd, sols = poly_solutions_exact(op, exact)
        exp_nd = [5 - j for j in range(1, NN + 1)]
        res['mode'] = 'exact'
        res['nulldims'] = nd
        res['ok_null'] = (nd == exp_nd)
        # exact gauge-fixed solutions must have exact degree d_j
        res['ok_deg'] = all(s is not None and s[op.d[j]] == 1 and all(c == 0 for c in s[op.d[j] + 1:]) for j, s in enumerate(sols)) if all(s is not None for s in sols) else False
        res['ok_deg'] = res['ok_deg'] and all(len(sols[j]) == op.d[j] + 1 for j in range(NN))
        if res['ok_null'] and res['ok_deg']:
            W = p_trim(wronskian_poly(sols, Fr(0)), Fr(0))
            T = wronskian_target(op, Fr(1))
            lead = W[-1] / T[-1]
            res['wr_ok'] = (len(W) == len(T) and all(W[i] == lead * T[i] for i in range(len(W))))
            res['wr_lead'] = lead
            res['wr_dev'] = mp.mpf(0) if res['wr_ok'] else mp.mpf(1)
            res['orders1'] = vanishing_orders_exact(sols, op.t)
            res['orders2'] = vanishing_orders_exact(sols, 0)
        else:
            res['wr_ok'] = False
            res['orders1'] = res['orders2'] = None
        res['gap'] = None
        res['sols'] = sols
        return res
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        data = poly_solutions_mp(op, v, dps)
        res['mode'] = 'mp'
        gaps = [x['gap'] for x in data if x['gap'] is not None]
        res['gap'] = min(gaps) if gaps else None
        res['null_rel'] = max(x['null_rel'] for x in data)
        res['kept_rel'] = min((x['kept_rel'] for x in data if x['kept_rel'] is not None), default=None)
        res['resid'] = max(x['resid'] for x in data)
        res['ok_null'] = (res['gap'] is None or res['gap'] >= mp.mpf(10) ** 20) and res['null_rel'] < mp.mpf(10) ** (-(dps // 2))
        sols = [x['sol'] for x in data]
        res['sols'] = sols
        W = p_trim(wronskian_poly(sols, mp.mpf(0)), mp.mpf(0))
        T = wronskian_target(op, mp.mpf(1))
        lead = W[-1] / T[-1]
        wm = max(abs(c) for c in W)
        dev = max(abs(W[i] - lead * T[i]) for i in range(min(len(W), len(T)))) / wm if len(W) == len(T) else mp.mpf(1)
        res['wr_dev'] = dev
        res['wr_ok'] = (len(W) == len(T)) and dev < mp.mpf(10) ** (-(dps // 2))
        res['wr_lead'] = lead
        res['orders1'] = vanishing_orders_mp(sols, op.t)
        res['orders2'] = vanishing_orders_mp(sols, 0)
        return res
    finally:
        mp.mp.dps = old


def degree_d1_residual(op, v, mats, d):
    """least squares residual  min_c || M_1(v) (u^{d} + sum_{i<d} c_i u^i) ||  (a polynomial solution of degree EXACTLY d_1 exists iff this is 0), relative to ||M_1(v)||_F;
       also the smallest singular value of M_1(v) relative to the largest (the brief's measure: it also sees the lower-degree solutions that survive a perturbation)"""
    A = op.matrix(v, d, mode='mp', mats=mats)
    ncol = A.cols
    B = mp.matrix(A.rows, ncol - 1)
    for m in range(A.rows):
        for i in range(ncol - 1):
            B[m, i] = A[m, i]
    a = mp.matrix([A[m, ncol - 1] for m in range(A.rows)])
    if ncol - 1 == 0:
        res = mp.norm(a)
    else:
        x, res = mp_lstsq(B, -a)
    fro = mp.sqrt(sum(A[m, i] ** 2 for m in range(A.rows) for i in range(ncol)))
    sv = mp_singular_values(A)
    return res / fro, sv[-1] / sv[0], sum(1 for s_ in sv if s_ / sv[0] < mp.mpf(10) ** (-20))


def d1_solution_exists_exact(op, v):
    """EXACT, threshold-free: does D(v) (v rational) have a polynomial solution of degree EXACTLY d_1?  (column u^{d_1} of the pencil lies in the span of the other columns)"""
    d = op.d[0]
    M = op.matrix(v, d)
    A = [row[:-1] for row in M]
    return mat_rank_exact(M) == mat_rank_exact(A)


def negative_tests(op, trip, dps=60, seed=0, delta=mp.mpf(10) ** (-3)):
    """(v) random triple and the known triple perturbed by delta = 1e-3 in each component: (a) least-squares residual of a degree-EXACTLY-d_1 solution relative to ||M_1||_F,
       (b) smallest singular value of M_1 relative to the largest, (c) number of singular values below 1e-20 relative.  Entries: (resid_rel, smin_rel, nsmall)"""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        d = op.d[0]
        mats = op.pencil(d)
        rng = random.Random(seed)
        scale = max(max(abs(x) for x in trip), mp.mpf(1))
        rv = [mp.mpf(rng.uniform(-1, 1)) * 3 * scale for _ in range(3)]
        out = dict(random=degree_d1_residual(op, rv, mats, d), perturbed=[], true=degree_d1_residual(op, list(trip), mats, d))
        for i in range(3):
            v = list(trip)
            v[i] = v[i] + delta
            out['perturbed'].append(degree_d1_residual(op, v, mats, d))
        return out
    finally:
        mp.mp.dps = old

# ====================================================================================================== 6b. G5.3: the independent solve (exact quotient algebra, subsets of conditions, comparison with G2, Jacobian, random starts)
def monomial_poly(ctx, mon):
    p = ctx.one
    for a, e in zip(ctx.v, mon):
        if e:
            p = p * a ** e
    return p


def quotient_basis(G):
    """standard monomials (exponent tuples) of a zero-dimensional Groebner basis G"""
    lms = [tuple(g.LM) for g in G]
    bounds = []
    for i in range(3):
        pures = [m[i] for m in lms if all((m[k] == 0) for k in range(3) if k != i) and m[i] > 0]
        bounds.append(min(pures))
    return [mon for mon in itertools.product(*[range(b) for b in bounds]) if not any(all(mon[i] >= m[i] for i in range(3)) for m in lms)]


def poly_coeff_fr(c):
    return Fr(int(c.numerator), int(c.denominator))


def quotient_mult_matrices(G, ctx):
    """exact multiplication matrices M_x, M_y, M_z (columns = coordinates of x_i * b_k in the basis of standard monomials) of Q[v]/<G>, and the basis"""
    B = quotient_basis(G)
    idx = {b: i for i, b in enumerate(B)}
    n = len(B)
    mats = []
    for a in range(3):
        M = [[Fr(0)] * n for _ in range(n)]
        for k, b in enumerate(B):
            mon = list(b)
            mon[a] += 1
            nf = monomial_poly(ctx, tuple(mon)).rem(G)
            for m_, c in nf.terms():
                M[idx[tuple(m_)]][k] = poly_coeff_fr(c)
        mats.append(M)
    return mats, B


def mult_of_monomial(mats, mon):
    n = len(mats[0])
    P = [[Fr(1) if i == j else Fr(0) for j in range(n)] for i in range(n)]
    for a in range(3):
        for _ in range(mon[a]):
            P = G2.mat_fr_mul(P, mats[a])
    return P


def distinct_points(mats, B):
    """number of distinct points of the zero-dimensional ideal = rank over Q of the trace form Tr(M_{b_a} M_{b_b}) on the quotient algebra"""
    n = len(B)
    Ms = [mult_of_monomial(mats, b) for b in B]
    T = [[sum((G2.mat_fr_mul(Ms[a], Ms[b])[i][i] for i in range(n)), Fr(0)) for b in range(n)] for a in range(n)]
    return mat_rank_exact(T)


def poly_eval_mp(poly, v):
    tot = mp.mpf(0)
    for mon, c in poly.terms():
        t = mp.mpf(int(c.numerator)) / int(c.denominator)
        for a in range(len(mon)):
            if mon[a]:
                t *= v[a] ** mon[a]
        tot += t
    return tot


def jacobian_mp(polys, ctx, v):
    rows = []
    for f in polys:
        rows.append([poly_eval_mp(f.diff(ctx.v[a]), v) for a in range(3)])
    return mp.matrix(rows)


def jacobian_rank(polys, ctx, v, tol=mp.mpf(10) ** (-20)):
    J = jacobian_mp(polys, ctx, v)
    sv = mp_singular_values(J)
    return sum(1 for s in sv if s > tol * sv[0]), sv


SUBSETS_J = [tuple(S) for r in range(1, 5) for S in itertools.combinations((4, 3, 2, 1), r)]


# ---- Groebner bases of some subsets of the conditions (single j = 1, 2 with degree-7 generators) can take hours in sympy: they are run in a worker process with a time-box
import multiprocessing as _mpc


def _poly_to_dict(f):
    return {tuple(mon): (int(c.numerator), int(c.denominator)) for mon, c in f.terms()}


def _subset_task(polydicts, m):
    ctx = RingCtx(0)
    polys = [ctx.R.from_dict({mon: _QQ(a, b) for mon, (a, b) in d.items()}) for d in polydicts]
    G, dim, cnt = gb_info(polys, ctx)
    entry = dict(dim=dim, cnt=cnt)
    if dim == 0:
        mS, BS = quotient_mult_matrices(G, ctx)
        entry['distinct'] = distinct_points(mS, BS)
        entry['equalsK'] = (cnt == m)
    else:
        entry['distinct'] = None
        entry['equalsK'] = False
    return entry


def _worker_loop(conn):
    while True:
        try:
            msg = conn.recv()
        except EOFError:
            break
        if msg is None:
            break
        polydicts, m = msg
        try:
            out = _subset_task(polydicts, m)
        except Exception as e:
            out = dict(error=repr(e))
        conn.send(out)


class Watchdog:
    """one persistent worker process (spawn) computing Groebner invariants of a list of polynomials; run() returns dict(timeout=True) and restarts the worker when the time-box is exceeded"""

    def __init__(self):
        self.proc = None
        self.conn = None

    def _start(self):
        c = _mpc.get_context('spawn')
        parent, child = c.Pipe()
        self.proc = c.Process(target=_worker_loop, args=(child,), daemon=True)
        self.proc.start()
        self.conn = parent

    def run(self, polys, m, timeout):
        if self.proc is None or not self.proc.is_alive():
            self._start()
        self.conn.send(([_poly_to_dict(f) for f in polys], m))
        if self.conn.poll(timeout):
            return self.conn.recv()
        self.proc.terminate()
        self.proc.join()
        self.proc = None
        return dict(timeout=True)

    def close(self):
        if self.proc is not None and self.proc.is_alive():
            try:
                self.conn.send(None)
            except Exception:
                pass
            self.proc.join(timeout=2)
            if self.proc.is_alive():
                self.proc.terminate()
        self.proc = None


_WATCHDOG = [None]
SUBSET_TIMEOUT = 5.0


def watchdog():
    if _WATCHDOG[0] is None:
        _WATCHDOG[0] = Watchdog()
    return _WATCHDOG[0]


def solve_sector(l1, l2, mu, dps=60, subsets=True, want_points=True, verbose=False):
    """G5.3 on one sector: gauge-fixed Frobenius conditions for j = 4..1 (polynomial in v), Groebner basis of the full system and of every subset of the four conditions,
       the quotient algebra's points vs G2's joint triples.  Returns a result dict."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        prod = Product.of(l1, l2)
        sec = prod.sector(mu)
        m = sec.m
        op = Oper(l1, l2, mu, Fr(1))
        ctx = RingCtx(0)
        t0 = time.perf_counter()
        conds = {}
        for j in (4, 3, 2, 1):
            conds[j], _ = frob_conditions(op, j, ctx)
        res = dict(mu=mu, m=m, d=op.d, ncond={j: len(c) for j, c in conds.items()}, maxdeg=max((c.degree() for cs in conds.values() for c in cs), default=0), nW=sec.hw.nW)
        Rm = sec.Rmats()
        # every generator vanishes on the commuting G2 matrices  => V(conditions) contains the m joint triples, ideal J subset of K
        res['all_vanish_on_G2'] = all(vanishes_on_matrices(c, Rm) for cs in conds.values() for c in cs)
        res['t_build'] = time.perf_counter() - t0
        t0 = time.perf_counter()
        full = [c for j in (4, 3, 2, 1) for c in conds[j]]
        G, dim, cnt = gb_info(full, ctx)
        res['full'] = dict(dim=dim, cnt=cnt)
        res['t_gb'] = time.perf_counter() - t0
        if dim == 0:
            mats, B = quotient_mult_matrices(G, ctx)
            res['distinct'] = distinct_points(mats, B)
            res['ideal_equals_K'] = (cnt == m and res['all_vanish_on_G2'])
            if want_points:
                jsq = G2.joint_spectrum(mats, G2.GENERIC_COMBO, dps)
                pts = [tuple(row) for row in jsq['trip']]
                g2 = [tuple(row) for row in sec.joint(dps)['trip']]
                dev = mp.mpf(0)
                used = set()
                for pg in g2:
                    best, bi = None, None
                    for i, pp in enumerate(pts):
                        if i in used:
                            continue
                        dd = max(abs(pg[a] - pp[a]) for a in range(3))
                        if best is None or dd < best:
                            best, bi = dd, i
                    if bi is not None:
                        used.add(bi)
                        dev = max(dev, best)
                res['points'] = pts
                res['dev'] = dev
                res['offmax'] = jsq['offmax']
                # Jacobian rank of the whole condition system at each G2 triple
                ranks = []
                for pg in g2:
                    rk, sv = jacobian_rank(full, ctx, list(pg))
                    ranks.append(rk)
                res['jac_ranks'] = ranks
        if subsets:
            sub = {}
            for S in SUBSETS_J:
                gens = [c for j in S for c in conds[j]]
                t1 = time.perf_counter()
                entry = watchdog().run(gens, m, SUBSET_TIMEOUT)
                entry['dt'] = time.perf_counter() - t1
                if entry.get('timeout'):
                    entry.update(dim=None, cnt=None, distinct=None, equalsK=False)
                assert 'error' not in entry, entry
                sub[S] = entry
            res['subsets'] = sub
        return res
    finally:
        mp.mp.dps = old


# ====================================================================================================== 7. G5.4: the July-2026 gl(4) Gaudin ODE of GaudinSolver.wb
import sympy as sp

GS_PATH = os.path.normpath(os.path.join(HERE, '..', '..', 'Mathematica', 'old-conventions', 'old', 'Gaudin', 'Clean', 'GaudinSolver.wb'))


def load_gs_cells():
    with open(GS_PATH, encoding='utf8') as f:
        return [c['value'] for c in json.load(f)['cells']]


def find_cell(cells, needle):
    for c in cells:
        if needle in c:
            return c
    raise KeyError(needle)


WB_SYMS = {n: sp.Symbol(n) for n in ('C2', 'C4', 'S', 'P', 'w1', 'w2', 'z')}


def transcribe_gaudinODE(cells=None):
    """read the cell defining gaudinODE from the notebook JSON and translate the Mathematica expression to sympy:
       Derivative[n][Q][z] -> Q_n, Q[z] -> Q_0, ^ -> **.  Returns (coeffs {n: c_n(z)} with ODE = sum_n c_n Q^{(n)}, source text)"""
    cells = cells or load_gs_cells()
    src = find_cell(cells, 'gaudinODE = ')
    body = src[src.index('gaudinODE = ') + len('gaudinODE = '):]
    body = body[:body.index(';')]
    txt = re.sub(r'Derivative\[(\d)\]\[Q\]\[z\]', r'Q_\1', body)
    txt = ' '.join(txt.replace('Q[z]', 'Q_0').replace('^', '**').split())          # one line (the Mathematica source is wrapped)
    loc = dict(WB_SYMS)
    for n in range(5):
        loc['Q_%d' % n] = sp.Symbol('Q_%d' % n)
    expr = sp.expand(sp.parse_expr(txt, local_dict=loc, transformations=sp.parsing.sympy_parser.standard_transformations + (sp.parsing.sympy_parser.implicit_multiplication_application,)))
    coeffs = {n: sp.expand(expr.coeff(loc['Q_%d' % n])) for n in range(5)}
    assert sp.expand(sum(coeffs[n] * loc['Q_%d' % n] for n in range(5)) - expr) == 0
    return coeffs, body


def wb_muhat(mu):
    """muhat[mu_] := Table[mu[[i]] + 4 - i, {i, 4}]"""
    return [mu[i] + 4 - (i + 1) for i in range(4)]


def wb_xw(mu):
    """xw[mu_] := Table[mu[[i]] + 5/2 - i, {i, 4}]"""
    return [sp.Rational(mu[i]) + sp.Rational(5, 2) - (i + 1) for i in range(4)]


def wb_Pw(mu, k):
    return sum(x ** k for x in wb_xw(mu))


def wb_Cas2(mu):
    """Cas2[mu_] := Pw[mu, 2] - 5"""
    return wb_Pw(mu, 2) - 5


def wb_Cas4(mu):
    """Cas4[mu_] := Pw[mu, 4] + 4 Pw[mu, 3] - Pw[mu, 1] Pw[mu, 2] + Pw[mu, 2]/2 - 2 Pw[mu, 1]^2 - 10 Pw[mu, 1] - 51/4"""
    P1, P2, P3, P4 = (wb_Pw(mu, k) for k in (1, 2, 3, 4))
    return P4 + 4 * P3 - P1 * P2 + P2 / 2 - 2 * P1 ** 2 - 10 * P1 - sp.Rational(51, 4)


def wb_degQ(S_, P_, mu):
    """degQ[SS_, PP_, mu_] := Table[(SS + PP + 3) - muhat[mu][[k]], {k, 4}]"""
    return [(S_ + P_ + 3) - h for h in wb_muhat(mu)]


def wb_irreps(S_, P_):
    """irreps of [S,S,0,0] x [P,P,0,0] (all multiplicity 1), via Racah-Speiser of G1 (the notebook's own `irreps` peels Schur polynomials)"""
    rs = G1.racah_speiser((S_, S_, 0, 0), (P_, P_, 0, 0))
    assert all(m == 1 for m in rs.values())
    return sorted(rs, reverse=True)


def wb_ode_operator_matrix(coeffs, vals, nmax):
    """matrix (nmax+1) x (nmax+1) of the ODE on monomials z^i (exact Fractions) with the numeric substitution vals = {C2:.., C4:.., S:.., P:.., w1:.., w2:..}:
       column i = coefficients of ODE[z^i] (degree <= i)"""
    z = WB_SYMS['z']
    cpoly = {}
    for n in range(5):
        pol = sp.Poly(coeffs[n].subs(vals), z)
        cpoly[n] = {mon[0]: Fr(int(c.p), int(c.q)) for mon, c in pol.terms()}
    M = [[Fr(0)] * (nmax + 1) for _ in range(nmax + 1)]
    for i in range(nmax + 1):
        for n in range(5):
            f = falling(i, n)
            if f == 0:
                continue
            for deg, c in cpoly[n].items():
                row = deg + i - n
                assert 0 <= row <= nmax or c == 0 or True
                if row <= nmax:
                    M[row][i] += c * f
                else:
                    raise AssertionError(('degree overflow', i, n, deg))
    return M


def wb_poly_solution_degrees(coeffs, vals, nmax):
    """degrees n <= nmax at which the dimension of the polynomial solutions of degree <= n jumps (exact), and the total dimension"""
    M = wb_ode_operator_matrix(coeffs, vals, nmax)
    dims = []
    for n in range(nmax + 1):
        sub = [row[:n + 1] for row in M[:n + 1]] if False else [row[:n + 1] for row in M]
        dims.append(n + 1 - mat_rank_exact(sub))
    jumps = [n for n in range(nmax + 1) if dims[n] > (dims[n - 1] if n > 0 else 0)]
    return jumps, dims[-1]


def adjoint_gauge_coeffs(s_exprs, p1, p2, w1, w2):
    """G[Q] := Pi^p D^dagger [Pi^-p Q],  Pi^p = (z-w1)^p1 (z-w2)^p2, for D = d^4 - s1 d^3 + s2 d^2 - s3 d + s4  (s_exprs = [s1, s2, s3, s4] sympy expressions in z; p1, p2 may be symbols).
       D^dagger = d^4 + d^3 s1 + d^2 s2 + d s3 + s4 (formal adjoint, d acting on everything to its right); conjugation by Pi^p replaces d by L = d - ell, ell = p1/(z-w1) + p2/(z-w2).
       Returns {n: coefficient of Q^{(n)}} (rational functions of z)"""
    z = WB_SYMS['z']
    s1, s2, s3, s4 = s_exprs
    Q = sp.Function('Q')
    ell = p1 / (z - w1) + p2 / (z - w2)

    def L(h):
        return sp.diff(h, z) - ell * h
    q = Q(z)
    expr = L(L(L(L(q)))) + L(L(L(s1 * q))) + L(L(s2 * q)) + L(s3 * q) + s4 * q
    expr = sp.expand(expr)
    out = {}
    for n in range(5):
        d = sp.Derivative(Q(z), (z, n)) if n > 0 else Q(z)
        out[n] = sp.cancel(sp.together(expr.coeff(d)))
    return out


def our_oper_exprs(S_, P_, w1, w2, v, tm, swap=False):
    """our oper D (the dual oper of the two-site model with [S,S,0,0] at theta_1 = w1 and [P,P,0,0] at theta_2 = w2, or swapped) as s1..s4 sympy expressions in z,
       through the DICTIONARY with symbolic/numeric v = (v3, v4a, v4b) and T_k(mu) = tm[k] (k = 2,3,4); returns [s1, s2, s3, s4]"""
    z = WB_SYMS['z']
    lam1, lam2 = ((S_, S_, 0, 0), (P_, P_, 0, 0)) if not swap else ((P_, P_, 0, 0), (S_, S_, 0, 0))
    t = w1 - w2

    def Tsym(nu, k):
        ys = [nu[c] + NN - 1 - c for c in range(NN)]
        th = sp.Symbol('th')
        F = sp.prod([th - y for y in ys])
        vals = [sp.expand(F.subs(th, i)) for i in range(NN + 1)]
        a = []
        d = vals[:]
        for kk in range(NN + 1):
            a.append(sp.expand(d[0] / sp.factorial(kk)))
            d = [d[i + 1] - d[i] for i in range(len(d) - 1)]
        return (-1) ** k * a[NN - k]
    T1 = {k: Tsym(lam1, k) for k in (2, 3, 4)}
    T2 = {k: Tsym(lam2, k) for k in (2, 3, 4)}
    n1, n2 = sum(lam1), sum(lam2)
    s = {1: n1 / (z - w1) + n2 / (z - w2)}
    for k in (2, 3, 4):
        s[k] = 0
    for (k, j, site), (coef, kind) in DICT_TABLE.items():
        if kind[0] == 'zero':
            sc = 0
        elif kind[0] == 'T':
            sc = T1[kind[2]] if kind[1] == 1 else T2[kind[2]]
        else:
            sc = kind[2] * (tm[kind[1]] - T1[kind[1]] - T2[kind[1]])
        val = coef[0] * v[0] + coef[1] * v[1] + coef[2] * v[2] + sc
        fac = (t if site == 1 else -t) ** (j - k)
        pos = w1 if site == 1 else w2
        s[k] = s[k] + fac * val / (z - pos) ** j
    return [s[1], s[2], s[3], s[4]]


def wb_lhs(S_, P_, p1, p2, w1, w2, swap=False):
    """LHS of the matching: for n = 0..4 the pair (numerator, denominator) of 8 Pi^2 G_n (G = Pi^p D^dagger Pi^-p of our oper with unknowns v3, v4a, v4b, tm2, tm3, tm4 = T_k(mu)),
       the numerator being a polynomial in z linear in the unknowns.  S_, P_, w1, w2 may be numbers or symbols."""
    z = WB_SYMS['z']
    unk = sp.symbols('v3 v4a v4b tm2 tm3 tm4')
    tm = {2: unk[3], 3: unk[4], 4: unk[5]}
    s_exprs = our_oper_exprs(S_, P_, w1, w2, unk[:3], tm, swap)
    G = adjoint_gauge_coeffs(s_exprs, p1, p2, w1, w2)
    z_ = z
    Pi2 = (z_ - w1) ** 2 * (z_ - w2) ** 2
    out = {}
    for n in range(5):
        num, den = sp.fraction(sp.together(8 * Pi2 * G[n]))
        out[n] = (sp.expand(num), sp.expand(den))
    return unk, out


def wb_match_from_lhs(unk, lhs, coeffs, S_, P_, w1, w2, Cv):
    """linear system num_n - rhs_n * den_n = 0 (coefficients of z^i) for the numeric (C2, C4) = Cv; returns the linsolve solution set (empty if inconsistent)"""
    z = WB_SYMS['z']
    subs = {WB_SYMS['S']: S_, WB_SYMS['P']: P_, WB_SYMS['w1']: w1, WB_SYMS['w2']: w2, WB_SYMS['C2']: Cv[0], WB_SYMS['C4']: Cv[1]}
    eqs = []
    for n in range(5):
        num, den = lhs[n]
        rhs = coeffs[n].subs(subs)
        poly = sp.Poly(sp.expand(num - rhs * den), z)
        eqs.extend(poly.coeffs())
    return sp.linsolve(eqs, unk)


def wb_match(S_, P_, p1, p2, w1, w2, Cv, swap=False, coeffs=None):
    """convenience: LHS + match for one (C2, C4)"""
    unk, lhs = wb_lhs(S_, P_, p1, p2, w1, w2, swap)
    coeffs = coeffs or transcribe_gaudinODE()[0]
    return unk, wb_match_from_lhs(unk, lhs, coeffs, S_, P_, w1, w2, Cv)


# ---- the closed-form dictionary (derived symbolically with S, P, w1, w2, C2, C4 all symbolic: stage wbsym; verified below on every irrep of every tested (S, P)), N = S + P:
#       v3 = (S+1) (2N(N+2) - C2)/2,   v4a = 0,   v4b = T_4(mu) = (C2^2 + 2(N+1) C2 - 2 C4 + 4 N (N+2)(N+3))/8,
#       T_2(mu) = (4N^2 + 6N - C2)/2,  T_3(mu) = (N+2)(2N(N+2) - C2)/2 = (N+2)/(S+1) v3
def wb_dictionary(S_, P_, C2, C4):
    N = S_ + P_
    v3 = sp.Rational(S_ + 1, 2) * (2 * N * (N + 2) - C2)
    v4b = (C2 ** 2 + 2 * (N + 1) * C2 - 2 * C4 + 4 * N * (N + 2) * (N + 3)) / 8
    T2 = (4 * N * N + 6 * N - C2) / 2
    T3 = sp.Rational(N + 2, 2) * (2 * N * (N + 2) - C2)
    return dict(v3=v3, v4a=sp.Integer(0), v4b=v4b, T2=T2, T3=T3, T4=v4b)


def sympy_poly_to_fr(poly_expr, z):
    pol = sp.Poly(poly_expr, z)
    return {mon[0]: Fr(int(c.p), int(c.q)) for mon, c in pol.terms()}


def fpoly_der(a, k):
    for _ in range(k):
        a = [a[i] * i for i in range(1, len(a))] or [Fr(0)]
    return a


def apply_ode_to_poly(coeffs, subs, Q):
    """exact: sum_n c_n(z) Q^{(n)}(z) for the numeric substitution subs and the polynomial Q (list of Fractions); returns the coefficient list (trimmed)"""
    z = WB_SYMS['z']
    tot = [Fr(0)]
    for n in range(5):
        cp = sympy_poly_to_fr(coeffs[n].subs(subs), z)
        carr = [cp.get(i, Fr(0)) for i in range(max(cp) + 1)] if cp else [Fr(0)]
        tot = p_add(tot, p_mul(carr, fpoly_der(list(Q), n), Fr(0)), Fr(0))
    return p_trim(tot, Fr(0))


def poly_divide_linear_power(P, root, k):
    """exact division of the polynomial P by (z - root)^k; returns (quotient, remainder-is-zero flag)"""
    Q = list(P)
    for _ in range(k):
        # synthetic division by (z - root)
        n = len(Q) - 1
        if n < 1:
            return Q, False
        b = [Fr(0)] * n
        carry = Fr(0)
        for i in range(n, 0, -1):
            carry = Q[i] + carry * root
            b[i - 1] = carry
        rem = Q[0] + carry * root
        if rem != 0:
            return b, False
        Q = b
    return Q, True


# ====================================================================================================== 8. stages
REQUIRED_PRODUCTS = [((2, 1, 0, 0), (2, 1, 0, 0)), ((2, 1, 1, 0), (2, 1, 1, 0)), ((3, 2, 1, 0), (3, 2, 1, 0)), ((2, 2, 1, 0), (2, 1, 1, 0)),
                     ((1, 1, 0, 0), (2, 1, 0, 0)), ((4, 2, 1, 0), (4, 2, 1, 0))]
EXTRA_DICT_PRODUCTS = [((3, 2, 0, 0), (2, 1, 1, 0)), ((2, 2, 1, 0), (2, 2, 1, 0)), ((3, 2, 0, 0), (3, 2, 0, 0))]
EXTRA_COPY_PRODUCTS = [((3, 2, 1, 0), (4, 2, 1, 0)), ((4, 2, 2, 0), (4, 2, 2, 0))]
DICT_PRODUCTS = REQUIRED_PRODUCTS + EXTRA_DICT_PRODUCTS
EXPS_PRODUCTS = REQUIRED_PRODUCTS + EXTRA_DICT_PRODUCTS + EXTRA_COPY_PRODUCTS
SOLVE_PRODUCTS = REQUIRED_PRODUCTS + EXTRA_DICT_PRODUCTS + EXTRA_COPY_PRODUCTS + [((4, 3, 1, 0), (4, 3, 1, 0))]
WB_SP = [(1, 1), (2, 1), (2, 2), (3, 2), (3, 3)]
DPS = 60


def pname(l1, l2):
    return lam_str(l1) + 'x' + lam_str(l2)


def fr_str(x):
    x = Fr(x)
    return str(x.numerator) if x.denominator == 1 else '%d/%d' % (x.numerator, x.denominator)


# ---------------------------------------------------------------------------------------------------- stage dict
def casimir_forms():
    held = [p for p in partitions4(6) if p[0] >= 5 and G1.weyl_dim(p) <= 3000]
    random.Random(1).shuffle(held)
    held = held[:20]
    return fit_T_in_casimirs(maxpart=4, heldout=held), held


def tau_cas(forms, l1, l2, mu, k):
    """tau_k from the Casimir closed forms (n, C2, C3, C4 of lam, lam', mu)"""
    def T_(nu):
        n, c2, c3, c4 = casimirs(nu)
        return sum((c * Fr(n) ** a * c2 ** b * c3 ** cc * c4 ** d for (a, b, cc, d), c in forms[k][0].items()), Fr(0))
    return T_(mu) - T_(l1) - T_(l2)


def stage_dict(write=True):
    hdr('G5.1  THE DICTIONARY: all 18 Laurent coefficients S^{(s)}_{k,j} = t^{j-k} (a V3 + b V4a + c V4b + scalar), the scalars in closed form')
    print('DERIVATION (exact, from G2 facts only).  Unknowns A_{k,j} = S^{(1)}_{k,j}, B_{k,j} = S^{(2)}_{k,j} (k = 2,3,4; j = 1..k); known: the top poles A_{k,k} = T_k(lam), B_{k,k} = T_k(lam\') (central scalars) and the')
    print('  expansion at infinity: the coefficient of u^{-r} in S_k(u) is  sum_s sum_{j<=r} C(r-1,j-1) theta_s^{r-j} S^{(s)}_{k,j}  = 0 for r < k (residue relations, G2) and = T_k(mu) for r = k (the total Casimir:')
    print('  at infinity L -> E^tot/u, so cdet(d - L) -> cdet(d - E^tot/u) = u^-N prod_c (theta - y_c(mu)), theta = u d_u).  Solving these k relations (theta_1 - theta_2 = t) for the Laurent coefficients leaves exactly')
    print('  the free parameters A_{3,2} = v3/t, A_{4,3} = v4a/t, A_{4,2} = v4b/t^2 (the definition of the triple) and gives:')
    print('     t   S^{(1)}_{2,1} = tau_2,   S^{(2)}_{2,1} = -tau_2 / t')
    print('     t^2 S^{(1)}_{3,1} = tau_3 - 2 v3,   t S^{(1)}_{3,2} = v3,   S^{(2)}_{3,1} = -S^{(1)}_{3,1},   t S^{(2)}_{3,2} = v3 - tau_3')
    print('     t^3 S^{(1)}_{4,1} = tau_4 - 3 v4a - 3 v4b,  t S^{(1)}_{4,3} = v4a,  t^2 S^{(1)}_{4,2} = v4b,  S^{(2)}_{4,1} = -S^{(1)}_{4,1},')
    print('     t^2 S^{(2)}_{4,2} = 3 v4a + 2 v4b - tau_4,   t S^{(2)}_{4,3} = 2 v4a + v4b - tau_4')
    print('  with T_k(nu) = (-1)^k a_{4-k}(nu) [prod_c (th - y_c) = sum_k a_k th^{falling k}, y_c = nu_c + 4 - c] and tau_k = T_k(mu) - T_k(lam) - T_k(lam\').')
    print('  This reproduces G2\'s integer table (section 6c) and fixes ALL scalars; below it is verified against the restricted operators on every sector.\n')
    print('%-6s %-6s %-6s | %-45s | scalar' % ('k', 'j', 'site', 'hat = a V3 + b V4a + c V4b + scalar'))
    for (k, j, s), (coef, kind) in sorted(DICT_TABLE.items()):
        sc = {'T': lambda: 'T_%d(lam%s)' % (kind[2], '' if kind[1] == 1 else "'"), 'tau': lambda: '%stau_%d' % ('' if kind[2] > 0 else '-', kind[1]), 'zero': lambda: '0'}[kind[0]]()
        print('%-6d %-6d %-6d | (%2d, %2d, %2d)  t-power %2d%s                  | %s' % (k, j, s, coef[0], coef[1], coef[2], j - k, '' if s == 1 else ' (with (-t) at theta_2)', sc))
    t0 = time.perf_counter()
    forms, held = casimir_forms()
    print('\nCLOSED FORMS OF THE CENTRAL CHARACTERS IN THE TRACE-FORM CASIMIRS (n = |nu|, C_k = sum e_{i2 i1} e_{i3 i2} ... e_{i1 i_k} on the HW vector of V_nu, GLNRep): exact fit on %d partitions (nu_1 <= 4, |nu| <= 12),' % forms[2][1])
    print('  verified on %d held-out partitions with nu_1 in {5, 6} (rank of the fit systems: %s):' % (forms[2][2], {k: forms[k][3] for k in forms}))
    for k in (2, 3, 4):
        print('   T_%d(nu) = %s      [fit ok on held-out: %s]' % (k, fmt_cas_poly(forms[k][0]), forms[k][4]))
    print('   => tau_2 = n1 n2 - Omega with Omega = (C2(mu) - C2(lam) - C2(lam\'))/2 (G2: t S_{2,1} = n1 n2 - Omega), and tau_3, tau_4 = the same differences of the T_k above (n_mu = n1 + n2).  [%.1fs]' % (time.perf_counter() - t0))
    assert all(forms[k][4] for k in forms)
    hdr('G5.1 verification: every Laurent coefficient restricted to M_mu, minus the integer combination of the restricted triple, is the SCALAR MATRIX of the closed form')
    tot_sectors = tot_ops = tot_bad = tot_copies = 0
    tau_bad = 0
    table = []
    for l1, l2 in DICT_PRODUCTS:
        prod = Product.of(l1, l2)
        t0 = time.perf_counter()
        nb = 0
        nsec = 0
        for mu, m in prod.irreps():
            sec = prod.sector(mu)
            assert sec.faithful_viol == 0 and sec.comm_bad == 0
            bad, found = check_dictionary_sector(sec, Fr(1))
            nb += len(bad)
            nsec += 1
            tot_ops += len(found)
            tot_copies += m
            # tau from the Casimir closed forms
            for k in (2, 3, 4):
                if tau_cas(forms, l1, l2, mu, k) != tau_of(l1, l2, mu, k):
                    tau_bad += 1
            n1, n2 = sum(l1), sum(l2)
            Om = (Fr(casimirs(mu)[1]) - casimirs(l1)[1] - casimirs(l2)[1]) / 2
            if tau_of(l1, l2, mu, 2) != n1 * n2 - Om:
                tau_bad += 1
            table.append((l1, l2, mu, m, found, sec))
        tot_sectors += nsec
        tot_bad += nb
        print('%-26s dim %6d: %2d sectors, Sum m = %2d: restricted 18 operators + triple per sector: %d violations of (remainder scalar and = closed form)  [%.1fs]' %
              (pname(l1, l2), prod.R1.dim * prod.R2.dim, nsec, sum(m for _, m in prod.irreps()), nb, time.perf_counter() - t0), flush=True)
    # second and third t (the t-power bookkeeping) on the small products
    for tv in (Fr(31, 35), Fr(-3, 7)):
        nb = ns = 0
        for l1, l2 in REQUIRED_PRODUCTS[:2] + REQUIRED_PRODUCTS[3:5]:
            prod = Product.of(l1, l2)
            for mu, m in prod.irreps():
                bad, found = check_dictionary_sector(prod.sector(mu), tv)
                nb += len(bad)
                ns += 1
        print('t = %-6s: %d sectors of the four small products, violations %d' % (fr_str(tv), ns, nb))
        tot_bad += nb
    print('\nTOTAL: %d sectors (Sum m = %d copies), %d (sector, operator) tests at t = 1: %d violations; tau_k from the Casimir closed forms and tau_2 = n1 n2 - Omega on all sectors: %d violations' %
          (tot_sectors, tot_copies, tot_ops, tot_bad, tau_bad))
    assert tot_bad == 0 and tau_bad == 0 and tot_sectors >= 60
    if write:
        write_dictionary_file(forms, table)
    return forms, table


def write_dictionary_file(forms, table, path=DICT_FILE):
    lines = [
        '# gl4-oper-dictionary.txt -- the dictionary between the accessory parameters (v3, v4a, v4b) of the two-site gl(4) Gaudin oper and ALL 18 Laurent coefficients (WP-G5, oper_gl4.py, 2026-10-22).  Exact.',
        '#',
        '# CONVENTIONS: those of ref/gl4-gaudin-ops.txt (G2).  [E_ij,E_kl] = d_jk E_il - d_li E_kj; two sites V_lam (site 1, pole theta_1) (x) V_lam\' (site 2, pole theta_2); t = theta_1 - theta_2;',
        '#   L(u)_ij = e1_ij/(u-theta_1) + e2_ij/(u-theta_2);  cdet(d_u - L) = d^4 - S_1 d^3 + S_2 d^2 - S_3 d + S_4 (column-ordered, derivative to the right);',
        '#   S_k(u) = sum_{s=1,2} sum_{j=1}^{k} S^{(s)}_{k,j} / (u - theta_s)^j  EXACTLY (poles of order <= k at theta_1, theta_2 only; S_k = O(u^-k) at infinity);  S_1 = |lam|/(u-theta_1) + |lam\'|/(u-theta_2).',
        '#   On a joint eigenvector with eigenvalues (v3, v4a, v4b) of V3 = t S_{3,2}, V4a = t S_{4,3}, V4b = t^2 S_{4,2} (G2: t-independent) this is the scalar 4th order Fuchsian oper (exponents at theta_1:',
        '#   {lam_c + 4 - c}, at theta_2: {lam\'_c + 4 - c}, at infinity: -(mu_c + 4 - c); polynomial solutions of degrees mu_c + 4 - c).',
        '#',
        '# FORMAT.  "DICT k j site tpow a b c scalar" : S^{(site)}_{k,j} = s^(tpow) * (a*v3 + b*v4a + c*v4b + scalar), s = t (site 1) or -t (site 2), tpow = j - k.',
        '#   scalar = T<k>(lam) / T<k>(lamp) (the central top poles), +-tau<k> (tau_k = T_k(mu) - T_k(lam) - T_k(lam\')) or 0.  [this table is an exact consequence of the residue relations and the Casimir values at infinity; it',
        '#   reproduces the integer table of RESULTS-G2 section 6c and was verified against the restricted operators on every sector listed below, at t = 1, 31/35, -3/7]',
        '#   "CASFORM k poly" : T_k(nu) = (-1)^k a_{4-k}(nu) as a polynomial in n = |nu| and the trace-form Casimir eigenvalues C2, C3, C4 of V_nu (C_k = sum e_{i2 i1} e_{i3 i2} ... e_{i1 i_k}).',
        '#       prod_c (theta - y_c) = sum_k a_k theta(theta-1)...(theta-k+1), y_c = nu_c + 4 - c, i.e. a_k = Delta^k F(0)/k!, F(theta) = prod_c (theta - y_c).',
        '#   "SECTOR lam lamp mu m | T2(lam) T3(lam) T4(lam) | T2(lamp) T3(lamp) T4(lamp) | T2(mu) T3(mu) T4(mu) | tau2 tau3 tau4 | triples" : per-sector table of the exact scalars; "triples" = the m joint',
        '#       eigenvalue triples (v3,v4a,v4b) of G2 (25 digits, ascending in the generic combination v3 + 7/19 v4a + 11/13 v4b), i.e. the data of the oper.',
        '#   closed forms: tau_2 = n1 n2 - Omega, Omega = (C2(mu) - C2(lam) - C2(lamp))/2;  tau_3, tau_4 from the CASFORM lines (n_mu = n1 + n2).',
        '#   multiplicity-free family [S,S,0,0] x [P,P,0,0] (N = S + P; Cas2, Cas4 = the su(4) Casimirs of the July-2026 GaudinSolver.wb): v3 = (S+1)(2N(N+2) - Cas2)/2, v4a = 0,',
        '#       v4b = T4(mu) = (Cas2^2 + 2(N+1) Cas2 - 2 Cas4 + 4N(N+2)(N+3))/8, T2(mu) = (4N^2 + 6N - Cas2)/2, T3(mu) = (N+2)(2N(N+2) - Cas2)/2  (RESULTS-G5 section 5; verified on 28 irreps).',
        '#',
    ]
    for (k, j, s), (coef, kind) in sorted(DICT_TABLE.items()):
        sc = {'T': lambda: 'T%d(%s)' % (kind[2], 'lam' if kind[1] == 1 else 'lamp'), 'tau': lambda: '%stau%d' % ('+' if kind[2] > 0 else '-', kind[1]), 'zero': lambda: '0'}[kind[0]]()
        lines.append('DICT\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s' % (k, j, s, j - k, coef[0], coef[1], coef[2], sc))
    for k in (2, 3, 4):
        lines.append('CASFORM\t%d\t%s' % (k, fmt_cas_poly(forms[k][0])))
    for l1, l2, mu, m, found, sec in table:
        js = sec.joint(DPS)
        trip = '; '.join('(' + ', '.join(mp.nstr(x, 25) for x in row) + ')' for row in js['trip'])
        vals = [fr_str(Tchar(l1, k)) for k in (2, 3, 4)] + [fr_str(Tchar(l2, k)) for k in (2, 3, 4)] + [fr_str(Tchar(mu, k)) for k in (2, 3, 4)] + [fr_str(tau_of(l1, l2, mu, k)) for k in (2, 3, 4)]
        lines.append('SECTOR\t%s\t%s\t%s\t%d\t|\t%s\t|\t%s\t|\t%s\t|\t%s\t|\t%s' % (lam_str(l1), lam_str(l2), '(' + ','.join(map(str, mu)) + ')', m, ' '.join(vals[0:3]), ' '.join(vals[3:6]), ' '.join(vals[6:9]), ' '.join(vals[9:12]), trip))
    with open(path, 'w', encoding='utf8', newline='\n') as f:
        f.write('\n'.join(lines) + '\n')
    print('written %s (%d lines)' % (os.path.relpath(path, os.path.join(HERE, '..', '..')), len(lines)))


# ---------------------------------------------------------------------------------------------------- stage exps
def stage_exps():
    hdr('G5.2  EXPONENTS, POLYNOMIAL SOLUTIONS OF EVERY COPY, WRONSKIAN, NEGATIVE TESTS  (mpmath %d digits, exact Fractions where the triple is rational)' % DPS)
    mp.mp.dps = DPS
    tot = Counter()
    worst = dict(gap=None, null_rel=mp.mpf(0), wr=mp.mpf(0), resid=mp.mpf(0))
    neg = dict(true=mp.mpf(0), resid_min=None, smin_min=None, below=0, ntests=0, rnd_resid_min=None, rnd_smin_min=None, rnd_n=0, exact_n=0, exact_rej=0, exact_true_ok=0, exact_true_n=0)
    fails = []
    for l1, l2 in EXPS_PRODUCTS:
        prod = Product.of(l1, l2)
        t0 = time.perf_counter()
        c = Counter()
        for mu, m in prod.irreps():
            sec = prod.sector(mu)
            op = Oper(l1, l2, mu, Fr(1))
            c['sectors'] += 1
            # (i) local indicial polynomials at theta_1, theta_2: v-independent and prod (e - y_c)
            ok_loc = True
            for s_, y_, nu_ in ((1, op.y1, l1), (2, op.y2, l2)):
                L = op.local_indicial(s_)
                a = falling_coeffs(nu_)
                ok_loc = ok_loc and all(L[i] == aff_const(a[NN - i]) for i in range(5))
            # (ii) indicial polynomial at infinity: v-independent and prod (n - (mu_c + 4 - c))
            I0 = op.I0_poly_coeffs()
            a = falling_coeffs(mu)
            ok_inf = all(I0[i] == aff_const(a[NN - i]) for i in range(5))
            # residue relations S_k = O(u^-k) as identities in v; S^{(2)}_{2,1} = -S^{(1)}_{2,1}
            ok_res = all(x == AFF0 for x in op.residue_defects().values()) and all(op.S[(2, 2, 1)][i] == -op.S[(1, 2, 1)][i] for i in range(4))     # sigma_{k,r} = 0 (r < k); S^{(2)}_{2,1} = -S^{(1)}_{2,1}
            c['ok_local'] += ok_loc
            c['ok_inf'] += ok_inf
            c['ok_res'] += ok_res
            if not (ok_loc and ok_inf and ok_res):
                fails.append(('indicial/residue', l1, l2, mu))
            js = sec.joint(DPS)
            ex = js['exact']
            if ex is not None and (l1, l2) in REQUIRED_PRODUCTS[:5]:         # t-independence: the same (t-free) triple at two other values of t = theta_1 - theta_2
                for tv in (Fr(31, 35), Fr(-3, 7)):
                    opt = Oper(l1, l2, mu, tv)
                    for exv in ex:
                        nd_t, _ = poly_solutions_exact(opt, exv)
                        c['tind'] += 1
                        c['tind_ok'] += (nd_t == [4, 3, 2, 1])
            for ci, trip in enumerate(js['trip']):
                c['copies'] += 1
                v = list(trip)
                r = copy_check(op, v, ex[ci] if ex is not None else None, DPS)
                ok_n = r['ok_null']
                ok_w = r['wr_ok']
                ok_e = (r['orders1'] == sorted(op.y1)) and (r['orders2'] == sorted(op.y2))
                c['ok_null'] += ok_n
                c['ok_wr'] += ok_w
                c['ok_orders'] += ok_e
                c['exactcopies'] += (r['mode'] == 'exact')
                if not (ok_n and ok_w and ok_e):
                    fails.append(('copy', l1, l2, mu, ci, ok_n, ok_w, ok_e))
                if r['mode'] == 'mp':
                    if r['gap'] is not None and (worst['gap'] is None or r['gap'] < worst['gap']):
                        worst['gap'] = r['gap']
                    worst['null_rel'] = max(worst['null_rel'], r['null_rel'])
                    worst['resid'] = max(worst['resid'], r['resid'])
                worst['wr'] = max(worst['wr'], r['wr_dev'])
                # negative tests
                ng = negative_tests(op, trip, DPS, seed=hash((mu, ci)) % 1000)
                neg['true'] = max(neg['true'], ng['true'][0])
                for (rr, sm, ns) in ng['perturbed']:
                    neg['ntests'] += 1
                    neg['resid_min'] = rr if neg['resid_min'] is None else min(neg['resid_min'], rr)
                    neg['smin_min'] = sm if neg['smin_min'] is None else min(neg['smin_min'], sm)
                    neg['below'] += (sm < mp.mpf(10) ** (-6))
                if ex is not None:                                   # exact, threshold-free negative tests on the rational triples
                    exv = ex[ci]
                    rng_ = random.Random(hash((mu, ci, 7)) % 1000)
                    tests = []
                    for i in range(3):
                        vv = list(exv)
                        vv[i] = vv[i] + Fr(1, 1000)
                        tests.append(vv)
                    tests.append([Fr(rng_.randint(-300, 300), rng_.randint(1, 9)) for _ in range(3)])
                    neg['exact_n'] += len(tests)
                    neg['exact_rej'] += sum(0 if d1_solution_exists_exact(op, vv) else 1 for vv in tests)
                    neg['exact_true_ok'] += (1 if d1_solution_exists_exact(op, list(exv)) else 0)
                    neg['exact_true_n'] += 1
                rr, sm, ns = ng['random']
                neg['rnd_n'] += 1
                neg['rnd_resid_min'] = rr if neg['rnd_resid_min'] is None else min(neg['rnd_resid_min'], rr)
                neg['rnd_smin_min'] = sm if neg['rnd_smin_min'] is None else min(neg['rnd_smin_min'], sm)
        tot.update(c)
        print('%-26s %2d sectors, %2d copies (%d with rational triple, exact): indicial theta_1,theta_2 v-independent & = prod(e-y): %d/%d; indicial at inf & = prod(n-d): %d/%d; residue relations: %d/%d; '
              'null dims = #{k: d_k <= d_j} for all j: %d/%d; Wronskian ~ (u-th1)^|lam| (u-th2)^|lam\'|: %d/%d; measured vanishing orders at theta_1, theta_2 = {lam_c+4-c}, {lam\'_c+4-c}: %d/%d  [%.1fs]' %
              (pname(l1, l2), c['sectors'], c['copies'], c['exactcopies'], c['ok_local'], c['sectors'], c['ok_inf'], c['sectors'], c['ok_res'], c['sectors'], c['ok_null'], c['copies'], c['ok_wr'], c['copies'], c['ok_orders'], c['copies'], time.perf_counter() - t0), flush=True)
    print('\nTOTAL: %d sectors, %d copies (%d exact/rational): null-space dimensions %d/%d, Wronskians %d/%d, vanishing orders %d/%d; indicial v-independent: local %d/%d, infinity %d/%d, residue relations %d/%d' %
          (tot['sectors'], tot['copies'], tot['exactcopies'], tot['ok_null'], tot['copies'], tot['ok_wr'], tot['copies'], tot['ok_orders'], tot['copies'], tot['ok_local'], tot['sectors'], tot['ok_inf'], tot['sectors'], tot['ok_res'], tot['sectors']))
    print('t-independence: the same exact triples at t = 31/35 and t = -3/7 (theta_2 = 0, theta_1 = t): null dimensions #{k: d_k <= d_j} in %d of %d tests (rational copies of the small products)' % (tot['tind_ok'], tot['tind']))
    assert tot['tind_ok'] == tot['tind']
    print('numerical copies: worst (smallest) gap between kept and dropped singular values %s (criterion >= 1e20), worst dropped-sv / scale %s, worst gauge-solution residual %s, worst Wronskian deviation (rel, incl. exact copies = 0) %s' %
          (mp.nstr(worst['gap'], 4) if worst['gap'] is not None else '-', mp.nstr(worst['null_rel'], 4), mp.nstr(worst['resid'], 4), mp.nstr(worst['wr'], 4)))
    print('NEGATIVE TESTS (v): true triples: least-squares residual of a degree-exactly-d_1 solution <= %s (relative to ||M_1||_F);' % mp.nstr(neg['true'], 3))
    print('   %d perturbations (1e-3 in one component of every known triple): smallest residual %s, smallest sigma_min/sigma_max of M_1 %s (the brief asked for >= 1e-6: %d of %d perturbations are BELOW 1e-6 in this measure;' %
          (neg['ntests'], mp.nstr(neg['resid_min'], 4), mp.nstr(neg['smin_min'], 4), neg['below'], neg['ntests']))
    print('   all %d perturbed triples are nevertheless rejected: no polynomial solution of degree d_1 (smallest residual %s vs %s at the true triples, a separation of > %d orders of magnitude);' %
          (neg['ntests'], mp.nstr(neg['resid_min'], 3), mp.nstr(neg['true'], 3), int(mp.floor(mp.log10(neg['resid_min'] / max(neg['true'], mp.mpf(10) ** (-80)))))))
    print('   %d random triples (uniform in the box |v_i| <= 3 max|true triple|): smallest residual %s, smallest sigma_min/sigma_max %s' % (neg['rnd_n'], mp.nstr(neg['rnd_resid_min'], 4), mp.nstr(neg['rnd_smin_min'], 4)))
    print('   EXACT threshold-free negative tests on the %d rational triples (Fractions: the column u^{d_1} of the pencil is NOT in the span of the other columns <=> no polynomial solution of degree exactly d_1): '
          'true triples have such a solution in %d of %d; the %d tests (each triple + 1/1000 in each component, and one random rational triple) are rejected in %d of %d' %
          (neg['exact_true_n'], neg['exact_true_ok'], neg['exact_true_n'], neg['exact_n'], neg['exact_rej'], neg['exact_n']))
    print('failures: %s' % (fails if fails else 'none'))
    assert not fails and neg['exact_rej'] == neg['exact_n'] and neg['exact_true_ok'] == neg['exact_true_n']
    return tot, neg


# ---------------------------------------------------------------------------------------------------- stage solve
def poly_arrays(polys, nvars=3):
    E, C, own = [], [], []
    for i, f in enumerate(polys):
        for mon, c in f.terms():
            E.append(tuple(mon))
            C.append(float(Fr(int(c.numerator), int(c.denominator))))
            own.append(i)
    return np.array(E, dtype=np.int64), np.array(C), np.array(own, dtype=np.int64), len(polys)


def numeric_search(l1, l2, mu, nstarts=200, seed=0, verbose=True):
    """independent numerical search for solutions of the full (gauge-fixed) condition system: Levenberg-Marquardt (float64) from nstarts random real starts in |v_i| <= 3 max|spectrum|;
       returns dict(found=[distinct solutions], matched=#G2 triples reproduced, extra=#solutions that are NOT G2 triples, converged=#starts that converged, ...)"""
    import scipy.optimize as so
    sec = Product.of(l1, l2).sector(mu)
    op = Oper(l1, l2, mu, Fr(1))
    ctx = RingCtx(0)
    polys = []
    for j in (4, 3, 2, 1):
        polys += frob_conditions(op, j, ctx)[0]
    E, C, own, npoly = poly_arrays(polys)
    g2 = [tuple(float(x) for x in row) for row in sec.joint(DPS)['trip']]
    scale = max(1.0, max(abs(x) for row in g2 for x in row))
    box = 3.0 * scale
    deg = E.sum(axis=1)
    norm = np.bincount(own, weights=np.abs(C) * (box ** deg), minlength=npoly)
    norm[norm == 0] = 1.0

    def F(x):
        mon = np.prod(np.power(x[None, :], E), axis=1)
        return np.bincount(own, weights=C * mon, minlength=npoly) / norm

    def Jf(x):
        J = np.zeros((npoly, 3))
        for a in range(3):
            Ea = E.copy()
            ca = C * E[:, a]
            Ea[:, a] = np.maximum(Ea[:, a] - 1, 0)
            mon = np.prod(np.power(x[None, :], Ea), axis=1)
            J[:, a] = np.bincount(own, weights=ca * mon, minlength=npoly) / norm
        return J
    rng = np.random.default_rng(seed)
    sols = []
    conv = 0
    for _ in range(nstarts):
        x0 = rng.uniform(-box, box, 3)
        try:
            r = so.least_squares(F, x0, jac=Jf, method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=400)
        except Exception:
            continue
        if np.max(np.abs(F(r.x))) < 1e-10:
            conv += 1
            if not any(np.max(np.abs(r.x - s_)) < 1e-5 * scale for s_ in sols):
                sols.append(r.x)
    matched = 0
    for pg in g2:
        if any(np.max(np.abs(np.array(pg) - s_)) < 1e-5 * scale for s_ in sols):
            matched += 1
    extra = [s_ for s_ in sols if not any(np.max(np.abs(np.array(pg) - s_)) < 1e-5 * scale for pg in g2)]
    return dict(m=sec.m, nstarts=nstarts, converged=conv, distinct=len(sols), matched=matched, extra=extra, box=box)


MIN_SUBSET_NOTE = ('{4}: the condition for P_4 alone (degree d_4 = mu_4, no free parameter); {3}: P_3 alone with its u^{d_4} coefficient gauge-fixed to 0 (echelon gauge, as the gl(3) Q0/J reports), etc.; '
                   'a subset J_S is the ideal generated by the gauge-fixed Frobenius conditions (resonance + bottom equations) of the listed j.')


def minimal_exact_subsets(sub):
    ex = [S for S, e in sub.items() if e['equalsK']]
    return [S for S in ex if not any(set(T) < set(S) for T in ex)]


def fmt_subsets(L):
    return ' '.join('{' + ','.join(map(str, sorted(S))) + '}' for S in sorted(L, key=lambda S: (len(S), S))) if L else '-'


SUBSET_PRODUCTS = REQUIRED_PRODUCTS[:5] + EXTRA_DICT_PRODUCTS


def fmt_single(e):
    if e.get('timeout'):
        return 'T/O'
    if e['dim'] is None:
        return '?'
    if e['dim'] == 1:
        return 'curve'
    if e['dim'] == 2:
        return 'surface'
    if e['dim'] == 3:
        return 'all'
    return '%d/%d' % (e['cnt'], e['distinct'])


def stage_solve(write_summary=True):
    hdr('G5.3  THE INDEPENDENT SOLVE: {v : the oper D(v) has polynomial solutions of degrees d_j = mu_j + 4 - j for all j}  =?=  G2\'s joint spectrum')
    print('METHOD (exact).  Frobenius recursion at infinity for each j (p_{d_j} = 1; the coefficient of u^{d_j-q} is determined by the equation c_q unless d_j - q is another exponent d_k, k > j: a RESONANCE, where')
    print('  the equation is a condition on v and p_{d_k} is free -> gauge-fixed to 0, valid when the solution of degree d_k exists), plus the 4 bottom equations (D~ p = 0 has d+5 equations): all polynomial in v = (v3, v4a, v4b);')
    print('  Groebner basis (sympy, grevlex) of the ideal J generated by all conditions of j = 4..1; dim_Q Q[v]/J = #standard monomials; every generator is evaluated EXACTLY on the commuting G2 matrices (V3, V4a, V4b): J subset K,')
    print('  K = ideal of the m joint triples (radical, m points, simple joint spectrum).  So  dim Q[v]/J = m  =>  J = K  =>  Sigma = V(J) = exactly the m triples of G2, each a simple (Jacobian rank 3) solution.  The points of Sigma are')
    print('  also computed from J alone (joint eigenvalues of the multiplication matrices of Q[v]/J, %d digits) and compared with G2.  Subsets of the four conditions (15, small products only): dimension / #solutions / #distinct solutions' % DPS)
    print('  (time-box %.0f s per subset, run in a worker process: T/O = not finished).' % SUBSET_TIMEOUT)
    print('  ' + MIN_SUBSET_NOTE + '\n')
    mp.mp.dps = DPS
    agg = Counter()
    sub_stats = {S: Counter() for S in SUBSETS_J}
    rows = []
    worstdev = mp.mpf(0)
    for l1, l2 in SOLVE_PRODUCTS:
        prod = Product.of(l1, l2)
        t0 = time.perf_counter()
        do_sub = (l1, l2) in SUBSET_PRODUCTS
        print('\n%s (dim %d), %d irreps, Sum m = %d%s' % (pname(l1, l2), prod.R1.dim * prod.R2.dim, len(prod.rs), sum(prod.rs.values()), '' if do_sub else '   [full system only]'), flush=True)
        for mu, m in prod.irreps():
            r = solve_sector(l1, l2, mu, DPS, subsets=do_sub)
            rows.append((l1, l2, r))
            agg['sectors'] += 1
            agg['copies'] += m
            agg['m>=2'] += (m >= 2)
            full = r['full']
            ok = (full['dim'] == 0 and full['cnt'] == m and r['all_vanish_on_G2'] and r['distinct'] == m and all(k == 3 for k in r['jac_ranks']))
            agg['exact_ok'] += ok
            agg['maxm'] = max(agg['maxm'], m)
            dev = r['dev']
            worstdev = max(worstdev, dev)
            digits = ('exact' if dev == 0 else '%d' % int(-mp.floor(mp.log10(dev))))
            extra = ''
            if do_sub:
                sub = r['subsets']
                for S, e in sub.items():
                    sub_stats[S]['sectors'] += 1
                    if e.get('timeout'):
                        sub_stats[S]['timeout'] += 1
                        continue
                    sub_stats[S]['exact'] += e['equalsK']
                    sub_stats[S]['dim_pos'] += (e['dim'] is not None and e['dim'] > 0)
                    sub_stats[S]['finite_extra'] += (e['dim'] == 0 and not e['equalsK'])
                    if m >= 2:
                        sub_stats[S]['sectors2'] += 1
                        sub_stats[S]['exact2'] += e['equalsK']
                extra = ' | exact subsets (minimal): %s | singles %s' % (fmt_subsets(minimal_exact_subsets(sub)), ' '.join('{%d}:%s' % (j, fmt_single(sub[(j,)])) for j in (4, 3, 2, 1)))
            print('  mu=%-13s m=%d nW=%-5d d=%-14s conds %s deg<=%d | Sigma: dim %d, dimQ[v]/J = %s (distinct %s) %s | dev %s digits | Jacobian rank %s%s | %.1fs' %
                  (str(mu), m, r['nW'], str(r['d']), '/'.join(str(r['ncond'][j]) for j in (4, 3, 2, 1)), r['maxdeg'], full['dim'], full['cnt'], r['distinct'], '== m' if ok else '!!', digits, set(r['jac_ranks']),
                   extra, r['t_build'] + r['t_gb']), flush=True)
            if not ok:
                agg['FAIL'] += 1
        print('  [%.1fs]' % (time.perf_counter() - t0))
    print('\nTOTAL: %d sectors (Sum m = %d copies, %d sectors with m >= 2, max m = %d): Sigma = G2\'s joint spectrum EXACTLY (dim_Q Q[v]/J = m = #distinct points, all conditions vanish on the G2 matrices, Jacobian rank 3) in %d of %d sectors; worst deviation of the points of Sigma (computed from J alone) from G2\'s triples: %s' %
          (agg['sectors'], agg['copies'], agg['m>=2'], agg['maxm'], agg['exact_ok'], agg['sectors'], mp.nstr(worstdev, 3)))
    print('\nTHE gl(4) PENCIL-BUG ANALYSIS (small products: %s): for each subset of the four conditions (gauge-fixed as above): in how many sectors it cuts out EXACTLY the spectrum (J_S = K) / sectors with m >= 2 where it does,' % ', '.join(pname(*x) for x in SUBSET_PRODUCTS))
    print('  in how many sectors its solution set has positive dimension, in how many it is finite but has MORE solutions than the spectrum, and in how many the Groebner computation hit the time-box:')
    print('  %-12s %9s %14s %12s %18s %8s' % ('subset', 'exact', 'exact (m>=2)', 'dim > 0', 'finite, too many', 'T/O'))
    for S in SUBSETS_J:
        st = sub_stats[S]
        print('  %-12s %5d/%-3d %8d/%-5d %12d %18d %8d' % ('{' + ','.join(map(str, S)) + '}', st['exact'], st['sectors'], st['exact2'], st['sectors2'], st['dim_pos'], st['finite_extra'], st['timeout']))
    if _WATCHDOG[0] is not None:
        _WATCHDOG[0].close()
    assert agg['FAIL'] == 0
    return rows, agg


def stage_search(sectors=None):
    hdr('G5.3 (numerical, independent of the Groebner computation): random-start Levenberg-Marquardt on the full condition system, >= 200 real starts in the box |v_i| <= 3 max|spectrum|')
    if sectors is None:
        sectors = []
        for l1, l2, thr in [((3, 2, 1, 0), (3, 2, 1, 0), 3), ((4, 2, 1, 0), (4, 2, 1, 0), 4)]:
            for mu, m in Product.of(l1, l2).irreps():
                if m >= thr or (l1 == (3, 2, 1, 0) and mu == (4, 4, 2, 2)):
                    sectors.append((l1, l2, mu))
        sectors.append(((6, 4, 2, 0), (6, 4, 2, 0), (9, 7, 5, 3)))
    out = []
    for l1, l2, mu in sectors:
        t0 = time.perf_counter()
        r = numeric_search(l1, l2, mu, nstarts=(2400 if Product.of(l1, l2).rs[mu] >= 8 else 240), seed=1)
        out.append((l1, l2, mu, r))
        print('  %s mu=%-13s m=%d: %d starts, %d converged (|residual| < 1e-10), %d distinct solutions found, %d of the %d G2 triples reproduced, %d solutions that are NOT G2 triples  [%.1fs]' %
              (pname(l1, l2), str(mu), r['m'], r['nstarts'], r['converged'], r['distinct'], r['matched'], r['m'], len(r['extra']), time.perf_counter() - t0), flush=True)
    return out


# ---------------------------------------------------------------------------------------------------- stage wb
def stage_wb(scan=True):
    hdr('G5.4  THE JULY-2026 gl(4) GAUDIN ODE (GaudinSolver.wb: gaudinODE, Cas2, Cas4, muhat, degQ, wsub) vs OUR OPER')
    cells = load_gs_cells()
    coeffs, body = transcribe_gaudinODE(cells)
    print('TRANSCRIPTION: the cell `ClearAll[gaudinODE]; gaudinODE = ...` of %s is read as text from the notebook JSON and translated mechanically (Derivative[n][Q][z] -> Q_n, ^ -> **, sympy parse);' % os.path.relpath(GS_PATH, os.path.join(HERE, '..', '..')))
    print('  Cas2, Cas4, muhat, xw, Pw, degQ, irreps are transcribed by hand (their Mathematica source is quoted in the docstrings of wb_Cas2, wb_Cas4, wb_muhat, wb_xw, wb_degQ):')
    print('     xw_i = mu_i + 5/2 - i,  Pw_k = sum_i xw_i^k,  Cas2 = Pw_2 - 5,  Cas4 = Pw_4 + 4 Pw_3 - Pw_1 Pw_2 + Pw_2/2 - 2 Pw_1^2 - 10 Pw_1 - 51/4,  muhat_k = mu_k + 4 - k,  degQ_k = (S+P+3) - muhat_k.')
    print('  ODE = sum_n c_n(z) Q^(n)(z) with (factored; S at the site w1, P at w2, leading coefficient 8 (z-w1)^2 (z-w2)^2):')
    for n in (4, 3, 2, 1, 0):
        print('     c_%d = %s' % (n, sp.factor(coeffs[n])))
    # (a) polynomial-solution degrees
    print('\n(a) the notebook\'s claim: polynomial solutions of degrees degQ[S,P,mu]_k exist for every irrep of [S,S,0,0] x [P,P,0,0], and only those degrees (exact: the dimension of the polynomial solutions of degree <= n for n = 0..2(S+P)+8):')
    totirr = okirr = 0
    for (S_, P_) in WB_SP:
        irr = wb_irreps(S_, P_)
        ok = 0
        for mu in irr:
            C2, C4 = wb_Cas2(mu), wb_Cas4(mu)
            for (w1, w2) in ((0, 1), (sp.Rational(2, 7), sp.Rational(-3, 5))):
                subs = {WB_SYMS['S']: S_, WB_SYMS['P']: P_, WB_SYMS['w1']: w1, WB_SYMS['w2']: w2, WB_SYMS['C2']: C2, WB_SYMS['C4']: C4}
                jumps, tot_ = wb_poly_solution_degrees(coeffs, subs, 2 * (S_ + P_) + 8)
                good = (jumps == sorted(wb_degQ(S_, P_, mu)) and tot_ == 4)
                ok += good
                totirr += 1
                okirr += good
        print('   (S,P) = (%d,%d): %d irreps (all multiplicity 1, Sum dim = dim^2 not rechecked here), polynomial solutions exactly at degQ with total dimension 4: %d of %d (w = (0,1) as the notebook and w = (2/7,-3/5))' % (S_, P_, len(irr), ok, 2 * len(irr)))
    assert okirr == totirr
    # (b) gauge scan
    if scan:
        print('\n(b) relation to OUR oper: the ODE must equal  8 (z-w1)^2 (z-w2)^2  Pi^p o D^dagger o Pi^-p  (D^dagger = formal adjoint of our D = d^4 - s1 d^3 + s2 d^2 - s3 d + s4; Pi^p = (z-w1)^p1 (z-w2)^p2), with our D built from the')
        print('    DICTIONARY with UNKNOWN accessory parameters (v3, v4a, v4b) and unknown T_k(mu) (k = 2,3,4); matching all coefficients of z^i in the 5 coefficient functions = a LINEAR system.  Scan of the gauge exponents (p1, p2) and of')
        print('    which of the two sites carries [S,S,0,0] (site 1 = pole w1 = theta_1 or swapped), w = (2/7, -3/5), the top irrep of each (S,P):')
        w1, w2 = sp.Rational(2, 7), sp.Rational(-3, 5)
        for (S_, P_) in [(1, 1), (2, 1)]:
            mu = wb_irreps(S_, P_)[0]
            C2, C4 = wb_Cas2(mu), wb_Cas4(mu)
            cons = []
            tested = 0
            t0 = time.perf_counter()
            for swap in ((False, True) if S_ != P_ else (False,)):
                for p1 in range(0, S_ + P_ + 2):
                    for p2 in range(0, S_ + P_ + 2):
                        unk, lhs = wb_lhs(S_, P_, p1, p2, w1, w2, swap)
                        sol = wb_match_from_lhs(unk, lhs, coeffs, S_, P_, w1, w2, (C2, C4))
                        tested += 1
                        if sol != sp.EmptySet:
                            cons.append((('swapped sites' if swap else 'S at w1'), p1, p2, list(sol)[0]))
            print('    (S,P) = (%d,%d), mu = %s: %d gauge/assignment candidates tested, consistent: %s  [%.1fs]' % (S_, P_, mu, tested, [(c[0], c[1], c[2]) for c in cons], time.perf_counter() - t0), flush=True)
            assert len(cons) == 1 and cons[0][1] == S_ and cons[0][2] == P_ and cons[0][0] == 'S at w1'
    # (c) all irreps: unique consistent solution with p = (S, P); compare with G2's triple and the T_k(mu)
    print('\n(c) all irreps of all five (S,P): with p = (S, P) (S at w1, P at w2) the linear system has a UNIQUE solution; it equals G2\'s exact triple (v3, v4a, v4b) of the (multiplicity-one) sector and T_k(mu) = the central characters;')
    print('    then the two operators coincide COEFFICIENT BY COEFFICIENT (exact polynomial identity in z):')
    w1, w2 = sp.Rational(2, 7), sp.Rational(-3, 5)
    rows = []
    nbad = 0
    for (S_, P_) in WB_SP:
        t0 = time.perf_counter()
        unk, lhs = wb_lhs(S_, P_, S_, P_, w1, w2, False)
        prod = Product.of((S_, S_, 0, 0), (P_, P_, 0, 0))
        for mu in wb_irreps(S_, P_):
            C2, C4 = wb_Cas2(mu), wb_Cas4(mu)
            sol = wb_match_from_lhs(unk, lhs, coeffs, S_, P_, w1, w2, (C2, C4))
            sol = list(sol)
            uniq = (len(sol) == 1 and all(not x.free_symbols for x in sol[0]))
            sec = prod.sector(mu)
            g2 = sec.exact_triples()[0]
            tT = tuple(Tchar(mu, k) for k in (2, 3, 4))
            vals_ = tuple(Fr(int(x.p), int(x.q)) for x in sol[0]) if uniq else None
            ok = uniq and vals_[:3] == g2 and vals_[3:] == tT
            # closed-form dictionary
            cf = wb_dictionary(S_, P_, C2, C4)
            ok_cf = uniq and Fr(int(cf['v3'].p), int(cf['v3'].q)) == g2[0] and cf['v4a'] == 0 == g2[1] and Fr(int(cf['v4b'].p), int(cf['v4b'].q)) == g2[2] and \
                [Fr(int(cf[k].p), int(cf[k].q)) for k in ('T2', 'T3', 'T4')] == list(tT)
            rows.append((S_, P_, mu, C2, C4, g2, tT, ok, ok_cf))
            nbad += (not ok) + (not ok_cf)
        print('   (S,P) = (%d,%d): %d irreps  [%.1fs]' % (S_, P_, len(wb_irreps(S_, P_)), time.perf_counter() - t0), flush=True)
    print('\n    %-7s %-14s %6s %9s | G2: (v3, v4a, v4b)   | T_2, T_3, T_4 (mu)   | unique solution = G2 and T_k: | closed-form dictionary:' % ('(S,P)', 'mu', 'Cas2', 'Cas4'))
    for (S_, P_, mu, C2, C4, g2, tT, ok, ok_cf) in rows:
        print('    (%d,%d)  %-14s %6s %9s | (%s, %s, %s) | (%s, %s, %s) | %-29s | %s' % (S_, P_, str(mu), C2, C4, fr_str(g2[0]), fr_str(g2[1]), fr_str(g2[2]), fr_str(tT[0]), fr_str(tT[1]), fr_str(tT[2]), 'yes' if ok else 'NO', 'yes' if ok_cf else 'NO'))
    print('\n    DICTIONARY (N = S + P; derived symbolically with S, P, w1, w2, C2, C4 all symbolic by `oper_gl4.py wbsym` (150 s), verified here on %d irreps of the five (S,P)):' % len(rows))
    print('        v3  = (S+1) (2N(N+2) - Cas2)/2,    v4a = 0,    v4b = T_4(mu) = (Cas2^2 + 2(N+1) Cas2 - 2 Cas4 + 4N(N+2)(N+3))/8,')
    print('        T_2(mu) = (4N^2 + 6N - Cas2)/2,    T_3(mu) = (N+2)(2N(N+2) - Cas2)/2 = ((N+2)/(S+1)) v3,    T_4(mu) = v4b')
    print('    (Cas2, Cas4 = the notebook\'s su(4) Casimirs: Cas2 = C2^{gl}(mu) - |mu|^2/4 with C2^{gl} = sum mu_i (mu_i + 5 - 2i); the family has mu_1 + mu_4 = mu_2 + mu_3 = N, so the exponents at infinity {mu_c + 4 - c}')
    print('    are symmetric about (N+3)/2, the odd characters are slaved to the even ones, and only (Cas2, Cas4) appear.)')
    print('    verification: unique consistent solution equals G2 and T_k(mu): %d of %d irreps; closed-form dictionary equals G2 triple and T_k: %d of %d' %
          (sum(1 for r in rows if r[7]), len(rows), sum(1 for r in rows if r[8]), len(rows)))
    assert nbad == 0
    # exact operator identity with G2's triple substituted
    nid = 0
    nsub = 0
    for (S_, P_, mu, C2, C4, g2, tT, ok, ok_cf) in rows:
        v = [sp.Rational(x.numerator, x.denominator) for x in g2]
        tm = {k: sp.Rational(tT[k - 2].numerator, tT[k - 2].denominator) for k in (2, 3, 4)}
        s_exprs = our_oper_exprs(S_, P_, w1, w2, v, tm, False)
        G = adjoint_gauge_coeffs(s_exprs, S_, P_, w1, w2)
        z = WB_SYMS['z']
        subs = {WB_SYMS['S']: S_, WB_SYMS['P']: P_, WB_SYMS['w1']: w1, WB_SYMS['w2']: w2, WB_SYMS['C2']: C2, WB_SYMS['C4']: C4}
        good = all(sp.simplify(8 * (z - w1) ** 2 * (z - w2) ** 2 * G[n] - coeffs[n].subs(subs)) == 0 for n in range(5))
        nid += good
        nsub += 1
    print('    exact identity  8 Pi^2 Pi^p D^dagger(v_G2) Pi^-p = gaudinODE(Cas2, Cas4)  as differential operators (all 5 coefficient functions, rational functions of z), G2\'s triple substituted: %d of %d irreps' % (nid, nsub))
    assert nid == nsub
    # (d) degrees and the Wronskian-minor construction of the adjoint solutions
    print('\n(d) the degrees map as polynomial solutions of the adjoint: our polynomial solutions P_j (degrees d_j = mu_j + 4 - j) give the adjoint solutions psi_j = W_3[P_k, k != j] / W_4 and Q_j = Pi^p psi_j =')
    print('    W_3[P_k, k != j] / ((z-w1)^S (z-w2)^P) (W_4 = c (z-w1)^{2S} (z-w2)^{2P}); each is a POLYNOMIAL of degree S+P+3-d_j = (S+P+3) - muhat_j = degQ_j and solves gaudinODE exactly (w = (1,0)):')
    nq = ntot = 0
    for (S_, P_) in WB_SP:
        prod = Product.of((S_, S_, 0, 0), (P_, P_, 0, 0))
        okc = totc = 0
        for mu in wb_irreps(S_, P_):
            sec = prod.sector(mu)
            op = Oper((S_, S_, 0, 0), (P_, P_, 0, 0), mu, Fr(1))
            nd, sols = poly_solutions_exact(op, sec.exact_triples()[0])
            C2, C4 = wb_Cas2(mu), wb_Cas4(mu)
            subs = {WB_SYMS['S']: S_, WB_SYMS['P']: P_, WB_SYMS['w1']: 1, WB_SYMS['w2']: 0, WB_SYMS['C2']: C2, WB_SYMS['C4']: C4}
            for j in range(4):
                W3 = p_trim(wronskian_poly([sols[k] for k in range(4) if k != j], Fr(0)), Fr(0))
                Qj, ok1 = poly_divide_linear_power(W3, Fr(1), S_)
                Qj, ok2 = poly_divide_linear_power(Qj, Fr(0), P_)
                good = ok1 and ok2 and len(Qj) - 1 == (S_ + P_ + 3) - op.d[j]
                if good:
                    res = apply_ode_to_poly(coeffs, subs, Qj)
                    good = all(c == 0 for c in res)
                okc += good
                totc += 1
        nq += okc
        ntot += totc
        print('   (S,P) = (%d,%d): %d of %d adjoint polynomials Q_j (4 per irrep) are polynomials of degree degQ_j and solve the notebook ODE exactly' % (S_, P_, okc, totc))
    assert nq == ntot
    return rows


def stage_wbsym():
    """the 150-second symbolic derivation of the Casimir dictionary (S, P, w1, w2, C2, C4 all symbolic)"""
    hdr('G5.4 (symbolic): solving for (v3, v4a, v4b, T_2(mu), T_3(mu), T_4(mu)) with S, P, w1, w2, C2, C4 symbolic, gauge p = (S, P)')
    S_, P_, w1, w2 = (WB_SYMS[n] for n in ('S', 'P', 'w1', 'w2'))
    coeffs = transcribe_gaudinODE()[0]
    t0 = time.perf_counter()
    unk, lhs = wb_lhs(S_, P_, S_, P_, w1, w2, False)
    z = WB_SYMS['z']
    eqs = []
    for n in range(5):
        num, den = lhs[n]
        poly = sp.Poly(sp.expand(num - coeffs[n] * den), z)
        eqs.extend(poly.coeffs())
    sol = sp.linsolve(eqs, unk)
    for s_ in sol:
        for u, val in zip(unk, s_):
            print('  %s = %s' % (u, sp.factor(val)))
    print('  [%.1fs]' % (time.perf_counter() - t0))
    return sol


def stage_big(nsel=5):
    hdr('G5.3 (extras): the largest multiplicity spaces of G2 (m up to 12): exact solve of the full condition system (no subsets), joint triples of G2 as reference')
    mp.mp.dps = DPS
    out = []
    for l1, l2 in [((5, 3, 1, 0), (5, 3, 1, 0)), ((5, 3, 2, 0), (5, 3, 2, 0)), ((6, 4, 2, 0), (6, 4, 2, 0))]:
        prod = Product.of(l1, l2)
        sel = sorted(prod.rs.items(), key=lambda kv: (-kv[1], kv[0]))[:nsel]
        print('\n%s (dim %d): the %d sectors with the largest multiplicity' % (pname(l1, l2), prod.R1.dim * prod.R2.dim, len(sel)), flush=True)
        for mu, m in sel:
            t0 = time.perf_counter()
            r = solve_sector(l1, l2, mu, DPS, subsets=False)
            ok = (r['full']['dim'] == 0 and r['full']['cnt'] == m and r['all_vanish_on_G2'] and r['distinct'] == m and all(k == 3 for k in r['jac_ranks']))
            dev = r['dev']
            digits = ('exact' if dev == 0 else '%d' % int(-mp.floor(mp.log10(dev))))
            print('  mu=%-14s m=%-2d nW=%-5d d=%-16s conds deg<=%d | Sigma: dim %d, dimQ[v]/J = %s (distinct %s) %s | dev %s digits | Jacobian rank %s | %.1fs' %
                  (str(mu), m, r['nW'], str(r['d']), r['maxdeg'], r['full']['dim'], r['full']['cnt'], r['distinct'], '== m' if ok else '!!', digits, set(r['jac_ranks']), time.perf_counter() - t0), flush=True)
            out.append((l1, l2, mu, m, ok))
    assert all(o[4] for o in out)
    return out


# ---------------------------------------------------------------------------------------------------- CLI
STAGES = ['dict', 'exps', 'solve', 'search', 'wb', 'wbsym', 'big', 'all']


def main(argv):
    cmd = argv[1] if len(argv) > 1 else 'all'
    t0 = time.perf_counter()
    import platform
    import scipy
    print('# oper_gl4.py stage=%s  (Python %s, numpy %s, scipy %s, sympy %s, mpmath %s; %s)' % (cmd, platform.python_version(), np.__version__, scipy.__version__, sp.__version__, mp.__version__, platform.platform()))
    mp.mp.dps = DPS
    if cmd not in STAGES:
        print('usage: python -X utf8 oper_gl4.py <stage>, stage in %s' % STAGES)
        return 1
    if cmd in ('dict', 'all'):
        stage_dict()
    if cmd in ('exps', 'all'):
        stage_exps()
    if cmd in ('solve', 'all'):
        stage_solve()
    if cmd in ('search', 'all'):
        stage_search()
    if cmd in ('wb', 'all'):
        stage_wb()
    if cmd == 'wbsym':
        stage_wbsym()
    if cmd in ('big', 'all'):
        stage_big()
    print('\ntotal time %.1fs' % (time.perf_counter() - t0))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
