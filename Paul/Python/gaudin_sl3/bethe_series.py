"""S1 (round 5, WP-C1): the eps-series of the six Laurent data (H_1,a_1,b_1,H_2,a_2,b_2) of an HW-limit branch from the BETHE SIDE ALONE
(the no-log conditions + the two global rules), including the degenerate cases (Jacobian rank 5: top state and octets) by Fredholm alternative.

Equations (units theta_1 - theta_2 = 1, zeta = eps*phi):  F(x; eps) = 0  with the eight equations
   three no-log conditions per site (nu_general.frobenius(lam)['nolog'] = A, B1, B2 at orders lam2+1, lam1+2, lam1+2; the regular parts of S1,S2,S3 at
   theta_s are written through the other site's data and sigma1, e2, e3 of zeta),
   H_1 + H_2 = sum_{i != j} zeta_i mu_j,   b_1 + b_2 = zeta_1 zeta_2 mu_3 + zeta_1 zeta_3 mu_2 + zeta_2 zeta_3 mu_1     (the 1/u and 1/u^0 matching of the
   exponentials e^{zeta_j u} u^{mu_j} at infinity: q(z) = sum_j mu_j prod_{k != j}(z - zeta_k) = 2|L| z^2 - (H1+H2) z + (b1+b2)).
Order by order  J x_n = -r_n(x_0..x_{n-1}),  J = dF/dx at (x_0, eps = 0) (8 x 6).  rank J = 6: unique.  rank J = 5: x_n = particular + c_n k (k = kernel), and
c_{n-1} is fixed by the solvability N r_n = 0 (N = left null vectors, 3 of them) of the NEXT order (second order if that is empty).
All arithmetic is exact in Q(sqrt D) (class QS); phi is a numeric rational direction; exact polynomial-in-phi forms are obtained by interpolation
on a triangular grid of rational directions (x_n is homogeneous of degree n in phi) and verified on extra directions.
"""
import sys, os, time, pickle, itertools
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sympy as sp

# --------------------------------------------------------------------------------------------- the number field Q(sqrt D)
class QS:
    __slots__ = ('a', 'b', 'D')
    def __init__(self, a=0, b=0, D=5):
        self.a = a if isinstance(a, Fr) else Fr(a)
        self.b = b if isinstance(b, Fr) else Fr(b)
        self.D = D
    def _c(self, o):
        return o if isinstance(o, QS) else QS(o, 0, self.D)
    def __add__(self, o):
        o = self._c(o); return QS(self.a + o.a, self.b + o.b, self.D)
    __radd__ = __add__
    def __sub__(self, o):
        o = self._c(o); return QS(self.a - o.a, self.b - o.b, self.D)
    def __rsub__(self, o):
        return self._c(o) - self
    def __neg__(self):
        return QS(-self.a, -self.b, self.D)
    def __mul__(self, o):
        if isinstance(o, QS):
            return QS(self.a * o.a + self.D * self.b * o.b, self.a * o.b + self.b * o.a, self.D)
        o = Fr(o)
        return QS(self.a * o, self.b * o, self.D)
    __rmul__ = __mul__
    def inv(self):
        n = self.a * self.a - self.D * self.b * self.b
        return QS(self.a / n, -self.b / n, self.D)
    def __truediv__(self, o):
        if isinstance(o, QS):
            return self * o.inv()
        return self * (1 / Fr(o))
    def __rtruediv__(self, o):
        return self.inv() * o
    def iszero(self):
        return self.a == 0 and self.b == 0
    def __eq__(self, o):
        o = self._c(o); return self.a == o.a and self.b == o.b
    def __hash__(self):
        return hash((self.a, self.b))
    def __repr__(self):
        return 'QS(%s,%s)' % (self.a, self.b) if self.b != 0 else '%s' % self.a
    def sym(self):
        return sp.Rational(self.a.numerator, self.a.denominator) + sp.Rational(self.b.numerator, self.b.denominator) * sp.sqrt(self.D)
    def __float__(self):
        return float(self.a) + float(self.b) * self.D ** 0.5

def qs_from_sympy(e, D):
    e = sp.expand(sp.sympify(e))
    s = sp.sqrt(D)
    b = e.coeff(s, 1)
    a = sp.expand(e - b * s)
    if not (a.is_Rational and b.is_Rational):
        e2 = sp.expand(sp.radsimp(e)); b = e2.coeff(s, 1); a = sp.expand(e2 - b * s)
    a = sp.Rational(a); b = sp.Rational(b)
    return QS(Fr(int(a.p), int(a.q)), Fr(int(b.p), int(b.q)), D)

# --------------------------------------------------------------------------------------------- the equations
from nu_general import frobenius, SYMS
_Hs, _As, _Bs, _Ht, _At, _Bt, _SG1, _E2, _E3, _DL = SYMS
HS, AS_, BS_, HT, AT, BT = SYMS[:6]
XS = sp.symbols('X0:6')            # H1, a1, b1, H2, a2, b2
CS = sp.symbols('Cs1 Ce2 Ce3 Ch1 Ct1')   # sigma1, e2, e3, rule constants h1 = sum_{i!=j} phi_i mu_j, t1
CEPS = (1, 2, 3, 1, 2)             # eps powers carried by the constants
_SYS = {}
def system_terms(lam):
    """the eight equations as lists of (Fraction coef, expo6, expo5); rules use mu through the constants Ch1, Ct1 (not here)"""
    lam = tuple(lam[:2])
    if lam in _SYS:
        return _SYS[lam]
    fr = frobenius(lam)
    eqs = []
    for s, t in ((0, 1), (1, 0)):
        d = 1 if s == 0 else -1
        sub = {_Hs: XS[3 * s], _As: XS[3 * s + 1], _Bs: XS[3 * s + 2], _Ht: XS[3 * t], _At: XS[3 * t + 1], _Bt: XS[3 * t + 2],
               _SG1: CS[0], _E2: CS[1], _E3: CS[2], _DL: d}
        for ex in fr['nolog']:
            eqs.append(sp.expand(ex.subs(sub)))
    eqs.append(XS[0] + XS[3] - CS[3])
    eqs.append(XS[2] + XS[5] - CS[4])
    out = []
    for e in eqs:
        P = sp.Poly(e, *XS, *CS)
        terms = []
        for mon, co in P.terms():
            terms.append((Fr(int(co.p), int(co.q)), tuple(mon[:6]), tuple(mon[6:])))
        out.append(terms)
    _SYS[lam] = out
    return out

def consts_of(phi, mu, D):
    """numeric constants (sigma1, e2, e3, h1, t1) for zeta = eps*phi (phi list of Fractions/QS)"""
    p = [QS(x, 0, D) if not isinstance(x, QS) else x for x in phi]
    m = [Fr(x) for x in mu]
    s1 = p[0] + p[1] + p[2]
    e2 = p[0] * p[1] + p[0] * p[2] + p[1] * p[2]
    e3 = p[0] * p[1] * p[2]
    tot = sum(m)
    h1 = p[0] * (tot - m[0]) + p[1] * (tot - m[1]) + p[2] * (tot - m[2])
    t1 = p[0] * p[1] * m[2] + p[0] * p[2] * m[1] + p[1] * p[2] * m[0]
    return [s1, e2, e3, h1, t1]

# --------------------------------------------------------------------------------------------- linear algebra over QS
def rref(M, D):
    """M: list of rows (QS). returns (T, R, pivcols): T invertible (row ops) with T M = R in reduced row echelon form"""
    n, m = len(M), len(M[0])
    A = [[x for x in row] + [QS(1 if i == j else 0, 0, D) for j in range(n)] for i, row in enumerate(M)]
    piv = []
    r = 0
    for c in range(m):
        pr = None
        for i in range(r, n):
            if not A[i][c].iszero():
                pr = i; break
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        iv = A[r][c].inv()
        A[r] = [x * iv for x in A[r]]
        for i in range(n):
            if i != r and not A[i][c].iszero():
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c)
        r += 1
        if r == n:
            break
    T = [row[m:] for row in A]
    R = [row[:m] for row in A]
    return T, R, piv

class Solver:
    """order-by-order solver for one branch at one numeric direction phi.  x0: list of 6 QS (the eps^0 point)."""
    def __init__(self, lam, mu, x0, phi, D, eqterms=None):
        self.lam = tuple(lam[:2]); self.mu = tuple(mu); self.D = D
        self.terms = eqterms or system_terms(self.lam)
        self.C = consts_of(phi, mu, D)
        self.x0 = x0
        self.Z = QS(0, 0, D)
        # constants powers cache: Cpow[c][k]
        self.cp = [[QS(1, 0, D)] for _ in range(5)]
        # pre-contracted terms: (coef QS, expo6, shift)
        self.pt = []
        for eq in self.terms:
            lst = []
            for co, e6, e5 in eq:
                val = QS(co, 0, D); sh = 0
                for c in range(5):
                    if e5[c]:
                        val = val * self.cpow(c, e5[c]); sh += e5[c] * CEPS[c]
                lst.append((val, e6, sh))
            self.pt.append(lst)
        self.xs = [list(x0)]                 # xs[n] = list of 6 QS
        self.mon = {}                         # exps6 -> list of QS coefficients
        # Jacobian at (x0, eps=0)
        self.J = []
        for lst in self.pt:
            row = [self.Z for _ in range(6)]
            for val, e6, sh in lst:
                if sh != 0: continue
                for v in range(6):
                    if e6[v] > 0:
                        e = list(e6); e[v] -= 1
                        t = val * e6[v]
                        for w in range(6):
                            t = t * self.pw(x0[w], e[w])
                        row[v] = row[v] + t
            self.J.append(row)
        self.T, self.R, self.piv = rref(self.J, D)
        self.rank = len(self.piv)
        self.free = [c for c in range(6) if c not in self.piv]
        self.nleft = 8 - self.rank
    def cpow(self, c, k):
        L = self.cp[c]
        while len(L) <= k:
            L.append(L[-1] * self.C[c])
        return L[k]
    def pw(self, x, k):
        r = QS(1, 0, self.D)
        for _ in range(k): r = r * x
        return r
    # --- monomial series
    def mcoef(self, key, n):
        L = self.mon.get(key)
        if L is None:
            L = self.mon[key] = []
        while len(L) <= n:
            m = len(L)
            if sum(key) == 0:
                L.append(QS(1 if m == 0 else 0, 0, self.D)); continue
            i = next(v for v in range(6) if key[v] > 0)
            k2 = list(key); k2[i] -= 1; k2 = tuple(k2)
            tot = self.Z
            for k in range(m + 1):
                xv = self.xs[m - k][i] if (m - k) < len(self.xs) else self.Z
                if xv.iszero(): continue
                tot = tot + self.mcoef(k2, k) * xv
            L.append(tot)
        return L[n]
    def Fn(self, n):
        out = []
        for lst in self.pt:
            tot = self.Z
            for val, e6, sh in lst:
                if sh > n: continue
                tot = tot + val * self.mcoef(e6, n - sh)
            out.append(tot)
        return out
    def invalidate(self, n):
        for L in self.mon.values():
            del L[n:]
    def kernel_basis(self):
        ks = []
        for f in self.free:
            k = [self.Z] * 6
            k[f] = QS(1, 0, self.D)
            for i, pc in enumerate(self.piv):
                k[pc] = -self.R[i][f]
            ks.append(k)
        return ks
    def step(self, n, cvec=None):
        """compute x_n given x_0..x_{n-1} (already in self.xs).  cvec: values of the free coordinates (list, len = 6-rank; default 0).
        returns (xn, g) with g = the solvability residual N r_n (list of nleft QS); x_n is appended to self.xs"""
        assert len(self.xs) == n
        self.invalidate(n)
        self.xs.append([self.Z] * 6)
        r = self.Fn(n)                                   # with x_n = 0
        Tr = [sum((self.T[i][j] * r[j] for j in range(8)), self.Z) for i in range(8)]
        g = Tr[self.rank:]
        xn = [self.Z] * 6
        cvec = cvec or [self.Z] * len(self.free)
        for f, c in zip(self.free, cvec):
            xn[f] = c if isinstance(c, QS) else QS(c, 0, self.D)
        for i, pc in enumerate(self.piv):
            v = -Tr[i]
            for f in self.free:
                v = v - self.R[i][f] * xn[f]
            xn[pc] = v
        self.xs[n] = xn
        self.invalidate(n)
        return xn, g
    def snapshot(self):
        return len(self.xs)
    def truncate(self, n):
        del self.xs[n:]
        self.invalidate(n)
    def resid(self, n):
        """full residual F_n with the current x_n (should vanish once x_n is the right solution)"""
        return self.Fn(n)

# --------------------------------------------------------------------------------------------- polynomial-in-c helpers (univariate over QS)
def polyfit_c(vals, pts):
    """interpolate the values (QS) at the points (QS) -> coefficient list (low to high), exact Lagrange via Newton"""
    n = len(pts)
    coef = list(vals)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i - 1]) / (pts[i] - pts[i - j])
    # convert Newton form to monomial
    poly = [coef[-1]]
    for i in range(n - 2, -1, -1):
        # poly = poly*(x - pts[i]) + coef[i]
        new = [poly[0] * (-pts[i])] + [poly[k - 1] - poly[k] * pts[i] if k < len(poly) else poly[k - 1] for k in range(1, len(poly) + 1)]
        new[0] = new[0] + coef[i]
        poly = new
    return poly

def ptrim(p):
    p = list(p)
    while p and p[-1].iszero():
        p.pop()
    return p

def pdivmod(a, b):
    a = ptrim(a); b = ptrim(b)
    q = [QS(0, 0, b[0].D) for _ in range(max(len(a) - len(b) + 1, 1))]
    while len(a) >= len(b) and a:
        c = a[-1] / b[-1]; d = len(a) - len(b)
        q[d] = c
        a = ptrim([x - (c * b[i - d] if 0 <= i - d < len(b) else 0) for i, x in enumerate(a)])
    return q, a

def pgcd(a, b):
    a = ptrim(a); b = ptrim(b)
    while b:
        _, r = pdivmod(a, b)
        a, b = b, r
    if a:
        iv = a[-1].inv()
        a = [x * iv for x in a]
    return a

def proots_linear(g):
    """g: monic poly of degree 1 -> its root"""
    return -g[0] / g[1]


# --------------------------------------------------------------------------------------------- the degenerate driver
class Degenerate(Exception):
    pass

def roots_in_K(G, D):
    """roots in Q(sqrt D) of the polynomial G (list of QS, low to high), as QS; multiplicities dropped"""
    c = sp.Symbol('c')
    expr = sum(g.sym() * c ** i for i, g in enumerate(G))
    P = sp.Poly(expr, c, extension=sp.sqrt(D)) if D not in (1, 4, 9) else sp.Poly(expr, c)
    out = []
    for f, mult in P.factor_list()[1]:
        if f.degree() == 1:
            a, b = f.all_coeffs()
            out.append(qs_from_sympy(sp.simplify(-b / a), D))
    return out

def _resolve_candidates(S, m, N, maxdeg, look, verbose, zero, pts):
    """candidate values (list of QS) of the first pending kernel coordinate c_{m+1}, and its order a = m+1"""
    D = S.D; nk = len(S.free)
    def rerun(nn, trials):
        S.truncate(m + 1)
        for k in range(m + 1, nn):
            S.step(k, trials.get(k, zero))
        return S.step(nn, zero)[1]
    logs = []
    for n in range(m + 1, m + 2 + look):
        pend = list(range(m + 1, n))
        if not pend:
            g = rerun(n, {})
            if not all(x.iszero() for x in g):
                raise Degenerate('solvability fails with no pending unknown at order %d' % n)
            continue
        a = pend[0]
        gv = [rerun(n, {a: [c]}) for c in pts]
        dep = False
        if len(pend) > 1:
            g1 = rerun(n, {a: [QS(3, 0, D)]})
            g2 = rerun(n, {a: [QS(3, 0, D)], **{k: [QS(7 + k, 0, D)] for k in pend[1:]}})
            dep = any(not (x - y).iszero() for x, y in zip(g1, g2))
        polys = [ptrim(polyfit_c([gv[i][j] for i in range(len(pts))], pts)) for j in range(len(gv[0]))]
        degs = [len(p) - 1 for p in polys]
        logs.append((n, tuple(pend), tuple(degs), dep))
        if verbose: print('  order %d pending %s: deg of g in c_%d = %s, depends on later c: %s' % (n, pend, a, degs, dep))
        if any(d > maxdeg for d in degs):
            raise Degenerate('degree too high %s' % (degs,))
        if dep:
            continue
        nz = [p for p in polys if p]
        if not nz:
            continue
        G = nz[0]
        for p in nz[1:]:
            G = pgcd(G, p)
        if len(G) < 2:
            raise Degenerate('no common root at order %d: degrees %s' % (n, degs))
        cands = roots_in_K(G, D)
        if not cands:
            raise Degenerate('no root in K at order %d (degree-%d factor)' % (n, len(G) - 1))
        return a, cands, logs
    raise Degenerate('could not resolve c_%d within lookahead %d' % (m + 1, look))

def solve_all(S, N, verbose=False, maxdeg=4, look=3, _m=0, _log=None):
    """all formal solutions x_0..x_N of the branch (list of (xs, log)); a branch is DFS-ed over the K-rational roots of the solvability polynomial.
    rank-6 branches have a single solution."""
    D = S.D; nk = len(S.free); Z = QS(0, 0, D)
    log = list(_log or [])
    if nk == 0:
        for n in range(1, N + 1):
            xn, g = S.step(n, None)
            if not all(x.iszero() for x in g):
                raise Degenerate('rank-6 branch: solvability fails at order %d' % n)
        return [(list(S.xs[:N + 1]), log)]
    zero = [Z] * nk
    pts = [QS(i, 0, D) for i in range(-2, maxdeg + 2)]
    m = _m
    if m >= N:
        return [(list(S.xs[:N + 1]), log)]
    a, cands, lg = _resolve_candidates(S, m, N, maxdeg, look, verbose, zero, pts)
    out = []
    for c in cands:
        S.truncate(m + 1)
        for k in range(m + 1, a):
            S.step(k, zero)
        S.step(a, [c] + [Z] * (nk - 1))
        if verbose: print('  -> c_%d = %s (candidate %d of %d)' % (a, c, cands.index(c) + 1, len(cands)))
        try:
            sub = solve_all(S, N, verbose, maxdeg, look, a, log + lg + [('c', a, c)])
        except Degenerate as ex:
            if verbose: print('     candidate c_%d = %s dies: %s' % (a, c, ex))
            continue
        out += sub
    if not out:
        raise Degenerate('all candidates for c_%d die' % a)
    return out

def solve_series(S, N, **kw):
    sols = solve_all(S, N, **kw)
    if len(sols) != 1:
        raise Degenerate('%d surviving solutions' % len(sols))
    xs, log = sols[0]
    S.truncate(N + 1)
    return xs, log

# --------------------------------------------------------------------------------------------- the eps = 0 system
def order0_equations(lam):
    """the eight equations at eps = 0 as sympy expressions in XS (H1,a1,b1,H2,a2,b2)"""
    eqs = []
    for lst in system_terms(lam):
        e = 0
        for co, e6, e5 in lst:
            if sum(e5) > 0: continue
            t = sp.Rational(co.numerator, co.denominator)
            for v in range(6): t *= XS[v] ** e6[v]
            e += t
        eqs.append(sp.expand(e))
    return eqs

def order0_points(lam, verbose=False):
    """all solutions of the eps = 0 system (exact, sympy).  Eliminate H2 = -H1, b2 = -b1 (the global rules) first; for lam2 = 1 also b_s from the first condition."""
    lam = tuple(lam[:2])
    eqs = order0_equations(lam)
    sub = {XS[3]: -XS[0], XS[5]: -XS[2]}
    eq6 = [sp.expand(e.subs(sub)) for e in eqs[:6]]
    vars_ = [XS[2], XS[4], XS[1], XS[0]]
    if lam[1] == 1:
        bsol = {XS[2]: sp.solve(eq6[0], XS[2])[0]}
        bs2 = sp.solve(eq6[3].subs(bsol).subs(XS[5], 0), XS[2]) if False else None
        # b2 = -b1 and A_2: c b2 = a2 (H2 - a2) with H2 = -H1 -> consistency a1(H - a1) + a2(-H - a2) = 0 ... keep as equation
        rest = [sp.expand(sp.numer(sp.together(e.subs(bsol)))) for e in eq6[1:]]
        vars_ = [XS[4], XS[1], XS[0]]
        eq6 = rest
    G = sp.groebner(eq6, *vars_, order='lex')
    if verbose:
        for g in G.exprs: print(sp.factor(g))
    sols = sp.solve(list(G.exprs), vars_, dict=True)
    out = []
    for s in sols:
        d = dict(s)
        if lam[1] == 1:
            d[XS[2]] = sp.simplify(bsol[XS[2]].subs(d))
        d[XS[3]] = -d[XS[0]]; d[XS[5]] = -d[XS[2]]
        out.append([sp.nsimplify(sp.simplify(d[v])) for v in XS])
    return out


def branches_w4(lam):
    """the eps = 0 branch data from W4's exact files (units theta12 = 1): dict (mu, copy) -> (D, x0 as sympy list [H1,a1,b1,H2,a2,b2]); copy 0 = smaller a_1.
    (They are checked against the eps = 0 equations in check_order0.)"""
    from nu_general import load_w4
    w = load_w4(tuple(lam[:2]) + (0,))
    out = {}
    for mu, s in w['sectors'].items():
        for b in s['branches']:
            D = b['D0'] or 5
            def cv(p):
                return sp.Rational(p[0].numerator if hasattr(p[0], 'numerator') else p[0], p[0].denominator if hasattr(p[0], 'denominator') else 1) + \
                    sp.Rational(p[1].numerator if hasattr(p[1], 'numerator') else p[1], p[1].denominator if hasattr(p[1], 'denominator') else 1) * sp.sqrt(D)
            x0 = [cv(b['eta'][0]), cv(b['alpha'][0]), cv(b['beta'][0]), cv(b['eta'][1]), cv(b['alpha'][1]), cv(b['beta'][1])]
            out[(mu, b['copy'] - 1)] = (D, x0)
    return out

def check_order0(lam, br):
    """residual of the eps = 0 equations and the rank of J for each branch"""
    eqs = order0_equations(tuple(lam[:2]))
    res = {}
    for key, (D, x0) in br.items():
        sub = dict(zip(XS, x0))
        r = [sp.simplify(e.subs(sub)) for e in eqs]
        J = sp.Matrix([[sp.diff(e, v) for v in XS] for e in eqs]).subs(sub).applyfunc(sp.simplify)
        res[key] = (all(x == 0 for x in r), J.rank())
    return res


# --------------------------------------------------------------------------------------------- comparison with W5 / W4 files
PHIA = (Fr(1, 3), Fr(-2, 7), Fr(5, 11))
PHIB = (Fr(-2, 5), Fr(3, 8), Fr(7, 11))
TH12 = Fr(31, 35)

def load_e4(lam, fname=None):
    """parse the E4_SERIES records of W5's laurent-series-<lam>.txt (or W4's untwisted-exact): dict (mu, copy, phi-tuple) -> (D0, H[s][n], a[s][n], b[s][n]) as QS lists"""
    from nu_general import _parse, W4DIR
    fname = fname or os.path.join(W4DIR, 'laurent-series-%d%d%d.txt' % (lam[0], lam[1], 0))
    out = {}
    for line in open(fname, encoding='utf8'):
        if not line.startswith('E4_SERIES'):
            continue
        f = line.rstrip('\n').split('\t')
        mu = tuple(_parse(f[1])); cp = int(f[2]); D0 = int(f[3]); phi = tuple(_parse(f[4]))
        D = D0 or 5
        rows = [_parse(f[6]), _parse(f[7]), _parse(f[8])]
        data = []
        for r in rows:
            data.append([[QS(Fr(p[0]), Fr(p[1]), D) for p in site] for site in r])
        out[(mu, cp - 1, tuple(Fr(x) for x in phi))] = (D, data[0], data[1], data[2])
    return out

def to_theta(xs, th12=TH12):
    """series in units theta12 = 1 -> series at general theta12: H,a ~ th12^(n-1), b ~ th12^(n-2)"""
    out = []
    for n, x in enumerate(xs):
        out.append([x[0] * th12 ** (n - 1), x[1] * th12 ** (n - 1), x[2] * th12 ** (n - 2), x[3] * th12 ** (n - 1), x[4] * th12 ** (n - 1), x[5] * th12 ** (n - 2)])
    return out

def compare_w5(lam, N, phis=(PHIA, PHIB), fname=None, verbose=True):
    """run the Bethe-side solver for every branch and phi in the file and compare ALL six data series (both sites) with the file exactly, orders 0..N"""
    ref = load_e4(lam, fname)
    br = branches_w4(lam)
    res = {}
    for (mu, cp), (D, x0) in br.items():
        for phi in phis:
            if (mu, cp, tuple(phi)) not in ref: continue
            S = Solver(lam, mu, [qs_from_sympy(v, D) for v in x0], list(phi), D)
            t = time.time()
            xs, log = solve_series(S, N)
            ser = to_theta(xs)
            D2, Hs_, As_, Bs_ = ref[(mu, cp, tuple(phi))]
            nm = min(N, len(Hs_[0]) - 1)
            bad = 0
            for n in range(nm + 1):
                for si in range(2):
                    for q, R in ((0, Hs_), (1, As_), (2, Bs_)):
                        if not (ser[n][3 * si + q] - R[si][n]).iszero(): bad += 1
            res[(mu, cp, phi)] = (bad, nm, time.time() - t)
            if verbose: print('  %s copy %d phi %s: rank %d, orders 0..%d, mismatching entries %d (%.1fs)' % (mu, cp, [str(x) for x in phi], S.rank, nm, bad, time.time() - t), flush=True)
    return res


# --------------------------------------------------------------------------------------------- exact polynomial-in-phi series (interpolation on a triangular grid)
P1, P2, P3 = sp.symbols('p1 p2 p3')
GRIDA = [Fr(1, 3) + Fr(k, 5) + Fr(1, 7 * (k + 2)) for k in range(14)]
GRIDB = [Fr(-2, 7) + Fr(k, 3) + Fr(1, 11 * (k + 3)) for k in range(14)]

def _monoms(n):
    return [(i, j, n - i - j) for i in range(n + 1) for j in range(n + 1 - i)]

def interp_homog(points, vals, n, D):
    """points: list of phi triples (phi1 = 1 normalised not required), vals: QS values of a degree-n homogeneous polynomial -> dict monomial -> QS"""
    mons = _monoms(n)
    assert len(points) >= len(mons)
    A = []
    for ph in points[:len(mons)]:
        A.append([QS(ph[0] ** m[0] * ph[1] ** m[1] * ph[2] ** m[2], 0, D) for m in mons])
    # solve A c = v by Gauss-Jordan with augmented column (use rref on [A | v])
    M = [row + [v] for row, v in zip(A, vals[:len(mons)])]
    T, R, piv = rref(M, D)
    assert piv == list(range(len(mons))), 'singular interpolation grid'
    return {m: R[i][len(mons)] for i, m in enumerate(mons)}

def poly_of(coefs):
    e = 0
    for (i, j, k), c in coefs.items():
        if c.iszero(): continue
        e += c.sym() * P1 ** i * P2 ** j * P3 ** k
    return sp.expand(e)

def symbolic_series(lam, mu, D, x0, N, verbose=False):
    """exact series x_n(phi), n = 0..N, as sympy polynomials in p1,p2,p3 (list over n of the six data), by interpolation over a triangular grid of rational directions
    (phi = (1, y, z)), checked at extra directions."""
    npts = len(_monoms(N))
    pts = []
    for i in range(N + 1):
        for j in range(N + 1 - i):
            pts.append((Fr(1), GRIDA[i], GRIDB[j]))
    runs = []
    t0 = time.time()
    for ph in pts:
        S = Solver(lam, mu, x0, list(ph), D)
        xs, log = solve_series(S, N)
        runs.append(xs)
    if verbose: print('    %d solver runs %.1fs' % (len(pts), time.time() - t0), flush=True)
    series = []
    for n in range(N + 1):
        # use the points of the first |monoms(n)| in a unisolvent subset: triangular grid of degree n sits inside the degree-N grid (i + j <= n)
        sub = [(ph, runs[k][n]) for k, (ph, _) in enumerate([(p, 0) for p in pts]) if True]
        idx = [k for k in range(len(pts)) if (lambda i, j: i + j <= n)(*divmod_idx(k, N))]
        row = []
        for q in range(6):
            cf = interp_homog([pts[k] for k in idx], [runs[k][n][q] for k in idx], n, D)
            row.append((cf, poly_of(cf)))
        series.append(row)
    # verification at extra directions
    extra = [(Fr(2, 3), Fr(-5, 4), Fr(3, 5)), (Fr(-1, 2), Fr(7, 3), Fr(11, 13))]
    for ph in extra:
        S = Solver(lam, mu, x0, list(ph), D)
        xs, log = solve_series(S, N)
        for n in range(N + 1):
            for q in range(6):
                cf = series[n][q][0]
                v = QS(0, 0, D)
                for (i, j, k), c in cf.items():
                    v = v + c * (ph[0] ** i * ph[1] ** j * ph[2] ** k)
                assert (v - xs[n][q]).iszero(), ('interpolation check fails', n, q)
    return [[poly for (cf, poly) in row] for row in series]

def divmod_idx(k, N):
    """inverse of the enumeration k <-> (i, j) with i in 0..N, j in 0..N-i used in symbolic_series"""
    c = 0
    for i in range(N + 1):
        for j in range(N + 1 - i):
            if c == k: return i, j
            c += 1


def stage_symbolic(lam, default_N, special=None, out=None):
    """exact symbolic series of every branch of lam; special = {(mu,copy): N}; saves bethe_series_<lam>_sym.pkl: {(mu,copy): {'H1':[poly...], ...}} (sympy strings, units theta12 = 1)"""
    special = special or {}
    br = branches_w4(lam)
    res = {}
    names = ('H1', 'a1', 'b1', 'H2', 'a2', 'b2')
    for key, (D, x0) in br.items():
        N = special.get(key, default_N)
        t = time.time()
        ser = symbolic_series(lam, key[0], D, [qs_from_sympy(v, D) for v in x0], N)
        res[key] = {nm: [str(ser[n][q]) for n in range(N + 1)] for q, nm in enumerate(names)}
        print('  %s copy %d: N = %d  (%.0fs)' % (key[0], key[1], N, time.time() - t), flush=True)
    fn = out or os.path.join(HERE, 'bethe_series_%d%d0_sym.pkl' % lam)
    with open(fn, 'wb') as f:
        pickle.dump(res, f)
    return res

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage', default='sym210')
    a = ap.parse_args()
    if a.stage == 'sym210':
        stage_symbolic((2, 1), 8)
    elif a.stage == 'sym310':
        stage_symbolic((3, 1), 7, {((6, 1, 1), 0): 9, ((5, 3, 0), 0): 9})
    elif a.stage == 'sym320':
        stage_symbolic((3, 2), 7)


def compare_symbolic(lam, pkl=None, fname=None):
    """compare my exact symbolic series (pkl of stage_symbolic) with W5's E4_SYMBOLIC records (orders 0..3, theta12 = 31/35): returns (#checked, #mismatch)"""
    import re
    from nu_general import W4DIR
    fname = fname or os.path.join(W4DIR, 'laurent-series-%d%d0.txt' % lam)
    mine = pickle.load(open(pkl or os.path.join(HERE, 'bethe_series_%d%d0_sym.pkl' % lam), 'rb'))
    names = ('H1', 'a1', 'b1', 'H2', 'a2', 'b2')
    chk = bad = 0
    th = sp.Rational(31, 35)
    for line in open(fname, encoding='utf8'):
        if not line.startswith('E4_SYMBOLIC'): continue
        f = line.rstrip('\n').split('\t')
        mu = tuple(int(x) for x in f[1].strip('{}').split(',')); cp = int(f[2]) - 1; nm = f[4]
        strs = re.findall(r'"([^"]*)"', f[5])
        key = (mu, cp)
        if key not in mine: continue
        for n, s in enumerate(strs):
            e = sp.sympify(s.replace('Sqrt', 'sqrt').replace('[', '(').replace(']', ')').replace('^', '**'))
            q = names.index(nm)
            ex = 1 if nm[0] in 'Ha' else 2
            mymine = sp.sympify(mine[key][nm][n]) * th ** (n - ex)
            chk += 1
            if sp.simplify(sp.expand(e - mymine)) != 0:
                bad += 1
                print('  mismatch', key, nm, n, e, mymine)
    return chk, bad
