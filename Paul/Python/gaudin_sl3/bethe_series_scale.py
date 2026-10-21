"""WP-S (round 6, agent S): the eps-series of the Laurent data (H_1,a_1,b_1,H_2,a_2,b_2) of the HW-limit branches of the two-site sl(3) Gaudin chain, Lambda = [lam1,lam2,0] on
both sites, from the BETHE SIDE ALONE, at scale: mpmath arithmetic (no Q(sqrt D)), series arithmetic on the Frobenius recursion (no expanded polynomials), order 0 from the untwisted
Q-system (exact pencil in b), degenerate (rank-5) orders by a numeric Fredholm/candidate-following solver.   Brief: SPEC-bethe-series-scale.md.   Report: RESULTS-bethe-series-scale.md.

Layers (all in this file):
  1. representation theory   decompose(lam)            irreps mu of V_Lambda x V_Lambda with multiplicities (Racah-Speiser on GT weights)
  2. order 0                 eta_of, q_of, pencil_poly, order0_branches, field_degrees
                             x_0(b) = (eta, q-b/2, b, -eta, -q-b/2, -b); b = common roots of the three echelon Baxter systems (exact rational arithmetic)
  3. series engine           Engine                    Frobenius recursion on truncated eps-series, order by order: the eight equations F_n, clone/truncate, nu-series
  4. numeric degenerate solver  BranchSolver           Jacobian + rank, rank 6 / rank 5 (Fredholm alternative, lag 1..lambda1, candidate following, field-membership tie-break)
  5. validation / runs       load_ref, validate, run_family, check_order0, test_homogeneity, timing_study
  6. depth analysis          kernel_T, G_terms_from_T, leading_analysis, kB_rule   (nu-series x exact kernel -> leading polynomial, N_B, depth)

Usage (from Paul/Python/gaudin_sl3):   python -X utf8 bethe_series_scale.py --stage STAGE --lam l1 l2 [--N N] [--dps D] [--phis A|B|AB] [--cache file.pkl]
  --stage order0    pencil roots vs W4 / agent R / the eight equations at eps = 0 (counts vs Racah-Speiser)
  --stage validate  acceptance 1: 210, 310, 320 against W5's exact series (all branches, both directions, orders 0..N)
  --stage family    all branches of Lambda x Lambda to order N (cache pickle, log lines); --stage depth  N_B and depth per branch (exact kernel: general_kernel_rep / block theorem on the one-row kernels / agent K's kernel_Theta_*.npz)
  --stage timing    cumulative time to fix order n on the top state / the largest-multiplicity irrep / a rank-6 irrep
Units: theta_1 - theta_2 = 1 (H, a at order n scale as theta12^(n-1), b as theta12^(n-2)); zeta = eps*phi.  Everything real (the physical branches are real).
Gotcha: importing the earlier modules (nu_general, general_kernel, ...) resets mp.mp.dps; every entry point here re-applies the working precision after such imports.
"""
import sys, os, time, pickle, itertools, random, argparse
from fractions import Fraction as Fr
from math import factorial
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
import mpmath as mp
from mpmath import mpf
import sympy as sp

DPS = 80          # working precision (decimal digits) of the series engine
TOL_EXP = 40      # minimum exponent of the zero tests |x| <= 10^-tol_exp * scale (BranchSolver default: max(TOL_EXP, 0.6 * dps))


def setdps(d):
    global DPS
    DPS = d
    mp.mp.dps = d


mp.mp.dps = DPS


def mpq(x):
    """exact rational (int/Fraction) -> mpf at the current precision"""
    if isinstance(x, Fr):
        return mpf(x.numerator) / x.denominator
    return mpf(x)


# ============================================================================================ 1. representation theory
def gt_weights(lam):
    """multiset of weights (w1,w2,w3) of the GL(3) irrep lam (GT patterns of gl3_reps)"""
    from gl3_reps import gt_patterns, _weight
    cnt = {}
    for p in gt_patterns(tuple(lam)):
        w = tuple(int(_weight(k, p)) for k in (1, 2, 3))
        cnt[w] = cnt.get(w, 0) + 1
    return cnt


_PERMS = [((0, 1, 2), 1), ((1, 0, 2), -1), ((0, 2, 1), -1), ((1, 2, 0), 1), ((2, 0, 1), 1), ((2, 1, 0), -1)]


_WMC = {}


def weights_cached(lam):
    lam = tuple(lam)
    if lam not in _WMC:
        _WMC[lam] = gt_weights(lam)
    return _WMC[lam]


def mult_rs(lam, mu, wm=None):
    """multiplicity of the irrep mu in V_lam x V_lam (Racah-Speiser / Klimyk: sum_s sgn(s) mult_lam[s(mu+rho) - rho - lam]); wm = weight multiplicities of V_lam"""
    wm = wm or gt_weights(lam)
    rho = (2, 1, 0)
    tot = 0
    mr = tuple(mu[i] + rho[i] for i in range(3))
    for perm, sg in _PERMS:
        w = tuple(mr[perm[i]] - rho[i] - lam[i] for i in range(3))
        tot += sg * wm.get(w, 0)
    return tot


def decompose(lam):
    """list of (mu, multiplicity) of V_lam x V_lam, mu dominant with |mu| = 2|lam|, sorted like W4: mu1 descending then mu2 descending"""
    lam = tuple(lam)
    wm = gt_weights(lam)
    tot = 2 * sum(lam)
    out = []
    for m1 in range(tot, -1, -1):
        for m2 in range(min(m1, tot - m1), -1, -1):
            m3 = tot - m1 - m2
            if m3 < 0 or m3 > m2:
                continue
            m = mult_rs(lam, (m1, m2, m3), wm)
            if m > 0:
                out.append(((m1, m2, m3), m))
    return out


# ============================================================================================ 2. order 0: closed forms and the Q-system pencil
def c2_cas(l):
    return l[0] * (l[0] + 2) + l[1] * l[1] + l[2] * (l[2] - 2)


def eta_of(lam, mu):
    """H_1^0 * theta12 = |Lambda|^2 - (C2(mu) - 2 C2(Lambda))/2   (C2(l) = sum l_i (l_i + 4 - 2i))"""
    L = sum(lam)
    return Fr(L * L) - Fr(c2_cas(mu) - 2 * c2_cas(lam), 2)


def q_of(mu):
    """(a_1 - a_2)/2 * theta12 = (mu1+2)(mu2+1) mu3 / 2   (product of the exponents at infinity / 2)"""
    return Fr((mu[0] + 2) * (mu[1] + 1) * mu[2], 2)


def pencil_coeffs(lam, mu):
    """exact coefficient functions of the three-term recursion alpha_{N-1} p_{N-1} + beta_N p_N + rho(N+1) p_{N+1} = 0 of the untwisted oper (units theta12 = 1),
    beta_N = beta0_N + b/2.  Returns (alpha, beta0, rho, exps) with exps = (mu3, mu2+1, mu1+2)."""
    l1, l2 = lam[0], lam[1]
    L = l1 + l2
    c2 = l2 * (l1 + 1)
    H = eta_of(lam, mu)
    q = q_of(mu)
    alpha = lambda j: Fr(j * (j - 1) * (j - 2) - 2 * L * j * (j - 1)) + (H + 2 * c2) * j - 2 * q
    beta0 = lambda j: Fr(-2 * j * (j - 1) * (j - 2) + 3 * L * j * (j - 1)) - (H + 2 * c2) * j + q
    rho = lambda j: Fr(j * (j - l2 - 1) * (j - l1 - 2))
    return alpha, beta0, rho, (mu[2], mu[1] + 1, mu[0] + 2)


def _fr_det(M):
    """exact determinant of a square list-of-lists of Fractions (Gaussian elimination)"""
    n = len(M)
    A = [list(r) for r in M]
    det = Fr(1)
    for c in range(n):
        p = None
        for i in range(c, n):
            if A[i][c] != 0:
                p = i
                break
        if p is None:
            return Fr(0)
        if p != c:
            A[c], A[p] = A[p], A[c]
            det = -det
        det *= A[c][c]
        inv = 1 / A[c][c]
        for i in range(c + 1, n):
            f = A[i][c] * inv
            if f != 0:
                A[i] = [a - f * b for a, b in zip(A[i], A[c])]
    return det


def _interp(xs, ys):
    """Newton interpolation (exact) -> coefficient list low..high of the polynomial through (xs, ys)"""
    n = len(xs)
    coef = list(ys)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i - 1]) / (xs[i] - xs[i - j])
    poly = [coef[-1]]
    for i in range(n - 2, -1, -1):
        new = [Fr(0)] * (len(poly) + 1)
        for k, c in enumerate(poly):
            new[k + 1] += c
            new[k] -= c * xs[i]
        new[0] += coef[i]
        poly = new
    return poly


_BSYM = sp.Symbol('b')


def pencil_poly(lam, mu, nproj=3, seed=12345, which=(1, 2, 3)):
    """the polynomial in b (sympy Poly over QQ) whose roots are the admissible values of b = b_1 for the sector mu: for each exponent d_j (j = 1,2,3) the echelon
    system (P_j monic of degree d_j, no u^{d_i} term for i < j) has d_j + 1 equations in d_j - (j-1) unknowns, affine in b; b is admissible iff [A|v] has a null vector
    with last entry 1.  Polynomial = gcd over `nproj` exact determinants of random integer projections; the result is the gcd over j in `which`."""
    alpha, beta0, rho, exps = pencil_coeffs(lam, mu)
    rng = random.Random(seed)
    polys = []
    for j in which:
        dj = exps[j - 1]
        fixed = set(exps[:j - 1])
        U = [i for i in range(dj) if i not in fixed]               # free unknown coefficients
        nr = dj + 1                                                  # equations N = 0..dj
        nc = len(U) + 1
        colidx = {i: c for c, i in enumerate(U)}
        colidx[dj] = nc - 1                                          # p_dj = 1 (the inhomogeneity)
        def aug(bv):
            M = [[Fr(0)] * nc for _ in range(nr)]
            for N in range(nr):
                if N >= 1 and (N - 1) in colidx:
                    M[N][colidx[N - 1]] += alpha(N - 1)
                if N in colidx:
                    M[N][colidx[N]] += beta0(N) + bv / 2
                if (N + 1) in colidx:
                    M[N][colidx[N + 1]] += rho(N + 1)
            return M
        gs = None
        for t in range(nproj):
            R = [[Fr(rng.randint(-9, 9)) for _ in range(nr)] for _ in range(nc)]
            xs = [Fr(k) for k in range(nc + 2)]
            ys = []
            for bv in xs:
                A = aug(bv)
                RA = [[sum(R[i][k] * A[k][c] for k in range(nr)) for c in range(nc)] for i in range(nc)]
                ys.append(_fr_det(RA))
            co = _interp(xs, ys)
            P = sp.Poly([sp.Rational(c.numerator, c.denominator) for c in reversed(co)], _BSYM, domain='QQ')
            gs = P if gs is None else sp.gcd(gs, P)
        polys.append(gs)
    g = polys[0]
    for P in polys[1:]:
        g = sp.gcd(g, P)
    return g


def order0_branches(lam, mu, dps=None, return_poly=False):
    """branches of the sector mu at eps = 0 (units theta12 = 1): list of x_0 = [H1,a1,b1,H2,a2,b2] (mpf) sorted by a_1 ascending (copy 1 = smaller a_1, as W4/W5)
    and the exact squarefree pencil polynomial"""
    dps = dps or DPS
    g = pencil_poly(lam, mu)
    gs = sp.Poly(sp.sqf_part(g.as_expr()), _BSYM, domain='QQ')
    eta = eta_of(lam, mu)
    q = q_of(mu)
    mp.mp.dps = dps + 20
    if gs.degree() == 0:
        roots = []
    else:
        roots = [r for r in gs.nroots(n=dps + 10, maxsteps=500)]
    mp.mp.dps = dps
    out = []
    for r in roots:
        im = abs(mpf(sp.im(r)))
        if im > mpf(10) ** (-(dps // 2)):
            raise ValueError('complex pencil root %s for %s' % (r, mu))
        b = mpf(str(sp.re(r)))
        H = mpq(eta)
        qq = mpq(q)
        out.append([H, qq - b / 2, b, -H, -qq - b / 2, -b])
    out.sort(key=lambda x: x[1])
    if return_poly:
        return out, gs
    return out


# ============================================================================================ 3. the series engine
def field_degrees(gs, roots_b):
    """for every pencil root b (mpf) the degree over Q of its minimal polynomial = degree of the irreducible factor of the (squarefree) pencil polynomial gs having b as a root:
    the copies of a multiplicity-m sector live in number fields K_b = Q(b) of this degree (1 for rational b)"""
    facs = [sp.Poly(f, _BSYM, domain='QQ') for f, _ in sp.factor_list(gs.as_expr())[1]]
    cf = [[mpq(Fr(int(c.p), int(c.q))) for c in f.all_coeffs()] for f in facs]
    out = []
    for b in roots_b:
        best = min(range(len(facs)), key=lambda i: abs(_peval(cf[i][::-1], b)) / max(1, sum(abs(c) * max(1, abs(b)) ** k for k, c in enumerate(cf[i][::-1]))))
        out.append(facs[best].degree())
    return out


def pairs_of(lam):
    """non-forced pairs (k,l), 0 <= k < l <= lam1+1, not both <= lam2 (same order as nu_general.pairs_of)"""
    l1, l2 = lam[0], lam[1]
    return [(k, l) for k in range(l1 + 2) for l in range(k + 1, l1 + 2) if not (l <= l2)]


def conv(P, Q, n):
    """order-n coefficient of the product of two series (lists of mpf, len >= n+1)"""
    return mp.fdot(P[:n + 1], Q[n::-1] if n > 0 else Q[:1])


def series_mul(P, Q, N):
    return [conv(P, Q, n) for n in range(N + 1)]


class Engine:
    """eps-series engine for one branch: the unknown series X[q] (q = H1,a1,b1,H2,a2,b2), the Frobenius solutions u, v (both sites) and the eight equations, computed ORDER BY ORDER.
    compute(n, xn): requires orders < n done; appends x_n, runs the order-n pass, returns F_n (8 numbers: A,B1,B2 of site 1, A,B1,B2 of site 2, rule H, rule b).
    truncate(n): forget all orders >= n.  All arithmetic mpf (current mp.dps)."""
    def __init__(self, lam, mu, phi, th12=1):
        self.lam = (lam[0], lam[1])
        self.l1, self.l2 = lam[0], lam[1]
        self.mu = tuple(mu)
        self.L = self.l1 + self.l2
        self.c2 = self.l2 * (self.l1 + 1)
        self.mres1, self.mres2 = self.l2 + 1, self.l1 + 2
        self.th12 = Fr(th12)
        p = [mpq(Fr(x)) if not isinstance(x, mpf) else x for x in phi]
        self.phi = p
        self.s1 = p[0] + p[1] + p[2]
        self.e2 = p[0] * p[1] + p[0] * p[2] + p[1] * p[2]
        self.e3 = p[0] * p[1] * p[2]
        tot = sum(mu)
        self.h1 = p[0] * (tot - mu[0]) + p[1] * (tot - mu[1]) + p[2] * (tot - mu[2])
        self.t1 = p[0] * p[1] * mu[2] + p[0] * p[2] * mu[1] + p[1] * p[2] * mu[0]
        K = self.mres2 + 2
        self.wg1, self.wg2, self.dl = [], [], []
        for s in (0, 1):
            d = mpq(self.th12) * (1 if s == 0 else -1)
            self.dl.append(d)
            self.wg1.append([(-1) ** k * d ** (-k - 1) for k in range(K)])
            self.wg2.append([(-1) ** k * (k + 1) * d ** (-k - 2) for k in range(K)])
        self.rho = [mpf(m * (m - self.l2 - 1) * (m - self.l1 - 2)) for m in range(self.mres2 + 1)]
        self.clear()

    def clear(self):
        self.X = [[] for _ in range(6)]
        self.n = 0
        # C[s][sol][m] = list over orders;  Wl/Xl/Yl[s][sol][m] = scalar-weighted sums (series)
        self.C = [[[[] for _ in range(self.mres2 + 1)] for _ in range(2)] for _ in range(2)]
        self.W = [[[[] for _ in range(self.mres2 + 1)] for _ in range(2)] for _ in range(2)]
        self.Xs = [[[[] for _ in range(self.mres2 + 1)] for _ in range(2)] for _ in range(2)]
        self.Y = [[[[] for _ in range(self.mres2 + 1)] for _ in range(2)] for _ in range(2)]
        self.res = [[[] for _ in range(3)] for _ in range(2)]      # A, B1, B2 per site
        self.F = []

    def truncate(self, n):
        assert n <= self.n
        for q in range(6):
            del self.X[q][n:]
        for s in (0, 1):
            for sol in (0, 1):
                for m in range(self.mres2 + 1):
                    del self.C[s][sol][m][n:]
                    del self.W[s][sol][m][n:]
                    del self.Xs[s][sol][m][n:]
                    del self.Y[s][sol][m][n:]
            for k in range(3):
                del self.res[s][k][n:]
        del self.F[n:]
        self.n = n

    def clone(self):
        """independent copy of the current state (parameters shared, arrays copied by slicing)"""
        e = Engine.__new__(Engine)
        e.__dict__.update(self.__dict__)
        e.X = [row[:] for row in self.X]
        e.C = [[[row[:] for row in self.C[s][sol]] for sol in range(2)] for s in range(2)]
        e.W = [[[row[:] for row in self.W[s][sol]] for sol in range(2)] for s in range(2)]
        e.Xs = [[[row[:] for row in self.Xs[s][sol]] for sol in range(2)] for s in range(2)]
        e.Y = [[[row[:] for row in self.Y[s][sol]] for sol in range(2)] for s in range(2)]
        e.res = [[row[:] for row in self.res[s]] for s in range(2)]
        e.F = self.F[:]
        return e

    def compute(self, n, xn):
        assert n == self.n, (n, self.n)
        for q in range(6):
            self.X[q].append(xn[q] if isinstance(xn[q], mpf) else mpf(xn[q]))
        zero = mpf(0)
        l2, mres1, mres2, L, c2 = self.l2, self.mres1, self.mres2, self.L, self.c2
        for s in (0, 1):
            t = 1 - s
            Hs, As, Bs, Ht, At, Bt = self.X[3 * s], self.X[3 * s + 1], self.X[3 * s + 2], self.X[3 * t], self.X[3 * t + 1], self.X[3 * t + 2]
            wg1, wg2 = self.wg1[s], self.wg2[s]
            for sol in (0, 1):
                Cs = self.C[s][sol]
                Ws, Xs_, Ys = self.W[s][sol], self.Xs[s][sol], self.Y[s][sol]
                # seeds
                Cs[0].append(mpf(1) if (sol == 0 and n == 0) else zero)
                for m in range(1, mres2 + 1):
                    if sol == 1 and m <= l2:
                        Cs[m].append(zero); Ws[m].append(zero); Xs_[m].append(zero); Ys[m].append(zero)
                        continue
                    if sol == 1 and m == mres1:
                        Cs[m].append(mpf(1) if n == 0 else zero)
                        Ws[m].append(zero); Xs_[m].append(zero); Ys[m].append(zero)
                        continue
                    # scalar-weighted sums at order n
                    S = zero; Wn = zero; Zn = zero; Xn = zero; Yn = zero
                    for i in range(2, m):
                        ci = Cs[i][n]
                        if ci != 0:
                            S += wg1[m - 1 - i] * (i * (i - 1)) * ci
                    for i in range(1, m - 1):
                        ci = Cs[i][n]
                        if ci != 0:
                            Wn += wg1[m - 2 - i] * i * ci
                            Zn += wg2[m - 2 - i] * i * ci
                    for i in range(0, m - 2):
                        ci = Cs[i][n]
                        if ci != 0:
                            Xn += wg1[m - 3 - i] * ci
                            Yn += wg2[m - 3 - i] * ci
                    Ws[m].append(Wn); Xs_[m].append(Xn); Ys[m].append(Yn)
                    R = -L * S + c2 * Zn
                    if n >= 1 and m >= 3:
                        R -= self.s1 * ((m - 1) * (m - 2)) * Cs[m - 1][n - 1]
                    if n >= 2 and m >= 3:
                        R += self.e2 * (m - 2) * Cs[m - 2][n - 2]
                    if n >= 3 and m >= 3:
                        R -= self.e3 * Cs[m - 3][n - 3]
                    if m >= 3:
                        R += conv(Ht, Ws[m], n) - conv(Bt, Xs_[m], n) - conv(At, Ys[m], n)
                    if m >= 2:
                        R -= conv(Bs, Cs[m - 2], n)
                    c1 = Cs[m - 1]
                    # ((m-1) H_s - a_s) c_{m-1}
                    R += (m - 1) * conv(Hs, c1, n) - conv(As, c1, n)
                    if m == mres1 and sol == 0:
                        self.res[s][0].append(R)
                        Cs[m].append(zero)
                    elif m == mres2:
                        self.res[s][1 + sol].append(R)
                        Cs[m].append(zero)
                    else:
                        Cs[m].append(-R / self.rho[m])
        X = self.X
        F = [self.res[0][0][n], self.res[0][1][n], self.res[0][2][n], self.res[1][0][n], self.res[1][1][n], self.res[1][2][n],
             X[0][n] + X[3][n] - (self.h1 if n == 1 else 0), X[2][n] + X[5][n] - (self.t1 if n == 2 else 0)]
        self.F.append(F)
        self.n = n + 1
        return F

    def nu_series(self, N=None):
        """nu_kl(theta_s) series for both sites: dict s -> list over pairs_of(lam) of lists (orders 0..N)"""
        N = (self.n - 1) if N is None else N
        out = []
        pr = pairs_of(self.lam)
        den = factorial(self.l2 + 1)
        for s in (0, 1):
            Cu, Cv = self.C[s][0], self.C[s][1]
            lst = []
            for (k, l) in pr:
                a = series_mul(Cu[k], Cv[l], N)
                b = series_mul(Cu[l], Cv[k], N)
                f = mpf(factorial(k) * factorial(l)) / den
                lst.append([f * (x - y) for x, y in zip(a, b)])
            out.append(lst)
        return out


# ============================================================================================ 4. the numeric degenerate solver
def rref_mp(M, tol):
    """M: list of rows (mpf).  Returns (T, R, piv): T invertible (row operations) with T M = R in reduced row echelon form; pivots chosen by partial pivoting,
    a column is 'free' if its best remaining pivot is <= tol * max|M|."""
    n, m = len(M), len(M[0])
    scale = max(abs(x) for row in M for x in row) or mpf(1)
    A = [list(row) + [mpf(1) if i == j else mpf(0) for j in range(n)] for i, row in enumerate(M)]
    piv = []
    r = 0
    for c in range(m):
        if r == n:
            break
        pr = max(range(r, n), key=lambda i: abs(A[i][c]))
        if abs(A[pr][c]) <= tol * scale:
            continue
        A[r], A[pr] = A[pr], A[r]
        iv = 1 / A[r][c]
        A[r] = [x * iv for x in A[r]]
        for i in range(n):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c)
        r += 1
    T = [row[m:] for row in A]
    R = [row[:m] for row in A]
    return T, R, piv


class Degenerate(Exception):
    pass


def _trim(p, tol):
    p = list(p)
    while p and abs(p[-1]) <= tol:
        p.pop()
    return p


def _peval(p, x):
    r = mpf(0)
    for c in reversed(p):
        r = r * x + c
    return r


def _pder(p):
    return [k * c for k, c in enumerate(p)][1:]


def _fit(pts, vals):
    """polynomial through the points (coefficients low..high) by Newton divided differences"""
    n = len(pts)
    coef = list(vals)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i - 1]) / (pts[i] - pts[i - j])
    poly = [coef[-1]]
    for i in range(n - 2, -1, -1):
        new = [mpf(0)] * (len(poly) + 1)
        for k, c in enumerate(poly):
            new[k + 1] += c
            new[k] -= c * pts[i]
        new[0] += coef[i]
        poly = new
    return poly


class BranchSolver:
    """order-by-order solver of one branch at one direction phi.  x0: the six eps^0 values (mpf).  Rank-5 Jacobians (one kernel direction) are handled by the Fredholm
    alternative with candidate following (SPEC section 5): the pending kernel coordinate c_{m+1} is fixed by the first later order whose solvability residual g_n(c_{m+1}) is
    not identically zero (lag = n - m - 1); the polynomial g_n(c) is reconstructed from trial runs on cloned engines, its real common roots are the candidates (followed depth first,
    pruned by the consistency tests of the later orders).  Rank < 5 raises NotImplementedError with the data."""
    def __init__(self, lam, mu, x0, phi, th12=1, tol_exp=None, maxdeg=None, look=None, verbose=False, field_deg=None):
        self.lam, self.mu, self.phi = tuple(lam[:2]), tuple(mu), tuple(phi)
        self.field_deg = field_deg
        self.selected_by = None
        self.eng = Engine(lam, mu, phi, th12)
        # zero tests |g| <= 10^-tol_exp * scale: the pruning tolerance must be MUCH tighter than the loosest value that works for few candidates -- at [14,7,0] (28,13,1)
        # a near-solution with relative inconsistency ~1e-41 survived the old 1e-40; default now 0.6 * dps digits (the true solution's noise is ~1e-(dps-10))
        self.tol_exp = tol_exp if tol_exp is not None else max(TOL_EXP, int(0.6 * mp.mp.dps))
        self.tol = mpf(10) ** (-self.tol_exp)
        self.maxdeg = maxdeg if maxdeg is not None else 2 * lam[0] + 4
        self.look = look if look is not None else 2 * lam[0] + 4
        self.verbose = verbose
        self.x0 = [mpf(x) for x in x0]
        self.Z6 = [mpf(0)] * 6
        F0 = self.eng.compute(0, self.x0)
        self.res0 = max(abs(f) for f in F0)
        # Jacobian J[:, j] = F_1(x_1 = e_j) - F_1(0)
        F1 = self.eng.compute(1, self.Z6)
        cols = []
        for j in range(6):
            self.eng.truncate(1)
            e = list(self.Z6)
            e[j] = mpf(1)
            Fj = self.eng.compute(1, e)
            cols.append([a - b for a, b in zip(Fj, F1)])
        self.eng.truncate(1)
        self.J = [[cols[j][i] for j in range(6)] for i in range(8)]
        self.T, self.R, self.piv = rref_mp(self.J, mpf(10) ** (-35))
        self.rank = len(self.piv)
        self.free = [c for c in range(6) if c not in self.piv]
        self.nk = 6 - self.rank
        self.nleft = 8 - self.rank
        self.passes = 0
        self.t_start = time.perf_counter()
        self.t_fix = {}
        self.lags = {}
        self.info = {}
        self.ndeaths = 0

    def svals(self):
        A = mp.matrix(self.J)
        return [mp.mpf(s) for s in mp.svd_r(A, compute_uv=False)]

    # ---- one order ---------------------------------------------------------------------------------------------
    def kernel_vec(self, f):
        k = [mpf(0)] * 6
        k[f] = mpf(1)
        for i, pc in enumerate(self.piv):
            k[pc] = -self.R[i][f]
        return k

    def _step_on(self, eng, n, cvec=None):
        """compute x_n on engine eng (free coordinates = cvec, default 0); returns (x_n, g, sc) with g the solvability residual N r_n (3 numbers), sc its natural scale"""
        assert eng.n == n, (eng.n, n)
        r = eng.compute(n, self.Z6)
        eng.truncate(n)
        Tr = [mp.fdot(self.T[i], r) for i in range(8)]
        g = Tr[self.rank:]
        sc = max(mpf(1), max(abs(x) for x in r))
        xn = [mpf(0)] * 6
        cvec = cvec or [mpf(0)] * self.nk
        for f, c in zip(self.free, cvec):
            xn[f] = mpf(c)
        for i, pc in enumerate(self.piv):
            v = -Tr[i]
            for f in self.free:
                v -= self.R[i][f] * xn[f]
            xn[pc] = v
        eng.compute(n, xn)
        self.passes += 2
        return xn, g, sc

    def step(self, n, cvec=None):
        return self._step_on(self.eng, n, cvec)

    def _g_from(self, m, n, trials):
        """solvability residual g_n for the base state x_0..x_m, orders m+1..n-1 stepped with trial free coordinates (dict order -> value, default 0); runs on a clone"""
        e = self.eng.clone()
        assert e.n == m + 1
        for k in range(m + 1, n):
            self._step_on(e, k, [mpf(trials.get(k, 0))])
        r = e.compute(n, self.Z6)
        self.passes += 1
        Tr = [mp.fdot(self.T[i], r) for i in range(8)]
        return Tr[self.rank:], max(mpf(1), max(abs(x) for x in r))

    # ---- candidates for the pending coordinate c_{m+1} ------------------------------------------------------------
    def resolve(self, m):
        if self.nk != 1:
            raise NotImplementedError('kernel dimension %d (rank %d)' % (self.nk, self.rank))
        tol = self.tol
        a = m + 1
        base = self.eng
        assert base.n == m + 1
        # (i) consistency at order a: no pending unknown
        e0 = base.clone()
        r = e0.compute(a, self.Z6)
        self.passes += 1
        Tr = [mp.fdot(self.T[i], r) for i in range(8)]
        g = Tr[self.rank:]
        sc = max(mpf(1), max(abs(x) for x in r))
        if any(abs(x) > tol * sc for x in g):
            raise Degenerate('solvability fails with no pending unknown at order %d (|g| = %s)' % (a, mp.nstr(max(abs(x) for x in g), 3)))
        # (ii) lockstep scan of three trial runs: first order n > a at which g_n is not identically zero
        pts = [mpf('0.37'), mpf('-1.21'), mpf('2.43')]
        runs = []
        for c in pts:
            e = base.clone()
            self._step_on(e, a, [c])
            runs.append(e)
        vals = [[] for _ in pts]           # vals[j][k - a - 1] = g_k at trial j
        scs = []
        n = a
        nstar_list = []
        while n < a + self.look:
            n += 1
            gs, ss = [], []
            for j, e in enumerate(runs):
                xn, g, sc = self._step_on(e, n, None)
                gs.append(g)
                ss.append(sc)
            sc = max(ss)
            for j in range(3):
                vals[j].append(gs[j])
            scs.append(ss)
            nonzero = any(abs(x) > tol * sc for g in gs for x in g)
            if self.verbose:
                print('   order %d (lag %d): g %s' % (n, n - a, 'nonzero' if nonzero else 'identically zero'), flush=True)
            if not nonzero:
                continue
            # dependence on the later pending coordinates c_{a+1}..c_{n-1}
            dep = False
            if n - a > 1:
                g1, _ = self._g_from(m, n, {a: mpf(3)})
                tr2 = {a: mpf(3)}
                tr2.update({k: mpf(7 + k) for k in range(a + 1, n)})
                g2, _ = self._g_from(m, n, tr2)
                dep = any(abs(x - y) > tol * sc for x, y in zip(g1, g2))
            if dep:
                if self.verbose:
                    print('   order %d: g depends on later pending coordinates, skipped' % n, flush=True)
                continue
            # polynomial reconstruction in c_a from the three scan points and extra trial points
            ptl = list(pts)
            vl = [vals[j][n - a - 1] for j in range(3)]
            sl = list(scs[n - a - 1])
            extra = []
            for kk in range(2, 40):
                extra += [mpf(kk), mpf(-kk)]
            polys = None
            for deg in range(1, self.maxdeg + 1):
                while len(ptl) < deg + 2:
                    c = extra[len(ptl) - 3]
                    g, s2 = self._g_from(m, n, {a: c})
                    ptl.append(c)
                    vl.append(g)
                    sl.append(s2)
                scl = max(sl[:deg + 2])
                ok = True
                cand = []
                for j in range(len(vl[0])):
                    p = _fit(ptl[:deg + 1], [v[j] for v in vl[:deg + 1]])
                    pred = _peval(p, ptl[deg + 1])
                    if abs(pred - vl[deg + 1][j]) > tol * max(scl, abs(pred)):
                        ok = False
                        break
                    cand.append(p)
                if ok:
                    polys = cand
                    break
            if polys is None:
                raise Degenerate('degree of g in c_%d exceeds %d at order %d' % (a, self.maxdeg, n))
            sc = max(sl)
            degs = [len(_trim(p, tol * sc)) - 1 for p in polys]
            if self.verbose:
                print('   order %d (lag %d): degrees of g in c_%d = %s' % (n, n - a, a, degs), flush=True)
            nzj = [(j, _trim(p, tol * sc)) for j, p in enumerate(polys)]
            nzj = [(j, p) for j, p in nzj if p]
            nz = [p for j, p in nzj]
            if any(len(p) == 1 for p in nz):
                raise Degenerate('no common root at order %d: a component is a non-zero constant' % n)
            cands = self._common_roots(nz, sc)
            if not cands:
                raise Degenerate('no real common root at order %d (degrees %s)' % (n, degs))
            cands = [self._refine(m, n, a, c, nzj) for c in cands]
            self.lags[a] = n - a
            self.info[a] = dict(lag=n - a, degs=degs, ncand=len(cands), mults=list(getattr(self, 'last_mults', [])))
            return a, cands
        raise Degenerate('could not resolve c_%d within lookahead %d' % (a, self.look))

    def _refine(self, m, n, a, c, nzj, iters=3):
        """Newton steps on the solvability residual g_n(c_a) itself (evaluated by trial runs, accurate at the root); the slope is taken from the fitted polynomials (nzj = list of
        (component index, fitted polynomial)).  The fit is spread over c ~ O(1..5) where g is huge, so the roots of the fitted polynomial lose digits at deep lag / high degree;
        the direct evaluation at the root does not."""
        tol = mpf(10) ** (-(mp.mp.dps - 10))
        for _ in range(iters):
            g, _s = self._g_from(m, n, {a: c})
            if max(abs(x) for x in g) <= mpf(10) ** (-(mp.mp.dps - 12)) * _s:      # already at the noise floor of the residual
                break
            jj, p = max(nzj, key=lambda t: abs(_peval(_pder(t[1]), c)))
            d = _peval(_pder(p), c)
            if d == 0:
                break
            dc = g[jj] / d
            c = c - dc
            if abs(dc) <= tol * max(1, abs(c)):
                break
        return c

    def _common_roots(self, nz, sc):
        """real common roots of the component polynomials nz (lists low..high).  Roots of the lowest-degree component (polyroots); clusters of nearby roots (a multiple root split by the
        rounding noise into points ~ noise^(1/m) apart) are replaced by their centroid, which is accurate to the noise itself; simple roots are polished by Newton on the component of
        largest slope; a root is kept if it is a root (to tolerance) of every component."""
        nz = sorted(nz, key=len)
        p0 = nz[0]
        mult = {}
        if len(p0) == 2:
            roots = [-p0[0] / p0[1]]
        else:
            rr = mp.polyroots(p0[::-1], maxsteps=5000, extraprec=4 * mp.mp.prec, error=False)
            cands = []
            for r in rr:
                r = mp.mpmathify(r)
                if abs(mp.im(r)) <= mpf(10) ** (-20) * max(1, abs(r)):
                    cands.append(mp.re(r))
            cands.sort()
            clusters = []
            for r in cands:
                if clusters and abs(r - clusters[-1][-1]) <= mpf(10) ** (-20) * max(1, abs(r)):
                    clusters[-1].append(r)
                else:
                    clusters.append([r])
            roots = [sum(c) / len(c) for c in clusters]
            mult = {i: len(c) for i, c in enumerate(clusters)}
        out = []
        self.last_mults = []
        conv_tol = mpf(10) ** (-(mp.mp.dps - 8))
        for idx, r in enumerate(roots):
            if mult.get(idx, 1) == 1:
                for _ in range(200):                                    # Newton on the component of largest slope until converged
                    best = max(nz, key=lambda p: abs(_peval(_pder(p), r)))
                    d = _peval(_pder(best), r)
                    if d == 0:
                        break
                    dr = _peval(best, r) / d
                    r = r - dr
                    if abs(dr) <= conv_tol * max(1, abs(r)):
                        break
            resid = max(abs(_peval(p, r)) / max(sc, sum(abs(c) * max(1, abs(r)) ** k for k, c in enumerate(p))) for p in nz)
            if resid <= mpf(10) ** (-(self.tol_exp - 10)):
                dup = [i for i, o in enumerate(out) if abs(r - o[0]) <= mpf(10) ** (-30) * max(1, abs(r))]
                if dup:
                    if resid < out[dup[0]][1]:                           # keep the more accurate of two coinciding roots
                        out[dup[0]] = (r, resid, mult.get(idx, 1))
                else:
                    out.append((r, resid, mult.get(idx, 1)))
        self.last_mults = [o[2] for o in out]
        out = [o[0] for o in out]
        return out

    def solve_all(self, N, m=0, log=None):
        """all formal solutions up to order N: list of (xs, log); DFS over candidate roots (real), pruned by the consistency tests of the next orders"""
        eng = self.eng
        log = list(log or [])
        if self.nk == 0:
            for n in range(1, N + 1):
                xn, g, sc = self.step(n)
                self.t_fix[n] = time.perf_counter() - self.t_start
                if any(abs(x) > self.tol * sc for x in g):
                    raise Degenerate('rank-6 branch: solvability fails at order %d' % n)
            return [([list(eng.X[q]) for q in range(6)], log)]
        if m >= N:
            return [([list(eng.X[q][:N + 1]) for q in range(6)], log)]
        a, cands = self.resolve(m)
        out = []
        for c in cands:
            eng.truncate(m + 1)
            self.step(a, [c])
            self.t_fix[a] = time.perf_counter() - self.t_start
            if self.verbose:
                print('  -> c_%d = %s (candidate %d of %d)' % (a, mp.nstr(c, 20), cands.index(c) + 1, len(cands)), flush=True)
            try:
                sub = self.solve_all(N, a, log + [(a, c, len(cands))])
            except Degenerate as ex:
                self.ndeaths += 1
                if self.verbose:
                    print('     candidate c_%d = %s dies: %s' % (a, mp.nstr(c, 20), ex), flush=True)
                continue
            out += sub
        if not out:
            raise Degenerate('all %d candidates for c_%d die' % (len(cands), a))
        return out

    def sols_dbg(self, sols):
        return [mp.nstr(l[0][1], 25) for c, l in sols]

    def in_field(self, val, maxcoeff=None):
        """is val in K_b = Q(b), b = b_1^0 (degree self.field_deg)?  PSLQ on (val, 1, b, ..., b^(d-1)) with coefficient bound ~ 10^(dps/(d+2)); True for val = 0"""
        if abs(val) < mpf(10) ** (-(mp.mp.dps // 2)):
            return True
        d = self.field_deg or 1
        b = self.x0[2]
        vec = [val] + ([b ** k for k in range(d)] if d > 1 else [mpf(1)])
        maxcoeff = maxcoeff or 10 ** max(8, mp.mp.dps // (d + 2) - 2)
        rel = mp.pslq(vec, tol=mpf(10) ** (-(mp.mp.dps - 20)), maxsteps=500000, maxcoeff=maxcoeff)
        return rel is not None and rel[0] != 0

    def series(self, N):
        """the unique formal solution to order N: returns (xs, log) with xs[n] = list of the 6 data at order n.  If several candidates survive all consistency tests (with the
        default tolerance 10^-(0.6 dps) this never happened in the final runs; with 1e-40 a near-solution survived at (28,13,1) of [14,7,0]), the physical one is selected by field
        membership: its coefficients lie in K_b = Q(b) (Galois stable: the chain is defined over Q) while spurious formal solutions are algebraic of higher degree.  Exactly one
        survivor must pass, else Degenerate is raised with the candidate data."""
        sols = self.solve_all(N)
        self.nsurvivors = len(sols)
        if len(sols) > 1 and self.field_deg is not None:
            keep = []
            for cols, log in sols:
                if all(self.in_field(mpf(l[1])) for l in log[:1]):
                    keep.append((cols, log))
            if len(keep) == 1:
                sols = keep
                self.selected_by = 'field'
            else:
                raise Degenerate('%d surviving solutions, %d of them with c_1 in K_b (c_1 of survivors: %s)' % (len(sols), len(keep), self.sols_dbg(sols)))
        if len(sols) != 1:
            raise Degenerate('%d surviving solutions' % len(sols))
        cols, log = sols[0]
        xs = [[cols[q][n] for q in range(6)] for n in range(N + 1)]
        return xs, log


def full_residual(lam, mu, phi, xs, th12=1):
    """recompute all eight equations from scratch on the series xs (list over orders of the six data); returns the list over orders of max |F_n|"""
    eng = Engine(lam, mu, phi, th12)
    return [max(abs(f) for f in eng.compute(n, x)) for n, x in enumerate(xs)], eng


# ============================================================================================ 5. validation against the exact reference series, runs at scale
PHIA = (Fr(1, 3), Fr(-2, 7), Fr(5, 11))
PHIB = (Fr(-2, 5), Fr(3, 8), Fr(7, 11))
TH12_REF = Fr(31, 35)


def load_ref(lam):
    """W5's exact series (ref/laurent-series-<lam>.txt, E4_SERIES): dict (mu, copy (1-based), phi tuple) -> [H, a, b] with H[s][n] mpf (theta12 = 31/35)"""
    from nu_general import _parse, W4DIR
    mp.mp.dps = DPS            # importing the earlier modules resets mp.dps
    fn = os.path.join(W4DIR, 'laurent-series-%d%d0.txt' % (lam[0], lam[1]))
    out = {}
    for line in open(fn, encoding='utf8'):
        if not line.startswith('E4_SERIES'):
            continue
        f = line.rstrip('\n').split('\t')
        mu = tuple(_parse(f[1]))
        cp = int(f[2])
        D0 = int(f[3])
        phi = tuple(Fr(x) for x in _parse(f[4]))
        D = D0 or 5
        sq = mp.sqrt(D)
        rows = []
        for fld in (f[6], f[7], f[8]):
            data = _parse(fld)
            rows.append([[mpq(Fr(p[0])) + sq * mpq(Fr(p[1])) for p in site] for site in data])
        out[(mu, cp, phi)] = rows
    return out


def to_theta_mp(xs, th12=TH12_REF):
    """series in units theta12 = 1 (list over orders of the six data) -> theta12: H, a ~ th12^(n-1), b ~ th12^(n-2) at order n"""
    t = mpq(Fr(th12))
    return [[x[0] * t ** (n - 1), x[1] * t ** (n - 1), x[2] * t ** (n - 2), x[3] * t ** (n - 1), x[4] * t ** (n - 1), x[5] * t ** (n - 2)] for n, x in enumerate(xs)]


def compare_series(ser_theta, ref_rows, N):
    """max abs and max relative (entries with |ref| > 1e-12) deviation of the six data of both sites, orders 0..N"""
    ea = mpf(0)
    er = mpf(0)
    nent = 0
    for n in range(N + 1):
        for s in range(2):
            for q in range(3):
                mine = ser_theta[n][3 * s + q]
                ref = ref_rows[q][s][n]
                d = abs(mine - ref)
                ea = max(ea, d)
                if abs(ref) > mpf(10) ** (-12):
                    er = max(er, d / abs(ref))
                nent += 1
    return ea, er, nent


def solve_sector_branches(lam, mu, phi, N, dps=None, verbose=False, th12=1):
    """all branches of the sector mu: returns list of dicts (copy, x0, rank, series xs (units th12=1), log, time, passes, residual per order)"""
    out = []
    if dps:
        setdps(dps)
    brs, gs = order0_branches(lam, mu, dps, return_poly=True)
    m_exp = mult_rs(lam, mu, weights_cached(lam))
    if len(brs) != m_exp:                      # never run a series on order-0 data with a wrong spectrum
        raise ValueError('pencil gives %d roots for the irrep %s, multiplicity (Racah-Speiser) is %d' % (len(brs), mu, m_exp))
    fdeg = field_degrees(gs, [x[2] for x in brs])
    for cp, x0 in enumerate(brs):
        t0 = time.perf_counter()
        S = BranchSolver(lam, mu, x0, phi, verbose=verbose, field_deg=fdeg[cp])
        rec = dict(mu=tuple(mu), copy=cp + 1, x0=x0, rank=S.rank, res0=S.res0)
        try:
            xs, log = S.series(N)
            rec.update(xs=xs, log=log, ok=True)
        except (Degenerate, NotImplementedError) as ex:
            rec.update(xs=None, log=None, ok=False, error=str(ex))
        rec['time'] = time.perf_counter() - t0
        rec['passes'] = S.passes
        rec['field_deg'] = fdeg[cp]
        rec['selected_by'] = S.selected_by
        rec['nsurv'] = getattr(S, 'nsurvivors', None)
        rec['S'] = S
        out.append(rec)
    return out


def validate(lam, N=10, phis=(PHIA, PHIB), dps=None, verbose=True):
    """acceptance 1: the numeric solver vs W5's exact series, all branches, both directions, orders 0..N"""
    lam = tuple(lam)
    if dps:
        setdps(dps)
    ref = load_ref(lam)
    mp.mp.dps = DPS
    rows = []
    for mu, m in decompose(lam):
        for phi in phis:
            for rec in solve_sector_branches(lam, mu, phi, N, dps):
                key = (rec['mu'], rec['copy'], tuple(phi))
                if not rec['ok']:
                    rows.append(dict(rec, ea=None, er=None, key=key))
                    if verbose:
                        print('  FAIL', key, rec['error'], flush=True)
                    continue
                ea, er, nent = compare_series(to_theta_mp(rec['xs']), ref[key], N)
                rows.append(dict(rec, ea=ea, er=er, key=key, nent=nent))
                if verbose:
                    print('  %s copy %d phi %s: rank %d, max abs dev %s, max rel dev %s, %.2fs' % (rec['mu'], rec['copy'], [str(x) for x in phi][0], rec['rank'], mp.nstr(ea, 3), mp.nstr(er, 3), rec['time']), flush=True)
    return rows


def run_family(lam, N, phis=(PHIA, PHIB), dps=None, verbose=True, cache=None, only=None, check_homog=True):
    """all branches of Lambda x Lambda (or the irreps in `only`) to order N for each direction; records per branch: rank, series (units theta12 = 1), residuals of all eight equations
    (recomputed from scratch), singular-value gap of J, time, passes; and, if check_homog, the homogeneity test x_n(2 phi) = 2^n x_n(phi) is NOT done here (see test_homogeneity).
    Saved to `cache` (pickle with strings) if given.  Returns the list of records."""
    lam = tuple(lam)
    recs = []
    for mu, m in decompose(lam):
        if only is not None and mu not in only:
            continue
        for phi in phis:
            for r in solve_sector_branches(lam, mu, phi, N, dps):
                r['mult'] = m
                r['phi'] = tuple(phi)
                S = r.pop('S')
                sv = S.svals()
                r['sv'] = sv
                r['kernel_dim'] = S.nk
                r['info'] = S.info
                r['ndeaths'] = S.ndeaths
                r['t_fix'] = dict(S.t_fix)
                if r['ok']:
                    res, _ = full_residual(lam, mu, phi, r['xs'])
                    r['resid'] = max(res)
                    r['resid_by_order'] = res
                else:
                    r['resid'] = None
                recs.append(r)
                if verbose:
                    print('%s mult %d copy %d phi %s rank %d  %s  resid %s  %.2fs' % (r['mu'], m, r['copy'], str(phi[0]), r['rank'], 'OK' if r['ok'] else 'FAIL: ' + r['error'],
                          mp.nstr(r['resid'], 3) if r['resid'] is not None else '-', r['time']), flush=True)
    if cache:
        dump = []
        for r in recs:
            d = {k: v for k, v in r.items() if k not in ('xs', 'x0', 'sv', 'resid', 'resid_by_order', 'log', 'res0')}
            d['x0'] = [mp.nstr(x, DPS) for x in r['x0']]
            d['xs'] = [[mp.nstr(x, DPS) for x in row] for row in r['xs']] if r['xs'] is not None else None
            d['sv'] = [mp.nstr(x, 20) for x in r['sv']]
            d['resid'] = mp.nstr(r['resid'], 5) if r['resid'] is not None else None
            d['res0'] = mp.nstr(r['res0'], 5)
            d['log'] = [(a, mp.nstr(c, 30), nc) for (a, c, nc) in r['log']] if r['log'] is not None else None
            dump.append(d)
        with open(cache, 'wb') as f:
            pickle.dump(dump, f)
    return recs


# ============================================================================================ 6. depth analysis: nu-series x exact kernel -> leading polynomial
def kernel_T(lam):
    """exact one-site-pair kernel Theta^{lam} in T-format {(i,j): {(n,p,q): Fraction}} (indices into pairs_of(lam)): directly from general_kernel_rep_<lam>.pkl if it exists, else
    by the block theorem Theta^{[l1,l2,0]}_ij = [l2!(l1-l2)!/l1!]^2 (a a')^(-l2) Theta^{[l1,0,0]}_ij from the one-row kernel (general_kernel_rep_<l1>00.pkl or master_kernel2_rep_<l1>00.pkl)"""
    import general_kernel as gk
    mp.mp.dps = DPS            # importing the earlier modules resets mp.dps
    lam = tuple(lam)
    l1, l2 = lam[0], lam[1]
    fn = os.path.join(HERE, 'general_kernel_rep_%d%d%d.pkl' % lam)
    if os.path.exists(fn):
        return gk.load_rep(lam)[0], 'general_kernel_rep_%d%d%d.pkl' % lam
    one = (l1, 0, 0)
    fn1 = os.path.join(HERE, 'general_kernel_rep_%d%d%d.pkl' % one)
    fn2 = os.path.join(HERE, 'master_kernel2_rep_%d%d%d.pkl' % one)
    if os.path.exists(fn1):
        T0 = gk.load_rep(one)[0]
        src = 'block theorem from general_kernel_rep_%d00.pkl' % l1
    elif os.path.exists(fn2):
        with open(fn2, 'rb') as f:
            T0 = pickle.load(f)
        src = 'block theorem from master_kernel2_rep_%d00.pkl' % l1
    else:
        fnK = os.path.join(HERE, 'kernel_Theta_%d%d%d.npz' % lam)
        if os.path.exists(fnK):                # agent K's exact kernels (Garner digits -> Fractions): [4,2,0], [8,4,0], [14,7,0]
            import kernel_L14 as KL
            T, arrs = KL.theta_T_from_file(fnK)
            mp.mp.dps = DPS
            return T, 'kernel_Theta_%d%d%d.npz (agent K, exact coefficients)' % lam
        raise FileNotFoundError('no one-row kernel for l1 = %d' % l1)
    p0 = pairs_of(one)
    idx0 = {pr: i for i, pr in enumerate(p0)}
    pr = pairs_of(lam)
    fac = Fr(factorial(l2) * factorial(l1 - l2), factorial(l1)) ** 2
    T = {}
    for i, pi in enumerate(pr):
        for j, pj in enumerate(pr):
            ent = T0.get((idx0[pi], idx0[pj]), {})
            T[(i, j)] = {(n, p - l2, q - l2): c * fac for (n, p, q), c in ent.items()}
    return T, src


def G_terms_from_T(lam, T, th12=TH12_REF):
    """{(i,j): [(alpha, beta, d, coef)]}: G_ij(eps w, eps w'; t) = sum coef eps^d w^alpha w'^beta, G = (w w')^(2 l1 - l2) Theta(1/w, 1/w', t), t = 1/theta12 (as nu_general.G_terms)"""
    l1, l2 = lam[0], lam[1]
    t = 1 / Fr(th12)
    out = {}
    for (i, j), ent in T.items():
        terms = {}
        for (n, p, q), c in ent.items():
            al, be = 2 * l1 - l2 - p, 2 * l1 - l2 - q
            assert al >= 0 and be >= 0, (lam, i, j, n, p, q)
            d = 4 * l1 - 2 * l2 - p - q
            key = (al, be, d)
            terms[key] = terms.get(key, 0) + Fr(c) * t ** n
        out[(i, j)] = [(al, be, d, c) for (al, be, d), c in terms.items() if c != 0]
    return out


def leading_analysis(lam, mu, phi, xs, Gt, tmax, thr=mpf(10) ** (-30), thr_depth=mpf(10) ** (-25)):
    """xs: series of the six data in units theta12 = 1 (list over orders).  Computes the nu-series at theta12 = 31/35 (replay on an engine with th12 = 31/35), the leading
    polynomials R_t(w, w') = sum_{d+n+m = t} G^[d]_ij nu_i^(n)(theta_1) nu_j^(m)(theta_2) for t <= tmax, and returns
      N_B measured (first t with max|R_t| > thr * (sum of |terms|)), the lower-order cancellation size, the depth (largest n or m of a contribution with |term| > thr_depth * max at t = N_B)."""
    N = len(xs) - 1
    th = to_theta_mp(xs)
    eng = Engine(lam, mu, phi, th12=TH12_REF)
    for n in range(N + 1):
        eng.compute(n, th[n])
    nus = eng.nu_series(N)
    prs = pairs_of(lam)
    npr = len(prs)
    Sij = {}
    for (i, j) in Gt:
        a, b = nus[0][i], nus[1][j]
        Sij[(i, j)] = [mp.fdot(a[:s + 1], b[s::-1]) for s in range(tmax + 1)]
    R = {}
    scale = {}
    for (i, j), terms in Gt.items():
        S = Sij[(i, j)]
        for (al, be, d, c0) in terms:
            c = mpq(c0)
            for s in range(0, tmax - d + 1):
                t = d + s
                v = c * S[s]
                Rt = R.setdefault(t, {})
                Rt[(al, be)] = Rt.get((al, be), mpf(0)) + v
    # scale of the contributions at each order: sum |c| |nu_i^(n)| |nu_j^(m)| over all terms (an upper estimate of the size of the individual contributions)
    absn = [[[abs(x) for x in nus[s][i]] for i in range(npr)] for s in (0, 1)]
    absS = {}
    for (i, j) in Gt:
        a, b = absn[0][i], absn[1][j]
        absS[(i, j)] = [mp.fdot(a[:s + 1], b[s::-1]) for s in range(tmax + 1)]
    for (i, j), terms in Gt.items():
        AS = absS[(i, j)]
        for (al, be, d, c0) in terms:
            c = abs(mpq(c0))
            for s in range(0, tmax - d + 1):
                scale[d + s] = scale.get(d + s, mpf(0)) + c * AS[s]
    rel = []
    for t in range(tmax + 1):
        mx = max((abs(v) for v in R.get(t, {}).values()), default=mpf(0))
        sc = scale.get(t, mpf(0))
        rel.append((t, mx, sc, mx / sc if sc > 0 else mpf(0)))
    NB = next((t for t, mx, sc, r in rel if r > thr), None)
    low = max([r for t, mx, sc, r in rel if NB is not None and t < NB], default=mpf(0))
    depth = None
    nmax = None
    if NB is not None:
        # contributions at t = NB
        contrib = {}
        for (i, j), terms in Gt.items():
            for (al, be, d, c0) in terms:
                if d > NB:
                    continue
                c = mpq(c0)
                for n in range(0, NB - d + 1):
                    m_ = NB - d - n
                    v = abs(c * nus[0][i][n] * nus[1][j][m_])
                    key = (i, n, j, m_)
                    contrib[key] = contrib.get(key, mpf(0)) + v
        mx = max(contrib.values())
        used = [k for k, v in contrib.items() if v > thr_depth * mx]
        depth = max(max(k[1], k[3]) for k in used)
        nmax = (max(k[1] for k in used), max(k[3] for k in used))
    return dict(NB=NB, low=low, depth=depth, nmax=nmax, rel=rel, mx=rel[NB][1] if NB is not None else None)


def kernel_summary(lam, T, th12=TH12_REF, pts=((Fr(3, 4), Fr(-12, 5)), (Fr(-5, 3), Fr(7, 2)))):
    """branch-independent summary of the kernel for the leading analysis: for every entry (i,j) and every eps-order d of G_ij(eps w, eps w'; t = 1/theta12) = sum_d eps^d P_{ij,d}(w, w'):
    A = sum |coef| over the monomials of that order and the values of P_{ij,d} at a few fixed points (w, w') (a polynomial that is non-zero is non-zero at a random point).
    Returns {(i,j): {d: (A, [P(w_k, w'_k)])}}.  Coefficients are combined over the powers of t first (as nu_general.G_terms)."""
    l1, l2 = lam[0], lam[1]
    tm = mpq(1 / Fr(th12))
    nmax = 0
    for ent in T.values():
        for (n, p, q) in ent:
            nmax = max(nmax, n)
    tpow = [mpf(1)]
    for n in range(nmax + 1):
        tpow.append(tpow[-1] * tm)
    amax = 2 * l1 - l2
    wt = [([mpq(w) ** a for a in range(amax + 1)], [mpq(w2) ** a for a in range(amax + 1)]) for (w, w2) in pts]
    KS = {}
    for (i, j), ent in T.items():
        acc = {}
        for (n, p, q), c in ent.items():
            al, be = 2 * l1 - l2 - p, 2 * l1 - l2 - q
            assert al >= 0 and be >= 0, (lam, i, j, n, p, q)
            key = (al, be, 4 * l1 - 2 * l2 - p - q)
            acc[key] = acc.get(key, mpf(0)) + mpq(c) * tpow[n]
        per_d = {}
        for (al, be, d), c in acc.items():
            if c == 0:
                continue
            ent_d = per_d.setdefault(d, [mpf(0)] + [mpf(0)] * len(pts))
            ent_d[0] += abs(c)
            for k, (wp, wq) in enumerate(wt):
                ent_d[1 + k] += c * wp[al] * wq[be]
        KS[(i, j)] = {d: (v[0], v[1:]) for d, v in per_d.items()}
    return KS


def leading_analysis_fast(lam, mu, phi, xs, KS, tmax, thr=mpf(10) ** (-30), thr_depth=mpf(10) ** (-25)):
    """same quantities as leading_analysis (N_B measured, lower-order cancellation, depth) from the branch-independent kernel summary; cost ~ (#entries) x tmax x (#d) per branch"""
    N = len(xs) - 1
    th = to_theta_mp(xs)
    eng = Engine(lam, mu, phi, th12=TH12_REF)
    for n in range(N + 1):
        eng.compute(n, th[n])
    nus = eng.nu_series(N)
    npts = len(next(iter(next(iter(KS.values())).values()))[1])
    R = [[mpf(0)] * npts for _ in range(tmax + 1)]
    scale = [mpf(0)] * (tmax + 1)
    absn = [[[abs(x) for x in nus[s][i]] for i in range(len(nus[s]))] for s in (0, 1)]
    for (i, j), per_d in KS.items():
        if not per_d:
            continue
        a, b = nus[0][i], nus[1][j]
        aa, bb = absn[0][i], absn[1][j]
        S = [mp.fdot(a[:s + 1], b[s::-1]) for s in range(tmax + 1)]
        SA = [mp.fdot(aa[:s + 1], bb[s::-1]) for s in range(tmax + 1)]
        for d, (A, vals) in per_d.items():
            if d > tmax:
                continue
            for s in range(0, tmax - d + 1):
                t = d + s
                scale[t] += A * SA[s]
                Ss = S[s]
                if Ss != 0:
                    for k in range(npts):
                        R[t][k] += vals[k] * Ss
    rel = []
    for t in range(tmax + 1):
        mx = max(abs(x) for x in R[t])
        rel.append((t, mx, scale[t], mx / scale[t] if scale[t] > 0 else mpf(0)))
    NB = next((t for t, mx, sc, r in rel if r > thr), None)
    low = max([r for t, mx, sc, r in rel if NB is not None and t < NB], default=mpf(0))
    depth = nmax = None
    if NB is not None:
        contrib = {}
        for (i, j), per_d in KS.items():
            for d, (A, vals) in per_d.items():
                if d > NB:
                    continue
                for n in range(0, NB - d + 1):
                    m_ = NB - d - n
                    v = A * absn[0][i][n] * absn[1][j][m_]
                    if v != 0:
                        contrib[(i, n, j, m_)] = v
        mxc = max(contrib.values())
        used = [k for k, v in contrib.items() if v > thr_depth * mxc]
        depth = max(max(k[1], k[3]) for k in used)
        nmax = (max(k[1] for k in used), max(k[3] for k in used))
    return dict(NB=NB, low=low, depth=depth, nmax=nmax, rel=rel)


def kB_rule(lam, mu):
    """k_B = 2 lam1 - (mu1 - mu3) (lam3 = 0);  N_B = 4 lam1 - 2 lam2 - k_B"""
    k = 2 * lam[0] - (mu[0] - mu[2])
    return k, 4 * lam[0] - 2 * lam[1] - k


# ============================================================================================ 7. checks of the order-0 route, timing study, CLI
def load_R_copies(lam):
    """agent R's ref/cg-slices-<lam>.txt: dict (mu, k) -> V_red eigenvalue (30 digits); a_1^0 = eigenvalue + eta/2"""
    from nu_general import _parse, W4DIR
    mp.mp.dps = DPS
    fn = os.path.join(W4DIR, 'cg-slices-%d%d%d.txt' % tuple(lam))
    out = {}
    for line in open(fn, encoding='utf8'):
        if line.startswith('COPY\t'):
            f = line.rstrip('\n').split('\t')
            out[(tuple(_parse(f[1])), int(f[2]))] = mpf(f[3])
    return out


def load_Q0_branches(lam):
    """agent Q0's ref/qsystem-spectrum-<lam>.txt (the CORRECTED pencil: GCD of all three consistency polynomials, equations themselves when a q_k has no unknowns):
    dict (mu, k) -> [V, H1, a1, b1, H2, a2, b2] (30 digits; theta12 = 1, sites at u = 1, 0; copies by increasing V = increasing a_1)"""
    from nu_general import _parse, W4DIR
    mp.mp.dps = DPS
    fn = os.path.join(W4DIR, 'qsystem-spectrum-%d%d%d.txt' % tuple(lam))
    out = {}
    num = lambda t: mpq(Fr(t)) if '/' in t else mpf(t.replace('*^', 'e'))
    for line in open(fn, encoding='utf8'):
        if line.startswith('BRANCH\t'):
            f = line.rstrip('\n').split('\t')
            out[(tuple(_parse(f[1])), int(f[2]))] = [num(t) for t in f[3:10]]
    return out


def check_order0(lam, dps=60, verbose=True):
    """order-0 route: pencil roots vs (a) W4's exact eta/alpha/beta (if untwisted-exact-<lam>.txt exists), (b) agent R's V_red eigenvalues (ref/cg-slices-<lam>.txt, 30 digits);
    also the residual of the eight equations at eps = 0 on x_0 (Engine) and the number of roots vs the Racah-Speiser multiplicity"""
    lam = tuple(lam)
    from nu_general import load_w4, W4DIR
    setdps(dps)
    dec = decompose(lam)
    out = dict(lam=lam, nirrep=len(dec), nbranch=sum(m for _, m in dec), count_mismatch=0, vs_w4=None, vs_R=None, vs_Q0=None, res0=mpf(0), time=0)
    t0 = time.time()
    w4 = None
    if os.path.exists(os.path.join(W4DIR, 'untwisted-exact-%d%d%d.txt' % lam)):
        w4 = load_w4(lam)
        out['vs_w4'] = mpf(0)
    cpR = None
    if os.path.exists(os.path.join(W4DIR, 'cg-slices-%d%d%d.txt' % lam)):
        cpR = load_R_copies(lam)
        out['vs_R'] = mpf(0)
    q0 = None
    if os.path.exists(os.path.join(W4DIR, 'qsystem-spectrum-%d%d%d.txt' % lam)):
        q0 = load_Q0_branches(lam)
        out['vs_Q0'] = mpf(0)
    setdps(dps)
    for mu, m in dec:
        br = order0_branches(lam, mu, dps)
        if len(br) != m:
            out['count_mismatch'] += 1
            continue
        eta = mpq(eta_of(lam, mu))
        for k, x0 in enumerate(br):
            eng = Engine(lam, mu, PHIA)
            F = eng.compute(0, x0)
            out['res0'] = max(out['res0'], max(abs(f) for f in F))
            if cpR is not None:
                out['vs_R'] = max(out['vs_R'], abs((x0[1] - eta / 2) - cpR[(mu, k + 1)]))
            if q0 is not None:
                out['vs_Q0'] = max(out['vs_Q0'], max(abs(x - y) for x, y in zip([x0[0], x0[1], x0[2], x0[3], x0[4], x0[5]], q0[(mu, k + 1)][1:7])), abs(x0[1] - eta / 2 - q0[(mu, k + 1)][0]))
            if w4 is not None:
                b = w4['sectors'][mu]['branches'][k]
                D = b['D0'] or 5
                pv = lambda p: mpq(Fr(p[0])) + mpq(Fr(p[1])) * mp.sqrt(D)
                ref = [pv(b['eta'][0]), pv(b['alpha'][0]), pv(b['beta'][0]), pv(b['eta'][1]), pv(b['alpha'][1]), pv(b['beta'][1])]
                out['vs_w4'] = max(out['vs_w4'], max(abs(x - y) for x, y in zip(x0, ref)))
    out['time'] = time.time() - t0
    if verbose:
        print('Lambda %s: %d irreps, %d branches, root-count mismatches %d, max eight-equation residual at eps=0 %s, vs W4 %s, vs agent R (30 digits) %s, vs agent Q0 (30 digits, all of V,H1,a1,b1,H2,a2,b2) %s, %.1fs' % (
            lam, out['nirrep'], out['nbranch'], out['count_mismatch'], mp.nstr(out['res0'], 3), mp.nstr(out['vs_w4'], 3) if out['vs_w4'] is not None else '-',
            mp.nstr(out['vs_R'], 3) if out['vs_R'] is not None else '-', mp.nstr(out['vs_Q0'], 3) if out['vs_Q0'] is not None else '-', out['time']), flush=True)
    return out


def test_homogeneity(lam, mu, N=8, t=Fr(2), phi=PHIA, dps=None):
    """x_n(t*phi) = t^n x_n(phi) for every copy of the sector (zeta = eps*phi: rescaling phi by t is the substitution eps -> t*eps); returns the max deviation"""
    dps = dps or DPS
    setdps(dps)
    worst = mpf(0)
    for x0 in order0_branches(lam, mu, dps):
        a = BranchSolver(lam, mu, x0, phi).series(N)[0]
        b = BranchSolver(lam, mu, x0, tuple(t * p for p in phi)).series(N)[0]
        for n in range(N + 1):
            for q in range(6):
                worst = max(worst, abs(b[n][q] - mpq(t) ** n * a[n][q]))
    return worst


def timing_study(lam, N=40, which=('top', 'oct', 'rank6'), dps=None, verbose=True):
    """cumulative time to fix order n (BranchSolver.t_fix) on representative branches: the top state (k_B = 0), the rank-5 branch with the largest multiplicity (first copy), a rank-6 branch"""
    lam = tuple(lam)
    dps = dps or DPS
    setdps(dps)
    dec = decompose(lam)
    cand = {}
    cand['top'] = (2 * lam[0], 2 * lam[1], 0)
    mx = max(dec, key=lambda t: (t[1], -t[0][2]))
    cand['oct'] = mx[0]
    cand['rank6'] = dec[-1][0]                    # the singlet-like irrep
    out = {}
    for key in which:
        mu = cand[key]
        t0 = time.perf_counter()
        br = order0_branches(lam, mu, dps)
        t_p = time.perf_counter() - t0
        S = BranchSolver(lam, mu, br[0], PHIA)
        t1 = time.perf_counter()
        xs, log = S.series(N)
        tt = time.perf_counter() - t1
        out[key] = dict(mu=mu, rank=S.rank, t_pencil=t_p, t_total=tt, t_fix=dict(S.t_fix), info=dict(S.info), passes=S.passes)
        if verbose:
            tf = S.t_fix
            sel = [n for n in (5, 10, 15, 20, 30, 40, 50) if n in tf and n <= N]
            print('  %s %s rank %d: pencil %.1fs; total to N=%d: %.1fs (%d passes); cumulative time to fix order n: %s' % (key, mu, S.rank, t_p, N, tt, S.passes,
                  ', '.join('n=%d: %.2fs' % (n, tf[n]) for n in sel)), flush=True)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--stage', required=True, help='order0 | validate | family | depth | timing')
    ap.add_argument('--lam', type=int, nargs=2, default=[2, 1])
    ap.add_argument('--N', type=int, default=10)
    ap.add_argument('--dps', type=int, default=80)
    ap.add_argument('--phis', default='AB')
    ap.add_argument('--cache', default=None)
    a = ap.parse_args()
    lam = (a.lam[0], a.lam[1], 0)
    setdps(a.dps)
    phis = tuple({'A': PHIA, 'B': PHIB}[c] for c in a.phis)
    if a.stage == 'order0':
        check_order0(lam, dps=min(a.dps, 60))
    elif a.stage == 'validate':
        for l in ((2, 1, 0), (3, 1, 0), (3, 2, 0)):
            t0 = time.time()
            rows = validate(l, a.N, phis, a.dps, verbose=False)
            nf = sum(1 for r in rows if not r['ok'])
            ea = max(r['ea'] for r in rows if r['ok'])
            er = max(r['er'] for r in rows if r['ok'])
            print('Lambda %s: %d (branch, direction) runs, %d failures, max abs dev %s, max rel dev %s, %.1fs' % (l, len(rows), nf, mp.nstr(ea, 3), mp.nstr(er, 3), time.time() - t0), flush=True)
    elif a.stage == 'family':
        t0 = time.time()
        recs = run_family(lam, a.N, phis, a.dps, cache=a.cache or os.path.join(HERE, 'bethe_series_scale_%d%d0.pkl' % lam[:2]))
        nf = sum(1 for r in recs if not r['ok'])
        print('TOTAL %s N=%d dps=%d: %.1fs, runs %d, failures %d, max resid %s' % (lam, a.N, a.dps, time.time() - t0, len(recs), nf,
              mp.nstr(max(r['resid'] for r in recs if r['ok']), 3)), flush=True)
    elif a.stage == 'depth':
        T, src = kernel_T(lam)
        setdps(a.dps)
        KS = kernel_summary(lam, T)
        del T
        print('Lambda %s, kernel: %s' % (lam, src), flush=True)
        for mu, m in decompose(lam):
            kB, NBp = kB_rule(lam, mu)
            for r in solve_sector_branches(lam, mu, phis[0], NBp + 2, a.dps):
                an = leading_analysis_fast(lam, mu, phis[0], r['xs'], KS, NBp + 2)
                print('  %s mult %d copy %d: k_B=%d N_B(rule)=%d N_B measured=%s lower-order residue %s depth=%s' % (mu, m, r['copy'], kB, NBp, an['NB'], mp.nstr(an['low'], 3), an['depth']), flush=True)
    elif a.stage == 'timing':
        timing_study(lam, a.N, dps=a.dps)


if __name__ == '__main__':
    main()
