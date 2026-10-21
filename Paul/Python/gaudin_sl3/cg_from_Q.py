"""F3 (round 5): CG coefficients of adj (x) adj = [2,1,0]^2 from Q-side data only.   python -X utf8 cg_from_Q.py   (about a minute)

WHAT IT DOES (equation numbers: Gaudin-sl3-CG-notes.md).
 1. Bethe-side series.  For each of the six HW-limit branches the exact eps-series of the Laurent data x = (H_s,a_s,b_s)_{s=1,2} (zeta = eps*phi, theta_1-theta_2 = 1) is solved
    from the six no-log conditions + the two global rules (eq. 21 and sect. 2.3), via bethe_series.py (degenerate Fredholm solver; NO chain, NO fit; the eps^0 point is
    itself solved from the eps = 0 equations, see bethe_series.branches_w4/order0_points).  Order 8.
 2. nu-formulas (sect. 7, eq. 70): nu(theta_s) = (nu_02,nu_03,nu_12,nu_13,nu_23) = (1, pi, kappa, kappa*pi, eps_s) are polynomials in the Laurent data of both sites
    (nu_general.nu_symbolic) -> exact eps-series of nu_i(theta_1), nu_j(theta_2) in Q(sqrt5).
 3. Kernel (sect. 8): G_ij(eps w1, eps w1'; t) = sum c t^n eps^d w1^al w1'^be from the committed kernel general_kernel_rep_210.pkl (exact), t = 1/theta12 = 1.
 4. Leading polynomial (sect. 9.3-9.4): L_B(w1,w1') = sum_{d+n+m = N_B} G^[d]_ij nu_i^(n)(theta_1) nu_j^(m)(theta_2); all lower total orders cancel EXACTLY (checked) and
    N_B = 6 - k_B.   Identity (70):  L_B = chat * <phi_u (x) e1|Psi_0> / X_B  (chat, Gauss parameters of phi_u: untwist_limit.gf_polys, exact).
 5. Dictionary (sect. 9.5): solve L_B = chat * sum_nu s_nu [e1^T rho(n')]_nu coefficientwise in (w1,w1') (exact linear algebra over Q(sqrt5), overdetermined, must be
    consistent) -> s_nu = C_nu / X_B.
 6. Absolute normalisation (Q1, RESULTS-gaudin-norm.md):  N(Psi) = <Psi|Psi>_G/<x0^G|Psi>^2 = (9/4)|det J6| |z12|^(mu3-6)|z13|^(mu2-6)|z23|^(mu1-6) (J6 = Jacobian of the six
    no-log conditions of untwisted_Q.frobenius_conditions w.r.t. (H1,a1,b1,H2,a2,b2); z_ij = eps(phi_i-phi_j)), and lim eps^(4+2k_B) N = 1/X_B^2, so |X_B| = lim^(-1/2)
    and C_nu = X_B s_nu with ||Psi_0||_G = 1.  det J6 needs the Laurent data to eps^(8-2k_B) (series of step 1, exact).  The overall sign of a slice is a phase
    convention (not fixed by the Q side): it is aligned to the exact table's first nonzero entry.
 Output: the six slices C_nu = <e_nu (x) e_1|Psi_0> vs W4's exact table (digits of agreement), k_B, N_B, X_B vs W4's closed forms.
"""
import sys, os, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sympy as sp
import mpmath as mp
import bethe_series as B
from bethe_series import QS, Solver, solve_series, qs_from_sympy, branches_w4, PHIA
import nu_general as NG
import untwist_limit as U
import general_kernel as gk
import untwisted_Q as UQ

mp.mp.dps = 90
LAM = (2, 1, 0)
LAM2 = (2, 1)
PHI = PHIA
M = 8                                  # series depth (orders 0..M)

# ---------------------------------------------------------------------------------- truncated series with QS coefficients
def Z(D): return QS(0, 0, D)
def s_const(c, D): return [QS(c, 0, D) if not isinstance(c, QS) else c] + [Z(D)] * M
def s_add(a, b): return [x + y for x, y in zip(a, b)]
def s_scale(a, c): return [x * c for x in a]
def s_mul(a, b, D):
    out = [Z(D) for _ in range(M + 1)]
    for i, x in enumerate(a):
        if x.iszero(): continue
        for j in range(M + 1 - i):
            y = b[j]
            if not y.iszero(): out[i + j] = out[i + j] + x * y
    return out
def s_shift(a, k, D):
    return ([Z(D)] * k + a)[:M + 1]

class SeriesEval:
    """evaluate polynomials (sympy) in named variables whose values are series"""
    def __init__(self, D):
        self.D = D; self.pc = {}
    def power(self, key, ser, k):
        L = self.pc.setdefault(key, [s_const(1, self.D)])
        while len(L) <= k:
            L.append(s_mul(L[-1], ser, self.D))
        return L[k]
    def eval(self, expr, vals, shifts=None):
        """expr: sympy polynomial in the symbols of vals (dict symbol -> (key, series)); returns a series"""
        syms = list(vals)
        P = sp.Poly(expr, *syms)
        out = [Z(self.D)] * (M + 1)
        for mon, co in P.terms():
            t = s_const(QS(Fr(int(co.p), int(co.q)), 0, self.D), self.D)
            for s, k in zip(syms, mon):
                if k:
                    key, ser = vals[s]
                    t = s_mul(t, self.power(key, ser, k), self.D)
            out = s_add(out, t)
        return out

def consts_series(phi, D):
    p = [QS(x, 0, D) for x in phi]
    s1 = p[0] + p[1] + p[2]; e2 = p[0] * p[1] + p[0] * p[2] + p[1] * p[2]; e3 = p[0] * p[1] * p[2]
    sig = s_shift(s_const(s1, D), 1, D); E2 = s_shift(s_const(e2, D), 2, D); E3 = s_shift(s_const(e3, D), 3, D)
    return sig, E2, E3

# ---------------------------------------------------------------------------------- the nu series (step 2)
def nu_series_of(xs, phi, D):
    """xs[n] = six data at order n (QS).  returns nus[site][i] = series of nu_i(theta_site), i over nu_general.pairs_of(2,1)"""
    prs, nus, nolog = NG.nu_symbolic(2, 1)
    Hs, As, Bs, Ht, At, Bt, SG1, E2s, E3s, DL = NG.SYMS
    sig, E2, E3 = consts_series(phi, D)
    X = [[x[q] for x in xs] + [Z(D)] * (M + 1 - len(xs)) for q in range(6)]
    SE = SeriesEval(D)
    out = []
    for s, (t, dl) in enumerate(((1, 1), (0, -1))):
        vals = {Hs: (('H', s), X[3 * s]), As: (('a', s), X[3 * s + 1]), Bs: (('b', s), X[3 * s + 2]),
                Ht: (('H', t), X[3 * t]), At: (('a', t), X[3 * t + 1]), Bt: (('b', t), X[3 * t + 2]),
                SG1: (('s1', 0), sig), E2s: (('e2', 0), E2), E3s: (('e3', 0), E3)}
        out.append([SE.eval(sp.expand(e.subs(DL, dl)), vals) for e in nus])
    return prs, out

# ---------------------------------------------------------------------------------- kernel pieces (step 3)
_GT = {}
def G_units():
    """{(i,j): {(al,be,d): Fraction}}: G_ij = sum c eps^d w1^al w1'^be at t = 1/theta12 = 1; G = (w w')^{2 l1 - l2} Theta (nu_general.G_terms with t = 1)"""
    if 'g' in _GT: return _GT['g']
    Th, rep = gk.load_rep(LAM)
    l1, l2 = LAM[0], LAM[1]
    out = {}
    for (i, j), ent in Th.items():
        terms = {}
        for (n, p, q), c in ent.items():
            al, be = 2 * l1 - l2 - p, 2 * l1 - l2 - q
            d = 4 * l1 - 2 * l2 - p - q
            terms[(al, be, d)] = terms.get((al, be, d), 0) + Fr(c)
        out[(i, j)] = {k: v for k, v in terms.items() if v != 0}
    _GT['g'] = out
    return out

def leading_poly(nus, D, tmax):
    """R[t] = {(a,b): QS} for total order t = d+n+m <= tmax"""
    G = G_units()
    R = {}
    for (i, j), terms in G.items():
        for (al, be, d), c in terms.items():
            for n in range(0, tmax - d + 1):
                x1 = nus[0][i][n]
                if x1.iszero(): continue
                for m in range(0, tmax - d - n + 1):
                    y = nus[1][j][m]
                    if y.iszero(): continue
                    t = d + n + m
                    Rt = R.setdefault(t, {})
                    Rt[(al, be)] = Rt.get((al, be), Z(D)) + x1 * y * c
    return R

# ---------------------------------------------------------------------------------- dictionary (step 5)
def dictionary(LB, phi, D):
    """solve LB = sum_nu s_nu * (chat [e1 rho]_nu) exactly; returns (s as list over nu of QS or None, number of equations, residual-consistent flag)"""
    R = NG.rep_of(LAM)
    polys, chat = U.gf_polys(R.ch, NG.PHI if tuple(phi) == tuple(NG.PHI) else tuple(phi))
    cols = {nu: B_poly(P, D) for nu, P in enumerate(polys) if P is not None}
    mons = sorted(set(m for c in cols.values() for m in c) | set(LB))
    nus = sorted(cols)
    A = [[cols[nu].get(m, Z(D)) for nu in nus] + [LB.get(m, Z(D))] for m in mons]
    T, Rr, piv = B.rref(A, D)
    k = len(nus)
    consistent = (k not in piv)
    s = [None] * R.d
    if consistent:
        assert piv == list(range(k)), 'dictionary columns dependent'
        for i, nu in enumerate(nus):
            s[nu] = Rr[i][k]
    return s, len(mons), consistent

def B_poly(P, D):
    return {m: QS(Fr(int(sp.Rational(c).p), int(sp.Rational(c).q)), 0, D) for m, c in zip(P.monoms(), P.coeffs())}

# ---------------------------------------------------------------------------------- det J6 (step 6)
_J6 = {}
def j6_exprs():
    if 'j' in _J6: return _J6['j']
    vars_, C = UQ.frobenius_conditions()
    Dv = sp.symbols('H1 a1 b1 H2 a2 b2')
    e1, e2, e3 = sp.symbols('e1 e2 e3')
    Fs = []
    for s in range(2):
        o = 1 - s
        ds = 1 if s == 0 else -1                 # d = theta_1 - theta_2 = 1
        reg = lambda k, m: (-1) ** k * sp.binomial(m + k - 1, k) * sp.Integer(ds) ** (-m - k)
        H, a, b = Dv[3 * s:3 * s + 3]; Ho, ao, bo = Dv[3 * o:3 * o + 3]
        sk = [3 * reg(k, 1) + (e1 if k == 0 else 0) for k in range(3)]
        hk = [Ho * reg(k, 1) + 3 * reg(k, 2) + (e2 if k == 0 else 0) for k in range(3)]
        tk = [bo * reg(k, 1) + ao * reg(k, 2) + (e3 if k == 0 else 0) for k in range(3)]
        sub = dict(zip(vars_, [H, a, b] + sk + hk + tk))
        for (n, idx, poly) in C:
            Fs.append(sp.expand(poly.subs(sub)))
    Fm = sp.Matrix(Fs)
    J = Fm.jacobian(sp.Matrix(Dv))
    _J6['j'] = (Fs, J, Dv, (e1, e2, e3))
    return _J6['j']

def det_series(Mx, D):
    n = len(Mx)
    memo = {}
    def minor(r, cols):
        if r == n: return s_const(1, D)
        key = (r, cols)
        if key in memo: return memo[key]
        tot = [Z(D)] * (M + 1)
        for pos, c in enumerate(cols):
            if all(x.iszero() for x in Mx[r][c]): continue
            sub = minor(r + 1, tuple(x for x in cols if x != c))
            term = s_mul(Mx[r][c], sub, D)
            tot = s_add(tot, term) if pos % 2 == 0 else s_add(tot, s_scale(term, -1))
        memo[key] = tot
        return tot
    return minor(0, tuple(range(n)))

def j6_series(xs, phi, D):
    Fs, J, Dv, (e1, e2, e3) = j6_exprs()
    sig, E2, E3 = consts_series(phi, D)
    X = [[x[q] for x in xs] + [Z(D)] * (M + 1 - len(xs)) for q in range(6)]
    vals = {Dv[q]: (('x', q), X[q]) for q in range(6)}
    vals[e1] = (('e1', 0), sig); vals[e2] = (('e2', 0), E2); vals[e3] = (('e3', 0), E3)
    SE = SeriesEval(D)
    Fres = [SE.eval(f, vals) for f in Fs]
    Jm = [[SE.eval(J[i, j], vals) for j in range(6)] for i in range(6)]
    return det_series(Jm, D), Fres

def qs_mp(q):
    return mp.mpf(q.a.numerator) / q.a.denominator + mp.sqrt(q.D) * mp.mpf(q.b.numerator) / q.b.denominator

def xb_closed(mu, cp, p, D):
    """|X_B| from W4's closed forms (notes eq. 73, theta12 = 1) for Psi_0 of G-norm 1 (the (4,1,1),(3,3,0),(3,2,1),(2,2,2) formulas are for the unnormalised kets: divide by sqrt(g0) = sqrt(2), sqrt(2), sqrt(1080(7-3 sqrt5)), sqrt(288))"""
    f = [mp.mpf(x.numerator) / x.denominator for x in p]
    d12, d13, d23 = f[0] - f[1], f[0] - f[2], f[1] - f[2]
    if mu == (4, 2, 0): return abs(d12 ** 2)
    if mu == (4, 1, 1): return abs(d12 * d13 * d23 / 2) / mp.sqrt(2)
    if mu == (3, 3, 0): return abs(d12 ** 3 / 2) / mp.sqrt(2)
    if mu == (2, 2, 2): return abs((d12 * d13 * d23) ** 2 / 16) / mp.sqrt(288)
    if mu == (3, 2, 1):
        sg = 1 if cp == 0 else -1
        s5 = mp.sqrt(5)
        return abs(5766 * (sg * s5 - 3) / 6125 * d12 ** 2 * d13 * d23) / mp.sqrt(1080 * (7 - sg * 3 * s5)) / (mp.mpf(31) / 35) ** 2   # (73) is written at theta12 = 31/35; X_B ~ theta12^2
    return None

# ---------------------------------------------------------------------------------- main
def run(verbose=True):
    R = NG.rep_of(LAM)
    br = branches_w4(LAM2)
    table = []
    t00 = time.time()
    for (mu, cp), (D, x0) in br.items():
        t0 = time.time()
        S = Solver(LAM2, mu, [qs_from_sympy(v, D) for v in x0], list(PHI), D)
        xs, log = solve_series(S, M)                      # step 1
        prs, nus = nu_series_of(xs, PHI, D)               # step 2
        Rl = leading_poly(nus, D, M)                      # steps 3-4
        nonzero = sorted(t for t, P in Rl.items() if any(not c.iszero() for c in P.values()))
        NB = nonzero[0]
        kB = 6 - NB
        LB = {m: c for m, c in Rl[NB].items() if not c.iszero()}
        s, neq, ok = dictionary(LB, PHI, D)               # step 5
        # step 6: norm from det J6
        det, Fres = j6_series(xs, PHI, D)
        Fzero = all(all(c.iszero() for c in f) for f in Fres)
        lead = next(n for n, c in enumerate(det) if not c.iszero())
        p = PHI
        d12 = abs(mp.mpf(p[0].numerator) / p[0].denominator - mp.mpf(p[1].numerator) / p[1].denominator)
        d13 = abs(mp.mpf(p[0].numerator) / p[0].denominator - mp.mpf(p[2].numerator) / p[2].denominator)
        d23 = abs(mp.mpf(p[1].numerator) / p[1].denominator - mp.mpf(p[2].numerator) / p[2].denominator)
        lim = mp.mpf(9) / 4 * abs(qs_mp(det[lead])) * d12 ** (mu[2] - 6) * d13 ** (mu[1] - 6) * d23 ** (mu[0] - 6)
        XB = 1 / mp.sqrt(lim)
        XBc = xb_closed(mu, cp, p, D)
        C = [None if q is None else XB * qs_mp(q) for q in s]
        ex = NG.w4_cg(R, mu, cp + 1)                       # exact W4 slice in the same (gl3_reps) one-site index
        # sign alignment to the exact table
        sg = 1
        for a_, b_ in zip(C, ex):
            if a_ is not None and abs(b_) > mp.mpf(10) ** -30:
                sg = 1 if a_ * b_ > 0 else -1; break
        C = [None if c is None else sg * c for c in C]
        dev = max(abs((c if c is not None else 0) - e) for c, e in zip(C, ex))
        rec = dict(XBclosed=XBc, mu=mu, cp=cp, kB=kB, NB=NB, detlead=lead, expected_detlead=8 - 2 * kB, consistent=ok, neq=neq, Fzero=Fzero, XB=XB, C=C, ex=ex, dev=dev, time=time.time() - t0)
        table.append(rec)
        if verbose:
            nzs = [(nu, mp.nstr(c, 14)) for nu, c in enumerate(C) if c is not None and abs(c) > mp.mpf(10) ** -30]
            print('%s copy %d: k_B = %d (N_B = %d, lower total orders vanish exactly), dictionary consistent: %s (%d coefficient equations), no-log conditions vanish along the series: %s' % (mu, cp, kB, NB, ok, neq, Fzero))
            print('     |X_B| from det J6 = %s vs W4 closed form (73) %s' % (mp.nstr(XB, 14), None if XBc is None else mp.nstr(XBc, 14)))
            print('     det J6 ~ eps^%d (expected %d), 1/X_B^2 = %s;  slice C_nu (gl3_reps index): %s' % (lead, 8 - 2 * kB, mp.nstr(lim, 12), nzs))
            print('     exact (W4): %s   max |Q - exact| = %s   [%.1fs]' % ([(nu, mp.nstr(e, 14)) for nu, e in enumerate(ex) if abs(e) > mp.mpf(10) ** -30], mp.nstr(dev, 3), rec['time']), flush=True)
    print('total %.0fs' % (time.time() - t00))
    return table

if __name__ == '__main__':
    if 'phi2' in sys.argv:
        PHI = B.PHIB                       # second direction phi' = (-2/5, 3/8, 7/11): the result is phi-independent
    print('direction phi =', [str(x) for x in PHI])
    run()
