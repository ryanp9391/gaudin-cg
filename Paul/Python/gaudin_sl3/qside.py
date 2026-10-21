"""Q side of the sl(3) Gaudin chain (L=2, Lambda=(2,1,0)): Q1 opers, Q2 Wronskian/QQ, Q3 brackets/det N.
Run top to bottom:  python -X utf8 qside.py     (prints the tables; caches in qside_cache.pkl)
"""
import sys, os, time, itertools, pickle
from fractions import Fraction as Fr
from math import comb, factorial
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
import sympy
from gaudin_chain import *

DPS = 60
mp.mp.dps = DPS
ZETA = [Fr(3, 2), Fr(5, 4), Fr(2, 3)]
OM = [Fr(5, 2), Fr(7, 3), Fr(3, 5)]
OMP = [Fr(7, 4), Fr(-2, 5), Fr(9, 7)]
TOL = mp.mpf(10) ** (-35)

def mpf(x):
    if isinstance(x, Fr):
        return mp.mpf(x.numerator) / x.denominator
    return mp.mpf(x)

# ------------------------------------------------------------------ rational-function eigenvalue functions
class PF:
    """c0 + sum_{(s,m)} c_{s,m} (u-theta_s)^(-m)"""
    def __init__(self, c0, terms, theta):
        self.c0, self.terms, self.theta = c0, terms, theta
    def __call__(self, u, k=0):
        tot = self.c0 if k == 0 else 0
        for (s, m), c in self.terms.items():
            if c == 0:
                continue
            f = 1
            for j in range(k):
                f *= -(m + j)
            tot = tot + c * f * (u - self.theta[s]) ** (-m - k)
        return tot

def pf_S2(lam2, theta):
    return PF(lam2[0], {(0, 1): lam2[1], (0, 2): lam2[2], (1, 1): lam2[3], (1, 2): lam2[4]}, theta)

def pf_S3(lam3, theta):
    return PF(lam3[0], {(0, 1): lam3[1], (0, 2): lam3[2], (0, 3): lam3[3],
                        (1, 1): lam3[4], (1, 2): lam3[5], (1, 3): lam3[6]}, theta)

def pf_S1(trg, lam1, theta):
    return PF(trg, {(0, 1): lam1, (1, 1): lam1}, theta)

# ------------------------------------------------------------------ ODE polynomial solver
def dpoly(j, kappa, n, u):
    """(d + kappa)^j u^n  (number-type generic; u may be Fraction or mpf)"""
    tot = 0
    for i in range(0, min(j, n) + 1):
        ff = 1
        for t in range(i):
            ff *= (n - t)
        tot = tot + comb(j, i) * kappa ** (j - i) * ff * u ** (n - i) if (n - i) >= 0 else tot
    return tot

def ser_mul(A, B, n=4):
    return [sum(A[i] * B[k - i] for i in range(k + 1)) for k in range(n)]

def ser_pow(x, e, n=4):
    """(x+delta)^e as Taylor series in delta, e integer (may be negative)"""
    out = []
    c = 1
    for k in range(n):
        out.append(c * x ** (e - k))
        c = c * Fr(e - k, k + 1) if isinstance(x, Fr) else c * (e - k) / (k + 1)
    return out

def ser_exp(kappa, n=4):
    return [kappa ** k / factorial(k) for k in range(n)]

def ser_pf(S, u, n=4):
    return [S(u, k) / factorial(k) for k in range(n)]

def ode_row(which, S1, S2, S3, kappa, d, u, gauge=0):
    """Values at u of the ODE applied to e^{kappa u} u^n (n = 0..d); common factor e^{kappa u} dropped.
    'D': psi''' - S1 psi'' + S2 psi' - S3 psi        (psi = Q)
    'A': -psi''' - (S1 psi)'' - (S2 psi)' - S3 psi   with psi = Pi^{-gauge} Q, Pi = prod_s (u - theta_s);
         the row is multiplied by Pi^{gauge} (so gauge = 0 is the plain formal adjoint)."""
    th = S1.theta
    row = []
    s1, s2, s3 = ser_pf(S1, u), ser_pf(S2, u), ser_pf(S3, u)
    w = 1
    for t in th:
        w = w * (u - t) ** gauge
    for n in range(d + 1):
        psi = ser_mul(ser_pow(u, n), ser_exp(kappa))
        if which == 'A':
            for t in th:
                psi = ser_mul(psi, ser_pow(u - t, -gauge))
        dv = [psi[k] * factorial(k) for k in range(4)]
        if which == 'D':
            v = dv[3] - s1[0] * dv[2] + s2[0] * dv[1] - s3[0] * dv[0]
        else:
            P1 = ser_mul(s1, psi)
            P2 = ser_mul(s2, psi)
            v = -dv[3] - 2 * P1[2] - P2[1] - s3[0] * psi[0]
            v = v * w
        row.append(v)
    return row

def sample_points(theta, npts, exact):
    pts = []
    for k in range(npts):
        x = Fr(-17, 10) + Fr(4 * k, max(npts - 1, 1)) + Fr(1, 97 + k)
        if any(x == t for t in theta):
            x += Fr(1, 1000)
        pts.append(x)
    return pts

def solve_poly_exact(which, S1, S2, S3, kappa, theta, dmax=8, gauge=0):
    """exact (Fraction) version: smallest degree d with a nonzero polynomial solution; returns (d, nullity, monic coeffs)"""
    for d in range(dmax + 1):
        pts = sample_points(theta, d + 1 + 6, True)
        rows = []
        for u in pts:
            w = 1
            for t in theta:
                w = w * (u - t) ** 3
            rows.append([w * z for z in ode_row(which, S1, S2, S3, kappa, d, u, gauge)])
        M = sympy.Matrix([[sympy.Rational(z.numerator, z.denominator) for z in r] for r in rows])
        ns = M.nullspace()
        if ns:
            v = ns[0]
            top = max(i for i in range(d + 1) if v[i] != 0)
            v = [v[i] / v[top] for i in range(top + 1)]
            return top, len(ns), [Fr(int(c.p), int(c.q)) for c in v]
    return None, 0, None

def solve_poly_mp(which, S1, S2, S3, kappa, theta, dmax=8, gauge=0):
    """mp version: smallest degree with a polynomial solution of the ODE (absolute singular-value test).
    returns (degree, nullity at that degree, sigma_min, monic coefficient list)"""
    th = [mpf(t) for t in theta]
    for d in range(dmax + 1):
        npts = d + 1 + 6
        pts = [mpf(x) for x in sample_points(theta, npts, False)]
        rows = []
        for u in pts:
            w = 1
            for t in th:
                w = w * (u - t) ** 3
            r = ode_row(which, S1, S2, S3, kappa, d, u, gauge)
            rows.append([w * z for z in r])
        scl = [max(abs(pts[k] ** n * (1 if True else 0)) * abs(rows[k][n] / (pts[k] ** n) if pts[k] != 0 else 0) * 0 + abs(pts[k]) ** n for k in range(npts)) for n in range(d + 1)]
        M = mp.matrix(npts, d + 1)
        for k in range(npts):
            for n in range(d + 1):
                M[k, n] = rows[k][n] / scl[n]
        U, Sg, V = mp.svd_r(M)
        nul = sum(1 for i in range(len(Sg)) if Sg[i] < TOL)
        if nul > 0:
            v = [V[V.rows - 1, i] / scl[i] for i in range(d + 1)]
            mx = max(abs(x) for x in v)
            top = d        # smallest d with a solution; the leading coefficient may be tiny (e.g. eps-suppressed in the untwisting regime)
            v = [x / v[top] for x in v[:top + 1]]
            return top, nul, Sg[len(Sg) - 1], v
    return None, 0, None, None

# ------------------------------------------------------------------ polynomial helpers (mp lists, index = power)
def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]

def pscale(a, c):
    return [c * x for x in a]

def pmul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out

def pder(a):
    return [i * a[i] for i in range(1, len(a))] or [0]

def pdivmod(a, b):
    a = list(a)
    while len(b) > 1 and b[-1] == 0:
        b = b[:-1]
    q = [0] * max(len(a) - len(b) + 1, 1)
    for i in range(len(a) - len(b), -1, -1):
        c = a[i + len(b) - 1] / b[-1]
        q[i] = c
        for j in range(len(b)):
            a[i + j] -= c * b[j]
    return q, a[:len(b) - 1] if len(b) > 1 else [0]

def peval(a, x):
    r = 0
    for c in reversed(a):
        r = r * x + c
    return r

def trim(a, tol=None):
    tol = tol or mp.mpf(10) ** (-30)
    mx = max(abs(x) for x in a) or 1
    a = list(a)
    while len(a) > 1 and abs(a[-1]) < tol * mx:
        a.pop()
    return a

class Qfun:
    """e^{kappa u} q(u)"""
    def __init__(self, kappa, q):
        self.kappa, self.q = kappa, q
    def __call__(self, u, k=0):
        # k-th derivative: e^{kappa u} (d + kappa)^k q
        q = self.q
        cur = q
        for _ in range(k):
            cur = padd(pder(cur), pscale(cur, self.kappa))
        return mp.exp(self.kappa * u) * peval(cur, u)
    def poly_der(self, k=1):
        cur = self.q
        for _ in range(k):
            cur = padd(pder(cur), pscale(cur, self.kappa))
        return cur


# ================================================================== setup / Q1
THETA_MP = None
GAUGE = 2          # direct Q_i = Pi^GAUGE * psi_i, psi_i solving the formal adjoint equation (found in Q1)

def build_states():
    ch = Chain()
    g = [[ZETA[i] if i == j else Fr(0) for j in range(3)] for i in range(3)]
    ops = Operators(ch, g)
    C = ops.fit()
    st = eigenstates(ch, ops, C, dps=DPS)
    res = eigenfunctions(ch, C, st)
    return ch, ops, C, st, res

def Qfun_from(kappa, coeffs):
    return Qfun(kappa, [mpf(c) for c in coeffs])

def ode_residual(which, S1, S2, S3, kappa, q, gauge, u):
    """|ODE[e^{kappa u} q]| at u, relative to the size of the individual terms"""
    d = len(q) - 1
    row = ode_row(which, S1, S2, S3, kappa, d, u, gauge)
    val = sum(c * r for c, r in zip(q, row))
    scale = sum(abs(c) * abs(u) ** n for n, c in enumerate(q)) or 1   # O(1) structural scale (terms of the ODE are O(1) at u=0.37)
    return abs(val) / scale

def q1_run(ch, st, th):
    S1 = pf_S1(sum(mpf(z) for z in ZETA), 3, th)
    out = []
    for k, s in enumerate(st):
        S2 = pf_S2(s['lam2'], th)
        S3 = pf_S3(s['lam3'], th)
        mu = s['sector']
        rec = dict(mu=mu, ev=s['ev'], dual=[], direct=[])
        for j in range(3):
            d, nul, sm, q = solve_poly_mp('D', S1, S2, S3, mpf(ZETA[j]), th)
            rec['dual'].append((d, nul, sm, q))
        for i in range(3):
            d, nul, sm, q = solve_poly_mp('A', S1, S2, S3, -mpf(ZETA[i]), th, gauge=GAUGE)
            rec['direct'].append((d, nul, sm, q))
        # residuals at u = 0.37 (away from poles)
        u = mp.mpf('0.37')
        rr = []
        for j in range(3):
            d, nul, sm, q = rec['dual'][j]
            rr.append(ode_residual('D', S1, S2, S3, mpf(ZETA[j]), q, 0, u))
        for i in range(3):
            d, nul, sm, q = rec['direct'][i]
            rr.append(ode_residual('A', S1, S2, S3, -mpf(ZETA[i]), q, GAUGE, u))
        rec['resid'] = max(rr)
        out.append(rec)
    return out

def q1_negative_controls(ch, st, th):
    """one state per sector: plain adjoint (gauge 0,1), D with -zeta, A with +zeta: smallest degree found (None = none up to 8)"""
    S1 = pf_S1(sum(mpf(z) for z in ZETA), 3, th)
    seen = set()
    res = {}
    for s in st:
        if s['sector'] in seen:
            continue
        seen.add(s['sector'])
        S2 = pf_S2(s['lam2'], th); S3 = pf_S3(s['lam3'], th)
        r = {}
        for label, which, sg, gg in (('D,-zeta', 'D', -1, 0), ('A g0,-zeta', 'A', -1, 0), ('A g1,-zeta', 'A', -1, 1),
                                     ('A g3,-zeta', 'A', -1, 3), ('A g2,+zeta', 'A', 1, 2)):
            r[label] = [solve_poly_mp(which, S1, S2, S3, sg * mpf(ZETA[i]), th, gauge=gg)[0] for i in range(3)]
        res[s['sector']] = r
    return res

def vac_data(ch, om):
    """exact eigenvalue functions of the bra vacuum for bra twist om, plus exact Q's (Fractions)"""
    cov = vacuum_covector(ch, ZETA, om)
    gp = gprime(ZETA, om, 0)
    ops = Operators(ch, gp)
    C = ops.fit()
    ok, res = covector_eigen_test(C, cov, ch)
    th = list(ch.theta)
    S1 = PF(sum(om), {(0, 1): Fr(3), (1, 1): Fr(3)}, th)
    S2 = pf_S2(res[2], th)
    S3 = pf_S3(res[3], th)
    Qd = []
    for i in range(3):
        d, nul, q = solve_poly_exact('A', S1, S2, S3, -om[i], th, gauge=GAUGE)
        Qd.append((d, nul, q))
    Pd = []
    for j in range(3):
        d, nul, q = solve_poly_exact('D', S1, S2, S3, om[j], th)
        Pd.append((d, nul, q))
    return dict(cov=cov, ok=ok, lam2=res[2], lam3=res[3], Qd=Qd, Pd=Pd, S1=S1, S2=S2, S3=S3, gp=gp)

def fmt_pf(lam, kind):
    names2 = ['c', 'H1', 'a1', 'H2', 'a2']
    names3 = ['c', 'b1', 'a1', 'T1', 'b2', 'a2', 'T2']
    return lam


# ================================================================== Q2
def poly_from(rec_item):
    d, nul, sm, q = rec_item
    return q

def wronskian_direct(qa, za, qb, zb):
    """W[Q_a,Q_b] = e^{-(za+zb)u} * w(u), Q_i = e^{-z_i u} q_i"""
    t1 = pmul(qa, padd(pder(qb), pscale(qb, -zb)))
    t2 = pmul(qb, padd(pder(qa), pscale(qa, -za)))
    return padd(t1, pscale(t2, -1))

def wronskian_dual(pa, za, pb, zb):
    """W[P^a,P^b] = e^{+(za+zb)u} * w(u), P^j = e^{z_j u} p_j"""
    t1 = pmul(pa, padd(pder(pb), pscale(pb, zb)))
    t2 = pmul(pb, padd(pder(pa), pscale(pa, za)))
    return padd(t1, pscale(t2, -1))

def q2_wronskians(recs, th):
    """returns dict with deviations; for every state and pair (a<b), c = complement"""
    z = [mpf(x) for x in ZETA]
    Pi = pmul([-th[0], 1], [-th[1], 1])
    rows = []
    for rec in recs:
        dq = [x[3] for x in rec['direct']]
        dp = [x[3] for x in rec['dual']]
        for (a, b, c) in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):
            # direct Wronskian / dual P^c
            w = wronskian_direct(dq[a], z[a], dq[b], z[b])
            quo, rem = pdivmod(trim(w), dp[c])
            remn = max(abs(x) for x in rem) / max(abs(x) for x in w)
            lead = quo[-1]
            F = [x / lead for x in quo]
            dev = max(abs(F[i] - Pi[i]) for i in range(3)) if len(F) == 3 else mp.mpf(1)
            # dual Wronskian / direct q_c
            w2 = wronskian_dual(dp[a], z[a], dp[b], z[b])
            quo2, rem2 = pdivmod(trim(w2), dq[c])
            remn2 = max(abs(x) for x in rem2) / max(abs(x) for x in w2)
            lead2 = quo2[-1]
            F2 = [x / lead2 for x in quo2]
            dev2 = max(abs(F2[i] - Pi[i]) for i in range(3)) if len(F2) == 3 else mp.mpf(1)
            rows.append(dict(mu=rec['mu'], pair=(a, b), c=c, deg_w=len(trim(w)) - 1, deg_pc=len(dp[c]) - 1, rem=remn, lead=lead, F_dev=dev,
                             deg_w2=len(trim(w2)) - 1, rem2=remn2, lead2=lead2, F2_dev=dev2,
                             cab=lead / (z[a] - z[b]), cab2=lead2 / (z[b] - z[a]), Fs=F, F2s=F2))
    return rows

def det3(M):
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))

def q2_weyl(recs, st, ch, th, upts=('0.37', '1.9', '-1.3')):
    """Cramer / Wronskian-determinant formulas for S1,S2,S3 from the three dual Q's; also direct side"""
    z = [mpf(x) for x in ZETA]
    S1f = pf_S1(sum(z), 3, th)
    out = []
    for rec, s in zip(recs, st):
        S2 = pf_S2(s['lam2'], th); S3 = pf_S3(s['lam3'], th)
        P = [Qfun(z[j], rec['dual'][j][3]) for j in range(3)]
        Qd = [Qfun(-z[i], rec['direct'][i][3]) for i in range(3)]
        dev = dict(S1=0, S2=0, S3=0, S1d=0, S2d=0, S3d=0, Wconst=None)
        for us in upts:
            u = mp.mpf(us)
            D = [[P[j](u, k) for k in range(4)] for j in range(3)]
            W = det3([[D[j][0], D[j][1], D[j][2]] for j in range(3)])
            s3 = det3([[D[j][3], D[j][1], D[j][2]] for j in range(3)]) / W
            s2 = -det3([[D[j][0], D[j][3], D[j][2]] for j in range(3)]) / W
            s1 = det3([[D[j][0], D[j][1], D[j][3]] for j in range(3)]) / W
            dev['S1'] = max(dev['S1'], abs(s1 - S1f(u)) / abs(S1f(u)))
            dev['S2'] = max(dev['S2'], abs(s2 - S2(u)) / abs(S2(u)))
            dev['S3'] = max(dev['S3'], abs(s3 - S3(u)) / abs(S3(u)))
            Pi3 = ((u - th[0]) * (u - th[1])) ** 3
            wc = W / (mp.exp(sum(z) * u) * Pi3)
            dev['Wconst'] = wc if dev['Wconst'] is None else max(dev['Wconst'], wc, key=abs) if False else wc
            # direct: psi_i = Pi^{-2} Q_i, monic operator d^3 - S1' d^2 + S2' d - S3' with S1'=-S1, S2'=2S1'+S2, S3'=-(S1''+S2'+S3)
            def psi(uu, k):
                # derivatives of psi = Pi^{-GAUGE} Q via series
                out = []
                return None
            # build series
            Dd = []
            for i in range(3):
                qs = [Qd[i](u, k) / factorial(k) for k in range(4)]
                ps = qs
                for t in th:
                    ps = ser_mul(ps, ser_pow(u - t, -GAUGE))
                Dd.append([ps[k] * factorial(k) for k in range(4)])
            Wd = det3([[Dd[i][0], Dd[i][1], Dd[i][2]] for i in range(3)])
            t3 = det3([[Dd[i][3], Dd[i][1], Dd[i][2]] for i in range(3)]) / Wd
            t2 = -det3([[Dd[i][0], Dd[i][3], Dd[i][2]] for i in range(3)]) / Wd
            t1 = det3([[Dd[i][0], Dd[i][1], Dd[i][3]] for i in range(3)]) / Wd
            e1 = -S1f(u)
            e2 = 2 * S1f(u, 1) + S2(u)
            e3 = -(S1f(u, 2) + S2(u, 1) + S3(u))
            dev['S1d'] = max(dev['S1d'], abs(t1 - e1) / abs(e1))
            dev['S2d'] = max(dev['S2d'], abs(t2 - e2) / abs(e2))
            dev['S3d'] = max(dev['S3d'], abs(t3 - e3) / abs(e3))
        out.append(dev)
    return out


# ================================================================== Q3
COLS = ['H1', 'H2', 'a1', 'a2', 'b1', 'b2']          # unknown Laurent coefficients of Delta S2 (simple poles), Delta S3 (double, simple)

def qtaylor(Qf, x, n=8):
    """Q^{(k)}(x)/k!  k = 0..n-1 as list"""
    return [Qf(x, k) / factorial(k) for k in range(n)]

class SiteData:
    """Taylor data of a set of Q functions at theta_s, with derivatives"""
    def __init__(self, Qfuns, th, nder=6):
        self.t = [[qtaylor(Q, th[s], nder + 4) for Q in Qfuns] for s in range(2)]   # t[s][idx][k]

def site_series(td, s, idx, a, n=4):
    """Taylor series (in x = u - theta_s) of the a-th derivative of function idx"""
    c = td.t[s][idx]
    # f^{(a)}(theta+x) = sum_k f^{(a+k)}/k! x^k ; stored c[m] = f^{(m)}/m!  ->  f^{(a+k)}/k! = c[a+k]*(a+k)!/k!
    return [c[a + k] * factorial(a + k) / factorial(k) for k in range(n)]

def res_functional(tdP, jdx, tdQ, idx, s, a, m, sp, th, n=4):
    """Res_{u=theta_s}[ Pi^{-2} * Q_idx * P_jdx^{(a)} / (u - theta_sp)^m ]  (sp = site of the pole factor)"""
    other = 1 - s
    Ps = site_series(tdP, s, jdx, a, n)
    Qs = site_series(tdQ, s, idx, 0, n)
    G = ser_mul(Ps, Qs, n)
    G = ser_mul(G, ser_pow(th[s] - th[other], -GAUGE, n), n)
    ms = 0
    if sp == s:
        ms = m
    else:
        G = ser_mul(G, ser_pow(th[s] - th[sp], -m, n), n)
    # integrand = x^{-GAUGE - ms} * G(x)  -> residue = coefficient of x^{GAUGE+ms-1}
    k = GAUGE + ms - 1
    return G[k]

def Nrows(tdP, tdQ, th):
    """dict (s,i,j) -> [N entries for COLS]; i indexes direct Q (bra), j dual P (ket)"""
    out = {}
    for s in range(2):
        for i in range(3):
            for j in range(3):
                row = []
                for col in COLS:
                    kind, sp = col[0], int(col[1]) - 1
                    if kind == 'H':
                        v = -res_functional(tdP, j, tdQ, i, s, 1, 1, sp, th)
                    elif kind == 'a':
                        v = res_functional(tdP, j, tdQ, i, s, 0, 2, sp, th)
                    else:
                        v = res_functional(tdP, j, tdQ, i, s, 0, 1, sp, th)
                    row.append(v)
                out[(s, i, j)] = row
    return out

def state_Q_objects(rec):
    z = [mpf(x) for x in ZETA]
    P = [Qfun(z[j], rec['dual'][j][3]) for j in range(3)]
    Q = [Qfun(-z[i], rec['direct'][i][3]) for i in range(3)]
    return P, Q

def delta_vec(sB, sA):
    """Delta = S^B - S^A Laurent data in COLS order"""
    l2B, l3B, l2A, l3A = sB['lam2'], sB['lam3'], sA['lam2'], sA['lam3']
    return [l2B[1] - l2A[1], l2B[3] - l2A[3], l3B[2] - l3A[2], l3B[5] - l3A[5], l3B[1] - l3A[1], l3B[4] - l3A[4]]


# ================================================================== reporting
def hdr(t):
    print('\n' + '=' * 78 + '\n' + t + '\n' + '=' * 78)

def nst(x, n=3):
    return mp.nstr(x, n)

def is_scalar_op(op, n):
    vals = {v for r_, row in op.r.items() for c_, v in row.items() if r_ == c_}
    offd = any(r_ != c_ for r_, row in op.r.items() for c_ in row)
    return (not offd) and len(vals) == 1 and len(op.r) == n, vals

def section0(ch):
    hdr('SECTION 0: conventions, symbolic cdet, chain checks')
    print('Lambda =', ch.lam, ' theta =', [str(t) for t in ch.theta], ' dim =', ch.dim)
    idx, ok = ch.check_hw()
    print("HW ket index in gl3_reps ordering (0-based) =", idx, "(NOT 0); GT pattern", ch.rep.patterns[idx],
          "; annihilated by E12,E13,E23:", ok)
    S1, S2, S3 = symbolic_cdet()
    print('symbolic cdet = sum_sigma sgn(sigma) M_{s(1)1} M_{s(2)2} M_{s(3)3}, M = d - L, L_ij = g_ij + sum_s E^(s)_ij/(u-theta_s), d f = f d + f\':')
    print('  S1 =', sym_str(S1))
    print('  S2 =', sym_str(S2))
    print('  S3 =', sym_str(S3))
    print("  (l_ij = L_ij(u), l'_ij = dL_ij/du, l''_ij = d^2L_ij/du^2; monomials are ordered products, left to right)")
    g = [[ZETA[i] if i == j else Fr(0) for j in range(3)] for i in range(3)]
    for (tr, ro) in ((False, False), (False, True), (True, False), (True, True)):
        ops = Operators(ch, g, transpose=tr, rowordered=ro)
        a = ops.S_at(Fr(3, 7))[1]
        b = ops.S_at(Fr(-2, 9))[2]
        print('  variant transposed-E=%s row-ordered=%s: [S2(3/7),S3(-2/9)] == 0 exactly: %s' % (tr, ro, comm(a, b).is_zero()))
    ops = Operators(ch, g)
    C = ops.fit()
    print('partial-fraction fit of S2 (basis 1, 1/x_s, 1/x_s^2) and S3 (1, 1/x_s .. 1/x_s^3), x_s = u - theta_s: mismatches at extra sample points:', C['res2'], C['res3'])
    print('S1 = tr g + 3/x_1 + 3/x_2 (times identity), exact:', (C['S1'] - SM.eye(ch.dim, S1_scalar(ch, ZETA, Fr(3, 7)))).is_zero())
    pts = [Fr(5, 3), Fr(-4, 5), Fr(7, 2)]
    allc = True
    for ua in pts:
        for ub in pts:
            A2 = ops.eval_S(C, 2, ua); B3 = ops.eval_S(C, 3, ub); B2 = ops.eval_S(C, 2, ub); A3 = ops.eval_S(C, 3, ua)
            allc = allc and comm(A2, B3).is_zero() and comm(A2, B2).is_zero() and comm(A3, B3).is_zero()
    print('all [S_k(u),S_l(v)] = 0 exactly on a 3x3 rational grid:', allc)
    wc = True
    for k in (1, 2, 3):
        Ek = ch.Es[(1, k, k)] + ch.Es[(2, k, k)]
        wc = wc and comm(ops.eval_S(C, 2, Fr(3, 7)), Ek).is_zero() and comm(ops.eval_S(C, 3, Fr(-2, 9)), Ek).is_zero()
    print('weight conservation [S2(u),E_kk^tot] = [S3(u),E_kk^tot] = 0 (k=1,2,3):', wc)
    e2 = ZETA[0] * ZETA[1] + ZETA[0] * ZETA[2] + ZETA[1] * ZETA[2]
    sc = is_scalar_op(C[2][0], ch.dim)
    print('S2 constant term is scalar =', sc[1], ' (e2(zeta) =', e2, ')')
    print('S2 double-pole coefficients scalar:', is_scalar_op(C[2][2], ch.dim)[0], is_scalar_op(C[2][4], ch.dim)[0], is_scalar_op(C[2][2], ch.dim)[1])
    sc = is_scalar_op(C[3][0], ch.dim)
    print('S3 constant term is scalar =', sc[1], ' (e3(zeta) =', ZETA[0] * ZETA[1] * ZETA[2], ')')
    print('S3 triple-pole coefficient operators (u-theta_s)^-3 are identically ZERO:', C[3][3].is_zero(), C[3][6].is_zero())

    def Hs(s):
        sp = 3 - s
        out = SM(ch.dim)
        for i in (1, 2, 3):
            for j in (1, 2, 3):
                if i < j:
                    out = out + ch.Es[(s, j, j)] * g[i - 1][i - 1] + ch.Es[(s, i, i)] * g[j - 1][j - 1]
                    out = out - ch.Es[(s, j, i)] * g[i - 1][j - 1] - ch.Es[(s, i, j)] * g[j - 1][i - 1]
                    cr = (ch.Es[(s, i, i)] * ch.Es[(sp, j, j)] + ch.Es[(s, j, j)] * ch.Es[(sp, i, i)]
                          - ch.Es[(s, i, j)] * ch.Es[(sp, j, i)] - ch.Es[(s, j, i)] * ch.Es[(sp, i, j)])
                    out = out + cr * (1 / (ch.theta[s - 1] - ch.theta[sp - 1]))
        return out
    print("simple-pole coefficient of S2 at theta_s equals H_s = sum_{i<j}[g_ii E_jj + g_jj E_ii - g_ij E_ji - g_ji E_ij]^(s)")
    print("   + 1/(theta_s - theta_s') sum_{i<j}[E^(s)_ii E^(s')_jj + E^(s)_jj E^(s')_ii - E^(s)_ij E^(s')_ji - E^(s)_ji E^(s')_ij]:",
          (C[2][1] - Hs(1)).is_zero(), (C[2][3] - Hs(2)).is_zero())
    return ops, C

def section_vacuum(ch):
    hdr("SECTION 0b: bra vacuum and g' orientation")
    r = gauss_r(ZETA, OM)
    print('Gauss parameters r_k(zeta,omega) =', [str(x) for x in r])
    out = {}
    for name, om in (('omega', OM), ('omega_prime', OMP)):
        cov = vacuum_covector(ch, ZETA, om)
        print('bra twist', name, [str(x) for x in om])
        for var, label in enumerate(['P gc P^-1', 'P^-1 gc P', '(P gc P^-1)^T', '(P^-1 gc P)^T', 'P^T gc P^-T', 'P^-T gc P^T']):
            ops = Operators(ch, gprime(ZETA, om, var))
            C = ops.fit()
            ok, res = covector_eigen_test(C, cov, ch)
            print("   variant %d  g' = %-14s : Omega x Omega is a left eigencovector of all S2,S3 coefficient operators: %s" % (var, label, ok))
        vd = vac_data(ch, om)
        out[name] = vd
        l2, l3 = vd['lam2'], vd['lam3']
        print('   exact eigenvalue functions (x1 = u-theta1, x2 = u-theta2):')
        print('   S1 = %s + 3/x1 + 3/x2' % sum(om))
        print('   S2 = %s + (%s)/x1 + (%s)/x1^2 + (%s)/x2 + (%s)/x2^2' % tuple(str(x) for x in l2))
        print('   S3 = %s + (%s)/x1 + (%s)/x1^2 + (%s)/x1^3 + (%s)/x2 + (%s)/x2^2 + (%s)/x2^3' % tuple(str(x) for x in l3))
        e2 = om[0] * om[1] + om[0] * om[2] + om[1] * om[2]
        print('   e2(omega) =', e2, ' e3(omega) =', om[0] * om[1] * om[2])
        print('   exact vacuum direct Q_i = e^{-omega_i u} q_i (gauge-%d adjoint): degrees' % GAUGE, [x[0] for x in vd['Qd']], 'nullities', [x[1] for x in vd['Qd']])
        print('   exact vacuum dual   P^j = e^{+omega_j u} p_j: degrees', [x[0] for x in vd['Pd']], 'nullities', [x[1] for x in vd['Pd']])
        for i in range(3):
            print('     q_%d =' % (i + 1), [str(c) for c in vd['Qd'][i][2]], '(coefficients of u^0, u^1, ...)')
    return out

def section1(ch, st, th, recs):
    hdr('SECTION 1 (Q1): opers and Q-functions, 64 states, 60 digits')
    print('  dual ODE  (D):  [d^3 - S1 d^2 + S2 d - S3] P = 0,     P^j = e^{+zeta_j u} p_j,  deg p_j = mu_j')
    print('  direct ODE:     Q_i = Pi^2 psi_i, A psi_i = 0 (A = formal adjoint of D), Q_i = e^{-zeta_i u} q_i, deg q_i = 4 - mu_i = 2 lam1 - mu_i, Pi = (u-th1)(u-th2)')
    bys = {}
    for r in recs:
        bys.setdefault(r['mu'], []).append(r)
    print('  %-10s %-4s %-18s %-20s %-8s %-10s %-9s' % ('sector mu', '#st', 'deg dual j=1,2,3', 'deg direct i=1,2,3', 'check', 'nullities', 'max resid'))
    allok = True
    for mu in sorted(bys):
        rs = bys[mu]
        dd = {tuple(x[0] for x in r['dual']) for r in rs}
        di = {tuple(x[0] for x in r['direct']) for r in rs}
        nl = {tuple(x[1] for x in r['dual'] + r['direct']) for r in rs}
        ok = dd == {tuple(mu)} and di == {tuple(4 - m for m in mu)} and nl == {(1,) * 6}
        allok = allok and ok
        print('  %-10s %-4d %-18s %-20s %-8s %-10s %-9s' % (mu, len(rs), sorted(dd)[0] if len(dd) == 1 else dd, sorted(di)[0] if len(di) == 1 else di,
                                                  'OK' if ok else 'FAIL', 'all 1' if nl == {(1,) * 6} else nl, nst(max(r['resid'] for r in rs))))
    print('  ALL 64 states: degrees as expected, kernel dimension 1 for each of the 6 exponentials:', allok)
    print('  max ODE residual (u=0.37, 60 digits) over 64 states x 6 Q:', nst(max(r['resid'] for r in recs)),
          '; largest smallest-singular-value at the degree found:', nst(max(x[2] for r in recs for x in r['dual'] + r['direct'])))

def section1_neg(ch, st, th):
    S = q1_negative_controls(ch, st, th)
    print('  negative controls (one state per sector, 19 states; smallest polynomial degree found for i=1,2,3; None = no solution up to degree 8):')
    names = {'D,-zeta': 'dual ODE D with e^{-zeta u}', 'A g0,-zeta': 'plain formal adjoint A (gauge 0), e^{-zeta u}',
             'A g1,-zeta': 'A with Pi^-1 gauge, e^{-zeta u}', 'A g3,-zeta': 'A with Pi^-3 gauge, e^{-zeta u}', 'A g2,+zeta': 'A with Pi^-2 gauge, e^{+zeta u}'}
    for l, nm in names.items():
        vals = {tuple(S[mu][l]) for mu in S}
        print('    %-45s : %s' % (nm, sorted(vals, key=str)[:4]))

def section2(ch, st, th, recs):
    hdr('SECTION 2 (Q2): Wronskian / QQ relation and quantum-Weyl formulas')
    rows = q2_wronskians(recs, th)
    z = [mpf(x) for x in ZETA]
    print('  (a) W[Q_a,Q_b] = Q_a Q_b\' - Q_b Q_a\' = c_ab e^{-(z_a+z_b)u} F(u) p_c  (c = complement of {a,b}; all p, q monic)  <=>  W/P^c = c_ab e^{-(z1+z2+z3)u} F(u)')
    print('      over 64 states x 3 pairs: max relative remainder of w_ab / p_c =', nst(max(r['rem'] for r in rows)))
    print('      quotient degree = deg w - deg p_c in', sorted({r['deg_w'] - r['deg_pc'] for r in rows}))
    print('      max |F - (u-th1)(u-th2)| (F = monic quotient) =', nst(max(r['F_dev'] for r in rows)))
    print('      max |c_ab/(zeta_a - zeta_b) - 1| =', nst(max(abs(r['cab'] - 1) for r in rows)))
    print('      F = u^2 + (%s) u + (%s);   (u-th1)(u-th2) = u^2 + %s u + %s' % (nst(rows[0]['Fs'][1], 12), nst(rows[0]['Fs'][0], 12),
                                                                           nst(-(th[0] + th[1]), 12), nst(th[0] * th[1], 12)))
    print("  (a') W[P^a,P^b] = c'_ab e^{+(z_a+z_b)u} F(u) q_c")
    print('      max relative remainder =', nst(max(r['rem2'] for r in rows)), '; max |F\' - Pi| =', nst(max(r['F2_dev'] for r in rows)),
          '; max |c\'_ab/(zeta_b - zeta_a) - 1| =', nst(max(abs(r['cab2'] - 1) for r in rows)))
    print('      sample rows (sector, pair (a,b), c, deg w, deg p_c, c_ab/(za-zb), rel. remainder):')
    for r in rows[:6]:
        print('        ', r['mu'], tuple(x + 1 for x in r['pair']), r['c'] + 1, r['deg_w'], r['deg_pc'], nst(r['cab'], 10), nst(r['rem']))
    wy = q2_weyl(recs, st, ch, th)
    print("  (b) quantum-Weyl formulas from the three dual Q's (the unique monic third-order operator killing P^1,P^2,P^3):")
    print("      W = det[P^j, P^j', P^j''],  S1 = det[P,P',P''']/W = (ln W)',  S2 = -det[P,P''',P'']/W,  S3 = det[P''',P',P'']/W  (rows j)")
    for k in ('S1', 'S2', 'S3'):
        print('      max relative deviation from the eigenvalue function %s (64 states x 3 u-values): %s' % (k, nst(max(w[k] for w in wy))))
    Wc = [w['Wconst'] for w in wy]
    vand = (z[1] - z[0]) * (z[2] - z[0]) * (z[2] - z[1])
    print('      W = c_W e^{(z1+z2+z3)u} Pi^3 with monic p_j; c_W state-independent: max |c_W/prod_{i<j}(zeta_j-zeta_i) - 1| =', nst(max(abs(c / vand - 1) for c in Wc)),
          '(prod =', nst(vand, 10), ')')
    print("      direct side: psi_i = Pi^-2 Q_i span ker(adjoint) = ker of the monic operator d^3 + (-S1) d^2 + (2S1'+S2) d - (S1''+S2'+S3):")
    for k in ('S1d', 'S2d', 'S3d'):
        print('      max relative deviation %s: %s' % (k, nst(max(w[k] for w in wy))))


def sv_ratio(N, ncols=6, keys=None):
    keys = keys or sorted(N)
    M = mp.matrix(len(keys), ncols)
    for r, key in enumerate(keys):
        for c in range(ncols):
            M[r, c] = N[key][c]
    for c in range(ncols):
        mx = max(abs(M[i, c]) for i in range(M.rows)) or 1
        for i in range(M.rows):
            M[i, c] /= mx
    U, S, V = mp.svd_r(M)
    return [S[i] / S[0] for i in range(ncols)]

def build_N_det(N, trip, cols=6):
    M = mp.matrix(2 * len(trip), cols)
    r = 0
    for s in range(2):
        for (i, j) in trip:
            for c in range(cols):
                M[r, c] = N[(s, i, j)][c]
            r += 1
    return M

def naive_rows(tdP, tdQ, th):
    """naive per-site brackets <f>_s = Res f/(u-th_s)^3 (Taylor coeff 2 of f), f = P Q, P Q', P' Q, with the raw direct Q (no Pi^-2)"""
    out = {}
    for s in range(2):
        for i in range(3):
            for j in range(3):
                row = []
                for s2 in range(2):
                    if s2 != s:
                        row += [0, 0, 0]
                        continue
                    P0 = site_series(tdP, s, j, 0, 4); P1 = site_series(tdP, s, j, 1, 4)
                    Q0 = site_series(tdQ, s, i, 0, 4); Q1 = site_series(tdQ, s, i, 1, 4)
                    row += [ser_mul(P0, Q0)[2], ser_mul(P0, Q1)[2], ser_mul(P1, Q0)[2]]
                out[(s, i, j)] = row
    return out

def section3(ch, st, th, recs):
    import random
    hdr('SECTION 3 (Q3): residue functionals, orthogonality, det N, x0-free overlap ratio')
    objs = [state_Q_objects(r) for r in recs]
    tdP = [SiteData(o[0], th) for o in objs]
    tdQ = [SiteData(o[1], th) for o in objs]
    rng = random.Random(7)
    pairs_ab = [(a, b) for a in range(64) for b in range(64) if a != b]
    sample = rng.sample(pairs_ab, 300)
    print('  Brackets used: <f>_s = Res_{u=theta_s} f/(u-theta_s)^3 = f\'\'(theta_s)/2 only enter through Taylor data; the exact identity below is derived, not fitted.')
    print("  Lagrange identity: psi D_B[P] - P D_B^dag[psi] = C(P,psi)', with P = P^j (ket B, D_B P = 0) and psi = Pi^-2 Q_i (bra A, A_A psi = 0):")
    print("     P (dS2 psi)' + dS3 P psi = C', dS_k = S_k^B - S_k^A  (S1 equal, constants equal)  =>  Res_{theta_s}[...] = 0 for s = 1, 2, giving")
    print("     sum_m Delta_m N[(s,i,j),m] = 0,  Delta = (dH1, dH2, da1, da2, db1, db2) = simple-pole coeffs of dS2 and double-/simple-pole coeffs of dS3 (state differences),")
    print("     N[(s,i,j),H_t] = -Res_{th_s}[ psi P^j' /(u-th_t) ],  N[.,a_t] = Res_{th_s}[ psi P^j /(u-th_t)^2 ],  N[.,b_t] = Res_{th_s}[ psi P^j /(u-th_t) ]   (t = 1,2).")
    mx = mp.mpf(0)
    t0 = time.time()
    rank_ne = []
    for (A, B) in sample:
        N = Nrows(tdP[B], tdQ[A], th)
        D = delta_vec(st[B], st[A])
        for key, row in N.items():
            v = abs(sum(r_ * d for r_, d in zip(row, D)))
            sc = sum(abs(r_ * d) for r_, d in zip(row, D))
            if sc > 0:
                mx = max(mx, v / sc)
    print('  (a) identity residual max |sum_m Delta_m N| / sum_m |Delta_m N| over 300 random ordered pairs A != B x 18 rows (%.0fs): %s' % (time.time() - t0, nst(mx)))
    for (A, B) in sample[:60]:
        S = sv_ratio(Nrows(tdP[B], tdQ[A], th))
        rank_ne.append(S)
    print('      18x6 stacked N (rows (s,i,j), 9 pairs x 2 sites; columns column-normalised), 60 pairs A != B:')
    print('        max sigma_6/sigma_1 =', nst(max(s_[5] for s_ in rank_ne)), '; min sigma_5/sigma_1 =', nst(min(s_[4] for s_ in rank_ne)), ' => rank 5')
    diag = [sv_ratio(Nrows(tdP[B], tdQ[B], th)) for B in range(64)]
    print('        A = B (64 states): min sigma_6/sigma_1 =', nst(min(s_[5] for s_ in diag)), ' => rank 6')
    naive = [sv_ratio(naive_rows(tdP[B], tdQ[A], th), 6) for (A, B) in sample[:20]]
    nv = sorted(s_[5] for s_ in naive)
    nfull = [x for x in nv if x > mp.mpf(10) ** -30]
    print('      naive per-site columns {<PQ>_s, <PQ_x>_s, <P_x Q>_s}, s=1,2 (raw Q, no Pi^-2): of 20 pairs A != B, %d have sigma_6/sigma_1 > 1e-30 (full rank, median %s) and %d are deficient (generic pairs: does NOT close)' % (len(nfull), nst(nfull[len(nfull) // 2]), len(nv) - len(nfull)))
    # (b) det N
    xxx = ((0, 2), (0, 1), (1, 2))      # (1,3),(1,2),(2,3) 0-based (i,j)
    print('  (b) det N: 6x6, rows (s,(i,j)), (i,j) in {(1,3),(1,2),(2,3)} x s = 1,2 ; columns (H1,H2,a1,a2,b1,b2); entries from the above functionals')
    def rel_det(N):
        M = build_N_det(N, xxx)
        nr = mp.mpf(1)
        for r in range(6):
            nr *= mp.sqrt(sum(abs(M[r, c]) ** 2 for c in range(6)))
        return abs(mp.det(M)) / nr
    ne = max(rel_det(Nrows(tdP[B], tdQ[A], th)) for (A, B) in sample[:100])
    dg = min(rel_det(Nrows(tdP[B], tdQ[B], th)) for B in range(64))
    print('      |det N| / prod(row norms): max over 100 pairs A != B =', nst(ne), '; min over the 64 diagonal pairs A = B =', nst(dg))
    # (c) ratio test
    hdr('SECTION 3c: overlap-ratio identity R_B(w)/R_B(w\') = [D_B(w)/D_B(w\')] [D_B0(w\')/D_B0(w)], B0 = e1 x e1 (top state), 64 states')
    vac = {}
    for name, om in (('w', OM), ('wp', OMP)):
        vd = vac_data(ch, om)
        Qs = [Qfun(-mpf(om[i]), [mpf(c) for c in vd['Qd'][i][2]]) for i in range(3)]
        Ps = [Qfun(mpf(om[i]), [mpf(c) for c in vd['Pd'][i][2]]) for i in range(3)]
        vac[name] = dict(cov=vd['cov'], td=SiteData(Qs, th), tdP=SiteData(Ps, th))
    def overlap(cov, s):
        return sum(mpf(c) * x for c, x in zip(cov, s['vec']) if c != 0)
    Rw = [overlap(vac['w']['cov'], s) for s in st]
    Rp = [overlap(vac['wp']['cov'], s) for s in st]
    top = [k for k, s in enumerate(st) if s['sector'] == (4, 2, 0)][0]
    print('  top state B0 (weight (4,2,0)) has index', top, '; overlaps R_B(w), R_B(w\') nonzero for all 64:', all(abs(x) > mp.mpf(10) ** -30 for x in Rw + Rp))
    NC = {}
    def Nvac(B, name):
        if (B, name) not in NC:
            NC[(B, name)] = Nrows(tdP[B], vac[name]['td'], th)
        return NC[(B, name)]
    def ratio_test(trip, Nfun, cols=6):
        dw = [mp.det(build_N_det(Nfun(B, 'w'), trip, cols)) for B in range(64)]
        dp = [mp.det(build_N_det(Nfun(B, 'wp'), trip, cols)) for B in range(64)]
        devs = []
        nd = 0
        for B in range(64):
            L = Rw[B] * dp[B] * dw[top]
            R = Rp[B] * dw[B] * dp[top]
            m = max(abs(L), abs(R))
            if m > mp.mpf(10) ** -45:
                nd += 1
                devs.append(abs(L - R) / m)
            else:
                devs.append(None)
        good = [d for d in devs if d is not None]
        return (max(good) if good else None), nd, devs, dw, dp
    mdev, nd, devs, dw, dp = ratio_test(xxx, Nvac)
    print('  main result, rows {(1,3),(1,2),(2,3)} x 2 sites, columns (H1,H2,a1,a2,b1,b2), 60-digit arithmetic:')
    print('     max over 64 states of |R(w)D(w\')D0(w) / (R(w\')D(w)D0(w\')) - 1|-type relative deviation =', nst(mdev), ' (non-degenerate states:', nd, ')')
    print('     the 5 worst states (index, sector, dev):', [(k, st[k]['sector'], nst(d)) for k, d in sorted([(k, d) for k, d in enumerate(devs) if d is not None], key=lambda t: -t[1])[:5]])
    print('     sample: R(w)/R(w\') vs [D(w)/D(w\')][D0(w\')/D0(w)] (first 4 states):')
    for B in range(4):
        print('        B=%d %s  LHS = %s  RHS = %s' % (B, st[B]['sector'], nst(Rw[B] / Rp[B], 15), nst((dw[B] / dp[B]) * (dp[top] / dw[top]), 15)))
    # all row triples
    allpairs = [(i, j) for i in range(3) for j in range(3)]
    ok_t, bad_t, deg_t = [], [], []
    for trip in itertools.combinations(allpairs, 3):
        md, nd2, _, _, _ = ratio_test(trip, Nvac)
        if nd2 < 64:
            deg_t.append((trip, nd2))
        elif md < mp.mpf(10) ** -30:
            ok_t.append(trip)
        else:
            bad_t.append((trip, md))
    print('  all C(9,3) = 84 row-triples (same triple on both sites): pass (all 64 states non-degenerate, dev < 1e-30): %d ; fail: %d ; degenerate (some det = 0): %d' % (len(ok_t), len(bad_t), len(deg_t)))
    print('     failing triples:', [(t, nst(m)) for t, m in bad_t][:6])
    # controls
    def Nmod(B, name):
        N = Nvac(B, name)
        N2 = {}
        for key, row in N.items():
            s, i, j = key
            row2 = list(row)
            # replace column a2 by Res[psi P / (u-th2)^3]
            row2[3] = res_functional(tdP[B], j, vac[name]['td'], i, s, 0, 3, 1, th, 6)
            N2[key] = row2
        return N2
    NCm = {}
    def Nmod_c(B, name):
        if (B, name) not in NCm:
            NCm[(B, name)] = Nmod(B, name)
        return NCm[(B, name)]
    md2, nd2, _, _, _ = ratio_test(xxx, Nmod_c)
    print('  control 1 (column a2 replaced by Res[psi P/(u-th2)^3]): max deviation =', nst(md2) if md2 is not None else None, ' non-degenerate', nd2)
    def Nswap(B, name):
        # opposite orientation: ket = state's DIRECT Q's, bra = vacuum's DUAL P's (roles of Q and P exchanged in the functionals)
        return Nrows(vac[name]['tdP'], tdQ[B], th)
    NCs = {}
    def Nswap_c(B, name):
        if (B, name) not in NCs:
            NCs[(B, name)] = Nswap(B, name)
        return NCs[(B, name)]
    md3, nd3, _, _, _ = ratio_test(xxx, Nswap_c)
    print('  control 2 (orientation swapped: vacuum dual P^j in the ket slot, state direct Q_i in the bra slot): max deviation =', nst(md3) if md3 is not None else None, ' non-degenerate', nd3)
    return dict(xxx_dev=mdev, ok_triples=ok_t, bad=bad_t, deg=deg_t)

# ================================================================== Q4: untwisting regime
PHI = [Fr(1, 3), Fr(-2, 7), Fr(5, 11)]
OMH = [Fr(3, 4), Fr(-1, 5), Fr(1, 2)]        # third component not fixed by the plan: chosen arbitrarily (documented)
OMPH = [Fr(-12, 5), Fr(1, 8), Fr(-1, 3)]
HW_SECTORS = {(4, 2, 0): 1, (4, 1, 1): 1, (3, 3, 0): 1, (3, 2, 1): 2, (2, 2, 2): 1}

def set_regime(zeta, om, omp, dps, tol):
    global ZETA, OM, OMP, DPS, TOL
    ZETA, OM, OMP, DPS, TOL = zeta, om, omp, dps, mp.mpf(10) ** (-tol)
    mp.mp.dps = dps

def raising_residual(ch, s):
    v = s['vec']
    nv = mp.sqrt(sum(abs(x) ** 2 for x in v))
    tot = 0
    for (i, j) in ((1, 2), (2, 3)):
        E = ch.Es[(1, i, j)] + ch.Es[(2, i, j)]
        w = E.apply(v)
        tot += sum(abs(x) ** 2 for x in w)
    return mp.sqrt(tot) / nv

def q4_point(eps):
    """one regime point: returns dict with branches, det N data at this eps"""
    zeta = [eps * p for p in PHI]
    om = [eps * x for x in OMH]
    omp = [eps * x for x in OMPH]
    set_regime(zeta, om, omp, 100, 50)
    ch, ops, C, st, res = build_states()
    th = [mpf(t) for t in ch.theta]
    # branches
    br = []
    bysec = {}
    for k, s in enumerate(st):
        s['rr'] = raising_residual(ch, s)
        bysec.setdefault(s['sector'], []).append(k)
    gaps = {}
    for mu, m in HW_SECTORS.items():
        ks = sorted(bysec[mu], key=lambda k: st[k]['rr'])
        br += ks[:m]
        gaps[mu] = [nst(st[k]['rr']) for k in ks[:m + 1]]
    S1 = pf_S1(sum(mpf(z) for z in zeta), 3, th)
    out = dict(eps=eps, res=res, gaps=gaps, branches=br, sectors=[st[k]['sector'] for k in br])
    recs = {}
    for k in br:
        s = st[k]
        S2 = pf_S2(s['lam2'], th); S3 = pf_S3(s['lam3'], th)
        rec = dict(mu=s['sector'], dual=[], direct=[])
        for j in range(3):
            rec['dual'].append(solve_poly_mp('D', S1, S2, S3, mpf(zeta[j]), th))
        for i in range(3):
            rec['direct'].append(solve_poly_mp('A', S1, S2, S3, -mpf(zeta[i]), th, gauge=GAUGE))
        recs[k] = rec
    out['degrees'] = {k: ([x[0] for x in recs[k]['dual']], [x[0] for x in recs[k]['direct']], [x[1] for x in recs[k]['dual'] + recs[k]['direct']]) for k in br}
    # vacua (exact)
    vac = {}
    for name, o in (('w', om), ('wp', omp)):
        vd = vac_data(ch, o)
        Qs = [Qfun(-mpf(o[i]), [mpf(c) for c in vd['Qd'][i][2]]) for i in range(3)]
        vac[name] = dict(cov=vd['cov'], td=SiteData(Qs, th), ok=vd['ok'], deg=[x[0] for x in vd['Qd']])
    out['vac_ok'] = (vac['w']['ok'], vac['wp']['ok'], vac['w']['deg'], vac['wp']['deg'])
    out['r_k'] = gauss_r(zeta, om)
    objs = {k: state_Q_objects(recs[k]) for k in br}
    tdP = {k: SiteData(objs[k][0], th) for k in br}
    xxx = ((0, 2), (0, 1), (1, 2))
    def overlap(cov, s):
        return sum(mpf(c) * x for c, x in zip(cov, s['vec']) if c != 0)
    data = {}
    for k in br:
        d = {}
        for name in ('w', 'wp'):
            N = Nrows(tdP[k], vac[name]['td'], th)
            M = build_N_det(N, xxx)
            d['det_' + name] = mp.det(M)
            d['maxentry_' + name] = max(abs(M[r, c]) for r in range(6) for c in range(6))
            nr = mp.mpf(1)
            for r in range(6):
                nr *= mp.sqrt(sum(abs(M[r, c]) ** 2 for c in range(6)))
            d['relDet_' + name] = abs(d['det_' + name]) / nr
            nc = mp.mpf(1)
            for c in range(6):
                nc *= mp.sqrt(sum(abs(M[r, c]) ** 2 for r in range(6)))
            d['relDetCol_' + name] = abs(d['det_' + name]) / nc
            Mc = M.copy()
            for c in range(6):
                cn = mp.sqrt(sum(abs(M[r, c]) ** 2 for r in range(6)))
                for r in range(6):
                    Mc[r, c] /= cn
            Sg = mp.svd_r(Mc, compute_uv=False)
            d['svratio_' + name] = Sg[5] / Sg[0]
            d['R_' + name] = overlap(vac[name]['cov'], st[k])
        data[k] = d
    top = [k for k in br if st[k]['sector'] == (4, 2, 0)][0]
    for k in br:
        d = data[k]
        L = d['R_w'] * d['det_wp'] * data[top]['det_w']
        R = d['R_wp'] * d['det_w'] * data[top]['det_wp']
        d['dev'] = abs(L - R) / max(abs(L), abs(R))
    out['data'] = data
    out['top'] = top
    return out

def section4():
    hdr('SECTION 4 (Q4): untwisting regime zeta = eps*phi, omega = eps*om^, omega\' = eps*om^\', HW-limit branches, 100 digits')
    print('  phi =', [str(x) for x in PHI], ' om^ =', [str(x) for x in OMH], " om^' =", [str(x) for x in OMPH], '(third components chosen arbitrarily: not fixed by the plan)')
    saved = (ZETA, OM, OMP, DPS, TOL)
    pts = {}
    try:
        for eps in (Fr(1, 10 ** 6), Fr(2, 10 ** 6)):
            t0 = time.time()
            pts[eps] = q4_point(eps)
            o = pts[eps]
            print('  eps = %s (%.0fs): eigen-residual %s; Gauss r_k(zeta,omega) = %s' % (str(eps), time.time() - t0, nst(o['res']), [nst(mpf(x), 8) for x in o['r_k']]))
            print('     vacuum eigencovector test (omega, omega\'):', o['vac_ok'][:2], ' vacuum direct Q degrees', o['vac_ok'][2], o['vac_ok'][3])
            print('     raising residuals (smallest first, m+1 per HW sector: gap shows the branch selection):', o['gaps'])
            for k in o['branches']:
                dd = o['degrees'][k]
                d = o['data'][k]
                print('     branch idx %2d sector %s : deg dual %s, deg direct %s, nullities %s | |det N(w)| = %s, det/prod(col norms) = %s, sigma_6/sigma_1 = %s | ratio-test deviation = %s' % (
                    k, o['sectors'][o['branches'].index(k)], dd[0], dd[1], dd[2], nst(abs(d['det_w'])), nst(d['relDetCol_w']), nst(d['svratio_w']), nst(d['dev'])))
        e1, e2 = Fr(1, 10 ** 6), Fr(2, 10 ** 6)
        print('  epsilon-exponents p in f ~ eps^p from the two points eps = 1e-6, 2e-6, p = log2(f(2e-6)/f(1e-6)):')
        a, b = pts[e1], pts[e2]
        for k1, k2 in zip(a['branches'], b['branches']):
            if a['sectors'][a['branches'].index(k1)] != b['sectors'][b['branches'].index(k2)]:
                continue
            da, db = a['data'][k1], b['data'][k2]
            ex = lambda f: nst(mp.log(abs(db[f]) / abs(da[f])) / mp.log(2), 4)
            print('     sector %s: |det N(w)| p = %s, |det N(wp)| p = %s, det/prod(row norms) p = %s, det/prod(col norms) p = %s, sigma_6/sigma_1 p = %s, max|N(w)| p = %s, |R(w)| p = %s' % (
                a['sectors'][a['branches'].index(k1)], ex('det_w'), ex('det_wp'), ex('relDet_w'), ex('relDetCol_w'), ex('svratio_w'), ex('maxentry_w'), ex('R_w')))
        print('  (note: Q normalisation is monic, so absolute exponents of det N include the eps-dependence of the monic normalisation; the ratio test is normalisation-free)')
    finally:
        set_regime(*saved[:4], 0)
        mp.mp.dps = saved[3]
        globals()['TOL'] = saved[4]
    return pts


def main():
    t0 = time.time()
    ch = Chain()
    ops, C = section0(ch)
    section_vacuum(ch)
    ch, ops, C, st, res = build_states()
    th = [mpf(t) for t in ch.theta]
    print('\n64 eigenstates found per weight sector (19 sectors); max eigen-residual of all coefficient operators on all states:', nst(res))
    cache = os.path.join(HERE, 'qside_cache.pkl')
    if os.path.exists(cache):
        with open(cache, 'rb') as f:
            recs = pickle.load(f)
        assert all(r['mu'] == s['sector'] and abs(r['ev'] - s['ev']) < mp.mpf(10) ** -40 for r, s in zip(recs, st)), 'stale cache'
    else:
        recs = q1_run(ch, st, th)
        with open(cache, 'wb') as f:
            pickle.dump(recs, f)
    section1(ch, st, th, recs)
    section1_neg(ch, st, th)
    section2(ch, st, th, recs)
    section3(ch, st, th, recs)
    section4()
    print('\ntotal time %.0fs' % (time.time() - t0))

if __name__ == '__main__':
    main()
