"""D2: the Gaudin two-vacua Q-side formula as the formal h -> 0 limit of the XXX ten-product bracket formula.

Run top to bottom:   python -X utf8 twovac_limit.py      (about 2-3 minutes; prints everything quoted in RESULTS-twovac-limit.md)

Sections
  L0  parse the closed forms of the XXX notebook (cell 62 of ref/xxx-210-infinf-formula-cells.txt) and verify them; exact top-state
      check of the whole XXX identity at finite h (conventions of minor, Bk2, ladder points, (e,e'), kappa, x0)
  L1  every h-dependent ingredient: kappa(z,w), x0(z), vacuum Gauss parameters, h-order table
  L2  the formal limit: exact ladder constraint minor_s(-1,0) = 0, leading order N = 10, the explicit Gaudin formula
      (Pluecker basis, twisted-jet form, direct Q_1-basis)
  L3  numerical test on the 64 Gaudin states (several (omega, omega'), theta, zeta points)

Only exact arithmetic (Fraction / sympy) until L3 and the h-scaling checks of L1/L2 (mpmath).  Reads (never edits) xxx_chain.py, gaudin_chain.py,
qside.py of the earlier rounds.
"""
import os
import sys
import time
import pickle
from fractions import Fraction as Fr
from math import factorial
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))
import sympy as sp
from sympy import QQ
from sympy.parsing.mathematica import parse_mathematica
from sympy.polys.rings import ring
import mpmath as mp

import xxx_chain
from xxx_chain import Ser, gl3_generators

T0 = time.time()
REF = os.path.normpath(os.path.join(HERE, "..", "..", "Mathematica", "Gaudin", "ref", "xxx-210-infinf-formula-cells.txt"))
sS, tt = sp.symbols("sS tt")
SLOTS = [(-3, -2), (-3, -1), (-2, -1), (-3, 0), (-2, 0)]          # ladder-1/2 slot pairs that occur in the kernel laws
EE = {(0, 0): (7, 5), (0, 1): (7, 4), (1, 0): (6, 5), (1, 1): (6, 4)}   # (e, e') per block (c, c')
MPTS = [-3, -2, -1, 0]                                              # ladder points theta_s + m h


def say(*a):
    print(*a, flush=True)


def hdr(t):
    say("\n" + "=" * 110)
    say(t)
    say("=" * 110)


def tm():
    return "[%.0fs]" % (time.time() - T0)


# ======================================================================================================================
# L0: parse and verify
# ======================================================================================================================
def _grab(cell, marker):
    k = cell.index(marker)
    j = cell.index("{", k + len(marker))
    depth = 0
    for p in range(j, len(cell)):
        if cell[p] == "{":
            depth += 1
        elif cell[p] == "}":
            depth -= 1
            if depth == 0:
                return cell[j:p + 1]
    raise ValueError("unbalanced")


def load_closed_forms():
    txt = open(REF, encoding="utf-8").read()
    i0, i1 = txt.index("===== cell 62"), txt.index("===== cell 70")
    cell = txt[i0:i1]
    m1 = _grab(cell, "\"kernel laws {(c,c'), slot (m m'), slot (n n'), kappa(s, t)} (tt = t = theta12/h, sS = s = a'/a)\",")
    m2 = _grab(cell, '(6,4)"')
    K = parse_mathematica(" ".join(m1.split()))
    P = parse_mathematica(" ".join(m2.split()))
    Kd = {(tuple(k[0]), tuple(k[1]), tuple(k[2])): k[3] for k in K}
    return Kd, list(P)


KD, PRODUCTS = load_closed_forms()


def wedge(F, G):
    """(F^G)_{mm'} = F_m G_m' - F_m' G_m, m < m' in ladder order (-3,-2,-1,0)  [notebook wedgeV]"""
    return {(MPTS[i], MPTS[j]): F[i] * G[j] - F[j] * G[i] for i in range(4) for j in range(i + 1, 4)}


def run_L0():
    hdr("L0  parse and verify the closed forms (cell 62)")
    say("conventions read off the setup cells (cells 17-19, 49):")
    say("  ladder points theta_s + m h, m = -3,-2,-1,0 (mPts = {-1-lambda1,..,0});  minor_s(m,m') = Q2(m') Q3(m) - Q3(m') Q2(m), m < m'")
    say("  <Q^j F>_s = sum_m F_m Q^j(theta_s + m h);  Bk2_s[F,G] = <Q3 F><Q2 G> - <Q2 F><Q3 G> = sum_{m<m'} (F^G)_{mm'} minor_s(m,m')")
    say("  (F^G)_{mm'} = F_m G_m' - F_m' G_m;   F_m = a^m f_m(s,t), G_m = a^(m+e) g_m(s,t)  =>  (F^G)_{mm'} = a^(m+m'+e) (f^g)_{mm'}")
    say("  M_vac = minor_s(-3,-2) Q3(theta_s - h)   (cell 49: mVac = (LV[2,2]LV[1,3]-LV[2,3]LV[1,2]) LV[3,3], LV[s,k,j] = Q^j(theta_s+mPts[k] h))")
    say("  identity: <Phi_c|Psi_c>/<x0_c|Psi_c> * M1_vac M2_vac = sum_{cc'} Q3(theta1-ch) Q3(theta2-c'h) sum_k Bk2^(1)[F_k,G_k] Bk2^(2)[F'_k,G'_k]")
    say("parsed: %d kernel laws, %d products (blocks (0,0),(0,1),(1,0),(1,1) with 3,3,2,2 products)" % (len(KD), len(PRODUCTS)))
    cnt = {}
    for p in PRODUCTS:
        cnt[(int(p[0]), int(p[1]))] = cnt.get((int(p[0]), int(p[1])), 0) + 1
    say("  products per block:", cnt)
    blocks = {}
    for p in PRODUCTS:
        c, cp, f, g, fp, gp = p
        w1, w2 = wedge(f, g), wedge(fp, gp)
        B = blocks.setdefault((int(c), int(cp)), {})
        for a, va in w1.items():
            for b, vb in w2.items():
                B[(a, b)] = B.get((a, b), 0) + va * vb
    bad = 0
    nchk = 0
    for cc, B in blocks.items():
        for (a, b), v in B.items():
            kv = KD.get((cc, a, b), 0)
            nchk += 1
            if sp.cancel(sp.together(v - kv)) != 0:
                bad += 1
                say("  MISMATCH", cc, a, b)
        for key in KD:
            if key[0] == cc and (key[1], key[2]) not in B:
                bad += 1
                say("  kernel law without wedge entry", key)
    say("RESULT L0: sum_k (f_k^g_k)_{mm'} (f'_k^g'_k)_{nn'} reproduces the printed kernel law of ALL %d block entries (4 blocks x 6x6 slots) "
        "identically in (s,t): %d mismatches" % (nchk, bad))
    say("  slots that occur: ladder-1/2 pairs", SLOTS, " -- the slot (-1,0) never occurs (see L2: minor_s(-1,0) = 0 exactly)")
    # structure: polynomial in s, degree in t
    degs = [(sp.Poly(sp.fraction(sp.together(v))[0], sS, tt).degree(sS)) for v in KD.values()]
    say("  kernel laws are polynomials in s of degree <= %d; denominators only (t-1) t (t+1) [t-poles at 0, +-1]" % max(degs))
    # wedge-entry h-orders (numerical, at t = theta12/h, s = exp(-h delta))
    mp.mp.dps = 50
    hs = [mp.mpf(10) ** -8, mp.mpf(10) ** -9]
    ords = {}
    for ip, p in enumerate(PRODUCTS):
        c, cp, f, g, fp, gp = p
        rec = []
        for vec in (f, g, fp, gp):
            os_ = []
            for e in vec:
                if e == 0:
                    continue
                fn = sp.lambdify((sS, tt), e, "mpmath")
                vals = [fn(mp.exp(-h * (mp.mpf(7) / 4 - mp.mpf(5) / 2)), mp.mpf(31) / 35 / h) for h in hs]
                if vals[0] == 0 or vals[1] == 0:
                    continue
                os_.append(round(float(mp.log(abs(vals[1] / vals[0])) / mp.log(hs[1] / hs[0]))))
            rec.append((min(os_), max(os_)) if os_ else None)
        ords[ip + 1] = rec
    say("  h-orders (min,max over nonzero entries) of the wedge vectors (f,g,f',g') per product, s=exp(-h*delta), t=theta12/h:")
    for ip in ords:
        say("    product %2d  (c,c')=%s  %s" % (ip, (int(PRODUCTS[ip - 1][0]), int(PRODUCTS[ip - 1][1])), ords[ip]))
    say("  -> individual wedge entries are Laurent in h (orders -2 .. +2: 1/(s-1) = O(1/h) in products 2 and 7, and t = theta12/h in product 2); the wedge sums, i.e. the 40 kernel laws, are polynomial in s with t-degree 0, hence O(1): no Laurent pole survives")
    return blocks


# ----------------------------------------------------------------------------------------------------------------------
# L0b: exact top-state check of the full XXX identity (finite h, exact rationals)
# ----------------------------------------------------------------------------------------------------------------------
def run_L0b():
    hdr("L0b  exact check of the whole XXX identity on the top state e1 (x) e1 (finite h = 1/2, exact rationals)")
    rep, E1 = gl3_generators((2, 1, 0))
    R = sp.Rational

    def tosp(op):
        M = sp.zeros(8, 8)
        for i, row in op.r.items():
            for j, v in row.items():
                M[i, j] = R(v.numerator, v.denominator)
        return M

    Ee = {k: tosp(v) for k, v in E1.items()}
    I8 = sp.eye(8)

    def nexp(A):
        return I8 + sum((A ** k / sp.factorial(k) for k in range(1, 5)), sp.zeros(8, 8))

    wts = [[Ee[(k, k)][p, p] for k in (1, 2, 3)] for p in range(8)]

    def rhodiag(d):
        return sp.diag(*[d[0] ** w[0] * d[1] ** w[1] * d[2] ** w[2] for w in wts])

    def Pphys(zv):
        z1, z2, z3 = zv
        sV = [z1 ** 2 / ((z1 - z2) * (z1 - z3)), -(z2 / (z1 * (z2 - z3))), 1 / (z1 * z2)]
        rV = [-(z3 / z2), -(z3 / z1), (z2 ** 2 * (z3 - z1)) / (z1 ** 2 * (z2 - z3))]
        return (nexp(rV[0] * Ee[(3, 2)]) * nexp(rV[1] * Ee[(3, 1)]) * nexp(rV[2] * Ee[(2, 1)]) * rhodiag(sV)
                * nexp(Ee[(1, 2)] / z1) * nexp(-Ee[(1, 3)] / (z1 * z2)) * nexp((z1 + z2) / (z1 * z2) * Ee[(2, 3)]))

    e1 = sp.Matrix([[1] + [0] * 7])

    def gaussCov(al, be, ga):
        return e1 * nexp(al * Ee[(1, 2)]) * nexp(be * Ee[(1, 3)]) * nexp(ga * Ee[(2, 3)])

    h = R(1, 2)
    zv = [R(3, 2), R(5, 4), R(2, 3)]
    thv = [R(2, 7), R(-3, 5)]
    Lam = (2, 1, 0)

    def nuf(n, x):
        r = 1
        for t in thv:
            r *= (x - t + h * Lam[n - 1])
        return r

    x = sp.Symbol("x")
    t1 = lambda u: sum(zv[i] * nuf(i + 1, u) for i in range(3))
    t2 = lambda u: sum(zv[i] * zv[j] * nuf(i + 1, u) * nuf(j + 1, u - h) for i in range(3) for j in range(i + 1, 3))
    t3 = lambda u: zv[0] * zv[1] * zv[2] * nuf(1, u) * nuf(2, u - h) * nuf(3, u - 2 * h)
    nu3 = lambda u: nuf(3, u)
    # dual Baxter (cell 7): nu3(x)nu3(x-h)nu3(x-2h) Q(x+h) - t1 nu3(x-h)nu3(x-2h) Q(x) + t2 nu3(x-2h) Q(x-h) - t3 Q(x-2h) = 0, Q = z_j^{x/h} q
    c = sp.symbols("c0:3")
    q2 = lambda y: c[0] + c[1] * y + c[2] * y ** 2
    zj = zv[1]
    expr = nu3(x) * nu3(x - h) * nu3(x - 2 * h) * zj * q2(x + h) - t1(x) * nu3(x - h) * nu3(x - 2 * h) * q2(x) + t2(x) * nu3(x - 2 * h) * q2(x - h) / zj - t3(x) * q2(x - 2 * h) / zj ** 2
    sol = sp.solve(sp.Poly(sp.expand(expr), x).all_coeffs() + [c[2] - 1], c, dict=True)[0]
    q2s = q2(x).subs(sol)
    say("top state (4,2,0): eigenvalues t1 = sum z_i nu_i(u), t2 = sum_{i<j} z_i z_j nu_i(u) nu_j(u-h), t3 = z1 z2 z3 nu1(u) nu2(u-h) nu3(u-2h)  (checked against the chain in the scratch run)")
    say("  dual Q's: q^3 = 1,  q^2 =", sp.expand(q2s))

    def ladder(ths):
        return {m: (zv[1] ** m * q2s.subs(x, ths + m * h), zv[2] ** m) for m in MPTS}

    def minor(L, m, mp_):
        return L[mp_][0] * L[m][1] - L[mp_][1] * L[m][0]

    Pinv = Pphys(zv).inv()
    col = Pinv[:, 0]
    x1 = (e1 * col)[0]
    L1l, L2l = ladder(thv[0]), ladder(thv[1])
    say("  minor_s(-1,0) at both sites (must vanish exactly):", sp.simplify(minor(L1l, -1, 0)), sp.simplify(minor(L2l, -1, 0)))
    # (rho, alpha, beta) of the ladder from Baxter (t-values) versus from the actual Q's
    def rab_from_t(ths):
        rho = t1(ths) * nu3(ths - h) / t2(ths)
        den = t3(ths - h)
        al = (nu3(ths - h) * nu3(ths - 2 * h) * nu3(ths - 3 * h) - rho * t1(ths - h) * nu3(ths - 2 * h) * nu3(ths - 3 * h)) / den
        be = t2(ths - h) * nu3(ths - 3 * h) / den
        return rho, al, be

    def rab_from_Q(L):
        D = minor(L, -2, 0)
        return L[-1][0] / L[0][0], -minor(L, -3, -2) / D, minor(L, -3, 0) / D

    say("  ladder ratios (rho, alpha, beta) from t-values (Baxter at x=theta_s, theta_s-h) vs from the Q's:  site 1 equal: %s, site 2 equal: %s ; values site 1: %s"
        % (all(sp.simplify(u - v) == 0 for u, v in zip(rab_from_t(thv[0]), rab_from_Q(L1l))), all(sp.simplify(u - v) == 0 for u, v in zip(rab_from_t(thv[1]), rab_from_Q(L2l))),
           [float(u) for u in rab_from_Q(L1l)]))
    ok_all = True
    for (a, ap) in ((R(3, 5), R(-2, 7)), (R(2, 3), R(5, 4)), (R(1, 3), R(1, 5))):
        s_, t_ = ap / a, (thv[0] - thv[1]) / h
        tot = 0
        for (cc, sl1, sl2), kv in KD.items():
            e, ep = EE[cc]
            tot += L1l[-cc[0]][1] * L2l[-cc[1]][1] * a ** (sl1[0] + sl1[1] + e) * ap ** (sl2[0] + sl2[1] + ep) * kv.subs({sS: s_, tt: t_}) * minor(L1l, *sl1) * minor(L2l, *sl2)
        M1 = minor(L1l, -3, -2) * L1l[-1][1]
        M2 = minor(L2l, -3, -2) * L2l[-1][1]
        lhs = (gaussCov(a, 0, a) * col)[0] * (gaussCov(ap, 0, ap) * col)[0] / x1 ** 2
        ok = sp.simplify(lhs * M1 * M2 - tot) == 0
        ok_all = ok_all and ok
        say("  (a,a') = (%s,%s):  <Phi_c|Psi_c>/<x0_c|Psi_c> M1 M2 = %.12g ; sum_{cc'} Q3 Q3 sum_k Bk2 Bk2 = %.12g ; exact equality: %s"
            % (a, ap, float(lhs * M1 * M2), float(tot), ok))
        # kappa check: <Omega_c(a)|P^-1 e1> = 1/kappa(z, w1 = 1/a)
        w = 1 / a
        kap = w ** 3 / ((zv[1] - w) * (zv[2] - w) ** 2 * (zv[0] - zv[1]))
        say("      kappa(z,w=(1/a,inf)) closed form (L1) * <Omega_c(a)|P^-1|e1> = %s (must be 1)" % ((gaussCov(a, 0, a) * col)[0] * kap))
    say("RESULT L0b: the ten-product / kernel-law formula, ladder points, minor/Bk2 conventions, (e,e'), M_vac, the companion ket Psi_c = (P^-1 (x) P^-1) psi_d and the "
        "x0 normalisation reproduce the exact overlap on the top state at finite h: %s" % ok_all)
    return dict(q2=q2s, zv=zv)


# ======================================================================================================================
# L1: ingredients
# ======================================================================================================================
def run_L1():
    hdr("L1  every h-dependent ingredient")
    z1, z2, z3, w1, h = sp.symbols("z1 z2 z3 w1 h")
    ze1, ze2, ze3, om = sp.symbols("zeta1 zeta2 zeta3 omega")

    def e3(i, j):
        M = sp.zeros(3)
        M[i - 1, j - 1] = 1
        return M

    def expN(M):
        return sp.eye(3) + M + M * M / 2

    sVec = [z1 ** 2 / ((z1 - z2) * (z1 - z3)), -(z2 / (z1 * (z2 - z3))), 1 / (z1 * z2)]
    rVec = [-(z3 / z2), -(z3 / z1), (z2 ** 2 * (z3 - z1)) / (z1 ** 2 * (z2 - z3))]
    P3 = (expN(rVec[0] * e3(3, 2)) * expN(rVec[1] * e3(3, 1)) * expN(rVec[2] * e3(2, 1)) * sp.diag(*sVec)
          * expN(e3(1, 2) / z1) * expN(-e3(1, 3) / (z1 * z2)) * expN((z1 + z2) / (z1 * z2) * e3(2, 3)))
    r1 = (w1 - z1) / (w1 - z2)
    r2 = -((w1 - z1) * (z2 - z3)) / ((z1 - z2) * (w1 - z3))          # w2 = infinity limit of the notebook rGauss
    r3 = ((w1 - z2) * (z1 - z3)) / ((z1 - z2) * (w1 - z3))
    n = sp.Matrix([[1, r1, r2 + r1 * r3], [0, 1, r3], [0, 0, 1]])
    g = (n * P3).applyfunc(sp.cancel)
    D1, D2 = sp.factor(g[0, 0]), sp.factor(g[:2, :2].det())
    say("(a) kappaTie(z,w) = Delta_1(n P3)^(lam1-lam2) Delta_2(n P3)^(lam2-lam3) Delta_3^lam3, Lambda = (2,1,0), w = (w1, infinity), n = nOf, P3 = Pphys3 (cell 252):")
    say("    Delta_1 =", D1, "   Delta_2 =", D2)
    kap = sp.factor(D1 * D2)
    say("    kappa(z,w) =", kap)
    kapz = kap.subs({z1: sp.exp(h * ze1), z2: sp.exp(h * ze2), z3: sp.exp(h * ze3), w1: sp.exp(h * om)})
    lead = sp.simplify(sp.limit(h ** 4 * kapz, h, 0))
    say("    z = exp(h zeta), w1 = exp(h omega):  kappa = h^-4 * [1/((zeta2-omega)(zeta3-omega)^2 (zeta1-zeta2))] (1 + O(h));  limit h^4 kappa =", sp.factor(lead))
    assert sp.simplify(lead - 1 / ((ze2 - om) * (ze3 - om) ** 2 * (ze1 - ze2))) == 0
    say("    => kappa kappa' = h^-8 / c(zeta,omega,omega'),  c := (zeta2-omega)(zeta3-omega)^2 (zeta2-omega')(zeta3-omega')^2 (zeta1-zeta2)^2")
    # (b) x0
    rep, E1 = gl3_generators((2, 1, 0))
    R = sp.Rational

    def tosp(op):
        M = sp.zeros(8, 8)
        for i, row in op.r.items():
            for j, v in row.items():
                M[i, j] = R(v.numerator, v.denominator)
        return M

    Ee = {k: tosp(v) for k, v in E1.items()}
    I8 = sp.eye(8)

    def nexp8(A):
        return I8 + sum((A ** k / sp.factorial(k) for k in range(1, 5)), sp.zeros(8, 8))

    wts = [[Ee[(k, k)][p, p] for k in (1, 2, 3)] for p in range(8)]
    ppinv = nexp8(-(z1 + z2) / (z1 * z2) * Ee[(2, 3)]) * nexp8(Ee[(1, 3)] / (z1 * z2)) * nexp8(-Ee[(1, 2)] / z1)
    dinv = sp.diag(*[sVec[0] ** (-w[0]) * sVec[1] ** (-w[1]) * sVec[2] ** (-w[2]) for w in wts])
    pminv = nexp8(-rVec[2] * Ee[(2, 1)]) * nexp8(-rVec[1] * Ee[(3, 1)]) * nexp8(-rVec[0] * Ee[(3, 2)])
    e1 = sp.Matrix([[1] + [0] * 7])
    xv = (e1 * ppinv * dinv * pminv).applyfunc(sp.cancel)
    say("(b) one-site x0(z) = e1 P_phys(z)^-1 (notebook order: GT patterns descending, index 0 = highest weight), EXACT components (all linear in z):")
    for i in range(8):
        say("      x0[%d] = %s" % (i, sp.factor(xv[i])))
    lin = [sp.expand(xv[i]).subs({z1: ze1, z2: ze2, z3: ze3}) for i in range(8)]
    say("    => x0(e^{h zeta}) = h x^G(zeta) + O(h^2) with x^G(zeta) = (%s) -- coefficient c0 = 1 exactly;  <x0^d|psi_d> = h^2 <x0^G|psi> (1+O(h)), x0^G = x^G (x) x^G" % ", ".join(str(sp.factor(l)) for l in lin))
    return dict(kappa=kap, xG=lin)


def run_L1_table():
    hdr("L1b  h-orders of every factor of RHS, kappa kappa' x0 x LHS (exponent e: factor ~ h^e)")
    mp.mp.dps = 60
    zeta = [mp.mpf(3) / 2, mp.mpf(5) / 4, mp.mpf(2) / 3]
    w1, w1p = mp.mpf(5) / 2, mp.mpf(7) / 4
    th12 = mp.mpf(31) / 35
    kf = {k: sp.lambdify((sS, tt), v, "mpmath") for k, v in KD.items()}

    def Kmat(h):
        a = mp.exp(-h * w1)
        ap = mp.exp(-h * w1p)
        s, t = ap / a, th12 / h
        M = mp.zeros(10, 10)
        for (cc, s1, s2), f in kf.items():
            e, ep = EE[cc]
            M[cc[0] * 5 + SLOTS.index(s1), cc[1] * 5 + SLOTS.index(s2)] = f(s, t) * a ** (s1[0] + s1[1] + e) * ap ** (s2[0] + s2[1] + ep)
        return M

    sv = {}
    for k in (3, 4):
        h = mp.mpf(10) ** -k
        sv[k] = mp.svd_r(Kmat(h), compute_uv=False)
    pw = [round(float(mp.log(sv[4][i] / sv[3][i]) / mp.log(mp.mpf(10) ** -1))) for i in range(6)]
    say("kernel matrix K(h) (10x10: (c,slot) x (c',slot'), exponentials a^(m+m'+e) a'^(n+n'+e') included): singular values at h=1e-4:",
        [mp.nstr(sv[4][i], 3) for i in range(7)], "; orders", pw, "(rank 6; 4 null directions)")
    say("  leading kernel (s=1, t=inf): K_0 = U (x) U with U = (u_0 | u_1), u_c = (c, 1-3c, -2+5c) on slots [(-3,-2),(-3,-1),(-2,-1)] (c = 0,1):")
    U = {0: (0, 1, -2), 1: (1, -2, 3)}
    # verify by limit of the exact laws
    ok = True
    for (cc, s1, s2), v in KD.items():
        lim = sp.limit(sp.limit(v.subs(sS, 1), tt, sp.oo), tt, sp.oo) if False else None
    h = mp.mpf(10) ** -9
    K9 = Kmat(h)
    err = 0
    for c in (0, 1):
        for cp in (0, 1):
            for i, s1 in enumerate(SLOTS[:3]):
                for j, s2 in enumerate(SLOTS[:3]):
                    err = max(err, abs(K9[c * 5 + i, cp * 5 + j] - U[c][i] * U[cp][j]))
    say("      max |K(h=1e-9) - U(x)U| over the 6x6 leading slots: %s ; slots (-3,0),(-2,0) are O(h) in K" % mp.nstr(err, 3))
    say("")
    say("%-58s %-8s %s" % ("factor", "order", "how"))
    rows = [
        ("kappa(z,w) = w1^3/((z2-w1)(z3-w1)^2(z1-z2))", "h^-4", "closed form (L1a); coefficient 1/((zeta2-omega)(zeta3-omega)^2(zeta1-zeta2))"),
        ("kappa kappa'", "h^-8", "product"),
        ("<x0^d|psi_d> = (x(z) (x) x(z)).psi", "h^2", "x(z) linear in z, x(e^{h zeta}) = h x^G + O(h^2) (L1b)"),
        ("<Phi_d|psi_d>, Phi_d = Omega_D[z,(w1,inf)] (x) Omega_D[z,(w1',inf)]", "h^0", "Gauss parameters r_k(z,w) = r_k(zeta,omega)+O(h) (ratios of O(h) quantities)"),
        ("Qnum_XXX = <Phi_d|psi_d>/(kappa kappa' <x0^d|psi_d>)", "h^6", "h^8 / h^2"),
        ("a^m, a'^n, s = a'/a, 1/t = h/theta12", "1, 1, 1+O(h), O(h)", "a = exp(-h omega) unit; (s-1) = O(h); t = theta12/h large"),
        ("kernel entries kappa^{cc'}_{(mm'),(nn')}(s,t)", "h^0", "polynomial in s, rational in t of degree 0: all 40 laws are O(1); singular values of K(h): 1, h^2, h^2, h^4, h^4, h^4"),
        ("wedge-vector entries f,g,f',g' (individually)", "h^-2..h^2", "1/(s-1) = O(1/h) in products 2, 7 (f'), t = theta12/h in product 2; all cancel in the kernel laws (L0 table)"),
        ("Q^j(theta_s+mh)  (jets J^j_{n,k} of P^j_n = e^{zeta_j u} q_n)", "h^0", "sum_{n,k} h^{n+k} m^k J_{n,k}/k!"),
        ("lone factors Q3(theta_s - c h), c = 0,1", "h^0", ""),
        ("minor_s(m,m') for m<m'", "h^2", "Wronskian zero: O(h^1) term = -(m'-m)W(theta_s) = 0; leading = (m-m')(m+m'+1) W'(theta_s)/2 (L2)"),
        ("M_vac = minor_s(-3,-2) Q3(theta_s-h)", "h^2 (per site)", "= 2 h^2 W'(theta_s) P3(theta_s)(1+O(h))   [NOT (5/2) h^2 W': see L2]"),
        ("M1 M2", "h^4", ""),
        ("naive (generic) order of RHS = sum_{cc'} Q3 Q3 sum A minor minor", "h^4", "kernel O(1) x minors h^2 x h^2"),
        ("exact leading order of RHS = Qnum_XXX M1 M2", "h^10", "h^6 x h^4; confirmed by the formal expansion (L2): cancellations through h^9"),
    ]
    for r in rows:
        say("%-58s %-8s %s" % r)
    return sv


def run_L1c():
    hdr("L1c  independent check of kappa kappa' and x0 on the top state: Qnum_XXX(h) = <Phi_c|Psi_c>/<x0_c|Psi_c> = F(h) -> h^6 c(zeta,omega,omega') Qnum_G")
    mp.mp.dps = 60
    zeta = [mp.mpf(3) / 2, mp.mpf(5) / 4, mp.mpf(2) / 3]
    th = [mp.mpf(2) / 7, -mp.mpf(3) / 5]
    w1, w1p = mp.mpf(5) / 2, mp.mpf(7) / 4
    Lam = (2, 1, 0)
    kf = {k: sp.lambdify((sS, tt), v, "mpmath") for k, v in KD.items()}

    def nuf(n, x, h):
        r = mp.mpf(1)
        for t_ in th:
            r *= (x - t_ + h * Lam[n - 1])
        return r

    def F_of(h):
        z = [mp.exp(h * zz) for zz in zeta]
        t1 = lambda u: sum(z[i] * nuf(i + 1, u, h) for i in range(3))
        t2 = lambda u: sum(z[i] * z[j] * nuf(i + 1, u, h) * nuf(j + 1, u - h, h) for i in range(3) for j in range(i + 1, 3))
        t3 = lambda u: z[0] * z[1] * z[2] * nuf(1, u, h) * nuf(2, u - h, h) * nuf(3, u - 2 * h, h)
        nu3 = lambda u: nuf(3, u, h)

        def rab(ths):
            rho = t1(ths) * nu3(ths - h) / t2(ths)
            den = t3(ths - h)
            al = (nu3(ths - h) * nu3(ths - 2 * h) * nu3(ths - 3 * h) - rho * t1(ths - h) * nu3(ths - 2 * h) * nu3(ths - 3 * h)) / den
            be = t2(ths - h) * nu3(ths - 3 * h) / den
            return rho, al, be

        r1, r2 = rab(th[0]), rab(th[1])
        a, ap = mp.exp(-h * w1), mp.exp(-h * w1p)
        s_, t_ = ap / a, (th[0] - th[1]) / h
        mino = lambda r: {(-3, -2): -r[1], (-3, -1): r[0] * r[2], (-3, 0): r[2], (-2, -1): r[0], (-2, 0): mp.mpf(1)}
        m1, m2 = mino(r1), mino(r2)
        tot = mp.mpf(0)
        for (cc, sl1, sl2), f in kf.items():
            e, ep = EE[cc]
            tot += r1[0] ** cc[0] * r2[0] ** cc[1] * a ** (sl1[0] + sl1[1] + e) * ap ** (sl2[0] + sl2[1] + ep) * f(s_, t_) * m1[sl1] * m2[sl2]
        return tot / ((-r1[1]) * r1[0] * (-r2[1]) * r2[0])

    z1, z2, z3 = zeta
    pred = (z2 - w1) * (z3 - w1) ** 2 * (z2 - w1p) * (z3 - w1p) ** 2          # c * Qnum_G(top), Qnum_G(top) = 1/(x^G_HW)^2 = 1/(zeta2-zeta1)^2
    vals = []
    for k in (3, 4, 5, 6):
        h = mp.mpf(10) ** -k
        vals.append(F_of(h) / h ** 6)
        say("  h = 1e-%d:  F(h)/h^6 = %s" % (k, mp.nstr(vals[-1], 14)))
    rich = (10 * vals[-1] - vals[-2]) / 9
    say("  Richardson (h=1e-5,1e-6): %s ; prediction c(zeta,omega,omega') Qnum_G(top) = (zeta2-w)(zeta3-w)^2 (zeta2-w')(zeta3-w')^2 = %s = 102245/41472"
        % (mp.nstr(rich, 14), mp.nstr(pred, 14)))
    say("  (F is built from the exact top-state t-eigenvalues through the Baxter ladder ratios (rho,alpha,beta), checked against the exact Q's in L0b)")


# ======================================================================================================================
# L2: formal limit
# ======================================================================================================================
def run_L2_ladder():
    hdr("L2a  the exact ladder constraint  minor_s(-1,0) = 0  (generalises W(theta_s) = 0)")
    x, h, th1, th2 = sp.symbols("x h theta1 theta2")
    z1, z2, z3 = sp.symbols("z1 z2 z3")
    Lam = (2, 1, 0)
    nuf = lambda n, y, th: (y - th[0] + h * Lam[n - 1]) * (y - th[1] + h * Lam[n - 1])
    th = (th1, th2)
    t3 = lambda y: z1 * z2 * z3 * nuf(1, y, th) * nuf(2, y - h, th) * nuf(3, y - 2 * h, th)
    nu3 = lambda y: nuf(3, y, th)
    say("dual Baxter: nu3(x)nu3(x-h)nu3(x-2h) Q(x+h) - t1(x) nu3(x-h)nu3(x-2h) Q(x) + t2(x) nu3(x-2h) Q(x-h) - t3(x) Q(x-2h) = 0,  nu3(x) = (x-theta1)(x-theta2), t3 = chi3 nu1(x) nu2(x-h) nu3(x-2h)")
    say("  at x = theta_1:  nu3(theta1) = %s,  t3(theta1) = %s   (nu2(x-h) = prod_s (x-theta_s))" % (sp.simplify(nu3(th1)), sp.simplify(t3(th1))))
    say("  => -t1 nu3(th-h) nu3(th-2h) Q(th) + t2 nu3(th-2h) Q(th-h) = 0:  Q^j(th_s - h) = rho_s Q^j(th_s), rho_s = t1(th_s) nu3(th_s-h)/t2(th_s) INDEPENDENT of j")
    say("  => minor_s(-1,0) = Q2(0)Q3(-1) - Q3(0)Q2(-1) = 0 EXACTLY (all orders in h, any XXX eigenstate; verified on the top state in L0b). With Q(-1) = rho Q(0) the ladder 4-vectors")
    say("     of Q2,Q3 lie in the 2-plane {v: v_-1 = rho v_0}; Cramer: v_-3 = alpha v_0 + beta v_-2; every minor_s(m,m') = (u_m' v_m - v_m' u_m) D, D = minor_s(-2,0):")
    say("     minor(-3,-2)/D = -alpha, minor(-3,-1)/D = rho beta, minor(-2,-1)/D = rho, minor(-3,0)/D = beta, minor(-2,0)/D = 1, minor(-1,0)/D = 0")
    say("     (leading values (alpha,beta,rho) = (-2,3,1); x=theta_s-h Baxter gives alpha, beta, rho through t1, t2 at theta_s, theta_s-h)")


def build_jets_ring(W=8):
    """symbolic generic jets with the ladder constraint imposed: J^j_{0,r} (r>=1) eliminated through Q(-1) = rho Q(0), rho = 1 + sum rho_i h^i"""
    h = sp.Symbol("h")
    rho = sp.symbols("rho1:%d" % (W + 1))
    J = {}
    for j in (2, 3):
        for n in range(W + 1):
            for k in range(W + 1 - n):
                J[(j, n, k)] = sp.Symbol("J%d_%d%d" % (j, n, k))
    sol = {}
    for j in (2, 3):
        for r in range(1, W + 1):
            rest = sum((-1) ** k * J[(j, r - k, k)] / sp.factorial(k) for k in range(1, r))
            rhs = sum(rho[i - 1] * J[(j, r - i, 0)] for i in range(1, r + 1))
            sol[J[(j, 0, r)]] = sp.expand(((-1) ** r * sp.factorial(r)) * (rhs - rest).subs(sol))
    return h, rho, J, sol


def run_L2_Z():
    hdr("L2b  site-local contraction Z_s = sum_c Q3(-c) u_c . minor(slots):  vanishes through weight 4 once minor(-1,0)=0 is imposed")
    W = 6
    h, rho, J, sol = build_jets_ring(W)

    def Qv(j, m):
        expr = sum(h ** (n + k) * m ** k * J[(j, n, k)] / sp.factorial(k) for n in range(W + 1) for k in range(W + 1 - n) if n + k <= W)
        return sp.expand(expr.subs(sol))

    def trunc(e, w):
        e = sp.expand(e)
        return sum(e.coeff(h, i) * h ** i for i in range(w + 1))

    Q = {(j, m): Qv(j, m) for j in (2, 3) for m in MPTS}
    mn = lambda m, mp_: trunc(Q[(2, mp_)] * Q[(3, m)] - Q[(3, mp_)] * Q[(2, m)], W)
    say("  minor(-1,0) with J^j_{0,r}, r=1..%d, eliminated through Q(-1)=rho Q(0) (rho_i free, j-independent):" % W, sp.simplify(mn(-1, 0)))
    M = {s: mn(*s) for s in SLOTS[:3]}
    u = {0: (0, 1, -2), 1: (1, -2, 3)}
    Z = sum(Q[(3, -c)] * sum(ui * M[s] for ui, s in zip(u[c], SLOTS[:3])) for c in (0, 1))
    Z = trunc(Z, W)
    for w in range(W + 1):
        cw = sp.expand(Z.coeff(h, w))
        say("  weight %d: %s" % (w, "0 (identically, in all jets and rho_i)" if cw == 0 else "%d terms%s" % (len(sp.Add.make_args(cw)), " (first nonzero weight)" if w == 5 else "")))
    say("  also generic weight-2 minors (no constraint): minor(m,m') = h^2 (m-m')[(m+m')p/2 + Q],  p = W'(theta_s) = J2_00 J3_02 - J3_00 J2_02,  Q = q11+q10 (q1-jets);")
    say("  minor(-1,0) = 0 forces Q = p/2:  minor(m,m') = h^2 (m-m')(m+m'+1) p/2:  (-3,-2): 2p, (-3,-1): 3p, (-2,-1): p, (-3,0): 3p, (-2,0): p, (-1,0): 0.")
    say("  => M_vac = 2 h^2 W'(theta_s) P3(theta_s) (1+O(h)); the spec's (5/2) h^2 W' used only the n=0 jets and is corrected by the n=1 jet Q = p/2.")


# ----------------------------------------------------------------------------------------------------------------------
# exact random-jet engine (Fractions, truncated series)
# ----------------------------------------------------------------------------------------------------------------------
NS = xxx_chain.NSER
x_ = sp.Symbol("x")


def kernel_series_data():
    out = {}
    for key, v in KD.items():
        num, den = sp.fraction(sp.together(v))
        pn = sp.Poly(num, sS, tt)
        pd = sp.Poly(den, tt)
        m = max(pn.degree(tt), pd.degree(tt))
        Nh = sp.Poly(sp.expand(sum(c * sS ** i * x_ ** (m - j) for (i, j), c in pn.terms())), sS, x_)
        Dh = sp.Poly(sp.expand(sum(c * x_ ** (m - j) for (j,), c in pd.terms())), x_)
        out[key] = (Nh, Dh)
    return out


KSD = kernel_series_data()


def kernel_series(key, delta, th12):
    """kappa(s,t) as Ser in h with s = exp(-h delta), t = th12/h"""
    Nh, Dh = KSD[key]
    sser = Ser.exp(-Fr(delta))
    xser = Ser([Fr(0), 1 / Fr(th12)])
    n = Ser.const(0)
    for (i, j), c in Nh.terms():
        n = n + (sser ** i) * (xser ** j) * Fr(int(c.p), int(c.q))
    d = Ser.const(0)
    for (j,), c in Dh.terms():
        d = d + (xser ** j) * Fr(int(c.p), int(c.q))
    return n / d


class JetEngine:
    """RHS of the XXX identity as an exact h-series for random jets obeying only the ladder constraint minor_s(-1,0)=0"""

    def __init__(self, w1, w1p, th12):
        self.w1, self.w1p, self.th12 = Fr(w1), Fr(w1p), Fr(th12)
        self.KER = {}
        for key in KD:
            cc, s1, s2 = key
            e, ep = EE[cc]
            self.KER[key] = (kernel_series(key, self.w1p - self.w1, self.th12) * Ser.exp(-self.w1 * (s1[0] + s1[1] + e))
                             * Ser.exp(-self.w1p * (s2[0] + s2[1] + ep)))
        self.W = NS - 1

    def randvars(self, rnd):
        W = self.W
        var = {}
        for r in range(1, W + 1):
            var[("rho", r)] = rnd()
        for j in (2, 3):
            for n in range(W + 1):
                for k in range(W + 1 - n):
                    if (n, k) == (0, 1) or (k == 1 and n >= 1):
                        continue
                    var[("J", j, n, k)] = rnd() if (n, k) != (0, 0) else abs(rnd()) + 1
        return var

    def build(self, var):
        W = self.W
        rho = [Fr(1)] + [var[("rho", r)] for r in range(1, W + 1)]
        J = {}
        for j in (2, 3):
            for n in range(W + 1):
                for k in range(W + 1 - n):
                    if (n, k) == (0, 1) or (k == 1 and n >= 1):
                        continue
                    J[(j, n, k)] = var[("J", j, n, k)]
            for r in range(1, W + 1):            # ladder relation at weight r solved for J_{r-1,1}
                rhs = sum(rho[i] * J[(j, r - i, 0)] for i in range(1, r + 1))
                rest = sum(Fr((-1) ** k, factorial(k)) * J[(j, r - k, k)] for k in range(2, r + 1) if (j, r - k, k) in J)
                J[(j, r - 1, 1)] = -(rhs - rest)
        Q = {}
        for j in (2, 3):
            for m in MPTS:
                c = [Fr(0)] * NS
                for n in range(NS):
                    for k in range(NS - n):
                        c[n + k] += Fr(m) ** k * J[(j, n, k)] / factorial(k)
                Q[(j, m)] = Ser(c)
        return Q, J

    @staticmethod
    def minor(Q, m, mp_):
        return Q[(2, mp_)] * Q[(3, m)] - Q[(3, mp_)] * Q[(2, m)]

    def rhs(self, Q1, Q2):
        V1 = {(c, s): Q1[(3, -c)] * self.minor(Q1, *s) for c in (0, 1) for s in SLOTS}
        V2 = {(c, s): Q2[(3, -c)] * self.minor(Q2, *s) for c in (0, 1) for s in SLOTS}
        tot = Ser.const(0)
        for (cc, s1, s2), ker in self.KER.items():
            tot = tot + ker * V1[(cc[0], s1)] * V2[(cc[1], s2)]
        return tot


def run_L2_random(Gs5):
    import random
    hdr("L2c  leading order from generic jets (exact Fractions): only minor_s(-1,0)=0 imposed, all other jets random")
    w1, w1p, th12 = Fr(5, 2), Fr(7, 4), Fr(31, 35)
    eng = JetEngine(w1, w1p, th12)
    random.seed(7)
    rnd = lambda: Fr(random.randint(-9, 9), random.randint(1, 5))
    Gn = [[Fr(int(sp.Rational(Gs5[i][j].subs({sp.Symbol("w1"): sp.Rational(5, 2), sp.Symbol("w1p"): sp.Rational(7, 4), sp.Symbol("thi"): sp.Rational(35, 31)})).p),
              int(sp.Rational(Gs5[i][j].subs({sp.Symbol("w1"): sp.Rational(5, 2), sp.Symbol("w1p"): sp.Rational(7, 4), sp.Symbol("thi"): sp.Rational(35, 31)})).q)) for j in range(5)] for i in range(5)]

    def Nvec(var):
        J = {}
        for j in (2, 3):
            for k in (0, 2, 3):
                J[(j, k)] = var[("J", j, 0, k)]
        kap = -var[("rho", 1)]
        for j in (2, 3):
            J[(j, 1)] = kap * J[(j, 0)]
        p = lambda k, l: J[(2, k)] * J[(3, l)] - J[(3, k)] * J[(2, l)]
        b0 = J[(3, 0)]
        return [b0 * p(0, 2), kap * b0 * p(0, 2), b0 * p(0, 3), kap * b0 * p(0, 3), b0 * p(2, 3)]

    for seed in range(3):
        random.seed(100 + seed)
        v1, v2 = eng.randvars(rnd), eng.randvars(rnd)
        Q1, _ = eng.build(v1)
        Q2, _ = eng.build(v2)
        assert eng.minor(Q1, -1, 0).val() is None and eng.minor(Q2, -1, 0).val() is None
        R = eng.rhs(Q1, Q2)
        N1, N2 = Nvec(v1), Nvec(v2)
        val = sum(Gn[i][j] * N1[i] * N2[j] for i in range(5) for j in range(5))
        say("  random jets #%d: val(minor(-3,-2)) = %d, valuation of RHS = %d (coefficients h^4..h^9 all zero: %s); coefficient of h^10 == sum_ij G_ij N_i N_j (exact): %s"
            % (seed, eng.minor(Q1, -3, -2).val(), R.val(), all(R.c[i] == 0 for i in range(10)), R.c[10] == val))
    # which jets enter at order 10 (perturbation of one jet at a time)
    random.seed(5)
    v1, v2 = eng.randvars(rnd), eng.randvars(rnd)
    Q1, _ = eng.build(v1)
    Q2, _ = eng.build(v2)
    base = eng.rhs(Q1, Q2)
    table = {}
    for key in sorted(v1.keys(), key=lambda k: (k[0], k[1:])):
        n_k = key[2:] if key[0] == "J" else None
        if key[0] == "J" and key[2] > 5:
            continue
        vv = dict(v1)
        vv[key] = vv[key] + Fr(1, 3)
        Qp, _ = eng.build(vv)
        d = eng.rhs(Qp, Q2) - base
        table[key] = d.val()
    say("  first order (valuation of the change of RHS) at which each site-1 variable enters when perturbed (None: beyond h^%d):" % (NS - 1))
    say("    free variables of the ladder-constrained jets: J^j_{n,k} (J^j_{0,1} = kappa J^j_{0,0}, J^j_{n>=1,1} fixed by the ladder relation), rho_r (j-independent, rho_1 = -kappa)")
    ent10 = sorted(set(("J%d_%d%d" % (k[1], k[2], k[3])) if k[0] == "J" else ("rho%d" % k[1]) for k, v in table.items() if v == 10))
    say("    variables that already enter at order 10:", ent10)
    nge1 = {(k[2], k[3]): v for k, v in table.items() if k[0] == "J" and k[1] == 3 and k[2] >= 1}
    say("    n >= 1 jets J3_{n,k} (k != 1) first enter at order:", {k: nge1[k] for k in sorted(nge1) if k[0] <= 4 and k[1] <= 4})
    say("    rho_r first enter at order:", {k[1]: v for k, v in table.items() if k[0] == "rho" and k[1] <= 6})
    say("    J^j_{0,k}, k>=4 first enter at order:", {k[3]: v for k, v in table.items() if k[0] == "J" and k[1] == 3 and k[2] == 0 and 4 <= k[3] <= 8})
    say("  => order-10 coefficient is a function of the n=0 jets J^j_{0,0}, J^j_{0,1}=kappa J^j_{0,0}, J^j_{0,2}, J^j_{0,3} ONLY (rho_1 = -kappa); q_1 (n=1), q_2 jets and rho_{r>=2} enter first at order 11 (shift n), J_{0,k>=4} at 7+k.")


def run_L2_symbolic():
    hdr("L2d  the order-10 coefficient: symbolic in (omega1, omega1', theta12) in the effective ladder model")
    names = ["h", "w1", "w1p", "thi"] + ["a0_1", "a2_1", "a3_1", "b0_1", "b2_1", "b3_1", "k_1"] + ["a0_2", "a2_2", "a3_2", "b0_2", "b2_2", "b3_2", "k_2"]
    R, *G = ring(",".join(names), QQ)
    h, w1, w1p, thi = G[:4]
    S1, S2 = G[4:11], G[11:18]
    NMAX = 10

    def trunc(p, M):
        return R.from_dict({m: c for m, c in p.items() if m[0] <= M})

    def expser(c, M=NMAX):
        out = R(1)
        term = R(1)
        for k in range(1, M + 1):
            term = trunc(term * (-h * c), M) * QQ(1, k)
            out = out + term
        return out

    def Qeff(j, m, S):
        a0, a2, a3, b0, b2, b3, k = S
        J0, J2, J3 = (a0, a2, a3) if j == 2 else (b0, b2, b3)
        return J0 * (1 + k * h * m) + J2 * QQ(1, 2) * h ** 2 * m * (m + 1) + J3 * QQ(1, 6) * h ** 3 * m * (m ** 2 - 1)

    def Vvec(S):
        Q = {(j, m): Qeff(j, m, S) for j in (2, 3) for m in MPTS}
        mn = lambda m, mp_: Q[(2, mp_)] * Q[(3, m)] - Q[(3, mp_)] * Q[(2, m)]
        return {(c, s): Q[(3, -c)] * mn(*s) for c in (0, 1) for s in SLOTS}, mn

    V1, mn1 = Vvec(S1)
    V2, mn2 = Vvec(S2)
    say("effective local ladder function (all that enters at order 10; J^j_{0,1} = kappa J^j_{0,0}, J_{1,1} = J_{0,2}/2, J_{2,1} = -J_{0,3}/6 from the ladder relation):")
    say("    Q^j(theta_s + m h) = J^j_0 (1 + kappa m h) + (J^j_2/2) h^2 m(m+1) + (J^j_3/6) h^3 m(m^2-1)      (all other jets set to 0: they do not enter at order 10, L2c)")
    say("  minor(-1,0) == 0 identically in the effective model:", mn1(-1, 0) == 0, mn2(-1, 0) == 0)
    x = sp.Symbol("x")
    delta = w1p - w1

    def kernel_ser(key, M):
        v = KD[key]
        cc, s1, s2 = key
        e, ep = EE[cc]
        Nh, Dh = KSD[key]
        xr = h * thi
        N = R(0)
        spow = {}
        for (i, j), c in Nh.terms():
            if i not in spow:
                spow[i] = expser(i * delta, M)
            N = N + trunc(spow[i] * (xr ** j), M) * QQ(int(c.p), int(c.q))
        dco = [(j, QQ(int(c.p), int(c.q))) for (j,), c in Dh.terms()]
        D = R(0)
        d0 = QQ(0)
        for j, c in dco:
            D = D + (xr ** j) * c
            if j == 0:
                d0 = c
        rest = D - R(d0)
        inv = R(1)
        term = R(1)
        for k in range(1, M + 1):
            term = trunc(term * (-rest) * (1 / d0), M)
            inv = inv + term
        inv = inv * (1 / d0)
        ker = trunc(N * inv, M)
        ker = trunc(ker * expser(w1 * (s1[0] + s1[1] + e), M), M)
        ker = trunc(ker * expser(w1p * (s2[0] + s2[1] + ep), M), M)
        return ker

    KER = {key: kernel_ser(key, NMAX - 4) for key in KD}

    def hparts(p):
        out = {}
        for m, c in p.items():
            out.setdefault(m[0], {})[(0,) + m[1:]] = c
        return {k: R.from_dict(v) for k, v in out.items()}

    V1p = {k: hparts(v) for k, v in V1.items()}
    V2p = {k: hparts(v) for k, v in V2.items()}
    Kp = {k: hparts(v) for k, v in KER.items()}

    def coefficient(N):
        tot = R(0)
        for (cc, s1, s2), kp in Kp.items():
            A = V1p[(cc[0], s1)]
            B = V2p[(cc[1], s2)]
            for wa, pa in A.items():
                for wb, pb in B.items():
                    m = N - wa - wb
                    if m < 0 or m not in kp:
                        continue
                    tot = tot + pa * pb * kp[m]
        return tot

    coeffs = {N: coefficient(N) for N in range(4, 11)}
    say("  coefficient of h^N of the full RHS, symbolic in all jets, omega1, omega1', 1/theta12: ", {N: ("0" if c == 0 else "%d terms" % len(c.terms())) for N, c in coeffs.items()})
    say("  => N = 10: the RHS vanishes identically at orders h^4 ... h^9 (cancellations between the ten products / four blocks / jets)")
    e10 = sp.sympify(str(coeffs[10].as_expr()))
    # convert to Pluecker variables
    s1v = sp.symbols("a0_1 a2_1 a3_1 b0_1 b2_1 b3_1 k_1")
    s2v = sp.symbols("a0_2 a2_2 a3_2 b0_2 b2_2 b3_2 k_2")
    sub = {}
    for tag in ("1", "2"):
        a0, a2, a3, b0, b2, b3, k = sp.symbols("a0_%s a2_%s a3_%s b0_%s b2_%s b3_%s k_%s" % ((tag,) * 7))
        p02, p03 = sp.symbols("p02_%s p03_%s" % (tag, tag))
        sub[tag] = {a2: (a0 * b2 - p02) / b0, a3: (a0 * b3 - p03) / b0}
    e2 = sp.expand(e10.subs(sub["1"]).subs(sub["2"]))
    a01, a02 = sp.symbols("a0_1 a0_2")
    say("  depends on a0 (i.e. on anything but the Pluecker coordinates p_kl = J2_k J3_l - J3_k J2_l, J3_k, kappa)?", sp.simplify(sp.diff(e2, a01)) != 0, sp.simplify(sp.diff(e2, a02)) != 0)
    v1 = sp.symbols("b0_1 b2_1 b3_1 p02_1 p03_1 k_1")
    v2 = sp.symbols("b0_2 b2_2 b3_2 p02_2 p03_2 k_2")
    P = sp.Poly(e2, *v1, *v2)
    C = defaultdict(dict)
    for m, c in P.terms():
        C[m[:6]][m[6:]] = c
    B0, B2, B3, P02, P03, K = v1
    C0, C2, C3, Q02, Q03, KK = v2
    mons = {"b0p02": B0 * P02, "kb0p02": K * B0 * P02, "b0p03": B0 * P03, "kb0p03": K * B0 * P03, "b2p03": B2 * P03, "b3p02": B3 * P02}
    mons2 = {"b0p02": C0 * Q02, "kb0p02": KK * C0 * Q02, "b0p03": C0 * Q03, "kb0p03": KK * C0 * Q03, "b2p03": C2 * Q03, "b3p02": C3 * Q02}
    nm = list(mons)
    Gm = sp.zeros(6, 6)
    for m1, d in C.items():
        mono1 = sp.Mul(*[s ** p for s, p in zip(v1, m1)])
        i = [n for n, b in mons.items() if sp.expand(b - mono1) == 0][0]
        for m2, val in d.items():
            mono2 = sp.Mul(*[s ** p for s, p in zip(v2, m2)])
            j = [n for n, b in mons2.items() if sp.expand(b - mono2) == 0][0]
            Gm[nm.index(i), nm.index(j)] = sp.factor(val)
    say("  order-10 coefficient as a bilinear form in the 6 site monomials (the b2p03 and b3p02 columns are exactly opposite:",
        all(sp.simplify(Gm[i, 4] + Gm[i, 5]) == 0 for i in range(6)) and all(sp.simplify(Gm[4, j] + Gm[5, j]) == 0 for j in range(6)), ")")
    say("  => with the Pluecker identity b0 p23 = b2 p03 - b3 p02 it is a 5x5 form in N = (P3 p02, P3' p02, P3 p03, P3' p03, P3 p23)")
    G5 = Gm[:5, :5]
    names5 = ["P3 p02", "P3' p02", "P3 p03", "P3' p03", "P3 p23"]
    w1s, w1ps, this = sp.symbols("w1 w1p thi")
    sym_ok = sp.simplify(G5.T.subs({w1s: w1ps, w1ps: w1s, this: -this}, simultaneous=True) - G5) == sp.zeros(5, 5)
    wt = [2, 3, 3, 4, 5]
    lam = sp.Symbol("lam")
    hom_ok = all(sp.simplify(sp.expand(G5[i, j].subs({w1s: lam * w1s, w1ps: lam * w1ps, this: lam * this}, simultaneous=True)) - lam ** (10 - wt[i] - wt[j]) * G5[i, j]) == 0
                 for i in range(5) for j in range(5))
    say("  rank of G5:", G5.rank(), "; site-exchange symmetry (omega<->omega', theta12 -> -theta12, transpose):", sym_ok, "; homogeneity deg 10 - w_i - w_j (weights 2,3,3,4,5):", hom_ok)
    # ---- twisted jets: jets of e^{-omega_s u} P^j at theta_s  (a^m Q(theta+mh) = twisted ladder function)
    def site_twist(tag, om):
        a0, a2, a3, b0, b2, b3, k = sp.symbols("a0_%s a2_%s a3_%s b0_%s b2_%s b3_%s k_%s" % ((tag,) * 7))
        A0, A2, A3, B0, B2, B3, Kt = sp.symbols("A0_%s A2_%s A3_%s B0_%s B2_%s B3_%s K_%s" % ((tag,) * 7))
        return {a0: A0, b0: B0, k: Kt + om, a2: A2 + (2 * om * Kt + om ** 2) * A0, b2: B2 + (2 * om * Kt + om ** 2) * B0,
                a3: A3 + 3 * om * A2 + 3 * om ** 2 * Kt * A0 + om ** 3 * A0, b3: B3 + 3 * om * B2 + 3 * om ** 2 * Kt * B0 + om ** 3 * B0}
    wA, wB = sp.symbols("w1 w1p")
    et = sp.expand(e10.subs(site_twist("1", wA)).subs(site_twist("2", wB)))
    for tag in ("1", "2"):
        A0, A2, A3, B0, B2, B3, Kt = sp.symbols("A0_%s A2_%s A3_%s B0_%s B2_%s B3_%s K_%s" % ((tag,) * 7))
        P02t, P03t = sp.symbols("P02_%s P03_%s" % (tag, tag))
        et = et.subs({A2: (A0 * B2 - P02t) / B0, A3: (A0 * B3 - P03t) / B0})
    et = sp.expand(et)
    t1v = sp.symbols("B0_1 B2_1 B3_1 P02_1 P03_1 K_1")
    t2v = sp.symbols("B0_2 B2_2 B3_2 P02_2 P03_2 K_2")
    Pt = sp.Poly(et, *t1v, *t2v)
    Ct = defaultdict(dict)
    for m, c in Pt.terms():
        Ct[m[:6]][m[6:]] = c
    TB0, TB2, TB3, TP02, TP03, TK = t1v
    UB0, UB2, UB3, UP02, UP03, UK = t2v
    tm1 = {"b0p02": TB0 * TP02, "kb0p02": TK * TB0 * TP02, "b0p03": TB0 * TP03, "kb0p03": TK * TB0 * TP03, "b2p03": TB2 * TP03, "b3p02": TB3 * TP02}
    tm2 = {"b0p02": UB0 * UP02, "kb0p02": UK * UB0 * UP02, "b0p03": UB0 * UP03, "kb0p03": UK * UB0 * UP03, "b2p03": UB2 * UP03, "b3p02": UB3 * UP02}
    Gt = sp.zeros(6, 6)
    for m1, d in Ct.items():
        mono1 = sp.Mul(*[s ** p for s, p in zip(t1v, m1)])
        i = [n for n, b in tm1.items() if sp.expand(b - mono1) == 0][0]
        for m2, val in d.items():
            mono2 = sp.Mul(*[s ** p for s, p in zip(t2v, m2)])
            j = [n for n, b in tm2.items() if sp.expand(b - mono2) == 0][0]
            Gt[nm.index(i), nm.index(j)] = val
    dl = sp.Symbol("delta")
    G5t = Gt[:5, :5].applyfunc(lambda e: sp.factor(e.subs(wB, wA + dl)))
    indep = all(sp.simplify(sp.diff(G5t[i, j], wA)) == 0 for i in range(5) for j in range(5))
    say("  twisted jets (J~_k = (d/du - omega_s)^k P at theta_s, kappa~ = kappa - omega_s): the kernel G~ depends on delta = omega' - omega and thi = 1/theta12 ONLY:", indep)
    return G5, names5, G5t


def run_L4():
    hdr("L4  naive companion-frame limit: a = 1/omega1, a' = 1/omega1' FIXED, t -> infinity, ladder brackets -> jet limits")
    w, wp = sp.symbols("w wp")
    lim = {}
    for key, v in KD.items():
        num, den = sp.fraction(sp.together(v))
        pn, pd = sp.Poly(num, tt), sp.Poly(den, tt)
        lim[key] = sp.factor(sp.simplify(pn.coeff_monomial(tt ** pd.degree()) / pd.LC())) if pn.degree() == pd.degree() else 0
    nz = {k: v for k, v in lim.items() if v != 0}
    U = {0: (0, 1, -2), 1: (1, -2, 3)}
    ok = all(sp.simplify(v - sS * U[k[0][0]][SLOTS.index(k[1])] * U[k[0][1]][SLOTS.index(k[2])]) == 0 for k, v in nz.items()) and len(nz) == 25
    say("lim_{t->inf} kappa^{cc'}_{(mm'),(nn')}(s,t): 25 nonzero laws (slots (-3,-2),(-3,-1),(-2,-1) only), all equal to s * u_c[mm'] u_c'[nn']: %s  (u_0 = (0,1,-2), u_1 = (1,-2,3); the same U (x) U as K(0), times s)" % ok)
    M = lambda m, mp_: sp.Rational((m - mp_) * (m + mp_ + 1), 2)
    tot = 0
    for (cc, s1, s2), v in nz.items():
        e, ep = EE[cc]
        tot += (1 / w) ** (s1[0] + s1[1] + e) * (1 / wp) ** (s2[0] + s2[1] + ep) * v.subs(sS, w / wp) * M(*s1) * M(*s2)
    say("naive companion-frame RHS/(M1_vac M2_vac) with minors -> h^2 (m-m')(m+m'+1) p/2:", sp.factor(sp.simplify(tot / 4)), " = (1-a)^3 (1-a')^3, a = 1/omega, STATE-INDEPENDENT")
    say("Gaudin companion tie: kappa^G(zeta,omega) = omega^3/((zeta2-omega)(zeta3-omega)^2(zeta1-zeta2)) (W2's Delta_1 Delta_2 at omega2 -> inf) = XXX kappa with z->zeta, w->omega;")
    say("  hence Qnum_G/(kappa^G kappa^G') = c Qnum_G/(omega^3 omega'^3), which VARIES with the state (L3), whereas the naive limit above is a state-independent constant:")
    say("  the two do NOT agree. The fixed-a limit is a different limit (a = e^{-h omega} -> 1 in the scaling family): (1-a)^3(1-a')^3 ~ h^6 omega^3 omega'^3 reproduces only the h^6 power of Qnum_XXX,")
    say("  all state dependence sits in the O(h^4)-relative cancellations of the diagonal-frame analysis (L2): the naive kernel s U(x)U is rank one and cannot produce it.")


# ======================================================================================================================
# L3: numerics on the 64 Gaudin states
# ======================================================================================================================
def run_L3(G5, G5t):
    hdr("L3  numerical test on the 64 Gaudin states of the sl(3) chain, Lambda=(2,1,0), L=2   (60 digits)")
    import qside
    from qside import Qfun, mpf, q1_run, Operators, eigenstates, eigenfunctions, Chain, exp_nil
    from gl3_reps import GL3Rep
    mp.mp.dps = 60
    rep = GL3Rep((2, 1, 0))
    order = sorted(rep.patterns, reverse=True)
    perm = [rep._index[p] for p in order]
    say("basis map: notebook index a (descending GT patterns, 0 = HW) -> gl3_reps index perm[a] =", perm, "(HW is index 7 in gl3_reps order)")
    w1s, w1ps, this = sp.symbols("w1 w1p thi")
    dls = sp.Symbol("delta")
    from math import comb
    zeta_saved = qside.ZETA

    def run_case(zeta_fr, theta_fr, oms, label, Q1check=False):
        t0 = time.time()
        qside.ZETA = [Fr(x) for x in zeta_fr]
        ch = Chain(theta=tuple(Fr(x) for x in theta_fr))
        g = [[qside.ZETA[i] if i == j else Fr(0) for j in range(3)] for i in range(3)]
        ops = Operators(ch, g)
        C = ops.fit()
        st = eigenstates(ch, ops, C, dps=60)
        res = eigenfunctions(ch, C, st)
        th = [mpf(t) for t in ch.theta]
        recs = q1_run(ch, st, th)
        zeta = [mpf(z) for z in qside.ZETA]
        z1, z2, z3 = zeta
        xn = [z2 - z1, z3 - z1, z2 - z1, 2 * z3 - z1 - z2, 3 * (z2 - z1), 3 * (z3 - z1), 2 * (z3 - z2), 3 * (z3 - z2)]
        xg = [0] * 8
        for a in range(8):
            xg[perm[a]] = xn[a]
        x0G = [a * b for a in xg for b in xg]
        sigma = z1 + z2 + z3
        say("[%s] zeta = %s, theta = %s: %d states, max eigen-residual %s, %.0fs" % (label, [str(x) for x in zeta_fr], [str(x) for x in theta_fr], len(st), mp.nstr(res, 3), time.time() - t0))
        for (o1, o2) in oms:
            def omsite(w):
                zz = qside.ZETA
                w = Fr(w)
                r1 = (w - zz[0]) / (w - zz[1])
                r2 = -(w - zz[0]) * (zz[1] - zz[2]) / ((zz[0] - zz[1]) * (w - zz[2]))
                r3 = (w - zz[1]) * (zz[0] - zz[2]) / ((zz[0] - zz[1]) * (w - zz[2]))
                e = [Fr(0)] * ch.d1
                e[ch.hw_index] = Fr(1)
                v = exp_nil(ch.E1[(1, 2)], r1).rapply(e)
                v = exp_nil(ch.E1[(1, 3)], r2).rapply(v)
                v = exp_nil(ch.E1[(2, 3)], r3).rapply(v)
                return [mpf(a) for a in v]

            Phi = [a * b for a in omsite(o1) for b in omsite(o2)]
            th12 = Fr(theta_fr[0]) - Fr(theta_fr[1])
            sub = {w1s: sp.Rational(o1.numerator, o1.denominator), w1ps: sp.Rational(o2.numerator, o2.denominator), this: sp.Rational(th12.denominator, th12.numerator)}
            G5n = [[(lambda q: mp.mpf(q.p) / mp.mpf(q.q))(sp.Rational(G5[i, j].subs(sub))) for j in range(5)] for i in range(5)]
            subt = {dls: sp.Rational((o2 - o1).numerator, (o2 - o1).denominator), this: sp.Rational(th12.denominator, th12.numerator)}
            G5tn = [[(lambda q: mp.mpf(q.p) / mp.mpf(q.q))(sp.Rational(G5t[i, j].subs(subt))) for j in range(5)] for i in range(5)]
            om1, om1p = mpf(o1), mpf(o2)
            omS = (om1, om1p)
            devtw = mp.mpf(0)
            cc = (z2 - om1) * (z3 - om1) ** 2 * (z2 - om1p) * (z3 - om1p) ** 2 * (z1 - z2) ** 2
            ratios, rat2, frob, wr, dq = [], [], mp.mpf(0), mp.mpf(0), mp.mpf(0)
            for s, rec in zip(st, recs):
                psi = s["vec"]
                qn = mp.fsum(a * b for a, b in zip(Phi, psi)) / mp.fsum(a * b for a, b in zip(x0G, psi))
                P = [Qfun(zeta[j], rec["dual"][j][3]) for j in range(3)]
                Qd1 = Qfun(-zeta[0], rec["direct"][0][3])
                Ns, nus, Ntw = [], [], []
                for ss in (0, 1):
                    u = th[ss]
                    J2 = [P[1](u, k) for k in range(4)]
                    J3 = [P[2](u, k) for k in range(4)]
                    kap = J3[1] / J3[0]
                    frob = max(frob, abs(J2[1] / J2[0] - kap))
                    p = lambda k, l: J2[k] * J3[l] - J3[k] * J2[l]
                    wr = max(wr, abs(p(0, 1)))
                    N = [J3[0] * p(0, 2), kap * J3[0] * p(0, 2), J3[0] * p(0, 3), kap * J3[0] * p(0, 3), J3[0] * p(2, 3)]
                    Ns.append((N, J3[0] * p(0, 2)))
                    # twisted jets (e^{-omega_s u} P^j)
                    tw = lambda J: [mp.fsum(comb(k, jj) * (-omS[ss]) ** (k - jj) * J[jj] for jj in range(k + 1)) for k in range(4)]
                    T2_, T3_ = tw(J2), tw(J3)
                    kapt = T3_[1] / T3_[0]
                    pt = lambda k, l: T2_[k] * T3_[l] - T3_[k] * T2_[l]
                    Nt = [T3_[0] * pt(0, 2), kapt * T3_[0] * pt(0, 2), T3_[0] * pt(0, 3), kapt * T3_[0] * pt(0, 3), T3_[0] * pt(2, 3)]
                    Ntw.append(Nt)
                    # direct-basis check: W'(theta_s) = (zeta3-zeta2) e^{sigma theta_s} Pi'(theta_s) Q1(theta_s),  Pi' = theta_s - theta_s'
                    Pip = th[ss] - th[1 - ss]
                    A = (z3 - z2) * mp.exp(sigma * u) * Pip
                    Q1v, Q1p = Qd1(u, 0), Qd1(u, 1)
                    dq = max(dq, abs(p(0, 2) - A * Q1v) / abs(p(0, 2)))
                    Lam = 2 * sigma + 2 / Pip
                    pi_ = (Lam - kap) + 2 * Q1p / Q1v                      # p03/p02
                    eps_ = (J3[2] / J3[0]) * pi_ - J3[3] / J3[0]         # p23/p02 = (b2 p03 - b3 p02)/(b0 p02)
                    nus.append([mp.mpf(1), kap, pi_, kap * pi_, eps_])
                R10 = mp.fsum(G5n[i][j] * Ns[0][0][i] * Ns[1][0][j] for i in range(5) for j in range(5))
                R10t = mp.fsum(G5tn[i][j] * Ntw[0][i] * Ntw[1][j] for i in range(5) for j in range(5))
                devtw = max(devtw, abs(R10t / R10 - 1))
                ratios.append(R10 / (cc * qn * 4 * Ns[0][1] * Ns[1][1]))
                Qformula = mp.fsum(G5n[i][j] * nus[0][i] * nus[1][j] for i in range(5) for j in range(5)) / (4 * cc)
                rat2.append(Qformula / qn)
            say("   (omega,omega') = (%s,%s): R10/(c Qnum 4 (P3 W')_1 (P3 W')_2) = %s on all 64 states, max rel. deviation %s;  direct-Q1 form: Qformula/Qnum = %s (dev %s); twisted-jet kernel vs untwisted: max rel dev %s"
                % (o1, o2, mp.nstr(ratios[0], 20), mp.nstr(max(abs(v / ratios[0] - 1) for v in ratios), 3), mp.nstr(rat2[0], 20), mp.nstr(max(abs(v / rat2[0] - 1) for v in rat2), 3), mp.nstr(devtw, 3)))
        say("   Frobenius max |P2'/P2 - P3'/P3| at theta_s: %s ; Wronskian max |W(theta_s)|: %s ; max rel dev of W'(theta_s) vs (zeta3-zeta2) e^(sigma theta) Pi' Q1(theta): %s"
            % (mp.nstr(frob, 3), mp.nstr(wr, 3), mp.nstr(dq, 3)))
        if len(oms):
            qn_all = []
        return None

    run_case([Fr(3, 2), Fr(5, 4), Fr(2, 3)], [Fr(2, 7), Fr(-3, 5)], [(Fr(5, 2), Fr(7, 4)), (Fr(-2, 5), Fr(3)), (Fr(9, 7), Fr(-5, 3))], "base")
    run_case([Fr(3, 2), Fr(5, 4), Fr(2, 3)], [Fr(1, 5), Fr(-4, 9)], [(Fr(5, 2), Fr(7, 4))], "theta variant")
    run_case([Fr(7, 5), Fr(-3, 4), Fr(1, 6)], [Fr(2, 7), Fr(-3, 5)], [(Fr(5, 2), Fr(7, 4))], "zeta variant")
    qside.ZETA = zeta_saved


# ======================================================================================================================
if __name__ == "__main__":
    run_L0()
    run_L0b()
    run_L1()
    run_L1_table()
    run_L1c()
    run_L2_ladder()
    run_L2_Z()
    G5, names5, G5t = run_L2_symbolic()
    say("")
    say("kernel G (untwisted jets), entries (w1 = omega1, w1p = omega1', thi = 1/theta12):  basis N =", names5)
    Gs5 = [[sp.factor(G5[i, j]) for j in range(5)] for i in range(5)]
    for i in range(5):
        for j in range(5):
            say("   G[%s | %s] = %s" % (names5[i], names5[j], Gs5[i][j]))
    say("")
    say("kernel G~ (twisted jets), entries (delta = omega1' - omega1, thi = 1/theta12):  basis N~ = (P~3 p~02, P~3' p~02, P~3 p~03, P~3' p~03, P~3 p~23)")
    for i in range(5):
        for j in range(5):
            say("   G~[%s | %s] = %s" % (names5[i], names5[j], G5t[i, j]))
    run_L2_random(Gs5)
    run_L4()
    run_L3(G5, G5t)
    hdr("SUMMARY: the Gaudin two-vacua formula (case (i): N = 10, only n = 0 jets)")
    say("  Qnum_G(B) := <Phi_inf|Psi_B>/<x0^G|Psi_B>,  Phi_inf = Omega_D^inf[zeta,omega1] (x) Omega_D^inf[zeta,omega1'],  x0^G = x^G(zeta) (x) x^G(zeta) (L1b)")
    say("  c(zeta,omega,omega') = (zeta2-omega1)(zeta3-omega1)^2 (zeta2-omega1')(zeta3-omega1')^2 (zeta1-zeta2)^2")
    say("  per site s=1,2:  P3^(k) = (d/du)^k P^3(theta_s), P^j = e^{zeta_j u} p_j,  kappa_s = P3'/P3 = P2'/P2,  p_kl = P2^(k) P3^(l) - P3^(k) P2^(l)  (p_01 = 0, p_12 = kappa p_02)")
    say("  4 c Qnum_G(B) (P3 p02)_1 (P3 p02)_2 = sum_ij G_ij(omega1,omega1',theta12) N_i(theta1) N_j(theta2),  N = (P3 p02, P3' p02, P3 p03, P3' p03, P3 p23)")
    say("  equivalently Qnum_G = (1/4c) sum_ij G_ij nu_i(theta1) nu_j(theta2), nu = (1, kappa, pi, kappa pi, eps),  pi = p03/p02,  eps = p23/p02 = (P3''/P3) pi - P3'''/P3")
    say("  direct Q1 form: p02 = W'(theta_s) = A_s Q1(theta_s), p03 = W'' - kappa W' = A_s[(Lambda_s - kappa) Q1 + 2 Q1'], A_s = (zeta3-zeta2) e^{sigma theta_s} Pi'(theta_s), Lambda_s = 2 sigma + 2/Pi'(theta_s)")
    say("  (sigma = zeta1+zeta2+zeta3, Pi'(theta_1) = theta12, Pi'(theta_2) = -theta12;  Q1 = e^{-zeta1 u} q_1 the direct Q of the first exponent)")
    fit = os.path.join(HERE, "RESULTS-twovac-fit.md")
    say("\nD1 results file present:", os.path.exists(fit), "(the comparison of the kernels is done in RESULTS-twovac-limit.md after reading it)")
    say("\ntotal time %.0fs" % (time.time() - T0))
