"""Q side of the sl(3) Gaudin chain, L = 2, Lambda = (3,2,0) on both sites (agent P4; dual of (3,1,0)); copy of qside_320.py with lambda_2-dependence fixed.
Run top to bottom:   python -X utf8 qside_320.py            (everything; caches in qside_320_*.pkl)
                     python -X utf8 qside_320.py 0 1 2 3    (only the named sections: 0 conventions, 1 Q1, 2 Q2, 3 Q3, 4 Q4, R)
Reuses qside.py / gaudin_chain.py / gl3_reps.py read-only; every function hard-wired to (2,1,0) is re-written here with
(Lambda, gauge power P, |Lambda|) as parameters.
"""
import sys, os, time, itertools, pickle, random
from fractions import Fraction as Fr
from math import factorial
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mpmath as mp
import sympy
import qside as q
from gaudin_chain import *

DPS = 60
mp.mp.dps = DPS
LAM = (3, 2, 0)
ABS = sum(LAM)                    # |Lambda| = residue of S1 at theta_s
ZETA, OM, OMP = list(q.ZETA), list(q.OM), list(q.OMP)
mpf, nst, hdr = q.mpf, q.nst, q.hdr
C2V = 8                           # scalar double pole of S2 (= lam2 (lam1+1))
PG = 3                            # gauge power of the direct oper (confirmed by the scan in Q1)
XXX = ((0, 2), (0, 1), (1, 2))    # row triple (i,j) = (1,3),(1,2),(2,3)

def cache(name, fn):
    path = os.path.join(HERE, name)
    if os.path.exists(path):
        with open(path, 'rb') as f:
            return pickle.load(f)
    val = fn()
    with open(path, 'wb') as f:
        pickle.dump(val, f)
    return val

def build_states(lam=LAM, zeta=None, dps=DPS):
    zeta = zeta or ZETA
    mp.mp.dps = dps
    ch = Chain(lam=lam)
    g = [[zeta[i] if i == j else Fr(0) for j in range(3)] for i in range(3)]
    ops = Operators(ch, g)
    C = ops.fit()
    st = eigenstates(ch, ops, C, dps=dps)
    for s in st:      # mp.eig returns mpc with ~1e-59 imaginary noise in some sectors: the problem is real, clean it
        mxim = max([abs(mp.im(x)) for x in s['vec']] + [abs(mp.im(s['ev']))])
        assert mxim < mp.mpf(10) ** (-dps // 2), mxim
        s['ev'] = mp.re(s['ev'])
        s['vec'] = [mp.re(x) for x in s['vec']]
    res = eigenfunctions(ch, C, st)
    return ch, ops, C, st, res

def th_mp(ch):
    return [mpf(t) for t in ch.theta]

def Pi_poly(th, k=1):
    out = [mp.mpf(1)]
    for _ in range(k):
        out = q.pmul(out, q.pmul([-th[0], 1], [-th[1], 1]))
    return out

def is_scalar_op(op, n):
    return q.is_scalar_op(op, n)

def Q_objects(rec, zeta):
    z = [mpf(x) for x in zeta]
    return ([q.Qfun(z[j], rec['dual'][j][3]) for j in range(3)], [q.Qfun(-z[i], rec['direct'][i][3]) for i in range(3)])

# =============================================================== Section 0
def vac_data(ch, om, zeta, gauge, lam1=ABS, variant=0, solve=True):
    cov = vacuum_covector(ch, zeta, om)
    gp = gprime(zeta, om, variant)
    ops = Operators(ch, gp)
    C = ops.fit()
    ok, res = covector_eigen_test(C, cov, ch)
    th = list(ch.theta)
    out = dict(cov=cov, ok=ok, lam2=res[2], lam3=res[3], gp=gp)
    if solve and ok:
        S1 = q.PF(sum(om), {(0, 1): Fr(lam1), (1, 1): Fr(lam1)}, th)
        S2 = q.pf_S2(res[2], th)
        S3 = q.pf_S3(res[3], th)
        out['Qd'] = [q.solve_poly_exact('A', S1, S2, S3, -om[i], th, gauge=gauge) for i in range(3)]
        out['Pd'] = [q.solve_poly_exact('D', S1, S2, S3, om[i], th) for i in range(3)]
    return out

def nilpotency(E):
    k, P = 1, E
    while not P.is_zero():
        P = P * E
        k += 1
    return k          # smallest k with E^k = 0

def section0(ch):
    hdr('SECTION 0: conventions for Lambda = %s' % (LAM,))
    print('Lambda =', ch.lam, ' |Lambda| =', ABS, ' theta =', [str(t) for t in ch.theta], ' one-site dim =', ch.d1, ' two-site dim =', ch.dim,
          ' weight sectors:', len(ch.sectors), ' largest sector:', max(len(v) for v in ch.sectors.values()))
    idx, ok = ch.check_hw()
    print('HW ket index in gl3_reps ordering (0-based) =', idx, ' GT pattern', ch.rep.patterns[idx], '; annihilated by E12, E13, E23:', ok)
    print('nilpotency indices (smallest k with E^k = 0 on the one-site rep): E12:', nilpotency(ch.E1[(1, 2)]), ' E23:', nilpotency(ch.E1[(2, 3)]), ' E13:', nilpotency(ch.E1[(1, 3)]))
    g = [[ZETA[i] if i == j else Fr(0) for j in range(3)] for i in range(3)]
    ops = Operators(ch, g)
    t0 = time.time()
    C = ops.fit()
    print('partial-fraction fit (S2: 1, 1/x, 1/x^2; S3: 1, 1/x..1/x^3 per site), mismatches at extra points:', C['res2'], C['res3'], '(%.1fs)' % (time.time() - t0))
    print('S1 = e1(zeta) + %d/x1 + %d/x2 (scalar), exact:' % (ABS, ABS), (C['S1'] - SM.eye(ch.dim, S1_scalar(ch, ZETA, Fr(3, 7)))).is_zero())
    a, b = ops.S_at(Fr(3, 7))[1], ops.S_at(Fr(-2, 9))[2]
    print('[S2(3/7), S3(-2/9)] = 0 exactly:', comm(a, b).is_zero())
    pts = [Fr(5, 3), Fr(-4, 5), Fr(7, 2)]
    allc = True
    for ua in pts:
        for ub in pts:
            A2 = ops.eval_S(C, 2, ua); B3 = ops.eval_S(C, 3, ub); B2 = ops.eval_S(C, 2, ub); A3 = ops.eval_S(C, 3, ua)
            allc = allc and comm(A2, B3).is_zero() and comm(A2, B2).is_zero() and comm(A3, B3).is_zero()
    print('all [S_k(u), S_l(v)] = 0 exactly on a 3x3 rational grid (via the fitted Laurent form):', allc)
    wc = True
    for k in (1, 2, 3):
        Ek = ch.Es[(1, k, k)] + ch.Es[(2, k, k)]
        wc = wc and comm(ops.eval_S(C, 2, Fr(3, 7)), Ek).is_zero() and comm(ops.eval_S(C, 3, Fr(-2, 9)), Ek).is_zero()
    print('weight conservation [S2(u), E_kk^tot] = [S3(u), E_kk^tot] = 0:', wc)
    e2 = ZETA[0] * ZETA[1] + ZETA[0] * ZETA[2] + ZETA[1] * ZETA[2]
    e3 = ZETA[0] * ZETA[1] * ZETA[2]
    print('S2 constant scalar =', is_scalar_op(C[2][0], ch.dim)[1], '(e2(zeta) =', e2, ');  S3 constant scalar =', is_scalar_op(C[3][0], ch.dim)[1], '(e3(zeta) =', e3, ')')
    sc1, sc2 = is_scalar_op(C[2][2], ch.dim), is_scalar_op(C[2][4], ch.dim)
    c2 = list(sc1[1])[0]
    print('S2 double-pole coefficients scalar: %s, %s; value c2(Lambda) = %s' % (sc1[0], sc2[0], sc1[1]))
    cas = sum(LAM[i] * (LAM[i] + 3 - 2 * (i + 1)) for i in range(3))
    print('   quadratic Casimir sum_i lam_i (lam_i + 3 - 2i) = %d, half = %s' % (cas, Fr(cas, 2)))
    print('S3 triple-pole coefficient operators identically ZERO:', C[3][3].is_zero(), C[3][6].is_zero())
    print('S3 double-pole coefficient scalar?', is_scalar_op(C[3][2], ch.dim)[0], is_scalar_op(C[3][5], ch.dim)[0], ' (non-scalar: state dependent)')
    print('non-scalar Laurent data: H_s (S2 simple pole):', [not is_scalar_op(C[2][i], ch.dim)[0] for i in (1, 3)], ' b_s (S3 simple):', [not is_scalar_op(C[3][i], ch.dim)[0] for i in (1, 4)],
          ' a_s (S3 double):', [not is_scalar_op(C[3][i], ch.dim)[0] for i in (2, 5)], ' => six non-scalar data (H1,H2,a1,a2,b1,b2)')
    e = sympy.Symbol('e')
    ind = sympy.expand(e * ((e - 1) * (e - 2) - ABS * (e - 1) + sympy.Rational(c2.numerator, c2.denominator)))
    print('indicial polynomial of the dual oper at theta_s: e(e-1)(e-2) - %d e(e-1) + c2 e = %s ; roots %s' % (ABS, sympy.factor(ind), sympy.roots(ind)))
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
    print('H_s = sum_{i<j}[g_ii E_jj + g_jj E_ii - g_ij E_ji - g_ji E_ij]^(s) + (1/(th_s - th_s\')) sum_{i<j}[...] (same formula as (2,1,0)), exact:',
          (C[2][1] - Hs(1)).is_zero(), (C[2][3] - Hs(2)).is_zero())
    return ops, C, c2

def section0_other_reps():
    print('c2(Lambda) for other reps (two-site fit, exact):')
    for lam in [(1, 0, 0), (1, 1, 0), (2, 0, 0), (2, 1, 0), (2, 2, 0), (3, 1, 0), (3, 0, 0), (2, 1, 1), (3, 2, 0)]:
        ch = Chain(lam=lam)
        g = [[ZETA[i] if i == j else Fr(0) for j in range(3)] for i in range(3)]
        C = Operators(ch, g).fit()
        sc = is_scalar_op(C[2][2], ch.dim)
        cas = sum(lam[i] * (lam[i] + 3 - 2 * (i + 1)) for i in range(3))
        print('    %s: c2 = %s   (Casimir/2 = %s, |Lambda| = %d, S3 triple pole zero: %s)' % (lam, sorted(sc[1]) if sc[0] else 'NONSCALAR', Fr(cas, 2), sum(lam), C[3][3].is_zero() and C[3][6].is_zero()))

def section_vacuum(ch):
    hdr("SECTION 0b: bra vacuum and g' orientation")
    r = gauss_r(ZETA, OM)
    print('Gauss parameters r_k(zeta,omega) =', [str(x) for x in r], '(rep-independent formulas)')
    out = {}
    for name, om in (('omega', OM), ('omega_prime', OMP)):
        print('bra twist', name, [str(x) for x in om])
        for var, label in ((0, 'P gc P^-1'), (1, 'P^-1 gc P'), (2, '(P gc P^-1)^T')):
            vd = vac_data(ch, om, ZETA, PG, variant=var, solve=False)
            print("   variant %d  g' = %-14s : Omega x Omega left eigencovector of all S2,S3 coefficient operators: %s" % (var, label, vd['ok']))
        vd = vac_data(ch, om, ZETA, PG)
        out[name] = vd
        l2, l3 = vd['lam2'], vd['lam3']
        print('   exact vacuum eigenvalue functions (x1 = u-theta1, x2 = u-theta2):')
        print('   S1 = %s + %d/x1 + %d/x2' % (sum(om), ABS, ABS))
        print('   S2 = %s + (%s)/x1 + (%s)/x1^2 + (%s)/x2 + (%s)/x2^2' % tuple(str(x) for x in l2))
        print('   S3 = %s + (%s)/x1 + (%s)/x1^2 + (%s)/x1^3 + (%s)/x2 + (%s)/x2^2 + (%s)/x2^3' % tuple(str(x) for x in l3))
        print('   e2(omega) = %s, e3(omega) = %s' % (om[0] * om[1] + om[0] * om[2] + om[1] * om[2], om[0] * om[1] * om[2]))
        print('   vacuum direct Q_i = e^{-om_i u} q_i (gauge %d): degrees %s nullities %s ; vacuum dual P^j: degrees %s nullities %s' % (
            PG, [x[0] for x in vd['Qd']], [x[1] for x in vd['Qd']], [x[0] for x in vd['Pd']], [x[1] for x in vd['Pd']]))
    return out

# =============================================================== Q1
def sector_reps(st):
    seen, out = set(), []
    for k, s in enumerate(st):
        if s['sector'] not in seen:
            seen.add(s['sector']); out.append(k)
    return out

def make_S(s, th, zeta):
    return q.pf_S1(sum(mpf(z) for z in zeta), ABS, th), q.pf_S2(s['lam2'], th), q.pf_S3(s['lam3'], th)

def q1_pscan(ch, st, th, zeta=None):
    zeta = zeta or ZETA
    res = {}
    for k in sector_reps(st):
        s = st[k]
        S1, S2, S3 = make_S(s, th, zeta)
        mu = s['sector']
        for p in range(6):
            out = [q.solve_poly_mp('A', S1, S2, S3, -mpf(zeta[i]), th, dmax=10, gauge=p) for i in range(3)]
            res[(mu, p)] = ([o[0] for o in out], [o[1] for o in out], [o[2] for o in out])
    return res

def q1_run(st, th, zeta=None, gauge=PG):
    zeta = zeta or ZETA
    out = []
    for s in st:
        S1, S2, S3 = make_S(s, th, zeta)
        rec = dict(mu=s['sector'], ev=s['ev'], dual=[], direct=[])
        for j in range(3):
            rec['dual'].append(q.solve_poly_mp('D', S1, S2, S3, mpf(zeta[j]), th, dmax=10))
        for i in range(3):
            rec['direct'].append(q.solve_poly_mp('A', S1, S2, S3, -mpf(zeta[i]), th, dmax=10, gauge=gauge))
        u = mp.mpf('0.37')
        rr = [q.ode_residual('D', S1, S2, S3, mpf(zeta[j]), rec['dual'][j][3], 0, u) for j in range(3)]
        rr += [q.ode_residual('A', S1, S2, S3, -mpf(zeta[i]), rec['direct'][i][3], gauge, u) for i in range(3)]
        rec['resid'] = max(rr)
        out.append(rec)
    return out

def jet_orders(Qfuns, x, K=9, tol=mp.mpf(10) ** -25):
    """vanishing orders at x of the span of the functions (relative Gaussian elimination on Taylor jets)"""
    M = [q.qtaylor(f, x, K) for f in Qfuns]
    for r in M:
        mx = max(abs(v) for v in r) or 1
        for k in range(K):
            r[k] = r[k] / mx
    rows = list(range(len(M)))
    orders = []
    for k in range(K):
        if not rows:
            break
        piv = max(rows, key=lambda r: abs(M[r][k]))
        if abs(M[piv][k]) > tol:
            orders.append(k)
            rows.remove(piv)
            for r in rows:
                f = M[r][k] / M[piv][k]
                for kk in range(k, K):
                    M[r][kk] -= f * M[piv][kk]
    return orders

def local_exponents(recs, th, zeta=None):
    zeta = zeta or ZETA
    tally = {}
    for rec in recs:
        P, Q = Q_objects(rec, zeta)
        for s in range(2):
            key = (tuple(jet_orders(P, th[s])), tuple(jet_orders(Q, th[s])))
            tally[key] = tally.get(key, 0) + 1
    return tally


# ---- Frobenius (no-logarithm) conditions of the dual oper at theta_s, exponents {0,3,5}, |Lambda| = 5, c2 = 8
_FROB = {}
def frobenius_polys(c2v=8):
    """symbolic: with x = u - theta_s, S1 = 4/x + sum s_k x^k, S2 = c2/x^2 + H/x + sum h_k x^k, S3 = a/x^2 + b/x + sum t_k x^k,
    the exponent-0 series needs C1 = 0 (order x^2) and C2 = C3 = 0 (order x^5, coefficients of the two free constants), the exponent-2 series needs C3 = 0 again.
    Returns (vars, [C1, C2, C3])"""
    if 'p' in _FROB:
        return _FROB['p']
    x = sympy.Symbol('x')
    H, a, b = sympy.symbols('H a b')
    N = 6
    sv = sympy.symbols('s0:%d' % N); hv = sympy.symbols('h0:%d' % N); tv = sympy.symbols('t0:%d' % N)
    S1 = ABS / x + sum(sv[k] * x ** k for k in range(N))
    S2 = c2v / x ** 2 + H / x + sum(hv[k] * x ** k for k in range(N))
    S3 = a / x ** 2 + b / x + sum(tv[k] * x ** k for k in range(N))
    kmax = 5
    c = sympy.symbols('k0:%d' % (kmax + 1))
    P = sum(c[k] * x ** k for k in range(kmax + 1))
    expr = sympy.expand((sympy.diff(P, x, 3) - S1 * sympy.diff(P, x, 2) + S2 * sympy.diff(P, x) - S3 * P) * x ** 3)
    eqs = [expr.coeff(x, n) for n in range(kmax + 1)]
    sol, conds = {}, []
    for n in range(kmax + 1):
        eq = sympy.expand(eqs[n].subs(sol))
        if sympy.simplify(eq.coeff(c[n])) == 0:
            conds.append((n, eq))
        else:
            sol[c[n]] = sympy.solve(eq, c[n])[0]
    C = []
    for n, eq in conds:
        if n == 0:
            continue
        poly = sympy.Poly(eq, *c)
        for mon, co in zip(poly.monoms(), poly.coeffs()):
            C.append((n, mon.index(1), sympy.expand(co * 192)))
    vars_ = (H, a, b) + tuple(sv[:3]) + tuple(hv[:3]) + tuple(tv[:3])
    _FROB['p'] = (vars_, C)
    return _FROB['p']

def frobenius_check(st, th):
    """evaluate the no-log conditions on every state and both sites; returns max relative violation per condition"""
    vars_, C = frobenius_polys()
    fns = []
    for (n, idx, poly) in C:
        P = sympy.Poly(poly, *vars_)
        fns.append((n, idx, [(mon, mp.mpf(int(co.p)) / int(co.q)) for mon, co in zip(P.monoms(), P.coeffs())]))
    e1 = sum(mpf(z) for z in ZETA)
    e2 = ZETA[0] * ZETA[1] + ZETA[0] * ZETA[2] + ZETA[1] * ZETA[2]
    e3 = ZETA[0] * ZETA[1] * ZETA[2]
    worst = [mp.mpf(0)] * len(fns)
    for s_ in st:
        l2, l3 = s_['lam2'], s_['lam3']
        data = [(l2[1], l3[2], l3[1]), (l2[3], l3[5], l3[4])]
        for s in range(2):
            o = 1 - s
            d = th[s] - th[o]
            def reg(k, m):
                return (-1) ** k * mp.binomial(m + k - 1, k) * d ** (-m - k)
            H, a, b = data[s]
            Ho, ao, bo = data[o]
            sk = [ABS * reg(k, 1) + (mpf(e1) if k == 0 else 0) for k in range(3)]
            hk = [Ho * reg(k, 1) + C2V * reg(k, 2) + (mpf(e2) if k == 0 else 0) for k in range(3)]
            tk = [bo * reg(k, 1) + ao * reg(k, 2) + (mpf(e3) if k == 0 else 0) for k in range(3)]
            vals = [H, a, b] + sk + hk + tk
            for i, (n, idx, terms) in enumerate(fns):
                tot, scl = mp.mpf(0), mp.mpf(0)
                for mon, co in terms:
                    t_ = co
                    for v, e in zip(vals, mon):
                        if e:
                            t_ = t_ * v ** e
                    tot += t_
                    scl += abs(t_)
                worst[i] = max(worst[i], abs(tot) / scl if scl != 0 else abs(tot))
    return [(n, idx, w) for (n, idx, _), w in zip(fns, worst)]

def section1(ch, st, th, res):
    hdr('SECTION 1 (Q1): opers and Q-functions for Lambda = %s, %d states, %d digits' % (LAM, len(st), DPS))
    t0 = time.time()
    scan = cache('qside_320_pscan.pkl', lambda: q1_pscan(ch, st, th))
    print('  gauge-power scan: direct ODE  Pi^p o A o Pi^-p, Q_i = e^{-zeta_i u} q_i, one state per weight sector (%.0fs)' % (time.time() - t0))
    print('  (smallest polynomial degree d_i, nullity at it; None = none up to degree 9); target degrees 2 lam1 - mu_i = 6 - mu_i')
    mus = sorted({mu for (mu, p) in scan})
    summ = {}
    for p in range(6):
        okdeg = okker = none = 0
        shift = {}
        for mu in mus:
            degs, nul, _ = scan[(mu, p)]
            exp = [6 - m for m in mu]
            okdeg += (degs == exp)
            okker += all(n == 1 for n in nul)
            if any(d is None for d in degs):
                none += 1
            else:
                dd = tuple(d - e for d, e in zip(degs, exp))
                shift[dd] = shift.get(dd, 0) + 1
        summ[p] = (okdeg, okker, none, shift)
        print('    p = %d: sectors with degrees == 6 - mu_i: %2d / %d ; kernel dim 1 for all three i: %2d ; sectors where some i has no solution: %2d ; (deg - (6-mu)) pattern counts: %s' % (
            p, okdeg, len(mus), okker, none, shift))
    good = [p for p in range(6) if summ[p][0] == len(mus) and summ[p][1] == len(mus)]
    print('  => gauge powers p giving ALL degrees 6 - mu_i with kernel dimension 1 in all %d sectors:' % len(mus), good)
    assert good == [PG], good
    t0 = time.time()
    recs = cache('qside_320_recs.pkl', lambda: q1_run(st, th))
    assert all(r['mu'] == s['sector'] and abs(r['ev'] - s['ev']) < mp.mpf(10) ** -40 for r, s in zip(recs, st)), 'stale cache'
    print('  full run on all %d states (p = %d), %.0fs; max eigen-residual of the eigenstates %s' % (len(st), PG, time.time() - t0, nst(res)))
    bys = {}
    for r in recs:
        bys.setdefault(r['mu'], []).append(r)
    print('  %-10s %-4s %-18s %-20s %-6s %-9s' % ('sector mu', '#st', 'deg dual p_j', 'deg direct q_i', 'check', 'max resid'))
    allok = True
    for mu in sorted(bys):
        rs = bys[mu]
        dd = {tuple(x[0] for x in r['dual']) for r in rs}
        di = {tuple(x[0] for x in r['direct']) for r in rs}
        nl = {tuple(x[1] for x in r['dual'] + r['direct']) for r in rs}
        ok = dd == {tuple(mu)} and di == {tuple(6 - m for m in mu)} and nl == {(1,) * 6}
        allok = allok and ok
        print('  %-10s %-4d %-18s %-20s %-6s %-9s' % (mu, len(rs), sorted(dd)[0] if len(dd) == 1 else dd, sorted(di)[0] if len(di) == 1 else di, 'OK' if ok else 'FAIL', nst(max(r['resid'] for r in rs))))
    print('  ALL %d states: dual degrees mu_j, direct degrees 6 - mu_i, kernel dimension 1 for each of 6 exponentials:' % len(st), allok)
    print('  max ODE residual (u=0.37) over all states x 6 Q:', nst(max(r['resid'] for r in recs)), '; largest smallest-singular-value at the degree found:', nst(max(x[2] for r in recs for x in r['dual'] + r['direct'])))
    tal = local_exponents(recs, th)
    print('  local exponents at theta_s (vanishing orders of the span of the 3 dual P^j | of the 3 direct Q_i; psi = Pi^-%d Q has orders minus %d), both sites, all states: counts' % (PG, PG))
    for k, v in sorted(tal.items()):
        print('    dual', k[0], ' direct Q', k[1], ' count', v)
    vars_, Cs = frobenius_polys()
    print('  no-logarithm (Frobenius) conditions of the dual oper at theta_s (exponents 0,3,5; Taylor data of S1,S2,S3 at theta_s: s_k, h_k, t_k from the other site, e_k(zeta)):')
    print('    exponent-0 series: c1 = a c0/8, c2 = (a^2 - H a + 8 b) c0/48 (so dP/P = a/8, ddP/P = (a^2 - H a + 8 b)/24); resonance at order x^3 (exponent 3) and at order x^5 (exponent 5) give:')
    for n, idx, poly in Cs:
        print('    order x^%d, coefficient of free constant c_%d (x192): %s' % (n, idx, poly))
    print('    (see the polynomials above: one condition at order x^3 in c_0, two at order x^5 in c_0 and c_3) => 3 conditions per site, 6 in total for the 6 unknowns (H1,H2,a1,a2,b1,b2)')
    fc = frobenius_check(st, th)
    for n, idx, w in fc:
        print('    max relative violation over all %d states x 2 sites of the condition (order %d, c_%d): %s' % (len(st), n, idx, nst(w)))
    mx1, mx2, rk3, rk5, rk6, s5, s6 = kappa_checks(recs, st, th)
    print('  kappa_s = dP(th_s)/P(th_s) = a_s/c2 = a_s/%d and kappa2_s = ddP(th_s)/P(th_s) = (a_s^2 - H_s a_s + %d b_s)/%d, all 3 P^j, all %d states, both sites: max rel dev %s, %s' % (C2V, C2V, 3 * C2V, len(st), nst(mx1), nst(mx2)))
    print('  jet matrix ranks [P^j(k)(th_s)/k!], j = 1..3 (900 / 450 cases... counted per (state,site)): k<=2 (3x3): %s ; k<=4 (3x5): %s ; k<=5 (3x6): %s ; min sigma_2/sigma_1 (3x5) %s, min sigma_3/sigma_1 (3x6) %s' % (rk3, rk5, rk6, nst(s5), nst(s6)))
    return recs

def section1_neg(ch, st, th):
    res = {}
    for k in sector_reps(st):
        s = st[k]
        S1, S2, S3 = make_S(s, th, ZETA)
        r = {}
        for label, which, sg, gg in (('D,-zeta', 'D', -1, 0), ('A g0,-zeta', 'A', -1, 0), ('A gP,+zeta', 'A', 1, PG)):
            r[label] = [q.solve_poly_mp(which, S1, S2, S3, sg * mpf(ZETA[i]), th, dmax=9, gauge=gg)[0] for i in range(3)]
        res[s['sector']] = r
    print('  negative controls (one state per sector, %d states; smallest degree for i = 1,2,3; None = no solution up to degree 9):' % len(res))
    for l, nm in (('D,-zeta', 'dual ODE D with e^{-zeta u}'), ('A g0,-zeta', 'plain formal adjoint (p = 0) with e^{-zeta u}'), ('A gP,+zeta', 'adjoint with p = %d and e^{+zeta u}' % PG)):
        vals = {tuple(res[mu][l]) for mu in res}
        print('    %-50s : %s' % (nm, sorted(vals, key=str)[:4]))


def kappa_checks(recs, st, th):
    """local exponents {0,3,5}: P(theta+x) = c0 (1 + kappa x + (kappa2/2) x^2) + O(x^3) for every dual solution.
    Indicial recursion at orders x^1, x^2 gives  kappa_s = a_s/c2,  kappa2_s = (a_s^2 - H_s a_s + c2 b_s)/(3 c2)   (c2 = 8: a/8, (a^2 - H a + 8 b)/24)."""
    mx1 = mx2 = mp.mpf(0)
    rk5 = {}; rk6 = {}; rk3 = {}
    sg5 = []; sg6 = []
    for rec, s in zip(recs, st):
        P, Q = Q_objects(rec, ZETA)
        l2, l3 = s['lam2'], s['lam3']
        data = [(l2[1], l3[2], l3[1]), (l2[3], l3[5], l3[4])]
        for site in range(2):
            H, a, b = data[site]
            k1 = a / C2V
            k2 = (a * a - H * a + C2V * b) / (3 * C2V)
            M = [[P[j](th[site], k) for k in range(6)] for j in range(3)]
            for j in range(3):
                mx1 = max(mx1, abs(M[j][1] / M[j][0] - k1) / (abs(k1) + 1))
                mx2 = max(mx2, abs(M[j][2] / M[j][0] - k2) / (abs(k2) + 1))
            def rank_of(cols):
                Mm = mp.matrix([[M[j][k] / mp.factorial(k) for k in cols] for j in range(3)])
                for j in range(3):
                    nr = max(abs(Mm[j, c]) for c in range(Mm.cols))
                    for c in range(Mm.cols):
                        Mm[j, c] /= nr
                Sg = mp.svd_r(Mm, compute_uv=False)
                return sum(1 for x in Sg if x / Sg[0] > mp.mpf(10) ** -30), [x / Sg[0] for x in Sg]
            r5, S5 = rank_of(range(5)); r6, S6 = rank_of(range(6)); r3, S3_ = rank_of(range(3))
            rk5[r5] = rk5.get(r5, 0) + 1; rk6[r6] = rk6.get(r6, 0) + 1; rk3[r3] = rk3.get(r3, 0) + 1
            sg5.append(S5[1]); sg6.append(S6[2])
    return mx1, mx2, rk3, rk5, rk6, min(sg5), min(sg6)

def vac_closed_forms(ch):
    hdr('SECTION 0c: vacuum Laurent data vs closed forms')
    th = [mpf(t) for t in ch.theta]
    th12 = ch.theta[0] - ch.theta[1]
    for name, om in (('omega', OM), ('omega_prime', OMP)):
        vd = vac_data(ch, om, ZETA, PG)
        l2, l3 = vd['lam2'], vd['lam3']
        Hs = [l2[1], l2[3]]; a = [l3[2], l3[5]]; b = [l3[1], l3[4]]
        base = sum((ABS - LAM[i]) * om[i] for i in range(3))
        for k in (12,):
            print('  %s: H_s - sum_i(|L|-lam_i) om_i = %s, %s ; times th12 = %s, %s  (predicted +-k = +-%d)' % (name, str(Hs[0] - base), str(Hs[1] - base), str((Hs[0] - base) * th12), str((Hs[1] - base) * th12), k))
        print('  %s: a_s = %s, %s ; a_s/om3 = %s, %s (predicted |Lambda| = %d, observed value = c2 = %d);  b_s - om3 (H_s - |L| om3) = %s, %s' % (
            name, str(a[0]), str(a[1]), str(a[0] / om[2]), str(a[1] / om[2]), ABS, C2V, str(b[0] - om[2] * (Hs[0] - ABS * om[2])), str(b[1] - om[2] * (Hs[1] - ABS * om[2]))))
        print('  %s: b_s - om3 (H_s - c2 om3) = %s, %s (for comparison)' % (name, str(b[0] - om[2] * (Hs[0] - C2V * om[2])), str(b[1] - om[2] * (Hs[1] - C2V * om[2]))))

# =============================================================== Q2
def pi_pow_dev(F, th):
    k = (len(F) - 1) // 2
    if len(F) - 1 != 2 * k:
        return None, mp.mpf(1)
    Pk = Pi_poly(th, k)
    return k, max(abs(F[i] - Pk[i]) for i in range(len(F)))

def q2_wronskians(recs, th, zeta=None):
    zeta = zeta or ZETA
    z = [mpf(x) for x in zeta]
    rows = []
    for rec in recs:
        dq = [x[3] for x in rec['direct']]
        dp = [x[3] for x in rec['dual']]
        for (a, b, c) in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):
            w = q.wronskian_direct(dq[a], z[a], dq[b], z[b])
            quo, rem = q.pdivmod(q.trim(w), dp[c])
            remn = max(abs(x) for x in rem) / max(abs(x) for x in w)
            lead = quo[-1]
            F = [x / lead for x in quo]
            k1, dev1 = pi_pow_dev(F, th)
            w2 = q.wronskian_dual(dp[a], z[a], dp[b], z[b])
            quo2, rem2 = q.pdivmod(q.trim(w2), dq[c])
            remn2 = max(abs(x) for x in rem2) / max(abs(x) for x in w2)
            lead2 = quo2[-1]
            F2 = [x / lead2 for x in quo2]
            k2, dev2 = pi_pow_dev(F2, th)
            rows.append(dict(mu=rec['mu'], pair=(a, b), c=c, deg_w=len(q.trim(w)) - 1, deg_pc=len(dp[c]) - 1, rem=remn, kF=k1, F_dev=dev1,
                             deg_w2=len(q.trim(w2)) - 1, deg_qc=len(dq[c]) - 1, rem2=remn2, kF2=k2, F2_dev=dev2,
                             cab=lead / (z[a] - z[b]), cab2=lead2 / (z[b] - z[a])))
    return rows

def q2_weyl(recs, st, th, gauge=PG, zeta=None, upts=('0.37', '1.9', '-1.3')):
    zeta = zeta or ZETA
    z = [mpf(x) for x in zeta]
    S1f = q.pf_S1(sum(z), ABS, th)
    det3 = q.det3
    out = []
    for rec, s in zip(recs, st):
        S2 = q.pf_S2(s['lam2'], th); S3 = q.pf_S3(s['lam3'], th)
        P, Qd = Q_objects(rec, zeta)
        dev = dict(S1=0, S2=0, S3=0, S1d=0, S2d=0, S3d=0, Wc=[])
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
            dev['Wc'].append(W / (mp.exp(sum(z) * u) * ((u - th[0]) * (u - th[1])) ** ABS))
            Dd = []
            for i in range(3):
                ps = [Qd[i](u, k) / factorial(k) for k in range(4)]
                for t in th:
                    ps = q.ser_mul(ps, q.ser_pow(u - t, -gauge))
                Dd.append([ps[k] * factorial(k) for k in range(4)])
            Wd = det3([[Dd[i][0], Dd[i][1], Dd[i][2]] for i in range(3)])
            t3 = det3([[Dd[i][3], Dd[i][1], Dd[i][2]] for i in range(3)]) / Wd
            t2 = -det3([[Dd[i][0], Dd[i][3], Dd[i][2]] for i in range(3)]) / Wd
            t1 = det3([[Dd[i][0], Dd[i][1], Dd[i][3]] for i in range(3)]) / Wd
            e1 = -S1f(u); e2 = 2 * S1f(u, 1) + S2(u); e3 = -(S1f(u, 2) + S2(u, 1) + S3(u))
            dev['S1d'] = max(dev['S1d'], abs(t1 - e1) / abs(e1))
            dev['S2d'] = max(dev['S2d'], abs(t2 - e2) / abs(e2))
            dev['S3d'] = max(dev['S3d'], abs(t3 - e3) / abs(e3))
        out.append(dev)
    return out

def section2(ch, st, th, recs):
    hdr('SECTION 2 (Q2): Wronskian / QQ relations and quantum-Weyl formulas, Lambda = %s' % (LAM,))
    rows = q2_wronskians(recs, th)
    z = [mpf(x) for x in ZETA]
    n = len(recs)
    print('  (a) direct: Q_a Q_b\' - Q_b Q_a\' = c_ab e^{-(z_a+z_b)u} w_ab(u),  w_ab = F p_c  (c = complement; all p, q monic)')
    print('      %d states x 3 pairs: max relative remainder of w_ab / p_c = %s' % (n, nst(max(r['rem'] for r in rows))))
    print('      deg F = deg w - deg p_c in', sorted({r['deg_w'] - r['deg_pc'] for r in rows}), '; F = Pi^k with k in', sorted({r['kF'] for r in rows}, key=str), ', max |F - Pi^k| =', nst(max(r['F_dev'] for r in rows)))
    print('      max |c_ab/(zeta_a - zeta_b) - 1| =', nst(max(abs(r['cab'] - 1) for r in rows)))
    print("  (a') dual: P^a P^b' - P^b P^a' = c'_ab e^{+(z_a+z_b)u} F'(u) q_c")
    print('      deg F\' in', sorted({r['deg_w2'] - r['deg_qc'] for r in rows}), '; F\' = Pi^k with k in', sorted({r['kF2'] for r in rows}, key=str), '; max relative remainder', nst(max(r['rem2'] for r in rows)),
          '; max |F\' - Pi^k| =', nst(max(r['F2_dev'] for r in rows)), '; max |c\'_ab/(zeta_b - zeta_a) - 1| =', nst(max(abs(r['cab2'] - 1) for r in rows)))
    wy = q2_weyl(recs, st, th)
    print('  (b) quantum Weyl (monic operator killing P^1,P^2,P^3): S1 = det[P,P\',P\'\'\']/W, S2 = -det[P,P\'\'\',P\'\']/W, S3 = det[P\'\'\',P\',P\'\']/W, W = det[P,P\',P\'\']')
    for k in ('S1', 'S2', 'S3'):
        print('      max relative deviation of %s from the eigenvalue function (%d states x 3 u): %s' % (k, n, nst(max(w[k] for w in wy))))
    vand = (z[1] - z[0]) * (z[2] - z[0]) * (z[2] - z[1])
    wc = [c for w in wy for c in w['Wc']]
    print('      W = c_W e^{(z1+z2+z3)u} Pi^%d (power = |Lambda|); max |c_W/prod_{i<j}(zeta_j - zeta_i) - 1| = %s (prod = %s)' % (ABS, nst(max(abs(c / vand - 1) for c in wc)), nst(vand, 10)))
    print('      direct side: psi_i = Pi^-%d Q_i span ker of the monic operator d^3 + (-S1) d^2 + (2S1\' + S2) d - (S1\'\' + S2\' + S3):' % PG)
    for k in ('S1d', 'S2d', 'S3d'):
        print('      max relative deviation %s: %s' % (k, nst(max(w[k] for w in wy))))

# =============================================================== Q3
COLS = q.COLS

def res_functional(tdP, jdx, tdQ, idx, s, a, m, sp, th, gauge, n=9):
    other = 1 - s
    Ps = q.site_series(tdP, s, jdx, a, n)
    Qs = q.site_series(tdQ, s, idx, 0, n)
    G = q.ser_mul(Ps, Qs, n)
    G = q.ser_mul(G, q.ser_pow(th[s] - th[other], -gauge, n), n)
    ms = 0
    if sp == s:
        ms = m
    else:
        G = q.ser_mul(G, q.ser_pow(th[s] - th[sp], -m, n), n)
    return G[gauge + ms - 1]

def Nrows(tdP, tdQ, th, gauge):
    out = {}
    for s in range(2):
        for i in range(3):
            for j in range(3):
                row = []
                for col in COLS:
                    kind, sp = col[0], int(col[1]) - 1
                    if kind == 'H':
                        v = -res_functional(tdP, j, tdQ, i, s, 1, 1, sp, th, gauge)
                    elif kind == 'a':
                        v = res_functional(tdP, j, tdQ, i, s, 0, 2, sp, th, gauge)
                    else:
                        v = res_functional(tdP, j, tdQ, i, s, 0, 1, sp, th, gauge)
                    row.append(v)
                out[(s, i, j)] = row
    return out

def ranks_from_sv(S, thr=mp.mpf(10) ** -40):
    return sum(1 for x in S if x > thr)

def rel_det(N, trip=XXX):
    M = q.build_N_det(N, trip)
    nr = mp.mpf(1)
    for r in range(M.rows):
        nr *= mp.sqrt(sum(abs(M[r, c]) ** 2 for c in range(6)))
    return abs(mp.det(M)) / nr

def section3(ch, st, th, recs):
    hdr('SECTION 3 (Q3): residue functionals, orthogonality, det N, x0-free overlap ratio, Lambda = %s, gauge p = %d' % (LAM, PG))
    nS = len(st)
    objs = [Q_objects(r, ZETA) for r in recs]
    tdP = [q.SiteData(o[0], th, 8) for o in objs]
    tdQ = [q.SiteData(o[1], th, 8) for o in objs]
    rng = random.Random(7)
    pairs = [(a, b) for a in range(nS) for b in range(nS) if a != b]
    sample = rng.sample(pairs, 600)
    print('  Lagrange identity (bra A: psi = Pi^-%d Q_i, ket B: P^j; Delta S_k = S_k^B - S_k^A):  P (dS2 psi)\' + dS3 P psi = C\'; the double pole of S2 is the scalar c2 (cancels), S3 has no triple pole,' % PG)
    print('  so Res_{theta_s} gives the SAME six columns: N[(s,i,j),H_t] = -Res[psi P\'/(u-th_t)], N[.,a_t] = Res[psi P/(u-th_t)^2], N[.,b_t] = Res[psi P/(u-th_t)], and sum_m Delta_m N_m = 0.')
    mx = mp.mpf(0)
    t0 = time.time()
    ranks_ne = {}
    nullcheck = []
    for (A, B) in sample:
        N = Nrows(tdP[B], tdQ[A], th, PG)
        D = q.delta_vec(st[B], st[A])
        for key, row in N.items():
            v = abs(sum(r_ * d for r_, d in zip(row, D)))
            sc = sum(abs(r_ * d) for r_, d in zip(row, D))
            if sc > 0:
                mx = max(mx, v / sc)
        S = q.sv_ratio(N)
        rk = ranks_from_sv(S)
        same = st[A]['sector'] == st[B]['sector']
        ranks_ne[(same, rk)] = ranks_ne.get((same, rk), 0) + 1
        nullcheck.append((A, B, rk, S))
    print('  (a) identity residual max |sum_m Delta_m N| / sum_m |Delta_m N| over %d random ordered pairs A != B x 18 rows (%.0fs): %s' % (len(sample), time.time() - t0, nst(mx)))
    print('      rank tally (stacked 18x6, columns normalised, rank = #{sigma_k/sigma_1 > 1e-40}) for A != B, key (same sector?, rank):', dict(sorted(ranks_ne.items())))
    nn = {}
    for (A, B, rk, S) in nullcheck:
        dm = tuple(b - a for a, b in zip(st[A]['sector'], st[B]['sector']))
        key = (sum(1 for x in dm if x > 0), sum(1 for x in dm if x < 0), rk)
        nn[key] = nn.get(key, 0) + 1
    print('      same pairs tallied by (#positive, #negative components of mu_B - mu_A, rank):', dict(sorted(nn.items())), '  [rank 4 <=> exactly two negative components, as in (2,1,0)]')
    r5 = [S[4] for (_, _, rk, S) in nullcheck if rk == 5]
    print('      max sigma_6/sigma_1 over A != B pairs =', nst(max(S[5] for (_, _, rk, S) in nullcheck)), '; min sigma_5/sigma_1 among rank-5 pairs =', nst(min(r5)) if r5 else None,
          '; max sigma_5/sigma_1 among rank-4 pairs =', nst(max([S[4] for (_, _, rk, S) in nullcheck if rk == 4] or [mp.mpf(0)])))
    diag = [q.sv_ratio(Nrows(tdP[B], tdQ[B], th, PG)) for B in range(nS)]
    dk = {}
    for S in diag:
        dk[ranks_from_sv(S)] = dk.get(ranks_from_sv(S), 0) + 1
    print('      A = B (%d states): rank tally' % nS, dk, '; min sigma_6/sigma_1 =', nst(min(S[5] for S in diag)))
    ne = max(rel_det(Nrows(tdP[B], tdQ[A], th, PG)) for (A, B) in sample[:150])
    dd_ = [rel_det(Nrows(tdP[B], tdQ[B], th, PG)) for B in range(nS)]
    print('  (b) det N (6x6, rows (s,(i,j)) with (i,j) in {(1,3),(1,2),(2,3)}, columns H1,H2,a1,a2,b1,b2): |det N|/prod(row norms): max over 150 pairs A != B = %s ; A = B over all %d states: min %s, max %s' % (nst(ne), nS, nst(min(dd_)), nst(max(dd_))))
    hdr('SECTION 3c: overlap-ratio identity R_B(w)/R_B(w\') = [D_B(w)/D_B(w\')][D_B0(w\')/D_B0(w)], B0 = e1 x e1, all %d states' % nS)
    vac = {}
    for name, om in (('w', OM), ('wp', OMP)):
        vd = vac_data(ch, om, ZETA, PG)
        Qs = [q.Qfun(-mpf(om[i]), [mpf(c) for c in vd['Qd'][i][2]]) for i in range(3)]
        Ps = [q.Qfun(mpf(om[i]), [mpf(c) for c in vd['Pd'][i][2]]) for i in range(3)]
        vac[name] = dict(cov=vd['cov'], td=q.SiteData(Qs, th, 8), tdP=q.SiteData(Ps, th, 8), ok=vd['ok'], degQ=[x[0] for x in vd['Qd']], degP=[x[0] for x in vd['Pd']])
        print('  vacuum', name, ': eigencovector ok', vd['ok'], ' direct Q degrees', vac[name]['degQ'], ' dual P degrees', vac[name]['degP'])
    def overlap(cov, s):
        return sum(mpf(c) * x for c, x in zip(cov, s['vec']) if c != 0)
    Rw = [overlap(vac['w']['cov'], s) for s in st]
    Rp = [overlap(vac['wp']['cov'], s) for s in st]
    top = [k for k, s in enumerate(st) if s['sector'] == tuple(2 * x for x in LAM)][0]
    print('  top state B0 (weight %s) has index %d; overlaps nonzero for all states: %s' % (tuple(2 * x for x in LAM), top, all(abs(x) > mp.mpf(10) ** -30 for x in Rw + Rp)))
    NC = {}
    def Nvac(B, name):
        if (B, name) not in NC:
            NC[(B, name)] = Nrows(tdP[B], vac[name]['td'], th, PG)
        return NC[(B, name)]
    def ratio_test(trip, Nfun):
        dw = [mp.det(q.build_N_det(Nfun(B, 'w'), trip)) for B in range(nS)]
        dp = [mp.det(q.build_N_det(Nfun(B, 'wp'), trip)) for B in range(nS)]
        devs, nd = [], 0
        for B in range(nS):
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
    mdev, nd, devs, dw, dp = ratio_test(XXX, Nvac)
    print('  main result (rows {(1,3),(1,2),(2,3)} x 2 sites, columns H1,H2,a1,a2,b1,b2, %d digits): max relative deviation over %d states = %s (non-degenerate: %d)' % (DPS, nS, nst(mdev), nd))
    print('     5 worst states (index, sector, dev):', [(k, st[k]['sector'], nst(d)) for k, d in sorted([(k, d) for k, d in enumerate(devs) if d is not None], key=lambda t: -t[1])[:5]])
    for B in (0, 1, 2, 3):
        print('     B=%d %s  LHS = %s  RHS = %s' % (B, st[B]['sector'], nst(Rw[B] / Rp[B], 15), nst((dw[B] / dp[B]) * (dp[top] / dw[top]), 15)))
    allpairs = [(i, j) for i in range(3) for j in range(3)]
    ok_t, bad_t, deg_t = [], [], []
    t0 = time.time()
    for trip in itertools.combinations(allpairs, 3):
        md, nd2, _, _, _ = ratio_test(trip, Nvac)
        if nd2 < nS:
            deg_t.append((trip, nd2))
        elif md < mp.mpf(10) ** -30:
            ok_t.append(trip)
        else:
            bad_t.append((trip, md))
    print('  84 row triples (same triple on both sites; %.0fs): pass (all %d states non-degenerate, dev < 1e-30): %d ; fail: %d ; degenerate: %d' % (time.time() - t0, nS, len(ok_t), len(bad_t), len(deg_t)))
    def is_perm(trip):
        return len({i for i, j in trip}) == 3 and len({j for i, j in trip}) == 3
    print('     failing triples:', [(tuple((i + 1, j + 1) for i, j in t), nst(m, 2)) for t, m in bad_t][:8], '; all failing are permutation triples:', all(is_perm(t) for t, _ in bad_t), '; number of permutation triples:', sum(1 for t in itertools.combinations(allpairs, 3) if is_perm(t)))
    print('     degenerate triples (first 8, with number of non-degenerate states):', [(tuple((i + 1, j + 1) for i, j in t), n_) for t, n_ in deg_t][:8],
          '; all degenerate are "same i or same j" triples:', all(len({i for i, j in t}) == 1 or len({j for i, j in t}) == 1 for t, _ in deg_t))
    NCm = {}
    def Nmod(B, name):
        if (B, name) not in NCm:
            N = Nvac(B, name)
            N2 = {}
            for key, row in N.items():
                s, i, j = key
                row2 = list(row)
                row2[3] = res_functional(tdP[B], j, vac[name]['td'], i, s, 0, 3, 1, th, PG)
                N2[key] = row2
            NCm[(B, name)] = N2
        return NCm[(B, name)]
    md2, nd2, _, _, _ = ratio_test(XXX, Nmod)
    print('  control 1 (a2 column replaced by Res[psi P/(u-th2)^3]): max deviation', nst(md2) if md2 is not None else None, '(non-degenerate', nd2, ')')
    NCs = {}
    def Nswap(B, name):
        if (B, name) not in NCs:
            NCs[(B, name)] = Nrows(vac[name]['tdP'], tdQ[B], th, PG)
        return NCs[(B, name)]
    md3, nd3, _, _, _ = ratio_test(XXX, Nswap)
    print('  control 2 (orientation swapped: vacuum dual in ket slot, state direct in bra slot): max deviation', nst(md3) if md3 is not None else None, '(non-degenerate', nd3, ')')
    for g_ in (2, 4):
        NCg = {}
        def Ng(B, name, g_=g_, NCg=NCg):
            if (B, name) not in NCg:
                NCg[(B, name)] = Nrows(tdP[B], vac[name]['td'], th, g_)
            return NCg[(B, name)]
        mdg, ndg, _, _, _ = ratio_test(XXX, Ng)
        print('  control 3 (weight Pi^-%d instead of Pi^-%d in the functionals, same Q): max deviation' % (g_, PG), nst(mdg) if mdg is not None else None, '(non-degenerate', ndg, ')')
    return dict(dev=mdev, ok=ok_t, bad=bad_t, deg=deg_t)

# =============================================================== Q4
PHI = [Fr(1, 3), Fr(-2, 7), Fr(5, 11)]
OMH = [Fr(3, 4), Fr(-1, 5), Fr(1, 2)]
OMPH = [Fr(-12, 5), Fr(1, 8), Fr(-1, 3)]
HW_SECTORS = {(6, 4, 0): 1, (5, 5, 0): 1, (6, 3, 1): 1, (5, 4, 1): 2, (6, 2, 2): 1, (5, 3, 2): 2, (4, 4, 2): 1, (4, 3, 3): 1}

def raising_residual(ch, s):
    v = s['vec']
    nv = mp.sqrt(sum(abs(x) ** 2 for x in v))
    tot = 0
    for (i, j) in ((1, 2), (2, 3)):
        E = ch.Es[(1, i, j)] + ch.Es[(2, i, j)]
        w = E.apply(v)
        tot += sum(abs(x) ** 2 for x in w)
    return mp.sqrt(tot) / nv

def q4_point(eps, dps=100):
    q.TOL = mp.mpf(10) ** -50          # kernel threshold for the 100-digit regime solves (as in round 1's set_regime)
    zeta = [eps * p for p in PHI]
    om = [eps * x for x in OMH]
    omp = [eps * x for x in OMPH]
    ch, ops, C, st, res = build_states(LAM, zeta, dps)
    th = th_mp(ch)
    bysec = {}
    for k, s in enumerate(st):
        s['rr'] = raising_residual(ch, s)
        bysec.setdefault(s['sector'], []).append(k)
    br, gaps = [], {}
    for mu, m in HW_SECTORS.items():
        ks = sorted(bysec[mu], key=lambda k: st[k]['rr'])
        br += ks[:m]
        gaps[mu] = [nst(st[k]['rr']) for k in ks[:m + 1]]
    S1 = q.pf_S1(sum(mpf(z) for z in zeta), ABS, th)
    recs = {}
    for k in br:
        s = st[k]
        S2 = q.pf_S2(s['lam2'], th); S3 = q.pf_S3(s['lam3'], th)
        rec = dict(mu=s['sector'], dual=[], direct=[])
        for j in range(3):
            rec['dual'].append(q.solve_poly_mp('D', S1, S2, S3, mpf(zeta[j]), th, dmax=10))
        for i in range(3):
            rec['direct'].append(q.solve_poly_mp('A', S1, S2, S3, -mpf(zeta[i]), th, dmax=10, gauge=PG))
        recs[k] = rec
    out = dict(eps=eps, res=res, gaps=gaps, branches=br, sectors=[st[k]['sector'] for k in br])
    out['degrees'] = {k: ([x[0] for x in recs[k]['dual']], [x[0] for x in recs[k]['direct']], [x[1] for x in recs[k]['dual'] + recs[k]['direct']]) for k in br}
    vac = {}
    for name, o in (('w', om), ('wp', omp)):
        vd = vac_data(ch, o, zeta, PG)
        Qs = [q.Qfun(-mpf(o[i]), [mpf(c) for c in vd['Qd'][i][2]]) for i in range(3)]
        vac[name] = dict(cov=vd['cov'], td=q.SiteData(Qs, th, 8), ok=vd['ok'], deg=[x[0] for x in vd['Qd']])
    out['vac_ok'] = (vac['w']['ok'], vac['wp']['ok'], vac['w']['deg'], vac['wp']['deg'])
    out['r_k'] = gauss_r(zeta, om)
    objs = {k: Q_objects(recs[k], zeta) for k in br}
    tdP = {k: q.SiteData(objs[k][0], th, 8) for k in br}
    def overlap(cov, s):
        return sum(mpf(c) * x for c, x in zip(cov, s['vec']) if c != 0)
    data = {}
    for k in br:
        d = {}
        for name in ('w', 'wp'):
            N = Nrows(tdP[k], vac[name]['td'], th, PG)
            M = q.build_N_det(N, XXX)
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
    top = [k for k in br if st[k]['sector'] == tuple(2 * x for x in LAM)][0]
    for k in br:
        d = data[k]
        L = d['R_w'] * d['det_wp'] * data[top]['det_w']
        R = d['R_wp'] * d['det_w'] * data[top]['det_wp']
        d['dev'] = abs(L - R) / max(abs(L), abs(R))
    out['data'] = data
    out['top'] = top
    return out

def section4():
    hdr('SECTION 4 (Q4): untwisting regime zeta = eps*phi, omega = eps*om^, omega\' = eps*om^\', Lambda = %s, 100 digits' % (LAM,))
    print('  phi =', [str(x) for x in PHI], ' om^ =', [str(x) for x in OMH], " om^' =", [str(x) for x in OMPH])
    pts = {}
    for eps in (Fr(1, 10 ** 6), Fr(2, 10 ** 6)):
        t0 = time.time()
        pts[eps] = cache('qside_320_q4_%s.pkl' % str(eps).replace('/', '_'), lambda eps=eps: q4_point(eps))
        o = pts[eps]
        print('  eps = %s (%.0fs): eigen-residual %s; Gauss r_k = %s' % (str(eps), time.time() - t0, nst(o['res']), [nst(mpf(x), 8) for x in o['r_k']]))
        print('     vacuum eigencovector test (omega, omega\'):', o['vac_ok'][:2], ' vacuum direct Q degrees', o['vac_ok'][2], o['vac_ok'][3])
        print('     raising residuals (smallest first, m+1 per HW sector):', o['gaps'])
        for k in o['branches']:
            dd = o['degrees'][k]
            d = o['data'][k]
            print('     branch idx %3d sector %s : deg dual %s, deg direct %s, nullities %s | |det N(w)| = %s, det/prod(col norms) = %s, sigma_6/sigma_1 = %s | ratio-test deviation = %s' % (
                k, o['sectors'][o['branches'].index(k)], dd[0], dd[1], dd[2], nst(abs(d['det_w'])), nst(d['relDetCol_w']), nst(d['svratio_w']), nst(d['dev'])))
    e1, e2 = Fr(1, 10 ** 6), Fr(2, 10 ** 6)
    a, b = pts[e1], pts[e2]
    print('  epsilon-exponents p in f ~ eps^p, p = log2(f(2e-6)/f(1e-6)):')
    for k1, k2 in zip(a['branches'], b['branches']):
        sa = a['sectors'][a['branches'].index(k1)]
        if sa != b['sectors'][b['branches'].index(k2)]:
            continue
        da, db = a['data'][k1], b['data'][k2]
        ex = lambda f: nst(mp.log(abs(db[f]) / abs(da[f])) / mp.log(2), 4)
        print('     sector %s: |det N(w)| p = %s, |det N(wp)| p = %s, det/prod(col norms) p = %s, sigma_6/sigma_1 p = %s, max|N(w)| p = %s, |R(w)| p = %s' % (
            sa, ex('det_w'), ex('det_wp'), ex('relDetCol_w'), ex('svratio_w'), ex('maxentry_w'), ex('R_w')))
    return pts

# =============================================================== R: the rank-4 puzzle of (2,1,0)
def R_setup():
    lam0 = (2, 1, 0)
    q.GAUGE = 2
    ch, ops, C, st, res = build_states(lam0, ZETA, DPS)
    th = th_mp(ch)
    with open(os.path.join(HERE, 'qside_cache.pkl'), 'rb') as f:
        recs = pickle.load(f)
    assert all(r['mu'] == s['sector'] and abs(r['ev'] - s['ev']) < mp.mpf(10) ** -40 for r, s in zip(recs, st))
    objs = [Q_objects(r, ZETA) for r in recs]
    tdP = [q.SiteData(o[0], th, 8) for o in objs]
    tdQ = [q.SiteData(o[1], th, 8) for o in objs]
    return ch, st, th, recs, tdP, tdQ

def R_compute(st, th, tdP, tdQ):
    out = {}
    for A in range(len(st)):
        for B in range(len(st)):
            if A == B:
                continue
            N = Nrows(tdP[B], tdQ[A], th, 2)
            keys = sorted(N)
            M = mp.matrix(18, 6)
            for r, key in enumerate(keys):
                for c in range(6):
                    M[r, c] = N[key][c]
            colscale = []
            for c in range(6):
                mx = max(abs(M[i, c]) for i in range(18))
                colscale.append(mx)
                if mx != 0:
                    for i in range(18):
                        M[i, c] /= mx
            U, S, V = mp.svd_r(M)
            out[(A, B)] = dict(S=[S[i] / S[0] for i in range(6)], V=[[V[i, j] for j in range(6)] for i in range(6)], colscale=colscale,
                               D=q.delta_vec(st[B], st[A]), M=[[M[i, j] for j in range(6)] for i in range(18)])
    return out

def sectionR():
    hdr('SECTION R: rank of the stacked 18x6 matrix for all ordered pairs A != B of Lambda = (2,1,0) (64 states, gauge p = 2)')
    from collections import Counter
    ch, st, th, recs, tdP, tdQ = R_setup()
    t0 = time.time()
    data = cache('qside_320_R.pkl', lambda: R_compute(st, th, tdP, tdQ))
    print('  computed/loaded %d ordered pairs (%.0fs)' % (len(data), time.time() - t0))
    thr = mp.mpf(10) ** -40
    mus = [s['sector'] for s in st]
    for key, d in data.items():
        d['rank'] = ranks_from_sv(d['S'], thr)
    print('  overall rank tally:', dict(sorted(Counter(d['rank'] for d in data.values()).items())))
    tab = {}
    for (A, B), d in data.items():
        tab.setdefault((mus[A], mus[B]), Counter())[d['rank']] += 1
    same = Counter()
    for (ma, mb), c in tab.items():
        if ma == mb:
            same.update(c)
    print('  same-sector pairs rank tally:', dict(same))
    low = {k: v for k, v in tab.items() if any(r < 5 for r in v)}
    print('  ordered sector pairs: %d total; with some rank<5 pair: %d; with ALL pairs rank<5: %d' % (len(tab), len(low), sum(1 for v in low.values() if all(r < 5 for r in v))))
    print('  rank<5 sector pairs (A weight, B weight, mu_B - mu_A, rank counts):')
    for k in sorted(low):
        print('     ', k[0], k[1], tuple(b - a for a, b in zip(*k)), dict(low[k]))
    # characterisation: rank is a function of Dmu = mu_B - mu_A only, via the number of NEGATIVE components
    bydm = {}
    for (A, B), d in data.items():
        dm = tuple(b - a for a, b in zip(mus[A], mus[B]))
        bydm.setdefault(dm, set()).add(d['rank'])
    print('  rank depends on Dmu = mu_B - mu_A only (every Dmu gives a single rank):', all(len(v) == 1 for v in bydm.values()), '; number of distinct Dmu:', len(bydm))
    tally = {}
    for dm, v in bydm.items():
        nneg = sum(1 for x in dm if x < 0)
        npos = sum(1 for x in dm if x > 0)
        tally.setdefault((npos, nneg), set()).update(v)
    print('  (number of positive Dmu components, number of negative components) -> ranks seen:', dict(sorted(tally.items())))
    print('  => rank 4  <=>  Dmu has exactly one positive and two negative components (sign pattern (+,-,-) in some order); otherwise rank 5')
    # mechanism: for each k with Dmu_k <= -1 the rational function psi_k P^k (exponentials cancel) has degree Dmu_k <= -1, so the sum of its residues
    # (with the weights 1/(u-th_t)^m of the six columns) vanishes: row(1,k,k) + row(2,k,k) = 0
    worst = mp.mpf(0)
    n_checked = 0
    worst_pos = mp.mpf(1)
    for (A, B), d in data.items():
        dm = tuple(b - a for a, b in zip(mus[A], mus[B]))
        M = d['M']
        for k in range(3):
            r1, r2 = M[4 * k], M[9 + 4 * k]
            nrm = max(abs(x) for x in r1 + r2)
            dev = max(abs(x + y) for x, y in zip(r1, r2)) / nrm if nrm > 0 else mp.mpf(0)
            if dm[k] <= -1:
                worst = max(worst, dev)
                n_checked += 1
            else:
                worst_pos = min(worst_pos, dev)
    print('  row relation row(1,k,k) + row(2,k,k) = 0 (all six columns) for every k with Dmu_k <= -1 (%d (pair,k) cases): max relative violation %s;' % (n_checked, nst(worst)),
          'for k with Dmu_k >= 0 the smallest violation is %s (relation does NOT hold)' % nst(worst_pos))
    # per-site row space dimension
    def rk_rows(rows):
        Mm = mp.matrix(rows)
        for c in range(Mm.cols):
            mx = max(abs(Mm[i, c]) for i in range(Mm.rows)) or 1
            for i in range(Mm.rows):
                Mm[i, c] /= mx
        Sg = mp.svd_r(Mm, compute_uv=False)
        return sum(1 for x in Sg if x / Sg[0] > thr)
    sdim = {}
    for key in list(data)[:400]:
        M = data[key]['M']
        sdim[(rk_rows(M[:9]), rk_rows(M[9:]), data[key]['rank'])] = sdim.get((rk_rows(M[:9]), rk_rows(M[9:]), data[key]['rank']), 0) + 1
    print('  (rank of site-1 rows, rank of site-2 rows, total rank) over 400 pairs ->', sdim, ' [each site contributes a 3-dim row space S_s; total rank = 6 - dim(S_1 n S_2)]')
    return ch, st, th, recs, tdP, tdQ, data

def main(which):
    t0 = time.time()
    ch = Chain(lam=LAM)
    if '0' in which:
        section0(ch)
        vac_closed_forms(ch)
        pass  # section0_other_reps() already done by P3 (c2 scan incl. (3,2,0) = 8)
        section_vacuum(ch)
    if {'1', '2', '3'} & set(which):
        t1 = time.time()
        def mk():
            ch_, ops, C, st, res = build_states()
            return st, res
        st, res = cache('qside_320_states.pkl', mk)
        ch = Chain(lam=LAM)
        th = th_mp(ch)
        print('\n%d eigenstates (%d sectors), %.0fs; max eigen-residual of all coefficient operators on all states: %s' % (len(st), len(ch.sectors), time.time() - t1, nst(res)))
        q.GAUGE = PG
        recs = section1(ch, st, th, res)
        section1_neg(ch, st, th)
        if '2' in which:
            section2(ch, st, th, recs)
        if '3' in which:
            section3(ch, st, th, recs)
    if '4' in which:
        section4()
    if 'R' in which:
        sectionR()
    print('\ntotal time %.0fs' % (time.time() - t0))

if __name__ == '__main__':
    which = sys.argv[1:] or ['0', '1', '2', '3', '4']
    main(which)
