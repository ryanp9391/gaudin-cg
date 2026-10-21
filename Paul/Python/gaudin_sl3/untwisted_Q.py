"""U1: the untwisted (zeta = 0) Q-system of the sl(3) Gaudin chain, L = 2, and the eps-expansion of the HW-limit branches.
Stages: python -X utf8 untwisted_Q.py --stages T0,T1,T2,T3   (see RESULTS-untwisted-Q.md).  Read-only reuse of qside / gaudin_chain / twovac_fit.
"""
import sys, os, time, itertools, pickle, argparse
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mpmath as mp
import sympy
from gaudin_chain import Chain, Operators, SM

def hdr(t):
    print('\n' + '=' * 100 + '\n' + t + '\n' + '=' * 100, flush=True)

def C2(l):
    return sum(l[i] * (l[i] + 4 - 2 * (i + 1)) for i in range(3))

def S(x):
    return sympy.Rational(x.numerator, x.denominator)

# --------------------------------------------------------------------------------------------- T0 : exact untwisted data
def hw_spaces(ch):
    """exact basis (Fractions, full length) of the total HW vectors, per weight mu"""
    out = {}
    for mu, idx in ch.sectors.items():
        E12 = ch.Es[(1, 1, 2)] + ch.Es[(2, 1, 2)]
        E23 = ch.Es[(1, 2, 3)] + ch.Es[(2, 2, 3)]
        rows = []
        for E in (E12, E23):
            for r, row in E.r.items():
                rows.append([S(row[j]) if j in row else sympy.Integer(0) for j in idx])
        if rows:
            ns = sympy.Matrix(rows).nullspace()
        else:
            ns = [sympy.eye(len(idx))[:, k] for k in range(len(idx))]
        if ns:
            out[mu] = [[(Fr(int(sympy.Rational(v[k]).p), int(sympy.Rational(v[k]).q)) if True else 0) for k in range(len(idx))] for v in ns]
    return out

def op_matrix_on_hw(op, ch, mu, basis):
    """matrix of op (SM, commuting with total gl3) restricted to the HW space of weight mu, in the given basis"""
    idx = ch.sectors[mu]
    n = len(basis)
    vecs = []
    for b in basis:
        full = [Fr(0)] * ch.dim
        for k, j in enumerate(idx):
            full[j] = b[k]
        vecs.append(full)
    B = sympy.Matrix([[S(x) for x in b] for b in basis]).T          # columns = basis (in sector coordinates)
    cols = []
    for full in vecs:
        w = op.apply(full)
        cols.append([S(w[j]) for j in idx])
        # w must lie in the sector
        assert all(w[j] == 0 for j in range(ch.dim) if j not in set(idx))
    W = sympy.Matrix(cols).T
    # solve B X = W
    X = (B.T * B).inv() * B.T * W
    assert B * X == W, 'op does not preserve the HW space'
    return X

def t0_data(lam, theta):
    """for one rep and theta: dict mu -> dict(mult, H1,a1,b1,H2,a2,b2 as sympy matrices on the HW space)"""
    ch = Chain(lam, theta)
    g = [[Fr(0)] * 3 for _ in range(3)]
    ops = Operators(ch, g)
    C = ops.fit()
    assert C['res2'] == 0 and C['res3'] == 0
    hw = hw_spaces(ch)
    out = {}
    names = {'H1': (2, 1), 'H2': (2, 3), 'b1': (3, 1), 'a1': (3, 2), 'b2': (3, 4), 'a2': (3, 5)}
    triple = [C[3][3], C[3][6]]
    assert all(t.is_zero() for t in triple)
    for mu, basis in hw.items():
        d = dict(mult=len(basis), basis=basis)
        for nm, (k, b) in names.items():
            d[nm] = op_matrix_on_hw(C[k][b], ch, mu, basis)
        out[mu] = d
    return ch, out

def simdiag(d):
    """simultaneous eigen-data on the HW space: returns list of dict(H1,a1,b1,H2,a2,b2) (sympy numbers) for multiplicity m (eigenvalues as sympy exprs)"""
    m = d['mult']
    if m == 1:
        return [{k: sympy.nsimplify(d[k][0, 0]) for k in ('H1', 'a1', 'b1', 'H2', 'a2', 'b2')}]
    keys = ('H1', 'a1', 'b1', 'H2', 'a2', 'b2')
    gen = d['a1'] + 3 * d['b1'] + 7 * d['a2'] + 11 * d['b2'] + 13 * d['H1']
    ev = gen.eigenvects()
    res = []
    for val, mult_, vecs in ev:
        for v in vecs:
            r = {}
            for k in keys:
                w = d[k] * v
                # eigenvalue
                i = [t for t in range(m) if v[t] != 0][0]
                r[k] = sympy.simplify(w[i] / v[i])
            res.append(r)
    return res


# --------------------------------------------------------------------------------------------- regime machinery (eps-branches)
LAM = (2, 1, 0)
THETA = (Fr(2, 7), Fr(-3, 5))
PHI = (Fr(1, 3), Fr(-2, 7), Fr(5, 11))
PHI2 = (Fr(-2, 5), Fr(3, 8), Fr(7, 11))
HW_SECT = {(4, 2, 0): 1, (4, 1, 1): 1, (3, 3, 0): 1, (3, 2, 1): 2, (2, 2, 2): 1}
BRANCHES = [((4, 2, 0), 0), ((4, 1, 1), 0), ((3, 3, 0), 0), ((3, 2, 1), 0), ((3, 2, 1), 1), ((2, 2, 2), 0)]
NEPS = 12
def eps_k(k):
    return Fr(1, 10 ** 6 * 2 ** k)

def mpf(x):
    if isinstance(x, Fr):
        return mp.mpf(x.numerator) / x.denominator
    return mp.mpf(x)

def real(x, tol):
    if isinstance(x, mp.mpc):
        assert abs(x.imag) < tol * max(1, abs(x.real)), x
        return x.real
    return x

def raising_residual(ch, s):
    v = s['vec']
    nv = mp.sqrt(sum(abs(x) ** 2 for x in v))
    tot = 0
    for (i, j) in ((1, 2), (2, 3)):
        E = ch.Es[(1, i, j)] + ch.Es[(2, i, j)]
        w = E.apply(v)
        tot += sum(abs(x) ** 2 for x in w)
    return mp.sqrt(tot) / nv

# polynomial helpers on mp lists (index = power)
def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]
def pmul(a, b):
    out = [mp.mpf(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out
def pder(a):
    return [i * a[i] for i in range(1, len(a))] or [mp.mpf(0)]
def pscale(a, c):
    return [c * x for x in a]
def peval(a, x):
    r = 0
    for c in reversed(a):
        r = r * x + c
    return r

def oper_coeffs(data, theta, zeta=(0, 0, 0)):
    """data = (H1,a1,b1,H2,a2,b2); returns coefficient lists of Pi^2 S1, Pi^2 S2, Pi^2 S3, Pi=(u-th1)(u-th2),
    S1 = e1 + 3/x1 + 3/x2, S2 = e2 + sum H/x + 3/x^2, S3 = e3 + sum b/x + a/x^2   (x_s = u - theta_s)"""
    H1, a1, b1, H2, a2, b2 = data
    t1, t2 = [mpf(t) for t in theta]
    z = [mpf(x) for x in zeta]
    e1 = z[0] + z[1] + z[2]; e2 = z[0] * z[1] + z[0] * z[2] + z[1] * z[2]; e3 = z[0] * z[1] * z[2]
    x1 = [-t1, mp.mpf(1)]; x2 = [-t2, mp.mpf(1)]
    Pi = pmul(x1, x2); Pi2 = pmul(Pi, Pi)
    def pf(c0, s1, s2):
        (A1, B1), (A2, B2) = s1, s2
        out = pscale(Pi2, c0)
        out = padd(out, pscale(pmul(x1, pmul(x2, x2)), A1))
        out = padd(out, pscale(pmul(x2, pmul(x1, x1)), A2))
        out = padd(out, pscale(pmul(x2, x2), B1))
        out = padd(out, pscale(pmul(x1, x1), B2))
        return out
    return pf(e1, (3, 0), (3, 0)), pf(e2, (H1, 3), (H2, 3)), pf(e3, (b1, a1), (b2, a2))

def dual_poly(cf, Pi2, kappa, d):
    """monic polynomial p (deg d) with D(e^{kappa u} p) = 0, i.e. Pi^2[(d+k)^3 - S1 (d+k)^2 + S2 (d+k) - S3] p = 0 (QR least squares)"""
    PS1, PS2, PS3 = cf
    cols = []
    for n in range(d + 1):
        D0 = [mp.mpf(0)] * n + [mp.mpf(1)]
        D1 = padd(pder(D0), pscale(D0, kappa))
        D2 = padd(pder(D1), pscale(D1, kappa))
        D3 = padd(pder(D2), pscale(D2, kappa))
        tot = pmul(Pi2, D3)
        tot = padd(tot, pscale(pmul(PS1, D2), -1))
        tot = padd(tot, pmul(PS2, D1))
        tot = padd(tot, pscale(pmul(PS3, D0), -1))
        cols.append(tot)
    R = max(len(c) for c in cols)
    M = mp.matrix(R, d + 1)
    for n, c in enumerate(cols):
        for m, v in enumerate(c):
            M[m, n] = v
    if d == 0:
        return [mp.mpf(1)], max(abs(M[m, 0]) for m in range(R))
    A = M[:, :d]
    b = -M[:, d]
    sc = [max(abs(A[m, n]) for m in range(R)) or mp.mpf(1) for n in range(d)]
    for n in range(d):
        for m in range(R):
            A[m, n] /= sc[n]
    x, resn = mp.qr_solve(A, b)
    p = [x[n] / sc[n] for n in range(d)] + [mp.mpf(1)]
    return p, resn

def pi2(theta):
    t1, t2 = [mpf(t) for t in theta]
    Pi = pmul([-t1, mp.mpf(1)], [-t2, mp.mpf(1)])
    return pmul(Pi, Pi)

def build_regime(eps, dps=200, phi=PHI, lam=LAM, theta=THETA, tag='', verbose=True):
    """eigenstates at zeta = eps*phi, the six HW-limit branches, Laurent data and dual monic polynomials p_j.  Cached per (eps, dps, tag)."""
    path = os.path.join(HERE, 'untwisted_Q_reg%s_%s_%d.pkl' % (tag, str(eps).replace('/', '_'), dps))
    if os.path.exists(path):
        with open(path, 'rb') as f:
            out = pickle.load(f)
        mp.mp.dps = dps
        return out
    t0 = time.time()
    mp.mp.dps = dps
    zeta = tuple(eps * p for p in phi)
    ch = Chain(lam, theta)
    from gaudin_chain import eigenstates, eigenfunctions
    g = [[zeta[i] if i == j else Fr(0) for j in range(3)] for i in range(3)]
    ops = Operators(ch, g)
    C = ops.fit()
    st = eigenstates(ch, ops, C, dps=dps)
    tol = mp.mpf(10) ** -(dps - 30)
    for x in st:
        x['ev'] = real(x['ev'], tol)
        x['vec'] = [real(c, tol) for c in x['vec']]
    res = eigenfunctions(ch, C, st)
    bysec = {}
    for k, s in enumerate(st):
        s['rr'] = raising_residual(ch, s)
        bysec.setdefault(s['sector'], []).append(k)
    PI2 = pi2(theta)
    recs = []
    for mu, m in HW_SECT.items():
        ks = sorted(bysec[mu], key=lambda k: st[k]['rr'])
        sel = sorted(ks[:m], key=lambda k: st[k]['ev'])
        gap = (st[ks[m - 1]]['rr'], st[ks[m]]['rr'] if len(ks) > m else mp.mpf(0))
        for c, k in enumerate(sel):
            s = st[k]
            l2 = [real(x, tol) for x in s['lam2']]; l3 = [real(x, tol) for x in s['lam3']]
            data = (l2[1], l3[2], l3[1], l2[3], l3[5], l3[4])
            cf = oper_coeffs(data, theta, zeta)
            ps = []; rs = []
            for j in range(3):
                p, r = dual_poly(cf, PI2, mpf(zeta[j]), mu[j])
                ps.append(p); rs.append(r)
            recs.append(dict(mu=mu, copy=c, ev=s['ev'], rr=s['rr'], gap=gap, data=data, lam2=l2, lam3=l3, p=ps, res=rs, vec=s['vec']))
    out = dict(eps=eps, zeta=zeta, dps=dps, eigres=res, branches=recs, time=time.time() - t0)
    with open(path, 'wb') as f:
        pickle.dump(out, f)
    if verbose:
        print('  built eps=%s (%.0fs) eigres %s' % (eps, time.time() - t0, mp.nstr(res, 3)), flush=True)
    return out


# --------------------------------------------------------------------------------------------- jets, Plucker monomials
PAIRS = [(0, 2), (1, 2), (0, 3), (1, 3), (2, 3)]       # the five adjoint pairs (kl) (and p01 = 0)

def poly_derivs_exp(kappa, q, n):
    """polynomials r_k with d^k (e^{kappa u} q) = e^{kappa u} r_k, k = 0..n"""
    out = [list(q)]
    cur = list(q)
    for _ in range(n):
        cur = padd(pder(cur), pscale(cur, kappa))
        out.append(cur)
    return out

def jets_at(kappa, q, x, n=4):
    e = mp.exp(kappa * x)
    return [e * peval(r, x) for r in poly_derivs_exp(kappa, q, n)]

def branch_jets(rec, zeta, theta, n=4):
    """J[s][j][k] = (d/du)^k P^j at theta_s"""
    z = [mpf(x) for x in zeta]
    t = [mpf(x) for x in theta]
    return [[jets_at(z[j], rec['p'][j], t[s], n) for j in range(3)] for s in range(2)]

def plucker_all(J, a=1, b=2, order=3):
    return {(k, l): J[a][k] * J[b][l] - J[b][k] * J[a][l] for k in range(order + 1) for l in range(k + 1, order + 1)}

def monomials(rec, zeta, theta):
    """per site: dict with p_kl (all k<l<=3), N_i = P3 p_kl for the 5 pairs, M = P3 p_02, kappa_s etc."""
    J = branch_jets(rec, zeta, theta, 4)
    out = []
    for s in range(2):
        pl = plucker_all(J[s])
        P3 = J[s][2][0]
        N = [P3 * pl[pr] for pr in PAIRS]
        out.append(dict(p=pl, P3=P3, N=N, M=P3 * pl[(0, 2)], p01=pl[(0, 1)], J=J[s]))
    return out


# --------------------------------------------------------------------------------------------- eps-series fits
def series_fit(eps_list, vals, order=None):
    """coefficients c_0..c_{K-1} of the polynomial interpolating vals at eps_list (exact Lagrange/Vandermonde solve at working precision);
    the variable is scaled x = eps/eps_0 for conditioning (coefficients returned for eps)."""
    K = len(eps_list)
    e0 = mpf(eps_list[0])
    xs = [mpf(e) / e0 for e in eps_list]
    A = mp.matrix(K, K)
    for i in range(K):
        for j in range(K):
            A[i, j] = xs[i] ** j
    c = mp.lu_solve(A, mp.matrix(vals))
    return [c[j] / e0 ** j for j in range(K)]

def fit_stability(eps_list, vals, nmax=10):
    """fit with all K points, and with the first K-1 and last K-1 points: returns (coeffs_all, [rel. difference per order])"""
    full = series_fit(eps_list, vals)
    a = series_fit(eps_list[:-1], vals[:-1])
    b = series_fit(eps_list[1:], vals[1:])
    rel = []
    for n in range(min(nmax + 1, len(a))):
        ref = abs(full[n]) or mp.mpf(1)
        rel.append(max(abs(full[n] - a[n]), abs(full[n] - b[n])) / ref if full[n] != 0 else max(abs(a[n]), abs(b[n])))
    return full, rel


# --------------------------------------------------------------------------------------------- no-log (Gaudin Bethe) conditions
_FROB = {}
def frobenius_conditions():
    """symbolic no-log conditions at theta_s for S1 = 3/x + sum s_k x^k, S2 = 3/x^2 + H/x + sum h_k x^k, S3 = a/x^2 + b/x + sum t_k x^k
    (local exponents 0,2,4).  Returns (vars, [(order, free-constant index, polynomial)])."""
    if 'p' in _FROB:
        return _FROB['p']
    x = sympy.Symbol('x')
    H, a, b = sympy.symbols('H a b')
    N = 6
    sv = sympy.symbols('s0:%d' % N); hv = sympy.symbols('h0:%d' % N); tv = sympy.symbols('t0:%d' % N)
    S1 = 3 / x + sum(sv[k] * x ** k for k in range(N))
    S2 = 3 / x ** 2 + H / x + sum(hv[k] * x ** k for k in range(N))
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
            C.append((n, mon.index(1), sympy.expand(co)))
    vars_ = (H, a, b) + tuple(sv[:3]) + tuple(hv[:3]) + tuple(tv[:3])
    _FROB['p'] = (vars_, C)
    return _FROB['p']

def nolog_values(data, theta, zeta):
    """data = (H1,a1,b1,H2,a2,b2) numbers/sympy; theta, zeta rational/sympy; returns the list of no-log conditions (both sites)"""
    vars_, C = frobenius_conditions()
    z = [sympy.nsimplify(q) for q in zeta]
    e1 = sum(z); e2 = z[0] * z[1] + z[0] * z[2] + z[1] * z[2]; e3 = z[0] * z[1] * z[2]
    out = []
    for s in range(2):
        o = 1 - s
        d = sympy.nsimplify(theta[s]) - sympy.nsimplify(theta[o])
        def reg(k, m):
            return (-1) ** k * sympy.binomial(m + k - 1, k) * d ** (-m - k)
        H, a, b = data[3 * s:3 * s + 3]
        Ho, ao, bo = data[3 * o:3 * o + 3]
        sk = [3 * reg(k, 1) + (e1 if k == 0 else 0) for k in range(3)]
        hk = [Ho * reg(k, 1) + 3 * reg(k, 2) + (e2 if k == 0 else 0) for k in range(3)]
        tk = [bo * reg(k, 1) + ao * reg(k, 2) + (e3 if k == 0 else 0) for k in range(3)]
        vals = [H, a, b] + sk + hk + tk
        sub = dict(zip(vars_, vals))
        for (n, idx, poly) in C:
            out.append(sympy.expand(poly.subs(sub)))
    return out



# --------------------------------------------------------------------------------------------- exact perturbation theory of the Laurent data
def zeta0_data(mu, copy):
    """exact zeta=0 Laurent data (H1,a1,b1,H2,a2,b2) in units theta1 - theta2 = 1 (theta=(1,0)), as sympy numbers"""
    ch, out = t0_data(LAM, (Fr(1), Fr(0)))
    ev = simdiag(out[mu])
    if len(ev) > 1:
        # copy 0 = the one with b1 > 0 (matches the regime branch ordering by eigenvalue)
        ev = sorted(ev, key=lambda e: -float(sympy.N(e['b1'])))
    e = ev[copy]
    return [sympy.nsimplify(e[k]) if not e[k].has(sympy.sqrt(5)) else sympy.simplify(e[k]) for k in ('H1', 'a1', 'b1', 'H2', 'a2', 'b2')]

def perturb_data(mu, copy, phi, N, e=None):
    """order-by-order exact solution of the six no-log conditions around the zeta = 0 solution, zeta = e*phi (phi sympy/rational), units theta=(1,0).
    returns list X_0..X_N of 6-vectors (sympy)."""
    e = e or sympy.Symbol('e')
    th = (sympy.Integer(1), sympy.Integer(0))
    X = [zeta0_data(mu, copy)]
    zt = [e * sympy.sympify(p) for p in phi]
    syms = sympy.symbols('X0:6')
    mu_ = [sympy.Integer(m) for m in mu]
    def extra(data, ee):
        H1, a1, b1, H2, a2, b2 = data
        p = [sympy.sympify(q) for q in phi]
        g1 = H1 + H2 - ee * (p[0] * (mu_[1] + mu_[2]) + p[1] * (mu_[0] + mu_[2]) + p[2] * (mu_[0] + mu_[1]))
        g2 = b1 + b2 - ee ** 2 * (p[0] * p[1] * mu_[2] + p[0] * p[2] * mu_[1] + p[1] * p[2] * mu_[0])
        return [g1, g2]
    def F_of(data):
        return nolog_values(data, th, zt) + extra(data, e)
    # Jacobian at zeta = 0
    F0 = F_of(list(syms[:]))
    F0 = [f.subs(dict(zip(syms, [0] * 6))) if False else f for f in F0]
    zt0 = [0, 0, 0]
    Fz = nolog_values(list(syms), th, zt0) + extra(list(syms), 0)
    sol0 = dict(zip(syms, X[0]))
    assert all(sympy.simplify(f.subs(sol0)) == 0 for f in Fz), 'zeta=0 data does not satisfy the no-log conditions'
    J = sympy.Matrix([[sympy.diff(f, s).subs(sol0) for s in syms] for f in Fz])
    J = J.applyfunc(sympy.simplify)
    JtJ = (J.T * J).applyfunc(sympy.simplify)
    if JtJ.det() == 0 or sympy.simplify(JtJ.det()) == 0:
        raise RuntimeError('rank-deficient Jacobian at zeta = 0')
    Jinv = (JtJ.inv() * J.T).applyfunc(sympy.simplify)
    Fe = F_of(list(syms))
    for n in range(1, N + 1):
        cur = [sum(X[k][i] * e ** k for k in range(n)) for i in range(6)]
        sub = dict(zip(syms, cur))
        rhs = []
        for f in Fe:
            ex = sympy.expand(f.subs(sub))
            rhs.append(ex.coeff(e, n))
        Xn = -(Jinv * sympy.Matrix(rhs))
        X.append([sympy.simplify(x) for x in Xn])
    return X


def perturb_data_deg(mu, copy, phi, N, e=None, verbose=False):
    """as perturb_data, but for a rank-deficient Jacobian (kernel dimension 1): the kernel component c_n of X_n is fixed by the
    solvability condition (left null vector) at order n+1 of the no-log + sum-rule system."""
    e = e or sympy.Symbol('e')
    th = (sympy.Integer(1), sympy.Integer(0))
    mu_ = [sympy.Integer(m) for m in mu]
    p = [sympy.sympify(q) for q in phi]
    zt = [e * q for q in p]
    syms = sympy.symbols('X0:6')
    def extra(data, ee):
        H1, a1, b1, H2, a2, b2 = data
        g1 = H1 + H2 - ee * (p[0] * (mu_[1] + mu_[2]) + p[1] * (mu_[0] + mu_[2]) + p[2] * (mu_[0] + mu_[1]))
        g2 = b1 + b2 - ee ** 2 * (p[0] * p[1] * mu_[2] + p[0] * p[2] * mu_[1] + p[1] * p[2] * mu_[0])
        return [g1, g2]
    X0 = zeta0_data(mu, copy)
    sol0 = dict(zip(syms, X0))
    Fz = nolog_values(list(syms), th, [0, 0, 0]) + extra(list(syms), 0)
    J = sympy.Matrix([[sympy.diff(f, s).subs(sol0) for s in syms] for f in Fz]).applyfunc(sympy.simplify)
    ker = J.nullspace()
    assert len(ker) == 1
    K = ker[0].applyfunc(sympy.simplify)
    left = J.T.nullspace()
    Fe = nolog_values(list(syms), th, zt) + extra(list(syms), e)
    X = [list(X0)]
    cs = []
    def Rn(Xl, n):
        cur = [sum(Xl[k][i] * e ** k for k in range(len(Xl))) for i in range(6)]
        sub = dict(zip(syms, cur))
        return [sympy.expand(sympy.expand(f.subs(sub)).coeff(e, n)) for f in Fe]
    # unknown parts: X_n = Xp_n + c_n K
    Xp = [list(X0)]
    csym = []
    for n in range(1, N + 2):
        # solvability at order n determines c_{n-1} (for n >= 2); then solve for X_n particular
        R = Rn(X, n)
        if n >= 2:
            eqs = []
            for l in left:
                eqs.append(sympy.simplify(sum(l[i] * R[i] for i in range(len(R)))))
            c = csym[-1]
            eqs = [q for q in eqs if q != 0]
            if eqs:
                s = sympy.solve(eqs[0], c, dict=True)
                if verbose: print('n', n, 'solvability eq', eqs[0], '->', s)
                # prefer the solution consistent with the numerics: take all, caller picks (store)
                s = s[0] if len(s) == 1 else s[-1]
                X = [[sympy.simplify(x.subs(s)) for x in v] for v in X]
                R = Rn(X, n)
        if n == N + 1:
            break
        Rm = sympy.Matrix(R)
        part, params = J.gauss_jordan_solve(-Rm)
        part = part.subs({q: 0 for q in params})
        c = sympy.Symbol('c%d' % n)
        csym.append(c)
        X.append([sympy.simplify(part[i] + c * K[i]) for i in range(6)])
    return X[:N + 1], csym

# --------------------------------------------------------------------------------------------- stages
V0_BASES = None
def v0_bases():
    """exact echelon bases of the untwisted solution space V0 (polynomials) in x = (u - theta2)/(theta1 - theta2), per branch (mu, copy)"""
    global V0_BASES
    if V0_BASES is None:
        xs = sympy.Symbol('x'); s5 = sympy.sqrt(5)
        V0_BASES = {((4, 2, 0), 0): [1, xs ** 2 * (2 * xs - 3) / 2, xs ** 4 * (2 * xs - 3) ** 2 / 4],
                    ((4, 1, 1), 0): [(2 * xs - 1) / 2, xs ** 2, xs ** 4 * (2 * xs ** 2 - 6 * xs + 5) / 2],
                    ((3, 3, 0), 0): [1, xs ** 2 * (xs - 1) ** 2, xs ** 2 * (4 * xs ** 3 - 10 * xs + 5) / 4],
                    ((3, 2, 1), 0): [(10 * xs - 5 + s5) / 10, xs ** 2 * (2 * xs - 3 - s5) / 2, xs ** 4 * (2 * xs - 5 + s5) / 2],
                    ((3, 2, 1), 1): [(10 * xs - 5 - s5) / 10, xs ** 2 * (2 * xs - 3 + s5) / 2, xs ** 4 * (2 * xs - 5 - s5) / 2],
                    ((2, 2, 2), 0): [(2 * xs - 1) ** 2 / 4, (2 * xs - 1) * (4 * xs ** 2 + 2 * xs - 1) / 8, xs ** 4]}
    return V0_BASES

def casimirs_total(ch):
    Et = {(i, j): ch.Es[(1, i, j)] + ch.Es[(2, i, j)] for i in (1, 2, 3) for j in (1, 2, 3)}
    C2o = SM(ch.dim); C3o = SM(ch.dim)
    for i in (1, 2, 3):
        for j in (1, 2, 3):
            C2o = C2o + Et[(i, j)] * Et[(j, i)]
            for k in (1, 2, 3):
                C3o = C3o + Et[(i, j)] * Et[(j, k)] * Et[(k, i)]
    return C2o, C3o

def eig_on(op, ch, mu, v):
    idx = ch.sectors[mu]
    full = [Fr(0)] * ch.dim
    for k, j in enumerate(idx):
        full[j] = v[k]
    w = op.apply(full)
    i = [k for k in range(len(idx)) if v[k] != 0][0]
    lam = w[idx[i]] / v[i]
    assert all(w[j] == lam * full[j] for j in range(ch.dim))
    return lam

def stage_T0():
    hdr('T0: untwisted (zeta = 0) Laurent data, exact; H1 = (|L|^2 - (C2(mu) - 2 C2(L))/2)/theta12; (a1-a2)/2 = (mu1+2)(mu2+1)mu3/(2 theta12)')
    for lam in [(2, 1, 0), (1, 1, 0), (2, 2, 0), (3, 1, 0), (3, 2, 0)]:
        for th in [(Fr(3), Fr(0)), (Fr(2, 7), Fr(-3, 5))]:
            ch, out = t0_data(lam, th)
            t12 = S(th[0] - th[1])
            g = [[Fr(0)] * 3 for _ in range(3)]
            Cf = Operators(ch, g).fit()
            glob = ((Cf[3][1] + Cf[3][4]).is_zero(), (Cf[3][2] + Cf[3][5] + Cf[3][1] * th[0] + Cf[3][4] * th[1]).is_zero(), (Cf[2][1] + Cf[2][3]).is_zero())
            C2o, C3o = casimirs_total(ch)
            L = sum(lam)
            print('Lambda %s theta %s: global operator identities b1+b2=0, a1+a2+th1 b1+th2 b2=0, H1+H2=0: %s' % (lam, tuple(map(str, th)), glob))
            for mu in sorted(out, reverse=True):
                d = out[mu]
                ev = simdiag(d)
                H1 = ev[0]['H1'] * t12
                pred = L ** 2 - sympy.Rational(1, 2) * (C2(mu) - 2 * C2(lam))
                q = [sympy.simplify((e['a1'] - e['a2']) / 2 * t12) for e in ev]
                qp = sympy.Rational((mu[0] + 2) * (mu[1] + 1) * mu[2], 2)
                s = [sympy.simplify((e['a1'] + e['a2']) * t12) for e in ev]
                b = [sympy.simplify(e['b1'] * t12 ** 2) for e in ev]
                c2 = eig_on(C2o, ch, mu, d['basis'][0]); c3 = eig_on(C3o, ch, mu, d['basis'][0])
                print('  mu=%s mult=%d C2=%s C3=%s | H1 th12 = %s (formula %s) | (a1-a2)/2 th12 = %s (formula %s) | (a1+a2) th12 = %s | b1 th12^2 = %s' % (
                    mu, d['mult'], c2, c3, H1, pred, q, qp, s, b))
    print('  b1 th12^2 squared (multiplicity-2 sectors):  beta = b1^2 th12^4: Lambda=210 mu=321: 45; Lambda=310: (5,2,1): 105, (4,3,1): 192; Lambda=320: (5,4,1): 105, (5,3,2): 192')
    print('  eigenvalues of a1 th12 on the multiplicity space: (a1-a2)/2 th12 +- sqrt(beta)/2 ; a2 = -a1 - b1 th12 ; b2 = -b1 ; H2 = -H1')

def stage_T1():
    hdr('T1: the untwisted oper D0 = d^3 - S1 d^2 + S2 d - S3 (adjoint, units theta1 - theta2 = 1, x = u - theta2)')
    xs = sympy.Symbol('x')
    ch, out = t0_data(LAM, (Fr(1), Fr(0)))
    k = sympy.Symbol('k')
    for (mu, cp), basis in v0_bases().items():
        ev = simdiag(out[mu])
        if len(ev) > 1:
            ev = sorted(ev, key=lambda e: -float(sympy.N(e['b1'])))
        e = ev[cp]
        x1 = xs - 1; x2 = xs
        S1 = 3 / x1 + 3 / x2
        S2 = e['H1'] / x1 + e['H2'] / x2 + 3 / x1 ** 2 + 3 / x2 ** 2
        S3 = e['b1'] / x1 + e['b2'] / x2 + e['a1'] / x1 ** 2 + e['a2'] / x2 ** 2
        D = lambda f: sympy.diff(f, xs, 3) - S1 * sympy.diff(f, xs, 2) + S2 * sympy.diff(f, xs) - S3 * f
        ok = [sympy.simplify(D(f)) == 0 for f in basis]
        W = sympy.factor(sympy.simplify(sympy.Matrix([[sympy.diff(f, xs, m) for m in range(3)] for f in basis]).det()))
        l = (mu[2], mu[1] + 1, mu[0] + 2)
        # exponent at infinity: u^k ; k^3 - 9 k^2 + (14 + h) k - 2 q = 0
        h = e['H1']; q = (e['a1'] - e['a2']) / 2
        roots = sympy.roots(sympy.Poly(k ** 3 - 9 * k ** 2 + (14 + h) * k - 2 * q, k))
        print('  branch %s copy %d: H1=%s a1=%s b1=%s | exponents at theta_s {0,2,4}; at infinity (growth u^k) k in %s (= (mu3, mu2+1, mu1+2) = %s); Fuchs: 12 - sum k = 3' % (mu, cp, h, e['a1'], e['b1'], sorted(roots), l))
        print('     basis of V0 (monic, echelon degrees %s): %s ; all solve D0: %s ; Wronskian = %s' % (l, [sympy.factor(f) for f in basis], ok, W))

def stage_T2(dps=200, K=12, phi=PHI, tag='', nmax=10):
    """series fits + export for U2"""
    hdr('T2: eps-series (dps %d, K = %d points eps_k = 1e-6 2^-k) of the six branches, phi = %s' % (dps, K, tuple(map(str, phi))))
    mp.mp.dps = dps
    eps = [eps_k(k) for k in range(K)]
    outs = [build_regime(eps_k(k), dps=dps, phi=phi, tag=tag) for k in range(K)]
    export = dict(meta=dict(phi=tuple(str(x) for x in phi), theta=tuple(str(x) for x in THETA), dps=dps, K=K, eps=[str(x) for x in eps], nmax=nmax,
                            branches=[(mu, c) for mu, c in BRANCHES]), branches={})
    for bi, br in enumerate(BRANCHES):
        ms = [monomials(o['branches'][bi], o['zeta'], THETA) for o in outs]
        ent = dict(N={}, M={}, P3jets={}, plucker={}, data={}, p={})
        def pack(vals, nz=None):
            c, rel = fit_stability(eps, vals, nmax)
            thr = mp.mpf(10) ** (-int(dps * 0.35))
            c = [mp.mpf(0) if abs(x) < thr else x for x in c[:nmax + 1]]
            n0 = next((i for i, x in enumerate(c) if x != 0), None)
            nrel = max([n for n, r in enumerate(rel) if r < mp.mpf(10) ** -10] + [-1])
            return dict(n0=n0, coef=[mp.nstr(x, dps) for x in c], relstab=[float(r) for r in rel], n_reliable=nrel)
        for s in range(2):
            for i, pr in enumerate(PAIRS):
                ent['N'][(s, pr)] = pack([m[s]['N'][i] for m in ms])
            ent['M'][s] = pack([m[s]['M'] for m in ms])
            for k in range(5):
                ent['P3jets'][(s, k)] = pack([m[s]['J'][2][k] for m in ms])
            for kl in [(k, l) for k in range(4) for l in range(k + 1, 4)]:
                ent['plucker'][(s, kl)] = pack([m[s]['p'][kl] for m in ms])
        for di, nm in enumerate(('H1', 'a1', 'b1', 'H2', 'a2', 'b2')):
            ent['data'][nm] = pack([o['branches'][bi]['data'][di] for o in outs])
        for j in range(3):
            for m in range(BRANCHES[bi][0][j] + 1):
                vals = [o['branches'][bi]['p'][j][m] for o in outs]
                # Laurent: find pole order e from the first two points, fit eps^e c_m
                v0 = abs(vals[-1]); v1 = abs(vals[-2])
                ex = int(round(float(mp.log(v1 / v0) / mp.log(2)))) if v0 != 0 and v1 != 0 else 0
                pole = max(0, -ex)
                sc = [v * mpf(e) ** pole for v, e in zip(vals, eps)]
                d = pack(sc)
                d['pole'] = pole
                ent['p'][(j, m)] = d
        export['branches'][br] = ent
        row = []
        for s in range(2):
            row.append(' '.join('%s:n0=%s,nrel=%d' % (pr, ent['N'][(s, pr)]['n0'], ent['N'][(s, pr)]['n_reliable']) for pr in PAIRS))
        print('  branch %s: first nonzero order n0 and reliable depth (rel.stab < 1e-10) of N_i(site1): %s' % (br, row[0]))
        print('      N_i(site2): %s' % row[1])
    return export

def stage_export():
    ex = stage_T2(200, 12, PHI, '', 10)
    with open(os.path.join(HERE, 'untwisted_Q_series_210.pkl'), 'wb') as f:
        pickle.dump(ex, f)
    ex2 = stage_T2(400, 24, PHI, '', 20)
    with open(os.path.join(HERE, 'untwisted_Q_series_210_ext.pkl'), 'wb') as f:
        pickle.dump(ex2, f)
    ex3 = stage_T2(200, 12, PHI2, '_phi2', 10)
    with open(os.path.join(HERE, 'untwisted_Q_series_210_phi2.pkl'), 'wb') as f:
        pickle.dump(ex3, f)

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--stages', default='T0,T1')
    args = ap.parse_args()
    st = args.stages.split(',')
    if 'T0' in st: stage_T0()
    if 'T1' in st: stage_T1()
    if 'EXPORT' in st: stage_export()
