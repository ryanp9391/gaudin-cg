"""M2: jet basis of one-row representations in the polynomial realisation (round 5, WP-A2).  Python only.
J0: dictionary GT basis <-> monomials x^alpha;  H0: e1^T rho(e^{aE12} e^{aE23}) = evaluation at (1,a,a^2).
"""
import sys, os, pickle, time, itertools
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))
import sympy as sp
from gl3_reps import GL3Rep
from gaudin_chain import SM, sm_from_sympy, comm, exp_nil


def monos(L):
    return [(a, b, L - a - b) for a in range(L, -1, -1) for b in range(L - a, -1, -1)]


def poly_E(L, i, j):
    """E_ij = x_i d_j on degree-L monomials, column convention M[target, source]; returns sympy Matrix and mono list"""
    ms = monos(L); ix = {m: k for k, m in enumerate(ms)}
    M = sp.zeros(len(ms), len(ms))
    for s, m in enumerate(ms):
        if m[j - 1] == 0:
            continue
        t = list(m); c = t[j - 1]; t[j - 1] -= 1; t[i - 1] += 1
        M[ix[tuple(t)], s] += c
    return M, ms


def dictionary(L):
    """diag scalars s_p with  v_p^{GT} = s_p x^{alpha(p)}, E^{GT} = D E^{poly} D^{-1}; returns (alpha list in GT order, s list, hw index)"""
    rep = GL3Rep((L, 0, 0)); d = rep.dim
    Eg = {(i, j): rep.generator(i, j) for (i, j) in [(1, 1), (2, 2), (3, 3), (1, 2), (2, 3), (2, 1), (3, 2)]}
    Eg[(1, 3)] = Eg[(1, 2)] * Eg[(2, 3)] - Eg[(2, 3)] * Eg[(1, 2)]
    Eg[(3, 1)] = Eg[(3, 2)] * Eg[(2, 1)] - Eg[(2, 1)] * Eg[(3, 2)]
    wt = [tuple(int(Eg[(k, k)][p, p]) for k in (1, 2, 3)) for p in range(d)]
    assert len(set(wt)) == d
    alpha = wt  # weights = exponents
    hw = alpha.index((L, 0, 0))
    s = [None] * d; s[hw] = sp.Integer(1)
    # propagate with E21, E32, E31 ... use lowering E21,E32 (and raising for robustness)
    ms = monos(L); ix = {m: k for k, m in enumerate(ms)}
    Ep = {k: poly_E(L, *k)[0] for k in Eg}
    changed = True
    while changed:
        changed = False
        for (i, j) in [(2, 1), (3, 2), (1, 2), (2, 3)]:
            for p in range(d):
                if s[p] is None:
                    continue
                for q in range(d):
                    if Eg[(i, j)][q, p] != 0 and s[q] is None:
                        # E^{GT}v_p = sum_q Eg[q,p] v_q ; E^{poly} s_p x^{alpha_p} = s_p * Ep[ix q, ix p] x^{alpha_q}
                        e_poly = Ep[(i, j)][ix[alpha[q]], ix[alpha[p]]]
                        s[q] = s[p] * e_poly / Eg[(i, j)][q, p]; changed = True
    assert all(v is not None for v in s)
    # verify all generators
    for k in Eg:
        for p in range(d):
            for q in range(d):
                lhs = Eg[k][q, p]
                rhs = s[p] / s[q] * Ep[k][ix[alpha[q]], ix[alpha[p]]]
                assert sp.simplify(lhs - rhs) == 0, (k, p, q, lhs, rhs)
    return alpha, s, hw, Eg


if __name__ == '__main__':
    for L in (1, 2, 3):
        alpha, s, hw, Eg = dictionary(L)
        print(L, 'hw', hw)
        for p, (a, sc) in enumerate(zip(alpha, s)):
            print('   GT', p, 'alpha', a, 's', sc)


# ================================================================ companion frame
import general_kernel as gk


def pinv_rows(R, zeta):
    """matrix M (list of rows, Fractions) with cov*M = cov * P_phys^{-1}; same recipe as gk.x_cov"""
    z1, z2, z3 = [Fr(x) for x in zeta]
    E = R.E
    s = [z1 ** 2 / ((z1 - z2) * (z1 - z3)), -(z2 / (z1 * (z2 - z3))), 1 / (z1 * z2)]
    r = [-(z3 / z2), -(z3 / z1), (z2 ** 2 * (z3 - z1)) / (z1 ** 2 * (z2 - z3))]
    rows = []
    for q in range(R.d):
        x = [Fr(0)] * R.d; x[q] = Fr(1)
        for (c, key) in ((-(z1 + z2) / (z1 * z2), (2, 3)), (1 / (z1 * z2), (1, 3)), (-(1 / z1), (1, 2))):
            x = exp_nil(E[key], c, 40).rapply(x)
        for p in range(R.d):
            f = Fr(1)
            for k in range(3):
                f *= s[k] ** (-R.wt[p][k])
            x[p] *= f
        for (c, key) in ((-r[2], (2, 1)), (-r[1], (3, 1)), (-r[0], (3, 2))):
            x = exp_nil(E[key], c, 40).rapply(x)
        rows.append(x)
    return rows


def vecmat(cov, M):
    n = len(M[0]); out = [Fr(0)] * n
    for i, c in enumerate(cov):
        if c != 0:
            for j in range(n):
                out[j] += c * M[i][j]
    return out


import mpmath as mp


def mp_to_fr(x, maxden=10 ** 28, tol=None):
    tol = tol or mp.mpf(10) ** -55
    if x == 0:
        return Fr(0)
    man, ex = mp.mpf(x).man, mp.mpf(x).exp
    f = Fr(int(man)) * (Fr(2) ** int(ex))
    if x < 0:
        f = -f
    r = f.limit_denominator(maxden)
    ref = max(abs(x), mp.mpf(1))
    if abs(x - mp.mpf(r.numerator) / r.denominator) < tol * ref:
        return r
    return None


def load_pt(lam, tag, th):
    return gk.get_point2(lam, tag, gk.ZETAS2[tag] if tag in gk.ZETAS2 else None, th)


def th_from_key(key):
    a, b = key.split('_')
    f = lambda s: Fr(int(s.split('.')[0]), int(s.split('.')[1]))
    return (f(a), f(b))


def comp_states(lam, pt, M):
    """companion kets Psi_c as mp matrix columns (d^2 x nstates) and x0^c.Psi_c"""
    d = pt['d1']
    Mm = [[mp.mpf(x.numerator) / x.denominator for x in row] for row in M]
    cols = []
    for st in pt['states']:
        vec = {i: v for i, v in zip(st['idx'], st['v'])}
        # (M x M) psi : psi_c[(p,q)] = sum_{a,b} M[p][a] M[q][b] psi[(a,b)]
        out = [mp.mpf(0)] * (d * d)
        for p in range(d):
            for q in range(d):
                t = mp.mpf(0)
                for (i, v) in vec.items():
                    a, b = divmod(i, d)
                    t += Mm[p][a] * Mm[q][b] * v
                out[p * d + q] = t
        cols.append(out)
    return cols


def jet_Y(lam, pt):
    """Y (npairs^2 x d^2) in companion frame, mp: rows (i,j) over pairs_of(lam) x pairs_of(lam)"""
    mp.mp.dps = pt['dps']
    R = gk.one(lam); d = pt['d1']
    M = pinv_rows(R, pt['zeta'])
    cols = comp_states(lam, pt, M)
    sm = gk.site_monos(lam, pt)
    prs = gk.pairs_of(lam); n = len(prs)
    assert len(cols) == d * d == n * n
    e = R.hw_index
    W = mp.zeros(n * n, d * d); Psi = mp.zeros(d * d, d * d)
    for al, (col, per) in enumerate(zip(cols, sm)):
        x0 = col[e * d + e]
        (m1, M1), (m2, M2) = per
        for i in range(n):
            for j in range(n):
                W[i * n + j, al] = m1[i] * m2[j] / (M1 * M2) * x0
        for r in range(d * d):
            Psi[r, al] = col[r]
    Y = W * mp.inverse(Psi)
    return Y, Psi


def recon(Y, maxden=10 ** 12):
    out = []
    bad = 0
    for i in range(Y.rows):
        row = []
        for j in range(Y.cols):
            f = mp_to_fr(Y[i, j], maxden=maxden)
            if f is None:
                bad += 1
            row.append(f)
        out.append(row)
    return out, bad


# ================================================================ t-structure
def Y_poly(lam, tag='P', keys=None, deg=None):
    """Y(t) = sum_n t^n Y_n (sympy Matrices, exact), from all cached theta points with this tag (dedupe in t)."""
    import glob, re
    L = lam[0]
    pat = 'general_kernel_pt_%d%d%d_%s_*_80.pkl' % (lam[0], lam[1], lam[2], tag)
    pts = {}
    for f in sorted(glob.glob(os.path.join(HERE, pat))):
        key = re.match(r'.*_%s_(.+)_80.pkl' % tag, os.path.basename(f)).group(1)
        th = th_from_key(key)
        pt = gk.get_point2(lam, tag, gk.ZETAS2[tag], th)
        t = 1 / (th[0] - th[1])
        if t in pts:
            continue
        Y, _ = jet_Y(lam, pt)
        Yr, bad = recon(Y)
        assert bad == 0, (key, bad)
        pts[t] = Yr
    ts = list(pts)
    n = len(pts[ts[0]]); m = len(pts[ts[0]][0])
    D = deg if deg else len(ts) - 3
    use = ts[:D + 1]; chk = ts[D + 1:]
    assert len(chk) >= 1
    tt = sp.symbols('t')
    coef = [sp.zeros(n, m) for _ in range(D + 1)]
    # Vandermonde solve (exact)
    V = sp.Matrix([[sp.Rational(t.numerator, t.denominator) ** k for k in range(D + 1)] for t in use])
    Vi = V.inv()
    for i in range(n):
        for j in range(m):
            ys = sp.Matrix([sp.Rational(pts[t][i][j].numerator, pts[t][i][j].denominator) for t in use])
            c = Vi * ys
            for k in range(D + 1):
                coef[k][i, j] = c[k]
    # check
    for t in chk:
        tv = sp.Rational(t.numerator, t.denominator)
        Ym = sum((coef[k] * tv ** k for k in range(D + 1)), sp.zeros(n, m))
        Yt = sp.Matrix(n, m, lambda i, j: sp.Rational(pts[t][i][j].numerator, pts[t][i][j].denominator))
        assert Ym == Yt, 'poly fit failed'
    while len(coef) > 1 and coef[-1] == sp.zeros(n, m):
        coef.pop()
    return coef, len(ts)


def omega12(lam):
    R = gk.one(lam); d = R.d
    from gaudin_chain import kron
    tot = SM(d * d)
    I = SM.eye(d)
    for i in (1, 2, 3):
        for k in (1, 2, 3):
            tot = tot + kron(R.E[(i, k)], R.E[(k, i)])
    return tot


a_s, ap_s, t_s = sp.symbols('a ap t')


def ev_cov(L):
    """symbolic evaluation covector (two-site, GT components of companion frame): Omega_c[a] (x) Omega_c[a']"""
    alpha, s, hw, Eg = dictionary(L); d = len(alpha)
    v1 = [sp.Rational(s[p]) * a_s ** (alpha[p][1] + 2 * alpha[p][2]) for p in range(d)]
    v2 = [sp.Rational(s[p]) * ap_s ** (alpha[p][1] + 2 * alpha[p][2]) for p in range(d)]
    return sp.Matrix([[v1[p] * v2[q] for p in range(d) for q in range(d)]])


def Theta_from_Y(L, coef):
    Yt = sum((coef[k] * t_s ** k for k in range(len(coef))), sp.zeros(*coef[0].shape))
    Yi = Yt.inv()  # rational in t
    Yi = Yi.applyfunc(sp.cancel)
    th = (ev_cov(L) * Yi).applyfunc(sp.factor)
    return th


# ================================================================ differential-operator fit of A_n (universal in L)
def mono_index(L):
    ms = monos(L); return ms, {m: k for k, m in enumerate(ms)}


def diffop_matrix(L, xs, ds):
    """single-site operator  x_{xs[0]}...x_{xs[p-1]} d_{ds[0]}...d_{ds[p-1]} on degree-L polynomials; returns dict (target,source)->Fraction (monomial basis)"""
    ms, ix = mono_index(L)
    out = {}
    for s, m in enumerate(ms):
        e = list(m); c = 1
        for k in ds:
            if e[k - 1] == 0:
                c = 0; break
            c *= e[k - 1]; e[k - 1] -= 1
        if c == 0:
            continue
        for k in xs:
            e[k - 1] += 1
        out[(ix[tuple(e)], s)] = out.get((ix[tuple(e)], s), 0) + Fr(c)
    return out


def op_basis(maxorder):
    """list of (name, xs, ds, gl3 weight shift) for orders 0..maxorder; xs,ds sorted tuples"""
    out = []
    for p in range(maxorder + 1):
        for xs in itertools.combinations_with_replacement((1, 2, 3), p):
            for ds in itertools.combinations_with_replacement((1, 2, 3), p):
                w = [0, 0, 0]
                for k in xs: w[k - 1] += 1
                for k in ds: w[k - 1] -= 1
                out.append(((xs, ds), tuple(w)))
    return out


def fit_A(An_by_L, shift, maxorder1, maxorder2, Ls):
    """An_by_L[L] = sympy matrix of A_n in the GT basis of the two-site space. Returns solution dict {((xs1,ds1),(xs2,ds2)): coeff} for weight shift `shift`
    (stacked over all L in Ls: universal L-independent coefficients), free variables = 0; plus rank/consistency."""
    base = op_basis(max(maxorder1, maxorder2))
    b1 = [o for o in base if len(o[0][0]) <= maxorder1]
    b2 = [o for o in base if len(o[0][0]) <= maxorder2]
    pairs = [(o1, o2) for o1 in b1 for o2 in b2 if tuple(a + b for a, b in zip(o1[1], o2[1])) == tuple(shift)]
    rows = []; rhs = []
    for L in Ls:
        alpha, s, hw, Eg = dictionary(L); d = len(alpha)
        ms, ix = mono_index(L)
        gt2mono = [ix[alpha[p]] for p in range(d)]
        sc = [Fr(int(sp.Rational(x).p), int(sp.Rational(x).q)) for x in s]
        An = An_by_L[L]
        # monomial-basis matrix of A_n :  c' = S A S^-1 c ; index two-site (p,q)->(gt2mono[p], gt2mono[q])
        # target data: entries with weight shift == shift
        w3 = [tuple(alpha[p][k] + alpha[q][k] for k in range(3)) for p in range(d) for q in range(d)]
        mats = []
        for (o1, o2) in pairs:
            m1 = diffop_matrix(L, *o1[0]); m2 = diffop_matrix(L, *o2[0])
            mats.append((m1, m2))
        # equations: for each (r,c) two-site GT index pair with right weight shift (others are automatically zero for both sides)
        # build in monomial indices
        mono_of = {}
        for p in range(d):
            mono_of[gt2mono[p]] = p
        for r in range(d * d):
            pr, qr = divmod(r, d)
            for c in range(d * d):
                if tuple(w3[r][k] - w3[c][k] for k in range(3)) != tuple(shift):
                    continue
                pc, qc = divmod(c, d)
                tr1, tr2, sr1, sr2 = gt2mono[pr], gt2mono[qr], gt2mono[pc], gt2mono[qc]
                row = [m1.get((tr1, sr1), 0) * m2.get((tr2, sr2), 0) for (m1, m2) in mats]
                val = An[r, c] * sc[pr] * sc[qr] / (sc[pc] * sc[qc])
                val = Fr(int(sp.Rational(val).p), int(sp.Rational(val).q))
                if any(v != 0 for v in row) or val != 0:
                    rows.append(row); rhs.append(val)
    x, rank, cons = gk.exact_solve(rows, rhs)
    sol = {pairs[i]: x[i] for i in range(len(pairs)) if x is not None and x[i] != 0}
    return sol, rank, cons, len(pairs), len(rows)


# ================================================================ Frobenius derivation of nu_kl from the Laurent data (one-row reps: exponents {0,1,L+2}, c2 = 0)
Hs, As, Bs, Ht, At, Bt, E1, E2, E3, DL, xx = sp.symbols('H_s a_s b_s H_t a_t b_t e1 e2 e3 delta x')


def laurent_series(L, N):
    """coefficient lists (Laurent in x, dict power->coeff) of S1,S2,S3 at the pole theta_s up to x^N (delta = theta_s - theta_t)"""
    def inv1(n):  # coefficient of x^n in 1/(x+delta)
        return (-1) ** n / DL ** (n + 1)
    def inv2(n):  # coefficient of x^n in 1/(x+delta)^2
        return (n + 1) * (-1) ** n / DL ** (n + 2)
    S1 = {-1: sp.Integer(L)}; S2 = {-1: Hs}; S3 = {-2: As, -1: Bs}
    for n in range(0, N + 1):
        S1[n] = L * inv1(n)
        S2[n] = Ht * inv1(n)
        S3[n] = Bt * inv1(n) + At * inv2(n)
    S1[0] += E1; S2[0] += E2; S3[0] += E3
    return S1, S2, S3


def frobenius_planes(L):
    """returns (gamma, d) lists (length L+2) of Frobenius coefficients of f0 = 1 + 0 x + ..., f1 = x + ... and the no-log condition(s)"""
    N = L + 3
    S1, S2, S3 = laurent_series(L, N)
    conds = []
    def solve(rho, c0, kmax, fixed):
        c = {0: sp.Integer(c0)}; c.update(fixed)
        unknown = {}
        for k in range(1, kmax + 1):
            if k in c:
                continue
            c[k] = sp.Symbol('c%d' % k)
            unknown[k] = c[k]
        # coefficient of x^{rho+m-3}, m = 0..kmax
        res = {}
        for m in range(0, kmax + 1):
            eq = 0
            for k in range(0, m + 1):
                n = rho + k
                # psi''' = n(n-1)(n-2) x^{n-3};  -S1 psi'' ; S2 psi' ; -S3 psi ; collect power x^{rho+m-3}
                if k == m:
                    eq += c[k] * n * (n - 1) * (n - 2)
                # S1 psi'' : x^{n-2} * S1[p] x^p -> power n-2+p = rho+m-3 -> p = m-k-1
                p = m - k - 1
                if p in S1:
                    eq += -S1[p] * c[k] * n * (n - 1)
                # S2 psi': x^{n-1} S2[p]: p = m-k-2
                p = m - k - 2
                if p in S2:
                    eq += S2[p] * c[k] * n
                # S3 psi: x^n S3[p]: p = m-k-3
                p = m - k - 3
                if p in S3:
                    eq += -S3[p] * c[k]
            res[m] = sp.expand(eq)
        return c, res, unknown
    # f0: rho = 0, gamma_1 = 0 gauge
    c0_, res0, unk0 = solve(0, 1, L + 1, {1: sp.Integer(0)})
    sol = {}
    cond0 = None
    for m in range(1, L + 2):
        eq = res0[m].subs(sol)
        if m == 1:
            cond0 = sp.simplify(eq)  # no-log condition (does not involve unknowns since gamma_1 gauge-fixed = 0 and indicial vanishes)
            continue
        ck = sp.Symbol('c%d' % m)
        s = sp.solve(eq, ck)
        sol[ck] = sp.simplify(s[0])
    gam = [sp.Integer(1), sp.Integer(0)] + [sp.simplify(sol[sp.Symbol('c%d' % k)]) for k in range(2, L + 2)]
    # f1: rho = 1
    c1_, res1, unk1 = solve(1, 1, L, {})
    sol1 = {}
    for m in range(1, L + 1):
        eq = res1[m].subs(sol1)
        ck = sp.Symbol('c%d' % m)
        s = sp.solve(eq, ck)
        sol1[ck] = sp.simplify(s[0])
    # f1 = sum c_k x^{1+k}: coefficient of x^j is c_{j-1}
    d = [sp.Integer(0), sp.Integer(1)] + [sp.simplify(sol1[sp.Symbol('c%d' % (j - 1))]) for j in range(2, L + 2)]
    return gam, d, cond0


def nu_formulas(L):
    gam, d, cond = frobenius_planes(L)
    nu = {}
    for k in range(L + 2):
        for l in range(k + 1, L + 2):
            nu[(k, l)] = sp.simplify(sp.factorial(k) * sp.factorial(l) * (gam[k] * d[l] - gam[l] * d[k]))
    return nu, cond


from gaudin_chain import Chain, Operators, gc


def companion_ops(lam, zeta, th):
    """Laurent-coefficient operators (companion frame, twist g_c(zeta)) as dict of SM: H1,H2,a1,a2,b1,b2 + data; exact"""
    ch = Chain(lam=lam, theta=th)
    gcm = gc(zeta)
    gmat = [[Fr(int(gcm[i, j].p), int(gcm[i, j].q)) for j in range(3)] for i in range(3)]
    ops = Operators(ch, gmat)
    C = ops.fit()
    assert C['res2'] == 0 and C['res3'] == 0
    return dict(H1=C[2][1], H2=C[2][3], a1=C[3][2], a2=C[3][5], b1=C[3][1], b2=C[3][4], c2a=C[2][2], t1=C[3][3], ch=ch)


def eval_poly_on_cov(poly, gens, ops_by_sym, cov, dim):
    """cov * poly(ops) with commuting ops; poly = sympy Poly in symbols `gens` with Fraction/rational coefficients"""
    out = [Fr(0)] * dim
    for monom, coeff in poly.terms():
        v = list(cov)
        for g, pw in zip(gens, monom):
            for _ in range(pw):
                v = ops_by_sym[g].rapply(v)
        cf = sp.Rational(coeff)
        cf = Fr(int(cf.p), int(cf.q))
        out = [o + cf * x for o, x in zip(out, v)]
    return out


def Y_from_operators(lam, zeta, th, nu=None):
    """Y rows (i,j) over pairs_of(lam)^2 as covectors over the GT product basis: <x0| nu_i(site1 data) nu_j(site2 data) (exact), companion frame"""
    L = lam[0]
    if nu is None:
        nu, cond = nu_formulas(L)
    O = companion_ops(lam, zeta, th)
    ch = O['ch']; dim = ch.dim; d = ch.d1
    z1, z2, z3 = [Fr(x) for x in zeta]
    e1v = z1 + z2 + z3; e2v = z1 * z2 + z1 * z3 + z2 * z3; e3v = z1 * z2 * z3
    delta = th[0] - th[1]
    R = gk.one(lam)
    x0 = [Fr(0)] * dim; x0[R.hw_index * d + R.hw_index] = Fr(1)
    prs = gk.pairs_of(lam); n = len(prs)
    gens = [Hs, As, Bs, Ht, At, Bt]
    rows = []
    # site 1: (s,t) = (1,2), delta = th1 - th2 ; site 2: (s,t)=(2,1), delta = th2-th1
    sub = lambda dd: {E1: sp.Rational(e1v.numerator, e1v.denominator), E2: sp.Rational(e2v.numerator, e2v.denominator),
                      E3: sp.Rational(e3v.numerator, e3v.denominator), DL: sp.Rational(dd.numerator, dd.denominator)}
    nu1 = {kl: sp.Poly(sp.expand(sp.together(v.subs(sub(delta))) ), *gens) for kl, v in nu.items()}
    nu2 = {kl: sp.Poly(sp.expand(sp.together(v.subs(sub(-delta))) ), *gens) for kl, v in nu.items()}
    ops1 = {Hs: O['H1'], As: O['a1'], Bs: O['b1'], Ht: O['H2'], At: O['a2'], Bt: O['b2']}
    ops2 = {Hs: O['H2'], As: O['a2'], Bs: O['b2'], Ht: O['H1'], At: O['a1'], Bt: O['b1']}
    for i in range(n):
        v1 = eval_poly_on_cov(nu1[prs[i]], gens, ops1, x0, dim)
        for j in range(n):
            rows.append(eval_poly_on_cov(nu2[prs[j]], gens, ops2, v1, dim))
    return rows


def Y_state_rational(lam, tag, th):
    pt = gk.get_point2(lam, tag, gk.ZETAS2[tag], th)
    Y, _ = jet_Y(lam, pt)
    return recon(Y)[0]


def poly_fit_in_t(res, D):
    """res: dict t -> rows (list of lists of Fractions). Fit entries by polynomials of degree D in t (exact), verify on remaining points. returns list of coefficient matrices (sympy)"""
    ts = list(res)
    use, chk = ts[:D + 1], ts[D + 1:]
    assert chk, 'need more points than D+1'
    n = len(res[ts[0]]); m = len(res[ts[0]][0])
    V = sp.Matrix([[sp.Rational(t.numerator, t.denominator) ** k for k in range(D + 1)] for t in use]).inv()
    coef = [sp.zeros(n, m) for _ in range(D + 1)]
    for i in range(n):
        for j in range(m):
            ys = sp.Matrix([sp.Rational(res[t][i][j].numerator, res[t][i][j].denominator) for t in use])
            if all(y == 0 for y in ys):
                continue
            c = V * ys
            for k in range(D + 1):
                coef[k][i, j] = c[k]
    for t in chk:
        tv = sp.Rational(t.numerator, t.denominator)
        for i in range(n):
            for j in range(m):
                v = sum(coef[k][i, j] * tv ** k for k in range(D + 1))
                assert v == sp.Rational(res[t][i][j].numerator, res[t][i][j].denominator), 'degree too low'
    while len(coef) > 1 and coef[-1] == sp.zeros(n, m):
        coef.pop()
    return coef


def theta_check(L, coef, vals_list):
    T, rep = gk.load_rep((L, 0, 0))
    prs = gk.pairs_of((L, 0, 0)); n = len(prs)
    evs = ev_cov(L)
    out = []
    for vals in vals_list:
        tv = vals[t_s]
        Yt = sum((coef[k] * tv ** k for k in range(len(coef))), sp.zeros(*coef[0].shape))
        ev = evs.subs(vals)
        th = Yt.T.LUsolve(ev.T)
        bad = 0
        for i in range(n):
            for j in range(n):
                ref = sum(sp.Rational(c.numerator, c.denominator) * tv ** nn * vals[a_s] ** p * vals[ap_s] ** q for (nn, p, q), c in T[(i, j)].items())
                if th[i * n + j] != ref:
                    bad += 1
        out.append(bad)
    return out
