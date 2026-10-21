"""Q1 (round 5, WP-D): the Gaudin norm from the Q side.  python -X utf8 gaudin_norm.py [--stages D0,G,...]
Reuses twovac_fit / gaudin_B read-only.  Basis: psi vectors in gl3_reps order (twovac_fit), companion frame / jet matrix in notebook order (gaudin_B)."""
import sys, os, time, pickle
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mpmath as mp
import twovac_fit as T
import gaudin_B as B

mp.mp.dps = 80
NB2G = B.NB2G
TH = T.THETA
ZET = B.ZETAS


def mpf(x):
    return T.mpf(x)


def cpath(name):
    return os.path.join(HERE, 'gaudin_norm_%s.pkl' % name)


def cached(name, fn):
    p = cpath(name)
    if os.path.exists(p):
        with open(p, 'rb') as f:
            return pickle.load(f)
    r = fn()
    with open(p, 'wb') as f:
        pickle.dump(r, f)
    return r


def gdiag_nb():
    """one-site invariant diagonal metric in notebook order (Fractions)"""
    d = T.gram_diag(T.chain())
    return [d[NB2G[a]] for a in range(8)]


def Gc_exact(tag):
    """metric G (x) G transported to the companion frame: Psi_d = (S (x) S) Psi_c, S = P_phys^{-1}(zeta)^{-1}; G_c = (S(x)S)^T GG (S(x)S) (exact 64x64, notebook order)"""
    zeta = ZET[tag]
    Pi = B.pinv_exact(zeta)
    S = B.frac_inverse(Pi)
    d = gdiag_nb()
    # one-site G_c1 = S^T diag(d) S; two-site metric is the Kronecker square
    n = 8
    G1 = [[sum(S[k][i] * d[k] * S[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    return G1, B.kron_mat(G1, G1)


def comp_Y(tag='P'):
    r = B.comp_jet_matrix('P')
    return r['Y'], r['Yinv']


def load(tag, th=TH):
    """all per-state data at a zeta point: kets (gl3 order), x0 overlaps, ground-truth norms, one-site monomials, Laurent data"""
    T.set_precision(80)
    mp.mp.dps = 80
    pt = T.get_point(tag, th=th)
    ch = T.chain(th)
    sj = [T.StateJets(st, pt['zeta'], th, n=6) for st in pt['states']]
    d = T.gram_diag(ch)
    dd = [mpf(d[a // 8] * d[a % 8]) for a in range(64)]
    x0 = T.x0_D2(ch, pt['zeta'])
    x0b = [a * b for a in T.x0_onesite(ch, pt['zeta']) for b in T.x0_onesite(ch, pt['zeta'])]
    out = []
    for st, J in zip(pt['states'], sj):
        v = st['vec']
        nrm = mp.fsum(dd[a] * v[a] ** 2 for a in range(64))
        ox = T.overlap(x0, v)
        oxb = T.overlap(x0b, v)
        mons = [B.mono_funcs([J], s)[0] for s in range(2)]
        out.append(dict(vec=v, norm=nrm, x0=ox, x0b=oxb, N=nrm / ox ** 2, Nb=nrm / oxb ** 2, mons=mons, lam2=st['lam2'], lam3=st['lam3'], sector=st['sector'], ev=st['ev']))
    return pt, sj, out



def G1_sym():
    """one-site companion-frame metric G_c1(zeta) = S^T D S, symbolic in (z1,z2,z3) (sympy 8x8; cached as srepr strings)"""
    import sympy as sp
    z1, z2, z3 = sp.symbols('z1 z2 z3')
    def build():
        def mat(key):
            Ek = B.CH.E1[key]
            return sp.Matrix(8, 8, lambda i, j: sp.Rational(Ek.r.get(i, {}).get(j, Fr(0)).numerator, Ek.r.get(i, {}).get(j, Fr(0)).denominator))
        E = {k: mat(k) for k in B.CH.E1}
        def expN(c, key):
            M = E[key]
            return sp.eye(8) + c * M + c ** 2 * M * M / 2
        s = [z1 ** 2 / ((z1 - z2) * (z1 - z3)), -(z2 / (z1 * (z2 - z3))), 1 / (z1 * z2)]
        r = [-(z3 / z2), -(z3 / z1), (z2 ** 2 * (z3 - z1)) / (z1 ** 2 * (z2 - z3))]
        pp = expN(-(z1 + z2) / (z1 * z2), (2, 3)) * expN(1 / (z1 * z2), (1, 3)) * expN(-(1 / z1), (1, 2))
        wts = [[int(B.CH.E1[(k, k)].r.get(p, {}).get(p, 0)) for k in (1, 2, 3)] for p in range(8)]
        dinv = sp.diag(*[s[0] ** (-w[0]) * s[1] ** (-w[1]) * s[2] ** (-w[2]) for w in wts])
        pm = expN(-r[2], (2, 1)) * expN(-r[1], (3, 1)) * expN(-r[0], (3, 2))
        S = sp.simplify((pp * dinv * pm).inv())
        d = gdiag_nb()
        D = sp.diag(*[sp.Rational(x.numerator, x.denominator) for x in d])
        G1 = (S.T * D * S).applyfunc(sp.factor)
        return [[sp.srepr(G1[i, j]) for j in range(8)] for i in range(8)]
    rows = cached('G1sym', build)
    return sp.Matrix(8, 8, lambda i, j: sp.sympify(eval(rows[i][j], sp.__dict__)))


_JAC = {}
def jac_funcs():
    """symbolic no-log conditions F(H1,a1,b1,H2,a2,b2; e1,e2,e3,d=th1-th2) (6 of them) and their 6x6 Jacobian, lambdified to mpmath"""
    if 'f' in _JAC:
        return _JAC['f']
    import sympy as sp
    import untwisted_Q as U
    vars_, C = U.frobenius_conditions()
    D = sp.symbols('H1 a1 b1 H2 a2 b2')
    e1, e2, e3, d = sp.symbols('e1 e2 e3 d')
    Fs = []
    for s in range(2):
        o = 1 - s
        ds = d if s == 0 else -d
        def reg(k, m):
            return (-1) ** k * sp.binomial(m + k - 1, k) * ds ** (-m - k)
        H, a, b = D[3 * s:3 * s + 3]
        Ho, ao, bo = D[3 * o:3 * o + 3]
        sk = [3 * reg(k, 1) + (e1 if k == 0 else 0) for k in range(3)]
        hk = [Ho * reg(k, 1) + 3 * reg(k, 2) + (e2 if k == 0 else 0) for k in range(3)]
        tk = [bo * reg(k, 1) + ao * reg(k, 2) + (e3 if k == 0 else 0) for k in range(3)]
        sub = dict(zip(vars_, [H, a, b] + sk + hk + tk))
        for (n, idx, poly) in C:
            Fs.append(sp.expand(poly.subs(sub)))
    Fm = sp.Matrix(Fs)
    J = Fm.jacobian(sp.Matrix(D))
    args = D + (e1, e2, e3, d)
    out = (Fm, J, sp.lambdify(args, Fm, 'mpmath'), sp.lambdify(args, J, 'mpmath'), D, (e1, e2, e3, d))
    _JAC['f'] = out
    return out


def laurent(st):
    return [st['lam2'][1], st['lam3'][2], st['lam3'][1], st['lam2'][3], st['lam3'][5], st['lam3'][4]]     # (H1,a1,b1,H2,a2,b2)


def jac_at(st, zeta, th=TH):
    Fm, J, Ff, Jf, D, _ = jac_funcs()
    z = [mpf(x) for x in zeta]
    e1 = sum(z); e2 = z[0] * z[1] + z[0] * z[2] + z[1] * z[2]; e3 = z[0] * z[1] * z[2]
    d = mpf(th[0] - th[1])
    args = laurent(st) + [e1, e2, e3, d]
    return Ff(*args), mp.matrix(Jf(*args))


def G1_exact_zeta(zeta):
    """exact one-site companion-frame metric at arbitrary rational zeta (8x8 Fractions)"""
    Pi = B.pinv_exact(tuple(Fr(x) for x in zeta))
    S = B.frac_inverse(Pi)
    d = gdiag_nb()
    return [[sum(S[k][i] * d[k] * S[k][j] for k in range(8)) for j in range(8)] for i in range(8)]


_YI = {}
def Yinv_mp():
    if 'y' not in _YI:
        Y, Yinv = comp_Y()
        _YI['y'] = (Y, Yinv)
    return _YI['y']


def Gamma_mp(zeta, t=None):
    """Gamma = Yinv(t)^T (G1 (x) G1) Yinv(t) as an mpmath 64x64 matrix (current mp.dps); Y(t) = t^L Y(1) t^-L, t = 1/theta12 (default value of the base theta);
    Yinv stored at theta = (2/7,-3/5) i.e. t = 35/31."""
    t0 = Fr(35, 31)
    Y, Yinv = Yinv_mp()
    t = t0 if t is None else Fr(t)
    LV = [0, 1, 1, 2, 2, 3, 3, 4]
    lev = [LV[i // 8] + LV[i % 8] for i in range(64)]
    G1 = G1_exact_zeta(zeta)
    G1m = mp.matrix(8, 8)
    for i in range(8):
        for j in range(8):
            G1m[i, j] = mpf(G1[i][j])
    G2 = mp.matrix(64, 64)
    for a in range(64):
        for b in range(64):
            G2[a, b] = G1m[a // 8, b // 8] * G1m[a % 8, b % 8]
    # Yinv(t)_{c a} = Yinv(t0)_{c a} (t/t0)^{lev_c - lev_a}
    r = t / t0
    Ym = mp.matrix(64, 64)
    for c in range(64):
        for a in range(64):
            if Yinv[c][a] != 0:
                Ym[c, a] = mpf(Yinv[c][a] * r ** (lev[c] - lev[a]))
    return Ym.T * G2 * Ym


def nu8(nu5):
    """the eight one-site nu's (m_i/M, notebook W_NAMES order) from (1, kappa, pi, kappa pi, eps): m4 = P3' p12 = kappa^2, m6 = kappa^2 pi, m7 = kappa eps"""
    _, k, p, kp, e = nu5
    return [mp.mpf(1), k, p, kp, k * k, e, k * k * p, k * e]


def nu2(n1, n2):
    return [a * b for a in n1 for b in n2]


KB = [0, 1, 1, 2, 2, 4]


def N_formula(b, zeta, th=TH):
    """N(Psi) = <Psi|Psi>_G/<x0|Psi>^2 = (9/4) |det J6| prod_c |zeta_a - zeta_b|^(mu_c - 6)  (pair (a,b) complementary to c; mu = sector); J6 = Jacobian of the six no-log conditions"""
    F, J = jac_at(b, zeta, th)
    z = [mpf(q) for q in zeta]
    d12, d13, d23 = abs(z[0] - z[1]), abs(z[0] - z[2]), abs(z[1] - z[2])
    m = b['sector']
    return mp.mpf(9) / 4 * abs(mp.det(J)) * d23 ** (m[0] - 6) * d13 ** (m[1] - 6) * d12 ** (m[2] - 6)


def cg_table(phi, tag, dps=200, verbose=True):
    """D2: Q-side fully normalised CG slices of adj (x) adj from Q data only:  Psi_0^Q = X_Q * v_Q,  X_Q = 1/sqrt(lim eps^(4+2kB) N_formula),
    v_Q = the HW ket (in H_mu) whose GF polynomial equals the leading nu-form polynomial L_B (U2's dictionary, shape only)."""
    import untwist_limit as UL
    mp.mp.dps = dps
    pts = UL.load_all(phi, tag)
    ch = T.chain()
    Hs = {mu: [[T.mpf(x) for x in h] for h in T.hw_space(ch, mu)] for mu in T.HW_SECT}
    Fpolys, chat = UL.gf_polys(ch, phi)
    Fd = [UL.poly_dict(P) if P is not None else {} for P in Fpolys]
    G = UL.G_pieces()
    out = []
    for bi in range(6):
        mu = UL.SECT[bi]
        # Q-side norm
        vals = []
        for p in pts:
            eps = mp.mpf(p['eps'].numerator) / p['eps'].denominator
            vals.append(N_formula(p['branches'][bi], p['zeta']) * eps ** (4 + 2 * KB[bi]))
        c0 = UL.series_fit(vals, 12)[0]
        XQ = 1 / mp.sqrt(c0)
        # leading polynomial from the nu-series
        nus = UL.nu_series(pts, bi)
        NBp = 6 - KB[bi]
        R, trip = UL.leading_poly(nus, NBp + 1, G)
        LB = R[NBp]
        mons = sorted(set(LB) | set(m_ for nu in range(8) for m_ in Fd[nu]))
        H = Hs[mu]
        Am = mp.matrix(len(mons), len(H))
        for i_, h in enumerate(H):
            for nu in range(8):
                for m_, c_ in Fd[nu].items():
                    Am[mons.index(m_), i_] += h[nu * 8 + 7] * UL.mpq(c_)
        l = mp.matrix([LB.get(m_, mp.mpf(0)) for m_ in mons])
        beta = mp.lu_solve(Am.T * Am, Am.T * l)
        resid = mp.norm(Am * beta - l) / mp.norm(l)
        vq = [mp.fsum(beta[i_] * H[i_][j] for i_ in range(len(H))) for j in range(64)]
        d = T.gram_diag(ch)
        nvq = mp.sqrt(mp.fsum(mpf(d[a // 8] * d[a % 8]) * vq[a] ** 2 for a in range(64)))
        # chain Psi_0 (G-norm 1, U2's sign convention)
        vs = UL.branch_psi(pts, bi, ch)
        PS = UL.comp_series(vs)
        psi0 = PS[0]
        v_Q = [XQ * x for x in vq]
        # sign: align with chain on the largest component
        j0 = max(range(64), key=lambda j: abs(psi0[j]))
        sg = 1 if v_Q[j0] * psi0[j0] > 0 else -1
        v_Q = [sg * x for x in v_Q]
        dev = max(abs(x - y) for x, y in zip(v_Q, psi0))
        sl = [(nu, v_Q[nu * 8 + 7], psi0[nu * 8 + 7]) for nu in range(8) if abs(psi0[nu * 8 + 7]) > mp.mpf(10) ** -30 or abs(v_Q[nu * 8 + 7]) > mp.mpf(10) ** -30]
        rec = dict(bi=bi, mu=mu, XQ=XQ, XQ_norm_vq=XQ * nvq, resid=resid, dev=dev, slice=sl, depth_series=NBp)
        out.append(rec)
        if verbose:
            print('%s: X_Q = 1/sqrt(lim eps^%d N) = %s ; ||X_Q v_q||_G = %s ; dictionary resid %s ; max|Psi0^Q - Psi0^chain| over 64 comps = %s' % (
                UL.BNAME[bi], 4 + 2 * KB[bi], mp.nstr(XQ, 20), mp.nstr(XQ * nvq, 20), mp.nstr(resid, 3), mp.nstr(dev, 3)))
            print('   slice C_nu (Q, chain): %s' % [(nu, mp.nstr(a, 15), mp.nstr(b, 15)) for nu, a, b in sl])
    return out


# =============================================================================================== stages
hdr = T.hdr


def stage_D0():
    hdr('D0: ground truth N = <Psi|Psi>_G/<x0|Psi>^2 for the 64 states at P, Pp, X1 (theta=(2/7,-3/5)); data table cached gaudin_norm_D0_<tag>.pkl')
    for tag in ('P', 'Pp', 'X1'):
        pt, sj, S = load(tag)
        dd = T.detdata(pt)
        recs = []
        for k, (s, st, J) in enumerate(zip(S, pt['states'], sj)):
            F, Jm = jac_at(st, pt['zeta'])
            recs.append(dict(sector=st['sector'], laurent=[str(x) for x in laurent(st)], m1=[str(x) for x in s['mons'][0]], m2=[str(x) for x in s['mons'][1]],
                             M1=str(s['mons'][0][0]), M2=str(s['mons'][1][0]), detN=str(dd['dets'][k]), Rv=str(dd['Rv'][k]), x0=str(s['x0']), norm=str(s['norm']),
                             N=str(s['N']), Nv=str(s['norm'] / dd['Rv'][k] ** 2), J6=str(mp.det(Jm))))
        with open(cpath('D0_' + tag), 'wb') as f:
            pickle.dump(recs, f)
        Ns = [s['N'] for s in S]
        Nv = [s['norm'] / dd['Rv'][k] ** 2 for k, s in enumerate(S)]
        rr = [(dd['Rv'][k] / s['x0']) * s['mons'][0][0] * s['mons'][1][0] / dd['dets'][k] for k, s in enumerate(S)]
        print('%s: zeta=%s  N range [%s, %s], N_Omega range [%s, %s]; (Rv/x0) M1 M2/detN = %s, spread %s' % (tag, [str(x) for x in pt['zeta']], mp.nstr(min(Ns), 5), mp.nstr(max(Ns), 5),
              mp.nstr(min(Nv), 5), mp.nstr(max(Nv), 5), mp.nstr(rr[0], 15), mp.nstr(max(abs(r / rr[0] - 1) for r in rr), 3)))
        print('   first rows (sector, N, |J6|, detN):')
        for r in recs[:4]:
            print('    ', r['sector'], mp.nstr(mp.mpf(r['N']), 10), mp.nstr(abs(mp.mpf(r['J6'])), 10), mp.nstr(mp.mpf(r['detN']), 10))


def stage_D1():
    import numpy as np
    hdr('D1: candidates (i)-(iii): log-linear regressions; (ii) Jacobian of the no-log conditions')
    for tag in ('P', 'X1'):
        pt, sj, S = load(tag)
        dd = T.detdata(pt)
        rows = []
        for k, (s, J, st) in enumerate(zip(S, sj, pt['states'])):
            f = {'detN': dd['dets'][k], 'M1': s['mons'][0][0], 'M2': s['mons'][1][0]}
            for si in range(2):
                for j in range(3):
                    f['P%d_%d' % (j + 1, si + 1)] = J.P[si][j][0]
                    f['Q%d_%d' % (j + 1, si + 1)] = J.Q[si][j][0]
            F, Jm = jac_at(st, pt['zeta'])
            f['J6'] = mp.det(Jm)
            rows.append((f, s, st))
        y = np.array([float(mp.log(s['N'])) for f, s, st in rows])

        def fit(fs, extra=False):
            A = np.array([[1.0] + ([st['sector'][0], st['sector'][1]] if extra else []) + [float(mp.log(abs(f[n]))) for n in fs] for f, s, st in rows])
            c, *_ = np.linalg.lstsq(A, y, rcond=None)
            return np.linalg.norm(A @ c - y)
        allq = ['detN', 'M1', 'M2'] + ['P%d_%d' % (j, s) for j in (1, 2, 3) for s in (1, 2)] + ['Q%d_%d' % (j, s) for j in (1, 2, 3) for s in (1, 2)]
        print('%s  regression of log N (64 states), residual norm (spread of log N: %.1f):' % (tag, np.linalg.norm(y - y.mean())))
        print('   (i)  detN_B alone %.2f; detN, M1, M2: %.2f' % (fit(['detN']), fit(['detN', 'M1', 'M2'])))
        print('   (i)+(iii) detN, M1, M2, P^j(theta_s), Q_j(theta_s) (15 params): %.2f; same + sector exponents mu1, mu2 (17 params): %.2f' % (fit(allq), fit(allq, True)))
        print('   (ii) J6 alone: %.2f ; J6 + sector exponents (4 params): %.2e' % (fit(['J6']), fit(['J6'], True)))
    print('\n(ii) exact check: c0 = N/(|J6| d23^mu1 d13^mu2 d12^mu3), state-independence over the 64 states, theta and zeta dependence')
    for tag in ('P', 'Pp', 'X1'):
        for th in ([TH, (Fr(1, 5), Fr(-4, 9)), (Fr(3, 7), Fr(-1, 4)), (Fr(1, 2), Fr(-1, 3)), (Fr(-1, 7), Fr(-6, 5))] if tag == 'P' else [TH]):
            pt, sj, S = load(tag, th)
            z = [mpf(x) for x in pt['zeta']]
            d12, d13, d23 = abs(z[0] - z[1]), abs(z[0] - z[2]), abs(z[1] - z[2])
            vals = []
            ns = 0
            for s, st in zip(S, pt['states']):
                F, Jm = jac_at(st, pt['zeta'], th)
                m = st['sector']
                dj = mp.det(Jm)
                ns += int(dj > 0)
                vals.append(s['N'] / abs(dj) / (d23 ** m[0] * d13 ** m[1] * d12 ** m[2]))
            print('  %-3s theta=%s  c0 = %s  spread over 64 states %s  (#det J6>0: %d)' % (tag, tuple(str(x) for x in th), mp.nstr(vals[0], 25), mp.nstr(max(abs(v / vals[0] - 1) for v in vals), 3), ns))
    p = os.path.join(HERE, 'gaudin_norm_c0.pkl')
    if os.path.exists(p):
        d = pickle.load(open(p, 'rb'))
        worst = mp.mpf(0)
        for k, (z, c) in d.items():
            zz = [mpf(x) for x in z]
            worst = max(worst, abs(mp.mpf(c) * (abs(zz[0] - zz[1]) * abs(zz[0] - zz[2]) * abs(zz[1] - zz[2])) ** 6 / (mp.mpf(9) / 4) - 1))
        print('c0 versus zeta: %d zeta points (8 named + 26 random), max |c0 (d12 d13 d23)^6/(9/4) - 1| = %s' % (len(d), mp.nstr(worst, 3)))


def stage_G():
    hdr('G / D1(iv): the SoV measure Gamma = Yinv^T (G1 (x) G1)_c Yinv; norm identity; completeness')
    for tag in ('P', 'Pp', 'X1'):
        pt, sj, S = load(tag)
        G = Gamma_mp(pt['zeta'])
        e = 0
        for s in S:
            nu = nu2([x / s['mons'][0][0] for x in s['mons'][0]], [x / s['mons'][1][0] for x in s['mons'][1]])
            e = max(e, abs(mp.fdot(nu, G * mp.matrix(nu)) / s['N'] - 1))
        print('%s: max |nu Gamma nu / N - 1| over 64 states = %s' % (tag, mp.nstr(e, 3)))
    for th in [(Fr(1, 5), Fr(-4, 9)), (Fr(1, 2), Fr(-1, 3)), (Fr(3, 7), Fr(-1, 4)), (Fr(-1, 7), Fr(-6, 5)), (Fr(4, 5), Fr(-2, 3))]:
        pt, sj, S = load('P', th)
        G = Gamma_mp(pt['zeta'], 1 / (th[0] - th[1]))
        e = 0
        for s in S:
            nu = nu2([x / s['mons'][0][0] for x in s['mons'][0]], [x / s['mons'][1][0] for x in s['mons'][1]])
            e = max(e, abs(mp.fdot(nu, G * mp.matrix(nu)) / s['N'] - 1))
        print('P, theta=%s: identity with Y(t) = t^L Y(1) t^-L: %s' % (tuple(str(x) for x in th), mp.nstr(e, 3)))
    pt, sj, S = load('P')
    Y, Yinv = comp_Y()
    Gi = B.frac_inverse(B.frac_matmul([list(r) for r in zip(*Yinv)], B.frac_matmul(Gc_exact('P')[1], Yinv)))
    nus = [nu2([x / s['mons'][0][0] for x in s['mons'][0]], [x / s['mons'][1][0] for x in s['mons'][1]]) for s in S]
    err = 0
    for a in range(64):
        for b in range(64):
            v = mp.fsum(nus[k][a] * nus[k][b] / S[k]['N'] for k in range(64))
            err = max(err, abs(v - mpf(Gi[a][b])))
    print('completeness: max |sum_A nu_a nu_b / N_A - (Gamma^-1)_ab| = %s' % mp.nstr(err, 3))


def stage_D2():
    hdr('D2: untwisting limit, regime caches, phi and phi2')
    import untwist_limit as UL
    for tag, phi in (('', UL.PHI), ('2', UL.PHI2)):
        print('### phi =', phi)
        cg_table(phi, tag)


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--stages', default='D0,D1,G,D2')
    a = ap.parse_args()
    for st in a.stages.split(','):
        {'D0': stage_D0, 'D1': stage_D1, 'G': stage_G, 'D2': stage_D2}[st]()
