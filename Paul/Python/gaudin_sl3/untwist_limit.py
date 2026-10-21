"""U2: analytic untwisting limit of the Gaudin jet formula (adjoint). Stage flags: python -X utf8 untwist_limit.py --stages R,S,...
R: build 12 regime points at 200 digits (cache untwist_limit_reg_<k>.pkl). Reuses twovac_fit (read-only)."""
import sys, os, time, pickle, argparse
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import mpmath as mp
import twovac_fit as T
import qside as Q
from gaudin_chain import *

DPS = 200
PHI = T.PHI
PHI2 = (Fr(-2, 5), Fr(3, 8), Fr(7, 11))
EPS = [Fr(1, 10 ** 6) / 2 ** k for k in range(12)]
GRID1 = [Fr(3, 4), Fr(-5, 3), Fr(1, 2), Fr(-7, 3)]
GRID2 = [Fr(-12, 5), Fr(7, 2), Fr(1, 8), Fr(5, 4)]
SECT = [(4, 2, 0), (4, 1, 1), (3, 3, 0), (3, 2, 1), (3, 2, 1), (2, 2, 2)]
BNAME = ['(4,2,0)', '(4,1,1)', '(3,3,0)', '(3,2,1)a', '(3,2,1)b', '(2,2,2)']


def build_regime(k, phi=PHI, dps=DPS, tag=''):
    eps = EPS[k]
    cache = os.path.join(HERE, 'untwist_limit_reg%s_%d.pkl' % (tag, k))
    if os.path.exists(cache):
        with open(cache, 'rb') as f:
            out = pickle.load(f)
        mp.mp.dps = dps
        return out
    T.set_precision(dps)
    zeta = tuple(eps * p for p in phi)
    Q.set_regime(list(zeta), Q.OM, Q.OMP, dps, 60)
    mp.mp.dps = dps
    ch = T.chain()
    g = [[zeta[i] if i == j else Fr(0) for j in range(3)] for i in range(3)]
    ops = Operators(ch, g)
    C = ops.fit()
    st = eigenstates(ch, ops, C, dps=dps)
    for x in st:
        x['ev'] = T.real(x['ev'], mp.mpf(10) ** -(dps - 30))
        x['vec'] = [T.real(c, mp.mpf(10) ** -(dps - 30)) for c in x['vec']]
    res = eigenfunctions(ch, C, st)
    thm = [T.mpf(t) for t in ch.theta]
    bysec = {}
    for kk, s in enumerate(st):
        s['rr'] = Q.raising_residual(ch, s)
        bysec.setdefault(s['sector'], []).append(kk)
    br, gaps = [], {}
    for mu, m in T.HW_SECT.items():
        ks = sorted(bysec[mu], key=lambda kk: st[kk]['rr'])
        sel = sorted(ks[:m], key=lambda kk: st[kk]['ev'])
        br += sel
        gaps[mu] = [st[kk]['rr'] for kk in ks[:m + 1]]
    S1 = Q.pf_S1(sum(T.mpf(z) for z in zeta), 3, thm)
    recs = []
    for kk in br:
        s = st[kk]
        S2 = Q.pf_S2(s['lam2'], thm)
        S3 = Q.pf_S3(s['lam3'], thm)
        dual = [Q.solve_poly_mp('D', S1, S2, S3, T.mpf(zeta[j]), thm) for j in range(3)]
        recs.append(dict(sector=s['sector'], ev=s['ev'], lam2=[T.real(x) for x in s['lam2']], lam3=[T.real(x) for x in s['lam3']], vec=s['vec'], rr=s['rr'],
                         dual=[[T.real(c) for c in x[3]] for x in dual], degs=[x[0] for x in dual], nul=[x[1] for x in dual]))
    out = dict(eps=eps, zeta=zeta, dps=dps, branches=recs, gaps=gaps, eigres=res)
    with open(cache, 'wb') as f:
        pickle.dump(out, f)
    return out


def stage_R(ks):
    for k in ks:
        t0 = time.time()
        o = build_regime(k)
        print('k=%d eps=%s built in %.0fs eigres %s gaps %s' % (k, o['eps'], time.time() - t0, mp.nstr(o['eigres'], 3),
              {m: [mp.nstr(x, 3) for x in g] for m, g in o['gaps'].items()}), flush=True)


import sympy as sp

# ------------------------------------------------------------------------------------------------ jets / series
def branch_N(b, zeta, th=T.THETA):
    """N^(s) = (P3 p02, P3' p02, P3 p03, P3' p03, P3 p23) at theta_s for a branch record (dual Q's only)."""
    z = [T.mpf(x) for x in zeta]
    t = [T.mpf(x) for x in th]
    out = []
    for s in range(2):
        P = [T.jets(z[j], b['dual'][j], t[s], 4) for j in range(3)]
        pl = T.plucker(P, 1, 2, 3)
        out.append([P[2][0] * pl[(0, 2)], P[2][1] * pl[(0, 2)], P[2][0] * pl[(0, 3)], P[2][1] * pl[(0, 3)], P[2][0] * pl[(2, 3)]])
    return out


def load_all(phi=PHI, tag='', dps=DPS):
    mp.mp.dps = dps
    return [build_regime(k, phi, dps, tag) for k in range(12)]


def series_fit(vals, npts=12, nmin=0):
    """fit f(eps_k) = sum_{n=nmin}^{nmin+npts-1} c_n eps_k^n through the first npts points (x_k = 2^-k scaling for conditioning); returns c_n (true eps)"""
    e0 = mp.mpf(EPS[0].numerator) / EPS[0].denominator
    xs = [mp.mpf(2) ** (-k) for k in range(npts)]
    A = mp.matrix(npts, npts)
    for i, x in enumerate(xs):
        for j in range(npts):
            A[i, j] = x ** (nmin + j)
    rhs = mp.matrix([vals[k] for k in range(npts)])
    c = mp.lu_solve(A, rhs)
    return [c[j] / e0 ** (nmin + j) for j in range(npts)]


# ------------------------------------------------------------------------------------------------ exact chain-side objects
def sym_E(ch):
    out = {}
    for key, sm in ch.E1.items():
        M = sp.zeros(ch.d1, ch.d1)
        for i, row in sm.r.items():
            for j, v in row.items():
                M[i, j] = sp.Rational(v.numerator, v.denominator)
        out[key] = M
    return out


def expm_nil(M, x, n=5):
    out = sp.eye(M.shape[0]); term = sp.eye(M.shape[0])
    for k in range(1, n + 1):
        term = term * M * x / k
        out = out + term
    return out


W1S, W2S = sp.symbols('W1 W2')


def gauss_closed(phi):
    p1, p2, p3 = [sp.Rational(x.numerator, x.denominator) for x in phi]
    a = (p1 - p2) * (W1S - W2S) / ((p2 - W1S) * (p2 - W2S))
    c = (p1 - p3) * (p2 - p3) * (W1S - W2S) / ((p1 - p2) * (p3 - W1S) * (p3 - W2S))
    g = (p3 - p1) * (p3 - p2) * (W1S - W2S) / ((p3 - W1S) * (p2 - W2S) * (p3 - W2S))
    return a, g, c


def chat_expr(phi):
    p1, p2, p3 = [sp.Rational(x.numerator, x.denominator) for x in phi]
    f = lambda w: (p2 - w) * (p3 - w) ** 2
    return f(W1S) * f(W2S) * (p1 - p2) ** 2


_GFC = {}
def gf_polys(ch, phi):
    """chat * [e1^T rho(n')]_nu as sympy polynomials in (W1, W2), nu = 0..7 (gl3_reps order). Also chat."""
    key = tuple(phi)
    if key in _GFC:
        return _GFC[key]
    E = sym_E(ch)
    a, g, c = gauss_closed(phi)
    ch_ = chat_expr(phi)
    e1 = sp.zeros(1, ch.d1); e1[0, ch.hw_index] = 1
    A, Gp, C = sp.symbols('A Gp C')
    v = e1 * expm_nil(E[(1, 2)], A) * expm_nil(E[(1, 3)], Gp) * expm_nil(E[(2, 3)], C)
    out = []
    for nu in range(ch.d1):
        ex = sp.together(sp.expand(v[0, nu]).subs({A: a, Gp: g, C: c}, simultaneous=True) * ch_)
        ex = sp.cancel(ex)
        num, den = sp.fraction(ex)
        out.append(sp.Poly(ex, W1S, W2S) if not den.has(W1S, W2S) else None)   # None: chat does not clear this component (weight never occurs in an octet/adjoint slice)
    _GFC[key] = (out, sp.Poly(sp.expand(ch_), W1S, W2S))
    return _GFC[key]


def poly_dict(P):
    return {m: sp.Rational(c) for m, c in zip(P.monoms(), P.coeffs())}


def mpq(q):
    q = sp.Rational(q)
    return mp.mpf(int(q.p)) / int(q.q)

# ------------------------------------------------------------------------------------------------ chain-side branch data
def branch_psi(pts, bi, ch):
    """G-normalised branch vectors at the 12 eps, sign made continuous (aligned with the previous/smaller eps)."""
    vs = []
    for p in pts:
        vs.append(T.gnormalise(ch, p['branches'][bi]['vec']))
    # sign continuity: align each with the next (smaller eps) one
    for k in range(10, -1, -1):
        d = mp.fsum(a * b for a, b in zip(vs[k], vs[k + 1]))
        if d < 0:
            vs[k] = [-x for x in vs[k]]
    return vs


def comp_series(vs, npts=12):
    """per-component eps series coefficients Psi^(n), n = 0..npts-1 (list of vectors)"""
    nc = len(vs[0])
    cols = [series_fit([v[j] for v in vs], npts) for j in range(nc)]
    return [[cols[j][n] for j in range(nc)] for n in range(npts)]


def lead_exp(c, tol=mp.mpf(10) ** -25, ref=None):
    ref = ref if ref is not None else max(abs(x) for x in c[:9])
    for n, x in enumerate(c):
        if abs(x) > tol * ref:
            return n
    return None


# ------------------------------------------------------------------------------------------------ kernel pieces, series, leading polynomial
_GP = {}
def G_pieces(thi=None):
    """spec kernel Gs = G_D2/4 split into (omega1, omega1')-degree pieces: {(i,j): {(a,b): mpf}} (monomial exponents of (w1, w1p); thi numeric)"""
    th = T.THETA
    thi = Fr(1) / (th[0] - th[1]) if thi is None else thi
    if thi in _GP:
        return _GP[thi]
    T.d2_G(Fr(1), Fr(1), thi)
    w, wp, t = sp.symbols('w1 w1p thi')
    out = {}
    for (i, j), ex in T._G2.items():
        P = sp.Poly(sp.expand(ex.subs(t, sp.Rational(thi.numerator, thi.denominator))), w, wp)
        out[(i, j)] = {m: mpq(c) / 4 for m, c in zip(P.monoms(), P.coeffs())}
    _GP[thi] = out
    return out


def N_series(pts, bi, nser=12):
    """Nser[s][i] = list of eps-series coefficients N_i^{(n)}(theta_s), n = 0..nser-1"""
    Ns = [branch_N(p['branches'][bi], p['zeta']) for p in pts]
    return [[series_fit([N[s][i] for N in Ns], nser) for i in range(5)] for s in range(2)]


def leading_poly(Nser, tmax, G=None):
    """R[t] = {(a,b): coeff} for t = total eps order = d + n + m, t <= tmax: sum_ij G_ij^[d] N_i^(n)(th1) N_j^(m)(th2); also triples contributions"""
    G = G or G_pieces()
    R = {}
    trip = {}
    for (i, j), pc in G.items():
        for (a, b), gc in pc.items():
            d = a + b
            for n in range(0, tmax - d + 1):
                for m in range(0, tmax - d - n + 1):
                    t = d + n + m
                    v = gc * Nser[0][i][n] * Nser[1][j][m]
                    R.setdefault(t, {})
                    R[t][(a, b)] = R[t].get((a, b), mp.mpf(0)) + v
                    trip.setdefault(t, {})
                    key = (d, n, m)
                    trip[t][key] = trip[t].get(key, mp.mpf(0)) + abs(v)
    return R, trip



# ------------------------------------------------------------------------------------------------ the analysis
def pe(P, w1, w2):
    """evaluate dict-polynomial {(a,b): mp} at rational (w1, w2)"""
    w1, w2 = T.mpf(w1), T.mpf(w2)
    return mp.fsum(c * w1 ** a * w2 ** b for (a, b), c in P.items())


def sdict(Pd, w1, w2):
    w1, w2 = sp.Rational(w1.numerator, w1.denominator), sp.Rational(w2.numerator, w2.denominator)
    return sum(c * w1 ** a * w2 ** b for (a, b), c in Pd.items())


def analyse(phi=PHI, tag='', verbose=True):
    pts = load_all(phi, tag)
    ch = T.chain()
    mp.mp.dps = DPS
    Hs = {mu: [[T.mpf(x) for x in h] for h in T.hw_space(ch, mu)] for mu in T.HW_SECT}
    Fpolys, chat = gf_polys(ch, phi)
    Fd = [poly_dict(P) if P is not None else {} for P in Fpolys]
    nonpoly = [nu for nu, P in enumerate(Fpolys) if P is None]
    xphi = T.x0_D2(ch, phi)                                   # x(phi) (x) x(phi) (mp list); x0(eps phi) = eps^2 x0(phi)
    # eps-free vacua
    covs = {}
    for a1 in GRID1:
        for a2 in GRID2:
            covs[(a1, a2)] = T.cov_to_mp(T.cov_inf(ch, phi, a1, a2))
    # check eps-freeness (exact) once
    e0 = EPS[0]
    c1 = T.cov_inf(ch, tuple(e0 * p for p in phi), e0 * GRID1[0], e0 * GRID2[0])
    c2 = T.cov_inf(ch, phi, GRID1[0], GRID2[0])
    epsfree = (c1 == c2)
    res = dict(epsfree=epsfree, br=[])
    for bi in range(6):
        mu = SECT[bi]
        vs = branch_psi(pts, bi, ch)
        PS = comp_series(vs)
        psi0 = PS[0]
        R = {}
        # <x0|Psi_B> series
        xs = series_fit([T.overlap(T.x0_D2(ch, p['zeta']), vs[k]) for k, p in enumerate(pts)])
        kx = lead_exp(xs)
        xphiPsi_n = [T.overlap(xphi, PS[n]) for n in range(6)]
        # M1M2
        Ns = [branch_N(p['branches'][bi], p['zeta']) for p in pts]
        M12 = series_fit([N[0][0] * N[1][0] for N in Ns])
        kM = lead_exp(M12)
        # Phi overlaps
        phiexp = []
        phirat = []
        C = [psi0[nu * 8 + 7] for nu in range(8)]
        if any(abs(C[nu]) > mp.mpf(10) ** -40 for nu in nonpoly):
            print('WARNING: nonzero slice on a component chat does not clear', nonpoly)
        for key, cov in covs.items():
            s_ = series_fit([T.overlap(cov, v) for v in vs])
            kP = lead_exp(s_)
            phiexp.append(kP)
            gf = mp.fsum(C[nu] * pe(Fd[nu], *key) for nu in range(8)) / pe(poly_dict(chat), *key)
            phirat.append(abs(s_[0] / gf - 1) if kP == 0 else None)
            cD = T.overlap(cov, psi0)
        # N series and exponents
        Nser = [[series_fit([N[s][i] for N in Ns]) for i in range(5)] for s in range(2)]
        nexp = [[lead_exp(Nser[s][i]) for i in range(5)] for s in range(2)]
        Rt, trip = leading_poly(Nser, 11)
        order = [(t, max(abs(v) for v in Rt[t].values())) for t in sorted(Rt)]
        NB = next(t for t, v in order if v > mp.mpf(10) ** -30)
        LB = Rt[NB]
        # exact GF with Psi0 = chain limit
        GFd = {}
        for nu in range(8):
            for m_, c_ in Fd[nu].items():
                GFd[m_] = GFd.get(m_, mp.mpf(0)) + C[nu] * mpq(c_)
        mons = sorted(set(GFd) | set(LB))
        l = mp.matrix([LB.get(m_, mp.mpf(0)) for m_ in mons])
        f = mp.matrix([GFd.get(m_, mp.mpf(0)) for m_ in mons])
        RB = (l.T * f)[0] / (f.T * f)[0]
        resid = mp.norm(l - RB * f) / mp.norm(l)
        pred = M12[kM] / xs[kx]
        # dictionary: solve for Psi0 in H_mu from L_B
        H = Hs[mu]
        Amat = mp.matrix(len(mons), len(H))
        for i, h in enumerate(H):
            for nu in range(8):
                for m_, c_ in Fd[nu].items():
                    Amat[mons.index(m_) if m_ in mons else 0, i] += 0 if m_ not in mons else h[nu * 8 + 7] * mpq(c_)
        beta = mp.lu_solve(Amat.T * Amat, Amat.T * l)
        vq = [mp.fsum(beta[i] * H[i][j] for i in range(len(H))) for j in range(64)]
        vqn = T.gnormalise(ch, vq)
        n2 = mp.sqrt(mp.fsum(T.mpf(T.gram_diag(ch)[a // 8] * T.gram_diag(ch)[a % 8]) * vq[a] ** 2 for a in range(64)))
        sgn = 1 if (vq[[NB2G_ for NB2G_ in [T.NB2G[a] * 8 + T.NB2G[b] for a in range(8) for b in range(8)] if abs(vq[NB2G_]) > mp.mpf(10) ** -30][0]] > 0) else -1
        dev = max(abs(x - y) for x, y in zip(vqn, psi0))
        slq = [vqn[nu * 8 + 7] for nu in range(8)]
        trs = sorted(trip[NB].items(), key=lambda kv: -kv[1])
        rec = dict(bi=bi, mu=mu, kx=kx, kM=kM, xs=xs[:8], M12=M12[:8], phiexp=phiexp, phirat=phirat, nexp=nexp, NB=NB, order=order, RB=RB, resid=resid, pred=pred,
                   dev=dev, slice_chain=C, slice_q=slq, trip=[(k, v) for k, v in trs], xphiPsi_n=xphiPsi_n, Nser=Nser, LB=LB, mons=mons, beta=beta, scaleQ=n2 * sgn)
        res['br'].append(rec)
        if verbose:
            print('=== %s: k_B=%s (<x0|Psi> ~ eps^%s, expected 2+k_B), M1M2 ~ eps^%s, Phi exps (16 pts) %s, max|Phi lead / GF(Psi0) - 1| %s' % (
                BNAME[bi], kx - 2, kx, kM, sorted(set(phiexp)), mp.nstr(max(x for x in phirat if x is not None), 3)))
            print('   <x(phi)x(phi)|Psi_n> n=0..5: %s' % [mp.nstr(abs(x), 3) for x in xphiPsi_n])
            print('   N_i exponents site1 %s site2 %s ; N_B = %d (order table max|R_t|: %s)' % (nexp[0], nexp[1], NB, [(t, mp.nstr(v, 2)) for t, v in order[:NB + 2]]))
            print('   identity L_B = R_B * chat*GF(Psi0): R_B = %s, relative residual %s; predicted R_B = m_B/X_B = %s, rel diff %s' % (mp.nstr(RB, 12), mp.nstr(resid, 3), mp.nstr(pred, 12), mp.nstr(abs(RB / pred - 1), 3)))
            print('   dictionary: Q-side G-normalised vector vs chain Psi0: max diff %s; slice(Q) %s' % (mp.nstr(dev, 3), [mp.nstr(x, 10) for x in slq if abs(x) > 1e-20]))
            print('   contributing (d,n,m) at order N_B (by |term|): %s' % [(k, mp.nstr(v, 3)) for k, v in trs[:12]])
    return res



# ------------------------------------------------------------------------------------------------ nu-form (Laurent data only)
import untwist_limit_frob as FR
_FNU = None
def nu_funcs():
    global _FNU
    if _FNU is None:
        r = FR.solve_nu()
        _FNU = {k: sp.lambdify((FR.s1_, FR.e2, FR.e3, FR.Hs, FR.Ht, FR.As, FR.At, FR.Bs, FR.Bt, FR.dl), r[k], 'mpmath') for k in ('kappa', 'pi', 'eps')}
    return _FNU


def nu_from_laurent(b, zeta, th=T.THETA):
    """nu_s = (1, kappa, pi, kappa pi, eps)(theta_s) from the Laurent data (lam2, lam3) of a branch and zeta only (Frobenius at theta_s of the dual oper)"""
    f = nu_funcs()
    z = [T.mpf(q) for q in zeta]
    e2v = z[0] * z[1] + z[0] * z[2] + z[1] * z[2]; e3v = z[0] * z[1] * z[2]; sg = sum(z)
    l2, l3 = b['lam2'], b['lam3']
    out = []
    for s in range(2):
        Hs, Ht = (l2[1], l2[3]) if s == 0 else (l2[3], l2[1])
        As, At = (l3[2], l3[5]) if s == 0 else (l3[5], l3[2])
        Bs, Bt = (l3[1], l3[4]) if s == 0 else (l3[4], l3[1])
        d = T.mpf(th[s] - th[1 - s])
        args = (sg, e2v, e3v, Hs, Ht, As, At, Bs, Bt, d)
        k_, p_ = f['kappa'](*args), f['pi'](*args)
        out.append([mp.mpf(1), k_, p_, k_ * p_, f['eps'](*args)])
    return out


def nu_series(pts, bi):
    Ns = [nu_from_laurent(p['branches'][bi], p['zeta']) for p in pts]
    return [[series_fit([N[s][i] for N in Ns]) for i in range(5)] for s in range(2)]



# ------------------------------------------------------------------------------------------------ nu-form analysis + extra checks
def analyse_nu(phi=PHI, tag='', verbose=True):
    pts = load_all(phi, tag)
    ch = T.chain()
    Hs = {mu: [[T.mpf(x) for x in h] for h in T.hw_space(ch, mu)] for mu in T.HW_SECT}
    Fpolys, chat = gf_polys(ch, phi)
    Fd = [poly_dict(P) if P is not None else {} for P in Fpolys]
    chatd = {m_: mpq(c_) for m_, c_ in poly_dict(chat).items()}
    xphi = T.x0_D2(ch, phi)
    # exact test: x(phi) (x) x(phi) on exact HW kets
    xtest = {}
    for mu, Hl in Hs.items():
        xtest[mu] = [abs(T.overlap(xphi, h)) for h in Hl]
    G = G_pieces()
    out = []
    for bi in range(6):
        mu = SECT[bi]
        vs = branch_psi(pts, bi, ch); PS = comp_series(vs)
        psi0 = PS[0]
        C = [psi0[nu * 8 + 7] for nu in range(8)]
        xn = [T.overlap(xphi, PS[n]) for n in range(8)]
        kB = next(n for n in range(8) if abs(xn[n]) > 1e-25)
        XB = xn[kB]
        nus = nu_series(pts, bi)
        NBp = 6 - kB
        R, trip = leading_poly(nus, NBp + 1)
        LB = R[NBp]
        GFd = {}
        for nu in range(8):
            for m_, c_ in Fd[nu].items():
                GFd[m_] = GFd.get(m_, mp.mpf(0)) + C[nu] * mpq(c_)
        mons = sorted(set(GFd) | set(LB))
        ident = max(abs(LB.get(m_, 0) - GFd.get(m_, 0) / XB) for m_ in mons) / max(abs(v) for v in LB.values())
        low = max([max(abs(v) for v in R[t].values()) for t in R if t < NBp] + [mp.mpf(0)])
        # rho_B at the 16 grid points relative to (3/4, -12/5)
        w0 = (GRID1[0], GRID2[0])
        rho_dev = mp.mpf(0)
        for a1 in GRID1:
            for a2 in GRID2:
                r_nu = (pe(chatd, *w0) / pe(chatd, a1, a2)) * pe(LB, a1, a2) / pe(LB, *w0)
                r_ex = (sum(C[nu] * pe(Fd[nu], a1, a2) for nu in range(8)) / pe(chatd, a1, a2)) / (sum(C[nu] * pe(Fd[nu], *w0) for nu in range(8)) / pe(chatd, *w0))
                rho_dev = max(rho_dev, abs(r_nu / r_ex - 1))
        # dictionary: solve for the HW ket from L'_B
        H = Hs[mu]
        Amat = mp.matrix(len(mons), len(H))
        for i, h in enumerate(H):
            for nu in range(8):
                for m_, c_ in Fd[nu].items():
                    Amat[mons.index(m_), i] += h[nu * 8 + 7] * mpq(c_)
        l = mp.matrix([LB.get(m_, mp.mpf(0)) for m_ in mons])
        beta = mp.lu_solve(Amat.T * Amat, Amat.T * l)
        vq = [mp.fsum(beta[i] * H[i][j] for i in range(len(H))) for j in range(64)]
        vqn = T.gnormalise(ch, vq)
        dev = max(abs(x - y) for x, y in zip(vqn, psi0))
        slq = [vqn[nu * 8 + 7] for nu in range(8)]
        d_ = T.gram_diag(ch)
        n2 = mp.sqrt(mp.fsum(T.mpf(d_[a // 8] * d_[a % 8]) * vq[a] ** 2 for a in range(64)))
        rec = dict(bi=bi, kB=kB, NBp=NBp, XB=XB, ident=ident, low=low, rho_dev=rho_dev, dev=dev, slice_q=slq, slice_chain=C, beta_norm=n2, xn=[abs(x) for x in xn[:7]],
                   trip=sorted(trip[NBp].items()), xtest=xtest[mu], LB=LB)
        out.append(rec)
        if verbose:
            print('%s: k_B=%d  N_B-2=%d  ||Sum G nu nu||: lower orders <= %s; L\'_B vs chat*GF(Psi0)/X_B max rel coeff dev %s; rho_B dev (16 pts) %s; dictionary |v_Q - Psi0| = %s' % (
                BNAME[bi], kB, NBp, mp.nstr(low, 3), mp.nstr(ident, 3), mp.nstr(rho_dev, 3), mp.nstr(dev, 3)))
            print('    X_B = <x x|Psi_kB> = %s ; |<xx|H_mu>| exact = %s ; slice(Q) %s' % (mp.nstr(XB, 12), [mp.nstr(v, 2) for v in xtest[mu]], [(nu, mp.nstr(slq[nu], 12)) for nu in range(8) if abs(slq[nu]) > 1e-20]))
    return out

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--stages', default='R,N,NU')
    a = ap.parse_args()
    st = a.stages.split(',')
    if 'R' in st:
        for tag, phi in (('', PHI), ('2', PHI2)):
            for k in range(12):
                build_regime(k, phi, DPS, tag)
        print('regime points built (200 digits, 12 eps, phi and phi2)')
    for tag, phi in (('', PHI), ('2', PHI2)):
        print('######## phi = %s' % (phi,))
        if 'N' in st:
            print('--- N-form (jet monomials N_i): leading polynomial L_B, identity L_B = R_B chat GF')
            analyse(phi, tag)
        if 'NU' in st:
            print('--- nu-form (nu = N/M from Laurent data): sum G nu nu = eps^(6-k_B) chat GF / X_B')
            analyse_nu(phi, tag)
