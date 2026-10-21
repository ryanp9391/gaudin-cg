"""D1: direct fit of the Q-side formula for the CG-complete two-vacua overlap of the sl(3) Gaudin chain (L=2, Lambda=(2,1,0)).

Run top to bottom:  python -X utf8 twovac_fit.py [--stages F0,F1,F2,F3,F4,F5]   (caches in twovac_fit_*.pkl)
Basis order everywhere: gl3_reps order (HW ket index 7), the order of gaudin_chain/qside (NOT the notebook order of xxx_limit).
"""
import sys, os, time, itertools, pickle, argparse
from fractions import Fraction as Fr
from math import factorial
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mpmath as mp
import sympy
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication
import qside as Q
from gaudin_chain import *          # Chain, Operators, vacuum_covector, gauss_r, exp_nil, SM, kron ...

DPS = 80
TOLQ = 50
mp.mp.dps = DPS

THETA = (Fr(2, 7), Fr(-3, 5))
ZETAS = {
    'P':  (Fr(3, 2), Fr(5, 4), Fr(2, 3)),
    'Pp': (Fr(2), Fr(5, 3), Fr(-1, 2)),
    'X1': (Fr(7, 5), Fr(-3, 4), Fr(1, 6)),
    'X2': (Fr(-4, 3), Fr(2, 9), Fr(5, 7)),
    'X3': (Fr(1, 3), Fr(-2, 7), Fr(5, 11)),
    'X4': (Fr(5, 3), Fr(-1, 2), Fr(3, 8)),
    'X5': (Fr(-3, 5), Fr(13, 7), Fr(2)),
    'X6': (Fr(4, 5), Fr(-5, 6), Fr(3, 4)),
    'X7': (Fr(-1, 4), Fr(6, 5), Fr(9, 7)),
    'X8': (Fr(3, 7), Fr(-8, 5), Fr(1, 2)),
    'X9': (Fr(5, 6), Fr(1, 7), Fr(-9, 4)),
}
OM0 = (Fr(5, 2), Fr(7, 3), Fr(3, 5))


def _candidate_zetas(n, seed=12345):
    """deterministic candidate rational zeta-points (distinct entries, separated by >= 0.15, |zeta| <= 3, avoiding every omega/theta value used)"""
    import random
    rnd = random.Random(seed)
    avoid = {Fr(5, 2), Fr(-2, 5), Fr(9, 7), Fr(7, 4), Fr(3), Fr(-5, 3), Fr(7, 3), Fr(3, 5), Fr(2, 7), Fr(-3, 5), Fr(1, 5), Fr(-4, 9), Fr(3, 7), Fr(-1, 4)}
    out = []
    while len(out) < n:
        z = tuple(Fr(rnd.randint(-17, 17), rnd.choice([2, 3, 4, 5, 6, 7, 8, 9, 11])) for _ in range(3))
        if any(abs(a - b) < Fr(15, 100) for a in z for b in z if a != b) or len(set(z)) < 3 or any(a in avoid for a in z) or any(abs(a) > 3 for a in z):
            continue
        out.append(z)
    return out


def _spectrum_is_real(z):
    ch = Chain(theta=THETA)
    g = [[z[i] if i == j else Fr(0) for j in range(3)] for i in range(3)]
    ops = Operators(ch, g)
    C = ops.fit()
    mp.mp.dps = 30
    st = eigenstates(ch, ops, C, dps=30)
    return all(abs(mp.im(s['ev'])) < mp.mpf(10) ** -20 for s in st)


def extra_zetas(n=60):
    """tags 'R1'.. -> zeta triples with a real S2+pi*S3 spectrum (cached)"""
    path = os.path.join(HERE, 'twovac_fit_zetalist.pkl')
    if os.path.exists(path):
        with open(path, 'rb') as f:
            return pickle.load(f)
    out = {}
    for z in _candidate_zetas(400):
        if _spectrum_is_real(z):
            out['R%d' % (len(out) + 1)] = z
        if len(out) >= n:
            break
    with open(path, 'wb') as f:
        pickle.dump(out, f)
    return out


def zeta_of(tag):
    if tag in ZETAS:
        return ZETAS[tag]
    return extra_zetas()[tag]
W1BASE = (Fr(5, 2), Fr(7, 4))
XXX_TRIP = ((0, 2), (0, 1), (1, 2))


def mpf(x):
    if isinstance(x, Fr):
        return mp.mpf(x.numerator) / x.denominator
    return mp.mpf(x)


def hdr(t):
    print('\n' + '=' * 100 + '\n' + t + '\n' + '=' * 100, flush=True)


def nst(x, n=3):
    return mp.nstr(x, n)


def real(x, tol=None):
    tol = tol or mp.mpf(10) ** (-(DPS - 20))
    if isinstance(x, mp.mpc):
        if abs(x.imag) > tol * max(1, abs(x.real)):
            raise ValueError('complex value with imag %s' % nst(x.imag))
        return x.real
    return x


# ------------------------------------------------------------------------------------------------ chain / point data
_CH = {}
def chain(th=THETA):
    key = tuple(Fr(t) for t in th)
    if key not in _CH:
        _CH[key] = Chain(theta=key)
    return _CH[key]


def set_precision(dps, tolq=None):
    global DPS
    DPS = dps
    mp.mp.dps = dps


def build_point(zeta, dps=None, th=THETA):
    """states, dual/direct Q coefficient lists, Laurent data, for one zeta point. Plain (picklable) record list."""
    dps = dps or DPS
    zeta = tuple(Fr(z) for z in zeta)
    Q.set_regime(list(zeta), Q.OM, Q.OMP, dps, max(TOLQ, 35))
    mp.mp.dps = dps
    ch = chain(th)
    g = [[zeta[i] if i == j else Fr(0) for j in range(3)] for i in range(3)]
    ops = Operators(ch, g)
    thf = set(Fr(t) for t in th)
    pts = [Fr(3, 7), Fr(-2, 9), Fr(5, 3), Fr(7, 2), Fr(-4, 5), Fr(11, 4), Fr(-7, 3), Fr(13, 5), Fr(1, 3)]
    if any(u in thf for u in pts):
        pts = [Fr(17, 13), Fr(-19, 11), Fr(23, 17), Fr(29, 7), Fr(-31, 13), Fr(37, 19), Fr(-41, 17), Fr(43, 23), Fr(47, 29)]
    C = ops.fit(pts)
    u0, u1 = (Fr(3, 7), Fr(-2, 9)) if Fr(3, 7) not in thf and Fr(-2, 9) not in thf else (Fr(17, 13), Fr(-19, 11))
    st = eigenstates(ch, ops, C, u0=u0, u1=u1, dps=dps)
    for x in st:
        x['ev'] = real(x['ev'], mp.mpf(10) ** -(dps - 30))
        x['vec'] = [real(c, mp.mpf(10) ** -(dps - 30)) for c in x['vec']]
    res = eigenfunctions(ch, C, st)
    thm = [mpf(t) for t in ch.theta]
    recs = Q.q1_run(ch, st, thm)
    out = []
    for s, r in zip(st, recs):
        assert tuple(len(x[3]) - 1 for x in r['dual']) == tuple(s['sector'])
        out.append(dict(sector=s['sector'], ev=real(s['ev']), lam2=[real(x) for x in s['lam2']], lam3=[real(x) for x in s['lam3']],
                        vec=[real(x) for x in s['vec']], pivot=s['pivot'],
                        dual=[[real(c) for c in x[3]] for x in r['dual']], direct=[[real(c) for c in x[3]] for x in r['direct']],
                        resid=r['resid']))
    return dict(zeta=zeta, theta=tuple(Fr(t) for t in th), dps=dps, states=out, eigres=res)


def thtag(th):
    return '' if tuple(Fr(t) for t in th) == THETA else '_th%s' % '_'.join('%d.%d' % (Fr(t).numerator, Fr(t).denominator) for t in th)


def get_point(tag, zeta=None, dps=None, th=THETA):
    dps = dps or DPS
    path = os.path.join(HERE, 'twovac_fit_pt_%s_%d%s.pkl' % (tag, dps, thtag(th)))
    if os.path.exists(path):
        with open(path, 'rb') as f:
            pt = pickle.load(f)
        mp.mp.dps = dps
        return pt
    t0 = time.time()
    pt = build_point(zeta if zeta is not None else zeta_of(tag), dps, th)
    with open(path, 'wb') as f:
        pickle.dump(pt, f)
    print('  built point %s (%d digits) in %.1fs; eigen residual %s' % (tag, dps, time.time() - t0, nst(pt['eigres'])), flush=True)
    return pt


# ------------------------------------------------------------------------------------------------ Q objects, jets
def poly_derivs(kappa, q, n):
    """polynomials r_k with (d^k)(e^{kappa u} q) = e^{kappa u} r_k, k = 0..n"""
    out = [list(q)]
    cur = list(q)
    for _ in range(n):
        cur = Q.padd(Q.pder(cur), Q.pscale(cur, kappa))
        out.append(cur)
    return out


def jets(kappa, q, x, n=4):
    """[ (d/du)^k (e^{kappa u} q(u)) at u = x ], k = 0..n"""
    e = mp.exp(kappa * x)
    return [e * Q.peval(r, x) for r in poly_derivs(kappa, q, n)]


class StateJets:
    """jets at theta_s, s = 0,1 of P^j (j = 0,1,2), direct Q_i (i = 0,1,2) of one state, zeta point"""
    def __init__(self, st, zeta, th=THETA, n=4):
        z = [mpf(x) for x in zeta]
        t = [mpf(x) for x in th]
        self.z, self.t = z, t
        self.P = [[jets(z[j], st['dual'][j], t[s], n) for j in range(3)] for s in range(2)]       # [s][j][k]
        self.Q = [[jets(-z[i], st['direct'][i], t[s], n) for i in range(3)] for s in range(2)]    # [s][i][k]


def plucker(J, a=1, b=2, order=3):
    """p_kl = J[a][k] J[b][l] - J[b][k] J[a][l], 0 <= k < l <= order"""
    return {(k, l): J[a][k] * J[b][l] - J[b][k] * J[a][l] for k in range(order + 1) for l in range(k + 1, order + 1)}


# ------------------------------------------------------------------------------------------------ chain-side objects
def onesite_inf(ch, zeta, w1):
    """Omega^inf[zeta, w1] = e_hw^T rho(exp(r1 E12) exp(r2 E13) exp(r3 E23)), one site, Fractions"""
    z1, z2, z3 = [Fr(x) for x in zeta]
    w1 = Fr(w1)
    r1 = (w1 - z1) / (w1 - z2)
    r2 = -(w1 - z1) * (z2 - z3) / ((z1 - z2) * (w1 - z3))
    r3 = (w1 - z2) * (z1 - z3) / ((z1 - z2) * (w1 - z3))
    e = [Fr(0)] * ch.d1
    e[ch.hw_index] = Fr(1)
    v = exp_nil(ch.E1[(1, 2)], r1).rapply(e)
    v = exp_nil(ch.E1[(1, 3)], r2).rapply(v)
    v = exp_nil(ch.E1[(2, 3)], r3).rapply(v)
    return v


def cov_inf(ch, zeta, w1, w1p):
    a = onesite_inf(ch, zeta, w1)
    b = onesite_inf(ch, zeta, w1p)
    return [x * y for x in a for y in b]


def x0_onesite(ch, zeta):
    """x0^G one site = e_hw^T P_phys^{-1}, P_phys^{-1} = ppinv dinv pminv (formula of xxx_limit.mp_x0_onesite), gl3_reps basis order."""
    z1, z2, z3 = [mpf(x) for x in zeta]
    d = ch.d1
    E = {k: v.to_mp() for k, v in ch.E1.items()}
    nil = ch.lam[0] - ch.lam[2]
    def expN(c, key):
        A = E[key] * c
        out = mp.eye(d)
        Pw = mp.eye(d)
        for k in range(1, nil + 1):
            Pw = Pw * A
            out = out + Pw / mp.factorial(k)
        return out
    s = [z1 ** 2 / ((z1 - z2) * (z1 - z3)), -(z2 / (z1 * (z2 - z3))), 1 / (z1 * z2)]
    r = [-(z3 / z2), -(z3 / z1), (z2 ** 2 * (z3 - z1)) / (z1 ** 2 * (z2 - z3))]
    ppinv = expN(-(z1 + z2) / (z1 * z2), (2, 3)) * expN(1 / (z1 * z2), (1, 3)) * expN(-(1 / z1), (1, 2))
    wts = [[E[(k, k)][p, p] for k in (1, 2, 3)] for p in range(d)]
    dinv = mp.diag([mp.fprod([s[k] ** (-wts[p][k]) for k in range(3)]) for p in range(d)])
    pminv = expN(-r[2], (2, 1)) * expN(-r[1], (3, 1)) * expN(-r[0], (3, 2))
    e1 = mp.zeros(1, d)
    e1[0, ch.hw_index] = 1
    x = e1 * ppinv * dinv * pminv
    return [x[0, j] for j in range(d)]


def cov_to_mp(c):
    return [mpf(x) for x in c]


def overlap(cov, vec):
    return mp.fsum(c * v for c, v in zip(cov, vec) if c != 0)


def S2S3_labels(st, th=THETA):
    thm = [mpf(t) for t in th]
    S2 = Q.pf_S2(st['lam2'], thm)
    S3 = Q.pf_S3(st['lam3'], thm)
    return S2(mp.mpf(3) / 7), S3(mp.mpf(-2) / 9)


def detN_for_states(pt, om0=OM0, th=THETA):
    """P1's D_B(om0): det N with the vacuum direct Q's of om0 (exact), rows (s,(i,j)) in XXX_TRIP; list over states"""
    ch = chain(th)
    zeta = pt['zeta']
    Q.set_regime(list(zeta), list(om0), list(om0), pt['dps'], 35)
    mp.mp.dps = pt['dps']
    thm = [mpf(t) for t in th]
    vd = Q.vac_data(ch, list(om0))
    assert vd['ok']
    Qs = [Q.Qfun(-mpf(om0[i]), [mpf(c) for c in vd['Qd'][i][2]]) for i in range(3)]
    tdV = Q.SiteData(Qs, thm)
    z = [mpf(x) for x in zeta]
    dets = []
    for st in pt['states']:
        P = [Q.Qfun(z[j], st['dual'][j]) for j in range(3)]
        tdP = Q.SiteData(P, thm)
        N = Q.Nrows(tdP, tdV, thm)
        dets.append(mp.det(Q.build_N_det(N, XXX_TRIP)))
    return dets, vd['cov']


# ------------------------------------------------------------------------------------------------ linear algebra (mpmath)
def dot(a, b):
    return mp.fdot(zip(a, b))


FP = 320        # fixed-point bits (about 96 digits) of the integer-arithmetic QR


def _toint(x):
    return int(mp.floor(x * mp.mpf(2) ** FP + mp.mpf(1) / 2))


def _fromint(i):
    return mp.mpf(i) / mp.mpf(2) ** FP


class PivQR:
    """column-pivoted modified Gram-Schmidt (two passes) in fixed-point integer arithmetic (FP bits) of the m x n matrix with
    columns `cols` (lists of mpf); each column is scaled to max-abs 1 internally (scales are undone in the solutions)."""
    def __init__(self, cols, rtol=None):
        n = len(cols)
        self.n = n
        self.m = len(cols[0])
        self.scale = []
        W = []
        for c in cols:
            mx = max(abs(x) for x in c) or mp.mpf(1)
            self.scale.append(mx)
            W.append([_toint(x / mx) for x in c])
        one = 1 << FP
        def idot(a, b):
            return sum(x * y for x, y in zip(a, b)) >> FP
        self.idot = idot
        nrm2 = [idot(c, c) for c in W]
        self.perm = []
        self.Qv = []
        self.diag = []
        self.R = []
        used = set()
        first = None
        floor = 10 ** (-(mp.mp.dps - 10)) if rtol is None else rtol
        for k in range(min(n, self.m)):
            cand = [j for j in range(n) if j not in used]
            if not cand:
                break
            j = max(cand, key=lambda t: nrm2[t])
            if nrm2[j] <= 0:
                break
            rkk = math_isqrt(nrm2[j] << FP)
            if first is None:
                first = rkk
            if rkk * 10 ** 80 <= first:
                break
            used.add(j)
            q = [(x << FP) // rkk for x in W[j]]
            self.perm.append(j)
            self.Qv.append(q)
            self.diag.append(_fromint(rkk))
            row = {j: _fromint(rkk)}
            for j2 in range(n):
                if j2 in used:
                    continue
                w = W[j2]
                r = idot(q, w)
                w = [a - ((r * b) >> FP) for a, b in zip(w, q)]
                r2 = idot(q, w)
                w = [a - ((r2 * b) >> FP) for a, b in zip(w, q)]
                row[j2] = _fromint(r + r2)
                W[j2] = w
                nrm2[j2] = idot(w, w)
            self.R.append(row)
        self.k = len(self.perm)
        self.rest = [j for j in range(n) if j not in used]

    def spectrum(self):
        return [d / self.diag[0] for d in self.diag]

    def rank(self, tol=mp.mpf(10) ** -40):
        return sum(1 for d in self.diag if d > tol * self.diag[0])

    def lstsq(self, y, tol=mp.mpf(10) ** -40):
        """basic least-squares solution in the ORIGINAL column scaling (non-basic unknowns 0), relative residual; (x, resid, rank)"""
        r = self.rank(tol)
        ys = max(abs(v) for v in y)
        yi = [_toint(v / ys) for v in y]
        idot = self.idot
        t = [idot(q, yi) for q in self.Qv[:r]]
        res = list(yi)
        for k in range(r):
            res = [a - ((t[k] * b) >> FP) for a, b in zip(res, self.Qv[k])]
        for k in range(r):
            c = idot(self.Qv[k], res)
            res = [a - ((c * b) >> FP) for a, b in zip(res, self.Qv[k])]
            t[k] += c
        ny = mp.sqrt(_fromint(idot(yi, yi)))
        resid = mp.sqrt(_fromint(idot(res, res))) / ny
        t = [_fromint(v) * ys for v in t]
        x = [mp.mpf(0)] * self.n
        for k in range(r - 1, -1, -1):
            s = t[k]
            for kk in range(k + 1, r):
                s -= self.R[k][self.perm[kk]] * x[self.perm[kk]]
            x[self.perm[k]] = s / self.R[k][self.perm[k]]
        return [xi / sc for xi, sc in zip(x, self.scale)], resid, r

    def nullspace(self, tol=mp.mpf(10) ** -40):
        """basis (list of n-vectors, ORIGINAL scaling) of the null space of A"""
        r = self.rank(tol)
        basic = self.perm[:r]
        bset = set(basic)
        nonbasic = [j for j in range(self.n) if j not in bset]
        out = []
        for j in nonbasic:
            v = [mp.mpf(0)] * self.n
            v[j] = mp.mpf(1)
            for k in range(r - 1, -1, -1):
                row = self.R[k]
                s = row.get(j, 0)
                for kk in range(k + 1, r):
                    s += row[self.perm[kk]] * v[self.perm[kk]]
                v[self.perm[k]] = -s / row[self.perm[k]]
            out.append([vi / sc for vi, sc in zip(v, self.scale)])
        return out


from math import isqrt as math_isqrt


def min_norm(x, Nb):
    """minimal-norm representative x - N (N^T N)^{-1} N^T x for null basis Nb (list of vectors)"""
    if not Nb:
        return list(x)
    G = mp.matrix(len(Nb), len(Nb))
    for a in range(len(Nb)):
        for b in range(a, len(Nb)):
            G[a, b] = G[b, a] = dot(Nb[a], Nb[b])
    rhs = mp.matrix([dot(v, x) for v in Nb])
    c = mp.lu_solve(G, rhs)
    out = list(x)
    for a, v in enumerate(Nb):
        out = [o - c[a] * w for o, w in zip(out, v)]
    return out


# ------------------------------------------------------------------------------------------------ F0
def F0(tag='P'):
    hdr('F0: setup and sanity (point %s)' % tag)
    ch = chain()
    pt = get_point(tag)
    zeta = pt['zeta']
    print('zeta =', [str(x) for x in zeta], ' theta =', [str(t) for t in THETA], ' digits', pt['dps'], ' HW index', ch.hw_index)
    print('states', len(pt['states']), ' max ODE residual of all 6 Q:', nst(max(s['resid'] for s in pt['states'])))
    # omega2 -> infinity check
    for w1 in (Fr(5, 2), Fr(7, 4)):
        big = Fr(10) ** 6
        v1 = vacuum_covector_1site(ch, zeta, (w1, big, Fr(3, 5)))
        v2 = onesite_inf(ch, zeta, w1)
        # projective: compare after normalising by the HW-component
        k = ch.hw_index
        d = max(abs(mpf(a / v1[k]) - mpf(b / v2[k])) for a, b in zip(v1, v2))
        print('  Omega^inf[zeta,%s] vs vacuum_covector_1site at omega=(%s,1e6,3/5): max projective diff %s (expected O(1e-6))' % (w1, w1, nst(d)))
    # x0 check
    x = x0_onesite(ch, zeta)
    NB2G = [7, 5, 6, 4, 2, 1, 3, 0]      # notebook index a -> gl3_reps index (sorted(patterns, reverse=True) permutation of xxx_chain)
    nbk = [mpf(a) for a in (Fr(1, 10), Fr(1, 3), Fr(1, 10), Fr(17, 30), Fr(3, 10), Fr(1), Fr(7, 15), Fr(7, 10))]
    dirx = [x[NB2G[a]] / x[NB2G[5]] for a in range(8)]
    print('  x0^G one-site direction mapped to notebook order (component 5 = 1):', [nst(a, 8) for a in dirx])
    print('  RESULTS-xxx-limit X1c direction (notebook order)                  :', [nst(a, 8) for a in nbk])
    if tag == 'P':
        print('  max diff:', nst(max(abs(a - b) for a, b in zip(dirx, nbk))))
    x0 = [a * b for a in x for b in x]
    # overlaps
    base = cov_inf(ch, zeta, *W1BASE)
    covb = cov_to_mp(base)
    rows = []
    for k, st in enumerate(pt['states']):
        ox = overlap(x0, st['vec'])
        op = overlap(covb, st['vec'])
        rows.append((ox, op))
    print('  <x0^G|Psi> for all 64 states (max-component normalised vectors): min |.| = %s, max |.| = %s' % (nst(min(abs(a) for a, b in rows)), nst(max(abs(a) for a, b in rows))))
    print('  <Phi_inf|Psi> at (w1,w1p) = (5/2,7/4): min |.| = %s, max |.| = %s; vanishing states (|.|<1e-30):' % (nst(min(abs(b) for a, b in rows)), nst(max(abs(b) for a, b in rows))),
          [(k, pt['states'][k]['sector']) for k, (a, b) in enumerate(rows) if abs(b) < mp.mpf(10) ** -30])
    # T / Qnum
    dets, vcov = detN_for_states(pt)
    mp.mp.dps = pt['dps']
    vc = cov_to_mp(vcov)
    Rv = [overlap(vc, st['vec']) for st in pt['states']]
    Qn = [b / a for a, b in rows]
    T = [b * d / r for (a, b), d, r in zip(rows, dets, Rv)]
    ratio = [t / q for t, q in zip(T, Qn)]
    r0 = ratio[0]
    print('  T(B)/Qnum(B): max_B |ratio/ratio_0 - 1| = %s   (ratio_0 = %s)' % (nst(max(abs(r / r0 - 1) for r in ratio)), nst(r0, 20)))
    print('  |D_B(omega0)| min/max = %s / %s ; |<Omega_D Omega_D|Psi>| min = %s' % (nst(min(abs(d) for d in dets)), nst(max(abs(d) for d in dets)), nst(min(abs(r) for r in Rv))))
    print('  the 64 Qnum(B) = <Phi_inf|Psi>/<x0^G|Psi> at P, (w1,w1p) = (5/2,7/4); labels (S2(3/7), S3(-2/9)), 30 digits:')
    for k, st in enumerate(pt['states']):
        l2, l3 = S2S3_labels(st)
        print('   %2d %s  S2=%s S3=%s  Qnum=%s' % (k, st['sector'], nst(l2, 12), nst(l3, 12), nst(Qn[k], 30)))
    return pt, Qn, T


# ------------------------------------------------------------------------------------------------ F1
def sv(M):
    S = mp.svd_r(M, compute_uv=False)
    return [S[i] for i in range(len(S))]


def kappa_of(st, s):
    """Frobenius constant kappa_s = a_s/3, a_s = coefficient of (u-theta_s)^-2 in S3 (derivation in F1)"""
    return (st['lam3'][2] if s == 0 else st['lam3'][5]) / 3


def F1(tag='P'):
    hdr('F1: jet data at theta_s, Wronskian facts, Frobenius constraint (point %s)' % tag)
    pt = get_point(tag)
    zeta = pt['zeta']
    th = THETA
    z = [mpf(x) for x in zeta]
    t = [mpf(x) for x in th]
    SJ = [StateJets(st, zeta) for st in pt['states']]
    dev_p01 = dev_p02 = dev_kap = mp.mpf(0)
    kap = []
    sv_rank4, sv_rank5 = [], []
    for k, (st, sj) in enumerate(zip(pt['states'], SJ)):
        for s in range(2):
            J = sj.P[s]
            pl = plucker(J)
            sc = abs(J[1][0] * J[2][1]) + abs(J[2][0] * J[1][1])
            dev_p01 = max(dev_p01, abs(pl[(0, 1)]) / sc)
            Pip = t[s] - t[1 - s]
            pred = (z[2] - z[1]) * mp.exp(sum(z) * t[s]) * Pip * sj.Q[s][0][0]
            dev_p02 = max(dev_p02, abs(pl[(0, 2)] - pred) / abs(pred))
            kappa = kappa_of(st, s)
            for j in range(3):
                d = abs(J[j][1] - kappa * J[j][0]) / (abs(J[j][1]) + abs(kappa * J[j][0]))
                dev_kap = max(dev_kap, d)
            kap.append((k, s, kappa))
            M4 = mp.matrix([[J[j][c] for c in range(4)] for j in range(3)])
            M5 = mp.matrix([[J[j][c] for c in range(5)] for j in range(3)])
            for M in (M4, M5):
                for c in range(M.cols):
                    nrm = mp.sqrt(sum(M[r, c] ** 2 for r in range(3))) or 1
                    for r in range(3):
                        M[r, c] /= nrm
            s4, s5 = sv(M4), sv(M5)
            sv_rank4.append(s4[2] / s4[0])
            sv_rank5.append(s5[2] / s5[0])
    print("  (a) Wronskian: max |p01|/(|J2_0 J3_1|+|J3_0 J2_1|) = %s (p01 = W(P2,P3)(theta_s)); max rel dev of p02 from (z3-z2) e^{sum z theta_s} Pi'(theta_s) Q1(theta_s) = %s" % (nst(dev_p01), nst(dev_p02)))
    print("  (b) Frobenius P^j'(theta_s) = kappa_s P^j(theta_s), kappa_s := a_s/3 (a_s = coefficient of (u-theta_s)^-2 in S3): max over states, sites, j of |J1 - kappa J0|/(|J1| + |kappa J0|) = %s" % nst(dev_kap))
    print('      derivation: indicial polynomial rho(rho-1)(rho-2) - 3 rho(rho-1) + 3 rho = rho(rho-2)(rho-4). Acting on c0 + c1 x + ... (x = u-theta_s), the x^-2 coefficient of the oper')
    print('      is I(1) c1 - a_s c0 with I(1) = 1(-1)(-3) = 3: only the leading singular parts 3/x (S1), 3/x^2 (S2) and the double pole a_s/x^2 of S3 enter;')
    print('      H_s, b_s and the regular parts of S1, S2 first appear at order x^-1 (c2). => kappa_s = c1/c0 = a_s/3.')
    print('      kappa_s (site 1 | site 2), states 0..7:', [nst(kk, 8) for (k, s, kk) in kap[0:16:2]], '|', [nst(kk, 8) for (k, s, kk) in kap[1:16:2]])
    print('  (c) jet matrix [J^j_k] j=1..3 (columns normalised): k=0..3: max over states, sites of sigma_3/sigma_1 = %s (rank 2); k=0..4: min sigma_3/sigma_1 = %s (rank 3)' % (nst(max(sv_rank4)), nst(min(sv_rank5))))
    print('      structure: exponent-4 solution has J_0..J_3 = 0, exponent-2 solution has J_0 = J_1 = 0; all three have J_1 = kappa_s J_0.')
    return pt, SJ


# ------------------------------------------------------------------------------------------------ fit data
def targets(pt, w1, w1p, th=THETA):
    ch = chain(th)
    zeta = pt['zeta']
    cov = cov_to_mp(cov_inf(ch, zeta, w1, w1p))
    x = x0_onesite(ch, zeta)
    x0 = [a * b for a in x for b in x]
    return [(overlap(cov, st['vec']), overlap(x0, st['vec'])) for st in pt['states']]


def detdata(pt, th=THETA):
    path = os.path.join(HERE, 'twovac_fit_detN_%s_%d%s.pkl' % (''.join(str(x) for x in pt['zeta']).replace('/', '_'), pt['dps'], thtag(th)))
    if os.path.exists(path):
        with open(path, 'rb') as f:
            return pickle.load(f)
    dets, vcov = detN_for_states(pt, th=th)
    mp.mp.dps = pt['dps']
    vc = cov_to_mp(vcov)
    Rv = [overlap(vc, st['vec']) for st in pt['states']]
    out = dict(dets=dets, Rv=Rv)
    with open(path, 'wb') as f:
        pickle.dump(out, f)
    return out


def kappaG(zeta, w1):
    """frame-tie scalar kappa^G(zeta, omega) for Lambda=(2,1,0) at omega_2 -> infinity (W2): Delta_1 Delta_2"""
    z1, z2, z3 = [Fr(x) for x in zeta]
    w1 = Fr(w1)
    D1 = w1 ** 2 / ((z2 - w1) * (z3 - w1))
    D2 = w1 / ((z1 - z2) * (z3 - w1))
    return mpf(D1 * D2)


_G2 = None
NB2G = [7, 5, 6, 4, 2, 1, 3, 0]      # notebook index a -> gl3_reps index


def x0_D2(ch, zeta):
    """D2/L1(b): x^G(zeta) = (z2-z1, z3-z1, z2-z1, 2z3-z1-z2, 3(z2-z1), 3(z3-z1), 2(z3-z2), 3(z3-z2)) in notebook order -> gl3_reps order, x (x) x"""
    z1, z2, z3 = [Fr(x) for x in zeta]
    nb = [z2 - z1, z3 - z1, z2 - z1, 2 * z3 - z1 - z2, 3 * (z2 - z1), 3 * (z3 - z1), 2 * (z3 - z2), 3 * (z3 - z2)]
    x = [None] * 8
    for a in range(8):
        x[NB2G[a]] = mpf(nb[a])
    return [a * b for a in x for b in x]


def cscal(zeta, w1, w1p):
    z1, z2, z3 = [Fr(x) for x in zeta]
    f = lambda w: (z2 - w) * (z3 - w) ** 2
    return mpf(f(Fr(w1)) * f(Fr(w1p)) * (z1 - z2) ** 2)


class Data:
    """rows over (zeta point, state) at fixed (w1, w1p, theta): jets, Qnum, T"""
    def __init__(self, tags, w1=W1BASE[0], w1p=W1BASE[1], th=THETA, with_T=True, dps=None):
        dps = dps or DPS
        mp.mp.dps = dps
        self.th = th
        self.rows = []
        x0Dx = {}
        for tag in tags:
            pt = get_point(tag, dps=dps, th=th)
            mp.mp.dps = dps
            tg = targets(pt, w1, w1p, th)
            x0Dx[tag] = x0_D2(chain(th), pt['zeta'])
            cs = cscal(pt['zeta'], w1, w1p)
            dd = detdata(pt, th) if with_T else None
            for k, st in enumerate(pt['states']):
                sj = StateJets(st, pt['zeta'], th)
                op, ox = tg[k]
                r_x0 = x0Dx[tag]
                kk = kappaG(pt['zeta'], w1) * kappaG(pt['zeta'], w1p)
                r = dict(tag=tag, k=k, sj=sj, st=st, Qnum=op / ox, Phi=op, X0=ox, Qc=op / ox / kk, kap=kk)
                r['Qd2'] = op / overlap(r_x0, st['vec'])
                r['Y0'] = 4 * cs * r['Qd2']
                r['lau'] = [st['lam2'][1], st['lam2'][3], st['lam3'][2], st['lam3'][5], st['lam3'][1], st['lam3'][4]]
                if with_T:
                    r['T'] = op * dd['dets'][k] / dd['Rv'][k]
                    r['Tc'] = r['T'] / kk
                    r['Deq'] = dd['dets'][k]
                self.rows.append(r)
        self.n = len(self.rows)


# per-site basis functions: take StateJets, s (site index) and return a list of numbers
def b_dir(sj, s, ks=(0, 1, 2), cs=(0, 1), lone=2):
    P = sj.P[s]
    return [sj.Q[s][0][k] * P[lone][c] for k in ks for c in cs]


def b_pl(sj, s, order=3, cs=(0, 1), lone=2, pair=(1, 2)):
    P = sj.P[s]
    pl = plucker(P, pair[0], pair[1], order)
    return [pl[kl] * P[lone][c] for kl in sorted(pl) for c in cs]


def b_mono(sj, s, order=3, a_set=None, lone=(1, 2), sym=True):
    """all J^2_a J^3_b J^3_c (a, b <= c <= order): trilinear monomials with one P^2 jet and two P^3 jets (what det N contains)"""
    P = sj.P[s]
    out = []
    for a in range(order + 1):
        for b in range(order + 1):
            for c in range(b, order + 1):
                out.append(P[1][a] * P[2][b] * P[2][c])
    return out


def mono_names(order=3):
    return [(a, b, c) for a in range(order + 1) for b in range(order + 1) for c in range(b, order + 1)]


def pw(sj, s):
    """P3 p02 at theta_s (the vacuum monomial of D2) and the five D2 monomials N = (P3 p02, P3' p02, P3 p03, P3' p03, P3 p23)"""
    P = sj.P[s]
    pl = plucker(P, 1, 2, 3)
    return P, pl


def b_N(sj, s):
    P, pl = pw(sj, s)
    return [P[2][0] * pl[(0, 2)], P[2][1] * pl[(0, 2)], P[2][0] * pl[(0, 3)], P[2][1] * pl[(0, 3)], P[2][0] * pl[(2, 3)]]


def d2_G(w1, w1p, thi):
    """D2's untwisted-jet kernel G_ij(w1, w1p, 1/theta12) parsed from RESULTS-twovac-limit.md (exact Fractions)"""
    global _G2
    if '_G2' not in globals() or _G2 is None:
        import re
        txt = open(os.path.join(HERE, 'RESULTS-twovac-limit.md'), encoding='utf8').read()
        blk = txt[txt.index('G[1,1] = -4 w1^2'):]
        blk = blk[:blk.index('```')]
        w, wp, t = sympy.symbols('w1 w1p thi')
        _G2 = {}
        for line in blk.splitlines():
            m = re.match(r'G\[(\d),(\d)\] = (.*)', line.strip())
            if m:
                _G2[(int(m.group(1)) - 1, int(m.group(2)) - 1)] = parse_expr(m.group(3).replace('^', '**'), local_dict=dict(w1=w, w1p=wp, thi=t), transformations=standard_transformations + (implicit_multiplication,))
    w, wp, t = sympy.symbols('w1 w1p thi')
    return [[sympy.Rational(_G2[(i, j)].subs({w: sympy.Rational(w1.numerator, w1.denominator), wp: sympy.Rational(w1p.numerator, w1p.denominator),
                                          t: sympy.Rational(thi.numerator, thi.denominator)})) for j in range(5)] for i in range(5)]


def AS(sj, s):
    z, t = sj.z, sj.t
    return (z[2] - z[1]) * mp.exp(sum(z) * t[s]) * (t[s] - t[1 - s])


def b_dirA(sj, s, ks=(0, 1), cs=(0, 1, 2, 3), lone=2):
    """A_s Q1^(k) P3^(c): the direct-side basis in the normalisation of D2 (A_s Q1 = p02)"""
    A = AS(sj, s)
    return [A * f for f in b_dir(sj, s, ks, cs, lone)]


def Mfac(sj, s):
    return sj.Q[s][0][0] * sj.P[s][2][0]


def reduce_basis(data, bas, s, tol=mp.mpf(10) ** -40):
    """indices of a maximal independent subset of the per-site features of `bas` at site s over all data rows"""
    F = [_clean(bas(r['sj'], s)) for r in data.rows]
    n = len(F[0])
    cols = []
    for j in range(n):
        c = [F[i][j] for i in range(len(F))]
        m = max(abs(x) for x in c)
        cols.append([x / m for x in c] if m != 0 else c)
    qr = PivQR(cols)
    r = qr.rank(tol)
    return sorted(qr.perm[:r]), qr


def sel_basis(bas, sel):
    return lambda sj, s: [bas(sj, s)[i] for i in sel]


def _clean(f):
    """set features that vanish identically (e.g. p01 = W(theta_s) = 0, 1e-78 noise) to exact zero"""
    mx = max(abs(a) for a in f)
    return [a if abs(a) > mx * mp.mpf(10) ** -45 else mp.mpf(0) for a in f]


def design(data, bas1, bas2, target='Qnum', Mf=True, rowscale=True, extra=None):
    """columns: all products b1_a b2_b (a-major); rows over data rows. Returns (cols, y)"""
    A = []
    y = []
    for r in data.rows:
        f1 = _clean(bas1(r['sj'], 0))
        f2 = _clean(bas2(r['sj'], 1))
        row = [a * b for a in f1 for b in f2]
        if extra is not None:                       # kernel linear in Laurent data: multiply by (1, Delta_m)
            row = [c * e for e in extra(r) for c in row]
        tv = r[target]
        if Mf == 'D2':
            tv = tv * pw(r['sj'], 0)[0][2][0] * pw(r['sj'], 0)[1][(0, 2)] * pw(r['sj'], 1)[0][2][0] * pw(r['sj'], 1)[1][(0, 2)]
        elif Mf:
            tv = tv * Mfac(r['sj'], 0) * Mfac(r['sj'], 1)
        if rowscale:
            nr = mp.sqrt(dot(row, row))
            row = [c / nr for c in row]
            tv = tv / nr
        A.append(row)
        y.append(tv)
    cols = [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]
    return cols, y


def build_all(tags, th=THETA, dps=None, workers=10, dets=None):
    """build points (+ det N data on the base theta only) for many zeta tags in parallel processes"""
    from concurrent.futures import ProcessPoolExecutor
    dets = (tuple(Fr(t) for t in th) == THETA) if dets is None else dets
    with ProcessPoolExecutor(max_workers=workers) as ex:
        list(ex.map(_build_one, [(t, zeta_of(t), th, dps or DPS, dets) for t in tags]))


def _build_one(a):
    tag, zeta, th, dps, dets = a
    mp.mp.dps = dps
    pt = get_point(tag, zeta=zeta, dps=dps, th=th)
    if dets:
        detdata(pt, th)
    return tag


# ------------------------------------------------------------------------------------------------ F0b: normalisations
def F0b(tags=('P', 'Pp', 'X1', 'X2', 'X3')):
    hdr('F0b: x0 normalisations, per-zeta scalars, T versus Qnum')
    D = Data(list(tags))
    print('  frame-tie scalar of W2: kappa^G(zeta,omega) = Delta_1 Delta_2 = w^3/((z2-w)(z3-w)^2(z1-z2)); c(zeta,w,w\') := (z2-w)(z3-w)^2(z2-w\')(z3-w\')^2(z1-z2)^2 = w^3 w\'^3/(kappa kappa\')')
    print('  Qnum_D2 := <Phi_inf|Psi>/<x_D2 (x) x_D2|Psi> with x_D2 = (z2-z1, z3-z1, z2-z1, 2z3-z1-z2, 3(z2-z1), 3(z3-z1), 2(z3-z2), 3(z3-z2)) (notebook order, mapped to gl3_reps order)')
    for tag in tags:
        rr = [r for r in D.rows if r['tag'] == tag]
        s1 = [r['Qnum'] / r['Qd2'] for r in rr]          # x_D2-overlap / my-x0-overlap ratio (state independent?)
        sp1 = max(abs(a / s1[0] - 1) for a in s1)
        cT = [r['T'] / (r['Qnum'] * Mfac(r['sj'], 0) * Mfac(r['sj'], 1)) for r in rr]
        cT2 = [r['T'] / (r['Qd2'] * pw(r['sj'], 0)[0][2][0] * pw(r['sj'], 0)[1][(0, 2)] * pw(r['sj'], 1)[0][2][0] * pw(r['sj'], 1)[1][(0, 2)]) for r in rr]
        spT = max(abs(a / cT[0] - 1) for a in cT)
        spT2 = max(abs(a / cT2[0] - 1) for a in cT2)
        spTQ = max(abs((r['T'] / r['Qnum']) / (rr[0]['T'] / rr[0]['Qnum']) - 1) for r in rr)
        print('  %-3s my-x0 vs x_D2 overlap ratio: %s (state spread %s) | T/(Qnum M1 M2), M=Q1 P3: %s (spread %s) | T/(Qnum_D2 (P3p02)_1 (P3p02)_2): %s (spread %s) | T/Qnum state spread %s' % (
            tag, nst(s1[0], 10), nst(sp1), nst(cT[0], 10), nst(spT), nst(cT2[0], 10), nst(spT2), nst(spTQ)))
    print('  => T(B) = c_T(zeta) Qnum(B) M1 M2 with a state-independent but zeta-DEPENDENT c_T; T/Qnum is NOT state independent (<x0|Psi>/n_B = const M1 M2, exponents 1 on Q1(th_s), P3(th_s), 0 on the rest: log-regression residual 2e-77).')
    return D


# ------------------------------------------------------------------------------------------------ F2
def rat_recon(x, maxden=10 ** 12, tol=mp.mpf(10) ** -40):
    """Fraction approximating mp number x, or None"""
    f = Fr(int(mp.nint(x * 10 ** 40)), 10 ** 40).limit_denominator(maxden)
    if abs(mpf(f) - x) < tol * max(1, abs(x)):
        return f
    return None


def fit_level(D, bas1, bas2, target, Mf, label, extra=None, reduce=True):
    s0 = list(range(len(bas1(D.rows[0]['sj'], 0))))
    s1 = list(range(len(bas2(D.rows[0]['sj'], 1))))
    raw = len(s0) * len(s1) * (len(extra(D.rows[0])) if extra is not None else 1)
    if reduce:
        s0, _ = reduce_basis(D, bas1, 0)
        s1, _ = reduce_basis(D, bas2, 1)
    b0 = sel_basis(bas1, s0)
    b1 = sel_basis(bas2, s1)
    cols, y = design(D, b0, b1, target, Mf, extra=extra)
    qr = PivQR(cols)
    x, res, rk = qr.lstsq(y)
    return dict(label=label, unknowns=raw if reduce else len(cols), site_rank=(len(s0), len(s1)), rank=rk, null=(raw if reduce else len(cols)) - rk, resid=res, x=x, qr=qr, s0=s0, s1=s1, ncols=len(cols))


def b_dirK(sj, s):
    """D2's nine direct-side monomials A_s x (P3 Q1, k P3 Q1, k^2 P3 Q1, P3 Q1', k P3 Q1', P3'' Q1, k P3'' Q1, P3'' Q1', P3''' Q1), k = P3'/P3"""
    P = sj.P[s][2]
    Q1 = sj.Q[s][0]
    A = AS(sj, s)
    k = P[1] / P[0]
    return [A * v for v in (P[0] * Q1[0], k * P[0] * Q1[0], k * k * P[0] * Q1[0], P[0] * Q1[1], k * P[0] * Q1[1], P[2] * Q1[0], k * P[2] * Q1[0], P[2] * Q1[1], P[3] * Q1[0])]


def b_Nlone(lone):
    def f(sj, s):
        P = sj.P[s]
        pl = plucker(P, 1, 2, 3)
        return [P[lone][0] * pl[(0, 2)], P[lone][1] * pl[(0, 2)], P[lone][0] * pl[(0, 3)], P[lone][1] * pl[(0, 3)], P[lone][0] * pl[(2, 3)]]
    return f


def b_zero_jets(Dzero):
    """basis evaluated from precomputed jets at u = 0 (control); Dzero[(tag,k)] -> StateJets at theta=(0,0)"""
    return None


def F2(tags=('P', 'Pp', 'X1', 'X2', 'X3'), tags_big=('P', 'Pp', 'X1', 'X2', 'X3', 'X4', 'X5', 'X6')):
    hdr('F2: the fit (base point (w1,w1p) = (5/2,7/4), theta = (2/7,-3/5)); target Y0 = 4 c Qnum_D2, with factor M1 M2 (M = P3 p02) or without')
    D5 = Data(list(tags))
    D8 = Data(list(tags_big))
    print('  rows: %d (5 zeta points), %d (8 zeta points); kernel ZETA-INDEPENDENT after the c(zeta,w,w\') scaling (stacked over all points)' % (D5.n, D8.n))
    ch = chain()
    levels = [
        ('a  N5 = {P3^(c) p_kl}: (02,c0),(02,c1),(03,c0),(03,c1),(23,c0)  [D2 basis]', b_N, D5),
        ('a  A_s Q1^(k) P3^(c), k<=2, c<=1  (direct basis 6)', lambda sj, s: b_dirA(sj, s, (0, 1, 2), (0, 1)), D5),
        ('a  A_s Q1^(k) P3^(c), k<=1, c<=3  (8)', lambda sj, s: b_dirA(sj, s, (0, 1), (0, 1, 2, 3)), D8),
        ('a  D2 direct monomials with kappa = P3\'/P3 (9)', b_dirK, D5),
        ('b  Plucker basis {p_kl P3^(c)}, k<l<=3, c<=1 (12)', b_pl, D5),
        ('c  Plucker k<l<=3, c<=3 (24)', lambda sj, s: b_pl(sj, s, 3, (0, 1, 2, 3)), D8),
        ('c  A_s Q1^(k) P3^(c), k<=3, c<=3 (16)', lambda sj, s: b_dirA(sj, s, (0, 1, 2, 3), (0, 1, 2, 3)), D8),
    ]
    print('  %-70s %-9s %-8s %-6s %-5s %-5s %-12s %-12s' % ('level (basis per site; same on both sites)', 'unknowns', 'sitrank', 'rank', 'null', 'rows', 'resid (M)', 'resid (no M)'))
    out = {}
    for name, b, D in levels:
        rM = fit_level(D, b, b, 'Y0', 'D2', name)
        rN = fit_level(D, b, b, 'Y0', False, name)
        out[name] = (rM, rN)
        print('  %-70s %-9d %-8s %-6d %-5d %-5d %-12s %-12s' % (name[:70], rM['unknowns'], '%d/%d' % rM['site_rank'], rM['rank'], rM['null'], D.n, nst(rM['resid']), nst(rN['resid'])), flush=True)
    print('  (T(B) as target: T = c_T(zeta) Qnum M1 M2 with a zeta-dependent constant, so T at stacked zeta points equals the M-column above only up to one free scalar per zeta point; not refit)')
    # level d: kernel linear in the six Laurent data
    ext = lambda r: [mp.mpf(1)] + [v / 100 for v in r['lau']]
    rd = fit_level(D8, b_N, b_N, 'Y0', 'D2', 'd', extra=ext)
    base = out[levels[0][0]][0]
    print('  d  N5 with K = K0 + sum_m Delta_m K_m (Laurent data H1,H2,a1,a2,b1,b2): unknowns %d rank %d null %d resid %s' % (rd['unknowns'], rd['rank'], rd['null'], nst(rd['resid'])))
    Nb = rd['qr'].nullspace()
    xm = min_norm(rd['x'], Nb)
    nK0 = mp.sqrt(sum(v ** 2 for v in xm[-25 * 0:][:25])) if False else None
    # in design(extra) the column order is: extra index major (7 groups of 25)
    grp = [mp.sqrt(sum(v ** 2 for v in xm[25 * g:25 * (g + 1)])) for g in range(7)]
    print('     minimal-norm representative: |K0| = %s, |K_m| (m=1..6) = %s  (nonzero only because the Laurent data are linearly tied to the jets, kappa_s = a_s/3 and H1 + H2 = const - W: the null space is 56-dim; the K0-only fit has residual 0 at rank 25, so no Laurent-dependent kernel is needed)' % (nst(grp[0]), [nst(g) for g in grp[1:]]))
    # kernel at the base point
    G = d2_G(W1BASE[0], W1BASE[1], Fr(1) / (THETA[0] - THETA[1]))
    x = base['x']
    print('  Kernel at the base point in basis N5 (fit over 5 zeta points, 320 rows; rank %d, null %d => unique):' % (base['rank'], base['null']))
    mx = mp.mpf(0)
    nrat = 0
    print('     i\\j     ' + ' '.join('%-24s' % j for j in range(1, 6)))
    for i in range(5):
        row = []
        for j in range(5):
            v = x[5 * i + j]
            g = mpf(Fr(int(G[i][j].p), int(G[i][j].q)))
            mx = max(mx, abs(v - g))
            f = rat_recon(v)
            nrat += f is not None
            row.append('%-24s' % (str(f) if f is not None else nst(v, 15)))
        print('     %d       ' % (i + 1) + ' '.join(row))
    print('     all 25 entries rational (denominators <= 1e12, 50-digit test): %d/25; max |fit - D2 G| = %s' % (nrat, nst(mx)))
    # Plucker 12x12 level: minimal norm and sparse representative
    rpl = fit_level(D5, b_pl, b_pl, 'Y0', 'D2', 'pl12 raw', reduce=False)
    Nb = rpl['qr'].nullspace()
    xm = min_norm(rpl['x'], Nb)
    print('  12x12 Plucker level: rank %d, null %d (of %d), residual %s' % (rpl['rank'], rpl['null'], rpl['unknowns'], nst(rpl['resid'])))
    print('     minimal-norm representative: |K|_F = %s, support (|entry|>1e-30) %d of %d, rational entries %d' % (
        nst(mp.sqrt(sum(v ** 2 for v in xm))), sum(1 for v in xm if abs(v) > 1e-30), len(xm), sum(1 for v in xm if rat_recon(v) is not None)))
    print('     sparse representative: the N5 kernel above (25 nonzero entries, support in {P3p02, P3\'p02, P3p03, P3\'p03, P3p23}^2, = 5 of the 12 Plucker-basis elements per site); it is minimal since the 25 products are independent (rank 25).')
    # embedding check: N5 inside the 12x12 basis, residual of the pl12 fit evaluated at the sparse kernel
    # controls
    print('  Controls (basis N5 shape, same target Y0 with M; must FAIL):')
    c1 = fit_level(D5, b_Nlone(1), b_Nlone(1), 'Y0', 'D2', 'lone P2')
    print('     lone index P2 instead of P3: rank %d resid %s' % (c1['rank'], nst(c1['resid'])))
    c1b = fit_level(D5, b_Nlone(0), b_Nlone(0), 'Y0', 'D2', 'lone P1')
    print('     lone index P1 instead of P3: rank %d resid %s' % (c1b['rank'], nst(c1b['resid'])))
    c2 = fit_level(D5, lambda sj, s: b_N(sj, 0), lambda sj, s: b_N(sj, 0), 'Y0', 'D2', 'both at theta1')
    print('     both factors from the site-1 jets (swapped/collapsed sites): rank %d resid %s' % (c2['rank'], nst(c2['resid'])))
    c2b = fit_level(D5, lambda sj, s: b_N(sj, 1), lambda sj, s: b_N(sj, 1), 'Y0', 'D2', 'both at theta2')
    print('     both factors from the site-2 jets: rank %d resid %s' % (c2b['rank'], nst(c2b['resid'])))
    # jets at u = 0
    for r in D5.rows:
        r['sj0'] = StateJets(r['st'], r['sj'].z, (Fr(0), Fr(1, 3)))
    # at wrong points: replace sj by jets at (0, 1/3) for the basis only
    class _W:
        pass
    Dw = Data.__new__(Data)
    Dw.rows = []
    for r in D5.rows:
        r2 = dict(r)
        r2['sj'] = r['sj0']
        Dw.rows.append(r2)
    Dw.n = D5.n
    # target needs the true M factor: precompute Y0*M_true
    for r, r2 in zip(D5.rows, Dw.rows):
        r2['Yw'] = r['Y0'] * pw(r['sj'], 0)[0][2][0] * pw(r['sj'], 0)[1][(0, 2)] * pw(r['sj'], 1)[0][2][0] * pw(r['sj'], 1)[1][(0, 2)]
    c3 = fit_level(Dw, b_N, b_N, 'Yw', False, 'u=0')
    print('     jets taken at the wrong points (u = 0, 1/3) with the true target: rank %d resid %s' % (c3['rank'], nst(c3['resid'])))
    c4 = fit_level(D5, b_N, b_N, 'Y0', False, 'noM')
    print('     right basis but target without the M1 M2 factor: rank %d resid %s' % (c4['rank'], nst(c4['resid'])))
    return D5, base


# ------------------------------------------------------------------------------------------------ F3
W_VALUES = [Fr(5, 2), Fr(-2, 5), Fr(9, 7), Fr(3), Fr(-5, 3), Fr(7, 4), Fr(11, 3), Fr(-7, 2), Fr(13, 5)]
THETA_VARS = [(Fr(1, 5), Fr(-4, 9)), (Fr(3, 7), Fr(-1, 4)), (Fr(1, 2), Fr(-1, 3)), (Fr(-1, 7), Fr(-6, 5)), (Fr(4, 5), Fr(-2, 3))]
Z5 = ('P', 'Pp', 'X1', 'X2', 'X3')


def fitG(tags, w1, w1p, th=THETA):
    """fit the 5x5 kernel in the N5 basis (25 unknowns) on the given zeta tags; returns (G as 5x5 mp, rank, resid)"""
    D = Data(list(tags), w1=w1, w1p=w1p, th=th, with_T=False)
    cols, y = design(D, b_N, b_N, 'Y0', 'D2', rowscale=False)
    qr = PivQR(cols)
    x, res, rk = qr.lstsq(y)
    return [[x[5 * i + j] for j in range(5)] for i in range(5)], rk, res


def G_vs_D2(G, w1, w1p, th):
    G2 = d2_G(w1, w1p, Fr(1) / (Fr(th[0]) - Fr(th[1])))
    return max(abs(G[i][j] - mpf(Fr(int(G2[i][j].p), int(G2[i][j].q)))) / max(1, abs(mpf(Fr(int(G2[i][j].p), int(G2[i][j].q))))) for i in range(5) for j in range(5))


def F3():
    hdr('F3: kernel laws over the (w1,w1p) grid, theta variants, scaling symmetry, per-zeta independence')
    tags = Z5
    # (0) zeta independence: per-zeta fits (rows = 64 per point, 25 unknowns) compared with the stacked kernel
    G0, rk0, res0 = fitG(tags, *W1BASE)
    print('  (0) stacked fit over 5 zeta points: rank %d, residual %s' % (rk0, nst(res0)))
    print('      per-zeta-point fits (64 rows, 25 unknowns, rank <= 25) versus the stacked kernel:')
    for t in tags:
        Gt, rk, res = fitG([t], *W1BASE)
        dev = max(abs(Gt[i][j] - G0[i][j]) / max(1, abs(G0[i][j])) for i in range(5) for j in range(5))
        print('        zeta %-3s rank %d residual %s  max |G_zeta - G_stacked| = %s' % (t, rk, nst(res), nst(dev)))
    # (1) 3x3 grid
    print('  (1) 3x3 grid w1 in {5/2,-2/5,9/7} x w1p in {7/4,3,-5/3} (theta base), kernel fit vs the analytic D2 kernel:')
    grid = {}
    for w1 in (Fr(5, 2), Fr(-2, 5), Fr(9, 7)):
        for w1p in (Fr(7, 4), Fr(3), Fr(-5, 3)):
            G, rk, res = fitG(tags, w1, w1p)
            grid[(w1, w1p)] = G
            print('        (%s, %s): rank %d residual %s  max |fit - D2 G| = %s  rational entries %d/25' % (w1, w1p, rk, nst(res), nst(G_vs_D2(G, w1, w1p, THETA)),
                  sum(1 for i in range(5) for j in range(5) if rat_recon(G[i][j]) is not None)), flush=True)
    # (2) polynomial dependence on (w1,w1p) at fixed theta, degree <= 5 total, 9x9 minus diagonal
    print('  (2) each kernel entry as a polynomial of total degree <= 6 in (w1,w1p) at theta base: fit on all ordered pairs of 9 values with w1 != w1p (72 points; 28 coefficients)')
    pts = [(a, b) for a in W_VALUES for b in W_VALUES if a != b]
    Gs = {}
    for (a, b) in pts:
        Gs[(a, b)] = fitG(tags, a, b)[0]
    mons = [(i, j) for d in range(7) for i in range(d + 1) for j in [d - i]]
    cols = [[mpf(a) ** i * mpf(b) ** j for (a, b) in pts] for (i, j) in mons]
    qr = PivQR(cols)
    wsym, wpsym, tsym = sympy.symbols('w1 w1p thi')
    th12 = THETA[0] - THETA[1]
    laws = {}
    maxres = mp.mpf(0)
    for i in range(5):
        for j in range(5):
            y = [Gs[p][i][j] for p in pts]
            x, res, rk = qr.lstsq(y)
            maxres = max(maxres, res)
            expr = 0
            ok = True
            for (e1, e2), c in zip(mons, x):
                f = rat_recon(c)
                if f is None and abs(c) > mp.mpf(10) ** -30:
                    ok = False
                if f is not None and f != 0:
                    expr += sympy.Rational(f.numerator, f.denominator) * wsym ** e1 * wpsym ** e2
            laws[(i, j)] = (expr, ok)
    print('      max residual of the 25 polynomial fits: %s; all coefficients rational: %s' % (nst(maxres), all(v[1] for v in laws.values())))
    G2sym = {}
    d2_G(Fr(1), Fr(2), Fr(1))      # populate the parse cache
    nbad = 0
    for i in range(5):
        for j in range(5):
            ex = laws[(i, j)][0]
            ref = _G2[(i, j)].subs(tsym, sympy.Rational(1, 1) / sympy.Rational(th12.numerator, th12.denominator))
            if sympy.expand(ex - ref) != 0:
                nbad += 1
    print('      fitted polynomial laws equal the D2 analytic G_ij(w1,w1p,1/theta12) (exact rational comparison): %d/25 identical' % (25 - nbad))
    print('      no poles: every entry is a polynomial in (w1,w1p) (fit residual above), none in w1 - w1p or w1 - zeta_j (the kernel is zeta-free).')
    # (3) theta dependence
    print('  (3) theta variants at (w1,w1p) = (5/2,7/4): kernel fit (5 zeta points each) vs D2')
    thetas = [THETA] + THETA_VARS
    Gth = {}
    for th in THETA_VARS:
        build_all(list(tags), th=th, workers=10)
        G, rk, res = fitG(tags, *W1BASE, th=th)
        Gth[th] = G
        print('        theta = (%s,%s) (theta12 = %s): rank %d residual %s  max |fit - D2 G| = %s' % (th[0], th[1], th[0] - th[1], rk, nst(res), nst(G_vs_D2(G, W1BASE[0], W1BASE[1], th))), flush=True)
    Gth[THETA] = G0
    # entry-wise polynomial in thi, degree <= 4 (weights: G has degree 10-w_i-w_j in (w, thi))
    ths = [(Fr(1) / (t[0] - t[1])) for t in thetas]
    mons2 = list(range(5))
    qr2 = PivQR([[mpf(x) ** e for x in ths] for e in mons2])
    nbad = 0
    print('      entry-wise polynomials in thi = 1/theta12 of degree <= 4 (6 theta values, 5 coefficients, 1 redundant check):')
    maxr = mp.mpf(0)
    for i in range(5):
        for j in range(5):
            y = [Gth[t][i][j] for t in thetas]
            x, res, rk = qr2.lstsq(y)
            maxr = max(maxr, res)
            ex = sum(sympy.Rational(rat_recon(c).numerator, rat_recon(c).denominator) * tsym ** e for e, c in enumerate(x) if rat_recon(c) not in (None, 0))
            ref = _G2[(i, j)].subs({wsym: sympy.Rational(5, 2), wpsym: sympy.Rational(7, 4)})
            if sympy.expand(ex - ref) != 0:
                nbad += 1
    print('        max residual %s; identical to the D2 polynomial in thi at (5/2,7/4): %d/25' % (nst(maxr), 25 - nbad))
    # (4) scaling symmetry c = 2: (theta, zeta, omega) -> (c theta, zeta/c, omega/c)
    print('  (4) scaling symmetry u -> c u, g -> g/c at c = 2: data at (2 theta, zeta/2, omega/2)')
    th2 = (2 * THETA[0], 2 * THETA[1])
    tagsS = []
    for t in tags:
        nt = 'S2' + t
        ZETAS[nt] = tuple(z / 2 for z in zeta_of(t))
        tagsS.append(nt)
    build_all(tagsS, th=th2, workers=10)
    w1s, w1ps = W1BASE[0] / 2, W1BASE[1] / 2
    G2s, rk, res = fitG(tagsS, w1s, w1ps, th=th2)
    print('        rank %d residual %s  max |fit - D2 G(w/2, 2 theta)| = %s' % (rk, nst(res), nst(G_vs_D2(G2s, w1s, w1ps, th2))))
    print('        entrywise ratio G\'_ij / G_ij (base point) and 2^-(10-w_i-w_j) (D2 homogeneity degree) times the jet-scaling factor:')
    wts = [2, 3, 3, 4, 5]
    rows = []
    for i in range(5):
        rows.append(' '.join('%9s' % nst(G2s[i][j] / G0[i][j], 5) for j in range(5)))
    print('\n'.join('          ' + r for r in rows))
    print('        D2 homogeneity check G(w/2, w\'/2, thi/2) = 2^-(10-w_i-w_j) G(w, w\', thi) exactly: ', end='')
    Ga = d2_G(Fr(5, 2), Fr(7, 4), Fr(35, 31))
    Gb = d2_G(Fr(5, 4), Fr(7, 8), Fr(35, 62))
    print(all(Gb[i][j] == Ga[i][j] / sympy.Integer(2) ** (10 - wts[i] - wts[j]) for i in range(5) for j in range(5)))
    # (5) epsilon-order of each entry: lowest total degree in (w1, w1p) at fixed thi
    print('  (5) (w1,w1p) -> (eps w1, eps w1p): epsilon-order of each kernel entry (lowest total degree of the fitted polynomial in (w1,w1p); thi fixed):')
    for i in range(5):
        row = []
        for j in range(5):
            ex = laws[(i, j)][0]
            pol = sympy.Poly(ex, wsym, wpsym)
            row.append(min(sum(m) for m in pol.monoms()) if not ex == 0 else None)
        print('          ' + ' '.join('%4s' % r for r in row))
    return laws


# ------------------------------------------------------------------------------------------------ F4
PHI = (Fr(1, 3), Fr(-2, 7), Fr(5, 11))
OMH1 = (Fr(3, 4), Fr(-5, 3))
OMH1P = (Fr(-12, 5), Fr(7, 2))
EPSS = (Fr(1, 10 ** 6), Fr(2, 10 ** 6), Fr(4, 10 ** 6))
HW_SECT = {(4, 2, 0): 1, (4, 1, 1): 1, (3, 3, 0): 1, (3, 2, 1): 2, (2, 2, 2): 1}


def build_regime(eps, dps=120):
    """states at zeta = eps*phi (120 digits), the six HW-limit branches (smallest raising residual per HW sector), their dual/direct Q's"""
    cache = os.path.join(HERE, 'twovac_fit_reg_%s_%d.pkl' % (str(eps).replace('/', '_'), dps))
    if os.path.exists(cache):
        with open(cache, 'rb') as f:
            out = pickle.load(f)
        mp.mp.dps = dps
        return out
    zeta = tuple(eps * p for p in PHI)
    Q.set_regime(list(zeta), Q.OM, Q.OMP, dps, 60)
    mp.mp.dps = dps
    ch = chain()
    g = [[zeta[i] if i == j else Fr(0) for j in range(3)] for i in range(3)]
    ops = Operators(ch, g)
    C = ops.fit()
    st = eigenstates(ch, ops, C, dps=dps)
    for x in st:
        x['ev'] = real(x['ev'], mp.mpf(10) ** -(dps - 30))
        x['vec'] = [real(c, mp.mpf(10) ** -(dps - 30)) for c in x['vec']]
    res = eigenfunctions(ch, C, st)
    thm = [mpf(t) for t in ch.theta]
    bysec = {}
    for k, s in enumerate(st):
        s['rr'] = Q.raising_residual(ch, s)
        bysec.setdefault(s['sector'], []).append(k)
    br = []
    gaps = {}
    for mu, m in HW_SECT.items():
        ks = sorted(bysec[mu], key=lambda k: st[k]['rr'])
        sel = sorted(ks[:m], key=lambda k: st[k]['ev'])
        br += sel
        gaps[mu] = [st[k]['rr'] for k in ks[:m + 1]]
    S1 = Q.pf_S1(sum(mpf(z) for z in zeta), 3, thm)
    recs = []
    for k in br:
        s = st[k]
        S2 = Q.pf_S2(s['lam2'], thm)
        S3 = Q.pf_S3(s['lam3'], thm)
        dual = [Q.solve_poly_mp('D', S1, S2, S3, mpf(zeta[j]), thm) for j in range(3)]
        direct = [Q.solve_poly_mp('A', S1, S2, S3, -mpf(zeta[i]), thm, gauge=Q.GAUGE) for i in range(3)]
        recs.append(dict(sector=s['sector'], ev=s['ev'], lam2=[real(x) for x in s['lam2']], lam3=[real(x) for x in s['lam3']], vec=s['vec'], rr=s['rr'],
                         dual=[[real(c) for c in x[3]] for x in dual], direct=[[real(c) for c in x[3]] for x in direct],
                         degs=([x[0] for x in dual], [x[0] for x in direct], [x[1] for x in dual + direct])))
    out = dict(eps=eps, zeta=zeta, dps=dps, branches=recs, gaps=gaps, eigres=res)
    with open(cache, 'wb') as f:
        pickle.dump(out, f)
    return out


def formula_terms(st, zeta, w1, w1p, th=THETA):
    """returns dict: N1[5], N2[5], G (5x5 mp), c, M1M2, terms matrix, Qform = sum G N N /(4 c M1 M2) for one branch record"""
    sj = StateJets(st, zeta, th)
    N1, N2 = b_N(sj, 0), b_N(sj, 1)
    G2 = d2_G(Fr(w1), Fr(w1p), Fr(1) / (Fr(th[0]) - Fr(th[1])))
    G = [[mpf(Fr(int(G2[i][j].p), int(G2[i][j].q))) for j in range(5)] for i in range(5)]
    c = cscal(zeta, w1, w1p)
    M1 = N1[0]
    M2 = N2[0]
    terms = [[G[i][j] * N1[i] * N2[j] for j in range(5)] for i in range(5)]
    tot = mp.fsum(t for row in terms for t in row)
    return dict(sj=sj, N1=N1, N2=N2, G=G, c=c, M1M2=M1 * M2, terms=terms, tot=tot, Qform=tot / (4 * c * M1 * M2))


def hw_space(ch, mu):
    """exact (Fraction) basis of the total-highest-weight vectors of weight mu in the 64-dim space"""
    idx = ch.sectors[mu]
    E12 = ch.Es[(1, 1, 2)] + ch.Es[(2, 1, 2)]
    E23 = ch.Es[(1, 2, 3)] + ch.Es[(2, 2, 3)]
    rows = []
    for E in (E12, E23):
        for r, row in E.r.items():
            rows.append([sympy.Rational(row.get(j, 0).numerator, row.get(j, 0).denominator) if j in row else sympy.Integer(0) for j in idx])
    M = sympy.Matrix(rows)
    ns = M.nullspace()
    out = []
    for v in ns:
        full = [Fr(0)] * ch.dim
        for k, j in enumerate(idx):
            q = sympy.Rational(v[k])
            full[j] = Fr(int(q.p), int(q.q))
        out.append(full)
    return out


def gram_diag(ch):
    """positive diagonal G (one site, gl3_reps order) with <a|E_ij b> = <E_ji a|b> for the plain GT generators: d_a E_ij[a,b] = E_ji[b,a] d_b"""
    d = [None] * ch.d1
    d[ch.hw_index] = Fr(1)
    changed = True
    while changed:
        changed = False
        for (i, j) in ((1, 2), (2, 3)):
            Eu = ch.E1[(i, j)]
            El = ch.E1[(j, i)]
            for a, row in Eu.r.items():
                for b, v in row.items():
                    lb = El.r.get(b, {}).get(a)
                    if d[b] is not None and d[a] is None and v != 0:
                        d[a] = lb * d[b] / v
                        changed = True
                    elif d[a] is not None and d[b] is None and lb is not None and lb != 0:
                        d[b] = v * d[a] / lb
                        changed = True
    return d


def gnormalise(ch, v):
    """G-norm 1 (Gram diag d (x) d), sign: first nonzero component in the notebook (GT descending) order positive"""
    d = gram_diag(ch)
    n2 = mp.fsum(mpf(d[a // 8] * d[a % 8]) * v[a] ** 2 for a in range(64))
    v = [x / mp.sqrt(n2) for x in v]
    nb_order = [NB2G[a] * 8 + NB2G[b] for a in range(8) for b in range(8)]
    for idx in nb_order:
        if abs(v[idx]) > mp.mpf(10) ** -30:
            if v[idx] < 0:
                v = [-x for x in v]
            break
    return v


def richardson3(f1, f2, f3):
    """eliminate a, b in f(eps) = f0 + a eps + b eps^2 from eps, 2eps, 4eps"""
    return (8 * f1 - 6 * f2 + f3) / 3


def F4(dps=120):
    hdr('F4: untwisting regime, zeta = eps phi, w = eps w^, 120 digits; Q-side formula vs the chain; Q-side CG')
    mp.mp.dps = dps
    ch = chain()
    pts = {}
    for eps in EPSS:
        t0 = time.time()
        pts[eps] = build_regime(eps, dps)
        o = pts[eps]
        print('  eps = %s: eigen residual %s (%.0fs); branch sectors %s; raising residuals (selected, next) %s' % (
            eps, nst(o['eigres']), time.time() - t0, [b['sector'] for b in o['branches']],
            {mu: [nst(x) for x in g] for mu, g in o['gaps'].items()}), flush=True)
        print('        Q degrees (dual; direct; nullities) per branch:', [(b['degs'][0], b['degs'][1], 'all1' if set(b['degs'][2]) == {1} else b['degs'][2]) for b in o['branches']])
    grid = [(a, b) for a in OMH1 for b in OMH1P]
    # formula versus chain
    data = {}
    print('  Q-side formula  Qform = (1/4c) sum_ij G_ij N_i(th1) N_j(th2) / (M1 M2)   versus chain  Qnum = <Phi|Psi>/<x0^G|Psi>,  max relative deviation over the 4 grid points, per branch:')
    for eps in EPSS:
        o = pts[eps]
        zeta = o['zeta']
        row = []
        for bi, b in enumerate(o['branches']):
            devs = []
            depth = []
            for (a1, a2) in grid:
                w1, w1p = eps * a1, eps * a2
                ft = formula_terms(b, zeta, w1, w1p)
                cov = cov_to_mp(cov_inf(ch, zeta, w1, w1p))
                chainQ = overlap(cov, b['vec']) / overlap(x0_D2(ch, zeta), b['vec'])
                data[(eps, bi, a1, a2)] = (ft, chainQ)
                devs.append(abs(ft['Qform'] / chainQ - 1))
                mx = max(abs(t) for r_ in ft['terms'] for t in r_)
                depth.append(mp.log10(mx / abs(ft['tot'])))
            row.append((b['sector'], max(devs), max(depth)))
        print('     eps = %s: ' % eps + '; '.join('%s dev %s (cancellation %s digits)' % (s_, nst(d_), nst(dp_, 3)) for s_, d_, dp_ in row))
    # exponents
    print('  epsilon-exponents p (quantity ~ eps^p; from eps = 1e-6 -> 2e-6 | 2e-6 -> 4e-6), first grid point (w^1,w^1\') = (3/4,-12/5), branch by branch:')
    e1, e2, e3 = EPSS
    a1, a2 = grid[0]
    def ex(f):
        x1, x2, x3 = f(e1), f(e2), f(e3)
        return '%s|%s' % (nst(mp.log(abs(x2 / x1)) / mp.log(2), 4), nst(mp.log(abs(x3 / x2)) / mp.log(2), 4))
    names = ['P3(th1)', 'kappa_1', 'p02(1)', 'p03(1)', 'p23(1)', 'N1(1)..: |N|max', 'sum G N N', '|G N N| max term', 'c', 'M1M2', 'Qnum(chain)', '<x0|Psi>']
    for bi in range(6):
        sec = pts[e1]['branches'][bi]['sector']
        fts = lambda e: data[(e, bi, a1, a2)][0]
        vals = [lambda e: fts(e)['sj'].P[0][2][0], lambda e: fts(e)['sj'].P[0][2][1] / fts(e)['sj'].P[0][2][0],
                lambda e: plucker(fts(e)['sj'].P[0], 1, 2, 3)[(0, 2)], lambda e: plucker(fts(e)['sj'].P[0], 1, 2, 3)[(0, 3)],
                lambda e: plucker(fts(e)['sj'].P[0], 1, 2, 3)[(2, 3)], lambda e: max(abs(x) for x in fts(e)['N1']),
                lambda e: fts(e)['tot'], lambda e: max(abs(t) for r_ in fts(e)['terms'] for t in r_), lambda e: fts(e)['c'], lambda e: fts(e)['M1M2'],
                lambda e: data[(e, bi, a1, a2)][1], lambda e: overlap(x0_D2(ch, pts[e]['zeta']), pts[e]['branches'][bi]['vec'])]
        print('     %s: ' % (sec,) + ' '.join('%s %s' % (n, ex(v)) for n, v in zip(names, vals)))
    print('     kernel entries G_ij: eps-order of the entry at (w1,w1p) -> eps (w^) : table of lowest total degree (F3 (5)); thi = 1/theta12 is eps-free.')
    # ---- CG from the Q side
    print('  CG from the Q side alone: for each branch the limit eps -> 0 of the ratios rho_a = Qform(grid a)/Qform(grid 0) (Richardson in eps from 1,2,4e-6; error O(eps^3));')
    print('  with Psi_B -> v_B in the exact highest-weight space H_mu (untwisting): <Phi_a|v_B> = Qnum_a <x0|v_B> => the linear system (f_a - rho_a f_0)(v) = 0, f_a(v) = <Phi_a|v> (Phi_a exact, eps-free).')
    gn = None
    for bi in range(6):
        sec = pts[e1]['branches'][bi]['sector']
        rho = []
        rho_chain = []
        for gi in range(1, 4):
            r_q = [data[(e, bi, *grid[gi])][0]['Qform'] / data[(e, bi, *grid[0])][0]['Qform'] for e in EPSS]
            r_c = [data[(e, bi, *grid[gi])][1] / data[(e, bi, *grid[0])][1] for e in EPSS]
            rho.append(richardson3(*r_q))
            rho_chain.append(richardson3(*r_c))
        H = hw_space(ch, sec)
        m = len(H)
        phis = [cov_to_mp(cov_inf(ch, PHI, a, b)) for (a, b) in grid]      # Phi_a is zeta-independent in the regime (r_k eps-free): use zeta=phi, w = w^
        f = lambda gi, v: overlap(phis[gi], [mpf(x) for x in v])
        Mx = mp.matrix(3, m)
        for gi in range(1, 4):
            for i in range(m):
                Mx[gi - 1, i] = f(gi, H[i]) - rho[gi - 1] * f(0, H[i])
        S = mp.svd_r(Mx, compute_uv=False)
        scl = max(abs(f(gi, H[i])) for gi in range(4) for i in range(m))
        Sv = [S[i] / scl for i in range(len(S))]
        U, S2_, V = mp.svd_r(Mx)
        beta = [V[V.rows - 1, i] for i in range(m)]
        v = [mp.mpf(0)] * 64
        for i in range(m):
            for j in range(64):
                v[j] += beta[i] * mpf(H[i][j])
        vn = gnormalise(ch, v)
        # chain-side direction: projection of the branch eigenvector (eps = 1e-6) onto H via least squares, then G-normalised
        vec = pts[e1]['branches'][bi]['vec']
        Hm = mp.matrix(64, m)
        for i in range(m):
            for j in range(64):
                Hm[j, i] = mpf(H[i][j])
        cc = mp.lu_solve(Hm, mp.matrix(vec))
        vc = [mp.mpf(0)] * 64
        for i in range(m):
            for j in range(64):
                vc[j] += cc[i] * mpf(H[i][j])
        vcn = gnormalise(ch, vc)
        sl = lambda w: [(a, nst(w[a * 8 + 7], 12)) for a in range(8) if abs(w[a * 8 + 7]) > 1e-20]      # slice: second factor = HW ket (index 7)
        dirdev = max(abs(x - y) for x, y in zip(vn, vcn))
        print('   branch %s: HW space dim %d; consistency singular values of the %dx%d system %s; rho_a(Q side, limit) = %s ; chain rho_a = %s' % (
            sec, m, 3, m, [nst(x) for x in Sv], [nst(r, 12) for r in rho], [nst(r, 12) for r in rho_chain]))
        print('       Q-side vector (G-norm 1, sign: first nonzero GT component > 0): slice <a, hw|v> (a = gl3_reps index of the first factor): %s' % sl(vn))
        print('       chain-side eigenvector projected on H (eps=1e-6), same normalisation: max |v_Q - v_chain| = %s  (O(eps) expected); |v - Pv|/|v| of the chain vector = %s' % (
            nst(dirdev), nst(mp.norm(mp.matrix(vec) - Hm * cc) / mp.norm(mp.matrix(vec)))))
    return pts, data


# ------------------------------------------------------------------------------------------------ F5
def F5():
    hdr('F5: cross-checks with D2 (RESULTS-twovac-limit.md) and W2 (Gaudin-L2-chain-2.wb, A5-1 table of the [2,1,0] run)')
    # D2: kernel at the base point and at the (-2/5,3) point from the fit vs its printed analytic G (done at length in F2/F3); here the Qnum comparison with W2
    import json
    nbp = os.path.join(HERE, '..', '..', 'Mathematica', 'Gaudin', 'Experiments', 'Gaudin-L2-chain-2.wb')
    d = json.load(open(nbp, encoding='utf8'))
    cell = None
    for c in d['cells']:
        if c['value'].lstrip().startswith('**[2,1,0] run**') and 'A5-1' in c['value']:
            cell = c['value']
            break
    if cell is None:
        print('  W2 [2,1,0] table cell not found')
        return
    rows = []
    for line in cell.splitlines():
        cells_ = [x.strip() for x in line.split('|')]
        if len(cells_) > 10 and cells_[1].isdigit():
            rows.append(cells_)
    print('  W2 A5-1 rows parsed: %d' % len(rows))
    pt = get_point('P')
    ch = chain()
    zeta = pt['zeta']
    pairs = [(Fr(5, 2), Fr(7, 4)), (Fr(-2, 5), Fr(3))]
    x0 = x0_D2(ch, zeta)
    mine = []
    for k, st in enumerate(pt['states']):
        l2, l3 = S2S3_labels(st)
        ox = overlap(x0, st['vec'])
        q = [overlap(cov_to_mp(cov_inf(ch, zeta, *pr)), st['vec']) / ox for pr in pairs]
        mine.append((l2, l3, q, ox / max(abs(v) for v in st['vec'])))
    worst = [mp.mpf(0), mp.mpf(0)]
    unmatched = 0
    for r in rows:
        s2 = mp.mpf(r[4].replace('e', 'e'))
        s3 = mp.mpf(r[5])
        qw = [mp.mpf(r[9]), mp.mpf(r[11])]
        best = min(mine, key=lambda m: abs(m[0] - s2) + abs(m[1] - s3))
        if abs(best[0] - s2) + abs(best[1] - s3) > mp.mpf(10) ** -20:
            unmatched += 1
            continue
        for i in range(2):
            worst[i] = max(worst[i], abs(best[2][i] / qw[i] - 1))
    print('  matching by the (S2(3/7), S3(-2/9)) labels (tolerance 1e-20): unmatched %d of %d (the unmatched rows belong to other zeta points of the W2 table, e.g. the second zeta point)' % (unmatched, len(rows)))
    print('  max relative difference of Qnum(B) between this code (x0 of D2/W2, Phi_inf) and W2: (5/2,7/4): %s ; (-2/5,3): %s  (30 digits printed by W2)' % (nst(worst[0]), nst(worst[1])))
    print('  D2: kernel G fitted here equals D2\'s analytic G at 9 + 72 (w1,w1p) points, 6 theta values, c = 2 scaling (F2, F3); D2\'s printed vacuum monomial P3 p02 and constant 1 reproduced (F2 level a).')


STAGES = ['F0', 'F1', 'F2', 'F3', 'F4', 'F5']

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stages', default=','.join(STAGES))
    ap.add_argument('--build', default='')
    ap.add_argument('--workers', type=int, default=10)
    a = ap.parse_args()
    stg = a.stages.split(',')
    if a.build:
        build_all(a.build.split(','), workers=a.workers)
        return
    t0 = time.time()
    if 'F0' in stg:
        F0()
    if 'F1' in stg:
        F1()
    if 'F0' in stg:
        F0b()
    if 'F2' in stg:
        F2()
    if 'F3' in stg:
        F3()
    if 'F4' in stg:
        F4()
    if 'F5' in stg:
        F5()
    print('\ntotal %.0fs' % (time.time() - t0))

if __name__ == '__main__':
    main()
