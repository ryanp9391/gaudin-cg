"""P4/J: the [3,2,0] jet formula by the round-2 fitting pipeline (D1's twovac_fit.py generalised; Lambda = (3,2,0), L = 2, dim 225).

Run top to bottom:  python -X utf8 twovac_fit_320.py [--stages T0,T1,T2,T3,T4]     (caches twovac_fit_320_*.pkl)
  T0  points, objects (Phi_inf, x0^G), sanity: <x0|Psi> != 0, state-independence of <x0|Psi>/(R/D)/(Q1 P3 Q1 P3) (frame test)
  T1  jet facts: kappa_s = a_s/8, kappa2_s, 3x5 jet rank 2, p01 = p02 = p12 = 0, per-site feature ranks
  T2  hypothesis ladder (jet-Plucker bases, vacuum-monomial candidates, c exponent scan)
  T3  nu-monomial scheme, sparse representative, kernel
Imports twovac_fit.py (D1), qside_320.py, qside.py, gaudin_chain.py read-only.
"""
import sys, os, time, itertools, pickle, argparse
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mpmath as mp
import twovac_fit as tf
import qside as Q
import qside_320 as q3
import qside_310 as q310
from gaudin_chain import Chain, Operators, eigenstates, eigenfunctions

LAM = (3, 2, 0)
DPS = 80
THETA = tf.THETA
W1BASE = tf.W1BASE
mpf, hdr, nst = tf.mpf, tf.hdr, tf.nst
ZETAS = dict(tf.ZETAS)


def qm():
    """the qside module of the current rep (LAM): qside_320 or qside_310 (for the [3,1,0] comparison)"""
    return {(3, 2, 0): q3, (3, 1, 0): q310}[LAM]


def set_lam(lam):
    global LAM
    LAM = tuple(lam)


def chain320(th=THETA):
    key = tuple(Fr(t) for t in th)
    if key not in tf._CH or tf._CH[key].lam != LAM:
        tf._CH[key] = Chain(lam=LAM, theta=key)
    return tf._CH[key]


def pt_path(tag, dps, th):
    return os.path.join(HERE, 'twovac_fit_320_pt%s_%s_%d%s.pkl' % ('' if LAM == (3, 2, 0) else '310', tag, dps, tf.thtag(th)))


def build_point(zeta, dps=DPS, th=THETA):
    """states (225), dual/direct Q coefficient lists, Laurent data for one zeta point; plain picklable record"""
    mp.mp.dps = dps
    Q.TOL = mp.mpf(10) ** -50
    zeta = tuple(Fr(z) for z in zeta)
    ch = chain320(th)
    g = [[zeta[i] if i == j else Fr(0) for j in range(3)] for i in range(3)]
    ops = Operators(ch, g)
    thf = set(Fr(t) for t in th)
    pts = [Fr(3, 7), Fr(-2, 9), Fr(5, 3), Fr(7, 2), Fr(-4, 5), Fr(11, 4), Fr(-7, 3), Fr(13, 5), Fr(1, 3)]
    if any(u in thf for u in pts):
        pts = [Fr(17, 13), Fr(-19, 11), Fr(23, 17), Fr(29, 7), Fr(-31, 13), Fr(37, 19), Fr(-41, 17), Fr(43, 23), Fr(47, 29)]
    C = ops.fit(pts)
    u0, u1 = (Fr(3, 7), Fr(-2, 9)) if Fr(3, 7) not in thf and Fr(-2, 9) not in thf else (Fr(17, 13), Fr(-19, 11))
    st = eigenstates(ch, ops, C, u0=u0, u1=u1, dps=dps)
    for s in st:
        mxim = max([abs(mp.im(x)) for x in s['vec']] + [abs(mp.im(s['ev']))])
        if mxim > mp.mpf(10) ** (-dps // 2):
            raise ValueError('complex spectrum (imag %s)' % nst(mxim))
        s['ev'] = mp.re(s['ev'])
        s['vec'] = [mp.re(x) for x in s['vec']]
    res = eigenfunctions(ch, C, st)
    thm = [mpf(t) for t in ch.theta]
    recs = qm().q1_run(st, thm, list(zeta))
    out = []
    for s, r in zip(st, recs):
        assert tuple(len(x[3]) - 1 for x in r['dual']) == tuple(s['sector']), (s['sector'], [len(x[3]) - 1 for x in r['dual']])
        out.append(dict(sector=s['sector'], ev=s['ev'], lam2=list(s['lam2']), lam3=list(s['lam3']), vec=list(s['vec']),
                        dual=[list(x[3]) for x in r['dual']], direct=[list(x[3]) for x in r['direct']], resid=r['resid']))
    return dict(zeta=zeta, theta=tuple(Fr(t) for t in th), dps=dps, states=out, eigres=res)


def _build_job(a):
    tag, zeta, dps, th = a
    t0 = time.time()
    mp.mp.dps = dps
    pt = build_point(zeta, dps, th)
    with open(pt_path(tag, dps, th), 'wb') as f:
        pickle.dump(pt, f)
    return tag, time.time() - t0, str(pt['eigres'])[:12]


def build_many(tags_zetas, dps=DPS, th=THETA, workers=10):
    """build (in parallel processes) the points not yet cached"""
    todo = [(tag, z, dps, th) for tag, z in tags_zetas if not os.path.exists(pt_path(tag, dps, th))]
    if not todo:
        return
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=min(workers, len(todo))) as ex:
        for tag, dt, er in ex.map(_build_job, todo):
            print('  built point %s in %.0fs (eigen residual %s)' % (tag, dt, er), flush=True)


def get_point(tag, zeta=None, dps=DPS, th=THETA):
    path = pt_path(tag, dps, th)
    if not os.path.exists(path):
        build_many([(tag, zeta if zeta is not None else ZETAS[tag])], dps, th)
    with open(path, 'rb') as f:
        pt = pickle.load(f)
    mp.mp.dps = dps
    return pt


# ------------------------------------------------------------------------------------------------ normalisation scalars
def Delta12(zeta, w):
    z1, z2, z3 = [Fr(x) for x in zeta]
    w = Fr(w)
    return w ** 2 / ((z2 - w) * (z3 - w)), w / ((z1 - z2) * (z3 - w))


def kappaG320(zeta, w):
    """kappa^G = Delta1^(l1-l2) Delta2^(l2-l3) Delta3^l3 = Delta1 Delta2^2 for (3,2,0) (W2's Delta's, rep independent)"""
    D1, D2 = Delta12(zeta, w)
    return mpf(D1 * D2 ** 2)


def cscal320(zeta, w1, w1p, a=1, b=3, m=2):
    """c = prod_w (z2-w)^a (z3-w)^b (z1-z2)^m  (derived: 1/(kappa kappa') up to w-powers = (z2-w)(z3-w)^3 (z2-w')(z3-w')^3 (z1-z2)^4)"""
    z1, z2, z3 = [Fr(x) for x in zeta]
    f = lambda w: (z2 - w) ** a * (z3 - w) ** b
    return mpf(f(Fr(w1)) * f(Fr(w1p)) * (z1 - z2) ** (2 * m))


def targets(pt, w1, w1p, th=THETA):
    ch = chain320(th)
    zeta = pt['zeta']
    cov = tf.cov_to_mp(tf.cov_inf(ch, zeta, w1, w1p))
    x = tf.x0_onesite(ch, zeta)
    x0 = [a * b for a in x for b in x]
    return [(tf.overlap(cov, st['vec']), tf.overlap(x0, st['vec'])) for st in pt['states']]


class Data:
    """rows over (zeta point, state) at fixed (w1, w1p, theta): jets (to order 4), Qnum"""
    def __init__(self, tags, w1=W1BASE[0], w1p=W1BASE[1], th=THETA, dps=DPS):
        mp.mp.dps = dps
        self.rows = []
        self.th = th
        for tag in tags:
            pt = get_point(tag, dps=dps, th=th)
            mp.mp.dps = dps
            tg = targets(pt, w1, w1p, th)
            for k, st in enumerate(pt['states']):
                sj = tf.StateJets(st, pt['zeta'], th, n=4)
                op, ox = tg[k]
                self.rows.append(dict(tag=tag, k=k, sj=sj, st=st, zeta=pt['zeta'], w=(w1, w1p), Qnum=op / ox, Phi=op, X0=ox))
        self.n = len(self.rows)


# ------------------------------------------------------------------------------------------------ per-site features
PAIRS = [(0, 3), (0, 4), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]      # Plucker p_kl with the forced p01 = p02 = p12 = 0 removed (and l <= 4)


def pl_(sj, s):
    """Plucker coordinates of (P^2, P^3) jets at theta_s up to order 4, forced zeros set to exact 0"""
    pl = tf.plucker(sj.P[s], 1, 2, 4)
    for kl in (((0, 1), (0, 2), (1, 2)) if LAM == (3, 2, 0) else ((0, 1),)):      # forced zeros: P' = kappa P (and P'' = kappa2 P for (3,2,0))
        pl[kl] = mp.mpf(0)
    return pl


def b_B1(sj, s, cs=(0, 1, 2), pairs=PAIRS):
    P3 = sj.P[s][2]
    pl = pl_(sj, s)
    return [P3[c] * pl[kl] for kl in pairs for c in cs]


def B1_names(cs=(0, 1, 2), pairs=PAIRS):
    return [(kl, c) for kl in pairs for c in cs]


def Mvac(sj, s, kind='P3p03'):
    P3 = sj.P[s][2]
    pl = pl_(sj, s)
    return {'P3p02': P3[0] * pl[(0, 2)], 'P3p03': P3[0] * pl[(0, 3)], 'P3p04': P3[0] * pl[(0, 4)], "P3'p03": P3[1] * pl[(0, 3)], 'P3p13': P3[0] * pl[(1, 3)]}[kind]


def nu_vars(sj, s):
    """kappa = P3'/P3, kappa2 = P3''/P3, pi4 = p04/p03, pi34 = p34/p03 at theta_s"""
    P3 = sj.P[s][2]
    pl = pl_(sj, s)
    return P3[1] / P3[0], P3[2] / P3[0], pl[(0, 4)] / pl[(0, 3)], pl[(3, 4)] / pl[(0, 3)]


WTS = (1, 2, 1, 4)


def nu_monos(maxw, maxdeg=None):
    out = []
    for e in itertools.product(range(maxw + 1), repeat=4):
        w = sum(a * b for a, b in zip(e, WTS))
        if w <= maxw:
            out.append(e)
    return sorted(out, key=lambda e: (sum(a * b for a, b in zip(e, WTS)), e))


def b_nu(sj, s, monos):
    v = nu_vars(sj, s)
    return [mp.fprod([x ** e for x, e in zip(v, m)]) for m in monos]


# ------------------------------------------------------------------------------------------------ design / fitting
def _clean(f):
    mx = max(abs(a) for a in f)
    return [a if abs(a) > mx * mp.mpf(10) ** -45 else mp.mpf(0) for a in f]


def target_of(r, cexp=(1, 3, 2), M='P3p03', useM=True, withM_in_feats=True):
    """Y = c * Qnum * M1 * M2 (the features then carry P3^(c) p_kl)"""
    a, b, m = cexp
    y = cscal320(r['zeta'], r['w'][0], r['w'][1], a, b, m) * r['Qnum']
    if useM:
        y = y * Mvac(r['sj'], 0, M) * Mvac(r['sj'], 1, M)
    return y


def design(data, bas1, bas2, cexp=(1, 3, 2), M='P3p03', useM=True, rowscale=True, rows=None):
    A, y = [], []
    for i, r in enumerate(data.rows):
        if rows is not None and i not in rows:
            continue
        f1 = _clean(bas1(r['sj'], 0))
        f2 = _clean(bas2(r['sj'], 1))
        row = [a * b for a in f1 for b in f2]
        tv = target_of(r, cexp, M, useM)
        if rowscale:
            nr = mp.sqrt(tf.dot(row, row))
            row = [c / nr for c in row]
            tv = tv / nr
        A.append(row)
        y.append(tv)
    cols = [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]
    return cols, y


def fit(data, bas1, bas2, label='', reduce=True, **kw):
    """returns (resid, rank, nunknowns, kernel list (full-length, unknown order a-major), qr)"""
    if reduce:
        sel1, _ = reduce_basis(data, bas1, 0)
        sel2, _ = reduce_basis(data, bas2, 1)
        b1 = lambda sj, s, sel=sel1: [bas1(sj, s)[i] for i in sel]
        b2 = lambda sj, s, sel=sel2: [bas2(sj, s)[i] for i in sel]
    else:
        sel1 = sel2 = None
        b1, b2 = bas1, bas2
    cols, y = design(data, b1, b2, **kw)
    t0 = time.time()
    qr = tf.PivQR(cols)
    x, resid, rank = qr.lstsq(y)
    print('   %-46s unknowns %4d (per site %s/%s) rows %4d rank %4d null %4d resid %s  (%.0fs)' % (
        label, len(cols), len(sel1) if sel1 is not None else '-', len(sel2) if sel2 is not None else '-', len(y), rank, len(cols) - rank, nst(resid, 3), time.time() - t0), flush=True)
    return dict(resid=resid, rank=rank, n=len(cols), x=x, qr=qr, sel1=sel1, sel2=sel2, y=y)


def reduce_basis(data, bas, s, tol=mp.mpf(10) ** -40):
    F = [_clean(bas(r['sj'], s)) for r in data.rows]
    n = len(F[0])
    cols = []
    for j in range(n):
        c = [F[i][j] for i in range(len(F))]
        m = max(abs(x) for x in c)
        cols.append([x / m for x in c] if m != 0 else c)
    qr = tf.PivQR(cols)
    r = qr.rank(tol)
    return sorted(qr.perm[:r]), qr


# ================================================================================================ stages
THETA_EXTRA = [(Fr(9, 7), Fr(2, 5)), (Fr(1, 3), Fr(-1, 2)), (Fr(2, 3), Fr(-1, 7)), (Fr(1, 9), Fr(-5, 4)), (Fr(3, 5), Fr(-1, 2)), (Fr(1, 4), Fr(-3, 4)), (Fr(5, 6), Fr(1, 6)), (Fr(-2, 5), Fr(-8, 5))]
OM0 = (Fr(5, 2), Fr(7, 3), Fr(3, 5))
TAGS5 = ('P', 'Pp', 'X1', 'X2', 'X3')
TAGS7 = ('P', 'Pp', 'X1', 'X2', 'X3', 'X4', 'X5')


def detN_for_states(pt, om0=OM0, th=THETA, gauge=3):
    """det N (rows (s,(i,j)) in XXX triple) of each state against the exact vacuum Q's of om0, and the vacuum covector"""
    ch = chain320(th)
    zeta = pt['zeta']
    mp.mp.dps = pt['dps']
    Q.TOL = mp.mpf(10) ** -50
    thm = [mpf(t) for t in ch.theta]
    vd = q3.vac_data(ch, list(om0), list(zeta), gauge)
    assert vd['ok']
    Qs = [Q.Qfun(-mpf(om0[i]), [mpf(c) for c in vd['Qd'][i][2]]) for i in range(3)]
    tdV = Q.SiteData(Qs, thm, 8)
    z = [mpf(x) for x in zeta]
    dets = []
    for st in pt['states']:
        P = [Q.Qfun(z[j], st['dual'][j]) for j in range(3)]
        tdP = Q.SiteData(P, thm, 8)
        N = q3.Nrows(tdP, tdV, thm, gauge)
        dets.append(mp.det(Q.build_N_det(N, q3.XXX)))
    return dets, vd['cov']


def T0(tag='P'):
    hdr('T0: objects and sanity (point %s, Lambda = %s, %d digits)' % (tag, LAM, DPS))
    D = Data([tag])
    ch = chain320()
    pt = get_point(tag)
    print('zeta =', [str(x) for x in pt['zeta']], ' theta =', [str(x) for x in THETA], ' HW index', ch.hw_index, ' states', len(pt['states']),
          ' max ODE residual', nst(max(s['resid'] for s in pt['states'])), ' eigen residual', nst(pt['eigres']))
    X0 = [abs(r['X0']) for r in D.rows]; PH = [abs(r['Phi']) for r in D.rows]
    print('<x0^G|Psi_B> (max-component normalised kets): min %s max %s ; <Phi_inf|Psi_B> at (5/2,7/4): min %s max %s' % (nst(min(X0)), nst(max(X0)), nst(min(PH)), nst(max(PH))))
    # Phi_inf vs the vacuum covector at omega_2 = 1e6
    zeta = pt['zeta']
    from gaudin_chain import vacuum_covector_1site
    for w1 in (Fr(5, 2), Fr(7, 4)):
        v1 = vacuum_covector_1site(ch, zeta, (w1, Fr(10) ** 6, Fr(3, 5)))
        v2 = tf.onesite_inf(ch, zeta, w1)
        k = ch.hw_index
        d = max(abs(mpf(a / v1[k]) - mpf(b / v2[k])) for a, b in zip(v1, v2))
        print('  Omega^inf[zeta,%s] vs vacuum_covector_1site(omega=(%s,1e6,3/5)): max projective difference %s (O(1e-6))' % (w1, w1, nst(d)))
    # the frame test: <x0|Psi_B> = n_B * const * Q1(th1) P3(th1) Q1(th2) P3(th2)
    dets, vcov = detN_for_states(pt)
    mp.mp.dps = DPS
    vc = tf.cov_to_mp(vcov)
    Rv = [tf.overlap(vc, st['vec']) for st in pt['states']]
    vals = []
    for r, d, R in zip(D.rows, dets, Rv):
        sj = r['sj']
        M1 = sj.Q[0][0][0] * sj.P[0][2][0]
        M2 = sj.Q[1][0][0] * sj.P[1][2][0]
        nB = R / d
        vals.append((r['X0'] / nB / (M1 * M2), r['X0'] / nB))
    v0 = vals[0][0]
    print('frame test: <x0|Psi_B>/n_B/(Q1P3)(th1)(Q1P3)(th2), n_B = R_B/D_B(omega0): max_B |./. (B=0) - 1| = %s ; without the M1 M2 factor: %s' % (
        nst(max(abs(v[0] / v0 - 1) for v in vals)), nst(max(abs(v[1] / vals[0][1] - 1) for v in vals))))
    # the 225 Qnum at P, (5/2,7/4), 30 digits, with labels
    return D


def T1(tags=TAGS5):
    hdr('T1: jet facts')
    D = Data(list(tags))
    # forced Plucker zeros, kappa formulas, jet ranks
    mxp = mx1 = mx2 = mp.mpf(0)
    rk = {}
    for r in D.rows:
        sj, st = r['sj'], r['st']
        for s in range(2):
            pl = tf.plucker(sj.P[s], 1, 2, 4)
            sc = max(abs(v) for v in pl.values())
            mxp = max(mxp, max(abs(pl[kl]) for kl in ((0, 1), (0, 2), (1, 2))) / sc)
            a = st['lam3'][2 + 3 * s]; H = st['lam2'][1 + 2 * s]; b = st['lam3'][1 + 3 * s]
            P3 = sj.P[s][2]
            mx1 = max(mx1, abs(P3[1] / P3[0] - a / 8) / (abs(a / 8) + 1))
            mx2 = max(mx2, abs(P3[2] / P3[0] - (a * a - H * a + 8 * b) / 24) / (abs((a * a - H * a + 8 * b) / 24) + 1))
            M = mp.matrix([[sj.P[s][j][k] / mp.factorial(k) for k in range(5)] for j in range(3)])
            for j in range(3):
                nr = max(abs(M[j, c]) for c in range(5))
                for c in range(5):
                    M[j, c] /= nr
            S = tf.sv(M)
            key = sum(1 for x in S if x / S[0] > mp.mpf(10) ** -30)
            rk[key] = rk.get(key, 0) + 1
    print('  points', tags, ' rows', D.n, '(state, zeta) x 2 sites')
    print('  forced Pluckers p01 = p02 = p12: max |p|/max|p_kl| = %s' % nst(mxp))
    print('  kappa_s = P3\'/P3 = a_s/8: max rel dev %s ; kappa2_s = P3\'\'/P3 = (a^2 - H a + 8 b)/24: %s' % (nst(mx1), nst(mx2)))
    print('  3x5 jet matrix rank tally', rk)
    # per-site feature ranks
    for lab, bas, nm in (('B1 = {P3^(c) p_kl, c<=2, 7 pairs}', b_B1, 21),):
        s1, _ = reduce_basis(D, bas, 0)
        s2, _ = reduce_basis(D, bas, 1)
        print('  %s: per-site rank %d / %d (site 1), %d / %d (site 2)' % (lab, len(s1), nm, len(s2), nm))
    return D


def _lone(idx):
    def f(sj, s, cs=(0, 1, 2), pairs=PAIRS):
        P = sj.P[s][idx]
        pl = pl_(sj, s)
        return [P[c] * pl[kl] for kl in pairs for c in cs]
    return f


def T2(tags=TAGS5):
    hdr('T2: hypothesis ladder (base point (w1,w1p) = (5/2,7/4), theta = (2/7,-3/5); target Y = c Qnum M1 M2, c = prod_w (z2-w)(z3-w)^3 (z1-z2)^4 (both w), M = P3 p03)')
    D = Data(list(tags))
    print('rows %d = 225 states x %d zeta points %s' % (D.n, len(tags), tags))
    print(' ladder in the per-site basis (same on both sites); resid = relative ||A k - y||/||y||, 0.0 = below the 1e-90 floor:')
    fit(D, b_B1, b_B1, 'B1 = P3^(c) p_kl, c<=2, 7 pairs (21)')
    fit(D, lambda sj, s: b_B1(sj, s, cs=(0, 1)), lambda sj, s: b_B1(sj, s, cs=(0, 1)), 'B1, lone c<=1 (14)')
    fit(D, lambda sj, s: b_B1(sj, s, cs=(0,)), lambda sj, s: b_B1(sj, s, cs=(0,)), 'B1, lone c = 0 only (7)')
    pr3 = [(0, 3), (1, 3), (2, 3)]
    fit(D, lambda sj, s: b_B1(sj, s, pairs=pr3), lambda sj, s: b_B1(sj, s, pairs=pr3), 'p_kl with l <= 3 only (9)')
    pr4 = [(0, 3), (0, 4), (1, 3), (1, 4), (2, 3), (2, 4)]
    fit(D, lambda sj, s: b_B1(sj, s, pairs=pr4), lambda sj, s: b_B1(sj, s, pairs=pr4), 'without p34 (18)')
    print(' vacuum-monomial candidates M (same B1 basis):')
    for M in ('P3p03', 'P3p04', "P3'p03"):
        fit(D, b_B1, b_B1, 'M = %s' % M, M=M)
    fit(D, b_B1, b_B1, 'no M (target c Qnum)', useM=False)
    print(' exponents of the frame-tie scalar c = prod_w (z2-w)^a (z3-w)^b (z1-z2)^(2m):')
    for ce in ((1, 3, 2), (1, 2, 2), (1, 3, 1), (0, 3, 2), (2, 3, 2), (1, 3, 3), (1, 4, 2)):
        fit(D, b_B1, b_B1, 'c exponents (a,b,m) = %s' % (ce,), cexp=ce)
    print(' controls (must fail):')
    fit(D, _lone(1), _lone(1), 'lone index P2 instead of P3')
    fit(D, _lone(0), _lone(0), 'lone index P1 instead of P3')
    fit(D, lambda sj, s: b_B1(sj, s), lambda sj, s: [x for x in b_B1(sj, 0)], 'both factors from the site-1 jets')
    fit(D, lambda sj, s: [x for x in b_B1(sj, 1)], lambda sj, s: [x for x in b_B1(sj, 1)], 'both factors from the site-2 jets')
    return D


N7NAMES = ['P3 p03', 'P3 p04', "P3' p03", "P3' p04", "P3'' p03", "P3'' p04", 'P3 p34']


def b_N7(sj, s):
    P3 = sj.P[s][2]
    pl = pl_(sj, s)
    return [P3[0] * pl[(0, 3)], P3[0] * pl[(0, 4)], P3[1] * pl[(0, 3)], P3[1] * pl[(0, 4)], P3[2] * pl[(0, 3)], P3[2] * pl[(0, 4)], P3[0] * pl[(3, 4)]]


def rat(x, maxden=10 ** 25, tol=mp.mpf(10) ** -55):
    f = Fr(int(mp.nint(x * 10 ** 70)), 10 ** 70).limit_denominator(maxden)
    return f if abs(mpf(f) - x) < tol * max(1, abs(x)) else None


def fitN7(D, **kw):
    cols, y = design(D, b_N7, b_N7, **kw)
    qr = tf.PivQR(cols)
    x, resid, rank = qr.lstsq(y)
    return x, resid, rank, qr


def kernel_matrix(x, n=7):
    return [[x[i * n + j] for j in range(n)] for i in range(n)]


def T3(tags=TAGS5):
    hdr('T3: the sparse formula: N7 = {P3^(c) p03, P3^(c) p04 (c = 0,1,2), P3 p34} per site (nu = (1, pi4) x (1, kappa, kappa2) + pi34)')
    D = Data(list(tags))
    x, resid, rank, qr = fitN7(D)
    print('  N7 x N7 fit at the base point (w1,w1p) = (5/2,7/4), theta12 = 31/35: %d rows, 49 unknowns, rank %d, resid %s' % (D.n, rank, nst(resid)))
    K = kernel_matrix(x)
    thi = mp.mpf(35) / 31
    print('  kernel K_ij in the basis (%s):' % ', '.join(N7NAMES))
    nonrat = 0
    for i in range(7):
        row = []
        for j in range(7):
            r = rat(K[i][j])
            nonrat += r is None
            row.append(str(r) if r is not None else nst(K[i][j], 12))
        print('   ', row)
    print('  rational reconstructions failed for %d of 49 entries' % nonrat)
    # the same kernel in the pure Plucker basis P3 p_kl (lone c = 0 suffices: P3' p0l = P3 p1l, P3'' p0l = P3 p2l, l = 3,4)
    PB = lambda sj, s: b_B1(sj, s, cs=(0,))
    cols, y = design(D, PB, PB)
    qrp = tf.PivQR(cols)
    xp, rp, rkp = qrp.lstsq(y)
    Kp = kernel_matrix(xp)
    names = ['P3 p%d%d' % kl for kl in PAIRS]
    print('  Plucker-basis fit ({P3 p_kl : (k,l) in %s}, 7 per site, lone c = 0 only): rank %d, resid %s; kernel in the basis (%s):' % (PAIRS, rkp, nst(rp), ', '.join(names)))
    for i in range(7):
        print('   ', [str(rat(Kp[i][j])) if rat(Kp[i][j]) is not None else nst(Kp[i][j], 12) for j in range(7)])
    print('  (identities used: dP3 p0l = P3 p1l and ddP3 p0l = P3 p2l for l = 3,4, exact to ~1e-72)')
    # per-zeta fits
    print('  per-zeta-point fits (225 rows, 49 unknowns) vs the stacked kernel: max |K_zeta - K| / |K|')
    Kn = max(abs(v) for v in x)
    for tg in tags:
        Dz = Data([tg])
        xz, rz, rkz, _ = fitN7(Dz)
        print('     %-3s rank %d resid %s max rel diff %s' % (tg, rkz, nst(rz), nst(max(abs(a - b) for a, b in zip(xz, x)) / Kn)))
    # minimality: drop one monomial (both sites)
    print('  minimality: drop one of the 7 per-site monomials (same on both sites) -> residual')
    for i in range(7):
        keep = [k for k in range(7) if k != i]
        bb = lambda sj, s, keep=keep: [b_N7(sj, s)[k] for k in keep]
        cols, y = design(D, bb, bb)
        q2 = tf.PivQR(cols)
        _, r2, rk2 = q2.lstsq(y)
        print('     without %-8s resid %s (rank %d / 36)' % (N7NAMES[i], nst(r2), rk2))
    # site exchange symmetry: K(w1,w1p) vs K(w1p,w1)^T
    D2 = Data(list(tags), w1=W1BASE[1], w1p=W1BASE[0])
    x2, r2, rk2, _ = fitN7(D2)
    K2 = kernel_matrix(x2)
    print('  site exchange: kernel at (w1,w1p) = (7/4,5/2): resid %s; max |K(w1,w1p)_ij - K(w1p,w1)_ji| / max|K| = %s' % (nst(r2), nst(max(abs(K[i][j] - K2[j][i]) for i in range(7) for j in range(7)) / Kn)))
    return D, K


def T4(tags=TAGS5):
    hdr('T4: kernel laws: the (w1,w1p) grid and polynomial dependence')
    gp = os.path.join(HERE, 'twovac_fit_320_grid.pkl')
    if os.path.exists(gp):
        with open(gp, 'rb') as f:
            Ks = pickle.load(f)
        print('  loaded the cached grid of %d kernel fits (delete twovac_fit_320_grid.pkl to recompute); each: 49 unknowns, rank 49, resid 0.0 at the original run' % len(Ks))
        return Ks
    ws = [Fr(5, 2), Fr(-2, 5), Fr(9, 7), Fr(3), Fr(-5, 3), Fr(7, 4), Fr(11, 3), Fr(-7, 2), Fr(13, 5)]
    pts, Ks = [], {}
    for w1 in ws:
        for w1p in ws:
            if w1 == w1p:
                continue
            D = Data(list(tags), w1=w1, w1p=w1p)
            x, resid, rank, _ = fitN7(D)
            Ks[(w1, w1p)] = kernel_matrix(x)
            print('   (w1,w1p) = (%s,%s): rank %d resid %s' % (w1, w1p, rank, nst(resid)), flush=True)
    with open(os.path.join(HERE, 'twovac_fit_320_grid.pkl'), 'wb') as f:
        pickle.dump(Ks, f)
    return Ks


def polyfit(xs, ys, monos):
    """least squares of ys on the monomials (functions of the point); returns (coeffs, relative residual, rank)"""
    cols = [[mp.fprod([xv ** e for xv, e in zip(x, m)]) for x in xs] for m in monos]
    qr = tf.PivQR(cols)
    c, resid, rank = qr.lstsq(ys)
    return c, resid, rank


def T5():
    hdr('T5: polynomial dependence of the kernel on (w1,w1p) (72 ordered pairs of 9 values; total-degree fits; fixed theta, zeta-independent)')
    with open(os.path.join(HERE, 'twovac_fit_320_grid.pkl'), 'rb') as f:
        Ks = pickle.load(f)
    keys = sorted(Ks)
    xs = [(mpf(a), mpf(b)) for a, b in keys]
    degs = {}
    allrat = True
    for i in range(7):
        for j in range(7):
            ys = [Ks[k][i][j] for k in keys]
            found = None
            for d in range(0, 11):
                monos = [(a, b) for a in range(d + 1) for b in range(d + 1 - a)]
                if len(monos) > len(keys) - 4:
                    break
                c, resid, rank = polyfit(xs, ys, monos)
                if resid < mp.mpf(10) ** -60:
                    found = (d, c, monos, resid)
                    break
            if found is None:
                degs[(i, j)] = None
                continue
            d, c, monos, resid = found
            rr = [rat(v) for v in c]
            allrat = allrat and all(r is not None or abs(v) < mp.mpf(10) ** -50 for r, v in zip(rr, c))
            degs[(i, j)] = (d, max(a for (a, b), v in zip(monos, c) if abs(v) > mp.mpf(10) ** -50) if any(abs(v) > 1e-50 for v in c) else 0,
                            max(b for (a, b), v in zip(monos, c) if abs(v) > mp.mpf(10) ** -50) if any(abs(v) > 1e-50 for v in c) else 0)
    print('  minimal total degree d_ij (7x7; None = not polynomial of degree <= 10 on these 72 points):')
    for i in range(7):
        print('    ', [degs[(i, j)][0] if degs[(i, j)] else None for j in range(7)])
    print('  max degree in w1 / in w1p per entry (7x7):')
    for i in range(7):
        print('    ', [(degs[(i, j)][1], degs[(i, j)][2]) if degs[(i, j)] else None for j in range(7)])
    print('  all coefficients rational (denominator <= 1e25):', allrat)
    return Ks, degs


def T6():
    hdr('T6: theta dependence of the kernel (single zeta point P per theta, base (w1,w1p) = (5/2,7/4))')
    ths = [THETA] + list(tf.THETA_VARS) + list(THETA_EXTRA)
    data = []
    for th in ths:
        if not os.path.exists(pt_path('P', DPS, th)):
            print('   theta', th, 'missing point; skipped')
            continue
        D = Data(['P'], th=th)
        x, resid, rank, _ = fitN7(D)
        t12 = Fr(th[0]) - Fr(th[1])
        data.append((th, t12, kernel_matrix(x)))
        print('   theta = (%s,%s) theta12 = %s: rank %d resid %s' % (th[0], th[1], t12, rank, nst(resid)), flush=True)
    # theta12-only dependence: pairs with equal theta12 and different theta
    byt = {}
    for th, t12, K in data:
        byt.setdefault(t12, []).append((th, K))
    for t12, v in byt.items():
        if len(v) > 1:
            d = max(abs(v[0][1][i][j] - v[1][1][i][j]) for i in range(7) for j in range(7)) / max(abs(v[0][1][i][j]) for i in range(7) for j in range(7))
            print('  theta12 = %s reached by %d thetas %s: max relative kernel difference %s' % (t12, len(v), [tuple(map(str, a)) for a, _ in v], nst(d)))
    # polynomial in thi = 1/theta12
    pts = {}
    for th, t12, K in data:
        pts.setdefault(t12, K)
    ts = sorted(pts)
    xs = [(mpf(1 / t),) for t in ts]
    print('  %d distinct theta12 values' % len(ts))
    for dmax in range(2, len(ts) - 1):
        monos = [(a,) for a in range(dmax + 1)]
        res = []
        for i in range(7):
            for j in range(7):
                ys = [pts[t][i][j] for t in ts]
                c, resid, rank = polyfit(xs, ys, monos)
                res.append(resid)
        print('   polynomial in thi = 1/theta12 of degree <= %d (%d coefficients, %d points): max resid over 49 entries %s' % (dmax, dmax + 1, len(ts), nst(max(res))))
    return data


def T7(tags=('P', 'Pp', 'X1')):
    """the same pipeline for [3,1,0] (lone c <= 1, M = P3 p02, c = (z2-w)^2 (z3-w)^3 (z1-z2)^2 per w): ladder + the N5-type and extended sparse bases"""
    set_lam((3, 1, 0))
    hdr('T7: [3,1,0] comparison (same pipeline; exponents {0,2,5}: only dP = kappa P forced, p01 = 0; first nonzero p02)')
    D = Data(list(tags))
    ce = (2, 3, 1)
    def mk(cs, pairs):
        return lambda sj, s: b_B1(sj, s, cs=cs, pairs=pairs)
    allp = [(0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
    print('  rows %d = 225 states x %d zeta points' % (D.n, len(tags)))
    r = fit(D, mk((0, 1), allp), mk((0, 1), allp), 'full: P3^(c) p_kl, c<=1, 9 pairs (18)', cexp=ce, M='P3p02')
    fit(D, mk((0, 1, 2), allp), mk((0, 1, 2), allp), 'full with lone c<=2 (27)', cexp=ce, M='P3p02')
    n5 = lambda sj, s: [sj.P[s][2][0] * pl_(sj, s)[(0, 2)], sj.P[s][2][1] * pl_(sj, s)[(0, 2)], sj.P[s][2][0] * pl_(sj, s)[(0, 3)], sj.P[s][2][1] * pl_(sj, s)[(0, 3)], sj.P[s][2][0] * pl_(sj, s)[(2, 3)]]
    fit(D, n5, n5, 'N5 = {P3^(c) p02, P3^(c) p03 (c=0,1), P3 p23}', reduce=False, cexp=ce, M='P3p02')
    n7 = lambda sj, s: n5(sj, s) + [sj.P[s][2][0] * pl_(sj, s)[(0, 4)], sj.P[s][2][1] * pl_(sj, s)[(0, 4)]]
    fit(D, n7, n7, 'N5 + {P3 p04, dP3 p04}', reduce=False, cexp=ce, M='P3p02')
    fit(D, n5, n5, 'N5, wrong c exponents (1,3,2)', reduce=False, cexp=(1, 3, 2), M='P3p02')
    # the 320-shaped hypothesis p_{0,e}, p_{0,e+1}, p_{e,e+1} with e = 2 is N5; test dropping each of the 5
    for i in range(5):
        keep = [k for k in range(5) if k != i]
        bb = lambda sj, s, keep=keep: [n5(sj, s)[k] for k in keep]
        cols, y = design(D, bb, bb, cexp=ce, M='P3p02')
        _, rr, rk = tf.PivQR(cols).lstsq(y)
        print('     N5 without element %d: resid %s (rank %d)' % (i, nst(rr), rk))
    cols, y = design(D, n5, n5, cexp=ce, M='P3p02')
    qr = tf.PivQR(cols)
    x, rr, rk = qr.lstsq(y)
    print('  N5 kernel (basis P3 p02, dP3 p02, P3 p03, dP3 p03, P3 p23), rational reconstruction:')
    K = kernel_matrix(x, 5)
    for i in range(5):
        print('   ', [str(rat(K[i][j])) if rat(K[i][j]) is not None else nst(K[i][j], 12) for j in range(5)])
    # the pure Plucker basis: P3 p_kl, all 9 non-forced pairs (k<l<=4), lone c = 0 only
    pb = lambda sj, s: b_B1(sj, s, cs=(0,), pairs=allp)
    cols, y = design(D, pb, pb, cexp=ce, M='P3p02')
    qr9 = tf.PivQR(cols)
    x9, r9, rk9 = qr9.lstsq(y)
    print('  pure Plucker basis {P3 p_kl: 9 pairs, k<l<=4, (0,1) forced zero}, lone c = 0: unknowns 81, rank %d, resid %s (unique kernel)' % (rk9, nst(r9)))
    pb5 = lambda sj, s: b_B1(sj, s, cs=(0,), pairs=[(0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])
    cols, y = design(D, pb5, pb5, cexp=ce, M='P3p02')
    _, r5, rk5 = tf.PivQR(cols).lstsq(y)
    print('  l <= 3 subset (5 Plucker coordinates = the adjoint pattern): rank %d resid %s -> FAILS' % (rk5, nst(r5)))
    for i in range(9):
        keep = [k for k in range(9) if k != i]
        bb = lambda sj, s, keep=keep: [pb(sj, s)[k] for k in keep]
        cols, y = design(D, bb, bb, cexp=ce, M='P3p02')
        _, rr, rk = tf.PivQR(cols).lstsq(y)
        print('     without P3 p%d%d: resid %s (rank %d / 64)' % (allp[i][0], allp[i][1], nst(rr), rk))
    K9 = kernel_matrix(x9, 9)
    nrat = sum(1 for i in range(9) for j in range(9) if rat(K9[i][j]) is None)
    print('  kernel entries with no rational reconstruction (denominator <= 1e25): %d of 81; K_99 = %s' % (nrat, rat(K9[8][8])))
    set_lam((3, 2, 0))
    return D


def T6b():
    hdr('T6b: site-exchange symmetry: is K_ij(w,w,thi) = K_ji(w,w,-thi) ? (w1 = w1p = 5/2; degree-7 polynomial in thi from 12 theta12 values)')
    ths = [THETA] + list(tf.THETA_VARS) + list(THETA_EXTRA)
    pts = {}
    for th in ths:
        if not os.path.exists(pt_path('P', DPS, th)):
            continue
        D = Data(['P'], w1=Fr(5, 2), w1p=Fr(5, 2), th=th)
        x, resid, rank, _ = fitN7(D)
        pts.setdefault(Fr(th[0]) - Fr(th[1]), kernel_matrix(x))
    ts = sorted(pts)
    monos = [(a,) for a in range(8)]
    xs = [(mpf(1 / t),) for t in ts]
    coef = {}
    worst = mp.mpf(0)
    for i in range(7):
        for j in range(7):
            c, resid, rank = polyfit(xs, [pts[t][i][j] for t in ts], monos)
            coef[(i, j)] = c
            worst = max(worst, resid)
    print('  %d theta12 values, degree-7 polynomial fit in thi: max resid %s' % (len(ts), nst(worst)))
    thi = mp.mpf(35) / 31
    ev = lambda c, x: sum(ci * x ** k for k, ci in enumerate(c))
    Kp = [[ev(coef[(i, j)], thi) for j in range(7)] for i in range(7)]
    Km = [[ev(coef[(i, j)], -thi) for j in range(7)] for i in range(7)]
    n = max(abs(Kp[i][j]) for i in range(7) for j in range(7))
    print('  max |K(thi)_ij - K(-thi)_ji| / max|K| = %s   (0 would mean exchange symmetric: K(w,w,thi)^T = K(w,w,-thi))' % nst(max(abs(Kp[i][j] - Km[j][i]) for i in range(7) for j in range(7)) / n))
    print('  max |K(thi)_ij + K(-thi)_ji| / max|K| = %s' % nst(max(abs(Kp[i][j] + Km[j][i]) for i in range(7) for j in range(7)) / n))
    return pts


STAGES = ['T0', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T6b', 'T7']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stages', default=','.join(STAGES))
    a = ap.parse_args()
    t0 = time.time()
    for st in a.stages.split(','):
        if st in globals():
            globals()[st]()
    print('\ntotal time %.0fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
