"""G1: general-Lambda structure of the kernel of the two-vacua jet formula (sl(3) Gaudin).  Python only, read-only imports of the earlier files.

Stages (run heavy ones separately):   python -X utf8 general_kernel.py --stages K1,K0,...
  K1  one-site kernel theta_i(a; Lambda), EXACT (Fractions), several Lambda
Notation: one-site monomials m_kl = P3 p_kl (p_kl = Plucker coordinate of (P2,P3) jets at theta), (kl) non-forced pairs
0 <= k < l <= lam1+1, not both <= lam2; M = P3 p_{0,lam2+1}.  The common exponential factor e^{(z2+2z3)theta} of all monomials is dropped
(the formula is homogeneous of degree 1 in the monomials: M on the left), so every quantity is an exact rational.
Jet weight of (kl): w = k+l ; e_kl = k + l - (lam2+1).  Canonical kernel Theta = K/(w w')^{2 lam1 - lam2}.
"""
import sys, os, time, itertools, pickle, argparse
from fractions import Fraction as Fr
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
import sympy
import mpmath as mp
import qside as Q
from gaudin_chain import SM, sm_from_sympy, comm, symbolic_cdet, exp_nil
from gl3_reps import GL3Rep


def hdr(t):
    print('\n' + '=' * 100 + '\n' + t + '\n' + '=' * 100, flush=True)


# ====================================================================================== exact linear algebra
def exact_solve(A, y):
    """solve A x = y over Q (list of rows of Fractions). returns (x or None, rank, consistent). Free variables set to 0."""
    n = len(A[0])
    M = [list(map(Fr, r)) + [Fr(v)] for r, v in zip(A, y)]
    piv = []
    r = 0
    for c in range(n):
        p = None
        for i in range(r, len(M)):
            if M[i][c] != 0:
                p = i
                break
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        inv = 1 / M[r][c]
        M[r] = [v * inv for v in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        piv.append(c)
        r += 1
        if r == len(M):
            break
    cons = all(M[i][n] == 0 for i in range(r, len(M)))
    x = [Fr(0)] * n
    for i, c in enumerate(piv):
        x[c] = M[i][n]
    return (x if cons else None), len(piv), cons


# ====================================================================================== one-site rep, operators, states
class One:
    """one-site gl(3) rep Lambda in the GT basis (rational, non-orthonormal): E dict (i,j)->SM, weights, hw index"""
    def __init__(self, lam):
        self.lam = tuple(lam)
        rep = GL3Rep(self.lam)
        self.rep = rep
        d = rep.dim
        self.d = d
        E = {(i, j): sm_from_sympy(rep.generator(i, j)) for (i, j) in [(1, 1), (2, 2), (3, 3), (1, 2), (2, 3), (2, 1), (3, 2)]}
        E[(1, 3)] = comm(E[(1, 2)], E[(2, 3)])
        E[(3, 1)] = comm(E[(3, 2)], E[(2, 1)])
        self.E = E
        hw = tuple(tuple(self.lam[:k]) for k in range(1, 4))
        self.hw_index = rep.patterns.index(hw)
        self.wt = [tuple(int(E[(k, k)].r.get(a, {}).get(a, 0)) for k in (1, 2, 3)) for a in range(d)]
        self.sectors = {}
        for a, w in enumerate(self.wt):
            self.sectors.setdefault(w, []).append(a)
        self.d1 = d
        self.E1 = E          # shim for twovac_fit.onesite_inf (needs .E1, .d1, .hw_index)


_ONE = {}
def one(lam):
    lam = tuple(lam)
    if lam not in _ONE:
        _ONE[lam] = One(lam)
    return _ONE[lam]


def S_at(R, g, theta, u):
    """(S1, S2, S3) as SM at rational u; one site L_ij = g_ij + E_ij/(u-theta)"""
    Sy = symbolic_cdet(False)
    u = Fr(u)
    cache = {}
    def L(i, j, dv):
        k = (i, j, dv)
        if k not in cache:
            c = Fr((-1) ** dv * (1 if dv == 0 else sympy.factorial(dv))) / (u - theta) ** (dv + 1)
            out = R.E[(i, j)] * c
            if dv == 0 and g[i - 1][j - 1] != 0:
                out = out + SM.eye(R.d, g[i - 1][j - 1])
            cache[k] = out
        return cache[k]
    pc = {(): None}
    def prod(mono):
        if mono in pc:
            return pc[mono]
        head = prod(mono[:-1])
        last = L(*mono[-1])
        r = last if head is None else head * last
        pc[mono] = r
        return r
    res = []
    for S in Sy:
        tot = SM(R.d)
        for mono, c in S.items():
            tot = tot + (SM.eye(R.d, c) if mono == () else prod(mono) * c)
        res.append(tot)
    return res


def fit_pf(R, zeta, theta, npts=(4, 5)):
    """exact partial-fraction coefficient operators: S2 = sum_b a^b C2[b] (b=0..2), S3 = sum_b a^b C3[b] (b=0..3), a = 1/(u-theta)"""
    g = [[Fr(zeta[i]) if i == j else Fr(0) for j in range(3)] for i in range(3)]
    pts = [Fr(3, 7), Fr(-2, 9), Fr(5, 3), Fr(7, 2), Fr(-4, 5), Fr(11, 4), Fr(-7, 3), Fr(13, 5)]
    pts = [p for p in pts if p != theta][:7]
    data = [S_at(R, g, theta, u) for u in pts]
    out = {}
    for k, nb in ((2, 3), (3, 4)):
        Phi = sympy.Matrix([[sympy.Rational((1 / (p - theta)) ** b) for b in range(nb)] for p in pts[:nb]])
        Pinv = Phi.inv()
        C = []
        for b in range(nb):
            op = SM(R.d)
            for kk in range(nb):
                w = Pinv[b, kk]
                if w != 0:
                    op = op + data[kk][k - 1] * Fr(int(w.p), int(w.q))
            C.append(op)
        bad = 0
        for kk in range(nb, len(pts)):
            pred = SM(R.d)
            for b in range(nb):
                pred = pred + C[b] * ((1 / (pts[kk] - theta)) ** b)
            if not (pred - data[kk][k - 1]).is_zero():
                bad += 1
        out[k] = C
        out['res%d' % k] = bad
    out['S1'] = data[0][0]
    return out


def diag_states(R, C):
    """GT basis vectors that are exact eigenvectors of all C2, C3 coefficient operators: dict idx -> (lam2 list, lam3 list)"""
    res = {}
    for a in range(R.d):
        ok = True
        lam = {2: [], 3: []}
        for k in (2, 3):
            for op in C[k]:
                col = {r: row[a] for r, row in op.r.items() if a in row}
                if any(r != a for r in col):
                    ok = False
                lam[k].append(col.get(a, Fr(0)))
        if ok:
            res[a] = (lam[2], lam[3])
    return res


def solve_dual(R, zeta, theta, lam2, lam3, nu):
    """exact dual Q polynomials p_j (coefficient lists, monic) for one-site eigenvalue functions"""
    th = [theta]
    S1 = Q.PF(sum(Fr(z) for z in zeta), {(0, 1): Fr(sum(R.lam))}, th)
    S2 = Q.PF(lam2[0], {(0, 1): lam2[1], (0, 2): lam2[2]}, th)
    S3 = Q.PF(lam3[0], {(0, 1): lam3[1], (0, 2): lam3[2], (0, 3): lam3[3]}, th)
    out = []
    for j in range(3):
        d, nul, q = Q.solve_poly_exact('D', S1, S2, S3, Fr(zeta[j]), th, dmax=sum(R.lam) + 2)
        out.append((d, nul, q))
    return out


def poly_jets(kappa, q, theta, n):
    """r_k(theta), (d/du)^k (e^{kappa u} q) = e^{kappa u} r_k"""
    cur = list(q)
    out = [Q.peval(cur, theta)]
    for _ in range(n):
        cur = Q.padd(Q.pder(cur), Q.pscale(cur, kappa))
        out.append(Q.peval(cur, theta))
    return out


def pairs_of(lam):
    l1, l2 = lam[0], lam[1]
    return [(k, l) for k in range(l1 + 2) for l in range(k + 1, l1 + 2) if not (l <= l2)]


def Mpair(lam):
    return (0, lam[1] + 1)


def mono_values(lam, J):
    """J[j][k] jets (j = 0,1,2 for P1,P2,P3). returns dict (kl) -> P3 p_kl (exponentials dropped), over non-forced pairs"""
    P2, P3 = J[1], J[2]
    return {(k, l): P3[0] * (P2[k] * P3[l] - P3[k] * P2[l]) for (k, l) in pairs_of(lam)}


# ---------- covectors (exact, Fractions)
def omega_inf(R, zeta, w):
    z1, z2, z3 = [Fr(x) for x in zeta]
    w = Fr(w)
    r1 = (w - z1) / (w - z2)
    r2 = -(w - z1) * (z2 - z3) / ((z1 - z2) * (w - z3))
    r3 = (w - z2) * (z1 - z3) / ((z1 - z2) * (w - z3))
    e = [Fr(0)] * R.d
    e[R.hw_index] = Fr(1)
    v = exp_nil(R.E[(1, 2)], r1, 40).rapply(e)
    v = exp_nil(R.E[(1, 3)], r2, 40).rapply(v)
    v = exp_nil(R.E[(2, 3)], r3, 40).rapply(v)
    return v


def x_cov(R, zeta):
    """x(zeta) = e_hw^T P_phys^{-1} (exact), same formula as twovac_fit.x0_onesite"""
    z1, z2, z3 = [Fr(x) for x in zeta]
    E = R.E
    s = [z1 ** 2 / ((z1 - z2) * (z1 - z3)), -(z2 / (z1 * (z2 - z3))), 1 / (z1 * z2)]
    r = [-(z3 / z2), -(z3 / z1), (z2 ** 2 * (z3 - z1)) / (z1 ** 2 * (z2 - z3))]
    x = [Fr(0)] * R.d
    x[R.hw_index] = Fr(1)
    for (c, key) in ((-(z1 + z2) / (z1 * z2), (2, 3)), (1 / (z1 * z2), (1, 3)), (-(1 / z1), (1, 2))):
        x = exp_nil(E[key], c, 40).rapply(x)
    for p in range(R.d):
        f = Fr(1)
        for k in range(3):
            f *= s[k] ** (-R.wt[p][k])
        x[p] *= f
    for (c, key) in ((-r[2], (2, 1)), (-r[1], (3, 1)), (-r[0], (3, 2))):
        x = exp_nil(E[key], c, 40).rapply(x)
    return x


def c1_scalar(lam, zeta, w):
    z1, z2, z3 = [Fr(x) for x in zeta]
    w = Fr(w)
    return (z2 - w) ** (lam[0] - lam[1]) * (z3 - w) ** lam[0] * (z1 - z2) ** lam[1]


# ====================================================================================== K1
ZETA1 = [(Fr(3, 2), Fr(5, 4), Fr(2, 3)), (Fr(2), Fr(5, 3), Fr(-1, 2)), (Fr(7, 5), Fr(-3, 4), Fr(1, 6)), (Fr(-4, 3), Fr(2, 9), Fr(5, 7)),
         (Fr(1, 3), Fr(-2, 7), Fr(5, 11)), (Fr(5, 3), Fr(-1, 2), Fr(3, 8)), (Fr(-3, 5), Fr(13, 7), Fr(2)), (Fr(4, 5), Fr(-5, 6), Fr(3, 4))]
THETA1 = Fr(2, 7)


def onesite_rows(lam, zeta, theta=THETA1, verbose=False):
    """for one zeta: list of dict(nu index, weight, monos, Qraw(omega-independent pieces)) over exactly diagonalised GT states"""
    R = one(lam)
    C = fit_pf(R, zeta, theta)
    assert C['res2'] == 0 and C['res3'] == 0, (C['res2'], C['res3'])
    ds = diag_states(R, C)
    xc = x_cov(R, zeta)
    out = []
    for a, (l2, l3) in sorted(ds.items()):
        sol = solve_dual(R, zeta, theta, l2, l3, R.wt[a])
        degs = tuple(s[0] for s in sol)
        if any(s[0] is None or s[1] != 1 for s in sol):
            out.append(dict(a=a, ok=False, degs=degs, nul=tuple(s[1] for s in sol)))
            continue
        J = [poly_jets(Fr(zeta[j]), sol[j][2], theta, lam[0] + 2) for j in range(3)]
        out.append(dict(a=a, ok=True, wt=R.wt[a], degs=degs, J=J, mono=mono_values(lam, J), Mv=None, x=xc[a], q=[s[2] for s in sol],
                        lam2=l2, lam3=l3))
    for r in out:
        if r['ok']:
            r['Mv'] = r['mono'][Mpair(lam)]
    return R, out, len(ds)


def fit_onesite(lam, w, zetas=None, theta=THETA1, cache=None):
    """exact fit of  c1 * Qnum * M = sum_kl K_kl m_kl  at fixed omega w, over all zeta points and diagonalised states; returns dict"""
    zetas = zetas or ZETA1
    cache = cache if cache is not None else {}
    prs = pairs_of(lam)
    A, y, tags = [], [], []
    for zi, z in enumerate(zetas):
        if zi not in cache:
            cache[zi] = onesite_rows(lam, z, theta)
        R, rows, nds = cache[zi]
        om = omega_inf(R, z, w)
        c1 = c1_scalar(lam, z, w)
        for r in rows:
            if not r['ok']:
                continue
            a = r['a']
            Qn = om[a] / r['x']
            A.append([r['mono'][kl] for kl in prs])
            y.append(c1 * Qn * r['Mv'])
            tags.append((zi, a))
    x, rk, cons = exact_solve(A, y)
    return dict(pairs=prs, K=x, rank=rk, consistent=cons, nrows=len(A))


def ckl_closed(lam, k, l):
    """closed form of the one-site coefficient c_kl(Lambda) (theta_kl(a) = c_kl a^(k+l-lam2-1)); lam = (lam1, lam2, 0)"""
    l1, l2 = lam[0], lam[1]
    return Fr((-1) ** (k + l - l2 - 1) * comb(l1 - l2, l - l2 - 1) * comb(l1 + 1, k) * (l - k), (l2 + 1) * comb(l, l2 + 1))


def stage_K1(lams=None):
    hdr('K1: one-site kernel, exact')
    lams = lams or [(1, 0, 0), (1, 1, 0), (2, 0, 0), (2, 1, 0), (2, 2, 0), (3, 1, 0), (3, 2, 0), (4, 1, 0), (4, 2, 0)]
    res = {}
    for lam in lams:
        t0 = time.time()
        cache = {}
        ws = [Fr(5, 2), Fr(7, 4), Fr(-2, 5), Fr(3)]
        out = {}
        for w in ws:
            f = fit_onesite(lam, w, cache=cache)
            out[w] = f
        R, rows, nds = cache[0]
        print('Lambda', lam, 'dim', R.d, 'diagonal GT states', nds, 'with unique dual solution', sum(1 for r in rows if r['ok']),
              ' time %.1fs' % (time.time() - t0))
        for w, f in out.items():
            print('  omega', w, 'rows', f['nrows'], 'unknowns', len(f['pairs']), 'rank', f['rank'], 'consistent', f['consistent'])
        # theta_kl(a) = K_kl / omega^(2 l1 - l2); coefficient c_kl = theta / a^e_kl
        coef = {}
        for w, f in out.items():
            a = 1 / w
            e = lambda kl: kl[0] + kl[1] - lam[1] - 1
            th = [f['K'][i] / w ** (2 * lam[0] - lam[1]) for i in range(len(f['pairs']))]
            coef[w] = [th[i] / a ** e(kl) for i, kl in enumerate(f['pairs'])]
        same = all(coef[w] == coef[ws[0]] for w in ws)
        print('  c_kl = theta_kl / a^e_kl omega-independent:', same)
        print('  ', [(kl, str(c)) for kl, c in zip(out[ws[0]]['pairs'], coef[ws[0]])])
        cf = all(ckl_closed(lam, *kl) == c for kl, c in zip(out[ws[0]]['pairs'], coef[ws[0]]))
        print('  closed form c_kl = (-1)^(k+l-l2-1) C(l1-l2,l-l2-1) C(l1+1,k) (l-k) / ((l2+1) C(l,l2+1)) reproduces all:', cf)
        res[lam] = dict(pairs=out[ws[0]]['pairs'], c=coef[ws[0]], closed=cf)
    with open(os.path.join(HERE, 'general_kernel_K1.pkl'), 'wb') as f:
        pickle.dump(res, f)
    return res



# ====================================================================================== two-site machinery (rep-generic)
from gaudin_chain import Chain, Operators, eigenstates, eigenfunctions
import twovac_fit as tf

DPS = 80
ZETAS2 = {'P': (Fr(3, 2), Fr(5, 4), Fr(2, 3)), 'Pp': (Fr(2), Fr(5, 3), Fr(-1, 2)), 'X1': (Fr(7, 5), Fr(-3, 4), Fr(1, 6)),
          'X2': (Fr(-4, 3), Fr(2, 9), Fr(5, 7)), 'X3': (Fr(1, 3), Fr(-2, 7), Fr(5, 11))}
TH0 = (Fr(2, 7), Fr(-3, 5))


def thkey(th):
    return '%d.%d_%d.%d' % (th[0].numerator, th[0].denominator, th[1].numerator, th[1].denominator)


def pt2_path(lam, tag, th, dps=DPS):
    return os.path.join(HERE, 'general_kernel_pt_%d%d%d_%s_%s_%d.pkl' % (lam[0], lam[1], lam[2], tag, thkey(th), dps))


def build_point2(lam, zeta, th=TH0, dps=DPS):
    """two-site states of Lambda at (zeta, theta): sparse eigenvectors, dual Q polynomials (monic), plain picklable record"""
    mp.mp.dps = dps
    Q.TOL = mp.mpf(10) ** -50
    zeta = tuple(Fr(z) for z in zeta)
    th = tuple(Fr(t) for t in th)
    t0 = time.time()
    ch = Chain(lam=lam, theta=th)
    g = [[zeta[i] if i == j else Fr(0) for j in range(3)] for i in range(3)]
    ops = Operators(ch, g)
    thf = set(th)
    pts = [Fr(3, 7), Fr(-2, 9), Fr(5, 3), Fr(7, 2), Fr(-4, 5), Fr(11, 4), Fr(-7, 3), Fr(13, 5), Fr(1, 3)]
    if any(u in thf for u in pts):
        pts = [Fr(17, 13), Fr(-19, 11), Fr(23, 17), Fr(29, 7), Fr(-31, 13), Fr(37, 19), Fr(-41, 17), Fr(43, 23), Fr(47, 29)]
    C = ops.fit(pts)
    assert C['res2'] == 0 and C['res3'] == 0, (C['res2'], C['res3'])
    t1 = time.time()
    u0, u1 = (Fr(3, 7), Fr(-2, 9)) if Fr(3, 7) not in thf and Fr(-2, 9) not in thf else (Fr(17, 13), Fr(-19, 11))
    st = eigenstates(ch, ops, C, u0=u0, u1=u1, dps=dps)
    for s in st:
        mxim = max([abs(mp.im(x)) for x in s['vec']] + [abs(mp.im(s['ev']))])
        if mxim > mp.mpf(10) ** (-dps // 2):
            raise ValueError('complex spectrum (imag %s)' % mp.nstr(mxim, 3))
        s['ev'] = mp.re(s['ev'])
        s['vec'] = [mp.re(x) for x in s['vec']]
    res = eigenfunctions(ch, C, st)
    t2 = time.time()
    thm = [tf.mpf(t) for t in th]
    ABS = Fr(sum(lam))
    S1 = Q.pf_S1(sum(tf.mpf(z) for z in zeta), ABS, thm)
    out = []
    for s in st:
        S2 = Q.pf_S2(s['lam2'], thm)
        S3 = Q.pf_S3(s['lam3'], thm)
        dual = []
        for j in range(3):
            d, nul, sm, q = Q.solve_poly_mp('D', S1, S2, S3, tf.mpf(zeta[j]), thm, dmax=2 * lam[0] + 3)
            dual.append((d, nul, q))
        mu = s['sector']
        ok = all(dual[j][0] == mu[j] and dual[j][1] == 1 for j in range(3))
        idx = s['idx']
        out.append(dict(sector=mu, ev=s['ev'], lam2=list(s['lam2']), lam3=list(s['lam3']), idx=list(idx), v=[s['vec'][i] for i in idx],
                        dual=[x[2] for x in dual], ok=ok, degs=tuple(x[0] for x in dual), nul=tuple(x[1] for x in dual)))
    t3 = time.time()
    return dict(lam=tuple(lam), zeta=zeta, theta=th, dps=dps, states=out, eigres=res, d1=ch.d1, hw=ch.hw_index,
                times=(t1 - t0, t2 - t1, t3 - t2))


def get_point2(lam, tag, zeta, th=TH0, dps=DPS):
    path = pt2_path(lam, tag, th, dps)
    if os.path.exists(path):
        with open(path, 'rb') as f:
            pt = pickle.load(f)
        mp.mp.dps = dps
        return pt
    pt = build_point2(lam, zeta, th, dps)
    with open(path + '.tmp', 'wb') as f:
        pickle.dump(pt, f)
    os.replace(path + '.tmp', path)
    print('  built point', lam, tag, 'theta', thkey(th), 'times (ops, eig, Q) %.0f/%.0f/%.0f s' % pt['times'], 'eigres', mp.nstr(pt['eigres'], 3),
          'bad states', sum(1 for s in pt['states'] if not s['ok']), flush=True)
    return pt


def _build_job(a):
    lam, tag, zeta, th = a
    get_point2(lam, tag, zeta, th)
    return tag, thkey(th)


def build_many2(jobs, workers=2):
    """jobs: list of (lam, tag, zeta, th); builds the missing points, `workers` processes at a time"""
    todo = [j for j in jobs if not os.path.exists(pt2_path(j[0], j[1], j[3]))]
    if not todo:
        return
    if workers <= 1 or len(todo) == 1:
        for j in todo:
            _build_job(j)
        return
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=workers) as ex:
        list(ex.map(_build_job, todo))


def poly_jets_mp(kappa, q, x, n):
    cur = list(q)
    out = [Q.peval(cur, x)]
    for _ in range(n):
        cur = Q.padd(Q.pder(cur), Q.pscale(cur, kappa))
        out.append(Q.peval(cur, x))
    return out


def site_monos(lam, pt):
    """for every state: [site] -> (list of monomials P3 p_kl over pairs_of(lam), M). Exponential factors dropped (common to a site)."""
    mp.mp.dps = pt['dps']
    z = [tf.mpf(x) for x in pt['zeta']]
    th = [tf.mpf(t) for t in pt['theta']]
    prs = pairs_of(lam)
    mpair = Mpair(lam)
    out = []
    for st in pt['states']:
        per = []
        for s in range(2):
            J = [poly_jets_mp(z[j], st['dual'][j], th[s], lam[0] + 2) for j in range(3)]
            mono = mono_values(lam, J)
            per.append(([mono[kl] for kl in prs], mono[mpair]))
        out.append(per)
    return out


def overlap2(pt, c1, c2):
    """<c1 (x) c2 | psi> for all states (sparse eigenvectors), covectors one-site (mp lists)"""
    d = pt['d1']
    out = []
    for st in pt['states']:
        tot = mp.mpf(0)
        for i, v in zip(st['idx'], st['v']):
            a, b = divmod(i, d)
            tot += c1[a] * c2[b] * v
        out.append(tot)
    return out


def cfac(lam, zeta, w1, w1p):
    z1, z2, z3 = [Fr(x) for x in zeta]
    f = lambda w: (z2 - Fr(w)) ** (lam[0] - lam[1]) * (z3 - Fr(w)) ** lam[0]
    return tf.mpf(f(w1) * f(w1p) * (z1 - z2) ** (2 * lam[1]))


class Fit2:
    """design matrix + QR of one (lam, zeta-points, theta) data set; kernels for many (w1, w1p) by re-using the factorisation"""
    def __init__(self, lam, pts):
        self.lam = tuple(lam)
        self.pts = pts
        self.prs = pairs_of(lam)
        self.sm = [site_monos(lam, pt) for pt in pts]
        self.index = [(p, k) for p, pt in enumerate(pts) for k in range(len(pt['states']))]
        self.F = []
        for p, k in self.index:
            f1, _ = self.sm[p][k][0]
            f2, _ = self.sm[p][k][1]
            self.F.append([a * b for a in f1 for b in f2])
        self.Mfac = [self.sm[p][k][0][1] * self.sm[p][k][1][1] for p, k in self.index]
        self.rn = [mp.sqrt(tf.dot(r, r)) for r in self.F]
        cols = [[self.F[i][j] / self.rn[i] for i in range(len(self.F))] for j in range(len(self.F[0]))]
        t0 = time.time()
        self.qr = tf.PivQR(cols)
        self.rank = self.qr.rank()
        self.nunk = len(cols)
        self.tqr = time.time() - t0
        R = one(lam)
        self.R = R
        self.x0 = [[tf.mpf(x) for x in x_cov(R, pt['zeta'])] for pt in pts]
        self._omc = {}

    def omega_cov(self, p, w):
        key = (p, Fr(w))
        if key not in self._omc:
            self._omc[key] = [tf.mpf(x) for x in omega_inf(self.R, self.pts[p]['zeta'], w)]
        return self._omc[key]

    def kernel(self, w1, w1p):
        """returns (K nested list n x n of mp numbers, relative residual, rank); K = fit of c * Qnum * M1 * M2"""
        ov = []
        for p, pt in enumerate(self.pts):
            ov.append((overlap2(pt, self.omega_cov(p, w1), self.omega_cov(p, w1p)), overlap2(pt, self.x0[p], self.x0[p])))
        y = []
        for i, (p, k) in enumerate(self.index):
            phi, x0 = ov[p]
            c = cfac(self.lam, self.pts[p]['zeta'], w1, w1p)
            y.append(c * (phi[k] / x0[k]) * self.Mfac[i] / self.rn[i])
        x, res, rk = self.qr.lstsq(y)
        n = len(self.prs)
        K = [[x[a * n + b] for b in range(n)] for a in range(n)]
        return K, res, rk


def theta_from_K(lam, K, w1, w1p):
    """canonical kernel Theta = K / (w1 w1p)^(2 l1 - l2), as nested list of mp numbers"""
    f = (tf.mpf(Fr(w1)) * tf.mpf(Fr(w1p))) ** (2 * lam[0] - lam[1])
    return [[x / f for x in row] for row in K]


def rat_mat(K, maxden=10 ** 14):
    out = []
    for row in K:
        r = []
        for x in row:
            f = tf.rat_recon(x, maxden=maxden)
            r.append(f)
        out.append(r)
    return out



# ====================================================================================== kernel on grids, polynomial reconstruction
W0 = Fr(5, 2)
WPS_ALL = [Fr(7, 4), Fr(3), Fr(-5, 3), Fr(11, 3), Fr(-7, 2), Fr(13, 5), Fr(9, 7), Fr(-2, 5), Fr(17, 6), Fr(-11, 4), Fr(11, 5), Fr(19, 7),
           Fr(-13, 5), Fr(7, 3), Fr(23, 9), Fr(-3, 2), Fr(10, 3)]
TH12_ALL = [Fr(31, 35), Fr(29, 45), Fr(5, 6), Fr(37, 35), Fr(22, 15), Fr(7, 9), Fr(13, 8), Fr(-9, 7), Fr(-17, 10), Fr(11, 6), Fr(3, 4), Fr(-5, 3),
            Fr(19, 13), Fr(-7, 11), Fr(9, 5), Fr(-13, 9), Fr(25, 14), Fr(-19, 10), Fr(15, 11), Fr(-11, 7)]


def theta_pair(th12):
    """theta = (theta1, theta2) with theta1 = 2/7 and theta12 = theta1 - theta2 (base: (2/7, -3/5), theta12 = 31/35)"""
    t1 = Fr(2, 7)
    return (t1, t1 - Fr(th12))


def grid_for(lam, nwp, nth, w0=W0, tags=('P',), workers=2, extra_w=(), th_list=None):
    """Theta on the grid (w1 = w0 fixed) x (w1p in WPS) x (th12 in TH12); zeta points `tags` are stacked (overdetermined if more than one);
    returns dict[(wp, th12)] = (Theta nested mp list, relative residual, rank) and diagnostics"""
    wps = WPS_ALL[:nwp]
    ths = th_list or TH12_ALL[:nth]
    jobs = [(lam, tg, ZETAS2[tg], theta_pair(t)) for t in ths for tg in tags]
    build_many2(jobs, workers)
    out = {}
    diag = {}
    for t in ths:
        pts = [get_point2(lam, tg, ZETAS2[tg], theta_pair(t)) for tg in tags]
        F = Fit2(lam, pts)
        diag[t] = dict(rank=F.rank, nunk=F.nunk, tqr=F.tqr, rows=len(F.index), bad=sum(1 for pt in pts for s in pt['states'] if not s['ok']))
        for wp in wps + list(extra_w):
            K, res, rk = F.kernel(w0, wp)
            out[(wp, t)] = (theta_from_K(lam, K, w0, wp), res, rk)
    return out, diag


def vander_pinv(xs, D):
    A = mp.matrix(len(xs), D + 1)
    for r, x in enumerate(xs):
        xm = tf.mpf(x)
        for p in range(D + 1):
            A[r, p] = xm ** p
    I = mp.eye(len(xs))
    P = mp.matrix(D + 1, len(xs))
    for c in range(len(xs)):
        sol, _ = mp.qr_solve(A, I[:, c])
        for p in range(D + 1):
            P[p, c] = sol[p]
    return A, P


def reconstruct(lam, grid, nwp, nth, Dx, Dy, w0=W0, th_list=None):
    """per entry: polynomial in (a' = 1/wp, t = 1/th12) of bidegree (Dx, Dy) at a = 1/w0 (tensor-product least squares with shared pseudo-inverses);
    returns dict[(i,j)] -> C[q][n] (Fractions or None) and the worst relative grid residual"""
    wps = WPS_ALL[:nwp]
    ths = th_list or TH12_ALL[:nth]
    n = len(pairs_of(lam))
    Ax, Px = vander_pinv([1 / w for w in wps], Dx)
    Ay, Py = vander_pinv([1 / t for t in ths], Dy)
    res = {}
    worst = mp.mpf(0)
    for i in range(n):
        for j in range(n):
            G = mp.matrix(len(wps), len(ths))
            for r, wp in enumerate(wps):
                for c, t in enumerate(ths):
                    G[r, c] = grid[(wp, t)][0][i][j]
            C = Px * G * Py.T
            Gp = Ax * C * Ay.T
            mx = max(abs(v) for v in G) or 1
            worst = max(worst, max(abs(Gp[r, c] - G[r, c]) for r in range(len(wps)) for c in range(len(ths))) / mx)
            res[(i, j)] = [[tf.rat_recon(C[q, m], maxden=10 ** 15, tol=mp.mpf(10) ** -30) for m in range(Dy + 1)] for q in range(Dx + 1)]
    return res, worst


# ====================================================================================== analysis of reconstructed kernels
asym, apsym, tsym = sympy.symbols('a ap t')


def e_of(lam, kl):
    return kl[0] + kl[1] - lam[1] - 1


def polys_from_rep(lam, rep, w0=W0):
    """T[(i,j)] = {(n, p, q): Fraction} with Theta_ij = sum T a^p ap^q t^n, deduced from the fixed-a fit assuming the homogeneity
    degree p + q = e_i + e_j + n (a-power p = e_i+e_j+n-q must be >= 0); returns (T, report)"""
    prs = pairs_of(lam)
    a0 = 1 / Fr(w0)
    T = {}
    nonhom = 0
    neg = 0
    for (i, j), C in rep['res'].items():
        ei, ej = e_of(lam, prs[i]), e_of(lam, prs[j])
        ent = {}
        for q, row in enumerate(C):
            for n, c in enumerate(row):
                if c is None:
                    nonhom += 1
                    continue
                if c == 0:
                    continue
                p = ei + ej + n - q
                if p < 0:
                    neg += 1
                    continue
                ent[(n, p, q)] = c / a0 ** p
        T[(i, j)] = ent
    return T, dict(nonrational_or_neg=nonhom + neg)


def sym_entry(ent):
    return sum(sympy.Rational(c.numerator, c.denominator) * tsym ** n * asym ** p * apsym ** q for (n, p, q), c in ent.items())


def eval_T(ent, a, ap, t):
    return sum(c * t ** n * a ** p * ap ** q for (n, p, q), c in ent.items())


def check_homogeneity_extra(lam, T, ths, w_alt=Fr(-2, 5), wps=(Fr(7, 4), Fr(3)), tag='P'):
    """independent check at a second value a = 1/w_alt (not used in the reconstruction): fit values vs polynomial prediction"""
    worst = mp.mpf(0)
    for t12 in ths:
        pt = get_point2(lam, tag, ZETAS2[tag], theta_pair(t12))
        F = Fit2(lam, [pt])
        for wp in wps:
            K, res, rk = F.kernel(w_alt, wp)
            Th = theta_from_K(lam, K, w_alt, wp)
            n = len(Th)
            for i in range(n):
                for j in range(n):
                    pred = eval_T(T[(i, j)], 1 / w_alt, 1 / wp, 1 / t12)
                    d = abs(Th[i][j] - tf.mpf(pred)) / max(abs(Th[i][j]), mp.mpf(10) ** -10)
                    worst = max(worst, d)
    return worst



def zeta_check(lam, tags=('P', 'Pp', 'X1'), th12=Fr(31, 35), w=(Fr(5, 2), Fr(7, 4)), workers=1):
    """stacked fit over several zeta points (overdetermined: tests that the formula holds and the kernel is zeta-independent), vs each single-zeta fit"""
    th = theta_pair(th12)
    build_many2([(lam, tg, ZETAS2[tg], th) for tg in tags], workers)
    pts = [get_point2(lam, tg, ZETAS2[tg], th) for tg in tags]
    Fs = Fit2(lam, pts)
    K, res, rk = Fs.kernel(*w)
    Th = theta_from_K(lam, K, *w)
    out = dict(stacked_resid=res, stacked_rank=rk, nunk=Fs.nunk, rows=len(Fs.index), per_zeta={})
    for tg, pt in zip(tags, pts):
        F1 = Fit2(lam, [pt])
        K1, r1, rk1 = F1.kernel(*w)
        T1 = theta_from_K(lam, K1, *w)
        d = max(abs(T1[i][j] - Th[i][j]) for i in range(len(Th)) for j in range(len(Th)))
        out['per_zeta'][tg] = (r1, rk1, d, len(F1.index))
    return out, Th



# ====================================================================================== top-state sum rule (K3c)
_U = sympy.Symbol('u')


def top_monos_exact(lam, z2, z3, thetas):
    """exact top-state data (rationals): P3 = e^{z3 u}, P2 = e^{z2 u} p2, p2 = sum_k (z3-z2)^(-k) (Pi^l2)^(k) (monic, degree l2 * L),
    returns [site] -> (dict kl -> P3 p_kl exp-free, M)"""
    l1, l2 = lam[0], lam[1]
    d = sympy.Rational(z3) - sympy.Rational(z2)
    Pi = sympy.Integer(1)
    for th in thetas:
        Pi *= (_U - sympy.Rational(th)) ** l2
    p2 = sympy.Integer(0)
    cur = Pi
    k = 0
    while cur != 0:
        p2 += cur / d ** k
        cur = sympy.diff(cur, _U)
        k += 1
    out = []
    for th in thetas:
        th = sympy.Rational(th)
        r = []
        cur = p2
        for k in range(l1 + 3):
            r.append(cur.subs(_U, th))
            cur = sympy.diff(cur, _U) + sympy.Rational(z2) * cur
        z3r = sympy.Rational(z3)
        mono = {}
        for (k, l) in pairs_of(lam):
            mono[(k, l)] = r[k] * z3r ** l - z3r ** k * r[l]
        out.append((mono, mono[Mpair(lam)]))
    return out


def top_sum_rule_check(lam, T, trials=4, seed=7):
    """checks  sum_ij Theta_ij m_i(th1) m_j(th2) = M(th1) M(th2) R(a) R(a'),  R(a) = (1 - z2 a)^(l1-l2) (1 - z3 a)^l1  (exact rationals)"""
    import random
    rnd = random.Random(seed)
    prs = pairs_of(lam)
    n = len(prs)
    worst = 0
    for _ in range(trials):
        z2, z3 = Fr(rnd.randint(-9, 9), rnd.randint(2, 7)), Fr(rnd.randint(-9, 9), rnd.randint(2, 7))
        if z2 == z3:
            continue
        th = (Fr(rnd.randint(-9, 9), rnd.randint(2, 7)), Fr(rnd.randint(-9, 9), rnd.randint(2, 7)))
        if th[0] == th[1]:
            continue
        a, ap = Fr(rnd.randint(-9, 9), rnd.randint(2, 7)), Fr(rnd.randint(-9, 9), rnd.randint(2, 7))
        t = 1 / (th[0] - th[1])
        (m1, M1), (m2, M2) = top_monos_exact(lam, z2, z3, th)
        lhs = sympy.Integer(0)
        for i in range(n):
            for j in range(n):
                th_ij = eval_T(T[(i, j)], a, ap, t)
                lhs += sympy.Rational(th_ij.numerator, th_ij.denominator) * m1[prs[i]] * m2[prs[j]]
        R = lambda x: (1 - sympy.Rational(z2) * sympy.Rational(x)) ** (lam[0] - lam[1]) * (1 - sympy.Rational(z3) * sympy.Rational(x)) ** lam[0]
        rhs = M1 * M2 * R(a) * R(ap)
        worst = max(worst, abs(lhs - rhs) / (abs(rhs) + 1))
    return worst



def block_check(l1, l2s, ws=((Fr(5, 2), Fr(7, 4)), (Fr(-2, 5), Fr(3))), th12s=(Fr(31, 35), Fr(29, 45)), tags0=('P', 'Pp'), workers=2):
    """direct fits of [l1,0,0] (stacked zeta points) and [l1,l2,0] at a few (w,w',theta12); tests  Theta^[l1,l2,0]_ij = N (a a')^(-l2) Theta^[l1,0,0]_ij
    (non-forced i,j), N = (l2! (l1-l2)! / l1!)^2.  Returns worst relative deviation per l2."""
    from math import factorial as fa
    lam0 = (l1, 0, 0)
    out = {}
    for t12 in th12s:
        th = theta_pair(t12)
        build_many2([(lam0, tg, ZETAS2[tg], th) for tg in tags0] + [((l1, l2, 0), 'P', ZETAS2['P'], th) for l2 in l2s], workers)
        F0 = Fit2(lam0, [get_point2(lam0, tg, ZETAS2[tg], th) for tg in tags0])
        prs0 = pairs_of(lam0)
        for l2 in l2s:
            lam = (l1, l2, 0)
            F1 = Fit2(lam, [get_point2(lam, 'P', ZETAS2['P'], th)])
            prs = pairs_of(lam)
            for (w, wp) in ws:
                K0, r0, _ = F0.kernel(w, wp)
                T0 = theta_from_K(lam0, K0, w, wp)
                K1, r1, _ = F1.kernel(w, wp)
                T1 = theta_from_K(lam, K1, w, wp)
                N = (Fr(fa(l2) * fa(l1 - l2), fa(l1))) ** 2
                pref = tf.mpf(N) * (tf.mpf(Fr(w)) * tf.mpf(Fr(wp))) ** l2
                worst = mp.mpf(0)
                for i, p in enumerate(prs):
                    for j, q in enumerate(prs):
                        v0 = T0[prs0.index(p)][prs0.index(q)] * pref
                        worst = max(worst, abs(v0 - T1[i][j]) / max(abs(T1[i][j]), mp.mpf(10) ** -12))
                out[(l2, t12, w, wp)] = (worst, r0, r1)
    return out


# ---------------------------------------------------------------------------------- closed form of the [l,l,0] kernel (conjecture C3)
def poch_fr(z, n):
    """rising factorial (z)_n, n >= 0 integer, or n = -1: 1/(z-1)  (Gamma continuation)"""
    z = Fr(z)
    if n == -1:
        return 1 / (z - 1)
    r = Fr(1)
    for i in range(n):
        r *= z + i
    return r


def F_mm(m, mp_, k, kp):
    """F_{m m'}(k,k') = (k+m+1)_{m'-1} (k'+m'+1)_{m-1} (k k' + m' k + m k')  (Gamma continuation for m or m' = 0)"""
    if m == 0 and mp_ == 0:
        return Fr(1)
    if mp_ == 0:
        return poch_fr(kp, m)
    if m == 0:
        return poch_fr(k, mp_)
    return poch_fr(k + m + 1, mp_ - 1) * poch_fr(kp + mp_ + 1, m - 1) * (k * kp + mp_ * k + m * kp)


def theta_single_closed(L, k, kp, a, ap, t):
    """Theta_{k k'} of Lambda = [L,L,0] (pairs (k, L+1)); exact Fractions, evaluated at (a, a', t)"""
    a, ap, t = Fr(a), Fr(ap), Fr(t)
    R = Fr(0)
    for m in range(L - k + 1):
        for mp_ in range(L - kp + 1):
            R += comb(L - k, m) * comb(L - kp, mp_) * F_mm(m, mp_, k, kp) * (a * t) ** m * (-ap * t) ** mp_
    return (-1) ** (k + kp) * comb(L, k) * comb(L, kp) * a ** k * ap ** kp * R


def single_check(L, th12s=(Fr(31, 35),), ws=((Fr(5, 2), Fr(7, 4)), (Fr(-2, 5), Fr(3)), (Fr(13, 5), Fr(-7, 2))), tags=('P',), workers=1):
    lam = (L, L, 0)
    out = {}
    for t12 in th12s:
        th = theta_pair(t12)
        build_many2([(lam, tg, ZETAS2[tg], th) for tg in tags], workers)
        F = Fit2(lam, [get_point2(lam, tg, ZETAS2[tg], th) for tg in tags])
        for (w, wp) in ws:
            K, res, rk = F.kernel(w, wp)
            Th = theta_from_K(lam, K, w, wp)
            worst = mp.mpf(0)
            for k in range(L + 1):
                for kp in range(L + 1):
                    c = tf.mpf(theta_single_closed(L, k, kp, 1 / Fr(w), 1 / Fr(wp), 1 / t12))
                    worst = max(worst, abs(Th[k][kp] - c) / max(abs(c), mp.mpf(10) ** -12))
            out[(t12, w, wp)] = (worst, res, rk, F.nunk, len(F.index))
    return out


# ====================================================================================== ANALYZE: K0/K3 verdicts from the reconstructed polynomial kernels
def load_rep(lam):
    with open(os.path.join(HERE, 'general_kernel_rep_%d%d%d.pkl' % tuple(lam)), 'rb') as f:
        rep = pickle.load(f)
    T, rpt = polys_from_rep(lam, rep)
    return T, rep


def theta_vec_sym(lam, av):
    prs = pairs_of(lam)
    out = []
    for (k, l) in prs:
        c = ckl_closed(lam, k, l)
        out.append(sympy.Rational(c.numerator, c.denominator) * av ** e_of(lam, (k, l)))
    return out


def stage_ANALYZE(lams):
    from math import factorial as fa
    hdr('ANALYZE: polynomial laws, exchange rule, t = 0 factorisation, expansion ranks, top-state sum rule, block theorem')
    Ts = {}
    for lam in lams:
        T, rep = load_rep(lam)
        Ts[lam] = T
        prs = pairs_of(lam)
        n = len(prs)
        Dt = max(k[0] for e in T.values() for k in e)
        Da = max(k[1] for e in T.values() for k in e)
        Dq = max(k[2] for e in T.values() for k in e)
        ex = all(sympy.expand(sym_entry(T[(i, j)]) - sym_entry(T[(j, i)]).subs({asym: apsym, apsym: asym, tsym: -tsym}, simultaneous=True)) == 0
                 for i in range(n) for j in range(n))
        th_a, th_b = theta_vec_sym(lam, asym), theta_vec_sym(lam, apsym)
        fact0 = all(sympy.expand(sum(c * asym ** p * apsym ** q for (nn, p, q), c in [(k, sympy.Rational(v.numerator, v.denominator)) for k, v in T[(i, j)].items() if k[0] == 0]) - th_a[i] * th_b[j]) == 0
                    for i in range(n) for j in range(n))
        av, apv, tv = sympy.Rational(2, 5), sympy.Rational(4, 7), sympy.Rational(3, 11)
        ranks = []
        for k in range(Dt + 1):
            Mk = sympy.Matrix(n, n, lambda i, j: sum(c * av ** p * apv ** q for (nn, p, q), c in T[(i, j)].items() if nn == k))
            ranks.append(Mk.rank())
        Mfull = sympy.Matrix(n, n, lambda i, j: sym_entry(T[(i, j)]).subs({asym: av, apsym: apv, tsym: tv}))
        sr = top_sum_rule_check(lam, T, trials=3)
        print('Lambda %s: pairs %d, a-degree %d, a\'-degree %d (2l1-l2 = %d), t-degree %d, exchange rule %s, Theta(t=0) = theta(a) x theta(a\') exactly %s, rank Theta(generic) %d, ranks of Theta_n %s, top-state sum rule worst dev %s'
              % (lam, n, Da, Dq, 2 * lam[0] - lam[1], Dt, ex, fact0, Mfull.rank(), ranks, sr), flush=True)
    # block theorem
    for l1 in sorted(set(l[0] for l in lams)):
        base = (l1, 0, 0)
        if base not in Ts:
            continue
        for lam in lams:
            if lam[0] != l1 or lam[1] == 0:
                continue
            N = sympy.Rational(fa(lam[1]) * fa(l1 - lam[1]), fa(l1)) ** 2
            prs = pairs_of(lam)
            prs0 = pairs_of(base)
            ok = True
            for i, p in enumerate(prs):
                for j, q in enumerate(prs):
                    e1 = sym_entry(Ts[lam][(i, j)])
                    e0 = sym_entry(Ts[base][(prs0.index(p), prs0.index(q))])
                    if sympy.expand(e1 - N * (asym * apsym) ** (-lam[1]) * e0) != 0:
                        ok = False
            print('block theorem: Theta^%s == N (a a\')^(-l2) Theta^%s restricted to non-forced pairs, N = %s : %s' % (lam, base, N, ok), flush=True)



def corner_check(lam, th12=Fr(31, 35), ws=((Fr(5, 2), Fr(7, 4)), (Fr(-2, 5), Fr(3))), tag='P', workers=1):
    """corner block of [l1,l2,0]: the pairs (k, l1+1), k <= l1, must equal  (N_{l1,l2} / N_{l1,l1}) (a a')^(l1-l2) Theta^{[l1,l1,0]}_{kk'} (closed form C3),
    N_{l1,l2} = (l2! (l1-l2)! / l1!)^2"""
    from math import factorial as fa
    l1, l2 = lam[0], lam[1]
    th = theta_pair(th12)
    build_many2([(lam, tag, ZETAS2[tag], th)], workers)
    F = Fit2(lam, [get_point2(lam, tag, ZETAS2[tag], th)])
    prs = pairs_of(lam)
    cidx = [i for i, p in enumerate(prs) if p[1] == l1 + 1]
    out = {}
    for (w, wp) in ws:
        K, res, rk = F.kernel(w, wp)
        Th = theta_from_K(lam, K, w, wp)
        a, ap = 1 / Fr(w), 1 / Fr(wp)
        ratio = Fr(fa(l2) * fa(l1 - l2), fa(l1)) ** 2 / Fr(fa(l1) * fa(0), fa(l1)) ** 2
        worst = mp.mpf(0)
        for i in cidx:
            for j in cidx:
                k, kp = prs[i][0], prs[j][0]
                pred = tf.mpf(ratio * (a * ap) ** (l1 - l2) * theta_single_closed(l1, k, kp, a, ap, 1 / Fr(th12)))
                worst = max(worst, abs(Th[i][j] - pred) / max(abs(pred), mp.mpf(10) ** -12))
        out[(w, wp)] = (worst, res, rk, F.nunk, len(F.index))
    return out


def theta1_closed(lam, p, q, a, ap, t=1):
    """first-order coefficient (times t) of Theta_{p q}, p = (k,l), q = (k',l'), exact Fractions (conjecture C4):
    t theta_p(a) theta_q(a') [ alpha a + beta a' + gamma a^2/a' + eps a'^2/a ],  n = l1+1, rho(p) = (n-k)(n-l), sigma(p) = k+l-1, N0 = n(n-1)"""
    n = lam[0] + 1
    N0 = n * (n - 1)
    rho = lambda r: (n - r[0]) * (n - r[1])
    sig = lambda r: r[0] + r[1] - 1
    f2 = lambda r: 2 * (n - 1) - sig(r)
    al = Fr(-3, n * n) * rho(p) * (rho(q) - N0) - Fr(2, n) * rho(p) * sig(q) + Fr(1, n) * f2(p) * (rho(q) - N0) + f2(p) * sig(q)
    be = -(Fr(-3, n * n) * rho(q) * (rho(p) - N0) - Fr(2, n) * rho(q) * sig(p) + Fr(1, n) * f2(q) * (rho(p) - N0) + f2(q) * sig(p))
    ga = Fr(q[0] * q[1] * rho(p), n * n)
    ep = -Fr(p[0] * p[1] * rho(q), n * n)
    a, ap = Fr(a), Fr(ap)
    cp, cq = ckl_closed(lam, *p), ckl_closed(lam, *q)
    th = cp * a ** e_of(lam, p) * cq * ap ** e_of(lam, q)
    return th * (al * a + be * ap + ga * a * a / ap + ep * ap * ap / a) * Fr(t)

def run_rep(lam, nwp, nth, Dx, Dy, workers=2, save=True, tags=('P',)):
    t0 = time.time()
    grid, diag = grid_for(lam, nwp, nth, workers=workers, tags=tags)
    print('Lambda', lam, 'grid', nwp, 'x', nth, 'time %.0fs' % (time.time() - t0), ' rows/unknowns/rank at theta0:', list(diag.values())[0])
    print('  max Q-side relative residual', mp.nstr(max(v[1] for v in grid.values()), 3), ' ranks', sorted(set(v[2] for v in grid.values())))
    res, worst = reconstruct(lam, grid, nwp, nth, Dx, Dy)
    nbad = sum(1 for v in res.values() for row in v for c in row if c is None)
    print('  polynomial (ap, t) fit of bidegree (%d,%d): worst relative residual %s, non-rational coefficients %d' % (Dx, Dy, mp.nstr(worst, 3), nbad))
    if save:
        with open(os.path.join(HERE, 'general_kernel_rep_%d%d%d.pkl' % tuple(lam)), 'wb') as f:
            pickle.dump(dict(lam=lam, nwp=nwp, nth=nth, Dx=Dx, Dy=Dy, res=res, worst=worst, grid={k: ([[str(x) for x in r] for r in v[0]], str(v[1]), v[2]) for k, v in grid.items()}), f)
    return grid, res, worst

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--stages', default='K1')
    ap.add_argument('--lams', default='')
    ap.add_argument('--lam', default='')
    ap.add_argument('--nwp', type=int, default=6)
    ap.add_argument('--nth', type=int, default=6)
    ap.add_argument('--dx', type=int, default=3)
    ap.add_argument('--dy', type=int, default=4)
    ap.add_argument('--workers', type=int, default=2)
    ap.add_argument('--tags', default='P')
    ap.add_argument('--l2s', default='1')
    ap.add_argument('--th12s', default='31/35')
    args = ap.parse_args()
    st = args.stages.split(',')
    if 'K1' in st:
        lams = [tuple(int(c) for c in s) for s in args.lams.split(',')] if args.lams else None
        stage_K1(lams)
    if 'REP' in st:
        run_rep(tuple(int(c) for c in args.lam), args.nwp, args.nth, args.dx, args.dy, args.workers, tags=tuple(args.tags.split(',')))
    if 'ZC' in st:
        for tag in args.lam.split(','):
            lam = tuple(int(c) for c in tag)
            out, Th = zeta_check(lam, workers=args.workers)
            print('Lambda', lam, 'stacked rows', out['rows'], 'unknowns', out['nunk'], 'rank', out['stacked_rank'], 'resid', mp.nstr(out['stacked_resid'], 3))
            for tg, v in out['per_zeta'].items():
                print('   zeta', tg, 'single-zeta resid', mp.nstr(v[0], 3), 'rank', v[1], 'rows', v[3], ' max |Theta_single - Theta_stacked|', mp.nstr(v[2], 3))
    if 'BLK' in st:
        l1 = int(args.lam)
        l2s = [int(c) for c in args.l2s.split(',')]
        res = block_check(l1, l2s, workers=args.workers)
        for k, v in res.items():
            print('block check [%d,%d,0] vs [%d,0,0]: theta12 %s w=(%s,%s)  max rel dev %s  (fit resid %s, %s)' % (l1, k[0], l1, k[1], k[2], k[3], mp.nstr(v[0], 3), mp.nstr(v[1], 2), mp.nstr(v[2], 2)))
    if 'SGL' in st:
        for L in [int(c) for c in args.lam.split(',')]:
            res = single_check(L, th12s=tuple(Fr(x) for x in args.th12s.split(',')), workers=args.workers)
            for k, v in res.items():
                print('closed form [%d,%d,0]: theta12 %s w=(%s,%s): max rel dev %s (fit resid %s, rank %d of %d unknowns, %d rows)' % (L, L, k[0], k[1], k[2], mp.nstr(v[0], 3), mp.nstr(v[1], 2), v[2], v[3], v[4]))
    if 'ANALYZE' in st:
        stage_ANALYZE([tuple(int(c) for c in tg) for tg in args.lam.split(',')])
    if 'CORNER' in st:
        for tg in args.lam.split(','):
            lam = tuple(int(c) for c in tg)
            res = corner_check(lam, workers=args.workers)
            for k, v in res.items():
                print('corner block of %s vs C3 closed form [%d,%d,0]: w=(%s,%s) max rel dev %s (fit resid %s, rank %d of %d unknowns, %d rows)' % (lam, lam[0], lam[0], k[0], k[1], mp.nstr(v[0], 3), mp.nstr(v[1], 2), v[2], v[3], v[4]))
