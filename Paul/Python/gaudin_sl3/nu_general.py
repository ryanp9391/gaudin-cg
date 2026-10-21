"""N1 (round 5, WP-B): the nu-formulas for general Lambda (B1) and the untwisting limit for [3,1,0], [3,2,0] (B2).
Part B1 (this block): symbolic Frobenius at theta_s of the dual oper  P''' - S1 P'' + S2 P' - S3 P = 0,
  S1 = |L|/x + sum s_k x^k, S2 = c2/x^2 + H/x + sum h_k x^k, S3 = a/x^2 + b/x + sum t_k x^k   (x = u - theta_s),
  exponents {0, lam2+1, lam1+2}, c2 = lam2 (lam1+1).  P = sum c_m x^m;  c_m = c_0 u_m + c_{lam2+1} v_m for m <= lam1+1.
  nu_kl = k! l! (u_k v_l - u_l v_k)/(lam2+1)!  (non-forced pairs k<l<=lam1+1, not both <= lam2), nu_{0,lam2+1} = 1.
"""
import sys, os, time, pickle
from fractions import Fraction as Fr
from math import factorial
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import sympy as sp
import mpmath as mp

Hs, As, Bs, Ht, At, Bt, SG1, E2, E3, DL = sp.symbols('H_s a_s b_s H_t a_t b_t sigma1 e2 e3 delta')
SYMS = (Hs, As, Bs, Ht, At, Bt, SG1, E2, E3, DL)


def pairs_of(lam):
    l1, l2 = lam[0], lam[1]
    return [(k, l) for k in range(l1 + 2) for l in range(k + 1, l1 + 2) if not (l <= l2)]


def rho(m, lam):
    L = lam[0] + lam[1]
    c2 = lam[1] * (lam[0] + 1)
    return m * (m - 1) * (m - 2) - L * m * (m - 1) + c2 * m


def reg_coeffs(lam, K):
    """Taylor coefficients (x^k, k<K) of the regular parts s_k, h_k, t_k at theta_s (other site at distance delta = theta_s-theta_t, constants sigma1,e2,e3)"""
    L = lam[0] + lam[1]
    c2 = lam[1] * (lam[0] + 1)
    d = lambda k, m: (-1) ** k * sp.binomial(m + k - 1, k) * DL ** (-m - k)     # coefficient of x^k in 1/(x+delta)^m
    s = [L * d(k, 1) + (SG1 if k == 0 else 0) for k in range(K)]
    h = [Ht * d(k, 1) + c2 * d(k, 2) + (E2 if k == 0 else 0) for k in range(K)]
    t = [Bt * d(k, 1) + At * d(k, 2) + (E3 if k == 0 else 0) for k in range(K)]
    return s, h, t


def resid_eq(lam, m, c, s, h, t):
    """coefficient of x^{m-3} of the ODE, EXCLUDING the c_m rho(m) term; c: dict index->expr (missing = 0)"""
    L = lam[0] + lam[1]
    tot = 0
    F = lambda n: c.get(n, 0) if n >= 0 else 0
    for k in range(len(s)):
        n = m - 1 - k
        if n >= 0:
            tot -= s[k] * n * (n - 1) * F(n)
        n = m - 2 - k
        if n >= 0:
            tot += h[k] * n * F(n)
        n = m - 3 - k
        if n >= 0:
            tot -= t[k] * F(n)
    tot += Hs * (m - 1) * F(m - 1) - As * F(m - 1) - Bs * F(m - 2)
    return tot


_CACHE = {}
def frobenius(lam):
    """returns dict(u=[u_0..u_{l1+1}], v=[...], nolog=[A, B1, B2], rho_res=...) as sympy expressions in SYMS.
    u: c_0 = 1, c_{l2+1} = 0 ;  v: c_0 = 0, c_{l2+1} = 1 ; both solve the recursion for m not in {l2+1, l1+2}."""
    lam = tuple(lam[:2])
    if lam in _CACHE:
        return _CACHE[lam]
    l1, l2 = lam
    mres1, mres2 = l2 + 1, l1 + 2
    K = mres2 + 2
    s, h, t = reg_coeffs(lam, K)
    out = {}
    for name, seed in (('u', {0: sp.Integer(1)}), ('v', {mres1: sp.Integer(1)})):
        c = dict(seed)
        for m in range(1, mres2 + 1):
            if m in (mres1, mres2):
                continue
            if m in seed:
                continue
            c[m] = sp.expand(-resid_eq(lam, m, c, s, h, t) / rho(m, lam))
        out[name] = c
    # no-log: equation at m = mres1 (only c_0 solution matters -> u), and at m = mres2 for u and v
    A = sp.expand(resid_eq(lam, mres1, out['u'], s, h, t))
    B1 = sp.expand(resid_eq(lam, mres2, out['u'], s, h, t))
    B2 = sp.expand(resid_eq(lam, mres2, out['v'], s, h, t))
    # v at m = mres1 is seed 1: the equation at mres1 for v must vanish identically (c_{m-1}=... =0 below)
    Av = sp.expand(resid_eq(lam, mres1, out['v'], s, h, t))
    assert Av == 0
    res = dict(u=[out['u'].get(m, sp.Integer(0)) for m in range(l1 + 2)], v=[out['v'].get(m, sp.Integer(0)) for m in range(l1 + 2)], nolog=[A, B1, B2])
    _CACHE[lam] = res
    return res


def nu_symbolic(lam1, lam2):
    """list of nu_kl over pairs_of((lam1,lam2)) as sympy expressions in (H_s,a_s,b_s,H_t,a_t,b_t,sigma1,e2,e3,delta); also the pair list and the no-log conditions"""
    lam = (lam1, lam2)
    fr = frobenius(lam)
    u, v = fr['u'], fr['v']
    prs = pairs_of(lam)
    nu = []
    for (k, l) in prs:
        nu.append(sp.expand(factorial(k) * factorial(l) * (u[k] * v[l] - u[l] * v[k]) / sp.Integer(factorial(lam2 + 1))))
    return prs, nu, fr['nolog']


_NUF = {}
def nu_numeric_funcs(lam):
    lam = tuple(lam[:2])
    if lam not in _NUF:
        prs, nu, nl = nu_symbolic(*lam)
        _NUF[lam] = (prs, [sp.lambdify(SYMS, e, 'mpmath') for e in nu], [sp.lambdify(SYMS, e, 'mpmath') for e in nl])
    return _NUF[lam]


def laurent_of(lam2, lam3, s):
    """(H_s, a_s, b_s) from the coefficient eigenvalue lists (gaudin_chain basis: lam2 = [e2,H1,c2,H2,c2], lam3 = [e3,b1,a1,0,b2,a2,0])"""
    if s == 0:
        return lam2[1], lam3[2], lam3[1]
    return lam2[3], lam3[5], lam3[4]


def nu_numeric(lam, laurent_s, laurent_t, zeta, theta_s, theta_t):
    """nu_kl at theta_s (mpmath) for laurent_* = (H, a, b); zeta = 3 numbers"""
    z = [mp.mpf(x.numerator) / x.denominator if isinstance(x, Fr) else mp.mpf(x) for x in zeta]
    e2 = z[0] * z[1] + z[0] * z[2] + z[1] * z[2]; e3 = z[0] * z[1] * z[2]; sg = sum(z)
    dl = _mp(theta_s) - _mp(theta_t)
    prs, fn, nl = nu_numeric_funcs(lam)
    args = (laurent_s[0], laurent_s[1], laurent_s[2], laurent_t[0], laurent_t[1], laurent_t[2], sg, e2, e3, dl)
    return [f(*args) for f in fn]


def nolog_numeric(lam, laurent_s, laurent_t, zeta, theta_s, theta_t):
    z = [mp.mpf(x.numerator) / x.denominator if isinstance(x, Fr) else mp.mpf(x) for x in zeta]
    e2 = z[0] * z[1] + z[0] * z[2] + z[1] * z[2]; e3 = z[0] * z[1] * z[2]; sg = sum(z)
    dl = _mp(theta_s) - _mp(theta_t)
    prs, fn, nl = nu_numeric_funcs(lam)
    args = (laurent_s[0], laurent_s[1], laurent_s[2], laurent_t[0], laurent_t[1], laurent_t[2], sg, e2, e3, dl)
    return [f(*args) for f in nl]


def _mp(x):
    return mp.mpf(x.numerator) / x.denominator if isinstance(x, Fr) else mp.mpf(x)


# ================================================================================================ B1 verification on actual states
def verify_B1(lam, tag='P', dps=80, th=None):
    import general_kernel as gk
    lam = tuple(lam)
    zeta = gk.ZETAS2[tag]
    th = th or gk.TH0
    pt = gk.get_point2(lam, tag, zeta, th, dps)
    mp.mp.dps = dps
    sm = gk.site_monos(lam, pt)
    prs = pairs_of(lam)
    worst = mp.mpf(0); worst_nl = [mp.mpf(0)] * 3; nbad = 0; nst = 0
    wsite = [mp.mpf(0), mp.mpf(0)]
    for st, per in zip(pt['states'], sm):
        L = [laurent_of(st['lam2'], st['lam3'], s) for s in (0, 1)]
        for s in (0, 1):
            mono, M = per[s]
            nu = nu_numeric(lam, L[s], L[1 - s], zeta, th[s], th[1 - s])
            for i, (m_, n_) in enumerate(zip(mono, nu)):
                d = abs(m_ / M - n_) / max(1, abs(n_))
                worst = max(worst, d); wsite[s] = max(wsite[s], d)
            nl = nolog_numeric(lam, L[s], L[1 - s], zeta, th[s], th[1 - s])
            for i in range(3):
                scale = max(1, abs(L[s][0]), abs(L[s][1]) ** 2, abs(L[s][2]))
                worst_nl[i] = max(worst_nl[i], abs(nl[i]) / scale)
        nst += 1
    return dict(nstates=nst, worst=worst, wsite=wsite, nolog=worst_nl, pairs=prs)

if __name__ == '__main__' and 'B1' in sys.argv:
    for lam in ((2, 1, 0), (3, 1, 0), (3, 2, 0)):
        t0 = time.time()
        r = verify_B1(lam)
        print(lam, 'states', r['nstates'], 'max |m_i/M - nu_i| = %s' % mp.nstr(r['worst'], 3), 'by site', [mp.nstr(x, 3) for x in r['wsite']], 'no-log residuals', [mp.nstr(x, 3) for x in r['nolog']], '%.0fs' % (time.time() - t0), flush=True)


# ================================================================================================ W4 exact data
import re
W4DIR = os.path.join(HERE, '..', '..', 'Mathematica', 'Gaudin', 'ref')


def _parse(s):
    s = re.sub(r'(-?\d+)/(\d+)', r'F(\1,\2)', s)
    s = s.replace('{', '[').replace('}', ']')
    return eval(s, {'F': Fr})


def load_w4(lam):
    """parse untwisted-exact-<lam>.txt -> dict(rep, dim, gtbasis, sectors{mu: dict(mult, hwbasis[list of {idx:val}], gram, branches[list])}, xvec)"""
    fn = os.path.join(W4DIR, 'untwisted-exact-%d%d%d.txt' % tuple(lam))
    out = dict(sectors={}, order=[])
    cur = None
    for line in open(fn, encoding='utf8'):
        line = line.rstrip('\n')
        if not line or line.startswith('#'):
            continue
        f = line.split('\t')
        tag = f[0]
        if tag in ('REP', 'DIM', 'GTBASIS', 'WEIGHTS1', 'GMET1', 'THETA'):
            out[tag.lower()] = _parse(f[1])
        elif tag == 'XVEC':
            out['xvec'] = _parse(f[1])
        elif tag == 'SECTOR':
            mu = tuple(_parse(f[1]))
            cur = out['sectors'].setdefault(mu, dict(mult=int(f[2]), wdim=int(f[3]), hwbasis=[], branches=[]))
            out['order'].append(mu)
        elif tag == 'HWBASIS':
            cur['hwbasis'].append({int(i): v for i, v in _parse(f[3])})
        elif tag == 'GRAM':
            cur['gram'] = _parse(f[2])
        elif tag == 'LAURMAT':
            cur['laurmat'] = [_parse(x) for x in f[2:]]
        elif tag == 'BRANCH':
            mu = tuple(_parse(f[1]))
            cp = int(f[2]); D0 = int(f[3])
            vec = _parse(f[4]); g2 = _parse(f[5]); cg = _parse(f[6])
            eta, alpha, beta = _parse(f[7]), _parse(f[8]), _parse(f[9])
            out['sectors'][mu]['branches'].append(dict(copy=cp, D0=D0, vec=vec, gnorm2=g2, cg=cg, eta=eta, alpha=alpha, beta=beta))
    return out


def pairval(p, D0):
    return mp.mpf(p[0].numerator if isinstance(p[0], Fr) else p[0]) / (p[0].denominator if isinstance(p[0], Fr) else 1) + \
        (mp.mpf(p[1].numerator if isinstance(p[1], Fr) else p[1]) / (p[1].denominator if isinstance(p[1], Fr) else 1)) * mp.sqrt(D0)


def w4_to_mine(lam, w4, ch):
    """map W4 one-site GT index (0-based) -> gl3_reps index"""
    pats = {tuple(map(tuple, p)): i for i, p in enumerate(ch.rep.patterns)}
    return [pats[tuple(tuple(r) for r in p)] for p in w4['gtbasis']]


# ================================================================================================ B2: regime points
from gaudin_chain import Chain, Operators
DPS = 200
THETA = (Fr(2, 7), Fr(-3, 5))
PHI = (Fr(1, 3), Fr(-2, 7), Fr(5, 11))
PHI2 = (Fr(-2, 5), Fr(3, 8), Fr(7, 11))
EPS0 = Fr(1, 10 ** 6)
GRID1 = [Fr(3, 4), Fr(-5, 3), Fr(1, 2), Fr(-7, 3)]
GRID2 = [Fr(-12, 5), Fr(7, 2), Fr(1, 8), Fr(5, 4)]
_CHS = {}


def chain_of(lam):
    lam = tuple(lam)
    if lam not in _CHS:
        _CHS[lam] = Chain(lam=lam, theta=THETA)
    return _CHS[lam]


def eps_k(k, eps0=EPS0):
    return eps0 / 2 ** k


def raising_res(ch, vec):
    nv = mp.sqrt(mp.fsum(abs(x) ** 2 for x in vec))
    tot = mp.mpf(0)
    for (i, j) in ((1, 2), (2, 3)):
        E = ch.Es[(1, i, j)] + ch.Es[(2, i, j)]
        w = E.apply(vec)
        tot += mp.fsum(abs(x) ** 2 for x in w)
    return mp.sqrt(tot) / nv


def eig_op(op, vec, piv):
    w = op.apply(vec)
    lam = w[piv] / vec[piv]
    return lam, max(abs(w[j] - lam * vec[j]) for j in range(len(vec)))


def build_regime(lam, k, phi=PHI, tag='', dps=DPS, eps0=EPS0, force=False):
    """at eps_k, zeta = eps phi: the HW-limit branches (smallest raising residual per HW sector; copies ordered by a_1 ascending)"""
    lam = tuple(lam)
    path = os.path.join(HERE, 'nu_general_reg_%d%d%d%s_%d%s.pkl' % (lam + (tag, k, '' if eps0 == EPS0 else '_e0')))
    if os.path.exists(path) and not force:
        with open(path, 'rb') as f:
            out = pickle.load(f)
        mp.mp.dps = dps
        return out
    mp.mp.dps = dps
    w4 = load_w4(lam)
    ch = chain_of(lam)
    eps = eps_k(k, eps0)
    zeta = tuple(eps * p for p in phi)
    g = [[zeta[i] if i == j else Fr(0) for j in range(3)] for i in range(3)]
    ops = Operators(ch, g)
    C = ops.fit()
    assert C['res2'] == 0 and C['res3'] == 0
    u0, u1 = Fr(3, 7), Fr(-2, 9)
    A2 = ops.eval_S(C, 2, u0); A3 = ops.eval_S(C, 3, u1)
    pi3 = mp.mpf(1) / 3 * mp.pi
    brs = []; gaps = {}
    for mu in w4['order']:
        mult = w4['sectors'][mu]['mult']
        idx = ch.sectors[mu]
        M = A2.to_mp(idx) + pi3 * A3.to_mp(idx)
        if len(idx) == 1:
            E = [M[0, 0]]; ER = mp.matrix([[1]])
        else:
            E, ER = mp.eig(M)
        cand = []
        for kk in range(len(idx)):
            v = [mp.mpf(0)] * ch.dim
            for a, j in enumerate(idx):
                v[j] = ER[a, kk]
            m_ = max(range(ch.dim), key=lambda j: abs(v[j]))
            v = [x / v[m_] for x in v]
            # real: the problem is real
            im = max(abs(mp.im(x)) for x in v)
            assert im < mp.mpf(10) ** (-dps // 2), im
            v = [mp.re(x) for x in v]
            cand.append(dict(ev=mp.re(E[kk]), vec=v, piv=m_, rr=raising_res(ch, v)))
        cand.sort(key=lambda c: c['rr'])
        gaps[mu] = [c['rr'] for c in cand[:mult + 1]]
        sel = cand[:mult]
        for c in sel:
            lam2 = [eig_op(op, c['vec'], c['piv']) for op in C[2]]
            lam3 = [eig_op(op, c['vec'], c['piv']) for op in C[3]]
            c['lam2'] = [x[0] for x in lam2]; c['lam3'] = [x[0] for x in lam3]
            c['eres'] = max(x[1] for x in lam2 + lam3)
            c['sector'] = mu
        sel.sort(key=lambda c: c['lam3'][2])          # copy 1 = smaller a_1
        for ci, c in enumerate(sel):
            c['copy'] = ci + 1
            c['idxs'] = idx
            c['vec'] = {j: c['vec'][j] for j in idx}
            brs.append(c)
    out = dict(lam=lam, k=k, eps=eps, zeta=zeta, phi=phi, dps=dps, branches=brs, gaps=gaps)
    with open(path + '.tmp', 'wb') as f:
        pickle.dump(out, f)
    os.replace(path + '.tmp', path)
    return out


# ================================================================================================ B2: analysis
import general_kernel as gk
from gaudin_chain import exp_nil

NPTS = 12
_SF = {}


class SeriesFitter:
    """f(eps_k) = sum_{n<npts} c_n eps_k^n through the first npts points, eps_k = eps0 2^-k (x = 2^-k scaling for conditioning); fit(vals) -> true-eps coefficients"""
    def __init__(self, npts=NPTS, eps0=EPS0):
        self.npts = npts
        e0 = mp.mpf(eps0.numerator) / eps0.denominator
        A = mp.matrix(npts, npts)
        for i in range(npts):
            x = mp.mpf(2) ** (-i)
            for j in range(npts):
                A[i, j] = x ** j
        self.Ainv = A ** -1
        self.scale = [e0 ** j for j in range(npts)]

    def fit(self, vals):
        npts = self.npts
        out = []
        for j in range(npts):
            s = mp.fsum(self.Ainv[j, i] * vals[i] for i in range(npts))
            out.append(s / self.scale[j])
        return out


def fitter(npts=NPTS):
    key = (npts, mp.mp.dps)
    if key not in _SF:
        _SF[key] = SeriesFitter(npts)
    return _SF[key]


def gram_diag(ch):
    d = [None] * ch.d1
    d[ch.hw_index] = Fr(1)
    changed = True
    while changed:
        changed = False
        for (i, j) in ((1, 2), (2, 3)):
            Eu = ch.E1[(i, j)]; El = ch.E1[(j, i)]
            for a, row in Eu.r.items():
                for b, v in row.items():
                    lb = El.r.get(b, {}).get(a)
                    if d[b] is not None and d[a] is None and v != 0:
                        d[a] = lb * d[b] / v; changed = True
                    elif d[a] is not None and d[b] is None and lb is not None and lb != 0:
                        d[b] = v * d[a] / lb; changed = True
    return d


class Rep:
    """per-rep data: chain, W4 data, index maps, gram metric, x(phi)"""
    def __init__(self, lam):
        self.lam = tuple(lam)
        self.ch = chain_of(lam)
        self.d = self.ch.d1
        self.w4 = load_w4(lam)
        self.w2m = w4_to_mine(lam, self.w4, self.ch)          # W4 one-site index (0-based) -> mine
        self.m2w = {m: w for w, m in enumerate(self.w2m)}
        self.gd = gram_diag(self.ch)
        self.gdm = [mp.mpf(g.numerator) / g.denominator for g in self.gd]
        self.one = gk.one(lam)
        self.hw = self.ch.hw_index
        # my two-site basis indices ordered by the W4 two-site index
        self.w4order = sorted(range(self.d * self.d), key=lambda ab: self.m2w[ab // self.d] * self.d + self.m2w[ab % self.d])
        self.prs = pairs_of(lam)

    def gnorm(self, v):
        d = self.d
        n2 = mp.fsum(self.gdm[a // d] * self.gdm[a % d] * v[a] ** 2 for a in range(d * d) if v[a] != 0)
        return [x / mp.sqrt(n2) for x in v]

    def full(self, br):
        v = [mp.mpf(0)] * (self.d * self.d)
        for j, x in br['vec'].items():
            v[j] = x
        return v

    def first_nonzero(self, v, tol):
        for ab in self.w4order:
            if abs(v[ab]) > tol:
                return ab
        return None

    def xcov2(self, phi):
        x = gk.x_cov(self.one, phi)
        xm = [_mp(q) for q in x]
        return [a * b for a in xm for b in xm]

    def phi_u(self, phi, w1, w2):
        """one-site covector e1^T exp(a E12) exp(g E13) exp(c E23) (exact Fractions) with the closed-form Gauss parameters of the plan, as mp list"""
        p1, p2, p3 = phi
        w1, w2 = Fr(w1), Fr(w2)
        a = (p1 - p2) * (w1 - w2) / ((p2 - w1) * (p2 - w2))
        c = (p1 - p3) * (p2 - p3) * (w1 - w2) / ((p1 - p2) * (p3 - w1) * (p3 - w2))
        g = (p3 - p1) * (p3 - p2) * (w1 - w2) / ((p3 - w1) * (p2 - w2) * (p3 - w2))
        e = [Fr(0)] * self.d
        e[self.hw] = Fr(1)
        v = exp_nil(self.ch.E1[(1, 2)], a, 10).rapply(e)
        v = exp_nil(self.ch.E1[(1, 3)], g, 10).rapply(v)
        v = exp_nil(self.ch.E1[(2, 3)], c, 10).rapply(v)
        return [_mp(x) for x in v]

    def chat(self, phi, w1, w2):
        l1, l2 = self.lam[0], self.lam[1]
        p1, p2, p3 = phi
        f = lambda w: (p2 - Fr(w)) ** (l1 - l2) * (p3 - Fr(w)) ** l1
        return _mp(f(w1) * f(w2) * (p1 - p2) ** (2 * l2))


_REPS = {}
def rep_of(lam):
    lam = tuple(lam)
    if lam not in _REPS:
        _REPS[lam] = Rep(lam)
    return _REPS[lam]


def load_pts(lam, tag='', phi=PHI, npts=NPTS):
    return [build_regime(lam, k, phi, tag) for k in range(npts)]


def branch_psi(R, pts, bi):
    vs = [R.gnorm(R.full(p['branches'][bi])) for p in pts]
    for k in range(len(vs) - 2, -1, -1):             # align each with the next (smaller eps) one
        if mp.fsum(a * b for a, b in zip(vs[k], vs[k + 1])) < 0:
            vs[k] = [-x for x in vs[k]]
    return vs


def comp_series(vs, npts=NPTS):
    F = fitter(npts)
    n = len(vs[0])
    z = [mp.mpf(0)] * npts
    cols = [F.fit([v[j] for v in vs[:npts]]) if any(v[j] != 0 for v in vs[:npts]) else z for j in range(n)]
    return [[cols[j][m] for j in range(n)] for m in range(npts)]


def nu_series(lam, pts, bi, npts=NPTS):
    """per site s, per pair i: eps-series coefficients of nu_i(theta_s), from the Laurent data of branch bi at the npts regime points"""
    F = fitter(npts)
    vals = []
    th = THETA
    for k, p in enumerate(pts[:npts]):
        b = p['branches'][bi]
        L = [laurent_of(b['lam2'], b['lam3'], s) for s in (0, 1)]
        vals.append([nu_numeric(lam, L[s], L[1 - s], p['zeta'], th[s], th[1 - s]) for s in (0, 1)])
    npr = len(vals[0][0])
    return [[F.fit([vals[k][s][i] for k in range(npts)]) for i in range(npr)] for s in (0, 1)]


def laurent_series(pts, bi, npts=NPTS):
    """eps-series of (H_s, a_s, b_s), s = 0, 1"""
    F = fitter(npts)
    out = []
    for s in (0, 1):
        out.append([F.fit([laurent_of(p['branches'][bi]['lam2'], p['branches'][bi]['lam3'], s)[q] for p in pts[:npts]]) for q in range(3)])
    return out


_GT = {}
def G_terms(lam):
    """{(i,j): [(alpha, beta, d, coef)]}: G_ij(eps w^1, eps w^1'; t) = sum coef eps^d w^1^alpha w^1'^beta, G = (w w')^{2l1-l2} Theta, t = 1/theta12"""
    lam = tuple(lam)
    if lam in _GT:
        return _GT[lam]
    T, rep = gk.load_rep(lam)
    l1, l2 = lam[0], lam[1]
    t = 1 / (THETA[0] - THETA[1])
    out = {}
    for (i, j), ent in T.items():
        terms = {}
        for (n, p, q), c in ent.items():
            al, be = 2 * l1 - l2 - p, 2 * l1 - l2 - q
            assert al >= 0 and be >= 0, (lam, i, j, n, p, q)
            d = 4 * l1 - 2 * l2 - p - q
            key = (al, be, d)
            terms[key] = terms.get(key, 0) + Fr(c) * t ** n
        out[(i, j)] = [(al, be, d, c) for (al, be, d), c in terms.items() if c != 0]
    _GT[lam] = out
    return out


def leading_poly(lam, nus, tmax):
    """R[t] = {(a,b): coeff}, t = d + n + m <= tmax: sum_ij G_ij^[d] nu_i^(n)(theta1) nu_j^(m)(theta2);
    contrib[t] = {(i,n,j,m): sum |term|}"""
    G = G_terms(lam)
    R = {}; contrib = {}
    for (i, j), terms in G.items():
        for (al, be, d, c0) in terms:
            c = _mp(c0)
            for n in range(0, tmax - d + 1):
                x1 = nus[0][i][n]
                for m in range(0, tmax - d - n + 1):
                    v = c * x1 * nus[1][j][m]
                    t = d + n + m
                    Rt = R.setdefault(t, {})
                    Rt[(al, be)] = Rt.get((al, be), mp.mpf(0)) + v
                    ct = contrib.setdefault(t, {})
                    ct[(i, n, j, m)] = ct.get((i, n, j, m), mp.mpf(0)) + abs(v)
    return R, contrib


def pe(P, w1, w2):
    w1, w2 = _mp(w1), _mp(w2)
    return mp.fsum(c * w1 ** a * w2 ** b for (a, b), c in P.items())


def lead_exp(c, tol=mp.mpf(10) ** -25):
    for n, x in enumerate(c):
        if abs(x) > tol:
            return n
    return None


# ---------------------------------------------------------------------------------------------- the analysis of one rep / phi
GRID8_1 = [Fr(3, 4), Fr(-5, 3), Fr(1, 2), Fr(-7, 3), Fr(9, 7), Fr(-4, 9), Fr(11, 5), Fr(-13, 6)]
GRID8_2 = [Fr(-12, 5), Fr(7, 2), Fr(1, 8), Fr(5, 4), Fr(-8, 3), Fr(13, 9), Fr(-3, 7), Fr(17, 6)]
_HWK = {}


def hw_kets(R, mu):
    """W4's exact HW kets of sector mu mapped to my two-site basis (list of lists of Fractions), checked exactly annihilated by E12, E23 (total)"""
    key = (R.lam, mu)
    if key in _HWK:
        return _HWK[key]
    d = R.d
    out = []
    for hb in R.w4['sectors'][mu]['hwbasis']:
        v = [Fr(0)] * (d * d)
        for idx, val in hb.items():
            i, j = divmod(idx - 1, d)
            v[R.w2m[i] * d + R.w2m[j]] = Fr(val)
        for (a, b) in ((1, 2), (2, 3)):
            E = R.ch.Es[(1, a, b)] + R.ch.Es[(2, a, b)]
            assert all(x == 0 for x in E.apply(v)), 'W4 HW ket not annihilated'
        out.append(v)
    _HWK[key] = out
    return out


def w4_cg(R, mu, copy):
    """exact normalised CG slice C_nu = <e_nu (x) e_1|Psi0> in MY one-site index (list of mp, length d) for W4 branch (mu, copy), with W4's sign"""
    b = R.w4['sectors'][mu]['branches'][copy - 1]
    D0 = b['D0']
    sq = mp.sqrt(D0) if D0 else mp.mpf(0)
    def pv(p):
        return _mp(Fr(p[0])) + sq * _mp(Fr(p[1]))
    g2 = pv(b['gnorm2'])
    out = [mp.mpf(0)] * R.d
    for (nu, ca, cb) in b['cg']:
        out[R.w2m[nu - 1]] = (_mp(Fr(ca)) + sq * _mp(Fr(cb))) / mp.sqrt(g2)
    return out


def w4_full(R, mu, copy):
    b = R.w4['sectors'][mu]['branches'][copy - 1]
    D0 = b['D0']
    sq = mp.sqrt(D0) if D0 else mp.mpf(0)
    d = R.d
    v = [mp.mpf(0)] * (d * d)
    for (idx, va, vb) in b['vec']:
        i, j = divmod(idx - 1, d)
        v[R.w2m[i] * d + R.w2m[j]] = _mp(Fr(va)) + sq * _mp(Fr(vb))
    return R.gnorm(v)


def analyse(lam, tag='', phi=PHI, npts=NPTS, bis=None, verbose=True, nmax_x=10, grid=None):
    lam = tuple(lam)
    R = rep_of(lam)
    mp.mp.dps = DPS
    pts = load_pts(lam, tag, phi, npts)
    l1, l2 = lam[0], lam[1]
    d = R.d; hw = R.hw
    w4 = R.w4
    blist = [(mu, cp) for mu in w4['order'] for cp in range(1, w4['sectors'][mu]['mult'] + 1)]
    xx = R.xcov2(phi)
    grid = grid or [(a, b) for a in GRID8_1 for b in GRID8_2]
    grid16 = [(a, b) for a in GRID1 for b in GRID2]
    cov_cache = {}
    # eps-freeness of the vacua (exact)
    e0 = eps_k(0)
    import twovac_fit as T
    c1 = T.cov_inf(R.ch, tuple(e0 * p for p in phi), e0 * GRID1[0], e0 * GRID2[0])
    c2 = T.cov_inf(R.ch, phi, GRID1[0], GRID2[0])
    epsfree = (c1 == c2)
    out = dict(lam=lam, tag=tag, phi=phi, epsfree=epsfree, br=[])
    for bi, (mu, cp) in enumerate(blist):
        if bis is not None and bi not in bis:
            continue
        t0 = time.time()
        vs = branch_psi(R, pts, bi)
        PS = comp_series(vs, npts)
        ab = R.first_nonzero(PS[0], mp.mpf(10) ** -25)
        if PS[0][ab] < 0:
            vs = [[-x for x in v] for v in vs]
            PS = [[-x for x in v] for v in PS]
        psi0 = PS[0]
        xn = [mp.fsum(a * b for a, b in zip(xx, PS[n]) if a != 0) for n in range(nmax_x)]
        kB = lead_exp(xn)
        XB = xn[kB]
        nus = nu_series(lam, pts, bi, npts)
        nexp = [[lead_exp(nus[s][i]) for i in range(len(nus[s]))] for s in (0, 1)]
        NBp = 4 * l1 - 2 * l2 - kB
        tmax = min(NBp + 1, npts - 1)
        Rr, contrib = leading_poly(lam, nus, tmax)
        order = [(t, max(abs(v) for v in Rr[t].values()), max(contrib[t].values())) for t in sorted(Rr)]
        NB = next((t for t, v, c in order if v > mp.mpf(10) ** -25), None)
        LB = Rr[NB]
        low = max([v / c for t, v, c in order if t < NB and c > mp.mpf(10) ** -30] + [mp.mpf(0)])
        low_abs = max([v for t, v, c in order if t < NB] + [mp.mpf(0)])
        # the identity at the 64 points
        rows = []
        for (w1, w2) in grid:
            Lp = pe(LB, w1, w2)
            cph = R.chat(phi, w1, w2)
            phu = R.phi_u(phi, w1, w2)
            gf = mp.fsum(phu[nu] * psi0[nu * d + hw] for nu in range(d))
            rows.append((w1, w2, Lp, cph * gf / XB))
        num = mp.fsum(a[2] * a[3] for a in rows); den = mp.fsum(a[3] * a[3] for a in rows)
        Cc = num / den
        resid = mp.sqrt(mp.fsum((a[2] - Cc * a[3]) ** 2 for a in rows)) / mp.sqrt(mp.fsum(a[2] ** 2 for a in rows))
        rat = [a[2] / a[3] for a in rows if abs(a[3]) > mp.mpf(10) ** -30]
        spread = max(abs(r / Cc - 1) for r in rat) if rat else None
        # dictionary: solve L_B/chat in span over HW kets
        H = hw_kets(R, mu)
        Hm = [[_mp(x) for x in h] for h in H]
        A = mp.matrix(len(grid), len(H))
        rhs = mp.matrix(len(grid), 1)
        for p_, (w1, w2) in enumerate(grid):
            cph = R.chat(phi, w1, w2)
            phu = R.phi_u(phi, w1, w2)
            rhs[p_] = pe(LB, w1, w2)
            for i, h in enumerate(Hm):
                A[p_, i] = cph * mp.fsum(phu[nu] * h[nu * d + hw] for nu in range(d))
        beta = mp.lu_solve(A.T * A, A.T * rhs) if len(H) > 1 else mp.matrix([(A.T * rhs)[0] / (A.T * A)[0]])
        dres = mp.norm(A * beta - rhs) / mp.norm(rhs)
        vq = [mp.fsum(beta[i] * Hm[i][j] for i in range(len(H))) for j in range(d * d)]
        vq = R.gnorm(vq)
        ab2 = R.first_nonzero(vq, mp.mpf(10) ** -30)
        if vq[ab2] < 0:
            vq = [-x for x in vq]
        slq = [vq[nu * d + hw] for nu in range(d)]
        ex = w4_cg(R, mu, cp)
        dev_cg = max(abs(a - b) for a, b in zip(slq, ex))
        dev_psi = max(abs(a - b) for a, b in zip(vq, psi0))
        exf = w4_full(R, mu, cp)
        dev_w4 = max(abs(a - b) for a, b in zip(psi0, exf))
        dev_w4q = max(abs(a - b) for a, b in zip(vq, exf))
        # Phi limit: <Phi_inf|Psi_B>(eps->0) vs <phi_u (x) e1|Psi0> at the 16 grid points
        phirat = []; phiexp = []
        for (w1, w2) in grid16:
            if (w1, w2) not in cov_cache:
                cov_cache[(w1, w2)] = [_mp(x) for x in T.cov_inf(R.ch, phi, w1, w2)]
            cov = cov_cache[(w1, w2)]
            s_ = fitter(npts).fit([mp.fsum(c * v for c, v in zip(cov, vk) if c != 0) for vk in vs[:npts]])
            phu = R.phi_u(phi, w1, w2)
            gf = mp.fsum(phu[nu] * psi0[nu * d + hw] for nu in range(d))
            phiexp.append(lead_exp(s_))
            phirat.append(abs(s_[0] / gf - 1))
        # depth: which (i, n) enter the leading polynomial
        cN = contrib[NB]
        thr = max(cN.values()) * mp.mpf(10) ** -25
        depth = [{}, {}]
        for (i, n, j, m), v in cN.items():
            if v > thr:
                depth[0].setdefault(i, set()).add(n); depth[1].setdefault(j, set()).add(m)
        rec = dict(bi=bi, mu=mu, copy=cp, kB=kB, XB=XB, xn=[abs(x) for x in xn[:8]], NB=NB, NBpred=NBp, order=order, low=low, low_abs=low_abs, C=Cc, resid=resid,
                   spread=spread, dres=dres, dev_cg=dev_cg, dev_psi=dev_psi, dev_w4=dev_w4, dev_w4q=dev_w4q, slq=slq, slex=ex, nexp=nexp, phirat=max(phirat),
                   phiexp=sorted(set(phiexp)), depth=[{i: sorted(v) for i, v in dd.items()} for dd in depth], LB=LB, psi0=psi0, time=time.time() - t0)
        out['br'].append(rec)
        if verbose:
            print('%s#%d: k_B=%d N_B=%d (pred %d)  C=%s  resid=%s spread=%s | lower orders max %s | dict resid %s, |slice-W4|=%s |v_Q-Psi0|=%s |Psi0-W4|=%s | Phi lead dev %s exps %s | %.0fs' % (
                mu, cp, kB, NB, NBp, mp.nstr(Cc, 12), mp.nstr(resid, 3), mp.nstr(spread, 3) if spread is not None else None, mp.nstr(low_abs, 3), mp.nstr(dres, 3),
                mp.nstr(dev_cg, 3), mp.nstr(dev_psi, 3), mp.nstr(dev_w4, 3), mp.nstr(max(phirat), 3), sorted(set(phiexp)), rec['time']), flush=True)
    return out


def run_all(lam, npts=24):
    lam = tuple(lam)
    res = {}
    for tag, phi in (('', PHI), ('2', PHI2)):
        path = os.path.join(HERE, 'nu_general_res_%d%d%d%s.pkl' % (lam + (tag,)))
        print('######## %s phi=%s npts=%d' % (lam, [str(x) for x in phi], npts), flush=True)
        o = analyse(lam, tag, phi, npts=npts)
        with open(path, 'wb') as f:
            pickle.dump(o, f)
        res[tag] = o
    return res


if __name__ == '__main__' and 'RUN' in sys.argv:
    run_all(tuple(int(c) for c in sys.argv[sys.argv.index('RUN') + 1]), int(sys.argv[sys.argv.index('RUN') + 2]) if len(sys.argv) > sys.argv.index('RUN') + 2 else 24)
