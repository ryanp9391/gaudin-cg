"""lattice_recursion_gl4.py -- G8a: block-Jacobi structure of the gl(4) triple (V3, V4a, V4b) with respect to lattice gradings, and the forward recursion that builds a slice
from its eigenvalue triple (2026-10-24).

Setting (see RESULTS-G8a-lattice-recursion.md and SPEC-G8a): the weight space W = V_lam[mu - lam'] of a sector is a region of the GT lattice with free coordinates
x = (m12, m13, m23).  A grading g = a m12 + b m13 + c m23 cuts it into levels L_d = {g = d}.  If every operator of the triple is block-banded with respect to the levels and the stacked
inter-level blocks F_d = [M_k[L_d, L_{d+B}]]_k (3 n_d x n_{d+B}) have full column rank (greedy for band B > 1), the eigenvector equations (V_k - v_k) psi = 0 determine psi level by level
from the bottom level L_0:   F_d psi_{d+1} = [ v_k psi_d - sum_{c known} M_k[L_d, c] psi_c ]_k .
The rows of the equations that are not used (3 n_d - n_{d+1} per level, and all 3 n_top at the top) are the CHARACTERISTIC EQUATIONS C_j(v) (psi_0) = 0.

Public functions: sector, analyze, build_plan, slice_from_triple, sym_recursion.  CLI: python -X utf8 lattice_recursion_gl4.py gradings|recursion|charpoly|n0|separability|all|run <stages>
"""
import os, sys, time, itertools, math, random, subprocess, json
from fractions import Fraction as Fr
from collections import defaultdict, Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import mpmath as mp
import glN_reps as gr
import slice_engine_gl4 as se
import slice_ops_gl4 as so
import structure_gl4 as sg

TRIPLE = so.TRIPLE
HERE = os.path.dirname(os.path.abspath(__file__))
TABLES_PATH = os.path.join(HERE, 'lattice_recursion_gl4_tables.txt')
OUTPUT_PATH = os.path.join(HERE, 'lattice_recursion_gl4_output.txt')

P_SMALL = [((2, 1, 0, 0), (2, 1, 0, 0)), ((2, 1, 1, 0), (2, 1, 1, 0)), ((3, 2, 1, 0), (3, 2, 1, 0)), ((4, 2, 0, 0), (4, 2, 0, 0)), ((4, 2, 1, 0), (4, 2, 1, 0)),
           ((3, 2, 1, 0), (4, 2, 1, 0)), ((4, 2, 1, 0), (3, 2, 1, 0)), ((4, 2, 2, 0), (4, 2, 2, 0)), ((4, 3, 1, 0), (4, 3, 1, 0))]
P_BIG = [((6, 4, 2, 0), (6, 4, 2, 0)), ((5, 3, 1, 0), (5, 3, 1, 0)), ((5, 3, 2, 1), (5, 3, 2, 1))]

_OUT = []
_TAB = []


def out(s=''):
    print(s, flush=True)
    _OUT.append(s)


def tab(s=''):
    _TAB.append(s)


def flush_tables(mode='a', title=None):
    with open(TABLES_PATH, mode, encoding='utf-8') as f:
        if title:
            f.write('\n' + '=' * 150 + '\n' + title + '\n' + '=' * 150 + '\n')
        f.write('\n'.join(_TAB) + '\n')
    _TAB.clear()


def hdr(s):
    out('')
    out('=' * 130)
    out(s)
    out('=' * 130)


def lam_str(l):
    return '[' + ','.join(str(x) for x in l) + ']'


def prod_str(lam, lamp):
    return lam_str(lam) + ('^2' if lam == lamp else 'x' + lam_str(lamp))


def mu_str(mu):
    return '(' + ','.join(str(x) for x in mu) + ')'


# ======================================================================================== 0. sector data
_OPS = {}
_SEC = {}


def reduced_ops(lamp):
    lamp = tuple(lamp)
    if lamp not in _OPS:
        _OPS[lamp] = {n: so.reduced_op(n, lamp) for n in TRIPLE}
    return _OPS[lamp]


class Space:
    """a coordinate subspace of a weight space: positions pos (indices into the sector's W), exact matrices M[k][r][c] restricted, lattice coordinates xs"""

    def __init__(self, sec, name, pos):
        self.sec, self.name, self.pos = sec, name, list(pos)
        self.n = len(self.pos)
        self.xs = [sec.xs[p] for p in self.pos]
        self.M = [[[sec.M[k][r][c] for c in self.pos] for r in self.pos] for k in range(3)]
        self.ents = sorted({(r, c) for k in range(3) for r in range(self.n) for c in range(self.n) if self.M[k][r][c] != 0})

    def grade(self, g):
        a, b, c = g
        return [a * x[0] + b * x[1] + c * x[2] for x in self.xs]


class Sector:
    def __init__(self, lam, lamp, mu):
        self.lam, self.lamp, self.mu = tuple(lam), tuple(lamp), tuple(mu)
        self.R = gr.GLNRep.of(self.lam)
        self.nu = tuple(mu[i] - lamp[i] for i in range(4))
        W = self.R.wsp.get(self.nu)
        self.ok = W is not None
        if not self.ok:
            return
        self.W = np.asarray(W)
        assert all(self.W[i] < self.W[i + 1] for i in range(len(self.W) - 1))
        self.w = len(self.W)
        mats = so.weight_space_mats(self.R, reduced_ops(lamp), self.W)
        self.M = [mats[n] for n in TRIPLE]
        self.pats = [self.R.pats[int(c)] for c in self.W]
        self.xs = [sg.coords(p) for p in self.pats]
        self.k1 = lamp[0] - lamp[1] + 1
        self.cut = mu[0] - lamp[1]                                   # E12^k1 kills e_x iff m12 <= mu1 - lam'_2
        self.in1 = [i for i, x in enumerate(self.xs) if x[0] <= self.cut]
        self.active = frozenset(self._active())
        self._m = None
        self._spaces = {}
        self._prv = None

    def _active(self):
        out_ = set()
        for i in range(1, 4):
            k = self.lamp[i - 1] - self.lamp[i] + 1
            tw = list(self.nu)
            tw[i - 1] += k
            tw[i] -= k
            if tuple(tw) in self.R.wsp:
                out_.add(i)
        return out_

    @property
    def prv(self):
        if self._prv is None:
            self._prv = se.prv_space(self.lam, self.lamp, self.mu)
        return self._prv

    @property
    def m(self):
        return self.prv.m

    def space(self, name):
        if name not in self._spaces:
            self._spaces[name] = Space(self, name, range(self.w) if name == 'W' else self.in1)
        return self._spaces[name]

    def w1_invariant(self):
        s1 = set(self.in1)
        return all(self.M[k][r][c] == 0 for k in range(3) for c in s1 for r in range(self.w) if r not in s1)

    def cut_is_kernel(self):
        """True iff {c : E12^k1 e_c = 0} equals the coordinate set m12 <= mu1 - lam'_2"""
        imgs = se.power_images(self.R, (1, 2), self.k1, self.W)
        killed = {i for i, v in enumerate(imgs) if not v}
        return killed == set(self.in1)


def sector(lam, lamp, mu):
    key = (tuple(lam), tuple(lamp), tuple(mu))
    if key not in _SEC:
        _SEC[key] = Sector(*key)
    return _SEC[key]


def product_sectors(lam, lamp, min_w=2):
    """all Sector objects (dominant mu, non-empty weight space of dimension >= min_w) of the product, descending in mu"""
    res = []
    for mu, mult in sg.sector_list(lam, lamp, 'all'):
        s = sector(lam, lamp, mu)
        if s.ok and s.w >= min_w:
            res.append(s)
    return res


# ======================================================================================== 1. exact rank, gradings, forward determinacy
def frank(rows):
    """exact rank over Q of a list of rows of Fractions"""
    M = [list(r) for r in rows if any(x != 0 for x in r)]
    if not M:
        return 0
    nc = len(M[0])
    rk = 0
    for c in range(nc):
        p = None
        for i in range(rk, len(M)):
            if M[i][c] != 0:
                p = i
                break
        if p is None:
            continue
        M[rk], M[p] = M[p], M[rk]
        pr = M[rk]
        inv = 1 / pr[c]
        for i in range(rk + 1, len(M)):
            f = M[i][c]
            if f != 0:
                f = f * inv
                row = M[i]
                M[i] = [x - f * y for x, y in zip(row, pr)]
        rk += 1
        if rk == len(M):
            break
    return rk


GRADINGS = []                       # 124 directed gradings (a, b, c), g and -g both
for _t in itertools.product(range(-2, 3), repeat=3):
    if _t != (0, 0, 0):
        GRADINGS.append(_t)


def pair_rep(g):
    """canonical representative of the pair {g, -g}"""
    return max(g, tuple(-x for x in g))


def levels_of(sp, g):
    gv = sp.grade(g)
    lv = sorted(set(gv))
    lid = {d: i for i, d in enumerate(lv)}
    li = [lid[x] for x in gv]
    groups = [[] for _ in lv]
    for i, l in enumerate(li):
        groups[l].append(i)
    return gv, lv, li, groups


def analyze(sp, g):
    """band (in g-values and in level steps), level sizes, forward determinacy (greedy; = 'every stacked block F_d has full column rank' for band 1).
    Returns dict(band_val, band, sizes, n0, ntop, ok, fail, deficit)"""
    gv, lv, li, groups = levels_of(sp, g)
    D = len(lv)
    sizes = [len(x) for x in groups]
    if sp.ents:
        band_val = max(abs(gv[r] - gv[c]) for r, c in sp.ents)
        B = max(abs(li[r] - li[c]) for r, c in sp.ents)
    else:
        band_val = B = 0
    det = set(groups[0])
    ok, fail, deficit = True, None, 0
    for d in range(D):
        U = [c for l in range(d + 1, min(d + B, D - 1) + 1) for c in groups[l] if c not in det]
        if U:
            rows = [[sp.M[k][r][c] for c in U] for k in range(3) for r in groups[d]]
            rk = frank(rows)
            if rk < len(U):
                ok, fail, deficit = False, d, len(U) - rk
                break
            det.update(U)
    return dict(band_val=band_val, band=B, sizes=sizes, n0=sizes[0], ntop=sizes[-1], ok=ok, fail=fail, deficit=deficit, D=D)


# ======================================================================================== 2. the recursion plan and its evaluation
class Plan:
    """exact data of the forward recursion on a Space for a grading g: init = positions of the bottom level, steps = list of dict(rows, U, known, coef[k][ri][ci], Lm, T)"""
    pass


def build_plan(sp, g):
    from slice_engine_gl4 import exact_left_inverse
    gv, lv, li, groups = levels_of(sp, g)
    D = len(lv)
    B = max([abs(li[r] - li[c]) for r, c in sp.ents] + [0])
    pl = Plan()
    pl.sp, pl.g, pl.init = sp, g, list(groups[0])
    pl.n0 = len(pl.init)
    pl.steps = []
    pl.ok = True
    pl.fail = None
    det = list(groups[0])
    detset = set(det)
    for d in range(D):
        rows = groups[d]
        U = [c for l in range(d + 1, min(d + B, D - 1) + 1) for c in groups[l] if c not in detset]
        known = sorted(det)
        coef = [[[sp.M[k][r][c] for c in known] for r in rows] for k in range(3)]
        st = dict(rows=rows, U=U, known=known, coef=coef)
        nr = 3 * len(rows)
        if U:
            A = [[sp.M[k][r][c] for c in U] for k in range(3) for r in rows]
            try:
                Lm, T = exact_left_inverse(A, len(U))
            except ValueError:
                pl.ok, pl.fail = False, d
                return pl
            st['Lm'], st['T'] = Lm, T
        else:
            st['Lm'] = []
            st['T'] = [[Fr(int(i == j)) for j in range(nr)] for i in range(nr)]
        pl.steps.append(st)
        det += U
        detset.update(U)
    pl.nunused = sum(len(st['T']) for st in pl.steps)
    pl._conv = {}
    return pl


def _prep(pl, key, conv):
    if key not in pl._conv:
        steps = []
        for st in pl.steps:
            steps.append(dict(rows=st['rows'], U=st['U'], known=st['known'],
                              coef=[[[conv(x) for x in row] for row in M] for M in st['coef']],
                              Lm=[[conv(x) for x in row] for row in st['Lm']], T=[[conv(x) for x in row] for row in st['T']]))
        pl._conv[key] = steps
    return pl._conv[key]


def run_plan(pl, v, psi0cols, conv, zero, key):
    """evaluate the recursion for scalars v = (v3, v4a, v4b) (numbers or polynomials) and psi0cols[i] = list (over ncol columns) of scalars at the init position i.
    Returns (psi: dict position -> list over columns, unused: list of lists [t][j])"""
    steps = _prep(pl, key, conv)
    ncol = len(psi0cols[0])
    psi = {p: list(c) for p, c in zip(pl.init, psi0cols)}
    unused = []
    for st in steps:
        rows, U, known = st['rows'], st['U'], st['known']
        rhs_all = []
        for j in range(ncol):
            rhs = []
            for k in range(3):
                Ck = st['coef'][k]
                for ri, r in enumerate(rows):
                    val = v[k] * psi[r][j]
                    row = Ck[ri]
                    for ci, c in enumerate(known):
                        x = row[ci]
                        if x != 0:
                            val = val - x * psi[c][j]
                    rhs.append(val)
            rhs_all.append(rhs)
        if U:
            sols = []
            for j in range(ncol):
                rhs = rhs_all[j]
                sols.append([_dot(Lrow, rhs, zero) for Lrow in st['Lm']])
            for ui, u in enumerate(U):
                psi[u] = [sols[j][ui] for j in range(ncol)]
        for trow in st['T']:
            unused.append([_dot(trow, rhs_all[j], zero) for j in range(ncol)])
    return psi, unused


def _dot(row, vec, zero):
    s = zero
    for x, y in zip(row, vec):
        if x != 0:
            s = s + x * y
    return s


_PLAN = {}


def get_plan(sec, space, g):
    key = (sec.lam, sec.lamp, sec.mu, space, tuple(g))
    if key not in _PLAN:
        _PLAN[key] = build_plan(sec.space(space), tuple(g))
    return _PLAN[key]


def best_gradings(sec, space, max_n0=None):
    """list of (n0, band, |g|_1, g) of the fully determined gradings of the space, sorted (smallest n0, then smallest band, then simplest)"""
    sp = sec.space(space)
    res = []
    for g in GRADINGS:
        a = analyze_cached(sec, space, g)
        if a['ok']:
            res.append((a['n0'], a['band'], sum(abs(x) for x in g), g))
    res.sort()
    if max_n0 is not None:
        res = [r for r in res if r[0] <= max_n0]
    return res


_AN = {}


def analyze_cached(sec, space, g):
    key = (sec.lam, sec.lamp, sec.mu, space, tuple(g))
    if key not in _AN:
        _AN[key] = analyze(sec.space(space), tuple(g))
    return _AN[key]


def mpq(x):
    return mp.mpf(x.numerator) / x.denominator


def slice_from_triple(lam, lamp, mu, trip, grading='best', space='W1', dps=50, psi0=None, info=None):
    """The forward lattice recursion.  Given the joint eigenvalue triple trip = (v3, v4a, v4b) of an eigenvector of the reduced operators on the coordinate subspace `space`
    ('W1' = the E12-cut subspace {m12 <= mu1 - lam'_2}, or 'W') of the weight space W = V_lam[mu - lam'], build the eigenvector level by level for the grading g = a m12 + b m13 + c m23:

        levels L_0 < L_1 < ... < L_D (distinct g-values ascending),  bottom component psi_0 = psi0 (the number 1 if the bottom level has n0 = 1 point, else a vector of length n0),
        F_d psi_{d+1} = [ v_k psi_d - sum_{c determined, c not in L_{d+1}} M_k[L_d, c] psi_c ]_{k = 3, 4a, 4b},    F_d = [M_k[L_d, L_{d+1}]]_k  (3 n_d x n_{d+1}, full column rank),

    solved with the EXACT left inverse of F_d (computed once over Q); for a band B > 1 the greedy version determines every not yet known unknown of the levels <= d + B at once.
    The 2 n + n0 - ... equations that are not needed are returned as `unused` (the characteristic equations evaluated at the triple: they vanish at a true joint eigenvalue).
    Output: dict(slice = {global GT index: mpf} in G3's normalisation (first nonzero component in ascending GT index = 1), psi = list over W positions, resid = max |unused|, grading).
    grading = 'best' (the determined grading of smallest n0 then band) or a tuple (a, b, c); raises ValueError if the grading is not fully determined or n0 > 1 and psi0 is None."""
    sec = sector(lam, lamp, mu)
    if grading == 'best':
        bg = best_gradings(sec, space)
        if not bg:
            raise ValueError('no fully determined grading')
        grading = bg[0][3]
    pl = get_plan(sec, space, grading)
    if not pl.ok:
        raise ValueError('grading %s not fully determined on %s' % (str(grading), space))
    old = mp.mp.dps
    mp.mp.dps = dps + 10
    try:
        if psi0 is None:
            if pl.n0 != 1:
                raise ValueError('n0 = %d > 1: give psi0' % pl.n0)
            psi0 = [mp.mpf(1)]
        psi, unused = run_plan(pl, [mp.mpf(x) for x in trip], [[mp.mpf(x)] for x in psi0], mpq, mp.mpf(0), ('mp', dps + 10))
        sp = sec.space(space)
        vec = [mp.mpf(0)] * sec.w
        for p, val in psi.items():
            vec[sp.pos[p]] = val[0]
        prv = sec.prv
        struct0 = [all(prv.basis[kk][a] == 0 for kk in range(prv.m)) for a in range(sec.w)]
        tolz = mp.mpf(10) ** (-(dps - 20))
        mx = max(abs(x) for x in vec)
        first = next((a for a in range(sec.w) if not struct0[a] and abs(vec[a]) > tolz * mx), None)
        resid = max([abs(u[0]) for u in unused] + [mp.mpf(0)])
        if first is None:
            return dict(slice={}, psi=vec, resid=resid, grading=grading, n0=pl.n0)
        sc = vec[first]
        nrm = [x / sc for x in vec]
        mx = max(abs(x) for x in nrm)
        nz = {int(sec.W[a]): nrm[a] for a in range(sec.w) if not struct0[a] and abs(nrm[a]) > tolz * mx}
        return dict(slice=nz, psi=vec, resid=resid, grading=grading, n0=pl.n0, first=first)
    finally:
        mp.mp.dps = old


# ======================================================================================== 3. symbolic recursion
def _ring():
    from sympy import QQ
    from sympy.polys.rings import ring
    R, v3, v4a, v4b = ring('v3,v4a,v4b', QQ)
    return R, (v3, v4a, v4b), QQ


def sym_recursion(sec, space, g, ncol=None):
    """psi = P(v3, v4a, v4b) psi0 exactly (sympy sparse polynomial ring over Q).  n0 = 1: psi0 = 1.  n0 > 1: ncol = n0 columns (psi0 = unit vectors), so P is a polynomial matrix.
    Returns dict(R, gens, psi {local position: list over columns of polynomials}, unused [[poly per column]], plan)"""
    pl = get_plan(sec, space, g)
    if not pl.ok:
        raise ValueError('not determined')
    R, gens, QQ = _ring()
    n0 = pl.n0
    psi0 = [[R.one if i == j else R.zero for j in range(n0)] for i in range(n0)]
    psi, unused = run_plan(pl, list(gens), psi0, lambda x: QQ(x.numerator, x.denominator), R.zero, 'sym')
    return dict(R=R, gens=gens, psi=psi, unused=unused, plan=pl)


def tdeg(p):
    return max((sum(m) for m in p.keys()), default=-1)


def peval(p, vals):
    s = mp.mpf(0)
    for mon, c in p.terms():
        t = mp.mpf(int(c.numerator)) / int(c.denominator)
        for e, x in zip(mon, vals):
            if e:
                t *= x ** e
        s += t
    return s


def prim_int(p):
    """make a polynomial primitive with integer coefficients and positive leading coefficient (sympy ring element -> same ring)"""
    from sympy import QQ, ZZ
    from math import gcd
    cs = [c for _, c in p.terms()]
    if not cs:
        return p
    den = 1
    for c in cs:
        den = den * int(c.denominator) // gcd(den, int(c.denominator))
    q = p * QQ(den)
    g = 0
    for _, c in q.terms():
        g = gcd(g, abs(int(c.numerator)))
    q = q * QQ(1, g)
    lm = q.LC
    if lm < 0:
        q = -q
    return q


def dgen(m):
    d = 0
    while math.comb(d + 3, 3) < m:
        d += 1
    return d


GROEBNER_CODE = r'''
import sys, json
from sympy import symbols, groebner, QQ, sympify
data = json.load(sys.stdin)
v3, v4a, v4b = symbols('v3 v4a v4b')
polys = [sympify(s) for s in data['polys']]
G = groebner(polys, v3, v4a, v4b, order='grevlex', domain=QQ)
lms = [tuple(__import__('sympy').Poly(g, v3, v4a, v4b).monoms(order='grevlex')[0]) for g in G.exprs]
if G.is_zero_dimensional:
    # count standard monomials by brute force
    def divisible(m):
        return any(all(m[i] >= l[i] for i in range(3)) for l in lms)
    seen = {(0, 0, 0)} if not divisible((0, 0, 0)) else set()
    stack = list(seen)
    while stack and len(seen) < 100000:
        m = stack.pop()
        for i in range(3):
            n = list(m); n[i] += 1; n = tuple(n)
            if n not in seen and not divisible(n):
                seen.add(n); stack.append(n)
    print(json.dumps(dict(dim=len(seen), zero_dim=True, ngens=len(lms))))
else:
    print(json.dumps(dict(dim=None, zero_dim=False, ngens=len(lms))))
'''


def groebner_dim(polys, timeout=120):
    """dimension of Q[v]/J over Q (None if not zero-dimensional or timeout) via a sympy subprocess; polys = list of sympy ring elements"""
    data = json.dumps(dict(polys=[str(p.as_expr()) for p in polys]))
    try:
        r = subprocess.run([sys.executable, '-c', GROEBNER_CODE], input=data, capture_output=True, text=True, timeout=timeout)
        if r.returncode != 0:
            return dict(dim=None, error=r.stderr[-200:])
        return json.loads(r.stdout.strip().splitlines()[-1])
    except subprocess.TimeoutExpired:
        return dict(dim=None, timeout=True)


# ======================================================================================== 4. A1: gradings
_RECS = None
STRATA = [('w2-4', 2, 4), ('w5-8', 5, 8), ('w9+', 9, 10 ** 9)]
SPACES = ('W', 'W1')


def all_recs():
    global _RECS
    if _RECS is None:
        _RECS = [(lam, lamp, s) for (lam, lamp) in P_SMALL + P_BIG for s in product_sectors(lam, lamp)]
    return _RECS


def stratum(w):
    for nm, lo, hi in STRATA:
        if lo <= w <= hi:
            return nm


def gstr(g):
    return '(%+d,%+d,%+d)' % g


NAMED = {'+m12': (1, 0, 0), '-m12': (-1, 0, 0), '+m13': (0, 1, 0), '-m13': (0, -1, 0), '+m23': (0, 0, 1), '-m23': (0, 0, -1), '+m33': (0, -1, -1), '-m33': (0, 1, 1),
         '+(m12+m13)': (1, 1, 0), '-(m12+m13)': (-1, -1, 0), '+(m12-m13)': (1, -1, 0), '-(m12-m13)': (-1, 1, 0), '+(m12+m13+m23)': (1, 1, 1), '-(m12+m13+m23)': (-1, -1, -1)}


def grading_stats():
    """stats[(space, g)] = Counter: N, det, n0one, worst deficit, per stratum  (det_<s>, n_<s>, n0one_<s>); band histogram"""
    stats = {}
    for space in SPACES:
        for g in GRADINGS:
            st = Counter()
            bands = Counter()
            for lam, lamp, s in all_recs():
                a = analyze_cached(s, space, g)
                sn = stratum(s.w)
                st['N'] += 1
                st['n_' + sn] += 1
                st['band%d' % a['band']] += 1
                if a['ok']:
                    st['det'] += 1
                    st['det_' + sn] += 1
                    if a['n0'] == 1:
                        st['n0one'] += 1
                        st['n0one_' + sn] += 1
                else:
                    st['worstdef'] = max(st['worstdef'], a['deficit'])
                st['band_val_max'] = max(st['band_val_max'], a['band_val'])
            stats[(space, g)] = st
    return stats


def stage_gradings():
    hdr('A1  LATTICE GRADINGS: band, level sizes, forward determinacy of the stacked inter-level blocks (exact ranks over Q)')
    recs = all_recs()
    prods = Counter()
    mzero = 0
    eq = 0
    inv_bad = 0
    cut_bad = 0
    for lam, lamp, s in recs:
        prods[prod_str(lam, lamp)] += 1
        mzero += (s.m == 0)
        eq += (len(s.in1) == s.w)
        inv_bad += (not s.w1_invariant())
        cut_bad += (not s.cut_is_kernel())
    out('weight spaces (w >= 2) of P_small + P_big: %d  (%d with m = 0; %d with W1 = W; W1 NOT invariant in %d; the cut m12 <= mu1 - lam\'_2 differs from ker E12^k1 in %d)' % (len(recs), mzero, eq, inv_bad, cut_bad))
    out('per product: ' + ', '.join('%s %d' % kv for kv in prods.items()))
    wsz = Counter(s.w for _, _, s in recs)
    out('weight-space sizes w: ' + ', '.join('%d:%d' % kv for kv in sorted(wsz.items())))
    t0 = time.perf_counter()
    stats = grading_stats()
    out('(%d directed gradings x %d spaces x %d weight spaces analysed in %.1f s; 62 grading pairs, both directions listed since g -> -g reverses the recursion)' % (len(GRADINGS), len(SPACES), len(recs), time.perf_counter() - t0))
    bound = lambda g: abs(g[0]) + max(abs(g[1]), abs(g[2]), abs(g[1] - g[2]))
    diffs = [(g, bound(g), stats[('W', g)]['band_val_max']) for g in GRADINGS if stats[('W', g)]['band_val_max'] != bound(g)]
    out('band (max |g(target) - g(source)| over the nonzero entries of the three matrices, all %d weight spaces): equals the structural bound |a| + max(|b|,|c|,|b-c|) for %d of the %d directed gradings; smaller for %d: %s' % (
        len(recs), len(GRADINGS) - len(diffs), len(GRADINGS), len(diffs), ' '.join('%s:%d<%d' % (gstr(g), d, b) for g, b, d in diffs[:8])))
    prim = [g for g in GRADINGS if math.gcd(math.gcd(abs(g[0]), abs(g[1])), abs(g[2])) == 1]
    out('%d of the 124 directed gradings are primitive (gcd 1); the others coincide with them up to scale (identical level structure, identical results): 49 genuinely different pairs' % len(prim))
    tab('A1 per-grading statistics (directed gradings g = a m12 + b m13 + c m23, ascending levels = recursion direction)')
    for space in SPACES:
        order = sorted(GRADINGS, key=lambda g: (-stats[(space, g)]['det_w5-8'] - stats[(space, g)]['det_w9+'], -stats[(space, g)]['n0one_w5-8'] - stats[(space, g)]['n0one_w9+'], -stats[(space, g)]['det'], g))
        tab('')
        tab('space %s: grading   N   det  n0=1 | det/n for w2-4, w5-8, w9+ | n0=1 for w2-4, w5-8, w9+ | max band (g-values) | worst deficit' % space)
        for g in order:
            st = stats[(space, g)]
            tab('  %s %4d %4d %4d | %s | %s | %d | %d' % (gstr(g), st['N'], st['det'], st['n0one'],
                ' '.join('%d/%d' % (st['det_' + nm], st['n_' + nm]) for nm, _, _ in STRATA), ' '.join('%d' % st['n0one_' + nm] for nm, _, _ in STRATA), st['band_val_max'], st['worstdef']))
        out('')
        out('space %s -- the 12 best directed gradings (ranked by full determinacy then n0 = 1 on the weight spaces with w >= 5; %d spaces in total, %d with w >= 5):' % (space, len(recs), sum(stats[(space, GRADINGS[0])]['n_' + nm] for nm in ('w5-8', 'w9+'))))
        out('   grading        det/all  n0=1   | w5-8: det/n  n0=1 | w9+: det/n  n0=1 | band')
        for g in order[:12]:
            st = stats[(space, g)]
            out('   %s  %4d/%d %4d   |  %3d/%3d %4d  |  %3d/%3d %4d  | %d' % (gstr(g), st['det'], st['N'], st['n0one'], st['det_w5-8'], st['n_w5-8'], st['n0one_w5-8'], st['det_w9+'], st['n_w9+'], st['n0one_w9+'], st['band_val_max']))
        out('named candidates:')
        for nm, g in NAMED.items():
            st = stats[(space, g)]
            out('   %-16s %s  det %4d/%d  n0=1 %4d  | w>=5: det %3d/%3d  n0=1 %3d  | worst deficit %d' % (nm, gstr(g), st['det'], st['N'], st['n0one'], st['det_w5-8'] + st['det_w9+'],
                st['n_w5-8'] + st['n_w9+'], st['n0one_w5-8'] + st['n0one_w9+'], st['worstdef']))
        # pairs fully determined in both directions
        both = [g for g in GRADINGS if g > tuple(-x for x in g) and stats[(space, g)]['det'] == stats[(space, g)]['N'] and stats[(space, tuple(-x for x in g))]['det'] == stats[(space, g)]['N']]
        full = [g for g in GRADINGS if stats[(space, g)]['det'] == stats[(space, g)]['N']]
        out('gradings (directed) fully determined on EVERY weight space of the list: %d  %s' % (len(full), ' '.join(gstr(g) for g in full[:12])))
    # per weight space: best gradings and n0
    tab('')
    tab('A1 per weight space: best fully determined gradings (smallest n0, then band, then |g|_1), P_small + P_big; columns: product mu  w |W1| m  A | space: n0 [gradings, up to 4]')
    minn0 = {sp: Counter() for sp in SPACES}
    nodet = {sp: 0 for sp in SPACES}
    for lam, lamp, s in recs:
        row = '%-22s %-12s w=%-2d w1=%-2d m=%-2d A=%-8s' % (prod_str(lam, lamp), mu_str(s.mu), s.w, len(s.in1), s.m, ''.join(map(str, sorted(s.active))) or '-')
        for space in SPACES:
            bg = best_gradings(s, space)
            if not bg:
                nodet[space] += 1
                row += ' | %s: NONE' % space
                continue
            n0 = bg[0][0]
            minn0[space][n0] += 1
            same = [r[3] for r in bg if r[0] == n0][:4]
            row += ' | %s: n0=%d %s' % (space, n0, ' '.join(gstr(g) for g in same))
        tab(row)
    for space in SPACES:
        out('')
        out('space %s: smallest n0 over all fully determined gradings, distribution over the %d weight spaces: %s ; no determined grading at all: %d' % (
            space, len(recs), ', '.join('n0=%d: %d' % kv for kv in sorted(minn0[space].items())), nodet[space]))
    stage_rule(stats)
    return stats


def stage_rule(stats):
    out('')
    out('RULE TEST: active strings A (string E_i^{k_i} non-void) -> natural gradings, space W1; entries det/n0=1 counts over the weight spaces of the class (w >= 2)')
    tab('')
    tab('A1 rule test: class A (active strings) x candidate grading on W1 and W: det/n  n0one')
    cls = defaultdict(list)
    for lam, lamp, s in all_recs():
        cls[s.active].append(s)
    for space in SPACES:
        for A in sorted(cls, key=lambda a: (len(a), sorted(a))):
            ss = cls[A]
            line = 'A=%-8s %s n=%3d (m=0: %d, mean w %.1f):' % (''.join(map(str, sorted(A))) or '-', space, len(ss), sum(1 for s in ss if s.m == 0), sum(s.w for s in ss) / len(ss))
            tab(line)
            res = []
            for g in GRADINGS:
                d = sum(1 for s in ss if analyze_cached(s, space, g)['ok'])
                o = sum(1 for s in ss if analyze_cached(s, space, g)['ok'] and analyze_cached(s, space, g)['n0'] == 1)
                res.append((-d, -o, g, d, o))
            res.sort()
            tab('    ' + '  '.join('%s %d/%d n0=1:%d' % (gstr(r[2]), r[3], len(ss), r[4]) for r in res[:6]))
            tab('    named: ' + '  '.join('%s %d/%d,%d' % (nm, sum(1 for s in ss if analyze_cached(s, space, g)['ok']), len(ss),
                sum(1 for s in ss if analyze_cached(s, space, g)['ok'] and analyze_cached(s, space, g)['n0'] == 1)) for nm, g in NAMED.items()))
            prim = [g for g in GRADINGS if math.gcd(math.gcd(abs(g[0]), abs(g[1])), abs(g[2])) == 1]
            fulls = [g for g in prim if all(analyze_cached(s, space, g)['ok'] for s in ss)]
            tab('    primitive directed gradings determined on ALL %d spaces of the class: %d  %s' % (len(ss), len(fulls), ' '.join(gstr(g) for g in fulls[:14])))
            if space == 'W1':
                out('      A=%-8s: %d of the 98 primitive directed gradings are determined on every W1 of the class, n0 = 1 on all: %d  %s' % (''.join(map(str, sorted(A))) or '-', len(fulls),
                    sum(1 for g in fulls if all(analyze_cached(s, space, g)['n0'] == 1 for s in ss)), ' '.join(gstr(g) for g in fulls[:6])))
            if space == 'W1':
                out('  A=%-8s n=%3d: best %s det %d n0=1 %d | +m12: %d/%d  -m12: %d/%d | +m13: %d  +m23: %d' % (''.join(map(str, sorted(A))) or '-', len(ss), gstr(res[0][2]), res[0][3], res[0][4],
                    sum(1 for s in ss if analyze_cached(s, space, (1, 0, 0))['ok']), len(ss), sum(1 for s in ss if analyze_cached(s, space, (-1, 0, 0))['ok']), len(ss),
                    sum(1 for s in ss if analyze_cached(s, space, (0, 1, 0))['ok']), sum(1 for s in ss if analyze_cached(s, space, (0, 0, 1))['ok'])))


# ---------------------------------------------------------------- Vandermonde side check
def vandermonde_sources(secs, rng, nsamples):
    """(sector, source position) pairs for which all three coherent Delta m12 = +1 root-move targets exist with nonzero V3 entry, in both orientation sets"""
    cand = []
    for s in secs:
        for c in range(s.w):
            pc = s.pats[c]
            tg = {}
            for r in range(s.w):
                pr = s.pats[r]
                if pr[1][0] - pc[1][0] != 1:
                    continue
                d = tuple(pr[2][i] - pc[2][i] for i in range(3))
                if sorted(d) != [-1, 0, 1] or s.M[0][r][c] == 0:
                    continue
                tg[(d.index(0), d)] = r
            for orient in (((0, -1, 1), (-1, 0, 1), (-1, 1, 0)), ((0, 1, -1), (1, 0, -1), (1, -1, 0))):
                if all((k, orient[k]) in tg for k in range(3)):
                    cand.append((s, c, orient, [tg[(k, orient[k])] for k in range(3)]))
    rng.shuffle(cand)
    return cand[:nsamples], len(cand)


def det3(M):
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))


def vandermonde_check(nsamples=50, seed=7):
    secs = [s for lam, lamp, s in all_recs() if s.w >= 3]
    rng = random.Random(seed)
    samples, ncand = vandermonde_sources(secs, rng, nsamples)
    ratios = Counter()
    for s, c, orient, tgt in samples:
        y = [s.pats[c][2][k] - (k + 1) for k in range(3)]
        Mx = [[s.M[t][r][c] for r in tgt] for t in range(3)]
        prodc = Mx[0][0] * Mx[0][1] * Mx[0][2]
        vd = (y[1] - y[0]) * (y[1] - y[2])
        ratios[(orient[0], det3(Mx) / (prodc * vd) if prodc * vd != 0 else 'zero')] += 1
    return ratios, len(samples), ncand


def stage_vandermonde():
    hdr('A1 side check: the Vandermonde determinant of the three Delta m12 = +1 root-move targets')
    ratios, n, ncand = vandermonde_check()
    out('%d candidate (sector, source) pairs with all three targets of a coherent orientation set (spectators 1, 2, 3) present and V3 != 0; %d sampled (seed 7).' % (ncand, n))
    out('ratio det[(V3, V4a, V4b) entries to the 3 targets] / (prod_alpha c_alpha (y2-y1)(y2-y3)), y_k = m_k3 - k of the SOURCE, per orientation set:')
    for (o, r), cnt in sorted(ratios.items(), key=lambda t: str(t)):
        out('   targets (moves of row 3) %s : ratio %s  x %d' % (o, r, cnt))
    # all mixed combos
    secs = [s for lam, lamp, s in all_recs() if s.w >= 3 and lam == lamp and lam in ((4, 2, 1, 0), (6, 4, 2, 0))]
    mixed = Counter()
    for s in secs:
        for c in range(s.w):
            pc = s.pats[c]
            tg = {}
            for r in range(s.w):
                pr = s.pats[r]
                if pr[1][0] - pc[1][0] != 1:
                    continue
                d = tuple(pr[2][i] - pc[2][i] for i in range(3))
                if sorted(d) != [-1, 0, 1] or s.M[0][r][c] == 0:
                    continue
                tg.setdefault(d.index(0), []).append((r, d))
            if len(tg) < 3:
                continue
            y = [pc[2][k] - (k + 1) for k in range(3)]
            for combo in itertools.product(*[tg[k] for k in range(3)]):
                Mx = [[s.M[t][r][c] for (r, d) in combo] for t in range(3)]
                pcc = Mx[0][0] * Mx[0][1] * Mx[0][2]
                vd = (y[1] - y[0]) * (y[1] - y[2])
                if pcc == 0 or vd == 0:
                    continue
                r_ = det3(Mx) / (pcc * vd)
                mixed[(tuple(d for _, d in combo), r_ in (1, -1))] += 1
    coh = sum(v for (k, ok), v in mixed.items() if ok)
    out('all target triples with spectators 1,2,3 over [4,2,1,0]^2 and [6,4,2,0]^2: %d with ratio exactly +-1, %d with another ratio (the mixed-orientation triples)' % (coh, sum(v for (k, ok), v in mixed.items() if not ok)))
    return ratios


# ======================================================================================== 5. A2: recursion validation
def slice_err(a, b):
    keys = set(a) | set(b)
    return max([abs(a.get(k, 0) - b.get(k, 0)) for k in keys] + [mp.mpf(0)])


def n0one_gradings(sec, space='W1', limit=8):
    return [r[3] for r in best_gradings(sec, space) if r[0] == 1][:limit]


def validate_sector(sec, space='W1', dps=50):
    """compare slice_from_triple with so.copies for every copy; returns dict"""
    gl = n0one_gradings(sec, space)
    if not gl:
        return None
    t0 = time.perf_counter()
    cps = so.copies(sec.lam, sec.lamp, sec.mu, dps=dps)
    t_eig = time.perf_counter() - t0
    res = dict(ncopies=len(cps), grading=gl[0], t_eig=t_eig, first_ok=0, any_ok=0, maxerr=mp.mpf(0), maxerr_any=mp.mpf(0), resid=mp.mpf(0), resid_pert=mp.mpf(1), t_rec=0.0, bad=0, fb=Counter())
    old = mp.mp.dps
    mp.mp.dps = dps + 10
    tol = mp.mpf(10) ** -40
    try:
        for cp in cps:
            best = None
            for gi, g in enumerate(gl):
                t1 = time.perf_counter()
                try:
                    sf = slice_from_triple(sec.lam, sec.lamp, sec.mu, cp['ev'], grading=g, space=space, dps=dps)
                except ValueError:
                    continue
                dt = time.perf_counter() - t1
                err = slice_err(sf['slice'], cp['slice'])
                if gi == 0:
                    res['t_rec'] += dt
                    res['maxerr'] = max(res['maxerr'], err)
                    res['resid'] = max(res['resid'], sf['resid'])
                    if err <= tol:
                        res['first_ok'] += 1
                if err <= tol:
                    res['any_ok'] += 1
                    res['fb'][gi] += 1
                    best = err
                    break
            if best is None:
                res['bad'] += 1
            else:
                res['maxerr_any'] = max(res['maxerr_any'], best)
            # perturbed triple
            ev = cp['ev']
            pert = tuple(x * (1 + mp.mpf(10) ** -3) + mp.mpf(10) ** -3 for x in ev)
            sf = slice_from_triple(sec.lam, sec.lamp, sec.mu, pert, grading=gl[0], space=space, dps=dps)
            res['resid_pert'] = min(res['resid_pert'], sf['resid']) if sf['resid'] != 0 else res['resid_pert']
    finally:
        mp.mp.dps = old
    return res


def stage_recursion():
    hdr('A2  slice_from_triple: forward lattice recursion validated against so.copies (50 digits)')
    prods = P_SMALL + [((6, 4, 2, 0), (6, 4, 2, 0))]
    tot = Counter()
    tab('A2 recursion validation: product mu w w1 m grading(n0=1)  copies  first-grading ok / any-grading ok  max err  resid(true triple)  resid(perturbed)  t_eig(s) t_rec/copy(ms)')
    maxerr = mp.mpf(0)
    maxresid = mp.mpf(0)
    minpert = mp.mpf(1)
    t_eig_tot = 0.0
    t_rec_tot = 0.0
    gr_used = Counter()
    for lam, lamp in prods:
        c = Counter()
        for sec in product_sectors(lam, lamp):
            if sec.m == 0:
                c['m0'] += 1
                continue
            r = validate_sector(sec)
            if r is None:
                c['no_n0_1'] += 1
                c['copies_skipped'] += sec.m
                tab('%-20s %-12s w=%-2d w1=%-2d m=%-2d  NO fully determined n0 = 1 grading' % (prod_str(lam, lamp), mu_str(sec.mu), sec.w, len(sec.in1), sec.m))
                continue
            c['sectors'] += 1
            c['copies'] += r['ncopies']
            c['first_ok'] += r['first_ok']
            c['any_ok'] += r['any_ok']
            c['bad'] += r['bad']
            gr_used[r['grading']] += 1
            maxerr = max(maxerr, r['maxerr_any'])
            maxresid = max(maxresid, r['resid'])
            minpert = min(minpert, r['resid_pert'])
            t_eig_tot += r['t_eig']
            t_rec_tot += r['t_rec']
            tab('%-20s %-12s w=%-2d w1=%-2d m=%-2d %s  %2d  %2d/%2d  err %s  resid %s  pert %s  %.3f s %.2f ms' % (prod_str(lam, lamp), mu_str(sec.mu), sec.w, len(sec.in1), sec.m, gstr(r['grading']),
                r['ncopies'], r['first_ok'], r['any_ok'], mp.nstr(r['maxerr_any'], 3), mp.nstr(r['resid'], 3), mp.nstr(r['resid_pert'], 3), r['t_eig'], 1000 * r['t_rec'] / max(1, r['ncopies'])))
        out('%-20s sectors validated %3d (m = 0: %d; no n0 = 1 grading on W1: %d [%d copies not tested]); copies %4d: reproduced with the FIRST grading %4d, with a fallback n0=1 grading %4d, NOT reproduced %d' % (
            prod_str(lam, lamp), c['sectors'], c['m0'], c['no_n0_1'], c['copies_skipped'], c['copies'], c['first_ok'], c['any_ok'], c['bad']))
        for k, v in c.items():
            tot[k] += v
    out('TOTAL: %d sectors, %d copies; reproduced (<= 1e-40 at 50 digits): %d (first grading) / %d (with fallback); NOT reproduced: %d; copies in sectors without an n0=1 grading (untested): %d' % (
        tot['sectors'], tot['copies'], tot['first_ok'], tot['any_ok'], tot['bad'], tot['copies_skipped']))
    out('max error over reproduced copies %s; max residual of the unused equations at the true triples %s; smallest residual at the perturbed triples (relative 1e-3) %s' % (
        mp.nstr(maxerr, 3), mp.nstr(maxresid, 3), mp.nstr(minpert, 3)))
    out('time: joint eigen-decomposition (so.copies, 50 digits) %.2f s total vs recursion %.3f s total over the copies (first grading)  ratio %.0fx' % (t_eig_tot, t_rec_tot, t_eig_tot / max(t_rec_tot, 1e-9)))
    out('gradings used first: ' + ', '.join('%s x%d' % (gstr(g), n) for g, n in gr_used.most_common(8)))
    return tot


# ======================================================================================== 6. A3: polynomials and characteristic equations
def mtilde_prv(sec):
    """exact basis of M~ = ker E12^k1 cap ker E23^k2 inside W (like se.prv_space with i in {1,2} only)"""
    R, lamp, nu = sec.R, sec.lamp, sec.nu
    rows = []
    for i in (1, 2):
        k = lamp[i - 1] - lamp[i] + 1
        if gr.wadd(nu, gr.dij(4, i, i + 1), k) not in R.wsp:
            continue
        imgs = se.power_images(R, (i, i + 1), k, sec.W)
        byrow = defaultdict(lambda: [Fr(0)] * sec.w)
        for a, vec in enumerate(imgs):
            for r, x in vec.items():
                byrow[r][a] = x
        rows.extend(byrow.values())
    basis, free = se._rref_nullspace(rows, sec.w)
    return se.PRV(sec.lam, sec.lamp, sec.mu, sec.W, basis, free)


def spec_count(trips, ref, tol):
    return sum(1 for t in trips if any(max(abs(a - b) for a, b in zip(t, r)) < tol for r in ref))


def charpoly_sector(sec, space='W1', dps=50, groebner=True, gtimeout=90):
    gl = [r for r in best_gradings(sec, space) if r[0] == 1]
    if not gl:
        return None
    band1 = [r for r in gl if r[1] == 1]
    g = (band1 or gl)[0][3]
    band = (band1 or gl)[0][1]
    sp = sec.space(space)
    sy = sym_recursion(sec, space, g)
    R, gens, pl = sy['R'], sy['gens'], sy['plan']
    gv, lv, li, groups = levels_of(sp, g)
    deg = {p: tdeg(v[0]) for p, v in sy['psi'].items()}
    maxdeg = max(deg.values())
    leveldeg = [max(deg[p] for p in grp) for grp in groups]
    exceeds = sum(1 for p in deg if deg[p] > li[p])
    # degrees on the PRV free coordinates and G4's Krylov degree
    prv = sec.prv
    loc = {gp: i for i, gp in enumerate(sec.in1)} if space == 'W1' else {i: i for i in range(sec.w)}
    free_loc = [loc[f] for f in prv.free if f in loc]
    degfree = max([deg[p] for p in free_loc] + [0])
    kry = None
    if prv.m >= 2:
        _, md = so.triple_mats(sec.lam, sec.lamp, sec.mu, prv)
        mats = [md[n] for n in TRIPLE]
        for x0 in range(prv.m):
            kr = sg.left_krylov(mats, x0)
            if kr['cyclic']:
                d_ = max(kr['deg'])
                kry = d_ if kry is None else min(kry, d_)
    # unused equations -> polynomials
    polys = []
    seen = set()
    for row in sy['unused']:
        p = row[0]
        if p == 0:
            continue
        q = prim_int(p)
        key = tuple(sorted(q.terms()))
        if key in seen:
            continue
        seen.add(key)
        polys.append(q)
    polys.sort(key=lambda p: (tdeg(p), len(p)))
    nuns = len(sy['unused'])
    # joint spectrum of the space
    old = mp.mp.dps
    mp.mp.dps = dps + 10
    try:
        mats_sp = [[[x for x in row] for row in sp.M[k]] for k in range(3)]
        jd = so.joint_diag(mats_sp, dps=dps)
        trips = jd['trip']
        init = pl.init[0]
        bottom_nz = [abs(jd['U'][k][init]) > mp.mpf(10) ** -20 for k in range(sp.n)]       # U columns are normalised to max component 1
        cps = so.copies(sec.lam, sec.lamp, sec.mu, dps=dps) if sec.m > 0 else []
        phys = [cp['ev'] for cp in cps]
        tol = mp.mpf(10) ** -30
        vmax_all = mp.mpf(0)
        for k, t in enumerate(trips):
            if bottom_nz[k]:
                for p in polys:
                    vmax_all = max(vmax_all, abs(peval(p, t)))
        vmax_phys = mp.mpf(0)
        for t in phys:
            for p in polys:
                vmax_phys = max(vmax_phys, abs(peval(p, t)))
        # which of the polynomials vanish at ALL W1 triples: all of those with bottom_nz.  Which are needed to cut out the spectrum of W1
        nbz = sum(1 for b in bottom_nz if not b)
        # triples of M~ and the physical ones among the W1 triples
        mt_info = None
        try:
            pt = mtilde_prv(sec)
            mt_info = dict(dim=pt.m)
            if pt.m > 0:
                R_ = sec.R
                mm = [so.vred_matrix(R_, reduced_ops(sec.lamp)[n], pt) for n in TRIPLE]
                jt = so.joint_diag(mm, dps=dps)
                mt_info['in_W1'] = spec_count(jt['trip'], trips, tol)
                mt_info['phys_in_Mt'] = spec_count(phys, jt['trip'], tol)
                # W1 triples that are not in M~
                mt_info['W1_notin_Mt'] = sp.n - spec_count(trips, jt['trip'], tol)
                # does some C_j fail to vanish at a non-M~ triple of W1 (bottom comp nonzero)? by construction no
        except Exception as e:
            mt_info = dict(error=str(e)[:80])
        dist_trips = 0
        for i, t in enumerate(trips):
            if all(max(abs(a - b) for a, b in zip(t, trips[j])) > tol for j in range(i)):
                dist_trips += 1
    finally:
        mp.mp.dps = old
    gb = None
    if groebner and polys:
        gb = groebner_dim(polys, timeout=gtimeout)
    return dict(grading=g, band=band, n=sp.n, m=sec.m, maxdeg=maxdeg, leveldeg=leveldeg, D=len(lv) - 1, exceeds=exceeds, degfree=degfree, kry=kry, dgen=dgen(sec.m), dgen_n=dgen(sp.n),
                nunused=nuns, npoly=len(polys), polydeg=[tdeg(p) for p in polys], vmax_all=vmax_all, vmax_phys=vmax_phys, n_bottom_zero=nbz, dist_trips=dist_trips, gb=gb, mt=mt_info,
                jd_ok=jd['ok'], polys=polys, sy=sy)


def stage_charpoly(groebner=True, gtimeout=90):
    hdr('A3  polynomials P_x(v3, v4a, v4b) and the characteristic equations (n0 = 1 sectors, space W1, band-1 gradings where available)')
    sel = []
    for lam, lamp in P_SMALL:
        cand = [s for s in product_sectors(lam, lamp) if s.m >= 3]
        cand.sort(key=lambda s: -s.m)
        sel += cand[:3]
    sel = [s for s in sel if best_gradings(s, 'W1') and best_gradings(s, 'W1')[0][0] == 1]
    sel.sort(key=lambda s: (-s.m, s.w))
    seen = set()
    sel2 = []
    for s in sel:
        key = (s.m, s.w, len(s.in1))
        if key in seen and len(sel2) >= 8:
            continue
        seen.add(key)
        sel2.append(s)
    sel = sel2[:14]
    big = [s for s in product_sectors((6, 4, 2, 0), (6, 4, 2, 0)) if s.m >= 6 and best_gradings(s, 'W1') and best_gradings(s, 'W1')[0][0] == 1]
    big.sort(key=lambda s: (-s.m, s.w))
    sel += big[:3]
    out('%d sectors selected (>= 3 with m >= 3 from P_small; 3 of [6,4,2,0]^2 with m >= 6)' % len(sel))
    tab('A3 per sector: product mu | w1=dim W1  m | grading band | max deg P_x, max deg per level | #x with deg>level | max deg on PRV coords, G4 Krylov deg, d_gen(m) | #unused eqs, #distinct nonzero C_j, their degrees | max|C_j| at W1 triples (bottom comp != 0), at physical triples | groebner dim Q[v]/J vs dim W1 vs #distinct W1 triples | M~')
    agg = Counter()
    for sec in sel:
        t0 = time.perf_counter()
        r = charpoly_sector(sec, groebner=groebner, gtimeout=gtimeout)
        if r is None:
            continue
        agg['sectors'] += 1
        agg['vanish_all'] += (r['vmax_all'] < mp.mpf(10) ** -30)
        agg['vanish_phys'] += (r['vmax_phys'] < mp.mpf(10) ** -30)
        agg['deg_le_level'] += (r['exceeds'] == 0)
        agg['krylov_equal'] += (r['kry'] is not None and r['degfree'] == r['kry'])
        agg['krylov_le'] += (r['kry'] is not None and r['degfree'] <= r['kry'])
        agg['band1'] += (r['band'] == 1)
        gbd = r['gb']['dim'] if r['gb'] else None
        if gbd is not None:
            agg['gb_done'] += 1
            agg['gb_eq_n'] += (gbd == r['n'])
            agg['gb_eq_dist'] += (gbd == r['dist_trips'])
        tab('%-18s %-12s | w1=%-2d m=%-2d | %s B=%d | maxdeg %d %s | exceeds-level %d | PRV-coord deg %d, krylov %s, dgen %d (dgen(w1) %d) | unused %d, C_j %d %s | vmax all %s phys %s (bottom comp = 0 at %d triples) | gb %s vs w1 %d vs distinct %d | M~ %s | %.1fs' % (
            prod_str(sec.lam, sec.lamp), mu_str(sec.mu), r['n'], r['m'], gstr(r['grading']), r['band'], r['maxdeg'], r['leveldeg'], r['exceeds'], r['degfree'], r['kry'], r['dgen'], r['dgen_n'],
            r['nunused'], r['npoly'], r['polydeg'], mp.nstr(r['vmax_all'], 3), mp.nstr(r['vmax_phys'], 3), r['n_bottom_zero'], r['gb'], r['n'], r['dist_trips'], r['mt'], time.perf_counter() - t0))
        out('%-18s %-12s w1=%-2d m=%-2d %s: deg P_x <= %d (levels 0..%d: %s), C_j: %d nonzero of %d unused, degrees %s; |C_j| at W1 triples %s, at phys %s; groebner dim %s (w1 = %d, distinct triples %d); M~ %s' % (
            prod_str(sec.lam, sec.lamp), mu_str(sec.mu), r['n'], r['m'], gstr(r['grading']), r['maxdeg'], r['D'], r['leveldeg'], r['npoly'], r['nunused'], sorted(set(r['polydeg'])),
            mp.nstr(r['vmax_all'], 2), mp.nstr(r['vmax_phys'], 2), r['gb'].get('dim') if r['gb'] else None, r['n'], r['dist_trips'], r['mt']))
        if sec is sel[0]:
            out('   example C_j polynomials (first sector), primitive integer form:')
            for p in r['polys'][:8]:
                out('      ' + str(p.as_expr()))
    out('AGGREGATE over %d sectors: all C_j vanish at every W1 joint triple (bottom component != 0) in %d, at all physical triples in %d; deg P_x <= level in %d; PRV-coordinate degree = G4 Krylov degree in %d (<= in %d); band-1 gradings %d; groebner finished %d, dim Q[v]/J = dim W1 in %d, = #distinct triples in %d' % (
        agg['sectors'], agg['vanish_all'], agg['vanish_phys'], agg['deg_le_level'], agg['krylov_equal'], agg['krylov_le'], agg['band1'], agg['gb_done'], agg['gb_eq_n'], agg['gb_eq_dist']))
    return agg


# ======================================================================================== 7. A4: n0 > 1
def min_n0_table():
    """per (space): rows (w, m, A) -> Counter of the minimal n0 over fully determined gradings"""
    res = {sp_: defaultdict(Counter) for sp_ in SPACES}
    for lam, lamp, s in all_recs():
        for space in SPACES:
            bg = best_gradings(s, space)
            res[space][(s.w, s.m, ''.join(map(str, sorted(s.active))) or '-')][bg[0][0] if bg else None] += 1
    return res


def check_n0_sector(sec, space='W1', dps=50):
    bg = best_gradings(sec, space)
    n0, band, _, g = bg[0]
    if n0 < 2:
        return None
    pl = get_plan(sec, space, g)
    sp = sec.space(space)
    cps = so.copies(sec.lam, sec.lamp, sec.mu, dps=dps)
    old = mp.mp.dps
    mp.mp.dps = dps + 10
    res = dict(grading=g, n0=n0, band=band, ncopies=len(cps), ok=0, rank_ok=0, maxker=mp.mpf(0), maxfull=mp.mpf(0), minsig_pert=mp.mpf(1), maxsig=mp.mpf(0), skipped=0)
    try:
        ident = [[mp.mpf(int(i == j)) for j in range(n0)] for i in range(n0)]
        for cp in cps:
            ev = cp['ev']
            b = [cp['slice'].get(int(sec.W[sp.pos[p]]), mp.mpf(0)) for p in pl.init]
            nb = mp.sqrt(sum(x * x for x in b))
            psi, unused = run_plan(pl, [mp.mpf(x) for x in ev], ident, mpq, mp.mpf(0), ('mp', dps + 10))
            C = mp.matrix(len(unused), n0)
            for i, row in enumerate(unused):
                for j in range(n0):
                    C[i, j] = row[j]
            U_, S, V = mp.svd_r(C)
            s1 = S[0]
            res['maxsig'] = max(res['maxsig'], S[n0 - 1] / s1)
            rank_ok = (S[n0 - 1] / s1 < mp.mpf(10) ** -40) and (n0 == 1 or S[n0 - 2] / s1 > mp.mpf(10) ** -15)
            res['rank_ok'] += rank_ok
            if nb < mp.mpf(10) ** -30:
                res['skipped'] += 1
                continue
            ker = [V[n0 - 1, j] for j in range(n0)]
            nk = mp.sqrt(sum(x * x for x in ker))
            cosv = abs(sum(x * y for x, y in zip(ker, b))) / (nk * nb)
            res['maxker'] = max(res['maxker'], 1 - cosv)
            # reconstruct the full vector from the kernel and compare with the copy's slice (scaled)
            sc = [sum(ker[j] * psi[p][j] for j in range(n0)) for p in range(sp.n)]
            vec = {int(sec.W[sp.pos[p]]): sc[p] for p in range(sp.n)}
            # scale: match at the first component of the copy
            k0 = min(cp['slice'])
            f = cp['slice'][k0] / vec[k0] if vec.get(k0, 0) != 0 else None
            if f is not None:
                err = slice_err({k: f * x for k, x in vec.items() if abs(f * x) > mp.mpf(10) ** -35}, cp['slice'])
                res['maxfull'] = max(res['maxfull'], err)
            res['ok'] += 1
            # a perturbed triple: C should have full rank n0
            pert = [x * (1 + mp.mpf(10) ** -3) + mp.mpf(10) ** -3 for x in ev]
            psi2, un2 = run_plan(pl, pert, ident, mpq, mp.mpf(0), ('mp', dps + 10))
            C2 = mp.matrix(len(un2), n0)
            for i, row in enumerate(un2):
                for j in range(n0):
                    C2[i, j] = row[j]
            S2 = mp.svd_r(C2, compute_uv=False)
            res['minsig_pert'] = min(res['minsig_pert'], S2[n0 - 1] / S2[0])
    finally:
        mp.mp.dps = old
    return res


def stage_n0(gtimeout=60):
    hdr('A4  n0 > 1: the recursion psi = T(v) psi0, characteristic equations C(v) psi0 = 0 with C a polynomial matrix; minimal n0 over gradings')
    tab_min = min_n0_table()
    tab('A4 minimal n0 over fully determined gradings, by (w, m, A): space -> {n0: number of weight spaces}')
    for space in SPACES:
        tab('space ' + space)
        for key in sorted(tab_min[space]):
            tab('   w=%-2d m=%-2d A=%-4s %s' % (key[0], key[1], key[2], dict(tab_min[space][key])))
    for space in SPACES:
        byw = defaultdict(Counter)
        for key, cnt in tab_min[space].items():
            for n0, k in cnt.items():
                byw[key[0]][n0] += k
        out('space %s: minimal n0 by weight-space dimension w: ' % space + '; '.join('w=%d: %s' % (w, dict(c)) for w, c in sorted(byw.items())))
        # relation n0 = 1 vs A
        byA = defaultdict(Counter)
        for key, cnt in tab_min[space].items():
            for n0, k in cnt.items():
                byA[key[2]][n0] += k
        out('space %s: minimal n0 by active strings A: ' % space + '; '.join('A=%s %s' % (a, dict(c)) for a, c in sorted(byA.items())))
    # products: a single grading with n0 = 1 everywhere?
    out('')
    out('Per product: the grading with the largest number of sectors (w >= 2) fully determined with n0 = 1 (space W1 / W):')
    tab('')
    tab('A4 per product: best single grading for n0 = 1 (det and n0=1 counts)')
    for lam, lamp in P_SMALL + P_BIG:
        ss = product_sectors(lam, lamp)
        row = '%-22s %3d sectors:' % (prod_str(lam, lamp), len(ss))
        for space in SPACES:
            best = None
            for g in GRADINGS:
                o = sum(1 for s in ss if analyze_cached(s, space, g)['ok'] and analyze_cached(s, space, g)['n0'] == 1)
                d = sum(1 for s in ss if analyze_cached(s, space, g)['ok'])
                if best is None or (o, d) > best[:2]:
                    best = (o, d, g)
            nmin1 = sum(1 for s in ss if best_gradings(s, space)[0][0] == 1)
            row += '  %s: %s n0=1 in %d (det %d); some grading gives n0=1 in %d' % (space, gstr(best[2]), best[0], best[1], nmin1)
        out('   ' + row)
        tab(row)
    # numeric verification
    out('')
    out('Numerical verification on sectors whose minimal n0 over the determined gradings (W1) is 2 or 3 (50 digits):')
    secs = []
    for lam, lamp in P_SMALL + P_BIG:
        for s in product_sectors(lam, lamp):
            if s.m == 0 or s.w > 15:
                continue
            bg = best_gradings(s, 'W1')
            if bg and bg[0][0] >= 2:
                secs.append(s)
    out('   %d candidate sectors (m >= 1)' % len(secs))
    tab('')
    tab('A4 verification: product mu w1 m n0 grading | copies, rank(C)=n0-1 at the copy triple, 1-cos(ker C, bottom comp) max, max error of T(v)ker vs copy, smallest sigma_n0/sigma_1 at a perturbed triple')
    tot = Counter()
    for s in secs:
        r = check_n0_sector(s)
        if r is None:
            continue
        tot['sectors'] += 1
        tot['n0_%d' % r['n0']] += 1
        tot['copies'] += r['ncopies']
        tot['rank_ok'] += r['rank_ok']
        tot['ok'] += r['ok']
        tot['skipped'] += r['skipped']
        tot['maxker'] = max(tot['maxker'], float(r['maxker']))
        tot['maxfull'] = max(tot['maxfull'], float(r['maxfull']))
        tab('%-20s %-12s w1=%-2d m=%-2d n0=%d %s B=%d | %d copies, rank n0-1 at %d, kernel check %d (bottom comp 0 at %d), 1-cos %s, full err %s, pert sigma %s' % (
            prod_str(s.lam, s.lamp), mu_str(s.mu), len(s.in1), s.m, r['n0'], gstr(r['grading']), r['band'], r['ncopies'], r['rank_ok'], r['ok'], r['skipped'], mp.nstr(r['maxker'], 3),
            mp.nstr(r['maxfull'], 3), mp.nstr(r['minsig_pert'], 3)))
    out('   verified: %d sectors (n0 = 2: %d, n0 = 3: %d), %d copies: rank C(v) = n0 - 1 at %d of them (sigma_n0/sigma_1 < 1e-40, sigma_{n0-1}/sigma_1 > 1e-15); kernel = bottom component of the copy (1 - cos <= %.1e) and T(v)ker = the copy (error <= %.1e) checked at %d (bottom component zero at %d)' % (
        tot['sectors'], tot['n0_2'], tot['n0_3'], tot['copies'], tot['rank_ok'], tot['maxker'], tot['maxfull'], tot['ok'], tot['skipped']))
    # symbolic C(v) for n0 = 2 sectors
    out('')
    out('Symbolic C(v) (polynomial matrix, n0 = 2) and its 2 x 2 minors:')
    nsym = 0
    for s in secs:
        bg = best_gradings(s, 'W1')
        if bg[0][0] != 2 or s.w > 8 or nsym >= 5:
            continue
        g = bg[0][3]
        sy = sym_recursion(s, 'W1', g)
        rows = [r for r in sy['unused'] if any(p != 0 for p in r)]
        minors = []
        seen = set()
        for i in range(len(rows)):
            for j in range(i + 1, len(rows)):
                p = rows[i][0] * rows[j][1] - rows[i][1] * rows[j][0]
                if p != 0:
                    q = prim_int(p)
                    key = tuple(sorted(q.terms()))
                    if key not in seen:
                        seen.add(key)
                        minors.append(q)
        minors.sort(key=lambda p: (tdeg(p), len(p)))
        degC = sorted(set(tdeg(p) for r in rows for p in r if p != 0))
        gb = groebner_dim(minors[:40], timeout=gtimeout) if minors else None
        sp = s.space('W1')
        out('   %s %s w1=%d m=%d g=%s: C has %d nonzero rows, entry degrees %s, %d distinct nonzero 2x2 minors (degrees %s); dim Q[v]/(minors) = %s' % (
            prod_str(s.lam, s.lamp), mu_str(s.mu), sp.n, s.m, gstr(g), len(rows), degC, len(minors), sorted(set(tdeg(p) for p in minors)), gb))
        nsym += 1
    return tot


# ======================================================================================== 8. A5: separability
def matrix_rank_mp(rows, thr=mp.mpf(10) ** -30):
    if not rows or not rows[0]:
        return 0
    A = mp.matrix(rows)
    S = mp.svd_r(A, compute_uv=False)
    s1 = S[0]
    if s1 == 0:
        return 0
    return sum(1 for k in range(len(S)) if S[k] / s1 > thr)


def sep_ranks(sec, cp):
    """ranks of Psi[kappa][m12], Psi[m13][(m23,m12)], Psi[m23][(m13,m12)] (full-range rows, all rows with missing -> 0)"""
    pats = {int(g): sec.pats[i] for i, g in enumerate(sec.W)}
    sl = cp['slice']
    allp = list(pats.values())
    m12s = sorted({p[1][0] for p in allp})
    kap = sorted({p[2] for p in allp})
    have = {(p[2], p[1][0]) for p in allp}
    val = {}
    for g, v in sl.items():
        p = pats[g]
        val[(p[2], p[1][0])] = v
    out_ = {}
    mat = [[val.get((k, a), mp.mpf(0)) for a in m12s] for k in kap]
    full = [[val.get((k, a), mp.mpf(0)) for a in m12s] for k in kap if all((k, a) in have for a in m12s)]
    def nzrc(M):
        return (sum(1 for r in M if any(abs(x) > mp.mpf(10) ** -30 for x in r)), sum(1 for j in range(len(M[0])) if any(abs(r[j]) > mp.mpf(10) ** -30 for r in M)) if M else 0)
    out_['kappa'] = (matrix_rank_mp(full) if full else 0, matrix_rank_mp(mat), len(full), len(kap), len(m12s)) + nzrc(mat)
    for name, idx, other in (('m13', 0, 1), ('m23', 1, 0)):
        rowsv = sorted({p[2][idx] for p in allp})
        cols = sorted({(p[2][other], p[1][0]) for p in allp})
        have2 = {(p[2][idx], (p[2][other], p[1][0])) for p in allp}
        val2 = {}
        for g, v in sl.items():
            p = pats[g]
            val2[(p[2][idx], (p[2][other], p[1][0]))] = v
        mat2 = [[val2.get((r, c), mp.mpf(0)) for c in cols] for r in rowsv]
        full2 = [[val2.get((r, c), mp.mpf(0)) for c in cols] for r in rowsv if all((r, c) in have2 for c in cols)]
        out_[name] = (matrix_rank_mp(full2) if full2 else 0, matrix_rank_mp(mat2), len(full2), len(rowsv), len(cols)) + nzrc(mat2)
    return out_


def stage_separability(dps=50):
    hdr('A5  separability diagnostic: rank of Psi[kappa][m12], Psi[m13][(m23,m12)], Psi[m23][(m13,m12)] for every physical copy of P_small (50 digits, threshold 1e-30)')
    old = mp.mp.dps
    mp.mp.dps = dps + 10
    try:
        dist = {k: Counter() for k in ('kappa', 'm13', 'm23')}
        nontriv = {k: Counter() for k in ('kappa', 'm13', 'm23')}
        ncop = 0
        tab('A5 per sector: product mu m | for kappa, m13, m23: rank over full-range rows / over all rows (shape full-rows x cols; all rows x cols), listed per copy')
        for lam, lamp in P_SMALL:
            for sec in product_sectors(lam, lamp, min_w=1):
                if sec.m == 0:
                    continue
                cps = so.copies(sec.lam, sec.lamp, sec.mu, dps=dps)
                line = []
                for cp in cps:
                    r = sep_ranks(sec, cp)
                    ncop += 1
                    for k in r:
                        dist[k][(r[k][0], r[k][1])] += 1
                        if r[k][5] >= 2 and r[k][6] >= 2:
                            nontriv[k][(r[k][0], r[k][1])] += 1
                    line.append('%d/%d,%d/%d,%d/%d' % (r['kappa'][0], r['kappa'][1], r['m13'][0], r['m13'][1], r['m23'][0], r['m23'][1]))
                tab('%-20s %-12s m=%-2d  %s' % (prod_str(lam, lamp), mu_str(sec.mu), sec.m, ' | '.join(line)))
        out('%d physical copies of P_small' % ncop)
        for k, nm in (('kappa', 'Psi[kappa][m12]'), ('m13', 'Psi[m13][(m23,m12)]'), ('m23', 'Psi[m23][(m13,m12)]')):
            out('%s: (rank over full-range rows, rank over all rows) distribution over copies: %s' % (nm, dict(sorted(dist[k].items()))))
            out('      among copies whose matrix has >= 2 nonzero rows and >= 2 nonzero columns: %s' % dict(sorted(nontriv[k].items())))
        return dist
    finally:
        mp.mp.dps = old


# ======================================================================================== CLI
STAGES = {'gradings': lambda: (stage_gradings(), stage_vandermonde()), 'recursion': stage_recursion, 'charpoly': stage_charpoly, 'n0': stage_n0, 'separability': stage_separability}
ORDER = ['gradings', 'recursion', 'charpoly', 'n0', 'separability']


def main(argv):
    mp.mp.dps = 60
    cmd = argv[1] if len(argv) > 1 else 'all'
    if cmd == 'all':
        stages = ORDER
    elif cmd == 'run':
        stages = argv[2].split(',')
    else:
        stages = [cmd]
    first = True
    for st in stages:
        t0 = time.perf_counter()
        STAGES[st]()
        out('[stage %s: %.1f s]' % (st, time.perf_counter() - t0))
        flush_tables('w' if (first and cmd in ('all',)) else 'a', title='stage ' + st)
        first = False
    if cmd == 'all':
        with open(OUTPUT_PATH, 'w', encoding='utf-8') as f:
            f.write('\n'.join(_OUT) + '\n')


if __name__ == '__main__':
    main(sys.argv)
