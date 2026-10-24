"""nested_gl3_gl4.py -- WP-G8b of gl4-jacobi-plan.md: nesting gl(3) in gl(4).  The weight space W = V_lam[nu] of gl(4) is a disjoint union of gl(3)-type paths P_kappa (one per row 3 of the
GT pattern); test whether the within-path block of the Gaudin triple is the gl(3) Jacobi chain (B1), build M~ = ker E12^k1 cap ker E23^k2 = (+)_kappa (gl(3) PRV slices) and the nested basis
(B2), solve the depth-1 sectors (B3).

CLI (run from this directory):  python -X utf8 nested_gl3_gl4.py ident | within | nested | depth | all      and      run <stage,...>
Everything exact (Fractions / exact rref) unless '50 digits' is said.  Read-only imports of G1-G5 (glN_reps, slice_engine_gl4, slice_ops_gl4, structure_gl4).
"""
import os, sys, time, itertools
from fractions import Fraction as Fr
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mpmath as mp
import glN_reps as gr
import slice_engine_gl4 as se
import slice_ops_gl4 as so
import structure_gl4 as st

TABLES_PATH = os.path.join(HERE, 'nested_gl3_gl4_tables.txt')
DATE = '2026-10-24'
TRIPLE = so.TRIPLE
_TABLE_BUF = []

P_SMALL = [((2, 1, 0, 0),) * 2, ((2, 1, 1, 0),) * 2, ((3, 2, 1, 0),) * 2, ((4, 2, 0, 0),) * 2, ((4, 2, 1, 0),) * 2, ((3, 2, 1, 0), (4, 2, 1, 0)), ((4, 2, 1, 0), (3, 2, 1, 0)),
           ((4, 2, 2, 0),) * 2, ((4, 3, 1, 0),) * 2]
P_BIG = [((6, 4, 2, 0),) * 2, ((5, 3, 1, 0),) * 2, ((5, 3, 2, 1),) * 2]


def tab(s=''):
    _TABLE_BUF.append(s)
    print(s)


def hdr(s):
    tab('')
    tab('=' * 8 + ' ' + s)


def flush_tables(path=TABLES_PATH, mode='w', title=None):
    with open(path, mode, encoding='utf-8') as f:
        if title:
            f.write(title + '\n')
        f.write('\n'.join(_TABLE_BUF) + '\n')
    _TABLE_BUF.clear()


def lstr(l):
    return '[' + ','.join(str(x) for x in l) + ']'


def pstr(lam, lamp):
    return lstr(lam) + ('^2' if tuple(lam) == tuple(lamp) else 'x' + lstr(lamp))


def fr(x):
    return str(x)


# ============================================================================================ paths
def sub_mu(mu, lamp):
    return tuple(mu[i] - lamp[i] for i in range(len(mu)))


def sectors(lam, lamp):
    """all dominant mu = lam + w (w a weight of V_lamp) with non-empty weight space V_lam[mu - lamp], with Racah-Speiser multiplicity m (0 = formal)"""
    R = gr.GLNRep.of(lam)
    rs = gr.racah_speiser(lam, lamp)
    out = []
    for mu in gr.candidates(lam, lamp):
        nu = sub_mu(mu, lamp)
        if nu in R.wsp:
            out.append((mu, rs.get(mu, 0), len(R.wsp[nu])))
    return out


def paths_of(R, W):
    """{kappa: [global indices ordered by m12]} for the weight space W (array of global indices) of the gl(4) rep R"""
    P = defaultdict(list)
    for c in W:
        p = R.pats[int(c)]
        P[p[2]].append(int(c))
    for k in P:
        P[k].sort(key=lambda c: R.pats[c][1][0])
    return dict(P)


def r3_index(R3, pat4):
    """index in the gl(3) rep R3 = GLNRep.of(kappa) of the gl(4) pattern (rows 1..3 = rows[:3], row 3 = kappa)"""
    return R3.index[tuple(pat4[:3])]


def path_map(R, R3, path):
    return [r3_index(R3, R.pats[c]) for c in path]


# ============================================================================================ B0
def ident_check(lam, verbose=True):
    """B0: for each kappa the generators E12,E21,E23,E32 (and the diagonal E11,E22,E33) restricted to the patterns with row 3 = kappa equal the gl(3) GLNRep of kappa EXACTLY, and
    never leave the row-3 = kappa set.  Returns dict(kappas, entries checked, violations)"""
    R = gr.GLNRep.of(lam)
    by = defaultdict(list)
    for c, p in enumerate(R.pats):
        by[p[2]].append(c)
    nent = 0
    bad = 0
    leave = 0
    for kappa, cols in by.items():
        R3 = gr.GLNRep.of(kappa)
        for ij in [(1, 2), (2, 1), (2, 3), (3, 2), (1, 1), (2, 2), (3, 3), (1, 3), (3, 1)]:
            E4, E3 = R.E[ij], R3.E[ij]
            for c in cols:
                img4 = {}
                for r, v in E4[c].items():
                    if R.pats[r][2] != kappa:
                        leave += 1
                    else:
                        img4[r3_index(R3, R.pats[r])] = v
                img3 = E3[r3_index(R3, R.pats[c])]
                nent += len(img3) + len(img4)
                if img3 != img4:
                    bad += 1
    if verbose:
        tab('  ident %-12s dim %5d  kappas %3d  columns-generators compared (matrix entries %7d)  mismatching columns %d  entries leaving the kappa-block %d' % (lstr(lam), R.dim, len(by), nent, bad, leave))
    return dict(kappas=len(by), entries=nent, bad=bad, leave=leave)


def stage_ident():
    hdr('B0: the kappa-block of the gl(4) GT generators equals the gl(3) GLNRep of kappa (E12,E21,E23,E32 and the Cartan E11,E22,E33 exact, plus E13,E31 = commutators; none leaves the kappa-block)')
    for lam in [(2, 1, 0, 0), (3, 2, 1, 0), (4, 2, 1, 0), (5, 3, 2, 1)]:
        ident_check(lam)


# ============================================================================================ operators on a weight space
_REP_CACHE = {}


def triple_on_W(lam, lamp, mu):
    """(R, W, paths, {name: exact w x w matrix on W in the order of W}), cached"""
    key = (tuple(lam), tuple(lamp), tuple(mu))
    if key not in _REP_CACHE:
        R = gr.GLNRep.of(lam)
        W = R.wsp[sub_mu(mu, lamp)]
        ops = {n: so.reduced_op(n, lamp) for n in TRIPLE}
        M = so.weight_space_mats(R, ops, W)
        _REP_CACHE[key] = (R, W, paths_of(R, W), M)
    return _REP_CACHE[key]


_V3_CACHE = {}


def gl3_V(kappa, lamp3, nu3):
    """(R3, W3, exact matrix of the gl(3) reduced operator 'V' = t S_{3,2} (lamp3) on V_kappa[nu3]), cached"""
    key = (tuple(kappa), tuple(lamp3), tuple(nu3))
    if key not in _V3_CACHE:
        R3 = gr.GLNRep.of(kappa)
        W3 = R3.wsp.get(tuple(nu3))
        M = so.weight_space_mats(R3, {'V': so.reduced_op('V', tuple(lamp3))}, W3)['V'] if W3 is not None else None
        _V3_CACHE[key] = (R3, W3, M)
    return _V3_CACHE[key]


def block(M, pos, rows, cols):
    return [[M[pos[r]][pos[c]] for c in cols] for r in rows]


def gl3_block_for_path(lam, lamp, mu, kappa):
    """the gl(3) operator V^(3)_kappa on V_kappa[nu123] REORDERED as the path P_kappa of the gl(4) weight space (same basis, identified exactly)"""
    R, W, P, M = triple_on_W(lam, lamp, mu)
    nu = sub_mu(mu, lamp)
    R3, W3, M3 = gl3_V(kappa, lamp[:3], nu[:3])
    pos3 = {int(c): a for a, c in enumerate(W3)}
    idx = [pos3[i] for i in path_map(R, R3, P[kappa])]
    return [[M3[r][c] for c in idx] for r in idx]


# ============================================================================================ B1
def tridiag_violations(Jfull, path_idx_in_W):
    """number of non-zero entries of the within-path block with |Delta m12| >= 2 (path ordered by m12, consecutive m12)"""
    n = len(path_idx_in_W)
    bad = 0
    for a in range(n):
        for b in range(n):
            if abs(a - b) >= 2 and Jfull[a][b] != 0:
                bad += 1
    return bad


def fit_ansatz(J, V, D=None):
    """exact test of J = a V + b 1 [+ c D].  Off-diagonal: all non-zero entries of V give the same ratio a (and J = 0 wherever V = 0); then the diagonal J_ii - a V_ii must be of the form
    b + c D_ii.  Returns dict(ok1, ok2, a, b, c) (ok1 = the 1-parameter ansatz, ok2 = with D); a None when undetermined"""
    n = len(J)
    ratios = set()
    pat_ok = True
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if V[i][j] == 0:
                if J[i][j] != 0:
                    pat_ok = False
            else:
                ratios.add(J[i][j] / V[i][j])
    out = dict(ok1=False, ok2=False, a=None, b=None, c=None, pat=pat_ok, nratios=len(ratios))
    if not pat_ok or len(ratios) > 1:
        return out
    a = next(iter(ratios)) if ratios else Fr(0)
    diff = [J[i][i] - a * V[i][i] for i in range(n)]
    out['a'] = a
    if len(set(diff)) == 1:
        out.update(ok1=True, ok2=True, b=diff[0], c=Fr(0))
        return out
    if D is not None and len(set(D)) > 1:
        # diff_i = b + c D_i: solve from the first two distinct D values, verify the rest
        i0 = 0
        i1 = next(i for i in range(n) if D[i] != D[0])
        c = (diff[i1] - diff[i0]) / (D[i1] - D[i0])
        b = diff[i0] - c * D[i0]
        if all(diff[i] == b + c * D[i] for i in range(n)):
            out.update(ok2=True, b=b, c=c)
    return out


def product_test(J, V):
    """(iii) the off-diagonal products J[k+1][k] J[k][k+1] = a^2 V[k+1][k] V[k][k+1] with ONE a^2 (ratio of the products); returns (ok, a2, is_square)"""
    n = len(J)
    rs = set()
    ok = True
    for k in range(n - 1):
        pj = J[k + 1][k] * J[k][k + 1]
        pv = V[k + 1][k] * V[k][k + 1]
        if pv == 0:
            if pj != 0:
                ok = False
        else:
            rs.add(pj / pv)
    if not ok or len(rs) > 1:
        return False, None, None
    a2 = next(iter(rs)) if rs else None
    sq = None
    if a2 is not None:
        import math
        if a2 >= 0:
            num, den = a2.numerator, a2.denominator
            sq = (math.isqrt(num) ** 2 == num and math.isqrt(den) ** 2 == den)
        else:
            sq = False
    return True, a2, sq


def span_fit(Jt, bases):
    """exact: is the matrix Jt (list of rows) a linear combination of the matrices in `bases`?  returns coefficients or None"""
    n = len(Jt)
    rows = []
    for i in range(n):
        for j in range(n):
            rows.append([B[i][j] for B in bases] + [Jt[i][j]])
    # solve rows * x = rhs by exact elimination
    ncol = len(bases)
    M = [list(r) for r in rows]
    piv = []
    r = 0
    for c in range(ncol):
        p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        inv = 1 / M[r][c]
        M[r] = [x * inv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        piv.append(c)
        r += 1
    for i in range(r, len(M)):
        if M[i][ncol] != 0:
            return None
    x = [Fr(0)] * ncol
    for k, c in enumerate(piv):
        x[c] = M[k][ncol]
    return x


def within_data(lam, lamp, mu, with_gl3=True):
    """per kappa (with |P_kappa| >= 2): dict(n, J = {name: block}, V3 = gl(3) block, D, tri = {name: violations}) for one weight space"""
    R, W, P, M = triple_on_W(lam, lamp, mu)
    pos = {int(c): a for a, c in enumerate(W)}
    out = {}
    for kappa, path in P.items():
        n = len(path)
        J = {nm: block(M[nm], pos, path, path) for nm in TRIPLE}
        D = [(R.pats[c][1][0] - R.pats[c][1][1]) for c in path]
        D = [Fr(x * (x + 2)) for x in D]
        d = dict(n=n, J=J, D=D, tri={nm: tridiag_violations(J[nm], path) for nm in TRIPLE}, path=path)
        if with_gl3 and n >= 1:
            d['V3'] = gl3_block_for_path(lam, lamp, mu, kappa)
        out[kappa] = d
    return out


def stage_within_product(lam, lamp, stats, examples, verbose=False):
    for mu, m, w in sectors(lam, lamp):
        if w < 2:
            continue
        data = within_data(lam, lamp, mu)
        for kappa, d in data.items():
            n = d['n']
            for nm in TRIPLE:
                stats['tri_bad'] += d['tri'][nm]
                stats['tri_entries'] += 1
            if n < 2:
                continue
            for nm in TRIPLE:
                f = fit_ansatz(d['J'][nm], d['V3'], d['D'])
                k = (nm,)
                stats['pairs2'][k] += 1
                stats['i'][k] += f['ok1']
                stats['ii'][k] += f['ok2']
                po, a2, sq = product_test(d['J'][nm], d['V3'])
                stats['iii'][k] += po
                if n >= 3:
                    stats['pairs3'][k] += 1
                    stats['i3'][k] += f['ok1']
                    stats['ii3'][k] += f['ok2']
                    stats['iii3'][k] += po
                if not f['ok2']:
                    ex = examples.setdefault(nm, [])
                    ex.append((n, lam, lamp, mu, kappa))
            # J^(4a) in span(J3, 1), J^(4b) in span(J3, 1, D)
            n_ = n
            one = [[Fr(int(i == j)) for j in range(n_)] for i in range(n_)]
            Dm = [[d['D'][i] if i == j else Fr(0) for j in range(n_)] for i in range(n_)]
            for nm, bases in (('V4a', [d['J']['V3'], one]), ('V4a+D', [d['J']['V3'], one, Dm]), ('V4b', [d['J']['V3'], one]), ('V4b+D', [d['J']['V3'], one, Dm])):
                nm0 = nm.split('+')[0]
                stats['span'][nm] += 1
                if span_fit(d['J'][nm0], bases) is not None:
                    stats['span_ok'][nm] += 1
                if n >= 3:
                    stats['span3'][nm] += 1
                    if span_fit(d['J'][nm0], bases) is not None:
                        stats['span3_ok'][nm] += 1


def new_stats():
    return dict(tri_bad=0, tri_entries=0, pairs2=Counter(), pairs3=Counter(), i=Counter(), ii=Counter(), iii=Counter(), i3=Counter(), ii3=Counter(), iii3=Counter(),
                span=Counter(), span_ok=Counter(), span3=Counter(), span3_ok=Counter())


def stage_within(products=None, verbose=True):
    hdr('B0b + B1: within-path blocks of the gl(4) triple vs the gl(3) reduced vertex operator V^(3)_red(kappa; lam123) on V_kappa[nu123]')
    products = products or (P_SMALL + P_BIG)
    tot = new_stats()
    examples = {}
    tab('%-26s %6s | %-13s | (i) a V+b: pairs>=2 / >=3 | (ii) +cD | (iii) products' % ('product', 'tri', 'op'))
    for lam, lamp in products:
        t0 = time.perf_counter()
        stats = new_stats()
        stage_within_product(lam, lamp, stats, examples)
        for k in stats:
            if isinstance(stats[k], Counter):
                tot[k].update(stats[k])
            else:
                tot[k] += stats[k]
        for nm in TRIPLE:
            k = (nm,)
            tab('%-26s tri_bad=%d | %-4s | (i) %d/%d  %d/%d | (ii) %d/%d  %d/%d | (iii) %d/%d  %d/%d   [%.1fs]' % (pstr(lam, lamp), stats['tri_bad'], nm, stats['i'][k], stats['pairs2'][k], stats['i3'][k],
                                                                                                                  stats['pairs3'][k], stats['ii'][k], stats['pairs2'][k], stats['ii3'][k], stats['pairs3'][k],
                                                                                                                  stats['iii'][k], stats['pairs2'][k], stats['iii3'][k], stats['pairs3'][k], time.perf_counter() - t0))
        for nm in ('V4a', 'V4a+D', 'V4b', 'V4b+D'):
            tab('%-26s span(J3,1[,D]) -> %-6s all pairs %d/%d, |P|>=3: %d/%d' % ('', nm, stats['span_ok'][nm], stats['span'][nm], stats['span3_ok'][nm], stats['span3'][nm]))
    tab('TOTAL: tridiagonality violations %d (over %d path-blocks x operators)' % (tot['tri_bad'], tot['tri_entries']))
    for nm in TRIPLE:
        k = (nm,)
        tab('TOTAL %-4s (i) %d/%d (|P|>=3: %d/%d)  (ii) %d/%d (%d/%d)  (iii) %d/%d (%d/%d)' % (nm, tot['i'][k], tot['pairs2'][k], tot['i3'][k], tot['pairs3'][k], tot['ii'][k], tot['pairs2'][k],
                                                                                              tot['ii3'][k], tot['pairs3'][k], tot['iii'][k], tot['pairs2'][k], tot['iii3'][k], tot['pairs3'][k]))
    for nm in ('V4a', 'V4a+D', 'V4b', 'V4b+D'):
        tab('TOTAL span %-6s all %d/%d, |P|>=3 %d/%d' % (nm, tot['span_ok'][nm], tot['span'][nm], tot['span3_ok'][nm], tot['span3'][nm]))
    return tot, examples


# ============================================================================================ B1b: the V4b within-path form
def within_V4b_fit(d):
    """J^{(4b)}_kappa = a V3 + e [xi, V3] + b + c xi + g xi^2 (xi = diag gl(2) Casimir x(x+2));  returns the coefficient list [a, e, b, c, g] or None"""
    n = d['n']
    V, D = d['V3'], d['D']
    comm = [[(D[i] - D[j]) * V[i][j] for j in range(n)] for i in range(n)]
    one = [[Fr(int(i == j)) for j in range(n)] for i in range(n)]
    Dm = [[D[i] if i == j else Fr(0) for j in range(n)] for i in range(n)]
    D2 = [[D[i] ** 2 if i == j else Fr(0) for j in range(n)] for i in range(n)]
    return span_fit(d['J']['V4b'], [V, comm, one, Dm, D2])


def stage_within_V4b(products=None):
    hdr('B1c: V4b within a path:  J = a V^(3) + e [xi, V^(3)] + b + c xi + g xi^2  (xi = x(x+2), x = m12 - m22); exact fits; g = -e/4 ?')
    products = products or (P_SMALL + P_BIG)
    tot = Counter()
    gcase = Counter()
    for lam, lamp in products:
        for mu, m, w in sectors(lam, lamp):
            if w < 2:
                continue
            for kappa, d in within_data(lam, lamp, mu).items():
                n = d['n']
                if n < 2:
                    continue
                tot['pairs'] += 1
                x = within_V4b_fit(d)
                if x is not None:
                    tot['fit'] += 1
                    if n >= 3:
                        tot['fit3'] += 1
                    if x[1] != 0 and n >= 3:
                        tot['e!=0'] += 1
                        gcase['g=-e/4'] += (x[4] == -x[1] / 4)
                if n >= 3:
                    tot['pairs3'] += 1
    tab('pairs (|P|>=2) %d: fit %d;  |P|>=3: %d/%d;  (|P|>=3) with e != 0 in %d, of which the fitted g = -e/4 in %d' % (tot['pairs'], tot['fit'], tot['fit3'], tot['pairs3'], tot['e!=0'], gcase['g=-e/4']))
    return tot


# ============================================================================================ held-out generic stage
def stage_heldout(seed=3, lams=((9, 6, 3, 0),), per=25):
    hdr("B1d: held-out: generic rational non-dominant lam-prime and random weight spaces of V_(9,6,3,0) (dim 4096); every (weight space, kappa) with |P_kappa| >= 2")
    import random
    tot = Counter()
    for lam in lams:
        R = gr.GLNRep.of(lam)
        for lamp in st.GENERIC_LAMPS:
            rng = random.Random(seed)
            wts = sorted(nu for nu in R.wsp if 6 <= len(R.wsp[nu]) <= 40)
            rng.shuffle(wts)
            for nu in wts[:per]:
                mu = tuple(nu[i] + lamp[i] for i in range(4))
                for kappa, d in within_data(lam, lamp, mu).items():
                    n = d['n']
                    if n < 2:
                        continue
                    tot[('pairs', n)] += 1
                    f3 = fit_ansatz(d['J']['V3'], d['V3'], d['D'])
                    f4 = fit_ansatz(d['J']['V4a'], d['V3'], d['D'])
                    tot[('V3 (i)', n)] += f3['ok1']
                    tot[('V4a (i)', n)] += f4['ok1']
                    tot[('V4b fit', n)] += within_V4b_fit(d) is not None
    for n in sorted({k[1] for k in tot}):
        tab('  |P|=%d: pairs %3d   V3 (i) %3d   V4a (i) %3d   V4b [a,e,b,c,xi^2] fit %3d' % (n, tot[('pairs', n)], tot[('V3 (i)', n)], tot[('V4a (i)', n)], tot[('V4b fit', n)]))
    tab('  total pairs %d, V3 (i) %d, V4a (i) %d, V4b fit %d' % (sum(v for k, v in tot.items() if k[0] == 'pairs'), sum(v for k, v in tot.items() if k[0] == 'V3 (i)'),
                                                             sum(v for k, v in tot.items() if k[0] == 'V4a (i)'), sum(v for k, v in tot.items() if k[0] == 'V4b fit')))
    return tot


# ============================================================================================ B2: M~ and the nested basis
def mpq(x):
    return mp.mpf(x.numerator) / x.denominator


def mtilde_space(lam, lamp, mu):
    """exact PRV-type object for M~ = ker E12^k1 cap ker E23^k2 in W (void strings skipped); same API as slice_engine.PRV"""
    lam, lamp, mu = tuple(lam), tuple(lamp), tuple(mu)
    R = gr.GLNRep.of(lam)
    N = 4
    w0 = sub_mu(mu, lamp)
    W = R.wsp[w0]
    rows = []
    for i in (1, 2):
        k = lamp[i - 1] - lamp[i] + 1
        if gr.wadd(w0, gr.dij(N, i, i + 1), k) not in R.wsp:
            continue
        imgs = se.power_images(R, (i, i + 1), k, W)
        byrow = defaultdict(lambda: [Fr(0)] * len(W))
        for a, vec in enumerate(imgs):
            for r, x in vec.items():
                byrow[r][a] = x
        rows.extend(byrow.values())
    basis, free = se._rref_nullspace(rows, len(W))
    return se.PRV(lam, lamp, mu, W, basis, free)


def kappa_prv(lam, lamp, mu, kappa):
    """the gl(3) PRV object of (kappa, lamp_123, mu_123) (an empty one if nu_123 is not a weight of V_kappa)"""
    return se.prv_space(kappa, tuple(lamp[:3]), tuple(mu[:3]))


def m3_of(kappa, lamp, mu):
    return gr.racah_speiser(kappa, tuple(lamp[:3])).get(tuple(mu[:3]), 0)


def embed3(R, R3, lam, vec3):
    """a vector of the gl(3) rep R3 {R3 index: x} -> the gl(4) global index {R index: x}  (pattern + (lam,))"""
    return {R.index[tuple(R3.pats[i]) + (tuple(lam),)]: x for i, x in vec3.items()}


def b2_exact(lam, lamp, mu):
    """exact B2(a),(b): M~ = (+)_kappa M^(3)_kappa embedded; dim = sum m3; M_mu in M~; the triple preserves M~.  Returns a dict of facts."""
    R, W, P, M = triple_on_W(lam, lamp, mu)
    pt = mtilde_space(lam, lamp, mu)
    pos = {int(c): a for a, c in enumerate(W)}
    out = dict(w=len(W), dimM=pt.m, kappas=len(P))
    dimsum = 0
    msum = 0
    ok_in = True
    for kappa in P:
        R3 = gr.GLNRep.of(kappa)
        pk = kappa_prv(lam, lamp, mu, kappa)
        dimsum += pk.m
        msum += m3_of(kappa, lamp, mu)
        for b in range(pk.m):
            v4 = embed3(R, R3, lam, pk.vec(b))
            vec = [Fr(0)] * len(W)
            for r, x in v4.items():
                vec[pos[r]] = x
            try:
                pt.coords(vec)
            except ValueError:
                ok_in = False
    out['dim_sum'] = dimsum
    out['rs_sum'] = msum
    out['emb_in_Mt'] = ok_in
    out['direct_sum'] = ok_in and dimsum == pt.m
    pres = True
    mats = {}
    for nm in TRIPLE:
        try:
            mats[nm] = so.vred_matrix(R, so.reduced_op(nm, lamp), pt)
        except (ValueError, KeyError):
            pres = False
            break
    out['preserved'] = pres
    out['mats'] = mats if pres else None
    pm = se.prv_space(lam, lamp, mu)
    out['m'] = pm.m
    inside = True
    for k in range(pm.m):
        try:
            pt.coords(list(pm.basis[k]))
        except ValueError:
            inside = False
    out['Mmu_in_Mt'] = inside
    out['pt'] = pt
    return out


def eigen_nested(lam, lamp, mu, dps=50):
    """the nested basis: for each kappa with m3 > 0 the gl(3) canonical copies (eigenvectors of V^(3)_red on M^(3)_kappa, ascending eigenvalue, 50 digits, first nonzero component of the
    W-vector = 1).  Returns list of dict(kappa, s, ev3, vec = {global index: mpf})"""
    R, W, P, M = triple_on_W(lam, lamp, mu)
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        out = []
        lamp3 = tuple(lamp[:3])
        tolz = mp.mpf(10) ** (-(dps - 20))
        for kappa in sorted(P):
            pk = kappa_prv(lam, lamp, mu, kappa)
            if pk.m == 0:
                continue
            R3 = gr.GLNRep.of(kappa)
            A = so.vred_matrix(R3, so.reduced_op('V', lamp3), pk)
            jd = so.joint_diag([A], dps=dps, combo=(Fr(1),))
            order = sorted(range(pk.m), key=lambda k: jd['trip'][k][0])
            Bm = [[mpq(x) for x in row] for row in pk.basis]
            for s, k in enumerate(order):
                col = jd['U'][k]
                sl = [mp.fsum(col[kk] * Bm[kk][a] for kk in range(pk.m)) for a in range(pk.w)]
                vec3 = {int(pk.W[a]): sl[a] for a in range(pk.w)}
                vec4 = embed3(R, R3, lam, vec3)
                mx0 = max(abs(x) for x in vec4.values())
                first = min(i for i, x in vec4.items() if abs(x) > tolz * mx0)
                sc = vec4[first]
                vec4 = {i: x / sc for i, x in vec4.items()}
                mx = max(abs(x) for x in vec4.values())
                vec4 = {i: x for i, x in vec4.items() if abs(x) > tolz * mx}
                out.append(dict(kappa=kappa, s=s + 1, ev3=jd['trip'][k][0], vec=vec4, offmax=jd['offmax'], ok=jd['ok']))
        return out
    finally:
        mp.mp.dps = old


def mp_matrix_from_vecs(vecs, W):
    pos = {int(c): a for a, c in enumerate(W)}
    B = mp.matrix(len(W), len(vecs))
    for j, v in enumerate(vecs):
        for r, x in v.items():
            B[pos[r], j] = x
    return B


def exact_to_mp(M):
    n = len(M)
    A = mp.matrix(n, len(M[0]) if n else 0)
    for i in range(n):
        for j in range(len(M[0])):
            A[i, j] = mpq(M[i][j])
    return A


def nested_triple(lam, lamp, mu, dps=50, nb=None):
    """the gl(4) triple in the nested basis: dict(nb, B = (w x mt) mp matrix, T = {name: (mt x mt) mp matrix}, res = max residual |M B - B T|, labels = [(kappa, s)])"""
    R, W, P, M = triple_on_W(lam, lamp, mu)
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        nb = nb or eigen_nested(lam, lamp, mu, dps)
        B = mp_matrix_from_vecs([x['vec'] for x in nb], W)
        T = {}
        res = mp.mpf(0)
        Bp = mp.inverse(B.T * B) * B.T
        for nm in TRIPLE:
            MB = exact_to_mp(M[nm]) * B
            Tm = Bp * MB
            res = max(res, max(abs(x) for x in (B * Tm - MB)))
            T[nm] = Tm
        return dict(nb=nb, B=B, Bp=Bp, T=T, res=res, labels=[(x['kappa'], x['s']) for x in nb])
    finally:
        mp.mp.dps = old


# ============================================================================================ B2 analysis
def mp_rank(A, tol):
    """numerical rank of the mp matrix A (singular values > tol * max(1, largest))"""
    if A.rows == 0 or A.cols == 0:
        return 0
    s = mp.svd_r(A, compute_uv=False)
    sv = [abs(s[i]) for i in range(len(s))]
    top = max(sv)
    if top == 0:
        return 0
    return sum(1 for x in sv if x > tol * max(top, 1))


def apply_mp(R, ij, k, vec):
    """E_ij^k on a vector {index: mpf} (the Fractions are converted)"""
    for _ in range(k):
        out = {}
        cols = R.E[ij]
        for c, v in vec.items():
            for r, e in cols[c].items():
                out[r] = out.get(r, 0) + v * mpq(e)
        vec = out
    return vec


def vnorm(vec):
    return max([abs(x) for x in vec.values()] + [mp.mpf(0)])


def nested_analysis(lam, lamp, mu, dps=50, tol_exp=30, detail=False):
    """everything B2 asks for one sector; returns a dict"""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        tol = mp.mpf(10) ** (-tol_exp)
        R, W, P, M = triple_on_W(lam, lamp, mu)
        pos = {int(c): a for a, c in enumerate(W)}
        b2 = b2_exact(lam, lamp, mu)
        out = dict(mu=mu, w=len(W), m=b2['m'], mt=b2['dimM'], kappas=len(P), b2={k: v for k, v in b2.items() if k not in ('mats', 'pt')})
        nt = nested_triple(lam, lamp, mu, dps)
        labels = nt['labels']
        out['labels'] = labels
        out['res'] = nt['res']
        mt = len(labels)
        # blocks of T in the nested basis
        kidx = defaultdict(list)
        for a, (k, s) in enumerate(labels):
            kidx[k].append(a)
        diag_off = {}
        for nm in TRIPLE:
            T = nt['T'][nm]
            worst = mp.mpf(0)
            for k, ids in kidx.items():
                for a in ids:
                    for b in ids:
                        if a != b:
                            worst = max(worst, abs(T[a, b]))
            diag_off[nm] = worst
        out['diag_off'] = diag_off
        # hopping blocks
        hop = []
        for nm in TRIPLE:
            T = nt['T'][nm]
            for k2, i2 in kidx.items():
                for k1, i1 in kidx.items():
                    if k1 == k2:
                        continue
                    A = mp.matrix(len(i2), len(i1))
                    for a, r in enumerate(i2):
                        for b, c in enumerate(i1):
                            A[a, b] = T[r, c]
                    nz = sum(1 for a in range(len(i2)) for b in range(len(i1)) if abs(A[a, b]) > tol)
                    if nz == 0:
                        continue
                    rk = mp_rank(A, tol)
                    hop.append(dict(op=nm, k2=k2, k1=k1, shape=(len(i2), len(i1)), nz=nz, rank=rk, A=A, i2=i2, i1=i1))
        out['hop'] = hop
        # lattice hopping blocks (exact) and their support relative to the path ends
        lat = []
        for nm in ('V3',):
            for k2, p2 in P.items():
                for k1, p1 in P.items():
                    if k1 == k2:
                        continue
                    H = block(M[nm], pos, p2, p1)
                    sup = [(a, b) for a in range(len(p2)) for b in range(len(p1)) if H[a][b] != 0]
                    if sup:
                        ends = all(a in (0, len(p2) - 1) and b in (0, len(p1) - 1) for a, b in sup)
                        lat.append(dict(k2=k2, k1=k1, n2=len(p2), n1=len(p1), nsup=len(sup), ends=ends, dm=sorted(set(R.pats[p2[a]][1][0] - R.pats[p1[b]][1][0] for a, b in sup))))
        out['lat'] = lat
        # physical copies in the nested basis
        cps = so.copies(lam, lamp, mu, dps=dps)
        Cm = []
        resc = mp.mpf(0)
        for cp in cps:
            sv = mp.matrix(len(W), 1)
            for r, x in cp['slice'].items():
                sv[pos[r]] = x
            c = nt['Bp'] * sv
            resc = max(resc, max(abs(x) for x in (nt['B'] * c - sv)))
            Cm.append([c[i] for i in range(mt)])
        out['cps'] = cps
        out['C'] = Cm
        out['res_copies'] = resc
        # joint spectrum of the triple on M~: contains Spec(M_mu); the others violate E34^k3
        mats = b2['mats']
        pt = b2['pt']
        if mats is not None and pt.m > 0:
            jd = so.joint_diag([mats[n] for n in TRIPLE], dps=dps)
            Bm = [[mpq(x) for x in row] for row in pt.basis]
            k3 = lamp[2] - lamp[3] + 1
            nu = sub_mu(mu, lamp)
            void3 = gr.wadd(nu, gr.dij(4, 3, 4), k3) not in R.wsp
            phys = 0
            viol = []
            for i, trip in enumerate(jd['trip']):
                match = any(max(abs(trip[t] - cp['ev'][t]) for t in range(3)) < tol for cp in cps)
                phys += match
                if not match:
                    col = jd['U'][i]
                    sl = [mp.fsum(col[kk] * Bm[kk][a] for kk in range(pt.m)) for a in range(pt.w)]
                    vec = {int(pt.W[a]): sl[a] for a in range(pt.w)}
                    if void3:
                        viol.append(mp.mpf(0))
                    else:
                        img = apply_mp(R, (3, 4), k3, vec)
                        viol.append(vnorm(img) / vnorm(vec))
            out['spec'] = dict(n_eig=len(jd['trip']), n_phys_matched=phys, n_other=len(viol), min_viol=min(viol) if viol else None, void3=void3, ok=jd['ok'], mmatch=(phys == len(cps)))
        else:
            out['spec'] = dict(n_eig=0, n_phys_matched=0, n_other=0, min_viol=None, void3=None, ok=True, mmatch=(len(cps) == 0))
        out['nt'] = nt
        out['P'] = P
        out['R'] = R
        return out
    finally:
        mp.mp.dps = old


def _prop(vec, ref, tol, pw=1):
    """is vec proportional to ref**pw (componentwise, zeros matched)?"""
    ratio = None
    for x, y in zip(vec, ref):
        if abs(y) < tol:
            if pw == -1:
                return False
            if abs(x) > tol:
                return False
            continue
        r = x / (y ** pw)
        if ratio is None:
            ratio = r
        elif abs(r - ratio) > tol * max(1, abs(ratio)):
            return False
    return ratio is not None and abs(ratio) > tol


def hop_rank_one_test(an, tol_exp=25):
    """for each rank-one hopping block A = f g^T (non-trivial only for the factor of length >= 2): is g proportional to the first / last component of the gl(3) copies of k1 (the end of the
    chain in the lattice basis) and f to (the end component of the copies of k2)^(+-1)?  counts: 'g tests' / 'f tests' (denominators) and successes per end"""
    tol = mp.mpf(10) ** (-tol_exp)
    nb = an['nt']['nb']
    P = an['P']
    res = Counter()
    for h in an['hop']:
        if h['rank'] != 1:
            continue
        A = h['A']
        imax, jmax = max(((a, b) for a in range(A.rows) for b in range(A.cols)), key=lambda ab: abs(A[ab[0], ab[1]]))
        f = [A[a, jmax] for a in range(A.rows)]
        g = [A[imax, b] for b in range(A.cols)]
        for which, vec, k, ids in (('g', g, h['k1'], h['i1']), ('f', f, h['k2'], h['i2'])):
            if len(vec) < 2:
                continue
            res[which + ' tests'] += 1
            path = P[k]
            for endn, c in (('first', path[0]), ('last', path[-1])):
                ref = [nb[i]['vec'].get(c, mp.mpf(0)) for i in ids]
                for pw in (1, -1):
                    if _prop(vec, ref, tol, pw):
                        res['%s ~ %s end^%+d' % (which, endn, pw)] += 1
    return res


def fm(x, d=3):
    return mp.nstr(x, d)


def stage_nested(products=None, detail_products=(), maxdetail=40):
    hdr('B2: M~ = (+)_kappa M^(3)_kappa, the nested basis |kappa,s>, the triple in it (50 digits; thresholds 1e-30)')
    products = products or [(((3, 2, 1, 0),) * 2)] + [((4, 2, 1, 0),) * 2] + [((6, 4, 2, 0),) * 2]
    glob = Counter()
    for lam, lamp in products:
        t0 = time.perf_counter()
        st_ = Counter()
        worst = {nm: mp.mpf(0) for nm in TRIPLE}
        offcount = {nm: [0, 0] for nm in TRIPLE}
        rkhist = Counter()
        fullhist = Counter()
        fracs = []
        latc = Counter()
        r1 = Counter()
        crk = Counter()
        maxres = mp.mpf(0)
        minviol = None
        lines = []
        for mu, m, w in sectors(lam, lamp):
            if m == 0:
                st_['formal sectors (m=0) skipped'] += 1
                continue
            an = nested_analysis(lam, lamp, mu)
            b2 = an['b2']
            st_['sectors'] += 1
            st_['direct sum & dim = RS'] += b2['direct_sum'] and b2['dim_sum'] == b2['rs_sum']
            st_['triple preserves M~'] += b2['preserved']
            st_['M_mu in M~'] += b2['Mmu_in_Mt']
            st_['spec(M_mu) in spec(M~)'] += an['spec']['mmatch']
            if an['spec']['n_other'] > 0:
                st_['sectors with spurious eigenvectors on M~'] += 1
                if not an['spec']['void3']:
                    mv = an['spec']['min_viol']
                    minviol = mv if minviol is None else min(minviol, mv)
                    st_['spurious violate E34^k3 (min ratio > 1e-10)'] += mv > mp.mpf(10) ** (-10)
                else:
                    st_['spurious with void E34 string (!)'] += 1
            maxres = max(maxres, an['res'], an['res_copies'])
            for nm in TRIPLE:
                worst[nm] = max(worst[nm], an['diag_off'][nm])
            kcount = Counter(k for k, s in an['labels'])
            nbig = sum(1 for k, c in kcount.items() if c >= 2)
            for nm in TRIPLE:
                offcount[nm][1] += nbig
                if nbig:
                    # per-kappa diagonality
                    T = an['nt']['T'][nm]
                    kidx = defaultdict(list)
                    for a, (k, s) in enumerate(an['labels']):
                        kidx[k].append(a)
                    for k, ids in kidx.items():
                        if len(ids) >= 2:
                            ow = max(abs(T[a, b]) for a in ids for b in ids if a != b)
                            offcount[nm][0] += ow < mp.mpf(10) ** (-30)
            for h in an['hop']:
                mn = min(h['shape'])
                rkhist[(h['op'], h['rank'], mn)] += 1
                if mn >= 2:
                    fullhist[(h['op'], h['rank'] == mn)] += 1
                fracs.append(h['nz'] / (h['shape'][0] * h['shape'][1]))
            for l in an['lat']:
                latc['pairs'] += 1
                latc['support only at path ends'] += l['ends']
            rr = hop_rank_one_test(an)
            r1.update(rr)
            # physical copies: c(kappa, s) arrays
            kap = sorted(set(k for k, s in an['labels']))
            smax = max(s for k, s in an['labels'])
            ncz = 0
            for c in an['C']:
                A = mp.matrix(len(kap), smax)
                for a, (k, s) in enumerate(an['labels']):
                    A[kap.index(k), s - 1] = c[a]
                    ncz += abs(c[a]) < mp.mpf(10) ** (-30)
                crk[(mp_rank(A, mp.mpf(10) ** (-30)), min(len(kap), smax))] += 1
            Cm = mp.matrix(len(an['C']), len(an['labels']))
            for i, c in enumerate(an['C']):
                for j in range(len(c)):
                    Cm[i, j] = c[j]
            st_['C (m x mt) full row rank m'] += mp_rank(Cm, mp.mpf(10) ** (-30)) == len(an['C'])
            st_['zero entries c(kappa,s) (below 1e-30)'] += ncz
            st_['entries c(kappa,s)'] += len(an['C']) * len(an['labels'])
            n4 = depth_of(mu, lamp, lam)
            lines.append('   mu=%-14s m=%2d w=%2d mt=%2d #kappa=%d  d=%d  offdiag(V3,V4a,V4b)=%s %s %s  res=%s' % (lstr(mu), m, w, len(an['labels']), len(kap), n4, fm(an['diag_off']['V3'], 2),
                                                                                                          fm(an['diag_off']['V4a'], 2), fm(an['diag_off']['V4b'], 2), fm(max(an['res'], an['res_copies']), 2)))
        tab('--- %s  (%.1fs)' % (pstr(lam, lamp), time.perf_counter() - t0))
        for k, v in st_.items():
            tab('   %-52s %d' % (k, v))
        tab('   largest nested-basis residual (|M B - B T|, |B c - slice|)          %s' % fm(maxres, 2))
        tab('   largest off-diagonal entry inside a kappa-block of T: V3 %s  V4a %s  V4b %s' % tuple(fm(worst[nm], 2) for nm in TRIPLE))
        tab('   kappa-blocks (|nested| >= 2) diagonal to 1e-30: V3 %d/%d  V4a %d/%d  V4b %d/%d' % tuple(x for nm in TRIPLE for x in offcount[nm]))
        if minviol is not None:
            tab('   smallest |E34^k3 v|/|v| over the spurious eigenvectors of M~ (non-void E34): %s' % fm(minviol, 3))
        tab('   hopping blocks T[k2,k1] (non-zero ones): (op, rank, min dim) histogram: ' + ', '.join('%s r%d/%d: %d' % (a, b, c, n) for (a, b, c), n in sorted(rkhist.items())))
        tab('   hopping blocks with min dim >= 2: full rank / total: ' + ', '.join('%s %d/%d' % (op, fullhist[(op, True)], fullhist[(op, True)] + fullhist[(op, False)]) for op in TRIPLE))
        if fracs:
            tab('   fraction of entries above 1e-30 in a non-zero hopping block: mean %.3f, min %.3f' % (sum(fracs) / len(fracs), min(fracs)))
        tab('   lattice-basis hopping blocks H[k2,k1] of V3 (non-zero): %d, supported only on the (first/last, first/last) corner entries: %d' % (latc['pairs'], latc['support only at path ends']))
        tab('   rank-one hopping blocks, factor tests (T = f g^T, only factors of length >= 2 counted): %s' % dict(sorted(r1.items())))
        tab('   physical copies: rank of the (#kappa x s_max) array c(kappa,s) per copy vs its max possible: ' + ', '.join('r%d/%d: %d' % (a, b, n) for (a, b), n in sorted(crk.items())))
        if (lam, lamp) in detail_products:
            for ln in lines[:maxdetail]:
                tab(ln)
        glob.update(st_)
    return glob


# ============================================================================================ B3: depth, the universal hopping formula, nested solution
def depth_of(mu, lamp, lam=None):
    """depth d = nu4 - lam4 = |lam_123| - |kappa| (the brief's n = nu4 coincides with d when lam4 = 0; for lam4 > 0 the true depth is nu4 - lam4: kappa = lam_123 is d = 0)"""
    return sub_mu(mu, lamp)[3] - (lam[3] if lam is not None else 0)


FCL = set(tuple(d) for d in st.f_classes())


def lattice_hop_check(lam, lamp, mu):
    """exact: every lattice matrix element of the triple between two DIFFERENT paths with |dm12| = 1 (the 12 universal classes) equals structure_gl4.universal_F; also counts of the other entries
    (dm12 = 0 hopping: the cubic Q classes, not re-tested here).  Returns Counter(tests, nonzero, bad, other_hops, within_path_entries)"""
    R, W, P, M = triple_on_W(lam, lamp, mu)
    pos = {int(c): a for a, c in enumerate(W)}
    cnt = Counter()
    for c in W:
        src = R.pats[int(c)]
        pt = st.pt_of_pattern(src, tuple(lam), tuple(lamp))
        ent = st.entries(src)
        for r in W:
            tgt = R.pats[int(r)]
            if tgt[2] == src[2]:
                cnt['within_path_entries'] += (M['V3'][pos[int(r)]][pos[int(c)]] != 0)
                continue
            d = tuple(b - a for a, b in zip(ent, st.entries(tgt)))
            if d not in FCL:
                cnt['other_hops (dm12=0 or outside the 12 classes)'] += (M['V3'][pos[int(r)]][pos[int(c)]] != 0)
                continue
            for nm in TRIPLE:
                try:
                    val = st.universal_F(nm, pt, d)
                except ZeroDivisionError:
                    cnt['pole'] += 1
                    continue
                exp_ = M[nm][pos[int(r)]][pos[int(c)]]
                cnt['tests'] += 1
                cnt['nonzero'] += exp_ != 0
                cnt['bad'] += val != exp_
    return cnt


def nested_copies(lam, lamp, mu, dps=50, an=None):
    """solve the triple in the NESTED coordinates: joint eigenvectors c of (T_V3, T_V4a, T_V4b) on M~ (size mt = sum_kappa m3(kappa)), keep those with E34^k3 (B c) = 0.  Returns a list of
    dict(ev, c, slice) sorted by ev; slice normalised like so.copies is NOT attempted (projective comparison in validate_nested_copies)"""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        R, W, P, M = triple_on_W(lam, lamp, mu)
        nt = an['nt'] if an is not None else nested_triple(lam, lamp, mu, dps)
        T = nt['T']
        mt = len(nt['labels'])
        tol = mp.mpf(10) ** (-30)
        combo = so.GENERIC
        Mc = mp.matrix(mt, mt)
        for cf, nm in zip(combo, TRIPLE):
            Mc += mpq(Fr(cf)) * T[nm]
        E, ER = mp.eig(Mc)
        k3 = lamp[2] - lamp[3] + 1
        nu = sub_mu(mu, lamp)
        void3 = gr.wadd(nu, gr.dij(4, 3, 4), k3) not in R.wsp
        out = []
        for i in range(mt):
            col = [ER[a, i] for a in range(mt)]
            big = max(range(mt), key=lambda a: abs(col[a]))
            col = [mp.re(x / col[big]) for x in col]
            cv = mp.matrix(col)
            ev = tuple((T[nm] * cv)[big] for nm in TRIPLE)
            wv = nt['B'] * cv
            vec = {int(W[a]): wv[a] for a in range(len(W)) if abs(wv[a]) > mp.mpf(10) ** (-45)}
            if void3:
                phys = True
            else:
                img = apply_mp(R, (3, 4), k3, vec)
                phys = vnorm(img) < tol * max(1, vnorm(vec))
            if phys:
                out.append(dict(ev=ev, c=col, slice=vec))
        out.sort(key=lambda x: tuple(x['ev']))
        return out
    finally:
        mp.mp.dps = old


def validate_nested_copies(lam, lamp, mu, dps=50):
    """compare nested_copies with so.copies: eigenvalues (max deviation after matching) and slices (projective residual)"""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        cps = so.copies(lam, lamp, mu, dps=dps)
        nc = nested_copies(lam, lamp, mu, dps)
        if len(cps) != len(nc):
            return dict(ok=False, n=len(cps), nn=len(nc), dev=None, sres=None)
        dev = mp.mpf(0)
        sres = mp.mpf(0)
        for cp in cps:
            best = min(nc, key=lambda x: max(abs(x['ev'][t] - cp['ev'][t]) for t in range(3)))
            dev = max(dev, max(abs(best['ev'][t] - cp['ev'][t]) for t in range(3)))
            j = max(cp['slice'], key=lambda r: abs(cp['slice'][r]))
            sc = best['slice'].get(j, mp.mpf(0)) / cp['slice'][j]
            if sc == 0:
                sres = mp.inf
                continue
            r = max(abs(best['slice'].get(i, 0) / sc - cp['slice'].get(i, 0)) for i in set(best['slice']) | set(cp['slice']))
            sres = max(sres, r)
        return dict(ok=(dev < mp.mpf(10) ** (-30) and sres < mp.mpf(10) ** (-30)), n=len(cps), nn=len(nc), dev=dev, sres=sres)
    finally:
        mp.mp.dps = old


def depth1_copies(lam, lamp, mu, dps=50):
    """the depth-1 sector (nu4 = mu4 - lamp4 = 1, at most 3 paths): the physical copies from three gl(3) problems (the nested basis: gl(3) copies of each kappa = lam_123 - e_i) coupled by the
    hopping blocks, solved as ONE mt x mt joint eigenproblem in the nested basis (mt = sum m3(kappa) <= 3 m_max); validated against slice_ops_gl4.copies by validate_nested_copies"""
    assert depth_of(mu, lamp, lam) == 1, 'not a depth-1 sector'
    return nested_copies(lam, lamp, mu, dps)


def stage_depth(products=None, detail=((((3, 2, 1, 0),) * 2), (((4, 2, 1, 0),) * 2)), validate=True):
    hdr('B3: depth d = nu4 - lam4 = |lam_123| - |kappa| (= nu4 = mu4 - lamp4 when lam4 = 0): paths, dim M~, n = 0 reduces to gl(3), the depth-1 sectors')
    products = products or (P_SMALL + P_BIG)
    tabn = {}
    n0 = Counter()
    n0par = []
    hopc = Counter()
    valc = Counter()
    for lam, lamp in products:
        t0 = time.perf_counter()
        byn = defaultdict(list)
        for mu, m, w in sectors(lam, lamp):
            if m == 0:
                continue
            R, W, P, M = triple_on_W(lam, lamp, mu)
            n = depth_of(mu, lamp, lam)
            mt = mtilde_space(lam, lamp, mu).m
            byn[n].append((mu, len(P), w, m, mt))
            if n == 0:
                n0['sectors'] += 1
                n0['one path'] += (len(P) == 1)
            if n >= 1:
                hc = lattice_hop_check(lam, lamp, mu)
                for k, v in hc.items():
                    hopc[k] += v
            if n == 1 and validate:
                v = validate_nested_copies(lam, lamp, mu)
                valc['depth-1 sectors'] += 1
                valc['nested solution = so.copies (1e-30)'] += v['ok']
                valc['copies'] += v['n']
        tab('--- %s  [%.1fs]' % (pstr(lam, lamp), time.perf_counter() - t0))
        tab('    d | sectors | #kappa (hist) | max w | max m | max mt | mt = m')
        for n in sorted(byn):
            L = byn[n]
            hist = Counter(x[1] for x in L)
            tab('    %d | %3d | %s | %2d | %2d | %2d | %d/%d' % (n, len(L), ' '.join('%d:%d' % kv for kv in sorted(hist.items())), max(x[2] for x in L), max(x[3] for x in L), max(x[4] for x in L),
                                                                  sum(1 for x in L if x[3] == x[4]), len(L)))
        tabn[(lam, lamp)] = byn
    tab('depth-0 sectors: %d, with exactly one path: %d' % (n0['sectors'], n0['one path']))
    tab('lattice hopping between different paths, |dm12| = 1 entries vs structure_gl4.universal_F (all sectors with n >= 1): %s' % dict(hopc))
    tab('depth-1 sectors solved in the nested basis: %s' % dict(valc))
    return tabn


def void3_of(lam, lamp, mu):
    R = gr.GLNRep.of(lam)
    nu = sub_mu(mu, lamp)
    k3 = lamp[2] - lamp[3] + 1
    return gr.wadd(nu, gr.dij(4, 3, 4), k3) not in R.wsp


def stage_depth_detail(prods=(((3, 2, 1, 0),) * 2, ((4, 2, 1, 0),) * 2, ((5, 3, 2, 1),) * 2), maxn=1, show_c=True):
    hdr('B3 detail: sectors of depth <= %d: paths, within-path parameters (J = a V^(3) + b), hopping, dim M~' % maxn)
    tot = Counter()
    for lam, lamp in prods:
        tab('--- %s' % pstr(lam, lamp))
        for mu, m, w in sectors(lam, lamp):
            if m == 0:
                continue
            n = depth_of(mu, lamp, lam)
            if n > maxn:
                continue
            R, W, P, M = triple_on_W(lam, lamp, mu)
            mt = mtilde_space(lam, lamp, mu).m
            data = within_data(lam, lamp, mu)
            parts = []
            for kappa in sorted(P):
                d = data[kappa]
                nn = d['n']
                m3 = m3_of(kappa, lamp, mu)
                if nn >= 2:
                    f3 = fit_ansatz(d['J']['V3'], d['V3'], d['D'])
                    f4 = fit_ansatz(d['J']['V4a'], d['V3'], d['D'])
                    fb = within_V4b_fit(d)
                    s = 'k=%s |P|=%d m3=%d  V3: a=%s b=%s | V4a: a=%s b=%s | V4b: %s' % (lstr(kappa), nn, m3, f3['a'], f3['b'], f4['a'], f4['b'], ('[a,e,b,c,g]=' + str([str(x) for x in fb])) if fb else 'NO FIT')
                else:
                    dv = [d['J'][nm][0][0] - d['V3'][0][0] for nm in TRIPLE]
                    s = 'k=%s |P|=1 m3=%d  J-V^(3)(1x1): V3 %s, V4a %s, V4b %s  (J = %s, %s, %s)' % (lstr(kappa), m3, dv[0], dv[1], dv[2], d['J']['V3'][0][0], d['J']['V4a'][0][0], d['J']['V4b'][0][0])
                parts.append(s)
            tab('  mu=%s  d=%d  w=%d m=%d mt=%d  E34 void: %s' % (lstr(mu), n, w, m, mt, void3_of(lam, lamp, mu)))
            for s in parts:
                tab('       ' + s)
            tot['sectors'] += 1
    return tot


def stage_mt_equals_m(products=None):
    hdr('M~ = M_mu exactly when the E34 string is void (n < lamp3 - lamp4 + 1), otherwise dim M~ - m extra spurious eigenvectors')
    products = products or (P_SMALL + P_BIG)
    c = Counter()
    for lam, lamp in products:
        for mu, m, w in sectors(lam, lamp):
            if m == 0:
                continue
            mt = mtilde_space(lam, lamp, mu).m
            v = void3_of(lam, lamp, mu)
            c[('void', v, 'mt=m' if mt == m else 'mt>m')] += 1
            # M_mu = M~ cap ker E34^k3 : dimension check is done in nested_analysis; here only counts
    for k, v in sorted(c.items(), key=str):
        tab('  E34 void=%s : %s : %d sectors' % (k[1], k[2], v))
    return c


# ============================================================================================ examples, c-tables, CLI
def stage_example():
    hdr('B1 failing example (smallest): V4b vs the gl(3) operator; the two tridiagonal matrices in full')
    lam = lamp = (2, 1, 1, 0)
    for mu, m, w in sectors(lam, lamp):
        if w < 2:
            continue
        for kappa, d in within_data(lam, lamp, mu).items():
            if d['n'] == 2:
                f = fit_ansatz(d['J']['V4b'], d['V3'], d['D'])
                if not f['ok1']:
                    R, W, P, M = triple_on_W(lam, lamp, mu)
                    tab('  product %s, mu = %s, kappa = %s, path (m12) = %s, xi = %s' % (pstr(lam, lamp), lstr(mu), lstr(kappa), [R.pats[c][1][0] for c in d['path']], [str(x) for x in d['D']]))
                    tab('  J^(4b)  = %s' % [[str(x) for x in r] for r in d['J']['V4b']])
                    tab('  V^(3)   = %s' % [[str(x) for x in r] for r in d['V3']])
                    tab('  J^(3)   = %s' % [[str(x) for x in r] for r in d['J']['V3']])
                    tab('  ratios of corresponding off-diagonal entries (J4b / V3): lower %s, upper %s  (a V + b would need them equal)' % (d['J']['V4b'][1][0] / d['V3'][1][0], d['J']['V4b'][0][1] / d['V3'][0][1]))
                    x = within_V4b_fit(d)
                    tab('  fit J4b = a V + e [xi,V] + b + c xi + g xi^2: [a,e,b,c,g] = %s (n = 2: g, c under-determined, a minimal solution is printed)' % ([str(t) for t in x] if x else None))
                    return


def stage_ctables(items=None, maxm=40):
    hdr('B2: physical copies in the nested basis, c(kappa, s)  (50 digits computed, 10 shown); rows = copies sorted by (v3, v4a, v4b); columns = nested vectors (kappa, s)')
    if items is None:
        items = [((3, 2, 1, 0),) * 2 + (None,), ((4, 2, 1, 0),) * 2 + (None,), ((6, 4, 2, 0),) * 2 + (6,)]
    for lam, lamp, mmin in items:
        for mu, m, w in sectors(lam, lamp):
            if m == 0 or (mmin is not None and m < mmin) or m < 2:
                continue
            an = nested_analysis(lam, lamp, mu)
            tab('--- %s  mu = %s  d=%d  m=%d  mt=%d  paths=%d' % (pstr(lam, lamp), lstr(mu), depth_of(mu, lamp, lam), m, len(an['labels']), an['kappas']))
            tab('    columns: ' + '  '.join('%s#%d' % (lstr(k), s) for k, s in an['labels']))
            for cp, c in zip(an['cps'], an['C']):
                tab('    v=(%s,%s,%s): ' % tuple(mp.nstr(x, 8) for x in cp['ev']) + '  '.join(mp.nstr(x, 10) for x in c))


def stage_depth0(products=None):
    hdr('B3: depth 0 (one path, kappa = lam_123): the gl(4) triple on the path in terms of the gl(3) operator V^(3)_red(lam_123; lamp_123)')
    products = products or (P_SMALL + P_BIG)
    c = Counter()
    ex = {}
    for lam, lamp in products:
        for mu, m, w in sectors(lam, lamp):
            if depth_of(mu, lamp, lam) != 0:
                continue
            d = within_data(lam, lamp, mu)
            assert len(d) == 1
            (kappa, dd), = d.items()
            n = dd['n']
            if lam[3] == 0:
                c['l4=0 sectors'] += 1
                c['l4=0: J3 = V^(3), J4a = 0, J4b = 0 exactly'] += (dd['J']['V3'] == dd['V3'] and all(x == 0 for r in dd['J']['V4a'] for x in r) and all(x == 0 for r in dd['J']['V4b'] for x in r))
            else:
                c['l4>0 sectors'] += 1
                if n >= 2:
                    c['l4>0, |P|>=2'] += 1
                    f3 = fit_ansatz(dd['J']['V3'], dd['V3'], dd['D'])
                    f4 = fit_ansatz(dd['J']['V4a'], dd['V3'], dd['D'])
                    fb = fit_ansatz(dd['J']['V4b'], dd['V3'], dd['D'])
                    c['l4>0: J3 = V + b3, J4a = V + b4a, J4b = -V + b4b (a = 1, 1, -1)'] += (f3['ok1'] and f3['a'] == 1 and f4['ok1'] and f4['a'] == 1 and fb['ok1'] and fb['a'] == -1)
                    ex[(lam, mu)] = (f3['b'], f4['b'], fb['b'], lam[3], lamp[3], sum(lam), sum(lamp))
    for k, v in c.items():
        tab('  %-80s %d' % (k, v))
    for k, v in sorted(ex.items())[:6]:
        tab('  example %s mu=%s: b3, b4a, b4b = %s, %s, %s' % (lstr(k[0]), lstr(k[1]), v[0], v[1], v[2]))
    return c


def main(argv):
    stages = ['ident', 'within', 'nested', 'depth']
    if len(argv) < 2:
        print(__doc__)
        return
    if argv[1] == 'all':
        run = stages
    elif argv[1] == 'run':
        run = argv[2].split(',')
    else:
        run = [argv[1]]
    t0 = time.perf_counter()
    mp.mp.dps = 50
    for s in run:
        t1 = time.perf_counter()
        if s == 'ident':
            stage_ident()
        elif s == 'within':
            stage_within()
            stage_within_V4b()
            stage_heldout()
            stage_example()
        elif s == 'nested':
            stage_nested(products=[((3, 2, 1, 0),) * 2, ((4, 2, 1, 0),) * 2, ((3, 2, 1, 0), (4, 2, 1, 0)), ((6, 4, 2, 0),) * 2, ((5, 3, 1, 0),) * 2, ((5, 3, 2, 1),) * 2])
            stage_ctables()
        elif s == 'depth':
            stage_mt_equals_m()
            stage_depth()
            stage_depth_detail()
            stage_depth0()
        else:
            print('unknown stage', s)
        tab('[stage %s: %.1fs]' % (s, time.perf_counter() - t1))
    tab('[total %.1fs]' % (time.perf_counter() - t0))
    flush_tables(mode='w' if argv[1] != 'run' else 'a', title='nested_gl3_gl4 tables (%s)' % DATE)


if __name__ == '__main__':
    main(sys.argv)
