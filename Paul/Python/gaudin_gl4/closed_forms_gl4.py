"""closed_forms_gl4.py -- G8d: closed forms of the gl(3)-in-gl(4) nesting scalars and the gl(4) triple assembled from labels alone (2026-10-24).

Brief: SPEC-G8d-closed-forms.md.  Results: RESULTS-G8d-closed-forms.md.  Everything is exact (Fractions) unless "50 digits" is said.

Public API
  a_of(kappa, lam), a4a_of(kappa, lam)                  the nesting scalars (D1)
  b_of(name, kappa, lam, nu, lamp)                      the shifts b3, b4a (D2), name in {'V3', 'V4a'}
  gl3_chain(kappa, lamp3, nu123)                        the gl(3) Jacobi chain (u, d, alpha) for GENERAL lamp3 (closed form)
  word_groups(...)                                      D3: kappa-blocks of the word groups alpha / beta
  closed_form_mats(lam, lamp, mu)                       the w x w matrices of V3, V4a on W = V_lam[mu - lamp] from labels alone (no representation matrix); V4b not included
  slices_closed(lam, lamp, mu)                          all copies of a sector from the closed-form matrices + the G8a forward recursion
CLI: python -X utf8 closed_forms_gl4.py a | b | words | gl3chain | mats | slices | all   (or: run <stage,...>)

Conventions (G4/G8b): GT pattern ((m11,), (m12, m22), (m13, m23, m33), lam); kappa = row 3 = (m13, m23, m33); l-coordinates y_i = kappa_i - i, Y_s = lam_s - s; d = |lam| - |kappa| (the depth);
nu = mu - lamp (the weight of W); matrices M[row][col] act on coordinate columns.
"""
import os, sys, time, itertools, random, math, pickle
from fractions import Fraction as Fr
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import mpmath as mp
import glN_reps as gr
import slice_engine_gl4 as se
import slice_ops_gl4 as so
import structure_gl4 as sg
import nested_gl3_gl4 as ng
import lattice_recursion_gl4 as lr

TABLES_PATH = os.path.join(HERE, 'closed_forms_gl4_tables.txt')
OUTPUT_PATH = os.path.join(HERE, 'closed_forms_gl4_output.txt')
CACHE_PATH = os.path.join(HERE, 'closed_forms_gl4_qcache.pkl')
_TAB = []
_OUT = []


def out(s=''):
    print(s, flush=True)
    _OUT.append(s)


def tab(s=''):
    _TAB.append(s)


def hdr(s):
    out('')
    out('=' * 150)
    out(s)
    out('=' * 150)
    tab('')
    tab('=' * 150)
    tab(s)
    tab('=' * 150)


def flush_tables(mode='w'):
    with open(TABLES_PATH, mode, encoding='utf-8') as f:
        f.write('\n'.join(_TAB) + '\n')
    with open(OUTPUT_PATH, mode, encoding='utf-8') as f:
        f.write('\n'.join(_OUT) + '\n')


def lstr(l):
    return '[' + ','.join(str(x) for x in l) + ']'


# ======================================================================================== D1: the nesting scalars
def ylist(kappa):
    return [kappa[i] - (i + 1) for i in range(3)]


def Ylist(lam):
    return [lam[s] - (s + 1) for s in range(4)]


def R_list(kappa, lam, kind='+'):
    """R_i = prod_{s=1..4} (y_i - Y_s) / prod_{j != i} (y_i - y_j)(y_i - y_j + 1)   (kind '+': 'raise kappa_i by E_{i4}, lower back')
       or     prod_s (y_i - Y_s - 1) / prod_{j != i} (y_i - y_j)(y_i - y_j - 1)       (kind '-': 'lower first').
    The term of a neighbour kappa +- e_i that is no dominant weight (pole of the denominator) is DROPPED (R_i := 0): the GT factor is a norm of an intermediate state that does not exist."""
    y = ylist(kappa)
    Y = Ylist(lam)
    sh = 0 if kind == '+' else 1
    out_ = []
    for i in range(3):
        den = 1
        for j in range(3):
            if j != i:
                den *= (y[i] - y[j]) * (y[i] - y[j] + (1 if kind == '+' else -1))
        if den == 0:
            out_.append(Fr(0))
            continue
        num = 1
        for s in range(4):
            num *= (y[i] - Y[s] - sh)
        out_.append(Fr(num, den))
    return out_


def a_of(kappa, lam):
    """a(kappa; lam) = 1 - sum_i R_i  (G8b's nesting scale of V3: J^{(3)}_kappa = a V^{(3)}_kappa + b3)"""
    return 1 - sum(R_list(kappa, lam))


def a4a_of(kappa, lam):
    """a_4a(kappa; lam) = d + sum_i R_i (2 y_i - |kappa| + 4),  d = |lam| - |kappa|   (J^{(4a)}_kappa = a_4a V^{(3)}_kappa + b4a)"""
    R = R_list(kappa, lam)
    y = ylist(kappa)
    return Fr(sum(lam) - sum(kappa)) + sum(R[i] * (2 * y[i] - sum(kappa) + 4) for i in range(3))


def a_data(lam, lamp, kappa, A=None, B=None, m11=None):
    """EXACT a, a4a, b3, b4a of the kappa-block of G3's operators from the two-point path (A,B),(A+1,B-1) (needs kappa1 > kappa2 > kappa3); the off-diagonal ratios of the lower and the upper entry must be equal"""
    R = gr.GLNRep.of(lam)
    k1, k2, k3 = kappa
    if A is None:
        A, B, m11 = k2, k2, k2
    pa = ((m11,), (A, B), tuple(kappa), tuple(lam))
    pb = ((m11,), (A + 1, B - 1), tuple(kappa), tuple(lam))
    ia, ib = R.index[pa], R.index[pb]
    nu = R.wt[ia]
    R3, W3, M3 = ng.gl3_V(tuple(kappa), tuple(lamp[:3]), tuple(nu[:3]))
    pos3 = {int(c): a for a, c in enumerate(W3)}
    j = [pos3[R3.index[tuple(R.pats[i][:3])]] for i in (ia, ib)]
    V = [[M3[r][c] for c in j] for r in j]
    res = {}
    for n in ('V3', 'V4a'):
        op = so.reduced_op(n, lamp)
        J = [[Fr(0)] * 2 for _ in range(2)]
        for cc, ic in enumerate((ia, ib)):
            img = op.apply(R, {ic: Fr(1)})
            for rr, ir in enumerate((ia, ib)):
                J[rr][cc] = img.get(ir, Fr(0))
        a1 = J[1][0] / V[1][0]
        a2 = J[0][1] / V[0][1]
        assert a1 == a2, (n, a1, a2)
        b = J[0][0] - a1 * V[0][0]
        assert J[1][1] - a1 * V[1][1] == b
        res[n] = (a1, b)
    return res, nu


def kappas_of(lam, strict=True):
    out_ = []
    for k1 in range(lam[1], lam[0] + 1):
        for k2 in range(lam[2], lam[1] + 1):
            for k3 in range(lam[3], lam[2] + 1):
                if not strict or k1 > k2 > k3:
                    out_.append((k1, k2, k3))
    return out_


# ======================================================================================== D2: the shifts b3, b4a (sparse tables in labels + invariants, see the end of the file)
_MON6 = [e for e in itertools.product(range(3), repeat=6) if sum(e) <= 2]      # monomials in (nu1, nu2, lam'1, lam'2, lam'3, lam'4)
_BCOEF_CACHE = {}


def _mono(vals, e):
    t = Fr(1)
    for v, k in zip(vals, e):
        if k:
            t *= v ** k
    return t


def b_coefficients(name, kappa, lam):
    """the 28 coefficients c_m(kappa, lam) of b = sum_m c_m nu^m  (m over _MON6): c_m = P_m(kappa, P1, P2, P3) + sum_i R_i S_m(y_i, e1(y), e2(y), P1, P2, P3), P_p = sum_s Y_s^p"""
    key = (name, tuple(kappa), tuple(lam))
    if key in _BCOEF_CACHE:
        return _BCOEF_CACHE[key]
    tab_ = B_TABLES[name]
    y = ylist(kappa)
    Y = Ylist(lam)
    P = [sum(Y), sum(v * v for v in Y), sum(v ** 3 for v in Y)]
    s1 = sum(y)
    e2 = y[0] * y[1] + y[0] * y[2] + y[1] * y[2]
    R = R_list(kappa, lam)
    pv = [Fr(kappa[0]), Fr(kappa[1]), Fr(kappa[2])] + [Fr(x) for x in P]
    out_ = {}
    for m, (pt, rt) in tab_.items():
        c = Fr(0)
        for e, co in pt:
            c += co * _mono(pv, e)
        for e, co in rt:
            sm = Fr(0)
            for i in range(3):
                if R[i] != 0:
                    sm += R[i] * _mono([Fr(y[i]), Fr(s1), Fr(e2)] + [Fr(x) for x in P], e)
            c += co * sm
        out_[m] = c
    _BCOEF_CACHE[key] = out_
    return out_


def b_of(name, kappa, lam, nu, lamp):
    """b3 (name 'V3') or b4a ('V4a') of the kappa-block, nu = weight of W (only nu1, nu2 enter), lamp = lam' (4 numbers)"""
    c = b_coefficients(name, kappa, lam)
    v = [Fr(nu[0]), Fr(nu[1])] + [Fr(x) for x in lamp]
    return sum(c[m] * _mono(v, m) for m in c)


# ======================================================================================== the gl(3) chain with general lam' (closed form)
def gl3_entries(kappa, q, t, s, m12):
    """closed forms of the gl(3) reduced operator V^{(3)}(q) = t S_{3,2} on V_kappa[(t, s - t, .)] in the GT basis, chain coordinate m1 = m12, m2 = s - m1 (m11 = t):
       up   = M[m1 + 1, m1] = -(k1-m1)(k2-m1-1)(k3-m1-2)(t+1-m2)(t-m1+q1-q2) / ((x+1)(x+2))
       down = M[m1 - 1, m1] = -(k2-m2)(t-m1)(k1-m2+1)(k3-m2-1)(t+1-m2+q1-q2) / (x(x+1))      (x = m1 - m2)
       alpha = M[m1, m1] = A0 + g xi + K'/xi,  xi = x(x+2),  g = 3s/16 - (k1+k2+k3)/8 + (q1+q2)/8 - q3/4 + 3/8,
               K' = (s-2k2)(s-2t)(s-2k1-2)(s-2k3+2)(2-s+2t+2q1-2q2)/16  (K'/xi := 0 at x = 0),  A0 = a cubic polynomial in (s, t, k, q) (table A0_TERMS)
    returns (alpha, up, down); up/down are None when the denominator vanishes (target not in the chain)"""
    k1, k2, k3 = kappa
    q1, q2, q3 = q
    m1 = m12
    m2 = s - m1
    x = m1 - m2
    up = down = None
    if (x + 1) * (x + 2) != 0:
        up = -Fr((k1 - m1) * (k2 - m1 - 1) * (k3 - m1 - 2) * (t + 1 - m2)) * (t - m1 + q1 - q2) / ((x + 1) * (x + 2))
    if x * (x + 1) != 0:
        down = -Fr((k2 - m2) * (t - m1) * (k1 - m2 + 1) * (k3 - m2 - 1)) * (t + 1 - m2 + q1 - q2) / (x * (x + 1))
    xi = x * (x + 2)
    vals = [Fr(s), Fr(t), Fr(k1), Fr(k2), Fr(k3), Fr(q1), Fr(q2), Fr(q3)]
    A0 = sum(co * _mono(vals, e) for e, co in A0_TERMS)
    g = Fr(3 * s, 16) - Fr(k1 + k2 + k3, 8) + (Fr(q1) + Fr(q2)) / 8 - Fr(q3) / 4 + Fr(3, 8)
    alpha = A0 + g * xi
    if xi != 0:
        Kp = Fr((s - 2 * k2) * (s - 2 * t) * (s - 2 * k1 - 2) * (s - 2 * k3 + 2)) * (2 - s + 2 * t + 2 * Fr(q1) - 2 * Fr(q2)) / 16
        alpha += Kp / xi
    return alpha, up, down


def chain_path(kappa, t, s):
    """the m12 values of the path P_kappa of the weight space (m11 = t, m12 + m22 = s)"""
    return [m for m in range(kappa[1], kappa[0] + 1) if kappa[2] <= s - m <= kappa[1] and s - m <= t <= m]


def gl3_chain(kappa, lamp3, nu123):
    """closed-form gl(3) chain for general lamp3: (path, alpha list, up list, down list) with up[i] = M[i+1, i], down[i] = M[i-1, i] (0 at the ends), no representation matrix"""
    t = nu123[0]
    s = nu123[0] + nu123[1]
    P = chain_path(kappa, t, s)
    al, up, dn = [], [], []
    for i, m in enumerate(P):
        a_, u_, d_ = gl3_entries(kappa, lamp3, t, s, m)
        al.append(a_)
        up.append(u_ if i + 1 < len(P) else Fr(0))
        dn.append(d_ if i > 0 else Fr(0))
    return P, al, up, dn


# ======================================================================================== weight spaces from labels
def pattern_key(p):
    """the order of glN_reps / gl3_reps.gt_patterns: (row 3, row 2, row 1) lexicographic (itertools.product order)"""
    return (p[2], p[1], p[0])


def weight_patterns(lam, nu):
    """all GT patterns of V_lam of weight nu (labels only), in the basis order of GLNRep"""
    tot = sum(lam)
    if sum(nu) != tot:
        return []
    out_ = []
    m11 = nu[0]
    for m13 in range(lam[1], lam[0] + 1):
        for m23 in range(lam[2], lam[1] + 1):
            m33 = tot - nu[3] - m13 - m23
            if not (lam[3] <= m33 <= lam[2]):
                continue
            for m12 in range(m23, m13 + 1):
                m22 = nu[0] + nu[1] - m12
                if m33 <= m22 <= m23 and m22 <= m11 <= m12:
                    out_.append(((m11,), (m12, m22), (m13, m23, m33), tuple(lam)))
    out_.sort(key=pattern_key)
    return out_


def entries6(p):
    return (p[0][0], p[1][0], p[1][1], p[2][0], p[2][1], p[2][2])


def pat_of(e, lam):
    return ((e[0],), (e[1], e[2]), (e[3], e[4], e[5]), tuple(lam))


def pt_of(p, lam, lamp):
    return sg.pt_of_pattern(p, lam, lamp)


# ======================================================================================== the Delta m12 = 0 classes (G4's polynomial Q, cached)
QCLASSES = [(0, 0, 0, 1, -1, 0), (0, 0, 0, 1, 0, -1), (0, 0, 0, -1, 1, 0), (0, 0, 0, 0, 1, -1), (0, 0, 0, -1, 0, 1), (0, 0, 0, 0, -1, 1)]
_QREC = {}


def _load_cache():
    if not _QREC and os.path.exists(CACHE_PATH):
        try:
            with open(CACHE_PATH, 'rb') as f:
                _QREC.update(pickle.load(f))
        except Exception:
            pass


def q_record(name, d):
    """slim record {forms, net, qdeg, Qexp, zero} of the class d (6-vector) of the operator `name` (G4's closed_form_one + reconstruct_q), cached in CACHE_PATH"""
    _load_cache()
    key = (name, tuple(d))
    if key not in _QREC:
        rec = sg.closed_form_one(name, tuple(d))
        if not rec['zero'] and rec['qdeg'] > 0:
            sg.reconstruct_q(rec)
        slim = dict(zero=rec['zero'])
        if not rec['zero']:
            slim.update(forms=rec['forms'], net=rec['net'], qdeg=rec['qdeg'], qaxis=rec['qaxis'], Qexp=rec.get('Qexp'), kappa=rec.get('kappa'))
        _QREC[key] = slim
        with open(CACHE_PATH, 'wb') as f:
            pickle.dump(_QREC, f)
    return _QREC[key]


def _lagr(xs, ys_, x0):
    tot = Fr(0)
    for i in range(len(xs)):
        term = ys_[i]
        for j in range(len(xs)):
            if j != i:
                term = term * (x0 - xs[j]) / (xs[i] - xs[j])
        tot += term
    return tot


def eval_qrec(rec, pt, rng):
    """exact value of the class: Q~(pt) * prod (net forms)^e, with the exact Laurent handling of removable poles along a random line (G4's eval_class, on the slim record)"""
    if rec['zero']:
        return Fr(0)
    qd = rec['qdeg']
    net = rec['net']
    if all(sg.form_value_f(f, pt) != 0 for f, e in net if e < 0):
        v = Fr(1)
        for f, e in net:
            v *= sg.form_value_f(f, pt) ** e
        return (sg.eval_Qexp(rec, pt) if qd else rec['kappa']) * v
    allforms = [f for f, e in net] + [(c0, items) for c0, items, m in rec['forms']]
    for _try in range(500):
        v = [Fr(rng.randint(1, 13)) for _ in range(14)]
        if all(sum(c * v[i] for i, c in f[1]) != 0 for f in allforms):
            break
    if qd == 0:
        qco = [Fr(rec['kappa'])]
    else:
        ss = [Fr(i) + Fr(1, 10007) for i in range(1, qd + 3)]
        vals = [sg.eval_Qexp(rec, [pt[i] + sx * v[i] for i in range(14)]) for sx in ss]
        qco = sg._poly_coeffs(ss[:qd + 1], vals[:qd + 1])
    pw = 0
    const = Fr(1)
    for f, e in net:
        fv = sg.form_value_f(f, pt)
        if fv == 0:
            lv = sum(c * v[i] for i, c in f[1])
            pw += e
            const *= lv ** e
        else:
            const *= fv ** e
    for j, a in enumerate(qco):
        if j + pw < 0 and a != 0:
            raise ZeroDivisionError('genuine pole')
    j0 = -pw
    if j0 < 0:
        return Fr(0)
    return (qco[j0] if j0 < len(qco) else Fr(0)) * const


# ======================================================================================== the closed-form matrices on a weight space
_CLOSED_OPS = ('V3', 'V4a')
_RNG = random.Random(12345)


def closed_form_mats(lam, lamp, mu, ops=_CLOSED_OPS, pats=None):
    """the w x w matrices (list of rows, M[row][col]) of V3 and V4a on W = V_lam[mu - lamp] from LABELS ALONE (no GLNRep, no representation matrix):
       diagonal                       a V^(3)_ii + b                      (D1, D2, gl3_entries),
       within a path (Dm12 = +-1)     a * (gl(3) chain entry up/down),
       12 hopping classes             structure_gl4.universal_F,
       6 classes with Dm12 = 0        G4's cubic Q x linear factors (q_record).
    V4b is NOT part of the closed forms (see RESULTS-G8d section 7).  Returns (pats, {name: matrix})."""
    lam, lamp, mu = tuple(lam), tuple(lamp), tuple(mu)
    nu = tuple(mu[i] - lamp[i] for i in range(4))
    if pats is None:
        pats = weight_patterns(lam, nu)
    w = len(pats)
    idx = {p: i for i, p in enumerate(pats)}
    M = {n: [[Fr(0)] * w for _ in range(w)] for n in ops}
    qrec = {(n, d): q_record(n, d) for n in ops for d in QCLASSES}
    hops = sg.f_classes()
    q3 = tuple(lamp[:3])
    ab = {}
    for c, p in enumerate(pats):
        e = entries6(p)
        kappa = p[2]
        t = e[0]
        s = e[1] + e[2]
        pt = None
        if kappa not in ab:
            ab[kappa] = {n: (a_of(kappa, lam) if n == 'V3' else a4a_of(kappa, lam)) for n in ops}
        al, up, dn = gl3_entries(kappa, q3, t, s, e[1])
        for n in ops:
            M[n][c][c] = ab[kappa][n] * al + b_of(n, kappa, lam, nu, lamp)
        # within path
        for dm, val in ((1, up), (-1, dn)):
            tg = (e[0], e[1] + dm, e[2] - dm, e[3], e[4], e[5])
            r = idx.get(pat_of(tg, lam))
            if r is not None:
                for n in ops:
                    M[n][r][c] = ab[kappa][n] * val
        # hopping and Delta m12 = 0 classes
        pt = pt_of(p, lam, lamp)
        for d in hops:
            tg = tuple(a + b for a, b in zip(e, d))
            r = idx.get(pat_of(tg, lam))
            if r is None:
                continue
            for n in ops:
                M[n][r][c] = sg.universal_F(n, pt, d)
        for d in QCLASSES:
            tg = tuple(a + b for a, b in zip(e, d))
            r = idx.get(pat_of(tg, lam))
            if r is None:
                continue
            for n in ops:
                M[n][r][c] = eval_qrec(qrec[(n, d)], pt, _RNG)
    return pats, M


def g3_mats(lam, lamp, mu, ops=_CLOSED_OPS):
    """G3's exact matrices on W (representation matrices) for the comparison; returns (pats, {name: matrix}) in the order of W"""
    R = gr.GLNRep.of(tuple(lam))
    nu = tuple(mu[i] - lamp[i] for i in range(4))
    W = R.wsp.get(nu)
    if W is None:
        return [], {n: [] for n in ops}
    mats = so.weight_space_mats(R, {n: so.reduced_op(n, tuple(lamp)) for n in ops}, W)
    return [R.pats[int(c)] for c in W], mats


# ==== GENERATED TABLES (do not edit by hand; produced by the fit of RESULTS-G8d section 3) ====
B_TABLES = {
    'V3': {
        (0, 0, 0, 0, 0, 0): ([((0, 0, 0, 0, 0, 0), Fr(-55, 2)), ((0, 0, 0, 0, 0, 1), Fr(-1, 3)), ((0, 0, 0, 0, 1, 0), Fr(7, 4)), ((0, 0, 0, 1, 0, 0), Fr(-10, 1)), ((0, 0, 0, 1, 1, 0), Fr(1, 2)), ((0, 0, 0, 2, 0, 0), Fr(-7, 4)), ((0, 0, 0, 3, 0, 0), Fr(-1, 6)), ((0, 0, 1, 0, 0, 0), Fr(9, 2)), ((0, 0, 1, 1, 0, 0), Fr(1, 2)), ((0, 0, 2, 0, 0, 0), Fr(-1, 2)), ((0, 1, 0, 0, 0, 0), Fr(6, 1)), ((0, 1, 0, 1, 0, 0), Fr(1, 2)), ((0, 1, 1, 0, 0, 0), Fr(-5, 2)), ((0, 2, 0, 0, 0, 0), Fr(-1, 2)), ((1, 0, 0, 0, 0, 0), Fr(19, 2)), ((1, 0, 0, 1, 0, 0), Fr(1, 2)), ((1, 0, 1, 0, 0, 0), Fr(-7, 2)), ((1, 1, 0, 0, 0, 0), Fr(-9, 2)), ((1, 1, 1, 0, 0, 0), Fr(1, 1)), ((2, 0, 0, 0, 0, 0), Fr(-1, 2))], [((0, 0, 0, 0, 0, 0), Fr(-45, 2)), ((0, 0, 1, 0, 0, 0), Fr(-9, 2)), ((0, 1, 0, 0, 0, 0), Fr(-15, 2)), ((1, 0, 0, 0, 0, 0), Fr(-23, 2)), ((1, 0, 1, 0, 0, 0), Fr(-1, 1)), ((2, 0, 0, 0, 0, 0), Fr(-9, 2))]),
        (0, 0, 0, 0, 0, 1): ([((0, 0, 1, 0, 0, 0), Fr(2, 1)), ((0, 1, 0, 0, 0, 0), Fr(1, 1)), ((0, 1, 1, 0, 0, 0), Fr(1, 1)), ((1, 0, 1, 0, 0, 0), Fr(1, 1)), ((1, 1, 0, 0, 0, 0), Fr(1, 1))], []),
        (0, 0, 0, 0, 0, 2): ([], []),
        (0, 0, 0, 0, 1, 0): ([((0, 0, 0, 0, 0, 0), Fr(45, 1)), ((0, 0, 0, 0, 1, 0), Fr(-1, 2)), ((0, 0, 0, 1, 0, 0), Fr(8, 1)), ((0, 0, 0, 2, 0, 0), Fr(1, 2)), ((0, 0, 1, 0, 0, 0), Fr(-11, 1)), ((0, 0, 1, 1, 0, 0), Fr(-1, 1)), ((0, 0, 2, 0, 0, 0), Fr(1, 1)), ((0, 1, 0, 0, 0, 0), Fr(-10, 1)), ((0, 1, 0, 1, 0, 0), Fr(-1, 1)), ((0, 1, 1, 0, 0, 0), Fr(1, 1)), ((0, 2, 0, 0, 0, 0), Fr(1, 1)), ((1, 0, 0, 0, 0, 0), Fr(-9, 1)), ((1, 0, 0, 1, 0, 0), Fr(-1, 1)), ((1, 0, 1, 0, 0, 0), Fr(1, 1)), ((1, 1, 0, 0, 0, 0), Fr(1, 1)), ((2, 0, 0, 0, 0, 0), Fr(1, 1))], [((0, 0, 0, 0, 0, 0), Fr(-6, 1)), ((1, 0, 0, 0, 0, 0), Fr(-5, 1)), ((2, 0, 0, 0, 0, 0), Fr(-1, 1))]),
        (0, 0, 0, 0, 1, 1): ([], []),
        (0, 0, 0, 0, 2, 0): ([], []),
        (0, 0, 0, 1, 0, 0): ([((0, 0, 0, 0, 0, 0), Fr(65, 2)), ((0, 0, 0, 0, 1, 0), Fr(-1, 4)), ((0, 0, 0, 1, 0, 0), Fr(5, 1)), ((0, 0, 0, 2, 0, 0), Fr(1, 4)), ((0, 0, 1, 0, 0, 0), Fr(7, 2)), ((0, 0, 1, 1, 0, 0), Fr(1, 2)), ((0, 0, 2, 0, 0, 0), Fr(-1, 2)), ((0, 1, 0, 0, 0, 0), Fr(4, 1)), ((0, 1, 0, 1, 0, 0), Fr(1, 2)), ((0, 1, 1, 0, 0, 0), Fr(-3, 2)), ((0, 2, 0, 0, 0, 0), Fr(-1, 2)), ((1, 0, 0, 0, 0, 0), Fr(9, 2)), ((1, 0, 0, 1, 0, 0), Fr(1, 2)), ((1, 0, 1, 0, 0, 0), Fr(-3, 2)), ((1, 1, 0, 0, 0, 0), Fr(-3, 2)), ((2, 0, 0, 0, 0, 0), Fr(-1, 2))], [((0, 0, 0, 0, 0, 0), Fr(25, 2)), ((0, 0, 1, 0, 0, 0), Fr(1, 2)), ((0, 1, 0, 0, 0, 0), Fr(5, 2)), ((1, 0, 0, 0, 0, 0), Fr(5, 2)), ((2, 0, 0, 0, 0, 0), Fr(1, 2))]),
        (0, 0, 0, 1, 0, 1): ([], []),
        (0, 0, 0, 1, 1, 0): ([], []),
        (0, 0, 0, 2, 0, 0): ([], []),
        (0, 0, 1, 0, 0, 0): ([((0, 0, 0, 0, 0, 0), Fr(65, 2)), ((0, 0, 0, 0, 1, 0), Fr(-1, 4)), ((0, 0, 0, 1, 0, 0), Fr(5, 1)), ((0, 0, 0, 2, 0, 0), Fr(1, 4)), ((0, 0, 1, 0, 0, 0), Fr(7, 2)), ((0, 0, 1, 1, 0, 0), Fr(1, 2)), ((0, 0, 2, 0, 0, 0), Fr(-1, 2)), ((0, 1, 0, 0, 0, 0), Fr(4, 1)), ((0, 1, 0, 1, 0, 0), Fr(1, 2)), ((0, 1, 1, 0, 0, 0), Fr(-3, 2)), ((0, 2, 0, 0, 0, 0), Fr(-1, 2)), ((1, 0, 0, 0, 0, 0), Fr(9, 2)), ((1, 0, 0, 1, 0, 0), Fr(1, 2)), ((1, 0, 1, 0, 0, 0), Fr(-3, 2)), ((1, 1, 0, 0, 0, 0), Fr(-3, 2)), ((2, 0, 0, 0, 0, 0), Fr(-1, 2))], [((0, 0, 0, 0, 0, 0), Fr(25, 2)), ((0, 0, 1, 0, 0, 0), Fr(1, 2)), ((0, 1, 0, 0, 0, 0), Fr(5, 2)), ((1, 0, 0, 0, 0, 0), Fr(5, 2)), ((2, 0, 0, 0, 0, 0), Fr(1, 2))]),
        (0, 0, 1, 0, 0, 1): ([], []),
        (0, 0, 1, 0, 1, 0): ([], []),
        (0, 0, 1, 1, 0, 0): ([], []),
        (0, 0, 2, 0, 0, 0): ([], []),
        (0, 1, 0, 0, 0, 0): ([((0, 0, 0, 0, 0, 0), Fr(-45, 2)), ((0, 0, 0, 0, 1, 0), Fr(1, 4)), ((0, 0, 0, 1, 0, 0), Fr(-4, 1)), ((0, 0, 0, 2, 0, 0), Fr(-1, 4)), ((0, 0, 1, 0, 0, 0), Fr(31, 2)), ((0, 0, 1, 1, 0, 0), Fr(3, 2)), ((0, 0, 2, 0, 0, 0), Fr(-3, 2)), ((0, 1, 0, 0, 0, 0), Fr(15, 1)), ((0, 1, 0, 1, 0, 0), Fr(3, 2)), ((0, 1, 1, 0, 0, 0), Fr(-5, 2)), ((0, 2, 0, 0, 0, 0), Fr(-3, 2)), ((1, 0, 0, 0, 0, 0), Fr(29, 2)), ((1, 0, 0, 1, 0, 0), Fr(3, 2)), ((1, 0, 1, 0, 0, 0), Fr(-5, 2)), ((1, 1, 0, 0, 0, 0), Fr(-5, 2)), ((2, 0, 0, 0, 0, 0), Fr(-3, 2))], [((0, 0, 0, 0, 0, 0), Fr(31, 2)), ((0, 0, 1, 0, 0, 0), Fr(1, 2)), ((0, 1, 0, 0, 0, 0), Fr(5, 2)), ((1, 0, 0, 0, 0, 0), Fr(13, 2)), ((2, 0, 0, 0, 0, 0), Fr(3, 2))]),
        (0, 1, 0, 0, 0, 1): ([], []),
        (0, 1, 0, 0, 1, 0): ([((0, 0, 0, 0, 0, 0), Fr(10, 1)), ((0, 0, 0, 1, 0, 0), Fr(1, 1)), ((0, 0, 1, 0, 0, 0), Fr(-1, 1)), ((0, 1, 0, 0, 0, 0), Fr(-1, 1)), ((1, 0, 0, 0, 0, 0), Fr(-1, 1))], [((0, 0, 0, 0, 0, 0), Fr(3, 1)), ((1, 0, 0, 0, 0, 0), Fr(1, 1))]),
        (0, 1, 0, 1, 0, 0): ([((0, 0, 0, 0, 0, 0), Fr(-10, 1)), ((0, 0, 0, 1, 0, 0), Fr(-1, 1)), ((0, 0, 1, 0, 0, 0), Fr(1, 1)), ((0, 1, 0, 0, 0, 0), Fr(1, 1)), ((1, 0, 0, 0, 0, 0), Fr(1, 1))], [((0, 0, 0, 0, 0, 0), Fr(-3, 1)), ((1, 0, 0, 0, 0, 0), Fr(-1, 1))]),
        (0, 1, 1, 0, 0, 0): ([], []),
        (0, 2, 0, 0, 0, 0): ([((0, 0, 0, 0, 0, 0), Fr(-10, 1)), ((0, 0, 0, 1, 0, 0), Fr(-1, 1)), ((0, 0, 1, 0, 0, 0), Fr(1, 1)), ((0, 1, 0, 0, 0, 0), Fr(1, 1)), ((1, 0, 0, 0, 0, 0), Fr(1, 1))], [((0, 0, 0, 0, 0, 0), Fr(-3, 1)), ((1, 0, 0, 0, 0, 0), Fr(-1, 1))]),
        (1, 0, 0, 0, 0, 0): ([((0, 0, 0, 0, 0, 0), Fr(-65, 2)), ((0, 0, 0, 0, 1, 0), Fr(1, 4)), ((0, 0, 0, 1, 0, 0), Fr(-5, 1)), ((0, 0, 0, 2, 0, 0), Fr(-1, 4)), ((0, 0, 1, 0, 0, 0), Fr(33, 2)), ((0, 0, 1, 1, 0, 0), Fr(3, 2)), ((0, 0, 2, 0, 0, 0), Fr(-3, 2)), ((0, 1, 0, 0, 0, 0), Fr(16, 1)), ((0, 1, 0, 1, 0, 0), Fr(3, 2)), ((0, 1, 1, 0, 0, 0), Fr(-5, 2)), ((0, 2, 0, 0, 0, 0), Fr(-3, 2)), ((1, 0, 0, 0, 0, 0), Fr(31, 2)), ((1, 0, 0, 1, 0, 0), Fr(3, 2)), ((1, 0, 1, 0, 0, 0), Fr(-5, 2)), ((1, 1, 0, 0, 0, 0), Fr(-5, 2)), ((2, 0, 0, 0, 0, 0), Fr(-3, 2))], [((0, 0, 0, 0, 0, 0), Fr(25, 2)), ((0, 0, 1, 0, 0, 0), Fr(1, 2)), ((0, 1, 0, 0, 0, 0), Fr(5, 2)), ((1, 0, 0, 0, 0, 0), Fr(11, 2)), ((2, 0, 0, 0, 0, 0), Fr(3, 2))]),
        (1, 0, 0, 0, 0, 1): ([], []),
        (1, 0, 0, 0, 1, 0): ([((0, 0, 0, 0, 0, 0), Fr(10, 1)), ((0, 0, 0, 1, 0, 0), Fr(1, 1)), ((0, 0, 1, 0, 0, 0), Fr(-1, 1)), ((0, 1, 0, 0, 0, 0), Fr(-1, 1)), ((1, 0, 0, 0, 0, 0), Fr(-1, 1))], [((0, 0, 0, 0, 0, 0), Fr(3, 1)), ((1, 0, 0, 0, 0, 0), Fr(1, 1))]),
        (1, 0, 0, 1, 0, 0): ([], []),
        (1, 0, 1, 0, 0, 0): ([((0, 0, 0, 0, 0, 0), Fr(-10, 1)), ((0, 0, 0, 1, 0, 0), Fr(-1, 1)), ((0, 0, 1, 0, 0, 0), Fr(1, 1)), ((0, 1, 0, 0, 0, 0), Fr(1, 1)), ((1, 0, 0, 0, 0, 0), Fr(1, 1))], [((0, 0, 0, 0, 0, 0), Fr(-3, 1)), ((1, 0, 0, 0, 0, 0), Fr(-1, 1))]),
        (1, 1, 0, 0, 0, 0): ([((0, 0, 0, 0, 0, 0), Fr(-10, 1)), ((0, 0, 0, 1, 0, 0), Fr(-1, 1)), ((0, 0, 1, 0, 0, 0), Fr(1, 1)), ((0, 1, 0, 0, 0, 0), Fr(1, 1)), ((1, 0, 0, 0, 0, 0), Fr(1, 1))], [((0, 0, 0, 0, 0, 0), Fr(-3, 1)), ((1, 0, 0, 0, 0, 0), Fr(-1, 1))]),
        (2, 0, 0, 0, 0, 0): ([((0, 0, 0, 0, 0, 0), Fr(-10, 1)), ((0, 0, 0, 1, 0, 0), Fr(-1, 1)), ((0, 0, 1, 0, 0, 0), Fr(1, 1)), ((0, 1, 0, 0, 0, 0), Fr(1, 1)), ((1, 0, 0, 0, 0, 0), Fr(1, 1))], [((0, 0, 0, 0, 0, 0), Fr(-3, 1)), ((1, 0, 0, 0, 0, 0), Fr(-1, 1))]),
    },
    'V4a': {
        (0, 0, 0, 0, 0, 0): ([((0, 0, 0, 0, 0, 0), Fr(2015, 12)), ((0, 0, 0, 0, 0, 1), Fr(1, 6)), ((0, 0, 0, 0, 1, 0), Fr(-95, 24)), ((0, 0, 0, 1, 0, 0), Fr(42, 1)), ((0, 0, 0, 1, 1, 0), Fr(-1, 4)), ((0, 0, 0, 2, 0, 0), Fr(95, 24)), ((0, 0, 0, 3, 0, 0), Fr(1, 12)), ((0, 0, 1, 0, 0, 0), Fr(-381, 4)), ((0, 0, 1, 0, 1, 0), Fr(3, 4)), ((0, 0, 1, 1, 0, 0), Fr(-185, 12)), ((0, 0, 1, 2, 0, 0), Fr(-3, 4)), ((0, 0, 2, 0, 0, 0), Fr(13, 4)), ((0, 0, 2, 0, 1, 0), Fr(1, 3)), ((0, 0, 2, 1, 0, 0), Fr(-5, 2)), ((0, 0, 2, 2, 0, 0), Fr(-1, 3)), ((0, 0, 3, 0, 0, 0), Fr(9, 2)), ((0, 0, 3, 1, 0, 0), Fr(2, 3)), ((0, 0, 4, 0, 0, 0), Fr(-2, 3)), ((0, 1, 0, 0, 0, 0), Fr(-649, 6)), ((0, 1, 0, 0, 1, 0), Fr(17, 12)), ((0, 1, 0, 1, 0, 0), Fr(-263, 12)), ((0, 1, 0, 2, 0, 0), Fr(-17, 12)), ((0, 1, 1, 0, 0, 0), Fr(341, 12)), ((0, 1, 1, 1, 0, 0), Fr(17, 6)), ((0, 1, 2, 0, 0, 0), Fr(5, 2)), ((0, 1, 2, 1, 0, 0), Fr(2, 3)), ((0, 1, 3, 0, 0, 0), Fr(-2, 3)), ((0, 2, 0, 0, 0, 0), Fr(131, 12)), ((0, 2, 0, 0, 1, 0), Fr(1, 3)), ((0, 2, 0, 1, 0, 0), Fr(-7, 6)), ((0, 2, 0, 2, 0, 0), Fr(-1, 3)), ((0, 2, 1, 0, 0, 0), Fr(19, 6)), ((0, 2, 1, 1, 0, 0), Fr(2, 3)), ((0, 2, 2, 0, 0, 0), Fr(-4, 3)), ((0, 3, 0, 0, 0, 0), Fr(5, 2)), ((0, 3, 0, 1, 0, 0), Fr(2, 3)), ((0, 3, 1, 0, 0, 0), Fr(-2, 3)), ((0, 4, 0, 0, 0, 0), Fr(-2, 3)), ((1, 0, 0, 0, 0, 0), Fr(-1465, 12)), ((1, 0, 0, 0, 1, 0), Fr(25, 12)), ((1, 0, 0, 1, 0, 0), Fr(-341, 12)), ((1, 0, 0, 2, 0, 0), Fr(-25, 12)), ((1, 0, 1, 0, 0, 0), Fr(455, 12)), ((1, 0, 1, 1, 0, 0), Fr(25, 6)), ((1, 0, 2, 0, 0, 0), Fr(1, 2)), ((1, 0, 2, 1, 0, 0), Fr(2, 3)), ((1, 0, 3, 0, 0, 0), Fr(-2, 3)), ((1, 1, 0, 0, 0, 0), Fr(505, 12)), ((1, 1, 0, 1, 0, 0), Fr(11, 2)), ((1, 1, 1, 0, 0, 0), Fr(-9, 2)), ((1, 1, 2, 0, 0, 0), Fr(-2, 3)), ((1, 2, 0, 0, 0, 0), Fr(-5, 6)), ((1, 2, 0, 1, 0, 0), Fr(2, 3)), ((1, 2, 1, 0, 0, 0), Fr(-2, 3)), ((1, 3, 0, 0, 0, 0), Fr(-2, 3)), ((2, 0, 0, 0, 0, 0), Fr(191, 12)), ((2, 0, 0, 0, 1, 0), Fr(1, 3)), ((2, 0, 0, 1, 0, 0), Fr(1, 6)), ((2, 0, 0, 2, 0, 0), Fr(-1, 3)), ((2, 0, 1, 0, 0, 0), Fr(11, 6)), ((2, 0, 1, 1, 0, 0), Fr(2, 3)), ((2, 0, 2, 0, 0, 0), Fr(-4, 3)), ((2, 1, 0, 0, 0, 0), Fr(-1, 6)), ((2, 1, 0, 1, 0, 0), Fr(2, 3)), ((2, 1, 1, 0, 0, 0), Fr(-2, 3)), ((2, 2, 0, 0, 0, 0), Fr(-4, 3)), ((3, 0, 0, 0, 0, 0), Fr(1, 2)), ((3, 0, 0, 1, 0, 0), Fr(2, 3)), ((3, 0, 1, 0, 0, 0), Fr(-2, 3)), ((3, 1, 0, 0, 0, 0), Fr(-2, 3)), ((4, 0, 0, 0, 0, 0), Fr(-2, 3))], [((0, 0, 0, 0, 0, 0), Fr(-169, 4)), ((0, 0, 1, 0, 0, 0), Fr(-115, 12)), ((0, 0, 2, 0, 0, 0), Fr(-10, 3)), ((0, 1, 0, 0, 0, 0), Fr(-159, 4)), ((0, 1, 1, 0, 0, 0), Fr(-2, 3)), ((0, 2, 0, 0, 0, 0), Fr(-22, 3)), ((0, 2, 1, 0, 0, 0), Fr(2, 3)), ((0, 3, 0, 0, 0, 0), Fr(-2, 3)), ((1, 0, 0, 0, 0, 0), Fr(123, 4)), ((1, 0, 1, 0, 0, 0), Fr(-5, 2)), ((1, 1, 1, 0, 0, 0), Fr(20, 3)), ((1, 3, 0, 0, 0, 0), Fr(-4, 3)), ((2, 0, 0, 0, 0, 0), Fr(91, 4)), ((2, 0, 1, 0, 0, 0), Fr(-9, 1)), ((3, 0, 0, 0, 0, 0), Fr(15, 2)), ((3, 1, 0, 0, 0, 0), Fr(5, 1)), ((4, 0, 0, 0, 0, 0), Fr(-3, 1))]),
        (0, 0, 0, 0, 0, 1): ([((0, 0, 1, 0, 0, 0), Fr(2, 1)), ((0, 1, 1, 0, 0, 0), Fr(2, 1)), ((1, 0, 1, 0, 0, 0), Fr(1, 1)), ((1, 1, 1, 0, 0, 0), Fr(1, 1))], []),
        (0, 0, 0, 0, 0, 2): ([], []),
        (0, 0, 0, 0, 1, 0): ([((0, 0, 0, 0, 0, 0), Fr(385, 4)), ((0, 0, 0, 0, 0, 1), Fr(1, 6)), ((0, 0, 0, 0, 1, 0), Fr(-17, 8)), ((0, 0, 0, 1, 0, 0), Fr(22, 1)), ((0, 0, 0, 1, 1, 0), Fr(-1, 4)), ((0, 0, 0, 2, 0, 0), Fr(17, 8)), ((0, 0, 0, 3, 0, 0), Fr(1, 12)), ((0, 0, 1, 0, 0, 0), Fr(-157, 4)), ((0, 0, 1, 0, 1, 0), Fr(1, 4)), ((0, 0, 1, 1, 0, 0), Fr(-23, 4)), ((0, 0, 1, 2, 0, 0), Fr(-1, 4)), ((0, 0, 2, 0, 0, 0), Fr(29, 4)), ((0, 0, 2, 1, 0, 0), Fr(1, 2)), ((0, 0, 3, 0, 0, 0), Fr(-1, 2)), ((0, 1, 0, 0, 0, 0), Fr(-65, 2)), ((0, 1, 0, 0, 1, 0), Fr(1, 4)), ((0, 1, 0, 1, 0, 0), Fr(-21, 4)), ((0, 1, 0, 2, 0, 0), Fr(-1, 4)), ((0, 1, 1, 0, 0, 0), Fr(27, 4)), ((0, 1, 1, 1, 0, 0), Fr(1, 2)), ((0, 1, 2, 0, 0, 0), Fr(-1, 2)), ((0, 2, 0, 0, 0, 0), Fr(25, 4)), ((0, 2, 0, 1, 0, 0), Fr(1, 2)), ((0, 2, 1, 0, 0, 0), Fr(-1, 2)), ((0, 3, 0, 0, 0, 0), Fr(-1, 2)), ((1, 0, 0, 0, 0, 0), Fr(-107, 4)), ((1, 0, 0, 0, 1, 0), Fr(1, 4)), ((1, 0, 0, 1, 0, 0), Fr(-19, 4)), ((1, 0, 0, 2, 0, 0), Fr(-1, 4)), ((1, 0, 1, 0, 0, 0), Fr(25, 4)), ((1, 0, 1, 1, 0, 0), Fr(1, 2)), ((1, 0, 2, 0, 0, 0), Fr(-1, 2)), ((1, 1, 0, 0, 0, 0), Fr(23, 4)), ((1, 1, 0, 1, 0, 0), Fr(1, 2)), ((1, 1, 1, 0, 0, 0), Fr(-1, 2)), ((1, 2, 0, 0, 0, 0), Fr(-1, 2)), ((2, 0, 0, 0, 0, 0), Fr(21, 4)), ((2, 0, 0, 1, 0, 0), Fr(1, 2)), ((2, 0, 1, 0, 0, 0), Fr(-1, 2)), ((2, 1, 0, 0, 0, 0), Fr(-1, 2)), ((3, 0, 0, 0, 0, 0), Fr(-1, 2))], [((0, 0, 0, 0, 0, 0), Fr(-57, 4)), ((0, 0, 1, 0, 0, 0), Fr(-9, 4)), ((0, 1, 0, 0, 0, 0), Fr(-15, 4)), ((1, 0, 0, 0, 0, 0), Fr(-21, 4)), ((1, 0, 1, 0, 0, 0), Fr(-1, 2)), ((2, 0, 0, 0, 0, 0), Fr(-1, 4)), ((3, 0, 0, 0, 0, 0), Fr(1, 2))]),
        (0, 0, 0, 0, 1, 1): ([], []),
        (0, 0, 0, 0, 2, 0): ([], []),
        (0, 0, 0, 1, 0, 0): ([((0, 0, 0, 0, 0, 0), Fr(815, 8)), ((0, 0, 0, 0, 0, 1), Fr(1, 12)), ((0, 0, 0, 0, 1, 0), Fr(-31, 16)), ((0, 0, 0, 1, 0, 0), Fr(45, 2)), ((0, 0, 0, 1, 1, 0), Fr(-1, 8)), ((0, 0, 0, 2, 0, 0), Fr(31, 16)), ((0, 0, 0, 3, 0, 0), Fr(1, 24)), ((0, 0, 1, 0, 0, 0), Fr(-111, 8)), ((0, 0, 1, 0, 1, 0), Fr(-1, 8)), ((0, 0, 1, 1, 0, 0), Fr(-5, 8)), ((0, 0, 1, 2, 0, 0), Fr(1, 8)), ((0, 0, 2, 0, 0, 0), Fr(-1, 8)), ((0, 0, 2, 1, 0, 0), Fr(-1, 4)), ((0, 0, 3, 0, 0, 0), Fr(1, 4)), ((0, 1, 0, 0, 0, 0), Fr(-35, 4)), ((0, 1, 0, 0, 1, 0), Fr(-1, 8)), ((0, 1, 0, 1, 0, 0), Fr(-3, 8)), ((0, 1, 0, 2, 0, 0), Fr(1, 8)), ((0, 1, 1, 0, 0, 0), Fr(-43, 8)), ((0, 1, 1, 1, 0, 0), Fr(-3, 4)), ((0, 1, 2, 0, 0, 0), Fr(3, 4)), ((0, 2, 0, 0, 0, 0), Fr(-1, 8)), ((0, 2, 0, 1, 0, 0), Fr(-1, 4)), ((0, 2, 1, 0, 0, 0), Fr(3, 4)), ((0, 3, 0, 0, 0, 0), Fr(1, 4)), ((1, 0, 0, 0, 0, 0), Fr(-33, 8)), ((1, 0, 0, 0, 1, 0), Fr(-1, 8)), ((1, 0, 0, 1, 0, 0), Fr(-1, 8)), ((1, 0, 0, 2, 0, 0), Fr(1, 8)), ((1, 0, 1, 0, 0, 0), Fr(-41, 8)), ((1, 0, 1, 1, 0, 0), Fr(-3, 4)), ((1, 0, 2, 0, 0, 0), Fr(3, 4)), ((1, 1, 0, 0, 0, 0), Fr(-39, 8)), ((1, 1, 0, 1, 0, 0), Fr(-3, 4)), ((1, 1, 1, 0, 0, 0), Fr(5, 4)), ((1, 2, 0, 0, 0, 0), Fr(3, 4)), ((2, 0, 0, 0, 0, 0), Fr(-1, 8)), ((2, 0, 0, 1, 0, 0), Fr(-1, 4)), ((2, 0, 1, 0, 0, 0), Fr(3, 4)), ((2, 1, 0, 0, 0, 0), Fr(3, 4)), ((3, 0, 0, 0, 0, 0), Fr(1, 4))], [((0, 0, 0, 0, 0, 0), Fr(189, 8)), ((0, 0, 1, 0, 0, 0), Fr(-3, 8)), ((0, 1, 0, 0, 0, 0), Fr(151, 8)), ((0, 1, 1, 0, 0, 0), Fr(1, 2)), ((0, 2, 0, 0, 0, 0), Fr(5, 2)), ((1, 0, 0, 0, 0, 0), Fr(-191, 8)), ((1, 0, 1, 0, 0, 0), Fr(-3, 4)), ((2, 0, 0, 0, 0, 0), Fr(-59, 8)), ((3, 0, 0, 0, 0, 0), Fr(-1, 4))]),
        (0, 0, 0, 1, 0, 1): ([], []),
        (0, 0, 0, 1, 1, 0): ([], []),
        (0, 0, 0, 2, 0, 0): ([], []),
        (0, 0, 1, 0, 0, 0): ([((0, 0, 0, 0, 0, 0), Fr(815, 8)), ((0, 0, 0, 0, 0, 1), Fr(1, 12)), ((0, 0, 0, 0, 1, 0), Fr(-31, 16)), ((0, 0, 0, 1, 0, 0), Fr(45, 2)), ((0, 0, 0, 1, 1, 0), Fr(-1, 8)), ((0, 0, 0, 2, 0, 0), Fr(31, 16)), ((0, 0, 0, 3, 0, 0), Fr(1, 24)), ((0, 0, 1, 0, 0, 0), Fr(-111, 8)), ((0, 0, 1, 0, 1, 0), Fr(-1, 8)), ((0, 0, 1, 1, 0, 0), Fr(-5, 8)), ((0, 0, 1, 2, 0, 0), Fr(1, 8)), ((0, 0, 2, 0, 0, 0), Fr(-1, 8)), ((0, 0, 2, 1, 0, 0), Fr(-1, 4)), ((0, 0, 3, 0, 0, 0), Fr(1, 4)), ((0, 1, 0, 0, 0, 0), Fr(-35, 4)), ((0, 1, 0, 0, 1, 0), Fr(-1, 8)), ((0, 1, 0, 1, 0, 0), Fr(-3, 8)), ((0, 1, 0, 2, 0, 0), Fr(1, 8)), ((0, 1, 1, 0, 0, 0), Fr(-43, 8)), ((0, 1, 1, 1, 0, 0), Fr(-3, 4)), ((0, 1, 2, 0, 0, 0), Fr(3, 4)), ((0, 2, 0, 0, 0, 0), Fr(-1, 8)), ((0, 2, 0, 1, 0, 0), Fr(-1, 4)), ((0, 2, 1, 0, 0, 0), Fr(3, 4)), ((0, 3, 0, 0, 0, 0), Fr(1, 4)), ((1, 0, 0, 0, 0, 0), Fr(-33, 8)), ((1, 0, 0, 0, 1, 0), Fr(-1, 8)), ((1, 0, 0, 1, 0, 0), Fr(-1, 8)), ((1, 0, 0, 2, 0, 0), Fr(1, 8)), ((1, 0, 1, 0, 0, 0), Fr(-41, 8)), ((1, 0, 1, 1, 0, 0), Fr(-3, 4)), ((1, 0, 2, 0, 0, 0), Fr(3, 4)), ((1, 1, 0, 0, 0, 0), Fr(-39, 8)), ((1, 1, 0, 1, 0, 0), Fr(-3, 4)), ((1, 1, 1, 0, 0, 0), Fr(5, 4)), ((1, 2, 0, 0, 0, 0), Fr(3, 4)), ((2, 0, 0, 0, 0, 0), Fr(-1, 8)), ((2, 0, 0, 1, 0, 0), Fr(-1, 4)), ((2, 0, 1, 0, 0, 0), Fr(3, 4)), ((2, 1, 0, 0, 0, 0), Fr(3, 4)), ((3, 0, 0, 0, 0, 0), Fr(1, 4))], [((0, 0, 0, 0, 0, 0), Fr(189, 8)), ((0, 0, 1, 0, 0, 0), Fr(-3, 8)), ((0, 1, 0, 0, 0, 0), Fr(151, 8)), ((0, 1, 1, 0, 0, 0), Fr(1, 2)), ((0, 2, 0, 0, 0, 0), Fr(5, 2)), ((1, 0, 0, 0, 0, 0), Fr(-191, 8)), ((1, 0, 1, 0, 0, 0), Fr(-3, 4)), ((2, 0, 0, 0, 0, 0), Fr(-59, 8)), ((3, 0, 0, 0, 0, 0), Fr(-1, 4))]),
        (0, 0, 1, 0, 0, 1): ([], []),
        (0, 0, 1, 0, 1, 0): ([], []),
        (0, 0, 1, 1, 0, 0): ([], []),
        (0, 0, 2, 0, 0, 0): ([], []),
        (0, 1, 0, 0, 0, 0): ([((0, 0, 0, 0, 0, 0), Fr(-135, 8)), ((0, 0, 0, 0, 0, 1), Fr(-1, 12)), ((0, 0, 0, 0, 1, 0), Fr(7, 16)), ((0, 0, 0, 1, 0, 0), Fr(-7, 2)), ((0, 0, 0, 1, 1, 0), Fr(1, 8)), ((0, 0, 0, 2, 0, 0), Fr(-7, 16)), ((0, 0, 0, 3, 0, 0), Fr(-1, 24)), ((0, 0, 1, 0, 0, 0), Fr(247, 8)), ((0, 0, 1, 0, 1, 0), Fr(-3, 8)), ((0, 0, 1, 1, 0, 0), Fr(45, 8)), ((0, 0, 1, 2, 0, 0), Fr(3, 8)), ((0, 0, 2, 0, 0, 0), Fr(-63, 8)), ((0, 0, 2, 1, 0, 0), Fr(-3, 4)), ((0, 0, 3, 0, 0, 0), Fr(3, 4)), ((0, 1, 0, 0, 0, 0), Fr(115, 4)), ((0, 1, 0, 0, 1, 0), Fr(-3, 8)), ((0, 1, 0, 1, 0, 0), Fr(43, 8)), ((0, 1, 0, 2, 0, 0), Fr(3, 8)), ((0, 1, 1, 0, 0, 0), Fr(-101, 8)), ((0, 1, 1, 1, 0, 0), Fr(-5, 4)), ((0, 1, 2, 0, 0, 0), Fr(5, 4)), ((0, 2, 0, 0, 0, 0), Fr(-55, 8)), ((0, 2, 0, 1, 0, 0), Fr(-3, 4)), ((0, 2, 1, 0, 0, 0), Fr(5, 4)), ((0, 3, 0, 0, 0, 0), Fr(3, 4)), ((1, 0, 0, 0, 0, 0), Fr(217, 8)), ((1, 0, 0, 0, 1, 0), Fr(-3, 8)), ((1, 0, 0, 1, 0, 0), Fr(41, 8)), ((1, 0, 0, 2, 0, 0), Fr(3, 8)), ((1, 0, 1, 0, 0, 0), Fr(-95, 8)), ((1, 0, 1, 1, 0, 0), Fr(-5, 4)), ((1, 0, 2, 0, 0, 0), Fr(5, 4)), ((1, 1, 0, 0, 0, 0), Fr(-89, 8)), ((1, 1, 0, 1, 0, 0), Fr(-5, 4)), ((1, 1, 1, 0, 0, 0), Fr(7, 4)), ((1, 2, 0, 0, 0, 0), Fr(5, 4)), ((2, 0, 0, 0, 0, 0), Fr(-47, 8)), ((2, 0, 0, 1, 0, 0), Fr(-3, 4)), ((2, 0, 1, 0, 0, 0), Fr(5, 4)), ((2, 1, 0, 0, 0, 0), Fr(5, 4)), ((3, 0, 0, 0, 0, 0), Fr(3, 4))], [((0, 0, 0, 0, 0, 0), Fr(251, 8)), ((0, 0, 1, 0, 0, 0), Fr(11, 8)), ((0, 1, 0, 0, 0, 0), Fr(161, 8)), ((0, 1, 1, 0, 0, 0), Fr(1, 2)), ((0, 2, 0, 0, 0, 0), Fr(5, 2)), ((1, 0, 0, 0, 0, 0), Fr(-129, 8)), ((1, 0, 1, 0, 0, 0), Fr(-1, 4)), ((2, 0, 0, 0, 0, 0), Fr(-53, 8)), ((3, 0, 0, 0, 0, 0), Fr(-3, 4))]),
        (0, 1, 0, 0, 0, 1): ([], []),
        (0, 1, 0, 0, 1, 0): ([((0, 0, 0, 0, 0, 0), Fr(45, 2)), ((0, 0, 0, 0, 1, 0), Fr(-1, 4)), ((0, 0, 0, 1, 0, 0), Fr(4, 1)), ((0, 0, 0, 2, 0, 0), Fr(1, 4)), ((0, 0, 1, 0, 0, 0), Fr(-11, 2)), ((0, 0, 1, 1, 0, 0), Fr(-1, 2)), ((0, 0, 2, 0, 0, 0), Fr(1, 2)), ((0, 1, 0, 0, 0, 0), Fr(-5, 1)), ((0, 1, 0, 1, 0, 0), Fr(-1, 2)), ((0, 1, 1, 0, 0, 0), Fr(1, 2)), ((0, 2, 0, 0, 0, 0), Fr(1, 2)), ((1, 0, 0, 0, 0, 0), Fr(-9, 2)), ((1, 0, 0, 1, 0, 0), Fr(-1, 2)), ((1, 0, 1, 0, 0, 0), Fr(1, 2)), ((1, 1, 0, 0, 0, 0), Fr(1, 2)), ((2, 0, 0, 0, 0, 0), Fr(1, 2))], [((0, 0, 0, 0, 0, 0), Fr(13, 2)), ((0, 0, 1, 0, 0, 0), Fr(1, 2)), ((0, 1, 0, 0, 0, 0), Fr(5, 2)), ((1, 0, 0, 0, 0, 0), Fr(-5, 2)), ((2, 0, 0, 0, 0, 0), Fr(-1, 2))]),
        (0, 1, 0, 1, 0, 0): ([((0, 0, 0, 0, 0, 0), Fr(-45, 2)), ((0, 0, 0, 0, 1, 0), Fr(1, 4)), ((0, 0, 0, 1, 0, 0), Fr(-4, 1)), ((0, 0, 0, 2, 0, 0), Fr(-1, 4)), ((0, 0, 1, 0, 0, 0), Fr(11, 2)), ((0, 0, 1, 1, 0, 0), Fr(1, 2)), ((0, 0, 2, 0, 0, 0), Fr(-1, 2)), ((0, 1, 0, 0, 0, 0), Fr(5, 1)), ((0, 1, 0, 1, 0, 0), Fr(1, 2)), ((0, 1, 1, 0, 0, 0), Fr(-1, 2)), ((0, 2, 0, 0, 0, 0), Fr(-1, 2)), ((1, 0, 0, 0, 0, 0), Fr(9, 2)), ((1, 0, 0, 1, 0, 0), Fr(1, 2)), ((1, 0, 1, 0, 0, 0), Fr(-1, 2)), ((1, 1, 0, 0, 0, 0), Fr(-1, 2)), ((2, 0, 0, 0, 0, 0), Fr(-1, 2))], [((0, 0, 0, 0, 0, 0), Fr(-13, 2)), ((0, 0, 1, 0, 0, 0), Fr(-1, 2)), ((0, 1, 0, 0, 0, 0), Fr(-5, 2)), ((1, 0, 0, 0, 0, 0), Fr(5, 2)), ((2, 0, 0, 0, 0, 0), Fr(1, 2))]),
        (0, 1, 1, 0, 0, 0): ([], []),
        (0, 2, 0, 0, 0, 0): ([((0, 0, 0, 0, 0, 0), Fr(-45, 2)), ((0, 0, 0, 0, 1, 0), Fr(1, 4)), ((0, 0, 0, 1, 0, 0), Fr(-4, 1)), ((0, 0, 0, 2, 0, 0), Fr(-1, 4)), ((0, 0, 1, 0, 0, 0), Fr(11, 2)), ((0, 0, 1, 1, 0, 0), Fr(1, 2)), ((0, 0, 2, 0, 0, 0), Fr(-1, 2)), ((0, 1, 0, 0, 0, 0), Fr(5, 1)), ((0, 1, 0, 1, 0, 0), Fr(1, 2)), ((0, 1, 1, 0, 0, 0), Fr(-1, 2)), ((0, 2, 0, 0, 0, 0), Fr(-1, 2)), ((1, 0, 0, 0, 0, 0), Fr(9, 2)), ((1, 0, 0, 1, 0, 0), Fr(1, 2)), ((1, 0, 1, 0, 0, 0), Fr(-1, 2)), ((1, 1, 0, 0, 0, 0), Fr(-1, 2)), ((2, 0, 0, 0, 0, 0), Fr(-1, 2))], [((0, 0, 0, 0, 0, 0), Fr(-13, 2)), ((0, 0, 1, 0, 0, 0), Fr(-1, 2)), ((0, 1, 0, 0, 0, 0), Fr(-5, 2)), ((1, 0, 0, 0, 0, 0), Fr(5, 2)), ((2, 0, 0, 0, 0, 0), Fr(1, 2))]),
        (1, 0, 0, 0, 0, 0): ([((0, 0, 0, 0, 0, 0), Fr(-315, 8)), ((0, 0, 0, 0, 0, 1), Fr(-1, 12)), ((0, 0, 0, 0, 1, 0), Fr(11, 16)), ((0, 0, 0, 1, 0, 0), Fr(-15, 2)), ((0, 0, 0, 1, 1, 0), Fr(1, 8)), ((0, 0, 0, 2, 0, 0), Fr(-11, 16)), ((0, 0, 0, 3, 0, 0), Fr(-1, 24)), ((0, 0, 1, 0, 0, 0), Fr(291, 8)), ((0, 0, 1, 0, 1, 0), Fr(-3, 8)), ((0, 0, 1, 1, 0, 0), Fr(49, 8)), ((0, 0, 1, 2, 0, 0), Fr(3, 8)), ((0, 0, 2, 0, 0, 0), Fr(-67, 8)), ((0, 0, 2, 1, 0, 0), Fr(-3, 4)), ((0, 0, 3, 0, 0, 0), Fr(3, 4)), ((0, 1, 0, 0, 0, 0), Fr(135, 4)), ((0, 1, 0, 0, 1, 0), Fr(-3, 8)), ((0, 1, 0, 1, 0, 0), Fr(47, 8)), ((0, 1, 0, 2, 0, 0), Fr(3, 8)), ((0, 1, 1, 0, 0, 0), Fr(-105, 8)), ((0, 1, 1, 1, 0, 0), Fr(-5, 4)), ((0, 1, 2, 0, 0, 0), Fr(5, 4)), ((0, 2, 0, 0, 0, 0), Fr(-59, 8)), ((0, 2, 0, 1, 0, 0), Fr(-3, 4)), ((0, 2, 1, 0, 0, 0), Fr(5, 4)), ((0, 3, 0, 0, 0, 0), Fr(3, 4)), ((1, 0, 0, 0, 0, 0), Fr(253, 8)), ((1, 0, 0, 0, 1, 0), Fr(-3, 8)), ((1, 0, 0, 1, 0, 0), Fr(45, 8)), ((1, 0, 0, 2, 0, 0), Fr(3, 8)), ((1, 0, 1, 0, 0, 0), Fr(-99, 8)), ((1, 0, 1, 1, 0, 0), Fr(-5, 4)), ((1, 0, 2, 0, 0, 0), Fr(5, 4)), ((1, 1, 0, 0, 0, 0), Fr(-93, 8)), ((1, 1, 0, 1, 0, 0), Fr(-5, 4)), ((1, 1, 1, 0, 0, 0), Fr(7, 4)), ((1, 2, 0, 0, 0, 0), Fr(5, 4)), ((2, 0, 0, 0, 0, 0), Fr(-51, 8)), ((2, 0, 0, 1, 0, 0), Fr(-3, 4)), ((2, 0, 1, 0, 0, 0), Fr(5, 4)), ((2, 1, 0, 0, 0, 0), Fr(5, 4)), ((3, 0, 0, 0, 0, 0), Fr(3, 4))], [((0, 0, 0, 0, 0, 0), Fr(199, 8)), ((0, 0, 1, 0, 0, 0), Fr(7, 8)), ((0, 1, 0, 0, 0, 0), Fr(141, 8)), ((0, 1, 1, 0, 0, 0), Fr(1, 2)), ((0, 2, 0, 0, 0, 0), Fr(5, 2)), ((1, 0, 0, 0, 0, 0), Fr(-109, 8)), ((1, 0, 1, 0, 0, 0), Fr(-1, 4)), ((2, 0, 0, 0, 0, 0), Fr(-49, 8)), ((3, 0, 0, 0, 0, 0), Fr(-3, 4))]),
        (1, 0, 0, 0, 0, 1): ([], []),
        (1, 0, 0, 0, 1, 0): ([((0, 0, 0, 0, 0, 0), Fr(45, 2)), ((0, 0, 0, 0, 1, 0), Fr(-1, 4)), ((0, 0, 0, 1, 0, 0), Fr(4, 1)), ((0, 0, 0, 2, 0, 0), Fr(1, 4)), ((0, 0, 1, 0, 0, 0), Fr(-11, 2)), ((0, 0, 1, 1, 0, 0), Fr(-1, 2)), ((0, 0, 2, 0, 0, 0), Fr(1, 2)), ((0, 1, 0, 0, 0, 0), Fr(-5, 1)), ((0, 1, 0, 1, 0, 0), Fr(-1, 2)), ((0, 1, 1, 0, 0, 0), Fr(1, 2)), ((0, 2, 0, 0, 0, 0), Fr(1, 2)), ((1, 0, 0, 0, 0, 0), Fr(-9, 2)), ((1, 0, 0, 1, 0, 0), Fr(-1, 2)), ((1, 0, 1, 0, 0, 0), Fr(1, 2)), ((1, 1, 0, 0, 0, 0), Fr(1, 2)), ((2, 0, 0, 0, 0, 0), Fr(1, 2))], [((0, 0, 0, 0, 0, 0), Fr(13, 2)), ((0, 0, 1, 0, 0, 0), Fr(1, 2)), ((0, 1, 0, 0, 0, 0), Fr(5, 2)), ((1, 0, 0, 0, 0, 0), Fr(-5, 2)), ((2, 0, 0, 0, 0, 0), Fr(-1, 2))]),
        (1, 0, 0, 1, 0, 0): ([], []),
        (1, 0, 1, 0, 0, 0): ([((0, 0, 0, 0, 0, 0), Fr(-45, 2)), ((0, 0, 0, 0, 1, 0), Fr(1, 4)), ((0, 0, 0, 1, 0, 0), Fr(-4, 1)), ((0, 0, 0, 2, 0, 0), Fr(-1, 4)), ((0, 0, 1, 0, 0, 0), Fr(11, 2)), ((0, 0, 1, 1, 0, 0), Fr(1, 2)), ((0, 0, 2, 0, 0, 0), Fr(-1, 2)), ((0, 1, 0, 0, 0, 0), Fr(5, 1)), ((0, 1, 0, 1, 0, 0), Fr(1, 2)), ((0, 1, 1, 0, 0, 0), Fr(-1, 2)), ((0, 2, 0, 0, 0, 0), Fr(-1, 2)), ((1, 0, 0, 0, 0, 0), Fr(9, 2)), ((1, 0, 0, 1, 0, 0), Fr(1, 2)), ((1, 0, 1, 0, 0, 0), Fr(-1, 2)), ((1, 1, 0, 0, 0, 0), Fr(-1, 2)), ((2, 0, 0, 0, 0, 0), Fr(-1, 2))], [((0, 0, 0, 0, 0, 0), Fr(-13, 2)), ((0, 0, 1, 0, 0, 0), Fr(-1, 2)), ((0, 1, 0, 0, 0, 0), Fr(-5, 2)), ((1, 0, 0, 0, 0, 0), Fr(5, 2)), ((2, 0, 0, 0, 0, 0), Fr(1, 2))]),
        (1, 1, 0, 0, 0, 0): ([((0, 0, 0, 0, 0, 0), Fr(-45, 2)), ((0, 0, 0, 0, 1, 0), Fr(1, 4)), ((0, 0, 0, 1, 0, 0), Fr(-4, 1)), ((0, 0, 0, 2, 0, 0), Fr(-1, 4)), ((0, 0, 1, 0, 0, 0), Fr(11, 2)), ((0, 0, 1, 1, 0, 0), Fr(1, 2)), ((0, 0, 2, 0, 0, 0), Fr(-1, 2)), ((0, 1, 0, 0, 0, 0), Fr(5, 1)), ((0, 1, 0, 1, 0, 0), Fr(1, 2)), ((0, 1, 1, 0, 0, 0), Fr(-1, 2)), ((0, 2, 0, 0, 0, 0), Fr(-1, 2)), ((1, 0, 0, 0, 0, 0), Fr(9, 2)), ((1, 0, 0, 1, 0, 0), Fr(1, 2)), ((1, 0, 1, 0, 0, 0), Fr(-1, 2)), ((1, 1, 0, 0, 0, 0), Fr(-1, 2)), ((2, 0, 0, 0, 0, 0), Fr(-1, 2))], [((0, 0, 0, 0, 0, 0), Fr(-13, 2)), ((0, 0, 1, 0, 0, 0), Fr(-1, 2)), ((0, 1, 0, 0, 0, 0), Fr(-5, 2)), ((1, 0, 0, 0, 0, 0), Fr(5, 2)), ((2, 0, 0, 0, 0, 0), Fr(1, 2))]),
        (2, 0, 0, 0, 0, 0): ([((0, 0, 0, 0, 0, 0), Fr(-45, 2)), ((0, 0, 0, 0, 1, 0), Fr(1, 4)), ((0, 0, 0, 1, 0, 0), Fr(-4, 1)), ((0, 0, 0, 2, 0, 0), Fr(-1, 4)), ((0, 0, 1, 0, 0, 0), Fr(11, 2)), ((0, 0, 1, 1, 0, 0), Fr(1, 2)), ((0, 0, 2, 0, 0, 0), Fr(-1, 2)), ((0, 1, 0, 0, 0, 0), Fr(5, 1)), ((0, 1, 0, 1, 0, 0), Fr(1, 2)), ((0, 1, 1, 0, 0, 0), Fr(-1, 2)), ((0, 2, 0, 0, 0, 0), Fr(-1, 2)), ((1, 0, 0, 0, 0, 0), Fr(9, 2)), ((1, 0, 0, 1, 0, 0), Fr(1, 2)), ((1, 0, 1, 0, 0, 0), Fr(-1, 2)), ((1, 1, 0, 0, 0, 0), Fr(-1, 2)), ((2, 0, 0, 0, 0, 0), Fr(-1, 2))], [((0, 0, 0, 0, 0, 0), Fr(-13, 2)), ((0, 0, 1, 0, 0, 0), Fr(-1, 2)), ((0, 1, 0, 0, 0, 0), Fr(-5, 2)), ((1, 0, 0, 0, 0, 0), Fr(5, 2)), ((2, 0, 0, 0, 0, 0), Fr(1, 2))]),
    },
}
A0_TERMS = [((0, 0, 0, 0, 1, 0, 1, 0), Fr(1, 1)), ((0, 0, 0, 0, 1, 1, 0, 0), Fr(1, 1)), ((0, 0, 0, 1, 0, 0, 0, 0), Fr(3, 2)), ((0, 0, 0, 1, 0, 0, 1, 0), Fr(1, 2)), ((0, 0, 0, 1, 0, 1, 0, 0), Fr(1, 2)), ((0, 0, 0, 1, 1, 0, 0, 0), Fr(-3, 2)), ((0, 0, 0, 1, 1, 0, 1, 0), Fr(1, 2)), ((0, 0, 0, 1, 1, 1, 0, 0), Fr(1, 2)), ((0, 0, 1, 0, 1, 0, 1, 0), Fr(1, 2)), ((0, 0, 1, 0, 1, 1, 0, 0), Fr(1, 2)), ((0, 0, 1, 1, 0, 0, 0, 0), Fr(3, 2)), ((0, 0, 1, 1, 0, 0, 1, 0), Fr(1, 2)), ((0, 0, 1, 1, 0, 1, 0, 0), Fr(1, 2)), ((0, 0, 1, 1, 1, 0, 0, 0), Fr(-3, 2)), ((0, 1, 0, 0, 0, 0, 0, 0), Fr(-1, 2)), ((0, 1, 0, 0, 0, 0, 1, 0), Fr(1, 2)), ((0, 1, 0, 0, 0, 1, 0, 0), Fr(-1, 2)), ((0, 1, 0, 0, 1, 0, 0, 0), Fr(-1, 2)), ((0, 1, 0, 0, 1, 0, 1, 0), Fr(1, 2)), ((0, 1, 0, 0, 1, 1, 0, 0), Fr(-1, 2)), ((0, 1, 0, 1, 0, 0, 0, 0), Fr(-1, 2)), ((0, 1, 0, 1, 0, 0, 1, 0), Fr(1, 2)), ((0, 1, 0, 1, 0, 1, 0, 0), Fr(-1, 2)), ((0, 1, 1, 0, 0, 0, 0, 0), Fr(-1, 2)), ((0, 1, 1, 0, 0, 0, 1, 0), Fr(1, 2)), ((0, 1, 1, 0, 0, 1, 0, 0), Fr(-1, 2)), ((0, 2, 0, 0, 0, 0, 0, 0), Fr(-1, 2)), ((0, 2, 0, 0, 1, 0, 0, 0), Fr(-1, 2)), ((0, 2, 0, 1, 0, 0, 0, 0), Fr(-1, 2)), ((0, 2, 1, 0, 0, 0, 0, 0), Fr(-1, 2)), ((1, 0, 0, 0, 0, 0, 0, 0), Fr(-1, 2)), ((1, 0, 0, 0, 0, 0, 0, 1), Fr(1, 2)), ((1, 0, 0, 0, 0, 0, 1, 0), Fr(-1, 2)), ((1, 0, 0, 0, 1, 0, 0, 0), Fr(1, 1)), ((1, 0, 0, 0, 1, 0, 1, 0), Fr(-1, 4)), ((1, 0, 0, 0, 1, 1, 0, 0), Fr(1, 4)), ((1, 0, 0, 1, 0, 0, 0, 0), Fr(3, 4)), ((1, 0, 0, 1, 0, 0, 1, 0), Fr(-1, 4)), ((1, 0, 0, 1, 0, 1, 0, 0), Fr(1, 4)), ((1, 0, 0, 1, 1, 0, 0, 0), Fr(1, 4)), ((1, 0, 1, 0, 0, 0, 0, 0), Fr(1, 2)), ((1, 0, 1, 0, 0, 0, 1, 0), Fr(-1, 4)), ((1, 0, 1, 0, 0, 1, 0, 0), Fr(1, 4)), ((1, 0, 1, 0, 1, 0, 0, 0), Fr(1, 4)), ((1, 0, 1, 1, 0, 0, 0, 0), Fr(1, 4)), ((1, 1, 0, 0, 0, 0, 0, 0), Fr(3, 4)), ((1, 1, 0, 0, 0, 0, 1, 0), Fr(-1, 4)), ((1, 1, 0, 0, 0, 1, 0, 0), Fr(1, 4)), ((1, 1, 0, 0, 1, 0, 0, 0), Fr(1, 2)), ((1, 1, 0, 1, 0, 0, 0, 0), Fr(1, 2)), ((1, 1, 1, 0, 0, 0, 0, 0), Fr(1, 2)), ((1, 2, 0, 0, 0, 0, 0, 0), Fr(1, 4)), ((2, 0, 0, 0, 0, 0, 0, 0), Fr(-1, 1)), ((2, 0, 0, 0, 0, 0, 0, 1), Fr(1, 4)), ((2, 0, 0, 0, 0, 1, 0, 0), Fr(-1, 4)), ((2, 1, 0, 0, 0, 0, 0, 0), Fr(-1, 4)), ((3, 0, 0, 0, 0, 0, 0, 0), Fr(-1, 8))]


# ==== PART B
LAMS6 = [(6, 4, 2, 0), (9, 6, 3, 0), (5, 3, 2, 1), (7, 4, 1, 0), (8, 5, 3, 1), (6, 6, 3, 0)]
HELD_A = [(12, 8, 5, 1), (11, 7, 4, 0), (10, 9, 4, 2), (9, 8, 6, 1), (12, 10, 6, 3), (10, 6, 5, 2), (8, 7, 4, 3), (11, 9, 3, 0), (10, 7, 4, 1), (9, 5, 2, 1)]
GENERIC_LAMP = (Fr(13, 3), Fr(-7, 5), Fr(5, 7), Fr(11, 13))


def solve_exact(rows, rhs):
    """exact solve (free columns = 0): returns (x, rank) or None if inconsistent"""
    n = len(rows)
    m = len(rows[0])
    M = [list(r) + [b] for r, b in zip(rows, rhs)]
    piv = []
    r = 0
    for c in range(m):
        p = next((i for i in range(r, n) if M[i][c] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        inv = 1 / M[r][c]
        M[r] = [x * inv for x in M[r]]
        for i in range(n):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        piv.append(c)
        r += 1
        if r == n:
            break
    for i in range(r, n):
        if M[i][m] != 0:
            return None
    x = [Fr(0)] * m
    for k, c in enumerate(piv):
        x[c] = M[k][m]
    return x, r


def a_samples(lam, rng=None, nrand=0):
    """exact (kappa, a, a4a, b3, b4a, nu, lamp) of the kappa-block for every kappa1 > kappa2 > kappa3 of lam (lam' = lam, or a generic rational lam' when the gl(3) off-diagonal V vanishes);
    with nrand > 0 additionally nrand random (lam', two-point path) samples per kappa"""
    out_ = []
    for k in kappas_of(lam):
        try:
            res, nu = a_data(lam, lam, k)
            lamp = lam
        except ZeroDivisionError:
            res, nu = a_data(lam, GENERIC_LAMP, k)
            lamp = GENERIC_LAMP
        out_.append((k, res['V3'][0], res['V4a'][0], res['V3'][1], res['V4a'][1], nu, lamp))
        for _ in range(nrand):
            for _try in range(20):
                A = rng.randint(k[1], k[0] - 1)
                B = rng.randint(k[2] + 1, k[1])
                m11 = rng.randint(B, A)
                lp = tuple(Fr(rng.randint(-20, 20), rng.choice([1, 2, 3, 5])) for _ in range(4))
                try:
                    res, nu = a_data(lam, lp, k, A, B, m11)
                except ZeroDivisionError:
                    continue
                out_.append((k, res['V3'][0], res['V4a'][0], res['V3'][1], res['V4a'][1], nu, lp))
                break
    return out_


def stage_a(verbose=True):
    hdr('D1: the nesting scalars a(kappa; lam), a_4a(kappa; lam)')
    rng = random.Random(21)
    tot = Counter()
    data = {}
    t0 = time.perf_counter()
    for lam in LAMS6:
        data[lam] = a_samples(lam, rng, nrand=0)
        tab('  lam = %s: kappa -> (a, a4a)   [exact, G8b two-point paths]' % lstr(lam))
        for k, a, a4, b3, b4, nu, lp in data[lam]:
            tab('    %-10s a = %-14s a4a = %-14s   formula: a = %-14s a4a = %-14s' % (k, a, a4, a_of(k, lam), a4a_of(k, lam)))
    out('  data: %d kappa values of 6 lam (kappa1 > kappa2 > kappa3), computed in %.1f s' % (sum(len(v) for v in data.values()), time.perf_counter() - t0))
    # fit on [6,4,2,0] (21 kappa), verify on the other five lam
    L0 = LAMS6[0]
    rows, rhs = [], []
    for k, a, a4, b3, b4, nu, lp in data[L0]:
        R = R_list(k, L0)
        rows.append([sum(R)])
        rhs.append(a - 1)
    sol = solve_exact(rows, rhs)
    out('  (i) ansatz a - 1 = c * sum_i R_i  (1 free parameter) fitted on lam = %s (%d kappa): c = %s' % (lstr(L0), len(rows), sol[0][0] if sol else None))
    c_a = sol[0][0]
    held = bad = 0
    for lam in LAMS6[1:]:
        for k, a, a4, b3, b4, nu, lp in data[lam]:
            held += 1
            bad += (a != 1 + c_a * sum(R_list(k, lam)))
    out('      held-out: %d kappa on the other 5 lam: %d mismatches' % (held, bad))
    # a4a
    rows, rhs = [], []
    for k, a, a4, b3, b4, nu, lp in data[L0]:
        R = R_list(k, L0)
        y = ylist(k)
        d = Fr(sum(L0) - sum(k))
        rows.append([Fr(1), sum(R[i] * y[i] for i in range(3)), sum(k) * sum(R), sum(R)])
        rhs.append(a4 - d)
    sol = solve_exact(rows, rhs)
    out('  (ii) ansatz a4a - d = c0 + c1 sum R_i y_i + c2 |kappa| sum R_i + c3 sum R_i  (d = |lam| - |kappa|; 4 free parameters) fitted on %s (%d kappa): c = %s' % (lstr(L0), len(rows), [str(x) for x in sol[0]] if sol else None))
    c4 = sol[0]
    held = bad = 0
    for lam in LAMS6[1:]:
        for k, a, a4, b3, b4, nu, lp in data[lam]:
            R = R_list(k, lam)
            y = ylist(k)
            d = Fr(sum(lam) - sum(k))
            held += 1
            bad += (a4 != d + c4[0] + c4[1] * sum(R[i] * y[i] for i in range(3)) + c4[2] * sum(k) * sum(R) + c4[3] * sum(R))
    out('      held-out: %d kappa on the other 5 lam: %d mismatches' % (held, bad))
    # final formulas, all 150 + independence of lam', nu, and further lam
    n_all = bad_a = bad_a4 = 0
    for lam in LAMS6:
        for k, a, a4, b3, b4, nu, lp in data[lam]:
            n_all += 1
            bad_a += (a != a_of(k, lam))
            bad_a4 += (a4 != a4a_of(k, lam))
    out('  final formulas  a = 1 - sum_i R_i,   a4a = d + sum_i R_i (2 y_i - |kappa| + 4):  %d kappa of the 6 lam: %d / %d mismatches' % (n_all, bad_a, bad_a4))
    n_h = bh_a = bh_4 = 0
    t0 = time.perf_counter()
    for lam in HELD_A:
        for k, a, a4, b3, b4, nu, lp in a_samples(lam, rng, nrand=2):
            n_h += 1
            bh_a += (a != a_of(k, lam))
            bh_4 += (a4 != a4a_of(k, lam))
    out('  HELD-OUT lam (10 further lam, none used in the fits), each kappa with 3 samples (lam\' = lam / two random generic lam\', random weight space nu): %d samples, %d / %d mismatches [%.1f s]' % (n_h, bh_a, bh_4, time.perf_counter() - t0))
    # degenerate kappa: the corner values a_dega: on a 1x1 path only a V + b is visible; the pole-dropping convention is tested in stage b
    out('  in the labels: R_i = prod_{s=1..4} (kappa_i - lam_s - i + s) / prod_{j != i} (kappa_i - kappa_j - i + j)(kappa_i - kappa_j - i + j + 1);  a = 1 - R_1 - R_2 - R_3;  a4a = |lam| - |kappa| + sum_i R_i (2 kappa_i - 2 i - |kappa| + 4)')
    tab('  in the labels: R_i = prod_s (kappa_i - lam_s - i + s) / prod_{j != i}(kappa_i - kappa_j - i + j)(kappa_i - kappa_j - i + j + 1); a = 1 - sum R_i; a4a = |lam| - |kappa| + sum R_i (2 kappa_i - 2i - |kappa| + 4)')
    tab('  in l-coordinates y_i = kappa_i - i, Y_s = lam_s - s:   R_i = prod_{s=1..4}(y_i - Y_s) / prod_{j != i}(y_i - y_j)(y_i - y_j + 1);   a = 1 - sum_i R_i;   a4a = d + sum_i R_i (2 y_i - |kappa| + 4),   d = |lam| - |kappa|')
    tab('  lowering-first version (equal sums):  R-_i = prod_s (y_i - Y_s - 1) / prod_{j != i}(y_i - y_j)(y_i - y_j - 1);  a neighbour kappa +- e_i that is no dominant weight contributes R_i := 0 (pole dropped)')
    # identities: sum R+ = sum R-
    nid = nbad = 0
    for lam in LAMS6 + HELD_A:
        for k in kappas_of(lam):
            nid += 1
            nbad += (sum(R_list(k, lam, '+')) != sum(R_list(k, lam, '-')))
    out('  identity  sum_i R+_i = sum_i R-_i (raise-then-lower = lower-then-raise), all kappa (kappa1 > kappa2 > kappa3) of 16 lam: %d checked, %d violations' % (nid, nbad))
    # 3 random lam' independence on [9,6,3,0]
    n_i = bad_i = 0
    lam = (9, 6, 3, 0)
    for k, a, a4, b3, b4, nu, lp in a_samples(lam, rng, nrand=3):
        n_i += 1
        bad_i += (a != a_of(k, lam) or a4 != a4a_of(k, lam))
    out('  independence of (lam\', nu): [9,6,3,0], every kappa with 3 random (lam\', weight space) samples + the dominant one: %d samples, %d deviations from a(kappa), a4a(kappa)' % (n_i, bad_i))
    return dict(held=n_h, bad=bh_a + bh_4 + bad_a + bad_a4 + bad_i)


def stage_gl3chain(nspaces=600, seed=5):
    hdr('gl(3) chain with general lam\': closed forms (u, d, alpha) vs the exact reduced operator V = t S_{3,2} on GLNRep(kappa)')
    rng = random.Random(seed)
    tot = Counter()
    ks = [(k1, k2, k3) for k1 in range(0, 9) for k2 in range(0, k1 + 1) for k3 in range(0, k2 + 1)]
    t0 = time.perf_counter()
    rng.shuffle(ks)
    for kappa in ks:
        if tot['spaces'] >= nspaces:
            break
        R3 = gr.GLNRep.of(kappa)
        if R3.dim > 400:
            continue
        for nu, W in R3.wsp.items():
            lamp = tuple(Fr(rng.randint(-25, 25), rng.choice([1, 2, 3, 4, 7])) for _ in range(3))
            M = so.weight_space_mats(R3, {'V': so.reduced_op('V', lamp)}, W)['V']
            pos = {int(c): i for i, c in enumerate(W)}
            P, al, up, dn = gl3_chain(kappa, lamp, nu)
            # order of W in the GT index = ascending (kappa=fixed) row 2 = ascending m12
            idxs = [pos[R3.index[((nu[0],), (m, nu[0] + nu[1] - m), tuple(kappa))]] for m in P]
            tot['spaces'] += 1
            tot['paths>=3'] += (len(P) >= 3)
            ok = len(P) == len(W)
            for i, ri in enumerate(idxs):
                for j, cj in enumerate(idxs):
                    ex = M[ri][cj]
                    if i == j:
                        val = al[i]
                    elif i == j + 1:
                        val = up[j]
                    elif i == j - 1:
                        val = dn[j]
                    else:
                        val = Fr(0)
                    tot['entries'] += 1
                    tot['nonzero'] += (ex != 0)
                    if val != ex:
                        ok = False
                        tot['bad_entries'] += 1
            tot['bad_spaces'] += (not ok)
    tab("  closed forms (m1 = m12, m2 = s - m1 = m22, x = m1 - m2, xi = x(x+2), t = m11, s = m12 + m22, kappa = (k1,k2,k3), lamprime = q; chain coordinate m1; GT basis of GLNRep(kappa)):")
    tab('    up    M[m1+1, m1] = -(k1-m1)(k2-m1-1)(k3-m1-2)(t+1-m2)(t-m1+q1-q2) / ((x+1)(x+2))')
    tab('    down  M[m1-1, m1] = -(k2-m2)(t-m1)(k1-m2+1)(k3-m2-1)(t+1-m2+q1-q2) / (x(x+1))')
    tab("    alpha M[m1, m1]   = A0 + g xi + Kp/xi,  g = 3s/16 - (k1+k2+k3)/8 + (q1+q2)/8 - q3/4 + 3/8,  Kp = (s-2k2)(s-2t)(s-2k1-2)(s-2k3+2)(2-s+2t+2q1-2q2)/16 (Kp/xi := 0 at x = 0)")
    tab('    A0 = sum of %d terms c * s^a t^b k1^c k2^d k3^e q1^f q2^g q3^h (exponent tuples (a,b,c,d,e,f,g,h)):' % len(A0_TERMS))
    for e, co in A0_TERMS:
        tab('       %s  %s' % (e, co))
    out('  random gl(3) irreps (dim <= 400) with generic rational lam\' (3 components, denominators 1,2,3,4,7): %d weight spaces (%d with a path of >= 3 points), %d matrix entries (%d nonzero), %d bad spaces, %d mismatching entries  [%.1f s]'
        % (tot['spaces'], tot['paths>=3'], tot['entries'], tot['nonzero'], tot['bad_spaces'], tot['bad_entries'], time.perf_counter() - t0))
    return tot


# ======================================================================================== D3: word groups
def word_group(w):
    """alpha: all indices in {1,2,3} or the Cartan E44;  beta: exactly one E_{i4} and one E_{4j} (i, j <= 3);  gamma: the rest"""
    if all((g[0] <= 3 and g[1] <= 3) or g == (4, 4) for g in w):
        return 'alpha'
    n4 = sum(1 for g in w if g[0] != g[1] and g[1] == 4)
    m4 = sum(1 for g in w if g[0] != g[1] and g[0] == 4)
    return 'beta' if (n4 == 1 and m4 == 1) else 'gamma'


def group_ops(name, lamp):
    """{group: ReducedOp} of the word groups of the reduced operator `name` at lam' = lamp"""
    univ = so.universal_reduction(name)
    parts = {'alpha': {}, 'beta': {}, 'gamma': {}}
    for w, p in univ.items():
        v = so.eval_poly(p, lamp)
        if v != 0:
            parts[word_group(w)][w] = v
    return {g: so.ReducedOp(d, name + ':' + g, lamp) for g, d in parts.items() if d}


def fit_v4b_block(J, V, D):
    n = len(J)
    comm = [[(D[i] - D[j]) * V[i][j] for j in range(n)] for i in range(n)]
    one = [[Fr(int(i == j)) for j in range(n)] for i in range(n)]
    Dm = [[D[i] if i == j else Fr(0) for j in range(n)] for i in range(n)]
    D2 = [[D[i] ** 2 if i == j else Fr(0) for j in range(n)] for i in range(n)]
    return ng.span_fit(J, [V, comm, one, Dm, D2])


def words_product(lam, lamp, names=('V3', 'V4a', 'V4b'), min_path=2):
    """per (operator, group): number of path blocks fitted exactly as a V^(3) + b (V3, V4a: and compared to the closed-form a of the group), per-kappa a values"""
    R = gr.GLNRep.of(lam)
    cnt = Counter()
    seen = {}
    gops = {n: group_ops(n, lamp) for n in names}
    for mu, m, w in ng.sectors(lam, lamp):
        if w < 2:
            continue
        nu = ng.sub_mu(mu, lamp)
        W = R.wsp[nu]
        P = ng.paths_of(R, W)
        pos = {int(c): a for a, c in enumerate(W)}
        for n in names:
            M = {g: so.weight_space_mats(R, {g: op}, W)[g] for g, op in gops[n].items()}
            for kappa, path in P.items():
                if len(path) < min_path:
                    continue
                V = ng.gl3_block_for_path(lam, lamp, mu, kappa)
                D = [Fr((R.pats[c][1][0] - R.pats[c][1][1]) * (R.pats[c][1][0] - R.pats[c][1][1] + 2)) for c in path]
                for g, Mg in M.items():
                    J = ng.block(Mg, pos, path, path)
                    if n in ('V3', 'V4a'):
                        f = ng.fit_ansatz(J, V, None)
                        cnt[(n, g, 'pairs')] += 1
                        if f['ok1']:
                            cnt[(n, g, 'fit')] += 1
                            seen.setdefault((n, g, kappa), set()).add(f['a'])
                    else:
                        f = fit_v4b_block(J, V, D)
                        cnt[(n, g, 'pairs')] += 1
                        if f is not None:
                            cnt[(n, g, 'fit')] += 1
                            seen.setdefault((n, g, kappa), set()).add((f[0], f[1]))
    return cnt, seen


def stage_words(products=(((4, 2, 1, 0),) * 2, ((6, 4, 2, 0),) * 2)):
    hdr('D3: the 39 / 78 / 85 universal reduced words of V3 / V4a / V4b split into the groups alpha (gl(3) + E44), beta (one E_{i4}, one E_{4j}), gamma (rest)')
    for n in ('V3', 'V4a', 'V4b'):
        univ = so.universal_reduction(n)
        c = Counter(word_group(w) for w in univ)
        out('  %-4s %d universal words (max length %d): alpha %d, beta %d, gamma %d' % (n, len(univ), max(len(w) for w in univ), c['alpha'], c['beta'], c['gamma']))
    tot = Counter()
    for lam, lamp in products:
        t0 = time.perf_counter()
        cnt, seen = words_product(lam, lamp)
        bad = Counter()
        n_k = Counter()
        for (n, g, kappa), vals in sorted(seen.items()):
            d = Fr(sum(lam) - sum(kappa))
            if n == 'V3':
                ref = {'alpha': Fr(1), 'beta': a_of(kappa, lam) - 1}[g]
            elif n == 'V4a':
                ref = {'alpha': d, 'beta': a4a_of(kappa, lam) - d}[g]
            else:
                ref = None
            if ref is not None:
                n_k[(n, g)] += 1
                bad[(n, g)] += (vals != {ref})
        out('  %s^2: path blocks (|P| >= 2) fitted exactly by  a_group V^(3) + b_group  (V4b: a V + e [xi,V] + b + c xi + g xi^2)  [%.1f s]' % (lstr(lam), time.perf_counter() - t0))
        for n in ('V3', 'V4a', 'V4b'):
            for g in ('alpha', 'beta', 'gamma'):
                if cnt[(n, g, 'pairs')]:
                    extra = ''
                    if n in ('V3', 'V4a') and n_k[(n, g)]:
                        extra = '   a_group = %s on all %d kappa: %d deviations' % ({('V3', 'alpha'): '1', ('V3', 'beta'): 'a - 1', ('V4a', 'alpha'): 'd = |lam|-|kappa|', ('V4a', 'beta'): 'a4a - d'}[(n, g)], n_k[(n, g)], bad[(n, g)])
                    out('      %-4s %-5s: %d / %d blocks fit%s' % (n, g, cnt[(n, g, 'fit')], cnt[(n, g, 'pairs')], extra))
                    tot[(n, g, 'fit')] += cnt[(n, g, 'fit')]
                    tot[(n, g, 'pairs')] += cnt[(n, g, 'pairs')]
                    tot[(n, g, 'bad')] += bad[(n, g)]
        # V4b group values: e of the alpha and beta parts
        e_al = {v for (n, g, k), vals in seen.items() if n == 'V4b' and g == 'alpha' for v in vals}
        e_be = {v[1] for (n, g, k), vals in seen.items() if n == 'V4b' and g == 'beta' for v in vals}
        e_alpha = {v[1] for (n, g, k), vals in seen.items() if n == 'V4b' and g == 'alpha' for v in vals}
        out('      V4b: e of the alpha part takes the values %s ; e of the beta part takes %d distinct values' % (sorted(str(x) for x in e_alpha)[:6], len(e_be)))
    return tot


# ======================================================================================== D2: b3, b4a held-out, degenerate kappa; the V4b within-path parameters (structure only)
def b_heldout(lams, nrand=3, seed=3):
    """b3, b4a from the closed-form tables vs the exact shifts J_ii - a V_ii (a from the data) on random (kappa, weight space, lam') samples of lam not used in the fits"""
    rng = random.Random(seed)
    tot = Counter()
    for lam in lams:
        for k, a, a4, b3, b4, nu, lp in a_samples(lam, rng, nrand=nrand):
            tot['samples'] += 1
            tot['bad3'] += (b_of('V3', k, lam, nu, lp) != b3)
            if 'V4a' in B_TABLES:
                tot['bad4a'] += (b_of('V4a', k, lam, nu, lp) != b4)
    return tot


def degenerate_check(lams, nper=6, seed=4):
    """kappa with kappa_i = kappa_{i+1} (1 x 1 paths, only a V + b is visible): the full diagonal entry J_ii of G3's operator vs a(kappa) alpha_ii + b with the pole-dropping convention"""
    rng = random.Random(seed)
    tot = Counter()
    for lam in lams:
        R = gr.GLNRep.of(lam)
        for k in kappas_of(lam, strict=False):
            if k[0] > k[1] > k[2]:
                continue
            for _ in range(nper):
                A = rng.randint(k[1], k[0])
                B = rng.randint(k[2], k[1])
                m11 = rng.randint(B, A)
                p = ((m11,), (A, B), tuple(k), tuple(lam))
                ip = R.index[p]
                nu = R.wt[ip]
                lp = tuple(Fr(rng.randint(-20, 20), rng.choice([1, 2, 3, 5])) for _ in range(4))
                al, up, dn = gl3_entries(k, lp[:3], m11, A + B, A)
                for n in _CLOSED_OPS:
                    if n not in B_TABLES:
                        continue
                    J = so.reduced_op(n, lp).apply(R, {ip: Fr(1)}).get(ip, Fr(0))
                    a = a_of(k, lam) if n == 'V3' else a4a_of(k, lam)
                    pred = a * al + b_of(n, k, lam, nu, lp)
                    tot['tests'] += 1
                    tot['bad'] += (pred != J)
                    tot['dropped'] += any(r == 0 for r in R_list(k, lam)) 
    return tot


def stage_b():
    hdr('D2: the shifts b3(kappa; nu, lam\'), b4a(kappa; nu, lam\') as polynomials in (nu1, nu2, lam\') (degree 2) with coefficients  P(kappa, p1, p2, p3) + sum_i R_i S(y_i, e1(y), e2(y), p1, p2, p3),  p_k = sum_s Y_s^k')
    tab('  B_TABLES (exponent tuples: P over (k1,k2,k3,p1,p2,p3), S over (y,e1,e2,p1,p2,p3)); monomials over (nu1, nu2, q1, q2, q3, q4)')
    for n, tb in B_TABLES.items():
        tab('  %s:' % n)
        nm = ['k1', 'k2', 'k3', 'p1', 'p2', 'p3']
        nr = ['y', 'e1', 'e2', 'p1', 'p2', 'p3']
        for m, (pt, rt) in tb.items():
            def mono_s(e, names):
                return '*'.join('%s^%d' % (x, k) if k > 1 else x for x, k in zip(names, e) if k) or '1'
            tab('    c[%s] = %s%s' % (m, ' + '.join('(%s)%s' % (co, mono_s(e, nm)) for e, co in pt) or '0', ''.join(' + R[(%s)%s]' % (co, mono_s(e, nr)) for e, co in rt)))
    for n in B_TABLES:
        nterm = sum(len(p) + len(r) for p, r in B_TABLES[n].values())
        out('  %s: 28 coefficient functions, %d terms in total (tables in closed_forms_gl4_tables.txt)' % (n, nterm))
    t0 = time.perf_counter()
    tot = b_heldout(HELD_A, nrand=3)
    out('  HELD-OUT lam (10 lam never used for the fits): %d samples (kappa x lam\' x weight space): b3 mismatches %d, b4a mismatches %d  [%.1f s]' % (tot['samples'], tot['bad3'], tot['bad4a'], time.perf_counter() - t0))
    dg = degenerate_check([(9, 6, 3, 0), (8, 5, 3, 1), (7, 4, 1, 0), (10, 7, 4, 1), (6, 6, 3, 0), (7, 7, 3, 3)])
    out('  degenerate kappa (kappa_i = kappa_{i+1}: 1x1 paths, no neighbour kappa + e_i): full diagonal J_ii = a alpha_ii + b with the pole-dropping convention R_i := 0: %d tests (%d with a dropped R_i), %d mismatches' % (dg['tests'], dg['dropped'], dg['bad']))
    return dict(held=tot, degenerate=dg)


# ======================================================================================== D4: the recursion from closed forms
class ClosedSector:
    """a G8a-style sector whose operator matrices come from the closed forms (V3, V4a) and V4b := 0 (the closed forms of V4b are not available: the recursion then uses the two equations
    (V3 - v3) psi = 0, (V4a - v4a) psi = 0 only).  Attributes as lattice_recursion_gl4.Sector."""

    def __init__(self, lam, lamp, mu):
        self.lam, self.lamp, self.mu = tuple(lam), tuple(lamp), tuple(mu)
        self.nu = tuple(mu[i] - lamp[i] for i in range(4))
        t0 = time.perf_counter()
        self.pats, M = closed_form_mats(self.lam, self.lamp, self.mu)
        self.t_mats = time.perf_counter() - t0
        self.w = len(self.pats)
        self.ok = self.w > 0
        if not self.ok:
            return
        zero = [[Fr(0)] * self.w for _ in range(self.w)]
        self.M = [M['V3'], M['V4a'], zero]
        self.xs = [sg.coords(p) for p in self.pats]
        self.k1 = lamp[0] - lamp[1] + 1
        self.cut = mu[0] - lamp[1]
        self.in1 = [i for i, x in enumerate(self.xs) if x[0] <= self.cut]
        self._spaces = {}

    def space(self, name):
        if name not in self._spaces:
            self._spaces[name] = lr.Space(self, name, range(self.w) if name == 'W' else self.in1)
        return self._spaces[name]


def best_gradings_closed(sec, space='W1', first_ok=False):
    """[(n0, band, |g|_1, g)] of the fully determined gradings of the space for the TWO closed-form operators (the V4b rows are zero), sorted"""
    sp = sec.space(space)
    res = []
    for g in lr.GRADINGS:
        a = lr.analyze(sp, g)
        if a['ok']:
            res.append((a['n0'], a['band'], sum(abs(x) for x in g), g))
            if first_ok and a['n0'] == 1:
                break
    res.sort()
    return res


def _null_space(C, n0, tol_exp=30):
    """null space of the (rows x n0) mp matrix C by Gaussian elimination with full pivoting (tolerance 10^-tol_exp relative): list of null vectors"""
    A = [list(r) for r in C]
    mx = max([abs(x) for r in A for x in r] + [mp.mpf(0)])
    if mx == 0:
        return [[mp.mpf(int(i == j)) for j in range(n0)] for i in range(n0)]
    tol = mx * mp.mpf(10) ** (-tol_exp)
    piv = []
    r = 0
    cols = list(range(n0))
    for _ in range(n0):
        best = None
        for i in range(r, len(A)):
            for j in cols:
                if j in piv:
                    continue
                if best is None or abs(A[i][j]) > abs(A[best[0]][best[1]]):
                    best = (i, j)
        if best is None or abs(A[best[0]][best[1]]) <= tol:
            break
        i, j = best
        A[r], A[i] = A[i], A[r]
        pv = A[r][j]
        A[r] = [x / pv for x in A[r]]
        for k in range(len(A)):
            if k != r and A[k][j] != 0:
                f = A[k][j]
                A[k] = [x - f * y for x, y in zip(A[k], A[r])]
        piv.append(j)
        r += 1
    free = [j for j in range(n0) if j not in piv]
    basis = []
    for f in free:
        x = [mp.mpf(0)] * n0
        x[f] = mp.mpf(1)
        for k, j in enumerate(piv):
            x[j] = -A[k][f]
        basis.append(x)
    return basis


def slice_from_closed(sec, trip, g, dps=50, plan=None):
    """one slice psi (dict pattern -> mpf, normalised: first non-negligible component in the GT order = 1) from the closed-form matrices by the forward lattice recursion of grading g on W1;
    trip = (v3, v4a); n0 > 1: psi0 = the null vector of the characteristic equations C(v) psi0 = 0"""
    old = mp.mp.dps
    mp.mp.dps = dps + 10
    try:
        sp = sec.space('W1')
        pl = plan if plan is not None else lr.build_plan(sp, g)
        assert pl.ok
        v = [mp.mpf(trip[0].numerator) / trip[0].denominator if isinstance(trip[0], Fr) else mp.mpf(trip[0]), mp.mpf(trip[1].numerator) / trip[1].denominator if isinstance(trip[1], Fr) else mp.mpf(trip[1]), mp.mpf(0)]
        key = ('mp', dps + 10)
        if pl.n0 == 1:
            psi, unused = lr.run_plan(pl, v, [[mp.mpf(1)]], lr.mpq, mp.mpf(0), key)
            res = max([abs(u[0]) for u in unused] + [mp.mpf(0)])
            vec = [mp.mpf(0)] * sec.w
            for p_, val in psi.items():
                vec[sp.pos[p_]] = val[0]
        else:
            n0 = pl.n0
            cols = [[mp.mpf(int(i == j)) for j in range(n0)] for i in range(n0)]
            psiT, unused = lr.run_plan(pl, v, cols, lr.mpq, mp.mpf(0), key)
            C = [[u[j] for j in range(n0)] for u in unused]
            basis = _null_space(C, n0)
            ambiguous = len(basis) != 1
            x0 = basis[0] if basis else [mp.mpf(1)] + [mp.mpf(0)] * (n0 - 1)
            res = max(abs(sum(row[j] * x0[j] for j in range(n0))) for row in C)
            vec = [mp.mpf(0)] * sec.w
            for p_, val in psiT.items():
                vec[sp.pos[p_]] = sum(val[j] * x0[j] for j in range(n0))
        mx = max(abs(x) for x in vec)
        tolz = mp.mpf(10) ** (-(dps - 20))
        first = next(a for a in range(sec.w) if abs(vec[a]) > tolz * mx)
        sc = vec[first]
        nz = {sec.pats[a]: vec[a] / sc for a in range(sec.w) if abs(vec[a]) > tolz * mx}
        if pl.n0 > 1 and ambiguous:
            res = mp.mpf('inf')          # the (V3, V4a) characteristic equations leave a null space of dimension != 1: the copy is not determined by two operators
        return nz, res
    finally:
        mp.mp.dps = old


def physical_triples(lam, lamp, mu, dps=50):
    """the joint eigenvalue triples (v3, v4a, v4b) of the physical copies of the sector from G3's m x m PRV matrices (this part still uses representation matrices of size w; the closed-form
    matrices of V3, V4a cannot see which eigenvalues of W1 are physical), ascending lexicographic as G3's copies"""
    prv = se.prv_space(tuple(lam), tuple(lamp), tuple(mu))
    _, mats = so.triple_mats(tuple(lam), tuple(lamp), tuple(mu), prv)
    jd = so.joint_diag([mats[n] for n in so.TRIPLE], dps=dps)
    return jd['trip']


def slices_closed(lam, lamp, mu, dps=50, trips=None, timing=None):
    """all copies of the sector mu in V_lam (x) V_lam' from the closed-form matrices (V3, V4a from labels alone) and the G8a recursion: list of (triple, slice dict pattern -> mpf).
    Eigenvalues: G3's PRV matrices (see physical_triples).  timing (dict) collects the breakdown."""
    timing = timing if timing is not None else {}
    t0 = time.perf_counter()
    sec = ClosedSector(lam, lamp, mu)
    timing['mats'] = timing.get('mats', 0) + sec.t_mats
    t1 = time.perf_counter()
    if trips is None:
        trips = sorted(physical_triples(lam, lamp, mu, dps), key=_g3_key(dps))
    timing['eigen'] = timing.get('eigen', 0) + time.perf_counter() - t1
    t1 = time.perf_counter()
    bg = best_gradings_closed(sec, 'W1', first_ok=True)
    if not bg:
        raise ValueError('no fully determined grading with (V3, V4a)')
    g = bg[0][3]
    pl = lr.build_plan(sec.space('W1'), g)
    timing['plan'] = timing.get('plan', 0) + time.perf_counter() - t1
    t1 = time.perf_counter()
    outl = []
    for trip in trips:
        nz, res = slice_from_closed(sec, (trip[0], trip[1]), g, dps, plan=pl)
        outl.append((trip, nz, res))
    timing['recursion'] = timing.get('recursion', 0) + time.perf_counter() - t1
    timing['n0'] = pl.n0
    timing['grading'] = g
    return outl


def compare_mats(lam, lamp, mu, ops=_CLOSED_OPS):
    """(w, closed matrices == G3's exactly for each op, t_closed, t_g3)"""
    t0 = time.perf_counter()
    pats, M = closed_form_mats(lam, lamp, mu, ops)
    t1 = time.perf_counter()
    pg, G = g3_mats(lam, lamp, mu, ops)
    t2 = time.perf_counter()
    assert pats == pg, 'weight-space patterns/order differ'
    return len(pats), {n: M[n] == G[n] for n in ops}, t1 - t0, t2 - t1


def stage_mats(heldout=True, seed=17):
    hdr('D4a: closed_form_mats (V3, V4a from labels alone) == G3 weight_space_mats, EXACTLY, on every weight space of the listed products')
    tot = Counter()
    for lam, lamp in ng.P_SMALL + ng.P_BIG:
        c = Counter()
        tc = tg = 0.0
        for mu, m, w in ng.sectors(lam, lamp):
            ww, eq, t1, t2 = compare_mats(lam, lamp, mu)
            c['spaces'] += 1
            c['entries'] += 2 * ww * ww
            c['w'] += ww
            c['wmax'] = max(c['wmax'], ww)
            for n, ok in eq.items():
                c['bad_' + n] += (not ok)
            tc += t1
            tg += t2
        out('  %-26s %3d weight spaces (w <= %2d, %5d matrix entries x 2 ops): V3 mismatching spaces %d, V4a %d   [closed %.2fs, G3 %.2fs]' % (lstr(lam) + ('x' + lstr(lamp) if lamp != lam else '^2'), c['spaces'], c['wmax'], c['entries'], c['bad_V3'], c['bad_V4a'], tc, tg))
        for k, v in c.items():
            if k != 'wmax':
                tot[k] += v
    out('  TOTAL products P_small + P_big: %d weight spaces, %d matrix entries; mismatching spaces V3 %d, V4a %d' % (tot['spaces'], tot['entries'], tot['bad_V3'], tot['bad_V4a']))
    sd = scale_demo()
    out('  beyond G3: lam = (30,20,10,0) (dim V_lam = %d, GLNRep infeasible): closed_form_mats on %d random weight spaces (w = %s) built in %s s; exact check [V3, V4a] = 0 on W: %s' % (
        gr.weyl_dim((30, 20, 10, 0)), len(sd), ','.join(str(x[1]) for x in sd), ','.join('%.2f' % x[2] for x in sd), 'all True' if all(x[3] for x in sd) else 'FAILED'))
    comm_n = comm_bad = 0
    for lam, lamp in ng.P_SMALL + ng.P_BIG:
        for mu, m, w in ng.sectors(lam, lamp):
            if w < 3:
                continue
            pats, M = closed_form_mats(lam, lamp, mu)
            n = len(pats)
            AB = [[sum(M['V3'][i][k] * M['V4a'][k][j] for k in range(n) if M['V3'][i][k] != 0) for j in range(n)] for i in range(n)]
            BA = [[sum(M['V4a'][i][k] * M['V3'][k][j] for k in range(n) if M['V4a'][i][k] != 0) for j in range(n)] for i in range(n)]
            comm_n += 1
            comm_bad += (AB != BA)
    out('  [V3, V4a] = 0 on the WHOLE weight space (not only on M_mu) for the closed-form matrices: %d weight spaces with w >= 3 of the 12 products, %d violations' % (comm_n, comm_bad))
    if heldout:
        rng = random.Random(seed)
        lam = (9, 6, 3, 0)
        R = gr.GLNRep.of(lam)
        nus = [nu for nu, W in R.wsp.items() if 2 <= len(W) <= 14]
        rng.shuffle(nus)
        c = Counter()
        t0 = time.perf_counter()
        for nu in nus[:20]:
            lamp = tuple(Fr(rng.randint(-25, 25), rng.choice([1, 2, 3, 5, 7])) for _ in range(4))
            mu = tuple(nu[i] + lamp[i] for i in range(4))
            ww, eq, t1, t2 = compare_mats(lam, lamp, mu)
            c['spaces'] += 1
            c['w'] += ww
            c['entries'] += 2 * ww * ww
            for n, ok in eq.items():
                c['bad_' + n] += (not ok)
        out('  [9,6,3,0]^2 with GENERIC non-dominant rational lam\' (random weight spaces, 2 <= w <= 14): %d weight spaces (total w %d, %d entries x 2 ops): V3 mismatches %d, V4a %d  [%.1f s]' % (c['spaces'], c['w'], c['entries'], c['bad_V3'], c['bad_V4a'], time.perf_counter() - t0))
        tot['held'] = c['spaces']
        tot['held_bad'] = c['bad_V3'] + c['bad_V4a']
    return tot


def scale_demo(lam=(30, 20, 10, 0), nspaces=4, seed=2):
    """closed_form_mats far beyond the reach of GLNRep(lam) (dim V_lam ~ 1e6): the exact matrices of V3, V4a on random weight spaces (20 <= w <= 60); independent exact check: [V3, V4a] = 0 on W"""
    rng = random.Random(seed)
    tot = sum(lam)
    done = 0
    t00 = time.perf_counter()
    res = []
    while done < nspaces:
        nu = sorted([rng.randint(8, 24) for _ in range(3)], reverse=True)
        nu.append(tot - sum(nu))
        if nu[3] < 0:
            continue
        nu = tuple(nu)
        pats = weight_patterns(lam, nu)
        if not (20 <= len(pats) <= 60):
            continue
        t0 = time.perf_counter()
        pats, M = closed_form_mats(lam, lam, tuple(nu[i] + lam[i] for i in range(4)))
        dt = time.perf_counter() - t0
        n = len(pats)

        def mm(A, B):
            return [[sum(A[i][k] * B[k][j] for k in range(n) if A[i][k] != 0) for j in range(n)] for i in range(n)]
        comm0 = mm(M['V3'], M['V4a']) == mm(M['V4a'], M['V3'])
        res.append((nu, n, dt, comm0))
        done += 1
    return res


def _g3_key(dps):
    return so.tol_cmp(mp.mpf(10) ** (-(dps // 2)))


def slices_product(lam, lamp, dps=50, mus=None, verbose=True):
    """validation + timing of slices_closed vs G3's copies and G8a's slice_from_triple on the sectors (m > 0) of the product"""
    old = mp.mp.dps
    mp.mp.dps = dps + 10
    R = gr.GLNRep.of(lam)
    tot = Counter()
    T = defaultdict(float)
    worst = mp.mpf(0)
    try:
        secs = [(mu, m) for mu, m, w in ng.sectors(lam, lamp) if m > 0 and (mus is None or mu in mus)]
        for mu, m in secs:
            t0 = time.perf_counter()
            cps = so.copies(lam, lamp, mu, dps=dps)
            tg3 = time.perf_counter() - t0
            T['copies'] += tg3
            timing = {}
            t0 = time.perf_counter()
            trips = physical_triples(lam, lamp, mu, dps)
            key = _g3_key(dps)
            trips = sorted(trips, key=key)
            te = time.perf_counter() - t0
            T['eigen'] += te
            res = slices_closed(lam, lamp, mu, dps, trips=trips, timing=timing)
            for k in ('mats', 'plan', 'recursion'):
                T[k] += timing.get(k, 0)
            T['closed_total'] += te + timing['mats'] + timing['plan'] + timing['recursion']
            tot['sectors'] += 1
            tot['copies'] += len(cps)
            tot['n0>1'] += (timing.get('n0', 1) > 1)
            assert len(res) == len(cps) == m
            for cp, (trip, nz, resid) in zip(cps, res):
                ref = {R.pats[i]: v for i, v in cp['slice'].items()}
                keys = set(ref) | set(nz)
                dv = max([abs(ref.get(k, 0) - nz.get(k, 0)) for k in keys] + [mp.mpf(0)])
                if resid == mp.inf:
                    tot['ambiguous'] += 1
                    continue
                worst = max(worst, dv)
                tot['bad'] += (dv > mp.mpf(10) ** -30)
        tot['worst'] = worst
        return tot, T
    finally:
        mp.mp.dps = old


def g8a_timing(lam, lamp, mus, dps=50):
    """G8a's slice_from_triple on the sectors mus whose best grading (with the triple) has n0 = 1, with G3's matrices: (time of building the sectors (weight_space_mats + grading search),
    time of the recursions, number of copies, number of sectors skipped for n0 > 1)"""
    tm = tr = 0.0
    n = skipped = 0
    for mu in mus:
        cps = so.copies(lam, lamp, mu, dps=dps)
        lr._SEC.pop((tuple(lam), tuple(lamp), tuple(mu)), None)
        t0 = time.perf_counter()
        sec = lr.sector(lam, lamp, mu)
        bg = lr.best_gradings(sec, 'W1')
        tm += time.perf_counter() - t0
        if not bg or bg[0][0] != 1:
            skipped += 1
            continue
        t0 = time.perf_counter()
        for cp in cps:
            lr.slice_from_triple(lam, lamp, mu, cp['ev'], dps=dps)
            n += 1
        tr += time.perf_counter() - t0
    return tm, tr, n, skipped


def two_op_determinacy(products=None):
    """fraction of weight spaces (w >= 2) of the products on which the forward lattice recursion is determined with (V3, V4a) alone (closed forms, best grading of the 124) vs with the triple (G8a: 222/222)"""
    products = products or (ng.P_SMALL + ng.P_BIG)
    tot = Counter()
    for lam, lamp in products:
        for mu, m, w in ng.sectors(lam, lamp):
            if w < 2:
                continue
            sec = ClosedSector(lam, lamp, mu)
            bg2 = best_gradings_closed(sec, 'W1')
            s3 = lr.sector(lam, lamp, mu)
            bg3 = lr.best_gradings(s3, 'W1')
            tot['spaces'] += 1
            tot['det2'] += bool(bg2)
            tot['det3'] += bool(bg3)
            tot['n0one2'] += bool(bg2 and bg2[0][0] == 1)
            tot['n0one3'] += bool(bg3 and bg3[0][0] == 1)
            tot['m>0'] += (m > 0)
            tot['m>0_det2'] += (m > 0 and bool(bg2))
    return tot


def stage_slices(products=(((6, 4, 2, 0),) * 2,), dps=50, extra=True):
    hdr('D4b: slices_closed (closed-form V3, V4a + forward lattice recursion) vs G3 copies (50 digits); timings vs G3 copies and G8a slice_from_triple')
    res = {}
    for lam, lamp in products:
        t0 = time.perf_counter()
        tot, T = slices_product(lam, lamp, dps)
        wall = time.perf_counter() - t0
        out('  %s^2: %d sectors (m > 0), %d copies; %d sectors with n0 > 1 (null-vector start); copies NOT determined by (V3, V4a) alone (null space of C(v) of dimension > 1, a shared (v3, v4a) pair): %d; of the %d determined copies: max slice deviation %s, copies deviating by more than 1e-30: %d' % (lstr(lam), tot['sectors'], tot['copies'], tot['n0>1'], tot['ambiguous'], tot['copies'] - tot['ambiguous'], mp.nstr(tot['worst'], 3), tot['bad']))
        out('      time, whole product (all copies): G3 `copies` %.2f s;  closed form: matrices %.2f s + PRV eigenvalues (still G3 matrices) %.2f s + grading/plan %.2f s + recursion %.2f s = %.2f s   [wall incl. validation %.1f s]' % (T['copies'], T['mats'], T['eigen'], T['plan'], T['recursion'], T['closed_total'], wall))
        mus = [mu for mu, m, w in ng.sectors(lam, lamp) if m > 0]
        tm, tr, n, sk = g8a_timing(lam, lamp, mus, dps)
        out('      G8a slice_from_triple (G3 matrices, triple of operators, sectors with an n0 = 1 grading only; %d sectors skipped for n0 > 1): sector matrices + grading search %.2f s, recursion %.2f s (%d copies)' % (sk, tm, tr, n))
        res[(lam, lamp)] = (tot, T, (tm, tr, n))
    return res


# ======================================================================================== V4b within a path: what is established (parameters of J^{(4b)} = a' V + e [xi, V] + b' + c xi - (e/4) xi^2)
def v4b_sample(lam, kappa, A, B, m11, lamp):
    """exact (a', e, b', c) of the two-point path (A,B),(A+1,B-1) of kappa (g = -e/4 assumed, verified on all |P| >= 3 blocks by G8b): from the lower/upper off-diagonal ratios and the two diagonal entries"""
    R = gr.GLNRep.of(lam)
    pa = ((m11,), (A, B), tuple(kappa), tuple(lam))
    pb = ((m11,), (A + 1, B - 1), tuple(kappa), tuple(lam))
    ia, ib = R.index[pa], R.index[pb]
    nu = R.wt[ia]
    al, up, dn = gl3_entries(kappa, tuple(lamp[:3]), m11, A + B, A)
    al2, up2, dn2 = gl3_entries(kappa, tuple(lamp[:3]), m11, A + B, A + 1)
    V = [[al, dn2], [up, al2]]            # V[row][col]: V[1][0] = M[m+1, m] = up(A), V[0][1] = M[m, m+1] = down(A+1)
    op = so.reduced_op('V4b', tuple(lamp))
    J = [[Fr(0)] * 2 for _ in range(2)]
    for cc, ic in enumerate((ia, ib)):
        img = op.apply(R, {ic: Fr(1)})
        for rr, ir in enumerate((ia, ib)):
            J[rr][cc] = img.get(ir, Fr(0))
    x = A - B
    xi = [Fr(x * (x + 2)), Fr((x + 2) * (x + 4))]
    D = xi[1] - xi[0]
    rl = J[1][0] / V[1][0]
    ru = J[0][1] / V[0][1]
    e = (rl - ru) / (2 * D)
    a = (rl + ru) / 2
    g = -e / 4
    r0 = J[0][0] - a * V[0][0] - g * xi[0] ** 2
    r1 = J[1][1] - a * V[1][1] - g * xi[1] ** 2
    c = (r1 - r0) / (xi[1] - xi[0])
    return dict(nu=nu, a=a, e=e, b=r0 - c * xi[0], c=c)


def stage_v4b(deep=False, seed=29):
    hdr('V4b within a path: J^{(4b)} = a\' V^{(3)} + e [xi, V^{(3)}] + b\' + c xi - (e/4) xi^2   (what is established; the full closed form of a\', b\', c is NOT identified, see section 7)')
    rng = random.Random(seed)
    n = bad = 0
    t0 = time.perf_counter()
    for lam in LAMS6[:5] + HELD_A[:4]:
        for k in kappas_of(lam):
            if k[0] - k[1] < 1 or k[1] - k[2] < 1:
                continue
            for _try in range(30):
                A = rng.randint(k[1], k[0] - 1)
                B = rng.randint(k[2] + 1, k[1])
                m11 = rng.randint(B, A)
                lp = tuple(Fr(rng.randint(-20, 20), rng.choice([1, 2, 3])) for _ in range(4))
                try:
                    s = v4b_sample(lam, k, A, B, m11, lp)
                except ZeroDivisionError:
                    continue
                n += 1
                bad += (s['e'] != (1 - a_of(k, lam)) / 8)
                break
    out('  e(kappa; lam) = (1 - a(kappa; lam))/8 = (sum_i R_i)/8, g = -e/4:  %d (lam, kappa, weight space, lam\') samples, %d mismatches  [%.1f s]' % (n, bad, time.perf_counter() - t0))
    lam, kappa = (9, 6, 3, 0), (8, 5, 2)
    S = []
    for _ in range(400):
        A = rng.randint(kappa[1], kappa[0] - 1)
        B = rng.randint(kappa[2] + 1, kappa[1])
        m11 = rng.randint(B, A)
        lp = tuple(Fr(rng.randint(-20, 20), rng.choice([1, 2, 3])) for _ in range(4))
        try:
            s = v4b_sample(lam, kappa, A, B, m11, lp)
        except ZeroDivisionError:
            continue
        s['lamp'] = lp
        S.append(s)
        if len(S) >= (260 if deep else 70):
            break
    for key, degs in (('e', (0,)), ('a', (1,)), ('c', (2,)), ('b', (3, 4) if deep else (3,))):
        found = None
        for deg in degs:
            ms = [e for e in itertools.product(range(deg + 1), repeat=6) if sum(e) <= deg]
            if len(ms) + 5 > len(S):
                continue
            rows = [[_mono([Fr(s['nu'][0]), Fr(s['nu'][1])] + list(s['lamp']), e) for e in ms] for s in S]
            sol = solve_exact(rows, [s[key] for s in S])
            if sol is not None:
                found = (deg, len(ms), sol[1])
                break
        out('  kappa = (8,5,2) of [9,6,3,0]: %s as a polynomial in (nu1, nu2, lam\'1..4) of degree <= %s: %s  (%d samples)' % (
            {'e': 'e', 'a': 'a\'', 'c': 'c', 'b': 'b\''}[key], degs, ('exact fit of degree %d (%d monomials, rank %d)' % found) if found else ('no exact fit' + ('' if deep else ' with the lowest degrees tried / not run (deep)')), len(S)))
    return dict(e_samples=n, e_bad=bad)


def stage_determinacy():
    hdr('D4c: forward determinacy of the lattice recursion with the two closed-form operators (V3, V4a) alone')
    t0 = time.perf_counter()
    tot = two_op_determinacy()
    out('  weight spaces with w >= 2 of the 12 products: %d;  determined with (V3, V4a): %d (n0 = 1: %d);  determined with the triple (G8a): %d (n0 = 1: %d)  [%.1f s]' % (
        tot['spaces'], tot['det2'], tot['n0one2'], tot['det3'], tot['n0one3'], time.perf_counter() - t0))
    out('  of which sectors with m > 0: %d, determined with (V3, V4a): %d' % (tot['m>0'], tot['m>0_det2']))
    return tot


# ======================================================================================== CLI
def stage_slices9():
    return stage_slices(products=(((9, 6, 3, 0),) * 2,))


STAGES = {'a': stage_a, 'b': stage_b, 'words': stage_words, 'gl3chain': stage_gl3chain, 'v4b': lambda: stage_v4b(deep=True), 'mats': stage_mats, 'det': stage_determinacy, 'slices': stage_slices, 'slices9': stage_slices9}


def main(argv):
    if not argv or argv[0] == 'all':
        names = ['a', 'b', 'words', 'gl3chain', 'v4b', 'mats', 'det', 'slices', 'slices9']
    elif argv[0] == 'run':
        names = argv[1].split(',')
    else:
        names = argv[0].split(',')
    for nm in names:
        t0 = time.perf_counter()
        STAGES[nm]()
        out('  [stage %s: %.1f s]' % (nm, time.perf_counter() - t0))
    flush_tables()


if __name__ == '__main__':
    main(sys.argv[1:])
