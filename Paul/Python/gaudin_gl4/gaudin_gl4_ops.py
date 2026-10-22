"""gaudin_gl4_ops.py -- WP-G2: the two-site gl(N) Gaudin algebra, exactly: cdet(d_u - L(u)) -> S_1..S_N -> Laurent coefficients at theta_1 -> term lists,
numerical builders (exact Fractions / float64 / mod p), HW spaces of V_lam (x) V_lam', the canonical triple (V3, V4a, V4b) for gl(4).

Conventions (binding, brief SPEC-G2-gaudin-gl4-ops.md section 1 = Gaudin-sl3-CG-notes.md section 1.1 generalised verbatim)
  * [E_ij, E_kl] = delta_jk E_il - delta_li E_kj;  two sites V_lam (site 1) (x) V_lam' (site 2), e^s_ij = E_ij on site s, e^1 and e^2 commute;
  * untwisted Lax matrix  L(u)_ij = e^1_ij/(u-theta_1) + e^2_ij/(u-theta_2),  t = theta_12 = theta_1 - theta_2;
  * cdet(d_u - L) = sum_sigma sgn(sigma) M_{sigma(1) 1} M_{sigma(2) 2} ... M_{sigma(N) N},  M = d_u 1 - L(u), derivative moved to the right (d f = f d + f'),
        = d^N - S_1 d^(N-1) + S_2 d^(N-2) - ... + (-1)^N S_N     (operator products kept in the order written);
  * Laurent data at theta_1: x = u - theta_1 (so u - theta_2 = x + t):  S_k(u) = sum_j S_{k,j} x^(-j) + regular;  theta_2 data: swap the sites 1<->2 and t -> -t.
  * GT basis exactly as gl3_reps.gt_patterns (HW pattern = last index), plain rational non-unitary GT matrices acting on kets as columns,
    two-site flat index nu*dim2 + q (site 1 slow).

Layout
  1. symbolic noncommutative differential-operator algebra: cdet_S(N) -> {k: {word: coef}}, word = tuple of factors (i, j, d) = d-th u-derivative of l_ij  (stage A)
  2. the gl(3) reference formula eq. (4) of the notes (convention test)
  3. Laurent expansion at theta_1 (stage B): laurent(N) -> {k: {j: terms}}, terms = {(w1, w2): {power: int}}, w_s = tuple of (i, j) = the site-s factors in product order
  4. term-list algebra, PBW normal ordering in U(gl_N) (x) U(gl_N), invariants (n_s, C2_s, Omega, T^(s1..sr), antisymmetrised forms)
  5. GT representations of gl(N) (private scaffold from gl3_reps), Racah-Speiser, HW spaces (exact / float)
  6. two-site operators: exact (Fractions), float64 (orthonormal basis), mod p builders; restriction to HW spaces
  7. the stages (CLI: python -X utf8 gaudin_gl4_ops.py <stage>):  sym (convention test, Laurent data, residue relations)  |  scaffold (cross-check vs glN_reps)  |  gl3 (G2.1)  |  comm (G2.2)  |
     triple (G2.3, products up to [4,3,1,0]^2)  |  struct (G2.4)  |  terms (G2.5: writes ../../Mathematica/Gaudin/ref/gl4-gaudin-ops.txt)  |  big (G2.3 for [5,3,1,0]^2, [5,3,2,0]^2, [6,4,2,0]^2)  |  all
"""
import os, sys, time, itertools, math
from fractions import Fraction as Fr
from collections import defaultdict, Counter
from math import comb, factorial
from functools import lru_cache

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))                 # gl3_reps (read-only)
import numpy as np
import scipy.sparse as sps
import scipy.linalg as sla
import scipy.sparse.linalg as spla
import mpmath as mp
from gl3_reps import gt_patterns, _raise, _lower, _weight   # general-N GT formulas (read-only reuse)

REFDIR = os.path.normpath(os.path.join(HERE, '..', '..', 'Mathematica', 'Gaudin', 'ref'))
TERMS_FILE = os.path.join(REFDIR, 'gl4-gaudin-ops.txt')


def perm_sign(p):
    s = 1
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


# ====================================================================================================== 1. symbolic engine (stage A)
# nc polynomial = dict {word: int}, word = tuple of factors (i, j, d);  differential operator = dict {k: poly}  meaning  sum_k poly_k d^k
def _clean(d):
    return {k: v for k, v in d.items() if v != 0}


def p_add(P, Q, c=1):
    for w, v in Q.items():
        P[w] = P.get(w, 0) + c * v


def p_mul(P, Q):
    out = {}
    for w1, c1 in P.items():
        for w2, c2 in Q.items():
            w = w1 + w2
            out[w] = out.get(w, 0) + c1 * c2
    return out


def p_d(P):
    """u-derivative of an nc polynomial (Leibniz, factor by factor, order kept)"""
    out = {}
    for w, c in P.items():
        for k in range(len(w)):
            i, j, d = w[k]
            w2 = w[:k] + ((i, j, d + 1),) + w[k + 1:]
            out[w2] = out.get(w2, 0) + c
    return _clean(out)


def op_mul(A, B):
    """(A_a d^a)(B_b d^b) = sum_j C(a,j) A_a B_b^{(j)} d^{a+b-j}   (derivative moved to the right)"""
    out = {}
    towers = {b: [Pb] for b, Pb in B.items()}
    for a, Pa in A.items():
        for b, Pb in B.items():
            tw = towers[b]
            for j in range(a + 1):
                while len(tw) <= j:
                    tw.append(p_d(tw[-1]))
                tgt = out.setdefault(a + b - j, {})
                p_add(tgt, p_mul(Pa, tw[j]), comb(a, j))
    res = {}
    for k, v in out.items():
        v = _clean(v)
        if v:
            res[k] = v
    return res


def op_matrix_M(N):
    """M_ij = delta_ij d - l_ij  as differential operators"""
    M = {}
    for i in range(1, N + 1):
        for j in range(1, N + 1):
            M[(i, j)] = {1: {(): 1}, 0: {((i, i, 0),): -1}} if i == j else {0: {((i, j, 0),): -1}}
    return M


def cdet_S(N, variant='col'):
    """expand cdet(d - L) for gl(N): returns {m: S_m (nc polynomial in l_ij^{(d)})}, m = 1..N.
       variant 'col' : sum_sigma sgn(sigma) M_{sigma(1)1} ... M_{sigma(N)N}   (the binding convention)
               'row' : sum_sigma sgn(sigma) M_{1 sigma(1)} ... M_{N sigma(N)}
               'col-left' : column-ordered, derivative moved to the LEFT (d f = f d + f' is replaced by f d = d f - f', i.e. coefficients written as d^k C_k) """
    M = op_matrix_M(N)
    tot = {}
    for sg in itertools.permutations(range(1, N + 1)):
        s = perm_sign(sg)
        if variant in ('col', 'col-left'):
            facs = [M[(sg[c], c + 1)] for c in range(N)]
        elif variant == 'row':
            facs = [M[(c + 1, sg[c])] for c in range(N)]
        else:
            raise ValueError(variant)
        prod = facs[0]
        for f in facs[1:]:
            prod = op_mul(prod, f)
        for k, P in prod.items():
            p_add(tot.setdefault(k, {}), P, s)
    tot = {k: _clean(v) for k, v in tot.items()}
    if variant == 'col-left':
        tot = to_left_form(tot)
    assert tot.get(N) == {(): 1}, ('leading coefficient of d^N is not 1', variant, tot.get(N))
    return {m: {w: (-1) ** m * c for w, c in tot.get(N - m, {}).items()} for m in range(1, N + 1)}


def to_left_form(op):
    """rewrite sum_k C_k d^k (derivative right) as sum_k d^k C'_k ... returned in the SAME dict layout {k: poly} but meaning sum_k d^k poly_k
       (used only for the 'col-left' variant): C d^k = sum_j (-1)^j C(k,j) d^(k-j) C^{(j)}"""
    out = {}
    for k, P in op.items():
        tw = [P]
        for j in range(k + 1):
            while len(tw) <= j:
                tw.append(p_d(tw[-1]))
            tgt = out.setdefault(k - j, {})
            p_add(tgt, tw[j], (-1) ** j * comb(k, j))
    return {k: _clean(v) for k, v in out.items() if _clean(v)}


# ====================================================================================================== 2. eq. (4) of the gl(3) notes
def eq4_gl3():
    def l(i, j, d=0):
        return ((i, j, d),)

    def mono(*factors):
        w = ()
        for f in factors:
            w += f
        return w
    S1 = {l(1, 1): 1, l(2, 2): 1, l(3, 3): 1}
    S2 = {}
    for i in range(1, 4):
        for j in range(i + 1, 4):
            p_add(S2, {mono(l(i, i), l(j, j)): 1, mono(l(j, i), l(i, j)): -1})
    p_add(S2, {l(2, 2, 1): -1, l(3, 3, 1): -2})
    S3 = {}
    terms = [
        (l(3, 3, 2), 1), (mono(l(1, 1), l(3, 3, 1)), -1), (mono(l(2, 2), l(3, 3, 1)), -1), (mono(l(2, 2, 1), l(3, 3)), -1),
        (mono(l(3, 1), l(1, 3, 1)), 1), (mono(l(3, 2), l(2, 3, 1)), 1), (mono(l(3, 2, 1), l(2, 3)), 1),
        (mono(l(1, 1), l(2, 2), l(3, 3)), 1), (mono(l(1, 1), l(3, 2), l(2, 3)), -1), (mono(l(2, 1), l(1, 2), l(3, 3)), -1),
        (mono(l(2, 1), l(3, 2), l(1, 3)), 1), (mono(l(3, 1), l(1, 2), l(2, 3)), 1), (mono(l(3, 1), l(2, 2), l(1, 3)), -1)]
    for w, c in terms:
        S3[w] = S3.get(w, 0) + c
    return {1: _clean(S1), 2: _clean(S2), 3: _clean(S3)}


def fmt_word_l(w):
    return ' '.join("l%s%s%s" % (i, j, "'" * d if d <= 2 else "^(%d)" % d) for (i, j, d) in w) if w else '1'


def fmt_poly_l(P):
    keys = sorted(P, key=lambda w: (len(w), w))
    return ' '.join(('%+d' % P[w]) + '*' + fmt_word_l(w) for w in keys)


# ====================================================================================================== 3. Laurent expansion at theta_1 (stage B)
def comb_neg(B, n):
    """binomial coefficient C(-B, n)  (the coefficient of x^n in (1+x)^(-B))"""
    if B == 0:
        return 1 if n == 0 else 0
    return (-1) ** n * comb(B + n - 1, n)


def laurent_poly(Sk, k, jmin=0):
    """Laurent coefficients at theta_1 of an nc polynomial Sk = S_k:  l_ij^{(d)} = (-1)^d d! [ e1_ij x^(-d-1) + e2_ij (x+t)^(-d-1) ],
       (x+t)^(-B) = sum_n C(-B,n) t^(-B-n) x^n.   Returns {j: {(w1, w2): {power: int}}}, j = k..jmin (x^(-j) coefficient), w_s = tuple of (i,j) site-s factors in order."""
    res = {}
    for word, coef in Sk.items():
        pre = coef
        for (_, _, d) in word:
            pre *= (-1) ** d * factorial(d)
        r = len(word)
        for sites in itertools.product((1, 2), repeat=r):
            A = sum(f[2] + 1 for f, s in zip(word, sites) if s == 1)
            B = sum(f[2] + 1 for f, s in zip(word, sites) if s == 2)
            w1 = tuple((f[0], f[1]) for f, s in zip(word, sites) if s == 1)
            w2 = tuple((f[0], f[1]) for f, s in zip(word, sites) if s == 2)
            for n in range(0, A - jmin + 1):
                cb = comb_neg(B, n)
                if cb == 0:
                    continue
                j = A - n
                power = -B - n
                dd = res.setdefault(j, {}).setdefault((w1, w2), {})
                dd[power] = dd.get(power, 0) + pre * cb
    out = {}
    for j, T in res.items():
        T2 = {}
        for ww, pw in T.items():
            pw = {p: c for p, c in pw.items() if c != 0}
            if pw:
                T2[ww] = pw
        out[j] = T2
    return out


_LAUR_CACHE = {}


def laurent(N, variant='col'):
    """{k: {j: terms}}, k = 1..N, j = 0..k (j = 0: the regular part x^0);  terms = {(w1, w2): {power: int}}"""
    key = (N, variant)
    if key not in _LAUR_CACHE:
        S = cdet_S(N, variant)
        _LAUR_CACHE[key] = {k: laurent_poly(S[k], k) for k in range(1, N + 1)}
    return _LAUR_CACHE[key]


# ====================================================================================================== 4. term-list algebra, PBW normal ordering, invariants
# 'LT'  = Laurent terms  {(w1, w2): {power: int}}   (stage B output)
# 'UP'  = t-free polynomial in U(gl_N) (x) U(gl_N):  {(w1, w2): Fraction},  w_s = tuple of (i, j) in product order (leftmost factor acts last)
def lt_hat(T):
    """strip the (single) power of t from an LT: returns (UP, power).  S_{k,j} is homogeneous: S_{k,j} = t^(j-k) * hat"""
    pows = {p for pw in T.values() for p in pw}
    if not pows:
        return {}, 0
    assert len(pows) == 1, ('not homogeneous in t', pows)
    p = next(iter(pows))
    return {ww: Fr(pw[p]) for ww, pw in T.items()}, p


def up_clean(P):
    return {k: v for k, v in P.items() if v != 0}


def up_add(P, Q, c=1):
    for k, v in Q.items():
        P[k] = P.get(k, 0) + c * v
    return P


def up_lin(*pairs):
    """sum of c*P over pairs (c, P)"""
    out = {}
    for c, P in pairs:
        up_add(out, P, c)
    return up_clean(out)


def up_mul(P, Q):
    out = {}
    for (a1, a2), c1 in P.items():
        for (b1, b2), c2 in Q.items():
            k = (a1 + b1, a2 + b2)
            out[k] = out.get(k, 0) + c1 * c2
    return up_clean(out)


def up_swap(P):
    """site swap 1 <-> 2 (the theta_2 rule, without the t -> -t)"""
    return {(w2, w1): c for (w1, w2), c in P.items()}


def up_scale(P, c):
    return {k: v * c for k, v in P.items()} if c != 0 else {}


@lru_cache(maxsize=None)
def pbw_word(word):
    """normal ordering of a word of gl(N) generators (i,j) in the lexicographic order of (i,j): {sorted word: int}"""
    for k in range(len(word) - 1):
        if word[k] > word[k + 1]:
            a, b = word[k], word[k + 1]
            res = {}
            for w, c in pbw_word(word[:k] + (b, a) + word[k + 2:]).items():
                res[w] = res.get(w, 0) + c
            (i1, j1), (i2, j2) = a, b                       # [e_{i1 j1}, e_{i2 j2}] = delta_{j1 i2} e_{i1 j2} - delta_{j2 i1} e_{i2 j1}
            if j1 == i2:
                for w, c in pbw_word(word[:k] + ((i1, j2),) + word[k + 2:]).items():
                    res[w] = res.get(w, 0) + c
            if j2 == i1:
                for w, c in pbw_word(word[:k] + ((i2, j1),) + word[k + 2:]).items():
                    res[w] = res.get(w, 0) - c
            return {w: c for w, c in res.items() if c != 0}
    return {word: 1}


def up_pbw(P):
    """PBW normal form (each site separately; e^1 and e^2 commute): the canonical representative"""
    out = {}
    for (w1, w2), c in P.items():
        d1 = pbw_word(w1)
        d2 = pbw_word(w2)
        for u1, c1 in d1.items():
            for u2, c2 in d2.items():
                k = (u1, u2)
                out[k] = out.get(k, 0) + c * c1 * c2
    return up_clean(out)


def fmt_up(P, maxterms=60):
    """human-readable PBW-ordered polynomial; e1[ij] = site-1 generator E_ij, e2[ij] = site 2"""
    def fw(w1, w2):
        return ' '.join(['e1[%d%d]' % g for g in w1] + ['e2[%d%d]' % g for g in w2]) or '1'
    keys = sorted(P, key=lambda k: (-(len(k[0]) + len(k[1])), -len(k[0]), k))
    s = []
    for k in keys[:maxterms]:
        c = P[k]
        s.append('%s%s*%s' % ('+' if c > 0 else '-', str(abs(c)), fw(*k)))
    if len(keys) > maxterms:
        s.append('... (%d more terms)' % (len(keys) - maxterms))
    return ' '.join(s) if s else '0'


# ---- standard invariants (explicit definitions; N = number of rows)
def inv_one():
    return {((), ()): Fr(1)}


def _site_word(s, word):
    return (word, ()) if s == 1 else ((), word)


def inv_n(s, N):
    return {_site_word(s, ((i, i),)): Fr(1) for i in range(1, N + 1)}


def inv_C2(s, N):
    out = {}
    for i in range(1, N + 1):
        for j in range(1, N + 1):
            k = _site_word(s, ((i, j), (j, i)))
            out[k] = out.get(k, 0) + Fr(1)
    return out


def inv_Omega(N):
    """Omega = sum_ij e1_ij e2_ji"""
    return {(((i, j),), ((j, i),)): Fr(1) for i in range(1, N + 1) for j in range(1, N + 1)}


def inv_T(seq, N):
    """polarised trace T^{(s1..sr)} = sum_{i1..ir} e^{s1}_{i2 i1} e^{s2}_{i3 i2} ... e^{sr}_{i1 ir}  (CPV polarised traces; e.g. T^{(1,1,2)} = sum e1_ji e1_kj e2_ik)"""
    r = len(seq)
    out = {}
    for idx in itertools.product(range(1, N + 1), repeat=r):
        w1, w2 = [], []
        for c in range(r):
            g = (idx[(c + 1) % r], idx[c])
            (w1 if seq[c] == 1 else w2).append(g)
        k = (tuple(w1), tuple(w2))
        out[k] = out.get(k, 0) + 1
    return {k: Fr(v) for k, v in out.items()}


def antisym_form(N, a, b, I=None):
    """D^{(a,b)}_I = 1/(a! b!) sum_{sigma,tau in S(I)} sgn(sigma) sgn(tau) e1_{s1 t1} ... e1_{sa ta} e2_{s(a+1) t(a+1)} ... e2_{s(a+b) t(a+b)},  |I| = a+b;
       the (a,b)-polarisation of det(X+Y) on the index set I (I = {1..N} if None); signs relative to the sorted order of I"""
    I = tuple(range(1, N + 1)) if I is None else tuple(I)
    r = a + b
    assert len(I) == r
    out = {}
    perms = list(itertools.permutations(range(r)))
    sg = {p: perm_sign(p) for p in perms}
    for s in perms:
        for t in perms:
            w1 = tuple((I[s[c]], I[t[c]]) for c in range(a))
            w2 = tuple((I[s[c]], I[t[c]]) for c in range(a, r))
            k = (w1, w2)
            out[k] = out.get(k, 0) + sg[s] * sg[t]
    f = Fr(1, factorial(a) * factorial(b))
    return {k: v * f for k, v in out.items() if v != 0}


def antisym_sum_subsets(N, a, b):
    """sum over the (a+b)-subsets I of {1..N} of D^{(a,b)}_I   (for N = 3 and (2,1) this is the gl(3) vertex operator V of the notes)"""
    out = {}
    for I in itertools.combinations(range(1, N + 1), a + b):
        up_add(out, antisym_form(N, a, b, I))
    return up_clean(out)


def solve_exact(cols, target):
    """exact rational solution c of sum_i c_i cols[i] = target  (cols, target: UP polynomials, already PBW-normalised).
       Returns (c, rank, ninconsistent): c = None if the system is inconsistent (ninconsistent = number of rows that cannot be matched)."""
    keys = sorted({k for P in cols for k in P} | set(target))
    ridx = {k: i for i, k in enumerate(keys)}
    n = len(cols)
    rows = [dict() for _ in keys]
    for j, P in enumerate(cols):
        for k, v in P.items():
            rows[ridx[k]][j] = Fr(v)
    for k, v in target.items():
        rows[ridx[k]][n] = Fr(v)
    piv = {}                                              # col -> fully reduced row with pivot coefficient 1
    incons = 0
    for r in rows:
        r = dict(r)
        while True:
            cs = [c for c in r if c < n and c in piv]
            if not cs:
                break
            c = cs[0]
            f = r[c]
            for cc, v in piv[c].items():
                r[cc] = r.get(cc, 0) - f * v
            r = {cc: v for cc, v in r.items() if v != 0}
        free = [c for c in r if c < n]
        if free:
            c = min(free)
            f = r[c]
            r = {cc: v / f for cc, v in r.items()}
            for c2 in list(piv):                          # keep the pivot rows fully reduced
                if c in piv[c2]:
                    f2 = piv[c2][c]
                    for cc, v in r.items():
                        piv[c2][cc] = piv[c2].get(cc, 0) - f2 * v
                    piv[c2] = {cc: v for cc, v in piv[c2].items() if v != 0}
            piv[c] = r
        elif r.get(n, 0) != 0:
            incons += 1
    rank = len(piv)
    if incons:
        return None, rank, incons
    sol = [Fr(0)] * n
    for c, r in piv.items():
        sol[c] = r.get(n, Fr(0))
    return sol, rank, 0


# ====================================================================================================== 5. GT representations of gl(N), Racah-Speiser, HW spaces
def _frac(c):
    return Fr(int(c.p), int(c.q))


def _compose(A, B):
    """(A o B)[c] = A applied to (B e_c);  A, B = lists of dicts row -> Fraction (columns)"""
    out = []
    for col in B:
        acc = {}
        for r, v in col.items():
            for r2, v2 in A[r].items():
                acc[r2] = acc.get(r2, 0) + v * v2
        out.append({r: v for r, v in acc.items() if v != 0})
    return out


def _sub(A, B):
    out = []
    for ca, cb in zip(A, B):
        acc = dict(ca)
        for r, v in cb.items():
            acc[r] = acc.get(r, 0) - v
        out.append({r: v for r, v in acc.items() if v != 0})
    return out


def _addm(A, B, c=1):
    out = []
    for ca, cb in zip(A, B):
        acc = dict(ca)
        for r, v in cb.items():
            acc[r] = acc.get(r, 0) + c * v
        out.append({r: v for r, v in acc.items() if v != 0})
    return out


def _lincomb_cols(pairs, n):
    """sum of c*M over pairs (c, M), M = list of column dicts of length n"""
    out = [dict() for _ in range(n)]
    for c, M in pairs:
        for j in range(n):
            o = out[j]
            for r, v in M[j].items():
                o[r] = o.get(r, 0) + c * v
    return [{r: v for r, v in o.items() if v != 0} for o in out]


def weight_mults(lam):
    N = len(lam)
    return Counter(tuple(_weight(k, p) for k in range(1, N + 1)) for p in gt_patterns(tuple(lam)))


def weyl_dim(lam):
    N = len(lam)
    d = Fr(1)
    for i in range(N):
        for j in range(i + 1, N):
            d *= Fr(lam[i] - lam[j] + j - i, j - i)
    assert d.denominator == 1
    return int(d)


def casimir2(lam):
    """C_2 = sum_ij E_ij E_ji = sum_i lam_i (lam_i + N + 1 - 2 i)  (scalar on V_lam)"""
    N = len(lam)
    return sum(lam[i] * (lam[i] + N + 1 - 2 * (i + 1)) for i in range(N))


def racah_speiser(lam, lam2=None):
    """decomposition {mu: multiplicity} of V_lam (x) V_lam2 (default lam2 = lam): Klimyk/Racah-Speiser summing over the GT weights of V_lam, rho = (N-1,...,0)"""
    lam = tuple(lam)
    lam2 = lam if lam2 is None else tuple(lam2)
    N = len(lam)
    rho = tuple(range(N - 1, -1, -1))
    wm = weight_mults(lam)
    acc = defaultdict(int)
    for w, m in wm.items():
        x = tuple(lam2[i] + w[i] + rho[i] for i in range(N))
        if len(set(x)) < N:
            continue
        order = sorted(range(N), key=lambda i: -x[i])
        sgn = perm_sign(order)
        mu = tuple(x[order[i]] - rho[i] for i in range(N))
        acc[mu] += sgn * m
    out = {}
    for mu, m in acc.items():
        if m < 0:
            raise AssertionError(('negative multiplicity', mu, m))
        if m > 0:
            out[mu] = m
    return out


class Rep:
    """exact sparse GT irrep of gl(N) in the gl3_reps order (HW pattern = last index): E[(i,j)] = list over columns of dicts {row: Fraction} for ALL i, j
       (diagonal included).  Private scaffold from gl3_reps._raise/_lower/_weight (general N), non-adjacent generators by commutators."""
    _cache = {}

    @classmethod
    def of(cls, lam):
        lam = tuple(int(x) for x in lam)
        if lam not in cls._cache:
            cls._cache[lam] = cls(lam)
        return cls._cache[lam]

    def __init__(self, lam):
        lam = tuple(int(x) for x in lam)
        N = len(lam)
        if not (all(lam[i] >= lam[i + 1] for i in range(N - 1)) and lam[-1] >= 0):
            raise ValueError(lam)
        self.lam, self.N = lam, N
        self.pats = gt_patterns(lam)
        self.dim = n = len(self.pats)
        self.index = {p: i for i, p in enumerate(self.pats)}
        self.wt = [tuple(_weight(k, p) for k in range(1, N + 1)) for p in self.pats]
        d = defaultdict(list)
        for i, w in enumerate(self.wt):
            d[w].append(i)
        self.wsp = {w: np.array(v, dtype=np.int64) for w, v in d.items()}
        self.top = self.index[tuple(tuple(lam[:k]) for k in range(1, N + 1))]
        assert self.top == n - 1, 'HW pattern must be the last index'
        E = {}
        for k in range(1, N):
            U, D = [], []
            for p in self.pats:
                u, dd = {}, {}
                for tg, c in _raise(k, p).items():
                    j = self.index.get(tg)               # raw GT formula can propose targets outside the interlacing range: drop (as GL3Rep does)
                    if j is not None:
                        u[j] = _frac(c)
                for tg, c in _lower(k, p).items():
                    j = self.index.get(tg)
                    if j is not None:
                        dd[j] = _frac(c)
                U.append(u)
                D.append(dd)
            E[(k, k + 1)] = U
            E[(k + 1, k)] = D
        for k in range(1, N + 1):
            E[(k, k)] = [({c: Fr(self.wt[c][k - 1])} if self.wt[c][k - 1] != 0 else {}) for c in range(n)]
        for span in range(2, N):
            for i in range(1, N - span + 1):
                j = i + span
                a = _compose(E[(i, j - 1)], E[(j - 1, j)])
                b = _compose(E[(j - 1, j)], E[(i, j - 1)])
                E[(i, j)] = _sub(a, b)                                    # [E_{i,j-1}, E_{j-1,j}] = E_ij
                a = _compose(E[(j, j - 1)], E[(j - 1, i)])
                b = _compose(E[(j - 1, i)], E[(j, j - 1)])
                E[(j, i)] = _sub(a, b)                                    # [E_{j,j-1}, E_{j-1,i}] = E_ji
        self.E = E
        self._g = None
        self._Eh = None
        self._wc = {(): [{c: Fr(1)} for c in range(n)]}
        self._wcsr = {}

    # ---- invariant metric: diagonal G with E_ij^T G = G E_ji, G[top] = 1, by a walk from the top through the simple raising operators (exact)
    @property
    def g(self):
        if self._g is None:
            n, N = self.dim, self.N
            g = [None] * n
            g[self.top] = Fr(1)
            rev = defaultdict(list)
            for k in range(1, N):
                U, D = self.E[(k, k + 1)], self.E[(k + 1, k)]
                for q in range(n):
                    for r, v in U[q].items():
                        rev[r].append((q, v, D[r].get(q, 0)))
            stack = [self.top]
            while stack:
                r = stack.pop()
                for q, v, w in rev[r]:
                    if g[q] is None and w != 0:
                        g[q] = g[r] * v / w
                        stack.append(q)
            assert all(x is not None and x > 0 for x in g)
            self._g = g
        return self._g

    def check_metric(self):
        """exact: g_r E_ij[r,q] = g_q E_ji[q,r] for ALL off-diagonal pairs; returns the number of violations"""
        g = self.g
        bad = 0
        N = self.N
        for i in range(1, N + 1):
            for j in range(i + 1, N + 1):
                U, D = self.E[(i, j)], self.E[(j, i)]
                for q in range(self.dim):
                    for r, v in U[q].items():
                        if g[r] * v != g[q] * D[r].get(q, 0):
                            bad += 1
                    for r, w in D[q].items():
                        if g[q] * U[r].get(q, 0) != g[r] * w:
                            bad += 1
        return bad

    def check_commutators(self):
        """exact: all N^4 relations [E_ij, E_kl] = delta_jk E_il - delta_li E_kj on the matrices; returns the number of violations"""
        N = self.N
        bad = 0
        for i, j, k, l in itertools.product(range(1, N + 1), repeat=4):
            lhs = _sub(_compose(self.E[(i, j)], self.E[(k, l)]), _compose(self.E[(k, l)], self.E[(i, j)]))
            rhs = [dict() for _ in range(self.dim)]
            if j == k:
                rhs = _addm(rhs, self.E[(i, l)], 1)
            if l == i:
                rhs = _addm(rhs, self.E[(k, j)], -1)
            if lhs != rhs:
                bad += 1
        return bad

    def sqrt_metric(self):
        if not hasattr(self, '_sg'):
            self._sg = np.array([math.sqrt(float(x)) for x in self.g])
        return self._sg

    def float_ops(self):
        """CSR matrices of all E_ij (incl. diagonal) in the ORTHONORMAL basis: Eh[r,q] = E[r,q] sqrt(g_r/g_q);  Eh_ji = Eh_ij^T"""
        if self._Eh is None:
            n = self.dim
            g = self.g
            Eh = {}
            for ij, cols in self.E.items():
                rows_, cols_, vals_ = [], [], []
                for q in range(n):
                    for r, v in cols[q].items():
                        rows_.append(r)
                        cols_.append(q)
                        vals_.append(float(v) * math.sqrt(float(g[r] / g[q])))
                Eh[ij] = sps.csr_matrix((vals_, (rows_, cols_)), shape=(n, n))
            self._Eh = Eh
        return self._Eh

    # ---- products of generators (word = tuple of (i,j); leftmost factor acts last)
    def word_cols(self, word):
        word = tuple(word)
        if word not in self._wc:
            self._wc[word] = _compose(self.E[word[0]], self.word_cols(word[1:]))
        return self._wc[word]

    def word_csr(self, word):
        word = tuple(word)
        if word not in self._wcsr:
            if not word:
                self._wcsr[word] = sps.identity(self.dim, format='csr')
            else:
                self._wcsr[word] = (self.float_ops()[word[0]] @ self.word_csr(word[1:])).tocsr()
        return self._wcsr[word]

    def apply(self, ij, vec):
        out = {}
        cols = self.E[ij]
        for c, v in vec.items():
            for r, e in cols[c].items():
                out[r] = out.get(r, 0) + v * e
        return {r: v for r, v in out.items() if v != 0}


def sum_dim_check(lam1, lam2, rs=None):
    rs = rs or racah_speiser(lam1, lam2)
    return sum(m * weyl_dim(mu) for mu, m in rs.items()) == weyl_dim(lam1) * weyl_dim(lam2)


def nW_of(lam1, lam2, mu):
    """dimension of the weight-mu subspace of V_lam1 (x) V_lam2"""
    N = len(mu)
    w1 = weight_mults(lam1)
    w2 = weight_mults(lam2)
    return sum(m * w2.get(tuple(mu[k] - w[k] for k in range(N)), 0) for w, m in w1.items())


# ---------------------------------------------------------------------------------------------------- HW spaces of V_lam (x) V_lam'
def hw_pairs(rep1, rep2, mu):
    """flat list of pairs (nu, q) (site 1 slow) of total weight mu"""
    N = rep1.N
    pairs = []
    for nu in range(rep1.dim):
        qw = tuple(mu[k] - rep1.wt[nu][k] for k in range(N))
        qs = rep2.wsp.get(qw)
        if qs is not None:
            pairs.extend((nu, int(q)) for q in qs)
    return pairs


def hw_rows(rep1, rep2, pairs, vals='exact'):
    """rows of the stacked total simple raising operators E_{i,i+1}^{(1)} + E_{i,i+1}^{(2)} on the weight-mu subspace: list of dicts {pair index: value}
       (vals = 'exact': plain GT Fractions; 'float': ORTHONORMAL-basis floats)"""
    N = rep1.N
    rows = {}
    if vals == 'exact':
        S1 = {i: rep1.E[(i, i + 1)] for i in range(1, N)}
        S2 = {i: rep2.E[(i, i + 1)] for i in range(1, N)}
    else:
        Eh1, Eh2 = rep1.float_ops(), rep2.float_ops()
        S1, S2 = {}, {}
        for i in range(1, N):
            for S, Eh, rep in ((S1, Eh1, rep1), (S2, Eh2, rep2)):
                C = Eh[(i, i + 1)].tocsc()
                S[i] = [{int(C.indices[k]): float(C.data[k]) for k in range(C.indptr[c], C.indptr[c + 1])} for c in range(rep.dim)]
    for ci, (nu, q) in enumerate(pairs):
        for i in range(1, N):
            for r, v in S1[i][nu].items():
                d = rows.setdefault((i, r, q), {})
                d[ci] = d.get(ci, 0) + v
            for r, v in S2[i][q].items():
                d = rows.setdefault((i, nu, r), {})
                d[ci] = d.get(ci, 0) + v
    return list(rows.values())


def nullspace_exact(rows, ncols):
    """exact sparse RREF nullspace: canonical basis (1 at a free column, 0 at the other free columns).  Returns (basis [dict col->Fr], free [cols])"""
    piv = {}
    for r in rows:
        r = dict(r)
        for c in [c for c in r if c in piv]:
            f = r.get(c, 0)
            if f != 0:
                for cc, v in piv[c].items():
                    r[cc] = r.get(cc, 0) - f * v
        r = {c: v for c, v in r.items() if v != 0}
        if not r:
            continue
        c = min(r)
        f = r[c]
        r = {cc: v / f for cc, v in r.items()}
        for c2 in list(piv):
            r2 = piv[c2]
            if c in r2:
                f2 = r2[c]
                for cc, v in r.items():
                    r2[cc] = r2.get(cc, 0) - f2 * v
                piv[c2] = {cc: v for cc, v in r2.items() if v != 0}
        piv[c] = r
    free = [c for c in range(ncols) if c not in piv]
    basis = []
    for f in free:
        v = {f: Fr(1)}
        for c, row in piv.items():
            x = row.get(f, 0)
            if x != 0:
                v[c] = -x
        basis.append(v)
    return basis, free


class HWExact:
    """exact HW space of weight mu in V_lam (x) V_lam' (plain GT basis): .pairs, .basis[k] = dict {(nu,q): Fr}, .free (pair indices), .m, .nW"""
    def __init__(self, rep1, rep2, mu):
        self.rep1, self.rep2, self.mu = rep1, rep2, tuple(mu)
        self.pairs = hw_pairs(rep1, rep2, mu)
        self.nW = len(self.pairs)
        rows = hw_rows(rep1, rep2, self.pairs, 'exact')
        b, self.free = nullspace_exact(rows, self.nW)
        self.m = len(b)
        self.basis = [{self.pairs[c]: v for c, v in vec.items()} for vec in b]

    def gram(self):
        """exact Gram matrix of the basis for the invariant form G (x) G (diagonal metric of each site)"""
        g1, g2 = self.rep1.g, self.rep2.g
        m = self.m
        G = [[Fr(0)] * m for _ in range(m)]
        for i in range(m):
            for j in range(i, m):
                s = Fr(0)
                bi, bj = self.basis[i], self.basis[j]
                if len(bi) > len(bj):
                    bi, bj = bj, bi
                for (nu, q), v in bi.items():
                    w = bj.get((nu, q))
                    if w is not None:
                        s += g1[nu] * g2[q] * v * w
                G[i][j] = G[j][i] = s
        return G


class HWFloat:
    """float HW space (orthonormal two-site basis): .pairs, .N (nW x m orthonormal null vectors), .m, .gap = (smallest non-null sv, largest null sv)"""
    def __init__(self, rep1, rep2, mu, m_expected):
        self.rep1, self.rep2, self.mu = rep1, rep2, tuple(mu)
        self.pairs = hw_pairs(rep1, rep2, mu)
        self.nW = nW = len(self.pairs)
        rows = hw_rows(rep1, rep2, self.pairs, 'float')
        A = np.zeros((len(rows), nW))
        for r, d in enumerate(rows):
            for c, v in d.items():
                A[r, c] = v
        if A.shape[0] < nW:
            A = np.vstack([A, np.zeros((nW - A.shape[0], nW))])
        U, s, Vt = sla.svd(A, full_matrices=False)
        m = m_expected
        self.m = m
        self.N = Vt[nW - m:, :].T.copy() if m > 0 else np.zeros((nW, 0))
        self.gap = (float(s[nW - m - 1]) if nW - m - 1 >= 0 else np.inf, float(s[nW - m]) if m > 0 and nW - m < len(s) else 0.0)
        self.s = s

    def dense(self, k):
        n1, n2 = self.rep1.dim, self.rep2.dim
        X = np.zeros((n1, n2))
        for c, (nu, q) in enumerate(self.pairs):
            X[nu, q] = self.N[c, k]
        return X


# ====================================================================================================== 6. two-site operators (exact / float64 / mod p), restriction to HW spaces
def lt_at(T, t):
    """numeric t (Fraction) in a Laurent-term list: {(w1,w2): Fraction}"""
    t = Fr(t)
    out = {}
    for ww, pw in T.items():
        v = sum(Fr(c) * t ** p for p, c in pw.items())
        if v != 0:
            out[ww] = v
    return out


def frac_mod(x, p):
    x = Fr(x)
    return x.numerator % p * pow(x.denominator % p, -1, p) % p


class Op2:
    """operator O = sum_{(w1,w2)} c_{w1 w2} (word w1 on site 1)(word w2 on site 2) on V_lam (x) V_lam' (c = Fraction, UP form, t already substituted),
       stored in Kronecker-sum form  O = sum_g P_g (x) Q_g  by grouping on the side with fewer distinct words (the other side's words are combined into one matrix per group)."""
    def __init__(self, up, rep1, rep2, scale=1):
        self.rep1, self.rep2 = rep1, rep2
        self.n1, self.n2 = rep1.dim, rep2.dim
        scale = Fr(scale)
        up = {k: v * scale for k, v in up.items() if v != 0}
        self.up = up
        d1 = {w1 for (w1, w2) in up}
        d2 = {w2 for (w1, w2) in up}
        self.group_site = 2 if len(d2) <= len(d1) else 1          # the site whose words label the groups (kept as plain word matrices)
        grp = defaultdict(dict)
        for (w1, w2), c in up.items():
            if self.group_site == 2:
                grp[w2][w1] = grp[w2].get(w1, 0) + c
            else:
                grp[w1][w2] = grp[w1].get(w2, 0) + c
        self.groups = grp
        self._exact = None
        self._float = None

    # ---- exact: list of (P_cols, Q_cols)
    def exact_pairs(self):
        if self._exact is None:
            out = []
            for kw, comb_ in self.groups.items():
                if self.group_site == 2:
                    P = _lincomb_cols([(c, self.rep1.word_cols(w)) for w, c in comb_.items()], self.n1)
                    Q = self.rep2.word_cols(kw)
                else:
                    P = self.rep1.word_cols(kw)
                    Q = _lincomb_cols([(c, self.rep2.word_cols(w)) for w, c in comb_.items()], self.n2)
                out.append((P, Q))
            self._exact = out
        return self._exact

    def apply_exact(self, X):
        """X = dict {(nu,q): Fraction} -> dict"""
        out = {}
        for P, Q in self.exact_pairs():
            for (nu, q), v in X.items():
                cp, cq = P[nu], Q[q]
                if not cp or not cq:
                    continue
                for r, a in cp.items():
                    va = v * a
                    for s, b in cq.items():
                        k = (r, s)
                        out[k] = out.get(k, 0) + va * b
        return {k: v for k, v in out.items() if v != 0}

    def matrix_exact(self):
        """dict-of-columns exact matrix on the flat index nu*n2 + q (list of dicts)"""
        n1, n2 = self.n1, self.n2
        cols = []
        for nu in range(n1):
            for q in range(n2):
                col = {}
                for P, Q in self.exact_pairs():
                    cp, cq = P[nu], Q[q]
                    for r, a in cp.items():
                        for s, b in cq.items():
                            k = r * n2 + s
                            col[k] = col.get(k, 0) + a * b
                cols.append({k: v for k, v in col.items() if v != 0})
        return cols

    # ---- float64 in the orthonormal basis: list of (P dense, Q csr)
    def float_pairs(self):
        if self._float is None:
            out = []
            for kw, comb_ in self.groups.items():
                if self.group_site == 2:
                    P = sum(float(c) * self.rep1.word_csr(w) for w, c in comb_.items()).toarray()
                    Q = self.rep2.word_csr(kw)
                else:
                    P = self.rep1.word_csr(kw).toarray()
                    Q = sum(float(c) * self.rep2.word_csr(w) for w, c in comb_.items()).tocsr()
                out.append((P, sps.csr_matrix(Q)))
            self._float = out
        return self._float

    def apply_float(self, X):
        """X dense (n1, n2) (orthonormal components) -> dense"""
        out = np.zeros((self.n1, self.n2))
        for P, Q in self.float_pairs():
            Y = (Q @ X.T).T                         # X Q^T
            out += P @ Y
        return out

    def matrix_float(self):
        tot = None
        for P, Q in self.float_pairs():
            term = sps.kron(sps.csr_matrix(P), Q, format='csr')
            tot = term if tot is None else tot + term
        return tot.tocsr() if tot is not None else sps.csr_matrix((self.n1 * self.n2,) * 2)

    # ---- mod p (dense int64), exact residues
    def matrix_modp(self, p):
        n1, n2 = self.n1, self.n2
        tot = np.zeros((n1 * n2, n1 * n2), dtype=np.int64)
        for P, Q in self.exact_pairs():
            Pm = np.zeros((n1, n1), dtype=np.int64)
            Qm = np.zeros((n2, n2), dtype=np.int64)
            for c, col in enumerate(P):
                for r, v in col.items():
                    Pm[r, c] = frac_mod(v, p)
            for c, col in enumerate(Q):
                for r, v in col.items():
                    Qm[r, c] = frac_mod(v, p)
            tot = (tot + np.kron(Pm, Qm) % p) % p
        return tot


def op_matrix(terms, rep1, rep2, t_value, kind='float'):
    """matrix of the operator given by a Laurent-term list `terms` ({(w1,w2): {power: int}}) at the value t = t_value on V_lam (x) V_lam' (flat index nu*n2+q).
       kind = 'float' (scipy CSR, orthonormal GT basis) | 'exact' (list of dicts, plain GT basis)"""
    op = Op2(lt_at(terms, t_value), rep1, rep2)
    return op.matrix_float() if kind == 'float' else op.matrix_exact()


# ---- total generators on the two-site space
def tot_generators_exact(rep1, rep2):
    n1, n2 = rep1.dim, rep2.dim
    N = rep1.N
    out = {}
    for i in range(1, N + 1):
        for j in range(1, N + 1):
            cols = []
            for nu in range(n1):
                for q in range(n2):
                    col = {}
                    for r, v in rep1.E[(i, j)][nu].items():
                        col[r * n2 + q] = col.get(r * n2 + q, 0) + v
                    for r, v in rep2.E[(i, j)][q].items():
                        col[nu * n2 + r] = col.get(nu * n2 + r, 0) + v
                    cols.append(col)
            out[(i, j)] = cols
    return out


def tot_generators_float(rep1, rep2):
    N = rep1.N
    E1, E2 = rep1.float_ops(), rep2.float_ops()
    I1, I2 = sps.identity(rep1.dim, format='csr'), sps.identity(rep2.dim, format='csr')
    out = {}
    for i in range(1, N + 1):
        for j in range(1, N + 1):
            A = E1.get((i, j))
            if A is None:                                     # diagonal generators are not stored by float_ops of E-dict? (they are, E has (k,k))
                raise KeyError((i, j))
            out[(i, j)] = (sps.kron(A, I2) + sps.kron(I1, E2[(i, j)])).tocsr()
    return out


def mat_equal_exact(A, B):
    return all(a == b for a, b in zip(A, B))


def comm_exact(A, B):
    """[A,B] exact (dict-of-columns); returns the commutator"""
    return _sub(_compose(A, B), _compose(B, A))


def nnz_exact(A):
    return sum(len(c) for c in A)


def rel_comm_float(A, B):
    C = (A @ B - B @ A)
    nC = spla.norm(C) if C.nnz else 0.0
    nA, nB = spla.norm(A), spla.norm(B)
    den = nA * nB
    return (nC / den) if den > 0 else nC


# ---- restriction to HW spaces
def restrict_exact(op, hw):
    """exact m x m matrix R (list of lists of Fractions, column j = coordinates of O b_j in the canonical HW basis) and the number of basis vectors for which O b_j
       is NOT exactly sum_i R_ij b_i (must be 0: the operator preserves the HW space)"""
    m = hw.m
    cols = []
    viol = 0
    frees = [hw.pairs[f] for f in hw.free]
    for b in hw.basis:
        y = op.apply_exact(b)
        co = [y.get(k, Fr(0)) for k in frees]
        rec = {}
        for c, bi in zip(co, hw.basis):
            if c != 0:
                for k, v in bi.items():
                    rec[k] = rec.get(k, 0) + c * v
        rec = {k: v for k, v in rec.items() if v != 0}
        if rec != y:
            viol += 1
        cols.append(co)
    return [[cols[j][i] for j in range(m)] for i in range(m)], viol


def restrict_float(op, hw):
    """float m x m matrix N^T O N (orthonormal HW basis) and the leakage = max ||O N_k - N N^T O N_k||"""
    m = hw.m
    R = np.zeros((m, m))
    leak = 0.0
    for k in range(m):
        Y = op.apply_float(hw.dense(k))
        y = np.array([Y[nu, q] for (nu, q) in hw.pairs])
        c = hw.N.T @ y
        R[:, k] = c
        tot = np.linalg.norm(Y)
        res = np.linalg.norm(y - hw.N @ c)
        outside = np.sqrt(max(tot ** 2 - np.linalg.norm(y) ** 2, 0.0))
        leak = max(leak, max(res, outside) / max(tot, 1e-300))
    return R, leak


def mat_fr_mul(A, B):
    n, k, m = len(A), len(B), len(B[0]) if B else 0
    return [[sum((A[i][l] * B[l][j] for l in range(k)), Fr(0)) for j in range(m)] for i in range(n)]


def mat_fr_comm(A, B):
    AB, BA = mat_fr_mul(A, B), mat_fr_mul(B, A)
    return [[AB[i][j] - BA[i][j] for j in range(len(AB[0]))] for i in range(len(AB))]


def mat_fr_zero(A):
    return all(x == 0 for r in A for x in r)


def mat_fr_transpose(A):
    return [list(r) for r in zip(*A)]


# ====================================================================================================== 7. operators of the two-site model on a product, sector analysis
def laurent_up(N, variant='col'):
    """{(k,j): (hat UP, power p)} for k = 1..N, j = 0..k :  S_{k,j} = t^p * hat  (theta_1 data)"""
    L = laurent(N, variant)
    out = {}
    for k in range(1, N + 1):
        for j in range(0, k + 1):
            if j in L[k]:
                out[(k, j)] = lt_hat(L[k][j])
            else:
                out[(k, j)] = ({}, j - k)
    return out


def build_ops(R1, R2, N, t, variant='col'):
    """all theta_1 / theta_2 Laurent coefficients S_{k,j}, k = 2..N, j = 1..k, as Op2 at the numeric t (theta_2: sites swapped, t -> -t); keys (k, j, s)"""
    t = Fr(t)
    LU = laurent_up(N, variant)
    ops = {}
    for k in range(2, N + 1):
        for j in range(1, k + 1):
            hat, p = LU[(k, j)]
            ops[(k, j, 1)] = Op2(hat, R1, R2, scale=t ** p)
            ops[(k, j, 2)] = Op2(up_swap(hat), R1, R2, scale=(-t) ** p)
    return ops


def triple_ops(R1, R2):
    """the canonical triple of gl(4): V3 = t S_{3,2}, V4a = t S_{4,3}, V4b = t^2 S_{4,2}  (t-free hats)"""
    LU = laurent_up(4)
    out = {}
    for name, (k, j) in (('V3', (3, 2)), ('V4a', (4, 3)), ('V4b', (4, 2))):
        hat, p = LU[(k, j)]
        assert p == j - k
        out[name] = Op2(hat, R1, R2, scale=1)
    return out


def nonzero_exact_ops(mats):
    return [k for k, M in mats.items() if nnz_exact(M)]


def exact_charpoly(R):
    """characteristic polynomial coefficients [1, c_{m-1}, ..., c_0] (Faddeev-LeVerrier, exact Fractions) of a rational square matrix"""
    m = len(R)
    if m == 0:
        return [Fr(1)]
    coeffs = [Fr(1)]
    Mk = [[Fr(0)] * m for _ in range(m)]
    I = [[Fr(1) if i == j else Fr(0) for j in range(m)] for i in range(m)]
    for k in range(1, m + 1):
        # M_k = R M_{k-1} + c_{m-k+1} I
        Mk = mat_fr_mul(R, Mk) if k > 1 else [[Fr(0)] * m for _ in range(m)]
        Mk = [[Mk[i][j] + coeffs[-1] * I[i][j] for j in range(m)] for i in range(m)]
        RM = mat_fr_mul(R, Mk)
        c = -sum((RM[i][i] for i in range(m)), Fr(0)) / k
        coeffs.append(c)
    return coeffs


def poly_disc_and_real(coeffs):
    """exact discriminant (sympy) and the number of distinct real roots (Sturm) of the monic polynomial with rational coefficients (high degree first)"""
    import sympy
    x = sympy.Symbol('x')
    P = sympy.Poly([sympy.Rational(c.numerator, c.denominator) for c in coeffs], x)
    if P.degree() <= 0:
        return sympy.Integer(1), 0
    return P.discriminant(), P.count_roots()


def fr_mp(R):
    return mp.matrix([[mp.mpf(x.numerator) / x.denominator for x in row] for row in R])


def joint_spectrum(Rs, combo, dps=40):
    """Rs = list of exact commuting m x m rational matrices; combo = list of rational weights (generic combination M = sum combo_i R_i).
       Returns dict: eigenvalues of M (sorted), the joint eigenvalue tuples (one per eigenvector of M), the max off-diagonal entry of U^-1 R_i U (joint diagonalisability),
       the max imaginary part, min gap of each R_i alone (from its own eigenvalues), the min joint distance."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        m = len(Rs[0])
        M = [[sum((c * R[i][j] for c, R in zip(combo, Rs)), Fr(0)) for j in range(m)] for i in range(m)]
        Mm = fr_mp(M)
        E, U = mp.eig(Mm)
        order = sorted(range(m), key=lambda k: mp.re(E[k]))
        U = mp.matrix([[U[i, k] for k in order] for i in range(m)])
        Uinv = mp.inverse(U)
        imag = max([abs(mp.im(E[k])) for k in range(m)] + [mp.mpf(0)])
        trip = [[None] * len(Rs) for _ in range(m)]
        offmax = mp.mpf(0)
        for a, R in enumerate(Rs):
            D = Uinv * fr_mp(R) * U
            for i in range(m):
                trip[i][a] = D[i, i]
                for j in range(m):
                    if i != j:
                        offmax = max(offmax, abs(D[i, j]))
        imag = max([imag] + [abs(mp.im(x)) for row in trip for x in row])
        gaps = []
        for a in range(len(Rs)):
            ev = sorted(mp.re(trip[i][a]) for i in range(m))
            gaps.append(min((ev[i + 1] - ev[i] for i in range(m - 1)), default=mp.inf))
        jd = min((max(abs(mp.re(trip[i][a]) - mp.re(trip[j][a])) for a in range(len(Rs))) for i in range(m) for j in range(i + 1, m)), default=mp.inf)
        return dict(trip=[[mp.re(x) for x in row] for row in trip], offmax=offmax, imag=imag, gaps=gaps, jgap=jd, m=m)
    finally:
        mp.mp.dps = old


def fmt_mp(x, n=12):
    if x == mp.inf:
        return 'inf'
    return mp.nstr(x, n)


GENERIC_COMBO = (Fr(1), Fr(7, 19), Fr(11, 13))
SUBCOMBOS = (GENERIC_COMBO, (Fr(1), Fr(5, 17), Fr(3, 11)))
SUBSETS = [('V3',), ('V4a',), ('V4b',), ('V3', 'V4a'), ('V3', 'V4b'), ('V4a', 'V4b'), ('V3', 'V4a', 'V4b')]


def analyse_sector(R1, R2, mu, m, ops, others=None, dps=40, quiet=False):
    """exact analysis of the multiplicity space of mu: restriction of the triple (ops = {'V3','V4a','V4b'} Op2's) to the HW space, exact commutation,
       exact simple-spectrum certificates (discriminant + Sturm), mpmath joint spectrum, self-adjointness for the invariant form, and (others) the other Laurent coefficients."""
    t0 = time.perf_counter()
    hw = HWExact(R1, R2, mu)
    assert hw.m == m, ('HW dimension != multiplicity', mu, hw.m, m)
    res = dict(mu=mu, m=m, nW=hw.nW)
    mats = {}
    viol = 0
    for name in ('V3', 'V4a', 'V4b'):
        R, v = restrict_exact(ops[name], hw)
        mats[name] = R
        viol += v
    res['faithful_viol'] = viol
    names = ['V3', 'V4a', 'V4b']
    comm_bad = sum(0 if mat_fr_zero(mat_fr_comm(mats[a], mats[b])) else 1 for a, b in itertools.combinations(names, 2))
    res['comm_bad'] = comm_bad
    G = hw.gram()
    sa = 0
    for name in names:
        R = mats[name]
        lhs = mat_fr_mul(mat_fr_transpose(R), G)
        rhs = mat_fr_mul(G, R)
        sa += 0 if lhs == rhs else 1
    res['selfadj_bad'] = sa
    # exact certificates
    cp3 = exact_charpoly(mats['V3'])
    d3, nr3 = poly_disc_and_real(cp3)
    res['V3_disc_nonzero'] = (d3 != 0)
    res['V3_real_simple'] = (d3 != 0 and nr3 == m)
    M = [[sum((c * mats[n][i][j] for c, n in zip(GENERIC_COMBO, names)), Fr(0)) for j in range(m)] for i in range(m)]
    cpM = exact_charpoly(M)
    dM, nrM = poly_disc_and_real(cpM)
    res['joint_disc_nonzero'] = (dM != 0)
    res['joint_real_simple'] = (dM != 0 and nrM == m)
    js = joint_spectrum([mats[n] for n in names], GENERIC_COMBO, dps)
    res['js'] = js
    res['mats'] = mats
    res['charpoly_V3'] = cp3
    # number of distinct V3 eigenvalues (clusters at 1e-20) and exact separation certificates for every subset of the triple
    ev3 = sorted(row[0] for row in js['trip'])
    res['V3_distinct'] = 1 + sum(1 for i in range(m - 1) if ev3[i + 1] - ev3[i] > mp.mpf(10) ** (-20))
    sub_cert = {}
    for sub in SUBSETS:
        idx = [names.index(x) for x in sub]
        flags = []
        for combo in SUBCOMBOS:
            Mc = [[sum((combo[a] * mats[names[a]][i][j] for a in idx), Fr(0)) for j in range(m)] for i in range(m)]
            d, nr = poly_disc_and_real(exact_charpoly(Mc))
            flags.append(d != 0)
        sub_cert[sub] = flags
    res['sub_cert'] = sub_cert
    if others:
        offs = []
        bad = 0
        for key, op in others.items():
            R, v = restrict_exact(op, hw)
            viol += v
            if not mat_fr_zero(mat_fr_comm(R, M)):
                bad += 1
            offs.append((key, R))
        res['others_noncomm'] = bad
        res['others_faithful_viol'] = viol - res['faithful_viol']
        # diagonal in the joint eigenbasis (mpmath), worst off-diagonal entry over all other Laurent coefficients
        old = mp.mp.dps
        mp.mp.dps = dps
        try:
            Mm = fr_mp(M)
            E, U = mp.eig(Mm)
            order = sorted(range(m), key=lambda k: mp.re(E[k]))
            U = mp.matrix([[U[i, k] for k in order] for i in range(m)])
            Uinv = mp.inverse(U)
            worst = mp.mpf(0)
            for key, R in offs:
                D = Uinv * fr_mp(R) * U
                for i in range(m):
                    for j in range(m):
                        if i != j:
                            worst = max(worst, abs(D[i, j]))
            res['others_offdiag'] = worst
        finally:
            mp.mp.dps = old
        res['others'] = dict(offs)
    res['hw'] = hw
    res['time'] = time.perf_counter() - t0
    return res


def sector_line(res):
    js = res['js']
    return ('  mu=%-14s m=%d nW=%-4d | faithful viol %d, triple non-commuting pairs %d, self-adjoint viol %d | V3 alone: exact disc!=0 %s, real&simple (Sturm) %s, min gap %s | '
            'joint (generic combo): disc!=0 %s, real&simple %s, min joint distance %s, joint offdiag %s, max imag %s%s  (%.1fs)' %
            (str(tuple(res['mu'])), res['m'], res['nW'], res['faithful_viol'], res['comm_bad'], res['selfadj_bad'], res['V3_disc_nonzero'], res['V3_real_simple'], fmt_mp(js['gaps'][0], 6),
             res['joint_disc_nonzero'], res['joint_real_simple'], fmt_mp(js['jgap'], 6), mp.nstr(js['offmax'], 2), mp.nstr(js['imag'], 2),
             (' | other Laurent coeffs: non-commuting with combo %d, worst off-diagonal in joint eigenbasis %s' % (res['others_noncomm'], mp.nstr(res['others_offdiag'], 2))) if 'others_noncomm' in res else '',
             res['time']))


# ====================================================================================================== 8. stages
def hdr(s):
    print('\n' + '=' * 150 + '\n' + s + '\n' + '=' * 150, flush=True)


def lam_str(lam):
    return '[' + ','.join(str(x) for x in lam) + ']'


def raw_word_count(Sk):
    return sum(2 ** len(w) for w in Sk)


def classical_laurent_check(N, seed=1, t=Fr(31, 35)):
    """stage-B check independent of representations: replace e^s_ij by commuting random integers E^s_ij; S_k(x) from the stage-A polynomials with l_ij^{(d)} = d-th derivative of
       E1_ij/x + E2_ij/(x+t) (sympy rational functions); its Laurent coefficients at x = 0 must equal the term lists (with t^power).  Returns the number of mismatching coefficients."""
    import sympy, random
    rng = random.Random(seed)
    E1 = [[rng.randint(-3, 3) for _ in range(N)] for _ in range(N)]
    E2 = [[rng.randint(-3, 3) for _ in range(N)] for _ in range(N)]
    x = sympy.Symbol('x')
    tt = sympy.Rational(t.numerator, t.denominator)
    S = cdet_S(N)
    LT = laurent(N)
    bad = 0
    ncheck = 0
    cache = {}

    def lder(i, j, d):
        key = (i, j, d)
        if key not in cache:
            f = sympy.Integer(E1[i - 1][j - 1]) / x + sympy.Integer(E2[i - 1][j - 1]) / (x + tt)
            cache[key] = sympy.diff(f, x, d)
        return cache[key]
    for k in range(1, N + 1):
        expr = 0
        for w, c in S[k].items():
            term = sympy.Integer(c)
            for (i, j, d) in w:
                term = term * lder(i, j, d)
            expr += term
        expr = sympy.together(expr)
        ser = sympy.series(expr, x, 0, 1).removeO()
        ser = sympy.expand(ser)
        for j in range(0, k + 1):
            cj = ser.coeff(x, -j)
            pred = sympy.Integer(0)
            for (w1, w2), pw in LT[k].get(j, {}).items():
                val = sympy.Integer(1)
                for (a, b) in w1:
                    val *= E1[a - 1][b - 1]
                for (a, b) in w2:
                    val *= E2[a - 1][b - 1]
                for p, cc in pw.items():
                    pred += cc * val * tt ** p
            ncheck += 1
            if sympy.simplify(cj - pred) != 0:
                bad += 1
    return bad, ncheck


def stage_sym():
    hdr('G2.1a  CONVENTION TEST: symbolic cdet(d - L) for gl(3) against eq. (4) of the notes (as noncommutative polynomials)')
    S3 = cdet_S(3)
    ref = eq4_gl3()
    for m in (1, 2, 3):
        print('S_%d  (engine, %d words): %s' % (m, len(S3[m]), fmt_poly_l(S3[m])))
        print('S_%d  (eq. (4), %d words): %s' % (m, len(ref[m]), fmt_poly_l(ref[m])))
        print('   equal as noncommutative polynomials: %s' % (S3[m] == ref[m]))
        assert S3[m] == ref[m]
    print('CONVENTION TEST PASSED: S_1, S_2, S_3 of the engine = eq. (4) exactly (N = 3)')
    hdr('G2.1b  the engine for N = 4: S_1..S_4 as nc polynomials in l_ij^(d) (column-ordered cdet, derivative to the right)')
    S4 = cdet_S(4)
    for m in (1, 2, 3, 4):
        print('S_%d: %d words' % (m, len(S4[m])))
    print('S_1 =', fmt_poly_l(S4[1]))
    print('S_2 =', fmt_poly_l(S4[2]))
    print('S_3 =', fmt_poly_l(S4[3]))
    print('S_4 =', fmt_poly_l(S4[4]))
    hdr('Stage B: Laurent expansion at theta_1 (x = u - theta_1, u - theta_2 = x + t):  term counts, homogeneity in t')
    for N in (3, 4):
        S = cdet_S(N)
        LT = laurent(N)
        print('N = %d: raw (word, site-assignment) counts per S_k: %s' % (N, {k: raw_word_count(S[k]) for k in S}))
        for k in range(1, N + 1):
            row = []
            for j in range(k, -1, -1):
                T = LT[k].get(j, {})
                pows = sorted({p for pw in T.values() for p in pw})
                row.append('j=%d: %d terms, t^%s' % (j, len(T), pows))
            print('  S_%d: %s' % (k, '; '.join(row)))
            for j in range(0, k + 1):
                T = LT[k].get(j, {})
                pows = {p for pw in T.values() for p in pw}
                assert pows <= {j - k}, ('S_{k,j} must be homogeneous of degree j-k in t', k, j, pows)
    print('every S_{k,j} is homogeneous: S_{k,j} = t^(j-k) * (t-independent operator); the canonical triple t*S_{3,2}, t*S_{4,3}, t^2*S_{4,2} is therefore t-independent')
    hdr('Classical-limit check of the Laurent expansion (commuting random integer matrices, sympy series), independent of any representation')
    for N in (3, 4):
        for seed in (1, 2):
            t0 = time.perf_counter()
            bad, nchk = classical_laurent_check(N, seed)
            print('  N=%d seed=%d: %d coefficients (k=1..N, j=0..k) compared, %d mismatches (%.1fs)' % (N, seed, nchk, bad, time.perf_counter() - t0), flush=True)
            assert bad == 0
    hdr('Residue relations (S_k = O(u^-k) at infinity): sum_s sum_{i<=j} C(j-1,i-1) theta_s^(j-i) S_{k,i}^(s) = 0 for j = 1..k-1, exact in U(gl_N) (x) U(gl_N)')
    for N in (3, 4):
        LU = laurent_up(N)
        for k in range(2, N + 1):
            for j in range(1, k):
                tot = {}
                for i in range(1, j + 1):
                    hat, p = LU[(k, i)]
                    up_add(tot, hat, comb(j - 1, i - 1))
                hatj, pj = LU[(k, j)]
                up_add(tot, up_swap(hatj), (-1) ** (j - k))
                res = up_pbw(up_clean(tot))
                print('  N=%d k=%d j=%d:  sum_i C(j-1,i-1) S^_{k,i}(theta_1) + (-1)^(j-k) swap(S^_{k,j}) = 0 in PBW normal form: residual terms %d' % (N, k, j, len(res)))
                assert len(res) == 0
    print('  (they use the theta_2 rule: site swap and t -> -t; all hold identically, so stage B and the swap rule are consistent)')
    hdr('Exact operator identities for N = 3 in U(gl_3) (x) U(gl_3) (PBW normal form)')
    N = 3
    LU = laurent_up(3)
    n1, n2, Om = inv_n(1, N), inv_n(2, N), inv_Omega(N)
    h21, _ = LU[(2, 1)]
    eta = up_lin((1, up_mul(n1, n2)), (-1, Om))
    print('  t*S_{2,1} - (n1 n2 - Omega)  = %s   [H1 = -(1/t) sum E1_ik E2_ki + |L|^2/t, |L|^2 = n1 n2]' % fmt_up(up_pbw(up_lin((1, h21), (-1, eta)))))
    h32, _ = LU[(3, 2)]
    V = antisym_sum_subsets(3, 2, 1)
    rem = up_pbw(up_lin((1, h32), (-1, V), (Fr(-1, 2), eta)))
    print('  t*S_{3,2} - V - (n1 n2 - Omega)/2 = %s   [V = 1/2 sum sgn sgn e1 e1 e2 (36 terms), the notes: t*a_1 = V + eta_1/2, eta_1 = n1 n2 - Omega]' % fmt_up(rem))
    assert len(rem) == 0


def site_scalar(rep, hat, site=1):
    """value of a one-site operator (terms with empty word on the other site) on V_lam if it is a scalar, None otherwise"""
    n = rep.dim
    val = None
    for c in range(n):
        acc = {}
        for (w1, w2), co in hat.items():
            if (w2 if site == 1 else w1):
                continue
            for r, v in rep.word_cols(w1 if site == 1 else w2)[c].items():
                acc[r] = acc.get(r, 0) + co * v
        acc = {r: v for r, v in acc.items() if v != 0}
        if any(r != c for r in acc):
            return None
        v = acc.get(c, Fr(0))
        if val is None:
            val = v
        elif v != val:
            return None
    return val if val is not None else Fr(0)


def newton_Smm(lam, m):
    """closed form of the pole of order m of S_m: S_{m,m} = (-1)^m a_{N-m}, prod_c (th - y_c) = sum_k a_k th(th-1)...(th-k+1), y_c = lam_c + N - c"""
    N = len(lam)
    y = [lam[c] + N - 1 - c for c in range(N)]

    def f(th):
        r = 1
        for yy in y:
            r *= (th - yy)
        return r
    k = N - m
    a = Fr(sum((-1) ** (k - j) * comb(k, j) * f(j) for j in range(k + 1)), factorial(k))
    return (-1) ** m * a


def stage_gl3():
    hdr('G2.1c  gl(3) numerical facts: [2,1,0]^2 (dim 64) and [3,1,0]^2 (dim 225), exact arithmetic (Fractions), all coefficients at theta_1 AND theta_2')
    sys.path.insert(0, os.path.join(HERE, '..', 'gaudin_sl3'))
    import slice_engine as se
    LU = laurent_up(3)
    for lam in [(2, 1, 0), (3, 1, 0)]:
        t0 = time.perf_counter()
        R = Rep.of(lam)
        S = se.Rep.of(lam)
        mism = sum(1 for ij, cols in S.E.items() for c in range(R.dim) if R.E[ij][c] != cols[c]) + sum(1 for a, b in zip(R.g, S.g) if a != b)
        print('\n%s^2  dim %d x %d = %d;  scaffold == slice_engine.Rep.of(lam): %d mismatches (generators + metric);  commutation relations on the scaffold: %d violations' %
              (lam_str(lam), R.dim, R.dim, R.dim ** 2, mism, R.check_commutators()))
        assert mism == 0
        for t in (Fr(31, 35), Fr(-3, 7)):
            ops = build_ops(R, R, 3, t)
            mats = {key: op.matrix_exact() for key, op in ops.items()}
            tot = tot_generators_exact(R, R)
            keys = list(mats)
            nz = [k for k in keys if nnz_exact(mats[k])]
            bad = 0
            npairs = 0
            for a, b in itertools.combinations(keys, 2):
                npairs += 1
                if nnz_exact(comm_exact(mats[a], mats[b])):
                    bad += 1
            badE = 0
            nE = 0
            for k in keys:
                for ij, E in tot.items():
                    nE += 1
                    if nnz_exact(comm_exact(mats[k], E)):
                        badE += 1
            print('  t = %s: %d operators (S_{2,1..2}, S_{3,1..3} at theta_1 and theta_2; nonzero: %d); pairwise commutators: %d pairs, %d violations; with the 9 E^tot_ij: %d tests, %d violations' %
                  (t, len(keys), len(nz), npairs, bad, nE, badE), flush=True)
            assert bad == 0 and badE == 0
        t = Fr(31, 35)
        ops = build_ops(R, R, 3, t)
        mats = {key: op.matrix_exact() for key, op in ops.items()}
        # double pole of S_2
        c2 = lam[1] * (lam[0] + 1)
        sc = site_scalar(R, LU[(2, 2)][0])
        sc2 = site_scalar(R, up_swap(LU[(2, 2)][0]), site=2)
        print('  S_{2,2} (theta_1) = %s * 1 (scalar on V_lam), formula lam_2 (lam_1 + 1) = %d; theta_2 (swapped) = %s' % (sc, c2, sc2))
        assert sc == c2 and sc2 == c2
        # H1
        h21 = {}
        for i in range(1, 4):
            for k in range(1, 4):
                h21[(((i, k),), ((k, i),))] = Fr(-1) / t
        h21[((), ())] = Fr(sum(lam) ** 2) / t
        opH = Op2(h21, R, R)
        MH = opH.matrix_exact()
        diff = nnz_exact(_sub(MH, mats[(2, 1, 1)]))
        print('  H_1 = S_{2,1}(theta_1) vs -(1/t) sum_ik E1_ik E2_ki + |Lambda|^2/t (|Lambda| = %d): %d differing entries (exact, whole space)' % (sum(lam), diff))
        assert diff == 0
        # triple pole
        s33 = site_scalar(R, LU[(3, 3)][0])
        print('  S_{3,3} (theta_1) (the triple pole of S_3) = %s as an operator on V_lam (x) V_lam (scalar: %s); lam_3 = %d%s' %
              (s33, s33 is not None, lam[2], ' -> vanishes identically, not only on HW vectors' if s33 == 0 else ''))
        assert s33 == 0 or lam[2] > 0
        # t S_{3,2} = V + eta_1/2 on every HW vector
        V = {}
        for coef, ij1, ij2, ij3 in se.VTERMS:
            V[((ij1, ij2), (ij3,))] = V.get(((ij1, ij2), (ij3,)), 0) + coef
        opV = Op2(V, R, R)
        op32 = Op2(LU[(3, 2)][0], R, R, scale=1)
        rs = racah_speiser(lam)
        rs_se = se.racah_speiser(lam)
        assert rs == rs_se, 'Racah-Speiser differs from slice_engine'
        nvec = 0
        bad = 0
        for mu, m in sorted(rs.items(), reverse=True):
            hw = HWExact(R, R, mu)
            assert hw.m == m
            eta = Fr(se.eta1(lam, mu))
            for b in hw.basis:
                nvec += 1
                y1 = op32.apply_exact(b)
                y2 = opV.apply_exact(b)
                res = dict(y1)
                for k, v in y2.items():
                    res[k] = res.get(k, 0) - v
                for k, v in b.items():
                    res[k] = res.get(k, 0) - eta / 2 * v
                if any(v != 0 for v in res.values()):
                    bad += 1
        print('  t*S_{3,2} = V + eta_1/2 on HW vectors: %d irreps (Racah-Speiser = slice_engine.racah_speiser), %d exact HW basis vectors (canonical RREF basis of each multiplicity space), %d violations' %
              (len(rs), nvec, bad))
        assert bad == 0
        # whole-space identity
        n1, n2, Om = inv_n(1, 3), inv_n(2, 3), inv_Omega(3)
        eta_up = up_lin((1, up_mul(n1, n2)), (-1, Om))
        full = up_lin((1, V), (Fr(1, 2), eta_up))
        Mfull = Op2(full, R, R).matrix_exact()
        M32 = op32.matrix_exact()
        print('  whole-space operator identity t*S_{3,2} = V + (n1 n2 - Omega)/2: %d differing entries (exact)   [%.1fs]' % (nnz_exact(_sub(M32, Mfull)), time.perf_counter() - t0), flush=True)
        assert nnz_exact(_sub(M32, Mfull)) == 0
    print('\ngl(3) scalar poles for lam_3 != 0 (the triple pole of S_3 is the central scalar S_{3,3} = (lam_1+2)(lam_2+1) lam_3, zero iff lam_3 = 0; S_{2,2} = lam_2(lam_1+1) + lam_3(lam_1+lam_2+2)):')
    for lam in [(2, 1, 1), (3, 2, 1), (4, 2, 2), (2, 2, 2)]:
        R = Rep.of(lam)
        s22 = site_scalar(R, LU[(2, 2)][0])
        s33 = site_scalar(R, LU[(3, 3)][0])
        f22 = lam[1] * (lam[0] + 1) + lam[2] * (lam[0] + lam[1] + 2)
        f33 = (lam[0] + 2) * (lam[1] + 1) * lam[2]
        print('  %s: S_{2,2} = %s (formula %d), S_{3,3} = %s (formula (l1+2)(l2+1)l3 = %d, Newton %s)' % (lam_str(lam), s22, f22, s33, f33, newton_Smm(lam, 3)))
        assert s22 == f22 and s33 == f33 == newton_Smm(lam, 3) and newton_Smm(lam, 2) == f22
    print('\nG2.1 acceptance: exact equalities, 0 violations.')


def comm_products():
    return [((1, 0, 0, 0), (1, 1, 0, 0)), ((1, 1, 0, 0), (2, 1, 0, 0)), ((2, 1, 0, 0), (2, 1, 0, 0)), ((2, 1, 1, 0), (2, 1, 1, 0))]


def float_matfree_comm(R1, R2, t, nvec=2, seed=7):
    """matrix-free float64 test on a product: residual of [S_a, S_b] X and [S_a, E_tot_ij] X on random vectors X (orthonormal two-site components), measured relative to the
       natural scale ||S_a X|| ||S_b X|| / ||X|| (resp. ||S_a X|| ||E X|| / ||X||); operators that vanish numerically (||S X|| < 1e-10 ||X||; identically zero in exact arithmetic)
       are skipped.  Returns (worst pair, worst E, #operators, #nonzero)"""
    rng = np.random.default_rng(seed)
    ops = build_ops(R1, R2, 4, t)
    keys = list(ops)
    Eh1, Eh2 = R1.float_ops(), R2.float_ops()
    worst_pair = 0.0
    worst_E = 0.0
    nz = None
    for _ in range(nvec):
        X = rng.standard_normal((R1.dim, R2.dim))
        nX = np.linalg.norm(X)
        SX = {k: ops[k].apply_float(X) for k in keys}
        nz = [k for k in keys if np.linalg.norm(SX[k]) > 1e-10 * nX]
        for a, b in itertools.combinations(nz, 2):
            ab = ops[a].apply_float(SX[b])
            ba = ops[b].apply_float(SX[a])
            scale = np.linalg.norm(SX[a]) * np.linalg.norm(SX[b]) / nX
            worst_pair = max(worst_pair, np.linalg.norm(ab - ba) / scale)
        for k in nz:
            for ij in itertools.product(range(1, 5), repeat=2):
                E1, E2 = Eh1[ij], Eh2[ij]
                EX = E1 @ X + (E2 @ X.T).T
                a = ops[k].apply_float(EX)
                b = E1 @ SX[k] + (E2 @ SX[k].T).T
                scale = max(np.linalg.norm(SX[k]) * np.linalg.norm(EX) / nX, 1e-300)
                worst_E = max(worst_E, np.linalg.norm(a - b) / scale)
    return worst_pair, worst_E, len(keys), len(nz)


def stage_comm():
    hdr('G2.2  COMMUTATIVITY FOR N = 4: all Laurent coefficients S_{k,j}, k = 2,3,4, j = 1..k, at theta_1 AND theta_2, pairwise and with every E^tot_ij')
    LU = laurent_up(4)
    print('\nConvention variants (information; the binding convention is the column-ordered cdet with the derivative to the right, which passes everything below):')
    for variant, descr in (('col', 'column-ordered, derivative right (binding)'), ('row', 'row-ordered  sum sgn M_{1 s(1)} ... M_{N s(N)}'), ('col-left', 'column-ordered cdet, coefficients read off in the derivative-LEFT normal form')):
        for N, (l1, l2) in ((3, ((2, 1, 0), (2, 1, 0))), (4, ((1, 1, 0, 0), (2, 1, 0, 0))), (4, ((2, 1, 1, 0), (2, 1, 1, 0)))):
            R1, R2 = Rep.of(l1), Rep.of(l2)
            ops = build_ops(R1, R2, N, Fr(31, 35), variant)
            mats = {k: o.matrix_exact() for k, o in ops.items()}
            nz = [k for k in mats if nnz_exact(mats[k])]
            npairs = nbad = 0
            for a, b in itertools.combinations(nz, 2):
                npairs += 1
                nbad += 1 if nnz_exact(comm_exact(mats[a], mats[b])) else 0
            tot = tot_generators_exact(R1, R2)
            nE = sum(1 for k in nz for ij, E in tot.items() if nnz_exact(comm_exact(mats[k], E)))
            print('  %-9s %-62s N=%d %s x %s: %3d nonzero operators, %3d pairs: %3d non-commuting; [S, E^tot] violations: %d' % (variant, descr, N, lam_str(l1), lam_str(l2), len(nz), npairs, nbad, nE))
    print('  => the row-ordered determinant is not even gl(N)-invariant with this Lax matrix (L_ij = E_ij/(u - theta)); the column-ordered one is the right one, the derivative-left re-normal-ordering is the same commutative family.')
    for l1, l2 in comm_products():
        R1, R2 = Rep.of(l1), Rep.of(l2)
        n = R1.dim * R2.dim
        print('\n%s (x) %s : dim %d x %d = %d (scaffold: commutation relations %d violations, metric %d violations)' %
              (lam_str(l1), lam_str(l2), R1.dim, R2.dim, n, R1.check_commutators() + R2.check_commutators(), R1.check_metric() + R2.check_metric()), flush=True)
        for t in (Fr(31, 35), Fr(-3, 7)):
            t0 = time.perf_counter()
            ops = build_ops(R1, R2, 4, t)
            mats = {key: op.matrix_exact() for key, op in ops.items()}
            keys = list(mats)
            nz = [k for k in keys if nnz_exact(mats[k])]
            bad = 0
            npairs = 0
            for a, b in itertools.combinations(nz, 2):
                npairs += 1
                if nnz_exact(comm_exact(mats[a], mats[b])):
                    bad += 1
            tot = tot_generators_exact(R1, R2)
            badE = nE = 0
            for k in nz:
                for ij, E in tot.items():
                    nE += 1
                    if nnz_exact(comm_exact(mats[k], E)):
                        badE += 1
            zero = [k for k in keys if k not in nz]
            print('  t = %-6s exact (Fractions, whole space): %2d operators (18 = 9 at theta_1 + 9 at theta_2), %2d nonzero, identically zero: %s;  %d pairs: %d violations;  with the 16 E^tot_ij: %d tests, %d violations  [%.1fs]' %
                  (t, len(keys), len(nz), [(k[0], k[1], k[2]) for k in zero], npairs, bad, nE, badE, time.perf_counter() - t0), flush=True)
            assert bad == 0 and badE == 0
        if n <= 400:
            ws, wE, nk, nnz_ = float_matfree_comm(R1, R2, Fr(31, 35))
            print('  float64 matrix-free cross-check (2 random vectors, %d nonzero operators): worst relative [S_a,S_b] residual %.2e, worst [S,E^tot] residual %.2e' % (nnz_, ws, wE), flush=True)
    for l1, l2 in [((3, 2, 0, 0), (3, 2, 0, 0)), ((3, 2, 1, 0), (3, 2, 1, 0))]:
        R1, R2 = Rep.of(l1), Rep.of(l2)
        t0 = time.perf_counter()
        ws, wE, nk, nnz_ = float_matfree_comm(R1, R2, Fr(31, 35))
        print('\n%s (x) %s : dim %d: float64 matrix-free test of all %d operators (%d nonzero; 2 random vectors): worst [S_a,S_b] residual %.2e, worst [S,E^tot] residual %.2e (relative to ||S_a X|| ||S_b X||/||X||)  [%.1fs]' %
              (lam_str(l1), lam_str(l2), R1.dim * R2.dim, nk, nnz_, ws, wE, time.perf_counter() - t0), flush=True)
        assert ws < 1e-12 and wE < 1e-12
    hdr('Scalar poles S_{k,k} (central elements of site 1): measured on V_lam vs closed forms')
    print('  S_{2,2} = e_2(lam) + sum_k (k-1) lam_k = lam_2 (lam_1+1) + lam_3 (lam_1+lam_2+2) + lam_4 (lam_1+lam_2+lam_3+3);  general closed form (Capelli / falling factorials):')
    print("  cdet(d - E/x) = x^-N prod_c (theta - y_c), theta = x d_x, y_c = lam_c + N - c  =>  S_{m,m} = (-1)^m a_{N-m},  prod_c(theta - y_c) = sum_k a_k theta(theta-1)...(theta-k+1)")
    reps = [(1, 0, 0, 0), (1, 1, 0, 0), (2, 1, 0, 0), (2, 1, 1, 0), (2, 2, 0, 0), (3, 1, 0, 0), (3, 2, 1, 0), (2, 2, 2, 0), (3, 3, 0, 0), (4, 2, 1, 0), (2, 2, 1, 1), (3, 1, 1, 1), (4, 3, 2, 1)]
    print('  %-12s %6s %6s %6s | %-18s | %-12s %-12s | %-12s %-12s | %s' % ('lam', 'S_{2,2}', 'formula', 'newton', 'lam2(l1+1)+l3(..)+l4(..)', 'S_{3,3}', 'newton', 'S_{4,4}', 'newton', 'S_{4,3}: scalar?'))
    nbad = 0
    for lam in reps:
        R = Rep.of(lam)
        N = 4
        s22 = site_scalar(R, LU[(2, 2)][0])
        s33 = site_scalar(R, LU[(3, 3)][0])
        s44 = site_scalar(R, LU[(4, 4)][0])
        s43 = site_scalar(R, LU[(4, 3)][0])
        form = lam[1] * (lam[0] + 1) + lam[2] * (lam[0] + lam[1] + 2) + lam[3] * (lam[0] + lam[1] + lam[2] + 3)
        e2 = sum(lam[i] * lam[j] for i in range(4) for j in range(i + 1, 4)) + sum(k * lam[k] for k in range(4))
        n22, n33, n44 = newton_Smm(lam, 2), newton_Smm(lam, 3), newton_Smm(lam, 4)
        ok = (s22 == form == e2 == n22) and (s33 == n33) and (s44 == n44)
        nbad += 0 if ok else 1
        # is S_{4,3} (site-1 part only, w2 = ()) a one-site operator?  (it needs a site-2 factor: no pure site-1 terms)
        two_site = all(w2 for (w1, w2) in LU[(4, 3)][0])
        print('  %-12s %6s %6d %6s | %-18d | %-12s %-12s | %-12s %-12s | %s' % (lam_str(lam), s22, form, n22, e2, s33, n33, s44, n44, 'no (every term of S_{4,3} has a site-2 factor: genuine two-site operator V4a/t)' if two_site else s43))
    print('  closed-form mismatches: %d of %d reps' % (nbad, len(reps)))
    assert nbad == 0
    print('\n  Higher poles: S_3 has the triple pole S_{3,3} = (lam_1+2)(lam_2+1) lam_3 [lam_4 = 0] (a central scalar; 0 for lam_3 = 0);  S_4 has S_{4,4} = (lam_1+3)(lam_2+2)(lam_3+1) lam_4 (central scalar, 0 for lam_4 = 0)')
    print('  and S_{4,3} = V4a/t, a genuine two-site operator (no pure one-site part); no pole of order > k occurs (weight counting: sum of (d+1) over the factors of every word = k).')
    print('\nG2.2 acceptance: 0 violations.')


# ---------------------------------------------------------------------------------------------------- G2.3: the canonical triple on multiplicity spaces
TRIPLE_PRODUCTS = [
    ((2, 1, 0, 0), (2, 1, 0, 0)),
    ((2, 1, 1, 0), (2, 1, 1, 0)),
    ((3, 2, 1, 0), (3, 2, 1, 0)),
    ((2, 2, 1, 0), (2, 1, 1, 0)),
    ((3, 2, 0, 0), (2, 1, 1, 0)),
    ((3, 2, 0, 0), (3, 2, 0, 0)),
    ((2, 2, 1, 0), (2, 2, 1, 0)),
]
EXTRA_PRODUCTS = [
    ((4, 2, 1, 0), (4, 2, 1, 0)),
    ((3, 2, 1, 0), (4, 2, 1, 0)),
    ((4, 2, 2, 0), (4, 2, 2, 0)),
    ((4, 3, 1, 0), (4, 3, 1, 0)),
]
SUBSETS = [('V3',), ('V4a',), ('V4b',), ('V3', 'V4a'), ('V3', 'V4b'), ('V4a', 'V4b'), ('V3', 'V4a', 'V4b')]


def min_joint_dist(trip, idx):
    """min over pairs of eigenvectors of the max-norm distance of their joint-eigenvalue tuples restricted to the operator indices idx"""
    m = len(trip)
    best = mp.inf
    for i in range(m):
        for j in range(i + 1, m):
            d = max(abs(trip[i][a] - trip[j][a]) for a in idx)
            if d < best:
                best = d
    return best


def multiplicity_table(l1, l2, rs):
    rows = []
    for mu, m in sorted(rs.items(), reverse=True):
        rows.append((mu, m, weyl_dim(mu), nW_of(l1, l2, mu)))
    return rows


def triple_product(l1, l2, others=True, verbose=True, t=Fr(31, 35), maxnW=None):
    """analyse every sector with m >= 2 of V_l1 (x) V_l2; returns a summary dict"""
    t0 = time.perf_counter()
    R1, R2 = Rep.of(l1), Rep.of(l2)
    rs = racah_speiser(l1, l2)
    assert sum_dim_check(l1, l2, rs)
    ops = triple_ops(R1, R2)
    oth = build_ops(R1, R2, 4, t) if others else None
    if verbose:
        print('\n%s (x) %s : dim %d x %d = %d;  %d irreps, sum of multiplicities %d, max m %d (sum m dim V_mu = dim x dim: True);  Racah-Speiser' %
              (lam_str(l1), lam_str(l2), R1.dim, R2.dim, R1.dim * R2.dim, len(rs), sum(rs.values()), max(rs.values())))
        print('   multiplicity table (mu : m, dim V_mu, nW):  ' + '  '.join('%s:%d,%d,%d' % (''.join(str(x) for x in mu) if max(mu) < 10 else str(mu), m, d, nw) for mu, m, d, nw in multiplicity_table(l1, l2, rs)))
        print('   t-independence: V3 = t*S_{3,2} built from the term list at t = 31/35 and t = -3/7 and scaled by t: ', end='')
        LTm = laurent(4)[3][2]
        ok = True
        for tv in (Fr(31, 35), Fr(-3, 7)):
            up = lt_at(LTm, tv)
            ok = ok and all(v * tv == ops['V3'].up.get(k, 0) for k, v in up.items()) and len(up) == len(ops['V3'].up)
        print('identical to the t-free operator: %s' % ok)
        assert ok
    rows = []
    for mu, m in sorted(rs.items(), reverse=True):
        if m < 2:
            continue
        if maxnW is not None and nW_of(l1, l2, mu) > maxnW:
            if verbose:
                print('  mu=%s m=%d nW=%d skipped (nW > %d)' % (mu, m, nW_of(l1, l2, mu), maxnW))
            continue
        res = analyse_sector(R1, R2, mu, m, ops, others=oth)
        rows.append(res)
        if verbose:
            print(sector_line(res), flush=True)
            for i, row in enumerate(res['js']['trip']):
                print('        eigen-triple (V3, V4a, V4b) #%d: (%s)' % (i + 1, ', '.join(mp.nstr(x, 18) for x in row)))
    summ = dict(l1=l1, l2=l2, rows=rows, nsectors=len(rows), time=time.perf_counter() - t0)
    summ['all_comm'] = all(r['comm_bad'] == 0 for r in rows)
    summ['all_faithful'] = all(r['faithful_viol'] == 0 for r in rows)
    summ['all_selfadj'] = all(r['selfadj_bad'] == 0 for r in rows)
    summ['V3_simple_all'] = all(r['V3_real_simple'] for r in rows)
    summ['joint_simple_all'] = all(r['joint_real_simple'] for r in rows)
    summ['V3_degenerate'] = [(r['mu'], r['m']) for r in rows if not r['V3_real_simple']]
    summ['min_V3_gap'] = min((r['js']['gaps'][0] for r in rows if r['V3_real_simple']), default=None)
    summ['min_joint'] = min((r['js']['jgap'] for r in rows), default=None)
    subs = {}
    for sub in SUBSETS:
        idx = [('V3', 'V4a', 'V4b').index(x) for x in sub]
        # exact certificate: a generic rational combination of the subset has nonzero discriminant on EVERY sector; failing sectors = both combinations have disc = 0
        failing = [(r['mu'], r['m']) for r in rows if not any(r['sub_cert'][sub])]
        subs[sub] = (len(failing) == 0, min((min_joint_dist(r['js']['trip'], idx) for r in rows), default=mp.inf), failing)
    summ['subsets'] = subs
    summ['V3_degenerate_detail'] = [(r['mu'], r['m'], r['V3_distinct']) for r in rows if not r['V3_real_simple']]
    if oth is not None:
        summ['others_noncomm'] = sum(r['others_noncomm'] for r in rows)
        summ['others_offdiag'] = max((r['others_offdiag'] for r in rows), default=mp.mpf(0))
    if verbose:
        print('  SUMMARY %s (x) %s : %d sectors with m >= 2 (sum m = %d); HW-space faithfulness violations %d, non-commuting triple pairs %d, non-self-adjoint operators %d;' %
              (lam_str(l1), lam_str(l2), len(rows), sum(r['m'] for r in rows), sum(r['faithful_viol'] for r in rows), sum(r['comm_bad'] for r in rows), sum(r['selfadj_bad'] for r in rows)))
        print('     V3 alone has real simple spectrum on every sector (exact: discriminant != 0 and Sturm count = m): %s; min gap (simple sectors) %s; V3-degenerate sectors (mu, m, #distinct V3 eigenvalues): %s' %
              (summ['V3_simple_all'], fmt_mp(summ['min_V3_gap'], 6) if summ['min_V3_gap'] is not None else '-', summ['V3_degenerate_detail'] if summ['V3_degenerate_detail'] else 'none'))
        print('     joint spectrum of the triple simple on every sector (exact certificate for a generic rational combination): %s; min joint distance %s' %
              (summ['joint_simple_all'], fmt_mp(summ['min_joint'], 6) if summ['min_joint'] is not None else '-'))
        print('     which subsets of {V3, V4a, V4b} separate every sector (exact: a generic rational combination has nonzero discriminant everywhere; in brackets the min max-norm distance of the joint eigenvalue tuples):  ' +
              '; '.join('{%s}: %s (%s)' % (','.join(s), 'yes' if v[0] else 'NO', fmt_mp(v[1], 5)) for s, v in subs.items()))
        for s, v in subs.items():
            if not v[0]:
                print('        {%s} fails (disc = 0 for two independent generic combinations) on %d sectors: %s' % (','.join(s), len(v[2]), v[2] if len(v[2]) <= 8 else v[2][:8] + ['...']))
        if oth is not None:
            print('     the other 15 Laurent coefficients (theta_1 and theta_2, S_2..S_4 with j = 1..k, 18 operators in all incl. the triple\'s own): non-commuting with the generic combination %d (exact); worst off-diagonal entry in the joint eigenbasis %s' %
                  (summ['others_noncomm'], mp.nstr(summ['others_offdiag'], 3)))
        print('     [%.1fs]' % summ['time'], flush=True)
    return summ


def print_triple_summary(sums, title):
    hdr(title)
    print('%-26s %6s %8s %5s %5s | %-9s %-9s %-9s %-9s | %-12s %-14s | %s' % ('product', 'dim', '#sect m>=2', 'sum m', 'maxm', 'comm', 'self-adj', 'V3 simple', 'joint simple', 'min V3 gap', 'min joint dist', 'V3-degenerate sectors (mu, m)'))
    for s in sums:
        l1, l2 = s['l1'], s['l2']
        rows = s['rows']
        print('%-26s %6d %8d %5d %5d | %-9s %-9s %-9s %-9s | %-12s %-14s | %s' % (lam_str(l1) + 'x' + lam_str(l2), weyl_dim(l1) * weyl_dim(l2), len(rows), sum(r['m'] for r in rows), max((r['m'] for r in rows), default=1),
              'yes' if s['all_comm'] else 'NO', 'yes' if s['all_selfadj'] else 'NO', 'yes' if s['V3_simple_all'] else 'NO', 'yes' if s['joint_simple_all'] else 'NO',
              fmt_mp(s['min_V3_gap'], 5) if s['min_V3_gap'] is not None else '-', fmt_mp(s['min_joint'], 5) if s['min_joint'] is not None else '-',
              ('%d of %d sectors: %s' % (len(s['V3_degenerate']), len(rows), s['V3_degenerate'] if len(s['V3_degenerate']) <= 4 else s['V3_degenerate'][:4] + ['...'])) if s['V3_degenerate'] else 'none'))
    print('\nsubsets of the triple that separate every sector (exact certificates; the number is the min max-norm distance of the joint eigenvalue tuples):')
    print('%-26s | %s' % ('product', ' | '.join('%-16s' % ('{' + ','.join(x) + '}') for x in SUBSETS)))
    for s in sums:
        print('%-26s | %s' % (lam_str(s['l1']) + 'x' + lam_str(s['l2']), ' | '.join('%-16s' % (('yes ' if s['subsets'][x][0] else 'NO  ') + fmt_mp(s['subsets'][x][1], 4)) for x in SUBSETS)))
    print('\ntotals: %d products, %d sectors with m >= 2, sum of multiplicities %d, joint spectrum simple on every sector: %s, V3 alone simple on every sector: %s' %
          (len(sums), sum(len(s['rows']) for s in sums), sum(sum(r['m'] for r in s['rows']) for s in sums), all(s['joint_simple_all'] for s in sums), all(s['V3_simple_all'] for s in sums)))


def stage_triple(extras=False):
    hdr('G2.3  THE CANONICAL TRIPLE V3 = t S_{3,2}, V4a = t S_{4,3}, V4b = t^2 S_{4,2} on the multiplicity spaces (exact: Fractions; certificates: exact discriminants + Sturm; eigenvalues: mpmath 40 digits)')
    print('HW spaces: exact null space (RREF over Q, canonical basis: 1 at a free pair, 0 at the other free pairs) of the total simple raising operators on the weight-mu subspace in the plain GT basis;')
    print('restriction: O b_j = sum_i R_ij b_i verified entry by entry (faithful); invariant form G (x) G (diagonal GT metric): O self-adjoint <=> R^T Gram = Gram R (exact).')
    print('generic combination for the joint certificate: M = V3 + (7/19) V4a + (11/13) V4b;  min gaps are in the eigenvalue units of the operators themselves.')
    sums = []
    for l1, l2 in TRIPLE_PRODUCTS + (EXTRA_PRODUCTS if extras else []):
        sums.append(triple_product(l1, l2))
    print_triple_summary(sums, 'G2.3 SUMMARY TABLE (per product)')
    hdr('G2.3 (reference values for WP-G5): the triple on the multiplicity-one sectors of [3,2,1,0]^2 (exact rationals; the m >= 2 sectors are listed above)')
    l1 = l2 = (3, 2, 1, 0)
    R1 = R2 = Rep.of(l1)
    ops = triple_ops(R1, R2)
    rs = racah_speiser(l1, l2)
    for mu, m in sorted(rs.items(), reverse=True):
        if m != 1:
            continue
        hw = HWExact(R1, R2, mu)
        vals = [restrict_exact(ops[n], hw)[0][0][0] for n in ('V3', 'V4a', 'V4b')]
        print('  mu = %-14s dim V_mu = %4d  nW = %3d :  (V3, V4a, V4b) = (%s, %s, %s)' % (str(mu), weyl_dim(mu), hw.nW, vals[0], vals[1], vals[2]))
    return sums


# ---------------------------------------------------------------------------------------------------- G2.4: structure of the triple
def necklaces(r):
    seen, out = set(), []
    for w in itertools.product((1, 2), repeat=r):
        c = min(w[i:] + w[:i] for i in range(r))
        if c not in seen:
            seen.add(c)
            out.append(c)
    return out


def tname(w):
    """display name of the trace generator T^{(w)}"""
    s = ''.join(str(x) for x in w)
    return {'1': 'n1', '2': 'n2', '11': 'C2_1', '22': 'C2_2', '12': 'Omega', '111': 'C3_1', '222': 'C3_2', '1111': 'C4_1', '2222': 'C4_2'}.get(s, 'T(' + ','.join(str(x) for x in w) + ')')


class TraceLib:
    """generators of the invariant algebra of the diagonal gl(N) in U(gl_N)(x)U(gl_N) up to degree 4 (polarised traces over the necklaces of {1,2}) and monomial libraries"""
    def __init__(self, N=4, maxlen=4):
        self.N = N
        self.gens = []                        # (display name, poly, bidegree, weight, is_scalar)
        for r in range(1, maxlen + 1):
            for w in necklaces(r):
                P = inv_T(w, N)
                scalar = (len(set(w)) == 1) or (w == (1, 2))
                self.gens.append((tname(w), P, (w.count(1), w.count(2)), r, scalar))
        self._pbw = {}

    def library(self, maxd1, maxd2, maxw, nonscalar=True, extra=None):
        """all monomials (multisets of generators) with bidegree <= (maxd1, maxd2) and weight <= maxw; entries (factors tuple of names, poly, bidegree, n_nonscalar)"""
        out = [((), inv_one(), (0, 0), 0)]
        gens = [g for g in self.gens if nonscalar or g[4]]

        def rec(start, names, poly, bd, w, ns):
            for i in range(start, len(gens)):
                nm, P, (a, b), r, sc = gens[i]
                nd = (bd[0] + a, bd[1] + b)
                nw = w + r
                if nd[0] > maxd1 or nd[1] > maxd2 or nw > maxw:
                    continue
                pp = up_mul(poly, P) if poly is not None else P
                nn = names + (nm,)
                nns = ns + (0 if sc else 1)
                out.append((nn, pp, nd, nns))
                rec(i, nn, pp, nd, nw, nns)
        rec(0, (), None, (0, 0), 0, 0)
        if extra:
            out = extra + out
        return out

    def pbw(self, key, poly):
        if key not in self._pbw:
            self._pbw[key] = up_pbw(poly)
        return self._pbw[key]

    def fit(self, target, lib):
        cols = [self.pbw(e[0] if isinstance(e[0], tuple) else (e[0],), e[1]) for e in lib]
        sol, rank, inc = solve_exact(cols, up_pbw(target))
        return sol, rank, inc


def fmt_mono(factors):
    if not factors:
        return '1'
    c = Counter(factors)
    return '*'.join(('%s^%d' % (k, v) if v > 1 else k) for k, v in sorted(c.items()))


def fmt_fit(lib, sol, scalar_names=None):
    """split the solution into the non-scalar part (grouped by the non-scalar generator) and the scalar part"""
    ns, sc = defaultdict(list), []
    for e, c in zip(lib, sol):
        if c == 0:
            continue
        names = e[0]
        nsn = [n for n in names if n.startswith('T(')]
        rest = [n for n in names if not n.startswith('T(')]
        if nsn:
            assert len(nsn) == 1
            ns[nsn[0]].append((c, rest))
        else:
            sc.append((c, rest))

    def fpoly(terms):
        return ' '.join('%s%s*%s' % ('+' if c > 0 else '-', abs(c), fmt_mono(r)) for c, r in terms) if terms else '0'
    parts = []
    for g_, terms in ns.items():
        parts.append('[%s] * %s' % (fpoly(terms), g_))
    return parts, fpoly(sc)


def eval_scalar_factors(names, ctx):
    v = Fr(1)
    for n in names:
        v *= ctx[n]
    return v


def sector_ctx(l1, l2, mu, cache={}):
    """values of the scalar generators on the sector (l1, l2, mu): n1 = |l1|, n2 = |l2|, C2_s = casimir2, Omega = (C2(mu) - C2(l1) - C2(l2))/2, C3_s (central cubic, evaluated on the HW vector of V_l)"""
    N = len(mu)

    def c3(lam, site):
        key = (lam, site)
        if key not in cache:
            P = inv_T((1, 1, 1) if site == 1 else (2, 2, 2), N)
            cache[key] = site_scalar(Rep.of(lam), P, site=site)
        return cache[key]
    return {'n1': Fr(sum(l1)), 'n2': Fr(sum(l2)), 'C2_1': Fr(casimir2(l1)), 'C2_2': Fr(casimir2(l2)),
            'Omega': Fr(casimir2(mu) - casimir2(l1) - casimir2(l2), 2), 'C3_1': c3(tuple(l1), 1), 'C3_2': c3(tuple(l2), 2)}


def mat_fr_sub_scalar(R, c):
    m = len(R)
    return [[R[i][j] - (c if i == j else 0) for j in range(m)] for i in range(m)]


def mat_fr_lin(terms, m):
    """sum of c*M"""
    out = [[Fr(0)] * m for _ in range(m)]
    for c, M in terms:
        for i in range(m):
            for j in range(m):
                out[i][j] += c * M[i][j]
    return out


def is_scalar_matrix(R):
    m = len(R)
    c = R[0][0]
    return all(R[i][j] == (c if i == j else 0) for i in range(m) for j in range(m)), c


def span_theorem(verbose=True):
    """PBW-level exact statement: every Laurent coefficient S_{k,j} (k = 2..4, j = 1..k, theta_1 and theta_2) is an exact linear combination of V3, V4a, V4b and a polynomial in the
       scalar invariants n_s, C2_s, C3_s, C4_s (Casimirs of the two sites), Omega and C3tot, C4tot (central elements of the total gl(4); scalars on every sector).  Returns the coefficient table."""
    N = 4
    LU = laurent_up(4)

    def Ctot(r):
        P = {}
        for s in itertools.product((1, 2), repeat=r):
            up_add(P, inv_T(s, N))
        return P
    gens = [('n1', inv_n(1, N), 1), ('n2', inv_n(2, N), 1), ('C2_1', inv_C2(1, N), 2), ('C2_2', inv_C2(2, N), 2), ('Omega', inv_Omega(N), 2),
            ('C3_1', inv_T((1, 1, 1), N), 3), ('C3_2', inv_T((2, 2, 2), N), 3), ('C3tot', Ctot(3), 3),
            ('C4_1', inv_T((1, 1, 1, 1), N), 4), ('C4_2', inv_T((2, 2, 2, 2), N), 4), ('C4tot', Ctot(4), 4)]
    lib = [((), inv_one())]

    def rec(start, names, poly, w):
        for i in range(start, len(gens)):
            nm, P, r = gens[i]
            if w + r > 4:
                continue
            pp = up_mul(poly, P) if poly is not None else P
            lib.append((names + (nm,), pp))
            rec(i, names + (nm,), pp, w + r)
    rec(0, (), None, 0)
    base = [(('V3',), LU[(3, 2)][0]), (('V4a',), LU[(4, 3)][0]), (('V4b',), LU[(4, 2)][0])]
    cols = [up_pbw(P) for _, P in base + lib]
    table = {}
    for (k, j) in [(2, 1), (3, 1), (3, 2), (4, 1), (4, 2), (4, 3)]:
        hat, p = LU[(k, j)]
        for side in (1, 2):
            P = hat if side == 1 else up_swap(hat)
            sol, rank, inc = solve_exact(cols, up_pbw(P))
            assert sol is not None, ('S_{%d,%d} at theta_%d is not in the span' % (k, j, side))
            table[(k, j, side)] = (sol[0], sol[1], sol[2], p)
    if verbose:
        print('  library: V3, V4a, V4b and %d monomials (weight <= 4) in the scalar invariants n_s, C2_s, C3_s, C4_s, Omega, C3tot, C4tot  (all independent: rank %d)' % (len(lib), rank))
        print('  every S_{k,j} at theta_1 and theta_2 (hat = the t-free operator: S^{(1)}_{k,j} = t^(j-k) hat, S^{(2)}_{k,j} = (-t)^(j-k) swap(hat)) lies in this span; coefficients (V3, V4a, V4b) of hat / swap(hat):')
        for (k, j, side), (a, b, c, p) in table.items():
            print('     S_{%d,%d} at theta_%d : (%s, %s, %s) * (V3, V4a, V4b) + [polynomial in the scalar invariants]' % (k, j, side, a, b, c))
    # the polarised traces over the module generated by the triple with coefficients in the scalars n_1, n_2 (V3 times n_s appears at degree 4)
    n1_, n2_ = inv_n(1, N), inv_n(2, N)
    base2 = base + [(('n1*V3',), up_mul(n1_, LU[(3, 2)][0])), (('n2*V3',), up_mul(n2_, LU[(3, 2)][0]))]
    cols2 = [up_pbw(P) for _, P in base2 + lib]
    names2 = [n[0] if n else '1' for n, _ in base2]
    Q = up_lin((1, inv_T((1, 1, 2, 2), N)), (Fr(1, 2), inv_T((1, 2, 1, 2), N)))
    expected = {'T(1,2,2)': {'V3': -1}, 'T(1,1,1,2)': {'V3': -3, 'V4a': -1, 'n1*V3': 1}, 'T(1,2,2,2)': {'V3': 3, 'V4a': 2, 'V4b': 1, 'n2*V3': -1},
                'T(1,1,2,2) + 1/2 T(1,2,1,2)': {'V4a': -1, 'V4b': -1, 'n1*V3': -1, 'n2*V3': 1}}
    targets = {'T(1,2,2)': inv_T((1, 2, 2), N), 'T(1,1,1,2)': inv_T((1, 1, 1, 2), N), 'T(1,2,2,2)': inv_T((1, 2, 2, 2), N), 'T(1,1,2,2) + 1/2 T(1,2,1,2)': Q,
               'T(1,1,2,2)': inv_T((1, 1, 2, 2), N), 'T(1,2,1,2)': inv_T((1, 2, 1, 2), N)}
    tab2 = {}
    for nm, P in targets.items():
        sol, rank2, inc = solve_exact(cols2, up_pbw(P))
        if sol is None:
            tab2[nm] = None
        else:
            tab2[nm] = {names2[i]: x for i, x in enumerate(sol[:len(base2)]) if x != 0}
    if verbose:
        print('  polarised traces modulo scalars, in the module spanned by V3, V4a, V4b, n1 V3, n2 V3 (scalar coefficients n_s):')
        for nm, v in tab2.items():
            print('     %-30s %s' % (nm, ('= ' + ' '.join('%s%s*%s' % ('+' if x > 0 else '-', abs(x), k) for k, x in v.items()) + ' + [scalar polynomial]') if v is not None else 'NOT in the span (not part of the commutative family)'))
    for nm, exp in expected.items():
        assert tab2[nm] is not None and {k: Fr(x) for k, x in exp.items()} == tab2[nm], (nm, tab2[nm])
    assert tab2['T(1,1,2,2)'] is None and tab2['T(1,2,1,2)'] is None
    return table


def naive_trace_commutation(verbose=True):
    """which naive polarised traces commute with the triple on the HW spaces (exact), on the sectors with m >= 2 of [3,2,1,0]^2 and [4,2,1,0]^2"""
    N = 4
    Tw = {'T(1,1,2)': (1, 1, 2), 'T(1,2,2)': (1, 2, 2), 'T(1,1,1,2)': (1, 1, 1, 2), 'T(1,1,2,2)': (1, 1, 2, 2), 'T(1,2,1,2)': (1, 2, 1, 2), 'T(1,2,2,2)': (1, 2, 2, 2)}
    Tpoly = {k: inv_T(w, N) for k, w in Tw.items()}
    Tpoly['T(1,1,2,2) + 1/2 T(1,2,1,2)'] = up_lin((1, Tpoly['T(1,1,2,2)']), (Fr(1, 2), Tpoly['T(1,2,1,2)']))
    cnt = {k: {'V3': [0, 0], 'V4a': [0, 0], 'V4b': [0, 0]} for k in Tpoly}
    nsect = 0
    for l1, l2 in [((3, 2, 1, 0), (3, 2, 1, 0)), ((4, 2, 1, 0), (4, 2, 1, 0))]:
        R1, R2 = Rep.of(l1), Rep.of(l2)
        ops = triple_ops(R1, R2)
        tops = {k: Op2(P, R1, R2) for k, P in Tpoly.items()}
        rs = racah_speiser(l1, l2)
        for mu, m in sorted(rs.items()):
            if m < 2:
                continue
            hw = HWExact(R1, R2, mu)
            nsect += 1
            Rv = {n: restrict_exact(o, hw)[0] for n, o in ops.items()}
            for k, op in tops.items():
                RT = restrict_exact(op, hw)[0]
                for n in ('V3', 'V4a', 'V4b'):
                    cnt[k][n][0] += 1
                    cnt[k][n][1] += 0 if mat_fr_zero(mat_fr_comm(RT, Rv[n])) else 1
    if verbose:
        print('  %d sectors with m >= 2 ([3,2,1,0]^2: 8, [4,2,1,0]^2: 14); entries = number of sectors where the commutator [T, V] is NOT zero (exact):' % nsect)
        print('     %-32s %6s %6s %6s' % ('trace', '[.,V3]', '[.,V4a]', '[.,V4b]'))
        for k in cnt:
            print('     %-32s %6d %6d %6d' % (k, cnt[k]['V3'][1], cnt[k]['V4a'][1], cnt[k]['V4b'][1]))
    return cnt


def cubic_general_N(Ns=(3, 4, 5, 6), verbose=True):
    """V3 = t S_{3,2} = sum_I D(2,1)_I + c_N eta, eta = n1 n2 - Omega, for gl(N): the exact coefficient is c_N = (N-2)/2 (N = 3, 4, 5, 6 solved exactly in PBW normal form)"""
    out = {}
    for N in Ns:
        LU = laurent_up(N)
        V3 = LU[(3, 2)][0]
        n1, n2, Om = inv_n(1, N), inv_n(2, N), inv_Omega(N)
        eta = up_lin((1, up_mul(n1, n2)), (-1, Om))
        D21 = antisym_sum_subsets(N, 2, 1)
        sol, rank, inc = solve_exact([up_pbw(D21), up_pbw(eta)], up_pbw(V3))
        assert sol is not None and sol[0] == 1 and sol[1] == Fr(N - 2, 2), (N, sol)
        out[N] = (sol[1], len(V3))
    if verbose:
        print('  general N (exact, PBW): V3 = t S_{3,2} = sum_{|I|=3} D(2,1)_I + c_N (n1 n2 - Omega) with c_N = %s for N = %s, i.e. c_N = (N-2)/2 (terms of V3: %s)' %
              (', '.join(str(out[N][0]) for N in Ns), ', '.join(str(N) for N in Ns), ', '.join(str(out[N][1]) for N in Ns)))
    return out


def exchange_parity(verbose=True):
    """site exchange P (nu, q) -> (q, nu) on V_lam (x) V_lam: P swap(O) P = O.  On every sector: P V3 P = -V3 + scalar, P V4a P = -2 V4a - V4b + scalar, P V4b P = 3 V4a + 2 V4b + scalar (exact);
       V3 - s/2 is P-odd, so on M_mu (p+ symmetric, p- antisymmetric copies) its spectrum is symmetric and has >= |p+ - p-| zero modes: V3 alone cannot separate the copies when |p+ - p-| >= 2."""
    lams = [(2, 1, 0, 0), (2, 1, 1, 0), (3, 2, 1, 0), (3, 2, 0, 0), (2, 2, 1, 0), (4, 2, 1, 0), (4, 2, 2, 0), (4, 3, 1, 0), (3, 3, 1, 0), (3, 3, 2, 0)]
    tot = defaultdict(int)
    detail = []
    for lam in lams:
        R = Rep.of(lam)
        rs = racah_speiser(lam, lam)
        ops = triple_ops(R, R)
        for mu, m in sorted(rs.items(), reverse=True):
            if m < 2:
                continue
            hw = HWExact(R, R, mu)
            frees = [hw.pairs[f] for f in hw.free]
            cols = []
            for b in hw.basis:
                y = {(q, nu): v for (nu, q), v in b.items()}
                cols.append([y.get(k, Fr(0)) for k in frees])
            RP =[[cols[j][i] for j in range(m)] for i in range(m)]
            # exact check that P b_j = sum_i RP_ij b_i  (P preserves the HW space)
            ok_pres = True
            for j, b in enumerate(hw.basis):
                y = {(q, nu): v for (nu, q), v in b.items()}
                rec = {}
                for i, bi in enumerate(hw.basis):
                    c = RP[i][j]
                    if c != 0:
                        for k, v in bi.items():
                            rec[k] = rec.get(k, 0) + c * v
                if {k: v for k, v in rec.items() if v != 0} != y:
                    ok_pres = False
            RP2 = mat_fr_mul(RP, RP)
            inv = all(RP2[i][j] == (1 if i == j else 0) for i in range(m) for j in range(m))
            trP = sum(RP[i][i] for i in range(m))
            pplus = int((m + trP) / 2)
            pminus = m - pplus
            Rm = {n: restrict_exact(ops[n], hw)[0] for n in ('V3', 'V4a', 'V4b')}

            def conj(M):
                return mat_fr_mul(mat_fr_mul(RP, M), RP)
            c3 = mat_fr_lin([(1, conj(Rm['V3'])), (1, Rm['V3'])], m)
            ca = mat_fr_lin([(1, conj(Rm['V4a'])), (2, Rm['V4a']), (1, Rm['V4b'])], m)
            cb = mat_fr_lin([(1, conj(Rm['V4b'])), (-3, Rm['V4a']), (-2, Rm['V4b'])], m)
            ok_id = is_scalar_matrix(c3)[0] and is_scalar_matrix(ca)[0] and is_scalar_matrix(cb)[0]
            res = analyse_sector(R, R, mu, m, ops)
            deg = res['V3_distinct'] < m
            forced = abs(pplus - pminus) >= 2
            tot['sectors'] += 1
            tot['ok_pres'] += int(ok_pres and inv)
            tot['ok_id'] += int(ok_id)
            tot['forced'] += int(forced)
            tot['forced_and_degenerate'] += int(forced and deg)
            tot['degenerate'] += int(deg)
            if deg:
                detail.append((lam, mu, m, pplus, pminus, res['V3_distinct'], forced))
    if verbose:
        print('  %d sectors with m >= 2 of the products V_lam (x) V_lam, lam in %s' % (tot['sectors'], [lam_str(l) for l in lams]))
        print('  P maps every HW space to itself and P^2 = 1 (exact): %d / %d sectors;  P V3 P = -V3 + c, P V4a P = -2 V4a - V4b + c\', P V4b P = 3 V4a + 2 V4b + c\'\' with scalars c (exact): %d / %d sectors' %
              (tot['ok_pres'], tot['sectors'], tot['ok_id'], tot['sectors']))
        print('  (so on the P-even / P-odd combinations: V4a + V4b is P-even, V3 and 3 V4a + V4b are P-odd, modulo scalars)')
        print('  V3-degenerate sectors: %d;  of them with |p+ - p-| >= 2 (forced by the exchange symmetry): %d of the %d sectors with |p+ - p-| >= 2 (all of them are degenerate);' % (tot['degenerate'], tot['forced_and_degenerate'], tot['forced']))
        print('  the other V3-degenerate sectors have |p+ - p-| <= 1 (an additional, non-forced degeneracy: the off-diagonal block of the odd operator is rank deficient).')
        print('     lam            mu              m  p+ p-  #distinct V3 eigenvalues  forced')
        for lam, mu, m, pp, pm, nd, forced in detail:
            print('     %-14s %-15s %2d  %d  %d   %2d                       %s' % (lam_str(lam), str(mu), m, pp, pm, nd, 'yes' if forced else 'no'))
    assert tot['ok_pres'] == tot['sectors'] and tot['ok_id'] == tot['sectors'] and tot['forced_and_degenerate'] == tot['forced']
    return tot, detail


def stage_struct(extras=False):
    hdr('G2.4  STRUCTURE OF THE TRIPLE: exact identities in U(gl_4) (x) U(gl_4) (PBW normal form, rational coefficients) and their numerical verification on the HW spaces')
    N = 4
    LU = laurent_up(4)
    V3, V4a, V4b = LU[(3, 2)][0], LU[(4, 3)][0], LU[(4, 2)][0]
    TL = TraceLib(4, 4)
    n1, n2, Om = inv_n(1, N), inv_n(2, N), inv_Omega(N)
    eta = up_lin((1, up_mul(n1, n2)), (-1, Om))
    D21, D12 = antisym_sum_subsets(N, 2, 1), antisym_sum_subsets(N, 1, 2)
    D31, D22 = antisym_form(N, 3, 1), antisym_form(N, 2, 2)
    print('Notation: n_s = sum_i e^s_ii, C2_s = sum e^s_ij e^s_ji, C3_s = sum e^s_{ji} e^s_{kj} e^s_{ik}, Omega = sum e^1_ij e^2_ji, eta = n1 n2 - Omega (= t S_{2,1}),')
    print('  T(s1..sr) = sum_{i1..ir} e^{s1}_{i2 i1} e^{s2}_{i3 i2} ... e^{sr}_{i1 ir} (polarised traces, CPV): T(1,1,2) = sum e1_ji e1_kj e2_ik, T(1,2,2) = sum e1_ji e2_kj e2_ik, T(1,1,1,2), T(1,1,2,2), T(1,2,1,2), T(1,2,2,2);')
    print('  D(a,b)_I = 1/(a! b!) sum_{sigma,tau in S(I)} sgn sgn e1_{s1 t1}..e1_{sa ta} e2_{s(a+1) t(a+1)}..e2_{s(a+b) t(a+b)} (first factor acts last); D(a,b) = D(a,b)_{1234}; sum_I D(2,1)_I over the 3-subsets (gl(3)-type vertex operators).')
    print('  Scalars on a sector (lam, lam\', mu): n1 = |lam|, n2 = |lam\'|, C2_s = C2(lam^s), Omega = (C2(mu) - C2(lam) - C2(lam\'))/2, C3_s = value of the central cubic element on V_{lam^s}.')
    hdr('(a) cubic:  V3 = t S_{3,2}')
    sol, rank, inc = TL.fit(V3, [(('D21sum',), D21, (2, 1), 0), (('Omega',), Om, (1, 1), 0), (('n1', 'n2'), up_mul(n1, n2), (1, 1), 0), (('n2',), n2, (0, 1), 0)])
    print('  V3 = sum_I D(2,1)_I + c_Omega*Omega + c_nn*n1*n2 + c_n2*n2  (4 unknowns, rank %d, inconsistent rows %d): (c_Omega, c_nn, c_n2) = (%s, %s, %s)' % (rank, inc, sol[1], sol[2], sol[3]))
    print('     i.e.  V3 = sum_I D(2,1)_I + (n1 n2 - Omega) = sum_I D(2,1)_I + t S_{2,1}  [gl(3): V3 = V + (n1 n2 - Omega)/2 = V + eta_1/2; here the coefficient is (N-2)/2 = 1]   (coefficient of sum_I D(2,1)_I: %s)' % sol[0])
    assert sol == [1, -1, 1, 0]
    delta = {}
    for s_ in itertools.permutations(range(1, 5)):
        for t_ in itertools.permutations(range(1, 5)):
            if s_[3] != t_[3]:
                continue
            k_ = (((s_[0], t_[0]), (s_[1], t_[1])), ((s_[2], t_[2]),))
            delta[k_] = delta.get(k_, 0) + Fr(perm_sign(s_) * perm_sign(t_), 2)
    print('     the brief\'s delta form: 1/2 sum_{sigma,tau in S_4} sgn(sigma) sgn(tau) e1_{s1 t1} e1_{s2 t2} e2_{s3 t3} delta_{s4 t4} == sum_I D(2,1)_I as polynomials (exact, %d terms): %s' % (len(D21), up_clean(delta) == D21))
    assert up_clean(delta) == D21
    cubic_general_N()
    lib3 = TL.library(2, 1, 3)
    sol, rank, inc = TL.fit(V3, lib3)
    parts, scal = fmt_fit(lib3, sol)
    print('  V3 = T(1,1,2) + remainder, remainder = %s   (library of %d invariants of bidegree <= (2,1), rank %d, inconsistent %d)' % (scal, len(lib3), rank, inc))
    print('     non-scalar part: %s' % parts)
    assert inc == 0 and parts == ['[+1*1] * T(1,1,2)']
    C3tot = {}
    for s in itertools.product((1, 2), repeat=3):
        up_add(C3tot, inv_T(s, N))
    C31, C32 = inv_T((1, 1, 1), N), inv_T((2, 2, 2), N)
    T112, T122 = inv_T((1, 1, 2), N), inv_T((1, 2, 2), N)
    extra = [(('C3tot',), C3tot, (3, 3), 0)]
    libz = extra + TL.library(3, 3, 3, nonscalar=False)
    sol, rank, inc = TL.fit(up_lin((1, T112), (1, T122)), libz)
    print('  T(1,1,2) + T(1,2,2) is a combination of Casimirs (C3tot = the cubic central element of the total gl(4): C3tot = sum over the 8 site assignments of T): consistency %s, rank %d' % ('exact' if inc == 0 else 'FAILED', rank))
    nz = [(c, e[0]) for e, c in zip(libz, sol) if c != 0] if sol is not None else None
    print('     T(1,1,2) + T(1,2,2) = ' + ' '.join('%s%s*%s' % ('+' if c > 0 else '-', abs(c), fmt_mono(list(n))) for c, n in nz))
    assert inc == 0
    lib3c = [(('T(1,2,2)',), T122, (1, 2), 1)] + libz
    sol, rank, inc = TL.fit(V3, lib3c)
    parts, scal = fmt_fit(lib3c, sol)
    print('  V3 in terms of T(1,2,2) and Casimirs (C3tot = central cubic of the total gl(4)): rank %d, inconsistent %d:  non-scalar part %s;  scalar part %s' % (rank, inc, parts, scal))
    assert inc == 0
    print('     (on a sector T(1,1,2) + T(1,2,2) is a scalar, so V3 - (a T(1,1,2) + b T(1,2,2)) is a scalar for every (a, b) with a - b = 1, e.g. (1,0) and (0,-1))')
    hdr('(b) quartic pair:  V4a = t S_{4,3},  V4b = t^2 S_{4,2}')
    sol, rank, inc = TL.fit(V4a, [(('D31',), D31, (3, 1), 0), (('D21sum',), D21, (2, 1), 0), (('n1', 'n2'), up_mul(n1, n2), (1, 1), 0), (('Omega',), Om, (1, 1), 0)])
    print('  V4a = D(3,1) + c1 sum_I D(2,1)_I + c2 n1 n2 + c3 Omega: rank %d, inconsistent %d, (c1, c2, c3) = (%s, %s, %s)  =>  V4a = D(3,1) + V3 - eta/3' % (rank, inc, sol[1], sol[2], sol[3]))
    assert sol == [1, 1, Fr(2, 3), Fr(-2, 3)]
    libb = [(('D22',), D22, (2, 2), 0), (('D31',), D31, (3, 1), 0), (('D12sum',), D12, (1, 2), 0), (('D21sum',), D21, (2, 1), 0), (('n1', 'n2'), up_mul(n1, n2), (1, 1), 0), (('Omega',), Om, (1, 1), 0)]
    sol, rank, inc = TL.fit(V4b, libb)
    print('  V4b = D(2,2) - D(3,1) + c1 sum_I D(1,2)_I + c2 sum_I D(2,1)_I + c3 n1 n2 + c4 Omega: rank %d, inconsistent %d, solution (D22, D31, D12sum, D21sum, n1n2, Omega) = %s' % (rank, inc, [str(x) for x in sol]))
    print('     (classically V4a = P(3,1), V4b = P(2,2) - P(3,1); the antisymmetrised forms are the gl(4)-analogues of the gl(3) vertex operator; quantum corrections are the lower antisymmetrised forms and eta.)')
    assert inc == 0
    lib4a = TL.library(3, 1, 4)
    sol, rank, inc = TL.fit(V4a, lib4a)
    parts, scal = fmt_fit(lib4a, sol)
    print('  V4a in the trace basis (library of %d invariants of bidegree <= (3,1), weight <= 4: rank %d, inconsistent %d):' % (len(lib4a), rank, inc))
    print('     non-scalar part: ' + ' + '.join(parts))
    print('     scalar part (a polynomial in n1, n2, C2_1, C3_1, Omega): ' + scal)
    assert inc == 0
    solA = (lib4a, sol)
    lib4b = TL.library(3, 2, 4)
    sol, rank, inc = TL.fit(V4b, lib4b)
    parts, scal = fmt_fit(lib4b, sol)
    print('  V4b in the trace basis (library of %d invariants of bidegree <= (3,2), weight <= 4: rank %d, inconsistent %d):' % (len(lib4b), rank, inc))
    print('     non-scalar part: ' + ' + '.join(parts))
    print('     scalar part (a polynomial in n1, n2, C2_1, C2_2, C3_1, Omega): ' + scal)
    assert inc == 0
    solB = (lib4b, sol)
    hdr('(c) the triple generates the whole family: every Laurent coefficient is a combination of V3, V4a, V4b and Casimir-type scalars (exact, PBW level)')
    span_theorem()
    hdr('(d) naive polarised traces versus the exact Gaudin combinations (commutation on the HW spaces, exact)')
    naive_trace_commutation()
    print('  => T(1,1,2), T(1,2,2), T(1,1,1,2), T(1,2,2,2) are (mod scalars) combinations of the triple (see (c)) and commute with it; T(1,1,2,2) and T(1,2,1,2) do NOT commute individually: only T(1,1,2,2) + 1/2 T(1,2,1,2) belongs to the commutative family (the quantum correction is the relative weight 1/2).')
    hdr('(e) site-exchange symmetry (lam = lam\'): why V3 alone can fail to separate the copies')
    exchange_parity()
    hdr('Numerical verification on the HW spaces (exact): remainders are multiples of the identity, with the predicted scalar')
    prods = [((2, 1, 0, 0), (2, 1, 0, 0)), ((2, 1, 1, 0), (2, 1, 1, 0)), ((3, 2, 1, 0), (3, 2, 1, 0)), ((2, 2, 1, 0), (2, 1, 1, 0)), ((3, 2, 0, 0), (2, 1, 1, 0)), ((3, 2, 0, 0), (3, 2, 0, 0)), ((4, 2, 1, 0), (4, 2, 1, 0))]
    if extras:
        prods += [((4, 3, 1, 0), (4, 3, 1, 0))]
    Tops = {}
    stats = defaultdict(lambda: [0, 0, 0])              # identity -> [sectors, violations, sectors with m >= 2]
    shown = []
    for l1, l2 in prods:
        t0 = time.perf_counter()
        R1, R2 = Rep.of(l1), Rep.of(l2)
        rs = racah_speiser(l1, l2)
        ops = triple_ops(R1, R2)
        nsgens = {}
        for nm, P, bd, r, sc in TL.gens:
            if not sc:
                nsgens[nm] = Op2(P, R1, R2)
        cubic_ops = {'D21sum': Op2(D21, R1, R2), 'D12sum': Op2(D12, R1, R2), 'D31': Op2(D31, R1, R2), 'D22': Op2(D22, R1, R2)}
        for mu, m in sorted(rs.items(), reverse=True):
            if nW_of(l1, l2, mu) > 260:
                continue
            hw = HWExact(R1, R2, mu)
            assert hw.m == m
            ctx = sector_ctx(l1, l2, mu)
            Rm = {}
            for nm, op in ops.items():
                Rm[nm] = restrict_exact(op, hw)[0]
            RT = {nm: restrict_exact(op, hw)[0] for nm, op in nsgens.items()}
            RD = {nm: restrict_exact(op, hw)[0] for nm, op in cubic_ops.items()}
            I = [[Fr(1) if i == j else Fr(0) for j in range(m)] for i in range(m)]
            eta_v = ctx['n1'] * ctx['n2'] - ctx['Omega']

            def check(label, target, lib, sol):
                """remainder = target - sum_{non-scalar} ... must equal the scalar part (value) times the identity"""
                nsm = [[Fr(0)] * m for _ in range(m)]
                sc = Fr(0)
                for e, c in zip(lib, sol):
                    if c == 0:
                        continue
                    names = e[0]
                    nsn = [n for n in names if n.startswith('T(')]
                    rest = [n for n in names if not n.startswith('T(')]
                    val = c * eval_scalar_factors(rest, ctx)
                    if nsn:
                        for i in range(m):
                            for j in range(m):
                                nsm[i][j] += val * RT[nsn[0]][i][j]
                    else:
                        sc += val
                diff = [[target[i][j] - nsm[i][j] - (sc if i == j else 0) for j in range(m)] for i in range(m)]
                ok = mat_fr_zero(diff)
                s = stats[label]
                s[0] += 1
                s[1] += 0 if ok else 1
                s[2] += 1 if m >= 2 else 0
                return ok
            # (1) V3 = sum_I D(2,1)_I + eta
            ok1 = mat_fr_zero(mat_fr_lin([(1, Rm['V3']), (-1, RD['D21sum']), (-eta_v, I)], m))
            stats['V3 = sum_I D(2,1)_I + eta'][0] += 1
            stats['V3 = sum_I D(2,1)_I + eta'][1] += 0 if ok1 else 1
            stats['V3 = sum_I D(2,1)_I + eta'][2] += 1 if m >= 2 else 0
            # (2) V3 = T112 + scalars
            check('V3 = T(1,1,2) + scalar polynomial', Rm['V3'], lib3, TL.fit(V3, lib3)[0])
            # (3) V4a = D31 + V3 - eta/3
            ok3 = mat_fr_zero(mat_fr_lin([(1, Rm['V4a']), (-1, RD['D31']), (-1, Rm['V3']), (eta_v / 3, I)], m))
            stats['V4a = D(3,1) + V3 - eta/3'][0] += 1
            stats['V4a = D(3,1) + V3 - eta/3'][1] += 0 if ok3 else 1
            stats['V4a = D(3,1) + V3 - eta/3'][2] += 1 if m >= 2 else 0
            # (4) V4b = D22 - D31 + (D12 - D21)/2 - eta/6   (solution of the library with a dependency: use the stored one)
            # (5),(6) trace-basis identities
            check('V4a = -T(1,1,1,2) + (n1-3) T(1,1,2) + scalar polynomial', Rm['V4a'], solA[0], solA[1])
            check('V4b = trace-basis non-scalar part + scalar polynomial', Rm['V4b'], solB[0], solB[1])
            # (7) sector scalar T112 + T122
            Rsum = mat_fr_lin([(1, RT['T(1,1,2)']), (1, RT['T(1,2,2)'])], m)
            isc, cval = is_scalar_matrix(Rsum)
            stats['T(1,1,2) + T(1,2,2) is a multiple of the identity on M_mu'][0] += 1
            stats['T(1,1,2) + T(1,2,2) is a multiple of the identity on M_mu'][1] += 0 if isc else 1
            stats['T(1,1,2) + T(1,2,2) is a multiple of the identity on M_mu'][2] += 1 if m >= 2 else 0
            # (8) with the central cubic of the total gl(4) on V_mu (only if V_mu is small)
            if weyl_dim(mu) <= 700:
                c3mu = site_scalar(Rep.of(mu), inv_T((1, 1, 1), N), 1)
                pred = (c3mu - ctx['C3_1'] - ctx['C3_2']) / 3 + Fr(2, 3) * ctx['n1'] * ctx['n2'] - Fr(8, 3) * ctx['Omega']
                stats['T(1,1,2) + T(1,2,2) = (C3(mu) - C3(l) - C3(l\'))/3 + 2/3 n1 n2 - 8/3 Omega'][0] += 1
                stats['T(1,1,2) + T(1,2,2) = (C3(mu) - C3(l) - C3(l\'))/3 + 2/3 n1 n2 - 8/3 Omega'][1] += 0 if (isc and cval == pred) else 1
                stats['T(1,1,2) + T(1,2,2) = (C3(mu) - C3(l) - C3(l\'))/3 + 2/3 n1 n2 - 8/3 Omega'][2] += 1 if m >= 2 else 0
            # (9) remainder values of V3 - T112 (for display)
            rem = mat_fr_lin([(1, Rm['V3']), (-1, RT['T(1,1,2)'])], m)
            isr, rv = is_scalar_matrix(rem)
            if m >= 2 and len(shown) < 40 and (l1, l2) in (((3, 2, 1, 0), (3, 2, 1, 0)), ((2, 1, 0, 0), (2, 1, 0, 0)), ((4, 2, 1, 0), (4, 2, 1, 0))):
                shown.append((l1, l2, mu, m, isr, rv, ctx))
        print('  %s (x) %s: %d sectors checked  [%.1fs]' % (lam_str(l1), lam_str(l2), len([1 for mu, m in rs.items() if nW_of(l1, l2, mu) <= 260]), time.perf_counter() - t0), flush=True)
    print('\n  identity                                                                                    | sectors | with m>=2 | violations')
    for k, (ns, nv, nm2) in stats.items():
        print('  %-92s | %7d | %9d | %d' % (k, ns, nm2, nv))
    print('\n  V3 - T(1,1,2) on M_mu (a multiple of the identity; value = 1/2 n1^2 n2 - 1/2 n2 C2_1 - n1 Omega + 1/2 n1 n2 + Omega):')
    for l1, l2, mu, m, isr, rv, ctx in shown:
        pred = Fr(1, 2) * ctx['n1'] ** 2 * ctx['n2'] - Fr(1, 2) * ctx['n2'] * ctx['C2_1'] - ctx['n1'] * ctx['Omega'] + Fr(1, 2) * ctx['n1'] * ctx['n2'] + ctx['Omega']
        print('     %s x %s  mu=%-14s m=%d : multiple of the identity: %s, value %s, predicted %s' % (lam_str(l1), lam_str(l2), str(mu), m, isr, rv, pred))
    assert all(nv == 0 for (_, nv, _) in stats.values())
    return stats


# ---------------------------------------------------------------------------------------------------- G2.5: term lists (file + loader + re-verification)
import re

TERMS_OPS = [('S_%d_%d' % (k, j), k, j) for k in (2, 3, 4) for j in range(1, k + 1)]


def lt_shift(T, s):
    """multiply a Laurent-term list by t^s"""
    return {ww: {p + s: c for p, c in pw.items()} for ww, pw in T.items()}


def terms_for_file():
    """ordered dict name -> LT (the theta_1 Laurent coefficients S_{k,j} of the N = 4 model and the triple V3 = t S_{3,2}, V4a = t S_{4,3}, V4b = t^2 S_{4,2})"""
    L = laurent(4)
    ops = {}
    for name, k, j in TERMS_OPS:
        ops[name] = L[k][j]
    ops['V3'] = lt_shift(L[3][2], 1)
    ops['V4a'] = lt_shift(L[4][3], 1)
    ops['V4b'] = lt_shift(L[4][2], 2)
    return ops


def fmt_coef(pw):
    return '{' + ','.join('{%d,%d/%d}' % (p, Fr(c).numerator, Fr(c).denominator) for p, c in sorted(pw.items())) + '}'


def fmt_word_file(w1, w2):
    return '{' + ','.join('{%d,%d,%d}' % (1, i, j) for (i, j) in w1) + (',' if w1 and w2 else '') + ','.join('{%d,%d,%d}' % (2, i, j) for (i, j) in w2) + '}'


def write_terms(path=TERMS_FILE):
    ops = terms_for_file()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    counts = {name: len(T) for name, T in ops.items()}
    raw = laurent(4)
    hdr_lines = [
        '# gl4-gaudin-ops.txt -- exact term lists of the two-site gl(4) Gaudin algebra (WP-G2, written by gaudin_gl4_ops.py, 2026-10-22).  All numbers are exact rationals.',
        '#',
        '# CONVENTIONS (the gl(3) notes Gaudin-sl3-CG-notes.md section 1.1, generalised verbatim to N = 4)',
        '#   [E_ij, E_kl] = delta_jk E_il - delta_li E_kj; two sites V_lam (site 1) (x) V_lam\' (site 2); e^s_ij = E_ij acting on site s (s = 1, 2); e^1 and e^2 COMMUTE.',
        '#   untwisted Lax matrix L(u)_ij = e^1_ij/(u - theta_1) + e^2_ij/(u - theta_2);  t := theta_12 = theta_1 - theta_2.',
        '#   column-ordered quantum determinant with the derivative moved to the right (d f = f d + f\'):',
        '#       cdet(d_u - L(u)) = sum_{sigma in S_4} sgn(sigma) M_{sigma(1) 1} M_{sigma(2) 2} M_{sigma(3) 3} M_{sigma(4) 4},  M = d_u 1 - L(u),',
        '#                        = d^4 - S_1 d^3 + S_2 d^2 - S_3 d + S_4  (operator products kept in the order written).',
        '#   Laurent data at theta_1:  x := u - theta_1 (so u - theta_2 = x + t):  S_k(u) = sum_j S_{k,j} x^(-j) + (regular at theta_1).  S_{k,j} (k = 2,3,4; j = 1..k) are operators on V_lam (x) V_lam\'.',
        '#   theta_2 data: the SAME formulas with the sites exchanged and t -> -t:  S_{k,j}^{(theta_2)}(t) = [S_{k,j} with every site label 1 <-> 2](t -> -t).',
        '#   Every S_{k,j} is homogeneous of degree j - k in t, so t^(k-j) S_{k,j} is t-independent; the canonical triple (all t-independent) is',
        '#       V3 = t S_{3,2},   V4a = t S_{4,3},   V4b = t^2 S_{4,2}   (classically P_{2,1}, P_{3,1}, P_{2,2} - P_{3,1}).',
        '#   The residue relations sum_s sum_{i<=j} C(j-1,i-1) theta_s^(j-i) S_{k,i}^{(s)} = 0 (j = 1..k-1, S_k = O(u^-k) at infinity) hold identically in U(gl_4) (x) U(gl_4).',
        '#',
        '# FORMAT.  Records are TAB-separated.  "OP<TAB>name" starts an operator (S_<k>_<j> = S_{k,j} at theta_1, including its power of t; V3, V4a, V4b = the t-scaled triple, power 0).',
        '#   "TERM<TAB>coefficient<TAB>word": coefficient = exact Laurent polynomial in t written {{power,p/q},...}; word = {{site,i,j},...} = the product of the generators e^site_ij in PRODUCT ORDER',
        '#   (the leftmost factor acts last on a ket; matrices act on kets as columns).  The empty word {} is the identity.  Indices 1-based.',
        '#   Words are stored with all site-1 factors first, then all site-2 factors (each site\'s internal order is the cdet order; e^1, e^2 commute, so this is the same operator and',
        '#   terms that differ only by the interleaving of the two sites have been merged).  WP-G3 may therefore normal-order the site-2 factors alone.',
        '#   Number of site-2 factors per term: V3 = t S_{3,2} and V4a = t S_{4,3}: exactly one; V4b = t^2 S_{4,2}: one or two (maximal word lengths per operator, site 1 / site 2: ' + '; '.join('%s: %d/%d' % (n, max((len(w1) for (w1, w2) in T), default=0), max((len(w2) for (w1, w2) in T), default=0)) for n, T in ops.items()) + ').',
        '#',
        '# CHECKS (gaudin_gl4_ops.py sym | gl3 | comm | triple | struct | terms):',
        '#   N = 3: the same engine reproduces eq. (4) of the notes exactly (S_1, S_2, S_3 as noncommutative polynomials in l_ij^(d)); t S_{3,2} = V + (n1 n2 - Omega)/2, t S_{2,1} = n1 n2 - Omega exactly in U(gl_3)(x)U(gl_3).',
        '#   N = 4: all 18 Laurent coefficients (theta_1, theta_2; k = 2..4, j = 1..k) commute pairwise and with all E^tot_ij exactly on [1,0,0,0]x[1,1,0,0], [1,1,0,0]x[2,1,0,0], [2,1,0,0]^2, [2,1,1,0]^2.',
        '#   V3 = sum_I D(2,1)_I + (n1 n2 - Omega), V4a = D(3,1) + V3 - (n1 n2 - Omega)/3  (D = antisymmetrised forms, see RESULTS-G2-gaudin-gl4-ops.md).',
        '#',
        '# COUNTS (number of TERM records per operator): ' + '; '.join('%s: %d' % (n, c) for n, c in counts.items()),
        '# Raw (word, site-assignment) counts of the cdet before merging: S_2: %d, S_3: %d, S_4: %d.' % (raw_word_count(cdet_S(4)[2]), raw_word_count(cdet_S(4)[3]), raw_word_count(cdet_S(4)[4])),
        '#',
    ]
    with open(path, 'w', encoding='utf8', newline='\n') as f:
        for line in hdr_lines:
            f.write(line + '\n')
        for name, T in ops.items():
            f.write('OP\t%s\n' % name)
            for ww in sorted(T, key=lambda k: (len(k[0]) + len(k[1]), k[0], k[1])):
                f.write('TERM\t%s\t%s\n' % (fmt_coef(T[ww]), fmt_word_file(*ww)))
    return counts


_COEF_RE = re.compile(r'\{(-?\d+),(-?\d+)/(\d+)\}')
_WORD_RE = re.compile(r'\{(\d+),(\d+),(\d+)\}')


def load_terms(path=TERMS_FILE):
    """parse gl4-gaudin-ops.txt: returns {name: {(w1, w2): {power: Fraction}}} with w_s = tuple of (i, j) (site-s factors in product order; each site's internal order is kept)"""
    ops = {}
    cur = None
    with open(path, encoding='utf8') as f:
        for line in f:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            fld = line.split('\t')
            if fld[0] == 'OP':
                cur = fld[1]
                ops[cur] = {}
            elif fld[0] == 'TERM':
                coef = {int(p): Fr(int(a), int(b)) for p, a, b in _COEF_RE.findall(fld[1])}
                w1, w2 = [], []
                for s, i, j in _WORD_RE.findall(fld[2]):
                    (w1 if int(s) == 1 else w2).append((int(i), int(j)))
                key = (tuple(w1), tuple(w2))
                d = ops[cur].setdefault(key, {})
                for p, c in coef.items():
                    d[p] = d.get(p, 0) + c
    return ops


def stage_terms(path=TERMS_FILE):
    hdr('G2.5  TERM LISTS: %s' % os.path.relpath(path, os.path.join(HERE, '..', '..')))
    counts = write_terms(path)
    print('written; terms per operator: %s' % counts)
    print('file size %d bytes, %d lines' % (os.path.getsize(path), sum(1 for _ in open(path, encoding='utf8'))))
    loaded = load_terms(path)
    gen = terms_for_file()
    ok = (set(loaded) == set(gen))
    for name in gen:
        g_ = {ww: {p: Fr(c) for p, c in pw.items()} for ww, pw in gen[name].items()}
        ok = ok and (loaded[name] == g_)
    print('loader round trip: %d operators, loaded term lists identical to the generated ones: %s' % (len(loaded), ok))
    assert ok
    # re-verification on a sector: the loaded V3, V4a, V4b reproduce the matrices; commutativity + exact simple spectrum
    for (l1, l2, mu) in [((3, 2, 1, 0), (3, 2, 1, 0), (5, 4, 2, 1)), ((4, 2, 1, 0), (4, 2, 1, 0), (6, 4, 3, 1))]:
        R1, R2 = Rep.of(l1), Rep.of(l2)
        m = racah_speiser(l1, l2)[mu]
        hw = HWExact(R1, R2, mu)
        direct = triple_ops(R1, R2)
        mats_d, mats_l = {}, {}
        for name in ('V3', 'V4a', 'V4b'):
            up = lt_at(loaded[name], 1)                       # power 0: the t-independent operator (t^0)
            opl = Op2(up, R1, R2)
            mats_l[name] = restrict_exact(opl, hw)[0]
            mats_d[name] = restrict_exact(direct[name], hw)[0]
        same = all(mats_l[n] == mats_d[n] for n in mats_l)
        comm = all(mat_fr_zero(mat_fr_comm(mats_l[a], mats_l[b])) for a, b in itertools.combinations(['V3', 'V4a', 'V4b'], 2))
        M = [[sum((c * mats_l[n][i][j] for c, n in zip(GENERIC_COMBO, ['V3', 'V4a', 'V4b'])), Fr(0)) for j in range(m)] for i in range(m)]
        dM, nrM = poly_disc_and_real(exact_charpoly(M))
        d3, nr3 = poly_disc_and_real(exact_charpoly(mats_l['V3']))
        js = joint_spectrum([mats_l[n] for n in ('V3', 'V4a', 'V4b')], GENERIC_COMBO, 30)
        print('  re-verification on %s x %s, mu = %s (m = %d, nW = %d): loaded V3/V4a/V4b restricted == directly built: %s; commute exactly: %s; joint disc != 0 and all roots real: %s; V3 alone simple: %s; min joint distance %s' %
              (lam_str(l1), lam_str(l2), mu, m, hw.nW, same, comm, dM != 0 and nrM == m, d3 != 0 and nr3 == m, fmt_mp(js['jgap'], 6)))
        assert same and comm and dM != 0 and nrM == m
    # all nine S_{k,j}(theta_1) from the file vs build_ops at t = 31/35
    R1, R2 = Rep.of((2, 1, 1, 0)), Rep.of((2, 1, 1, 0))
    t = Fr(31, 35)
    ops_d = build_ops(R1, R2, 4, t)
    bad = 0
    for name, k, j in TERMS_OPS:
        Mf = Op2(lt_at(loaded[name], t), R1, R2).matrix_exact()
        Md = ops_d[(k, j, 1)].matrix_exact()
        bad += 0 if Mf == Md else 1
    print('  the nine S_{k,j}(theta_1) rebuilt from the FILE at t = 31/35 on [2,1,1,0]^2 equal build_ops: %d mismatches' % bad)
    assert bad == 0
    return counts


# ---------------------------------------------------------------------------------------------------- CLI
def stage_scaffold():
    hdr("Representation scaffold: the private gl3_reps-based construction of this module vs G1's glN_reps.GLNRep (read-only cross-check, if importable)")
    lams = [(1, 0, 0, 0), (1, 1, 0, 0), (2, 1, 0, 0), (2, 1, 1, 0), (2, 2, 1, 0), (3, 2, 0, 0), (3, 2, 1, 0), (4, 2, 1, 0)]
    try:
        import glN_reps as G
    except Exception as e:
        print('glN_reps not importable (%r): cross-check skipped; this module uses its own scaffold' % (e,))
        return
    tot = 0
    for lam in lams:
        R, S = Rep.of(lam), G.GLNRep.of(lam)
        bad = int(R.pats != S.pats) + int(R.wt != S.wt) + int(R.top != S.top)
        for ij in itertools.product(range(1, 5), repeat=2):
            bad += int(R.E[ij] != S.E[ij])
        bad += int(list(R.g) != list(S.g))
        tot += bad
        print('  %s: dim %d, patterns/weights/top, all 16 generators, metric: %d mismatches' % (lam_str(lam), R.dim, bad))
    prods = [((2, 1, 0, 0), (2, 1, 0, 0)), ((3, 2, 1, 0), (3, 2, 1, 0)), ((3, 2, 0, 0), (2, 1, 1, 0)), ((4, 2, 1, 0), (3, 2, 1, 0)), ((5, 3, 1, 0), (5, 3, 1, 0))]
    for l1, l2 in prods:
        same = racah_speiser(l1, l2) == G.racah_speiser(l1, l2)
        tot += 0 if same else 1
        print('  Racah-Speiser %s x %s identical to glN_reps: %s' % (lam_str(l1), lam_str(l2), same))
    print('  total mismatches: %d (this module uses its OWN scaffold; G1\'s glN_reps passes its own tests at the time of writing: 10 tests, 0 failed)' % tot)
    assert tot == 0


def stage_big():
    hdr('G2.3 (large multiplicity spaces, extras): [5,3,1,0]^2 (dim 129600, m <= 7), [5,3,2,0]^2 (m <= 8), [6,4,2,0]^2 (dim 531441, m <= 12)')
    sums = [triple_product((5, 3, 1, 0), (5, 3, 1, 0)), triple_product((5, 3, 2, 0), (5, 3, 2, 0)), triple_product((6, 4, 2, 0), (6, 4, 2, 0), others=False)]
    print_triple_summary(sums, 'G2.3 SUMMARY TABLE (large products)')
    return sums


STAGES = ['sym', 'scaffold', 'gl3', 'comm', 'triple', 'struct', 'terms', 'big', 'all']


def main(argv):
    cmd = argv[1] if len(argv) > 1 else 'all'
    t0 = time.perf_counter()
    import platform, sympy
    print('# gaudin_gl4_ops.py stage=%s  (Python %s, numpy %s, scipy %s, sympy %s, mpmath %s; %s)' % (cmd, platform.python_version(), np.__version__, __import__('scipy').__version__, sympy.__version__, mp.__version__, platform.platform()))
    if cmd in ('sym', 'all'):
        stage_sym()
    if cmd in ('scaffold', 'all'):
        stage_scaffold()
    if cmd in ('gl3', 'all'):
        stage_gl3()
    if cmd in ('comm', 'all'):
        stage_comm()
    if cmd in ('triple', 'all'):
        stage_triple(extras=True)
    if cmd in ('struct', 'all'):
        stage_struct()
    if cmd in ('terms', 'all'):
        stage_terms()
    if cmd in ('big', 'all'):
        stage_big()
    if cmd not in STAGES:
        print('usage: python -X utf8 gaudin_gl4_ops.py <stage>, stage in %s' % STAGES)
        return 1
    print('\ntotal time %.1fs' % (time.perf_counter() - t0))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
