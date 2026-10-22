"""slice_ops_gl4.py -- WP-G3: the two-site gl(N) Gaudin missing-label operators reduced to the PRV slice space M_mu, and the canonical copies.

THEOREM (reduction rule; SPEC-G3-slice-ops-gl4.md section 1).  Psi = sum_q psi_q (x) e_q a highest-weight (HW) vector of weight mu in V_lam (x) V_lam' (site 1 = V_lam,
site 2 = V_lam'), slice psi_top = (id (x) <top|) Psi in M_mu.  For a term c W1 W2 (W1 a word of site-1 generators, W2 a word of site-2 generators; e1, e2 commute):
normal-order W2 in U(gl_N) as sum_w c_w L H R (L lowering, H Cartan, R raising factors, class order lowering < Cartan < raising, fixed order inside a class); then
      slice(c W1 L H R Psi) = 0                                                         if L is not empty   (<top| E_ab = 0 for a > b),
                            = c * prod_{E_aa in H} lam'_a * (-1)^|R| * W1 e1_{r_|R|} ... e1_{r_1} psi_top      for R = e2_{r_1} ... e2_{r_|R|}
(<top| H = prod lam'_a <top|; for a < b e2_ab Psi = -e1_ab Psi since E^tot_ab Psi = 0, and e2 commutes with e1, so the rightmost raising factor is transferred first and the
transferred factors stand in REVERSED order to the right of W1).  Hence every operator O of the two-site algebra acts on the slice as the universal polynomial
Red(O) in the site-1 generators (coefficients in Q[lam']); for O in the commutant of the diagonal gl_N it maps M_mu to M_mu and O -> Red(O)|_{M_mu} is the matrix of O on
the multiplicity space in slice coordinates (algebra homomorphism: slice(O1 O2 Psi) = Red(O1) slice(O2 Psi) because O2 Psi is again HW).

Public API (lam, lamp, mu = tuples of ints; exact = Fractions unless stated)
  normal_order(word)                       {word: coefficient}, word = tuple of (i, j); class order lowering < Cartan < raising, lexicographic inside a class
  reduce_terms_poly(terms, N)              terms {(w1, w2): c} -> {site-1 word: {lam'-monomial (exponent tuple): coefficient}}   (universal, lam'-independent)
  reduce_terms(terms, lamp)                -> {site-1 word: Fraction}   (the polynomial evaluated at lam')
  vred_apply(rep1, red, v)                 exact action of a reduced term list (dict word -> coefficient, or ReducedOp) on a slice v = {GT index: Fraction}
  op_terms(name, N)                        t-free term list {(w1, w2): Fraction} of 'V3', 'V4a', 'V4b', 'S31', 'S41', 'S_k_j' (file ref/gl4-gaudin-ops.txt for N = 4, G2's symbolic engine otherwise)
  V_red(lam, lamp, mu, name, prv=None)     exact m x m matrix in G1's PRV basis (column j = V_red b_j in PRV coordinates; PRV.coords raises if the image leaves M_mu)
  triple_mats(lam, lamp, mu)               (prv, {'V3','V4a','V4b': matrix})
  joint_diag / copies(lam, lamp, mu, dps)  joint eigen-decomposition (mpmath) of the commuting triple on M_mu: the canonical copies (joint_diag_float: float64 variant)
  direct_sector / direct_scale / swapped_check   the direct route (G2's full two-site operators on the exact HW vectors); theta_2 data via names 'S_k_j_sw'
  c_closed / cprime_closed                 closed forms (fitted, not proved) of the scalars of t^2 S_{3,1} = -2 V3 + c and t^3 S_{4,1} = -3 (V4a + V4b) + c'
  structure_checks / generic_lamp_commutation   BGG intertwining E_i^k V(lam') = V(s_i.lam') E_i^k, commutation on whole weight spaces, simplicity of the spectrum there
  exchange_analysis / separation_counts    site exchange P on the copies (lam = lam'), separation of the copies by subsets of the triple
  write_copies_ref / parse_copies_ref / verify_copies_ref / mutation_test_copies   the reference files ref/gl4-cg-copies-*.txt and their independent verification
CLI:  python -X utf8 slice_ops_gl4.py <stage> [ext],  stage in  order | reduce | gl3 | direct | copies | bench | ref | verify | all   (RESULTS-G3-slice-ops-gl4.md describes every stage)
"""
import os, sys, time, itertools, math, re, random, platform
from fractions import Fraction as Fr
from functools import lru_cache, cmp_to_key
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import scipy.sparse as sps
import mpmath as mp
import glN_reps as gr
import slice_engine_gl4 as se
import gaudin_gl4_ops as go

REFDIR = se.REFDIR
DATE = '2026-10-22'
TRIPLE = ('V3', 'V4a', 'V4b')
GENERIC = tuple(go.GENERIC_COMBO)                 # (1, 7/19, 11/13): the generic rational combination of G2 (M = V3 + 7/19 V4a + 11/13 V4b)
NONADJ = [(1, 3), (2, 4), (1, 4), (3, 1), (4, 2), (4, 1)]
GL3DIR = os.path.normpath(os.path.join(HERE, '..', 'gaudin_sl3'))

# the products of G3.3 (lam, lamp)
G33_PRODUCTS = [((2, 1, 0, 0), (2, 1, 0, 0)), ((2, 1, 1, 0), (2, 1, 1, 0)), ((2, 2, 1, 0), (2, 1, 1, 0)), ((1, 1, 0, 0), (2, 1, 0, 0)),
                ((3, 2, 0, 0), (2, 1, 1, 0)), ((3, 2, 1, 0), (3, 2, 1, 0))]


def hdr(s):
    print('\n' + '=' * 150 + '\n' + s + '\n' + '=' * 150, flush=True)


def lam_str(lam):
    return '[' + ','.join(str(x) for x in lam) + ']'


def prod_str(lam, lamp):
    return lam_str(lam) + ('^2' if tuple(lam) == tuple(lamp) else 'x' + lam_str(lamp))


fq = se.fq
fmu = se.fmu


# ======================================================================================== 1. normal ordering (G3.1)
def _cls(g):
    """0 = lowering (i > j), 1 = Cartan, 2 = raising (i < j)"""
    return 0 if g[0] > g[1] else (1 if g[0] == g[1] else 2)


def gen_key(g):
    return (_cls(g), g[0], g[1])


def _comm(a, b):
    """[E_ij, E_kl] = delta_jk E_il - delta_li E_kj  as a list of (generator, coefficient)"""
    (i, j), (k, l) = a, b
    out = []
    if j == k:
        out.append(((i, l), 1))
    if l == i:
        out.append(((k, j), -1))
    return out


@lru_cache(maxsize=None)
def _normal_order(word):
    for k in range(len(word) - 1):
        a, b = word[k], word[k + 1]
        if gen_key(a) > gen_key(b):                 # E_a E_b = E_b E_a + [E_a, E_b]
            res = {}
            for w, c in _normal_order(word[:k] + (b, a) + word[k + 2:]):
                res[w] = res.get(w, 0) + c
            for g, cg in _comm(a, b):
                for w, c in _normal_order(word[:k] + (g,) + word[k + 2:]):
                    res[w] = res.get(w, 0) + cg * c
            return tuple((w, c) for w, c in res.items() if c != 0)
    return ((word, 1),)


def normal_order(word):
    """PBW reordering in U(gl(N)): word = tuple of generators (i, j) (product order, leftmost factor acts last) -> {word': integer coefficient} with every word' of the
    form (lowering factors)(Cartan factors)(raising factors), each class sorted lexicographically.  Memoised."""
    return dict(_normal_order(tuple(word)))


def clear_normal_order_cache():
    _normal_order.cache_clear()


def is_normal_ordered(w):
    return all(gen_key(w[k]) <= gen_key(w[k + 1]) for k in range(len(w) - 1))


def word_cols(rep, word):
    """exact matrix of the product of generators `word` (leftmost factor acts last) as a list of column dicts (GLNRep)"""
    cols = [{c: Fr(1)} for c in range(rep.dim)]
    for g in reversed(word):
        E = rep.E[g]
        new = []
        for col in cols:
            acc = {}
            for r, v in col.items():
                for r2, v2 in E[r].items():
                    acc[r2] = acc.get(r2, 0) + v * v2
            new.append({r: v for r, v in acc.items() if v != 0})
        cols = new
    return cols


def cols_lincomb(pairs, n):
    out = [dict() for _ in range(n)]
    for c, M in pairs:
        for j in range(n):
            o = out[j]
            for r, v in M[j].items():
                o[r] = o.get(r, 0) + c * v
    return [{r: v for r, v in o.items() if v != 0} for o in out]


# ======================================================================================== 2. the reduction (G3.2)
def split_word(w):
    """a normal-ordered word -> (number of lowering factors, tuple of Cartan indices a (E_aa), tuple of raising generators)"""
    nl = sum(1 for g in w if g[0] > g[1])
    car = tuple(g[0] for g in w if g[0] == g[1])
    rai = tuple(g for g in w if g[0] < g[1])
    return nl, car, rai


def reduce_terms_poly(terms, N, variant='rule'):
    """the universal reduction: terms {(w1, w2): coefficient} -> {site-1 word: {monomial: coefficient}}, monomial = exponent tuple (e_1..e_N) of lam'_1..lam'_N.
    Normal-order w2; words with a lowering factor die; Cartan E_aa -> lam'_a; raising factors move to site 1 in reversed order with the sign (-1)^|R|.
    variant != 'rule' is for the sensitivity test of the verification only (deliberately WRONG rules): 'no_sign' (omit (-1)^|R|), 'no_reverse' (transfer in the original order),
    'naive_order' (no normal ordering: a lowering factor anywhere in the word as it stands kills the term)."""
    out = {}
    for (w1, w2), c in terms.items():
        for w, cw in ([(w2, 1)] if variant == 'naive_order' else normal_order(w2).items()):
            nl, car, rai = split_word(w) if variant != 'naive_order' else (sum(1 for g in w if g[0] > g[1]), tuple(g[0] for g in w if g[0] == g[1]), tuple(g for g in w if g[0] < g[1]))
            if nl:
                continue
            mono = [0] * N
            for a in car:
                mono[a - 1] += 1
            nw = w1 + (tuple(rai) if variant == 'no_reverse' else tuple(reversed(rai)))
            sgn = 1 if variant == 'no_sign' else (-1 if len(rai) % 2 else 1)
            d = out.setdefault(nw, {})
            key = tuple(mono)
            d[key] = d.get(key, 0) + sgn * Fr(c) * cw
    res = {}
    for w, d in out.items():
        d = {k: v for k, v in d.items() if v != 0}
        if d:
            res[w] = d
    return res


def eval_poly(p, lamp):
    s = Fr(0)
    for mono, c in p.items():
        t = Fr(c)
        for a, e in enumerate(mono):
            if e:
                t *= Fr(lamp[a]) ** e
        s += t
    return s


def reduce_terms(terms, lamp):
    """site-1 term list {word: Fraction} of the reduction of `terms` at lam' = lamp (see reduce_terms_poly)"""
    lamp = tuple(int(x) for x in lamp)
    out = {}
    for w, p in reduce_terms_poly(terms, len(lamp)).items():
        v = eval_poly(p, lamp)
        if v != 0:
            out[w] = v
    return out


class ReducedOp:
    """a reduced operator {site-1 word: coefficient} stored as a prefix trie for the exact evaluation  W_p v = c_p v + sum_g E_g W_{pg} v  (word = g_1 g_2 ... g_k acts
    right to left; the node of the prefix p collects  sum_s c_{p s} (s v)  over the suffixes s that act first): every trie node costs ONE generator application,
    sums are formed before the generator is applied, branches that die (zero vector) are cut."""

    def __init__(self, red, name=None, lamp=None):
        self.red, self.name, self.lamp = red, name, lamp
        self.root = [Fr(0), {}]
        for w, c in red.items():
            node = self.root
            for g in w:
                node = node[1].setdefault(g, [Fr(0), {}])
            node[0] += c
        self.nwords = len(red)
        self.maxlen = max((len(w) for w in red), default=0)
        self.nnodes = self._count(self.root)

    @staticmethod
    def _count(node):
        return sum(1 + ReducedOp._count(ch) for ch in node[1].values())

    def apply(self, rep, v, stats=None):
        """exact action on the slice v ({GT index: Fraction}) in the GLNRep `rep` of site 1"""
        Eb = {}
        wt = rep.wt

        def gen(g, vec):
            if g[0] == g[1]:
                k = g[0] - 1
                return {r: x * wt[r][k] for r, x in vec.items() if wt[r][k] != 0}
            return rep.apply(g, vec)

        def ev(node):
            c, ch = node
            out = {}
            if c != 0:
                for r, x in v.items():
                    out[r] = c * x
            for g, child in ch.items():
                cc, cch = child
                if not cch:
                    e = Eb.get(g)
                    if e is None:
                        e = Eb[g] = gen(g, v)
                    if stats is not None:
                        stats['nodes'] += 1
                        if not e:
                            stats['dead'] += 1
                    for r, x in e.items():
                        out[r] = out.get(r, 0) + cc * x
                else:
                    sub = ev(child)
                    if stats is not None:
                        stats['nodes'] += 1
                    if sub:
                        e = gen(g, sub)
                        if stats is not None and not e:
                            stats['dead'] += 1
                        for r, x in e.items():
                            out[r] = out.get(r, 0) + x
                    elif stats is not None:
                        stats['dead'] += 1
            return {r: x for r, x in out.items() if x != 0}
        return ev(self.root)

    def apply_naive(self, rep, v):
        """word by word (reference implementation for the tests)"""
        out = {}
        for w, c in self.red.items():
            vec = v
            for g in reversed(w):
                if g[0] == g[1]:
                    k = g[0] - 1
                    vec = {r: x * rep.wt[r][k] for r, x in vec.items() if rep.wt[r][k] != 0}
                else:
                    vec = rep.apply(g, vec)
                if not vec:
                    break
            for r, x in vec.items():
                out[r] = out.get(r, 0) + c * x
        return {r: x for r, x in out.items() if x != 0}


def vred_apply(rep1, red, v):
    """exact action of the reduced operator `red` (a ReducedOp or a dict {site-1 word: coefficient}) on the slice v = {GT index: Fraction} (site-1 irrep rep1 = GLNRep of lam)"""
    op = red if isinstance(red, ReducedOp) else ReducedOp(red)
    return op.apply(rep1, v)


# ---- term lists
_FILE_TERMS = [None]
_ALIAS = {'V3': 'V3', 'V4A': 'V4a', 'V4B': 'V4b', 'S31': 'S_3_1', 'S41': 'S_4_1', 'S32': 'S_3_2', 'S42': 'S_4_2', 'S43': 'S_4_3', 'S21': 'S_2_1', 'S22': 'S_2_2',
          'S33': 'S_3_3', 'S44': 'S_4_4', 'V': 'S_3_2'}


def canon(name):
    return _ALIAS.get(name.upper(), name) if '_' not in name else name


def file_terms():
    if _FILE_TERMS[0] is None:
        _FILE_TERMS[0] = go.load_terms()
    return _FILE_TERMS[0]


def op_terms(name, N=4):
    """t-free term list {(w1, w2): Fraction} of the operator `name`.  N = 4: 'V3', 'V4a', 'V4b' (= t S_{3,2}, t S_{4,3}, t^2 S_{4,2}, power 0 in the file) and 'S_k_j' / 'S31' / 'S41'
    (S_{k,j} evaluated at t = 1, i.e. t^(k-j) S_{k,j}).  Other N: S_k_j from G2's symbolic engine (name 'V' = t S_{3,2} for N = 3).
    A name with the suffix '_sw' (e.g. 'S_4_2_sw') is the same term list with the two sites exchanged: the theta_2 Laurent coefficient up to the sign (-1)^(j-k) (G2: S^{(2)}_{k,j}(t) =
    (-t)^(j-k) swap(hat)); its site-2 words are the former site-1 words (products of up to three generators, still in the column order of the cdet)."""
    nm = canon(name)
    if nm.endswith('_sw'):                                   # the theta_2 data: sites exchanged (t-free hat swapped, G2's rule S^{(2)}_{k,j} = (-t)^(j-k) swap(hat))
        return go.up_swap(op_terms(nm[:-3], N))
    if N == 4 and nm in file_terms():
        return go.lt_at(file_terms()[nm], 1)
    m = re.fullmatch(r'S_(\d+)_(\d+)', nm)
    if not m:
        raise KeyError((name, N))
    k, j = int(m.group(1)), int(m.group(2))
    return go.lt_at(go.laurent(N)[k][j], 1)


_UNIV = {}
_RED = {}


def universal_reduction(name, N=4):
    """{word: {monomial: coeff}}: the lam'-independent reduction of the operator `name` (cached)"""
    key = (canon(name), N)
    if key not in _UNIV:
        _UNIV[key] = reduce_terms_poly(op_terms(name, N), N)
    return _UNIV[key]


def reduced_op(name, lamp):
    """ReducedOp of the operator `name` at lam' = lamp (cached).  lamp may be ANY weight (integers, Fractions, non-dominant): the reduction is a polynomial in lam'"""
    lamp = tuple(x if isinstance(x, int) else Fr(x) for x in lamp)
    key = (canon(name), lamp)
    if key not in _RED:
        red = {}
        for w, p in universal_reduction(name, len(lamp)).items():
            v = eval_poly(p, lamp)
            if v != 0:
                red[w] = v
        _RED[key] = ReducedOp(red, canon(name), lamp)
    return _RED[key]


def clear_caches():
    """cold start for timings: G1's caches, the normal-ordering memo, the universal / evaluated reductions, G2's scaffold reps"""
    se._clear_caches()
    clear_normal_order_cache()
    _UNIV.clear()
    _RED.clear()
    _OP2.clear()
    _SLICEOP.clear()
    go.Rep._cache.clear()


def vred_matrix(rep, op, prv, stats=None):
    """exact m x m matrix of the ReducedOp `op` on M_mu = span(prv) in PRV coordinates: column j = op b_j (PRV.coords raises ValueError if the image leaves M_mu,
    KeyError if it leaves the weight space)"""
    m = prv.m
    pos = {int(c): a for a, c in enumerate(prv.W)}
    cols = []
    for k in range(m):
        img = op.apply(rep, prv.vec(k), stats)
        vec = [Fr(0)] * prv.w
        for r, x in img.items():
            vec[pos[r]] = x
        cols.append(prv.coords(vec))
    return [[cols[j][i] for j in range(m)] for i in range(m)]


def V_red(lam, lamp, mu, name, prv=None):
    """exact m x m matrix (list of rows of Fractions) of the reduced operator `name` on M_mu in G1's PRV basis; V_red(lam, lamp, mu, 'V3') etc."""
    lam, lamp, mu = tuple(lam), tuple(lamp), tuple(mu)
    prv = prv or se.prv_space(lam, lamp, mu)
    if prv.m == 0:
        return []
    return vred_matrix(gr.GLNRep.of(lam), reduced_op(name, lamp), prv)


def triple_mats(lam, lamp, mu, prv=None, names=TRIPLE):
    lam, lamp, mu = tuple(lam), tuple(lamp), tuple(mu)
    prv = prv or se.prv_space(lam, lamp, mu)
    R = gr.GLNRep.of(lam)
    return prv, {n: (vred_matrix(R, reduced_op(n, lamp), prv) if prv.m else []) for n in names}


# ---- small exact matrix helpers
def mat_mul(A, B):
    return go.mat_fr_mul(A, B)


def mat_comm(A, B):
    return go.mat_fr_comm(A, B)


def mat_is_zero(A):
    return all(x == 0 for r in A for x in r)


def mat_scalar(A):
    """(is_scalar, c) for an m x m exact matrix"""
    m = len(A)
    c = A[0][0]
    return all(A[i][j] == (c if i == j else 0) for i in range(m) for j in range(m)), c


def mat_lin(terms, m):
    out = [[Fr(0)] * m for _ in range(m)]
    for c, M in terms:
        for i in range(m):
            for j in range(m):
                out[i][j] += c * M[i][j]
    return out


def mat_trace(A):
    return sum((A[i][i] for i in range(len(A))), Fr(0))


def charpoly_asc(A):
    """exact characteristic polynomial det(x - A), coefficients in ASCENDING powers [c0, ..., c_{m-1}, 1]"""
    return list(reversed(go.exact_charpoly(A)))


# ======================================================================================== 3. G3.1 stage: normal ordering tests
def test_normal_order(nrandom=400, seed=11, verbose=True):
    """G3.1: (a) exhaustive length-2 words and `nrandom` random words of length <= 4 in gl(4): the normal-ordered sum equals the original word as an EXACT matrix on [2,1,0,0] and
    [1,1,0,0] (GLNRep) and in the PBW normal form of G2 (representation-free, a different total order); every output word is (lowering)(Cartan)(raising) sorted; normal_order is
    idempotent on its own output.  (b) all distinct length-3 site-2 words occurring in S_{4,1} (up to three site-2 factors)."""
    rng = random.Random(seed)
    gens = [(i, j) for i in range(1, 5) for j in range(1, 5)]
    reps = [gr.GLNRep.of((2, 1, 0, 0)), gr.GLNRep.of((1, 1, 0, 0))]
    wc = {}

    def wcols(ri, w):
        k = (ri, w)
        if k not in wc:
            wc[k] = word_cols(reps[ri], w)
        return wc[k]

    def check(word):
        no = normal_order(word)
        bad = dict(rep=0, pbw=0, order=0, idem=0)
        for ri, R in enumerate(reps):
            lhs = wcols(ri, word)
            rhs = cols_lincomb([(c, wcols(ri, w)) for w, c in no.items()], R.dim)
            if lhs != rhs:
                bad['rep'] += 1
        acc = {}
        for w, c in no.items():
            for u, d in go.pbw_word(w).items():
                acc[u] = acc.get(u, 0) + c * d
        acc = {u: c for u, c in acc.items() if c != 0}
        if acc != go.pbw_word(tuple(word)):
            bad['pbw'] += 1
        for w in no:
            if not is_normal_ordered(w):
                bad['order'] += 1
            if normal_order(w) != {w: 1}:
                bad['idem'] += 1
        return bad, len(no)

    words = [tuple(p) for p in itertools.product(gens, repeat=2)]
    n_ex = len(words)
    for _ in range(nrandom):
        words.append(tuple(rng.choice(gens) for _ in range(rng.randint(1, 4))))
    tot = Counter()
    maxterms = 0
    for w in words:
        b, nt = check(w)
        tot.update(b)
        maxterms = max(maxterms, nt)
    s41 = sorted({w2 for (w1, w2) in file_terms()['S_4_1'] if len(w2) == 3})
    tot3 = Counter()
    nt3 = 0
    for w in s41:
        b, nt = check(w)
        tot3.update(b)
        nt3 = max(nt3, nt)
    # the example of the brief
    ex = normal_order(((1, 2), (2, 1)))
    ex_ok = (ex == {((2, 1), (1, 2)): 1, ((1, 1),): 1, ((2, 2),): -1})
    ex2 = normal_order(((2, 1), (1, 2)))
    ex2_ok = (ex2 == {((2, 1), (1, 2)): 1})
    res = dict(n_words=len(words), n_exhaustive2=n_ex, bad=dict(tot), maxterms=maxterms, n_s41=len(s41), bad_s41=dict(tot3), maxterms_s41=nt3, example_ok=ex_ok and ex2_ok)
    if verbose:
        print('  words tested: %d (all %d words of length 2 + %d random words of length 1..4, seed %d); exact-matrix mismatches on [2,1,0,0] and [1,1,0,0]: %d, PBW (representation-free, G2) mismatches: %d, '
              'class-order violations: %d, idempotence violations: %d; max number of terms of a normal-ordered word: %d' %
              (len(words), n_ex, nrandom, seed, tot['rep'], tot['pbw'], tot['order'], tot['idem'], maxterms))
        print('  the %d distinct length-3 site-2 words of S_{4,1}: exact-matrix mismatches %d, PBW mismatches %d, order violations %d, idempotence violations %d; max terms of a normal form %d' %
              (len(s41), tot3['rep'], tot3['pbw'], tot3['order'], tot3['idem'], nt3))
        print('  e12 e21 = e21 e12 + e11 - e22 (raising left of lowering): %s;  e21 e12 is already normal ordered: %s' % (ex_ok, ex2_ok))
    return res


def stage_order():
    hdr('G3.1  NORMAL ORDERING in U(gl(4)): class order lowering < Cartan < raising, PBW reordering by commutators (memoised)')
    t0 = time.perf_counter()
    r = test_normal_order()
    nbad = sum(r['bad'].values()) + sum(r['bad_s41'].values()) + (0 if r['example_ok'] else 1)
    print('  normal-ordering failures in total: %d   (%.2fs)' % (nbad, time.perf_counter() - t0))
    return r


# ======================================================================================== 4. G3.2: reduction tables, closed forms, sparsity
def fmt_gen(g):
    return 'E%d%d' % g


def fmt_word(w):
    return ' '.join(fmt_gen(g) for g in w) if w else '1'


def fmt_monopoly(p, N=4):
    """polynomial in lam'_1..lam'_N from {monomial: coefficient}"""
    parts = []
    for mono in sorted(p, key=lambda m_: (-sum(m_), m_), reverse=False):
        c = p[mono]
        mon = '*'.join(('l%d' % (a + 1) if e == 1 else 'l%d^%d' % (a + 1, e)) for a, e in enumerate(mono) if e)
        cs = fq(c)
        parts.append(('%s%s' % ('+' if c > 0 else '-', fq(abs(c)))) + ('*' + mon if mon else ''))
    return ' '.join(parts) if parts else '0'


def hand_reduce_single(terms, lamp):
    """the ONE-LINE rule for terms with at most one site-2 factor (no normal ordering, independent of normal_order): lowering -> 0, Cartan E_aa -> lam'_a, raising E_ab -> -E_ab
    appended to the right of the site-1 word.  Raises ValueError on a term with more than one site-2 factor."""
    out = {}
    for (w1, w2), c in terms.items():
        c = Fr(c)
        if len(w2) == 0:
            nw, cc = w1, c
        elif len(w2) == 1:
            g = w2[0]
            if g[0] > g[1]:
                continue
            if g[0] == g[1]:
                nw, cc = w1, c * lamp[g[0] - 1]
            else:
                nw, cc = w1 + (g,), -c
        else:
            raise ValueError('more than one site-2 factor')
        if cc != 0:
            out[nw] = out.get(nw, 0) + cc
    return {w: c for w, c in out.items() if c != 0}


def closed_forms(N=4):
    """G2's exact identities in U(gl4)(x)U(gl4) (RESULTS-G2 section 6): V3 = sum_I D(2,1)_I + eta, V4a = D(3,1) + V3 - eta/3, eta = n1 n2 - Omega, as UP term lists"""
    eta = go.up_lin((1, go.up_mul(go.inv_n(1, N), go.inv_n(2, N))), (-1, go.inv_Omega(N)))
    V3 = go.up_lin((1, go.antisym_sum_subsets(N, 2, 1)), (1, eta))
    V4a = go.up_lin((1, go.antisym_form(N, 3, 1)), (1, V3), (Fr(-1, 3), eta))
    return dict(V3=V3, V4a=V4a, eta=eta)


def reduction_table(lamps, verbose=True):
    rows = []
    names = ['V3', 'V4a', 'V4b', 'S31', 'S41']
    for name in names:
        terms = op_terms(name, 4)
        w2s = {w2 for (w1, w2) in terms}
        maxlen2 = max(len(w2) for w2 in w2s)
        nl = 0
        for (w1, w2) in terms:
            if all(split_word(w)[0] for w in normal_order(w2)):
                nl += 1
        clear_normal_order_cache()
        t0 = time.perf_counter()
        U = reduce_terms_poly(terms, 4)
        t_poly = time.perf_counter() - t0
        row = dict(name=name, nterms=len(terms), nw2=len(w2s), maxlen2=maxlen2, nkilled=nl, nuniv=len(U), univmaxlen=max(len(w) for w in U), t_poly=t_poly, per_lamp={})
        for lamp in lamps:
            t0 = time.perf_counter()
            op = reduced_op(name, lamp)
            row['per_lamp'][lamp] = (op.nwords, op.nnodes, op.maxlen, time.perf_counter() - t0)
        rows.append(row)
        if verbose:
            print('  %-4s terms %4d | distinct site-2 words %3d (max length %d) | terms killed outright (every normal-ordered site-2 word has a lowering factor) %4d | universal reduced words %3d (max length %d, %.1f ms) | %s' %
                  (name, row['nterms'], row['nw2'], row['maxlen2'], row['nkilled'], row['nuniv'], row['univmaxlen'], 1e3 * t_poly,
                   ' | '.join('lam\'=%s: %d words, %d trie nodes, max len %d' % (lam_str(lp), *v[:3]) for lp, v in row['per_lamp'].items())), flush=True)
    return rows


def closed_form_check(lam, lamp, verbose=True):
    """V_red,3 and V_red,4a (engine: the 88 / 152-term lists of G2's file, reduced by the engine) against the hand-reduced closed forms V3 = sum_I D(2,1)_I + eta,
    V4a = D(3,1) + V3 - eta/3 (G2 section 6; single site-2 factor: the one-line rule, no normal ordering) on every PRV basis vector of every irrep"""
    cf = closed_forms()
    R = gr.GLNRep.of(lam)
    hand = {n: ReducedOp(hand_reduce_single(cf[n], lamp)) for n in ('V3', 'V4a')}
    rs = gr.racah_speiser(lam, lamp)
    nvec = nbad = 0
    for mu in sorted(rs, reverse=True):
        prv = se.prv_space(lam, lamp, mu)
        for n in ('V3', 'V4a'):
            op = reduced_op(n, lamp)
            for k in range(prv.m):
                v = prv.vec(k)
                nvec += 1
                if op.apply(R, v) != hand[n].apply(R, v):
                    nbad += 1
    if verbose:
        print('  %-18s closed forms (sum_I D(2,1)_I + eta ; D(3,1) + V3 - eta/3, hand-reduced: %d / %d words) vs engine: %d (operator, PRV vector) pairs, mismatches %d' %
              (prod_str(lam, lamp), len(hand['V3'].red), len(hand['V4a'].red), nvec, nbad), flush=True)
    return nvec, nbad


def closed_form_V4b_check(lam, lamp, verbose=True):
    """where the NORMAL ORDERING really matters: G2's antisymmetrised representative of the other quartic operator,
        V4b = D(2,2) - D(3,1) + 1/2 sum_I D(1,2)_I - 1/2 sum_I D(2,1)_I - eta/6      (G2 section 6b; exact in U(gl4)(x)U(gl4)),
    is a DIFFERENT polynomial from the cdet term list of the file; its site-2 words (D(2,2): e2_{s3 t3} e2_{s4 t4} with arbitrary indices) contain raising factors to the left of lowering factors with
    nonvanishing commutators (e.g. e2_12 e2_21).  Reduced with the rule (normal ordering) it must give the same operator on every M_mu as the file's V4b; reduced with the naive rule (no normal
    ordering) it must not.  Returns (#PRV vectors, #images differing with the rule, #images differing with the naive rule)"""
    lam, lamp = tuple(lam), tuple(lamp)
    N = len(lam)
    eta = go.up_lin((1, go.up_mul(go.inv_n(1, N), go.inv_n(2, N))), (-1, go.inv_Omega(N)))
    V4b = go.up_lin((1, go.antisym_form(N, 2, 2)), (-1, go.antisym_form(N, 3, 1)), (Fr(1, 2), go.antisym_sum_subsets(N, 1, 2)), (Fr(-1, 2), go.antisym_sum_subsets(N, 2, 1)), (Fr(-1, 6), eta))
    R = gr.GLNRep.of(lam)
    ref = reduced_op('V4b', lamp)
    ops = {}
    for var in ('rule', 'naive_order'):
        U = reduce_terms_poly(V4b, N, var)
        ops[var] = ReducedOp({w: v for w, v in ((w, eval_poly(p, lamp)) for w, p in U.items()) if v != 0})
    nv = bad_rule = bad_naive = 0
    for mu in sorted(gr.racah_speiser(lam, lamp), reverse=True):
        prv = se.prv_space(lam, lamp, mu)
        for k in range(prv.m):
            v = prv.vec(k)
            r0 = ref.apply(R, v)
            nv += 1
            bad_rule += (ops['rule'].apply(R, v) != r0)
            bad_naive += (ops['naive_order'].apply(R, v) != r0)
    if verbose:
        print('   %-20s G2\'s antisymmetrised V4b (%d terms; site-2 words of length 2 in arbitrary index order): %d PRV vectors: images differ from the file\'s V4b with the RULE (normal ordering): %d, with the NAIVE rule: %d' %
              (prod_str(lam, lamp), len(V4b), nv, bad_rule, bad_naive), flush=True)
    return nv, bad_rule, bad_naive


def sparsity_report(lam, lamp, verbose=True, names=TRIPLE):
    """(a) trie nodes visited / dead (zero vector, cut) while building the triple on every M_mu; (b) (reduced word, weight-space basis vector) pairs with w e_c = 0 ('word acts trivially
    on the weight space V_lam[mu - lam'])"""
    R = gr.GLNRep.of(lam)
    rs = gr.racah_speiser(lam, lamp)
    out = {}
    for name in names:
        op = reduced_op(name, lamp)
        st = dict(nodes=0, dead=0)
        pairs = zero_pairs = 0
        nsec = 0
        for mu in sorted(rs, reverse=True):
            prv = se.prv_space(lam, lamp, mu)
            for k in range(prv.m):
                op.apply(R, prv.vec(k), st)
            nsec += 1
            for w in op.red:
                for c in prv.W:
                    pairs += 1
                    vec = {int(c): Fr(1)}
                    for g in reversed(w):
                        if g[0] == g[1]:
                            kk = g[0] - 1
                            vec = {r: x * R.wt[r][kk] for r, x in vec.items() if R.wt[r][kk] != 0}
                        else:
                            vec = R.apply(g, vec)
                        if not vec:
                            break
                    if not vec:
                        zero_pairs += 1
        out[name] = dict(nodes=st['nodes'], dead=st['dead'], pairs=pairs, zero_pairs=zero_pairs, nwords=op.nwords, nsectors=nsec)
        if verbose:
            print('  %-18s %-4s %3d reduced words: trie nodes evaluated %6d, of which dead (zero vector / zero image) %6d (%.0f%%); (word, weight-space basis vector) pairs %7d, with w e_c = 0: %7d (%.0f%%)' %
                  (prod_str(lam, lamp), name, op.nwords, st['nodes'], st['dead'], 100.0 * st['dead'] / max(st['nodes'], 1), pairs, zero_pairs, 100.0 * zero_pairs / max(pairs, 1)), flush=True)
    return out


def variant_sensitivity(lam, lamp, variants=('no_sign', 'no_reverse', 'naive_order'), names=('V3', 'V4a', 'V4b', 'S31', 'S41')):
    """a test that cannot fail proves nothing: reduce with a deliberately WRONG rule (see reduce_terms_poly) and run the same two tests as the real ones on every sector of the product:
    (i) does the image of a PRV basis vector leave M_mu (PRV.coords raises)?  (ii) does the slice differ from the slice read off the FULL two-site operator on the exact HW vector?
    Returns {variant: {operator: (#vectors with image outside M_mu, #vectors with a wrong slice, #vectors)}}"""
    lam, lamp = tuple(lam), tuple(lamp)
    N = len(lam)
    R = gr.GLNRep.of(lam)
    rs = gr.racah_speiser(lam, lamp)
    out = {}
    for var in variants:
        ops = {}
        for n in names:
            U = reduce_terms_poly(op_terms(n, N), N, var)
            ops[n] = ReducedOp({w: v for w, v in ((w, eval_poly(p, lamp)) for w, p in U.items()) if v != 0})
        res = {n: [0, 0, 0] for n in names}
        for mu in sorted(rs, reverse=True):
            prv = se.prv_space(lam, lamp, mu)
            psis = [se.hw_from_slice_exact(lam, lamp, mu, prv.basis[k])[0] for k in range(prv.m)]
            pos = {int(c): a for a, c in enumerate(prv.W)}
            for n in names:
                sop = slice_op(n, lam, lamp)
                for j in range(prv.m):
                    img = ops[n].apply(R, prv.vec(j))
                    res[n][2] += 1
                    try:
                        vec = [Fr(0)] * prv.w
                        for r, x in img.items():
                            vec[pos[r]] = x
                        prv.coords(vec)
                    except (ValueError, KeyError):
                        res[n][0] += 1
                    res[n][1] += (sop.apply(by_q(psis[j])) != img)
        out[var] = {n: tuple(v) for n, v in res.items()}
    return out


def weight_space_mats(R, ops, W):
    """exact w x w matrices of the ReducedOps `ops` (dict name -> ReducedOp) on the WHOLE weight space spanned by the GT basis vectors W (global indices) of the GLNRep R"""
    pos = {int(c): a for a, c in enumerate(W)}
    out = {}
    for n, op in ops.items():
        cols = []
        for c in W:
            img = op.apply(R, {int(c): Fr(1)})
            col = [Fr(0)] * len(W)
            for r, x in img.items():
                col[pos[r]] = x
            cols.append(col)
        out[n] = [[cols[j][i] for j in range(len(W))] for i in range(len(W))]
    return out


def power_apply(R, ij, k, vec):
    for _ in range(k):
        vec = R.apply(ij, vec)
        if not vec:
            break
    return vec


def structure_checks(lam, lamp, verbose=True):
    """structure of the reduced operators, beyond the brief (exact arithmetic):
    (1) they commute on the WHOLE weight space V_lam[mu - lam'] of every candidate mu, not only on M_mu;
    (2) the intertwining relation  E_i^(k_i) V_red(lam') = V_red(s_i.lam') E_i^(k_i)  on V_lam[mu - lam'] for i = 1..3, k_i = lam'_i - lam'_{i+1} + 1, s_i.lam' = lam' - k_i alpha_i
        (the Verma embedding M(s_i.lam') -> M(lam'); it implies that M_mu = intersection of the kernels is invariant);
    (3) their joint spectrum on the whole weight space (50 digits): simple?  and the joint eigenvalue triples of M_mu are among them."""
    lam, lamp = tuple(lam), tuple(lamp)
    N = len(lam)
    R = gr.GLNRep.of(lam)
    rs = gr.racah_speiser(lam, lamp)
    opsA = {n: reduced_op(n, lamp) for n in TRIPLE}
    st = Counter()
    for mu in gr.candidates(lam, lamp, both=True):
        nu = tuple(mu[a] - lamp[a] for a in range(N))
        W = R.wsp.get(nu)
        if W is None:
            continue
        st['spaces'] += 1
        w = len(W)
        m = rs.get(mu, 0)
        M = weight_space_mats(R, opsA, W)
        st['comm_bad'] += sum(0 if mat_is_zero(mat_comm(M[a], M[b])) else 1 for a, b in itertools.combinations(TRIPLE, 2))
        st['comm_pairs'] += 3
        for i in range(1, N):
            k = lamp[i - 1] - lamp[i] + 1
            slamp = list(lamp)
            slamp[i - 1] -= k
            slamp[i] += k
            opsB = {n: reduced_op(n, tuple(slamp)) for n in TRIPLE}
            for n in TRIPLE:
                for c in W:
                    v = {int(c): Fr(1)}
                    lhs = power_apply(R, (i, i + 1), k, opsA[n].apply(R, v))
                    rhs = opsB[n].apply(R, power_apply(R, (i, i + 1), k, v))
                    st['inter_tests'] += 1
                    st['inter_nontrivial'] += bool(lhs or rhs)
                    st['inter_bad'] += (lhs != rhs)
        jd = joint_diag([M[n] for n in TRIPLE], dps=50)
        gap = min_joint_gap(jd['trip'])
        st['simple'] += (w == 1) or (gap > mp.mpf(10) ** -20)
        st['w_total'] += w
        st['m_total'] += m
        if m:
            cps = copies(lam, lamp, mu)
            st['contained'] += all(any(max(abs(a - b) for a, b in zip(cp['ev'], t)) < mp.mpf(10) ** -30 for t in jd['trip']) for cp in cps)
        else:
            st['contained'] += 1
    if verbose:
        print('   %-20s %3d weight spaces V_lam[mu-lam\']: reduced operators commute on the whole space: %d of %d pairs fail | intertwining E_i^k V(lam\') = V(s_i.lam\') E_i^k: %d tests (%d nonzero), %d fail | joint spectrum on the whole weight space simple: %d of %d spaces; '
              'M_mu triples contained in it: %d of %d; dim sum w = %d, of which %d not in M_mu (sum m)' %
              (prod_str(lam, lamp), st['spaces'], st['comm_bad'], st['comm_pairs'], st['inter_tests'], st['inter_nontrivial'], st['inter_bad'], st['simple'], st['spaces'], st['contained'], st['spaces'], st['w_total'],
               st['w_total'] - st['m_total']), flush=True)
    return st


def generic_lamp_commutation(lam, lamps, verbose=True):
    """the three reduced operators commute on every weight space V_lam[nu] of dimension >= 2 for GENERIC (rational, non-dominant) lam' (the reduction is a polynomial in lam')"""
    R = gr.GLNRep.of(tuple(lam))
    tested = bad = 0
    for lamp in lamps:
        ops = {n: reduced_op(n, lamp) for n in TRIPLE}
        for w_, W in sorted(R.wsp.items()):
            if len(W) < 2:
                continue
            M = weight_space_mats(R, ops, W)
            for a, b in itertools.combinations(TRIPLE, 2):
                tested += 1
                bad += (not mat_is_zero(mat_comm(M[a], M[b])))
    if verbose:
        print('   V_lam = %-10s generic lam\' (%d values, e.g. %s): commutator tests (pair x weight space of dim >= 2 x lam\'): %d, failures %d' %
              (lam_str(lam), len(lamps), tuple(str(x) for x in lamps[0]), tested, bad), flush=True)
    return tested, bad


def stage_reduce():
    hdr('G3.2  THE REDUCTION: universal polynomials in the site-1 generators, term counts, closed-form cross-check, sparsity')
    lamps = [(2, 1, 1, 0), (3, 2, 1, 0), (6, 4, 2, 0), (12, 6, 0, 0), (8, 8, 0, 0)]
    print('  term counts of the reduction (original term list -> normal ordering of the site-2 words -> universal polynomial in the site-1 generators with coefficients in Q[lam\'])')
    rows = reduction_table(lamps)
    for nm_ in ('V3', 'V4a', 'V4b'):
        U = universal_reduction(nm_)
        print('\n  universal V_red,%s (lam\'-independent): %s reduces to %d words; the polynomial (word = product of site-1 generators E_ij, leftmost acts last; coefficient = polynomial in lam\'_1..lam\'_4 written l1..l4):' %
              (nm_[1:], nm_, len(U)))
        for w in sorted(U, key=lambda w_: (len(w_), w_)):
            print('     %-34s %s' % (fmt_word(w), fmt_monopoly(U[w])))
    print('\n  closed-form cross-check (independent of the engine for the one-site-2-factor operators):')
    tot = [0, 0]
    for lam, lamp in G33_PRODUCTS + [((4, 2, 0, 0), (4, 2, 0, 0)), ((4, 2, 1, 0), (4, 2, 1, 0)), ((6, 4, 2, 0), (6, 4, 2, 0))]:
        a, b = closed_form_check(lam, lamp)
        tot[0] += a
        tot[1] += b
    print('  closed-form cross-check: %d (operator, vector) pairs, %d mismatches' % tuple(tot))
    print('\n  where the normal ordering really matters: G2\'s antisymmetrised representative of V4b (a different polynomial from the file\'s cdet list, site-2 words in arbitrary index order):')
    t4 = [0, 0, 0]
    for lam, lamp in G33_PRODUCTS + [((4, 2, 0, 0), (4, 2, 0, 0)), ((4, 2, 1, 0), (4, 2, 1, 0)), ((6, 4, 2, 0), (6, 4, 2, 0))]:
        a, b, c = closed_form_V4b_check(lam, lamp)
        t4[0] += a
        t4[1] += b
        t4[2] += c
    print('  antisymmetrised V4b: %d PRV vectors; rule (normal ordering) differs from the file\'s V4b on %d, naive rule (no normal ordering) on %d' % tuple(t4))
    print('\n  sparsity (how often the reduced words act trivially):')
    sp = {}
    for lam, lamp in [((3, 2, 1, 0), (3, 2, 1, 0)), ((4, 2, 0, 0), (4, 2, 0, 0)), ((6, 4, 2, 0), (6, 4, 2, 0))]:
        sp[(lam, lamp)] = sparsity_report(lam, lamp)
    print('\n  sensitivity of the tests to the rule (deliberately WRONG rules, same two tests as for the real rule: image outside M_mu / slice differs from the full two-site operator; entries = failing vectors of all).')
    print('  (the rule "naive_order" = no normal ordering is NOT detected on the cdet term lists: in a column-ordered cdet word a raising factor never stands left of a lowering factor it fails to commute with,')
    print('   because a raising factor E_{a c} (a < c) and a lowering factor E_{r c\'} with c\' > c have r > c\' > c > a, so [E_{ac}, E_{rc\'}] = 0 -- the same holds for the theta_2 words (sites exchanged, same word order);')
    print('   it is detected by G2\'s antisymmetrised V4b above and by the random words of G3.1)')
    sens = {}
    for lam, lamp in [((2, 1, 0, 0), (2, 1, 0, 0)), ((3, 2, 1, 0), (3, 2, 1, 0)), ((2, 2, 1, 0), (2, 1, 1, 0))]:
        sens[(lam, lamp)] = variant_sensitivity(lam, lamp)
        for var, d in sens[(lam, lamp)].items():
            print('   %-20s rule %-12s ' % (prod_str(lam, lamp), var) + '; '.join('%s: %d/%d outside M_mu, %d/%d wrong slice' % (n, v[0], v[2], v[1], v[2]) for n, v in d.items()), flush=True)
    print('\n  structure of the reduced operators beyond the brief (exact arithmetic; on the WHOLE weight space V_lam[mu - lam\'] of every candidate mu, not only on M_mu):')
    stc = {}
    for lam, lamp in G33_PRODUCTS + [((4, 2, 0, 0), (4, 2, 0, 0)), ((4, 2, 1, 0), (4, 2, 1, 0)), ((6, 4, 2, 0), (6, 4, 2, 0))]:
        stc[(lam, lamp)] = structure_checks(lam, lamp)
    gen = []
    for lam in [(3, 2, 1, 0), (4, 2, 0, 0), (2, 1, 1, 0)]:
        lamps = [(Fr(1, 3), Fr(5, 7), Fr(-2), Fr(4, 5)), (Fr(7, 2), Fr(-1, 3), Fr(2), Fr(11, 4)), (Fr(0), Fr(0), Fr(0), Fr(0)), (Fr(-3), Fr(5), Fr(1, 2), Fr(-7, 3))]
        gen.append(generic_lamp_commutation(lam, lamps))
    return rows, tot, sp, stc, gen


# ======================================================================================== 5. G3.2 gl(3): the convention test
def _se3():
    if GL3DIR not in sys.path:
        sys.path.insert(0, GL3DIR)
    import slice_engine as s3
    return s3


def gl3_convention_test(lams=((2, 1, 0), (3, 1, 0), (4, 2, 0)), verbose=True):
    """N = 3: the term list of t S_{3,2} (G2's symbolic engine, N = 3) reduced by THIS engine == slice_engine.vred_apply + (eta_1/2) 1 exactly on every PRV basis vector of every irrep;
    S_{3,1} (two-site-2-factor words) satisfies t^2 S_{3,1,red} = -2 (reduced t S_{3,2}) + c 1 on every M_mu, c computed per sector; the integer -2 is re-verified by the direct route
    (the full two-site S_{3,1}, S_{3,2} on the exact HW vectors)."""
    s3 = _se3()
    out = dict(vec=0, bad_conv=0, nsec=0, bad_rel=0, cs={}, direct_bad=0, direct_n=0)
    for lam in lams:
        lam = tuple(lam)
        R = gr.GLNRep.of(lam)
        rep3 = s3.Rep.of(lam)
        op32 = reduced_op('S32', lam)
        op31 = reduced_op('S31', lam)
        rs = gr.racah_speiser(lam, lam)
        nvec = nbad = nrel = 0
        t32 = op_terms('S32', 3)
        t31 = op_terms('S31', 3)
        R1, R2 = go.Rep.of(lam), go.Rep.of(lam)
        o32, o31 = go.Op2(t32, R1, R2), go.Op2(t31, R1, R2)
        top = R2.top
        cs = {}
        for mu in sorted(rs, reverse=True):
            prv = se.prv_space(lam, lam, mu)
            eta = s3.eta1(lam, mu)
            out['nsec'] += 1
            A32 = vred_matrix(R, op32, prv)
            A31 = vred_matrix(R, op31, prv)
            m = prv.m
            for k in range(m):
                v = prv.vec(k)
                lhs = op32.apply(R, v)
                rhs = dict(s3.vred_apply(rep3, v))
                for r, x in v.items():
                    rhs[r] = rhs.get(r, 0) + eta / 2 * x
                rhs = {r: x for r, x in rhs.items() if x != 0}
                nvec += 1
                nbad += (lhs != rhs)
            D = mat_lin([(1, A31), (2, A32)], m)             # S31 + 2 (V + eta/2) = c 1 ?
            isc, c = mat_scalar(D)
            nrel += (not isc)
            cs[mu] = (m, c, c - eta, isc)
            # direct route for the integer: full two-site operators on the exact HW vectors
            psis = [se.hw_from_slice_exact(lam, lam, mu, prv.basis[k])[0] for k in range(m)]
            for A, o in ((A32, o32), (A31, o31)):
                for j in range(m):
                    y = o.apply_exact(psis[j])
                    rec = {}
                    for i in range(m):
                        if A[i][j] != 0:
                            for kk, vv in psis[i].items():
                                rec[kk] = rec.get(kk, 0) + A[i][j] * vv
                    rec = {kk: vv for kk, vv in rec.items() if vv != 0}
                    out['direct_n'] += 1
                    out['direct_bad'] += (rec != y)
        out['vec'] += nvec
        out['bad_conv'] += nbad
        out['bad_rel'] += nrel
        out['cs'][lam] = cs
        if verbose:
            print('  N=3  %-10s: %2d irreps, %2d PRV basis vectors: reduced tS32 == slice_engine.vred_apply + eta1/2: mismatches %d;  t^2 S31_red + 2 (tS32_red) = c 1 on every M_mu: failures %d' %
                  (prod_str(lam, lam), len(rs), nvec, nbad, nrel), flush=True)
            print('       c per irrep (m): ' + ', '.join('%s: %s (m=%d, c-eta1 = %s)' % (fmu(mu), fq(v[1]), v[0], fq(v[2])) for mu, v in sorted(cs.items(), reverse=True)))
    if verbose:
        print('  N=3 totals: %d PRV vectors, convention-test mismatches %d; scalar-relation failures %d (on %d sectors); direct route (full two-site tS32 and S31 on exact HW vectors, slice AND full vector) mismatches %d of %d' %
              (out['vec'], out['bad_conv'], out['bad_rel'], out['nsec'], out['direct_bad'], out['direct_n']))
    return out


def stage_gl3():
    hdr('G3.2  gl(3) CONVENTION TEST: the reduced N = 3 t S_{3,2} equals slice_engine.vred_apply + eta_1/2 exactly; S_{3,1}')
    t0 = time.perf_counter()
    r = gl3_convention_test()
    print('  (%.2fs)' % (time.perf_counter() - t0))
    return r


# ======================================================================================== 6. joint diagonalisation and the canonical copies (G3.4)
def mpq(x):
    return mp.mpf(x.numerator) / x.denominator


def tol_cmp(tol):
    """lexicographic comparison of tuples of mpf with a tolerance (values closer than tol are treated as equal)"""
    def cmp(a, b):
        for x, y in zip(a, b):
            if abs(x - y) > tol:
                return -1 if x < y else 1
        return 0
    return cmp_to_key(cmp)


def joint_diag(mats, dps=50, combo=GENERIC, maxtries=3):
    """joint eigen-decomposition of the commuting exact m x m matrices `mats` = [A_1..A_r] (PRV coordinates, column j = image of b_j) by mpmath (dps digits) of the generic rational
    combination sum combo_i A_i; then EVERY A_i is conjugated to the eigenbasis and the largest off-diagonal entry measured (acceptance 10^-(dps-10) relative to max(1, |A|)); if it is
    not met the precision is raised by 30 digits (at most maxtries times).  Returns dict(trip = [tuple of r real mpf per eigenvector], U = [columns: lists of m real mpf, each normalised
    to its largest component 1], offmax, imag (max imaginary part met), dps (used), ok)"""
    m = len(mats[0])
    r = len(mats)
    if m == 1:
        old = mp.mp.dps
        mp.mp.dps = dps
        try:
            return dict(trip=[tuple(mpq(A[0][0]) for A in mats)], U=[[mp.mpf(1)]], offmax=mp.mpf(0), imag=mp.mpf(0), dps=dps, ok=True)
        finally:
            mp.mp.dps = old
    scale = max([1] + [abs(float(x)) for A in mats for row in A for x in row])
    for attempt in range(maxtries):
        d = dps + 30 * attempt
        old = mp.mp.dps
        mp.mp.dps = d
        try:
            M = [[sum((Fr(c) * A[i][j] for c, A in zip(combo, mats)), Fr(0)) for j in range(m)] for i in range(m)]
            E, ER = mp.eig(go.fr_mp(M))
            imag = max([abs(mp.im(E[k])) for k in range(m)] + [mp.mpf(0)])
            cols = []
            for k in range(m):
                big = max(range(m), key=lambda i: abs(ER[i, k]))
                ref = ER[big, k]
                col = [ER[i, k] / ref for i in range(m)]
                imag = max(imag, max(abs(mp.im(x)) for x in col))
                cols.append([mp.re(x) for x in col])
            Um = mp.matrix(m, m)
            for k in range(m):
                for i in range(m):
                    Um[i, k] = cols[k][i]
            Uinv = mp.inverse(Um)
            trip = [[None] * r for _ in range(m)]
            offmax = mp.mpf(0)
            for t, A in enumerate(mats):
                D = Uinv * go.fr_mp(A) * Um
                for i in range(m):
                    trip[i][t] = D[i, i]
                    for j in range(m):
                        if i != j:
                            offmax = max(offmax, abs(D[i, j]))
            ok = offmax < mp.mpf(10) ** (-(dps - 10)) * scale and imag < mp.mpf(10) ** (-(dps - 10)) * scale
            res = dict(trip=[tuple(row) for row in trip], U=cols, offmax=offmax, imag=imag, dps=d, ok=ok)
        finally:
            mp.mp.dps = old
        if ok:
            break
    return res


def joint_diag_float(mats, combo=GENERIC):
    """float64 (numpy) version of joint_diag, for the timing of a double-precision copies pipeline: returns dict(trip = (m x r) array of joint eigenvalues, U = (m x m) real eigenvector
    columns scaled to unit max-norm, offmax = largest off-diagonal entry of U^-1 A_i U relative to max(1, |A_i|), imag = largest imaginary part met)"""
    m = len(mats[0])
    Af = [np.array([[float(x) for x in row] for row in A], dtype=float).reshape(m, m) for A in mats]
    if m == 1:
        return dict(trip=np.array([[A[0, 0] for A in Af]]), U=np.ones((1, 1)), offmax=0.0, imag=0.0)
    M = sum(float(c) * A for c, A in zip(combo, Af))
    E, U = np.linalg.eig(M)
    imag = float(np.abs(E.imag).max())
    if np.iscomplexobj(U):
        big = np.abs(U).argmax(axis=0)
        U = U / U[big, np.arange(m)]
        imag = max(imag, float(np.abs(U.imag).max()))
        U = U.real
    U = U / np.abs(U).max(axis=0)
    Uinv = np.linalg.inv(U)
    trip = np.zeros((m, len(Af)))
    off = 0.0
    for t, A in enumerate(Af):
        D = Uinv @ A @ U
        trip[:, t] = np.diag(D)
        off = max(off, float(np.abs(D - np.diag(np.diag(D))).max()) / max(1.0, float(np.abs(A).max())))
    return dict(trip=trip, U=U, offmax=off, imag=imag)


def min_joint_gap(trips):
    """smallest max-norm distance between two joint eigenvalue tuples (inf for m = 1)"""
    best = mp.inf
    for i in range(len(trips)):
        for j in range(i + 1, len(trips)):
            best = min(best, max(abs(a - b) for a, b in zip(trips[i], trips[j])))
    return best


SEP_SUBSETS = [('V3', (0,)), ('V4a', (1,)), ('V4b', (2,)), ('V3,V4a', (0, 1)), ('V3,V4b', (0, 2)), ('V4a,V4b', (1, 2)), ('V3,V4a,V4b', (0, 1, 2))]


def separation_counts(cps_by_mu, tol_digits=20):
    """over all sectors with m >= 2: in how many is the spectrum of the subset of the commuting triple DEGENERATE (two copies whose eigenvalues of the subset agree to 10^-tol_digits)?
    returns (number of sectors with m >= 2, {subset name: number of sectors where it fails to separate the copies})"""
    tol = mp.mpf(10) ** (-tol_digits)
    cnt = Counter()
    n = 0
    for mu, cps in cps_by_mu.items():
        m = len(cps)
        if m < 2:
            continue
        n += 1
        for name, idx in SEP_SUBSETS:
            deg = False
            for i in range(m):
                for j in range(i + 1, m):
                    if max(abs(cps[i]['ev'][t] - cps[j]['ev'][t]) for t in idx) < tol:
                        deg = True
                        break
                if deg:
                    break
            cnt[name] += deg
    return n, {name: cnt[name] for name, idx in SEP_SUBSETS}


class Copies(list):
    """list of canonical copies (dicts: ev = (ev3, ev4a, ev4b), coords = free entries of the slice (PRV coordinates), slice = {global GT index: mpf}, first nonzero component 1);
    .info = diagnostics of the joint diagonalisation, .prv, .mats"""
    pass


def copies(lam, lamp, mu, dps=50, prv=None, mats=None):
    """the canonical copies of the irrep mu in V_lam (x) V_lam': joint eigenvectors of (V_red,3, V_red,4a, V_red,4b) on M_mu, ordered lexicographically by (ev3, ev4a, ev4b) ascending
    (ties of one eigenvalue are resolved by the next one, with tolerance 10^-(dps/2)), slice normalised to first nonzero component (ascending GT index) = 1.  Everything exact until the
    mpmath eigen-decomposition (`joint_diag`, `dps` digits; the precision is raised automatically if the diagonality test 10^-(dps-10) fails)."""
    lam, lamp, mu = tuple(lam), tuple(lamp), tuple(mu)
    if prv is None:
        prv = se.prv_space(lam, lamp, mu)
    if mats is None:
        _, mats = triple_mats(lam, lamp, mu, prv)
    m = prv.m
    jd = joint_diag([mats[n] for n in TRIPLE], dps=dps)
    old = mp.mp.dps
    mp.mp.dps = max(dps, jd['dps'])
    try:
        key = tol_cmp(mp.mpf(10) ** (-(dps // 2)))
        order = sorted(range(m), key=lambda k: key(jd['trip'][k]))
        Bm = [[mpq(x) for x in row] for row in prv.basis]
        tolz = mp.mpf(10) ** (-(dps - 20))
        out = Copies()
        smallest = mp.inf
        struct0 = [all(prv.basis[kk][a] == 0 for kk in range(m)) for a in range(prv.w)]       # components that vanish on the whole of M_mu
        for k in order:
            col = jd['U'][k]
            sl = [mp.fsum(col[kk] * Bm[kk][a] for kk in range(m)) for a in range(prv.w)]
            mx = max(abs(x) for x in sl)
            first = next(a for a in range(prv.w) if not struct0[a] and abs(sl[a]) > tolz * mx)
            sc = sl[first]
            sl = [x / sc for x in sl]
            mx = max(abs(x) for x in sl)
            nz = {int(prv.W[a]): sl[a] for a in range(prv.w) if not struct0[a] and abs(sl[a]) > tolz * mx}
            if nz:
                smallest = min(smallest, min(abs(x) for x in nz.values()) / mx)
            coords = [sl[f] for f in prv.free]
            out.append(dict(ev=tuple(jd['trip'][k]), coords=coords, slice=nz))
        out.info = dict(offmax=jd['offmax'], imag=jd['imag'], dps=jd['dps'], ok=jd['ok'], jgap=min_joint_gap(jd['trip']), smallest=smallest)
        out.prv, out.mats = prv, mats
        return out
    finally:
        mp.mp.dps = old


def copy_float_slices(prv, cps):
    """(w, m) float array of the copies' slices on the weight space (columns in the order of cps)"""
    pos = {int(c): a for a, c in enumerate(prv.W)}
    F = np.zeros((prv.w, len(cps)))
    for j, cp in enumerate(cps):
        for c, x in cp['slice'].items():
            F[pos[c], j] = float(x)
    return F


# ======================================================================================== 7. direct route (G3.3): the full two-site operators on the reconstructed HW vectors
_OP2 = {}


def two_site_op(name, lam, lamp):
    """G2's Kronecker-sum builder Op2 of the operator `name` (t-free term list at t = 1) on V_lam (x) V_lam' (G2's private GT scaffold, entry-by-entry identical to GLNRep)"""
    key = (canon(name), tuple(lam), tuple(lamp))
    if key not in _OP2:
        _OP2[key] = go.Op2(op_terms(name, len(lam)), go.Rep.of(tuple(lam)), go.Rep.of(tuple(lamp)))
    return _OP2[key]


class SliceOp:
    """G2's Op2 reduced to what the slice read-off needs: for every Kronecker group (P, Q) of O = sum_g P_g (x) Q_g the site-1 matrix P (columns) and Qtop = {q: <top|Q|q>};
    slice(O X)[r] = sum_g sum_{(nu,q) in X} P_g[r,nu] <top|Q_g|q> X[nu,q]  (the same numbers as the components (r, top) of Op2.apply_exact(X), without forming the image)"""

    def __init__(self, op, top):
        self.pairs = []
        for P, Q in op.exact_pairs():
            qt = {q: col[top] for q, col in enumerate(Q) if top in col}
            if qt:
                self.pairs.append((P, qt))

    def apply(self, X_byq):
        out = {}
        for P, qt in self.pairs:
            for q, c in qt.items():
                lst = X_byq.get(q)
                if not lst:
                    continue
                for nu, v in lst:
                    vq = v * c
                    for r, a in P[nu].items():
                        out[r] = out.get(r, 0) + vq * a
        return {r: x for r, x in out.items() if x != 0}


def by_q(X):
    """two-site vector {(nu, q): v} -> {q: [(nu, v)]} (zeros dropped)"""
    d = {}
    for (nu, q), v in X.items():
        if v != 0:
            d.setdefault(q, []).append((nu, v))
    return d


_SLICEOP = {}


def slice_op(name, lam, lamp):
    key = (canon(name), tuple(lam), tuple(lamp))
    if key not in _SLICEOP:
        _SLICEOP[key] = SliceOp(two_site_op(name, lam, lamp), go.Rep.of(tuple(lamp)).top)
    return _SLICEOP[key]


def op_slice_exact(name, lam, lamp, X):
    """exact slice {nu: Fraction} of the full two-site operator `name` applied to the exact two-site vector X = {(nu, q): Fraction}"""
    return slice_op(name, lam, lamp).apply(by_q(X))


DIRECT_OPS = ('V3', 'V4a', 'V4b', 'S31', 'S41')


def direct_sector(lam, lamp, mu, names=DIRECT_OPS, prv=None, mats=None, full=True, psis=None):
    """G3.3 for ONE sector, EXACT: the exact HW vectors Psi_k of the PRV basis (G1's exact twin), the full two-site operator (G2's Op2, term lists from the file) applied exactly
    (full=True: the whole image vector; else only its slice components via op_slice_exact); compared with the reduced operator:
        slice_bad = number of basis vectors j for which slice(O Psi_j) != sum_i A_ij b_i   (A = V_red matrix),
        full_bad  = number of j for which O Psi_j != sum_i A_ij Psi_i as whole two-site vectors (only if full=True; implies that O maps the HW space into itself).
    Returns dict(name -> (slice_bad, full_bad)), the matrices A and the exact HW vectors."""
    lam, lamp, mu = tuple(lam), tuple(lamp), tuple(mu)
    prv = prv or se.prv_space(lam, lamp, mu)
    R = gr.GLNRep.of(lam)
    top = gr.GLNRep.of(lamp).top
    if mats is None:
        mats = {n: vred_matrix(R, reduced_op(n, lamp), prv) for n in names}
    if psis is None:
        psis = [se.hw_from_slice_exact(lam, lamp, mu, prv.basis[k])[0] for k in range(prv.m)]
    pos = {int(c): a for a, c in enumerate(prv.W)}
    res = {}
    for n in names:
        op = two_site_op(n, lam, lamp)
        A = mats[n]
        sb = fb = 0
        for j in range(prv.m):
            if full:
                y = op.apply_exact(psis[j])
                sl = {r[0]: v for r, v in y.items() if r[1] == top}
            else:
                y = None
                sl = op_slice_exact(n, lam, lamp, psis[j])
            img = {}
            for i in range(prv.m):
                if A[i][j] != 0:
                    for c, x in prv.vec(i).items():
                        img[c] = img.get(c, 0) + A[i][j] * x
            img = {c: x for c, x in img.items() if x != 0}
            if sl != img:
                sb += 1
            if full:
                rec = {}
                for i in range(prv.m):
                    if A[i][j] != 0:
                        for kk, vv in psis[i].items():
                            rec[kk] = rec.get(kk, 0) + A[i][j] * vv
                rec = {kk: vv for kk, vv in rec.items() if vv != 0}
                if rec != y:
                    fb += 1
        res[n] = (sb, fb)
    return res, mats, psis


def direct_sector_float(lam, lamp, mu, names=TRIPLE, prv=None, mats=None):
    """float direct route: float HW vectors of the PRV basis (G1's recursion), G2's float operators (orthonormal basis, Op2.apply_float) applied to the whole two-site vector;
    returns {name: max over j of ||slice(O Psi_j) - V_red b_j|| / max(||V_red b_j||, ||b_j||)} in the orthonormal basis of site 1 (V_red b_j from the exact reduction, converted to float)"""
    lam, lamp, mu = tuple(lam), tuple(lamp), tuple(mu)
    prv = prv or se.prv_space(lam, lamp, mu)
    R = gr.GLNRep.of(lam)
    if mats is None:
        mats = {n: vred_matrix(R, reduced_op(n, lamp), prv) for n in names}
    hw = se.hw_from_slice(lam, lamp, mu, prv.float_basis())
    sg1, sg2 = hw.sg1, hw.sg2
    top = gr.GLNRep.of(lamp).top
    W = np.asarray(prv.W)
    out = {}
    for n in names:
        op = two_site_op(n, lam, lamp)
        A = np.array([[float(x) for x in row] for row in mats[n]], dtype=float).reshape(prv.m, prv.m)
        worst = 0.0
        for j in range(prv.m):
            X = np.zeros((hw.dim1, hw.dim2))
            X[hw.rows, hw.cols] = hw.vals_hat[:, j]
            Y = op.apply_float(X)
            sl_hat = Y[W, top]                                  # orthonormal components of the slice of the image (site 1 orthonormal, site 2 = top)
            ref = (prv.float_basis() @ A[:, j]) * sg1[W] * sg2[top]
            nb = np.linalg.norm(prv.float_basis()[:, j] * sg1[W] * sg2[top])       # relative to the norm of the slice itself (V4a vanishes identically when lam_3 = 0)
            worst = max(worst, np.linalg.norm(sl_hat - ref) / max(np.linalg.norm(ref), nb))
        out[n] = worst
    return out


# ---- G2's route: exact HW basis of the sector in the two-site space, full operators restricted exactly
def g2_sector(lam, lamp, mu, dps=40, names=TRIPLE):
    """G2's independent route for one sector: exact HW basis (HWExact), the triple as full two-site operators restricted exactly (restrict_exact) -> exact matrices in G2's basis, their
    characteristic polynomials (also of the generic combination), joint eigenvalue triples (mpmath, dps digits) and the eigenvectors' SLICES (Psi[(nu, top)], normalised to first
    nonzero component 1) of the joint eigenbasis."""
    lam, lamp, mu = tuple(lam), tuple(lamp), tuple(mu)
    R1, R2 = go.Rep.of(lam), go.Rep.of(lamp)
    hw = go.HWExact(R1, R2, mu)
    m = hw.m
    mats = {}
    viol = 0
    for n in names:
        M, v = go.restrict_exact(two_site_op(n, lam, lamp), hw)
        mats[n] = M
        viol += v
    comb = mat_lin([(c, mats[n]) for c, n in zip(GENERIC, names)], m)
    out = dict(m=m, mats=mats, viol=viol, cp={n: charpoly_asc(mats[n]) for n in names}, cpJ=charpoly_asc(comb), hw=hw)
    jd = joint_diag([mats[n] for n in names], dps=dps)
    out['jd'] = jd
    # slices of the eigenvectors: Psi = sum_i u_i hw.basis[i]  (u in G2's canonical HW basis coordinates)
    old = mp.mp.dps
    mp.mp.dps = max(dps, jd['dps'])
    try:
        key = tol_cmp(mp.mpf(10) ** (-(dps // 2)))
        order = sorted(range(m), key=lambda k: key(jd['trip'][k]))
        W = sorted(nu for nu in range(R1.dim) if tuple(R1.wt[nu][a] + R2.wt[R2.top][a] for a in range(len(mu))) == mu)
        sls, trips = [], []
        for k in order:
            col = jd['U'][k]
            sl = {}
            for i in range(m):
                for (nu, q), v in hw.basis[i].items():
                    if q == R2.top:
                        sl[nu] = sl.get(nu, 0) + col[i] * mpq(v)
            mx = max([abs(x) for x in sl.values()] + [mp.mpf(0)])
            nzs = {nu: x for nu, x in sl.items() if abs(x) > mp.mpf(10) ** (-(dps - 20)) * mx}
            first = min(nzs)
            sc = nzs[first]
            sls.append({nu: x / sc for nu, x in nzs.items()})
            trips.append(tuple(jd['trip'][k]))
        out['slices'] = sls
        out['trips'] = trips
    finally:
        mp.mp.dps = old
    return out


def relations(mats, m):
    """the scalar relations of G2 section 6c on one sector: S31 + 2 V3 = c 1, S41 + 3 (V4a + V4b) = c' 1  (t^2 S_{3,1} = -2 V3 + c, t^3 S_{4,1} = -3 (V4a + V4b) + c').  Returns (c, ok3, c', ok4)"""
    D3 = mat_lin([(1, mats['S31']), (2, mats['V3'])], m)
    D4 = mat_lin([(1, mats['S41']), (3, mats['V4a']), (3, mats['V4b'])], m)
    ok3, c3 = mat_scalar(D3)
    ok4, c4 = mat_scalar(D4)
    return c3, ok3, c4, ok4


def run_direct_product(lam, lamp, names=DIRECT_OPS, g2=True, full=True, verbose=True, listing=True, maxlines=40):
    """G3.3 for a product: every irrep, every PRV basis vector, every operator in `names`: exact direct route; M_mu membership; exact commutation of the triple; the scalar relations;
    spectrum against G2 (exact characteristic polynomials, 40-digit joint eigenvalues)."""
    lam, lamp = tuple(lam), tuple(lamp)
    R = gr.GLNRep.of(lam)
    rs = gr.racah_speiser(lam, lamp)
    st = Counter()
    maxdiff = mp.mpf(0)
    cs = []
    lines = []
    t0 = time.perf_counter()
    for mu in sorted(rs, reverse=True):
        prv = se.prv_space(lam, lamp, mu)
        m = prv.m
        mats = {}
        for n in names:
            try:
                mats[n] = vred_matrix(R, reduced_op(n, lamp), prv)
            except (ValueError, KeyError):
                st['leave_M_' + n] += 1
                mats[n] = None
        if any(mats[n] is None for n in names):
            continue
        st['sectors'] += 1
        st['vectors'] += m
        dres, _, psis = direct_sector(lam, lamp, mu, names, prv=prv, mats=mats, full=full)
        for n in names:
            st['slice_bad_' + n] += dres[n][0]
            st['full_bad_' + n] += dres[n][1]
        comm_bad = 0
        for a, b in itertools.combinations(names, 2):
            if not mat_is_zero(mat_comm(mats[a], mats[b])):
                comm_bad += 1
        st['comm_bad'] += comm_bad
        st['comm_pairs'] += len(list(itertools.combinations(names, 2)))
        c3, ok3, c4, ok4 = relations(mats, m) if all(n in mats for n in ('S31', 'S41', 'V3', 'V4a', 'V4b')) else (None, True, None, True)
        st['rel_bad'] += (not ok3) + (not ok4)
        cs.append((mu, m, c3, c4))
        line = None
        if g2:
            ref = g2_sector(lam, lamp, mu)
            mine = [mats[n] for n in TRIPLE]
            cpbad = sum(1 for n in TRIPLE if charpoly_asc(mats[n]) != ref['cp'][n])
            combo = mat_lin([(c, mats[n]) for c, n in zip(GENERIC, TRIPLE)], m)
            cpbad += (charpoly_asc(combo) != ref['cpJ'])
            jd = joint_diag(mine, dps=50)
            key = tol_cmp(mp.mpf(10) ** (-25))
            order = sorted(range(m), key=lambda k: key(jd['trip'][k]))
            mine_t = [jd['trip'][k] for k in order]
            dmax = max([abs(a - b) for t1, t2 in zip(mine_t, ref['trips']) for a, b in zip(t1, t2)] + [mp.mpf(0)])
            maxdiff = max(maxdiff, dmax)
            st['cp_bad'] += cpbad
            st['g2_viol'] += ref['viol']
            st['g2_sectors'] += 1
            if m >= 2:
                line = '    mu=%-14s m=%2d  V3/V4a/V4b/generic exact char. polys equal G2: %s   joint-triple max |diff| vs G2 (40 digits): %s   offdiag(mine) %s' % (
                    str(mu), m, cpbad == 0, mp.nstr(dmax, 3), mp.nstr(jd['offmax'], 3))
        if line:
            lines.append(line)
    if verbose:
        print('  %-20s sectors %3d, PRV vectors %3d | direct route, exact: slice mismatches ' % (prod_str(lam, lamp), st['sectors'], st['vectors']) +
              ', '.join('%s %d' % (n, st['slice_bad_' + n]) for n in names) +
              (' | full-vector mismatches ' + ', '.join('%s %d' % (n, st['full_bad_' + n]) for n in names) if full else ' (slice read-off only)') +
              ' | M_mu left: %d | non-commuting pairs: %d of %d | scalar-relation failures: %d' % (sum(v for k, v in st.items() if k.startswith('leave_M')), st['comm_bad'], st['comm_pairs'], st['rel_bad']) +
              ((' | vs G2: %d sectors, charpoly mismatches %d, restriction violations %d, max joint-eigenvalue diff %s' % (st['g2_sectors'], st['cp_bad'], st['g2_viol'], mp.nstr(maxdiff, 3))) if g2 else '') +
              '  (%.1fs)' % (time.perf_counter() - t0), flush=True)
        if listing:
            for ln in lines[:maxlines]:
                print(ln, flush=True)
    return dict(stats=st, maxdiff=maxdiff, cs=cs)


def pp_casimir(mu, r):
    """Perelomov-Popov: eigenvalue of the chain trace sum E_{i1 i2} E_{i2 i3} ... E_{ir i1} on V_mu: sum_i l_i^r prod_{j != i} (1 - 1/(l_i - l_j)), l_i = mu_i + N - i"""
    N = len(mu)
    l = [Fr(mu[i] + N - 1 - i) for i in range(N)]
    s = Fr(0)
    for i in range(N):
        t = l[i] ** r
        for j in range(N):
            if j != i:
                t *= 1 - Fr(1) / (l[i] - l[j])
        s += t
    return s


def c_closed(lam, lamp, mu):
    """closed form (exact polynomial fit, found by this stage, NOT proved) of the scalar c of  t^2 S_{3,1} = -2 V3 + c  on the irrep mu of V_lam (x) V_lam' (gl(4)); n1 = |lam|, n2 = |lam'|,
    C_r = Perelomov-Popov chain Casimirs:
      c = 1/2 n1^2 n2 + 7/3 n1 n2 + 1/2 n1 n2^2 + 1/2 n1 C2(lam) + 1/2 n2 C2(lam') + 5/3 (C2(lam) + C2(lam')) - 1/3 (C3(lam) + C3(lam')) - (5/3 + (n1 + n2)/2) C2(mu) + 1/3 C3(mu)"""
    n1, n2 = Fr(sum(lam)), Fr(sum(lamp))
    C2a, C2b, C2m = pp_casimir(lam, 2), pp_casimir(lamp, 2), pp_casimir(mu, 2)
    C3a, C3b, C3m = pp_casimir(lam, 3), pp_casimir(lamp, 3), pp_casimir(mu, 3)
    return (Fr(1, 2) * n1 * n1 * n2 + Fr(7, 3) * n1 * n2 + Fr(1, 2) * n1 * n2 * n2 + Fr(1, 2) * n1 * C2a + Fr(1, 2) * n2 * C2b + Fr(5, 3) * (C2a + C2b) - Fr(1, 3) * (C3a + C3b)
            - (Fr(5, 3) + (n1 + n2) / 2) * C2m + Fr(1, 3) * C3m)


_CP_TERMS = [(Fr(1, 6), ('n1', 'n1', 'n1', 'n2')), (Fr(5, 4), ('n1', 'n1', 'n2')), (Fr(1, 4), ('n1', 'n1', 'n2', 'n2')), (Fr(1, 4), ('n1', 'n1', 'C21')), (Fr(-1, 4), ('n1', 'n1', 'C2m')),
             (Fr(11, 4), ('n1', 'n2')), (Fr(5, 4), ('n1', 'n2', 'n2')), (Fr(1, 6), ('n1', 'n2', 'n2', 'n2')), (Fr(-1, 2), ('n1', 'n2', 'C2m')), (Fr(5, 3), ('n1', 'C21')),
             (Fr(-5, 3), ('n1', 'C2m')), (Fr(-1, 3), ('n1', 'C31')), (Fr(1, 3), ('n1', 'C3m')), (Fr(1, 4), ('n2', 'n2', 'C22')), (Fr(-1, 4), ('n2', 'n2', 'C2m')), (Fr(5, 3), ('n2', 'C22')),
             (Fr(-5, 3), ('n2', 'C2m')), (Fr(-1, 3), ('n2', 'C32')), (Fr(1, 3), ('n2', 'C3m')), (Fr(11, 4), ('C21',)), (Fr(-1, 8), ('C21', 'C21')), (Fr(11, 4), ('C22',)),
             (Fr(-1, 8), ('C22', 'C22')), (Fr(-11, 4), ('C2m',)), (Fr(1, 8), ('C2m', 'C2m')), (Fr(-3, 2), ('C31',)), (Fr(-3, 2), ('C32',)), (Fr(3, 2), ('C3m',)), (Fr(1, 4), ('C41',)),
             (Fr(1, 4), ('C42',)), (Fr(-1, 4), ('C4m',))]


def cprime_closed(lam, lamp, mu):
    """closed form (exact fit of weight <= 4 in n1, n2 and the Perelomov-Popov Casimirs C2, C3, C4 of lam, lam', mu; found by this stage, NOT proved) of the scalar c' of
    t^3 S_{4,1} = -3 (V4a + V4b) + c' on the irrep mu of V_lam (x) V_lam' (gl(4)); the 31 monomials are listed in _CP_TERMS"""
    ctx = dict(n1=Fr(sum(lam)), n2=Fr(sum(lamp)), C21=pp_casimir(lam, 2), C22=pp_casimir(lamp, 2), C2m=pp_casimir(mu, 2), C31=pp_casimir(lam, 3), C32=pp_casimir(lamp, 3),
               C3m=pp_casimir(mu, 3), C41=pp_casimir(lam, 4), C42=pp_casimir(lamp, 4), C4m=pp_casimir(mu, 4))
    s = Fr(0)
    for c, mo in _CP_TERMS:
        v = c
        for g in mo:
            v *= ctx[g]
        s += v
    return s


def closed_scalar_check(lam, lamp, verbose=True):
    """compare the computed scalars c, c' (from the reduced operators, every irrep) with the closed forms c_closed / cprime_closed"""
    lam, lamp = tuple(lam), tuple(lamp)
    R = gr.GLNRep.of(lam)
    rs = gr.racah_speiser(lam, lamp)
    n = bad = 0
    for mu in sorted(rs, reverse=True):
        prv = se.prv_space(lam, lamp, mu)
        mats = {nm: vred_matrix(R, reduced_op(nm, lamp), prv) for nm in ('V3', 'V4a', 'V4b', 'S31', 'S41')}
        c3, ok3, c4, ok4 = relations(mats, prv.m)
        n += 1
        bad += (c3 != c_closed(lam, lamp, mu)) + (c4 != cprime_closed(lam, lamp, mu)) + (not ok3) + (not ok4)
    if verbose:
        print('   %-22s %3d irreps: closed forms c_closed, cprime_closed differ from the computed scalars in %d (irrep, scalar) checks' % (prod_str(lam, lamp), n, bad), flush=True)
    return n, bad


def direct_scale(lam, lamp, names=TRIPLE, verbose=True):
    """the direct route at scale: every irrep of lam (x) lam', every PRV basis vector: exact HW vector (G1's exact recursion), the three full two-site operators (G2's Op2) applied exactly, only the
    slice components read off (SliceOp), versus the slice-reduced operators.  Returns (#vectors x operators, #mismatches)"""
    lam, lamp = tuple(lam), tuple(lamp)
    rs = gr.racah_speiser(lam, lamp)
    t0 = time.perf_counter()
    nv = nb = 0
    for mu in sorted(rs, reverse=True):
        prv = se.prv_space(lam, lamp, mu)
        res, _, _ = direct_sector(lam, lamp, mu, names, prv=prv, full=False)
        nv += prv.m * len(names)
        nb += sum(v[0] for v in res.values())
    if verbose:
        print('   %-22s %3d irreps, %4d PRV vectors x %d operators: slice read-off of the FULL two-site operators on the exact HW vectors vs the reduced operators: mismatches %d  (%.1fs)' %
              (prod_str(lam, lamp), len(rs), nv // len(names), len(names), nb, time.perf_counter() - t0), flush=True)
    return nv, nb


# G2 section 6c, theta_2 column: swap(hat S_{k,j}) = a V3 + b V4a + c V4b + scalar  (hat = t^(k-j) S_{k,j}; the theta_2 coefficient is S^{(2)}_{k,j} = (-t)^(j-k) swap(hat))
SWAPPED_G2 = [('S_3_1', (2, 0, 0)), ('S_3_2', (-1, 0, 0)), ('S_4_1', (0, -3, -3)), ('S_4_2', (0, 3, 2)), ('S_4_3', (0, -2, -1))]


def swapped_check(lam, lamp, verbose=True):
    """the theta_2 Laurent coefficients (sites exchanged: the former site-1 words, products of up to three generators in cdet column order, now stand on site 2; the theta_1 data have at most
    two site-2 factors in V3, V4a, V4b):
    (1) the reduced swapped operators preserve M_mu and agree EXACTLY with the slice read-off of the FULL swapped two-site operators on the exact HW vectors (G2's Op2, exact);
    (2) G2's relations (section 6c) swap(hat S_{k,j}) = a V3 + b V4a + c V4b + scalar, (k,j) = (3,1), (3,2), (4,1), (4,2), (4,3), hold on every M_mu with the integers of G2's table;
    (3) the same test with the naive rule (no normal ordering): counts of wrong slices (0 expected: the cdet word order never puts a raising factor left of a lowering factor it does not commute with)."""
    lam, lamp = tuple(lam), tuple(lamp)
    N = len(lam)
    R = gr.GLNRep.of(lam)
    rs = gr.racah_speiser(lam, lamp)
    st = Counter()
    names = [n for n, _ in SWAPPED_G2]
    for mu in sorted(rs, reverse=True):
        prv = se.prv_space(lam, lamp, mu)
        m = prv.m
        st['sectors'] += 1
        st['vectors'] += m
        base = {n: vred_matrix(R, reduced_op(n, lamp), prv) for n in TRIPLE}
        psis = [se.hw_from_slice_exact(lam, lamp, mu, prv.basis[k])[0] for k in range(m)]
        for n, (a, b, c) in SWAPPED_G2:
            try:
                A = vred_matrix(R, reduced_op(n + '_sw', lamp), prv)
            except (ValueError, KeyError):
                st['leave_M'] += 1
                continue
            sop = slice_op(n + '_sw', lam, lamp)
            for j in range(m):
                img = {}
                for i in range(m):
                    if A[i][j] != 0:
                        for cc, x in prv.vec(i).items():
                            img[cc] = img.get(cc, 0) + A[i][j] * x
                st['direct_bad'] += (sop.apply(by_q(psis[j])) != {cc: x for cc, x in img.items() if x != 0})
                st['direct_n'] += 1
            D = mat_lin([(1, A), (-a, base['V3']), (-b, base['V4a']), (-c, base['V4b'])], m)
            ok, sc = mat_scalar(D)
            st['rel_bad'] += (not ok)
            st['rel_n'] += 1
    # the naive rule on the same operators
    sens = variant_sensitivity(lam, lamp, variants=('naive_order',), names=tuple(n + '_sw' for n in names))['naive_order']
    wrong = sum(v[1] for v in sens.values())
    tot = sum(v[2] for v in sens.values())
    if verbose:
        print('   %-20s %3d irreps, %3d PRV vectors | theta_2 operators S_{3,1}, S_{3,2}, S_{4,1}, S_{4,2}, S_{4,3} (swapped): image left M_mu %d times; slice vs FULL swapped two-site operator (exact): %d mismatches of %d; '
              'G2 relations swap(hat) = aV3 + bV4a + cV4b + scalar: %d failures of %d | naive rule (no normal ordering): %d of %d slices wrong' %
              (prod_str(lam, lamp), st['sectors'], st['vectors'], st['leave_M'], st['direct_bad'], st['direct_n'], st['rel_bad'], st['rel_n'], wrong, tot), flush=True)
    return st, wrong, tot


def stage_direct(extras=True):
    hdr('G3.3  DIRECT ROUTE on gl(4): the full two-site operators on the reconstructed exact HW vectors, slice read off  vs  the slice-reduced operators (exact)')
    results = {}
    for lam, lamp in G33_PRODUCTS:
        results[(lam, lamp)] = run_direct_product(lam, lamp)
    print('\n  scalars c, c\' of  S31 + 2 V3 = c 1  (t^2 S_{3,1} = -2 V3 + c)  and  S41 + 3 (V4a + V4b) = c\' 1  (t^3 S_{4,1} = -3 (V4a + V4b) + c\')  per irrep (m):')
    for (lam, lamp), r in results.items():
        print('   %s: ' % prod_str(lam, lamp) + ';  '.join('%s (m=%d): c=%s, c\'=%s' % (fmu(mu), m, fq(c3), fq(c4)) for mu, m, c3, c4 in r['cs']))
    if extras:
        print('\n  closed forms of the scalars (exact polynomial fits in n1 = |lam|, n2 = |lam\'| and the Perelomov-Popov chain Casimirs C2, C3, C4 of lam, lam\', mu; fitted on random small products,')
        print('  see RESULTS-G3; C_r = Perelomov-Popov chain Casimir sum_i l_i^r prod_{j != i} (1 - 1/(l_i - l_j)), l_i = mu_i + 4 - i; C2x = C2(lam) (x = 1), C2(lam\') (x = 2), C2(mu) (x = m), likewise C3x, C4x):')
        print('   c  = 1/2 n1^2 n2 + 7/3 n1 n2 + 1/2 n1 n2^2 + 1/2 n1 C21 + 1/2 n2 C22 + 5/3 (C21 + C22) - 1/3 (C31 + C32) - (5/3 + (n1 + n2)/2) C2m + 1/3 C3m')
        print("   c' = " + ' '.join('%s%s*%s' % ('+' if cc > 0 else '-', fq(abs(cc)), '*'.join(mo)) for cc, mo in _CP_TERMS))
        print('  check against the computed scalars on the six products, and on larger products NOT used for the fits (dim x dim\' > 2500):')
        tot = [0, 0]
        for lam, lamp in G33_PRODUCTS + [((4, 2, 1, 0), (4, 2, 1, 0)), ((3, 2, 1, 0), (4, 2, 1, 0)), ((4, 3, 1, 0), (4, 3, 1, 0)), ((5, 3, 2, 0), (5, 3, 2, 0)), ((6, 4, 2, 0), (6, 4, 2, 0)),
                                         ((5, 2, 1, 1), (3, 2, 2, 1))]:
            a, b = closed_scalar_check(lam, lamp)
            tot[0] += a
            tot[1] += b
        print('  closed-form scalar check total: %d irreps, %d mismatches' % tuple(tot))
        print('\n  the direct route at scale (exact; slice read-off of the full two-site operators V3, V4a, V4b of G2 on the exact HW vectors of the PRV basis):')
        for lam, lamp in [((4, 2, 0, 0), (4, 2, 0, 0)), ((4, 2, 1, 0), (4, 2, 1, 0)), ((3, 2, 1, 0), (4, 2, 1, 0)), ((4, 2, 1, 0), (3, 2, 1, 0)), ((6, 3, 0, 0), (6, 3, 0, 0)),
                          ((6, 4, 2, 0), (6, 4, 2, 0))] + [((2, 1, 0, 0), (1, 1, 0, 0)), ((2, 1, 1, 0), (2, 2, 1, 0)), ((2, 1, 1, 0), (3, 2, 0, 0)), ((2, 1, 0, 0), (3, 2, 1, 0)),
                                                           ((3, 2, 1, 0), (2, 1, 0, 0))]:
            direct_scale(lam, lamp)
        print('\n  the theta_2 Laurent coefficients (sites exchanged: the former site-1 words, up to THREE generators in cdet column order, now on site 2):')
        for lam, lamp in G33_PRODUCTS + [((4, 2, 1, 0), (4, 2, 1, 0)), ((3, 2, 1, 0), (4, 2, 1, 0)), ((6, 4, 2, 0), (6, 4, 2, 0))]:
            swapped_check(lam, lamp)
        print('\n  float direct route (G1 float recursion + G2 float operators on the whole two-site vector; max over irreps of ||slice(O Psi_j) - V_red b_j|| / max(||V_red b_j||, ||b_j||), orthonormal basis of site 1):')
        for lam, lamp in [((4, 2, 0, 0), (4, 2, 0, 0)), ((4, 2, 1, 0), (4, 2, 1, 0))]:
            rs = gr.racah_speiser(lam, lamp)
            worst = {n: 0.0 for n in TRIPLE}
            t0 = time.perf_counter()
            for mu in sorted(rs, reverse=True):
                r = direct_sector_float(lam, lamp, mu)
                for n in TRIPLE:
                    worst[n] = max(worst[n], r[n])
            print('   %-22s %3d irreps: worst deviation V3 %.1e, V4a %.1e, V4b %.1e  (%.1fs)' % (prod_str(lam, lamp), len(rs), worst['V3'], worst['V4a'], worst['V4b'], time.perf_counter() - t0), flush=True)
    return results


# ======================================================================================== 8. G3.4: checks of the canonical copies
def prv_gram_exact(lam, lamp, mu, prv, psis=None):
    """exact G (x) G Gram matrix of the full HW vectors of the PRV basis (invariant form: diagonal GT metrics of the two sites) and the exact vectors"""
    if psis is None:
        psis = [se.hw_from_slice_exact(lam, lamp, mu, prv.basis[k])[0] for k in range(prv.m)]
    g1, g2 = gr.GLNRep.of(lam).g, gr.GLNRep.of(lamp).g
    m = prv.m
    G = [[Fr(0)] * m for _ in range(m)]
    for i in range(m):
        for j in range(i, m):
            a, b = psis[i], psis[j]
            if len(a) > len(b):
                a, b = b, a
            s = Fr(0)
            for (nu, q), v in a.items():
                if v != 0:
                    w = b.get((nu, q))
                    if w is not None and w != 0:
                        s += g1[nu] * g2[q] * v * w
            G[i][j] = G[j][i] = s
    return G, psis


def selfadjoint_bad(mats, G):
    """number of operators A in `mats` (dict) with A^T G != G A (exact; PRV coordinates, column j = image of b_j)"""
    bad = 0
    for n, A in mats.items():
        if mat_mul(go.mat_fr_transpose(A), G) != mat_mul(G, A):
            bad += 1
    return bad


def gram_offdiag_rel(Gm):
    """max_{k != l} |G_kl| / sqrt(G_kk G_ll) of a symmetric (m x m) array/list of floats or mpf"""
    m = len(Gm)
    worst = 0
    for k in range(m):
        for l in range(k + 1, m):
            worst = max(worst, abs(Gm[k][l]) / mp.sqrt(abs(Gm[k][k]) * abs(Gm[l][l])))
    return worst


def float_gram_of_copies(lam, lamp, mu, prv, cps):
    """float G (x) G Gram of the copies' full vectors from the FLOAT recursion applied to the float slices; returns (Gram (m x m), HW object)"""
    F = copy_float_slices(prv, cps)
    hw = se.hw_from_slice(lam, lamp, mu, F)
    V = hw.vals_hat
    return V.T @ V, hw


def exchange_analysis(lam, mu, cps, mats, hw=None, tol_digits=25):
    """lam = lam': the site exchange P (nu, q) -> (q, nu) preserves the HW space and (G2 section 6e) conjugates P V3 P = -V3 + c, P V4a P = -2 V4a - V4b + c',
    P V4b P = 3 V4a + 2 V4b - c' (c, c' fixed by the traces).  Eigenvalue side: the affine map T(e3, e4a, e4b) = (c - e3, -2 e4a - e4b + c', 3 e4a + 2 e4b - c') must permute the joint
    eigenvalue triples (partner sigma(k)); with the float full vectors (hw): overlap matrix <Psi_j | P Psi_k> / norms must be the signed permutation matrix of sigma.
    Returns dict(c, cp, spectrum_symmetric_dev, T_ok, sigma, nfixed, npairs, parity, overlap_dev)"""
    m = len(cps)
    old = mp.mp.dps
    mp.mp.dps = 50
    try:
        t3 = mpq(mat_trace(mats['V3']))
        t4a = mpq(mat_trace(mats['V4a']))
        t4b = mpq(mat_trace(mats['V4b']))
        c = 2 * t3 / m
        cp = (3 * t4a + t4b) / m
        trips = [cp_['ev'] for cp_ in cps]
        ev3 = [t[0] for t in trips]
        sym_dev = max(abs(ev3[k] + ev3[m - 1 - k] - c) for k in range(m))
        tol = mp.mpf(10) ** (-tol_digits)
        sigma = []
        Tdev = mp.mpf(0)
        for k in range(m):
            e3, e4a, e4b = trips[k]
            img = (c - e3, -2 * e4a - e4b + cp, 3 * e4a + 2 * e4b - cp)
            best = min(range(m), key=lambda j: max(abs(a - b) for a, b in zip(img, trips[j])))
            dev = max(abs(a - b) for a, b in zip(img, trips[best]))
            Tdev = max(Tdev, dev)
            sigma.append(best)
        T_ok = (Tdev < tol) and sorted(sigma) == list(range(m))
        nfixed = sum(1 for k in range(m) if sigma[k] == k)
        npairs = sum(1 for k in range(m) if sigma[k] > k)
        v3deg = any(abs(ev3[k + 1] - ev3[k]) < mp.mpf(10) ** -20 for k in range(m - 1))
    finally:
        mp.mp.dps = old
    out = dict(c=c, cp=cp, sym_dev=sym_dev, T_ok=T_ok, Tdev=Tdev, sigma=sigma, nfixed=nfixed, npairs=npairs, v3_degenerate=v3deg)
    if hw is not None:
        n1, n2 = hw.dim1, hw.dim2
        Xs = [sps.csr_matrix((hw.vals_hat[:, k], (hw.rows, hw.cols)), shape=(n1, n2)) for k in range(m)]
        nrm = [math.sqrt(float(X.multiply(X).sum())) for X in Xs]
        O = np.zeros((m, m))
        for j in range(m):
            for k in range(m):
                O[j, k] = float(Xs[j].multiply(Xs[k].T).sum()) / (nrm[j] * nrm[k])
        pred = np.zeros((m, m))
        for k in range(m):
            pred[sigma[k], k] = 1.0
        dev = float(np.max(np.abs(np.abs(O) - pred)))
        par = []
        for k in range(m):
            if sigma[k] == k:
                par.append(1 if O[k, k] > 0 else -1)
        out.update(overlap_dev=dev, parity=par, O=O, forced=(abs(sum(par)) >= 2))
    return out


def copies_report(lam, lamp, mu, dps=50, g2=False, exch=False):
    """all checks of G3.4 for one sector; returns a dict"""
    lam, lamp, mu = tuple(lam), tuple(lamp), tuple(mu)
    prv = se.prv_space(lam, lamp, mu)
    _, mats = triple_mats(lam, lamp, mu, prv)
    cps = copies(lam, lamp, mu, dps=dps, prv=prv, mats=mats)
    m = prv.m
    G, psis = prv_gram_exact(lam, lamp, mu, prv)
    sab = selfadjoint_bad(mats, G)
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        # orthogonality from the exact Gram of the PRV basis and the mp coordinates of the copies
        Gc = [[mp.fsum(cps[k]['coords'][i] * mpq(G[i][j]) * cps[l]['coords'][j] for i in range(m) for j in range(m)) for l in range(m)] for k in range(m)]
        orth_exact = gram_offdiag_rel(Gc) if m > 1 else mp.mpf(0)
        # orthogonality from the float recursion
        Gf, hw = float_gram_of_copies(lam, lamp, mu, prv, cps)
        orth_float = float(gram_offdiag_rel(Gf.tolist())) if m > 1 else 0.0
        res_hw, wdev = se.hw_residual(hw, extra=True)
    finally:
        mp.mp.dps = old
    out = dict(mu=mu, m=m, cps=cps, prv=prv, mats=mats, G=G, selfadj_bad=sab, orth_exact=orth_exact, orth_float=orth_float, hwres=res_hw, wdev=wdev, info=cps.info)
    if g2:
        ref = g2_sector(lam, lamp, mu)
        key = tol_cmp(mp.mpf(10) ** (-(dps // 2)))
        dtrip = mp.mpf(0)
        dsl = mp.mpf(0)
        for cp_, t2, s2 in zip(cps, ref['trips'], ref['slices']):
            dtrip = max(dtrip, max(abs(a - b) for a, b in zip(cp_['ev'], t2)))
            keys = set(cp_['slice']) | set(s2)
            dsl = max(dsl, max(abs(cp_['slice'].get(kk, 0) - s2.get(kk, 0)) for kk in keys))
        out['g2_dtrip'] = dtrip
        out['g2_dslice'] = dsl
    if exch and lam == lamp:
        ex = exchange_analysis(lam, mu, cps, mats, hw=hw)
        ex['c_matches_S31'] = bool(abs(ex['c'] - mpq(c_closed(lam, lamp, mu))) < mp.mpf(10) ** -30 and abs(ex['cp'] - mpq(cprime_closed(lam, lamp, mu))) < mp.mpf(10) ** -30)
        out['exch'] = ex
    return out


def run_copies_product(lam, lamp, dps=50, g2=False, exch=False, listing=True, verbose=True, min_m_list=2, maxlines=60):
    lam, lamp = tuple(lam), tuple(lamp)
    rs = gr.racah_speiser(lam, lamp)
    t0 = time.perf_counter()
    S = Counter()
    worst = dict(offmax=mp.mpf(0), imag=mp.mpf(0), orth_exact=mp.mpf(0), orth_float=0.0, hwres=0.0, dtrip=mp.mpf(0), dslice=mp.mpf(0), exch=0.0, Tdev=mp.mpf(0))
    jgap = mp.inf
    lines = []
    exch_rows = []
    for mu in sorted(rs, reverse=True):
        r = copies_report(lam, lamp, mu, dps=dps, g2=g2, exch=exch)
        m = r['m']
        S['sectors'] += 1
        S['copies'] += m
        S['diag_fail'] += (not r['info']['ok'])
        S['selfadj_bad'] += r['selfadj_bad']
        worst['offmax'] = max(worst['offmax'], r['info']['offmax'])
        worst['imag'] = max(worst['imag'], r['info']['imag'])
        worst['orth_exact'] = max(worst['orth_exact'], r['orth_exact'])
        worst['orth_float'] = max(worst['orth_float'], r['orth_float'])
        worst['hwres'] = max(worst['hwres'], r['hwres'])
        jgap = min(jgap, r['info']['jgap'])
        if g2:
            worst['dtrip'] = max(worst['dtrip'], r['g2_dtrip'])
            worst['dslice'] = max(worst['dslice'], r['g2_dslice'])
        if 'exch' in r:
            ex = r['exch']
            S['exch_sectors'] += 1
            S['exch_T_bad'] += (not ex['T_ok'])
            S['exch_fixed'] += ex['nfixed']
            S['exch_pairs'] += ex['npairs']
            S['exch_c_match'] += ex['c_matches_S31']
            worst['exch'] = max(worst['exch'], ex['overlap_dev'])
            worst['Tdev'] = max(worst['Tdev'], ex['Tdev'])
            if m >= 2:
                exch_rows.append((mu, m, ex))
                S['exch_ms2'] += 1
                S['exch_forced'] += ex['forced']
                S['exch_v3deg'] += ex['v3_degenerate']
                S['exch_accidental'] += (ex['v3_degenerate'] and not ex['forced'])
                S['exch_forced_nondeg'] += (ex['forced'] and not ex['v3_degenerate'])
        if listing and m >= min_m_list:
            lines.append('    mu=%-14s m=%2d | joint-diag offmax %s imag %s (dps %d) | selfadj viol %d | G-orth: exact-Gram %s, float %s | HW resid %.1e | joint gap %s%s' %
                         (str(mu), m, mp.nstr(r['info']['offmax'], 3), mp.nstr(r['info']['imag'], 3), r['info']['dps'], r['selfadj_bad'], mp.nstr(r['orth_exact'], 3), '%.1e' % r['orth_float'],
                          r['hwres'], mp.nstr(r['info']['jgap'], 4), (' | vs G2: triples %s, slices %s' % (mp.nstr(r['g2_dtrip'], 3), mp.nstr(r['g2_dslice'], 3))) if g2 else ''))
    if verbose:
        print('  %-20s sectors %3d, copies %4d | joint diag failed: %d, worst off-diagonal %s, worst imaginary part %s, smallest joint eigenvalue distance %s | self-adjointness violations (exact, G-Gram of the PRV basis): %d | '
              'G-orthogonality: exact-Gram %s, float full vectors %.1e | worst HW residual %.1e%s%s  (%.1fs)' %
              (prod_str(lam, lamp), S['sectors'], S['copies'], S['diag_fail'], mp.nstr(worst['offmax'], 3), mp.nstr(worst['imag'], 3), mp.nstr(jgap, 4), S['selfadj_bad'], mp.nstr(worst['orth_exact'], 3),
               worst['orth_float'], worst['hwres'],
               (' | vs G2 joint eigenbasis: max triple diff %s, max slice diff %s' % (mp.nstr(worst['dtrip'], 3), mp.nstr(worst['dslice'], 3))) if g2 else '',
               (' | exchange: %d sectors, T-map failures %d, fixed copies %d, 2-cycles %d, max |overlap - signed permutation| %.1e, exchange scalars (c, c\') = the closed-form scalars of the S31 / S41 relations on %d of %d sectors; of the %d sectors with m >= 2: |p+ - p-| >= 2 (V3 degeneracy forced): %d, V3 degenerate: %d (not forced: %d), forced but V3 simple: %d' %
                (S['exch_sectors'], S['exch_T_bad'], S['exch_fixed'], S['exch_pairs'], worst['exch'], S['exch_c_match'], S['exch_sectors'], S['exch_ms2'], S['exch_forced'], S['exch_v3deg'], S['exch_accidental'], S['exch_forced_nondeg'])) if exch and lam == lamp else '',
               time.perf_counter() - t0), flush=True)
        if listing:
            for ln in lines[:maxlines]:
                print(ln, flush=True)
    return dict(S=S, worst=worst, jgap=jgap, exch_rows=exch_rows)


def fmt_ev(t, n=8):
    return '(' + ', '.join(mp.nstr(x, n) for x in t) + ')'


def exchange_table(rows, title):
    """tabulate the exchange structure: per sector with m >= 2: c, c', symmetry of the V3 spectrum about c/2, T-map, fixed copies with their P-parity, 2-cycles"""
    print('  exchange structure (lam = lam\'), %s: P V3 P = -V3 + c,  P V4a P = -2 V4a - V4b + c\',  P V4b P = 3 V4a + 2 V4b - c\'' % title)
    for mu, m, ex in rows:
        par = ex.get('parity', [])
        cyc = [(k + 1, s + 1) for k, s in enumerate(ex['sigma']) if s > k]
        fixed = [k + 1 for k, s in enumerate(ex['sigma']) if s == k]
        print('    mu=%-14s m=%2d  c = %s (spectrum of V3 symmetric about c/2 to %s), c\' = %s | T permutes the joint triples: %s | fixed copies %s (P-parities %s; p+ - p- = %d) | 2-cycles %s | max |overlap - signed perm| %.1e' %
              (str(mu), m, mp.nstr(ex['c'], 12), mp.nstr(ex['sym_dev'], 2), mp.nstr(ex['cp'], 12), ex['T_ok'], fixed, ''.join('+' if p > 0 else '-' for p in par) or '-', sum(par), cyc, ex.get('overlap_dev', float('nan'))))


def stage_copies(extras=True):
    hdr('G3.4  CANONICAL COPIES: joint eigenvectors of (V_red,3, V_red,4a, V_red,4b) on M_mu; G-orthogonality, self-adjointness, agreement with G2\'s eigenbasis, exchange structure')
    out = {}
    listed = [((3, 2, 1, 0), (3, 2, 1, 0)), ((2, 1, 0, 0), (2, 1, 0, 0))]
    for lam, lamp in G33_PRODUCTS:
        out[(lam, lamp)] = run_copies_product(lam, lamp, g2=True, exch=(lam == lamp), listing=((lam, lamp) in listed))
    print()
    for lam, lamp in G33_PRODUCTS:
        r = out[(lam, lamp)]
        if lam == lamp and r['exch_rows']:
            exchange_table(r['exch_rows'], prod_str(lam, lamp))
    if extras:
        print('\n  larger products (copies, G-orthogonality, self-adjointness, exchange; no G2 comparison):')
        for lam, lamp in [((4, 2, 0, 0), (4, 2, 0, 0)), ((4, 2, 1, 0), (4, 2, 1, 0)), ((6, 4, 2, 0), (6, 4, 2, 0))]:
            out[(lam, lamp)] = run_copies_product(lam, lamp, g2=False, exch=True, listing=False)
        for lam, lamp in [((4, 2, 0, 0), (4, 2, 0, 0)), ((6, 4, 2, 0), (6, 4, 2, 0))]:
            r = out[(lam, lamp)]
            exchange_table([x for x in r['exch_rows'] if x[1] >= (6 if lam[0] == 6 else 2)], prod_str(lam, lamp) + (' (m >= 6 only)' if lam[0] == 6 else ''))
    return out


# ======================================================================================== 9. G3.5: timing of the n_W-free part vs the full reconstruction
def bench_one(lam, lamp=None, reps=3, do_b=True, dps=50, resid=True):
    """G3.5 for the product lam (x) lamp, min over `reps` cold runs (all caches cleared before each run, single process).
    t_a = the n_W-free part: GT generators of both sites (adjacent + diagonal: t_gen) + the 6 non-adjacent generators of site 1 that the reduced words use (t_gen_na) + Racah-Speiser (t_rs)
          + exact PRV spaces of all irreps (t_prv)  [= G1's t_a, plus t_gen_na]  +  the reduction of the three operators for this lam' (t_red: normal ordering + evaluation at lam', once
          per lam')  +  the three reduced operators applied exactly to the PRV basis of every M_mu (t_ops)  +  the joint eigen-decomposition (mpmath, dps digits) and the slices of all copies (t_eig).
    t_b = t_a + the float setup of the level recursion (t_setup) + the float reconstruction of EVERY full copy, batched per irrep (t_recon).
    resid: untimed pass for the worst HW residual of the reconstructed copies and the diagnostics of the joint diagonalisation (taken from the best run)."""
    lam = tuple(lam)
    lamp = lam if lamp is None else tuple(lamp)
    best = None
    for _ in range(reps):
        clear_caches()
        T = {}
        t0 = time.perf_counter()
        R1 = gr.GLNRep.of(lam)
        R2 = gr.GLNRep.of(lamp)
        T['gen'] = time.perf_counter() - t0
        t0 = time.perf_counter()
        for ij in NONADJ:
            R1.E[ij]
        T['gen_na'] = time.perf_counter() - t0
        t0 = time.perf_counter()
        rs = gr.racah_speiser(lam, lamp)
        mus = sorted(rs, reverse=True)
        T['rs'] = time.perf_counter() - t0
        t0 = time.perf_counter()
        prvs = {mu: se.prv_space(lam, lamp, mu) for mu in mus}
        T['prv'] = time.perf_counter() - t0
        t0 = time.perf_counter()
        ops = [reduced_op(n, lamp) for n in TRIPLE]
        T['red'] = time.perf_counter() - t0
        t0 = time.perf_counter()
        allm = {mu: {n: vred_matrix(R1, op, prvs[mu]) for n, op in zip(TRIPLE, ops)} for mu in mus}
        T['ops'] = time.perf_counter() - t0
        t0 = time.perf_counter()
        cps = {mu: copies(lam, lamp, mu, dps=dps, prv=prvs[mu], mats=allm[mu]) for mu in mus}
        T['eig'] = time.perf_counter() - t0
        T['a'] = sum(T[k] for k in ('gen', 'gen_na', 'rs', 'prv', 'red', 'ops', 'eig'))
        T['a_g1'] = T['gen'] + T['rs'] + T['prv']
        # supplementary: the same joint eigen-decomposition in float64 (numpy) instead of mpmath: copies' slices as (w, m) float arrays, not normalised
        t0 = time.perf_counter()
        fl = {}
        for mu in mus:
            jf = joint_diag_float([allm[mu][n] for n in TRIPLE])
            fl[mu] = (jf, prvs[mu].float_basis() @ jf['U'])
        T['eig_f'] = time.perf_counter() - t0
        T['a_f'] = T['a'] - T['eig'] + T['eig_f']
        if do_b:
            t0 = time.perf_counter()
            rec = se.Recursion.of(lam, lamp)
            T['setup'] = time.perf_counter() - t0
            t0 = time.perf_counter()
            nent = 0
            for mu in mus:
                F = copy_float_slices(prvs[mu], cps[mu])
                hw = se.hw_from_slice(lam, lamp, mu, F)
                nent += hw.nW * hw.m
            T['recon'] = time.perf_counter() - t0
            T['nent'] = nent
            T['b'] = T['a'] + T['setup'] + T['recon']
            T['b_f'] = T['a_f'] + T['setup'] + T['recon']
        else:
            T['b'] = None
            T['b_f'] = None
        if best is None or (T['b'] if do_b else T['a']) < (best['b'] if do_b else best['a']):
            best = T
            best_cps, best_prvs, best_mus, best_fl = cps, prvs, mus, fl
    # untimed diagnostics
    R1 = gr.GLNRep.of(lam)
    worst, wmu = 0.0, None
    if resid and do_b:
        for mu in best_mus:
            hw = se.hw_from_slice(lam, lamp, mu, copy_float_slices(best_prvs[mu], best_cps[mu]))
            r, wd = se.hw_residual(hw, extra=(R1.dim <= 400))
            if r > worst:
                worst, wmu = r, mu
    # accuracy of the float64 joint eigen-decomposition against the 50-digit one (sorted lexicographically; joint eigenvalue distances are >= 6 here)
    dev_f = 0.0
    off_f = 0.0
    for mu in best_mus:
        jf = best_fl[mu][0]
        tr = jf['trip']
        trr = np.round(tr, 6)                                     # float ties (V3 degenerate sectors) must not be broken by rounding noise
        order = np.lexsort(tuple(trr[:, t] for t in reversed(range(tr.shape[1]))))
        ref = np.array([[float(x) for x in cp['ev']] for cp in best_cps[mu]])
        dev_f = max(dev_f, float(np.abs(tr[order] - ref).max() / max(1.0, np.abs(ref).max())))
        off_f = max(off_f, jf['offmax'])
    offmax = max(float(c.info['offmax']) for c in best_cps.values())
    diag_fail = sum(1 for c in best_cps.values() if not c.info['ok'])
    dps_raised = sum(1 for c in best_cps.values() if c.info['dps'] > dps)
    jg = min((c.info['jgap'] for c in best_cps.values() if len(c) > 1), default=mp.inf)
    best.update(lam=lam, lamp=lamp, dim=R1.dim, nirrep=len(rs), ncopies=sum(rs.values()), maxm=max(rs.values()), maxw=max(p.w for p in best_prvs.values()),
                maxnW=max(gr.nW_of(lam, lamp, mu) for mu in best_mus), resid=worst, worst_mu=wmu, offmax=offmax, diag_fail=diag_fail, dps_raised=dps_raised, jgap=jg,
                nwords=[op.nwords for op in ops], nnodes=[op.nnodes for op in ops], ta=best['a'], tb=best['b'], reps=reps, dev_f=dev_f, off_f=off_f,
                sep=separation_counts(best_cps))
    return best


# [2a,a,0,0]^2: the scaling family (m = a + 1);  [3a,2a,a,0]^2: the stress family (max m = 4, 12, 27, 51);  [S,S,0,0]^2: multiplicity-free
BENCH_FAMILIES_G3 = [('[2a,a,0,0]^2', lambda a: ((2 * a, a, 0, 0),) * 2, range(1, 7)),
                     ('[3a,2a,a,0]^2', lambda a: ((3 * a, 2 * a, a, 0),) * 2, range(1, 4))]


def bench_reps(dim):
    return 5 if dim <= 500 else (3 if dim <= 5000 else 2)


def bench_family_rows(name, f, rng_, verbose=True, reps=None, do_b=True):
    rows = []
    if verbose:
        print('  family %s' % name)
        print('   %2s %-14s %6s %6s %6s %5s %5s %8s | %8s | %8s | %11s | %8s | breakdown t_gen / t_gen_na / t_rs / t_prv / t_red / t_ops / t_eig | t_setup / t_recon | words (V3,V4a,V4b) | worst joint-diag offdiag, min joint dist |'
              ' float64 copies: t_eig_f, t_a_f, t_b_f, max |ev - ev_mp|, max offdiag' %
              ('a', 'lam', 'dimV', 'irreps', 'copies', 'max m', 'max w', 'max nW', 't_a', 't_b', 'entries', 'resid'), flush=True)
    for a in rng_:
        lam, lamp = f(a)
        r = reps or bench_reps(gr.weyl_dim(lam))
        T = bench_one(lam, lamp, reps=r, do_b=do_b)
        T['family'], T['a'] = name, a
        rows.append(T)
        if verbose:
            print('   %2d %-14s %6d %6d %6d %5d %5d %8d | %8.4f | %8s | %11s | %8.1e | %.4f / %.4f / %.4f / %.4f / %.4f / %.4f / %.4f | %s | %s | %.1e, %s | %.4f, %.4f, %s, %.1e, %.1e   [min of %d cold runs%s]' %
                  (a, lam_str(T['lam']), T['dim'], T['nirrep'], T['ncopies'], T['maxm'], T['maxw'], T['maxnW'], T['ta'], ('%.4f' % T['tb']) if do_b else '-', T['nent'] if do_b else '-', T['resid'],
                   T['gen'], T['gen_na'], T['rs'], T['prv'], T['red'], T['ops'], T['eig'], ('%.4f / %.4f' % (T['setup'], T['recon'])) if do_b else '-', '/'.join(str(x) for x in T['nwords']),
                   T['offmax'], mp.nstr(T['jgap'], 4), T['eig_f'], T['a_f'], ('%.4f' % T['b_f']) if do_b else '-', T['dev_f'], T['off_f'], T['reps'],
                   '' if not T['diag_fail'] and not T['dps_raised'] else '; diag failures %d, precision raised in %d sectors' % (T['diag_fail'], T['dps_raised'])), flush=True)
    if verbose:
        print('   separation of the copies by subsets of the commuting triple (sectors with m >= 2; number of sectors where the spectrum of the subset is DEGENERATE, 50-digit eigenvalues, tolerance 1e-20):')
        for T in rows:
            n, d = T['sep']
            print('     a = %d %-14s %3d sectors with m >= 2: ' % (T['a'], lam_str(T['lam']), n) + ', '.join('{%s}: %d' % (k, v) for k, v in d.items()), flush=True)
    return rows


def rect_scalars(S, verbose=True):
    """[S,S,0,0]^2 is multiplicity-free: the reduced operators are scalars on every M_mu.  Returns [(mu, V3, V4a, V4b)] (exact) and prints the V3 table"""
    lam = (S, S, 0, 0)
    rs = gr.racah_speiser(lam, lam)
    rows = []
    for mu in sorted(rs, reverse=True):
        prv, mats = triple_mats(lam, lam, mu)
        assert prv.m == 1
        rows.append((mu, mats['V3'][0][0], mats['V4a'][0][0], mats['V4b'][0][0]))
    if verbose:
        print('  [%d,%d,0,0]^2 (m = 1 everywhere, %d irreps): the scalars (V_red,3, V_red,4a, V_red,4b) on M_mu:' % (S, S, len(rows)))
        print('     ' + '; '.join('%s: (%s, %s, %s)' % (fmu(mu), fq(a), fq(b), fq(c)) for mu, a, b, c in rows))
    return rows


def rect_closed_check(S, verbose=False):
    """closed forms on the multiplicity-free family [S,S,0,0]^2 (m = 1: the reduced operators are scalars): with C_r = Perelomov-Popov chain Casimirs and top = (2S,2S,0,0),
         V_red,3 = (S+1)/2 (C2(top) - C2(mu)) = 4 S (S+1)^2 - (S+1)/2 C2(mu),   V_red,4a = 0,   V_red,4b = F(mu) - F(top),   F(mu) = (2S+1)/4 C2(mu) - C4(mu)/4 + C2(mu)^2/8.
    Returns (#irreps, #mismatches) against the scalars computed by the reduction."""
    lam = (S, S, 0, 0)
    top = (2 * S, 2 * S, 0, 0)
    C2 = lambda mu: pp_casimir(mu, 2)
    C4 = lambda mu: pp_casimir(mu, 4)
    F = lambda mu: Fr(2 * S + 1, 4) * C2(mu) - C4(mu) / 4 + C2(mu) ** 2 / 8
    rows = rect_scalars(S, verbose=False)
    bad = 0
    for mu, v3, v4a, v4b in rows:
        bad += (v3 != Fr(S + 1, 2) * (C2(top) - C2(mu))) + (v4a != 0) + (v4b != F(mu) - F(top))
    return len(rows), bad


def gl3_cost_reference(rows4, verbose=True):
    """where the reduction cost sits relative to gl(3) (same machine): slice_engine.bench_one(a) for [2a,a,0]^2 (ONE operator, the 24-term V_red, exact, float64 eigen) against the gl(4) rows of the
    scaling family [2a,a,0,0]^2 (THREE operators with 36 / 61 / 67 reduced words, exact; t_ops = their application to every PRV basis vector of every irrep)"""
    s3 = _se3()
    out = []
    if verbose:
        print('\n  cost of the reduction relative to gl(3) (same machine, same cold-run protocol; gl(3): slice_engine.bench_one, ONE operator = the 24-term V_red, exact application + float64 eigen-decomposition;')
        print('  gl(4): t_ops = the THREE reduced operators (V3, V4a, V4b: 36 / 61 / 67 words at lam\' = [2a,a,0,0]) applied exactly to every PRV basis vector of every irrep):')
    for a in (4, 7, 10, 14):
        T = s3.bench_one(a, reps=3)
        out.append(('gl3', a, T))
        if verbose:
            print('   gl(3) [%2d,%2d,0]^2: dim %5d, irreps %3d, copies %4d, max w %2d | t_Vred %.4f s = %4.0f us per copy (%4.0f us per irrep) | t_a %.4f s, t_b %.4f s' %
                  (2 * a, a, T['dim'], T['nirrep'], T['ncopies'], T['maxw'], T['vred'], 1e6 * T['vred'] / T['ncopies'], 1e6 * T['vred'] / T['nirrep'], T['ta'], T['tb']), flush=True)
    for T in rows4:
        out.append(('gl4', T['a'], T))
        if verbose:
            print('   gl(4) [%2d,%2d,0,0]^2: dim %5d, irreps %3d, copies %4d, max w %2d | t_ops %.4f s = %4.0f us per copy for the triple = %4.0f us per operator and copy (%4.0f us per irrep) | t_a %.4f s (t_a_f %.4f s with float64 eig), t_b %.4f s' %
                  (2 * T['a'], T['a'], T['dim'], T['nirrep'], T['ncopies'], T['maxw'], T['ops'], 1e6 * T['ops'] / T['ncopies'], 1e6 * T['ops'] / T['ncopies'] / 3, 1e6 * T['ops'] / T['nirrep'], T['ta'], T['a_f'], T['tb']), flush=True)
    return out


def stage_bench(extended=False):
    hdr('G3.5  TIMING: n_W-free part t_a (PRV spaces + the three reduced operators on every M_mu + joint eigen-decomposition + copies\' slices) vs t_b (+ full reconstruction of every copy)')
    print('  machine: %s; Python %s, numpy %s, scipy %s, mpmath %s; single process; min over cold runs (caches cleared); perf_counter' % (platform.platform(), platform.python_version(), np.__version__,
                                                                                                                       __import__('scipy').__version__, mp.__version__))
    allrows = {}
    for name, f, rng_ in BENCH_FAMILIES_G3:
        allrows[name] = bench_family_rows(name, f, rng_)
    allrows['gl3ref'] = gl3_cost_reference(allrows[BENCH_FAMILIES_G3[0][0]])
    if extended:
        print('\n  extension (slices only, t_a; one cold run each): ')
        allrows['ext1'] = bench_family_rows('[2a,a,0,0]^2 (extension, a = 7, 8)', BENCH_FAMILIES_G3[0][1], range(7, 9), reps=1, do_b=False)
        allrows['ext2'] = bench_family_rows('[3a,2a,a,0]^2 (extension, a = 4)', BENCH_FAMILIES_G3[1][1], range(4, 5), reps=1, do_b=False)
    print('\n  multiplicity-free family [S,S,0,0]^2:')
    for S in (4, 8):
        T = bench_one((S, S, 0, 0), reps=3)
        print('   S = %d: dim V %d, irreps %d, max w %d, max nW %d | t_a %.4f s (gen %.4f, gen_na %.4f, rs %.4f, prv %.4f, red %.4f, ops %.4f, eig %.4f), t_b %.4f s | worst HW resid %.1e' %
              (S, T['dim'], T['nirrep'], T['maxw'], T['maxnW'], T['ta'], T['gen'], T['gen_na'], T['rs'], T['prv'], T['red'], T['ops'], T['eig'], T['tb'], T['resid']), flush=True)
        allrows['rect%d' % S] = T
    allrows['rect_scalars'] = {}
    for S in (4, 8):
        allrows['rect_scalars'][S] = rect_scalars(S)
    tot = [0, 0]
    for S in range(1, 11):
        a, b = rect_closed_check(S)
        tot[0] += a
        tot[1] += b
    print('  closed forms on [S,S,0,0]^2 (m = 1): V_red,3 = (S+1)/2 (C2(top) - C2(mu)), V_red,4a = 0, V_red,4b = F(mu) - F(top), F = (2S+1)/4 C2 - C4/4 + C2^2/8 (Perelomov-Popov chain Casimirs, top = (2S,2S,0,0)):'
          ' checked on S = 1..10: %d irreps, %d mismatches' % tuple(tot))
    return allrows


# ======================================================================================== 10. G3.6: reference files ref/gl4-cg-copies-*.txt
def copies_ref_name(lam, lamp):
    s = ''.join(str(x) for x in lam)
    if tuple(lam) != tuple(lamp):
        s += '-' + ''.join(str(x) for x in lamp)
    return 'gl4-cg-copies-%s.txt' % s


def fd(x, nd=30):
    """decimal with nd significant digits (mpf / mpc real part); exact zero -> '0'"""
    x = mp.mpf(x)
    if x == 0:
        return '0'
    return mp.nstr(x, nd, strip_zeros=False, min_fixed=-8, max_fixed=24)


def f_mat(A):
    return '{%s}' % ','.join('{%s}' % ','.join(fq(x) for x in row) for row in A)


def f_poly(c):
    return '{%s}' % ','.join(fq(x) for x in c)


def _copies_hdr(lam, lamp, has_vec, name, date):
    N = len(lam)
    L = []
    A = L.append
    A('# %s   (WP-G3 of the gl(4) slice-engine programme; generator Paul/Python/gaudin_gl4/slice_ops_gl4.py `python -X utf8 slice_ops_gl4.py ref`; %s)' % (name, date))
    A('# Canonical multiplicity-space data of gl(%d), V_Lambda (x) V_Lambda\' with Lambda = {%s} (REP, dim DIM) and Lambda\' = {%s} (REP2, dim DIM2): for EVERY irrep mu the PRV slice space M_mu (exact basis), the three two-site Gaudin missing-label operators V3, V4a, V4b '
      'REDUCED TO THE SLICE (exact m x m matrices in the PRV basis), their exact characteristic polynomials, the exact G-Gram of the PRV basis, and the canonical copies (joint eigenvectors; 30 digits): slice, eigenvalue triple, HW residual, G-norm^2%s.' %
      (N, ','.join(str(x) for x in lam), ','.join(str(x) for x in lamp), ', full vector (VEC)' if has_vec else ''))
    A('# FORMAT: one record per line, TAB-separated fields, first field = tag, nested braces = lists (Wolfram InputForm style). p/q = exact rationals (integers without denominator); the COPY / GNORM2 / VEC numbers are decimals with 30 significant digits; HWRESID are Python floats. Lines starting with # are comments. '
      'Same conventions as the G1 files ref/gl4-cg-slices-*.txt (GT order, 1-based indices, two sites, slice).')
    A('# ORDER OF THE GT BASES, TWO SITES, SLICE: exactly as in ref/gl4-cg-slices-*.txt (order of gl3_reps.gt_patterns, HW state LAST; site 1 = first factor Lambda (nu, slow), site 2 = Lambda\' (q, fast); flat index (nu-1)*DIM2 + q; slice = Psi[., q_top] in the weight space V_Lambda[mu-Lambda\'] of site 1; '
      'GMET1/GMET2 = invariant diagonal metrics, G[HW] = 1).')
    A('# THE OPERATORS (WP-G2, ref/gl4-gaudin-ops.txt; the t-free triple of the two-site gl(4) Gaudin model with L(u)_ij = e1_ij/(u-theta1) + e2_ij/(u-theta2), t = theta1 - theta2, cdet(d_u - L) = d^4 - S1 d^3 + S2 d^2 - S3 d + S4 column-ordered): '
      'V3 = t S_{3,2}, V4a = t S_{4,3}, V4b = t^2 S_{4,2} (S_{k,j} = Laurent coefficient of S_k at theta1 of order x^-j, x = u - theta1); exactly t-independent, mutually commuting, self-adjoint for G (x) G, simple joint spectrum on every sector tested (G2).')
    A('# THE REDUCTION RULE (this file\'s VRED matrices): a term c W1 W2 of an operator (W1 = word of site-1 generators, W2 = word of site-2 generators) acts on the slice psi_top of a HW vector as follows: normal-order W2 in U(gl_4) as sum c_w L H R (L lowering, H Cartan, R raising factors); '
      'a term with L != 1 contributes 0; otherwise c c_w prod_{E_aa in H} Lambda\'_a (-1)^|R| W1 E_{r_|R|} ... E_{r_1} psi_top for R = e2_{r_1} ... e2_{r_|R|} (raising factors move to site 1 with the sign (-1)^|R|, in REVERSED order, to the right of W1; every site-2 Cartan factor is evaluated to Lambda\'). '
      'So each operator acts on M_mu as a polynomial in the site-1 generators of V_Lambda: no site-2 vector, no n_W-sized object.')
    A('# RECORDS: REP | REP2 | DIM | DIM2 | GTBASIS | GTBASIS2 | WEIGHTS1 | WEIGHTS2 | GMET1 | GMET2 | NSECTORS | NCOPIES | GENERIC | then per sector (descending lexicographic mu): SECTOR, WSPACE, PRVBASIS x m, PRVFREE, PRVGRAM, VRED3, VRED4A, VRED4B, CHARPOLY3, CHARPOLYJ, then per copy k = 1..m: COPY, HWRESID, GNORM2%s |' %
      (', VEC' if has_vec else ''))
    A('#   SECTOR mu m w nW | WSPACE mu {nu,...} (1-based GT indices of V_Lambda[mu-Lambda\'], ascending; slice and PRVBASIS entries are indexed by them) | PRVBASIS mu {{nu,val},...} (m lines, exact canonical RREF basis of M_mu as in the G1 files: basis vector k = 1 at the k-th free position, 0 at the other free positions) |')
    A('#   PRVFREE mu {nu_1,...,nu_m} : the free positions (1-based GT indices): the PRV coordinates u_k of a slice s in M_mu are u_k = s[nu_k] (s = sum_k u_k PRVBASIS_k) | PRVGRAM mu {{g_kl}} : EXACT m x m Gram matrix of the full HW vectors of the PRV basis for the invariant form G (x) G: g_kl = sum g1[nu] g2[q] Psi_k[nu,q] Psi_l[nu,q] |')
    A('#   VRED3 / VRED4A / VRED4B mu {{a_ij}} : EXACT m x m matrices (rows i, columns j) of the slice-reduced V3, V4a, V4b in the PRV basis: V_red PRVBASIS_j = sum_i a_ij PRVBASIS_i (column j = image of the j-th basis vector in PRV coordinates) |')
    A('#   CHARPOLY3 mu {c_0,...,c_m} : exact characteristic polynomial det(x - A3) of VRED3, ASCENDING powers (c_m = 1) | CHARPOLYJ mu {c_0,...,c_m} : the same for the GENERIC rational combination A3 + 7/19 A4a + 11/13 A4b of G2 (the combination used for the eigen-decomposition) | GENERIC {1,7/19,11/13} |')
    A('#   COPY mu k ev3 ev4a ev4b {{nu,val},...} : the k-th canonical copy: joint eigenvector of (VRED3, VRED4A, VRED4B) with eigenvalues (ev3, ev4a, ev4b) (30 digits; copies ordered lexicographically ascending by (ev3, ev4a, ev4b), ties of ev3 resolved by ev4a, then ev4b); slice = {{nu,val}} on the weight space, '
      'normalised so that the first nonzero component (ascending GT index nu) is exactly 1 (nonzero components only); computed by mpmath at 50 digits from the exact matrices, every V_red,k diagonal in the eigenbasis to 1e-40 |')
    A('#   HWRESID mu k residual : max over i < j of ||(E_ij (x) 1 + 1 (x) E_ij) Psi|| / ||Psi|| for the full vector of copy k reconstructed by the FLOAT level recursion of slice_engine_gl4 (double precision) from the copy\'s slice | GNORM2 mu k value : G-norm^2 sum g1 g2 Psi^2 of the full vector of copy k '
      '(= u^T PRVGRAM u, exact Gram, 30 digits) |')
    if has_vec:
        A('#   VEC mu k {{idx,val},...} : (DIM*DIM2 <= 2000) the full HW vector of copy k, idx = (nu-1)*DIM2 + q (1-based), 30 digits, same scale as the slice of COPY (= sum_i u_i of the exact full vectors of the PRV basis) |')
    else:
        A('#   (no VEC records in this file: DIM*DIM2 > 2000; regenerate with slice_engine_gl4.hw_from_slice[_exact] from the COPY slices (linearity: full vector = sum_i u_i * full vector of PRVBASIS_i))')
    return L


def write_copies_ref(lam, lamp=None, vec=None, dps=50, date=DATE, verbose=True, path=None):
    """write ref/gl4-cg-copies-<lam>[-<lamp>].txt (see the header of the file for the record formats)"""
    lam = tuple(lam)
    lamp = lam if lamp is None else tuple(lamp)
    t0 = time.perf_counter()
    R1, R2 = gr.GLNRep.of(lam), gr.GLNRep.of(lamp)
    n1, n2 = R1.dim, R2.dim
    if vec is None:
        vec = (n1 * n2 <= 2000)
    rs = gr.racah_speiser(lam, lamp)
    mus = sorted(rs, reverse=True)
    body = []
    st = Counter()
    worst = dict(res=0.0, offmax=mp.mpf(0))
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        for mu in mus:
            prv = se.prv_space(lam, lamp, mu)
            m = prv.m
            assert m == rs[mu]
            _, mats = triple_mats(lam, lamp, mu, prv)
            cps = copies(lam, lamp, mu, dps=dps, prv=prv, mats=mats)
            G, psis = prv_gram_exact(lam, lamp, mu, prv)
            st['exact_vectors'] += m
            body.append('SECTOR\t%s\t%d\t%d\t%d' % (fmu(mu), m, prv.w, gr.nW_of(lam, lamp, mu)))
            body.append('WSPACE\t%s\t{%s}' % (fmu(mu), ','.join(str(int(c) + 1) for c in prv.W)))
            for k in range(m):
                body.append('PRVBASIS\t%s\t{%s}' % (fmu(mu), ','.join('{%d,%s}' % (c + 1, fq(x)) for c, x in sorted(prv.vec(k).items()))))
            body.append('PRVFREE\t%s\t{%s}' % (fmu(mu), ','.join(str(int(prv.W[f]) + 1) for f in prv.free)))
            body.append('PRVGRAM\t%s\t%s' % (fmu(mu), f_mat(G)))
            for tag, n in (('VRED3', 'V3'), ('VRED4A', 'V4a'), ('VRED4B', 'V4b')):
                body.append('%s\t%s\t%s' % (tag, fmu(mu), f_mat(mats[n])))
            body.append('CHARPOLY3\t%s\t%s' % (fmu(mu), f_poly(charpoly_asc(mats['V3']))))
            comb = mat_lin([(c, mats[n]) for c, n in zip(GENERIC, TRIPLE)], m)
            body.append('CHARPOLYJ\t%s\t%s' % (fmu(mu), f_poly(charpoly_asc(comb))))
            hw = se.hw_from_slice(lam, lamp, mu, copy_float_slices(prv, cps))
            res, wd = se.hw_residual(hw, extra=True, per_vector=True)
            assert wd == 0
            worst['res'] = max(worst['res'], max(res))
            worst['offmax'] = max(worst['offmax'], cps.info['offmax'])
            for k, cp in enumerate(cps):
                body.append('COPY\t%s\t%d\t%s\t%s\t%s\t{%s}' % (fmu(mu), k + 1, fd(cp['ev'][0]), fd(cp['ev'][1]), fd(cp['ev'][2]),
                                                                ','.join('{%d,%s}' % (c + 1, fd(x)) for c, x in sorted(cp['slice'].items()))))
                body.append('HWRESID\t%s\t%d\t%.2e' % (fmu(mu), k + 1, res[k]))
                u = cp['coords']
                gn = mp.fsum(u[i] * mpq(G[i][j]) * u[j] for i in range(m) for j in range(m))
                body.append('GNORM2\t%s\t%d\t%s' % (fmu(mu), k + 1, fd(gn)))
                if vec:
                    keys = set()
                    for ps in psis:
                        keys.update(kk for kk, v in ps.items() if v != 0)
                    ent = {}
                    for kk in keys:
                        ent[kk] = mp.fsum(u[i] * mpq(psis[i].get(kk, Fr(0))) for i in range(m))
                    mx = max(abs(x) for x in ent.values())
                    items = sorted(((kk[0] * n2 + kk[1] + 1, x) for kk, x in ent.items() if abs(x) > mp.mpf(10) ** (-(dps - 15)) * mx))
                    body.append('VEC\t%s\t%d\t{%s}' % (fmu(mu), k + 1, ','.join('{%d,%s}' % (i, fd(x)) for i, x in items)))
                    st['vec'] += 1
                st['copies'] += 1
    finally:
        mp.mp.dps = old
    name = copies_ref_name(lam, lamp)
    head = _copies_hdr(lam, lamp, vec, name, date)
    head.append('# SUMMARY (self-check of this file): %d sectors, %d copies, exact full vectors used %d, VEC records %d; max HWRESID (float recursion, all E_ij) %.2e; worst off-diagonal entry of the three V_red,k in the joint eigenbasis %s' %
                (len(mus), st['copies'], st['exact_vectors'], st['vec'], worst['res'], mp.nstr(worst['offmax'], 3)))
    head.append('REP\t{%s}' % ','.join(str(x) for x in lam))
    head.append('REP2\t{%s}' % ','.join(str(x) for x in lamp))
    head.append('DIM\t%d' % n1)
    head.append('DIM2\t%d' % n2)
    gtb = lambda R: '{%s}' % ','.join('{%s}' % ','.join('{%s}' % ','.join(str(x) for x in row) for row in pat) for pat in R.pats)
    head.append('GTBASIS\t%s' % gtb(R1))
    head.append('GTBASIS2\t%s' % gtb(R2))
    head.append('WEIGHTS1\t{%s}' % ','.join('{%s}' % ','.join(str(x) for x in w) for w in R1.wt))
    head.append('WEIGHTS2\t{%s}' % ','.join('{%s}' % ','.join(str(x) for x in w) for w in R2.wt))
    head.append('GMET1\t{%s}' % ','.join(fq(x) for x in R1.g))
    head.append('GMET2\t{%s}' % ','.join(fq(x) for x in R2.g))
    head.append('NSECTORS\t%d' % len(mus))
    head.append('NCOPIES\t%d' % sum(rs.values()))
    head.append('GENERIC\t{%s}' % ','.join(fq(c) for c in GENERIC))
    path = path or os.path.join(REFDIR, name)
    with open(path, 'w', encoding='utf8', newline='\n') as f:
        f.write('\n'.join(head + body) + '\n')
    st['t'] = time.perf_counter() - t0
    st['bytes'] = os.path.getsize(path)
    st['path'] = path
    st['nsec'] = len(mus)
    st['maxres'] = worst['res']
    st['offmax'] = worst['offmax']
    if verbose:
        print('  wrote %s: %d sectors, %d copies, VEC %d, max HWRESID %.2e, worst joint-eigenbasis off-diagonal %s, %d bytes, %.1fs' % (name, len(mus), st['copies'], st['vec'], worst['res'], mp.nstr(worst['offmax'], 3), st['bytes'], st['t']), flush=True)
    return st


COPIES_REF_PRODUCTS = [((2, 1, 0, 0), (2, 1, 0, 0)), ((2, 1, 1, 0), (2, 1, 1, 0)), ((3, 2, 1, 0), (3, 2, 1, 0)), ((2, 2, 1, 0), (2, 1, 1, 0)), ((4, 2, 0, 0), (4, 2, 0, 0)), ((6, 4, 2, 0), (6, 4, 2, 0))]


def stage_ref():
    hdr('G3.6  REFERENCE FILES  ref/gl4-cg-copies-*.txt')
    out = {}
    for lam, lamp in COPIES_REF_PRODUCTS:
        out[(lam, lamp)] = write_copies_ref(lam, lamp)
    return out


# ======================================================================================== 11. independent re-verification of the reference files
_TOKX = re.compile(r'\{|\}|,|-?\d+/\d+|-?\d+\.\d+(?:[eE][-+]?\d+)?|-?\d+[eE][-+]?\d+|-?\d+')


def parse_nested_exact(s):
    """'{{1,2/3},{4,1.5e-3}}' -> nested lists; integer tokens -> int, rationals / decimals -> Fraction (EXACT, also for the 30-digit decimals); no eval"""
    toks = _TOKX.findall(s)
    pos = [0]

    def val():
        t = toks[pos[0]]
        pos[0] += 1
        if t == '{':
            out = []
            while toks[pos[0]] != '}':
                if toks[pos[0]] == ',':
                    pos[0] += 1
                    continue
                out.append(val())
            pos[0] += 1
            return out
        if re.fullmatch(r'-?\d+', t):
            return int(t)
        return Fr(t)
    return val()


def parse_copies_ref(path):
    """parse a gl4-cg-copies-*.txt: dict(lam, lamp, dim, dim2, gtbasis, gtbasis2, weights, weights2, g, g2, nsectors, ncopies, generic,
    sectors = {mu: dict(m, w, nW, wspace, prv=[{nu: Fraction}], free=[nu], gram, vred3, vred4a, vred4b, cp3, cpJ, copies=[dict(k, ev=(3 Fractions), slice={nu: Fraction}, hwres, gnorm2, vec)])})"""
    out = dict(sectors={})
    cur = None
    for line in open(path, encoding='utf8'):
        line = line.rstrip('\n')
        if not line or line.startswith('#'):
            continue
        f = line.split('\t')
        tag = f[0]
        if tag == 'REP':
            out['lam'] = tuple(parse_nested_exact(f[1]))
        elif tag == 'REP2':
            out['lamp'] = tuple(parse_nested_exact(f[1]))
        elif tag == 'DIM':
            out['dim'] = int(f[1])
        elif tag == 'DIM2':
            out['dim2'] = int(f[1])
        elif tag == 'GTBASIS':
            out['gtbasis'] = parse_nested_exact(f[1])
        elif tag == 'GTBASIS2':
            out['gtbasis2'] = parse_nested_exact(f[1])
        elif tag == 'WEIGHTS1':
            out['weights'] = [tuple(w) for w in parse_nested_exact(f[1])]
        elif tag == 'WEIGHTS2':
            out['weights2'] = [tuple(w) for w in parse_nested_exact(f[1])]
        elif tag == 'GMET1':
            out['g'] = [Fr(x) for x in parse_nested_exact(f[1])]
        elif tag == 'GMET2':
            out['g2'] = [Fr(x) for x in parse_nested_exact(f[1])]
        elif tag == 'NSECTORS':
            out['nsectors'] = int(f[1])
        elif tag == 'NCOPIES':
            out['ncopies'] = int(f[1])
        elif tag == 'GENERIC':
            out['generic'] = tuple(Fr(x) for x in parse_nested_exact(f[1]))
        elif tag == 'SECTOR':
            mu = tuple(parse_nested_exact(f[1]))
            cur = out['sectors'][mu] = dict(m=int(f[2]), w=int(f[3]), nW=int(f[4]), prv=[], copies=[])
        elif tag == 'WSPACE':
            cur['wspace'] = parse_nested_exact(f[2])
        elif tag == 'PRVBASIS':
            cur['prv'].append({nu: Fr(x) for nu, x in parse_nested_exact(f[2])})
        elif tag == 'PRVFREE':
            cur['free'] = parse_nested_exact(f[2])
        elif tag == 'PRVGRAM':
            cur['gram'] = [[Fr(x) for x in row] for row in parse_nested_exact(f[2])]
        elif tag in ('VRED3', 'VRED4A', 'VRED4B'):
            cur[tag.lower()] = [[Fr(x) for x in row] for row in parse_nested_exact(f[2])]
        elif tag == 'CHARPOLY3':
            cur['cp3'] = [Fr(x) for x in parse_nested_exact(f[2])]
        elif tag == 'CHARPOLYJ':
            cur['cpJ'] = [Fr(x) for x in parse_nested_exact(f[2])]
        elif tag == 'COPY':
            cur['copies'].append(dict(k=int(f[2]), ev=tuple(Fr(x) for x in f[3:6]), slice={nu: Fr(x) for nu, x in parse_nested_exact(f[6])}, hwres=None, gnorm2=None, vec=None))
        elif tag == 'HWRESID':
            cur['copies'][int(f[2]) - 1]['hwres'] = float(f[3])
        elif tag == 'GNORM2':
            cur['copies'][int(f[2]) - 1]['gnorm2'] = Fr(f[3])
        elif tag == 'VEC':
            cur['copies'][int(f[2]) - 1]['vec'] = {i: Fr(x) for i, x in parse_nested_exact(f[3])}
    return out


class SparseOp:
    """G2's Op2 as sparse CSR Kronecker groups in the ORTHONORMAL basis (G2's scaffold word_csr): apply(X) = sum_g P_g X Q_g^T for a sparse X (n1 x n2)"""

    def __init__(self, op):
        self.n1, self.n2 = op.n1, op.n2
        self.pairs = []
        for kw, comb in op.groups.items():
            if op.group_site == 2:
                P = sum(float(c) * op.rep1.word_csr(w) for w, c in comb.items())
                Q = op.rep2.word_csr(kw)
            else:
                P = op.rep1.word_csr(kw)
                Q = sum(float(c) * op.rep2.word_csr(w) for w, c in comb.items())
            self.pairs.append((sps.csr_matrix(P), sps.csr_matrix(Q).T.tocsr()))

    def apply(self, X):
        out = None
        for P, QT in self.pairs:
            Y = P @ (X @ QT)
            out = Y if out is None else out + Y
        return out


def rank_frac(rows):
    """exact rank of a list of rows (lists of Fractions)"""
    M = [list(r) for r in rows if any(x != 0 for x in r)]
    rk = 0
    ncols = len(M[0]) if M else 0
    for c in range(ncols):
        p = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
        if p is None:
            continue
        M[rk], M[p] = M[p], M[rk]
        inv = 1 / M[rk][c]
        for i in range(rk + 1, len(M)):
            if M[i][c] != 0:
                f = M[i][c] * inv
                M[i] = [x - f * y for x, y in zip(M[i], M[rk])]
        rk += 1
    return rk


def verify_copies_ref(path, full=True, verbose=True, max_full_sectors=None):
    """independent re-verification of a gl4-cg-copies-*.txt (see _verify_copies_ref for the list of checks); restores the mpmath precision afterwards"""
    old = mp.mp.dps
    try:
        return _verify_copies_ref(path, full, verbose, max_full_sectors)
    finally:
        mp.mp.dps = old


def _verify_copies_ref(path, full=True, verbose=True, max_full_sectors=None):
    """independent re-verification of a gl4-cg-copies-*.txt, parsed back from the TEXT (the reduction code is not used):
      (1) GTBASIS/WEIGHTS/GMET/DIM = glN_reps; the set of sectors and m = Racah-Speiser, w = |WSPACE| = weight space computed from WEIGHTS1, nW;
      (2) PRVBASIS: exact string conditions E_i^(lam'_i - lam'_{i+1} + 1) v = 0, rank m, canonical RREF structure at PRVFREE, span = the exact nullspace of the string conditions;
      (3) VRED3/4A/4B vs the FULL two-site operators (G2's Op2, term lists of ref/gl4-gaudin-ops.txt) applied exactly to the exact HW vectors of the PRV basis (G1's exact recursion): slice read-off
          sum_i A_ij b_i, exact; PRVGRAM vs the exact Gram of those vectors; CHARPOLY3 / CHARPOLYJ recomputed exactly; the three matrices commute exactly and are G-self-adjoint (exact);
      (4) COPY: slice in span(PRVBASIS), eigen-equation A_k u = ev_k u (u = slice at PRVFREE) to 1e-24, first nonzero component exactly 1, ordering, G-orthogonality of the copies and GNORM2 from PRVGRAM;
      (5) full = True: the FLOAT full vector of every copy (level recursion from the COPY slice): HW residual of all E_ij <= 1e-12, GNORM2 vs the float norm, VEC (if present) vs the recursion, and the
          FULL two-site operators V3, V4a, V4b (G2's scaffold, sparse) applied to it: ||O Psi - ev Psi|| / (||Psi|| (1 + |ev|)) <= 1e-12.
    Returns (total number of failed checks, Counter of failures per category, dict of the largest deviations)."""
    D = parse_copies_ref(path)
    lam, lamp = D['lam'], D['lamp']
    N = len(lam)
    R1, R2 = gr.GLNRep.of(lam), gr.GLNRep.of(lamp)
    bad = Counter()
    dev = defaultdict(float)
    t0 = time.perf_counter()
    # (1) header
    bad['hdr'] += (D['dim'] != R1.dim) + (D['dim2'] != R2.dim)
    bad['hdr'] += ([tuple(tuple(r) for r in p) for p in D['gtbasis']] != [tuple(tuple(r) for r in p) for p in R1.pats])
    bad['hdr'] += ([tuple(tuple(r) for r in p) for p in D['gtbasis2']] != [tuple(tuple(r) for r in p) for p in R2.pats])
    bad['hdr'] += (D['weights'] != R1.wt) + (D['weights2'] != R2.wt) + (D['g'] != list(R1.g)) + (D['g2'] != list(R2.g))
    bad['hdr'] += (D['generic'] != GENERIC)
    rs = gr.racah_speiser(lam, lamp)
    bad['sector'] += (set(D['sectors']) != set(rs)) + (D['nsectors'] != len(rs)) + (D['ncopies'] != sum(rs.values()))
    ops2 = {n: two_site_op(n, lam, lamp) for n in TRIPLE}
    sops = {n: SparseOp(ops2[n]) for n in TRIPLE} if full else None
    top1 = R1.top
    top2 = R2.top
    nfull = 0
    mp.mp.dps = 60
    for mu in sorted(D['sectors'], reverse=True):
        S = D['sectors'][mu]
        m = S['m']
        # (1) sector fields
        w0 = tuple(mu[a] - lamp[a] for a in range(N))
        Wt = [int(c) for c in R1.wsp.get(w0, [])]
        bad['sector'] += (m != rs.get(mu, 0)) + (S['wspace'] != [c + 1 for c in Wt]) + (S['w'] != len(Wt)) + (S['nW'] != gr.nW_of(lam, lamp, mu)) + (len(S['prv']) != m) + (len(S['copies']) != m)
        W0 = [c - 1 for c in S['wspace']]
        # (2) PRVBASIS
        vecs = [{nu - 1: x for nu, x in b.items()} for b in S['prv']]
        for v in vecs:
            for i in range(1, N):
                kk = lamp[i - 1] - lamp[i] + 1
                u = dict(v)
                for _ in range(kk):
                    u = R1.apply((i, i + 1), u)
                    if not u:
                        break
                bad['prv_string'] += bool(u)
        Bm = [[b.get(c + 1, Fr(0)) for c in W0] for b in S['prv']]
        bad['prv_rank'] += (rank_frac(Bm) != m)
        fr = S['free']
        pos1 = {c + 1: a for a, c in enumerate(W0)}
        for k, b in enumerate(S['prv']):
            for l, f_ in enumerate(fr):
                bad['prv_free'] += (b.get(f_, Fr(0)) != (1 if k == l else 0))
        # independent nullspace dimension
        rows = []
        for i in range(1, N):
            kk = lamp[i - 1] - lamp[i] + 1
            if tuple(w0[a] + (kk if a == i - 1 else 0) - (kk if a == i else 0) for a in range(N)) not in R1.wsp:
                continue
            imgs = se.power_images(R1, (i, i + 1), kk, W0)
            byrow = defaultdict(lambda: [Fr(0)] * len(W0))
            for a, vec in enumerate(imgs):
                for r, x in vec.items():
                    byrow[r][a] = x
            rows.extend(byrow.values())
        nb, _ = se._rref_nullspace(rows, len(W0))
        bad['prv_null'] += (len(nb) != m)
        # exact HW vectors of the PRV basis (G1's exact recursion; the slice is the text's PRVBASIS entry)
        psis = []
        for b in S['prv']:
            psi, nchk, viol = se.hw_from_slice_exact(lam, lamp, mu, [b.get(c + 1, Fr(0)) for c in W0])
            bad['exact_twin_viol'] += sum(viol)
            psis.append(psi)
        # (3) VRED vs the full two-site operators
        mats = {}
        for n, tg in (('V3', 'vred3'), ('V4a', 'vred4a'), ('V4b', 'vred4b')):
            A = S[tg]
            mats[n] = A
            bad['vred_shape'] += (len(A) != m) or any(len(r) != m for r in A)
            sop = slice_op(n, lam, lamp)
            for j in range(m):
                sl = sop.apply(by_q(psis[j]))
                img = {}
                for i in range(m):
                    if A[i][j] != 0:
                        for nu, x in S['prv'][i].items():
                            img[nu - 1] = img.get(nu - 1, 0) + A[i][j] * x
                img = {c: x for c, x in img.items() if x != 0}
                bad['vred_vs_full_' + n] += (sl != img)
        # PRVGRAM
        g1, g2 = R1.g, R2.g
        G = [[Fr(0)] * m for _ in range(m)]
        for i in range(m):
            for j in range(i, m):
                a_, b_ = psis[i], psis[j]
                if len(a_) > len(b_):
                    a_, b_ = b_, a_
                s_ = Fr(0)
                for (nu, q), v in a_.items():
                    if v != 0:
                        wv = b_.get((nu, q))
                        if wv is not None and wv != 0:
                            s_ += g1[nu] * g2[q] * v * wv
                G[i][j] = G[j][i] = s_
        bad['gram'] += (G != S['gram'])
        bad['cp3'] += (charpoly_asc(mats['V3']) != S['cp3'])
        comb = mat_lin([(c, mats[n]) for c, n in zip(GENERIC, TRIPLE)], m)
        bad['cpJ'] += (charpoly_asc(comb) != S['cpJ'])
        bad['comm'] += sum(0 if mat_is_zero(mat_comm(mats[a], mats[b])) else 1 for a, b in itertools.combinations(TRIPLE, 2))
        bad['selfadj'] += selfadjoint_bad(mats, S['gram'])
        # (4) COPY
        Bmp = [[mpq(x) for x in row] for row in Bm]
        Gm = [[mpq(x) for x in row] for row in S['gram']]
        Am = {n: [[mpq(x) for x in row] for row in mats[n]] for n in TRIPLE}
        us = []
        prev = None
        for k, cp in enumerate(S['copies']):
            bad['copy_k'] += (cp['k'] != k + 1)
            sl = cp['slice']
            nzs = sorted(nu for nu, x in sl.items() if x != 0)
            bad['copy_norm'] += (not nzs) or (sl[nzs[0]] != 1) or any(nu not in pos1 for nu in nzs)
            u = [mpq(sl.get(f_, Fr(0))) for f_ in fr]
            us.append(u)
            # slice in the span of the PRV basis
            rec = [mp.fsum(u[i] * Bmp[i][a] for i in range(m)) for a in range(len(W0))]
            dd = max([abs(rec[a] - mpq(sl.get(c + 1, Fr(0)))) for a, c in enumerate(W0)] + [mp.mpf(0)])
            dev['copy_span'] = max(dev['copy_span'], float(dd))
            bad['copy_span'] += (dd > mp.mpf(10) ** -27)
            for t, n in enumerate(TRIPLE):
                ev = mpq(cp['ev'][t])
                Au = [mp.fsum(Am[n][i][j] * u[j] for j in range(m)) for i in range(m)]
                r_ = max(abs(Au[i] - ev * u[i]) for i in range(m)) / (max(abs(x) for x in u) * (1 + abs(ev)))
                dev['copy_eig'] = max(dev['copy_eig'], float(r_))
                bad['copy_eig'] += (r_ > mp.mpf(10) ** -24)
            gn = mp.fsum(u[i] * Gm[i][j] * u[j] for i in range(m) for j in range(m))
            bad['gnorm2'] += (abs(gn - mpq(cp['gnorm2'])) > mp.mpf(10) ** -24 * abs(gn))
            if prev is not None:
                key = tol_cmp(mp.mpf(10) ** -20)
                bad['copy_order'] += (key(tuple(mpq(x) for x in prev)) > key(tuple(mpq(x) for x in cp['ev'])))
            prev = cp['ev']
        for k in range(m):
            for l in range(k + 1, m):
                o = mp.fsum(us[k][i] * Gm[i][j] * us[l][j] for i in range(m) for j in range(m))
                gk = mp.fsum(us[k][i] * Gm[i][j] * us[k][j] for i in range(m) for j in range(m))
                gl = mp.fsum(us[l][i] * Gm[i][j] * us[l][j] for i in range(m) for j in range(m))
                rel = abs(o) / mp.sqrt(gk * gl)
                dev['copy_orth'] = max(dev['copy_orth'], float(rel))
                bad['copy_orth'] += (rel > mp.mpf(10) ** -24)
        # (5) float full vectors and the full two-site operators
        if full and (max_full_sectors is None or nfull < max_full_sectors):
            nfull += 1
            F = np.zeros((len(W0), m))
            for j, cp in enumerate(S['copies']):
                for nu, x in cp['slice'].items():
                    F[pos1[nu], j] = float(x)
            hw = se.hw_from_slice(lam, lamp, mu, F)
            res, wd = se.hw_residual(hw, extra=True, per_vector=True)
            bad['hw_resid'] += sum(1 for r in res if r > 1e-12) + (wd != 0)
            dev['hw_resid'] = max([dev['hw_resid']] + list(res))
            n1, n2 = hw.dim1, hw.dim2
            for j, cp in enumerate(S['copies']):
                X = sps.csr_matrix((hw.vals_hat[:, j], (hw.rows, hw.cols)), shape=(n1, n2))
                nx = math.sqrt(float(X.multiply(X).sum()))
                gf = float(nx ** 2)
                bad['gnorm2_float'] += (abs(gf - float(cp['gnorm2'])) > 1e-11 * gf)
                dev['gnorm2_float'] = max(dev['gnorm2_float'], abs(gf - float(cp['gnorm2'])) / gf)
                for t, n in enumerate(TRIPLE):
                    Y = sops[n].apply(X) - float(cp['ev'][t]) * X
                    r_ = math.sqrt(float(Y.multiply(Y).sum())) / (nx * (1 + abs(float(cp['ev'][t]))))
                    dev['full_eig'] = max(dev['full_eig'], r_)
                    bad['full_eig'] += (r_ > 1e-12)
                if cp['vec'] is not None:
                    vals = hw.vals[:, j]
                    flat = hw.flat() + 1
                    d = {int(i): float(v) for i, v in zip(flat, vals)}
                    mx = max(abs(float(x)) for x in cp['vec'].values())
                    dv = max([abs(d.get(i, 0.0) - float(x)) for i, x in cp['vec'].items()] + [abs(v) for i, v in d.items() if i not in cp['vec'] and abs(v) > 1e-14 * mx]) / mx
                    dev['vec'] = max(dev['vec'], dv)
                    bad['vec'] += (dv > 1e-12)
    total = sum(bad.values())
    if verbose:
        print('  %s: %d sectors, %d copies | TOTAL MISMATCHES %d %s | max deviations: copy-in-span %.1e, eigen-equation %.1e, G-orthogonality %.1e%s  (%.1fs)' %
              (os.path.basename(path), len(D['sectors']), sum(len(s['copies']) for s in D['sectors'].values()), total, dict(bad) if total else '',
               dev['copy_span'], dev['copy_eig'], dev['copy_orth'],
               (', HW residual %.1e, full-operator eigen-equation %.1e, GNORM2 %.1e, VEC %.1e (full vectors checked in %d sectors)' % (dev['hw_resid'], dev['full_eig'], dev['gnorm2_float'], dev['vec'], nfull)) if full else '',
               time.perf_counter() - t0), flush=True)
    return total, bad, dict(dev)


def mutation_test_copies(name='gl4-cg-copies-2100.txt', verbose=True):
    """a verifier that cannot fail proves nothing: corrupt ONE record at a time in a copy of the file (in the m = 2 sector mu = (3,2,1,0), or in the header); every corruption must be
    detected (a nonzero number of failed checks, or an exception of the verifier), the pristine copy must pass"""
    import tempfile
    src = os.path.join(REFDIR, name)
    lines = open(src, encoding='utf8').read().split('\n')
    tmp = os.path.join(tempfile.gettempdir(), 'g3_mut_%d.txt' % os.getpid())
    res = []
    sec = '{3,2,1,0}'

    def find(tag, sector=True):
        for i, ln in enumerate(lines):
            f = ln.split('\t')
            if f[0] == tag and (not sector or f[1] == sec):
                return i
        raise KeyError(tag)

    def run(label, tag, field, fn, sector=True):
        L = list(lines)
        if tag is not None:
            i = find(tag, sector)
            f = L[i].split('\t')
            f[field] = fn(f[field])
            L[i] = '\t'.join(f)
        with open(tmp, 'w', encoding='utf8', newline='\n') as fh:
            fh.write('\n'.join(L))
        try:
            tot, bad, dev = verify_copies_ref(tmp, full=True, verbose=False)
            note = ''
        except Exception as e:                                # a corrupted record that breaks the parser / the checks also counts as detected
            tot, bad, note = 1, {}, ' (verifier raised %s)' % type(e).__name__
        res.append((label, tot, dict(bad), note))

    def bump_mat(s, d=1):
        A = parse_nested_exact(s)
        A[0][0] = Fr(A[0][0]) + d
        return f_mat(A)

    def bump_poly(s):
        c = parse_nested_exact(s)
        c[1] = Fr(c[1]) + 1
        return f_poly(c)

    def bump_pairs(s, which=0, d=1):
        P = parse_nested_exact(s)
        P[which][1] = Fr(P[which][1]) + d
        return '{%s}' % ','.join('{%d,%s}' % (nu, fq(Fr(x))) for nu, x in P)

    eps = Fr(1, 10 ** 6)
    run('pristine copy', None, 0, None)
    run('VRED3 entry changed', 'VRED3', 2, bump_mat)
    run('VRED4A entry changed', 'VRED4A', 2, bump_mat)
    run('VRED4B entry changed', 'VRED4B', 2, bump_mat)
    run('CHARPOLY3 coefficient changed', 'CHARPOLY3', 2, bump_poly)
    run('CHARPOLYJ coefficient changed', 'CHARPOLYJ', 2, bump_poly)
    run('PRVGRAM entry changed', 'PRVGRAM', 2, bump_mat)
    run('PRVBASIS value changed', 'PRVBASIS', 2, bump_pairs)
    run('PRVFREE index changed', 'PRVFREE', 2, lambda s: '{%d%s' % (int(s[1:s.index(',')]) + 1, s[s.index(','):]))
    run('COPY eigenvalue ev3 changed by 1e-6', 'COPY', 3, lambda s: str(Fr(s) + eps))
    run('COPY eigenvalue ev4b changed by 1e-6', 'COPY', 5, lambda s: str(Fr(s) + eps))
    run('COPY slice value changed by 1e-6', 'COPY', 6, lambda s: bump_pairs(s, which=1, d=eps))
    run('GNORM2 changed by 1e-6 (relative)', 'GNORM2', 3, lambda s: str(Fr(s) * (1 + eps)))
    run('VEC value changed by 1e-6', 'VEC', 3, lambda s: bump_pairs(s, which=0, d=eps))
    run('SECTOR multiplicity changed', 'SECTOR', 2, lambda s: str(int(s) + 1))
    run('WSPACE index changed', 'WSPACE', 2, lambda s: '{%d%s' % (int(s[1:s.index(',')]) + 1, s[s.index(','):]))
    run('GMET1 entry changed', 'GMET1', 1, lambda s: '{%s%s' % (str(Fr(s[1:s.index(',')]) + 1), s[s.index(','):]), sector=False)
    run('GENERIC changed', 'GENERIC', 1, lambda s: '{1,7/19,12/13}', sector=False)
    try:
        os.remove(tmp)
    except OSError:
        pass
    if verbose:
        for label, tot, bad, note in res:
            ok = (tot == 0) if label == 'pristine copy' else (tot > 0)
            nz = {k: v for k, v in bad.items() if v}
            print('    %-44s -> %3d failed checks %s %s%s' % (label, tot, (str(nz) if nz else ''), ('OK (passes)' if label == 'pristine copy' and ok else ('DETECTED' if ok else 'NOT DETECTED / UNEXPECTED')), note))
    return res


def compare_with_g1(lam, lamp, verbose=True):
    """the G3 file against G1's ref/gl4-cg-slices-*.txt of the same product (if it exists): SECTOR fields, WSPACE and PRVBASIS identical (exact), and on every m = 1 sector the COPY slice (30 digits
    vs G1's exact rationals, 1e-28), the GNORM2 (1e-27 relative) and, where both files have it, the full vector VEC (1e-28)"""
    p3 = os.path.join(REFDIR, copies_ref_name(lam, lamp))
    p1 = os.path.join(REFDIR, se.ref_name(lam, lamp))
    if not os.path.exists(p1):
        if verbose:
            print('   %-22s no G1 file (%s) to compare with' % (prod_str(lam, lamp), os.path.basename(p1)))
        return None
    D3 = parse_copies_ref(p3)
    D1 = se.parse_ref(p1)
    nbad = n1 = nvec = 0
    tolc = Fr(1, 10 ** 28)
    for mu, S3 in D3['sectors'].items():
        S1 = D1['sectors'][mu]
        nbad += (S3['m'] != S1['m']) + (S3['w'] != S1['w']) + (S3['nW'] != S1['nW']) + (S3['wspace'] != S1['wspace']) + (S3['prv'] != S1['prv'])
        if S3['m'] == 1:
            n1 += 1
            c1, c3 = S1['copy'], S3['copies'][0]
            nbad += not (set(c1) == set(c3['slice']) and all(abs(c3['slice'][k] - c1[k]) < tolc for k in c1))
            g1_ = Fr(S1['gnorm2']) if not isinstance(S1['gnorm2'], float) else Fr(S1['gnorm2'])
            nbad += abs(c3['gnorm2'] - g1_) > Fr(1, 10 ** 27) * g1_
            if S1.get('vec') is not None and c3['vec'] is not None:
                nvec += 1
                v1, v3 = S1['vec'], c3['vec']
                nbad += not (set(v1) == set(v3) and all(abs(v3[k] - v1[k]) < tolc for k in v1))
    if verbose:
        print('   %-22s vs %s: %d sectors (SECTOR fields, WSPACE, PRVBASIS exact), %d m = 1 copies (slice, GNORM2), %d full vectors: %d mismatches' % (prod_str(lam, lamp), os.path.basename(p1), len(D3['sectors']), n1, nvec, nbad))
    return len(D3['sectors']), n1, nvec, nbad


def stage_verify(full_big=True):
    hdr('G3.6  INDEPENDENT RE-VERIFICATION of the reference files (parsed from the text; G1/G2 builders: full two-site operators, exact recursion, float recursion)')
    out = {}
    for lam, lamp in COPIES_REF_PRODUCTS:
        path = os.path.join(REFDIR, copies_ref_name(lam, lamp))
        out[(lam, lamp)] = verify_copies_ref(path, full=True)
    print('\n  consistency with the G1 slice files ref/gl4-cg-slices-*.txt (same products where G1 wrote a file):')
    for lam, lamp in COPIES_REF_PRODUCTS:
        out[('g1', lam, lamp)] = compare_with_g1(lam, lamp)
    print('\n  mutation test of the verifier (ONE record corrupted at a time in a copy of gl4-cg-copies-2100.txt, in the m = 2 sector (3,2,1,0) or the header):')
    out['mutation'] = mutation_test_copies()
    return out


# ======================================================================================== CLI
STAGES = ['order', 'reduce', 'gl3', 'direct', 'copies', 'bench', 'ref', 'verify', 'all']


def main(argv):
    cmd = argv[1] if len(argv) > 1 else 'all'
    ext = 'ext' in argv[2:]
    if cmd not in STAGES:
        print('usage: python -X utf8 slice_ops_gl4.py <stage> [ext], stage in %s   (bench ext: also the extension rows a = 7, 8 / [12,8,4,0]^2, slices only)' % STAGES)
        return 1
    t0 = time.perf_counter()
    print('# slice_ops_gl4.py stage=%s%s  (Python %s, numpy %s, scipy %s, mpmath %s; %s; %s)' % (cmd, ' ext' if ext else '', platform.python_version(), np.__version__, __import__('scipy').__version__,
                                                                                             mp.__version__, platform.platform(), time.strftime('%Y-%m-%d %H:%M')), flush=True)
    if cmd in ('order', 'all'):
        stage_order()
    if cmd in ('reduce', 'all'):
        stage_reduce()
    if cmd in ('gl3', 'all'):
        stage_gl3()
    if cmd in ('direct', 'all'):
        stage_direct()
    if cmd in ('copies', 'all'):
        stage_copies()
    if cmd in ('bench', 'all'):
        stage_bench(extended=ext)
    if cmd in ('ref', 'all'):
        stage_ref()
    if cmd in ('verify', 'all'):
        stage_verify()
    print('\ntotal time %.1fs' % (time.perf_counter() - t0))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
