"""D (round 6): the Gauss dictionary for general Lambda -- c-hat, phi_u = e1^T rho_Lambda(n'), the PRV slice space M_mu and the (m+1)-point dictionary.
Brief: SPEC-gauss-dictionary.md (this file implements it); report: RESULTS-gauss-dictionary.md.  Python only; run from Paul/Python/gaudin_sl3 with `python -X utf8`.
CLI: `python -X utf8 gauss_dictionary.py [--stages polys,prv,timing,branches,cond,points,scale,misc,crossR] [--amax 14]` (everything printed is in gauss_dictionary_output.txt).

The identity (plan round 6, section 1; constant exactly 1 on all 26 known branches):

    L_B(w1, w1p) = c-hat(phi; w1, w1p) * <phi_u(w1, w1p) (x) e1|Psi_0> / X_B ,
    <phi_u (x) e1|Psi_0> = sum_{nu in W_mu} [e1^T rho_Lambda(n')]_nu C_nu ,   W_mu = V_Lambda[mu - Lambda] ,   C_nu = <e_nu (x) e1|Psi_0>,
    n' = exp(a' E12) exp(g' E13) exp(c' E23)   (Lambda-independent group data, `gauss_params`),
    c-hat = (phi1-phi2)^(2 lam2) prod_{w in {w1,w1p}} (phi2-w)^(lam1-lam2) (phi3-w)^lam1   (`chat`, lam3 = 0).

The slice (C_nu) lies in the PRV space M_mu = {v in W_mu : E12^(l1-l2+1) v = E23^(l2-l3+1) v = 0} (dim m = multiplicity), so L_B at m+1 generic points fixes it:
    L_B(w_k)/c-hat(w_k) = sum_j x_j <phi_u(w_k)|b_j>,   b_j = prv_basis(lam, mu)   (least squares; the residual is the consistency check).
Everything lives in the ONE-SITE space of dimension dim V_Lambda (512 at [14,7,0], 3375 at [28,14,0]); no two-site object is ever built.

Conventions: one-site GT basis = gl3_reps.GL3Rep order ("my" order, HW = the last index); the basis is non-orthonormal.  E12, E23 = gl3_reps raising matrices,
E13 := [E12, E23] (as gaudin_chain.Chain).  Matrices act on kets as columns (M[row][col]); covectors act from the right: (e1^T M)_c = M[hw][c].
Numbers: Fraction (exact) for int/Fraction/sympy-Rational inputs, mpmath mpf for mpf/float inputs (global mp.mp.dps; the entry points take `dps`).
In the PRV basis (RREF nullspace form) the coordinates of any v in M_mu are x_j = v[f_j], f_j = prv_free_indices(lam, mu)[j] (`coords_of`).

PUBLIC INTERFACE (what round 7 calls)
  slice_from_points(lam, mu, phi, points, values, dps=None) -> (x, residual, slice_vector)      [dps default: max(mp.mp.dps, 30); pass dps >= the digits of your data]
      lam = (l1, l2, 0); mu = (m1, m2, m3) dominant with M_mu != 0; phi = (phi1, phi2, phi3) (Fraction or mpf; the phi of the eps-regime that produced L_B);
      points = [(w1_k, w1p_k)], >= m+1 pairs (Fraction/int exact, or mpf), w1 != w1p, never equal to a phi_i; values = [L_B(w1_k, w1p_k)] (mpf/Fraction/float, one per point);
      x = list of m mpf: coordinates in prv_basis(lam, mu) of RAW = sum x_j b_j = C/X_B (the overall scale 1/X_B is kept; X_B * RAW = the G-normalised slice, verified to 1e-53...1e-121 on the 26 branches);
      residual = ||D(Ax - y)||/||D y|| (mpf; D = row scaling 1/||A_k||, y_k = values[k]/c-hat(w_k); exact squared residual on the exact path) -- the consistency check;
      slice_vector = list of dim V_Lambda mpf in my GT order, zero outside W_mu, scaled so that its first nonzero entry (|.| > 10^(-dps/2) max|.|) is 1.
      Exact path: all of phi, points, values rational -> exact rational least squares (x_exact in slice_from_points_ex).
  slices_from_points(lam, mu, phi, points, values_per_branch, dps=None) -> [(x, residual, slice_vector), ...]   (all copies of one irrep, one dictionary matrix)
  DictionaryProblem(lam, mu, phi, points, dps=None) with .solve(values) -> dict(x, residual, slice, raw, pivot, x_exact, exact) and .cond();
  slice_from_points_ex(...) -> dict (+ A, chat, sv, cond, timing).
  points: default_points(phi, n); shared_points(lam, phi) (RECOMMENDED: shared_points(lam, phi)[:m+1] for an irrep of multiplicity m; one nested list for the whole Lambda, worst cond
          43 at [14,7,0], 3e3 at [28,14,0]); choose_points(lam, mu, phi, n) (per-irrep greedy, worst cond 8 / 92).
  Gauss data: gauss_params(phi, w1, w1p) -> (a', g', c'); chat(lam, phi, w1, w1p); gauss_covector(lam, phi, w1, w1p, idx=None, method='terms'|'sweep');
          gauss_covector_poly(lam, phi) -> (list of Poly or None, c-hat Poly) (= untwist_limit.gf_polys for the adjoint); uncleared_components(lam, phi); pole_order(...).
  PRV: prv_basis(lam, mu), prv_free_indices(lam, mu), coords_of(lam, mu, vec), multiplicity(lam, mu) (Racah-Speiser), irreps(lam), dim_irrep(mu), weight_space(lam, nu).
  The identity itself: synthetic_values(lam, mu, phi, points, slice) = c-hat <phi_u|slice> (exact for exact input) -- for self-consistency tests and held-out checks of L_B.
  rep1(lam) -> Rep1: sparse exact one-site E12, E23, E13, GT patterns/weights, the Lambda-only tensor T (phi_nu = sum T a'^p g'^q c'^r).
"""
import sys, os, time, pickle, argparse, itertools
from fractions import Fraction as Fr
from math import factorial
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))
import sympy as sp
import mpmath as mp
from gl3_reps import GL3Rep, _raise, _weight
from gaudin_chain import SM, comm

W1, W2 = sp.symbols('W1 W2')           # same Symbol objects (name-equal) as untwist_limit.W1S / W2S


# ------------------------------------------------------------------------------------------------------------- numbers
def _num(x):
    """exact -> Fraction, numeric -> mpf"""
    if isinstance(x, Fr):
        return x
    if isinstance(x, int):
        return Fr(x)
    if isinstance(x, sp.Rational):
        return Fr(int(x.p), int(x.q))
    if isinstance(x, (mp.mpf, float)):
        return mp.mpf(x)
    if isinstance(x, str):
        return Fr(x)
    raise TypeError('unsupported number type %r' % (type(x),))


def to_mp(x):
    """Fraction / int / mpf / float -> mpf at the current precision"""
    if isinstance(x, Fr):
        return mp.mpf(x.numerator) / x.denominator
    if isinstance(x, int):
        return mp.mpf(x)
    return mp.mpf(x)


def _unify(*xs):
    """list of numbers all Fraction (if every input is exact) or all mpf"""
    ys = [_num(x) for x in xs]
    if all(isinstance(y, Fr) for y in ys):
        return ys
    return [to_mp(y) for y in ys]


# ------------------------------------------------------------------------------------------------------------- the Gauss data
def gauss_params(phi, w1, w1p):
    """(a', g', c') of n' = exp(a'E12) exp(g'E13) exp(c'E23)  (round-4 plan section 0; Lambda-independent). Exact for exact inputs, mpf otherwise."""
    p1, p2, p3, w, wp = _unify(phi[0], phi[1], phi[2], w1, w1p)
    d = w - wp
    a = (p1 - p2) * d / ((p2 - w) * (p2 - wp))
    c = (p1 - p3) * (p2 - p3) * d / ((p1 - p2) * (p3 - w) * (p3 - wp))
    g = (p3 - p1) * (p3 - p2) * d / ((p3 - w) * (p2 - wp) * (p3 - wp))
    return a, g, c


def chat(lam, phi, w1, w1p):
    """c-hat = (phi1-phi2)^(2 lam2) prod_{w in {w1,w1p}} (phi2-w)^(lam1-lam2) (phi3-w)^lam1  (lam3 = 0; RESULTS-nu-general section 0)"""
    l1, l2, l3 = lam
    if l3 != 0:
        raise ValueError('c-hat is only defined for lam3 = 0 (lam3 > 0 not attempted)')
    p1, p2, p3, w, wp = _unify(phi[0], phi[1], phi[2], w1, w1p)
    return (p1 - p2) ** (2 * l2) * ((p2 - w) ** (l1 - l2) * (p3 - w) ** l1) * ((p2 - wp) ** (l1 - l2) * (p3 - wp) ** l1)


# ------------------------------------------------------------------------------------------------------------- sparse one-site machinery
def _apply_cov(rows, cov):
    """covector (dict r -> val) times sparse matrix (rows[r] = {c: val})"""
    out = {}
    for r, x in cov.items():
        row = rows.get(r)
        if row:
            for c, v in row.items():
                out[c] = out.get(c, 0) + x * v
    return out


def _apply_ket(cols, ket):
    """sparse matrix (cols[c] = {r: val}) times ket (dict c -> val)"""
    out = {}
    for c, x in ket.items():
        col = cols.get(c)
        if col:
            for r, v in col.items():
                out[r] = out.get(r, 0) + v * x
    return {r: v for r, v in out.items() if v != 0}


def _cols_of(sm):
    cols = {}
    for r, row in sm.r.items():
        for c, v in row.items():
            cols.setdefault(c, {})[r] = v
    return cols


def _exp_apply_cov(rows, cov, x):
    """covector * exp(x E) for a nilpotent sparse E (series terminates structurally)"""
    tot = dict(cov)
    term = cov
    k = 0
    while term:
        k += 1
        term = {c: v * x / k for c, v in _apply_cov(rows, term).items()}
        for c, v in term.items():
            tot[c] = tot.get(c, 0) + v
    return tot


class Rep1:
    """One-site irrep Lambda in the gl3_reps GT basis: patterns, weights, sparse exact E12, E23, E13 (rows and columns), the HW index and the
    Lambda-only tensor T = {nu: [(p, q, r, T)]}, T = <hw|E12^p E13^q E23^r|e_nu>/(p! q! r!) (phi_nu(n') = sum T a'^p g'^q c'^r)."""

    def __init__(self, lam):
        lam = tuple(int(x) for x in lam)
        if len(lam) != 3 or not (lam[0] >= lam[1] >= lam[2] >= 0):
            raise ValueError('lam must be a decreasing triple of non-negative integers: %r' % (lam,))
        self.lam = lam
        t0 = time.time()
        g = GL3Rep(lam)                               # only patterns + index (no dense generators are built)
        self.patterns = g.patterns
        self.index = g._index
        self.dim = len(self.patterns)
        self.hw = self.index[tuple(tuple(lam[:k]) for k in (1, 2, 3))]
        self.weights = [tuple(_weight(k, pat) for k in (1, 2, 3)) for pat in self.patterns]
        self.wspace = {}
        for i, w in enumerate(self.weights):
            self.wspace.setdefault(w, []).append(i)
        self.E12 = self._raising(1)
        self.E23 = self._raising(2)
        self.E13 = comm(self.E12, self.E23)            # [E12, E23] = E13, as gaudin_chain.Chain
        self.C12, self.C23, self.C13 = _cols_of(self.E12), _cols_of(self.E23), _cols_of(self.E13)
        self.t_setup = time.time() - t0
        self._T = None
        self._Tmp = {}
        self.t_terms = None

    def _raising(self, k):
        rows = {}
        for col, pat in enumerate(self.patterns):
            for target, coeff in _raise(k, pat).items():
                r = self.index.get(target)               # boundary terms of the raw formula are discarded, as in GL3Rep
                if r is not None:
                    q = sp.Rational(coeff)
                    rows.setdefault(r, {})[col] = Fr(int(q.p), int(q.q))
        return SM(self.dim, rows)

    # -------- the Lambda-only tensor
    @property
    def T(self):
        if self._T is None:
            t0 = time.time()
            T = {}
            cov0 = {self.hw: Fr(1)}
            covp = cov0; p = 0
            while covp:
                covq = covp; q = 0
                while covq:
                    covr = covq; r = 0
                    while covr:
                        for nu, v in covr.items():
                            if v != 0:
                                T.setdefault(nu, []).append((p, q, r, v))
                        r += 1
                        covr = {c: v / r for c, v in _apply_cov(self.E23.r, covr).items()}
                    q += 1
                    covq = {c: v / q for c, v in _apply_cov(self.E13.r, covq).items()}
                p += 1
                covp = {c: v / p for c, v in _apply_cov(self.E12.r, covp).items()}
            self._T = T
            self.t_terms = time.time() - t0
        return self._T

    def _Tm(self, exact):
        """T converted to the working number type (cached per precision)"""
        if exact:
            return self.T
        key = mp.mp.dps
        if key not in self._Tmp:
            self._Tmp[key] = {nu: [(p, q, r, to_mp(v)) for (p, q, r, v) in lst] for nu, lst in self.T.items()}
        return self._Tmp[key]

    # -------- covectors
    def covector_terms(self, a, g, c, idx=None):
        """phi_nu = sum T a^p g^q c^r for nu in idx (default all): dict nu -> number"""
        exact = isinstance(a, Fr) and isinstance(g, Fr) and isinstance(c, Fr)
        Tm = self._Tm(exact)
        one = Fr(1) if exact else mp.mpf(1)
        zero = Fr(0) if exact else mp.mpf(0)
        nmax = max(self.lam[0] + 2, 3)
        ap, gq, cr = [one], [one], [one]
        for _ in range(nmax):
            ap.append(ap[-1] * a); gq.append(gq[-1] * g); cr.append(cr[-1] * c)
        out = {}
        for nu in (range(self.dim) if idx is None else idx):
            s = zero
            for (p, q, r, t) in Tm.get(nu, ()):
                s += t * ap[p] * gq[q] * cr[r]
            out[nu] = s
        return out

    def covector_sweep(self, a, g, c):
        """phi = e_hw^T exp(a E12) exp(g E13) exp(c E23) by direct application (dict nu -> number; only the nonzero support)"""
        exact = isinstance(a, Fr) and isinstance(g, Fr) and isinstance(c, Fr)
        cov = {self.hw: Fr(1) if exact else mp.mpf(1)}
        for rows, x in ((self.E12.r, a), (self.E13.r, g), (self.E23.r, c)):
            cov = _exp_apply_cov(rows, cov, x)
        return cov

    def covector_alt(self, x12, x13, x23):
        """phi = e_hw^T exp(x23 E23) exp(x13 E13) exp(x12 E12): the OTHER Gauss ordering, in which (x12, x13, x23) are exactly the matrix entries of the 3x3 unipotent
        n' (no cross terms: E23E12 = E23E13 = E13E12 = 0); equals the (a', g', c') ordering with x12 = a', x23 = c', x13 = g' + a'c'."""
        exact = isinstance(x12, Fr) and isinstance(x13, Fr) and isinstance(x23, Fr)
        cov = {self.hw: Fr(1) if exact else mp.mpf(1)}
        for rows, x in ((self.E23.r, x23), (self.E13.r, x13), (self.E12.r, x12)):
            cov = _exp_apply_cov(rows, cov, x)
        return cov

    def check_algebra(self):
        """[E13,E31] = E11 - E33 (E31 := [E32,E21]) on this rep -- needs the lowering operators; only for small reps"""
        from gl3_reps import _lower
        low = {}
        for k in (1, 2):
            rows = {}
            for col, pat in enumerate(self.patterns):
                for target, coeff in _lower(k, pat).items():
                    r = self.index.get(target)
                    if r is not None:
                        q = sp.Rational(coeff)
                        rows.setdefault(r, {})[col] = Fr(int(q.p), int(q.q))
            low[k] = SM(self.dim, rows)
        E21, E32 = low[1], low[2]
        E31 = comm(E32, E21)
        diag = SM(self.dim, {i: {i: Fr(self.weights[i][0] - self.weights[i][2])} for i in range(self.dim) if self.weights[i][0] != self.weights[i][2]})
        return (comm(self.E13, E31) - diag).is_zero() and comm(self.E12, E21).__sub__(
            SM(self.dim, {i: {i: Fr(self.weights[i][0] - self.weights[i][1])} for i in range(self.dim) if self.weights[i][0] != self.weights[i][1]})).is_zero()


_REP1 = {}
def rep1(lam):
    lam = tuple(int(x) for x in lam)
    if lam not in _REP1:
        _REP1[lam] = Rep1(lam)
    return _REP1[lam]


def gauss_covector(lam, phi, w1, w1p, idx=None, method='terms'):
    """phi_u = e1^T rho_Lambda(n') as a list over my GT index (zeros included). Exact Fractions for exact inputs, mpf for numeric inputs.
    idx: restrict to these GT indices (then the list is {nu: value}-dict instead). method: 'terms' (the Lambda-only tensor T; default) or 'sweep' (direct nilpotent exponentials)."""
    R = rep1(lam)
    a, g, c = gauss_params(phi, w1, w1p)
    if idx is not None:
        if method != 'terms':
            raise ValueError('idx only with method=terms')
        return R.covector_terms(a, g, c, idx)
    if method == 'terms':
        d = R.covector_terms(a, g, c)
    elif method == 'sweep':
        d = R.covector_sweep(a, g, c)
    else:
        raise ValueError(method)
    zero = Fr(0) if isinstance(a, Fr) and isinstance(g, Fr) and isinstance(c, Fr) else mp.mpf(0)
    return [d.get(i, zero) for i in range(R.dim)]


# ------------------------------------------------------------------------------------------------------------- symbolic polynomials
def _sym_gauss(phi):
    p1, p2, p3 = [sp.Rational(Fr(x).numerator, Fr(x).denominator) for x in phi]
    a = (p1 - p2) * (W1 - W2) / ((p2 - W1) * (p2 - W2))
    c = (p1 - p3) * (p2 - p3) * (W1 - W2) / ((p1 - p2) * (p3 - W1) * (p3 - W2))
    g = (p3 - p1) * (p3 - p2) * (W1 - W2) / ((p3 - W1) * (p2 - W2) * (p3 - W2))
    return a, g, c, (p1, p2, p3)


def chat_poly(lam, phi):
    l1, l2, l3 = lam
    _, _, _, (p1, p2, p3) = _sym_gauss(phi)
    f = lambda w: (p2 - w) ** (l1 - l2) * (p3 - w) ** l1
    return sp.Poly(sp.expand(f(W1) * f(W2) * (p1 - p2) ** (2 * l2)), W1, W2)


def phi_component_expr(lam, phi, nu):
    """phi_nu(w1, w2) as a sympy rational function in (W1, W2), from the Lambda-only tensor T (NOT multiplied by c-hat)"""
    R = rep1(lam)
    a, g, c, _ = _sym_gauss(phi)
    return sum(sp.Rational(t.numerator, t.denominator) * a ** p * g ** q * c ** r for (p, q, r, t) in R.T.get(nu, ()))


def gauss_covector_poly(lam, phi, comps=None):
    """c-hat * phi_nu as sympy Poly in (W1, W2) for every nu (index = my GT index), or None where c-hat does not clear the denominator
    (`untwist_limit.gf_polys` for the adjoint, generalised). Returns (list over nu, c-hat Poly)."""
    R = rep1(lam)
    ch = chat_poly(lam, phi)
    chx = ch.as_expr()
    out = []
    for nu in (range(R.dim) if comps is None else comps):
        ex = sp.cancel(sp.together(phi_component_expr(lam, phi, nu) * chx))
        num, den = sp.fraction(ex)
        out.append(sp.Poly(ex, W1, W2) if not den.has(W1, W2) else None)
    return out, ch


def uncleared_denominator(lam, phi, nu):
    """the leftover denominator (sympy factor list) of c-hat * phi_nu, or None if it is a polynomial"""
    ex = sp.cancel(sp.together(phi_component_expr(lam, phi, nu) * chat_poly(lam, phi).as_expr()))
    num, den = sp.fraction(ex)
    if not den.has(W1, W2):
        return None
    return sp.factor_list(den)


def _lin_pow_coeffs(c0, c1, k, K):
    """coefficients of t^0..t^(K-1) in (c0 + c1 t)^k (exact)"""
    out = []
    binom = 1
    for j in range(K):
        if j > k:
            out.append(Fr(0))
            continue
        out.append(binom * c0 ** (k - j) * c1 ** j)
        binom = binom * (k - j) // (j + 1)
    return out


def _poly_mul_trunc(a, b, K):
    out = [Fr(0)] * K
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            if i + j >= K:
                break
            out[i + j] += x * y
    return out


def pole_order(lam, phi, nu, w1=Fr(7, 3)):
    """order of the pole of c-hat * phi_nu at omega1' = phi2 (the ONLY place where c-hat can fail to clear a denominator; see RESULTS), exactly, from the
    Lambda-only tensor T. Structure: with n1 = p+q, n2 = q+r (fixed by the weight Lambda - nu = n1 a1 + n2 a2),
      c-hat phi_nu = sum_terms T K (W1-W2)^(p+q+r) (phi2-W1)^(l1-l2-p) (phi3-W1)^(l1-q-r) (phi2-W2)^(l1-l2-n1) (phi3-W2)^(l1-q-r),
    K = (phi1-phi2)^(2 l2) A^p Gg^q Cc^r, A = phi1-phi2, Gg = (phi3-phi1)(phi3-phi2), Cc = (phi1-phi3)(phi2-phi3)/(phi1-phi2);
    all exponents are >= 0 except (phi2-W2)^(l1-l2-n1) (p <= l1-l2 and n2 <= l1 hold automatically), so the pole order is
    max(0, n1-(l1-l2) - v), v = vanishing order at W2 = phi2 of the remaining sum (Taylor expansion in t = W2-phi2 at the generic rational W1)."""
    R = rep1(lam)
    l1, l2, l3 = lam
    terms = R.T.get(nu, ())
    if not terms:
        return 0
    p1, p2, p3 = [Fr(x) for x in phi]
    n1 = terms[0][0] + terms[0][1]
    excess = n1 - (l1 - l2)
    if excess <= 0:
        return 0
    A = p1 - p2
    Gg = (p3 - p1) * (p3 - p2)
    Cc = (p1 - p3) * (p2 - p3) / (p1 - p2)
    w1 = Fr(w1)
    for cand in (Fr(11, 5), Fr(13, 7), Fr(17, 3)):                      # the generic evaluation point must avoid the phi_i
        if w1 in (p1, p2, p3):
            w1 = cand
    K = excess
    tot = [Fr(0)] * K
    for (p, q, r, T) in terms:
        const = T * A ** (2 * l2) * A ** p * Gg ** q * Cc ** r * (p2 - w1) ** (l1 - l2 - p) * (p3 - w1) ** (l1 - q - r)
        poly = _poly_mul_trunc(_lin_pow_coeffs(w1 - p2, Fr(-1), p + q + r, K), _lin_pow_coeffs(p3 - p2, Fr(-1), l1 - q - r, K), K)
        for i in range(K):
            tot[i] += const * poly[i]
    v = next((i for i in range(K) if tot[i] != 0), K)
    return max(0, excess - v)


def uncleared_components(lam, phi, idx=None):
    """{nu: pole order at omega1' = phi2} for the components where c-hat does not clear the denominator"""
    R = rep1(lam)
    out = {}
    for nu in (range(R.dim) if idx is None else idx):
        o = pole_order(lam, phi, nu)
        if o > 0:
            out[nu] = o
    return out


# ------------------------------------------------------------------------------------------------------------- exact linear algebra (small)
def nullspace_exact(M, ncols):
    return nullspace_exact_free(M, ncols)[0]


def nullspace_exact_free(M, ncols):
    """right nullspace of the Fraction matrix M (list of rows), RREF form: one basis vector per free column (increasing column order),
    entry 1 at the free column, minus the RREF entries at the pivot columns. Returns list of coefficient lists (length ncols)."""
    A = [[Fr(x) for x in row] for row in M if any(x != 0 for x in row)]
    piv = []
    r = 0
    for c in range(ncols):
        pr = None
        for i in range(r, len(A)):
            if A[i][c] != 0:
                pr = i; break
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c)
        r += 1
        if r == len(A):
            break
    free = [c for c in range(ncols) if c not in piv]
    out = []
    for f in free:
        v = [Fr(0)] * ncols
        v[f] = Fr(1)
        for i, pc in enumerate(piv):
            v[pc] = -A[i][f]
        out.append(v)
    return out, free


def rank_exact(rows, ncols):
    A = [[Fr(x) for x in row] for row in rows if any(x != 0 for x in row)]
    return ncols - len(nullspace_exact(A, ncols)) if A else 0


# ------------------------------------------------------------------------------------------------------------- the PRV slice space
def weight_space(lam, nu):
    return list(rep1(lam).wspace.get(tuple(nu), []))


def prv_matrix(lam, mu):
    """(idx, rows): the PRV constraints on the weight space W_mu = V_Lambda[mu - Lambda]: rows of the (exact) matrix of (E12^(l1-l2+1), E23^(l2-l3+1)) restricted to W_mu"""
    R = rep1(lam)
    nu = tuple(int(mu[i] - lam[i]) for i in range(3))
    idx = R.wspace.get(nu, [])
    k1 = lam[0] - lam[1] + 1
    k2 = lam[1] - lam[2] + 1
    cons = {}
    for pos, j in enumerate(idx):
        for which, cols, k in ((0, R.C12, k1), (1, R.C23, k2)):
            ket = {j: Fr(1)}
            for _ in range(k):
                ket = _apply_ket(cols, ket)
                if not ket:
                    break
            for r_, v in ket.items():
                cons.setdefault((which, r_), {})[pos] = v
    rows = [[row.get(pos, Fr(0)) for pos in range(len(idx))] for row in cons.values()]
    return idx, rows


_PRV = {}
_PRV_FREE = {}
def prv_basis(lam, mu):
    """exact basis of M_mu = {v in V_Lambda[mu-Lambda]: E12^(l1-l2+1) v = E23^(l2-l3+1) v = 0}: list of {my GT index: Fraction}, RREF nullspace form
    (one vector per free column of the constraint matrix, free columns in increasing my-index order, entry 1 at its free index)."""
    key = (tuple(lam), tuple(mu))
    if key not in _PRV:
        idx, rows = prv_matrix(lam, mu)
        if not idx:
            _PRV[key] = []
            _PRV_FREE[key] = []
        else:
            ns, free = nullspace_exact_free(rows, len(idx))
            _PRV[key] = [{idx[pos]: v for pos, v in enumerate(vec) if v != 0} for vec in ns]
            _PRV_FREE[key] = [idx[f] for f in free]
    return _PRV[key]


def _clear_prv_cache(lam):
    for k in [k for k in _PRV if k[0] == tuple(lam)]:
        _PRV.pop(k, None)
        _PRV_FREE.pop(k, None)


def _sgn_perm(p):
    s = 1
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


def multiplicity(lam, mu):
    """Racah-Speiser multiplicity of V_mu in V_Lambda (x) V_Lambda: sum_sigma sgn(sigma) dim V_Lambda[sigma.mu - Lambda], sigma.mu = sigma(mu + rho) - rho, rho = (2,1,0)"""
    R = rep1(lam)
    rho = (2, 1, 0)
    v = [mu[i] + rho[i] for i in range(3)]
    tot = 0
    for perm in itertools.permutations(range(3)):
        wv = [v[perm[i]] for i in range(3)]
        sh = tuple(wv[i] - rho[i] - lam[i] for i in range(3))
        tot += _sgn_perm(perm) * len(R.wspace.get(sh, ()))
    return tot


def dim_irrep(mu):
    return (mu[0] - mu[1] + 1) * (mu[1] - mu[2] + 1) * (mu[0] - mu[2] + 2) // 2


def irreps(lam):
    """all dominant mu with multiplicity m >= 1 in V_Lambda (x) V_Lambda: list of (mu, m), sorted decreasing (dominance-ish lexicographic)"""
    R = rep1(lam)
    out = []
    for nu in R.wspace:
        mu = tuple(lam[i] + nu[i] for i in range(3))
        if mu[0] >= mu[1] >= mu[2]:
            m = multiplicity(lam, mu)
            if m > 0:
                out.append((mu, m))
    out.sort(reverse=True)
    return out


# ------------------------------------------------------------------------------------------------------------- the dictionary
_PTS = [Fr(3, 4), Fr(-5, 3), Fr(1, 2), Fr(-7, 3), Fr(9, 7), Fr(-4, 9), Fr(11, 5), Fr(-13, 6), Fr(-12, 5), Fr(7, 2), Fr(1, 8), Fr(5, 4),
        Fr(-8, 3), Fr(13, 9), Fr(-3, 7), Fr(17, 6), Fr(2, 9), Fr(-9, 8), Fr(19, 7), Fr(-6, 11), Fr(15, 4), Fr(-17, 10), Fr(5, 13), Fr(-21, 8)]


def default_points(phi, n, offset=0):
    """n deterministic generic rational pairs (w1, w1p): w1 != w1p, no coincidence with phi_i (poles of the Gauss data). Pair k = (P[k'], P[(7k'+3) mod len])
    with k' = k + offset over the fixed list P of small-height rationals; extends by a second pass with shifted partners when n exceeds the list."""
    ph = [x if isinstance(x, Fr) else (Fr(int(x)) if isinstance(x, int) else Fr(float(x))) for x in phi]
    out = []
    L = len(_PTS)
    k = offset
    shift = 3
    while len(out) < n:
        w1 = _PTS[k % L]
        w2 = _PTS[(7 * k + shift + (k // L)) % L]
        k += 1
        if w1 == w2 or w1 in ph or w2 in ph or (w1, w2) in out:
            continue
        out.append((w1, w2))
    return out


def _fix_scale(vec, order=None, tol=None, pivot=None):
    """divide by the first nonzero entry (|x| > tol * max|x|) in `order` (default: my index order); returns (vector, pivot index)"""
    n = len(vec)
    if tol is None:
        tol = mp.mpf(10) ** (-(mp.mp.dps // 2))
    if pivot is None:
        mx = max(abs(x) for x in vec)
        for i in (range(n) if order is None else order):
            if abs(vec[i]) > tol * mx:
                pivot = i
                break
    p = vec[pivot]
    return [x / p for x in vec], pivot


def dictionary_matrix(lam, mu, phi, points, basis=None):
    """A[k][j] = sum_nu [phi_u(w_k)]_nu b_j[nu] (exact Fractions for exact inputs, mpf else) and the c-hat(w_k) list"""
    R = rep1(lam)
    basis = prv_basis(lam, mu) if basis is None else basis
    nu = tuple(int(mu[i] - lam[i]) for i in range(3))
    idx = R.wspace.get(nu, [])
    A = []
    chs = []
    for (w1, w1p) in points:
        a, g, c = gauss_params(phi, w1, w1p)
        cov = R.covector_terms(a, g, c, idx)
        exact = isinstance(a, Fr)
        row = []
        for b in basis:
            s = Fr(0) if exact else mp.mpf(0)
            for j, bv in b.items():
                s += cov[j] * (bv if exact else to_mp(bv))
            row.append(s)
        A.append(row)
        chs.append(chat(lam, phi, w1, w1p))
    return A, chs


def _is_exact_num(x):
    return isinstance(x, (int, Fr, sp.Rational))


def _solve_exact_ls(A, y):
    """exact rational least squares (normal equations, Gauss-Jordan) for the Fraction matrix A (list of rows) and vector y; returns (x, squared relative residual as Fraction)"""
    m = len(A[0])
    M = [[sum(A[k][i] * A[k][j] for k in range(len(A))) for j in range(m)] + [sum(A[k][i] * y[k] for k in range(len(A)))] for i in range(m)]
    for c in range(m):
        pr = next((i for i in range(c, m) if M[i][c] != 0), None)
        if pr is None:
            raise ValueError('singular dictionary matrix (the points do not separate the PRV basis)')
        M[c], M[pr] = M[pr], M[c]
        pv = M[c][c]
        M[c] = [v / pv for v in M[c]]
        for i in range(m):
            if i != c and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[c])]
    x = [M[i][m] for i in range(m)]
    r2 = sum((sum(A[k][j] * x[j] for j in range(m)) - y[k]) ** 2 for k in range(len(A)))
    y2 = sum(v * v for v in y)
    return x, (r2 / y2 if y2 != 0 else Fr(0))


def prv_free_indices(lam, mu):
    """the GT index f_j carried by each PRV basis vector (b_j[f_j] = 1, b_k[f_j] = 0 for k != j): the coordinates of ANY v in M_mu in this basis are x_j = v[f_j]"""
    prv_basis(lam, mu)
    return list(_PRV_FREE[(tuple(lam), tuple(mu))])


def coords_of(lam, mu, vec):
    """coordinates of a vector of M_mu (list over my GT index, or dict) in prv_basis(lam, mu); checks membership (returns (x, max |v - sum x_j b_j|))"""
    basis = prv_basis(lam, mu)
    fr = prv_free_indices(lam, mu)
    get = (lambda i: vec.get(i, 0)) if isinstance(vec, dict) else (lambda i: vec[i])
    x = [get(f) for f in fr]
    n = rep1(lam).dim
    rec = {}
    for xj, b in zip(x, basis):
        for i, v in b.items():
            rec[i] = rec.get(i, 0) + xj * (v if isinstance(xj, Fr) else to_mp(v))
    err = max([abs(to_mp(get(i)) - to_mp(rec.get(i, 0))) for i in range(n)] + [mp.mpf(0)])
    return x, err


class DictionaryProblem:
    """The (m+1)-point dictionary of ONE irrep mu at FIXED points: A[k][j] = <phi_u(w_k)|b_j> (sum over W_mu only), c-hat(w_k).  All copies (branches) of the
    irrep share it: solve(values) gives (x, residual, slice) for any list of values L_B(w_k) (one list per copy).
    row_scale='norm': every equation is divided by ||A_k||_2 before the (QR) least-squares solve; residual = ||D(Ax - y)||/||D y||."""

    def __init__(self, lam, mu, phi, points, dps=None, basis=None, row_scale='norm'):
        self.lam = tuple(int(x) for x in lam)
        self.mu = tuple(int(x) for x in mu)
        self.phi = tuple(phi)
        self.points = [tuple(p) for p in points]
        self.dps = dps or max(mp.mp.dps, 30)          # never below 30 digits; pass dps explicitly for more (the precision of the DATA is not inferred)
        self.row_scale = row_scale
        t0 = time.time()
        with mp.workdps(self.dps):
            self.R = rep1(self.lam)
            self.basis = prv_basis(self.lam, self.mu) if basis is None else basis
            self.m = len(self.basis)
            if self.m == 0:
                raise ValueError('M_mu = 0 for mu = %s in V_%s x V_%s' % (self.mu, self.lam, self.lam))
            if len(self.points) < self.m:
                raise ValueError('need at least m = %d points (got %d)' % (self.m, len(self.points)))
            for (w1, w1p) in self.points:
                if w1 == w1p:
                    raise ValueError("degenerate point w1 == w1p = %s (n_prime = 1: the covector is e_hw)" % (w1,))
                for p_ in self.phi:
                    if w1 == p_ or w1p == p_:
                        raise ValueError('point (%s, %s) coincides with a phi_i = %s (pole of the Gauss data)' % (w1, w1p, p_))
            self.A, self.chat = dictionary_matrix(self.lam, self.mu, self.phi, self.points, self.basis)
            self.exact_A = all(isinstance(a, Fr) for row in self.A for a in row) and all(isinstance(c, Fr) for c in self.chat)
            n = len(self.points)
            self.Am = mp.matrix(n, self.m)
            for k in range(n):
                for j in range(self.m):
                    self.Am[k, j] = to_mp(self.A[k][j])
            if row_scale == 'norm':
                self.d = [1 / mp.sqrt(mp.fsum(self.Am[k, j] ** 2 for j in range(self.m))) for k in range(n)]
            else:
                self.d = [mp.mpf(1)] * n
            self.As = mp.matrix(n, self.m)
            for k in range(n):
                for j in range(self.m):
                    self.As[k, j] = self.Am[k, j] * self.d[k]
        self.t_build = time.time() - t0

    def cond(self):
        """cond(D A) (2-norm) at the working precision"""
        with mp.workdps(self.dps):
            sv = mp.svd_r(self.As, compute_uv=False)
            svl = [sv[i] for i in range(len(sv))]
            return svl[0] / svl[-1] if svl[-1] != 0 else mp.inf

    def solve(self, values, pivot_tol=None, exact=None, pivot=None):
        """-> dict(x, residual, slice (normalised), raw (= sum x_j b_j, scale 1/X_B), pivot, x_exact, exact)"""
        n = len(self.points)
        if len(values) != n:
            raise ValueError('need one value per point (%d points, %d values)' % (n, len(values)))
        t0 = time.time()
        with mp.workdps(self.dps):
            exact_ok = self.exact_A and all(_is_exact_num(v) for v in values)
            if exact is None:
                exact = exact_ok
            if exact and not exact_ok:
                raise ValueError('exact=True needs exact phi, points and values')
            x_exact = None
            if exact:
                yq = [_num(v) / self.chat[k] for k, v in enumerate(values)]
                x_exact, r2 = _solve_exact_ls(self.A, yq)
                x = [to_mp(v) for v in x_exact]
                res = mp.sqrt(to_mp(r2))
            else:
                y = [to_mp(_num(v)) / to_mp(self.chat[k]) for k, v in enumerate(values)]
                ys = mp.matrix([y[k] * self.d[k] for k in range(n)])
                if n == self.m:
                    xs = mp.lu_solve(self.As, ys)
                    res = mp.mpf(0)
                else:
                    xs, _r = mp.qr_solve(self.As, ys)
                    res = mp.norm(self.As * xs - ys) / mp.norm(ys)
                x = [xs[j] for j in range(self.m)]
            raw = [mp.mpf(0)] * self.R.dim
            for j, b in enumerate(self.basis):
                for i, bv in b.items():
                    raw[i] += x[j] * to_mp(bv)
            vec, piv = _fix_scale(raw, tol=pivot_tol, pivot=pivot)
        return dict(x=x, x_exact=x_exact, residual=res, slice=vec, raw=raw, pivot=piv, exact=bool(exact), t_solve=time.time() - t0)


def slice_from_points_ex(lam, mu, phi, points, values, dps=None, pivot_tol=None, want_sv=False, basis=None, row_scale='norm', exact=None, pivot=None):
    """the (m+1)-point dictionary, full output (see slice_from_points).
    Equations (one per point k): sum_j x_j A[k][j] = y_k, A[k][j] = <phi_u(w_k)|b_j> over W_mu, y_k = values[k]/c-hat(w_k).
    Numeric path: least squares by QR in mpmath after scaling every equation by 1/||A_k||_2 (row_scale='norm'; None = no scaling); residual = ||D(Ax-y)||/||D y||.
    Exact path (all of phi, points, values exact rationals; exact=None/True): exact rational least squares, residual from the exact squared residual."""
    if len(points) != len(values):
        raise ValueError('need as many values as points (got %d points, %d values)' % (len(points), len(values)))
    t0 = time.time()
    P = DictionaryProblem(lam, mu, phi, points, dps=dps, basis=basis, row_scale=row_scale)
    out = P.solve(values, pivot_tol=pivot_tol, exact=exact, pivot=pivot)
    out.update(m=P.m, A=P.Am, chat=P.chat, npts=len(points), t_rows=P.t_build, t_total=time.time() - t0)
    if want_sv:
        with mp.workdps(P.dps):
            sv = mp.svd_r(P.As, compute_uv=False)
            svl = [sv[i] for i in range(len(sv))]
        out['sv'] = svl
        out['cond'] = svl[0] / svl[-1] if svl[-1] != 0 else mp.inf
    return out


def slice_from_points(lam, mu, phi, points, values, dps=None):
    """THE INTERFACE: (x, residual, slice_vector) -- see the module docstring"""
    o = slice_from_points_ex(lam, mu, phi, points, values, dps=dps)
    return o['x'], o['residual'], o['slice']


def slices_from_points(lam, mu, phi, points, values_per_branch, dps=None):
    """all copies of one irrep at once (one dictionary matrix): values_per_branch = [values of copy 1, values of copy 2, ...] -> [(x, residual, slice_vector), ...]"""
    P = DictionaryProblem(lam, mu, phi, points, dps=dps)
    out = []
    for vals in values_per_branch:
        o = P.solve(vals)
        out.append((o['x'], o['residual'], o['slice']))
    return out


_SHARED = {}
def shared_points(lam, phi, pool_n=None, mmin=None):
    """ONE ordered list of points for ALL irreps of V_Lambda (x) V_Lambda: the pool of `pool_n` default points ordered by column-pivoted QR (float64) of the row-normalised dictionary
    matrices of the irreps with m >= mmin (default max(2, m_max - 2): the hardest ones) stacked side by side; use points = shared_points(lam, phi)[:m+1] for an irrep of
    multiplicity m (nested prefixes: the kernel pieces G^[d](w_k) are then needed at the first m_max+1 points only). Measured worst cond(A) over ALL irreps: see RESULTS."""
    import numpy as np
    from scipy.linalg import qr
    irr = irreps(lam)
    mmax = max(m for _, m in irr)
    if pool_n is None:
        pool_n = max(40, 3 * (mmax + 1))
    key = (tuple(lam), tuple(Fr(x) if isinstance(x, (int, Fr)) else float(x) for x in phi), pool_n, mmin)
    if key not in _SHARED:
        pool = default_points(phi, pool_n)
        lo = max(2, mmax - 2) if mmin is None else mmin
        mats = []
        for mu, m in irr:
            if m >= lo:
                A, ch = dictionary_matrix(lam, mu, phi, pool)
                M = np.array([[float(v) for v in row] for row in A], dtype=float)
                M = M / np.linalg.norm(M, axis=1)[:, None]
                mats.append(M)
        if not mats:
            _SHARED[key] = pool
        else:
            _, _, P = qr(np.hstack(mats).T, pivoting=True)
            _SHARED[key] = [pool[i] for i in P]
    return list(_SHARED[key])


def synthetic_values(lam, mu, phi, points, slice_vec):
    """L(w_k) = c-hat(w_k) <phi_u(w_k)|slice> for a slice given as {my GT index: number} or a list over my GT index (the identity itself; exact for exact input)"""
    R = rep1(lam)
    sd = slice_vec if isinstance(slice_vec, dict) else {i: v for i, v in enumerate(slice_vec) if v != 0}
    sd = {i: _num(v) for i, v in sd.items()}
    out = []
    for (w1, w1p) in points:
        a, g, c = gauss_params(phi, w1, w1p)
        cov = R.covector_terms(a, g, c, list(sd))
        ch = chat(lam, phi, w1, w1p)
        if isinstance(a, Fr) and all(isinstance(v, Fr) for v in sd.values()):
            out.append(ch * sum(cov[i] * v for i, v in sd.items()))
        else:
            out.append(to_mp(ch) * mp.fsum(to_mp(cov[i]) * to_mp(v) for i, v in sd.items()))
    return out


def choose_points(lam, mu, phi, n, pool=None, basis=None):
    """greedy well-conditioned choice of n points out of a pool (default: 40 `default_points`) by column-pivoted QR (float64) on the row-normalised
    dictionary matrix; returns the list of chosen points."""
    import numpy as np
    from scipy.linalg import qr
    pool = pool or default_points(phi, 40)
    A, chs = dictionary_matrix(lam, mu, phi, pool, basis)
    M = np.array([[float(v) for v in row] for row in A], dtype=float)
    nr = np.linalg.norm(M, axis=1)
    nr[nr == 0] = 1.0
    M = M / nr[:, None]
    _, _, P = qr(M.T, pivoting=True)
    return [pool[i] for i in P[:n]]


# ============================================================================================================== validation (lazy imports)
PHI = (Fr(1, 3), Fr(-2, 7), Fr(5, 11))       # N1/U2's two directions
PHI2 = (Fr(-2, 5), Fr(3, 8), Fr(7, 11))


def _gt_str(pat):
    return '%d|%d,%d' % (pat[0][0], pat[1][0], pat[1][1])


def sym_pole_order(lam, phi, nu):
    """sympy truth: order of the leftover (W2 - phi2) in the denominator of c-hat phi_nu (0 if cleared), or None if other factors survive"""
    fl = uncleared_denominator(lam, phi, nu)
    if fl is None:
        return 0
    p2 = sp.Rational(Fr(phi[1]).numerator, Fr(phi[1]).denominator)
    tot = 0
    for f, e in fl[1]:
        Pf = sp.Poly(f, W1, W2)
        if Pf.degree(W1) == 0 and Pf.degree(W2) == 1 and sp.simplify(f.subs(W2, p2)) == 0:
            tot += e
        else:
            return None
    return tot


def _w4_slices(lam):
    """all exact W4 slices at site 2 = HW, in MY one-site index: list of (mu, tag, {my index: value}) for every HWBASIS ket (rational) and every BRANCH (mp, G-normalised)"""
    import nu_general as N
    R = N.rep_of(lam)
    d = R.d
    hw = R.hw
    out = []
    for mu in R.w4['order']:
        for k, v in enumerate(N.hw_kets(R, mu)):
            out.append((mu, 'hw%d' % (k + 1), {nu: v[nu * d + hw] for nu in range(d) if v[nu * d + hw] != 0}))
        for b in R.w4['sectors'][mu]['branches']:
            cg = N.w4_cg(R, mu, b['copy'])
            out.append((mu, 'br%d' % b['copy'], {nu: x for nu, x in enumerate(cg) if x != 0}))
    return R, out


def stage_polys():
    import untwist_limit as U
    from gaudin_chain import Chain
    import nu_general as N
    mp.mp.dps = 50
    print('=== D1(c): untwist_limit.gf_polys reproduced exactly (adjoint, sympy Poly equality component by component)')
    for phi in (PHI, PHI2):
        th, thc = U.gf_polys(Chain(lam=(2, 1, 0)), phi)
        mine, minec = gauss_covector_poly((2, 1, 0), phi)
        ok = all(((P is None) == (Q is None)) and (P is None or (P - Q).is_zero) for P, Q in zip(mine, th)) and (minec - thc).is_zero
        print('  phi=%s: %d components, polynomial on %d, None on %s, c-hat equal: %s -> reproduced: %s' % (
            [str(x) for x in phi], len(mine), sum(P is not None for P in mine), [nu for nu, P in enumerate(mine) if P is None], (minec - thc).is_zero, ok))
        assert ok
    print('=== D1(d): components where c-hat * phi_nu is NOT a polynomial; the leftover is (omega1p - phi2)^order; assertion: C_nu = 0 on every W4 HW ket')
    for lam in [(2, 1, 0), (3, 1, 0), (3, 2, 0)]:
        R, sl = _w4_slices(lam)
        R1 = rep1(lam)
        for phi in (PHI, PHI2):
            polys, _ = gauss_covector_poly(lam, phi)
            unc = [nu for nu, P in enumerate(polys) if P is None]
            fast = uncleared_components(lam, phi)
            symord = {nu: sym_pole_order(lam, phi, nu) for nu in unc}
            assert set(unc) == set(fast) and all(symord[nu] == fast[nu] for nu in unc), (lam, unc, fast, symord)
            cmax = mp.mpf(0)
            for (mu, tag, s) in sl:
                for nu in unc:
                    cmax = max(cmax, abs(to_mp(s.get(nu, 0))))
            print('  %s phi=%s: %d W4 slices (%d HWBASIS kets + branches); uncleared nu: %s; max |C_nu| over all of them = %s' % (
                lam, [str(x) for x in phi], len(sl), sum(1 for t in sl if t[1].startswith('hw')),
                [(nu, _gt_str(R1.patterns[nu]), R1.weights[nu], 'pole^%d' % fast[nu]) for nu in unc], mp.nstr(cmax, 3)))
            assert cmax == 0
    print('=== fast exact pole-order analysis (Lambda-only tensor T, Taylor expansion at omega1p = phi2) vs the sympy truth (cancel of c-hat * phi_nu): same uncleared set and same orders, both phi')
    for lam in [(4, 1, 0), (4, 2, 0), (4, 3, 0), (5, 2, 0)]:
        for phi in (PHI, PHI2):
            polys, _ = gauss_covector_poly(lam, phi)
            unc = [nu for nu, P in enumerate(polys) if P is None]
            fast = uncleared_components(lam, phi)
            symord = {nu: sym_pole_order(lam, phi, nu) for nu in unc}
            ok = set(unc) == set(fast) and all(symord[nu] == fast[nu] for nu in unc)
            print('  %s phi=%s: uncleared (sympy) %s ; fast %s ; agree: %s' % (lam, [str(x) for x in phi], {nu: symord[nu] for nu in unc}, fast, ok))
            assert ok
    print('=== structure: uncleared components vs the PRV supports (all irreps), fast exact pole-order analysis, both phi')
    for lam in [(2, 1, 0), (3, 1, 0), (3, 2, 0), (4, 2, 0), (6, 3, 0), (8, 4, 0), (14, 7, 0), (28, 14, 0)]:
        t0 = time.time()
        R1 = rep1(lam)
        uu = [uncleared_components(lam, phi) for phi in (PHI, PHI2)]
        sup = set()
        for mu, m in irreps(lam):
            for b in prv_basis(lam, mu):
                sup |= set(b)
        inter = sup & set(uu[0])
        gt = [nu for nu in range(R1.dim) if lam[0] - R1.weights[nu][0] > lam[0] - lam[1]]
        print('  %-10s dim %5d: #uncleared %4d (max pole order %2d; phi-independent: %s), #nu with n1 > l1-l2: %4d, |union of PRV supports| %4d, intersection with uncleared: %d  (%.1fs)' % (
            lam, R1.dim, len(uu[0]), max(uu[0].values()) if uu[0] else 0, uu[0] == uu[1], len(gt), len(sup), len(inter), time.time() - t0))
        assert not inter and uu[0] == uu[1]
    print('=== D1(e): gauss_covector (Fractions) vs N1 nu_general.Rep.phi_u (Chain-based, mp) and the other Gauss ordering (matrix-entry coordinates), exact')
    for lam in [(2, 1, 0), (3, 1, 0), (3, 2, 0)]:
        RN = N.rep_of(lam)
        devs = []
        for (w1, w2) in default_points(PHI, 4) + default_points(PHI2, 2):
            for phi in (PHI, PHI2):
                c1 = gauss_covector(lam, phi, w1, w2)
                c2 = RN.phi_u(phi, w1, w2)
                devs.append(max(abs(to_mp(a) - b) for a, b in zip(c1, c2)))
        print('  %s: max |gauss_covector - Rep.phi_u| over %d (phi, point) combinations = %s' % (lam, len(devs), mp.nstr(max(devs), 3)))
    for lam in [(2, 1, 0), (3, 1, 0), (3, 2, 0), (4, 2, 0), (8, 4, 0), (14, 7, 0)]:
        R1 = rep1(lam)
        okc = okd = True
        for (w1, w2) in default_points(PHI, 3):
            a, g, c = gauss_params(PHI, w1, w2)
            d1 = R1.covector_terms(a, g, c)
            d2 = R1.covector_sweep(a, g, c)
            d3 = R1.covector_alt(a, g + a * c, c)
            okc &= all(d1[i] == d2.get(i, 0) for i in range(R1.dim))
            okd &= all(d1[i] == d3.get(i, 0) for i in range(R1.dim))
        extra = R1.check_algebra() if R1.dim <= 130 else None
        print('  %-10s terms == sweep (exact): %s ; (a,g,c)-ordering == (x23,x13,x12)-ordering with x13 = g+ac (exact): %s ; [E13,E31] = E11-E33 and [E12,E21] = E11-E22: %s' % (lam, okc, okd, extra))
        assert okc and okd


def stage_timing(amax=14):
    print('=== D1(f): timing, Lambda = [2a,a,0] (dim = (a+1)^3); times in seconds (exact = Fractions, mp = 50 digits); medians over 5 points')
    print('  %3s %6s %8s %8s %9s %9s %9s %8s | %5s %6s %6s %7s' % ('a', 'dim', 'setup', 'T-tens', 'cov exact', 'cov mp50', 'chat', 'T size', '#irr', 'maxm', 'maxw', 'PRVall'))
    rows = []
    for a in range(1, amax + 1):
        lam = (2 * a, a, 0)
        t0 = time.time(); R = Rep1(lam); t_setup = time.time() - t0
        t0 = time.time(); T = R.T; t_T = time.time() - t0
        _REP1[lam] = R
        pts = default_points(PHI, 5)
        te, tm, tc = [], [], []
        mp.mp.dps = 50
        phim = tuple(to_mp(x) for x in PHI)
        for (w1, w2) in pts:
            t0 = time.time(); gauss_covector(lam, PHI, w1, w2); te.append(time.time() - t0)
            w1m, w2m = to_mp(w1), to_mp(w2)
            t0 = time.time(); gauss_covector(lam, phim, w1m, w2m); tm.append(time.time() - t0)
            t0 = time.time(); chat(lam, PHI, w1, w2); tc.append(time.time() - t0)
        med = lambda l: sorted(l)[len(l) // 2]
        _clear_prv_cache(lam)
        t0 = time.time()
        irr = irreps(lam)
        for mu, m in irr:
            prv_basis(lam, mu)
        t_prv = time.time() - t0
        maxw = max(len(v) for v in R.wspace.values())
        row = (a, R.dim, t_setup, t_T, med(te), med(tm), med(tc), sum(len(v) for v in T.values()), len(irr), max(m for _, m in irr), maxw, t_prv)
        rows.append(row)
        print('  %3d %6d %8.3f %8.3f %9.4f %9.4f %9.5f %8d | %5d %6d %6d %7.2f' % row, flush=True)
    return rows


def stage_prv():
    mp.mp.dps = 60
    print('=== D2: dim M_mu = Racah-Speiser multiplicity for every irrep; sum m dim V_mu = dim^2')
    for lam in [(2, 1, 0), (3, 1, 0), (3, 2, 0), (4, 2, 0), (6, 3, 0), (8, 4, 0), (14, 7, 0), (28, 14, 0)]:
        _clear_prv_cache(lam)
        t0 = time.time()
        R = rep1(lam)
        irr = irreps(lam)
        bad = [(mu, m, len(prv_basis(lam, mu))) for mu, m in irr if len(prv_basis(lam, mu)) != m]
        tot = sum(m * dim_irrep(mu) for mu, m in irr)
        dt = time.time() - t0
        print('  %-10s dim %5d: %3d irreps, %4d copies, max m %2d, max w %2d; dim M_mu == m for all: %s (mismatches %s); sum m dim V_mu = %d = dim^2: %s; PRV for all irreps %.2fs (%.1f ms/irrep)' % (
            lam, R.dim, len(irr), sum(m for _, m in irr), max(m for _, m in irr), max(len(v) for v in R.wspace.values()), not bad, bad, tot, tot == R.dim ** 2, dt, 1e3 * dt / len(irr)))
        assert not bad and tot == R.dim ** 2
    print('=== D2: W4 exact slices (every HWBASIS ket and every BRANCH slice) satisfy the PRV conditions and span M_mu (exact ranks; surd branches checked in mp)')
    for lam in [(2, 1, 0), (3, 1, 0), (3, 2, 0)]:
        R, sl = _w4_slices(lam)
        R1 = rep1(lam)
        k1 = lam[0] - lam[1] + 1
        k2 = lam[1] - lam[2] + 1
        for mu in R.w4['order']:
            basis = prv_basis(lam, mu)
            idx = weight_space(lam, tuple(mu[i] - lam[i] for i in range(3)))
            mine = [s for (m_, tag, s) in sl if m_ == mu and tag.startswith('hw')]
            viol = 0
            for s in mine:
                for cols, k in ((R1.C12, k1), (R1.C23, k2)):
                    ket = dict(s)
                    for _ in range(k):
                        ket = _apply_ket(cols, ket)
                    viol += sum(1 for v in ket.values() if v != 0)
            outside = sum(1 for s in mine for j in s if j not in idx)
            rs = rank_exact([[s.get(j, 0) for j in idx] for s in mine], len(idx))
            rb = rank_exact([[b.get(j, 0) for j in idx] for b in basis], len(idx))
            rsb = rank_exact([[s.get(j, 0) for j in idx] for s in mine] + [[b.get(j, 0) for j in idx] for b in basis], len(idx))
            bres = mp.mpf(0)
            for (m_, tag, s) in sl:
                if m_ == mu and tag.startswith('br'):
                    A = mp.matrix(len(idx), len(basis))
                    for j, b in enumerate(basis):
                        for pos, i in enumerate(idx):
                            A[pos, j] = to_mp(b.get(i, 0))
                    y = mp.matrix([to_mp(s.get(i, 0)) for i in idx])
                    if len(basis) == len(idx):
                        x = mp.lu_solve(A, y); r = mp.mpf(0)
                    else:
                        x, r = mp.qr_solve(A, y)
                    bres = max(bres, r / mp.norm(y))
            print('  %s mu=%s: m=%d  HW-ket slices %d: PRV violations %d, outside W_mu %d, rank(slices)=%d rank(basis)=%d rank(both)=%d; branch slices in span(M_mu): rel. residual %s' % (
                lam, mu, len(basis), len(mine), viol, outside, rs, rb, rsb, mp.nstr(bres, 3)))
            assert viol == 0 and outside == 0 and rs == rb == rsb == len(basis)


# ------------------------------------------------------------------------------------------------------------- D3: the 26 known branches
W4REPS = [(2, 1, 0), (3, 1, 0), (3, 2, 0)]


def _load_records(lam, tag):
    fn = os.path.join(HERE, 'nu_general_res_%d%d%d%s.pkl' % (tuple(lam) + (tag,)))
    with open(fn, 'rb') as f:
        return pickle.load(f)


def _pe(LB, w1, w2):
    w1, w2 = to_mp(Fr(w1)), to_mp(Fr(w2))
    return mp.fsum(c * w1 ** a * w2 ** b for (a, b), c in LB.items())


def _first_nz(v, tol=mp.mpf(10) ** -20):
    return next(i for i, x in enumerate(v) if abs(x) > tol)


def slice_dev(sl, ex):
    """max |sl/sl[piv] - ex/ex[piv]|, piv = first nonzero of the exact slice ex (my order); also the max of sl outside supp(ex), normalised"""
    piv = _first_nz(ex)
    a = [x / sl[piv] for x in sl]
    b = [x / ex[piv] for x in ex]
    return max(abs(x - y) for x, y in zip(a, b))


def _digits(x):
    """number of correct digits (floor(-log10 x)); 'exact' if x == 0"""
    if x == 0:
        return 'exact'
    return int(mp.floor(-mp.log10(abs(x))))


def stage_branches(dps=120):
    import nu_general as N
    mp.mp.dps = dps
    grid = [(a, b) for a in N.GRID8_1 for b in N.GRID8_2]                   # N1's 64-point grid (held-out test points)
    print('=== D3: the (m+1)-point dictionary on the leading polynomials L_B of N1 (records nu_general_res_<lam>[2].pkl, 200-digit data); dps = %d' % dps)
    print('    per branch record: rep, phi (1/2), mu, copy | m, w | m+1 points: residual, correct digits of the slice | m+2 pts | 2m+2 pts | abs. scale max|raw*X_B - exact slice|/max | held-out: L_B at the 64 grid points of N1 vs prediction c-hat<phi_u|raw> | identity-only: digits | cond(A)')
    print('    (correct digits: componentwise vs the exact W4 slice after normalising both on the first nonzero entry of the exact slice; "exact" = identical at working precision; for m = 1 the slice is the PRV vector itself)')
    allrows = []
    for lam in W4REPS:
        R = N.rep_of(lam)
        for tag in ('', '2'):
            o = _load_records(lam, tag)
            phi = tuple(Fr(x) for x in o['phi'])
            for rec in o['br']:
                mu, cp = rec['mu'], rec['copy']
                basis = prv_basis(lam, mu)
                m = len(basis)
                w = len(weight_space(lam, tuple(mu[i] - lam[i] for i in range(3))))
                ex = N.w4_cg(R, mu, cp)                                    # exact G-normalised W4 slice, my index
                assert max(abs(a - b) for a, b in zip(ex, rec['slex'])) < mp.mpf(10) ** -90
                row = dict(lam=lam, tag=tag, mu=mu, copy=cp, m=m, w=w)
                for key, npts in (('a', m + 1), ('b', m + 2), ('c', 2 * m + 2)):
                    pts = default_points(phi, npts)
                    vals = [_pe(rec['LB'], *p) for p in pts]
                    o1 = slice_from_points_ex(lam, mu, phi, pts, vals, dps=dps, want_sv=(key == 'a'))
                    row['res_' + key] = o1['residual']
                    row['dev_' + key] = slice_dev(o1['slice'], ex)
                    if key == 'a':
                        row['cond'] = o1['cond']
                        row['scale_dev'] = max(abs(o1['raw'][nu] * rec['XB'] - rec['slex'][nu]) for nu in range(R.d)) / max(abs(x) for x in rec['slex'])
                        Lh = [_pe(rec['LB'], *p) for p in grid]
                        Ph = synthetic_values(lam, mu, phi, grid, o1['raw'])
                        row['held'] = max(abs(a - b) for a, b in zip(Lh, Ph)) / max(abs(a) for a in Lh)
                # identity-only: L from the exact slice through the identity
                pts = default_points(phi, m + 1)
                vals = synthetic_values(lam, mu, phi, pts, ex)
                o2 = slice_from_points_ex(lam, mu, phi, pts, vals, dps=dps)
                row['dev_id'] = slice_dev(o2['slice'], ex)
                row['res_id'] = o2['residual']
                allrows.append(row)
                print('  %s %d %s #%d | m=%d w=%d | res %s dig %s | dig %s | dig %s | scale %s | held-out %s | id-only dig %s | cond %s' % (
                    lam, 1 if tag == '' else 2, mu, cp, m, w, mp.nstr(row['res_a'], 2), _digits(row['dev_a']), _digits(row['dev_b']), _digits(row['dev_c']),
                    mp.nstr(row['scale_dev'], 2), mp.nstr(row['held'], 2), _digits(row['dev_id']), mp.nstr(row['cond'], 3)), flush=True)
    n = len(allrows)
    print('  SUMMARY: %d branch records (26 branches x 2 phi); worst slice error: m+1 points %s, m+2 %s, 2m+2 %s; worst residual (m+1) %s; worst abs-scale deviation %s; worst held-out prediction error %s; identity-only worst %s' % (
        n, mp.nstr(max(r['dev_a'] for r in allrows), 3), mp.nstr(max(r['dev_b'] for r in allrows), 3), mp.nstr(max(r['dev_c'] for r in allrows), 3),
        mp.nstr(max(r['res_a'] for r in allrows), 3), mp.nstr(max(r['scale_dev'] for r in allrows), 3), mp.nstr(max(r['held'] for r in allrows), 3), mp.nstr(max(r['dev_id'] for r in allrows), 3)))
    for lam in W4REPS:
        sub = [r for r in allrows if r['lam'] == lam]
        print('    %s: %d records, worst residual %s, worst slice error %s, worst held-out %s' % (lam, len(sub), mp.nstr(max(r['res_a'] for r in sub), 3), mp.nstr(max(r['dev_a'] for r in sub), 3), mp.nstr(max(r['held'] for r in sub), 3)))
    m2 = [r for r in allrows if r['m'] == 2]
    print('  multiplicity-2 sectors: %d branch records (both copies of each sector, both phi): worst slice error (m+1 points) %s, worst cond(A) %s; every copy matches ITS OWN exact slice (the two copies are different vectors of the 2-dim M_mu)' % (
        len(m2), mp.nstr(max(r['dev_a'] for r in m2), 3), mp.nstr(max(r['cond'] for r in m2), 3)))
    return allrows


# ------------------------------------------------------------------------------------------------------------- conditioning and noise
def _rand_slice(lam, mu, rng):
    """random exact element of M_mu: (coefficients in the PRV basis, slice dict)"""
    basis = prv_basis(lam, mu)
    cs = [Fr(rng.choice([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5])) for _ in basis]
    s = {}
    for c, b in zip(cs, basis):
        for i, v in b.items():
            s[i] = s.get(i, 0) + c * v
    return cs, {i: v for i, v in s.items() if v != 0}


def _noisy(vals, delta, rng):
    return [v * (1 + delta * to_mp(Fr(rng.randint(-10 ** 6, 10 ** 6), 10 ** 6))) for v in vals]


def _dense(sd, n):
    v = [mp.mpf(0)] * n
    for i, x in sd.items():
        v[i] = to_mp(x)
    return v


def dir_err(u, v):
    """scale-free direction error of u relative to v: ||u - c v||/||u||, c = (u.v)/(v.v) (Euclidean norm in the GT coordinates)"""
    uv = mp.fsum(a * b for a, b in zip(u, v))
    vv = mp.fsum(b * b for b in v)
    c = uv / vv
    return mp.sqrt(mp.fsum((a - c * b) ** 2 for a, b in zip(u, v)) / mp.fsum(a * a for a in u))


def stage_cond(dps=60):
    import random
    import nu_general as N
    mp.mp.dps = dps
    rng = random.Random(20261021)
    print('=== D3(i): conditioning of the dictionary matrix (rows scaled to unit norm) with the default points vs greedily chosen points (column-pivoted QR over a pool of 40), real branches (copy 1 of each sector, phi)')
    print('    %-10s %-10s %2s %2s | cond default (m+1 pts) | cond chosen (m+1 pts) | cond default (2m+2 pts)' % ('rep', 'mu', 'm', 'w'))
    for lam in W4REPS:
        o = _load_records(lam, '')
        phi = tuple(Fr(x) for x in o['phi'])
        for rec in o['br']:
            if rec['copy'] != 1:
                continue
            mu = rec['mu']
            m = len(prv_basis(lam, mu))
            if m == 1:
                continue
            w = len(weight_space(lam, tuple(mu[i] - lam[i] for i in range(3))))
            row = []
            for pts in (default_points(phi, m + 1), choose_points(lam, mu, phi, m + 1), default_points(phi, 2 * m + 2)):
                row.append(DictionaryProblem(lam, mu, phi, pts, dps=dps).cond())
            print('    %-10s %-10s %2d %2d | %-22s | %-21s | %s' % (lam, mu, m, w, mp.nstr(row[0], 3), mp.nstr(row[1], 3), mp.nstr(row[2], 3)))
    print('    (m = 1: cond = 1 trivially; the slice is the PRV vector itself and the data only fix the scale and test consistency)')
    print('=== D3(ii): noise amplification on the real branches with m = 2: relative noise delta on every L_B value (uniform in [-delta, delta], 4 trials); amplification = worst direction error / delta')
    print('    direction error = ||s - c e||/||s|| (Euclidean in the GT coordinates, scale free); phi and phi\'; the 12 branch records of the 6 sectors with m = 2 (both copies)')
    amp = {}
    for lam in W4REPS:
        for tag in ('', '2'):
            o = _load_records(lam, tag)
            phi = tuple(Fr(x) for x in o['phi'])
            for rec in o['br']:
                mu, cp = rec['mu'], rec['copy']
                m = len(prv_basis(lam, mu))
                if m != 2:
                    continue
                ex = rec['slex']
                for npts_rule, pts in (('m+1 default', default_points(phi, m + 1)), ('m+1 chosen', choose_points(lam, mu, phi, m + 1)), ('2m+2 default', default_points(phi, 2 * m + 2))):
                    P = DictionaryProblem(lam, mu, phi, pts, dps=dps)
                    vals = [_pe(rec['LB'], *p) for p in pts]
                    for delta in (mp.mpf(10) ** -10, mp.mpf(10) ** -20, mp.mpf(10) ** -40):
                        worst = mp.mpf(0)
                        for trial in range(4):
                            o1 = P.solve(_noisy(vals, delta, rng))
                            worst = max(worst, dir_err(o1['raw'], ex))
                        key = (npts_rule, int(-mp.log10(delta)))
                        amp[key] = max(amp.get(key, mp.mpf(0)), worst / delta)
    for key in sorted(amp):
        print('    m=2, points %-12s delta=1e-%d: worst amplification %s' % (key[0], key[1], mp.nstr(amp[key], 3)))
    print('=== D3(i,ii) at scale: synthetic data (random exact element of M_mu; L = c-hat <phi_u|s> exactly), one irrep per multiplicity m (the one with the largest w), delta = 1e-20 relative noise (3 trials)')
    print('    amplification = worst direction error/delta (and, in brackets, the same for the componentwise error relative to the first nonzero entry)')
    print('    %-10s %-12s %2s %2s | cond: default, chosen, default 2m+2 | amplification: m+1 default, m+1 chosen, 2m+2 default' % ('rep', 'mu', 'm', 'w'))
    for lam in [(8, 4, 0), (14, 7, 0), (28, 14, 0)]:
        R = rep1(lam)
        phi = PHI
        best = {}
        for mu, m in irreps(lam):
            w = len(weight_space(lam, tuple(mu[i] - lam[i] for i in range(3))))
            if m not in best or w > best[m][1]:
                best[m] = (mu, w)
        for m in sorted(best):
            if m == 1:
                continue
            mu, w = best[m]
            cs, sd = _rand_slice(lam, mu, rng)
            ex = _dense(sd, R.dim)
            conds, amps, ampp = [], [], []
            for pts in (default_points(phi, m + 1), choose_points(lam, mu, phi, m + 1), default_points(phi, 2 * m + 2)):
                P = DictionaryProblem(lam, mu, phi, pts, dps=dps)
                vals = [to_mp(v) for v in synthetic_values(lam, mu, phi, pts, sd)]
                conds.append(P.cond())
                delta = mp.mpf(10) ** -20
                worst = worstp = mp.mpf(0)
                for trial in range(3):
                    o2 = P.solve(_noisy(vals, delta, rng))
                    worst = max(worst, dir_err(o2['raw'], ex))
                    worstp = max(worstp, slice_dev(o2['slice'], ex))
                amps.append(worst / delta)
                ampp.append(worstp / delta)
            print('    %-10s %-12s %2d %2d | %-9s %-9s %-9s | %-9s (%-9s) %-9s (%-9s) %-9s (%s)' % (
                lam, mu, m, w, mp.nstr(conds[0], 3), mp.nstr(conds[1], 3), mp.nstr(conds[2], 3), mp.nstr(amps[0], 3), mp.nstr(ampp[0], 3), mp.nstr(amps[1], 3), mp.nstr(ampp[1], 3),
                mp.nstr(amps[2], 3), mp.nstr(ampp[2], 3)), flush=True)


# ------------------------------------------------------------------------------------------------------------- D3(iii): synthetic scale test
def stage_scale(dps=50):
    import random
    mp.mp.dps = dps
    rng = random.Random(777)
    print('=== D3(iii): synthetic scale test (SYNTHETIC data: random exact element of M_mu -> L(w_k) = c-hat <phi_u|s> -> recover), every irrep; exact path = Fractions, numeric path = mp %d digits' % dps)
    print('    columns: rep, dim, #irreps, #copies, max m, max w | exact recovery (m+1 default points) | numeric, default m+1: worst direction error, worst residual, worst cond | numeric, chosen m+1 (pool 40): worst direction error, worst cond | time/irrep (mean, max): exact, numeric')
    summary = []
    for lam in [(4, 2, 0), (6, 3, 0), (8, 4, 0), (14, 7, 0), (28, 14, 0)]:
        R = rep1(lam)
        phi = PHI
        irr = irreps(lam)
        ok = 0
        wd = wr = wc = wd2 = wc2 = mp.mpf(0)
        te, tn = [], []
        for mu, m in irr:
            cs, sd = _rand_slice(lam, mu, rng)
            ex = _dense(sd, R.dim)
            pts = default_points(phi, m + 1)
            vals = synthetic_values(lam, mu, phi, pts, sd)                   # exact Fractions
            t0 = time.time()
            o1 = slice_from_points_ex(lam, mu, phi, pts, vals, dps=dps)       # exact path (all rational)
            te.append(time.time() - t0)
            ok += int(o1['exact'] and o1['x_exact'] == cs and o1['residual'] == 0)
            valsm = [to_mp(v) for v in vals]
            t0 = time.time()
            P = DictionaryProblem(lam, mu, phi, pts, dps=dps)
            o2 = P.solve(valsm)
            tn.append(time.time() - t0)
            wd = max(wd, dir_err(o2['raw'], ex))
            wr = max(wr, o2['residual'])
            if m > 1:
                wc = max(wc, P.cond())
                pts2 = choose_points(lam, mu, phi, m + 1)
                P2 = DictionaryProblem(lam, mu, phi, pts2, dps=dps)
                o3 = P2.solve([to_mp(v) for v in synthetic_values(lam, mu, phi, pts2, sd)])
                wd2 = max(wd2, dir_err(o3['raw'], ex))
                wc2 = max(wc2, P2.cond())
        print('    %-10s %5d %4d %5d %3d %3d | %3d/%-3d | %s %s %s | %s %s | %.4f %.4f s, %.4f %.4f s' % (
            lam, R.dim, len(irr), sum(m for _, m in irr), max(m for _, m in irr), max(len(v) for v in R.wspace.values()), ok, len(irr),
            mp.nstr(wd, 3), mp.nstr(wr, 3), mp.nstr(wc, 3), mp.nstr(wd2, 3), mp.nstr(wc2, 3), sum(te) / len(te), max(te), sum(tn) / len(tn), max(tn)), flush=True)
        summary.append((lam, len(irr), ok, wd, wr, wc, wd2, wc2))
    return summary


def stage_points(dps=40):
    print('=== D3(i) which points: worst cond(D A) over ALL irreps of multiplicity m (m >= 2), for four point rules (dps %d):' % dps)
    print('    default = the first m+1 `default_points`; chosen = per-irrep greedy choice (pivoted QR, pool of 40); shared = prefix of length m+1 of ONE ordered list `shared_points(lam, phi)` for the whole Lambda; default2 = first 2m+2 default points')
    for lam in [(8, 4, 0), (14, 7, 0), (28, 14, 0)]:
        phi = PHI
        t0 = time.time()
        shared = shared_points(lam, phi)
        t_sh = time.time() - t0
        worst = {}
        cnt = {}
        for mu, m in irreps(lam):
            if m == 1:
                continue
            basis = prv_basis(lam, mu)
            cs = [DictionaryProblem(lam, mu, phi, pts, dps=dps, basis=basis).cond() for pts in (
                default_points(phi, m + 1), choose_points(lam, mu, phi, m + 1, basis=basis), shared[:m + 1], default_points(phi, 2 * m + 2))]
            w = worst.setdefault(m, [mp.mpf(0)] * 4)
            for i in range(4):
                w[i] = max(w[i], cs[i])
            cnt[m] = cnt.get(m, 0) + 1
        print('  %s (shared list of %d points built in %.1fs)' % (lam, len(shared), t_sh))
        print('    %3s %5s | %-10s %-10s %-10s %-10s' % ('m', '#irr', 'default', 'chosen', 'shared', 'default2'))
        for m in sorted(worst):
            print('    %3d %5d | %-10s %-10s %-10s %-10s' % (m, cnt[m], mp.nstr(worst[m][0], 3), mp.nstr(worst[m][1], 3), mp.nstr(worst[m][2], 3), mp.nstr(worst[m][3], 3)), flush=True)
        print('    worst over all m: default %s, chosen %s, shared %s, default2 %s' % tuple(mp.nstr(max(worst[m][i] for m in worst), 3) for i in range(4)))


def stage_misc(dps=60):
    """interface robustness: mp inputs, double-precision data, all copies at once, wrong-irrep negative control, coordinates"""
    import nu_general as N
    mp.mp.dps = dps
    print('=== misc: robustness of the interface (real branches, N1 records, phi)')
    # (a) mp inputs (points and phi as mpf) agree with exact inputs
    lam = (3, 1, 0)
    o = _load_records(lam, '')
    phi = tuple(Fr(x) for x in o['phi'])
    rec = [r for r in o['br'] if r['mu'] == (5, 2, 1) and r['copy'] == 2][0]
    pts = default_points(phi, 3)
    vals = [_pe(rec['LB'], *p) for p in pts]
    x1, r1, s1 = slice_from_points(lam, rec['mu'], phi, pts, vals, dps=dps)
    phim = tuple(to_mp(v) for v in phi)
    ptsm = [(to_mp(a), to_mp(b)) for a, b in pts]
    x2, r2, s2 = slice_from_points(lam, rec['mu'], phim, ptsm, vals, dps=dps)
    print('  (a) [3,1,0] (5,2,1)#2: exact-rational inputs vs mpf inputs (phi and points as mpf): max |slice difference| = %s, residuals %s / %s' % (
        mp.nstr(max(abs(a - b) for a, b in zip(s1, s2)), 3), mp.nstr(r1, 3), mp.nstr(r2, 3)))
    # (b) double-precision data
    print('  (b) values rounded to double precision (float64) -- what a machine-precision kernel (WP-K) would deliver; slice error vs exact (componentwise, first-nonzero normalisation) and direction error:')
    for lam in W4REPS:
        for tag in ('',):
            o = _load_records(lam, tag)
            phi = tuple(Fr(x) for x in o['phi'])
            for rec in o['br']:
                mu = rec['mu']
                m = len(prv_basis(lam, mu))
                if m != 2 or rec['copy'] != 1:
                    continue
                for rule, pts in (('m+1 default', default_points(phi, m + 1)), ('m+1 chosen', choose_points(lam, mu, phi, m + 1)), ('2m+2 default', default_points(phi, 2 * m + 2))):
                    vals = [mp.mpf(float(_pe(rec['LB'], *p))) for p in pts]
                    o1 = slice_from_points_ex(lam, mu, phi, pts, vals, dps=dps)
                    print('      %s %s %-12s: slice error %s, direction error %s, residual %s' % (lam, mu, rule, mp.nstr(slice_dev(o1['slice'], rec['slex']), 3), mp.nstr(dir_err(o1['raw'], rec['slex']), 3), mp.nstr(o1['residual'], 3)))
    # (c) all copies at once; coordinates
    print('  (c) all copies of an irrep at once (slices_from_points) and coordinates x_j = v[f_j]:')
    for lam, mu in (((3, 1, 0), (5, 2, 1)), ((3, 2, 0), (5, 4, 1)), ((2, 1, 0), (3, 2, 1))):
        for tag in ('',):
            o = _load_records(lam, tag)
            phi = tuple(Fr(x) for x in o['phi'])
            recs = [r for r in o['br'] if r['mu'] == mu]
            pts = default_points(phi, 3)
            res = slices_from_points(lam, mu, phi, pts, [[_pe(r['LB'], *p) for p in pts] for r in recs], dps=dps)
            for r, (x, resid, sl) in zip(recs, res):
                xs, err = coords_of(lam, mu, r['slex'])
                scale = x[0] / xs[0]
                print('      %s %s copy %d: residual %s, slice error %s; coordinates of the exact slice %s (membership error %s), x/(coords) ratio spread %s, free indices %s' % (
                    lam, mu, r['copy'], mp.nstr(resid, 3), mp.nstr(slice_dev(sl, r['slex']), 3), [mp.nstr(v, 6) for v in xs], mp.nstr(err, 3),
                    mp.nstr(max(abs(x[j] / xs[j] / scale - 1) for j in range(len(x))), 3), prv_free_indices(lam, mu)))
    # (d) wrong-irrep / wrong-branch negative control: the residual is a real consistency check
    print('  (d) negative control: feed the L_B of branch B into the dictionary of ANOTHER irrep with the same m, w (m+2 points); residuals right vs wrong:')
    for lam in W4REPS:
        o = _load_records(lam, '')
        phi = tuple(Fr(x) for x in o['phi'])
        recs = o['br']
        ok_res, bad_res = [], []
        for ra in recs:
            ma = len(prv_basis(lam, ra['mu']))
            pts = default_points(phi, ma + 2)
            va = [_pe(ra['LB'], *p) for p in pts]
            for rb in recs:
                if rb['copy'] != 1 and len(prv_basis(lam, rb['mu'])) == 1:
                    continue
                mb = len(prv_basis(lam, rb['mu']))
                if mb != ma:
                    continue
                if ra['mu'] == rb['mu']:
                    ok_res.append(slice_from_points_ex(lam, rb['mu'], phi, pts, va, dps=dps)['residual'])
                else:
                    try:
                        bad_res.append(slice_from_points_ex(lam, rb['mu'], phi, pts, va, dps=dps)['residual'])
                    except Exception:
                        pass
        print('      %s: right irrep: worst residual %s (%d cases); wrong irrep: smallest residual %s (%d cases)' % (
            lam, mp.nstr(max(ok_res), 3), len(ok_res), mp.nstr(min(bad_res), 3) if bad_res else None, len(bad_res)))


# ------------------------------------------------------------------------------------------------------------- cross-check against agent R (slice_engine.py / ref/cg-slices-*.txt)
REFDIR = os.path.join(HERE, '..', '..', 'Mathematica', 'Gaudin', 'ref')


def _parse_R(s):
    import re
    s = re.sub(r'(-?\d+)/(\d+)', r'F(\1,\2)', s)
    s = re.sub(r'(?<![\w.])(-?\d+\.\d+)', r"M('\1')", s)
    return eval(s.replace('{', '[').replace('}', ']'), {'F': Fr, 'M': mp.mpf})


def load_R(lam):
    """parse WP-R's ref/cg-slices-<lam>.txt: {mu: dict(m, w, wspace (0-based my indices), prv [dict nu -> Fraction], copies [(eigenvalue, dict nu -> mpf)])}"""
    fn = os.path.join(REFDIR, 'cg-slices-%d%d%d.txt' % tuple(lam))
    out = {}
    for line in open(fn, encoding='utf8'):
        line = line.rstrip('\n')
        if not line or line.startswith('#'):
            continue
        f = line.split('\t')
        tag = f[0]
        if tag == 'SECTOR':
            mu = tuple(_parse_R(f[1]))
            out[mu] = dict(m=int(f[2]), w=int(f[3]), prv=[], copies=[])
        elif tag == 'WSPACE':
            out[tuple(_parse_R(f[1]))]['wspace'] = [i - 1 for i in _parse_R(f[2])]
        elif tag == 'PRVBASIS':
            out[tuple(_parse_R(f[1]))]['prv'].append({nu - 1: v for nu, v in _parse_R(f[2])})
        elif tag == 'COPY':
            out[tuple(_parse_R(f[1]))]['copies'].append((mp.mpf(f[3]), {nu - 1: v for nu, v in _parse_R(f[4])}))
    return out


def stage_crossR(dps=40):
    """independent cross-check against WP-R: PRV bases identical, R's V_red eigen-copies lie in M_mu and are recovered by the dictionary from the identity at the shared points"""
    mp.mp.dps = dps
    print('=== cross-check against WP-R (slice_engine.py, ref/cg-slices-<lam>.txt): PRV bases (exact), COPY slices (30 digits) in M_mu, dictionary round trip on the copies of R (identity -> L at shared_points[:m+1] -> slice)')
    for lam in [(4, 2, 0), (6, 3, 0), (8, 4, 0), (14, 7, 0)]:
        fn = os.path.join(REFDIR, 'cg-slices-%d%d%d.txt' % lam)
        if not os.path.exists(fn):
            print('  %s: %s not found, skipped' % (lam, fn))
            continue
        t0 = time.time()
        Rd = load_R(lam)
        irr = dict(irreps(lam))
        phi = PHI
        sh = shared_points(lam, phi)
        same_set = set(Rd) == set(irr)
        nbasis = nbad_dim = nbad_w = nbad_prv = ncopy = 0
        worst_mem = worst_rt = worst_res = mp.mpf(0)
        for mu, rd in Rd.items():
            wsp = weight_space(lam, tuple(mu[i] - lam[i] for i in range(3)))
            basis = prv_basis(lam, mu)
            m = len(basis)
            nbad_dim += int(m != rd['m'] or m != irr.get(mu))
            nbad_w += int(len(wsp) != rd['w'] or wsp != rd['wspace'])
            nbasis += m
            nbad_prv += int(basis != rd['prv'])
            fr = prv_free_indices(lam, mu)
            P = DictionaryProblem(lam, mu, phi, sh[:m + 1], dps=dps)
            n = P.R.dim
            for ev, sl in rd['copies']:
                ncopy += 1
                x = [sl.get(f, mp.mpf(0)) for f in fr]
                rec = {}
                for xj, b in zip(x, basis):
                    for i, v in b.items():
                        rec[i] = rec.get(i, 0) + xj * to_mp(v)
                worst_mem = max(worst_mem, max([abs(sl.get(i, mp.mpf(0)) - rec.get(i, mp.mpf(0))) for i in set(sl) | set(rec)] + [mp.mpf(0)]))
                vals = synthetic_values(lam, mu, phi, sh[:m + 1], sl)
                o1 = P.solve(vals)
                ref = [sl.get(i, mp.mpf(0)) for i in range(n)]
                worst_rt = max(worst_rt, max(abs(a - b) for a, b in zip(o1['slice'], ref)))
                worst_res = max(worst_res, o1['residual'])
        print('  %-10s: %3d sectors (same irrep set as mine: %s), %3d PRV basis vectors, %3d COPY slices | dim M_mu, w, WSPACE mismatches: %d, %d | PRVBASIS != mine in %d sectors | COPY in M_mu: max |v - sum v[f_j] b_j| = %s | dictionary round trip on the copies of R: max |slice - R| = %s, worst residual %s | %.1fs' % (
            lam, len(Rd), same_set, nbasis, ncopy, nbad_dim, nbad_w, nbad_prv, mp.nstr(worst_mem, 3), mp.nstr(worst_rt, 3), mp.nstr(worst_res, 3), time.time() - t0), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stages', default='polys,prv,timing,branches,cond,points,scale,misc,crossR')
    ap.add_argument('--amax', type=int, default=14)
    a = ap.parse_args()
    for st in a.stages.split(','):
        t0 = time.time()
        if st == 'timing':
            stage_timing(a.amax)
        else:
            globals()['stage_' + st]()
        print('--- stage %s finished in %.1fs' % (st, time.time() - t0), flush=True)


if __name__ == '__main__':
    main()
