"""jacobi_cg.py -- round 7, WP-J: the Jacobi (three-term-recurrence) structure of the reduced Gaudin vertex operator V_red on the weight spaces of V_Lambda,
the eigenvector formula  slice = (p_1(V), p_2(V), ...)  and its Q-side meaning.   Brief: SPEC-jacobi-cg.md.   Report: RESULTS-jacobi-cg.md.

Maths in one paragraph.  Psi -> slice psi_top = Psi[., q_top] identifies the HW vectors of weight mu in V_Lam (x) V_Lam with M_mu = {v in W : E12^(l1-l2+1) v = E23^(l2-l3+1) v = 0},
W = V_Lam[mu-Lam] (dimension w), and the Gaudin vertex operator acts on M_mu as V_red = 1/2 sum sgn(s) sgn(t) E_{s1 t1} E_{s2 t2} Phi_{s3 t3} (Phi_cc = l_c, Phi_ab = -E_ab (a<b), 0 (a>b)),
a weight-zero operator on W.  In the GT basis {t; m1, m2} a weight space has t = mu1-l1 and s = m1+m2 = (mu1-l1)+(mu2-l2) FIXED, so m1 runs along a chain (m2 = s-m1).
Every term of V_red contains exactly one generator with first index 3 and one with second index 3 (or one E33), hence moves (m1,m2) by 0 or +-(1,-1): V_red is TRIDIAGONAL on the chain
(structural proof, `check_index_structure`).  The closed forms (derived symbolically in `derive_symbolic`, not fitted), x = m1-m2, xi = x(x+2), t = mu1-l1, s = m1+m2:
    u(m1) = M[m1+1,m1] = -(l1-m1)(l2-m1-1)(l3-m1-2)(t+1-m2)(mu1-l2-m1) / ((x+1)(x+2))
    d(m1) = M[m1-1,m1] = -(l2-m2)(t-m1)(l1-m2+1)(l3-m2-1)(mu1-l2-m2+1) / (x(x+1))
    alpha(m1) = A0(l,s,t) + 3/16 (s-2 l3+2) xi - K/(16 xi),   K = (s-2l2)(s-2t)(s-2l1-2)(s-2l3+2)(s-2mu1+2l2-2)     (K/xi := 0 at x = 0, where K = 0)
The factor (mu1-l2-m1) of u vanishes at m1 = mu1-l2, exactly where the E12 constraint cuts the chain (E12 is diagonal in the chain coordinate): the initial segment W1 = {m1 <= mu1-l2}
is V_red-invariant, V_red|W1 is a Jacobi matrix (all off-diagonal products > 0), M_mu lies in W1.  The E23 constraint is a coordinate cut in a SECOND chain: V' = U V_red U^-1 (U = the Weyl element w0,
E_ij -> E_{4-i,4-j}) on W' = V_Lam[(e3,e2,e1)], where E23^(l2-l3+1) becomes E21^(l2-l3+1), diagonal in the same GT basis: cut m1 <= mu2-l3, again Jacobi (closed forms of chain 2 differ only in the last
factor of u, d, in the pole residue and in the cubic polynomial).  At most one constraint is non-void (the image weights need mu2 >= l1+l3+1 resp. mu2 <= l1+l3-1): for mu2 >= l1+l3, M_mu = W1 (m = w1);
for mu2 <= l1+l3, U M_mu = W2 (m = w2).  So m = w1 or w2 is a closed-form multiplicity (`mult_closed`) and the physical V-spectrum (= the Q-system polynomial) is the spectrum of an explicit
m x m Jacobi matrix (`physical_V`).  The slice of a copy with Q-system eigenvalue V is the forward three-term recurrence p_1 = 1, p_{i+1} = ((V - alpha_i) p_i - u_{i-1} p_{i-1}) / d_i on W1 (zero
beyond), O(w) work (`slice_from_V`).

Public API (everything exact unless it says mp/float; lam = (l1,l2,l3), mu = (mu1,mu2,mu3) tuples of ints)
  chain(lam, mu, kind=1|2) -> Chain   closed-form chain: t, s, lo, hi, w, cut, w1, idx (0-based GL3Rep indices), alpha/u/d (Fractions); matrix(), products(), contin(n); NO representation matrices are built
  slice_from_V(lam, mu, V, dps=30)          mp slice (list of mpf aligned with Chain.idx, zeros beyond the cut, first component 1);  slice_from_V_float(ch, V, twisted=False|True)  python floats
  slices_from_V, slice_dict(ch, p) -> {0-based global GT index: value}
  physical_V(lam, mu, dps, method='auto'|'top')   Q-free physical V-spectrum (chain 1 on W1 / chain 2 on W2;  or the top m eigenvalues of chain 1 on W1);  jacobi_eigs
  rs_mult / mult_closed / irreps / wt_mult / gt_index   Racah-Speiser, closed-form (PRV) multiplicity, irreps of V_lam (x) V_lam, weight multiplicity, GL3Rep index of {t;m1,m2} -- all closed form
  pipeline(lam, dps)                        all irreps, all copies, from the chains alone (Q-free)
Verification stages (CLI: python -X utf8 jacobi_cg.py <stage> [args]):  derive | closed | tables l1 l2 [l3] | dump_tables | struct | chain2 | scan [lmax] | prv [lmax] | recur | qside | qdata | systems [l1 l2 l3] | bench | all
"""
import os, sys, time, itertools, math, re, pickle, random
from fractions import Fraction as Fr
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
import numpy as np
import mpmath as mp
import sympy as sp

REFDIR = os.path.normpath(os.path.join(HERE, '..', '..', 'Mathematica', 'Gaudin', 'ref'))
MAIN_REPS = [(4, 2, 0), (6, 3, 0), (8, 4, 0), (14, 7, 0)]


# ======================================================================================== 1. combinatorics: weights, GL3Rep index map, Racah-Speiser, irreps
def chain_bounds(lam, t, s):
    """valid m1 of the GT patterns {t; m1, s-m1; lam}:  lam2 <= m1 <= lam1, lam3 <= s-m1 <= lam2, s-m1 <= t <= m1   ->   lo <= m1 <= hi"""
    l1, l2, l3 = lam
    return max(l2, s - l2, t, s - t), min(l1, s - l3)


def wt_mult(lam, wt):
    """multiplicity of the weight wt = (E11,E22,E33) in V_lam = number of valid m1 (closed form; = number of GT patterns of that weight)"""
    if wt[0] + wt[1] + wt[2] != sum(lam):
        return 0
    lo, hi = chain_bounds(lam, wt[0], wt[0] + wt[1])
    return max(0, hi - lo + 1)


_PERMS = [((0, 1, 2), 1), ((1, 0, 2), -1), ((0, 2, 1), -1), ((1, 2, 0), 1), ((2, 0, 1), 1), ((2, 1, 0), -1)]


def rs_mult(lam, mu, lam2=None):
    """Racah-Speiser/Klimyk multiplicity of V_mu in V_lam (x) V_lam2 (default lam2 = lam): sum_s sgn(s) mult_lam[s(mu+rho) - rho - lam2]  (rho = (2,1,0)), closed form (no GT patterns built)"""
    lam2 = lam if lam2 is None else lam2
    rho = (2, 1, 0)
    mr = tuple(mu[i] + rho[i] for i in range(3))
    tot = 0
    for perm, sg in _PERMS:
        wt = tuple(mr[perm[i]] - rho[i] - lam2[i] for i in range(3))
        tot += sg * wt_mult(lam, wt)
    return tot


def irreps(lam):
    """{mu: multiplicity} of V_lam (x) V_lam for all dominant mu with multiplicity > 0, sorted descending lexicographically (closed form)"""
    lam = tuple(lam)
    l1, l3 = lam[0], lam[2]
    L = sum(lam)
    out = {}
    for e1 in range(l3, l1 + 1):
        for e2 in range(l3, l1 + 1):
            e3 = L - e1 - e2
            if not (l3 <= e3 <= l1):
                continue
            mu = (lam[0] + e1, lam[1] + e2, lam[2] + e3)
            if mu[0] >= mu[1] >= mu[2]:
                m = rs_mult(lam, mu)
                if m > 0:
                    out[mu] = m
    return dict(sorted(out.items(), reverse=True))


_OFF = {}


def _offsets(lam):
    """offsets of the (m1,m2) blocks in the GL3Rep order (ascending lexicographic in (m1,m2,t)); returns (dict, dim)"""
    lam = tuple(lam)
    if lam not in _OFF:
        l1, l2, l3 = lam
        off, c = {}, 0
        for m1 in range(l2, l1 + 1):
            for m2 in range(l3, l2 + 1):
                off[(m1, m2)] = c
                c += m1 - m2 + 1
        _OFF[lam] = (off, c)
    return _OFF[lam]


def gt_index(lam, t, m1, m2):
    """0-based index of the GT pattern {t; m1,m2; lam} in the order of gl3_reps.GL3Rep(lam).patterns"""
    return _offsets(lam)[0][(m1, m2)] + (t - m2)


def dim_gl3(lam):
    return (lam[0] - lam[1] + 1) * (lam[1] - lam[2] + 1) * (lam[0] - lam[2] + 2) // 2


# ======================================================================================== 2. the closed forms (exact)
def up_coef(lam, t, m1, m2):
    """u = coefficient of e(m1+1, m2-1) in V_red e(m1, m2) at fixed t:  -(l1-m1)(l2-m1-1)(l3-m1-2)(t+1-m2)(l1-l2-m1+t) / ((x+1)(x+2)),  x = m1-m2"""
    l1, l2, l3 = lam
    x = m1 - m2
    return -Fr((l1 - m1) * (l2 - m1 - 1) * (l3 - m1 - 2) * (t + 1 - m2) * (l1 - l2 - m1 + t), (x + 1) * (x + 2))


def dn_coef(lam, t, m1, m2):
    """d = coefficient of e(m1-1, m2+1) in V_red e(m1, m2):  -(l2-m2)(t-m1)(l1-m2+1)(l3-m2-1)(l1-l2-m2+t+1) / (x(x+1)),  x = m1-m2 >= 1 (x = 0 never occurs for a valid down step)"""
    l1, l2, l3 = lam
    x = m1 - m2
    return -Fr((l2 - m2) * (t - m1) * (l1 - m2 + 1) * (l3 - m2 - 1) * (l1 - l2 - m2 + t + 1), x * (x + 1))


# 8 A0, the x-independent part of the diagonal entry: a cubic polynomial in (l1,l2,l3,s,t) (text shared by the exact Python function and the sympy comparison in derive_symbolic)
_A0_POLY = ("4*l1**2*l2 + 4*l1**2*l3 + 2*l1**2*s - 4*l1**2*t - 4*l1**2 + 4*l1*l2**2 - 4*l1*l2*l3 + 2*l1*l2*s + 12*l1*l2 + 4*l1*l3*s - 4*l1*l3*t + 8*l1*l3 - 2*l1*s**2 "
            "+ 6*l1*s*t - 4*l1*t**2 - 4*l1*t - 8*l1 + 4*l2**2*l3 - 2*l2**2*s + 4*l2**2*t + 4*l2*l3*t - 4*l2*l3 + 2*l2*s*t + 2*l2*s - 4*l2*t**2 - 4*l2*t + 8*l2 "
            "+ 2*l3*s**2 + 4*l3*s*t + 4*l3*s - 4*l3*t**2 - 4*l3*t - s**3 - 2*s**2*t - 4*s**2 + 2*s*t**2 + 2*s*t")
_A0_CODE = compile(_A0_POLY, '<A0>', 'eval')


def A0(lam, s, t):
    """A0(l; s, t) = (_A0_POLY)/8, exact Fraction"""
    l1, l2, l3 = lam
    return Fr(eval(_A0_CODE, {'__builtins__': {}}, {'l1': l1, 'l2': l2, 'l3': l3, 's': s, 't': t}), 8)


def Kpole(lam, s, t):
    """K = (s-2l2)(s-2t)(s-2l1-2)(s-2l3+2)(s-2l1+2l2-2t-2) (the pole residue is K/32; the pole term of the diagonal is -K/(16 xi))"""
    l1, l2, l3 = lam
    return (s - 2 * l2) * (s - 2 * t) * (s - 2 * l1 - 2) * (s - 2 * l3 + 2) * (s - 2 * l1 + 2 * l2 - 2 * t - 2)


def diag_coef(lam, t, m1, m2):
    """alpha = A0 + 3/16 (s - 2 l3 + 2) xi - K/(16 xi),  xi = x(x+2) = (m1-m2)(m1-m2+2) = 4 j(j+1) (j = (m1-m2)/2 the gl2 spin); the pole term is 0 at x = 0 (there K = 0: m1 = m2 = l2)"""
    l1, l2, l3 = lam
    x = m1 - m2
    s = m1 + m2
    xi = x * (x + 2)
    val = A0(lam, s, t) + Fr(3 * (s - 2 * l3 + 2) * xi, 16)
    if xi != 0:
        val -= Fr(Kpole(lam, s, t), 16 * xi)
    else:
        assert Kpole(lam, s, t) == 0
    return val


# ---- SECOND CHAIN: V'_red = U V_red U^{-1} on the (12)-adapted GT basis, U = the Weyl element w0 (E_ij -> E_{4-i,4-j}; it maps the (12)-adapted GT basis to the (23)-adapted one).  V' has the same form with
#      Phi'_cc = l_{4-c}, Phi'_ab = -E_ab for a > b (LOWERING), 0 for a < b; it acts on W' = V_lam[(e3,e2,e1)], (e1,e2,e3) = mu - lam (so t' = mu3 - l3, s' = (mu3-l3)+(mu2-l2)), and U M_mu =
#      {v' in W' : E32^(l1-l2+1) v' = E21^(l2-l3+1) v' = 0}.  E21 is DIAGONAL in the chain (t -> t-1): E21^(l2-l3+1) keeps the components with m2 >= t' - (l2-l3), i.e. m1 <= mu2 - l3 (the second cut).
#      The closed forms differ from chain 1 only in the last factor of u, d, in the pole residue and in the cubic polynomial (derived symbolically, `derive_symbolic(second=True)`).
_A0P_POLY = ("4*l1*l2**2 - 4*l1*l2*l3 + 4*l1*l2*t + 4*l1*l2 + 4*l1*l3**2 + 4*l1*l3*s - 4*l1*l3*t - 8*l1*l3 + 2*l1*s**2 + 4*l1*s*t - 4*l1*s - 4*l1*t**2 + 4*l1*t + 4*l2**2*l3 - 2*l2**2*s "
             "+ 4*l2**2*t + 4*l2*l3**2 + 2*l2*l3*s - 12*l2*l3 + 2*l2*s*t - 2*l2*s - 4*l2*t**2 + 4*l2*t + 8*l2 + 2*l3**2*s - 4*l3**2*t + 4*l3**2 - 2*l3*s**2 + 6*l3*s*t - 4*l3*t**2 "
             "+ 4*l3*t - 8*l3 - s**3 - 2*s**2*t + 4*s**2 + 2*s*t**2 - 2*s*t")
_A0P_CODE = compile(_A0P_POLY, '<A0p>', 'eval')


def A0p(lam, s, t):
    l1, l2, l3 = lam
    return Fr(eval(_A0P_CODE, {'__builtins__': {}}, {'l1': l1, 'l2': l2, 'l3': l3, 's': s, 't': t}), 8)


def up2_coef(lam, t, m1, m2):
    l1, l2, l3 = lam
    x = m1 - m2
    return -Fr((l1 - m1) * (l2 - m1 - 1) * (l3 - m1 - 2) * (t + 1 - m2) * (t + l3 - l2 - m2), (x + 1) * (x + 2))


def dn2_coef(lam, t, m1, m2):
    l1, l2, l3 = lam
    x = m1 - m2
    return -Fr((l2 - m2) * (t - m1) * (l1 - m2 + 1) * (l3 - m2 - 1) * (t + l3 - l2 - m1 - 1), x * (x + 1))


def Kpole2(lam, s, t):
    l1, l2, l3 = lam
    return (s - 2 * l2) * (s - 2 * t) * (s - 2 * l1 - 2) * (s - 2 * l3 + 2) * (s + 2 * l2 - 2 * l3 - 2 * t + 2)


def diag2_coef(lam, t, m1, m2):
    l1, l2, l3 = lam
    x = m1 - m2
    s = m1 + m2
    xi = x * (x + 2)
    val = A0p(lam, s, t) + Fr(3 * (s - 2 * l1 - 2) * xi, 16)
    if xi != 0:
        val -= Fr(Kpole2(lam, s, t), 16 * xi)
    else:
        assert Kpole2(lam, s, t) == 0
    return val


class Chain:
    """closed-form chain of V_red on W = V_lam[mu-lam]:  positions k = 0..w-1 <-> m1 = lo+k (m2 = s-m1, t fixed).  M[i,j] = coefficient of e_i in V_red e_j:
       M[k,k] = alpha[k], M[k+1,k] = u[k] (up step m1 -> m1+1), M[k,k+1] = d[k] (down step m1+1 -> m1), all other entries 0 (tridiagonal).
       cut = mu1 - l2 (E12 constraint: components with m1 > cut vanish on M_mu), w1 = #{k : lo+k <= cut}, idx = 0-based GL3Rep indices of the chain states (ascending)."""

    def __init__(self, lam, mu, kind=1):
        lam = tuple(int(x) for x in lam)
        mu = tuple(int(x) for x in mu)
        assert sum(mu) == 2 * sum(lam), (lam, mu)
        l1, l2, l3 = lam
        self.lam, self.mu, self.kind = lam, mu, kind
        if kind == 1:                                   # the chain of V_red on W = V[mu - lam], cut m1 <= mu1 - l2 (E12 constraint)
            self.t = t = mu[0] - l1
            self.s = s = (mu[0] - l1) + (mu[1] - l2)
            self.cut = mu[0] - l2
            fu, fd, fa = up_coef, dn_coef, diag_coef
        else:                                           # the second chain of V' = U V_red U^-1 on W' = V[(e3,e2,e1)], cut m1 <= mu2 - l3 (E23 constraint seen through E21)
            self.t = t = mu[2] - l3
            self.s = s = (mu[2] - l3) + (mu[1] - l2)
            self.cut = mu[1] - l3
            fu, fd, fa = up2_coef, dn2_coef, diag2_coef
        self.lo, self.hi = lo, hi = chain_bounds(lam, t, s)
        self.w = w = max(0, hi - lo + 1)
        self.w1 = max(0, min(hi, self.cut) - lo + 1) if w else 0
        self.m1s = list(range(lo, lo + w))
        self.idx = [gt_index(lam, t, m1, s - m1) for m1 in self.m1s]
        self.alpha = [fa(lam, t, m1, s - m1) for m1 in self.m1s]
        self.u = [fu(lam, t, m1, s - m1) for m1 in self.m1s[:-1]]
        self.d = [fd(lam, t, m1 + 1, s - m1 - 1) for m1 in self.m1s[:-1]]
        self._mp = {}
        self._fl = None

    # ---- exact objects
    def matrix(self, n=None):
        n = self.w if n is None else n
        M = [[Fr(0)] * n for _ in range(n)]
        for k in range(n):
            M[k][k] = self.alpha[k]
            if k + 1 < n:
                M[k + 1][k] = self.u[k]
                M[k][k + 1] = self.d[k]
        return M

    def products(self, n=None):
        """off-diagonal products M[k+1,k] M[k,k+1] = u[k] d[k] (invariant under diagonal similarity), k = 0..n-2"""
        n = self.w if n is None else n
        return [self.u[k] * self.d[k] for k in range(n - 1)]

    def contin(self, n=None):
        """characteristic polynomial det(x - M_n) of the leading n x n block (n = w1 by default) as ascending Fraction coefficients, by the continuant recursion
           D_0 = 1, D_1 = x - a_0, D_{k+1} = (x - a_k) D_k - u_{k-1} d_{k-1} D_{k-1}"""
        n = self.w1 if n is None else n
        Dm, Dk = [Fr(1)], [-self.alpha[0], Fr(1)] if n > 0 else [Fr(1)]
        if n == 0:
            return Dm
        for k in range(1, n):
            pr = self.u[k - 1] * self.d[k - 1]
            new = [Fr(0)] * (len(Dk) + 1)
            for i, c in enumerate(Dk):
                new[i + 1] += c
                new[i] -= self.alpha[k] * c
            for i, c in enumerate(Dm):
                new[i] -= pr * c
            Dm, Dk = Dk, new
        return Dk

    # ---- numeric coefficient caches
    def mp_coeffs(self, dps):
        if dps not in self._mp:
            with mp.workdps(dps):
                f = lambda x: mp.mpf(x.numerator) / x.denominator
                self._mp[dps] = ([f(x) for x in self.alpha], [f(x) for x in self.u], [f(x) for x in self.d])
        return self._mp[dps]

    def f_coeffs(self):
        if self._fl is None:
            self._fl = ([float(x) for x in self.alpha], [float(x) for x in self.u], [float(x) for x in self.d])
        return self._fl


_CH = {}


def chain(lam, mu, kind=1):
    key = (tuple(lam), tuple(mu), kind)
    if key not in _CH:
        _CH[key] = Chain(lam, mu, kind)
    return _CH[key]


def mult_closed(lam, mu):
    """CLOSED-FORM multiplicity of V_mu in V_lam (x) V_lam from the PRV theorem and the two chains (no Racah-Speiser sum): the E23 constraint is void for mu2 >= l1+l3 and then M_mu = W1, so
       m = w1 = #{m1 : max(l2, s-l2, t, s-t) <= m1 <= min(l1, s-l3, mu1-l2)}  (t = mu1-l1, s = (mu1-l1)+(mu2-l2));  the E12 constraint is void for mu2 <= l1+l3 and then U M_mu = W2,
       m = w2 = #{m1 : max(l2, s'-l2, t', s'-t') <= m1 <= min(l1, s'-l3, mu2-l3)}  (t' = mu3-l3, s' = (mu3-l3)+(mu2-l2))."""
    if mu[1] >= lam[0] + lam[2]:
        return chain(lam, mu, 1).w1
    return chain(lam, mu, 2).w1


# ======================================================================================== 3. the eigenvector formula
def slice_from_V(lam, mu, V, dps=30, ch=None, residual=False):
    """THE FORMULA.  Slice of the copy of V_mu in V_lam (x) V_lam with vertex-operator eigenvalue V (= the Q-system's pencil root), from the forward three-term recurrence on W1
         p_0 = 1,  p_1 = (V - alpha_0) / d_0,  p_{k+1} = ((V - alpha_k) p_k - u_{k-1} p_{k-1}) / d_k      (k = 1..w1-2),
       returned as a list of mpf aligned with Chain.idx (length w; zero beyond the E12 cut; first component exactly 1 = R's COPY normalisation).  residual=True also returns the residual
       of the last row of the truncated chain, (alpha_{w1-1} - V) p_{w1-1} + u_{w1-2} p_{w1-2} (0 iff V is an eigenvalue of V_red|W1)."""
    ch = ch or chain(lam, mu)
    assert ch.kind == 1, 'slice_from_V needs chain 1 (the second chain lives in the w0-image of the weight space)'
    n = ch.w1
    with mp.workdps(dps):
        V = mp.mpf(V.numerator) / V.denominator if isinstance(V, Fr) else mp.mpf(V)
        al, uu, dd = ch.mp_coeffs(dps)
        p = [mp.mpf(0)] * ch.w
        res = mp.mpf(0)
        if n > 0:
            p[0] = mp.mpf(1)
            if n > 1:
                p[1] = (V - al[0]) / dd[0]
                for k in range(1, n - 1):
                    p[k + 1] = ((V - al[k]) * p[k] - uu[k - 1] * p[k - 1]) / dd[k]
            if residual:
                res = (al[n - 1] - V) * p[n - 1] + (uu[n - 2] * p[n - 2] if n > 1 else 0)
    return (p, res) if residual else p


def slice_from_V_float(ch, V, twisted=False):
    """python-float version of the same recurrence (speed).  twisted=False: the plain forward recurrence (the formula; accuracy degrades with the chain length because the rounding error of V is
       amplified, see the RESULTS).  twisted=True: forward recurrence from the bottom (f) and backward recurrence from the top (g) of W1, glued at the twist index r with the smallest
       relative residual |rho_r|/|f_r| of the row r (Parlett-Dhillon): same O(w) cost, accuracy ~ eps/gap instead of the exponential loss."""
    al, uu, dd = ch.f_coeffs()
    n = ch.w1
    p = [0.0] * ch.w
    if n == 0:
        return p
    p[0] = 1.0
    if n == 1:
        return p
    f = [0.0] * n
    f[0] = 1.0
    f[1] = (V - al[0]) / dd[0]
    for k in range(1, n - 1):
        f[k + 1] = ((V - al[k]) * f[k] - uu[k - 1] * f[k - 1]) / dd[k]
    if not twisted:
        p[:n] = f
        return p
    g = [0.0] * n
    g[n - 1] = 1.0
    g[n - 2] = (V - al[n - 1]) / uu[n - 2]
    for k in range(n - 2, 0, -1):
        g[k - 1] = ((V - al[k]) * g[k] - dd[k] * g[k + 1]) / uu[k - 1]
    best, r = None, 0
    for k in range(n):
        gl = g[k + 1] * f[k] / g[k] if k + 1 < n else 0.0                      # glued upper neighbour
        fl = f[k - 1] if k > 0 else 0.0
        rho = (uu[k - 1] * fl if k > 0 else 0.0) + (al[k] - V) * f[k] + (dd[k] * gl if k + 1 < n else 0.0)
        crit = abs(rho) / abs(f[k]) if f[k] != 0 else float('inf')
        if best is None or crit < best:
            best, r = crit, k
    scale = f[r] / g[r] if g[r] != 0 else 0.0
    for k in range(n):
        p[k] = f[k] if k <= r else scale * g[k]
    return p


def slices_from_V(lam, mu, Vs, dps=30):
    ch = chain(lam, mu)
    return [slice_from_V(lam, mu, V, dps, ch) for V in Vs]


def slice_dict(ch, p, tol=None, one_based=False):
    """{global GT index: value} of a slice (entries with |value| <= tol dropped; the zero tail beyond the cut is dropped)"""
    out = {}
    for k, x in enumerate(p):
        if x == 0 or (tol is not None and abs(x) <= tol):
            continue
        out[ch.idx[k] + (1 if one_based else 0)] = x
    return out


def jacobi_eigs(lam, mu, dps=30, ch=None):
    """ALL w1 eigenvalues (ascending, mpf) of V_red|W1: the truncated chain is a Jacobi matrix (u_k d_k > 0), symmetrised by a positive diagonal similarity (off-diagonals sqrt(u d)), mpmath eigsy"""
    ch = ch or chain(lam, mu)
    n = ch.w1
    if n == 0:
        return []
    with mp.workdps(dps + 10):
        al, uu, dd = ch.mp_coeffs(dps + 10)
        if n == 1:
            ev = [al[0]]
        else:
            A = mp.zeros(n, n)
            for k in range(n):
                A[k, k] = al[k]
                if k + 1 < n:
                    b = mp.sqrt(uu[k] * dd[k])
                    A[k, k + 1] = b
                    A[k + 1, k] = b
            ev = sorted(mp.eigsy(A, eigvals_only=True))
    with mp.workdps(dps):
        return [+x for x in ev]


def physical_V(lam, mu, dps=30, ch=None, m=None, method='auto'):
    """Q-FREE physical V-spectrum of the sector (ascending list of mpf), two equivalent rules:
         method = 'auto' (rigorous, from the PRV theorem): if w1 = m (always when mu2 >= l1+l3) all eigenvalues of the Jacobi matrix V_red|W1 of chain 1; otherwise (w1 > m, then mu2 < l1+l3 and the E12
                  constraint is void) U M_mu = W2 is a coordinate segment of chain 2 (w2 = m): all eigenvalues of the Jacobi matrix V'_red|W2;
         method = 'top'  (empirical rule, checked on thousands of sectors): the TOP m eigenvalues of V_red|W1 of chain 1."""
    ch = ch or chain(lam, mu)
    m = rs_mult(lam, mu) if m is None else m
    if m == 0:
        return []
    if method == 'top' or ch.w1 == m:
        ev = jacobi_eigs(lam, mu, dps, ch)
        return ev[len(ev) - m:]
    ch2 = chain(lam, mu, 2)
    assert ch2.w1 == m, ('second chain', lam, mu, ch2.w1, m)
    return jacobi_eigs(lam, mu, dps, ch2)


def refine_V(lam, mu, V, dps=50, ch=None, iters=4):
    """Newton polish of an approximate V-root (e.g. Q0's 30-digit pencil root) on the exact continuant of the Jacobi matrix V_red|W1 (D_{k+1} = (V - alpha_k) D_k - u_{k-1} d_{k-1} D_{k-1}, with its derivative):
       every physical eigenvalue is a simple root of D_{w1}, so 2-3 steps take a 30-digit value to `dps` digits.  The Q-system only has to tell WHICH root; the digits come from the chain."""
    ch = ch or chain(lam, mu)
    n = ch.w1
    with mp.workdps(dps + 15):
        al, uu, dd = ch.mp_coeffs(dps + 15)
        x = mp.mpf(V.numerator) / V.denominator if isinstance(V, Fr) else mp.mpf(V)
        pr = [uu[k] * dd[k] for k in range(max(n - 1, 0))]
        for _ in range(iters):
            d0, d1 = mp.mpf(1), x - al[0]
            e0, e1 = mp.mpf(0), mp.mpf(1)
            for k in range(1, n):
                d0, d1, e0, e1 = d1, (x - al[k]) * d1 - pr[k - 1] * d0, e1, d1 + (x - al[k]) * e1 - pr[k - 1] * e0
            x = x - d1 / e1
    with mp.workdps(dps):
        return +x


def sector_copies(lam, mu, dps=30, V=None, one_based=False, tol=None):
    """all copies of the sector mu as R-style records: [(V, {GL3Rep index: value})], first nonzero component 1, ascending V.  V = list of Q-system eigenvalues (any number of digits; they are polished on the
       chain with `refine_V` to dps digits) or None (Q-free: physical_V)."""
    ch = chain(lam, mu)
    m = rs_mult(lam, mu)
    if V is None:
        Vs = physical_V(lam, mu, dps, ch, m)
    else:
        Vs = [refine_V(lam, mu, v, dps, ch) for v in sorted(V)]
    return [(v, slice_dict(ch, slice_from_V(lam, mu, v, dps, ch), tol, one_based)) for v in Vs]


def pipeline(lam, dps=30, vsource=None):
    """all irreps and all copies of V_lam (x) V_lam from the chain alone: {mu: dict(m, w, w1, V=[...], slices=[...])}.  vsource(mu) -> list of V (Q-system numbers) replaces physical_V."""
    lam = tuple(lam)
    out = {}
    for mu, m in irreps(lam).items():
        ch = chain(lam, mu)
        Vs = vsource(mu) if vsource is not None else physical_V(lam, mu, dps, ch, m)
        out[mu] = dict(m=m, w=ch.w, w1=ch.w1, V=Vs, slices=[slice_from_V(lam, mu, V, dps, ch) for V in Vs])
    return out


# ======================================================================================== 4. symbolic derivation of the closed forms (generic GT labels) and the structural argument
def derive_symbolic(verbose=True, npoints=200, seed=7, second=False):
    """Re-derive u, d, alpha from scratch with sympy, for GENERIC labels (t, m1, m2, l1, l2, l3 symbols): the GT formulas of gl3_reps (E12: t->t+1, coefficient -(t-m1)(t-m2+1); E21: t->t-1, coefficient 1;
       E23/E32: m1 or m2 +-1 with the rational coefficients of gl3_reps._raise(2,.)/_lower(2,.); E13 = [E12,E23], E31 = [E32,E21]; E11 = t, E22 = m1+m2-t, E33 = |l|-m1-m2) are applied to the state
       {t; m1, m2} through the 36-term V_red, ALL intermediate states kept (no boundary treatment); the three weight-zero components (m1 shifted by +1, 0, -1) are factored and compared with the
       closed forms (exact sympy equality for u, d; for alpha: exact equality of the rational functions).  Also compares the exact-Fraction Python functions with the symbolic result at `npoints`
       random generic integer points (guards against transcription errors between the two implementations).  Returns dict of checks.
       second=True: the same for the second chain, V' = U V_red U^-1 with Phi'_cc = l_{4-c}, Phi'_ab = -E_ab (a > b), Phi'_ab = 0 (a < b)."""
    t, m1, m2, l1, l2, l3 = sp.symbols('t m1 m2 l1 l2 l3')
    L = l1 + l2 + l3
    lamv = (l1, l2, l3)

    def add(out, key, val):
        out[key] = out.get(key, 0) + val

    def E12(vec):
        out = {}
        for (a, b, c), co in vec.items():
            T, M1, M2 = t + a, m1 + b, m2 + c
            add(out, (a + 1, b, c), co * (-(T - M1) * (T - M2 + 1)))
        return out

    def E21(vec):
        out = {}
        for (a, b, c), co in vec.items():
            add(out, (a - 1, b, c), co)
        return out

    def E23(vec):
        out = {}
        for (a, b, c), co in vec.items():
            M1, M2 = m1 + b, m2 + c
            add(out, (a, b + 1, c), co * (-(M1 - l1) * (M1 - l2 + 1) * (M1 - l3 + 2) / (M1 - M2 + 1)))
            add(out, (a, b, c + 1), co * (-(M2 - 1 - l1) * (M2 - l2) * (M2 + 1 - l3) / (M2 - M1 - 1)))
        return out

    def E32(vec):
        out = {}
        for (a, b, c), co in vec.items():
            T, M1, M2 = t + a, m1 + b, m2 + c
            add(out, (a, b - 1, c), co * ((M1 - T) / (M1 - M2 + 1)))
            add(out, (a, b, c - 1), co * ((M2 - T - 1) / (M2 - M1 - 1)))
        return out

    def Ed(i, vec):
        out = {}
        for (a, b, c), co in vec.items():
            T, M1, M2 = t + a, m1 + b, m2 + c
            out[(a, b, c)] = co * {1: T, 2: M1 + M2 - T, 3: L - M1 - M2}[i]
        return out

    def lin(v1, v2, sgn=1):
        out = dict(v1)
        for k, x in v2.items():
            add(out, k, sgn * x)
        return out

    def Eij(i, j, vec):
        if i == j:
            return Ed(i, vec)
        return {(1, 2): E12, (2, 1): E21, (2, 3): E23, (3, 2): E32,
                (1, 3): lambda v: lin(E12(E23(v)), E23(E12(v)), -1), (3, 1): lambda v: lin(E32(E21(v)), E21(E32(v)), -1)}[(i, j)](vec)

    def Phi(a, b, vec):
        if second:
            if a < b:
                return {}
            if a == b:
                return {k: x * lamv[3 - a] for k, x in vec.items()}
            return {k: -x for k, x in Eij(a, b, vec).items()}
        if a > b:
            return {}
        if a == b:
            return {k: x * lamv[a - 1] for k, x in vec.items()}
        return {k: -x for k, x in Eij(a, b, vec).items()}

    t0 = time.perf_counter()
    out = {}
    v0 = {(0, 0, 0): sp.Integer(1)}
    for s_ in itertools.permutations((1, 2, 3)):
        for t_ in itertools.permutations((1, 2, 3)):
            coef = sp.Rational(perm_sign(s_) * perm_sign(t_), 2)
            ph = Phi(s_[2], t_[2], v0)
            if not ph:
                continue
            r = Eij(s_[0], t_[0], Eij(s_[1], t_[1], ph))
            for k, x in r.items():
                add(out, k, coef * x)
    comps = {k: sp.cancel(sp.together(v)) for k, v in out.items()}
    # closed forms (sympy versions)
    x, s = m1 - m2, m1 + m2
    if not second:
        up_s = -(l1 - m1) * (l2 - m1 - 1) * (l3 - m1 - 2) * (t + 1 - m2) * (l1 - l2 - m1 + t) / ((x + 1) * (x + 2))
        dn_s = -(l2 - m2) * (t - m1) * (l1 - m2 + 1) * (l3 - m2 - 1) * (l1 - l2 - m2 + t + 1) / (x * (x + 1))
        A0_s = sp.sympify(_A0_POLY, locals=dict(l1=l1, l2=l2, l3=l3, s=s, t=t)) / 8
        K_s = (s - 2 * l2) * (s - 2 * t) * (s - 2 * l1 - 2) * (s - 2 * l3 + 2) * (s - 2 * l1 + 2 * l2 - 2 * t - 2)
        al_s = A0_s + sp.Rational(3, 16) * (s - 2 * l3 + 2) * x * (x + 2) - K_s / (16 * x * (x + 2))
    else:
        up_s = -(l1 - m1) * (l2 - m1 - 1) * (l3 - m1 - 2) * (t + 1 - m2) * (t + l3 - l2 - m2) / ((x + 1) * (x + 2))
        dn_s = -(l2 - m2) * (t - m1) * (l1 - m2 + 1) * (l3 - m2 - 1) * (t + l3 - l2 - m1 - 1) / (x * (x + 1))
        A0_s = sp.sympify(_A0P_POLY, locals=dict(l1=l1, l2=l2, l3=l3, s=s, t=t)) / 8
        K_s = (s - 2 * l2) * (s - 2 * t) * (s - 2 * l1 - 2) * (s - 2 * l3 + 2) * (s + 2 * l2 - 2 * l3 - 2 * t + 2)
        al_s = A0_s + sp.Rational(3, 16) * (s - 2 * l1 - 2) * x * (x + 2) - K_s / (16 * x * (x + 2))
    res = dict(keys=sorted(comps.keys()))
    res['up'] = sp.cancel(sp.together(comps[(0, 1, -1)] - up_s)) == 0
    res['dn'] = sp.cancel(sp.together(comps[(0, -1, 1)] - dn_s)) == 0
    res['diag'] = sp.cancel(sp.together(comps[(0, 0, 0)] - al_s)) == 0
    res['only_three_components'] = sorted(comps.keys()) == [(0, -1, 1), (0, 0, 0), (0, 1, -1)]
    # python functions vs the symbolic result at random generic integer points
    f_up = sp.lambdify((t, m1, m2, l1, l2, l3), comps[(0, 1, -1)], 'sympy')
    f_dn = sp.lambdify((t, m1, m2, l1, l2, l3), comps[(0, -1, 1)], 'sympy')
    f_dg = sp.lambdify((t, m1, m2, l1, l2, l3), comps[(0, 0, 0)], 'sympy')
    rng = random.Random(seed)
    bad = 0
    for _ in range(npoints):
        while True:
            vals = [rng.randint(-9, 25) for _ in range(6)]
            tt, a1, a2, b1, b2, b3 = vals
            if abs(a1 - a2) < 4:
                continue                                 # keep away from the poles of the intermediate-state denominators
            break
        lam_ = (b1, b2, b3)
        for fs, fp in ((f_up, up2_coef if second else up_coef), (f_dn, dn2_coef if second else dn_coef), (f_dg, diag2_coef if second else diag_coef)):
            sv = Fr(str(fs(*[sp.Integer(v) for v in (tt, a1, a2, b1, b2, b3)])))     # sympy Integers: exact rational arithmetic (python ints would give floats)
            if sv != fp(lam_, tt, a1, a2):
                bad += 1
    res['random_point_mismatches'] = bad
    res['time'] = time.perf_counter() - t0
    if verbose:
        print('  symbolic derivation of chain %d (generic labels, %d random-point comparisons): components %s, only the three weight-zero components: %s' % (2 if second else 1, npoints, res['keys'], res['only_three_components']))
        print('  u = closed form: %s,  d = closed form: %s,  alpha = A0 + 3/16 (s-2l%d%s) xi - K/(16 xi): %s,  python-vs-symbolic random-point mismatches: %d  (%.1fs)' %
              (res['up'], res['dn'], 1 if second else 3, '-2' if second else '+2', res['diag'], bad, res['time']))
    return res


def perm_sign(p):
    s = 1
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


def check_index_structure(verbose=True, second=False):
    """STRUCTURAL PROOF of tridiagonality.  V_red = 1/2 sum sgn(s) sgn(t) E_{s1 t1} E_{s2 t2} Phi_{s3 t3}; the three factors carry the index pairs (s1,t1),(s2,t2),(s3,t3), s and t permutations.
       Phi_{s3 t3} = 0 unless s3 <= t3 (then lambda_c or -E_ab, a < b).  Claim: in every surviving term exactly one factor has first index 3 and exactly one has second index 3, so the term is either
       (i) E33 (diagonal) x two generators of the gl2^(12) (they change t only) or (ii) E_{a3} (raises m1+m2 by 1: E23 or E13 = [E12,E23]) x E_{3b} (lowers it by 1: E32 or E31) x one gl2^(12) generator.
       Each of E23, E32 moves exactly one of m1, m2 by +-1 and E12, E21, E11, E22 do not touch (m1,m2), so the net shift of (m1,m2) is 0 or +-(1,-1): tridiagonal on the chain (m1 = running coordinate).
       Returns counts of the term types."""
    n_terms = n_nonzero = 0
    types = Counter()
    for s_ in itertools.permutations((1, 2, 3)):
        for t_ in itertools.permutations((1, 2, 3)):
            n_terms += 1
            if (s_[2] < t_[2]) if second else (s_[2] > t_[2]):
                continue
            n_nonzero += 1
            facs = [(s_[0], t_[0]), (s_[1], t_[1]), (s_[2], t_[2])]
            n1 = sum(1 for a, b in facs if a == 3)
            n2 = sum(1 for a, b in facs if b == 3)
            assert n1 == 1 and n2 == 1
            if (3, 3) in facs:
                types['E33 x gl2^(12) x gl2^(12) (shift 0)'] += 1
            else:
                assert any(a == 3 for a, b in facs) and any(b == 3 for a, b in facs)
                types['E_a3 x E_3b x gl2^(12) (net shift 0 or +-(1,-1))'] += 1
    if verbose:
        print("  %s: %d of 36 terms have Phi != 0; in each of them exactly one factor has first index 3 and one has second index 3: %s" % ("chain 2 (Phi' lowering)" if second else 'chain 1 (Phi raising)', n_nonzero, dict(types)))
    return dict(n_terms=n_terms, n_nonzero=n_nonzero, types=dict(types))


# ======================================================================================== 5. data loaders (R's slice engine, R's/Q0's files, S's pencil)
def _se():
    import slice_engine as se
    return se


def rname(lam):
    return '%d%d%d' % tuple(lam) if max(lam) < 10 else None


def load_R_file(lam):
    """R's reference file parsed (30-digit data): {mu: dict(m, w, charpoly (exact), copies [(ev, {1-based idx: mpf})], ...)} (slice_engine.parse_ref); files exist for 420, 630, 840, 1470"""
    se = _se()
    fn = os.path.join(REFDIR, 'cg-slices-%d%d0.txt' % (lam[0], lam[1]))
    with mp.workdps(50):
        return se.parse_ref(fn, dps=50)['sectors']


def _cache_name(prefix, lam, extra=''):
    return os.path.join(HERE, 'jacobi_cg_%s_%d_%d_%d%s.pkl' % ((prefix,) + tuple(lam) + (extra,)))


def load_R_mem(lam, dps=50, regenerate=False):
    """R's engine run in memory (exact PRV basis + exact V_red + mpmath eigen-decomposition at `dps` digits): {mu: dict(m, w, ev [mpf], sl [{0-based idx: mpf}])}.  Cached (gitignored pickle,
       numbers stored as strings); this is the 'slice_engine.py data' for reps without a ref file ([20,10,0], [28,14,0])."""
    lam = tuple(lam)
    fn = _cache_name('Rdata', lam, '_%d' % dps)
    if os.path.exists(fn) and not regenerate:
        raw = pickle.load(open(fn, 'rb'))
    else:
        se = _se()
        raw = {}
        for mu in sorted(se.racah_speiser(lam), reverse=True):
            p = se.prv_space(lam, mu)
            Vm = se.V_red(lam, mu, p)
            cps = se.copies(p, Vm, dps)
            raw[mu] = dict(m=p.m, w=p.w, ev=[mp.nstr(c['ev'], dps + 5) for c in cps],
                           sl=[{int(k): mp.nstr(v, dps + 5) for k, v in c['slice'].items()} for c in cps])
        pickle.dump(raw, open(fn, 'wb'))
    out = {}
    with mp.workdps(dps + 10):
        for mu, r in raw.items():
            out[mu] = dict(m=r['m'], w=r['w'], ev=[mp.mpf(x) for x in r['ev']], sl=[{k: mp.mpf(v) for k, v in d.items()} for d in r['sl']])
    return out


def _mpf(sx):
    return mp.mpf(sx.replace('*^', 'e'))


def _mu_of(sx):
    return tuple(int(x) for x in sx.strip('{}').split(','))


def load_Q0_spectrum(lam):
    """Q0's qsystem-spectrum-<lam>.txt (30 digits): {'branch': {(mu,k): [V,H1,a1,b1,H2,a2,b2]}, 'qpoly': {(mu,k): [q1,q2,q3 coefficient lists (ascending)]}}"""
    fn = os.path.join(REFDIR, 'qsystem-spectrum-%d%d0.txt' % (lam[0], lam[1]))
    out = dict(branch={}, qpoly={})
    with mp.workdps(50):
        for line in open(fn, encoding='utf8'):
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            f = line.split('\t')
            if f[0] == 'BRANCH':
                out['branch'][(_mu_of(f[1]), int(f[2]))] = [_mpf(x) for x in f[3:10]]
            elif f[0] == 'QPOLY':
                out['qpoly'][(_mu_of(f[1]), int(f[2]))] = [[_mpf(x) for x in s.strip('{}').split(',')] for s in f[3:6]]
    return out


def load_Q0_roots(a):
    """Q0's qsystem-roots-a7-a10-a14.txt: {mu: [V_1 < ... < V_m]} (qsSolve roots, 30 digits) for Lambda = [2a,a,0], a in (7, 10, 14)"""
    fn = os.path.join(REFDIR, 'qsystem-roots-a7-a10-a14.txt')
    out = {}
    cur = None
    with mp.workdps(50):
        for line in open(fn, encoding='utf8'):
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            f = line.split('\t')
            if f[0] == 'LAMBDA':
                cur = int(f[1])
            elif f[0] == 'ROOTS' and cur == a:
                out[_mu_of(f[1])] = [_mpf(x) for x in f[9].strip('{}').split(',')]
    return out


def load_S_pencil(lam):
    """S's order-0 pencil (generated by `gen_S_pencil`, cached): {mu: dict(mult, gdeg, gsdeg, g [exact coefficients, descending, strings], gs, V [strings, 50 digits, ascending], t, eta, q)}"""
    fn = _cache_name('Spencil', lam)
    if not os.path.exists(fn):
        gen_S_pencil(lam)
    return pickle.load(open(fn, 'rb'))


def gen_S_pencil(lam, dps=50):
    """run S's `bethe_series_scale.pencil_poly` for every irrep of V_lam (x) V_lam: the exact consistency polynomial g(b) (gcd of the three echelon systems), its squarefree part, the roots,
       V = q - b/2 - eta/2.  Slow (S: 0.5 s / 3 s / 12 s / 238 s for 420 / 630 / 840 / 1470); cached."""
    import bethe_series_scale as bs
    lam = tuple(lam)
    t0 = time.perf_counter()
    out = {}
    for mu, mult in bs.decompose(lam):
        t1 = time.perf_counter()
        g = bs.pencil_poly(lam, mu)
        gs = sp.Poly(sp.sqf_part(g.as_expr()), bs._BSYM, domain='QQ')
        eta, q = bs.eta_of(lam, mu), bs.q_of(mu)
        mp.mp.dps = dps + 20
        roots = [] if gs.degree() == 0 else list(gs.nroots(n=dps + 10, maxsteps=500))
        Vs = sorted(mp.mpf(q.numerator) / q.denominator - mp.mpf(str(sp.re(r))) / 2 - mp.mpf(eta.numerator) / eta.denominator / 2 for r in roots)
        mp.mp.dps = dps
        out[mu] = dict(mult=mult, gdeg=g.degree(), gsdeg=gs.degree(), g=[str(c) for c in g.all_coeffs()], gs=[str(c) for c in gs.all_coeffs()],
                       V=[mp.nstr(v, dps) for v in Vs], t=time.perf_counter() - t1, eta=str(eta), q=str(q))
    pickle.dump(out, open(_cache_name('Spencil', lam), 'wb'))
    return out


def gen_S_pencil_sample(lam=(28, 14, 0), mus=((42, 28, 14), (52, 17, 15), (51, 18, 15), (41, 39, 4), (41, 38, 5), (49, 28, 7)), dps=50):
    """S's pencil on a SAMPLE of irreps (the full [28,14,0] run is hours; 25-70 s per irrep here): the m = 15 irrep (42,28,14), two E23-spectral, two E12-cut and the m = 8 irrep (49,28,7).
       Cached as jacobi_cg_Spencil_sample_28_14_0.pkl (picked up by stage_recur as the source 'S pencil (sample ...)')."""
    import bethe_series_scale as bs
    lam = tuple(lam)
    irr = irreps(lam)
    out = {}
    for mu in mus:
        t1 = time.perf_counter()
        g = bs.pencil_poly(lam, mu)
        gs = sp.Poly(sp.sqf_part(g.as_expr()), bs._BSYM, domain='QQ')
        eta, q = bs.eta_of(lam, mu), bs.q_of(mu)
        mp.mp.dps = dps + 20
        roots = [] if gs.degree() == 0 else list(gs.nroots(n=dps + 10, maxsteps=500))
        Vs = sorted(mp.mpf(q.numerator) / q.denominator - mp.mpf(str(sp.re(r))) / 2 - mp.mpf(eta.numerator) / eta.denominator / 2 for r in roots)
        mp.mp.dps = dps
        out[mu] = dict(mult=irr[mu], V=[mp.nstr(v, dps) for v in Vs], t=time.perf_counter() - t1, gdeg=g.degree(), gsdeg=gs.degree(), gs=[str(c) for c in gs.all_coeffs()], eta=str(eta), q=str(q))
        pickle.dump(out, open(_cache_name('Spencil_sample', lam), 'wb'))
    return out


# ======================================================================================== 6. verification stages: J1 / J2 (closed forms, bandwidth)
def exact_full_matrix(R, lam, mu):
    """V_red on the full weight space from R's operator (slice_engine.vred_apply, built from the exact GT generators): (W = global indices, M as Fractions, M[row][col])"""
    se = _se()
    W = R.wsp.get((mu[0] - lam[0], mu[1] - lam[1], mu[2] - lam[2]))
    if W is None:
        return [], []
    pos = {int(c): a for a, c in enumerate(W)}
    w = len(W)
    M = [[Fr(0)] * w for _ in range(w)]
    for j, c in enumerate(W):
        for r, x in se.vred_apply(R, {int(c): Fr(1)}).items():
            M[pos[r]][j] = x
    return [int(c) for c in W], M


def stage_closed(reps=None, others=None, verbose=True):
    """J1/J2: bandwidth of V_red on the full weight space from R's operator, and exact equality with the closed-form chain (indices + the whole matrix), for ALL dominant candidates mu = Lam + wt
       (a superset of the irreps: includes weights of multiplicity zero) of every rep."""
    se = _se()
    reps = reps or (MAIN_REPS + [(28, 14, 0)])
    others = others if others is not None else [(20, 10, 0), (10, 3, 0), (9, 1, 0), (5, 5, 0), (7, 7, 0), (3, 2, 1), (5, 3, 2), (6, 4, 1), (6, 4, 2), (7, 5, 3), (11, 6, 3), (9, 5, 2), (12, 5, 0), (13, 13, 0)]
    tot = Counter()
    rows = []
    for lam in reps + others:
        t0 = time.perf_counter()
        R = se.Rep.of(lam)
        cands = se.candidates(lam)
        irr = irreps(lam)
        irr_ok = (se.racah_speiser(lam) == irr)
        nstate = maxband = ntri = badW = badM = 0
        wmax = 0
        for mu in cands:
            ch = chain(lam, mu)
            W, M = exact_full_matrix(R, lam, mu)
            w = len(W)
            nstate += w
            wmax = max(wmax, w)
            band = max([abs(i - j) for i in range(w) for j in range(w) if M[i][j] != 0] + [0])
            maxband = max(maxband, band)
            ntri += (band <= 1)
            if W != ch.idx:
                badW += 1
                continue
            if M != ch.matrix():
                badM += 1
        rows.append((lam, len(cands), len(irr), sum(irr.values()), nstate, wmax, maxband, ntri, badW, badM, irr_ok, time.perf_counter() - t0))
        tot['sectors'] += len(cands)
        tot['states'] += nstate
        tot['bad'] += badW + badM
        if verbose:
            print('  %-11s dominant candidates %4d (irreps %3d, copies %4d, irreps = closed-form irreps: %s), states (sum w) %6d, max w %2d: max band of V_red (R operator) = %d, tridiagonal in %d/%d, '
                  'chain-index mismatches %d, closed-form matrix mismatches %d  (%.1fs)' % (str(list(lam)), len(cands), len(irr), sum(irr.values()), irr_ok, nstate, wmax, maxband, ntri, len(cands), badW, badM,
                                                                                         rows[-1][-1]), flush=True)
    if verbose:
        print('  TOTAL: %d sectors, %d weight-space states, %d reps, closed-form mismatches %d' % (tot['sectors'], tot['states'], len(reps + others), tot['bad']))
    return rows, tot


def fmt_fr(x):
    return str(x.numerator) if x.denominator == 1 else '%d/%d' % (x.numerator, x.denominator)


def stage_tables(lam, sectors=None, verbose=True, kind=1):
    """exact tables of the chain: per irrep (descending mu) m, w, w1, t, s, m1 range, alpha_k, u_k = M[k+1,k], d_k = M[k,k+1], the products u_k d_k (k < w1-1 are > 0); kind=2: the second chain"""
    lam = tuple(lam)
    irr = irreps(lam)
    for mu, m in irr.items():
        if sectors is not None and mu not in sectors:
            continue
        ch = chain(lam, mu, kind)
        print('  mu %-12s m %d  w %d  w%d %d  t %d  s %d  m1 in [%d,%d]  cut m1 <= %d' % (str(list(mu)), m, ch.w, kind, ch.w1, ch.t, ch.s, ch.lo, ch.hi, ch.cut))
        print('      alpha = %s' % [fmt_fr(x) for x in ch.alpha])
        if ch.w > 1:
            print('      u     = %s     (M[k+1,k], up step m1 -> m1+1)' % [fmt_fr(x) for x in ch.u])
            print('      d     = %s     (M[k,k+1], down step m1+1 -> m1)' % [fmt_fr(x) for x in ch.d])
            print('      u d   = %s' % [fmt_fr(x) for x in ch.products()])


TABLES_HEADER = """# jacobi_cg_tables.txt  (agent J, round 7; generator: Paul/Python/gaudin_sl3/jacobi_cg.py  `python -X utf8 jacobi_cg.py dump_tables`)
# Exact tridiagonal chains of the reduced Gaudin vertex operator V_red on the weight spaces V_Lambda[mu-Lambda] (chain 1) and of V' = U V_red U^-1 on V_Lambda[(e3,e2,e1)] (chain 2), (e1,e2,e3) = mu - Lambda,
# for EVERY irrep mu of V_Lambda (x) V_Lambda, Lambda = {4,2,0}, {6,3,0}, {8,4,0}, {14,7,0}.  Closed forms of RESULTS-jacobi-cg.md section 3; verified entry for entry against slice_engine's exact V_red.
# RECORDS: REP {l1,l2,l3} | SECTOR mu m | CHAIN1 mu k t s lo hi w w1 cut | ALPHA1 mu {a_0,..,a_{w-1}} | UP1 mu {u_0,..,u_{w-2}} | DOWN1 mu {d_0,..,d_{w-2}} | the same with 2 for chain 2.
#   k = chain kind; t, s: fixed GT labels (t = E11 eigenvalue, s = m1+m2); position j <-> m1 = lo + j, m2 = s - m1; hi = last m1; w = hi-lo+1 = dimension of the weight space; w1 = number of positions with m1 <= cut (the invariant initial
#   segment: cut = mu1-l2 for chain 1, mu2-l3 for chain 2); the PRV subspace M_mu is the whole W1 (chain 1, mu2 >= l1+l3) or corresponds to the whole W2 (chain 2, mu2 <= l1+l3); m = multiplicity.
#   Matrix: M[j,j] = ALPHA[j], M[j+1,j] = UP[j] (m1 -> m1+1), M[j,j+1] = DOWN[j] (m1+1 -> m1), all other entries 0 (M[row,col] = coefficient of e_row in V e_col, GL3Rep plain GT basis).  Exact rationals p/q.
#   Eigenvector of eigenvalue V on W1: p_0 = 1, p_1 = (V-alpha_0)/d_0, p_{j+1} = ((V-alpha_j) p_j - u_{j-1} p_{j-1})/d_j (chain 1); chain 2 gives the physical eigenvalues when w1(chain 1) > m.
"""


def write_tables(lams=((4, 2, 0), (6, 3, 0), (8, 4, 0), (14, 7, 0)), path=None):
    """exact chains of every irrep (chain 1 and chain 2) in a TAB-separated record file (header TABLES_HEADER)"""
    path = path or os.path.join(HERE, 'jacobi_cg_tables.txt')
    TAB, NL = chr(9), chr(10)
    fl = lambda xs: '{' + ','.join(fmt_fr(x) for x in xs) + '}'
    with open(path, 'w', encoding='utf8') as f:
        f.write(TABLES_HEADER)
        for lam in lams:
            f.write(TAB.join(['REP', '{%d,%d,%d}' % tuple(lam)]) + NL)
            for mu, m in irreps(lam).items():
                mus = '{%d,%d,%d}' % mu
                f.write(TAB.join(['SECTOR', mus, str(m)]) + NL)
                for kind in (1, 2):
                    ch = chain(lam, mu, kind)
                    f.write(TAB.join(['CHAIN%d' % kind, mus] + [str(x) for x in (kind, ch.t, ch.s, ch.lo, ch.hi, ch.w, ch.w1, ch.cut)]) + NL)
                    f.write(TAB.join(['ALPHA%d' % kind, mus, fl(ch.alpha)]) + NL)
                    f.write(TAB.join(['UP%d' % kind, mus, fl(ch.u)]) + NL)
                    f.write(TAB.join(['DOWN%d' % kind, mus, fl(ch.d)]) + NL)
    return path


def verify_tables(path=None):
    """re-parse jacobi_cg_tables.txt and check, for every sector of its reps, that the continuant of the STORED chain (chain 1 on W1 if w1 = m, else chain 2 on W2) equals R's exact CHARPOLY (ref files)"""
    path = path or os.path.join(HERE, 'jacobi_cg_tables.txt')
    cur = None
    data = {}
    for line in open(path, encoding='utf8'):
        line = line.rstrip(chr(10))
        if not line or line.startswith('#'):
            continue
        f = line.split(chr(9))
        mu = tuple(int(x) for x in f[1].strip('{}').split(',')) if f[0] != 'REP' else None
        if f[0] == 'REP':
            cur = tuple(int(x) for x in f[1].strip('{}').split(','))
        elif f[0] == 'SECTOR':
            data[(cur, mu)] = {'m': int(f[2])}
        elif f[0].startswith('CHAIN'):
            data[(cur, mu)]['w1_%s' % f[0][-1]] = int(f[8])
        elif f[0].startswith(('ALPHA', 'UP', 'DOWN')):
            body = f[2].strip('{}')
            data[(cur, mu)][f[0]] = [Fr(x) for x in body.split(',')] if body else []
    ok = tot = 0
    for lam in sorted(set(k[0] for k in data)):
        R = load_R_file(lam)
        for mu, sec in R.items():
            d = data[(lam, mu)]
            kind = 1 if d['w1_1'] == d['m'] else 2
            al, up, dn = d['ALPHA%d' % kind], d['UP%d' % kind], d['DOWN%d' % kind]
            n = d['w1_%d' % kind]
            D0, D1 = [Fr(1)], [-al[0], Fr(1)]
            for k in range(1, n):
                pr = up[k - 1] * dn[k - 1]
                new = [Fr(0)] * (len(D1) + 1)
                for i, c in enumerate(D1):
                    new[i + 1] += c
                    new[i] -= al[k] * c
                for i, c in enumerate(D0):
                    new[i] -= pr * c
                D0, D1 = D1, new
            tot += 1
            ok += (D1 == sec['charpoly'])
    return ok, tot


# ======================================================================================== 7. J3: the PRV subspace inside the chain, (w, w1, m), the top-m rule
def poly_divmod(a, b):
    """exact division of polynomials with Fraction coefficients (ascending lists): (quotient, remainder)"""
    a = list(a)
    qd = len(a) - len(b)
    if qd < 0:
        return [Fr(0)], a
    q = [Fr(0)] * (qd + 1)
    for k in range(qd, -1, -1):
        c = a[k + len(b) - 1] / b[-1]
        q[k] = c
        for i, bi in enumerate(b):
            a[k + i] -= c * bi
    r = a[:len(b) - 1]
    while len(r) > 1 and r[-1] == 0:
        r.pop()
    return q, r


def poly_roots(poly, dps=60):
    """real roots (mpf, ascending) of an exact polynomial with only real roots (the continuant of a Jacobi matrix) by mpmath"""
    with mp.workdps(dps):
        cs = [mp.mpf(c.numerator) / c.denominator for c in poly][::-1]
        if len(cs) == 1:
            return []
        rts = mp.polyroots(cs, maxsteps=5000, extraprec=400)
        return sorted(mp.re(r) for r in rts)


def sign_changes(vals, rel=mp.mpf(10) ** (-25)):
    big = max([abs(v) for v in vals] + [mp.mpf(0)])
    sg = [1 if v > 0 else -1 for v in vals if abs(v) > rel * big]
    return sum(1 for a, b in zip(sg, sg[1:]) if a != b)


def constraint_targets(lam, mu):
    """does the target weight of the E12 / E23 PRV constraint exist?  E12^(l1-l2+1) maps W = V[mu-lam] to weight (t + l1-l2+1, e2 - (l1-l2+1), e3), E23^(l2-l3+1) to (t, e2 + l2-l3+1, e3 - (l2-l3+1))
       (t, e2, e3) = mu - lam.  When the target weight is not a weight of V_lam the constraint is void."""
    t, e2, e3 = mu[0] - lam[0], mu[1] - lam[1], mu[2] - lam[2]
    a1, b1 = lam[0] - lam[1] + 1, lam[1] - lam[2] + 1
    return wt_mult(lam, (t + a1, e2 - a1, e3)) > 0, wt_mult(lam, (t, e2 + b1, e3 - b1)) > 0


def analyze_rep(lam, check_prv=True, check_rule=True, nrand=100, seed=3, Rmem=None, Rfile=None, verbose=False):
    """J3 data of one rep: for every irrep mu: (w, w1, m); u at the cut = 0; products > 0 on W1 (and signs beyond the cut); M_mu inside W1 (exact PRV basis from slice_engine); the target-weight
       criteria; exact CHARPOLY relation (R's file: equal to the continuant on W1 if w1 = m, divisor of it if w1 > m); top-m rule against R's eigenvalues; Sturm oscillation of the physical
       eigenvectors; sign changes of random vectors of M_mu (variation-diminishing mechanism).  Returns a dict of counters and the distribution."""
    se = _se()
    lam = tuple(lam)
    irr = irreps(lam)
    C = Counter()
    dist = Counter()
    gaps = []
    for mu, m in irr.items():
        ch = chain(lam, mu)
        w, w1 = ch.w, ch.w1
        C['sectors'] += 1
        typ = 'E12-cut' if w > w1 else ('E23-spectral' if w1 > m else 'free')
        dist[(w, w1, m, typ)] += 1
        C['type:' + typ] += 1
        if w1 < w:
            C['cut_u_zero'] += (ch.u[w1 - 1] == 0)
            neg = sum(1 for k in range(w1, w - 1) if ch.u[k] * ch.d[k] < 0)
            zer = sum(1 for k in range(w1, w - 1) if ch.u[k] * ch.d[k] == 0)
            C['beyond_cut_products_neg'] += neg
            C['beyond_cut_products_zero'] += zer
            C['beyond_cut_products_pos'] += (w - 1 - w1 - neg - zer)
            C['cut_sectors'] += 1
        C['products_pos_on_W1'] += all(p > 0 for p in ch.products(w1))
        C['products_total_W1'] += max(w1 - 1, 0)
        # signs of u and d individually on W1 (all negative: Perron-Frobenius / Sturm in the plain basis)
        C['u_neg_W1'] += all(x < 0 for x in ch.u[:max(w1 - 1, 0)])
        C['d_neg_W1'] += all(x < 0 for x in ch.d[:max(w1 - 1, 0)])
        ex12, ex23 = constraint_targets(lam, mu)
        C['e12_iff_target'] += ((w > w1) == ex12)
        C['e23_implies_target'] += ((not (w1 > m)) or ex23)
        C['never_both'] += (not (w > w1 and w1 > m))
        C['m_le_w1'] += (m <= w1)
        # regimes: the E23 image weight needs mu2 <= l1+l3-1 and the E12 image weight needs mu2 >= l1+l3+1 (E22 of the image must lie in [l3, l1]), so for mu2 >= l1+l3 only E12 can act
        # (M_mu = W1, m = w1) and for mu2 <= l1+l3 only E23 can act (W1 = W): the two constraints are never both active
        C['regA_m_eq_w1'] += ((mu[1] < lam[0] + lam[2]) or (m == w1))
        C['regB_w_eq_w1'] += ((mu[1] > lam[0] + lam[2]) or (w == w1))
        if check_prv:
            p = se.prv_space(lam, mu)
            C['prv_dim_eq_m'] += (p.m == m)
            C['W_eq_chain'] += ([int(c) for c in p.W] == ch.idx)
            inW1 = all(all(x == 0 for x in b[w1:]) for b in p.basis)
            C['prv_in_W1'] += inW1
            if w1 > m and p.m > 0:
                # random vectors of M_mu: sign changes >= w1 - m ?
                rng = random.Random(seed + w + 31 * m)
                ok = True
                minsc = 10**9
                for _ in range(nrand):
                    cs = [rng.randint(-6, 6) for _ in range(p.m)]
                    if all(c == 0 for c in cs):
                        continue
                    v = [sum(cs[k] * p.basis[k][i] for k in range(p.m)) for i in range(w1)]
                    sc = sign_changes([mp.mpf(x.numerator) / x.denominator for x in v])
                    minsc = min(minsc, sc)
                    ok &= (sc >= w1 - m)
                C['kernel_vectors_signchanges_ok'] += ok
                C['kernel_sectors'] += 1
                C['kernel_min_sc_margin'] = min(C.get('kernel_min_sc_margin', 10**9), minsc - (w1 - m))
        if Rfile is not None and mu in Rfile:
            P = Rfile[mu]['charpoly']
            C1 = ch.contin()
            if w1 == m:
                C['charpoly_equal_w1_eq_m'] += (C1 == P)
                C['charpoly_checked_w1_eq_m'] += 1
            else:
                qq, rr = poly_divmod(C1, P)
                C['charpoly_divides_w1_gt_m'] += (all(r == 0 for r in rr))
                C['charpoly_checked_w1_gt_m'] += 1
                rp, rq = poly_roots(P), poly_roots(qq)
                if rq:
                    gaps.append(min(rp) - max(rq))
        if check_rule and Rmem is not None and mu in Rmem:
            ev = Rmem[mu]['ev']
            phys = physical_V(lam, mu, dps=50, ch=ch, m=m, method='top')
            dev = max([abs(a - b) for a, b in zip(phys, ev)] + [mp.mpf(0)])
            C['top_rule_checked'] += 1
            C['top_rule_ok'] += (dev < mp.mpf(10) ** (-35))
            C['top_rule_maxdev_log10'] = max(C.get('top_rule_maxdev_log10', -999), float(mp.log10(dev)) if dev > 0 else -999)
            phys2 = physical_V(lam, mu, dps=50, ch=ch, m=m, method='auto')
            dev2 = max([abs(a - b) for a, b in zip(phys2, ev)] + [mp.mpf(0)])
            C['auto_rule_ok'] += (dev2 < mp.mpf(10) ** (-35))
            C['auto_rule_maxdev_log10'] = max(C.get('auto_rule_maxdev_log10', -999), float(mp.log10(dev2)) if dev2 > 0 else -999)
            if w1 > m:
                # Sturm oscillation: physical eigenvector number i (0-based among the w1 eigenvalues) has exactly i sign changes
                allev = jacobi_eigs(lam, mu, 50, ch)
                for V, sl in zip(ev, Rmem[mu]['sl']):
                    i = min(range(len(allev)), key=lambda k: abs(allev[k] - V))
                    pvec = slice_from_V(lam, mu, V, 50, ch)[:w1]
                    C['sturm_checked'] += 1
                    C['sturm_ok'] += (sign_changes(pvec) == i)
                    C['phys_index_is_top'] += (i >= w1 - m)
    return C, dist, gaps


def stage_struct(reps=None, verbose=True):
    se = _se()
    reps = reps or (MAIN_REPS + [(28, 14, 0)])
    allrows = {}
    for lam in reps:
        t0 = time.perf_counter()
        Rmem = load_R_mem(lam, 50)
        Rfile = load_R_file(lam) if lam in MAIN_REPS else None
        C, dist, gaps = analyze_rep(lam, Rmem=Rmem, Rfile=Rfile)
        allrows[lam] = (C, dist, gaps)
        if verbose:
            print('  %s: %d irreps; types: free (w = w1 = m) %d, E12-cut (w > w1 = m) %d, E23-spectral (w = w1 > m) %d; never both: %d/%d' %
                  (str(list(lam)), C['sectors'], C['type:free'], C['type:E12-cut'], C['type:E23-spectral'], C['never_both'], C['sectors']))
            print('     u at the cut = 0 in %d/%d cut sectors (products beyond the cut: %d negative, %d zero, %d positive); off-diagonal products on W1 all > 0 in %d/%d sectors (%d products); '
                  'u<0 and d<0 on W1 in %d/%d' % (C['cut_u_zero'], C['cut_sectors'], C['beyond_cut_products_neg'], C['beyond_cut_products_zero'], C['beyond_cut_products_pos'],
                                                  C['products_pos_on_W1'], C['sectors'], C['products_total_W1'], min(C['u_neg_W1'], C['d_neg_W1']), C['sectors']))
            print('     PRV: dim M = m in %d/%d, weight space = chain in %d/%d, M_mu inside W1 in %d/%d; m <= w1: %d/%d; E12 effective <=> target weight exists: %d/%d; E23 effective => target weight exists: %d/%d; '
                  'mu2 >= l1+l3 => m = w1: %d/%d; mu2 <= l1+l3 => w = w1: %d/%d' %
                  (C['prv_dim_eq_m'], C['sectors'], C['W_eq_chain'], C['sectors'], C['prv_in_W1'], C['sectors'], C['m_le_w1'], C['sectors'], C['e12_iff_target'], C['sectors'],
                   C['e23_implies_target'], C['sectors'], C['regA_m_eq_w1'], C['sectors'], C['regB_w_eq_w1'], C['sectors']))
            if Rfile is not None:
                print('     exact CHARPOLY (R file): continuant on W1 == CHARPOLY in %d/%d sectors with w1 = m; CHARPOLY divides the continuant in %d/%d sectors with w1 > m; smallest gap (min physical root - max other root) %s'
                      % (C['charpoly_equal_w1_eq_m'], C['charpoly_checked_w1_eq_m'], C['charpoly_divides_w1_gt_m'], C['charpoly_checked_w1_gt_m'], mp.nstr(min(gaps), 6) if gaps else '-'))
            print("     rigorous rule (chain 1 on W1 if w1 = m, else chain 2 on W2) vs R's 50-digit eigenvalues: %d/%d sectors agree (max deviation 1e%d); empirical top-m rule of chain 1: %d/%d (max deviation 1e%d); "
                  "Sturm oscillation (physical eigenvector i has i sign changes) %d/%d; physical eigenvalues are the top ones: %d/%d; "
                  "random vectors of M_mu have >= w1-m sign changes in %d/%d E23-sectors (min margin %s)  [%.1fs]" %
                  (C['auto_rule_ok'], C['top_rule_checked'], C['auto_rule_maxdev_log10'], C['top_rule_ok'], C['top_rule_checked'], C['top_rule_maxdev_log10'], C['sturm_ok'], C['sturm_checked'],
                   C['phys_index_is_top'], C['sturm_checked'], C['kernel_vectors_signchanges_ok'], C['kernel_sectors'], C.get('kernel_min_sc_margin'), time.perf_counter() - t0), flush=True)
    # the (w, w1, m) table of [14,7,0]
    lam = (14, 7, 0)
    if lam in allrows and verbose:
        dist = allrows[lam][1]
        print('  (w, w1, m) distribution of the 113 irreps of [14,7,0]^2 (count; type):')
        for (w, w1, m, typ), n in sorted(dist.items()):
            print('     (%d, %d, %d)  x %2d  %s' % (w, w1, m, n, typ))
        wm = Counter()
        for (w, w1, m, typ), n in dist.items():
            wm[(w, m)] += n
        print('  (w, m) distribution: %s' % sorted(wm.items()))
    return allrows


def stage_scan(lmax=9, verbose=True):
    """broad exact scan: every rep [l1,l2,l3] with l1 <= lmax (l1 >= 1): closed-form chain 1 = R's exact V_red matrix and closed-form chain 2 = the exact V' matrix on every dominant candidate, bandwidth 1,
       products > 0 on W1 / W2, u(cut) = 0, M_mu inside W1, E12 effective <=> target exists, never both constraints effective, the regime statements (mu2 >= l1+l3 => w1 = m, mu2 <= l1+l3 => w2 = m),
       and both selection rules (rigorous chain-1/chain-2 rule, empirical top-m rule) against the eigenvalues of R's exact V_red on M_mu (float64)."""
    se = _se()
    tot = Counter()
    t0 = time.perf_counter()
    worst = worst_top = 0.0
    nreps = 0
    for l1 in range(1, lmax + 1):
        for l2 in range(0, l1 + 1):
            for l3 in range(0, l2 + 1):
                lam = (l1, l2, l3)
                nreps += 1
                R = se.Rep.of(lam)
                irr = irreps(lam)
                tot['irreps_eq'] += (se.racah_speiser(lam) == irr)
                for mu in se.candidates(lam):
                    ch = chain(lam, mu)
                    W, M = exact_full_matrix(R, lam, mu)
                    tot['candidates'] += 1
                    tot['closed_ok'] += (W == ch.idx and M == ch.matrix())
                    ch2 = chain(lam, mu, 2)
                    W2, M2 = exact_full_matrix2(R, lam, mu)
                    tot['closed2_ok'] += (W2 == ch2.idx and M2 == ch2.matrix())
                    tot['mult_closed_ok'] += (mult_closed(lam, mu) == rs_mult(lam, mu))
                for mu, m in irr.items():
                    ch = chain(lam, mu)
                    ch2 = chain(lam, mu, 2)
                    w, w1, w2 = ch.w, ch.w1, ch2.w1
                    tot['irrep_sectors'] += 1
                    tot['copies'] += m
                    tot['prod_pos'] += all(p > 0 for p in ch.products(w1))
                    tot['prod2_pos'] += all(p > 0 for p in ch2.products(w2))
                    tot['cut_ok'] += (w1 == w or ch.u[w1 - 1] == 0)
                    tot['cut2_ok'] += (w2 == ch2.w or ch2.u[w2 - 1] == 0)
                    p = se.prv_space(lam, mu)
                    tot['prv_ok'] += (p.m == m and all(all(x == 0 for x in b[w1:]) for b in p.basis))
                    ex12, ex23 = constraint_targets(lam, mu)
                    tot['e12_iff'] += ((w > w1) == ex12)
                    tot['both_ok'] += (not (w > w1 and w1 > m))
                    tot['m_le_w1'] += (m <= w1)
                    tot['regA'] += ((mu[1] < lam[0] + lam[2]) or (m == w1))
                    tot['regB'] += ((mu[1] > lam[0] + lam[2]) or (m == w2))
                    tot['type_e12'] += (w > w1)
                    tot['type_e23'] += (w1 > m)
                    Vm = se.V_red(lam, mu, p)
                    ev = np.sort(np.linalg.eigvals(np.array([[float(x) for x in row] for row in Vm])).real)
                    phys = [float(x) for x in physical_V(lam, mu, dps=30, ch=ch, m=m)]
                    dev = max(abs(a - b) / max(1.0, abs(a)) for a, b in zip(ev, phys))
                    worst = max(worst, dev)
                    tot['auto_ok'] += (dev < 1e-8)
                    phys = [float(x) for x in physical_V(lam, mu, dps=30, ch=ch, m=m, method='top')]
                    dev = max(abs(a - b) / max(1.0, abs(a)) for a, b in zip(ev, phys))
                    worst_top = max(worst_top, dev)
                    tot['top_ok'] += (dev < 1e-8)
    if verbose:
        print("  scan of all %d reps [l1,l2,l3], 1 <= l1 <= %d, l3 >= 0 (%.1fs): irreps = closed-form irreps in %d/%d reps; closed-form chain 1 == R exact V_red matrix: %d/%d dominant candidates, chain 2 == exact V' matrix: %d/%d; "
              "%d irrep sectors (%d copies): products on W1 > 0: %d, on W2 > 0: %d; u(cut1) = 0: %d, u'(cut2) = 0: %d; PRV dim = m and M_mu inside W1: %d; E12 effective <=> target weight exists: %d; never both effective: %d; "
              "m <= w1: %d; mu2 >= l1+l3 => m = w1: %d; mu2 <= l1+l3 => m = w2: %d; sectors with E12-cut %d, with E23-spectral %d; selection vs eigenvalues of R's V_red on M_mu: rigorous rule %d/%d (max rel. dev. %.1e), "
              "top-m rule %d/%d (%.1e); closed-form multiplicity (w1 or w2) == Racah-Speiser on %d/%d dominant candidates" %
              (nreps, lmax, time.perf_counter() - t0, tot['irreps_eq'], nreps, tot['closed_ok'], tot['candidates'], tot['closed2_ok'], tot['candidates'], tot['irrep_sectors'], tot['copies'], tot['prod_pos'],
               tot['prod2_pos'], tot['cut_ok'], tot['cut2_ok'], tot['prv_ok'], tot['e12_iff'], tot['both_ok'], tot['m_le_w1'], tot['regA'], tot['regB'], tot['type_e12'], tot['type_e23'],
               tot['auto_ok'], tot['irrep_sectors'], worst, tot['top_ok'], tot['irrep_sectors'], worst_top, tot['mult_closed_ok'], tot['candidates']))
    return tot


# ======================================================================================== 8. J4: the recurrence against R's slices, from every V source
def r_slices_file(lam):
    """R's COPY records from the ref file as {mu: dict(ev [mpf], sl [{0-based idx: mpf}])} (30 digits; the file indices are 1-based)"""
    sec = load_R_file(lam)
    return {mu: dict(ev=[c[0] for c in s['copies']], sl=[{int(k) - 1: v for k, v in c[1].items()} for c in s['copies']]) for mu, s in sec.items()}


def compare_slices(lam, Rdat, Vsrc, dps, ev_ref=None):
    """slice_from_V(lam, mu, V) for every copy against the reference slices Rdat (ev, sl): Vsrc = {mu: [V,..]} (ascending V per sector, same copies as Rdat).
       Returns a dict: ncopies, max absolute deviation (componentwise, all entries), max deviation relative to max|ref slice|, worst copy, support violations (reference entries outside the chain/W1),
       max |V - V_ref| (if ev_ref given) and max amplification (relative slice deviation / |V - V_ref|), time of the recurrences (mp, dps) in seconds."""
    nc = 0
    maxabs = maxrel = mp.mpf(0)
    worst = None
    supp_bad = 0
    maxdV = mp.mpf(0)
    maxamp = mp.mpf(0)
    ttot = 0.0
    with mp.workdps(dps + 10):
        for mu, rd in Rdat.items():
            ch = chain(lam, mu)
            Vs = Vsrc[mu]
            assert len(Vs) == len(rd['sl']), (lam, mu, len(Vs), len(rd['sl']))
            for k, (V, sl) in enumerate(zip(Vs, rd['sl'])):
                t0 = time.perf_counter()
                p = slice_from_V(lam, mu, V, dps, ch)
                ttot += time.perf_counter() - t0
                nc += 1
                pos = {i: j for j, i in enumerate(ch.idx)}
                dev = mp.mpf(0)
                for j, i in enumerate(ch.idx):
                    dev = max(dev, abs(p[j] - sl.get(i, 0)))
                supp_bad += sum(1 for i, v in sl.items() if (i not in pos or pos[i] >= ch.w1) and v != 0)
                scale = max(abs(v) for v in sl.values())
                rel = dev / scale
                if rel > maxrel:
                    maxrel, worst = rel, (mu, k + 1)
                maxabs = max(maxabs, dev)
                if ev_ref is not None:
                    dV = abs(mp.mpf(V) - ev_ref[mu][k])
                    maxdV = max(maxdV, dV)
                    if dV > 0:
                        maxamp = max(maxamp, rel / dV)
    return dict(ncopies=nc, maxabs=maxabs, maxrel=maxrel, worst=worst, supp_bad=supp_bad, maxdV=maxdV, maxamp=maxamp, time=ttot)


def _sources(lam):
    """the V sources available for lam: dict name -> {mu: [V ascending]} (mpf)"""
    src = {}
    a = lam[1] if (lam[0] == 2 * lam[1] and lam[2] == 0) else None
    if lam in [(6, 3, 0), (8, 4, 0), (14, 7, 0)] and os.path.exists(os.path.join(REFDIR, 'qsystem-spectrum-%d%d0.txt' % (lam[0], lam[1]))):
        q0 = load_Q0_spectrum(lam)
        d = defaultdict(list)
        for (mu, k), vals in sorted(q0['branch'].items(), key=lambda kv: (kv[0][0], kv[0][1])):
            d[mu].append((k, vals[0]))
        src['Q0 spectrum file'] = {mu: [v for k, v in sorted(lst)] for mu, lst in d.items()}
    if a in (7, 10, 14):
        src['Q0 roots file'] = load_Q0_roots(a)
    if os.path.exists(_cache_name('Spencil', lam)):
        sp_ = load_S_pencil(lam)
        with mp.workdps(60):
            src["S pencil (b -> V)"] = {mu: [mp.mpf(x) for x in r['V']] for mu, r in sp_.items()}
    fs = _cache_name('Spencil_sample', lam)
    if os.path.exists(fs):
        sp_ = pickle.load(open(fs, 'rb'))
        with mp.workdps(60):
            src["S pencil (sample of %d irreps)" % len(sp_)] = {mu: [mp.mpf(x) for x in r['V']] for mu, r in sp_.items()}
    return src


def stage_recur(reps=None, verbose=True, hi=True):
    reps = reps or [(4, 2, 0), (6, 3, 0), (8, 4, 0), (14, 7, 0), (20, 10, 0), (28, 14, 0)]
    rows = []
    for lam in reps:
        irr = irreps(lam)
        Rmem = load_R_mem(lam, 50)
        ncop = sum(irr.values())
        # chains (cold build)
        _CH.clear()
        t0 = time.perf_counter()
        for mu in irr:
            chain(lam, mu)
        t_chain = time.perf_counter() - t0
        for mu in irr:
            chain(lam, mu).mp_coeffs(50)
        # sources, reference = R in memory at 50 digits
        ev_ref = {mu: r['ev'] for mu, r in Rmem.items()}
        srcs = {'R eigenvalues (50 digits)': {mu: r['ev'] for mu, r in Rmem.items()}}
        t0 = time.perf_counter()
        srcs['chain only: chain 1 on W1 / chain 2 on W2 (Q-free)'] = {mu: physical_V(lam, mu, 50, chain(lam, mu), m) for mu, m in irr.items()}
        t_phys = time.perf_counter() - t0
        srcs['chain only: top-m of V_red|W1 (Q-free)'] = {mu: physical_V(lam, mu, 50, chain(lam, mu), m, method='top') for mu, m in irr.items()}
        qsrc = _sources(lam)
        for nm in list(qsrc):
            if nm.startswith('Q0'):
                srcs[nm + ' + Newton polish on the chain (50 digits)'] = {mu: [refine_V(lam, mu, V, 50, chain(lam, mu)) for V in Vs] for mu, Vs in qsrc[nm].items()}
        srcs.update(qsrc)
        if verbose:
            print('  %s: %d irreps, %d copies; chain coefficients (closed forms, all irreps): %.3f s; eigenvalues of V_red|W1 for all irreps (mpmath eigsy, 60 digits): %.2f s' % (str(list(lam)), len(irr), ncop, t_chain, t_phys))
        for name, vs in srcs.items():
            r = compare_slices(lam, {mu: Rmem[mu] for mu in vs}, vs, 50, ev_ref)
            rows.append((lam, 'R mem (50 digits)', name, r))
            if verbose:
                print('     V from %-40s vs R in memory (50 digits): %4d copies, max |dev| %s, max relative dev %s (worst %s), support violations %d, max |V-V_R| %s, max amplification %s, recurrence time %.3f s (%.0f us/copy at 50 digits)' %
                      (name, r['ncopies'], mp.nstr(r['maxabs'], 3), mp.nstr(r['maxrel'], 3), r['worst'], r['supp_bad'], mp.nstr(r['maxdV'], 3), mp.nstr(r['maxamp'], 3), r['time'], 1e6 * r['time'] / r['ncopies']), flush=True)
        if lam in MAIN_REPS:
            Rf = r_slices_file(lam)
            ev_f = {mu: r['ev'] for mu, r in Rf.items()}
            srcs_f = {'R file eigenvalues (30 digits)': ev_f}
            for nm, vs in srcs.items():
                if nm.startswith('Q0') or nm.startswith('S pencil'):
                    srcs_f[nm] = vs
            for name, vs in srcs_f.items():
                r = compare_slices(lam, {mu: Rf[mu] for mu in vs}, vs, 50, ev_f)
                rows.append((lam, 'R file (30 digits)', name, r))
                if verbose:
                    print('     V from %-40s vs R FILE slices (30 digits):     %4d copies, max |dev| %s, max relative dev %s (worst %s), support violations %d, recurrence time %.3f s' %
                          (name, r['ncopies'], mp.nstr(r['maxabs'], 3), mp.nstr(r['maxrel'], 3), r['worst'], r['supp_bad'], r['time']), flush=True)
    # other shapes and lambda3 > 0 against R's engine in memory (no ref files exist for them), V from the Q-free rule
    if hi:
        tot_c = 0
        worst_rel = mp.mpf(0)
        for lam in [(10, 3, 0), (9, 1, 0), (12, 5, 0), (3, 2, 1), (5, 3, 2), (6, 4, 1), (6, 4, 2), (7, 5, 3), (11, 6, 3), (9, 5, 2), (12, 7, 3), (15, 8, 4)]:
            irr = irreps(lam)
            Rm = load_R_mem(lam, 50)
            vs = {mu: physical_V(lam, mu, 50, chain(lam, mu), m) for mu, m in irr.items()}
            r = compare_slices(lam, Rm, vs, 50, {mu: x['ev'] for mu, x in Rm.items()})
            tot_c += r['ncopies']
            worst_rel = max(worst_rel, r['maxrel'])
            rows.append((lam, 'R mem (50 digits)', 'chain only (Q-free)', r))
            if verbose:
                print('  %-10s (other shape / lambda3 > 0): %4d copies, chain-only V: max relative dev %s, max |V-V_R| %s, support violations %d' % (str(list(lam)), r['ncopies'], mp.nstr(r['maxrel'], 3), mp.nstr(r['maxdV'], 3), r['supp_bad']), flush=True)
        if verbose:
            print('  other shapes / lambda3 > 0 together: %d copies, worst relative deviation %s' % (tot_c, mp.nstr(worst_rel, 3)))
    # 80-digit check: the 1e-42 floor of the 50-digit comparison at a = 14 is R's mp.eig accuracy, not the recurrence
    if hi:
        for lam in [(14, 7, 0), (28, 14, 0)]:
            irr = irreps(lam)
            Rm = load_R_mem(lam, 80)
            vs = {mu: physical_V(lam, mu, 80, chain(lam, mu), m) for mu, m in irr.items()}
            r = compare_slices(lam, Rm, vs, 80, {mu: x['ev'] for mu, x in Rm.items()})
            rows.append((lam, 'R mem (80 digits)', 'chain only (Q-free), 80 digits', r))
            if verbose:
                print('  %s at 80 digits (R in memory 80 digits, V = physical eigenvalues of chain 1 / chain 2 at 80 digits): %d copies, max |dev| %s, max relative dev %s, max |V-V_R| %s' %
                      (str(list(lam)), r['ncopies'], mp.nstr(r['maxabs'], 3), mp.nstr(r['maxrel'], 3), mp.nstr(r['maxdV'], 3)), flush=True)
    return rows


# ======================================================================================== 9. J5: the Q side (S's Baxter pencil, Q0's Q polynomials)
def poly_of_V_from_b(gb, q, eta):
    """exact polynomial in V (ascending Fractions) from a polynomial in b (descending coefficient strings 'p/q'), via b = 2q - 2V - eta (V = q - b/2 - eta/2); returned monic"""
    A = 2 * Fr(q) - Fr(eta)
    # Horner in the variable V with b = A - 2V
    res = [Fr(0)]
    for c in gb:
        # res = res * (A - 2V) + c
        new = [Fr(0)] * (len(res) + 1)
        for i, r in enumerate(res):
            new[i] += r * A
            new[i + 1] -= 2 * r
        new[0] += Fr(c)
        res = new
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    lead = res[-1]
    return [x / lead for x in res]


def baxter_vs_chain_roots(lam, mu):
    """roots (as integers/half-integers, with multiplicity) of the off-diagonal PRODUCTS of the two three-term recursions for the sector mu:
       Baxter (S): alpha_N rho_{N+1} = (N-mu3)(N-mu2-1)(N-mu1-2) (N+1)(N-l2)(N-l1-1)   (polynomial in N, chain N = 0..d_j, d_j = (mu3, mu2+1, mu1+2));
       GT chain (here): u(m1) d(m1+1) = N10(m1) / ((x+1)(x+2)^2(x+3)),  N10 = prod of the ten linear factors listed below, x = 2 m1 - s."""
    l1, l2, l3 = lam
    ch = chain(lam, mu)
    t, s = ch.t, ch.s
    bax = sorted([mu[2], mu[1] + 1, mu[0] + 2, -1, l2, l1 + 1])
    gt_num = sorted([l1, l2 - 1, l3 - 2, s - t - 1, mu[0] - l2, s - l2 - 1, t - 1, s - l1 - 2, s - l3, s - mu[0] + l2 - 2])
    gt_den = sorted([Fr(s - 1, 2), Fr(s - 2, 2), Fr(s - 2, 2), Fr(s - 3, 2)])
    return bax, gt_num, gt_den


def _np(x):
    return float(x)


def stage_qside(verbose=True):
    se = _se()
    out = {}
    # ---- (1) the two recurrences side by side: sizes and off-diagonal products
    if verbose:
        print('  (1) the two three-term recursions, sector examples (chain sizes and the roots of the off-diagonal products; Baxter: roots in N, GT: roots in m1):')
    lam = (14, 7, 0)
    for mu in [(28, 14, 0), (21, 14, 7), (22, 10, 10), (14, 14, 14), (24, 12, 6)]:
        ch = chain(lam, mu)
        bax, gnum, gden = baxter_vs_chain_roots(lam, mu)
        if verbose:
            print('     mu %-12s m %d: GT chain length w = %d (W1: %d); Baxter chains d_j + 1 = %s;  Baxter product roots %s;  GT numerator roots %s, denominator roots %s' %
                  (str(list(mu)), irreps(lam)[mu], ch.w, ch.w1, [mu[2] + 1, mu[1] + 2, mu[0] + 3], bax, gnum, [str(x) for x in gden]))
    sizes = [(chain(lam, mu).w1, mu[2] + 1, mu[1] + 2, mu[0] + 3) for mu in irreps(lam)]
    if verbose:
        print('     [14,7,0]: max GT chain length w1 = %d, Baxter chain lengths d_1+1 up to %d, d_2+1 up to %d, d_3+1 up to %d; w1 = d_1+1 in %d of 113 sectors' %
              (max(x[0] for x in sizes), max(x[1] for x in sizes), max(x[2] for x in sizes), max(x[3] for x in sizes), sum(1 for x in sizes if x[0] == x[1])))
    # ---- (2) exact polynomial identity: S's consistency polynomial (in V) == R's CHARPOLY (and the continuant on W1 when w1 = m)
    if verbose:
        print("  (2) S's echelon-pencil consistency polynomial g(b) (gcd of the three exact systems), mapped to V by b = 2q - 2V - eta, vs R's exact CHARPOLY and the chain continuant:")
    for lam in MAIN_REPS:
        fn = _cache_name('Spencil', lam)
        if not os.path.exists(fn):
            continue
        Sp = load_S_pencil(lam)
        Rf = load_R_file(lam)
        neq = nsqf = ndiv = n1 = n2 = 0
        tot_t = 0.0
        for mu, m in irreps(lam).items():
            r = Sp[mu]
            P = poly_of_V_from_b(r['gs'], r['q'], r['eta'])
            neq += (P == Rf[mu]['charpoly'])
            nsqf += (r['gdeg'] == r['gsdeg'])
            ch = chain(lam, mu)
            C1 = ch.contin()
            if ch.w1 == m:
                n1 += (poly_of_V_from_b(r['gs'], r['q'], r['eta']) == C1)
            else:
                qq, rr = poly_divmod(C1, P)
                n2 += all(x == 0 for x in rr)
            tot_t += r['t']
        nw1 = sum(1 for mu, m in irreps(lam).items() if chain(lam, mu).w1 == m)
        if verbose:
            print('     %s: g -> V polynomial == R CHARPOLY (exact, monic) in %d/%d sectors; g already squarefree in %d/%d; continuant on W1 == Q-side polynomial in %d/%d sectors with w1 = m; Q-side polynomial divides the '
                  'continuant in %d/%d sectors with w1 > m (S pencil time %.1f s)' % (str(list(lam)), neq, len(irreps(lam)), nsqf, len(irreps(lam)), n1, nw1, n2, len(irreps(lam)) - nw1, tot_t))
    # ---- (3) slice vs Q-polynomial coefficients (Q0): is slice_j^(c) = f_j * q_i^(c) (or f_j * q_i^(c)/q_i0^(c)) with a copy-independent factor f_j, over the copies c of a sector?
    if verbose:
        print("  (3) regression of R's slices (R file, 30 digits) against Q0's QPOLY coefficients (q1,q2,q3): ratio constancy over the copies c = 1..m of each sector with m >= 2")
        print("      (test A: z_c = s_j^(c) / q_i^(c) constant;  test B: z_c = s_j^(c) q_i0^(c) / q_i^(c) constant;  'constant' = (max z - min z)/max|z| < 1e-7;  j >= 2 slice components; the Q columns q_i are: the monomial coefficients "
              "(= Taylor coefficients at the site u=0), the Taylor coefficients at the site u=1, and the values at u = -1, 2, ..., 10, of q1, q2, q3)")
    TOL = 1e-7
    for lam in [(6, 3, 0), (8, 4, 0), (14, 7, 0)]:
        Rf = r_slices_file(lam)
        Q0 = load_Q0_spectrum(lam)
        nA = nB = hA = hB = nsec = ncop = 0
        hits = []
        for mu, rd in Rf.items():
            m = len(rd['sl'])
            if m < 2:
                continue
            ch = chain(lam, mu)
            nsec += 1
            ncop += m
            S = np.array([[float(rd['sl'][c].get(ch.idx[j], 0)) for c in range(m)] for j in range(1, ch.w1)])
            qs = [Q0['qpoly'][(mu, c + 1)] for c in range(m)]
            Qrows, labels = [], []
            for k in range(3):
                for i in range(len(qs[0][k])):
                    col = np.array([float(qs[c][k][i]) for c in range(m)])
                    if np.max(np.abs(col)) > 0:
                        Qrows.append(col)
                        labels.append(('c', k + 1, i))                                   # coefficient of u^i = Taylor coefficient at the site u = 0
            with mp.workdps(40):
                for k in range(3):
                    n_ = len(qs[0][k])
                    for i in range(n_):                                                  # Taylor coefficients at the site u = 1: coefficient of v^i in q(1+v)
                        col = np.array([float(mp.fsum(qs[c][k][j] * math.comb(j, i) for j in range(i, n_))) for c in range(m)])
                        if np.max(np.abs(col)) > 0:
                            Qrows.append(col)
                            labels.append(('t1', k + 1, i))
                    for N in (-1, 2, 3, 4, 5, 6, 7, 8, 9, 10):                          # values at integer points
                        col = np.array([float(mp.fsum(qs[c][k][j] * mp.mpf(N) ** j for j in range(n_))) for c in range(m)])
                        if np.max(np.abs(col)) > 0:
                            Qrows.append(col)
                            labels.append(('v', k + 1, N))
            Qm = np.array(Qrows)
            okq = np.where(np.all(np.abs(Qm) > 1e-12 * np.max(np.abs(Qm), axis=1)[:, None], axis=1))[0]
            for j in range(len(S)):
                sv = S[j]
                if np.max(np.abs(sv)) == 0 or np.any(np.abs(sv) < 1e-12 * np.max(np.abs(sv))):
                    continue
                Qok = Qm[okq]
                labok = [labels[i] for i in okq]
                Z = sv[None, :] / Qok                                                     # test A, vectorised over the Q columns
                spread = (Z.max(axis=1) - Z.min(axis=1)) / np.abs(Z).max(axis=1)
                nA += len(okq)
                for ii in np.where(spread < TOL)[0]:
                    hA += 1
                    hits.append(('A', mu, m, j + 2, labok[ii], None, j + 2 == ch.w1))
                for ii0 in range(len(okq)):                                               # test B: s_j q_i0 / q_i over all i != i0
                    Z = (sv * Qok[ii0])[None, :] / Qok
                    spread = (Z.max(axis=1) - Z.min(axis=1)) / np.abs(Z).max(axis=1)
                    spread[ii0] = np.inf
                    nB += len(okq) - 1
                    for ii in np.where(spread < TOL)[0]:
                        hB += 1
                        hits.append(('B', mu, m, j + 2, labok[ii], labok[ii0], j + 2 == ch.w1))
        hit_secs = sorted(set((h[1], h[2]) for h in hits))
        n_last = sum(1 for h in hits if h[6])
        n_other = len(hits) - n_last
        hit_sec_other = sorted(set((h[1], h[2], h[3]) for h in hits if not h[6]))
        out[lam] = dict(nsec=nsec, ncop=ncop, testsA=nA, hitsA=hA, testsB=nB, hitsB=hB, hit_sectors=hit_secs, hits_last=n_last, hits_other=n_other)
        if verbose:
            print('     %s: %d sectors with m >= 2 (%d copies): test A %d pairs, %d hits; test B %d pairs, %d hits; hits in %d sectors; hits at the LAST slice component (j = w1): %d, at any other component: %d %s' %
                  (str(list(lam)), nsec, ncop, nA, hA, nB, hB, len(hit_secs), n_last, n_other, hit_sec_other[:6]))
    # ---- (4) the last component: ratio of two Taylor coefficients of ONE q_k at ONE site (u = 0 or u = 1)
    if verbose:
        print("  (4) the LAST slice component s_w1 (largest m1) against ratios s_w1 * q_k^[a](site) / q_k^[b](site) = const over the copies (Taylor coefficients of orders a != b of one q_k at the site u = 0 or u = 1):")
    for lam in [(6, 3, 0), (8, 4, 0), (14, 7, 0)]:
        Rf = r_slices_file(lam)
        Q0 = load_Q0_spectrum(lam)
        nsec = ncov = 0
        fam = Counter()
        unc = []
        with mp.workdps(40):
            for mu, rd in Rf.items():
                m = len(rd['sl'])
                if m < 2:
                    continue
                ch = chain(lam, mu)
                nsec += 1
                sl_last = np.array([float(rd['sl'][c].get(ch.idx[ch.w1 - 1], 0)) for c in range(m)])
                if np.any(sl_last == 0):
                    unc.append(mu)
                    continue
                found = []
                for site in (0, 1):
                    for k in range(3):
                        polys = [Q0['qpoly'][(mu, c + 1)][k] for c in range(m)]
                        n_ = len(polys[0])
                        T = np.array([[float(mp.fsum(polys[c][j] * math.comb(j, i) * mp.mpf(site) ** (j - i) for j in range(i, n_))) for i in range(n_)] for c in range(m)])   # (m, n)
                        ok = np.where((np.max(np.abs(T), axis=0) > 0) & np.all(np.abs(T) > 1e-25 * np.maximum(np.max(np.abs(T), axis=0), 1e-300)[None, :], axis=0))[0]
                        if len(ok) < 2:
                            continue
                        Tn = T[:, ok]
                        Z = sl_last[:, None, None] * Tn[:, :, None] / Tn[:, None, :]                  # (m, na, nb)
                        spread = (Z.max(axis=0) - Z.min(axis=0)) / np.abs(Z).max(axis=0)
                        for ia in range(len(ok)):
                            spread[ia, ia] = np.inf
                        for ia, ib in zip(*np.where(spread < 1e-7)):
                            found.append((site, k + 1, int(ok[ia]), int(ok[ib])))
                if found:
                    ncov += 1
                    for p in found:
                        fam[p] += 1
                else:
                    unc.append(mu)
        out[lam]['last_component'] = dict(nsec=nsec, covered=ncov, families=fam.most_common(6))
        if verbose:
            print('     %s: %d sectors with m >= 2: an exact relation of this type exists for the last component in %d (%d without); most frequent (site, k, a, b): %s' %
                  (str(list(lam)), nsec, ncov, len(unc), fam.most_common(5)))
    return out


# ======================================================================================== 10. J6: complexity
def float_physical_V(ch, m):
    """float64 version of physical_V ('auto'): all eigenvalues of the Jacobi matrix of chain 1 on W1 if w1 = m, else of chain 2 on W2 (scipy LAPACK stebz on the symmetrised tridiagonal matrix; top m of W1 would also do)"""
    from scipy.linalg import eigvalsh_tridiagonal
    if ch.w1 != m:
        ch = chain(ch.lam, ch.mu, 2)
    al, uu, dd = ch.f_coeffs()
    n = ch.w1
    if n == 1:
        return [al[0]]
    off = [math.sqrt(uu[k] * dd[k]) for k in range(n - 1)]
    return list(eigvalsh_tridiagonal(np.array(al[:n]), np.array(off), select='i', select_range=(n - m, n - 1)))


def physical_V_fast(lam, mu, dps=30, ch=None, m=None):
    """high-precision Q-free V-spectrum without an mpmath eigen-solver: float64 eigenvalues of the m x m Jacobi matrix (chain 1 on W1 or chain 2 on W2, scipy stebz) polished by Newton steps on the exact
       continuant (`refine_V`); the roots are simple and separated, so the polish converges to the right root.  Ascending list of mpf at `dps` digits."""
    ch = ch or chain(lam, mu)
    m = rs_mult(lam, mu) if m is None else m
    use = ch if ch.w1 == m else chain(lam, mu, 2)
    return [refine_V(lam, mu, v, dps, use) for v in float_physical_V(ch, m)]


def stage_bench(alist=None, verbose=True, with_R=True, Rreps=3):
    se = _se()
    alist = alist or [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 20, 30]
    rows = []
    if verbose:
        print('  %3s %6s %6s %6s %5s | %8s %8s | %8s %8s %8s | %8s %8s | %9s %9s %9s | %8s %8s | %s' %
              ('a', 'dimV', 'irreps', 'copies', 'max w', 't_irr ms', 't_chn ms', 't_eigF ms', 't_recF ms', 'us/cp F', 't_rec30', 'us/cp 30', 't_float', 't_Qfree30', 't_fast30', 'R t_a s', 'R t_b s', 'float64 dev vs mp: forward / twisted (us/copy twisted)'))
    for a in alist:
        lam = (2 * a, a, 0)
        # best of 5 for the pure closed-form parts
        best = None
        for rep in range(5):
            T = {}
            _CH.clear()
            _OFF.clear()
            t0 = time.perf_counter()
            irr = irreps(lam)
            T['irr'] = time.perf_counter() - t0
            t0 = time.perf_counter()
            chains = {mu: chain(lam, mu) for mu in irr}
            T['chain'] = time.perf_counter() - t0
            for ch in chains.values():
                ch.f_coeffs()
            t0 = time.perf_counter()
            evs = {mu: float_physical_V(chains[mu], m) for mu, m in irr.items()}
            T['eig_f'] = time.perf_counter() - t0
            t0 = time.perf_counter()
            slf = {mu: [slice_from_V_float(chains[mu], V) for V in evs[mu]] for mu in irr}
            T['rec_f'] = time.perf_counter() - t0
            t0 = time.perf_counter()
            slt = {mu: [slice_from_V_float(chains[mu], V, twisted=True) for V in evs[mu]] for mu in irr}
            T['rec_t'] = time.perf_counter() - t0
            if best is None or sum(T.values()) < sum(best.values()):
                best = T
        T = best
        ncop = sum(irr.values())
        # mp: recurrence only (V given, 30 digits), then the whole Q-free pipeline at 30 digits
        for ch in chains.values():
            ch.mp_coeffs(30)
        t0 = time.perf_counter()
        phys = {mu: physical_V(lam, mu, 30, chains[mu], m) for mu, m in irr.items()}
        T['eig_mp'] = time.perf_counter() - t0
        t0 = time.perf_counter()
        sl30 = {mu: [slice_from_V(lam, mu, V, 30, chains[mu]) for V in phys[mu]] for mu in irr}
        T['rec30'] = time.perf_counter() - t0
        # fast high-precision Q-free path: float eigenvalues + Newton polish, then the mp recurrence
        t0 = time.perf_counter()
        fastV = {mu: physical_V_fast(lam, mu, 30, chains[mu], m) for mu, m in irr.items()}
        T['eig_fast'] = time.perf_counter() - t0
        fdevV = max([abs(a - b) / max(1, abs(b)) for mu in irr for a, b in zip(fastV[mu], phys[mu])] + [0])
        # accuracy of the float pipeline vs the mp pipeline
        dev = devt = 0.0
        for mu in irr:
            for pf, pt, pm in zip(slf[mu], slt[mu], sl30[mu]):
                sc = max(abs(float(x)) for x in pm)
                dev = max(dev, max(abs(float(x) - y) for x, y in zip(pm, pf)) / sc)
                devt = max(devt, max(abs(float(x) - y) for x, y in zip(pm, pt)) / sc)
        ta = tb = float('nan')
        if with_R and a in (1, 2, 4, 7, 10, 14):
            Tr = se.bench_one(a, Rreps)
            ta, tb = Tr['ta'], Tr['tb']
        maxw = max(ch.w for ch in chains.values())
        t_float = T['irr'] + T['chain'] + T['eig_f'] + T['rec_f']
        row = dict(a=a, dim=dim_gl3(lam), nirr=len(irr), ncop=ncop, maxw=maxw, t_irr=T['irr'], t_chain=T['chain'], t_eig_f=T['eig_f'], t_rec_f=T['rec_f'], t_rec30=T['rec30'], t_eig_mp=T['eig_mp'],
                   t_float=t_float, t_q30=T['irr'] + T['chain'] + T['eig_mp'] + T['rec30'], t_fast30=T['irr'] + T['chain'] + T['eig_fast'] + T['rec30'], fdevV=float(fdevV), fdev=dev, fdev_t=devt, t_rec_t=T['rec_t'], R_ta=ta, R_tb=tb)
        rows.append(row)
        if verbose:
            print('  %3d %6d %6d %6d %5d | %8.2f %8.2f | %8.2f %8.2f %8.2f | %8.2f %8.1f | %9.2f %9.1f %9.1f | %8.3f %8.3f | %.1e' %
                  (a, dim_gl3(lam), len(irr), ncop, maxw, 1e3 * T['irr'], 1e3 * T['chain'], 1e3 * T['eig_f'], 1e3 * T['rec_f'], 1e6 * T['rec_f'] / ncop, 1e3 * T['rec30'], 1e6 * T['rec30'] / ncop,
                   1e3 * t_float, 1e3 * row['t_q30'], 1e3 * row['t_fast30'], ta, tb, dev), end='', flush=True)
            print(' / %.1e (%.2f)' % (devt, 1e6 * T['rec_t'] / ncop), flush=True)
    return rows


def system_polys(lam, mu):
    """S's three echelon Baxter systems of the sector, one by one: {j: monic squarefree polynomial in V (ascending Fractions)} for j = 1, 2, 3 (consistency polynomial of the j-th system alone,
       b -> V by b = 2q - 2V - eta).  Slow (exact determinants): S's pencil costs ~1 s per system per sector at a = 7."""
    import bethe_series_scale as bs
    eta, q = bs.eta_of(lam, mu), bs.q_of(mu)
    out = {}
    for j in (1, 2, 3):
        g = bs.pencil_poly(lam, mu, which=(j,))
        gs = sp.Poly(sp.sqf_part(g.as_expr()), bs._BSYM, domain='QQ')
        out[j] = poly_of_V_from_b([str(c) for c in gs.all_coeffs()], str(q), str(eta))
    return out


def _pgcd(a, b):
    """monic gcd of two Fraction polynomials (ascending)"""
    a, b = list(a), list(b)
    while len(b) > 1 or (len(b) == 1 and b[0] != 0):
        _, r = poly_divmod(a, b)
        a, b = b, r
        while len(b) > 1 and b[-1] == 0:
            b.pop()
        if len(b) == 1 and b[0] == 0:
            break
    return [x / a[-1] for x in a]


def stage_systems(reps=None, verbose=True, special_only_for=((14, 7, 0),)):
    """which information sits in which of S's three Baxter systems?  Per sector: degrees of the single-system polynomials P1, P2, P3 (in V), of gcd(P1,P2), and whether gcd(P1,P2) is the continuant of
       V_red|W1 (the full chain spectrum, physical + non-physical) and gcd(P1,P2,P3) the physical polynomial."""
    reps = reps or [(4, 2, 0), (6, 3, 0), (8, 4, 0), (14, 7, 0)]
    for lam in reps:
        t0 = time.perf_counter()
        irr = irreps(lam)
        n = Counter()
        rows = []
        cache_fn = _cache_name('systems', lam)
        cache = pickle.load(open(cache_fn, 'rb')) if os.path.exists(cache_fn) else {}
        for mu, m in irr.items():
            ch = chain(lam, mu)
            typ = 'E12-cut' if ch.w > ch.w1 else ('E23-spectral' if ch.w1 > m else 'free')
            if lam in special_only_for and typ == 'free' and mu[2] > 0 and mu[0] != 2 * lam[0]:
                continue                                                   # at [14,7,0] only the special sectors (+ the M1 = 0 ones, mu1 = 2 l1, + mu3 = 0) are recomputed
            if mu not in cache:
                cache[mu] = {j: [str(c) for c in p] for j, p in system_polys(lam, mu).items()}
                pickle.dump(cache, open(cache_fn, 'wb'))
            P = {j: [Fr(c) for c in p] for j, p in cache[mu].items()}
            C1 = ch.contin()
            C1m = [x / C1[-1] for x in C1]
            g12 = _pgcd(P[1], P[2])
            g13 = _pgcd(P[1], P[3])
            g23 = _pgcd(P[2], P[3])
            g123 = _pgcd(g12, P[3])
            n['sectors'] += 1
            n['g123_deg_eq_m'] += (len(g123) - 1 == m)
            n['g12_eq_W1_cont'] += (g12 == C1m)
            n['g12_deg_eq_w1'] += (len(g12) - 1 == ch.w1)
            n['g12_deg_gt_m'] += (len(g12) - 1 > m)
            n['P3_deg_eq_m'] += (len(P[3]) - 1 == m)
            n['P1_deg_gt_m'] += (len(P[1]) - 1 > m)
            n['P2_deg_gt_m'] += (len(P[2]) - 1 > m)
            n['typ:' + typ] += 1
            if typ == 'E23-spectral':
                n['e23_total'] += 1
                n['e23_P2_eq_m'] += (len(P[2]) - 1 == m)
                n['e23_g12_eq_m'] += (len(g12) - 1 == m)
                n['e23_g13_gt_m'] += (len(g13) - 1 > m)
            if typ == 'E12-cut':
                n['e12_total'] += 1
                n['e12_P3_eq_m'] += (len(P[3]) - 1 == m)
                n['e12_g13_eq_m'] += (len(g13) - 1 == m)
                n['e12_g12_gt_m'] += (len(g12) - 1 > m)
            rows.append((mu, m, ch.w, ch.w1, typ, len(P[1]) - 1, len(P[2]) - 1, len(P[3]) - 1, len(g12) - 1, len(g13) - 1, len(g23) - 1, len(g123) - 1, g12 == C1m))
        if verbose:
            print('  %s: %d sectors checked: gcd of all three systems has degree m in %d; gcd(P1,P2) is exactly the continuant of V_red|W1 in %d (degree w1 in %d, degree > m in %d); P3 alone has degree m in %d; '
                  'P1 alone has degree > m in %d, P2 alone in %d  [%.1f s]' % (str(list(lam)), n['sectors'], n['g123_deg_eq_m'], n['g12_eq_W1_cont'], n['g12_deg_eq_w1'], n['g12_deg_gt_m'], n['P3_deg_eq_m'],
                                                                                 n['P1_deg_gt_m'], n['P2_deg_gt_m'], time.perf_counter() - t0))
            print("       E23-spectral sectors (%d): the system P2 ALONE has degree m in %d, gcd(P1,P2) has degree m in %d, gcd(P1,P3) has degree > m in %d;  E12-cut sectors (%d): the system P3 ALONE has degree m in %d, "
                  "gcd(P1,P3) has degree m in %d, gcd(P1,P2) has degree > m in %d" % (n['e23_total'], n['e23_P2_eq_m'], n['e23_g12_eq_m'], n['e23_g13_gt_m'], n['e12_total'], n['e12_P3_eq_m'], n['e12_g13_eq_m'], n['e12_g12_gt_m']))
            for r in rows:
                if r[4] != 'free' and (lam == (14, 7, 0) or len(rows) < 12):
                    print('       mu %-12s m %d w %d w1 %d %-12s deg P1 %2d P2 %2d P3 %2d | gcd12 %2d gcd13 %2d gcd23 %2d gcd123 %2d | gcd12 == continuant(W1): %s' %
                          ((str(list(r[0])),) + r[1:]))
    return None


# ======================================================================================== 11. the second chain (E23 constraint as a coordinate cut)
def vred2_apply(rep, v):
    """V'_red = U V_red U^-1 = 1/2 sum sgn(s) sgn(t) E_{s1 t1} E_{s2 t2} Phi'_{s3 t3}, Phi'_cc = l_{4-c}, Phi'_ab = -E_ab (a > b), 0 (a < b)  (slice_engine generators; the weight-space image of V_red under w0)"""
    se = _se()
    lam = rep.lam
    phi = {}
    out = {}
    for coef, ij1, ij2, ij3 in se.VTERMS:
        a, b = ij3
        if a < b:
            continue
        if ij3 not in phi:
            if a == b:
                phi[ij3] = {c: x * lam[3 - a] for c, x in v.items() if lam[3 - a] != 0}
            else:
                phi[ij3] = {c: -x for c, x in rep.apply(ij3, v).items()}
        ph = phi[ij3]
        if not ph:
            continue
        se._vadd(out, rep.apply(ij1, rep.apply(ij2, ph)), coef)
    return {r: x for r, x in out.items() if x != 0}


def exact_full_matrix2(R, lam, mu):
    """V'_red on W' = V[(e3,e2,e1)] from the exact GT generators (independent of the closed forms of chain 2)"""
    W = R.wsp.get((mu[2] - lam[2], mu[1] - lam[1], mu[0] - lam[0]))
    if W is None:
        return [], []
    pos = {int(c): a for a, c in enumerate(W)}
    w = len(W)
    M = [[Fr(0)] * w for _ in range(w)]
    for j, c in enumerate(W):
        for r, x in vred2_apply(R, {int(c): Fr(1)}).items():
            M[pos[r]][j] = x
    return [int(c) for c in W], M


def analyze_chain2(lam, Rfile=None, Rmem=None, verbose=False):
    """chain 2 for every irrep: closed forms vs the exact V' matrix, invariant segment W2 (m1 <= mu2 - l3), positivity of products on W2, w2 vs m, and the physical polynomial = continuant of V'|W2
       (exact against R's CHARPOLY when a file exists, else against the 50-digit eigenvalues): per type (free / E12-cut / E23-spectral)."""
    se = _se()
    lam = tuple(lam)
    R = se.Rep.of(lam)
    C = Counter()
    for mu, m in irreps(lam).items():
        ch1, ch2 = chain(lam, mu, 1), chain(lam, mu, 2)
        typ = 'E12-cut' if ch1.w > ch1.w1 else ('E23-spectral' if ch1.w1 > m else 'free')
        C['sectors'] += 1
        C['type:' + typ] += 1
        W, M = exact_full_matrix2(R, lam, mu)
        C['closed2_ok'] += (W == ch2.idx and M == ch2.matrix())
        band = max([abs(i - j) for i in range(len(W)) for j in range(len(W)) if M[i][j] != 0] + [0])
        C['tridiag2'] += (band <= 1)
        w2 = ch2.w1
        C['w2_eq_m:' + typ] += (w2 == m)
        C['w2_ge_m'] += (w2 >= m)
        C['cut2_u_zero'] += (ch2.w1 == ch2.w or ch2.u[ch2.w1 - 1] == 0)
        C['prod2_pos'] += all(p > 0 for p in ch2.products(w2))
        C['u2d2_neg'] += all(x < 0 for x in ch2.u[:max(w2 - 1, 0)]) and all(x < 0 for x in ch2.d[:max(w2 - 1, 0)])
        C['w2_lt_w:' + typ] += (ch2.w1 < ch2.w)
        if Rfile is not None and mu in Rfile:
            P = Rfile[mu]['charpoly']
            C2 = ch2.contin()
            if typ in ('free', 'E23-spectral'):
                C['charpoly2_eq:' + typ] += (C2 == P)
                C['charpoly2_chk:' + typ] += 1
            else:
                qq, rr = poly_divmod(C2, P)
                C['charpoly2_div:E12-cut'] += all(r == 0 for r in rr)
                C['charpoly2_chk:E12-cut'] += 1
            if typ == 'free':
                C['cont1_eq_cont2:free'] += (ch1.contin() == ch2.contin())
        elif Rmem is not None and mu in Rmem:
            ev = Rmem[mu]['ev']
            if typ in ('free', 'E23-spectral'):
                allev = jacobi_eigs(lam, mu, 50, ch2)
                dev = max([abs(a - b) for a, b in zip(allev, ev)] + [mp.mpf(0)]) if len(allev) == len(ev) else mp.mpf(1)
                C['eig2_eq:' + typ] += (dev < mp.mpf(10) ** (-35))
                C['eig2_chk:' + typ] += 1
    return C


def stage_chain2(reps=None, verbose=True):
    reps = reps or (MAIN_REPS + [(28, 14, 0)])
    for lam in reps:
        t0 = time.perf_counter()
        Rfile = load_R_file(lam) if lam in MAIN_REPS else None
        Rmem = None if Rfile is not None else load_R_mem(lam, 50)
        C = analyze_chain2(lam, Rfile, Rmem)
        if verbose:
            print('  %s: %d irreps: closed-form chain 2 == exact V\' matrix (R generators) in %d, tridiagonal %d, u\'(cut2) = 0 in %d, products > 0 on W2 in %d, u\',d\' < 0 on W2 in %d; w2 >= m in %d; w2 = m: free %d/%d, E23-spectral %d/%d, E12-cut %d/%d (there w2 > m is allowed)' %
                  (str(list(lam)), C['sectors'], C['closed2_ok'], C['tridiag2'], C['cut2_u_zero'], C['prod2_pos'], C['u2d2_neg'], C['w2_ge_m'], C['w2_eq_m:free'], C['type:free'],
                   C['w2_eq_m:E23-spectral'], C['type:E23-spectral'], C['w2_eq_m:E12-cut'], C['type:E12-cut']))
            if Rfile is not None:
                print('       physical polynomial (R CHARPOLY, exact) == continuant of V\'|W2: free %d/%d, E23-spectral %d/%d; in E12-cut sectors CHARPOLY divides it: %d/%d; free sectors: continuant(chain 1, W1) == continuant(chain 2, W2): %d/%d' %
                      (C['charpoly2_eq:free'], C['charpoly2_chk:free'], C['charpoly2_eq:E23-spectral'], C['charpoly2_chk:E23-spectral'], C['charpoly2_div:E12-cut'], C['charpoly2_chk:E12-cut'],
                       C['cont1_eq_cont2:free'], C['type:free']))
            else:
                print('       eigenvalues of V\'|W2 == R 50-digit eigenvalues: free %d/%d, E23-spectral %d/%d  [%.1fs]' % (C['eig2_eq:free'], C['eig2_chk:free'], C['eig2_eq:E23-spectral'], C['eig2_chk:E23-spectral'], time.perf_counter() - t0))
    return None


# ======================================================================================== 12. direct PRV test of the recurrence slices (independent of R's eigenvectors)
def stage_prv_check(lmax=8, dps=50, verbose=True):
    """For every rep [l1,l2,l3] with l1 <= lmax and EVERY copy: take V from the Q-free rule (physical_V), build the slice with slice_from_V (dps digits), convert the components to exact Fractions and apply the
       PRV constraints E12^(l1-l2+1) v and E23^(l2-l3+1) v with R's exact GT generators: the residual (relative to max|v|) must vanish up to the working precision.  Also: the V_red eigen-equation on W
       (exact V_red applied to v) and a NEGATIVE control: the same test with V shifted by 1/3 of the nearest eigenvalue gap must fail."""
    se = _se()
    t0 = time.perf_counter()
    worst_prv = worst_eig = mp.mpf(0)
    ncopies = nbad = nneg = nnegfail = 0
    nreps = 0
    for l1 in range(1, lmax + 1):
        for l2 in range(0, l1 + 1):
            for l3 in range(0, l2 + 1):
                lam = (l1, l2, l3)
                nreps += 1
                R = se.Rep.of(lam)
                a1, b1 = l1 - l2 + 1, l2 - l3 + 1
                for mu, m in irreps(lam).items():
                    ch = chain(lam, mu)
                    Vs = physical_V(lam, mu, dps, ch, m)
                    for V in Vs:
                        p = slice_from_V(lam, mu, V, dps, ch)
                        with mp.workdps(dps + 10):
                            vec = {ch.idx[k]: Fr(int(mp.floor(p[k] * mp.mpf(10) ** (dps + 5))), 10 ** (dps + 5)) for k in range(ch.w) if p[k] != 0}
                            scale = max(abs(x) for x in p)
                            r_prv = mp.mpf(0)
                            for ij, kk in (((1, 2), a1), ((2, 3), b1)):
                                x = dict(vec)
                                for _ in range(kk):
                                    x = R.apply(ij, x)
                                    if not x:
                                        break
                                if x:
                                    r_prv = max(r_prv, max(abs(mp.mpf(c.numerator) / c.denominator) for c in x.values()))
                            # eigen-equation: V_red v - V v = 0 (exact V_red on the rounded vector)
                            img = se.vred_apply(R, vec)
                            r_eig = mp.mpf(0)
                            for k in range(ch.w):
                                i = ch.idx[k]
                                val = mp.mpf(img.get(i, Fr(0)).numerator) / img.get(i, Fr(1)).denominator - mp.mpf(V) * (mp.mpf(vec.get(i, Fr(0)).numerator) / vec.get(i, Fr(1)).denominator)
                                r_eig = max(r_eig, abs(val))
                            worst_prv = max(worst_prv, r_prv / scale)
                            worst_eig = max(worst_eig, r_eig / scale)
                            ncopies += 1
                            nbad += (r_prv / scale > mp.mpf(10) ** (-30)) or (r_eig / scale > mp.mpf(10) ** (-25))
                    # negative control for sectors with w1 > m (E23-spectral): a NON-physical eigenvalue of V_red|W1 must violate the PRV constraint
                    if ch.w1 > m:
                        allev = jacobi_eigs(lam, mu, dps, ch)
                        phys = set(float(v) for v in Vs)
                        for V in allev:
                            if float(V) in phys:
                                continue
                            p = slice_from_V(lam, mu, V, dps, ch)
                            with mp.workdps(dps + 10):
                                vec = {ch.idx[k]: Fr(int(mp.floor(p[k] * mp.mpf(10) ** (dps + 5))), 10 ** (dps + 5)) for k in range(ch.w) if p[k] != 0}
                                scale = max(abs(x) for x in p)
                                x = dict(vec)
                                for _ in range(b1):
                                    x = R.apply((2, 3), x)
                                    if not x:
                                        break
                                r = max([abs(mp.mpf(c.numerator) / c.denominator) for c in x.values()] + [mp.mpf(0)]) / scale
                            nneg += 1
                            nnegfail += (r > mp.mpf(10) ** (-10))
    if verbose:
        print("  PRV test of the recurrence slices at V = physical eigenvalues (Q-free rule), all reps with l1 <= %d (%d reps), %d copies, %d digits: worst relative residual of E12^(l1-l2+1) v, E23^(l2-l3+1) v: %s, "
              "of the exact eigen-equation V_red v = V v: %s, failures: %d;  negative control: %d non-physical eigenvalues of V_red|W1 (w1 > m sectors), %d of them violate E23^(l2-l3+1) v = 0 (> 1e-10)  [%.1fs]" %
              (lmax, nreps, ncopies, dps, mp.nstr(worst_prv, 3), mp.nstr(worst_eig, 3), nbad, nneg, nnegfail, time.perf_counter() - t0))
    return ncopies, nbad, nneg, nnegfail


# ======================================================================================== 13. closed-form Q-system data from the chain (V = alpha - eta/2 for m = 1, sum rule over the copies)
def stage_qdata(lmax=15, dps=50, verbose=True):
    """Q0's closed forms re-derived as consequences of the chains (lambda3 = 0 reps, l1 <= lmax): for m = 1 the V of the sector equals alpha - eta/2 (alpha = (mu1+2)(mu2+1)mu3/2, eta = |L|^2 + cas2(L) - cas2(mu)/2);
       for every m the MEAN of V over the copies equals alpha - eta/2 (Q0: 'mean of a1 over copies = alpha, numerical identity, no proof'; here it is the trace of the explicit m x m Jacobi matrix, a finite sum
       identity in the closed forms), equivalently sum_c b1^(c) = 0 (b1 = 2(alpha - a1), a1 = V + eta/2)."""
    def cas2(mu):
        return sum(mu[i] * (mu[i] + 4 - 2 * (i + 1)) for i in range(3))
    mp.mp.dps = dps
    n1 = b1 = nm = bm = 0
    t0 = time.perf_counter()
    for l1 in range(1, lmax + 1):
        for l2 in range(0, l1 + 1):
            lam = (l1, l2, 0)
            for mu, m in irreps(lam).items():
                alpha = Fr((mu[0] + 2) * (mu[1] + 1) * mu[2], 2)
                eta = Fr(sum(lam) ** 2) + Fr(cas2(lam)) - Fr(cas2(mu), 2)
                tgt = mp.mpf((alpha - eta / 2).numerator) / (alpha - eta / 2).denominator
                Vs = physical_V(lam, mu, dps, m=m)
                if m == 1:
                    n1 += 1
                    b1 += (abs(Vs[0] - tgt) > mp.mpf(10) ** (-(dps - 10)))
                nm += 1
                bm += (abs(mp.fsum(Vs) / m - tgt) > mp.mpf(10) ** (-(dps - 10)))
    if verbose:
        print('  lambda3 = 0 reps with l1 <= %d (%.1fs): m = 1 sectors %d: V = alpha - eta/2 violated in %d; all %d sectors: mean of V over the copies = alpha - eta/2 violated in %d' % (lmax, time.perf_counter() - t0, n1, b1, nm, bm))
    return n1, b1, nm, bm


# @@STAGES@@


# ======================================================================================== CLI
def main(argv):
    cmd = argv[1] if len(argv) > 1 else 'all'
    if cmd in ('derive', 'all'):
        print('== J1/J2 structural argument and symbolic derivation of the closed forms')
        check_index_structure()
        check_index_structure(second=True)
        derive_symbolic()
        derive_symbolic(second=True)
    if cmd in ('closed', 'all'):
        print('== J1/J2: bandwidth of V_red (R operator) and closed-form chain vs R exact matrices, all dominant candidates')
        stage_closed()
    if cmd == 'tables':
        lam = tuple(int(x) for x in argv[2:5]) if len(argv) > 4 else (int(argv[2]), int(argv[3]), 0)
        print('== chain tables of %s' % str(list(lam)))
        stage_tables(lam)
    if cmd == 'all':
        print('== J1: exact chain tables, all irreps of [4,2,0] (alpha_k, u_k = M[k+1,k], d_k = M[k,k+1], u_k d_k; positions k = 0..w-1 ascending m1)')
        stage_tables((4, 2, 0))
        print('== J1: exact chain tables, selected irreps of [14,7,0] (free: (21,14,7), (14,14,14), (28,14,0), (24,12,6); E12-cut: (20,15,7), (19,17,6); E23-spectral: (22,10,10), (23,10,9))')
        stage_tables((14, 7, 0), sectors={(21, 14, 7), (14, 14, 14), (28, 14, 0), (24, 12, 6), (20, 15, 7), (19, 17, 6), (22, 10, 10), (23, 10, 9)})
        print('== J3: exact tables of the SECOND chain (W2 of m1 <= mu2 - l3), E23-spectral irreps: [4,2,0] (6,3,3); [14,7,0] (22,10,10), (23,10,9); and a free one (21,14,7) (cont2 = cont1 there)')
        stage_tables((4, 2, 0), sectors={(6, 3, 3)}, kind=2)
        stage_tables((14, 7, 0), sectors={(22, 10, 10), (23, 10, 9), (21, 14, 7)}, kind=2)
    if cmd in ('struct', 'all'):
        print('== J3: the PRV subspace inside the chain')
        stage_struct()
    if cmd in ('chain2', 'all'):
        print("== J3: the second chain (V' = U V_red U^-1; the E23 constraint is a coordinate cut there)")
        stage_chain2()
    if cmd in ('scan', 'all'):
        print('== J3/J1 broad scan over all small reps')
        stage_scan(int(argv[2]) if (cmd == 'scan' and len(argv) > 2) else 16)
    if cmd in ('dump_tables', 'all'):
        print('== exact chain tables of all irreps of [4,2,0], [6,3,0], [8,4,0], [14,7,0] written to jacobi_cg_tables.txt')
        print('  wrote', write_tables())
        ok, tot = verify_tables()
        print('  re-parsed from the file: continuant of the stored chain (chain 1 on W1 if w1 = m, else chain 2 on W2) == R exact CHARPOLY in %d of %d sectors' % (ok, tot))
    if cmd == 'spencil_sample':
        gen_S_pencil_sample()
    if cmd in ('prv', 'all'):
        print("== J4: direct PRV test of the recurrence slices (exact GT generators), all reps with l1 <= lmax, V from the Q-free rule")
        stage_prv_check(int(argv[2]) if (cmd == 'prv' and len(argv) > 2) else 10)
    if cmd in ('recur', 'all'):
        print("== J4: slice_from_V (three-term recurrence) vs R's slices, V from every source")
        stage_recur()
    if cmd in ('qside', 'all'):
        print("== J5: the Q side (S's Baxter pencil, Q0's Q polynomials)")
        stage_qside()
    if cmd in ('qdata', 'all'):
        print("== J5: Q0's closed forms (V = alpha - eta/2 for m = 1; mean of V over the copies) as consequences of the chains")
        stage_qdata()
    if cmd in ('systems', 'all'):
        print("== J5/J3: which of S's three echelon Baxter systems carries what (single-system polynomials and their gcds)")
        stage_systems(reps=[(4, 2, 0), (6, 3, 0), (8, 4, 0), (14, 7, 0)] if cmd == 'all' or len(argv) < 3 else [tuple(int(x) for x in argv[2:5])])
    if cmd in ('bench', 'all'):
        print('== J6: complexity (Q-free chain pipeline and recurrence vs R slice engine; Lambda = [2a,a,0])')
        stage_bench()


if __name__ == '__main__':
    main(sys.argv)
