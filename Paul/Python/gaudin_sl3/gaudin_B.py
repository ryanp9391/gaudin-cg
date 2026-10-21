"""D3: the Gaudin B operator as the h -> 0 expansion of the XXX SoV operator B(u), the jet basis, the Gaudin wavefunction rule,
the Gaudin K tensor, and the kernel G of round 2 from first principles.   (sl(3), L = 2, Lambda = (2,1,0), theta = (2/7,-3/5))

Run top to bottom:   python -X utf8 gaudin_B.py [--stages B0,B1,B2,B2b,B2c,B3,B3c,B3d,B3e,B4,B4b,B4c,B4d,B4e,B5,B1b]   (caches gaudin_B_*.pkl; delete to recompute;
  about 6 min with warm caches, about 15 min cold; GB_THETA="1/5,-4/9" repeats everything at another theta)

Conventions (new conventions; xxx_chain.py): Lax L_ij(x) = x delta_ij + h E_ij, T(u) = L^(1)(u-theta_1) L^(2)(u-theta_2) UNTWISTED,
  B(u) = T11(u) qmin_{12,12}(u+h) + T21(u) qmin_{12,13}(u+h)   (ref-L2-rotation-sov-all.txt cell 16, Bop2).
Basis order: xxx_chain / notebook order (descending GT patterns, HW ket = index 0) for everything XXX / covector side.  The Gaudin
  eigenstates of P1 (gaudin_chain/qside) are in gl3_reps order (HW index 7); NB2G maps notebook index a -> gl3 index NB2G[a].
Exact Fractions throughout (sparse SpOp), mpmath only for the Gaudin states (60-80 digits).
"""
import os
import sys
import time
import pickle
import argparse
import itertools
from fractions import Fraction as Fr
from math import factorial

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))
import xxx_chain as X  # noqa: E402
from xxx_chain import Chain, SpOp, Ser, op_commutator  # noqa: E402

NB2G = [7, 5, 6, 4, 2, 1, 3, 0]          # notebook (descending GT) index a  ->  gl3_reps index
U0, U1 = Fr(3, 7), Fr(-2, 9)
_TH_ENV = os.environ.get("GB_THETA")          # e.g. GB_THETA="1/5,-4/9" runs the whole pipeline at another (theta_1, theta_2); caches get a suffix
if _TH_ENV:
    CH = Chain(thetas=[Fr(x) for x in _TH_ENV.split(",")])
    SUF = "_th" + _TH_ENV.replace("/", "o").replace(",", "_")
    if U0 in CH.thetas:
        U0 = Fr(5, 11)
    if U1 in CH.thetas:
        U1 = Fr(-7, 13)
else:
    CH = Chain()
    SUF = ""
TH = CH.thetas
N64 = CH.dim

# one-site labels (a|b c), roots {-1-a, -b, -c}; lam2 <= a <= lam1, 0 <= c <= lam2, c <= b <= a (Lambda = (2,1,0))
LABELS1 = [(a, b, c) for a in (2, 1) for c in (0, 1) for b in range(c, a + 1)]
# order used in the notebooks (trip list): sorted root triples
TRIP_ORDER = [(2, 2, 1), (2, 1, 1), (1, 1, 1), (2, 2, 0), (2, 1, 0), (2, 0, 0), (1, 1, 0), (1, 0, 0)]
assert sorted(LABELS1) == sorted(TRIP_ORDER)
LABELS1 = TRIP_ORDER
LABELS2 = [(l1, l2) for l1 in LABELS1 for l2 in LABELS1]           # 64 two-site labels
VAC = (2, 2, 1)                                                     # x0: roots {-3,-2,-1} on both sites


def roots(lab):
    a, b, c = lab
    return (-1 - a, -b, -c)


def hdr(t):
    print("\n" + "=" * 110 + "\n" + t + "\n" + "=" * 110, flush=True)


def cache_path(name):
    return os.path.join(HERE, "gaudin_B_%s%s.pkl" % (name, SUF))


def cached(name, fn):
    p = cache_path(name)
    if os.path.exists(p):
        with open(p, "rb") as f:
            return pickle.load(f)
    r = fn()
    with open(p, "wb") as f:
        pickle.dump(r, f)
    return r


# ---------------------------------------------------------------------------------------------- B(u) and exact algebra
def Bser(u):
    """B(u) with entries Fractions / Ser in h (polynomial in h of degree 6)."""
    h = Ser.h()
    Tof = lambda x: CH.mono_untwisted(x, h)
    T = Tof(u)
    q1 = CH.qmin((1, 2), (1, 2), u + h, h, Tof)
    q2 = CH.qmin((1, 2), (1, 3), u + h, h, Tof)
    return T[0][0] * q1 + T[1][0] * q2


def Bcoeffs(u):
    """[B_0 .. B_6] as SpOp (Fraction), B(u) = sum_k h^k B_k."""
    Bs = Bser(u)
    out = []
    for k in range(7):
        r = {}
        for i, row in Bs.r.items():
            d = {j: v.c[k] for j, v in row.items() if v.c[k] != 0}
            if d:
                r[i] = d
        out.append(SpOp(N64, r))
    for i, row in Bs.r.items():
        for v in row.values():
            assert all(x == 0 for x in v.c[7:]), "B has degree > 6 in h"
    return out


def Bat(Bk, h):
    """B(u) at a rational h from its coefficient matrices."""
    tot = None
    for k, M in enumerate(Bk):
        t = M.scale(h ** k)
        tot = t if tot is None else tot + t
    return tot


def nu(u):
    return (u - TH[0]) * (u - TH[1])


def lam_x(l1, l2, u, h):
    v = Fr(1)
    for s, lab in enumerate((l1, l2)):
        for m in roots(lab):
            v *= u - TH[s] - m * h
    return v


def left_null(M, n=N64):
    """exact basis of {v : v M = 0} (M a SpOp): dense Fraction Gauss-Jordan on M^T."""
    A = [[Fr(0)] * n for _ in range(n)]
    for i, row in M.r.items():
        for j, v in row.items():
            A[j][i] = Fr(v)
    piv, r = [], 0
    for c in range(n):
        p = None
        for i in range(r, n):
            if A[i][c] != 0:
                p = i
                break
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        inv = 1 / A[r][c]
        A[r] = [x * inv for x in A[r]]
        for i in range(n):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c)
        r += 1
    free = [c for c in range(n) if c not in piv]
    sols = []
    for f in free:
        v = [Fr(0)] * n
        v[f] = Fr(1)
        for k, c in enumerate(piv):
            v[c] = -A[k][f]
        sols.append(v)
    return sols


def dense_vec_times(v, M):
    """covector v (list) times SpOp M."""
    return M.vec_times(v)


def solve_lin(cols, target):
    """express sparse vector `target` (dict) as a combination of sparse `cols` (dicts): (coeffs|None, nullity, residual size)."""
    basis, pivs, null = [], {}, 0
    for j, v in enumerate(cols):
        v = dict(v)
        comb = {j: Fr(1)}
        while True:
            ks = [k for k in v if k in pivs]
            if not ks:
                break
            k = ks[0]
            bv, bc = basis[pivs[k]]
            f = v[k] / bv[k]
            for kk, x in bv.items():
                nv = v.get(kk, 0) - f * x
                if nv == 0:
                    v.pop(kk, None)
                else:
                    v[kk] = nv
            for kk, x in bc.items():
                comb[kk] = comb.get(kk, 0) - f * x
        if not v:
            null += 1
            continue
        k = next(iter(v))
        pivs[k] = len(basis)
        basis.append((v, comb))
    t, comb = dict(target), {}
    while True:
        ks = [k for k in t if k in pivs]
        if not ks:
            break
        k = ks[0]
        bv, bc = basis[pivs[k]]
        f = t[k] / bv[k]
        for kk, x in bv.items():
            nv = t.get(kk, 0) - f * x
            if nv == 0:
                t.pop(kk, None)
            else:
                t[kk] = nv
        for kk, x in bc.items():
            comb[kk] = comb.get(kk, 0) + f * x
    return (comb if not t else None), null, len(t)


# ---------------------------------------------------------------------------------------------- operator-valued polynomials in h
NH = 9


class OPol:
    """sum_k h^k M_k with M_k SpOp, truncated at degree NH-1"""
    __slots__ = ("c",)

    def __init__(self, c=None):
        self.c = c if c is not None else [None] * NH

    @staticmethod
    def scalar(x):
        o = OPol()
        o.c[0] = SpOp.identity(N64).scale(x)
        return o

    def __add__(self, o):
        out = OPol()
        for k in range(NH):
            a, b = self.c[k], o.c[k]
            out.c[k] = b if a is None else (a if b is None else a + b)
        return out

    def __sub__(self, o):
        return self + o.scale(-1)

    def scale(self, s):
        out = OPol()
        for k in range(NH):
            if self.c[k] is not None:
                out.c[k] = self.c[k].scale(s)
        return out

    def __mul__(self, o):
        out = OPol()
        for i in range(NH):
            if self.c[i] is None:
                continue
            for j in range(NH - i):
                if o.c[j] is None:
                    continue
                t = self.c[i] * o.c[j]
                out.c[i + j] = t if out.c[i + j] is None else out.c[i + j] + t
        return out

    def shift(self, k):
        out = OPol()
        for i in range(NH - k):
            out.c[i + k] = self.c[i]
        return out


def Lder(i, j, d, u):
    """d-th u-derivative of the Gaudin Lax entry L_ij(u) = sum_s E^(s)_ij/(u-theta_s)."""
    tot = None
    for s in range(2):
        f = Fr((-1) ** d * factorial(d)) / (u - TH[s]) ** (d + 1)
        t = CH.E[s][(i, j)].scale(f)
        tot = t if tot is None else tot + t
    return tot


def Cder(i, j, d, u):
    """d-th derivative of C_ij(u) = sum_k E1_ik E2_kj /((u-theta1)(u-theta2)) (the h^2 term of T(u)/nu(u))."""
    cross = SpOp(N64)
    for k in range(1, 4):
        cross = cross + CH.E[0][(i, k)] * CH.E[1][(k, j)]
    g = Fr((-1) ** d * factorial(d)) / (TH[0] - TH[1]) * ((u - TH[0]) ** (-d - 1) - (u - TH[1]) ** (-d - 1))
    return cross.scale(g)


def tmat_poly(i, j, n, u):
    """t_ij(u + n h) = delta_ij + h L_ij(u+nh) + h^2 C_ij(u+nh) as an OPol (Taylor in h, exact once resummed)."""
    o = OPol()
    o.c[0] = SpOp.identity(N64) if i == j else SpOp(N64)
    for d in range(NH - 1):
        t = Lder(i, j, d, u).scale(Fr(n) ** d / factorial(d))
        o.c[1 + d] = t if o.c[1 + d] is None else o.c[1 + d] + t
    for d in range(NH - 2):
        t = Cder(i, j, d, u).scale(Fr(n) ** d / factorial(d))
        o.c[2 + d] = t if o.c[2 + d] is None else o.c[2 + d] + t
    return o


def b_from_lax(u):
    """b_k(u), k = 0..8, from the exact formula
       B/nu^3 = (1 + h s1 + h^2 s2) * { t11(u)[t11(u+h) t22(u) - t21(u+h) t12(u)] + t21(u)[t11(u+h) t23(u) - t21(u+h) t13(u)] },
       t = T/nu = 1 + h L + h^2 C,  s1 = sum_s 1/(u-theta_s), s2 = 1/((u-theta_1)(u-theta_2))   (no reference to Chain.mono)."""
    t = lambda i, j, n: tmat_poly(i, j, n, u)
    br = t(1, 1, 0) * (t(1, 1, 1) * t(2, 2, 0) - t(2, 1, 1) * t(1, 2, 0)) + t(2, 1, 0) * (t(1, 1, 1) * t(2, 3, 0) - t(2, 1, 1) * t(1, 3, 0))
    s1 = 1 / (u - TH[0]) + 1 / (u - TH[1])
    s2 = 1 / ((u - TH[0]) * (u - TH[1]))
    pref = OPol.scalar(1) + OPol.scalar(s1).shift(1) + OPol.scalar(s2).shift(2)
    return pref * br


def eigenclass_blocks(diag):
    out = {}
    for i, v in enumerate(diag):
        out.setdefault(v, []).append(i)
    return out


def charpoly_block(M, idx):
    import sympy
    A = sympy.Matrix([[sympy.Rational(M.r.get(i, {}).get(j, 0).numerator, M.r.get(i, {}).get(j, 0).denominator) if M.r.get(i, {}).get(j, 0) != 0 else 0
                       for j in idx] for i in idx])
    x = sympy.Symbol("x")
    return sympy.factor_list(A.charpoly(x).as_expr())


# ================================================================================================ B0
def stage_B0():
    hdr("B0: port and verify B(u)  (exact; untwisted T, new conventions)")
    print("check of the gl(3) relations on the one-site generators: violations =", CH.check_commutation())
    h = Fr(1, 2)
    Bk0 = Bcoeffs(U0)
    Bk1 = Bcoeffs(U1)
    Bh0, Bh1 = Bat(Bk0, h), Bat(Bk1, h)
    print("B(u) is a polynomial of degree 6 in h; nnz of B_0..B_6 at u0 = 3/7:", [M.nnz() for M in Bk0])
    print("B_0 = nu(u)^3 * 1 :", all(v == nu(U0) ** 3 for v in Bk0[0].flat().values()) and Bk0[0].nnz() == 64)
    comm = op_commutator(Bh0, Bh1)
    print("[B(3/7), B(-2/9)] = 0 at h = 1/2 exactly:", comm.is_zero())
    # degree in u from 7 points + an 8th point
    us = [Fr(1, 7) + Fr(k, 3) for k in range(7)]
    Bs = [Bat(Bcoeffs(u), h) for u in us]
    lead, tot8 = None, None
    u8 = Fr(19, 23)
    for k in range(7):
        den, w = Fr(1), Fr(1)
        for l in range(7):
            if l != k:
                den *= us[k] - us[l]
                w *= (u8 - us[l]) / (us[k] - us[l])
        t = Bs[k].scale(1 / den)
        lead = t if lead is None else lead + t
        t8 = Bs[k].scale(w)
        tot8 = t8 if tot8 is None else tot8 + t8
    print("degree 6 in u: leading coefficient is the identity:", lead.nnz() >= 64 and set(lead.flat().values()) == {Fr(1)} and len(lead.flat()) == 64,
          "; interpolation from 7 points reproduces B(19/23) exactly:", (tot8 - Bat(Bcoeffs(u8), h)).is_zero())
    print("one-site labels (a|bc) -> roots {-1-a,-b,-c}:")
    for lab in LABELS1:
        print("   ", lab, roots(lab))
    # spectrum, nullspaces, u independence
    t0 = time.time()
    covs, ok_sp, ok_u = {}, 0, 0
    evs = {}
    for lab in LABELS2:
        lam0 = lam_x(*lab, U0, h)
        M = Bh0 - SpOp.identity(N64).scale(lam0)
        ns = left_null(M)
        if len(ns) != 1:
            print("  label", lab, "nullity", len(ns))
            continue
        v = ns[0]
        ok_sp += 1
        r = Bh1.vec_times(v)
        lam1 = lam_x(*lab, U1, h)
        if all(r[i] == lam1 * v[i] for i in range(N64)):
            ok_u += 1
        covs[lab] = v
        evs[lab] = lam0
    print("the 64 label-rule eigenvalues each have a 1-dim exact left nullspace of B(3/7)-lambda_x:", ok_sp, "/ 64;",
          "and the covector is an eigencovector of B(-2/9) with the label-rule eigenvalue (u independence):", ok_u, "/ 64   (%.1fs)" % (time.time() - t0))
    vals = list(evs.values())
    print("the 64 eigenvalues are pairwise distinct (simple spectrum):", len(set(vals)) == 64)
    # trace check: sum of eigenvalues = trace B
    trB = sum(Bh0.r.get(i, {}).get(i, 0) for i in range(N64))
    print("trace B(3/7) = sum of the label-rule eigenvalues:", trB == sum(vals))
    print("the 64 pairs (l1,l2) (first site slow): first 8 =", LABELS2[:8], "...")
    return Bk0, Bk1


# ================================================================================================ B1
def stage_B1():
    hdr("B1: the h-expansion B = nu^3 [1 + h b1 + h^2 b2 + ...] in Gaudin Lax terms")
    res = {}
    # b_k(u) matrices
    for u in (U0, U1, Fr(5, 3)):
        Bk = Bcoeffs(u)
        bk = [M.scale(1 / nu(u) ** 3) for M in Bk]
        bl = b_from_lax(u)
        agree = all((bk[k] - bl.c[k]).is_zero() if bl.c[k] is not None else bk[k].is_zero() for k in range(7))
        zero78 = all((bl.c[k] is None) or bl.c[k].is_zero() for k in (7, 8))
        print("u = %s: b_k(u) = Taylor/Lax formula (t = 1 + hL + h^2 C, k = 0..6): %s ; orders 7, 8 of the formula vanish: %s" % (u, agree, zero78))
    u = U0
    Bk0, Bk1 = Bcoeffs(U0), Bcoeffs(U1)
    b0 = [M.scale(1 / nu(U0) ** 3) for M in Bk0]
    b1_ = [M.scale(1 / nu(U1) ** 3) for M in Bk1]
    # explicit b1, b2
    s1 = lambda u: 1 / (u - TH[0]) + 1 / (u - TH[1])
    s2 = lambda u: 1 / ((u - TH[0]) * (u - TH[1]))
    for uu, bk in ((U0, b0), (U1, b1_), (Fr(5, 3), [M.scale(1 / nu(Fr(5, 3)) ** 3) for M in Bcoeffs(Fr(5, 3))])):
        L = lambda i, j, d=0: Lder(i, j, d, uu)
        C = lambda i, j, d=0: Cder(i, j, d, uu)
        b1x = L(1, 1).scale(2) + L(2, 2) + SpOp.identity(N64).scale(s1(uu))
        b2x = (C(1, 1).scale(2) + C(2, 2) + L(1, 1, 1) + L(1, 1) * L(1, 1) + (L(1, 1) * L(2, 2)).scale(2) - L(2, 1) * L(1, 2) + L(2, 1) * L(2, 3)
               + (L(1, 1).scale(2) + L(2, 2)).scale(s1(uu)) + SpOp.identity(N64).scale(s2(uu)))
        print("u = %s: b1 = 2 L11 + L22 + nu'/nu: %s ; b2 = [2C11 + C22 + L11' + L11^2 + 2 L11 L22 - L21 L12 + L21 L23] + (nu'/nu)(2L11+L22) + s2: %s" % (
            uu, (bk[1] - b1x).is_zero(), (bk[2] - b2x).is_zero()))
    # is b2 in the span of {L_ij L_kl, L'_ij, L''_ij, L_ij, 1} (Gaudin Lax entries only)?
    Id = SpOp.identity(N64)
    Ld = {(i, j, d_): Lder(i, j, d_, U0) for i in range(1, 4) for j in range(1, 4) for d_ in range(3)}
    cols_ = [op.flat() for op in Ld.values()] + [Id.flat()]
    base_n = len(cols_)
    cols_ += [(Ld[(i, j, 0)] * Ld[(k, l, 0)]).flat() for i in range(1, 4) for j in range(1, 4) for k in range(1, 4) for l in range(1, 4)]
    co, nul, res_ = solve_lin(cols_, b0[2].flat())
    print("b2(u0) in span{ L_ij L_kl (81 ordered products), L_ij, L'_ij, L''_ij, 1 } :", co is not None, "(residual entries %d)" % res_)
    cross = []
    for i in range(1, 4):
        for j in range(1, 4):
            cr = SpOp(N64)
            for k in range(1, 4):
                cr = cr + CH.E[0][(i, k)] * CH.E[1][(k, j)]
            cross.append(cr.flat())
    co2, nul2, res2 = solve_lin(cols_ + cross, b0[2].flat())
    print("b2(u0) in the span after adding the 9 cross-site entries C_ij = sum_k E1_ik E2_kj (h^2 term of T(u)/nu(u)):", co2 is not None, "(residual entries %d)" % res2)
    # b1 is Cartan: diagonal
    d1 = [b0[1].r.get(i, {}).get(i, Fr(0)) for i in range(N64)]
    print("b1(u0) diagonal in the GT product basis:", b0[1].nnz() == 64 and all(j == i for i, row in b0[1].r.items() for j in row))
    print("b1(u0) = sum_s (2 E^s_11 + E^s_22 + 1)/(u-theta_s): eigenvalue w_s = 1 + a + b + c per site; explicit check of the multiset:")
    pred = sorted(sum(Fr(1 + sum(l), 1) / (U0 - TH[s]) for s, l in enumerate(lab)) for lab in LABELS2)
    print("   multiset of diag(b1) == multiset of sum_s (1+a+b+c)_s/(u0-theta_s) over the 64 labels:", sorted(d1) == pred)
    blocks = eigenclass_blocks(d1)
    print("   b1(u0) eigenspaces (generic u0): number %d, dimensions %s" % (len(blocks), sorted(len(v) for v in blocks.values())))
    from collections import Counter
    print("   dimension tally:", dict(Counter(len(v) for v in blocks.values())))
    # commutators
    hdr("B1: commutators  [b_i(u), b_j(v)],  u = 3/7, v = -2/9")
    zero = lambda M: M.is_zero()
    tab = {}
    for i in range(1, 7):
        for j in range(1, 7):
            tab[(i, j)] = not zero(op_commutator(b0[i], b1_[j]))
    print("nonzero [b_i(u), b_j(v)] (1 = nonzero), rows i, cols j = 1..6:")
    for i in range(1, 7):
        print("   i=%d  " % i + " ".join("1" if tab[(i, j)] else "0" for j in range(1, 7)))
    for n in range(2, 9):
        tot = None
        for i in range(1, n):
            j = n - i
            if i <= 6 and j <= 6:
                t = op_commutator(b0[i], b1_[j])
                tot = t if tot is None else tot + t
        print("   order h^%d: sum_{i+j=%d} [b_i(u), b_j(v)] = 0: %s" % (n, n, zero(tot) if tot is not None else True))
    same = [(i, j) for i in range(1, 5) for j in range(i + 1, 5) if not zero(op_commutator(b0[i], b0[j]))]
    print("same-u commutators [b_i(u), b_j(u)] nonzero for (i,j) in", same)
    print("b1 vs the Lax entries: [b1(u), L_ij(v)] etc. are not zero (b1 Cartan only); [b1(u), b1(v)]:", zero(op_commutator(b0[1], b1_[1])))
    # spectral data of b2, b3 compressed on b1 eigenspaces
    hdr("B1: b2 compressed to the b1-eigenspaces (generic u0 = 3/7) vs the label-rule prediction")
    u = U0
    t1, t2 = 1 / (u - TH[0]), 1 / (u - TH[1])
    def esym(lab):
        m = roots(lab)
        return -sum(m), m[0] * m[1] + m[0] * m[2] + m[1] * m[2], -m[0] * m[1] * m[2]   # e1' = -e1 (so b1 eigenvalue = e1'), e2, e3'
    # label-rule coefficient series of  lambda_x(h)/nu^3 = prod_s (1 - e1 t_s + e2 t_s^2 - e3 t_s^3), t_s = h/(u-theta_s)
    def coef_series(lab):
        pol = [Fr(1)]
        for s, l in enumerate(lab):
            m = roots(l)
            ts = (t1, t2)[s]
            f = [Fr(1)]
            for mk in m:
                f2 = [Fr(0)] * (len(f) + 1)
                for k, c in enumerate(f):
                    f2[k] += c
                    f2[k + 1] += -mk * ts * c
                f = f2
            new = [Fr(0)] * (len(pol) + len(f) - 1)
            for i, a in enumerate(pol):
                for j, b_ in enumerate(f):
                    new[i + j] += a * b_
            pol = new
        return pol  # coefficients of h^0..h^6
    series = {lab: coef_series(lab) for lab in LABELS2}
    ok_all, n_split = True, 0
    for key, idx in sorted(blocks.items(), key=lambda kv: kv[0]):
        labs = [lab for lab in LABELS2 if series[lab][1] == key]
        assert len(labs) == len(idx)
        pred = sorted(series[lab][2] for lab in labs)
        if len(idx) == 1:
            got = [b0[2].r.get(idx[0], {}).get(idx[0], Fr(0))]
        else:
            fl = charpoly_block(b0[2], idx)
            import sympy
            got = []
            for f, mult in fl[1]:
                if sympy.degree(f, sympy.Symbol("x")) != 1:
                    got = None
                    break
                r = sympy.solve(f, sympy.Symbol("x"))[0]
                got += [Fr(int(r.p), int(r.q))] * mult
            if got is not None:
                got.sort()
        # b1 diagonal => the 1x1 case needs b2 diagonal entry; for blocks the compressed b2 eigenvalues
        match = got is not None and got == pred
        ok_all &= match
        n_split += len(set(got)) if got is not None else 0
    print("compressed-b2 eigenvalues on every b1-eigenspace = the h^2 coefficients of the label-rule eigenvalue series  (second-order degenerate PT; all", len(blocks), "eigenspaces):", ok_all)
    print("  number of distinct (b1, compressed-b2) pairs =", n_split, "of 64  (b2 splits every b1-eigenspace completely: %s)" % (n_split == 64))
    hdr("B1: full characteristic polynomials of b_k(u0) (generalised eigenspaces of single b_k) ")
    from sympy.polys.matrices import DomainMatrix
    from sympy import QQ
    import sympy
    xs = sympy.Symbol("x")
    def dm(M):
        rows = [[QQ(M.r.get(i, {}).get(j, Fr(0)).numerator, M.r.get(i, {}).get(j, Fr(0)).denominator) for j in range(N64)] for i in range(N64)]
        return DomainMatrix(rows, (N64, N64), QQ)
    from collections import Counter
    for k in range(1, 7):
        cp = dm(b0[k]).charpoly()
        P = sympy.Poly([QQ.to_sympy(c) for c in cp], xs)
        fl = sympy.factor_list(P.as_expr())
        mult = Counter()
        for f, m in fl[1]:
            mult[(sympy.degree(f, xs), m)] += 1
        print("   b%d(u0): factor types (degree, multiplicity): count" % k, dict(mult))
        if k == 2:
            # Jordan check for repeated roots
            for f, m in fl[1]:
                if m > 1 and sympy.degree(f, xs) == 1:
                    r = sympy.solve(f, xs)[0]
                    lam = Fr(int(r.p), int(r.q))
                    Mx = b0[2] - SpOp.identity(N64).scale(lam)
                    n1 = len(left_null(Mx))
                    n2 = len(left_null(Mx * Mx))
                    print("      repeated root %s of b2: dim ker(b2-lam) = %d, dim ker (b2-lam)^2 = %d" % (lam, n1, n2))
    return b0, b1_



# ================================================================================================ B2 helpers: exact rational reconstruction
def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]


def pmul_(a, b):
    out = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x != 0:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def ptrim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def pdivmod_(a, b):
    a = ptrim(a)
    b = ptrim(b)
    q = [Fr(0)] * max(len(a) - len(b) + 1, 1)
    a = list(a)
    while len(a) >= len(b) and a:
        c = a[-1] / b[-1]
        k = len(a) - len(b)
        q[k] = c
        for i, y in enumerate(b):
            a[i + k] -= c * y
        a = ptrim(a)
    return ptrim(q), a


def peval_(a, x):
    r = Fr(0)
    for c in reversed(a):
        r = r * x + c
    return r


def newton_interp(xs, ys):
    """interpolating polynomial (coefficient list low->high) through the points."""
    n = len(xs)
    cf = list(ys)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            cf[i] = (cf[i] - cf[i - 1]) / (xs[i] - xs[i - j])
    # Newton -> monomial
    poly = [Fr(0)]
    for i in range(n - 1, -1, -1):
        poly = padd(pmul_(poly, [-xs[i], Fr(1)]), [cf[i]])
    return ptrim(poly)


def rat_recon(xs, ys, nmax=None):
    """rational function p/q through the points, minimal degree via the Euclid algorithm on (prod(h-h_i), interpolant);
    returns (p, q) (low->high), q monic-normalised q[0]!=0 scaling free, or None if no (dn+dd < len-2) solution."""
    n = len(xs)
    I = newton_interp(xs, ys)
    M = [Fr(1)]
    for x in xs:
        M = pmul_(M, [-x, Fr(1)])
    r0, r1 = M, I
    t0, t1 = [Fr(0)], [Fr(1)]
    while len(ptrim(r1)) - 1 >= n // 2 + (n % 2) - 1 + 0 and ptrim(r1):
        if len(ptrim(r1)) - 1 < (n - 1) // 2:
            break
        q, r = pdivmod_(r0, r1)
        r0, r1 = r1, r
        t0, t1 = t1, padd(t0, [-c for c in pmul_(q, t1)]) if q else t0
    p, qd = ptrim(r1), ptrim(t1)
    if not qd:
        return None
    # verify
    for x, y in zip(xs, ys):
        if peval_(qd, x) == 0 or peval_(p, x) / peval_(qd, x) != y:
            return None
    if len(p) + len(qd) - 2 >= n - 2:
        return None
    # normalise so that qd(0) = 1 if possible
    c = qd[0] if qd[0] != 0 else qd[-1]
    return [x / c for x in p], [x / c for x in qd]


def null_dense(rows, ncols):
    """right nullspace of a dense Fraction matrix (list of rows)."""
    A = [list(r) for r in rows]
    piv, r = [], 0
    for c in range(ncols):
        p = None
        for i in range(r, len(A)):
            if A[i][c] != 0:
                p = i
                break
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        inv = 1 / A[r][c]
        A[r] = [x * inv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c)
        r += 1
    free = [c for c in range(ncols) if c not in piv]
    out = []
    for f in free:
        v = [Fr(0)] * ncols
        v[f] = Fr(1)
        for k, c in enumerate(piv):
            v[c] = -A[k][f]
        out.append(v)
    return out


RVEC = [Fr(1, i + 2) for i in range(N64)]


def hpoints(n):
    """n rational values of h = 1/m (m skipping multiples of 5 and 7, where the ladders of the two sites collide)"""
    out, m = [], 2
    while len(out) < n:
        if m % 5 and m % 7:
            out.append(Fr(1, m))
        m += 1
    return out


def cov_at(Bk, lab, h):
    M = Bat(Bk, h) - SpOp.identity(N64).scale(lam_x(*lab, U0, h))
    ns = left_null(M)
    assert len(ns) == 1, (lab, h, len(ns))
    v = ns[0]
    n = sum(a * b for a, b in zip(v, RVEC))
    return [x / n for x in v]


def recon_covector(Bk, lab, npts=40):
    """normalised (v.RVEC = 1) eigencovector as a vector of rational functions with common denominator: (numerators, den)."""
    hs = hpoints(npts)
    vals = [cov_at(Bk, lab, h) for h in hs]
    # generic combination to find the common denominator
    w = [Fr((7 * i * i + 3 * i + 5) % 31 + 1) for i in range(N64)]
    ys = [sum(a * b for a, b in zip(w, v)) for v in vals]
    rr = rat_recon(hs, ys)
    if rr is None:
        return None
    p, q = rr
    nums = []
    for c in range(N64):
        yc = [v[c] * peval_(q, h) for v, h in zip(vals, hs)]
        P = newton_interp(hs, yc)
        if len(P) > len(hs) - 4:
            return None
        nums.append(P)
    return nums, q

# ================================================================================================ B2
def get_covs():
    Bk = Bcoeffs(U0)
    def build():
        out = {}
        for lab in LABELS2:
            r = recon_covector(Bk, lab, npts=34)
            assert r is not None, lab
            out[lab] = r
        return out
    return cached("cov", build)


def cov_series(nums, q, K):
    """Taylor coefficients 0..K of v(h) = nums(h)/q(h) (q[0] != 0), vector valued: list of K+1 vectors."""
    qs = list(q) + [Fr(0)] * (K + 1)
    out = []
    # series of 1/q
    inv = [Fr(0)] * (K + 1)
    inv[0] = 1 / qs[0]
    for n in range(1, K + 1):
        inv[n] = -sum(qs[k] * inv[n - k] for k in range(1, n + 1)) / qs[0]
    for k in range(K + 1):
        vec = []
        for c in range(N64):
            P = nums[c]
            tot = Fr(0)
            for i in range(min(len(P), k + 1)):
                tot += P[i] * inv[k - i]
            vec.append(tot)
        out.append(vec)
    return out


def rank_of(rows):
    return len(rows[0]) - len(null_dense(rows, len(rows[0]))) if rows else 0


def proj_key(v):
    for x in v:
        if x != 0:
            return tuple(y / x for y in v)
    return None


def stage_B2():
    hdr("B2: SoV covectors <x|(h) as exact rational functions of h")
    covs = get_covs()
    Bk = Bcoeffs(U0)
    # degrees
    degs = {lab: (len(q) - 1, max(len(P) - 1 for P in nums)) for lab, (nums, q) in covs.items()}
    from collections import Counter
    print("(deg den, max deg num) of the r-normalised covector components, tally over the 64 labels:", dict(Counter(degs.values())))
    print("maximal degrees: den %d, num %d" % (max(d[0] for d in degs.values()), max(d[1] for d in degs.values())))
    bad = [lab for lab, (nums, q) in covs.items() if q[0] != 1]
    print("q(0) = 1 for all labels (the covector is analytic at h = 0 with <x|.r = 1):", not bad)
    # held-out validation: exact left-eigen relation at unseen rational h
    ok = 0
    held = [Fr(3, 11), Fr(4, 13), Fr(9, 17), Fr(2, 19), Fr(-5, 23), Fr(7, 3)]
    for lab, (nums, q) in covs.items():
        good = True
        for h in held:
            qv = peval_(q, h)
            v = [peval_(P, h) / qv for P in nums]
            if sum(a * b for a, b in zip(v, RVEC)) != 1:
                good = False
            r = Bat(Bk, h).vec_times(v)
            lam = lam_x(*lab, U0, h)
            good = good and all(r[i] == lam * v[i] for i in range(N64))
        ok += good
    print("held-out exact validation (6 unseen h, incl. h=7/3 and h<0): v(h)B(h) = lambda_x(h) v(h), v.r = 1:", ok, "/ 64")
    # leading covectors
    K = 12
    ser = {lab: cov_series(nums, q, K) for lab, (nums, q) in covs.items()}
    lead = {lab: ser[lab][0] for lab in LABELS2}
    classes = {}
    for lab in LABELS2:
        classes.setdefault(proj_key(lead[lab]), []).append(lab)
    print("distinct leading covectors <x|_0 (projectively): %d  among 64 labels" % len(classes))
    print("class sizes:", dict(Counter(len(v) for v in classes.values())))
    # filtration ranks (union of Taylor coefficients up to order k)
    for kmax in range(0, 7):
        rows = [ser[lab][k] for lab in LABELS2 for k in range(kmax + 1)]
        print("   rank span{<x|_k : k <= %d} = %d" % (kmax, rank_of(rows)))
    return covs, ser, classes

# ================================================================================================ B2: diagonal frame, lattice, jets
KSER = 40
ZETAS = {'P': (Fr(3, 2), Fr(5, 4), Fr(2, 3)), 'Pp': (Fr(2), Fr(5, 3), Fr(-1, 2)), 'X1': (Fr(7, 5), Fr(-3, 4), Fr(1, 6))}
ZETA_P = ZETAS['P']


def set_kser(k=KSER):
    X.NSER = k


def ser_shift(a, e):
    """h^e * a  (e >= 0)"""
    return Ser([Fr(0)] * e + list(a.c))


def ser_zero():
    return Ser([Fr(0)])


def sp_dense(M):
    """SpOp (Ser/Fraction entries) -> dense list of lists of Ser"""
    n = M.n
    out = [[Ser([Fr(0)]) for _ in range(n)] for _ in range(n)]
    for i, row in M.r.items():
        for j, v in row.items():
            out[i][j] = v if isinstance(v, Ser) else Ser([Fr(v)])
    return out


def pinv_series(zeta):
    """one-site P_phys(z)^{-1} (notebook order, acts on kets as a column matrix), z_j = exp(h zeta_j), as an 8x8 matrix of Ser
       (formula of xxx_limit.mp_x0_onesite: ppinv . dinv . pminv).  Every entry is O(h)."""
    z = [Ser.exp(Fr(x)) for x in zeta]
    z1, z2, z3 = z
    one = Ser.const(1)
    d12, d13, d23 = (z1 - z2).divh(1), (z1 - z3).divh(1), (z2 - z3).divh(1)
    E = CH.E1
    def expN(c, key):
        Ek = E[key]
        return SpOp.identity(8) + Ek.scale(c) + (Ek * Ek).scale(c * c * Fr(1, 2))
    z1i, z2i = z1.inv(), z2.inv()
    pp = expN(-(z1 + z2) * z1i * z2i, (2, 3)) * expN(z1i * z2i, (1, 3)) * expN(-z1i, (1, 2))
    # s1 = h^-2 z1^2/(d12 d13), s2 = h^-1 (-z2/(z1 d23)), s3 = 1/(z1 z2)
    u1 = (d12 * d13) * (z1i * z1i)               # (h^2 s1)^-1
    u2 = (-(z1 * d23)) * z2i                       # (h s2)^-1 = -z1 d23 / z2
    u3 = z1 * z2                                   # s3^-1
    wts = [[CH.E1[(k, k)].r.get(p, {}).get(p, Fr(0)) for k in (1, 2, 3)] for p in range(8)]
    dr = {}
    for p, (w1, w2, w3) in enumerate(wts):
        w1, w2, w3 = int(w1), int(w2), int(w3)
        un = (u1 ** w1) * (u2 ** w2) * (u3 ** w3)
        dr[p] = {p: ser_shift(un, 2 * w1 + w2)}
    dinv = SpOp(8, dr)
    r1 = -(z3 * z2i)
    r2 = -(z3 * z1i)
    r3 = (z2 * z2 * (z3 - z1).divh(1)) * (z1i * z1i) * (z2 - z3).divh(1).inv()
    pm = expN(-r3, (2, 1)) * expN(-r2, (3, 1)) * expN(-r1, (3, 2))
    return pp * dinv * pm


def ser_val(a):
    return a.val()


def diag_rows(covs, P):
    """f^d_x(h) = nums_x(h) (P (x) P) for every label: dict lab -> list of 64 Ser (nums = polynomial numerator of the r-normalised
       covector; the dropped unit 1/q(h) is irrelevant for lattices)."""
    cache = {}
    def outer(i, k):
        if (i, k) not in cache:
            cache[(i, k)] = [P[i][j] * P[k][l] for j in range(8) for l in range(8)]
        return cache[(i, k)]
    rows = {}
    for lab, (nums, q) in covs.items():
        row = [Ser([Fr(0)]) for _ in range(N64)]
        for c, Pc in enumerate(nums):
            if not any(x != 0 for x in Pc):
                continue
            i, k = divmod(c, 8)
            ps = Ser(Pc)
            o = outer(i, k)
            for m in range(N64):
                if o[m].val() is not None:
                    row[m] = row[m] + ps * o[m]
        rows[lab] = row
    return rows


def depth1(lab):
    return 6 - (1 + sum(lab))


def depth2(lab2):
    return depth1(lab2[0]) + depth1(lab2[1])


def row_val(row):
    vs = [x.val() for x in row if x.val() is not None]
    return min(vs) if vs else None


def row_shift(row, v):
    return [x.divh(v) if v > 0 else x for x in row]


def saturate(rows, order, verbose=False):
    """Saturate the lattice generated over Q[[h]] by the given rows (after dividing each by its own h-valuation): echelon reduction on the
       leading (h^0) coefficients with constant multipliers; a row whose reduced leading term vanishes is divided by h.  Rows are processed in `order`
       (this fixes the gauge: each new covector is a divided difference with respect to the earlier ones).  Returns list of dicts."""
    acc = []
    K = X.NSER
    for lab in order:
        row = list(rows[lab])
        v0 = row_val(row)
        assert v0 is not None
        d = v0                       # total number of divisions by h applied to this generator (valuation of the generator)
        row = row_shift(row, v0)
        prec = K - v0
        base_val = v0
        ndiv = 0                      # divisions beyond the normalisation (the jet defect)
        while True:
            lead = [x.c[0] for x in row]
            for a in acc:
                p = a["piv"]
                if lead[p] != 0:
                    c = lead[p]
                    row = [x - a["row"][m] * c for m, x in enumerate(row)]
                    lead = [x.c[0] for x in row]
            if any(x != 0 for x in lead):
                break
            v = row_val(row)
            assert v is not None and v < prec - 1, ("precision exhausted", lab, prec)
            row = row_shift(row, v)
            prec -= v
            ndiv += v
        piv = next(i for i, x in enumerate(lead) if x != 0)
        c = 1 / lead[piv]
        row = [x * c for x in row]
        acc.append(dict(lab=lab, row=row, piv=piv, defect=ndiv, val0=base_val, prec=prec, lead=[x.c[0] for x in row]))
        if verbose:
            print("   %s: valuation %d, jet defect %d, pivot %d, prec %d" % (lab, base_val, ndiv, piv, prec))
    return acc


def stage_B2b(tag="P"):
    hdr("B2 (diagonal frame): P_phys^{-1}(e^{h zeta}) as an exact h-series and the frame change of the SoV covectors")
    set_kser(KSER)
    zeta = ZETAS[tag]
    print("zeta point", tag, zeta)
    covs = get_covs()
    t0 = time.time()
    P = sp_dense(pinv_series(zeta))
    print("P_phys^{-1}(e^{h zeta}) series built (%.1fs); valuation of every entry:" % (time.time() - t0), sorted({x.val() for row in P for x in row if x.val() is not None}))
    z = [Ser.exp(Fr(x)) for x in zeta]
    lin = [z[1] - z[0], z[2] - z[0], z[1] - z[0], 2 * z[2] - z[0] - z[1], 3 * (z[1] - z[0]), 3 * (z[2] - z[0]), 2 * (z[2] - z[1]), 3 * (z[2] - z[1])]
    print("first row e1 P^{-1} equals D2's linear forms (z2-z1, z3-z1, z2-z1, 2z3-z1-z2, 3(z2-z1), 3(z3-z1), 2(z3-z2), 3(z3-z2)):",
          all(P[0][j] == lin[j] for j in range(8)))
    t0 = time.time()
    rows = diag_rows(covs, P)
    print("diagonal-frame covectors f^d_x(h) = f_x(h)(P (x) P) built (%.1fs)" % (time.time() - t0))
    vals = {lab: row_val(r) for lab, r in rows.items()}
    from collections import Counter
    print("h-valuation v_x of f^d_x (vacuum = %d):" % vals[(VAC, VAC)], dict(sorted(Counter(vals.values()).items())))
    print("  v_x - v_vac against the depth D(x) = sum_s (6 - (1+a+b+c)) :", all(vals[lab] - vals[(VAC, VAC)] >= 0 for lab in LABELS2),
          "; equal to depth for all labels:", all(vals[lab] - vals[(VAC, VAC)] == depth2(lab) for lab in LABELS2))
    # leading terms
    leads = {lab: [x.c[0] for x in row_shift(rows[lab], vals[lab])] for lab in LABELS2}
    print("rank of the 64 normalised leading covectors (h^0 term of f^d_x / h^{v_x}):", rank_of([leads[lab] for lab in LABELS2]))
    xg = [z_.c[1] for z_ in lin]
    xG2 = [a * b for a in xg for b in xg]
    print("leading covector of the vacuum label = x^G (x) x^G (D2's x0^G, per site (z2-z1,...)'): ",
          proj_key(leads[(VAC, VAC)]) == proj_key(xG2))
    ok = all(proj_key(leads[lab]) == proj_key(xG2) for lab in LABELS2 if True)
    print("all 64 leading covectors are proportional to x^G (x) x^G (single class):", ok)
    order = sorted(LABELS2, key=lambda l: (depth2(l), l))
    t0 = time.time()
    def build_acc():
        a_ = saturate(rows, order)
        return [dict(lab=a["lab"], lead=a["lead"], defect=a["defect"], piv=a["piv"]) for a in a_]
    acc = cached("jets_" + tag, build_acc)
    print("lattice saturation done / loaded (%.1fs)" % (time.time() - t0))
    print("jet defects d_y (number of extra divisions by h), tally:", dict(sorted(Counter(a["defect"] for a in acc).items())))
    print("rank of the 64 jet leading covectors:", rank_of([a["lead"] for a in acc]))
    return acc

# ================================================================================================ B3
W_NAMES = [(0, (0, 2)), (1, (0, 2)), (0, (0, 3)), (1, (0, 3)), (1, (1, 2)), (0, (2, 3)), (1, (1, 3)), (1, (2, 3))]   # (c, (k,l)): P3^(c) p_kl
W_WEIGHT = [c + k + l for c, (k, l) in W_NAMES]                                                                  # 2,3,3,4,4,5,5,6


def mono_name(i):
    c, (k, l) = W_NAMES[i]
    return "P3%s p%d%d" % ("'" * c, k, l)


_PT = {}


def gaudin_point(tag="P"):
    if tag not in _PT:
        import twovac_fit as T
        import mpmath as mp
        T.set_precision(80)
        mp.mp.dps = 80
        thp = tuple(CH.thetas)
        pt = T.get_point(tag, th=thp)
        _PT[tag] = (T, pt, [T.StateJets(st, pt["zeta"], thp, n=6) for st in pt["states"]])
    return _PT[tag]


def mpf(x):
    import mpmath as mp
    if isinstance(x, Fr):
        return mp.mpf(x.numerator) / x.denominator
    return mp.mpf(x)


def psi_nb(pt):
    """Gaudin eigenkets (diagonal frame) in the notebook (xxx_chain) basis order, list of 64 lists of mp numbers."""
    out = []
    for st in pt["states"]:
        v = st["vec"]
        out.append([v[NB2G[a] * 8 + NB2G[k]] for a in range(8) for k in range(8)])
    return out


def mono_funcs(sj, s):
    """per state: list of the 8 one-site jet monomials P3^(c) p_kl at theta_s"""
    out = []
    for J in sj:
        P = J.P[s]
        pl = {}
        for k in range(0, 4):
            for l in range(k + 1, 4):
                pl[(k, l)] = P[1][k] * P[2][l] - P[2][k] * P[1][l]
        out.append([P[2][c] * pl[kl] for (c, kl) in W_NAMES])
    return out


def solve_cols(A, b):
    import mpmath as mp
    return mp.lu_solve(A, b)


def stage_B3(acc=None, tag="P"):
    import mpmath as mp
    hdr("B3: the Gaudin wavefunction rule in the jet basis (diagonal frame, P1 states, 80 digits)")
    if acc is None:
        acc = stage_B2b(tag)
    T, pt, sj = gaudin_point(tag)
    mp.mp.dps = 80
    psis = psi_nb(pt)
    zeta = pt["zeta"]
    z = [Ser.exp(Fr(x)) for x in zeta]
    xg = [mpf(c) for c in [(zeta[1] - zeta[0]), (zeta[2] - zeta[0]), (zeta[1] - zeta[0]), 2 * zeta[2] - zeta[0] - zeta[1], 3 * (zeta[1] - zeta[0]),
                           3 * (zeta[2] - zeta[0]), 2 * (zeta[2] - zeta[1]), 3 * (zeta[2] - zeta[1])]]
    xG2 = [a * b for a in xg for b in xg]
    nst = len(psis)
    den = [sum(a * b for a, b in zip(xG2, ps)) for ps in psis]
    m1, m2 = mono_funcs(sj, 0), mono_funcs(sj, 1)
    M1M2 = [m1[B][0] * m2[B][0] for B in range(nst)]
    # check x0^G of D1 (gl3 order) agrees with the nb-order x^G
    xD1 = T.x0_onesite(T.chain(tuple(CH.thetas)), zeta)
    xnb_as_gl3 = [None] * 8
    for a in range(8):
        xnb_as_gl3[NB2G[a]] = xg[a]
    r0 = [xD1[i] / xnb_as_gl3[i] for i in range(8)]
    print("x^G (notebook order, mapped by NB2G) equals D1's gl3-order x0_onesite up to one constant:", max(abs(r / r0[0] - 1) for r in r0) < mp.mpf(10) ** -60)
    # product matrix
    A = mp.matrix(nst, 64)
    for B in range(nst):
        for i in range(8):
            for j in range(8):
                A[B, i * 8 + j] = m1[B][i] * m2[B][j] / M1M2[B]
    sv = mp.svd_r(A, compute_uv=False)
    print("the 64 products m_i(theta_1) m_j(theta_2) (8 one-site jets per site) as functions on the 64 states: smallest/largest singular value %s (rank 64: %s)" % (
        mp.nstr(sv[63] / sv[0], 3), sv[63] / sv[0] > mp.mpf(10) ** -40))
    # solve for W
    Ainv_cache = mp.lu_solve  # noqa
    rows_out = []
    Wm = []
    t0 = time.time()
    # LU once
    for ent in acc:
        y = [mpf(c) for c in ent["lead"]]
        r = mp.matrix(nst, 1)
        for B in range(nst):
            num = sum(a * b for a, b in zip(y, psis[B]))
            r[B] = num / den[B]
        w = mp.lu_solve(A, r)
        Wm.append([w[i] for i in range(64)])
    print("W solved (%.1fs)" % (time.time() - t0))
    return acc, Wm, A, (m1, m2, psis, den, M1M2)


def mp_to_fr(x, maxden=10 ** 28, tol=None):
    """rational reconstruction of an mpf (80 digits): Fraction with denominator <= maxden if the residual is < tol, else None."""
    import mpmath as mp
    tol = tol or mp.mpf(10) ** -55
    if x == 0:
        return Fr(0)
    man, ex = mp.mpf(x).man, mp.mpf(x).exp
    f = Fr(int(man)) * (Fr(2) ** int(ex))
    if x < 0:
        f = -f
    r = f.limit_denominator(maxden)
    ref = max(abs(x), mp.mpf(1))
    if abs(x - mpf(r)) < tol * ref:
        return r
    return None


def compute_W(tag="P"):
    """exact rational 64x64 matrix W (rows: jet covectors y in the saturation basis, columns: products m_i(theta_1) m_j(theta_2), i,j = 0..7):
       <y|psi_B>/<x0^G|psi_B> = sum_ij W_{y,ij} m_i(th1) m_j(th2)/(M1 M2)  for all 64 Gaudin states B."""
    def build():
        acc, Wm, A, aux = stage_B3(tag=tag)
        W = []
        bad = 0
        for row in Wm:
            r = []
            for x in row:
                f = mp_to_fr(x)
                if f is None:
                    bad += 1
                    f = Fr(0)
                r.append(f)
            W.append(r)
        return dict(W=W, bad=bad, defects=[a["defect"] for a in acc], labs=[a["lab"] for a in acc])
    return cached("W_" + tag, build)


def frac_inverse(M):
    """exact inverse of a square Fraction matrix (list of lists)."""
    n = len(M)
    A = [list(r) + [Fr(int(i == j)) for j in range(n)] for i, r in enumerate(M)]
    for c in range(n):
        p = next(i for i in range(c, n) if A[i][c] != 0)
        A[c], A[p] = A[p], A[c]
        inv = 1 / A[c][c]
        A[c] = [x * inv for x in A[c]]
        for i in range(n):
            if i != c and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    return [r[n:] for r in A]


def frac_matmul(A, B):
    Bt = list(zip(*B))
    return [[sum((a * b for a, b in zip(r, c) if a != 0 and b != 0), Fr(0)) for c in Bt] for r in A]


def stage_B3b(tag="P"):
    hdr("B3 (result): exact wavefunction matrix W of the jet basis  [zeta point %s]" % tag)
    d = compute_W(tag)
    W, defects, labs = d["W"], d["defects"], d["labs"]
    print("rationalised entries: %d failures; nonzero entries: %d of 4096" % (d["bad"], sum(1 for r in W for x in r if x != 0)))
    wt2 = [W_WEIGHT[i] + W_WEIGHT[j] - 4 for i in range(8) for j in range(8)]       # product jet weight 0..8
    viol = sum(1 for y, r in enumerate(W) for m, x in enumerate(r) if x != 0 and wt2[m] > defects[y])
    print("filtration: W[y, (i,j)] = 0 whenever w_i + w_j - 4 > d_y (jet defect of y):", viol == 0, "(%d violations)" % viol)
    diag_ok = []
    from collections import Counter
    for lev in range(9):
        ys = [y for y in range(64) if defects[y] == lev]
        ms = [m for m in range(64) if wt2[m] == lev]
        blk = [[W[y][m] for m in ms] for y in ys]
        rk = rank_of(blk) if blk and blk[0] else 0
        diag_ok.append((lev, len(ys), len(ms), rk))
    print("graded blocks (level d: #jet covectors, #products of weight 4+d, rank of the leading block):", diag_ok)
    print("block W^(d,d) invertible at every level (the leading part of <y|psi> is a pure weight-(4+d) combination):", all(a[1] == a[2] == a[3] for a in diag_ok))
    return d


def stage_B3c():
    """compare W across zeta points"""
    hdr("B3 (zeta dependence): W at P, Pp, X1")
    ds = {}
    for tag in ("P", "Pp", "X1"):
        if os.path.exists(cache_path("jets_" + tag)):
            ds[tag] = compute_W(tag)
    base = ds["P"]["W"]
    for tag, d in ds.items():
        print(tag, "max |W - W_P| =", max(abs(a - b) for r1, r2 in zip(d["W"], base) for a, b in zip(r1, r2)), " defects equal:", d["defects"] == ds["P"]["defects"])
    return ds


def pinv_exact(z):
    """one-site P_phys^{-1} at the exact numbers z = (z1,z2,z3) (Fractions), notebook order, as a dense 8x8 Fraction matrix (column action on kets)."""
    z1, z2, z3 = [Fr(x) for x in z]
    E = CH.E1
    def expN(c, key):
        Ek = E[key]
        return SpOp.identity(8) + Ek.scale(c) + (Ek * Ek).scale(c * c * Fr(1, 2))
    s = [z1 ** 2 / ((z1 - z2) * (z1 - z3)), -(z2 / (z1 * (z2 - z3))), 1 / (z1 * z2)]
    r = [-(z3 / z2), -(z3 / z1), (z2 ** 2 * (z3 - z1)) / (z1 ** 2 * (z2 - z3))]
    pp = expN(-(z1 + z2) / (z1 * z2), (2, 3)) * expN(1 / (z1 * z2), (1, 3)) * expN(-(1 / z1), (1, 2))
    wts = [[int(CH.E1[(k, k)].r.get(p, {}).get(p, 0)) for k in (1, 2, 3)] for p in range(8)]
    dinv = SpOp(8, {p: {p: s[0] ** (-w[0]) * s[1] ** (-w[1]) * s[2] ** (-w[2])} for p, w in enumerate(wts)})
    pm = expN(-r[2], (2, 1)) * expN(-r[1], (3, 1)) * expN(-r[0], (3, 2))
    return (pp * dinv * pm).dense_list()


def kron_mat(A, B):
    n, m = len(A), len(B)
    return [[A[i // m][j // m] * B[i % m][j % m] for j in range(n * m)] for i in range(n * m)]


def jet_matrix(tag):
    acc = cached("jets_" + tag, None)
    return [a["lead"] for a in acc]


def adapted_covectors(tag):
    """exact covectors <mu| (rows indexed by mu = 8 i + j, one-site jet monomials m_i m_j) with <mu|psi>/<x0^G|psi> = m_i m_j/(M1 M2):
       <mu| = sum_y (W^{-1})_{mu y} <y|  (diagonal frame, notebook order)"""
    W = compute_W(tag)["W"]
    Y = jet_matrix(tag)
    Wi = frac_inverse(W)
    return frac_matmul(Wi, Y), Wi


def stage_B3d():
    hdr("B3 (frames): the monomial covectors <i,j| in the Gaudin companion frame at several zeta (is the rule zeta independent there?)")
    res = {}
    for tag in ("P", "Pp", "X1"):
        if not os.path.exists(cache_path("W_" + tag)):
            continue
        Ya, Wi = adapted_covectors(tag)
        zeta = ZETAS[tag]
        PG = pinv_exact(zeta)
        PGi = frac_inverse(PG)
        Mi = kron_mat(PGi, PGi)                    # (P (x) P)^-1
        Ycomp = frac_matmul(Ya, Mi)                # <mu|^c = <mu|^d (P_G^{-1} (x) P_G^{-1})^{-1}
        res[tag] = Ycomp
        # vacuum <x0|^c should be e1 (x) e1 up to scalar
        x0c = Ycomp[0]
        print(tag, "<0,0|^c (the x0 covector in the Gaudin companion frame) support:", [i for i, v in enumerate(x0c) if v != 0], "value", x0c[0])
    base = res["P"]
    for tag, Yc in res.items():
        if tag == "P":
            continue
        ratios = []
        okp = 0
        for mu in range(64):
            a, b = base[mu], Yc[mu]
            ia = next(i for i, v in enumerate(a) if v != 0)
            sc = b[ia] / a[ia]
            if all(b[i] == sc * a[i] for i in range(64)):
                okp += 1
            ratios.append(sc)
        print(tag, "monomial covectors <mu|^c proportional to those at P (exactly):", okp, "/ 64")
        print("   scalars (first 8):", [str(r) for r in ratios[:8]])
    return res


def comp_jet_matrix(tag="P"):
    """the Gaudin companion-frame monomial covectors <i,j|^c as an exact 64x64 rational matrix (rows mu = 8 i + j), and its inverse"""
    def build():
        Ya, Wi = adapted_covectors(tag)
        zeta = ZETAS[tag]
        PGi = frac_inverse(pinv_exact(zeta))
        Yc = frac_matmul(Ya, kron_mat(PGi, PGi))
        return dict(Y=Yc, Yinv=frac_inverse(Yc), W=compute_W(tag)["W"], Wi=Wi)
    return cached("Ycomp_" + tag, build)


def omega_comp_sym(a):
    """Omega_c^inf[a] = e1^T rho(e^{a E12} e^{a E23}) (notebook order), sympy vector, a a sympy symbol or number"""
    import sympy
    def mat(key):
        return sympy.Matrix(8, 8, lambda i, j: sympy.Rational(CH.E1[key].r.get(i, {}).get(j, Fr(0)).numerator, CH.E1[key].r.get(i, {}).get(j, Fr(0)).denominator))
    E12, E23 = mat((1, 2)), mat((2, 3))
    I8 = sympy.eye(8)
    e1 = sympy.Matrix([[1] + [0] * 7])
    return list(e1 * (I8 + a * E12 + a ** 2 * E12 ** 2 / 2) * (I8 + a * E23 + a ** 2 * E23 ** 2 / 2))


N5_IDX = [0, 1, 2, 3, 5]          # W_NAMES indices of N = (P3 p02, P3' p02, P3 p03, P3' p03, P3 p23)


def stage_B4(tag="P"):
    import sympy
    hdr("B4: the Gaudin K tensor: Phi_c^inf = Omega_c^inf[a] (x) Omega_c^inf[a'] expanded in the monomial (jet) basis")
    d = comp_jet_matrix(tag)
    Y, Yinv = d["Y"], d["Yinv"]
    nzY = sum(1 for r in Y for x in r if x != 0)
    nzI = sum(1 for r in Yinv for x in r if x != 0)
    print("companion-frame monomial covectors <i,j|^c (zeta independent, exact): %d nonzero entries of 4096; inverse: %d nonzero" % (nzY, nzI))
    a, ap, thi = sympy.symbols("a ap thi")
    Om, Omp = omega_comp_sym(a), omega_comp_sym(ap)
    print("Omega_c^inf[a] components (notebook order):", [sympy.factor(x) for x in Om])
    Yi = sympy.Matrix(64, 64, lambda i, j: sympy.Rational(Yinv[i][j].numerator, Yinv[i][j].denominator))
    Phi = sympy.Matrix([[Om[i] * Omp[k] for i in range(8) for k in range(8)]])
    Theta = (Phi * Yi).applyfunc(sympy.expand)
    nz = [(m // 8, m % 8) for m in range(64) if Theta[m] != 0]
    print("nonzero K-tensor entries Theta_(i,j) (jet monomials m_i(th1) m_j(th2)): %d of 64" % len(nz))
    print("   support (i,j) with i,j in W_NAMES order:", nz)
    print("   in terms of N5 = (P3p02, P3'p02, P3p03, P3'p03, P3p23) ~ W indices", N5_IDX, ": all support inside N5 x N5:", all(i in N5_IDX and j in N5_IDX for i, j in nz))
    # polynomial structure
    thi_val = Fr(1) / (CH.thetas[0] - CH.thetas[1])
    print("1/theta12 =", thi_val)
    bad = compare_with_G(Theta, a, ap, thi_val)
    print("Theta_(N5_i, N5_j)(a, a') == G_ij(1/a, 1/a', 1/theta12) a^3 a'^3 / 4  (= G/(4 w1^3 w1'^3)), exact polynomial identity in (a, a'), entries failing:", bad, "(of 25)")
    return d, Theta, nz


_G2SYM = None


def d2_G_sym():
    """D2's kernel G_ij(w1, w1p, thi) (symbolic, from RESULTS-twovac-limit.md, read only): dict (i,j) -> sympy expr, plus the symbols"""
    global _G2SYM
    import re
    import sympy
    from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication
    w, wp, t = sympy.symbols("w1 w1p thi")
    if _G2SYM is None:
        txt = open(os.path.join(HERE, "RESULTS-twovac-limit.md"), encoding="utf8").read()
        blk = txt[txt.index("G[1,1] = -4 w1^2"):]
        blk = blk[:blk.index("```")]
        G = {}
        for line in blk.splitlines():
            m = re.match(r"G\[(\d),(\d)\] = (.*)", line.strip())
            if m:
                G[(int(m.group(1)) - 1, int(m.group(2)) - 1)] = parse_expr(m.group(3).replace("^", "**"), local_dict=dict(w1=w, w1p=wp, thi=t),
                                                                          transformations=standard_transformations + (implicit_multiplication,))
        _G2SYM = G
    return _G2SYM, (w, wp, t)


def compare_with_G(Theta, a, ap, thi_val, label=""):
    import sympy
    G, (w, wp, t) = d2_G_sym()
    bad = []
    for i in range(5):
        for j in range(5):
            mu = N5_IDX[i] * 8 + N5_IDX[j]
            gij = G[(i, j)].subs({w: 1 / a, wp: 1 / ap, t: sympy.Rational(thi_val.numerator, thi_val.denominator)})
            diff = sympy.simplify(Theta[mu] - gij * a ** 3 * ap ** 3 / 4)
            if diff != 0:
                bad.append((i, j))
    return bad


def comp_kets(tag):
    """Gaudin companion-frame kets Psi_c = (P_G^-1 (x) P_G^-1) psi_d (notebook order), P_G^-1 = P_phys^{-1} at z = zeta; list over the 64 states"""
    import mpmath as mp
    T, pt, sj = gaudin_point(tag)
    mp.mp.dps = 80
    PG = pinv_exact(ZETAS[tag])
    PGm = mp.matrix(8, 8)
    for i in range(8):
        for j in range(8):
            PGm[i, j] = mpf(PG[i][j])
    out = []
    for ps in psi_nb(pt):
        M = mp.matrix(8, 8)
        for i in range(8):
            for k in range(8):
                M[i, k] = ps[8 * i + k]
        C = PGm * M * PGm.T
        out.append([C[i, k] for i in range(8) for k in range(8)])
    return out


def stage_B4b():
    import mpmath as mp
    hdr("B3/B4 verification: the zeta-independent companion-frame monomial covectors vs the Gaudin states at P, Pp, X1 (80 digits)")
    d = comp_jet_matrix("P")
    Y = d["Y"]
    mp.mp.dps = 80
    Ymp = [[mpf(x) for x in r] for r in Y]
    for tag in ("P", "Pp", "X1"):
        T, pt, sj = gaudin_point(tag)
        kets = comp_kets(tag)
        m1, m2 = mono_funcs(sj, 0), mono_funcs(sj, 1)
        worst = mp.mpf(0)
        worstx0 = mp.mpf(0)
        for B, ket in enumerate(kets):
            x0 = ket[0]                       # <x0^c| = e1 (x) e1
            M12 = m1[B][0] * m2[B][0]
            for mu in range(64):
                lhs = sum((y * k for y, k in zip(Ymp[mu], ket) if y != 0), mp.mpf(0)) / x0
                rhs = m1[B][mu // 8] * m2[B][mu % 8] / M12
                worst = max(worst, abs(lhs - rhs) / max(abs(rhs), mp.mpf(10) ** -30))
        print("zeta point %s (%s): max over 64 states x 64 jets of | <mu|^c Psi_c/<x0|Psi_c - m_i(th1) m_j(th2)/(M1 M2) | / | . | = %s" % (tag, ZETAS[tag], mp.nstr(worst, 3)))
    return d


def stage_B4c():
    import sympy
    import mpmath as mp
    hdr("B4: homogeneity, structure of the K tensor entries Theta_ij(a, a', 1/theta12)")
    d, Theta, nz = stage_B4()
    a, ap = sympy.symbols("a ap")
    thi = Fr(1) / (CH.thetas[0] - CH.thetas[1])
    W5 = [W_WEIGHT[i] for i in N5_IDX]
    print("N5 jet weights:", W5)
    # degree structure: each Theta_ij is a polynomial in (a, a') with thi numeric; total degree and lowest degree
    rows = []
    for i in range(5):
        for j in range(5):
            mu = N5_IDX[i] * 8 + N5_IDX[j]
            P = sympy.Poly(Theta[mu], a, ap)
            degs = [m[0] + m[1] for m in P.monoms()]
            rows.append((i, j, min(degs), max(degs), len(P.terms())))
    print("entries (i,j) of the N5 x N5 block: lowest/highest total degree in (a,a') and number of monomials (thi = 1/theta12 numeric):")
    for r in rows:
        print("   (%d,%d)  deg_in_(a,a'): %d..%d   #terms %d   [w_i+w_j-4 = %d]" % (r[0], r[1], r[2], r[3], r[4], W5[r[0]] + W5[r[1]] - 4))
    return Theta


OM_EQ = (Fr(5, 2), Fr(7, 3), Fr(3, 5))


def omega_D_nb(zeta, om):
    """one-site diagonal-frame vacuum covector Omega_D[zeta, om] = e1^T rho(exp(r1 E12) exp(r2 E13) exp(r3 E23)) (notebook order, exact)"""
    r = X.gauss_r(tuple(Fr(z) for z in zeta), tuple(Fr(w) for w in om))
    return X.covec_Omega(CH, r)


def theta_of_cov(tag, cov2):
    """jet-monomial coefficients vartheta (64, exact) of a two-site diagonal-frame covector (notebook order, Fractions):
       cov2 (P_G (x) P_G)^-1 = sum_mu vartheta_mu <mu|^c   (so that <cov2|psi>/<x0^G|psi> = sum_mu vartheta_mu m_mu/(M1 M2))"""
    d = comp_jet_matrix("P")
    Yinv = d["Yinv"]
    PGi = frac_inverse(pinv_exact(ZETAS[tag]))
    c = [Fr(0)] * 64
    # c = cov2 (PG^-1 (x) PG^-1)^{-1}... (P_G^{-1}) is pinv_exact; the companion covector is cov2 . (pinv_exact(zeta) (x) pinv_exact(zeta))^{-1}
    Mi = kron_mat(PGi, PGi)
    cc = [sum((cov2[i] * Mi[i][j] for i in range(64) if cov2[i] != 0), Fr(0)) for j in range(64)]
    return [sum((cc[i] * Yinv[i][mu] for i in range(64) if cc[i] != 0), Fr(0)) for mu in range(64)]


def stage_B4d(tag="P"):
    import mpmath as mp
    hdr("B4: equal-twist family Phi = Omega_D[zeta,w] (x) Omega_D[zeta,w] (finite w2): K tensor and the Q-number against P1's det N")
    T, pt, sj = gaudin_point(tag)
    mp.mp.dps = 80
    zeta = ZETAS[tag]
    om0 = OM_EQ
    o1 = omega_D_nb(zeta, om0)
    cov2 = [a * b for a in o1 for b in o1]
    th = theta_of_cov(tag, cov2)
    nz = [mu for mu in range(64) if th[mu] != 0]
    print("equal twist w =", om0, ": nonzero jet coefficients vartheta_(i,j):", len(nz), "of 64")
    lev = {}
    for mu in nz:
        lev.setdefault(W_WEIGHT[mu // 8] + W_WEIGHT[mu % 8] - 4, []).append(mu)
    print("   by jet level d = w_i + w_j - 4:", {k: len(v) for k, v in sorted(lev.items())})
    print("   one-site monomials involved on site 1:", sorted({mu // 8 for mu in nz}), "site 2:", sorted({mu % 8 for mu in nz}))
    # numerical check against the states
    psis = psi_nb(pt)
    xg = omega_x0_nb(zeta)
    m1, m2 = mono_funcs(sj, 0), mono_funcs(sj, 1)
    qn, qj = [], []
    for B, ps in enumerate(psis):
        num = sum(mpf(c) * v for c, v in zip(cov2, ps) if c != 0)
        den = sum(mpf(c) * v for c, v in zip(xg, ps))
        q = num / den
        M12 = m1[B][0] * m2[B][0]
        qj_ = sum(mpf(th[mu]) * m1[B][mu // 8] * m2[B][mu % 8] for mu in nz) / M12
        qn.append(q)
        qj.append(qj_)
    print("Q-number from the pairing vs from the jet expansion sum_mu vartheta_mu m_mu/(M1 M2): max rel. deviation over 64 states: %s" % mp.nstr(max(abs(a / b - 1) for a, b in zip(qn, qj)), 3))
    D, _ = T.detN_for_states(pt, om0, tuple(CH.thetas))
    M12s = [m1[B][0] * m2[B][0] for B in range(64)]
    ratio = [qn[B] * M12s[B] / D[B] for B in range(64)]
    r0 = ratio[0]
    dev = max(abs(r / r0 - 1) for r in ratio)
    print("Qnum_eq(B) M1 M2 / det N_B(w): state-independent? max |ratio/ratio_0 - 1| = %s ; ratio_0 = %s" % (mp.nstr(dev, 3), mp.nstr(r0, 15)))
    return th, ratio


def omega_x0_nb(zeta):
    z1, z2, z3 = [Fr(x) for x in zeta]
    x = [z2 - z1, z3 - z1, z2 - z1, 2 * z3 - z1 - z2, 3 * (z2 - z1), 3 * (z3 - z1), 2 * (z3 - z2), 3 * (z3 - z2)]
    return [a * b for a in x for b in x]


# ================================================================================================ B1b / B2c / B5: grading, Lie algebra, leading covectors, candidates
def htot_depth():
    """one-site depth 5 - (2 n1 + n2) of the 8 GT states (notebook order) and the two-site depth of the 64 components"""
    H = [int(2 * CH.E1[(1, 1)].r.get(p, {}).get(p, 0) + CH.E1[(2, 2)].r.get(p, {}).get(p, 0)) for p in range(8)]
    dep = [5 - h for h in H]
    return dep, [dep[c // 8] + dep[c % 8] for c in range(64)]


def mat_dense(M, idx):
    return [[M.r.get(i, {}).get(j, Fr(0)) for j in idx] for i in idx]


def mbr(A, B):
    n = len(A)
    AB = frac_matmul(A, B)
    BA = frac_matmul(B, A)
    return [[AB[i][j] - BA[i][j] for j in range(n)] for i in range(n)]


def flat(A):
    return [x for r in A for x in r]


def lie_closure(gens, maxdim=400):
    """Lie algebra generated by square Fraction matrices: returns the list of independent matrices (basis), by left-normed brackets with the generators"""
    n = len(gens[0])
    basis, vecs, piv = [], [], {}
    def add(A):
        v = flat(A)
        v = list(v)
        for k, (pv, pi) in enumerate(vecs):
            if v[pi] != 0:
                f = v[pi] / pv[pi]
                v = [a - f * b for a, b in zip(v, pv)]
        for i, x in enumerate(v):
            if x != 0:
                vecs.append((v, i))
                basis.append(A)
                return True
        return False
    front = []
    for g in gens:
        if add(g):
            front.append(g)
    while front and len(basis) < maxdim:
        new = []
        for x in front:
            for g in gens:
                c = mbr(x, g)
                if add(c):
                    new.append(c)
        front = new
    return basis


def derived(basis):
    out = []
    vecs = []
    def add(A):
        v = list(flat(A))
        for pv, pi in vecs:
            if v[pi] != 0:
                f = v[pi] / pv[pi]
                v = [a - f * b for a, b in zip(v, pv)]
        for i, x in enumerate(v):
            if x != 0:
                vecs.append((v, i))
                out.append(A)
                return
    for i in range(len(basis)):
        for j in range(i + 1, len(basis)):
            add(mbr(basis[i], basis[j]))
    return out


def stage_B1b():
    hdr("B1: the H_tot = 2E11 + E22 grading of the b_k and simultaneous triangularisability (Lie algebra generated by the b_k(u))")
    dep, cdep = htot_depth()
    levels = sorted(set(cdep))
    idxs = {D: [c for c in range(64) if cdep[c] == D] for D in levels}
    print("one-site H depths (5 - (2n1+n2)) in notebook order:", dep, "; two-site level sizes:", {D: len(v) for D, v in idxs.items()})
    ok = True
    for u in (U0, U1):
        Bk = Bcoeffs(u)
        for k, M in enumerate(Bk):
            for i, row in M.r.items():
                for j in row:
                    if cdep[i] != cdep[j]:
                        ok = False
    print("every B_k(u) is block diagonal in the two-site H_tot-grading (B commutes with H1 + H2 = 2E11 + E22 summed over sites):", ok)
    # coefficient matrices of u^m in B_k(u) (B_k is a polynomial of degree <= 6 in u): span_u {b_k(u)} = span_m {B_k^{(m)}}
    us = [Fr(1, 7) + Fr(m, 3) for m in range(7)]
    Bus = [Bcoeffs(u) for u in us]
    # Lagrange -> monomial coefficients
    def lag_coeffs(i):
        poly = [Fr(1)]
        den = Fr(1)
        for l in range(7):
            if l != i:
                poly = pmul_(poly, [-us[l], Fr(1)])
                den *= us[i] - us[l]
        return [c / den for c in poly]
    LC = [lag_coeffs(i) for i in range(7)]
    coef = {}
    for k in range(1, 7):
        for m in range(7):
            tot = None
            for i in range(7):
                t = Bus[i][k].scale(LC[i][m])
                tot = t if tot is None else tot + t
            coef[(k, m)] = tot
    res = []
    for D in levels:
        idx = idxs[D]
        n = len(idx)
        gens2 = []
        for key, M in coef.items():
            g = mat_dense(M, idx)
            tr = sum(g[i][i] for i in range(n)) / n
            g2 = [[g[i][j] - (tr if i == j else 0) for j in range(n)] for i in range(n)]
            if any(x != 0 for r in g2 for x in r):
                gens2.append(g2)
        if not gens2:
            print("   level D=%d (block size %d): all B_k^(m) are scalar on this block" % (D, n))
            res.append((D, n, 0, 0, 0, 0, 0))
            continue
        L = lie_closure(gens2)
        d1 = derived(L)
        d2 = derived(d1) if d1 else []
        d3 = derived(d2) if d2 else []
        d4 = derived(d3) if d3 else []
        res.append((D, n, len(L), len(d1), len(d2), len(d3), len(d4)))
        print("   level D=%d (block size %d): dim of the Lie algebra generated by all B_k^(m) (k=1..6, u^m coefficients, traceless parts) = %d; derived series dims: %d, %d, %d, %d" % (D, n, len(L), len(d1), len(d2), len(d3), len(d4)), flush=True)
    return res


def stage_B2c():
    hdr("B2: grading, leading covectors as joint eigenvectors, filtration ranks")
    covs = get_covs()
    dep, cdep = htot_depth()
    # each f_x supported in a single graded level D(x)
    ok = True
    for lab, (nums, q) in covs.items():
        sup = {cdep[c] for c, P in enumerate(nums) if any(x != 0 for x in P)}
        if sup != {depth2(lab)}:
            ok = False
    print("each SoV covector f_x(h) is supported on a single H_tot level, equal to the label depth D(x) = sum_s (6 - (1+a+b+c)):", ok)
    # leading covectors: left eigenvectors of b1 and of the b1-compressed b2
    Bk = Bcoeffs(U0)
    b = [M.scale(1 / nu(U0) ** 3) for M in Bk]
    d1 = [b[1].r.get(i, {}).get(i, Fr(0)) for i in range(N64)]
    blocks = eigenclass_blocks(d1)
    cls_of = {}
    for key, idx in blocks.items():
        for i in idx:
            cls_of[i] = key
    # compressed b2: entries (i,j) with equal b1 class
    comp = SpOp(N64, {i: {j: v for j, v in row.items() if cls_of[i] == cls_of[j]} for i, row in b[2].r.items()})
    K = 4
    ser = {lab: cov_series(nums, q, K) for lab, (nums, q) in covs.items()}
    ok1 = ok2 = 0
    for lab in LABELS2:
        v0 = ser[lab][0]
        r1 = b[1].vec_times(v0)
        l1 = sum(Fr(1 + sum(l)) / (U0 - TH[s]) for s, l in enumerate(lab))
        ok1 += all(r1[i] == l1 * v0[i] for i in range(N64))
        r2 = comp.vec_times(v0)
        i0 = next(i for i in range(N64) if v0[i] != 0)
        lam = r2[i0] / v0[i0]
        ok2 += all(r2[i] == lam * v0[i] for i in range(N64))
    print("leading covectors <x|_0 (companion frame): left eigenvectors of b1(u0): %d/64 ; of the b1-block-compressed b2 (sum_S P_S b2 P_S): %d/64" % (ok1, ok2))
    # diagonal frame filtration ranks of the Taylor coefficients of R_x = f^d_x / h^2
    set_kser(KSER)
    P = sp_dense(pinv_series(ZETAS["P"]))
    rows = diag_rows(covs, P)
    R = {lab: row_shift(rows[lab], 2) for lab in LABELS2}
    out = []
    for kmax in range(0, 10):
        mat = [[x.c[k] for x in R[lab]] for lab in LABELS2 for k in range(kmax + 1)]
        out.append(rank_of(mat))
    print("diagonal frame, zeta = P: rank span{ R_x,k : k <= K } for K = 0..9 :", out, "(R_x = f^d_x/h^2; lattice-saturated cumulative jet counts are 1,5,13,25,39,51,59,63,64)")
    return out


def stage_B5():
    hdr("B5: operators b_k(u) in the Gaudin companion-frame jet basis  (b_k acting on covectors from the right; matrix B_k = Y b_k Y^-1)")
    d = comp_jet_matrix("P")
    Y, Yinv = d["Y"], d["Yinv"]
    dep, cdep = htot_depth()
    lev = [W_WEIGHT[mu // 8] + W_WEIGHT[mu % 8] - 4 for mu in range(64)]
    cls = [(dep_i, dep_j) for dep_i in [0] for dep_j in [0]]  # placeholder
    onesite_depth = [w - 2 for w in W_WEIGHT]       # jet depth of the 8 one-site monomials (0,1,1,2,2,3,3,4)
    cls = [(onesite_depth[mu // 8], onesite_depth[mu % 8]) for mu in range(64)]
    for u, tag in ((U0, "u0"), (U1, "u1")):
        Bk = Bcoeffs(u)
        bk = [M.scale(1 / nu(u) ** 3) for M in Bk]
        print("u =", u)
        for k in (1, 2, 3):
            # jet-basis matrix: B = Y b Y^-1
            bd = bk[k].dense_list()
            Bm = frac_matmul(frac_matmul(Y, bd), Yinv)
            tri_lev = all(Bm[m][n] == 0 for m in range(64) for n in range(64) if lev[n] > lev[m])
            nz = sum(1 for r in Bm for x in r if x != 0)
            diag = [Bm[m][m] for m in range(64)]
            extra = ""
            if k == 1:
                lam = lambda c: sum(Fr(6 - c[s_]) / (u - TH[s_]) for s_ in range(2))
                # eigenvalue per site: w_s = 1 + a + b + c = 6 - depth
                extra = "; diagonal = sum_s (6 - depth_s)/(u-theta_s) of the monomial class: %s" % all(diag[m] == lam(cls[m]) for m in range(64))
                # class triangularity (partial order): entry (m,n) nonzero only if cls[n] <= cls[m] componentwise
                tri_cls = all(Bm[m][n] == 0 or (cls[n][0] <= cls[m][0] and cls[n][1] <= cls[m][1]) for m in range(64) for n in range(64))
                extra += "; triangular w.r.t. the (depth_1, depth_2) class partial order: %s" % tri_cls
            print("   b_%d: nonzero %d of 4096; block lower-triangular w.r.t. the jet level (b_k F_d subset F_d): %s%s" % (k, nz, tri_lev, extra))
    return d


def stage_B4e(tag="P"):
    import sympy
    hdr("B4: frame tie scalar kappa^G (origin of c) and the K tensor in the XXX-derived (saturation) jet basis")
    zeta = ZETAS[tag]
    PG = pinv_exact(zeta)
    z1, z2, z3 = zeta
    for w in (Fr(5, 2), Fr(-2, 5), Fr(9, 7)):
        a = 1 / w
        omc = X.covec_Omega(CH, (a, Fr(0), a))                     # Omega_c^inf[a]: e1^T rho(e^{a E12} e^{a E23})
        omd = omega_D_nb(zeta, (w, Fr(10) ** 12, Fr(3, 5)))        # placeholder, replaced below by the exact infinity formulas
        r1 = (w - z1) / (w - z2)
        r2 = -(w - z1) * (z2 - z3) / ((z1 - z2) * (w - z3))
        r3 = (w - z2) * (z1 - z3) / ((z1 - z2) * (w - z3))
        omd = X.covec_Omega(CH, (r1, r2, r3))                      # Omega_D^inf[zeta, w]
        lhs = [sum((omc[i] * PG[i][j] for i in range(8)), Fr(0)) for j in range(8)]      # Omega_c (P_G^-1): the diagonal-frame image of the companion covector
        ratios = {lhs[j] / omd[j] for j in range(8) if omd[j] != 0}
        kap = w ** 3 / ((z2 - w) * (z3 - w) ** 2 * (z1 - z2))
        print("   w = %s: Omega_c^inf[1/w] P_G^-1 = rho * Omega_D^inf[zeta,w] with the single scalar rho = %s ; D2/W2 kappa^G = w^3/((z2-w)(z3-w)^2(z1-z2)) = %s ; rho * kappa^G == 1: %s" % (
            w, sorted(ratios), kap, ratios == {1 / kap}))
    d0 = comp_jet_matrix(tag)
    rk = []
    for mu in range(64):
        M8 = [[d0["Y"][mu][8 * i + k] for k in range(8)] for i in range(8)]
        rk.append(rank_of(M8))
    from collections import Counter as _C
    print("   companion-frame monomial covectors <i,j|^c viewed as 8x8 matrices (site 1 x site 2): rank tally", dict(sorted(_C(rk).items())),
          "(rank 1 = product covector; rank >= 2 = the covector couples the two sites)")
    print("   rank-1 (decoupled) jets:", [(mu // 8, mu % 8) for mu in range(64) if rk[mu] == 1][:40], "...")
    # K tensor in the saturation basis: phi^sat = Theta W^{-1}
    d = comp_jet_matrix(tag)
    W, Wi = d["W"], d["Wi"]
    a, ap = Fr(2, 3), Fr(-5, 7)
    o1 = X.covec_Omega(CH, (a, Fr(0), a))
    o2 = X.covec_Omega(CH, (ap, Fr(0), ap))
    Phi = [x * y for x in o1 for y in o2]
    Theta = [sum((Phi[i] * d["Yinv"][i][mu] for i in range(64) if Phi[i] != 0), Fr(0)) for mu in range(64)]
    phi = [sum((Theta[mu] * Wi[mu][y] for mu in range(64) if Theta[mu] != 0), Fr(0)) for y in range(64)]
    defects = compute_W(tag)["defects"]
    print("   (a,a') = (2/3,-5/7): nonzero Theta_mu (adapted monomial basis): %d ; nonzero phi_y in the saturation (XXX-derived) jet basis: %d of 64" % (
        sum(1 for x in Theta if x != 0), sum(1 for x in phi if x != 0)))
    from collections import Counter
    print("   phi_y nonzero by jet defect level:", dict(sorted(Counter(defects[y] for y in range(64) if phi[y] != 0).items())))
    return phi


def stage_B3e():
    hdr("B3: table of the saturation jet basis y -> products m_i(theta_1) m_j(theta_2)   (zeta = P; all 64 rows)")
    d = compute_W("P")
    W, defects, labs = d["W"], d["defects"], d["labs"]
    print("one-site monomials m_i:", ", ".join("%d: %s (w=%d)" % (i, mono_name(i), W_WEIGHT[i]) for i in range(8)))
    for y in range(64):
        terms = [(m // 8, m % 8, W[y][m]) for m in range(64) if W[y][m] != 0]
        top = [t for t in terms if W_WEIGHT[t[0]] + W_WEIGHT[t[1]] - 4 == defects[y]]
        print("y=%2d %s d=%d  top-level part: %s   (%d terms in all)" % (y, labs[y], defects[y], " + ".join("(%s)m%d m%d" % (t[2], t[0], t[1]) for t in top), len(terms)))


STAGES = {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stages", default="B0,B1,B2,B2b,B2c,B3,B3c,B3d,B3e,B4,B4b,B4c,B4d,B4e,B5,B1b")
    a = ap.parse_args()
    st = a.stages.split(",")
    if "B0" in st:
        stage_B0()
    if "B1" in st:
        stage_B1()
    if "B2" in st:
        stage_B2()
    if "B2b" in st:
        stage_B2b()
    if "B3" in st:
        stage_B3b("P")
    if "B3c" in st:
        stage_B3c()
    if "B3d" in st:
        stage_B3d()
    if "B4" in st:
        stage_B4()
    if "B4b" in st:
        stage_B4b()
    if "B4c" in st:
        stage_B4c()
    if "B4d" in st:
        stage_B4d()
    if "B1b" in st:
        stage_B1b()
    if "B5" in st:
        stage_B5()
    if "B4e" in st:
        stage_B4e()
    if "B3e" in st:
        stage_B3e()
    if "B2c" in st:
        stage_B2c()


if __name__ == "__main__":
    main()
