"""XXX gl(3) spin-chain port (new conventions) for the Gaudin-limit study (agent P2).

Conventions (Paul/Mathematica/new-conventions):
  Lax   L_ij(u) = u delta_ij + h E_ij ;  T(u) = L^(1)(u-theta_1) L^(2)(u-theta_2) ... ;  twist on the RIGHT: T.G
  quantum minors qmin_{I,J}(u) = sum_sigma sgn(sigma) prod_k T_{sigma(i_k), j_k}(u-(k-1)h)  (rows antisymmetrised)
  t_a(u) = sum_{|I|=a} qmin_{I,I}(u),  nu_n(u) = prod_s (u - theta_s + h Lam_n).
Generators: plain (non-Hermitian) GT generators from gl3_reps.py (E_13, E_31 by commutators, as in the notebook JJ).
Everything is exact: Fractions in a sparse operator class (SpOp), and truncated power series in h (Ser).
"""
import os
import sys
from fractions import Fraction as Fr
from itertools import permutations, combinations

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from gl3_reps import GL3Rep  # noqa: E402

# ----------------------------------------------------------------------------------------------
# truncated power series in h with Fraction coefficients (valuation >= 0), absolute order NSER
# ----------------------------------------------------------------------------------------------
NSER = 16


class Ser:
    __slots__ = ("c",)

    def __init__(self, c):
        c = list(c)[:NSER]
        c += [Fr(0)] * (NSER - len(c))
        self.c = c

    @staticmethod
    def const(x):
        return Ser([Fr(x)])

    @staticmethod
    def h():
        return Ser([Fr(0), Fr(1)])

    @staticmethod
    def exp(a):
        """exp(a*h) as a series (a Fraction)."""
        out, term = [], Fr(1)
        for k in range(NSER):
            out.append(term)
            term = term * Fr(a) / (k + 1)
        return Ser(out)

    def _co(self, o):
        return o if isinstance(o, Ser) else Ser([Fr(o)])

    def __add__(self, o):
        o = self._co(o)
        return Ser([a + b for a, b in zip(self.c, o.c)])

    __radd__ = __add__

    def __sub__(self, o):
        o = self._co(o)
        return Ser([a - b for a, b in zip(self.c, o.c)])

    def __rsub__(self, o):
        return self._co(o) - self

    def __neg__(self):
        return Ser([-a for a in self.c])

    def __mul__(self, o):
        if not isinstance(o, Ser):
            o = Fr(o)
            return Ser([a * o for a in self.c])
        out = [Fr(0)] * NSER
        for i, a in enumerate(self.c):
            if a == 0:
                continue
            for j in range(NSER - i):
                b = o.c[j]
                if b != 0:
                    out[i + j] += a * b
        return Ser(out)

    __rmul__ = __mul__

    def val(self):
        for i, a in enumerate(self.c):
            if a != 0:
                return i
        return None

    def divh(self, k=1):
        """exact division by h^k (the k lowest coefficients must vanish); top k orders become unknown (set 0)."""
        assert all(a == 0 for a in self.c[:k]), "divh: not divisible"
        return Ser(self.c[k:])

    def inv(self):
        a0 = self.c[0]
        assert a0 != 0
        out = [Fr(0)] * NSER
        out[0] = 1 / a0
        for n in range(1, NSER):
            s = sum(self.c[k] * out[n - k] for k in range(1, n + 1))
            out[n] = -s / a0
        return Ser(out)

    def __truediv__(self, o):
        if not isinstance(o, Ser):
            return self * (1 / Fr(o))
        return self * o.inv()

    def __rtruediv__(self, o):
        return self._co(o) * self.inv()

    def __pow__(self, n):
        r = Ser.const(1)
        for _ in range(n):
            r = r * self
        return r

    def __eq__(self, o):
        o = self._co(o)
        return all(a == b for a, b in zip(self.c, o.c))

    def __hash__(self):
        return hash(tuple(self.c))

    def __repr__(self):
        return "Ser(" + ", ".join(str(a) for a in self.c[:6]) + ", ...)"


def ratio_h(a, b):
    """a/b for series that vanish to the same (or higher) order at h = 0: divide out h^val(b) first."""
    va, vb = a.val(), b.val()
    if va is None:
        return Ser.const(0)
    assert vb is not None and va >= vb, "ratio_h: denominator more singular than numerator"
    return a.divh(vb) / b.divh(vb)


# ----------------------------------------------------------------------------------------------
# sparse operators with exact entries (Fraction or Ser)
# ----------------------------------------------------------------------------------------------
class SpOp:
    __slots__ = ("n", "r")

    def __init__(self, n, r=None):
        self.n = n
        self.r = r if r is not None else {}

    @staticmethod
    def identity(n):
        return SpOp(n, {i: {i: Fr(1)} for i in range(n)})

    @staticmethod
    def from_dense(M):
        n = len(M)
        r = {}
        for i in range(n):
            row = {j: Fr(M[i][j]) for j in range(n) if M[i][j] != 0}
            if row:
                r[i] = row
        return SpOp(n, r)

    def copy(self):
        return SpOp(self.n, {i: dict(row) for i, row in self.r.items()})

    def scale(self, s):
        if isinstance(s, (int, Fr)) and s == 0:
            return SpOp(self.n)
        return SpOp(self.n, {i: {j: v * s for j, v in row.items()} for i, row in self.r.items()})

    def __add__(self, o):
        if isinstance(o, (int, Fr, Ser)):
            o = SpOp.identity(self.n).scale(o)
        out = {i: dict(row) for i, row in self.r.items()}
        for i, row in o.r.items():
            d = out.setdefault(i, {})
            for j, v in row.items():
                d[j] = d[j] + v if j in d else v
        return SpOp(self.n, out)

    __radd__ = __add__

    def __neg__(self):
        return self.scale(-1)

    def __sub__(self, o):
        if isinstance(o, (int, Fr, Ser)):
            o = SpOp.identity(self.n).scale(o)
        return self + o.scale(-1)

    def __mul__(self, o):
        if not isinstance(o, SpOp):
            return self.scale(o)
        out = {}
        for i, row in self.r.items():
            acc = {}
            for k, a in row.items():
                brow = o.r.get(k)
                if not brow:
                    continue
                for j, b in brow.items():
                    t = a * b
                    acc[j] = acc[j] + t if j in acc else t
            if acc:
                out[i] = acc
        return SpOp(self.n, out)

    def __rmul__(self, s):
        return self.scale(s)

    def kron(self, o):
        n2 = o.n
        out = {}
        for i, row in self.r.items():
            for i2, row2 in o.r.items():
                d = {}
                for j, a in row.items():
                    for j2, b in row2.items():
                        d[j * n2 + j2] = a * b
                out[i * n2 + i2] = d
        return SpOp(self.n * n2, out)

    def dense_list(self):
        M = [[Fr(0)] * self.n for _ in range(self.n)]
        for i, row in self.r.items():
            for j, v in row.items():
                M[i][j] = v
        return M

    def vec_times(self, v):
        """row vector v (list) times operator -> list (covector action)."""
        out = [Fr(0)] * self.n
        for i, row in self.r.items():
            vi = v[i]
            if _is_zero(vi):
                continue
            for j, a in row.items():
                out[j] = out[j] + vi * a
        return out

    def times_vec(self, v):
        out = [Fr(0)] * self.n
        for i, row in self.r.items():
            s = Fr(0)
            for j, a in row.items():
                s = s + a * v[j]
            out[i] = s
        return out

    def is_zero(self):
        return all(all(_is_zero(v) for v in row.values()) for row in self.r.values())

    def nnz(self):
        return sum(len(r) for r in self.r.values())

    def flat(self):
        """dict {(i,j): v} of nonzero entries."""
        return {(i, j): v for i, row in self.r.items() for j, v in row.items() if not _is_zero(v)}


def op_commutator(a, b):
    return a * b - b * a


# ----------------------------------------------------------------------------------------------
# gl(3) one-site generators, exact
# ----------------------------------------------------------------------------------------------
def _frac(x):
    return Fr(int(x.p), int(x.q))


def gl3_generators(lam):
    rep = GL3Rep(lam)
    d = rep.dim
    # gl3_reps lists GT patterns in ASCENDING lexicographic order (index 0 = lowest weight); the Mathematica notebooks
    # (ps = Reverse[Sort[..]]) list them DESCENDING (index 0 = highest weight vector e1). Reverse so that e1 is the HW ket.
    order = sorted(rep.patterns, reverse=True)
    perm = [rep._index[p] for p in order]
    base = {}
    for (i, j) in [(1, 1), (2, 2), (3, 3), (1, 2), (2, 3), (2, 1), (3, 2)]:
        M = rep.generator(i, j)
        base[(i, j)] = SpOp.from_dense([[_frac(M[perm[a], perm[b]]) for b in range(d)] for a in range(d)])
    base[(1, 3)] = op_commutator(base[(1, 2)], base[(2, 3)])
    base[(3, 1)] = op_commutator(base[(3, 2)], base[(2, 1)])
    return rep, base


class Chain:
    """gl(3) XXX chain, L = len(thetas) sites, all in the irrep lam."""

    def __init__(self, lam=(2, 1, 0), thetas=(Fr(2, 7), Fr(-3, 5))):
        self.lam = tuple(lam)
        self.thetas = [Fr(t) for t in thetas]
        self.L = len(thetas)
        self.rep, self.E1 = gl3_generators(self.lam)
        self.d1 = self.rep.dim
        self.dim = self.d1 ** self.L
        I1 = SpOp.identity(self.d1)
        self.E = []  # E[s][(i,j)]
        for s in range(self.L):
            Es = {}
            for key, op in self.E1.items():
                M = None
                for t in range(self.L):
                    f = op if t == s else I1
                    M = f if M is None else M.kron(f)
                Es[key] = M
            self.E.append(Es)
        self.id = SpOp.identity(self.dim)

    # --- checks -------------------------------------------------------------------------------
    def check_commutation(self):
        bad = 0
        for i in range(1, 4):
            for j in range(1, 4):
                for k in range(1, 4):
                    for l in range(1, 4):
                        lhs = op_commutator(self.E1[(i, j)], self.E1[(k, l)])
                        rhs = SpOp(self.d1)
                        if j == k:
                            rhs = rhs + self.E1[(i, l)]
                        if l == i:
                            rhs = rhs - self.E1[(k, j)]
                        if not (lhs - rhs).is_zero():
                            bad += 1
        return bad

    # --- monodromy ----------------------------------------------------------------------------
    def lax(self, s, x, h, transpose=False):
        """L^(s)_ij(x) = (x - theta_s) delta_ij + h E^(s)_ij  (transpose: E_ji)."""
        th = self.thetas[s]
        out = [[None] * 3 for _ in range(3)]
        for i in range(3):
            for j in range(3):
                e = self.E[s][(j + 1, i + 1)] if transpose else self.E[s][(i + 1, j + 1)]
                op = e.scale(h)
                if i == j:
                    op = op + self.id.scale(x - th)
                out[i][j] = op
        return out

    def mono_untwisted(self, x, h, sites=None, transpose=False):
        sites = list(range(self.L)) if sites is None else sites
        T = None
        for s in sites:
            Ls = self.lax(s, x, h, transpose)
            if T is None:
                T = Ls
            else:
                T = [[sum_ops([T[i][k] * Ls[k][j] for k in range(3)]) for j in range(3)] for i in range(3)]
        return T

    def mono(self, x, h, G=None, sites=None, transpose=False):
        """T(x).G (twist on the right); G a 3x3 nested list (entries Fraction/Ser) or None."""
        T = self.mono_untwisted(x, h, sites, transpose)
        if G is None:
            return T
        out = [[None] * 3 for _ in range(3)]
        for i in range(3):
            for j in range(3):
                terms = [T[i][l].scale(G[l][j]) for l in range(3) if not _is_zero(G[l][j])]
                out[i][j] = sum_ops(terms) if terms else SpOp(T[0][0].n)
        return out

    def qmin(self, rows, cols, u, h, Tof):
        """qmin_{rows,cols}(u): rows antisymmetrised, k-th factor at u-(k-1)h. Tof(x) -> 3x3 monodromy at x."""
        a = len(cols)
        tot = None
        Ts = [Tof(u - k * h) for k in range(a)]
        for perm in permutations(range(a)):
            sgn = perm_sign(perm)
            prod = None
            for k in range(a):
                Tk = Ts[k][rows[perm[k]] - 1][cols[k] - 1]
                prod = Tk if prod is None else prod * Tk
            term = prod.scale(sgn)
            tot = term if tot is None else tot + term
        return tot

    def t_a(self, a, u, h, G=None, sites=None, transpose=False):
        cache = {}

        def Tof(x):
            if x not in cache:
                cache[x] = self.mono(x, h, G, sites, transpose)
            return cache[x]

        tot = None
        for I in combinations((1, 2, 3), a):
            q = self.qmin(I, I, u, h, Tof)
            tot = q if tot is None else tot + q
        return tot

    def nu(self, n, x, h):
        out = Fr(1) if not isinstance(h, Ser) else Ser.const(1)
        for th in self.thetas:
            out = out * (x - th + h * self.lam[n - 1])
        return out


def _is_zero(x):
    if isinstance(x, Ser):
        return all(a == 0 for a in x.c)
    return x == 0


def sum_ops(lst):
    tot = None
    for o in lst:
        tot = o if tot is None else tot + o
    return tot


def perm_sign(p):
    s, p = 1, list(p)
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]
            p[i], p[j] = p[j], p[i]
            s = -s
    return s


def _fact(k):
    out = 1
    for i in range(2, k + 1):
        out *= i
    return out


# ----------------------------------------------------------------------------------------------
# frames and vacuum covectors
# ----------------------------------------------------------------------------------------------
def gauss_r(z, w):
    """Gauss parameters r_1..r_3 of the vacuum n(z,w) (plan section 1, notebook rGauss); w = (w1, w2, ...)."""
    z1, z2, z3 = z
    w1, w2 = w[0], w[1]
    r1 = (w1 - z1) / (w1 - z2)
    r2 = (w1 - z1) * (z1 - w2) * (z2 - z3) / ((z1 - z2) * (w1 - z3) * (w2 - z3))
    r3 = (w1 - z2) * (z2 - w2) * (z1 - z3) / ((z1 - z2) * (w1 - z3) * (z3 - w2))
    return r1, r2, r3


def covec_Omega(chain, r):
    """e1^T exp(r1 E12) exp(r2 E13) exp(r3 E23) on ONE site; r entries Fraction or Ser."""
    d = chain.d1
    v = [Fr(0)] * d
    v[0] = Fr(1)
    nil = chain.lam[0] - chain.lam[2]
    for (key, c) in (((1, 2), r[0]), ((1, 3), r[1]), ((2, 3), r[2])):
        E = chain.E1[key]
        term = list(v)
        acc = list(v)
        for k in range(1, nil + 1):
            term = E.vec_times(term)
            fac = (c ** k) * Fr(1, _fact(k))
            acc = [a + fac * t for a, t in zip(acc, term)]
        v = acc
    return v


def kron_vec(a, b):
    return [x * y for x in a for y in b]
