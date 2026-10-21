"""P2: the XXX -> Gaudin limit for gl(3), L = 2, Lam = (2,1,0)  (tasks X1-X4 of SPEC-xxx-limit.md).

Run top to bottom:   python -X utf8 xxx_limit.py        (writes nothing; prints; RESULTS-xxx-limit.md is the written record)
Exactness policy (stated once):
  * X2/X3 operator coefficients are EXACT: the whole XXX chain is evaluated with h = the formal series variable
    (class Ser, truncated at h^16, Fraction coefficients) and z = exp(h zeta) the exact series, so every h-expansion
    coefficient is an exact rational (no Richardson needed).  Cross-checks at the three rational h = 1e-2, 1e-3, 1e-4
    use the exponential series truncated at h^15 (rational arithmetic), i.e. z = exp(h zeta) up to O(h^16).
  * X1 uses mpmath at 80 digits with exact exponentials (for rates) and h = 1e-30 (for the limits).
  * X4 is sympy-symbolic in h (scalars only).
"""
import sys
import time
from fractions import Fraction as Fr
from itertools import permutations
from math import comb, factorial

sys.path.insert(0, __file__.rsplit("\\", 1)[0] if "\\" in __file__ else __file__.rsplit("/", 1)[0])
from xxx_chain import (NSER, Chain, Ser, SpOp, gauss_r, covec_Omega, kron_vec, ratio_h, perm_sign, sum_ops,  # noqa: E402
                       op_commutator, _is_zero)

ZETA = [Fr(3, 2), Fr(5, 4), Fr(2, 3)]
OMEGA = [Fr(5, 2), Fr(7, 3), Fr(3, 5)]
THETA = [Fr(2, 7), Fr(-3, 5)]
HS = [Fr(1, 100), Fr(1, 1000), Fr(1, 10000)]
LAM = (2, 1, 0)
REPORT = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    REPORT.append(s)


def fs(x):
    """compact string of a Fraction / list / nested list."""
    if isinstance(x, (list, tuple)):
        return "[" + ", ".join(fs(y) for y in x) + "]"
    return str(x)


# ==============================================================================================
# series / operator helpers
# ==============================================================================================
def coeff_op(op, k):
    """h^k coefficient of an operator with Ser entries -> SpOp with Fraction entries."""
    out = {}
    for i, row in op.r.items():
        d = {}
        for j, v in row.items():
            c = v.c[k] if isinstance(v, Ser) else (v if k == 0 else Fr(0))
            if c != 0:
                d[j] = c
        if d:
            out[i] = d
    return SpOp(op.n, out)


def maxabs(op):
    m = Fr(0)
    for row in op.r.values():
        for v in row.values():
            if abs(v) > m:
                m = abs(v)
    return m


def fl(x):
    return float(x)


def diag3(z):
    return [[z[0], 0, 0], [0, z[1], 0], [0, 0, z[2]]]


# ==============================================================================================
# minimal Gaudin model (independent of P1's gaudin_chain.py)
#   L_g(u)_ij = g_ij 1 + sum_s E^(s)_ij/(u-theta_s)            (transpose=True: E_ji)
#   cdet(d - L_g) = sum_sigma sgn(sigma) M_{s1,1} M_{s2,2} M_{s3,3}   (column-ordered; row_ordered=True: M_{1,s1}...)
#   = d^3 - S1 d^2 + S2 d - S3,  d moved to the right (d f = f d + f')
# differential operators with SpOp-valued coefficients are stored as {n: jets}, jets = [c, c', c'', ...]
# ==============================================================================================
JET = 5


def _jmul(a, b):
    J = min(len(a), len(b))
    out = []
    for j in range(J):
        acc = None
        for i in range(j + 1):
            t = (a[i] * b[j - i]).scale(comb(j, i))
            acc = t if acc is None else acc + t
        out.append(acc)
    return out


def _jadd(a, b):
    J = min(len(a), len(b))
    return [a[i] + b[i] for i in range(J)]


def do_mul(A, B):
    out = {}
    for m, a in A.items():
        for n, b in B.items():
            for k in range(m + 1):
                bk = b[k:]
                if not bk:
                    continue
                prod = [x.scale(comb(m, k)) for x in _jmul(a, bk)]
                key = m + n - k
                out[key] = prod if key not in out else _jadd(out[key], prod)
    return out


def do_add(A, B, sgn=1):
    out = {k: [x for x in v] for k, v in A.items()}
    for k, v in B.items():
        v = [x.scale(sgn) for x in v]
        out[k] = v if k not in out else _jadd(out[k], v)
    return out


def gaudin_L_jets(ch, u, g, transpose=False):
    """L_ij jets [L, L', L'', ...] (JET entries) as SpOp, i,j = 0..2."""
    L = [[None] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            jets = []
            for k in range(JET):
                acc = ch.id.scale(g[i][j]) if (k == 0 and not _is_zero(g[i][j])) else SpOp(ch.dim)
                for s, th in enumerate(ch.thetas):
                    e = ch.E[s][(j + 1, i + 1)] if transpose else ch.E[s][(i + 1, j + 1)]
                    acc = acc + e.scale(Fr((-1) ** k * factorial(k)) / (u - th) ** (k + 1))
                jets.append(acc)
            L[i][j] = jets
    return L


def gaudin_cdet(ch, u, g, transpose=False, row_ordered=False):
    """returns (S1, S2, S3) as SpOp; S1 should be scalar."""
    Ljets = gaudin_L_jets(ch, u, g, transpose)
    zero = SpOp(ch.dim)
    Z = [zero] * JET
    M = [[None] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            d1 = [ch.id] + [zero] * (JET - 1) if i == j else Z
            M[i][j] = {1: d1, 0: [x.scale(-1) for x in Ljets[i][j]]}
    tot = {}
    for p in permutations(range(3)):
        sgn = perm_sign(p)
        if row_ordered:
            facs = [M[k][p[k]] for k in range(3)]
        else:
            facs = [M[p[k]][k] for k in range(3)]
        prod = facs[0]
        for f in facs[1:]:
            prod = do_mul(prod, f)
        prod = {k: [x.scale(sgn) for x in v] for k, v in prod.items()}
        tot = do_add(tot, prod)
    assert (tot[3][0] - ch.id).is_zero()
    S1 = tot[2][0].scale(-1)
    S2 = tot[1][0]
    S3 = tot[0][0].scale(-1)
    return S1, S2, S3


def is_scalar_op(op, n):
    f = op.flat()
    if any(i != j for (i, j) in f):
        return None
    vals = {f.get((i, i), Fr(0)) for i in range(n)}
    return vals.pop() if len(vals) == 1 else None


# ==============================================================================================
# generic exact linear algebra over Q for operators (express target in span of basis)
# ==============================================================================================
def express(target, basis, names=None):
    """Solve target = sum_k c_k basis_k exactly (Fractions). Returns (coeffs or None, residual flag)."""
    keys = set(target.flat())
    for b in basis:
        keys |= set(b.flat())
    keys = sorted(keys)
    rows = []
    cols = [b.flat() for b in basis]
    tf = target.flat()
    nb = len(basis)
    for key in keys:
        rows.append([c.get(key, Fr(0)) for c in cols] + [tf.get(key, Fr(0))])
    # gaussian elimination
    piv = []
    r = 0
    for c in range(nb):
        p = None
        for i in range(r, len(rows)):
            if rows[i][c] != 0:
                p = i
                break
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        pv = rows[r][c]
        rows[r] = [x / pv for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                f = rows[i][c]
                rows[i] = [x - f * y for x, y in zip(rows[i], rows[r])]
        piv.append(c)
        r += 1
    inconsistent = any(all(x == 0 for x in rows[i][:nb]) and rows[i][nb] != 0 for i in range(r, len(rows)))
    if inconsistent:
        return None, piv
    sol = [Fr(0)] * nb
    for i, c in enumerate(piv):
        sol[c] = rows[i][nb]
    return sol, piv


# ==============================================================================================
# small exact 3x3 / vector helpers
# ==============================================================================================
def mat3_mul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def mat3_inv(A):
    a, b, c = A[0]
    d, e, f = A[1]
    g, hh, i = A[2]
    det = a * (e * i - f * hh) - b * (d * i - f * g) + c * (d * hh - e * g)
    cof = [[e * i - f * hh, -(d * i - f * g), d * hh - e * g],
           [-(b * i - c * hh), a * i - c * g, -(a * hh - b * g)],
           [b * f - c * e, -(a * f - c * d), a * e - b * d]]
    return [[cof[j][i_] / det for j in range(3)] for i_ in range(3)]  # adjugate / det


def mat3_T(A):
    return [[A[j][i] for j in range(3)] for i in range(3)]


def esym(t):
    return t[0] + t[1] + t[2], t[0] * t[1] + t[0] * t[2] + t[1] * t[2], t[0] * t[1] * t[2]


def Paux_exact(z):
    e1, e2, _ = esym(z)
    return [[Fr(1), zi - e1, zi * zi - e1 * zi + e2] for zi in z]


def Gc_exact(tw):
    e1, e2, e3 = esym(tw)
    return [[e1, -e2, e3], [Fr(1), Fr(0), Fr(0)], [Fr(0), Fr(1), Fr(0)]]


def gprime(zeta, omega):
    """g' = P_aux(zeta) g_c(omega) P_aux(zeta)^-1  (Lie-algebra companion frame, exact)."""
    P = Paux_exact(zeta)
    return mat3_mul(mat3_mul(P, Gc_exact(omega)), mat3_inv(P))


def vec_scale(v, c):
    return [c * x for x in v]


def vec_add(a, b):
    return [x + y for x, y in zip(a, b)]


def vec_sub(a, b):
    return [x - y for x, y in zip(a, b)]


def vec_dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def vec_norm2(a):
    return vec_dot(a, a)


def ls_fit(basis, target):
    """exact least squares target ~ sum c_k basis_k (basis: list of vectors; dependent ones dropped).
    returns (coeff dict {index: Fraction}, residual vector, indices used)."""
    # independent subset by greedy Gram-Schmidt in exact arithmetic
    ortho, used = [], []
    for idx, b in enumerate(basis):
        r = list(b)
        for q in ortho:
            r = vec_sub(r, vec_scale(q, vec_dot(r, q) / vec_dot(q, q)))
        if vec_norm2(r) != 0:
            ortho.append(r)
            used.append(idx)
    # solve normal equations on used columns
    B = [basis[i] for i in used]
    n = len(B)
    M = [[vec_dot(B[i], B[j]) for j in range(n)] + [vec_dot(B[i], target)] for i in range(n)]
    for c in range(n):
        p = next(i for i in range(c, n) if M[i][c] != 0)
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for i in range(n):
            if i != c and M[i][c] != 0:
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[c])]
    coef = {used[i]: M[i][n] for i in range(n)}
    res = list(target)
    for i in range(n):
        res = vec_sub(res, vec_scale(B[i], M[i][n]))
    return coef, res, used


def matrix_rank_exact(M):
    M = [list(r) for r in M]
    rk, ncol = 0, len(M[0])
    for c in range(ncol):
        p = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
        if p is None:
            continue
        M[rk], M[p] = M[p], M[rk]
        for i in range(rk + 1, len(M)):
            if M[i][c] != 0:
                f = M[i][c] / M[rk][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[rk])]
        rk += 1
    return rk


def fvec(v, nd=6):
    return "[" + ", ".join(("%.*g" % (nd, float(x))) for x in v) + "]"


def svec(v):
    return "[" + ", ".join(str(x) for x in v) + "]"


def laurent_fit(ops_at, us, kmax, thetas):
    """operator-valued rational function f(u) = c0 + sum_s sum_{k=1..kmax} c_{s,k}/(u-theta_s)^k from samples.
    ops_at: list of SpOp at points us (len >= 1+len(thetas)*kmax+1, last one used as a check).
    returns dict {(s,k): SpOp} with (None,0) = constant, and a bool for the consistency check."""
    nunk = 1 + len(thetas) * kmax
    keys = [(None, 0)] + [(s, k) for s in range(len(thetas)) for k in range(1, kmax + 1)]
    rows = []
    for u in us[:nunk]:
        rows.append([Fr(1)] + [1 / (u - thetas[s]) ** k for (s, k) in keys[1:]])
    # invert rows (nunk x nunk) exactly
    n = nunk
    M = [rows[i] + [Fr(1) if i == j else Fr(0) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(i for i in range(c, n) if M[i][c] != 0)
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for i in range(n):
            if i != c and M[i][c] != 0:
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[c])]
    Minv = [row[n:] for row in M]
    out = {}
    for ki, key in enumerate(keys):
        acc = SpOp(ops_at[0].n)
        for j in range(n):
            if Minv[ki][j] != 0:
                acc = acc + ops_at[j].scale(Minv[ki][j])
        out[key] = acc
    # check at the extra points
    ok = True
    for u, op in zip(us[nunk:], ops_at[nunk:]):
        pred = out[(None, 0)]
        for (s, k) in keys[1:]:
            pred = pred + out[(s, k)].scale(1 / (u - thetas[s]) ** k)
        ok = ok and (pred - op).is_zero()
    return out, ok


# ==============================================================================================
# small exact 3x3 / vector helpers
# ==============================================================================================
def mat3_mul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def mat3_inv(A):
    a, b, c = A[0]
    d, e, f = A[1]
    g, hh, i = A[2]
    det = a * (e * i - f * hh) - b * (d * i - f * g) + c * (d * hh - e * g)
    cof = [[e * i - f * hh, -(d * i - f * g), d * hh - e * g],
           [-(b * i - c * hh), a * i - c * g, -(a * hh - b * g)],
           [b * f - c * e, -(a * f - c * d), a * e - b * d]]
    return [[cof[j][i_] / det for j in range(3)] for i_ in range(3)]  # adjugate / det


def mat3_T(A):
    return [[A[j][i] for j in range(3)] for i in range(3)]


def esym(t):
    return t[0] + t[1] + t[2], t[0] * t[1] + t[0] * t[2] + t[1] * t[2], t[0] * t[1] * t[2]


def Paux_exact(z):
    e1, e2, _ = esym(z)
    return [[Fr(1), zi - e1, zi * zi - e1 * zi + e2] for zi in z]


def Gc_exact(tw):
    e1, e2, e3 = esym(tw)
    return [[e1, -e2, e3], [Fr(1), Fr(0), Fr(0)], [Fr(0), Fr(1), Fr(0)]]


def gprime(zeta, omega):
    """g' = P_aux(zeta) g_c(omega) P_aux(zeta)^-1  (Lie-algebra companion frame, exact)."""
    P = Paux_exact(zeta)
    return mat3_mul(mat3_mul(P, Gc_exact(omega)), mat3_inv(P))


def vec_scale(v, c):
    return [c * x for x in v]


def vec_add(a, b):
    return [x + y for x, y in zip(a, b)]


def vec_sub(a, b):
    return [x - y for x, y in zip(a, b)]


def vec_dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def vec_norm2(a):
    return vec_dot(a, a)


def ls_fit(basis, target):
    """exact least squares target ~ sum c_k basis_k (dependent basis vectors dropped).
    returns (coeff dict {index: Fraction}, residual vector, used indices)."""
    ortho, used = [], []
    for idx, b in enumerate(basis):
        r = list(b)
        for q in ortho:
            r = vec_sub(r, vec_scale(q, vec_dot(r, q) / vec_dot(q, q)))
        if vec_norm2(r) != 0:
            ortho.append(r)
            used.append(idx)
    B = [basis[i] for i in used]
    n = len(B)
    M = [[vec_dot(B[i], B[j]) for j in range(n)] + [vec_dot(B[i], target)] for i in range(n)]
    for c in range(n):
        p = next(i for i in range(c, n) if M[i][c] != 0)
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for i in range(n):
            if i != c and M[i][c] != 0:
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[c])]
    coef = {used[i]: M[i][n] for i in range(n)}
    res = list(target)
    for i in range(n):
        res = vec_sub(res, vec_scale(B[i], M[i][n]))
    return coef, res, used


def fvec(v, nd=6):
    return "[" + ", ".join(("%.*g" % (nd, float(x))) for x in v) + "]"


def svec(v):
    return "[" + ", ".join(str(x) for x in v) + "]"


def laurent_fit(ops_at, us, kmax, thetas):
    """operator-valued f(u) = c0 + sum_s sum_{k=1..kmax} c_{s,k}/(u-theta_s)^k from samples (extra samples = check).
    returns ({(None,0): c0, (s,k): c_{s,k}}, consistent?)"""
    nunk = 1 + len(thetas) * kmax
    keys = [(None, 0)] + [(s, k) for s in range(len(thetas)) for k in range(1, kmax + 1)]
    rows = []
    for u in us[:nunk]:
        rows.append([Fr(1)] + [1 / (u - thetas[s]) ** k for (s, k) in keys[1:]])
    n = nunk
    M = [rows[i] + [Fr(1) if i == j else Fr(0) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(i for i in range(c, n) if M[i][c] != 0)
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for i in range(n):
            if i != c and M[i][c] != 0:
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[c])]
    Minv = [row[n:] for row in M]
    out = {}
    for ki, key in enumerate(keys):
        acc = SpOp(ops_at[0].n)
        for j in range(n):
            if Minv[ki][j] != 0:
                acc = acc + ops_at[j].scale(Minv[ki][j])
        out[key] = acc
    ok = True
    for u, op in zip(us[nunk:], ops_at[nunk:]):
        pred = out[(None, 0)]
        for (s, k) in keys[1:]:
            pred = pred + out[(s, k)].scale(1 / (u - thetas[s]) ** k)
        ok = ok and (pred - op).is_zero()
    return out, ok


# ==============================================================================================
# X1  frames and the vacuum
# ==============================================================================================
def gauss_r_series(zeta, omega):
    """exact h-series of r_k(z,w) at z = exp(h zeta), w = exp(h omega) (all differences divided by h first)."""
    z = [Ser.exp(c) for c in zeta]
    w = [Ser.exp(c) for c in omega]

    def d(a, b):
        return (a - b).divh()

    r1 = d(w[0], z[0]) / d(w[0], z[1])
    r2 = d(w[0], z[0]) * d(z[0], w[1]) * d(z[1], z[2]) / (d(z[0], z[1]) * d(w[0], z[2]) * d(w[1], z[2]))
    r3 = d(w[0], z[1]) * d(z[1], w[1]) * d(z[0], z[2]) / (d(z[0], z[1]) * d(w[0], z[2]) * d(z[2], w[1]))
    return [r1, r2, r3]


def _mpf(x):
    import mpmath as mp
    return mp.mpf(x.numerator) / x.denominator


def mp_from_op(op):
    import mpmath as mp
    n = op.n
    M = mp.zeros(n, n)
    for (i, j), v in op.flat().items():
        M[i, j] = _mpf(v)
    return M


def mp_Paux_nb(z):
    """notebook definition of P_aux (product of exponentials), mpmath."""
    import mpmath as mp

    def e3(i, j):
        M = mp.zeros(3, 3)
        M[i - 1, j - 1] = 1
        return M

    def ex(c, i, j):
        return mp.eye(3) + c * e3(i, j)

    s = [z[0] ** 2 / ((z[0] - z[1]) * (z[0] - z[2])), -(z[1] / (z[0] * (z[1] - z[2]))), 1 / (z[0] * z[1])]
    rh = [-(z[2] / z[1]), -(z[2] / z[0]), (z[1] ** 2 * (z[2] - z[0])) / (z[0] ** 2 * (z[1] - z[2]))]
    D = mp.diag([1 / s[0], 1 / s[1], 1 / s[2]])
    return (ex(-rh[0], 2, 3) * ex(-rh[1], 1, 3) * ex(-rh[2], 1, 2) * D * ex(-1 / z[0], 2, 1)
            * ex(1 / (z[0] * z[1]), 3, 1) * ex(-(z[0] + z[1]) / (z[0] * z[1]), 3, 2))


def mp_Gc(tw):
    import mpmath as mp
    e1, e2, e3 = esym(tw)
    return mp.matrix([[e1, -e2, e3], [1, 0, 0], [0, 1, 0]])


def mp_x0_onesite(chain, z):
    """one-site x = e1^T P_phys(z)^-1 via the factorised inverse (notebook x0Of): e1 . ppinv . dinv . pminv."""
    import mpmath as mp
    E = {k: mp_from_op(v) for k, v in chain.E1.items()}
    d = chain.d1
    nil = chain.lam[0] - chain.lam[2]

    def expN(c, key):
        A = E[key] * c
        out = mp.eye(d)
        P = mp.eye(d)
        for k in range(1, nil + 1):
            P = P * A
            out = out + P / mp.factorial(k)
        return out

    z1, z2, z3 = z
    s = [z1 ** 2 / ((z1 - z2) * (z1 - z3)), -(z2 / (z1 * (z2 - z3))), 1 / (z1 * z2)]
    r = [-(z3 / z2), -(z3 / z1), (z2 ** 2 * (z3 - z1)) / (z1 ** 2 * (z2 - z3))]
    ppinv = expN(-(z1 + z2) / (z1 * z2), (2, 3)) * expN(1 / (z1 * z2), (1, 3)) * expN(-(1 / z1), (1, 2))
    wts = [[E[(k, k)][p, p] for k in (1, 2, 3)] for p in range(d)]
    dinv = mp.diag([mp.fprod([s[k] ** (-wts[p][k]) for k in range(3)]) for p in range(d)])
    pminv = expN(-r[2], (2, 1)) * expN(-r[1], (3, 1)) * expN(-r[0], (3, 2))
    e1 = mp.zeros(1, d)
    e1[0, 0] = 1
    return e1 * ppinv * dinv * pminv


def run_X1(chain1):
    import mpmath as mp
    mp.mp.dps = 80
    say("=" * 100)
    say("X1  frames and the vacuum   (mpmath, 80 digits, exact exponentials z = exp(h zeta), w = exp(h omega))")
    say("=" * 100)
    zeta, omega = ZETA, OMEGA
    zm = [_mpf(c) for c in zeta]
    om = [_mpf(c) for c in omega]

    # ---- (a) Gauss parameters -----------------------------------------------------------------
    say("\n(a) Gauss parameters r_k(z,w) -> r_k(zeta,omega)")
    r0e = gauss_r(zeta, omega)
    r0 = [_mpf(x) for x in r0e]
    say("  r_k(zeta,omega) (exact) =", fs(list(r0e)), "=", [mp.nstr(x, 12) for x in r0])
    prev = None
    say("  h        |r_k(z,w)-r_k(zeta,omega)|, k=1,2,3                   ; ratio to previous h (O(h) <=> ratio ~ 10)")
    for h in HS:
        hm = _mpf(h)
        z = [mp.exp(hm * c) for c in zm]
        w = [mp.exp(hm * c) for c in om]
        r = gauss_r(z, w)
        dev = [abs(r[k] - r0[k]) for k in range(3)]
        rate = "" if prev is None else "   ratios " + "  ".join(mp.nstr(prev[k] / dev[k], 5) for k in range(3))
        say("  %-8s %s%s" % (str(h), "  ".join(mp.nstr(x, 6) for x in dev), rate))
        prev = dev
    rs = gauss_r_series(zeta, omega)
    say("  exact series (Fractions):  r_k(h) = r_k(0) + h r_k^(1) + h^2 r_k^(2) + ...")
    for k in range(3):
        say("    r_%d = %s + h*(%s) + h^2*(%s) + h^3*(%s) + ..." % (k + 1, rs[k].c[0], rs[k].c[1], rs[k].c[2], rs[k].c[3]))
    say("  check r_k(h=0) == r_k(zeta,omega) exactly:", all(rs[k].c[0] == r0e[k] for k in range(3)))

    # ---- (b) frames -----------------------------------------------------------------------------
    say("\n(b) G'(z,w) = P_aux(z) G_c(w) P_aux(z)^-1 ; limit of (G'-1)/h")
    Pnb = mp_Paux_nb(zm)
    Prow = mp.matrix([[_mpf(x) for x in row] for row in Paux_exact(zeta)])
    say("  notebook P_aux (exponential product)  vs  rows (1, z_i-chi1, z_i^2-chi1 z_i+chi2): max |diff| =",
        mp.nstr(mp.mpf(max(abs(Pnb[i, j] - Prow[i, j]) for i in range(3) for j in range(3))), 5))
    dd = Prow * mp_Gc(zm) * Prow ** -1 - mp.diag(zm)
    say("  P_aux G_c(z) P_aux^-1 - diag(z), max abs:", mp.nstr(mp.mpf(max(abs(dd[i, j]) for i in range(3) for j in range(3))), 5))
    g1 = gprime(zeta, omega)
    say("  g' = P_aux(zeta) g_c(omega) P_aux(zeta)^-1 (exact) =", fs(g1))
    G1m = mp.matrix([[_mpf(x) for x in row] for row in g1])
    PG = Prow
    variants = {"g'=P gc P^-1": G1m, "(g')^T": G1m.T, "(g')^-1": G1m ** -1,
                "P^-1 gc P": PG ** -1 * mp_Gc(om) * PG, "-g'": -G1m}
    names = list(variants)
    say("  max-entry distance ||(G'-1)/h - X|| for X in:   " + " | ".join(names))
    hlist = [(str(HS[0]), _mpf(HS[0])), (str(HS[1]), _mpf(HS[1])), (str(HS[2]), _mpf(HS[2])), ("1e-30", mp.mpf("1e-30"))]
    lim = None
    for hstr, hm in hlist:
        mp.mp.dps = 250 if hstr == "1e-30" else 80   # G' = P G_c P^-1 loses ~3 log10(1/h) digits: 90 digits at h = 1e-30
        z = [mp.exp(hm * c) for c in zm]
        w = [mp.exp(hm * c) for c in om]
        c1 = z[0] + z[1] + z[2]
        c2 = z[0] * z[1] + z[0] * z[2] + z[1] * z[2]
        P = mp.matrix([[1, zi - c1, zi ** 2 - c1 * zi + c2] for zi in z])
        G = P * mp_Gc(w) * P ** -1
        D = (G - mp.eye(3)) / hm
        row = [mp.nstr(mp.mpf(max(abs(D[i, j] - variants[nm][i, j]) for i in range(3) for j in range(3))), 5) for nm in names]
        say("  h=%-8s %s" % (hstr, " | ".join(row)))
        lim = D
    mp.mp.dps = 80
    say("  (G'-1)/h at h=1e-30 =", [[mp.nstr(lim[i, j], 12) for j in range(3)] for i in range(3)])
    say("  => (G'-1)/h CONVERGES, error O(h) (distance drops ~10x per decade), to g' = P_aux(zeta) g_c(omega) P_aux(zeta)^-1 itself "
        "(no transpose, no inverse).  G_c(1,1,1) is a Jordan block but G' = 1 + h g' + O(h^2) anyway.")

    # ---- (c) x0 -------------------------------------------------------------------------------
    say("\n(c) x0(z) = e1^T P_phys(z)^-1 (one site; the two-site x0 is x (x) x so it is parallel iff the one-site vectors are)")
    xs = {}
    for hstr in ["1e-10", "1e-20", "1e-30"]:
        mp.mp.dps = 200
        hm = mp.mpf(hstr)
        z = [mp.exp(hm * _mpf(c)) for c in zeta]
        xs[hstr] = mp_x0_onesite(chain1, z)
    mp.mp.dps = 80
    n1, n2 = mp.norm(xs["1e-20"], 2), mp.norm(xs["1e-30"], 2)
    p_lead = mp.log(n2 / n1) / mp.log(mp.mpf("1e-30") / mp.mpf("1e-20"))
    say("  leading h-power of ||x0(e^{h zeta})|| (from h=1e-20 and 1e-30):  ||x0|| ~ h^(%s)" % mp.nstr(p_lead, 8))
    xv = xs["1e-30"]
    k0 = max(range(chain1.d1), key=lambda j: abs(xv[0, j]))
    dirx = [xv[0, j] / xv[0, k0] for j in range(chain1.d1)]
    xG = mp_x0_onesite(chain1, zm)
    dirG = [xG[0, j] / xG[0, k0] for j in range(chain1.d1)]
    say("  XXX x0 direction as h->0 (projective, component %d = 1):" % k0, [mp.nstr(c, 10) for c in dirx])
    say("  Gaudin x0^G(zeta) direction (same normalisation)       :", [mp.nstr(c, 10) for c in dirG])
    dist = max(abs(dirx[j] - dirG[j]) for j in range(chain1.d1))
    say("  projective distance max_j |dir_XXX - dir_G| = %s  -> %s" % (mp.nstr(dist, 6), "PARALLEL" if dist < mp.mpf("1e-10") else "NOT PARALLEL"))
    dx2 = xs["1e-20"]
    dir2 = [dx2[0, j] / dx2[0, k0] for j in range(chain1.d1)]
    say("  XXX direction stable between h=1e-20 and 1e-30: max diff = %s" % mp.nstr(max(abs(dir2[j] - dirx[j]) for j in range(chain1.d1)), 5))
    return dict(g1=g1, p_lead=p_lead, dirx=dirx, dirG=dirG)


# ==============================================================================================
# port checks (XXX chain) and the Gaudin bra-vacuum eigencovector condition
# ==============================================================================================
def run_checks(ch):
    say("=" * 100)
    say("PORT CHECKS (XXX chain at the notebook point h = 1/2, theta = (2/7,-3/5), z = (3/2,5/4,2/3), w = (5/2,7/3,3/5))")
    say("=" * 100)
    say("  basis: gl3_reps.py lists the GT patterns in a different order (index 0 = LOWEST weight (0,1,2)); xxx_chain re-orders them")
    say("         to the notebooks' order (patterns sorted descending), so basis vector 0 is the HW ket e1 (E12 e1 = E13 e1 = E23 e1 = 0).")
    e1 = [Fr(1)] + [Fr(0)] * (ch.d1 - 1)
    hw = all(all(x == 0 for x in ch.E1[k].times_vec(e1)) for k in [(1, 2), (1, 3), (2, 3)])
    say("  HW ket e1 annihilated by E12, E13, E23:", hw, "; weights of the 8 basis vectors (E11,E22,E33):",
        [tuple(int(ch.E1[(k, k)].r.get(p, {}).get(p, 0)) for k in (1, 2, 3)) for p in range(ch.d1)])
    say("  gl(3) commutation relations [E_ij,E_kl] = d_jk E_il - d_li E_kj violated in", ch.check_commutation(), "of 81 cases")
    h = Fr(1, 2)
    z = [Fr(3, 2), Fr(5, 4), Fr(2, 3)]
    w = [Fr(5, 2), Fr(7, 3), Fr(3, 5)]
    G = diag3(z)
    t0 = time.time()
    u = Fr(1, 5)
    t3 = ch.t_a(3, u, h, G)
    qdet = z[0] * z[1] * z[2] * ch.nu(1, u, h) * ch.nu(2, u - h, h) * ch.nu(3, u - 2 * h, h)
    say("  quantum determinant t3(u) = chi3 nu1(u) nu2(u-h) nu3(u-2h) * 1 :", (t3 - ch.id.scale(qdet)).is_zero(),
        "   (one 64x64 exact quantum minor t3: %.3f s)" % (time.time() - t0))
    t1 = ch.t_a(1, u, h, G)
    t2 = ch.t_a(2, Fr(2, 3), h, G)
    say("  [t1(1/5), t2(2/3)] = 0 :", op_commutator(t1, t2).is_zero())
    # bra vacuum eigencovector with the bra twist G' = P_aux(z) G_c(w) P_aux(z)^-1
    Gp = mat3_mul(mat3_mul(Paux_exact(z), Gc_exact(w)), mat3_inv(Paux_exact(z)))
    Om1 = covec_Omega(ch_one(ch), gauss_r(z, w))
    O2 = kron_vec(Om1, Om1)
    u = Fr(2, 9)
    res = {}
    for nm, GG in (("G'(w) = P_aux G_c(w) P_aux^-1", Gp), ("control diag(w)", diag3(w))):
        out = []
        for a in (1, 2, 3):
            ta = ch.t_a(a, u, h, GG)
            c = ta.vec_times(O2) if False else vec_times_op(ta, O2)
            out.append(all(x == 0 for x in vec_sub(c, vec_scale(O2, c[0] / O2[0]))))
        res[nm] = out
    say("  Omega_zw (x) Omega_zw is a LEFT eigencovector of t_a(u=2/9), a=1,2,3, twist G' :", res["G'(w) = P_aux G_c(w) P_aux^-1"],
        "; control twist diag(w):", res["control diag(w)"])


def ch_one(ch):
    """a one-site view of the chain object (only E1, lam, d1 used by covec_Omega)."""
    return ch


def vec_times_op(op, v):
    return op.vec_times(v)


def run_gaudin_eigencovector(ch):
    say("=" * 100)
    say("GAUDIN: bra-vacuum eigencovector condition  Omega_{zeta omega} (x) Omega_{zeta omega} . S_k(u; g') = lambda . Omega (x) Omega")
    say("=" * 100)
    zeta, omega = ZETA, OMEGA
    Om1 = covec_Omega(ch, gauss_r(zeta, omega))
    O2 = kron_vec(Om1, Om1)
    g1 = gprime(zeta, omega)
    P = Paux_exact(zeta)
    variants = {
        "g' = P gc P^-1": g1,
        "(g')^T": mat3_T(g1),
        "(g')^-1": mat3_inv(g1),
        "P^-1 gc P": mat3_mul(mat3_mul(mat3_inv(P), Gc_exact(omega)), P),
        "diag(omega) (control)": diag3(omega),
    }
    u = Fr(2, 9)
    for tr, ro in ((False, False), (True, False), (False, True), (True, True)):
        row = []
        for nm, gg in variants.items():
            S1, S2, S3 = gaudin_cdet(ch, u, gg, transpose=tr, row_ordered=ro)
            ok = []
            for S in (S2, S3):
                c = S.vec_times(O2)
                ok.append(all(x == 0 for x in vec_sub(c, vec_scale(O2, c[0] / O2[0]))))
            row.append("%s: S2 %s S3 %s" % (nm, ok[0], ok[1]))
        say("  Lax E_%s, cdet %s-ordered :: %s" % ("ji" if tr else "ij", "row" if ro else "column", " | ".join(row)))


# ==============================================================================================
# X2  the Hamiltonian dictionary
# ==============================================================================================
def S1_of(u, s=None):
    s = sum(ZETA) if s is None else s
    return s + 3 * sum(1 / (u - th) for th in THETA)


def S1p_of(u):
    return -3 * sum(1 / (u - th) ** 2 for th in THETA)


def S1pp_of(u):
    return 6 * sum(1 / (u - th) ** 3 for th in THETA)


def truncated_exp(a, h):
    """exp(a h) truncated at order NSER-1 (rational) -- used only for the finite-h cross-checks."""
    out, term = Fr(0), Fr(1)
    for k in range(NSER):
        out += term
        term = term * a * h / (k + 1)
    return out


def run_X2(ch):
    say("=" * 100)
    say("X2  the Hamiltonian dictionary: cdet(d - L_g) coefficients vs h-expansion of t1/nu, t2/(nu nu^[-1])")
    say("=" * 100)
    s = sum(ZETA)
    hS = Ser.h()
    G = diag3([Ser.exp(c) for c in ZETA])
    g = diag3(ZETA)
    US = [Fr(3, 7), Fr(-2, 9), Fr(5, 3), Fr(-7, 4)]
    say("  exact h-series: h is the formal variable (Ser, order %d), z = exp(h zeta); the h-coefficients below are exact rationals." % NSER)
    say("  s = sum zeta =", s, "; theta =", fs(THETA), "; zeta =", fs(ZETA))
    store = {}
    for u in US:
        t0 = time.time()
        nu0 = (u - THETA[0]) * (u - THETA[1])
        t1 = ch.t_a(1, u, hS, G)
        t2 = ch.t_a(2, u, hS, G)
        nuh = (Ser.const(u) - hS - THETA[0]) * (Ser.const(u) - hS - THETA[1])
        t2n = t2.scale((nuh * nu0).inv())
        A = [coeff_op(t1, k).scale(1 / nu0) for k in range(6)]
        B = [coeff_op(t2n, k) for k in range(6)]
        S1, S2, S3 = gaudin_cdet(ch, u, g)
        S1v = is_scalar_op(S1, ch.dim)
        store[u] = dict(A=A, B=B, S=(S1, S2, S3), S1v=S1v)
        say("\n  --- u = %s   (series build %.1f s)" % (u, time.time() - t0))
        say("  A0 = %s*1 ; B0 = %s*1 (scalars: %s, %s)" % (is_scalar_op(A[0], ch.dim), is_scalar_op(B[0], ch.dim),
                                                          is_scalar_op(A[0], ch.dim) is not None, is_scalar_op(B[0], ch.dim) is not None))
        a1 = is_scalar_op(A[1], ch.dim)
        b1 = is_scalar_op(B[1], ch.dim)
        say("  A1 = %s*1 (scalar), S1 = %s*1 (scalar), S1 formula s + 3 sum 1/(u-theta_s) = %s  ->  A1 == S1: %s" %
            (a1, S1v, S1_of(u), a1 == S1v == S1_of(u)))
        say("  B1 = %s*1 (scalar) ; B1 == 2 A1: %s" % (b1, b1 == 2 * a1))
        S1u, S1p, S1pp = S1_of(u), S1p_of(u), S1pp_of(u)
        c2 = S1u ** 2 / 2 + S1p / 2
        c3 = S1u ** 3 / 6 + S1u * S1p / 2 + S1pp / 2
        r2 = (S2 + A[2] - ch.id.scale(c2))
        r3 = (S3 - A[3] + B[3] - ch.id.scale(c3))
        say("  S2 = -A2 + c2(u)*1 with c2 = S1^2/2 + S1'/2 = %s : residual zero: %s" % (c2, r2.is_zero()))
        say("  S3 = A3 - B3 + c3(u)*1 with c3 = S1^3/6 + S1 S1'/2 + S1''/2 = %s : residual zero: %s" % (c3, r3.is_zero()))
        nup = 2 * u - THETA[0] - THETA[1]
        nuu = (u - THETA[0]) * (u - THETA[1])
        b2c = s * s / 2 + 3 * s * nup / nuu + 9 / nuu
        r4 = B[2] - A[2] + ch.id.scale(S1p - b2c)
        say("  B2 = A2 - S1' + (s^2/2 + 3 s nu'/nu + 9/nu) = A2 + %s : residual zero: %s" % (-S1p + b2c, r4.is_zero()))
        # generic fit (exact linear algebra): no ansatz on the c-numbers
        sol2, _ = express(S2, [A[2], ch.id, A[1] * A[1], B[2]])
        sol3, _ = express(S3, [A[3], B[3], ch.id, A[2], B[2]])
        say("  free exact fit S2 ~ alpha1 A2 + alpha0 1 (+A1^2, B2 (dependent/zero)): coefficients (A2, 1, A1^2, B2) =", None if sol2 is None else [str(x) for x in sol2])
        say("  free exact fit S3 ~ (A3, B3, 1, A2, B2): coefficients =", None if sol3 is None else [str(x) for x in sol3])
    # ---- zeta-independence of the c-numbers and explicit rational functions -------------------------
    say("\n  the c-numbers depend on zeta only through s = sum zeta (inside S1): same dictionary at other twists (u = 3/7, -2/9):")
    for zt in ([Fr(1, 3), Fr(-2, 5), Fr(7, 4)], [Fr(2), Fr(-1, 2), Fr(3, 7)]):
        st = sum(zt)
        Gz = diag3([Ser.exp(c) for c in zt])
        flags = []
        for u in (Fr(3, 7), Fr(-2, 9)):
            nu0 = (u - THETA[0]) * (u - THETA[1])
            nuh = (Ser.const(u) - hS - THETA[0]) * (Ser.const(u) - hS - THETA[1])
            t1 = ch.t_a(1, u, hS, Gz)
            t2n = ch.t_a(2, u, hS, Gz).scale((nuh * nu0).inv())
            A_ = [coeff_op(t1, k).scale(1 / nu0) for k in range(4)]
            B_ = [coeff_op(t2n, k) for k in range(4)]
            _, S2z, S3z = gaudin_cdet(ch, u, diag3(zt))
            S1u = S1_of(u, st)
            S1p, S1pp = S1p_of(u), S1pp_of(u)
            f2 = (S2z + A_[2] - ch.id.scale(S1u ** 2 / 2 + S1p / 2)).is_zero()
            f3 = (S3z - A_[3] + B_[3] - ch.id.scale(S1u ** 3 / 6 + S1u * S1p / 2 + S1pp / 2)).is_zero()
            f1 = is_scalar_op(A_[1], ch.dim) == S1u and is_scalar_op(B_[1], ch.dim) == 2 * S1u
            flags.append((f1, f2, f3))
        say("   zeta = %s (s = %s): [A1=S1,B1=2S1 | S2 dictionary | S3 dictionary] at the two u: %s" % (fs(zt), st, flags))
    import sympy as sp
    uu = sp.symbols("u")
    S1s = sp.Rational(s.numerator, s.denominator) + 3 * sum(1 / (uu - sp.Rational(t.numerator, t.denominator)) for t in THETA)
    c2s = S1s ** 2 / 2 + sp.diff(S1s, uu) / 2
    c3s = S1s ** 3 / 6 + S1s * sp.diff(S1s, uu) / 2 + sp.diff(S1s, uu, 2) / 2
    say("  explicit c-numbers for zeta = %s, theta = %s (s = %s):" % (fs(ZETA), fs(THETA), s))
    say("   S1(u)  = %s" % sp.apart(sp.together(S1s), uu))
    say("   c2(u)  = S1^2/2 + S1'/2        = %s" % sp.apart(sp.together(c2s), uu))
    say("   c3(u)  = S1^3/6 + S1 S1'/2 + S1''/2 = %s" % sp.apart(sp.together(c3s), uu))
    # ---- finite-h cross-check -----------------------------------------------------------------
    say("\n  finite-h cross-check at u = 3/7 (rational h, z = exp(h zeta) truncated at O(h^%d)):" % NSER)
    u = US[0]
    nu0 = (u - THETA[0]) * (u - THETA[1])
    A, B = store[u]["A"], store[u]["B"]
    prev = None
    say("   h        max|(t1/nu-3)/h - (A1 + hA2 + h^2 A3)|    ratio     max|(t2/(nu nu^[-1])-3)/h - (B1+hB2+h^2B3)|   ratio")
    for h in HS:
        zt = [truncated_exp(c, h) for c in ZETA]
        Gt = diag3(zt)
        t1 = ch.t_a(1, u, h, Gt)
        t2 = ch.t_a(2, u, h, Gt)
        nuh = (u - h - THETA[0]) * (u - h - THETA[1])
        lhsA = (t1.scale(1 / nu0) - ch.id.scale(3)).scale(1 / h)
        rhsA = A[1] + A[2].scale(h) + A[3].scale(h * h)
        lhsB = (t2.scale(1 / (nu0 * nuh)) - ch.id.scale(3)).scale(1 / h)
        rhsB = B[1] + B[2].scale(h) + B[3].scale(h * h)
        dA = float(maxabs(lhsA - rhsA))
        dB = float(maxabs(lhsB - rhsB))
        rt = "" if prev is None else ("%.3g   %.3g" % (prev[0] / dA, prev[1] / dB))
        say("   %-8s %.4e                                    %s" % (h, dA, ""))
        say("            B: %.4e  %s" % (dB, rt))
        prev = (dA, dB)
    say("   (expected: residual ~ h^3 A4, i.e. ratio ~ 10^3 per decade)")
    # ---- variants ------------------------------------------------------------------------------
    say("\n  cdet convention variants at u = 3/7 (do S1 scalar, [S2,S3] = 0 and the dictionary survive?)")
    u = US[0]
    A, B = store[u]["A"], store[u]["B"]
    S1u, S1p, S1pp = S1_of(u), S1p_of(u), S1pp_of(u)
    c2 = S1u ** 2 / 2 + S1p / 2
    c3 = S1u ** 3 / 6 + S1u * S1p / 2 + S1pp / 2
    for tr, ro in ((False, False), (False, True), (True, False), (True, True)):
        S1, S2, S3 = gaudin_cdet(ch, u, g, transpose=tr, row_ordered=ro)
        S1b, S2b, S3b = gaudin_cdet(ch, Fr(-2, 9), g, transpose=tr, row_ordered=ro)
        comm = op_commutator(S2, S3b).is_zero() and op_commutator(S3, S3b).is_zero() and op_commutator(S2, S2b).is_zero()
        d2 = (S2 + A[2] - ch.id.scale(c2)).is_zero()
        d3 = (S3 - A[3] + B[3] - ch.id.scale(c3)).is_zero()
        say("   Lax %s, %s-ordered cdet: S1 scalar: %s ; commuting at 2 points: %s ; S2 = -A2 + c2: %s ; S3 = A3 - B3 + c3: %s" %
            ("E_ji" if tr else "E_ij", "row" if ro else "column", is_scalar_op(S1, ch.dim) is not None, comm, d2, d3))
    # large-u asymptotics: S_k -> e_k(zeta)
    big = Fr(10) ** 9
    S1, S2, S3 = gaudin_cdet(ch, big, g)
    e1_, e2_, e3_ = esym(ZETA)
    say("   u -> infinity (u = 1e9): max|S2 - e2(zeta)| = %.3g ; max|S3 - e3(zeta)| = %.3g  (e2 = %s, e3 = %s)" %
        (float(maxabs(S2 - ch.id.scale(e2_))), float(maxabs(S3 - ch.id.scale(e3_))), e2_, e3_))
    return store


# ==============================================================================================
# X3  localised covectors
# ==============================================================================================
def to_ser(x):
    return x if isinstance(x, Ser) else Ser.const(x)


def vcoef(v, k):
    return [to_ser(x).c[k] for x in v]


def vval(v):
    vals = [to_ser(x).val() for x in v]
    vals = [x for x in vals if x is not None]
    return min(vals) if vals else None


def T_fused(ch, lam, v, h, G):
    """one-site fused transfer matrices: T_(1,1)(v) = t2(v); T_(2,1)(v) = t2(v) t1(v+h) - t3(v+h)."""
    if lam == "11":
        return ch.t_a(2, v, h, G)
    return ch.t_a(2, v, h, G) * ch.t_a(1, v + h, h, G) - ch.t_a(3, v + h, h, G)


def op_apply_chain(v, ops):
    for o in ops:
        v = o.vec_times(v)
    return v


def gaudin_onesite_laurent(ch1, zeta):
    """Laurent coefficients at theta of the one-site Gaudin S2 (kmax 2) and S3 (kmax 3), g = diag(zeta)."""
    th = ch1.thetas[0]
    g = diag3(zeta)
    us = [Fr(3, 7), Fr(-2, 9), Fr(5, 3), Fr(-7, 4), Fr(1, 2), Fr(11, 5), Fr(-9, 8)]
    S2s, S3s = [], []
    for u in us:
        _, S2, S3 = gaudin_cdet(ch1, u, g)
        S2s.append(S2)
        S3s.append(S3)
    L2, ok2 = laurent_fit(S2s, us, 2, [th])
    L3, ok3 = laurent_fit(S3s, us, 3, [th])
    return L2, ok2, L3, ok3


def decompose_report(name, target, Om, ops, opnames, ch1, label):
    """target (vector) modulo Omega: least-squares onto {Omega X : X in ops} (X applied from the right)."""
    basis = [op.vec_times(Om) if False else op.vec_times(Om) for op in ops]
    coef, res, used = ls_fit(basis, target)
    say("    %s = sum c_k * Omega.X_k   (exact LS; independent X_k used: %s)" % (name, [opnames[i] for i in used]))
    for i in used:
        say("        c[%s] = %s" % (opnames[i], coef[i]))
    r2 = vec_norm2(res)
    say("        residual ||res||^2 = %s   (||target||^2 = %s) %s" % (float(r2), float(vec_norm2(target)), "-> EXACT" if r2 == 0 else ""))
    return coef, res


def run_X3(ch1a, ch1b, ch2):
    say("=" * 100)
    say("X3  localised covectors (one-site chain at theta_1 / theta_2, exact h-series)")
    say("=" * 100)
    zeta, omega = ZETA, OMEGA
    th1, th2 = THETA
    th12 = th1 - th2
    hS = Ser.h()
    G = diag3([Ser.exp(c) for c in zeta])
    rs = gauss_r_series(zeta, omega)
    OmS = [to_ser(x) for x in covec_Omega(ch1a, rs)]  # Omega_zw as series
    Om0 = covec_Omega(ch1a, gauss_r(zeta, omega))      # Omega_{zeta omega}
    say("\n  Omega_{zeta omega} = e1^T rho(n(zeta,omega)) =", svec(Om0))
    say("  check the series Omega_zw(h=0) equals it:", [x.c[0] for x in OmS] == Om0)
    g1 = gprime(zeta, omega)
    gd = diag3(zeta)
    # operators used in the decomposition (one-site, theta_1 chain)
    def build_ops(ch1):
        L2, ok2, L3, ok3 = gaudin_onesite_laurent(ch1, zeta)
        C = sum_ops([ch1.E1[(k + 1, k + 1)].scale(zeta[k]) for k in range(3)])
        Y = sum_ops([ch1.E1[(j + 1, i + 1)].scale(gd[i][j] - g1[i][j]) for i in range(3) for j in range(3)
                     if gd[i][j] - g1[i][j] != 0])
        ops = [ch1.id, C, L3[(0, 1)], L3[(0, 2)], L3[(0, 3)], Y]
        names = ["1", "C_zeta", "S3[-1]", "S3[-2]", "S3[-3]", "Y=sum(g-g')_ij E_ji"]
        ops2 = [L2[(0, 1)], L2[(0, 2)], L3[(None, 0)]]
        names2 = ["S2[-1]", "S2[-2]", "S3[0]"]
        return ops, names, ops2, names2, ok2 and ok3, L2, L3

    ops, names, ops2, names2, okL, L2, L3 = build_ops(ch1a)
    say("  one-site Gaudin Laurent fits consistent with extra sample points:", okL)
    say("  S3^{1-site}(u) - const: coefficient ops are 8x8; scalar parts: S3[0] scalar =", is_scalar_op(L3[(None, 0)], 8),
        "; S2[-2] scalar =", is_scalar_op(L2[(0, 2)], 8), "; S3[-3] scalar =", is_scalar_op(L3[(0, 3)], 8))
    say("  g' =", fs(g1), " ; g - g' =", fs([[gd[i][j] - g1[i][j] for j in range(3)] for i in range(3)]))
    results = {}
    for lam, absl in (("11", 2), ("21", 3)):
        say("\n  ---- T_(%s,%s)(theta_1 - 2h) on the one-site chain at theta_1 ----" % (lam[0], lam[1]))
        v = Ser.const(th1) - hS * 2
        Tl = T_fused(ch1a, lam, v, hS, G)
        vl = op_apply_chain(OmS, [Tl])
        val = vval(vl)
        say("  valuation of <Omega_zw|T_lambda(theta_1-2h): h^%d (expected h^%d = h^|lambda|)" % (val, absl))
        vs = [vcoef(vl, absl + k) for k in range(4)]
        kappa = vs[0][0] / Om0[0]
        say("  v0 = %s" % svec(vs[0]))
        say("  v0 == kappa * Omega_{zeta omega}: %s  (kappa = %s)" % (all(a == kappa * b for a, b in zip(vs[0], Om0)), kappa))
        say("  v1 = %s" % svec(vs[1]))
        say("  v2 = %s" % svec(vs[2]))
        say("  v3 = %s" % svec(vs[3]))
        # v1 modulo Omega
        t1 = vec_sub(vs[1], vec_scale(Om0, vs[1][0] / Om0[0]))
        say("  v1 mod Omega (component 0 removed): %s" % svec(t1))
        say("  decomposition of v1 in {Omega.1, Omega.C_zeta, Omega.S3[-1], Omega.S3[-2], Omega.S3[-3], Omega.Y}:")
        c1, r1 = decompose_report("v1", vs[1], Om0, ops, names, ch1a, lam)
        say("  same with the larger basis (+ S2[-1], S2[-2], S3[0] Laurent coefficients):")
        c1b, r1b = decompose_report("v1", vs[1], Om0, ops + ops2, names + names2, ch1a, lam)
        # drop 1 and C: only S3 Laurent coefficients and Y
        say("  v1 using ONLY {1, Y, S3[-1], S3[-2], S3[-3]} (no C_zeta):")
        idx = [0, 5, 2, 3, 4]
        decompose_report("v1", vs[1], Om0, [ops[i] for i in idx], [names[i] for i in idx], ch1a, lam)
        # v2: degree <= 2 products
        gens = ops[1:]
        gnames = names[1:]
        prod_ops = [ch1a.id] + gens
        prod_names = ["1"] + gnames
        for i in range(len(gens)):
            for j in range(i, len(gens)):
                prod_ops.append(gens[i] * gens[j])
                prod_names.append(gnames[i] + "*" + gnames[j])
                if i != j:
                    prod_ops.append(gens[j] * gens[i])
                    prod_names.append(gnames[j] + "*" + gnames[i])
        say("  decomposition of v2 in Omega.{1, X_a, X_a X_b} (X in C_zeta, S3[-1..-3], Y; %d vectors):" % len(prod_ops))
        basis = [op.vec_times(Om0) for op in prod_ops]
        coef2, res2, used2 = ls_fit(basis, vs[2])
        say("      rank of the span of those vectors: %d (of 8); residual ||res||^2 = %s %s" %
            (len(used2), float(vec_norm2(res2)), "-> v2 lies in the span" if vec_norm2(res2) == 0 else "-> v2 NOT in span"))
        # also the pure linear span for v2
        c2l, r2l, u2l = ls_fit([op.vec_times(Om0) for op in ops], vs[2])
        say("      linear-only basis: rank %d, residual ||res||^2 = %s" % (len(u2l), float(vec_norm2(r2l))))
        results[lam] = dict(vs=vs, val=val)
        # finite-h numerical check
        say("  finite-h check (rational h, z,w = exp series truncated): max|<Omega_zw|T(theta1-2h)/h^%d - (v0+h v1+h^2 v2+h^3 v3)|" % absl)
        prev = None
        for h in HS:
            zt = [truncated_exp(c, h) for c in zeta]
            wt = [truncated_exp(c, h) for c in omega]
            Omt = covec_Omega(ch1a, gauss_r(zt, wt))
            Tn = T_fused(ch1a, lam, th1 - 2 * h, h, diag3(zt))
            vn = Tn.vec_times(Omt)
            pred = [sum(vs[k][i] * h ** k for k in range(4)) for i in range(8)]
            dev = max(abs(vn[i] / h ** absl - pred[i]) for i in range(8))
            say("     h=%-7s dev = %.4e%s" % (h, float(dev), "" if prev is None else "   ratio %.3g (expect ~10^4 for O(h^4))" % (prev / float(dev))))
            prev = float(dev)
    # ---- spectator factors ----------------------------------------------------------------
    say("\n  spectator factors f_(1,1) = (th12-h)(th12-3h), f_(2,1) = (th12-h)^2(th12-3h), th12 = %s" % th12)
    f11 = [th12 ** 2, -4 * th12, Fr(3)]
    f21 = [th12 ** 3, -5 * th12 ** 2, 7 * th12, Fr(-3)]
    say("    f11 = %s + h*(%s) + h^2*(%s)  -> th12^2 = %s ; relative first-order coeff -4/th12 = %s" % (f11[0], f11[1], f11[2], th12 ** 2, f11[1] / f11[0]))
    say("    f21 = %s + h*(%s) + h^2*(%s) + h^3*(%s) -> th12^3 = %s ; relative first-order coeff -5/th12 = %s" % (f21[0], f21[1], f21[2], f21[3], th12 ** 3, f21[1] / f21[0]))
    # ---- site 2 -----------------------------------------------------------------------------
    say("\n  ---- site 2: one-site chain at theta_2: <Omega|t2(theta_2+h) and <Omega|t2(theta_2+h) t1(theta_2-h) ----")
    ops_b, names_b, ops2_b, names2_b, okLb, L2b, L3b = build_ops(ch1b)
    say("  one-site Gaudin Laurent fits at theta_2 consistent:", okLb)
    v1p = Ser.const(th2) + hS
    v1m = Ser.const(th2) - hS
    wA = op_apply_chain(OmS, [ch1b.t_a(2, v1p, hS, G)])
    wB = op_apply_chain(OmS, [ch1b.t_a(2, v1p, hS, G), ch1b.t_a(1, v1m, hS, G)])
    # effective site-2 covector of the second identity: -f21' w21 + 6 h^4 chi3 f11 Omega
    chi3 = Ser.exp(sum(zeta))
    f11S = Ser(f11)
    f21S = Ser(f21)
    wBeff = [(-f21S) * x + Ser([Fr(0)] * 4 + [Fr(6)]) * chi3 * f11S * y for x, y in zip(wB, OmS)]
    for nm, wv in (("<Omega|t2(th2+h)", wA), ("<Omega|t2(th2+h)t1(th2-h)", wB), ("effective (2,1)-site-2 covector -f21 w + 6h^4 chi3 f11 Omega", wBeff)):
        val = vval(wv)
        ws = [vcoef(wv, val + k) for k in range(4)]
        kappa = ws[0][0] / Om0[0] if Om0[0] != 0 else None
        say("\n  %s : valuation h^%d" % (nm, val))
        say("    w0 = %s" % svec(ws[0]))
        say("    w0 == kappa*Omega_{zeta omega}: %s (kappa = %s)" % (all(a == kappa * b for a, b in zip(ws[0], Om0)), kappa))
        say("    decomposition of w0:")
        decompose_report("w0", ws[0], Om0, ops_b, names_b, ch1b, nm)
        say("    w1 = %s" % svec(ws[1]))
        say("    w2 = %s" % svec(ws[2]))
        say("    decomposition of w1:")
        c, r = decompose_report("w1", ws[1], Om0, ops_b, names_b, ch1b, nm)
    # ---- two-site sanity and X1 operator ------------------------------------------------------
    say("\n  ---- two-site sanity at small h (exact rationals, 64x64): <Om (x) Om|T_lambda(th1-2h) = f_lambda <v_lambda| (x) <Om| etc. ----")
    for h in (Fr(1, 1000),):
        zt = [truncated_exp(c, h) for c in zeta]
        wt = [truncated_exp(c, h) for c in omega]
        Gt = diag3(zt)
        Omt = covec_Omega(ch2, gauss_r(zt, wt))
        O2 = kron_vec(Omt, Omt)
        f11n = (th12 - h) * (th12 - 3 * h)
        f21n = (th12 - h) ** 2 * (th12 - 3 * h)
        u1 = th1 - 2 * h
        t0 = time.time()
        lhs11 = op_apply_chain(O2, [T_fused(ch2, "11", u1, h, Gt)])
        lhs21 = op_apply_chain(O2, [T_fused(ch2, "21", u1, h, Gt)])
        v11 = T_fused(ch1a, "11", u1, h, Gt).vec_times(Omt)
        v21 = T_fused(ch1a, "21", u1, h, Gt).vec_times(Omt)
        r11 = vec_sub(lhs11, vec_scale(kron_vec(v11, Omt), f11n))
        r21 = vec_sub(lhs21, vec_scale(kron_vec(v21, Omt), f21n))
        # site 2
        t2a = ch2.t_a(2, th2 + h, h, Gt)
        lhs2a = t2a.vec_times(O2)
        w2a = ch1b.t_a(2, th2 + h, h, Gt).vec_times(Omt)
        r2a = vec_sub(lhs2a, vec_scale(kron_vec(Omt, w2a), f11n))
        lhs2b = ch2.t_a(1, th2 - h, h, Gt).vec_times(lhs2a) if False else None
        # <Om (x) Om| t2(th2+h) t1(th2-h)
        lhs2b = op_apply_chain(O2, [t2a, ch2.t_a(1, th2 - h, h, Gt)])
        w2b = op_apply_chain(Omt, [ch1b.t_a(2, th2 + h, h, Gt), ch1b.t_a(1, th2 - h, h, Gt)])
        chi3n = zt[0] * zt[1] * zt[2]
        rhs2b = kron_vec(Omt, vec_add(vec_scale(w2b, -f21n), vec_scale(Omt, 6 * h ** 4 * chi3n * f11n)))
        r2b = vec_sub(lhs2b, rhs2b)
        mx = lambda r: float(max(abs(x) for x in r))
        say("    h = %s: max abs residuals  site1 T11: %.3g | site1 T21: %.3g | site2 t2: %.3g | site2 t2 t1: %.3g   (sizes: %.3g %.3g %.3g %.3g; %.1f s)" %
            (h, mx(r11), mx(r21), mx(r2a), mx(r2b), mx(lhs11), mx(lhs21), mx(lhs2a), mx(lhs2b), time.time() - t0))
    return results


def run_X3_X1op(ch2):
    th1, th2 = THETA
    th12 = th1 - th2
    zeta = ZETA
    hS = Ser.h()
    G = diag3([Ser.exp(c) for c in zeta])
    say("\n  ---- leading two-site operator X1 = lim (T_(1,1)(th1-2h)/(h^2 th12^2) - c)/h  (64x64) ----")
    t0 = time.time()
    T11 = ch2.t_a(2, Ser.const(th1) - hS * 2, hS, G)
    lowest = [coeff_op(T11, k) for k in range(5)]
    say("  coefficients of h^0, h^1 vanish:", lowest[0].is_zero(), lowest[1].is_zero(), "; h^2 coefficient is scalar:",
        is_scalar_op(lowest[2], ch2.dim), "= 3 th12^2 =", 3 * th12 ** 2, "(%.1f s)" % (time.time() - t0))
    c_lead = is_scalar_op(lowest[2], ch2.dim) / th12 ** 2
    say("  c = %s ;  X1 = (h^3 coefficient)/th12^2  has %d nonzero entries; scalar part (trace/64) = %s" %
        (c_lead, lowest[3].nnz(), sum(lowest[3].flat().get((i, i), 0) for i in range(ch2.dim)) / ch2.dim / th12 ** 2))
    X1 = lowest[3].scale(1 / th12 ** 2)
    # Gaudin Laurent coefficients, two sites
    g = diag3(zeta)
    us = [Fr(3, 7), Fr(-2, 9), Fr(5, 3), Fr(-7, 4), Fr(1, 2), Fr(11, 5), Fr(-9, 8), Fr(13, 3), Fr(-4, 3), Fr(7, 6)]
    S2s, S3s = [], []
    for u in us:
        _, S2, S3 = gaudin_cdet(ch2, u, g)
        S2s.append(S2)
        S3s.append(S3)
    L2, ok2 = laurent_fit(S2s, us, 2, ch2.thetas)
    L3, ok3 = laurent_fit(S3s, us, 3, ch2.thetas)
    say("  two-site Gaudin Laurent fits (S2: poles up to 2, S3: up to 3) consistent at extra points:", ok2, ok3)
    W = sum_ops([(ch2.E[0][(k + 1, k + 1)] + ch2.E[1][(k + 1, k + 1)]).scale(zeta[k]) for k in range(3)])
    H1 = L2[(0, 1)]
    H2 = L2[(1, 1)]
    solHW, _ = express(H1 + H2, [ch2.id, W])
    say("  residue theorem: H1 + H2 in span{1, W}:", None if solHW is None else "yes, H1+H2 = %s*1 + (%s)*W" % (solHW[0], solHW[1]),
        "; W is scalar:", is_scalar_op(W, ch2.dim) is not None)
    lin = [("1", ch2.id), ("W", W), ("H1", H1)] + [("S3[th%d,-%d]" % (s + 1, k), L3[(s, k)]) for s in range(2) for k in (1, 2, 3)]
    sol, piv = express(X1, [o for _, o in lin])
    ranks = len(piv)
    say("  LINEAR span test: X1 in span{1, W, H1, S3 Laurent coefficients at th1, th2} (rank of the basis = %d of %d): %s" %
        (ranks, len(lin), "YES" if sol is not None else "NO"))
    if sol is not None:
        for (nm, _), c in zip(lin, sol):
            if c != 0:
                say("      %s : %s" % (nm, c))
    # quadratic extension
    gens = [o for nm, o in lin[1:]]
    gn = [nm for nm, o in lin[1:]]
    quad = [("1", ch2.id)] + lin[1:]
    for i in range(len(gens)):
        for j in range(i, len(gens)):
            quad.append((gn[i] + "*" + gn[j], gens[i] * gens[j]))
            if i != j:
                quad.append((gn[j] + "*" + gn[i], gens[j] * gens[i]))
    sol2, piv2 = express(X1, [o for _, o in quad])
    say("  QUADRATIC span test: X1 in span{1, gens, products of two gens} (%d operators, rank %d): %s" %
        (len(quad), len(piv2), "YES" if sol2 is not None else "NO"))
    if sol2 is not None:
        nz = [(quad[i][0], c) for i, c in enumerate(sol2) if c != 0]
        say("      (one solution, free columns set to 0) nonzero coefficients:", nz[:30])
    # ---- localisation content of X1: formal-series identity and rank of Omega0(x)Omega0 . X1 ------
    hS2 = hS
    rsS = gauss_r_series(zeta, OMEGA)
    OmS = [to_ser(x) for x in covec_Omega(ch2, rsS)]
    Om0 = covec_Omega(ch2, gauss_r(zeta, OMEGA))
    OO = [x * y for x in OmS for y in OmS]
    lhs_series = T11.vec_times(OO)           # <Omega_zw (x) Omega_zw| T_(1,1)(th1-2h), exact series, 64 components
    ch1a = Chain(LAM, [th1])
    v11 = ch1a.t_a(2, Ser.const(th1) - hS2 * 2, hS2, G).vec_times(OmS)
    f11 = Ser([th12 ** 2, -4 * th12, Fr(3)])
    rhs_series = [f11 * (x * y) for x in v11 for y in OmS]
    ident = all(all(a == b for a, b in zip(to_ser(l).c[:12], to_ser(r).c[:12])) for l, r in zip(lhs_series, rhs_series))
    say("  formal-series identity <Om(x)Om|T11(th1-2h) = f11 <v11|(x)<Om| (all 64 components, h^0..h^11, exact): %s" % ident)
    OO0 = kron_vec(Om0, Om0)
    Y = X1.vec_times(OO0)
    Mrows = [[Y[8 * i + j] for j in range(8)] for i in range(8)]
    rk = matrix_rank_exact(Mrows)
    say("  Omega0 (x) Omega0 . X1 reshaped to an 8x8 matrix (site1 x site2): exact rank = %d  (rank 1 would mean a pure product v (x) w)" % rk)
    if rk == 1:
        uvec = [Mrows[i][0] / Om0[0] for i in range(8)]
        prod_ok = all(Mrows[i][j] == uvec[i] * Om0[j] for i in range(8) for j in range(8))
        say("  => Omega0(x)Omega0 . X1 = u (x) Omega0 exactly (%s) with the one-site covector (site 1) u = %s" % (prod_ok, svec(uvec)))
        L2a, ok2a, L3a, ok3a = gaudin_onesite_laurent(ch1a, zeta)
        Cz = sum_ops([ch1a.E1[(k + 1, k + 1)].scale(zeta[k]) for k in range(3)])
        ops_u = [ch1a.id, Cz, L3a[(0, 1)], L3a[(0, 2)], L3a[(0, 3)]]
        nms_u = ["1", "C_zeta", "S3[-1]", "S3[-2]", "S3[-3]"]
        decompose_report("u", uvec, Om0, ops_u, nms_u, ch1a, "u")
    return X1


# ==============================================================================================
# X4  Baxter -> oper (sympy, scalars only)
# ==============================================================================================
def run_X4():
    import sympy as sp
    say("=" * 100)
    say("X4  h-expansion of the XXX direct / dual Baxter equations (sympy, symbolic h, L = 2 adjoint, nu(x) = (x-a)(x-b))")
    say("=" * 100)
    x, h, s, a, b = sp.symbols("x h s a b")
    ORD = 3
    Qf = sp.Function("Q")(x)
    Af = {k: sp.Function("A%d" % k)(x) for k in (1, 2, 3)}
    Bf = {k: sp.Function("B%d" % k)(x) for k in (1, 2, 3)}
    A1, A2, A3, B1, B2, B3 = Af[1], Af[2], Af[3], Bf[1], Bf[2], Bf[3]

    def taylor(f, k, n=ORD + 1):
        return sum((k * h) ** m / sp.factorial(m) * sp.diff(f, x, m) for m in range(n + 1))

    nu = lambda y: (y - a) * (y - b)
    lam = {1: 2, 2: 1, 3: 0}
    nun = lambda n, y: (y + h * lam[n] - a) * (y + h * lam[n] - b)
    chi3 = sum((h * s) ** m / sp.factorial(m) for m in range(ORD + 2))   # chi3 = exp(h s), s = zeta_1+zeta_2+zeta_3

    def T1(k):
        y = x + k * h
        return nu(y) * (3 + h * taylor(A1, k) + h ** 2 * taylor(A2, k) + h ** 3 * taylor(A3, k))

    def T2(k):
        y = x + k * h
        return nu(y) * nu(y - h) * (3 + h * taylor(B1, k) + h ** 2 * taylor(B2, k) + h ** 3 * taylor(B3, k))

    def T3(k):
        y = x + k * h
        return chi3 * nun(1, y) * nun(2, y - h) * nun(3, y - 2 * h)

    direct = (nun(1, x) * nun(1, x + h) * nun(1, x + 2 * h) * taylor(Qf, -1) - T1(0) * nun(1, x + h) * nun(1, x + 2 * h) * Qf
              + T2(1) * nun(1, x + 2 * h) * taylor(Qf, 1) - T3(2) * taylor(Qf, 2))
    dual = (nun(3, x) * nun(3, x - h) * nun(3, x - 2 * h) * taylor(Qf, 1) - T1(0) * nun(3, x - h) * nun(3, x - 2 * h) * Qf
            + T2(0) * nun(3, x - 2 * h) * taylor(Qf, -1) - T3(0) * taylor(Qf, -2))
    say("  direct: nu1(x)nu1(x+h)nu1(x+2h)Q(x-h) - t1(x)nu1(x+h)nu1(x+2h)Q(x) + t2(x+h)nu1(x+2h)Q(x+h) - t3(x+2h)Q(x+2h) = 0")
    say("  dual  : nu3(x)nu3(x-h)nu3(x-2h)Q(x+h) - t1(x)nu3(x-h)nu3(x-2h)Q(x) + t2(x)nu3(x-2h)Q(x-h) - t3(x)Q(x-2h) = 0")
    say("  with t1(x) = nu(x)[3 + hA1 + h^2 A2 + h^3 A3], t2(x) = nu(x)nu(x-h)[3 + hB1 + h^2 B2 + h^3 B3], chi3 = exp(h s), "
        "nu_n(x) = nu(x + h Lam_n), Lam = (2,1,0); A_k, B_k unknown functions of x (shifted arguments Taylor-expanded).")
    q = sp.symbols("q0:6")

    def orders(expr):
        ex = sp.expand(expr)
        return [ex.coeff(h, k) for k in range(ORD + 1)]

    def decompose(expr, maxd=4):
        ex = sp.expand(expr)
        for n in range(maxd, 0, -1):
            ex = ex.subs(sp.diff(Qf, x, n), q[n])
        ex = sp.expand(ex.subs(Qf, q[0]))
        return [ex.coeff(q[n]) for n in range(maxd + 1)]

    nuf = (x - a) * (x - b)
    dnu = sp.diff(nuf, x)
    S1 = s + 3 * dnu / nuf
    S1p, S1pp = sp.diff(S1, x), sp.diff(S1, x, 2)
    res = {}
    for name, eq in (("direct", direct), ("dual", dual)):
        t0 = time.time()
        o = orders(eq)
        say("\n  ===== %s Baxter equation =====" % name)
        say("  order h^0: %s" % sp.simplify(o[0]))
        d1 = decompose(o[1])
        c1 = sp.factor(d1[0])
        say("  order h^1: only Q appears; coefficient = %s ; derivative-of-Q coefficients: %s" % (c1, [sp.simplify(c) for c in d1[1:]]))
        say("      => vanishes iff  B1 = A1 + s + 3 nu'/nu  (one constraint tying B1 to A1; A1 itself still free)")
        B1sol = A1 + s + 3 * dnu / nuf
        o2 = o[2].subs(B1, B1sol).doit()
        d2 = decompose(o2)
        cq1 = sp.factor(d2[1])
        say("  order h^2 (after B1 = A1 + s + 3 nu'/nu): Q' coefficient = %s" % cq1)
        say("      => vanishes iff  A1 = s + 3 nu'/nu = S1(x) (the Gaudin tr L_g), hence B1 = 2 A1 = 2 S1")
        o2b = o[2].subs({A1: S1, B1: 2 * S1}).doit()
        d2b = decompose(o2b)
        say("      with A1 = S1, B1 = 2 S1: Q' coefficient = %s ; Q'', Q''' coefficients = %s" % (sp.simplify(d2b[1]), [sp.simplify(c) for c in d2b[2:]]))
        B2sol = sp.solve(d2b[0], B2)[0]
        B2dict = A2 - S1p + (s ** 2 / 2 + 3 * s * dnu / nuf + 9 / nuf)
        say("      Q coefficient vanishes iff B2 = A2 - S1' + (s^2/2 + 3 s nu'/nu + 9/nu)  [= A2 + S1^2/2 + S1'/2]: check difference = %s" %
            sp.simplify(B2sol - B2dict))
        say("      (so order h^2 imposes: B1 = 2A1, A1 = S1, B2 = A2 + S1^2/2 + S1'/2; no constraint on A2, A3, B3)")
        o3 = o[3].subs({A1: S1, B1: 2 * S1}).doit().subs(B2, B2sol).doit()
        d3 = decompose(o3)
        lead = sp.factor(d3[3])
        say("  order h^3 ODE: [%s] * ( Q''' + p2 Q'' + p1 Q' + p0 Q ) = 0" % lead)
        mon = [sp.cancel(sp.together(d3[n] / d3[3])) for n in range(3)]
        res[name] = mon
        say("      p2 = %s" % sp.factor(mon[2]))
        say("      p1 = %s" % sp.factor(mon[1]))
        say("      p0 = (long; verified below)  [%.0f s]" % (time.time() - t0))
    # ---- compare with the oper -----------------------------------------------------------------
    c2 = S1 ** 2 / 2 + S1p / 2
    c3 = S1 ** 3 / 6 + S1 * S1p / 2 + S1pp / 2
    S2d = -A2 + c2          # X2 dictionary
    S3d = A3 - B3 + c3
    say("\n  X2 dictionary used: S1 = A1 = s + 3 nu'/nu ; S2 = -A2 + S1^2/2 + S1'/2 ; S3 = A3 - B3 + S1^3/6 + S1 S1'/2 + S1''/2")
    mon = res["dual"]
    chk = [sp.simplify(mon[2] + S1), sp.simplify(mon[1] - S2d), sp.simplify(mon[0] + S3d)]
    say("  DUAL: monic ODE minus (Q''' - S1 Q'' + S2 Q' - S3 Q) coefficients (Q'', Q', Q):", chk)
    say("        => dual Baxter at order h^3 is EXACTLY (d^3 - S1 d^2 + S2 d - S3) Q = 0, the Talalaev oper cdet(d - L_g), acting on Q^j (no gauge factor)")
    mon = res["direct"]
    # Q = nu^p R : R''' + (p2 + 3 l1) R'' + (p1 + 2 p2 l1 + 3 l2) R' + (p0 + p1 l1 + p2 l2 + l3) R = 0
    adj = [S1pp + sp.diff(S2d, x) + S3d, 2 * S1p + S2d, S1]
    for pw in (0, 1, 2):
        phi = nuf ** pw
        l1, l2, l3 = [sp.diff(phi, x, k) / phi for k in (1, 2, 3)]
        cn = [mon[0] + mon[1] * l1 + mon[2] * l2 + l3, mon[1] + 2 * mon[2] * l1 + 3 * l2, mon[2] + 3 * l1]
        dif = [sp.simplify(cn[n] - adj[n]) for n in range(3)]
        zero = all(d == 0 for d in dif)
        say("  DIRECT with Q = nu(x)^%d R: (R coefficients) - (formal-adjoint-oper coefficients [(S1''+S2'+S3), (2S1'+S2), S1]) = zero: %s" % (pw, zero))
    say("        => direct Baxter at order h^3 is the formal adjoint (-d)^3 - (-d)^2 S1 + (-d) S2 - S3 acting on R = Q/nu^2, i.e. on "
        "Q = nu(x)^2 R with nu = prod_s (x - theta_s).")
    return res


# ==============================================================================================
def main():
    t0 = time.time()
    ch2 = Chain(LAM, THETA)
    ch1a = Chain(LAM, [THETA[0]])
    ch1b = Chain(LAM, [THETA[1]])
    run_checks(ch2)
    run_X1(ch1a)
    run_gaudin_eigencovector(ch2)
    run_X2(ch2)
    run_X3(ch1a, ch1b, ch2)
    run_X3_X1op(ch2)
    run_X4()
    say("\nTOTAL RUNTIME %.0f s" % (time.time() - t0))


if __name__ == "__main__":
    main()
