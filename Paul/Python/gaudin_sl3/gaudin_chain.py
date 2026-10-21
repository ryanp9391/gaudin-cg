"""sl(3)/gl(3) Gaudin chain, L = 2 sites, exact sparse Fraction algebra + mpmath eigenstates.

Independent Python implementation (agent P1).  Conventions (see RESULTS-Qside.md section 0):
  L_g(u)_ij = g_ij + sum_s E^(s)_ij/(u - theta_s),   E^(1)=E x 1, E^(2)=1 x E (Kronecker order)
  cdet(d - L_g) = sum_sigma sgn(sigma) M_{s(1)1} M_{s(2)2} M_{s(3)3},  M = d*1 - L_g, d f = f d + f'
                = d^3 - S1 d^2 + S2 d - S3
Matrices act on kets as columns: E e_col = sum_row M[row,col] e_row; covectors act from the right.
"""
import sys, os, itertools, pickle
from fractions import Fraction as Fr
from math import comb, factorial
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import sympy
import mpmath as mp
from gl3_reps import GL3Rep

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------- sparse matrices
class SM:
    """sparse matrix: rows dict r -> {c: Fraction}"""
    __slots__ = ('n', 'r')
    def __init__(self, n, r=None):
        self.n = n
        self.r = r if r is not None else {}
    @staticmethod
    def eye(n, c=1):
        return SM(n, {i: {i: Fr(c)} for i in range(n)})
    def copy(self):
        return SM(self.n, {i: dict(row) for i, row in self.r.items()})
    def __mul__(self, o):
        if isinstance(o, SM):
            out = {}
            for i, row in self.r.items():
                acc = {}
                for k, a in row.items():
                    orow = o.r.get(k)
                    if orow:
                        for c, b in orow.items():
                            acc[c] = acc.get(c, 0) + a * b
                acc = {c: v for c, v in acc.items() if v != 0}
                if acc:
                    out[i] = acc
            return SM(self.n, out)
        o = Fr(o)
        if o == 0:
            return SM(self.n)
        return SM(self.n, {i: {c: v * o for c, v in row.items()} for i, row in self.r.items()})
    __rmul__ = lambda self, o: self.__mul__(o)
    def __add__(self, o):
        out = {i: dict(row) for i, row in self.r.items()}
        for i, row in o.r.items():
            d = out.setdefault(i, {})
            for c, v in row.items():
                nv = d.get(c, 0) + v
                if nv == 0:
                    d.pop(c, None)
                else:
                    d[c] = nv
            if not d:
                out.pop(i)
        return SM(self.n, out)
    def __neg__(self):
        return self * (-1)
    def __sub__(self, o):
        return self + (-o)
    def to_mp(self, idx=None):
        idx = list(range(self.n)) if idx is None else idx
        pos = {j: k for k, j in enumerate(idx)}
        A = mp.zeros(len(idx), len(idx))
        for i in idx:
            for c, v in self.r.get(i, {}).items():
                if c in pos:
                    A[pos[i], pos[c]] = mp.mpf(v.numerator) / v.denominator
        return A
    def apply(self, vec):
        """M * ket (list of numbers)"""
        out = [0] * self.n
        for i, row in self.r.items():
            s = 0
            for c, v in row.items():
                s += v * vec[c]
            out[i] = s
        return out
    def rapply(self, cov):
        """covector * M"""
        out = [0] * self.n
        for i, row in self.r.items():
            ci = cov[i]
            if ci == 0:
                continue
            for c, v in row.items():
                out[c] += ci * v
        return out
    def is_zero(self):
        return not self.r

def sm_from_sympy(M):
    n = M.shape[0]
    r = {}
    for i in range(n):
        for j in range(n):
            v = M[i, j]
            if v != 0:
                q = sympy.Rational(v)
                r.setdefault(i, {})[j] = Fr(int(q.p), int(q.q))
    return SM(n, r)

def kron(A, B):
    n = A.n * B.n
    out = {}
    for i, ra in A.r.items():
        for j, rb in B.r.items():
            row = {}
            for c, a in ra.items():
                for d, b in rb.items():
                    row[c * B.n + d] = a * b
            out[i * B.n + j] = row
    return SM(n, out)

def comm(A, B):
    return A * B - B * A

# ---------------------------------------------------------------- the rep
LAM = (2, 1, 0)
THETA = (Fr(2, 7), Fr(-3, 5))

class Chain:
    def __init__(self, lam=LAM, theta=THETA):
        self.lam = tuple(lam)
        self.theta = tuple(theta)
        rep = GL3Rep(self.lam)
        self.rep = rep
        d = rep.dim
        self.d1 = d
        self.dim = d * d
        self.E1 = {(i, j): sm_from_sympy(rep.generator(i, j)) for (i, j) in
                   [(1, 1), (2, 2), (3, 3), (1, 2), (2, 3), (2, 1), (3, 2)]}
        self.E1[(1, 3)] = comm(self.E1[(1, 2)], self.E1[(2, 3)])   # [E12,E23] = E13
        self.E1[(3, 1)] = comm(self.E1[(3, 2)], self.E1[(2, 1)])   # [E32,E21] = E31
        # sanity: [E13,E31] = E11 - E33
        assert (comm(self.E1[(1, 3)], self.E1[(3, 1)]) - (self.E1[(1, 1)] - self.E1[(3, 3)])).is_zero()
        I1 = SM.eye(d)
        self.Es = {}   # (s, i, j) -> 64x64
        for i in (1, 2, 3):
            for j in (1, 2, 3):
                self.Es[(1, i, j)] = kron(self.E1[(i, j)], I1)
                self.Es[(2, i, j)] = kron(I1, self.E1[(i, j)])
        hw = tuple(tuple(self.lam[:k]) for k in range(1, 4))
        self.hw_index = rep.patterns.index(hw)
        # weights of two-site basis
        w1 = [tuple(int(self.E1[(k, k)].r.get(a, {}).get(a, 0)) for k in (1, 2, 3)) for a in range(d)]
        self.w1 = w1
        self.wt = [tuple(w1[a][k] + w1[b][k] for k in range(3)) for a in range(d) for b in range(d)]
        self.sectors = {}
        for idx, w in enumerate(self.wt):
            self.sectors.setdefault(w, []).append(idx)

    def check_hw(self):
        e = self.hw_index
        ok = {}
        for (i, j) in [(1, 2), (1, 3), (2, 3)]:
            col = [self.E1[(i, j)].r.get(r, {}).get(e, 0) for r in range(self.d1)]
            ok[(i, j)] = all(c == 0 for c in col)
        return self.hw_index, ok

# ---------------------------------------------------------------- symbolic normal ordering
# algebra element: dict {(m, mono): coeff}; mono = tuple of (i, j, deriv-order); m = power of d placed right
def _dmono(mono):
    out = []
    for p in range(len(mono)):
        i, j, k = mono[p]
        out.append(mono[:p] + ((i, j, k + 1),) + mono[p + 1:])
    return out

def _dk(mono_dict, k):
    """k-th derivative of a dict {mono: coeff}"""
    cur = mono_dict
    for _ in range(k):
        nxt = {}
        for mono, c in cur.items():
            for m2 in _dmono(mono):
                nxt[m2] = nxt.get(m2, 0) + c
        cur = nxt
    return cur

def alg_mul(A, B):
    out = {}
    for (m, ma), ca in A.items():
        for (n, mb), cb in B.items():
            for k in range(m + 1):
                for mono, c in _dk({mb: 1}, k).items():
                    key = (m - k + n, ma + mono)
                    out[key] = out.get(key, 0) + ca * cb * comb(m, k) * c
    return {k: v for k, v in out.items() if v != 0}

def sgn(p):
    s = 1
    p = list(p)
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s

_SYM_CACHE = {}
def symbolic_cdet(rowordered=False):
    """returns (S1, S2, S3) as dicts {mono: coeff}, with cdet = d^3 - S1 d^2 + S2 d - S3.
    rowordered: product over columns of M_{1 s(1)} M_{2 s(2)} M_{3 s(3)} (the variant)."""
    key = rowordered
    if key in _SYM_CACHE:
        return _SYM_CACHE[key]
    def M(i, j):
        e = {}
        if i == j:
            e[(1, ())] = 1
        e[(0, ((i, j, 0),))] = -1
        return e
    tot = {}
    for p in itertools.permutations((1, 2, 3)):
        el = {(0, ()): 1}
        for k in (1, 2, 3):
            el = alg_mul(el, M(p[k - 1], k) if not rowordered else M(k, p[k - 1]))
        for kk, v in el.items():
            tot[kk] = tot.get(kk, 0) + sgn(p) * v
    tot = {k: v for k, v in tot.items() if v != 0}
    assert all(m <= 3 for (m, _) in tot)
    S1 = {mono: -c for (m, mono), c in tot.items() if m == 2}
    S2 = {mono: c for (m, mono), c in tot.items() if m == 1}
    S3 = {mono: -c for (m, mono), c in tot.items() if m == 0}
    assert tot.get((3, ()), 0) == 1
    _SYM_CACHE[key] = (S1, S2, S3)
    return _SYM_CACHE[key]

def sym_str(S):
    def f(s):
        i, j, k = s
        return 'l%s_%d%d' % ("'" * k, i, j)
    terms = []
    for mono, c in sorted(S.items(), key=lambda t: (len(t[0]), t[0])):
        terms.append('%+d %s' % (c, ' '.join(f(s) for s in mono) if mono else '1'))
    return ' '.join(terms)

# ---------------------------------------------------------------- numeric instantiation
class Operators:
    """S_k(u) for twist matrix g (3x3 Fractions), exact; plus partial-fraction coefficient operators."""
    def __init__(self, chain, g, transpose=False, rowordered=False):
        self.ch = chain
        self.g = [[Fr(x) for x in row] for row in g]
        self.transpose = transpose
        self.rowordered = rowordered

    def Lmat(self, u, i, j, deriv=0):
        ch = self.ch
        ie, je = (j, i) if self.transpose else (i, j)
        out = SM(ch.dim)
        for s in (1, 2):
            c = Fr((-1) ** deriv * factorial(deriv)) / (Fr(u) - ch.theta[s - 1]) ** (deriv + 1)
            out = out + ch.Es[(s, ie, je)] * c
        if deriv == 0 and self.g[i - 1][j - 1] != 0:
            out = out + SM.eye(ch.dim, self.g[i - 1][j - 1])
        return out

    def S_at(self, u):
        """(S1, S2, S3) at rational u as SM"""
        ch = self.ch
        Sy = symbolic_cdet(self.rowordered)
        cache = {}
        def L(i, j, d):
            k = (i, j, d)
            if k not in cache:
                cache[k] = self.Lmat(u, i, j, d)
            return cache[k]
        pcache = {(): None}
        def prod(mono):
            if mono in pcache:
                return pcache[mono]
            head = prod(mono[:-1])
            last = L(*mono[-1])
            r = last if head is None else head * last
            pcache[mono] = r
            return r
        res = []
        for S in Sy:
            tot = SM(ch.dim)
            for mono, c in S.items():
                if mono == ():
                    tot = tot + SM.eye(ch.dim, c)
                else:
                    tot = tot + prod(mono) * c
            res.append(tot)
        return res

    # partial fraction bases
    def basis(self, k, u):
        a = [1 / (Fr(u) - t) for t in self.ch.theta]
        if k == 2:
            return [Fr(1), a[0], a[0] ** 2, a[1], a[1] ** 2]
        if k == 3:
            return [Fr(1), a[0], a[0] ** 2, a[0] ** 3, a[1], a[1] ** 2, a[1] ** 3]
        raise ValueError

    def fit(self, upts=None):
        """exact coefficient operators C2[b], C3[b] with S_k(u) = sum_b basis_b(u) C_k[b]; returns also check residual"""
        pts = upts or [Fr(3, 7), Fr(-2, 9), Fr(5, 3), Fr(7, 2), Fr(-4, 5), Fr(11, 4), Fr(-7, 3), Fr(13, 5), Fr(1, 3)]
        data = [self.S_at(u) for u in pts]
        out = {}
        for k, nb in ((2, 5), (3, 7)):
            Phi = sympy.Matrix([[sympy.Rational(x.numerator, x.denominator) for x in self.basis(k, u)] for u in pts[:nb]])
            Pinv = Phi.inv()
            C = []
            for b in range(nb):
                op = SM(self.ch.dim)
                for kk in range(nb):
                    w = Pinv[b, kk]
                    if w != 0:
                        op = op + data[kk][k - 1] * Fr(int(w.p), int(w.q))
                C.append(op)
            # check at remaining points
            maxres = 0
            for kk in range(nb, len(pts)):
                pred = SM(self.ch.dim)
                for b in range(nb):
                    pred = pred + C[b] * self.basis(k, pts[kk])[b]
                diff = pred - data[kk][k - 1]
                if not diff.is_zero():
                    maxres += 1
            out[k] = C
            out['res%d' % k] = maxres
        # S1: scalar
        S1 = data[0][0]
        out['S1'] = S1
        out['data0'] = data[0]
        return out

    def eval_S(self, C, k, u):
        b = self.basis(k, u)
        op = SM(self.ch.dim)
        for bb in range(len(b)):
            op = op + C[k][bb] * b[bb]
        return op

def S1_scalar(chain, zeta, u):
    lam1 = sum(chain.lam)
    return sum(Fr(z) for z in zeta) + sum(Fr(lam1) / (Fr(u) - t) for t in chain.theta)

# ---------------------------------------------------------------- eigenstates
def eigenstates(chain, ops, C, u0=Fr(3, 7), u1=Fr(-2, 9), pi=mp.mpf(1) / 3 * mp.pi, dps=60):
    """per weight sector, right eigenvectors of S2(u0) + pi*S3(u1). Returns list of dicts:
    {sector, idx(list of basis indices), ev, vec (full-length mp list), lam2[b], lam3[b], resid}"""
    mp.mp.dps = dps
    A2 = ops.eval_S(C, 2, u0)
    A3 = ops.eval_S(C, 3, u1)
    states = []
    for w, idx in sorted(chain.sectors.items()):
        # restrict the combined operator
        M = A2.to_mp(idx) + pi * A3.to_mp(idx)
        if len(idx) == 1:
            E = [M[0, 0]]
            ER = mp.matrix([[1]])
        else:
            E, ER = mp.eig(M)
        for k in range(len(idx)):
            v = [mp.mpf(0)] * chain.dim
            for kk, j in enumerate(idx):
                v[j] = ER[kk, k]
            # normalise: largest component = 1
            m = max(range(chain.dim), key=lambda j: abs(v[j]))
            v = [x / v[m] for x in v]
            st = dict(sector=w, idx=idx, ev=E[k], vec=v, pivot=m)
            states.append(st)
    return states

def op_on_vec(op, v):
    return op.apply(v)

def eigenvalue_of(op, st):
    v = st['vec']
    w = op.apply(v)
    m = st['pivot']
    lam = w[m] / v[m]
    res = max(abs(w[j] - lam * v[j]) for j in range(len(v)))
    return lam, res

def eigenfunctions(chain, C, states):
    """attach lam2 (5), lam3 (7) coefficient eigenvalues to each state; returns max residual"""
    mx = 0
    for st in states:
        st['lam2'] = []
        st['lam3'] = []
        for k, key in ((2, 'lam2'), (3, 'lam3')):
            for op in C[k]:
                lam, res = eigenvalue_of(op, st)
                st[key].append(lam)
                mx = max(mx, res)
    return mx

def S2_eval(st, chain, u):
    a = [1 / (mp.mpf(u.numerator) / u.denominator - mp.mpf(t.numerator) / t.denominator) for t in chain.theta]
    b = [1, a[0], a[0] ** 2, a[1], a[1] ** 2]
    return sum(l * x for l, x in zip(st['lam2'], b))

# ---------------------------------------------------------------- vacuum
def gauss_r(zeta, om):
    z1, z2, z3 = [Fr(x) for x in zeta]
    w1, w2, w3 = [Fr(x) for x in om]
    r1 = (w1 - z1) / (w1 - z2)
    r2 = (w1 - z1) * (z1 - w2) * (z2 - z3) / ((z1 - z2) * (w1 - z3) * (w2 - z3))
    r3 = (w1 - z2) * (z2 - w2) * (z1 - z3) / ((z1 - z2) * (w1 - z3) * (z3 - w2))
    return r1, r2, r3

def exp_nil(E, r, n=6):
    out = SM.eye(E.n)
    term = SM.eye(E.n)
    for k in range(1, n + 1):
        term = term * E * (Fr(r) / k)
        if term.is_zero():
            break
        out = out + term
    return out

def vacuum_covector_1site(chain, zeta, om):
    r1, r2, r3 = gauss_r(zeta, om)
    e = [Fr(0)] * chain.d1
    e[chain.hw_index] = Fr(1)
    v = exp_nil(chain.E1[(1, 2)], r1).rapply(e)
    v = exp_nil(chain.E1[(1, 3)], r2).rapply(v)
    v = exp_nil(chain.E1[(2, 3)], r3).rapply(v)
    return v

def vacuum_covector(chain, zeta, om):
    v = vacuum_covector_1site(chain, zeta, om)
    return [a * b for a in v for b in v]

def Pmat(zeta):
    z = [Fr(x) for x in zeta]
    e1 = sum(z)
    e2 = z[0] * z[1] + z[0] * z[2] + z[1] * z[2]
    return sympy.Matrix([[1, zi - e1, zi ** 2 - e1 * zi + e2] for zi in
                         [sympy.Rational(x.numerator, x.denominator) for x in z]]).applyfunc(
        lambda x: sympy.nsimplify(x)) if False else sympy.Matrix(
        [[1, sympy.Rational(zi.numerator, zi.denominator) - sympy.Rational(e1.numerator, e1.denominator),
          sympy.Rational((zi ** 2 - e1 * zi + e2).numerator, (zi ** 2 - e1 * zi + e2).denominator)] for zi in z])

def gc(om):
    w = [Fr(x) for x in om]
    e1 = sum(w)
    e2 = w[0] * w[1] + w[0] * w[2] + w[1] * w[2]
    e3 = w[0] * w[1] * w[2]
    R = lambda x: sympy.Rational(x.numerator, x.denominator)
    return sympy.Matrix([[R(e1), -R(e2), R(e3)], [1, 0, 0], [0, 1, 0]])

def gprime(zeta, om, variant=0):
    P = Pmat(zeta)
    G = gc(om)
    if variant == 0:
        M = P * G * P.inv()
    elif variant == 1:
        M = P.inv() * G * P
    elif variant == 2:
        M = (P * G * P.inv()).T
    elif variant == 3:
        M = (P.inv() * G * P).T
    elif variant == 4:
        M = P.T * G * P.T.inv()
    elif variant == 5:
        M = P.T.inv() * G * P.T
    else:
        raise ValueError
    return [[Fr(int(M[i, j].p), int(M[i, j].q)) for j in range(3)] for i in range(3)]

def covector_eigen_test(ops_C, cov, chain):
    """is cov a left eigencovector of the coefficient operators? return (ok, {name: eigenvalue}) exact"""
    res = {}
    ok = True
    n0 = next(j for j in range(len(cov)) if cov[j] != 0)
    for k in (2, 3):
        lams = []
        for op in ops_C[k]:
            w = op.rapply(cov)
            lam = w[n0] / cov[n0]
            if any(w[j] != lam * cov[j] for j in range(len(cov))):
                ok = False
            lams.append(lam)
        res[k] = lams
    return ok, res

# ---------------------------------------------------------------- caching helper
def cached(name, fn):
    path = os.path.join(HERE, name)
    if os.path.exists(path):
        with open(path, 'rb') as f:
            return pickle.load(f)
    val = fn()
    with open(path, 'wb') as f:
        pickle.dump(val, f)
    return val
