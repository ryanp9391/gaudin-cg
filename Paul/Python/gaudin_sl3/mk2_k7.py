"""K7: span tests in X,Y (symbol) coordinates. Q^{(lam)}_n as vectors in Sym^lam(X) (x) Sym^lam(Y), homogenised (X0,X1,X2)=(1,E1,E2)."""
import sys, itertools, pickle
import sympy as sp
from sympy import Rational as R
import master_kernel2 as m
X = sp.symbols('X0 X1 X2'); Yv = sp.symbols('Y0 Y1 Y2')

def monos_of(lam):
    ms = [e for e in itertools.product(range(lam + 1), repeat=3) if sum(e) == lam]
    return ms

def to_vec(poly_expr, lam, idx):
    """poly in E1,E2,H1,H2 (numeric a,ap already substituted) -> homogenise degree lam in X (E1->X1,E2->X2, X0 pad)"""
    p = sp.Poly(sp.expand(poly_expr), m.E1, m.E2, m.H1, m.H2)
    v = {}
    for (i, j, k, l), c in p.terms():
        assert i + j <= lam and k + l <= lam, (i, j, k, l)
        key = ((lam - i - j, i, j), (lam - k - l, k, l))
        v[key] = v.get(key, 0) + c
    return v

def build(lam, a, ap):
    Q = m.symbol_orders_exact(lam)
    out = {}
    for n, q in Q.items():
        out[n] = to_vec(q.subs({m.a_: a, m.ap_: ap}), lam, None)
    return out

def space(lam):
    ms = monos_of(lam)
    idx = {(u, v): k for k, (u, v) in enumerate(itertools.product(ms, ms))}
    return ms, idx

def op_matrix(lam, site, i, j):
    """X_i d_{X_j} on site (1 or 2) as sparse dict-of-dicts {col: {row: coeff}} over basis (u,v)"""
    ms, idx = space(lam)
    mat = sp.zeros(len(idx), len(idx))
    for (u, v), c in idx.items():
        w = u if site == 1 else v
        if w[j] == 0: continue
        w2 = list(w); coef = w2[j]; w2[j] -= 1; w2[i] += 1; w2 = tuple(w2)
        key = (w2, v) if site == 1 else (u, w2)
        mat[idx[key], c] += coef
    return mat

def casimir(lam, ops=None):
    ms, idx = space(lam)
    n = len(idx)
    E1 = {(i, j): op_matrix(lam, 1, i, j) for i in range(3) for j in range(3)}
    E2 = {(i, j): op_matrix(lam, 2, i, j) for i in range(3) for j in range(3)}
    Om = sp.zeros(n, n)
    for i in range(3):
        for j in range(3):
            Om += E1[(i, j)] * E2[(j, i)]
    return E1, E2, Om

def irrep_omega(lam):
    c = lambda mu: sum(mu[i] * (mu[i] + 4 - 2 * (i + 1)) for i in range(3))
    out = {}
    for k in range(lam + 1):
        mu = (2 * lam - k, k, 0)
        out[mu] = R(c(mu) - 2 * c((lam, 0, 0)), 2)
    return out

def projectors(lam, Om):
    ws = irrep_omega(lam)
    n = Om.shape[0]
    I = sp.eye(n)
    P = {}
    for mu, w in ws.items():
        M = I
        for nu, w2 in ws.items():
            if nu != mu:
                M = M * (Om - w2 * I) / (w - w2)
        P[mu] = M
    return P

def vec(d, lam):
    ms, idx = space(lam)
    v = sp.zeros(len(idx), 1)
    for k, c in d.items():
        v[idx[k]] += c
    return v

if __name__ == '__main__':
    lam = int(sys.argv[1]); a = R(1, 3); ap = R(1, 5)
    Qn = build(lam, a, ap)
    E1, E2, Om = casimir(lam)
    P = projectors(lam, Om)
    print('P sums to I:', sum(P.values(), sp.zeros(*Om.shape)) == sp.eye(Om.shape[0]), 'dims', {mu: p.rank() for mu, p in P.items()})
    vs = {n: vec(d, lam) for n, d in Qn.items()}
    # check Q_1 = z (lam^2 - Om) v ; with t-normalisation: Q_1 = lam^2 Pi^(lam-1) K
    v0 = vs[0]
    chk = (lam ** 2 * sp.eye(Om.shape[0]) - Om) * v0 * (a * ap / (ap - a))
    print('Q1 == aa\'/(a\'-a) (lam^2-Om) v:', chk == vs[1])
    for mu, p in P.items():
        M = sp.Matrix.hstack(*[p * vs[n] for n in sorted(vs)])
        print(mu, 'component span dim over n:', M.rank())
    pickle.dump((Qn,), open('master_kernel2_k7_Q%d.pkl' % lam, 'wb'))
