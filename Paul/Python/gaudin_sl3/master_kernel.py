"""M1: master kernel R^{l1}_{(kl),(k'l')}(x,y), Theta^{[l1,0,0]} = theta_kl(a) theta_k'l'(a') R, R = sum F_{mm'} x^m y^m', x=a t, y=a' t.
Python only. Stage flags: python -X utf8 master_kernel.py --stages TAB"""
import sys, os, pickle
from fractions import Fraction as Fr
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import general_kernel as G

def tabulate(l1):
    lam = (l1, 0, 0)
    T, rep = G.load_rep(lam)
    prs = G.pairs_of(lam)
    R = {}
    for i, p in enumerate(prs):
        for j, q in enumerate(prs):
            ep, eq = G.e_of(lam, p), G.e_of(lam, q)
            cc = G.ckl_closed(lam, *p) * G.ckl_closed(lam, *q)
            for (n, a, b), c in T[(i, j)].items():
                m, mp_ = a - ep, b - eq
                assert m + mp_ == n
                R[(p[0], p[1], q[0], q[1], m, mp_)] = Fr(c) / cc
    return R

def get_tab():
    f = os.path.join(HERE, 'master_kernel_tab.pkl')
    if os.path.exists(f):
        return pickle.load(open(f, 'rb'))
    tab = {l1: tabulate(l1) for l1 in (1, 2, 3, 4)}
    pickle.dump(tab, open(f, 'wb'))
    return tab

if __name__ == '__main__':
    tab = get_tab()
    for l1, R in tab.items():
        ms = sorted({(k[4], k[5]) for k in R})
        print(l1, len(R), 'support (m,mp):', ms)

# ============================================================================ closed orders 0,1,2 (symbol form), all lambda
import sympy as sp
_x1, _x2, _y1, _y2 = sp.symbols('x1 x2 y1 y2')

def symbol_orders(l1, a, ap, upto=2):
    """returns dict n -> sympy Poly(x1,x2,y1,y2) of Q_n, the t^n coefficient of the symbol Q^{(l1)}(xi;eta) defined by
       sum_{pq} Theta_pq (x1^k x2^l - x1^l x2^k)(y1^k' y2^l' - ...) = (x2-x1)(y2-y1) Q ;  a,ap exact rationals.
       Closed forms (verified l1<=4):  Pi=F_a(xi)G_ap(eta);  K = aa'/(a'-a) (Pi - F_ap(xi) G_a(eta));
       Q_0 = Pi^l, Q_1 = l^2 Pi^(l-1) K, Q_2 = l^2 (l-1) Pi^(l-2) (A + l K^2/2),
       A = -aa'(a-a')(e1-h1) Pi + K (K/2 + (a-a')(aa'(e1 h1 - e2 - h2) - a e1 - a' h1 + 2)),  e = sym fns of xi, h of eta."""
    a, ap = sp.Rational(a.numerator, a.denominator), sp.Rational(ap.numerator, ap.denominator)
    e1, e2, h1, h2 = _x1 + _x2, _x1 * _x2, _y1 + _y2, _y1 * _y2
    F = lambda c: (1 - c * _x1) * (1 - c * _x2)
    G = lambda c: (1 - c * _y1) * (1 - c * _y2)
    Pi = F(a) * G(ap)
    K = a * ap / (ap - a) * (Pi - F(ap) * G(a))
    A = -a * ap * (a - ap) * (e1 - h1) * Pi + K * (K / 2 + (a - ap) * (a * ap * (e1 * h1 - e2 - h2) - a * e1 - ap * h1 + 2))
    l = l1
    out = {0: Pi ** l}
    if upto >= 1:
        out[1] = l ** 2 * Pi ** (l - 1) * K
    if upto >= 2:
        out[2] = l ** 2 * (l - 1) * Pi ** (l - 2) * (A + sp.Rational(l, 2) * K ** 2)
    return {n: sp.Poly(sp.expand(v), _x1, _x2, _y1, _y2) for n, v in out.items()}

def theta_from_symbol_orders(l1, a, ap, t, upto=2):
    """Theta_pq(a,a',t) truncated at order `upto` in t, as dict {(p,q): Fraction}"""
    Qn = symbol_orders(l1, a, ap, upto)
    tot = None
    for n, P in Qn.items():
        term = P * sp.Poly((_x2 - _x1) * (_y2 - _y1), _x1, _x2, _y1, _y2) * sp.Rational(t.numerator, t.denominator) ** n
        tot = term if tot is None else tot + term
    n = l1 + 1
    prs = [(k, l) for k in range(n + 1) for l in range(k + 1, n + 1)]
    co = dict(tot.terms())
    res = {}
    for p in prs:
        for q in prs:
            c = co.get((p[0], p[1], q[0], q[1]), 0)
            res[(p, q)] = Fr(int(c.p), int(c.q)) if c != 0 else Fr(0)
    return res

def theta_from_pkl(l1, a, ap, t, maxorder=None):
    """Theta from the exact polynomial kernel (l1 <= 4), optionally truncated at t-order maxorder"""
    lam = (l1, 0, 0)
    T, rep = G.load_rep(lam)
    prs = G.pairs_of(lam)
    out = {}
    for i, p in enumerate(prs):
        for j, q in enumerate(prs):
            out[(p, q)] = sum(c * t ** n * a ** pp * ap ** qq for (n, pp, qq), c in T[(i, j)].items() if maxorder is None or n <= maxorder)
    return out

# ============================================================================ validation on states (no fit)
import mpmath as mp, time
import twovac_fit as tf

def load_pt(lam, tag, th, dps):
    path = G.pt2_path(lam, tag, th, dps)
    with open(path, 'rb') as f:
        pt = pickle.load(f)
    mp.mp.dps = dps
    return pt

def state_data(lam, pt, w1, w1p):
    """per state: F (list of m_i m_j products over pairs x pairs, row-major), y = c Qnum M1 M2 (exact jet LHS)"""
    sm = G.site_monos(lam, pt)
    R = G.one(lam)
    zeta = pt['zeta']
    om1 = [tf.mpf(x) for x in G.omega_inf(R, zeta, w1)]
    om2 = [tf.mpf(x) for x in G.omega_inf(R, zeta, w1p)]
    x0 = [tf.mpf(x) for x in G.x_cov(R, zeta)]
    phi = G.overlap2(pt, om1, om2)
    xx = G.overlap2(pt, x0, x0)
    c = G.cfac(lam, zeta, w1, w1p)
    F, y = [], []
    for k in range(len(pt['states'])):
        f1, M1 = sm[k][0]
        f2, M2 = sm[k][1]
        F.append([u * v for u in f1 for v in f2])
        y.append(c * (phi[k] / xx[k]) * M1 * M2)
    return F, y

def residuals(lam, pt, w1, w1p, theta_dict, F=None, y=None):
    """relative residual |y - sum Theta_ij F_ij| / |y| per state; Theta_dict {(p,q): Fraction}, canonical Theta; K = Theta (w1 w1p)^(2 l1)"""
    if F is None:
        F, y = state_data(lam, pt, w1, w1p)
    prs = G.pairs_of(lam)
    fac = (tf.mpf(Fr(w1)) * tf.mpf(Fr(w1p))) ** (2 * lam[0] - lam[1])
    Kv = [tf.mpf(theta_dict[(p, q)]) * fac for p in prs for q in prs]
    out = []
    for f, yy in zip(F, y):
        pred = mp.fsum(a * b for a, b in zip(Kv, f))
        out.append(abs(pred - yy) / abs(yy))
    return out, F, y
