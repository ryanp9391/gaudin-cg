"""M3: exact Theta^{[5,0,0]} from M2's operator-derived Y(t); grading Y(t)=t^rho Y(1) t^-kappa.  python -X utf8 master_kernel2.py inv5"""
import sys, os, pickle, time
from fractions import Fraction as Fr
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy import QQ, ZZ
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import jet_basis_sym as J

def load_Y(L):
    return pickle.load(open(os.path.join(HERE, 'jet_basis_sym_Ycoef_ops%d.pkl' % L), 'rb'))

def grading(L, coef):
    """kappa_c (ket weight alpha2+2alpha3 summed over two sites), rho_r from Y_k pattern: rho_r - kappa_c = k"""
    alpha, s, hw, Eg = J.dictionary(L); d = len(alpha)
    w = [alpha[p][1] + 2 * alpha[p][2] for p in range(d)]
    kappa = [w[p] + w[q] for p in range(d) for q in range(d)]
    return kappa

if __name__ == '__main__':
    pass

def theta_T(L):
    """exact Theta^{[L,0,0]} as general_kernel-format dict T[(i,j)] = {(n,p,q): Fraction} (n=t power, p=a power, q=a' power), from Y1^-1 and the grading"""
    D = pickle.load(open(os.path.join(HERE, 'master_kernel2_Y1inv%d.pkl' % L), 'rb'))
    alpha, s, hw, Eg = J.dictionary(L); d = len(alpha)
    w = [alpha[p][1] + 2 * alpha[p][2] for p in range(d)]
    nprs = (L + 2) * (L + 1) // 2
    assert d * d == nprs * nprs
    out = {}
    for (c, r), v in D['Y1inv'].items():
        p, q = divmod(c, d)
        ev = sp.Rational(s[p]) * sp.Rational(s[q])
        coef = ev * v
        key = (D['kappa'][c] - D['rho'][r], w[p], w[q])
        i, j = divmod(r, nprs)
        dd = out.setdefault((i, j), {})
        dd[key] = dd.get(key, 0) + coef
    return {k: {kk: Fr(int(vv.p), int(vv.q)) for kk, vv in v.items() if vv != 0} for k, v in out.items()}

# ============================================================ symbol in symmetric variables
import general_kernel as G
E1, E2, H1, H2, a_, ap_ = sp.symbols('E1 E2 H1 H2 a ap')

def hcomp(m, e1, e2):
    """complete homogeneous h_m(x1,x2) in e1,e2:  h_m = e1 h_{m-1} - e2 h_{m-2}"""
    h = [sp.Integer(1), e1]
    for k in range(2, m + 1):
        h.append(sp.expand(e1 * h[-1] - e2 * h[-2]))
    return h[m] if m >= 0 else sp.Integer(0)

def get_T(L):
    if L == 5:
        return pickle.load(open(os.path.join(HERE, 'master_kernel2_rep_500.pkl'), 'rb'))
    if L == 6:
        return pickle.load(open(os.path.join(HERE, 'master_kernel2_rep_600.pkl'), 'rb'))
    return G.load_rep((L, 0, 0))[0]

def symbol_orders_exact(L, T=None, nmax=None):
    """dict n -> sympy expr in E1,E2,H1,H2,a,ap : Q_n, with  sum Theta_pq (x^k..)  = (x2-x1)(y2-y1) Q"""
    T = T or get_T(L)
    prs = G.pairs_of((L, 0, 0))
    hx = [hcomp(m, E1, E2) for m in range(L + 2)]
    hy = [hcomp(m, H1, H2) for m in range(L + 2)]
    Q = {}
    for i, p in enumerate(prs):
        fx = E2 ** p[0] * hx[p[1] - p[0] - 1]
        for j, q in enumerate(prs):
            fy = H2 ** q[0] * hy[q[1] - q[0] - 1]
            for (n, pp, qq), c in T[(i, j)].items():
                if nmax is not None and n > nmax:
                    continue
                Q[n] = Q.get(n, 0) + sp.Rational(c.numerator, c.denominator) * a_ ** pp * ap_ ** qq * fx * fy
    return Q

Pi_ = (1 - a_ * E1 + a_ ** 2 * E2) * (1 - ap_ * H1 + ap_ ** 2 * H2)
K_ = a_ * ap_ / (ap_ - a_) * (Pi_ - (1 - ap_ * E1 + ap_ ** 2 * E2) * (1 - a_ * H1 + a_ ** 2 * H2))
