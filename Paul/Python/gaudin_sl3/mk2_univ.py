import sys, pickle, time
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy import QQ
import general_kernel as G
import master_kernel2 as m
def Mmat(L):
    Y = m.load_Y(L)
    n = Y[0].shape[0]
    Y1 = sum(Y, sp.zeros(n, n))
    d0 = DomainMatrix.from_Matrix(Y[0]).convert_to(QQ)
    d1 = DomainMatrix.from_Matrix(Y1).convert_to(QQ)
    M = d1 * d0.inv()
    return M
def restr(M, L, Lp):
    """restrict M^{(L)} (pairs of labels<=L+1, row-major pairs x pairs) to pairs with labels <= Lp+1"""
    prs = G.pairs_of((L, 0, 0)); n = len(prs)
    keep = [i for i, p in enumerate(prs) if p[1] <= Lp + 1]
    sel = [i * n + j for i in keep for j in keep]
    return sel
if __name__ == '__main__':
    Ms = {}
    for L in (3, 4, 5):
        t = time.time(); Ms[L] = Mmat(L); print(L, 'M built', time.time() - t, flush=True)
    for L in (3, 4, 5):
        for Lp in range(3, L):
            sel = restr(Ms[L], L, Lp)
            A = Ms[L].to_Matrix()[sel, sel]
            B = Ms[Lp].to_Matrix()
            ok = (A == B)
            Ai = Ms[L].inv().to_Matrix()[sel, sel]; Bi = Ms[Lp].inv().to_Matrix()
            print('M^(%d)|_%d == M^(%d):' % (L, Lp, Lp), ok, ' inverse restricts:', Ai == Bi, flush=True)
