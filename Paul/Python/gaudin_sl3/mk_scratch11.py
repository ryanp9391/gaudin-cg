import sys
sys.path.insert(0,'.')
import sympy as sp
from fractions import Fraction as Fr
import general_kernel as G
from mk_lib import Nmat, wedge2
a,ap,t=sp.symbols('a ap t')
for l1 in (1,2,3):
    n=l1+1; lam=(l1,0,0)
    T,rep=G.load_rep(lam); prs=G.pairs_of(lam)
    N2=wedge2(Nmat(n),prs)
    Th=sp.Matrix(len(prs),len(prs),lambda i,j: sum(sp.Rational(c.numerator,c.denominator)*t**nn*a**p*ap**q for (nn,p,q),c in T[(i,j)].items()))
    V=N2*Th+Th*N2.T-Th.diff(a)-Th.diff(ap)
    W=Th.diff(t)
    # ratio test
    ratios=set()
    for i in range(V.shape[0]):
        for j in range(V.shape[1]):
            if W[i,j]!=0:
                ratios.add(sp.simplify(V[i,j]/W[i,j]))
            elif V[i,j]!=0: ratios.add('inf')
    print(l1,len(ratios),list(ratios)[:4],flush=True)
