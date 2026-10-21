import sys; sys.path.insert(0,'.')
import sympy as sp
from fractions import Fraction as Fr
import general_kernel as G
from mk_lib import wedge2
def Fmat(n,c=1):
    M=sp.zeros(n+1,n+1)
    for k in range(1,n+1): M[k-1,k]=c*k
    return M
def Tnum(T,prs,a,ap,t):
    def ev(ent):
        v=G.eval_T(ent,Fr(a),Fr(ap),Fr(t)); return sp.Rational(v.numerator,v.denominator)
    return sp.Matrix(len(prs),len(prs),lambda i,j: ev(T[(i,j)]))
