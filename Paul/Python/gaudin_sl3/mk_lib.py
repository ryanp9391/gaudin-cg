import sys
sys.path.insert(0,'.')
import sympy as sp
from fractions import Fraction as Fr
import general_kernel as G
from master_kernel import *
def Nmat(n):
    M=sp.zeros(n+1,n+1)
    for k in range(n): M[k+1,k]=-(n-k)
    return M
def wedge2(M,prs):
    idx={p:i for i,p in enumerate(prs)}
    out=sp.zeros(len(prs),len(prs))
    for (k,l),j in idx.items():
        for r in range(M.shape[0]):
            if M[r,k]!=0 and r!=l:
                s=1; a,b=r,l
                if a>b: a,b=b,a; s=-1
                out[idx[(a,b)],j]+=s*M[r,k]
            if M[r,l]!=0 and r!=k:
                s=1; a,b=k,r
                if a>b: a,b=b,a; s=-1
                out[idx[(a,b)],j]+=s*M[r,l]
    return out
