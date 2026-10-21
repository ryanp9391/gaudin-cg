"""Universal operator N (inverse of M(1)=Y(1)Y(0)^-1) in the symbol-monomial basis."""
import pickle, sys, itertools
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy import QQ
import general_kernel as G
import master_kernel2 as m
E1,E2,H1,H2=m.E1,m.E2,m.H1,m.H2

def Nmat(L):
    Y = m.load_Y(L)
    n = Y[0].shape[0]
    Y1 = sum(Y, sp.zeros(n, n))
    d0 = DomainMatrix.from_Matrix(Y[0]).convert_to(QQ)
    d1 = DomainMatrix.from_Matrix(Y1).convert_to(QQ)
    return (d1 * d0.inv()).inv().to_Matrix()   # N = M(1)^-1 in jet-pair basis, rows (r,s), cols (p,q)

def iota(k,l,e1,e2):
    return sp.expand(e2**k*m.hcomp(l-k-1,e1,e2))

def sym_op(L, N, grade):
    """dict: input monomial (i,j,k,l) -> output poly (E1^i E2^j H1^k H2^l basis) for the grade-n part of N^T acting on symbols; rows (r,s) -> sum_pq N_{rs,pq} iota(p)iota(q)"""
    prs=G.pairs_of((L,0,0)); npr=len(prs)
    lev=[p[0]+p[1]-1 for p in prs]
    # monomial of pair r: iota(r) = sum monomials; we need N acting on iota(r)iota(s) (basis elements) -> convert to monomial basis by solving
    ex={}
    out={}
    for r in range(npr):
        for s in range(npr):
            row=r*npr+s
            acc=0
            for c in range(npr*npr):
                v=N[row,c]
                if v==0: continue
                p,q=divmod(c,npr)
                if lev[r]+lev[s]-lev[p]-lev[q]!=grade: continue
                acc+=v*iota(*prs[p],E1,E2)*iota(*prs[q],H1,H2)
            out[(r,s)]=sp.expand(acc)
    return out,prs
