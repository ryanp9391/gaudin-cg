import sys
sys.path.insert(0,'.')
from master_kernel import *
import sympy as sp
tab=get_tab()
l1=int(sys.argv[1]); n=l1+1
R=tab[l1]
prs=[(k,l) for k in range(n+1) for l in range(k+1,n+1)]
ms=sorted({(k[4],k[5]) for k in R},key=lambda x:(x[0]+x[1],x))
byN={}
for (m,mp_) in ms:
    M=sp.Matrix(len(prs),len(prs),lambda i,j: sp.Rational(R.get((prs[i][0],prs[i][1],prs[j][0],prs[j][1],m,mp_),0).numerator, R.get((prs[i][0],prs[i][1],prs[j][0],prs[j][1],m,mp_),Fr(1)).denominator))
    byN.setdefault(m+mp_,[]).append(((m,mp_),M.rank()))
for N in sorted(byN): print(N,byN[N])
