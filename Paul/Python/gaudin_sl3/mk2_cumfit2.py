import pickle, sympy as sp, sys
import master_kernel2 as m
gens=(m.E1,m.E2,m.H1,m.H2)
d=pickle.load(open(sys.argv[1],'rb'))
Ls=sorted(d)
for n in range(1,max(d[Ls[0]])+1):
    ks=list(range(2,n+2))
    if len(Ls)<len(ks): continue
    use=Ls[:len(ks)]; chk=Ls[len(ks):]
    Vi=sp.Matrix([[L**k for k in ks] for L in use]).inv()
    co=[sp.expand(sum(Vi[ki,i]*d[use[i]][n] for i in range(len(use)))) for ki in range(len(ks))]
    ok=[sp.expand(sum(co[ki]*L**ks[ki] for ki in range(len(ks)))-d[L][n])==0 for L in chk]
    print('n',n,'unknown lam^2..lam^%d'%(n+1),'checks on',chk,ok)
