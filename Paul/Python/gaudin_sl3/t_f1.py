from jet_basis_sym import *
import sys
n=int(sys.argv[1]); Ls=[int(x) for x in sys.argv[2].split(',')]
A_by={}
for L in Ls:
    coef,A=pickle.load(open('jet_basis_sym_A%d.pkl'%L,'rb'))
    Y0i=coef[0].inv()
    A_by[L]=Y0i*coef[n]
shifts={1:[(-1,1,0),(0,-1,1)],2:[(-2,2,0),(-1,0,1),(0,-2,2)]}[n]
mo=int(sys.argv[3]) if len(sys.argv)>3 else n
for sh in shifts:
    t0=time.time()
    sol,rank,cons,npairs,neq=fit_A(A_by,sh,mo,mo,Ls)
    print('shift',sh,'unknowns',npairs,'eqs',neq,'rank',rank,'consistent',cons,'nterms',len(sol),'%.0fs'%(time.time()-t0),flush=True)
    if cons and len(sol)<40:
        for k,v in sol.items(): print('   ',v,k)
