import sys, pickle, time
import sympy as sp
from mk2_fast import *
from mk2_opQ import Dlist, exp_series
n=5; L=int(sys.argv[1]); dmin=int(sys.argv[2]); dmax=int(sys.argv[3])
D=Dlist(4)
D[5]=sp.Integer(0)
Phi=exp_series(D,5)
rest_terms=sym_to_terms(Phi[5])
terms=terms_general(n,dmin,dmax)
print('unknowns',len(terms),flush=True)
if L==6:
    d=pickle.load(open('master_kernel2_N6.pkl','rb'))
    N=sp.zeros(784,784)
    for (i,j),v in d.items(): N[i,j]=v
else: N=None
t=time.time()
res=fit_general(n,L,rest_terms,terms,N)
print('total',time.time()-t)
if res[0] is None: print('NO SOLUTION')
else:
    sol,free,terms=res
    print('free columns',len(free))
    D5={}
    for pj,(v,deps) in sol.items():
        if v!=0: D5[terms[pj]]=sp.Rational(v.numerator,v.denominator)
    print('nonzero D5 terms',len(D5),'dependents?',any(deps for v,deps in sol.values()))
    pickle.dump(D5,open('master_kernel2_D5_gen_L%d.pkl'%L,'wb'))
    degs={}
    for (e,d_),v in D5.items(): degs[sum(d_)]=degs.get(sum(d_),0)+1
    print('terms by p-degree',degs)
