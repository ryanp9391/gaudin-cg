import sys, pickle, time
import sympy as sp
from mk2_fast import *
n=int(sys.argv[1]); L=int(sys.argv[2]); dmin=int(sys.argv[3]); dmax=int(sys.argv[4])
D=X  # symbols
D1=X.D1; D2=X.getD2(); D3=X.sym_from_terms(pickle.load(open('master_kernel2_D3.pkl','rb')))
if n==4:
    rest=D1**4/24+D1**2*D2/2+D2**2/2+D1*D3
rest_terms=sym_to_terms(rest)
terms=terms_general(n,dmin,dmax)
print('unknowns',len(terms),flush=True)
t=time.time()
res=fit_general(n,L,rest_terms,terms)
print('total',time.time()-t)
if res[0] is None: print('NO SOLUTION')
else:
    sol,free,terms=res
    print('free columns',len(free))
    D4={}
    for pj,(v,deps) in sol.items():
        if v!=0: D4[terms[pj]]=sp.Rational(v.numerator,v.denominator)
    print('nonzero D%d terms'%n,len(D4),'dependents with free?',any(deps for v,deps in sol.values()))
    pickle.dump(D4,open('master_kernel2_D%d_gen.pkl'%n,'wb'))
    degs={}
    for (e,d),v in D4.items(): degs[sum(d)]=degs.get(sum(d),0)+1
    print('terms by p-degree',degs)
