from jet_basis_sym import *
import sys
L=int(sys.argv[1]); lam=(L,0,0)
zero=(Fr(0),Fr(0),Fr(0))
nu,cond=nu_formulas(L)
print('cond',cond, flush=True)
pickle.dump(nu,open('jet_basis_sym_nu%d.pkl'%L,'wb'))
ths=[(Fr(2,7),Fr(2,7)-Fr(35,k)) for k in (31,42,53,64,75,86,97)]  # theta12 = 35/k -> t = k/35 ... all distinct
res={}
for th in ths[:L+3]:
    t0=time.time()
    rows=Y_from_operators(lam,zero,th,nu=nu)
    t=1/(th[0]-th[1])
    res[t]=rows
    print('t',t,'%.0fs'%(time.time()-t0),flush=True)
    pickle.dump(res,open('jet_basis_sym_Yops%d.pkl'%L,'wb'))
