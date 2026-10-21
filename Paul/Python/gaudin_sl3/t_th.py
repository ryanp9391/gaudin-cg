from jet_basis_sym import *
import sys
L=int(sys.argv[1]); lam=(L,0,0)
coef,A=pickle.load(open('jet_basis_sym_A%d.pkl'%L,'rb'))
th=Theta_from_Y(L,coef)
T,rep=gk.load_rep(lam)
prs=gk.pairs_of(lam); n=len(prs)
bad=0
for i in range(n):
    for j in range(n):
        ref=sum(sp.Rational(c.numerator,c.denominator)*t_s**nn*a_s**p*ap_s**q for (nn,p,q),c in T[(i,j)].items())
        if sp.expand(th[i*n+j]-ref)!=0: bad+=1; 
        if i*n+j<12 or bad<3 and sp.expand(th[i*n+j]-ref)!=0: print((i,j),th[i*n+j],'|',ref)
print('mismatches',bad,'of',n*n)
pickle.dump(th,open('jet_basis_sym_Theta%d.pkl'%L,'wb'))
