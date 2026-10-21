from jet_basis_sym import *
import sys
L=int(sys.argv[1]); lam=(L,0,0)
coef,nt=Y_poly(lam, deg=int(sys.argv[2]) if len(sys.argv)>2 else None)
print('t-degree',len(coef)-1,'npts',nt)
Y0=coef[0]; print('det Y0',Y0.det() if Y0.shape[0]<=9 else Y0.rank())
A=Y0.inv()*coef[1]
d=gk.one(lam).d
Om=omega12(lam)
Omd=sp.zeros(d*d,d*d)
for i,row in Om.r.items():
    for c,v in row.items(): Omd[i,c]=sp.Rational(v.numerator,v.denominator)
print('A==c*Omega12?', [ (A[i,j]/Omd[i,j]) for i in range(d*d) for j in range(d*d) if Omd[i,j]!=0][:12])
print('A nonzero', sum(1 for x in A if x!=0), 'Omega nonzero', sum(1 for x in Omd if x!=0))
if len(coef)>2:
    print('exp test Y2==Y0 A^2/2', coef[2]==Y0*A*A/2)
pickle.dump((coef,A),open('jet_basis_sym_A%d.pkl'%L,'wb'))
