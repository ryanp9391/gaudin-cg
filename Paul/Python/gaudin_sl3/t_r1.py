from jet_basis_sym import *
import sys
L=int(sys.argv[1]); lam=(L,0,0)
coef,A=pickle.load(open('jet_basis_sym_A%d.pkl'%L,'rb'))
th=pickle.load(open('jet_basis_sym_Theta%d.pkl'%L,'rb'))
R=gk.one(lam); d=R.d
Om=omega12(lam)
Omd=sp.zeros(d*d,d*d)
for i,row in Om.r.items():
    for c,v in row.items(): Omd[i,c]=sp.Rational(v.numerator,v.denominator)
C2=lambda l: sum(l[i]*(l[i]+4-2*(i+1)) for i in range(3))
mus=[(2*L-k,k,0) for k in range(L+1)]
eig=[sp.Rational(C2(m)-2*C2((L,0,0)),2) for m in mus]
print('eigs',eig)
I=sp.eye(d*d)
P=[]
for k in range(L+1):
    M=I
    for j in range(L+1):
        if j!=k: M=M*(Omd-eig[j]*I)/(eig[k]-eig[j])
    P.append(M)
print('projector check', [ (P[k]*P[k]-P[k]).is_zero_matrix for k in range(L+1)], sum(P,sp.zeros(d*d,d*d))==I)
Y0=coef[0]
nn=len(gk.pairs_of(lam))
vals={a_s:sp.Rational(3,7),ap_s:sp.Rational(-5,4),t_s:sp.Rational(2,9)}
Ehat=(th*Y0).subs(vals)
E=ev_cov(L).subs(vals)
cols=sp.Matrix.hstack(*[ (E*P[k]).T for k in range(L+1)])
sol=sp.linsolve((cols,Ehat.T),*sp.symbols('r0:%d'%(L+1)))
print('Ehat in span E P_k?',sol)
