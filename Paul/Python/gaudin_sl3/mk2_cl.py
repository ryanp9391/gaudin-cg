import pickle, sympy as sp, master_kernel2 as m
from mk2_coord import *
Dm=ap-a
W3=pickle.load(open('master_kernel2_W3.pkl','rb'))   # times Dm^3, poly in lam
p3=sp.Poly(W3,lam)
K=m.K_; Pi=m.Pi_
C=sp.cancel((-a*ap*(a-ap)*(E1-H1)*Pi + K*((a-ap)*(a*ap*(E1*H1-E2-H2)-a*E1-ap*H1+2))))
Dd=sp.cancel(p3.coeff_monomial(lam**4)/Dm**3)
Kc=tocoord(K); Cc=tocoord(C); Dc=tocoord(Dd); Pc=tocoord(Pi)
print('Pi',Pc)
z=a*ap/(ap-a)
# numeric a,a' to inspect
sub={a:sp.Rational(1,3),ap:sp.Rational(1,5)}
for nm,e in (('K',Kc),('C',Cc),('D',Dc)):
    print(nm,sp.factor(e.subs(sub)))
k1=Kc/Pc; 
# exp(S)/Pi coefficients
k=Kc; c=Cc; d=Dc; P=Pc
c2=sp.cancel((c/P**2+(k/P)**2/2).subs(sub)); c3=sp.cancel((d/P**3+k*c/P**3+(k/P)**3/6).subs(sub))
print('tau^2 coeff of exp(S)/Pi :',sp.factor(c2)); print('tau^3:',sp.factor(c3))
