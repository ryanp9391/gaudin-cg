import pickle, sympy as sp, master_kernel2 as m
gens=(m.E1,m.E2,m.H1,m.H2,m.a_,m.ap_)
lam=sp.Symbol('lam')
c=pickle.load(open('master_kernel2_P3coef.pkl','rb'))
K=sp.cancel(m.K_); a,ap=m.a_,m.ap_
D=ap-a
Kd=sp.expand(K*D)   # polynomial
e1,e2,h1,h2=m.E1,m.E2,m.H1,m.H2
Pi=m.Pi_
A=-a*ap*(a-ap)*(e1-h1)*Pi + K*(K/2 + (a-ap)*(a*ap*(e1*h1-e2-h2)-a*e1-ap*h1+2))
A=sp.cancel(A)
P1=lam**2*K
P2=lam**2*(lam-1)*(A+lam*K**2/2)
P3=lam**2*(lam-1)*sum(c[k]*lam**k for k in range(4))
# cumulants in units of (t/Pi)^n : W_n*Pi^n
W1=P1
W2=P2-P1**2/2
W3=P3-P1*P2+P1**3/3
for n,W in ((1,W1),(2,W2),(3,W3)):
    W=sp.expand(sp.cancel(sp.together(W*D**3)))   # clear denominators
    p=sp.Poly(W,lam)
    print(n,'lam-degree',p.degree(),'coeff factors:',[sp.factor(cc) for cc in p.all_coeffs()][:0])
    pickle.dump(W,open('master_kernel2_W%d.pkl'%n,'wb'))
    print([ (k, sp.Poly(p.coeff_monomial(lam**k),*gens).is_zero) for k in range(p.degree()+1)])
