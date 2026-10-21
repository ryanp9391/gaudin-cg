import pickle, sympy as sp, master_kernel2 as m
res=pickle.load(open('master_kernel2_Pn.pkl','rb'))
gens=(m.E1,m.E2,m.H1,m.H2,m.a_,m.ap_)
lam=sp.Symbol('lam')
P3={}
for L in (3,4,5): P3[L]=res[(L,3)][0]
P3[2]=sp.expand(res[(2,3)][0]*m.Pi_)
# P3 = lam^2 (lam-1) (c0 + c1 lam + c2 lam^2 + c3 lam^3)
Ls=[2,3,4,5]
V=sp.Matrix([[L**k for k in range(4)] for L in Ls])
rhs=[sp.expand(P3[L]/(L**2*(L-1))) for L in Ls]
Vi=V.inv()
c=[sp.expand(sum(Vi[k,i]*rhs[i] for i in range(4))) for k in range(4)]
K=m.K_
Kn=sp.cancel(K)
print('c3 - K^3/6 =', sp.simplify(sp.expand(c[3]*(m.ap_-m.a_)**3 - sp.expand(sp.cancel(K*(m.ap_-m.a_))**3)/6)))
pickle.dump(c,open('master_kernel2_P3coef.pkl','wb'))
for k in range(4): print(k,len(sp.Poly(c[k],*gens).terms()))
