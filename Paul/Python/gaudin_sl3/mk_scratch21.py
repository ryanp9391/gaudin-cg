import pickle, sympy as sp
a,ap,e1,e2,h1,h2=sp.symbols('a ap e1 e2 h1 h2')
Ae,Ke=pickle.load(open('mk_AK_e.pkl','rb'))
sol=sp.solve(Ke,h2)[0]
Ar=sp.factor(sp.together(Ae.subs(h2,sol)))
print(Ar)
Pi=(1-a*e1+a**2*e2)*(1-ap*h1+ap**2*h2)
print('Pi on K=0:',sp.factor(sp.together(Pi.subs(h2,sol))))
