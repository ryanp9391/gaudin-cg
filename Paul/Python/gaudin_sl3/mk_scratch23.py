import pickle, sympy as sp
a,ap,e1,e2,h1,h2=sp.symbols('a ap e1 e2 h1 h2')
u,up,v,vp=sp.symbols('u up v vp')
Ae,Ke=pickle.load(open('mk_AK_e.pkl','rb'))
s1=sp.solve([1-a*e1+a**2*e2-u,1-ap*e1+ap**2*e2-up],[e1,e2],dict=True)[0]
s2=sp.solve([1-ap*h1+ap**2*h2-v,1-a*h1+a**2*h2-vp],[h1,h2],dict=True)[0]
sub={**s1,**s2}
conv=lambda X: sp.factor(sp.together(X.subs(sub,simultaneous=True)))
print('K =',conv(Ke))
print('A =',conv(Ae))
pickle.dump((conv(Ae),conv(Ke)),open('mk_AK_u.pkl','wb'))
