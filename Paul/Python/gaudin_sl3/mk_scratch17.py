import pickle, sympy as sp, random
a,ap,t,x1,x2,y1,y2=sp.symbols('a ap t x1 x2 y1 y2')
vars_=(x1,x2,y1,y2,a,ap)
P=lambda n,l: sp.Poly(pickle.load(open('mk_Pn_%d_%d.pkl'%(n,l),'rb')),*vars_)
ps=[P(2,l) for l in (2,3,4)]
pts=[[sp.Rational(random.randint(-9,9),random.randint(1,7)) for _ in range(6)] for _ in range(6)]
M=sp.Matrix([[p.eval(dict(zip(vars_,pt))) if False else p.as_expr().subs(dict(zip(vars_,pt))) for p in ps] for pt in pts])
print(M.nullspace())
ps3=[P(3,l) for l in (3,4)]
