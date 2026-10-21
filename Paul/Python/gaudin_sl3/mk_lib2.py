import sympy as sp
from sympy.polys.polyfuncs import symmetrize
a,ap,t,x1,x2,y1,y2=sp.symbols('a ap t x1 x2 y1 y2')
e1,e2,h1,h2=sp.symbols('e1 e2 h1 h2')  # e: sym fns of x ; h: sym fns of y
def to_e(expr):
    r=symmetrize(sp.expand(expr),[x1,x2],formal=True)
    assert r[1]==0
    s1,s2=[s for s,_ in r[2]]
    ex=r[0].subs({s1:e1,s2:e2})
    ex=sp.expand(ex)
    r2=symmetrize(ex.subs({e1:e1,e2:e2}) if False else sp.expand(ex.subs({}) ) ,[y1,y2],formal=True)
    return r2
