import pickle, sympy as sp
from sympy.polys.polyfuncs import symmetrize
a,ap,t,x1,x2,y1,y2=sp.symbols('a ap t x1 x2 y1 y2')
e1,e2,h1,h2=sp.symbols('e1 e2 h1 h2')
def to_e(expr):
    expr=sp.expand(expr)
    r=symmetrize(expr,[x1,x2],formal=True); assert r[1]==0
    (s1,_),(s2,_)=r[2]
    ex=sp.expand(r[0].subs({s1:e1,s2:e2},simultaneous=True))
    r=symmetrize(ex,[y1,y2],formal=True); assert r[1]==0
    (s1,_),(s2,_)=r[2]
    return sp.expand(r[0].subs({s1:h1,s2:h2},simultaneous=True))
if __name__=='__main__':
    A,B=pickle.load(open('mk_AB2.pkl','rb'))
    K=pickle.load(open('mk_Pn_1_1.pkl','rb'))
    Ae,Ke=to_e(A),to_e(K)
    print('K=',sp.factor(Ke)); print('terms A',len(Ae.args))
    pickle.dump((Ae,Ke),open('mk_AK_e.pkl','wb'))
    print(sp.factor(Ae))
