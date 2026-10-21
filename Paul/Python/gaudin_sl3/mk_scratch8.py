import sys, pickle
import sympy as sp
l1=int(sys.argv[1])
Q=pickle.load(open('mk_symbol_%d.pkl'%l1,'rb'))
a,ap,t,x1,x2,y1,y2,e=sp.symbols('a ap t x1 x2 y1 y2 e')
P=sp.Poly(Q,t)
for n in range(P.degree()+1):
    c=P.coeff_monomial(t**n)
    # vanishing order along x1=y1
    def order(expr, sub):
        ex=sp.expand(expr.subs(sub,simultaneous=True))
        pe=sp.Poly(ex,e)
        ms=[m[0] for m in pe.monoms()]
        return min(ms) if ms else None
    o11=order(c,{x1:y1+e})
    o_all=order(c,{x1:y1+e}) 
    # both x1=y1 and x2=y2 
    o2=order(c,{x1:y1+e,x2:y2+e})
    o3=order(c,{x1:y1+e,x2:y1+e}) # x1=x2=y1
    print(l1,n,'ord(x1=y1):',o11,' ord(x1=y1,x2=y2):',o2,' ord(x1=x2=y1):',o3,flush=True)
