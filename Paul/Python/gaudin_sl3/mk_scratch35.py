import pickle, sympy as sp, sys
e1,e2,h1,h2,a,ap,t=sp.symbols('e1 e2 h1 h2 a ap t')
x=[sp.Integer(1),e1,e2]; y=[sp.Integer(1),h1,h2]
# NOTE: x_0 = 1 is dehomogenised: restore via Euler: use homogeneous x0,y0 symbols
x0,y0=sp.symbols('x0 y0')
X=[x0,e1,e2]; Y=[y0,h1,h2]
def Omega(f):
    out=0
    for i in range(3):
        for j in range(3):
            out+=X[i]*sp.diff(sp.diff(f,X[j]),Y[i])*Y[j] if False else X[i]*Y[j]*sp.diff(f,X[j],Y[i])
    return sp.expand(out)
l=int(sys.argv[1])
av,apv=sp.Rational(2,5),sp.Rational(3,7)
ell=lambda c,V: V[0]-c*V[1]+c**2*V[2]
v=sp.expand(ell(av,X)**l*ell(apv,Y)**l)
vecs=[v]
for k in range(l+1): vecs.append(Omega(vecs[-1]))
def dehom(f): return sp.expand(f.subs({x0:1,y0:1}))
# Q_n from data (symbol in xi,eta) -> e,h variables
from mk_scratch20 import to_e as to_e_
x1,x2,y1,y2=sp.symbols('x1 x2 y1 y2')
Q=pickle.load(open('mk_symbol_%d.pkl'%l,'rb'))
Qe=to_e_(Q.subs({sp.Symbol('a'):av,sp.Symbol('ap'):apv}))
Qe=sp.expand(Qe)
Qn=[sp.expand(Qe.coeff(t,n)) for n in range(0,4*(l-1)+2)] if l>1 else [sp.expand(Qe.coeff(t,n)) for n in range(2)]
# span test
basis=[dehom(w) for w in vecs[:l+1]]
mons=sorted({m for b in basis+Qn for m in sp.Poly(b,e1,e2,h1,h2).monoms()})
def vec(b):
    P=sp.Poly(b,e1,e2,h1,h2); return [P.coeff_monomial(m) for m in mons]
B=sp.Matrix([vec(b) for b in basis]).T
print('rank basis',B.rank(),'of',l+1)
for n,qn in enumerate(Qn):
    M=B.row_join(sp.Matrix(vec(qn)))
    print(n,'in span' if M.rank()==B.rank() else 'NOT in span')
