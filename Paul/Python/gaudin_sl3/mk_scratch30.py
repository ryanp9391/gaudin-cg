import pickle, sympy as sp, sys
z,w,wp=sp.symbols('z w wp'); p0,p1,q0,q1=sp.symbols('p0 p1 q0 q1')
l=int(sys.argv[1])
out=pickle.load(open('mk_Qw_%d.pkl'%l,'rb'))
Q=0
for (a,b,c,d),co in out.items():
    Q+=sp.factor(co.subs(w,wp+1))*(p0+p1)**a*p0**b*q0**c*(q0+q1)**d
Q=sp.expand(Q)
P=sp.Poly(Q,p0,p1,q0,q1)
print('nterms',len(P.terms()))
pickle.dump(Q,open('mk_Qp_%d.pkl'%l,'wb'))
if l<=2:
    for m,c in sorted(P.terms()):
        print(m,sp.factor(c))
