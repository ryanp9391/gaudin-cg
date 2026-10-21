import pickle, sympy as sp
a,ap,t,x1,x2,y1,y2=sp.symbols('a ap t x1 x2 y1 y2')
vars_=(x1,x2,y1,y2,a,ap)
P=lambda n,l: sp.Poly(pickle.load(open('mk_Pn_%d_%d.pkl'%(n,l),'rb')),*vars_)
K=P(1,1)  # lam=1: Q_1 = K
p2,p3,p4=P(2,2),P(2,3),P(2,4)
u2=p2*sp.Rational(1,4); u3=p3*sp.Rational(1,18); u4=p4*sp.Rational(1,48)
B=u3-u2; A=u2-2*B
print('check', (u4-(A+4*B)).is_zero)
print('terms A,B,K',len(A.terms()),len(B.terms()),len(K.terms()))
q,r=sp.div(B,K*K); print('B/K^2 rem zero?',r.is_zero)
q,r=sp.div(A,K*K); print('A/K^2 rem zero?',r.is_zero)
Pi=sp.Poly(sp.expand((1-a*x1)*(1-a*x2)*(1-ap*y1)*(1-ap*y2)),*vars_)
for nm,X in (('A',A),('B',B)):
    # try X = c K^2 + Pi * M
    pass
pickle.dump((A.as_expr(),B.as_expr()),open('mk_AB2.pkl','wb'))
print(sp.factor(B.as_expr()) if len(B.terms())<80 else 'B big')
