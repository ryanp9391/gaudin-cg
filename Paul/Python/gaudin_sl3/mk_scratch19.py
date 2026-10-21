import pickle, sympy as sp, random
a,ap,t,x1,x2,y1,y2=sp.symbols('a ap t x1 x2 y1 y2')
vars_=(x1,x2,y1,y2,a,ap)
A,B=pickle.load(open('mk_AB2.pkl','rb'))
K=pickle.load(open('mk_Pn_1_1.pkl','rb'))
rs=set()
for _ in range(6):
    vals={x2:sp.Rational(random.randint(-9,9),random.randint(1,7)),y1:sp.Rational(random.randint(-9,9),random.randint(1,7)),y2:sp.Rational(random.randint(-9,9),random.randint(1,7)),a:sp.Rational(random.randint(1,9),random.randint(1,7)),ap:sp.Rational(random.randint(1,9),random.randint(1,7))}
    vals[x1]=1/vals[a]
    rs.add(sp.nsimplify(A.subs(vals)/K.subs(vals)**2))
print(rs)
