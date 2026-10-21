import pickle, sympy as sp, itertools, random
a,ap,t,x1,x2,y1,y2=sp.symbols('a ap t x1 x2 y1 y2')
Q={0:sp.Integer(1)}
for l in (1,2,3,4):
    Q[l]=pickle.load(open('mk_symbol_%d.pkl'%l,'rb'))
vars_={'a':a,'ap':ap,'t':t}
def val(expr,pt): return expr.subs(pt)
random.seed(3)
def rpt():
    return {a:sp.Rational(random.randint(1,9),random.randint(2,7)),ap:sp.Rational(random.randint(1,9),random.randint(2,7)),t:sp.Rational(random.randint(1,9),random.randint(2,7)),
            x1:sp.Rational(random.randint(-9,9),random.randint(1,5)),x2:sp.Rational(random.randint(-9,9),random.randint(1,5)),y1:sp.Rational(random.randint(-9,9),random.randint(1,5)),y2:sp.Rational(random.randint(-9,9),random.randint(1,5))}
p1=rpt(); p2=dict(p1); 
for s in (x1,x2,y1,y2): p2[s]=sp.Rational(random.randint(-9,9),random.randint(1,5))
ops=[('a','ap'),('t','t'),('a','a'),('ap','ap'),('a','t'),('ap','t')]
for (d1,d2) in ops:
    for l in (1,2,3):
        res=[]
        for pt in (p1,p2):
            T=Q[l]; D1=sp.diff(T,vars_[d1]); D2=sp.diff(T,vars_[d2]); D12=sp.diff(T,vars_[d1],vars_[d2])
            num=(T*D12-D1*D2).subs(pt)
            den=(Q[l+1]*Q[l-1]).subs(pt)
            res.append(sp.nsimplify(num/den))
        print(d1,d2,l,res,flush=True)
