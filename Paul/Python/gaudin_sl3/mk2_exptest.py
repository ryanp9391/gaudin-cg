import sympy as sp, pickle
from mk2_expfit import *
D2=getD2()
x=D1
D2pred=sp.Rational(1,2)*sum(sp.diff(x,p)*sp.diff(x,q) for p,q in zip(PV,V))
print('D2 - 1/2 sum x_p x_q =', sp.expand(D2-D2pred))
# compare also with D2 alone structure
print(sp.factor(sp.expand(D2)))
