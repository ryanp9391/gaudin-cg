import pickle, sympy as sp
a,ap,t,x1,x2,y1,y2=sp.symbols('a ap t x1 x2 y1 y2')
P=lambda n,l: pickle.load(open('mk_Pn_%d_%d.pkl'%(n,l),'rb'))
print('n=1 ratios', sp.cancel(P(1,2)/P(1,1)), sp.cancel(P(1,3)/P(1,1)))
r=sp.cancel(P(2,3)/P(2,2)); print('n=2 ratio 3/2', r)
