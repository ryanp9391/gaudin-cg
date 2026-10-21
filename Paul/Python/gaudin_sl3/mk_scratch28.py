import pickle, sympy as sp
a,ap,u,up,v,vp=sp.symbols('a ap u up v vp')
w,wp,U,Up,V,Vp,nu=sp.symbols('w wp U Up V Vp nu')
A,K=pickle.load(open('mk_AK_u.pkl','rb'))
def conv(X):
    X=X.subs({a:1/w,ap:1/wp,u:U/w**2,up:Up/wp**2,v:V/wp**2,vp:Vp/w**2},simultaneous=True)
    return sp.factor(sp.simplify(X*(w*wp)**4))
At=conv(A)
print('A~ =',At)
# subtract: write A~ = alpha*K~^2 + ...
Kt=sp.factor(sp.simplify(K.subs({a:1/w,ap:1/wp,u:U/w**2,up:Up/wp**2,v:V/wp**2,vp:Vp/w**2},simultaneous=True)*(w*wp)**2))
print('K~',Kt)
