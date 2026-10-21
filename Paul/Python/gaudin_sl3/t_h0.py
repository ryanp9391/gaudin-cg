from jet_basis_sym import *
for L in (1,2,3,4):
    alpha,s,hw,Eg = dictionary(L)
    d=len(alpha)
    E12=sm_from_sympy(Eg[(1,2)]); E23=sm_from_sympy(Eg[(2,3)])
    for a in (Fr(2,5),Fr(-3,7)):
        e=[Fr(0)]*d; e[hw]=Fr(1)
        v=exp_nil(E12,a,40).rapply(e); v=exp_nil(E23,a,40).rapply(v)
        pred=[Fr(int(sp.Rational(s[p]).p),int(sp.Rational(s[p]).q))*a**(alpha[p][1]+2*alpha[p][2]) for p in range(d)]
        print(L,a,v==pred)
