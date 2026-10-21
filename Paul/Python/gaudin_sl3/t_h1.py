from jet_basis_sym import *
for L in (1,2):
    lam=(L,0,0); R=gk.one(lam); zeta=(Fr(3,2),Fr(5,4),Fr(2,3))
    M=pinv_rows(R,zeta)
    e=[Fr(0)]*R.d; e[R.hw_index]=Fr(1)
    print('x', vecmat(e,M)==gk.x_cov(R,zeta))
    alpha,s,hw,Eg=dictionary(L); assert hw==R.hw_index
    for w in (Fr(5,2),Fr(7,4)):
        a=1/w
        sf=[Fr(int(sp.Rational(x).p),int(sp.Rational(x).q)) for x in s]
        Oc=[sf[p]*a**(alpha[p][1]+2*alpha[p][2]) for p in range(R.d)]
        lhs=vecmat(Oc,M); rhs=gk.omega_inf(R,zeta,w)
        ratios=set(l/r for l,r in zip(lhs,rhs) if r!=0)
        z1,z2,z3=zeta
        kap=(w**2/((z2-w)*(z3-w)))**L
        print(L,w,ratios, 'zero pattern ok', all((l==0)==(r==0) for l,r in zip(lhs,rhs)), 'Delta1^L', kap, 1/kap)
