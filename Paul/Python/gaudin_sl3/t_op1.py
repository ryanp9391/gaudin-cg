from jet_basis_sym import *
from gaudin_chain import Chain, Operators, gc
import mpmath as mp
L=1; lam=(1,0,0)
tag='P'; th=th_from_key('2.7_-3.5'); zeta=gk.ZETAS2[tag]
pt=gk.get_point2(lam,tag,zeta,th)
mp.mp.dps=pt['dps']
R=gk.one(lam); d=R.d
M=pinv_rows(R,zeta)
cols=comp_states(lam,pt,M)
ch=Chain(lam=lam,theta=th)
gcm=gc(zeta)   # sympy matrix companion of roots zeta
for variant,gmat in (('gc',[[Fr(int(gcm[i,j].p),int(gcm[i,j].q)) for j in range(3)] for i in range(3)]),):
    for transpose in (False,True):
        ops=Operators(ch,gmat,transpose=transpose)
        C=ops.fit()
        print(variant,'transpose',transpose,'fit res',C['res2'],C['res3'])
        # H1 operator = C2[1], a1=C3[2], b1=C3[1]; check eigen on companion kets
        for name,(k,b) in (('H1',(2,1)),('H2',(2,3)),('a1',(3,2)),('b1',(3,1))):
            op=C[k][b]
            worst=0
            for st,col in zip(pt['states'],cols):
                v=[mp.mpf(x) for x in col]
                w=op.apply(v)
                # eigenvalue from Laurent data
                lam_ev={('H1'):st['lam2'][1],'H2':st['lam2'][3],'a1':st['lam3'][2],'b1':st['lam3'][1]}[name]
                num=max(abs(w[i]-lam_ev*v[i]) for i in range(len(v)))
                worst=max(worst,num)
            print('  ',name,'eigen residual',mp.nstr(worst,3))
