from jet_basis_sym import *
import sys
L=int(sys.argv[1]); lam=(L,0,0)
for tag,key in (('P','2.7_-3.5'),('Pp','2.7_-3.5'),('X1','2.7_-3.5')):
    if L>=3 and tag!='P': 
        pass
    th=th_from_key(key); zeta=gk.ZETAS2[tag]
    t0=time.time()
    rows=Y_from_operators(lam,zeta,th)
    Yr=Y_state_rational(lam,tag,th)
    ok=all(rows[i][j]==Yr[i][j] for i in range(len(rows)) for j in range(len(rows[0])))
    print('L',L,tag,key,'Y == <x0|nu nu(ops):',ok,'%.0fs'%(time.time()-t0),flush=True)
