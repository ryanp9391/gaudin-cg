from jet_basis_sym import *
import sys
L=int(sys.argv[1]); lam=(L,0,0)
th=th_from_key('2.7_-3.5')
Yr=Y_state_rational(lam,'P',th)
t0=time.time()
rows=Y_from_operators(lam,(Fr(0),Fr(0),Fr(0)),th)
ok=all(rows[i][j]==Yr[i][j] for i in range(len(rows)) for j in range(len(rows[0])))
print('L',L,'zeta=0 (nilpotent companion twist): Y == <x0|nu nu(ops):',ok,'%.0fs'%(time.time()-t0))
