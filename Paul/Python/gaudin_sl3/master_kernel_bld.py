import sys, time, inspect, pickle, os
import general_kernel as gk
from fractions import Fraction as Fr
# patched build_point2: SVD null threshold scaled to the working precision (G1's 1e-50 is for 80 digits)
src = inspect.getsource(gk.build_point2).replace('Q.TOL = mp.mpf(10) ** -50', 'Q.TOL = mp.mpf(10) ** -(5 * dps // 8)')
exec(src, gk.__dict__)
lam=tuple(int(c) for c in sys.argv[1]); tag=sys.argv[2]; dps=int(sys.argv[3])
th=gk.TH0
if len(sys.argv)>4: th=gk.theta_pair(Fr(sys.argv[4]))
path=gk.pt2_path(lam,tag,th,dps)
if os.path.exists(path):
    print('exists',path); sys.exit()
t=time.time()
pt=gk.build_point2(lam,gk.ZETAS2[tag],th,dps)
pickle.dump(pt,open(path+'.tmp','wb')); os.replace(path+'.tmp',path)
print('done',lam,tag,dps,th,'bad',sum(1 for s in pt['states'] if not s['ok']),'times',pt['times'],time.time()-t,flush=True)
