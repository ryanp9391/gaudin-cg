import sys
sys.path.insert(0,'.')
from master_kernel import *
import sympy
tab=get_tab()
l1=int(sys.argv[1]); m=int(sys.argv[2]); mp_=int(sys.argv[3])
R=tab[l1]; n=l1+1
prs=[(k,l) for k in range(n+1) for l in range(k+1,n+1)]
print('rows p, cols q')
for p in prs:
    print(p, [str(R.get((p[0],p[1],q[0],q[1],m,mp_),0)) for q in prs])
