import sys, time, pickle
import master_kernel2 as m, sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy import QQ
L = int(sys.argv[1])
Y = m.load_Y(L)
n = Y[0].shape[0]
kappa = m.grading(L, Y)
rho = [None]*n; bad = 0
for k, M in enumerate(Y):
    for (i, j), v in M.todok().items():
        if v == 0: continue
        r = k + kappa[j]
        if rho[i] is None: rho[i] = r
        elif rho[i] != r: bad += 1
print('rho conflicts', bad, 'unset', sum(1 for r in rho if r is None), flush=True)
S = sum(Y, sp.zeros(n, n))
dm = DomainMatrix.from_Matrix(S).convert_to(QQ)
t0 = time.time()
inv = dm.inv()
print('inv time', time.time()-t0, flush=True)
Yi = inv.to_Matrix()
pickle.dump({'kappa': kappa, 'rho': rho, 'Y1inv': {(i, j): v for (i, j), v in Yi.todok().items()}}, open('master_kernel2_Y1inv%d.pkl' % L, 'wb'))
