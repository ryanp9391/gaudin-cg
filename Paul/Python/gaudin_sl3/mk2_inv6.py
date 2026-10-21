import sys, time, pickle
import master_kernel2 as m, sympy as sp
import general_kernel as G
from sympy.polys.matrices import DomainMatrix
from sympy import QQ
L = 6
Y1 = pickle.load(open('master_kernel2_Y1pts6.pkl', 'rb'))[0]
n = 784
prs = G.pairs_of((L, 0, 0)); npr = len(prs)
lev = [p[0] + p[1] - 1 for p in prs]
rho = [lev[i] + lev[j] for i in range(npr) for j in range(npr)]
kappa = m.grading(L, None)
dm = DomainMatrix({i: {} for i in range(n)}, (n, n), QQ)
rows = {}
for (r, c), v in Y1.items():
    rows.setdefault(r, {})[c] = QQ(v.numerator, v.denominator)
dm = DomainMatrix(rows, (n, n), QQ)
t0 = time.time()
inv = dm.inv()
print('inv time', time.time() - t0, flush=True)
Yi = inv.to_Matrix()
pickle.dump({'kappa': kappa, 'rho': rho, 'Y1inv': {(i, j): v for (i, j), v in Yi.todok().items()}}, open('master_kernel2_Y1inv6.pkl', 'wb'))
T = m.theta_T(6)
pickle.dump(T, open('master_kernel2_rep_600.pkl', 'wb'))
print('rep_600 entries', len(T), sum(len(v) for v in T.values()), 'max t-power', max(k[0] for v in T.values() for k in v))
