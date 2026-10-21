import sys, time, pickle
sys.path.insert(0, '.')
from jet_basis_sym import *
import master_kernel2 as m
L = int(sys.argv[1]); lam = (L, 0, 0)
zero = (Fr(0), Fr(0), Fr(0))
t0 = time.time()
nu, cond = nu_formulas(L)
print('nu', time.time() - t0, 'cond', cond, flush=True)
prs = gk.pairs_of(lam); n = len(prs)
lev = [p[0] + p[1] - 1 for p in prs]
rho = [lev[i] + lev[j] for i in range(n) for j in range(n)]
kappa = m.grading(L, None)
Y1s = []
for k in (31, 42):
    th = (Fr(2, 7), Fr(2, 7) - Fr(35, k))
    t = 1 / (th[0] - th[1])
    t1 = time.time()
    rows = Y_from_operators(lam, zero, th, nu=nu)
    print('t', t, time.time() - t1, flush=True)
    Y1 = {}
    for r, row in enumerate(rows):
        for c, v in enumerate(row):
            if v != 0:
                Y1[(r, c)] = v * t ** (kappa[c] - rho[r])
    Y1s.append(Y1)
    pickle.dump(Y1s, open('master_kernel2_Y1pts%d.pkl' % L, 'wb'))
print('Y1 agree at two t:', Y1s[0] == Y1s[1], len(Y1s[0]), len(Y1s[1]), flush=True)
