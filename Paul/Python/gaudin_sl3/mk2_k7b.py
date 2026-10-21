import sys, pickle, itertools
import sympy as sp
from sympy import Rational as R
from mk2_k7 import *
lam = int(sys.argv[1]); a = R(1, 3); ap = R(1, 5)
Qn = build(lam, a, ap)
E1, E2, Om = casimir(lam)
n = Om.shape[0]
vs = {k: vec(d, lam) for k, d in Qn.items()}
# W_Q = span of Omega^k Q_n
cols = []
for k in sorted(vs):
    w = vs[k]
    for j in range(lam + 2):
        cols.append(w); w = Om * w
WQ = sp.Matrix.hstack(*cols)
Wb = WQ.columnspace()
print('dim W_Q', len(Wb))
B = sp.Matrix.hstack(*Wb)
# operator family
fam = {}
fam['1'] = sp.eye(n)
for i in range(3):
    for j in range(3):
        fam['E1_%d%d' % (i, j)] = E1[(i, j)]
        fam['E2_%d%d' % (i, j)] = E2[(i, j)]
for i, j, k, l in itertools.product(range(3), repeat=4):
    fam['E1_%d%d*E2_%d%d' % (i, j, k, l)] = E1[(i, j)] * E2[(k, l)]
names = list(fam)
# condition: (sum c_f F) B has rank unchanged when appended: i.e. in column space of B. Use complement projector: N^T (F B)=0 with N = nullspace basis of B^T
N = sp.Matrix.hstack(*(B.T.nullspace()))
print('complement dim', N.shape[1])
rows = sp.Matrix.hstack(*[ (N.T * fam[nm] * B).reshape(N.shape[1] * B.shape[1], 1) for nm in names])
ns = rows.nullspace()
print('operators in family preserving W_Q:', len(ns))
sols = []
for v in ns:
    d = {names[k]: v[k] for k in range(len(names)) if v[k] != 0}
    sols.append(d)
    pass
pickle.dump(sols, open('master_kernel2_k7_stab%d.pkl' % lam, 'wb'))
# --- modulo operator relations
def flat(M): return M.reshape(M.shape[0]*M.shape[1],1)
allm = sp.Matrix.hstack(*[flat(fam[nm]) for nm in names])
print('family dim', len(names), 'rank as operators on V(x)V', allm.rank())
stab = [sum((v[k]*fam[names[k]] for k in range(len(names)) if v[k]!=0), sp.zeros(n,n)) for v in ns]
Sm = sp.Matrix.hstack(*[flat(s) for s in stab])
print('stabiliser operators (independent as operators):', Sm.rank())
lie_names=[nm for nm in names if '*' not in nm]
lie_idx=[names.index(nm) for nm in lie_names]
lie = [v for v in ns]
# stabiliser within Lie part only
rowsL = sp.Matrix.hstack(*[ (N.T * fam[nm] * B).reshape(N.shape[1]*B.shape[1],1) for nm in lie_names])
nsL = rowsL.nullspace()
LM = sp.Matrix.hstack(*[flat(sum((v[k]*fam[lie_names[k]] for k in range(len(lie_names)) if v[k]!=0), sp.zeros(n,n))) for v in nsL])
print('Lie part (single-site + 1) stabiliser rank (as operators):', LM.rank(), 'of', sp.Matrix.hstack(*[flat(fam[nm]) for nm in lie_names]).rank())
# mod out by span{1,Omega, Lie stabiliser, and products Om*...}
base = [sp.eye(n), Om] + [sum((v[k]*fam[lie_names[k]] for k in range(len(lie_names)) if v[k]!=0), sp.zeros(n,n)) for v in nsL]
Bm = sp.Matrix.hstack(*[flat(x) for x in base])
print('rank of {1,Omega,Lie stab}', Bm.rank(), ' rank with all stabiliser', sp.Matrix.hstack(Bm, Sm).rank())
pickle.dump((names, ns, nsL), open('master_kernel2_k7_stab%d.pkl' % lam, 'wb'))
