"""Span check of ndsu3lib HW-row coefficients (site 2 at its highest weight) against the PRV space M_mu of cg-slices-1470.txt.
usage: python -I span_check.py dump_hwrow.txt   (run with Paul/Python on sys.path via PYTHONPATH set below)"""
import sys, itertools, math
sys.path.insert(0, '/home/paul/git_projects/gaudin-CG/Paul/Python/gaudin_sl3')
sys.path.insert(0, '/home/paul/git_projects/gaudin-CG/Paul/Python')
import numpy as np
from sympy.physics.wigner import clebsch_gordan
from sympy import Rational
import slice_engine as se

LAM = (14, 7, 0)
rep = se.Rep.of(LAM)
sg = rep.sqrt_metric()

def parse(fn):
    out = []; cur = None
    for l in open(fn):
        f = l.split()
        if f[0] == 'IRREP':
            cur = dict(mu=tuple(map(int, f[1:4])), rho=int(f[4]), rows=[]); out.append(cur)
        else:
            cur['rows'].append((int(f[1]), int(f[2]), [float(x) for x in f[3:]]))
    return out

def cg(twoL1, twoM1, twoL2, twoM2, twoL3, twoM3):
    return float(clebsch_gordan(Rational(twoL1, 2), Rational(twoL2, 2), Rational(twoL3, 2), Rational(twoM1, 2), Rational(twoM2, 2), Rational(twoM3, 2)))

def analyse(ir, signfun=None, search=False):
    mu = ir['mu']; m = ir['rho']
    w0 = (mu[0]-LAM[0], mu[1]-LAM[1], mu[2]-LAM[2])
    W = rep.wsp[w0]
    prv = se.prv_space(LAM, mu)
    assert prv.m == m, (mu, prv.m, m)
    B = np.array([[float(x) for x in b] for b in prv.basis]).T          # w x m, columns span M_mu (aligned with W)
    Q, _ = np.linalg.qr(B)
    twoL3 = mu[0]-mu[1]
    table = {(e, L): c for e, L, c in ir['rows']}
    V = np.zeros((len(W), m)); keys = []
    for a, nu in enumerate(W):
        p = rep.pats[nu]                      # ((k,),(m12,m22),(l1,l2,l3))
        k = p[0][0]; m12, m22 = p[1]
        eps1 = 3*(m12+m22) - 2*sum(LAM); twoL1 = m12-m22; twoM1 = 2*k-(m12+m22)
        keys.append((eps1, twoL1))
        c = table.get((eps1, twoL1))
        if c is None:
            continue
        # site 2 at HW: Lambda2=7/2, M2=7/2 ; coupled top state Lambda3=L3/2, M3=Lambda3
        assert twoM1 == twoL3 - 7, (twoM1, twoL3)
        g = cg(twoL1, twoM1, 7, 7, twoL3, twoL3)
        V[a, :] = np.array(c) * g / sg[nu]
    return W, keys, V, Q

def resid(V, Q, s=None):
    if s is not None: V = V * s[:, None]
    R = V - Q @ (Q.T @ V)
    nv = np.linalg.norm(V, axis=0)
    return max(np.linalg.norm(R[:, j]) / nv[j] for j in range(V.shape[1])) if V.shape[1] else 0.0, V

if __name__ == '__main__':
    data = parse(sys.argv[1])
    worst0 = 0; worst1 = 0; nfail0 = 0; res = []
    for ir in data:
        W, keys, V, Q = analyse(ir)
        r0, _ = resid(V, Q)
        worst0 = max(worst0, r0); nfail0 += r0 > 1e-8
        res.append((ir['mu'], r0, len(W)))
    print('plain signs (+1 everywhere): max relative residual %.3e, irreps failing (>1e-8): %d / %d' % (worst0, nfail0, len(data)))
    # ---- strength of the test
    nontriv = [(ir['mu'], len(W_), ir['rho']) for ir in data for W_ in [rep.wsp[tuple(ir['mu'][i]-LAM[i] for i in range(3))]] if len(W_) > ir['rho']]
    print('irreps with dim(weight space) w > m (non-vacuous test):', len(nontriv), '; max (w-m) =', max(w-m for _, w, m in nontriv))
    miss = 0; sv_min = 1.0; sv_proj_min = 1.0
    for ir in data:
        W, keys, V, Q = analyse(ir)
        tab = {(e, L) for e, L, _ in ir['rows']}
        L3 = ir['mu'][0]-ir['mu'][1]
        # states of the slice weight space absent from the dump must violate the SU(2) triangle |L1-7/2|<=L3<=L1+7/2 (CG = 0)
        absent = set(keys) - tab
        ok_abs = all(not (abs(L-7) <= L3 <= L+7 and (L+7+L3) % 2 == 0) for _, L in absent)
        if (tab - set(keys)) or not ok_abs: miss += 1
        nabs = locals().get('nabs', 0) + len(absent)
        Vn = V / np.linalg.norm(V, axis=0)
        sv = np.linalg.svd(Vn, compute_uv=False)
        sv_min = min(sv_min, sv[-1])
        # span equality: principal angles between span(V) and M_mu (both dim m) via projector difference
        Qv, _ = np.linalg.qr(V)
        sv_proj_min = min(sv_proj_min, np.linalg.svd(Q.T @ Qv, compute_uv=False)[-1])
    print('dump keys subset of slice weight-space states, and absent states all forbidden by the SU(2) triangle rule:', miss == 0, '(violating irreps: %d; %d absent states in total)' % (miss, nabs))
    print('min singular value of the normalised m ndsu3lib vectors per irrep (independence): %.4f' % sv_min)
    print('min cos(principal angle) span(ndsu3lib) vs M_mu: %.15f  (1 = identical subspaces)' % sv_proj_min)
    # ---- negative controls
    rng = np.random.default_rng(1)
    for name in ('flip one random component sign',):
        mx = []
        for ir in data:
            W, keys, V, Q = analyse(ir)
            if len(W) <= ir['rho']: continue
            s = np.ones(len(W))
            if name.startswith('flip'): s[rng.integers(len(W))] = -1
            r, _ = resid(V, Q, s); mx.append(r)
        print('control "%s": over the %d non-vacuous irreps residual>1e-8 in %d, median %.2e' % (name, len(mx), sum(x > 1e-8 for x in mx), float(np.median(mx))))
