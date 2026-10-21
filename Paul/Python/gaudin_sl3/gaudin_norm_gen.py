"""Q2 (round 5): the Gaudin norm formula N = <Psi|Psi>_G/<x0^G|Psi>^2 = C |det J6| prod (zeta_i-zeta_j)^(p_ij(mu)) for general Lambda ([2,1,0],[3,1,0],[3,2,0]).
python -X utf8 gaudin_norm_gen.py [stages]  (stages: TAB, FIT, SIGN, REG)   Reuses nu_general / general_kernel read-only."""
import sys, os, time, pickle
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import sympy as sp
import mpmath as mp
import nu_general as NG
import general_kernel as gk

DPS = 80
mp.mp.dps = DPS
ZET = gk.ZETAS2
TH0 = gk.TH0


def cpath(name):
    return os.path.join(HERE, 'gaudin_norm_gen_%s.pkl' % name)


def cached(name, fn):
    p = cpath(name)
    if os.path.exists(p):
        with open(p, 'rb') as f:
            return pickle.load(f)
    r = fn()
    with open(p + '.tmp', 'wb') as f:
        pickle.dump(r, f)
    os.replace(p + '.tmp', p)
    return r


from math import factorial as _f
def Cpred(l1, l2):
    return Fr(_f(l2) * _f(l1 - l2) * _f(l1 + 1), (l2 + 1) * (l1 - l2 + 1)) ** 2


def mpf(x):
    return mp.mpf(x.numerator) / x.denominator if isinstance(x, Fr) else mp.mpf(x)


_REP = {}
class MiniRep:
    """rep data without W4 (chain, one-site rep, Gram metric)"""
    def __init__(self, lam):
        self.lam = tuple(lam)
        self.ch = NG.chain_of(lam)
        self.d = self.ch.d1
        self.gd = NG.gram_diag(self.ch)
        self.gdm = [mp.mpf(q.numerator) / q.denominator for q in self.gd]
        self.one = gk.one(lam)


def rep(lam):
    if lam not in _REP:
        _REP[lam] = MiniRep(lam)
    return _REP[lam]


# ------------------------------------------------------------------ Jacobian of the six no-log conditions (N1's raw normalisation)
_J = {}
def jac_funcs(lam):
    """F = (A_0,B1_0,B2_0,A_1,B1_1,B2_1) as functions of (H1,a1,b1,H2,a2,b2,e1,e2,e3,d), d = theta1-theta2, and the 6x6 Jacobian; lambdified"""
    lam = tuple(lam[:2])
    if lam in _J:
        return _J[lam]
    nl = NG.frobenius(lam)['nolog']
    H1, a1, b1, H2, a2, b2, e1, e2, e3, d = sp.symbols('H1 a1 b1 H2 a2 b2 e1 e2 e3 d')
    Fs = []
    for s in range(2):
        if s == 0:
            sub = {NG.Hs: H1, NG.As: a1, NG.Bs: b1, NG.Ht: H2, NG.At: a2, NG.Bt: b2, NG.DL: d}
        else:
            sub = {NG.Hs: H2, NG.As: a2, NG.Bs: b2, NG.Ht: H1, NG.At: a1, NG.Bt: b1, NG.DL: -d}
        sub.update({NG.SG1: e1, NG.E2: e2, NG.E3: e3})
        for f in nl:
            Fs.append(sp.expand(f.subs(sub, simultaneous=True)))
    Fm = sp.Matrix(Fs)
    D = sp.Matrix([H1, a1, b1, H2, a2, b2])
    J = Fm.jacobian(D)
    args = (H1, a1, b1, H2, a2, b2, e1, e2, e3, d)
    out = (sp.lambdify(args, Fm, 'mpmath'), sp.lambdify(args, J, 'mpmath'))
    _J[lam] = out
    return out


def laurent6(st):
    L0 = NG.laurent_of(st['lam2'], st['lam3'], 0)
    L1 = NG.laurent_of(st['lam2'], st['lam3'], 1)
    return list(L0) + list(L1)


def jac_at(lam, lau, zeta, th):
    Ff, Jf = jac_funcs(lam)
    z = [mpf(x) for x in zeta]
    e1 = sum(z); e2 = z[0] * z[1] + z[0] * z[2] + z[1] * z[2]; e3 = z[0] * z[1] * z[2]
    args = list(lau) + [e1, e2, e3, mpf(th[0] - th[1])]
    return mp.matrix(Ff(*args)), mp.matrix(Jf(*args))


# ------------------------------------------------------------------ ground truth table
def table(lam, tag, th=TH0, zeta=None):
    """per-state records at (lam, tag, th): sector, laurent (mp), N = <Psi|Psi>_G/<x0^G|Psi>^2, det J6, det of the per-site 3x3 diagonal blocks"""
    lam = tuple(lam)
    zeta = zeta or ZET[tag]
    key = '%d%d%d_%s_%s' % (lam[0], lam[1], lam[2], tag, gk.thkey(th))
    def build():
        pt = gk.get_point2(lam, tag, zeta, th, DPS)
        mp.mp.dps = DPS
        R = rep(lam)
        d = R.d
        x = [mpf(q) for q in gk.x_cov(R.one, zeta)]
        gd = R.gdm
        out = []
        for st in pt['states']:
            nrm = mp.fsum(gd[i // d] * gd[i % d] * v ** 2 for i, v in zip(st['idx'], st['v']))
            ox = mp.fsum(x[i // d] * x[i % d] * v for i, v in zip(st['idx'], st['v']))
            lau = laurent6(st)
            F, J = jac_at(lam, lau, zeta, th)
            out.append(dict(sector=tuple(st['sector']), laurent=[str(q) for q in lau], N=str(nrm / ox ** 2), x0=str(ox), norm=str(nrm), Fres=[str(q) for q in F],
                            J=[[str(J[i, j]) for j in range(6)] for i in range(6)], ev=str(st['ev'])))
        return out
    recs = cached('tab_' + key, build)
    return recs


def unpack(r):
    J = mp.matrix(6, 6)
    for i in range(6):
        for j in range(6):
            J[i, j] = mp.mpf(r['J'][i][j])
    return mp.mpf(r['N']), J


def copyset(recs):
    from collections import Counter
    return Counter(r['sector'] for r in recs)


# extra zeta points (own tags 'q2a'...) for the zeta-dependence of the constant
ZQ = {'q2a': (Fr(5, 2), Fr(1, 3), Fr(-3, 4)), 'q2b': (Fr(-2, 5), Fr(7, 3), Fr(4, 9)), 'q2c': (Fr(9, 4), Fr(-5, 6), Fr(2, 7)),
      'q2d': (Fr(1, 2), Fr(-7, 5), Fr(11, 6)), 'q2e': (Fr(-8, 3), Fr(3, 5), Fr(1, 4)), 'q2f': (Fr(6, 5), Fr(2, 3), Fr(-9, 7)),
      'q2g': (Fr(-1, 3), Fr(8, 7), Fr(5, 2)), 'q2h': (Fr(13, 6), Fr(-4, 7), Fr(3, 8))}
ZALL = dict(ZET); ZALL.update(ZQ)
ZET.update(ZQ)

if __name__ == '__main__' and 'BUILD' in sys.argv:
    lam = tuple(int(c) for c in sys.argv[2])
    for tag in sys.argv[3:]:
        t = time.time()
        table(lam, tag)
        print(lam, tag, 'table built %.0fs' % (time.time() - t), flush=True)
