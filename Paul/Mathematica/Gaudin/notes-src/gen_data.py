# -*- coding: utf-8 -*-
"""Generate appendix tables (markup fragments) from the exact data files / RESULTS files."""
import re, os, math, sys
from fractions import Fraction as Fr
from math import comb

ROOT = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get('GAUDIN_REPO') or os.path.abspath(os.path.join(ROOT, '..', '..', '..', '..'))  # notes-src/ lives in Paul/Mathematica/Gaudin/
REF = os.path.join(REPO, 'Paul', 'Mathematica', 'Gaudin', 'ref')
OUT = os.path.join(ROOT, 'gen')
os.makedirs(OUT, exist_ok=True)


def parse_nested(s):
    s = s.strip()
    pos = [0]

    def rec():
        if s[pos[0]] == '{':
            pos[0] += 1
            lst = []
            while s[pos[0]] != '}':
                lst.append(rec())
                if s[pos[0]] == ',':
                    pos[0] += 1
            pos[0] += 1
            return lst
        m = re.match(r'-?\d+(?:/\d+)?', s[pos[0]:])
        tok = m.group(0)
        pos[0] += len(tok)
        return Fr(tok)
    return rec()


def ftex(x, paren=False):
    x = Fr(x)
    if x.denominator == 1:
        return str(x.numerator)
    sgn = '-' if x < 0 else ''
    return sgn + '\\frac{%d}{%d}' % (abs(x.numerator), x.denominator)


def pair_tex(p, q, D0):
    p = Fr(p); q = Fr(q)
    if D0 == 0 or q == 0:
        return ftex(p)
    sq = '\\sqrt{%d}' % D0
    if p == 0:
        if q == 1: return sq
        if q == -1: return '-' + sq
        return ftex(q) + sq
    qs = ''
    if q == 1: qs = '+' + sq
    elif q == -1: qs = '-' + sq
    elif q > 0: qs = '+' + (ftex(q) if q.denominator != 1 else str(q)) + sq
    else: qs = '-' + (ftex(-q) if q.denominator != 1 else str(-q)) + sq
    if q.denominator != 1 and q > 0:
        qs = '+' + ftex(q) + sq
    return ftex(p) + qs


def pair_val(p, q, D0):
    return float(p) + (float(q) * math.sqrt(D0) if D0 else 0.0)


def cell(s):
    return '$' + s + '$'


REPS = ['110', '210', '220', '310', '320', '410']

branches = []   # dicts
for rep in REPS:
    lam = tuple(int(c) for c in rep)
    w1 = None
    for ln in open(os.path.join(REF, 'untwisted-exact-%s.txt' % rep), encoding='utf8'):
        ln = ln.rstrip('\n')
        if ln.startswith('WEIGHTS1'):
            w1 = parse_nested(ln.split('\t')[1])
        if ln.startswith('SECTOR'):
            f = ln.split('\t')
            cur_sector = (parse_nested(f[1]), int(f[2]))
        if ln.startswith('BRANCH'):
            f = ln.split('\t')
            mu = tuple(int(x) for x in parse_nested(f[1]))
            cp = int(f[2]); D0 = int(f[3])
            vec = parse_nested(f[4]); gn = parse_nested(f[5]); cg = parse_nested(f[6])
            eta = parse_nested(f[7]); al = parse_nested(f[8]); be = parse_nested(f[9])
            branches.append(dict(lam=lam, mu=mu, copy=cp, D0=D0, gn=gn, cg=cg, eta=eta, al=al, be=be, w1=w1, mult=cur_sector[1]))

# ------------------------------------------------------------------ D1: CG slices
rows = ['| $\\Lambda$ | $\\mu$ | copy | $\\nu$ (site-1 weight) | unnormalised $\\langle e_\\nu\\otimes e_1|\\Psi_0\\rangle$ | $\\|\\Psi_0\\|^2_G$ | $C_\\nu$ |']
lines = []
for b in branches:
    D0 = b['D0']
    g0, g1 = b['gn']
    gn_tex = pair_tex(g0, g1, D0)
    gval = pair_val(g0, g1, D0)
    for k, (nu, cA, cB) in enumerate(b['cg']):
        nu = int(nu)
        wt = b['w1'][nu - 1]
        val = pair_val(cA, cB, D0) / math.sqrt(gval)
        first = (k == 0)
        lamS = '(%d,%d,%d)' % b['lam'] if first else ''
        muS = '(%d,%d,%d)' % b['mu'] if first else ''
        cpS = (str(b['copy']) if b['mult'] == 2 else '') if first else ''
        lines.append('| %s | %s | %s | %d, (%d,%d,%d) | %s | %s | %.10f |' % (
            cell(lamS) if lamS else '', cell(muS) if muS else '', cpS, nu, wt[0], wt[1], wt[2],
            cell(pair_tex(cA, cB, D0)), cell(gn_tex) if first else '', val))
open(os.path.join(OUT, 'tab-cg.src'), 'w', encoding='utf8').write(
    '@table tab:cg | Exact CG slices $C_\\nu=\\langle e_\\nu\\otimes e_1|\\Psi_0\\rangle/\\|\\Psi_0\\|_G$ of the highest-weight kets of $V_\\Lambda\\otimes V_\\Lambda$, six representations (45 HW branches); $\\nu$ is the site-1 Gelfand--Tsetlin index in notebook order (index 1 = highest weight) with its weight; the two copies in a multiplicity-2 sector are the eigenvectors of $\\theta_{12}\\,a_1$, copy 1 having the smaller eigenvalue; $C_\\nu=(c_A+c_B\\sqrt{D_0})/\\sqrt{g_A+g_B\\sqrt{D_0}}$. Source: `ref/untwisted-exact-<Lambda>.txt`, records `BRANCH`.\n@size scriptsize\n@cols llcccc>{\\raggedleft\\arraybackslash}p{2.6cm}\n' +
    rows[0].replace('| $\\Lambda$ | $\\mu$ | copy | $\\nu$ (site-1 weight) | unnormalised $\\langle e_\\nu\\otimes e_1|\\Psi_0\\rangle$ | $\\|\\Psi_0\\|^2_G$ | $C_\\nu$ |',
                    '| $\\Lambda$ | $\\mu$ | copy | $\\nu$, weight | $c_A+c_B\\sqrt{D_0}$ | $g_A+g_B\\sqrt{D_0}$ | $C_\\nu$ |') + '\n|---|---|---|---|---|---|---|\n' + '\n'.join(lines) + '\n@endtable\n')

# ------------------------------------------------------------------ E1: sectors / Laurent data / oper data
def threeF2(lam, mu):
    l1, l2, l3 = lam
    exps_th = [0, l2 + 1, l1 + 2]
    exps_inf = [mu[2], mu[1] + 1, mu[0] + 2]
    def consec(e):
        e = sorted(e)
        return any(e[i + 1] - e[i] == 1 for i in range(2))
    return consec(exps_th) or consec(exps_inf)

rows = []
n3f2 = 0
nbr = 0
n_dir_exc = 0
for b in branches:
    lam, mu = b['lam'], b['mu']
    D0 = b['D0']
    nbr += 1
    eta = b['eta'][0][0]
    al = pair_tex(b['al'][0][0], b['al'][0][1], D0)
    be = pair_tex(b['be'][0][0], b['be'][0][1], D0)
    dinf = '(%d,%d,%d)' % (mu[0] + 2, mu[1] + 1, mu[2])
    dd = (2 * lam[0] + 2 - mu[2], 2 * lam[0] + 1 - mu[1], 2 * lam[0] - mu[0])
    if (lam, mu) == ((1, 1, 0), (2, 2, 0)):
        dd = (0, 2, 4)
    if (lam, mu) == ((2, 2, 0), (4, 4, 0)):
        dd = (0, 2, 6)
    f = threeF2(lam, mu)
    n3f2 += f
    rows.append('| %s | %s | %s | %s | %s | %s | %s | %s | %s |' % (
        cell('(%d,%d,%d)' % lam), cell('(%d,%d,%d)' % mu), (str(b['copy']) if b['mult'] == 2 else ''),
        cell(ftex(eta)), cell(al), cell(be), cell(dinf), cell('(%d,%d,%d)' % dd), 'yes' if f else 'no'))
open(os.path.join(OUT, 'tab-untw.src'), 'w', encoding='utf8').write(
    '@table tab:untw | Untwisted Laurent data and oper data on all 45 highest-weight branches: $H_1=-H_2=\\eta/\\theta_{12}$, $a_1=-a_2=\\alpha/\\theta_{12}$, $b_1=-b_2=\\beta/\\theta_{12}^2$; exponents of the dual oper at infinity $(\\mu_1+2,\\mu_2+1,\\mu_3)$ (the three polynomial solutions have exactly these degrees); degrees of the three polynomial solutions of the untwisted direct oper (see text for the two exceptions); last column: whether the oper is a ${}_3F_2$ equation (consecutive exponent pair at a singular point). Source: `ref/untwisted-exact-<Lambda>.txt`.\n@size scriptsize\n@cols lllrccccc\n| $\\Lambda$ | $\\mu$ | copy | $\\eta$ | $\\alpha$ | $\\beta$ | at $\\infty$ | direct deg. | ${}_3F_2$ |\n|---|---|---|---|---|---|---|---|---|\n' + '\n'.join(rows) + '\n@endtable\n')
print('branches', nbr, '3F2', n3f2)

# ------------------------------------------------------------------ D2: adjoint series
ser = {}
for ln in open(os.path.join(REF, 'untwisted-exact-210.txt'), encoding='utf8'):
    if ln.startswith('E4_SERIES'):
        f = ln.rstrip('\n').split('\t')
        mu = tuple(int(x) for x in parse_nested(f[1])); cp = int(f[2]); D0 = int(f[3]); phi = parse_nested(f[4])
        if phi[0] != Fr(1, 3):
            continue
        H = parse_nested(f[6]); A = parse_nested(f[7]); B = parse_nested(f[8])
        ser[(mu, cp)] = (D0, H, A, B)
names = {(4, 2, 0): 'top', (4, 1, 1): '(4,1,1)', (3, 3, 0): '(3,3,0)', (3, 2, 1): 'octet', (2, 2, 2): 'singlet'}
rows_ex = []
rows_num = []
for (mu, cp), (D0, H, A, B) in ser.items():
    lab = '(%d,%d,%d)' % mu + ((' copy %d' % cp) if mu == (3, 2, 1) else '')
    for nm, X in (('H_1', H[0]), ('a_1', A[0]), ('b_1', B[0])):
        ex = []
        for n in range(3):
            p, q = X[n]
            ex.append(cell(pair_tex(p, q, D0)) if (p != 0 or q != 0) else '0')
        rows_ex.append('| %s | %s | %s |' % (cell(lab), cell(nm), ' | '.join(ex)))
        nums = []
        for n in range(8):
            p, q = X[n]
            v = pair_val(p, q, D0)
            nums.append('%.6g' % v if v != 0 else '0')
        rows_num.append('| %s | %s | %s |' % (cell(lab), cell(nm), ' | '.join(nums)))
open(os.path.join(OUT, 'tab-series-exact.src'), 'w', encoding='utf8').write(
    '@table tab:series-ex | Exact $\\epsilon$-series coefficients of the Laurent data of the six adjoint branches at $\\phi=(1/3,-2/7,5/11)$, $\\theta_{12}=31/35$ (so $H_1^{(0)}=4/\\theta_{12}$ etc.), orders $n=0,1,2$ ($H_s=\\sum_n\\epsilon^nH_s^{(n)}$). Source: `ref/untwisted-exact-210.txt`, record `E4_SERIES` (agent W4; exact Rayleigh--Schr\\"odinger).\n@size scriptsize\n@cols llccc\n| branch | series | $n=0$ | $n=1$ | $n=2$ |\n|---|---|---|---|---|\n' + '\n'.join(rows_ex) + '\n@endtable\n')
open(os.path.join(OUT, 'tab-series-num.src'), 'w', encoding='utf8').write(
    '@table tab:series-num | The same series, numerically (6 significant digits), orders $n=0,\\dots,7$; exact values to order 7 are in the data file. Odd orders vanish identically on the algebraic branches $(4,1,1)$ and $(3,3,0)$ beyond $n=1$ and on the singlet for $H_1$ (see text).\n@size tiny\n@cols llrrrrrrrr\n| branch | series | $n=0$ | $n=1$ | $n=2$ | $n=3$ | $n=4$ | $n=5$ | $n=6$ | $n=7$ |\n|---|---|---|---|---|---|---|---|---|---|\n' + '\n'.join(rows_num) + '\n@endtable\n')

# ------------------------------------------------------------------ C: one-site kernel c_kl
def ckl(l1, l2, k, l):
    n = l1 + 1
    sgn = (-1) ** (k + l - l2 - 1)
    return Fr(sgn * comb(l1 - l2, l - l2 - 1) * comb(n, k) * (l - k), (l2 + 1) * comb(l, l2 + 1))

rows = []
for lam in [(1, 0), (1, 1), (2, 0), (2, 1), (2, 2), (3, 0), (3, 1), (3, 2), (3, 3), (4, 0), (4, 1), (4, 2), (4, 3), (4, 4)]:
    l1, l2 = lam
    prs = [(k, l) for k in range(l1 + 2) for l in range(k + 1, l1 + 2) if l >= l2 + 1]
    ent = ', '.join('(%d%d): %s' % (k, l, '$' + ftex(ckl(l1, l2, k, l)) + '$') for k, l in prs)
    rows.append('| %s | %d | %s |' % (cell('[%d,%d,0]' % (l1, l2)), len(prs), ent))
open(os.path.join(OUT, 'tab-onesite.src'), 'w', encoding='utf8').write(
    '@table tab:onesite | The one-site kernel $\\theta_{kl}(a)=c_{kl}\\,a^{k+l-\\lambda_2-1}$: coefficients $c_{kl}$ from the closed form (eq.~\\ref{eq:ckl}) for all two-row (and one-row) representations up to $\\lambda_1=4$. Pairs $(kl)$ in lexicographic order, $0\\le k<l\\le\\lambda_1+1$, $l\\ge\\lambda_2+1$. The entries for $\\Lambda\\le[4,2,0]$ coincide with the exact one-site fits of G1; the rest were predicted by the formula and then confirmed (17 representations in all).\n@size scriptsize\n@cols lcp{11cm}\n| $\\Lambda$ | # pairs | $c_{kl}$ |\n|---|---|---|\n' + '\n'.join(rows) + '\n@endtable\n')

# ------------------------------------------------------------------ B: adjoint kernel G (D2) and twisted G~
txt = open(os.path.join(REPO, 'Paul', 'Python', 'gaudin_sl3', 'RESULTS-twovac-limit.md'), encoding='utf8').read()


def conv(expr):
    e = expr.strip()
    e = e.replace('w1p', 'W1P').replace('w1', 'W1').replace('thi', 'T')
    e = re.sub(r'\bd\b', 'D', e)
    e = e.replace('*', ' ')
    e = re.sub(r'\^(\d+)', r'^{\1}', e)
    e = e.replace('W1P', "\\omega_1'").replace('W1', '\\omega_1').replace('T', 't').replace('D', '\\delta')
    e = re.sub(r"\\omega_1'\s+", "\\\\omega_1' ", e)
    # tidy: join digits followed by space and letter
    return e


Gs = {}
for m in re.finditer(r'^G\[(\d),(\d)\] = (.*)$', txt, flags=re.M):
    Gs[(int(m.group(1)), int(m.group(2)))] = m.group(3)
print('G entries', len(Gs))
labs = ['P^3p_{02}', "P^{3\\prime}p_{02}", 'P^3p_{03}', "P^{3\\prime}p_{03}", 'P^3p_{23}']
# D2 order: N = (P3 p02, P3' p02, P3 p03, P3' p03, P3 p23)
lines = []
for i in range(1, 6):
    for j in range(1, 6):
        lines.append('| %d,%d | %s |' % (i, j, cell(conv(Gs[(i, j)]))))
open(os.path.join(OUT, 'tab-G210.src'), 'w', encoding='utf8').write(
    '@table tab:G210 | The adjoint kernel $G_{ij}(\\omega_1,\\omega_1\',t)$, $t=1/\\theta_{12}$, in the basis $N=(P^3p_{02},\\,P^{3\\prime}p_{02},\\,P^3p_{03},\\,P^{3\\prime}p_{03},\\,P^3p_{23})$ (weights $2,3,3,4,5$); $G=4K$. Source: `RESULTS-twovac-limit.md` section L2(d), parsed from the agent\'s printout (D2); equal to the independent fit of D1 in 25/25 entries.\n@size scriptsize\n@cols cp{13.5cm}\n| $(i,j)$ | $G_{ij}$ |\n|---|---|\n' + '\n'.join(lines) + '\n@endtable\n')

Gt = {}
block = txt.split('**Twisted-jet form**')[1].split('(d = omega1')[0]
for m in re.finditer(r'G~\[(\d),(\d)\] =\s*(.*?)(?=\s{2,}G~\[|\s*$)', block, flags=re.M):
    pass
for m in re.finditer(r'G~\[(\d),(\d)\]\s*=\s*([^G]*?)(?=\s{2,}G~\[|$)', block, flags=re.M):
    Gt[(int(m.group(1)), int(m.group(2)))] = m.group(3).strip()
print('Gt entries', len(Gt))


def convt(e):
    e = e.strip().replace('*', ' ')
    e = re.sub(r'\bthi\b', 't', e)
    e = re.sub(r'\bd\b', '\\\\delta', e)
    e = re.sub(r'\^(\d+)', r'^{\1}', e)
    return e


lines = []
for i in range(1, 6):
    for j in range(1, 6):
        lines.append('| %d,%d | %s |' % (i, j, cell(convt(Gt[(i, j)]))))
open(os.path.join(OUT, 'tab-Gtw.src'), 'w', encoding='utf8').write(
    '@table tab:Gtw | The adjoint kernel in twisted jets $\\tilde N$ built from $e^{-\\omega_su}P^j$ ($\\omega_s=\\omega_1$ on site 1, $\\omega_1\'$ on site 2): it depends only on $\\delta=\\omega_1\'-\\omega_1$ and $t=1/\\theta_{12}$. Same basis order as Table~\\ref{tab:G210}. Source: `RESULTS-twovac-limit.md`.\n@size scriptsize\n@cols cp{13.5cm}\n| $(i,j)$ | $\\tilde G_{ij}$ |\n|---|---|\n' + '\n'.join(lines) + '\n@endtable\n')

# ------------------------------------------------------------------ B: [3,1,0] K at base point (and [3,2,0] sub-block)
import pickle


def kernel_K(lam, w1, w1p, th12, w0=Fr(5, 2)):
    """exact K = Theta*(w1 w1')^(2 l1 - l2) from G1's stored exact polynomial laws (general_kernel_rep_<lam>.pkl);
    the fit stored C[q][n] = coefficient of a'^q t^n at fixed a = 1/w0, entry (i,j); homogeneity p + q - n = e_i + e_j restores the a-power."""
    l1, l2, _ = lam
    rep = pickle.load(open(os.path.join(REPO, 'Paul', 'Python', 'gaudin_sl3', 'general_kernel_rep_%d%d%d.pkl' % lam), 'rb'))
    prs = [(k, l) for k in range(l1 + 2) for l in range(k + 1, l1 + 2) if l >= l2 + 1]
    a0 = 1 / Fr(w0)
    a, ap, t = 1 / Fr(w1), 1 / Fr(w1p), 1 / Fr(th12)
    n = len(prs)
    M = [[Fr(0)] * n for _ in range(n)]
    for (i, j), C in rep['res'].items():
        ei, ej = prs[i][0] + prs[i][1] - l2 - 1, prs[j][0] + prs[j][1] - l2 - 1
        s = Fr(0)
        for q, row in enumerate(C):
            for nn, c in enumerate(row):
                if c is None or c == 0:
                    continue
                p = ei + ej + nn - q
                if p < 0:
                    continue
                s += (c / a0 ** p) * t ** nn * a ** p * ap ** q
        M[i][j] = s * (Fr(w1) * Fr(w1p)) ** (2 * l1 - l2)
    return prs, M


prs, K = kernel_K((3, 1, 0), Fr(5, 2), Fr(7, 4), Fr(31, 35))
assert K[8][8] == Fr(1, 9)
lines = []
for i, pi in enumerate(prs):
    ent = ', '.join('$' + ftex(K[i][j]) + '$' for j in range(9))
    lines.append('- Row $(%d%d)$: %s.' % (pi[0], pi[1], ent))
open(os.path.join(OUT, 'K310.src'), 'w', encoding='utf8').write('\n'.join(lines) + '\n')
print('done')
