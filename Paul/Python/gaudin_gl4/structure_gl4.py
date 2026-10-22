"""structure_gl4.py -- WP-G4: structure of the slice-reduced gl(4) missing-label operators V_red,3, V_red,4a, V_red,4b on the three-parameter GT lattice of a weight space
(the gl(4) analogue of the Jacobi chain of gaudin_sl3/jacobi_cg.py).   Brief: SPEC-G4-structure-gl4.md.   Report: RESULTS-G4-structure-gl4.md.

Setting.  Site 1 = V_lam in the GT basis (rows (m11), (m12,m22), (m13,m23,m33), lam).  A weight space W = V_lam[nu], nu = mu - lam', fixes the row sums and leaves THREE free lattice
coordinates (m12, m13, m23) (then m11 = nu_1, m22 = (nu_1+nu_2) - m12, m33 = (nu_1+nu_2+nu_3) - m13 - m23): the region R_nu of the GT lattice, w = |R_nu|.  The reduced operators are
weight-zero polynomials in the site-1 generators with coefficients in Q[lam'] (G3: universal_reduction); M_mu = intersection of ker E_i^{k_i} in W, k_i = lam'_i - lam'_{i+1} + 1.

Stages (CLI: python -X utf8 structure_gl4.py <stage> [args];  stage in  stencil | symbolic | symq | univ | weyl | rule | chain | poly | cycW | all;
        `run stage1,stage2,..` writes one part file per stage, `merge` assembles structure_gl4_tables.txt):
  stencil   Q1: lattice displacements (structural from the universal words, confirmed on data), tables of stencil sizes
  symbolic  Q2: closed-form matrix elements as rational functions of the GT labels and lam' (generic GT formulas through the universal words, term sums simplified by exact line reconstruction),
            verification on held-out reps;  symq: the explicit non-linear polynomials Q of the row-3 moves (exact interpolation);  univ: the ONE product formula of the 12 factorising classes
            of V3, V4a, V4b (M^V4a = (y_k + 3) M^V3, M^V4b = -(y_a - y_i + 2 y_k - q3 + 5) M^V3) and its verification
  weyl      Q3 + Q6: the Weyl/BGG characterisation of the physical spectrum: sum_w sgn(w) mult_tau(V_red(w.lam') on V_lam[mu - w.lam']) = [tau in Spec(M_mu)];
            one-step version, truncation depth, signed lattice count, gl(3) sanity check
  rule      Q4: rank of the physical triples in the orderings by V3, V4a, V4b and G2's generic combination; candidate rules; held-out products
  chain     Q5(a): sparsity graph of V_red,3, bandwidth, tridiagonality for m = w
  poly      Q5(b,c): polynomiality of the slice components in the triple, cyclic coordinate vectors, graded Krylov degree profile;  cycW: cyclic lattice points on the whole weight space W
Every table printed is also appended to structure_gl4_tables.txt (see `tab()`).
Imports G1 (glN_reps, slice_engine_gl4), G2 (gaudin_gl4_ops) and G3 (slice_ops_gl4) read-only.
"""
import os, sys, time, itertools, math, random, re
from fractions import Fraction as Fr
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import mpmath as mp
import sympy as sp
import glN_reps as gr
import slice_engine_gl4 as se
import slice_ops_gl4 as so
import gaudin_gl4_ops as go

TABLES_PATH = os.path.join(HERE, 'structure_gl4_tables.txt')
DATE = '2026-10-22'
_TABLE_BUF = []


def tab(s=''):
    """print a line AND record it for structure_gl4_tables.txt"""
    print(s, flush=True)
    _TABLE_BUF.append(s)


def hdr(s):
    tab('\n' + '=' * 150 + '\n' + s + '\n' + '=' * 150)


def flush_tables(path=TABLES_PATH, mode='w', title=None):
    with open(path, mode, encoding='utf-8') as f:
        if title:
            f.write(title + '\n')
        f.write('\n'.join(_TABLE_BUF) + '\n')
    _TABLE_BUF.clear()


def lam_str(lam):
    return '[' + ','.join(str(x) for x in lam) + ']'


def prod_str(lam, lamp):
    return lam_str(lam) + ('^2' if tuple(lam) == tuple(lamp) else 'x' + lam_str(lamp))


def mu_str(mu):
    return '(' + ','.join(str(x) for x in mu) + ')'


# the family of operators per N: (names, generic combination)
FAM = {4: dict(names=so.TRIPLE, combo=so.GENERIC), 3: dict(names=('V',), combo=(Fr(1),))}


# ======================================================================================== 1. Weyl group, dot action
def weyl_elements(N):
    """[(p, sign, length)] for all permutations p of range(N); w acts on x = lam'+rho by (w x)_i = x_{p(i)}"""
    out = []
    for p in itertools.permutations(range(N)):
        inv = sum(1 for i in range(N) for j in range(i + 1, N) if p[i] > p[j])
        out.append((p, -1 if inv % 2 else 1, inv))
    out.sort(key=lambda t: (t[2], t[0]))
    return out


WEYL = {N: weyl_elements(N) for N in (3, 4)}


def dot_action(p, lamp):
    """w . lam' = w(lam' + rho) - rho,  rho = (N-1, ..., 0)"""
    N = len(lamp)
    x = [lamp[i] + (N - 1 - i) for i in range(N)]
    return tuple(x[p[i]] - (N - 1 - i) for i in range(N))


def simple_refl(N, i):
    """permutation p of the simple reflection s_i (i = 1..N-1): swap the entries i, i+1"""
    p = list(range(N))
    p[i - 1], p[i] = p[i], p[i - 1]
    return tuple(p)


def p_name(p):
    """short name of a permutation: its one-line notation in 1-based form"""
    return ''.join(str(x + 1) for x in p)


# ======================================================================================== 2. joint spectrum on a whole weight space (complex-safe, defect-safe)
def _mpf_fr(x):
    return mp.mpf(x.numerator) / x.denominator


def jspec(mats, combo, dps=50, tries=3):
    """joint spectrum of the commuting exact w x w matrices `mats` (rows of Fractions; column j = image of basis vector j) by an mpmath eigen-decomposition of the generic
    combination sum combo_t A_t; every A_t is conjugated to the eigenbasis (U^-1 A_t U), the largest off-diagonal entry measured (accept < 10^-(dps-10) * max(1,|A|), else retry with
    +30 digits, `tries` times).  Entries are mpf when real (imag < 10^-(dps-10)) else mpc.
    Returns dict(trip = list of tuples (joint eigenvalues, one per eigenvector), U = mp.matrix (right eigenvectors as columns, scaled to largest component 1) or None for w = 1,
    ok, offmax, gap = smallest max-norm distance between two joint triples (inf for w = 1), cplx = largest imaginary part, dps = precision used)."""
    w = len(mats[0])
    r = len(mats)
    if w == 1:
        old = mp.mp.dps
        mp.mp.dps = dps
        try:
            t = tuple(_mpf_fr(A[0][0]) for A in mats)
        finally:
            mp.mp.dps = old
        return dict(trip=[t], U=None, ok=True, offmax=mp.mpf(0), gap=mp.inf, cplx=mp.mpf(0), dps=dps, w=1)
    scale = max([1.0] + [abs(float(x)) for A in mats for row in A for x in row])
    res = None
    for att in range(tries):
        d = dps + 30 * att
        old = mp.mp.dps
        mp.mp.dps = d
        try:
            M = [[sum((Fr(c) * A[i][j] for c, A in zip(combo, mats)), Fr(0)) for j in range(w)] for i in range(w)]
            E, U = mp.eig(go.fr_mp(M))
            for k in range(w):
                big = max(range(w), key=lambda i: abs(U[i, k]))
                ref = U[big, k]
                for i in range(w):
                    U[i, k] = U[i, k] / ref
            Uinv = mp.inverse(U)
            trip = [[None] * r for _ in range(w)]
            offmax = mp.mpf(0)
            for t, A in enumerate(mats):
                D = Uinv * go.fr_mp(A) * U
                for i in range(w):
                    trip[i][t] = D[i, i]
                    for j in range(w):
                        if i != j:
                            offmax = max(offmax, abs(D[i, j]))
            thr = mp.mpf(10) ** (-(dps - 10)) * scale
            cplx = max(abs(mp.im(x)) for row in trip for x in row)
            ok = offmax < thr
            tr = [tuple((mp.re(x) if abs(mp.im(x)) < thr else x) for x in row) for row in trip]
            gap = min((max(abs(a - b) for a, b in zip(tr[i], tr[j])) for i in range(w) for j in range(i + 1, w)), default=mp.inf)
            res = dict(trip=tr, U=U, ok=ok, offmax=offmax, gap=gap, cplx=cplx, dps=d, w=w)
        finally:
            mp.mp.dps = old
        if ok:
            break
    return res


def group_triples(items, tol):
    """items = [(triple, tag)]; groups = [dict(trip = representative, members = [tags], dev = max deviation inside)] (triples closer than tol, max-norm, are the same)"""
    groups = []
    for trip, tag in items:
        for g in groups:
            if max(abs(a - b) for a, b in zip(trip, g['trip'])) < tol:
                g['members'].append(tag)
                g['dev'] = max(g['dev'], max(abs(a - b) for a, b in zip(trip, g['trip'])))
                break
        else:
            groups.append(dict(trip=trip, members=[tag], dev=mp.mpf(0)))
    return groups


def min_group_sep(groups):
    best = mp.inf
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            best = min(best, max(abs(a - b) for a, b in zip(groups[i]['trip'], groups[j]['trip'])))
    return best


# ======================================================================================== 3. sectors and the physical spectrum
def dominant(mu):
    return all(mu[i] >= mu[i + 1] for i in range(len(mu) - 1))


def sector_list(lam, lamp, kind='rs'):
    """kind 'rs': the irreps of V_lam (x) V_lam' (Racah-Speiser, m > 0), descending;  kind 'all': every dominant mu with |mu| = |lam| + |lam'| for which SOME Weyl translate
    V_lam[mu - w.lam'] is non-empty (includes the RS irreps, the multiplicity-0 candidates, and the formal sectors with mu_N < 0).  Returns [(mu, m)]"""
    N = len(lam)
    rs = gr.racah_speiser(lam, lamp)
    if kind == 'rs':
        return sorted(rs.items(), reverse=True)
    R = gr.GLNRep.of(lam)
    mus = set()
    for p, sg, ln in WEYL[N]:
        wd = dot_action(p, lamp)
        for nu in R.wsp:
            mu = tuple(nu[i] + wd[i] for i in range(N))
            if dominant(mu):
                mus.add(mu)
    return sorted(((mu, rs.get(mu, 0)) for mu in mus), reverse=True)


def physical_spec(lam, lamp, mu, fam, dps=50):
    """(prv, joint spectrum on M_mu): the m joint triples of the reduced operators restricted to M_mu (empty if m = 0)"""
    prv = se.prv_space(lam, lamp, mu)
    if prv.m == 0:
        return prv, dict(trip=[], U=None, ok=True, offmax=mp.mpf(0), gap=mp.inf, cplx=mp.mpf(0), dps=dps, w=0)
    R = gr.GLNRep.of(lam)
    mats = [so.vred_matrix(R, so.reduced_op(n, lamp), prv) for n in fam['names']]
    return prv, jspec(mats, fam['combo'], dps)


def string_images(R, lamp, nu, W):
    """for i = 1..N-1: exact images E_i^{k_i} e_c, c in W, as dicts {target index: Fraction}  (k_i = lam'_i - lam'_{i+1} + 1; void if the image weight is not a weight of V_lam)"""
    N = R.N
    out = []
    for i in range(1, N):
        k = lamp[i - 1] - lamp[i] + 1
        out.append([so.power_apply(R, (i, i + 1), k, {int(c): Fr(1)}) for c in W])
    return out


def string_violations(U, imgs, w, thr):
    """for every eigenvector k (column of U) and string i: is E_i^{k_i} u_k != 0 ?  returns [set of violated strings (1-based) per eigenvector], and the largest violation sizes"""
    res = []
    for k in range(w):
        u = [U[c, k] for c in range(w)] if U is not None else [mp.mpf(1)]
        nrm = max(abs(x) for x in u)
        viol = set()
        for i, imi in enumerate(imgs, start=1):
            acc = {}
            for c in range(w):
                for r, v in imi[c].items():
                    acc[r] = acc.get(r, 0) + u[c] * _mpf_fr(v)
            if acc and max(abs(x) for x in acc.values()) > thr * nrm:
                viol.add(i)
        res.append(viol)
    return res


# ======================================================================================== 4. THE WEYL / BGG TEST (Q3, Q6)
def weyl_sector(lam, lamp, mu, fam=None, dps=50, onestep=True, maxlen=None, variant='rule'):
    old = mp.mp.dps
    mp.mp.dps = dps + 10
    try:
        return _weyl_sector(lam, lamp, mu, fam, dps, onestep, maxlen, variant)
    finally:
        mp.mp.dps = old


def _weyl_sector(lam, lamp, mu, fam, dps, onestep, maxlen, variant='rule'):
    """all 24 (N!) Weyl translates of the sector mu: joint spectra of V_red(w.lam') on V_lam[mu - w.lam'] (50 digits), triples matched across spaces (tolerance 10^-(dps-20)),
    signed counts, comparison with the physical spectrum, depth-truncated variants, one-step analysis (strings violated by each eigenvector of the identity space)."""
    lam, lamp, mu = tuple(lam), tuple(lamp), tuple(mu)
    N = len(lam)
    fam = fam or FAM[N]
    R = gr.GLNRep.of(lam)
    tol = mp.mpf(10) ** (-(dps - 20))
    spaces = []
    items = []
    t0 = time.perf_counter()
    for p, sg, ln in WEYL[N]:
        wd = dot_action(p, lamp) if variant != 'plain' else tuple(lamp[p[i]] for i in range(N))     # 'plain': the ordinary Weyl action w(lam') instead of the dot action (sensitivity test)
        if variant == 'nosign':
            sg = 1                                                                                   # 'nosign': all translates counted with +1 (sensitivity test)
        nu = tuple(mu[i] - wd[i] for i in range(N))
        W = R.wsp.get(nu)
        if W is None:
            continue
        ops = {n: so.reduced_op(n, wd) for n in fam['names']}
        mats = so.weight_space_mats(R, ops, W)
        js = jspec([mats[n] for n in fam['names']], fam['combo'], dps)
        sp_ = dict(p=p, sgn=sg, len=ln, wd=wd, nu=nu, W=W, w=len(W), js=js)
        spaces.append(sp_)
        for k, t in enumerate(js['trip']):
            items.append((t, (len(spaces) - 1, k)))
    groups = group_triples(items, tol)
    prv, pj = physical_spec(lam, lamp, mu, fam, dps)
    m = prv.m
    # signed counts and the indicator
    for g in groups:
        g['signed'] = sum(spaces[si]['sgn'] for si, k in g['members'])
        g['phys'] = sum(1 for t in pj['trip'] if max(abs(a - b) for a, b in zip(t, g['trip'])) < tol)
        g['lens'] = sorted(spaces[si]['len'] for si, k in g['members'])
    phys_found = sum(1 for g in groups if g['phys'])
    fails = [g for g in groups if g['signed'] != g['phys']]
    # depth-truncated variants: keep only translates of length <= L
    maxL = max((sp_['len'] for sp_ in spaces), default=0)
    depth_ok = {}
    for L in range(0, (maxlen if maxlen is not None else N * (N - 1) // 2) + 1):
        ok = True
        for g in groups:
            s = sum(spaces[si]['sgn'] for si, k in g['members'] if spaces[si]['len'] <= L)
            if s != g['phys']:
                ok = False
                break
        depth_ok[L] = ok
    # signed lattice count (Q6)
    dims = [(sp_['len'], sp_['sgn'], sp_['w']) for sp_ in spaces]
    mult_signed = sum(sg * w for ln, sg, w in dims)
    out = dict(mu=mu, m=m, spaces=spaces, groups=groups, nphys=len(pj['trip']), phys_found=phys_found, fails=fails, depth_ok=depth_ok, mult_signed=mult_signed,
               sep=min_group_sep(groups), maxdev=max((g['dev'] for g in groups), default=mp.mpf(0)), jok=all(s['js']['ok'] for s in spaces) and pj['ok'],
               jgap=min((s['js']['gap'] for s in spaces), default=mp.inf), time=time.perf_counter() - t0, pj=pj, prv=prv)
    # one-step analysis on the identity space
    out['onestep'] = None
    if onestep:
        ident = next((si for si, sp_ in enumerate(spaces) if sp_['len'] == 0), None)
        out['onestep'] = None
        if ident is not None:
            sp0 = spaces[ident]
            imgs = string_images(R, lamp, sp0['nu'], sp0['W'])
            viol = string_violations(sp0['js']['U'], imgs, sp0['w'], tol)
            refl = {i: simple_refl(N, i) for i in range(1, N)}
            rows = []
            for k, t in enumerate(sp0['js']['trip']):
                g = next(g for g in groups if (ident, k) in g['members'])
                S = viol[k]
                T = set()
                for i in range(1, N):
                    for si, kk in g['members']:
                        if spaces[si]['p'] == refl[i]:
                            T.add(i)
                rows.append(dict(S=S, T=T, phys=bool(g['phys']), signed=g['signed'], lens=g['lens']))
            out['onestep'] = rows
    return out


def weyl_product(lam, lamp, kind='rs', fam=None, dps=50, verbose=True, listing=True, maxlist=100, onestep=True, variant='rule'):
    """run weyl_sector over the sectors of a product; returns (list of sector results, summary dict)"""
    lam, lamp = tuple(lam), tuple(lamp)
    N = len(lam)
    fam = fam or FAM[N]
    secs = sector_list(lam, lamp, kind)
    rs = gr.racah_speiser(lam, lamp)
    res = []
    t0 = time.perf_counter()
    for mu, m in secs:
        res.append(weyl_sector(lam, lamp, mu, fam, dps, onestep=onestep, variant=variant))
    T = time.perf_counter() - t0
    S = summarise_weyl(res, lam, lamp, rs)
    S['time'] = T
    S['kind'] = kind
    if verbose:
        print_weyl_summary(S, res, lam, lamp, listing=listing, maxlist=maxlist)
    return res, S


def summarise_weyl(res, lam, lamp, rs):
    N = len(lam)
    S = Counter()
    S['sectors'] = len(res)
    depth = Counter()
    for r in res:
        S['m_total'] += r['m']
        S['spaces'] += len(r['spaces'])
        S['groups'] += len(r['groups'])
        S['w_ident'] += next((sp_['w'] for sp_ in r['spaces'] if sp_['len'] == 0), 0)
        ok = (not r['fails']) and r['phys_found'] == r['nphys'] == r['m']
        S['ok'] += ok
        S['fail_sectors'] += (not ok)
        S['fail_groups'] += len(r['fails'])
        S['jfail'] += (not r['jok'])
        S['mult_ok'] += (r['mult_signed'] == r['m'] == rs.get(r['mu'], 0))
        S['nonsimple'] += sum(1 for sp_ in r['spaces'] if sp_['js']['gap'] < mp.mpf(10) ** -20)
        for L, v in r['depth_ok'].items():
            depth[L] += bool(v)
        # minimal depth needed
        Lmin = next((L for L in sorted(r['depth_ok']) if r['depth_ok'][L]), None)
        S['Lmin_%s' % Lmin] += 1
        if r['onestep'] is not None:
            for row in r['onestep']:
                S['os_tau'] += 1
                S['os_S_eq_T'] += (row['S'] == row['T'])
                S['os_S_sub_T'] += (row['S'] <= row['T'])
                S['os_multi'] += (len(row['S']) >= 2)
                S['os_phys_consistent'] += (row['phys'] == (len(row['S']) == 0))
                if not row['phys']:
                    S['os_spurious'] += 1
                    S['os_spur_single'] += (len(row['S']) == 1)
                    S['os_spur_multi'] += (len(row['S']) >= 2)
    S['depth'] = dict(depth)
    S['maxlen'] = max((sp_['len'] for r in res for sp_ in r['spaces']), default=0)
    S['maxsep'] = min((r['sep'] for r in res), default=mp.inf)
    S['maxdev'] = max((r['maxdev'] for r in res), default=mp.mpf(0))
    return S


def fmt_S(S):
    return ('sectors %d (sum m = %d, %d translate spaces, %d distinct triples): signed count = indicator on %d, FAIL on %d sectors (%d triples); eigendecomposition not ok in %d sectors; '
            'non-simple spaces %d; signed lattice count = m on %d/%d; closest distinct triples %s, worst identification deviation %s' %
            (S['sectors'], S['m_total'], S['spaces'], S['groups'], S['ok'], S['fail_sectors'], S['fail_groups'], S['jfail'], S['nonsimple'], S['mult_ok'], S['sectors'],
             mp.nstr(S['maxsep'], 4), mp.nstr(S['maxdev'], 3)))


def print_weyl_summary(S, res, lam, lamp, listing=True, maxlist=100):
    N = len(lam)
    tab('  %s  [%s sectors]: %s  (%.1f s)' % (prod_str(lam, lamp), S['kind'], fmt_S(S), S['time']))
    dd = S['depth']
    tab('     truncation depth (sum over translates of length <= L only): sectors OK at L = ' + ', '.join('%d: %d' % (L, dd[L]) for L in sorted(dd)) + ' (of %d)' % S['sectors'])
    if S['os_tau']:
        tab('     one-step analysis (identity-space eigenvectors, %d triples): physical <=> no string violated: %d/%d; violated set S subset of inherited set T (BGG): %d/%d; S = T: %d/%d;'
            ' spurious triples %d: one string violated %d, two or more %d' %
            (S['os_tau'], S['os_phys_consistent'], S['os_tau'], S['os_S_sub_T'], S['os_tau'], S['os_S_eq_T'], S['os_tau'], S['os_spurious'], S['os_spur_single'], S['os_spur_multi']))
    if listing:
        tab('     per sector:  mu  m  w(identity)  #translate spaces (lengths: dims)  #triples  signed-count failures  min depth  [fail details]')
        for r in res[:maxlist]:
            lens = Counter()
            for sp_ in r['spaces']:
                lens[sp_['len']] += sp_['w']
            Lmin = next((L for L in sorted(r['depth_ok']) if r['depth_ok'][L]), None)
            w0 = next((sp_['w'] for sp_ in r['spaces'] if sp_['len'] == 0), 0)
            line = '     %-14s m=%-3d w=%-3d spaces=%-2d dims by length %-34s triples=%-3d fails=%d depth=%s' % (mu_str(r['mu']), r['m'], w0, len(r['spaces']),
                   ','.join('%d:%d' % (L, lens[L]) for L in sorted(lens)), len(r['groups']), len(r['fails']), Lmin)
            if r['fails']:
                line += '   FAIL: ' + '; '.join('tau=(%s) signed=%d phys=%d lens=%s' % (','.join(mp.nstr(mp.re(x), 8) for x in g['trip']), g['signed'], g['phys'], g['lens']) for g in r['fails'][:3])
            tab(line)
        if len(res) > maxlist:
            tab('     ... (%d more sectors)' % (len(res) - maxlist))


def bgg_signature_table(res, title):
    """for each triple of each sector: the multiset of lengths of the translates it occurs in, and whether it is physical -> Counter"""
    c = Counter()
    for r in res:
        for g in r['groups']:
            c[(tuple(g['lens']), bool(g['phys']), g['signed'])] += 1
    tab('  %s: occurrence pattern of a triple over the Weyl translates (lengths of the translates containing it), physical?, signed count -> number of triples' % title)
    for (lens, ph, sg), n in sorted(c.items(), key=lambda kv: (-kv[1], kv[0])):
        tab('     lengths %-22s physical=%-5s signed=%+d : %d' % (lens, ph, sg, n))
    return c


# ======================================================================================== 5. STENCILS (Q1)
# GT entries of a gl(4) pattern in the order (m11, m12, m22, m13, m23, m33); row 4 = lam never moves.  Free lattice coordinates of a weight space: (m12, m13, m23).
ENT = {(1, 1): 0, (1, 2): 1, (2, 2): 2, (1, 3): 3, (2, 3): 4, (3, 3): 5}
ROWS = {1: (0,), 2: (1, 2), 3: (3, 4, 5)}


def coords(pat):
    """free lattice coordinates (m12, m13, m23) of a gl(4) GT pattern ((m11,), (m12, m22), (m13, m23, m33), lam)"""
    return (pat[1][0], pat[2][0], pat[2][1])


def entries(pat):
    return (pat[0][0], pat[1][0], pat[1][1], pat[2][0], pat[2][1], pat[2][2])


def _unit(k, sign):
    v = [0] * 6
    v[k] = sign
    return tuple(v)


def _sumset(A, B):
    return {tuple(a + b for a, b in zip(x, y)) for x in A for y in B}


_GM = {}


def gen_moves(g):
    """the set of 6-vectors (changes of the entries (m11,m12,m22,m13,m23,m33)) with which the generator E_ij (acting on a gl(4) GT pattern, GLNRep's definition: adjacent = GT
    formulas, E_ij = [E_{i,j-1}, E_{j-1,j}] for j > i+1, E_ij = [E_{i,i-1}, E_{i-1,j}] for i > j+1) can change a pattern: E_{k,k+1} raises ONE entry of row k, E_{k+1,k} lowers ONE entry of
    row k, the diagonal ones do not move, the long ones move one entry in each of the rows i..j-1 (a sumset of the moves of their two factors)"""
    if g in _GM:
        return _GM[g]
    i, j = g
    if i == j:
        r = {(0,) * 6}
    elif j == i + 1:
        r = {_unit(ENT[(a, i)], +1) for a in range(1, i + 1)}
    elif i == j + 1:
        r = {_unit(ENT[(a, j)], -1) for a in range(1, j + 1)}
    elif j > i + 1:
        r = _sumset(gen_moves((i, j - 1)), gen_moves((j - 1, j)))
    else:
        r = _sumset(gen_moves((i, i - 1)), gen_moves((i - 1, j)))
    _GM[g] = r
    return r


def word_moves(word):
    r = {(0,) * 6}
    for g in word:
        r = _sumset(r, gen_moves(g))
    return r


def weight_zero(v):
    """the row sums (m11), (m12+m22), (m13+m23+m33) do not change"""
    return v[0] == 0 and v[1] + v[2] == 0 and v[3] + v[4] + v[5] == 0


def project(v):
    """6-vector of entry changes (weight zero) -> change of the free coordinates (dm12, dm13, dm23)"""
    return (v[1], v[3], v[4])


def structural_stencil(name):
    """{displacement (dm12, dm13, dm23): number of universal reduced words of the operator that can produce it} (weight-zero part of the sumset of the words' moves)"""
    U = so.universal_reduction(name)
    cnt = Counter()
    for w in U:
        ds = {project(v) for v in word_moves(w) if weight_zero(v)}
        for d in ds:
            cnt[d] += 1
    return cnt


def data_stencil(lam, lamps, names=so.TRIPLE, min_w=2):
    """nonzero matrix elements of the reduced operators on EVERY weight space of V_lam, for each lam' in `lamps`: {name: Counter(displacement -> number of nonzero entries)}, and the
    set of weight-space sizes met"""
    lam = tuple(lam)
    R = gr.GLNRep.of(lam)
    out = {n: Counter() for n in names}
    for lamp in lamps:
        ops = {n: so.reduced_op(n, lamp) for n in names}
        for nu, W in R.wsp.items():
            if len(W) < min_w:
                continue
            mats = so.weight_space_mats(R, ops, W)
            cs = [coords(R.pats[int(c)]) for c in W]
            for n in names:
                M = mats[n]
                for a in range(len(W)):
                    for b in range(len(W)):
                        if M[a][b] != 0:
                            out[n][tuple(x - y for x, y in zip(cs[a], cs[b]))] += 1
    return out


GENERIC_LAMPS = [(Fr(1, 3), Fr(5, 7), Fr(-2), Fr(4, 5)), (Fr(7, 2), Fr(-1, 3), Fr(2), Fr(11, 4)), (Fr(-3), Fr(5), Fr(1, 2), Fr(-7, 3))]


def stage_stencil(verbose=True):
    hdr('Q1  STENCILS: lattice displacements (dm12, dm13, dm23) of the nonzero matrix elements of V_red,3 / V_red,4a / V_red,4b on the weight-space lattice')
    tab('  Structural stencil = weight-zero part of the Minkowski sum of the moves of the generators of each universal reduced word (E_{k,k+1}: raises ONE entry of row k; E_{k+1,k}: lowers one entry')
    tab('  of row k; E_ii: no move; E_ij, |i-j| > 1: one entry in each of the rows min..max-1, as GLNRep builds them by commutators).  A displacement is the change of (m12, m13, m23); the other entries')
    tab('  follow from the weight (m22 = const - m12, m33 = const - m13 - m23).  gl(3) analogue: {0, +-(1,-1)} in (m1, m2).')
    st = {}
    for n in so.TRIPLE:
        st[n] = structural_stencil(n)
    tab('')
    tab('  %-6s %-9s %-14s %-12s %s' % ('op', 'words', 'stencil size', 'max |d|_1', 'symmetric under d -> -d'))
    for n in so.TRIPLE:
        c = st[n]
        tab('  %-6s %-9d %-14d %-12d %s' % (n, len(so.universal_reduction(n)), len(c), max(sum(abs(x) for x in d) for d in c), all((tuple(-x for x in d) in c) for d in c)))
    tab('')
    tab('  structural stencil of V3 (displacement : number of words that can produce it):')
    for d, k in sorted(st['V3'].items()):
        tab('     %-14s %d' % (d, k))
    # data
    lams = [(2, 1, 0, 0), (3, 2, 1, 0), (4, 2, 1, 0), (4, 3, 1, 0), (4, 2, 2, 0), (5, 3, 2, 1), (6, 4, 2, 0)]
    dstat = {}
    tab('')
    tab('  Data stencil: union over all weight spaces of dim >= 2 of V_lam, (a) GENERIC lam\' (3 rational non-dominant weights, e.g. (1/3, 5/7, -2, 4/5)), (b) the dominant lam\' = lam of the product')
    tab('  (the stencil of a dominant lam\' can be smaller: coefficients that vanish identically at special lam\').  Sizes (number of displacements incl. 0) and containment in the structural stencil:')
    tab('  %-12s %-5s | %-26s | %-26s | %s' % ('V_lam', 'op', 'generic lam\' (size, in struct.)', 'dominant lam\' = lam (size)', 'max nonzero entries per displacement'))
    allgen = {n: Counter() for n in so.TRIPLE}
    alldom = {n: Counter() for n in so.TRIPLE}
    for lam in lams:
        dg = data_stencil(lam, GENERIC_LAMPS)
        dd = data_stencil(lam, [lam])
        for n in so.TRIPLE:
            allgen[n].update(dg[n])
            alldom[n].update(dd[n])
            inside = all(d in st[n] for d in dg[n]) and all(d in st[n] for d in dd[n])
            tab('  %-12s %-5s | %3d  %-20s | %3d                        | %d' % (lam_str(lam), n, len(dg[n]), 'contained' if inside else 'NOT CONTAINED', len(dd[n]), max(dg[n].values(), default=0)))
        dstat[lam] = (dg, dd)
    tab('')
    tab('  Union over the reps above:')
    for n in so.TRIPLE:
        g, dm = allgen[n], alldom[n]
        sym = all((tuple(-x for x in d) in g) for d in g)
        tab('   %-4s structural %3d | generic-lam\' data %3d (equal to structural: %s; symmetric: %s) | dominant-lam\' data %3d' % (n, len(st[n]), len(g), set(g) == set(st[n]), sym, len(dm)))
        miss = sorted(set(st[n]) - set(g))
        if miss:
            tab('        structural displacements never seen with generic lam\': %s' % miss)
    tab('')
    tab('  Stencils by operator (data, generic lam\'): displacement lists')
    for n in so.TRIPLE:
        tab('   %s: %s' % (n, sorted(allgen[n])))
    return st, allgen, alldom


# ======================================================================================== 6. SYMBOLIC CLOSED FORMS (Q2)
# Generic GT formulas in "l-coordinates" y_ik = m_ik - i (rows 1..3), y_j4 = lam_j - j, applied to a generic pattern through the universal reduced words.  A coefficient is a sum of TERMS
# c * prod_f f^e over affine forms f (GT factors y_a - y_b + const, Cartan sums, lam'-monomials); terms with equal exponent vectors are merged, nothing is expanded.  The simplification into
# closed forms is done by exact line reconstruction: N(t) = M(P + t v) * D(P + t v) is a polynomial in t (D = product of the denominator forms), interpolated exactly, factored over Q;
# along the 14 axes the linear factors give their roots, which are clustered into affine forms (a form with coefficient vector c and value F at P has the roots -F/c_i).
SYM_VARS = ['y11', 'y12', 'y22', 'y13', 'y23', 'y33', 'y14', 'y24', 'y34', 'y44', 'q1', 'q2', 'q3', 'q4']
SYM_LABEL = ['m11', 'm12', 'm22', 'm13', 'm23', 'm33', 'l1', 'l2', 'l3', 'l4', 'q1', 'q2', 'q3', 'q4']       # y_v = label_v - SYM_SHIFT[v]
SYM_SHIFT = [1, 1, 2, 1, 2, 3, 1, 2, 3, 4, 0, 0, 0, 0]
_EV = {(1, 1): 0, (1, 2): 1, (2, 2): 2, (1, 3): 3, (2, 3): 4, (3, 3): 5, (1, 4): 6, (2, 4): 7, (3, 4): 8, (4, 4): 9}
_DI = {(1, 1): 0, (1, 2): 1, (2, 2): 2, (1, 3): 3, (2, 3): 4, (3, 3): 5}
_FORM_ID = {}
_FORMS = []


def _fid(form):
    i = _FORM_ID.get(form)
    if i is None:
        i = _FORM_ID[form] = len(_FORMS)
        _FORMS.append(form)
    return i


def _mkform(coefs, c0):
    """affine form sum coefs[v] y_v + c0 -> (id of the canonical form (first nonzero coefficient positive), sign s with form = s * canonical)"""
    items = tuple(sorted((v, c) for v, c in coefs.items() if c != 0))
    sgn = 1
    if items[0][1] < 0:
        sgn = -1
        items = tuple((v, -c) for v, c in items)
        c0 = -c0
    return _fid((c0, items)), sgn


def _ydiff(a, da, b, db):
    return _mkform({_EV[a]: 1, _EV[b]: -1}, da - db)


def _shift(disp, e):
    return 0 if e[1] == 4 else disp[_DI[e]]


def _key_mul(key, factors):
    d = dict(key)
    for fid, e in factors:
        v = d.get(fid, 0) + e
        if v == 0:
            d.pop(fid, None)
        else:
            d[fid] = v
    return tuple(sorted(d.items()))


_ADJC = {}


def _adj_terms(k, up, disp):
    kk = (k, up, disp)
    if kk in _ADJC:
        return _ADJC[kk]
    out = []
    for i in range(1, k + 1):
        e_i = (i, k)
        di = _shift(disp, e_i)
        coef = Fr(-1) if up else Fr(1)
        fac = []
        rng_num = range(1, k + 2) if up else range(1, k)
        rown = (k + 1) if up else (k - 1)
        for j in rng_num:
            e_j = (j, rown)
            fid, s = _ydiff(e_i, di, e_j, _shift(disp, e_j))
            fac.append((fid, 1))
            coef *= s
        for j in range(1, k + 1):
            if j != i:
                e_j = (j, k)
                fid, s = _ydiff(e_i, di, e_j, _shift(disp, e_j))
                fac.append((fid, -1))
                coef *= s
        nd = list(disp)
        nd[_DI[e_i]] += 1 if up else -1
        out.append((tuple(nd), coef, fac))
    _ADJC[kk] = out
    return out


def _cartan(k, disp):
    lst = [(1, _EV[(i, k)]) for i in range(1, k + 1)]
    const = sum(_shift(disp, (i, k)) + i for i in range(1, k + 1))
    if k > 1:
        lst += [(-1, _EV[(i, k - 1)]) for i in range(1, k)]
        const -= sum(_shift(disp, (i, k - 1)) + i for i in range(1, k))
    return lst, const


def _clean(st):
    return {d: {k: c for k, c in t.items() if c != 0} for d, t in st.items() if any(c != 0 for c in t.values())}


def _add_into(acc, st, s=1):
    for d, terms in st.items():
        o = acc.setdefault(d, defaultdict(Fr))
        for key, c in terms.items():
            o[key] += s * c


def _apply_adj(st, k, up):
    out = {}
    for d, terms in st.items():
        for nd, coef, fac in _adj_terms(k, up, d):
            o = out.setdefault(nd, defaultdict(Fr))
            for key, c in terms.items():
                o[_key_mul(key, fac)] += c * coef
    return _clean(out)


def _apply_gen(st, g):
    i, j = g
    if i == j:
        out = {}
        for d, terms in st.items():
            lst, const = _cartan(i, d)
            o = out.setdefault(d, defaultdict(Fr))
            for key, c in terms.items():
                if const:
                    o[key] += c * const
                for cv, v in lst:
                    o[_key_mul(key, [(_fid((0, ((v, 1),))), 1)])] += c * cv
        return _clean(out)
    if j == i + 1:
        return _apply_adj(st, i, True)
    if i == j + 1:
        return _apply_adj(st, j, False)
    A, B = ((i, j - 1), (j - 1, j)) if j > i + 1 else ((i, i - 1), (i - 1, j))
    out = {}
    _add_into(out, _apply_gen(_apply_gen(st, B), A))
    _add_into(out, _apply_gen(_apply_gen(st, A), B), -1)
    return _clean(out)


def _apply_word(word, st):
    for g in reversed(word):
        st = _apply_gen(st, g)
        if not st:
            break
    return st


_SYMT = {}


def sym_terms(name):
    """exact generic matrix elements of the universal reduced operator `name` on a generic GT pattern: {displacement (6-vector of entry changes): {term key: Fraction}}; a term key is a
    sorted tuple of (form id, exponent); lam'-monomials enter as single-variable forms (q_a)^e.  Only terms that did not cancel are kept."""
    if name in _SYMT:
        return _SYMT[name]
    U = so.universal_reduction(name)
    start = {(0,) * 6: {(): Fr(1)}}
    total = {}
    for w, poly in U.items():
        res = _apply_word(w, start)
        for mono, coef in poly.items():
            qf = [(_fid((0, ((10 + a, 1),))), e) for a, e in enumerate(mono) if e]
            for d, terms in res.items():
                o = total.setdefault(d, defaultdict(Fr))
                for key, c in terms.items():
                    o[_key_mul(key, qf)] += c * Fr(coef)
    _SYMT[name] = _clean(total)
    return _SYMT[name]


def form_value(fid, pt):
    c0, items = _FORMS[fid]
    s = c0
    for v, c in items:
        s += c * pt[v]
    return s


def eval_terms(terms, pt, cache=None):
    """exact value of a term sum at the point pt (14 numbers: y_11 .. y_44, q_1..q_4); ZeroDivisionError if a denominator form vanishes"""
    fv = {} if cache is None else cache
    tot = Fr(0)
    for key, c in terms.items():
        p = c
        for fid, e in key:
            val = fv.get(fid)
            if val is None:
                val = fv[fid] = form_value(fid, pt)
            p = p * (val ** e)
        tot += p
    return tot


def form_str(f, mode='label'):
    """affine form (c0, ((var, coef), ...)) in the y/q variables -> text in the GT labels (m_ik, l_j = lam_j, q_j = lam'_j)"""
    c0, items = f
    parts = []
    const = Fr(c0)
    for v, c in items:
        const -= c * SYM_SHIFT[v]
        parts.append((c, SYM_LABEL[v]))
    s = ''
    for c, nm in parts:
        c = Fr(c)
        sg = '-' if c < 0 else '+'
        a = abs(c)
        s += ' %s %s%s' % (sg, ('' if a == 1 else '%s*' % a), nm)
    if const != 0:
        s += ' %s %s' % ('-' if const < 0 else '+', abs(const))
    s = s.strip()
    if s.startswith('+ '):
        s = s[2:]
    return '(' + s + ')'


def pt_of_pattern(pat, lam, lamp):
    """the 14 numbers of the generic point for a concrete GT pattern ((m11,), (m12, m22), (m13, m23, m33), lam), site-1 weights lam and lam' = lamp"""
    m = entries(pat)
    y = [m[0] - 1, m[1] - 1, m[2] - 2, m[3] - 1, m[4] - 2, m[5] - 3]
    y += [lam[j] - (j + 1) for j in range(4)]
    return [Fr(x) for x in y] + [Fr(x) for x in lamp]


def term_stats(name):
    T = sym_terms(name)
    return {d: len(t) for d, t in T.items()}


def interp_newton(xs, ys):
    """exact Newton interpolation -> sympy Poly in t"""
    t = sp.Symbol('t')
    n = len(xs)
    coef = list(ys)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i - 1]) / (xs[i] - xs[i - j])
    p = sp.Rational(coef[-1].numerator, coef[-1].denominator)
    for i in range(n - 2, -1, -1):
        p = p * (t - sp.Rational(xs[i].numerator, xs[i].denominator)) + sp.Rational(coef[i].numerator, coef[i].denominator)
    return sp.Poly(sp.expand(p), t)


class LineProbe:
    """N(x) = M(x) D(x) for a term sum M (D = product of the denominator forms with their maximal exponents); exact line reconstruction"""

    def __init__(self, terms):
        self.terms = terms
        ids = set()
        for key in terms:
            for fid, e in key:
                ids.add(fid)
        self.ids = sorted(ids)
        self.emin = {fid: min(dict(key).get(fid, 0) for key in terms) for fid in self.ids}
        self.den = {fid: -e for fid, e in self.emin.items() if e < 0}
        self.degb = max(sum(dict(key).get(fid, 0) + max(0, -self.emin[fid]) for fid in self.ids) for key in terms)
        self.t = sp.Symbol('t')

    def N(self, pt):
        val = eval_terms(self.terms, pt)
        D = Fr(1)
        for fid, e in self.den.items():
            D *= form_value(fid, pt) ** e
        return val * D

    def line(self, P, v, extra=2):
        xs = [Fr(i) for i in range(self.degb + 1 + extra)]
        ys = [self.N([P[i] + x * v[i] for i in range(14)]) for x in xs]
        poly = interp_newton(xs[:self.degb + 1], ys[:self.degb + 1])
        ok = all(poly.eval(sp.Rational(x.numerator, x.denominator)) == sp.Rational(y.numerator, y.denominator) for x, y in zip(xs[self.degb + 1:], ys[self.degb + 1:]))
        return poly, ok

    def factor_line(self, P, v):
        poly, ok = self.line(P, v)
        if poly.is_zero:
            return None, [], [], ok, -1
        c, fl = sp.factor_list(poly.as_expr(), self.t)
        lin, nonlin = [], []
        for f, m in fl:
            pf = sp.Poly(f, self.t)
            if pf.degree() == 1:
                a, b = pf.all_coeffs()
                lin.append((Fr(int(-b.p) * int(a.q), int(b.q) * int(a.p)), m))
            else:
                nonlin.append((pf.degree(), m))
        return c, lin, nonlin, ok, poly.degree()


DEBUG_ID = False


def rand_point(rng):
    """random base point whose coordinates have the fractional parts 2^i / 10007 (a Sidon-type choice: no form with small integer coefficients vanishes or coincides in value with another)"""
    return [Fr(rng.randint(-60, 60)) + Fr(2 ** i, 10007) for i in range(14)]


def identify_forms(pr, rng, maxcoef=4):
    """clusters the axis roots of N (the numerator of M_Delta over the maximal denominator) into candidate affine forms, keeps those that are genuine factors on two generic lines, and
    measures what is left.  Returns dict(zero, forms = [(c0, items, mult)] genuine linear factors (items = ((var, int coef), ...), first coefficient positive), degN, qdeg = degree of the
    non-linear remainder Q (0: N is a product of linear forms), qaxis = degree of Q along each axis, incomplete = linear roots on generic lines not explained by the forms, kappa (if qdeg == 0),
    check)"""
    P = rand_point(rng)
    axis, axdeg = [], []
    for i in range(14):
        v = [Fr(0)] * 14
        v[i] = Fr(1)
        c, lin, nonlin, ok, deg = pr.factor_line(P, v)
        if c is None:
            return dict(zero=True)
        axis.append(lin)
        axdeg.append(deg)
    pool = sorted(((abs(r), i, r, m) for i, lst in enumerate(axis) for r, m in lst), key=lambda x: -x[0])
    used = [False] * len(pool)
    cands = []
    for a0, e0 in enumerate(pool):
        if used[a0]:
            continue
        _, i0, r0, m0 = e0
        used[a0] = True
        coefs = {i0: Fr(1)}
        for a1 in range(a0 + 1, len(pool)):
            if used[a1]:
                continue
            _, i1, r1, m1 = pool[a1]
            if i1 in coefs or m1 != m0:
                continue
            ratio = r0 / r1
            if ratio.denominator == 1 and abs(ratio) <= maxcoef and ratio != 0:
                coefs[i1] = ratio
                used[a1] = True
        c0 = -r0 - sum(c * P[i] for i, c in coefs.items())
        den = 1
        for c in list(coefs.values()) + [c0]:
            den = den * c.denominator // math.gcd(den, c.denominator)
        ic = {i: int(c * den) for i, c in coefs.items()}
        ic0 = c0 * den
        g = 0
        for c in ic.values():
            g = math.gcd(g, abs(c))
        ic = {i: c // g for i, c in ic.items()}
        ic0 = ic0 / g
        first = min(ic)
        if ic[first] < 0:
            ic = {i: -c for i, c in ic.items()}
            ic0 = -ic0
        ic0 = int(ic0) if Fr(ic0).denominator == 1 else Fr(ic0)
        cands.append((ic0, tuple(sorted(ic.items())), m0))
    lines = []
    for _ in range(2):
        P2 = rand_point(rng)
        for _try in range(200):                       # small slopes (cheap factorisation) such that no candidate form is constant along the line
            v2 = [Fr(rng.randint(1, 13)) for _ in range(14)]
            if all(sum(c_ * v2[i] for i, c_ in items) != 0 for c0, items, m in cands):
                break
        c, lin, nonlin, ok, deg = pr.factor_line(P2, v2)
        lines.append((P2, v2, Counter({r: m for r, m in lin}), deg, ok))
    genuine = []
    for c0, items, m in cands:
        mm = m
        for P2, v2, act, deg, ok in lines:
            Lp = c0 + sum(c_ * P2[i] for i, c_ in items)
            lv = sum(c_ * v2[i] for i, c_ in items)
            if lv == 0:
                mm = 0
                break
            mm = min(mm, act.get(-Lp / lv, 0))
        if mm > 0:
            genuine.append((c0, items, mm))
    incomplete = 0
    for P2, v2, act, deg, ok in lines:
        pred = Counter()
        for c0, items, m in genuine:
            Lp = c0 + sum(c_ * P2[i] for i, c_ in items)
            lv = sum(c_ * v2[i] for i, c_ in items)
            pred[-Lp / lv] += m
        incomplete += sum((act - pred).values())
    degN = lines[0][3]
    nlin = sum(m for _, _, m in genuine)
    qdeg = degN - nlin
    qaxis = [axdeg[i] - sum(m for c0, items, m in genuine for v_, c_ in items if v_ == i) for i in range(14)]
    res = dict(zero=False, forms=genuine, degN=degN, qdeg=qdeg, qaxis=qaxis, incomplete=incomplete, lines_ok=all(l[4] for l in lines) and lines[0][3] == lines[1][3])
    res['check'] = res['lines_ok'] and incomplete == 0
    if qdeg == 0:
        ks = []
        for _ in range(3):
            P3 = rand_point(rng)
            prod = Fr(1)
            for c0, items, m in genuine:
                prod *= (c0 + sum(c_ * P3[i] for i, c_ in items)) ** m
            ks.append(pr.N(P3) / prod)
        res['kappa'] = ks[0]
        res['check'] = res['check'] and ks[0] == ks[1] == ks[2]
    return res


def classify_displacement(d):
    """(row-2 move, row-3 move) of the displacement 6-vector"""
    return (d[1], (d[3], d[4], d[5]))


def stage_symbolic_probe(names=('V3',), verbose=True, seed=11):
    rng = random.Random(seed)
    out = {}
    for name in names:
        T = sym_terms(name)
        tab('  %s: %d displacement classes, %d terms in all (before any simplification)' % (name, len(T), sum(len(t) for t in T.values())))
        for d in sorted(T):
            t0 = time.perf_counter()
            pr = LineProbe(T[d])
            res = identify_forms(pr, rng)
            res['terms'] = len(T[d])
            res['den'] = dict(pr.den)
            res['disp'] = d
            out[(name, d)] = res
            if verbose:
                if res['zero']:
                    tab('    %-4s d=%-12s terms %4d  ZERO' % (name, project(d), len(T[d])))
                else:
                    tab('    %-4s d=%-12s terms %4d  den: %2d forms (deg %2d)  N: degree %2d = %2d linear factors (%2d distinct) + non-linear part of degree %d  [check %s]  %.1fs' %
                        (name, project(d), len(T[d]), len(pr.den), sum(pr.den.values()), res['degN'], sum(m for _, _, m in res['forms']), len(res['forms']), res['qdeg'], res['check'], time.perf_counter() - t0))
    return out


_CF = {}


def closed_forms(name, seed=11):
    """for every displacement class d (6-vector) of the universal operator `name`: the identification record (forms, kappa, qdeg, qaxis, ...), the net factor list (numerator forms minus the
    denominator forms), and the LineProbe (for the evaluation of the non-factorising classes).  Cached."""
    if name in _CF:
        return _CF[name]
    rng = random.Random(seed)
    T = sym_terms(name)
    out = {}
    for d in sorted(T):
        pr = LineProbe(T[d])
        res = identify_forms(pr, rng)
        for _retry in range(3):                    # a failed identification (coincidences of the random base point) is retried with fresh points
            if res['zero'] or res['check']:
                break
            res = identify_forms(pr, rng)
        res['terms'] = len(T[d])
        res['pr'] = pr
        if not res['zero']:
            net = {}
            for c0, items, m in res['forms']:
                net[(c0, items)] = net.get((c0, items), 0) + m
            for fid, e in pr.den.items():
                f = _FORMS[fid]
                net[f] = net.get(f, 0) - e
            res['net'] = [(f, e) for f, e in net.items() if e != 0]
        out[d] = res
    _CF[name] = out
    return out


def closed_form_one(name, d, seed=11):
    """identification record (as in closed_forms) for the single displacement class d (6-vector) of `name`"""
    rng = random.Random(seed)
    T = sym_terms(name)
    pr = LineProbe(T[d])
    res = identify_forms(pr, rng)
    for _retry in range(3):
        if res['zero'] or res['check']:
            break
        res = identify_forms(pr, rng)
    res['terms'] = len(T[d])
    res['pr'] = pr
    if not res['zero']:
        net = {}
        for c0, items, m in res['forms']:
            net[(c0, items)] = net.get((c0, items), 0) + m
        for fid, e in pr.den.items():
            f = _FORMS[fid]
            net[f] = net.get(f, 0) - e
        res['net'] = [(f, e) for f, e in net.items() if e != 0]
    return res


def form_value_f(f, pt):
    c0, items = f
    s_ = Fr(c0)
    for v, c in items:
        s_ += c * pt[v]
    return s_


def eval_net(rec, pt):
    """kappa * prod f^e of a fully factorising class at the point pt (ZeroDivisionError at a pole)"""
    v = rec['kappa']
    for f, e in rec['net']:
        v = v * (form_value_f(f, pt) ** e)
    return v


def _lagrange(xs, ys, x0):
    tot = Fr(0)
    for i in range(len(xs)):
        term = ys[i]
        for j in range(len(xs)):
            if j != i:
                term = term * (x0 - xs[j]) / (xs[i] - xs[j])
        tot += term
    return tot


def _poly_coeffs(xs, ys):
    """ascending coefficients (Fractions) of the interpolating polynomial through the points (Newton form expanded)"""
    n = len(xs)
    c = list(ys)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            c[i] = (c[i] - c[i - 1]) / (xs[i] - xs[i - j])
    poly = [c[-1]]
    for i in range(n - 2, -1, -1):
        new = [Fr(0)] * (len(poly) + 1)
        for k, a in enumerate(poly):
            new[k + 1] += a
            new[k] -= a * xs[i]
        new[0] += c[i]
        poly = new
    return poly


def eval_class(rec, pt, rng):
    """exact value of the matrix element class at pt: product form if the class factorises, else Q~(pt) * prod_net f^e with the polynomial Q~ = N / prod L^m of degree qdeg
    reconstructed on a generic line pt + s v (no zero denominators at s != 0; the extra point checks the degree).  Forms with f(pt) = 0 (removable poles at boundary patterns) are
    handled by the exact Laurent series in s: M(pt) = lim_{s -> 0}."""
    if rec['zero']:
        return Fr(0)
    pr = rec.get('pr')
    qd = rec['qdeg']
    allforms = [f for f, e in rec['net']] + [(c0, items) for c0, items, m in rec['forms']] + ([_FORMS[fid] for fid in pr.den] if pr is not None else [])
    for _try in range(500):                       # a direction along which no form is constant
        v = [Fr(rng.randint(1, 13)) for _ in range(14)]
        if all(sum(c * v[i] for i, c in f[1]) != 0 for f in allforms):
            break
    if qd == 0:
        qco = [Fr(rec['kappa'])]
    else:
        ss = [Fr(i) + Fr(1, 10007) for i in range(1, qd + 3)]
        vals = []
        for sx in ss:
            q = [pt[i] + sx * v[i] for i in range(14)]
            prod = Fr(1)
            for c0, items, m in rec['forms']:
                prod *= form_value_f((c0, items), q) ** m
            vals.append(eval_Qexp(rec, q) if 'Qexp' in rec else pr.N(q) / prod)
        if _lagrange(ss[:qd + 1], vals[:qd + 1], ss[qd + 1]) != vals[qd + 1]:
            raise ValueError('degree of the non-linear part is not %d' % qd)
        qco = _poly_coeffs(ss[:qd + 1], vals[:qd + 1])
        # Q~(pt + s v) = sum qco[j] s^j  (the interpolation variable IS s)
    # the product of the net forms along the line: forms vanishing at pt contribute (s * l)^e
    pw = 0
    const = Fr(1)
    for f, e in rec['net']:
        fv = form_value_f(f, pt)
        if fv == 0:
            lv = sum(c * v[i] for i, c in f[1])
            pw += e
            const *= lv ** e
        else:
            const *= fv ** e
    # sum_j qco[j] s^(j+pw): finite limit needs qco[j] = 0 for j + pw < 0
    for j, a in enumerate(qco):
        if j + pw < 0 and a != 0:
            raise ZeroDivisionError('genuine pole')
    j0 = -pw
    if j0 < 0:
        return Fr(0)
    return (qco[j0] if j0 < len(qco) else Fr(0)) * const


def reconstruct_q(rec, seed=3):
    """explicit non-linear part Q~ = N / prod L^m of a class with small qdeg: exact interpolation of the polynomial in its support variables (monomials of total degree <= qdeg and per-variable degree
    <= qaxis), solved exactly over Q (DomainMatrix rref).  Stores rec['Qexp'] = (support variable indices, {exponent tuple: Fraction}) in the y/q variables and rec['Qstr'] (GT labels)."""
    import math
    from sympy import QQ
    from sympy.polys.matrices import DomainMatrix
    rng = random.Random(seed)
    pr = rec['pr']
    supp = [i for i, q in enumerate(rec['qaxis']) if q]
    qd = rec['qdeg']
    mons = [e for e in itertools.product(*[range(rec['qaxis'][i] + 1) for i in supp]) if sum(e) <= qd]
    base = rand_point(rng)
    rows, rhs = [], []
    for _ in range(len(mons) + 8):
        pt = list(base)
        for i in supp:
            pt[i] = Fr(rng.randint(-30, 30)) + Fr(2 ** i, 10007)
        prod = Fr(1)
        for c0, items, m in rec['forms']:
            prod *= form_value_f((c0, items), pt) ** m
        rhs.append(pr.N(pt) / prod)
        rows.append([math.prod([pt[i] ** k for i, k in zip(supp, e)]) for e in mons])
    A = DomainMatrix([[QQ(x.numerator, x.denominator) for x in r] + [QQ(y.numerator, y.denominator)] for r, y in zip(rows, rhs)], (len(rows), len(mons) + 1), QQ)
    R_, piv = A.rref()
    if len(mons) in piv:
        raise ValueError('inconsistent: wrong degree / support')
    Rl = R_.to_Matrix()
    coef = {}
    for r_, c in enumerate(piv):
        v = Rl[r_, len(mons)]
        if v != 0:
            coef[mons[c]] = Fr(int(v.p), int(v.q))
    rec['Qexp'] = (supp, coef)
    ys = sp.symbols([SYM_VARS[i] for i in supp])
    lab = {ys[k]: sp.Symbol(SYM_LABEL[i]) - SYM_SHIFT[i] for k, i in enumerate(supp)}
    Qy = sum(sp.Rational(c.numerator, c.denominator) * sp.Mul(*[y ** e for y, e in zip(ys, ex)]) for ex, c in coef.items())
    rec['Qstr'] = str(sp.expand(Qy.subs(lab)))
    rec['Qterms'] = len(sp.Poly(sp.expand(Qy.subs(lab)), *[sp.Symbol(SYM_LABEL[i]) for i in supp]).terms())
    return rec


def eval_Qexp(rec, pt):
    supp, coef = rec['Qexp']
    tot = Fr(0)
    for ex, c in coef.items():
        t = c
        for i, k in zip(supp, ex):
            if k:
                t *= pt[i] ** k
        tot += t
    return tot


def print_closed_forms(name, CF):
    tab("  %s: closed forms by displacement class d = (dm12, dm13, dm23); M_d(x) = kappa * prod (forms)^e   [GT labels m_ik, l_j = lam_j, q_j = lam'_j]" % name)
    fs = lambda L: ' '.join('%s%s' % (f, ('^%d' % e if e > 1 else '')) for f, e in L) or '1'
    for d in sorted(CF, key=lambda k: (abs(k[1]), project(k))):
        rec = CF[d]
        if rec['zero']:
            tab('   d=%-12s  = 0 identically' % (project(d),))
            continue
        num = [(form_str(f), e) for f, e in rec['net'] if e > 0]
        den = [(form_str(f), -e) for f, e in rec['net'] if e < 0]
        if rec['qdeg'] == 0:
            tab('   d=%-12s  kappa=%+d  numerator: %s' % (project(d), rec['kappa'], fs(num)))
            tab('   %-16s  denominator: %s' % ('', fs(den)))
        else:
            sup = ','.join('%s^%d' % (SYM_LABEL[i], q) for i, q in enumerate(rec['qaxis']) if q)
            tab('   d=%-12s  NOT a product of linear forms: [linear numerator factors] x Q, Q a polynomial of total degree %d (degrees per variable: %s)' % (project(d), rec['qdeg'], sup))
            tab('   %-16s  linear numerator: %s' % ('', fs(num)))
            tab('   %-16s  denominator: %s' % ('', fs(den)))


def verify_closed(name, lam, lamps, rng, which='all', verbose=False, maxdim=None, CF=None):
    """held-out test of the closed forms against G3's exact matrices (GLNRep + reduced_op): every weight space W of V_lam, every lam' in lamps, every source pattern x in W and every
    displacement class d of the stencil: closed form (product form for the factorising classes, `eval_class` for the others) == G3's entry (0 if the target is no GT pattern); also every
    nonzero entry of G3's matrices must lie in the stencil.  which = 'F' (factorising classes only) or 'all'."""
    lam = tuple(lam)
    R = gr.GLNRep.of(lam)
    CF = CF if CF is not None else closed_forms(name)
    stats = Counter()
    for lamp in lamps:
        op = so.reduced_op(name, lamp)
        for nu, W in sorted(R.wsp.items()):
            if maxdim and len(W) > maxdim:
                continue
            M = so.weight_space_mats(R, {name: op}, W)[name]
            pos = {int(c): a for a, c in enumerate(W)}
            stats['spaces'] += 1
            covered = set()
            for xi, cx in enumerate(W):
                pat = R.pats[int(cx)]
                ent = entries(pat)
                pt = pt_of_pattern(pat, lam, lamp)
                for d, rec in CF.items():
                    if which == 'F' and (not rec['zero'] and rec['qdeg'] != 0):
                        continue
                    tgt = tuple(a + b for a, b in zip(ent, d))
                    tp = ((tgt[0],), (tgt[1], tgt[2]), (tgt[3], tgt[4], tgt[5]), lam)
                    ti = R.index.get(tp)
                    if ti is None or ti not in pos:
                        stats['invalid_target'] += 1               # the target is no GT pattern: the matrix has no such entry (the formal coefficient need not vanish there)
                        continue
                    expected = M[pos[ti]][xi]
                    covered.add((pos[ti], xi))
                    try:
                        val = eval_class(rec, pt, rng)
                    except ZeroDivisionError:
                        stats['zerodiv'] += 1
                        stats['zerodiv_expected_zero'] += (expected == 0)
                        continue
                    stats['tests'] += 1
                    stats['nonzero_tests'] += (expected != 0)
                    if val != expected:
                        stats['bad'] += 1
                        if verbose and stats['bad'] < 5:
                            print('   MISMATCH', name, lam, lamp, d, pat, val, expected)
            if which == 'all':
                for a in range(len(W)):
                    for b in range(len(W)):
                        if M[a][b] != 0 and (a, b) not in covered:
                            stats['uncovered_nonzero'] += 1
    return stats


def explicit_q_block(names, verify=False, rng=None):
    """prints the explicit non-linear parts Q~ of the row-3 moves with dm12 = 0 (V3: all six; V4a, V4b: the class (0,-1,0) as an example) and optionally verifies them against G3
    (the verification then uses the explicit polynomial, not the term-sum oracle)"""
    tab('')
    tab('  Explicit non-linear parts Q~ (the polynomial left after dividing N by its linear factors), reconstructed exactly by interpolation in the support variables (V3: all six row-3 moves with dm12 = 0;')
    tab("  V4a, V4b: the class (0,-1,0) as an example); printed in the GT labels m_ik, l_j = lam_j, q_j = lam'_j:")
    sel = {}
    for name in names:
        if name in _CF:
            CF = _CF[name]
        else:
            T = sym_terms(name)
            CF = {}
            for d in sorted(T, key=lambda k: project(k)):
                if name == 'V3' or project(d) == (0, -1, 0):
                    if project(d)[0] == 0 and project(d) != (0, 0, 0):
                        CF[d] = closed_form_one(name, d)
        for d in sorted(CF, key=lambda k: project(k)):
            rec = CF[d]
            if rec['zero'] or rec['qdeg'] == 0 or rec['qdeg'] > 4 or project(d)[0] != 0 or project(d) == (0, 0, 0):
                continue
            if name != 'V3' and project(d) != (0, -1, 0):
                continue
            t0 = time.perf_counter()
            reconstruct_q(rec)
            sel.setdefault(name, {})[d] = rec
            tab('   %s d=%s: Q~ has %d terms, total degree %d: %s   [%.0f s]' % (name, project(d), rec['Qterms'], rec['qdeg'], rec['Qstr'], time.perf_counter() - t0))
    if verify:
        rng = rng or random.Random(9)
        tab('')
        tab("  Verification of the classes with an explicit polynomial (closed form = linear factors x explicit Q~, evaluated with the exact Laurent handling of removable poles) against G3, all weight spaces:")
        for name, CFs in sel.items():
            for lam in ((3, 2, 1, 0), (4, 3, 1, 0), (5, 3, 2, 1)):
                for lamp_set, tag in (([lam], "dominant lam'=lam"), ([(Fr(1, 3), Fr(5, 7), Fr(-2), Fr(4, 5))], 'generic rational')):
                    t0 = time.perf_counter()
                    st = verify_closed(name, lam, lamp_set, rng, which='all', CF=CFs)
                    tab('   %-4s %-12s %-20s %d classes: %d tests (%d nonzero entries), %d mismatches, %d pole skips  [%.1fs]' % (name, lam_str(lam), tag, len(CFs), st['tests'], st['nonzero_tests'], st['bad'], st['zerodiv'], time.perf_counter() - t0))
    return sel


def f_classes():
    """the 12 displacement 6-vectors (m11, m12, m22, m13, m23, m33) with dm12 = +-1 and a non-zero row-3 move (one entry of row 3 raised, one lowered)"""
    out = []
    for s in (1, -1):
        for i in (1, 2, 3):
            for j in (1, 2, 3):
                if i != j:
                    d3 = [0, 0, 0]
                    d3[i - 1] += 1
                    d3[j - 1] -= 1
                    out.append((0, s, -s, d3[0], d3[1], d3[2]))
    return out


def universal_F(op, pt, d):
    """THE closed form of the 12 fully factorising classes in l-coordinates y_ik = m_ik - i, Y_s = lam_s - s (pt = the 14 numbers y11, y12, y22, y13, y23, y33, Y1..Y4, q1..q4):
       M^{V3}_d = (y11 - y_b)(y11 - y_a + q1 - q2)(y_a - y_j)(y_a - y_j + 1)(y_a - y_k) prod_s (y_i - Y_s) / [(y_a - y_b)(y_a - y_b + 1)(y_i - y_j)(y_i - y_j + 1)(y_j - y_k)(y_i - y_k)],
       M^{V4a}_d = (y_k + 3) M^{V3}_d,   M^{V4b}_d = -(y_a - y_i + 2 y_k - q3 + 5) M^{V3}_d,
       a, b = the raised / lowered entry of row 2, i, j, k = the raised / lowered / spectator entry of row 3"""
    row3 = {1: 3, 2: 4, 3: 5}
    a, b = (1, 2) if d[1] == 1 else (2, 1)
    d3 = {1: d[3], 2: d[4], 3: d[5]}
    i = next(k for k, v in d3.items() if v == 1)
    j = next(k for k, v in d3.items() if v == -1)
    k = next(k for k, v in d3.items() if v == 0)
    ya, yb, yi, yj, yk = pt[a], pt[b], pt[row3[i]], pt[row3[j]], pt[row3[k]]
    num = (pt[0] - yb) * (pt[0] - ya + pt[10] - pt[11]) * (ya - yj) * (ya - yj + 1) * (ya - yk)
    for s_ in range(4):
        num *= (yi - pt[6 + s_])
    den = (ya - yb) * (ya - yb + 1) * (yi - yj) * (yi - yj + 1) * (yj - yk) * (yi - yk)
    m3 = num / den
    if op == 'V3':
        return m3
    if op == 'V4a':
        return (yk + 3) * m3
    if op == 'V4b':
        return -(ya - yi + 2 * yk - pt[12] + 5) * m3
    raise KeyError(op)


def verify_universal(lam, lamps, ops=('V3', 'V4a', 'V4b')):
    """the universal formula against G3's exact matrices: every weight space of V_lam, every lam', every source pattern and every one of the 12 classes with a valid target pattern"""
    lam = tuple(lam)
    R = gr.GLNRep.of(lam)
    stats = Counter()
    for lamp in lamps:
        opsd = {n: so.reduced_op(n, lamp) for n in ops}
        for nu, W in sorted(R.wsp.items()):
            mats = so.weight_space_mats(R, opsd, W)
            pos = {int(c): a for a, c in enumerate(W)}
            for xi, cx in enumerate(W):
                pat = R.pats[int(cx)]
                ent = entries(pat)
                pt = pt_of_pattern(pat, lam, lamp)
                for d in f_classes():
                    tgt = tuple(a + b for a, b in zip(ent, d))
                    ti = R.index.get(((tgt[0],), (tgt[1], tgt[2]), (tgt[3], tgt[4], tgt[5]), lam))
                    if ti is None or ti not in pos:
                        continue
                    for n in ops:
                        try:
                            val = universal_F(n, pt, d)
                        except ZeroDivisionError:
                            stats['zerodiv'] += 1
                            continue
                        exp_ = mats[n][pos[ti]][xi]
                        stats['tests'] += 1
                        stats['nonzero'] += (exp_ != 0)
                        stats['bad'] += (val != exp_)
    return stats


def stage_univ():
    hdr('Q2 (continued)  ONE FORMULA FOR THE 12 FACTORISING CLASSES OF V3, V4a, V4b')
    tab('  In l-coordinates y_ik = m_ik - i (rows 1..3), Y_s = lam_s - s, q = lam\' ; a, b = the RAISED / LOWERED entry of row 2, i, j, k = the RAISED / LOWERED / SPECTATOR entry of row 3 of the displacement:')
    tab('    M^{V3}  = (y11 - y_b)(y11 - y_a + q1 - q2)(y_a - y_j)(y_a - y_j + 1)(y_a - y_k) prod_{s=1..4}(y_i - Y_s) / [ (y_a - y_b)(y_a - y_b + 1)(y_i - y_j)(y_i - y_j + 1)(y_j - y_k)(y_i - y_k) ]')
    tab('    M^{V4a} = (y_k + 3) M^{V3},      M^{V4b} = - (y_a - y_i + 2 y_k - q3 + 5) M^{V3}')
    tab('  (found from the 36 factor lists of the symbolic stage; the sign of V3 is +1 for all 12 classes in this orientation; the relations to V3 hold exactly for all 12 classes, 120 random-point tests of the term sums)')
    lamps_gen = [(Fr(1, 3), Fr(5, 7), Fr(-2), Fr(4, 5)), (Fr(-3), Fr(5), Fr(1, 2), Fr(-7, 3))]
    tot = Counter()
    tab("  %-12s %-34s %s" % ('V_lam', "lam'", 'tests (valid source and target) / nonzero entries / mismatches / pole skips   (V3, V4a, V4b together)'))
    for lam in ((3, 2, 1, 0), (4, 3, 1, 0), (5, 3, 2, 1), (6, 4, 2, 0), (4, 2, 2, 0), (5, 4, 2, 1), (6, 3, 2, 1)):
        for lamp_set, tag in (([lam], "dominant lam' = lam"), (lamps_gen, 'generic rational (2 values)')):
            t0 = time.perf_counter()
            st = verify_universal(lam, lamp_set)
            tab("  %-12s %-34s %d / %d / %d / %d   [%.1fs]" % (lam_str(lam), tag, st['tests'], st['nonzero'], st['bad'], st['zerodiv'], time.perf_counter() - t0))
            for k_, v_ in st.items():
                tot[k_] += v_
    tab('  TOTAL: %d tests (%d nonzero entries), %d mismatches, %d pole skips' % (tot['tests'], tot['nonzero'], tot['bad'], tot['zerodiv']))
    return tot


def stage_symq():
    hdr('Q2 (continued)  EXPLICIT NON-LINEAR PARTS OF THE ROW-3 MOVES (dm12 = 0)')
    explicit_q_block(('V3', 'V4a', 'V4b'), verify=True)


def stage_symbolic(names=('V3', 'V4a', 'V4b'), heldout=((3, 2, 1, 0), (4, 3, 1, 0), (5, 3, 2, 1), (6, 4, 2, 0)), nonf_rep=(3, 2, 1, 0), seed=5):
    hdr('Q2  CLOSED-FORM MATRIX ELEMENTS (symbolic GT formulas through the universal reduced words; exact line reconstruction; held-out verification)')
    tab("  Method: the universal reduced operator (G3: 39 / 78 / 85 words with coefficients in Q[lam']) is applied to a GENERIC GT pattern with symbolic labels; every intermediate state is kept (no")
    tab('  boundary treatment).  Coefficients are sums of terms c * prod(affine forms)^e (GT factors y_a - y_b + c, Cartan sums, lam\'-monomials), never expanded.  The sum over the terms of one')
    tab('  displacement class is simplified by exact line reconstruction (see the module docstring): N = M * D (D = maximal denominator) is a polynomial in the line parameter, interpolated exactly and')
    tab('  factored over Q; the axis roots are clustered into affine forms, kept if they are genuine factors on generic lines; the remainder is the non-linear part Q.')
    tab("  Variables: y_ik = m_ik - i (rows 1..3), y_j4 = lam_j - j; forms are printed in the GT labels m_ik, l_j = lam_j (site-1 irrep) and q_j = lam'_j (site-2 highest weight).")
    rng = random.Random(seed)
    summ = {}
    for name in names:
        t0 = time.perf_counter()
        CF = closed_forms(name)
        tt = time.perf_counter() - t0
        nF = sum(1 for r in CF.values() if not r['zero'] and r['qdeg'] == 0)
        nZ = sum(1 for r in CF.values() if r['zero'])
        nN = len(CF) - nF - nZ
        ok_id = all(r['check'] for r in CF.values() if not r['zero'])
        tab('')
        tab('  %s: %d displacement classes with terms (%d terms in all) [%.0f s]: %d identically zero after summation (|dm12| = 2), %d fully factorising, %d with a non-linear part; identification checks all True: %s' %
            (name, len(CF), sum(r['terms'] for r in CF.values()), tt, nZ, nF, nN, ok_id))
        summ[name] = dict(F=nF, N=nN, Z=nZ)
        print_closed_forms(name, CF)
    explicit_q_block(names)
    tab('')
    tab("  Held-out verification against G3 (exact Fractions, all weight spaces of V_lam, every source pattern, every displacement class; lam' = lam (dominant) and generic rational lam'):")
    tab("  %-6s %-12s %-8s %-30s %s" % ('op', 'V_lam', 'classes', "lam'", 'tests (valid source and valid target) / of which nonzero entries / mismatches / pole skips'))
    lamps_gen = [(Fr(1, 3), Fr(5, 7), Fr(-2), Fr(4, 5)), (Fr(-3), Fr(5), Fr(1, 2), Fr(-7, 3))]
    tot = Counter()
    for name in names:
        for lam in heldout:
            for lamp_set, tag in (([lam], "dominant lam'=lam"), (lamps_gen, 'generic rational (2 values)')):
                t0 = time.perf_counter()
                st = verify_closed(name, lam, lamp_set, rng, which='F')
                tab('  %-6s %-12s %-8s %-30s %d / %d / %d / %d   [%d weight spaces, %.1fs]' % (name, lam_str(lam), 'F only', tag, st['tests'], st['nonzero_tests'], st['bad'], st['zerodiv'], st['spaces'], time.perf_counter() - t0))
                tot['tests'] += st['tests']
                tot['bad'] += st['bad']
                tot['nonzero'] += st['nonzero_tests']
                tot['zerodiv'] += st['zerodiv']
    tab('  factorising classes, all held-out runs: %d tests (%d nonzero entries), %d mismatches, %d pole skips' % (tot['tests'], tot['nonzero'], tot['bad'], tot['zerodiv']))
    tab('')
    tab('  All classes (non-factorising ones via the polynomial extrapolation of Q, completeness of the stencil), rep %s:' % lam_str(nonf_rep))
    tab("  %-6s %-12s %-8s %-30s %s" % ('op', 'V_lam', 'classes', "lam'", 'tests / nonzero entries / mismatches / uncovered nonzero entries of G3 (not in the stencil)'))
    for name in names:
        for lamp_set, tag in (([nonf_rep], "dominant lam'=lam"), (lamps_gen[:1], 'generic rational')):
            t0 = time.perf_counter()
            st = verify_closed(name, nonf_rep, lamp_set, rng, which='all')
            tab('  %-6s %-12s %-8s %-30s %d / %d / %d / %d   [%d weight spaces, zero-denominator skips %d (of which the true entry is 0: %d), %.1fs]' % (name, lam_str(nonf_rep), 'all', tag, st['tests'], st['nonzero_tests'], st['bad'], st['uncovered_nonzero'], st['spaces'], st['zerodiv'], st['zerodiv_expected_zero'], time.perf_counter() - t0))
    return summ


# ======================================================================================== 7. per-sector data of the identity space (rules, chains, polynomiality)
def ident_sector(lam, lamp, mu, fam=None, dps=50):
    """the whole weight space W = V_lam[mu - lam'] of a sector: exact matrices of the family, joint spectrum (50 digits), which eigenvectors are physical (joint triple in Spec(M_mu)),
    and which strings E_i^{k_i} each eigenvector violates"""
    lam, lamp, mu = tuple(lam), tuple(lamp), tuple(mu)
    N = len(lam)
    fam = fam or FAM[N]
    old = mp.mp.dps
    mp.mp.dps = dps + 10
    try:
        R = gr.GLNRep.of(lam)
        nu = tuple(mu[i] - lamp[i] for i in range(N))
        W = R.wsp.get(nu)
        if W is None:
            return None
        ops = {n: so.reduced_op(n, lamp) for n in fam['names']}
        mats = so.weight_space_mats(R, ops, W)
        js = jspec([mats[n] for n in fam['names']], fam['combo'], dps)
        prv, pj = physical_spec(lam, lamp, mu, fam, dps)
        tol = mp.mpf(10) ** (-(dps - 20))
        phys = [any(max(abs(a - b) for a, b in zip(t, tp)) < tol for tp in pj['trip']) for t in js['trip']]
        imgs = string_images(R, lamp, nu, W)
        S = string_violations(js['U'], imgs, len(W), tol)
        # which strings are non-void (the image weight is a weight of V_lam)
        active = set()
        for i in range(1, N):
            k = lamp[i - 1] - lamp[i] + 1
            tw = list(nu)
            tw[i - 1] += k
            tw[i] -= k
            if tuple(tw) in R.wsp:
                active.add(i)
        return dict(lam=lam, lamp=lamp, mu=mu, R=R, W=W, w=len(W), m=prv.m, mats=mats, js=js, phys=phys, S=S, active=active, prv=prv, pj=pj, names=fam['names'], combo=fam['combo'])
    finally:
        mp.mp.dps = old


def gen_value(t, combo):
    return sum(mp.mpf(c.numerator) / c.denominator * x for c, x in zip(combo, t))


# ======================================================================================== 8. Q4: which eigenlines are physical?
def rank_orderings(sec):
    """for every ordering (V3, V4a, V4b, G2's generic combination): the values per eigenvector"""
    trip = sec['js']['trip']
    vals = {}
    for k, nm in enumerate(sec['names']):
        vals[nm] = [mp.re(t[k]) for t in trip]
    if len(sec['names']) == 3:
        vals['G'] = [mp.re(gen_value(t, sec['combo'])) for t in trip]
    return vals


def sep_direction(vals, phys, eps=mp.mpf(10) ** -20):
    """'top' if every physical value exceeds every spurious value, 'bottom' if every physical value is below every spurious one, 'none' otherwise; 'tie' if the ordering has no strict gap"""
    ph = [v for v, p in zip(vals, phys) if p]
    sp_ = [v for v, p in zip(vals, phys) if not p]
    if not sp_ or not ph:
        return 'all'
    if min(ph) - max(sp_) > eps:
        return 'top'
    if min(sp_) - max(ph) > eps:
        return 'bottom'
    if abs(min(ph) - max(sp_)) <= eps or abs(max(ph) - min(sp_)) <= eps:
        return 'tie'
    return 'none'


RULE_PRODUCTS = [((2, 1, 0, 0), (2, 1, 0, 0)), ((2, 1, 1, 0), (2, 1, 1, 0)), ((3, 2, 1, 0), (3, 2, 1, 0)), ((4, 2, 0, 0), (4, 2, 0, 0)), ((4, 2, 1, 0), (4, 2, 1, 0)), ((6, 4, 2, 0), (6, 4, 2, 0))]
RULE_HELD = [((3, 2, 1, 0), (4, 2, 1, 0)), ((4, 2, 2, 0), (4, 2, 2, 0)), ((4, 3, 1, 0), (4, 3, 1, 0))]


def rule_records(lam, lamp, dps=50, only_spurious=True):
    recs = []
    for mu, m in sector_list(lam, lamp, 'rs'):
        sec = ident_sector(lam, lamp, mu, dps=dps)
        if sec is None or (only_spurious and sec['w'] == sec['m']):
            continue
        vals = rank_orderings(sec)
        # classification: strings violated by the spurious eigenvectors
        Sset = set()
        for S_, p in zip(sec['S'], sec['phys']):
            if not p:
                Sset |= S_
        dirs = {nm: sep_direction(v, sec['phys']) for nm, v in vals.items()}
        recs.append(dict(prod=(lam, lamp), mu=mu, m=sec['m'], w=sec['w'], vals=vals, phys=sec['phys'], S=sec['S'], Sset=frozenset(Sset), dirs=dirs, active=sec['active'], trip=sec['js']['trip']))
    return recs


def rank_of(vals, k):
    """rank (1 = largest) of eigenvector k among vals (ties get the best rank)"""
    return 1 + sum(1 for v in vals if v - vals[k] > mp.mpf(10) ** -20)


def stage_rule(products=RULE_PRODUCTS, held=RULE_HELD, dps=50, verbose=True):
    hdr('Q4  WHICH EIGENLINES ARE PHYSICAL?  ranks of the physical triples in the orderings by V3, V4a, V4b and G2\'s generic combination, candidate rules, held-out products')
    tab('  Sectors with w > m only (w = m: every eigenline is physical).  rank 1 = largest eigenvalue.  S = strings E_i^{k_i} violated by the spurious eigenvectors of the sector (BGG: inherited from s_i.lam\').')
    tab("  direction of an ordering on a sector: top = every physical value above every spurious one, bottom = every physical below every spurious one, none = interleaved, tie = no strict gap.")
    allrecs = {}
    for lam, lamp in products + held:
        t0 = time.perf_counter()
        recs = rule_records(lam, lamp, dps)
        allrecs[(lam, lamp)] = recs
        tag = 'HELD-OUT' if (lam, lamp) in held else 'training'
        tab('')
        tab('  %s [%s]: %d sectors with w > m  (%.1f s)' % (prod_str(lam, lamp), tag, len(recs), time.perf_counter() - t0))
        if verbose and len(recs) <= 40:
            tab('     mu            m  w  active strings  S(spurious)   ranks of the PHYSICAL triples in V3 | V4a | V4b | G   (directions V3/V4a/V4b/G)')
            for r in recs:
                ph = [k for k, p in enumerate(r['phys']) if p]
                rk = {nm: sorted(rank_of(r['vals'][nm], k) for k in ph) for nm in r['vals']}
                tab('     %-13s %-2d %-2d %-15s %-13s %s | %s | %s | %s   (%s)' % (mu_str(r['mu']), r['m'], r['w'], sorted(r['active']), sorted(r['Sset']), rk['V3'], rk['V4a'], rk['V4b'], rk['G'],
                    '/'.join(r['dirs'][nm] for nm in ('V3', 'V4a', 'V4b', 'G'))))
    # aggregate by direction and by class (active strings A = intrinsic; violated strings S of the spurious modes = from the data)
    def agg(prods, label, keyf, kname):
        tab('')
        tab('  Direction statistics, %s products (sectors with w > m), per ordering and per %s:' % (label, kname))
        cls = defaultdict(lambda: defaultdict(Counter))
        tot = defaultdict(Counter)
        for pr in prods:
            for r in allrecs[pr]:
                for nm, d in r['dirs'].items():
                    cls[keyf(r)][nm][d] += 1
                    tot[nm][d] += 1
        tab('     %-12s | %s' % (kname, ' | '.join('%-34s' % nm for nm in ('V3', 'V4a', 'V4b', 'G'))))
        for c in sorted(cls):
            tab('     %-12s | %s' % (c, ' | '.join('%-34s' % ', '.join('%s %d' % (d, n) for d, n in sorted(cls[c][nm].items())) for nm in ('V3', 'V4a', 'V4b', 'G'))))
        tab('     %-12s | %s' % ('ALL', ' | '.join('%-34s' % ', '.join('%s %d' % (d, n) for d, n in sorted(tot[nm].items())) for nm in ('V3', 'V4a', 'V4b', 'G'))))
    keyA = lambda r: tuple(sorted(r['active']))
    keyS = lambda r: tuple(sorted(r['Sset']))
    agg(products, 'TRAINING', keyA, 'active strings A')
    agg(held, 'HELD-OUT', keyA, 'active strings A')
    agg(products, 'TRAINING', keyS, 'violated strings S')
    # linear separability: physical = top m of  l = (a, b, c) . (v3, v4a, v4b)  ?   LP per class: maximise the margin of  l . (tau_phys - tau_spur) >= margin * |tau_phys - tau_spur|
    from scipy.optimize import linprog

    def feats(t, deg):
        v = [float(mp.re(t[0])), float(mp.re(t[1])), float(mp.re(t[2]))]
        f = list(v)
        if deg >= 2:
            f += [v[i] * v[j] for i in range(3) for j in range(i, 3)]
        return np.array(f)

    def separate(recs, deg=1):
        rows = []
        for r in recs:
            ph = [k for k, p in enumerate(r['phys']) if p]
            sp_ = [k for k, p in enumerate(r['phys']) if not p]
            for a in ph:
                for b in sp_:
                    d = feats(r['trip'][a], deg) - feats(r['trip'][b], deg)
                    nd = np.linalg.norm(d)
                    if nd > 1e-12:
                        rows.append(d / nd)
        if not rows:
            return None
        A = np.array(rows)
        nf = A.shape[1]
        res = linprog(c=[0] * nf + [-1], A_ub=np.hstack([-A, np.ones((len(A), 1))]), b_ub=np.zeros(len(A)), bounds=[(-1, 1)] * nf + [(None, 1)], method='highs')
        if res.status != 0:
            return None
        return res.x[:nf], res.x[nf], len(rows)

    def apply_rule(l, recs, deg=1):
        ok = tot_ = 0
        for r in recs:
            vals = [float(np.dot(l, feats(t, deg))) for t in r['trip']]
            order = sorted(range(len(vals)), key=lambda k: -vals[k])
            top = set(order[:r['m']])
            tot_ += 1
            ok += (top == {k for k, p in enumerate(r['phys']) if p})
        return ok, tot_
    tab('')
    tab('  Linear separability (LP, margin maximised): is there a linear functional l = a V3 + b V4a + c V4b with  physical = top-m of l  on EVERY training sector of the class?  Then the same l on the held-out sectors')
    tab('  of the class (exactly the physical set picked?).  margin > 0: separable (the value is the smallest normalised gap l . (tau_phys - tau_spur) / |tau_phys - tau_spur|).')
    recs_tr = [r for pr in products for r in allrecs[pr]]
    recs_he = [r for pr in held for r in allrecs[pr]]
    for kname, keyf in (('active strings A', keyA), ('violated strings S', keyS), ('ALL sectors', lambda r: ())):
        tab('   classes by %s:' % kname)
        classes = sorted({keyf(r) for r in recs_tr + recs_he})
        for cl in classes:
            tr = [r for r in recs_tr if keyf(r) == cl]
            he = [r for r in recs_he if keyf(r) == cl]
            sep = separate(tr) if tr else None
            if sep is None:
                tab('     %-10s training %3d sectors: %s;  held-out sectors %d' % (cl, len(tr), 'no data' if not tr else 'LP failed', len(he)))
                continue
            l, mg, npairs = sep
            n = np.max(np.abs(l)) or 1
            l = l / n
            okh = apply_rule(l, he) if he else (0, 0)
            okt = apply_rule(l, tr)
            line = '     %-10s training %3d sectors (%4d pairs): %s  margin %+.4f  l = (%+.3f, %+.3f, %+.3f);  linear rule picks the physical set on training %d/%d, on HELD-OUT %d/%d' % (
                cl, len(tr), npairs, 'SEPARABLE' if mg > 1e-9 else 'not separable', mg, l[0], l[1], l[2], okt[0], okt[1], okh[0], okh[1])
            tab(line)
            if mg <= 1e-9:
                sq = separate(tr, 2)
                if sq is not None and sq[1] > 1e-9:
                    lq = sq[0] / (np.max(np.abs(sq[0])) or 1)
                    okq = apply_rule(lq, tr, 2)
                    okqh = apply_rule(lq, he, 2) if he else (0, 0)
                    tab('                quadratic functional (9 monomials of degree <= 2 in the triple): SEPARABLE, margin %+.4f; picks the physical set on training %d/%d, on HELD-OUT %d/%d' % (sq[1], okq[0], okq[1], okqh[0], okqh[1]))
                else:
                    tab('                quadratic functional (9 monomials of degree <= 2 in the triple): not separable')
    return allrecs


# ======================================================================================== 9. Q5(a): the sparsity graph
def graph_stats(M, names=None):
    """sparsity graph of an exact matrix (vertices = lattice points, edges = nonzero off-diagonal entries, symmetrised): components, path/tree, degrees, RCM bandwidth"""
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components, reverse_cuthill_mckee
    n = len(M)
    rows, cols = [], []
    asym = 0
    for a in range(n):
        for b in range(n):
            if a != b and M[a][b] != 0:
                rows.append(a)
                cols.append(b)
                if M[b][a] == 0:
                    asym += 1
    A = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    A = ((A + A.T) > 0).astype(np.int8).tocsr()
    nc, lab = connected_components(A, directed=False)
    deg = np.asarray(A.sum(axis=1)).ravel()
    edges = int(A.sum() // 2)
    comp_sizes = Counter(lab.tolist())
    # components are paths iff max degree <= 2 and edges = n - ncomp
    is_paths = (deg.max() <= 2 if n else True) and edges == n - nc
    perm = reverse_cuthill_mckee(A, symmetric_mode=True) if n > 1 else np.arange(n)
    inv = np.empty(n, dtype=int)
    inv[perm] = np.arange(n)
    bw = 0
    for a, b in zip(rows, cols):
        bw = max(bw, abs(int(inv[a]) - int(inv[b])))
    return dict(n=n, edges=edges, ncomp=nc, maxdeg=int(deg.max()) if n else 0, paths=is_paths, tree=(edges == n - nc), bw=bw, asym=asym)


def stage_chain(products=RULE_PRODUCTS + RULE_HELD, dps=50, verbose=True):
    hdr('Q5(a)  SPARSITY GRAPH OF THE REDUCED OPERATORS ON THE LATTICE R_nu (vertices = lattice points of W, edges = nonzero off-diagonal matrix elements)')
    tab('  For every weight space W = V_lam[mu - lam\'] of every irrep (identity translate): connected?  path (then tridiagonal in the path order)?  tree?  maximal degree, bandwidth after reverse Cuthill-McKee.')
    tot = defaultdict(Counter)
    for lam, lamp in products:
        R = gr.GLNRep.of(lam)
        N = 4
        agg = defaultdict(Counter)
        t0 = time.perf_counter()
        for mu, m in sector_list(lam, lamp, 'rs'):
            nu = tuple(mu[i] - lamp[i] for i in range(N))
            W = R.wsp.get(nu)
            if W is None or len(W) < 2:
                continue
            ops = {n: so.reduced_op(n, lamp) for n in so.TRIPLE}
            mats = so.weight_space_mats(R, ops, W)
            w = len(W)
            union = [[(mats['V3'][a][b] != 0 or mats['V4b'][a][b] != 0 or mats['V4a'][a][b] != 0) for b in range(w)] for a in range(w)]
            for tag, M in (('V3', mats['V3']), ('V4b', mats['V4b']), ('V3|V4b', [[mats['V3'][a][b] if mats['V3'][a][b] != 0 else mats['V4b'][a][b] for b in range(w)] for a in range(w)]), ('V3|V4a|V4b', [[Fr(1) if union[a][b] else Fr(0) for b in range(w)] for a in range(w)])):
                g = graph_stats(M)
                kind = 'm=w' if m == w else 'm<w'
                for key in ('conn', 'paths', 'tree'):
                    pass
                agg[(tag, kind)]['sectors'] += 1
                agg[(tag, kind)]['connected'] += (g['ncomp'] == 1)
                agg[(tag, kind)]['paths'] += g['paths']
                agg[(tag, kind)]['tree'] += g['tree']
                agg[(tag, kind)]['asym'] += (g['asym'] > 0)
                agg[(tag, kind)]['bw_sum'] += g['bw']
                agg[(tag, kind)]['bw_max'] = max(agg[(tag, kind)]['bw_max'], g['bw'])
                agg[(tag, kind)]['deg_max'] = max(agg[(tag, kind)]['deg_max'], g['maxdeg'])
                agg[(tag, kind)]['w_max'] = max(agg[(tag, kind)]['w_max'], w)
        tab('')
        tab('  %s  (%.1f s)   [weight spaces of dimension >= 2]' % (prod_str(lam, lamp), time.perf_counter() - t0))
        tab('     %-12s %-5s %8s %10s %7s %6s %8s %9s %8s %8s %9s' % ('operator', 'kind', 'spaces', 'connected', 'paths', 'trees', 'asym', 'mean bw', 'max bw', 'max deg', 'max w'))
        for (tag, kind), c in sorted(agg.items()):
            tab('     %-12s %-5s %8d %10d %7d %6d %8d %9.2f %8d %8d %9d' % (tag, kind, c['sectors'], c['connected'], c['paths'], c['tree'], c['asym'], c['bw_sum'] / max(1, c['sectors']), c['bw_max'], c['deg_max'], c['w_max']))
            for k_, v_ in c.items():
                tot[(tag, kind)][k_] += v_ if k_ not in ('bw_max', 'deg_max', 'w_max') else 0
                if k_ in ('bw_max', 'deg_max', 'w_max'):
                    tot[(tag, kind)][k_] = max(tot[(tag, kind)][k_], v_)
    tab('')
    tab('  TOTAL over all products above:')
    tab('     %-12s %-5s %8s %10s %7s %6s %8s %9s %8s %8s %9s' % ('operator', 'kind', 'spaces', 'connected', 'paths', 'trees', 'asym', 'mean bw', 'max bw', 'max deg', 'max w'))
    for (tag, kind), c in sorted(tot.items()):
        tab('     %-12s %-5s %8d %10d %7d %6d %8d %9.2f %8d %8d %9d' % (tag, kind, c['sectors'], c['connected'], c['paths'], c['tree'], c['asym'], c['bw_sum'] / max(1, c['sectors']), c['bw_max'], c['deg_max'], c['w_max']))
    # gl(3) comparison
    tab('')
    tab('  gl(3) comparison (the chain of jacobi_cg): V on V_Lam[mu - Lam], Lam = [8,4,0] and [6,4,2]:')
    for lam in ((8, 4, 0), (6, 4, 2)):
        R = gr.GLNRep.of(lam)
        c = Counter()
        for mu, m in sector_list(lam, lam, 'rs'):
            nu = tuple(mu[i] - lam[i] for i in range(3))
            W = R.wsp.get(nu)
            if W is None or len(W) < 2:
                continue
            M = so.weight_space_mats(R, {'V': so.reduced_op('V', lam)}, W)['V']
            g = graph_stats(M)
            c['spaces'] += 1
            c['connected'] += (g['ncomp'] == 1)
            c['paths'] += g['paths']
            c['bw_max'] = max(c['bw_max'], g['bw'])
            c['w_max'] = max(c['w_max'], len(W))
            c['m=w paths'] += (g['paths'] and m == len(W))
            c['m=w'] += (m == len(W))
        tab('     %s: %d spaces (dim >= 2), connected %d, paths (tridiagonal) %d, max bandwidth after RCM %d, max w %d; sectors with m = w: %d, of which paths %d' %
            (lam_str(lam), c['spaces'], c['connected'], c['paths'], c['bw_max'], c['w_max'], c['m=w'], c['m=w paths']))


# ======================================================================================== 10. Q5(b,c): cyclic coordinate vectors, graded Krylov degrees, polynomiality of the copies
def left_krylov(mats, x0):
    """row-vector Krylov space of e_{x0}^T under the commuting exact m x m matrices `mats` (e^T -> e^T A): graded by total degree.  Returns dict(cyclic, deg[y], poly[y] = {alpha: coeff}) with
    e_y^T = sum_alpha poly[y][alpha] * e_{x0}^T prod_t A_t^alpha_t  and deg[y] = minimal total degree of such an expression (None if e_y is not in the Krylov space)"""
    m = len(mats[0])
    r = len(mats)
    e0 = [Fr(0)] * m
    e0[x0] = Fr(1)
    vecs = {(0,) * r: e0}
    E = []                                    # echelon basis: (pivot, vector, combo {alpha: coeff})

    def reduce_vec(vec, combo):
        vec = list(vec)
        combo = dict(combo)
        for piv, bv, bc in E:
            f = vec[piv]
            if f != 0:
                vec = [a - f * b for a, b in zip(vec, bv)]
                for al, c in bc.items():
                    combo[al] = combo.get(al, 0) - f * c
        return vec, combo

    def insert(vec, combo):
        vec, combo = reduce_vec(vec, combo)
        piv = next((i for i, x in enumerate(vec) if x != 0), None)
        if piv is None:
            return False
        inv = 1 / vec[piv]
        E.append((piv, [x * inv for x in vec], {al: c * inv for al, c in combo.items()}))
        return True

    def rmul(vec, A):
        return [sum((vec[i] * A[i][j] for i in range(m) if vec[i] != 0), Fr(0)) for j in range(m)]

    deg = [None] * m
    poly = [None] * m
    insert(e0, {(0,) * r: Fr(1)})
    # E rows are reduced against earlier rows only; membership test below reduces fully
    def check_members(d):
        for y in range(m):
            if deg[y] is not None:
                continue
            ey = [Fr(0)] * m
            ey[y] = Fr(1)
            res, combo = reduce_vec(ey, {})
            if all(x == 0 for x in res):
                # e_y = - sum (coefficients picked up in `combo`) ... reduce_vec subtracts f*basis: ey - sum f_k b_k = 0  => ey = sum f_k b_k ; combo accumulated -f*c, so ey = -combo
                deg[y] = d
                poly[y] = {al: -c for al, c in combo.items() if c != 0}
    check_members(0)
    d = 0
    while any(x is None for x in deg):
        d += 1
        new = {}
        grew = False
        for al_prev, v_prev in vecs.items():
            if sum(al_prev) != d - 1:
                continue
            for t in range(r):
                al = tuple(al_prev[i] + (1 if i == t else 0) for i in range(r))
                if al in new:
                    continue
                new[al] = rmul(v_prev, mats[t])
        for al, v in new.items():
            if insert(v, {al: Fr(1)}):
                grew = True
        vecs.update(new)
        check_members(d)
        if not grew:
            break
    return dict(cyclic=all(x is not None for x in deg), deg=deg, poly=poly, dim=len(E))


def d_gen(m):
    """smallest total degree d such that the polynomials of degree <= d in three variables can interpolate any function on m generic points: C(d+3, 3) >= m"""
    d = 0
    while math.comb(d + 3, 3) < m:
        d += 1
    return d


def poly_eval(poly, trip):
    return sum((mp.mpf(c.numerator) / c.denominator * mp.fprod([t ** a for t, a in zip(trip, al)]) for al, c in poly.items()), mp.mpf(0))


BYM = defaultdict(Counter)


def stage_poly(products=RULE_PRODUCTS + RULE_HELD + [((5, 3, 1, 0), (5, 3, 1, 0)), ((5, 3, 2, 0), (5, 3, 2, 0))], dps=50, verbose=True):
    hdr('Q5(b,c)  POLYNOMIALITY OF THE SLICE COMPONENTS IN THE TRIPLE, CYCLIC COORDINATE VECTORS, GRADED KRYLOV DEGREES')
    tab('  On M_mu (m x m exact matrices in the PRV coordinates, m >= 2).  A coordinate x0 is CYCLIC iff e_x0^T generates the dual space under the commutative algebra of the triple, equivalently every')
    tab('  copy has a nonzero component x0.  Then every other component of every copy is a polynomial in the copy\'s triple with RATIONAL coefficients, u_j[y] / u_j[x0] = P_y(v3, v4a, v4b),')
    tab('  of minimal total degree deg(y) (graded Krylov degree, exact); also reported for the single generic operator G (degree <= m-1, always).  The polynomials are checked at 50 digits against G3\'s copies.')
    tot = Counter()
    for lam, lamp in products:
        t0 = time.perf_counter()
        c = Counter()
        rows = []
        for mu, m in sector_list(lam, lamp, 'rs'):
            if m < 2:
                continue
            prv, mats_d = so.triple_mats(lam, lamp, mu)
            mats = [mats_d[n] for n in so.TRIPLE]
            Gm = [[sum((Fr(cf) * A[i][j] for cf, A in zip(so.GENERIC, mats)), Fr(0)) for j in range(m)] for i in range(m)]
            cps = so.copies(lam, lamp, mu, dps=dps, prv=prv, mats=mats_d)
            c['sectors'] += 1
            c['m_total'] += m
            cyc = []
            best = None
            for x0 in range(m):
                kr = left_krylov(mats, x0)
                krG = left_krylov([Gm], x0)
                cyc.append((kr['cyclic'], krG['cyclic']))
                if kr['cyclic']:
                    dmax = max(kr['deg'])
                    if best is None or dmax < best[0]:
                        best = (dmax, x0, kr, krG)
            c['cyclic_coords'] += sum(1 for a, b in cyc if a)
            c['coords'] += m
            c['all_cyclic'] += all(a for a, b in cyc)
            c['some_cyclic'] += (best is not None)
            c['G_cyclic_equal'] += all(a == b for a, b in cyc)
            if best is None:
                rows.append((mu, m, 'NO cyclic coordinate', None))
                continue
            dmax, x0, kr, krG = best
            BYM[m][(dmax, max(krG['deg']))] += 1
            c['dmax_sum'] += dmax
            c['dmax_le_m-1'] += (dmax <= m - 1)
            dg = d_gen(m)
            c['dmax_eq_gen'] += (dmax == dg)
            c['dmax_gt_gen'] += (dmax > dg)
            c['dmax_lt_gen'] += (dmax < dg)
            c['ncopy_gt4'] += (m > 4)
            c['dmax_le1_gt4'] += (m > 4 and dmax <= 1)
            c['dG_max_sum'] += max(krG['deg'])
            # numerical check at 50 digits against the copies
            err = mp.mpf(0)
            old = mp.mp.dps
            mp.mp.dps = dps + 10
            try:
                for cp in cps:
                    u = cp['coords']
                    trip = cp['ev']
                    for y in range(m):
                        pred = poly_eval(kr['poly'][y], trip) * u[x0]
                        err = max(err, abs(pred - u[y]) / max(mp.mpf(1), abs(u[y])))
            finally:
                mp.mp.dps = old
            c['err_max'] = max(c['err_max'], float(err))
            c['err_bad'] += (err > mp.mpf(10) ** -(dps - 20))
            # does every coordinate have a cyclic first coordinate?  (x0 = 0 is G3's normalisation point only when its coordinate is nonzero)
            rows.append((mu, m, x0, dmax, max(krG['deg']), d_gen(m), sorted(Counter(kr['deg']).items())))
        c['time'] = time.perf_counter() - t0
        tab('')
        tab('  %s: %d sectors with m >= 2 (sum m = %d): coordinates cyclic for the triple %d of %d (all coordinates cyclic in %d sectors; some cyclic coordinate in %d); cyclic for the triple <=> cyclic for G in %d sectors;' %
            (prod_str(lam, lamp), c['sectors'], c['m_total'], c['cyclic_coords'], c['coords'], c['all_cyclic'], c['some_cyclic'], c['G_cyclic_equal']))
        if c['some_cyclic']:
            tab('     best cyclic coordinate: mean minimal max-degree in the triple %.2f (<= m-1 in %d sectors; = / > / < the generic dimension-counting bound d_gen(m) in %d / %d / %d sectors), mean for G alone %.2f;  exact polynomials vs G3 copies at 50 digits: max relative error %.2e, failures %d  (%.1f s)' %
                (c['dmax_sum'] / c['some_cyclic'], c['dmax_le_m-1'], c['dmax_eq_gen'], c['dmax_gt_gen'], c['dmax_lt_gen'], c['dG_max_sum'] / c['some_cyclic'], c['err_max'], c['err_bad'], c['time']))
        if verbose and len(rows) <= 40:
            tab('     mu            m  x0(best)  max total degree (triple)  max degree (G)  generic bound  degree profile {degree: #coordinates}')
            for r_ in rows:
                if r_[2] == 'NO cyclic coordinate':
                    tab('     %-13s %-2d  NO cyclic coordinate' % (mu_str(r_[0]), r_[1]))
                else:
                    tab('     %-13s %-2d  %-8d  %-26d %-15d %-13d %s' % (mu_str(r_[0]), r_[1], r_[2], r_[3], r_[4], r_[5], r_[6]))
        for k_, v_ in c.items():
            if k_ not in ('err_max', 'time'):
                tot[k_] += v_
            else:
                tot[k_] = max(tot[k_], v_) if k_ == 'err_max' else tot[k_] + v_
    tab('')
    tab('  TOTAL: %d sectors with m >= 2; cyclic coordinates %d of %d; all coordinates cyclic in %d sectors; a cyclic coordinate exists in %d sectors; polynomial fits checked at 50 digits: worst error %.2e, failures %d' %
        (tot['sectors'], tot['cyclic_coords'], tot['coords'], tot['all_cyclic'], tot['some_cyclic'], tot['err_max'], tot['err_bad']))
    tab('         by multiplicity m: number of sectors with (best total degree in the triple, degree for the generic single operator G) = ...  and d_gen(m):')
    for m_ in sorted(BYM):
        tab('            m = %2d  d_gen = %d : %s' % (m_, d_gen(m_), ', '.join('(%d, %d): %d' % (a, b, n) for (a, b), n in sorted(BYM[m_].items()))))
    tab('         degree of the best cyclic coordinate compared with the dimension-counting bound d_gen(m) (smallest d with C(d+3,3) >= m): equal in %d, larger in %d, smaller in %d sectors; sectors with m > 4: %d, of which affine (degree 1) %d' %
        (tot['dmax_eq_gen'], tot['dmax_gt_gen'], tot['dmax_lt_gen'], tot['ncopy_gt4'], tot['dle1_gt4'] if False else tot['dmax_le1_gt4']))


def stage_cycW(products=RULE_PRODUCTS + RULE_HELD, verbose=True):
    """Q5(c) on the WHOLE weight space: which lattice points x of R_nu are cyclic coordinate vectors for the triple on W (every eigenvector of W, physical or spurious, has a nonzero component at x)?"""
    hdr('Q5(c)  CYCLIC LATTICE POINTS ON THE WHOLE WEIGHT SPACE W = V_lam[mu - lam\'] (exact graded Krylov spaces of the coordinate vectors under the triple)')
    tab('  For every weight space with w >= 2 (all irreps of the product, identity translate): the number of lattice points that are cyclic; whether the first / last GT-index point is cyclic; the best graded')
    tab('  degree (max over lattice points y of the minimal total degree of e_y^T in the Krylov space of e_x^T, minimised over the cyclic x) against d_gen(w).')
    tot = Counter()
    for lam, lamp in products:
        R = gr.GLNRep.of(lam)
        c = Counter()
        ops = {n: so.reduced_op(n, lamp) for n in so.TRIPLE}
        t0 = time.perf_counter()
        for mu, m in sector_list(lam, lamp, 'rs'):
            nu = tuple(mu[i] - lamp[i] for i in range(4))
            W = R.wsp.get(nu)
            if W is None or len(W) < 2:
                continue
            w = len(W)
            mats_d = so.weight_space_mats(R, ops, W)
            mats = [mats_d[n] for n in so.TRIPLE]
            cyc = []
            best = None
            for x in range(w):
                kr = left_krylov(mats, x)
                cyc.append(kr['cyclic'])
                if kr['cyclic']:
                    d = max(kr['deg'])
                    best = d if best is None else min(best, d)
            c['spaces'] += 1
            c['points'] += w
            c['cyclic_points'] += sum(cyc)
            c['some'] += any(cyc)
            c['all'] += all(cyc)
            c['first'] += cyc[0]
            c['last'] += cyc[-1]
            if best is not None:
                c['best_eq_gen'] += (best == d_gen(w))
                c['best_gt_gen'] += (best > d_gen(w))
            c['m<w'] += (m < w)
            c['m<w_some'] += (m < w and any(cyc))
        tab('  %-22s %3d spaces (w >= 2, %d with m < w): cyclic lattice points %d of %d; a cyclic point exists in %d, all points cyclic in %d; first GT-index point cyclic in %d, last in %d; best degree = d_gen(w) in %d, > in %d  (%.1f s)' %
            (prod_str(lam, lamp), c['spaces'], c['m<w'], c['cyclic_points'], c['points'], c['some'], c['all'], c['first'], c['last'], c['best_eq_gen'], c['best_gt_gen'], time.perf_counter() - t0))
        for k_, v_ in c.items():
            tot[k_] += v_
    tab('  TOTAL: %d spaces: cyclic lattice points %d of %d (%.1f %%); a cyclic point exists in %d of %d, in %d of the %d spaces with m < w; first point cyclic in %d, last in %d' %
        (tot['spaces'], tot['cyclic_points'], tot['points'], 100.0 * tot['cyclic_points'] / max(1, tot['points']), tot['some'], tot['spaces'], tot['m<w_some'], tot['m<w'], tot['first'], tot['last']))


# ======================================================================================== 11. the Weyl stage driver (Q3, Q6) incl. gl(3) sanity check and sensitivity tests
WEYL_MAIN = [((2, 1, 0, 0), (2, 1, 0, 0)), ((2, 1, 1, 0), (2, 1, 1, 0)), ((3, 2, 1, 0), (3, 2, 1, 0)), ((4, 2, 0, 0), (4, 2, 0, 0)), ((4, 2, 1, 0), (4, 2, 1, 0)), ((6, 4, 2, 0), (6, 4, 2, 0))]
WEYL_EXTRA = [((3, 2, 1, 0), (4, 2, 1, 0)), ((4, 2, 1, 0), (3, 2, 1, 0)), ((4, 2, 2, 0), (4, 2, 2, 0)), ((4, 3, 1, 0), (4, 3, 1, 0)), ((3, 2, 0, 0), (2, 1, 1, 0)), ((2, 1, 1, 0), (2, 2, 1, 0)),
              ((5, 3, 1, 0), (5, 3, 1, 0)), ((5, 3, 2, 0), (5, 3, 2, 0)), ((6, 3, 0, 0), (6, 3, 0, 0)), ((5, 2, 1, 1), (3, 2, 2, 1)),
              ((5, 4, 2, 1), (5, 4, 2, 1)), ((6, 3, 2, 1), (6, 3, 2, 1)), ((8, 4, 0, 0), (8, 4, 0, 0)), ((4, 4, 2, 0), (4, 3, 2, 1)), ((7, 4, 1, 0), (7, 4, 1, 0)), ((6, 4, 2, 0), (5, 3, 1, 0))]
WEYL_GL3 = [((4, 2, 0), (4, 2, 0)), ((6, 3, 0), (6, 3, 0)), ((8, 4, 0), (8, 4, 0)), ((6, 4, 2), (6, 4, 2)), ((14, 7, 0), (14, 7, 0)), ((4, 2, 0), (3, 1, 0)), ((5, 2, 1), (4, 3, 0))]


def gl3_jacobi_crosscheck(lams=((4, 2, 0), (6, 3, 0), (6, 4, 2), (8, 4, 0))):
    """the N = 3 family used here (G3's reduction of t S_{3,2}) equals the Jacobi operator V_red of gaudin_sl3/jacobi_cg plus the sector constant eta_1/2, eta_1 = |Lam|^2 + cas2(Lam) - cas2(mu)/2:
    compares the 50-digit physical spectra of ident_sector with jacobi_cg.physical_V on every sector"""
    if so.GL3DIR not in sys.path:
        sys.path.insert(0, so.GL3DIR)
    import jacobi_cg as jc
    out = []
    old = mp.mp.dps
    mp.mp.dps = 60
    try:
        for lam in lams:
            nsec = nok = 0
            worst = mp.mpf(0)
            for mu, m in sector_list(lam, lam, 'rs'):
                sec = ident_sector(lam, lam, mu, FAM[3])
                ours = sorted(mp.re(t[0]) for t, p_ in zip(sec['js']['trip'], sec['phys']) if p_)
                theirs = jc.physical_V(lam, mu, dps=50)
                cas = lambda nu: sum(nu[i] * (nu[i] + 4 - 2 * (i + 1)) for i in range(3))
                eta1 = sum(lam) ** 2 + cas(lam) - mp.mpf(cas(mu)) / 2
                dev = max(abs(a - b - eta1 / 2) for a, b in zip(ours, theirs))
                nsec += 1
                worst = max(worst, dev)
                nok += dev < mp.mpf(10) ** -30
            out.append((lam, nsec, nok, worst))
    finally:
        mp.mp.dps = old
    return out


def stage_weyl(main=WEYL_MAIN, extra=WEYL_EXTRA, gl3=WEYL_GL3, sens=True, dps=50):
    hdr('Q3 + Q6  THE WEYL / BGG CHARACTERISATION OF THE PHYSICAL SPECTRUM')
    tab('  Conjecture: for every joint triple tau,  sum_{w in S_4} sgn(w) mult_tau( V_red(w.lam\') on V_lam[mu - w.lam\'] ) = [tau in Spec(M_mu)],  w.lam\' = w(lam\' + rho) - rho,  rho = (3,2,1,0).')
    tab('  Test: for every sector mu, all 24 Weyl translates; the operators V_red,k at the (generally non-dominant) weight w.lam\' on the WHOLE weight space V_lam[mu - w.lam\'] (G3 reduced operators, universal')
    tab('  polynomials evaluated at w.lam\'), joint spectra at 50 digits (mpmath eigen-decomposition of the generic combination V3 + 7/19 V4a + 11/13 V4b, every operator conjugated to the eigenbasis),')
    tab('  triples matched across the translate spaces (tolerance 1e-30), signed counts compared with the indicator of the m physical triples (joint eigenvalues of the reduced operators on M_mu, PRV basis).')
    tab('  Sectors: kind "all" = every dominant mu with |mu| = |lam| + |lam\'| for which some translate space is non-empty (the irreps of Racah-Speiser AND formal/multiplicity-0 candidates).')
    results = {}
    rows = []
    for lam, lamp in main + extra:
        small = (lam, lamp) in main[:5]
        t0 = time.perf_counter()
        res, S = weyl_product(lam, lamp, kind='all', listing=small, maxlist=40, dps=dps)
        results[(lam, lamp)] = (res, S)
        bgg_signature_table(res, prod_str(lam, lamp))
        # multi-string triples: which translates contain them
        ms = []
        for r in res:
            for g in r['groups']:
                if len(g['lens']) >= 3:
                    ms.append((r['mu'], tuple(sorted(p_name(r['spaces'][si]['p']) for si, k in g['members'])), g['signed'], g['phys']))
        if ms:
            tab('     triples occurring in three or more translate spaces (sector, permutations of the translates containing them, signed count, physical): %s' % ms[:6])
        rows.append((lam, lamp, S))
    tab('')
    tab('  MAIN TABLE (one row per product; sectors = all dominant mu with a non-empty translate space):')
    tab('  %-24s %7s %6s %8s %9s | %-10s %-10s | %-12s | %-34s | %s | %s' % ('product', 'sectors', 'sum m', 'spaces', 'triples', 'signed=ind.', 'failures', 'lattice cnt', 'one-step: S<=T / S=T / multi (of triples)', 'sectors needing depth L = 0 / 1 / 2 / >2', 'longest translate'))
    for lam, lamp, S in rows:
        Lc = [S['Lmin_%s' % L] for L in (0, 1, 2)] + [sum(v for k, v in S.items() if k.startswith('Lmin_') and k.split('_')[1] not in ('0', '1', '2'))]
        tab('  %-24s %7d %6d %8d %9d | %-10s %-10d | %-12s | %-34s | %s' % (prod_str(lam, lamp), S['sectors'], S['m_total'], S['spaces'], S['groups'], '%d/%d' % (S['ok'], S['sectors']), S['fail_sectors'],
            '%d/%d' % (S['mult_ok'], S['sectors']), '%d / %d / %d (of %d)' % (S['os_S_sub_T'], S['os_S_eq_T'], S['os_multi'], S['os_tau']), ' / '.join(str(x) for x in Lc) + '   | %d' % S['maxlen']))
    tot = Counter()
    for lam, lamp, S in rows:
        for k_ in ('sectors', 'm_total', 'spaces', 'groups', 'ok', 'fail_sectors', 'fail_groups', 'mult_ok', 'os_tau', 'os_S_sub_T', 'os_S_eq_T', 'os_multi', 'os_spurious', 'os_spur_single', 'os_spur_multi'):
            tot[k_] += S[k_]
    tab('  %-24s %7d %6d %8d %9d | %-10s %-10d | %-12s | %-34s |' % ('TOTAL', tot['sectors'], tot['m_total'], tot['spaces'], tot['groups'], '%d/%d' % (tot['ok'], tot['sectors']), tot['fail_sectors'], '%d/%d' % (tot['mult_ok'], tot['sectors']),
        '%d / %d / %d (of %d)' % (tot['os_S_sub_T'], tot['os_S_eq_T'], tot['os_multi'], tot['os_tau'])))
    tab('  spurious triples (w > m): %d; one violated string %d, two or more %d (these need the length-2 translate of the inclusion-exclusion)' % (tot['os_spurious'], tot['os_spur_single'], tot['os_spur_multi']))
    # sensitivity
    if sens:
        tab('')
        tab('  Sensitivity (the test must be able to fail): the same signed count with deliberately wrong rules, per product (sectors where the identity "signed count = indicator" still holds):')
        for lam, lamp in [((3, 2, 1, 0), (3, 2, 1, 0)), ((4, 2, 1, 0), (4, 2, 1, 0)), ((6, 4, 2, 0), (6, 4, 2, 0))]:
            out = []
            for var in ('rule', 'nosign', 'plain'):
                res, S = weyl_product(lam, lamp, kind='rs', onestep=False, listing=False, verbose=False, variant=var, dps=dps)
                out.append('%s: %d/%d' % ({'rule': 'dot action + signs (the conjecture)', 'nosign': 'all signs +1', 'plain': 'plain Weyl action w(lam\') without the rho-shift'}[var], S['ok'], S['sectors']))
            tab('     %-14s %s' % (prod_str(lam, lamp), ';  '.join(out)))
    # gl(3) sanity check
    tab('')
    tab('  gl(3) sanity check of the same statement (S_3, rho = (2,1,0), the single operator V = t S_{3,2} reduced by the same engine, N = 3):')
    g3rows = []
    for lam, lamp in gl3:
        res, S = weyl_product(lam, lamp, kind='all', fam=FAM[3], listing=False, verbose=False, dps=dps)
        # top-m rule for the physical eigenvalue(s) of V on the whole weight space
        nd = Counter()
        for mu, m in sector_list(lam, lamp, 'rs'):
            sec = ident_sector(lam, lamp, mu, FAM[3], dps)
            if sec is None or sec['w'] == sec['m']:
                continue
            vals = [mp.re(t[0]) for t in sec['js']['trip']]
            nd[(tuple(sorted(sec['active'])), sep_direction(vals, sec['phys']))] += 1
        g3rows.append((lam, lamp, S, nd))
        tab('     %-18s sectors %4d (sum m %4d, %4d spaces, %4d triples): signed count = indicator %d/%d, failures %d; depth needed L=0/1/2/3: %s; one-step S<=T %d/%d, S=T %d/%d, multi %d; lattice count %d/%d; '
            'ordering of the physical eigenvalues of V on W, sectors with w > m, by (active string, direction): %s' % (prod_str(lam, lamp), S['sectors'], S['m_total'], S['spaces'], S['groups'], S['ok'], S['sectors'], S['fail_sectors'],
            '/'.join(str(S['Lmin_%s' % L]) for L in (0, 1, 2, 3)), S['os_S_sub_T'], S['os_tau'], S['os_S_eq_T'], S['os_tau'], S['os_multi'], S['mult_ok'], S['sectors'], dict(nd)))
    for lam, nsec, nok, worst in gl3_jacobi_crosscheck():
        tab('     cross-check with gaudin_sl3/jacobi_cg: %s: the 50-digit physical spectra of the N = 3 family used here equal jacobi_cg.physical_V + eta_1/2 on %d/%d sectors (worst deviation %s)' % (lam_str(lam), nok, nsec, mp.nstr(worst, 3)))
    return results, g3rows


# ======================================================================================== quick driver for development
def main(argv):
    stage = argv[1] if len(argv) > 1 else 'all'
    if stage == 'stencil':
        stage_stencil()
        return
    if stage == 'symq':
        stage_symq()
        return
    if stage == 'univ':
        stage_univ()
        return
    if stage == 'symbolic':
        stage_symbolic(tuple(argv[2].split(',')) if len(argv) > 2 else ('V3', 'V4a', 'V4b'))
        return
    if stage == 'rule':
        stage_rule(verbose=False)
        return
    if stage == 'chain':
        stage_chain()
        return
    if stage == 'poly':
        stage_poly()
        return
    if stage == 'cycW':
        stage_cycW()
        return
    if stage == 'rule_dev':
        stage_rule(products=[((2, 1, 0, 0), (2, 1, 0, 0)), ((3, 2, 1, 0), (3, 2, 1, 0))], held=[((4, 2, 1, 0), (4, 2, 1, 0))])
        return
    if stage == 'chain_dev':
        stage_chain(products=[((3, 2, 1, 0), (3, 2, 1, 0)), ((4, 2, 1, 0), (4, 2, 1, 0))])
        return
    if stage == 'poly_dev':
        stage_poly(products=[((3, 2, 1, 0), (3, 2, 1, 0)), ((4, 2, 1, 0), (4, 2, 1, 0))])
        return
    if stage == 'sym_probe':
        stage_symbolic_probe(tuple(argv[2].split(',')) if len(argv) > 2 else ('V3',))
        return
    if stage == 'weyl':
        stage_weyl()
        return
    if stage == 'run':                       # python structure_gl4.py run stage1,stage2,...  -> structure_gl4_tables_<stage>.part for every stage
        fns = {'stencil': stage_stencil, 'symbolic': stage_symbolic, 'symq': stage_symq, 'univ': stage_univ, 'weyl': stage_weyl, 'rule': lambda: stage_rule(verbose=True), 'chain': stage_chain, 'poly': stage_poly, 'cycW': stage_cycW}
        for nm in argv[2].split(','):
            t0 = time.perf_counter()
            fns[nm]()
            tab('  [stage %s finished in %.1f s]' % (nm, time.perf_counter() - t0))
            flush_tables(os.path.join(HERE, 'structure_gl4_tables_%s.part' % nm), 'w')
        return
    if stage == 'merge':
        title = 'structure_gl4_tables.txt: every table of structure_gl4.py (WP-G4), %s' % DATE
        with open(TABLES_PATH, 'w', encoding='utf-8') as f:
            f.write(title + chr(10))
            for nm in ('stencil', 'symbolic', 'symq', 'univ', 'weyl', 'rule', 'chain', 'poly', 'cycW'):
                pth = os.path.join(HERE, 'structure_gl4_tables_%s.part' % nm)
                if os.path.exists(pth):
                    f.write(open(pth, encoding='utf-8').read())
                    os.remove(pth)
        return
    if stage == 'all':
        title = 'structure_gl4_tables.txt: every table of structure_gl4.py (WP-G4), %s' % DATE
        first = True
        for nm, fn in (('stencil', stage_stencil), ('symbolic', stage_symbolic), ('weyl', stage_weyl), ('rule', lambda: stage_rule(verbose=True)), ('chain', stage_chain), ('poly', stage_poly), ('cycW', stage_cycW)):
            t0 = time.perf_counter()
            fn()
            tab('  [stage %s finished in %.1f s]' % (nm, time.perf_counter() - t0))
            flush_tables(TABLES_PATH, 'w' if first else 'a', title if first else None)
            first = False
        return
    if stage == 'weyl_dev':
        lam = tuple(int(x) for x in argv[2].split(','))
        lamp = tuple(int(x) for x in argv[3].split(',')) if len(argv) > 3 else lam
        kind = argv[4] if len(argv) > 4 else 'rs'
        res, S = weyl_product(lam, lamp, kind=kind)
        bgg_signature_table(res, prod_str(lam, lamp))
        return
    raise SystemExit('unknown stage %s' % stage)


if __name__ == '__main__':
    main(sys.argv)
