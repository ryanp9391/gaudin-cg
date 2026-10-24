"""relabelled_gl4.py -- WP-G8c: the 24 relabelled Gelfand-Tsetlin chains of gl(4); the commuting triple (V3, V4a, V4b) and the three PRV strings in every sigma-GT basis.

CONVENTION (fixed and checked, section C0).  For sigma in S_4 (tuple sigma = (s1, s2, s3, s4), 1-based values) the sigma-chain is gl(4) > gl({s1,s2,s3}) > gl({s1,s2}) > gl({s1}) and
the sigma-GT basis b_p of V_lam is the GT basis adapted to it.  By definition E_{s_a s_b} acts on b_p by the STANDARD GT formulas of E_ab.  Hence, with tau = sigma^-1,
      matrix of pi(E_cd) in the sigma-basis  =  standard-basis matrix of E_{tau(c) tau(d)},
      pi(X) for a word X  ->  the word relabelled E_cd -> E_{tau(c) tau(d)}, evaluated in the plain GLNRep.of(lam),
      the sigma-basis vector b_p has E_cc eigenvalue nu_c  <=>  the standard pattern has weight  wt_a = nu_{sigma(a)}   (so it lives in R.wsp[(nu_{s1},..,nu_{s4})]),
      the string E_{i,i+1}^{k_i}  ->  E_{tau(i) tau(i+1)}^{k_i}  (acts on W^sigma = R.wsp[nu^sigma]).
(The brief's phrase 'E_ab -> E_{sigma(a)sigma(b)}' is the same statement for the representation pi o phi_sigma; for involutions sigma = tau.)  A string E_{i,i+1} is a coordinate cut in the sigma-basis
exactly when {tau(i), tau(i+1)} = {1,2}  <=>  {i, i+1} = {s1, s2}.  `conv='dir'` uses tau = sigma instead (the WRONG convention for non-involutions; used only by the check).

Public API
  sigma_mats(lam, lamp, mu, sigma, conv='inv', names=TRIPLE)  dict: W (global GT indices of W^sigma), pats, coords (m12,m13,m23 of the sigma-patterns), mats {name: w x w exact}, strings [...], A (active strings)
  restrict(sm, ...), kernel info, graph statistics, the CLI stages  check | bandwidth | cuts | small | parity | all   (python -X utf8 relabelled_gl4.py <stage>)
"""
import os, sys, time, itertools, pickle
from fractions import Fraction as Fr
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mpmath as mp
import glN_reps as gr
import slice_engine_gl4 as se
import slice_ops_gl4 as so
import structure_gl4 as sg

TRIPLE = so.TRIPLE
SIGMAS = list(itertools.permutations((1, 2, 3, 4)))
ID = (1, 2, 3, 4)
W0 = (4, 3, 2, 1)
TABLES_PATH = os.path.join(HERE, 'relabelled_gl4_tables.txt')
CACHE_PATH = os.path.join(HERE, 'relabelled_gl4_records.pkl')
_TAB = []


def tab(s=''):
    _TAB.append(s)
    print(s, flush=True)


def hdr(s):
    tab('')
    tab('=' * 140)
    tab(s)
    tab('=' * 140)


def sstr(s):
    return ''.join(str(x) for x in s)


def inv(s):
    t = [0] * 4
    for a, x in enumerate(s):
        t[x - 1] = a + 1
    return tuple(t)


def pstr(lam, lamp):
    return so.prod_str(lam, lamp)


# ======================================================================================== 1. relabelled operators
_OPS = {}


def relabelled_op(name, lamp, tau):
    """ReducedOp of `name` at lam' with every generator E_cd relabelled to E_{tau(c) tau(d)} (tau: tuple, tau[c-1] = tau(c))"""
    key = (name, tuple(lamp), tuple(tau))
    if key not in _OPS:
        base = so.reduced_op(name, lamp)
        if tuple(tau) == ID:
            _OPS[key] = base
        else:
            red = {tuple((tau[a - 1], tau[b - 1]) for (a, b) in w): c for w, c in base.red.items()}
            _OPS[key] = so.ReducedOp(red, name, tuple(lamp))
    return _OPS[key]


def conv_maps(sigma, conv='inv'):
    """(tau, nu-permutation): tau relabels the generators, nu^sigma_a = nu_{pi(a)} the weight"""
    sigma = tuple(sigma)
    if conv == 'inv':
        return inv(sigma), sigma
    if conv == 'dir':                                   # consistent pair for the chain of sigma^-1 (NOT wrong: = 'inv' at sigma^-1)
        return sigma, inv(sigma)
    if conv == 'badw':                                  # WRONG: generators relabelled by tau = sigma^-1 but the weight permuted by sigma^-1 too (should be sigma)
        return inv(sigma), inv(sigma)
    if conv == 'badt':                                  # WRONG: generators relabelled by sigma, weight permuted by sigma (should be sigma^-1 for the generators... pair mismatch)
        return sigma, sigma
    raise ValueError(conv)


def string_info(R, tau, lamp, nup, W, void_ok=True):
    """per string i = 1..3 of the sigma-chain: k_i, generator (tau(i), tau(i+1)), void?, exact kernel rows (list of rows over the columns of W), zero columns, kernel dimension, coordinate?"""
    out = []
    w = len(W)
    for i in range(1, 4):
        k = lamp[i - 1] - lamp[i] + 1
        g = (tau[i - 1], tau[i])
        tw = gr.wadd(nup, gr.dij(4, g[0], g[1]), k)
        if tw not in R.wsp:
            out.append(dict(i=i, k=k, gen=g, void=True, rows=[], zero=set(range(w)), nullity=w, coord=True))
            continue
        imgs = se.power_images(R, g, k, W)
        byrow = defaultdict(lambda: [Fr(0)] * w)
        for a, vec in enumerate(imgs):
            for r, x in vec.items():
                byrow[r][a] = x
        rows = list(byrow.values())
        basis, free = se._rref_nullspace(rows, w)
        zero = {j for j in range(w) if all(r[j] == 0 for r in rows)}
        out.append(dict(i=i, k=k, gen=g, void=False, rows=rows, zero=zero, nullity=len(basis), coord=(len(basis) == len(zero))))
    return out


def sigma_mats(lam, lamp, mu, sigma, conv='inv', names=TRIPLE, with_strings=True):
    """the sigma-basis data of the sector (lam, lamp, mu): W^sigma, the sigma-patterns and coordinates, exact matrices of the triple, and the strings"""
    lam, lamp, mu = tuple(lam), tuple(lamp), tuple(mu)
    tau, perm = conv_maps(sigma, conv)
    R = gr.GLNRep.of(lam)
    nu = tuple(mu[a] - lamp[a] for a in range(4))
    nup = tuple(nu[perm[a] - 1] for a in range(4))
    W = R.wsp.get(nup)
    if W is None:
        return None
    W = [int(c) for c in W]
    ops = {n: relabelled_op(n, lamp, tau) for n in names}
    mats = so.weight_space_mats(R, ops, W)
    pats = [R.pats[c] for c in W]
    d = dict(lam=lam, lamp=lamp, mu=mu, sigma=tuple(sigma), tau=tau, R=R, nu=nu, nup=nup, W=W, w=len(W), pats=pats, coords=[sg.coords(p) for p in pats], mats=mats)
    if with_strings:
        d['strings'] = string_info(R, tau, lamp, nup, W)
        d['A'] = active_strings(R, lamp, nu)
    return d


def active_strings(R, lamp, nu):
    """identity-frame active strings A = {i : E_{i,i+1}^{k_i} maps W = V_lam[nu] to an existing weight space} (Weyl invariant notion: the same set in every sigma-frame, with nu^sigma)"""
    A = set()
    for i in range(1, 4):
        k = lamp[i - 1] - lamp[i] + 1
        if gr.wadd(nu, gr.dij(4, i, i + 1), k) in R.wsp:
            A.add(i)
    return frozenset(A)


# ======================================================================================== 2. exact linear algebra helpers
def mat_vec(M, v):
    return [sum((M[i][j] * v[j] for j in range(len(v)) if v[j] != 0), Fr(0)) for i in range(len(M))]


def charpoly(M):
    """exact characteristic polynomial coefficients [1, c1, ..., cn] of det(x - M) (Faddeev-LeVerrier, Fractions)"""
    n = len(M)
    if n == 0:
        return [Fr(1)]
    cs = [Fr(1)]
    Mk = [[Fr(i == j) for j in range(n)] for i in range(n)]
    for k in range(1, n + 1):
        AM = [[sum((M[i][t] * Mk[t][j] for t in range(n) if M[i][t] != 0), Fr(0)) for j in range(n)] for i in range(n)]
        ck = -sum((AM[i][i] for i in range(n)), Fr(0)) / k
        cs.append(ck)
        Mk = [[AM[i][j] + (ck if i == j else 0) for j in range(n)] for i in range(n)]
    return cs


def lincomb(mats, combo):
    n = len(mats[0])
    return [[sum((Fr(c) * A[i][j] for c, A in zip(combo, mats)), Fr(0)) for j in range(n)] for i in range(n)]


def submatrix(M, idx):
    return [[M[a][b] for b in idx] for a in idx]


def kernel_basis(sm, which=None):
    """exact basis (list of w-vectors, RREF free-position form) and free positions of M^sigma = W^sigma intersected with the kernels of the non-void strings in `which` (default all)"""
    rows = []
    for s in sm['strings']:
        if which is None or s['i'] in which:
            rows.extend(s['rows'])
    return se._rref_nullspace(rows, sm['w']) if rows else ([[Fr(int(a == b)) for b in range(sm['w'])] for a in range(sm['w'])], list(range(sm['w'])))


def restrict_to(M, basis, free):
    """matrix of M on span(basis) (basis vectors have a 1 at their own free position and 0 at the other free positions); raises ValueError if the span is not M-invariant"""
    m = len(basis)
    w = len(M)
    cols = []
    for k in range(m):
        img = mat_vec(M, basis[k])
        co = [img[f] for f in free]
        rec = [Fr(0)] * w
        for kk, c in enumerate(co):
            if c != 0:
                for i, x in enumerate(basis[kk]):
                    if x != 0:
                        rec[i] += c * x
        if rec != img:
            raise ValueError('span not invariant')
        cols.append(co)
    return [[cols[j][i] for j in range(m)] for i in range(m)]


def graph(M):
    return sg.graph_stats(M) if len(M) >= 2 else dict(n=len(M), edges=0, ncomp=len(M), maxdeg=0, paths=True, tree=True, bw=0, asym=0)


def union_matrix(mats):
    w = len(mats['V3'])
    return [[Fr(1) if (mats['V3'][a][b] != 0 or mats['V4a'][a][b] != 0 or mats['V4b'][a][b] != 0) else Fr(0) for b in range(w)] for a in range(w)]


def stencil_of(sm, name='V3'):
    M = sm['mats'][name]
    C = sm['coords']
    S = set()
    for a in range(sm['w']):
        for b in range(sm['w']):
            if a != b and M[a][b] != 0:
                S.add(tuple(C[a][t] - C[b][t] for t in range(3)))
    return S


# ======================================================================================== 3. C0: convention check
def check_sector(lam, lamp, mu, sigmas=SIGMAS, conv='inv', verbose=False):
    """exact: charpolys of V3, V4a, V4b and of the generic combination on M^sigma equal those on M_mu (identity chain, G3's PRV matrices); returns dict sigma -> status"""
    prv = se.prv_space(lam, lamp, mu)
    R = gr.GLNRep.of(lam)
    ref = {n: so.vred_matrix(R, so.reduced_op(n, lamp), prv) for n in TRIPLE}
    gen = lincomb([ref[n] for n in TRIPLE], so.GENERIC)
    refcp = {n: charpoly(ref[n]) for n in TRIPLE}
    refcp['G'] = charpoly(gen)
    res = {}
    for s in sigmas:
        try:
            sm = sigma_mats(lam, lamp, mu, s, conv)
            if sm is None:
                res[s] = 'no weight space'
                continue
            basis, free = kernel_basis(sm)
            if len(basis) != prv.m:
                res[s] = 'dim M^sigma = %d != m = %d' % (len(basis), prv.m)
                continue
            rm = {n: restrict_to(sm['mats'][n], basis, free) for n in TRIPLE}
            cp = {n: charpoly(rm[n]) for n in TRIPLE}
            cp['G'] = charpoly(lincomb([rm[n] for n in TRIPLE], so.GENERIC))
            # commutation of the restricted triple (exact)
            comm = all(so.mat_mul(rm[a], rm[b]) == so.mat_mul(rm[b], rm[a]) for a, b in (('V3', 'V4a'), ('V3', 'V4b'), ('V4a', 'V4b')))
            res[s] = 'ok' if (cp == refcp and comm) else 'charpoly/commutation mismatch'
        except ValueError as e:
            res[s] = 'M^sigma not invariant under the triple'
    return res, prv.m


def stage_check(verbose=True):
    hdr('C0  CONVENTION CHECK: joint spectrum of the relabelled triple on M^sigma = Spec(M_mu)  (exact charpolys of V3, V4a, V4b and the generic combination + exact commutation)')
    lam = lamp = (3, 2, 1, 0)
    tot = Counter()
    wrong = Counter()
    nsec = 0
    for mu, m in sg.sector_list(lam, lamp, 'rs'):
        nsec += 1
        res, mm = check_sector(lam, lamp, mu)
        for s, st in res.items():
            tot[st] += 1
        res2, _ = check_sector(lam, lamp, mu, conv='badw')
        for s, st in res2.items():
            wrong[(s == inv(s), st)] += 1
    tab('  [3,2,1,0]^2: %d sectors x 24 sigma, convention tau = sigma^-1 (E_cd -> E_{tau c tau d}, weight nu_{sigma(a)}):' % nsec)
    for st, c in tot.items():
        tab('      %-48s %5d' % (st, c))
    tab('  the WRONG convention (tau = sigma^-1 but weight nu_{sigma^-1(a)} instead of nu_{sigma(a)}), split by involution (sigma = sigma^-1) / non-involution:')
    for (invol, st), c in sorted(wrong.items(), key=lambda t: (not t[0][0], t[0][1])):
        tab('      %-14s %-48s %5d' % ('involution' if invol else 'non-involution', st, c))
    return tot, wrong


# ======================================================================================== 4. products, records, bandwidth
P_SMALL = [((2, 1, 0, 0), (2, 1, 0, 0)), ((2, 1, 1, 0), (2, 1, 1, 0)), ((3, 2, 1, 0), (3, 2, 1, 0)), ((4, 2, 0, 0), (4, 2, 0, 0)), ((4, 2, 1, 0), (4, 2, 1, 0)),
           ((3, 2, 1, 0), (4, 2, 1, 0)), ((4, 2, 1, 0), (3, 2, 1, 0)), ((4, 2, 2, 0), (4, 2, 2, 0)), ((4, 3, 1, 0), (4, 3, 1, 0))]
P_BIG = [((6, 4, 2, 0), (6, 4, 2, 0)), ((5, 3, 1, 0), (5, 3, 1, 0)), ((5, 3, 2, 1), (5, 3, 2, 1))]


def sector_record(lam, lamp, mu, m, sigmas=SIGMAS):
    """one record per weight space (w >= 2) and per sigma: graph statistics of V3, of the triple's union graph, stencil, string classification, coordinate subspace data"""
    R = gr.GLNRep.of(lam)
    nu = tuple(mu[a] - lamp[a] for a in range(4))
    W = R.wsp.get(nu)
    if W is None or len(W) < 2:
        return None
    rec = dict(lam=lam, lamp=lamp, mu=mu, m=m, w=len(W), A=active_strings(R, lamp, nu), per={})
    for s in sigmas:
        sm = sigma_mats(lam, lamp, mu, s)
        g3 = graph(sm['mats']['V3'])
        gu = graph(union_matrix(sm['mats']))
        st = sm['strings']
        info = dict(g3=g3, gu=gu, stencil=frozenset(stencil_of(sm)), coord=tuple(x['coord'] for x in st), void=tuple(x['void'] for x in st),
                    nullity=tuple(x['nullity'] for x in st), zero=tuple(len(x['zero']) for x in st))
        # M^sigma: coordinate iff every string kernel is coordinate (a coordinate kernel intersection is a coordinate subspace); then the triple restricted to it
        allc = all(x['coord'] for x in st)
        info['Mcoord'] = allc
        info['coords'] = tuple(sm['coords'])
        info['zero'] = tuple(tuple(sorted(x['zero'])) for x in st)
        # the coordinate sub-lattice cut out by the strings that ARE coordinate in this sigma (non-void ones): invariant under the triple (BGG), contains M_mu
        cut = [x for x in st if x['coord'] and not x['void']]
        info['cuts'] = tuple(x['i'] for x in cut)
        if cut:
            cidx = sorted(set.intersection(*[x['zero'] for x in cut]))
            info['cdim'] = len(cidx)
            sub = {n: submatrix(sm['mats'][n], cidx) for n in TRIPLE}
            cin = set(cidx)
            info['cinv'] = all(sm['mats'][n][a][b] == 0 for n in TRIPLE for b in cidx for a in range(sm['w']) if a not in cin)
            info['cg3'] = graph(sub['V3'])
        if allc:
            idx = sorted(set.intersection(*[x['zero'] for x in st]))
            info['Mdim'] = len(idx)
            info['Midx'] = tuple(idx)
            sub = {n: submatrix(sm['mats'][n], idx) for n in TRIPLE}
            inside = set(idx)
            ok = all(sm['mats'][n][a][b] == 0 for n in TRIPLE for b in idx for a in range(sm['w']) if a not in inside)
            info['Minv'] = ok
            info['gM3'] = graph(sub['V3'])
            info['gMu'] = graph(union_matrix(sub)) if len(idx) >= 1 else None
            info['Mcoords'] = tuple(sm['coords'][a] for a in idx)
        rec['per'][s] = info
    return rec


def all_sectors(products):
    for lam, lamp in products:
        for mu, m in sg.sector_list(lam, lamp, 'rs'):
            yield lam, lamp, mu, m


def _cache_file(lam, lamp):
    return os.path.join(HERE, 'relabelled_gl4_rec_%s_%s.pkl' % (''.join(str(x) for x in lam), ''.join(str(x) for x in lamp)))


def compute_records(products, cache=True, verbose=True):
    out = []
    for lam, lamp in products:
        f = _cache_file(lam, lamp)
        if cache and os.path.exists(f):
            try:
                out.extend(pickle.load(open(f, 'rb')))
                continue
            except Exception:
                pass
        t0 = time.perf_counter()
        recs = []
        for lam_, lamp_, mu, m in all_sectors([(lam, lamp)]):
            r = sector_record(lam_, lamp_, mu, m)
            if r is not None:
                recs.append(r)
        out.extend(recs)
        if verbose:
            print('  records %s: %d weight spaces (%.1f s)' % (pstr(lam, lamp), len(recs), time.perf_counter() - t0), flush=True)
        if cache:
            pickle.dump(recs, open(f, 'wb'))
    return out


def A_str(A):
    return '{' + ','.join(str(i) for i in sorted(A)) + '}' if A else '{}'


def kind_of(r):
    return 'm=w' if r['m'] == r['w'] else 'm<w'


def stage_bandwidth(recs, label=''):
    hdr('C1  BANDWIDTH AND CHAIN STRUCTURE PER SIGMA  %s' % label)
    kinds = ('m=w', 'm<w')
    n = Counter(kind_of(r) for r in recs)
    tab('  weight spaces (w >= 2): m = w: %d, m < w: %d, total %d' % (n['m=w'], n['m<w'], len(recs)))
    tab('')
    tab('  A. per sigma: number of weight spaces where the sparsity graph of V3 is a PATH (union graph of the triple in brackets), mean RCM bandwidth of V3, max bandwidth; connected counts')
    tab('     %-6s %-4s | %-12s %-12s | %-12s %-12s | %-9s %-9s | %-9s' % ('sigma', 'chn', 'paths m=w', '(triple)', 'paths m<w', '(triple)', 'mean bw', 'max bw', 'connected'))
    rows = []
    for s in SIGMAS:
        c = Counter()
        for r in recs:
            k = kind_of(r)
            p = r['per'][s]
            c['p3', k] += p['g3']['paths']
            c['pu', k] += p['gu']['paths']
            c['bw'] += p['g3']['bw']
            c['bwmax'] = max(c['bwmax'], p['g3']['bw'])
            c['conn'] += (p['g3']['ncomp'] == 1)
        rows.append((s, c))
    for s, c in rows:
        chain = '{%d%d%d}{%d%d}{%d}' % (s[0], s[1], s[2], s[0], s[1], s[0])
        tab('     %-6s %-4s | %5d/%-5d %5d/%-5d | %5d/%-5d %5d/%-5d | %9.2f %9d | %9d' % (sstr(s), 'id' if s == ID else ('w0' if s == W0 else ''), c['p3', 'm=w'], n['m=w'], c['pu', 'm=w'], n['m=w'],
                                                                                    c['p3', 'm<w'], n['m<w'], c['pu', 'm<w'], n['m<w'], c['bw'] / max(1, len(recs)), c['bwmax'], c['conn']))
    # B. some sigma
    tab('')
    tab('  B. weight spaces where SOME sigma makes V3 (resp. the whole triple) a path / lowers the bandwidth:')
    for k in kinds + ('all',):
        sel = [r for r in recs if k == 'all' or kind_of(r) == k]
        if not sel:
            continue
        some3 = sum(any(r['per'][s]['g3']['paths'] for s in SIGMAS) for r in sel)
        someu = sum(any(r['per'][s]['gu']['paths'] for s in SIGMAS) for r in sel)
        idp = sum(r['per'][ID]['g3']['paths'] for r in sel)
        better = sum(min(r['per'][s]['g3']['bw'] for s in SIGMAS) < r['per'][ID]['g3']['bw'] for r in sel)
        nopath = [r for r in sel if not any(r['per'][s]['g3']['paths'] for s in SIGMAS)]
        tab('     %-5s: %3d spaces; V3 path at sigma = id: %3d; at SOME sigma: %3d (triple union: %3d); bandwidth strictly smaller at some sigma than at id: %3d; NO sigma gives a path: %3d; min over sigma of bw: mean %.2f, max %d (id: mean %.2f, max %d)' %
            (k, len(sel), idp, some3, someu, better, len(nopath), sum(min(r['per'][s]['g3']['bw'] for s in SIGMAS) for r in sel) / len(sel), max(min(r['per'][s]['g3']['bw'] for s in SIGMAS) for r in sel),
             sum(r['per'][ID]['g3']['bw'] for r in sel) / len(sel), max(r['per'][ID]['g3']['bw'] for r in sel)))
    # C. by class A
    tab('')
    tab('  C. by active-string class A (strings E_{i,i+1}^{k_i} that are non-void on the weight space) and kind: number of spaces; at id / at some sigma path; best sigma (max #paths, ties -> smaller mean bw)')
    classes = sorted({(r['A'], kind_of(r)) for r in recs}, key=lambda t: (len(t[0]), sorted(t[0]), t[1]))
    best = {}
    for A, k in classes:
        sel = [r for r in recs if r['A'] == A and kind_of(r) == k]
        sc = []
        for s in SIGMAS:
            sc.append((-sum(r['per'][s]['g3']['paths'] for r in sel), sum(r['per'][s]['g3']['bw'] for r in sel) / len(sel), s))
        sc.sort()
        top = [(sstr(s), -a) for a, b, s in sc[:3]]
        some = sum(any(r['per'][s]['g3']['paths'] for s in SIGMAS) for r in sel)
        idp = sum(r['per'][ID]['g3']['paths'] for r in sel)
        best[(A, k)] = sc[0][2]
        tab('     A = %-9s %-4s: %3d spaces; id: %3d paths; some sigma: %3d; best sigmas (sigma, #paths): %s; all-sigma-path sigmas: %s' % (A_str(A), k, len(sel), idp, some, top,
            [sstr(s) for s in SIGMAS if all(r['per'][s]['g3']['paths'] for r in sel)]))
    return best


# ======================================================================================== 5. coordinate cuts
def cut_inequality(sm, i):
    """for a string whose generator is E_21 in the standard frame (sigma-chain with {sigma1,sigma2} = {i,i+1} and the order making it a lowering E_21): the kernel is {m22 > nu'_1 - k}"""
    s = sm['strings'][i - 1]
    return s


def stage_cuts(recs, label=''):
    hdr('C2  COORDINATE CUTS  %s' % label)
    # which (sigma, string) give a coordinate kernel, split by void / non-void; against the prediction {sigma_1, sigma_2} = {i, i+1}
    tab('  A. for every sigma and string i: is ker E_{tau(i)tau(i+1)}^{k_i} on W^sigma a coordinate subspace?  counts over weight spaces where the string is NON-VOID; "pred" = prediction ({sigma1,sigma2} = {i,i+1})')
    tab('     %-6s | %s' % ('sigma', ' | '.join('i=%d: coord/nonvoid (pred) triv' % i for i in (1, 2, 3))))
    viol = 0
    acc = Counter()
    for s in SIGMAS:
        cells = []
        for i in (1, 2, 3):
            nv = co = triv = 0
            for r in recs:
                p = r['per'][s]
                if not p['void'][i - 1]:
                    nv += 1
                    co += p['coord'][i - 1]
                    # trivial: kernel is everything (nullity = w) or zero (nullity = 0): coordinate for a degenerate reason
                    triv += (p['coord'][i - 1] and (p['nullity'][i - 1] in (0, r['w'])))
            pred = {s[0], s[1]} == {i, i + 1}
            cells.append('%4d/%-4d (%s) %3d' % (co, nv, 'yes' if pred else 'no ', triv))
            if pred:
                viol += (co != nv)
            else:
                acc[(i, s)] += co - triv
        tab('     %-6s | %s' % (sstr(s), ' | '.join(cells)))
    tab('  predicted cuts that fail to be coordinate (non-void): %d;  accidental non-trivial coordinate kernels outside the prediction (sum over spaces): %d' % (viol, sum(acc.values())))
    # B. two active strings
    tab('')
    tab('  B. sectors with >= 2 active strings: is there a sigma in which ALL their (non-void) kernels are coordinate subspaces?  (non-trivially: kernel not 0 and not everything)')
    c = Counter()
    ex = []
    for r in recs:
        if len(r['A']) >= 2:
            good = [s for s in SIGMAS if all(r['per'][s]['coord'][i - 1] and r['per'][s]['nullity'][i - 1] not in (0, r['w']) for i in r['A'])]
            good_triv = [s for s in SIGMAS if all(r['per'][s]['coord'][i - 1] for i in r['A'])]
            c[(A_str(r['A']), 'nontriv', bool(good))] += 1
            c[(A_str(r['A']), 'any', bool(good_triv))] += 1
            if good:
                ex.append((r['lam'], r['lamp'], r['mu'], r['w'], r['m'], [sstr(s) for s in good]))
    for k, v in sorted(c.items()):
        tab('     A = %-9s %-8s: some sigma makes all strings coordinate: %-5s in %d spaces' % (k[0], k[1], k[2], v))
    tab('     counter-examples (non-trivial, two active strings both coordinate in one sigma): %d' % len(ex))
    for e in ex[:10]:
        tab('        %s' % (e,))
    # C. one active string: M^sigma shape
    tab('')
    tab('  C. sectors with ONE active string i: in the sigma with {s1,s2} = {i,i+1}, M^sigma is a coordinate subspace; is the invariant, and is the triple a path on it?')
    for i in (1, 2, 3):
        sel = [r for r in recs if r['A'] == frozenset([i])]
        if not sel:
            continue
        cs = Counter()
        for r in sel:
            for s in SIGMAS:
                if {s[0], s[1]} == {i, i + 1}:
                    p = r['per'][s]
                    cs[(sstr(s), 'M coordinate')] += p['Mcoord']
                    cs[(sstr(s), 'dim M = m')] += (p.get('Mdim') == r['m'])
                    cs[(sstr(s), 'invariant')] += bool(p.get('Minv'))
                    cs[(sstr(s), 'V3|M path')] += bool(p.get('gM3') and p['gM3']['paths'])
                    cs[(sstr(s), 'triple|M path')] += bool(p.get('gMu') and p['gMu']['paths'])
                    cs[(sstr(s), 'V3|M connected')] += bool(p.get('gM3') and p['gM3']['ncomp'] == 1)
        tab('     A = {%d}: %d spaces' % (i, len(sel)))
        for s in sorted({k[0] for k in cs}):
            tab('        sigma = %s: M coordinate %d, dim M^sigma = m %d, triple-invariant %d, V3|M path %d, triple|M path %d, V3|M connected %d' %
                (s, cs[(s, 'M coordinate')], cs[(s, 'dim M = m')], cs[(s, 'invariant')], cs[(s, 'V3|M path')], cs[(s, 'triple|M path')], cs[(s, 'V3|M connected')]))
    return ex


def w0_cut_check(lam, lamp, verbose=True):
    """C0: for sigma = w0 the string E_34^{k_3} becomes E_21^{k_3}; its kernel on W^{w0} is exactly {patterns with m22 > nu'_1 - k_3} (nu'_1 = nu_4); verified on every sector of the product"""
    R = gr.GLNRep.of(lam)
    bad = checked = 0
    for mu, m in sg.sector_list(lam, lamp, 'rs'):
        sm = sigma_mats(lam, lamp, mu, W0)
        if sm is None:
            continue
        s = sm['strings'][2]
        assert s['gen'] == (2, 1), s['gen']
        if s['void']:
            continue
        k = s['k']
        pred = {j for j in range(sm['w']) if sm['pats'][j][1][1] > sm['nup'][0] - k}
        checked += 1
        if pred != s['zero'] or not s['coord']:
            bad += 1
    return checked, bad


def swap12(s):
    return (s[1], s[0], s[2], s[3])


def cut_sigmas(i):
    """the 4 sigma whose innermost gl(2) is gl({i,i+1}) (the sigma in which the string E_{i,i+1} is a coordinate cut)"""
    return [s for s in SIGMAS if {s[0], s[1]} == {i, i + 1}]


def stage_extra(recs, label=''):
    hdr('C1b  STENCIL, sigma-EQUIVALENCE, RULES  %s' % label)
    bad = sum((r['per'][s]['g3'], r['per'][s]['gu']) != (r['per'][swap12(s)]['g3'], r['per'][swap12(s)]['gu']) for r in recs for s in SIGMAS)
    tab('  (a) sigma and sigma o (12) (swap of the first two chain entries) have identical graph statistics of V3 and of the triple in ALL %d x 24 cases: %d differences (the gl(2)_{sigma1 sigma2} Weyl reflection is a monomial map of the GT basis);' % (len(recs), bad))
    tab('      so only the 12 ordered pairs (sigma4, sigma3) give different bases.')
    tab('')
    tab('  (b) stencil of V3 in the sigma-patterns: set of displacements (dm12, dm13, dm23) with a non-zero entry, union over all weight spaces:')
    st = {}
    for s in SIGMAS:
        U = set()
        for r in recs:
            U |= set(r['per'][s]['stencil'])
        st[s] = U
    base = st[ID]
    offd = {(a, b, c) for a in (-1, 0, 1) for (b, c) in [(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]} - {(0, 0, 0)}
    prod = {(a, b, c) for a in (-1, 0, 1) for (b, c) in [(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]}
    tab('      sigma = id: %d off-diagonal displacements; equal to the 20 off-diagonal points of the 21-point set {0,+-1} x ({0} u gl3 roots}: %s' % (len(base), base == offd))
    tab('      all 24 sigma: sizes %s; every stencil inside the 21-point set: %s; every stencil equal to the 20 off-diagonal points: %s' % (sorted(Counter(len(U) for U in st.values()).items()), all(U <= prod for U in st.values()), all(U == offd for U in st.values())))
    tab('')
    nd = [r for r in recs if len({r['per'][s]['g3']['paths'] for s in SIGMAS}) > 1]
    tab('  (c) the path property of V3 is the same in all 24 sigma-bases in %d of %d weight spaces; exceptions:' % (len(recs) - len(nd), len(recs)))
    for r in nd:
        tab('      %s mu=%s w=%d m=%d: path in sigma = %s; at id: edges %d, maxdeg %d' % (pstr(r['lam'], r['lamp']), r['mu'], r['w'], r['m'], [sstr(s) for s in SIGMAS if r['per'][s]['g3']['paths']],
            r['per'][ID]['g3']['edges'], r['per'][ID]['g3']['maxdeg']))
    c = Counter()
    for r in recs:
        e = {r['per'][s]['g3']['edges'] for s in SIGMAS}
        c[(kind_of(r), len(e) == 1)] += 1
    tab('      the number of EDGES of the V3 graph is sigma-independent in: %s' % dict(c))
    tab('')
    tab("  (d) candidate rules 'active strings A -> sigma': number of spaces where V3 is a path / mean RCM bandwidth of V3 at the rule's sigma (id and the optimum over sigma for comparison)")
    for A in sorted({r['A'] for r in recs}, key=lambda a: (len(a), sorted(a))):
        sel = [r for r in recs if r['A'] == A]

        def ev(fn):
            ss = [fn(r) for r in sel]
            return (sum(r['per'][x]['g3']['paths'] for r, x in zip(sel, ss)), sum(r['per'][x]['g3']['bw'] for r, x in zip(sel, ss)) / len(sel))
        pid = ev(lambda r: ID)
        popt = (sum(any(r['per'][x]['g3']['paths'] for x in SIGMAS) for r in sel), sum(min(r['per'][x]['g3']['bw'] for x in SIGMAS) for r in sel) / len(sel))
        rules = []
        for i in sorted(A):
            rules.append(('cut sigma of string %d (%s)' % (i, sstr(cut_sigmas(i)[0])), ev(lambda r, i=i: cut_sigmas(i)[0])))
        if A:
            rules.append(('w0', ev(lambda r: W0)))
        tab('     A = %-9s (%3d spaces): id: %3d paths, bw %.2f | optimum over sigma: %3d paths, bw %.2f | %s' % (A_str(A), len(sel), pid[0], pid[1], popt[0], popt[1],
            ' | '.join('%s: %d paths, bw %.2f' % (nm, v[0], v[1]) for nm, v in rules)))


def stage_cuts2(recs, label=''):
    hdr('C2b  CUT INEQUALITY AND THE CUT SUB-LATTICE  %s' % label)
    chk = bad = 0
    for r in recs:
        for i in (1, 2, 3):
            if i not in r['A']:
                continue
            thr = r['mu'][i - 1] - r['lamp'][i]
            for s in cut_sigmas(i):
                p = r['per'][s]
                pred = tuple(j for j, c in enumerate(p['coords']) if c[0] <= thr)
                chk += 1
                bad += (pred != p['zero'][i - 1] or not p['coord'][i - 1])
    tab("  (a) cut inequality in the sigma-coordinates:  for {sigma1, sigma2} = {i, i+1}:   ker E_{i,i+1}^{k_i} on W^sigma = { patterns with m12 <= mu_i - lam'_{i+1} }   (m12 = highest weight of the gl({i,i+1}) irrep)")
    tab("      (equivalently m22 >= mu_{i+1} - lam'_i); same threshold for E_12-type (sigma1 = i) and E_21-type (sigma1 = i+1) chains.  Checked on %d (space, string, sigma) triples: %d failures" % (chk, bad))
    tab('')
    tab('  (b) coordinate sub-lattice cut out by the strings that are coordinate in sigma (all non-void such strings), compared with M_mu: invariance under the triple; V3 path on it; dimension = m?')
    tab('      %-6s | %-17s | %-9s %-11s %-6s %-10s | %s' % ('sigma', 'spaces with a cut', 'invariant', 'V3|cut path', 'dim=m', 'mean dim/w', 'whole-W V3 path (same spaces)'))
    for s in SIGMAS:
        if s[0] > s[1]:
            continue
        sel = [r for r in recs if 'cdim' in r['per'][s]]
        if not sel:
            continue
        tab('      %-6s | %-17d | %-9d %-11d %-6d %-10.2f | %d' % (sstr(s), len(sel), sum(r['per'][s]['cinv'] for r in sel), sum(r['per'][s]['cg3']['paths'] for r in sel), sum(r['per'][s]['cdim'] == r['m'] for r in sel),
            sum(r['per'][s]['cdim'] / r['w'] for r in sel) / len(sel), sum(r['per'][s]['g3']['paths'] for r in sel)))
    tab('')
    tab('  (c) sectors with exactly ONE active string i (sigma = first cut sigma of i): M^sigma = {m12 <= thr}; triple restricted to it a path?')
    for i in (1, 2, 3):
        sel = [r for r in recs if r['A'] == frozenset([i])]
        if not sel:
            continue
        s = cut_sigmas(i)[0]
        tab('      A = {%d} (sigma %s): %d spaces; M coordinate: %d; dim M = m: %d; invariant: %d; V3|M path: %d; triple|M path: %d; V3|M connected: %d; whole W V3 path at id: %d; both paths: %d; W path & M not: %d; M path & W not: %d' % (
            i, sstr(s), len(sel), sum(r['per'][s]['Mcoord'] for r in sel), sum(r['per'][s].get('Mdim') == r['m'] for r in sel), sum(bool(r['per'][s].get('Minv')) for r in sel),
            sum(bool(r['per'][s].get('gM3') and r['per'][s]['gM3']['paths']) for r in sel), sum(bool(r['per'][s].get('gMu') and r['per'][s]['gMu']['paths']) for r in sel),
            sum(bool(r['per'][s].get('gM3') and r['per'][s]['gM3']['ncomp'] == 1) for r in sel), sum(r['per'][ID]['g3']['paths'] for r in sel),
            sum(bool(r['per'][s]['gM3']['paths'] and r['per'][ID]['g3']['paths']) for r in sel), sum(bool(not r['per'][s]['gM3']['paths'] and r['per'][ID]['g3']['paths']) for r in sel),
            sum(bool(r['per'][s]['gM3']['paths'] and not r['per'][ID]['g3']['paths']) for r in sel)))


import sympy as sp


def spmat(M):
    return sp.Matrix([[sp.Rational(v.numerator, v.denominator) for v in r] for r in M])


def mat_str(M, width=9):
    return '\n'.join('        [' + ', '.join('%*s' % (width, ('%s' % x) if x.denominator == 1 else '%d/%d' % (x.numerator, x.denominator)) for x in r) + ' ]' for r in M)


def fac_str(p):
    return str(sp.factor(p))


def galois_blocks(mats, combo=so.GENERIC):
    """rational invariant decomposition: the irreducible factors over Q of the charpoly of the generic combination G of the (restricted) triple; for each factor f the
    rational subspace ker f(G) (basis, exact) and the restricted exact matrices.  Returns list of dict(f, deg, basis (sympy Matrix w x d), R {name: d x d sympy Matrix})"""
    x = sp.symbols('x')
    G = spmat(lincomb(mats, combo))
    cp = G.charpoly(x).as_expr()
    out = []
    for f, mult in sp.factor_list(cp)[1]:
        fm = sp.zeros(G.shape[0], G.shape[0])
        P = sp.Poly(f, x)
        for k, c in enumerate(reversed(P.all_coeffs())):
            fm += c * G ** k
        ns = fm.nullspace()
        B = sp.Matrix.hstack(*ns)
        out.append(dict(f=f, mult=mult, deg=P.degree(), basis=B))
    return out


def block_matrices(B, mats):
    """exact d x d matrices R_k with  V_k B = B R_k  (B = rational basis of an invariant subspace)"""
    out = []
    Bt = B.T
    for M in mats:
        V = spmat(M)
        R = (Bt * B).inv() * Bt * V * B
        assert V * B == B * R
        out.append(R)
    return out


def eigvecs_2(R3, R4a, R4b):
    """exact joint eigenvectors of commuting 2 x 2 matrices over their common quadratic field (radical form)"""
    x = sp.symbols('x')
    res = []
    cpoly = R3.charpoly(x).as_expr()
    # use the matrix with distinct eigenvalues among the three
    for R in (R3, R4a, R4b):
        if len(set(R.eigenvals())) == 2:
            ev = list(R.eigenvals())
            break
    else:
        return None
    for e in ev:
        v = (R - e * sp.eye(2)).nullspace()[0]
        v = sp.simplify(v / v[0]) if v[0] != 0 else v
        trip = []
        for Rk in (R3, R4a, R4b):
            w = Rk * v
            idx = 0 if v[0] != 0 else 1
            trip.append(sp.simplify(sp.radsimp(sp.simplify(w[idx] / v[idx]))))
        res.append((v, trip))
    return res


def nested_basis_graph(sm, tol=mp.mpf(10) ** -30, dps=50):
    """change to the within-path eigenbasis: group the lattice points by the sigma-row-3 (kappa = (m13, m23)), diagonalise the kappa-diagonal block of V3 (50 digits), conjugate the triple;
    returns (graph stats of V3 and of the triple union in the new basis with entries < tol treated as zero, number of nonzero off-diagonal entries)"""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        w = sm['w']
        groups = defaultdict(list)
        for a, c in enumerate(sm['coords']):
            groups[(c[1], c[2])].append(a)
        U = mp.zeros(w, w)
        for kap, idx in groups.items():
            if len(idx) == 1:
                U[idx[0], idx[0]] = 1
                continue
            blk = mp.matrix([[mp.mpf(sm['mats']['V3'][a][b].numerator) / sm['mats']['V3'][a][b].denominator for b in idx] for a in idx])
            E, Uk = mp.eig(blk)
            for jj, b in enumerate(idx):
                for ii, a in enumerate(idx):
                    U[a, b] = Uk[ii, jj]
        Ui = mp.inverse(U)
        newm = {}
        for n in TRIPLE:
            M = mp.matrix([[mp.mpf(v.numerator) / v.denominator for v in r] for r in sm['mats'][n]])
            newm[n] = Ui * M * U
        def to_graph(names):
            Mx = [[Fr(0)] * w for _ in range(w)]
            nz = 0
            for a in range(w):
                for b in range(w):
                    if a != b and any(abs(newm[n][a, b]) > tol for n in names):
                        Mx[a][b] = Fr(1)
                        nz += 1
            return graph(Mx), nz
        g3, nz3 = to_graph(('V3',))
        gu, nzu = to_graph(TRIPLE)
        return g3, nz3, gu, nzu, len(groups), max(len(i) for i in groups.values())
    finally:
        mp.mp.dps = old


def anatomy(lam, lamp, mu, m, full_mats=True):
    """the whole anatomy of one sector (C3): patterns, exact matrices, spectra over Q, graphs in all sigma-bases, nested-basis graph"""
    tab('')
    tab('-' * 140)
    tab('SECTOR %s  mu = %s   (w = dim W, m = dim M_mu)' % (pstr(lam, lamp), tuple(mu)))
    sm = sigma_mats(lam, lamp, mu, ID)
    w = sm['w']
    prv = se.prv_space(lam, lamp, mu)
    tab('  w = %d, m = %d, active strings A = %s' % (w, prv.m, A_str(sm['A'])))
    tab('  lattice points (id chain), coordinates (m12, m13, m23) and full patterns:')
    for a in range(w):
        tab('     %d: coords %s   pattern %s' % (a + 1, sm['coords'][a], sm['pats'][a][:3]))
    if full_mats and w <= 6:
        for n in TRIPLE:
            tab('  %s on W (column j = image of lattice point j):' % n)
            tab(mat_str(sm['mats'][n]))
    # physical copies (matrices on M_mu, PRV coordinates)
    R = sm['R']
    if prv.m == w:
        mm = [sm['mats'][n] for n in TRIPLE]
    else:
        mm = [so.vred_matrix(R, so.reduced_op(n, lamp), prv) for n in TRIPLE]
        tab('  M_mu (PRV basis, free positions %s of W): V3|M = ' % (list(prv.free),))
        tab(mat_str(mm[0]))
    x = sp.symbols('x')
    for n, M in zip(TRIPLE, mm):
        tab('  charpoly of %s on M_mu = %s' % (n, fac_str(spmat(M).charpoly(x).as_expr())))
    blocks = galois_blocks(mm)
    tab('  generic combination G = V3 + 7/19 V4a + 11/13 V4b: irreducible factors over Q of its charpoly: degrees %s' % [b['deg'] for b in blocks])
    for b in blocks:
        tab('     factor %s' % sp.factor(b['f']))
    # joint eigenvectors exactly in the quadratic blocks
    for b in blocks:
        if b['deg'] == 2:
            Rs = block_matrices(b['basis'], mm)
            ev = eigvecs_2(*Rs)
            tab('     rational 2-dim invariant block, basis (in M_mu coordinates) columns: %s ; the triple on it:' % (b['basis'].T.tolist(),))
            for n, Rk in zip(TRIPLE, Rs):
                tab('        %s|block = %s' % (n, Rk.tolist()))
            if ev:
                for v, trip in ev:
                    tab('        joint eigenvector (block coords) %s   joint eigenvalues (V3, V4a, V4b) = %s' % (list(v), trip))
        elif b['deg'] == 1:
            Rs = block_matrices(b['basis'], mm)
            tab('     rational line: vector %s, joint eigenvalues %s' % (b['basis'].T.tolist(), [Rk[0, 0] for Rk in Rs]))
        else:
            tab('     degree-%d block: joint eigenvalues are roots of %s (numbers not written)' % (b['deg'], sp.factor(b['f'])))
    # graphs in all sigma bases
    tab('  sparsity graph of V3 (and of the triple) in the 24 sigma-bases [grouped by sigma ~ sigma o (12)]:')
    rows = []
    allsm = {}
    for s in SIGMAS:
        if s[0] > s[1]:
            continue
        sms = sigma_mats(lam, lamp, mu, s)
        allsm[s] = sms
        g3 = graph(sms['mats']['V3'])
        gu = graph(union_matrix(sms['mats']))
        ng = nested_basis_graph(sms)
        rows.append((s, g3, gu, ng))
        tab('     sigma %s: V3 edges %d maxdeg %d path %-5s bw %d | triple edges %d path %-5s bw %d | cut strings %s | nested (kappa-eigenbasis of V3): %d kappa-groups (max size %d), V3 edges %d path %s, triple edges %d path %s' %
            (sstr(s), g3['edges'], g3['maxdeg'], g3['paths'], g3['bw'], gu['edges'], gu['paths'], gu['bw'], sms.get('A') and [x_['i'] for x_ in sms['strings'] if x_['coord'] and not x_['void']],
             ng[4], ng[5], ng[0]['edges'], ng[0]['paths'], ng[2]['edges'], ng[2]['paths']))
    best = min(rows, key=lambda t: (not t[1]['paths'], t[1]['bw'], t[2]['bw']))
    tab('  best sigma (path first, then smaller bandwidth of V3): %s (V3 path %s, bw %d); number of the 12 classes with V3 path: %d' % (sstr(best[0]), best[1]['paths'], best[1]['bw'], sum(t[1]['paths'] for t in rows)))
    return dict(sm=sm, rows=rows, blocks=blocks)


def pick_sectors(recs, products=P_SMALL):
    """K4 first, then the next 3 smallest m = w spaces whose V3 graph at id is not a path, and the smallest m < w sector with a cycle in the V3 graph"""
    nonpath = [r for r in recs if r['m'] == r['w'] and not r['per'][ID]['g3']['paths']]
    nonpath.sort(key=lambda r: (r['w'], r['per'][ID]['g3']['edges'], r['lam'], r['mu']))
    cyc = [r for r in recs if r['m'] < r['w'] and r['per'][ID]['g3']['edges'] > r['w'] - r['per'][ID]['g3']['ncomp']]
    cyc.sort(key=lambda r: (r['w'], r['m'], r['lam'], r['mu']))
    return nonpath, cyc


def stage_small():
    hdr('C3  THE SMALLEST NON-CHAIN SECTORS BY HAND')
    recs = compute_records(P_SMALL, verbose=False)
    nonpath, cyc = pick_sectors(recs)
    tab('  m = w spaces (P_small) whose V3 graph at sigma = id is not a path: %d; the smallest ones (w, edges, product, mu):' % len(nonpath))
    for r in nonpath[:10]:
        tab('     w=%d edges=%d %s mu=%s (A = %s)' % (r['w'], r['per'][ID]['g3']['edges'], pstr(r['lam'], r['lamp']), r['mu'], A_str(r['A'])))
    tab('  m < w spaces with a cycle: %d; the smallest (w, m, edges, product, mu):' % len(cyc))
    for r in cyc[:6]:
        tab('     w=%d m=%d edges=%d %s mu=%s (A = %s)' % (r['w'], r['m'], r['per'][ID]['g3']['edges'], pstr(r['lam'], r['lamp']), r['mu'], A_str(r['A'])))
    k4 = [r for r in nonpath if r['lam'] == (3, 2, 1, 0) and r['lamp'] == (3, 2, 1, 0) and r['mu'] == (5, 4, 2, 1)]
    assert k4
    chosen = [k4[0]]
    seen = {(3, 2, 1, 0, 3, 2, 1, 0, 5, 4, 2, 1)}
    seen = {((3, 2, 1, 0), (3, 2, 1, 0), (5, 4, 2, 1))}
    for r in nonpath:
        key = (r['lam'], r['lamp'], r['mu'])
        if key not in seen and len(chosen) < 4:
            chosen.append(r)
            seen.add(key)
    chosen.append(cyc[0])
    out = []
    for r in chosen:
        out.append(anatomy(r['lam'], r['lamp'], r['mu'], r['m']))
    return out


def stage_galois(products=None):
    hdr('C3c  RATIONAL (Galois) BLOCK STRUCTURE OF THE COPIES: degrees of the irreducible factors over Q of the charpoly of the generic combination G on M_mu')
    products = products or P_SMALL + [P_BIG[1], P_BIG[2]]
    x = sp.symbols('x')
    for lam, lamp in products:
        R = gr.GLNRep.of(lam)
        c = Counter()
        ex = []
        for mu, m in sg.sector_list(lam, lamp, 'rs'):
            if m < 2:
                continue
            prv = se.prv_space(lam, lamp, mu)
            mm = [so.vred_matrix(R, so.reduced_op(n, lamp), prv) for n in TRIPLE]
            G = spmat(lincomb(mm, so.GENERIC))
            degs = sorted(sp.Poly(f, x).degree() * 1 for f, mult in sp.factor_list(G.charpoly(x).as_expr())[1] for _ in range(mult))
            kind = 'irreducible' if degs == [m] else ('all rational' if max(degs) == 1 else 'split')
            c[kind] += 1
            if kind == 'split' and len(ex) < 8:
                ex.append((mu, m, degs))
        tab('  %s: sectors with m >= 2: %s; examples of split sectors (mu, m, degrees): %s' % (pstr(lam, lamp), dict(c), ex))


def exchange_targets(mats):
    """for lam = lam': the affine images demanded of the triple by the site exchange P:  P V3 P = -V3 + c,  P V4a P = -2 V4a - V4b + c',  P V4b P = 3 V4a + 2 V4b - c'  (G3 section 8.1; c, c' from the traces)"""
    m = len(mats['V3'])
    tr = {n: sum((mats[n][a][a] for a in range(m)), Fr(0)) for n in TRIPLE}
    c = 2 * tr['V3'] / m
    cp = (3 * tr['V4a'] + tr['V4b']) / m
    I = [[Fr(int(a == b)) for b in range(m)] for a in range(m)]
    def comb(cs, ci):
        return [[sum((k * mats[n][a][b] for k, n in cs), Fr(0)) + ci * I[a][b] for b in range(m)] for a in range(m)]
    return {'V3': comb([(-1, 'V3')], c), 'V4a': comb([(-2, 'V4a'), (-1, 'V4b')], cp), 'V4b': comb([(3, 'V4a'), (2, 'V4b')], -cp)}, c, cp


def monomial_intertwiners(mats, targets, limit=50):
    """all monomial matrices S (S e_a = d_a e_{pi(a)}) with S V S^-1 = target(V) for the three operators (backtracking on the diagonal and zero pattern, then scales along a spanning tree);
    returns list of (pi, d) with d[0] = 1"""
    m = len(mats['V3'])
    ops = [(mats[n], targets[n]) for n in TRIPLE]
    cand = [[b for b in range(m) if all(V[a][a] == Tg[b][b] for V, Tg in ops)] for a in range(m)]
    res = []
    pi = [None] * m
    used = set()

    def scales(pi):
        d = [None] * m
        d[0] = Fr(1)
        stack = [0]
        while stack:
            a = stack.pop()
            for b in range(m):
                if b == a:
                    continue
                for V, Tg in ops:
                    if V[a][b] != 0:
                        # d_a V_ab / d_b = T_{pi a, pi b}
                        t = Tg[pi[a]][pi[b]]
                        if t == 0:
                            return None
                        db = d[a] * V[a][b] / t
                        if d[b] is None:
                            d[b] = db
                            stack.append(b)
                        elif d[b] != db:
                            return None
        if any(x is None for x in d):
            return None
        for a in range(m):
            for b in range(m):
                for V, Tg in ops:
                    if d[a] * V[a][b] != Tg[pi[a]][pi[b]] * d[b]:
                        return None
        return d

    def rec(a):
        if len(res) >= limit:
            return
        if a == m:
            d = scales(pi)
            if d is not None:
                res.append((tuple(pi), tuple(d)))
            return
        for b in cand[a]:
            if b in used:
                continue
            # zero-pattern prefilter against earlier assignments
            ok = True
            for a2 in range(a):
                for V, Tg in ops:
                    if (V[a2][a] == 0) != (Tg[pi[a2]][b] == 0) or (V[a][a2] == 0) != (Tg[b][pi[a2]] == 0):
                        ok = False
                        break
                if not ok:
                    break
            if not ok:
                continue
            pi[a] = b
            used.add(b)
            rec(a + 1)
            used.discard(b)
            pi[a] = None
    rec(0)
    return res


def is_involutive(pi_d):
    pi, d = pi_d
    m = len(pi)
    if any(pi[pi[a]] != a for a in range(m)):
        return False
    k = d[0] * d[pi[0]]
    return all(d[a] * d[pi[a]] == k for a in range(m))


def stage_parity(products=(((3, 2, 1, 0), (3, 2, 1, 0)), ((4, 2, 1, 0), (4, 2, 1, 0)))):
    hdr('C4  THE EXCHANGE PARITY AS A LATTICE SYMMETRY (lam = lam\')')
    tab('  for each sector with m >= 2 and each sigma in which M^sigma is a coordinate sub-lattice: does a MONOMIAL matrix S (S e_a = d_a e_pi(a), exact) exist with S V3 S^-1 = -V3 + c, S V4a S^-1 = -2 V4a - V4b + c\', S V4b S^-1 = 3 V4a + 2 V4b - c\' on M^sigma?')
    tab('  (S then maps the joint eigenlines like P does; it equals P only up to a function of the triple, so this is a necessary condition for P to be a lattice reflection.)')
    for lam, lamp in products:
        tot = Counter()
        rows = []
        for mu, m in sg.sector_list(lam, lamp, 'rs'):
            if m < 2:
                continue
            found = {}
            ncoord = 0
            for s in SIGMAS:
                if s[0] > s[1]:
                    continue
                sm = sigma_mats(lam, lamp, mu, s)
                if not all(x['coord'] for x in sm['strings']):
                    continue
                ncoord += 1
                idx = sorted(set.intersection(*[x['zero'] for x in sm['strings']]))
                sub = {n: submatrix(sm['mats'][n], idx) for n in TRIPLE}
                tg, c, cp = exchange_targets(sub)
                sols = monomial_intertwiners(sub, tg)
                found[s] = (len(sols), sum(is_involutive(x) for x in sols), [x[0] for x in sols[:1]])
            ids = [s for s, v in found.items() if v[0] > 0]
            tot['sectors'] += 1
            tot['with coordinate M^sigma in some sigma'] += bool(found)
            tot['monomial S exists in some sigma'] += bool(ids)
            tot['involutive monomial S in some sigma'] += any(v[1] > 0 for v in found.values())
            rows.append((mu, m, ncoord, [(sstr(s), v[0], v[1]) for s, v in found.items() if v[0] > 0]))
        tab('  %s: %s' % (pstr(lam, lamp), dict(tot)))
        for mu, m, nc, fnd in rows:
            tab('      mu=%s m=%d: sigma-classes with coordinate M: %d; monomial S found in (sigma, #solutions, #involutive): %s' % (mu, m, nc, fnd[:6]))


# ======================================================================================== main
def main(argv):
    stages = argv[1:] or ['all']
    t0 = time.perf_counter()
    if 'run' in stages:
        stages = [x for x in ','.join(argv[argv.index('run') + 1:]).split(',') if x]
    if 'all' in stages:
        stages = ['check', 'bandwidth', 'cuts', 'small', 'parity']
    state = {}
    if 'check' in stages:
        stage_check()
    if 'small' in stages:
        stage_small()
        stage_galois()
    if 'bandwidth' in stages or 'cuts' in stages:
        small = compute_records(P_SMALL)
        big = compute_records(P_BIG)
        allr = small + big
        stage_bandwidth(small, '[P_small]')
        stage_bandwidth(allr, '[P_small + P_big]')
        stage_extra(allr, '[P_small + P_big]')
        stage_cuts(allr, '[P_small + P_big]')
        stage_cuts2(allr, '[P_small + P_big]')
    if 'parity' in stages:
        stage_parity()
    tab('')
    tab('(total %.1f s)' % (time.perf_counter() - t0))
    with open(TABLES_PATH, 'w', encoding='utf-8') as f:
        f.write('\n'.join(_TAB) + '\n')


if __name__ == '__main__':
    main(sys.argv)
