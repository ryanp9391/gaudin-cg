"""G7: literal comparison of Munich's HW Clebsch-Gordan coefficients (ClebschGordan.cpp via Paul/Cpp/AlexCGC/hwcgc.sh) with the slice engine's.

  python -X utf8 munich_bridge.py selftest       synthetic Munich outputs from our own vectors (random rotation of the copies for m >= 2): accepted; corrupted: rejected
  python -X utf8 munich_bridge.py compare [2a|2b|all]   the whole D1 comparison (raw Munich outputs go to Paul/Cpp/AlexCGC/results/su4-same-machine/)
  python -X utf8 munich_bridge.py sector LAM LAMP MU   one product/irrep (labels as 4 comma-free ints each, OUR gl(4) labels), e.g.  sector 3,2,1,0 3,2,1,0 5,4,2,1

Conventions (RESULTS-G1 section 9): Munich Q(M) = our GT index + 1; Munich's orthonormal basis = ours rescaled by 1/sqrt(g), trivial sign gauge;
Munich's SU(4) label = ours minus mu_4 (1,1,1,1).  Site 1 = first factor (index Q(M)), site 2 = second factor (index Q(M')).
"""
import os
import re
import subprocess
import sys
import time

import numpy as np

import glN_reps as gr
import slice_engine_gl4 as se
import slice_ops_gl4 as so

HERE = os.path.dirname(os.path.abspath(__file__))
CGC_DIR = os.path.normpath(os.path.join(HERE, '..', '..', 'Cpp', 'AlexCGC'))
OUT_DIR = os.path.join(CGC_DIR, 'results', 'su4-same-machine')
TOL_O = 1e-9      # Munich prints 15 decimals
TOL_SV = 1e-8


def _s(x):
    return ' '.join(str(int(v)) for v in x)


def munich_label(mu):
    return tuple(int(x) - int(mu[3]) for x in mu)


def run_hwcgc(lam, lamp, mu_m, alpha, fname, timeout=1800, reuse=True):
    """run hwcgc.sh (N=4) and save the output; returns (text, seconds, status)"""
    path = os.path.join(OUT_DIR, fname)
    if reuse and os.path.exists(path) and os.path.getsize(path) > 0:
        return open(path).read(), 0.0, 'cached'
    if reuse and os.path.exists(path + '.TIMEOUT'):
        return '', 0.0, 'TIMEOUT'
    env = dict(os.environ, LC_ALL='C')
    cmd = ['bash', 'hwcgc.sh', '4', _s(lam), _s(lamp), _s(mu_m), str(alpha)]
    t0 = time.perf_counter()
    try:
        r = subprocess.run(cmd, cwd=CGC_DIR, capture_output=True, text=True, timeout=timeout, env=env)
        text, st = r.stdout, 'ok' if r.returncode == 0 else 'rc=%d' % r.returncode
    except subprocess.TimeoutExpired:
        text, st = '', 'TIMEOUT'
        open(path + '.TIMEOUT', 'w').write('hwcgc.sh exceeded %d s\n' % timeout)
    dt = time.perf_counter() - t0
    if st == 'ok':
        with open(path, 'w') as f:
            f.write(text)
    return text, dt, st


def bench_texts(lam, lamp, mu_m, mode='eig', fname=None, timeout=1800):
    """all m alpha copies at once from the repaired profiler (cgc_hw_prof_ob, HW_MODE=mode; prints the HW coefficients per alpha when not asked for 'machine' output,
    12 decimals).  Returns (list of hwcgc.sh-style texts, seconds, status)"""
    path = os.path.join(OUT_DIR, fname) if fname else None
    if path and os.path.exists(path) and os.path.getsize(path) > 0:
        out, dt, st = open(path).read(), 0.0, 'cached'
    else:
        env = dict(os.environ, LC_ALL='C', HW_MODE=mode, OPENBLAS_NUM_THREADS='1')
        cmd = ['./cgc_hw_prof_ob', '4'] + [str(x) for x in list(lam) + list(lamp) + list(mu_m)] + ['1']
        t0 = time.perf_counter()
        try:
            r = subprocess.run(cmd, cwd=CGC_DIR, capture_output=True, text=True, timeout=timeout, env=env)
            out, st = r.stdout, 'ok' if r.returncode == 0 else 'rc=%d' % r.returncode
        except subprocess.TimeoutExpired:
            return [], timeout, 'TIMEOUT'
        dt = time.perf_counter() - t0
        if st == 'ok' and path:
            open(path, 'w').write(out)
    dimmu = gr.weyl_dim(gr.munich_label_to_ours(lam, lamp, mu_m))
    texts, cur = [], None
    for line in out.splitlines():
        if re.match(r'\s*alpha = \d+:', line):
            cur = []
            texts.append(cur)
        elif cur is not None:
            f = line.split()
            if len(f) == 3 and all(re.fullmatch(r'-?\d+', x) for x in f[:2]):
                cur.append('%s\t%s\t%d\t%s' % (f[0], f[1], dimmu, f[2]))
    return ['\n'.join(t) + '\n' for t in texts], dt, st


def hw_rows(text, dimmu):
    rows = se.parse_munich_rows(text)
    return {(a, b): v for (a, b, c), v in rows.items() if c == dimmu}


# ---------------------------------------------------------------------------------------------------------------- 2a
def compare_mf(lam, lamp, mu, fname=None, tol=TOL_O):
    """multiplicity-free: one vector, up to a global sign.  mu = OUR label.  returns dict"""
    mm = munich_label(mu)
    fname = fname or 'hw_%s_%s_%s.txt' % ('x'.join(map(str, lam)), 'x'.join(map(str, lamp)), '_'.join(map(str, mm)))
    text, dt, st = run_hwcgc(lam, lamp, mm, 1, fname)
    if st not in ('ok', 'cached'):
        if st == 'TIMEOUT':     # hwcgc.sh builds the full descendant table of mu (too big): fall back to the repaired profiler, which prints only the HW vector (HW_MODE=eig, same algorithm class)
            r = compare_sector(lam, lamp, mu, source='bench')
            if 'sv' in r:
                return dict(status='MATCH(prof)' if r['status'] == 'OK' else 'FAIL(prof)', n_o=r['n_ours'], n_m=r['n_mun'], mism=r['n_ours'] - r['n_mun'], dev=1 - float(r['sv'][0]) if False else r['maxdev'],
                            sigma=float(np.sign(r['Q'][0, 0])), secs=dt)
        return dict(status=st, secs=dt)
    n_o, n_m, mism, dev, sigma = se.munich_compare(lam, lamp, mu, text, tol=tol)
    return dict(status='MATCH' if (mism == 0 and dev <= tol) else 'FAIL', n_o=n_o, n_m=n_m, mism=mism, dev=dev, sigma=sigma, secs=dt)


# ---------------------------------------------------------------------------------------------------------------- 2b
def sparse_from_rows(rows, index):
    """{(Q,Q'): v} -> dense vector over `index` ({(Q,Q'): position}); keys outside index are returned as leftover"""
    v = np.zeros(len(index))
    left = 0.0
    for k, x in rows.items():
        if k in index:
            v[index[k]] = x
        else:
            left += x * x
    return v, left


def compare_sector(lam, lamp, mu, dir=None, tol_sv=TOL_SV, alpha_text=None, verbose=False, source='hwcgc'):
    """subspace comparison for the irrep mu (OUR label) of V_lam (x) V_lamp.  alpha_text (selftest): list of Munich outputs per alpha instead of running hwcgc.sh.
    returns dict(m, sv, maxdev, orth_munich, orth_ours, Q, signed_perm, secs, status)"""
    lam, lamp, mu = tuple(lam), tuple(lamp), tuple(mu)
    mult = gr.racah_speiser(lam, lamp)[mu]
    mm = munich_label(mu)
    dimmu = gr.weyl_dim(mu)
    tag = '%s__%s__%s' % ('_'.join(map(str, lam)), '_'.join(map(str, lamp)), '_'.join(map(str, mm)))
    texts, secs = [], 0.0
    if source == 'bench':
        texts, secs, st = bench_texts(lam, lamp, mm, 'eig', 'bench_eig_%s.txt' % tag)
        if st not in ('ok', 'cached') or len(texts) != mult:
            return dict(m=mult, status=st if st not in ('ok', 'cached') else 'alpha count %d != m=%d' % (len(texts), mult), secs=secs)
    for a in range(1, mult + 1):
        if source == 'bench':
            break
        if alpha_text is not None:
            texts.append(alpha_text[a - 1])
            continue
        text, dt, st = run_hwcgc(lam, lamp, mm, a, 'hw_%s_a%d.txt' % (tag, a))
        secs += dt
        if st not in ('ok', 'cached'):
            return dict(m=mult, status=st, secs=secs)
        texts.append(text)
    Mrows = [hw_rows(t, dimmu) for t in texts]
    # ours
    prv = se.prv_space(lam, lamp, mu)
    assert prv.m == mult
    cps = so.copies(lam, lamp, mu, dps=30, prv=prv)
    hw = se.hw_from_slice(lam, lamp, mu, so.copy_float_slices(prv, cps))
    O = [se.to_munich(lam, lamp, hw, k) for k in range(mult)]
    keys = {}
    for d in Mrows + O:
        for k in d:
            keys.setdefault(k, len(keys))
    Mmat = np.zeros((len(keys), mult))
    Omat = np.zeros((len(keys), mult))
    for j in range(mult):
        Mmat[:, j] = sparse_from_rows(Mrows[j], keys)[0]
        Omat[:, j] = sparse_from_rows(O[j], keys)[0]
    gm = Mmat.T @ Mmat
    orth_m = float(np.abs(gm - np.eye(mult)).max())
    Onorm = Omat / np.linalg.norm(Omat, axis=0)
    go_ = Onorm.T @ Onorm
    orth_o = float(np.abs(go_ - np.eye(mult)).max()) if mult > 1 else 0.0
    Q = Onorm.T @ Mmat                          # rotation: Q[k, a] = <ours_k | Munich_alpha>
    sv = np.linalg.svd(Q, compute_uv=False)
    maxdev = float(np.abs(1 - sv).max())
    # residual of Munich span outside the span of ours (in the full space, catches support differences)
    resid = float(np.linalg.norm(Mmat - Onorm @ Q))
    signed_perm = bool(np.allclose(np.sort(np.abs(Q).ravel())[-mult:], 1, atol=1e-6) and np.allclose(np.abs(Q).sum(1), 1, atol=1e-6)) if mult > 1 else True
    ok = maxdev <= tol_sv and orth_m <= TOL_O and resid <= 1e-7
    return dict(m=mult, status='OK' if ok else 'FAIL', sv=sv, maxdev=maxdev, orth_munich=orth_m, orth_ours=orth_o, resid=resid, Q=Q, signed_perm=signed_perm,
                n_ours=sum(1 for x in O[0].values() if abs(x) > TOL_O), n_mun=sum(1 for x in Mrows[0].values() if abs(x) > TOL_O), secs=secs, tag=tag)


def fmt_Q(Q):
    return '\n'.join('      [' + ' '.join('%8.3f' % x for x in row) + ']' for row in Q)


# ---------------------------------------------------------------------------------------------------------------- driver lists
def mf_list():
    L = []
    for S in (1, 2, 3, 4):
        lam = (S, S, 0, 0)
        for mu in sorted(gr.racah_speiser(lam, lam), reverse=True):
            L.append((lam, lam, mu))
    for S, mus in ((6, ['0 0 0 0', '2 1 1 0', '6 4 2 0', '12 12 0 0']), (8, ['0 0 0 0', '2 1 1 0', '8 6 2 0', '16 16 0 0'])):
        lam = (S, S, 0, 0)
        for m in mus:
            mm = tuple(int(x) for x in m.split())
            mu = gr.munich_label_to_ours(lam, lam, mm)
            L.append((lam, lam, mu))
    return L


PRODUCTS = [((2, 1, 0, 0), (2, 1, 0, 0)), ((2, 1, 1, 0), (2, 1, 1, 0)), ((2, 1, 0, 0), (1, 1, 0, 0)), ((1, 1, 0, 0), (2, 1, 0, 0)),
            ((2, 2, 1, 0), (2, 1, 1, 0)), ((3, 2, 1, 0), (3, 2, 1, 0)), ((4, 2, 0, 0), (4, 2, 0, 0)), ((4, 2, 1, 0), (4, 2, 1, 0))]


def run_2a(verbose=True):
    fails = 0
    print('== 2a: multiplicity-free irreps of [S,S,0,0]^2 vector by vector (up to one global sign), tolerance %g' % TOL_O, flush=True)
    print('  %-12s %-14s %6s %6s %5s %10s %5s  %7s %s' % ('product', 'mu (Munich)', 'n_ours', 'n_mun', 'mism', 'max dev', 'sign', 'secs', 'status'), flush=True)
    res = []
    skip = [x for x in os.environ.get('G7_MFSKIP', '').split(',') if x]      # e.g. '8:0000,8:2110' (memory: run those alone)
    only = os.environ.get('G7_MFONLY', '')
    for lam, lamp, mu in mf_list():
        S = lam[0]
        mm = munich_label(mu)
        key = '%d:%s' % (S, ''.join(map(str, mm)))
        if key in skip or (only and key not in only.split(',')):
            continue
        fname = 'hw_%d_%s.txt' % (S, '_'.join(map(str, mm)))
        r = compare_mf(lam, lamp, mu, fname)
        res.append((lam, mu, r))
        if 'dev' in r:
            print('  S=%-10d %-14s %6d %6d %5d %10.2e %5s  %7.2f %s' % (S, _s(mm), r['n_o'], r['n_m'], r['mism'], r['dev'], '%+d' % r['sigma'], r['secs'], r['status']), flush=True)
        else:
            print('  S=%-10d %-14s %s (%.0f s)' % (S, _s(mm), r['status'], r['secs']), flush=True)
        fails += not r['status'].startswith('MATCH')
    n = sum(1 for _, _, r in res if r['status'].startswith('MATCH'))
    print('  2a: %d MATCH of %d; max deviation over matches %.2e' % (n, len(res), max([r['dev'] for _, _, r in res if 'dev' in r] + [0])), flush=True)
    return res


BIG = [(9, 7, 5, 3), (10, 7, 5, 2), (10, 7, 4, 3), (9, 8, 5, 2), (10, 8, 4, 2), (10, 6, 5, 3)]


def run_6420():
    lam = (6, 4, 2, 0)
    print('\n== 2b: [6,4,2,0] (x) [6,4,2,0]: six largest-m sectors from the repaired profiler (HW_MODE=eig, prints all m alpha copies, 12 decimals)', flush=True)
    out = []
    for mu in BIG:
        m = gr.racah_speiser(lam, lam)[mu]
        t0 = time.perf_counter()
        r = compare_sector(lam, lam, mu, source='bench')
        out.append((lam, lam, mu, r))
        if 'sv' not in r:
            print('  mu=%s m=%d %s' % (mu, m, r['status']), flush=True)
            continue
        print('  mu=%s (Munich %s) m=%d n_W=%d  maxdev|1-sv|=%.2e  orthM=%.1e resid=%.1e  supp ours/Munich=%d/%d  signed-perm=%s  %s  (munich %.2f s, total %.1f s)'
              % (mu, _s(munich_label(mu)), m, gr.nW_of(lam, lam, mu), r['maxdev'], r['orth_munich'], r['resid'], r['n_ours'], r['n_mun'], r['signed_perm'], r['status'], r['secs'], time.perf_counter() - t0), flush=True)
        print('    sv = %s' % np.array2string(r['sv'], precision=12), flush=True)
    return out


def run_2b(products=None, verbose=True):
    products = products or PRODUCTS
    allres = []
    for lam, lamp in products:
        irr = sorted(gr.racah_speiser(lam, lamp).items(), reverse=True)
        print('\n== 2b: %s (x) %s : %d irreps, sum m = %d, max m = %d' % (lam, lamp, len(irr), sum(m for _, m in irr), max(m for _, m in irr)), flush=True)
        for mu, m in irr:
            r = compare_sector(lam, lamp, mu)
            allres.append((lam, lamp, mu, r))
            if 'sv' not in r:
                print('  mu=%s (Munich %s) m=%d  %s  (%.0f s)' % (mu, _s(munich_label(mu)), m, r['status'], r['secs']), flush=True)
                continue
            print('  mu=%s (Munich %s) m=%d  maxdev|1-sv|=%.2e  orthM=%.1e resid=%.1e  supp ours/Munich=%d/%d  signed-perm=%s  %s  (%.1f s)'
                  % (mu, _s(munich_label(mu)), m, r['maxdev'], r['orth_munich'], r['resid'], r['n_ours'], r['n_mun'], r['signed_perm'], r['status'], r['secs']), flush=True)
            if m >= 2:
                print('    sv = %s' % np.array2string(r['sv'], precision=12), flush=True)
                print('    Q (rows = our canonical copies k, cols = Munich alpha):\n' + fmt_Q(r['Q']), flush=True)
    ok = [r for *_, r in allres if r.get('status') == 'OK']
    sectors = len(allres)
    print('\n  2b summary: %d sectors, %d OK, sum m = %d; worst max|1-sv| = %.2e; worst Munich orthonormality = %.2e'
          % (sectors, len(ok), sum(r['m'] for *_, r in allres), max([r['maxdev'] for r in ok] + [0]), max([r['orth_munich'] for r in ok] + [0])), flush=True)
    return allres


# ---------------------------------------------------------------------------------------------------------------- selftest
def synth_text(lam, lamp, mu, vec, extra_rows=True):
    """text as hwcgc.sh would print it for a Munich vector {(Q,Q'): v}"""
    dimmu = gr.weyl_dim(mu)
    lines = ['Highest-weight-state CGCs of ...', 'Q(M)\tQ(M\')\tQ(M\'\')\tCGC']
    for (a, b), v in sorted(vec.items()):
        lines.append('%d\t%d\t%d\t%.15f' % (a, b, dimmu, v))
    return '\n'.join(lines) + '\n'


def selftest():
    rng = np.random.default_rng(5)
    cases = [((2, 1, 0, 0), (2, 1, 0, 0)), ((3, 2, 1, 0), (3, 2, 1, 0)), ((2, 2, 1, 0), (2, 1, 1, 0))]
    nacc = nrej = ntot_acc = ntot_rej = 0
    for lam, lamp in cases:
        for mu, m in sorted(gr.racah_speiser(lam, lamp).items(), reverse=True):
            prv = se.prv_space(lam, lamp, mu)
            cps = so.copies(lam, lamp, mu, dps=30, prv=prv)
            hw = se.hw_from_slice(lam, lamp, mu, so.copy_float_slices(prv, cps))
            O = [se.to_munich(lam, lamp, hw, k) for k in range(m)]
            keys = {k: i for i, k in enumerate(sorted(set().union(*O)))}
            Om = np.stack([sparse_from_rows(o, keys)[0] for o in O], axis=1)
            Qr, _ = np.linalg.qr(rng.normal(size=(m, m)))
            Mm = Om @ Qr                         # Munich copies = random orthogonal rotation of ours
            inv = {i: k for k, i in keys.items()}
            texts = [synth_text(lam, lamp, mu, {inv[i]: Mm[i, a] for i in range(len(keys)) if abs(Mm[i, a]) > 1e-14}) for a in range(m)]
            r = compare_sector(lam, lamp, mu, alpha_text=texts)
            ntot_acc += 1
            nacc += r['status'] == 'OK'
            # corrupted: perturb one component of the first copy
            bad = [dict(hw_rows(t, gr.weyl_dim(mu))) for t in texts]
            k0 = max(bad[0], key=lambda k: abs(bad[0][k]))
            bad[0][k0] *= 1.001
            btexts = [synth_text(lam, lamp, mu, b) for b in bad]
            r2 = compare_sector(lam, lamp, mu, alpha_text=btexts)
            ntot_rej += 1
            nrej += r2['status'] != 'OK'
            # corrupted: for m >= 2 replace the second copy by a vector leaving the span (mix in a foreign direction)
            if m >= 2:
                bad = [dict(hw_rows(t, gr.weyl_dim(mu))) for t in texts]
                a0, a1 = bad[0], bad[1]
                for k in a1:
                    a1[k] = a1[k]
                k1 = max(a0, key=lambda k: abs(a0[k]))
                a1[k1] = a1.get(k1, 0) + 0.3        # leaves the span and breaks orthonormality
                r3 = compare_sector(lam, lamp, mu, alpha_text=[synth_text(lam, lamp, mu, b) for b in bad])
                ntot_rej += 1
                nrej += r3['status'] != 'OK'
    print('munich_bridge selftest: %d/%d rotated-copy outputs accepted, %d/%d corrupted outputs rejected' % (nacc, ntot_acc, nrej, ntot_rej))
    return nacc == ntot_acc and nrej == ntot_rej


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 0
    cmd = argv[1]
    if cmd == 'selftest':
        return 0 if selftest() else 1
    if cmd == 'compare':
        which = argv[2] if len(argv) > 2 else 'all'
        os.makedirs(OUT_DIR, exist_ok=True)
        if which in ('2a', 'all'):
            run_2a()
        if which in ('2b', 'all'):
            run_2b()
        if which in ('6420', 'all'):
            run_6420()
        return 0
    if cmd == 'sector':
        lam, lamp, mu = [tuple(int(x) for x in a.split(',')) for a in argv[2:5]]
        run_2b_one = compare_sector(lam, lamp, mu)
        print({k: v for k, v in run_2b_one.items() if k != 'Q'})
        print(fmt_Q(run_2b_one['Q']))
        return 0
    print(__doc__)
    return 1


if __name__ == '__main__':
    sys.exit(main(sys.argv))
