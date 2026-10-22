"""G7 3c: Munich hw-solve times on multiplicity products, original (cgc_hw_bench) vs repaired (cgc_hw_prof_ob HW_MODE=eig).  Run from anywhere:
   ../.venv/bin/python run_products.py   (single thread, LC_ALL=C; memory guard n_W dim^2 8 B <= 10 GB for the original)"""
import os, re, subprocess, sys, time
sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', 'Python', 'gaudin_gl4')))
import glN_reps as gr
CGC = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
R = os.path.join(CGC, 'results', 'su4-same-machine')
ENV = dict(os.environ, LC_ALL='C', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
PRODS = [(2, 1, 1, 0), (3, 2, 1, 0), (4, 2, 0, 0), (4, 2, 1, 0), (4, 3, 1, 0), (6, 3, 0, 0), (6, 4, 2, 0), (8, 4, 0, 0), (9, 6, 3, 0)]
NO_ORIG = [(6, 4, 2, 0), (9, 6, 3, 0)]     # brief: repaired variants only
NWMAX_PROF = 6000                         # n_W^2 * 8 B (A A^T) and n_W^3 time guard for eig/compact
NWMAX_COMPACT = 1500
def reps_for(nW):
    return 5 if nW < 300 else (3 if nW < 1000 else 1)
LIMIT = 10e9

def run(cmd, mode=None):
    env = dict(ENV)
    if mode:
        env['HW_MODE'] = mode
    t0 = time.perf_counter()
    r = subprocess.run(['/usr/bin/time', '-v', 'timeout', '1800'] + cmd, cwd=CGC, capture_output=True, text=True, env=env)
    m = re.search(r'hw_solve_s=(\S+)', r.stdout)
    rss = re.search(r'Maximum resident set size \(kbytes\): (\d+)', r.stderr)
    return (float(m.group(1)) if m else None), (int(rss.group(1)) if rss else None), time.perf_counter() - t0

def main(prods):
    for lam in prods:
        dim = gr.weyl_dim(lam)
        dec = subprocess.run(['./cgc_hw_bench', 'decompose', '4'] + [str(x) for x in lam] * 2, cwd=CGC, capture_output=True, text=True, env=ENV).stdout
        irr = []
        for line in dec.splitlines():
            f = line.split()
            if f and f[0] == 'MU':
                irr.append((tuple(int(x) for x in f[1:5]), int(f[6]), int(f[8])))
        out = open(os.path.join(R, 'products_%s.txt' % ''.join(map(str, lam))), 'w')
        print('# product %s^2: %d irreps, sum m = %d, max m = %d; cmds: ./cgc_hw_bench 4 L L MU reps machine (reps 5 if nW<300, 3 if <1000, else 1) | HW_MODE=eig ./cgc_hw_prof_ob ... ; load %s' % (lam, len(irr), sum(m for _, m, _ in irr), max(m for _, m, _ in irr), open('/proc/loadavg').read().strip()), file=out, flush=True)
        print('# mu_munich  m  dimmu  nW  reps  mem_orig_GB  orig_hw_s  orig_rss_MB  eig_hw_s  eig_rss_MB  compact_hw_s  compact_rss_MB   (NA = not run: memory/time guard)', file=out, flush=True)
        for mu, m, dmu in irr:
            ours = gr.munich_label_to_ours(lam, lam, mu)
            nW = gr.nW_of(lam, lam, ours)
            mem = nW * dim * dim * 8
            rp = reps_for(nW)
            args = [str(x) for x in lam + lam + mu] + [str(rp), 'machine']
            to = ro = te = re_ = tc = rc = None
            if mem <= LIMIT and lam not in NO_ORIG:
                to, ro, _ = run(['./cgc_hw_bench', '4'] + args)
            if nW <= NWMAX_PROF:
                te, re_, _ = run(['./cgc_hw_prof_ob', '4'] + args, mode='eig')
            if nW <= NWMAX_COMPACT:
                tc, rc, _ = run(['./cgc_hw_prof_ob', '4'] + args, mode='compact')
            f = lambda x: 'NA' if x is None else '%.6g' % x
            g = lambda x: 'NA' if x is None else str(x // 1024)
            print('%s %d %d %d %d %.3f %s %s %s %s %s %s' % (' '.join(map(str, mu)), m, dmu, nW, rp, mem / 1e9, f(to), g(ro), f(te), g(re_), f(tc), g(rc)), file=out, flush=True)
        out.close()

if __name__ == '__main__':
    main([tuple(int(c) for c in a) for a in sys.argv[1:]] or PRODS)
