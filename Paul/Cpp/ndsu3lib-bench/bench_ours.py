"""Wrapper: time our Python engines at one a, sequentially.  usage: python -I bench_ours.py {se|jc} a reps
se: slice_engine.bench_one(a,1) (cold caches) -> t_a, t_b, #irreps/copies/entries, HW residual, float-eigen residual
jc: jacobi_cg.stage_bench([a], verbose=False, with_R=False) -> t_float and float64-vs-mpmath(30 digits) deviation, forward / twisted recurrence"""
import sys
sys.path.insert(0, '/home/paul/git_projects/gaudin-CG/Paul/Python/gaudin_sl3')
sys.path.insert(0, '/home/paul/git_projects/gaudin-CG/Paul/Python')
kind, a, reps = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
open('/proc/loadavg')
for i in range(reps):
    load = open('/proc/loadavg').read().split()[0]
    if kind == 'se':
        import slice_engine as se
        T = se.bench_one(a, 1)
        print('se a=%d rep %d load %s  t_a %.4f t_b %.4f  irreps %d copies %d entries %d maxw %d  hw_resid %.2e eigres %.2e mingap %.2e smin %.2e' %
              (a, i, load, T['ta'], T['tb'], T['nirrep'], T['ncopies'], T['nent'], T['maxw'], T['resid'], T['eigres'], T['mingap'], T['smin']), flush=True)
    else:
        import jacobi_cg as jc
        r = jc.stage_bench([a], verbose=False, with_R=False)[0]
        print('jc a=%d rep %d load %s  t_float %.4f (irr %.4f chain %.4f eig_f %.4f rec_f %.4f rec_t %.4f)  irreps %d copies %d  float64 dev vs mp30: forward %.2e twisted %.2e  fdevV %.2e' %
              (a, i, load, r['t_float'], r['t_irr'], r['t_chain'], r['t_eig_f'], r['t_rec_f'], r['t_rec_t'], r['nirr'], r['ncop'], r['fdev'], r['fdev_t'], r['fdevV']), flush=True)
