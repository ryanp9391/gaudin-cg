"""G7 3d: our t_a / t_b (slice_ops_gl4.bench_one, min of reps cold runs, single thread) on the multiplicity products.  Run from Paul/Python/gaudin_gl4 with the venv."""
import sys, time
import glN_reps as gr, slice_ops_gl4 as so
PRODS = [(2,1,1,0),(3,2,1,0),(4,2,0,0),(4,2,1,0),(4,3,1,0),(6,3,0,0),(6,4,2,0),(8,4,0,0),(9,6,3,0)]
reps = {(9,6,3,0): 2, (8,4,0,0): 2}
print('# product  dimV  irreps  copies(sum m)  max m  max nW | t_a[s]  t_b[s]  reps | worst HW resid')
for lam in PRODS:
    t0 = time.perf_counter()
    T = so.bench_one(lam, lam, reps=reps.get(lam, 3 if gr.weyl_dim(lam) < 500 else 2))
    print('%s %d %d %d %d %d | %.4f %.4f %d | %.1e   (wall %.0f s)' % (lam, T['dim'], T['nirrep'], T['ncopies'], T['maxm'], T['maxnW'], T['ta'], T['tb'], T['reps'], T['resid'], time.perf_counter()-t0), flush=True)
