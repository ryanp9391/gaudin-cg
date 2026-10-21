"""Write irreps_a<a>.txt (GT mu, m) for [2a,a,0]x[2a,a,0] from slice_engine.racah_speiser, cross-checked against jacobi_cg.irreps (closed form).
usage: python -I gen_irreps.py a [a ...]"""
import sys
sys.path.insert(0, '/home/paul/git_projects/gaudin-CG/Paul/Python/gaudin_sl3')
sys.path.insert(0, '/home/paul/git_projects/gaudin-CG/Paul/Python')
import slice_engine as se, jacobi_cg as jc
for a in map(int, sys.argv[1:]):
    lam = (2*a, a, 0)
    A = se.racah_speiser(lam) if a <= 20 else None
    B = jc.irreps(lam)
    if A is not None:
        assert dict(A) == dict(B), 'MISMATCH'
    with open('irreps_a%d.txt' % a, 'w') as f:
        f.write('# GT mu1 mu2 mu3, m ; [%d,%d,0]x[%d,%d,0]\n' % (2*a, a, 2*a, a))
        for mu, m in sorted(B.items(), reverse=True):
            f.write('%d %d %d %d\n' % (*mu, m))
    print(a, 'irreps', len(B), 'copies', sum(B.values()), 'RS==closed form:', A is not None)
