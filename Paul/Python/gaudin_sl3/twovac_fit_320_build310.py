"""helper: build [3,1,0] zeta points (for the comparison pipeline) in parallel"""
import sys, time, os, pickle
import twovac_fit_320 as t
from concurrent.futures import ProcessPoolExecutor

def job(a):
    tag, z = a
    t.set_lam((3, 1, 0))
    try:
        t0 = time.time()
        pt = t.build_point(z)
        with open(t.pt_path(tag, t.DPS, t.THETA), 'wb') as f:
            pickle.dump(pt, f)
        return tag, 'ok %.0fs' % (time.time() - t0)
    except Exception as e:
        return tag, 'FAIL %r' % (e,)

if __name__ == '__main__':
    t.set_lam((3, 1, 0))
    tags = sys.argv[1:]
    todo = [(g, t.ZETAS[g]) for g in tags if not os.path.exists(t.pt_path(g, t.DPS, t.THETA))]
    with ProcessPoolExecutor(max_workers=6) as ex:
        for tag, msg in ex.map(job, todo):
            print(tag, msg, flush=True)
