"""helper: build zeta points at theta variants in parallel"""
import sys, time, os, pickle
import twovac_fit_320 as t
from concurrent.futures import ProcessPoolExecutor

def job(a):
    tag, th = a
    try:
        t0 = time.time()
        pt = t.build_point(t.ZETAS[tag], t.DPS, th)
        with open(t.pt_path(tag, t.DPS, th), 'wb') as f:
            pickle.dump(pt, f)
        return tag, th, 'ok %.0fs' % (time.time() - t0)
    except Exception as e:
        return tag, th, 'FAIL %r' % (e,)

if __name__ == '__main__':
    tags = sys.argv[1:]
    todo = [(g, th) for th in t.tf.THETA_VARS for g in tags if not os.path.exists(t.pt_path(g, t.DPS, th))]
    with ProcessPoolExecutor(max_workers=14) as ex:
        for r in ex.map(job, todo):
            print(*r, flush=True)
