"""helper: extra theta variants (single zeta point P)"""
import sys, time, os, pickle
import twovac_fit_320 as t
from concurrent.futures import ProcessPoolExecutor

def job(th):
    try:
        t0 = time.time()
        pt = t.build_point(t.ZETAS['P'], t.DPS, th)
        with open(t.pt_path('P', t.DPS, th), 'wb') as f:
            pickle.dump(pt, f)
        return th, 'ok %.0fs' % (time.time() - t0)
    except Exception as e:
        return th, 'FAIL %r' % (e,)

if __name__ == '__main__':
    todo = [th for th in t.THETA_EXTRA if not os.path.exists(t.pt_path('P', t.DPS, th))]
    with ProcessPoolExecutor(max_workers=8) as ex:
        for r in ex.map(job, todo):
            print(*r, flush=True)
