"""test_cg_fast_gl4.py -- acceptance tests of cg_fast_gl4.py (WP-G9).  Run:  python -X utf8 test_cg_fast_gl4.py   (or pytest); total < 90 s.

Each test compares the float64 engine with the Windows-verified exact engine (G3 `copies`, 30 digits here for speed), with the exact weight-space matrices, or with the cached Munich output of G7.
"""
import os, sys, time, unittest
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import scipy.sparse as sps
import glN_reps as gr
import slice_engine_gl4 as se
import slice_ops_gl4 as so
import cg_fast_gl4 as cg

MUNICH_DIR = os.path.normpath(os.path.join(HERE, '..', '..', 'Cpp', 'AlexCGC', 'results', 'su4-same-machine'))


def exact(lam, lamp, mu, dps=30):
    cps = so.copies(lam, lamp, mu, dps=dps)
    F = so.copy_float_slices(cps.prv, cps)
    trip = np.array([[float(x) for x in cp['ev']] for cp in cps])
    return F, trip


def slice_dev(sec, F):
    return float((np.abs(sec['slices'] - F).max(axis=0) / np.abs(F).max(axis=0)).max())


def trip_dev(sec, trip):
    return float((np.abs(sec['trip'] - trip) / np.maximum(1.0, np.abs(trip))).max())


class TestFastCG(unittest.TestCase):

    def test_01_slices_hw_orthogonality_multiplicity(self):
        """T1 on [4,2,1,0]^2 (56 copies): every sector: #physical = m (no fallback), slices and triples = G3's copies, HW residual, G(x)G orthogonality of the copies"""
        lam = (4, 2, 1, 0)
        secs = cg.hw_vectors(lam, lam, strict=True)
        rs = gr.racah_speiser(lam, lam)
        self.assertEqual(set(secs), set(rs))
        worst = 0.0
        for mu, sec in secs.items():
            self.assertEqual(sec['m'], rs[mu])
            self.assertIsNone(sec['fallback'], mu)
            self.assertEqual(sec['info']['nphys'], rs[mu], mu)
            F, trip = exact(lam, lam, mu)
            self.assertLess(slice_dev(sec, F), 1e-11, mu)
            self.assertLess(trip_dev(sec, trip), 1e-11, mu)
            res, wd, off = cg.hw_checks(lam, lam, sec)
            self.assertLess(res, 1e-11, mu)
            self.assertEqual(wd, 0)
            self.assertLess(off, 1e-11, mu)
            worst = max(worst, slice_dev(sec, F))
        self.assertLess(worst, 1e-11)

    def test_02_route_A_equals_route_B(self):
        """route A (closed forms of G8d, float) == route B (reduced words) in the orthonormal basis, 1e-12 relative to max |M| (to 1 if M = 0)"""
        import closed_forms_gl4 as cf
        if not os.path.exists(cf.CACHE_PATH):
            self.skipTest('closed_forms_gl4_qcache.pkl absent: its one-off symbolic build takes ~5 min (run `python -X utf8 closed_forms_gl4.py mats` or the validate stage once)')
        lam = (3, 2, 1, 0)
        S = cg.Setup.of(lam, lam)
        for mu in sorted(gr.racah_speiser(lam, lam)):
            W, pos, B, leak = cg.operators_B(S, mu)
            _, _, A = cg.operators_A(S, mu)
            for n in A:
                self.assertLess(float(np.abs(A[n] - B[n]).max() / max(np.abs(B[n]).max(), 1.0)), 1e-12, (mu, n))

    def test_03_operators_vs_exact_weight_space_mats(self):
        """exact twin: route B on W1 == G3's exact weight_space_mats conjugated by sqrt(g), and the exact matrices leave W1 invariant (leak exactly 0)"""
        for lam, lamp in (((4, 2, 1, 0),) * 2, ((3, 2, 1, 0), (4, 2, 1, 0))):
            for mu in sorted(gr.racah_speiser(lam, lamp), reverse=True)[::3]:
                d = cg.mats_vs_exact(lam, lamp, mu, 'B')
                for n, (dev, leak) in d.items():
                    self.assertLess(dev, 1e-12, (lam, lamp, mu, n))
                    self.assertEqual(leak, 0.0)
                self.assertTrue(cg.check_cut(lam, lamp, mu)[0], mu)

    def test_04_string_test_gap(self):
        """physical eigenvectors violate no string (s <= 1e-9), spurious ones violate by >= 0.5: gap >= 1e8 over [4,3,1,0]^2"""
        lam = (4, 3, 1, 0)
        sp, ph = [], []
        for mu in gr.racah_speiser(lam, lam):
            sec = cg.physical_slices(lam, lam, mu, strict=True)
            ph.append(sec['s_phys'])
            if np.isfinite(sec['s_spur']):
                sp.append(sec['s_spur'])
        self.assertTrue(sp)
        self.assertLess(max(ph), 1e-9)
        self.assertGreater(min(sp), 0.5)
        self.assertGreater(min(sp) / max(max(ph), 1e-300), 1e8)

    def test_05_munich_multiplicity_free_match(self):
        """cached Munich hwcgc.sh outputs of G7: [S,S,0,0]^2 irreps match after one global sign to 1e-9 (S = 2, 3)"""
        n = 0
        for S in (2, 3):
            lam = (S, S, 0, 0)
            secs = cg.hw_vectors(lam, lam)
            for mu, sec in secs.items():
                path = os.path.join(MUNICH_DIR, 'hw_%d_%s.txt' % (S, '_'.join(map(str, [x - mu[3] for x in mu]))))
                if sec['m'] != 1 or not os.path.exists(path):
                    continue
                ok, dev, nsup = cg.compare_mf_file(lam, lam, mu, sec, path)
                self.assertTrue(ok, (S, mu, dev, nsup))
                n += 1
        self.assertGreaterEqual(n, 3, 'no cached Munich file found (G7 cache missing?)')

    def test_06_lambda_ne_lambdap(self):
        """[3,2,1,0] x [4,2,1,0] and the reversed order: all sectors vs the exact copies"""
        for lam, lamp in (((3, 2, 1, 0), (4, 2, 1, 0)), ((4, 2, 1, 0), (3, 2, 1, 0))):
            secs = cg.hw_vectors(lam, lamp)
            for mu, sec in secs.items():
                self.assertIsNone(sec['fallback'])
                F, trip = exact(lam, lamp, mu)
                self.assertLess(slice_dev(sec, F), 1e-11, (lam, lamp, mu))
                self.assertLess(trip_dev(sec, trip), 1e-11, (lam, lamp, mu))
                self.assertLess(cg.hw_checks(lam, lamp, sec)[0], 1e-11)

    def test_07_lambda4_positive(self):
        """lambda_4 > 0: [4,3,2,1] x [3,2,1,1] (all sectors) and the sectors of [5,3,2,1]^2 with m >= 3"""
        lam, lamp = (4, 3, 2, 1), (3, 2, 1, 1)
        for mu, sec in cg.hw_vectors(lam, lamp).items():
            F, trip = exact(lam, lamp, mu)
            self.assertLess(slice_dev(sec, F), 1e-11, mu)
            self.assertLess(trip_dev(sec, trip), 1e-11, mu)
        lam = (5, 3, 2, 1)
        rs = gr.racah_speiser(lam, lam)
        mus = [mu for mu, m in rs.items() if m >= 3]
        self.assertTrue(mus)
        for mu, sec in cg.hw_vectors(lam, lam, mus=mus).items():
            F, trip = exact(lam, lam, mu)
            self.assertLess(slice_dev(sec, F), 1e-11, mu)
            self.assertIsNone(sec['fallback'])

    def test_08_degenerate_sector_6420(self):
        """(9,7,5,3) of [6,4,2,0]^2 (m = 12; V3 and V4a do not separate the copies): separated by the triple, slices and triples = exact, no fallback"""
        lam = (6, 4, 2, 0)
        mu = (9, 7, 5, 3)
        sec = cg.physical_slices(lam, lam, mu, strict=True)
        self.assertEqual(sec['m'], 12)
        self.assertIsNone(sec['fallback'])
        d34 = min(max(abs(a[0] - b[0]), abs(a[1] - b[1])) for i, a in enumerate(sec['trip']) for b in sec['trip'][:i])
        self.assertLess(d34, 1e-8)                                  # (V3, V4a) alone do not separate: the sector is the degenerate one of G8d
        F, trip = exact(lam, lam, mu)
        self.assertLess(slice_dev(sec, F), 1e-11)
        self.assertLess(trip_dev(sec, trip), 1e-11)

    def test_09_fallback_path(self):
        """the logged exact fallback (G3 copies) agrees with the float path, and a sector whose physical count is wrong is routed to it"""
        lam = (4, 2, 1, 0)
        mu = max(gr.racah_speiser(lam, lam), key=lambda k: gr.racah_speiser(lam, lam)[k])
        sec = cg.physical_slices(lam, lam, mu)
        prv, F, trip, info = cg._fallback_sector(lam, lam, mu, dps=30)
        self.assertEqual(F.shape, sec['slices'].shape)
        self.assertLess(slice_dev(sec, F), 1e-11)
        self.assertLess(trip_dev(sec, trip), 1e-11)
        old = cg._slices_float
        try:
            cg._slices_float = lambda *a, **k: (False, dict(w1=0, ncomplex=0, s_phys=0.0, s_spur=np.inf, nphys=0, ritz='forced'))
            sec2 = cg.physical_slices(lam, lam, mu)
        finally:
            cg._slices_float = old
        self.assertEqual(sec2['fallback'], 'forced')
        self.assertLess(slice_dev(sec2, F), 1e-11)

    def test_10_whole_product_6420(self):
        """[6,4,2,0]^2 end to end (76 irreps, 297 copies): no fallback, slices / triples = G3's exact copies (30 digits), HW residual, orthogonality, Racah-Speiser multiplicities"""
        lam = (6, 4, 2, 0)
        secs = cg.hw_vectors(lam, lam)
        rs = gr.racah_speiser(lam, lam)
        self.assertEqual(sum(s['m'] for s in secs.values()), 297)
        worst = 0.0
        for mu, sec in secs.items():
            self.assertEqual(sec['m'], rs[mu])
            self.assertIsNone(sec['fallback'], mu)
            F, trip = exact(lam, lam, mu)
            worst = max(worst, slice_dev(sec, F))
            self.assertLess(trip_dev(sec, trip), 1e-11, mu)
            res, wd, off = cg.hw_checks(lam, lam, sec)
            self.assertLess(res, 1e-11, mu)
            self.assertLess(off, 1e-11, mu)
        self.assertLess(worst, 1e-11)

    def test_11_descendants_unitarity(self):
        """S4: the descendants of all HW vectors of [2,1,0,0]^2 form an orthogonal matrix (full CG table): unitarity <= 1e-12"""
        if not hasattr(cg, 'cg_table'):
            self.skipTest('descendants not implemented')
        lam = (2, 1, 0, 0)
        C = cg.cg_table(lam, lam)
        n = C.shape[0]
        self.assertEqual(C.shape[0], C.shape[1])
        self.assertLess(float(abs(C @ C.T - sps.identity(n)).max()), 1e-12)


if __name__ == '__main__':
    t0 = time.perf_counter()
    r = unittest.main(exit=False, verbosity=2)
    print('total %.1f s' % (time.perf_counter() - t0))
    sys.exit(0 if r.result.wasSuccessful() else 1)
