"""test_structure_gl4.py -- quick regression tests for structure_gl4.py (WP-G4); pytest-style functions AND runnable as `python -X utf8 test_structure_gl4.py` (a few seconds).
Exact (Fractions) except the explicitly marked 50-digit (mpmath) checks.  The heavy runs (all products, all displacement classes of V4a/V4b, the rule/chain/poly tables) are the CLI stages."""
import os, sys, time, random
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mpmath as mp
import glN_reps as gr
import slice_ops_gl4 as so
import structure_gl4 as sg


def test_stencil_v3():
    st = sg.structural_stencil('V3')
    assert len(st) == 21 and all(tuple(-x for x in d) in st for d in st)
    prod = {(a, b, c) for a in (-1, 0, 1) for (b, c) in [(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]}      # row 2: 0 or +-one root of gl(2); row 3: 0 or one root of gl(3)
    assert set(st) == prod
    for lam, lamp in (((4, 3, 1, 0), (Fr(1, 3), Fr(5, 7), Fr(-2), Fr(4, 5))),):
        dat = sg.data_stencil(lam, [lamp], names=('V3',))['V3']
        assert set(dat) == prod, sorted(set(dat) ^ prod)
    # V4a / V4b: the structural bound has the |dm12| = 2 classes, the data (and the exact summation) do not
    assert len(sg.structural_stencil('V4a')) == 35
    return 'V3 stencil: 21 displacements = {0, +-(gl2 root)} x {0, 6 gl3 roots}, structural = data (generic lam\'); V4a/V4b structural bound 35'


def test_symbolic_class_heldout():
    """one fully factorising displacement class of V3 (dm12, dm13, dm23) = (1, 1, 0): identified symbolically, then compared with G3's exact matrices on the held-out rep [4,3,1,0]"""
    T = sg.sym_terms('V3')
    d = next(k for k in T if sg.project(k) == (1, 1, 0))
    rec = sg.closed_form_one('V3', d)
    assert rec['check'] and rec['qdeg'] == 0 and sum(m for _, _, m in rec['forms']) == 9
    CF = {d: rec}
    rng = random.Random(1)
    lam = (4, 3, 1, 0)
    st = sg.verify_closed('V3', lam, [lam, (Fr(1, 3), Fr(5, 7), Fr(-2), Fr(4, 5))], rng, which='F', CF=CF)
    assert st['bad'] == 0 and st['tests'] >= 30 and st['nonzero_tests'] >= 30 and st['zerodiv'] == 0, dict(st)
    nm = [sg.form_str(f) for f, e in rec['net'] if e > 0]
    assert '(m11 - m12 + q1 - q2)' in nm                      # the only lam'-dependence of the Delta m12 = +1 classes (through mu_1 - lam'_2 - ...)
    return 'V3 class (+1,+1,0): numerator of 9 linear factors verified on [4,3,1,0] (%d tests, %d nonzero entries, 0 mismatches)' % (st['tests'], st['nonzero_tests'])


def test_symbolic_nonfactorising_class():
    """the pure row-2 move (+1, 0, 0) of V3 is NOT a product of linear forms: 8 linear factors x a polynomial of degree 6 in (m13, m23, m33, l1..l4); extrapolation reproduces G3's entries on [3,2,1,0]"""
    T = sg.sym_terms('V3')
    d = next(k for k in T if sg.project(k) == (1, 0, 0))
    rec = sg.closed_form_one('V3', d)
    assert rec['check'] and rec['qdeg'] == 6 and sum(m for _, _, m in rec['forms']) == 8
    assert {sg.SYM_LABEL[i] for i, q in enumerate(rec['qaxis']) if q} == {'m13', 'm23', 'm33', 'l1', 'l2', 'l3', 'l4'}
    return 'V3 class (+1,0,0): 8 linear factors x Q (degree 6, support m13, m23, m33, l1..l4)'


def test_universal_formula():
    """the single product formula of the 12 factorising classes (M^V4a = (y_k + 3) M^V3, M^V4b = -(y_a - y_i + 2 y_k - q3 + 5) M^V3) against G3's matrices on held-out reps, dominant and generic lam'"""
    assert len(sg.f_classes()) == 12
    tot = 0
    for lam in ((4, 3, 1, 0), (5, 4, 2, 1)):
        st = sg.verify_universal(lam, [lam, (Fr(1, 3), Fr(5, 7), Fr(-2), Fr(4, 5))])
        assert st['bad'] == 0 and st['zerodiv'] == 0 and st['tests'] > 500 and st['nonzero'] > 400, dict(st)
        tot += st['tests']
    return 'universal formula for the 12 factorising classes of V3, V4a, V4b: %d exact tests on [4,3,1,0], [5,4,2,1], 0 mismatches' % tot


def test_weyl_sector_main():
    """the Weyl/BGG identity on one sector of [4,2,0,0]^2 with a triple occurring in FOUR translates (e, s1, s3, s1 s3): signed count = indicator"""
    lam = lamp = (4, 2, 0, 0)
    mu = (5, 5, 1, 1)
    r = sg.weyl_sector(lam, lamp, mu)
    assert r['m'] == 1 and not r['fails'] and r['phys_found'] == r['nphys'] == 1 and r['jok'] and r['mult_signed'] == 1
    lens = sorted(tuple(g['lens']) for g in r['groups'])
    assert (0, 1, 1, 2) in lens, lens
    assert r['depth_ok'][0] is False and r['depth_ok'][1] is False and r['depth_ok'][2] is True
    # one-step: the spurious triple violates two strings
    assert any(len(row['S']) == 2 for row in r['onestep'])
    # another sector, [3,2,1,0]^2 mu = (5,4,2,1) (m = 4, V3 degenerate): no spurious triple, identity only
    r2 = sg.weyl_sector((3, 2, 1, 0), (3, 2, 1, 0), (5, 4, 2, 1))
    assert not r2['fails'] and r2['m'] == 4 and len(r2['spaces']) == 1
    return 'Weyl/BGG: [4,2,0,0]^2 mu = (5,5,1,1): spurious triple in e, s1, s3, s1 s3 (signed 1 - 2 + 1 = 0), depth 2 needed; [3,2,1,0]^2 (5,4,2,1): identity only'


def test_weyl_sensitivity():
    """the test can fail: all signs +1, or the plain Weyl action instead of the dot action, break the identity on [3,2,1,0]^2 sectors that need a translate"""
    lam = lamp = (3, 2, 1, 0)
    ok_rule = sum(1 for mu, m in sg.sector_list(lam, lamp) if not sg.weyl_sector(lam, lamp, mu, onestep=False)['fails'])
    ok_nosign = sum(1 for mu, m in sg.sector_list(lam, lamp) if not sg.weyl_sector(lam, lamp, mu, onestep=False, variant='nosign')['fails'])
    ok_plain = sum(1 for mu, m in sg.sector_list(lam, lamp) if not sg.weyl_sector(lam, lamp, mu, onestep=False, variant='plain')['fails'])
    assert ok_rule == 16 and ok_nosign == 11 and ok_plain < 16, (ok_rule, ok_nosign, ok_plain)
    return 'sensitivity on [3,2,1,0]^2: dot action + signs 16/16, all signs +1 %d/16, plain Weyl action %d/16' % (ok_nosign, ok_plain)


def test_weyl_gl3():
    """the gl(3) sanity check of the same statement ([4,2,0]^2, all sectors) and the two-sided ordering rule of V (E12-type: bottom m, E23-type: top m)"""
    lam = lamp = (4, 2, 0)
    res, S = sg.weyl_product(lam, lamp, kind='all', fam=sg.FAM[3], verbose=False)
    assert S['ok'] == S['sectors'] == 13 and S['mult_ok'] == 13 and S['os_S_eq_T'] == S['os_tau'] == 21
    dirs = set()
    for mu, m in sg.sector_list(lam, lamp):
        sec = sg.ident_sector(lam, lamp, mu, sg.FAM[3])
        if sec['w'] > sec['m']:
            dirs.add(sg.sep_direction([mp.re(t[0]) for t in sec['js']['trip']], sec['phys']))
    assert dirs == {'top', 'bottom'}
    return 'gl(3) [4,2,0]^2: signed count = indicator on 13/13 sectors; physical = top m (E23-type) or bottom m (E12-type) of V'


def test_krylov_polynomials():
    """graded Krylov degrees on M_mu: copies of [3,2,1,0]^2, mu = (5,4,2,1) (m = 4): every coordinate cyclic, components affine in the triple, exact polynomials reproduce the 50-digit copies"""
    lam = lamp = (3, 2, 1, 0)
    mu = (5, 4, 2, 1)
    prv, md = so.triple_mats(lam, lamp, mu)
    mats = [md[n] for n in so.TRIPLE]
    kr = sg.left_krylov(mats, 0)
    assert kr['cyclic'] and max(kr['deg']) == 1
    cps = so.copies(lam, lamp, mu, dps=50, prv=prv, mats=md)
    old = mp.mp.dps
    mp.mp.dps = 60
    try:
        err = max(abs(sg.poly_eval(kr['poly'][y], cp['ev']) * cp['coords'][0] - cp['coords'][y]) for cp in cps for y in range(4))
    finally:
        mp.mp.dps = old
    assert err < mp.mpf(10) ** -40, err
    assert sg.d_gen(4) == 1 and sg.d_gen(5) == 2 and sg.d_gen(10) == 2 and sg.d_gen(11) == 3
    return 'Krylov: all 4 coordinates of (5,4,2,1) of [3,2,1,0]^2 cyclic, affine polynomials in the triple, error %s' % mp.nstr(err, 3)


def test_graph_stats():
    M = [[Fr(1), Fr(2), Fr(0)], [Fr(3), Fr(1), Fr(5)], [Fr(0), Fr(7), Fr(1)]]
    g = sg.graph_stats(M)
    assert g['paths'] and g['ncomp'] == 1 and g['bw'] == 1 and g['maxdeg'] == 2
    M2 = [[Fr(0), Fr(1), Fr(1)], [Fr(1), Fr(0), Fr(1)], [Fr(1), Fr(1), Fr(0)]]
    assert not sg.graph_stats(M2)['paths']
    return 'graph statistics (path / triangle)'


def run_all():
    tests = [test_stencil_v3, test_symbolic_class_heldout, test_symbolic_nonfactorising_class, test_universal_formula, test_weyl_sector_main, test_weyl_sensitivity, test_weyl_gl3, test_krylov_polynomials, test_graph_stats]
    nfail = 0
    t00 = time.perf_counter()
    for t in tests:
        t0 = time.perf_counter()
        try:
            msg = t()
            print('PASS  %-40s %.1fs  %s' % (t.__name__, time.perf_counter() - t0, msg), flush=True)
        except Exception as e:
            nfail += 1
            import traceback
            print('FAIL  %-40s %s' % (t.__name__, repr(e)), flush=True)
            traceback.print_exc()
    print('%d tests, %d failed (%.1f s)' % (len(tests), nfail, time.perf_counter() - t00))
    return nfail


if __name__ == '__main__':
    sys.exit(1 if run_all() else 0)
