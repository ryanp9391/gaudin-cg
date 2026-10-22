"""test_slice_ops_gl4.py -- quick regression tests for slice_ops_gl4.py (pytest-style functions AND runnable as `python -X utf8 test_slice_ops_gl4.py`).
Everything is exact (Fractions) except the explicitly marked mpmath / float checks; the heavy runs (all sectors of all products, benchmark, reference files) are the CLI stages of the module."""
import os, sys, time, itertools, random
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mpmath as mp
import glN_reps as gr
import slice_engine_gl4 as se
import gaudin_gl4_ops as go
import slice_ops_gl4 as so


def test_normal_order_basic():
    assert so.normal_order(((1, 2), (2, 1))) == {((2, 1), (1, 2)): 1, ((1, 1),): 1, ((2, 2),): -1}      # e12 e21 = e21 e12 + e11 - e22
    assert so.normal_order(((2, 1), (1, 2))) == {((2, 1), (1, 2)): 1}
    assert so.normal_order(((1, 1), (2, 2))) == {((1, 1), (2, 2)): 1}
    for w in [((1, 3), (3, 1), (2, 2)), ((1, 4), (4, 1), (1, 2), (2, 1))]:
        for u in so.normal_order(w):
            assert so.is_normal_ordered(u)
    return 'e12 e21 = e21 e12 + e11 - e22 and the class order lowering < Cartan < raising'


def test_normal_order_random():
    r = so.test_normal_order(nrandom=60, seed=5, verbose=False)
    assert sum(r['bad'].values()) == 0 and sum(r['bad_s41'].values()) == 0 and r['example_ok']
    return '%d random + exhaustive words and the %d length-3 site-2 words of S_{4,1}: exact matrices on [2,1,0,0], [1,1,0,0] and the PBW normal form agree' % (r['n_words'], r['n_s41'])


def test_n3_convention():
    r = so.gl3_convention_test(lams=((2, 1, 0), (3, 1, 0)), verbose=False)
    assert r['bad_conv'] == 0 and r['bad_rel'] == 0 and r['direct_bad'] == 0
    return 'N = 3: reduced t S_{3,2} = slice_engine.vred_apply + eta_1/2 on %d PRV vectors; t^2 S_{3,1} + 2 tS_{3,2} = c 1 (integer c) on %d sectors; direct route %d checks' % (r['vec'], r['nsec'], r['direct_n'])


def test_hand_rule_matches_engine():
    lam = lamp = (2, 1, 0, 0)
    n, bad = so.closed_form_check(lam, lamp, verbose=False)
    assert bad == 0
    # one-factor terms: the hand rule is literally the engine's output
    for name in ('V3', 'V4a'):
        T = so.op_terms(name, 4)
        assert so.hand_reduce_single(T, lamp) == so.reduced_op(name, lamp).red
    return 'V_red,3 / V_red,4a: engine == one-line rule == hand-reduced closed forms (sum_I D(2,1)_I + eta, D(3,1) + V3 - eta/3) on all PRV vectors of [2,1,0,0]^2'


def test_gl4_direct_route_one_sector():
    lam = lamp = (3, 2, 1, 0)
    mu = (5, 4, 2, 1)                                         # m = 4, V3 degenerate (57 twice)
    prv = se.prv_space(lam, lamp, mu)
    assert prv.m == 4
    res, mats, psis = so.direct_sector(lam, lamp, mu, so.DIRECT_OPS, prv=prv, full=True)
    for n, (sb, fb) in res.items():
        assert sb == 0 and fb == 0, (n, sb, fb)
    for a, b in itertools.combinations(so.DIRECT_OPS, 2):
        assert so.mat_is_zero(so.mat_comm(mats[a], mats[b])), (a, b)
    c3, ok3, c4, ok4 = so.relations(mats, 4)
    assert ok3 and ok4 and c3 == 114 and c4 == 144 == so.cprime_closed(lam, lamp, mu) and c3 == so.c_closed(lam, lamp, mu)
    return '[3,2,1,0]^2, mu = (5,4,2,1), m = 4: slice AND full-vector direct route exact for V3, V4a, V4b, S31, S41; the 5 operators commute; c = %s, c\' = %s equal the closed forms' % (c3, c4)


def test_m_mu_preserved_and_g2_spectrum():
    lam = lamp = (2, 1, 0, 0)
    mu = (3, 2, 1, 0)
    prv, mats = so.triple_mats(lam, lamp, mu)                 # PRV.coords inside raises if the image leaves M_mu
    ref = so.g2_sector(lam, lamp, mu)
    for n in so.TRIPLE:
        assert so.charpoly_asc(mats[n]) == ref['cp'][n]
    cps = so.copies(lam, lamp, mu)
    assert all(abs(a - b) < mp.mpf(10) ** -30 for c1, t2 in zip(cps, ref['trips']) for a, b in zip(c1['ev'], t2))
    return 'sector (3,2,1,0) of [2,1,0,0]^2: reduced matrices preserve M_mu; exact characteristic polynomials and 40-digit joint eigenvalues equal G2\'s full-space result'


def test_copies_one_sector():
    lam = lamp = (3, 2, 1, 0)
    mu = (5, 4, 2, 1)
    r = so.copies_report(lam, lamp, mu, dps=50, g2=True, exch=True)
    cps = r['cps']
    assert r['info']['ok'] and r['info']['offmax'] < mp.mpf(10) ** -40 and r['info']['imag'] == 0
    assert r['selfadj_bad'] == 0 and r['orth_exact'] < mp.mpf(10) ** -40 and r['orth_float'] < 1e-12 and r['hwres'] < 1e-12
    assert r['g2_dtrip'] < mp.mpf(10) ** -30 and r['g2_dslice'] < mp.mpf(10) ** -30
    ev = [tuple(float(x) for x in c['ev']) for c in cps]
    assert abs(ev[0][0] - 46.6076951546) < 1e-9 and abs(ev[1][0] - 57) < 1e-12 and abs(ev[2][0] - 57) < 1e-12 and abs(ev[3][0] - 67.3923048454) < 1e-9     # G2's table
    assert abs(ev[1][1] - 21.3030615433) < 1e-9 and abs(ev[2][1] - 50.6969384567) < 1e-9                                                                      # V3 degenerate: ordered by V4a
    ex = r['exch']
    assert ex['T_ok'] and ex['nfixed'] == 0 and ex['sigma'] == [3, 2, 1, 0] and ex['overlap_dev'] < 1e-12
    return 'copies of (5,4,2,1) in [3,2,1,0]^2: joint-diag off-diagonal %s, self-adjoint, G-orthogonal (%s exact Gram / %.0e float), triples and slices equal G2\'s, exchange pairs (1,4), (2,3)' % (
        mp.nstr(r['info']['offmax'], 2), mp.nstr(r['orth_exact'], 2), r['orth_float'])


def test_closed_forms_and_structure():
    lam = lamp = (3, 2, 1, 0)
    n, bad = so.closed_scalar_check(lam, lamp, verbose=False)
    assert bad == 0
    st = so.structure_checks((2, 1, 1, 0), (2, 1, 1, 0), verbose=False)
    assert st['comm_bad'] == 0 and st['inter_bad'] == 0 and st['simple'] == st['spaces'] and st['contained'] == st['spaces']
    nt, nb = so.generic_lamp_commutation((2, 1, 1, 0), [(Fr(1, 3), Fr(5, 7), Fr(-2), Fr(4, 5))], verbose=False)
    assert nb == 0 and nt > 0
    tot = [so.rect_closed_check(S) for S in (1, 2, 3)]
    assert all(b == 0 for a, b in tot)
    return 'c, c\' closed forms (%d irreps); commutation on whole weight spaces, BGG intertwining E_i^k V(lam\') = V(s_i.lam\') E_i^k, generic non-dominant lam\'; [S,S,0,0]^2 scalar closed forms' % n


def test_reference_file_roundtrip():
    import tempfile
    path = os.path.join(tempfile.gettempdir(), 'g3_test_copies_%d.txt' % os.getpid())
    so.write_copies_ref((2, 1, 0, 0), (2, 1, 0, 0), path=path, verbose=False)
    tot, bad, dev = so.verify_copies_ref(path, full=True, verbose=False)
    os.remove(path)
    assert tot == 0, dict(bad)
    return 'write -> parse -> independent verification of gl4-cg-copies-2100 (VRED vs full operators, char. polys, copies, HW residuals, full-operator eigen-equations): 0 mismatches'


TESTS = [test_normal_order_basic, test_normal_order_random, test_n3_convention, test_hand_rule_matches_engine, test_gl4_direct_route_one_sector, test_m_mu_preserved_and_g2_spectrum,
         test_copies_one_sector, test_closed_forms_and_structure, test_reference_file_roundtrip]


def main():
    nfail = 0
    for f in TESTS:
        t0 = time.perf_counter()
        try:
            msg = f()
            print('  PASS %-36s (%5.2fs)  %s' % (f.__name__, time.perf_counter() - t0, msg), flush=True)
        except Exception as e:
            nfail += 1
            import traceback
            print('  FAIL %-36s (%5.2fs)  %r' % (f.__name__, time.perf_counter() - t0, e), flush=True)
            traceback.print_exc()
    print('%d tests, %d failed' % (len(TESTS), nfail))
    return nfail


if __name__ == '__main__':
    sys.exit(1 if main() else 0)
