from untwist_limit import *
import pickle
def run(phi, tag):
    pts = load_all(phi, tag); ch = T.chain()
    Fpolys, chat = gf_polys(ch, phi)
    Fd = [poly_dict(P) if P is not None else {} for P in Fpolys]
    out = []
    for bi in range(6):
        vs = branch_psi(pts, bi, ch); PS = comp_series(vs)
        C = [PS[0][nu*8+7] for nu in range(8)]
        xphi = T.x0_D2(ch, phi)
        xn = [T.overlap(xphi, PS[n]) for n in range(8)]
        kB = next(n for n in range(8) if abs(xn[n]) > 1e-25)
        XB = xn[kB]
        nus = nu_series(pts, bi)
        nexp = [[lead_exp(nus[s][i]) for i in range(5)] for s in range(2)]
        NBp = 6 - kB
        R, trip = leading_poly(nus, NBp + 2)
        order = [(t, max(abs(v) for v in R[t].values())) for t in sorted(R)]
        LB = R[NBp]
        GFd = {}
        for nu in range(8):
            for m_, c_ in Fd[nu].items():
                GFd[m_] = GFd.get(m_, mp.mpf(0)) + C[nu]*mpq(c_)
        mons = sorted(set(GFd)|set(LB))
        l = mp.matrix([LB.get(m_, mp.mpf(0)) for m_ in mons]); f = mp.matrix([GFd.get(m_, mp.mpf(0)) for m_ in mons])
        rat = (l.T*f)[0]/(f.T*f)[0]; res = mp.norm(l - rat*f)/mp.norm(l)
        # lower orders vanish?
        low = max([v for t, v in order if t < NBp] + [mp.mpf(0)])
        # triples + needed (i, n)
        G = G_pieces()
        used = {0: set(), 1: set()}; triples = {}
        for (i, j), pc in G.items():
            for (a, b_), gc in pc.items():
                d = a + b_
                for n in range(NBp - d + 1):
                    m = NBp - d - n
                    v = abs(gc * nus[0][i][n] * nus[1][j][m])
                    if v > mp.mpf(10)**-30:
                        used[0].add((i, n)); used[1].add((j, m)); triples[(d, n, m)] = triples.get((d, n, m), 0) + v
        out.append(dict(bi=bi, kB=kB, NBp=NBp, nexp=nexp, low=low, ratio=rat, XB=XB, rel=abs(rat*XB-1), res=res, used=used, triples=sorted(triples), maxn=max(n for (i,n) in used[0]|used[1])))
        print(BNAME[bi], 'k_B', kB, 'nu exps', nexp, 'lower-order max', mp.nstr(low,3), '| L\'_B = r*chat*GF: r*X_B-1 =', mp.nstr(abs(rat*XB-1),3), 'resid', mp.nstr(res,3), '| #triples', len(triples), 'max n', out[-1]['maxn'], 'used(i,n) site1', sorted(used[0]), 'site2', sorted(used[1]))
    return out
a = run(PHI, ''); b = run(PHI2, '2')
