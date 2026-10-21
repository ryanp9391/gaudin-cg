# Q0 report: the untwisted Q-system of two-site sl(3) Gaudin for every irrep of [2a,a,0]^2, a = 1..14 (round 6, 2026-10-21)

Agent Q0 (Sonnet), kernel K1, Mathematica 15.0.0 (Windows, 12 threads, shared with the other round-6 agents). Brief: `Gaudin-Qsystem-scale-spec.md`; plan:
`../Gaudin-sl3-round6-plan.md` (row Q0). Notation: Lam = [2a,a,0] on both sites; the benchmark [14,7,0] is **a = 7**, [28,14,0] is a = 14.
Q-system = the three polynomial solutions q_1,q_2,q_3 of the untwisted third-order Baxter equation of `computeV` (sites u = 1 and u = 0, theta12 = 1)
and the pencil root V (one per multiplicity copy).

## 0. Headline

1. **The shipped `pencilSolve`/`computeV` does not return `multiplicity` roots on all irreps.** It is wrong on **C(a,2) + a + 1** irreps of [2a,a,0]^2
   (a = 7: 29 of 113, a = 10: 56 of 221, a = 14: 106 of 421; a = 1..14 all checked), e.g. 29 roots for the multiplicity-1 irrep (56,14,14) at a = 14. Two mechanisms (Section 2):
   (i) the a+1 irreps with mu1 = 4a (q_1 == 1, no unknowns): `CoefficientArrays[eqs, {}]` silently takes V itself as the unknown, P1 == 0, the GCD is P2 alone;
   (ii) the C(a,2) irreps with mu3 >= a+1 and mu2 <= 2a-1: P1 and P2 share roots at which the degree-M3 equation is inconsistent (the GCD of P1, P2 is not enough; P3 is needed).
   The claim "exactly multiplicity roots for every irrep tested at a <= 7 (113 irreps at a = 7)" does not hold for the shipped files as they are now.
2. **Corrected solver `qsSolve` (exact GCD of the consistency polynomials P1,P2,P3): #roots = multiplicity on every irrep for a = 1..14** — 113/113 irreps (344 copies) at a = 7,
   221/221 (891) at a = 10, **421/421 (2255) at a = 14** — roots real and simple (min spacing 27.7 / 38.1 / 51.9), all three Baxter systems satisfied at every root
   (relative residual <= 1e-62 at 60 digits), identical root sets between the exact (P1,P2,P3 always) and the adaptive solver (1e-55 deviation).
3. **Independent confirmation**: V-spectra equal agent R's slice-engine eigenvalues (`ref/cg-slices-{420,630,840,1470}.txt`; Python, exact PRV slice, no Baxter equation) on **all 492 branches of
   [4,2,0],[6,3,0],[8,4,0],[14,7,0]** (19+44+85+344 = 492 branches, the [14,7,0] ones being the 344 copies of the benchmark) to all 30 printed digits (exact CHARPOLY residual 1e-30);
   equal to `fastEig`'s V-spectrum for a = 1..7 (a = 7: 6.5e-12 absolute, 4e-15 relative);
   reproduce all 45 exact branches of W4 (a_s, b_s, H_s for [1,1,0],[2,1,0],[2,2,0],[3,1,0],[3,2,0],[4,1,0]) to 1e-58.
4. **Laurent data are fixed by V alone, in closed form** (polynomial identities in the sector weight, proved with `Simplify`, checked on every branch): H1 = eta exactly (all sectors),
   a1 = V + eta/2, **b1 = 2(alpha - a1), a2 = a1 - 2 alpha, b2 = -b1**, alpha = (mu1+2)(mu2+1)mu3/2. Consequence: for every sector the **mean of a1 over the copies is alpha**
   (sum_k b1 = 0, verified to 1e-54 on all 2255 copies at a = 14), W4's product formula is the multiplicity-average; multiplicity 1: V = alpha - eta/2 (28/40/56 sectors at a = 7/10/14).
5. **Degrees**: deg q_k = M_k = lam1 + mu1 - E_k + (k-1) (E = sector weight), i.e. (4a-mu1, 4a-mu2+1, 4a-mu3+2), M1+M2+M3 = 6a+3; this is the indicial equation at infinity,
   I_inf(n) = (n-M1)(n-M2)(n-M3), an exact identity. The exponents at u = 0 and u = 1 are {0, a+1, 2a+2} for every sector.
6. **Cost**: the whole Q-system of all 421 irreps of [28,14,0]^2 takes **38 s CPU (48-49 s wall)** with the corrected default solver (0.12 s per irrep on average, max 0.28 s), versus 200 s CPU (216-223 s wall) for the
   shipped code, which also returns wrong roots; a = 7: 3.0 s CPU for all 113 irreps. Per-irrep time spans 0.05-0.28 s while n_W spans 1 .. 25,425: no n_W dependence.
   Scaling in a: time of the full set ~ a^3.7 (corrected, CPU) / a^4.0 (wall, a = 6..14), a^5.2 (shipped, wall; dominated by least-squares solves for spurious roots).
7. Reference files: `ref/qsystem-spectrum-{630,840,1470}.txt` (per sector: V, H1,a1,b1,H2,a2,b2, q-polynomials, 30 digits; validated independently in Python/mpmath: all q_k solve the
   adjoint oper built from the file's own Laurent data to <= 5e-25), and `ref/qsystem-roots-a7-a10-a14.txt` (roots, counts, times of every irrep at a = 7, 10, 14).

## 1. Files and reproduction

| file | content |
|---|---|
| `Experiments/Gaudin-Qsystem-scale-spec.md` | the brief (written first) |
| `Experiments/Gaudin-Qsystem-scale.wl` | the script (prefix `qs`; never edits the two inputs) |
| `Experiments/Gaudin-Qsystem-scale.md` | this report |
| `ref/qsystem-spectrum-630.txt`, `-840.txt`, `-1470.txt` | the reference files ([6,3,0], [8,4,0], [14,7,0]); 30 KB, 65 KB, 390 KB |
| `ref/qsystem-roots-a7-a10-a14.txt` | per-irrep table for a = 7, 10, 14 (755 irreps; roots to 30 digits, counts, times, M, n_W); 150 KB |

Reproduce (fresh kernel, tested: 17 s wall):
```
Get[".../Paul/Mathematica/Gaudin/Experiments/gl3-GT-tools.wl"]; Get[".../Minimal_complete_working-fast.wl"];
Get[".../Paul/Mathematica/Gaudin/Experiments/Gaudin-Qsystem-scale.wl"]; runAll[7]
```
`runAll[a]` runs, for every irrep of [2a,a,0]^2: the shipped `computeV` (count, time), the corrected `qsSolve`, the Laurent checks, `fastEig` (a <= 7), and returns the table
(`"table"`, `"summary"`, `"totals"`); options `"Shipped" -> False`, `"FastEig" -> False`, `"SolverOptions" -> {...}`. `qsWriteRef[Lam, file, 100]` writes a reference file
(it refuses to write if a check fails), `qsWriteTable[{runAll[7], runAll[10], runAll[14]}, file]` the roots table, `qsCompareW4[Lam]` the W4 comparison, `qsShippedTimed[lam,lam,nu]` the step timing of the shipped pipeline, `qsMarkdownTable[res]` the table of Section 10, `qsPencilSolvePatched` is the minimal patch of the
shipped solver (Section 2.4). Direct use for round 7: `sol = qsSolve[lam, lam, nu]` with `nu = mu - mu[[3]] {1,1,1}`; `sol["sols"][[k]]` has `"V"`, `"q"` ({q1,q2,q3} coefficient lists, ascending powers,
monic echelon) and `"resid"`; `qsLaurent[lam, lam, nu, V]` gives `{H1, a1, b1, H2, a2, b2}`.

## 2. The shipped solver is wrong: details

### 2.1 Counts (shipped `computeV` as is; irreps with #roots != multiplicity)

| a | Lam | #irreps | #copies | wrong irreps | = a+1 (M1 = 0) + C(a,2) (band) |
|---|---|---|---|---|---|
| 1 | [2,1,0] | 5 | 6 | 2 | 2 + 0 |
| 2 | [4,2,0] | 13 | 19 | 4 | 3 + 1 |
| 3 | [6,3,0] | 25 | 44 | 7 | 4 + 3 |
| 4 | [8,4,0] | 41 | 85 | 11 | 5 + 6 |
| 5 | [10,5,0] | 61 | 146 | 16 | 6 + 10 |
| 6 | [12,6,0] | 85 | 231 | 22 | 7 + 15 |
| 7 | [14,7,0] | 113 | 344 | 29 | 8 + 21 |
| 8 | [16,8,0] | 145 | 489 | 37 | 9 + 28 |
| 9 | [18,9,0] | 181 | 670 | 46 | 10 + 36 |
| 10 | [20,10,0] | 221 | 891 | 56 | 11 + 45 |
| 11 | [22,11,0] | 265 | 1156 | 67 | 12 + 55 |
| 12 | [24,12,0] | 313 | 1469 | 79 | 13 + 66 |
| 13 | [26,13,0] | 365 | 1834 | 92 | 14 + 78 |
| 14 | [28,14,0] | 421 | 2255 | 106 | 15 + 91 |

(All irreps enumerated by Racah-Speiser; sum m dim = dim^2 checked for every a; the counts agree with `multRS`.) At a = 7 the 29 wrong irreps are, as mu : mult/#roots:
(28,14,0):1/8 (28,13,1):1/9 (28,12,2):1/10 (28,11,3):1/11 (28,10,4):1/12 (28,9,5):1/13 (28,8,6):1/14 (28,7,7):1/15 [M1 = 0] and
(26,8,8):1/3 (25,9,8):2/4 (24,10,8):3/5 (24,9,9):1/5 (23,11,8):4/6 (23,10,9):2/6 (22,12,8):5/7 (22,11,9):3/7 (22,10,10):1/7 (21,13,8):6/8 (21,12,9):4/8 (21,11,10):2/8 (20,13,9):5/7
(20,12,10):3/7 (20,11,11):1/7 (19,13,10):4/6 (19,12,11):2/6 (18,13,11):3/5 (18,12,12):1/5 (17,13,12):2/4 (16,13,13):1/3 [band].

### 2.2 Mechanism (i): M1 = 0 (the a+1 irreps with mu1 = 4a)

q_1 = u^0 = 1 has no unknowns; the first Baxter system is the single equation `V + const = 0`. `pencilSolve` calls `CoefficientArrays[eqs, {}]`; with an empty variable list Mathematica
treats V itself as the unknown: for the adjoint top state `CoefficientArrays[{98 + V}, {}]` returns `{{98}, {{1}}}`, `charPoly` of that is identically 0, and
`PolynomialGCD[0, P2] = P2`: all roots of P2 (a+1+mu3 of them, 8..15 at a = 7) are returned, one genuine. (The `Thread::tdlen` message printed for these irreps is the symptom.)

### 2.3 Mechanism (ii): the band mu3 >= a+1, mu2 <= 2a-1 (C(a,2) irreps)

There P1 and P2 have extra common roots: deg GCD(P1,P2) = min(M1+1, 2a+1-M1) > mult (checked for a = 1..14: the set of irreps with deg GCD(P1,P2) > mult is exactly this band, 0,1,3,6,..,91 irreps).
Example mu = (16,13,13) (mult 1, M = (12,16,17)): GCD(P1,P2) has the three roots V = 1383, 1431, 1479; the genuine one is 1431 (`fastEig`: 1431.0000000000002). At 1383 and 1479 q_1, q_2 exist
but the degree-M3 equation is inconsistent (residual of the third system 191.9 and 48.0 versus an exact zero at 1431); since every polynomial solution has degree M1, M2 or M3 (Section 5), the ODE has exactly two
polynomial solutions there, not three (so these V are not eigenvalues on the HW space: `fastEig` and agent R have only the genuine roots). The exact
fix is the third consistency polynomial: GCD(P1,P2,P3) has degree = multiplicity on every irrep (a = 1..14).

### 2.4 Minimal patch of `pencilSolve` and what it costs

Two changes (see `qsPencilSolvePatched`): GCD over all three systems (`sys` instead of `sys[[;; 2]]`), and `{eqs, ConstantArray[0, {Length[eqs], 0}]}` instead of `CoefficientArrays` when there are no unknowns.
With `Block[{pencilSolve = qsPencilSolvePatched}, computeV[lam, mu][nu]]` all irreps are correct for a = 1..7 (0 wrong of 5,13,25,41,61,85,113), but the symbolic third charPoly makes it slow:
total 0.06 / 0.29 / 0.81 / 2.6 / 6.6 / 17.7 / 32.5 s for a = 1..7 (shipped: 0.06, 0.25, 0.52, 1.4, ..., 5.5 s at a = 7; `qsSolve`: 3.6 s). Use `qsSolve`. NOTE for agent S: its Python port of
the "linear-pencil idea" must include P3 (or a residual filter) and the empty-unknown guard, or it inherits the bug; `ref/qsystem-spectrum-*.txt` is the safe order-0 source.

## 3. The corrected solver

`qsPencil`: with Q = l q, l = (u-1)^-lam1 u^-mu1, the Baxter operator times u^3 (u-1)^3 / l is p3 q''' + p2 q'' + p1 q' + (p0 + V u(u-1)) q with polynomial p_j (degrees 6,5,4,3); V enters only as the multiplication
by u(u-1). Each q_k (monic echelon, computeV's gauge) gives an exact integer system (A0 + V A1) c + (b0 + V b1) = 0 assembled directly on the monomials (banded; no `Together`/`CoefficientArrays`; coefficients
identical to computeV's: checked coefficientwise on all irreps of [6,3,0], deviation 2e-59, and the roots/q's of all shipped-correct irreps).
`qsConsPoly`: P_k(V) = interpolated determinant of a random integer square projection of [A|b](V) (exact, n+1 points), GCD over projections/k. `qsSolve` (default): P1, P2, GCD; P3 and the
final GCD only if deg GCD(P1,P2) != Racah-Speiser multiplicity (91 of 421 irreps at a = 14, exactly the band; no retry was ever needed); roots by `NSolve` at 60 digits; the three linear systems at each root by QR +
back substitution; relative residuals of all three recorded. Exactness: the claim #roots = multiplicity is a statement about exact polynomial GCD degrees, it was also checked with P3 always computed.

### 3.1 Map Baxter coefficients -> Laurent data (state of record)

`computeV`'s operator is -(formal adjoint of d^3 - S1 d^2 + S2 d - S3) acting on Q = l q, hence h1 = S1, h2 = S2 + S1', h3 = S3 + S2'/2 + S1''/6, i.e. **S2 = h2 - h1', S3 = h3 - h2'/2 + h1''/3**.
With S2 = sum_s [H_s/(u-th_s) + c2/(u-th_s)^2], S3 = sum_s [b_s/(u-th_s) + a_s/(u-th_s)^2], th = (1,0), c2 = lam2(lam1+1), and
h2 = hh220/(u-1)^2 + hh211/((u-1)u) + hh020/u^2, h3 = hh330/(u-1)^3 + hh321/((u-1)^2 u) + hh312/((u-1)u^2) + hh303/u^3, hh321 = V (E = findweight, sector weight):
**H1 = -H2 = hh211, a1 = hh321 + hh211/2 = V + H1/2, a2 = -hh312 - hh211/2, b1 = -b2 = hh312 - hh321.**
Checked for every irrep at a = 1,2,3,7,10: the double pole of S2 is c2 at both sites, no triple pole in S2 or S3, and S2, S3 equal their principal parts (so S2 = O(u^-2), S3 = O(u^-3) at infinity).
Closed forms (equal sites, lam3 = 0; `Simplify` gives 0 for both differences, symbolic in lam1, lam2, E1, E2): **hh211 = eta = |Lam|^2 - (C2(E) - 2 C2(Lam))/2**, **hh312 + V = 2 alpha - eta**,
alpha = (E1+2)(E2+1)E3/2, C2(l) = sum_i l_i (l_i + 4 - 2i). I.e. hh312 = 2 alpha - eta - V (the explicit polynomial expressions of hh211 and hh312 in E are `qsHH` in the script).
Hence the relations of Section 0.4, and "H = eta" on all sectors is an identity, not a fit.

## 4. Cross-checks

| check | scope | result |
|---|---|---|
| #roots = multiplicity | all irreps, a = 1..14 (3 sweeps: a = 7,10,14 with P3 always, with P3 adaptive) | 113/113, 221/221, 421/421; deg of the final GCD = mult everywhere; roots real, simple |
| residuals of the 3 Baxter systems | every root, 60 digits | <= 1.9e-63 (a = 7), 1.8e-63 (a = 10), 8.0e-64 (a = 14) with the default solver (<= 8.3e-62 with `LeastSquares`) |
| precision | a = 7, wp = 100 versus 160 (V and all q coefficients) | no difference at 100 digits (relative to max coefficient); coefficient dynamic range per branch <= 3.8e7 |
| agent R (slice engine) | [4,2,0],[6,3,0],[8,4,0],[14,7,0]: 13+25+41+113 sectors, 492 branches | max relative |V_Q0 - V_R| = 0 at the 30 printed digits; R's exact CHARPOLY at my V: 1e-30 |
| `fastEig` V-spectrum | all irreps, a = 1..7; max abs deviation | 1.8e-15, 1.4e-14, 1.1e-13, 3.4e-13, 9.7e-13, 2.2e-12, 6.5e-12 (a = 7; |V| <= 1470) |
| W4 exact branches | 45 branches, 6 representations | count mismatches 0; max |dH| = 0 (exact), |da|, |db| < 1e-58 |
| H1 = eta | all irreps | exact, 113/221/421 of 113/221/421 (a = 7/10/14) |
| V = alpha - eta/2 | multiplicity-1 sectors: 28 / 40 / 56 (a = 7/10/14) | < 1e-55 |
| mean of a1 over copies = alpha | all sectors a = 7, 10, 14 | < 1e-54 |
| reference files (Python/mpmath, independent of the Mathematica hh formulas) | 630: 44 branches, 840: 85, 1470: 344 (3 solutions x 5 points each) | structure/closed-form problems 0; worst relative ODE residual 1.3e-28, 2.2e-27, 5.4e-25 (30-digit data) |

## 5. Degrees and local exponents

Rule (equal sites lam = mu = [lam1,lam2,0], here [2a,a,0]): **deg q_k = M_k = lam1 + mu1 - E_k + (k-1)**, E = findweight[lam,mu][nu] = the sector weight (sum E = |lam|+|mu|):
M = (4a - mu1, 4a - mu2 + 1, 4a - mu3 + 2); M1 < M2 < M3 strictly; M1 + M2 + M3 = 6a + 3 (45, 63, 87 for a = 7, 10, 14); M1 = 0 iff mu1 = 4a; the largest single degree is M3 = 4a+2 (mu3 = 0, e.g. (56,28,0): (0,29,58)),
the singlet has M = (2a, 2a+1, 2a+2) = (28,29,30) at a = 14. **Proof**: u^n -> I_inf(n) u^(n+3) + lower, and I_inf(n) = (n-M1)(n-M2)(n-M3) identically (symbolic in lam1, lam2, E1, E2; also checked on all irreps for a = 1,2,3,4,7,10).
The indicial polynomials at u = 0 and u = 1 have the sector-independent roots {0, lam1-lam2+1, lam1+2} = {0, a+1, 2a+2} (a = 2, 3, 7, all irreps), the exponents of the direct oper of the notes.
(The plan's expectation "near-singlet irreps have the largest Q-degrees (M ~ 42)" is not what happens: near the singlet the three degrees are balanced at ~2a.)

## 6. Timing

All times are wall seconds on a shared 12-thread machine with other heavy jobs running (noise +-50%); CPU times of the kernel (single-threaded) were re-measured for a = 7, 10, 14 with each irrep timed back to back.

### 6.1 Totals

No extrapolation was needed: a = 14 was run completely for all 421 irreps (the shipped baseline alone takes ~4 min, the corrected solver under a minute).

| a | shipped (wrong) CPU / wall | corrected default CPU / wall | corrected median / max per irrep (wall) | M1+M2+M3 | max n_W |
|---|---|---|---|---|---|
| 7 | 4.3 / 5.2 s | 3.0 / 3.4 s | 0.031 / 0.057 s | 45 | 2080 |
| 10 | 18.1 / 22.6 s | 11.3 / 14.4 s | 0.063 / 0.142 s | 63 | 7381 |
| 14 | 200 / 216 s | 38.3 / 47.7 s | 0.110 / 0.281 s | 87 | 25,425 |

Wall totals over all a (shipped / corrected default, s): a = 1..14: 0.035/0.030, 0.12/0.093, 0.31/0.23, 0.72/0.51, 1.44/0.98, 2.8/1.9, 5.5/3.6, 9.4/6.1, 15.1/9.7, 23.6/14.9, 41/24.6, 90/38, 167/43.5, 223/48.9.
Log-log fit over a >= 6: corrected t ~ a^4.0 (wall), a^3.7 (CPU, three points); shipped a^5.2. The corrected median cost per irrep grows ~ a^1.8 (0.031 -> 0.110 s from a = 7 to 14), times the number of irreps ~ a^2/2.
M1+M2+M3 = 6a+3 is the same for every irrep of a given a, so the per-irrep variation is not a function of it; the explanatory variables are the multiplicity m (time ~ linear in m: 0.063 s at m = 1,
0.15 s at m = 8, 0.20 s at m = 14 and 15; a = 14) and the band (needs P3).

### 6.2 Versus mu and n_W (a = 14, corrected default)

| class | #irreps | mean / max corrected (s) | mean / max shipped (s) |
|---|---|---|---|
| M1 = 0 (mu1 = 4a) | 15 | 0.056 / 0.087 | 3.9 / 6.4 |
| band (P3 used) | 91 | 0.137 / 0.281 | 0.76 / 1.4 |
| rest | 315 | 0.113 / 0.207 | 0.30 / 1.7 |

The singlet (n_W = 25,425, the largest) takes 0.10 s, the top irrep (56,28,0) (n_W = 1) 0.08 s; the slowest corrected irreps are band irreps with mu3 = 15..17 and m = 6..12, e.g. (43,26,15): 0.28 s. The
log-log correlation of time with n_W is 0.55 only because n_W and m are both large in the middle of the weight diagram (correlation with m: 0.82); n_W ranges over four orders of magnitude, time over a factor 5.
The plan's n_W ~ 1e5 at a = 14 is 25,425 (n_W,max = C((a+1)^2+1, 2)).

### 6.3 Five slowest irreps, step by step

Shipped pipeline (`qsShippedTimed`, verbatim algorithm with timers; seconds):

| a | irrep (mult, #roots) | total | setup+Baxter assembly | CoefficientArrays | Det interpolation (charPoly) | NSolve | LeastSquares at the roots |
|---|---|---|---|---|---|---|---|
| 7 | (28,7,7) (1, 15) | 0.119 | 0.008 | 0.001 | 0.057 | 0.003 | 0.045 |
| 7 | (21,12,9) (4, 8) | 0.068 | 0.011 | 0.002 | 0.032 | 0.002 | 0.018 |
| 10 | (40,10,10) (1, 21) | 0.264 | 0.008 | 0.001 | 0.144 | 0.005 | 0.097 |
| 10 | (28,16,16) (1, 9) | 0.188 | 0.008 | 0.001 | 0.089 | 0.002 | 0.081 |
| 14 | (56,18,10) (1, 25) | 5.78 | 0.009 | 0.002 | 0.133 | 0.006 | **5.62 (97%)** |
| 14 | (56,22,6) (1, 21) | 5.73 | 0.009 | 0.002 | 0.101 | 0.004 | **5.61 (98%)** |
| 14 | (56,20,8) (1, 23) | 5.68 | 0.010 | 0.002 | 0.117 | 0.005 | **5.53 (97%)** |
| 14 | (56,16,12) (1, 27) | 5.37 | 0.009 | 0.002 | 0.152 | 0.012 | **5.18 (96%)** |
| 14 | (56,19,9) (1, 24) | 5.56 | 0.009 | 0.002 | 0.123 | 0.010 | **5.40 (97%)** |

(The five slowest shipped irreps at a = 14 are all M1 = 0 irreps, mu1 = 56; the table shows two rows each for a = 7, 10.) At a = 14 the per-root `LeastSquares` (50 digits, 3 systems of size up to 56 x 52)
for the 21-27 spurious roots dominates; the genuine root needs only the exact polynomials. Gain: the same five irreps take 0.043-0.056 s with `qsSolve` (P1 is linear when M1 = 0, GCD degree 1 = mult, one root): **a factor 100-140**.
Exact GCD of the consistency polynomials is never the problem (< 1 ms; the plan's concern about degree-40 GCDs is moot); the Det interpolation is next (0.1-0.15 s here).

Corrected solver, five slowest per a (default; seconds, build / det / gcd / nsolve / solve):

| a | irrep (mult) | total | build | det | gcd | nsolve | solve |
|---|---|---|---|---|---|---|---|
| 7 | (21,13,8) (6) | 0.057 | 0.006 | 0.025 | 0.0002 | 0.002 | 0.021 |
| 10 | (37,12,11) (2) | 0.142 | 0.008 | 0.109 | 0.0003 | 0.002 | 0.013 |
| 10 | (30,19,11) (9) | 0.137 | 0.007 | 0.072 | 0.0004 | 0.004 | 0.049 |
| 14 | (43,26,15) (12) | 0.281 | 0.008 | 0.142 | 0.0004 | 0.004 | 0.114 |
| 14 | (41,26,17) (10) | 0.256 | 0.008 | 0.131 | 0.0005 | 0.003 | 0.100 |

### 6.4 Improvements of the corrected solver (a = 14, 26-irrep subset: the 12 slowest plus 14 random; configurations run interleaved, median of 3; seconds)

| configuration | total | build | Det | gcd | nsolve | solve |
|---|---|---|---|---|---|---|
| A: P1,P2,P3 always, 2 projections, `LeastSquares` | 45.1 | 0.37 | 10.9 | 0.04 | 0.10 | 32.8 |
| B: A + QR/back-substitution solves | 15.6 | 0.36 | 10.7 | 0.04 | 0.10 | 3.5 |
| C: B + one projection per P_k | 9.8 | 0.36 | 5.4 | 0.015 | 0.10 | 3.5 |
| D: C + P3 only when deg GCD(P1,P2) != mult (default) | 6.6 | 0.35 | 2.3 | 0.008 | 0.10 | 3.5 |

Gains: QR instead of `LeastSquares` 9.3x on the solve step (26x on a 48 x 43 system, equal to 1e-68); one projection 2x on the Det step (the artefact roots of different k do not coincide, so the GCD over k removes them);
adaptive P3 2.3x; total 6.9x. Full sweep at a = 14: baseline config A 188 s, default D 85 s wall in the first run and 49 s in the clean run; the remaining cost is the QR solves (~50%) and the Det interpolation (~40%).
`Det` of the 57 x 57 integer matrices takes 1.5 ms (already integral; modular `Det` would be ~2x faster; not pursued). Nothing was redesigned.

## 7. Reference files

`ref/qsystem-spectrum-<l1><l2>0.txt`: header documents every convention (theta = (1,0), theta12 = 1, V = hh321 of computeV = eigenvalue of the notebook vertex operator on HW kets, E = mu, gauge of the q's, the map to the Laurent data,
the closed forms, number format: 30 significant digits, `*^` exponent only for |x| < 1e-4 or >= 1e10, integers/rationals exact). Records `REP`, `THETA`, `NSECTORS`, `NBRANCHES`, then per sector (descending lexicographic mu)
`SECTOR mu mult`, and per copy k (increasing V) `BRANCH mu k V H1 a1 b1 H2 a2 b2` and `QPOLY mu k {q1} {q2} {q3}` (ascending powers; last entry 1; exact zeros are gauge or valuation).
[6,3,0]: 25 sectors, 44 branches; [8,4,0]: 41, 85; [14,7,0]: 113, 344. The V values of agent R's `cg-slices-*.txt` COPY records coincide with these. Generated at 100 digits, written only after the checks pass.

## 8. Surprises

1. The shipped `pencilSolve` is wrong on a quarter to two fifths of all irreps (Section 2: 40% at a = 1, 26% at a = 7, 25% at a = 14), including 2 of the 5 irreps of the adjoint [2,1,0]^2, the earliest test case.
2. The failing set has a clean structure: M1 = 0 (a+1 irreps) and the band mu3 >= a+1, mu2 <= 2a-1 (C(a,2) irreps) with deg GCD(P1,P2) = min(M1+1, 2a+1-M1) there; I did not find the representation-theoretic reason (two of the
   three polynomial solutions exist at the spurious V).
3. The Laurent data collapse to one number: b1 = 2(alpha - a1), a2 = a1 - 2 alpha, H1 = eta; W4's "alpha = (mu1+2)(mu2+1)mu3/2" (found by fit on multiplicity-1 sectors) is the **mean** of a1 over the copies on every sector,
   and in multiplicity-1 sectors V = K/2 follows from b1 = 0 alone (K polynomial in E).
4. The plan's cost worry (exact GCD of degree-40 polynomials, large near-singlet degrees) did not materialise: the singlet is among the cheapest irreps; the whole [28,14,0]^2 Q-system takes under a minute.
5. Timing on the shared machine fluctuates by up to a factor 2 between identical runs (the Det step of an identical configuration: 78 s and 134 s); the CPU-time and the interleaved ablation numbers are the reliable ones.

## 9. Not done

- `fastEig` cross-check beyond a = 7 (the chain-side null space at a = 8..14 has n_W up to 25,425); the independent check at a = 7 is R's slice engine (all 344 copies) plus `fastEig` (a <= 7).
- No explanation of the band (mechanism ii) beyond the exact pattern; no proof that sum_k b1 = 0 (verified numerically for a <= 14, expected from the site-exchange symmetry of equal sites).
- The Laurent map is verified for lam3 = mu3 = 0 and equal sites only; unequal sites (`qsSolve` accepts them; the closed forms eta/alpha do not).
- The patched shipped `pencilSolve` was verified for a <= 7 only (it is too slow beyond), the corrected `qsSolve` for a <= 14.
- The ref file for a = 14 ([28,14,0], 2255 branches, ~5 MB) was not written (not requested); `qsWriteRef[{28,14,0}, file, 100]` does it in a few minutes.

## 10. Table a = 7 (shipped baseline vs corrected; generated from `runAll[7]` by `qsMarkdownTable[runAll[7]]`; times in ms; **x** = shipped root count wrong)

Columns: mu | mult | M = (M1,M2,M3) | n_W | #roots shipped | #roots corrected | t shipped (ms) | t corrected (ms) | V roots (increasing; 30 digits in `ref/qsystem-spectrum-1470.txt`).

| mu | m | M | n_W | shipped | corrected | t_ship | t_corr | V |
|---|---|---|---|---|---|---|---|---|
| (28,14,0) | 1 | 0,15,30 | 1 | 8 **x** | 1 | 50 | 20 | -98 |
| (28,13,1) | 1 | 0,16,29 | 2 | 9 **x** | 1 | 61 | 20 | 105 |
| (28,12,2) | 1 | 0,17,28 | 3 | 10 **x** | 1 | 68 | 25 | 279 |
| (28,11,3) | 1 | 0,18,27 | 4 | 11 **x** | 1 | 74 | 23 | 424 |
| (28,10,4) | 1 | 0,19,26 | 5 | 12 **x** | 1 | 81 | 26 | 540 |
| (28,9,5) | 1 | 0,20,25 | 6 | 13 **x** | 1 | 106 | 27 | 627 |
| (28,8,6) | 1 | 0,21,24 | 7 | 14 **x** | 1 | 108 | 29 | 685 |
| (28,7,7) | 1 | 0,22,23 | 8 | 15 **x** | 1 | 124 | 32 | 714 |
| (27,15,0) | 1 | 1,14,30 | 2 | 1 | 1 | 28 | 20 | -105 |
| (27,14,1) | 2 | 1,15,29 | 6 | 2 | 2 | 37 | 24 | 64.6112639, 145.388736 |
| (27,13,2) | 2 | 1,16,28 | 10 | 2 | 2 | 39 | 26 | 231, 343 |
| (27,12,3) | 2 | 1,17,27 | 14 | 2 | 2 | 48 | 27 | 374.102691, 507.897309 |
| (27,11,4) | 2 | 1,18,26 | 18 | 2 | 2 | 48 | 28 | 492.060024, 641.939976 |
| (27,10,5) | 2 | 1,19,25 | 22 | 2 | 2 | 56 | 34 | 584.183851, 745.816149 |
| (27,9,6) | 2 | 1,20,24 | 26 | 2 | 2 | 69 | 32 | 650.147186, 819.852814 |
| (27,8,7) | 2 | 1,21,23 | 30 | 2 | 2 | 71 | 35 | 689.780449, 864.219551 |
| (26,16,0) | 1 | 2,13,30 | 3 | 1 | 1 | 25 | 17 | -111 |
| (26,15,1) | 2 | 2,14,29 | 10 | 2 | 2 | 32 | 21 | 49, 161 |
| (26,14,2) | 3 | 2,15,28 | 20 | 3 | 3 | 37 | 26 | 200.085145, 294, 387.914855 |
| (26,13,3) | 3 | 2,16,27 | 30 | 3 | 3 | 41 | 28 | 337.661502, 456, 574.338498 |
| (26,12,4) | 3 | 2,17,26 | 40 | 3 | 3 | 46 | 29 | 454.823644, 591, 727.176356 |
| (26,11,5) | 3 | 2,18,25 | 50 | 3 | 3 | 50 | 31 | 549.601205, 699, 848.398795 |
| (26,10,6) | 3 | 2,19,24 | 60 | 3 | 3 | 57 | 31 | 621.141572, 780, 938.858428 |
| (26,9,7) | 3 | 2,20,23 | 70 | 3 | 3 | 65 | 35 | 669.015152, 834, 998.984848 |
| (26,8,8) | 1 | 2,21,22 | 74 | 3 **x** | 1 | 71 | 50 | 861 |
| (25,17,0) | 1 | 3,12,30 | 4 | 1 | 1 | 24 | 17 | -116 |
| (25,16,1) | 2 | 3,13,29 | 14 | 2 | 2 | 36 | 26 | 38.1026906, 171.897309 |
| (25,15,2) | 3 | 3,14,28 | 30 | 3 | 3 | 34 | 27 | 181.661502, 300, 418.338498 |
| (25,14,3) | 4 | 3,15,27 | 50 | 4 | 4 | 45 | 32 | 312.424889, 426.967458, 511.032542, 625.575111 |
| (25,13,4) | 4 | 3,16,26 | 70 | 4 | 4 | 44 | 31 | 426.942168, 556.451833, 667.548167, 797.057832 |
| (25,12,5) | 4 | 3,17,25 | 90 | 4 | 4 | 55 | 35 | 522.415545, 664.615506, 793.384494, 935.584455 |
| (25,11,6) | 4 | 3,18,24 | 110 | 4 | 4 | 59 | 37 | 597.495657, 749.460525, 890.539475, 1042.50434 |
| (25,10,7) | 4 | 3,19,23 | 130 | 4 | 4 | 74 | 39 | 651.482738, 810.326120, 959.673880, 1118.51726 |
| (25,9,8) | 2 | 3,20,22 | 142 | 4 **x** | 2 | 66 | 51 | 846.928605, 1001.07140 |
| (24,18,0) | 1 | 4,11,30 | 5 | 1 | 1 | 22 | 16 | -120 |
| (24,17,1) | 2 | 4,12,29 | 18 | 2 | 2 | 28 | 20 | 30.0600240, 179.939976 |
| (24,16,2) | 3 | 4,13,28 | 40 | 3 | 3 | 33 | 26 | 168.823644, 305, 441.176356 |
| (24,15,3) | 4 | 4,14,27 | 70 | 4 | 4 | 37 | 28 | 294.942168, 424.451833, 535.548167, 665.057832 |
| (24,14,4) | 5 | 4,15,26 | 105 | 5 | 5 | 42 | 32 | 406.478435, 538.380624, 630, 721.619376, 853.521565 |
| (24,13,5) | 5 | 4,16,25 | 140 | 5 | 5 | 45 | 35 | 501.481279, 642.286832, 755, 867.713168, 1008.51872 |
| (24,12,6) | 5 | 4,17,24 | 175 | 5 | 5 | 51 | 36 | 578.562100, 727.822613, 855, 982.177387, 1131.43790 |
| (24,11,7) | 5 | 4,18,23 | 210 | 5 | 5 | 60 | 37 | 636.894824, 792.838213, 930, 1067.16179, 1223.10518 |
| (24,10,8) | 3 | 4,19,22 | 235 | 5 **x** | 3 | 61 | 48 | 836.486934, 980, 1123.51307 |
| (24,9,9) | 1 | 4,20,21 | 244 | 5 **x** | 1 | 68 | 45 | 1005 |
| (23,19,0) | 1 | 5,10,30 | 6 | 1 | 1 | 23 | 15 | -123 |
| (23,18,1) | 2 | 5,11,29 | 22 | 2 | 2 | 28 | 21 | 24.1838506, 185.816149 |
| (23,17,2) | 3 | 5,12,28 | 50 | 3 | 3 | 32 | 24 | 159.601205, 309, 458.398795 |
| (23,16,3) | 4 | 5,13,27 | 90 | 4 | 4 | 37 | 30 | 282.415545, 424.615506, 553.384494, 695.584455 |
| (23,15,4) | 5 | 5,14,26 | 140 | 5 | 5 | 40 | 33 | 391.481279, 532.286832, 645, 757.713168, 898.518721 |
| (23,14,5) | 6 | 5,15,25 | 196 | 6 | 6 | 49 | 40 | 485.566634, 629.988930, 737.048926, 816.951074, 924.011070, 1068.43337 |
| (23,13,6) | 6 | 5,16,24 | 252 | 6 | 6 | 51 | 37 | 563.629506, 713.899202, 833.410052, 936.589948, 1056.10080, 1206.37049 |
| (23,12,7) | 6 | 5,17,23 | 308 | 6 | 6 | 60 | 40 | 624.932551, 780.784684, 910.294502, 1027.70550, 1157.21532, 1313.06745 |
| (23,11,8) | 4 | 5,18,22 | 352 | 6 **x** | 4 | 62 | 51 | 829.199208, 965.757264, 1092.24274, 1228.80079 |
| (23,10,9) | 2 | 5,19,21 | 376 | 6 **x** | 2 | 73 | 46 | 999.187767, 1130.81223 |
| (22,20,0) | 1 | 6,9,30 | 7 | 1 | 1 | 22 | 15 | -125 |
| (22,19,1) | 2 | 6,10,29 | 26 | 2 | 2 | 24 | 19 | 20.1471863, 189.852814 |
| (22,18,2) | 3 | 6,11,28 | 60 | 3 | 3 | 28 | 23 | 153.141572, 312, 470.858428 |
| (22,17,3) | 4 | 6,12,27 | 110 | 4 | 4 | 34 | 27 | 273.495657, 425.460525, 566.539475, 718.504343 |
| (22,16,4) | 5 | 6,13,26 | 175 | 5 | 5 | 39 | 32 | 380.562100, 529.822613, 657, 784.177387, 933.437900 |
| (22,15,5) | 6 | 6,14,25 | 252 | 6 | 6 | 43 | 37 | 473.629506, 623.899202, 743.410052, 846.589948, 966.100798, 1116.37049 |
| (22,14,6) | 7 | 6,15,24 | 336 | 7 | 7 | 49 | 40 | 552.044766, 705.683208, 825.983937, 910, 994.016063, 1114.31679, 1267.95523 |
| (22,13,7) | 7 | 6,16,23 | 420 | 7 | 7 | 54 | 43 | 615.291186, 773.031545, 900.227018, 1002, 1103.77298, 1230.96845, 1388.70881 |
| (22,12,8) | 5 | 6,17,22 | 490 | 7 **x** | 5 | 60 | 55 | 824.441877, 957.893448, 1071, 1184.10655, 1317.55812 |
| (22,11,9) | 3 | 6,18,21 | 536 | 7 **x** | 3 | 66 | 52 | 996.933352, 1117, 1237.06665 |
| (22,10,10) | 1 | 6,19,20 | 552 | 7 **x** | 1 | 72 | 41 | 1140 |
| (21,21,0) | 1 | 7,8,30 | 8 | 1 | 1 | 22 | 16 | -126 |
| (21,20,1) | 2 | 7,9,29 | 30 | 2 | 2 | 27 | 22 | 17.7804494, 192.219551 |
| (21,19,2) | 3 | 7,10,28 | 70 | 3 | 3 | 29 | 23 | 149.015152, 314, 478.984848 |
| (21,18,3) | 4 | 7,11,27 | 130 | 4 | 4 | 36 | 29 | 267.482738, 426.326120, 575.673880, 734.517262 |
| (21,17,4) | 5 | 7,12,26 | 210 | 5 | 5 | 36 | 31 | 372.894824, 528.838213, 666, 803.161787, 959.105176 |
| (21,16,5) | 6 | 7,13,25 | 308 | 6 | 6 | 49 | 37 | 464.932551, 620.784684, 750.294502, 867.705498, 997.215316, 1153.06745 |
| (21,15,6) | 7 | 7,14,24 | 420 | 7 | 7 | 52 | 47 | 543.291186, 701.031545, 828.227018, 930, 1031.77298, 1158.96845, 1316.70881 |
| (21,14,7) | 8 | 7,15,23 | 540 | 8 | 8 | 54 | 45 | 607.713066, 768.307900, 897.790929, 993.417302, 1064.58270, 1160.20907, 1289.69210, 1450.28693 |
| (21,13,8) | 6 | 7,16,22 | 644 | 8 **x** | 6 | 59 | 57 | 821.539096, 955.216688, 1060.71880, 1151.28120, 1256.78331, 1390.46090 |
| (21,12,9) | 4 | 7,17,21 | 720 | 8 **x** | 4 | 62 | 52 | 997.435797, 1110.27277, 1211.72723, 1324.56420 |
| (21,11,10) | 2 | 7,18,20 | 760 | 8 **x** | 2 | 65 | 43 | 1140.33437, 1247.66563 |
| (20,20,2) | 1 | 8,9,28 | 74 | 1 | 1 | 22 | 16 | 315 |
| (20,19,3) | 2 | 8,10,27 | 142 | 2 | 2 | 31 | 20 | 426.928605, 581.071395 |
| (20,18,4) | 3 | 8,11,26 | 235 | 3 | 3 | 30 | 24 | 528.486934, 672, 815.513066 |
| (20,17,5) | 4 | 8,12,25 | 352 | 4 | 4 | 36 | 26 | 619.199208, 755.757264, 882.242736, 1018.80079 |
| (20,16,6) | 5 | 8,13,24 | 490 | 5 | 5 | 41 | 31 | 698.441877, 831.893448, 945, 1058.10655, 1191.55812 |
| (20,15,7) | 6 | 8,14,23 | 644 | 6 | 6 | 44 | 35 | 765.539096, 899.216688, 1004.71880, 1095.28120, 1200.78331, 1334.46090 |
| (20,14,8) | 7 | 8,15,22 | 791 | 7 | 7 | 55 | 39 | 819.881676, 955.791618, 1061.25764, 1134, 1206.74236, 1312.20838, 1448.11832 |
| (20,13,9) | 5 | 8,16,21 | 906 | 7 **x** | 5 | 59 | 55 | 999.636971, 1110.05269, 1197, 1283.94731, 1394.36303 |
| (20,12,10) | 3 | 8,17,20 | 979 | 7 **x** | 3 | 58 | 41 | 1143.92109, 1239, 1334.07891 |
| (20,11,11) | 1 | 8,18,19 | 1004 | 7 **x** | 1 | 65 | 35 | 1260 |
| (19,19,4) | 1 | 9,10,26 | 244 | 1 | 1 | 26 | 18 | 675 |
| (19,18,5) | 2 | 9,11,25 | 376 | 2 | 2 | 32 | 21 | 759.187767, 890.812233 |
| (19,17,6) | 3 | 9,12,24 | 536 | 3 | 3 | 35 | 23 | 834.933352, 955, 1075.06665 |
| (19,16,7) | 4 | 9,13,23 | 720 | 4 | 4 | 41 | 29 | 901.435797, 1014.27277, 1115.72723, 1228.56420 |
| (19,15,8) | 5 | 9,14,22 | 906 | 5 | 5 | 54 | 31 | 957.636971, 1068.05269, 1155, 1241.94731, 1352.36303 |
| (19,14,9) | 6 | 9,15,21 | 1070 | 6 | 6 | 52 | 39 | 1002.44775, 1114.25268, 1195.41491, 1254.58509, 1335.74732, 1447.55225 |
| (19,13,10) | 4 | 9,16,20 | 1186 | 6 **x** | 4 | 59 | 48 | 1149.48527, 1237.82402, 1312.17598, 1400.51473 |
| (19,12,11) | 2 | 9,17,19 | 1246 | 6 **x** | 2 | 66 | 41 | 1264.16803, 1345.83197 |
| (18,18,6) | 1 | 10,11,24 | 552 | 1 | 1 | 32 | 22 | 960 |
| (18,17,7) | 2 | 10,12,23 | 760 | 2 | 2 | 46 | 26 | 1020.33437, 1127.66563 |
| (18,16,8) | 3 | 10,13,22 | 979 | 3 | 3 | 42 | 31 | 1073.92109, 1169, 1264.07891 |
| (18,15,9) | 4 | 10,14,21 | 1186 | 4 | 4 | 44 | 32 | 1119.48527, 1207.82402, 1282.17598, 1370.51473 |
| (18,14,10) | 5 | 10,15,20 | 1356 | 5 | 5 | 50 | 34 | 1155.34758, 1243.33940, 1302, 1360.66060, 1448.65242 |
| (18,13,11) | 3 | 10,16,19 | 1462 | 5 **x** | 3 | 53 | 42 | 1271.03624, 1340, 1408.96376 |
| (18,12,12) | 1 | 10,17,18 | 1498 | 5 **x** | 1 | 57 | 35 | 1359 |
| (17,17,8) | 1 | 11,12,22 | 1004 | 1 | 1 | 35 | 19 | 1176 |
| (17,16,9) | 2 | 11,13,21 | 1246 | 2 | 2 | 41 | 26 | 1216.16803, 1297.83197 |
| (17,15,10) | 3 | 11,14,20 | 1462 | 3 | 3 | 49 | 28 | 1251.03624, 1320, 1388.96376 |
| (17,14,11) | 4 | 11,15,19 | 1626 | 4 | 4 | 57 | 35 | 1278.50578, 1342.70090, 1387.29910, 1451.49422 |
| (17,13,12) | 2 | 11,16,18 | 1710 | 4 **x** | 2 | 54 | 38 | 1364.50455, 1419.49545 |
| (16,16,10) | 1 | 12,13,20 | 1498 | 1 | 1 | 39 | 21 | 1329 |
| (16,15,11) | 2 | 12,14,19 | 1710 | 2 | 2 | 46 | 26 | 1352.50455, 1407.49545 |
| (16,14,12) | 3 | 12,15,18 | 1855 | 3 | 3 | 48 | 27 | 1371.80995, 1414, 1456.19005 |
| (16,13,13) | 1 | 12,16,17 | 1904 | 3 **x** | 1 | 54 | 31 | 1431 |
| (15,15,12) | 1 | 13,14,18 | 1904 | 1 | 1 | 52 | 24 | 1425 |
| (15,14,13) | 2 | 13,15,17 | 2016 | 2 | 2 | 56 | 30 | 1435.17068, 1462.82932 |
| (14,14,14) | 1 | 14,15,16 | 2080 | 1 | 1 | 50 | 24 | 1470 |
