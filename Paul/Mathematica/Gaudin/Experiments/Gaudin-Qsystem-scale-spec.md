# Spec Q0: the untwisted Q-system of two-site sl(3) Gaudin at scale — every irrep of [2a,a,0]x[2a,a,0] at a = 7, 10, 14 (round 6, WP-Q0; 2026-10-21)

Agent Q0 (SONNET, general-purpose) — the ONLY agent on the wolfbook kernel K1 this round. Planned by the main session in
`../Gaudin-sl3-round6-plan.md` (§0 target, §1 ingredient (a), §3 row Q0, §5 risk Q0, §6 conventions). This file is the self-contained
brief; the plan is the authority if they disagree. No git commits. Existing `.wl` files and notebooks are never edited (variants go into
`Gaudin-Qsystem-scale.wl`).

## 0. Why (one paragraph)

The Q-route of the project (Q-functions -> Laurent data -> jet formula -> the slice of the HW vector) starts from the ε = 0 Q-system of
every branch B = (mu, copy): the polynomial Q-functions of the untwisted third-order Baxter equation, the pencil spectrum V (the vertex
operator = leading-order Hamiltonian, one root per multiplicity copy) and the untwisted Laurent data (H_s, a_s, b_s) of the oper.
Round 5/today established that `computeV` + `pencilSolve` give exactly `multiplicity` roots for every irrep tested at a <= 7 in 0.05 s
(a = 7, nu = lam) and 0.2 s (a = 14, nu = lam). Round 6 needs this for EVERY irrep of [14,7,0]^2 (and [6,3,0]^2, [8,4,0]^2 as test
cases), with 30-digit data written to reference files for agent S (order-0 points of the eps-series) and agent R (spectrum
cross-check), plus a complexity measurement: how does the cost scale with a and with the sector mu — this is the "poly(|lambda|),
independent of n_W" claim of the Q-route (plan §0).

## 1. Inputs (read their comments; load IN THIS ORDER with `Get`)

1. `Paul/Mathematica/Gaudin/Experiments/gl3-GT-tools.wl` — the notebook's GT machinery (`ps`, `ee`, `V`, `findweight`, `computeV`).
   `computeV[lam,mu][nu]` = the untwisted Q-system of the two-site Gaudin model, sites at u = 1 (rep lam) and u = 0 (rep mu), theta12 = 1;
   returns one `{V, q1(u), q2(u), q3(u)}` per solution; calls `pencilSolve` from (2).
2. `Minimal_complete_working-fast.wl` (repo root) — `pencilSolve`/`charPoly` (the three Baxter equations are linear in the Q coefficients
   c[k,n] and affine in V; the spectrum = roots of the GCD of two exact consistency polynomials), `fastEig[lam,mu][nu, MachinePrecision]`
   (sparse HW null space + V eigenvalues, cross-check at a <= 7), `multRS[lam,mu,wt]` (Racah–Speiser multiplicity).

Conventions fixed by those files (do not change): sites u = 1 (lam) and u = 0 (mu), theta12 = 1; GT order `Reverse@Sort`, index 1 = HW;
V = (1/2) Sum_{s,t in S3} sgn s sgn t e1_{s1 t1} e1_{s2 t2} e2_{s3 t3}; weights: `findweight[lam,mu][nu]` = nu - |nu|/3 + |lam+mu|/3
(central shift, so calling with nu = mu_sector - mu_sector[[3]]*{1,1,1} gives the sector weight itself); Baxter equation
`Q''' + h1 Q'' + (h2 + h1') Q' + (h3 + h2'/2 + h1''/3) Q = 0`, Q_k = l(u) q_k(u), l = 1/((u-1)^lam1 u^mu1), q_k monic of degree
M = (lam1+mu1-E1, lam1+mu1-E2+1, lam1+mu1-E3+2) in echelon form (E = findweight).
W4's closed forms for the untwisted Laurent data of sector mu (header and BRANCH records of `ref/untwisted-exact-310.txt`):
H1 = -H2 = eta/theta12, eta = |Lambda|^2 - (C2(mu) - 2 C2(Lambda))/2, C2(l) = Sum l_i (l_i + 4 - 2i); a1 = -a2 = alpha/theta12,
alpha = (mu1+2)(mu2+1)mu3/2 on multiplicity-1 sectors; in general a1 + a2 + theta1 b1 + theta2 b2 = 0, b1 + b2 = 0.
Laurent convention: S2 = Sum_s [H_s/(u-theta_s) + c2/(u-theta_s)^2], S3 = Sum_s [b_s/(u-theta_s) + a_s/(u-theta_s)^2], c2 = lam2 (lam1+1).

## 2. Tasks

1. **Enumerate and solve.** Enumerate the irreps of [2a,a,0] x [2a,a,0] (Racah–Speiser multiplicities; assert Sum m dim(mu) = dim^2). Run
   `computeV` for EVERY irrep at a = 7, then a = 10, then a = 14. If a = 14 would take more than ~30 min in total, do the 20 slowest
   irreps and extrapolate (say so). Per irrep record: mu, multiplicity, #roots, the roots to 30 digits, the time. Assert #roots = multiplicity
   everywhere; any mismatch is a headline finding (repeated roots? the random projection of `charPoly`? precision?).
2. **Cross-checks.** (i) a <= 7: the roots equal the eigenvalues of `fastEig[lam,lam][nu, MachinePrecision][[1]]` for all irreps (max
   deviation). (ii) Derive the untwisted Laurent data (H_s, a_s, b_s) at both sites from the solved Q's/Baxter coefficients (state the map
   from computeV's h[2], h[3], the E's and V explicitly) and compare H with the eta closed form on all sectors and the pencil root V with
   W4's alpha on >= 20 multiplicity-1 sectors (state the relation as a formula). (iii) State the degrees of q1, q2, q3 as a rule in lam, mu.
3. **Timing.** Per-irrep time vs M1+M2+M3 and vs mu; totals at a = 7, 10, 14; instrument `pencilSolve`'s steps (Det interpolation, GCD,
   per-root LeastSquares) for the five slowest irreps; if one step dominates, try an obvious improvement (e.g. exact sparse
   LinearSolve/LU for the Dets) and report the gain. Do not redesign the algorithm.
4. **Reference files** `Paul/Mathematica/Gaudin/ref/qsystem-spectrum-{630,840,1470}.txt` for [6,3,0], [8,4,0], [14,7,0]: TAB-separated,
   `#` header documenting every record and convention (theta12 = 1, sites at u = 1 (site 1, rep lam), 0 (site 2), V's normalisation,
   weights). Records: `SECTOR mu mult`; `BRANCH mu k V H1 a1 b1 H2 a2 b2` (30 digits, k by increasing V); `QPOLY mu k {coeffs q1} {coeffs q2}
   {coeffs q3}` (ascending powers, 30 digits). Consumers: agent S (order-0 points of the eps-series), agent R (spectrum cross-check).

## 3. Deliverables and acceptance (plan §3 row Q0)

Files (all in `Paul/Mathematica/Gaudin/Experiments/` unless stated): this spec; `Gaudin-Qsystem-scale.wl` (loadable with `Get` after the two
inputs; entry point `runAll[a]` reproducing the a-table); `Gaudin-Qsystem-scale.md` (report: files, conventions, headline, tables,
surprises, not done); `../ref/qsystem-spectrum-{630,840,1470}.txt`; plus (added during the run, Task 1's "per irrep record" for a = 10, 14 and the
a = 7 timings) `../ref/qsystem-roots-a7-a10-a14.txt`. Result in one line: the shipped `pencilSolve` is wrong on C(a,2) + a + 1 irreps of [2a,a,0]^2 (see the report);
all numbers below refer to the corrected solver `qsSolve` of `Gaudin-Qsystem-scale.wl`.
Acceptance: #roots = multiplicity for all 113 irreps at a = 7 and all irreps at a = 10, 14; roots = V-spectrum of `fastEig` (a <= 7) and
H0 = W4's eta/theta12 closed form; timing vs a and vs mu. Risk (plan §5): near-singlet irreps at a = 14 may have large Q-degrees; the exact GCD of two
degree-~40 polynomials over Q may be slow — fall back to numeric roots with a tolerance and report.

## 4. Kernel and style rules

Use K1 via `mcp__wolfbook__wolfbook_evaluateExpression` (pass a large `timeoutSeconds` for runs > 30 s; if it says still running call
`wolfbook_waitEvaluation`); separate statements with `;`; absolute Windows paths with doubled backslashes inside Mathematica strings;
memoised definitions cache bad results (re-run downstream after fixing). Report failures with numbers, not adjectives.
