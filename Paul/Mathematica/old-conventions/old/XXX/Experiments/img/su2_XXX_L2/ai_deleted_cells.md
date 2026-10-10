## Deleted Cell 38 — 2026-07-14 18:22:20
**Kind:** code
```wolfram
(*confirm which (M1,M2) states are "clean" (exact naive product formula) vs genuinely mixed, for (λ1,λ2)=(1,1) and (1,2)*)
cleanCheck11 = Table[{M1,M2,Simplify[buildPsi2Raw[1,1][{M1,M2}][[1]]-tEigenvalue2[1,1,M1,M2][u]]===0},{M1,0,1},{M2,0,1}];
cleanCheck12 = Table[{M1,M2,Simplify[buildPsi2Raw[1,2][{M1,M2}][[1]]-tEigenvalue2[1,2,M1,M2][u]]===0},{M1,0,1},{M2,0,2}];
{cleanCheck11, cleanCheck12}
```

## Deleted Cell 30 — 2026-07-14 18:22:20
**Kind:** code
```wolfram
(*direct diagonalization confirms this for λ1=λ2=1*)
MT=Transpose[B2[1,1][u]];
{evals,evecs}=Simplify[{Eigenvalues[MT],Eigenvectors[MT]}]
```

## Deleted Cell 29 — 2026-07-14 18:22:20
**Kind:** markdown
```markdown
(*RESULT of the naive-tensor-product test: diagonalizing B2[1,1][u] directly (Eigensystem of Transpose[B2[u]], symbolic in u) confirms the eigenvalue guess (u-θ1-s1 h)(u-θ2-s2 h) is exactly right for every state. But the naive tensor covector x2[{s1,s2}]=xSingle[λ1,θ1][{s1}]⊗xSingle[λ2,θ2][{s2}] is only an exact eigenvector at the boundary of the (s1,s2) grid (s2=0 or s1=λ1); interior states pick up genuine mixing with x2[{s1+1,s2-1}] under B2(u).

CORRECTION (superseding the ad hoc fix originally in this section): rather than hand-solving for the admixture coefficients state by state, the correct closed-form SoV covector is obtained exactly as in the L=1 case (x[λ][{s}]=x[λ][{0}].t[λ][1,s][θ]), but acting with the TWO-SITE fused transfer matrices t2[λ1,λ2][1,s][u] (built in the Hirota hierarchy above) on the two-site vacuum, fused to level s1 at θ1 and level s2 at θ2. This gives an EXACT eigenvector directly, with no ad hoc correction term needed — see xCorrected below.*)
```

## Deleted Cell 28 — 2026-07-14 18:22:20
**Kind:** code
```wolfram
(*direct eigenvector check for λ1=λ2=1: all 4 candidate states {s1,s2}∈{0,1}x{0,1}*)
case11=Table[{s1,s2,Simplify[x2[1,1][{s1,s2}].B2[1,1][u]-(u-θ1-s1 h)(u-θ2-s2 h) x2[1,1][{s1,s2}]]//Flatten//Union},{s1,0,1},{s2,0,1}];
case11
```

## Deleted Cell 27 — 2026-07-14 18:22:20
**Kind:** code
```wolfram
ClearAll[x2];
x2[λ1_,λ2_][{s1_,s2_}]:=KroneckerProduct[xSingle[λ1,θ1][{s1}],xSingle[λ2,θ2][{s2}]];
```

## Deleted Cell 26 — 2026-07-14 18:22:20
**Kind:** markdown
```markdown
(*Step 9 — THE KEY CHECK: candidate SoV covector basis x2[λ1,λ2][{s1,s2}] built as the naive tensor product of the single-site fused-transfer-matrix covectors, x2[λ1,λ2][{s1,s2}] = xSingle[λ1,θ1][{s1}] ⊗ xSingle[λ2,θ2][{s2}]. Test directly whether x2[{s1,s2}].B2[λ1,λ2][u] is proportional to x2[{s1,s2}], and if so with what eigenvalue — compare against the naive product-structure guess (u-θ1-s1 h)(u-θ2-s2 h).*)
```

## Deleted Cell 36 — 2026-07-14 20:28:22
**Kind:** code
```wolfram
(*VALIDATION 5c: confirm the fixed n-ordering actually holds — Re[J_n] must be NON-DECREASING in n
  within each sector. Spot-checked on sectors with M>=2 (where >1 solution exists to order):
  (λ1,λ2,M) = (1,1,2), (2,1,2), (2,1,3).*)
order5c = Table[
   With[{sols = BaxterSolve[lp[[1]], lp[[2]], lp[[3]]]},
    {lp, N[sols[[All, 1]], 6], OrderedQ[Re[sols[[All, 1]]]]}],
   {lp, {{1, 1, 2}, {2, 1, 2}, {2, 1, 3}}}];
Print["5c: Re[J_n] non-decreasing in n for every checked sector? ", And @@ (Last /@ order5c)];
Print[If[And @@ (Last /@ order5c), "5c PASSED", "5c FAILED"]];
TableForm[order5c, TableHeadings -> {None, {"(λ1,λ2,M)", "J_n  (n=0..M)", "Re non-decreasing?"}}]
```

## Deleted Cell 35 — 2026-07-14 20:28:22
**Kind:** code
```wolfram
(*VALIDATION 5b: genuine generalization to a NEW, previously-untested representation λ1=2, λ2=1
  (asymmetric; λ1+λ2=3, so M=0..3). Same end-to-end wrapper residual check, plus the full
  (M,n,J_n) table for this case.*)
sampleU5b = {0, 1, 3/2, -2, 10};
res5b = Flatten[Table[
    Table[{M, n,
      Max[Abs[Table[
         aFun[2, 1][uu] Q1[2, 1][M, n][uu + h]
         - τ[2, 1][M, n][uu] Q1[2, 1][M, n][uu]
         + Qθ[uu] Q1[2, 1][M, n][uu - h],
         {uu, sampleU5b}]]]},
     {n, 0, M}],
    {M, 0, 3}], 1];
max5b = Max[res5b[[All, 3]]];
Print["5b max residual (via Q1/τ wrappers, λ1=2,λ2=1) = ", max5b];
Print[If[max5b < 10^-18, "5b PASSED", "5b FAILED"]];
JnTable5b = Flatten[Table[
    Table[{M, n, N[BaxterSolve[2, 1, M][[n + 1, 1]], 10]}, {n, 0, M}],
    {M, 0, 3}], 1];
Column[{
  TableForm[res5b, TableHeadings -> {None, {"M", "n", "max|residual|"}}],
  "(M,n,J_n) table for λ1=2, λ2=1:",
  TableForm[JnTable5b, TableHeadings -> {None, {"M", "n", "J_n"}}]}]
```

## Deleted Cell 34 — 2026-07-14 20:28:22
**Kind:** code
```wolfram
(*VALIDATION 5a: end-to-end check of the ORIGINAL (pre-ansatz-substitution) Baxter equation
    a(u) Q(u+h) - τ(u) Q(u) + d(u) Q(u-h) = 0
  assembled ENTIRELY from the wrapper functions Q1 and τ (NOT the raw q_n/J_n as Cell 31 did),
  for λ1=λ2=1, every (M,n), at several sample u. Exact-rational u so only κ1 carries finite
  (30-digit) precision. This exercises the wrappers — memoization, sorting, and the κ1^(u/h)
  prefactor — top to bottom, independently of the internal matrix-building logic.*)
sampleU5a = {0, 1, 3/2, -2, 10};
res5a = Flatten[Table[
    Table[{M, n,
      Max[Abs[Table[
         aFun[1, 1][uu] Q1[1, 1][M, n][uu + h]
         - τ[1, 1][M, n][uu] Q1[1, 1][M, n][uu]
         + Qθ[uu] Q1[1, 1][M, n][uu - h],
         {uu, sampleU5a}]]]},
     {n, 0, M}],
    {M, 0, 2}], 1];
max5a = Max[res5a[[All, 3]]];
Print["5a max residual (via Q1/τ wrappers, λ1=λ2=1) = ", max5a];
Print[If[max5a < 10^-18, "5a PASSED", "5a FAILED"]];
TableForm[res5a, TableHeadings -> {None, {"M", "n", "max|residual|"}}]
```

## Deleted Cell 32 — 2026-07-14 20:28:22
**Kind:** code
```wolfram
(*Test run: full (M,n,J_n,q_n(u)) table for λ1=λ2=1 (so λ1+λ2=2, M=0,1,2). Values rounded to 10
  digits for display; full 30-digit values are available from BaxterSolve[1,1,M].*)
resultTable = Flatten[
   Table[
    Module[{sols = BaxterSolve[1, 1, M]},
     Table[{M, n - 1, N[sols[[n, 1]], 10], N[qpoly[sols[[n, 2]], u], 10]}, {n, 1, M + 1}]],
    {M, 0, 2}],
   1];
TableForm[resultTable, TableHeadings -> {None, {"M", "n", "J_n", "q_n(u)"}}]
```

## Deleted Cell 31 — 2026-07-14 20:28:22
**Kind:** code
```wolfram
(*INDEPENDENT residual check: re-evaluate the FULL Baxter equation from scratch (using τ=τ0+J_n,
  NOT any BaxterMatrix internals) at several u that are not roots of a(u) or d(u).
  a(u) roots: 4/3, 8/7 ; d(u) roots: 1/3, 1/7 ; sample u in {0,1,3/2,-2,10} avoid all of these.
  Exact rationals are used for u so only κ1 carries finite (30-digit) precision.*)
sampleU = {0, 1, 3/2, -2, 10};
resTable = Flatten[
   Table[
    Module[{sols = BaxterSolve[1, 1, M]},
     Table[
      With[{Jn = sols[[n, 1]], cl = sols[[n, 2]]},
       {M, n - 1,
        Max[Abs[Table[
            κ1 aFun[1, 1][uu] qpoly[cl, uu + h]
            - (τ0[1, 1, M][uu] + Jn) qpoly[cl, uu]
            + κ2 Qθ[uu] qpoly[cl, uu - h],
            {uu, sampleU}]]]}],
      {n, 1, M + 1}]],
    {M, 0, 2}],
   1];
maxResidual = Max[resTable[[All, 3]]];
Print["max residual over all (M,n) and sample u = ", maxResidual];
Print[If[maxResidual < 10^-18, "RESIDUAL CHECK PASSED", "RESIDUAL CHECK FAILED"]];
TableForm[resTable, TableHeadings -> {None, {"M", "n", "max|residual|"}}]
```

## Deleted Cell 29 — 2026-07-14 20:28:22
**Kind:** code
```wolfram
(*Explicit display of the high-degree (u^(M+1), u^(M+2)) coefficients for λ1=λ2=1, M=0,1,2 —
  confirms the cancellation the assertion relies on (should all be ~10^-30).*)
Table[
 Module[{clists},
  clists = Table[
     PadRight[CoefficientList[Expand[κ1 aFun[1, 1][u] (u + h)^k + κ2 Qθ[u] (u - h)^k - τ0[1, 1, M][u] u^k], u], M + 3],
     {k, 0, M}];
  {M, "max|deg M+1,M+2 coeff|" -> Max[Abs[Flatten[clists[[All, {M + 2, M + 3}]]]]]}],
 {M, 0, 2}] // Column
```

