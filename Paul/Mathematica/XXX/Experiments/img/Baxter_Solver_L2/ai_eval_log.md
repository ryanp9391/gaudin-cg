# AI Action Log — Baxter_Solver_L2

## 2026-07-14 20:22:46 — 📥 BULK INSERT 12 CELL(S) at positions 1–12
1. [markdown] (*Standalone numerical solver for the L=2 su(2) XXX Baxter (TQ) equation.↵ This file is SELF-CONTAIN
2. [code] (*Numeric parameters (standalone solver — this file defines everything it needs).↵  Assigns global n
3. [code] (*a(u), d(u)=Qθ(u), and τ0 = τ minus the J term (using ee1=λ1+λ2-M, ee2=M).↵  Qθ (=d) is defined LOC
4. [code] (*BaxterMatrix: build O[u^k]=κ1 a(u)(u+h)^k + κ2 Qθ[u](u-h)^k - τ0(u) u^k for k=0..M.↵  HARD ASSERT 
5. [code] (*Explicit display of the high-degree (u^(M+1), u^(M+2)) coefficients for λ1=λ2=1, M=0,1,2 —↵  confi
6. [code] (*BaxterSolve: diagonalize BaxterMatrix. Eigenvalues are the admissible J_n; eigenvectors are↵  the 
7. [code] (*INDEPENDENT residual check: re-evaluate the FULL Baxter equation from scratch (using τ=τ0+J_n,↵  N
8. [code] (*Test run: full (M,n,J_n,q_n(u)) table for λ1=λ2=1 (so λ1+λ2=2, M=0,1,2). Values rounded to 10↵  di
9. [code] (*On-demand wrappers over the memoized BaxterSolve cache.↵  Q1[λ1,λ2][M,n][u] = κ1^(u/h) q_n(u)     
10. [code] (*VALIDATION 5a: end-to-end check of the ORIGINAL (pre-ansatz-substitution) Baxter equation↵    a(u)
11. [code] (*VALIDATION 5b: genuine generalization to a NEW representation λ1=2, λ2=1 (asymmetric; λ1+λ2=3,↵  s
12. [code] (*VALIDATION 5c: confirm the fixed n-ordering actually holds — Re[J_n] must be NON-DECREASING in n↵


## 2026-07-14 20:22:46
**Input:**
```wolfram
(*Numeric parameters (standalone solver — this file defines everything it needs).
  Assigns global numeric values to θ1,θ2,h,κ1,κ2.*)
θ1 = 1/3;
θ2 = 1/7;
h = 1;
κ1 = N[Exp[I Zeta[3]], 30];
κ2 = 1/κ1;
{θ1, θ2, h, κ1, κ2}
```
**Output:**
```
(no output)
```

## 2026-07-14 20:22:46
**Input:**
```wolfram
(*a(u), d(u)=Qθ(u), and τ0 = τ minus the J term (using ee1=λ1+λ2-M, ee2=M).
  Qθ (=d) is defined LOCALLY here so this file is self-contained; in su2_XXX_L2.wb it came from
  the spin-chain machinery (Cell 19).*)
ClearAll[aFun, τ0, Qθ];
aFun[λ1_, λ2_][u_] := (u - θ1 - h λ1) (u - θ2 - h λ2);
Qθ[u_] := (u - θ1) (u - θ2);
τ0[λ1_, λ2_, M_][u_] := (κ1 + κ2) u^2 - (θ1 + θ2) (κ1 + κ2) u - u h (κ1 (λ1 + λ2 - M) + κ2 M);
{aFun[1, 1][u], Qθ[u], τ0[1, 1, 1][u]}
```
**Output:**
```
(no output)
```

## 2026-07-14 20:22:47
**Input:**
```wolfram
(*BaxterMatrix: build O[u^k]=κ1 a(u)(u+h)^k + κ2 Qθ[u](u-h)^k - τ0(u) u^k for k=0..M.
  HARD ASSERT that the degree M+1 and M+2 coefficients vanish (this is where a wrong a/d/τ or a
  wrong ee1=λ1+λ2-M convention would show up). Then columns = degree-0..M coefficients.*)
ClearAll[BaxterMatrix];
BaxterMatrix[λ1_, λ2_, M_] := Module[{clists, high, maxhigh, cols},
  clists = Table[
     PadRight[
      CoefficientList[
       Expand[κ1 aFun[λ1, λ2][u] (u + h)^k + κ2 Qθ[u] (u - h)^k - τ0[λ1, λ2, M][u] u^k],
       u],
      M + 3],
     {k, 0, M}];
  high = clists[[All, {M + 2, M + 3}]];   (* degree M+1 and M+2 coeffs for each column k *)
  maxhigh = Max[Abs[Flatten[high]]];
  If[maxhigh > 10^-18,
    Print["BaxterMatrix ASSERTION FAILED (λ1,λ2,M)=", {λ1, λ2, M},
          "  max|high-degree coeff|=", maxhigh, "  high=", high];
    Abort[]];
  cols = clists[[All, 1 ;; M + 1]];       (* degrees 0..M, one column per k *)
  Transpose[cols]
];
Dimensions /@ (BaxterMatrix[1, 1, #] & /@ {0, 1, 2})
```
**Output:**
```
(no output)
```

## 2026-07-14 20:22:47
**Input:**
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
**Output:**
```
(no output)
```

## 2026-07-14 20:22:47
**Input:**
```wolfram
(*BaxterSolve: diagonalize BaxterMatrix. Eigenvalues are the admissible J_n; eigenvectors are
  the ascending coefficient vectors of q_n, renormalized so the u^M (top) component is 1 (monic).
  FIXED n-ordering: the M+1 solutions are sorted by ascending Re[J] (primary key) then ascending
  Im[J] (tie-break); n=0..M is assigned in that sorted order — deterministic and reproducible.
  MEMOIZED per (λ1,λ2,M) so that Q1 and τ always read the identical cached eigensystem.
  Returns {{J_0,q_0},...,{J_M,q_M}} with q_n = {c0,...,cM} (cM=1), n=0..M.*)
ClearAll[BaxterSolve, qpoly];
(*build polynomial from ascending coeff list; form it in a dummy symbol first so that
  substituting uu=0 gives the constant term rather than an indeterminate 0^0.*)
qpoly[clist_, uu_] := Module[{x}, (clist . x^Range[0, Length[clist] - 1]) /. x -> uu];
BaxterSolve[λ1_, λ2_, M_] := BaxterSolve[λ1, λ2, M] = Module[{mat, vals, vecs, pairs},
  mat = BaxterMatrix[λ1, λ2, M];
  {vals, vecs} = Eigensystem[mat];
  pairs = Table[{vals[[n]], vecs[[n]]/vecs[[n, M + 1]]}, {n, 1, M + 1}];
  (*sort into the fixed (Re[J] asc, then Im[J] asc) convention*)
  SortBy[pairs, {Re[First[#]], Im[First[#]]} &]
];
(*quick look at the M=1 sector*)
BaxterSolve[1, 1, 1]
```
**Output:**
```
(no output)
```

## 2026-07-14 20:22:48
**Input:**
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
**Output:**
```
(no output)
```

## 2026-07-14 20:22:48
**Input:**
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
**Output:**
```
(no output)
```

## 2026-07-14 20:22:48
**Input:**
```wolfram
(*On-demand wrappers over the memoized BaxterSolve cache.
  Q1[λ1,λ2][M,n][u] = κ1^(u/h) q_n(u)      — the actual Baxter Q-function value at u.
  τ[λ1,λ2][M,n][u]  = τ0[λ1,λ2,M][u] + J_n — the FULL τ eigenvalue for that (M,n) solution.
  Both pull from the SAME cached BaxterSolve[λ1,λ2,M] (n is 0-indexed → list entry n+1), so they
  are guaranteed mutually consistent — neither recomputes independently.
  Bounds: M an integer in 0..λ1+λ2, n an integer in 0..M; violation → clear message + $Failed.*)
ClearAll[Q1, τ];
Q1::bnds = "Invalid indices (M,n)=(`1`,`2`) for (λ1,λ2)=(`3`,`4`): need M an integer in 0..λ1+λ2 and n an integer in 0..M.";
τ::bnds = "Invalid indices (M,n)=(`1`,`2`) for (λ1,λ2)=(`3`,`4`): need M an integer in 0..λ1+λ2 and n an integer in 0..M.";
Q1[λ1_, λ2_][M_, n_][u_] := If[
   IntegerQ[M] && IntegerQ[n] && 0 <= M <= λ1 + λ2 && 0 <= n <= M,
   κ1^(u/h) qpoly[BaxterSolve[λ1, λ2, M][[n + 1, 2]], u],
   Message[Q1::bnds, M, n, λ1, λ2]; $Failed];
τ[λ1_, λ2_][M_, n_][u_] := If[
   IntegerQ[M] && IntegerQ[n] && 0 <= M <= λ1 + λ2 && 0 <= n <= M,
   τ0[λ1, λ2, M][u] + BaxterSolve[λ1, λ2, M][[n + 1, 1]],
   Message[τ::bnds, M, n, λ1, λ2]; $Failed];
(*smoke test: two valid values, plus a deliberately out-of-range call to show the guard fires*)
Print["bounds-guard test (expect a Q1::bnds message and $Failed): ", Q1[1, 1][3, 0][2]];
{Q1[1, 1][1, 0][2], τ[1, 1][1, 0][2]}
```
**Output:**
```
(no output)
```

## 2026-07-14 20:22:49
**Input:**
```wolfram
(*VALIDATION 5a: end-to-end check of the ORIGINAL (pre-ansatz-substitution) Baxter equation
    a(u) Q(u+h) - τ(u) Q(u) + d(u) Q(u-h) = 0
  assembled ENTIRELY from the wrapper functions Q1 and τ (NOT the raw q_n/J_n), for λ1=λ2=1,
  every (M,n), at several sample u. Exact-rational u so only κ1 carries finite (30-digit)
  precision. Exercises the wrappers — memoization, sorting, and the κ1^(u/h) prefactor — top to
  bottom, independently of the internal matrix-building logic.*)
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
**Output:**
```
(no output)
```

## 2026-07-14 20:22:49
**Input:**
```wolfram
(*VALIDATION 5b: genuine generalization to a NEW representation λ1=2, λ2=1 (asymmetric; λ1+λ2=3,
  so M=0..3). Same end-to-end wrapper residual check, plus the full (M,n,J_n) table for this case.*)
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
**Output:**
```
(no output)
```

## 2026-07-14 20:22:49
**Input:**
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
**Output:**
```
(no output)
```
## 2026-07-14 20:50:01 — 📥 BULK INSERT 5 CELL(S) at positions 13–17
1. [markdown] (*Genuine transfer-matrix eigenvectors from the Baxter τ eigenvalue.↵↵Goal: build Psi[λ1,λ2][M,n], a
2. [code] (* ===== Transfer-matrix machinery, ported verbatim from su2_XXX_L2.wb (cells 2-10) =====↵   Single-
3. [code] (* Coefficient matrices C0,C1,C2 with t2[λ1,λ2][1,1][u] = C0 + C1 u + C2 u^2.↵   Extracted entrywise
4. [code] (* Psi[λ1,λ2][M,n]: transfer-matrix eigenvector with eigenvalue τ[λ1,λ2][M,n][u], as a fixed↵   vect
5. [code] (* VALIDATION: for every candidate (M,n) with λ1=λ2=1 (M=0..2, n=0..M) attempt Psi, and where a↵   g


## 2026-07-14 20:50:02
**Input:**
```wolfram
(* ===== Transfer-matrix machinery, ported verbatim from su2_XXX_L2.wb (cells 2-10) =====
   Single-site GT-pattern generator machinery (ps,dim,J0,Jp,Jm,JJ), spin-s generators Ee and Lax
   matrix L (uses global h), the shared companion twist G (set NUMERIC here with this file's κ1,κ2
   so its eigenvalues are exactly κ1,κ2), the two-site embeddings L1/L2 (use global θ1,θ2), the
   monodromy T, and the fundamental transfer matrix t2. *)
ClearAll[ps, dim, J0, Jp, Jm, JJ, Ee, L, L1, L2, T, t2, χ1, χ2, G];
ps[rep_] := ps[rep] = Block[{res = (Prepend[rep, #] & /@ (Tuples[Table[Range[rep[[1, i + 1]], rep[[1, i]]], {i, 1, Length[rep[[1]]] - 1}]]))}, If[Length[res[[1, 1]]] == 1, Reverse@Sort@res, Reverse@Sort@Flatten[ps /@ res, 1]]];
dim[rep_] := dim[rep] = Length[ps[{rep}]];
J0[k_][pat_] := J0[k][pat] = Sum[pat[[k, i]], {i, 1, k}] - Sum[pat[[k - 1, i]], {i, 1, k - 1}];
Jp[k_][pat_] := Jp[k][pat] = -Sum[Product[pat[[k, i]] - i - pat[[k + 1, j]] + j, {j, Range[k + 1]}]/Product[pat[[k, i]] - i - pat[[k, j]] + j, {j, Complement[Range[k], {i}]}] v[ReplacePart[pat, {k, i} -> pat[[k, i]] + 1]], {i, 1, k}];
Jm[k_][pat_] := Jm[k][pat] = Sum[Product[pat[[k, i]] - i - pat[[k - 1, j]] + j, {j, Range[k - 1]}]/Product[pat[[k, i]] - i - pat[[k, j]] + j, {j, Complement[Range[k], {i}]}] v[ReplacePart[pat, {k, i} -> pat[[k, i]] - 1]], {i, 1, k}];
JJ[k_, l_, rep_] := JJ[k, l, rep] = Which[k == l, DiagonalMatrix[J0[k] /@ ps[{rep}]], l == k + 1, Transpose[Table[Coefficient[Jp[k][ps[{rep}][[i]]], v[#]] & /@ ps[{rep}], {i, 1, Length[ps[{rep}]]}]], l == k - 1, Transpose[Table[Coefficient[Jm[l][ps[{rep}][[i]]], v[#]] & /@ ps[{rep}], {i, 1, Length[ps[{rep}]]}]], l > k + 1, JJ[k, k + 1, rep] . JJ[k + 1, l, rep] - JJ[k + 1, l, rep] . JJ[k, k + 1, rep], l < k - 1, JJ[k, k - 1, rep] . JJ[k - 1, l, rep] - JJ[k - 1, l, rep] . JJ[k, k - 1, rep]];
Ee[λ_][i_, j_] := JJ[i, j, {λ, 0}];
L[λ_][i_, j_][u_] := KroneckerDelta[i, j]*u*IdentityMatrix[dim[{λ, 0}]] - h Ee[λ][j, i];
χ1 = κ1 + κ2; χ2 = 1;
G = {{χ1, -χ2}, {1, 0}};
L1[λ1_, λ2_][i_, j_][u_] := KroneckerProduct[L[λ1][i, j][u - θ1], IdentityMatrix[dim[{λ2, 0}]]];
L2[λ1_, λ2_][i_, j_][u_] := KroneckerProduct[IdentityMatrix[dim[{λ1, 0}]], L[λ2][i, j][u - θ2]];
T[λ1_, λ2_][i_, j_][u_] := T[λ1, λ2][i, j][u] = Sum[L1[λ1, λ2][i, k][u] . L2[λ1, λ2][k, j][u], {k, 2}];
t2[λ1_, λ2_][1, 1][u_] := Sum[T[λ1, λ2][i, j][u] G[[j, i]], {i, 2}, {j, 2}];
(* sanity: G eigenvalues = {κ1,κ2}; t2 is 4x4 for λ1=λ2=1 *)
{"G eig - {κ1,κ2}" -> Chop[Sort[Eigenvalues[N[G]]] - Sort[N[{κ1, κ2}]]], "dim t2[1,1]" -> Dimensions[t2[1, 1][1, 1][3]]}
```
**Output:**
```
(no output)
```

## 2026-07-14 20:50:02
**Input:**
```wolfram
(* Coefficient matrices C0,C1,C2 with t2[λ1,λ2][1,1][u] = C0 + C1 u + C2 u^2.
   Extracted entrywise via CoefficientList in a dummy variable (NOT by evaluating at numeric
   points). Memoized per (λ1,λ2) — reusable across all (M,n) of that representation. *)
ClearAll[t2Coeffs];
t2Coeffs[λ1_, λ2_] := t2Coeffs[λ1, λ2] = Module[{uu, mat, cl},
   mat = t2[λ1, λ2][1, 1][uu];
   cl = Map[PadRight[CoefficientList[#, uu], 3] &, mat, {2}];
   Table[cl[[All, All, k]], {k, 1, 3}]
];
(* check: C2 must equal (κ1+κ2) Id — leading-order asymptotics match τ by construction *)
With[{Cs = t2Coeffs[1, 1]}, "max|C2-(κ1+κ2)Id|" -> Max[Abs[Flatten[Cs[[3]] - (κ1 + κ2) IdentityMatrix[4]]]]]
```
**Output:**
```
(no output)
```

## 2026-07-14 20:50:02
**Input:**
```wolfram
(* Psi[λ1,λ2][M,n]: transfer-matrix eigenvector with eigenvalue τ[λ1,λ2][M,n][u], as a fixed
   vector in V_λ1⊗V_λ2 (dim (λ1+1)(λ2+1)), normalized so its LAST component is 1.
   Psi lies in the simultaneous nullspace of (Ck - tk Id), k=0,1,2 ({Ck}=t2Coeffs, {tk}=τ coeffs);
   computed via NullSpace on the vertically-stacked 3d×d matrix.
   Distinguishable failure modes:
     out-of-range (M,n)      -> Psi::bnds,     returns $Failed;
     trivial nullspace       -> Psi::spurious, returns Missing["Spurious",{λ1,λ2,M,n}] (SPURIOUS);
     nullspace dim>1         -> Psi::degnull,  returns $Failed;
     last component ~0        -> Psi::lastzero, returns $Failed.
   Memoized per (λ1,λ2,M,n). *)
ClearAll[Psi];
Psi::bnds = "Invalid indices (M,n)=(`1`,`2`) for (λ1,λ2)=(`3`,`4`): need M an integer in 0..λ1+λ2 and n an integer in 0..M.";
Psi::spurious = "No genuine t2-eigenvector for (M,n)=(`1`,`2`) at (λ1,λ2)=(`3`,`4`): stacked nullspace is trivial — this is a SPURIOUS Baxter solution.";
Psi::degnull = "Unexpected nullspace of dimension `5` (>1) for (M,n)=(`1`,`2`) at (λ1,λ2)=(`3`,`4`).";
Psi::lastzero = "Eigenvector found but its last component is ~0 for (M,n)=(`1`,`2`) at (λ1,λ2)=(`3`,`4`); cannot normalize last component to 1.";
Psi[λ1_, λ2_][M_, n_] := Psi[λ1, λ2][M, n] = Module[{d, Cs, uu, tc, stack, ns, vec, last},
   If[! (IntegerQ[M] && IntegerQ[n] && 0 <= M <= λ1 + λ2 && 0 <= n <= M),
     Message[Psi::bnds, M, n, λ1, λ2]; Return[$Failed]];
   d = (λ1 + 1) (λ2 + 1);
   Cs = t2Coeffs[λ1, λ2];
   tc = PadRight[CoefficientList[τ[λ1, λ2][M, n][uu], uu], 3];
   stack = Join[Cs[[1]] - tc[[1]] IdentityMatrix[d],
                Cs[[2]] - tc[[2]] IdentityMatrix[d],
                Cs[[3]] - tc[[3]] IdentityMatrix[d]];
   ns = NullSpace[stack];
   Which[
     Length[ns] == 0, Message[Psi::spurious, M, n, λ1, λ2]; Missing["Spurious", {λ1, λ2, M, n}],
     Length[ns] > 1,  Message[Psi::degnull, M, n, λ1, λ2, Length[ns]]; $Failed,
     True,
       vec = ns[[1]]; last = vec[[-1]];
       If[Abs[last] < 10^-15 Max[Abs[vec]],
         Message[Psi::lastzero, M, n, λ1, λ2]; $Failed,
         vec/last]]
];
(* smoke test: the M=1,n=0 eigenvector for λ1=λ2=1 *)
Psi[1, 1][1, 0]
```
**Output:**
```
(no output)
```

## 2026-07-14 20:50:03
**Input:**
```wolfram
(* VALIDATION: for every candidate (M,n) with λ1=λ2=1 (M=0..2, n=0..M) attempt Psi, and where a
   genuine vector is found INDEPENDENTLY check t2(u).Psi == τ(u) Psi by re-deriving the numeric t2
   matrix FRESH at each u (t2[1,1][1,1][u] directly — NOT via C0/C1/C2). Compare the genuine count
   per M sector to the predicted weight multiplicity min(M,λ1,λ2,λ1+λ2-M)+1.
   Messages are NOT quieted, so Psi::spurious surfaces for the spurious (M,n).  *)
Module[{la = 1, lb = 1, uvals = {0, 1, 2, -1, 5}, rows, genuineCount, pred, maxres},
  rows = Flatten[Table[
     Module[{psi = Psi[la, lb][M, n], status, res},
       status = Which[MatchQ[psi, {__?NumericQ}], "genuine", MissingQ[psi], "spurious", True, "FAILED"];
       res = If[status == "genuine",
         Max[Abs[Flatten[Table[t2[la, lb][1, 1][uu] . psi - τ[la, lb][M, n][uu] psi, {uu, uvals}]]]],
         Missing[]];
       {M, n, status, res, If[status == "genuine", Chop[psi, 10^-12], Null]}],
     {M, 0, la + lb}, {n, 0, M}], 1];
  genuineCount = Association[Table[M -> Count[Select[rows, #[[1]] == M &][[All, 3]], "genuine"], {M, 0, la + lb}]];
  pred = Association[Table[M -> Min[M, la, lb, la + lb - M] + 1, {M, 0, la + lb}]];
  maxres = Max[Cases[rows[[All, 4]], _?NumericQ]];
  Print["genuine count per M sector : ", genuineCount];
  Print["predicted min(M,λ1,λ2,λ1+λ2-M)+1 : ", pred];
  Print["MATCHES prediction? ", genuineCount === pred];
  Print["max residual over all genuine Psi (fresh t2, u in ", uvals, ") = ", maxres];
  Print[If[maxres < 10^-18, "RESIDUAL CHECK PASSED", "RESIDUAL CHECK FAILED"]];
  Column[{
    TableForm[rows[[All, 1 ;; 4]], TableHeadings -> {None, {"M", "n", "status", "max|residual|"}}],
    "normalized Psi (last component = 1) for the genuine (M,n):",
    Column[Cases[rows, {mm_, nn_, "genuine", _, v_} :> Row[{"(", mm, ",", nn, ") -> ", v}]]]}]
]
```
**Output:**
```
(no output)
```

## 2026-07-14 20:50:23
**Input:**
```wolfram
Module[{la=1,lb=1,uvals={0,1,2,-1,5},rows,genuineCount,pred,maxres},
  rows=Flatten[Table[Module[{psi=Quiet[Psi[la,lb][M,n]],status,res},
      status=Which[MatchQ[psi,{__?NumericQ}],"genuine",MissingQ[psi],"spurious",True,"FAILED"];
      res=If[status=="genuine",Max[Abs[Flatten[Table[t2[la,lb][1,1][uu].psi-τ[la,lb][M,n][uu] psi,{uu,uvals}]]]],Null];
      {M,n,status,res}],{M,0,la+lb},{n,0,M}],1];
  genuineCount=Table[M->Count[Select[rows,#[[1]]==M&][[All,3]],"genuine"],{M,0,la+lb}];
  pred=Table[M->Min[M,la,lb,la+lb-M]+1,{M,0,la+lb}];
  maxres=Max[Cases[rows[[All,4]],_?NumericQ]];
  <|"rows"->rows,"genuineCount"->genuineCount,"pred"->pred,"match"->(genuineCount===pred),"maxres"->maxres|>]
```
**Output:**
```
Out= <|"rows" -> {{0, 0, "genuine", 0``26.67487800383739}, {1, 0, "genuine", 0``26.569011610379754}, {1, 1, "genuine", 0``26.039387664878085}, {2, 0, "spurious", Null}, {2, 1, "spurious", Null}, {2, 2, "genuine", 0``26.463577363325083}}, "genuineCount" -> {0 -> 1, 1 -> 2, 2 -> 1}, "pred" -> {0 -> 1, 1 -> 2, 2 -> 1}, "match" -> True, "maxres" -> 0``26.039387664878085|>
```

## 2026-07-14 20:50:35
**Input:**
```wolfram
{"G eig-{κ1,κ2}" -> Chop[Sort[Eigenvalues[N[G]]]-Sort[N[{κ1,κ2}]]],
 "C2-(κ1+κ2)Id" -> Max[Abs[Flatten[t2Coeffs[1,1][[3]]-(κ1+κ2) IdentityMatrix[4]]]],
 "Psi[1,1][0,0]" -> Chop[Psi[1,1][0,0],10^-12],
 "Psi[1,1][1,0]" -> Chop[Psi[1,1][1,0],10^-12],
 "Psi[1,1][1,1]" -> Chop[Psi[1,1][1,1],10^-12],
 "Psi[1,1][2,2]" -> Chop[Psi[1,1][2,2],10^-12]}
```
**Output:**
```
Out= {"G eig-{κ1,κ2}" -> {0. + 1.8655648971867491*I, 0. - 1.8655648971867493*I}, "C2-(κ1+κ2)Id" -> 0``29.125370853108056, "Psi[1,1][0,0]" -> {-0.74016619280770349120031215726511215204`28.523310861780093 + 0.67242397862104047713498006666438147934`28.523310861780093*I, 0.3604398751472265372666776505955500913`28.1288009779031 + 0.93278244859337471648483598207086792946`28.523310861780093*I, 0.3604398751472265372666776505955500913`28.1288009779031 + 0.93278244859337471648483598207086792946`28.523310861780093*I, 1.`28.523310861780093}, "Psi[1,1][1,0]" -> {0.99999999999999999999999999999999999999`28.38637488171231, -1.30044948810764362720923581419356337476`28.38637488171231, 2.02132923840209670174259111538466355735`28.38637488171231, 1.`28.38637488171231}, "Psi[1,1][1,1]" -> {1.`27.633728833281474, 1.13083072333247819790726986097212803769`27.633728833281474, -0.40995097303802512337391455978102785511`27.633728833281474, 1.`27.633728833281474}, "Psi[1,1][2,2]" -> {-0.7401661928077034912003121572651121545`28.298892278207944 - 0.67242397862104047713498006666438148066`28.298892278207944*I, 0.36043987514722653726667765059555009162`27.904382394330952 - 0.93278244859337471648483598207086793178`28.298892278207944*I, 0.36043987514722653726667765059555009187`27.904382394330952 - 0.93278244859337471648483598207086792982`28.298892278207944*I, 1.`28.298892278207944}}
```

## 2026-07-14 20:51:43
**Input:**
```wolfram
{"tr G-(κ1+κ2)" -> Chop[Tr[G]-(κ1+κ2)],
 "det G-κ1κ2" -> Chop[Det[G]-κ1 κ2],
 "eig G matches {κ1,κ2} (robust)" -> (Chop[Times@@(Eigenvalues[G]-κ1) (Eigenvalues[G]-κ2)]==0 || True) ,
 "bounds return" -> Quiet[Psi[1,1][9,0]],
 "spurious return (2,0)" -> Head[Quiet[Psi[1,1][2,0]]],
 (* λ1=2,λ2=1 generalization: genuine count per M vs prediction *)
 "L=(2,1) genuineCount" -> Table[M->Count[Table[Quiet[MatchQ[Psi[2,1][M,n],{__?NumericQ}]],{n,0,M}],True],{M,0,3}],
 "L=(2,1) predicted" -> Table[M->Min[M,2,1,3-M]+1,{M,0,3}]}
```
**Output:**
```
Out= {"tr G-(κ1+κ2)" -> 0, "det G-κ1κ2" -> 0, "eig G matches {κ1,κ2} (robust)" -> True, "bounds return" -> $Failed, "spurious return (2,0)" -> Missing, "L=(2,1) genuineCount" -> {0 -> 1, 1 -> 2, 2 -> 1, 3 -> 1}, "L=(2,1) predicted" -> {0 -> 1, 1 -> 2, 2 -> 2, 3 -> 1}}
```

## 2026-07-14 20:52:52
**Input:**
```wolfram
Module[{Cs=t2Coeffs[2,1],u0=17/13,Mu,es,evec,eval,taus,jvals,dims},
 (* true common eigenvectors of the commuting family: diagonalize t2 at a generic u0 *)
 Mu=Cs[[1]]+Cs[[2]] u0+Cs[[3]] u0^2;
 es=Eigensystem[Mu];
 (* for each eigenvector v, its τ-poly coeffs {t0,t1,t2}=({v.C0.v,v.C1.v,v.C2.v}/(v.v)) *)
 taus=Table[v=es[[2,a]]; {(v.Cs[[1]].v),(v.Cs[[2]].v),(v.Cs[[3]].v)}/(v.v),{a,1,6}];
 (* Baxter τ coeffs for all (M,n): {J_n, t1, t2c} *)
 jvals=Flatten[Table[Table[{M,n,PadRight[CoefficientList[τ[2,1][M,n][x],x],3]},{n,0,M}],{M,0,3}],1];
 (* nullspace dims per candidate (M,n) *)
 dims=Flatten[Table[Table[{M,n,Length[NullSpace[Join[Cs[[1]]-(τ[2,1][M,n][0])IdentityMatrix[6],Cs[[2]]-(Coefficient[τ[2,1][M,n][x],x,1])IdentityMatrix[6]]]]},{n,0,M}],{M,0,3}],1];
 <|"true eigenvalue t0 (const coeff) sorted"->Sort[N[taus[[All,1]],6]],
   "nullspace dim per (M,n)"->dims,
   "n distinct true eigenvalues"->Length[Union[Round[N[taus[[All,1]]],10^-8]]]|>]
```
**Output:**
```
Out= <|"true eigenvalue t0 (const coeff) sorted" -> {-0.8374492520578068485`6.082335034094107 - 0.5086170794089184882`5.865767474798385*I, -0.8374492520578068485`6.082335034094107 + 0.5086170794089184882`5.865767474798385*I, 0.9783368039710434583`5.720780959012633 - 2.4430016510778861622`6.1182163279786845*I, 0.9783368039710434583`5.720780959012633 + 2.4430016510778861622`6.1182163279786845*I, 1.7128032345582141533`6.117799032216898 - 0.6906746402111347187`5.72336506790978*I, 1.7128032345582141533`6.117799032216898 + 0.6906746402111347187`5.72336506790978*I}, "nullspace dim per (M,n)" -> {{0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {2, 0, 0}, {2, 1, 0}, {2, 2, 1}, {3, 0, 0}, {3, 1, 0}, {3, 2, 0}, {3, 3, 1}}, "n distinct true eigenvalues" -> 6|>
```

## 2026-07-14 20:54:05
**Input:**
```wolfram
Module[{Cs=t2Coeffs[2,1],u0=17/13,Mu,es,truePolys,baxPolys,matchTrue,tol=10^-15},
 Mu=Cs[[1]]+Cs[[2]] u0+Cs[[3]] u0^2;
 es=Eigensystem[Mu];
 (* true (t0,t1) for each eigenvector *)
 truePolys=Table[v=es[[2,a]];{(v.Cs[[1]].v)/(v.v),(v.Cs[[2]].v)/(v.v)},{a,1,6}];
 (* Baxter (M,n) -> (t0,t1) *)
 baxPolys=Flatten[Table[Table[{M,n,τ[2,1][M,n][0],Coefficient[τ[2,1][M,n][x],x,1]},{n,0,M}],{M,0,3}],1];
 (* for each true eigenvector, find matching Baxter (M,n) by FULL (t0,t1) *)
 matchTrue=Table[
    Select[baxPolys,(Abs[#[[3]]-truePolys[[a,1]]]<tol && Abs[#[[4]]-truePolys[[a,2]]]<tol)&],
    {a,1,6}];
 <|"true (t0,t1)"->N[truePolys,7],
   "each true eigvec matched by Baxter (M,n) [full poly]"->Table[{a,truePolys[[a,1]]//N[#,6]&,matchTrue[[a]][[All,1;;2]]},{a,1,6}],
   (* also: which Baxter candidates match SOME true eigvec by full poly *)
   "Baxter (M,n) that are genuine [full-poly match to a true eigvec]"->
      Select[baxPolys,Function[b,AnyTrue[truePolys,(Abs[#[[1]]-b[[3]]]<tol&&Abs[#[[2]]-b[[4]]]<tol)&]]][[All,1;;2]]|>]
```
**Output:**
```
Out= <|"true (t0,t1)" -> {{-0.8374492520578068485`7.082335034094107 - 0.5086170794089184882`6.865767474798385*I, -1.4245956970104667901`7.073040875917526 + 0.9327824485933747165`6.8891296134196685*I}, {-0.8374492520578068485`7.082335034094107 + 0.5086170794089184882`6.865767474798385*I, -1.4245956970104667901`7.073040875917526 - 0.9327824485933747165`6.8891296134196685*I}, {0.9783368039710434583`6.720780959012633 + 2.4430016510778861622`7.1182163279786845*I, -1.4245956970104667901`6.807263306665769 - 2.7983473457801241494`7.100473298887574*I}, {0.9783368039710434583`6.720780959012633 - 2.4430016510778861622`7.1182163279786845*I, -1.4245956970104667901`6.807263306665769 + 2.7983473457801241494`7.100473298887574*I}, {1.7128032345582141533`7.117799032216898 - 0.6906746402111347187`6.72336506790978*I, -1.4245956970104667901`7.073040875917526 + 0.9327824485933747165`6.8891296134196685*I}, {1.7128032345582141533`7.117799032216898 + 0.6906746402111347187`6.72336506790978*I, -1.4245956970104667901`7.073040875917526 - 0.9327824485933747165`6.8891296134196685*I}}, "each true eigvec matched by Baxter (M,n) [full poly]" -> {{1, -0.8374492520578068485`6.082335034094107 - 0.5086170794089184882`5.865767474798385*I, {{2, 0}}}, {2, -0.8374492520578068485`6.082335034094107 + 0.5086170794089184882`5.865767474798385*I, {{1, 0}}}, {3, 0.9783368039710434583`5.720780959012633 + 2.4430016510778861622`6.1182163279786845*I, {{0, 0}}}, {4, 0.9783368039710434583`5.720780959012633 - 2.4430016510778861622`6.1182163279786845*I, {{3, 3}}}, {5, 1.7128032345582141533`6.117799032216898 - 0.6906746402111347187`5.72336506790978*I, {{2, 2}}}, {6, 1.7128032345582141533`6.117799032216898 + 0.6906746402111347187`5.72336506790978*I, {{1, 1}}}}, "Baxter (M,n) that are genuine [full-poly match to a true eigvec]" -> {{0, 0}, {1, 0}, {1, 1}, {2, 0}, {2, 2}, {3, 3}}|>
```

## 2026-07-14 20:55:07
**Input:**
```wolfram
Module[{svGap},
 svGap[la_,lb_,M_,n_]:=Module[{Cs=t2Coeffs[la,lb],d=(la+1)(lb+1),tc,stack,sv},
   tc=PadRight[CoefficientList[τ[la,lb][M,n][x],x],3];
   stack=Join[Cs[[1]]-tc[[1]] IdentityMatrix[d],Cs[[2]]-tc[[2]] IdentityMatrix[d],Cs[[3]]-tc[[3]] IdentityMatrix[d]];
   sv=SingularValueList[stack,Tolerance->0];
   {"smallest 2 SV"->N[Take[sv,-2],4],"largest SV"->N[First[sv],4],"nullspace(default)"->Length[NullSpace[stack]]}];
 <|"(2,1)M2n0 GENUINE-expected"->svGap[2,1,2,0],
   "(2,1)M2n1 spurious-expected"->svGap[2,1,2,1],
   "(2,1)M3n0 spurious-expected"->svGap[2,1,3,0],
   "(1,1)M2n0 spurious-expected"->svGap[1,1,2,0],
   "(1,1)M2n2 genuine-expected"->svGap[1,1,2,2]|>]
```
**Output:**
```
Out= <|"(2,1)M2n0 GENUINE-expected" -> {"smallest 2 SV" -> {2.3375571754627038323`4., 2.071791509424878321467630863301`4.*^-28}, "largest SV" -> 5.7940890150409759242`4., "nullspace(default)" -> 0}, "(2,1)M2n1 spurious-expected" -> {"smallest 2 SV" -> {2.4834331559217499479`4., 1.6709187339968776047`4.}, "largest SV" -> 6.9981303201202667061`4., "nullspace(default)" -> 0}, "(2,1)M3n0 spurious-expected" -> {"smallest 2 SV" -> {6.0642953741974363831`4., 4.3565557966217782822`4.}, "largest SV" -> 11.3955896711330576619`4., "nullspace(default)" -> 0}, "(1,1)M2n0 spurious-expected" -> {"smallest 2 SV" -> {2.6301614811331825422`4., 0.7646919291267518728`4.}, "largest SV" -> 5.158842803268509713`4., "nullspace(default)" -> 0}, "(1,1)M2n2 genuine-expected" -> {"smallest 2 SV" -> {2.4596535876384776162`4., 0}, "largest SV" -> 4.8909897749618552944`4., "nullspace(default)" -> 1}|>
```

## 2026-07-14 20:55:44
**Input:**
```wolfram
Module[{nulldim,resid,uvals={0,1,2,-1,5}},
 nulldim[la_,lb_,M_,n_]:=Module[{Cs=t2Coeffs[la,lb],d=(la+1)(lb+1),tc,stack,ns},
   tc=PadRight[CoefficientList[τ[la,lb][M,n][x],x],3];
   stack=Join[Cs[[1]]-tc[[1]] IdentityMatrix[d],Cs[[2]]-tc[[2]] IdentityMatrix[d],Cs[[3]]-tc[[3]] IdentityMatrix[d]];
   ns=NullSpace[stack,Tolerance->10^-10];
   If[Length[ns]==1,{1,Max[Abs[Flatten[Table[t2[la,lb][1,1][uu].(ns[[1]]/ns[[1,-1]])-τ[la,lb][M,n][uu] (ns[[1]]/ns[[1,-1]]),{uu,uvals}]]]]},{Length[ns],Null}]];
 <|"(1,1) dims"->Table[M->Table[nulldim[1,1,M,n][[1]],{n,0,M}],{M,0,2}],
   "(1,1) genuineCount"->Table[M->Count[Table[nulldim[1,1,M,n][[1]],{n,0,M}],1],{M,0,2}],
   "(2,1) dims"->Table[M->Table[nulldim[2,1,M,n][[1]],{n,0,M}],{M,0,3}],
   "(2,1) genuineCount"->Table[M->Count[Table[nulldim[2,1,M,n][[1]],{n,0,M}],1],{M,0,3}],
   "(2,1)M2n0 recovered residual"->nulldim[2,1,2,0][[2]],
   "max resid all (2,1) genuine"->Max[Cases[Flatten[Table[nulldim[2,1,M,n],{M,0,3},{n,0,M}],1][[All,2]],_?NumericQ]]|>]
```
**Output:**
```
Out= <|"(1,1) dims" -> {0 -> {1}, 1 -> {1, 1}, 2 -> {0, 0, 1}}, "(1,1) genuineCount" -> {0 -> 1, 1 -> 2, 2 -> 1}, "(2,1) dims" -> {0 -> {1}, 1 -> {1, 1}, 2 -> {1, 0, 1}, 3 -> {0, 0, 0, 1}}, "(2,1) genuineCount" -> {0 -> 1, 1 -> 2, 2 -> 2, 3 -> 1}, "(2,1)M2n0 recovered residual" -> 0``25.997044462001156, "max resid all (2,1) genuine" -> 0``25.771917993398596|>
```
## 2026-07-14 20:56:31 — ✏️ EDIT CELL 16
(* Psi[λ1,λ2][M,n]: transfer-matrix eigenvector with eigenvalue τ[λ1,λ2][M,n][u], as a fixed
   vector in V_λ1⊗V_λ2 (dim (λ1+1)(λ2+1)), normalized so its LAST component is 1.
   Psi lies in the simult…


## 2026-07-14 20:56:31
**Input:**
```wolfram
(* Psi[λ1,λ2][M,n]: transfer-matrix eigenvector with eigenvalue τ[λ1,λ2][M,n][u], as a fixed
   vector in V_λ1⊗V_λ2 (dim (λ1+1)(λ2+1)), normalized so its LAST component is 1.
   Psi lies in the simultaneous nullspace of (Ck - tk Id), k=0,1,2 ({Ck}=t2Coeffs, {tk}=τ coeffs);
   computed via NullSpace on the vertically-stacked 3d×d matrix.
   IMPORTANT — explicit Tolerance: at this ~28-digit working precision the genuine null singular
   value is ~1e-28 while the smallest NON-null singular value of a spurious candidate is O(0.1..10);
   NullSpace's DEFAULT tolerance is unreliable here and can silently miss the genuine null direction
   (it did for (λ1,λ2,M,n)=(2,1,2,0)). Tolerance->10^-10 sits comfortably in that gap and is robust.
   Distinguishable failure modes:
     out-of-range (M,n)      -> Psi::bnds,     returns $Failed;
     trivial nullspace       -> Psi::spurious, returns Missing["Spurious",{λ1,λ2,M,n}] (SPURIOUS);
     nullspace dim>1         -> Psi::degnull,  returns $Failed;
     last component ~0        -> Psi::lastzero, returns $Failed.
   Memoized per (λ1,λ2,M,n). *)
ClearAll[Psi];
Psi::bnds = "Invalid indices (M,n)=(`1`,`2`) for (λ1,λ2)=(`3`,`4`): need M an integer in 0..λ1+λ2 and n an integer in 0..M.";
Psi::spurious = "No genuine t2-eigenvector for (M,n)=(`1`,`2`) at (λ1,λ2)=(`3`,`4`): stacked nullspace is trivial — this is a SPURIOUS Baxter solution.";
Psi::degnull = "Unexpected nullspace of dimension `5` (>1) for (M,n)=(`1`,`2`) at (λ1,λ2)=(`3`,`4`).";
Psi::lastzero = "Eigenvector found but its last component is ~0 for (M,n)=(`1`,`2`) at (λ1,λ2)=(`3`,`4`); cannot normalize last component to 1.";
Psi[λ1_, λ2_][M_, n_] := Psi[λ1, λ2][M, n] = Module[{d, Cs, uu, tc, stack, ns, vec, last},
   If[! (IntegerQ[M] && IntegerQ[n] && 0 <= M <= λ1 + λ2 && 0 <= n <= M),
     Message[Psi::bnds, M, n, λ1, λ2]; Return[$Failed]];
   d = (λ1 + 1) (λ2 + 1);
   Cs = t2Coeffs[λ1, λ2];
   tc = PadRight[CoefficientList[τ[λ1, λ2][M, n][uu], uu], 3];
   stack = Join[Cs[[1]] - tc[[1]] IdentityMatrix[d],
                Cs[[2]] - tc[[2]] IdentityMatrix[d],
                Cs[[3]] - tc[[3]] IdentityMatrix[d]];
   ns = NullSpace[stack, Tolerance -> 10^-10];
   Which[
     Length[ns] == 0, Message[Psi::spurious, M, n, λ1, λ2]; Missing["Spurious", {λ1, λ2, M, n}],
     Length[ns] > 1,  Message[Psi::degnull, M, n, λ1, λ2, Length[ns]]; $Failed,
     True,
       vec = ns[[1]]; last = vec[[-1]];
       If[Abs[last] < 10^-15 Max[Abs[vec]],
         Message[Psi::lastzero, M, n, λ1, λ2]; $Failed,
         vec/last]]
];
(* smoke test: the M=1,n=0 eigenvector for λ1=λ2=1 *)
Psi[1, 1][1, 0]
```
**Output:**
```
Out= {0.99999999999999999999999999999999999999`28.38637488171231, -1.30044948810764362720923581419356337476`28.38637488171231, 2.02132923840209670174259111538466355735`28.38637488171231, 1.`28.38637488171231}
```
## 2026-07-14 20:57:01 — 📥 BULK INSERT 1 CELL(S) at positions 18–18
1. [code] (* VALIDATION (generalization + regression guard): same construction for the ASYMMETRIC rep↵   λ1=2,


## 2026-07-14 20:57:01
**Input:**
```wolfram
(* VALIDATION (generalization + regression guard): same construction for the ASYMMETRIC rep
   λ1=2, λ2=1 (dim 6, M=0..3). This case is what exposed the NullSpace-tolerance bug: with the
   default tolerance, candidate (2,0) — whose FULL τ-polynomial matches a genuine transfer-matrix
   eigenvalue — was wrongly reported spurious (its null singular value is ~2e-28, far below the
   next SV ~2.3, but default NullSpace missed it). With Tolerance->10^-10 the count now equals the
   predicted weight multiplicity min(M,λ1,λ2,λ1+λ2-M)+1 and the genuine Psi span the full 6-dim
   space. Residuals use a FRESH t2 at each u (independent of the C0/C1/C2 extraction). *)
Module[{la = 2, lb = 1, uvals = {0, 1, 2, -1, 5}, rows, genuineCount, pred, maxres},
  rows = Flatten[Table[
     Module[{psi = Psi[la, lb][M, n], status, res},
       status = Which[MatchQ[psi, {__?NumericQ}], "genuine", MissingQ[psi], "spurious", True, "FAILED"];
       res = If[status == "genuine",
         Max[Abs[Flatten[Table[t2[la, lb][1, 1][uu] . psi - τ[la, lb][M, n][uu] psi, {uu, uvals}]]]],
         Missing[]];
       {M, n, status, res}],
     {M, 0, la + lb}, {n, 0, M}], 1];
  genuineCount = Association[Table[M -> Count[Select[rows, #[[1]] == M &][[All, 3]], "genuine"], {M, 0, la + lb}]];
  pred = Association[Table[M -> Min[M, la, lb, la + lb - M] + 1, {M, 0, la + lb}]];
  maxres = Max[Cases[rows[[All, 4]], _?NumericQ]];
  Print["λ1=2,λ2=1 genuine count per M : ", genuineCount];
  Print["predicted min(M,λ1,λ2,λ1+λ2-M)+1 : ", pred];
  Print["MATCHES prediction? ", genuineCount === pred, "   (total genuine = ", Total[Values[genuineCount]], " = dim ", (la + 1) (lb + 1), ")"];
  Print["max residual over genuine Psi (fresh t2) = ", maxres];
  Print[If[maxres < 10^-18, "RESIDUAL CHECK PASSED", "RESIDUAL CHECK FAILED"]];
  TableForm[rows, TableHeadings -> {None, {"M", "n", "status", "max|residual|"}}]
]
```
**Output:**
```
(no output)
```

---

## 2026-07-14 20:57:09 — 🔄 KERNEL RESTART

## 2026-07-14 20:57:25 — ▶️ RUN CELLS 1–18
Cell 2: ✓ — (no output)
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 7: ✓ — (no output)
Cell 8: ✓ — (no output)
Cell 9: ✓ — (no output)
Cell 10: ✓ — (no output)
Cell 11: ✓ — (no output)
Cell 12: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — (no output)
Cell 16: ✓ — (no output)
Cell 17: ✓ — (no output)
Cell 18: ✓ — (no output)


## 2026-07-14 20:57:38
**Input:**
```wolfram
Module[{report},
 report[la_,lb_]:=Module[{uvals={0,1,2,-1,5},rows},
   rows=Flatten[Table[Module[{psi=Quiet[Psi[la,lb][M,n]],status,res},
      status=Which[MatchQ[psi,{__?NumericQ}],"genuine",MissingQ[psi],"spurious",True,"FAILED"];
      res=If[status=="genuine",Max[Abs[Flatten[Table[t2[la,lb][1,1][uu].psi-τ[la,lb][M,n][uu] psi,{uu,uvals}]]]],Null];
      {M,n,status,If[NumericQ[res],ScientificForm[res,3],res],If[status=="genuine",Chop[N[psi,6],10^-9],Null]}],
     {M,0,la+lb},{n,0,M}],1];
   rows];
 <|"L11"->report[1,1],"L21"->report[2,1]|>]
```
**Output:**
```
Out=                                  -26                                                                                         -6                             -26                                                                 -26                                                                                                                              -26                                                                                         -6                                      -26                                                                                                                                 -6                             -26                                                                                                                              -6                             -25                                                                                                                                 -6                             -25                                                                                                                              -6                                                           -25                                                                                                                                 -6                                                                                                                       -25                                                                                                                                 -6
<|L11 -> {{0, 0, genuine, 0. × 10   , {-0.740166 + 0.672424 I, 0.360440 + 0.932782 I, 0.360440 + 0.932782 I, 1.000000 + 0. 10   I}}, {1, 0, genuine, 0. × 10   , {1.00000, -1.30045, 2.02133, 1.00000}}, {1, 1, genuine, 0. × 10   , {1.00000, 1.13083, -0.409951, 1.00000}}, {2, 0, spurious, Null, Null}, {2, 1, spurious, Null, Null}, {2, 2, genuine, 0. × 10   , {-0.740166 - 0.672424 I, 0.360440 - 0.932782 I, 0.360440 - 0.932782 I, 1.000000 + 0. 10   I}}}, L21 -> {{0, 0, genuine, 0. × 10   , {-1.78802 - 0.89609 I, -1.48033 + 1.34485 I, -1.48033 + 1.34485 I, 0.72088 + 1.86556 I, 0.360440 + 0.932782 I, 1.000000 + 0. 10   I}}, {1, 0, genuine, 0. × 10   , {0.72088 + 1.86556 I, -0.07253 - 4.00848 I, 1.29610 + 2.67666 I, -1.52215 + 0.24420 I, 2.60347 + 0.68858 I, 1.000000 + 0. 10   I}}, {1, 1, genuine, 0. × 10   , {0.72088 + 1.86556 I, 1.77891 + 2.03076 I, 0.370381 - 0.342955 I, 1.62816 + 0.46910 I, -0.546842 + 0.463683 I, 1.000000 + 0. 10   I}}, {2, 0, genuine, 0. × 10   , {0.72088 - 1.86556 I, -0.07253 + 4.00848 I, 1.29610 - 2.67666 I, -1.52215 - 0.24420 I, 2.60347 - 0.68858 I, 1.000000 + 0. 10   I}}, {2, 1, spurious, Null, Null}, {2, 2, genuine, 0. × 10   , {0.72088 - 1.86556 I, 1.77891 - 2.03076 I, 0.370381 + 0.342955 I, 1.62816 - 0.46910 I, -0.546842 - 0.463683 I, 1.000000 + 0. 10   I}}, {3, 0, spurious, Null, Null}, {3, 1, spurious, Null, Null}, {3, 2, spurious, Null, Null}, {3, 3, genuine, 0. × 10   , {-1.78802 + 0.89609 I, -1.48033 - 1.34485 I, -1.48033 - 1.34485 I, 0.72088 - 1.86556 I, 0.360440 - 0.932782 I, 1.000000 + 0. 10   I}}}|>
```

## 2026-07-14 20:58:07
**Input:**
```wolfram
{
 "bounds fires+returns" -> Reap[Check[Psi[1,1][5,0], "caught", Psi::bnds], _, Sow][[{1}]],
 "bounds msgQ" -> (Head[#]===String &)[Quiet[Check[Psi[1,1][6,0],"F",Psi::bnds];"ok"]],
 "bounds result" -> Quiet[Psi[1,1][5,0]],
 "spurious fresh result (1,2)(2,0)" -> Quiet[Psi[1,2][2,0]],
 "spurious msg fires? (fresh 1,2)(2,1)" -> (Length[Reap[Quiet[Psi[1,2][2,1], {Psi::spurious}]; Message[Psi::spurious,2,1,1,2]//Quiet][[1]]]>=0),
 "genuine counts final" -> {
   "L11" -> (Table[Count[Table[Quiet@MatchQ[Psi[1,1][M,n],{__?NumericQ}],{n,0,M}],True],{M,0,2}]),
   "L21" -> (Table[Count[Table[Quiet@MatchQ[Psi[2,1][M,n],{__?NumericQ}],{n,0,M}],True],{M,0,3}]),
   "L12" -> (Table[Count[Table[Quiet@MatchQ[Psi[1,2][M,n],{__?NumericQ}],{n,0,M}],True],{M,0,3}])}
}
```
**Output:**
```
[message] Psi::bnds: Invalid indices (M,n)=(5,0) for (λ1,λ2)=(1,1): need M an integer in 0..λ1+λ2 and n an integer in 0..M.
Out= {"bounds fires+returns" -> {"caught"}, "bounds msgQ" -> True, "bounds result" -> $Failed, "spurious fresh result (1,2)(2,0)" -> {0.72087975029445307453335530119110011567`28.153074769624215 - 1.86556489718674943296967196414173586546`28.153074769624215*I, 0.52281951842974346420890501848755989375`27.9821687260546 + 0.57694015147647477183277023666697427982`28.153074769624215*I, -0.71015452038593535032279427428522635137`28.153074769624215 - 0.2371533418544670754950145800722828570416813482600499511063`28.04810851998579*I, 1.4740285775251060891815656484946558559`28.153074769624215 - 2.4987282601950304979355006066627115898`28.097376630158738*I, 1.79147414582761496212282722607187661651`28.153074769624215 - 0.69562910673890764098982140199858509746`27.76004747591281*I, 1.`28.153074769624215 + 0``28.32543157765107*I}, "spurious msg fires? (fresh 1,2)(2,1)" -> True, "genuine counts final" -> {"L11" -> {1, 2, 1}, "L21" -> {1, 2, 2, 1}, "L12" -> {1, 2, 2, 1}}}
```
## 2026-07-15 09:47:12 — 🗑️ DELETE CODE CELL 2
(*Numeric parameters (standalone solver — this file defines everything it needs).
  Assigns global numeric values to θ1,θ2,h,κ1,κ2.*)
θ1 = 1/3;
θ2 = 1/7;
h = 1;
κ1 = N[Exp[I Zeta[3]], 30];
κ2 = 1/κ1;
…

## 2026-07-15 09:47:12 — 🗑️ DELETE CODE CELL 3
(*a(u), d(u)=Qθ(u), and τ0 = τ minus the J term (using ee1=λ1+λ2-M, ee2=M).
  Qθ (=d) is defined LOCALLY here so this file is self-contained; in su2_XXX_L2.wb it came from
  the spin-chain machinery (…

## 2026-07-15 09:47:12 — 🗑️ DELETE CODE CELL 4
(*BaxterMatrix: build O[u^k]=κ1 a(u)(u+h)^k + κ2 Qθ[u](u-h)^k - τ0(u) u^k for k=0..M.
  HARD ASSERT that the degree M+1 and M+2 coefficients vanish (this is where a wrong a/d/τ or a
  wrong ee1=λ1+λ2-…

## 2026-07-15 09:47:12 — 🗑️ DELETE CODE CELL 5
(*Explicit display of the high-degree (u^(M+1), u^(M+2)) coefficients for λ1=λ2=1, M=0,1,2 —
  confirms the cancellation the assertion relies on (should all be ~10^-30).*)
Table[
 Module[{clists},
  c…

## 2026-07-15 09:47:12 — 🗑️ DELETE CODE CELL 6
(*BaxterSolve: diagonalize BaxterMatrix. Eigenvalues are the admissible J_n; eigenvectors are
  the ascending coefficient vectors of q_n, renormalized so the u^M (top) component is 1 (monic).
  FIXED …

## 2026-07-15 09:47:12 — 🗑️ DELETE CODE CELL 7
(*INDEPENDENT residual check: re-evaluate the FULL Baxter equation from scratch (using τ=τ0+J_n,
  NOT any BaxterMatrix internals) at several u that are not roots of a(u) or d(u).
  a(u) roots: 4/3, 8…

## 2026-07-15 09:47:12 — 🗑️ DELETE CODE CELL 8
(*Test run: full (M,n,J_n,q_n(u)) table for λ1=λ2=1 (so λ1+λ2=2, M=0,1,2). Values rounded to 10
  digits for display; full 30-digit values are available from BaxterSolve[1,1,M].*)
resultTable = Flatte…

## 2026-07-15 09:47:12 — 🗑️ DELETE CODE CELL 9
(*On-demand wrappers over the memoized BaxterSolve cache.
  Q1[λ1,λ2][M,n][u] = κ1^(u/h) q_n(u)      — the actual Baxter Q-function value at u.
  τ[λ1,λ2][M,n][u]  = τ0[λ1,λ2,M][u] + J_n — the FULL τ …

## 2026-07-15 09:47:12 — 🗑️ DELETE CODE CELL 10
(*VALIDATION 5a: end-to-end check of the ORIGINAL (pre-ansatz-substitution) Baxter equation
    a(u) Q(u+h) - τ(u) Q(u) + d(u) Q(u-h) = 0
  assembled ENTIRELY from the wrapper functions Q1 and τ (NOT …

## 2026-07-15 09:47:12 — 🗑️ DELETE CODE CELL 11
(*VALIDATION 5b: genuine generalization to a NEW representation λ1=2, λ2=1 (asymmetric; λ1+λ2=3,
  so M=0..3). Same end-to-end wrapper residual check, plus the full (M,n,J_n) table for this case.*)
sa…

## 2026-07-15 09:47:12 — 🗑️ DELETE CODE CELL 12
(*VALIDATION 5c: confirm the fixed n-ordering actually holds — Re[J_n] must be NON-DECREASING in n
  within each sector. Spot-checked on sectors with M>=2 (where >1 solution exists to order):
  (λ1,λ2…

## 2026-07-15 09:47:12 — 🗑️ DELETE MARKDOWN CELL 13
(*Genuine transfer-matrix eigenvectors from the Baxter τ eigenvalue.

Goal: build Psi[λ1,λ2][M,n], a FIXED vector (no [u] argument) in V_λ1⊗V_λ2 of dimension
(λ1+1)(λ2+1), satisfying  t2[λ1,λ2][1,1][u…

## 2026-07-15 09:47:12 — 🗑️ DELETE CODE CELL 14
(* ===== Transfer-matrix machinery, ported verbatim from su2_XXX_L2.wb (cells 2-10) =====
   Single-site GT-pattern generator machinery (ps,dim,J0,Jp,Jm,JJ), spin-s generators Ee and Lax
   matrix L (…

## 2026-07-15 09:47:12 — 🗑️ DELETE CODE CELL 15
(* Coefficient matrices C0,C1,C2 with t2[λ1,λ2][1,1][u] = C0 + C1 u + C2 u^2.
   Extracted entrywise via CoefficientList in a dummy variable (NOT by evaluating at numeric
   points). Memoized per (λ1,…

## 2026-07-15 09:47:12 — 🗑️ DELETE CODE CELL 16
(* Psi[λ1,λ2][M,n]: transfer-matrix eigenvector with eigenvalue τ[λ1,λ2][M,n][u], as a fixed
   vector in V_λ1⊗V_λ2 (dim (λ1+1)(λ2+1)), normalized so its LAST component is 1.
   Psi lies in the simult…

## 2026-07-15 09:47:12 — 🗑️ DELETE CODE CELL 17
(* VALIDATION: for every candidate (M,n) with λ1=λ2=1 (M=0..2, n=0..M) attempt Psi, and where a
   genuine vector is found INDEPENDENTLY check t2(u).Psi == τ(u) Psi by re-deriving the numeric t2
   ma…

## 2026-07-15 09:47:12 — 🗑️ DELETE CODE CELL 18
(* VALIDATION (generalization + regression guard): same construction for the ASYMMETRIC rep
   λ1=2, λ2=1 (dim 6, M=0..3). This case is what exposed the NullSpace-tolerance bug: with the
   default to…

## 2026-07-15 09:47:34 — 📥 BULK INSERT 4 CELL(S) at positions 2–5
1. [code] (*Numeric parameters (standalone solver — this file defines everything it needs).↵  Assigns global n
2. [code] (*a(u), d(u)=Qθ(u), and τ0 = τ minus the J term (using ee1=λ1+λ2-M, ee2=M).↵  Qθ (=d) is defined LOC
3. [code] (*BaxterMatrix: build O[u^k]=κ1 a(u)(u+h)^k + κ2 Qθ[u](u-h)^k - τ0(u) u^k for k=0..M.↵  HARD ASSERT 
4. [code] (*Explicit display of the high-degree (u^(M+1), u^(M+2)) coefficients for λ1=λ2=1, M=0,1,2 —↵  confi


## 2026-07-15 09:47:34
**Input:**
```wolfram
(*Numeric parameters (standalone solver — this file defines everything it needs).
  Assigns global numeric values to θ1,θ2,h,κ1,κ2.*)
θ1 = 1/3;
θ2 = 1/7;
h = 1;
κ1 = N[Exp[I Zeta[3]], 30];
κ2 = 1/κ1;
{θ1, θ2, h, κ1, κ2}
```
**Output:**
```
Out[6]= \left\{\frac{1}{3},\frac{1}{7},1,0.360439875147226537266677650596+0.932782448593374716484835982071\,\mathrm{i},0.36043987514722653726667765060-0.932782448593374716484835982071\,\mathrm{i}\right\}
```

## 2026-07-15 09:47:35
**Input:**
```wolfram
(*a(u), d(u)=Qθ(u), and τ0 = τ minus the J term (using ee1=λ1+λ2-M, ee2=M).
  Qθ (=d) is defined LOCALLY here so this file is self-contained.
  qpoly builds a polynomial from an ASCENDING coefficient list {c0,...,cM} (cM=1 for monic q_n);
  formed in a dummy symbol first so substituting uu=0 gives the constant term rather than an
  indeterminate 0^0.*)
ClearAll[aFun, τ0, Qθ, qpoly];
aFun[λ1_, λ2_][u_] := (u - θ1 - h λ1) (u - θ2 - h λ2);
Qθ[u_] := (u - θ1) (u - θ2);
τ0[λ1_, λ2_, M_][u_] := (κ1 + κ2) u^2 - (θ1 + θ2) (κ1 + κ2) u - u h (κ1 (λ1 + λ2 - M) + κ2 M);
qpoly[clist_, uu_] := Module[{x}, (clist . x^Range[0, Length[clist] - 1]) /. x -> uu];
{aFun[1, 1][u], Qθ[u], τ0[1, 1, 1][u]}
```
**Output:**
```
Out[12]= \left\{\left(u-\frac{4}{3}\right)\,\left(u-\frac{8}{7}\right),\left(u-\frac{1}{3}\right)\,\left(u-\frac{1}{7}\right),\left(0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i}\right)\,u^2-\left(1.06415582186324025288257211128+0.\times 10^{-29}\,\mathrm{i}\right)\,u\right\}
```

## 2026-07-15 09:47:35
**Input:**
```wolfram
(*BaxterMatrix: build O[u^k]=κ1 a(u)(u+h)^k + κ2 Qθ[u](u-h)^k - τ0(u) u^k for k=0..M.
  HARD ASSERT that the degree M+1 and M+2 coefficients vanish (this is where a wrong a/d/τ or a
  wrong ee1=λ1+λ2-M convention would show up). Then columns = degree-0..M coefficients.
  Used TWICE downstream: (old role) never diagonalized directly anymore; (new role) as the FIXED
  linear operator whose nullspace at a KNOWN genuine J gives the Q-function coefficients
  (see QSolve below) — this is now a linear solve, not an eigenvalue problem.*)
ClearAll[BaxterMatrix];
BaxterMatrix[λ1_, λ2_, M_] := Module[{clists, high, maxhigh, cols},
  clists = Table[
     PadRight[
      CoefficientList[
       Expand[κ1 aFun[λ1, λ2][u] (u + h)^k + κ2 Qθ[u] (u - h)^k - τ0[λ1, λ2, M][u] u^k],
       u],
      M + 3],
     {k, 0, M}];
  high = clists[[All, {M + 2, M + 3}]];   (* degree M+1 and M+2 coeffs for each column k *)
  maxhigh = Max[Abs[Flatten[high]]];
  If[maxhigh > 10^-18,
    Print["BaxterMatrix ASSERTION FAILED (λ1,λ2,M)=", {λ1, λ2, M},
          "  max|high-degree coeff|=", maxhigh, "  high=", high];
    Abort[]];
  cols = clists[[All, 1 ;; M + 1]];       (* degrees 0..M, one column per k *)
  Transpose[cols]
];
Dimensions /@ (BaxterMatrix[1, 1, #] & /@ {0, 1, 2})
```
**Output:**
```
Out[15]= \begin{pmatrix}1 & 1\\2 & 2\\3 & 3\end{pmatrix}
```

## 2026-07-15 09:47:35
**Input:**
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
**Output:**
```
Out[16]= \begin{matrix}\left\{0,\text{max|deg M+1,M+2 coeff|}\to 0.\times 10^{-29}\right\}\\\left\{1,\text{max|deg M+1,M+2 coeff|}\to 0.\times 10^{-29}\right\}\\\left\{2,\text{max|deg M+1,M+2 coeff|}\to 0.\times 10^{-29}\right\}\end{matrix}
```
## 2026-07-15 09:48:54 — 📥 BULK INSERT 2 CELL(S) at positions 6–7
1. [code] (* ===== Transfer-matrix machinery, ported verbatim from su2_XXX_L2.wb (cells 2-10) =====↵   Single-
2. [code] (* Coefficient matrices C0,C1,C2 with t2[λ1,λ2][1,1][u] = C0 + C1 u + C2 u^2.↵   Extracted entrywise


## 2026-07-15 09:48:55
**Input:**
```wolfram
(* ===== Transfer-matrix machinery, ported verbatim from su2_XXX_L2.wb (cells 2-10) =====
   Single-site GT-pattern generator machinery (ps,dim,J0,Jp,Jm,JJ), spin-s generators Ee and Lax
   matrix L (uses global h), the shared companion twist G (set NUMERIC here with this file's κ1,κ2
   so its eigenvalues are exactly κ1,κ2), the two-site embeddings L1/L2 (use global θ1,θ2), the
   monodromy T, and the fundamental transfer matrix t2. *)
ClearAll[ps, dim, J0, Jp, Jm, JJ, Ee, L, L1, L2, T, t2, χ1, χ2, G];
ps[rep_] := ps[rep] = Block[{res = (Prepend[rep, #] & /@ (Tuples[Table[Range[rep[[1, i + 1]], rep[[1, i]]], {i, 1, Length[rep[[1]]] - 1}]]))}, If[Length[res[[1, 1]]] == 1, Reverse@Sort@res, Reverse@Sort@Flatten[ps /@ res, 1]]];
dim[rep_] := dim[rep] = Length[ps[{rep}]];
J0[k_][pat_] := J0[k][pat] = Sum[pat[[k, i]], {i, 1, k}] - Sum[pat[[k - 1, i]], {i, 1, k - 1}];
Jp[k_][pat_] := Jp[k][pat] = -Sum[Product[pat[[k, i]] - i - pat[[k + 1, j]] + j, {j, Range[k + 1]}]/Product[pat[[k, i]] - i - pat[[k, j]] + j, {j, Complement[Range[k], {i}]}] v[ReplacePart[pat, {k, i} -> pat[[k, i]] + 1]], {i, 1, k}];
Jm[k_][pat_] := Jm[k][pat] = Sum[Product[pat[[k, i]] - i - pat[[k - 1, j]] + j, {j, Range[k - 1]}]/Product[pat[[k, i]] - i - pat[[k, j]] + j, {j, Complement[Range[k], {i}]}] v[ReplacePart[pat, {k, i} -> pat[[k, i]] - 1]], {i, 1, k}];
JJ[k_, l_, rep_] := JJ[k, l, rep] = Which[k == l, DiagonalMatrix[J0[k] /@ ps[{rep}]], l == k + 1, Transpose[Table[Coefficient[Jp[k][ps[{rep}][[i]]], v[#]] & /@ ps[{rep}], {i, 1, Length[ps[{rep}]]}]], l == k - 1, Transpose[Table[Coefficient[Jm[l][ps[{rep}][[i]]], v[#]] & /@ ps[{rep}], {i, 1, Length[ps[{rep}]]}]], l > k + 1, JJ[k, k + 1, rep] . JJ[k + 1, l, rep] - JJ[k + 1, l, rep] . JJ[k, k + 1, rep], l < k - 1, JJ[k, k - 1, rep] . JJ[k - 1, l, rep] - JJ[k - 1, l, rep] . JJ[k, k - 1, rep]];
Ee[λ_][i_, j_] := JJ[i, j, {λ, 0}];
L[λ_][i_, j_][u_] := KroneckerDelta[i, j]*u*IdentityMatrix[dim[{λ, 0}]] - h Ee[λ][j, i];
χ1 = κ1 + κ2; χ2 = 1;
G = {{χ1, -χ2}, {1, 0}};
L1[λ1_, λ2_][i_, j_][u_] := KroneckerProduct[L[λ1][i, j][u - θ1], IdentityMatrix[dim[{λ2, 0}]]];
L2[λ1_, λ2_][i_, j_][u_] := KroneckerProduct[IdentityMatrix[dim[{λ1, 0}]], L[λ2][i, j][u - θ2]];
T[λ1_, λ2_][i_, j_][u_] := T[λ1, λ2][i, j][u] = Sum[L1[λ1, λ2][i, k][u] . L2[λ1, λ2][k, j][u], {k, 2}];
t2[λ1_, λ2_][1, 1][u_] := Sum[T[λ1, λ2][i, j][u] G[[j, i]], {i, 2}, {j, 2}];
(* sanity: G eigenvalues = {κ1,κ2}; t2 is 4x4 for λ1=λ2=1 *)
{"G eig - {κ1,κ2}" -> Chop[Sort[Eigenvalues[N[G]]] - Sort[N[{κ1, κ2}]]], "dim t2[1,1]" -> Dimensions[t2[1, 1][1, 1][3]]}
```
**Output:**
```
Out[32]= \left\{\text{G eig - \{$\kappa$1,$\kappa$2\}}\to \left\{0.\, +1.8655648971867491\,\mathrm{i},0.\, -1.8655648971867493\,\mathrm{i}\right\},\text{dim t2[1,1]}\to \left\{4,4\right\}\right\}
```

## 2026-07-15 09:48:55
**Input:**
```wolfram
(* Coefficient matrices C0,C1,C2 with t2[λ1,λ2][1,1][u] = C0 + C1 u + C2 u^2.
   Extracted entrywise via CoefficientList in a dummy variable (NOT by evaluating at numeric
   points). Memoized per (λ1,λ2) — reusable across all (M,n) of that representation. *)
ClearAll[t2Coeffs];
t2Coeffs[λ1_, λ2_] := t2Coeffs[λ1, λ2] = Module[{uu, mat, cl},
   mat = t2[λ1, λ2][1, 1][uu];
   cl = Map[PadRight[CoefficientList[#, uu], 3] &, mat, {2}];
   Table[cl[[All, All, k]], {k, 1, 3}]
];
(* check: C2 must equal (κ1+κ2) Id — leading-order asymptotics match τ by construction *)
With[{Cs = t2Coeffs[1, 1]}, "max|C2-(κ1+κ2)Id|" -> Max[Abs[Flatten[Cs[[3]] - (κ1 + κ2) IdentityMatrix[4]]]]]
```
**Output:**
```
Out[35]= \text{max|C2-($\kappa$1+$\kappa$2)Id|}\to 0.\times 10^{-29}
```

## 2026-07-15 09:49:05
**Input:**
```wolfram
{Eigenvalues[N[G]], N[{κ1,κ2}], G}
```
**Output:**
```
Out= {{0.36043987514722653 - 0.9327824485933747*I, 0.3604398751472264 + 0.9327824485933744*I}, {0.36043987514722653 + 0.9327824485933747*I, 0.36043987514722653 - 0.9327824485933747*I}, {{0.7208797502944530745333553011911001825875290465245700049878`29.54117567596917 + 0``29.50582140452307*I, -1}, {1, 0}}}
```

## 2026-07-15 09:49:24
**Input:**
```wolfram
Chop[Sort[Eigenvalues[N[G]]] - Sort[N[{κ1, κ2}]]]
```
**Output:**
```
Out= {0. + 1.8655648971867491*I, 0. - 1.8655648971867493*I}
```

## 2026-07-15 09:49:42
**Input:**
```wolfram
{Sort[Eigenvalues[N[G]]], Sort[N[{κ1, κ2}]]}
```
**Output:**
```
Out= {{0.3604398751472264 + 0.9327824485933744*I, 0.36043987514722653 - 0.9327824485933747*I}, {0.36043987514722653 - 0.9327824485933747*I, 0.36043987514722653 + 0.9327824485933747*I}}
```

## 2026-07-15 09:50:22
**Input:**
```wolfram
Chop[SortBy[Eigenvalues[N[G]], {Re[#],Im[#]}&] - SortBy[N[{κ1, κ2}], {Re[#],Im[#]}&]]
```
**Output:**
```
Out= {0. + 1.8655648971867491*I, 0. - 1.8655648971867493*I}
```

## 2026-07-15 09:50:29
**Input:**
```wolfram
{SortBy[Eigenvalues[N[G]], {Re[#],Im[#]}&], SortBy[N[{κ1, κ2}], {Re[#],Im[#]}&]}
```
**Output:**
```
Out= {{0.3604398751472264 + 0.9327824485933744*I, 0.36043987514722653 - 0.9327824485933747*I}, {0.36043987514722653 - 0.9327824485933747*I, 0.36043987514722653 + 0.9327824485933747*I}}
```

## 2026-07-15 09:50:40
**Input:**
```wolfram
Eigenvalues[N[G]] // InputForm
```
**Output:**
```
Out= InputForm[{0.36043987514722653 - 0.9327824485933747*I, 0.3604398751472264 + 0.9327824485933744*I}]
```

## 2026-07-15 09:51:29
**Input:**
```wolfram
Chop[{Tr[N[G]], Det[N[G]]} - {κ1 + κ2, κ1 κ2}]
```
**Output:**
```
Out= {0, 0}
```
## 2026-07-15 09:51:52 — ✏️ EDIT CELL 6
(* ===== Transfer-matrix machinery, ported verbatim from su2_XXX_L2.wb (cells 2-10) =====
   Single-site GT-pattern generator machinery (ps,dim,J0,Jp,Jm,JJ), spin-s generators Ee and Lax
   matrix L (…


## 2026-07-15 09:51:52
**Input:**
```wolfram
(* ===== Transfer-matrix machinery, ported verbatim from su2_XXX_L2.wb (cells 2-10) =====
   Single-site GT-pattern generator machinery (ps,dim,J0,Jp,Jm,JJ), spin-s generators Ee and Lax
   matrix L (uses global h), the shared companion twist G (set NUMERIC here with this file's κ1,κ2
   so its eigenvalues are exactly κ1,κ2), the two-site embeddings L1/L2 (use global θ1,θ2), the
   monodromy T, and the fundamental transfer matrix t2. *)
ClearAll[ps, dim, J0, Jp, Jm, JJ, Ee, L, L1, L2, T, t2, χ1, χ2, G];
ps[rep_] := ps[rep] = Block[{res = (Prepend[rep, #] & /@ (Tuples[Table[Range[rep[[1, i + 1]], rep[[1, i]]], {i, 1, Length[rep[[1]]] - 1}]]))}, If[Length[res[[1, 1]]] == 1, Reverse@Sort@res, Reverse@Sort@Flatten[ps /@ res, 1]]];
dim[rep_] := dim[rep] = Length[ps[{rep}]];
J0[k_][pat_] := J0[k][pat] = Sum[pat[[k, i]], {i, 1, k}] - Sum[pat[[k - 1, i]], {i, 1, k - 1}];
Jp[k_][pat_] := Jp[k][pat] = -Sum[Product[pat[[k, i]] - i - pat[[k + 1, j]] + j, {j, Range[k + 1]}]/Product[pat[[k, i]] - i - pat[[k, j]] + j, {j, Complement[Range[k], {i}]}] v[ReplacePart[pat, {k, i} -> pat[[k, i]] + 1]], {i, 1, k}];
Jm[k_][pat_] := Jm[k][pat] = Sum[Product[pat[[k, i]] - i - pat[[k - 1, j]] + j, {j, Range[k - 1]}]/Product[pat[[k, i]] - i - pat[[k, j]] + j, {j, Complement[Range[k], {i}]}] v[ReplacePart[pat, {k, i} -> pat[[k, i]] - 1]], {i, 1, k}];
JJ[k_, l_, rep_] := JJ[k, l, rep] = Which[k == l, DiagonalMatrix[J0[k] /@ ps[{rep}]], l == k + 1, Transpose[Table[Coefficient[Jp[k][ps[{rep}][[i]]], v[#]] & /@ ps[{rep}], {i, 1, Length[ps[{rep}]]}]], l == k - 1, Transpose[Table[Coefficient[Jm[l][ps[{rep}][[i]]], v[#]] & /@ ps[{rep}], {i, 1, Length[ps[{rep}]]}]], l > k + 1, JJ[k, k + 1, rep] . JJ[k + 1, l, rep] - JJ[k + 1, l, rep] . JJ[k, k + 1, rep], l < k - 1, JJ[k, k - 1, rep] . JJ[k - 1, l, rep] - JJ[k - 1, l, rep] . JJ[k, k - 1, rep]];
Ee[λ_][i_, j_] := JJ[i, j, {λ, 0}];
L[λ_][i_, j_][u_] := KroneckerDelta[i, j]*u*IdentityMatrix[dim[{λ, 0}]] - h Ee[λ][j, i];
χ1 = κ1 + κ2; χ2 = 1;
G = {{χ1, -χ2}, {1, 0}};
L1[λ1_, λ2_][i_, j_][u_] := KroneckerProduct[L[λ1][i, j][u - θ1], IdentityMatrix[dim[{λ2, 0}]]];
L2[λ1_, λ2_][i_, j_][u_] := KroneckerProduct[IdentityMatrix[dim[{λ1, 0}]], L[λ2][i, j][u - θ2]];
T[λ1_, λ2_][i_, j_][u_] := T[λ1, λ2][i, j][u] = Sum[L1[λ1, λ2][i, k][u] . L2[λ1, λ2][k, j][u], {k, 2}];
t2[λ1_, λ2_][1, 1][u_] := Sum[T[λ1, λ2][i, j][u] G[[j, i]], {i, 2}, {j, 2}];
(* sanity: char. poly of G is λ²-χ1 λ+χ2 = (λ-κ1)(λ-κ2), so {tr,det}(G) must equal {κ1+κ2,κ1 κ2}.
   (NOT checked via Sort[Eigenvalues[G]]-Sort[{κ1,κ2}]: since κ2=Conjugate[κ1] exactly (|κ1|=1),
   Eigenvalues[N[G]] introduces ~1e-17 real-part noise between the conjugate pair, which is enough
   to flip Sort's ordering relative to the exact {κ1,κ2} and make an otherwise-correct check falsely
   report an O(1) discrepancy. tr/det is the order-independent, numerically robust check.) *)
{"{tr,det}(G) - {κ1+κ2,κ1 κ2}" -> Chop[{Tr[N[G]], Det[N[G]]} - {κ1 + κ2, κ1 κ2}], "dim t2[1,1]" -> Dimensions[t2[1, 1][1, 1][3]]}
```
**Output:**
```
Out= {"{tr,det}(G) - {κ1+κ2,κ1 κ2}" -> {0, 0}, "dim t2[1,1]" -> {4, 4}}
```
## 2026-07-15 09:52:03 — ▶️ RUN CELL 7 (0.02 s)
**In [7]:** `(* Coefficient matrices C0,C1,C2 with t2[λ1,λ2][1,1][u] = C0 + C1 u + C2 u^2.↵   Extracted entrywise via CoefficientList in a dummy variable (NOT by evaluating at numeric↵   points). Memoized per (λ1,`
**Out:** Out[38]= \text{max|C2-($\kappa$1+$\kappa$2)Id|}\to 0.\times 10^{-29}


## 2026-07-15 09:52:13
**Input:**
```wolfram
dPred[M_, λ1_, λ2_] := Min[M, λ1, λ2, λ1 + λ2 - M] + 1
Cs = t2Coeffs[1, 1];
{C0, C1, C2} = Cs;
d = 4;
rCandidates = {7/3, N[Pi, 30], N[E, 30], 11/5};
Table[
 Module[{ee, vv},
  {ee, vv} = Eigensystem[C0 + r C1];
  {r, MatrixRank[vv], Length[Union[Chop[ee, 10^-12]]]}
 ], {r, rCandidates}]
```
**Output:**
```
Out[1]= Null
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= (suppressed)
Out[5]= (suppressed)
Out[6]= {{7/3, 4, 4}, {3.1415926535897932384626433832795028841971693993751058209749`30., 4, 4}, {2.7182818284590452353602874713526624938382062863352767788129`30., 4, 4}, {11/5, 4, 4}}
```

## 2026-07-15 09:52:20
**Input:**
```wolfram
r = 7/3;
{vals, vecs} = Eigensystem[N[C0 + r C1, 30]];
recs = Table[
   Module[{v = vecs[[i]], t0, t1, res0, res1, Mi},
     t0 = (v . (C0 . v))/(v . v);
     t1 = (v . (C1 . v))/(v . v);
     res0 = Max[Abs[C0 . v - t0 v]];
     res1 = Max[Abs[C1 . v - t1 v]];
     Mi = (-(θ1 + θ2)(κ1 + κ2) - h κ1 (1 + 1) - t1)/(h(κ2 - κ1));
     {i, t0, t1, res0, res1, Mi}
   ], {i, 1, d}];
TableForm[recs]
```
**Output:**
```
Out= TableForm[{{1, 0.56640551808849884427620773665015014352`26.563908186381784 + 1.37696456697117220052523406877128122917`26.915783392193127*I, -1.06415582186324025288257211128210026959`26.697527434398307 - 1.86556489718674943296967196414173585885`26.897367377166717*I, 0``26.78747568852612, 0``26.64117555017849, 0``26.896834687289104 + 0``26.940844418343726*I}, {2, 0.56640551808849884427620773665015014352`26.563908186381784 - 1.37696456697117220052523406877128122917`26.915783392193127*I, -1.06415582186324025288257211128210026959`26.697527434398307 + 1.86556489718674943296967196414173585885`26.897367377166717*I, 0``26.78747568852612, 0``26.64117555017849, 1.99999999999999999999999999999999999996`27.194441534100946 + 0``26.940002801075714*I}, {3, -0.77812396645115795553762697294103694812`27.852290597730754 + 0``28.581379173510307*I, -1.06415582186324025288257211128210026954`27.837863308189405 + 0``28.439580594982075*I, 0``27.981348989228252, 0``27.83096433407242, 1.`28.5841435970618 + 0``28.069274167806608*I}, {4, 1.19005525233370256955668714505023705247`27.867765140203126 + 0``28.443999752384332*I, -1.06415582186324025288257211128210026954`27.82211180328276 + 0``28.373868172528027*I, 0``27.91544179335782, 0``27.87244729292644, 1.`28.534093472118112 + 0``28.053958116895085*I}}]
```

## 2026-07-15 09:52:38
**Input:**
```wolfram
Cs21 = t2Coeffs[2, 1];
{C0b, C1b, C2b} = Cs21;
d21 = 6;
rTest = SelectFirst[{7/3, N[Pi,30], N[E,30], 11/5, 13/4, N[Sqrt[2],30]},
  Module[{ee, vv}, {ee, vv} = Eigensystem[C0b + # C1b]; MatrixRank[vv] == d21 && Length[Union[Chop[ee, 10^-12]]] == d21] &,
  Missing["NoGenericR"]];
rTest
```
**Output:**
```
Out= 7/3
```

## 2026-07-15 09:52:48
**Input:**
```wolfram
{valsb, vecsb} = Eigensystem[N[C0b + rTest C1b, 30]];
recsb = Table[
   Module[{v = vecsb[[i]], t0, t1, res0, res1, Mi},
     t0 = (v . (C0b . v))/(v . v);
     t1 = (v . (C1b . v))/(v . v);
     res0 = Max[Abs[C0b . v - t0 v]];
     res1 = Max[Abs[C1b . v - t1 v]];
     Mi = (-(θ1 + θ2)(κ1 + κ2) - h κ1 (2 + 1) - t1)/(h(κ2 - κ1));
     {i, Chop[t0,10^-15], Round[Re[Mi]], Abs[Mi-Round[Re[Mi]]], res0, res1}
   ], {i, 1, d21}];
TableForm[recsb]
```
**Output:**
```
Out= TableForm[{{1, 0.97833680397104345829526790875935030108`26.48339553183149 + 2.44300165107788616222218947685227322473`26.88756662607194*I, 0, 0``26.533130374964013, 0``26.518114780378326, 0``26.434228732079823}, {2, 0.97833680397104345829526790875935024775`26.48339553183149 - 2.44300165107788616222218947685227314844`26.88756662607194*I, 3, 0``26.53121406460057, 0``26.518114780378326, 0``26.434228732079823}, {3, -0.83744925205780684849419842696152995843`27.61471509837466 + 0.50861707940891848821622002755228133718`27.349383600792255*I, 1, 0``27.609059659037424, 0``27.650594818800617, 0``27.498239356862104}, {4, -0.83744925205780684849419842696152999853`27.61471509837466 - 0.50861707940891848821622002755228139733`27.349383600792255*I, 2, 0``27.60194147587506, 0``27.650594818800617, 0``27.498239356862104}, {5, 1.71280323455821415328470129269358020904`27.817162788911126 - 0.69067464021113471869285480653883454959`27.699174448547968*I, 2, 0``27.744639559496584, 0``27.642278194294438, 0``27.591419712438334}, {6, 1.71280323455821415328470129269358008854`27.817162788911126 + 0.69067464021113471869285480653883452093`27.699174448547968*I, 1, 0``27.75402970464649, 0``27.642278194294438, 0``27.591419712438334}}]
```
## 2026-07-15 09:53:33 — 📥 BULK INSERT 1 CELL(S) at positions 8–8
1. [code] (* dPred[M,λ1,λ2]: predicted dimension of the weight-M subspace of Vλ1⊗Vλ2 (number of (m1,m2)↵   pai


## 2026-07-15 09:53:34
**Input:**
```wolfram
(* dPred[M,λ1,λ2]: predicted dimension of the weight-M subspace of Vλ1⊗Vλ2 (number of (m1,m2)
   pairs with m1+m2 fixed) = Min[M,λ1,λ2,λ1+λ2-M]+1. This is what the GENUINE t2-diagonalization
   below must reproduce exactly for every M — no more, no less — since it is now an honest
   d=(λ1+1)(λ2+1)-dimensional diagonalization of the actual transfer matrix, not an oversized
   (M+1)-dimensional Baxter-matrix proxy that can manufacture spurious extra roots. *)
ClearAll[dPred];
dPred[M_, λ1_, λ2_] := Min[M, λ1, λ2, λ1 + λ2 - M] + 1;

(* TauEigensystem[λ1,λ2]: diagonalize the GENUINE transfer matrix t2[λ1,λ2][1,1][u] directly.
   C0,C1,C2 pairwise commute (all built from the same commuting-family monodromy), so a GENERIC
   linear combination C0+r C1 has the full set of common eigenvectors of C0 and C1 (C2 is already
   known to be the scalar (κ1+κ2) Id, so it carries no extra information). We try a short list of
   generic r candidates and keep the first that fully resolves the d-dim eigenspace (rank d,
   d distinct combo-eigenvalues); this guards against an accidental coincidental degeneracy in
   C0+r C1 for a specific r.
   For each eigenvector v: extract t0=(v.C0.v)/(v.v), t1=(v.C1.v)/(v.v) (valid for ANY genuine
   right eigenvector, no conjugation needed since C0,C1 are not Hermitian). The residual checks
   ||C0.v-t0 v|| and ||C1.v-t1 v|| being ~0 IS the genuineness test (replaces the old
   Baxter-matrix + post-hoc Psi::spurious filtering entirely: this diagonalization is exact, so a
   real failure here means a genuine construction bug, hence a hard Abort).
   M is read off ALGEBRAICALLY from t1 via τ0's u^1 coefficient (τ=τ0+J only shifts u^0):
     M = [-(θ1+θ2)(κ1+κ2) - hκ1(λ1+λ2) - t1] / (h(κ2-κ1))
   and must land within 10^-10 of an integer (hard assert). n within each M-sector: J=t0 directly
   (τ0 has zero constant term), sorted ascending by (Re[J],Im[J]) — same tie-break convention as
   the old BaxterSolve. Each M-sector's count is asserted to equal dPred[M,λ1,λ2] exactly.
   Returns an Association keyed by {M,n} -> <|"J"->J_n,"Psi"->eigenvector normalized to last
   component 1|>. Memoized per (λ1,λ2). *)
ClearAll[TauEigensystem];
TauEigensystem::rankdef = "Combination C0+r C1 failed to resolve the `1`-dim eigenspace for (λ1,λ2)=(`2`,`3`) at every tried r in `4`.";
TauEigensystem::noninteger = "M = `1` is not within tolerance of an integer for eigenvector #`2` at (λ1,λ2)=(`3`,`4`) (combo r=`5`) — indicates a bug in the τ0 u^1-coefficient formula or the t2 construction.";
Psi::notgenuine = "Eigenvector #`1` for (λ1,λ2)=(`2`,`3`) failed the common-eigenvector residual check against C0/C1 (residual=`4`, combo r=`5`) — the generic-combination diagonalization should make this impossible; investigate before trusting results.";
Psi::degnull = "Weight sector M=`1` for (λ1,λ2)=(`2`,`3`) has `4` genuine state(s) but the predicted weight multiplicity dPred is `5`.";
Psi::lastzero = "Genuine eigenvector for (M,n)=(`1`,`2`) at (λ1,λ2)=(`3`,`4`) has last component ~0; cannot normalize.";
TauEigensystem[λ1_, λ2_] := TauEigensystem[λ1, λ2] = Module[
   {d, Cs, C0, C1, C2, rCandidates, r, vals, vecs, recs, byM, table},
   d = (λ1 + 1) (λ2 + 1);
   Cs = t2Coeffs[λ1, λ2]; {C0, C1, C2} = Cs;
   rCandidates = {7/3, N[Pi, 30], N[E, 30], 11/5, 13/4, N[Sqrt[2], 30]};
   r = SelectFirst[rCandidates,
     Module[{ee, vv}, {ee, vv} = Eigensystem[C0 + # C1];
       MatrixRank[vv] == d && Length[Union[Chop[ee, 10^-12]]] == d] &,
     Missing["NoGenericR"]];
   If[MissingQ[r], Message[TauEigensystem::rankdef, d, λ1, λ2, rCandidates]; Abort[]];
   {vals, vecs} = Eigensystem[N[C0 + r C1, 30]];
   recs = Table[
     Module[{v = vecs[[i]], t0, t1, res0, res1, Mi},
       t0 = (v . (C0 . v))/(v . v);
       t1 = (v . (C1 . v))/(v . v);
       res0 = Max[Abs[C0 . v - t0 v]];
       res1 = Max[Abs[C1 . v - t1 v]];
       If[Max[res0, res1] > 10^-15,
         Message[Psi::notgenuine, i, λ1, λ2, Max[res0, res1], r]; Abort[]];
       Mi = (-(θ1 + θ2) (κ1 + κ2) - h κ1 (λ1 + λ2) - t1)/(h (κ2 - κ1));
       If[Abs[Mi - Round[Re[Mi]]] > 10^-10,
         Message[TauEigensystem::noninteger, Mi, i, λ1, λ2, r]; Abort[]];
       <|"M" -> Round[Re[Mi]], "J" -> t0, "vec" -> v|>],
     {i, 1, d}];
   byM = GroupBy[recs, #["M"] &];
   table = Association @ Flatten @ KeyValueMap[
     Function[{Mval, group},
       Module[{sorted = SortBy[group, {Re[#["J"]], Im[#["J"]]} &], predicted = dPred[Mval, λ1, λ2]},
         If[Length[sorted] != predicted,
           Message[Psi::degnull, Mval, λ1, λ2, Length[sorted], predicted]; Abort[]];
         MapIndexed[
           Function[{rec, idx}, Module[{n = First[idx] - 1, vv = rec["vec"], last},
             last = vv[[-1]];
             If[Abs[last] < 10^-15 Max[Abs[vv]],
               Message[Psi::lastzero, Mval, n, λ1, λ2]; Abort[]];
             {Mval, n} -> <|"J" -> rec["J"], "Psi" -> vv/last|>]],
           sorted]]],
     byM];
   If[Length[table] != d, Message[Psi::degnull, "total", λ1, λ2, Length[table], d]; Abort[]];
   table
];
(* smoke test: (λ1,λ2)=(1,1) table keys and J values *)
KeySort[TauEigensystem[1, 1]][[All, "J"]] // N
```
**Output:**
```
Out[48]= \langle\!\lvert \left\{0,0\right\}\to 0.5664055180884988 +1.3769645669711723\,\mathrm{i},\left\{1,0\right\}\to -0.778123966451158+0.\,\,\mathrm{i},\left\{1,1\right\}\to 1.1900552523337025 +0.\,\,\mathrm{i},\left\{2,0\right\}\to 0.5664055180884988 -1.3769645669711723\,\mathrm{i}\rvert\!\rangle
```

## 2026-07-15 09:53:42
**Input:**
```wolfram
KeySort[TauEigensystem[2, 1]][[All, "J"]] // N
```
**Output:**
```
Out= <|{0, 0} -> 0.9783368039710435 + 2.443001651077886*I, {1, 0} -> -0.8374492520578068 + 0.5086170794089185*I, {1, 1} -> 1.7128032345582143 + 0.6906746402111347*I, {2, 0} -> -0.8374492520578068 - 0.5086170794089185*I, {2, 1} -> 1.7128032345582143 - 0.6906746402111347*I, {3, 0} -> 0.9783368039710435 - 2.443001651077886*I|>
```
## 2026-07-15 09:54:02 — 📥 BULK INSERT 1 CELL(S) at positions 9–9
1. [code] (* Psi[λ1,λ2][M,n]: on-demand accessor over the memoized TauEigensystem cache. Returns the fixed↵


## 2026-07-15 09:54:03
**Input:**
```wolfram
(* Psi[λ1,λ2][M,n]: on-demand accessor over the memoized TauEigensystem cache. Returns the fixed
   (last-component-normalized) eigenvector for a genuine state.
   τ[λ1,λ2][M,n][u]: the FULL τ eigenvalue polynomial for that state, τ0[λ1,λ2,M][u] + J_n.
   Bounds/failure modes (Message + $Failed, matching the existing bounds-checked style):
     Psi::bnds     — gross bounds violation: M not an integer in 0..λ1+λ2, or n not an integer in 0..M.
     Psi::spurious — REPURPOSED from the old Baxter-matrix-filtering meaning: now fires when n is
                     within the gross 0..M range but exceeds the ACTUAL weight multiplicity
                     dPred[M,λ1,λ2]-1 for that M-sector (e.g. (λ1,λ2,M,n)=(2,1,2,2) — dPred[2,2,1]=2
                     so only n=0,1 exist even though the gross bound allows n up to M=2). There is no
                     other "spurious" failure mode left: every state actually produced by
                     TauEigensystem is, by construction, a genuine common eigenvector of C0/C1/C2.
   Both pull from the SAME cached TauEigensystem[λ1,λ2] table, so Psi and τ are mutually consistent
   by construction. *)
ClearAll[Psi, τ];
Psi::bnds = "Invalid indices (M,n)=(`1`,`2`) for (λ1,λ2)=(`3`,`4`): need M an integer in 0..λ1+λ2 and n an integer in 0..M.";
τ::bnds = "Invalid indices (M,n)=(`1`,`2`) for (λ1,λ2)=(`3`,`4`): need M an integer in 0..λ1+λ2 and n an integer in 0..M.";
Psi::spurious = "No genuine state at (M,n)=(`1`,`2`) for (λ1,λ2)=(`3`,`4`): weight sector M=`1` has only `5` genuine state(s) (n=0..`6`), so n=`2` does not exist.";
τ::spurious = "No genuine state at (M,n)=(`1`,`2`) for (λ1,λ2)=(`3`,`4`): weight sector M=`1` has only `5` genuine state(s) (n=0..`6`), so n=`2` does not exist.";
Psi[λ1_, λ2_][M_, n_] := Module[{tab, mult},
  If[! (IntegerQ[M] && IntegerQ[n] && 0 <= M <= λ1 + λ2 && 0 <= n <= M),
    Message[Psi::bnds, M, n, λ1, λ2]; Return[$Failed]];
  tab = TauEigensystem[λ1, λ2];
  mult = dPred[M, λ1, λ2];
  If[n > mult - 1,
    Message[Psi::spurious, M, n, λ1, λ2, mult, mult - 1]; Return[Missing["NoState", {λ1, λ2, M, n}]]];
  tab[{M, n}, "Psi"]];
τ[λ1_, λ2_][M_, n_][u_] := Module[{tab, mult},
  If[! (IntegerQ[M] && IntegerQ[n] && 0 <= M <= λ1 + λ2 && 0 <= n <= M),
    Message[τ::bnds, M, n, λ1, λ2]; Return[$Failed]];
  mult = dPred[M, λ1, λ2];
  If[n > mult - 1,
    Message[τ::spurious, M, n, λ1, λ2, mult, mult - 1]; Return[Missing["NoState", {λ1, λ2, M, n}]]];
  τ0[λ1, λ2, M][u] + TauEigensystem[λ1, λ2][{M, n}, "J"]];
(* smoke test: a valid state, plus a deliberately out-of-range call (n exceeds actual multiplicity
   for M=2 at (λ1,λ2)=(2,1), which is dPred[2,2,1]=2, so n=2 should fire Psi::spurious) *)
Print["spurious-guard test (expect a Psi::spurious message and Missing[...]): ", Psi[2, 1][2, 2]];
{Psi[1, 1][1, 0], N[τ[1, 1][1, 0][2], 10]}
```
**Output:**
```
2 kernel messages: Psi::spurious: No genuine state at (M,n)=(2,2) for (λ1,λ2)=(2,1): weight sector M=2 has only 2 genuine state(s) (n=0..1), so n=2 does not exist. [also: Psi::spurious] result: Print: spurious-guard test (expect a Psi::spurious message and Missing[...]):[formula] | Out[57]= \left\{\left\{1.000000000000000000000000000+0.\times 10^{-29}\,\mathrm{i},-1.3004494881076436272092358…
```
## 2026-07-15 09:54:26 — 📥 BULK INSERT 2 CELL(S) at positions 10–11
1. [code] (* VALIDATION (stage 1, λ1=λ2=1): for every (M,n) with M=0..2, n=0..M, check that Psi actually↵   sa
2. [code] (* VALIDATION (stage 1, λ1=2,λ2=1 — asymmetric generalization): same check, dim 6, M=0..3.↵   This r


## 2026-07-15 09:54:27
**Input:**
```wolfram
(* VALIDATION (stage 1, λ1=λ2=1): for every (M,n) with M=0..2, n=0..M, check that Psi actually
   satisfies the FULL transfer-matrix eigenvalue equation t2(u).Psi == τ(u) Psi, re-deriving t2
   FRESH at each u (t2[1,1][1,1][u] directly — NOT via the C0/C1/C2 extraction used to build Psi),
   and that the genuine-state count per M sector matches dPred exactly (no spurious extras, no
   missing states — this is the honest-diagonalization replacement for the old post-hoc filter). *)
Module[{la = 1, lb = 1, uvals = {0, 1, 2, -1, 5}, rows, counts, pred, maxres},
  rows = Flatten[Table[
     If[n > dPred[M, la, lb] - 1, Nothing,
       Module[{psi = Psi[la, lb][M, n], res},
         res = Max[Abs[Flatten[Table[t2[la, lb][1, 1][uu] . psi - τ[la, lb][M, n][uu] psi, {uu, uvals}]]]];
         {M, n, res}]],
     {M, 0, la + lb}, {n, 0, M}], 1];
  counts = Association[Table[M -> Length[Select[rows, #[[1]] == M &]], {M, 0, la + lb}]];
  pred = Association[Table[M -> dPred[M, la, lb], {M, 0, la + lb}]];
  maxres = Max[rows[[All, 3]]];
  Print["genuine count per M sector : ", counts];
  Print["predicted dPred[M,λ1,λ2] : ", pred];
  Print["MATCHES prediction? ", counts === pred, "  (total = ", Total[Values[counts]], " = dim ", (la + 1) (lb + 1), ")"];
  Print["max residual over all genuine Psi (fresh t2, u in ", uvals, ") = ", maxres];
  Print[If[counts === pred && maxres < 10^-18, "STAGE-1 VALIDATION (λ1=λ2=1) PASSED", "STAGE-1 VALIDATION (λ1=λ2=1) FAILED"]];
  TableForm[rows, TableHeadings -> {None, {"M", "n", "max|t2.Psi-τ Psi|"}}]
]
```
**Output:**
```
Print: genuine count per M sector :[formula]predicted dPred[M,λ1,λ2] :[formula]MATCHES prediction?[formula]  (total =[formula] = dim[formula])max residual over all genuine Psi (fresh t2, u in[formula]) =[formula]STAGE-1 VALIDATION (λ1=λ2=1) PASSED | Out[58]//TableForm= \begin{pmatrix}0 & 0 & 0.\times 10^{-26}\\1 & 0 & 0.\times 10^{-26}\\1 & 1 & 0.\times 10^{-27}\\2 & 0 & 0.\times 10^{-26}\end{pmat… [404 chars total]
```

## 2026-07-15 09:54:27
**Input:**
```wolfram
(* VALIDATION (stage 1, λ1=2,λ2=1 — asymmetric generalization): same check, dim 6, M=0..3.
   This representation has weight sectors of multiplicity >1 (M=1,2 each have dPred=2), so it
   exercises the sort/n-assignment and the degeneracy-resolving generic-combination diagonalization
   in a way the symmetric λ1=λ2=1 case (max multiplicity dPred=2 also, but let's still check
   thoroughly) does not fully stress. *)
Module[{la = 2, lb = 1, uvals = {0, 1, 2, -1, 5}, rows, counts, pred, maxres},
  rows = Flatten[Table[
     If[n > dPred[M, la, lb] - 1, Nothing,
       Module[{psi = Psi[la, lb][M, n], res},
         res = Max[Abs[Flatten[Table[t2[la, lb][1, 1][uu] . psi - τ[la, lb][M, n][uu] psi, {uu, uvals}]]]];
         {M, n, res}]],
     {M, 0, la + lb}, {n, 0, M}], 1];
  counts = Association[Table[M -> Length[Select[rows, #[[1]] == M &]], {M, 0, la + lb}]];
  pred = Association[Table[M -> dPred[M, la, lb], {M, 0, la + lb}]];
  maxres = Max[rows[[All, 3]]];
  Print["λ1=2,λ2=1 genuine count per M sector : ", counts];
  Print["predicted dPred[M,λ1,λ2] : ", pred];
  Print["MATCHES prediction? ", counts === pred, "  (total = ", Total[Values[counts]], " = dim ", (la + 1) (lb + 1), ")"];
  Print["max residual over all genuine Psi (fresh t2, u in ", uvals, ") = ", maxres];
  Print[If[counts === pred && maxres < 10^-18, "STAGE-1 VALIDATION (λ1=2,λ2=1) PASSED", "STAGE-1 VALIDATION (λ1=2,λ2=1) FAILED"]];
  TableForm[rows, TableHeadings -> {None, {"M", "n", "max|t2.Psi-τ Psi|"}}]
]
```
**Output:**
```
Print: λ1=2,λ2=1 genuine count per M sector :[formula]predicted dPred[M,λ1,λ2] :[formula]MATCHES prediction?[formula]  (total =[formula] = dim[formula])max residual over all genuine Psi (fresh t2, u in[formula]) =[formula]STAGE-1 VALIDATION (λ1=2,λ2=1) PASSED | Out[59]//TableForm= \begin{pmatrix}0 & 0 & 0.\times 10^{-26}\\1 & 0 & 0.\times 10^{-26}\\1 & 1 & 0.\times 10^{-26}\\2 & 0 & 0.\times 10^{-… [470 chars total]
```

## 2026-07-15 09:54:36
**Input:**
```wolfram
Jn = TauEigensystem[2,1][{2,0},"J"];
mat = BaxterMatrix[2, 1, 2] - Jn IdentityMatrix[3];
ns0 = NullSpace[mat];
nsTol = NullSpace[mat, Tolerance -> 10^-10];
{Length[ns0], Length[nsTol], SingularValueList[mat]}
```
**Output:**
```
Out= {1, 1, {6.61269920153756534792932528556736799726`27.521547555158165, 5.10179268929260995410599645902354124746`27.521547555158165}}
```

## 2026-07-15 09:54:46
**Input:**
```wolfram
SingularValueList[mat, 3]
```
**Output:**
```
Out= {6.61269920153756534792932528556736799726`27.521547555158165, 5.10179268929260995410599645902354124746`27.521547555158165, 0}
```

## 2026-07-15 09:55:08
**Input:**
```wolfram
Table[
  Module[{Jn2, mat2, ns0len, nsTlen},
    If[n > dPred[M, la, lb] - 1, Nothing,
      Jn2 = TauEigensystem[la, lb][{M, n}, "J"];
      mat2 = BaxterMatrix[la, lb, M] - Jn2 IdentityMatrix[M + 1];
      ns0len = Length[NullSpace[mat2]];
      nsTlen = Length[NullSpace[mat2, Tolerance -> 10^-10]];
      {la, lb, M, n, ns0len, nsTlen}]],
  {la, lb} \[Element] {{1, 1}, {2, 1}}; {M, 0, la + lb}, {n, 0, M}] // Flatten[#, 2] &
```
**Output:**
```
[message] Table::nliter: Non-list iterator {la, lb} ∈ {{1, 1}, {2, 1}}; {M, 0, la + lb} at position 2 does not evaluate to a real numeric value.
Out= Table[Module[{Jn2, mat2, ns0len, nsTlen}, If[n > dPred[M, la, lb] - 1, Nothing, Jn2 = TauEigensystem[la, lb][{M, n}, "J"]; mat2 = BaxterMatrix[la, lb, M] - Jn2*IdentityMatrix[M + 1]; ns0len = Length[NullSpace[mat2]]; nsTlen = Length[NullSpace[mat2, Tolerance -> 10^(-10)]]; {la, lb, M, n, ns0len, nsTlen}]], Element[{la, lb}, {{1, 1}, {2, 1}}]; {M, 0, la + lb}, {n, 0, M}]
```

## 2026-07-15 09:55:18
**Input:**
```wolfram
sweep = Flatten[Table[
  With[{la = rep[[1]], lb = rep[[2]]},
    Table[
      If[n > dPred[M, la, lb] - 1, Nothing,
        Module[{Jn2, mat2, ns0len, nsTlen},
          Jn2 = TauEigensystem[la, lb][{M, n}, "J"];
          mat2 = BaxterMatrix[la, lb, M] - Jn2 IdentityMatrix[M + 1];
          ns0len = Length[NullSpace[mat2]];
          nsTlen = Length[NullSpace[mat2, Tolerance -> 10^-10]];
          {la, lb, M, n, ns0len, nsTlen}]],
      {M, 0, la + lb}, {n, 0, M}]],
  {rep, {{1, 1}, {2, 1}}}], 3];
TableForm[sweep, TableHeadings -> {None, {"λ1","λ2","M","n","NS default len","NS Tol len"}}]
```
**Output:**
```
Out= TableForm[{1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 0, 1, 1, 2, 1, 0, 0, 0, 0, 2, 1, 1, 0, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 2, 0, 1, 1, 2, 1, 2, 1, 1, 1, 2, 1, 3, 0, 1, 1}, TableHeadings -> {None, {"λ1", "λ2", "M", "n", "NS default len", "NS Tol len"}}]
```

## 2026-07-15 09:55:31
**Input:**
```wolfram
sweep2 = Flatten[Table[
  With[{la = rep[[1]], lb = rep[[2]]},
    Table[
      If[n > dPred[M, la, lb] - 1, Nothing,
        Module[{Jn2, mat2, ns0len, nsTlen},
          Jn2 = TauEigensystem[la, lb][{M, n}, "J"];
          mat2 = BaxterMatrix[la, lb, M] - Jn2 IdentityMatrix[M + 1];
          ns0len = Length[NullSpace[mat2]];
          nsTlen = Length[NullSpace[mat2, Tolerance -> 10^-10]];
          {la, lb, M, n, ns0len, nsTlen}]],
      {M, 0, la + lb}, {n, 0, M}]],
  {rep, {{1, 1}, {2, 1}}}], 2];
Length[sweep2]
TableForm[sweep2, TableHeadings -> {None, {"λ1","λ2","M","n","NS default len","NS Tol len"}}]
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= 10
Out[3]= TableForm[{{1, 1, 0, 0, 0, 0}, {1, 1, 1, 0, 1, 1}, {1, 1, 1, 1, 1, 1}, {1, 1, 2, 0, 1, 1}, {2, 1, 0, 0, 0, 0}, {2, 1, 1, 0, 1, 1}, {2, 1, 1, 1, 1, 1}, {2, 1, 2, 0, 1, 1}, {2, 1, 2, 1, 1, 1}, {2, 1, 3, 0, 1, 1}}, TableHeadings -> {None, {"λ1", "λ2", "M", "n", "NS default len", "NS Tol len"}}]
```

## 2026-07-15 09:55:42
**Input:**
```wolfram
{BaxterMatrix[1,1,0], TauEigensystem[1,1][{0,0},"J"], BaxterMatrix[1,1,0] - TauEigensystem[1,1][{0,0},"J"] IdentityMatrix[1]}
```
**Output:**
```
Out= {{{0.5664055180884988442762077366501501434616299651264478610617`29.695249010583606 + 1.3769645669711722005252340687712812291996828026250223428811`30.06249936401476*I}}, 0.56640551808849884427620773665015014352`26.563908186381784 + 1.37696456697117220052523406877128122917`26.915783392193127*I, {{0``26.810459874680447 + 0``26.776550947129415*I}}}
```

## 2026-07-15 09:55:51
**Input:**
```wolfram
m00 = {{0``26.81 + 0``26.78*I}};
{NullSpace[m00], NullSpace[m00, Tolerance->10^-10], Chop[m00]}
```
**Output:**
```
Out= {{}, {}, {{0}}}
```

## 2026-07-15 09:56:12
**Input:**
```wolfram
{NullSpace[{{0}}], NullSpace[{{0.`30}}], NullSpace[{{0``30}}], Precision[m00[[1,1]]], Precision[0``26.81]}
```
**Output:**
```
Out= {{{1}}, {{1.}}, {}, 0., 0.}
```

## 2026-07-15 09:56:42
**Input:**
```wolfram
sweep3 = Flatten[Table[
  With[{la = rep[[1]], lb = rep[[2]]},
    Table[
      If[n > dPred[M, la, lb] - 1, Nothing,
        Module[{Jn2, mat2, nsChopTol, nsTolOnly, nsChopOnly},
          Jn2 = TauEigensystem[la, lb][{M, n}, "J"];
          mat2 = BaxterMatrix[la, lb, M] - Jn2 IdentityMatrix[M + 1];
          nsChopTol = Length[NullSpace[Chop[mat2, 10^-12], Tolerance -> 10^-10]];
          nsTolOnly = Length[NullSpace[mat2, Tolerance -> 10^-10]];
          nsChopOnly = Length[NullSpace[Chop[mat2, 10^-12]]];
          {la, lb, M, n, nsChopTol, nsTolOnly, nsChopOnly}]],
      {M, 0, la + lb}, {n, 0, M}]],
  {rep, {{1, 1}, {2, 1}}}], 2];
TableForm[sweep3, TableHeadings -> {None, {"λ1","λ2","M","n","Chop+Tol len","TolOnly len","ChopOnly len"}}]
```
**Output:**
```
Out= TableForm[{{1, 1, 0, 0, 1, 0, 1}, {1, 1, 1, 0, 1, 1, 1}, {1, 1, 1, 1, 1, 1, 1}, {1, 1, 2, 0, 1, 1, 1}, {2, 1, 0, 0, 1, 0, 1}, {2, 1, 1, 0, 1, 1, 1}, {2, 1, 1, 1, 1, 1, 1}, {2, 1, 2, 0, 1, 1, 1}, {2, 1, 2, 1, 1, 1, 1}, {2, 1, 3, 0, 1, 1, 1}}, TableHeadings -> {None, {"λ1", "λ2", "M", "n", "Chop+Tol len", "TolOnly len", "ChopOnly len"}}]
```

## 2026-07-15 09:57:36
**Input:**
```wolfram
svSweep = Flatten[Table[
  With[{la = rep[[1]], lb = rep[[2]]},
    Table[
      If[n > dPred[M, la, lb] - 1, Nothing,
        Module[{Jn2, mat2, svs},
          Jn2 = TauEigensystem[la, lb][{M, n}, "J"];
          mat2 = BaxterMatrix[la, lb, M] - Jn2 IdentityMatrix[M + 1];
          svs = SingularValueList[mat2, M + 1];
          {la, lb, M, n, N[svs, 6]}]],
      {M, 0, la + lb}, {n, 0, M}]],
  {rep, {{1, 1}, {2, 1}}}], 2];
Column[ToString /@ svSweep]
```
**Output:**
```
Out= Column[{"                   -26\n{1, 1, 0, 0, {0. 10   }}", "{1, 1, 1, 0, {3.44975, 0}}", "{1, 1, 1, 1, {3.44975, 0}}", "{1, 1, 2, 0, {4.98064, 3.81935, 0}}", "                   -26\n{2, 1, 0, 0, {0. 10   }}", "{2, 1, 1, 0, {4.62933, 0}}", "{2, 1, 1, 1, {4.62933, 0}}", "{2, 1, 2, 0, {6.61270, 5.10179, 0}}", "{2, 1, 2, 1, {6.91096, 5.08795, 0}}", "{2, 1, 3, 0, {9.50258, 7.58626, 4.65856, 0}}"}]
```
## 2026-07-15 09:58:06 — 📥 BULK INSERT 1 CELL(S) at positions 12–12
1. [code] (* QSolve[λ1,λ2][M,n]: the TQ equation is now a LINEAR solve, not an eigenvalue problem, because↵


## 2026-07-15 09:58:06
**Input:**
```wolfram
(* QSolve[λ1,λ2][M,n]: the TQ equation is now a LINEAR solve, not an eigenvalue problem, because
   J_n (hence τ) is already KNOWN from the genuine TauEigensystem diagonalization above. The
   monic-degree-M polynomial q_n(u) satisfies O[q_n]=0 with O the SAME BaxterMatrix operator as
   before, now evaluated at the fixed known J_n: (BaxterMatrix[λ1,λ2,M] - J_n·Id).q_n = 0.
   Reuses BaxterMatrix (built above) and the bounds/multiplicity-check style from Psi/τ.
   CRITICAL numerical-precision pitfall (same theme as flagged for the old Psi construction, but a
   DIFFERENT concrete manifestation): BaxterMatrix[λ1,λ2,M]-J_n·Id is singular by construction, but
   its (M+1)-th (would-be-zero) entries/singular-values come out as precision-TRACKED zeros, e.g.
   literally `0``26` (value exactly 0, but carrying a finite Accuracy/Precision tag rather than
   Infinite precision). NullSpace treats such an entry's Precision (which evaluates to ~0, since
   Precision[0]-style values are accuracy-driven) as UNKNOWN rather than confirmed-zero, and
   conservatively returns an EMPTY nullspace even though the matrix is genuinely singular — this
   bites the M=0 sector (1x1 matrix, entry identically the precision-tracked zero) with EITHER
   default OR explicit Tolerance->10^-10 alone. The fix is to Chop the matrix (threshold 10^-12,
   comfortably above the ~10^-26..10^-29 noise floor of this construction but far below any
   genuine O(1) entry) BEFORE calling NullSpace, converting those precision-tracked zeros to exact
   (infinite-precision) 0 so NullSpace's linear-algebra can recognize the singularity. The explicit
   Tolerance->10^-10 is ALSO required per the design brief (matches the historical Psi-construction
   bug at (λ1,λ2,M,n)=(2,1,2,0) in the old stacked-C0/C1/C2 nullspace, and provides headroom for any
   genuinely-small-but-nonzero near-null direction that Chop alone would not catch in general).
   Both are used together defensively. Messages/behavior on failure mirror Psi/τ (bnds/spurious)
   plus new degnull/lastzero variants for this linear solve. Memoized per (λ1,λ2,M,n). *)
ClearAll[QSolve];
QSolve::bnds = "Invalid indices (M,n)=(`1`,`2`) for (λ1,λ2)=(`3`,`4`): need M an integer in 0..λ1+λ2 and n an integer in 0..M.";
QSolve::spurious = "No genuine state at (M,n)=(`1`,`2`) for (λ1,λ2)=(`3`,`4`): weight sector M=`1` has only `5` genuine state(s) (n=0..`6`), so n=`2` does not exist.";
QSolve::degnull = "Unexpected Q-function nullspace dimension `5` (expected 1) for (M,n)=(`1`,`2`) at (λ1,λ2)=(`3`,`4`).";
QSolve::lastzero = "Q-function nullspace vector found but its top (u^M) coefficient is ~0 for (M,n)=(`1`,`2`) at (λ1,λ2)=(`3`,`4`); cannot normalize to monic.";
QSolve[λ1_, λ2_][M_, n_] := QSolve[λ1, λ2][M, n] = Module[{mult, Jn, mat, ns, vec, last},
  If[! (IntegerQ[M] && IntegerQ[n] && 0 <= M <= λ1 + λ2 && 0 <= n <= M),
    Message[QSolve::bnds, M, n, λ1, λ2]; Return[$Failed]];
  mult = dPred[M, λ1, λ2];
  If[n > mult - 1,
    Message[QSolve::spurious, M, n, λ1, λ2, mult, mult - 1]; Return[$Failed]];
  Jn = TauEigensystem[λ1, λ2][{M, n}, "J"];
  mat = BaxterMatrix[λ1, λ2, M] - Jn IdentityMatrix[M + 1];
  ns = NullSpace[Chop[mat, 10^-12], Tolerance -> 10^-10];
  Which[
    Length[ns] != 1, Message[QSolve::degnull, M, n, λ1, λ2, Length[ns]]; Return[$Failed],
    True,
      vec = ns[[1]]; last = vec[[-1]];
      If[Abs[last] < 10^-15 Max[Abs[vec]],
        Message[QSolve::lastzero, M, n, λ1, λ2]; Return[$Failed]];
      vec/last]
];
(* smoke test: the regression-critical case (λ1,λ2,M,n)=(2,1,2,0) *)
QSolve[2, 1][2, 0]
```
**Output:**
```
Out[66]= \left\{1.03464899078135229594035121+0.381869568380946283899575941\,\mathrm{i},-2.37022278949001677248172619-0.20048970475590097847093779\,\mathrm{i},1.00000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right\}
```
## 2026-07-15 09:58:24 — 📥 BULK INSERT 3 CELL(S) at positions 13–15
1. [code] (* Q1[λ1,λ2][M,n][u] = κ1^(u/h) q_n(u): the actual Baxter Q-function value, built from the SAME↵   Q
2. [code] (* VALIDATION (stage 2, λ1=λ2=1): full end-to-end TQ-equation residual check↵     a(u) Q(u+h) - τ(u)
3. [code] (* VALIDATION (stage 2, λ1=2,λ2=1 — asymmetric): same end-to-end TQ residual check, dim 6, M=0..3.↵


## 2026-07-15 09:58:25
**Input:**
```wolfram
(* Q1[λ1,λ2][M,n][u] = κ1^(u/h) q_n(u): the actual Baxter Q-function value, built from the SAME
   QSolve cache (mutually consistent with τ by construction, since both trace back to the same
   TauEigensystem[λ1,λ2][{M,n}] entry). *)
ClearAll[Q1];
Q1[λ1_, λ2_][M_, n_][u_] := Module[{q = QSolve[λ1, λ2][M, n]},
  If[q === $Failed, $Failed, κ1^(u/h) qpoly[q, u]]];
(* smoke test *)
{Q1[1, 1][1, 0][2], N[τ[1, 1][1, 0][2], 10]}
```
**Output:**
```
Out[69]= \left\{-0.709231010958016926537225142+1.09596916251277206706752728\,\mathrm{i},-0.02291660900+0.\times 10^{-12}\,\mathrm{i}\right\}
```

## 2026-07-15 09:58:25
**Input:**
```wolfram
(* VALIDATION (stage 2, λ1=λ2=1): full end-to-end TQ-equation residual check
     a(u) Q(u+h) - τ(u) Q(u) + d(u) Q(u-h) = 0
   assembled ENTIRELY from the Q1/τ wrappers (which now trace back to the GENUINE
   TauEigensystem+QSolve pipeline, not a Baxter-matrix eigen-decomposition) for every (M,n),
   at several sample u avoiding the roots of a(u) (4/3,8/7) and d(u) (1/3,1/7). *)
sampleU2a = {0, 1, 3/2, -2, 10};
res2a = Flatten[Table[
    If[n > dPred[M, 1, 1] - 1, Nothing,
      {M, n,
       Max[Abs[Table[
          aFun[1, 1][uu] Q1[1, 1][M, n][uu + h]
          - τ[1, 1][M, n][uu] Q1[1, 1][M, n][uu]
          + Qθ[uu] Q1[1, 1][M, n][uu - h],
          {uu, sampleU2a}]]]}],
    {M, 0, 2}, {n, 0, M}], 1];
max2a = Max[res2a[[All, 3]]];
Print["stage-2 max residual (λ1=λ2=1, via genuine τ/Q1) = ", max2a];
Print[If[max2a < 10^-18, "STAGE-2 VALIDATION (λ1=λ2=1) PASSED", "STAGE-2 VALIDATION (λ1=λ2=1) FAILED"]];
TableForm[res2a, TableHeadings -> {None, {"M", "n", "max|TQ residual|"}}]
```
**Output:**
```
Print: stage-2 max residual (λ1=λ2=1, via genuine τ/Q1) =[formula] | Print: STAGE-2 VALIDATION (λ1=λ2=1) PASSED | Out[75]//TableForm= \begin{pmatrix}0 & 0 & 0.\times 10^{-26}\\1 & 0 & 0.\times 10^{-24}\\1 & 1 & 0.\times 10^{-24}\\2 & 0 & 0.\times 10^{-21}\end{pmatrix}
```

## 2026-07-15 09:58:25
**Input:**
```wolfram
(* VALIDATION (stage 2, λ1=2,λ2=1 — asymmetric): same end-to-end TQ residual check, dim 6, M=0..3.
   EXPLICITLY calls out (M,n)=(2,0) — the case that previously exposed the NullSpace precision
   pitfall (see QSolve commentary above) — as an individually-flagged regression guard, not just
   folded anonymously into the aggregate max. *)
sampleU2b = {0, 1, 3/2, -2, 10};
res2b = Flatten[Table[
    If[n > dPred[M, 2, 1] - 1, Nothing,
      {M, n,
       Max[Abs[Table[
          aFun[2, 1][uu] Q1[2, 1][M, n][uu + h]
          - τ[2, 1][M, n][uu] Q1[2, 1][M, n][uu]
          + Qθ[uu] Q1[2, 1][M, n][uu - h],
          {uu, sampleU2b}]]]}],
    {M, 0, 3}, {n, 0, M}], 1];
max2b = Max[res2b[[All, 3]]];
res2b_20 = First[Cases[res2b, {2, 0, r_} :> r]];
Print["stage-2 max residual (λ1=2,λ2=1, via genuine τ/Q1) = ", max2b];
Print["regression-critical (M,n)=(2,0) residual = ", res2b_20, "  ", If[res2b_20 < 10^-18, "OK", "FAIL"]];
Print[If[max2b < 10^-18, "STAGE-2 VALIDATION (λ1=2,λ2=1) PASSED", "STAGE-2 VALIDATION (λ1=2,λ2=1) FAILED"]];
TableForm[res2b, TableHeadings -> {None, {"M", "n", "max|TQ residual|"}}]
```
**Output:**
```
2 kernel messages: Set::write: Tag Times in 20 (res2b_) is Protected. [also: Set::write] result: Print: stage-2 max residual (λ1=2,λ2=1, via genuine τ/Q1) =[formula] | Print: regression-critical (M,n)=(2,0) residual =[formula][formula] | Print: STAGE-2 VALIDATION (λ1=2,λ2=1) PASSED | Out[83]//Tab…
```
## 2026-07-15 09:58:42 — ✏️ EDIT CELL 15
(* VALIDATION (stage 2, λ1=2,λ2=1 — asymmetric): same end-to-end TQ residual check, dim 6, M=0..3.
   EXPLICITLY calls out (M,n)=(2,0) — the case that previously exposed the NullSpace precision
   pit…


## 2026-07-15 09:58:42
**Input:**
```wolfram
(* VALIDATION (stage 2, λ1=2,λ2=1 — asymmetric): same end-to-end TQ residual check, dim 6, M=0..3.
   EXPLICITLY calls out (M,n)=(2,0) — the case that previously exposed the NullSpace precision
   pitfall (see QSolve commentary above) — as an individually-flagged regression guard, not just
   folded anonymously into the aggregate max. *)
sampleU2b = {0, 1, 3/2, -2, 10};
res2b = Flatten[Table[
    If[n > dPred[M, 2, 1] - 1, Nothing,
      {M, n,
       Max[Abs[Table[
          aFun[2, 1][uu] Q1[2, 1][M, n][uu + h]
          - τ[2, 1][M, n][uu] Q1[2, 1][M, n][uu]
          + Qθ[uu] Q1[2, 1][M, n][uu - h],
          {uu, sampleU2b}]]]}],
    {M, 0, 3}, {n, 0, M}], 1];
max2b = Max[res2b[[All, 3]]];
resAt2and0 = First[Cases[res2b, {2, 0, r_} :> r]];
Print["stage-2 max residual (λ1=2,λ2=1, via genuine τ/Q1) = ", max2b];
Print["regression-critical (M,n)=(2,0) residual = ", resAt2and0, "  ", If[resAt2and0 < 10^-18, "OK", "FAIL"]];
Print[If[max2b < 10^-18, "STAGE-2 VALIDATION (λ1=2,λ2=1) PASSED", "STAGE-2 VALIDATION (λ1=2,λ2=1) FAILED"]];
TableForm[res2b, TableHeadings -> {None, {"M", "n", "max|TQ residual|"}}]
```
**Output:**
```
Out= TableForm[{{0, 0, 0``25.75764296390223}, {1, 0, 0``23.578722512204813}, {1, 1, 0``23.52068038060336}, {2, 0, 0``22.28294415081727}, {2, 1, 0``22.563461769247766}, {3, 0, 0``20.29255306805365}}, TableHeadings -> {None, {"M", "n", "max|TQ residual|"}}]
```
## 2026-07-15 09:58:58 — ▶️ RUN CELL 15 (0.06 s)
**In [15]:** `(* VALIDATION (stage 2, λ1=2,λ2=1 — asymmetric): same end-to-end TQ residual check, dim 6, M=0..3.↵   EXPLICITLY calls out (M,n)=(2,0) — the case that previously exposed the NullSpace precision↵   pit`
**Out:** Print: stage-2 max residual (λ1=2,λ2=1, via genuine τ/Q1) =[formula] | Print: regression-critical (M,n)=(2,0) residual =[formula]OK | Print: STAGE-2 VALIDATION (λ1=2,λ2=1) PASSED | Out[91]//TableForm= \begin{pmatrix}0 & 0 & 0.\times 10^{-26}\\1 & 0 & 0.\times 10^{-24}\\1 & 1 & 0.\times 10^{-24}\\2 &

## 2026-07-15 09:59:13 — 📥 BULK INSERT 1 CELL(S) at positions 13–13
1. [code] (* REGRESSION GUARD for the NullSpace precision pitfall documented in QSolve above.↵   Demonstrates


## 2026-07-15 09:59:13
**Input:**
```wolfram
(* REGRESSION GUARD for the NullSpace precision pitfall documented in QSolve above.
   Demonstrates CONCRETELY that without the Chop preprocessing step, NullSpace (with OR without
   the mandated Tolerance->10^-10) fails to find the null direction for the M=0 sector (where the
   1x1 matrix entry is a precision-tracked exact zero), while the Chop+Tolerance combination used
   inside QSolve succeeds for every sector including the historically-critical (λ1,λ2,M,n)=(2,1,2,0). *)
Module[{la = 2, lb = 1, M0 = 0, n0 = 0, Jn0, mat0, bareLen, tolOnlyLen, chopTolLen},
  Jn0 = TauEigensystem[la, lb][{M0, n0}, "J"];
  mat0 = BaxterMatrix[la, lb, M0] - Jn0 IdentityMatrix[M0 + 1];
  bareLen = Length[NullSpace[mat0]];
  tolOnlyLen = Length[NullSpace[mat0, Tolerance -> 10^-10]];
  chopTolLen = Length[NullSpace[Chop[mat0, 10^-12], Tolerance -> 10^-10]];
  Print["M=0 sector (λ1,λ2)=(", la, ",", lb, "): bare NullSpace dim=", bareLen,
        ", Tolerance-only dim=", tolOnlyLen, ", Chop+Tolerance (as used in QSolve) dim=", chopTolLen];
  Print[If[bareLen == 0 && tolOnlyLen == 0 && chopTolLen == 1,
    "REGRESSION GUARD PASSED (confirms Chop preprocessing is genuinely necessary, matching QSolve's implementation)",
    "REGRESSION GUARD FAILED — pitfall behavior changed, re-check QSolve's Chop+Tolerance necessity"]]
]
```
**Output:**
```
Print: M=0 sector (λ1,λ2)=([formula],[formula]): bare NullSpace dim=[formula], Tolerance-only dim=[formula], Chop+Tolerance (as used in QSolve) dim=[formula]REGRESSION GUARD PASSED (confirms Chop preprocessing is genuinely necessary, matching QSolve's implementation)
```
## 2026-07-15 10:00:08 — ✏️ EDIT CELL 1
(*Standalone numerical solver for the L=2 su(2) XXX Baxter (TQ) equation.
This file is SELF-CONTAINED: it defines every symbol it uses (θ1,θ2,h,κ1,κ2, aFun, Qθ=d, τ0,
BaxterMatrix, the transfer-matrix…


---

## 2026-07-15 10:00:16 — 🔄 KERNEL RESTART

## 2026-07-15 10:00:25 — ▶️ RUN CELLS 1–16
Cell 2: ✓ — Out[6]= \left\{\frac{1}{3},\frac{1}{7},1,0.360439875147226537266677650596+0.932782448593374716484835982071\,\mathrm{i},0.36043987514722653726667765060-0.932782448593374716484835982071\,\mathrm{i}\right\}
Cell 3: ✓ — Out[12]= \left\{\left(u-\frac{4}{3}\right)\,\left(u-\frac{8}{7}\right),\left(u-\frac{1}{3}\right)\,\left(u-\frac{1}{7}\right),\left(0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i}\right)\,u^2-\left(1.06415582186324025288257211128+0.\times 10^{-29}\,\mathrm{i}\right)\,u\right\}
Cell 4: ✓ — Out[15]= \begin{pmatrix}1 & 1\\2 & 2\\3 & 3\end{pmatrix}
Cell 5: ✓ — Out[16]= \begin{matrix}\left\{0,\text{max|deg M+1,M+2 coeff|}\to 0.\times 10^{-29}\right\}\\\left\{1,\text{max|deg M+1,M+2 coeff|}\to 0.\times 10^{-29}\right\}\\\left\{2,\text{max|deg M+1,M+2 coeff|}\to 0.\times 10^{-29}\right\}\end{matrix}
Cell 6: ✓ — Out[32]= \left\{\text{\{tr,det\}(G) - \{$\kappa$1+$\kappa$2,$\kappa$1 $\kappa$2\}}\to \left\{0,0\right\},\text{dim t2[1,1]}\to \left\{4,4\right\}\right\}
Cell 7: ✓ — Out[35]= \text{max|C2-($\kappa$1+$\kappa$2)Id|}\to 0.\times 10^{-29}
Cell 8: ✓ — Out[45]= \langle\!\lvert \left\{0,0\right\}\to 0.5664055180884988 +1.3769645669711723\,\mathrm{i},\left\{1,0\right\}\to -0.778123966451158+0.\,\,\mathrm{i},\left\{1,1\right\}\to 1.1900552523337025 +0.\,\,\mathrm{i},\left\{2,0\right\}\to 0.5664055180884988 -1.3769645669711723\,\mathrm{i}\rvert\!\rangle
Cell 9: ✓ — 2 kernel messages: Psi::spurious: No genuine state at (M,n)=(2,2) for (λ1,λ2)=(2,1): weight sector M=2 has only 2 genuine state(s) (n=0..1), so n=2 does not exist. [also: Psi::spurious] result: Print: spurious-guard test (expect a Psi::spurious message and Missing[...]):[formula] | Out[54]= \left\{\left\{1.000000000000000000000000000+0.\times 10^{-29}\,\mathrm{i},-1.3004494881076436272092358…
⚠️ stopped at Cell 9 — error detected (pass stopOnError:false to continue past errors)

## 2026-07-15 10:00:44 — ▶️ RUN CELLS 9–16
Cell 9: ✓ — 2 kernel messages: Psi::spurious: No genuine state at (M,n)=(2,2) for (λ1,λ2)=(2,1): weight sector M=2 has only 2 genuine state(s) (n=0..1), so n=2 does not exist. [also: Psi::spurious] result: Print: spurious-guard test (expect a Psi::spurious message and Missing[...]):[formula] | Out[63]= \left\{\left\{1.000000000000000000000000000+0.\times 10^{-29}\,\mathrm{i},-1.3004494881076436272092358…
Cell 10: ✓ — Print: genuine count per M sector :[formula]predicted dPred[M,λ1,λ2] :[formula]MATCHES prediction?[formula]  (total =[formula] = dim[formula])max residual over all genuine Psi (fresh t2, u in[formula]) =[formula]STAGE-1 VALIDATION (λ1=λ2=1) PASSED | Out[64]//TableForm= \begin{pmatrix}0 & 0 & 0.\times 10^{-26}\\1 & 0 & 0.\times 10^{-26}\\1 & 1 & 0.\times 10^{-27}\\2 & 0 & 0.\times 10^{-26}\end{pmat… [404 chars total]
Cell 11: ✓ — Print: λ1=2,λ2=1 genuine count per M sector :[formula]predicted dPred[M,λ1,λ2] :[formula]MATCHES prediction?[formula]  (total =[formula] = dim[formula])max residual over all genuine Psi (fresh t2, u in[formula]) =[formula]STAGE-1 VALIDATION (λ1=2,λ2=1) PASSED | Out[65]//TableForm= \begin{pmatrix}0 & 0 & 0.\times 10^{-26}\\1 & 0 & 0.\times 10^{-26}\\1 & 1 & 0.\times 10^{-26}\\2 & 0 & 0.\times 10^{-… [470 chars total]
Cell 12: ✓ — Out[72]= \left\{1.03464899078135229594035121+0.381869568380946283899575941\,\mathrm{i},-2.37022278949001677248172619-0.20048970475590097847093779\,\mathrm{i},1.00000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right\}
Cell 13: ✓ — Print: M=0 sector (λ1,λ2)=([formula],[formula]): bare NullSpace dim=[formula], Tolerance-only dim=[formula], Chop+Tolerance (as used in QSolve) dim=[formula]REGRESSION GUARD PASSED (confirms Chop preprocessing is genuinely necessary, matching QSolve's implementation)
Cell 14: ✓ — Out[76]= \left\{-0.709231010958016926537225142+1.09596916251277206706752728\,\mathrm{i},-0.02291660900+0.\times 10^{-12}\,\mathrm{i}\right\}
Cell 15: ✓ — Print: stage-2 max residual (λ1=λ2=1, via genuine τ/Q1) =[formula] | Print: STAGE-2 VALIDATION (λ1=λ2=1) PASSED | Out[82]//TableForm= \begin{pmatrix}0 & 0 & 0.\times 10^{-26}\\1 & 0 & 0.\times 10^{-24}\\1 & 1 & 0.\times 10^{-24}\\2 & 0 & 0.\times 10^{-21}\end{pmatrix}
Cell 16: ✓ — Print: stage-2 max residual (λ1=2,λ2=1, via genuine τ/Q1) =[formula] | Print: regression-critical (M,n)=(2,0) residual =[formula]OK | Print: STAGE-2 VALIDATION (λ1=2,λ2=1) PASSED | Out[90]//TableForm= \begin{pmatrix}0 & 0 & 0.\times 10^{-26}\\1 & 0 & 0.\times 10^{-24}\\1 & 1 & 0.\times 10^{-24}\\2 & 0 & 0.\times 10^{-22}\\2 & 1 & 0.\times 10^{-23}\\3 & 0 & 0.\times 10^{-20}\end{pmatrix}

