## Deleted Cell 18 — 2026-07-15 09:47:12
**Kind:** code
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

## Deleted Cell 17 — 2026-07-15 09:47:12
**Kind:** code
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

## Deleted Cell 16 — 2026-07-15 09:47:12
**Kind:** code
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

## Deleted Cell 15 — 2026-07-15 09:47:12
**Kind:** code
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

## Deleted Cell 14 — 2026-07-15 09:47:12
**Kind:** code
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

## Deleted Cell 13 — 2026-07-15 09:47:12
**Kind:** markdown
```markdown
(*Genuine transfer-matrix eigenvectors from the Baxter τ eigenvalue.

Goal: build Psi[λ1,λ2][M,n], a FIXED vector (no [u] argument) in V_λ1⊗V_λ2 of dimension
(λ1+1)(λ2+1), satisfying  t2[λ1,λ2][1,1][u] . Psi == τ[λ1,λ2][M,n][u] * Psi  for ALL u,
where τ is the eigenvalue already produced by the Baxter solver above.

Twist: the companion matrix G={{χ1,-χ2},{1,0}} with χ1=κ1+κ2, χ2=1 has characteristic polynomial
λ²-χ1 λ+χ2 = (λ-κ1)(λ-κ2), i.e. eigenvalues EXACTLY κ1,κ2 (since χ2=κ1 κ2=1 here). This matches
the diagonal twist implicit in the Baxter a(u)/d(u)/τ(u), so the τ eigenvalues are genuine t2
eigenvalues in the twisted transfer matrix.

Method (coefficient-matrix / simultaneous-nullspace): t2[λ1,λ2][1,1][u] = C0 + C1 u + C2 u² with
numeric d×d coefficient matrices (T=L1·L2, each linear in u ⇒ degree ≤2). τ[λ1,λ2][M,n][u] =
t0 + t1 u + t2c u² with scalar coefficients. The eigenvalue equation for all u is equivalent to
the three matrix equations Ck.Psi = tk Psi, so Psi lies in the simultaneous nullspace of
(Ck - tk Id), k=0,1,2. We stack these into a 3d×d matrix and take NullSpace. (Note C2 = tr(G) Id
= (κ1+κ2) Id = t2c Id exactly, so the u² block is ~0 and the real constraints are C0,C1.)

Spurious solutions: the Baxter matrix yields M+1 candidate (M,n) for each M=0..λ1+λ2, but the true
weight-space multiplicity for magnon number M is min(M,λ1,λ2,λ1+λ2-M)+1, which is < M+1 once
M>min(λ1,λ2). For those extra candidates the stacked nullspace is trivial (no genuine eigenvector)
— flagged via Psi::spurious. We attempt ALL (M,n) and report genuine vs spurious empirically.

Normalization: Psi is scaled so its LAST component (in the KroneckerProduct basis ordering used by
L1/L2/t2) is exactly 1; a ~0 last component is flagged (Psi::lastzero) rather than silently giving
ComplexInfinity.

The t2 machinery below (ps,dim,J0,Jp,Jm,JJ,Ee,L,G,L1,L2,T,t2) is ported verbatim from
su2_XXX_L2.wb cells 2-10 to keep this file self-contained.*)
```

## Deleted Cell 12 — 2026-07-15 09:47:12
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

## Deleted Cell 11 — 2026-07-15 09:47:12
**Kind:** code
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

## Deleted Cell 10 — 2026-07-15 09:47:12
**Kind:** code
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

## Deleted Cell 9 — 2026-07-15 09:47:12
**Kind:** code
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

## Deleted Cell 8 — 2026-07-15 09:47:12
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

## Deleted Cell 7 — 2026-07-15 09:47:12
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

## Deleted Cell 6 — 2026-07-15 09:47:12
**Kind:** code
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

## Deleted Cell 5 — 2026-07-15 09:47:12
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

## Deleted Cell 4 — 2026-07-15 09:47:12
**Kind:** code
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

## Deleted Cell 3 — 2026-07-15 09:47:12
**Kind:** code
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

## Deleted Cell 2 — 2026-07-15 09:47:12
**Kind:** code
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

