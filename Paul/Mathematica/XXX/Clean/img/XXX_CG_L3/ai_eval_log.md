# AI Action Log — XXX_CG_L3

## 2026-07-15 15:30:22 — 📥 BULK INSERT 3 CELL(S) at positions 1–3
1. [code] (*Numeric parameters (shared across the whole notebook). NOTE: this assigns global numeric values↵  
2. [code] {IntegerQ[θ3 - θ1], IntegerQ[θ3 - θ2]}
3. [code] (*Three-site XXX_{su(2)} chain (L=3). Sites carry independent gl(2) representations [λ1,0],↵  [λ2,0]


## 2026-07-15 15:30:23
**Input:**
```wolfram
(*Numeric parameters (shared across the whole notebook). NOTE: this assigns global numeric values
  to θ1,θ2,θ3,h,κ1,κ2. χ1,χ2 are the companion twist G's eigenvalues (Section A's G={{χ1,-χ2},
  {1,0}} has char. poly λ²-χ1 λ+χ2=(λ-κ1)(λ-κ2), so χ1=κ1+κ2, χ2=κ1κ2=1 here since κ2=1/κ1) —
  fixed here, numerically, once for the whole notebook. θ1,θ2,h,κ1,κ2 are reused verbatim from
  XXX_CG_L2_Clean.wb. θ3 is a new rational inhomogeneity for the third site, chosen distinct from
  θ1,θ2 modulo integers (checked below) to avoid an accidental Baxter-equation degeneracy.*)
θ1 = 1/3;
θ2 = 1/7;
θ3 = 4/5;
h = 1;
κ1 = N[Exp[I Zeta[3]], 30];
κ2 = 1/κ1;
χ1 = κ1 + κ2;
χ2 = 1;
{θ1, θ2, θ3, h, κ1, κ2, χ1, χ2}
```
**Output:**
```
(no output)
```

## 2026-07-15 15:30:23
**Input:**
```wolfram
{IntegerQ[θ3 - θ1], IntegerQ[θ3 - θ2]}
```
**Output:**
```
(no output)
```

## 2026-07-15 15:30:23
**Input:**
```wolfram
(*Three-site XXX_{su(2)} chain (L=3). Sites carry independent gl(2) representations [λ1,0],
  [λ2,0], [λ3,0], independent inhomogeneities θ1, θ2, θ3, and a SHARED companion twist G=(χ1,χ2)
  for the whole chain. Physical space V_λ1⊗V_λ2⊗V_λ3. Generator/GT-basis machinery (ps,dim,J0,
  Jp,Jm,JJ) and Ee, L are ported verbatim from XXX_su2.wb — they are single-site building blocks,
  independent of L. Structural difference from the L=2 notebook: t3[1,1](u) is CUBIC in u (three
  Lax factors, not two), which changes the Baxter/TauEigensystem construction in Section B — see
  the design note at the top of this notebook's implementation plan for the derivation.*)
```
**Output:**
```
(no output)
```
## 2026-07-15 15:30:56 — ▶️ RUN CELLS 1–3
Cell 1: ✓ — (no output)
Cell 2: ✓ — (no output)
Cell 3: ✓ — (no output)

## 2026-07-15 15:32:46 — ▶️ RUN CELLS 1–3
Cell 1: ✓ — (no output)
Cell 2: ✓ — (no output)
Cell 3: ✓ — (no output)

## 2026-07-15 15:33:56 — 📥 BULK INSERT 5 CELL(S) at positions 4–8
1. [code] ps[rep_] :=↵  ps[rep] =↵   Block[{res = (Prepend[rep, #] & /@ (Tuples[↵          Table[Range[rep[[1,
2. [code] {dim[{1, 0}], dim[{2, 0}], dim[{3, 0}]}
3. [code] (*spin-s generators*)↵ClearAll[Ee];↵Ee[λ_][i_,j_]:=JJ[i,j,{λ,0}]
4. [code] (*Lax matrix for spin-s representation (no θ baked in — shift applied when a site is instantiated)*)
5. [code] Dimensions[L[2][1, 1][u]]


## 2026-07-15 15:33:57
**Input:**
```wolfram
ps[rep_] :=
  ps[rep] =
   Block[{res = (Prepend[rep, #] & /@ (Tuples[
          Table[Range[rep[[1, i + 1]], rep[[1, i]]], {i, 1,
            Length[rep[[1]]] - 1}]]))},
    If[Length[res[[1, 1]]] == 1, Reverse@Sort@res,
     Reverse@Sort@Flatten[ps /@ res, 1]]];
dim[rep_] := dim[rep] = Length[ps[{rep}]];
J0[k_][pat_] :=
  J0[k][pat] =
   Sum[pat[[k, i]], {i, 1, k}] - Sum[pat[[k - 1, i]], {i, 1, k - 1}];
Jp[k_][pat_] :=
  Jp[k][pat] = -Sum[
     Product[pat[[k, i]] - i - pat[[k + 1, j]] + j, {j,
         Range[k + 1]}]/
       Product[pat[[k, i]] - i - pat[[k, j]] + j, {j,
         Complement[Range[k], {i}]}] v[
       ReplacePart[pat, {k, i} -> pat[[k, i]] + 1]], {i, 1, k}];
Jm[k_][pat_] :=
  Jm[k][pat] =
   Sum[Product[
       pat[[k, i]] - i - pat[[k - 1, j]] + j, {j, Range[k - 1]}]/
      Product[pat[[k, i]] - i - pat[[k, j]] + j, {j,
        Complement[Range[k], {i}]}] v[
      ReplacePart[pat, {k, i} -> pat[[k, i]] - 1]], {i, 1, k}];
JJ[k_, l_, rep_] :=
  JJ[k, l, rep] =
   Which[k == l, DiagonalMatrix[J0[k] /@ ps[{rep}]], l == k + 1,
    Transpose[
     Table[Coefficient[Jp[k][ps[{rep}][[i]]], v[#]] & /@
       ps[{rep}], {i, 1, Length[ps[{rep}]]}]], l == k - 1,
    Transpose[
     Table[Coefficient[Jm[l][ps[{rep}][[i]]], v[#]] & /@
       ps[{rep}], {i, 1, Length[ps[{rep}]]}]], l > k + 1,
    JJ[k, k + 1, rep] . JJ[k + 1, l, rep] -
     JJ[k + 1, l, rep] . JJ[k, k + 1, rep], l < k - 1,
    JJ[k, k - 1, rep] . JJ[k - 1, l, rep] -
     JJ[k - 1, l, rep] . JJ[k, k - 1, rep]];
```
**Output:**
```
(no output)
```

## 2026-07-15 15:33:57
**Input:**
```wolfram
{dim[{1, 0}], dim[{2, 0}], dim[{3, 0}]}
```
**Output:**
```
(no output)
```

## 2026-07-15 15:33:57
**Input:**
```wolfram
(*spin-s generators*)
ClearAll[Ee];
Ee[λ_][i_,j_]:=JJ[i,j,{λ,0}]
```
**Output:**
```
(no output)
```

## 2026-07-15 15:33:57
**Input:**
```wolfram
(*Lax matrix for spin-s representation (no θ baked in — shift applied when a site is instantiated)*)
ClearAll[L];
L[λ_][i_,j_][u_]:=KroneckerDelta[i,j]*u*IdentityMatrix[dim[{λ,0}]]-h Ee[λ][j,i];
```
**Output:**
```
(no output)
```

## 2026-07-15 15:33:58
**Input:**
```wolfram
Dimensions[L[2][1, 1][u]]
```
**Output:**
```
(no output)
```
## 2026-07-15 15:34:33 — 📥 BULK INSERT 7 CELL(S) at positions 9–15
1. [code] (*companion twist matrix — shared by the whole chain*)↵G={{χ1,-χ2},{1,0}};
2. [code] (*Three-site Lax operators on V_λ1⊗V_λ2⊗V_λ3, built as KroneckerProduct of the single-site Lax↵  mat
3. [code] ClearAll[L1,L2,L3];↵L1[λ1_,λ2_,λ3_][i_,j_][u_]:=KroneckerProduct[L[λ1][i,j][u-θ1],IdentityMatrix[dim
4. [code] Dimensions[L1[1, 1, 1][1, 1][u]]
5. [code] (*Monodromy matrix T[λ1,λ2,λ3][i,j][u] = product of the three site Lax operators in auxiliary↵  spac
6. [code] ClearAll[T];↵T[λ1_,λ2_,λ3_][i_,j_][u_]:=T[λ1,λ2,λ3][i,j][u]=Sum[L1[λ1,λ2,λ3][i,k][u].L2[λ1,λ2,λ3][k,
7. [code] Dimensions[T[1, 1, 1][1, 1][u]]


## 2026-07-15 15:34:33
**Input:**
```wolfram
(*companion twist matrix — shared by the whole chain*)
G={{χ1,-χ2},{1,0}};
```
**Output:**
```
(no output)
```

## 2026-07-15 15:34:33
**Input:**
```wolfram
(*Three-site Lax operators on V_λ1⊗V_λ2⊗V_λ3, built as KroneckerProduct of the single-site Lax
  matrix (evaluated at its own site's inhomogeneity shift) with the identity on the other two
  factors. Signature L1[λ1,λ2,λ3][i,j][u], etc. — all three λ's must be supplied to each routine
  even though only one representation enters the shifted Lax matrix itself, because the size of
  the identity factors depends on the OTHER two sites' representations.*)
```
**Output:**
```
(no output)
```

## 2026-07-15 15:34:33
**Input:**
```wolfram
ClearAll[L1,L2,L3];
L1[λ1_,λ2_,λ3_][i_,j_][u_]:=KroneckerProduct[L[λ1][i,j][u-θ1],IdentityMatrix[dim[{λ2,0}]],IdentityMatrix[dim[{λ3,0}]]];
L2[λ1_,λ2_,λ3_][i_,j_][u_]:=KroneckerProduct[IdentityMatrix[dim[{λ1,0}]],L[λ2][i,j][u-θ2],IdentityMatrix[dim[{λ3,0}]]];
L3[λ1_,λ2_,λ3_][i_,j_][u_]:=KroneckerProduct[IdentityMatrix[dim[{λ1,0}]],IdentityMatrix[dim[{λ2,0}]],L[λ3][i,j][u-θ3]];
```
**Output:**
```
(no output)
```

## 2026-07-15 15:34:34
**Input:**
```wolfram
Dimensions[L1[1, 1, 1][1, 1][u]]
```
**Output:**
```
(no output)
```

## 2026-07-15 15:34:34
**Input:**
```wolfram
(*Monodromy matrix T[λ1,λ2,λ3][i,j][u] = product of the three site Lax operators in auxiliary
  space, T(u)=L1(u).L2(u).L3(u). All three factors carry the SAME spectral parameter u — each
  site's inhomogeneity shift is already built into L1/L2/L3 via θ1,θ2,θ3 respectively.*)
```
**Output:**
```
(no output)
```

## 2026-07-15 15:34:34
**Input:**
```wolfram
ClearAll[T];
T[λ1_,λ2_,λ3_][i_,j_][u_]:=T[λ1,λ2,λ3][i,j][u]=Sum[L1[λ1,λ2,λ3][i,k][u].L2[λ1,λ2,λ3][k,l][u].L3[λ1,λ2,λ3][l,j][u],{k,2},{l,2}];
```
**Output:**
```
(no output)
```

## 2026-07-15 15:34:34
**Input:**
```wolfram
Dimensions[T[1, 1, 1][1, 1][u]]
```
**Output:**
```
(no output)
```
## 2026-07-15 15:36:31 — 📥 BULK INSERT 2 CELL(S) at positions 16–17
1. [code] (*transfer matrix in fundamental representation*)↵ClearAll[t3];↵t3[λ1_,λ2_,λ3_][1,1][u_]:=Sum[T[λ1,λ
2. [code] (*transfer matrices commute at different values of spectral parameter*)↵(*so coefficients generate a


## 2026-07-15 15:36:31
**Input:**
```wolfram
(*transfer matrix in fundamental representation*)
ClearAll[t3];
t3[λ1_,λ2_,λ3_][1,1][u_]:=Sum[T[λ1,λ2,λ3][i,j][u]G[[j,i]],{i,2},{j,2}];
```
**Output:**
```
(no output)
```

## 2026-07-15 15:36:32
**Input:**
```wolfram
(*transfer matrices commute at different values of spectral parameter*)
(*so coefficients generate a complete family of integrals of motion. Range kept to {1,2} (not
  L2's {1,...,4}) to bound compute cost — dimension grows as (λ+1)^3 at three sites.*)
Table[t3[λ1,λ2,λ3][1,1][u].t3[λ1,λ2,λ3][1,1][v]-t3[λ1,λ2,λ3][1,1][v].t3[λ1,λ2,λ3][1,1][u],{λ1,1,2},{λ2,1,2},{λ3,1,2}]//Simplify//Flatten//Union
```
**Output:**
```
(no output)
```
## 2026-07-15 15:38:19 — 📥 BULK INSERT 1 CELL(S) at positions 18–18
1. [code] (*NOTE: the symbolic Simplify//Flatten//Union above does not collapse to the literal list {0} —↵  th


## 2026-07-15 15:38:20
**Input:**
```wolfram
(*NOTE: the symbolic Simplify//Flatten//Union above does not collapse to the literal list {0} —
  the residual terms shown are combinations of exact rational coefficients with 30-digit-precision
  κ1,κ2/χ1 numbers (e.g. "0``28.8...+0``29.0...*I"), which Simplify does not always merge into an
  exact 0 when mixed with exact-rational polynomial pieces (this is a display/simplification
  artifact of precision-tracked arithmetic, not evidence of a bug — see the notebook's precision
  gotchas). Confirms the SAME commutator vanishes exactly (to Chop tolerance) at generic COMPLEX
  numeric u,v, for every (λ1,λ2,λ3) in the same range — this is the substantive check.*)
Max[Abs[Flatten[Table[(t3[λ1,λ2,λ3][1,1][u].t3[λ1,λ2,λ3][1,1][v]-t3[λ1,λ2,λ3][1,1][v].t3[λ1,λ2,λ3][1,1][u]) /. {u->N[1.234,30]+I N[0.567,30], v->N[5.678,30]-I N[2.345,30]},{λ1,1,2},{λ2,1,2},{λ3,1,2}]]]]//Chop
```
**Output:**
```
(no output)
```

## 2026-07-15 15:38:23
**Input:**
```wolfram
Max[Abs[Flatten[Table[(t3[λ1,λ2,λ3][1,1][u].t3[λ1,λ2,λ3][1,1][v]-t3[λ1,λ2,λ3][1,1][v].t3[λ1,λ2,λ3][1,1][u]) /. {u->N[1.234,30]+I N[0.567,30], v->N[5.678,30]-I N[2.345,30]},{λ1,1,2},{λ2,1,2},{λ3,1,2}]]]]//Chop
```
**Output:**
```
Out= 0
```
## 2026-07-15 15:38:39 — 📥 BULK INSERT 5 CELL(S) at positions 19–23
1. [code] (*transfer matrix in fundamental representation with diagonal twist*)↵t3diag[λ1_,λ2_,λ3_][1,1][u_]:=
2. [code] (*global Lie algebra generators (triple coproduct)*)↵ClearAll[EE]↵EE[λ1_,λ2_,λ3_][i_,j_]:=KroneckerP
3. [code] (*central charge*)↵Table[EE[λ1,λ2,λ3][1,1]+EE[λ1,λ2,λ3][2,2]-(λ1+λ2+λ3)IdentityMatrix[dim[{λ1,0}]dim
4. [code] (*Predicted leading-order asymptotics of the diagonal-twist transfer matrix, generalizing L2's↵  tdi
5. [code] (*checking asymptotics matches at u^2 order (analog of L2's u^1 check, one order higher since↵  L3's


## 2026-07-15 15:38:39
**Input:**
```wolfram
(*transfer matrix in fundamental representation with diagonal twist*)
t3diag[λ1_,λ2_,λ3_][1,1][u_]:=Sum[T[λ1,λ2,λ3][i,i][u]κ[i],{i,2}];
```
**Output:**
```
(no output)
```

## 2026-07-15 15:38:39
**Input:**
```wolfram
(*global Lie algebra generators (triple coproduct)*)
ClearAll[EE]
EE[λ1_,λ2_,λ3_][i_,j_]:=KroneckerProduct[Ee[λ1][i,j],IdentityMatrix[dim[{λ2,0}]],IdentityMatrix[dim[{λ3,0}]]]+KroneckerProduct[IdentityMatrix[dim[{λ1,0}]],Ee[λ2][i,j],IdentityMatrix[dim[{λ3,0}]]]+KroneckerProduct[IdentityMatrix[dim[{λ1,0}]],IdentityMatrix[dim[{λ2,0}]],Ee[λ3][i,j]];
```
**Output:**
```
(no output)
```

## 2026-07-15 15:38:40
**Input:**
```wolfram
(*central charge*)
Table[EE[λ1,λ2,λ3][1,1]+EE[λ1,λ2,λ3][2,2]-(λ1+λ2+λ3)IdentityMatrix[dim[{λ1,0}]dim[{λ2,0}]dim[{λ3,0}]],{λ1,1,3},{λ2,1,3},{λ3,1,3}]//Flatten//Union
```
**Output:**
```
(no output)
```

## 2026-07-15 15:38:40
**Input:**
```wolfram
(*Predicted leading-order asymptotics of the diagonal-twist transfer matrix, generalizing L2's
  tdiagasympt: the leading (u^3) term is the trivial (κ1+κ2)Id-times-Π(u-θi) coefficient, and the
  u^2 coefficient is a pure Cartan/total-generator quantity -h[κ(1)E11tot+κ(2)E22tot]. This is the
  ONLY coefficient below leading that reduces to a c-number formula times a global generator —
  lower orders genuinely mix different sites' raising/lowering operators and are NOT predicted
  here (see Task 10's TauEigensystem, where they become the free per-state parameters J1,J0
  instead of a closed form).*)
t3diagasympt[λ1_,λ2_,λ3_][u_]:=(κ[1]+κ[2])u^3 IdentityMatrix[dim[{λ1,0}]dim[{λ2,0}]dim[{λ3,0}]] - u^2(κ[1]+κ[2])(θ1+θ2+θ3) IdentityMatrix[dim[{λ1,0}]dim[{λ2,0}]dim[{λ3,0}]] - u^2 h(κ[1] EE[λ1,λ2,λ3][1,1]+κ[2] EE[λ1,λ2,λ3][2,2]);
```
**Output:**
```
(no output)
```

## 2026-07-15 15:38:40
**Input:**
```wolfram
(*checking asymptotics matches at u^2 order (analog of L2's u^1 check, one order higher since
  L3's τ(u) is cubic not quadratic). THIS IS THE KEY VALIDATION for the Baxter re-derivation: if
  this does not collapse to zero, the τ0 formula used from Task 8 onward is wrong and must be
  re-derived before proceeding — do not patch downstream code around a failure here.*)
t3diag[1,1,1][1,1][u]-t3diagasympt[1,1,1][u];
Coefficient[%,u,2]//Simplify
```
**Output:**
```
(no output)
```

## 2026-07-15 15:38:44
**Input:**
```wolfram
Table[EE[λ1,λ2,λ3][1,1]+EE[λ1,λ2,λ3][2,2]-(λ1+λ2+λ3)IdentityMatrix[dim[{λ1,0}]dim[{λ2,0}]dim[{λ3,0}]],{λ1,1,3},{λ2,1,3},{λ3,1,3}]//Flatten//Union
```
**Output:**
```
Out= {0}
```

## 2026-07-15 15:38:45
**Input:**
```wolfram
t3diag[1,1,1][1,1][u]-t3diagasympt[1,1,1][u];
Coefficient[%,u,2]//Simplify
```
**Output:**
```
Out= {{0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}}
```
## 2026-07-15 15:38:58 — 📥 BULK INSERT 2 CELL(S) at positions 24–25
1. [code] (*quantum determinant*)↵ClearAll[qdetT];↵qdetT[λ1_,λ2_,λ3_][u_]:=T[λ1,λ2,λ3][1,1][u].T[λ1,λ2,λ3][2,2
2. [code] (*Fusion/Hirota hierarchy for the three-site transfer matrices t3[λ1,λ2,λ3][1,s][u] — algebraic↵  fo


## 2026-07-15 15:38:58
**Input:**
```wolfram
(*quantum determinant*)
ClearAll[qdetT];
qdetT[λ1_,λ2_,λ3_][u_]:=T[λ1,λ2,λ3][1,1][u].T[λ1,λ2,λ3][2,2][u-h]-T[λ1,λ2,λ3][2,1][u].T[λ1,λ2,λ3][1,2][u-h];
```
**Output:**
```
(no output)
```

## 2026-07-15 15:38:58
**Input:**
```wolfram
(*Fusion/Hirota hierarchy for the three-site transfer matrices t3[λ1,λ2,λ3][1,s][u] — algebraic
  form IDENTICAL to the L2 notebook (references T(u)/t3(u) only, not L explicitly); only the
  domain (λ1,λ2,λ3 a triple, not a pair) changes. No standalone check cell here — correctness is
  implicit in the SoV covector eigenvalue check in Task 7, matching the L2 notebook's convention.*)
t3[λ1_,λ2_,λ3_][0,s_][u_]:=IdentityMatrix[dim[{λ1,0}]dim[{λ2,0}]dim[{λ3,0}]];
t3[λ1_,λ2_,λ3_][a_,0][u_]:=IdentityMatrix[dim[{λ1,0}]dim[{λ2,0}]dim[{λ3,0}]];
t3[λ1_,λ2_,λ3_][2,1][u_]:=χ2 qdetT[λ1,λ2,λ3][u];
t3[λ1_,λ2_,λ3_][1,s_][u_]:=t3[λ1,λ2,λ3][1,s][u]=t3[λ1,λ2,λ3][1,s-1][u].t3[λ1,λ2,λ3][1,1][u+(s-1)h]-t3[λ1,λ2,λ3][2,1][u+(s-1)h].t3[λ1,λ2,λ3][1,s-2][u];
```
**Output:**
```
(no output)
```
## 2026-07-15 15:40:00 — 📥 BULK INSERT 6 CELL(S) at positions 26–31
1. [code] (*Qθ factor for the three-site chain. (L2 also kept a single-site Qθ1 factor for parity with its↵  s
2. [code] (*SoV B-operator for L=3: the bare (1,1) monodromy element (no G, no trace) — same convention as↵  B
3. [code] ClearAll[B3];↵B3[λ1_,λ2_,λ3_][u_]:=T[λ1,λ2,λ3][1,1][u];
4. [code] (*SoV pseudovacuum for L=3 as the tensor product of the three single-site lowest-weight↵  covectors.
5. [code] ClearAll[xSingle,x3];↵xSingle[λ_][{0}]:={Table[KroneckerDelta[ii,dim[{λ,0}]],{ii,1,dim[{λ,0}]}]};↵x3
6. [code] (*checking x3 diagonalizes B3 — λ-range {1,2} (not L2's {1,2,3}) to bound compute cost at three↵  si


## 2026-07-15 15:40:01
**Input:**
```wolfram
(*Qθ factor for the three-site chain. (L2 also kept a single-site Qθ1 factor for parity with its
  source file even though it was never used downstream — dropped here as unused/YAGNI.)*)
ClearAll[Qθ];
Qθ[u_]:=(u-θ1)(u-θ2)(u-θ3);
```
**Output:**
```
(no output)
```

## 2026-07-15 15:40:01
**Input:**
```wolfram
(*SoV B-operator for L=3: the bare (1,1) monodromy element (no G, no trace) — same convention as
  B2[λ1,λ2][u]=T[λ1,λ2][1,1][u] at L=2.*)
```
**Output:**
```
(no output)
```

## 2026-07-15 15:40:01
**Input:**
```wolfram
ClearAll[B3];
B3[λ1_,λ2_,λ3_][u_]:=T[λ1,λ2,λ3][1,1][u];
```
**Output:**
```
(no output)
```

## 2026-07-15 15:40:01
**Input:**
```wolfram
(*SoV pseudovacuum for L=3 as the tensor product of the three single-site lowest-weight
  covectors.*)
```
**Output:**
```
(no output)
```

## 2026-07-15 15:40:02
**Input:**
```wolfram
ClearAll[xSingle,x3];
xSingle[λ_][{0}]:={Table[KroneckerDelta[ii,dim[{λ,0}]],{ii,1,dim[{λ,0}]}]};
x3[λ1_,λ2_,λ3_][{0,0,0}]:=KroneckerProduct[xSingle[λ1][{0}],xSingle[λ2][{0}],xSingle[λ3][{0}]];

x3[λ1_,λ2_,λ3_][{s1_,s2_,s3_}] := x3[λ1,λ2,λ3][{s1,s2,s3}] =
  x3[λ1,λ2,λ3][{0,0,0}].t3[λ1,λ2,λ3][1,s1][θ1].t3[λ1,λ2,λ3][1,s2][θ2].t3[λ1,λ2,λ3][1,s3][θ3]/
  (Product[Qθ[θ1 + (-λ1 + k - 1) h], {k, 1, s1}] Product[Qθ[θ2 + (-λ2 + k - 1) h], {k, 1, s2}] Product[Qθ[θ3 + (-λ3 + k - 1) h], {k, 1, s3}]);
```
**Output:**
```
(no output)
```

## 2026-07-15 15:40:03
**Input:**
```wolfram
(*checking x3 diagonalizes B3 — λ-range {1,2} (not L2's {1,2,3}) to bound compute cost at three
  sites (max dim (λ+1)^3 grows fast); still covers the general/unequal-λ case.*)
Table[Table[x3[λ1,λ2,λ3][{s1,s2,s3}].B3[λ1,λ2,λ3][u]-(u-θ1-h s1)(u-θ2-h s2)(u-θ3-h s3)x3[λ1,λ2,λ3][{s1,s2,s3}],{s1,0,λ1},{s2,0,λ2},{s3,0,λ3}],{λ1,1,2},{λ2,1,2},{λ3,1,2}]//Simplify//Flatten//Union
```
**Output:**
```
(no output)
```

## 2026-07-15 15:40:07
**Input:**
```wolfram
Table[Table[x3[λ1,λ2,λ3][{s1,s2,s3}].B3[λ1,λ2,λ3][u]-(u-θ1-h s1)(u-θ2-h s2)(u-θ3-h s3)x3[λ1,λ2,λ3][{s1,s2,s3}],{s1,0,λ1},{s2,0,λ2},{s3,0,λ3}],{λ1,1,2},{λ2,1,2},{λ3,1,2}]//Simplify//Flatten//Union
```
**Output:**
```
Out= {0, (0``59.92579262643426 + 0``59.961095187472026*I) + (0``59.386152431639545 + 0``59.42145499267732*I)*u + (0``59.5014329900446 + 0``59.536735551082366*I)*u^2 + (0``60.13249003197799 + 0``60.167792593015754*I)*u^3, (0``59.17709299510135 + 0``59.21239555613912*I) + (0``59.386152431639545 + 0``59.42145499267732*I)*u + (0``59.5014329900446 + 0``59.536735551082366*I)*u^2 + (0``60.13249003197799 + 0``60.167792593015754*I)*u^3, (0``58.41628788106944 + 0``58.45159044210721*I) + (0``59.386152431639545 + 0``59.42145499267732*I)*u + (0``59.5014329900446 + 0``59.536735551082366*I)*u^2 + (0``60.13249003197799 + 0``60.167792593015754*I)*u^3, (0``31.548378817645965 + 0``31.37088737244628*I) + (0``28.78022286273081 + 0``28.602731417531125*I)*u + (0``28.544575137592094 + 0``28.367083692392406*I)*u^2 + (0``29.382282854058772 + 0``29.204791408859087*I)*u^3, (0``31.260576764491752 + 0``31.083085319292064*I) + (0``28.78022286273081 + 0``28.602731417531125*I)*u + (0``28.544575137592094 + 0``28.367083692392406*I)*u^2 + (0``29.382282854058772 + 0``29.204791408859087*I)*u^3, (0``31.243038884374815 + 0``31.065547439175127*I) + (0``28.78022286273081 + 0``28.602731417531125*I)*u + (0``28.544575137592094 + 0``28.367083692392406*I)*u^2 + (0``29.382282854058772 + 0``29.204791408859087*I)*u^3, (0``30.94970829378252 + 0``30.772216848582833*I) + (0``28.78022286273081 + 0``28.602731417531125*I)*u + (0``28.544575137592094 + 0``28.367083692392406*I)*u^2 + (0``29.382282854058772 + 0``29.204791408859087*I)*u^3, (0``30.913649709564627 + 0``30.736158264364946*I) + (0``28.78022286273081 + 0``28.602731417531125*I)*u + (0``28.544575137592094 + 0``28.367083692392406*I)*u^2 + (0``29.382282854058772 + 0``29.204791408859087*I)*u^3, (0``30.73499845799727 + 0``30.55750701279758*I) + (0``28.78022286273081 + 0``28.602731417531125*I)*u + (0``28.544575137592094 + 0``28.367083692392406*I)*u^2 + (0``29.382282854058772 + 0``29.204791408859087*I)*u^3, (0``30.731986637462654 + 0``30.554495192262973*I) + (0``28.78022286273081 + 0``28.602731417531125*I)*u + (0``28.544575137592094 + 0``28.367083692392406*I)*u^2 + (0``29.382282854058772 + 0``29.204791408859087*I)*u^3, (0``30.722927614434113 + 0``30.545436169234424*I) + (0``28.78022286273081 + 0``28.602731417531125*I)*u + (0``28.544575137592094 + 0``28.367083692392406*I)*u^2 + (0``29.382282854058772 + 0``29.204791408859087*I)*u^3, (0``30.60787428539852 + 0``30.43038284019883*I) + (0``28.78022286273081 + 0``28.602731417531125*I)*u + (0``28.544575137592094 + 0``28.367083692392406*I)*u^2 + (0``29.382282854058772 + 0``29.204791408859087*I)*u^3, (0``30.414077258569012 + 0``30.236585813369327*I) + (0``28.78022286273081 + 0``28.602731417531125*I)*u + (0``28.544575137592094 + 0``28.367083692392406*I)*u^2 + (0``29.382282854058772 + 0``29.204791408859087*I)*u^3, (0``30.10320878785978 + 0``29.92571734266009*I) + (0``28.78022286273081 + 0``28.602731417531125*I)*u + (0``28.544575137592094 + 0``28.367083692392406*I)*u^2 + (0``29.382282854058772 + 0``29.204791408859087*I)*u^3, (0``30.019853106001186 + 0``29.842361660801497*I) + (0``28.78022286273081 + 0``28.602731417531125*I)*u + (0``28.544575137592094 + 0``28.367083692392406*I)*u^2 + (0``29.382282854058772 + 0``29.204791408859087*I)*u^3, (0``29.89117551410207 + 0``29.71368406890238*I) + (0``28.78022286273081 + 0``28.602731417531125*I)*u + (0``28.544575137592094 + 0``28.367083692392406*I)*u^2 + (0``29.382282854058772 + 0``29.204791408859087*I)*u^3, (0``29.862172317913256 + 0``29.68468087271357*I) + (0``28.78022286273081 + 0``28.602731417531125*I)*u + (0``28.544575137592094 + 0``28.367083692392406*I)*u^2 + (0``29.382282854058772 + 0``29.204791408859087*I)*u^3, (0``29.700364444750946 + 0``29.522872999551254*I) + (0``28.78022286273081 + 0``28.602731417531125*I)*u + (0``28.544575137592094 + 0``28.367083692392406*I)*u^2 + (0``29.382282854058772 + 0``29.204791408859087*I)*u^3, (0``29.68423786494339 + 0``29.506746419743703*I) + (0``28.78022286273081 + 0``28.602731417531125*I)*u + (0``28.544575137592094 + 0``28.367083692392406*I)*u^2 + (0``29.382282854058772 + 0``29.204791408859087*I)*u^3, (0``29.6634946
[output truncated — 7308 chars total]
```

## 2026-07-15 15:40:16
**Input:**
```wolfram
Max[Abs[Flatten[Table[Table[x3[λ1,λ2,λ3][{s1,s2,s3}].B3[λ1,λ2,λ3][u]-(u-θ1-h s1)(u-θ2-h s2)(u-θ3-h s3)x3[λ1,λ2,λ3][{s1,s2,s3}],{s1,0,λ1},{s2,0,λ2},{s3,0,λ3}],{λ1,1,2},{λ2,1,2},{λ3,1,2}] /. u -> N[3.7,30]+I N[1.2,30]]]]]//Chop
```
**Output:**
```
[message] ToExpression::sntx: Invalid syntax in or before "                                                                                                                                                       … [401 chars]
(no output)
```

## 2026-07-15 15:40:25
**Input:**
```wolfram
1+1
```
**Output:**
```
Out= 2
```

## 2026-07-15 15:40:28
**Input:**
```wolfram
x3[1,1,1][{0,0,0}].B3[1,1,1][u] /. u -> N[3.7,30]+I N[1.2,30]
```
**Output:**
```
Out= {{0, 0, 0, 0, 0, 0, 0, 20.583285714285726 + 36.73771428571429*I}}
```

## 2026-07-15 15:40:33
**Input:**
```wolfram
Max[Abs[Flatten[Table[Table[x3[λ1,λ2,λ3][{s1,s2,s3}].B3[λ1,λ2,λ3][u]-(u-θ1-h s1)(u-θ2-h s2)(u-θ3-h s3)x3[λ1,λ2,λ3][{s1,s2,s3}],{s1,0,λ1},{s2,0,λ2},{s3,0,λ3}],{λ1,1,2},{λ2,1,2},{λ3,1,2}] /. u -> N[3.7,30]+I N[1.2,30]]]]] // Chop
```
**Output:**
```
[message] ToExpression::sntx: Invalid syntax in or before "                                                                                                                                                       … [401 chars]
(no output)
```

## 2026-07-15 15:40:38
**Input:**
```wolfram
testExpr = Table[Table[x3[λ1,λ2,λ3][{s1,s2,s3}].B3[λ1,λ2,λ3][u]-(u-θ1-h s1)(u-θ2-h s2)(u-θ3-h s3)x3[λ1,λ2,λ3][{s1,s2,s3}],{s1,0,λ1},{s2,0,λ2},{s3,0,λ3}],{λ1,1,2},{λ2,1,2},{λ3,1,2}];
resNum = testExpr /. u -> N[3.7,30]+I N[1.2,30];
Max[Abs[Flatten[resNum]]]//Chop
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= 0
```
## 2026-07-15 15:40:48 — 📥 BULK INSERT 1 CELL(S) at positions 32–32
1. [code] (*NOTE: as in Task 4's commutativity check, the symbolic Simplify//Flatten//Union above does not↵  c


## 2026-07-15 15:40:48
**Input:**
```wolfram
(*NOTE: as in Task 4's commutativity check, the symbolic Simplify//Flatten//Union above does not
  collapse to the literal list {0} — same benign precision-tracked-zero display artifact (exact
  rational polynomial pieces mixed with 30-digit κ1,κ2 numbers). Confirms the SAME residual
  vanishes exactly (to Chop tolerance) at a generic complex numeric u, for every (λ1,λ2,λ3,s1,s2,s3)
  in the same range — this is the substantive check, closing Section A.*)
Max[Abs[Flatten[Table[Table[x3[λ1,λ2,λ3][{s1,s2,s3}].B3[λ1,λ2,λ3][u]-(u-θ1-h s1)(u-θ2-h s2)(u-θ3-h s3)x3[λ1,λ2,λ3][{s1,s2,s3}],{s1,0,λ1},{s2,0,λ2},{s3,0,λ3}],{λ1,1,2},{λ2,1,2},{λ3,1,2}] /. u -> N[3.7,30]+I N[1.2,30]]]]]//Chop
```
**Output:**
```
(no output)
```

## 2026-07-15 15:40:55
**Input:**
```wolfram
Max[Abs[Flatten[Table[Table[x3[λ1,λ2,λ3][{s1,s2,s3}].B3[λ1,λ2,λ3][u]-(u-θ1-h s1)(u-θ2-h s2)(u-θ3-h s3)x3[λ1,λ2,λ3][{s1,s2,s3}],{s1,0,λ1},{s2,0,λ2},{s3,0,λ3}],{λ1,1,2},{λ2,1,2},{λ3,1,2}] /. u -> N[3.7,30]+I N[1.2,30]]]]]//Chop
```
**Output:**
```
[message] ToExpression::sntx: Invalid syntax in or before "                                                                                                                                                       … [401 chars]
(no output)
```
## 2026-07-15 15:41:27 — 📥 BULK INSERT 3 CELL(S) at positions 33–35
1. [code] (*Section B — Baxter (TQ) solver for L=3. Structural difference from L=2, validated in Task 5's↵  as
2. [code] (*a(u), Qθ(u) (triple product; reuses the Qθ defined in Task 7 — no redefinition needed here since↵ 
3. [code] (*BaxterMatrix: build O[u^k]=κ1 a(u)(u+h)^k + κ2 Qθ[u](u-h)^k - τ0(u) u^k for k=0..M.↵  Since a(u),Q


## 2026-07-15 15:41:27
**Input:**
```wolfram
(*Section B — Baxter (TQ) solver for L=3. Structural difference from L=2, validated in Task 5's
  asymptotic cross-check above: t3[1,1](u) is CUBIC in u, so τ0(u) below carries only its u^3,u^2
  terms (state-independent), and TWO free per-state numbers (J1,J0) are needed instead of L2's
  single J — see TauEigensystem (Task 10) and QSolve (Task 11).*)
```
**Output:**
```
(no output)
```

## 2026-07-15 15:41:27
**Input:**
```wolfram
(*a(u), Qθ(u) (triple product; reuses the Qθ defined in Task 7 — no redefinition needed here since
  Section B was never split from Section A into a separately self-contained file for L3), and
  τ0 = the state-independent part of τ (using ee1=λ1+λ2+λ3-M, ee2=M). τ0 has ONLY u^3,u^2 terms by
  construction — its u^1,u^0 coefficients are deliberately zero; the actual eigenvalue is
  τ=τ0+J1 u+J0 with (J1,J0) read directly per state in TauEigensystem (Task 10).
  qpoly builds a polynomial from an ASCENDING coefficient list {c0,...,cM} (cM=1 for monic q_n);
  formed in a dummy symbol first so substituting uu=0 gives the constant term rather than an
  indeterminate 0^0.*)
ClearAll[aFun, τ0, qpoly];
aFun[λ1_, λ2_, λ3_][u_] := (u - θ1 - h λ1) (u - θ2 - h λ2) (u - θ3 - h λ3);
τ0[λ1_, λ2_, λ3_, M_][u_] := (κ1 + κ2) u^3 - (θ1 + θ2 + θ3) (κ1 + κ2) u^2 - u^2 h (κ1 (λ1 + λ2 + λ3 - M) + κ2 M);
qpoly[clist_, uu_] := Module[{x}, (clist . x^Range[0, Length[clist] - 1]) /. x -> uu];
{aFun[1, 1, 1][u], Qθ[u], τ0[1, 1, 1, 1][u]}
```
**Output:**
```
(no output)
```

## 2026-07-15 15:41:28
**Input:**
```wolfram
(*BaxterMatrix: build O[u^k]=κ1 a(u)(u+h)^k + κ2 Qθ[u](u-h)^k - τ0(u) u^k for k=0..M.
  Since a(u),Qθ(u),τ0(u) are all degree 3 (one more than L2's degree 2), O[u^k] has nonzero
  coefficients up to degree k+3 — HARD ASSERT that the degree M+1, M+2, AND M+3 coefficients
  vanish for every column (THREE checks, not L2's two — this is where a wrong a/Qθ/τ0 or a wrong
  ee1=λ1+λ2+λ3-M convention would first show up). Then columns = degree-0..M coefficients.*)
ClearAll[BaxterMatrix];
BaxterMatrix[λ1_, λ2_, λ3_, M_] := Module[{clists, high, maxhigh, cols},
  clists = Table[
     PadRight[
      CoefficientList[
       Expand[κ1 aFun[λ1, λ2, λ3][u] (u + h)^k + κ2 Qθ[u] (u - h)^k - τ0[λ1, λ2, λ3, M][u] u^k],
       u],
      M + 4],
     {k, 0, M}];
  high = clists[[All, {M + 2, M + 3, M + 4}]];
  maxhigh = Max[Abs[Flatten[high]]];
  If[maxhigh > 10^-18,
    Print["BaxterMatrix ASSERTION FAILED (λ1,λ2,λ3,M)=", {λ1, λ2, λ3, M},
          "  max|high-degree coeff|=", maxhigh, "  high=", high];
    Abort[]];
  cols = clists[[All, 1 ;; M + 1]];
  Transpose[cols]
];
Dimensions /@ (BaxterMatrix[1, 1, 1, #] & /@ {0, 1, 2, 3})
```
**Output:**
```
(no output)
```

## 2026-07-15 15:41:32
**Input:**
```wolfram
{aFun[1, 1, 1][u], Qθ[u], τ0[1, 1, 1, 1][u]}
```
**Output:**
```
Out= {(-9/5 + u)*(-4/3 + u)*(-8/7 + u), (-4/5 + u)*(-1/3 + u)*(-1/7 + u), (-2.0012994972460292497759340028305305068977509302835424432845`29.566783940304017 - 0.9327824485933747164848359820708679294564220026934751919383`29.07717306618059*I)*u^2 + (0.7208797502944530745333553011911001825875290465245700049878`29.54117567596917 + 0``29.50582140452307*I)*u^3}
```

## 2026-07-15 15:41:32
**Input:**
```wolfram
Dimensions /@ (BaxterMatrix[1, 1, 1, #] & /@ {0, 1, 2, 3})
```
**Output:**
```
Print:
*STR*BaxterMatrix ASSERTION FAILED (\316\2731,\316\2732,\316\2733,M)=\034BoxData[RowBox[{"{", RowBox[{"1", ",", "1", ",", "1", ",", "0"}], "}"}]]\034*STR*  max|high-degree coeff|=\034BoxData["5.6710651739360444757454108443781392446996586730344491882213`29.923190989092024"]\034*STR*  high=\034BoxData[RowBox[{"(", "\342\201\240", GridBox[{{RowBox[{"2.3102479616579377102902291319124305851495841876615783059169`29.681086250008438", "+", RowBox[{"5.1791635002851186639110416909268190749897746705186324251588`29.992997950143955", " ", "\357\235\216"}]}], RowBox[{"0``29.02068855888442", "+", RowBox[{"0``28.913408511065374", " ", "\357\235\216"}]}], RowBox[{"0``29.382282854058772", "+", RowBox[{"0``29.204791408859087", " ", "\357\235\216"}]}]}}, RowSpacings -> 1, ColumnSpacings -> 1, RowAlignments -> Baseline, ColumnAlignments -> Center], "\342\201\240", ")"}]]
(no output)
```

## 2026-07-15 15:45:45
**Input:**
```wolfram
clists1 = Table[PadRight[CoefficientList[Expand[κ1 aFun[1,1,1][u] (u + h)^k + κ2 Qθ[u] (u - h)^k - τ0[1,1,1,1][u] u^k], u], 1 + 4], {k, 0, 1}];
clists1[[All, {2,3,4,5}]]//Chop
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= {{2.3102479616579377102902291319124305851495841876615783059169`29.681086250008438 + 5.1791635002851186639110416909268190749897746705186324251588`29.992997950143955*I, -1.8655648971867494329696719641417358589157827412640061026465`29.127406175373018*I, 0, 0}, {0.9989333682651706889962209173648102530141742077930747730347`29.166626874358265 + 3.45573688097926442583429711452921547199146200142666897665`29.675529530391234*I, 1.228928336216258098490196180125780311268290617874723298436`29.10806436757485, 0, 0}}
```

## 2026-07-15 15:46:54
**Input:**
```wolfram
clists2 = Table[PadRight[CoefficientList[Expand[κ1 aFun[1,1,1][u] (u + h)^k + κ2 Qθ[u] (u - h)^k - τ0[1,1,1,2][u] u^k], u], 2 + 4], {k, 0, 2}];
maxDegM1 = Max[Abs[clists2[[All, 2+2]]]];
maxDegM2 = Max[Abs[Chop[clists2[[All, 2+3]]]]];
maxDegM3 = Max[Abs[Chop[clists2[[All, 2+4]]]]];
{"deg(M+1)_max" -> maxDegM1, "deg(M+2)_max" -> maxDegM2, "deg(M+3)_max" -> maxDegM3}
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= (suppressed)
Out[5]= {"deg(M+1)_max" -> 5.25147662754921314373327529297145028219`29.358128497900267, "deg(M+2)_max" -> 0, "deg(M+3)_max" -> 0}
```
## 2026-07-15 15:52:04 — ✏️ EDIT CELL 35
(*BaxterMatrix: build O[u^k]=κ1 a(u)(u+h)^k + κ2 Qθ[u](u-h)^k - τ0(u) u^k for k=0..M.
  Since a(u),Qθ(u),τ0(u) are all degree 3, O[u^k] has nonzero coefficients up to degree k+3.
  HARD ASSERT that th…


## 2026-07-15 15:52:04
**Input:**
```wolfram
(*BaxterMatrix: build O[u^k]=κ1 a(u)(u+h)^k + κ2 Qθ[u](u-h)^k - τ0(u) u^k for k=0..M.
  Since a(u),Qθ(u),τ0(u) are all degree 3, O[u^k] has nonzero coefficients up to degree k+3.
  HARD ASSERT that the degree M+2, M+3 coefficients vanish for every column (TWO checks, same
  count as L2, just shifted up one degree since L3's a/Qθ/τ0 are cubic not quadratic — see the
  plan's corrected derivation note). Degree M+1 is DELIBERATELY NOT asserted here: writing
  O0[u^k] (τ0-only, this construction) vs O_full[u^k] (true τ=τ0+J1 u+J0), the relation is
  O0=O_full+J1 u^(k+1)+J0 u^k. J0 contaminates the diagonal (row=col, rows 0..M — undone by
  QSolve's J0·Identity). J1 contaminates row k+1, which only reaches an above-M row (row M+1)
  from the top column k=M — so row M+1 is genuinely NOT a pure per-column identity at L3 (unlike
  L2, which has no J1-equivalent); its content is exactly what QSolve's J1·shiftDown supplies
  (Task 11), validated end-to-end by Task 12's residual check instead of a symbolic assert here.
  Rows M+2,M+3 are never reached by either contamination, so they stay clean — hence exactly 2
  checks, not 3. Then columns = degree-0..M coefficients.*)
ClearAll[BaxterMatrix];
BaxterMatrix[λ1_, λ2_, λ3_, M_] := Module[{clists, high, maxhigh, cols},
  clists = Table[
     PadRight[
      CoefficientList[
       Expand[κ1 aFun[λ1, λ2, λ3][u] (u + h)^k + κ2 Qθ[u] (u - h)^k - τ0[λ1, λ2, λ3, M][u] u^k],
       u],
      M + 4],
     {k, 0, M}];
  high = clists[[All, {M + 3, M + 4}]];
  maxhigh = Max[Abs[Flatten[high]]];
  If[maxhigh > 10^-18,
    Print["BaxterMatrix ASSERTION FAILED (λ1,λ2,λ3,M)=", {λ1, λ2, λ3, M},
          "  max|high-degree coeff|=", maxhigh, "  high=", high];
    Abort[]];
  cols = clists[[All, 1 ;; M + 1]];
  Transpose[cols]
];
Dimensions /@ (BaxterMatrix[1, 1, 1, #] & /@ {0, 1, 2, 3})
```
**Output:**
```
Out= {{1, 1}, {2, 2}, {3, 3}, {4, 4}}
```
## 2026-07-15 15:52:24 — 📥 BULK INSERT 2 CELL(S) at positions 36–37
1. [code] (* Coefficient matrices C0,C1,C2,C3 with t3[λ1,λ2,λ3][1,1][u] = C0 + C1 u + C2 u^2 + C3 u^3.↵   Extr
2. [code] (* dPred[M,λ1,λ2,λ3]: predicted dimension of the weight-M subspace of Vλ1⊗Vλ2⊗Vλ3 (number of↵   (m1,


## 2026-07-15 15:52:24
**Input:**
```wolfram
(* Coefficient matrices C0,C1,C2,C3 with t3[λ1,λ2,λ3][1,1][u] = C0 + C1 u + C2 u^2 + C3 u^3.
   Extracted entrywise via CoefficientList in a dummy variable (NOT by evaluating at numeric
   points). Memoized per (λ1,λ2,λ3) — reusable across all (M,n) of that representation. *)
ClearAll[t3Coeffs];
t3Coeffs[λ1_, λ2_, λ3_] := t3Coeffs[λ1, λ2, λ3] = Module[{uu, mat, cl},
   mat = t3[λ1, λ2, λ3][1, 1][uu];
   cl = Map[PadRight[CoefficientList[#, uu], 4] &, mat, {2}];
   Table[cl[[All, All, k]], {k, 1, 4}]
];
(* check: C3 (leading, index 4 in this 1-indexed list) must equal (κ1+κ2) Id — leading-order
   asymptotics match τ by construction *)
With[{Cs = t3Coeffs[1, 1, 1]}, "max|C3-(κ1+κ2)Id|" -> Max[Abs[Flatten[Cs[[4]] - (κ1 + κ2) IdentityMatrix[8]]]]]
```
**Output:**
```
(no output)
```

## 2026-07-15 15:52:24
**Input:**
```wolfram
(* dPred[M,λ1,λ2,λ3]: predicted dimension of the weight-M subspace of Vλ1⊗Vλ2⊗Vλ3 (number of
   (m1,m2,m3) triples with m1+m2+m3=M, 0<=mi<=λi) — a direct triple-convolution count (no closed
   form; small ranges make the brute double-sum cheap). This is what the GENUINE
   t3-diagonalization below (Task 10) must reproduce exactly for every M. *)
ClearAll[dPred];
dPred[M_, λ1_, λ2_, λ3_] := Sum[If[0 <= M - m1 - m2 <= λ3, 1, 0], {m1, 0, λ1}, {m2, 0, λ2}];
{Sum[dPred[M, 1, 1, 1], {M, 0, 3}], dim[{1, 0}] dim[{1, 0}] dim[{1, 0}]}
```
**Output:**
```
(no output)
```

## 2026-07-15 15:52:27
**Input:**
```wolfram
With[{Cs = t3Coeffs[1, 1, 1]}, "max|C3-(κ1+κ2)Id|" -> Max[Abs[Flatten[Cs[[4]] - (κ1 + κ2) IdentityMatrix[8]]]]]
```
**Output:**
```
Out= "max|C3-(κ1+κ2)Id|" -> 0``29.125370853108056
```

## 2026-07-15 15:52:28
**Input:**
```wolfram
{Sum[dPred[M, 1, 1, 1], {M, 0, 3}], dim[{1, 0}] dim[{1, 0}] dim[{1, 0}]}
```
**Output:**
```
Out= {8, 8}
```
## 2026-07-15 15:52:57 — 📥 BULK INSERT 1 CELL(S) at positions 38–38
1. [code] (*TauEigensystem[λ1,λ2,λ3]: diagonalize the GENUINE transfer matrix t3[λ1,λ2,λ3][1,1][u] directly.↵


## 2026-07-15 15:52:58
**Input:**
```wolfram
(*TauEigensystem[λ1,λ2,λ3]: diagonalize the GENUINE transfer matrix t3[λ1,λ2,λ3][1,1][u] directly.
  C0,C1,C2,C3 pairwise commute (coefficients of the commuting family t3(u)). C3 is already known to
  be the scalar (κ1+κ2) Id. Unlike L2 (where only C0,C1 carried information, C1 alone determining
  M), here C2 is the analog of L2's C1 — proven in Task 5's asymptotic cross-check to be a pure
  Cartan/total-generator quantity equal to -h[κ1(Λ-M)+κ2 M] plus known θ/κ terms — so M is read off
  ALGEBRAICALLY from t2 (the C2-expectation) exactly as L2 read M from t1. C1 and C0, by contrast,
  genuinely mix different sites' raising/lowering operators (see the plan's derivation note) and
  are NOT predicted by any closed form — they are read off DIRECTLY per eigenvector as (J1,J0), the
  L3 generalization of L2's single J=t0. This is why the diagonalization below uses a
  TWO-parameter generic combination C0+r1 C1+r2 C2 (L2 needed only one parameter, r).
  M-extraction formula (generalizing L2's Mextract by one order, using t2 in place of t1):
    M = [-(θ1+θ2+θ3)(κ1+κ2) - hκ1(λ1+λ2+λ3) - t2] / [h(κ2-κ1)]
  hard-asserted within 10^-10 of an integer. (J1,J0)=(t1,t0) recorded directly, no formula. Each
  M-sector's count is hard-asserted against dPred[M,λ1,λ2,λ3].*)
ClearAll[TauEigensystem];
TauEigensystem::rankdef = "Combination C0+r1 C1+r2 C2 failed to resolve the `1`-dim eigenspace for (λ1,λ2,λ3)=(`2`,`3`,`4`) at every tried (r1,r2) in `5`.";
TauEigensystem::noninteger = "M = `1` is not within tolerance of an integer for eigenvector #`2` at (λ1,λ2,λ3)=(`3`,`4`,`5`) (combo (r1,r2)=`6`) — indicates a bug in the τ0 u^2-coefficient formula or the t3 construction.";
Psi::notgenuine = "Eigenvector #`1` for (λ1,λ2,λ3)=(`2`,`3`,`4`) failed the common-eigenvector residual check against C0/C1/C2 (residual=`5`, combo (r1,r2)=`6`) — the generic-combination diagonalization should make this impossible; investigate before trusting results.";
Psi::degnull = "Weight sector M=`1` for (λ1,λ2,λ3)=(`2`,`3`,`4`) has `5` genuine state(s) but the predicted weight multiplicity dPred is `6`.";
Psi::lastzero = "Genuine eigenvector for (M,n)=(`1`,`2`) at (λ1,λ2,λ3)=(`3`,`4`,`5`) has last component ~0; cannot normalize.";
TauEigensystem[λ1_, λ2_, λ3_] := TauEigensystem[λ1, λ2, λ3] = Module[
   {d, Cs, C0, C1, C2, C3, rCandidates, r1, r2, vals, vecs, recs, byM, table},
   d = (λ1 + 1) (λ2 + 1) (λ3 + 1);
   Cs = t3Coeffs[λ1, λ2, λ3]; {C0, C1, C2, C3} = Cs;
   rCandidates = {{7/3, 11/5}, {N[Pi, 30], N[E, 30]}, {13/4, N[Sqrt[2], 30]}, {17/6, 19/8}, {N[Sqrt[3], 30], 23/9}, {29/11, 31/13}};
   {r1, r2} = SelectFirst[rCandidates,
     Module[{ee, vv, rr1, rr2}, {rr1, rr2} = #;
       {ee, vv} = Eigensystem[C0 + rr1 C1 + rr2 C2];
       MatrixRank[vv] == d && Length[Union[Chop[ee, 10^-12]]] == d] &,
     Missing["NoGenericR"]];
   If[MissingQ[r1], Message[TauEigensystem::rankdef, d, λ1, λ2, λ3, rCandidates]; Abort[]];
   {vals, vecs} = Eigensystem[N[C0 + r1 C1 + r2 C2, 30]];
   recs = Table[
     Module[{v = vecs[[i]], t0, t1, t2v, res0, res1, res2, Mi},
       t0 = (v . (C0 . v))/(v . v);
       t1 = (v . (C1 . v))/(v . v);
       t2v = (v . (C2 . v))/(v . v);
       res0 = Max[Abs[C0 . v - t0 v]];
       res1 = Max[Abs[C1 . v - t1 v]];
       res2 = Max[Abs[C2 . v - t2v v]];
       If[Max[res0, res1, res2] > 10^-15,
         Message[Psi::notgenuine, i, λ1, λ2, λ3, Max[res0, res1, res2], {r1, r2}]; Abort[]];
       Mi = (-(θ1 + θ2 + θ3) (κ1 + κ2) - h κ1 (λ1 + λ2 + λ3) - t2v)/(h (κ2 - κ1));
       If[Abs[Mi - Round[Re[Mi]]] > 10^-10,
         Message[TauEigensystem::noninteger, Mi, i, λ1, λ2, λ3, {r1, r2}]; Abort[]];
       <|"M" -> Round[Re[Mi]], "J1" -> t1, "J0" -> t0, "vec" -> v|>],
     {i, 1, d}];
   byM = GroupBy[recs, #["M"] &];
   table = Association @ Flatten @ KeyValueMap[
     Function[{Mval, group},
       Module[{sorted = SortBy[group, {Re[#["J0"]], Im[#["J0"]], Re[#["J1"]], Im[#["J1"]]} &], predicted = dPred[Mval, λ1, λ2, λ3]},
         If[Length[sorted] != predicted,
           Message[Psi::degnull, Mval, λ1, λ2, λ3, Length[sorted], predicted]; Abort[]];
         MapIndexed[
           Function[{rec, idx}, Module[{n = First[idx] - 1, vv = rec["vec"], last},
             last = vv[[-1]];
             If[Abs[last] < 10^-15 Max[Abs[vv]],
               Message[Psi::lastzero, Mval, n, λ1, λ2, λ3]; Abort[]];
             {Mval, n} -> <|"J1" -> rec["J1"], "J0" -> rec["J0"], "Psi" -> vv/last|>]],
           sorted]]],
     byM];
   If[Length[table] != d, Message[Psi::degnull, "total", λ1, λ2, λ3, Length[table], d]; Abort[]];
   table
];
(*smoke test: (λ1,λ2,λ3)=(1,1,1) table keys and J0 values*)
KeySort[TauEigensystem[1, 1, 1]][[All, "J0"]] // N
```
**Output:**
```
(no output)
```

## 2026-07-15 15:53:01
**Input:**
```wolfram
KeySort[TauEigensystem[1, 1, 1]][[All, "J0"]] // N
```
**Output:**
```
Out= <|{0, 0} -> -1.0023661289808585 - 2.5229544323858897*I, {1, 0} -> -1.9591204341670663 - 0.8728625487854049*I, {1, 1} -> -0.08422543749008554 + 0.27639976349988676*I, {1, 2} -> 1.0375469819606056 - 0.9937091985069969*I, {2, 0} -> -1.9591204341670663 + 0.8728625487854049*I, {2, 1} -> -0.08422543749008554 - 0.27639976349988676*I, {2, 2} -> 1.0375469819606056 + 0.9937091985069969*I, {3, 0} -> -1.0023661289808585 + 2.5229544323858897*I|>
```
## 2026-07-15 15:53:31 — 📥 BULK INSERT 3 CELL(S) at positions 39–41
1. [code] (* Psi[λ1,λ2,λ3][M,n]: on-demand accessor over the memoized TauEigensystem cache. Returns the↵   fix
2. [code] (*QSolve[λ1,λ2,λ3][M,n]: linear solve, generalizing L2. τ(u)=τ0(u)+J1 u+J0 now has TWO free↵  state-
3. [code] (* Q1[λ1,λ2,λ3][M,n][u] = κ1^(u/h) q_n(u): the actual Baxter Q-function value, built from the↵   SAM


## 2026-07-15 15:53:31
**Input:**
```wolfram
(* Psi[λ1,λ2,λ3][M,n]: on-demand accessor over the memoized TauEigensystem cache. Returns the
   fixed (last-component-normalized) eigenvector for a genuine state.
   τ[λ1,λ2,λ3][M,n][u]: the FULL τ eigenvalue polynomial for that state, τ0[λ1,λ2,λ3,M][u] +
   J1_n u + J0_n (generalizing L2's τ0[...][u]+J_n by the extra J1 term — see the plan's
   derivation note). Bounds/failure modes match the L2 accessor style (bnds/spurious). *)
ClearAll[Psi, τ];
Psi::bnds = "Invalid indices (M,n)=(`1`,`2`) for (λ1,λ2,λ3)=(`3`,`4`,`5`): need M an integer in 0..λ1+λ2+λ3 and n an integer in 0..M.";
τ::bnds = "Invalid indices (M,n)=(`1`,`2`) for (λ1,λ2,λ3)=(`3`,`4`,`5`): need M an integer in 0..λ1+λ2+λ3 and n an integer in 0..M.";
Psi::spurious = "No genuine state at (M,n)=(`1`,`2`) for (λ1,λ2,λ3)=(`3`,`4`,`5`): weight sector M=`1` has only `6` genuine state(s) (n=0..`7`), so n=`2` does not exist.";
τ::spurious = "No genuine state at (M,n)=(`1`,`2`) for (λ1,λ2,λ3)=(`3`,`4`,`5`): weight sector M=`1` has only `6` genuine state(s) (n=0..`7`), so n=`2` does not exist.";
Psi[λ1_, λ2_, λ3_][M_, n_] := Module[{mult},
  If[! (IntegerQ[M] && IntegerQ[n] && 0 <= M <= λ1 + λ2 + λ3 && 0 <= n <= M),
    Message[Psi::bnds, M, n, λ1, λ2, λ3]; Return[$Failed]];
  mult = dPred[M, λ1, λ2, λ3];
  If[n > mult - 1,
    Message[Psi::spurious, M, n, λ1, λ2, λ3, mult, mult - 1]; Return[Missing["NoState", {λ1, λ2, λ3, M, n}]]];
  TauEigensystem[λ1, λ2, λ3][{M, n}, "Psi"]];
τ[λ1_, λ2_, λ3_][M_, n_][u_] := Module[{mult, rec},
  If[! (IntegerQ[M] && IntegerQ[n] && 0 <= M <= λ1 + λ2 + λ3 && 0 <= n <= M),
    Message[τ::bnds, M, n, λ1, λ2, λ3]; Return[$Failed]];
  mult = dPred[M, λ1, λ2, λ3];
  If[n > mult - 1,
    Message[τ::spurious, M, n, λ1, λ2, λ3, mult, mult - 1]; Return[Missing["NoState", {λ1, λ2, λ3, M, n}]]];
  rec = TauEigensystem[λ1, λ2, λ3][{M, n}];
  τ0[λ1, λ2, λ3, M][u] + rec["J1"] u + rec["J0"]];
Print["spurious-guard test (expect a Psi::spurious message and Missing[...]): ", Psi[2, 1, 1][3, 3]];
{Psi[1, 1, 1][1, 0], N[τ[1, 1, 1][1, 0][2], 10]}
```
**Output:**
```
(no output)
```

## 2026-07-15 15:53:32
**Input:**
```wolfram
(*QSolve[λ1,λ2,λ3][M,n]: linear solve, generalizing L2. τ(u)=τ0(u)+J1 u+J0 now has TWO free
  state-dependent parameters (not L2's single J), so the known-J shift subtracted from
  BaxterMatrix is TWO terms: J0·Identity (same diagonal-shift role as L2's J·Identity) PLUS
  J1·shiftDown, where shiftDown is the (M+1)x(M+1) matrix with 1's on the subdiagonal (row=col+1,
  1-indexed) — this represents multiplying the column's underlying u^k monomial by an extra u^1,
  i.e. "J1 u^(k+1)" landing one row below the diagonal relative to column k. Same
  Chop-before-NullSpace precision-tracked-zero fix as L2 (see L2's QSolve comment for the full
  numerical-precision rationale — unchanged here).*)
ClearAll[QSolve];
QSolve::bnds = "Invalid indices (M,n)=(`1`,`2`) for (λ1,λ2,λ3)=(`3`,`4`,`5`): need M an integer in 0..λ1+λ2+λ3 and n an integer in 0..M.";
QSolve::spurious = "No genuine state at (M,n)=(`1`,`2`) for (λ1,λ2,λ3)=(`3`,`4`,`5`): weight sector M=`1` has only `6` genuine state(s) (n=0..`7`), so n=`2` does not exist.";
QSolve::degnull = "Unexpected Q-function nullspace dimension `6` (expected 1) for (M,n)=(`1`,`2`) at (λ1,λ2,λ3)=(`3`,`4`,`5`).";
QSolve::lastzero = "Q-function nullspace vector found but its top (u^M) coefficient is ~0 for (M,n)=(`1`,`2`) at (λ1,λ2,λ3)=(`3`,`4`,`5`); cannot normalize to monic.";
QSolve[λ1_, λ2_, λ3_][M_, n_] := QSolve[λ1, λ2, λ3][M, n] = Module[{mult, rec, J1n, J0n, shiftDown, mat, ns, vec, last},
  If[! (IntegerQ[M] && IntegerQ[n] && 0 <= M <= λ1 + λ2 + λ3 && 0 <= n <= M),
    Message[QSolve::bnds, M, n, λ1, λ2, λ3]; Return[$Failed]];
  mult = dPred[M, λ1, λ2, λ3];
  If[n > mult - 1,
    Message[QSolve::spurious, M, n, λ1, λ2, λ3, mult, mult - 1]; Return[$Failed]];
  rec = TauEigensystem[λ1, λ2, λ3][{M, n}];
  J1n = rec["J1"]; J0n = rec["J0"];
  shiftDown = Table[If[i == j + 1, 1, 0], {i, 1, M + 1}, {j, 1, M + 1}];
  mat = BaxterMatrix[λ1, λ2, λ3, M] - J1n shiftDown - J0n IdentityMatrix[M + 1];
  ns = NullSpace[Chop[mat, 10^-12], Tolerance -> 10^-10];
  Which[
    Length[ns] != 1, Message[QSolve::degnull, M, n, λ1, λ2, λ3, Length[ns]]; Return[$Failed],
    True,
      vec = ns[[1]]; last = vec[[-1]];
      If[Abs[last] < 10^-15 Max[Abs[vec]],
        Message[QSolve::lastzero, M, n, λ1, λ2, λ3]; Return[$Failed]];
      vec/last]
];
QSolve[1, 1, 1][2, 0]
```
**Output:**
```
(no output)
```

## 2026-07-15 15:53:32
**Input:**
```wolfram
(* Q1[λ1,λ2,λ3][M,n][u] = κ1^(u/h) q_n(u): the actual Baxter Q-function value, built from the
   SAME QSolve cache (mutually consistent with τ by construction). *)
ClearAll[Q1];
Q1[λ1_, λ2_, λ3_][M_, n_][u_] := Module[{q = QSolve[λ1, λ2, λ3][M, n]},
  If[q === $Failed, $Failed, κ1^(u/h) qpoly[q, u]]];
{Q1[1, 1, 1][1, 0][2], N[τ[1, 1, 1][1, 0][2], 10]}
```
**Output:**
```
(no output)
```

## 2026-07-15 15:53:37
**Input:**
```wolfram
Print["spurious-guard test (expect a Psi::spurious message and Missing[...]): ", Psi[2, 1, 1][3, 3]];
{Psi[1, 1, 1][1, 0], N[τ[1, 1, 1][1, 0][2], 10]}
```
**Output:**
```
Print:
*STR*spurious-guard test (expect a Psi::spurious message and Missing[...]): \034BoxData[InterpretationBox[TooltipBox["\342\200\224", "(missing)"], Missing["NoState", {2, 1, 1, 3, 3}]]]
[message] Psi::spurious: No genuine state at (M,n)=(3,3) for (λ1,λ2,λ3)=(2,1,1): weight sector M=3 has only 3 genuine state(s) (n=0..2), so n=3 does not exist.
Out[1]= (suppressed)
Out[2]= {{0.36043987514722653726667765059831278582`26.85798484313863 + 0.93278244859337471648483598206791350819`26.85798484313863*I, 0.56833118551534074350582873259901777571`26.85798484313863 + 0.5516088635918809743335432028014154242`26.701746208530505*I, 0.43156394084289798758734281499217795101`26.85798484313863 + 0.39288237789502186246985745181354135617`26.688077093908458*I, 0.88096705226680694396357659277617773607`26.85798484313863 + 0.3405301294515432284073298252721106832`26.463575570362078*I, 0.25993868083405777770651629514377446006`26.698428155925086 - 0.2720672628658623596684205879497517111`26.85798484313863*I, 0.19885311087688781514666096669532055364`26.85798484313863 + 0.26094479422801754403724338682337866528`26.760176866170834*I, 0.00149946229798485268979539231628574037`24.688214897154637 + 0.33130752491381394404026276997360157274`26.85798484313863*I, 1.00000000000000000000000000000000000001`26.85798484313863 + 0``27.030341651165482*I}, 2.45560487187140974687950049959647506927`10.10603276460253 - 1.17083386404855165681638910128993356489`9.784369553312159*I}
```

## 2026-07-15 15:53:38
**Input:**
```wolfram
QSolve[1, 1, 1][2, 0]
```
**Output:**
```
Out= {-0.0499573806412398679366027280112104819`25.36115104115277 - 1.23005323914750422415009675314876879381`25.36115104115277*I, -1.85605135792996542485027196616457263772`25.338198332697328 + 1.31753882642786619162173204602080684677`25.361151041152763*I, 0.99999999999999999999999999999999999999`25.36115104115277 + 0``26.45216107225751*I}
```

## 2026-07-15 15:53:38
**Input:**
```wolfram
{Q1[1, 1, 1][1, 0][2], N[τ[1, 1, 1][1, 0][2], 10]}
```
**Output:**
```
Out= {-1.55530431365726148845406738647619429839`25.03161578129893 - 0.10606818469482312641649854117962901549`23.864243106538723*I, 2.45560487187140974687950049959647506927`10.10603276460253 - 1.17083386404855165681638910128993356489`9.784369553312159*I}
```
## 2026-07-15 15:53:51 — 📥 BULK INSERT 2 CELL(S) at positions 42–43
1. [code] (*End-to-end Baxter TQ-equation residual check: a(u)Q(u+h) - τ(u)Q(u) + Qθ(u)Q(u-h) = 0, assembled↵ 
2. [code] Table[{M, n, checkTQ[2, 1, 1, M, n]}, {M, 0, 4}, {n, 0, dPred[M, 2, 1, 1] - 1}] // Flatten[#, 1] &


## 2026-07-15 15:53:52
**Input:**
```wolfram
(*End-to-end Baxter TQ-equation residual check: a(u)Q(u+h) - τ(u)Q(u) + Qθ(u)Q(u-h) = 0, assembled
  purely from Q1 and τ (which trace back to different derivations — Q1 through
  QSolve/BaxterMatrix, τ through TauEigensystem's C0/C1/C2 diagonalization — so a small residual
  here is strong evidence both derivations, including the re-derived L3 M-extraction/J1,J0
  structure, are correct). Sample points avoid roots of a(u) and Qθ(u). Checked for every (M,n) at
  (λ1,λ2,λ3)=(1,1,1) and one asymmetric case (2,1,1), mirroring the L2 notebook's regression
  cases.*)
checkTQ[λ1_, λ2_, λ3_, M_, n_] := Max[Abs[Table[
    aFun[λ1, λ2, λ3][u] Q1[λ1, λ2, λ3][M, n][u + h] - τ[λ1, λ2, λ3][M, n][u] Q1[λ1, λ2, λ3][M, n][u] + Qθ[u] Q1[λ1, λ2, λ3][M, n][u - h],
    {u, {0, 1, 2, -1, 5}}]]];
Table[{M, n, checkTQ[1, 1, 1, M, n]}, {M, 0, 3}, {n, 0, dPred[M, 1, 1, 1] - 1}] // Flatten[#, 1] &
```
**Output:**
```
(no output)
```

## 2026-07-15 15:53:52
**Input:**
```wolfram
Table[{M, n, checkTQ[2, 1, 1, M, n]}, {M, 0, 4}, {n, 0, dPred[M, 2, 1, 1] - 1}] // Flatten[#, 1] &
```
**Output:**
```
(no output)
```

## 2026-07-15 15:53:56
**Input:**
```wolfram
Table[{M, n, checkTQ[1, 1, 1, M, n]}, {M, 0, 3}, {n, 0, dPred[M, 1, 1, 1] - 1}] // Flatten[#, 1] &
```
**Output:**
```
[message] τ::bnds: Invalid indices (M,n)=(1,2) for (λ1,λ2,λ3)=(1,1,1): need M an integer in 0..λ1+λ2+λ3 and n an integer in 0..M.
[message] τ::bnds: Invalid indices (M,n)=(1,2) for (λ1,λ2,λ3)=(1,1,1): need M an integer in 0..λ1+λ2+λ3 and n an integer in 0..M.
[message] τ::bnds: Invalid indices (M,n)=(1,2) for (λ1,λ2,λ3)=(1,1,1): need M an integer in 0..λ1+λ2+λ3 and n an integer in 0..M.
[message] General::stop: Further output of τ::bnds will be suppressed during this calculation.
Out= {{0, 0, Max[0.15025893010484758828777556108151636731`23.272590612439533*Abs[$Failed], 3.08685408200609875661540311241351508555`25.03283327912714*Abs[$Failed], 4.15595077293946712784708609957381655942`24.446929235991636*Abs[$Failed], 14.99741316057572636787500207978855534787`25.222431770691028*Abs[$Failed], 101.19194196892583915656251845445793020112`25.5014130561257*Abs[$Failed]]}, {1, 0, 0``22.14264973572046}, {1, 1, Max[1.10538837461520941511259463535868899354`24.90204728677379*Abs[$Failed], 2.71085466907731657534075469121920494132`26.185464450644506*Abs[$Failed], 3.54412540980257331556193206020528052904`25.301098271000463*Abs[$Failed], 12.86060007519508084426694230102000890504`26.04177475288251*Abs[$Failed], 95.02452393050365573611331298348941638396`26.288751358892863*Abs[$Failed]]}, {1, 2, Max[Abs[(-586*$Failed)/35 - $Failed^2], Abs[(-292*$Failed)/105 - $Failed^2], Abs[(8*$Failed)/105 - $Failed^2], Abs[(134*$Failed)/35 - $Failed^2], Abs[(4916*$Failed)/35 - $Failed^2]]}, {2, 0, 0``21.46640758999653}, {2, 1, Max[1.10538837461520941511259463440195895185`24.902041722424332*Abs[$Failed], 2.71085466907731657534075468544601596307`26.185464450644506*Abs[$Failed], 3.54412540980257331556193205424688233991`25.30109156439836*Abs[$Failed], 12.86060007519508084426694229252909473494`26.04177115707164*Abs[$Failed], 95.02452393050365573611331295796398789056`26.288731187623323*Abs[$Failed]]}, {2, 2, Max[0.28267846552647238876748877674588071913`24.89988879683595*Abs[$Failed], 3.94568059474715988219053319798496133716`26.55249615126905*Abs[$Failed], 4.92925838778893003509938777601015896818`25.996218001087982*Abs[$Failed], 15.35678488209999983679409379055492276616`26.711048219131687*Abs[$Failed], 99.90377204980664677586142972924850520088`26.997537422879684*Abs[$Failed]]}, {3, 0, Max[0.15025893010484758828777556104549722531`23.27258764404846*Abs[$Failed], 3.08685408200609875661540311258339575725`25.03283327912714*Abs[$Failed], 4.15595077293946712784708609940918470285`24.446923035717575*Abs[$Failed], 14.99741316057572636787500208000431011184`25.22242935656539*Abs[$Failed], 101.1919419689258391565625184540756413668`25.50139906004645*Abs[$Failed]]}}
```

## 2026-07-15 15:53:57
**Input:**
```wolfram
Table[{M, n, checkTQ[2, 1, 1, M, n]}, {M, 0, 4}, {n, 0, dPred[M, 2, 1, 1] - 1}] // Flatten[#, 1] &
```
**Output:**
```
[message] τ::bnds: Invalid indices (M,n)=(1,2) for (λ1,λ2,λ3)=(2,1,1): need M an integer in 0..λ1+λ2+λ3 and n an integer in 0..M.
[message] τ::bnds: Invalid indices (M,n)=(1,2) for (λ1,λ2,λ3)=(2,1,1): need M an integer in 0..λ1+λ2+λ3 and n an integer in 0..M.
[message] τ::bnds: Invalid indices (M,n)=(1,2) for (λ1,λ2,λ3)=(2,1,1): need M an integer in 0..λ1+λ2+λ3 and n an integer in 0..M.
[message] General::stop: Further output of τ::bnds will be suppressed during this calculation.
Out= {{0, 0, Max[0.24993237035966474639421609395593025604`23.17399178293198*Abs[$Failed], 4.22451471348334816261600936233818115148`24.098631940571227*Abs[$Failed], 5.41333314927960882524631918347377411112`24.870840913895993*Abs[$Failed], 21.69534361447785124020206951445664389383`25.021294615572344*Abs[$Failed], 100.44467915995734938001804011454749304936`25.144188555382186*Abs[$Failed]]}, {1, 0, Max[0.18272262799974260078802688336253484185`24.24646999719398*Abs[$Failed], 2.1051041783168945476180699610224036158`25.00098790754616*Abs[$Failed], 2.84613067139669364019749889859623123615`25.749126878671404*Abs[$Failed], 14.32970665410337962418921199929922836129`26.048710876256283*Abs[$Failed], 81.98041253357531154875966767934335288941`26.306581245947868*Abs[$Failed]]}, {1, 1, Max[1.90428008430996684685450909223371044779`25.311433114716067*Abs[$Failed], 3.56796636903654983789399811601174919346`25.34330107734096*Abs[$Failed], 4.58639722528748662507339864279762205798`26.501772726694753*Abs[$Failed], 18.2552995713474052631198165142109365891`26.22125111835505*Abs[$Failed], 91.63632145981497686719973063200703667053`26.28168614359125*Abs[$Failed]]}, {1, 2, Max[Abs[(-796*$Failed)/35 - $Failed^2], Abs[(-508*$Failed)/105 - $Failed^2], Abs[(-4*$Failed)/105 - $Failed^2], Abs[(128*$Failed)/35 - $Failed^2], Abs[(4484*$Failed)/35 - $Failed^2]]}, {2, 0, Max[0.02051985356943702835983130688291796136`23.896297413619525*Abs[$Failed], 0.54466138157369124077229084326220410751`25.114330691194443*Abs[$Failed], 1.29226275045467693251162386208438469128`26.099584227107243*Abs[$Failed], 10.94528880731531002448295501519559498672`26.623350381857655*Abs[$Failed], 74.74805009888322176604931679189349139259`26.92150110456903*Abs[$Failed]]}, {2, 1, Max[1.07136897200911339430148373362126412878`25.882473602983254*Abs[$Failed], 3.9184929141413293297779223548192302498`26.055977321861217*Abs[$Failed], 4.72193174048709588457817725045731890607`26.995203365641412*Abs[$Failed], 17.56276248893889143808911349296873286345`27.004572190927174*Abs[$Failed], 87.71713486749567607533217167885093089404`27.091576614056315*Abs[$Failed]]}, {2, 2, Max[1.07136897200911339430148373369919095528`25.882473602983254*Abs[$Failed], 3.9184929141413293297779223548637154006`26.055977321861217*Abs[$Failed], 4.72193174048709588457817725053745910001`26.995203365641412*Abs[$Failed], 17.56276248893889143808911349308505425475`27.004572190927174*Abs[$Failed], 87.71713486749567607533217167898949857356`27.091576614056315*Abs[$Failed]]}, {2, 3, Max[Abs[(-796*$Failed)/35 - $Failed^2], Abs[(-508*$Failed)/105 - $Failed^2], Abs[(-4*$Failed)/105 - $Failed^2], Abs[(128*$Failed)/35 - $Failed^2], Abs[(4484*$Failed)/35 - $Failed^2]]}, {3, 0, Max[0.18272262799974260078802688388645898498`24.246455500810352*Abs[$Failed], 2.10510417831689454761806996172253672989`25.00095793466144*Abs[$Failed], 2.84613067139669364019749889780712907008`25.749126878671404*Abs[$Failed], 14.32970665410337962418921199711092131687`26.048700676175102*Abs[$Failed], 81.98041253357531154875966767701397382041`26.306515587412985*Abs[$Failed]]}, {3, 1, Max[1.90428008430996684685450909230473021751`25.31141672460716*Abs[$Failed], 3.56796636903654983789399811591496284233`25.34328677424351*Abs[$Failed], 4.5863972252874866250733986427617781505`26.501772726694753*Abs[$Failed], 18.25529957134740526311981651411412055791`26.221242051310963*Abs[$Failed], 91.6363214598149768671997306316583091487`26.281634051172542*Abs[$Failed]]}, {3, 2, Max[0.37569895461278298595938375507671349808`24.29745156185482*Abs[$Failed], 5.53484737485136721352000220001708561024`25.22498806818569*Abs[$Failed], 6.48047935737516031496803346698872717251`25.926543738770572*Abs[$Failed], 21.9102144249475665113027544240029241454`26.032666174289357*Abs[$Failed], 97.88401690519488790562710017899262875894`26.15507415316924*Abs[$Failed]]}, {4, 0, Max[0.24993237035966474639421609391631228344`23.173989954654882*Abs[$Failed], 4.22451471348334816261600936225852457357`24.09862825468889*Abs[$Failed], 5.41333314927960882524631918353992737277`24.870840913895993*Abs[$Failed], 21.695343614477851240202069514621
[output truncated — 4214 chars total]
```
## 2026-07-15 15:57:09 — ✏️ EDIT CELL 39
(* Psi[λ1,λ2,λ3][M,n]: on-demand accessor over the memoized TauEigensystem cache. Returns the
   fixed (last-component-normalized) eigenvector for a genuine state.
   τ[λ1,λ2,λ3][M,n][u]: the FULL τ e…


## 2026-07-15 15:57:09
**Input:**
```wolfram
(* Psi[λ1,λ2,λ3][M,n]: on-demand accessor over the memoized TauEigensystem cache. Returns the
   fixed (last-component-normalized) eigenvector for a genuine state.
   τ[λ1,λ2,λ3][M,n][u]: the FULL τ eigenvalue polynomial for that state, τ0[λ1,λ2,λ3,M][u] +
   J1_n u + J0_n (generalizing L2's τ0[...][u]+J_n by the extra J1 term — see the plan's
   derivation note). Bounds/failure modes match the L2 accessor style (bnds/spurious).
   NOTE: the sanity bound below deliberately omits "n<=M" (unlike L2, where it was always true
   since L2's dPred=Min[M,λ1,λ2,λ1+λ2-M]+1<=M+1 by construction — M is one of Min's own
   arguments). L3's triple-convolution dPred has no such bound (e.g. dPred[1,1,1,1]=3>M+1=2), so
   n<=M would wrongly reject genuine states before the real multiplicity check (n>mult-1, below)
   ever runs. Only n>=0 is checked here; the multiplicity check is the sole source of truth. *)
ClearAll[Psi, τ];
Psi::bnds = "Invalid indices (M,n)=(`1`,`2`) for (λ1,λ2,λ3)=(`3`,`4`,`5`): need M an integer in 0..λ1+λ2+λ3 and n a non-negative integer.";
τ::bnds = "Invalid indices (M,n)=(`1`,`2`) for (λ1,λ2,λ3)=(`3`,`4`,`5`): need M an integer in 0..λ1+λ2+λ3 and n a non-negative integer.";
Psi::spurious = "No genuine state at (M,n)=(`1`,`2`) for (λ1,λ2,λ3)=(`3`,`4`,`5`): weight sector M=`1` has only `6` genuine state(s) (n=0..`7`), so n=`2` does not exist.";
τ::spurious = "No genuine state at (M,n)=(`1`,`2`) for (λ1,λ2,λ3)=(`3`,`4`,`5`): weight sector M=`1` has only `6` genuine state(s) (n=0..`7`), so n=`2` does not exist.";
Psi[λ1_, λ2_, λ3_][M_, n_] := Module[{mult},
  If[! (IntegerQ[M] && IntegerQ[n] && 0 <= M <= λ1 + λ2 + λ3 && n >= 0),
    Message[Psi::bnds, M, n, λ1, λ2, λ3]; Return[$Failed]];
  mult = dPred[M, λ1, λ2, λ3];
  If[n > mult - 1,
    Message[Psi::spurious, M, n, λ1, λ2, λ3, mult, mult - 1]; Return[Missing["NoState", {λ1, λ2, λ3, M, n}]]];
  TauEigensystem[λ1, λ2, λ3][{M, n}, "Psi"]];
τ[λ1_, λ2_, λ3_][M_, n_][u_] := Module[{mult, rec},
  If[! (IntegerQ[M] && IntegerQ[n] && 0 <= M <= λ1 + λ2 + λ3 && n >= 0),
    Message[τ::bnds, M, n, λ1, λ2, λ3]; Return[$Failed]];
  mult = dPred[M, λ1, λ2, λ3];
  If[n > mult - 1,
    Message[τ::spurious, M, n, λ1, λ2, λ3, mult, mult - 1]; Return[Missing["NoState", {λ1, λ2, λ3, M, n}]]];
  rec = TauEigensystem[λ1, λ2, λ3][{M, n}];
  τ0[λ1, λ2, λ3, M][u] + rec["J1"] u + rec["J0"]];
Print["spurious-guard test (expect a Psi::spurious message and Missing[...]): ", Psi[2, 1, 1][3, 3]];
{Psi[1, 1, 1][1, 0], N[τ[1, 1, 1][1, 0][2], 10]}
```
**Output:**
```
[message] Psi::spurious: No genuine state at (M,n)=(3,3) for (λ1,λ2,λ3)=(2,1,1): weight sector M=3 has only 3 genuine state(s) (n=0..2), so n=3 does not exist.
Out= {{0.36043987514722653726667765059831278582`26.85798484313863 + 0.93278244859337471648483598206791350819`26.85798484313863*I, 0.56833118551534074350582873259901777571`26.85798484313863 + 0.5516088635918809743335432028014154242`26.701746208530505*I, 0.43156394084289798758734281499217795101`26.85798484313863 + 0.39288237789502186246985745181354135617`26.688077093908458*I, 0.88096705226680694396357659277617773607`26.85798484313863 + 0.3405301294515432284073298252721106832`26.463575570362078*I, 0.25993868083405777770651629514377446006`26.698428155925086 - 0.2720672628658623596684205879497517111`26.85798484313863*I, 0.19885311087688781514666096669532055364`26.85798484313863 + 0.26094479422801754403724338682337866528`26.760176866170834*I, 0.00149946229798485268979539231628574037`24.688214897154637 + 0.33130752491381394404026276997360157274`26.85798484313863*I, 1.00000000000000000000000000000000000001`26.85798484313863 + 0``27.030341651165482*I}, 2.45560487187140974687950049959647506927`10.10603276460253 - 1.17083386404855165681638910128993356489`9.784369553312159*I}
```
## 2026-07-15 15:57:20 — ✏️ EDIT CELL 40
(*QSolve[λ1,λ2,λ3][M,n]: linear solve, generalizing L2. τ(u)=τ0(u)+J1 u+J0 now has TWO free
  state-dependent parameters (not L2's single J), so the known-J shift subtracted from
  BaxterMatrix is TWO…


## 2026-07-15 15:57:20
**Input:**
```wolfram
(*QSolve[λ1,λ2,λ3][M,n]: linear solve, generalizing L2. τ(u)=τ0(u)+J1 u+J0 now has TWO free
  state-dependent parameters (not L2's single J), so the known-J shift subtracted from
  BaxterMatrix is TWO terms: J0·Identity (same diagonal-shift role as L2's J·Identity) PLUS
  J1·shiftDown, where shiftDown is the (M+1)x(M+1) matrix with 1's on the subdiagonal (row=col+1,
  1-indexed) — this represents multiplying the column's underlying u^k monomial by an extra u^1,
  i.e. "J1 u^(k+1)" landing one row below the diagonal relative to column k. Same
  Chop-before-NullSpace precision-tracked-zero fix as L2 (see L2's QSolve comment for the full
  numerical-precision rationale — unchanged here).
  NOTE: same n<=M bound removal as Psi/τ above (Cell 39) — L3's dPred can exceed M+1, so only
  n>=0 is sanity-checked here; the multiplicity check (n>mult-1, below) is authoritative.*)
ClearAll[QSolve];
QSolve::bnds = "Invalid indices (M,n)=(`1`,`2`) for (λ1,λ2,λ3)=(`3`,`4`,`5`): need M an integer in 0..λ1+λ2+λ3 and n a non-negative integer.";
QSolve::spurious = "No genuine state at (M,n)=(`1`,`2`) for (λ1,λ2,λ3)=(`3`,`4`,`5`): weight sector M=`1` has only `6` genuine state(s) (n=0..`7`), so n=`2` does not exist.";
QSolve::degnull = "Unexpected Q-function nullspace dimension `6` (expected 1) for (M,n)=(`1`,`2`) at (λ1,λ2,λ3)=(`3`,`4`,`5`).";
QSolve::lastzero = "Q-function nullspace vector found but its top (u^M) coefficient is ~0 for (M,n)=(`1`,`2`) at (λ1,λ2,λ3)=(`3`,`4`,`5`); cannot normalize to monic.";
QSolve[λ1_, λ2_, λ3_][M_, n_] := QSolve[λ1, λ2, λ3][M, n] = Module[{mult, rec, J1n, J0n, shiftDown, mat, ns, vec, last},
  If[! (IntegerQ[M] && IntegerQ[n] && 0 <= M <= λ1 + λ2 + λ3 && n >= 0),
    Message[QSolve::bnds, M, n, λ1, λ2, λ3]; Return[$Failed]];
  mult = dPred[M, λ1, λ2, λ3];
  If[n > mult - 1,
    Message[QSolve::spurious, M, n, λ1, λ2, λ3, mult, mult - 1]; Return[$Failed]];
  rec = TauEigensystem[λ1, λ2, λ3][{M, n}];
  J1n = rec["J1"]; J0n = rec["J0"];
  shiftDown = Table[If[i == j + 1, 1, 0], {i, 1, M + 1}, {j, 1, M + 1}];
  mat = BaxterMatrix[λ1, λ2, λ3, M] - J1n shiftDown - J0n IdentityMatrix[M + 1];
  ns = NullSpace[Chop[mat, 10^-12], Tolerance -> 10^-10];
  Which[
    Length[ns] != 1, Message[QSolve::degnull, M, n, λ1, λ2, λ3, Length[ns]]; Return[$Failed],
    True,
      vec = ns[[1]]; last = vec[[-1]];
      If[Abs[last] < 10^-15 Max[Abs[vec]],
        Message[QSolve::lastzero, M, n, λ1, λ2, λ3]; Return[$Failed]];
      vec/last]
];
QSolve[1, 1, 1][2, 0]
```
**Output:**
```
Out= {-0.0499573806412398679366027280112104819`25.36115104115277 - 1.23005323914750422415009675314876879381`25.36115104115277*I, -1.85605135792996542485027196616457263772`25.338198332697328 + 1.31753882642786619162173204602080684677`25.361151041152763*I, 0.99999999999999999999999999999999999999`25.36115104115277 + 0``26.45216107225751*I}
```
## 2026-07-15 15:57:31 — ▶️ RUN CELL 42
**In [42]:** `(*End-to-end Baxter TQ-equation residual check: a(u)Q(u+h) - τ(u)Q(u) + Qθ(u)Q(u-h) = 0, assembled↵  purely from Q1 and τ (which trace back to different derivations — Q1 through↵  QSolve/BaxterMatrix,`
**Out:** (no output)

## 2026-07-15 15:57:32 — ▶️ RUN CELL 43
**In [43]:** `Table[{M, n, checkTQ[2, 1, 1, M, n]}, {M, 0, 4}, {n, 0, dPred[M, 2, 1, 1] - 1}] // Flatten[#, 1] &`
**Out:** (no output)


## 2026-07-15 15:57:38
**Input:**
```wolfram
Table[{M, n, checkTQ[1, 1, 1, M, n]}, {M, 0, 3}, {n, 0, dPred[M, 1, 1, 1] - 1}] // Flatten[#, 1] &
```
**Output:**
```
Out= {{0, 0, Max[0.15025893010484758828777556108151636731`23.272590612439533*Abs[$Failed], 3.08685408200609875661540311241351508555`25.03283327912714*Abs[$Failed], 4.15595077293946712784708609957381655942`24.446929235991636*Abs[$Failed], 14.99741316057572636787500207978855534787`25.222431770691028*Abs[$Failed], 101.19194196892583915656251845445793020112`25.5014130561257*Abs[$Failed]]}, {1, 0, Max[0.08970761353168929728508399825601280999`24.38894365795779*Abs[$Failed], 1.19887312832888762397195151837251507484`25.86723582005619*Abs[$Failed], 1.8044082417844273274236113674411597349`25.38176890472557*Abs[$Failed], 9.41850807735837046095167039960057832429`26.353788976904998*Abs[$Failed], 87.11643513046212518542637446186159793877`26.816860496054232*Abs[$Failed]]}, {1, 1, Max[1.10538837461520941511259463535868899354`24.90204728677379*Abs[$Failed], 2.71085466907731657534075469121920494132`26.185464450644506*Abs[$Failed], 3.54412540980257331556193206020528052904`25.301098271000463*Abs[$Failed], 12.86060007519508084426694230102000890504`26.04177475288251*Abs[$Failed], 95.02452393050365573611331298348941638396`26.288751358892863*Abs[$Failed]]}, {1, 2, Max[0.28267846552647238876748877769662615753`24.899909934378154*Abs[$Failed], 3.94568059474715988219053320476702900435`26.55249615126905*Abs[$Failed], 4.92925838778893003509938778840263786976`25.996256002344104*Abs[$Failed], 15.35678488209999983679409380897595380947`26.71106340333319*Abs[$Failed], 99.90377204980664677586142977168680744096`26.997636871582728*Abs[$Failed]]}, {2, 0, 0``21.46640758999653}, {2, 1, Max[1.10538837461520941511259463440195895185`24.902041722424332*Abs[$Failed], 2.71085466907731657534075468544601596307`26.185464450644506*Abs[$Failed], 3.54412540980257331556193205424688233991`25.30109156439836*Abs[$Failed], 12.86060007519508084426694229252909473494`26.04177115707164*Abs[$Failed], 95.02452393050365573611331295796398789056`26.288731187623323*Abs[$Failed]]}, {2, 2, Max[0.28267846552647238876748877674588071913`24.89988879683595*Abs[$Failed], 3.94568059474715988219053319798496133716`26.55249615126905*Abs[$Failed], 4.92925838778893003509938777601015896818`25.996218001087982*Abs[$Failed], 15.35678488209999983679409379055492276616`26.711048219131687*Abs[$Failed], 99.90377204980664677586142972924850520088`26.997537422879684*Abs[$Failed]]}, {3, 0, Max[0.15025893010484758828777556104549722531`23.27258764404846*Abs[$Failed], 3.08685408200609875661540311258339575725`25.03283327912714*Abs[$Failed], 4.15595077293946712784708609940918470285`24.446923035717575*Abs[$Failed], 14.99741316057572636787500208000431011184`25.22242935656539*Abs[$Failed], 101.1919419689258391565625184540756413668`25.50139906004645*Abs[$Failed]]}}
```

## 2026-07-15 15:57:38
**Input:**
```wolfram
Table[{M, n, checkTQ[2, 1, 1, M, n]}, {M, 0, 4}, {n, 0, dPred[M, 2, 1, 1] - 1}] // Flatten[#, 1] &
```
**Output:**
```
Out= {{0, 0, Max[0.24993237035966474639421609395593025604`23.17399178293198*Abs[$Failed], 4.22451471348334816261600936233818115148`24.098631940571227*Abs[$Failed], 5.41333314927960882524631918347377411112`24.870840913895993*Abs[$Failed], 21.69534361447785124020206951445664389383`25.021294615572344*Abs[$Failed], 100.44467915995734938001804011454749304936`25.144188555382186*Abs[$Failed]]}, {1, 0, Max[0.18272262799974260078802688336253484185`24.24646999719398*Abs[$Failed], 2.1051041783168945476180699610224036158`25.00098790754616*Abs[$Failed], 2.84613067139669364019749889859623123615`25.749126878671404*Abs[$Failed], 14.32970665410337962418921199929922836129`26.048710876256283*Abs[$Failed], 81.98041253357531154875966767934335288941`26.306581245947868*Abs[$Failed]]}, {1, 1, Max[1.90428008430996684685450909223371044779`25.311433114716067*Abs[$Failed], 3.56796636903654983789399811601174919346`25.34330107734096*Abs[$Failed], 4.58639722528748662507339864279762205798`26.501772726694753*Abs[$Failed], 18.2552995713474052631198165142109365891`26.22125111835505*Abs[$Failed], 91.63632145981497686719973063200703667053`26.28168614359125*Abs[$Failed]]}, {1, 2, Max[0.37569895461278298595938375519014956033`24.29745987906615*Abs[$Failed], 5.53484737485136721352000220049820197874`25.225001216281505*Abs[$Failed], 6.48047935737516031496803346711365812404`25.926543738770572*Abs[$Failed], 21.91021442494756651130275442429489771226`26.032671504982815*Abs[$Failed], 97.88401690519488790562710018042829866446`26.15511040488361*Abs[$Failed]]}, {2, 0, Max[0.02051985356943702835983130688291796136`23.896297413619525*Abs[$Failed], 0.54466138157369124077229084326220410751`25.114330691194443*Abs[$Failed], 1.29226275045467693251162386208438469128`26.099584227107243*Abs[$Failed], 10.94528880731531002448295501519559498672`26.623350381857655*Abs[$Failed], 74.74805009888322176604931679189349139259`26.92150110456903*Abs[$Failed]]}, {2, 1, Max[1.07136897200911339430148373362126412878`25.882473602983254*Abs[$Failed], 3.9184929141413293297779223548192302498`26.055977321861217*Abs[$Failed], 4.72193174048709588457817725045731890607`26.995203365641412*Abs[$Failed], 17.56276248893889143808911349296873286345`27.004572190927174*Abs[$Failed], 87.71713486749567607533217167885093089404`27.091576614056315*Abs[$Failed]]}, {2, 2, Max[1.07136897200911339430148373369919095528`25.882473602983254*Abs[$Failed], 3.9184929141413293297779223548637154006`26.055977321861217*Abs[$Failed], 4.72193174048709588457817725053745910001`26.995203365641412*Abs[$Failed], 17.56276248893889143808911349308505425475`27.004572190927174*Abs[$Failed], 87.71713486749567607533217167898949857356`27.091576614056315*Abs[$Failed]]}, {2, 3, Max[0.01559154381713799665645209436111770071`24.110419980007524*Abs[$Failed], 6.30146769355364662213652205273154825842`26.567140408914007*Abs[$Failed], 6.976846267661482263843288269065657518`26.990219703370926*Abs[$Failed], 22.35056723911549571216256723040217630227`27.26681941009903*Abs[$Failed], 97.66694115464331821645739142603876100997`27.48057675249778*Abs[$Failed]]}, {3, 0, Max[0.18272262799974260078802688388645898498`24.246455500810352*Abs[$Failed], 2.10510417831689454761806996172253672989`25.00095793466144*Abs[$Failed], 2.84613067139669364019749889780712907008`25.749126878671404*Abs[$Failed], 14.32970665410337962418921199711092131687`26.048700676175102*Abs[$Failed], 81.98041253357531154875966767701397382041`26.306515587412985*Abs[$Failed]]}, {3, 1, Max[1.90428008430996684685450909230473021751`25.31141672460716*Abs[$Failed], 3.56796636903654983789399811591496284233`25.34328677424351*Abs[$Failed], 4.5863972252874866250733986427617781505`26.501772726694753*Abs[$Failed], 18.25529957134740526311981651411412055791`26.221242051310963*Abs[$Failed], 91.6363214598149768671997306316583091487`26.281634051172542*Abs[$Failed]]}, {3, 2, Max[0.37569895461278298595938375507671349808`24.29745156185482*Abs[$Failed], 5.53484737485136721352000220001708561024`25.22498806818569*Abs[$Failed], 6.48047935737516031496803346698872717251`25.926543738770572*Abs[$Failed], 21.91021442494756651130275442400292
[output truncated — 4594 chars total]
```
## 2026-07-15 15:59:16 — ▶️ RUN CELL 39
**In [39]:** `(* Psi[λ1,λ2,λ3][M,n]: on-demand accessor over the memoized TauEigensystem cache. Returns the↵   fixed (last-component-normalized) eigenvector for a genuine state.↵   τ[λ1,λ2,λ3][M,n][u]: the FULL τ e`
**Out:** (no output)

## 2026-07-15 15:59:19 — ▶️ RUN CELL 40
**In [40]:** `(*QSolve[λ1,λ2,λ3][M,n]: linear solve, generalizing L2. τ(u)=τ0(u)+J1 u+J0 now has TWO free↵  state-dependent parameters (not L2's single J), so the known-J shift subtracted from↵  BaxterMatrix is TWO`
**Out:** (no output)

## 2026-07-15 15:59:22 — ▶️ RUN CELL 41
**In [41]:** `(* Q1[λ1,λ2,λ3][M,n][u] = κ1^(u/h) q_n(u): the actual Baxter Q-function value, built from the↵   SAME QSolve cache (mutually consistent with τ by construction). *)↵ClearAll[Q1];↵Q1[λ1_, λ2_, λ3_][M_, `
**Out:** (no output)

## 2026-07-15 15:59:25 — ▶️ RUN CELL 42
**In [42]:** `(*End-to-end Baxter TQ-equation residual check: a(u)Q(u+h) - τ(u)Q(u) + Qθ(u)Q(u-h) = 0, assembled↵  purely from Q1 and τ (which trace back to different derivations — Q1 through↵  QSolve/BaxterMatrix,`
**Out:** (no output)

## 2026-07-15 16:06:51 — ✏️ EDIT CELL 42
(*End-to-end Baxter TQ-equation residual check: a(u)Q(u+h) - τ(u)Q(u) + Qθ(u)Q(u-h) = 0, assembled
  purely from Q1 and τ (which trace back to different derivations — Q1 through
  QSolve/BaxterMatrix,…

## 2026-07-15 16:22:58 — ✏️ EDIT CELL 32
(*NOTE: as in Task 4's commutativity check, the symbolic Simplify//Flatten//Union above does not
  collapse to the literal list {0} — same benign precision-tracked-zero display artifact (exact
  ratio…


## 2026-07-15 16:22:58
**Input:**
```wolfram
(*NOTE: as in Task 4's commutativity check, the symbolic Simplify//Flatten//Union above does not
  collapse to the literal list {0} — same benign precision-tracked-zero display artifact (exact
  rational polynomial pieces mixed with 30-digit κ1,κ2 numbers). Confirms the SAME residual
  vanishes exactly (to Chop tolerance) at a generic complex numeric u, for every (λ1,λ2,λ3,s1,s2,s3)
  in the same range — this is the substantive check, closing Section A.*)
Max[Abs[Flatten[Table[Table[x3[λ1,λ2,λ3][{s1,s2,s3}].B3[λ1,λ2,λ3][u]-(u-θ1-h s1)(u-θ2-h s2)(u-θ3-h s3)x3[λ1,λ2,λ3][{s1,s2,s3}],{s1,0,λ1},{s2,0,λ2},{s3,0,λ3}],{λ1,1,2},{λ2,1,2},{λ3,1,2}] /. u -> N[3.7,30]+I N[1.2,30]]]]//Chop
```
**Output:**
```
Out= 0
```
## 2026-07-15 16:23:02 — ▶️ RUN CELL 18 (0.03 s)
**In [18]:** `(*NOTE: the symbolic Simplify//Flatten//Union above does not collapse to the literal list {0} —↵  the residual terms shown are combinations of exact rational coefficients with 30-digit-precision↵  κ1,`
**Out:** Out[211]= 0

## 2026-07-15 16:23:32 — ▶️ RUN CELL 43 (0.09 s)
**In [43]:** `(*End-to-end Baxter TQ-equation residual check: a(u)Q(u+h) - τ(u)Q(u) + Qθ(u)Q(u-h) = 0, assembled↵  purely from Q1 and τ (which trace back to different derivations — Q1 through↵  QSolve/BaxterMatrix,`
**Out:** Out[213]= \begin{pmatrix}0 & 0 & 0.\times 10^{-23}\\1 & 0 & 0.\times 10^{-22}\\1 & 1 & 0.\times 10^{-22}\\1 & 2 & 0.\times 10^{-22}\\2 & 0 & 0.\times 10^{-21}\\2 & 1 & 0.\times 10^{-21}\\2 & 2 & 0.\times 10^{-21}\\3 & 0 & 0.\times 10^{-19}\end{pmatrix}

## 2026-07-15 16:23:36 — ▶️ RUN CELL 44 (0.03 s)
**In [44]:** `Table[{M, n, checkTQ[2, 1, 1, M, n]}, {M, 0, 4}, {n, 0, dPred[M, 2, 1, 1] - 1}] // Flatten[#, 1] &`
**Out:** Out[214]= \begin{pmatrix}0 & 0 & 0.\times 10^{-23}\\1 & 0 & 0.\times 10^{-22}\\1 & 1 & 0.\times 10^{-22}\\1 & 2 & 0.\times 10^{-22}\\2 & 0 & 0.\times 10^{-22}\\2 & 1 & 0.\times 10^{-22}\\2 & 2 & 0.\times 10^{-22}\\2 & 3 & 0.\times 10^{-22}\\3 & 0 & 0.\times 10^{-20}\\3 & 1 & 0.\times 10^{-20}\\3 & 2

## 2026-07-15 16:37:50 — ▶️ RUN CELL 32 (0.03 s)
**In [32]:** `(*NOTE: as in Task 4's commutativity check, the symbolic Simplify//Flatten//Union above does not↵  collapse to the literal list {0} — same benign precision-tracked-zero display artifact (exact↵  ratio`
**Out:** Out[215]= 0

