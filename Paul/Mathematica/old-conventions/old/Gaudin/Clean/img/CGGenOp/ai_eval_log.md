# AI Action Log — CGGenOp

## 2026-07-25 14:39:09 — 📥 BULK INSERT 4 CELL(S) at positions 1–4
1. [markdown] ## Weight decomposition of $\langle HWS| e^{a_{13}E_{13}} e^{a_{14}E_{14}} e^{a_{23}E_{23}} e^{a_{24
2. [markdown] ### Weight grading↵↵A covector $\langle HWS|E_{13}$ pairs nontrivially only with kets of weight $\la
3. [code] (* --- gl(4) generators in the Gelfand-Tsetlin basis (ported verbatim from Hamiltonians_SU4.wb) --- 
4. [code] (* Sanity: the four generators pairwise commute, and each annihilates |HWS> (they are raising ops).


## 2026-07-25 14:39:09
**Input:**
```wolfram
(* --- gl(4) generators in the Gelfand-Tsetlin basis (ported verbatim from Hamiltonians_SU4.wb) --- *)
ps[rep_] := ps[rep] = Block[{res = (Prepend[rep, #] & /@ (Tuples[Table[Range[rep[[1, i + 1]], rep[[1, i]]], {i, 1, Length[rep[[1]]] - 1}]]))}, If[Length[res[[1, 1]]] == 1, Reverse@Sort@res, Reverse@Sort@Flatten[ps /@ res, 1]]];
dim[rep_] := dim[rep] = Length[ps[{rep}]];
J0[k_][pat_] := J0[k][pat] = Sum[pat[[k, i]], {i, 1, k}] - Sum[pat[[k - 1, i]], {i, 1, k - 1}];
Jp[k_][pat_] := Jp[k][pat] = -Sum[Product[pat[[k, i]] - i - pat[[k + 1, j]] + j, {j, Range[k + 1]}]/Product[pat[[k, i]] - i - pat[[k, j]] + j, {j, Complement[Range[k], {i}]}] v[ReplacePart[pat, {k, i} -> pat[[k, i]] + 1]], {i, 1, k}];
Jm[k_][pat_] := Jm[k][pat] = Sum[Product[pat[[k, i]] - i - pat[[k - 1, j]] + j, {j, Range[k - 1]}]/Product[pat[[k, i]] - i - pat[[k, j]] + j, {j, Complement[Range[k], {i}]}] v[ReplacePart[pat, {k, i} -> pat[[k, i]] - 1]], {i, 1, k}];
JJ[k_, l_, rep_] := JJ[k, l, rep] = Which[k == l, DiagonalMatrix[J0[k] /@ ps[{rep}]], l == k + 1, Transpose[Table[Coefficient[Jp[k][ps[{rep}][[i]]], v[#]] & /@ ps[{rep}], {i, 1, Length[ps[{rep}]]}]], l == k - 1, Transpose[Table[Coefficient[Jm[l][ps[{rep}][[i]]], v[#]] & /@ ps[{rep}], {i, 1, Length[ps[{rep}]]}]], l > k + 1, JJ[k, k + 1, rep] . JJ[k + 1, l, rep] - JJ[k + 1, l, rep] . JJ[k, k + 1, rep], l < k - 1, JJ[k, k - 1, rep] . JJ[k - 1, l, rep] - JJ[k - 1, l, rep] . JJ[k, k - 1, rep]];
EE[rep_][i_, j_] := JJ[i, j, rep];
```
**Output:**
```
(no output)
```

## 2026-07-25 14:39:10
**Input:**
```wolfram
(* Sanity: the four generators pairwise commute, and each annihilates |HWS> (they are raising ops). *)
testrep = {3, 2, 1, 0};
idx = {{1, 3}, {1, 4}, {2, 3}, {2, 4}};
comms = Table[EE[testrep] @@ p . EE[testrep] @@ q - EE[testrep] @@ q . EE[testrep] @@ p, {p, idx}, {q, idx}];
hws = UnitVector[dim[testrep], 1];  (* ps[{rep}] is reverse-sorted, so pattern 1 is the HWS *)
{Max[Abs[Flatten[comms]]], Max[Abs[Flatten[Table[EE[testrep] @@ p . hws, {p, idx}]]]], J0[#][ps[{testrep}][[1]]] & /@ Range[4]}
```
**Output:**
```
(no output)
```

## 2026-07-25 14:39:17
**Input:**
```wolfram
testrep = {3, 2, 1, 0}; idx = {{1, 3}, {1, 4}, {2, 3}, {2, 4}}; comms = Table[EE[testrep] @@ p . EE[testrep] @@ q - EE[testrep] @@ q . EE[testrep] @@ p, {p, idx}, {q, idx}]; hws = UnitVector[dim[testrep], 1]; {dim[testrep], Max[Abs[Flatten[comms]]], Max[Abs[Flatten[Table[EE[testrep] @@ p . hws, {p, idx}]]]], J0[#][ps[{testrep}][[1]]] & /@ Range[4]}
```
**Output:**
```
Out= {64, 0, 0, {3, 2, 1, 0}}
```
## 2026-07-25 14:39:35 — 📥 BULK INSERT 2 CELL(S) at positions 5–6
1. [markdown] ### The grading routines↵↵`nList[d]` enumerates exponent tuples $n=(n_{13},n_{14},n_{23},n_{24})$ of
2. [code] (* Exponent tuples of total degree d, the weight they carry, and the multinomial coefficient. *)↵nLi


## 2026-07-25 14:39:35
**Input:**
```wolfram
(* Exponent tuples of total degree d, the weight they carry, and the multinomial coefficient. *)
nList[d_] := Flatten[Table[{n13, n14, n23, d - n13 - n14 - n23}, {n13, 0, d}, {n14, 0, d - n13}, {n23, 0, d - n13 - n14}], 2];
wt[lam_, {n13_, n14_, n23_, n24_}] := lam + {-(n13 + n14), -(n23 + n24), n13 + n23, n14 + n24};
coef[{n13_, n14_, n23_, n24_}] := a[1, 3]^n13 a[1, 4]^n14 a[2, 3]^n23 a[2, 4]^n24/(n13! n14! n23! n24!);

(* Weight sectors at fixed total degree d: association weight -> list of exponent tuples. *)
sectors[lam_, d_] := GroupBy[nList[d], wt[lam, #] &];

(* Sector labelled directly by the margins m = {m1,m2,m3,m4} (must satisfy m1+m2==m3+m4). *)
marginSector[{m1_, m2_, m3_, m4_}] := If[m1 + m2 =!= m3 + m4, {}, Table[{t, m1 - t, m3 - t, m2 - m3 + t}, {t, Max[0, m3 - m2], Min[m1, m3]}];
```
**Output:**
```
(no output)
```
## 2026-07-25 14:39:50 — ✏️ EDIT CELL 6
(* Exponent tuples of total degree d, the weight they carry, and the multinomial coefficient. *)
nList[d_] := Flatten[Table[{n13, n14, n23, d - n13 - n14 - n23}, {n13, 0, d}, {n14, 0, d - n13}, {n23, …


## 2026-07-25 14:39:50
**Input:**
```wolfram
(* Exponent tuples of total degree d, the weight they carry, and the multinomial coefficient. *)
nList[d_] := Flatten[Table[{n13, n14, n23, d - n13 - n14 - n23}, {n13, 0, d}, {n14, 0, d - n13}, {n23, 0, d - n13 - n14}], 2];
wt[lam_, {n13_, n14_, n23_, n24_}] := lam + {-(n13 + n14), -(n23 + n24), n13 + n23, n14 + n24};
coef[{n13_, n14_, n23_, n24_}] := a[1, 3]^n13 a[1, 4]^n14 a[2, 3]^n23 a[2, 4]^n24/(n13! n14! n23! n24!);

(* Weight sectors at fixed total degree d: association weight -> list of exponent tuples. *)
sectors[lam_, d_] := GroupBy[nList[d], wt[lam, #] &];

(* Same sector labelled directly by the margins m={m1,m2,m3,m4} (needs m1+m2==m3+m4);
   the free parameter t=n13 runs over min(m)+1 values. *)
marginSector[{m1_, m2_, m3_, m4_}] := If[m1 + m2 =!= m3 + m4, {},
   Table[{t, m1 - t, m3 - t, m2 - m3 + t}, {t, Max[0, m3 - m2], Min[m1, m3]}]];

(* Symbolic bra expansion truncated at total degree dmax, organised by weight sector:
   each entry is weight -> sum over that sector of coef[n] * bra[n]. *)
braExp[lam_, dmax_] := Association @@ Flatten[Table[
    Normal@Map[Total[coef[#] bra @@ # & /@ #] &, sectors[lam, d]] , {d, 0, dmax}]];
```
**Output:**
```
Out= Null
```
## 2026-07-25 14:40:09 — 📥 BULK INSERT 4 CELL(S) at positions 7–10
1. [markdown] ### Numerical verification in an explicit irrep↵↵Three things to check in $\lambda=\{3,2,1,0\}$ (dim
2. [code] (* (i) The four exponentials really do merge (nilpotent => MatrixExp is polynomial). *)↵mex[m_] := M
3. [code] (* (ii) Grading check: weights of the GT patterns supporting each bra monomial. *)↵patWt[pat_] := J0
4. [code] (* (iii) Sector sizes: every weight sector has exactly min(m1,m2,m3,m4)+1 monomials,↵   and grouping


## 2026-07-25 14:45:09
**Input:**
```wolfram
(* (i) The four exponentials really do merge (nilpotent => MatrixExp is polynomial). *)
mex[m_] := MatrixExp[m];
lhs = mex[a13 EE[testrep][1, 3]] . mex[a14 EE[testrep][1, 4]] . mex[a23 EE[testrep][2, 3]] . mex[a24 EE[testrep][2, 4]];
rhs = mex[a13 EE[testrep][1, 3] + a14 EE[testrep][1, 4] + a23 EE[testrep][2, 3] + a24 EE[testrep][2, 4]];
Max[Abs[Expand[lhs - rhs]]] /. {a13 -> 7/5, a14 -> -3/2, a23 -> 11/7, a24 -> 2/9}
```
**Output:**
```
(no output)
```
## 2026-07-25 14:45:19 — ✏️ EDIT CELL 8
(* (i) The four exponentials really do merge into one (nilpotent => MatrixExp is polynomial).
   Checked at rational parameter values to keep the 64x64 exponentials cheap. *)
avals = {7/5, -3/2, 11/7,…


## 2026-07-25 14:45:22
**Input:**
```wolfram
(* (i) The four exponentials really do merge into one (nilpotent => MatrixExp is polynomial).
   Checked at rational parameter values to keep the 64x64 exponentials cheap. *)
avals = {7/5, -3/2, 11/7, 2/9};
gens = EE[testrep] @@@ {{1, 3}, {1, 4}, {2, 3}, {2, 4}};
lhs = Dot @@ (MatrixExp[#1 #2] & @@@ Transpose[{avals, gens}]);
rhs = MatrixExp[avals . gens];
Max[Abs[lhs - rhs]]
```
**Output:**
```
Out= 0
```
## 2026-07-25 14:45:28 — ▶️ RUN CELLS 9–10
Cell 9: ✓ — (no output)
Cell 10: ✓ — (no output)


## 2026-07-25 14:45:39
**Input:**
```wolfram
{And @@ Flatten[Table[gradeCheck[testrep, testrep, n], {d, 0, 6}, {n, nList[d]}]], Table[With[{s = sectors[{l1, l2, l3, l4}, d]}, And[And @@ KeyValueMap[Length[#2] === Min[margins[First[#2]]] + 1 &, s], And @@ (Sort[#] === Sort[marginSector[margins[First[#]]]] & /@ Values[s])]], {d, 0, 8}]}
```
**Output:**
```
[message] MatrixPower::sing: Matrix {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, … [23031 chars]
[message] [MatrixPower::sing ×2 more suppressed — first instance shown above]
[message] General::stop: Further output of MatrixPower::sing will be suppressed during this calculation.
[message] Pick::incomp: Expressions {{3, 2, 1, 0}, {3, 2, 0, 1}, {3, 1, 2, 0}, {3, 1, 1, 1}, {3, 1, 1, 1}, {3, 1, 0, 2}, {3, 0, 2, 1}, {3, 0, 1, 2}, {2, 3, 1, 0}, {2, 3, 0, 1}, {2, 2, 2, 0}, {2, 2, 1, 1}, {2, 2… [101706 chars]
[message] [Pick::incomp ×2 more suppressed — first instance shown above]
[message] [General::stop ×1 more suppressed — first instance shown above]
Out= {False, {True, True, True, True, True, True, True, True, True}}
```
## 2026-07-25 14:45:54 — ✏️ EDIT CELL 9
(* (ii) Grading check: the GT patterns supporting each bra monomial all carry weight nu(n).
   mpow avoids MatrixPower[nilpotent, 0], which throws MatrixPower::sing. *)
mpow[m_, 0] := IdentityMatrix[L…


## 2026-07-25 14:45:59
**Input:**
```wolfram
(* (ii) Grading check: the GT patterns supporting each bra monomial all carry weight nu(n).
   mpow avoids MatrixPower[nilpotent, 0], which throws MatrixPower::sing. *)
mpow[m_, 0] := IdentityMatrix[Length[m]];
mpow[m_, n_] := MatrixPower[m, n];
patWt[pat_] := J0[#][pat] & /@ Range[4];
allWts = patWt /@ ps[{testrep}];
braVec[rep_, {n13_, n14_, n23_, n24_}] := UnitVector[dim[rep], 1] . mpow[EE[rep][1, 3], n13] . mpow[EE[rep][1, 4], n14] . mpow[EE[rep][2, 3], n23] . mpow[EE[rep][2, 4], n24];
gradeCheck[rep_, lam_, n_] := Module[{v = braVec[rep, n], supp},
   supp = DeleteDuplicates[Pick[allWts, Unitize[v], 1]];
   supp === {} || supp === {wt[lam, n]}];
(* nonzero monomials, and whether every one of them is correctly graded *)
{Count[Flatten[Table[braVec[testrep, n] =!= 0 Range[dim[testrep]], {d, 0, 6}, {n, nList[d]}], 1], True],
 And @@ Flatten[Table[gradeCheck[testrep, testrep, n], {d, 0, 6}, {n, nList[d]}]]}
```
**Output:**
```
Out= {34, True}
```
## 2026-07-25 14:46:17 — 📥 BULK INSERT 3 CELL(S) at positions 11–13
1. [markdown] ### Using it↵↵`sectorTable[\lambda, d]` lists, for total degree $d$, each weight $\nu$ together with
2. [code] (* Weight sectors at fixed total degree, with multiplicity and collected coefficient. *)↵sectorTable
3. [code] (* The weight-nu component of <HWS| Exp[Sum a_ij E_ij] in a concrete irrep:↵   sum over the sector o


## 2026-07-25 14:46:17
**Input:**
```wolfram
(* Weight sectors at fixed total degree, with multiplicity and collected coefficient. *)
sectorTable[lam_, d_] := TableForm[
   KeyValueMap[{#1, Length[#2], Total[coef /@ #2]} &, KeySort@sectors[lam, d]],
   TableHeadings -> {None, {"weight nu", "mult", "sum of coef[n]"}}];

sectorTable[{l1, l2, l3, l4}, 4]
```
**Output:**
```
(no output)
```

## 2026-07-25 14:46:18
**Input:**
```wolfram
(* The weight-nu component of <HWS| Exp[Sum a_ij E_ij] in a concrete irrep:
   sum over the sector of coef[n] * (bra monomial). Truncation is automatic --
   monomials outside the rep vanish. *)
sectorVec[rep_, lam_, nu_] := Module[{m = {lam[[1]] - nu[[1]], lam[[2]] - nu[[2]], nu[[3]] - lam[[3]], nu[[4]] - lam[[4]]}},
   Total[coef[#] braVec[rep, #] & /@ marginSector[m]]];

(* Example: two-step-down weight sectors in lambda = {3,2,1,0}. Check each is genuinely
   supported on that weight only, and record which sectors survive in the rep. *)
with2 = Select[Keys@sectors[testrep, 2], AllTrue[sectorVec[testrep, testrep, #], PossibleZeroQ] === False &];
Table[{nu, Min[{testrep[[1]] - nu[[1]], testrep[[2]] - nu[[2]], nu[[3]] - testrep[[3]], nu[[4]] - testrep[[4]]}] + 1,
   DeleteDuplicates[Pick[allWts, Unitize[sectorVec[testrep, testrep, nu] /. a[__] -> 1], 1]] === {nu}}, {nu, with2}] // TableForm
```
**Output:**
```
(no output)
```

## 2026-07-25 14:46:21
**Input:**
```wolfram
sectorTable[{l1, l2, l3, l4}, 3]
```
**Output:**
```
Out=   weight nu   mult   sum of coef[n]
  -3 + l1                   3
  l2                 a[1, 4]
  l3                 --------
  3 + l4      1         6

  -3 + l1                           2
  l2                 a[1, 3] a[1, 4]
  1 + l3             ----------------
  2 + l4      1             2

  -3 + l1                   2
  l2                 a[1, 3]  a[1, 4]
  2 + l3             ----------------
  1 + l4      1             2

  -3 + l1                   3
  l2                 a[1, 3]
  3 + l3             --------
  l4          1         6

  -2 + l1                   2
  -1 + l2            a[1, 4]  a[2, 4]
  l3                 ----------------
  3 + l4      1             2

  -2 + l1                   2
  -1 + l2            a[1, 4]  a[2, 3]
  1 + l3             ---------------- + a[1, 3] a[1, 4] a[2, 4]
  2 + l4      2             2

  -2 + l1                                             2
  -1 + l2                                      a[1, 3]  a[2, 4]
  2 + l3             a[1, 3] a[1, 4] a[2, 3] + ----------------
  1 + l4      2                                       2

  -2 + l1                   2
  -1 + l2            a[1, 3]  a[2, 3]
  3 + l3             ----------------
  l4          1             2

  -1 + l1                           2
  -2 + l2            a[1, 4] a[2, 4]
  l3                 ----------------
  3 + l4      1             2

  -1 + l1                                                     2
  -2 + l2                                      a[1, 3] a[2, 4]
  1 + l3             a[1, 4] a[2, 3] a[2, 4] + ----------------
  2 + l4      2                                       2

  -1 + l1                           2
  -2 + l2            a[1, 4] a[2, 3]
  2 + l3             ---------------- + a[1, 3] a[2, 3] a[2, 4]
  1 + l4      2             2

  -1 + l1                           2
  -2 + l2            a[1, 3] a[2, 3]
  3 + l3             ----------------
  l4          1             2

  l1                        3
  -3 + l2            a[2, 4]
  l3                 --------
  3 + l4      1         6

  l1                                2
  -3 + l2            a[2, 3] a[2, 4]
  1 + l3             ----------------
  2 + l4      1             2

  l1                        2
  -3 + l2            a[2, 3]  a[2, 4]
  2 + l3             ----------------
  1 + l4      1             2

  l1                        3
  -3 + l2            a[2, 3]
  3 + l3             --------
  l4          1         6
```

## 2026-07-25 14:46:28
**Input:**
```wolfram
with2 = Select[Keys@sectors[testrep, 2], AllTrue[sectorVec[testrep, testrep, #], PossibleZeroQ] === False &]; Column[{Length[with2], Table[{nu, Min[{testrep[[1]] - nu[[1]], testrep[[2]] - nu[[2]], nu[[3]] - testrep[[3]], nu[[4]] - testrep[[4]]}] + 1, DeleteDuplicates[Pick[allWts, Unitize[sectorVec[testrep, testrep, nu] /. a[__] -> 1], 1]] === {nu}}, {nu, with2}]}]
```
**Output:**
```
Out= Column[{8, {{{3, 0, 1, 2}, 1, True}, {{3, 0, 2, 1}, 1, True}, {{2, 1, 1, 2}, 1, True}, {{2, 1, 2, 1}, 2, True}, {{1, 2, 1, 2}, 1, True}, {{2, 1, 3, 0}, 1, True}, {{1, 2, 2, 1}, 1, True}, {{1, 2, 3, 0}, 1, True}}}]
```
## 2026-07-25 14:55:06 — 📥 BULK INSERT 2 CELL(S) at positions 14–15
1. [markdown] ## Overlap with a fixed-weight ket: $\langle HWS|e^{\sum a_{ij}E_{ij}}|\Psi\rangle$↵↵If $|\Psi\rangl
2. [code] (* Margins forced by the target weight, and the contributing exponent tuples. *)↵overlapMargins[lam_


## 2026-07-25 14:55:06
**Input:**
```wolfram
(* Margins forced by the target weight, and the contributing exponent tuples. *)
overlapMargins[lam_, mu_] := {lam[[1]] - mu[[1]], lam[[2]] - mu[[2]], mu[[3]] - lam[[3]], mu[[4]] - lam[[4]]};

(* The exponent tuples n contributing to <HWS|exp(..)|Psi_mu>; {} if the overlap vanishes
   on weight grounds (some margin negative, or total weight not conserved). *)
overlapTerms[lam_, mu_] := Module[{m = overlapMargins[lam, mu]},
   If[Total[lam] =!= Total[mu] || Min[m] < 0, {}, marginSector[m]]];

(* Symbolic overlap: sum of coef[n] * <HWS|E^n|Psi>, with the bra monomial left as bra[n]. *)
overlapSymbolic[lam_, mu_] := Total[coef[#] bra @@ # & /@ overlapTerms[lam, mu]];

(* Concrete overlap in a given irrep, for an explicit weight-mu ket Psi. *)
overlapVal[rep_, lam_, mu_, Psi_] := Total[coef[#] (braVec[rep, #] . Psi) & /@ overlapTerms[lam, mu]];

{overlapTerms[{l1, l2, l3, l4}, {l1 - 2, l2 - 2, l3 + 2, l4 + 2}], overlapSymbolic[{l1, l2, l3, l4}, {l1 - 2, l2 - 2, l3 + 2, l4 + 2}]}
```
**Output:**
```
(no output)
```

## 2026-07-25 14:55:10
**Input:**
```wolfram
{overlapTerms[{l1, l2, l3, l4}, {l1 - 2, l2 - 2, l3 + 2, l4 + 2}], overlapSymbolic[{l1, l2, l3, l4}, {l1 - 2, l2 - 2, l3 + 2, l4 + 2}]}
```
**Output:**
```
Out= {{{0, 2, 2, 0}, {1, 1, 1, 1}, {2, 0, 0, 2}}, (a[1, 4]^2*a[2, 3]^2*bra[0, 2, 2, 0])/4 + a[1, 3]*a[1, 4]*a[2, 3]*a[2, 4]*bra[1, 1, 1, 1] + (a[1, 3]^2*a[2, 4]^2*bra[2, 0, 0, 2])/4}
```
## 2026-07-25 14:55:28 — 📥 BULK INSERT 2 CELL(S) at positions 16–17
1. [markdown] ### Verification against the brute-force exponential↵↵For $\lambda=\{3,2,1,0\}$: build $M=e^{\sum a_
2. [code] (* Brute force: full matrix exponential at rational a-values. *)↵SeedRandom[42];↵avals = {7/5, -3/2,


## 2026-07-25 14:55:30
**Input:**
```wolfram
(* Brute force: full matrix exponential at rational a-values. *)
SeedRandom[42];
avals = {7/5, -3/2, 11/7, 2/9};
numRule = Thread[{a[1, 3], a[1, 4], a[2, 3], a[2, 4]} -> avals];
Mfull = MatrixExp[avals . gens];
braFull = UnitVector[dim[testrep], 1] . Mfull;

(* every weight occurring in the irrep, with the positions of its GT basis vectors *)
wtPos = GroupBy[Range[dim[testrep]], allWts[[#]] &];

check = Table[
   Module[{Psi = Normal@SparseArray[Thread[wtPos[mu] -> RandomInteger[{-4, 4}, Length[wtPos[mu]]]], dim[testrep]]},
    {mu, Length[wtPos[mu]], Length[overlapTerms[testrep, mu]],
     braFull . Psi - (overlapVal[testrep, testrep, mu, Psi] /. numRule)}],
   {mu, Keys[wtPos]}];

{Length[check], Max[Abs[check[[All, 4]]]], Total[check[[All, 3]]]}
```
**Output:**
```
(no output)
```

## 2026-07-25 14:55:35
**Input:**
```wolfram
{Length[check], Max[Abs[check[[All, 4]]]], Total[check[[All, 3]]], Count[check[[All,3]], 0], Max[check[[All,3]]]}
```
**Output:**
```
Out= {38, 0, 36, 12, 3}
```
## 2026-07-25 15:00:14 — 📥 BULK INSERT 3 CELL(S) at positions 19–21
1. [markdown] ### Worked example: $\lambda=\{16,16,0,0\}$, $\mu=\{8,8,8,8\}$↵↵The rectangular case relevant to the
2. [code] lamEx = {16, 16, 0, 0}; muEx = {8, 8, 8, 8};↵{overlapMargins[lamEx, muEx], Length[overlapTerms[lamEx
3. [code] (* The nine contributing monomials, with their multinomial prefactors. *)↵TableForm[↵  {#, Total[#],


## 2026-07-25 15:00:14
**Input:**
```wolfram
lamEx = {16, 16, 0, 0}; muEx = {8, 8, 8, 8};
{overlapMargins[lamEx, muEx], Length[overlapTerms[lamEx, muEx]]}
```
**Output:**
```
Out[100]= \left\{\left\{8,8,8,8\right\},9\right\}
```

## 2026-07-25 15:00:14
**Input:**
```wolfram
(* The nine contributing monomials, with their multinomial prefactors. *)
TableForm[
  {#, Total[#], coef[#], 1/(Times @@ (#!))} & /@ overlapTerms[lamEx, muEx],
  TableHeadings -> {None, {"{n13,n14,n23,n24}", "deg", "coef[n]", "prefactor"}}]
```
**Output:**
```
Out[101]//TableForm= \begin{pmatrix}\left\{0,8,8,0\right\} & 16 & \frac{a\left(1,4\right)^8\,a\left(2,3\right)^8}{1625702400} & \frac{1}{1625702400}\\\left\{1,7,7,1\right\} & 16 & \frac{a\left(1,3\right)\,a\left(1,4\right)^7\,a\left(2,3\right)^7\,a\left(2,4\right)}{25401600} & \frac{1}{25401600}\\\left\{2,6,6,2\right\} & 16 & \frac{a\left(1,3\right)^2\,a\left(1,4\right)^6\,a\left(2,3\right)^6\,a\l… [1258 chars total]
```
## 2026-07-25 15:00:30 — 📥 BULK INSERT 3 CELL(S) at positions 22–24
1. [code] (* The middle-weight overlap, written out. *)↵overlapSymbolic[lamEx, muEx]
2. [markdown] The $\lambda=\{16,16,0,0\}$ irrep has dimension 8721, too big for a comfortable brute-force $e^{M}$.
3. [code] (* Brute-force check of the rectangular middle-weight overlap at S=4. *)↵SeedRandom[7];↵lamR = {4, 4


## 2026-07-25 15:00:30
**Input:**
```wolfram
(* The middle-weight overlap, written out. *)
overlapSymbolic[lamEx, muEx]
```
**Output:**
```
Out[102]= \frac{a\left(1,4\right)^8\,a\left(2,3\right)^8\,\mathrm{bra}\left(0,8,8,0\right)}{1625702400}+\frac{a\left(1,3\right)\,a\left(1,4\right)^7\,a\left(2,4\right)\,a\left(2,3\right)^7\,\mathrm{bra}\left(1,7,7,1\right)}{25401600}+\frac{a\left(1,3\right)^2\,a\left(1,4\right)^6\,a\left(2,4\right)^2\,a\left(2,3\right)^6\,\mathrm{bra}\left(2,6,6,2\right)}{2073600}+\frac{a\left(1,3\right)^3\,a\left… [1119 chars total]
```

## 2026-07-25 15:00:40
**Input:**
```wolfram
(* Brute-force check of the rectangular middle-weight overlap at S=4. *)
SeedRandom[7];
lamR = {4, 4, 0, 0}; muR = {2, 2, 2, 2};
gensR = EE[lamR] @@@ {{1, 3}, {1, 4}, {2, 3}, {2, 4}};
wtsR = (J0[#1][#2] & @@@ Thread[{Range[4], #}] &) /@ ps[{lamR}];
wtsR = Table[J0[k][pat] & /@ Range[4] // Quiet, {pat, ps[{lamR}]}];
posR = Flatten@Position[wtsR, muR];
PsiR = Normal@SparseArray[Thread[posR -> RandomInteger[{-4, 4}, Length[posR]]], dim[lamR]];
braR = UnitVector[dim[lamR], 1] . MatrixExp[avals . gensR];
{dim[lamR], Length[posR], overlapTerms[lamR, muR],
 braR . PsiR - (overlapVal[lamR, lamR, muR, PsiR] /. numRule)}
```
**Output:**
```
4 kernel messages: Part::partd1: Depth of object {4} is not sufficient for the given part specification.
Part::partd1: Depth of object {4, 4} is not sufficient for the given part specification.
Part::partd1: Depth of ob… [also: Part::pkspec1, Part::pkspec1, Part::partd1] result: Out[111]= \left\{105,0,\begin{pmatrix}0 & 2 & 2 & 0\\1 & 1 & 1 & 1\\2 & 0 & 0 & 2\end{pmatrix},0\right\}
```

## 2026-07-25 15:00:52
**Input:**
```wolfram
lamR = {4, 4, 0, 0}; wtsR = patWt /@ ps[{lamR}]; {dim[lamR], Length[wtsR], First[wtsR], Count[wtsR, {2, 2, 2, 2}]}
```
**Output:**
```
Out= {105, 105, {4, 4, 0, 0}, 3}
```
## 2026-07-25 15:00:59 — ✏️ EDIT CELL 24
(* Brute-force check of the rectangular middle-weight overlap one rung down (S=4). *)
SeedRandom[7];
lamR = {4, 4, 0, 0}; muR = {2, 2, 2, 2};
gensR = EE[lamR] @@@ {{1, 3}, {1, 4}, {2, 3}, {2, 4}};
wts…


## 2026-07-25 15:01:05
**Input:**
```wolfram
(* Brute-force check of the rectangular middle-weight overlap one rung down (S=4). *)
SeedRandom[7];
lamR = {4, 4, 0, 0}; muR = {2, 2, 2, 2};
gensR = EE[lamR] @@@ {{1, 3}, {1, 4}, {2, 3}, {2, 4}};
wtsR = patWt /@ ps[{lamR}];
posR = Flatten@Position[wtsR, muR];
PsiR = Normal@SparseArray[Thread[posR -> RandomInteger[{-4, 4}, Length[posR]]], dim[lamR]];
braR = UnitVector[dim[lamR], 1] . MatrixExp[avals . gensR];
{dim[lamR], Length[posR], Length[overlapTerms[lamR, muR]],
 braR . PsiR - (overlapVal[lamR, lamR, muR, PsiR] /. numRule)}
```
**Output:**
```
Out= {105, 3, 3, 0}
```

## 2026-07-25 15:01:09
**Input:**
```wolfram
Column[{#, Total[#], coef[#]} & /@ overlapTerms[{16,16,0,0},{8,8,8,8}]]
```
**Output:**
```
Out= Column[{{{0, 8, 8, 0}, 16, (a[1, 4]^8*a[2, 3]^8)/1625702400}, {{1, 7, 7, 1}, 16, (a[1, 3]*a[1, 4]^7*a[2, 3]^7*a[2, 4])/25401600}, {{2, 6, 6, 2}, 16, (a[1, 3]^2*a[1, 4]^6*a[2, 3]^6*a[2, 4]^2)/2073600}, {{3, 5, 5, 3}, 16, (a[1, 3]^3*a[1, 4]^5*a[2, 3]^5*a[2, 4]^3)/518400}, {{4, 4, 4, 4}, 16, (a[1, 3]^4*a[1, 4]^4*a[2, 3]^4*a[2, 4]^4)/331776}, {{5, 3, 3, 5}, 16, (a[1, 3]^5*a[1, 4]^3*a[2, 3]^3*a[2, 4]^5)/518400}, {{6, 2, 2, 6}, 16, (a[1, 3]^6*a[1, 4]^2*a[2, 3]^2*a[2, 4]^6)/2073600}, {{7, 1, 1, 7}, 16, (a[1, 3]^7*a[1, 4]*a[2, 3]*a[2, 4]^7)/25401600}, {{8, 0, 0, 8}, 16, (a[1, 3]^8*a[2, 4]^8)/1625702400}}]
```

## 2026-07-25 15:07:39
**Input:**
```wolfram
rect[S_] := {S, S, 0, 0};
braKill[rep_] := {Max@Abs[UnitVector[dim[rep], 1] . EE[rep][1, 2]], Max@Abs[UnitVector[dim[rep], 1] . EE[rep][3, 4]], Max@Abs[UnitVector[dim[rep], 1] . EE[rep][2, 1]], Max@Abs[UnitVector[dim[rep], 1] . EE[rep][4, 3]]};
{braKill[rect[1]], braKill[rect[2]], braKill[rect[3]], braKill[{3, 2, 1, 0}], braKill[{2, 1, 1, 0}]}
```
**Output:**
```
Out= {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {1, 1, 0, 0}, {1, 1, 0, 0}}
```

## 2026-07-25 15:08:14
**Input:**
```wolfram
spanTest[S_] := Module[{rep = rect[S], all, mat},
  all = Flatten[Table[nList[d], {d, 0, 4 S}], 1];
  mat = braVec[rep, #] & /@ all;
  {S, dim[rep], Length[all], MatrixRank[mat], Count[mat, x_ /; ! AllTrue[x, PossibleZeroQ]]}];
TableForm[spanTest /@ {1, 2, 3}, TableHeadings -> {None, {"S", "dim rep", "#monomials", "rank", "#nonzero monomials"}}]
```
**Output:**
```
Out= TableForm[{{1, 6, 70, 6, 7}, {2, 20, 495, 20, 26}, {3, 50, 1820, 50, 70}}, TableHeadings -> {None, {"S", "dim rep", "#monomials", "rank", "#nonzero monomials"}}]
```

## 2026-07-25 15:08:37
**Input:**
```wolfram
detOp[rep_] := EE[rep][1, 3] . EE[rep][2, 4] - EE[rep][1, 4] . EE[rep][2, 3];
Table[{S, Max@Abs[UnitVector[dim[rect[S]], 1] . detOp[rect[S]]], Max@Abs[UnitVector[dim[rect[S]], 1] . MatrixPower[detOp[rect[S]], S + 1]]}, {S, 1, 3}]
```
**Output:**
```
Out= {{1, 8, 0}, {2, 36, 0}, {3, 96, 0}}
```

## 2026-07-25 15:08:56
**Input:**
```wolfram
wtMultTable[S_] := Module[{rep = rect[S], wts, g},
  wts = patWt /@ ps[{rep}]; g = Counts[wts];
  Table[{mu, g[mu], Length[overlapTerms[rep, mu]], g[mu] === Length[overlapTerms[rep, mu]]}, {mu, Keys[g]}]];
tab3 = wtMultTable[3];
{Length[tab3], And @@ tab3[[All, 4]], Select[tab3, ! #[[4]] &][[1 ;; UpTo[8]]]}
```
**Output:**
```
Out= {44, False, {{{2, 0, 3, 1}, 1, 2, False}, {{2, 0, 2, 2}, 1, 2, False}, {{2, 0, 1, 3}, 1, 2, False}, {{1, 1, 3, 1}, 1, 2, False}, {{1, 1, 2, 2}, 2, 3, False}, {{1, 1, 1, 3}, 1, 2, False}, {{1, 0, 3, 2}, 1, 3, False}, {{1, 0, 2, 3}, 1, 3, False}}}
```

## 2026-07-25 15:09:06
**Input:**
```wolfram
multFormula[S_, mu_] := Module[{m = overlapMargins[rect[S], mu]}, If[Min[m] < 0 || Total[mu] =!= 2 S, 0, Min[Min[m] + 1, S + 1 - Max[m]]];
checkFormula[S_] := Module[{rep = rect[S], g}, g = Counts[patWt /@ ps[{rep}]];
  {S, Length[g], And @@ (g[#] === multFormula[S, #] & /@ Keys[g]),
   Total[Length[overlapTerms[rep, #]] & /@ Keys[g]], Total[Values[g]]}];
TableForm[checkFormula /@ {1, 2, 3, 4}, TableHeadings -> {None, {"S", "#weights", "mult formula holds", "total #monomials", "dim"}}]
```
**Output:**
```
[message] ToExpression::sntx: Invalid syntax in or before "TableForm[checkFormula /@ {1, 2, 3, 4}, TableHeadings -> {None, {"S", "#weights", "mult formula holds", "total #monomials", "dim"}}]".
(no output)
```

## 2026-07-25 15:09:12
**Input:**
```wolfram
multFormula[S_, mu_] := Module[{m = overlapMargins[rect[S], mu]}, If[Min[m] < 0 || Total[mu] =!= 2 S, 0, Min[Min[m] + 1, S + 1 - Max[m]]]];
checkFormula[S_] := Module[{rep = rect[S], g}, g = Counts[patWt /@ ps[{rep}]]; {S, Length[g], And @@ (g[#] === multFormula[S, #] & /@ Keys[g]), Total[Length[overlapTerms[rep, #]] & /@ Keys[g]], Total[Values[g]]}];
checkFormula /@ {1, 2, 3, 4}
```
**Output:**
```
Out= {{1, 6, True, 7, 6}, {2, 19, True, 26, 20}, {3, 44, True, 70, 50}, {4, 85, True, 155, 105}}
```

## 2026-07-25 15:09:30
**Input:**
```wolfram
basisSel[S_] := Module[{rep = rect[S], g, sel}, g = Counts[patWt /@ ps[{rep}]];
  sel = Flatten[Table[Take[overlapTerms[rep, mu], multFormula[S, mu]], {mu, Keys[g]}], 1];
  {S, Length[sel], dim[rep], MatrixRank[braVec[rep, #] & /@ sel]}];
basisSel /@ {1, 2, 3}
```
**Output:**
```
Out= {{1, 6, 6, 6}, {2, 20, 20, 20}, {3, 50, 50, 50}}
```
## 2026-07-25 15:09:56 — 📥 BULK INSERT 5 CELL(S) at positions 25–29
1. [markdown] ## Rectangular reps $[S,S,0,0]$: monomial labels as CG indices↵↵For $\lambda=[S,S,0,0]$ the weights 
2. [code] rect[S_] := {S, S, 0, 0};↵↵(* (a) The extra annihilation conditions, absent for non-rectangular reps
3. [code] (* (b) The four commuting generators really do span the whole dual module: rank == dim.↵   But the n
4. [code] (* (c) Weight multiplicity vs number of monomials in the sector, and the basis truncation. *)↵multFo
5. [code] (* (d) The determinant operator D = E13 E24 - E14 E23 is nonzero on the bra, but nilpotent↵   of ord


## 2026-07-25 15:09:57
**Input:**
```wolfram
rect[S_] := {S, S, 0, 0};

(* (a) The extra annihilation conditions, absent for non-rectangular reps.
   Columns: <HWS|E12, <HWS|E34, <HWS|E21, <HWS|E43 (max abs entry). *)
braKill[rep_] := Max@Abs[UnitVector[dim[rep], 1] . EE[rep] @@ #] & /@ {{1, 2}, {3, 4}, {2, 1}, {4, 3}};
TableForm[{#, braKill[#]} & /@ {rect[1], rect[2], rect[3], {3, 2, 1, 0}, {2, 1, 1, 0}},
  TableHeadings -> {None, {"rep", "{E12, E34, E21, E43}"}}]
```
**Output:**
```
Out[114]//TableForm= \begin{pmatrix}\left\{1,1,0,0\right\} & \left\{0,0,0,0\right\}\\\left\{2,2,0,0\right\} & \left\{0,0,0,0\right\}\\\left\{3,3,0,0\right\} & \left\{0,0,0,0\right\}\\\left\{3,2,1,0\right\} & \left\{1,1,0,0\right\}\\\left\{2,1,1,0\right\} & \left\{1,1,0,0\right\}\end{pmatrix}
```

## 2026-07-25 15:10:20
**Input:**
```wolfram
(* (b) The four commuting generators really do span the whole dual module: rank == dim.
   But the nonzero monomials outnumber dim => linear relations. *)
spanTest[S_] := Module[{rep = rect[S], all, mat},
   all = Flatten[Table[nList[d], {d, 0, 4 S}], 1];
   mat = braVec[rep, #] & /@ all;
   {S, dim[rep], Length[all], MatrixRank[mat], Count[mat, x_ /; ! AllTrue[x, PossibleZeroQ]]}];
TableForm[spanTest /@ {1, 2, 3},
  TableHeadings -> {None, {"S", "dim", "#monomials", "rank", "#nonzero monomials"}}]
```
**Output:**
```
Out[116]//TableForm= \begin{pmatrix}1 & 6 & 70 & 6 & 7\\2 & 20 & 495 & 20 & 26\\3 & 50 & 1820 & 50 & 70\end{pmatrix}
```

## 2026-07-25 15:10:20
**Input:**
```wolfram
(* (c) Weight multiplicity vs number of monomials in the sector, and the basis truncation. *)
multFormula[S_, mu_] := Module[{m = overlapMargins[rect[S], mu]},
   If[Min[m] < 0 || Total[mu] =!= 2 S, 0, Min[Min[m] + 1, S + 1 - Max[m]]]];

checkFormula[S_] := Module[{rep = rect[S], g}, g = Counts[patWt /@ ps[{rep}]];
   {S, Length[g], And @@ (g[#] === multFormula[S, #] & /@ Keys[g]),
    Total[Length[overlapTerms[rep, #]] & /@ Keys[g]], Total[Values[g]]}];

(* keeping the first mult(mu) values of t per sector gives an exact basis *)
basisSel[S_] := Module[{rep = rect[S], g, sel}, g = Counts[patWt /@ ps[{rep}]];
   sel = Flatten[Table[Take[overlapTerms[rep, mu], multFormula[S, mu]], {mu, Keys[g]}], 1];
   {S, Length[sel], dim[rep], MatrixRank[braVec[rep, #] & /@ sel]}];

Column[{TableForm[checkFormula /@ {1, 2, 3, 4}, TableHeadings -> {None, {"S", "#weights", "mult formula ok", "#monomials", "dim"}}],
  TableForm[basisSel /@ {1, 2, 3}, TableHeadings -> {None, {"S", "#selected", "dim", "rank"}}]}]
```
**Output:**
```
Out[120]= \begin{matrix}\begin{aligned}\text{S} & \text{\#weights} & \text{mult formula ok} & \text{\#monomials} & \text{dim}\\1 & 6 & \mathrm{True} & 7 & 6\\2 & 19 & \mathrm{True} & 26 & 20\\3 & 44 & \mathrm{True} & 70 & 50\\4 & 85 & \mathrm{True} & 155 & 105\end{aligned}\\\begin{aligned}\text{S} & \text{\#selected} & \text{dim} & \text{rank}\\1 & 6 & 6 & 6\\2 & 20 & 20 & 20\\3 & 50 & 50 & 50\end… [421 chars total]
```

## 2026-07-25 15:10:21
**Input:**
```wolfram
(* (d) The determinant operator D = E13 E24 - E14 E23 is nonzero on the bra, but nilpotent
   of order S+1 -- the source of the relations. *)
detOp[rep_] := EE[rep][1, 3] . EE[rep][2, 4] - EE[rep][1, 4] . EE[rep][2, 3];
TableForm[Table[{S, Max@Abs[UnitVector[dim[rect[S]], 1] . detOp[rect[S]]],
    Max@Abs[UnitVector[dim[rect[S]], 1] . MatrixPower[detOp[rect[S]], S + 1]]}, {S, 1, 3}],
  TableHeadings -> {None, {"S", "|<HWS|D|", "|<HWS|D^(S+1)|"}}]
```
**Output:**
```
Out[122]//TableForm= \begin{pmatrix}1 & 8 & 0\\2 & 36 & 0\\3 & 96 & 0\end{pmatrix}
```

## 2026-07-25 15:14:57
**Input:**
```wolfram
rankOne[rep_, {u1_, u2_}, {v3_, v4_}] := u1 v3 EE[rep][1, 3] + u1 v4 EE[rep][1, 4] + u2 v3 EE[rep][2, 3] + u2 v4 EE[rep][2, 4];
Table[With[{rep = rect[S], X = rankOne[rect[S], {u1, u2}, {v3, v4}]},
  {S, Max@Abs@Expand[UnitVector[dim[rep], 1] . MatrixPower[X, S + 1]],
      Max@Abs@Expand[UnitVector[dim[rep], 1] . MatrixPower[X, S]]}], {S, 1, 3}]
```
**Output:**
```
Out= {{1, 0, Max[0, Abs[u1*v3], Abs[u2*v3], 2*Abs[u1*v4], 2*Abs[u2*v4]]}, {2, 0, Max[0, 8*Abs[u1^2*v3^2], 8*Abs[u1*u2*v3^2], 4*Abs[u2^2*v3^2], 24*Abs[u1^2*v3*v4], 24*Abs[u1*u2*v3*v4], 12*Abs[u2^2*v3*v4], 48*Abs[u1^2*v4^2], 48*Abs[u1*u2*v4^2], 24*Abs[u2^2*v4^2]]}, {3, 0, Max[0, 216*Abs[u1^3*v3^3], 216*Abs[u1^2*u2*v3^3], 108*Abs[u1*u2^2*v3^3], 36*Abs[u2^3*v3^3], 864*Abs[u1^3*v3^2*v4], 864*Abs[u1^2*u2*v3^2*v4], 432*Abs[u1*u2^2*v3^2*v4], 144*Abs[u2^3*v3^2*v4], 2592*Abs[u1^3*v3*v4^2], 2592*Abs[u1^2*u2*v3*v4^2], 1296*Abs[u1*u2^2*v3*v4^2], 432*Abs[u2^3*v3*v4^2], 5184*Abs[u1^3*v4^3], 5184*Abs[u1^2*u2*v4^3], 2592*Abs[u1*u2^2*v4^3], 864*Abs[u2^3*v4^3]]}}
```

## 2026-07-25 15:15:06
**Input:**
```wolfram
graded[S_] := Module[{rep = rect[S], prev = {}, res = {}, cur, r0, r1},
  Do[cur = Join[prev, braVec[rep, #] & /@ nList[k]];
     r0 = MatrixRank[prev]; r1 = MatrixRank[cur];
     AppendTo[res, {k, Length[nList[k]], r1 - r0, Length[NullSpace[Transpose[braVec[rep, #] & /@ nList[k]]]]}];
     prev = cur, {k, 0, 2 S + 2}];
  res];
TableForm[graded[2], TableHeadings -> {None, {"deg k", "#monomials", "new dim", "dim K_k"}}]
```
**Output:**
```
[message] MatrixRank::matrix: Argument {} at position 1 is not a non-empty rectangular matrix.
Out= TableForm[{{0, 1, 1 - MatrixRank[{}], 0}, {1, 4, 4, 0}, {2, 10, 10, 0}, {3, 20, 4, 16}, {4, 35, 1, 34}, {5, 56, 0, 56}, {6, 84, 0, 84}}, TableHeadings -> {None, {"deg k", "#monomials", "new dim", "dim K_k"}}]
```

## 2026-07-25 15:15:30
**Input:**
```wolfram
gradedDims[S_] := Module[{rep = rect[S], prev = {ConstantArray[0, dim[rect[S]]]}, res = {}, cur, r0, r1},
  Do[cur = Join[prev, braVec[rep, #] & /@ nList[k]]; r0 = MatrixRank[prev]; r1 = MatrixRank[cur];
     AppendTo[res, r1 - r0]; prev = cur, {k, 0, 2 S + 2}]; res];
predDims[S_] := Table[Total[Table[If[0 <= k - 2 j <= S && j <= S - (k - 2 j), (k - 2 j + 1)^2, 0], {j, 0, Floor[k/2]}]], {k, 0, 2 S + 2}];
Table[{S, gradedDims[S], predDims[S], gradedDims[S] === predDims[S], Total[gradedDims[S]], dim[rect[S]]}, {S, 1, 3}]
```
**Output:**
```
Out= {{1, {1, 4, 1, 0, 0}, {1, 4, 1, 0, 0}, True, 6, 6}, {2, {1, 4, 10, 4, 1, 0, 0}, {1, 4, 10, 4, 1, 0, 0}, True, 20, 20}, {3, {1, 4, 10, 20, 10, 4, 1, 0, 0}, {1, 4, 10, 20, 10, 4, 1, 0, 0}, True, 50, 50}}
```

## 2026-07-25 15:15:43
**Input:**
```wolfram
monBasis[k_] := nList[k];
kerAt[S_, k_] := NullSpace[Transpose[braVec[rect[S], #] & /@ nList[k]]];
(* multiply a degree-k kernel vector by each x_ij, expressed in the degree-(k+1) monomial basis *)
shift[k_, vec_, i_] := Module[{tgt = nList[k + 1], src = nList[k], pos},
   pos = (Position[tgt, # + UnitVector[4, i]][[1, 1]] &) /@ src;
   Normal@SparseArray[Thread[pos -> vec], Length[tgt]]];
genCheck[S_, k_] := Module[{Kk = kerAt[S, k], Kk1 = kerAt[S, k + 1], gen},
   gen = Flatten[Table[shift[k, v, i], {v, Kk}, {i, 4}], 1];
   {k, Length[Kk], Length[Kk1], If[gen === {}, 0, MatrixRank[gen]]}];
Table[{S, (S + 2)^2, genCheck[S, #] & /@ Range[S + 1, S + 3]}, {S, 1, 2}]
```
**Output:**
```
Out= {{1, 9, {{2, 9, 20, 20}, {3, 20, 35, 35}, {4, 35, 56, 56}}}, {2, 16, {{3, 16, 34, 34}, {4, 34, 56, 56}, {5, 56, 84, 84}}}}
```

## 2026-07-25 15:16:02
**Input:**
```wolfram
Rrel[S_, m_] := Total[(S + 1)!/(Times @@ (#!)) braVec[rect[S], #] & /@ marginSector[m]];
Table[{S, Max[Abs[Flatten[Table[Rrel[S, {m1, S + 1 - m1, m3, S + 1 - m3}], {m1, 0, S + 1}, {m3, 0, S + 1}]]]],
      Max[Abs[Flatten[Table[Rrel[S - 1 /. x_ :> S, {m1, S - m1, m3, S - m3}] , {m1, 0, S}, {m3, 0, S}]]]]}, {S, 1, 3}]
```
**Output:**
```
Out= {{1, 0, 4}, {2, 0, 144}, {3, 0, 20736}}
```
## 2026-07-25 15:16:41 — 📥 BULK INSERT 5 CELL(S) at positions 30–34
1. [markdown] ## Exact classification of the relations↵↵Write $x_{ij}=E_{ij}$ for the four commuting nilradical ge
2. [code] (* (e) The relations are the rank-one degenerations: <HWS|(sum u_i v_j E_ij)^(S+1) == 0↵   identical
3. [code] (* (f) Explicit sector-by-sector relations R_m, and sharpness of the degree S+1. *)↵Rrel[S_, m_] := 
4. [code] (* (g) Graded Levi decomposition: degree-k dimensions vs the D^j Sym^a (x) Sym^a prediction. *)↵grad
5. [code] (* (h) No generators below degree S+1, and none above: rank(x . K_k) == dim K_(k+1). *)↵kerAt[S_, k_


## 2026-07-25 15:16:46
**Input:**
```wolfram
(* (e) The relations are the rank-one degenerations: <HWS|(sum u_i v_j E_ij)^(S+1) == 0
   identically in u,v, while the S-th power does not vanish. *)
rankOne[rep_, {u1_, u2_}, {v3_, v4_}] := u1 v3 EE[rep][1, 3] + u1 v4 EE[rep][1, 4] + u2 v3 EE[rep][2, 3] + u2 v4 EE[rep][2, 4];
TableForm[Table[With[{rep = rect[S], X = rankOne[rect[S], {u1, u2}, {v3, v4}]},
    {S, Max@Abs@Expand[UnitVector[dim[rep], 1] . MatrixPower[X, S + 1]],
     Simplify[Max@Abs@Expand[UnitVector[dim[rep], 1] . MatrixPower[X, S]]] =!= 0}], {S, 1, 3}],
  TableHeadings -> {None, {"S", "|<HWS|X^(S+1)|", "X^S nonzero"}}]
```
**Output:**
```
Out[124]//TableForm= \begin{pmatrix}1 & 0 & \mathrm{True}\\2 & 0 & \mathrm{True}\\3 & 0 & \mathrm{True}\end{pmatrix}
```

## 2026-07-25 15:16:47
**Input:**
```wolfram
(* (f) Explicit sector-by-sector relations R_m, and sharpness of the degree S+1. *)
Rrel[S_, m_] := Total[(S + 1)!/(Times @@ (#!)) braVec[rect[S], #] & /@ marginSector[m]];
TableForm[Table[{S,
    Max[Abs[Flatten[Table[Rrel[S, {m1, S + 1 - m1, m3, S + 1 - m3}], {m1, 0, S + 1}, {m3, 0, S + 1}]]]],
    Max[Abs[Flatten[Table[Rrel[S, {m1, S - m1, m3, S - m3}], {m1, 0, S}, {m3, 0, S}]]]]}, {S, 1, 3}],
  TableHeadings -> {None, {"S", "max|R_m|, deg S+1", "max|R_m|, deg S"}}]
```
**Output:**
```
Out[126]//TableForm= \begin{pmatrix}1 & 0 & 4\\2 & 0 & 144\\3 & 0 & 20736\end{pmatrix}
```

## 2026-07-25 15:17:00
**Input:**
```wolfram
(* (g) Graded Levi decomposition: degree-k dimensions vs the D^j Sym^a (x) Sym^a prediction. *)
gradedDims[S_] := Module[{rep = rect[S], prev = {ConstantArray[0, dim[rect[S]]]}, res = {}, cur, r0, r1},
   Do[cur = Join[prev, braVec[rep, #] & /@ nList[k]]; r0 = MatrixRank[prev]; r1 = MatrixRank[cur];
      AppendTo[res, r1 - r0]; prev = cur, {k, 0, 2 S + 2}]; res];
predDims[S_] := Table[Total[Table[If[0 <= k - 2 j <= S && j <= S - (k - 2 j), (k - 2 j + 1)^2, 0], {j, 0, Floor[k/2]}]], {k, 0, 2 S + 2}];
TableForm[Table[{S, gradedDims[S], gradedDims[S] === predDims[S], Total[gradedDims[S]], dim[rect[S]]}, {S, 1, 3}],
  TableHeadings -> {None, {"S", "graded dims", "matches prediction", "total", "dim"}}]
```
**Output:**
```
Out[129]//TableForm= \begin{pmatrix}1 & \left\{1,4,1,0,0\right\} & \mathrm{True} & 6 & 6\\2 & \left\{1,4,10,4,1,0,0\right\} & \mathrm{True} & 20 & 20\\3 & \left\{1,4,10,20,10,4,1,0,0\right\} & \mathrm{True} & 50 & 50\end{pmatrix}
```

## 2026-07-25 15:17:01
**Input:**
```wolfram
(* (h) No generators below degree S+1, and none above: rank(x . K_k) == dim K_(k+1). *)
kerAt[S_, k_] := NullSpace[Transpose[braVec[rect[S], #] & /@ nList[k]]];
shift[k_, vec_, i_] := Module[{tgt = nList[k + 1], src = nList[k], pos},
   pos = (Position[tgt, # + UnitVector[4, i]][[1, 1]] &) /@ src;
   Normal@SparseArray[Thread[pos -> vec], Length[tgt]]];
genCheck[S_, k_] := Module[{Kk = kerAt[S, k], Kk1 = kerAt[S, k + 1], gen},
   gen = Flatten[Table[shift[k, v, i], {v, Kk}, {i, 4}], 1];
   {k, Length[Kk], Length[Kk1], If[gen === {}, 0, MatrixRank[gen]]}];
TableForm[Flatten[Table[Prepend[genCheck[S, k], S], {S, 1, 2}, {k, S, S + 3}], 1],
  TableHeadings -> {None, {"S", "k", "dim K_k", "dim K_(k+1)", "rank(x.K_k)"}}]
```
**Output:**
```
Out[133]//TableForm= \begin{pmatrix}1 & 1 & 0 & 9 & 0\\1 & 2 & 9 & 20 & 20\\1 & 3 & 20 & 35 & 35\\1 & 4 & 35 & 56 & 56\\2 & 2 & 0 & 16 & 0\\2 & 3 & 16 & 34 & 34\\2 & 4 & 34 & 56 & 56\\2 & 5 & 56 & 84 & 84\end{pmatrix}
```

## 2026-07-25 15:21:03
**Input:**
```wolfram
plueckerLabels[S_] := Select[Flatten[Table[{e13, e14, e23, e24, e34}, {e13, 0, S}, {e14, 0, S}, {e23, 0, S}, {e24, 0, S}, {e34, 0, S}], 4], Total[#] <= S && Min[#[[2]], #[[3]]] == 0 &];
pluVec[rep_, {e13_, e14_, e23_, e24_, e34_}] := UnitVector[dim[rep], 1] . mpow[EE[rep][1, 3], e13] . mpow[EE[rep][1, 4], e14] . mpow[EE[rep][2, 3], e23] . mpow[EE[rep][2, 4], e24] . mpow[detOp[rep], e34];
Table[With[{L = plueckerLabels[S]}, {S, Length[L], dim[rect[S]], MatrixRank[pluVec[rect[S], #] & /@ L]}], {S, 1, 3}]
```
**Output:**
```
Out= {{1, 6, 6, 6}, {2, 20, 20, 20}, {3, 50, 50, 50}}
```

## 2026-07-25 15:21:13
**Input:**
```wolfram
symVec[rep_, {a_, b_, c_}] := UnitVector[dim[rep], 1] . mpow[EE[rep][1, 2], a] . mpow[EE[rep][1, 3], b] . mpow[EE[rep][1, 4], c];
symTest[S_] := Module[{rep = {S, 0, 0, 0}, L},
  L = Select[Flatten[Table[{a, b, c}, {a, 0, S}, {b, 0, S}, {c, 0, S}], 2], Total[#] <= S &];
  {S, dim[rep], Binomial[S + 3, 3], Length[L], MatrixRank[symVec[rep, #] & /@ L], Max@Abs@Flatten[{UnitVector[dim[rep],1].EE[rep][2,3], UnitVector[dim[rep],1].EE[rep][2,4], UnitVector[dim[rep],1].EE[rep][3,4]}]}];
Column[{symTest /@ {2, 3, 4},
  Table[{S, (S + 1) (S + 2)^2 (S + 3)/12, Binomial[S + 5, 5] - Binomial[S + 3, 5]}, {S, 1, 6}]}]
```
**Output:**
```
Out= Column[{{{2, 10, 10, 10, 10, 0}, {3, 20, 20, 20, 20, 0}, {4, 35, 35, 35, 35, 0}}, {{1, 6, 6}, {2, 20, 20}, {3, 50, 50}, {4, 105, 105}, {5, 196, 196}, {6, 336, 336}}}]
```
## 2026-07-25 15:21:44 — 📥 BULK INSERT 4 CELL(S) at positions 35–38
1. [markdown] ## A genuine basis: Plücker / standard-monomial labels↵↵The fix is to add the determinant as a *fift
2. [code] (* (i) Pluecker / standard-monomial basis: sum e <= S with min(e14,e23)==0. *)↵plueckerLabels[S_] :=
3. [markdown] ### How rich are rectangular reps really?↵↵Contrast with the symmetric powers $[S,0,0,0]$: there the
4. [code] (* (j) Symmetric powers [S,0,0,0] for comparison: 3 commuting generators, monomials of↵   degree <=


## 2026-07-25 15:21:44
**Input:**
```wolfram
(* (i) Pluecker / standard-monomial basis: sum e <= S with min(e14,e23)==0. *)
plueckerLabels[S_] := Select[Flatten[Table[{e13, e14, e23, e24, e34}, {e13, 0, S}, {e14, 0, S}, {e23, 0, S}, {e24, 0, S}, {e34, 0, S}], 4], Total[#] <= S && Min[#[[2]], #[[3]]] == 0 &];
pluVec[rep_, {e13_, e14_, e23_, e24_, e34_}] := UnitVector[dim[rep], 1] . mpow[EE[rep][1, 3], e13] . mpow[EE[rep][1, 4], e14] . mpow[EE[rep][2, 3], e23] . mpow[EE[rep][2, 4], e24] . mpow[detOp[rep], e34];
TableForm[Table[With[{L = plueckerLabels[S]}, {S, Length[L], dim[rect[S]], MatrixRank[pluVec[rect[S], #] & /@ L]}], {S, 1, 3}],
  TableHeadings -> {None, {"S", "#labels", "dim", "rank"}}]
```
**Output:**
```
Out[136]//TableForm= \begin{pmatrix}1 & 6 & 6 & 6\\2 & 20 & 20 & 20\\3 & 50 & 50 & 50\end{pmatrix}
```

## 2026-07-25 15:21:45
**Input:**
```wolfram
(* (j) Symmetric powers [S,0,0,0] for comparison: 3 commuting generators, monomials of
   degree <= S are already a basis (no relations), and E23/E24/E34 also kill the bra.
   Then the so(6) harmonic dimension identity for the rectangular case. *)
symVec[rep_, {a_, b_, c_}] := UnitVector[dim[rep], 1] . mpow[EE[rep][1, 2], a] . mpow[EE[rep][1, 3], b] . mpow[EE[rep][1, 4], c];
symTest[S_] := Module[{rep = {S, 0, 0, 0}, L},
   L = Select[Flatten[Table[{a, b, c}, {a, 0, S}, {b, 0, S}, {c, 0, S}], 2], Total[#] <= S &];
   {S, dim[rep], Binomial[S + 3, 3], Length[L], MatrixRank[symVec[rep, #] & /@ L],
    Max@Abs@Flatten[UnitVector[dim[rep], 1] . EE[rep] @@ # & /@ {{2, 3}, {2, 4}, {3, 4}}]}];
Column[{TableForm[symTest /@ {2, 3, 4}, TableHeadings -> {None, {"S", "dim", "C(S+3,3)", "#monomials", "rank", "E23/E24/E34 on bra"}}],
  TableForm[Table[{S, (S + 1) (S + 2)^2 (S + 3)/12, Binomial[S + 5, 5] - Binomial[S + 3, 5]}, {S, 1, 6}],
   TableHeadings -> {None, {"S", "dim[S,S,0,0]", "so(6) harmonics"}}]}]
```
**Output:**
```
Out[139]= \begin{matrix}\begin{aligned}\text{S} & \text{dim} & \text{C(S+3,3)} & \text{\#monomials} & \text{rank} & \text{E23/E24/E34 on bra}\\2 & 10 & 10 & 10 & 10 & 0\\3 & 20 & 20 & 20 & 20 & 0\\4 & 35 & 35 & 35 & 35 & 0\end{aligned}\\\begin{aligned}\text{S} & \text{dim[S,S,0,0]} & \text{so(6) harmonics}\\1 & 6 & 6\\2 & 20 & 20\\3 & 50 & 50\\4 & 105 & 105\\5 & 196 & 196\\6 & 336 & 336\end{aligne… [414 chars total]
```

## 2026-07-25 17:56:02
**Input:**
```wolfram
symWt[S_] := Module[{rep = {S, 0, 0, 0}, g, L, wmap}, g = Counts[patWt /@ ps[{rep}]];
  L = Select[Flatten[Table[{a, b, c}, {a, 0, S}, {b, 0, S}, {c, 0, S}], 2], Total[#] <= S &];
  wmap = ({S, 0, 0, 0} + {-Total[#], #[[1]], #[[2]], #[[3]]}) & /@ L;
  {S, dim[rep], Length[g], Max[Values[g]], Length[L], Length[DeleteDuplicates[wmap]], Sort[DeleteDuplicates[wmap]] === Sort[Keys[g]]}];
TableForm[symWt /@ {2, 3, 4, 5}, TableHeadings -> {None, {"S", "dim", "#weights", "max mult", "#monomials", "#distinct monomial weights", "weights match"}}]
```
**Output:**
```
Out= TableForm[{{2, 10, 10, 1, 10, 10, True}, {3, 20, 20, 1, 20, 20, True}, {4, 35, 35, 1, 35, 35, True}, {5, 56, 56, 1, 56, 56, True}}, TableHeadings -> {None, {"S", "dim", "#weights", "max mult", "#monomials", "#distinct monomial weights", "weights match"}}]
```

## 2026-07-25 17:56:14
**Input:**
```wolfram
Table[With[{g = Counts[patWt /@ ps[{rect[S]}]]}, {S, Max[Values[g]], Floor[S/2] + 1, S + 1, First[Keys[Select[g, # === Max[Values[g]] &]]]}], {S, 1, 5}]
```
**Output:**
```
Out= {{1, 1, 1, 2, {1, 1, 0, 0}}, {2, 2, 2, 3, {1, 1, 1, 1}}, {3, 2, 2, 4, {2, 2, 1, 1}}, {4, 3, 3, 5, {2, 2, 2, 2}}, {5, 3, 3, 6, {3, 3, 2, 2}}}
```
