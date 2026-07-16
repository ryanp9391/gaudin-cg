# AI Action Log — XXX_CG_L2_V2

## 2026-07-16 08:25:23 — ▶️ RUN CELL 1
**In [1]:** `(*Numeric parameters (shared across the whole notebook). NOTE: this assigns global numeric values to θ1,θ2,h,κ1,κ2.↵  χ1,χ2 are the companion twist G's eigenvalues (Section A's G={{χ1,-χ2},{1,0}} has `
**Out:** (no output)


## 2026-07-16 08:25:34
**Input:**
```wolfram
{θ1,θ2,h}
```
**Output:**
```
Out= {1/3, 1/7, 1}
```

## 2026-07-16 08:25:41
**Input:**
```wolfram
{ValueQ[mu],ValueQ[Nn],ValueQ[nn],ValueQ[br],ValueQ[det]}
```
**Output:**
```
Out= {False, False, False, False, False}
```
## 2026-07-16 08:25:59 — 📥 BULK INSERT 6 CELL(S) at positions 35–40
1. [code] (*Section C — FSoV scalar product / orthogonality check. mu/Nn/nn/br/det ported verbatim from↵  Expe
2. [code] (*FSoV measures*)↵ClearAll[mu];↵mu[λ1_,λ2_][1][u_]:=Nn[λ1,λ2][1]/(Product[(u-θ1-h k),{k,0,λ1}]Produc
3. [code] (*Normalization factors in measure*)↵Nn[λ1_,λ2_][1]:=1/(2 π I)/(1-Exp[2 π I/h(θ1-θ2)]) nn[λ1,λ2][1];
4. [code] (*brackets of function f, given by contour integral around all poles in measure*)↵ClearAll[br];↵br[λ
5. [code] ClearAll[det];↵det[λ1_,λ2_][f_]:=1/(f/.u:>θ1)/(f/.u:>θ2)/(θ2-θ1)({{br[λ1,λ2][1][f],br[λ1,λ2][1][u f]
6. [code] (*Port-fidelity check: fsovL2CachedFormula is FSoV_L2.wb's own cached Out[215] for↵  det[1,1][f[u] z


## 2026-07-16 08:25:59
**Input:**
```wolfram
(*Section C — FSoV scalar product / orthogonality check. mu/Nn/nn/br/det ported verbatim from
  Experiments/FSoV_L2.wb (that file's own cells 0-4): mu is the FSoV measure per site (supported
  on that site's θ-tower), br is the bracket pairing (residue sum of f·mu over a site's poles),
  det is the 2x2 Wronskian-type determinant built from br[f] and br[u f] at both sites. For
  f=Q_A(u)Q_B(u), det[λ1,λ2][f] is the SoV scalar product <Q_A,Q_B>: it must vanish for A≠B and
  be nonzero for A=B. See docs/superpowers/specs/2026-07-16-fsov-scalar-product-design.md.*)
```
**Output:**
```
(no output)
```

## 2026-07-16 08:25:59
**Input:**
```wolfram
(*FSoV measures*)
ClearAll[mu];
mu[λ1_,λ2_][1][u_]:=Nn[λ1,λ2][1]/(Product[(u-θ1-h k),{k,0,λ1}]Product[(u-θ2-h k),{k,0,λ2}])(1-Exp[2 π I/h (u-θ2)]);
mu[λ1_,λ2_][2][u_]:=Nn[λ1,λ2][2]/(Product[(u-θ1-h k),{k,0,λ1}]Product[(u-θ2-h k),{k,0,λ2}])(1-Exp[2 π I/h (u-θ1)]);
```
**Output:**
```
(no output)
```

## 2026-07-16 08:26:00
**Input:**
```wolfram
(*Normalization factors in measure*)
Nn[λ1_,λ2_][1]:=1/(2 π I)/(1-Exp[2 π I/h(θ1-θ2)]) nn[λ1,λ2][1];
Nn[λ1_,λ2_][2]:=1/(2 π I)/(1-Exp[2 π I/h(θ2-θ1)]) nn[λ1,λ2][2];

nn[λ1_,λ2_][1]:= (-1)^(λ1+λ2-1)h^(λ1)Gamma[λ1+1] Product[θ2+k h-θ1,{k,0,λ2}]//Simplify;
nn[λ1_,λ2_][2]:= (-1)^(λ1+λ2-1)h^(λ2)Gamma[λ2+1] Product[θ1+k h-θ2,{k,0,λ1}]//Simplify;
```
**Output:**
```
(no output)
```

## 2026-07-16 08:26:00
**Input:**
```wolfram
(*brackets of function f, given by contour integral around all poles in measure*)
ClearAll[br];
br[λ1_,λ2_][1][f_]:=2 π I Plus@@(Residue[f mu[λ1,λ2][1][u],{u,#}]&/@Table[θ1+h k,{k,0,λ1}])//Expand//Simplify;
br[λ1_,λ2_][2][f_]:=2 π I Plus@@(Residue[f mu[λ1,λ2][2][u],{u,#}]&/@Table[θ2+h k,{k,0,λ2}])//Expand//Simplify;
```
**Output:**
```
(no output)
```

## 2026-07-16 08:26:00
**Input:**
```wolfram
ClearAll[det];
det[λ1_,λ2_][f_]:=1/(f/.u:>θ1)/(f/.u:>θ2)/(θ2-θ1)({{br[λ1,λ2][1][f],br[λ1,λ2][1][u f]},{br[λ1,λ2][2][f],br[λ1,λ2][2][u f]}}//Det//Expand//Simplify)
```
**Output:**
```
(no output)
```

## 2026-07-16 08:26:01
**Input:**
```wolfram
(*Port-fidelity check: fsovL2CachedFormula is FSoV_L2.wb's own cached Out[215] for
  det[1,1][f[u] z^(u/h)], transcribed verbatim with θ1,θ2,h,z,f still symbolic (as they were in
  that standalone notebook). Here θ1,θ2,h are already bound (cell 0), so comparing the freshly
  computed det[...] against this formula with the SAME substitution applied is a direct
  transcription check, independent of re-deriving the algebra.*)
ClearAll[fsovL2CachedFormula];
fsovL2CachedFormula[θ1_,θ2_,h_,z_,fn_] := -(z^((θ1+θ2)/h-θ1/h-θ2/h) (z fn[h+θ1] (fn[θ2] (h-θ1+θ2)+z (θ1-θ2) fn[h+θ2])+fn[θ1] (fn[θ2] (θ1-θ2)-z fn[h+θ2] (h+θ1-θ2))))/(fn[θ1] fn[θ2] (θ2-θ1));
portCheckDiff = Simplify[det[1,1][f[u] z^(u/h)] - fsovL2CachedFormula[θ1,θ2,h,z,f]]
```
**Output:**
```
(no output)
```

## 2026-07-16 08:26:12
**Input:**
```wolfram
portCheckDiff
```
**Output:**
```
Out= 0
```
## 2026-07-16 08:26:19 — ▶️ RUN CELL 40
**In [40]:** `(*Port-fidelity check: fsovL2CachedFormula is FSoV_L2.wb's own cached Out[215] for↵  det[1,1][f[u] z^(u/h)], transcribed verbatim with θ1,θ2,h,z,f still symbolic (as they were in↵  that standalone not`
**Out:** (no output)

## 2026-07-16 08:26:57 — ▶️ RUN CELL 40
**In [40]:** `(*Port-fidelity check: fsovL2CachedFormula is FSoV_L2.wb's own cached Out[215] for↵  det[1,1][f[u] z^(u/h)], transcribed verbatim with θ1,θ2,h,z,f still symbolic (as they were in↵  that standalone not`
**Out:** (no output)

## 2026-07-16 08:27:05 — 📥 BULK INSERT 1 CELL(S) at positions 41–41
1. [code] 1+1


## 2026-07-16 08:27:05
**Input:**
```wolfram
1+1
```
**Output:**
```
(no output)
```

## 2026-07-16 08:27:21
**Input:**
```wolfram
$Line
```
**Output:**
```
Out= 246
```

## 2026-07-16 08:27:21
**Input:**
```wolfram
1+1
```
**Output:**
```
Out= 2
```
## 2026-07-16 08:27:46 — 🗑️ DELETE CODE CELL 41
1+1


## 2026-07-16 08:27:52
**Input:**
```wolfram
{θ1,θ2,h}
{mu[1,1][1][u],mu[1,1][2][u]}
Nn[1,1][1]
nn[1,1][1]
br[1,1][1][1]
det[1,1][f[u] z^(u/h)] - fsovL2CachedFormula[θ1,θ2,h,z,f] // Simplify
portCheckDiff
```
**Output:**
```
Out[1]= {1/3, 1/7, 1}
Out[2]= {(((-34*I)/441)*(1 - E^((2*I)*Pi*(-1/7 + u))))/((1 - E^(((8*I)/21)*Pi))*Pi*(-4/3 + u)*(-8/7 + u)*(-1/3 + u)*(-1/7 + u)), (((50*I)/441)*(1 - E^((2*I)*Pi*(-1/3 + u))))/((1 - E^(((-8*I)/21)*Pi))*Pi*(-4/3 + u)*(-8/7 + u)*(-1/3 + u)*(-1/7 + u))}
Out[3]= ((-34*I)/441)/((1 - E^(((8*I)/21)*Pi))*Pi)
Out[4]= 68/441
Out[5]= 42/25
Out[6]= 0
Out[7]= 0
```

## 2026-07-16 08:30:56
**Input:**
```wolfram
{Length[DownValues[TauEigensystem]], Length[DownValues[Q1]], Length[DownValues[det]], Keys[TauEigensystem[2,1]]}
```
**Output:**
```
[message] Keys::invrl: The argument TauEigensystem[2, 1] is not a valid Association or a list of rules.
Out= {0, 0, 0, Keys[TauEigensystem[2, 1]]}
```

## 2026-07-16 08:31:04
**Input:**
```wolfram
{Context[TauEigensystem], Context[Q1], Context[det], Context[dPred], Context[Psi]}
```
**Output:**
```
Out= {"Global`", "Global`", "Global`", "Global`", "Global`"}
```

## 2026-07-16 08:31:23
**Input:**
```wolfram
{ValueQ[TauEigensystem], ValueQ[Q1], ValueQ[det], ValueQ[dPred], ValueQ[mu], ValueQ[br]}
```
**Output:**
```
Out= {False, False, False, False, False, False}
```

## 2026-07-16 08:32:18
**Input:**
```wolfram
{ValueQ[t2], ValueQ[T], ValueQ[L1], ValueQ[ps], ValueQ[dim], ValueQ[JJ], ValueQ[B2], ValueQ[x2], ValueQ[EE]}
```
**Output:**
```
Out= {False, False, False, False, False, False, False, False, False}
```

## 2026-07-16 08:32:31
**Input:**
```wolfram
{Length[DownValues[t2]], Length[DownValues[T]], Length[DownValues[L1]], Length[DownValues[ps]], Length[DownValues[dim]], Length[DownValues[JJ]], Length[DownValues[B2]], Length[DownValues[x2]], Length[DownValues[EE]], Length[DownValues[TauEigensystem]], Length[DownValues[Q1]], Length[DownValues[Psi]], Length[SubValues[det]]}
```
**Output:**
```
Out= {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}
```
## 2026-07-16 08:33:31 — ▶️ RUN CELLS 1–34
Cell 1: ✓ — (no output)
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 8: ✓ — (no output)
Cell 10: ✓ — (no output)
Cell 11: ✓ — (no output)
Cell 12: ✓ — (no output)
Cell 13: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — (no output)
Cell 16: ✓ — (no output)
Cell 17: ✓ — (no output)
Cell 18: ✓ — (no output)
Cell 20: ✓ — (no output)
Cell 21: ✓ — (no output)
Cell 23: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 26: ✓ — (no output)
Cell 28: ✓ — (no output)
Cell 29: ✓ — (no output)
Cell 30: ✓ — (no output)
Cell 31: ✓ — (no output)
Cell 32: ✓ — (no output)
Cell 33: ✓ — (no output)
Cell 34: ✓ — (no output)


## 2026-07-16 08:33:38
**Input:**
```wolfram
{Length[DownValues[TauEigensystem]], Length[DownValues[Q1]], Length[DownValues[t2]]}
```
**Output:**
```
Out= {3, 0, 0}
```

## 2026-07-16 08:33:38
**Input:**
```wolfram
KeySort[TauEigensystem[1,1]][[All,"J"]]//N
```
**Output:**
```
Out= <|{0, 0} -> 0.5664055180884988 + 1.3769645669711723*I, {1, 0} -> -0.778123966451158 + 0.*I, {1, 1} -> 1.1900552523337025 + 0.*I, {2, 0} -> 0.5664055180884988 - 1.3769645669711723*I|>
```

## 2026-07-16 08:33:48
**Input:**
```wolfram
{Length[SubValues[Q1]], Q1[1,1][1,0][2]}
```
**Output:**
```
Out= {1, -0.70923101095801692653722514224657913666`26.84522127451995 + 1.09596916251277206706752728049497901496`27.042725760166295*I}
```
## 2026-07-16 08:33:57 — 📥 BULK INSERT 3 CELL(S) at positions 41–43
1. [code] (*statesOf: the sorted {M,n} key list for a representation, read off TauEigensystem's memoized↵  cac
2. [code] states21 = statesOf[2, 1]
3. [code] G21 = GramMatrix[2, 1];↵Dimensions[G21]


## 2026-07-16 08:33:58
**Input:**
```wolfram
(*statesOf: the sorted {M,n} key list for a representation, read off TauEigensystem's memoized
  cache (built once in Section B) rather than re-deriving dPred here.
  GramMatrix: the full Gram matrix G_{(M,n),(M',n')} = det[λ1,λ2][Q1(M,n)(u) Q1(M',n')(u)] over
  ALL states of the representation at once — this covers same-M degeneracy AND cross-M pairs in
  a single object, per the design.*)
ClearAll[statesOf, GramMatrix];
statesOf[λ1_,λ2_] := SortBy[Keys[TauEigensystem[λ1,λ2]], {First,Last}];
GramMatrix[λ1_,λ2_] := Module[{s = statesOf[λ1,λ2]},
  Table[
    det[λ1,λ2][Q1[λ1,λ2][s1[[1]],s1[[2]]][u] Q1[λ1,λ2][s2[[1]],s2[[2]]][u]],
    {s1, s}, {s2, s}]
];
```
**Output:**
```
(no output)
```

## 2026-07-16 08:33:58
**Input:**
```wolfram
states21 = statesOf[2, 1]
```
**Output:**
```
(no output)
```

## 2026-07-16 08:33:59
**Input:**
```wolfram
G21 = GramMatrix[2, 1];
Dimensions[G21]
```
**Output:**
```
(no output)
```

## 2026-07-16 08:34:04
**Input:**
```wolfram
states21
```
**Output:**
```
Out= {{0, 0}, {1, 0}, {1, 1}, {2, 0}, {2, 1}, {3, 0}}
```

## 2026-07-16 08:34:10
**Input:**
```wolfram
Dimensions[G21]
```
**Output:**
```
Out= {6, 6}
```

## 2026-07-16 08:34:18
**Input:**
```wolfram
{FreeQ[G21, $Failed], FreeQ[G21, _Residue], Head /@ Flatten[G21] // Union // (Complement[#,{Complex,Integer,Real}]&), Diagonal[G21]}
```
**Output:**
```
Out= {True, True, {}, {2.90906428931508326308379330009016616504`28.70086473314406 - 5.80462006356593384304255385623629945047`28.813651008705538*I, 7.3660840301643179118107272726855156597`20.537831121691273 - 0.96804629979446140026588742680878951036`20.541470560591492*I, 2.68702160063169134011892018475885623259`20.583969362984988 - 1.34561532605578582573985646197209785482`20.305173588935705*I, 7.3660840301643179118107272726855165054`23.190597579888355 + 0.96804629979446140026588742680878801124`22.420796492243838*I, 2.68702160063169134011892018475885605506`23.500210129608277 + 1.34561532605578582573985646197209772463`23.182166798943467*I, 2.90906428931508326308379330009016636418`21.699567121386053 + 5.80462006356593384304255385623629952445`21.677262557344527*I}}
```

## 2026-07-16 08:36:25
**Input:**
```wolfram
{ValueQ[G21], ValueQ[states21], If[ValueQ[G21], Dimensions[G21], "N/A"]}
```
**Output:**
```
Out= {True, True, {6, 6}}
```

## 2026-07-16 08:36:29
**Input:**
```wolfram
{states21, FreeQ[G21, $Failed], FreeQ[G21, _Residue]}
```
**Output:**
```
Out= {{{0, 0}, {1, 0}, {1, 1}, {2, 0}, {2, 1}, {3, 0}}, True, True}
```
## 2026-07-16 08:36:40 — 📥 BULK INSERT 1 CELL(S) at positions 44–44
1. [code] (*Verification: G21 should be numerically diagonal. maxOffDiag is the largest |off-diagonal|↵  entry


## 2026-07-16 08:36:40
**Input:**
```wolfram
(*Verification: G21 should be numerically diagonal. maxOffDiag is the largest |off-diagonal|
  entry over all 30 off-diagonal pairs (same-M AND cross-M, since states21 spans all of M=0..3);
  minDiagAbs is the smallest |diagonal| entry, confirming non-degeneracy (a genuine scalar
  product, not one that degenerates to 0 on some state). Threshold 10^-15 matches this notebook's
  existing hard-assert convention (e.g. Section B's TauEigensystem residual check, cell 30) — well
  above the ~10^-20..10^-27 noise floor expected at this notebook's 30-digit working precision,
  and far below any genuine O(1)-scale diagonal entry.*)
FSoVGram::notorthogonal = "Off-diagonal Gram entries do not vanish for (λ1,λ2)=(2,1): max|off-diag|=`1`, min|diag|=`2`.";
FSoVGram::degenerate = "A diagonal Gram entry is suspiciously small for (λ1,λ2)=(2,1) (min|diag|=`1`) — scalar product may be degenerate.";
offDiagVals = Flatten[Table[If[i != j, G21[[i, j]], Nothing], {i, 1, 6}, {j, 1, 6}]];
diagVals = Table[G21[[i, i]], {i, 1, 6}];
maxOffDiag = Max[Abs[N[offDiagVals]]];
minDiagAbs = Min[Abs[N[diagVals]]];
If[maxOffDiag > 10^-15,
  Message[FSoVGram::notorthogonal, maxOffDiag, minDiagAbs];
  Abort[]];
If[minDiagAbs < 10^-15,
  Message[FSoVGram::degenerate, minDiagAbs];
  Abort[]];
<|
  "states" -> states21,
  "max|off-diagonal|" -> maxOffDiag,
  "digits vanished (-Log10)" -> N[-Log10[maxOffDiag]],
  "min|diagonal|" -> minDiagAbs,
  "diagonal values" -> diagVals
|>
```
**Output:**
```
(no output)
```

## 2026-07-16 08:36:46
**Input:**
```wolfram
{maxOffDiag, minDiagAbs, N[-Log10[maxOffDiag]], diagVals, states21}
```
**Output:**
```
Out= {5.065342385207116*^-34, 3.0051232400647923, 33.29539119374288, {2.90906428931508326308379330009016616504`28.70086473314406 - 5.80462006356593384304255385623629945047`28.813651008705538*I, 7.3660840301643179118107272726855156597`20.537831121691273 - 0.96804629979446140026588742680878951036`20.541470560591492*I, 2.68702160063169134011892018475885623259`20.583969362984988 - 1.34561532605578582573985646197209785482`20.305173588935705*I, 7.3660840301643179118107272726855165054`23.190597579888355 + 0.96804629979446140026588742680878801124`22.420796492243838*I, 2.68702160063169134011892018475885605506`23.500210129608277 + 1.34561532605578582573985646197209772463`23.182166798943467*I, 2.90906428931508326308379330009016636418`21.699567121386053 + 5.80462006356593384304255385623629952445`21.677262557344527*I}, {{0, 0}, {1, 0}, {1, 1}, {2, 0}, {2, 1}, {3, 0}}}
```

## 2026-07-16 08:36:52
**Input:**
```wolfram
Quiet[Check[
  <|
    "states" -> states21,
    "max|off-diagonal|" -> maxOffDiag,
    "digits vanished (-Log10)" -> N[-Log10[maxOffDiag]],
    "min|diagonal|" -> minDiagAbs,
    "diagonal values" -> diagVals
  |>,
  "MESSAGES_FIRED", {FSoVGram::notorthogonal, FSoVGram::degenerate}
], {FSoVGram::notorthogonal, FSoVGram::degenerate}]
```
**Output:**
```
Out= <|"states" -> {{0, 0}, {1, 0}, {1, 1}, {2, 0}, {2, 1}, {3, 0}}, "max|off-diagonal|" -> 5.065342385207116*^-34, "digits vanished (-Log10)" -> 33.29539119374288, "min|diagonal|" -> 3.0051232400647923, "diagonal values" -> {2.90906428931508326308379330009016616504`28.70086473314406 - 5.80462006356593384304255385623629945047`28.813651008705538*I, 7.3660840301643179118107272726855156597`20.537831121691273 - 0.96804629979446140026588742680878951036`20.541470560591492*I, 2.68702160063169134011892018475885623259`20.583969362984988 - 1.34561532605578582573985646197209785482`20.305173588935705*I, 7.3660840301643179118107272726855165054`23.190597579888355 + 0.96804629979446140026588742680878801124`22.420796492243838*I, 2.68702160063169134011892018475885605506`23.500210129608277 + 1.34561532605578582573985646197209772463`23.182166798943467*I, 2.90906428931508326308379330009016636418`21.699567121386053 + 5.80462006356593384304255385623629952445`21.677262557344527*I}|>
```
## 2026-07-16 08:37:07 — ▶️ RUN CELL 44
**In [44]:** `(*Verification: G21 should be numerically diagonal. maxOffDiag is the largest |off-diagonal|↵  entry over all 30 off-diagonal pairs (same-M AND cross-M, since states21 spans all of M=0..3);↵  minDiagA`
**Out:** (no output)

## 2026-07-16 08:41:35 — ▶️ RUN CELL 42
**In [42]:** `states21 = statesOf[2, 1]`
**Out:** (no output)


## 2026-07-16 08:41:48
**Input:**
```wolfram
states21
```
**Output:**
```
Out= {{0, 0}, {1, 0}, {1, 1}, {2, 0}, {2, 1}, {3, 0}}
```
## 2026-07-16 08:42:14 — ▶️ RUN CELL 42
**In [42]:** `states21 = statesOf[2, 1]`
**Out:** (no output)

## 2026-07-16 08:42:28 — 📥 BULK INSERT 1 CELL(S) at positions 45–45
1. [code] 1+1


## 2026-07-16 08:42:28
**Input:**
```wolfram
1+1
```
**Output:**
```
(no output)
```
## 2026-07-16 08:42:52 — 🗑️ DELETE CODE CELL 45
1+1


## 2026-07-16 08:43:02
**Input:**
```wolfram
{Dimensions[G21], Max[Abs[N[Flatten[Table[If[i!=j,G21[[i,j]],Nothing],{i,1,6},{j,1,6}]]]]], Min[Abs[N[Table[G21[[i,i]],{i,1,6}]]]]}
```
**Output:**
```
Out= {{6, 6}, 5.065342385207116*^-34, 3.0051232400647923}
```
## 2026-07-16 08:46:40 — 📥 BULK INSERT 1 CELL(S) at positions 45–45
1. [code] 1+1


## 2026-07-16 08:46:40
**Input:**
```wolfram
1+1
```
**Output:**
```
(no output)
```
## 2026-07-16 08:46:52 — 🗑️ DELETE CODE CELL 45
1+1


---

## 2026-07-16 08:59:42 — 🔄 KERNEL RESTART

## 2026-07-16 09:00:13 — 📥 BULK INSERT 1 CELL(S) at positions 45–45
1. [code] 1+1


## 2026-07-16 09:00:13
**Input:**
```wolfram
1+1
```
**Output:**
```
(no output)
```
## 2026-07-16 09:00:25 — 🗑️ DELETE CODE CELL 45
1+1

## 2026-07-16 09:01:39 — 📥 BULK INSERT 1 CELL(S) at positions 45–45
1. [code] 1+1


## 2026-07-16 09:01:40
**Input:**
```wolfram
1+1
```
**Output:**
```
(no output)
```
## 2026-07-16 09:02:16 — 🗑️ DELETE CODE CELL 45
1+1


---

## 2026-07-16 09:03:45 — 🔄 KERNEL RESTART

## 2026-07-16 09:58:22 — 🗑️ DELETE CODE CELL 45
*(empty)*

## 2026-07-16 09:58:22 — 🗑️ DELETE CODE CELL 47
*(empty)*

## 2026-07-16 09:58:59 — 📥 BULK INSERT 1 CELL(S) at positions 46–46
1. [code] 1+1


## 2026-07-16 09:59:00
**Input:**
```wolfram
1+1
```
**Output:**
```
Out[1]= 2
```
## 2026-07-16 10:15:09 — 📥 BULK INSERT 1 CELL(S) at positions 46–46
1. [code] (*Section D — left eigenvectors of the companion-twisted transfer matrix, paired with the↵  existing


## 2026-07-16 10:15:09
**Input:**
```wolfram
(*Section D — left eigenvectors of the companion-twisted transfer matrix, paired with the
  existing right eigenvectors (Psi, Section B) to form <Psi|Psi>, and its ratio against the
  FSoV scalar product det[Q,Q] (Section C). TauEigensystemLeft is an INDEPENDENT construction
  (does not edit TauEigensystem): it diagonalizes Transpose[C0+r C1] for its own generic r, then
  matches each left eigenvector to an existing {M,n} state by comparing eigenvalues (since
  t2(r).v=τ(r) v for the right eigenvector and C2=(κ1+κ2)Id exactly, the target combo-eigenvalue
  for state {M,n} is τ[λ1,λ2][M,n][r]-(κ1+κ2) r^2). Exploratory: no predicted formula for the
  ratio exists yet in paper/xxxCG.tex (still single-site only) — this reports values, it doesn't
  assert against one. See docs/superpowers/specs/2026-07-16-fsov-left-eigenvector-ratio-design.md.*)
```
**Output:**
```
(no output)
```
## 2026-07-16 10:15:29 — 📥 BULK INSERT 3 CELL(S) at positions 47–49
1. [code] ClearAll[TauEigensystemLeft];↵TauEigensystemLeft::rankdef = "Left combination C0+r C1 failed to reso
2. [code] (*smoke test: (λ1,λ2)=(2,1) must produce the same 6 keys as statesOf[2,1] (already verified in↵  Sec
3. [code] (*PsiLeft[λ1,λ2][M,n]: on-demand accessor over the memoized TauEigensystemLeft cache, same↵  bounds/


## 2026-07-16 10:15:30
**Input:**
```wolfram
ClearAll[TauEigensystemLeft];
TauEigensystemLeft::rankdef = "Left combination C0+r C1 failed to resolve the `1`-dim eigenspace for (λ1,λ2)=(`2`,`3`) at every tried r in `4`.";
TauEigensystemLeft::nomatch = "No left eigenvalue matched target for state `1` at (λ1,λ2)=(`2`,`3`) within tolerance (closest distance `4`).";
TauEigensystemLeft::dupmatch = "Left eigenvector index `4` matched more than one state (latest: `1`) at (λ1,λ2)=(`2`,`3`) — combo eigenvalues are not resolving distinctly.";
TauEigensystemLeft::lastzero = "Matched left eigenvector for state `1` at (λ1,λ2)=(`2`,`3`) has last component ~0; cannot normalize.";
TauEigensystemLeft::countmismatch = "TauEigensystemLeft produced `3` states for (λ1,λ2)=(`1`,`2`), expected `4`.";
TauEigensystemLeft[λ1_, λ2_] := TauEigensystemLeft[λ1, λ2] = Module[
   {d, Cs, C0, C1, C2, rCandidatesL, r, valsL, vecsL, states, usedIdx, table},
   d = (λ1 + 1) (λ2 + 1);
   Cs = t2Coeffs[λ1, λ2]; {C0, C1, C2} = Cs;
   rCandidatesL = {17/6, N[Pi, 30] + 1, N[E, 30] + 1, 9/4, 19/6, N[Sqrt[3], 30]};
   r = SelectFirst[rCandidatesL,
     Module[{ee, vv}, {ee, vv} = Eigensystem[Transpose[C0 + # C1]];
       MatrixRank[vv] == d && Length[Union[Chop[ee, 10^-12]]] == d] &,
     Missing["NoGenericR"]];
   If[MissingQ[r], Message[TauEigensystemLeft::rankdef, d, λ1, λ2, rCandidatesL]; Abort[]];
   {valsL, vecsL} = Eigensystem[Transpose[N[C0 + r C1, 30]]];
   states = statesOf[λ1, λ2];
   usedIdx = {};
   table = Association @ Table[
      Module[{Mn = states[[k]], target, dists, idx, w, last},
        target = τ[λ1, λ2][Mn[[1]], Mn[[2]]][r] - (κ1 + κ2) r^2;
        dists = Abs[valsL - target];
        idx = First[Ordering[dists]];
        If[dists[[idx]] > 10^-10,
          Message[TauEigensystemLeft::nomatch, Mn, λ1, λ2, dists[[idx]]]; Abort[]];
        If[MemberQ[usedIdx, idx],
          Message[TauEigensystemLeft::dupmatch, Mn, λ1, λ2, idx]; Abort[]];
        AppendTo[usedIdx, idx];
        w = vecsL[[idx]];
        last = w[[-1]];
        If[Abs[last] < 10^-15 Max[Abs[w]],
          Message[TauEigensystemLeft::lastzero, Mn, λ1, λ2]; Abort[]];
        Mn -> w/last],
      {k, 1, Length[states]}];
   If[Length[table] != d,
     Message[TauEigensystemLeft::countmismatch, λ1, λ2, Length[table], d]; Abort[]];
   table
];
```
**Output:**
```
(no output)
```

## 2026-07-16 10:15:30
**Input:**
```wolfram
(*smoke test: (λ1,λ2)=(2,1) must produce the same 6 keys as statesOf[2,1] (already verified in
  Section C), just via the independent left-eigenvector construction.*)
KeySort[TauEigensystemLeft[2, 1]] // Keys
```
**Output:**
```
2 kernel messages: Set::shape: Lists {C0$2666, C1$2666, C2$2666} and t2Coeffs[2, 1] are not the same shape.
Set::shape: Lists {ee$2673, vv$2673} and Eigensystem[Transpose[C0$2666 + ----------]] are not the same shape.  … [also: Set::shape] result: Out[9]= \text{\$Aborted}
```

## 2026-07-16 10:15:30
**Input:**
```wolfram
(*PsiLeft[λ1,λ2][M,n]: on-demand accessor over the memoized TauEigensystemLeft cache, same
  bounds/spurious-guard pattern as Psi (Section B).*)
ClearAll[PsiLeft];
PsiLeft::bnds = "Invalid indices (M,n)=(`1`,`2`) for (λ1,λ2)=(`3`,`4`): need M an integer in 0..λ1+λ2 and n an integer in 0..M.";
PsiLeft::spurious = "No genuine state at (M,n)=(`1`,`2`) for (λ1,λ2)=(`3`,`4`): weight sector M=`1` has only `5` genuine state(s) (n=0..`6`), so n=`2` does not exist.";
PsiLeft[λ1_, λ2_][M_, n_] := Module[{tab, mult},
  If[! (IntegerQ[M] && IntegerQ[n] && 0 <= M <= λ1 + λ2 && 0 <= n <= M),
    Message[PsiLeft::bnds, M, n, λ1, λ2]; Return[$Failed]];
  tab = TauEigensystemLeft[λ1, λ2];
  mult = dPred[M, λ1, λ2];
  If[n > mult - 1,
    Message[PsiLeft::spurious, M, n, λ1, λ2, mult, mult - 1]; Return[Missing["NoState", {λ1, λ2, M, n}]]];
  tab[{M, n}]];
(*smoke test: a valid state's left eigenvector, plus confirming its last component is exactly 1*)
{PsiLeft[2, 1][1, 0][[-1]], Length[PsiLeft[2, 1][1, 0]]}
```
**Output:**
```
2 kernel messages: Set::shape: Lists {C0$2696, C1$2696, C2$2696} and t2Coeffs[2, 1] are not the same shape.
Set::shape: Lists {ee$2697, vv$2697} and Eigensystem[Transpose[C0$2696 + ----------]] are not the same shape.  … [also: Set::shape] result: Out[14]= \text{\$Aborted}
```

## 2026-07-16 10:15:42
**Input:**
```wolfram
t2Coeffs[2,1]//Dimensions
```
**Output:**
```
Out= {2}
```

## 2026-07-16 10:16:01
**Input:**
```wolfram
ClearAll[t2CoeffsTest]; t2CoeffsTest[λ1_,λ2_]:=Module[{uu,mat,cl},mat=t2[λ1,λ2][1,1][uu];cl=Map[PadRight[CoefficientList[#,uu],3]&,mat,{2}];Table[cl[[All,All,k]],{k,1,3}]]; t2CoeffsTest[2,1]//Dimensions
```
**Output:**
```
[message] Symbol::argx: Symbol called with 0 arguments; 1 argument is expected.
[message] [Symbol::argx ×2 more suppressed — first instance shown above]
[message] General::stop: Further output of Symbol::argx will be suppressed during this calculation.
Out= {3}
```

## 2026-07-16 10:16:11
**Input:**
```wolfram
t2[2,1][1,1][uu]//Dimensions
t2Coeffs[2,1]
```
**Output:**
```
Out[1]= {1}
Out[2]= t2Coeffs[2, 1]
```

## 2026-07-16 10:16:18
**Input:**
```wolfram
t2[2,1][1,1][uu]
```
**Output:**
```
Out= t2[2, 1][1, 1][uu]
```

## 2026-07-16 10:16:31
**Input:**
```wolfram
?t2Coeffs
DownValues[t2Coeffs]
```
**Output:**
```
Out[1]= InformationData[<|"ObjectType" -> "Symbol", "Usage" -> "Global`t2Coeffs", "Documentation" -> None, "OwnValues" -> None, "UpValues" -> None, "DownValues" -> None, "SubValues" -> None, "DefaultValues" -> None, "NValues" -> None, "FormatValues" -> None, "Options" -> None, "Attributes" -> {}, "FullName" -> "Global`t2Coeffs"|>, False]
Out[2]= {}
```

## 2026-07-16 10:16:31
**Input:**
```wolfram
T[2,1][1,1][uu]//Dimensions
```
**Output:**
```
Out= {1}
```
## 2026-07-16 10:17:08 — ▶️ RUN CELLS 1–45
Cell 1: ✓ — Out[22]= \left\{\frac{1}{3},\frac{1}{7},1,0.360439875147226537266677650596+0.932782448593374716484835982071\,\mathrm{i},0.36043987514722653726667765060-0.932782448593374716484835982071\,\mathrm{i},0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i},1\right\}
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 8: ✓ — (no output)
Cell 10: ✓ — (no output)
Cell 11: ✓ — (no output)
Cell 12: ✓ — Out[41]= \left\{0\right\}
Cell 13: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — Out[45]= \left\{0\right\}
Cell 16: ✓ — (no output)
Cell 17: ✓ — Out[48]= \begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{pmatrix}
Cell 18: ✓ — (no output)
Cell 20: ✓ — (no output)
Cell 21: ✓ — (no output)
Cell 23: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 26: ✓ — Out[64]= \left\{0\right\}
Cell 28: ✓ — Out[70]= \left\{\left(u-\frac{4}{3}\right)\,\left(u-\frac{8}{7}\right),\left(u-\frac{1}{3}\right)\,\left(u-\frac{1}{7}\right),\left(0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i}\right)\,u^2-\left(1.06415582186324025288257211128+0.\times 10^{-29}\,\mathrm{i}\right)\,u\right\}
Cell 29: ✓ — Out[73]= \begin{pmatrix}1 & 1\\2 & 2\\3 & 3\end{pmatrix}
Cell 30: ✓ — Out[76]= \text{max|C2-($\kappa$1+$\kappa$2)Id|}\to 0.\times 10^{-29}
Cell 31: ✓ — Out[86]= \langle\!\lvert \left\{0,0\right\}\to 0.5664055180884988 +1.3769645669711723\,\mathrm{i},\left\{1,0\right\}\to -0.778123966451158+0.\,\,\mathrm{i},\left\{1,1\right\}\to 1.1900552523337025 +0.\,\,\mathrm{i},\left\{2,0\right\}\to 0.5664055180884988 -1.3769645669711723\,\mathrm{i}\rvert\!\rangle
Cell 32: ✓ — 2 kernel messages: Psi::spurious: No genuine state at (M,n)=(2,2) for (λ1,λ2)=(2,1): weight sector M=2 has only 2 genuine state(s) (n=0..1), so n=2 does not exist. [also: Psi::spurious] result: Print: spurious-guard test (expect a Psi::spurious message and Missing[...]):[formula] | Out[95]= \left\{\left\{1.000000000000000000000000000+0.\times 10^{-29}\,\mathrm{i},-1.3004494881076436272092358…
⚠️ stopped at Cell 32 — error detected (pass stopOnError:false to continue past errors)

## 2026-07-16 10:17:23 — ▶️ RUN CELLS 33–45
Cell 33: ✓ — Out[102]= \left\{1.03464899078135229594035121+0.381869568380946283899575941\,\mathrm{i},-2.37022278949001677248172619-0.20048970475590097847093779\,\mathrm{i},1.00000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right\}
Cell 34: ✓ — Out[105]= \left\{-0.709231010958016926537225142+1.09596916251277206706752728\,\mathrm{i},-0.02291660900+0.\times 10^{-12}\,\mathrm{i}\right\}
Cell 35: ✓ — (no output)
Cell 36: ✓ — (no output)
Cell 37: ✓ — (no output)
Cell 38: ✓ — (no output)
Cell 39: ✓ — (no output)
Cell 40: ✓ — Out[120]= 0
Cell 41: ✓ — (no output)
Cell 42: ✓ — Out[124]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 43: ✓ — Out[126]= \left\{6,6\right\}
Cell 44: ✓ — Out[127]= \begin{pmatrix}2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i} & 0 & 0 & 0 & 0 & 0\\0 & 7.36608403016431791181-0.968046299794461400266\,\mathrm{i} & 0 & 0 & 0 & 0\\0 & 0 & 2.68702160063169134012-1.3456153260557858257\,\mathrm{i} & 0 & 0 & 0\\0 & 0 & 0 & 7.3660840301643179118107+0.9680462997944614002659\,\mathrm{i} & 0 & 0\\0 & 0 & 0 & 0 & 2.687021600631691340118… [537 chars total]
Cell 45: ✓ — Out[136]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{max|off-diagonal|}\to 5.065342385207116\times 10^{-34},\text{digits vanished (-Log10)}\to 33.29539119374288,\text{min|diagonal|}\to 3.0051232400647923,\text{diagonal values}\to \left\{2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i},7.36608403016431791181-0… [699 chars total]

## 2026-07-16 10:17:30 — ▶️ RUN CELLS 46–49
Cell 46: ✓ — (no output)
Cell 47: ✓ — (no output)
Cell 48: ✓ — Out[144]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 49: ✓ — Out[149]= \left\{1.00000000000000000000000000000,6\right\}

## 2026-07-16 10:18:53 — 📥 BULK INSERT 1 CELL(S) at positions 50–50
1. [code] (*Biorthogonality: a genuine left/right eigenvector pair of a diagonalizable operator with↵  simple


## 2026-07-16 10:18:53
**Input:**
```wolfram
(*Biorthogonality: a genuine left/right eigenvector pair of a diagonalizable operator with
  simple spectrum must satisfy PsiLeft[M,n].Psi[M',n']=0 for {M,n}!={M',n'}. This is the
  correctness check on TauEigensystemLeft's construction itself (independent of Section C's
  det[Q,Q] machinery) — the natural companion to Section C's off-diagonal Gram check. Threshold
  10^-15 matches this notebook's existing hard-assert convention.*)
BiorthoCheck::fail = "Off-diagonal <PsiLeft|Psi> pairing does not vanish for (λ1,λ2)=(2,1): max|off-diag|=`1`.";
offBiortho = Flatten[Table[
   If[i != j, PsiLeft[2, 1] @@ states21[[i]] . (Psi[2, 1] @@ states21[[j]]), Nothing],
   {i, 1, 6}, {j, 1, 6}]];
maxOffBiortho = Max[Abs[N[offBiortho]]];
If[maxOffBiortho > 10^-15,
  Message[BiorthoCheck::fail, maxOffBiortho];
  Abort[]];
<|"max|<PsiLeft|Psi>| off-diagonal" -> maxOffBiortho,
  "digits vanished (-Log10)" -> N[-Log10[maxOffBiortho]]|>
```
**Output:**
```
Out[154]= \langle\!\lvert \text{max|<PsiLeft|Psi>| off-diagonal}\to 0.\,,\text{digits vanished (-Log10)}\to \mathrm{Indeterminate}\rvert\!\rangle
```

## 2026-07-16 10:19:03
**Input:**
```wolfram
{maxOffBiortho, InputForm[maxOffBiortho]}
```
**Output:**
```
Out= {0., InputForm[0.]}
```
## 2026-07-16 10:20:35 — 📥 BULK INSERT 2 CELL(S) at positions 51–52
1. [code] (*braket = <Psi|Psi>: plain bilinear pairing (no conjugation), matching this notebook's existing↵  n
2. [code] (*Final report: braket, G21 diagonal entry, and normRatio for all 6 states of (λ1,λ2)=(2,1).↵  Hard-


## 2026-07-16 10:20:35
**Input:**
```wolfram
(*braket = <Psi|Psi>: plain bilinear pairing (no conjugation), matching this notebook's existing
  non-Hermitian convention. normRatio divides by G21's diagonal entry for the SAME state (G21 is
  Section C's already-verified Gram matrix over states21) — only valid for (λ1,λ2)=(2,1) since
  G21/states21 are fixed globals for that representation, matching this plan's fixed test case.*)
ClearAll[braket, normRatio];
braket[λ1_, λ2_][M_, n_] := PsiLeft[λ1, λ2][M, n] . Psi[λ1, λ2][M, n];
normRatio[λ1_, λ2_][M_, n_] := Module[{i = First[FirstPosition[states21, {M, n}]]},
  braket[λ1, λ2][M, n]/G21[[i, i]]];
(*smoke test*)
{braket[2, 1][1, 0], normRatio[2, 1][1, 0]}
```
**Output:**
```
Out[158]= \left\{7.366084030164317911810727273-0.9680462997944614002658874268\,\mathrm{i},1.00000000000000000000+0.\times 10^{-21}\,\mathrm{i}\right\}
```

## 2026-07-16 10:20:35
**Input:**
```wolfram
(*Final report: braket, G21 diagonal entry, and normRatio for all 6 states of (λ1,λ2)=(2,1).
  Hard-assert every braket is nonzero (non-degenerate pairing) before reporting — a zero braket
  here would mean PsiLeft and Psi are somehow orthogonal even on the diagonal, which should be
  impossible given Task 2's biorthogonality check already passed.*)
BraketReport::degenerate = "braket vanished for state `1` at (λ1,λ2)=(2,1) (value `2`) — pairing may be degenerate.";
reportRows = Table[
   Module[{Mn = states21[[i]], b, g, ratio},
     b = braket[2, 1] @@ Mn;
     g = G21[[i, i]];
     ratio = normRatio[2, 1] @@ Mn;
     If[Abs[N[b]] < 10^-15,
       Message[BraketReport::degenerate, Mn, b];
       Abort[]];
     <|"M" -> Mn[[1]], "n" -> Mn[[2]], "braket" -> N[b, 10], "G21 diagonal" -> N[g, 10],
       "normRatio" -> N[ratio, 10]|>],
   {i, 1, 6}];
reportRows // Column
```
**Output:**
```
Out[161]= \begin{matrix}\langle\!\lvert \text{M}\to 0,\text{n}\to 0,\text{braket}\to 2.909064289-5.804620064\,\mathrm{i},\text{G21 diagonal}\to 2.909064289-5.804620064\,\mathrm{i},\text{normRatio}\to 1.000000000+0.\times 10^{-10}\,\mathrm{i}\rvert\!\rangle \\\langle\!\lvert \text{M}\to 1,\text{n}\to 0,\text{braket}\to 7.366084030-0.968046300\,\mathrm{i},\text{G21 diagonal}\to 7.366084030-0.9680463… [1247 chars total]
```

## 2026-07-16 10:20:41
**Input:**
```wolfram
reportRows
```
**Output:**
```
Out=                                                                                                                                         -10                                                                                                                                            -10                                                                                                                                            -10                                                                                                                                             -10                                                                                                                                            -10                                                                                                                                            -10
{<|M -> 0, n -> 0, braket -> 2.909064289 - 5.804620064 I, G21 diagonal -> 2.909064289 - 5.804620064 I, normRatio -> 1.0000000000 + 0. 10    I|>, <|M -> 1, n -> 0, braket -> 7.366084030 - 0.968046300 I, G21 diagonal -> 7.366084030 - 0.968046300 I, normRatio -> 1.000000000 + 0. 10    I|>, <|M -> 1, n -> 1, braket -> 2.687021601 - 1.345615326 I, G21 diagonal -> 2.687021601 - 1.345615326 I, normRatio -> 1.000000000 + 0. 10    I|>, <|M -> 2, n -> 0, braket -> 7.366084030 + 0.968046300 I, G21 diagonal -> 7.366084030 + 0.968046300 I, normRatio -> 1.0000000000 + 0. 10    I|>, <|M -> 2, n -> 1, braket -> 2.687021601 + 1.345615326 I, G21 diagonal -> 2.687021601 + 1.345615326 I, normRatio -> 1.000000000 + 0. 10    I|>, <|M -> 3, n -> 0, braket -> 2.909064289 + 5.804620064 I, G21 diagonal -> 2.909064289 + 5.804620064 I, normRatio -> 1.000000000 + 0. 10    I|>}
```

## 2026-07-16 10:26:47
**Input:**
```wolfram
1+1
```
**Output:**
```
Out= 2
```

## 2026-07-16 10:26:55
**Input:**
```wolfram
states11 = statesOf[1, 1];
G11 = GramMatrix[1, 1];
states11
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= {{0, 0}, {1, 0}, {1, 1}, {2, 0}}
```

## 2026-07-16 10:27:00
**Input:**
```wolfram
report11 = Table[
  Module[{M = states11[[i, 1]], n = states11[[i, 2]], braket, ratio},
    braket = PsiLeft[1, 1][M, n] . Psi[1, 1][M, n];
    ratio = braket/G11[[i, i]];
    <|"i" -> i, "M" -> M, "n" -> n, "braket" -> N[braket, 10], "ratio" -> N[ratio, 12]|>
  ],
  {i, Length[states11]}
];
report11 // TableForm
```
**Output:**
```
Out= TableForm[{<|"i" -> 1, "M" -> 0, "n" -> 0, "braket" -> 2.57602437156630798107022350951211819847`10.019844242584108 - 2.34025894965916922651994955262289010989`9.978158189805137*I, "ratio" -> 0.99999999999999999999999999999632540872`12.15051499783199 + 0``12.15051499783199*I|>, <|"i" -> 2, "M" -> 1, "n" -> 0, "braket" -> 5.4917989214912120519`10.15051499783199 + 0``9.410800370375338*I, "ratio" -> 1.00000000000000000000000000000197478268`12.15051499783199 + 0``12.15051499783199*I|>, <|"i" -> 3, "M" -> 1, "n" -> 1, "braket" -> 2.54732899300362789298831911184368934496`10.15051499783199 + 0``9.744429959182062*I, "ratio" -> 1.00000000000000000000000000000683224991`12.15051499783199 + 0``12.15051499783199*I|>, <|"i" -> 4, "M" -> 2, "n" -> 0, "braket" -> 2.57602437156630798107022350951211819847`10.019844242584108 + 2.34025894965916922651994955262289010989`9.978158189805137*I, "ratio" -> 0.99999999999999999999999999999632540802`12.15051499783199 + 0``12.15051499783199*I|>}]
```

## 2026-07-16 10:27:03
**Input:**
```wolfram
states12 = statesOf[1, 2];
G12 = GramMatrix[1, 2];
states12
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= {{0, 0}, {1, 0}, {1, 1}, {2, 0}, {2, 1}, {3, 0}}
```

## 2026-07-16 10:27:09
**Input:**
```wolfram
report12 = Table[
  Module[{M = states12[[i, 1]], n = states12[[i, 2]], braket, ratio},
    braket = PsiLeft[1, 2][M, n] . Psi[1, 2][M, n];
    ratio = braket/G12[[i, i]];
    <|"i" -> i, "M" -> M, "n" -> n, "braket" -> N[braket, 10], "ratio" -> N[ratio, 12]|>
  ],
  {i, Length[states12]}
];
report12 // TableForm
```
**Output:**
```
Out= TableForm[{<|"i" -> 1, "M" -> 0, "n" -> 0, "braket" -> 2.9090642893150832631`9.801837232770875 - 5.804620063565933843`10.101857712291416*I, "ratio" -> 0.99999999999999999999999999999607423623`12.15051499783199 + 0``12.15051499783199*I|>, <|"i" -> 2, "M" -> 1, "n" -> 0, "braket" -> 6.7869002179603894678`10.134093454342494 - 1.9022334194371847523`9.581685797487046*I, "ratio" -> 1.00000000000000000000000000000310507803`12.15051499783199 + 0``12.15051499783199*I|>, <|"i" -> 3, "M" -> 1, "n" -> 1, "braket" -> 2.85504176788371417162457950778356869939`10.112777142088998 - 1.24384437218714642763845194655571017999`9.75193072146676*I, "ratio" -> 1.00000000000000000000000000000059941782`12.15051499783199 + 0``12.15051499783199*I|>, <|"i" -> 4, "M" -> 2, "n" -> 0, "braket" -> 6.7869002179603894678`10.134093454342494 + 1.9022334194371847523`9.581685797487046*I, "ratio" -> 0.99999999999999999999999999999608040815`12.15051499783199 + 0``12.15051499783199*I|>, <|"i" -> 5, "M" -> 2, "n" -> 1, "braket" -> 2.85504176788371417162457950779800089846`10.112777142088998 + 1.24384437218714642763845194655171215948`9.75193072146676*I, "ratio" -> 1.00000000000000000000000000000433525262`12.15051499783199 + 0``12.15051499783199*I|>, <|"i" -> 6, "M" -> 3, "n" -> 0, "braket" -> 2.9090642893150832631`9.801837232770875 + 5.804620063565933843`10.101857712291416*I, "ratio" -> 1.00000000000000000000000000000043063875`12.15051499783199 + 0``12.15051499783199*I|>}]
```

## 2026-07-16 10:27:29
**Input:**
```wolfram
te12 = TauEigensystem[1, 2];
te21 = TauEigensystem[2, 1];
Keys[te12] === Keys[te21]
Table[N[te12[k]["J"] - te21[k]["J"], 5], {k, Keys[te12]}]
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= False
Out[4]= {0.0686552143137574357`4.707347828414425 - 0.1776728473511189936`5.120295363558723*I, 0.0686552143137574357`4.707347828414425 + 0.1776728473511189936`5.120295363558723*I, -0.0528247787173994827`4.497877766604335 - 0.2314450220493373158`5.139487924487642*I, 0.1901352073449143541`5.073662230860259 - 0.1239006726529006715`4.887673352373281*I, -0.0528247787173994827`4.497877766604335 + 0.2314450220493373158`5.139487924487642*I, 0.1901352073449143541`5.073662230860259 + 0.1239006726529006715`4.887673352373281*I}
```

## 2026-07-16 10:27:48
**Input:**
```wolfram
G12[[1,1]]
G21fromKernel = {{2.90906428931508326308379330009016616504`28.70086473314406 - 5.80462006356593384304255385623629945047`28.813651008705538*I}};
N[G12[[1,1]] - G21fromKernel[[1,1]], 20]
```
**Output:**
```
Out[1]= 2.90906428931508326308379330009016616504`28.672611766720735 - 5.80462006356593384304255385623629945046`28.776121292226968*I
Out[2]= (suppressed)
Out[3]= 0``27.921725226565176 + 0``27.72967708922445*I
```

## 2026-07-16 10:28:03
**Input:**
```wolfram
states22 = statesOf[2, 2];
G22 = GramMatrix[2, 2];
states22
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= {{0, 0}, {1, 0}, {1, 1}, {2, 0}, {2, 1}, {2, 2}, {3, 0}, {3, 1}, {4, 0}}
```

## 2026-07-16 10:28:07
**Input:**
```wolfram
report22 = Table[
  Module[{M = states22[[i, 1]], n = states22[[i, 2]], braket, ratio},
    braket = PsiLeft[2, 2][M, n] . Psi[2, 2][M, n];
    ratio = braket/G22[[i, i]];
    <|"i" -> i, "M" -> M, "n" -> n, "braket" -> N[braket, 10], "ratio" -> N[ratio, 12]|>
  ],
  {i, Length[states22]}
];
report22 // TableForm
```
**Output:**
```
Out= TableForm[{<|"i" -> 1, "M" -> 0, "n" -> 0, "braket" -> 1.1590896114437539988`9.131390565671127 - 12.0571281801963787844`10.148517430221233*I, "ratio" -> 0.99999999999999999999999999999655966034`12.15051499783199 + 0``12.15051499783199*I|>, <|"i" -> 2, "M" -> 1, "n" -> 0, "braket" -> 7.0247675875539993138`10.019844242584108 - 6.3818399381264284014`9.978158189805137*I, "ratio" -> 0.99999999999999999999999999999719068173`12.15051499783199 + 0``12.15051499783199*I|>, <|"i" -> 3, "M" -> 1, "n" -> 1, "braket" -> 3.29158322996769424627653763404389053449`10.019844242584108 - 2.99032773039959929504928368149130335496`9.978158189805137*I, "ratio" -> 1.00000000000000000000000000004304330883`12.15051499783199 + 0``12.15051499783199*I|>, <|"i" -> 4, "M" -> 2, "n" -> 0, "braket" -> 13.6441172259518654224`10.15051499783199 + 0``29.003230470623212*I, "ratio" -> 1.00000000000000000000000000000079645552`12.15051499783199 + 0``12.15051499783199*I|>, <|"i" -> 5, "M" -> 2, "n" -> 1, "braket" -> 7.0673384996419235617`10.15051499783199 + 0``28.217706837367327*I, "ratio" -> 1.00000000000000000000000000001156647757`12.15051499783199 + 0``12.15051499783199*I|>, <|"i" -> 6, "M" -> 2, "n" -> 2, "braket" -> 3.56430424546688441692640983426976220017`10.15051499783199 + 0``28.02364567690563*I, "ratio" -> 1.00000000000000000000000000001603808513`12.15051499783199 + 0``12.15051499783199*I|>, <|"i" -> 7, "M" -> 3, "n" -> 0, "braket" -> 7.0247675875539993138`10.019844242584108 + 6.3818399381264284014`9.978158189805137*I, "ratio" -> 1.00000000000000000000000000000511590873`12.15051499783199 + 0``12.15051499783199*I|>, <|"i" -> 8, "M" -> 3, "n" -> 1, "braket" -> 3.29158322996769424627653763404389053449`10.019844242584108 + 2.99032773039959929504928368149130335496`9.978158189805137*I, "ratio" -> 1.0000000000000000000000000000121658494`12.15051499783199 + 0``12.15051499783199*I|>, <|"i" -> 9, "M" -> 4, "n" -> 0, "braket" -> 1.1590896114437539988`9.131390565671127 + 12.0571281801963787844`10.148517430221233*I, "ratio" -> 1.00000000000000000000000000000005932109`12.15051499783199 + 0``12.15051499783199*I|>}]
```

## 2026-07-16 10:28:24
**Input:**
```wolfram
Check[{TauEigensystem[1,1];TauEigensystem[1,2];TauEigensystem[2,2];TauEigensystemLeft[1,1];TauEigensystemLeft[1,2];TauEigensystemLeft[2,2];}, "MESSAGES_FIRED", {General::"*"}]
```
**Output:**
```
Out= {Null}
```
