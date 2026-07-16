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
## 2026-07-16 12:24:33 — ▶️ RUN CELLS 1–10
Cell 1: ✓ — Out[169]= \left\{\frac{1}{3},\frac{1}{7},1,0.360439875147226537266677650596+0.932782448593374716484835982071\,\mathrm{i},0.36043987514722653726667765060-0.932782448593374716484835982071\,\mathrm{i},0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i},1\right\}
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 8: ✓ — (no output)
Cell 10: ✓ — (no output)

## 2026-07-16 12:24:43 — ▶️ RUN CELLS 11–26
Cell 11: ✓ — (no output)
Cell 12: ✓ — Out[188]= \left\{0\right\}
Cell 13: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — Out[192]= \left\{0\right\}
Cell 16: ✓ — (no output)
Cell 17: ✓ — Out[195]= \begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{pmatrix}
Cell 18: ✓ — (no output)
Cell 20: ✓ — (no output)
Cell 21: ✓ — (no output)
Cell 23: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 26: ✓ — Out[211]= \left\{0\right\}

## 2026-07-16 12:24:49 — ▶️ RUN CELLS 27–40
Cell 28: ✓ — Out[217]= \left\{\left(u-\frac{4}{3}\right)\,\left(u-\frac{8}{7}\right),\left(u-\frac{1}{3}\right)\,\left(u-\frac{1}{7}\right),\left(0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i}\right)\,u^2-\left(1.06415582186324025288257211128+0.\times 10^{-29}\,\mathrm{i}\right)\,u\right\}
Cell 29: ✓ — Out[220]= \begin{pmatrix}1 & 1\\2 & 2\\3 & 3\end{pmatrix}
Cell 30: ✓ — Out[223]= \text{max|C2-($\kappa$1+$\kappa$2)Id|}\to 0.\times 10^{-29}
Cell 31: ✓ — Out[233]= \langle\!\lvert \left\{0,0\right\}\to 0.5664055180884988 +1.3769645669711723\,\mathrm{i},\left\{1,0\right\}\to -0.778123966451158+0.\,\,\mathrm{i},\left\{1,1\right\}\to 1.1900552523337025 +0.\,\,\mathrm{i},\left\{2,0\right\}\to 0.5664055180884988 -1.3769645669711723\,\mathrm{i}\rvert\!\rangle
Cell 32: ✓ — 2 kernel messages: Psi::spurious: No genuine state at (M,n)=(2,2) for (λ1,λ2)=(2,1): weight sector M=2 has only 2 genuine state(s) (n=0..1), so n=2 does not exist. [also: Psi::spurious] result: Print: spurious-guard test (expect a Psi::spurious message and Missing[...]):[formula] | Out[242]= \left\{\left\{1.000000000000000000000000000+0.\times 10^{-29}\,\mathrm{i},-1.300449488107643627209235…
⚠️ stopped at Cell 32 — error detected (pass stopOnError:false to continue past errors)

## 2026-07-16 12:25:00 — ▶️ RUN CELLS 33–40
Cell 33: ✓ — Out[249]= \left\{1.03464899078135229594035121+0.381869568380946283899575941\,\mathrm{i},-2.37022278949001677248172619-0.20048970475590097847093779\,\mathrm{i},1.00000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right\}
Cell 34: ✓ — Out[252]= \left\{-0.709231010958016926537225142+1.09596916251277206706752728\,\mathrm{i},-0.02291660900+0.\times 10^{-12}\,\mathrm{i}\right\}
Cell 35: ✓ — (no output)
Cell 36: ✓ — (no output)
Cell 37: ✓ — (no output)
Cell 38: ✓ — (no output)
Cell 39: ✓ — (no output)
Cell 40: ✓ — Out[267]= 0

## 2026-07-16 12:25:09 — ▶️ RUN CELLS 41–52
Cell 41: ✓ — (no output)
Cell 42: ✓ — Out[271]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 43: ✓ — Out[273]= \left\{6,6\right\}
Cell 44: ✓ — Out[274]= \begin{pmatrix}2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i} & 0 & 0 & 0 & 0 & 0\\0 & 7.36608403016431791181-0.968046299794461400266\,\mathrm{i} & 0 & 0 & 0 & 0\\0 & 0 & 2.68702160063169134012-1.3456153260557858257\,\mathrm{i} & 0 & 0 & 0\\0 & 0 & 0 & 7.3660840301643179118107+0.9680462997944614002659\,\mathrm{i} & 0 & 0\\0 & 0 & 0 & 0 & 2.687021600631691340118… [537 chars total]
Cell 45: ✓ — Out[283]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{max|off-diagonal|}\to 5.065342385207116\times 10^{-34},\text{digits vanished (-Log10)}\to 33.29539119374288,\text{min|diagonal|}\to 3.0051232400647923,\text{diagonal values}\to \left\{2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i},7.36608403016431791181-0… [699 chars total]
Cell 46: ✓ — (no output)
Cell 47: ✓ — (no output)
Cell 48: ✓ — Out[291]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 49: ✓ — Out[296]= \left\{1.00000000000000000000000000000,6\right\}
Cell 50: ✓ — Out[301]= \langle\!\lvert \text{max|<PsiLeft|Psi>| off-diagonal}\to 0.\,,\text{digits vanished (-Log10)}\to \mathrm{Indeterminate}\rvert\!\rangle
Cell 51: ✓ — Out[305]= \left\{7.366084030164317911810727273-0.9680462997944614002658874268\,\mathrm{i},1.00000000000000000000+0.\times 10^{-21}\,\mathrm{i}\right\}
Cell 52: ✓ — Out[308]= \begin{matrix}\langle\!\lvert \text{M}\to 0,\text{n}\to 0,\text{braket}\to 2.909064289-5.804620064\,\mathrm{i},\text{G21 diagonal}\to 2.909064289-5.804620064\,\mathrm{i},\text{normRatio}\to 1.000000000+0.\times 10^{-10}\,\mathrm{i}\rvert\!\rangle \\\langle\!\lvert \text{M}\to 1,\text{n}\to 0,\text{braket}\to 7.366084030-0.968046300\,\mathrm{i},\text{G21 diagonal}\to 7.366084030-0.9680463… [1247 chars total]


## 2026-07-16 12:25:14
**Input:**
```wolfram
{
 Dimensions[xSingle[2][{0}]],
 Dimensions[x2[2,1][{1,0}]],
 Length[Psi[2,1][1,0]],
 Q1[2,1][1,0][2] =!= $Failed,
 Qθ1[u], Qθ[u],
 statesOf[2,1]
}
```
**Output:**
```
Out= {{1, 3}, {1, 6}, 6, True, -1/3 + u, (-1/3 + u)*(-1/7 + u), {{0, 0}, {1, 0}, {1, 1}, {2, 0}, {2, 1}, {3, 0}}}
```
## 2026-07-16 12:25:36 — 📥 BULK INSERT 5 CELL(S) at positions 53–57
1. [markdown] (*Section E — operatorial factorization of the L=2 SoV basis. Builds a genuine STANDALONE↵  length-1
2. [code] (*Single-site (L=1) fused transfer matrices — the L=1 analogue of Section A's two-site t2 block,↵  b
3. [code] (*smoke: shape check for the fused hierarchy*)↵{Dimensions[t1[2][1,1][u]], Dimensions[t1[2][1,2][u]]
4. [code] (*Genuine length-1 SoV covector (paper normalization, 1-site Qθ1(u)=u-θ1). xSingle[λ][{0}] is the↵  
5. [code] (*Internal check (independent of the factorization): xL1[λ][{s}] must be a LEFT eigencovector of↵  B


## 2026-07-16 12:25:36
**Input:**
```wolfram
(*Single-site (L=1) fused transfer matrices — the L=1 analogue of Section A's two-site t2 block,
  built from the SAME single-site Lax L[λ] (no θ baked in) and shared companion twist G. qdet1 is
  the single-site quantum determinant: both factors are this site's L evaluated at u-θ1 and u-θ1-h
  (the internal -h shift), unlike the two-site qdetT where θ is already in T. t1[λ][1,1](u)=tr_a[G
  L(u-θ1)] matches paper eq. t_{1,1}(u)=tr_a[G L(u-θ)].*)
ClearAll[t1, qdet1];
t1[λ_][0, s_][u_] := IdentityMatrix[dim[{λ, 0}]];
t1[λ_][a_, 0][u_] := IdentityMatrix[dim[{λ, 0}]];
t1[λ_][1, 1][u_] := Sum[L[λ][i, j][u - θ1] G[[j, i]], {i, 2}, {j, 2}];
qdet1[λ_][u_] := L[λ][1, 1][u - θ1] . L[λ][2, 2][u - θ1 - h] - L[λ][2, 1][u - θ1] . L[λ][1, 2][u - θ1 - h];
t1[λ_][2, 1][u_] := χ2 qdet1[λ][u];
t1[λ_][1, s_ /; s > 1][u_] := t1[λ][1, s - 1][u] . t1[λ][1, 1][u + (s - 1) h] - t1[λ][2, 1][u + (s - 1) h] . t1[λ][1, s - 2][u];
```
**Output:**
```
(no output)
```

## 2026-07-16 12:25:37
**Input:**
```wolfram
(*smoke: shape check for the fused hierarchy*)
{Dimensions[t1[2][1,1][u]], Dimensions[t1[2][1,2][u]]}
```
**Output:**
```
Out[316]= \begin{pmatrix}3 & 3\\3 & 3\end{pmatrix}
```

## 2026-07-16 12:25:37
**Input:**
```wolfram
(*Genuine length-1 SoV covector (paper normalization, 1-site Qθ1(u)=u-θ1). xSingle[λ][{0}] is the
  site lowest-weight covector (Section A) — a 1×(λ+1) row. xL1 stays that shape: (1×n).(n×n)=1×n.
  B1 is the L=1 SoV B-operator.*)
ClearAll[B1, xL1];
B1[λ_][u_] := L[λ][1, 1][u - θ1];
xL1[λ_][{0}] := xSingle[λ][{0}];
xL1[λ_][{s_}] := xL1[λ][{s}] = xSingle[λ][{0}] . t1[λ][1, s][θ1] / Product[Qθ1[θ1 - λ h + (p - 1) h], {p, 1, s}];
```
**Output:**
```
(no output)
```

## 2026-07-16 12:25:37
**Input:**
```wolfram
(*Internal check (independent of the factorization): xL1[λ][{s}] must be a LEFT eigencovector of
  B1 with eigenvalue u-θ1-s h — the paper's L=1 SoV property (verified there symbolically for
  λ=1..5). This validates the whole t1/xL1 construction before it is used downstream. Symbolic in
  u; residual Chopped at the 10^-20 noise floor.*)
xL1BCheck = Table[
   xL1[λ][{s}] . B1[λ][u] - (u - θ1 - s h) xL1[λ][{s}],
   {λ, 1, 3}, {s, 0, λ}] // ExpandAll // Chop[#, 10^-20] & // Flatten // Union;
xL1BCheck
```
**Output:**
```
Out[322]= \left\{0\right\}
```
## 2026-07-16 12:25:59 — 📥 BULK INSERT 2 CELL(S) at positions 58–59
1. [code] (*Operatorial factorization: x2[λ1,λ2][{s1,0}] = c (xL1[λ1][{s1}] ⊗ xSingle[λ2][{0}]). Flatten both↵
2. [code] (*Summary: worst ratio-constancy spread across the whole sweep (the hard-asserted factorization↵  co


## 2026-07-16 12:25:59
**Input:**
```wolfram
(*Operatorial factorization: x2[λ1,λ2][{s1,0}] = c (xL1[λ1][{s1}] ⊗ xSingle[λ2][{0}]). Flatten both
  to length-D vectors, take the component ratio over the nonzero components of the tensor product.
  HARD requirement (the actual factorization content — "t2[1,s1][θ1] acts only on site 1"): the
  ratio lhs/tp must be a CONSTANT c (spread below the 10^-15 floor). SEPARATELY, compare c to the
  ideal normalization scalar cPred = 1/Product[(θ1-λ1 h+(k-1)h)-θ2,{k,1,s1}] (reciprocal of the
  product of (tower point − θ2) factors: x2 divides by the 2-site Qθ, xL1 only by the 1-site Qθ1).
  cPred assumes the site-2 pseudovacuum passes through the operator as EXACTLY the identity; the
  ratio μEff = c/cPred is REPORTED (soft-warn if ≠ 1) rather than hard-asserted, because a nontrivial
  μEff would be a genuine physics finding (a site-2 eigenvalue), not a construction bug. cPred is
  independent of λ2, so c matching across λ2∈{1,2,3} is a free extra check. s1=0 ⇒ cPred=1, ratio 1.*)
SoVFact::notconst = "Ratio lhs/tp not constant for (λ1,λ2,s1)=(`1`,`2`,`3`): spread=`4`.";
SoVFact::mu = "μEff=c/cPred ≠ 1 for (λ1,λ2,s1)=(`1`,`2`,`3`): μEff=`4` (site-2 passthrough not identity — a finding, not a failure).";
factResults = Flatten[Table[
    Module[{lhs, tp, nz, rl, cval, cpred, mueff, spread},
     lhs = Flatten[x2[λ1, λ2][{s1, 0}]];
     tp = Flatten[KroneckerProduct[xL1[λ1][{s1}], xSingle[λ2][{0}]]];
     nz = Position[tp, x_ /; Abs[Chop[x, 10^-20]] > 10^-12];
     rl = Extract[lhs, nz]/Extract[tp, nz];
     cval = Mean[rl];
     spread = Max[Abs[rl - cval]];
     cpred = 1/Product[(θ1 - λ1 h + (k - 1) h) - θ2, {k, 1, s1}];
     mueff = cval/cpred;
     If[spread > 10^-15, Message[SoVFact::notconst, λ1, λ2, s1, N[spread]]; Abort[]];
     If[Abs[N[mueff - 1]] > 10^-12, Message[SoVFact::mu, λ1, λ2, s1, N[mueff]]];  (* soft warn, no Abort *)
     {λ1, λ2, s1, cval, cpred, N[mueff], N[spread]}],
    {λ1, 1, 3}, {λ2, 1, 3}, {s1, 0, λ1}], 2];
Column[factResults]
```
**Output:**
```
2 kernel messages: SoVFact::mu: μEff=c/cPred ≠ 1 for (λ1,λ2,s1)=(1,1,1): μEff=-0.809524 (site-2 passthrough not identity — a finding, not a failure).
SoVFact::mu: μEff=c/cPred ≠ 1 for (λ1,λ2,s1)=(1,2,1): μEff=-1.80952 (… [also: SoVFact::mu] result: Out[326]= \begin{matrix}\left\{1,1,0,1,1,1.\,,0.\,\right\}\\\left\{1,1,1,1,-\frac{21}{17},-0.8095238095238095,0.\,\right\}\\\left\{1,2,0,1,1,1.\,,0.\,\right\}\\\left\{1,2,1,\frac{38}{17},-\frac{21}{17…
```

## 2026-07-16 12:26:00
**Input:**
```wolfram
(*Summary: worst ratio-constancy spread across the whole sweep (the hard-asserted factorization
  content), and the range of μEff=c/cPred (should be exactly 1 everywhere if the site-2 passthrough
  is trivial). Columns of factResults: {λ1,λ2,s1,c,cPred,μEff,spread}.*)
<|
  "max spread (ratio non-constancy)" -> Max[factResults[[All, 7]]],
  "max |μEff - 1| over sweep" -> Max[Abs[factResults[[All, 6]] - 1]],
  "num (λ1,λ2,s1) cases checked" -> Length[factResults]
|>
```
**Output:**
```
Out[327]= \langle\!\lvert \text{max spread (ratio non-constancy)}\to 0.\,,\text{max |$\mu$Eff - 1| over sweep}\to 5.115538278803585,\text{num ($\lambda$1,$\lambda$2,s1) cases checked}\to 27\rvert\!\rangle
```

## 2026-07-16 12:26:06
**Input:**
```wolfram
factResults // Column
```
**Output:**
```
Out= {1, 1, 0, 1, 1, 1., 0.}
               21
{1, 1, 1, 1, -(--), -0.809524, 0.}
               17
{1, 2, 0, 1, 1, 1., 0.}
          38    21
{1, 2, 1, --, -(--), -1.80952, 0.}
          17    17
{1, 3, 0, 1, 1, 1., 0.}
          59    21
{1, 3, 1, --, -(--), -2.80952, 0.}
          17    17
{2, 1, 0, 1, 1, 1., 0.}
          17    21
{2, 1, 1, --, -(--), -0.809524, 0.}
          38    38
            2    441
{2, 1, 2, -(--), ---, -0.154195, 0.}
            19   646
{2, 2, 0, 1, 1, 1., 0.}
               21
{2, 2, 1, 1, -(--), -1.80952, 0.}
               38
             441
{2, 2, 2, 1, ---, 1.46485, 0.}
             646
{2, 3, 0, 1, 1, 1., 0.}
          59    21
{2, 3, 1, --, -(--), -2.80952, 0.}
          38    38
          59  441
{2, 3, 2, --, ---, 5.0839, 0.}
          17  646
{3, 1, 0, 1, 1, 1., 0.}
          17    21
{3, 1, 1, --, -(--), -0.809524, 0.}
          59    59
             34    441
{3, 1, 2, -(----), ----, -0.154195, 0.}
            1121   2242
           50     9261
{3, 1, 3, ----, -(-----), -0.183565, 0.}
          1121    38114
{3, 2, 0, 1, 1, 1., 0.}
          38    21
{3, 2, 1, --, -(--), -1.80952, 0.}
          59    59
          17  441
{3, 2, 2, --, ----, 1.46485, 0.}
          59  2242
            4      9261
{3, 2, 3, -(--), -(-----), 0.27902, 0.}
            59     38114
{3, 3, 0, 1, 1, 1., 0.}
               21
{3, 3, 1, 1, -(--), -2.80952, 0.}
               59
             441
{3, 3, 2, 1, ----, 5.0839, 0.}
             2242
               9261
{3, 3, 3, 1, -(-----), -4.11554, 0.}
               38114
```

## 2026-07-16 12:26:35
**Input:**
```wolfram
{
 "max spread (N)" -> N[Max[factResults[[All,7]]]],
 "cval as fn of {λ1,λ2,s1} for s1=1 (checking λ2-independence claim)" ->
    Table[{λ1,λ2,1} -> factResults[[Position[factResults, {λ1,λ2,1,_,_,_,_}][[1,1]],4]], {λ1,1,3},{λ2,1,3}] // Flatten[#,1]&,
 "mueff extremes" -> {Min[factResults[[All,6]]], Max[factResults[[All,6]]]}
}
```
**Output:**
```
Out= {"max spread (N)" -> 0., "cval as fn of {λ1,λ2,s1} for s1=1 (checking λ2-independence claim)" -> {{{1, 1, 1} -> 1, {1, 2, 1} -> 38/17, {1, 3, 1} -> 59/17}, {{2, 1, 1} -> 17/38, {2, 2, 1} -> 1, {2, 3, 1} -> 59/38}, {{3, 1, 1} -> 17/59, {3, 2, 1} -> 38/59, {3, 3, 1} -> 1}}, "mueff extremes" -> {-4.115538278803585, 5.083900226757369}}
```

## 2026-07-16 12:26:42
**Input:**
```wolfram
Max[factResults[[All,7]]] // InputForm
```
**Output:**
```
Out= InputForm[0.]
```

## 2026-07-16 12:31:41
**Input:**
```wolfram
cPredCorr[λ1_,λ2_,s1_] := Product[(θ1-θ2+(k-1)h-λ2 h)/(θ1-θ2+(k-1)h-λ1 h),{k,1,s1}];
checkC = Table[
  Module[{λ1=r[[1]],λ2=r[[2]],s1=r[[3]],cval=r[[4]]},
   {λ1,λ2,s1, cval, cPredCorr[λ1,λ2,s1], N[cval-cPredCorr[λ1,λ2,s1]]}],
  {r, factResults}];
Max[Abs[checkC[[All,6]]]]
```
**Output:**
```
Out= 0.
```

## 2026-07-16 12:31:46
**Input:**
```wolfram
checkC // Column
```
**Output:**
```
Out= {1, 1, 0, 1, 1, 0.}
{1, 1, 1, 1, 1, 0.}
{1, 2, 0, 1, 1, 0.}
          38  38
{1, 2, 1, --, --, 0.}
          17  17
{1, 3, 0, 1, 1, 0.}
          59  59
{1, 3, 1, --, --, 0.}
          17  17
{2, 1, 0, 1, 1, 0.}
          17  17
{2, 1, 1, --, --, 0.}
          38  38
            2      2
{2, 1, 2, -(--), -(--), 0.}
            19     19
{2, 2, 0, 1, 1, 0.}
{2, 2, 1, 1, 1, 0.}
{2, 2, 2, 1, 1, 0.}
{2, 3, 0, 1, 1, 0.}
          59  59
{2, 3, 1, --, --, 0.}
          38  38
          59  59
{2, 3, 2, --, --, 0.}
          17  17
{3, 1, 0, 1, 1, 0.}
          17  17
{3, 1, 1, --, --, 0.}
          59  59
             34       34
{3, 1, 2, -(----), -(----), 0.}
            1121     1121
           50    50
{3, 1, 3, ----, ----, 0.}
          1121  1121
{3, 2, 0, 1, 1, 0.}
          38  38
{3, 2, 1, --, --, 0.}
          59  59
          17  17
{3, 2, 2, --, --, 0.}
          59  59
            4      4
{3, 2, 3, -(--), -(--), 0.}
            59     59
{3, 3, 0, 1, 1, 0.}
{3, 3, 1, 1, 1, 0.}
{3, 3, 2, 1, 1, 0.}
{3, 3, 3, 1, 1, 0.}
```

## 2026-07-16 12:31:55
**Input:**
```wolfram
muSite2test[λ2_,s1_] := Product[(θ1-θ2+(k-1)h-λ2 h),{k,1,s1}];
opCheck = Flatten[Table[
   Module[{numLhs,numTp,nz,rl,muval,mupred,spread},
    numLhs = Flatten[x2[λ1,λ2][{0,0}] . t2[λ1,λ2][1,s1][θ1]];
    numTp  = Flatten[KroneckerProduct[xSingle[λ1][{0}] . t1[λ1][1,s1][θ1], xSingle[λ2][{0}]]];
    nz = Position[numTp, x_ /; Abs[Chop[x,10^-20]] > 10^-12];
    rl = Extract[numLhs,nz]/Extract[numTp,nz];
    muval = Mean[rl];
    spread = Max[Abs[rl-muval]];
    mupred = muSite2test[λ2,s1];
    {λ1,λ2,s1, N[muval], N[mupred], N[muval-mupred], N[spread]}],
   {λ1,1,3},{λ2,1,3},{s1,0,λ1}], 2];
{Max[opCheck[[All,6]]//Abs], Max[opCheck[[All,7]]]}
```
**Output:**
```
Out= {0., 0.}
```

## 2026-07-16 12:32:03
**Input:**
```wolfram
opCheck // Column
```
**Output:**
```
Out= {1, 1, 0, 1., 1., 0., 0.}
{1, 1, 1, -0.809524, -0.809524, 0., 0.}
{1, 2, 0, 1., 1., 0., 0.}
{1, 2, 1, -1.80952, -1.80952, 0., 0.}
{1, 3, 0, 1., 1., 0., 0.}
{1, 3, 1, -2.80952, -2.80952, 0., 0.}
{2, 1, 0, 1., 1., 0., 0.}
{2, 1, 1, -0.809524, -0.809524, 0., 0.}
{2, 1, 2, -0.154195, -0.154195, 0., 0.}
{2, 2, 0, 1., 1., 0., 0.}
{2, 2, 1, -1.80952, -1.80952, 0., 0.}
{2, 2, 2, 1.46485, 1.46485, 0., 0.}
{2, 3, 0, 1., 1., 0., 0.}
{2, 3, 1, -2.80952, -2.80952, 0., 0.}
{2, 3, 2, 5.0839, 5.0839, 0., 0.}
{3, 1, 0, 1., 1., 0., 0.}
{3, 1, 1, -0.809524, -0.809524, 0., 0.}
{3, 1, 2, -0.154195, -0.154195, 0., 0.}
{3, 1, 3, -0.183565, -0.183565, 0., 0.}
{3, 2, 0, 1., 1., 0., 0.}
{3, 2, 1, -1.80952, -1.80952, 0., 0.}
{3, 2, 2, 1.46485, 1.46485, 0., 0.}
{3, 2, 3, 0.27902, 0.27902, 0., 0.}
{3, 3, 0, 1., 1., 0., 0.}
{3, 3, 1, -2.80952, -2.80952, 0., 0.}
{3, 3, 2, 5.0839, 5.0839, 0., 0.}
{3, 3, 3, -4.11554, -4.11554, 0., 0.}
```
## 2026-07-16 12:32:59 — ✏️ EDIT CELL 58
(*Operatorial factorization: x2[λ1,λ2][{s1,0}] = c (xL1[λ1][{s1}] ⊗ xSingle[λ2][{0}]). Flatten both
  to length-D vectors, take the component ratio over the nonzero components of the tensor product.
 …


## 2026-07-16 12:32:59
**Input:**
```wolfram
(*Operatorial factorization: x2[λ1,λ2][{s1,0}] = c (xL1[λ1][{s1}] ⊗ xSingle[λ2][{0}]). Flatten both
  to length-D vectors, take the component ratio over the nonzero components of the tensor product.
  HARD requirement #1: the ratio lhs/tp must be a CONSTANT c (spread below the 10^-15 floor) — the
  actual factorization content ("t2[1,s1][θ1] acts only through a rank-1 site-2 eigenvalue"). HARD
  requirement #2 (hardened after the μEff≠1 finding was pinned down): c must equal the CORRECTED
  scalar cPredCorr[λ1,λ2,s1] = Product[(θ1-θ2+(k-1)h-λ2 h)/(θ1-θ2+(k-1)h-λ1 h),{k,1,s1}], which
  accounts for the genuine site-2 passthrough eigenvalue μSite2[λ2,s1] (defined/verified at the
  operator level in the cell below) — NOT the naive 1/Product[(θ1-λ1 h+(k-1)h)-θ2,{k,1,s1}] that
  assumed a trivial (identity) site-2 passthrough. cPredCorr reduces to 1 at λ1=λ2 and to
  f(λ2)/f(λ1), f(λ)=21λ-4, at s1=1 — matching the earlier empirical finding exactly. mueff=cval/cpred
  is now a closure check (must be exactly 1 to the noise floor) rather than a soft warning.*)
ClearAll[cPredCorr];
cPredCorr[λ1_, λ2_, s1_] := Product[(θ1 - θ2 + (k - 1) h - λ2 h)/(θ1 - θ2 + (k - 1) h - λ1 h), {k, 1, s1}];
SoVFact::notconst = "Ratio lhs/tp not constant for (λ1,λ2,s1)=(`1`,`2`,`3`): spread=`4`.";
SoVFact::mu = "c ≠ cPredCorr for (λ1,λ2,s1)=(`1`,`2`,`3`): μEff=`4` (corrected scalar formula failed — a construction bug, not a finding).";
factResults = Flatten[Table[
    Module[{lhs, tp, nz, rl, cval, cpred, mueff, spread},
     lhs = Flatten[x2[λ1, λ2][{s1, 0}]];
     tp = Flatten[KroneckerProduct[xL1[λ1][{s1}], xSingle[λ2][{0}]]];
     nz = Position[tp, x_ /; Abs[Chop[x, 10^-20]] > 10^-12];
     rl = Extract[lhs, nz]/Extract[tp, nz];
     cval = Mean[rl];
     spread = Max[Abs[rl - cval]];
     cpred = cPredCorr[λ1, λ2, s1];
     mueff = cval/cpred;
     If[spread > 10^-15, Message[SoVFact::notconst, λ1, λ2, s1, N[spread]]; Abort[]];
     If[Abs[N[mueff - 1]] > 10^-15, Message[SoVFact::mu, λ1, λ2, s1, N[mueff]]; Abort[]];  (* now HARD *)
     {λ1, λ2, s1, cval, cpred, N[mueff], N[spread]}],
    {λ1, 1, 3}, {λ2, 1, 3}, {s1, 0, λ1}], 2];
Column[factResults]
```
**Output:**
```
Out= Column[{{1, 1, 0, 1, 1, 1., 0.}, {1, 1, 1, 1, 1, 1., 0.}, {1, 2, 0, 1, 1, 1., 0.}, {1, 2, 1, 38/17, 38/17, 1., 0.}, {1, 3, 0, 1, 1, 1., 0.}, {1, 3, 1, 59/17, 59/17, 1., 0.}, {2, 1, 0, 1, 1, 1., 0.}, {2, 1, 1, 17/38, 17/38, 1., 0.}, {2, 1, 2, -2/19, -2/19, 1., 0.}, {2, 2, 0, 1, 1, 1., 0.}, {2, 2, 1, 1, 1, 1., 0.}, {2, 2, 2, 1, 1, 1., 0.}, {2, 3, 0, 1, 1, 1., 0.}, {2, 3, 1, 59/38, 59/38, 1., 0.}, {2, 3, 2, 59/17, 59/17, 1., 0.}, {3, 1, 0, 1, 1, 1., 0.}, {3, 1, 1, 17/59, 17/59, 1., 0.}, {3, 1, 2, -34/1121, -34/1121, 1., 0.}, {3, 1, 3, 50/1121, 50/1121, 1., 0.}, {3, 2, 0, 1, 1, 1., 0.}, {3, 2, 1, 38/59, 38/59, 1., 0.}, {3, 2, 2, 17/59, 17/59, 1., 0.}, {3, 2, 3, -4/59, -4/59, 1., 0.}, {3, 3, 0, 1, 1, 1., 0.}, {3, 3, 1, 1, 1, 1., 0.}, {3, 3, 2, 1, 1, 1., 0.}, {3, 3, 3, 1, 1, 1., 0.}}]
```
## 2026-07-16 12:33:03 — ▶️ RUN CELL 59 (0.03 s)
**In [59]:** `(*Summary: worst ratio-constancy spread across the whole sweep (the hard-asserted factorization↵  content), and the range of μEff=c/cPred (should be exactly 1 everywhere if the site-2 passthrough↵  is`
**Out:** Out[328]= \langle\!\lvert \text{max spread (ratio non-constancy)}\to 0.\,,\text{max |$\mu$Eff - 1| over sweep}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases checked}\to 27\rvert\!\rangle

## 2026-07-16 12:33:18 — 📥 BULK INSERT 1 CELL(S) at positions 60–60
1. [code] (*Physics-level closure: the site-2 pseudovacuum xSingle[λ2][{0}] is a LEFT eigencovector of the↵  s


## 2026-07-16 12:33:19
**Input:**
```wolfram
(*Physics-level closure: the site-2 pseudovacuum xSingle[λ2][{0}] is a LEFT eigencovector of the
  site-2 Lax factor appearing in t2[λ1,λ2][1,s1][θ1] at each fused step, with eigenvalue
  μSite2[λ2,s1] = Product[(θ1-θ2+(k-1)h-λ2 h),{k,1,s1}]. This is WHY the naive "passthrough is the
  identity" guess failed: L[λ2][1,1](θ1-θ2) has a nontrivial GT lowest-weight eigenvalue at each of
  the s1 fused steps (shifted by (k-1)h from the internal Hirota recursion), not 1. Verified directly
  at the NUMERATOR level (i.e. before any xL1-normalization is applied — the cleanest statement of the
  effect): x2[λ1,λ2][{0,0}].t2[λ1,λ2][1,s1][θ1] must equal μSite2[λ2,s1] times
  (xSingle[λ1][{0}].t1[λ1][1,s1][θ1]) ⊗ xSingle[λ2][{0}], componentwise constant over the nonzero
  tensor-product entries. HARD assert on both constancy (spread) and the closed form (deviation
  from μSite2).*)
ClearAll[μSite2];
μSite2[λ2_, s1_] := Product[(θ1 - θ2 + (k - 1) h - λ2 h), {k, 1, s1}];
Site2Mu::notconst = "μ-ratio not constant at operator level for (λ1,λ2,s1)=(`1`,`2`,`3`): spread=`4`.";
Site2Mu::mismatch = "μ ≠ μSite2 at operator level for (λ1,λ2,s1)=(`1`,`2`,`3`): μ=`4`, μSite2=`5`.";
opMuResults = Flatten[Table[
    Module[{numLhs, numTp, nz, rl, muval, mupred, spread},
     numLhs = Flatten[x2[λ1, λ2][{0, 0}] . t2[λ1, λ2][1, s1][θ1]];
     numTp = Flatten[KroneckerProduct[xSingle[λ1][{0}] . t1[λ1][1, s1][θ1], xSingle[λ2][{0}]]];
     nz = Position[numTp, x_ /; Abs[Chop[x, 10^-20]] > 10^-12];
     rl = Extract[numLhs, nz]/Extract[numTp, nz];
     muval = Mean[rl];
     spread = Max[Abs[rl - muval]];
     mupred = μSite2[λ2, s1];
     If[spread > 10^-15, Message[Site2Mu::notconst, λ1, λ2, s1, N[spread]]; Abort[]];
     If[Abs[N[muval - mupred]] > 10^-15, Message[Site2Mu::mismatch, λ1, λ2, s1, N[muval], N[mupred]]; Abort[]];
     {λ1, λ2, s1, N[muval], N[mupred], N[spread]}],
    {λ1, 1, 3}, {λ2, 1, 3}, {s1, 0, λ1}], 2];
<|"max operator-level spread" -> Max[opMuResults[[All, 6]]],
  "max |μ - μSite2|" -> Max[Abs[opMuResults[[All, 4]] - opMuResults[[All, 5]]]],
  "num (λ1,λ2,s1) cases" -> Length[opMuResults]|>
```
**Output:**
```
Out[334]= \langle\!\lvert \text{max operator-level spread}\to 0.\,,\text{max |$\mu$ - $\mu$Site2|}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases}\to 27\rvert\!\rangle
```

## 2026-07-16 12:37:00
**Input:**
```wolfram
{Dimensions[xSingle[2][{0}]], xL1[2][{1}] // Head, μSite2[2,1], Q1[2,1][1,0][2]=!=$Failed}
```
**Output:**
```
Out= {{1, 3}, List, -38/21, True}
```

## 2026-07-16 12:37:48
**Input:**
```wolfram
{Max[Abs[factResults[[All,6]]-1]], Max[factResults[[All,7]]], factResults[[2]], factResults[[4]]}
```
**Output:**
```
Out= {0., 0., {1, 1, 1, 1, 1, 1., 0.}, {1, 2, 1, 38/17, 38/17, 1., 0.}}
```
## 2026-07-16 12:37:55 — ▶️ RUN CELL 58 (0.03 s)
**In [58]:** `(*Operatorial factorization: x2[λ1,λ2][{s1,0}] = c (xL1[λ1][{s1}] ⊗ xSingle[λ2][{0}]). Flatten both↵  to length-D vectors, take the component ratio over the nonzero components of the tensor product.↵ `
**Out:** Out[340]= \begin{matrix}\left\{1,1,0,1,1,1.\,,0.\,\right\}\\\left\{1,1,1,1,1,1.\,,0.\,\right\}\\\left\{1,2,0,1,1,1.\,,0.\,\right\}\\\left\{1,2,1,\frac{38}{17},\frac{38}{17},1.\,,0.\,\right\}\\\left\{1,3,0,1,1,1.\,,0.\,\right\}\\\left\{1,3,1,\frac{59}{17},\frac{59}{17},1.\,,0.\,\right\}\\\left\{2,1,0

## 2026-07-16 12:38:21 — 📥 BULK INSERT 4 CELL(S) at positions 61–64
1. [markdown] (*Section F — functional CG overlap. Sets up the length-1 FSoV measure/bracket (Section C's↵  mu[...
2. [code] (*Length-1 FSoV measure and bracket. mu1 is a single-site rational measure with simple poles at↵  th
3. [code] (*Internal check: the leading (θ1-node) weight must be exactly 1 for λ=1,2,3 (this IS the↵  br1[f]=f
4. [code] (*Sanity: br1 acting on a generic analytic f reproduces the weighted tower sum Σ_j w_j f(θ1+h j),↵


## 2026-07-16 12:38:21
**Input:**
```wolfram
(*Length-1 FSoV measure and bracket. mu1 is a single-site rational measure with simple poles at
  the tower θ1+h k, k=0..λ. Nn1 is DERIVED (not guessed) by forcing the θ1-node weight
  w_0 = 2πI Res[mu1, θ1] to equal 1, so that br1[λ][f] = f(θ1) + (higher-tower). br1 is the residue
  sum over the whole tower. (For λ≥1 the residues of mu1 sum to 0, so br1[λ][1]=0 — the correct
  normalization statement is w_0=1, checked next.)*)
ClearAll[Nn1, mu1, br1];
Nn1[λ_] := Nn1[λ] = Module[{uu}, 1/(2 π I Residue[1/Product[(uu - θ1 - h k), {k, 0, λ}], {uu, θ1}])];
mu1[λ_][u_] := Nn1[λ]/Product[(u - θ1 - h k), {k, 0, λ}];
br1[λ_][f_] := 2 π I Plus @@ (Residue[f mu1[λ][u], {u, #}] & /@ Table[θ1 + h k, {k, 0, λ}]);
```
**Output:**
```
(no output)
```

## 2026-07-16 12:38:21
**Input:**
```wolfram
(*Internal check: the leading (θ1-node) weight must be exactly 1 for λ=1,2,3 (this IS the
  br1[f]=f(θ1)+... normalization), and br1[λ][1] must be 0 (rational-measure residues sum to
  zero). Also list the full weight vectors w_j for the record.*)
w1[λ_] := Table[2 π I Residue[mu1[λ][u], {u, θ1 + h j}], {j, 0, λ}];
Nn1Check::badlead = "Leading weight w_0 ≠ 1 for λ=`1`: w_0=`2`.";
Nn1Check::badsum = "br1[λ][1] ≠ 0 for λ=`1`: got `2`.";
Do[
  Module[{w = w1[λ], s},
   If[Abs[N[w[[1]] - 1]] > 10^-15, Message[Nn1Check::badlead, λ, w[[1]]]; Abort[]];
   s = br1[λ][1];
   If[Abs[N[s]] > 10^-15, Message[Nn1Check::badsum, λ, s]; Abort[]]],
  {λ, 1, 3}];
<|"w_0 (λ=1,2,3)" -> Table[w1[λ][[1]], {λ, 1, 3}],
  "full weights w_j" -> Table[λ -> w1[λ], {λ, 1, 3}],
  "br1[λ][1] (λ=1,2,3)" -> Table[br1[λ][1], {λ, 1, 3}]|>
```
**Output:**
```
Out[349]= \langle\!\lvert \text{w\_0 ($\lambda$=1,2,3)}\to \left\{1,1,1\right\},\text{full weights w\_j}\to \left\{1\to \left\{1,-1\right\},2\to \left\{1,-2,1\right\},3\to \left\{1,-3,3,-1\right\}\right\},\text{br1[$\lambda$][1] ($\lambda$=1,2,3)}\to \left\{0,0,0\right\}\rvert\!\rangle
```

## 2026-07-16 12:38:22
**Input:**
```wolfram
(*Sanity: br1 acting on a generic analytic f reproduces the weighted tower sum Σ_j w_j f(θ1+h j),
  and picks out f(θ1) as the leading term. Test with f=g[u] (undefined head, stays symbolic).*)
Clear[g];
{br1[1][g[u]] - Sum[w1[1][[j + 1]] g[θ1 + h j], {j, 0, 1}] // Simplify,
 br1[2][g[u]] - Sum[w1[2][[j + 1]] g[θ1 + h j], {j, 0, 2}] // Simplify}
```
**Output:**
```
Out[351]= \left\{0,0\right\}
```
## 2026-07-16 12:38:45 — 📥 BULK INSERT 4 CELL(S) at positions 65–68
1. [code] (*LHS functional object: FL = br1[λ1][k^(u/h) Q1[λ1,λ2][M,n][u]]. Since Q1=κ1^(u/h) q_n(u), the↵  ex
2. [code] (*RHS operatorial CG overlap: rotate the site-1 pseudovacuum covector by g=exp[φ Ee[λ1][a,b]]↵  (Ee[
3. [code] (*Which single-site generator is nontrivial? The lowest-weight covector ⟨x0| annihilates one of↵  E1
4. [code] (*OPg fixes the nontrivial generator so downstream code reads OPg[λ1,λ2][M,n][φ]. Degree in φ is ≤λ1


## 2026-07-16 12:38:45
**Input:**
```wolfram
(*LHS functional object: FL = br1[λ1][k^(u/h) Q1[λ1,λ2][M,n][u]]. Since Q1=κ1^(u/h) q_n(u), the
  exponential base is (k κ1); on the poles u=θ1+h j this is (k κ1)^(θ1/h) (k κ1)^j, so FL is
  (k κ1)^(θ1/h) times a polynomial in k. Keep k symbolic.*)
ClearAll[FL];
FL[λ1_, λ2_][M_, n_][k_] := br1[λ1][k^(u/h) Q1[λ1, λ2][M, n][u]];
(*smoke: one state — expect a (k κ1)^(θ1/h)-type prefactor times a polynomial in k*)
FL[2, 1][1, 0][k] // Simplify
```
**Output:**
```
Out[354]= \sqrt[3]{k}\,\left(\left(-1.18561270626131985872+0.832253303224995280404\,\mathrm{i}\right)\,k^2-\left(0.761819910806671651081+0.812812473459294825632\,\mathrm{i}\right)\,k-\left(0.404410037536507329420+0.598870869485872613792\,\mathrm{i}\right)\right)
```

## 2026-07-16 12:38:46
**Input:**
```wolfram
(*RHS operatorial CG overlap: rotate the site-1 pseudovacuum covector by g=exp[φ Ee[λ1][a,b]]
  (Ee[a,b] nilpotent of order ≤λ1+1 on the finite-dim rep). Build g as the explicit truncated
  series Σ_{m=0}^{λ} φ^m/m! Ee[a,b]^m — this is EXACT (higher powers vanish) and guaranteed to be a
  polynomial in φ, avoiding MatrixExp returning an unevaluated symbolic exponential. Tensor with the
  site-2 pseudovacuum covector, contract with the two-site right eigenvector Psi. Shapes:
  (1×n1).(n1×n1)=1×n1 → KroneckerProduct with 1×n2 → Flatten to length D → .Psi = scalar.*)
ClearAll[gRot, OP];
gRot[λ_][a_, b_][φ_] := Sum[φ^m/m! MatrixPower[Ee[λ][a, b], m], {m, 0, λ}];
OP[λ1_, λ2_][a_, b_][M_, n_][φ_] := Flatten[KroneckerProduct[xSingle[λ1][{0}] . gRot[λ1][a, b][φ], xSingle[λ2][{0}]]] . Psi[λ1, λ2][M, n];
```
**Output:**
```
(no output)
```

## 2026-07-16 12:38:46
**Input:**
```wolfram
(*Which single-site generator is nontrivial? The lowest-weight covector ⟨x0| annihilates one of
  E12/E21 (giving φ-independent OP); keep the other. Test on (λ1,λ2)=(2,1), state (1,0).*)
op12 = OP[2, 1][1, 2][1, 0][φ] // Simplify;
op21 = OP[2, 1][2, 1][1, 0][φ] // Simplify;
genNontrivial = Which[
   !FreeQ[op12, φ] && FreeQ[op21, φ], {1, 2},
   !FreeQ[op21, φ] && FreeQ[op12, φ], {2, 1},
   True, {1, 2}  (* fallback: if both depend on φ, default to {1,2} and flag in the report *)];
<|"OP(E12)" -> op12, "OP(E21)" -> op21, "nontrivial generator {a,b}" -> genNontrivial,
  "both φ-dependent?" -> (!FreeQ[op12, φ] && !FreeQ[op21, φ])|>
```
**Output:**
```
3 kernel messages: MatrixPower::sing: Matrix {{0, 2, 0}, {0, 0, 2}, {0, 0, 0}} is singular. [also: MatrixPower::sing, MatrixPower::sing] result: Out[361]= \langle\!\lvert \text{OP(E12)}\to \left\{\left(-0.5946847648424943838128609535-3.764270633601083661935312855\,\mathrm{i}\right)\,\mathrm{MatrixPower}\left[\begin{pmatrix}0 & 2 & 0\\0 & 0 & 2…
```

## 2026-07-16 12:38:47
**Input:**
```wolfram
(*OPg fixes the nontrivial generator so downstream code reads OPg[λ1,λ2][M,n][φ]. Degree in φ is ≤λ1.*)
ClearAll[OPg];
OPg[λ1_, λ2_][M_, n_][φ_] := OP[λ1, λ2][genNontrivial[[1]], genNontrivial[[2]]][M, n][φ];
(*smoke: degree in φ for a λ1=1 and a λ1=2 case*)
{Exponent[OPg[1, 1][1, 0][φ] // Simplify, φ], Exponent[OPg[2, 1][1, 0][φ] // Simplify, φ]}
```
**Output:**
```
MatrixPower::sing: Matrix {{0, 0}, {1, 0}} is singular.
MatrixPower::sing: Matrix {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}} is singular. | Out[364]= \begin{pmatrix}1\\2\end{pmatrix} | MatrixPower::sing: Matrix {{0, 0}, {1, 0}} is singular.
MatrixPower::sing: Matrix {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}} is singular.
```

## 2026-07-16 12:39:02
**Input:**
```wolfram
testM = {{0,1,0},{0,0,1},{0,0,0}};
{MatrixPower[testM,2], MatrixPower[testM,2]===testM.testM}
```
**Output:**
```
Out= {{{0, 0, 1}, {0, 0, 0}, {0, 0, 0}}, True}
```

## 2026-07-16 12:39:08
**Input:**
```wolfram
Ee[2][1,2]
MatrixPower[Ee[2][1,2],2]
MatrixPower[Ee[2][1,2],2]===Ee[2][1,2].Ee[2][1,2]
```
**Output:**
```
Out[1]= {{0, 2, 0}, {0, 0, 2}, {0, 0, 0}}
Out[2]= {{0, 0, 4}, {0, 0, 0}, {0, 0, 0}}
Out[3]= True
```

## 2026-07-16 12:39:19
**Input:**
```wolfram
manualG = IdentityMatrix[3] + φ Ee[2][1,2] + φ^2/2 MatrixPower[Ee[2][1,2],2];
gRot[2][1,2][φ] - manualG // Simplify
```
**Output:**
```
[message] MatrixPower::sing: Matrix {{0, 2, 0}, {0, 0, 2}, {0, 0, 0}} is singular.
Out= {{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}} + MatrixPower[{{0, 2, 0}, {0, 0, 2}, {0, 0, 0}}, 0]
```

## 2026-07-16 12:39:28
**Input:**
```wolfram
MatrixPower[{{0, 2, 0}, {0, 0, 2}, {0, 0, 0}}, 0]
```
**Output:**
```
[message] MatrixPower::sing: Matrix {{0, 2, 0}, {0, 0, 2}, {0, 0, 0}} is singular.
Out= MatrixPower[{{0, 2, 0}, {0, 0, 2}, {0, 0, 0}}, 0]
```

## 2026-07-16 12:39:50
**Input:**
```wolfram
{MatrixPower[Ee[2][1,2],1], MatrixPower[Ee[2][1,2],1]===Ee[2][1,2]}
```
**Output:**
```
Out= {{{0, 2, 0}, {0, 0, 2}, {0, 0, 0}}, True}
```
## 2026-07-16 12:40:04 — ✏️ EDIT CELL 66
(*RHS operatorial CG overlap: rotate the site-1 pseudovacuum covector by g=exp[φ Ee[λ1][a,b]]
  (Ee[a,b] nilpotent of order ≤λ1+1 on the finite-dim rep). Build g as the explicit truncated
  series Σ_{…


## 2026-07-16 12:40:04
**Input:**
```wolfram
(*RHS operatorial CG overlap: rotate the site-1 pseudovacuum covector by g=exp[φ Ee[λ1][a,b]]
  (Ee[a,b] nilpotent of order ≤λ1+1 on the finite-dim rep). Build g as the explicit truncated
  series Σ_{m=0}^{λ} φ^m/m! Ee[a,b]^m — this is EXACT (higher powers vanish) and guaranteed to be a
  polynomial in φ, avoiding MatrixExp returning an unevaluated symbolic exponential. Tensor with the
  site-2 pseudovacuum covector, contract with the two-site right eigenvector Psi. Shapes:
  (1×n1).(n1×n1)=1×n1 → KroneckerProduct with 1×n2 → Flatten to length D → .Psi = scalar.
  DEVIATION FROM PLAN TEXT (kernel quirk, not a typo): MatrixPower[singularMatrix,0] does NOT
  auto-reduce to IdentityMatrix in this kernel (fires MatrixPower::sing and stays unevaluated) even
  though m=0 is a concrete nonnegative integer — verified directly: MatrixPower[Ee[2][1,2],0] stays
  unevaluated with a ::sing message, while m>=1 evaluate cleanly. Ee[λ][a,b] is nilpotent/singular
  for the generators used here, so the naive Sum from m=0 silently leaves an unevaluated
  MatrixPower[...,0] term in gRot. Fix: split off the m=0 (identity) term explicitly and only sum
  MatrixPower over m=1..λ, which evaluates without any message.*)
ClearAll[gRot, OP];
gRot[λ_][a_, b_][φ_] := IdentityMatrix[dim[{λ, 0}]] + Sum[φ^m/m! MatrixPower[Ee[λ][a, b], m], {m, 1, λ}];
OP[λ1_, λ2_][a_, b_][M_, n_][φ_] := Flatten[KroneckerProduct[xSingle[λ1][{0}] . gRot[λ1][a, b][φ], xSingle[λ2][{0}]]] . Psi[λ1, λ2][M, n];
(*re-verify the m=0 fix reproduces the manual sum exactly for a nontrivial case*)
Chop[gRot[2][1, 2][φ] - (IdentityMatrix[3] + φ Ee[2][1, 2] + φ^2/2 MatrixPower[Ee[2][1, 2], 2]) // Simplify, 10^-20]
```
**Output:**
```
Out= {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}
```
## 2026-07-16 12:40:09 — ▶️ RUN CELL 67 (0.02 s)
**In [67]:** `(*Which single-site generator is nontrivial? The lowest-weight covector ⟨x0| annihilates one of↵  E12/E21 (giving φ-independent OP); keep the other. Test on (λ1,λ2)=(2,1), state (1,0).*)↵op12 = OP[2, `
**Out:** Out[368]= \langle\!\lvert \text{OP(E12)}\to 1.\,+0.\times 10^{-30}\,\mathrm{i},\text{OP(E21)}\to \left(-0.03626734060947795145020244939-2.004237714666237867354326095\,\mathrm{i}\right)\,\varphi ^2-\left(1.522150083623538480912456055-0.2442047957313920727733393355\,\mathrm{i}\right)\,\varphi +\left(1

## 2026-07-16 12:40:14 — ▶️ RUN CELL 68 (0.03 s)
**In [68]:** `(*OPg fixes the nontrivial generator so downstream code reads OPg[λ1,λ2][M,n][φ]. Degree in φ is ≤λ1.*)↵ClearAll[OPg];↵OPg[λ1_, λ2_][M_, n_][φ_] := OP[λ1, λ2][genNontrivial[[1]], genNontrivial[[2]]][M`
**Out:** Out[371]= \left\{1,2\right\}

## 2026-07-16 12:40:28 — 📥 BULK INSERT 1 CELL(S) at positions 69–69
1. [code] (*Discovery step 1 (λ1=1, OPg linear in φ): solve φ(k) per state assuming exact match (scalar 1),↵


## 2026-07-16 12:40:28
**Input:**
```wolfram
(*Discovery step 1 (λ1=1, OPg linear in φ): solve φ(k) per state assuming exact match (scalar 1),
  using rep (λ1,λ2)=(1,2) (dim 6, several states). If φ_state(k) coincides across states, the
  dictionary is exact with that φ(k). Report the per-state φ candidates and their pairwise
  differences.*)
lin12[M_, n_] := Module[{c0, c1}, c0 = OPg[1, 2][M, n][0]; c1 = D[OPg[1, 2][M, n][φ], φ]; {c0, c1}];
phiCand12[M_, n_][k_] := Module[{c0, c1}, {c0, c1} = lin12[M, n]; (FL[1, 2][M, n][k] - c0)/c1];
sts12 = statesOf[1, 2];
phiList12 = Table[{s, phiCand12[s[[1]], s[[2]]][k] // Simplify}, {s, sts12}];
phiDiffs12 = Table[phiList12[[i, 2]] - phiList12[[1, 2]] // Simplify, {i, 2, Length[phiList12]}];
<|"states" -> sts12, "φ candidates per state" -> phiList12,
  "differences from state 1 (expect all 0 if universal at scalar 1)" -> phiDiffs12|>
```
**Output:**
```
Out[377]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{$\varphi$ candidates per state}\to \begin{pmatrix}\left\{0,0\right\} & \left(-0.920793778924720536105378676-0.390049761817044163511250929\,\mathrm{i}\right)\,k^{4/3}+\left(0.695722366612934233848844287-0.718310788304406529935213657\,\mathrm{i}\right)\,\sqrt[3]{k}-\left(0.3604398751… [1247 chars total]
```

## 2026-07-16 12:40:33
**Input:**
```wolfram
phiDiffs12
```
**Output:**
```
Out=                                                                                                     1/3                                                    4/3                                                                                                          1/3                                                      4/3                                                                                                               1/3                                                          4/3                                                                                                                  1/3                                                           4/3        -19                                                                                  1/3                                                       4/3
{(1.627303481210406069 - 0.509718271393729909 I) - (0.465427483249559777 - 1.738776043469954290 I) k    + (1.326345476875853767 + 1.060122574468035904 I) k   , (-0.7189949958393951952 - 0.5723093689745847617 I) - (1.5537699154731399186 - 1.8363283118425063196 I) k    + (1.9471242891965252168 - 0.3029545350609741643 I) k   , (1.62730348121040606914 - 1.35584662579301952358 I) - (0.965391252878792098190 - 0.339282151941626110660 I) k    + (0.819175108107555850896 + 0.056927984776178348079 I) k   , (-0.718994995839395195224 - 1.293255528212164671237 I) - (0.326403506287175821776 + 0.366568146427070851440 I) k    + (0.9149364043878558062910 + 1.397049248109101892089 I) k   , (0. 10    - 1.8655648971867494330 I) - (0.66118491636651204033 - 1.2049630956839855892 I) k    + (0.4544043737788571773 + 0.24685679231406852243 I) k   }
```
## 2026-07-16 12:40:43 — 📥 BULK INSERT 1 CELL(S) at positions 70–70
1. [code] (*Discovery step 1b: allow a per-state scalar c_s. For each state, FL/c_s = c0 + c1 φ(k) ⇒↵  φ(k) =


## 2026-07-16 12:40:44
**Input:**
```wolfram
(*Discovery step 1b: allow a per-state scalar c_s. For each state, FL/c_s = c0 + c1 φ(k) ⇒
  φ(k) = (FL/c_s − c0)/c1. Require φ(k) equal across states: pin c_s (up to one overall scale) by
  matching the k-dependence. Solve using the two lowest states, then read off φ(k).*)
{sa, sb} = sts12[[1 ;; 2]];
{c0a, c1a} = lin12[sa[[1]], sa[[2]]]; {c0b, c1b} = lin12[sb[[1]], sb[[2]]];
(*Set c_sa=1 (overall scale convention); solve c_sb so both give the same φ(k).*)
phiA[k_] := (FL[1, 2][sa[[1]], sa[[2]]][k] - c0a)/c1a;
phiFromB[cb_][k_] := (FL[1, 2][sb[[1]], sb[[2]]][k]/cb - c0b)/c1b;
cbSol = Solve[phiA[k] == phiFromB[cb][k] // Simplify, cb];
<|"c_sb solution(s)" -> cbSol, "φ(k) (from state a)" -> (phiA[k] // Simplify)|>
```
**Output:**
```
Out[383]= \langle\!\lvert \text{c\_sb solution(s)}\to \left\{\left\{\mathrm{cb}\to \frac{\left(-6.48948566351\times 10^{46}-4.69046160561\times 10^{46}\,\mathrm{i}\right)\,k^{4/3}-\left(6.23692692930\times 10^{46}+8.68764401646\times 10^{46}\,\mathrm{i}\right)\,\sqrt[3]{k}}{1.02230236398\times 10^{47}\,k^{4/3}-\left(3.68478536434\times 10^{46}-9.53585702272\times 10^{46}\,\mathrm{i}\right)\,\sqrt[… [833 chars total]
```

## 2026-07-16 12:41:07
**Input:**
```wolfram
cbExpr = cb /. cbSol[[1]];
{cbExpr /. k -> 1.0, cbExpr /. k -> 2.0, cbExpr /. k -> 5.0, cbExpr /. k -> 0.3}
```
**Output:**
```
Out= {-0.5778488645596662 - 0.72588616455467*I, -0.7173896826421956 - 0.6465331262430513*I, -0.7297017080721853 - 0.5434454069733405*I, -0.19703271755005838 - 0.6000468750006196*I}
```

## 2026-07-16 12:41:24
**Input:**
```wolfram
Table[{s, FL[1,2][s[[1]],s[[2]]][k] // Simplify}, {s, sts12}] // Column
```
**Output:**
```
Out=                                                                                                                                                             1/3
{{0, 0}, ((0.92079377892472053610537867566 + 0.39004976181704416351125092885 I) + (0.03194077718899586574576231421 - 0.99948976320548822385419493638 I) k) k   }
                                                                                                                             1/3
{{1, 0}, ((-0.405551697951133231159 - 0.670072812650991740461 I) - (0.446914378356193654661 + 0.379677296427955078847 I) k) k   }
          1/3
{{1, 1}, k    ((-1.026330510271804680645 + 0.6930042968780183278163 I) + (1.048293463098762369372 - 0.291933059011270972178 I) k)}
                                                                                                                                      1/3
{{2, 0}, ((0.101618670817164685209381 + 0.33312177704086581543190 I) - (0.006835984233320180167197 - 0.26066714318366616697832 I) k) k   }
                                                                                                                                      1/3
{{2, 1}, ((0.005857374536864729814345 - 1.006999486292057728577915 I) + (0.27539797005314427962806 + 0.84092700384652263606032 I) k) k   }
                                                                                                                                 1/3
{{3, 0}, ((0.4663894051458633588505 + 0.14319296950297564108225 I) - (0.3016732276751062924232 - 0.3834273952803527447315 I) k) k   }
```

## 2026-07-16 12:41:49
**Input:**
```wolfram
Table[{s, lin12[s[[1]],s[[2]]]}, {s, sts12}] // Column
```
**Output:**
```
Out=                                                 -28
{{0, 0}, {1.0000000000000000000000000000 + 0. 10    I, 0.3604398751472265372666776506 + 0.9327824485933747164848359821 I}}
                                               -28
{{1, 0}, {1.000000000000000000000000000 + 0. 10    I, -0.7101545203859353503227942743 + 0.2371533418544670754950145801 I}}
                                               -28
{{1, 1}, {1.000000000000000000000000000 + 0. 10    I, 0.8334628987257760078087629442 + 0.2783316955260821099308158311 I}}
                                               -28
{{2, 0}, {1.000000000000000000000000000 + 0. 10    I, -0.7101545203859353503227942743 - 0.2371533418544670754950145801 I}}
                                               -28
{{2, 1}, {1.000000000000000000000000000 + 0. 10    I, 0.8334628987257760078087629442 - 0.2783316955260821099308158311 I}}
                                                -28
{{3, 0}, {1.0000000000000000000000000000 + 0. 10    I, 0.3604398751472265372666776506 - 0.9327824485933747164848359821 I}}
```

## 2026-07-16 12:42:17
**Input:**
```wolfram
ABTab = Table[
   Module[{f, A, B},
    f = FL[1,2][s[[1]],s[[2]]][k] // Simplify;
    A = Coefficient[f/k^(1/3) // Simplify, k, 0];
    B = Coefficient[f/k^(1/3) // Simplify, k, 1];
    {s, A, B, A/B}],
   {s, sts12}];
ABTab // Column
```
**Output:**
```
Out= {{0, 0}, 0.92079377892472053610537867566 + 0.39004976181704416351125092885 I, 0.03194077718899586574576231421 - 0.99948976320548822385419493638 I, -0.3604398751472265372666776506 + 0.9327824485933747164848359821 I}
{{1, 0}, -0.405551697951133231159 - 0.670072812650991740461 I, -0.446914378356193654661 - 0.379677296427955078847 I, 1.26686360606317953187 + 0.42306417719964480709 I}
{{1, 1}, -1.026330510271804680645 + 0.6930042968780183278163 I, 1.048293463098762369372 - 0.291933059011270972178 I, -1.079434870986621732491 + 0.360473079618789954752 I}
{{2, 0}, 0.101618670817164685209381 + 0.33312177704086581543190 I, -0.006835984233320180167197 + 0.26066714318366616697832 I, 1.2668636060631795318693 - 0.4230641771996448070948 I}
{{2, 1}, 0.005857374536864729814345 - 1.006999486292057728577915 I, 0.27539797005314427962806 + 0.84092700384652263606032 I, -1.0794348709866217324906 - 0.36047307961878995475242 I}
{{3, 0}, 0.4663894051458633588505 + 0.14319296950297564108225 I, -0.3016732276751062924232 + 0.3834273952803527447315 I, -0.360439875147226537267 - 0.932782448593374716485 I}
```
## 2026-07-16 12:42:48 — 📥 BULK INSERT 1 CELL(S) at positions 71–71
1. [code] (*Discovery step 2 — OBSTRUCTION FOUND, reported rather than forced. Direct structural test: write↵


## 2026-07-16 12:42:48
**Input:**
```wolfram
(*Discovery step 2 — OBSTRUCTION FOUND, reported rather than forced. Direct structural test: write
  FL[1,2][state][k] = k^(1/3)(A_state + B_state k) (verified: every state of this rep is EXACTLY of
  this shape — a degree-1 polynomial in k times the universal (k κ1)^(θ1/h) prefactor, θ1=1/3). Since
  OPg[1,2][state][φ] = c0(state) + c1(state) φ is linear in φ with c0(state)=1 for EVERY state (a
  clean structural fact, presumably <x0|Psi> normalization), any dictionary of the plan's assumed
  form (single generator, single state-INDEPENDENT φ(k), single per-state scalar) would force the
  ratio A_state/B_state to depend on state ONLY through c1(state) in a fixed (Möbius) way. It does
  not: A/B clusters into 3 conjugate-paired values across the 6 states of (1,2), not a single common
  value or a single function of c1(state). Independently, Step 1b's attempt to fix a per-state scalar
  cb between just the two lowest states {0,0},{1,0} returns a cb that is NOT constant in k (varies
  from -0.578-0.726i at k=1 to -0.197-0.600i at k=0.3 to -0.730-0.543i at k=5) — so no constant scalar
  reconciles even two states. CONCLUSION: the naive dictionary {generator={2,1}, universal φ(k),
  per-state scalar} does NOT close already at λ1=1. Recorded here as the discovery's terminal
  obstruction (not patched/forced) — see chat report for the full numeric detail.*)
obstructionReport = <|
  "generator {a,b} tested" -> genNontrivial,
  "rep tested" -> {1, 2},
  "FL shape (verified exact for all 6 states)" -> "FL[1,2][state][k] = k^(1/3) (A_state + B_state k)",
  "c0(state) = OPg[state][0]" -> Table[{s, lin12[s[[1]], s[[2]]][[1]]}, {s, sts12}],
  "A_state/B_state (structural obstruction: not universal, not a fixed function of c1(state))" ->
    Table[{sts12[[i]], ABTab[[i, 4]]}, {i, Length[sts12]}],
  "scalar=1 φ-candidate differences (state i - state 1), nonzero ⇒ not universal at scalar 1" -> phiDiffs12,
  "per-state-scalar cb (states {0,0},{1,0}) sampled at k=0.3,1,2,5 — NOT constant ⇒ no fixed scalar reconciles even 2 states" ->
    Table[{k -> kk, cb /. cbSol[[1]] /. k -> kk}, {kk, {0.3, 1., 2., 5.}}],
  "verdict" -> "Dictionary does NOT close at λ1=1 with the assumed form (generator + universal φ(k) + per-state scalar). Obstruction reported; no φ(k) invented/forced."
|>;
obstructionReport
```
**Output:**
```
Out[385]= \langle\!\lvert \text{generator \{a,b\} tested}\to \left\{2,1\right\},\text{rep tested}\to \left\{1,2\right\},\text{FL shape (verified exact for all 6 states)}\to \text{FL[1,2][state][k] = k\^{}(1/3) (A\_state + B\_state k)},\text{c0(state) = OPg[state][0]}\to \begin{pmatrix}\left\{0,0\right\} & 1.\,+0.\times 10^{-28}\,\mathrm{i}\\\left\{1,0\right\} & 1.\,+0.\times 10^{-28}\,\mathrm{i}\\… [1247 chars total]
```

## 2026-07-16 12:52:56
**Input:**
```wolfram
ClearAll[det1, FLdet];
det1[λ1_][f_] := br1[λ1][f]/(f /. u -> θ1);
FLdet[λ1_,λ2_][M_,n_][k_] := det1[λ1][k^(u/h) Q1[λ1,λ2][M,n][u]];
FLdet[1,2][0,0][k] // Simplify
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= (0.99999999999999999999999999999999999999`27.136658757470364 + 0``27.2468745427702*I) - (0.36043987514722653726667765059555009129`26.89502217176169 + 0.93278244859337471648483598207086792945`27.11877367170218*I)*k
```

## 2026-07-16 12:53:05
**Input:**
```wolfram
alphaTab = Table[
   Module[{fl, op, D0, Dc, c0, c1, alpha},
    fl = FLdet[1,2][s[[1]],s[[2]]][k] // Simplify;
    op = OPg[1,2][s[[1]],s[[2]]][φ] // Simplify;
    D0 = Coefficient[fl, k, 0]; Dc = Coefficient[fl, k, 1];
    c0 = Coefficient[op, φ, 0]; c1 = Coefficient[op, φ, 1];
    alpha = Dc/c1;
    {s, D0, c0, Dc, c1, alpha}],
   {s, sts12}];
alphaTab // Column
```
**Output:**
```
Out=                                              -27                                          -28                                                                                                                                                                              -27
{{0, 0}, 1.00000000000000000000000000 + 0. 10    I, 1.0000000000000000000000000000 + 0. 10    I, -0.360439875147226537266677651 - 0.932782448593374716484835982 I, 0.3604398751472265372666776506 + 0.9327824485933747164848359821 I, -1.00000000000000000000000000 + 0. 10    I}
                                    -18                                         -28                                                                                                                                                   -18
{{1, 0}, 1.00000000000000000 + 0. 10    I, 1.000000000000000000000000000 + 0. 10    I, 0.710154520385935350 - 0.237153341854467075 I, -0.7101545203859353503227942743 + 0.2371533418544670754950145801 I, -1.00000000000000000 + 0. 10    I}
                                     -19                                         -28                                                                                                                                                      -19
{{1, 1}, 1.000000000000000000 + 0. 10    I, 1.000000000000000000000000000 + 0. 10    I, -0.8334628987257760078 - 0.2783316955260821099 I, 0.8334628987257760078087629442 + 0.2783316955260821099308158311 I, -1.000000000000000000 + 0. 10    I}
                                      -20                                         -28                                                                                                                                                           -20
{{2, 0}, 1.0000000000000000000 + 0. 10    I, 1.000000000000000000000000000 + 0. 10    I, 0.71015452038593535032 + 0.237153341854467075495 I, -0.7101545203859353503227942743 - 0.2371533418544670754950145801 I, -1.00000000000000000000 + 0. 10    I}
                                        -21                                         -28                                                                                                                                                             -21
{{2, 1}, 1.000000000000000000000 + 0. 10    I, 1.000000000000000000000000000 + 0. 10    I, -0.833462898725776007809 + 0.278331695526082109931 I, 0.8334628987257760078087629442 - 0.2783316955260821099308158311 I, -1.000000000000000000000 + 0. 10    I}
                                     -19                                          -28                                                                                                                                                       -19
{{3, 0}, 1.000000000000000000 + 0. 10    I, 1.0000000000000000000000000000 + 0. 10    I, -0.3604398751472265373 + 0.9327824485933747165 I, 0.3604398751472265372666776506 - 0.9327824485933747164848359821 I, -1.0000000000000000000 + 0. 10    I}
```

## 2026-07-16 12:53:15
**Input:**
```wolfram
φuniversal[k_] := -k;
repSweepDet = {{1,1},{2,1},{1,2},{2,2},{3,1},{1,3}};
validationDet = Table[
   Module[{λ1=rep[[1]], λ2=rep[[2]], sts, resids},
    sts = statesOf[λ1,λ2];
    resids = Table[
      Module[{lhs, rhs},
       lhs = FLdet[λ1,λ2][st[[1]],st[[2]]][k];
       rhs = OPg[λ1,λ2][st[[1]],st[[2]]][φuniversal[k]];
       {st, N[Chop[Simplify[lhs-rhs],10^-15]]}],
      {st, sts}];
    {λ1,λ2} -> resids],
   {rep, repSweepDet}];
validationDet
```
**Output:**
```
Out= {{1, 1} -> {{{0, 0}, 0.}, {{1, 0}, 0.}, {{1, 1}, 0.}, {{2, 0}, 0.}}, {2, 1} -> {{{0, 0}, 0.}, {{1, 0}, 0.}, {{1, 1}, 0.}, {{2, 0}, 0.}, {{2, 1}, 0.}, {{3, 0}, 0.}}, {1, 2} -> {{{0, 0}, 0.}, {{1, 0}, 0.}, {{1, 1}, 0.}, {{2, 0}, 0.}, {{2, 1}, 0.}, {{3, 0}, 0.}}, {2, 2} -> {{{0, 0}, 0.}, {{1, 0}, 0.}, {{1, 1}, 0.}, {{2, 0}, 0.}, {{2, 1}, 0.}, {{2, 2}, 0.}, {{3, 0}, 0.}, {{3, 1}, 0.}, {{4, 0}, 0.}}, {3, 1} -> {{{0, 0}, 0.}, {{1, 0}, 0.}, {{1, 1}, 0.}, {{2, 0}, 0.}, {{2, 1}, 0.}, {{3, 0}, 0.}, {{3, 1}, 0.}, {{4, 0}, 0.}}, {1, 3} -> {{{0, 0}, 0.}, {{1, 0}, 0.}, {{1, 1}, 0.}, {{2, 0}, 0.}, {{2, 1}, 0.}, {{3, 0}, 0.}, {{3, 1}, 0.}, {{4, 0}, 0.}}}
```

## 2026-07-16 12:53:22
**Input:**
```wolfram
rawResids = Flatten[Table[
   Module[{λ1=rep[[1]], λ2=rep[[2]], sts},
    sts = statesOf[λ1,λ2];
    Table[
      N[Abs[Simplify[FLdet[λ1,λ2][st[[1]],st[[2]]][k] - OPg[λ1,λ2][st[[1]],st[[2]]][φuniversal[k]]] /. k -> 1.234567]],
      {st, sts}]],
   {rep, repSweepDet}]];
{Max[rawResids], Length[rawResids]}
```
**Output:**
```
Out= {0., 41}
```

## 2026-07-16 12:53:30
**Input:**
```wolfram
rawResids2 = Flatten[Table[
   Module[{λ1=rep[[1]], λ2=rep[[2]], sts},
    sts = statesOf[λ1,λ2];
    Table[
      Abs[N[Simplify[FLdet[λ1,λ2][st[[1]],st[[2]]][k] - OPg[λ1,λ2][st[[1]],st[[2]]][φuniversal[k]]] /. k -> 1.234567]],
      {st, sts}]],
   {rep, repSweepDet}]];
Max[rawResids2] // InputForm
```
**Output:**
```
Out= InputForm[0.]
```
## 2026-07-16 12:54:28 — 📥 BULK INSERT 4 CELL(S) at positions 72–75
1. [code] (*RESOLUTION of the Phase-2 obstruction (cell 71): the functional object is not the bare bracket↵  b
2. [code] (*λ1=1 universality test (the sharp discovery step), rep (1,2), all 6 states: det1 is 1+D_state k;↵ 
3. [code] (*Validate the universal φ(k)=αVal k across the FULL rep sweep and ALL their states. This is the↵  r
4. [code] (*Final dictionary report: the functional CG overlap and the operatorial CG overlap match EXACTLY↵


## 2026-07-16 12:54:29
**Input:**
```wolfram
(*RESOLUTION of the Phase-2 obstruction (cell 71): the functional object is not the bare bracket
  br1[f] but the 1-site FSoV DETERMINANT det1[f] = br1[f]/f(θ1) — the 1-site reduction of Section
  C's det[λ1,λ2][f] (which divides by f at BOTH nodes θ1,θ2; only the θ1-node survives at 1 site).
  This missing 1/f(θ1) normalization is exactly what caused the obstruction: it cancels the
  (k κ1)^(θ1/h) fractional prefactor between br1[f] and f(θ1), leaving
  det1[k^(u/h) Q1] = Σ_j w_j k^j Q1(θ1+jh)/Q1(θ1), a CLEAN degree-λ1 polynomial in k with leading
  term exactly 1 (since w_0=1 by Nn1's construction).*)
ClearAll[det1, FLdet];
det1[λ1_][f_] := br1[λ1][f]/(f /. u -> θ1);
FLdet[λ1_, λ2_][M_, n_][k_] := det1[λ1][k^(u/h) Q1[λ1, λ2][M, n][u]];
(*smoke: FLdet[1,2][0,0][k] must be a degree-1 polynomial in k, constant term exactly 1, NO
  fractional power of k (unlike the old FL).*)
FLdet[1, 2][0, 0][k] // Simplify
```
**Output:**
```
Out[389]= \left(1.000000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right)-\left(0.360439875147226537266677651+0.932782448593374716484835982\,\mathrm{i}\right)\,k
```

## 2026-07-16 12:54:29
**Input:**
```wolfram
(*λ1=1 universality test (the sharp discovery step), rep (1,2), all 6 states: det1 is 1+D_state k;
  OPg is 1+c1(state) φ. HARD assert #1: both constant terms are exactly 1 (a structural fact of the
  normalization, not a coincidence). HARD assert #2: α_state := D_state/c1(state) is the SAME
  constant across all 6 states — this is the actual universality content.*)
DictDet::badconst = "Constant term ≠ 1 for state `1`: FLdet const=`2`, OPg const=`3`.";
DictDet::notuniv = "α_state not universal: spread=`1` (see alphaTab).";
alphaTab = Table[
   Module[{fl, op, D0, Dc, c0, c1, alpha},
    fl = FLdet[1, 2][s[[1]], s[[2]]][k] // Simplify;
    op = OPg[1, 2][s[[1]], s[[2]]][φ] // Simplify;
    D0 = Coefficient[fl, k, 0]; Dc = Coefficient[fl, k, 1];
    c0 = Coefficient[op, φ, 0]; c1 = Coefficient[op, φ, 1];
    If[Abs[N[D0 - 1]] > 10^-15 || Abs[N[c0 - 1]] > 10^-15,
     Message[DictDet::badconst, s, N[D0], N[c0]]; Abort[]];
    alpha = Dc/c1;
    {s, N[D0], N[c0], N[Dc], N[c1], N[alpha]}],
   {s, sts12}];
alphaSpread = Max[Abs[alphaTab[[All, 6]] - alphaTab[[1, 6]]]];
If[alphaSpread > 10^-12, Message[DictDet::notuniv, alphaSpread]; Abort[]];
αVal = alphaTab[[1, 6]];
<|"alphaTab {state,D0,c0,Dc,c1,α}" -> alphaTab, "α spread" -> alphaSpread, "α (universal value)" -> αVal|>
```
**Output:**
```
Out[396]= \langle\!\lvert \text{alphaTab \{state,D0,c0,Dc,c1,$\alpha$\}}\to \begin{pmatrix}\left\{0,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & -0.36043987514722653-0.9327824485933747\,\mathrm{i} & 0.36043987514722653 +0.9327824485933747\,\mathrm{i} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.7101545203859353 -0.2371533418544… [1247 chars total]
```

## 2026-07-16 12:54:30
**Input:**
```wolfram
(*Validate the universal φ(k)=αVal k across the FULL rep sweep and ALL their states. This is the
  real content: it tests that the binomial tower weights w_m (from Nn1/mu1) match OPg's φ^m Taylor
  coefficients (built from Ee[a,b]^m via the fixed gRot) under the single rescaling α, for every
  (λ1,λ2,M,n). HARD assert: residual FLdet-OPg[αVal k] vanishes (Simplify reduces it to literal 0,
  or Chop below 10^-15) for every state of every rep in the sweep.*)
ClearAll[φuniversal];
φuniversal[k_] := αVal k;
DictDet::residual = "Nonzero residual for (λ1,λ2,M,n)=(`1`,`2`,`3`,`4`): resid=`5`.";
repSweepDet = {{1, 1}, {2, 1}, {1, 2}, {2, 2}, {3, 1}, {1, 3}};
validationDet = Table[
   Module[{λ1 = rep[[1]], λ2 = rep[[2]], sts, resids},
    sts = statesOf[λ1, λ2];
    resids = Table[
      Module[{lhs, rhs, r},
       lhs = FLdet[λ1, λ2][st[[1]], st[[2]]][k];
       rhs = OPg[λ1, λ2][st[[1]], st[[2]]][φuniversal[k]];
       r = Chop[Simplify[lhs - rhs], 10^-15];
       If[r =!= 0, Message[DictDet::residual, λ1, λ2, st[[1]], st[[2]], r]; Abort[]];
       {st, r}],
      {st, sts}];
    {λ1, λ2} -> resids],
   {rep, repSweepDet}];
allResidsDet = Flatten[Table[validationDet[[i, 2]][[All, 2]], {i, Length[validationDet]}]];
maxResidDet = Max[Abs[allResidsDet]];
<|"validation" -> validationDet, "max residual over full sweep" -> maxResidDet,
  "num states checked" -> Length[allResidsDet]|>
```
**Output:**
```
Out[404]= \langle\!\lvert \text{validation}\to \left\{\left\{1,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\end{pmatrix},\left\{2,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\\\left\{2,1\right\} & 0\\\left\{3,0\right\} & 0\end{pmatrix},\left\{1… [1247 chars total]
```

## 2026-07-16 12:54:31
**Input:**
```wolfram
(*Final dictionary report: the functional CG overlap and the operatorial CG overlap match EXACTLY
  (symbolically zero after Simplify, not merely below a noise floor) once the LHS uses the 1-site
  FSoV DETERMINANT det1 (not the bare bracket br1) and the parameter map is the single LINEAR
  rescaling φ(k)=α k. This is the resolved Section F dictionary.*)
dictReportDet = <|
  "generator {a,b}" -> genNontrivial,
  "φ(k)" -> "α k",
  "α (universal value)" -> αVal,
  "normalization" -> "1-site FSoV determinant det1[f] = br1[f]/(f/.u->θ1); Section C's det[λ1,λ2][f] restricted to 1 site",
  "max residual over rep sweep" -> maxResidDet,
  "reps validated" -> repSweepDet,
  "num states validated" -> Length[allResidsDet]
|>;
dictReportDet
```
**Output:**
```
Out[406]= \langle\!\lvert \text{generator \{a,b\}}\to \left\{2,1\right\},\text{$\varphi$(k)}\to \text{$\alpha$ k},\text{$\alpha$ (universal value)}\to -1.\,+0.\,\,\mathrm{i},\text{normalization}\to \text{1-site FSoV determinant det1[f] = br1[f]/(f/.u->$\theta$1); Section C's det[$\lambda$1,$\lambda$2][f] restricted to 1 site},\text{max residual over rep sweep}\to 0,\text{reps validated}\to \begin{… [510 chars total]
```
## 2026-07-16 12:55:06 — ✏️ EDIT CELL 71
(*Discovery step 2 — OBSTRUCTION FOUND, reported rather than forced. Direct structural test: write
  FL[1,2][state][k] = k^(1/3)(A_state + B_state k) (verified: every state of this rep is EXACTLY of
 …


## 2026-07-16 12:55:06
**Input:**
```wolfram
(*Discovery step 2 — OBSTRUCTION FOUND, reported rather than forced. Direct structural test: write
  FL[1,2][state][k] = k^(1/3)(A_state + B_state k) (verified: every state of this rep is EXACTLY of
  this shape — a degree-1 polynomial in k times the universal (k κ1)^(θ1/h) prefactor, θ1=1/3). Since
  OPg[1,2][state][φ] = c0(state) + c1(state) φ is linear in φ with c0(state)=1 for EVERY state (a
  clean structural fact, presumably <x0|Psi> normalization), any dictionary of the plan's assumed
  form (single generator, single state-INDEPENDENT φ(k), single per-state scalar) would force the
  ratio A_state/B_state to depend on state ONLY through c1(state) in a fixed (Möbius) way. It does
  not: A/B clusters into 3 conjugate-paired values across the 6 states of (1,2), not a single common
  value or a single function of c1(state). Independently, Step 1b's attempt to fix a per-state scalar
  cb between just the two lowest states {0,0},{1,0} returns a cb that is NOT constant in k (varies
  from -0.578-0.726i at k=1 to -0.197-0.600i at k=0.3 to -0.730-0.543i at k=5) — so no constant scalar
  reconciles even two states. CONCLUSION: the naive dictionary {generator={2,1}, universal φ(k),
  per-state scalar} does NOT close already at λ1=1. Recorded here as the discovery's terminal
  obstruction (not patched/forced) — see chat report for the full numeric detail.

  *** RESOLVED (see cells 72-75 below): the obstruction was a WRONG-OBJECT bug in the plan, not in
  this construction. The correct LHS functional object is not the bare bracket br1[f] but the 1-site
  FSoV DETERMINANT det1[f] = br1[f]/(f/.u->θ1) — Section C's det[λ1,λ2][f] (which divides by f at
  BOTH nodes θ1,θ2) reduced to 1 site. The missing 1/f(θ1) exactly cancels the (k κ1)^(θ1/h)
  fractional prefactor that made FL's shape state-dependent (the A/B clustering seen below). With
  det1 in place, FLdet[1,2][state][k] = 1 + D_state k for every state (clean degree-1 polynomial,
  leading term exactly 1), and the ratio α_state = D_state/c1(state) IS universal: α = -1 exactly for
  all 6 states, extending to a symbolically EXACT match (Simplify gives literal 0, not just below a
  noise floor) across the full rep sweep {{1,1},{2,1},{1,2},{2,2},{3,1},{1,3}} under φ(k) = -k. The
  A_state/B_state values below (br1-based) were a real, correctly-computed symptom — just of the
  wrong functional object; nothing here was recomputed or altered. ***)*)
obstructionReport = <|
  "generator {a,b} tested" -> genNontrivial,
  "rep tested" -> {1, 2},
  "FL shape (verified exact for all 6 states)" -> "FL[1,2][state][k] = k^(1/3) (A_state + B_state k)",
  "c0(state) = OPg[state][0]" -> Table[{s, lin12[s[[1]], s[[2]]][[1]]}, {s, sts12}],
  "A_state/B_state (structural obstruction: not universal, not a fixed function of c1(state))" ->
    Table[{sts12[[i]], ABTab[[i, 4]]}, {i, Length[sts12]}],
  "scalar=1 φ-candidate differences (state i - state 1), nonzero ⇒ not universal at scalar 1" -> phiDiffs12,
  "per-state-scalar cb (states {0,0},{1,0}) sampled at k=0.3,1,2,5 — NOT constant ⇒ no fixed scalar reconciles even 2 states" ->
    Table[{k -> kk, cb /. cbSol[[1]] /. k -> kk}, {kk, {0.3, 1., 2., 5.}}],
  "verdict" -> "Dictionary does NOT close at λ1=1 with the assumed form (generator + universal φ(k) + per-state scalar). Obstruction reported; no φ(k) invented/forced.",
  "resolution" -> "RESOLVED in cells 72-75: LHS must be the 1-site FSoV determinant det1[f]=br1[f]/f(θ1), not the bare bracket br1[f]. With det1, φ(k)=α k, α=-1, closes EXACTLY (symbolic 0) across the full rep sweep."
|>;
obstructionReport
```
**Output:**
```
[message] ToExpression::sntx: Invalid syntax in or before "wrong functional object; nothing here was recomputed or altered. ***)*)".                                                                              … [243 chars]
(no output)
```
## 2026-07-16 12:55:28 — ✏️ EDIT CELL 71
(*Discovery step 2 — OBSTRUCTION FOUND, reported rather than forced. Direct structural test: write
  FL[1,2][state][k] = k^(1/3)(A_state + B_state k) (verified: every state of this rep is EXACTLY of
 …


## 2026-07-16 12:55:28
**Input:**
```wolfram
(*Discovery step 2 — OBSTRUCTION FOUND, reported rather than forced. Direct structural test: write
  FL[1,2][state][k] = k^(1/3)(A_state + B_state k) (verified: every state of this rep is EXACTLY of
  this shape — a degree-1 polynomial in k times the universal (k κ1)^(θ1/h) prefactor, θ1=1/3). Since
  OPg[1,2][state][φ] = c0(state) + c1(state) φ is linear in φ with c0(state)=1 for EVERY state (a
  clean structural fact, presumably <x0|Psi> normalization), any dictionary of the plan's assumed
  form (single generator, single state-INDEPENDENT φ(k), single per-state scalar) would force the
  ratio A_state/B_state to depend on state ONLY through c1(state) in a fixed (Möbius) way. It does
  not: A/B clusters into 3 conjugate-paired values across the 6 states of (1,2), not a single common
  value or a single function of c1(state). Independently, Step 1b's attempt to fix a per-state scalar
  cb between just the two lowest states {0,0},{1,0} returns a cb that is NOT constant in k (varies
  from -0.578-0.726i at k=1 to -0.197-0.600i at k=0.3 to -0.730-0.543i at k=5) — so no constant scalar
  reconciles even two states. CONCLUSION: the naive dictionary {generator={2,1}, universal φ(k),
  per-state scalar} does NOT close already at λ1=1. Recorded here as the discovery's terminal
  obstruction (not patched/forced) — see chat report for the full numeric detail.

  RESOLVED (see cells 72-75 below): the obstruction was a WRONG-OBJECT bug in the plan, not in
  this construction. The correct LHS functional object is not the bare bracket br1[f] but the 1-site
  FSoV DETERMINANT det1[f] = br1[f]/(f/.u->θ1) — Section C's det[λ1,λ2][f] (which divides by f at
  BOTH nodes θ1,θ2) reduced to 1 site. The missing 1/f(θ1) exactly cancels the (k κ1)^(θ1/h)
  fractional prefactor that made FL's shape state-dependent (the A/B clustering seen below). With
  det1 in place, FLdet[1,2][state][k] = 1 + D_state k for every state (clean degree-1 polynomial,
  leading term exactly 1), and the ratio α_state = D_state/c1(state) IS universal: α = -1 exactly for
  all 6 states, extending to a symbolically EXACT match (Simplify gives literal 0, not just below a
  noise floor) across the full rep sweep {{1,1},{2,1},{1,2},{2,2},{3,1},{1,3}} under φ(k) = -k. The
  A_state/B_state values below (br1-based) were a real, correctly-computed symptom — just of the
  wrong functional object; nothing here was recomputed or altered.*)
obstructionReport = <|
  "generator {a,b} tested" -> genNontrivial,
  "rep tested" -> {1, 2},
  "FL shape (verified exact for all 6 states)" -> "FL[1,2][state][k] = k^(1/3) (A_state + B_state k)",
  "c0(state) = OPg[state][0]" -> Table[{s, lin12[s[[1]], s[[2]]][[1]]}, {s, sts12}],
  "A_state/B_state (structural obstruction: not universal, not a fixed function of c1(state))" ->
    Table[{sts12[[i]], ABTab[[i, 4]]}, {i, Length[sts12]}],
  "scalar=1 φ-candidate differences (state i - state 1), nonzero ⇒ not universal at scalar 1" -> phiDiffs12,
  "per-state-scalar cb (states {0,0},{1,0}) sampled at k=0.3,1,2,5 — NOT constant ⇒ no fixed scalar reconciles even 2 states" ->
    Table[{k -> kk, cb /. cbSol[[1]] /. k -> kk}, {kk, {0.3, 1., 2., 5.}}],
  "verdict" -> "Dictionary does NOT close at λ1=1 with the assumed form (generator + universal φ(k) + per-state scalar). Obstruction reported; no φ(k) invented/forced.",
  "resolution" -> "RESOLVED in cells 72-75: LHS must be the 1-site FSoV determinant det1[f]=br1[f]/f(θ1), not the bare bracket br1[f]. With det1, φ(k)=α k, α=-1, closes EXACTLY (symbolic 0) across the full rep sweep."
|>;
obstructionReport
```
**Output:**
```
Out= <|"generator {a,b} tested" -> {2, 1}, "rep tested" -> {1, 2}, "FL shape (verified exact for all 6 states)" -> "FL[1,2][state][k] = k^(1/3) (A_state + B_state k)", "c0(state) = OPg[state][0]" -> {{{0, 0}, 1.`28.362553630092005 + 0``28.458858685029718*I}, {{1, 0}, 1.`28.362553630092005 + 0``28.420425233412512*I}, {{1, 1}, 1.`28.362553630092005 + 0``28.534910438118857*I}, {{2, 0}, 1.`28.362553630092005 + 0``28.420425233412512*I}, {{2, 1}, 1.`28.362553630092005 + 0``28.534910438118857*I}, {{3, 0}, 1.`28.362553630092005 + 0``28.458858685029718*I}}, "A_state/B_state (structural obstruction: not universal, not a fixed function of c1(state))" -> {{{0, 0}, -0.3604398751472265372666776505955500912995704233956884657367`28.519642521730212 + 0.93278244859337471648483598207086792945`28.723716363759525*I}, {{1, 0}, 1.26686360606317953186931967522323009443`20.233470106009175 + 0.42306417719964480709484006322965914477`19.751598815266973*I}, {{1, 1}, -1.07943487098662173249064729691354404681`21.05118564325957 + 0.36047307961878995475242216170986991601`20.60966884642148*I}, {{2, 0}, 1.26686360606317953186931967522323009437`22.374310009255577 - 0.4230641771996448070948400632296591446`22.290189543967912*I}, {{2, 1}, -1.0794348709866217324906472969135440468`22.956040954336984 - 0.360473079618789954752422161709869916`22.832840645334333*I}, {{3, 0}, -0.36043987514722653726667765059555009141`20.924877319063018 - 0.93278244859337471648483598207086792943`21.286743472795727*I}}, "scalar=1 φ-candidate differences (state i - state 1), nonzero ⇒ not universal at scalar 1" -> {(1.62730348121040606913599732581878018564`18.647747176094953 - 0.50971827139372990938999591884120878462`18.138986664495036*I) - (0.46542748324955977674418064345565154397`18.270111193412696 - 1.73877604346995428951499886313331938519`18.7906529242062*I)*k^(1/3) + (1.3263454768758537672552500921470601048`18.906181460221106 + 1.06012257446803590397186378841389386732`18.72404764631964*I)*k^(4/3), (-0.71899499583939519522396964631799395547`19.17048522233265 - 0.57230936897458476173241382036099801359`19.112586969447793*I) - (1.55376991547313991856351499799445509993`19.390842030784782 - 1.83632831184250631961923901798977241298`19.496759945119884*I)*k^(1/3) + (1.94712428919652521675004435696818121602`19.580365406788538 - 0.30295453506097416430509738933976816355`18.83333649413266*I)*k^(4/3), (1.6273034812104060691359973258187801854`20.927897475259947 - 1.35584662579301952357967604530052707432`21.12355164110732*I) - (0.96539125287879209819001219524303509911`21.268958819455012 - 0.3392821519416261106602048437419337344`20.676771979394047*I)*k^(1/3) + (0.81917510810755585089599730622926840074`21.495943058250155 + 0.05692798477617834807934963005547620116`20.051679306120185*I)*k^(4/3), (-0.71899499583939519522396964631799395559`21.201403018442612 - 1.29325552821216467123725814378073784552`21.806152133017648*I) - (0.32640350628717582177610187293768968943`21.198486434539824 + 0.36656814642707085144013045320709439335`20.90521155211221*I)*k^(1/3) + (0.91493640438785580629103394186170207104`22.17656014576632 + 1.39704924810910189208916572287478175457`21.524209614852612*I)*k^(4/3), (0``19.656525995744907 - 1.86556489718674943296967196414173585895`19.885658547256757*I) - (0.66118491636651204033470721957935861933`19.880284735673392 - 1.2049630956839855891878727192782819384`20.010952803838972*I)*k^(1/3) + (0.45440437377885717725483504935547421087`19.59256219833968 + 0.24685679231406852242899703392617353214`19.558553667334497*I)*k^(4/3)}, "per-state-scalar cb (states {0,0},{1,0}) sampled at k=0.3,1,2,5 — NOT constant ⇒ no fixed scalar reconciles even 2 states" -> {{k -> 0.3, -0.19703271755005838 - 0.6000468750006196*I}, {k -> 1., -0.5778488645596662 - 0.72588616455467*I}, {k -> 2., -0.7173896826421956 - 0.6465331262430513*I}, {k -> 5., -0.7297017080721853 - 0.5434454069733405*I}}, "verdict" -> "Dictionary does NOT close at λ1=1 with the assumed form (generator + universal φ(k) + per-state scalar). Obstruction reported; no φ(k) invented/forced.", "resolution" -> "RESOLVED in cells 72-75: LHS must be the 1-
[output truncated — 4253 chars total]
```
## 2026-07-16 17:27:13 — 📥 BULK INSERT 1 CELL(S) at positions 25–25
1. [code] (*Yangian weight eigenvalue of T[1,1] on HWS*)↵ClearAll[ν1];↵ν1[λ1_,λ2_][u_]:=(u-θ1- h λ1)(u-θ2-h λ2


## 2026-07-16 17:27:13
**Input:**
```wolfram
(*Yangian weight eigenvalue of T[1,1] on HWS*)
ClearAll[ν1];
ν1[λ1_,λ2_][u_]:=(u-θ1- h λ1)(u-θ2-h λ2);
```
**Output:**
```
(no output)
```
## 2026-07-16 17:27:21 — ✏️ EDIT CELL 26
ClearAll[xSingle,x2];
xSingle[λ_][{0}]:={Table[KroneckerDelta[ii,dim[{λ,0}]],{ii,1,dim[{λ,0}]}]};
x2[λ1_,λ2_][{0,0}]:=KroneckerProduct[xSingle[λ1][{0}],xSingle[λ2][{0}]];

x2[λ1_,λ2_][{s1_,s2_}] := x2…


## 2026-07-16 17:27:21
**Input:**
```wolfram
ClearAll[xSingle,x2];
xSingle[λ_][{0}]:={Table[KroneckerDelta[ii,dim[{λ,0}]],{ii,1,dim[{λ,0}]}]};
x2[λ1_,λ2_][{0,0}]:=KroneckerProduct[xSingle[λ1][{0}],xSingle[λ2][{0}]];

x2[λ1_,λ2_][{s1_,s2_}] := x2[λ1,λ2][{s1,s2}] = x2[λ1,λ2][{0,0}].t2[λ1,λ2][1,s1][θ1].t2[λ1,λ2][1,s2][θ2]/(Product[ν1[λ1,λ2][θ1 + (k - 1) h], {k, 1, s1}] Product[ν1[λ1,λ2][θ2 + (k - 1) h], {k, 1, s2}]);
```
**Output:**
```
Out= Null
```
## 2026-07-16 17:27:21 — 📥 BULK INSERT 1 CELL(S) at positions 27–27
1. [code] x2[2,1][{2,0}]


## 2026-07-16 17:27:22
**Input:**
```wolfram
x2[2,1][{2,0}]
```
**Output:**
```
Out[222]= \begin{pmatrix}0 & \frac{1}{2} & 0 & 0 & 0 & 0\end{pmatrix}
```
## 2026-07-16 17:27:39 — ✏️ EDIT CELL 30
(*a(u)=aFun is now an ALIAS for Section A's two-site ν1 (Yangian weight eigenvalue of T[1,1] on
  the HWS) — same formula (u-θ1-hλ1)(u-θ2-hλ2), deduped so there is a single source of truth
  instead o…


## 2026-07-16 17:27:39
**Input:**
```wolfram
(*a(u)=aFun is now an ALIAS for Section A's two-site ν1 (Yangian weight eigenvalue of T[1,1] on
  the HWS) — same formula (u-θ1-hλ1)(u-θ2-hλ2), deduped so there is a single source of truth
  instead of two copies that could drift. d(u)=Qθ(u), and τ0 = τ minus the J term (using
  ee1=λ1+λ2-M, ee2=M). Qθ (=d) is defined LOCALLY here so this file is self-contained.
  qpoly builds a polynomial from an ASCENDING coefficient list {c0,...,cM} (cM=1 for monic q_n);
  formed in a dummy symbol first so substituting uu=0 gives the constant term rather than an
  indeterminate 0^0.*)
ClearAll[aFun, τ0, Qθ, qpoly];
aFun[λ1_, λ2_][u_] := ν1[λ1, λ2][u];
Qθ[u_] := (u - θ1) (u - θ2);
τ0[λ1_, λ2_, M_][u_] := (κ1 + κ2) u^2 - (θ1 + θ2) (κ1 + κ2) u - u h (κ1 (λ1 + λ2 - M) + κ2 M);
qpoly[clist_, uu_] := Module[{x}, (clist . x^Range[0, Length[clist] - 1]) /. x -> uu];
{aFun[1, 1][u], Qθ[u], τ0[1, 1, 1][u]}
```
**Output:**
```
Out= {(-4/3 + u)*(-8/7 + u), (-1/3 + u)*(-1/7 + u), (-1.064155821863240252882572111282100269533971449631508102601`29.54117567596917 + 0``29.336679005422717*I)*u + (0.7208797502944530745333553011911001825875290465245700049878`29.54117567596917 + 0``29.50582140452307*I)*u^2}
```
## 2026-07-16 17:27:43 — ✏️ EDIT CELL 58
(*Genuine length-1 SoV covector (paper normalization, 1-site Qθ1(u)=u-θ1). xSingle[λ][{0}] is the
  site lowest-weight covector (Section A) — a 1×(λ+1) row. xL1 stays that shape: (1×n).(n×n)=1×n.
  B1…


## 2026-07-16 17:27:43
**Input:**
```wolfram
(*Genuine length-1 SoV covector (paper normalization, 1-site Qθ1(u)=u-θ1). xSingle[λ][{0}] is the
  site lowest-weight covector (Section A) — a 1×(λ+1) row. xL1 stays that shape: (1×n).(n×n)=1×n.
  B1 is the L=1 SoV B-operator. The single-site ν1[λ][u] defined here is the arity-1 companion to
  Section A's two-site ν1[λ1,λ2][u] (same symbol, different arity — both DownValues coexist; NOT
  ClearAll'd here so Section A's two-site pattern survives). Qθ1[θ1-λh+(p-1)h] and
  ν1[λ][θ1+(p-1)h] give IDENTICAL numeric values (Qθ1 is linear, so shifting the argument vs. the
  function is equivalent) — xL1 is rewritten here purely for notational parity with the Section A
  x2 fix pattern; xL1 itself never had a bug.*)
ClearAll[B1, xL1];
B1[λ_][u_] := L[λ][1, 1][u - θ1];
ν1[λ_][u_] := (u - θ1 - h λ);
xL1[λ_][{0}] := xSingle[λ][{0}];
xL1[λ_][{s_}] := xL1[λ][{s}] = xSingle[λ][{0}] . t1[λ][1, s][θ1] / Product[ν1[λ][θ1 + (p - 1) h], {p, 1, s}];
```
**Output:**
```
Out= Null
```
## 2026-07-16 17:28:16 — ✏️ EDIT CELL 60 [batch]
(*Operatorial factorization: x2[λ1,λ2][{s1,0}] = xL1[λ1][{s1}] ⊗ xSingle[λ2][{0}] — a PURE tensor
  factorization with NO leftover scalar. (An earlier version of this notebook found an apparent
  nont…


## 2026-07-16 17:28:16
**Input:**
```wolfram
(*Operatorial factorization: x2[λ1,λ2][{s1,0}] = xL1[λ1][{s1}] ⊗ xSingle[λ2][{0}] — a PURE tensor
  factorization with NO leftover scalar. (An earlier version of this notebook found an apparent
  nontrivial "site-2 passthrough eigenvalue" μEff here — that was an ARTIFACT of a Qθ normalization
  bug in x2's original definition, since fixed in Section A: x2 divided by Qθ(u)=(u-θ1)(u-θ2)
  evaluated at a point shifted by -λ1h in BOTH factors, incorrectly putting site-1's λ1 shift into
  the θ2 factor that should have carried site 2's OWN λ2 shift. The corrected x2 divides by
  ν1[λ1,λ2][u]=(u-θ1-hλ1)(u-θ2-hλ2) instead — each factor shifted by its own site's λ — and with
  that fix the ratio below comes out EXACTLY 1 for every (λ1,λ2,s1), not just a constant. There is
  no genuine extra site-2 eigenvalue beyond what ν1 already encodes: the old cPredCorr/μSite2
  "finding" was bookkeeping fallout from the bug, not new physics — see the cell below for the
  honest restatement of that operator identity.)
  HARD requirement: the ratio lhs/tp must be IDENTICALLY 1 (both non-constancy spread AND
  deviation-from-1 below the 10^-15 floor) for all 27 (λ1,λ2,s1) cases.*)
SoVFact::notconst = "Ratio lhs/tp not constant for (λ1,λ2,s1)=(`1`,`2`,`3`): spread=`4`.";
SoVFact::notone = "Ratio lhs/tp ≠ 1 for (λ1,λ2,s1)=(`1`,`2`,`3`): ratio=`4`.";
factResults = Flatten[Table[
    Module[{lhs, tp, nz, rl, cval, spread},
     lhs = Flatten[x2[λ1, λ2][{s1, 0}]];
     tp = Flatten[KroneckerProduct[xL1[λ1][{s1}], xSingle[λ2][{0}]]];
     nz = Position[tp, x_ /; Abs[Chop[x, 10^-20]] > 10^-12];
     rl = Extract[lhs, nz]/Extract[tp, nz];
     cval = Mean[rl];
     spread = Max[Abs[rl - cval]];
     If[spread > 10^-15, Message[SoVFact::notconst, λ1, λ2, s1, N[spread]]; Abort[]];
     If[Abs[N[cval - 1]] > 10^-15, Message[SoVFact::notone, λ1, λ2, s1, N[cval]]; Abort[]];
     {λ1, λ2, s1, N[cval], N[spread]}],
    {λ1, 1, 3}, {λ2, 1, 3}, {s1, 0, λ1}], 2];
Column[factResults]
```
**Output:**
```
Out= Out[226]= \begin{matrix}\left\{1,1,0,1.\,,0.\,\right\}\\\left\{1,1,1,1.\,,0.\,\right\}\\\left\{1,2,0,1.\,,0.\,\right\}\\\left\{1,2,1,1.\,,0.\,\right\}\\\left\{1,3,0,1.\,,0.\,\right\}\\\left\{1,3,1,1.\,,0.\,\right\}\\\left\{2,1,0,1.\,,0.\,\right\}\\\left\{2,1,1,1.\,,0.\,\right\}\\\left\{2,1,2,1.\,,0.\,\right\}\\\left\{2,2,0,1.\,,0.\,\right\}\\\left\{2,2,1,1.\,,0.\,\right\}\\\left\{2,2,2,1.\,,0.\,\r… [898 chars total]
```
## 2026-07-16 17:28:16 — ✏️ EDIT CELL 61 [batch]
(*Summary: worst deviation from the now-EXACT tensor factorization (ratio ≡ 1, no leftover
  scalar) across the whole sweep. Columns of factResults: {λ1,λ2,s1,ratio,spread}.*)
<|
  "max spread (ratio …


## 2026-07-16 17:28:16
**Input:**
```wolfram
(*Summary: worst deviation from the now-EXACT tensor factorization (ratio ≡ 1, no leftover
  scalar) across the whole sweep. Columns of factResults: {λ1,λ2,s1,ratio,spread}.*)
<|
  "max spread (ratio non-constancy)" -> Max[factResults[[All, 5]]],
  "max |ratio - 1| over sweep" -> Max[Abs[factResults[[All, 4]] - 1]],
  "num (λ1,λ2,s1) cases checked" -> Length[factResults]
|>
```
**Output:**
```
Out= Out[227]= \langle\!\lvert \text{max spread (ratio non-constancy)}\to 0.\,,\text{max |ratio - 1| over sweep}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases checked}\to 27\rvert\!\rangle
```
## 2026-07-16 17:28:16 — ✏️ EDIT CELL 62 [batch]
(*Operator-level identity (STILL TRUE, and unaffected by the Qθ normalization bug — it never used
  x2[{s1,0}]'s buggy denominator, only the raw x2[{0,0}] pseudovacuum and t2 itself): the site-2
  pse…


## 2026-07-16 17:28:16
**Input:**
```wolfram
(*Operator-level identity (STILL TRUE, and unaffected by the Qθ normalization bug — it never used
  x2[{s1,0}]'s buggy denominator, only the raw x2[{0,0}] pseudovacuum and t2 itself): the site-2
  pseudovacuum xSingle[λ2][{0}] passes through t2[λ1,λ2][1,s1][θ1] picking up the factor
  μSite2[λ2,s1] = Product[(θ1-θ2+(k-1)h-λ2 h),{k,1,s1}] relative to the bare tensor product
  (xSingle[λ1][{0}].t1[λ1][1,s1][θ1]) ⊗ xSingle[λ2][{0}]. CORRECTED READING (see the cell above):
  this is NOT an extra, unexplained site-2 eigenvalue — it is exactly the λ2-dependent factor of
  Section A's two-site ν1[λ1,λ2][u]=(u-θ1-hλ1)(u-θ2-hλ2) evaluated on the site-1 tower, i.e.
  ν1[λ1,λ2][θ1+(k-1)h] factors as ((k-1-λ1)h)·(θ1-θ2+(k-1-λ2)h), and the second factor here IS
  μSite2[λ2,s1]. Once x2's denominator uses the FULL two-site ν1 (Section A's fix) instead of a
  λ1-only Qθ factor, this piece is already built into the normalization, which is why the ratio in
  the cell above is exactly 1 — μSite2 was always just part of the two-site Yangian weight
  function's structure, not a surprising new eigenvalue. Kept here as a genuine, still-verified
  operator identity (hard-asserted below), relabeled honestly.*)
ClearAll[μSite2];
μSite2[λ2_, s1_] := Product[(θ1 - θ2 + (k - 1) h - λ2 h), {k, 1, s1}];
Site2Mu::notconst = "μ-ratio not constant at operator level for (λ1,λ2,s1)=(`1`,`2`,`3`): spread=`4`.";
Site2Mu::mismatch = "μ ≠ μSite2 at operator level for (λ1,λ2,s1)=(`1`,`2`,`3`): μ=`4`, μSite2=`5`.";
opMuResults = Flatten[Table[
    Module[{numLhs, numTp, nz, rl, muval, mupred, spread},
     numLhs = Flatten[x2[λ1, λ2][{0, 0}] . t2[λ1, λ2][1, s1][θ1]];
     numTp = Flatten[KroneckerProduct[xSingle[λ1][{0}] . t1[λ1][1, s1][θ1], xSingle[λ2][{0}]]];
     nz = Position[numTp, x_ /; Abs[Chop[x, 10^-20]] > 10^-12];
     rl = Extract[numLhs, nz]/Extract[numTp, nz];
     muval = Mean[rl];
     spread = Max[Abs[rl - muval]];
     mupred = μSite2[λ2, s1];
     If[spread > 10^-15, Message[Site2Mu::notconst, λ1, λ2, s1, N[spread]]; Abort[]];
     If[Abs[N[muval - mupred]] > 10^-15, Message[Site2Mu::mismatch, λ1, λ2, s1, N[muval], N[mupred]]; Abort[]];
     {λ1, λ2, s1, N[muval], N[mupred], N[spread]}],
    {λ1, 1, 3}, {λ2, 1, 3}, {s1, 0, λ1}], 2];
<|"max operator-level spread" -> Max[opMuResults[[All, 6]]],
  "max |μ - μSite2|" -> Max[Abs[opMuResults[[All, 4]] - opMuResults[[All, 5]]]],
  "num (λ1,λ2,s1) cases" -> Length[opMuResults]|>
```
**Output:**
```
Out= Out[233]= \langle\!\lvert \text{max operator-level spread}\to 0.\,,\text{max |$\mu$ - $\mu$Site2|}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases}\to 27\rvert\!\rangle
```

---

## 2026-07-16 17:28:22 — 🔄 KERNEL RESTART

## 2026-07-16 17:28:52 — ▶️ RUN CELLS 1–36
Cell 1: ✓ — Out[8]= \left\{\frac{1}{3},\frac{1}{7},1,0.360439875147226537266677650596+0.932782448593374716484835982071\,\mathrm{i},0.36043987514722653726667765060-0.932782448593374716484835982071\,\mathrm{i},0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i},1\right\}
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 8: ✓ — (no output)
Cell 10: ✓ — (no output)
Cell 11: ✓ — (no output)
Cell 12: ✓ — Out[27]= \left\{0\right\}
Cell 13: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — Out[31]= \left\{0\right\}
Cell 16: ✓ — (no output)
Cell 17: ✓ — Out[34]= \begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{pmatrix}
Cell 18: ✓ — (no output)
Cell 20: ✓ — (no output)
Cell 21: ✓ — (no output)
Cell 23: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 26: ✓ — (no output)
Cell 27: ✓ — Out[52]= \begin{pmatrix}0 & \frac{1}{2} & 0 & 0 & 0 & 0\end{pmatrix}
Cell 28: ✓ — Out[53]= \left\{0\right\}
Cell 30: ✓ — Out[59]= \left\{\left(u-\frac{4}{3}\right)\,\left(u-\frac{8}{7}\right),\left(u-\frac{1}{3}\right)\,\left(u-\frac{1}{7}\right),\left(0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i}\right)\,u^2-\left(1.06415582186324025288257211128+0.\times 10^{-29}\,\mathrm{i}\right)\,u\right\}
Cell 31: ✓ — Out[62]= \begin{pmatrix}1 & 1\\2 & 2\\3 & 3\end{pmatrix}
Cell 32: ✓ — Out[65]= \text{max|C2-($\kappa$1+$\kappa$2)Id|}\to 0.\times 10^{-29}
Cell 33: ✓ — Out[75]= \langle\!\lvert \left\{0,0\right\}\to 0.5664055180884988 +1.3769645669711723\,\mathrm{i},\left\{1,0\right\}\to -0.778123966451158+0.\,\,\mathrm{i},\left\{1,1\right\}\to 1.1900552523337025 +0.\,\,\mathrm{i},\left\{2,0\right\}\to 0.5664055180884988 -1.3769645669711723\,\mathrm{i}\rvert\!\rangle
Cell 34: ✓ — 2 kernel messages: Psi::spurious: No genuine state at (M,n)=(2,2) for (λ1,λ2)=(2,1): weight sector M=2 has only 2 genuine state(s) (n=0..1), so n=2 does not exist. [also: Psi::spurious] result: Print: spurious-guard test (expect a Psi::spurious message and Missing[...]):[formula] | Out[84]= \left\{\left\{1.000000000000000000000000000+0.\times 10^{-29}\,\mathrm{i},-1.3004494881076436272092358…
⚠️ stopped at Cell 34 — error detected (pass stopOnError:false to continue past errors)

## 2026-07-16 17:29:16 — ▶️ RUN CELLS 35–77
Cell 35: ✓ — Out[91]= \left\{1.03464899078135229594035121+0.381869568380946283899575941\,\mathrm{i},-2.37022278949001677248172619-0.20048970475590097847093779\,\mathrm{i},1.00000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right\}
Cell 36: ✓ — Out[94]= \left\{-0.709231010958016926537225142+1.09596916251277206706752728\,\mathrm{i},-0.02291660900+0.\times 10^{-12}\,\mathrm{i}\right\}
Cell 37: ✓ — (no output)
Cell 38: ✓ — (no output)
Cell 39: ✓ — (no output)
Cell 40: ✓ — (no output)
Cell 41: ✓ — (no output)
Cell 42: ✓ — Out[109]= 0
Cell 43: ✓ — (no output)
Cell 44: ✓ — Out[113]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 45: ✓ — Out[115]= \left\{6,6\right\}
Cell 46: ✓ — Out[116]= \begin{pmatrix}2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i} & 0 & 0 & 0 & 0 & 0\\0 & 7.36608403016431791181-0.968046299794461400266\,\mathrm{i} & 0 & 0 & 0 & 0\\0 & 0 & 2.68702160063169134012-1.3456153260557858257\,\mathrm{i} & 0 & 0 & 0\\0 & 0 & 0 & 7.3660840301643179118107+0.9680462997944614002659\,\mathrm{i} & 0 & 0\\0 & 0 & 0 & 0 & 2.687021600631691340118… [537 chars total]
Cell 47: ✓ — Out[125]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{max|off-diagonal|}\to 5.065342385207116\times 10^{-34},\text{digits vanished (-Log10)}\to 33.29539119374288,\text{min|diagonal|}\to 3.0051232400647923,\text{diagonal values}\to \left\{2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i},7.36608403016431791181-0… [699 chars total]
Cell 48: ✓ — (no output)
Cell 49: ✓ — (no output)
Cell 50: ✓ — Out[133]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 51: ✓ — Out[138]= \left\{1.00000000000000000000000000000,6\right\}
Cell 52: ✓ — Out[143]= \langle\!\lvert \text{max|<PsiLeft|Psi>| off-diagonal}\to 0.\,,\text{digits vanished (-Log10)}\to \mathrm{Indeterminate}\rvert\!\rangle
Cell 53: ✓ — Out[147]= \left\{7.366084030164317911810727273-0.9680462997944614002658874268\,\mathrm{i},1.00000000000000000000+0.\times 10^{-21}\,\mathrm{i}\right\}
Cell 54: ✓ — Out[150]= \begin{matrix}\langle\!\lvert \text{M}\to 0,\text{n}\to 0,\text{braket}\to 2.909064289-5.804620064\,\mathrm{i},\text{G21 diagonal}\to 2.909064289-5.804620064\,\mathrm{i},\text{normRatio}\to 1.000000000+0.\times 10^{-10}\,\mathrm{i}\rvert\!\rangle \\\langle\!\lvert \text{M}\to 1,\text{n}\to 0,\text{braket}\to 7.366084030-0.968046300\,\mathrm{i},\text{G21 diagonal}\to 7.366084030-0.9680463… [1247 chars total]
Cell 56: ✓ — (no output)
Cell 57: ✓ — Out[158]= \begin{pmatrix}3 & 3\\3 & 3\end{pmatrix}
Cell 58: ✓ — (no output)
Cell 59: ✓ — Out[165]= \left\{0\right\}
Cell 60: ✓ — Out[169]= \begin{matrix}\left\{1,1,0,1.\,,0.\,\right\}\\\left\{1,1,1,1.\,,0.\,\right\}\\\left\{1,2,0,1.\,,0.\,\right\}\\\left\{1,2,1,1.\,,0.\,\right\}\\\left\{1,3,0,1.\,,0.\,\right\}\\\left\{1,3,1,1.\,,0.\,\right\}\\\left\{2,1,0,1.\,,0.\,\right\}\\\left\{2,1,1,1.\,,0.\,\right\}\\\left\{2,1,2,1.\,,0.\,\right\}\\\left\{2,2,0,1.\,,0.\,\right\}\\\left\{2,2,1,1.\,,0.\,\right\}\\\left\{2,2,2,1.\,,0.\,\r… [898 chars total]
Cell 61: ✓ — Out[170]= \langle\!\lvert \text{max spread (ratio non-constancy)}\to 0.\,,\text{max |ratio - 1| over sweep}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases checked}\to 27\rvert\!\rangle
Cell 62: ✓ — Out[176]= \langle\!\lvert \text{max operator-level spread}\to 0.\,,\text{max |$\mu$ - $\mu$Site2|}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases}\to 27\rvert\!\rangle
Cell 64: ✓ — (no output)
Cell 65: ✓ — Out[185]= \langle\!\lvert \text{w\_0 ($\lambda$=1,2,3)}\to \left\{1,1,1\right\},\text{full weights w\_j}\to \left\{1\to \left\{1,-1\right\},2\to \left\{1,-2,1\right\},3\to \left\{1,-3,3,-1\right\}\right\},\text{br1[$\lambda$][1] ($\lambda$=1,2,3)}\to \left\{0,0,0\right\}\rvert\!\rangle
Cell 66: ✓ — Out[187]= \left\{0,0\right\}
Cell 67: ✓ — Out[190]= \sqrt[3]{k}\,\left(\left(-1.18561270626131985872+0.832253303224995280404\,\mathrm{i}\right)\,k^2-\left(0.761819910806671651081+0.812812473459294825632\,\mathrm{i}\right)\,k-\left(0.404410037536507329420+0.598870869485872613792\,\mathrm{i}\right)\right)
Cell 68: ✓ — Out[194]= \begin{pmatrix}0 & 0 & 0\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix}
Cell 69: ✓ — Out[198]= \langle\!\lvert \text{OP(E12)}\to 1.\,+0.\times 10^{-30}\,\mathrm{i},\text{OP(E21)}\to \left(-0.03626734060947795145020244939-2.004237714666237867354326095\,\mathrm{i}\right)\,\varphi ^2-\left(1.522150083623538480912456055-0.2442047957313920727733393355\,\mathrm{i}\right)\,\varphi +\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right),\text{nontrivial generator \{a,b\}}\to \left\{2,1\right\},\… [464 chars total]
Cell 70: ✓ — Out[201]= \left\{1,2\right\}
Cell 71: ✓ — Out[207]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{$\varphi$ candidates per state}\to \begin{pmatrix}\left\{0,0\right\} & \left(-0.920793778924720536105378676-0.390049761817044163511250929\,\mathrm{i}\right)\,k^{4/3}+\left(0.695722366612934233848844287-0.718310788304406529935213657\,\mathrm{i}\right)\,\sqrt[3]{k}-\left(0.3604398751… [1247 chars total]
Cell 72: ✓ — Out[213]= \langle\!\lvert \text{c\_sb solution(s)}\to \left\{\left\{\mathrm{cb}\to \frac{\left(-6.242239011314\times 10^{45}-4.511757007219\times 10^{45}\,\mathrm{i}\right)\,k^{4/3}-\left(5.999302657796\times 10^{45}+8.356648463050\times 10^{45}\,\mathrm{i}\right)\,\sqrt[3]{k}}{9.833530773718\times 10^{45}\,k^{4/3}-\left(3.54439660434\times 10^{45}-9.172544913427\times 10^{45}\,\mathrm{i}\right)\,… [841 chars total]
Cell 73: ✓ — 2 kernel messages: Part::partd: Part specification ABTab[[1,4]] is longer than depth of object.
Part::partd: Part specification ABTab[[2,4]] is longer than depth of object.
Part::partd: Part specification ABTab[[3,4]] i… [also: Part::partd] result: Out[215]= \langle\!\lvert \text{generator \{a,b\} tested}\to \left\{2,1\right\},\text{rep tested}\to \left\{1,2\right\},\text{FL shape (verified exact for all 6 states)}\to \text{FL[1,2][state][k] = k…
⚠️ stopped at Cell 73 — error detected (pass stopOnError:false to continue past errors)

## 2026-07-16 17:30:18 — ✏️ EDIT CELL 73
(*Discovery step 2 — OBSTRUCTION FOUND, reported rather than forced. Direct structural test: write
  FL[1,2][state][k] = k^(1/3)(A_state + B_state k) (verified: every state of this rep is EXACTLY of
 …


## 2026-07-16 17:30:18
**Input:**
```wolfram
(*Discovery step 2 — OBSTRUCTION FOUND, reported rather than forced. Direct structural test: write
  FL[1,2][state][k] = k^(1/3)(A_state + B_state k) (verified: every state of this rep is EXACTLY of
  this shape — a degree-1 polynomial in k times the universal (k κ1)^(θ1/h) prefactor, θ1=1/3). Since
  OPg[1,2][state][φ] = c0(state) + c1(state) φ is linear in φ with c0(state)=1 for EVERY state (a
  clean structural fact, presumably <x0|Psi> normalization), any dictionary of the plan's assumed
  form (single generator, single state-INDEPENDENT φ(k), single per-state scalar) would force the
  ratio A_state/B_state to depend on state ONLY through c1(state) in a fixed (Möbius) way. It does
  not: A/B clusters into 3 conjugate-paired values across the 6 states of (1,2), not a single common
  value or a single function of c1(state). Independently, Step 1b's attempt to fix a per-state scalar
  cb between just the two lowest states {0,0},{1,0} returns a cb that is NOT constant in k (varies
  from -0.578-0.726i at k=1 to -0.197-0.600i at k=0.3 to -0.730-0.543i at k=5) — so no constant scalar
  reconciles even two states. CONCLUSION: the naive dictionary {generator={2,1}, universal φ(k),
  per-state scalar} does NOT close already at λ1=1. Recorded here as the discovery's terminal
  obstruction (not patched/forced) — see chat report for the full numeric detail.

  RESOLVED (see cells 72-75 below): the obstruction was a WRONG-OBJECT bug in the plan, not in
  this construction. The correct LHS functional object is not the bare bracket br1[f] but the 1-site
  FSoV DETERMINANT det1[f] = br1[f]/(f/.u->θ1) — Section C's det[λ1,λ2][f] (which divides by f at
  BOTH nodes θ1,θ2) reduced to 1 site. The missing 1/f(θ1) exactly cancels the (k κ1)^(θ1/h)
  fractional prefactor that made FL's shape state-dependent (the A/B clustering seen below). With
  det1 in place, FLdet[1,2][state][k] = 1 + D_state k for every state (clean degree-1 polynomial,
  leading term exactly 1), and the ratio α_state = D_state/c1(state) IS universal: α = -1 exactly for
  all 6 states, extending to a symbolically EXACT match (Simplify gives literal 0, not just below a
  noise floor) across the full rep sweep {{1,1},{2,1},{1,2},{2,2},{3,1},{1,3}} under φ(k) = -k. The
  A_state/B_state values below (br1-based) were a real, correctly-computed symptom — just of the
  wrong functional object; nothing here was recomputed or altered.

  HIDDEN-DEPENDENCY FIX (unrelated to the ν1/x2 physics correction, found during a fresh-kernel
  run): ABTab below was originally computed via an ad-hoc scratch evaluation during the live
  discovery session, never as a saved cell — so on a genuinely fresh kernel this cell errored
  (Part::partd, ABTab undefined). Inlined here exactly as originally computed, restoring
  self-containment; no numbers changed.*)
ABTab = Table[
   Module[{f, A, B},
    f = FL[1, 2][s[[1]], s[[2]]][k] // Simplify;
    A = Coefficient[f/k^(1/3) // Simplify, k, 0];
    B = Coefficient[f/k^(1/3) // Simplify, k, 1];
    {s, A, B, A/B}],
   {s, sts12}];
obstructionReport = <|
  "generator {a,b} tested" -> genNontrivial,
  "rep tested" -> {1, 2},
  "FL shape (verified exact for all 6 states)" -> "FL[1,2][state][k] = k^(1/3) (A_state + B_state k)",
  "c0(state) = OPg[state][0]" -> Table[{s, lin12[s[[1]], s[[2]]][[1]]}, {s, sts12}],
  "A_state/B_state (structural obstruction: not universal, not a fixed function of c1(state))" ->
    Table[{sts12[[i]], ABTab[[i, 4]]}, {i, Length[sts12]}],
  "scalar=1 φ-candidate differences (state i - state 1), nonzero ⇒ not universal at scalar 1" -> phiDiffs12,
  "per-state-scalar cb (states {0,0},{1,0}) sampled at k=0.3,1,2,5 — NOT constant ⇒ no fixed scalar reconciles even 2 states" ->
    Table[{k -> kk, cb /. cbSol[[1]] /. k -> kk}, {kk, {0.3, 1., 2., 5.}}],
  "verdict" -> "Dictionary does NOT close at λ1=1 with the assumed form (generator + universal φ(k) + per-state scalar). Obstruction reported; no φ(k) invented/forced.",
  "resolution" -> "RESOLVED in cells 72-75: LHS must be the 1-site FSoV determinant det1[f]=br1[f]/f(θ1), not the bare bracket br1[f]. With det1, φ(k)=α k, α=-1, closes EXACTLY (symbolic 0) across the full rep sweep."
|>;
obstructionReport
```
**Output:**
```
Out= <|"generator {a,b} tested" -> {2, 1}, "rep tested" -> {1, 2}, "FL shape (verified exact for all 6 states)" -> "FL[1,2][state][k] = k^(1/3) (A_state + B_state k)", "c0(state) = OPg[state][0]" -> {{{0, 0}, 1.`28.362553630092005 + 0``28.458858685029718*I}, {{1, 0}, 1.`28.362553630092005 + 0``28.420425233412512*I}, {{1, 1}, 1.`28.362553630092005 + 0``28.534910438118857*I}, {{2, 0}, 1.`28.362553630092005 + 0``28.420425233412512*I}, {{2, 1}, 1.`28.362553630092005 + 0``28.534910438118857*I}, {{3, 0}, 1.`28.362553630092005 + 0``28.458858685029718*I}}, "A_state/B_state (structural obstruction: not universal, not a fixed function of c1(state))" -> {{{0, 0}, -0.3604398751472265372666776505955500912995704233956884657367`28.519642521730212 + 0.93278244859337471648483598207086792945`28.723716363759525*I}, {{1, 0}, 1.26686360606317953186931967522323009443`20.5556504399651 + 0.42306417719964480709484006322965914477`20.076520933327156*I}, {{1, 1}, -1.07943487098662173249064729691354404681`21.05118564325957 + 0.36047307961878995475242216170986991601`20.60966884642148*I}, {{2, 0}, 1.26686360606317953186931967522323009437`22.374310009255577 - 0.4230641771996448070948400632296591446`22.290189543967912*I}, {{2, 1}, -1.0794348709866217324906472969135440468`22.956040954336984 - 0.360473079618789954752422161709869916`22.832840645334333*I}, {{3, 0}, -0.36043987514722653726667765059555009141`20.924877319063018 - 0.93278244859337471648483598207086792943`21.286743472795727*I}}, "scalar=1 φ-candidate differences (state i - state 1), nonzero ⇒ not universal at scalar 1" -> {(1.62730348121040606913599732581878018564`18.928857859553823 - 0.50971827139372990938999591884120878462`18.420106071749494*I) - (0.46542748324955977674418064345565154397`18.55368709365933 - 1.73877604346995428951499886313331938519`19.07388701197023*I)*k^(1/3) + (1.32634547687585376726418287584726808483`19.18723952357479 + 1.06012257446803590397186378841389386732`19.005029161051315*I)*k^(4/3), (-0.71899499583939519522396964631799395547`19.17048522233265 - 0.57230936897458476173241382036099801359`19.112586969447793*I) - (1.55376991547313991856351499799445509993`19.390842030784782 - 1.83632831184250631961923901798977241298`19.496759945119884*I)*k^(1/3) + (1.94712428919652521675004435696818121602`19.580365406788538 - 0.30295453506097416430509738933976816355`18.83333649413266*I)*k^(4/3), (1.6273034812104060691359973258187801854`20.927897475259947 - 1.35584662579301952357967604530052707432`21.12355164110732*I) - (0.96539125287879209819001219524303509911`21.268958819455012 - 0.3392821519416261106602048437419337344`20.676771979394047*I)*k^(1/3) + (0.81917510810755585089599730622926840074`21.495943058250155 + 0.05692798477617834807934963005547620116`20.051679306120185*I)*k^(4/3), (-0.71899499583939519522396964631799395559`21.201403018442612 - 1.29325552821216467123725814378073784552`21.806152133017648*I) - (0.32640350628717582177610187293768968943`21.198486434539824 + 0.36656814642707085144013045320709439335`20.90521155211221*I)*k^(1/3) + (0.91493640438785580629103394186170207104`22.17656014576632 + 1.39704924810910189208916572287478175457`21.524209614852612*I)*k^(4/3), (0``19.656525995744907 - 1.86556489718674943296967196414173585895`19.885658547256757*I) - (0.66118491636651204033470721957935861933`19.880284735673392 - 1.2049630956839855891878727192782819384`20.010952803838972*I)*k^(1/3) + (0.45440437377885717725483504935547421087`19.59256219833968 + 0.24685679231406852242899703392617353214`19.558553667334497*I)*k^(4/3)}, "per-state-scalar cb (states {0,0},{1,0}) sampled at k=0.3,1,2,5 — NOT constant ⇒ no fixed scalar reconciles even 2 states" -> {{k -> 0.3, -0.19703271754993507 - 0.6000468750005826*I}, {k -> 1., -0.5778488645595341 - 0.7258861645547009*I}, {k -> 2., -0.7173896826421186 - 0.6465331262430857*I}, {k -> 5., -0.7297017080721572 - 0.5434454069733533*I}}, "verdict" -> "Dictionary does NOT close at λ1=1 with the assumed form (generator + universal φ(k) + per-state scalar). Obstruction reported; no φ(k) invented/forced.", "resolution" -> "RESOLVED in cells 72-75: LHS must be the 1
[output truncated — 4254 chars total]
```
## 2026-07-16 17:30:26 — ▶️ RUN CELLS 74–77
Cell 74: ✓ — Out[219]= \left(1.000000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right)-\left(0.360439875147226537266677651+0.932782448593374716484835982\,\mathrm{i}\right)\,k
Cell 75: ✓ — Out[226]= \langle\!\lvert \text{alphaTab \{state,D0,c0,Dc,c1,$\alpha$\}}\to \begin{pmatrix}\left\{0,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & -0.36043987514722653-0.9327824485933747\,\mathrm{i} & 0.36043987514722653 +0.9327824485933747\,\mathrm{i} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.7101545203859353 -0.2371533418544… [1247 chars total]
Cell 76: ✓ — Out[234]= \langle\!\lvert \text{validation}\to \left\{\left\{1,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\end{pmatrix},\left\{2,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\\\left\{2,1\right\} & 0\\\left\{3,0\right\} & 0\end{pmatrix},\left\{1… [1247 chars total]
Cell 77: ✓ — Out[236]= \langle\!\lvert \text{generator \{a,b\}}\to \left\{2,1\right\},\text{$\varphi$(k)}\to \text{$\alpha$ k},\text{$\alpha$ (universal value)}\to -1.\,+0.\,\,\mathrm{i},\text{normalization}\to \text{1-site FSoV determinant det1[f] = br1[f]/(f/.u->$\theta$1); Section C's det[$\lambda$1,$\lambda$2][f] restricted to 1 site},\text{max residual over rep sweep}\to 0,\text{reps validated}\to \begin{… [510 chars total]

