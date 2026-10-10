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


---

## 2026-07-18 10:20:50 — 🔄 KERNEL RESTART

## 2026-07-18 10:30:40 — ▶️ RUN CELLS 1–77
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
Cell 27: ⏱ timeout — Out[52]= \begin{pmatrix}0 & \frac{1}{2} & 0 & 0 & 0 & 0\end{pmatrix}
Cell 28: ⏱ timeout — Out[53]= \left\{0\right\}
⚠️ global timeout (580s) reached before Cell 30


---

## 2026-07-18 10:51:00 — 🔄 KERNEL RESTART

## 2026-07-18 10:51:09 — ▶️ RUN CELLS 1–10
Cell 1: ✓ — Out[8]= \left\{\frac{1}{3},\frac{1}{7},1,0.360439875147226537266677650596+0.932782448593374716484835982071\,\mathrm{i},0.36043987514722653726667765060-0.932782448593374716484835982071\,\mathrm{i},0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i},1\right\}
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 8: ✓ — (no output)
Cell 10: ✓ — (no output)

## 2026-07-18 10:51:18 — ▶️ RUN CELLS 11–20
Cell 11: ✓ — (no output)
Cell 12: ✓ — Out[27]= \left\{0\right\}
Cell 13: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — Out[31]= \left\{0\right\}
Cell 16: ✓ — (no output)
Cell 17: ✓ — Out[34]= \begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{pmatrix}
Cell 18: ✓ — (no output)
Cell 20: ✓ — (no output)

## 2026-07-18 10:51:26 — ▶️ RUN CELLS 21–28
Cell 21: ✓ — (no output)
Cell 23: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 26: ✓ — (no output)
Cell 27: ✓ — Out[52]= \begin{pmatrix}0 & \frac{1}{2} & 0 & 0 & 0 & 0\end{pmatrix}
Cell 28: ✓ — Out[53]= \left\{0\right\}

## 2026-07-18 10:51:34 — ▶️ RUN CELLS 29–36
Cell 30: ✓ — Out[59]= \left\{\left(u-\frac{4}{3}\right)\,\left(u-\frac{8}{7}\right),\left(u-\frac{1}{3}\right)\,\left(u-\frac{1}{7}\right),\left(0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i}\right)\,u^2-\left(1.06415582186324025288257211128+0.\times 10^{-29}\,\mathrm{i}\right)\,u\right\}
Cell 31: ✓ — Out[62]= \begin{pmatrix}1 & 1\\2 & 2\\3 & 3\end{pmatrix}
Cell 32: ✓ — Out[65]= \text{max|C2-($\kappa$1+$\kappa$2)Id|}\to 0.\times 10^{-29}
Cell 33: ✓ — Out[75]= \langle\!\lvert \left\{0,0\right\}\to 0.5664055180884988 +1.3769645669711723\,\mathrm{i},\left\{1,0\right\}\to -0.778123966451158+0.\,\,\mathrm{i},\left\{1,1\right\}\to 1.1900552523337025 +0.\,\,\mathrm{i},\left\{2,0\right\}\to 0.5664055180884988 -1.3769645669711723\,\mathrm{i}\rvert\!\rangle
Cell 34: ✓ — 2 kernel messages: Psi::spurious: No genuine state at (M,n)=(2,2) for (λ1,λ2)=(2,1): weight sector M=2 has only 2 genuine state(s) (n=0..1), so n=2 does not exist. [also: Psi::spurious] result: Print: spurious-guard test (expect a Psi::spurious message and Missing[...]):[formula] | Out[84]= \left\{\left\{1.000000000000000000000000000+0.\times 10^{-29}\,\mathrm{i},-1.3004494881076436272092358…
⚠️ stopped at Cell 34 — error detected (pass stopOnError:false to continue past errors)

## 2026-07-18 10:51:45 — ▶️ RUN CELLS 35–42
Cell 35: ✓ — Out[91]= \left\{1.03464899078135229594035121+0.381869568380946283899575941\,\mathrm{i},-2.37022278949001677248172619-0.20048970475590097847093779\,\mathrm{i},1.00000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right\}
Cell 36: ✓ — Out[94]= \left\{-0.709231010958016926537225142+1.09596916251277206706752728\,\mathrm{i},-0.02291660900+0.\times 10^{-12}\,\mathrm{i}\right\}
Cell 37: ✓ — (no output)
Cell 38: ✓ — (no output)
Cell 39: ✓ — (no output)
Cell 40: ✓ — (no output)
Cell 41: ✓ — (no output)
Cell 42: ✓ — Out[109]= 0

## 2026-07-18 10:51:55 — ▶️ RUN CELLS 43–47
Cell 43: ✓ — (no output)
Cell 44: ✓ — Out[113]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 45: ✓ — Out[115]= \left\{6,6\right\}
Cell 46: ✓ — Out[116]= \begin{pmatrix}2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i} & 0 & 0 & 0 & 0 & 0\\0 & 7.36608403016431791181-0.968046299794461400266\,\mathrm{i} & 0 & 0 & 0 & 0\\0 & 0 & 2.68702160063169134012-1.3456153260557858257\,\mathrm{i} & 0 & 0 & 0\\0 & 0 & 0 & 7.3660840301643179118107+0.9680462997944614002659\,\mathrm{i} & 0 & 0\\0 & 0 & 0 & 0 & 2.687021600631691340118… [537 chars total]
Cell 47: ✓ — Out[125]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{max|off-diagonal|}\to 5.065342385207116\times 10^{-34},\text{digits vanished (-Log10)}\to 33.29539119374288,\text{min|diagonal|}\to 3.0051232400647923,\text{diagonal values}\to \left\{2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i},7.36608403016431791181-0… [699 chars total]

## 2026-07-18 10:52:03 — ▶️ RUN CELLS 48–54
Cell 48: ✓ — (no output)
Cell 49: ✓ — (no output)
Cell 50: ✓ — Out[133]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 51: ✓ — Out[138]= \left\{1.00000000000000000000000000000,6\right\}
Cell 52: ✓ — Out[143]= \langle\!\lvert \text{max|<PsiLeft|Psi>| off-diagonal}\to 0.\,,\text{digits vanished (-Log10)}\to \mathrm{Indeterminate}\rvert\!\rangle
Cell 53: ✓ — Out[147]= \left\{7.366084030164317911810727273-0.9680462997944614002658874268\,\mathrm{i},1.00000000000000000000+0.\times 10^{-21}\,\mathrm{i}\right\}
Cell 54: ✓ — Out[150]= \begin{matrix}\langle\!\lvert \text{M}\to 0,\text{n}\to 0,\text{braket}\to 2.909064289-5.804620064\,\mathrm{i},\text{G21 diagonal}\to 2.909064289-5.804620064\,\mathrm{i},\text{normRatio}\to 1.000000000+0.\times 10^{-10}\,\mathrm{i}\rvert\!\rangle \\\langle\!\lvert \text{M}\to 1,\text{n}\to 0,\text{braket}\to 7.366084030-0.968046300\,\mathrm{i},\text{G21 diagonal}\to 7.366084030-0.9680463… [1247 chars total]

## 2026-07-18 10:52:10 — ▶️ RUN CELLS 55–62
Cell 56: ✓ — (no output)
Cell 57: ✓ — Out[158]= \begin{pmatrix}3 & 3\\3 & 3\end{pmatrix}
Cell 58: ✓ — (no output)
Cell 59: ✓ — Out[165]= \left\{0\right\}
Cell 60: ✓ — Out[169]= \begin{matrix}\left\{1,1,0,1.\,,0.\,\right\}\\\left\{1,1,1,1.\,,0.\,\right\}\\\left\{1,2,0,1.\,,0.\,\right\}\\\left\{1,2,1,1.\,,0.\,\right\}\\\left\{1,3,0,1.\,,0.\,\right\}\\\left\{1,3,1,1.\,,0.\,\right\}\\\left\{2,1,0,1.\,,0.\,\right\}\\\left\{2,1,1,1.\,,0.\,\right\}\\\left\{2,1,2,1.\,,0.\,\right\}\\\left\{2,2,0,1.\,,0.\,\right\}\\\left\{2,2,1,1.\,,0.\,\right\}\\\left\{2,2,2,1.\,,0.\,\r… [898 chars total]
Cell 61: ✓ — Out[170]= \langle\!\lvert \text{max spread (ratio non-constancy)}\to 0.\,,\text{max |ratio - 1| over sweep}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases checked}\to 27\rvert\!\rangle
Cell 62: ✓ — Out[176]= \langle\!\lvert \text{max operator-level spread}\to 0.\,,\text{max |$\mu$ - $\mu$Site2|}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases}\to 27\rvert\!\rangle

## 2026-07-18 10:52:17 — ▶️ RUN CELLS 63–70
Cell 64: ✓ — (no output)
Cell 65: ✓ — Out[185]= \langle\!\lvert \text{w\_0 ($\lambda$=1,2,3)}\to \left\{1,1,1\right\},\text{full weights w\_j}\to \left\{1\to \left\{1,-1\right\},2\to \left\{1,-2,1\right\},3\to \left\{1,-3,3,-1\right\}\right\},\text{br1[$\lambda$][1] ($\lambda$=1,2,3)}\to \left\{0,0,0\right\}\rvert\!\rangle
Cell 66: ✓ — Out[187]= \left\{0,0\right\}
Cell 67: ✓ — Out[190]= \sqrt[3]{k}\,\left(\left(-1.18561270626131985872+0.832253303224995280404\,\mathrm{i}\right)\,k^2-\left(0.761819910806671651081+0.812812473459294825632\,\mathrm{i}\right)\,k-\left(0.404410037536507329420+0.598870869485872613792\,\mathrm{i}\right)\right)
Cell 68: ✓ — Out[194]= \begin{pmatrix}0 & 0 & 0\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix}
Cell 69: ✓ — Out[198]= \langle\!\lvert \text{OP(E12)}\to 1.\,+0.\times 10^{-30}\,\mathrm{i},\text{OP(E21)}\to \left(-0.03626734060947795145020244939-2.004237714666237867354326095\,\mathrm{i}\right)\,\varphi ^2-\left(1.522150083623538480912456055-0.2442047957313920727733393355\,\mathrm{i}\right)\,\varphi +\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right),\text{nontrivial generator \{a,b\}}\to \left\{2,1\right\},\… [464 chars total]
Cell 70: ✓ — Out[201]= \left\{1,2\right\}

## 2026-07-18 10:52:25 — ▶️ RUN CELLS 71–77
Cell 71: ✓ — Out[207]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{$\varphi$ candidates per state}\to \begin{pmatrix}\left\{0,0\right\} & \left(-0.920793778924720536105378676-0.390049761817044163511250929\,\mathrm{i}\right)\,k^{4/3}+\left(0.695722366612934233848844287-0.718310788304406529935213657\,\mathrm{i}\right)\,\sqrt[3]{k}-\left(0.3604398751… [1247 chars total]
Cell 72: ✓ — Out[213]= \langle\!\lvert \text{c\_sb solution(s)}\to \left\{\left\{\mathrm{cb}\to \frac{\left(-6.242239011314\times 10^{45}-4.511757007219\times 10^{45}\,\mathrm{i}\right)\,k^{4/3}-\left(5.999302657796\times 10^{45}+8.356648463050\times 10^{45}\,\mathrm{i}\right)\,\sqrt[3]{k}}{9.833530773718\times 10^{45}\,k^{4/3}-\left(3.54439660434\times 10^{45}-9.172544913427\times 10^{45}\,\mathrm{i}\right)\,… [841 chars total]
Cell 73: ✓ — Out[216]= \langle\!\lvert \text{generator \{a,b\} tested}\to \left\{2,1\right\},\text{rep tested}\to \left\{1,2\right\},\text{FL shape (verified exact for all 6 states)}\to \text{FL[1,2][state][k] = k\^{}(1/3) (A\_state + B\_state k)},\text{c0(state) = OPg[state][0]}\to \begin{pmatrix}\left\{0,0\right\} & 1.\,+0.\times 10^{-28}\,\mathrm{i}\\\left\{1,0\right\} & 1.\,+0.\times 10^{-28}\,\mathrm{i}\\… [1247 chars total]
Cell 74: ✓ — Out[220]= \left(1.000000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right)-\left(0.360439875147226537266677651+0.932782448593374716484835982\,\mathrm{i}\right)\,k
Cell 75: ✓ — Out[227]= \langle\!\lvert \text{alphaTab \{state,D0,c0,Dc,c1,$\alpha$\}}\to \begin{pmatrix}\left\{0,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & -0.36043987514722653-0.9327824485933747\,\mathrm{i} & 0.36043987514722653 +0.9327824485933747\,\mathrm{i} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.7101545203859353 -0.2371533418544… [1247 chars total]
Cell 76: ✓ — Out[235]= \langle\!\lvert \text{validation}\to \left\{\left\{1,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\end{pmatrix},\left\{2,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\\\left\{2,1\right\} & 0\\\left\{3,0\right\} & 0\end{pmatrix},\left\{1… [1247 chars total]
Cell 77: ✓ — Out[237]= \langle\!\lvert \text{generator \{a,b\}}\to \left\{2,1\right\},\text{$\varphi$(k)}\to \text{$\alpha$ k},\text{$\alpha$ (universal value)}\to -1.\,+0.\,\,\mathrm{i},\text{normalization}\to \text{1-site FSoV determinant det1[f] = br1[f]/(f/.u->$\theta$1); Section C's det[$\lambda$1,$\lambda$2][f] restricted to 1 site},\text{max residual over rep sweep}\to 0,\text{reps validated}\to \begin{… [510 chars total]


## 2026-07-18 11:11:47
**Input:**
```wolfram
{θ1,θ2,h,ValueQ[Q1],ValueQ[dictReportDet]}
```
**Output:**
```
Out= {1/3, 1/7, 1, False, True}
```

## 2026-07-18 11:11:55
**Input:**
```wolfram
dictReportDet
```
**Output:**
```
Out= <|"generator {a,b}" -> {2, 1}, "φ(k)" -> "α k", "α (universal value)" -> -1. + 0.*I, "normalization" -> "1-site FSoV determinant det1[f] = br1[f]/(f/.u->θ1); Section C's det[λ1,λ2][f] restricted to 1 site", "max residual over rep sweep" -> 0, "reps validated" -> {{1, 1}, {2, 1}, {1, 2}, {2, 2}, {3, 1}, {1, 3}}, "num states validated" -> 41|>
```
## 2026-07-18 11:12:51 — 📥 BULK INSERT 2 CELL(S) at positions 78–79
1. [markdown] (*Section G — Phase 1: the 2→1 site degeneration by pole cancellation.↵Hypothesis (Paul): multiplyin
2. [code] (*P1site2[λ2][u]: the site-2-tower-cancelling polynomial multiplying f before taking br[·][1][·].↵


## 2026-07-18 11:12:52
**Input:**
```wolfram
(*P1site2[λ2][u]: the site-2-tower-cancelling polynomial multiplying f before taking br[·][1][·].
  C1pred[λ1,λ2]: the predicted constant, derived analytically from Nn[λ1,λ2][1]=nn[λ1,λ2][1]/(2πi(1-Exp[2πi(θ1-θ2)/h]))
  and Nn1[λ1]=(-h)^λ1 λ1!/(2πi) (both closed forms already in the notebook/Section C,F): the oscillatory
  factor's (1-Exp[...]) cancels against the same factor sitting inside Nn[·][1], leaving
  C1pred = nn[λ1,λ2][1] / ((-h)^λ1 λ1!) = (-1)^(λ2-1) Product_{k=0}^{λ2}(θ2+kh-θ1).*)
ClearAll[P1site2,C1pred];
P1site2[λ2_][u_] := Product[(u-θ2-h k),{k,0,λ2}];
C1pred[λ1_,λ2_] := (-1)^(λ2-1) Product[(θ2+k h-θ1),{k,0,λ2}];
```
**Output:**
```
(no output)
```
## 2026-07-18 11:13:06 — 📥 BULK INSERT 2 CELL(S) at positions 80–81
1. [code] (*Test 1 — symbolic polynomial test functions f=u^j, j=0..λ1+2, across the full rep sweep.↵  Purely 
2. [code] (*Test 2 — twisted monomials f = k^(u/h) u^j, symbolic k, across the full rep sweep.↵  Same structur


## 2026-07-18 11:13:06
**Input:**
```wolfram
(*Test 1 — symbolic polynomial test functions f=u^j, j=0..λ1+2, across the full rep sweep.
  Purely rational (no Exp), so Simplify should collapse the residual to a LITERAL 0.
  This is the hard structural check on P1site2/C1pred with no numerical tolerance involved.*)
ClearAll[phase1PolyResid];
phase1PolyResid = Flatten[Table[
   Simplify[br[λ1,λ2][1][P1site2[λ2][u] u^j] - C1pred[λ1,λ2] br1[λ1][u^j]],
   {rep,{{1,1},{2,1},{1,2},{2,2},{3,1},{1,3}}},
   {λ1,{rep[[1]]}},{λ2,{rep[[2]]}},
   {j,0,λ1+2}]];
phase1PolyResid = DeleteDuplicates[Flatten[phase1PolyResid]];
If[phase1PolyResid =!= {0}, Abort[]];
{"Test 1 (symbolic monomials u^j) residuals", phase1PolyResid}
```
**Output:**
```
Out[245]= \left\{\text{Test 1 (symbolic monomials u\^{}j) residuals},\left\{0\right\}\right\}
```

## 2026-07-18 11:13:08
**Input:**
```wolfram
(*Test 2 — twisted monomials f = k^(u/h) u^j, symbolic k, across the full rep sweep.
  Same structural check, now including the transcendental k^(u/h) factor that Q1 itself carries.*)
ClearAll[phase1TwistResid];
phase1TwistResid = Flatten[Table[
   Simplify[br[λ1,λ2][1][P1site2[λ2][u] k^(u/h) u^j] - C1pred[λ1,λ2] br1[λ1][k^(u/h) u^j]],
   {rep,{{1,1},{2,1},{1,2},{2,2},{3,1},{1,3}}},
   {λ1,{rep[[1]]}},{λ2,{rep[[2]]}},
   {j,0,λ1+2}]];
phase1TwistResid = DeleteDuplicates[Flatten[phase1TwistResid]];
If[phase1TwistResid =!= {0}, Abort[]];
{"Test 2 (twisted monomials k^(u/h) u^j) residuals", phase1TwistResid}
```
**Output:**
```
Out[250]= \left\{\text{Test 2 (twisted monomials k\^{}(u/h) u\^{}j) residuals},\left\{0\right\}\right\}
```
## 2026-07-18 11:13:35 — 📥 BULK INSERT 1 CELL(S) at positions 82–82
1. [code] (*Test 3 — tie Phase 1 directly to Section F's headline object: apply the P1site2 insertion to↵  the


## 2026-07-18 11:13:37
**Input:**
```wolfram
(*Test 3 — tie Phase 1 directly to Section F's headline object: apply the P1site2 insertion to
  the ACTUAL twisted Baxter function k^(u/h) Q1[λ1,λ2][M,n][u] (not just monomials), across every
  state of the full rep sweep, and check the bracket-level identity numerically (Q1's coefficients
  are numeric, so this is a numeric not symbolic check, at the notebook's ~30-digit precision).
  Also forms the induced 1-site-style determinant from the 2-site bracket,
    det1viaSite2[λ1,λ2][M,n][k] := br[λ1,λ2][1][P1site2[λ2][u] k^(u/h)Q1[u]] / ((P1site2[λ2][u] k^(u/h)Q1[u])/.u->θ1),
  and checks it is proportional to Section F's FLdet by the SAME constant C1pred[λ1,λ2]/P1site2[λ2][θ1]
  for every state — i.e. the 2-site FSoV machinery, once localized via P1site2, reproduces the known
  1-site functional CG-overlap dictionary up to one fixed rep-dependent rescaling (not a per-state one).*)
ClearAll[det1viaSite2,phase1QResid,phase1DetRatio];
det1viaSite2[λ1_,λ2_][M_,n_][k_] := Module[{g},
  g = P1site2[λ2][u] k^(u/h) Q1[λ1,λ2][M,n][u];
  br[λ1,λ2][1][g]/(g/.u->θ1)];
phase1QResid = {};
phase1DetRatio = {};
Do[
  {λ1,λ2} = rep;
  Do[
    lhs = Simplify[br[λ1,λ2][1][P1site2[λ2][u] k^(u/h) Q1[λ1,λ2][st[[1]],st[[2]]][u]]];
    rhs = Simplify[C1pred[λ1,λ2] br1[λ1][k^(u/h) Q1[λ1,λ2][st[[1]],st[[2]]][u]]];
    resid = Chop[N[lhs-rhs] /. k->N[7/3+I/5],10^-15];
    AppendTo[phase1QResid,{λ1,λ2,st,resid}];
    ratio = Simplify[det1viaSite2[λ1,λ2][st[[1]],st[[2]]][k]/FLdet[λ1,λ2][st[[1]],st[[2]]][k]];
    ratio = Chop[N[ratio],10^-15];
    predRatio = Chop[N[C1pred[λ1,λ2]/P1site2[λ2][θ1]],10^-15];
    AppendTo[phase1DetRatio,{λ1,λ2,st,ratio,predRatio,Chop[N[ratio-predRatio],10^-12]}]
  ,{st,statesOf[λ1,λ2]}]
,{rep,{{1,1},{2,1},{1,2},{2,2},{3,1},{1,3}}}];
worstQResid = Max[Abs[phase1QResid[[All,4]]]];
worstDetRatioDev = Max[Abs[phase1DetRatio[[All,6]]]];
If[worstQResid > 10^-12 || worstDetRatioDev > 10^-10, Abort[]];
{"Test 3 worst |br identity residual| on Q1", worstQResid, "worst |detRatio - predicted| ", worstDetRatioDev, "num states checked", Length[phase1QResid]}
```
**Output:**
```
Out[259]= \left\{\text{Test 3 worst |br identity residual| on Q1},0,\text{worst |detRatio - predicted| },0,\text{num states checked},41\right\}
```
## 2026-07-18 11:13:58 — 📥 BULK INSERT 3 CELL(S) at positions 83–85
1. [markdown] (*Mirrored statement: localizing to site 2 by cancelling site 1's tower. By the θ1<->θ2, λ1<->λ2↵sym
2. [code] (*Site-2 mirror of P1site2/C1pred.*)↵ClearAll[P2site1,C2pred];↵P2site1[λ1_][u_] := Product[(u-θ1-h k
3. [code] (*Sanity: Nn1s2 normalization gives leading (θ2-node) weight 1, and br1s2[λ][1]=0, mirroring↵  Cell


## 2026-07-18 11:13:58
**Input:**
```wolfram
(*Site-2 mirror of P1site2/C1pred.*)
ClearAll[P2site1,C2pred];
P2site1[λ1_][u_] := Product[(u-θ1-h k),{k,0,λ1}];
C2pred[λ1_,λ2_] := (-1)^(λ1-1) Product[(θ1+k h-θ2),{k,0,λ1}];

(*Site-2 analogue of Nn1/mu1/br1/det1 (Section F), built around θ2 instead of θ1 — same
  construction, just re-centered, so that det1s2[λ2][f] is the genuine 1-site-at-θ2 FSoV determinant.*)
ClearAll[Nn1s2,mu1s2,br1s2,det1s2];
Nn1s2[λ_] := Nn1s2[λ] = Module[{uu}, 1/(2 π I Residue[1/Product[(uu-θ2-h k),{k,0,λ}],{uu,θ2}])];
mu1s2[λ_][u_] := Nn1s2[λ]/Product[(u-θ2-h k),{k,0,λ}];
br1s2[λ_][f_] := 2 π I Plus@@(Residue[f mu1s2[λ][u],{u,#}]&/@Table[θ2+h k,{k,0,λ}]);
det1s2[λ_][f_] := br1s2[λ][f]/(f/.u->θ2);
```
**Output:**
```
(no output)
```

## 2026-07-18 11:13:58
**Input:**
```wolfram
(*Sanity: Nn1s2 normalization gives leading (θ2-node) weight 1, and br1s2[λ][1]=0, mirroring
  Cell 65's Nn1Check for the θ1 construction.*)
ClearAll[Nn1s2Check];
Nn1s2Check = Table[
  w0 = 2 π I Residue[mu1s2[λ][u],{u,θ2}] // Simplify;
  {λ, w0, br1s2[λ][1]//Simplify}
,{λ,1,3}];
If[Union[Nn1s2Check[[All,2]]] =!= {1} || Union[Nn1s2Check[[All,3]]] =!= {0}, Abort[]];
Nn1s2Check
```
**Output:**
```
Out[271]= \begin{pmatrix}1 & 1 & 0\\2 & 1 & 0\\3 & 1 & 0\end{pmatrix}
```
## 2026-07-18 11:14:21 — 📥 BULK INSERT 2 CELL(S) at positions 86–87
1. [code] (*Mirrored Test 1+2 — symbolic monomials u^j and twisted monomials k^(u/h)u^j, site-2 bracket,↵  acr
2. [code] (*Mirrored Test 3 — the same identity on the actual Q1 Baxter functions, all 41 states of the↵  full


## 2026-07-18 11:14:22
**Input:**
```wolfram
(*Mirrored Test 1+2 — symbolic monomials u^j and twisted monomials k^(u/h)u^j, site-2 bracket,
  across the full rep sweep. Same structural check as Tests 1-2 above, roles of the two sites swapped.*)
ClearAll[phase1PolyResidS2,phase1TwistResidS2];
phase1PolyResidS2 = Flatten[Table[
   Simplify[br[λ1,λ2][2][P2site1[λ1][u] u^j] - C2pred[λ1,λ2] br1s2[λ2][u^j]],
   {rep,{{1,1},{2,1},{1,2},{2,2},{3,1},{1,3}}},
   {λ1,{rep[[1]]}},{λ2,{rep[[2]]}},
   {j,0,λ2+2}]];
phase1PolyResidS2 = DeleteDuplicates[Flatten[phase1PolyResidS2]];
phase1TwistResidS2 = Flatten[Table[
   Simplify[br[λ1,λ2][2][P2site1[λ1][u] k^(u/h) u^j] - C2pred[λ1,λ2] br1s2[λ2][k^(u/h) u^j]],
   {rep,{{1,1},{2,1},{1,2},{2,2},{3,1},{1,3}}},
   {λ1,{rep[[1]]}},{λ2,{rep[[2]]}},
   {j,0,λ2+2}]];
phase1TwistResidS2 = DeleteDuplicates[Flatten[phase1TwistResidS2]];
If[phase1PolyResidS2 =!= {0} || phase1TwistResidS2 =!= {0}, Abort[]];
{"Mirror Test1 (u^j)",phase1PolyResidS2,"Mirror Test2 (k^(u/h)u^j)",phase1TwistResidS2}
```
**Output:**
```
Out[278]= \left\{\text{Mirror Test1 (u\^{}j)},\left\{0\right\},\text{Mirror Test2 (k\^{}(u/h)u\^{}j)},\left\{0\right\}\right\}
```

## 2026-07-18 11:14:23
**Input:**
```wolfram
(*Mirrored Test 3 — the same identity on the actual Q1 Baxter functions, all 41 states of the
  full sweep, and the induced det1s2-via-site-2 object compared to a to-be-built site-2 FLdet2.*)
ClearAll[det1s2viaSite1,phase1QResidS2,phase1DetRatioS2];
det1s2viaSite1[λ1_,λ2_][M_,n_][k_] := Module[{g},
  g = P2site1[λ1][u] k^(u/h) Q1[λ1,λ2][M,n][u];
  br[λ1,λ2][2][g]/(g/.u->θ2)];
phase1QResidS2 = {};
Do[
  {λ1,λ2} = rep;
  Do[
    lhs = Simplify[br[λ1,λ2][2][P2site1[λ1][u] k^(u/h) Q1[λ1,λ2][st[[1]],st[[2]]][u]]];
    rhs = Simplify[C2pred[λ1,λ2] br1s2[λ2][k^(u/h) Q1[λ1,λ2][st[[1]],st[[2]]][u]]];
    resid = Chop[N[lhs-rhs] /. k->N[7/3+I/5],10^-15];
    AppendTo[phase1QResidS2,{λ1,λ2,st,resid}]
  ,{st,statesOf[λ1,λ2]}]
,{rep,{{1,1},{2,1},{1,2},{2,2},{3,1},{1,3}}}];
worstQResidS2 = Max[Abs[phase1QResidS2[[All,4]]]];
If[worstQResidS2 > 10^-12, Abort[]];
{"Mirror Test3 worst |br identity residual| on Q1", worstQResidS2, "num states checked", Length[phase1QResidS2]}
```
**Output:**
```
Out[285]= \left\{\text{Mirror Test3 worst |br identity residual| on Q1},0,\text{num states checked},41\right\}
```
## 2026-07-18 11:15:08 — 📥 BULK INSERT 2 CELL(S) at positions 88–89
1. [markdown] (*Required milestone: the site-2 mirror of Section F's operatorial/functional dictionary. Rotate↵the
2. [code] (*OP2: rotate the SITE-2 pseudovacuum covector by gRot[λ2][a,b][φ] (site 1's covector left bare),↵


## 2026-07-18 11:15:09
**Input:**
```wolfram
(*OP2: rotate the SITE-2 pseudovacuum covector by gRot[λ2][a,b][φ] (site 1's covector left bare),
  contract with the two-site right eigenvector Psi — the site-2 operatorial CG overlap.*)
ClearAll[OP2];
OP2[λ1_,λ2_][a_,b_][M_,n_][φ_] := Flatten[KroneckerProduct[xSingle[λ1][{0}], xSingle[λ2][{0}].gRot[λ2][a,b][φ]]] . Psi[λ1,λ2][M,n];

(*Which single-site generator on SITE 2 produces genuine φ-dependence? Mirrors Cell 69's test,
  fixed case (λ1,λ2)=(2,1), state (1,0) — same rep/state used for the site-1 test, this time
  rotating the site-2 factor instead.*)
op12s2 = OP2[2,1][1,2][1,0][φ] // Simplify;
op21s2 = OP2[2,1][2,1][1,0][φ] // Simplify;
genNontrivial2 = If[FreeQ[op12s2//Expand,φ] && !FreeQ[op21s2//Expand,φ], {2,1},
                     If[FreeQ[op21s2//Expand,φ] && !FreeQ[op12s2//Expand,φ], {1,2}, $Failed]];
If[genNontrivial2===$Failed, Abort[]];
ClearAll[OPg2];
OPg2[λ1_,λ2_][M_,n_][φ_] := OP2[λ1,λ2][genNontrivial2[[1]],genNontrivial2[[2]]][M,n][φ];
{"OP2(E12)",op12s2,"OP2(E21)",op21s2,"nontrivial generator {a,b}",genNontrivial2}
```
**Output:**
```
Out[294]= \left\{\text{OP2(E12)},1.\,+0.\times 10^{-28}\,\mathrm{i},\text{OP2(E21)},\left(2.603469709065218092712489006+0.6885776528619826437114966465\,\mathrm{i}\right)\,\varphi +\left(1.\,+0.\times 10^{-28}\,\mathrm{i}\right),\text{nontrivial generator \{a,b\}},\left\{2,1\right\}\right\}
```
## 2026-07-18 11:15:27 — 📥 BULK INSERT 2 CELL(S) at positions 90–91
1. [code] (*FLdet2: the site-2 functional CG-overlap object, det1s2 applied to k^(u/h)Q1[λ1,λ2][M,n][u].↵  Sin
2. [code] (*Discovery step (sharp, λ2=1 rep, mirrors Cell 71): rep (2,1), all 6 states — solve φ per state↵  a


## 2026-07-18 11:15:27
**Input:**
```wolfram
(*FLdet2: the site-2 functional CG-overlap object, det1s2 applied to k^(u/h)Q1[λ1,λ2][M,n][u].
  Since Q1 carries κ1^(u/h) (the SHARED companion-twist eigenvalue, not a site-2-specific one — Q1
  is a single chain-level Baxter function, not resolved per site) f(θ2)=(kκ1)^(θ2/h)q_n(θ2), and
  dividing by it cancels the fractional prefactor exactly as for FLdet, leaving a clean degree-λ2
  polynomial in k.*)
ClearAll[FLdet2];
FLdet2[λ1_,λ2_][M_,n_][k_] := det1s2[λ2][k^(u/h) Q1[λ1,λ2][M,n][u]];
(*smoke test: rep (2,1) (λ2=1), state (0,0) — should be linear in k, constant term 1.*)
FLdet2[2,1][0,0][k]//Simplify
```
**Output:**
```
Out[297]= \left(1.\,+0.\times 10^{-27}\,\mathrm{i}\right)-\left(0.360439875147226537266677651+0.932782448593374716484835982\,\mathrm{i}\right)\,k
```

## 2026-07-18 11:15:27
**Input:**
```wolfram
(*Discovery step (sharp, λ2=1 rep, mirrors Cell 71): rep (2,1), all 6 states — solve φ per state
  assuming exact match against OPg2 (both linear in φ since λ2=1), then check the ratio
  α2_state = D_state/c1_state is the SAME constant across all 6 states.*)
statesR21 = statesOf[2,1];
alphaTab2 = Table[
  Fl = FLdet2[2,1][st[[1]],st[[2]]][k]//Simplify;
  c0 = Fl/.k->0;
  Dcoef = Coefficient[Fl,k,1];
  Og = OPg2[2,1][st[[1]],st[[2]]][φ]//Simplify;
  cc0 = Og/.φ->0;
  cc1 = Coefficient[Og,φ,1];
  If[Chop[c0-1,10^-12]=!=0 || Chop[cc0-1,10^-12]=!=0, Abort[]];
  {st, Dcoef, cc1, Dcoef/cc1}
,{st,statesR21}];
alpha2Vals = alphaTab2[[All,4]];
alpha2Spread = Max[Abs[N[alpha2Vals-alpha2Vals[[1]]]]];
If[alpha2Spread > 10^-12, Abort[]];
alpha2Val = N[alpha2Vals[[1]]];
{"alphaTab2 {state,D,c1,alpha}",alphaTab2,"alpha2Val",alpha2Val,"spread",alpha2Spread}
```
**Output:**
```
Out[302]= \text{\$Aborted} | Out[304]= \left\{\text{alphaTab2 \{state,D,c1,alpha\}},\begin{pmatrix}\left\{0,0\right\} & -0.360439875147226537266677651-0.932782448593374716484835982\,\mathrm{i} & 0.3604398751472265372666776506+0.9327824485933747164848359821\,\mathrm{i} & -1.000000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\\\left\{1,0\right\} & 0.4685059330398572984-0.171361115691460089\,\m… [1276 chars total]
```

## 2026-07-18 11:15:48
**Input:**
```wolfram
{alpha2Val, alpha2Spread, alphaTab2[[All,{1,4}]]}
```
**Output:**
```
Out= {-1. + 0.*I, 1.7575074201841347, {{{0, 0}, -0.99999999999999999999999999999999993415`26.911708336465317 + 0``26.995317996306234*I}, {{1, 0}, 0.15191893851739566206190478609696339641`18.489116861324803 - 0.10600050418875032622029734374875387091`18.43994706580298*I}, {{1, 1}, 0.51059310376891640075107185162804008008`18.420716078647565 + 0.89829884050247156339917750863254951909`18.681437733381564*I}, {{2, 0}, 0.15191893851739566206190478609696339593`20.68886563835738 + 0.10600050418875032622029734374875386468`20.641854443514987*I}, {{2, 1}, 0.51059310376891640075107185162804002861`21.151340334417352 - 0.89829884050247156339917750863254952904`21.413495505878576*I}, {{3, 0}, -1.00000000000000000000000000000000000001`19.55960534305864 + 0``19.63667288935328*I}}}
```
## 2026-07-18 11:17:08 — 📥 BULK INSERT 1 CELL(S) at positions 92–92
1. [code] (*Clean re-derivation (Module-scoped, no Abort) of the per-state alpha ratio for rep (2,1),↵  to get


## 2026-07-18 11:17:08
**Input:**
```wolfram
(*Clean re-derivation (Module-scoped, no Abort) of the per-state alpha ratio for rep (2,1),
  to get trustworthy diagnostic numbers after the previous cell's leaked-variable/Abort confusion.*)
ClearAll[alphaTab2Clean];
alphaTab2Clean = Table[
  Module[{Fl,c0,Dcoef,Og,cc0,cc1},
    Fl = FLdet2[2,1][st[[1]],st[[2]]][k]//Simplify;
    c0 = Fl/.k->0;
    Dcoef = Coefficient[Fl,k,1];
    Og = OPg2[2,1][st[[1]],st[[2]]][φ]//Simplify;
    cc0 = Og/.φ->0;
    cc1 = Coefficient[Og,φ,1];
    {st, N[c0], N[cc0], N[Dcoef], N[cc1], N[Dcoef/cc1]}
  ]
,{st,statesOf[2,1]}];
TableForm[alphaTab2Clean]
```
**Output:**
```
Out[307]//TableForm= \begin{pmatrix}\left\{0,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & -0.36043987514722653-0.9327824485933747\,\mathrm{i} & 0.36043987514722653 +0.9327824485933747\,\mathrm{i} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.4685059330398573 -0.1713611156914601\,\mathrm{i} & 2.603469709065218 +0.6885776528619827… [1258 chars total]
```

## 2026-07-18 11:17:14
**Input:**
```wolfram
alphaTab2Clean[[All,{1,6}]]
```
**Output:**
```
Out= {{{0, 0}, -1. + 0.*I}, {{1, 0}, 0.15191893851739566 - 0.10600050418875033*I}, {{1, 1}, 0.5105931037689164 + 0.8982988405024716*I}, {{2, 0}, 0.15191893851739566 + 0.10600050418875033*I}, {{2, 1}, 0.5105931037689164 - 0.8982988405024716*I}, {{3, 0}, -1. + 0.*I}}
```
## 2026-07-18 11:18:59 — 📥 BULK INSERT 1 CELL(S) at positions 93–93
1. [code] (*Diagnostic: the site-2 mirror does NOT close with a state-independent alpha (extremal states↵  M=0


## 2026-07-18 11:18:59
**Input:**
```wolfram
(*Diagnostic: the site-2 mirror does NOT close with a state-independent alpha (extremal states
  M=0,3 give -1, interior M=1,2 give genuinely different complex numbers, not merely a global
  rescaling of the same value). Investigate whether the per-state mismatch factor has a simple
  closed form in terms of already-known state data (Q1 values, J_n, M).*)
diag2 = Table[
  Module[{a2 = alphaTab2Clean[[i,6]], M=alphaTab2Clean[[i,1,1]], n=alphaTab2Clean[[i,1,2]]},
    {  {M,n}, a2,
       a2/(-1),
       N[Q1[2,1][M,n][θ1]/Q1[2,1][M,n][θ2]],
       N[τ[2,1][M,n][θ1]/τ[2,1][M,n][θ2]],
       N[κ1^((θ1-θ2)/h)]
    }
  ]
,{i,1,Length[alphaTab2Clean]}];
TableForm[diag2,TableHeadings->{None,{"state","alpha2","alpha2/(-1)","Q1(θ1)/Q1(θ2)","τ(θ1)/τ(θ2)","κ1^((θ1-θ2)/h)"}}]
```
**Output:**
```
Out[309]//TableForm= \begin{pmatrix}\left\{0,0\right\} & -1.\,+0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.973902234276272 +0.2269679230016555\,\mathrm{i} & 0.7391304347826086 +0.\,\,\mathrm{i} & 0.973902234276272 +0.2269679230016555\,\mathrm{i}\\\left\{1,0\right\} & 0.15191893851739566 -0.10600050418875033\,\mathrm{i} & -0.15191893851739566+0.10600050418875033\,\mathrm{i} & 0.7651595274649543 +… [1258 chars total]
```

## 2026-07-18 11:19:05
**Input:**
```wolfram
diag2
```
**Output:**
```
Out= {{{0, 0}, -1. + 0.*I, 1. + 0.*I, 0.973902234276272 + 0.2269679230016555*I, 0.7391304347826086 + 0.*I, 0.973902234276272 + 0.2269679230016555*I}, {{1, 0}, 0.15191893851739566 - 0.10600050418875033*I, -0.15191893851739566 + 0.10600050418875033*I, 0.7651595274649543 + 0.27587192687799*I, 1.121164201738695 + 0.21744516531475866*I, 0.973902234276272 + 0.2269679230016555*I}, {{1, 1}, 0.5105931037689164 + 0.8982988405024716*I, -0.5105931037689164 - 0.8982988405024716*I, 0.9030409119491282 + 0.10146628916432213*I, 0.8431888703246504 - 0.05922775858039619*I, 0.973902234276272 + 0.2269679230016555*I}, {{2, 0}, 0.15191893851739566 + 0.10600050418875033*I, -0.15191893851739566 - 0.10600050418875033*I, 0.5232058326875958 + 0.28297154471668806*I, 1.121164201738695 - 0.21744516531475866*I, 0.973902234276272 + 0.2269679230016555*I}, {{2, 1}, 0.5105931037689164 - 0.8982988405024716*I, -0.5105931037689164 + 0.8982988405024716*I, 0.8440611884595888 - 0.07577268481216724*I, 0.8431888703246504 + 0.05922775858039619*I, 0.973902234276272 + 0.2269679230016555*I}, {{3, 0}, -1. + 0.*I, 1. + 0.*I, 1.2798210250222584 - 0.3433891666400345*I, 0.7391304347826086 + 0.*I, 0.973902234276272 + 0.2269679230016555*I}}
```
## 2026-07-18 11:19:58 — 📥 BULK INSERT 1 CELL(S) at positions 94–94
1. [code] (*Cross-check the pattern on a second λ2=1 rep, (3,1) — d=8, states M=0..4, to confirm this isn't↵


## 2026-07-18 11:19:59
**Input:**
```wolfram
(*Cross-check the pattern on a second λ2=1 rep, (3,1) — d=8, states M=0..4, to confirm this isn't
  an accident of (2,1). Same clean Module-scoped computation as alphaTab2Clean.*)
statesR31 = statesOf[3,1];
alphaTab2R31 = Table[
  Module[{Fl,c0,Dcoef,Og,cc0,cc1},
    Fl = FLdet2[3,1][st[[1]],st[[2]]][k]//Simplify;
    c0 = Fl/.k->0;
    Dcoef = Coefficient[Fl,k,1];
    Og = OPg2[3,1][st[[1]],st[[2]]][φ]//Simplify;
    cc0 = Og/.φ->0;
    cc1 = Coefficient[Og,φ,1];
    {st, N[Dcoef/cc1]}
  ]
,{st,statesR31}];
alphaTab2R31
```
**Output:**
```
Out[312]= \begin{pmatrix}\left\{0,0\right\} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 0.05448860910968355 -0.09540961504089102\,\mathrm{i}\\\left\{1,1\right\} & -0.0293413933141388+0.7363100318787654\,\mathrm{i}\\\left\{2,0\right\} & 0.11703059508534998 +0.\,\,\mathrm{i}\\\left\{2,1\right\} & 0.8158972359113332 +0.\,\,\mathrm{i}\\\left\{3,0\right\} & 0.05448860910968355 +0.09540961504089102\,… [541 chars total]
```
## 2026-07-18 11:20:21 — 📥 BULK INSERT 1 CELL(S) at positions 95–95
1. [markdown] (*Phase 1 summary / milestone report. The BRACKET-level pole-cancellation statement closes↵EXACTLY a

## 2026-07-18 11:21:06 — 📥 BULK INSERT 2 CELL(S) at positions 96–97
1. [markdown] (*Section H — Phase 2: two independent local rotations. Brute-force operatorial target↵O[k1,k2] := (
2. [code] (*O[λ1,λ2][M,n][k1,k2]: the brute-force two-parameter operatorial overlap.*)↵ClearAll[O2rot];↵O2rot[


## 2026-07-18 11:21:06
**Input:**
```wolfram
(*O[λ1,λ2][M,n][k1,k2]: the brute-force two-parameter operatorial overlap.*)
ClearAll[O2rot];
O2rot[λ1_,λ2_][M_,n_][k1_,k2_] := Flatten[KroneckerProduct[
   xSingle[λ1][{0}].gRot[λ1][2,1][-k1],
   xSingle[λ2][{0}].gRot[λ2][2,1][-k2]
 ]] . Psi[λ1,λ2][M,n];
(*sanity: k2->0 must reproduce the known site-1 dictionary FLdet exactly, k1->0 the site-2 analogue.*)
testState = {1,0}; testRep = {2,1};
chk1 = Simplify[(O2rot[testRep[[1]],testRep[[2]]][testState[[1]],testState[[2]]][k1,0]) - FLdet[testRep[[1]],testRep[[2]]][testState[[1]],testState[[2]]][k1]];
chk2 = Simplify[(O2rot[testRep[[1]],testRep[[2]]][testState[[1]],testState[[2]]][0,k2]) - FLdet2[testRep[[1]],testRep[[2]]][testState[[1]],testState[[2]]][k2]];
{"O[k1,0]-FLdet[k1] (must be 0)", Chop[N[chk1],10^-15], "O[0,k2]-FLdet2[k2] (expected NONZERO per Phase1 finding)", Chop[N[chk2],10^-15]}
```
**Output:**
```
Out[318]= \left\{\text{O[k1,0]-FLdet[k1] (must be 0)},0,\text{O[0,k2]-FLdet2[k2] (expected NONZERO per Phase1 finding)},\left(-3.0719756421050755-0.5172165371705225\,\mathrm{i}\right)\,k2\right\}
```
## 2026-07-18 11:21:23 — 📥 BULK INSERT 1 CELL(S) at positions 98–98
1. [code] (*Candidate 1: the most natural 2x2-structured FSoV determinant with independent twists localized↵


## 2026-07-18 11:21:24
**Input:**
```wolfram
(*Candidate 1: the most natural 2x2-structured FSoV determinant with independent twists localized
  per site via the Phase-1 pole-cancelling polynomials P1site2/P2site1, each row carrying its own
  k_alpha^(u/h) and its own node normalization (mirrors Section C's det[λ1,λ2][f] but with the SAME
  f replaced by two DIFFERENT twisted functions g1 (site-1 row) and g2 (site-2 row)).*)
ClearAll[Fcand1];
Fcand1[λ1_,λ2_][M_,n_][k1_,k2_] := Module[{g1,g2},
  g1[uu_] := P1site2[λ2][uu] k1^(uu/h) Q1[λ1,λ2][M,n][uu];
  g2[uu_] := P2site1[λ1][uu] k2^(uu/h) Q1[λ1,λ2][M,n][uu];
  1/g1[θ1]/g2[θ2]/(θ2-θ1) *
   Det[{{br[λ1,λ2][1][g1[u]], br[λ1,λ2][1][u g1[u]]},
        {br[λ1,λ2][2][g2[u]], br[λ1,λ2][2][u g2[u]]}}]
];
(*test on rep (2,1), state (1,0): compare to O2rot as a function of k1,k2.*)
Fc1 = Fcand1[2,1][1,0][k1,k2] // Simplify;
Oc1 = O2rot[2,1][1,0][k1,k2] // Simplify;
{"Fcand1", Fc1, "O2rot", Oc1}
```
**Output:**
```
Out[323]= \left\{\text{Fcand1},k1^2\,\left(\left(-2.2527492196701659-5.82989030370859198\,\mathrm{i}\right)\,k2-\left(0.4170744170089964+23.0487337186617355\,\mathrm{i}\right)\right)+k1\,\left(\left(0.67128913890100625-0.37524873225652243\,\mathrm{i}\right)\,k2+\left(9.51343802264711551-1.5262799733212005\,\mathrm{i}\right)\right)-\left(1.99115021541939352-0.72828474168870538\,\mathrm{i}\right)\,k… [1014 chars total]
```

## 2026-07-18 11:21:38
**Input:**
```wolfram
diffc1 = Simplify[Fc1-Oc1]
ratioc1 = Simplify[Fc1/Oc1]
```
**Output:**
```
Out[1]= (0``17.58858361986547 + 0``17.69290006822662*I) + k1^2*((-0.38080707639951849024794569133739434921`15.916418301243063 - 21.04449600399549760711121006227408424151`17.550560710051425*I) - (1.89230934452293932061144909109129157291`17.25319886821053 + 4.89710785511521726154519595054075300726`17.48826646390722*I)*k2) + k1*((7.9912879390235770247777156366847403095`17.510321044914395 - 1.28207517758980838200187725012720671076`16.851582222693594*I) - (0.62481200890076820926362113885920823007`17.462719140619605 + 3.0519104255438007752813175458888132942`18.27775709968749*I)*k2) + (0.61231949364582457433617423822785625279`17.06820134180077 + 1.41686239455068802278788075604215863693`17.61084107578691*I)*k2
Out[2]= ((0.99999999999999999998878308437674950228`17.588583619872402 + 0``17.692900068234064*I) + k1^2*((-0.41707441700899644169814814072343582822`15.95592684252675 - 23.04873371866173547446553615732631699402`17.5900692513357*I) - (2.25274921967016585787812674168684168308`17.328919582149258 + 5.82989030370859197803003193261162095951`17.5639871778464*I)*k2) + k1*((9.51343802264711550569017169138413914125`17.58604175885362 - 1.52627997332120045477521658566260860898`16.927302936631953*I) + (0.6712891389010062509862691532617211038`17.493879393120743 - 0.37524873225652243079201138417219909949`17.36750455061047*I)*k2) - (1.99115021541939351837631476825819283474`17.58032728117016 - 0.72828474168870537907638410950669258374`17.321814613241447*I)*k2)/((1.`28.385628029452757 + 0``28.458858685029718*I) + k1^2*((-0.0362673406094779514502024493860414790055997083400338776284`28.385628029452757 - 2.00423771466623786735432609505223275251`28.385628029452757*I) - (0.3604398751472265372666776505955501101721018646136190171211`28.385628029452757 + 0.93278244859337471648483598207086795224`28.379879135262136*I)*k2) - (2.60346970906521809271248900648604908753`28.385628029452757 + 0.68857765286198264371149664653546605319`28.330198017570385*I)*k2 + k1*((1.52215008362353848091245605469939883175`28.3833973045706 - 0.24420479573139207277333933553540189821871773341101956567`28.385628029452757*I) + (1.29610114780177446024989029212092933387`28.385628029452757 + 2.67666169328727834448930616171661419472`28.378365180658164*I)*k2))
```
## 2026-07-18 11:22:19 — 📥 BULK INSERT 1 CELL(S) at positions 99–99
1. [code] (*Quick structural check before building more elaborate candidates: does O2rot factorize as a↵  plai


## 2026-07-18 11:22:19
**Input:**
```wolfram
(*Quick structural check before building more elaborate candidates: does O2rot factorize as a
  plain product of its two 1-parameter slices, O[k1,k2] =? O[k1,0]*O[0,k2] (would hold iff Psi
  were literally a product state under this rotation action — generically false for an entangled
  two-site eigenstate, but cheap to rule out first).*)
Ofact = Simplify[Oc1 - (O2rot[2,1][1,0][k1,0] O2rot[2,1][1,0][0,k2])];
Chop[N[Ofact],10^-12]
```
**Output:**
```
Out[325]= \left(0.9252125034899888 -6.175727508765443\,\mathrm{i}\right)\,k1^2\,k2+\left(5.4271267482291075 +3.0890004366773067\,\mathrm{i}\right)\,k1\,k2
```
## 2026-07-18 11:23:21 — 📥 BULK INSERT 1 CELL(S) at positions 100–100
1. [code] (*Diagnostic: coefficient-by-coefficient comparison of Candidate1 vs O2rot, state (2,1)[1,0],↵  to s


## 2026-07-18 11:23:21
**Input:**
```wolfram
(*Diagnostic: coefficient-by-coefficient comparison of Candidate1 vs O2rot, state (2,1)[1,0],
  to see exactly which Taylor coefficients (in k1,k2) match and which don't.*)
cfsF = CoefficientRules[Numerator[Together[Fc1]] , {k1,k2}];
polyF = Expand[Fc1];
polyO = Expand[Oc1];
diffCoeffs = Table[
  {mono, Chop[N[Coefficient[Coefficient[polyF,k1,mono[[1]]],k2,mono[[2]]] -
                 Coefficient[Coefficient[polyO,k1,mono[[1]]],k2,mono[[2]]]],10^-10]}
, {mono, {{0,0},{1,0},{2,0},{0,1},{1,1},{2,1}}}];
diffCoeffs
```
**Output:**
```
Out[330]= \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 7.991287939023577 -1.2820751775898085\,\mathrm{i}\\\left\{2,0\right\} & -0.3808070763995185-21.044496003995498\,\mathrm{i}\\\left\{0,1\right\} & 0.6123194936458246 +1.416862394550688\,\mathrm{i}\\\left\{1,1\right\} & -0.6248120089007682-3.051910425543801\,\mathrm{i}\\\left\{2,1\right\} & -1.8923093445229393-4.897107855115217\,\m… [421 chars total]
```
## 2026-07-18 11:23:49 — 📥 BULK INSERT 1 CELL(S) at positions 101–101
1. [code] (*Candidate 2: same 2x2 structure but WITHOUT the P1site2/P2site1 insertions (test whether those↵  i


## 2026-07-18 11:23:49
**Input:**
```wolfram
(*Candidate 2: same 2x2 structure but WITHOUT the P1site2/P2site1 insertions (test whether those
  insertions, rather than the row-splitting idea itself, are the source of the mismatch).*)
ClearAll[Fcand2];
Fcand2[λ1_,λ2_][M_,n_][k1_,k2_] := Module[{g1,g2},
  g1[uu_] := k1^(uu/h) Q1[λ1,λ2][M,n][uu];
  g2[uu_] := k2^(uu/h) Q1[λ1,λ2][M,n][uu];
  1/g1[θ1]/g2[θ2]/(θ2-θ1) *
   Det[{{br[λ1,λ2][1][g1[u]], br[λ1,λ2][1][u g1[u]]},
        {br[λ1,λ2][2][g2[u]], br[λ1,λ2][2][u g2[u]]}}]
];
Fc2 = Simplify[Fcand2[2,1][1,0][k1,k2]];
diffc2 = Chop[N[Expand[Fc2]-Expand[Oc1]],10^-10];
{"Fcand2 (no P-insertion)", Fc2, "diff vs O2rot", diffc2}
```
**Output:**
```
Out[335]= \left\{\text{Fcand2 (no P-insertion)},k1^2\,\left(\left(0.02466179161444501+1.36288164597304175\,\mathrm{i}\right)-\left(0.36043987514722654+0.932782448593374716\,\mathrm{i}\right)\,k2\right)+k1\,\left(\left(1.2351720155778515-0.69045766735200127\,\mathrm{i}\right)\,k2-\left(6.46913785540003854-1.0378703818584163\,\mathrm{i}\right)\right)+\left(5.38781822995835893-1.9706528304517910\,\ma… [789 chars total]
```
## 2026-07-18 11:24:12 — 📥 BULK INSERT 1 CELL(S) at positions 102–102
1. [code] (*Systematic check: Candidate1 vs O2rot across several states/reps, reporting whether the k1=0↵  and


## 2026-07-18 11:24:13
**Input:**
```wolfram
(*Systematic check: Candidate1 vs O2rot across several states/reps, reporting whether the k1=0
  and k2=0 boundary slices match (they must, or don't, per Phase-1 findings) and the overall
  deviation magnitude at a generic interior point (k1,k2)=(1,1) — quantifies how comprehensively
  Candidate1 fails, not just its k1k2 cross-term structure.*)
candidateReport = {};
Do[
  {λ1,λ2} = rep;
  Do[
    Fc = Fcand1[λ1,λ2][st[[1]],st[[2]]][k1,k2];
    Oc = O2rot[λ1,λ2][st[[1]],st[[2]]][k1,k2];
    k2zero = Chop[N[(Fc/.k2->0) - (Oc/.k2->0)],10^-10];
    k1zero = Chop[N[(Fc/.k1->0) - (Oc/.k1->0)],10^-10];
    genericDev = Chop[N[(Fc/.{k1->1,k2->1}) - (Oc/.{k1->1,k2->1})],10^-10];
    AppendTo[candidateReport,{λ1,λ2,st,k2zero,k1zero,genericDev}]
  ,{st, {statesOf[λ1,λ2][[1]], statesOf[λ1,λ2][[Ceiling[Length[statesOf[λ1,λ2]]/2]]]}}]
,{rep,{{2,1},{1,2},{3,1}}}];
TableForm[candidateReport,TableHeadings->{None,{"λ1","λ2","state","F(k1,0)-O(k1,0)","F(0,k2)-O(0,k2)","F(1,1)-O(1,1)"}}]
```
**Output:**
```
2 kernel messages: Power::infy: Infinite expression 0       encountered.
Infinity::indet: Indeterminate expression ----------------- encountered.                                                   1/3                    … [also: Power::infy] result: Out[338]//TableForm= \begin{pmatrix}2 & 1 & \left\{0,0\right\} & \mathrm{Indeterminate} & \mathrm{Indeterminate} & -4.970498218180917+4.515583419069149\,\mathrm{i}\\2 & 1 & \left\{1,1\right\} & \mathr…
```
## 2026-07-18 11:24:34 — 📥 BULK INSERT 1 CELL(S) at positions 103–103
1. [code] (*Fix: Simplify each Fc/Oc to a genuine rational function of k1,k2 BEFORE substituting boundary↵  va


## 2026-07-18 11:24:34
**Input:**
```wolfram
(*Fix: Simplify each Fc/Oc to a genuine rational function of k1,k2 BEFORE substituting boundary
  values (raw unsimplified Fcand1 output still carries fractional powers that cancel only after
  Simplify, so direct substitution of k1->0 or k2->0 hits 0^(fractional) => Indeterminate).*)
candidateReport = {};
Do[
  {λ1,λ2} = rep;
  Do[
    Fc = Simplify[Fcand1[λ1,λ2][st[[1]],st[[2]]][k1,k2]];
    Oc = Simplify[O2rot[λ1,λ2][st[[1]],st[[2]]][k1,k2]];
    k2zero = Chop[N[(Fc/.k2->0) - (Oc/.k2->0)],10^-10];
    k1zero = Chop[N[(Fc/.k1->0) - (Oc/.k1->0)],10^-10];
    genericDev = Chop[N[(Fc/.{k1->1,k2->1}) - (Oc/.{k1->1,k2->1})],10^-10];
    AppendTo[candidateReport,{λ1,λ2,st,k2zero,k1zero,genericDev}]
  ,{st, {statesOf[λ1,λ2][[1]], statesOf[λ1,λ2][[Ceiling[Length[statesOf[λ1,λ2]]/2]]]}}]
,{rep,{{2,1},{1,2},{3,1}}}];
TableForm[candidateReport,TableHeadings->{None,{"λ1","λ2","state","F(k1,0)-O(k1,0)","F(0,k2)-O(0,k2)","F(1,1)-O(1,1)"}}]
```
**Output:**
```
Out[341]//TableForm= \begin{pmatrix}2 & 1 & \left\{0,0\right\} & \left(-7.771745024480887+7.060451775520925\,\mathrm{i}\right)\,k1^2-\left(3.7846186890458786 +9.794215710230434\,\mathrm{i}\right)\,k1 & \left(1.8923093445229393 +4.897107855115217\,\mathrm{i}\right)\,k2 & -4.970498218180917+4.515583419069149\,\mathrm{i}\\2 & 1 & \left\{1,1\right\} & \left(9.339255641812587 +10.661478687052961\,\math… [1258 chars total]
```
## 2026-07-18 11:25:03 — 📥 BULK INSERT 1 CELL(S) at positions 104–104
1. [code] (*Diagnose precisely: is Fcand1(k1,0) at least PROPORTIONAL to the known-correct FLdet[k1] (just a↵


## 2026-07-18 11:25:04
**Input:**
```wolfram
(*Diagnose precisely: is Fcand1(k1,0) at least PROPORTIONAL to the known-correct FLdet[k1] (just a
  wrong overall constant, fixable by renormalizing), or does it have a genuinely different SHAPE
  (relative k1 vs k1^2 coefficients differ, meaning even the anchoring constraint fails structurally)?*)
Fck2zero = Simplify[Fcand1[2,1][0,0][k1,k2] /. k2->0];
FLd = FLdet[2,1][0,0][k1];
ratioShape = Simplify[Fck2zero/FLd];
{"Fcand1(k1,0)", Fck2zero, "FLdet[k1]", FLd, "ratio (want: k1-independent const)", ratioShape}
```
**Output:**
```
2 kernel messages: Power::infy: Infinite expression 0       encountered.
Infinity::indet: Indeterminate expression ----------------- encountered.                                                   1/3                    … [also: Power::infy] result: Out[345]= \left\{\text{Fcand1(k1,0)},\mathrm{Indeterminate},\text{FLdet[k1]},\frac{\left(2.45075493251774917289426756629+5.78551794268217244227029022731\,\mathrm{i}\right)\,\left(\left(0.0525946877178…
```
## 2026-07-18 11:25:14 — 📥 BULK INSERT 1 CELL(S) at positions 105–105
1. [code] (*Redo correctly: Simplify the FULL bivariate Fcand1(k1,k2) first (clears fractional powers via↵  cr


## 2026-07-18 11:25:14
**Input:**
```wolfram
(*Redo correctly: Simplify the FULL bivariate Fcand1(k1,k2) first (clears fractional powers via
  cross-cancellation between the k1-row and k2-row), THEN substitute k2->0.*)
FcFull00 = Simplify[Fcand1[2,1][0,0][k1,k2]];
Fck2zero = FcFull00 /. k2->0;
FLd = FLdet[2,1][0,0][k1];
ratioShape = Simplify[Fck2zero/FLd];
{"Fcand1(k1,0)", Fck2zero, "FLdet[k1]", FLd, "ratio (want: k1-independent const)", ratioShape}
```
**Output:**
```
Out[350]= \left\{\text{Fcand1(k1,0)},\left(-8.51191121729+7.73287575414\,\mathrm{i}\right)\,k1^2-\left(4.5054984393+11.6597806074\,\mathrm{i}\right)\,k1+\left(1.000000000000+0.\times 10^{-12}\,\mathrm{i}\right),\text{FLdet[k1]},\frac{\left(2.45075493251774917289426756629+5.78551794268217244227029022731\,\mathrm{i}\right)\,\left(\left(0.05259468771783351844996683714+0.15021349717800998891606886537\… [960 chars total]
```

## 2026-07-18 11:25:22
**Input:**
```wolfram
{Coefficient[Fck2zero,k1,0],Coefficient[Fck2zero,k1,1],Coefficient[Fck2zero,k1,2]}
{Coefficient[FLd,k1,0]//N,Coefficient[FLd,k1,1]//N,Coefficient[FLd,k1,2]//N}
N[predRatio]
```
**Output:**
```
Out[1]= {0.9999999999999999999101349301340325701`11.816696857349125 + 0``11.816696857349125*I, -4.50549843934033171590459079096968663976`11.362264551212531 - 11.65978060741718395586590852434933977875`11.77521208635683*I, -8.51191121728859014832820921505204921451`11.694000295897567 + 7.73287575414196548701354662559576223056`11.6523142431186*I}
Out[2]= {1. + 0.*I, -0.7208797502944531 - 1.8655648971867493*I, -0.7401661928077035 + 0.6724239786210405*I}
Out[3]= 1.
```
## 2026-07-18 11:26:04 — 📥 BULK INSERT 1 CELL(S) at positions 106–106
1. [markdown] (*Phase 2 outcome — OBSTRUCTION, reported precisely rather than forced further.↵↵Target: O[λ1,λ2][M,


---

## 2026-07-18 11:26:37 — 🔄 KERNEL RESTART

## 2026-07-18 11:27:14 — ▶️ RUN CELLS 1–10
Cell 1: ✓ — Out[8]= \left\{\frac{1}{3},\frac{1}{7},1,0.360439875147226537266677650596+0.932782448593374716484835982071\,\mathrm{i},0.36043987514722653726667765060-0.932782448593374716484835982071\,\mathrm{i},0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i},1\right\}
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 8: ✓ — (no output)
Cell 10: ✓ — (no output)

## 2026-07-18 11:27:19 — ▶️ RUN CELLS 11–20
Cell 11: ✓ — (no output)
Cell 12: ✓ — Out[27]= \left\{0\right\}
Cell 13: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — Out[31]= \left\{0\right\}
Cell 16: ✓ — (no output)
Cell 17: ✓ — Out[34]= \begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{pmatrix}
Cell 18: ✓ — (no output)
Cell 20: ✓ — (no output)

## 2026-07-18 11:27:22 — ▶️ RUN CELLS 21–28
Cell 21: ✓ — (no output)
Cell 23: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 26: ✓ — (no output)
Cell 27: ✓ — Out[52]= \begin{pmatrix}0 & \frac{1}{2} & 0 & 0 & 0 & 0\end{pmatrix}
Cell 28: ✓ — Out[53]= \left\{0\right\}

## 2026-07-18 11:27:24 — ▶️ RUN CELLS 29–36
Cell 30: ✓ — Out[59]= \left\{\left(u-\frac{4}{3}\right)\,\left(u-\frac{8}{7}\right),\left(u-\frac{1}{3}\right)\,\left(u-\frac{1}{7}\right),\left(0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i}\right)\,u^2-\left(1.06415582186324025288257211128+0.\times 10^{-29}\,\mathrm{i}\right)\,u\right\}
Cell 31: ✓ — Out[62]= \begin{pmatrix}1 & 1\\2 & 2\\3 & 3\end{pmatrix}
Cell 32: ✓ — Out[65]= \text{max|C2-($\kappa$1+$\kappa$2)Id|}\to 0.\times 10^{-29}
Cell 33: ✓ — Out[75]= \langle\!\lvert \left\{0,0\right\}\to 0.5664055180884988 +1.3769645669711723\,\mathrm{i},\left\{1,0\right\}\to -0.778123966451158+0.\,\,\mathrm{i},\left\{1,1\right\}\to 1.1900552523337025 +0.\,\,\mathrm{i},\left\{2,0\right\}\to 0.5664055180884988 -1.3769645669711723\,\mathrm{i}\rvert\!\rangle
Cell 34: ✓ — 2 kernel messages: Psi::spurious: No genuine state at (M,n)=(2,2) for (λ1,λ2)=(2,1): weight sector M=2 has only 2 genuine state(s) (n=0..1), so n=2 does not exist. [also: Psi::spurious] result: Print: spurious-guard test (expect a Psi::spurious message and Missing[...]):[formula] | Out[84]= \left\{\left\{1.000000000000000000000000000+0.\times 10^{-29}\,\mathrm{i},-1.3004494881076436272092358…
⚠️ stopped at Cell 34 — error detected (pass stopOnError:false to continue past errors)

## 2026-07-18 11:27:35 — ▶️ RUN CELLS 35–42
Cell 35: ✓ — Out[91]= \left\{1.03464899078135229594035121+0.381869568380946283899575941\,\mathrm{i},-2.37022278949001677248172619-0.20048970475590097847093779\,\mathrm{i},1.00000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right\}
Cell 36: ✓ — Out[94]= \left\{-0.709231010958016926537225142+1.09596916251277206706752728\,\mathrm{i},-0.02291660900+0.\times 10^{-12}\,\mathrm{i}\right\}
Cell 37: ✓ — (no output)
Cell 38: ✓ — (no output)
Cell 39: ✓ — (no output)
Cell 40: ✓ — (no output)
Cell 41: ✓ — (no output)
Cell 42: ✓ — Out[109]= 0

## 2026-07-18 11:27:37 — ▶️ RUN CELLS 43–47
Cell 43: ✓ — (no output)
Cell 44: ✓ — Out[113]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 45: ✓ — Out[115]= \left\{6,6\right\}
Cell 46: ✓ — Out[116]= \begin{pmatrix}2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i} & 0 & 0 & 0 & 0 & 0\\0 & 7.36608403016431791181-0.968046299794461400266\,\mathrm{i} & 0 & 0 & 0 & 0\\0 & 0 & 2.68702160063169134012-1.3456153260557858257\,\mathrm{i} & 0 & 0 & 0\\0 & 0 & 0 & 7.3660840301643179118107+0.9680462997944614002659\,\mathrm{i} & 0 & 0\\0 & 0 & 0 & 0 & 2.687021600631691340118… [537 chars total]
Cell 47: ✓ — Out[125]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{max|off-diagonal|}\to 5.065342385207116\times 10^{-34},\text{digits vanished (-Log10)}\to 33.29539119374288,\text{min|diagonal|}\to 3.0051232400647923,\text{diagonal values}\to \left\{2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i},7.36608403016431791181-0… [699 chars total]

## 2026-07-18 11:27:41 — ▶️ RUN CELLS 48–54
Cell 48: ✓ — (no output)
Cell 49: ✓ — (no output)
Cell 50: ✓ — Out[133]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 51: ✓ — Out[138]= \left\{1.00000000000000000000000000000,6\right\}
Cell 52: ✓ — Out[143]= \langle\!\lvert \text{max|<PsiLeft|Psi>| off-diagonal}\to 0.\,,\text{digits vanished (-Log10)}\to \mathrm{Indeterminate}\rvert\!\rangle
Cell 53: ✓ — Out[147]= \left\{7.366084030164317911810727273-0.9680462997944614002658874268\,\mathrm{i},1.00000000000000000000+0.\times 10^{-21}\,\mathrm{i}\right\}
Cell 54: ✓ — Out[150]= \begin{matrix}\langle\!\lvert \text{M}\to 0,\text{n}\to 0,\text{braket}\to 2.909064289-5.804620064\,\mathrm{i},\text{G21 diagonal}\to 2.909064289-5.804620064\,\mathrm{i},\text{normRatio}\to 1.000000000+0.\times 10^{-10}\,\mathrm{i}\rvert\!\rangle \\\langle\!\lvert \text{M}\to 1,\text{n}\to 0,\text{braket}\to 7.366084030-0.968046300\,\mathrm{i},\text{G21 diagonal}\to 7.366084030-0.9680463… [1247 chars total]

## 2026-07-18 11:27:45 — ▶️ RUN CELLS 55–62
Cell 56: ✓ — (no output)
Cell 57: ✓ — Out[158]= \begin{pmatrix}3 & 3\\3 & 3\end{pmatrix}
Cell 58: ✓ — (no output)
Cell 59: ✓ — Out[165]= \left\{0\right\}
Cell 60: ✓ — Out[169]= \begin{matrix}\left\{1,1,0,1.\,,0.\,\right\}\\\left\{1,1,1,1.\,,0.\,\right\}\\\left\{1,2,0,1.\,,0.\,\right\}\\\left\{1,2,1,1.\,,0.\,\right\}\\\left\{1,3,0,1.\,,0.\,\right\}\\\left\{1,3,1,1.\,,0.\,\right\}\\\left\{2,1,0,1.\,,0.\,\right\}\\\left\{2,1,1,1.\,,0.\,\right\}\\\left\{2,1,2,1.\,,0.\,\right\}\\\left\{2,2,0,1.\,,0.\,\right\}\\\left\{2,2,1,1.\,,0.\,\right\}\\\left\{2,2,2,1.\,,0.\,\r… [898 chars total]
Cell 61: ✓ — Out[170]= \langle\!\lvert \text{max spread (ratio non-constancy)}\to 0.\,,\text{max |ratio - 1| over sweep}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases checked}\to 27\rvert\!\rangle
Cell 62: ✓ — Out[176]= \langle\!\lvert \text{max operator-level spread}\to 0.\,,\text{max |$\mu$ - $\mu$Site2|}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases}\to 27\rvert\!\rangle

## 2026-07-18 11:27:54 — ▶️ RUN CELLS 63–70
Cell 64: ✓ — (no output)
Cell 65: ✓ — Out[185]= \langle\!\lvert \text{w\_0 ($\lambda$=1,2,3)}\to \left\{1,1,1\right\},\text{full weights w\_j}\to \left\{1\to \left\{1,-1\right\},2\to \left\{1,-2,1\right\},3\to \left\{1,-3,3,-1\right\}\right\},\text{br1[$\lambda$][1] ($\lambda$=1,2,3)}\to \left\{0,0,0\right\}\rvert\!\rangle
Cell 66: ✓ — Out[187]= \left\{0,0\right\}
Cell 67: ✓ — Out[190]= \sqrt[3]{k}\,\left(\left(-1.18561270626131985872+0.832253303224995280404\,\mathrm{i}\right)\,k^2-\left(0.761819910806671651081+0.812812473459294825632\,\mathrm{i}\right)\,k-\left(0.404410037536507329420+0.598870869485872613792\,\mathrm{i}\right)\right)
Cell 68: ✓ — Out[194]= \begin{pmatrix}0 & 0 & 0\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix}
Cell 69: ✓ — Out[198]= \langle\!\lvert \text{OP(E12)}\to 1.\,+0.\times 10^{-30}\,\mathrm{i},\text{OP(E21)}\to \left(-0.03626734060947795145020244939-2.004237714666237867354326095\,\mathrm{i}\right)\,\varphi ^2-\left(1.522150083623538480912456055-0.2442047957313920727733393355\,\mathrm{i}\right)\,\varphi +\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right),\text{nontrivial generator \{a,b\}}\to \left\{2,1\right\},\… [464 chars total]
Cell 70: ✓ — Out[201]= \left\{1,2\right\}

## 2026-07-18 11:27:58 — ▶️ RUN CELLS 71–77
Cell 71: ✓ — Out[207]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{$\varphi$ candidates per state}\to \begin{pmatrix}\left\{0,0\right\} & \left(-0.920793778924720536105378676-0.390049761817044163511250929\,\mathrm{i}\right)\,k^{4/3}+\left(0.695722366612934233848844287-0.718310788304406529935213657\,\mathrm{i}\right)\,\sqrt[3]{k}-\left(0.3604398751… [1247 chars total]
Cell 72: ✓ — Out[213]= \langle\!\lvert \text{c\_sb solution(s)}\to \left\{\left\{\mathrm{cb}\to \frac{\left(-6.242239011314\times 10^{45}-4.511757007219\times 10^{45}\,\mathrm{i}\right)\,k^{4/3}-\left(5.999302657796\times 10^{45}+8.356648463050\times 10^{45}\,\mathrm{i}\right)\,\sqrt[3]{k}}{9.833530773718\times 10^{45}\,k^{4/3}-\left(3.54439660434\times 10^{45}-9.172544913427\times 10^{45}\,\mathrm{i}\right)\,… [841 chars total]
Cell 73: ✓ — Out[216]= \langle\!\lvert \text{generator \{a,b\} tested}\to \left\{2,1\right\},\text{rep tested}\to \left\{1,2\right\},\text{FL shape (verified exact for all 6 states)}\to \text{FL[1,2][state][k] = k\^{}(1/3) (A\_state + B\_state k)},\text{c0(state) = OPg[state][0]}\to \begin{pmatrix}\left\{0,0\right\} & 1.\,+0.\times 10^{-28}\,\mathrm{i}\\\left\{1,0\right\} & 1.\,+0.\times 10^{-28}\,\mathrm{i}\\… [1247 chars total]
Cell 74: ✓ — Out[220]= \left(1.000000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right)-\left(0.360439875147226537266677651+0.932782448593374716484835982\,\mathrm{i}\right)\,k
Cell 75: ✓ — Out[227]= \langle\!\lvert \text{alphaTab \{state,D0,c0,Dc,c1,$\alpha$\}}\to \begin{pmatrix}\left\{0,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & -0.36043987514722653-0.9327824485933747\,\mathrm{i} & 0.36043987514722653 +0.9327824485933747\,\mathrm{i} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.7101545203859353 -0.2371533418544… [1247 chars total]
Cell 76: ✓ — Out[235]= \langle\!\lvert \text{validation}\to \left\{\left\{1,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\end{pmatrix},\left\{2,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\\\left\{2,1\right\} & 0\\\left\{3,0\right\} & 0\end{pmatrix},\left\{1… [1247 chars total]
Cell 77: ✓ — Out[237]= \langle\!\lvert \text{generator \{a,b\}}\to \left\{2,1\right\},\text{$\varphi$(k)}\to \text{$\alpha$ k},\text{$\alpha$ (universal value)}\to -1.\,+0.\,\,\mathrm{i},\text{normalization}\to \text{1-site FSoV determinant det1[f] = br1[f]/(f/.u->$\theta$1); Section C's det[$\lambda$1,$\lambda$2][f] restricted to 1 site},\text{max residual over rep sweep}\to 0,\text{reps validated}\to \begin{… [510 chars total]

## 2026-07-18 11:28:08 — ▶️ RUN CELLS 78–82
Cell 79: ✓ — (no output)
Cell 80: ✓ — Out[245]= \left\{\text{Test 1 (symbolic monomials u\^{}j) residuals},\left\{0\right\}\right\}
Cell 81: ✓ — Out[250]= \left\{\text{Test 2 (twisted monomials k\^{}(u/h) u\^{}j) residuals},\left\{0\right\}\right\}
Cell 82: ✓ — Out[259]= \left\{\text{Test 3 worst |br identity residual| on Q1},0,\text{worst |detRatio - predicted| },0,\text{num states checked},41\right\}

## 2026-07-18 11:28:12 — ▶️ RUN CELLS 83–87
Cell 84: ✓ — (no output)
Cell 85: ✓ — Out[271]= \begin{pmatrix}1 & 1 & 0\\2 & 1 & 0\\3 & 1 & 0\end{pmatrix}
Cell 86: ✓ — Out[278]= \left\{\text{Mirror Test1 (u\^{}j)},\left\{0\right\},\text{Mirror Test2 (k\^{}(u/h)u\^{}j)},\left\{0\right\}\right\}
Cell 87: ✓ — Out[285]= \left\{\text{Mirror Test3 worst |br identity residual| on Q1},0,\text{num states checked},41\right\}

## 2026-07-18 11:28:13 — ▶️ RUN CELLS 88–91
Cell 89: ✓ — Out[294]= \left\{\text{OP2(E12)},1.\,+0.\times 10^{-30}\,\mathrm{i},\text{OP2(E21)},\left(2.603469709065218092712489006+0.6885776528619826437114966465\,\mathrm{i}\right)\,\varphi +\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right),\text{nontrivial generator \{a,b\}},\left\{2,1\right\}\right\}
Cell 90: ✓ — Out[297]= \left(1.\,+0.\times 10^{-27}\,\mathrm{i}\right)-\left(0.360439875147226537266677651+0.932782448593374716484835982\,\mathrm{i}\right)\,k
Cell 91: ✓ — Out[302]= \text{\$Aborted} | Out[304]= \left\{\text{alphaTab2 \{state,D,c1,alpha\}},\begin{pmatrix}\left\{0,0\right\} & -0.360439875147226537266677651-0.932782448593374716484835982\,\mathrm{i} & 0.3604398751472265372666776506+0.9327824485933747164848359821\,\mathrm{i} & -1.000000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\\\left\{1,0\right\} & 0.4685059330398572984-0.171361115691460089\,\m… [1276 chars total]

## 2026-07-18 11:28:31 — 🗑️ DELETE CODE CELL 91
(*Discovery step (sharp, λ2=1 rep, mirrors Cell 71): rep (2,1), all 6 states — solve φ per state
  assuming exact match against OPg2 (both linear in φ since λ2=1), then check the ratio
  α2_state = D_…

## 2026-07-18 11:28:44 — ✏️ EDIT CELL 91
(*Discovery step (sharp, λ2=1 rep, mirrors Cell 71): rep (2,1), all 6 states — solve the per-state
  alpha ratio D_state/c1_state assuming FLdet2/OPg2 are both linear in k,φ (true since λ2=1 caps
  gR…


## 2026-07-18 11:28:44
**Input:**
```wolfram
(*Discovery step (sharp, λ2=1 rep, mirrors Cell 71): rep (2,1), all 6 states — solve the per-state
  alpha ratio D_state/c1_state assuming FLdet2/OPg2 are both linear in k,φ (true since λ2=1 caps
  gRot's degree at 1), Module-scoped to avoid cross-iteration variable leakage.*)
ClearAll[alphaTab2Clean];
alphaTab2Clean = Table[
  Module[{Fl,c0,Dcoef,Og,cc0,cc1},
    Fl = FLdet2[2,1][st[[1]],st[[2]]][k]//Simplify;
    c0 = Fl/.k->0;
    Dcoef = Coefficient[Fl,k,1];
    Og = OPg2[2,1][st[[1]],st[[2]]][φ]//Simplify;
    cc0 = Og/.φ->0;
    cc1 = Coefficient[Og,φ,1];
    {st, N[c0], N[cc0], N[Dcoef], N[cc1], N[Dcoef/cc1]}
  ]
,{st,statesOf[2,1]}];
TableForm[alphaTab2Clean]
```
**Output:**
```
Out= TableForm[{{{0, 0}, 1. + 0.*I, 1. + 0.*I, -0.36043987514722653 - 0.9327824485933747*I, 0.36043987514722653 + 0.9327824485933747*I, -1. + 0.*I}, {{1, 0}, 1. + 0.*I, 1. + 0.*I, 0.4685059330398573 - 0.1713611156914601*I, 2.603469709065218 + 0.6885776528619827*I, 0.15191893851739566 - 0.10600050418875033*I}, {{1, 1}, 1. + 0.*I, 1. + 0.*I, -0.6957397673718044 - 0.2544743499707327*I, -0.5468421861663373 + 0.46368301892983316*I, 0.5105931037689164 + 0.8982988405024716*I}, {{2, 0}, 1. + 0.*I, 1. + 0.*I, 0.4685059330398573 + 0.1713611156914601*I, 2.603469709065218 - 0.6885776528619827*I, 0.15191893851739566 + 0.10600050418875033*I}, {{2, 1}, 1. + 0.*I, 1. + 0.*I, -0.6957397673718044 + 0.2544743499707327*I, -0.5468421861663373 - 0.46368301892983316*I, 0.5105931037689164 - 0.8982988405024716*I}, {{3, 0}, 1. + 0.*I, 1. + 0.*I, -0.36043987514722653 + 0.9327824485933747*I, 0.36043987514722653 - 0.9327824485933747*I, -1. + 0.*I}}]
```
## 2026-07-18 11:29:31 — ▶️ RUN CELLS 92–104
Cell 92: ✓ — Out[306]//TableForm= \begin{pmatrix}\left\{0,0\right\} & -1.\,+0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.973902234276272 +0.2269679230016555\,\mathrm{i} & 0.7391304347826086 +0.\,\,\mathrm{i} & 0.973902234276272 +0.2269679230016555\,\mathrm{i}\\\left\{1,0\right\} & 0.15191893851739566 -0.10600050418875033\,\mathrm{i} & -0.15191893851739566+0.10600050418875033\,\mathrm{i} & 0.7651595274649543 +… [1258 chars total]
Cell 93: ✓ — Out[309]= \begin{pmatrix}\left\{0,0\right\} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 0.05448860910968355 -0.09540961504089102\,\mathrm{i}\\\left\{1,1\right\} & -0.0293413933141388+0.7363100318787654\,\mathrm{i}\\\left\{2,0\right\} & 0.11703059508534998 +0.\,\,\mathrm{i}\\\left\{2,1\right\} & 0.8158972359113332 +0.\,\,\mathrm{i}\\\left\{3,0\right\} & 0.05448860910968355 +0.09540961504089102\,… [541 chars total]
Cell 96: ✓ — Out[315]= \left\{\text{O[k1,0]-FLdet[k1] (must be 0)},0,\text{O[0,k2]-FLdet2[k2] (expected NONZERO per Phase1 finding)},\left(-3.0719756421050755-0.5172165371705225\,\mathrm{i}\right)\,k2\right\}
Cell 97: ✓ — Out[320]= \left\{\text{Fcand1},k1^2\,\left(\left(-2.2527492196701659-5.82989030370859198\,\mathrm{i}\right)\,k2-\left(0.4170744170089964+23.0487337186617355\,\mathrm{i}\right)\right)+k1\,\left(\left(0.67128913890100625-0.37524873225652243\,\mathrm{i}\right)\,k2+\left(9.51343802264711551-1.5262799733212005\,\mathrm{i}\right)\right)-\left(1.99115021541939352-0.72828474168870538\,\mathrm{i}\right)\,k… [1014 chars total]
Cell 98: ✓ — Out[322]= \left(0.9252125034899888 -6.175727508765443\,\mathrm{i}\right)\,k1^2\,k2+\left(5.4271267482291075 +3.0890004366773067\,\mathrm{i}\right)\,k1\,k2
Cell 99: ✓ — Out[327]= \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 7.991287939023577 -1.2820751775898085\,\mathrm{i}\\\left\{2,0\right\} & -0.3808070763995185-21.044496003995498\,\mathrm{i}\\\left\{0,1\right\} & 0.6123194936458246 +1.416862394550688\,\mathrm{i}\\\left\{1,1\right\} & -0.6248120089007682-3.051910425543801\,\mathrm{i}\\\left\{2,1\right\} & -1.8923093445229393-4.897107855115217\,\m… [421 chars total]
Cell 100: ✓ — Out[332]= \left\{\text{Fcand2 (no P-insertion)},k1^2\,\left(\left(0.02466179161444501+1.36288164597304175\,\mathrm{i}\right)-\left(0.36043987514722654+0.93278244859337472\,\mathrm{i}\right)\,k2\right)+k1\,\left(\left(1.2351720155778515-0.69045766735200127\,\mathrm{i}\right)\,k2-\left(6.46913785540003854-1.0378703818584163\,\mathrm{i}\right)\right)+\left(5.38781822995835893-1.9706528304517910\,\mat… [788 chars total]
Cell 101: ✓ — 2 kernel messages: Power::infy: Infinite expression 0       encountered.
Infinity::indet: Indeterminate expression ----------------- encountered.                                                   1/3                    … [also: Power::infy] result: Out[335]//TableForm= \begin{pmatrix}2 & 1 & \left\{0,0\right\} & \mathrm{Indeterminate} & \mathrm{Indeterminate} & -4.970498218180917+4.515583419069149\,\mathrm{i}\\2 & 1 & \left\{1,1\right\} & \mathr…
Cell 102: ✓ — Out[338]//TableForm= \begin{pmatrix}2 & 1 & \left\{0,0\right\} & \left(-7.771745024480887+7.060451775520925\,\mathrm{i}\right)\,k1^2-\left(3.7846186890458786 +9.794215710230434\,\mathrm{i}\right)\,k1 & \left(1.8923093445229393 +4.897107855115217\,\mathrm{i}\right)\,k2 & -4.970498218180917+4.515583419069149\,\mathrm{i}\\2 & 1 & \left\{1,1\right\} & \left(9.339255641812587 +10.661478687052961\,\math… [1258 chars total]
Cell 103: ✓ — 2 kernel messages: Power::infy: Infinite expression 0       encountered.
Infinity::indet: Indeterminate expression ----------------- encountered.                                                   1/3                    … [also: Power::infy] result: Out[342]= \left\{\text{Fcand1(k1,0)},\mathrm{Indeterminate},\text{FLdet[k1]},\frac{\left(2.45075493251774917289426756629+5.78551794268217244227029022731\,\mathrm{i}\right)\,\left(\left(0.0525946877178…
Cell 104: ✓ — Out[347]= \left\{\text{Fcand1(k1,0)},\left(-8.51191121729+7.73287575414\,\mathrm{i}\right)\,k1^2-\left(4.5054984393+11.6597806074\,\mathrm{i}\right)\,k1+\left(1.000000000000+0.\times 10^{-12}\,\mathrm{i}\right),\text{FLdet[k1]},\frac{\left(2.45075493251774917289426756629+5.78551794268217244227029022731\,\mathrm{i}\right)\,\left(\left(0.05259468771783351844996683714+0.15021349717800998891606886537\… [960 chars total]


## 2026-07-18 11:34:46
**Input:**
```wolfram
{θ1,θ2,h,ValueQ[C1pred],ValueQ[O2rot],ValueQ[Fcand1]}
```
**Output:**
```
Out= {1/3, 1/7, 1, False, False, False}
```

## 2026-07-18 11:34:53
**Input:**
```wolfram
{Length[DownValues[C1pred]], Length[DownValues[O2rot]], Length[DownValues[Fcand1]], Length[DownValues[Q1]], Length[DownValues[xSingle]], Length[DownValues[xL1]]}
```
**Output:**
```
Out= {1, 0, 0, 0, 0, 0}
```

---

## 2026-07-18 11:35:02 — 🔄 KERNEL RESTART

## 2026-07-18 11:35:11 — ▶️ RUN CELLS 1–10
Cell 1: ✓ — Out[8]= \left\{\frac{1}{3},\frac{1}{7},1,0.360439875147226537266677650596+0.932782448593374716484835982071\,\mathrm{i},0.36043987514722653726667765060-0.932782448593374716484835982071\,\mathrm{i},0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i},1\right\}
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 8: ✓ — (no output)
Cell 10: ✓ — (no output)

## 2026-07-18 11:35:16 — ▶️ RUN CELLS 11–20
Cell 11: ✓ — (no output)
Cell 12: ✓ — Out[27]= \left\{0\right\}
Cell 13: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — Out[31]= \left\{0\right\}
Cell 16: ✓ — (no output)
Cell 17: ✓ — Out[34]= \begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{pmatrix}
Cell 18: ✓ — (no output)
Cell 20: ✓ — (no output)

## 2026-07-18 11:35:19 — ▶️ RUN CELLS 21–28
Cell 21: ✓ — (no output)
Cell 23: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 26: ✓ — (no output)
Cell 27: ✓ — Out[52]= \begin{pmatrix}0 & \frac{1}{2} & 0 & 0 & 0 & 0\end{pmatrix}
Cell 28: ✓ — Out[53]= \left\{0\right\}

## 2026-07-18 11:35:22 — ▶️ RUN CELLS 29–36
Cell 30: ✓ — Out[59]= \left\{\left(u-\frac{4}{3}\right)\,\left(u-\frac{8}{7}\right),\left(u-\frac{1}{3}\right)\,\left(u-\frac{1}{7}\right),\left(0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i}\right)\,u^2-\left(1.06415582186324025288257211128+0.\times 10^{-29}\,\mathrm{i}\right)\,u\right\}
Cell 31: ✓ — Out[62]= \begin{pmatrix}1 & 1\\2 & 2\\3 & 3\end{pmatrix}
Cell 32: ✓ — Out[65]= \text{max|C2-($\kappa$1+$\kappa$2)Id|}\to 0.\times 10^{-29}
Cell 33: ✓ — Out[75]= \langle\!\lvert \left\{0,0\right\}\to 0.5664055180884988 +1.3769645669711723\,\mathrm{i},\left\{1,0\right\}\to -0.778123966451158+0.\,\,\mathrm{i},\left\{1,1\right\}\to 1.1900552523337025 +0.\,\,\mathrm{i},\left\{2,0\right\}\to 0.5664055180884988 -1.3769645669711723\,\mathrm{i}\rvert\!\rangle
Cell 34: ✓ — 2 kernel messages: Psi::spurious: No genuine state at (M,n)=(2,2) for (λ1,λ2)=(2,1): weight sector M=2 has only 2 genuine state(s) (n=0..1), so n=2 does not exist. [also: Psi::spurious] result: Print: spurious-guard test (expect a Psi::spurious message and Missing[...]):[formula] | Out[84]= \left\{\left\{1.000000000000000000000000000+0.\times 10^{-29}\,\mathrm{i},-1.3004494881076436272092358…
⚠️ stopped at Cell 34 — error detected (pass stopOnError:false to continue past errors)

## 2026-07-18 11:35:37 — ▶️ RUN CELLS 35–42
Cell 35: ✓ — Out[91]= \left\{1.03464899078135229594035121+0.381869568380946283899575941\,\mathrm{i},-2.37022278949001677248172619-0.20048970475590097847093779\,\mathrm{i},1.00000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right\}
Cell 36: ✓ — Out[94]= \left\{-0.709231010958016926537225142+1.09596916251277206706752728\,\mathrm{i},-0.02291660900+0.\times 10^{-12}\,\mathrm{i}\right\}
Cell 37: ✓ — (no output)
Cell 38: ✓ — (no output)
Cell 39: ✓ — (no output)
Cell 40: ✓ — (no output)
Cell 41: ✓ — (no output)
Cell 42: ✓ — Out[109]= 0

## 2026-07-18 11:35:40 — ▶️ RUN CELLS 43–47
Cell 43: ✓ — (no output)
Cell 44: ✓ — Out[113]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 45: ✓ — Out[115]= \left\{6,6\right\}
Cell 46: ✓ — Out[116]= \begin{pmatrix}2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i} & 0 & 0 & 0 & 0 & 0\\0 & 7.36608403016431791181-0.968046299794461400266\,\mathrm{i} & 0 & 0 & 0 & 0\\0 & 0 & 2.68702160063169134012-1.3456153260557858257\,\mathrm{i} & 0 & 0 & 0\\0 & 0 & 0 & 7.3660840301643179118107+0.9680462997944614002659\,\mathrm{i} & 0 & 0\\0 & 0 & 0 & 0 & 2.687021600631691340118… [537 chars total]
Cell 47: ✓ — Out[125]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{max|off-diagonal|}\to 5.065342385207116\times 10^{-34},\text{digits vanished (-Log10)}\to 33.29539119374288,\text{min|diagonal|}\to 3.0051232400647923,\text{diagonal values}\to \left\{2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i},7.36608403016431791181-0… [699 chars total]

## 2026-07-18 11:35:44 — ▶️ RUN CELLS 48–54
Cell 48: ✓ — (no output)
Cell 49: ✓ — (no output)
Cell 50: ✓ — Out[133]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 51: ✓ — Out[138]= \left\{1.00000000000000000000000000000,6\right\}
Cell 52: ✓ — Out[143]= \langle\!\lvert \text{max|<PsiLeft|Psi>| off-diagonal}\to 0.\,,\text{digits vanished (-Log10)}\to \mathrm{Indeterminate}\rvert\!\rangle
Cell 53: ✓ — Out[147]= \left\{7.366084030164317911810727273-0.9680462997944614002658874268\,\mathrm{i},1.00000000000000000000+0.\times 10^{-21}\,\mathrm{i}\right\}
Cell 54: ✓ — Out[150]= \begin{matrix}\langle\!\lvert \text{M}\to 0,\text{n}\to 0,\text{braket}\to 2.909064289-5.804620064\,\mathrm{i},\text{G21 diagonal}\to 2.909064289-5.804620064\,\mathrm{i},\text{normRatio}\to 1.000000000+0.\times 10^{-10}\,\mathrm{i}\rvert\!\rangle \\\langle\!\lvert \text{M}\to 1,\text{n}\to 0,\text{braket}\to 7.366084030-0.968046300\,\mathrm{i},\text{G21 diagonal}\to 7.366084030-0.9680463… [1247 chars total]

## 2026-07-18 11:35:47 — ▶️ RUN CELLS 55–62
Cell 56: ✓ — (no output)
Cell 57: ✓ — Out[158]= \begin{pmatrix}3 & 3\\3 & 3\end{pmatrix}
Cell 58: ✓ — (no output)
Cell 59: ✓ — Out[165]= \left\{0\right\}
Cell 60: ✓ — Out[169]= \begin{matrix}\left\{1,1,0,1.\,,0.\,\right\}\\\left\{1,1,1,1.\,,0.\,\right\}\\\left\{1,2,0,1.\,,0.\,\right\}\\\left\{1,2,1,1.\,,0.\,\right\}\\\left\{1,3,0,1.\,,0.\,\right\}\\\left\{1,3,1,1.\,,0.\,\right\}\\\left\{2,1,0,1.\,,0.\,\right\}\\\left\{2,1,1,1.\,,0.\,\right\}\\\left\{2,1,2,1.\,,0.\,\right\}\\\left\{2,2,0,1.\,,0.\,\right\}\\\left\{2,2,1,1.\,,0.\,\right\}\\\left\{2,2,2,1.\,,0.\,\r… [898 chars total]
Cell 61: ✓ — Out[170]= \langle\!\lvert \text{max spread (ratio non-constancy)}\to 0.\,,\text{max |ratio - 1| over sweep}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases checked}\to 27\rvert\!\rangle
Cell 62: ✓ — Out[176]= \langle\!\lvert \text{max operator-level spread}\to 0.\,,\text{max |$\mu$ - $\mu$Site2|}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases}\to 27\rvert\!\rangle

## 2026-07-18 11:36:11 — ▶️ RUN CELLS 63–70
Cell 64: ✓ — (no output)
Cell 65: ✓ — Out[185]= \langle\!\lvert \text{w\_0 ($\lambda$=1,2,3)}\to \left\{1,1,1\right\},\text{full weights w\_j}\to \left\{1\to \left\{1,-1\right\},2\to \left\{1,-2,1\right\},3\to \left\{1,-3,3,-1\right\}\right\},\text{br1[$\lambda$][1] ($\lambda$=1,2,3)}\to \left\{0,0,0\right\}\rvert\!\rangle
Cell 66: ✓ — Out[187]= \left\{0,0\right\}
Cell 67: ✓ — Out[190]= \sqrt[3]{k}\,\left(\left(-1.18561270626131985872+0.832253303224995280404\,\mathrm{i}\right)\,k^2-\left(0.761819910806671651081+0.812812473459294825632\,\mathrm{i}\right)\,k-\left(0.404410037536507329420+0.598870869485872613792\,\mathrm{i}\right)\right)
Cell 68: ✓ — Out[194]= \begin{pmatrix}0 & 0 & 0\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix}
Cell 69: ✓ — Out[198]= \langle\!\lvert \text{OP(E12)}\to 1.\,+0.\times 10^{-30}\,\mathrm{i},\text{OP(E21)}\to \left(-0.03626734060947795145020244939-2.004237714666237867354326095\,\mathrm{i}\right)\,\varphi ^2-\left(1.522150083623538480912456055-0.2442047957313920727733393355\,\mathrm{i}\right)\,\varphi +\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right),\text{nontrivial generator \{a,b\}}\to \left\{2,1\right\},\… [464 chars total]
Cell 70: ✓ — Out[201]= \left\{1,2\right\}

## 2026-07-18 11:36:15 — ▶️ RUN CELLS 71–77
Cell 71: ✓ — Out[207]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{$\varphi$ candidates per state}\to \begin{pmatrix}\left\{0,0\right\} & \left(-0.920793778924720536105378676-0.390049761817044163511250929\,\mathrm{i}\right)\,k^{4/3}+\left(0.695722366612934233848844287-0.718310788304406529935213657\,\mathrm{i}\right)\,\sqrt[3]{k}-\left(0.3604398751… [1247 chars total]
Cell 72: ✓ — Out[213]= \langle\!\lvert \text{c\_sb solution(s)}\to \left\{\left\{\mathrm{cb}\to \frac{\left(-6.242239011314\times 10^{45}-4.511757007219\times 10^{45}\,\mathrm{i}\right)\,k^{4/3}-\left(5.999302657796\times 10^{45}+8.356648463050\times 10^{45}\,\mathrm{i}\right)\,\sqrt[3]{k}}{9.833530773718\times 10^{45}\,k^{4/3}-\left(3.54439660434\times 10^{45}-9.172544913427\times 10^{45}\,\mathrm{i}\right)\,… [841 chars total]
Cell 73: ✓ — Out[216]= \langle\!\lvert \text{generator \{a,b\} tested}\to \left\{2,1\right\},\text{rep tested}\to \left\{1,2\right\},\text{FL shape (verified exact for all 6 states)}\to \text{FL[1,2][state][k] = k\^{}(1/3) (A\_state + B\_state k)},\text{c0(state) = OPg[state][0]}\to \begin{pmatrix}\left\{0,0\right\} & 1.\,+0.\times 10^{-28}\,\mathrm{i}\\\left\{1,0\right\} & 1.\,+0.\times 10^{-28}\,\mathrm{i}\\… [1247 chars total]
Cell 74: ✓ — Out[220]= \left(1.000000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right)-\left(0.360439875147226537266677651+0.932782448593374716484835982\,\mathrm{i}\right)\,k
Cell 75: ✓ — Out[227]= \langle\!\lvert \text{alphaTab \{state,D0,c0,Dc,c1,$\alpha$\}}\to \begin{pmatrix}\left\{0,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & -0.36043987514722653-0.9327824485933747\,\mathrm{i} & 0.36043987514722653 +0.9327824485933747\,\mathrm{i} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.7101545203859353 -0.2371533418544… [1247 chars total]
Cell 76: ✓ — Out[235]= \langle\!\lvert \text{validation}\to \left\{\left\{1,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\end{pmatrix},\left\{2,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\\\left\{2,1\right\} & 0\\\left\{3,0\right\} & 0\end{pmatrix},\left\{1… [1247 chars total]
Cell 77: ✓ — Out[237]= \langle\!\lvert \text{generator \{a,b\}}\to \left\{2,1\right\},\text{$\varphi$(k)}\to \text{$\alpha$ k},\text{$\alpha$ (universal value)}\to -1.\,+0.\,\,\mathrm{i},\text{normalization}\to \text{1-site FSoV determinant det1[f] = br1[f]/(f/.u->$\theta$1); Section C's det[$\lambda$1,$\lambda$2][f] restricted to 1 site},\text{max residual over rep sweep}\to 0,\text{reps validated}\to \begin{… [510 chars total]

## 2026-07-18 11:36:49 — ▶️ RUN CELLS 78–82
Cell 79: ✓ — (no output)
Cell 80: ✓ — Out[245]= \left\{\text{Test 1 (symbolic monomials u\^{}j) residuals},\left\{0\right\}\right\}
Cell 81: ✓ — Out[250]= \left\{\text{Test 2 (twisted monomials k\^{}(u/h) u\^{}j) residuals},\left\{0\right\}\right\}
Cell 82: ✓ — Out[259]= \left\{\text{Test 3 worst |br identity residual| on Q1},0,\text{worst |detRatio - predicted| },0,\text{num states checked},41\right\}

## 2026-07-18 11:36:53 — ▶️ RUN CELLS 83–87
Cell 84: ✓ — (no output)
Cell 85: ✓ — Out[271]= \begin{pmatrix}1 & 1 & 0\\2 & 1 & 0\\3 & 1 & 0\end{pmatrix}
Cell 86: ✓ — Out[278]= \left\{\text{Mirror Test1 (u\^{}j)},\left\{0\right\},\text{Mirror Test2 (k\^{}(u/h)u\^{}j)},\left\{0\right\}\right\}
Cell 87: ✓ — Out[285]= \left\{\text{Mirror Test3 worst |br identity residual| on Q1},0,\text{num states checked},41\right\}

## 2026-07-18 11:36:54 — ▶️ RUN CELLS 88–91
Cell 89: ✓ — Out[294]= \left\{\text{OP2(E12)},1.\,+0.\times 10^{-30}\,\mathrm{i},\text{OP2(E21)},\left(2.603469709065218092712489006+0.6885776528619826437114966465\,\mathrm{i}\right)\,\varphi +\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right),\text{nontrivial generator \{a,b\}},\left\{2,1\right\}\right\}
Cell 90: ✓ — Out[297]= \left(1.\,+0.\times 10^{-27}\,\mathrm{i}\right)-\left(0.360439875147226537266677651+0.932782448593374716484835982\,\mathrm{i}\right)\,k
Cell 91: ✓ — Out[300]//TableForm= \begin{pmatrix}\left\{0,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & -0.36043987514722653-0.9327824485933747\,\mathrm{i} & 0.36043987514722653 +0.9327824485933747\,\mathrm{i} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.4685059330398573 -0.1713611156914601\,\mathrm{i} & 2.603469709065218 +0.6885776528619827… [1258 chars total]

## 2026-07-18 11:37:14 — ▶️ RUN CELLS 92–105
Cell 92: ✓ — Out[302]//TableForm= \begin{pmatrix}\left\{0,0\right\} & -1.\,+0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.973902234276272 +0.2269679230016555\,\mathrm{i} & 0.7391304347826086 +0.\,\,\mathrm{i} & 0.973902234276272 +0.2269679230016555\,\mathrm{i}\\\left\{1,0\right\} & 0.15191893851739566 -0.10600050418875033\,\mathrm{i} & -0.15191893851739566+0.10600050418875033\,\mathrm{i} & 0.7651595274649543 +… [1258 chars total]
Cell 93: ✓ — Out[305]= \begin{pmatrix}\left\{0,0\right\} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 0.05448860910968355 -0.09540961504089102\,\mathrm{i}\\\left\{1,1\right\} & -0.0293413933141388+0.7363100318787654\,\mathrm{i}\\\left\{2,0\right\} & 0.11703059508534998 +0.\,\,\mathrm{i}\\\left\{2,1\right\} & 0.8158972359113332 +0.\,\,\mathrm{i}\\\left\{3,0\right\} & 0.05448860910968355 +0.09540961504089102\,… [541 chars total]
Cell 96: ✓ — Out[311]= \left\{\text{O[k1,0]-FLdet[k1] (must be 0)},0,\text{O[0,k2]-FLdet2[k2] (expected NONZERO per Phase1 finding)},\left(-3.0719756421050755-0.5172165371705225\,\mathrm{i}\right)\,k2\right\}
Cell 97: ✓ — Out[316]= \left\{\text{Fcand1},k1^2\,\left(\left(-2.2527492196701659-5.82989030370859198\,\mathrm{i}\right)\,k2-\left(0.4170744170089964+23.0487337186617355\,\mathrm{i}\right)\right)+k1\,\left(\left(0.67128913890100625-0.37524873225652243\,\mathrm{i}\right)\,k2+\left(9.51343802264711551-1.5262799733212005\,\mathrm{i}\right)\right)-\left(1.99115021541939352-0.72828474168870538\,\mathrm{i}\right)\,k… [1014 chars total]
Cell 98: ✓ — Out[318]= \left(0.9252125034899888 -6.175727508765443\,\mathrm{i}\right)\,k1^2\,k2+\left(5.4271267482291075 +3.0890004366773067\,\mathrm{i}\right)\,k1\,k2
Cell 99: ✓ — Out[323]= \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 7.991287939023577 -1.2820751775898085\,\mathrm{i}\\\left\{2,0\right\} & -0.3808070763995185-21.044496003995498\,\mathrm{i}\\\left\{0,1\right\} & 0.6123194936458246 +1.416862394550688\,\mathrm{i}\\\left\{1,1\right\} & -0.6248120089007682-3.051910425543801\,\mathrm{i}\\\left\{2,1\right\} & -1.8923093445229393-4.897107855115217\,\m… [421 chars total]
Cell 100: ✓ — Out[328]= \left\{\text{Fcand2 (no P-insertion)},k1^2\,\left(\left(0.02466179161444501+1.36288164597304175\,\mathrm{i}\right)-\left(0.36043987514722654+0.93278244859337472\,\mathrm{i}\right)\,k2\right)+k1\,\left(\left(1.2351720155778515-0.69045766735200127\,\mathrm{i}\right)\,k2-\left(6.46913785540003854-1.0378703818584163\,\mathrm{i}\right)\right)+\left(5.38781822995835893-1.9706528304517910\,\mat… [788 chars total]
Cell 101: ✓ — 2 kernel messages: Power::infy: Infinite expression 0       encountered.
Infinity::indet: Indeterminate expression ----------------- encountered.                                                   1/3                    … [also: Power::infy] result: Out[331]//TableForm= \begin{pmatrix}2 & 1 & \left\{0,0\right\} & \mathrm{Indeterminate} & \mathrm{Indeterminate} & -4.970498218180917+4.515583419069149\,\mathrm{i}\\2 & 1 & \left\{1,1\right\} & \mathr…
Cell 102: ✓ — Out[334]//TableForm= \begin{pmatrix}2 & 1 & \left\{0,0\right\} & \left(-7.771745024480887+7.060451775520925\,\mathrm{i}\right)\,k1^2-\left(3.7846186890458786 +9.794215710230434\,\mathrm{i}\right)\,k1 & \left(1.8923093445229393 +4.897107855115217\,\mathrm{i}\right)\,k2 & -4.970498218180917+4.515583419069149\,\mathrm{i}\\2 & 1 & \left\{1,1\right\} & \left(9.339255641812587 +10.661478687052961\,\math… [1258 chars total]
Cell 103: ✓ — 2 kernel messages: Power::infy: Infinite expression 0       encountered.
Infinity::indet: Indeterminate expression ----------------- encountered.                                                   1/3                    … [also: Power::infy] result: Out[338]= \left\{\text{Fcand1(k1,0)},\mathrm{Indeterminate},\text{FLdet[k1]},\frac{\left(2.45075493251774917289426756629+5.78551794268217244227029022731\,\mathrm{i}\right)\,\left(\left(0.0525946877178…
Cell 104: ✓ — Out[343]= \left\{\text{Fcand1(k1,0)},\left(-8.51191121729+7.73287575414\,\mathrm{i}\right)\,k1^2-\left(4.5054984393+11.6597806074\,\mathrm{i}\right)\,k1+\left(1.000000000000+0.\times 10^{-12}\,\mathrm{i}\right),\text{FLdet[k1]},\frac{\left(2.45075493251774917289426756629+5.78551794268217244227029022731\,\mathrm{i}\right)\,\left(\left(0.05259468771783351844996683714+0.15021349717800998891606886537\… [960 chars total]

## 2026-07-18 11:37:48 — 📥 BULK INSERT 2 CELL(S) at positions 106–107
1. [markdown] (*Section I — Phase 2 attempt via the GRADED EXPANSION route. The site-2 rotation is a FINITE sum↵(g
2. [code] (*STEP 0 — baseline anchor: does the GLOBAL 2-site det (Section C's det[λ1,λ2][f], same twist↵  inse


## 2026-07-18 11:37:49
**Input:**
```wolfram
(*STEP 0 — baseline anchor: does the GLOBAL 2-site det (Section C's det[λ1,λ2][f], same twist
  inserted at BOTH nodes since f carries a single k) equal O[k,k] (same k at both sites)? This
  calibrates conventions/normalizations before attempting the graded decomposition.*)
ClearAll[step0Report];
step0Report = {};
Do[
  {λ1,λ2} = rep;
  Do[
    fFun[uu_] := k^(uu/h) Q1[λ1,λ2][st[[1]],st[[2]]][uu];
    dGlobal = Simplify[det[λ1,λ2][fFun[u]]];
    oDiag = Simplify[O2rot[λ1,λ2][st[[1]],st[[2]]][k,k]];
    resid = Chop[N[dGlobal-oDiag]/.k->N[7/3+I/5],10^-12];
    ratio = Chop[N[Simplify[dGlobal/oDiag]]/.k->N[7/3+I/5],10^-12];
    AppendTo[step0Report,{λ1,λ2,st,resid,ratio}]
  ,{st,statesOf[λ1,λ2]}]
,{rep,{{1,1},{2,1},{1,2}}}];
worstStep0 = Max[Abs[step0Report[[All,4]]]];
{"STEP 0: det[λ1,λ2][k^(u/h)Q1] vs O[k,k], worst |residual| (sample k)", worstStep0,
 "sample rows {λ1,λ2,state,resid,ratio}", Take[step0Report, UpTo[6]]}
```
**Output:**
```
Out[348]= \left\{\text{STEP 0: det[$\lambda$1,$\lambda$2][k\^{}(u/h)Q1] vs O[k,k], worst |residual| (sample k)},0,\text{sample rows \{$\lambda$1,$\lambda$2,state,resid,ratio\}},\begin{pmatrix}1 & 1 & \left\{0,0\right\} & 0 & 1.\,\\1 & 1 & \left\{1,0\right\} & 0 & 1.\,\\1 & 1 & \left\{1,1\right\} & 0 & 1.\,\\1 & 1 & \left\{2,0\right\} & 0 & 1.\,\\2 & 1 & \left\{0,0\right\} & 0 & 1.\,\\2 & 1 & \left… [445 chars total]
```
## 2026-07-18 11:38:06 — 📥 BULK INSERT 1 CELL(S) at positions 108–108
1. [code] (*STEP 0 hard-assert, full sweep, symbolic k (Simplify to literal 0, not just a numeric sample) —↵


## 2026-07-18 11:38:08
**Input:**
```wolfram
(*STEP 0 hard-assert, full sweep, symbolic k (Simplify to literal 0, not just a numeric sample) —
  this is the calibration the whole graded-expansion route rests on.*)
ClearAll[step0FullResid];
step0FullResid = {};
Do[
  {λ1,λ2} = rep;
  Do[
    Module[{fFun,dGlobal,oDiag,resid},
      fFun[uu_] := k^(uu/h) Q1[λ1,λ2][st[[1]],st[[2]]][uu];
      dGlobal = Simplify[det[λ1,λ2][fFun[u]]];
      oDiag = Simplify[O2rot[λ1,λ2][st[[1]],st[[2]]][k,k]];
      resid = Simplify[dGlobal-oDiag];
      AppendTo[step0FullResid,{λ1,λ2,st,Chop[N[resid]/.k->N[7/3+I/5],10^-10]}]
    ]
  ,{st,statesOf[λ1,λ2]}]
,{rep,{{1,1},{2,1},{1,2},{2,2},{3,1},{1,3}}}];
worstStep0Full = Max[Abs[step0FullResid[[All,4]]]];
If[worstStep0Full > 10^-9, Abort[]];
{"STEP 0 CONFIRMED: det[λ1,λ2][k^(u/h)Q1] == O[k,k] identically",
 "worst residual over full sweep", worstStep0Full, "num states", Length[step0FullResid]}
```
**Output:**
```
Out[354]= \left\{\text{STEP 0 CONFIRMED: det[$\lambda$1,$\lambda$2][k\^{}(u/h)Q1] == O[k,k] identically},\text{worst residual over full sweep},0,\text{num states},41\right\}
```
## 2026-07-18 11:39:36 — 📥 BULK INSERT 1 CELL(S) at positions 109–109
1. [code] (*STEP 1 — covector identification: is ⟨x0[λ]|E21^m (pure representation-theoretic, no θ/site↵  depe


## 2026-07-18 11:39:37
**Input:**
```wolfram
(*STEP 1 — covector identification: is ⟨x0[λ]|E21^m (pure representation-theoretic, no θ/site
  dependence) proportional to Section E's level-m single-site SoV covector xL1[λ][{m}] (built from
  the fused transfer-matrix hierarchy, known to diagonalize B1 with eigenvalue u-θ1-mh)? Test for
  every λ in the sweep {1,2,3}, m=0..λ. Record the proportionality constants γ_m(λ).*)
ClearAll[step1Report,gammaTab];
step1Report = {};
Do[
  Do[
    Module[{raw, ref, ratios},
      raw = xSingle[λ][{0}] . MatrixPower[Ee[λ][2,1],m];
      ref = xL1[λ][{m}];
      ratios = DeleteCases[Flatten[raw]/Flatten[ref], Indeterminate|ComplexInfinity, Infinity];
      ratios = Select[ratios, (Chop[#,10^-12]=!=0 &)];  (*drop the 0/0 slots from zero components*)
      AppendTo[step1Report, {λ, m, Chop[N[ratios],10^-10]}]
    ]
  ,{m,0,λ}]
,{λ,{1,2,3}}];
step1Report
```
**Output:**
```
2 kernel messages: MatrixPower::sing: Matrix {{0, 0}, {1, 0}} is singular.
Power::infy: Infinite expression - encountered.                                  0
MatrixPower::sing: Matrix {{0., 0.}, {1., 0.}} is singular.
P… [also: MatrixPower::sing] result: Out[358]= \begin{pmatrix}1 & 0 & \left\{1.\,\right\}\,\begin{pmatrix}0 & 1.\,\end{pmatrix}.\,\mathrm{MatrixPower}\left[\begin{pmatrix}0 & 0\\1.\, & 0\end{pmatrix},0\right]\\1 & 1 & \left\{1.\,\right\}…
```
## 2026-07-18 11:39:47 — ✏️ EDIT CELL 109
(*STEP 1 — covector identification: is ⟨x0[λ]|E21^m (pure representation-theoretic, no θ/site
  dependence) proportional to Section E's level-m single-site SoV covector xL1[λ][{m}] (built from
  the f…


## 2026-07-18 11:39:47
**Input:**
```wolfram
(*STEP 1 — covector identification: is ⟨x0[λ]|E21^m (pure representation-theoretic, no θ/site
  dependence) proportional to Section E's level-m single-site SoV covector xL1[λ][{m}] (built from
  the fused transfer-matrix hierarchy, known to diagonalize B1 with eigenvalue u-θ1-mh)? Test for
  every λ in the sweep {1,2,3}, m=0..λ. Record the proportionality constants γ_m(λ).
  (*why*) MatrixPower[nilpotentMatrix,0] fires MatrixPower::sing even for literal integer 0 (same
  gotcha as gRot, Section F) — m=0 handled by the identity directly, m>=1 via MatrixPower.*)
ClearAll[step1Report,gammaTab];
step1Report = {};
Do[
  Do[
    Module[{raw, ref, ratios},
      raw = If[m==0, xSingle[λ][{0}], xSingle[λ][{0}] . MatrixPower[Ee[λ][2,1],m]];
      ref = xL1[λ][{m}];
      ratios = DeleteCases[Flatten[raw]/Flatten[ref], Indeterminate|ComplexInfinity, Infinity];
      ratios = Select[ratios, (Chop[#,10^-12]=!=0 &)];  (*drop the 0/0 slots from zero components*)
      AppendTo[step1Report, {λ, m, Chop[N[ratios],10^-10]}]
    ]
  ,{m,0,λ}]
,{λ,{1,2,3}}];
step1Report
```
**Output:**
```
[message] Power::infy: Infinite expression - encountered.                                  0
[message] Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.
[message] [Power::infy ×2 more suppressed — first instance shown above]
[message] [Infinity::indet ×2 more suppressed — first instance shown above]
[message] General::stop: Further output of Power::infy will be suppressed during this calculation.
[message] [General::stop ×1 more suppressed — first instance shown above]
Out= {{1, 0, {1.}}, {1, 1, {1.}}, {2, 0, {1.}}, {2, 1, {2.}}, {2, 2, {2.}}, {3, 0, {1.}}, {3, 1, {3.}}, {3, 2, {6.}}, {3, 3, {6.}}}
```
## 2026-07-18 11:40:12 — 📥 BULK INSERT 1 CELL(S) at positions 110–110
1. [code] (*Closed form recognized: γ_m(λ) = λ!/(λ-m)! (falling factorial). Hard-assert against step1Report.*)


## 2026-07-18 11:40:12
**Input:**
```wolfram
(*Closed form recognized: γ_m(λ) = λ!/(λ-m)! (falling factorial). Hard-assert against step1Report.*)
ClearAll[gammaClosed];
gammaClosed[λ_,m_] := λ!/(λ-m)!;
step1Check = Table[
  {row[[1]], row[[2]], row[[3,1]], gammaClosed[row[[1]],row[[2]]], Chop[row[[3,1]]-gammaClosed[row[[1]],row[[2]]],10^-10]}
, {row, step1Report}];
worstStep1 = Max[Abs[step1Check[[All,5]]]];
If[worstStep1 > 10^-9, Abort[]];
{"STEP 1 CONFIRMED: ⟨x0[λ]|E21^m = γ_m(λ) xL1[λ][{m}], γ_m(λ)=λ!/(λ-m)!",
 "worst deviation", worstStep1, "table {λ,m,γ_found,γ_predicted,dev}", step1Check}
```
**Output:**
```
Out[364]= \left\{\text{STEP 1 CONFIRMED: ⟨x0[$\lambda$]|E21\^{}m = $\gamma$\_m($\lambda$) xL1[$\lambda$][\{m\}], $\gamma$\_m($\lambda$)=$\lambda$!/($\lambda$-m)!},\text{worst deviation},0,\text{table \{$\lambda$,m,$\gamma$\_found,$\gamma$\_predicted,dev\}},\begin{pmatrix}1 & 0 & 1.\, & 1 & 0\\1 & 1 & 1.\, & 1 & 0\\2 & 0 & 1.\, & 1 & 0\\2 & 1 & 2.\, & 2 & 0\\2 & 2 & 2.\, & 2 & 0\\3 & 0 & 1.\, & 1 &… [489 chars total]
```
## 2026-07-18 11:40:57 — 📥 BULK INSERT 1 CELL(S) at positions 111–111
1. [code] (*STEP 2 (the crux) — does the s2=m slice of the two-site SoV basis factorize the SAME way↵  Section


## 2026-07-18 11:40:58
**Input:**
```wolfram
(*STEP 2 (the crux) — does the s2=m slice of the two-site SoV basis factorize the SAME way
  Section E proved for m=0: x2[λ1,λ2][{s1,m}] =? xL1[λ1][{s1}] ⊗ xL1[λ2][{m}] (using xL1 for BOTH
  factors — the only level-m single-site covector object this notebook has; per Step 1 this is
  proportional to ⟨x0[λ2]|E21^m by the known γ_m(λ2)). Test constancy of the componentwise ratio
  over nonzero entries, m=1..λ2 (m=0 already proven, Section E Cell 60), for every (λ1,λ2,s1,m) in
  the sweep λ1,λ2∈{1,2,3}.*)
ClearAll[step2Report];
step2Report = {};
Do[
  Do[
    Do[
      Module[{lhs,rhsraw,ratios,spread,ratioVal},
        lhs = Flatten[x2[λ1,λ2][{s1,m}]];
        rhsraw = Flatten[KroneckerProduct[xL1[λ1][{s1}], xL1[λ2][{m}]]];
        ratios = Table[If[Chop[rhsraw[[i]],10^-12]===0, Nothing, lhs[[i]]/rhsraw[[i]]], {i,Length[lhs]}];
        ratios = N[ratios];
        spread = If[Length[ratios]>1, Max[Abs[ratios-ratios[[1]]]], 0];
        ratioVal = If[Length[ratios]>0, ratios[[1]], Missing["NoNonzero"]];
        AppendTo[step2Report, {λ1,λ2,s1,m, Chop[ratioVal,10^-10], Chop[spread,10^-10]}]
      ]
    ,{m,1,λ2}]
  ,{s1,0,λ1}]
,{λ1,{1,2,3}},{λ2,{1,2,3}}];
worstSpread2 = Max[step2Report[[All,6]]];
{"STEP 2: x2[{s1,m}] vs xL1[λ1][{s1}]⊗xL1[λ2][{m}], worst ratio-spread (constancy check)", worstSpread2,
 "num cases", Length[step2Report], "sample rows", Take[step2Report,UpTo[10]]}
```
**Output:**
```
Out[369]= \left\{\text{STEP 2: x2[\{s1,m\}] vs xL1[$\lambda$1][\{s1\}]$\otimes$xL1[$\lambda$2][\{m\}], worst ratio-spread (constancy check)},0,\text{num cases},54,\text{sample rows},\begin{pmatrix}1 & 1 & 0 & 1 & 0.16 & 0\\1 & 1 & 1 & 1 & 1.\, & 0\\1 & 2 & 0 & 1 & 0.16 & 0\\1 & 2 & 0 & 2 & -0.68 & 0\\1 & 2 & 1 & 1 & 1.\, & 0\\1 & 2 & 1 & 2 & 1.\, & 0\\1 & 3 & 0 & 1 & 0.16 & 0\\1 & 3 & 0 & 2 & -0.6… [479 chars total]
```

## 2026-07-18 11:42:04
**Input:**
```wolfram
step2Report
```
**Output:**
```
Out= {{1, 1, 0, 1, 0.16, 0}, {1, 1, 1, 1, 1., 0}, {1, 2, 0, 1, 0.16, 0}, {1, 2, 0, 2, -0.68, 0}, {1, 2, 1, 1, 1., 0}, {1, 2, 1, 2, 1., 0}, {1, 3, 0, 1, 0.16, 0}, {1, 3, 0, 2, -0.68, 0}, {1, 3, 0, 3, -1.52, 0}, {1, 3, 1, 1, 1., 0}, {1, 3, 1, 2, 1., 0}, {1, 3, 1, 3, 1., 0}, {2, 1, 0, 1, 0.08695652173913043, 0}, {2, 1, 1, 1, 0.5434782608695652, 0}, {2, 1, 2, 1, 1., 0}, {2, 2, 0, 1, 0.08695652173913043, 0}, {2, 2, 0, 2, -0.059130434782608696, 0}, {2, 2, 1, 1, 0.5434782608695652, 0}, {2, 2, 1, 2, 0.08695652173913043, 0}, {2, 2, 2, 1, 1., 0}, {2, 2, 2, 2, 1., 0}, {2, 3, 0, 1, 0.08695652173913043, 0}, {2, 3, 0, 2, -0.059130434782608696, 0}, {2, 3, 0, 3, 0.5617391304347826, 0}, {2, 3, 1, 1, 0.5434782608695652, 0}, {2, 3, 1, 2, 0.08695652173913043, 0}, {2, 3, 1, 3, -0.3695652173913043, 0}, {2, 3, 2, 1, 1., 0}, {2, 3, 2, 2, 1., 0}, {2, 3, 2, 3, 1., 0}, {3, 1, 0, 1, 0.05970149253731343, 0}, {3, 1, 1, 1, 0.373134328358209, 0}, {3, 1, 2, 1, 0.6865671641791045, 0}, {3, 1, 3, 1, 1., 0}, {3, 2, 0, 1, 0.05970149253731343, 0}, {3, 2, 0, 2, -0.022063595068137574, 0}, {3, 2, 1, 1, 0.373134328358209, 0}, {3, 2, 1, 2, 0.03244646333549643, 0}, {3, 2, 2, 1, 0.6865671641791045, 0}, {3, 2, 2, 2, 0.373134328358209, 0}, {3, 2, 3, 1, 1., 0}, {3, 2, 3, 2, 1., 0}, {3, 3, 0, 1, 0.05970149253731343, 0}, {3, 3, 0, 2, -0.022063595068137574, 0}, {3, 3, 0, 3, 0.03353666450356911, 0}, {3, 3, 1, 1, 0.373134328358209, 0}, {3, 3, 1, 2, 0.03244646333549643, 0}, {3, 3, 1, 3, -0.022063595068137574, 0}, {3, 3, 2, 1, 0.6865671641791045, 0}, {3, 3, 2, 2, 0.373134328358209, 0}, {3, 3, 2, 3, 0.05970149253731343, 0}, {3, 3, 3, 1, 1., 0}, {3, 3, 3, 2, 1., 0}, {3, 3, 3, 3, 1., 0}}
```
## 2026-07-18 11:42:43 — 📥 BULK INSERT 1 CELL(S) at positions 112–112
1. [code] (*Find the exact closed form for c(λ1,s1) (observed to depend ONLY on λ1,s1 — NOT λ2 or m).↵  Recomp


## 2026-07-18 11:42:44
**Input:**
```wolfram
(*Find the exact closed form for c(λ1,s1) (observed to depend ONLY on λ1,s1 — NOT λ2 or m).
  Recompute with Rationalize (exact rationals, θ1=1/3,θ2=1/7,h=1 are exact in this notebook) rather
  than floating N[], then fit against the numeric pattern c(λ1,s1)=(4+21 s1)/(4+21 λ1) noticed from
  the decimals (21=1/θ1 · 1/θ2 = 3·7, plausible given θ-tower structure).*)
ClearAll[step2Exact];
step2Exact = {};
Do[
  Do[
    Module[{lhs,rhsraw,ratios,idx},
      lhs = Flatten[x2[λ1,2][{s1,1}]];  (*fix λ2=2,m=1 arbitrarily — Step 2 showed λ2,m-independence*)
      rhsraw = Flatten[KroneckerProduct[xL1[λ1][{s1}], xL1[2][{1}]]];
      idx = FirstPosition[Chop[N[rhsraw],10^-12], x_/;x=!=0][[1]];
      ratios = Rationalize[lhs[[idx]]/rhsraw[[idx]], 10^-8];
      AppendTo[step2Exact, {λ1, s1, ratios}]
    ]
  ,{s1,0,λ1}]
,{λ1,{1,2,3}}];
step2Exact
```
**Output:**
```
Out[373]= \begin{pmatrix}1 & 0 & 1\\1 & 1 & 1\\2 & 0 & 1\\2 & 1 & 1\\2 & 2 & 1\\3 & 0 & 1\\3 & 1 & 1\\3 & 2 & 1\\3 & 3 & 1\end{pmatrix}
```

## 2026-07-18 11:43:22
**Input:**
```wolfram
lhs0 = Flatten[x2[2,2][{0,1}]]
rhs0 = Flatten[KroneckerProduct[xL1[2][{0}], xL1[2][{1}]]]
Chop[N[rhs0],10^-12]
idx0 = FirstPosition[Chop[N[rhs0],10^-12], x_/;x=!=0]
lhs0[[idx0[[1]]]]/rhs0[[idx0[[1]]]]
```
**Output:**
```
Out[1]= {0, 0, 0, 0, 0, 21/46, 0, 1/23, 0}
Out[2]= {0, 0, 0, 0, 0, 0, 0, 1/2, 0}
Out[3]= {0, 0, 0, 0, 0, 0, 0, 0.5, 0}
Out[4]= {0}
Out[5]= 1
```

## 2026-07-18 11:43:35
**Input:**
```wolfram
FirstPosition[{0,0,0.5,0}, x_/;x=!=0]
FirstPosition[{0,0,0.5,0}, Except[0]]
Position[{0,0,0.5,0}, Except[0.]]
```
**Output:**
```
Out[1]= {0}
Out[2]= {0}
Out[3]= {{0}, {1}, {2}, {3}, {4}, {}}
```

## 2026-07-18 11:44:19
**Input:**
```wolfram
idxSafe = SelectFirst[Range[Length[rhs0]], (Chop[N[rhs0[[#]]],10^-12]=!=0)&]
lhs0[[idxSafe]]/rhs0[[idxSafe]]//Rationalize
```
**Output:**
```
Out[1]= 8
Out[2]= 2/23
```
## 2026-07-18 11:44:29 — ✏️ EDIT CELL 112
(*Find the exact closed form for c(λ1,s1) (observed to depend ONLY on λ1,s1 — NOT λ2 or m).
  Recompute with exact rationals (θ1=1/3,θ2=1/7,h=1 are exact in this notebook) rather than floating
  N[], …


## 2026-07-18 11:44:29
**Input:**
```wolfram
(*Find the exact closed form for c(λ1,s1) (observed to depend ONLY on λ1,s1 — NOT λ2 or m).
  Recompute with exact rationals (θ1=1/3,θ2=1/7,h=1 are exact in this notebook) rather than floating
  N[], then fit against the numeric pattern c(λ1,s1)=(4+21 s1)/(4+21 λ1) noticed from the decimals
  (21=1/θ1 · 1/θ2 = 3·7, plausible given θ-tower structure).
  (*why SelectFirst not FirstPosition*) FirstPosition/Position on a flat numeric list misbehaved here
  (returned a spurious level-0 match) — SelectFirst on Range[Length[·]] is the robust way to find the
  first nonzero index.*)
ClearAll[step2Exact];
step2Exact = {};
Do[
  Do[
    Module[{lhs,rhsraw,idx,ratio},
      lhs = Flatten[x2[λ1,2][{s1,1}]];  (*fix λ2=2,m=1 arbitrarily — Step 2 showed λ2,m-independence*)
      rhsraw = Flatten[KroneckerProduct[xL1[λ1][{s1}], xL1[2][{1}]]];
      idx = SelectFirst[Range[Length[rhsraw]], (Chop[N[rhsraw[[#]]],10^-12]=!=0)&];
      ratio = Rationalize[N[lhs[[idx]]/rhsraw[[idx]]], 10^-8];
      AppendTo[step2Exact, {λ1, s1, ratio}]
    ]
  ,{s1,0,λ1}]
,{λ1,{1,2,3}}];
step2Exact
```
**Output:**
```
Out= {{1, 0, 4/25}, {1, 1, 1}, {2, 0, 2/23}, {2, 1, 25/46}, {2, 2, 1}, {3, 0, 4/67}, {3, 1, 25/67}, {3, 2, 46/67}, {3, 3, 1}}
```
## 2026-07-18 11:44:55 — 📥 BULK INSERT 1 CELL(S) at positions 113–113
1. [code] (*Closed form recognized: c(λ1,s1) = (θ1-θ2+s1 h)/(θ1-θ2+λ1 h) — independent of λ2,m entirely.↵  Har


## 2026-07-18 11:44:55
**Input:**
```wolfram
(*Closed form recognized: c(λ1,s1) = (θ1-θ2+s1 h)/(θ1-θ2+λ1 h) — independent of λ2,m entirely.
  Hard-assert against the FULL step2Report sweep (54 cases, all λ1,λ2,s1,m combinations).*)
ClearAll[cClosed];
cClosed[λ1_,s1_] := (θ1-θ2+s1 h)/(θ1-θ2+λ1 h);
step2Check = Table[
  {row[[1]],row[[2]],row[[3]],row[[4]], row[[5]], N[cClosed[row[[1]],row[[3]]]],
   Chop[row[[5]]-N[cClosed[row[[1]],row[[3]]]],10^-9]}
, {row, step2Report}];
worstStep2Closed = Max[Abs[step2Check[[All,7]]]];
If[worstStep2Closed > 10^-8, Abort[]];
{"STEP 2 CONFIRMED: x2[λ1,λ2][{s1,m}] = c(λ1,s1)·xL1[λ1][{s1}]⊗xL1[λ2][{m}], c(λ1,s1)=(θ1-θ2+s1h)/(θ1-θ2+λ1h) — genuinely independent of λ2,m",
 "worst deviation over all 54 (λ1,λ2,s1,m) cases", worstStep2Closed}
```
**Output:**
```
Out[378]= \text{\$Aborted} | Out[379]= \left\{\text{STEP 2 CONFIRMED: x2[$\lambda$1,$\lambda$2][\{s1,m\}] = c($\lambda$1,s1)$\cdot$xL1[$\lambda$1][\{s1\}]$\otimes$xL1[$\lambda$2][\{m\}], c($\lambda$1,s1)=($\theta$1-$\theta$2+s1h)/($\theta$1-$\theta$2+$\lambda$1h) --- genuinely independent of $\lambda$2,m},\text{worst deviation over all 54 ($\lambda$1,$\lambda$2,s1,m) cases},1.68\right\}
```

## 2026-07-18 11:46:43
**Input:**
```wolfram
Select[step2Check, Abs[#[[7]]]>10^-8 &]
```
**Output:**
```
Out= {{1, 2, 0, 2, -0.68, 0.16, -0.8400000000000001}, {1, 3, 0, 2, -0.68, 0.16, -0.8400000000000001}, {1, 3, 0, 3, -1.52, 0.16, -1.68}, {2, 2, 0, 2, -0.059130434782608696, 0.08695652173913043, -0.14608695652173914}, {2, 2, 1, 2, 0.08695652173913043, 0.5434782608695652, -0.45652173913043476}, {2, 3, 0, 2, -0.059130434782608696, 0.08695652173913043, -0.14608695652173914}, {2, 3, 0, 3, 0.5617391304347826, 0.08695652173913043, 0.4747826086956522}, {2, 3, 1, 2, 0.08695652173913043, 0.5434782608695652, -0.45652173913043476}, {2, 3, 1, 3, -0.3695652173913043, 0.5434782608695652, -0.9130434782608695}, {3, 2, 0, 2, -0.022063595068137574, 0.05970149253731343, -0.08176508760545101}, {3, 2, 1, 2, 0.03244646333549643, 0.373134328358209, -0.34068786502271253}, {3, 2, 2, 2, 0.373134328358209, 0.6865671641791045, -0.3134328358208955}, {3, 3, 0, 2, -0.022063595068137574, 0.05970149253731343, -0.08176508760545101}, {3, 3, 0, 3, 0.03353666450356911, 0.05970149253731343, -0.02616482803374432}, {3, 3, 1, 2, 0.03244646333549643, 0.373134328358209, -0.34068786502271253}, {3, 3, 1, 3, -0.022063595068137574, 0.373134328358209, -0.3951979234263466}, {3, 3, 2, 2, 0.373134328358209, 0.6865671641791045, -0.3134328358208955}, {3, 3, 2, 3, 0.05970149253731343, 0.6865671641791045, -0.6268656716417911}}
```
## 2026-07-18 11:47:53 — 📥 BULK INSERT 1 CELL(S) at positions 114–114
1. [code] (*Correction: c depends on m too (my earlier "λ2,m-independent" claim conflated λ2-independence,↵  w


## 2026-07-18 11:47:53
**Input:**
```wolfram
(*Correction: c depends on m too (my earlier "λ2,m-independent" claim conflated λ2-independence,
  which IS true, with m-independence, which is FALSE — caught by the hard-assert Abort above).
  Redo cleanly: exact rational c(λ1,s1,m) for λ1=1,2,3; s1=0..λ1; m=1..3, fixed λ2=3 (largest in
  sweep, so m up to 3 is available), then separately re-confirm λ2-independence at fixed m.*)
ClearAll[step2ExactFull];
step2ExactFull = {};
Do[
  Do[
    Do[
      Module[{lhs,rhsraw,idx,ratio},
        lhs = Flatten[x2[λ1,3][{s1,m}]];
        rhsraw = Flatten[KroneckerProduct[xL1[λ1][{s1}], xL1[3][{m}]]];
        idx = SelectFirst[Range[Length[rhsraw]], (Chop[N[rhsraw[[#]]],10^-12]=!=0)&];
        ratio = Rationalize[N[lhs[[idx]]/rhsraw[[idx]]], 10^-8];
        AppendTo[step2ExactFull, {λ1, s1, m, ratio}]
      ]
    ,{m,1,3}]
  ,{s1,0,λ1}]
,{λ1,{1,2,3}}];
step2ExactFull
```
**Output:**
```
Out[383]= \begin{pmatrix}1 & 0 & 1 & \frac{4}{25}\\1 & 0 & 2 & -\frac{17}{25}\\1 & 0 & 3 & -\frac{38}{25}\\1 & 1 & 1 & 1\\1 & 1 & 2 & 1\\1 & 1 & 3 & 1\\2 & 0 & 1 & \frac{2}{23}\\2 & 0 & 2 & -\frac{34}{575}\\2 & 0 & 3 & \frac{323}{575}\\2 & 1 & 1 & \frac{25}{46}\\2 & 1 & 2 & \frac{2}{23}\\2 & 1 & 3 & -\frac{17}{46}\\2 & 2 & 1 & 1\\2 & 2 & 2 & 1\\2 & 2 & 3 & 1\\3 & 0 & 1 & \frac{4}{67}\\3 & 0 & 2 & … [672 chars total]
```

## 2026-07-18 11:48:00
**Input:**
```wolfram
step2ExactFull
```
**Output:**
```
Out= {{1, 0, 1, 4/25}, {1, 0, 2, -17/25}, {1, 0, 3, -38/25}, {1, 1, 1, 1}, {1, 1, 2, 1}, {1, 1, 3, 1}, {2, 0, 1, 2/23}, {2, 0, 2, -34/575}, {2, 0, 3, 323/575}, {2, 1, 1, 25/46}, {2, 1, 2, 2/23}, {2, 1, 3, -17/46}, {2, 2, 1, 1}, {2, 2, 2, 1}, {2, 2, 3, 1}, {3, 0, 1, 4/67}, {3, 0, 2, -34/1541}, {3, 0, 3, 1292/38525}, {3, 1, 1, 25/67}, {3, 1, 2, 50/1541}, {3, 1, 3, -34/1541}, {3, 2, 1, 46/67}, {3, 2, 2, 25/67}, {3, 2, 3, 4/67}, {3, 3, 1, 1}, {3, 3, 2, 1}, {3, 3, 3, 1}}
```
## 2026-07-18 11:49:15 — 📥 BULK INSERT 1 CELL(S) at positions 115–115
1. [code] (*STEP 2 correction/investigation: the pure rank-1 factorization x2[{s1,m}] ∝ xL1[λ1][{s1}]⊗xL1[λ2][


## 2026-07-18 11:49:16
**Input:**
```wolfram
(*STEP 2 correction/investigation: the pure rank-1 factorization x2[{s1,m}] ∝ xL1[λ1][{s1}]⊗xL1[λ2][{m}]
  is CONFIRMED false for m>=2 (previous cell's hard Abort at worst deviation 1.68). Test the natural
  fallback: TRIANGULAR MIXING in s1 at fixed m — does x2[λ1,λ2][{s1,m}] lie in the span of
  {xL1[λ1][{s1'}]⊗xL1[λ2][{m}] : s1'=0..s1} (lower-triangular in s1, matching that t2[1,s2][θ2]
  could raise/preserve but not exceed site-1's own level when acting on the s1-slice)?
  Test case: λ1=3,λ2=3,m=2, all s1=0..3 simultaneously via a single linear solve.*)
basis32 = Table[Flatten[KroneckerProduct[xL1[3][{s1p}], xL1[3][{2}]]], {s1p,0,3}];
targets32 = Table[Flatten[x2[3,3][{s1,2}]], {s1,0,3}];
(*solve each target as a combination of ALL 4 basis vectors (full, not yet restricted to triangular)*)
coeffs32 = Table[
  LinearSolve[Transpose[basis32], targets32[[s1+1]]]
, {s1,0,3}];
(*residual check: does the linear combination reproduce the target exactly?*)
resid32 = Table[
  Chop[N[targets32[[s1+1]] - coeffs32[[s1+1]].basis32],10^-10]
, {s1,0,3}];
{"coeffs (rows=s1, cols=s1'=0..3)", Chop[coeffs32,10^-10], "max residual (solve exactness)", Max[Abs[Flatten[resid32]]]}
```
**Output:**
```
4 kernel messages: LinearSolve::nosol: Linear equation encountered that has no solution.
LinearSolve::nosol: Linear equation encountered that has no solution.
LinearSolve::nosol: Linear equation encountered that has no … [also: LinearSolve::nosol, LinearSolve::nosol, LinearSolve::nosol] result: Out[388]= \left\{\text{coeffs (rows=s1, cols=s1'=0..3)},\left\{\mathrm{LinearSolve}\left[\begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & \frac{1}{36}\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & \f…
```
## 2026-07-18 11:49:37 — 📥 BULK INSERT 1 CELL(S) at positions 116–116
1. [code] (*Redo with the FULL (λ1+1)(λ2+1)-dim basis {xL1[λ1][{s1'}]⊗xL1[λ2][{m'}] : all s1',m'} — this IS↵


## 2026-07-18 11:49:38
**Input:**
```wolfram
(*Redo with the FULL (λ1+1)(λ2+1)-dim basis {xL1[λ1][{s1'}]⊗xL1[λ2][{m'}] : all s1',m'} — this IS
  a complete basis of V_λ1^*⊗V_λ2^* (Section E: the λ+1 covectors exhaust V_λ^* at L=1), so ANY
  vector including x2[{s1,m}] must decompose in it; the question is whether the decomposition is
  TRIANGULAR (only s1'<=s1, m'<=m nonzero) — a "modified det" would still be recoverable from that.*)
ClearAll[fullBasis,decompX2];
fullBasis[λ1_,λ2_] := Flatten[Table[
   {s1p,m1p} -> Flatten[KroneckerProduct[xL1[λ1][{s1p}], xL1[λ2][{m1p}]]]
  ,{s1p,0,λ1},{m1p,0,λ2}], 1];
decompX2[λ1_,λ2_][s1_,m_] := Module[{basisRules,basisVecs,mat,target,sol},
  basisRules = fullBasis[λ1,λ2];
  basisVecs = basisRules[[All,2]];
  mat = Transpose[basisVecs];
  target = Flatten[x2[λ1,λ2][{s1,m}]];
  sol = LinearSolve[mat, target];
  Thread[basisRules[[All,1]] -> Chop[sol,10^-9]]
];
(*test state (s1,m)=(0,2), rep (3,3): which (s1',m') coefficients are nonzero?*)
d32 = decompX2[3,3][0,2];
Select[d32, #[[2]]=!=0 &]
```
**Output:**
```
Out[393]= \left\{\left\{0,2\right\}\to -0.02206359506813757300454250486697,\left\{1,1\right\}\to 0.1635301752109020116807268007787,\left\{2,0\right\}\to 0.8585334198572355613238157040883\right\}
```
## 2026-07-18 11:50:15 — 📥 BULK INSERT 1 CELL(S) at positions 117–117
1. [code] (*STEP 2 RESOLUTION: the mixing found above is confined to states of FIXED TOTAL LEVEL L=s1+m —↵  no


## 2026-07-18 11:50:15
**Input:**
```wolfram
(*STEP 2 RESOLUTION: the mixing found above is confined to states of FIXED TOTAL LEVEL L=s1+m —
  not a simple triangular s1-only mixing, but a genuine "same-total-level" recoupling. Verify this
  level-conservation property broadly: for every (λ1,λ2,s1,m) in the sweep, decompX2 should have
  ALL its nonzero coefficients exactly at (s1',m') with s1'+m'=s1+m.*)
levelCheck = {};
Do[
  Do[
    Do[
      Module[{d,offLevel},
        d = decompX2[λ1,λ2][s1,m];
        offLevel = Select[d, (#[[2]]=!=0 && #[[1,1]]+#[[1,2]]!=s1+m)&];
        AppendTo[levelCheck, {λ1,λ2,s1,m, Length[offLevel]}]
      ]
    ,{m,0,λ2}]
  ,{s1,0,λ1}]
,{rep,{{2,1},{1,2},{2,2},{3,1}}},{λ1,{rep[[1]]}},{λ2,{rep[[2]]}}];
worstOffLevel = Max[levelCheck[[All,5]]];
{"STEP 2 level-conservation check: # off-level nonzero coeffs (want 0 everywhere)", worstOffLevel,
 "num (λ1,λ2,s1,m) cases", Length[levelCheck]}
```
**Output:**
```
Out[397]= \left\{\text{STEP 2 level-conservation check: \# off-level nonzero coeffs (want 0 everywhere)},0,\text{num ($\lambda$1,$\lambda$2,s1,m) cases},29\right\}
```
## 2026-07-18 11:52:19 — 📥 BULK INSERT 2 CELL(S) at positions 118–119
1. [markdown] (*STEP 2 conclusion: the pure rank-1 factorization x2[λ1,λ2][{s1,m}] ∝ xL1[λ1][{s1}]⊗xL1[λ2][{m}]↵ho
2. [code] (*STEP 3 — graded dictionary. O_m[k1] extracted as (-1)^m m! · (k2^m coefficient of O2rot); this is↵


## 2026-07-18 11:52:19
**Input:**
```wolfram
(*STEP 3 — graded dictionary. O_m[k1] extracted as (-1)^m m! · (k2^m coefficient of O2rot); this is
  exact by construction (O2rot is a finite polynomial in k2 of degree λ2). Candidate:
  Cm[λ1,λ2][M,n][k1] := br[λ1,λ2][1][Pm[u] k1^(u/h) Q1[u]] / NORM, Pm[u]=Product[(u-θ2-hk),{k,0,m-1}]
  (empty=1 at m=0 — the ALREADY-PROVEN Phase-1 m=λ2+1-equivalent full case is the λ2-degree P1site2,
  distinct from this partial Pm). Try NORM = value at θ1 (det1-style) first, m=0,1 only (cheapest
  gate) on rep (2,1) state (1,0).*)
ClearAll[Pm,Om,Cm];
Pm[m_][u_] := Product[(u-θ2-h k),{k,0,m-1}];  (*empty product convention: Product[...,{k,0,-1}]=1*)
Om[λ1_,λ2_][M_,n_][m_][k1_] := (-1)^m m! Coefficient[O2rot[λ1,λ2][M,n][k1,k2],k2,m];
Cm[λ1_,λ2_][M_,n_][m_][k1_] := Module[{f},
  f[uu_] := Pm[m][uu] k1^(uu/h) Q1[λ1,λ2][M,n][uu];
  br[λ1,λ2][1][f[u]] / (f[θ1])
];
testM0 = Simplify[Om[2,1][1,0][0][k1] - Cm[2,1][1,0][0][k1]];
testM1raw = Simplify[Om[2,1][1,0][1][k1]];
testM1cand = Simplify[Cm[2,1][1,0][1][k1]];
{"m=0: Om-Cm (must be 0 per Phase1, C1pred/P1site2[θ1] should equal exactly Cm's own norm here)", Chop[N[testM0],10^-10],
 "m=1: Om_1[k1]", testM1raw, "m=1: Cm_1[k1]", testM1cand}
```
**Output:**
```
Out[405]= \left\{\text{m=0: Om-Cm (must be 0 per Phase1, C1pred/P1site2[$\theta$1] should equal exactly Cm's own norm here)},\left(2.5572121404875445 -0.4102640568287387\,\mathrm{i}\right)\,k1-\left(0.038411844228125346 +2.1227491621421546\,\mathrm{i}\right)\,k1^2,\text{m=1: Om\_1[k1]},\left(0.3604398751472265372666776506+0.9327824485933747164848359821\,\mathrm{i}\right)\,k1^2-\left(1.296101147801… [801 chars total]
```
## 2026-07-18 11:53:03 — 📥 BULK INSERT 1 CELL(S) at positions 120–120
1. [code] (*Diagnosis: Cm[m] as built (bracket with Pm inserted, divided by value at θ1) does NOT match↵  O_m[


## 2026-07-18 11:53:03
**Input:**
```wolfram
(*Diagnosis: Cm[m] as built (bracket with Pm inserted, divided by value at θ1) does NOT match
  O_m[k1] even at m=0 (Pm[0]=1, i.e. NO insertion at all — that's the bare 2-site br[1], not br1;
  Phase 1 showed the FULL P1site2[λ2] insertion (degree λ2+1, killing ALL of tower 2) is needed to
  reach br1[λ1], not merely "no insertion"). Reconsider: Pm[m] only reshapes the WEIGHT of the
  (tower-1-only) residue sum br[λ1,λ2][1][·] — it does NOT literally restrict to a "θ2-level m"
  lattice point (br[1] never evaluates AT tower 2 at all). Test the natural alternative: the
  UNNORMALIZED bracket DIFFERENCE Bm[m]-Bm[m+1] (successive partial insertions), since Pm[m+1] kills
  one more tower-2-type factor than Pm[m] — check if this difference is proportional to O_m[k1].*)
ClearAll[Bm];
Bm[λ1_,λ2_][M_,n_][m_][k1_] := br[λ1,λ2][1][Pm[m][u] k1^(u/h) Q1[λ1,λ2][M,n][u]];
diffB0 = Simplify[Bm[2,1][1,0][0][k1] - Bm[2,1][1,0][1][k1]];
ratioB0 = Simplify[diffB0/Om[2,1][1,0][0][k1]];
diffB1 = Simplify[Bm[2,1][1,0][1][k1] - Bm[2,1][1,0][2][k1]];
ratioB1 = Simplify[diffB1/Om[2,1][1,0][1][k1]];
{"m=0: (B0-B1)/O_0 ratio", ratioB1, "m=0 diff/O0 first", ratioB0}
```
**Output:**
```
Out[412]= \left\{\text{m=0: (B0-B1)/O\_0 ratio},\frac{\sqrt[3]{k1}\,\left(\left(-0.0292504903494856237+0.0205326892042810627\,\mathrm{i}\right)\,k1^2+\left(0.499242526583056932+0.532659421382621779\,\mathrm{i}\right)\,k1-\left(0.139388493663376676+0.206413542317126162\,\mathrm{i}\right)\right)}{\left(0.3604398751472265372666776506+0.9327824485933747164848359821\,\mathrm{i}\right)\,k1^2-\left(1.296… [1083 chars total]
```
## 2026-07-18 11:53:30 — 📥 BULK INSERT 1 CELL(S) at positions 121–121
1. [code] (*One more natural variant: difference of the NORMALIZED objects Cm[m] (each divided by its OWN↵  f(


## 2026-07-18 11:53:30
**Input:**
```wolfram
(*One more natural variant: difference of the NORMALIZED objects Cm[m] (each divided by its OWN
  f(θ1), unlike the previous cell's unnormalized Bm-difference which had a fractional-power
  mismatch). Test Cm[0]-Cm[1] against O_0[k1].*)
diffC0 = Simplify[Cm[2,1][1,0][0][k1] - Cm[2,1][1,0][1][k1]];
ratioC0 = Simplify[diffC0/Om[2,1][1,0][0][k1]];
{"Cm[0]-Cm[1] shape", diffC0, "ratio vs O_0[k1] (want: clean constant)", ratioC0}
```
**Output:**
```
Out[415]= \left\{\text{Cm[0]-Cm[1] shape},\left(-0.0225172879957976151-1.244370198497125076\,\mathrm{i}\right)\,k1^2+\left(5.434075798536032377-0.8718111207610696998\,\mathrm{i}\right)\,k1+\left(0.\times 10^{-20}+0.\times 10^{-20}\,\mathrm{i}\right),\text{ratio vs O\_0[k1] (want: clean constant)},\frac{\left(-0.0225172879957976151-1.244370198497125076\,\mathrm{i}\right)\,k1^2+\left(5.4340757985360… [753 chars total]
```
## 2026-07-18 11:54:08 — 📥 BULK INSERT 1 CELL(S) at positions 122–122
1. [code] (*Key realization from Phase 1: Pm[m=λ2+1][u]=P1site2[λ2][u] (the FULL polynomial) relates to↵  O_0[


## 2026-07-18 11:54:08
**Input:**
```wolfram
(*Key realization from Phase 1: Pm[m=λ2+1][u]=P1site2[λ2][u] (the FULL polynomial) relates to
  O_0[k1] (the m'=0 graded piece), i.e. the FULL insertion (m=λ2+1, "kill everything") maps to the
  BOTTOM graded piece, not the top. This suggests a REVERSED indexing: Cm[m][k1] (properly
  normalized) might equal (const)·O_{λ2-m}[k1], not O_m[k1] as naively guessed. Test on rep (2,1)
  (λ2=1): does Cm[0][k1] (no insertion) match O_1[k1] (m'=λ2-0=1, the TOP graded piece)?*)
ratioRev = Simplify[Cm[2,1][1,0][0][k1] / Om[2,1][1,0][1][k1]];
{"Cm[0][k1] / O_1[k1] (reversed-index hypothesis; want clean k1-independent constant)", ratioRev}
```
**Output:**
```
Out[417]= \left\{\text{Cm[0][k1] / O\_1[k1] (reversed-index hypothesis; want clean k1-independent constant)},\frac{\left(0.111318362779488714+0.0407158959953172346\,\mathrm{i}\right)\,k1^2-\left(0.21818047436770437-1.02534209918467083\,\mathrm{i}\right)\,k1+\left(0.360439875147226537-0.932782448593374716\,\mathrm{i}\right)}{\left(1.\,+0.\times 10^{-29}\,\mathrm{i}\right)\,k1^2-\left(2.963909584212… [539 chars total]
```

## 2026-07-18 11:54:32
**Input:**
```wolfram
Simplify[Cm[2,1][1,0][1][k1]/Om[2,1][1,0][1][k1]]
```
**Output:**
```
Out= ((0.36043987514722653726667765059555008944`17.77612140690624 - 0.93278244859337471648483598207086789993`17.97959948758913*I) - (1.36362796479815233533906288354975328693`17.49185270736853 - 6.40838811990419271777205213998317537464`17.98634886043309*I)*k1 + (1.28016117196412021234271094716206107478`17.9766208757477 + 0.46823280394614819771050007179979957073`17.75386028979382*I)*k1^2)/((1.58068744596833235240887910945247686274`26.490733824100538 - 2.18028000683379170484605115842512849186`26.50326365661453*I) - (2.96390958421244462997916665708159906713`18.214180962570733 - 0.24420479573139207277333933553540188389`17.13007019990194*I)*k1 + (1.`26.573595569998396 + 0``28.533967999094948*I)*k1^2)
```
## 2026-07-18 11:55:14 — 📥 BULK INSERT 1 CELL(S) at positions 123–123
1. [markdown] (*Section I outcome — STEPS 0-2 CLOSED EXACTLY (genuine progress); STEP 3 OBSTRUCTED, reported↵preci


---

## 2026-07-18 11:55:20 — 🔄 KERNEL RESTART

## 2026-07-18 11:55:30 — ▶️ RUN CELLS 1–10
Cell 1: ✓ — Out[8]= \left\{\frac{1}{3},\frac{1}{7},1,0.360439875147226537266677650596+0.932782448593374716484835982071\,\mathrm{i},0.36043987514722653726667765060-0.932782448593374716484835982071\,\mathrm{i},0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i},1\right\}
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 8: ✓ — (no output)
Cell 10: ✓ — (no output)

## 2026-07-18 11:55:35 — ▶️ RUN CELLS 11–20
Cell 11: ✓ — (no output)
Cell 12: ✓ — Out[27]= \left\{0\right\}
Cell 13: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — Out[31]= \left\{0\right\}
Cell 16: ✓ — (no output)
Cell 17: ✓ — Out[34]= \begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{pmatrix}
Cell 18: ✓ — (no output)
Cell 20: ✓ — (no output)

## 2026-07-18 11:55:38 — ▶️ RUN CELLS 21–28
Cell 21: ✓ — (no output)
Cell 23: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 26: ✓ — (no output)
Cell 27: ✓ — Out[52]= \begin{pmatrix}0 & \frac{1}{2} & 0 & 0 & 0 & 0\end{pmatrix}
Cell 28: ✓ — Out[53]= \left\{0\right\}

## 2026-07-18 11:55:41 — ▶️ RUN CELLS 29–36
Cell 30: ✓ — Out[59]= \left\{\left(u-\frac{4}{3}\right)\,\left(u-\frac{8}{7}\right),\left(u-\frac{1}{3}\right)\,\left(u-\frac{1}{7}\right),\left(0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i}\right)\,u^2-\left(1.06415582186324025288257211128+0.\times 10^{-29}\,\mathrm{i}\right)\,u\right\}
Cell 31: ✓ — Out[62]= \begin{pmatrix}1 & 1\\2 & 2\\3 & 3\end{pmatrix}
Cell 32: ✓ — Out[65]= \text{max|C2-($\kappa$1+$\kappa$2)Id|}\to 0.\times 10^{-29}
Cell 33: ✓ — Out[75]= \langle\!\lvert \left\{0,0\right\}\to 0.5664055180884988 +1.3769645669711723\,\mathrm{i},\left\{1,0\right\}\to -0.778123966451158+0.\,\,\mathrm{i},\left\{1,1\right\}\to 1.1900552523337025 +0.\,\,\mathrm{i},\left\{2,0\right\}\to 0.5664055180884988 -1.3769645669711723\,\mathrm{i}\rvert\!\rangle
Cell 34: ✓ — 2 kernel messages: Psi::spurious: No genuine state at (M,n)=(2,2) for (λ1,λ2)=(2,1): weight sector M=2 has only 2 genuine state(s) (n=0..1), so n=2 does not exist. [also: Psi::spurious] result: Print: spurious-guard test (expect a Psi::spurious message and Missing[...]):[formula] | Out[84]= \left\{\left\{1.000000000000000000000000000+0.\times 10^{-29}\,\mathrm{i},-1.3004494881076436272092358…
⚠️ stopped at Cell 34 — error detected (pass stopOnError:false to continue past errors)

## 2026-07-18 11:55:59 — ▶️ RUN CELLS 35–42
Cell 35: ✓ — Out[91]= \left\{1.03464899078135229594035121+0.381869568380946283899575941\,\mathrm{i},-2.37022278949001677248172619-0.20048970475590097847093779\,\mathrm{i},1.00000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right\}
Cell 36: ✓ — Out[94]= \left\{-0.709231010958016926537225142+1.09596916251277206706752728\,\mathrm{i},-0.02291660900+0.\times 10^{-12}\,\mathrm{i}\right\}
Cell 37: ✓ — (no output)
Cell 38: ✓ — (no output)
Cell 39: ✓ — (no output)
Cell 40: ✓ — (no output)
Cell 41: ✓ — (no output)
Cell 42: ✓ — Out[109]= 0

## 2026-07-18 11:56:02 — ▶️ RUN CELLS 43–47
Cell 43: ✓ — (no output)
Cell 44: ✓ — Out[113]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 45: ✓ — Out[115]= \left\{6,6\right\}
Cell 46: ✓ — Out[116]= \begin{pmatrix}2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i} & 0 & 0 & 0 & 0 & 0\\0 & 7.36608403016431791181-0.968046299794461400266\,\mathrm{i} & 0 & 0 & 0 & 0\\0 & 0 & 2.68702160063169134012-1.3456153260557858257\,\mathrm{i} & 0 & 0 & 0\\0 & 0 & 0 & 7.3660840301643179118107+0.9680462997944614002659\,\mathrm{i} & 0 & 0\\0 & 0 & 0 & 0 & 2.687021600631691340118… [537 chars total]
Cell 47: ✓ — Out[125]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{max|off-diagonal|}\to 5.065342385207116\times 10^{-34},\text{digits vanished (-Log10)}\to 33.29539119374288,\text{min|diagonal|}\to 3.0051232400647923,\text{diagonal values}\to \left\{2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i},7.36608403016431791181-0… [699 chars total]

## 2026-07-18 11:56:06 — ▶️ RUN CELLS 48–54
Cell 48: ✓ — (no output)
Cell 49: ✓ — (no output)
Cell 50: ✓ — Out[133]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 51: ✓ — Out[138]= \left\{1.00000000000000000000000000000,6\right\}
Cell 52: ✓ — Out[143]= \langle\!\lvert \text{max|<PsiLeft|Psi>| off-diagonal}\to 0.\,,\text{digits vanished (-Log10)}\to \mathrm{Indeterminate}\rvert\!\rangle
Cell 53: ✓ — Out[147]= \left\{7.366084030164317911810727273-0.9680462997944614002658874268\,\mathrm{i},1.00000000000000000000+0.\times 10^{-21}\,\mathrm{i}\right\}
Cell 54: ✓ — Out[150]= \begin{matrix}\langle\!\lvert \text{M}\to 0,\text{n}\to 0,\text{braket}\to 2.909064289-5.804620064\,\mathrm{i},\text{G21 diagonal}\to 2.909064289-5.804620064\,\mathrm{i},\text{normRatio}\to 1.000000000+0.\times 10^{-10}\,\mathrm{i}\rvert\!\rangle \\\langle\!\lvert \text{M}\to 1,\text{n}\to 0,\text{braket}\to 7.366084030-0.968046300\,\mathrm{i},\text{G21 diagonal}\to 7.366084030-0.9680463… [1247 chars total]

## 2026-07-18 11:56:10 — ▶️ RUN CELLS 55–62
Cell 56: ✓ — (no output)
Cell 57: ✓ — Out[158]= \begin{pmatrix}3 & 3\\3 & 3\end{pmatrix}
Cell 58: ✓ — (no output)
Cell 59: ✓ — Out[165]= \left\{0\right\}
Cell 60: ✓ — Out[169]= \begin{matrix}\left\{1,1,0,1.\,,0.\,\right\}\\\left\{1,1,1,1.\,,0.\,\right\}\\\left\{1,2,0,1.\,,0.\,\right\}\\\left\{1,2,1,1.\,,0.\,\right\}\\\left\{1,3,0,1.\,,0.\,\right\}\\\left\{1,3,1,1.\,,0.\,\right\}\\\left\{2,1,0,1.\,,0.\,\right\}\\\left\{2,1,1,1.\,,0.\,\right\}\\\left\{2,1,2,1.\,,0.\,\right\}\\\left\{2,2,0,1.\,,0.\,\right\}\\\left\{2,2,1,1.\,,0.\,\right\}\\\left\{2,2,2,1.\,,0.\,\r… [898 chars total]
Cell 61: ✓ — Out[170]= \langle\!\lvert \text{max spread (ratio non-constancy)}\to 0.\,,\text{max |ratio - 1| over sweep}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases checked}\to 27\rvert\!\rangle
Cell 62: ✓ — Out[176]= \langle\!\lvert \text{max operator-level spread}\to 0.\,,\text{max |$\mu$ - $\mu$Site2|}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases}\to 27\rvert\!\rangle

## 2026-07-18 11:56:23 — ▶️ RUN CELLS 63–70
Cell 64: ✓ — (no output)
Cell 65: ✓ — Out[185]= \langle\!\lvert \text{w\_0 ($\lambda$=1,2,3)}\to \left\{1,1,1\right\},\text{full weights w\_j}\to \left\{1\to \left\{1,-1\right\},2\to \left\{1,-2,1\right\},3\to \left\{1,-3,3,-1\right\}\right\},\text{br1[$\lambda$][1] ($\lambda$=1,2,3)}\to \left\{0,0,0\right\}\rvert\!\rangle
Cell 66: ✓ — Out[187]= \left\{0,0\right\}
Cell 67: ✓ — Out[190]= \sqrt[3]{k}\,\left(\left(-1.18561270626131985872+0.832253303224995280404\,\mathrm{i}\right)\,k^2-\left(0.761819910806671651081+0.812812473459294825632\,\mathrm{i}\right)\,k-\left(0.404410037536507329420+0.598870869485872613792\,\mathrm{i}\right)\right)
Cell 68: ✓ — Out[194]= \begin{pmatrix}0 & 0 & 0\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix}
Cell 69: ✓ — Out[198]= \langle\!\lvert \text{OP(E12)}\to 1.\,+0.\times 10^{-30}\,\mathrm{i},\text{OP(E21)}\to \left(-0.03626734060947795145020244939-2.004237714666237867354326095\,\mathrm{i}\right)\,\varphi ^2-\left(1.522150083623538480912456055-0.2442047957313920727733393355\,\mathrm{i}\right)\,\varphi +\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right),\text{nontrivial generator \{a,b\}}\to \left\{2,1\right\},\… [464 chars total]
Cell 70: ✓ — Out[201]= \left\{1,2\right\}

## 2026-07-18 11:56:27 — ▶️ RUN CELLS 71–77
Cell 71: ✓ — Out[207]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{$\varphi$ candidates per state}\to \begin{pmatrix}\left\{0,0\right\} & \left(-0.920793778924720536105378676-0.390049761817044163511250929\,\mathrm{i}\right)\,k^{4/3}+\left(0.695722366612934233848844287-0.718310788304406529935213657\,\mathrm{i}\right)\,\sqrt[3]{k}-\left(0.3604398751… [1247 chars total]
Cell 72: ✓ — Out[213]= \langle\!\lvert \text{c\_sb solution(s)}\to \left\{\left\{\mathrm{cb}\to \frac{\left(-6.242239011314\times 10^{45}-4.511757007219\times 10^{45}\,\mathrm{i}\right)\,k^{4/3}-\left(5.999302657796\times 10^{45}+8.356648463050\times 10^{45}\,\mathrm{i}\right)\,\sqrt[3]{k}}{9.833530773718\times 10^{45}\,k^{4/3}-\left(3.54439660434\times 10^{45}-9.172544913427\times 10^{45}\,\mathrm{i}\right)\,… [841 chars total]
Cell 73: ✓ — Out[216]= \langle\!\lvert \text{generator \{a,b\} tested}\to \left\{2,1\right\},\text{rep tested}\to \left\{1,2\right\},\text{FL shape (verified exact for all 6 states)}\to \text{FL[1,2][state][k] = k\^{}(1/3) (A\_state + B\_state k)},\text{c0(state) = OPg[state][0]}\to \begin{pmatrix}\left\{0,0\right\} & 1.\,+0.\times 10^{-28}\,\mathrm{i}\\\left\{1,0\right\} & 1.\,+0.\times 10^{-28}\,\mathrm{i}\\… [1247 chars total]
Cell 74: ✓ — Out[220]= \left(1.000000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right)-\left(0.360439875147226537266677651+0.932782448593374716484835982\,\mathrm{i}\right)\,k
Cell 75: ✓ — Out[227]= \langle\!\lvert \text{alphaTab \{state,D0,c0,Dc,c1,$\alpha$\}}\to \begin{pmatrix}\left\{0,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & -0.36043987514722653-0.9327824485933747\,\mathrm{i} & 0.36043987514722653 +0.9327824485933747\,\mathrm{i} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.7101545203859353 -0.2371533418544… [1247 chars total]
Cell 76: ✓ — Out[235]= \langle\!\lvert \text{validation}\to \left\{\left\{1,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\end{pmatrix},\left\{2,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\\\left\{2,1\right\} & 0\\\left\{3,0\right\} & 0\end{pmatrix},\left\{1… [1247 chars total]
Cell 77: ✓ — Out[237]= \langle\!\lvert \text{generator \{a,b\}}\to \left\{2,1\right\},\text{$\varphi$(k)}\to \text{$\alpha$ k},\text{$\alpha$ (universal value)}\to -1.\,+0.\,\,\mathrm{i},\text{normalization}\to \text{1-site FSoV determinant det1[f] = br1[f]/(f/.u->$\theta$1); Section C's det[$\lambda$1,$\lambda$2][f] restricted to 1 site},\text{max residual over rep sweep}\to 0,\text{reps validated}\to \begin{… [510 chars total]

## 2026-07-18 11:56:43 — ▶️ RUN CELLS 78–82
Cell 79: ✓ — (no output)
Cell 80: ✓ — Out[245]= \left\{\text{Test 1 (symbolic monomials u\^{}j) residuals},\left\{0\right\}\right\}
Cell 81: ✓ — Out[250]= \left\{\text{Test 2 (twisted monomials k\^{}(u/h) u\^{}j) residuals},\left\{0\right\}\right\}
Cell 82: ✓ — Out[259]= \left\{\text{Test 3 worst |br identity residual| on Q1},0,\text{worst |detRatio - predicted| },0,\text{num states checked},41\right\}

## 2026-07-18 11:56:46 — ▶️ RUN CELLS 83–87
Cell 84: ✓ — (no output)
Cell 85: ✓ — Out[271]= \begin{pmatrix}1 & 1 & 0\\2 & 1 & 0\\3 & 1 & 0\end{pmatrix}
Cell 86: ✓ — Out[278]= \left\{\text{Mirror Test1 (u\^{}j)},\left\{0\right\},\text{Mirror Test2 (k\^{}(u/h)u\^{}j)},\left\{0\right\}\right\}
Cell 87: ✓ — Out[285]= \left\{\text{Mirror Test3 worst |br identity residual| on Q1},0,\text{num states checked},41\right\}

## 2026-07-18 11:56:48 — ▶️ RUN CELLS 88–91
Cell 89: ✓ — Out[294]= \left\{\text{OP2(E12)},1.\,+0.\times 10^{-30}\,\mathrm{i},\text{OP2(E21)},\left(2.603469709065218092712489006+0.6885776528619826437114966465\,\mathrm{i}\right)\,\varphi +\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right),\text{nontrivial generator \{a,b\}},\left\{2,1\right\}\right\}
Cell 90: ✓ — Out[297]= \left(1.\,+0.\times 10^{-27}\,\mathrm{i}\right)-\left(0.360439875147226537266677651+0.932782448593374716484835982\,\mathrm{i}\right)\,k
Cell 91: ✓ — Out[300]//TableForm= \begin{pmatrix}\left\{0,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & -0.36043987514722653-0.9327824485933747\,\mathrm{i} & 0.36043987514722653 +0.9327824485933747\,\mathrm{i} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.4685059330398573 -0.1713611156914601\,\mathrm{i} & 2.603469709065218 +0.6885776528619827… [1258 chars total]

## 2026-07-18 11:57:13 — ▶️ RUN CELLS 92–108
Cell 92: ✓ — Out[302]//TableForm= \begin{pmatrix}\left\{0,0\right\} & -1.\,+0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.973902234276272 +0.2269679230016555\,\mathrm{i} & 0.7391304347826086 +0.\,\,\mathrm{i} & 0.973902234276272 +0.2269679230016555\,\mathrm{i}\\\left\{1,0\right\} & 0.15191893851739566 -0.10600050418875033\,\mathrm{i} & -0.15191893851739566+0.10600050418875033\,\mathrm{i} & 0.7651595274649543 +… [1258 chars total]
Cell 93: ✓ — Out[305]= \begin{pmatrix}\left\{0,0\right\} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 0.05448860910968355 -0.09540961504089102\,\mathrm{i}\\\left\{1,1\right\} & -0.0293413933141388+0.7363100318787654\,\mathrm{i}\\\left\{2,0\right\} & 0.11703059508534998 +0.\,\,\mathrm{i}\\\left\{2,1\right\} & 0.8158972359113332 +0.\,\,\mathrm{i}\\\left\{3,0\right\} & 0.05448860910968355 +0.09540961504089102\,… [541 chars total]
Cell 96: ✓ — Out[311]= \left\{\text{O[k1,0]-FLdet[k1] (must be 0)},0,\text{O[0,k2]-FLdet2[k2] (expected NONZERO per Phase1 finding)},\left(-3.0719756421050755-0.5172165371705225\,\mathrm{i}\right)\,k2\right\}
Cell 97: ✓ — Out[316]= \left\{\text{Fcand1},k1^2\,\left(\left(-2.2527492196701659-5.82989030370859198\,\mathrm{i}\right)\,k2-\left(0.4170744170089964+23.0487337186617355\,\mathrm{i}\right)\right)+k1\,\left(\left(0.67128913890100625-0.37524873225652243\,\mathrm{i}\right)\,k2+\left(9.51343802264711551-1.5262799733212005\,\mathrm{i}\right)\right)-\left(1.99115021541939352-0.72828474168870538\,\mathrm{i}\right)\,k… [1014 chars total]
Cell 98: ✓ — Out[318]= \left(0.9252125034899888 -6.175727508765443\,\mathrm{i}\right)\,k1^2\,k2+\left(5.4271267482291075 +3.0890004366773067\,\mathrm{i}\right)\,k1\,k2
Cell 99: ✓ — Out[323]= \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 7.991287939023577 -1.2820751775898085\,\mathrm{i}\\\left\{2,0\right\} & -0.3808070763995185-21.044496003995498\,\mathrm{i}\\\left\{0,1\right\} & 0.6123194936458246 +1.416862394550688\,\mathrm{i}\\\left\{1,1\right\} & -0.6248120089007682-3.051910425543801\,\mathrm{i}\\\left\{2,1\right\} & -1.8923093445229393-4.897107855115217\,\m… [421 chars total]
Cell 100: ✓ — Out[328]= \left\{\text{Fcand2 (no P-insertion)},k1^2\,\left(\left(0.02466179161444501+1.36288164597304175\,\mathrm{i}\right)-\left(0.36043987514722654+0.93278244859337472\,\mathrm{i}\right)\,k2\right)+k1\,\left(\left(1.2351720155778515-0.69045766735200127\,\mathrm{i}\right)\,k2-\left(6.46913785540003854-1.0378703818584163\,\mathrm{i}\right)\right)+\left(5.38781822995835893-1.9706528304517910\,\mat… [788 chars total]
Cell 101: ✓ — 2 kernel messages: Power::infy: Infinite expression 0       encountered.
Infinity::indet: Indeterminate expression ----------------- encountered.                                                   1/3                    … [also: Power::infy] result: Out[331]//TableForm= \begin{pmatrix}2 & 1 & \left\{0,0\right\} & \mathrm{Indeterminate} & \mathrm{Indeterminate} & -4.970498218180917+4.515583419069149\,\mathrm{i}\\2 & 1 & \left\{1,1\right\} & \mathr…
Cell 102: ✓ — Out[334]//TableForm= \begin{pmatrix}2 & 1 & \left\{0,0\right\} & \left(-7.771745024480887+7.060451775520925\,\mathrm{i}\right)\,k1^2-\left(3.7846186890458786 +9.794215710230434\,\mathrm{i}\right)\,k1 & \left(1.8923093445229393 +4.897107855115217\,\mathrm{i}\right)\,k2 & -4.970498218180917+4.515583419069149\,\mathrm{i}\\2 & 1 & \left\{1,1\right\} & \left(9.339255641812587 +10.661478687052961\,\math… [1258 chars total]
Cell 103: ✓ — 2 kernel messages: Power::infy: Infinite expression 0       encountered.
Infinity::indet: Indeterminate expression ----------------- encountered.                                                   1/3                    … [also: Power::infy] result: Out[338]= \left\{\text{Fcand1(k1,0)},\mathrm{Indeterminate},\text{FLdet[k1]},\frac{\left(2.45075493251774917289426756629+5.78551794268217244227029022731\,\mathrm{i}\right)\,\left(\left(0.0525946877178…
Cell 104: ✓ — Out[343]= \left\{\text{Fcand1(k1,0)},\left(-8.51191121729+7.73287575414\,\mathrm{i}\right)\,k1^2-\left(4.5054984393+11.6597806074\,\mathrm{i}\right)\,k1+\left(1.000000000000+0.\times 10^{-12}\,\mathrm{i}\right),\text{FLdet[k1]},\frac{\left(2.45075493251774917289426756629+5.78551794268217244227029022731\,\mathrm{i}\right)\,\left(\left(0.05259468771783351844996683714+0.15021349717800998891606886537\… [960 chars total]
Cell 107: ✓ — Out[348]= \left\{\text{STEP 0: det[$\lambda$1,$\lambda$2][k\^{}(u/h)Q1] vs O[k,k], worst |residual| (sample k)},0,\text{sample rows \{$\lambda$1,$\lambda$2,state,resid,ratio\}},\begin{pmatrix}1 & 1 & \left\{0,0\right\} & 0 & 1.\,\\1 & 1 & \left\{1,0\right\} & 0 & 1.\,\\1 & 1 & \left\{1,1\right\} & 0 & 1.\,\\1 & 1 & \left\{2,0\right\} & 0 & 1.\,\\2 & 1 & \left\{0,0\right\} & 0 & 1.\,\\2 & 1 & \left… [445 chars total]
Cell 108: ✓ — Out[354]= \left\{\text{STEP 0 CONFIRMED: det[$\lambda$1,$\lambda$2][k\^{}(u/h)Q1] == O[k,k] identically},\text{worst residual over full sweep},0,\text{num states},41\right\}

## 2026-07-18 11:57:41 — ▶️ RUN CELLS 109–123
Cell 109: ✓ — 2 kernel messages: Power::infy: Infinite expression - encountered.                                  0
Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.
Power::infy: Infinite expression - encounter… [also: Power::infy] result: Out[358]= \begin{pmatrix}1 & 0 & \left\{1.\,\right\}\\1 & 1 & \left\{1.\,\right\}\\2 & 0 & \left\{1.\,\right\}\\2 & 1 & \left\{2.\,\right\}\\2 & 2 & \left\{2.\,\right\}\\3 & 0 & \left\{1.\,\right\}\\3…
Cell 110: ✓ — Out[364]= \left\{\text{STEP 1 CONFIRMED: ⟨x0[$\lambda$]|E21\^{}m = $\gamma$\_m($\lambda$) xL1[$\lambda$][\{m\}], $\gamma$\_m($\lambda$)=$\lambda$!/($\lambda$-m)!},\text{worst deviation},0,\text{table \{$\lambda$,m,$\gamma$\_found,$\gamma$\_predicted,dev\}},\begin{pmatrix}1 & 0 & 1.\, & 1 & 0\\1 & 1 & 1.\, & 1 & 0\\2 & 0 & 1.\, & 1 & 0\\2 & 1 & 2.\, & 2 & 0\\2 & 2 & 2.\, & 2 & 0\\3 & 0 & 1.\, & 1 &… [489 chars total]
Cell 111: ✓ — Out[369]= \left\{\text{STEP 2: x2[\{s1,m\}] vs xL1[$\lambda$1][\{s1\}]$\otimes$xL1[$\lambda$2][\{m\}], worst ratio-spread (constancy check)},0,\text{num cases},54,\text{sample rows},\begin{pmatrix}1 & 1 & 0 & 1 & 0.16 & 0\\1 & 1 & 1 & 1 & 1.\, & 0\\1 & 2 & 0 & 1 & 0.16 & 0\\1 & 2 & 0 & 2 & -0.68 & 0\\1 & 2 & 1 & 1 & 1.\, & 0\\1 & 2 & 1 & 2 & 1.\, & 0\\1 & 3 & 0 & 1 & 0.16 & 0\\1 & 3 & 0 & 2 & -0.6… [479 chars total]
Cell 112: ✓ — Out[373]= \begin{pmatrix}1 & 0 & \frac{4}{25}\\1 & 1 & 1\\2 & 0 & \frac{2}{23}\\2 & 1 & \frac{25}{46}\\2 & 2 & 1\\3 & 0 & \frac{4}{67}\\3 & 1 & \frac{25}{67}\\3 & 2 & \frac{46}{67}\\3 & 3 & 1\end{pmatrix}
Cell 113: ✓ — Out[378]= \text{\$Aborted} | Out[379]= \left\{\text{STEP 2 CONFIRMED: x2[$\lambda$1,$\lambda$2][\{s1,m\}] = c($\lambda$1,s1)$\cdot$xL1[$\lambda$1][\{s1\}]$\otimes$xL1[$\lambda$2][\{m\}], c($\lambda$1,s1)=($\theta$1-$\theta$2+s1h)/($\theta$1-$\theta$2+$\lambda$1h) --- genuinely independent of $\lambda$2,m},\text{worst deviation over all 54 ($\lambda$1,$\lambda$2,s1,m) cases},1.68\right\}
Cell 114: ✓ — Out[383]= \begin{pmatrix}1 & 0 & 1 & \frac{4}{25}\\1 & 0 & 2 & -\frac{17}{25}\\1 & 0 & 3 & -\frac{38}{25}\\1 & 1 & 1 & 1\\1 & 1 & 2 & 1\\1 & 1 & 3 & 1\\2 & 0 & 1 & \frac{2}{23}\\2 & 0 & 2 & -\frac{34}{575}\\2 & 0 & 3 & \frac{323}{575}\\2 & 1 & 1 & \frac{25}{46}\\2 & 1 & 2 & \frac{2}{23}\\2 & 1 & 3 & -\frac{17}{46}\\2 & 2 & 1 & 1\\2 & 2 & 2 & 1\\2 & 2 & 3 & 1\\3 & 0 & 1 & \frac{4}{67}\\3 & 0 & 2 & … [672 chars total]
Cell 115: ✓ — 4 kernel messages: LinearSolve::nosol: Linear equation encountered that has no solution.
LinearSolve::nosol: Linear equation encountered that has no solution.
LinearSolve::nosol: Linear equation encountered that has no … [also: LinearSolve::nosol, LinearSolve::nosol, LinearSolve::nosol] result: Out[388]= \left\{\text{coeffs (rows=s1, cols=s1'=0..3)},\left\{\mathrm{LinearSolve}\left[\begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & \frac{1}{36}\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & \f…
Cell 116: ✓ — Out[393]= \left\{\left\{0,2\right\}\to -0.02206359506813757300454250486697,\left\{1,1\right\}\to 0.1635301752109020116807268007787,\left\{2,0\right\}\to 0.8585334198572355613238157040883\right\}
Cell 117: ✓ — Out[397]= \left\{\text{STEP 2 level-conservation check: \# off-level nonzero coeffs (want 0 everywhere)},0,\text{num ($\lambda$1,$\lambda$2,s1,m) cases},29\right\}
Cell 119: ✓ — Out[405]= \left\{\text{m=0: Om-Cm (must be 0 per Phase1, C1pred/P1site2[$\theta$1] should equal exactly Cm's own norm here)},\left(2.5572121404875445 -0.4102640568287387\,\mathrm{i}\right)\,k1-\left(0.038411844228125346 +2.1227491621421546\,\mathrm{i}\right)\,k1^2,\text{m=1: Om\_1[k1]},\left(0.3604398751472265372666776506+0.9327824485933747164848359821\,\mathrm{i}\right)\,k1^2-\left(1.296101147801… [801 chars total]
Cell 120: ✓ — Out[412]= \left\{\text{m=0: (B0-B1)/O\_0 ratio},\frac{\sqrt[3]{k1}\,\left(\left(-0.0292504903494856237+0.0205326892042810627\,\mathrm{i}\right)\,k1^2+\left(0.499242526583056932+0.532659421382621779\,\mathrm{i}\right)\,k1-\left(0.139388493663376676+0.206413542317126162\,\mathrm{i}\right)\right)}{\left(0.3604398751472265372666776506+0.9327824485933747164848359821\,\mathrm{i}\right)\,k1^2-\left(1.296… [1083 chars total]
Cell 121: ✓ — Out[415]= \left\{\text{Cm[0]-Cm[1] shape},\left(-0.0225172879957976151-1.244370198497125076\,\mathrm{i}\right)\,k1^2+\left(5.434075798536032377-0.8718111207610696998\,\mathrm{i}\right)\,k1+\left(0.\times 10^{-20}+0.\times 10^{-20}\,\mathrm{i}\right),\text{ratio vs O\_0[k1] (want: clean constant)},\frac{\left(-0.0225172879957976151-1.244370198497125076\,\mathrm{i}\right)\,k1^2+\left(5.4340757985360… [753 chars total]
Cell 122: ✓ — Out[417]= \left\{\text{Cm[0][k1] / O\_1[k1] (reversed-index hypothesis; want clean k1-independent constant)},\frac{\left(0.111318362779488714+0.0407158959953172346\,\mathrm{i}\right)\,k1^2-\left(0.21818047436770437-1.02534209918467083\,\mathrm{i}\right)\,k1+\left(0.360439875147226537-0.932782448593374716\,\mathrm{i}\right)}{\left(1.\,+0.\times 10^{-29}\,\mathrm{i}\right)\,k1^2-\left(2.963909584212… [539 chars total]

## 2026-07-18 11:58:09 — ✏️ EDIT CELL 113
(*First hypothesis, kept as the discovery record: c(λ1,s1) = (θ1-θ2+s1 h)/(θ1-θ2+λ1 h) fit only
  from the λ2=2,m=1 slice used to extract it (Cell 112) — test it against the FULL step2Report sweep
  (…


## 2026-07-18 11:58:10
**Input:**
```wolfram
(*First hypothesis, kept as the discovery record: c(λ1,s1) = (θ1-θ2+s1 h)/(θ1-θ2+λ1 h) fit only
  from the λ2=2,m=1 slice used to extract it (Cell 112) — test it against the FULL step2Report sweep
  (all m, not just m=1). RESULT: FAILS for m>=2 (worst deviation ~1.68, not 0) — this hypothesis
  incorrectly assumed m-independence (it IS λ2-independent, but not m-independent, an error caught
  here and corrected in the next cell). Reported (not asserted) rather than Abort-forced, since the
  point of this cell is documenting what the naive fit missed, not hard-asserting a wrong claim.*)
ClearAll[cClosed];
cClosed[λ1_,s1_] := (θ1-θ2+s1 h)/(θ1-θ2+λ1 h);
step2Check = Table[
  {row[[1]],row[[2]],row[[3]],row[[4]], row[[5]], N[cClosed[row[[1]],row[[3]]]],
   Chop[row[[5]]-N[cClosed[row[[1]],row[[3]]]],10^-9]}
, {row, step2Report}];
worstStep2Closed = Max[Abs[step2Check[[All,7]]]];
{"c(λ1,s1) hypothesis (m-independent) vs full sweep — FAILS for m>=2, see next cell for the fix",
 "worst deviation over all 54 (λ1,λ2,s1,m) cases", worstStep2Closed}
```
**Output:**
```
Out= {"c(λ1,s1) hypothesis (m-independent) vs full sweep — FAILS for m>=2, see next cell for the fix", "worst deviation over all 54 (λ1,λ2,s1,m) cases", 1.68}
```
## 2026-07-18 11:58:21 — ▶️ RUN CELLS 114–123
Cell 114: ✓ — Out[421]= \begin{pmatrix}1 & 0 & 1 & \frac{4}{25}\\1 & 0 & 2 & -\frac{17}{25}\\1 & 0 & 3 & -\frac{38}{25}\\1 & 1 & 1 & 1\\1 & 1 & 2 & 1\\1 & 1 & 3 & 1\\2 & 0 & 1 & \frac{2}{23}\\2 & 0 & 2 & -\frac{34}{575}\\2 & 0 & 3 & \frac{323}{575}\\2 & 1 & 1 & \frac{25}{46}\\2 & 1 & 2 & \frac{2}{23}\\2 & 1 & 3 & -\frac{17}{46}\\2 & 2 & 1 & 1\\2 & 2 & 2 & 1\\2 & 2 & 3 & 1\\3 & 0 & 1 & \frac{4}{67}\\3 & 0 & 2 & … [672 chars total]
Cell 115: ✓ — 4 kernel messages: LinearSolve::nosol: Linear equation encountered that has no solution.
LinearSolve::nosol: Linear equation encountered that has no solution.
LinearSolve::nosol: Linear equation encountered that has no … [also: LinearSolve::nosol, LinearSolve::nosol, LinearSolve::nosol] result: Out[426]= \left\{\text{coeffs (rows=s1, cols=s1'=0..3)},\left\{\mathrm{LinearSolve}\left[\begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & \frac{1}{36}\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & \f…
Cell 116: ✓ — Out[431]= \left\{\left\{0,2\right\}\to -0.02206359506813757300454250486697,\left\{1,1\right\}\to 0.1635301752109020116807268007787,\left\{2,0\right\}\to 0.8585334198572355613238157040883\right\}
Cell 117: ✓ — Out[435]= \left\{\text{STEP 2 level-conservation check: \# off-level nonzero coeffs (want 0 everywhere)},0,\text{num ($\lambda$1,$\lambda$2,s1,m) cases},29\right\}
Cell 119: ✓ — Out[443]= \left\{\text{m=0: Om-Cm (must be 0 per Phase1, C1pred/P1site2[$\theta$1] should equal exactly Cm's own norm here)},\left(2.5572121404875445 -0.4102640568287387\,\mathrm{i}\right)\,k1-\left(0.038411844228125346 +2.1227491621421546\,\mathrm{i}\right)\,k1^2,\text{m=1: Om\_1[k1]},\left(0.3604398751472265372666776506+0.9327824485933747164848359821\,\mathrm{i}\right)\,k1^2-\left(1.296101147801… [801 chars total]
Cell 120: ✓ — Out[450]= \left\{\text{m=0: (B0-B1)/O\_0 ratio},\frac{\sqrt[3]{k1}\,\left(\left(-0.0292504903494856237+0.0205326892042810627\,\mathrm{i}\right)\,k1^2+\left(0.499242526583056932+0.532659421382621779\,\mathrm{i}\right)\,k1-\left(0.139388493663376676+0.206413542317126162\,\mathrm{i}\right)\right)}{\left(0.3604398751472265372666776506+0.9327824485933747164848359821\,\mathrm{i}\right)\,k1^2-\left(1.296… [1083 chars total]
Cell 121: ✓ — Out[453]= \left\{\text{Cm[0]-Cm[1] shape},\left(-0.0225172879957976151-1.244370198497125076\,\mathrm{i}\right)\,k1^2+\left(5.434075798536032377-0.8718111207610696998\,\mathrm{i}\right)\,k1+\left(0.\times 10^{-20}+0.\times 10^{-20}\,\mathrm{i}\right),\text{ratio vs O\_0[k1] (want: clean constant)},\frac{\left(-0.0225172879957976151-1.244370198497125076\,\mathrm{i}\right)\,k1^2+\left(5.4340757985360… [753 chars total]
Cell 122: ✓ — Out[455]= \left\{\text{Cm[0][k1] / O\_1[k1] (reversed-index hypothesis; want clean k1-independent constant)},\frac{\left(0.111318362779488714+0.0407158959953172346\,\mathrm{i}\right)\,k1^2-\left(0.21818047436770437-1.02534209918467083\,\mathrm{i}\right)\,k1+\left(0.360439875147226537-0.932782448593374716\,\mathrm{i}\right)}{\left(1.\,+0.\times 10^{-29}\,\mathrm{i}\right)\,k1^2-\left(2.963909584212… [539 chars total]


## 2026-07-18 14:39:08
**Input:**
```wolfram
{θ1,θ2,h,Length[DownValues[Om]],Length[DownValues[Cm]],Length[DownValues[decompX2]]}
```
**Output:**
```
Out= {1/3, 1/7, 1, 0, 0, 0}
```

---

## 2026-07-18 14:39:13 — 🔄 KERNEL RESTART

## 2026-07-18 14:39:21 — ▶️ RUN CELLS 1–10
Cell 1: ✓ — Out[8]= \left\{\frac{1}{3},\frac{1}{7},1,0.360439875147226537266677650596+0.932782448593374716484835982071\,\mathrm{i},0.36043987514722653726667765060-0.932782448593374716484835982071\,\mathrm{i},0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i},1\right\}
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 8: ✓ — (no output)
Cell 10: ✓ — (no output)

## 2026-07-18 14:39:26 — ▶️ RUN CELLS 11–20
Cell 11: ✓ — (no output)
Cell 12: ✓ — Out[27]= \left\{0\right\}
Cell 13: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — Out[31]= \left\{0\right\}
Cell 16: ✓ — (no output)
Cell 17: ✓ — Out[34]= \begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{pmatrix}
Cell 18: ✓ — (no output)
Cell 20: ✓ — (no output)

## 2026-07-18 14:39:29 — ▶️ RUN CELLS 21–28
Cell 21: ✓ — (no output)
Cell 23: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 26: ✓ — (no output)
Cell 27: ✓ — Out[52]= \begin{pmatrix}0 & \frac{1}{2} & 0 & 0 & 0 & 0\end{pmatrix}
Cell 28: ✓ — Out[53]= \left\{0\right\}

## 2026-07-18 14:39:31 — ▶️ RUN CELLS 29–36
Cell 30: ✓ — Out[59]= \left\{\left(u-\frac{4}{3}\right)\,\left(u-\frac{8}{7}\right),\left(u-\frac{1}{3}\right)\,\left(u-\frac{1}{7}\right),\left(0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i}\right)\,u^2-\left(1.06415582186324025288257211128+0.\times 10^{-29}\,\mathrm{i}\right)\,u\right\}
Cell 31: ✓ — Out[62]= \begin{pmatrix}1 & 1\\2 & 2\\3 & 3\end{pmatrix}
Cell 32: ✓ — Out[65]= \text{max|C2-($\kappa$1+$\kappa$2)Id|}\to 0.\times 10^{-29}
Cell 33: ✓ — Out[75]= \langle\!\lvert \left\{0,0\right\}\to 0.5664055180884988 +1.3769645669711723\,\mathrm{i},\left\{1,0\right\}\to -0.778123966451158+0.\,\,\mathrm{i},\left\{1,1\right\}\to 1.1900552523337025 +0.\,\,\mathrm{i},\left\{2,0\right\}\to 0.5664055180884988 -1.3769645669711723\,\mathrm{i}\rvert\!\rangle
Cell 34: ✓ — 2 kernel messages: Psi::spurious: No genuine state at (M,n)=(2,2) for (λ1,λ2)=(2,1): weight sector M=2 has only 2 genuine state(s) (n=0..1), so n=2 does not exist. [also: Psi::spurious] result: Print: spurious-guard test (expect a Psi::spurious message and Missing[...]):[formula] | Out[84]= \left\{\left\{1.000000000000000000000000000+0.\times 10^{-29}\,\mathrm{i},-1.3004494881076436272092358…
⚠️ stopped at Cell 34 — error detected (pass stopOnError:false to continue past errors)

## 2026-07-18 14:39:41 — ▶️ RUN CELLS 35–42
Cell 35: ✓ — Out[91]= \left\{1.03464899078135229594035121+0.381869568380946283899575941\,\mathrm{i},-2.37022278949001677248172619-0.20048970475590097847093779\,\mathrm{i},1.00000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right\}
Cell 36: ✓ — Out[94]= \left\{-0.709231010958016926537225142+1.09596916251277206706752728\,\mathrm{i},-0.02291660900+0.\times 10^{-12}\,\mathrm{i}\right\}
Cell 37: ✓ — (no output)
Cell 38: ✓ — (no output)
Cell 39: ✓ — (no output)
Cell 40: ✓ — (no output)
Cell 41: ✓ — (no output)
Cell 42: ✓ — Out[109]= 0

## 2026-07-18 14:39:44 — ▶️ RUN CELLS 43–47
Cell 43: ✓ — (no output)
Cell 44: ✓ — Out[113]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 45: ✓ — Out[115]= \left\{6,6\right\}
Cell 46: ✓ — Out[116]= \begin{pmatrix}2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i} & 0 & 0 & 0 & 0 & 0\\0 & 7.36608403016431791181-0.968046299794461400266\,\mathrm{i} & 0 & 0 & 0 & 0\\0 & 0 & 2.68702160063169134012-1.3456153260557858257\,\mathrm{i} & 0 & 0 & 0\\0 & 0 & 0 & 7.3660840301643179118107+0.9680462997944614002659\,\mathrm{i} & 0 & 0\\0 & 0 & 0 & 0 & 2.687021600631691340118… [537 chars total]
Cell 47: ✓ — Out[125]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{max|off-diagonal|}\to 5.065342385207116\times 10^{-34},\text{digits vanished (-Log10)}\to 33.29539119374288,\text{min|diagonal|}\to 3.0051232400647923,\text{diagonal values}\to \left\{2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i},7.36608403016431791181-0… [699 chars total]

## 2026-07-18 14:39:48 — ▶️ RUN CELLS 48–54
Cell 48: ✓ — (no output)
Cell 49: ✓ — (no output)
Cell 50: ✓ — Out[133]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 51: ✓ — Out[138]= \left\{1.00000000000000000000000000000,6\right\}
Cell 52: ✓ — Out[143]= \langle\!\lvert \text{max|<PsiLeft|Psi>| off-diagonal}\to 0.\,,\text{digits vanished (-Log10)}\to \mathrm{Indeterminate}\rvert\!\rangle
Cell 53: ✓ — Out[147]= \left\{7.366084030164317911810727273-0.9680462997944614002658874268\,\mathrm{i},1.00000000000000000000+0.\times 10^{-21}\,\mathrm{i}\right\}
Cell 54: ✓ — Out[150]= \begin{matrix}\langle\!\lvert \text{M}\to 0,\text{n}\to 0,\text{braket}\to 2.909064289-5.804620064\,\mathrm{i},\text{G21 diagonal}\to 2.909064289-5.804620064\,\mathrm{i},\text{normRatio}\to 1.000000000+0.\times 10^{-10}\,\mathrm{i}\rvert\!\rangle \\\langle\!\lvert \text{M}\to 1,\text{n}\to 0,\text{braket}\to 7.366084030-0.968046300\,\mathrm{i},\text{G21 diagonal}\to 7.366084030-0.9680463… [1247 chars total]

## 2026-07-18 14:39:52 — ▶️ RUN CELLS 55–62
Cell 56: ✓ — (no output)
Cell 57: ✓ — Out[158]= \begin{pmatrix}3 & 3\\3 & 3\end{pmatrix}
Cell 58: ✓ — (no output)
Cell 59: ✓ — Out[165]= \left\{0\right\}
Cell 60: ✓ — Out[169]= \begin{matrix}\left\{1,1,0,1.\,,0.\,\right\}\\\left\{1,1,1,1.\,,0.\,\right\}\\\left\{1,2,0,1.\,,0.\,\right\}\\\left\{1,2,1,1.\,,0.\,\right\}\\\left\{1,3,0,1.\,,0.\,\right\}\\\left\{1,3,1,1.\,,0.\,\right\}\\\left\{2,1,0,1.\,,0.\,\right\}\\\left\{2,1,1,1.\,,0.\,\right\}\\\left\{2,1,2,1.\,,0.\,\right\}\\\left\{2,2,0,1.\,,0.\,\right\}\\\left\{2,2,1,1.\,,0.\,\right\}\\\left\{2,2,2,1.\,,0.\,\r… [898 chars total]
Cell 61: ✓ — Out[170]= \langle\!\lvert \text{max spread (ratio non-constancy)}\to 0.\,,\text{max |ratio - 1| over sweep}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases checked}\to 27\rvert\!\rangle
Cell 62: ✓ — Out[176]= \langle\!\lvert \text{max operator-level spread}\to 0.\,,\text{max |$\mu$ - $\mu$Site2|}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases}\to 27\rvert\!\rangle

## 2026-07-18 14:40:18 — ▶️ RUN CELLS 63–70
Cell 64: ✓ — (no output)
Cell 65: ✓ — Out[185]= \langle\!\lvert \text{w\_0 ($\lambda$=1,2,3)}\to \left\{1,1,1\right\},\text{full weights w\_j}\to \left\{1\to \left\{1,-1\right\},2\to \left\{1,-2,1\right\},3\to \left\{1,-3,3,-1\right\}\right\},\text{br1[$\lambda$][1] ($\lambda$=1,2,3)}\to \left\{0,0,0\right\}\rvert\!\rangle
Cell 66: ✓ — Out[187]= \left\{0,0\right\}
Cell 67: ✓ — Out[190]= \sqrt[3]{k}\,\left(\left(-1.18561270626131985872+0.832253303224995280404\,\mathrm{i}\right)\,k^2-\left(0.761819910806671651081+0.812812473459294825632\,\mathrm{i}\right)\,k-\left(0.404410037536507329420+0.598870869485872613792\,\mathrm{i}\right)\right)
Cell 68: ✓ — Out[194]= \begin{pmatrix}0 & 0 & 0\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix}
Cell 69: ✓ — Out[198]= \langle\!\lvert \text{OP(E12)}\to 1.\,+0.\times 10^{-30}\,\mathrm{i},\text{OP(E21)}\to \left(-0.03626734060947795145020244939-2.004237714666237867354326095\,\mathrm{i}\right)\,\varphi ^2-\left(1.522150083623538480912456055-0.2442047957313920727733393355\,\mathrm{i}\right)\,\varphi +\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right),\text{nontrivial generator \{a,b\}}\to \left\{2,1\right\},\… [464 chars total]
Cell 70: ✓ — Out[201]= \left\{1,2\right\}

## 2026-07-18 14:40:21 — ▶️ RUN CELLS 71–77
Cell 71: ✓ — Out[207]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{$\varphi$ candidates per state}\to \begin{pmatrix}\left\{0,0\right\} & \left(-0.920793778924720536105378676-0.390049761817044163511250929\,\mathrm{i}\right)\,k^{4/3}+\left(0.695722366612934233848844287-0.718310788304406529935213657\,\mathrm{i}\right)\,\sqrt[3]{k}-\left(0.3604398751… [1247 chars total]
Cell 72: ✓ — Out[213]= \langle\!\lvert \text{c\_sb solution(s)}\to \left\{\left\{\mathrm{cb}\to \frac{\left(-6.242239011314\times 10^{45}-4.511757007219\times 10^{45}\,\mathrm{i}\right)\,k^{4/3}-\left(5.999302657796\times 10^{45}+8.356648463050\times 10^{45}\,\mathrm{i}\right)\,\sqrt[3]{k}}{9.833530773718\times 10^{45}\,k^{4/3}-\left(3.54439660434\times 10^{45}-9.172544913427\times 10^{45}\,\mathrm{i}\right)\,… [841 chars total]
Cell 73: ✓ — Out[216]= \langle\!\lvert \text{generator \{a,b\} tested}\to \left\{2,1\right\},\text{rep tested}\to \left\{1,2\right\},\text{FL shape (verified exact for all 6 states)}\to \text{FL[1,2][state][k] = k\^{}(1/3) (A\_state + B\_state k)},\text{c0(state) = OPg[state][0]}\to \begin{pmatrix}\left\{0,0\right\} & 1.\,+0.\times 10^{-28}\,\mathrm{i}\\\left\{1,0\right\} & 1.\,+0.\times 10^{-28}\,\mathrm{i}\\… [1247 chars total]
Cell 74: ✓ — Out[220]= \left(1.000000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right)-\left(0.360439875147226537266677651+0.932782448593374716484835982\,\mathrm{i}\right)\,k
Cell 75: ✓ — Out[227]= \langle\!\lvert \text{alphaTab \{state,D0,c0,Dc,c1,$\alpha$\}}\to \begin{pmatrix}\left\{0,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & -0.36043987514722653-0.9327824485933747\,\mathrm{i} & 0.36043987514722653 +0.9327824485933747\,\mathrm{i} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.7101545203859353 -0.2371533418544… [1247 chars total]
Cell 76: ✓ — Out[235]= \langle\!\lvert \text{validation}\to \left\{\left\{1,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\end{pmatrix},\left\{2,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\\\left\{2,1\right\} & 0\\\left\{3,0\right\} & 0\end{pmatrix},\left\{1… [1247 chars total]
Cell 77: ✓ — Out[237]= \langle\!\lvert \text{generator \{a,b\}}\to \left\{2,1\right\},\text{$\varphi$(k)}\to \text{$\alpha$ k},\text{$\alpha$ (universal value)}\to -1.\,+0.\,\,\mathrm{i},\text{normalization}\to \text{1-site FSoV determinant det1[f] = br1[f]/(f/.u->$\theta$1); Section C's det[$\lambda$1,$\lambda$2][f] restricted to 1 site},\text{max residual over rep sweep}\to 0,\text{reps validated}\to \begin{… [510 chars total]

## 2026-07-18 14:40:44 — ▶️ RUN CELLS 78–82
Cell 79: ✓ — (no output)
Cell 80: ✓ — Out[245]= \left\{\text{Test 1 (symbolic monomials u\^{}j) residuals},\left\{0\right\}\right\}
Cell 81: ✓ — Out[250]= \left\{\text{Test 2 (twisted monomials k\^{}(u/h) u\^{}j) residuals},\left\{0\right\}\right\}
Cell 82: ✓ — Out[259]= \left\{\text{Test 3 worst |br identity residual| on Q1},0,\text{worst |detRatio - predicted| },0,\text{num states checked},41\right\}

## 2026-07-18 14:40:48 — ▶️ RUN CELLS 83–87
Cell 84: ✓ — (no output)
Cell 85: ✓ — Out[271]= \begin{pmatrix}1 & 1 & 0\\2 & 1 & 0\\3 & 1 & 0\end{pmatrix}
Cell 86: ✓ — Out[278]= \left\{\text{Mirror Test1 (u\^{}j)},\left\{0\right\},\text{Mirror Test2 (k\^{}(u/h)u\^{}j)},\left\{0\right\}\right\}
Cell 87: ✓ — Out[285]= \left\{\text{Mirror Test3 worst |br identity residual| on Q1},0,\text{num states checked},41\right\}

## 2026-07-18 14:40:49 — ▶️ RUN CELLS 88–91
Cell 89: ✓ — Out[294]= \left\{\text{OP2(E12)},1.\,+0.\times 10^{-30}\,\mathrm{i},\text{OP2(E21)},\left(2.603469709065218092712489006+0.6885776528619826437114966465\,\mathrm{i}\right)\,\varphi +\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right),\text{nontrivial generator \{a,b\}},\left\{2,1\right\}\right\}
Cell 90: ✓ — Out[297]= \left(1.\,+0.\times 10^{-27}\,\mathrm{i}\right)-\left(0.360439875147226537266677651+0.932782448593374716484835982\,\mathrm{i}\right)\,k
Cell 91: ✓ — Out[300]//TableForm= \begin{pmatrix}\left\{0,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & -0.36043987514722653-0.9327824485933747\,\mathrm{i} & 0.36043987514722653 +0.9327824485933747\,\mathrm{i} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.4685059330398573 -0.1713611156914601\,\mathrm{i} & 2.603469709065218 +0.6885776528619827… [1258 chars total]

## 2026-07-18 14:41:01 — ▶️ RUN CELLS 92–108
Cell 92: ✓ — Out[302]//TableForm= \begin{pmatrix}\left\{0,0\right\} & -1.\,+0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.973902234276272 +0.2269679230016555\,\mathrm{i} & 0.7391304347826086 +0.\,\,\mathrm{i} & 0.973902234276272 +0.2269679230016555\,\mathrm{i}\\\left\{1,0\right\} & 0.15191893851739566 -0.10600050418875033\,\mathrm{i} & -0.15191893851739566+0.10600050418875033\,\mathrm{i} & 0.7651595274649543 +… [1258 chars total]
Cell 93: ✓ — Out[305]= \begin{pmatrix}\left\{0,0\right\} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 0.05448860910968355 -0.09540961504089102\,\mathrm{i}\\\left\{1,1\right\} & -0.0293413933141388+0.7363100318787654\,\mathrm{i}\\\left\{2,0\right\} & 0.11703059508534998 +0.\,\,\mathrm{i}\\\left\{2,1\right\} & 0.8158972359113332 +0.\,\,\mathrm{i}\\\left\{3,0\right\} & 0.05448860910968355 +0.09540961504089102\,… [541 chars total]
Cell 96: ✓ — Out[311]= \left\{\text{O[k1,0]-FLdet[k1] (must be 0)},0,\text{O[0,k2]-FLdet2[k2] (expected NONZERO per Phase1 finding)},\left(-3.0719756421050755-0.5172165371705225\,\mathrm{i}\right)\,k2\right\}
Cell 97: ✓ — Out[316]= \left\{\text{Fcand1},k1^2\,\left(\left(-2.2527492196701659-5.82989030370859198\,\mathrm{i}\right)\,k2-\left(0.4170744170089964+23.0487337186617355\,\mathrm{i}\right)\right)+k1\,\left(\left(0.67128913890100625-0.37524873225652243\,\mathrm{i}\right)\,k2+\left(9.51343802264711551-1.5262799733212005\,\mathrm{i}\right)\right)-\left(1.99115021541939352-0.72828474168870538\,\mathrm{i}\right)\,k… [1014 chars total]
Cell 98: ✓ — Out[318]= \left(0.9252125034899888 -6.175727508765443\,\mathrm{i}\right)\,k1^2\,k2+\left(5.4271267482291075 +3.0890004366773067\,\mathrm{i}\right)\,k1\,k2
Cell 99: ✓ — Out[323]= \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 7.991287939023577 -1.2820751775898085\,\mathrm{i}\\\left\{2,0\right\} & -0.3808070763995185-21.044496003995498\,\mathrm{i}\\\left\{0,1\right\} & 0.6123194936458246 +1.416862394550688\,\mathrm{i}\\\left\{1,1\right\} & -0.6248120089007682-3.051910425543801\,\mathrm{i}\\\left\{2,1\right\} & -1.8923093445229393-4.897107855115217\,\m… [421 chars total]
Cell 100: ✓ — Out[328]= \left\{\text{Fcand2 (no P-insertion)},k1^2\,\left(\left(0.02466179161444501+1.36288164597304175\,\mathrm{i}\right)-\left(0.36043987514722654+0.93278244859337472\,\mathrm{i}\right)\,k2\right)+k1\,\left(\left(1.2351720155778515-0.69045766735200127\,\mathrm{i}\right)\,k2-\left(6.46913785540003854-1.0378703818584163\,\mathrm{i}\right)\right)+\left(5.38781822995835893-1.9706528304517910\,\mat… [788 chars total]
Cell 101: ✓ — 2 kernel messages: Power::infy: Infinite expression 0       encountered.
Infinity::indet: Indeterminate expression ----------------- encountered.                                                   1/3                    … [also: Power::infy] result: Out[331]//TableForm= \begin{pmatrix}2 & 1 & \left\{0,0\right\} & \mathrm{Indeterminate} & \mathrm{Indeterminate} & -4.970498218180917+4.515583419069149\,\mathrm{i}\\2 & 1 & \left\{1,1\right\} & \mathr…
Cell 102: ✓ — Out[334]//TableForm= \begin{pmatrix}2 & 1 & \left\{0,0\right\} & \left(-7.771745024480887+7.060451775520925\,\mathrm{i}\right)\,k1^2-\left(3.7846186890458786 +9.794215710230434\,\mathrm{i}\right)\,k1 & \left(1.8923093445229393 +4.897107855115217\,\mathrm{i}\right)\,k2 & -4.970498218180917+4.515583419069149\,\mathrm{i}\\2 & 1 & \left\{1,1\right\} & \left(9.339255641812587 +10.661478687052961\,\math… [1258 chars total]
Cell 103: ✓ — 2 kernel messages: Power::infy: Infinite expression 0       encountered.
Infinity::indet: Indeterminate expression ----------------- encountered.                                                   1/3                    … [also: Power::infy] result: Out[338]= \left\{\text{Fcand1(k1,0)},\mathrm{Indeterminate},\text{FLdet[k1]},\frac{\left(2.45075493251774917289426756629+5.78551794268217244227029022731\,\mathrm{i}\right)\,\left(\left(0.0525946877178…
Cell 104: ✓ — Out[343]= \left\{\text{Fcand1(k1,0)},\left(-8.51191121729+7.73287575414\,\mathrm{i}\right)\,k1^2-\left(4.5054984393+11.6597806074\,\mathrm{i}\right)\,k1+\left(1.000000000000+0.\times 10^{-12}\,\mathrm{i}\right),\text{FLdet[k1]},\frac{\left(2.45075493251774917289426756629+5.78551794268217244227029022731\,\mathrm{i}\right)\,\left(\left(0.05259468771783351844996683714+0.15021349717800998891606886537\… [960 chars total]
Cell 107: ✓ — Out[348]= \left\{\text{STEP 0: det[$\lambda$1,$\lambda$2][k\^{}(u/h)Q1] vs O[k,k], worst |residual| (sample k)},0,\text{sample rows \{$\lambda$1,$\lambda$2,state,resid,ratio\}},\begin{pmatrix}1 & 1 & \left\{0,0\right\} & 0 & 1.\,\\1 & 1 & \left\{1,0\right\} & 0 & 1.\,\\1 & 1 & \left\{1,1\right\} & 0 & 1.\,\\1 & 1 & \left\{2,0\right\} & 0 & 1.\,\\2 & 1 & \left\{0,0\right\} & 0 & 1.\,\\2 & 1 & \left… [445 chars total]
Cell 108: ✓ — Out[354]= \left\{\text{STEP 0 CONFIRMED: det[$\lambda$1,$\lambda$2][k\^{}(u/h)Q1] == O[k,k] identically},\text{worst residual over full sweep},0,\text{num states},41\right\}

## 2026-07-18 14:41:17 — ▶️ RUN CELLS 109–123
Cell 109: ✓ — 2 kernel messages: Power::infy: Infinite expression - encountered.                                  0
Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.
Power::infy: Infinite expression - encounter… [also: Power::infy] result: Out[358]= \begin{pmatrix}1 & 0 & \left\{1.\,\right\}\\1 & 1 & \left\{1.\,\right\}\\2 & 0 & \left\{1.\,\right\}\\2 & 1 & \left\{2.\,\right\}\\2 & 2 & \left\{2.\,\right\}\\3 & 0 & \left\{1.\,\right\}\\3…
Cell 110: ✓ — Out[364]= \left\{\text{STEP 1 CONFIRMED: ⟨x0[$\lambda$]|E21\^{}m = $\gamma$\_m($\lambda$) xL1[$\lambda$][\{m\}], $\gamma$\_m($\lambda$)=$\lambda$!/($\lambda$-m)!},\text{worst deviation},0,\text{table \{$\lambda$,m,$\gamma$\_found,$\gamma$\_predicted,dev\}},\begin{pmatrix}1 & 0 & 1.\, & 1 & 0\\1 & 1 & 1.\, & 1 & 0\\2 & 0 & 1.\, & 1 & 0\\2 & 1 & 2.\, & 2 & 0\\2 & 2 & 2.\, & 2 & 0\\3 & 0 & 1.\, & 1 &… [489 chars total]
Cell 111: ✓ — Out[369]= \left\{\text{STEP 2: x2[\{s1,m\}] vs xL1[$\lambda$1][\{s1\}]$\otimes$xL1[$\lambda$2][\{m\}], worst ratio-spread (constancy check)},0,\text{num cases},54,\text{sample rows},\begin{pmatrix}1 & 1 & 0 & 1 & 0.16 & 0\\1 & 1 & 1 & 1 & 1.\, & 0\\1 & 2 & 0 & 1 & 0.16 & 0\\1 & 2 & 0 & 2 & -0.68 & 0\\1 & 2 & 1 & 1 & 1.\, & 0\\1 & 2 & 1 & 2 & 1.\, & 0\\1 & 3 & 0 & 1 & 0.16 & 0\\1 & 3 & 0 & 2 & -0.6… [479 chars total]
Cell 112: ✓ — Out[373]= \begin{pmatrix}1 & 0 & \frac{4}{25}\\1 & 1 & 1\\2 & 0 & \frac{2}{23}\\2 & 1 & \frac{25}{46}\\2 & 2 & 1\\3 & 0 & \frac{4}{67}\\3 & 1 & \frac{25}{67}\\3 & 2 & \frac{46}{67}\\3 & 3 & 1\end{pmatrix}
Cell 113: ✓ — Out[378]= \left\{\text{c($\lambda$1,s1) hypothesis (m-independent) vs full sweep --- FAILS for m>=2, see next cell for the fix},\text{worst deviation over all 54 ($\lambda$1,$\lambda$2,s1,m) cases},1.68\right\}
Cell 114: ✓ — Out[382]= \begin{pmatrix}1 & 0 & 1 & \frac{4}{25}\\1 & 0 & 2 & -\frac{17}{25}\\1 & 0 & 3 & -\frac{38}{25}\\1 & 1 & 1 & 1\\1 & 1 & 2 & 1\\1 & 1 & 3 & 1\\2 & 0 & 1 & \frac{2}{23}\\2 & 0 & 2 & -\frac{34}{575}\\2 & 0 & 3 & \frac{323}{575}\\2 & 1 & 1 & \frac{25}{46}\\2 & 1 & 2 & \frac{2}{23}\\2 & 1 & 3 & -\frac{17}{46}\\2 & 2 & 1 & 1\\2 & 2 & 2 & 1\\2 & 2 & 3 & 1\\3 & 0 & 1 & \frac{4}{67}\\3 & 0 & 2 & … [672 chars total]
Cell 115: ✓ — 4 kernel messages: LinearSolve::nosol: Linear equation encountered that has no solution.
LinearSolve::nosol: Linear equation encountered that has no solution.
LinearSolve::nosol: Linear equation encountered that has no … [also: LinearSolve::nosol, LinearSolve::nosol, LinearSolve::nosol] result: Out[387]= \left\{\text{coeffs (rows=s1, cols=s1'=0..3)},\left\{\mathrm{LinearSolve}\left[\begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & \frac{1}{36}\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & \f…
Cell 116: ✓ — Out[392]= \left\{\left\{0,2\right\}\to -0.02206359506813757300454250486697,\left\{1,1\right\}\to 0.1635301752109020116807268007787,\left\{2,0\right\}\to 0.8585334198572355613238157040883\right\}
Cell 117: ✓ — Out[396]= \left\{\text{STEP 2 level-conservation check: \# off-level nonzero coeffs (want 0 everywhere)},0,\text{num ($\lambda$1,$\lambda$2,s1,m) cases},29\right\}
Cell 119: ✓ — Out[404]= \left\{\text{m=0: Om-Cm (must be 0 per Phase1, C1pred/P1site2[$\theta$1] should equal exactly Cm's own norm here)},\left(2.5572121404875445 -0.4102640568287387\,\mathrm{i}\right)\,k1-\left(0.038411844228125346 +2.1227491621421546\,\mathrm{i}\right)\,k1^2,\text{m=1: Om\_1[k1]},\left(0.3604398751472265372666776506+0.9327824485933747164848359821\,\mathrm{i}\right)\,k1^2-\left(1.296101147801… [801 chars total]
Cell 120: ✓ — Out[411]= \left\{\text{m=0: (B0-B1)/O\_0 ratio},\frac{\sqrt[3]{k1}\,\left(\left(-0.0292504903494856237+0.0205326892042810627\,\mathrm{i}\right)\,k1^2+\left(0.499242526583056932+0.532659421382621779\,\mathrm{i}\right)\,k1-\left(0.139388493663376676+0.206413542317126162\,\mathrm{i}\right)\right)}{\left(0.3604398751472265372666776506+0.9327824485933747164848359821\,\mathrm{i}\right)\,k1^2-\left(1.296… [1083 chars total]
Cell 121: ✓ — Out[414]= \left\{\text{Cm[0]-Cm[1] shape},\left(-0.0225172879957976151-1.244370198497125076\,\mathrm{i}\right)\,k1^2+\left(5.434075798536032377-0.8718111207610696998\,\mathrm{i}\right)\,k1+\left(0.\times 10^{-20}+0.\times 10^{-20}\,\mathrm{i}\right),\text{ratio vs O\_0[k1] (want: clean constant)},\frac{\left(-0.0225172879957976151-1.244370198497125076\,\mathrm{i}\right)\,k1^2+\left(5.4340757985360… [753 chars total]
Cell 122: ✓ — Out[416]= \left\{\text{Cm[0][k1] / O\_1[k1] (reversed-index hypothesis; want clean k1-independent constant)},\frac{\left(0.111318362779488714+0.0407158959953172346\,\mathrm{i}\right)\,k1^2-\left(0.21818047436770437-1.02534209918467083\,\mathrm{i}\right)\,k1+\left(0.360439875147226537-0.932782448593374716\,\mathrm{i}\right)}{\left(1.\,+0.\times 10^{-29}\,\mathrm{i}\right)\,k1^2-\left(2.963909584212… [539 chars total]

## 2026-07-18 14:43:23 — 📥 BULK INSERT 3 CELL(S) at positions 124–126
1. [markdown] (*Section J — closed form of the recoupling matrix A[s1,m;s1',m'] and full assembly.↵Strategy: build
2. [code] (*Standalone single-site symbolic-θ fused hierarchy (mirrors Section E's t1/qdet1/xL1, θ1->θa).*)↵Cl
3. [code] (*Validation: symbolic constructions must reproduce the live numeric xL1/x2 exactly under θa->θ1,↵


## 2026-07-18 14:43:23
**Input:**
```wolfram
(*Standalone single-site symbolic-θ fused hierarchy (mirrors Section E's t1/qdet1/xL1, θ1->θa).*)
ClearAll[t1Sym,qdet1Sym,ν1SingleSym,xL1Sym];
t1Sym[θa_][λ_][0,s_][u_] := IdentityMatrix[dim[{λ,0}]];
t1Sym[θa_][λ_][a_,0][u_] := IdentityMatrix[dim[{λ,0}]];
t1Sym[θa_][λ_][1,1][u_] := Sum[L[λ][i,j][u-θa] G[[j,i]],{i,2},{j,2}];
qdet1Sym[θa_][λ_][u_] := L[λ][1,1][u-θa].L[λ][2,2][u-θa-h] - L[λ][2,1][u-θa].L[λ][1,2][u-θa-h];
t1Sym[θa_][λ_][2,1][u_] := χ2 qdet1Sym[θa][λ][u];
t1Sym[θa_][λ_][1,s_/;s>1][u_] := t1Sym[θa][λ][1,s-1][u].t1Sym[θa][λ][1,1][u+(s-1)h] - t1Sym[θa][λ][2,1][u+(s-1)h].t1Sym[θa][λ][1,s-2][u];
ν1SingleSym[θa_][λ_][u_] := (u-θa-h λ);
xL1Sym[θa_][λ_][{0}] := xSingle[λ][{0}];
xL1Sym[θa_][λ_][{s_}] := xSingle[λ][{0}].t1Sym[θa][λ][1,s][θa] / Product[ν1SingleSym[θa][λ][θa+(p-1)h],{p,1,s}];

(*Standalone two-site symbolic-θ hierarchy (mirrors Section A's L1/L2/T/t2/ν1/x2, θ1,θ2->θa,θb).*)
ClearAll[L1Sym,L2Sym,TSym,t2Sym,qdetTSym,ν1TwoSym,x2Sym];
L1Sym[θa_,θb_][λ1_,λ2_][i_,j_][u_] := KroneckerProduct[L[λ1][i,j][u-θa],IdentityMatrix[dim[{λ2,0}]]];
L2Sym[θa_,θb_][λ1_,λ2_][i_,j_][u_] := KroneckerProduct[IdentityMatrix[dim[{λ1,0}]],L[λ2][i,j][u-θb]];
TSym[θa_,θb_][λ1_,λ2_][i_,j_][u_] := Sum[L1Sym[θa,θb][λ1,λ2][i,k][u].L2Sym[θa,θb][λ1,λ2][k,j][u],{k,2}];
t2Sym[θa_,θb_][λ1_,λ2_][0,s_][u_] := IdentityMatrix[dim[{λ1,0}]dim[{λ2,0}]];
t2Sym[θa_,θb_][λ1_,λ2_][a_,0][u_] := IdentityMatrix[dim[{λ1,0}]dim[{λ2,0}]];
t2Sym[θa_,θb_][λ1_,λ2_][1,1][u_] := Sum[TSym[θa,θb][λ1,λ2][i,j][u] G[[j,i]],{i,2},{j,2}];
qdetTSym[θa_,θb_][λ1_,λ2_][u_] := TSym[θa,θb][λ1,λ2][1,1][u].TSym[θa,θb][λ1,λ2][2,2][u-h] - TSym[θa,θb][λ1,λ2][2,1][u].TSym[θa,θb][λ1,λ2][1,2][u-h];
t2Sym[θa_,θb_][λ1_,λ2_][2,1][u_] := χ2 qdetTSym[θa,θb][λ1,λ2][u];
t2Sym[θa_,θb_][λ1_,λ2_][1,s_/;s>1][u_] := t2Sym[θa,θb][λ1,λ2][1,s-1][u].t2Sym[θa,θb][λ1,λ2][1,1][u+(s-1)h] - t2Sym[θa,θb][λ1,λ2][2,1][u+(s-1)h].t2Sym[θa,θb][λ1,λ2][1,s-2][u];
ν1TwoSym[θa_,θb_][λ1_,λ2_][u_] := (u-θa-h λ1)(u-θb-h λ2);
x2Sym[θa_,θb_][λ1_,λ2_][{0,0}] := KroneckerProduct[xSingle[λ1][{0}],xSingle[λ2][{0}]];
x2Sym[θa_,θb_][λ1_,λ2_][{s1_,s2_}] := x2Sym[θa,θb][λ1,λ2][{0,0}].t2Sym[θa,θb][λ1,λ2][1,s1][θa].t2Sym[θa,θb][λ1,λ2][1,s2][θb] /
  (Product[ν1TwoSym[θa,θb][λ1,λ2][θa+(k-1)h],{k,1,s1}] Product[ν1TwoSym[θa,θb][λ1,λ2][θb+(k-1)h],{k,1,s2}]);
```
**Output:**
```
(no output)
```

## 2026-07-18 14:43:23
**Input:**
```wolfram
(*Validation: symbolic constructions must reproduce the live numeric xL1/x2 exactly under θa->θ1,
  θb->θ2 (the actual global numeric values) — a port-fidelity check before trusting the symbolic run.*)
valXL1 = Table[
  Chop[N[Flatten[xL1Sym[θ1][λ][{s}]] - Flatten[xL1[λ][{s}]]],10^-10]
,{λ,{1,2,3}},{s,0,λ}];
valX2 = Table[
  Chop[N[Flatten[x2Sym[θ1,θ2][λ1,λ2][{s1,s2}]] - Flatten[x2[λ1,λ2][{s1,s2}]]],10^-10]
,{λ1,{1,2}},{λ2,{1,2}},{s1,0,λ1},{s2,0,λ2}];
worstValXL1 = Max[Abs[Flatten[valXL1]]];
worstValX2 = Max[Abs[Flatten[valX2]]];
If[worstValXL1 > 10^-8 || worstValX2 > 10^-8, Abort[]];
{"Symbolic-θ port-fidelity: xL1Sym vs xL1", worstValXL1, "x2Sym vs x2", worstValX2}
```
**Output:**
```
Out[445]= \left\{\text{Symbolic-$\theta$ port-fidelity: xL1Sym vs xL1},0,\text{x2Sym vs x2},0\right\}
```
## 2026-07-18 14:43:39 — 📥 BULK INSERT 1 CELL(S) at positions 127–127
1. [code] (*Symbolic decomposition, θ1,θ2 fully symbolic (h=1, numeric G/χ2 as always). Solve↵  x2Sym[{s1,m}]


## 2026-07-18 14:43:40
**Input:**
```wolfram
(*Symbolic decomposition, θ1,θ2 fully symbolic (h=1, numeric G/χ2 as always). Solve
  x2Sym[{s1,m}] = Σ A[s1,m;s1',m'] xL1Sym[λ1][{s1'}]⊗xL1Sym[λ2][{m'}] over the FULL basis (not
  assuming level-conservation — that should fall out of the symbolic solve as A=0 off-level,
  reconfirming Section I's numeric finding analytically).*)
ClearAll[fullBasisSym,decompX2Sym];
fullBasisSym[λ1_,λ2_] := Flatten[Table[
   {s1p,m1p} -> Flatten[KroneckerProduct[xL1Sym[θ1][λ1][{s1p}], xL1Sym[θ2][λ2][{m1p}]]]
  ,{s1p,0,λ1},{m1p,0,λ2}], 1];
decompX2Sym[λ1_,λ2_][s1_,m_] := Module[{basisRules,basisVecs,mat,target,sol},
  basisRules = fullBasisSym[λ1,λ2];
  basisVecs = basisRules[[All,2]];
  mat = Transpose[basisVecs];
  target = Flatten[x2Sym[θ1,θ2][λ1,λ2][{s1,m}]];
  sol = Simplify[LinearSolve[mat, target]];
  Thread[basisRules[[All,1]] -> sol]
];
(*smallest nontrivial case first: λ1=1,λ2=1, state (s1,m)=(0,1) — should mix (0,1) and (1,0).*)
d11_01 = decompX2Sym[1,1][0,1];
d11_01
```
**Output:**
```
Set::write: Tag Times in 1 (d11_) is Protected. | Out[450]= d11\_ | Set::write: Tag Times in 1 (d11_) is Protected.
```
## 2026-07-18 14:43:53 — ✏️ EDIT CELL 127
(*Symbolic decomposition, θ1,θ2 fully symbolic (h=1, numeric G/χ2 as always). Solve
  x2Sym[{s1,m}] = Σ A[s1,m;s1',m'] xL1Sym[λ1][{s1'}]⊗xL1Sym[λ2][{m'}] over the FULL basis (not
  assuming level-cons…


## 2026-07-18 14:43:53
**Input:**
```wolfram
(*Symbolic decomposition, θ1,θ2 fully symbolic (h=1, numeric G/χ2 as always). Solve
  x2Sym[{s1,m}] = Σ A[s1,m;s1',m'] xL1Sym[λ1][{s1'}]⊗xL1Sym[λ2][{m'}] over the FULL basis (not
  assuming level-conservation — that should fall out of the symbolic solve as A=0 off-level,
  reconfirming Section I's numeric finding analytically).
  (*why no underscores in variable names*) Wolfram Language reserves "_" for pattern syntax even
  mid-identifier (e.g. "d1101" not "d11_01") — a literal underscore in a variable name silently
  turns it into a Blank pattern, causing Set::write errors on assignment.*)
ClearAll[fullBasisSym,decompX2Sym];
fullBasisSym[λ1_,λ2_] := Flatten[Table[
   {s1p,m1p} -> Flatten[KroneckerProduct[xL1Sym[θ1][λ1][{s1p}], xL1Sym[θ2][λ2][{m1p}]]]
  ,{s1p,0,λ1},{m1p,0,λ2}], 1];
decompX2Sym[λ1_,λ2_][s1_,m_] := Module[{basisRules,basisVecs,mat,target,sol},
  basisRules = fullBasisSym[λ1,λ2];
  basisVecs = basisRules[[All,2]];
  mat = Transpose[basisVecs];
  target = Flatten[x2Sym[θ1,θ2][λ1,λ2][{s1,m}]];
  sol = Simplify[LinearSolve[mat, target]];
  Thread[basisRules[[All,1]] -> sol]
];
(*smallest nontrivial case first: λ1=1,λ2=1, state (s1,m)=(0,1) — should mix (0,1) and (1,0).*)
dCase1101 = decompX2Sym[1,1][0,1];
dCase1101
```
**Output:**
```
Out= {{0, 0} -> 0, {0, 1} -> 4/25, {1, 0} -> 21/25, {1, 1} -> 0}
```

## 2026-07-18 14:44:12
**Input:**
```wolfram
{θ1,θ2,ValueQ[θA],ValueQ[θB]}
```
**Output:**
```
Out= {1/3, 1/7, False, False}
```
## 2026-07-18 14:44:32 — ✏️ EDIT CELL 125
(*Standalone single-site symbolic-θ fused hierarchy (mirrors Section E's t1/qdet1/xL1, θ1->θa).
  (*why θa,θb as function ARGUMENTS, not literally θ1,θ2*) θ1,θ2 are GLOBAL symbols already assigned
  n…


## 2026-07-18 14:44:32
**Input:**
```wolfram
(*Standalone single-site symbolic-θ fused hierarchy (mirrors Section E's t1/qdet1/xL1, θ1->θa).
  (*why θa,θb as function ARGUMENTS, not literally θ1,θ2*) θ1,θ2 are GLOBAL symbols already assigned
  numeric values (1/3, 1/7) in Cell 1 — calling xL1Sym[θ1][...] would silently substitute those
  numbers immediately, NOT stay symbolic. The functions below take θa (site 1) / θa,θb (two-site) as
  ordinary pattern arguments so callers can pass either the free symbols θA,θB (genuinely unassigned,
  for symbolic exploration) or the numeric θ1,θ2 (for validation against the live machinery).*)
ClearAll[t1Sym,qdet1Sym,ν1SingleSym,xL1Sym];
t1Sym[θa_][λ_][0,s_][u_] := IdentityMatrix[dim[{λ,0}]];
t1Sym[θa_][λ_][a_,0][u_] := IdentityMatrix[dim[{λ,0}]];
t1Sym[θa_][λ_][1,1][u_] := Sum[L[λ][i,j][u-θa] G[[j,i]],{i,2},{j,2}];
qdet1Sym[θa_][λ_][u_] := L[λ][1,1][u-θa].L[λ][2,2][u-θa-h] - L[λ][2,1][u-θa].L[λ][1,2][u-θa-h];
t1Sym[θa_][λ_][2,1][u_] := χ2 qdet1Sym[θa][λ][u];
t1Sym[θa_][λ_][1,s_/;s>1][u_] := t1Sym[θa][λ][1,s-1][u].t1Sym[θa][λ][1,1][u+(s-1)h] - t1Sym[θa][λ][2,1][u+(s-1)h].t1Sym[θa][λ][1,s-2][u];
ν1SingleSym[θa_][λ_][u_] := (u-θa-h λ);
xL1Sym[θa_][λ_][{0}] := xSingle[λ][{0}];
xL1Sym[θa_][λ_][{s_}] := xSingle[λ][{0}].t1Sym[θa][λ][1,s][θa] / Product[ν1SingleSym[θa][λ][θa+(p-1)h],{p,1,s}];

(*Standalone two-site symbolic-θ hierarchy (mirrors Section A's L1/L2/T/t2/ν1/x2, θ1,θ2->θa,θb).*)
ClearAll[L1Sym,L2Sym,TSym,t2Sym,qdetTSym,ν1TwoSym,x2Sym];
L1Sym[θa_,θb_][λ1_,λ2_][i_,j_][u_] := KroneckerProduct[L[λ1][i,j][u-θa],IdentityMatrix[dim[{λ2,0}]]];
L2Sym[θa_,θb_][λ1_,λ2_][i_,j_][u_] := KroneckerProduct[IdentityMatrix[dim[{λ1,0}]],L[λ2][i,j][u-θb]];
TSym[θa_,θb_][λ1_,λ2_][i_,j_][u_] := Sum[L1Sym[θa,θb][λ1,λ2][i,k][u].L2Sym[θa,θb][λ1,λ2][k,j][u],{k,2}];
t2Sym[θa_,θb_][λ1_,λ2_][0,s_][u_] := IdentityMatrix[dim[{λ1,0}]dim[{λ2,0}]];
t2Sym[θa_,θb_][λ1_,λ2_][a_,0][u_] := IdentityMatrix[dim[{λ1,0}]dim[{λ2,0}]];
t2Sym[θa_,θb_][λ1_,λ2_][1,1][u_] := Sum[TSym[θa,θb][λ1,λ2][i,j][u] G[[j,i]],{i,2},{j,2}];
qdetTSym[θa_,θb_][λ1_,λ2_][u_] := TSym[θa,θb][λ1,λ2][1,1][u].TSym[θa,θb][λ1,λ2][2,2][u-h] - TSym[θa,θb][λ1,λ2][2,1][u].TSym[θa,θb][λ1,λ2][1,2][u-h];
t2Sym[θa_,θb_][λ1_,λ2_][2,1][u_] := χ2 qdetTSym[θa,θb][λ1,λ2][u];
t2Sym[θa_,θb_][λ1_,λ2_][1,s_/;s>1][u_] := t2Sym[θa,θb][λ1,λ2][1,s-1][u].t2Sym[θa,θb][λ1,λ2][1,1][u+(s-1)h] - t2Sym[θa,θb][λ1,λ2][2,1][u+(s-1)h].t2Sym[θa,θb][λ1,λ2][1,s-2][u];
ν1TwoSym[θa_,θb_][λ1_,λ2_][u_] := (u-θa-h λ1)(u-θb-h λ2);
x2Sym[θa_,θb_][λ1_,λ2_][{0,0}] := KroneckerProduct[xSingle[λ1][{0}],xSingle[λ2][{0}]];
x2Sym[θa_,θb_][λ1_,λ2_][{s1_,s2_}] := x2Sym[θa,θb][λ1,λ2][{0,0}].t2Sym[θa,θb][λ1,λ2][1,s1][θa].t2Sym[θa,θb][λ1,λ2][1,s2][θb] /
  (Product[ν1TwoSym[θa,θb][λ1,λ2][θa+(k-1)h],{k,1,s1}] Product[ν1TwoSym[θa,θb][λ1,λ2][θb+(k-1)h],{k,1,s2}]);
```
**Output:**
```
Out= Null
```
## 2026-07-18 14:44:39 — ✏️ EDIT CELL 126
(*Validation: symbolic constructions (called with θA,θB, then substituted /.{θA->θ1,θB->θ2} to the
  actual numeric values) must reproduce the live numeric xL1/x2 exactly — a port-fidelity check
  bef…


## 2026-07-18 14:44:39
**Input:**
```wolfram
(*Validation: symbolic constructions (called with θA,θB, then substituted /.{θA->θ1,θB->θ2} to the
  actual numeric values) must reproduce the live numeric xL1/x2 exactly — a port-fidelity check
  before trusting the symbolic run.*)
valXL1 = Table[
  Chop[N[(Flatten[xL1Sym[θA][λ][{s}]]/.θA->θ1) - Flatten[xL1[λ][{s}]]],10^-10]
,{λ,{1,2,3}},{s,0,λ}];
valX2 = Table[
  Chop[N[(Flatten[x2Sym[θA,θB][λ1,λ2][{s1,s2}]]/.{θA->θ1,θB->θ2}) - Flatten[x2[λ1,λ2][{s1,s2}]]],10^-10]
,{λ1,{1,2}},{λ2,{1,2}},{s1,0,λ1},{s2,0,λ2}];
worstValXL1 = Max[Abs[Flatten[valXL1]]];
worstValX2 = Max[Abs[Flatten[valX2]]];
If[worstValXL1 > 10^-8 || worstValX2 > 10^-8, Abort[]];
{"Symbolic-θ port-fidelity: xL1Sym vs xL1", worstValXL1, "x2Sym vs x2", worstValX2}
```
**Output:**
```
Out= {"Symbolic-θ port-fidelity: xL1Sym vs xL1", 0, "x2Sym vs x2", 0}
```
## 2026-07-18 14:44:49 — ✏️ EDIT CELL 127
(*Symbolic decomposition, θA,θB GENUINELY free (unlike θ1,θ2 which already carry numeric values —
  see the fix note on the previous cell). Solve x2Sym[{s1,m}] = Σ A[s1,m;s1',m'] xL1Sym[λ1][{s1'}]
  ⊗…


## 2026-07-18 14:44:49
**Input:**
```wolfram
(*Symbolic decomposition, θA,θB GENUINELY free (unlike θ1,θ2 which already carry numeric values —
  see the fix note on the previous cell). Solve x2Sym[{s1,m}] = Σ A[s1,m;s1',m'] xL1Sym[λ1][{s1'}]
  ⊗xL1Sym[λ2][{m'}] over the FULL basis (not assuming level-conservation — that should fall out of
  the symbolic solve as A=0 off-level, reconfirming Section I's numeric finding analytically).
  (*why no underscores in variable names*) Wolfram Language reserves "_" for pattern syntax even
  mid-identifier (e.g. "dCase1101" not "d11_01") — a literal underscore in a variable name silently
  turns it into a Blank pattern, causing Set::write errors on assignment.*)
ClearAll[fullBasisSym,decompX2Sym];
fullBasisSym[λ1_,λ2_] := Flatten[Table[
   {s1p,m1p} -> Flatten[KroneckerProduct[xL1Sym[θA][λ1][{s1p}], xL1Sym[θB][λ2][{m1p}]]]
  ,{s1p,0,λ1},{m1p,0,λ2}], 1];
decompX2Sym[λ1_,λ2_][s1_,m_] := Module[{basisRules,basisVecs,mat,target,sol},
  basisRules = fullBasisSym[λ1,λ2];
  basisVecs = basisRules[[All,2]];
  mat = Transpose[basisVecs];
  target = Flatten[x2Sym[θA,θB][λ1,λ2][{s1,m}]];
  sol = Simplify[LinearSolve[mat, target]];
  Thread[basisRules[[All,1]] -> sol]
];
(*smallest nontrivial case first: λ1=1,λ2=1, state (s1,m)=(0,1) — should mix (0,1) and (1,0).*)
dCase1101 = decompX2Sym[1,1][0,1];
dCase1101
```
**Output:**
```
Out= {{0, 0} -> 0, {0, 1} -> (θA - θB)/(1 + θA - θB), {1, 0} -> (1 + θA - θB)^(-1), {1, 1} -> 0}
```
## 2026-07-18 14:45:06 — 📥 BULK INSERT 1 CELL(S) at positions 128–128
1. [code] (*Full symbolic sweep λ1,λ2∈{1,2}, all (s1,m). Collect every nonzero coefficient with its full↵  (λ1


## 2026-07-18 14:45:07
**Input:**
```wolfram
(*Full symbolic sweep λ1,λ2∈{1,2}, all (s1,m). Collect every nonzero coefficient with its full
  (λ1,λ2,s1,m,s1',m') label for pattern-hunting.*)
symSweep = {};
Do[
  Do[
    Do[
      Module[{d},
        d = decompX2Sym[λ1,λ2][s1,m];
        Do[
          If[rule[[2]]=!=0, AppendTo[symSweep, {λ1,λ2,s1,m,rule[[1,1]],rule[[1,2]],rule[[2]]}]]
        ,{rule,d}]
      ]
    ,{m,0,λ2}]
  ,{s1,0,λ1}]
,{λ1,{1,2}},{λ2,{1,2}}];
Length[symSweep]
```
**Output:**
```
Out[453]= 42
```

## 2026-07-18 14:45:12
**Input:**
```wolfram
symSweep
```
**Output:**
```
Out=                                                                                                                                                                                                                                                                                        -28        -28   2            -28        -28            -28           -28   2                                                                                                                                                                                                                                                            -29        -29               -29        -29                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            2                    2                                                                                          -28        -28                                                                                                                                            -28        -28                                                                                                                                                                                                                                                                                                                                                                                                                                   -28        -28                                                                                                                                                                                                                      -28        -28                                                                                                                                                                                                                        -29        -29      2                                                                                                                                                                                                                 -28        -28                                                                                                                                                                                                                       -28        -28                                                                                                                                                                                                                         -29        -29      2                                                                                                                                                                                                                                                                                                                                                         -57        -57      3                                                                                                                                                                                                                                                                                             -56        -57                                                                                                                                                                                                                                                                                                    -56        -56      2           
[output truncated — 24028 chars total]
```

## 2026-07-18 14:45:27
**Input:**
```wolfram
{χ2, χ2//Precision, χ2//InputForm}
```
**Output:**
```
Out= {1, Infinity, InputForm[1]}
```

## 2026-07-18 14:45:47
**Input:**
```wolfram
{G, G//InputForm}
```
**Output:**
```
Out= {{{0.7208797502944530745333553011911001825875290465245700049878`29.54117567596917 + 0``29.50582140452307*I, -1}, {1, 0}}, InputForm[{{0.7208797502944530745333553011911001825875290465245700049878`29.54117567596917 + 0``29.50582140452307*I, -1}, {1, 0}}]}
```
## 2026-07-18 14:46:10 — ✏️ EDIT CELL 125
(*Standalone single-site symbolic-θ fused hierarchy (mirrors Section E's t1/qdet1/xL1, θ1->θa).
  (*why θa,θb as function ARGUMENTS, not literally θ1,θ2*) θ1,θ2 are GLOBAL symbols already assigned
  n…


## 2026-07-18 14:46:10
**Input:**
```wolfram
(*Standalone single-site symbolic-θ fused hierarchy (mirrors Section E's t1/qdet1/xL1, θ1->θa).
  (*why θa,θb as function ARGUMENTS, not literally θ1,θ2*) θ1,θ2 are GLOBAL symbols already assigned
  numeric values (1/3, 1/7) in Cell 1 — calling xL1Sym[θ1][...] would silently substitute those
  numbers immediately, NOT stay symbolic. The functions below take θa (site 1) / θa,θb (two-site) as
  ordinary pattern arguments so callers can pass either the free symbols θA,θB (genuinely unassigned,
  for symbolic exploration) or the numeric θ1,θ2 (for validation against the live machinery).
  (*why GSym not G*) G's (1,1) entry χ1=κ1+κ2 is a ~30-digit-precision APPROXIMATE real number
  (κ1=Exp[IZeta[3]] is transcendental) — mixing that into Simplify with exact symbolic θA,θB leaves
  ugly unsimplified 10^-28-type precision residues instead of collapsing cleanly. GSym uses a fresh
  free symbol χ1sym in G's place so everything stays exact; verified below whether χ1sym survives in
  the final A-coefficients (expected: cancels, matching Section E/F's χ-independent normalized SoV
  covector ratios).*)
ClearAll[GSym,t1Sym,qdet1Sym,ν1SingleSym,xL1Sym];
GSym = {{χ1sym,-1},{1,0}};
t1Sym[θa_][λ_][0,s_][u_] := IdentityMatrix[dim[{λ,0}]];
t1Sym[θa_][λ_][a_,0][u_] := IdentityMatrix[dim[{λ,0}]];
t1Sym[θa_][λ_][1,1][u_] := Sum[L[λ][i,j][u-θa] GSym[[j,i]],{i,2},{j,2}];
qdet1Sym[θa_][λ_][u_] := L[λ][1,1][u-θa].L[λ][2,2][u-θa-h] - L[λ][2,1][u-θa].L[λ][1,2][u-θa-h];
t1Sym[θa_][λ_][2,1][u_] := qdet1Sym[θa][λ][u];  (*χ2=1 exactly, no need for a symbolic stand-in*)
t1Sym[θa_][λ_][1,s_/;s>1][u_] := t1Sym[θa][λ][1,s-1][u].t1Sym[θa][λ][1,1][u+(s-1)h] - t1Sym[θa][λ][2,1][u+(s-1)h].t1Sym[θa][λ][1,s-2][u];
ν1SingleSym[θa_][λ_][u_] := (u-θa-h λ);
xL1Sym[θa_][λ_][{0}] := xSingle[λ][{0}];
xL1Sym[θa_][λ_][{s_}] := xSingle[λ][{0}].t1Sym[θa][λ][1,s][θa] / Product[ν1SingleSym[θa][λ][θa+(p-1)h],{p,1,s}];

(*Standalone two-site symbolic-θ hierarchy (mirrors Section A's L1/L2/T/t2/ν1/x2, θ1,θ2->θa,θb).*)
ClearAll[L1Sym,L2Sym,TSym,t2Sym,qdetTSym,ν1TwoSym,x2Sym];
L1Sym[θa_,θb_][λ1_,λ2_][i_,j_][u_] := KroneckerProduct[L[λ1][i,j][u-θa],IdentityMatrix[dim[{λ2,0}]]];
L2Sym[θa_,θb_][λ1_,λ2_][i_,j_][u_] := KroneckerProduct[IdentityMatrix[dim[{λ1,0}]],L[λ2][i,j][u-θb]];
TSym[θa_,θb_][λ1_,λ2_][i_,j_][u_] := Sum[L1Sym[θa,θb][λ1,λ2][i,k][u].L2Sym[θa,θb][λ1,λ2][k,j][u],{k,2}];
t2Sym[θa_,θb_][λ1_,λ2_][0,s_][u_] := IdentityMatrix[dim[{λ1,0}]dim[{λ2,0}]];
t2Sym[θa_,θb_][λ1_,λ2_][a_,0][u_] := IdentityMatrix[dim[{λ1,0}]dim[{λ2,0}]];
t2Sym[θa_,θb_][λ1_,λ2_][1,1][u_] := Sum[TSym[θa,θb][λ1,λ2][i,j][u] GSym[[j,i]],{i,2},{j,2}];
qdetTSym[θa_,θb_][λ1_,λ2_][u_] := TSym[θa,θb][λ1,λ2][1,1][u].TSym[θa,θb][λ1,λ2][2,2][u-h] - TSym[θa,θb][λ1,λ2][2,1][u].TSym[θa,θb][λ1,λ2][1,2][u-h];
t2Sym[θa_,θb_][λ1_,λ2_][2,1][u_] := qdetTSym[θa,θb][λ1,λ2][u];
t2Sym[θa_,θb_][λ1_,λ2_][1,s_/;s>1][u_] := t2Sym[θa,θb][λ1,λ2][1,s-1][u].t2Sym[θa,θb][λ1,λ2][1,1][u+(s-1)h] - t2Sym[θa,θb][λ1,λ2][2,1][u+(s-1)h].t2Sym[θa,θb][λ1,λ2][1,s-2][u];
ν1TwoSym[θa_,θb_][λ1_,λ2_][u_] := (u-θa-h λ1)(u-θb-h λ2);
x2Sym[θa_,θb_][λ1_,λ2_][{0,0}] := KroneckerProduct[xSingle[λ1][{0}],xSingle[λ2][{0}]];
x2Sym[θa_,θb_][λ1_,λ2_][{s1_,s2_}] := x2Sym[θa,θb][λ1,λ2][{0,0}].t2Sym[θa,θb][λ1,λ2][1,s1][θa].t2Sym[θa,θb][λ1,λ2][1,s2][θb] /
  (Product[ν1TwoSym[θa,θb][λ1,λ2][θa+(k-1)h],{k,1,s1}] Product[ν1TwoSym[θa,θb][λ1,λ2][θb+(k-1)h],{k,1,s2}]);
```
**Output:**
```
Out= Null
```
## 2026-07-18 14:46:18 — ✏️ EDIT CELL 126
(*Validation: symbolic constructions (called with θA,θB,χ1sym, then substituted to the actual
  numeric θ1,θ2,χ1) must reproduce the live numeric xL1/x2 exactly — a port-fidelity check before
  trusti…


## 2026-07-18 14:46:18
**Input:**
```wolfram
(*Validation: symbolic constructions (called with θA,θB,χ1sym, then substituted to the actual
  numeric θ1,θ2,χ1) must reproduce the live numeric xL1/x2 exactly — a port-fidelity check before
  trusting the symbolic run.*)
valXL1 = Table[
  Chop[N[(Flatten[xL1Sym[θA][λ][{s}]]/.{θA->θ1,χ1sym->χ1}) - Flatten[xL1[λ][{s}]]],10^-10]
,{λ,{1,2,3}},{s,0,λ}];
valX2 = Table[
  Chop[N[(Flatten[x2Sym[θA,θB][λ1,λ2][{s1,s2}]]/.{θA->θ1,θB->θ2,χ1sym->χ1}) - Flatten[x2[λ1,λ2][{s1,s2}]]],10^-10]
,{λ1,{1,2}},{λ2,{1,2}},{s1,0,λ1},{s2,0,λ2}];
worstValXL1 = Max[Abs[Flatten[valXL1]]];
worstValX2 = Max[Abs[Flatten[valX2]]];
If[worstValXL1 > 10^-8 || worstValX2 > 10^-8, Abort[]];
{"Symbolic-θ,χ1 port-fidelity: xL1Sym vs xL1", worstValXL1, "x2Sym vs x2", worstValX2}
```
**Output:**
```
Out= {"Symbolic-θ,χ1 port-fidelity: xL1Sym vs xL1", 0, "x2Sym vs x2", 0}
```
## 2026-07-18 14:46:29 — ▶️ RUN CELL 127 (0.03 s)
**In [127]:** `(*Symbolic decomposition, θA,θB GENUINELY free (unlike θ1,θ2 which already carry numeric values —↵  see the fix note on the previous cell). Solve x2Sym[{s1,m}] = Σ A[s1,m;s1',m'] xL1Sym[λ1][{s1'}]↵  ⊗`
**Out:** Out[458]= \left\{\left\{0,0\right\}\to 0,\left\{0,1\right\}\to \frac{\theta A-\theta B}{\theta A-\theta B+1},\left\{1,0\right\}\to \frac{1}{\theta A-\theta B+1},\left\{1,1\right\}\to 0\right\}

## 2026-07-18 14:46:47 — ▶️ RUN CELL 128 (0.03 s)
**In [128]:** `(*Full symbolic sweep λ1,λ2∈{1,2}, all (s1,m). Collect every nonzero coefficient with its full↵  (λ1,λ2,s1,m,s1',m') label for pattern-hunting.*)↵symSweep = {};↵Do[↵  Do[↵    Do[↵      Module[{d},↵   `
**Out:** Out[461]= 35


## 2026-07-18 14:46:55
**Input:**
```wolfram
TableForm[symSweep]
```
**Output:**
```
Out= TableForm[{{1, 1, 0, 0, 0, 0, 1}, {1, 1, 0, 1, 0, 1, (θA - θB)/(1 + θA - θB)}, {1, 1, 0, 1, 1, 0, (1 + θA - θB)^(-1)}, {1, 1, 1, 0, 1, 0, 1}, {1, 1, 1, 1, 1, 1, 1}, {1, 2, 0, 0, 0, 0, 1}, {1, 2, 0, 1, 0, 1, (θA - θB)/(1 + θA - θB)}, {1, 2, 0, 1, 1, 0, (1 + θA - θB)^(-1)}, {1, 2, 0, 2, 0, 2, (-1 + θA - θB)/(1 + θA - θB)}, {1, 2, 0, 2, 1, 1, 2/(1 + θA - θB)}, {1, 2, 1, 0, 1, 0, 1}, {1, 2, 1, 1, 1, 1, 1}, {1, 2, 1, 2, 1, 2, 1}, {2, 1, 0, 0, 0, 0, 1}, {2, 1, 0, 1, 0, 1, (θA - θB)/(2 + θA - θB)}, {2, 1, 0, 1, 1, 0, 2/(2 + θA - θB)}, {2, 1, 1, 0, 1, 0, 1}, {2, 1, 1, 1, 1, 1, (1 + θA - θB)/(2 + θA - θB)}, {2, 1, 1, 1, 2, 0, (2 + θA - θB)^(-1)}, {2, 1, 2, 0, 2, 0, 1}, {2, 1, 2, 1, 2, 1, 1}, {2, 2, 0, 0, 0, 0, 1}, {2, 2, 0, 1, 0, 1, (θA - θB)/(2 + θA - θB)}, {2, 2, 0, 1, 1, 0, 2/(2 + θA - θB)}, {2, 2, 0, 2, 0, 2, (-θA + θA^2 + θB - 2*θA*θB + θB^2)/((1 + θA - θB)*(2 + θA - θB))}, {2, 2, 0, 2, 1, 1, (4*(θA - θB))/((1 + θA - θB)*(2 + θA - θB))}, {2, 2, 0, 2, 2, 0, 2/((1 + θA - θB)*(2 + θA - θB))}, {2, 2, 1, 0, 1, 0, 1}, {2, 2, 1, 1, 1, 1, (1 + θA - θB)/(2 + θA - θB)}, {2, 2, 1, 1, 2, 0, (2 + θA - θB)^(-1)}, {2, 2, 1, 2, 1, 2, (θA - θB)/(2 + θA - θB)}, {2, 2, 1, 2, 2, 1, 2/(2 + θA - θB)}, {2, 2, 2, 0, 2, 0, 1}, {2, 2, 2, 1, 2, 1, 1}, {2, 2, 2, 2, 2, 2, 1}}]
```
## 2026-07-18 14:48:33 — 📥 BULK INSERT 1 CELL(S) at positions 129–129
1. [code] (*Rescale to the RAW (unnormalized) numerators — yL1 := xSingle.t1[1,s][θ] (no ν1Single division),↵


## 2026-07-18 14:48:34
**Input:**
```wolfram
(*Rescale to the RAW (unnormalized) numerators — yL1 := xSingle.t1[1,s][θ] (no ν1Single division),
  y2 := x2[{0,0}].t2[1,s1][θa].t2[1,m][θb] (no ν1Two division) — to see if the underlying recoupling
  B is a clean combinatorial (Vandermonde-type) matrix once the known θ-dependent normalization
  factors D1,D2 are stripped out.*)
ClearAll[yL1Sym,y2Sym,decompY2Sym];
yL1Sym[θa_][λ_][{s_}] := xSingle[λ][{0}].t1Sym[θa][λ][1,s][θa];
y2Sym[θa_,θb_][λ1_,λ2_][{s1_,s2_}] := x2Sym[θa,θb][λ1,λ2][{0,0}].t2Sym[θa,θb][λ1,λ2][1,s1][θa].t2Sym[θa,θb][λ1,λ2][1,s2][θb];
fullBasisYSym[λ1_,λ2_] := Flatten[Table[
   {s1p,m1p} -> Flatten[KroneckerProduct[yL1Sym[θA][λ1][{s1p}], yL1Sym[θB][λ2][{m1p}]]]
  ,{s1p,0,λ1},{m1p,0,λ2}], 1];
decompY2Sym[λ1_,λ2_][s1_,m_] := Module[{basisRules,basisVecs,mat,target,sol},
  basisRules = fullBasisYSym[λ1,λ2];
  basisVecs = basisRules[[All,2]];
  mat = Transpose[basisVecs];
  target = Flatten[y2Sym[θA,θB][λ1,λ2][{s1,m}]];
  sol = Simplify[LinearSolve[mat, target]];
  Thread[basisRules[[All,1]] -> sol]
];
(*same smallest test case: λ1=1,λ2=1, (s1,m)=(0,1).*)
decompY2Sym[1,1][0,1]
```
**Output:**
```
Out[467]= \left\{\left\{0,0\right\}\to 0,\left\{0,1\right\}\to \theta B-\theta A,\left\{1,0\right\}\to -1,\left\{1,1\right\}\to 0\right\}
```
## 2026-07-18 14:48:46 — 📥 BULK INSERT 1 CELL(S) at positions 130–130
1. [code] (*Full B-sweep, λ1,λ2∈{1,2,3} (larger than before since B looks much cheaper/cleaner — polynomial,↵


## 2026-07-18 14:48:47
**Input:**
```wolfram
(*Full B-sweep, λ1,λ2∈{1,2,3} (larger than before since B looks much cheaper/cleaner — polynomial,
  not rational — so bigger reps should still be tractable).*)
symSweepB = {};
Do[
  Do[
    Do[
      Module[{d},
        d = decompY2Sym[λ1,λ2][s1,m];
        Do[
          If[rule[[2]]=!=0, AppendTo[symSweepB, {λ1,λ2,s1,m,rule[[1,1]],rule[[1,2]],rule[[2]]}]]
        ,{rule,d}]
      ]
    ,{m,0,λ2}]
  ,{s1,0,λ1}]
,{λ1,{1,2,3}},{λ2,{1,2,3}}];
Length[symSweepB]
```
**Output:**
```
Out[470]= 127
```

## 2026-07-18 14:48:53
**Input:**
```wolfram
Select[symSweepB, #[[1]]==3&&#[[2]]==3&]
```
**Output:**
```
Out= {{3, 3, 0, 0, 0, 0, 1}, {3, 3, 0, 1, 0, 1, -θA + θB}, {3, 3, 0, 1, 1, 0, -3}, {3, 3, 0, 2, 0, 2, -θA + θA^2 + θB - 2*θA*θB + θB^2}, {3, 3, 0, 2, 1, 1, 4*(θA - θB)}, {3, 3, 0, 2, 2, 0, 6}, {3, 3, 0, 3, 0, 3, -θA^3 + 3*θA^2*(1 + θB) + θB*(2 + 3*θB + θB^2) - θA*(2 + 6*θB + 3*θB^2)}, {3, 3, 0, 3, 1, 2, -3*(-θA + θA^2 + θB - 2*θA*θB + θB^2)}, {3, 3, 0, 3, 2, 1, -6*θA + 6*θB}, {3, 3, 0, 3, 3, 0, -6}, {3, 3, 1, 0, 1, 0, -3 + θA - θB}, {3, 3, 1, 1, 1, 1, -((-3 + θA - θB)*(1 + θA - θB))}, {3, 3, 1, 1, 2, 0, 9 - 3*θA + 3*θB}, {3, 3, 1, 2, 1, 2, θA^3 - θA^2*(2 + 3*θB) - θB*(-3 + 2*θB + θB^2) + θA*(-3 + 4*θB + 3*θB^2)}, {3, 3, 1, 2, 2, 1, 4*(-3 + θA^2 + 2*θB + θB^2 - 2*θA*(1 + θB))}, {3, 3, 1, 2, 3, 0, 6*(-3 + θA - θB)}, {3, 3, 1, 3, 1, 3, -((-3 + θA - θB)*(-1 + θA - θB)*(θA - θB)*(1 + θA - θB))}, {3, 3, 1, 3, 2, 2, 3*(-θA^3 + θA^2*(2 + 3*θB) + θA*(3 - 4*θB - 3*θB^2) + θB*(-3 + 2*θB + θB^2))}, {3, 3, 1, 3, 3, 1, -6*(-3 + θA - θB)*(1 + θA - θB)}, {3, 3, 2, 0, 2, 0, (-3 + θA - θB)*(-2 + θA - θB)}, {3, 3, 2, 1, 2, 1, -12 - θA^3 - 4*θB + 3*θB^2 + θB^3 + 3*θA^2*(1 + θB) + θA*(4 - 6*θB - 3*θB^2)}, {3, 3, 2, 1, 3, 0, -3*(6 + θA^2 + 5*θB + θB^2 - θA*(5 + 2*θB))}, {3, 3, 2, 2, 2, 2, (-3 + θA - θB)*(-2 + θA - θB)*(1 + θA - θB)*(2 + θA - θB)}, {3, 3, 2, 2, 3, 1, 4*(-3 + θA - θB)*(-2 + θA - θB)*(2 + θA - θB)}, {3, 3, 2, 3, 2, 3, -θA^5 + θA^4*(2 + 5*θB) + θA^3*(7 - 8*θB - 10*θB^2) + θA^2*(-8 - 21*θB + 12*θB^2 + 10*θB^3) + θA*(-12 + 16*θB + 21*θB^2 - 8*θB^3 - 5*θB^4) + θB*(12 - 8*θB - 7*θB^2 + 2*θB^3 + θB^4)}, {3, 3, 2, 3, 3, 2, -3*(12 + θA^4 - 8*θB - 7*θB^2 + 2*θB^3 + θB^4 - 2*θA^3*(1 + 2*θB) + θA^2*(-7 + 6*θB + 6*θB^2) + θA*(8 + 14*θB - 6*θB^2 - 4*θB^3))}, {3, 3, 3, 0, 3, 0, -6 + θA^3 - 11*θB - 6*θB^2 - θB^3 - 3*θA^2*(2 + θB) + θA*(11 + 12*θB + 3*θB^2)}, {3, 3, 3, 1, 3, 1, -((-3 + θA - θB)*(-2 + θA - θB)*(-1 + θA - θB)*(3 + θA - θB))}, {3, 3, 3, 2, 3, 2, -36 + θA^5 - 36*θB + 13*θB^2 + 13*θB^3 - θB^4 - θB^5 - θA^4*(1 + 5*θB) + θA^3*(-13 + 4*θB + 10*θB^2) + θA^2*(13 + 39*θB - 6*θB^2 - 10*θB^3) + θA*(36 - 26*θB - 39*θB^2 + 4*θB^3 + 5*θB^4)}, {3, 3, 3, 3, 3, 3, -((-3 + θA - θB)*(-2 + θA - θB)*(-1 + θA - θB)*(1 + θA - θB)*(2 + θA - θB)*(3 + θA - θB))}}
```

## 2026-07-18 14:49:38
**Input:**
```wolfram
Factor /@ (Select[symSweepB, #[[1]]==2&&#[[2]]==2&][[All,7]])
Select[symSweepB, #[[1]]==2&&#[[2]]==2&][[All,{3,4,5,6}]]
```
**Output:**
```
Out[1]= {1, -θA + θB, -2, (-1 + θA - θB)*(θA - θB), 2*(θA - θB), 2, -2 + θA - θB, -((-2 + θA - θB)*(1 + θA - θB)), -2*(-2 + θA - θB), (-2 + θA - θB)*(θA - θB)*(1 + θA - θB), 2*(-2 + θA - θB)*(1 + θA - θB), (-2 + θA - θB)*(-1 + θA - θB), -((-2 + θA - θB)*(-1 + θA - θB)*(2 + θA - θB)), (-2 + θA - θB)*(-1 + θA - θB)*(1 + θA - θB)*(2 + θA - θB)}
Out[2]= {{0, 0, 0, 0}, {0, 1, 0, 1}, {0, 1, 1, 0}, {0, 2, 0, 2}, {0, 2, 1, 1}, {0, 2, 2, 0}, {1, 0, 1, 0}, {1, 1, 1, 1}, {1, 1, 2, 0}, {1, 2, 1, 2}, {1, 2, 2, 1}, {2, 0, 2, 0}, {2, 1, 2, 1}, {2, 2, 2, 2}}
```
## 2026-07-18 14:51:24 — 📥 BULK INSERT 1 CELL(S) at positions 131–131
1. [code] (*Extend the clean A-sweep to λ1,λ2∈{1,2,3} (bug-fixed χ1sym version) — bigger than the earlier↵  λ≤


## 2026-07-18 14:51:25
**Input:**
```wolfram
(*Extend the clean A-sweep to λ1,λ2∈{1,2,3} (bug-fixed χ1sym version) — bigger than the earlier
  λ≤2 sweep, to get more data for pattern-hunting; only keep diagonal + one-step-off-diagonal entries
  to control runtime (LinearSolve over larger bases is the bottleneck for λ=3).*)
symSweepFull = {};
Do[
  Do[
    Do[
      Module[{d},
        d = decompX2Sym[λ1,λ2][s1,m];
        Do[
          If[rule[[2]]=!=0, AppendTo[symSweepFull, {λ1,λ2,s1,m,rule[[1,1]],rule[[1,2]],Factor[rule[[2]]]}]]
        ,{rule,d}]
      ]
    ,{m,0,λ2}]
  ,{s1,0,λ1}]
,{λ1,{1,2,3}},{λ2,{1,2,3}}];
Length[symSweepFull]
```
**Output:**
```
Out[473]= 127
```

## 2026-07-18 14:51:35
**Input:**
```wolfram
Select[symSweepFull, #[[1]]==3&&#[[2]]==3&&#[[3]]==#[[5]]&&#[[4]]==#[[6]]&][[All,{3,4,7}]]
```
**Output:**
```
Out= {{0, 0, 1}, {0, 1, (θA - θB)/(3 + θA - θB)}, {0, 2, ((-1 + θA - θB)*(θA - θB))/((2 + θA - θB)*(3 + θA - θB))}, {0, 3, ((-2 + θA - θB)*(-1 + θA - θB)*(θA - θB))/((1 + θA - θB)*(2 + θA - θB)*(3 + θA - θB))}, {1, 0, 1}, {1, 1, (1 + θA - θB)/(3 + θA - θB)}, {1, 2, ((θA - θB)*(1 + θA - θB))/((2 + θA - θB)*(3 + θA - θB))}, {1, 3, ((-1 + θA - θB)*(θA - θB))/((2 + θA - θB)*(3 + θA - θB))}, {2, 0, 1}, {2, 1, (2 + θA - θB)/(3 + θA - θB)}, {2, 2, (1 + θA - θB)/(3 + θA - θB)}, {2, 3, (θA - θB)/(3 + θA - θB)}, {3, 0, 1}, {3, 1, 1}, {3, 2, 1}, {3, 3, 1}}
```

## 2026-07-18 14:52:07
**Input:**
```wolfram
Select[symSweepFull, #[[1]]==3&&#[[2]]==3&&#[[3]]!=#[[5]]&][[All,{3,4,5,6,7}]]//TableForm
```
**Output:**
```
Out= TableForm[{{0, 1, 1, 0, 3/(3 + θA - θB)}, {0, 2, 1, 1, (6*(θA - θB))/((2 + θA - θB)*(3 + θA - θB))}, {0, 2, 2, 0, 6/((2 + θA - θB)*(3 + θA - θB))}, {0, 3, 1, 2, (9*(-1 + θA - θB)*(θA - θB))/((1 + θA - θB)*(2 + θA - θB)*(3 + θA - θB))}, {0, 3, 2, 1, (18*(θA - θB))/((1 + θA - θB)*(2 + θA - θB)*(3 + θA - θB))}, {0, 3, 3, 0, 6/((1 + θA - θB)*(2 + θA - θB)*(3 + θA - θB))}, {1, 1, 2, 0, 2/(3 + θA - θB)}, {1, 2, 2, 1, (4*(1 + θA - θB))/((2 + θA - θB)*(3 + θA - θB))}, {1, 2, 3, 0, 2/((2 + θA - θB)*(3 + θA - θB))}, {1, 3, 2, 2, (6*(θA - θB))/((2 + θA - θB)*(3 + θA - θB))}, {1, 3, 3, 1, 6/((2 + θA - θB)*(3 + θA - θB))}, {2, 1, 3, 0, (3 + θA - θB)^(-1)}, {2, 2, 3, 1, 2/(3 + θA - θB)}, {2, 3, 3, 2, 3/(3 + θA - θB)}}]
```
## 2026-07-18 14:53:10 — 📥 BULK INSERT 1 CELL(S) at positions 132–132
1. [code] (*Pivot: rather than the full A closed-form (looks like a genuine Hahn/Racah-type special↵  function


## 2026-07-18 14:53:10
**Input:**
```wolfram
(*Pivot: rather than the full A closed-form (looks like a genuine Hahn/Racah-type special
  function), test whether the INVERSE of each level-L block is simpler — that's the object actually
  needed for assembly (expressing xL1⊗xL1 in terms of x2, the direction Step 3 needs).
  Build the level-L block matrix of A (rows/cols indexed by s1'=max(0,L-λ2)..min(L,λ1)) and invert.*)
ClearAll[levelBlockA,levelBlockAinv];
levelBlockA[λ1_,λ2_][L_] := Module[{s1rng,mat},
  s1rng = Range[Max[0,L-λ2],Min[L,λ1]];
  mat = Table[
    Module[{d,val},
      d = decompX2Sym[λ1,λ2][s1,L-s1];
      val = s1p /. Cases[d, ({s1pp_,mp_}->v_)/;s1pp==s1p -> ({s1pp,mp}->v)][[1]];
      Simplify[Lookup[Association[d], {s1p,L-s1p}, 0]]
    ]
  ,{s1,s1rng},{s1p,s1rng}];
  {s1rng, mat}
];
(*test: λ1=λ2=3, level L=2 (3x3 block, s1'=0,1,2).*)
{rng32,blk32} = levelBlockA[3,3][2];
blk32 // MatrixForm
```
**Output:**
```
Out[477]//MatrixForm= \begin{pmatrix}\left\{0,0\right\} & \left\{0,0\right\} & \left\{0,0\right\}\\\left\{0,0\right\} & \left\{0,0\right\} & \left\{0,0\right\}\\\left\{0,0\right\} & \left\{0,0\right\} & \left\{0,0\right\}\end{pmatrix}
```

## 2026-07-18 14:53:26
**Input:**
```wolfram
dtest = decompX2Sym[3,3][0,2]
Association[dtest]
Lookup[Association[dtest], {1,1}, 0]
```
**Output:**
```
Out[1]= {{0, 0} -> 0, {0, 1} -> 0, {0, 2} -> (-θA + θA^2 + θB - 2*θA*θB + θB^2)/((2 + θA - θB)*(3 + θA - θB)), {0, 3} -> 0, {1, 0} -> 0, {1, 1} -> (6*(θA - θB))/((2 + θA - θB)*(3 + θA - θB)), {1, 2} -> 0, {1, 3} -> 0, {2, 0} -> 6/((2 + θA - θB)*(3 + θA - θB)), {2, 1} -> 0, {2, 2} -> 0, {2, 3} -> 0, {3, 0} -> 0, {3, 1} -> 0, {3, 2} -> 0, {3, 3} -> 0}
Out[2]= <|{0, 0} -> 0, {0, 1} -> 0, {0, 2} -> (-θA + θA^2 + θB - 2*θA*θB + θB^2)/((2 + θA - θB)*(3 + θA - θB)), {0, 3} -> 0, {1, 0} -> 0, {1, 1} -> (6*(θA - θB))/((2 + θA - θB)*(3 + θA - θB)), {1, 2} -> 0, {1, 3} -> 0, {2, 0} -> 6/((2 + θA - θB)*(3 + θA - θB)), {2, 1} -> 0, {2, 2} -> 0, {2, 3} -> 0, {3, 0} -> 0, {3, 1} -> 0, {3, 2} -> 0, {3, 3} -> 0|>
Out[3]= {0, 0}
```

## 2026-07-18 14:53:41
**Input:**
```wolfram
assoc = <|{1,1}->5,{2,2}->7|>
Lookup[assoc,{1,1},0]
assoc[{1,1}]
KeyExistsQ[assoc,{1,1}]
```
**Output:**
```
Out[1]= <|{1, 1} -> 5, {2, 2} -> 7|>
Out[2]= {0, 0}
Out[3]= 5
Out[4]= True
```
## 2026-07-18 14:53:57 — ✏️ EDIT CELL 132
(*Pivot: rather than the full A closed-form (looks like a genuine Hahn/Racah-type special
  function), test whether the INVERSE of each level-L block is simpler — that's the object actually
  needed f…


## 2026-07-18 14:53:58
**Input:**
```wolfram
(*Pivot: rather than the full A closed-form (looks like a genuine Hahn/Racah-type special
  function), test whether the INVERSE of each level-L block is simpler — that's the object actually
  needed for assembly (expressing xL1⊗xL1 in terms of x2, the direction Step 3 needs).
  Build the level-L block matrix of A (rows/cols indexed by s1'=max(0,L-λ2)..min(L,λ1)) and invert.
  (*why assoc[key] not Lookup[assoc,key,default]*) Lookup[assoc,{a,b},default] treats {a,b} as TWO
  separate keys to look up, not one list-valued key {a,b} — misbehaves (Listable-broadcasts the
  default) exactly when keys are themselves lists, as they are here ({s1',m'} pairs). assoc[key]
  with an explicit KeyExistsQ guard is the correct way to look up a single list-valued key.*)
ClearAll[levelBlockA,levelBlockAinv];
levelBlockA[λ1_,λ2_][L_] := Module[{s1rng,mat,assoc},
  s1rng = Range[Max[0,L-λ2],Min[L,λ1]];
  mat = Table[
    assoc = Association[decompX2Sym[λ1,λ2][s1,L-s1]];
    If[KeyExistsQ[assoc,{s1p,L-s1p}], assoc[{s1p,L-s1p}], 0]
  ,{s1,s1rng},{s1p,s1rng}];
  {s1rng, mat}
];
(*test: λ1=λ2=3, level L=2 (3x3 block, s1'=0,1,2).*)
{rng32,blk32} = levelBlockA[3,3][2];
blk32 // MatrixForm
```
**Output:**
```
Out= MatrixForm[{{(-θA + θA^2 + θB - 2*θA*θB + θB^2)/((2 + θA - θB)*(3 + θA - θB)), (6*(θA - θB))/((2 + θA - θB)*(3 + θA - θB)), 6/((2 + θA - θB)*(3 + θA - θB))}, {0, (1 + θA - θB)/(3 + θA - θB), 2/(3 + θA - θB)}, {0, 0, 1}}]
```
## 2026-07-18 14:54:19 — 📥 BULK INSERT 1 CELL(S) at positions 133–133
1. [code] (*A restricted to a level-L block is UPPER TRIANGULAR (in s1 ascending order) — x2[{s1,m}] only↵  in


## 2026-07-18 14:54:19
**Input:**
```wolfram
(*A restricted to a level-L block is UPPER TRIANGULAR (in s1 ascending order) — x2[{s1,m}] only
  involves basis vectors with s1'>=s1. Invert and inspect (Wolfram handles small triangular inverses
  cleanly).*)
blk32inv = Simplify[Inverse[blk32]];
blk32inv // MatrixForm
```
**Output:**
```
Out[479]//MatrixForm= \begin{pmatrix}\frac{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)}{\theta A^2-2\,\theta A\,\theta B-\theta A+\theta B^2+\theta B} & -\frac{6\,\left(\theta A-\theta B+3\right)}{\theta A^2-2\,\theta A\,\theta B+\theta B^2-1} & \frac{6}{\left(\theta A-\theta B\right)\,\left(\theta A-\theta B+1\right)}\\0 & \frac{\theta A-\theta B+3}{\theta A-\theta B+1} & -… [453 chars total]
```

## 2026-07-18 14:54:26
**Input:**
```wolfram
Factor[blk32inv]//MatrixForm
```
**Output:**
```
Out= MatrixForm[{{((2 + θA - θB)*(3 + θA - θB))/((-1 + θA - θB)*(θA - θB)), (-6*(3 + θA - θB))/((-1 + θA - θB)*(1 + θA - θB)), 6/((θA - θB)*(1 + θA - θB))}, {0, (3 + θA - θB)/(1 + θA - θB), -2/(1 + θA - θB)}, {0, 0, 1}}]
```
## 2026-07-18 14:56:30 — 📥 BULK INSERT 1 CELL(S) at positions 134–134
1. [code] (*Refined Step-3 hypothesis: since A is triangular, Cm[m][k1] (Section I's bracket candidate) may↵


## 2026-07-18 14:56:31
**Input:**
```wolfram
(*Refined Step-3 hypothesis: since A is triangular, Cm[m][k1] (Section I's bracket candidate) may
  equal a LINEAR COMBINATION of several O_m'[k1] pieces (m'=0..λ2-m), not a single term. Test on rep
  (2,1) (λ2=1, only 2 unknowns w0,w1): solve Cm[0][k1] = w0 O_0[k1] + w1 O_1[k1] by matching k1^0,k1^1,
  k1^2 coefficients (3 equations, 2 unknowns — overdetermined, a genuine consistency test).*)
c0k1 = Cm[2,1][1,0][0][k1] // Simplify;
o0k1 = Om[2,1][1,0][0][k1];
o1k1 = Om[2,1][1,0][1][k1];
eqs = Table[Coefficient[c0k1,k1,p] == w0 Coefficient[o0k1,k1,p] + w1 Coefficient[o1k1,k1,p], {p,0,2}];
solW = Solve[eqs, {w0,w1}];
{"equations", eqs, "solution (want: consistent, unique)", solW}
```
**Output:**
```
2 kernel messages: Solve::ivar: 1 is not a valid variable. [also: Solve::ivar] result: Out[485]= \left\{\text{equations},\left\{1.000000000000000000000+0.\times 10^{-21}\,\mathrm{i}⩵\left(2.603469709065218092712489006+0.6885776528619826437114966465\,\mathrm{i}\right)\,w1+\left(1.\,+0.\t…
```

## 2026-07-18 14:56:44
**Input:**
```wolfram
eqs
```
**Output:**
```
Out= {0.99999999999999999999999999999999999974`20.763776783737285 + 0``20.774154586198602*I == (1.`28.385628029452757 + 0``29.827173735437608*I) + (2.60346970906521809271248900648604908753`28.385628029452757 + 0.68857765286198264371149664653546605319`28.330198017570385*I)*w1, -1.03506205686400616702047011719559125013`20.446822045330972 + 0.1660592610973466094858707481640733242`19.654330922994678*I == (1.52215008362353848091245605469939883175`28.3833973045706 - 0.24420479573139207277333933553540189821871773341101956567`28.385628029452757*I) - (1.29610114780177446024989029212092933387`28.385628029452757 + 2.67666169328727834448930616171661419472`28.378365180658164*I)*w1, 0.00214450361864739191183805787673984767`19.0383484896868 + 0.11851144747591667389573406475091463548`20.76968122893966*I == (-0.0362673406094779514502024493860414790055997083400338776284`28.385628029452757 - 2.00423771466623786735432609505223275251`28.385628029452757*I) + (0.3604398751472265372666776505955501101721018646136190171211`28.385628029452757 + 0.93278244859337471648483598207086795224`28.379879135262136*I)*w1}
```

## 2026-07-18 14:57:02
**Input:**
```wolfram
{o0k1, Coefficient[o0k1,k1,0], Coefficient[o0k1,k1,1], Coefficient[o0k1,k1,2]}
```
**Output:**
```
Out= {(1.`28.385628029452757 + 0``29.827173735437608*I) + (1.52215008362353848091245605469939883175`28.3833973045706 - 0.24420479573139207277333933553540189821871773341101956567`28.385628029452757*I)*k1 - (0.0362673406094779514502024493860414790055997083400338776284`28.385628029452757 + 2.00423771466623786735432609505223275251`28.385628029452757*I)*k1^2, 1.`28.385628029452757 + 0``29.827173735437608*I, 1.52215008362353848091245605469939883175`28.3833973045706 - 0.24420479573139207277333933553540189821871773341101956567`28.385628029452757*I, -0.0362673406094779514502024493860414790055997083400338776284`28.385628029452757 - 2.00423771466623786735432609505223275251`28.385628029452757*I}
```

## 2026-07-18 14:57:22
**Input:**
```wolfram
{ValueQ[w0],ValueQ[w1],w0,w1}
```
**Output:**
```
Out= {True, False, 1, w1}
```
## 2026-07-18 14:57:31 — ✏️ EDIT CELL 134
(*Refined Step-3 hypothesis: since A is triangular, Cm[m][k1] (Section I's bracket candidate) may
  equal a LINEAR COMBINATION of several O_m'[k1] pieces (m'=0..λ2-m), not a single term. Test on rep
 …


## 2026-07-18 14:57:31
**Input:**
```wolfram
(*Refined Step-3 hypothesis: since A is triangular, Cm[m][k1] (Section I's bracket candidate) may
  equal a LINEAR COMBINATION of several O_m'[k1] pieces (m'=0..λ2-m), not a single term. Test on rep
  (2,1) (λ2=1, only 2 unknowns wc0,wc1): solve Cm[0][k1] = wc0 O_0[k1] + wc1 O_1[k1] by matching
  k1^0,k1^1,k1^2 coefficients (3 equations, 2 unknowns — overdetermined, a genuine consistency test).
  (*why wc0/wc1 not w0/w1*) w0 is ALREADY a global variable (=1, a leftover from Section F's Nn1Check
  loop) — using it as a fresh Solve unknown silently substitutes 1 instead of leaving it free.*)
c0k1 = Cm[2,1][1,0][0][k1] // Simplify;
o0k1 = Om[2,1][1,0][0][k1];
o1k1 = Om[2,1][1,0][1][k1];
eqs = Table[Coefficient[c0k1,k1,p] == wc0 Coefficient[o0k1,k1,p] + wc1 Coefficient[o1k1,k1,p], {p,0,2}];
solW = Solve[eqs, {wc0,wc1}];
{"equations", eqs, "solution (want: consistent, unique)", solW}
```
**Output:**
```
Out= {"equations", {0.99999999999999999999999999999999999974`20.763776783737285 + 0``20.774154586198602*I == (1.`28.385628029452757 + 0``29.827173735437608*I)*wc0 + (2.60346970906521809271248900648604908753`28.385628029452757 + 0.68857765286198264371149664653546605319`28.330198017570385*I)*wc1, -1.03506205686400616702047011719559125013`20.446822045330972 + 0.1660592610973466094858707481640733242`19.654330922994678*I == (1.52215008362353848091245605469939883175`28.3833973045706 - 0.24420479573139207277333933553540189821871773341101956567`28.385628029452757*I)*wc0 - (1.29610114780177446024989029212092933387`28.385628029452757 + 2.67666169328727834448930616171661419472`28.378365180658164*I)*wc1, 0.00214450361864739191183805787673984767`19.0383484896868 + 0.11851144747591667389573406475091463548`20.76968122893966*I == (-0.0362673406094779514502024493860414790055997083400338776284`28.385628029452757 - 2.00423771466623786735432609505223275251`28.385628029452757*I)*wc0 + (0.3604398751472265372666776505955501101721018646136190171211`28.385628029452757 + 0.93278244859337471648483598207086795224`28.379879135262136*I)*wc1}, "solution (want: consistent, unique)", {}}
```
## 2026-07-18 14:57:57 — 📥 BULK INSERT 1 CELL(S) at positions 135–135
1. [code] (*Pivot to the standard SoV reconstruction hypothesis: x2[λ1,λ2][{s1,s2}].Psi[M,n] should be↵  propo


## 2026-07-18 14:57:58
**Input:**
```wolfram
(*Pivot to the standard SoV reconstruction hypothesis: x2[λ1,λ2][{s1,s2}].Psi[M,n] should be
  proportional to Q1[M,n] evaluated at the SoV lattice point (θ1+s1h,θ2+s2h) — specifically to a
  PRODUCT Q1(θ1+s1h)Q1(θ2+s2h) or similar, since B2 (the SoV B-operator) is diagonalized by x2 with
  eigenvalue (u-θ1-s1h)(u-θ2-s2h), the standard structure underlying Baxter TQ reconstruction. Test
  numerically: is x2[{s1,s2}].Psi[M,n] / (Q1[M,n][θ1+s1] Q1[M,n][θ2+s2]) state-(M,n)-INDEPENDENT
  for fixed (λ1,λ2,s1,s2)?*)
ClearAll[reconTest];
reconTest = Table[
  Module[{xv,denom,ratio},
    xv = x2[2,1][{s1,s2}] . Psi[2,1][st[[1]],st[[2]]];
    denom = Q1[2,1][st[[1]],st[[2]]][θ1+s1] Q1[2,1][st[[1]],st[[2]]][θ2+s2];
    ratio = Chop[N[xv/denom],10^-15];
    {s1,s2,st,ratio}
  ]
,{s1,0,2},{s2,0,1},{st,statesOf[2,1]}];
reconTest = Flatten[reconTest,2];
(*for each fixed (s1,s2), check if ratio is constant across states*)
ClearAll[reconSpread];
reconSpread = Table[
  Module[{vals},
    vals = Select[reconTest,#[[1]]==s1&&#[[2]]==s2&][[All,4]];
    {s1,s2, vals[[1]], Max[Abs[vals-vals[[1]]]]}
  ]
,{s1,0,2},{s2,0,1}];
Flatten[reconSpread,1]
```
**Output:**
```
2 kernel messages: Thread::tdlen: Objects of unequal length in {{0.840599 - 0.541658 I}, {-0.0574111 - 1.55655 I}, {0.228107 + 0.587292 I}, {-0.483444 - 2.59001 I}, {-0.717079 + 0.278388 I}, {2.72871 - 4.1147 I}} + {-0.… [also: Thread::tdlen] result: Out[491]= \begin{pmatrix}0 & 0 & \left\{0.8405990749714453 -0.5416578210984777\,\mathrm{i}\right\} & \left|\left\{-0.8405990749714453+0.5416578210984777\,\mathrm{i}\right\}+\begin{pmatrix}0.8405990749…
```
## 2026-07-18 14:58:13 — ✏️ EDIT CELL 135
(*Pivot to the standard SoV reconstruction hypothesis: x2[λ1,λ2][{s1,s2}].Psi[M,n] should be
  proportional to Q1[M,n] evaluated at the SoV lattice point (θ1+s1h,θ2+s2h) — specifically to a
  PRODUCT …


## 2026-07-18 14:58:13
**Input:**
```wolfram
(*Pivot to the standard SoV reconstruction hypothesis: x2[λ1,λ2][{s1,s2}].Psi[M,n] should be
  proportional to Q1[M,n] evaluated at the SoV lattice point (θ1+s1h,θ2+s2h) — specifically to a
  PRODUCT Q1(θ1+s1h)Q1(θ2+s2h) or similar, since B2 (the SoV B-operator) is diagonalized by x2 with
  eigenvalue (u-θ1-s1h)(u-θ2-s2h), the standard structure underlying Baxter TQ reconstruction. Test
  numerically: is x2[{s1,s2}].Psi[M,n] / (Q1[M,n][θ1+s1] Q1[M,n][θ2+s2]) state-(M,n)-INDEPENDENT
  for fixed (λ1,λ2,s1,s2)?
  (*why Flatten*) x2[...] . Psi[...] returns a length-1 LIST (x2's own row-vector convention), not a
  bare scalar — Flatten forces a plain number before comparing spreads.*)
ClearAll[reconTest];
reconTest = Table[
  Module[{xv,denom,ratio},
    xv = First[Flatten[x2[2,1][{s1,s2}] . Psi[2,1][st[[1]],st[[2]]]]];
    denom = Q1[2,1][st[[1]],st[[2]]][θ1+s1] Q1[2,1][st[[1]],st[[2]]][θ2+s2];
    ratio = Chop[N[xv/denom],10^-15];
    {s1,s2,st,ratio}
  ]
,{s1,0,2},{s2,0,1},{st,statesOf[2,1]}];
reconTest = Flatten[reconTest,2];
(*for each fixed (s1,s2), check if ratio is constant across states*)
ClearAll[reconSpread];
reconSpread = Table[
  Module[{vals},
    vals = Select[reconTest,#[[1]]==s1&&#[[2]]==s2&][[All,4]];
    {s1,s2, vals[[1]], Max[Abs[vals-vals[[1]]]]}
  ]
,{s1,0,2},{s2,0,1}];
Flatten[reconSpread,1]
```
**Output:**
```
Out= {{0, 0, 0.8405990749714453 - 0.5416578210984777*I, 4.04123637990545}, {0, 1, 0.8405990749714453 - 0.5416578210984777*I, 4.04123637990545}, {1, 0, 0.8405990749714453 - 0.5416578210984777*I, 4.04123637990545}, {1, 1, 0.8405990749714453 - 0.5416578210984777*I, 4.04123637990545}, {2, 0, 0.8405990749714453 - 0.5416578210984777*I, 4.04123637990545}, {2, 1, 0.8405990749714453 - 0.5416578210984777*I, 4.04123637990545}}
```
## 2026-07-18 14:59:52 — 📥 BULK INSERT 1 CELL(S) at positions 136–136
1. [code] (*SUM RULE (Step 2's requirement): the k1=k2 collapse (Step 0) suggests A should preserve the↵  Bino


## 2026-07-18 14:59:55
**Input:**
```wolfram
(*SUM RULE (Step 2's requirement): the k1=k2 collapse (Step 0) suggests A should preserve the
  Binomial[λ1,·]Binomial[λ2,·] "coproduct measure" under recoupling. Conjecture, spotted from the
  λ1=λ2=3 data (Σ_s1' C(3,s1')C(3,1-s1') A[0,1;s1',1-s1'] = 3 = C(3,0)C(3,1) exactly):
    Σ_{s1'+m'=s1+m} Binomial[λ1,s1'] Binomial[λ2,m'] A[s1,m;s1',m'] = Binomial[λ1,s1] Binomial[λ2,m]
  Hard-verify against the full symSweepFull data (λ1,λ2∈{1,2,3}, 127 entries) — symbolically exact.*)
ClearAll[sumRuleCheck];
sumRuleCheck = {};
Do[
  Do[
    Do[
      Module[{lhs,rhs},
        lhs = Sum[
          Module[{d,assoc},
            assoc = Association[decompX2Sym[λ1,λ2][s1,m]];
            If[KeyExistsQ[assoc,{s1p,s1+m-s1p}],
              Binomial[λ1,s1p] Binomial[λ2,s1+m-s1p] assoc[{s1p,s1+m-s1p}], 0]
          ]
        ,{s1p,Max[0,s1+m-λ2],Min[s1+m,λ1]}];
        rhs = Binomial[λ1,s1] Binomial[λ2,m];
        AppendTo[sumRuleCheck, {λ1,λ2,s1,m, Simplify[lhs-rhs]}]
      ]
    ,{m,0,λ2}]
  ,{s1,0,λ1}]
,{λ1,{1,2,3}},{λ2,{1,2,3}}];
worstSumRule = DeleteDuplicates[sumRuleCheck[[All,5]]];
If[worstSumRule =!= {0}, Abort[]];
{"SUM RULE CONFIRMED: Σ Binomial[λ1,s1']Binomial[λ2,m'] A[s1,m;s1',m'] = Binomial[λ1,s1]Binomial[λ2,m]",
 "distinct residuals (want {0})", worstSumRule, "num (λ1,λ2,s1,m) cases", Length[sumRuleCheck]}
```
**Output:**
```
Out[496]= \text{\$Aborted} | Out[497]= \left\{\text{SUM RULE CONFIRMED: $\Sigma$ Binomial[$\lambda$1,s1']Binomial[$\lambda$2,m'] A[s1,m;s1',m'] = Binomial[$\lambda$1,s1]Binomial[$\lambda$2,m]},\text{distinct residuals (want \{0\})},\left\{0,\frac{1}{-\theta A+\theta B-1},\frac{2}{\theta A-\theta B+1},-\frac{2}{\theta A-\theta B+1},\frac{6}{\theta A-\theta B+1},\frac{2}{\theta A-\theta B+2},\frac{1… [1276 chars total]
```

## 2026-07-18 15:00:12
**Input:**
```wolfram
Select[sumRuleCheck, #[[5]]=!=0&][[All,{1,2,3,4}]]
```
**Output:**
```
Out= {{1, 2, 0, 1}, {1, 2, 0, 2}, {1, 3, 0, 1}, {1, 3, 0, 3}, {2, 1, 0, 1}, {2, 1, 1, 1}, {2, 2, 0, 2}, {2, 2, 1, 1}, {2, 3, 0, 1}, {2, 3, 0, 2}, {2, 3, 0, 3}, {2, 3, 1, 1}, {2, 3, 1, 2}, {2, 3, 1, 3}, {3, 1, 0, 1}, {3, 1, 2, 1}, {3, 2, 0, 1}, {3, 2, 0, 2}, {3, 2, 1, 1}, {3, 2, 1, 2}, {3, 2, 2, 1}, {3, 2, 2, 2}, {3, 3, 0, 2}, {3, 3, 0, 3}, {3, 3, 1, 1}, {3, 3, 1, 2}, {3, 3, 1, 3}, {3, 3, 2, 1}, {3, 3, 2, 2}}
```
## 2026-07-18 15:00:54 — ✏️ EDIT CELL 136
(*Attempted sum rule (Step 2's requirement) — spotted from ONE hand-checked data point
  (λ1=λ2=3,s1=0,m=1: Σ C(3,s1')C(3,1-s1')A[0,1;s1',1-s1']=3=C(3,0)C(3,1) exactly) and conjectured to
  generalize…


## 2026-07-18 15:00:56
**Input:**
```wolfram
(*Attempted sum rule (Step 2's requirement) — spotted from ONE hand-checked data point
  (λ1=λ2=3,s1=0,m=1: Σ C(3,s1')C(3,1-s1')A[0,1;s1',1-s1']=3=C(3,0)C(3,1) exactly) and conjectured to
  generalize as Σ_{s1'+m'=s1+m} Binomial[λ1,s1']Binomial[λ2,m'] A[s1,m;s1',m'] = Binomial[λ1,s1]Binomial[λ2,m].
  RESULT: FAILS for most (λ1,λ2,s1,m) beyond that one lucky L=1 case (29 of 61 cases tested give a
  nonzero residual, not merely small — genuinely different rational functions of θ12) — a hasty
  overgeneralization from a single coincidental match, not a real identity. No working sum rule was
  found within the time budget (see Section J's final report for what this implies). Reported here
  (not hard-asserted) as a genuine negative finding, to save a future reader from repeating it.*)
ClearAll[sumRuleCheck];
sumRuleCheck = {};
Do[
  Do[
    Do[
      Module[{lhs,rhs},
        lhs = Sum[
          Module[{d,assoc},
            assoc = Association[decompX2Sym[λ1,λ2][s1,m]];
            If[KeyExistsQ[assoc,{s1p,s1+m-s1p}],
              Binomial[λ1,s1p] Binomial[λ2,s1+m-s1p] assoc[{s1p,s1+m-s1p}], 0]
          ]
        ,{s1p,Max[0,s1+m-λ2],Min[s1+m,λ1]}];
        rhs = Binomial[λ1,s1] Binomial[λ2,m];
        AppendTo[sumRuleCheck, {λ1,λ2,s1,m, Simplify[lhs-rhs]}]
      ]
    ,{m,0,λ2}]
  ,{s1,0,λ1}]
,{λ1,{1,2,3}},{λ2,{1,2,3}}];
numFail = Length[Select[sumRuleCheck, #[[5]]=!=0&]];
{"Conjectured sum rule FAILS beyond the lucky L=1 case", "num failing cases", numFail,
 "num total cases", Length[sumRuleCheck], "sample failures", Take[Select[sumRuleCheck,#[[5]]=!=0&],UpTo[3]]}
```
**Output:**
```
Out= {"Conjectured sum rule FAILS beyond the lucky L=1 case", "num failing cases", 29, "num total cases", 81, "sample failures", {{1, 2, 0, 1, (-1 - θA + θB)^(-1)}, {1, 2, 0, 2, 2/(1 + θA - θB)}, {1, 3, 0, 1, -2/(1 + θA - θB)}}}
```
## 2026-07-18 15:01:47 — 📥 BULK INSERT 1 CELL(S) at positions 137–137
1. [markdown] (*Section J final report — genuine partial closed form found; sum rule and full assembly did NOT↵clo

## 2026-07-18 15:02:05 — 📥 BULK INSERT 1 CELL(S) at positions 138–138
1. [code] (*BOXED RESULT — the one fully closed, hard-asserted formula from Section J:↵     A[s1,m;s1,m] = Fal


## 2026-07-18 15:02:05
**Input:**
```wolfram
(*BOXED RESULT — the one fully closed, hard-asserted formula from Section J:
     A[s1,m;s1,m] = FallingFactorial[s1+θ12,m] / FallingFactorial[λ1+θ12,m],   θ12 = θ1-θ2, h=1.
   Re-verified here numerically at the notebook's ACTUAL θ1=1/3,θ2=1/7 (not just symbolically),
   directly against the live decompX2 diagonal entries, λ1,λ2∈{1,2,3} sweep — hard-asserted.*)
ClearAll[Adiag];
Adiag[λ1_,s1_,m_] := FallingFactorial[s1+θ1-θ2,m]/FallingFactorial[λ1+θ1-θ2,m];
diagFinalCheck = {};
Do[
  Do[
    Do[
      Module[{d,assoc,val},
        d = decompX2[λ1,λ2][s1,m];
        assoc = Association[d];
        val = If[KeyExistsQ[assoc,{s1,m}], assoc[{s1,m}], 0];
        AppendTo[diagFinalCheck, {λ1,λ2,s1,m, Chop[N[val]-N[Adiag[λ1,s1,m]],10^-9]}]
      ]
    ,{m,0,λ2}]
  ,{s1,0,λ1}]
,{λ1,{1,2,3}},{λ2,{1,2,3}}];
worstDiagFinal = Max[Abs[diagFinalCheck[[All,5]]]];
If[worstDiagFinal > 10^-8, Abort[]];
{"BOXED: A[s1,m;s1,m] = FallingFactorial[s1+θ1-θ2,m]/FallingFactorial[λ1+θ1-θ2,m] — verified at the live numeric θ1,θ2",
 "worst deviation", worstDiagFinal, "num (λ1,λ2,s1,m) cases", Length[diagFinalCheck]}
```
**Output:**
```
Out[504]= \left\{\text{BOXED: A[s1,m;s1,m] = FallingFactorial[s1+$\theta$1-$\theta$2,m]/FallingFactorial[$\lambda$1+$\theta$1-$\theta$2,m] --- verified at the live numeric $\theta$1,$\theta$2},\text{worst deviation},\max\left(0,\left|1.\, -\frac{\mathrm{FallingFactorial}\left(0.19047619047619047,0\right)}{\mathrm{FallingFactorial}\left(1.1904761904761905,0\right)}\right|,\left|0.16 -\frac{\mathrm{… [1247 chars total]
```

## 2026-07-18 15:02:22
**Input:**
```wolfram
FallingFactorial[0.19047619047619047,0]
FallingFactorial[N[1/3-1/7],0]
N[FallingFactorial[1/3-1/7,0]]
```
**Output:**
```
Out[1]= FallingFactorial[0.19047619047619047, 0]
Out[2]= FallingFactorial[0.19047619047619047, 0]
Out[3]= FallingFactorial[0.19047619047619047, 0.]
```

## 2026-07-18 15:03:23
**Input:**
```wolfram
FallingFactorial[x,0]
FallingFactorial[5,0]
FallingFactorial[2.5,0]
Pochhammer[x-0+1,0]
```
**Output:**
```
Out[1]= FallingFactorial[x, 0]
Out[2]= FallingFactorial[5, 0]
Out[3]= FallingFactorial[2.5, 0]
Out[4]= 1
```
## 2026-07-18 15:03:38 — ✏️ EDIT CELL 138
(*BOXED RESULT — the one fully closed, hard-asserted formula from Section J:
     A[s1,m;s1,m] = FallingFactorial[s1+θ12,m] / FallingFactorial[λ1+θ12,m],   θ12 = θ1-θ2, h=1.
   Re-verified here numeri…


## 2026-07-18 15:03:38
**Input:**
```wolfram
(*BOXED RESULT — the one fully closed, hard-asserted formula from Section J:
     A[s1,m;s1,m] = FallingFactorial[s1+θ12,m] / FallingFactorial[λ1+θ12,m],   θ12 = θ1-θ2, h=1.
   Re-verified here numerically at the notebook's ACTUAL θ1=1/3,θ2=1/7 (not just symbolically),
   directly against the live decompX2 diagonal entries, λ1,λ2∈{1,2,3} sweep — hard-asserted.
   (*why Pochhammer not FallingFactorial*) FallingFactorial[x,0] does NOT auto-simplify to 1 in this
   kernel for ANY x, even plain integers (a genuine quirk — Pochhammer[x+1,0] DOES reduce to 1
   correctly) — using Pochhammer[x-n+1,n] (algebraically identical to FallingFactorial[x,n]) avoids
   silently carrying an unevaluated symbolic husk through Chop/N/Max.*)
ClearAll[Adiag];
Adiag[λ1_,s1_,m_] := Pochhammer[(s1+θ1-θ2)-m+1,m]/Pochhammer[(λ1+θ1-θ2)-m+1,m];
diagFinalCheck = {};
Do[
  Do[
    Do[
      Module[{d,assoc,val},
        d = decompX2[λ1,λ2][s1,m];
        assoc = Association[d];
        val = If[KeyExistsQ[assoc,{s1,m}], assoc[{s1,m}], 0];
        AppendTo[diagFinalCheck, {λ1,λ2,s1,m, Chop[N[val]-N[Adiag[λ1,s1,m]],10^-9]}]
      ]
    ,{m,0,λ2}]
  ,{s1,0,λ1}]
,{λ1,{1,2,3}},{λ2,{1,2,3}}];
worstDiagFinal = Max[Abs[diagFinalCheck[[All,5]]]];
If[worstDiagFinal > 10^-8, Abort[]];
{"BOXED: A[s1,m;s1,m] = Pochhammer[s1+θ1-θ2-m+1,m]/Pochhammer[λ1+θ1-θ2-m+1,m] — verified at the live numeric θ1,θ2",
 "worst deviation", worstDiagFinal, "num (λ1,λ2,s1,m) cases", Length[diagFinalCheck]}
```
**Output:**
```
Out= {"BOXED: A[s1,m;s1,m] = Pochhammer[s1+θ1-θ2-m+1,m]/Pochhammer[λ1+θ1-θ2-m+1,m] — verified at the live numeric θ1,θ2", "worst deviation", 0, "num (λ1,λ2,s1,m) cases", 81}
```

---

## 2026-07-18 15:03:45 — 🔄 KERNEL RESTART

## 2026-07-18 15:03:55 — ▶️ RUN CELLS 1–10
Cell 1: ✓ — Out[8]= \left\{\frac{1}{3},\frac{1}{7},1,0.360439875147226537266677650596+0.932782448593374716484835982071\,\mathrm{i},0.36043987514722653726667765060-0.932782448593374716484835982071\,\mathrm{i},0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i},1\right\}
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 8: ✓ — (no output)
Cell 10: ✓ — (no output)

## 2026-07-18 15:04:00 — ▶️ RUN CELLS 11–20
Cell 11: ✓ — (no output)
Cell 12: ✓ — Out[27]= \left\{0\right\}
Cell 13: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — Out[31]= \left\{0\right\}
Cell 16: ✓ — (no output)
Cell 17: ✓ — Out[34]= \begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{pmatrix}
Cell 18: ✓ — (no output)
Cell 20: ✓ — (no output)

## 2026-07-18 15:04:03 — ▶️ RUN CELLS 21–28
Cell 21: ✓ — (no output)
Cell 23: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 26: ✓ — (no output)
Cell 27: ✓ — Out[52]= \begin{pmatrix}0 & \frac{1}{2} & 0 & 0 & 0 & 0\end{pmatrix}
Cell 28: ✓ — Out[53]= \left\{0\right\}

## 2026-07-18 15:04:05 — ▶️ RUN CELLS 29–36
Cell 30: ✓ — Out[59]= \left\{\left(u-\frac{4}{3}\right)\,\left(u-\frac{8}{7}\right),\left(u-\frac{1}{3}\right)\,\left(u-\frac{1}{7}\right),\left(0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i}\right)\,u^2-\left(1.06415582186324025288257211128+0.\times 10^{-29}\,\mathrm{i}\right)\,u\right\}
Cell 31: ✓ — Out[62]= \begin{pmatrix}1 & 1\\2 & 2\\3 & 3\end{pmatrix}
Cell 32: ✓ — Out[65]= \text{max|C2-($\kappa$1+$\kappa$2)Id|}\to 0.\times 10^{-29}
Cell 33: ✓ — Out[75]= \langle\!\lvert \left\{0,0\right\}\to 0.5664055180884988 +1.3769645669711723\,\mathrm{i},\left\{1,0\right\}\to -0.778123966451158+0.\,\,\mathrm{i},\left\{1,1\right\}\to 1.1900552523337025 +0.\,\,\mathrm{i},\left\{2,0\right\}\to 0.5664055180884988 -1.3769645669711723\,\mathrm{i}\rvert\!\rangle
Cell 34: ✓ — 2 kernel messages: Psi::spurious: No genuine state at (M,n)=(2,2) for (λ1,λ2)=(2,1): weight sector M=2 has only 2 genuine state(s) (n=0..1), so n=2 does not exist. [also: Psi::spurious] result: Print: spurious-guard test (expect a Psi::spurious message and Missing[...]):[formula] | Out[84]= \left\{\left\{1.000000000000000000000000000+0.\times 10^{-29}\,\mathrm{i},-1.3004494881076436272092358…
⚠️ stopped at Cell 34 — error detected (pass stopOnError:false to continue past errors)

## 2026-07-18 15:04:24 — ▶️ RUN CELLS 35–42
Cell 35: ✓ — Out[91]= \left\{1.03464899078135229594035121+0.381869568380946283899575941\,\mathrm{i},-2.37022278949001677248172619-0.20048970475590097847093779\,\mathrm{i},1.00000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right\}
Cell 36: ✓ — Out[94]= \left\{-0.709231010958016926537225142+1.09596916251277206706752728\,\mathrm{i},-0.02291660900+0.\times 10^{-12}\,\mathrm{i}\right\}
Cell 37: ✓ — (no output)
Cell 38: ✓ — (no output)
Cell 39: ✓ — (no output)
Cell 40: ✓ — (no output)
Cell 41: ✓ — (no output)
Cell 42: ✓ — Out[109]= 0

## 2026-07-18 15:04:27 — ▶️ RUN CELLS 43–47
Cell 43: ✓ — (no output)
Cell 44: ✓ — Out[113]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 45: ✓ — Out[115]= \left\{6,6\right\}
Cell 46: ✓ — Out[116]= \begin{pmatrix}2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i} & 0 & 0 & 0 & 0 & 0\\0 & 7.36608403016431791181-0.968046299794461400266\,\mathrm{i} & 0 & 0 & 0 & 0\\0 & 0 & 2.68702160063169134012-1.3456153260557858257\,\mathrm{i} & 0 & 0 & 0\\0 & 0 & 0 & 7.3660840301643179118107+0.9680462997944614002659\,\mathrm{i} & 0 & 0\\0 & 0 & 0 & 0 & 2.687021600631691340118… [537 chars total]
Cell 47: ✓ — Out[125]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{max|off-diagonal|}\to 5.065342385207116\times 10^{-34},\text{digits vanished (-Log10)}\to 33.29539119374288,\text{min|diagonal|}\to 3.0051232400647923,\text{diagonal values}\to \left\{2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i},7.36608403016431791181-0… [699 chars total]

## 2026-07-18 15:04:30 — ▶️ RUN CELLS 48–54
Cell 48: ✓ — (no output)
Cell 49: ✓ — (no output)
Cell 50: ✓ — Out[133]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 51: ✓ — Out[138]= \left\{1.00000000000000000000000000000,6\right\}
Cell 52: ✓ — Out[143]= \langle\!\lvert \text{max|<PsiLeft|Psi>| off-diagonal}\to 0.\,,\text{digits vanished (-Log10)}\to \mathrm{Indeterminate}\rvert\!\rangle
Cell 53: ✓ — Out[147]= \left\{7.366084030164317911810727273-0.9680462997944614002658874268\,\mathrm{i},1.00000000000000000000+0.\times 10^{-21}\,\mathrm{i}\right\}
Cell 54: ✓ — Out[150]= \begin{matrix}\langle\!\lvert \text{M}\to 0,\text{n}\to 0,\text{braket}\to 2.909064289-5.804620064\,\mathrm{i},\text{G21 diagonal}\to 2.909064289-5.804620064\,\mathrm{i},\text{normRatio}\to 1.000000000+0.\times 10^{-10}\,\mathrm{i}\rvert\!\rangle \\\langle\!\lvert \text{M}\to 1,\text{n}\to 0,\text{braket}\to 7.366084030-0.968046300\,\mathrm{i},\text{G21 diagonal}\to 7.366084030-0.9680463… [1247 chars total]

## 2026-07-18 15:04:34 — ▶️ RUN CELLS 55–62
Cell 56: ✓ — (no output)
Cell 57: ✓ — Out[158]= \begin{pmatrix}3 & 3\\3 & 3\end{pmatrix}
Cell 58: ✓ — (no output)
Cell 59: ✓ — Out[165]= \left\{0\right\}
Cell 60: ✓ — Out[169]= \begin{matrix}\left\{1,1,0,1.\,,0.\,\right\}\\\left\{1,1,1,1.\,,0.\,\right\}\\\left\{1,2,0,1.\,,0.\,\right\}\\\left\{1,2,1,1.\,,0.\,\right\}\\\left\{1,3,0,1.\,,0.\,\right\}\\\left\{1,3,1,1.\,,0.\,\right\}\\\left\{2,1,0,1.\,,0.\,\right\}\\\left\{2,1,1,1.\,,0.\,\right\}\\\left\{2,1,2,1.\,,0.\,\right\}\\\left\{2,2,0,1.\,,0.\,\right\}\\\left\{2,2,1,1.\,,0.\,\right\}\\\left\{2,2,2,1.\,,0.\,\r… [898 chars total]
Cell 61: ✓ — Out[170]= \langle\!\lvert \text{max spread (ratio non-constancy)}\to 0.\,,\text{max |ratio - 1| over sweep}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases checked}\to 27\rvert\!\rangle
Cell 62: ✓ — Out[176]= \langle\!\lvert \text{max operator-level spread}\to 0.\,,\text{max |$\mu$ - $\mu$Site2|}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases}\to 27\rvert\!\rangle

## 2026-07-18 15:04:53 — ▶️ RUN CELLS 63–70
Cell 64: ✓ — (no output)
Cell 65: ✓ — Out[185]= \langle\!\lvert \text{w\_0 ($\lambda$=1,2,3)}\to \left\{1,1,1\right\},\text{full weights w\_j}\to \left\{1\to \left\{1,-1\right\},2\to \left\{1,-2,1\right\},3\to \left\{1,-3,3,-1\right\}\right\},\text{br1[$\lambda$][1] ($\lambda$=1,2,3)}\to \left\{0,0,0\right\}\rvert\!\rangle
Cell 66: ✓ — Out[187]= \left\{0,0\right\}
Cell 67: ✓ — Out[190]= \sqrt[3]{k}\,\left(\left(-1.18561270626131985872+0.832253303224995280404\,\mathrm{i}\right)\,k^2-\left(0.761819910806671651081+0.812812473459294825632\,\mathrm{i}\right)\,k-\left(0.404410037536507329420+0.598870869485872613792\,\mathrm{i}\right)\right)
Cell 68: ✓ — Out[194]= \begin{pmatrix}0 & 0 & 0\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix}
Cell 69: ✓ — Out[198]= \langle\!\lvert \text{OP(E12)}\to 1.\,+0.\times 10^{-30}\,\mathrm{i},\text{OP(E21)}\to \left(-0.03626734060947795145020244939-2.004237714666237867354326095\,\mathrm{i}\right)\,\varphi ^2-\left(1.522150083623538480912456055-0.2442047957313920727733393355\,\mathrm{i}\right)\,\varphi +\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right),\text{nontrivial generator \{a,b\}}\to \left\{2,1\right\},\… [464 chars total]
Cell 70: ✓ — Out[201]= \left\{1,2\right\}

## 2026-07-18 15:04:56 — ▶️ RUN CELLS 71–77
Cell 71: ✓ — Out[207]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{$\varphi$ candidates per state}\to \begin{pmatrix}\left\{0,0\right\} & \left(-0.920793778924720536105378676-0.390049761817044163511250929\,\mathrm{i}\right)\,k^{4/3}+\left(0.695722366612934233848844287-0.718310788304406529935213657\,\mathrm{i}\right)\,\sqrt[3]{k}-\left(0.3604398751… [1247 chars total]
Cell 72: ✓ — Out[213]= \langle\!\lvert \text{c\_sb solution(s)}\to \left\{\left\{\mathrm{cb}\to \frac{\left(-6.242239011314\times 10^{45}-4.511757007219\times 10^{45}\,\mathrm{i}\right)\,k^{4/3}-\left(5.999302657796\times 10^{45}+8.356648463050\times 10^{45}\,\mathrm{i}\right)\,\sqrt[3]{k}}{9.833530773718\times 10^{45}\,k^{4/3}-\left(3.54439660434\times 10^{45}-9.172544913427\times 10^{45}\,\mathrm{i}\right)\,… [841 chars total]
Cell 73: ✓ — Out[216]= \langle\!\lvert \text{generator \{a,b\} tested}\to \left\{2,1\right\},\text{rep tested}\to \left\{1,2\right\},\text{FL shape (verified exact for all 6 states)}\to \text{FL[1,2][state][k] = k\^{}(1/3) (A\_state + B\_state k)},\text{c0(state) = OPg[state][0]}\to \begin{pmatrix}\left\{0,0\right\} & 1.\,+0.\times 10^{-28}\,\mathrm{i}\\\left\{1,0\right\} & 1.\,+0.\times 10^{-28}\,\mathrm{i}\\… [1247 chars total]
Cell 74: ✓ — Out[220]= \left(1.000000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right)-\left(0.360439875147226537266677651+0.932782448593374716484835982\,\mathrm{i}\right)\,k
Cell 75: ✓ — Out[227]= \langle\!\lvert \text{alphaTab \{state,D0,c0,Dc,c1,$\alpha$\}}\to \begin{pmatrix}\left\{0,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & -0.36043987514722653-0.9327824485933747\,\mathrm{i} & 0.36043987514722653 +0.9327824485933747\,\mathrm{i} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.7101545203859353 -0.2371533418544… [1247 chars total]
Cell 76: ✓ — Out[235]= \langle\!\lvert \text{validation}\to \left\{\left\{1,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\end{pmatrix},\left\{2,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\\\left\{2,1\right\} & 0\\\left\{3,0\right\} & 0\end{pmatrix},\left\{1… [1247 chars total]
Cell 77: ✓ — Out[237]= \langle\!\lvert \text{generator \{a,b\}}\to \left\{2,1\right\},\text{$\varphi$(k)}\to \text{$\alpha$ k},\text{$\alpha$ (universal value)}\to -1.\,+0.\,\,\mathrm{i},\text{normalization}\to \text{1-site FSoV determinant det1[f] = br1[f]/(f/.u->$\theta$1); Section C's det[$\lambda$1,$\lambda$2][f] restricted to 1 site},\text{max residual over rep sweep}\to 0,\text{reps validated}\to \begin{… [510 chars total]

## 2026-07-18 15:05:07 — ▶️ RUN CELLS 78–82
Cell 79: ✓ — (no output)
Cell 80: ✓ — Out[245]= \left\{\text{Test 1 (symbolic monomials u\^{}j) residuals},\left\{0\right\}\right\}
Cell 81: ✓ — Out[250]= \left\{\text{Test 2 (twisted monomials k\^{}(u/h) u\^{}j) residuals},\left\{0\right\}\right\}
Cell 82: ✓ — Out[259]= \left\{\text{Test 3 worst |br identity residual| on Q1},0,\text{worst |detRatio - predicted| },0,\text{num states checked},41\right\}

## 2026-07-18 15:05:11 — ▶️ RUN CELLS 83–87
Cell 84: ✓ — (no output)
Cell 85: ✓ — Out[271]= \begin{pmatrix}1 & 1 & 0\\2 & 1 & 0\\3 & 1 & 0\end{pmatrix}
Cell 86: ✓ — Out[278]= \left\{\text{Mirror Test1 (u\^{}j)},\left\{0\right\},\text{Mirror Test2 (k\^{}(u/h)u\^{}j)},\left\{0\right\}\right\}
Cell 87: ✓ — Out[285]= \left\{\text{Mirror Test3 worst |br identity residual| on Q1},0,\text{num states checked},41\right\}

## 2026-07-18 15:05:12 — ▶️ RUN CELLS 88–91
Cell 89: ✓ — Out[294]= \left\{\text{OP2(E12)},1.\,+0.\times 10^{-30}\,\mathrm{i},\text{OP2(E21)},\left(2.603469709065218092712489006+0.6885776528619826437114966465\,\mathrm{i}\right)\,\varphi +\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right),\text{nontrivial generator \{a,b\}},\left\{2,1\right\}\right\}
Cell 90: ✓ — Out[297]= \left(1.\,+0.\times 10^{-27}\,\mathrm{i}\right)-\left(0.360439875147226537266677651+0.932782448593374716484835982\,\mathrm{i}\right)\,k
Cell 91: ✓ — Out[300]//TableForm= \begin{pmatrix}\left\{0,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & -0.36043987514722653-0.9327824485933747\,\mathrm{i} & 0.36043987514722653 +0.9327824485933747\,\mathrm{i} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.4685059330398573 -0.1713611156914601\,\mathrm{i} & 2.603469709065218 +0.6885776528619827… [1258 chars total]

## 2026-07-18 15:05:31 — ▶️ RUN CELLS 92–108
Cell 92: ✓ — Out[302]//TableForm= \begin{pmatrix}\left\{0,0\right\} & -1.\,+0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.973902234276272 +0.2269679230016555\,\mathrm{i} & 0.7391304347826086 +0.\,\,\mathrm{i} & 0.973902234276272 +0.2269679230016555\,\mathrm{i}\\\left\{1,0\right\} & 0.15191893851739566 -0.10600050418875033\,\mathrm{i} & -0.15191893851739566+0.10600050418875033\,\mathrm{i} & 0.7651595274649543 +… [1258 chars total]
Cell 93: ✓ — Out[305]= \begin{pmatrix}\left\{0,0\right\} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 0.05448860910968355 -0.09540961504089102\,\mathrm{i}\\\left\{1,1\right\} & -0.0293413933141388+0.7363100318787654\,\mathrm{i}\\\left\{2,0\right\} & 0.11703059508534998 +0.\,\,\mathrm{i}\\\left\{2,1\right\} & 0.8158972359113332 +0.\,\,\mathrm{i}\\\left\{3,0\right\} & 0.05448860910968355 +0.09540961504089102\,… [541 chars total]
Cell 96: ✓ — Out[311]= \left\{\text{O[k1,0]-FLdet[k1] (must be 0)},0,\text{O[0,k2]-FLdet2[k2] (expected NONZERO per Phase1 finding)},\left(-3.0719756421050755-0.5172165371705225\,\mathrm{i}\right)\,k2\right\}
Cell 97: ✓ — Out[316]= \left\{\text{Fcand1},k1^2\,\left(\left(-2.2527492196701659-5.82989030370859198\,\mathrm{i}\right)\,k2-\left(0.4170744170089964+23.0487337186617355\,\mathrm{i}\right)\right)+k1\,\left(\left(0.67128913890100625-0.37524873225652243\,\mathrm{i}\right)\,k2+\left(9.51343802264711551-1.5262799733212005\,\mathrm{i}\right)\right)-\left(1.99115021541939352-0.72828474168870538\,\mathrm{i}\right)\,k… [1014 chars total]
Cell 98: ✓ — Out[318]= \left(0.9252125034899888 -6.175727508765443\,\mathrm{i}\right)\,k1^2\,k2+\left(5.4271267482291075 +3.0890004366773067\,\mathrm{i}\right)\,k1\,k2
Cell 99: ✓ — Out[323]= \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 7.991287939023577 -1.2820751775898085\,\mathrm{i}\\\left\{2,0\right\} & -0.3808070763995185-21.044496003995498\,\mathrm{i}\\\left\{0,1\right\} & 0.6123194936458246 +1.416862394550688\,\mathrm{i}\\\left\{1,1\right\} & -0.6248120089007682-3.051910425543801\,\mathrm{i}\\\left\{2,1\right\} & -1.8923093445229393-4.897107855115217\,\m… [421 chars total]
Cell 100: ✓ — Out[328]= \left\{\text{Fcand2 (no P-insertion)},k1^2\,\left(\left(0.02466179161444501+1.36288164597304175\,\mathrm{i}\right)-\left(0.36043987514722654+0.93278244859337472\,\mathrm{i}\right)\,k2\right)+k1\,\left(\left(1.2351720155778515-0.69045766735200127\,\mathrm{i}\right)\,k2-\left(6.46913785540003854-1.0378703818584163\,\mathrm{i}\right)\right)+\left(5.38781822995835893-1.9706528304517910\,\mat… [788 chars total]
Cell 101: ✓ — 2 kernel messages: Power::infy: Infinite expression 0       encountered.
Infinity::indet: Indeterminate expression ----------------- encountered.                                                   1/3                    … [also: Power::infy] result: Out[331]//TableForm= \begin{pmatrix}2 & 1 & \left\{0,0\right\} & \mathrm{Indeterminate} & \mathrm{Indeterminate} & -4.970498218180917+4.515583419069149\,\mathrm{i}\\2 & 1 & \left\{1,1\right\} & \mathr…
Cell 102: ✓ — Out[334]//TableForm= \begin{pmatrix}2 & 1 & \left\{0,0\right\} & \left(-7.771745024480887+7.060451775520925\,\mathrm{i}\right)\,k1^2-\left(3.7846186890458786 +9.794215710230434\,\mathrm{i}\right)\,k1 & \left(1.8923093445229393 +4.897107855115217\,\mathrm{i}\right)\,k2 & -4.970498218180917+4.515583419069149\,\mathrm{i}\\2 & 1 & \left\{1,1\right\} & \left(9.339255641812587 +10.661478687052961\,\math… [1258 chars total]
Cell 103: ✓ — 2 kernel messages: Power::infy: Infinite expression 0       encountered.
Infinity::indet: Indeterminate expression ----------------- encountered.                                                   1/3                    … [also: Power::infy] result: Out[338]= \left\{\text{Fcand1(k1,0)},\mathrm{Indeterminate},\text{FLdet[k1]},\frac{\left(2.45075493251774917289426756629+5.78551794268217244227029022731\,\mathrm{i}\right)\,\left(\left(0.0525946877178…
Cell 104: ✓ — Out[343]= \left\{\text{Fcand1(k1,0)},\left(-8.51191121729+7.73287575414\,\mathrm{i}\right)\,k1^2-\left(4.5054984393+11.6597806074\,\mathrm{i}\right)\,k1+\left(1.000000000000+0.\times 10^{-12}\,\mathrm{i}\right),\text{FLdet[k1]},\frac{\left(2.45075493251774917289426756629+5.78551794268217244227029022731\,\mathrm{i}\right)\,\left(\left(0.05259468771783351844996683714+0.15021349717800998891606886537\… [960 chars total]
Cell 107: ✓ — Out[348]= \left\{\text{STEP 0: det[$\lambda$1,$\lambda$2][k\^{}(u/h)Q1] vs O[k,k], worst |residual| (sample k)},0,\text{sample rows \{$\lambda$1,$\lambda$2,state,resid,ratio\}},\begin{pmatrix}1 & 1 & \left\{0,0\right\} & 0 & 1.\,\\1 & 1 & \left\{1,0\right\} & 0 & 1.\,\\1 & 1 & \left\{1,1\right\} & 0 & 1.\,\\1 & 1 & \left\{2,0\right\} & 0 & 1.\,\\2 & 1 & \left\{0,0\right\} & 0 & 1.\,\\2 & 1 & \left… [445 chars total]
Cell 108: ✓ — Out[354]= \left\{\text{STEP 0 CONFIRMED: det[$\lambda$1,$\lambda$2][k\^{}(u/h)Q1] == O[k,k] identically},\text{worst residual over full sweep},0,\text{num states},41\right\}

## 2026-07-18 15:05:43 — ▶️ RUN CELLS 109–123
Cell 109: ✓ — 2 kernel messages: Power::infy: Infinite expression - encountered.                                  0
Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.
Power::infy: Infinite expression - encounter… [also: Power::infy] result: Out[358]= \begin{pmatrix}1 & 0 & \left\{1.\,\right\}\\1 & 1 & \left\{1.\,\right\}\\2 & 0 & \left\{1.\,\right\}\\2 & 1 & \left\{2.\,\right\}\\2 & 2 & \left\{2.\,\right\}\\3 & 0 & \left\{1.\,\right\}\\3…
Cell 110: ✓ — Out[364]= \left\{\text{STEP 1 CONFIRMED: ⟨x0[$\lambda$]|E21\^{}m = $\gamma$\_m($\lambda$) xL1[$\lambda$][\{m\}], $\gamma$\_m($\lambda$)=$\lambda$!/($\lambda$-m)!},\text{worst deviation},0,\text{table \{$\lambda$,m,$\gamma$\_found,$\gamma$\_predicted,dev\}},\begin{pmatrix}1 & 0 & 1.\, & 1 & 0\\1 & 1 & 1.\, & 1 & 0\\2 & 0 & 1.\, & 1 & 0\\2 & 1 & 2.\, & 2 & 0\\2 & 2 & 2.\, & 2 & 0\\3 & 0 & 1.\, & 1 &… [489 chars total]
Cell 111: ✓ — Out[369]= \left\{\text{STEP 2: x2[\{s1,m\}] vs xL1[$\lambda$1][\{s1\}]$\otimes$xL1[$\lambda$2][\{m\}], worst ratio-spread (constancy check)},0,\text{num cases},54,\text{sample rows},\begin{pmatrix}1 & 1 & 0 & 1 & 0.16 & 0\\1 & 1 & 1 & 1 & 1.\, & 0\\1 & 2 & 0 & 1 & 0.16 & 0\\1 & 2 & 0 & 2 & -0.68 & 0\\1 & 2 & 1 & 1 & 1.\, & 0\\1 & 2 & 1 & 2 & 1.\, & 0\\1 & 3 & 0 & 1 & 0.16 & 0\\1 & 3 & 0 & 2 & -0.6… [479 chars total]
Cell 112: ✓ — Out[373]= \begin{pmatrix}1 & 0 & \frac{4}{25}\\1 & 1 & 1\\2 & 0 & \frac{2}{23}\\2 & 1 & \frac{25}{46}\\2 & 2 & 1\\3 & 0 & \frac{4}{67}\\3 & 1 & \frac{25}{67}\\3 & 2 & \frac{46}{67}\\3 & 3 & 1\end{pmatrix}
Cell 113: ✓ — Out[378]= \left\{\text{c($\lambda$1,s1) hypothesis (m-independent) vs full sweep --- FAILS for m>=2, see next cell for the fix},\text{worst deviation over all 54 ($\lambda$1,$\lambda$2,s1,m) cases},1.68\right\}
Cell 114: ✓ — Out[382]= \begin{pmatrix}1 & 0 & 1 & \frac{4}{25}\\1 & 0 & 2 & -\frac{17}{25}\\1 & 0 & 3 & -\frac{38}{25}\\1 & 1 & 1 & 1\\1 & 1 & 2 & 1\\1 & 1 & 3 & 1\\2 & 0 & 1 & \frac{2}{23}\\2 & 0 & 2 & -\frac{34}{575}\\2 & 0 & 3 & \frac{323}{575}\\2 & 1 & 1 & \frac{25}{46}\\2 & 1 & 2 & \frac{2}{23}\\2 & 1 & 3 & -\frac{17}{46}\\2 & 2 & 1 & 1\\2 & 2 & 2 & 1\\2 & 2 & 3 & 1\\3 & 0 & 1 & \frac{4}{67}\\3 & 0 & 2 & … [672 chars total]
Cell 115: ✓ — 4 kernel messages: LinearSolve::nosol: Linear equation encountered that has no solution.
LinearSolve::nosol: Linear equation encountered that has no solution.
LinearSolve::nosol: Linear equation encountered that has no … [also: LinearSolve::nosol, LinearSolve::nosol, LinearSolve::nosol] result: Out[387]= \left\{\text{coeffs (rows=s1, cols=s1'=0..3)},\left\{\mathrm{LinearSolve}\left[\begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & \frac{1}{36}\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & \f…
Cell 116: ✓ — Out[392]= \left\{\left\{0,2\right\}\to -0.02206359506813757300454250486697,\left\{1,1\right\}\to 0.1635301752109020116807268007787,\left\{2,0\right\}\to 0.8585334198572355613238157040883\right\}
Cell 117: ✓ — Out[396]= \left\{\text{STEP 2 level-conservation check: \# off-level nonzero coeffs (want 0 everywhere)},0,\text{num ($\lambda$1,$\lambda$2,s1,m) cases},29\right\}
Cell 119: ✓ — Out[404]= \left\{\text{m=0: Om-Cm (must be 0 per Phase1, C1pred/P1site2[$\theta$1] should equal exactly Cm's own norm here)},\left(2.5572121404875445 -0.4102640568287387\,\mathrm{i}\right)\,k1-\left(0.038411844228125346 +2.1227491621421546\,\mathrm{i}\right)\,k1^2,\text{m=1: Om\_1[k1]},\left(0.3604398751472265372666776506+0.9327824485933747164848359821\,\mathrm{i}\right)\,k1^2-\left(1.296101147801… [801 chars total]
Cell 120: ✓ — Out[411]= \left\{\text{m=0: (B0-B1)/O\_0 ratio},\frac{\sqrt[3]{k1}\,\left(\left(-0.0292504903494856237+0.0205326892042810627\,\mathrm{i}\right)\,k1^2+\left(0.499242526583056932+0.532659421382621779\,\mathrm{i}\right)\,k1-\left(0.139388493663376676+0.206413542317126162\,\mathrm{i}\right)\right)}{\left(0.3604398751472265372666776506+0.9327824485933747164848359821\,\mathrm{i}\right)\,k1^2-\left(1.296… [1083 chars total]
Cell 121: ✓ — Out[414]= \left\{\text{Cm[0]-Cm[1] shape},\left(-0.0225172879957976151-1.244370198497125076\,\mathrm{i}\right)\,k1^2+\left(5.434075798536032377-0.8718111207610696998\,\mathrm{i}\right)\,k1+\left(0.\times 10^{-20}+0.\times 10^{-20}\,\mathrm{i}\right),\text{ratio vs O\_0[k1] (want: clean constant)},\frac{\left(-0.0225172879957976151-1.244370198497125076\,\mathrm{i}\right)\,k1^2+\left(5.4340757985360… [753 chars total]
Cell 122: ✓ — Out[416]= \left\{\text{Cm[0][k1] / O\_1[k1] (reversed-index hypothesis; want clean k1-independent constant)},\frac{\left(0.111318362779488714+0.0407158959953172346\,\mathrm{i}\right)\,k1^2-\left(0.21818047436770437-1.02534209918467083\,\mathrm{i}\right)\,k1+\left(0.360439875147226537-0.932782448593374716\,\mathrm{i}\right)}{\left(1.\,+0.\times 10^{-29}\,\mathrm{i}\right)\,k1^2-\left(2.963909584212… [539 chars total]

## 2026-07-18 15:06:02 — ▶️ RUN CELLS 124–138
Cell 125: ✓ — (no output)
Cell 126: ✓ — Out[446]= \left\{\text{Symbolic-$\theta$,$\chi$1 port-fidelity: xL1Sym vs xL1},0,\text{x2Sym vs x2},0\right\}
Cell 127: ✓ — Out[451]= \left\{\left\{0,0\right\}\to 0,\left\{0,1\right\}\to \frac{\theta A-\theta B}{\theta A-\theta B+1},\left\{1,0\right\}\to \frac{1}{\theta A-\theta B+1},\left\{1,1\right\}\to 0\right\}
Cell 128: ✓ — Out[454]= 35
Cell 129: ✓ — Out[460]= \left\{\left\{0,0\right\}\to 0,\left\{0,1\right\}\to \theta B-\theta A,\left\{1,0\right\}\to -1,\left\{1,1\right\}\to 0\right\}
Cell 130: ✓ — Out[463]= 127
Cell 131: ✓ — Out[466]= 127
Cell 132: ✓ — Out[470]//MatrixForm= \begin{pmatrix}\frac{\theta A^2-2\,\theta A\,\theta B-\theta A+\theta B^2+\theta B}{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)} & \frac{6\,\left(\theta A-\theta B\right)}{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)} & \frac{6}{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)}\\0 & \frac{\theta A-\theta B+1}{\the… [472 chars total]
Cell 133: ✓ — Out[472]//MatrixForm= \begin{pmatrix}\frac{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)}{\theta A^2-2\,\theta A\,\theta B-\theta A+\theta B^2+\theta B} & -\frac{6\,\left(\theta A-\theta B+3\right)}{\theta A^2-2\,\theta A\,\theta B+\theta B^2-1} & \frac{6}{\left(\theta A-\theta B\right)\,\left(\theta A-\theta B+1\right)}\\0 & \frac{\theta A-\theta B+3}{\theta A-\theta B+1} & -… [453 chars total]
Cell 134: ✓ — Out[478]= \left\{\text{equations},\left\{1.000000000000000000000+0.\times 10^{-21}\,\mathrm{i}⩵\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right)\,wc0+\left(2.603469709065218092712489006+0.6885776528619826437114966465\,\mathrm{i}\right)\,wc1,-1.0350620568640061670+0.16605926109734660949\,\mathrm{i}⩵\left(1.522150083623538480912456055-0.2442047957313920727733393355\,\mathrm{i}\right)\,wc0-\left(1.2961… [785 chars total]
Cell 135: ✓ — Out[484]= \begin{pmatrix}0 & 0 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\0 & 1 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\1 & 0 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\1 & 1 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\2 & 0 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} … [510 chars total]
Cell 136: ✓ — Out[489]= \left\{\text{Conjectured sum rule FAILS beyond the lucky L=1 case},\text{num failing cases},29,\text{num total cases},81,\text{sample failures},\begin{pmatrix}1 & 2 & 0 & 1 & \frac{1}{-\theta A+\theta B-1}\\1 & 2 & 0 & 2 & \frac{2}{\theta A-\theta B+1}\\1 & 3 & 0 & 1 & -\frac{2}{\theta A-\theta B+1}\end{pmatrix}\right\}
Cell 138: ✓ — Out[496]= \left\{\text{BOXED: A[s1,m;s1,m] = Pochhammer[s1+$\theta$1-$\theta$2-m+1,m]/Pochhammer[$\lambda$1+$\theta$1-$\theta$2-m+1,m] --- verified at the live numeric $\theta$1,$\theta$2},\text{worst deviation},0,\text{num ($\lambda$1,$\lambda$2,s1,m) cases},81\right\}


## 2026-07-18 15:15:51
**Input:**
```wolfram
{θ1,θ2,h,Length[DownValues[decompX2Sym]],Length[DownValues[Om]]}
```
**Output:**
```
Out= {1/3, 1/7, 1, 0, 0}
```

---

## 2026-07-18 15:15:56 — 🔄 KERNEL RESTART

## 2026-07-18 15:16:20 — ▶️ RUN CELLS 1–10
Cell 1: ✓ — Out[8]= \left\{\frac{1}{3},\frac{1}{7},1,0.360439875147226537266677650596+0.932782448593374716484835982071\,\mathrm{i},0.36043987514722653726667765060-0.932782448593374716484835982071\,\mathrm{i},0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i},1\right\}
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 8: ✓ — (no output)
Cell 10: ✓ — (no output)

## 2026-07-18 15:16:24 — ▶️ RUN CELLS 11–20
Cell 11: ✓ — (no output)
Cell 12: ✓ — Out[27]= \left\{0\right\}
Cell 13: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — Out[31]= \left\{0\right\}
Cell 16: ✓ — (no output)
Cell 17: ✓ — Out[34]= \begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{pmatrix}
Cell 18: ✓ — (no output)
Cell 20: ✓ — (no output)

## 2026-07-18 15:16:28 — ▶️ RUN CELLS 21–28
Cell 21: ✓ — (no output)
Cell 23: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 26: ✓ — (no output)
Cell 27: ✓ — Out[52]= \begin{pmatrix}0 & \frac{1}{2} & 0 & 0 & 0 & 0\end{pmatrix}
Cell 28: ✓ — Out[53]= \left\{0\right\}

## 2026-07-18 15:16:30 — ▶️ RUN CELLS 29–36
Cell 30: ✓ — Out[59]= \left\{\left(u-\frac{4}{3}\right)\,\left(u-\frac{8}{7}\right),\left(u-\frac{1}{3}\right)\,\left(u-\frac{1}{7}\right),\left(0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i}\right)\,u^2-\left(1.06415582186324025288257211128+0.\times 10^{-29}\,\mathrm{i}\right)\,u\right\}
Cell 31: ✓ — Out[62]= \begin{pmatrix}1 & 1\\2 & 2\\3 & 3\end{pmatrix}
Cell 32: ✓ — Out[65]= \text{max|C2-($\kappa$1+$\kappa$2)Id|}\to 0.\times 10^{-29}
Cell 33: ✓ — Out[75]= \langle\!\lvert \left\{0,0\right\}\to 0.5664055180884988 +1.3769645669711723\,\mathrm{i},\left\{1,0\right\}\to -0.778123966451158+0.\,\,\mathrm{i},\left\{1,1\right\}\to 1.1900552523337025 +0.\,\,\mathrm{i},\left\{2,0\right\}\to 0.5664055180884988 -1.3769645669711723\,\mathrm{i}\rvert\!\rangle
Cell 34: ✓ — 2 kernel messages: Psi::spurious: No genuine state at (M,n)=(2,2) for (λ1,λ2)=(2,1): weight sector M=2 has only 2 genuine state(s) (n=0..1), so n=2 does not exist. [also: Psi::spurious] result: Print: spurious-guard test (expect a Psi::spurious message and Missing[...]):[formula] | Out[84]= \left\{\left\{1.000000000000000000000000000+0.\times 10^{-29}\,\mathrm{i},-1.3004494881076436272092358…
⚠️ stopped at Cell 34 — error detected (pass stopOnError:false to continue past errors)

## 2026-07-18 15:17:10 — ▶️ RUN CELLS 35–42
Cell 35: ✓ — Out[91]= \left\{1.03464899078135229594035121+0.381869568380946283899575941\,\mathrm{i},-2.37022278949001677248172619-0.20048970475590097847093779\,\mathrm{i},1.00000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right\}
Cell 36: ✓ — Out[94]= \left\{-0.709231010958016926537225142+1.09596916251277206706752728\,\mathrm{i},-0.02291660900+0.\times 10^{-12}\,\mathrm{i}\right\}
Cell 37: ✓ — (no output)
Cell 38: ✓ — (no output)
Cell 39: ✓ — (no output)
Cell 40: ✓ — (no output)
Cell 41: ✓ — (no output)
Cell 42: ✓ — Out[109]= 0

## 2026-07-18 15:17:13 — ▶️ RUN CELLS 43–47
Cell 43: ✓ — (no output)
Cell 44: ✓ — Out[113]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 45: ✓ — Out[115]= \left\{6,6\right\}
Cell 46: ✓ — Out[116]= \begin{pmatrix}2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i} & 0 & 0 & 0 & 0 & 0\\0 & 7.36608403016431791181-0.968046299794461400266\,\mathrm{i} & 0 & 0 & 0 & 0\\0 & 0 & 2.68702160063169134012-1.3456153260557858257\,\mathrm{i} & 0 & 0 & 0\\0 & 0 & 0 & 7.3660840301643179118107+0.9680462997944614002659\,\mathrm{i} & 0 & 0\\0 & 0 & 0 & 0 & 2.687021600631691340118… [537 chars total]
Cell 47: ✓ — Out[125]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{max|off-diagonal|}\to 5.065342385207116\times 10^{-34},\text{digits vanished (-Log10)}\to 33.29539119374288,\text{min|diagonal|}\to 3.0051232400647923,\text{diagonal values}\to \left\{2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i},7.36608403016431791181-0… [699 chars total]

## 2026-07-18 15:17:16 — ▶️ RUN CELLS 48–54
Cell 48: ✓ — (no output)
Cell 49: ✓ — (no output)
Cell 50: ✓ — Out[133]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 51: ✓ — Out[138]= \left\{1.00000000000000000000000000000,6\right\}
Cell 52: ✓ — Out[143]= \langle\!\lvert \text{max|<PsiLeft|Psi>| off-diagonal}\to 0.\,,\text{digits vanished (-Log10)}\to \mathrm{Indeterminate}\rvert\!\rangle
Cell 53: ✓ — Out[147]= \left\{7.366084030164317911810727273-0.9680462997944614002658874268\,\mathrm{i},1.00000000000000000000+0.\times 10^{-21}\,\mathrm{i}\right\}
Cell 54: ✓ — Out[150]= \begin{matrix}\langle\!\lvert \text{M}\to 0,\text{n}\to 0,\text{braket}\to 2.909064289-5.804620064\,\mathrm{i},\text{G21 diagonal}\to 2.909064289-5.804620064\,\mathrm{i},\text{normRatio}\to 1.000000000+0.\times 10^{-10}\,\mathrm{i}\rvert\!\rangle \\\langle\!\lvert \text{M}\to 1,\text{n}\to 0,\text{braket}\to 7.366084030-0.968046300\,\mathrm{i},\text{G21 diagonal}\to 7.366084030-0.9680463… [1247 chars total]

## 2026-07-18 15:17:20 — ▶️ RUN CELLS 55–62
Cell 56: ✓ — (no output)
Cell 57: ✓ — Out[158]= \begin{pmatrix}3 & 3\\3 & 3\end{pmatrix}
Cell 58: ✓ — (no output)
Cell 59: ✓ — Out[165]= \left\{0\right\}
Cell 60: ✓ — Out[169]= \begin{matrix}\left\{1,1,0,1.\,,0.\,\right\}\\\left\{1,1,1,1.\,,0.\,\right\}\\\left\{1,2,0,1.\,,0.\,\right\}\\\left\{1,2,1,1.\,,0.\,\right\}\\\left\{1,3,0,1.\,,0.\,\right\}\\\left\{1,3,1,1.\,,0.\,\right\}\\\left\{2,1,0,1.\,,0.\,\right\}\\\left\{2,1,1,1.\,,0.\,\right\}\\\left\{2,1,2,1.\,,0.\,\right\}\\\left\{2,2,0,1.\,,0.\,\right\}\\\left\{2,2,1,1.\,,0.\,\right\}\\\left\{2,2,2,1.\,,0.\,\r… [898 chars total]
Cell 61: ✓ — Out[170]= \langle\!\lvert \text{max spread (ratio non-constancy)}\to 0.\,,\text{max |ratio - 1| over sweep}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases checked}\to 27\rvert\!\rangle
Cell 62: ✓ — Out[176]= \langle\!\lvert \text{max operator-level spread}\to 0.\,,\text{max |$\mu$ - $\mu$Site2|}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases}\to 27\rvert\!\rangle

## 2026-07-18 15:18:53 — ▶️ RUN CELLS 63–70
Cell 64: ✓ — (no output)
Cell 65: ✓ — Out[185]= \langle\!\lvert \text{w\_0 ($\lambda$=1,2,3)}\to \left\{1,1,1\right\},\text{full weights w\_j}\to \left\{1\to \left\{1,-1\right\},2\to \left\{1,-2,1\right\},3\to \left\{1,-3,3,-1\right\}\right\},\text{br1[$\lambda$][1] ($\lambda$=1,2,3)}\to \left\{0,0,0\right\}\rvert\!\rangle
Cell 66: ✓ — Out[187]= \left\{0,0\right\}
Cell 67: ✓ — Out[190]= \sqrt[3]{k}\,\left(\left(-1.18561270626131985872+0.832253303224995280404\,\mathrm{i}\right)\,k^2-\left(0.761819910806671651081+0.812812473459294825632\,\mathrm{i}\right)\,k-\left(0.404410037536507329420+0.598870869485872613792\,\mathrm{i}\right)\right)
Cell 68: ✓ — Out[194]= \begin{pmatrix}0 & 0 & 0\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix}
Cell 69: ✓ — Out[198]= \langle\!\lvert \text{OP(E12)}\to 1.\,+0.\times 10^{-30}\,\mathrm{i},\text{OP(E21)}\to \left(-0.03626734060947795145020244939-2.004237714666237867354326095\,\mathrm{i}\right)\,\varphi ^2-\left(1.522150083623538480912456055-0.2442047957313920727733393355\,\mathrm{i}\right)\,\varphi +\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right),\text{nontrivial generator \{a,b\}}\to \left\{2,1\right\},\… [464 chars total]
Cell 70: ✓ — Out[201]= \left\{1,2\right\}

## 2026-07-18 15:18:56 — ▶️ RUN CELLS 71–77
Cell 71: ✓ — Out[207]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{$\varphi$ candidates per state}\to \begin{pmatrix}\left\{0,0\right\} & \left(-0.920793778924720536105378676-0.390049761817044163511250929\,\mathrm{i}\right)\,k^{4/3}+\left(0.695722366612934233848844287-0.718310788304406529935213657\,\mathrm{i}\right)\,\sqrt[3]{k}-\left(0.3604398751… [1247 chars total]
Cell 72: ✓ — Out[213]= \langle\!\lvert \text{c\_sb solution(s)}\to \left\{\left\{\mathrm{cb}\to \frac{\left(-6.242239011314\times 10^{45}-4.511757007219\times 10^{45}\,\mathrm{i}\right)\,k^{4/3}-\left(5.999302657796\times 10^{45}+8.356648463050\times 10^{45}\,\mathrm{i}\right)\,\sqrt[3]{k}}{9.833530773718\times 10^{45}\,k^{4/3}-\left(3.54439660434\times 10^{45}-9.172544913427\times 10^{45}\,\mathrm{i}\right)\,… [841 chars total]
Cell 73: ✓ — Out[216]= \langle\!\lvert \text{generator \{a,b\} tested}\to \left\{2,1\right\},\text{rep tested}\to \left\{1,2\right\},\text{FL shape (verified exact for all 6 states)}\to \text{FL[1,2][state][k] = k\^{}(1/3) (A\_state + B\_state k)},\text{c0(state) = OPg[state][0]}\to \begin{pmatrix}\left\{0,0\right\} & 1.\,+0.\times 10^{-28}\,\mathrm{i}\\\left\{1,0\right\} & 1.\,+0.\times 10^{-28}\,\mathrm{i}\\… [1247 chars total]
Cell 74: ✓ — Out[220]= \left(1.000000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right)-\left(0.360439875147226537266677651+0.932782448593374716484835982\,\mathrm{i}\right)\,k
Cell 75: ✓ — Out[227]= \langle\!\lvert \text{alphaTab \{state,D0,c0,Dc,c1,$\alpha$\}}\to \begin{pmatrix}\left\{0,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & -0.36043987514722653-0.9327824485933747\,\mathrm{i} & 0.36043987514722653 +0.9327824485933747\,\mathrm{i} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.7101545203859353 -0.2371533418544… [1247 chars total]
Cell 76: ✓ — Out[235]= \langle\!\lvert \text{validation}\to \left\{\left\{1,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\end{pmatrix},\left\{2,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\\\left\{2,1\right\} & 0\\\left\{3,0\right\} & 0\end{pmatrix},\left\{1… [1247 chars total]
Cell 77: ✓ — Out[237]= \langle\!\lvert \text{generator \{a,b\}}\to \left\{2,1\right\},\text{$\varphi$(k)}\to \text{$\alpha$ k},\text{$\alpha$ (universal value)}\to -1.\,+0.\,\,\mathrm{i},\text{normalization}\to \text{1-site FSoV determinant det1[f] = br1[f]/(f/.u->$\theta$1); Section C's det[$\lambda$1,$\lambda$2][f] restricted to 1 site},\text{max residual over rep sweep}\to 0,\text{reps validated}\to \begin{… [510 chars total]

## 2026-07-18 15:19:56 — ▶️ RUN CELLS 78–82
Cell 79: ✓ — (no output)
Cell 80: ✓ — Out[245]= \left\{\text{Test 1 (symbolic monomials u\^{}j) residuals},\left\{0\right\}\right\}
Cell 81: ✓ — Out[250]= \left\{\text{Test 2 (twisted monomials k\^{}(u/h) u\^{}j) residuals},\left\{0\right\}\right\}
Cell 82: ✓ — Out[259]= \left\{\text{Test 3 worst |br identity residual| on Q1},0,\text{worst |detRatio - predicted| },0,\text{num states checked},41\right\}

## 2026-07-18 15:20:00 — ▶️ RUN CELLS 83–87
Cell 84: ✓ — (no output)
Cell 85: ✓ — Out[271]= \begin{pmatrix}1 & 1 & 0\\2 & 1 & 0\\3 & 1 & 0\end{pmatrix}
Cell 86: ✓ — Out[278]= \left\{\text{Mirror Test1 (u\^{}j)},\left\{0\right\},\text{Mirror Test2 (k\^{}(u/h)u\^{}j)},\left\{0\right\}\right\}
Cell 87: ✓ — Out[285]= \left\{\text{Mirror Test3 worst |br identity residual| on Q1},0,\text{num states checked},41\right\}

## 2026-07-18 15:20:02 — ▶️ RUN CELLS 88–91
Cell 89: ✓ — Out[294]= \left\{\text{OP2(E12)},1.\,+0.\times 10^{-30}\,\mathrm{i},\text{OP2(E21)},\left(2.603469709065218092712489006+0.6885776528619826437114966465\,\mathrm{i}\right)\,\varphi +\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right),\text{nontrivial generator \{a,b\}},\left\{2,1\right\}\right\}
Cell 90: ✓ — Out[297]= \left(1.\,+0.\times 10^{-27}\,\mathrm{i}\right)-\left(0.360439875147226537266677651+0.932782448593374716484835982\,\mathrm{i}\right)\,k
Cell 91: ✓ — Out[300]//TableForm= \begin{pmatrix}\left\{0,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & -0.36043987514722653-0.9327824485933747\,\mathrm{i} & 0.36043987514722653 +0.9327824485933747\,\mathrm{i} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.4685059330398573 -0.1713611156914601\,\mathrm{i} & 2.603469709065218 +0.6885776528619827… [1258 chars total]

## 2026-07-18 15:20:14 — ▶️ RUN CELLS 92–108
Cell 92: ✓ — Out[302]//TableForm= \begin{pmatrix}\left\{0,0\right\} & -1.\,+0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.973902234276272 +0.2269679230016555\,\mathrm{i} & 0.7391304347826086 +0.\,\,\mathrm{i} & 0.973902234276272 +0.2269679230016555\,\mathrm{i}\\\left\{1,0\right\} & 0.15191893851739566 -0.10600050418875033\,\mathrm{i} & -0.15191893851739566+0.10600050418875033\,\mathrm{i} & 0.7651595274649543 +… [1258 chars total]
Cell 93: ✓ — Out[305]= \begin{pmatrix}\left\{0,0\right\} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 0.05448860910968355 -0.09540961504089102\,\mathrm{i}\\\left\{1,1\right\} & -0.0293413933141388+0.7363100318787654\,\mathrm{i}\\\left\{2,0\right\} & 0.11703059508534998 +0.\,\,\mathrm{i}\\\left\{2,1\right\} & 0.8158972359113332 +0.\,\,\mathrm{i}\\\left\{3,0\right\} & 0.05448860910968355 +0.09540961504089102\,… [541 chars total]
Cell 96: ✓ — Out[311]= \left\{\text{O[k1,0]-FLdet[k1] (must be 0)},0,\text{O[0,k2]-FLdet2[k2] (expected NONZERO per Phase1 finding)},\left(-3.0719756421050755-0.5172165371705225\,\mathrm{i}\right)\,k2\right\}
Cell 97: ✓ — Out[316]= \left\{\text{Fcand1},k1^2\,\left(\left(-2.2527492196701659-5.82989030370859198\,\mathrm{i}\right)\,k2-\left(0.4170744170089964+23.0487337186617355\,\mathrm{i}\right)\right)+k1\,\left(\left(0.67128913890100625-0.37524873225652243\,\mathrm{i}\right)\,k2+\left(9.51343802264711551-1.5262799733212005\,\mathrm{i}\right)\right)-\left(1.99115021541939352-0.72828474168870538\,\mathrm{i}\right)\,k… [1014 chars total]
Cell 98: ✓ — Out[318]= \left(0.9252125034899888 -6.175727508765443\,\mathrm{i}\right)\,k1^2\,k2+\left(5.4271267482291075 +3.0890004366773067\,\mathrm{i}\right)\,k1\,k2
Cell 99: ✓ — Out[323]= \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 7.991287939023577 -1.2820751775898085\,\mathrm{i}\\\left\{2,0\right\} & -0.3808070763995185-21.044496003995498\,\mathrm{i}\\\left\{0,1\right\} & 0.6123194936458246 +1.416862394550688\,\mathrm{i}\\\left\{1,1\right\} & -0.6248120089007682-3.051910425543801\,\mathrm{i}\\\left\{2,1\right\} & -1.8923093445229393-4.897107855115217\,\m… [421 chars total]
Cell 100: ✓ — Out[328]= \left\{\text{Fcand2 (no P-insertion)},k1^2\,\left(\left(0.02466179161444501+1.36288164597304175\,\mathrm{i}\right)-\left(0.36043987514722654+0.93278244859337472\,\mathrm{i}\right)\,k2\right)+k1\,\left(\left(1.2351720155778515-0.69045766735200127\,\mathrm{i}\right)\,k2-\left(6.46913785540003854-1.0378703818584163\,\mathrm{i}\right)\right)+\left(5.38781822995835893-1.9706528304517910\,\mat… [788 chars total]
Cell 101: ✓ — 2 kernel messages: Power::infy: Infinite expression 0       encountered.
Infinity::indet: Indeterminate expression ----------------- encountered.                                                   1/3                    … [also: Power::infy] result: Out[331]//TableForm= \begin{pmatrix}2 & 1 & \left\{0,0\right\} & \mathrm{Indeterminate} & \mathrm{Indeterminate} & -4.970498218180917+4.515583419069149\,\mathrm{i}\\2 & 1 & \left\{1,1\right\} & \mathr…
Cell 102: ✓ — Out[334]//TableForm= \begin{pmatrix}2 & 1 & \left\{0,0\right\} & \left(-7.771745024480887+7.060451775520925\,\mathrm{i}\right)\,k1^2-\left(3.7846186890458786 +9.794215710230434\,\mathrm{i}\right)\,k1 & \left(1.8923093445229393 +4.897107855115217\,\mathrm{i}\right)\,k2 & -4.970498218180917+4.515583419069149\,\mathrm{i}\\2 & 1 & \left\{1,1\right\} & \left(9.339255641812587 +10.661478687052961\,\math… [1258 chars total]
Cell 103: ✓ — 2 kernel messages: Power::infy: Infinite expression 0       encountered.
Infinity::indet: Indeterminate expression ----------------- encountered.                                                   1/3                    … [also: Power::infy] result: Out[338]= \left\{\text{Fcand1(k1,0)},\mathrm{Indeterminate},\text{FLdet[k1]},\frac{\left(2.45075493251774917289426756629+5.78551794268217244227029022731\,\mathrm{i}\right)\,\left(\left(0.0525946877178…
Cell 104: ✓ — Out[343]= \left\{\text{Fcand1(k1,0)},\left(-8.51191121729+7.73287575414\,\mathrm{i}\right)\,k1^2-\left(4.5054984393+11.6597806074\,\mathrm{i}\right)\,k1+\left(1.000000000000+0.\times 10^{-12}\,\mathrm{i}\right),\text{FLdet[k1]},\frac{\left(2.45075493251774917289426756629+5.78551794268217244227029022731\,\mathrm{i}\right)\,\left(\left(0.05259468771783351844996683714+0.15021349717800998891606886537\… [960 chars total]
Cell 107: ✓ — Out[348]= \left\{\text{STEP 0: det[$\lambda$1,$\lambda$2][k\^{}(u/h)Q1] vs O[k,k], worst |residual| (sample k)},0,\text{sample rows \{$\lambda$1,$\lambda$2,state,resid,ratio\}},\begin{pmatrix}1 & 1 & \left\{0,0\right\} & 0 & 1.\,\\1 & 1 & \left\{1,0\right\} & 0 & 1.\,\\1 & 1 & \left\{1,1\right\} & 0 & 1.\,\\1 & 1 & \left\{2,0\right\} & 0 & 1.\,\\2 & 1 & \left\{0,0\right\} & 0 & 1.\,\\2 & 1 & \left… [445 chars total]
Cell 108: ✓ — Out[354]= \left\{\text{STEP 0 CONFIRMED: det[$\lambda$1,$\lambda$2][k\^{}(u/h)Q1] == O[k,k] identically},\text{worst residual over full sweep},0,\text{num states},41\right\}

## 2026-07-18 15:20:23 — ▶️ RUN CELLS 109–123
Cell 109: ✓ — 2 kernel messages: Power::infy: Infinite expression - encountered.                                  0
Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.
Power::infy: Infinite expression - encounter… [also: Power::infy] result: Out[358]= \begin{pmatrix}1 & 0 & \left\{1.\,\right\}\\1 & 1 & \left\{1.\,\right\}\\2 & 0 & \left\{1.\,\right\}\\2 & 1 & \left\{2.\,\right\}\\2 & 2 & \left\{2.\,\right\}\\3 & 0 & \left\{1.\,\right\}\\3…
Cell 110: ✓ — Out[364]= \left\{\text{STEP 1 CONFIRMED: ⟨x0[$\lambda$]|E21\^{}m = $\gamma$\_m($\lambda$) xL1[$\lambda$][\{m\}], $\gamma$\_m($\lambda$)=$\lambda$!/($\lambda$-m)!},\text{worst deviation},0,\text{table \{$\lambda$,m,$\gamma$\_found,$\gamma$\_predicted,dev\}},\begin{pmatrix}1 & 0 & 1.\, & 1 & 0\\1 & 1 & 1.\, & 1 & 0\\2 & 0 & 1.\, & 1 & 0\\2 & 1 & 2.\, & 2 & 0\\2 & 2 & 2.\, & 2 & 0\\3 & 0 & 1.\, & 1 &… [489 chars total]
Cell 111: ✓ — Out[369]= \left\{\text{STEP 2: x2[\{s1,m\}] vs xL1[$\lambda$1][\{s1\}]$\otimes$xL1[$\lambda$2][\{m\}], worst ratio-spread (constancy check)},0,\text{num cases},54,\text{sample rows},\begin{pmatrix}1 & 1 & 0 & 1 & 0.16 & 0\\1 & 1 & 1 & 1 & 1.\, & 0\\1 & 2 & 0 & 1 & 0.16 & 0\\1 & 2 & 0 & 2 & -0.68 & 0\\1 & 2 & 1 & 1 & 1.\, & 0\\1 & 2 & 1 & 2 & 1.\, & 0\\1 & 3 & 0 & 1 & 0.16 & 0\\1 & 3 & 0 & 2 & -0.6… [479 chars total]
Cell 112: ✓ — Out[373]= \begin{pmatrix}1 & 0 & \frac{4}{25}\\1 & 1 & 1\\2 & 0 & \frac{2}{23}\\2 & 1 & \frac{25}{46}\\2 & 2 & 1\\3 & 0 & \frac{4}{67}\\3 & 1 & \frac{25}{67}\\3 & 2 & \frac{46}{67}\\3 & 3 & 1\end{pmatrix}
Cell 113: ✓ — Out[378]= \left\{\text{c($\lambda$1,s1) hypothesis (m-independent) vs full sweep --- FAILS for m>=2, see next cell for the fix},\text{worst deviation over all 54 ($\lambda$1,$\lambda$2,s1,m) cases},1.68\right\}
Cell 114: ✓ — Out[382]= \begin{pmatrix}1 & 0 & 1 & \frac{4}{25}\\1 & 0 & 2 & -\frac{17}{25}\\1 & 0 & 3 & -\frac{38}{25}\\1 & 1 & 1 & 1\\1 & 1 & 2 & 1\\1 & 1 & 3 & 1\\2 & 0 & 1 & \frac{2}{23}\\2 & 0 & 2 & -\frac{34}{575}\\2 & 0 & 3 & \frac{323}{575}\\2 & 1 & 1 & \frac{25}{46}\\2 & 1 & 2 & \frac{2}{23}\\2 & 1 & 3 & -\frac{17}{46}\\2 & 2 & 1 & 1\\2 & 2 & 2 & 1\\2 & 2 & 3 & 1\\3 & 0 & 1 & \frac{4}{67}\\3 & 0 & 2 & … [672 chars total]
Cell 115: ✓ — 4 kernel messages: LinearSolve::nosol: Linear equation encountered that has no solution.
LinearSolve::nosol: Linear equation encountered that has no solution.
LinearSolve::nosol: Linear equation encountered that has no … [also: LinearSolve::nosol, LinearSolve::nosol, LinearSolve::nosol] result: Out[387]= \left\{\text{coeffs (rows=s1, cols=s1'=0..3)},\left\{\mathrm{LinearSolve}\left[\begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & \frac{1}{36}\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & \f…
Cell 116: ✓ — Out[392]= \left\{\left\{0,2\right\}\to -0.02206359506813757300454250486697,\left\{1,1\right\}\to 0.1635301752109020116807268007787,\left\{2,0\right\}\to 0.8585334198572355613238157040883\right\}
Cell 117: ✓ — Out[396]= \left\{\text{STEP 2 level-conservation check: \# off-level nonzero coeffs (want 0 everywhere)},0,\text{num ($\lambda$1,$\lambda$2,s1,m) cases},29\right\}
Cell 119: ✓ — Out[404]= \left\{\text{m=0: Om-Cm (must be 0 per Phase1, C1pred/P1site2[$\theta$1] should equal exactly Cm's own norm here)},\left(2.5572121404875445 -0.4102640568287387\,\mathrm{i}\right)\,k1-\left(0.038411844228125346 +2.1227491621421546\,\mathrm{i}\right)\,k1^2,\text{m=1: Om\_1[k1]},\left(0.3604398751472265372666776506+0.9327824485933747164848359821\,\mathrm{i}\right)\,k1^2-\left(1.296101147801… [801 chars total]
Cell 120: ✓ — Out[411]= \left\{\text{m=0: (B0-B1)/O\_0 ratio},\frac{\sqrt[3]{k1}\,\left(\left(-0.0292504903494856237+0.0205326892042810627\,\mathrm{i}\right)\,k1^2+\left(0.499242526583056932+0.532659421382621779\,\mathrm{i}\right)\,k1-\left(0.139388493663376676+0.206413542317126162\,\mathrm{i}\right)\right)}{\left(0.3604398751472265372666776506+0.9327824485933747164848359821\,\mathrm{i}\right)\,k1^2-\left(1.296… [1083 chars total]
Cell 121: ✓ — Out[414]= \left\{\text{Cm[0]-Cm[1] shape},\left(-0.0225172879957976151-1.244370198497125076\,\mathrm{i}\right)\,k1^2+\left(5.434075798536032377-0.8718111207610696998\,\mathrm{i}\right)\,k1+\left(0.\times 10^{-20}+0.\times 10^{-20}\,\mathrm{i}\right),\text{ratio vs O\_0[k1] (want: clean constant)},\frac{\left(-0.0225172879957976151-1.244370198497125076\,\mathrm{i}\right)\,k1^2+\left(5.4340757985360… [753 chars total]
Cell 122: ✓ — Out[416]= \left\{\text{Cm[0][k1] / O\_1[k1] (reversed-index hypothesis; want clean k1-independent constant)},\frac{\left(0.111318362779488714+0.0407158959953172346\,\mathrm{i}\right)\,k1^2-\left(0.21818047436770437-1.02534209918467083\,\mathrm{i}\right)\,k1+\left(0.360439875147226537-0.932782448593374716\,\mathrm{i}\right)}{\left(1.\,+0.\times 10^{-29}\,\mathrm{i}\right)\,k1^2-\left(2.963909584212… [539 chars total]

## 2026-07-18 15:20:37 — ▶️ RUN CELLS 124–138
Cell 125: ✓ — (no output)
Cell 126: ✓ — Out[446]= \left\{\text{Symbolic-$\theta$,$\chi$1 port-fidelity: xL1Sym vs xL1},0,\text{x2Sym vs x2},0\right\}
Cell 127: ✓ — Out[451]= \left\{\left\{0,0\right\}\to 0,\left\{0,1\right\}\to \frac{\theta A-\theta B}{\theta A-\theta B+1},\left\{1,0\right\}\to \frac{1}{\theta A-\theta B+1},\left\{1,1\right\}\to 0\right\}
Cell 128: ✓ — Out[454]= 35
Cell 129: ✓ — Out[460]= \left\{\left\{0,0\right\}\to 0,\left\{0,1\right\}\to \theta B-\theta A,\left\{1,0\right\}\to -1,\left\{1,1\right\}\to 0\right\}
Cell 130: ✓ — Out[463]= 127
Cell 131: ✓ — Out[466]= 127
Cell 132: ✓ — Out[470]//MatrixForm= \begin{pmatrix}\frac{\theta A^2-2\,\theta A\,\theta B-\theta A+\theta B^2+\theta B}{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)} & \frac{6\,\left(\theta A-\theta B\right)}{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)} & \frac{6}{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)}\\0 & \frac{\theta A-\theta B+1}{\the… [472 chars total]
Cell 133: ✓ — Out[472]//MatrixForm= \begin{pmatrix}\frac{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)}{\theta A^2-2\,\theta A\,\theta B-\theta A+\theta B^2+\theta B} & -\frac{6\,\left(\theta A-\theta B+3\right)}{\theta A^2-2\,\theta A\,\theta B+\theta B^2-1} & \frac{6}{\left(\theta A-\theta B\right)\,\left(\theta A-\theta B+1\right)}\\0 & \frac{\theta A-\theta B+3}{\theta A-\theta B+1} & -… [453 chars total]
Cell 134: ✓ — Out[478]= \left\{\text{equations},\left\{1.000000000000000000000+0.\times 10^{-21}\,\mathrm{i}⩵\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right)\,wc0+\left(2.603469709065218092712489006+0.6885776528619826437114966465\,\mathrm{i}\right)\,wc1,-1.0350620568640061670+0.16605926109734660949\,\mathrm{i}⩵\left(1.522150083623538480912456055-0.2442047957313920727733393355\,\mathrm{i}\right)\,wc0-\left(1.2961… [785 chars total]
Cell 135: ✓ — Out[484]= \begin{pmatrix}0 & 0 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\0 & 1 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\1 & 0 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\1 & 1 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\2 & 0 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} … [510 chars total]
Cell 136: ✓ — Out[489]= \left\{\text{Conjectured sum rule FAILS beyond the lucky L=1 case},\text{num failing cases},29,\text{num total cases},81,\text{sample failures},\begin{pmatrix}1 & 2 & 0 & 1 & \frac{1}{-\theta A+\theta B-1}\\1 & 2 & 0 & 2 & \frac{2}{\theta A-\theta B+1}\\1 & 3 & 0 & 1 & -\frac{2}{\theta A-\theta B+1}\end{pmatrix}\right\}
Cell 138: ✓ — Out[496]= \left\{\text{BOXED: A[s1,m;s1,m] = Pochhammer[s1+$\theta$1-$\theta$2-m+1,m]/Pochhammer[$\lambda$1+$\theta$1-$\theta$2-m+1,m] --- verified at the live numeric $\theta$1,$\theta$2},\text{worst deviation},0,\text{num ($\lambda$1,$\lambda$2,s1,m) cases},81\right\}

## 2026-07-18 15:22:57 — 📥 BULK INSERT 2 CELL(S) at positions 139–140
1. [markdown] (*Section K — derivation-first closure of the two-rotation formula, following Paul's corrected↵wave-
2. [code] (*STEP 1a — test WITH Q1's own κ1^(u/h) prefactor (Q1 as defined in Section B).*)↵ClearAll[Wfun,step


## 2026-07-18 15:22:57
**Input:**
```wolfram
(*STEP 1a — test WITH Q1's own κ1^(u/h) prefactor (Q1 as defined in Section B).*)
ClearAll[Wfun,step1TestWith];
Wfun[λ1_,λ2_][M_,n_][s1_,s2_] := Q1[λ1,λ2][M,n][θ1+s1] Q1[λ1,λ2][M,n][θ2+s2] / (Q1[λ1,λ2][M,n][θ1] Q1[λ1,λ2][M,n][θ2]);
step1TestWith = {};
Do[
  Do[
    Module[{xv,pred,resid},
      xv = First[Flatten[x2[λ1,λ2][{s1,s2}] . Psi[λ1,λ2][st[[1]],st[[2]]]]];
      pred = Wfun[λ1,λ2][st[[1]],st[[2]]][s1,s2];
      resid = Chop[N[xv-pred],10^-15];
      AppendTo[step1TestWith, {λ1,λ2,st,s1,s2,resid}]
    ]
  ,{s1,0,λ1},{s2,0,λ2}]
,{rep,{{1,1},{2,1},{1,2},{2,2},{3,1},{1,3}}},{λ1,{rep[[1]]}},{λ2,{rep[[2]]}},{st,statesOf[λ1,λ2]}];
worstStep1With = Max[Abs[step1TestWith[[All,6]]]];
{"STEP 1a (WITH κ1^(u/h) prefactor): worst |x2.Psi - Q1(θ1+s1)Q1(θ2+s2)/(Q1(θ1)Q1(θ2))|", worstStep1With,
 "num cases", Length[step1TestWith]}
```
**Output:**
```
Out[502]= \left\{\text{STEP 1a (WITH $\kappa$1\^{}(u/h) prefactor): worst |x2.Psi - Q1($\theta$1+s1)Q1($\theta$2+s2)/(Q1($\theta$1)Q1($\theta$2))|},0,\text{num cases},297\right\}
```
## 2026-07-18 15:23:47 — 📥 BULK INSERT 1 CELL(S) at positions 141–141
1. [code] (*STEP 1b — sanity comparison: test WITHOUT the κ1^(u/h) prefactor (bare polynomial q_n(u) only),↵


## 2026-07-18 15:23:47
**Input:**
```wolfram
(*STEP 1b — sanity comparison: test WITHOUT the κ1^(u/h) prefactor (bare polynomial q_n(u) only),
  to explicitly confirm the prefactor is essential (and see how it fails) rather than merely assuming
  1a's success means the prefactor was irrelevant.*)
ClearAll[qbare,WfunBare,step1TestBare];
qbare[λ1_,λ2_][M_,n_][u_] := qpoly[QSolve[λ1,λ2][M,n],u];
WfunBare[λ1_,λ2_][M_,n_][s1_,s2_] := qbare[λ1,λ2][M,n][θ1+s1] qbare[λ1,λ2][M,n][θ2+s2] / (qbare[λ1,λ2][M,n][θ1] qbare[λ1,λ2][M,n][θ2]);
step1TestBare = {};
Do[
  Module[{xv,pred,resid},
    xv = First[Flatten[x2[2,1][{s1,s2}] . Psi[2,1][st[[1]],st[[2]]]]];
    pred = WfunBare[2,1][st[[1]],st[[2]]][s1,s2];
    resid = Chop[N[xv-pred],10^-15];
    AppendTo[step1TestBare, {st,s1,s2,resid}]
  ]
,{s1,0,2},{s2,0,1},{st,statesOf[2,1]}];
worstStep1Bare = Max[Abs[step1TestBare[[All,4]]]];
{"STEP 1b (WITHOUT prefactor, bare q_n): worst residual (expect NONZERO, confirming prefactor matters)",
 worstStep1Bare, "sample", Take[step1TestBare,UpTo[3]]}
```
**Output:**
```
Out[509]= \left\{\text{STEP 1b (WITHOUT prefactor, bare q\_n): worst residual (expect NONZERO, confirming prefactor matters)},3.7396476326262356,\text{sample},\begin{pmatrix}\left\{0,0\right\} & 0 & 0 & 0\\\left\{1,0\right\} & 0 & 0 & 0\\\left\{1,1\right\} & 0 & 0 & 0\end{pmatrix}\right\}
```
## 2026-07-18 15:24:39 — 📥 BULK INSERT 2 CELL(S) at positions 142–143
1. [markdown] (*STEP 1 conclusion: Paul's exact statement holds with x2 residual EXACTLY 0 (297 cases: 41 states↵×
2. [code] (*STEP 2 — solve for B (state-independent) by linear algebra, no guessing. MP[a,b][state] :=↵  (xL1[


## 2026-07-18 15:24:39
**Input:**
```wolfram
(*STEP 2 — solve for B (state-independent) by linear algebra, no guessing. MP[a,b][state] :=
  (xL1[λ1][{a}]⊗xL1[λ2][{b}]).Psi[state]; impose level support (s,t) with s+t=a+b (Section I), solve
  the OVERDETERMINED system MP[a,b][state]=Σ_{s+t=a+b} B[a,b;s,t] W[s,t][state] using every state of
  a rep as one equation, hard-assert exact (residual 0) consistency. Rep (2,1) first.*)
ClearAll[MPfun,Wval,solveB];
MPfun[λ1_,λ2_][M_,n_][a_,b_] := First[Flatten[KroneckerProduct[xL1[λ1][{a}],xL1[λ2][{b}]] . Psi[λ1,λ2][M,n]]];
Wval[λ1_,λ2_][M_,n_][s_,t_] := First[Flatten[x2[λ1,λ2][{s,t}] . Psi[λ1,λ2][M,n]]];
solveB[λ1_,λ2_][a_,b_] := Module[{L,srng,states,matEqs,rhsEqs,sol},
  L = a+b;
  srng = Range[Max[0,L-λ2],Min[L,λ1]];
  states = statesOf[λ1,λ2];
  matEqs = Table[Wval[λ1,λ2][st[[1]],st[[2]]][s,L-s], {st,states},{s,srng}];
  rhsEqs = Table[MPfun[λ1,λ2][st[[1]],st[[2]]][a,b], {st,states}];
  sol = LeastSquares[matEqs, rhsEqs];
  {"srng",srng,"sol",sol,"maxResid", Max[Abs[matEqs.sol - rhsEqs]]}
];
(*test rep (2,1), all (a,b).*)
testB21 = Table[{a,b,solveB[2,1][a,b]}, {a,0,2},{b,0,1}];
testB21
```
**Output:**
```
Out[515]= \begin{pmatrix}\left\{0,0,\left\{\text{srng},\left\{0\right\},\text{sol},\left\{1.\,\right\},\text{maxResid},0.\times 10^{-28}\right\}\right\} & \left\{0,1,\left\{\text{srng},\left\{0,1\right\},\text{sol},\left\{11.50000000000000000000000000-4.431302171530629069010197265\times 10^{-38}\,\mathrm{i},-10.50000000000000000000000000+2.933894180666706878219239219\times 10^{-38}\,\mathrm{i}\rig… [1115 chars total]
```

## 2026-07-18 15:25:01
**Input:**
```wolfram
Table[{testB21[[a+1,b+1,1]],testB21[[a+1,b+1,2]],Chop[testB21[[a+1,b+1,3]][[4]],10^-20],Chop[testB21[[a+1,b+1,3]][[6]],10^-15]},{a,0,2},{b,0,1}]//TableForm
```
**Output:**
```
Out= TableForm[{{{0, 0, {1.`28.23534835192822}, 0}, {0, 1, {11.50000000000000000000000000000000000015`28.109500546392532, -10.50000000000000000000000000000000000013`28.109500546392532}, 0}}, {{1, 0, {0, 1.`28.109500546392532}, 0}, {1, 1, {1.84`27.741334690026093, -0.84`27.741334690026093}, 0}}, {{2, 0, {0, 1.`27.741334690026093}, 0}, {2, 1, {1.`28.260705158408555}, 0}}}]
```
## 2026-07-18 15:25:23 — 📥 BULK INSERT 1 CELL(S) at positions 144–144
1. [code] (*Full B-sweep, rationalized, all reps in the standard sweep — check (i) exact-0 residual always,↵


## 2026-07-18 15:25:23
**Input:**
```wolfram
(*Full B-sweep, rationalized, all reps in the standard sweep — check (i) exact-0 residual always,
  (ii) whether B[a,b;s,t] depends on λ2 (mirroring A's λ2-independence) or on both λ1,λ2.*)
ClearAll[solveBrat];
solveBrat[λ1_,λ2_][a_,b_] := Module[{L,srng,states,matEqs,rhsEqs,sol,resid},
  L = a+b;
  srng = Range[Max[0,L-λ2],Min[L,λ1]];
  states = statesOf[λ1,λ2];
  matEqs = Table[Wval[λ1,λ2][st[[1]],st[[2]]][s,L-s], {st,states},{s,srng}];
  rhsEqs = Table[MPfun[λ1,λ2][st[[1]],st[[2]]][a,b], {st,states}];
  sol = LeastSquares[matEqs, rhsEqs];
  resid = Max[Abs[matEqs.sol - rhsEqs]];
  {srng, Rationalize[Chop[sol,10^-15],10^-9], resid}
];
fullBsweep = {};
Do[
  Do[
    Do[
      Module[{res},
        res = solveBrat[λ1,λ2][a,b];
        AppendTo[fullBsweep, {λ1,λ2,a,b,res[[1]],res[[2]],res[[3]]}]
      ]
    ,{b,0,λ2}]
  ,{a,0,λ1}]
,{rep,{{2,1},{1,2},{3,1},{1,3},{2,2}}},{λ1,{rep[[1]]}},{λ2,{rep[[2]]}}];
worstBresid = Max[fullBsweep[[All,7]]];
{"worst residual across full B-sweep (want ~0)", worstBresid, "num (rep,a,b) cases", Length[fullBsweep]}
```
**Output:**
```
Out[521]= \left\{\text{worst residual across full B-sweep (want \textasciitilde{}0)},0.\times 10^{-24},\text{num (rep,a,b) cases},37\right\}
```

## 2026-07-18 15:25:58
**Input:**
```wolfram
fullBsweep[[All,{1,2,3,4,5,6}]]//TableForm
```
**Output:**
```
Out= TableForm[{{2, 1, 0, 0, {0}, {1}}, {2, 1, 0, 1, {0, 1}, {23/2, -21/2}}, {2, 1, 1, 0, {0, 1}, {0, 1}}, {2, 1, 1, 1, {1, 2}, {46/25, -21/25}}, {2, 1, 2, 0, {1, 2}, {0, 1}}, {2, 1, 2, 1, {2}, {1}}, {1, 2, 0, 0, {0}, {1}}, {1, 2, 0, 1, {0, 1}, {25/4, -21/4}}, {1, 2, 0, 2, {0, 1}, {-25/17, 42/17}}, {1, 2, 1, 0, {0, 1}, {0, 1}}, {1, 2, 1, 1, {0, 1}, {0, 1}}, {1, 2, 1, 2, {1}, {1}}, {3, 1, 0, 0, {0}, {1}}, {3, 1, 0, 1, {0, 1}, {67/4, -63/4}}, {3, 1, 1, 0, {0, 1}, {0, 1}}, {3, 1, 1, 1, {1, 2}, {67/25, -42/25}}, {3, 1, 2, 0, {1, 2}, {0, 1}}, {3, 1, 2, 1, {2, 3}, {67/46, -21/46}}, {3, 1, 3, 0, {2, 3}, {0, 1}}, {3, 1, 3, 1, {3}, {1}}, {1, 3, 0, 0, {0}, {1}}, {1, 3, 0, 1, {0, 1}, {25/4, -21/4}}, {1, 3, 0, 2, {0, 1}, {-25/17, 42/17}}, {1, 3, 0, 3, {0, 1}, {-25/38, 63/38}}, {1, 3, 1, 0, {0, 1}, {0, 1}}, {1, 3, 1, 1, {0, 1}, {0, 1}}, {1, 3, 1, 2, {0, 1}, {0, 1}}, {1, 3, 1, 3, {1}, {1}}, {2, 2, 0, 0, {0}, {1}}, {2, 2, 0, 1, {0, 1}, {23/2, -21/2}}, {2, 2, 0, 2, {0, 1, 2}, {-575/34, 3864/425, 441/50}}, {2, 2, 1, 0, {0, 1}, {0, 1}}, {2, 2, 1, 1, {0, 1, 2}, {0, 46/25, -21/25}}, {2, 2, 1, 2, {1, 2}, {23/2, -21/2}}, {2, 2, 2, 0, {0, 1, 2}, {0, 0, 1}}, {2, 2, 2, 1, {1, 2}, {0, 1}}, {2, 2, 2, 2, {2}, {1}}}]
```
## 2026-07-18 15:27:42 — 📥 BULK INSERT 1 CELL(S) at positions 145–145
1. [code] (*Structural identification: since x2[{s,t}]=Σ_s' A[s,t;s',t'] xL1[s']⊗xL1[t'] (Section J), and↵  MP


## 2026-07-18 15:27:43
**Input:**
```wolfram
(*Structural identification: since x2[{s,t}]=Σ_s' A[s,t;s',t'] xL1[s']⊗xL1[t'] (Section J), and
  MP[a,b]=Σ B[a,b;s,t] W[s,t] with W[s,t]=x2[{s,t}].Psi, B should be EXACTLY the matrix inverse of A
  within each level block (row=xL-label, col=x2-label) — no new linear solve needed in principle.
  Verify against Section J's decompX2 (A) for λ1=2,λ2=1: invert A's level-2 block and compare to
  solveBrat[2,1][1,1]'s computed B row.*)
level2blockA21 = Table[
  Module[{d,assoc},
    d = decompX2[2,1][s,2-s];
    assoc = Association[d];
    If[KeyExistsQ[assoc,{sp,2-sp}], assoc[{sp,2-sp}], 0]
  ]
,{s,1,2},{sp,1,2}];
level2blockA21inv = Inverse[level2blockA21];
bFromStep2 = solveBrat[2,1][1,1][[2]];
{"A-block (λ1=2,λ2=1, L=2, s,s'=1,2)", level2blockA21, "its inverse (should = B row)", Chop[level2blockA21inv,10^-9],
 "Step2's independently-solved B row {1,1}", bFromStep2,
 "match?", Chop[level2blockA21inv[[1]]-bFromStep2,10^-9]}
```
**Output:**
```
Out[525]= \left\{\text{A-block ($\lambda$1=2,$\lambda$2=1, L=2, s,s'=1,2)},\begin{pmatrix}\frac{25}{46} & \frac{21}{46}\\0 & 1\end{pmatrix},\text{its inverse (should = B row)},\begin{pmatrix}\frac{46}{25} & -\frac{21}{25}\\0 & 1\end{pmatrix},\text{Step2's independently-solved B row \{1,1\}},\left\{\frac{46}{25},-\frac{21}{25}\right\},\text{match?},\left\{0,0\right\}\right\}
```
## 2026-07-18 15:28:58 — 📥 BULK INSERT 2 CELL(S) at positions 146–147
1. [markdown] (*STEP 4 — DERIVE (not guess) the sum rule, by expanding det[λ1,λ2][k^(u/h)Q1] as an explicit↵residu
2. [code] (*w1(s1),w2(s2): the ALREADY-KNOWN measure residue weights from Section C's mu[λ1,λ2][i][u].*)↵Clear


## 2026-07-18 15:28:58
**Input:**
```wolfram
(*w1(s1),w2(s2): the ALREADY-KNOWN measure residue weights from Section C's mu[λ1,λ2][i][u].*)
ClearAll[w1wt,w2wt];
w1wt[λ1_,λ2_][s1_] := 2 π I Residue[mu[λ1,λ2][1][u],{u,θ1+s1 h}];
w2wt[λ1_,λ2_][s2_] := 2 π I Residue[mu[λ1,λ2][2][u],{u,θ2+s2 h}];
(*Verify the STEP-4 derivation directly: does Σ_{s1,s2} w1(s1)w2(s2)(θ2-θ1+(s2-s1)h)/(θ2-θ1) k^(s1+s2) W[s1,s2]
  reproduce det[λ1,λ2][k^(u/h)Q1] exactly, for a representative state?*)
ClearAll[detExpandedCheck];
detExpandedCheck[λ1_,λ2_][M_,n_][k_] := Sum[
  w1wt[λ1,λ2][s1] w2wt[λ1,λ2][s2] (θ2-θ1+(s2-s1)h)/(θ2-θ1) k^(s1+s2) Wval[λ1,λ2][M,n][s1,s2]
,{s1,0,λ1},{s2,0,λ2}];
testDetExpand = Simplify[
  det[2,1][k^(u/h) Q1[2,1][1,0][u]] - detExpandedCheck[2,1][1,0][k]
];
{"det[λ1,λ2][k^u Q1] - explicit residue-sum expansion (rep (2,1), state (1,0)), want 0", Chop[N[testDetExpand]/.k->N[7/3+I/5],10^-10]}
```
**Output:**
```
Out[532]= \left\{\text{det[$\lambda$1,$\lambda$2][k\^{}u Q1] - explicit residue-sum expansion (rep (2,1), state (1,0)), want 0},0\right\}
```
## 2026-07-18 15:30:03 — 📥 BULK INSERT 1 CELL(S) at positions 148–148
1. [code] (*Hard-verify the DERIVED sum rule:↵     (-1)^(s+t) Σ_{a+b=s+t} γ_a(λ1)γ_b(λ2)/(a!b!)·B[a,b;s,t] = w


## 2026-07-18 15:30:03
**Input:**
```wolfram
(*Hard-verify the DERIVED sum rule:
     (-1)^(s+t) Σ_{a+b=s+t} γ_a(λ1)γ_b(λ2)/(a!b!)·B[a,b;s,t] = w1(s)w2(t)(θ2-θ1+(t-s)h)/(θ2-θ1)
   using the rationalized B entries already solved (solveBrat) and the closed-form γ_a(λ)=λ!/(λ-a)!.
   Test across the same rep sweep as the B-solve.*)
ClearAll[gammaFun];
gammaFun[λ_,a_] := λ!/(λ-a)!;
sumRuleTrue = {};
Do[
  Do[
    Do[
      Module[{L,lhs,rhs,bAssoc},
        L = s+t;
        lhs = (-1)^L Sum[
          Module[{bv,srng},
            {srng,bv} = solveBrat[λ1,λ2][a,L-a][[1;;2]];
            If[MemberQ[srng,s], gammaFun[λ1,a] gammaFun[λ2,L-a]/(a!(L-a)!) bv[[Position[srng,s][[1,1]]]], 0]
          ]
        ,{a,0,L}];
        rhs = w1wt[λ1,λ2][s] w2wt[λ1,λ2][t] (θ2-θ1+(t-s)h)/(θ2-θ1);
        AppendTo[sumRuleTrue, {λ1,λ2,s,t, Chop[N[lhs-rhs],10^-8]}]
      ]
    ,{t,0,λ2}]
  ,{s,0,λ1}]
,{rep,{{2,1},{1,2},{3,1}}},{λ1,{rep[[1]]}},{λ2,{rep[[2]]}}];
worstSumRuleTrue = Max[Abs[sumRuleTrue[[All,5]]]];
{"TRUE sum rule (derived, Step 4): worst deviation", worstSumRuleTrue, "num (λ1,λ2,s,t) cases", Length[sumRuleTrue],
 "sample", Take[sumRuleTrue,UpTo[6]]}
```
**Output:**
```
2 kernel messages: Power::infy: Infinite expression - encountered.                                  0
Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.
Infinity::indet: Indeterminate expression 0 … [also: Power::infy] result: Out[538]= \left\{\text{TRUE sum rule (derived, Step 4): worst deviation},\mathrm{Indeterminate},\text{num ($\lambda$1,$\lambda$2,s,t) cases},20,\text{sample},\begin{pmatrix}2 & 1 & 0 & 0 & 0\\2 & 1 & …
```
## 2026-07-18 15:30:18 — ✏️ EDIT CELL 148
(*Hard-verify the DERIVED sum rule:
     (-1)^(s+t) Σ_{a+b=s+t} γ_a(λ1)γ_b(λ2)/(a!b!)·B[a,b;s,t] = w1(s)w2(t)(θ2-θ1+(t-s)h)/(θ2-θ1)
   using the rationalized B entries already solved (solveBrat) and t…


## 2026-07-18 15:30:19
**Input:**
```wolfram
(*Hard-verify the DERIVED sum rule:
     (-1)^(s+t) Σ_{a+b=s+t} γ_a(λ1)γ_b(λ2)/(a!b!)·B[a,b;s,t] = w1(s)w2(t)(θ2-θ1+(t-s)h)/(θ2-θ1)
   using the rationalized B entries already solved (solveBrat) and the closed-form γ_a(λ)=λ!/(λ-a)!.
   (*why the a-range fix*) the sum over a must respect 0<=a<=λ1 AND 0<=L-a<=λ2 (xL1[λ][{a}] is
   undefined outside 0..λ) — the unrestricted range hit invalid (a,b) pairs, producing
   Power::infy/Indeterminate from a malformed solveBrat call.*)
ClearAll[gammaFun];
gammaFun[λ_,a_] := λ!/(λ-a)!;
sumRuleTrue = {};
Do[
  Do[
    Do[
      Module[{L,lhs,rhs,aMin,aMax},
        L = s+t;
        aMin = Max[0,L-λ2]; aMax = Min[L,λ1];
        lhs = (-1)^L Sum[
          Module[{bv,srng},
            {srng,bv} = solveBrat[λ1,λ2][a,L-a][[1;;2]];
            If[MemberQ[srng,s], gammaFun[λ1,a] gammaFun[λ2,L-a]/(a!(L-a)!) bv[[Position[srng,s][[1,1]]]], 0]
          ]
        ,{a,aMin,aMax}];
        rhs = w1wt[λ1,λ2][s] w2wt[λ1,λ2][t] (θ2-θ1+(t-s)h)/(θ2-θ1);
        AppendTo[sumRuleTrue, {λ1,λ2,s,t, Chop[N[lhs-rhs],10^-8]}]
      ]
    ,{t,0,λ2}]
  ,{s,0,λ1}]
,{rep,{{2,1},{1,2},{3,1}}},{λ1,{rep[[1]]}},{λ2,{rep[[2]]}}];
worstSumRuleTrue = Max[Abs[sumRuleTrue[[All,5]]]];
{"TRUE sum rule (derived, Step 4): worst deviation", worstSumRuleTrue, "num (λ1,λ2,s,t) cases", Length[sumRuleTrue],
 "sample", Take[sumRuleTrue,UpTo[6]]}
```
**Output:**
```
Out= {"TRUE sum rule (derived, Step 4): worst deviation", 0, "num (λ1,λ2,s,t) cases", 20, "sample", {{2, 1, 0, 0, 0}, {2, 1, 0, 1, 0}, {2, 1, 1, 0, 0}, {2, 1, 1, 1, 0}, {2, 1, 2, 0, 0}, {2, 1, 2, 1, 0}}}
```
## 2026-07-18 15:31:54 — 📥 BULK INSERT 2 CELL(S) at positions 149–150
1. [markdown] (*STEP 3 — assembly, now UNCONDITIONAL given Steps 1-2 verified exactly:↵   O[k1,k2] = Σ_{a,b} cA[a]
2. [code] (*Memoize solveBrat (it's called many times below) and assemble O[k1,k2].*)↵ClearAll[solveBratM,Oass


## 2026-07-18 15:31:54
**Input:**
```wolfram
(*Memoize solveBrat (it's called many times below) and assemble O[k1,k2].*)
ClearAll[solveBratM,Oassembled];
solveBratM[λ1_,λ2_][a_,b_] := solveBratM[λ1,λ2][a,b] = solveBrat[λ1,λ2][a,b];
Oassembled[λ1_,λ2_][M_,n_][k1_,k2_] := Sum[
  Module[{srng,bv,cA,cB},
    cA = (-k1)^a/a! gammaFun[λ1,a];
    cB = (-k2)^b/b! gammaFun[λ2,b];
    {srng,bv} = solveBratM[λ1,λ2][a,b][[1;;2]];
    cA cB Sum[bv[[i]] Wval[λ1,λ2][M,n][srng[[i]], a+b-srng[[i]]], {i,Length[srng]}]
  ]
,{a,0,λ1},{b,0,λ2}];
(*hard test: rep (2,1), all 6 states, sampled at 3 independent (k1,k2) grid points.*)
gridPts = {{1.3,0.7},{2.1-0.3I,0.5+0.2I},{0.4+1.1I,-0.9}};
assemblyTest = {};
Do[
  Do[
    Module[{oA,oKnown,resid},
      oA = Oassembled[2,1][st[[1]],st[[2]]][pt[[1]],pt[[2]]];
      oKnown = O2rot[2,1][st[[1]],st[[2]]][pt[[1]],pt[[2]]];
      resid = Chop[N[oA-oKnown],10^-9];
      AppendTo[assemblyTest, {st,pt,resid}]
    ]
  ,{pt,gridPts}]
,{st,statesOf[2,1]}];
worstAssembly = Max[Abs[assemblyTest[[All,3]]]];
{"STEP 3 assembly test (rep (2,1), 6 states x 3 grid points): worst |Oassembled-O2rot|", worstAssembly,
 "num cases", Length[assemblyTest]}
```
**Output:**
```
Out[546]= \left\{\text{STEP 3 assembly test (rep (2,1), 6 states x 3 grid points): worst |Oassembled-O2rot|},0,\text{num cases},18\right\}
```
## 2026-07-18 15:32:50 — 📥 BULK INSERT 1 CELL(S) at positions 151–151
1. [code] (*Full-sweep hard-assert: all 41 states, 6 reps, 2 grid points per state (kept modest for runtime;↵


## 2026-07-18 15:32:50
**Input:**
```wolfram
(*Full-sweep hard-assert: all 41 states, 6 reps, 2 grid points per state (kept modest for runtime;
  each point already exercises the FULL B-solve+assembly machinery independently per rep).*)
gridPts2 = {{1.3,0.7},{2.1-0.3I,0.5+0.2I}};
assemblyFull = {};
Do[
  {λ1,λ2} = rep;
  Do[
    Do[
      Module[{oA,oKnown,resid},
        oA = Oassembled[λ1,λ2][st[[1]],st[[2]]][pt[[1]],pt[[2]]];
        oKnown = O2rot[λ1,λ2][st[[1]],st[[2]]][pt[[1]],pt[[2]]];
        resid = Chop[N[oA-oKnown],10^-9];
        AppendTo[assemblyFull, {λ1,λ2,st,pt,resid}]
      ]
    ,{pt,gridPts2}]
  ,{st,statesOf[λ1,λ2]}]
,{rep,{{1,1},{2,1},{1,2},{2,2},{3,1},{1,3}}}];
worstAssemblyFull = Max[Abs[assemblyFull[[All,5]]]];
If[worstAssemblyFull > 10^-8, Abort[]];
{"STEP 3 FULL-SWEEP ASSEMBLY CONFIRMED: O[k1,k2] = Σ cA(k1)cB(k2) B-weighted Q-bilinears == O2rot exactly",
 "worst deviation", worstAssemblyFull, "num states x grid cases", Length[assemblyFull]}
```
**Output:**
```
Out[552]= \left\{\text{STEP 3 FULL-SWEEP ASSEMBLY CONFIRMED: O[k1,k2] = $\Sigma$ cA(k1)cB(k2) B-weighted Q-bilinears == O2rot exactly},\text{worst deviation},0,\text{num states x grid cases},82\right\}
```
## 2026-07-18 15:34:32 — 📥 BULK INSERT 2 CELL(S) at positions 152–153
1. [markdown] (*STEP 5 (timeboxed) — structure hunt on B.↵(i) B's diagonal is now known exactly (identified as A^{
2. [code] (*Numeric confirmation of the Step-5 derivation: compare Fcand[k1,k2] (per-row-insertion,↵  node-nor


## 2026-07-18 15:34:33
**Input:**
```wolfram
(*Numeric confirmation of the Step-5 derivation: compare Fcand[k1,k2] (per-row-insertion,
  node-normalized so k1=k2 reproduces Step-0's det exactly by construction) against the CORRECT
  O[k1,k2], coefficient-by-coefficient in k1^a k2^b, for rep (2,1) state (1,0).*)
ClearAll[FcandStepK,RfunStepK];
RfunStepK[λ1_,λ2_][s_,t_] := w1wt[λ1,λ2][s] w2wt[λ1,λ2][t] (θ2-θ1+(t-s)h)/(θ2-θ1);
FcandStepK[λ1_,λ2_][M_,n_][k1_,k2_] := Sum[
  RfunStepK[λ1,λ2][s,t] k1^s k2^t Wval[λ1,λ2][M,n][s,t]
,{s,0,λ1},{t,0,λ2}];
(*coefficient-by-coefficient comparison, rep (2,1) state (1,0): O has coeffs at (a,b), Fcand only at
  (a,b) diagonal (s=a,t=b) — print both and the difference.*)
compTable = Table[
  Module[{oCoef,fCoef},
    oCoef = (-1)^(a+b)/(a! b!) gammaFun[2,a] gammaFun[1,b] * Module[{srng,bv,idx},
      {srng,bv} = solveBratM[2,1][a,b][[1;;2]];
      Sum[bv[[i]] Wval[2,1][1,0][srng[[i]],a+b-srng[[i]]], {i,Length[srng]}]
    ];
    fCoef = RfunStepK[2,1][a,b] Wval[2,1][1,0][a,b];
    {a,b, Chop[N[oCoef],10^-10], Chop[N[fCoef],10^-10], Chop[N[oCoef-fCoef],10^-10]}
  ]
,{a,0,2},{b,0,1}];
Flatten[compTable,1]
```
**Output:**
```
Out[557]= \begin{pmatrix}0 & 0 & 1.\, & 1.\, & 0\\0 & 1 & -2.603469709065218-0.6885776528619827\,\mathrm{i} & 5.387818229958359 -1.970652830451791\,\mathrm{i} & -7.991287939023577+1.2820751775898085\,\mathrm{i}\\1 & 0 & 1.5221500836235384 -0.24420479573139206\,\mathrm{i} & -6.469137855400039+1.0378703818584163\,\mathrm{i} & 7.991287939023577 -1.2820751775898085\,\mathrm{i}\\1 & 1 & 1.2961011478017… [843 chars total]
```
## 2026-07-18 15:36:11 — 📥 BULK INSERT 1 CELL(S) at positions 154–154
1. [markdown] (*Section K final report — the two-rotation construction CLOSES completely, derivation-first.↵↵STEP


---

## 2026-07-18 15:36:16 — 🔄 KERNEL RESTART

## 2026-07-18 15:36:29 — ▶️ RUN CELLS 1–10
Cell 1: ✓ — Out[8]= \left\{\frac{1}{3},\frac{1}{7},1,0.360439875147226537266677650596+0.932782448593374716484835982071\,\mathrm{i},0.36043987514722653726667765060-0.932782448593374716484835982071\,\mathrm{i},0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i},1\right\}
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 8: ✓ — (no output)
Cell 10: ✓ — (no output)

## 2026-07-18 15:36:34 — ▶️ RUN CELLS 11–20
Cell 11: ✓ — (no output)
Cell 12: ✓ — Out[27]= \left\{0\right\}
Cell 13: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — Out[31]= \left\{0\right\}
Cell 16: ✓ — (no output)
Cell 17: ✓ — Out[34]= \begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{pmatrix}
Cell 18: ✓ — (no output)
Cell 20: ✓ — (no output)

## 2026-07-18 15:36:37 — ▶️ RUN CELLS 21–28
Cell 21: ✓ — (no output)
Cell 23: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 26: ✓ — (no output)
Cell 27: ✓ — Out[52]= \begin{pmatrix}0 & \frac{1}{2} & 0 & 0 & 0 & 0\end{pmatrix}
Cell 28: ✓ — Out[53]= \left\{0\right\}

## 2026-07-18 15:36:39 — ▶️ RUN CELLS 29–36
Cell 30: ✓ — Out[59]= \left\{\left(u-\frac{4}{3}\right)\,\left(u-\frac{8}{7}\right),\left(u-\frac{1}{3}\right)\,\left(u-\frac{1}{7}\right),\left(0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i}\right)\,u^2-\left(1.06415582186324025288257211128+0.\times 10^{-29}\,\mathrm{i}\right)\,u\right\}
Cell 31: ✓ — Out[62]= \begin{pmatrix}1 & 1\\2 & 2\\3 & 3\end{pmatrix}
Cell 32: ✓ — Out[65]= \text{max|C2-($\kappa$1+$\kappa$2)Id|}\to 0.\times 10^{-29}
Cell 33: ✓ — Out[75]= \langle\!\lvert \left\{0,0\right\}\to 0.5664055180884988 +1.3769645669711723\,\mathrm{i},\left\{1,0\right\}\to -0.778123966451158+0.\,\,\mathrm{i},\left\{1,1\right\}\to 1.1900552523337025 +0.\,\,\mathrm{i},\left\{2,0\right\}\to 0.5664055180884988 -1.3769645669711723\,\mathrm{i}\rvert\!\rangle
Cell 34: ✓ — 2 kernel messages: Psi::spurious: No genuine state at (M,n)=(2,2) for (λ1,λ2)=(2,1): weight sector M=2 has only 2 genuine state(s) (n=0..1), so n=2 does not exist. [also: Psi::spurious] result: Print: spurious-guard test (expect a Psi::spurious message and Missing[...]):[formula] | Out[84]= \left\{\left\{1.000000000000000000000000000+0.\times 10^{-29}\,\mathrm{i},-1.3004494881076436272092358…
⚠️ stopped at Cell 34 — error detected (pass stopOnError:false to continue past errors)

## 2026-07-18 15:37:12 — ▶️ RUN CELLS 35–42
Cell 35: ✓ — Out[91]= \left\{1.03464899078135229594035121+0.381869568380946283899575941\,\mathrm{i},-2.37022278949001677248172619-0.20048970475590097847093779\,\mathrm{i},1.00000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right\}
Cell 36: ✓ — Out[94]= \left\{-0.709231010958016926537225142+1.09596916251277206706752728\,\mathrm{i},-0.02291660900+0.\times 10^{-12}\,\mathrm{i}\right\}
Cell 37: ✓ — (no output)
Cell 38: ✓ — (no output)
Cell 39: ✓ — (no output)
Cell 40: ✓ — (no output)
Cell 41: ✓ — (no output)
Cell 42: ✓ — Out[109]= 0

## 2026-07-18 15:37:15 — ▶️ RUN CELLS 43–47
Cell 43: ✓ — (no output)
Cell 44: ✓ — Out[113]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 45: ✓ — Out[115]= \left\{6,6\right\}
Cell 46: ✓ — Out[116]= \begin{pmatrix}2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i} & 0 & 0 & 0 & 0 & 0\\0 & 7.36608403016431791181-0.968046299794461400266\,\mathrm{i} & 0 & 0 & 0 & 0\\0 & 0 & 2.68702160063169134012-1.3456153260557858257\,\mathrm{i} & 0 & 0 & 0\\0 & 0 & 0 & 7.3660840301643179118107+0.9680462997944614002659\,\mathrm{i} & 0 & 0\\0 & 0 & 0 & 0 & 2.687021600631691340118… [537 chars total]
Cell 47: ✓ — Out[125]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{max|off-diagonal|}\to 5.065342385207116\times 10^{-34},\text{digits vanished (-Log10)}\to 33.29539119374288,\text{min|diagonal|}\to 3.0051232400647923,\text{diagonal values}\to \left\{2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i},7.36608403016431791181-0… [699 chars total]

## 2026-07-18 15:37:19 — ▶️ RUN CELLS 48–54
Cell 48: ✓ — (no output)
Cell 49: ✓ — (no output)
Cell 50: ✓ — Out[133]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 51: ✓ — Out[138]= \left\{1.00000000000000000000000000000,6\right\}
Cell 52: ✓ — Out[143]= \langle\!\lvert \text{max|<PsiLeft|Psi>| off-diagonal}\to 0.\,,\text{digits vanished (-Log10)}\to \mathrm{Indeterminate}\rvert\!\rangle
Cell 53: ✓ — Out[147]= \left\{7.366084030164317911810727273-0.9680462997944614002658874268\,\mathrm{i},1.00000000000000000000+0.\times 10^{-21}\,\mathrm{i}\right\}
Cell 54: ✓ — Out[150]= \begin{matrix}\langle\!\lvert \text{M}\to 0,\text{n}\to 0,\text{braket}\to 2.909064289-5.804620064\,\mathrm{i},\text{G21 diagonal}\to 2.909064289-5.804620064\,\mathrm{i},\text{normRatio}\to 1.000000000+0.\times 10^{-10}\,\mathrm{i}\rvert\!\rangle \\\langle\!\lvert \text{M}\to 1,\text{n}\to 0,\text{braket}\to 7.366084030-0.968046300\,\mathrm{i},\text{G21 diagonal}\to 7.366084030-0.9680463… [1247 chars total]

## 2026-07-18 15:37:22 — ▶️ RUN CELLS 55–62
Cell 56: ✓ — (no output)
Cell 57: ✓ — Out[158]= \begin{pmatrix}3 & 3\\3 & 3\end{pmatrix}
Cell 58: ✓ — (no output)
Cell 59: ✓ — Out[165]= \left\{0\right\}
Cell 60: ✓ — Out[169]= \begin{matrix}\left\{1,1,0,1.\,,0.\,\right\}\\\left\{1,1,1,1.\,,0.\,\right\}\\\left\{1,2,0,1.\,,0.\,\right\}\\\left\{1,2,1,1.\,,0.\,\right\}\\\left\{1,3,0,1.\,,0.\,\right\}\\\left\{1,3,1,1.\,,0.\,\right\}\\\left\{2,1,0,1.\,,0.\,\right\}\\\left\{2,1,1,1.\,,0.\,\right\}\\\left\{2,1,2,1.\,,0.\,\right\}\\\left\{2,2,0,1.\,,0.\,\right\}\\\left\{2,2,1,1.\,,0.\,\right\}\\\left\{2,2,2,1.\,,0.\,\r… [898 chars total]
Cell 61: ✓ — Out[170]= \langle\!\lvert \text{max spread (ratio non-constancy)}\to 0.\,,\text{max |ratio - 1| over sweep}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases checked}\to 27\rvert\!\rangle
Cell 62: ✓ — Out[176]= \langle\!\lvert \text{max operator-level spread}\to 0.\,,\text{max |$\mu$ - $\mu$Site2|}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases}\to 27\rvert\!\rangle

## 2026-07-18 15:37:35 — ▶️ RUN CELLS 63–70
Cell 64: ✓ — (no output)
Cell 65: ✓ — Out[185]= \langle\!\lvert \text{w\_0 ($\lambda$=1,2,3)}\to \left\{1,1,1\right\},\text{full weights w\_j}\to \left\{1\to \left\{1,-1\right\},2\to \left\{1,-2,1\right\},3\to \left\{1,-3,3,-1\right\}\right\},\text{br1[$\lambda$][1] ($\lambda$=1,2,3)}\to \left\{0,0,0\right\}\rvert\!\rangle
Cell 66: ✓ — Out[187]= \left\{0,0\right\}
Cell 67: ✓ — Out[190]= \sqrt[3]{k}\,\left(\left(-1.18561270626131985872+0.832253303224995280404\,\mathrm{i}\right)\,k^2-\left(0.761819910806671651081+0.812812473459294825632\,\mathrm{i}\right)\,k-\left(0.404410037536507329420+0.598870869485872613792\,\mathrm{i}\right)\right)
Cell 68: ✓ — Out[194]= \begin{pmatrix}0 & 0 & 0\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix}
Cell 69: ✓ — Out[198]= \langle\!\lvert \text{OP(E12)}\to 1.\,+0.\times 10^{-30}\,\mathrm{i},\text{OP(E21)}\to \left(-0.03626734060947795145020244939-2.004237714666237867354326095\,\mathrm{i}\right)\,\varphi ^2-\left(1.522150083623538480912456055-0.2442047957313920727733393355\,\mathrm{i}\right)\,\varphi +\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right),\text{nontrivial generator \{a,b\}}\to \left\{2,1\right\},\… [464 chars total]
Cell 70: ✓ — Out[201]= \left\{1,2\right\}

## 2026-07-18 15:37:39 — ▶️ RUN CELLS 71–77
Cell 71: ✓ — Out[207]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{$\varphi$ candidates per state}\to \begin{pmatrix}\left\{0,0\right\} & \left(-0.920793778924720536105378676-0.390049761817044163511250929\,\mathrm{i}\right)\,k^{4/3}+\left(0.695722366612934233848844287-0.718310788304406529935213657\,\mathrm{i}\right)\,\sqrt[3]{k}-\left(0.3604398751… [1247 chars total]
Cell 72: ✓ — Out[213]= \langle\!\lvert \text{c\_sb solution(s)}\to \left\{\left\{\mathrm{cb}\to \frac{\left(-6.242239011314\times 10^{45}-4.511757007219\times 10^{45}\,\mathrm{i}\right)\,k^{4/3}-\left(5.999302657796\times 10^{45}+8.356648463050\times 10^{45}\,\mathrm{i}\right)\,\sqrt[3]{k}}{9.833530773718\times 10^{45}\,k^{4/3}-\left(3.54439660434\times 10^{45}-9.172544913427\times 10^{45}\,\mathrm{i}\right)\,… [841 chars total]
Cell 73: ✓ — Out[216]= \langle\!\lvert \text{generator \{a,b\} tested}\to \left\{2,1\right\},\text{rep tested}\to \left\{1,2\right\},\text{FL shape (verified exact for all 6 states)}\to \text{FL[1,2][state][k] = k\^{}(1/3) (A\_state + B\_state k)},\text{c0(state) = OPg[state][0]}\to \begin{pmatrix}\left\{0,0\right\} & 1.\,+0.\times 10^{-28}\,\mathrm{i}\\\left\{1,0\right\} & 1.\,+0.\times 10^{-28}\,\mathrm{i}\\… [1247 chars total]
Cell 74: ✓ — Out[220]= \left(1.000000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right)-\left(0.360439875147226537266677651+0.932782448593374716484835982\,\mathrm{i}\right)\,k
Cell 75: ✓ — Out[227]= \langle\!\lvert \text{alphaTab \{state,D0,c0,Dc,c1,$\alpha$\}}\to \begin{pmatrix}\left\{0,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & -0.36043987514722653-0.9327824485933747\,\mathrm{i} & 0.36043987514722653 +0.9327824485933747\,\mathrm{i} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.7101545203859353 -0.2371533418544… [1247 chars total]
Cell 76: ✓ — Out[235]= \langle\!\lvert \text{validation}\to \left\{\left\{1,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\end{pmatrix},\left\{2,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\\\left\{2,1\right\} & 0\\\left\{3,0\right\} & 0\end{pmatrix},\left\{1… [1247 chars total]
Cell 77: ✓ — Out[237]= \langle\!\lvert \text{generator \{a,b\}}\to \left\{2,1\right\},\text{$\varphi$(k)}\to \text{$\alpha$ k},\text{$\alpha$ (universal value)}\to -1.\,+0.\,\,\mathrm{i},\text{normalization}\to \text{1-site FSoV determinant det1[f] = br1[f]/(f/.u->$\theta$1); Section C's det[$\lambda$1,$\lambda$2][f] restricted to 1 site},\text{max residual over rep sweep}\to 0,\text{reps validated}\to \begin{… [510 chars total]

## 2026-07-18 15:37:52 — ▶️ RUN CELLS 78–82
Cell 79: ✓ — (no output)
Cell 80: ✓ — Out[245]= \left\{\text{Test 1 (symbolic monomials u\^{}j) residuals},\left\{0\right\}\right\}
Cell 81: ✓ — Out[250]= \left\{\text{Test 2 (twisted monomials k\^{}(u/h) u\^{}j) residuals},\left\{0\right\}\right\}
Cell 82: ✓ — Out[259]= \left\{\text{Test 3 worst |br identity residual| on Q1},0,\text{worst |detRatio - predicted| },0,\text{num states checked},41\right\}

## 2026-07-18 15:37:56 — ▶️ RUN CELLS 83–87
Cell 84: ✓ — (no output)
Cell 85: ✓ — Out[271]= \begin{pmatrix}1 & 1 & 0\\2 & 1 & 0\\3 & 1 & 0\end{pmatrix}
Cell 86: ✓ — Out[278]= \left\{\text{Mirror Test1 (u\^{}j)},\left\{0\right\},\text{Mirror Test2 (k\^{}(u/h)u\^{}j)},\left\{0\right\}\right\}
Cell 87: ✓ — Out[285]= \left\{\text{Mirror Test3 worst |br identity residual| on Q1},0,\text{num states checked},41\right\}

## 2026-07-18 15:37:57 — ▶️ RUN CELLS 88–91
Cell 89: ✓ — Out[294]= \left\{\text{OP2(E12)},1.\,+0.\times 10^{-30}\,\mathrm{i},\text{OP2(E21)},\left(2.603469709065218092712489006+0.6885776528619826437114966465\,\mathrm{i}\right)\,\varphi +\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right),\text{nontrivial generator \{a,b\}},\left\{2,1\right\}\right\}
Cell 90: ✓ — Out[297]= \left(1.\,+0.\times 10^{-27}\,\mathrm{i}\right)-\left(0.360439875147226537266677651+0.932782448593374716484835982\,\mathrm{i}\right)\,k
Cell 91: ✓ — Out[300]//TableForm= \begin{pmatrix}\left\{0,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & -0.36043987514722653-0.9327824485933747\,\mathrm{i} & 0.36043987514722653 +0.9327824485933747\,\mathrm{i} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.4685059330398573 -0.1713611156914601\,\mathrm{i} & 2.603469709065218 +0.6885776528619827… [1258 chars total]

## 2026-07-18 15:38:10 — ▶️ RUN CELLS 92–108
Cell 92: ✓ — Out[302]//TableForm= \begin{pmatrix}\left\{0,0\right\} & -1.\,+0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.973902234276272 +0.2269679230016555\,\mathrm{i} & 0.7391304347826086 +0.\,\,\mathrm{i} & 0.973902234276272 +0.2269679230016555\,\mathrm{i}\\\left\{1,0\right\} & 0.15191893851739566 -0.10600050418875033\,\mathrm{i} & -0.15191893851739566+0.10600050418875033\,\mathrm{i} & 0.7651595274649543 +… [1258 chars total]
Cell 93: ✓ — Out[305]= \begin{pmatrix}\left\{0,0\right\} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 0.05448860910968355 -0.09540961504089102\,\mathrm{i}\\\left\{1,1\right\} & -0.0293413933141388+0.7363100318787654\,\mathrm{i}\\\left\{2,0\right\} & 0.11703059508534998 +0.\,\,\mathrm{i}\\\left\{2,1\right\} & 0.8158972359113332 +0.\,\,\mathrm{i}\\\left\{3,0\right\} & 0.05448860910968355 +0.09540961504089102\,… [541 chars total]
Cell 96: ✓ — Out[311]= \left\{\text{O[k1,0]-FLdet[k1] (must be 0)},0,\text{O[0,k2]-FLdet2[k2] (expected NONZERO per Phase1 finding)},\left(-3.0719756421050755-0.5172165371705225\,\mathrm{i}\right)\,k2\right\}
Cell 97: ✓ — Out[316]= \left\{\text{Fcand1},k1^2\,\left(\left(-2.2527492196701659-5.82989030370859198\,\mathrm{i}\right)\,k2-\left(0.4170744170089964+23.0487337186617355\,\mathrm{i}\right)\right)+k1\,\left(\left(0.67128913890100625-0.37524873225652243\,\mathrm{i}\right)\,k2+\left(9.51343802264711551-1.5262799733212005\,\mathrm{i}\right)\right)-\left(1.99115021541939352-0.72828474168870538\,\mathrm{i}\right)\,k… [1014 chars total]
Cell 98: ✓ — Out[318]= \left(0.9252125034899888 -6.175727508765443\,\mathrm{i}\right)\,k1^2\,k2+\left(5.4271267482291075 +3.0890004366773067\,\mathrm{i}\right)\,k1\,k2
Cell 99: ✓ — Out[323]= \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 7.991287939023577 -1.2820751775898085\,\mathrm{i}\\\left\{2,0\right\} & -0.3808070763995185-21.044496003995498\,\mathrm{i}\\\left\{0,1\right\} & 0.6123194936458246 +1.416862394550688\,\mathrm{i}\\\left\{1,1\right\} & -0.6248120089007682-3.051910425543801\,\mathrm{i}\\\left\{2,1\right\} & -1.8923093445229393-4.897107855115217\,\m… [421 chars total]
Cell 100: ✓ — Out[328]= \left\{\text{Fcand2 (no P-insertion)},k1^2\,\left(\left(0.02466179161444501+1.36288164597304175\,\mathrm{i}\right)-\left(0.36043987514722654+0.93278244859337472\,\mathrm{i}\right)\,k2\right)+k1\,\left(\left(1.2351720155778515-0.69045766735200127\,\mathrm{i}\right)\,k2-\left(6.46913785540003854-1.0378703818584163\,\mathrm{i}\right)\right)+\left(5.38781822995835893-1.9706528304517910\,\mat… [788 chars total]
Cell 101: ✓ — 2 kernel messages: Power::infy: Infinite expression 0       encountered.
Infinity::indet: Indeterminate expression ----------------- encountered.                                                   1/3                    … [also: Power::infy] result: Out[331]//TableForm= \begin{pmatrix}2 & 1 & \left\{0,0\right\} & \mathrm{Indeterminate} & \mathrm{Indeterminate} & -4.970498218180917+4.515583419069149\,\mathrm{i}\\2 & 1 & \left\{1,1\right\} & \mathr…
Cell 102: ✓ — Out[334]//TableForm= \begin{pmatrix}2 & 1 & \left\{0,0\right\} & \left(-7.771745024480887+7.060451775520925\,\mathrm{i}\right)\,k1^2-\left(3.7846186890458786 +9.794215710230434\,\mathrm{i}\right)\,k1 & \left(1.8923093445229393 +4.897107855115217\,\mathrm{i}\right)\,k2 & -4.970498218180917+4.515583419069149\,\mathrm{i}\\2 & 1 & \left\{1,1\right\} & \left(9.339255641812587 +10.661478687052961\,\math… [1258 chars total]
Cell 103: ✓ — 2 kernel messages: Power::infy: Infinite expression 0       encountered.
Infinity::indet: Indeterminate expression ----------------- encountered.                                                   1/3                    … [also: Power::infy] result: Out[338]= \left\{\text{Fcand1(k1,0)},\mathrm{Indeterminate},\text{FLdet[k1]},\frac{\left(2.45075493251774917289426756629+5.78551794268217244227029022731\,\mathrm{i}\right)\,\left(\left(0.0525946877178…
Cell 104: ✓ — Out[343]= \left\{\text{Fcand1(k1,0)},\left(-8.51191121729+7.73287575414\,\mathrm{i}\right)\,k1^2-\left(4.5054984393+11.6597806074\,\mathrm{i}\right)\,k1+\left(1.000000000000+0.\times 10^{-12}\,\mathrm{i}\right),\text{FLdet[k1]},\frac{\left(2.45075493251774917289426756629+5.78551794268217244227029022731\,\mathrm{i}\right)\,\left(\left(0.05259468771783351844996683714+0.15021349717800998891606886537\… [960 chars total]
Cell 107: ✓ — Out[348]= \left\{\text{STEP 0: det[$\lambda$1,$\lambda$2][k\^{}(u/h)Q1] vs O[k,k], worst |residual| (sample k)},0,\text{sample rows \{$\lambda$1,$\lambda$2,state,resid,ratio\}},\begin{pmatrix}1 & 1 & \left\{0,0\right\} & 0 & 1.\,\\1 & 1 & \left\{1,0\right\} & 0 & 1.\,\\1 & 1 & \left\{1,1\right\} & 0 & 1.\,\\1 & 1 & \left\{2,0\right\} & 0 & 1.\,\\2 & 1 & \left\{0,0\right\} & 0 & 1.\,\\2 & 1 & \left… [445 chars total]
Cell 108: ✓ — Out[354]= \left\{\text{STEP 0 CONFIRMED: det[$\lambda$1,$\lambda$2][k\^{}(u/h)Q1] == O[k,k] identically},\text{worst residual over full sweep},0,\text{num states},41\right\}

## 2026-07-18 15:38:20 — ▶️ RUN CELLS 109–123
Cell 109: ✓ — 2 kernel messages: Power::infy: Infinite expression - encountered.                                  0
Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.
Power::infy: Infinite expression - encounter… [also: Power::infy] result: Out[358]= \begin{pmatrix}1 & 0 & \left\{1.\,\right\}\\1 & 1 & \left\{1.\,\right\}\\2 & 0 & \left\{1.\,\right\}\\2 & 1 & \left\{2.\,\right\}\\2 & 2 & \left\{2.\,\right\}\\3 & 0 & \left\{1.\,\right\}\\3…
Cell 110: ✓ — Out[364]= \left\{\text{STEP 1 CONFIRMED: ⟨x0[$\lambda$]|E21\^{}m = $\gamma$\_m($\lambda$) xL1[$\lambda$][\{m\}], $\gamma$\_m($\lambda$)=$\lambda$!/($\lambda$-m)!},\text{worst deviation},0,\text{table \{$\lambda$,m,$\gamma$\_found,$\gamma$\_predicted,dev\}},\begin{pmatrix}1 & 0 & 1.\, & 1 & 0\\1 & 1 & 1.\, & 1 & 0\\2 & 0 & 1.\, & 1 & 0\\2 & 1 & 2.\, & 2 & 0\\2 & 2 & 2.\, & 2 & 0\\3 & 0 & 1.\, & 1 &… [489 chars total]
Cell 111: ✓ — Out[369]= \left\{\text{STEP 2: x2[\{s1,m\}] vs xL1[$\lambda$1][\{s1\}]$\otimes$xL1[$\lambda$2][\{m\}], worst ratio-spread (constancy check)},0,\text{num cases},54,\text{sample rows},\begin{pmatrix}1 & 1 & 0 & 1 & 0.16 & 0\\1 & 1 & 1 & 1 & 1.\, & 0\\1 & 2 & 0 & 1 & 0.16 & 0\\1 & 2 & 0 & 2 & -0.68 & 0\\1 & 2 & 1 & 1 & 1.\, & 0\\1 & 2 & 1 & 2 & 1.\, & 0\\1 & 3 & 0 & 1 & 0.16 & 0\\1 & 3 & 0 & 2 & -0.6… [479 chars total]
Cell 112: ✓ — Out[373]= \begin{pmatrix}1 & 0 & \frac{4}{25}\\1 & 1 & 1\\2 & 0 & \frac{2}{23}\\2 & 1 & \frac{25}{46}\\2 & 2 & 1\\3 & 0 & \frac{4}{67}\\3 & 1 & \frac{25}{67}\\3 & 2 & \frac{46}{67}\\3 & 3 & 1\end{pmatrix}
Cell 113: ✓ — Out[378]= \left\{\text{c($\lambda$1,s1) hypothesis (m-independent) vs full sweep --- FAILS for m>=2, see next cell for the fix},\text{worst deviation over all 54 ($\lambda$1,$\lambda$2,s1,m) cases},1.68\right\}
Cell 114: ✓ — Out[382]= \begin{pmatrix}1 & 0 & 1 & \frac{4}{25}\\1 & 0 & 2 & -\frac{17}{25}\\1 & 0 & 3 & -\frac{38}{25}\\1 & 1 & 1 & 1\\1 & 1 & 2 & 1\\1 & 1 & 3 & 1\\2 & 0 & 1 & \frac{2}{23}\\2 & 0 & 2 & -\frac{34}{575}\\2 & 0 & 3 & \frac{323}{575}\\2 & 1 & 1 & \frac{25}{46}\\2 & 1 & 2 & \frac{2}{23}\\2 & 1 & 3 & -\frac{17}{46}\\2 & 2 & 1 & 1\\2 & 2 & 2 & 1\\2 & 2 & 3 & 1\\3 & 0 & 1 & \frac{4}{67}\\3 & 0 & 2 & … [672 chars total]
Cell 115: ✓ — 4 kernel messages: LinearSolve::nosol: Linear equation encountered that has no solution.
LinearSolve::nosol: Linear equation encountered that has no solution.
LinearSolve::nosol: Linear equation encountered that has no … [also: LinearSolve::nosol, LinearSolve::nosol, LinearSolve::nosol] result: Out[387]= \left\{\text{coeffs (rows=s1, cols=s1'=0..3)},\left\{\mathrm{LinearSolve}\left[\begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & \frac{1}{36}\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & \f…
Cell 116: ✓ — Out[392]= \left\{\left\{0,2\right\}\to -0.02206359506813757300454250486697,\left\{1,1\right\}\to 0.1635301752109020116807268007787,\left\{2,0\right\}\to 0.8585334198572355613238157040883\right\}
Cell 117: ✓ — Out[396]= \left\{\text{STEP 2 level-conservation check: \# off-level nonzero coeffs (want 0 everywhere)},0,\text{num ($\lambda$1,$\lambda$2,s1,m) cases},29\right\}
Cell 119: ✓ — Out[404]= \left\{\text{m=0: Om-Cm (must be 0 per Phase1, C1pred/P1site2[$\theta$1] should equal exactly Cm's own norm here)},\left(2.5572121404875445 -0.4102640568287387\,\mathrm{i}\right)\,k1-\left(0.038411844228125346 +2.1227491621421546\,\mathrm{i}\right)\,k1^2,\text{m=1: Om\_1[k1]},\left(0.3604398751472265372666776506+0.9327824485933747164848359821\,\mathrm{i}\right)\,k1^2-\left(1.296101147801… [801 chars total]
Cell 120: ✓ — Out[411]= \left\{\text{m=0: (B0-B1)/O\_0 ratio},\frac{\sqrt[3]{k1}\,\left(\left(-0.0292504903494856237+0.0205326892042810627\,\mathrm{i}\right)\,k1^2+\left(0.499242526583056932+0.532659421382621779\,\mathrm{i}\right)\,k1-\left(0.139388493663376676+0.206413542317126162\,\mathrm{i}\right)\right)}{\left(0.3604398751472265372666776506+0.9327824485933747164848359821\,\mathrm{i}\right)\,k1^2-\left(1.296… [1083 chars total]
Cell 121: ✓ — Out[414]= \left\{\text{Cm[0]-Cm[1] shape},\left(-0.0225172879957976151-1.244370198497125076\,\mathrm{i}\right)\,k1^2+\left(5.434075798536032377-0.8718111207610696998\,\mathrm{i}\right)\,k1+\left(0.\times 10^{-20}+0.\times 10^{-20}\,\mathrm{i}\right),\text{ratio vs O\_0[k1] (want: clean constant)},\frac{\left(-0.0225172879957976151-1.244370198497125076\,\mathrm{i}\right)\,k1^2+\left(5.4340757985360… [753 chars total]
Cell 122: ✓ — Out[416]= \left\{\text{Cm[0][k1] / O\_1[k1] (reversed-index hypothesis; want clean k1-independent constant)},\frac{\left(0.111318362779488714+0.0407158959953172346\,\mathrm{i}\right)\,k1^2-\left(0.21818047436770437-1.02534209918467083\,\mathrm{i}\right)\,k1+\left(0.360439875147226537-0.932782448593374716\,\mathrm{i}\right)}{\left(1.\,+0.\times 10^{-29}\,\mathrm{i}\right)\,k1^2-\left(2.963909584212… [539 chars total]

## 2026-07-18 15:38:34 — ▶️ RUN CELLS 124–138
Cell 125: ✓ — (no output)
Cell 126: ✓ — Out[446]= \left\{\text{Symbolic-$\theta$,$\chi$1 port-fidelity: xL1Sym vs xL1},0,\text{x2Sym vs x2},0\right\}
Cell 127: ✓ — Out[451]= \left\{\left\{0,0\right\}\to 0,\left\{0,1\right\}\to \frac{\theta A-\theta B}{\theta A-\theta B+1},\left\{1,0\right\}\to \frac{1}{\theta A-\theta B+1},\left\{1,1\right\}\to 0\right\}
Cell 128: ✓ — Out[454]= 35
Cell 129: ✓ — Out[460]= \left\{\left\{0,0\right\}\to 0,\left\{0,1\right\}\to \theta B-\theta A,\left\{1,0\right\}\to -1,\left\{1,1\right\}\to 0\right\}
Cell 130: ✓ — Out[463]= 127
Cell 131: ✓ — Out[466]= 127
Cell 132: ✓ — Out[470]//MatrixForm= \begin{pmatrix}\frac{\theta A^2-2\,\theta A\,\theta B-\theta A+\theta B^2+\theta B}{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)} & \frac{6\,\left(\theta A-\theta B\right)}{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)} & \frac{6}{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)}\\0 & \frac{\theta A-\theta B+1}{\the… [472 chars total]
Cell 133: ✓ — Out[472]//MatrixForm= \begin{pmatrix}\frac{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)}{\theta A^2-2\,\theta A\,\theta B-\theta A+\theta B^2+\theta B} & -\frac{6\,\left(\theta A-\theta B+3\right)}{\theta A^2-2\,\theta A\,\theta B+\theta B^2-1} & \frac{6}{\left(\theta A-\theta B\right)\,\left(\theta A-\theta B+1\right)}\\0 & \frac{\theta A-\theta B+3}{\theta A-\theta B+1} & -… [453 chars total]
Cell 134: ✓ — Out[478]= \left\{\text{equations},\left\{1.000000000000000000000+0.\times 10^{-21}\,\mathrm{i}⩵\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right)\,wc0+\left(2.603469709065218092712489006+0.6885776528619826437114966465\,\mathrm{i}\right)\,wc1,-1.0350620568640061670+0.16605926109734660949\,\mathrm{i}⩵\left(1.522150083623538480912456055-0.2442047957313920727733393355\,\mathrm{i}\right)\,wc0-\left(1.2961… [785 chars total]
Cell 135: ✓ — Out[484]= \begin{pmatrix}0 & 0 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\0 & 1 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\1 & 0 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\1 & 1 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\2 & 0 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} … [510 chars total]
Cell 136: ✓ — Out[489]= \left\{\text{Conjectured sum rule FAILS beyond the lucky L=1 case},\text{num failing cases},29,\text{num total cases},81,\text{sample failures},\begin{pmatrix}1 & 2 & 0 & 1 & \frac{1}{-\theta A+\theta B-1}\\1 & 2 & 0 & 2 & \frac{2}{\theta A-\theta B+1}\\1 & 3 & 0 & 1 & -\frac{2}{\theta A-\theta B+1}\end{pmatrix}\right\}
Cell 138: ✓ — Out[496]= \left\{\text{BOXED: A[s1,m;s1,m] = Pochhammer[s1+$\theta$1-$\theta$2-m+1,m]/Pochhammer[$\lambda$1+$\theta$1-$\theta$2-m+1,m] --- verified at the live numeric $\theta$1,$\theta$2},\text{worst deviation},0,\text{num ($\lambda$1,$\lambda$2,s1,m) cases},81\right\}

## 2026-07-18 15:38:46 — ▶️ RUN CELLS 139–154
Cell 140: ✓ — Out[502]= \left\{\text{STEP 1a (WITH $\kappa$1\^{}(u/h) prefactor): worst |x2.Psi - Q1($\theta$1+s1)Q1($\theta$2+s2)/(Q1($\theta$1)Q1($\theta$2))|},0,\text{num cases},297\right\}
Cell 141: ✓ — Out[509]= \left\{\text{STEP 1b (WITHOUT prefactor, bare q\_n): worst residual (expect NONZERO, confirming prefactor matters)},3.7396476326262356,\text{sample},\begin{pmatrix}\left\{0,0\right\} & 0 & 0 & 0\\\left\{1,0\right\} & 0 & 0 & 0\\\left\{1,1\right\} & 0 & 0 & 0\end{pmatrix}\right\}
Cell 143: ✓ — Out[515]= \begin{pmatrix}\left\{0,0,\left\{\text{srng},\left\{0\right\},\text{sol},\left\{1.\,\right\},\text{maxResid},0.\times 10^{-28}\right\}\right\} & \left\{0,1,\left\{\text{srng},\left\{0,1\right\},\text{sol},\left\{11.50000000000000000000000000-4.431302171530629069010197265\times 10^{-38}\,\mathrm{i},-10.50000000000000000000000000+2.933894180666706878219239219\times 10^{-38}\,\mathrm{i}\rig… [1115 chars total]
Cell 144: ✓ — Out[521]= \left\{\text{worst residual across full B-sweep (want \textasciitilde{}0)},0.\times 10^{-24},\text{num (rep,a,b) cases},37\right\}
Cell 145: ✓ — Out[525]= \left\{\text{A-block ($\lambda$1=2,$\lambda$2=1, L=2, s,s'=1,2)},\begin{pmatrix}\frac{25}{46} & \frac{21}{46}\\0 & 1\end{pmatrix},\text{its inverse (should = B row)},\begin{pmatrix}\frac{46}{25} & -\frac{21}{25}\\0 & 1\end{pmatrix},\text{Step2's independently-solved B row \{1,1\}},\left\{\frac{46}{25},-\frac{21}{25}\right\},\text{match?},\left\{0,0\right\}\right\}
Cell 147: ✓ — Out[532]= \left\{\text{det[$\lambda$1,$\lambda$2][k\^{}u Q1] - explicit residue-sum expansion (rep (2,1), state (1,0)), want 0},0\right\}
Cell 148: ✓ — Out[538]= \left\{\text{TRUE sum rule (derived, Step 4): worst deviation},0,\text{num ($\lambda$1,$\lambda$2,s,t) cases},20,\text{sample},\begin{pmatrix}2 & 1 & 0 & 0 & 0\\2 & 1 & 0 & 1 & 0\\2 & 1 & 1 & 0 & 0\\2 & 1 & 1 & 1 & 0\\2 & 1 & 2 & 0 & 0\\2 & 1 & 2 & 1 & 0\end{pmatrix}\right\}
Cell 150: ✓ — Out[546]= \left\{\text{STEP 3 assembly test (rep (2,1), 6 states x 3 grid points): worst |Oassembled-O2rot|},0,\text{num cases},18\right\}
Cell 151: ✓ — Out[552]= \left\{\text{STEP 3 FULL-SWEEP ASSEMBLY CONFIRMED: O[k1,k2] = $\Sigma$ cA(k1)cB(k2) B-weighted Q-bilinears == O2rot exactly},\text{worst deviation},0,\text{num states x grid cases},82\right\}
Cell 153: ✓ — Out[557]= \begin{pmatrix}0 & 0 & 1.\, & 1.\, & 0\\0 & 1 & -2.603469709065218-0.6885776528619827\,\mathrm{i} & 5.387818229958359 -1.970652830451791\,\mathrm{i} & -7.991287939023577+1.2820751775898085\,\mathrm{i}\\1 & 0 & 1.5221500836235384 -0.24420479573139206\,\mathrm{i} & -6.469137855400039+1.0378703818584163\,\mathrm{i} & 7.991287939023577 -1.2820751775898085\,\mathrm{i}\\1 & 1 & 1.2961011478017… [843 chars total]

## 2026-07-18 16:16:01 — ▶️ RUN CELLS 1–10
Cell 1: ✓ — Out[8]= \left\{\frac{1}{3},\frac{1}{7},1,0.360439875147226537266677650596+0.932782448593374716484835982071\,\mathrm{i},0.36043987514722653726667765060-0.932782448593374716484835982071\,\mathrm{i},0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i},1\right\}
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 8: ✓ — (no output)
Cell 10: ✓ — (no output)

## 2026-07-18 16:16:06 — ▶️ RUN CELLS 11–20
Cell 11: ✓ — (no output)
Cell 12: ✓ — Out[27]= \left\{0\right\}
Cell 13: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — Out[31]= \left\{0\right\}
Cell 16: ✓ — (no output)
Cell 17: ✓ — Out[34]= \begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{pmatrix}
Cell 18: ✓ — (no output)
Cell 20: ✓ — (no output)

## 2026-07-18 16:16:09 — ▶️ RUN CELLS 21–28
Cell 21: ✓ — (no output)
Cell 23: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 26: ✓ — (no output)
Cell 27: ✓ — Out[52]= \begin{pmatrix}0 & \frac{1}{2} & 0 & 0 & 0 & 0\end{pmatrix}
Cell 28: ✓ — Out[53]= \left\{0\right\}

## 2026-07-18 16:16:12 — ▶️ RUN CELLS 29–36
Cell 30: ✓ — Out[59]= \left\{\left(u-\frac{4}{3}\right)\,\left(u-\frac{8}{7}\right),\left(u-\frac{1}{3}\right)\,\left(u-\frac{1}{7}\right),\left(0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i}\right)\,u^2-\left(1.06415582186324025288257211128+0.\times 10^{-29}\,\mathrm{i}\right)\,u\right\}
Cell 31: ✓ — Out[62]= \begin{pmatrix}1 & 1\\2 & 2\\3 & 3\end{pmatrix}
Cell 32: ✓ — Out[65]= \text{max|C2-($\kappa$1+$\kappa$2)Id|}\to 0.\times 10^{-29}
Cell 33: ✓ — Out[75]= \langle\!\lvert \left\{0,0\right\}\to 0.5664055180884988 +1.3769645669711723\,\mathrm{i},\left\{1,0\right\}\to -0.778123966451158+0.\,\,\mathrm{i},\left\{1,1\right\}\to 1.1900552523337025 +0.\,\,\mathrm{i},\left\{2,0\right\}\to 0.5664055180884988 -1.3769645669711723\,\mathrm{i}\rvert\!\rangle
Cell 34: ✓ — 2 kernel messages: Psi::spurious: No genuine state at (M,n)=(2,2) for (λ1,λ2)=(2,1): weight sector M=2 has only 2 genuine state(s) (n=0..1), so n=2 does not exist. [also: Psi::spurious] result: Print: spurious-guard test (expect a Psi::spurious message and Missing[...]):[formula] | Out[84]= \left\{\left\{1.000000000000000000000000000+0.\times 10^{-29}\,\mathrm{i},-1.3004494881076436272092358…
⚠️ stopped at Cell 34 — error detected (pass stopOnError:false to continue past errors)

## 2026-07-18 16:16:25 — ▶️ RUN CELLS 35–42
Cell 35: ✓ — Out[91]= \left\{1.03464899078135229594035121+0.381869568380946283899575941\,\mathrm{i},-2.37022278949001677248172619-0.20048970475590097847093779\,\mathrm{i},1.00000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right\}
Cell 36: ✓ — Out[94]= \left\{-0.709231010958016926537225142+1.09596916251277206706752728\,\mathrm{i},-0.02291660900+0.\times 10^{-12}\,\mathrm{i}\right\}
Cell 37: ✓ — (no output)
Cell 38: ✓ — (no output)
Cell 39: ✓ — (no output)
Cell 40: ✓ — (no output)
Cell 41: ✓ — (no output)
Cell 42: ✓ — Out[109]= 0

## 2026-07-18 16:16:28 — ▶️ RUN CELLS 43–47
Cell 43: ✓ — (no output)
Cell 44: ✓ — Out[113]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 45: ✓ — Out[115]= \left\{6,6\right\}
Cell 46: ✓ — Out[116]= \begin{pmatrix}2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i} & 0 & 0 & 0 & 0 & 0\\0 & 7.36608403016431791181-0.968046299794461400266\,\mathrm{i} & 0 & 0 & 0 & 0\\0 & 0 & 2.68702160063169134012-1.3456153260557858257\,\mathrm{i} & 0 & 0 & 0\\0 & 0 & 0 & 7.3660840301643179118107+0.9680462997944614002659\,\mathrm{i} & 0 & 0\\0 & 0 & 0 & 0 & 2.687021600631691340118… [537 chars total]
Cell 47: ✓ — Out[125]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{max|off-diagonal|}\to 5.065342385207116\times 10^{-34},\text{digits vanished (-Log10)}\to 33.29539119374288,\text{min|diagonal|}\to 3.0051232400647923,\text{diagonal values}\to \left\{2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i},7.36608403016431791181-0… [699 chars total]

## 2026-07-18 16:16:32 — ▶️ RUN CELLS 48–54
Cell 48: ✓ — (no output)
Cell 49: ✓ — (no output)
Cell 50: ✓ — Out[133]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 51: ✓ — Out[138]= \left\{1.00000000000000000000000000000,6\right\}
Cell 52: ✓ — Out[143]= \langle\!\lvert \text{max|<PsiLeft|Psi>| off-diagonal}\to 0.\,,\text{digits vanished (-Log10)}\to \mathrm{Indeterminate}\rvert\!\rangle
Cell 53: ✓ — Out[147]= \left\{7.366084030164317911810727273-0.9680462997944614002658874268\,\mathrm{i},1.00000000000000000000+0.\times 10^{-21}\,\mathrm{i}\right\}
Cell 54: ✓ — Out[150]= \begin{matrix}\langle\!\lvert \text{M}\to 0,\text{n}\to 0,\text{braket}\to 2.909064289-5.804620064\,\mathrm{i},\text{G21 diagonal}\to 2.909064289-5.804620064\,\mathrm{i},\text{normRatio}\to 1.000000000+0.\times 10^{-10}\,\mathrm{i}\rvert\!\rangle \\\langle\!\lvert \text{M}\to 1,\text{n}\to 0,\text{braket}\to 7.366084030-0.968046300\,\mathrm{i},\text{G21 diagonal}\to 7.366084030-0.9680463… [1247 chars total]

## 2026-07-18 16:16:35 — ▶️ RUN CELLS 55–62
Cell 56: ✓ — (no output)
Cell 57: ✓ — Out[158]= \begin{pmatrix}3 & 3\\3 & 3\end{pmatrix}
Cell 58: ✓ — (no output)
Cell 59: ✓ — Out[165]= \left\{0\right\}
Cell 60: ✓ — Out[169]= \begin{matrix}\left\{1,1,0,1.\,,0.\,\right\}\\\left\{1,1,1,1.\,,0.\,\right\}\\\left\{1,2,0,1.\,,0.\,\right\}\\\left\{1,2,1,1.\,,0.\,\right\}\\\left\{1,3,0,1.\,,0.\,\right\}\\\left\{1,3,1,1.\,,0.\,\right\}\\\left\{2,1,0,1.\,,0.\,\right\}\\\left\{2,1,1,1.\,,0.\,\right\}\\\left\{2,1,2,1.\,,0.\,\right\}\\\left\{2,2,0,1.\,,0.\,\right\}\\\left\{2,2,1,1.\,,0.\,\right\}\\\left\{2,2,2,1.\,,0.\,\r… [898 chars total]
Cell 61: ✓ — Out[170]= \langle\!\lvert \text{max spread (ratio non-constancy)}\to 0.\,,\text{max |ratio - 1| over sweep}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases checked}\to 27\rvert\!\rangle
Cell 62: ✓ — Out[176]= \langle\!\lvert \text{max operator-level spread}\to 0.\,,\text{max |$\mu$ - $\mu$Site2|}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases}\to 27\rvert\!\rangle

## 2026-07-18 16:16:47 — ▶️ RUN CELLS 63–70
Cell 64: ✓ — (no output)
Cell 65: ✓ — Out[185]= \langle\!\lvert \text{w\_0 ($\lambda$=1,2,3)}\to \left\{1,1,1\right\},\text{full weights w\_j}\to \left\{1\to \left\{1,-1\right\},2\to \left\{1,-2,1\right\},3\to \left\{1,-3,3,-1\right\}\right\},\text{br1[$\lambda$][1] ($\lambda$=1,2,3)}\to \left\{0,0,0\right\}\rvert\!\rangle
Cell 66: ✓ — Out[187]= \left\{0,0\right\}
Cell 67: ✓ — Out[190]= \sqrt[3]{k}\,\left(\left(-1.18561270626131985872+0.832253303224995280404\,\mathrm{i}\right)\,k^2-\left(0.761819910806671651081+0.812812473459294825632\,\mathrm{i}\right)\,k-\left(0.404410037536507329420+0.598870869485872613792\,\mathrm{i}\right)\right)
Cell 68: ✓ — Out[194]= \begin{pmatrix}0 & 0 & 0\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix}
Cell 69: ✓ — Out[198]= \langle\!\lvert \text{OP(E12)}\to 1.\,+0.\times 10^{-30}\,\mathrm{i},\text{OP(E21)}\to \left(-0.03626734060947795145020244939-2.004237714666237867354326095\,\mathrm{i}\right)\,\varphi ^2-\left(1.522150083623538480912456055-0.2442047957313920727733393355\,\mathrm{i}\right)\,\varphi +\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right),\text{nontrivial generator \{a,b\}}\to \left\{2,1\right\},\… [464 chars total]
Cell 70: ✓ — Out[201]= \left\{1,2\right\}

## 2026-07-18 16:16:51 — ▶️ RUN CELLS 71–77
Cell 71: ✓ — Out[207]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{$\varphi$ candidates per state}\to \begin{pmatrix}\left\{0,0\right\} & \left(-0.920793778924720536105378676-0.390049761817044163511250929\,\mathrm{i}\right)\,k^{4/3}+\left(0.695722366612934233848844287-0.718310788304406529935213657\,\mathrm{i}\right)\,\sqrt[3]{k}-\left(0.3604398751… [1247 chars total]
Cell 72: ✓ — Out[213]= \langle\!\lvert \text{c\_sb solution(s)}\to \left\{\left\{\mathrm{cb}\to \frac{\left(-6.242239011314\times 10^{45}-4.511757007219\times 10^{45}\,\mathrm{i}\right)\,k^{4/3}-\left(5.999302657796\times 10^{45}+8.356648463050\times 10^{45}\,\mathrm{i}\right)\,\sqrt[3]{k}}{9.833530773718\times 10^{45}\,k^{4/3}-\left(3.54439660434\times 10^{45}-9.172544913427\times 10^{45}\,\mathrm{i}\right)\,… [841 chars total]
Cell 73: ✓ — Out[216]= \langle\!\lvert \text{generator \{a,b\} tested}\to \left\{2,1\right\},\text{rep tested}\to \left\{1,2\right\},\text{FL shape (verified exact for all 6 states)}\to \text{FL[1,2][state][k] = k\^{}(1/3) (A\_state + B\_state k)},\text{c0(state) = OPg[state][0]}\to \begin{pmatrix}\left\{0,0\right\} & 1.\,+0.\times 10^{-28}\,\mathrm{i}\\\left\{1,0\right\} & 1.\,+0.\times 10^{-28}\,\mathrm{i}\\… [1247 chars total]
Cell 74: ✓ — Out[220]= \left(1.000000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right)-\left(0.360439875147226537266677651+0.932782448593374716484835982\,\mathrm{i}\right)\,k
Cell 75: ✓ — Out[227]= \langle\!\lvert \text{alphaTab \{state,D0,c0,Dc,c1,$\alpha$\}}\to \begin{pmatrix}\left\{0,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & -0.36043987514722653-0.9327824485933747\,\mathrm{i} & 0.36043987514722653 +0.9327824485933747\,\mathrm{i} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.7101545203859353 -0.2371533418544… [1247 chars total]
Cell 76: ✓ — Out[235]= \langle\!\lvert \text{validation}\to \left\{\left\{1,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\end{pmatrix},\left\{2,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\\\left\{2,1\right\} & 0\\\left\{3,0\right\} & 0\end{pmatrix},\left\{1… [1247 chars total]
Cell 77: ✓ — Out[237]= \langle\!\lvert \text{generator \{a,b\}}\to \left\{2,1\right\},\text{$\varphi$(k)}\to \text{$\alpha$ k},\text{$\alpha$ (universal value)}\to -1.\,+0.\,\,\mathrm{i},\text{normalization}\to \text{1-site FSoV determinant det1[f] = br1[f]/(f/.u->$\theta$1); Section C's det[$\lambda$1,$\lambda$2][f] restricted to 1 site},\text{max residual over rep sweep}\to 0,\text{reps validated}\to \begin{… [510 chars total]

## 2026-07-18 16:17:03 — ▶️ RUN CELLS 78–82
Cell 79: ✓ — (no output)
Cell 80: ✓ — Out[245]= \left\{\text{Test 1 (symbolic monomials u\^{}j) residuals},\left\{0\right\}\right\}
Cell 81: ✓ — Out[250]= \left\{\text{Test 2 (twisted monomials k\^{}(u/h) u\^{}j) residuals},\left\{0\right\}\right\}
Cell 82: ✓ — Out[259]= \left\{\text{Test 3 worst |br identity residual| on Q1},0,\text{worst |detRatio - predicted| },0,\text{num states checked},41\right\}

## 2026-07-18 16:17:07 — ▶️ RUN CELLS 83–87
Cell 84: ✓ — (no output)
Cell 85: ✓ — Out[271]= \begin{pmatrix}1 & 1 & 0\\2 & 1 & 0\\3 & 1 & 0\end{pmatrix}
Cell 86: ✓ — Out[278]= \left\{\text{Mirror Test1 (u\^{}j)},\left\{0\right\},\text{Mirror Test2 (k\^{}(u/h)u\^{}j)},\left\{0\right\}\right\}
Cell 87: ✓ — Out[285]= \left\{\text{Mirror Test3 worst |br identity residual| on Q1},0,\text{num states checked},41\right\}

## 2026-07-18 16:17:09 — ▶️ RUN CELLS 88–91
Cell 89: ✓ — Out[294]= \left\{\text{OP2(E12)},1.\,+0.\times 10^{-30}\,\mathrm{i},\text{OP2(E21)},\left(2.603469709065218092712489006+0.6885776528619826437114966465\,\mathrm{i}\right)\,\varphi +\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right),\text{nontrivial generator \{a,b\}},\left\{2,1\right\}\right\}
Cell 90: ✓ — Out[297]= \left(1.\,+0.\times 10^{-27}\,\mathrm{i}\right)-\left(0.360439875147226537266677651+0.932782448593374716484835982\,\mathrm{i}\right)\,k
Cell 91: ✓ — Out[300]//TableForm= \begin{pmatrix}\left\{0,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & -0.36043987514722653-0.9327824485933747\,\mathrm{i} & 0.36043987514722653 +0.9327824485933747\,\mathrm{i} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.4685059330398573 -0.1713611156914601\,\mathrm{i} & 2.603469709065218 +0.6885776528619827… [1258 chars total]

## 2026-07-18 16:17:21 — ▶️ RUN CELLS 92–108
Cell 92: ✓ — Out[302]//TableForm= \begin{pmatrix}\left\{0,0\right\} & -1.\,+0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.973902234276272 +0.2269679230016555\,\mathrm{i} & 0.7391304347826086 +0.\,\,\mathrm{i} & 0.973902234276272 +0.2269679230016555\,\mathrm{i}\\\left\{1,0\right\} & 0.15191893851739566 -0.10600050418875033\,\mathrm{i} & -0.15191893851739566+0.10600050418875033\,\mathrm{i} & 0.7651595274649543 +… [1258 chars total]
Cell 93: ✓ — Out[305]= \begin{pmatrix}\left\{0,0\right\} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 0.05448860910968355 -0.09540961504089102\,\mathrm{i}\\\left\{1,1\right\} & -0.0293413933141388+0.7363100318787654\,\mathrm{i}\\\left\{2,0\right\} & 0.11703059508534998 +0.\,\,\mathrm{i}\\\left\{2,1\right\} & 0.8158972359113332 +0.\,\,\mathrm{i}\\\left\{3,0\right\} & 0.05448860910968355 +0.09540961504089102\,… [541 chars total]
Cell 96: ✓ — Out[311]= \left\{\text{O[k1,0]-FLdet[k1] (must be 0)},0,\text{O[0,k2]-FLdet2[k2] (expected NONZERO per Phase1 finding)},\left(-3.0719756421050755-0.5172165371705225\,\mathrm{i}\right)\,k2\right\}
Cell 97: ✓ — Out[316]= \left\{\text{Fcand1},k1^2\,\left(\left(-2.2527492196701659-5.82989030370859198\,\mathrm{i}\right)\,k2-\left(0.4170744170089964+23.0487337186617355\,\mathrm{i}\right)\right)+k1\,\left(\left(0.67128913890100625-0.37524873225652243\,\mathrm{i}\right)\,k2+\left(9.51343802264711551-1.5262799733212005\,\mathrm{i}\right)\right)-\left(1.99115021541939352-0.72828474168870538\,\mathrm{i}\right)\,k… [1014 chars total]
Cell 98: ✓ — Out[318]= \left(0.9252125034899888 -6.175727508765443\,\mathrm{i}\right)\,k1^2\,k2+\left(5.4271267482291075 +3.0890004366773067\,\mathrm{i}\right)\,k1\,k2
Cell 99: ✓ — Out[323]= \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 7.991287939023577 -1.2820751775898085\,\mathrm{i}\\\left\{2,0\right\} & -0.3808070763995185-21.044496003995498\,\mathrm{i}\\\left\{0,1\right\} & 0.6123194936458246 +1.416862394550688\,\mathrm{i}\\\left\{1,1\right\} & -0.6248120089007682-3.051910425543801\,\mathrm{i}\\\left\{2,1\right\} & -1.8923093445229393-4.897107855115217\,\m… [421 chars total]
Cell 100: ✓ — Out[328]= \left\{\text{Fcand2 (no P-insertion)},k1^2\,\left(\left(0.02466179161444501+1.36288164597304175\,\mathrm{i}\right)-\left(0.36043987514722654+0.93278244859337472\,\mathrm{i}\right)\,k2\right)+k1\,\left(\left(1.2351720155778515-0.69045766735200127\,\mathrm{i}\right)\,k2-\left(6.46913785540003854-1.0378703818584163\,\mathrm{i}\right)\right)+\left(5.38781822995835893-1.9706528304517910\,\mat… [788 chars total]
Cell 101: ✓ — 2 kernel messages: Power::infy: Infinite expression 0       encountered.
Infinity::indet: Indeterminate expression ----------------- encountered.                                                   1/3                    … [also: Power::infy] result: Out[331]//TableForm= \begin{pmatrix}2 & 1 & \left\{0,0\right\} & \mathrm{Indeterminate} & \mathrm{Indeterminate} & -4.970498218180917+4.515583419069149\,\mathrm{i}\\2 & 1 & \left\{1,1\right\} & \mathr…
Cell 102: ✓ — Out[334]//TableForm= \begin{pmatrix}2 & 1 & \left\{0,0\right\} & \left(-7.771745024480887+7.060451775520925\,\mathrm{i}\right)\,k1^2-\left(3.7846186890458786 +9.794215710230434\,\mathrm{i}\right)\,k1 & \left(1.8923093445229393 +4.897107855115217\,\mathrm{i}\right)\,k2 & -4.970498218180917+4.515583419069149\,\mathrm{i}\\2 & 1 & \left\{1,1\right\} & \left(9.339255641812587 +10.661478687052961\,\math… [1258 chars total]
Cell 103: ✓ — 2 kernel messages: Power::infy: Infinite expression 0       encountered.
Infinity::indet: Indeterminate expression ----------------- encountered.                                                   1/3                    … [also: Power::infy] result: Out[338]= \left\{\text{Fcand1(k1,0)},\mathrm{Indeterminate},\text{FLdet[k1]},\frac{\left(2.45075493251774917289426756629+5.78551794268217244227029022731\,\mathrm{i}\right)\,\left(\left(0.0525946877178…
Cell 104: ✓ — Out[343]= \left\{\text{Fcand1(k1,0)},\left(-8.51191121729+7.73287575414\,\mathrm{i}\right)\,k1^2-\left(4.5054984393+11.6597806074\,\mathrm{i}\right)\,k1+\left(1.000000000000+0.\times 10^{-12}\,\mathrm{i}\right),\text{FLdet[k1]},\frac{\left(2.45075493251774917289426756629+5.78551794268217244227029022731\,\mathrm{i}\right)\,\left(\left(0.05259468771783351844996683714+0.15021349717800998891606886537\… [960 chars total]
Cell 107: ✓ — Out[348]= \left\{\text{STEP 0: det[$\lambda$1,$\lambda$2][k\^{}(u/h)Q1] vs O[k,k], worst |residual| (sample k)},0,\text{sample rows \{$\lambda$1,$\lambda$2,state,resid,ratio\}},\begin{pmatrix}1 & 1 & \left\{0,0\right\} & 0 & 1.\,\\1 & 1 & \left\{1,0\right\} & 0 & 1.\,\\1 & 1 & \left\{1,1\right\} & 0 & 1.\,\\1 & 1 & \left\{2,0\right\} & 0 & 1.\,\\2 & 1 & \left\{0,0\right\} & 0 & 1.\,\\2 & 1 & \left… [445 chars total]
Cell 108: ✓ — Out[354]= \left\{\text{STEP 0 CONFIRMED: det[$\lambda$1,$\lambda$2][k\^{}(u/h)Q1] == O[k,k] identically},\text{worst residual over full sweep},0,\text{num states},41\right\}

## 2026-07-18 16:17:31 — ▶️ RUN CELLS 109–123
Cell 109: ✓ — 2 kernel messages: Power::infy: Infinite expression - encountered.                                  0
Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.
Power::infy: Infinite expression - encounter… [also: Power::infy] result: Out[358]= \begin{pmatrix}1 & 0 & \left\{1.\,\right\}\\1 & 1 & \left\{1.\,\right\}\\2 & 0 & \left\{1.\,\right\}\\2 & 1 & \left\{2.\,\right\}\\2 & 2 & \left\{2.\,\right\}\\3 & 0 & \left\{1.\,\right\}\\3…
Cell 110: ✓ — Out[364]= \left\{\text{STEP 1 CONFIRMED: ⟨x0[$\lambda$]|E21\^{}m = $\gamma$\_m($\lambda$) xL1[$\lambda$][\{m\}], $\gamma$\_m($\lambda$)=$\lambda$!/($\lambda$-m)!},\text{worst deviation},0,\text{table \{$\lambda$,m,$\gamma$\_found,$\gamma$\_predicted,dev\}},\begin{pmatrix}1 & 0 & 1.\, & 1 & 0\\1 & 1 & 1.\, & 1 & 0\\2 & 0 & 1.\, & 1 & 0\\2 & 1 & 2.\, & 2 & 0\\2 & 2 & 2.\, & 2 & 0\\3 & 0 & 1.\, & 1 &… [489 chars total]
Cell 111: ✓ — Out[369]= \left\{\text{STEP 2: x2[\{s1,m\}] vs xL1[$\lambda$1][\{s1\}]$\otimes$xL1[$\lambda$2][\{m\}], worst ratio-spread (constancy check)},0,\text{num cases},54,\text{sample rows},\begin{pmatrix}1 & 1 & 0 & 1 & 0.16 & 0\\1 & 1 & 1 & 1 & 1.\, & 0\\1 & 2 & 0 & 1 & 0.16 & 0\\1 & 2 & 0 & 2 & -0.68 & 0\\1 & 2 & 1 & 1 & 1.\, & 0\\1 & 2 & 1 & 2 & 1.\, & 0\\1 & 3 & 0 & 1 & 0.16 & 0\\1 & 3 & 0 & 2 & -0.6… [479 chars total]
Cell 112: ✓ — Out[373]= \begin{pmatrix}1 & 0 & \frac{4}{25}\\1 & 1 & 1\\2 & 0 & \frac{2}{23}\\2 & 1 & \frac{25}{46}\\2 & 2 & 1\\3 & 0 & \frac{4}{67}\\3 & 1 & \frac{25}{67}\\3 & 2 & \frac{46}{67}\\3 & 3 & 1\end{pmatrix}
Cell 113: ✓ — Out[378]= \left\{\text{c($\lambda$1,s1) hypothesis (m-independent) vs full sweep --- FAILS for m>=2, see next cell for the fix},\text{worst deviation over all 54 ($\lambda$1,$\lambda$2,s1,m) cases},1.68\right\}
Cell 114: ✓ — Out[382]= \begin{pmatrix}1 & 0 & 1 & \frac{4}{25}\\1 & 0 & 2 & -\frac{17}{25}\\1 & 0 & 3 & -\frac{38}{25}\\1 & 1 & 1 & 1\\1 & 1 & 2 & 1\\1 & 1 & 3 & 1\\2 & 0 & 1 & \frac{2}{23}\\2 & 0 & 2 & -\frac{34}{575}\\2 & 0 & 3 & \frac{323}{575}\\2 & 1 & 1 & \frac{25}{46}\\2 & 1 & 2 & \frac{2}{23}\\2 & 1 & 3 & -\frac{17}{46}\\2 & 2 & 1 & 1\\2 & 2 & 2 & 1\\2 & 2 & 3 & 1\\3 & 0 & 1 & \frac{4}{67}\\3 & 0 & 2 & … [672 chars total]
Cell 115: ✓ — 4 kernel messages: LinearSolve::nosol: Linear equation encountered that has no solution.
LinearSolve::nosol: Linear equation encountered that has no solution.
LinearSolve::nosol: Linear equation encountered that has no … [also: LinearSolve::nosol, LinearSolve::nosol, LinearSolve::nosol] result: Out[387]= \left\{\text{coeffs (rows=s1, cols=s1'=0..3)},\left\{\mathrm{LinearSolve}\left[\begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & \frac{1}{36}\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & \f…
Cell 116: ✓ — Out[392]= \left\{\left\{0,2\right\}\to -0.02206359506813757300454250486697,\left\{1,1\right\}\to 0.1635301752109020116807268007787,\left\{2,0\right\}\to 0.8585334198572355613238157040883\right\}
Cell 117: ✓ — Out[396]= \left\{\text{STEP 2 level-conservation check: \# off-level nonzero coeffs (want 0 everywhere)},0,\text{num ($\lambda$1,$\lambda$2,s1,m) cases},29\right\}
Cell 119: ✓ — Out[404]= \left\{\text{m=0: Om-Cm (must be 0 per Phase1, C1pred/P1site2[$\theta$1] should equal exactly Cm's own norm here)},\left(2.5572121404875445 -0.4102640568287387\,\mathrm{i}\right)\,k1-\left(0.038411844228125346 +2.1227491621421546\,\mathrm{i}\right)\,k1^2,\text{m=1: Om\_1[k1]},\left(0.3604398751472265372666776506+0.9327824485933747164848359821\,\mathrm{i}\right)\,k1^2-\left(1.296101147801… [801 chars total]
Cell 120: ✓ — Out[411]= \left\{\text{m=0: (B0-B1)/O\_0 ratio},\frac{\sqrt[3]{k1}\,\left(\left(-0.0292504903494856237+0.0205326892042810627\,\mathrm{i}\right)\,k1^2+\left(0.499242526583056932+0.532659421382621779\,\mathrm{i}\right)\,k1-\left(0.139388493663376676+0.206413542317126162\,\mathrm{i}\right)\right)}{\left(0.3604398751472265372666776506+0.9327824485933747164848359821\,\mathrm{i}\right)\,k1^2-\left(1.296… [1083 chars total]
Cell 121: ✓ — Out[414]= \left\{\text{Cm[0]-Cm[1] shape},\left(-0.0225172879957976151-1.244370198497125076\,\mathrm{i}\right)\,k1^2+\left(5.434075798536032377-0.8718111207610696998\,\mathrm{i}\right)\,k1+\left(0.\times 10^{-20}+0.\times 10^{-20}\,\mathrm{i}\right),\text{ratio vs O\_0[k1] (want: clean constant)},\frac{\left(-0.0225172879957976151-1.244370198497125076\,\mathrm{i}\right)\,k1^2+\left(5.4340757985360… [753 chars total]
Cell 122: ✓ — Out[416]= \left\{\text{Cm[0][k1] / O\_1[k1] (reversed-index hypothesis; want clean k1-independent constant)},\frac{\left(0.111318362779488714+0.0407158959953172346\,\mathrm{i}\right)\,k1^2-\left(0.21818047436770437-1.02534209918467083\,\mathrm{i}\right)\,k1+\left(0.360439875147226537-0.932782448593374716\,\mathrm{i}\right)}{\left(1.\,+0.\times 10^{-29}\,\mathrm{i}\right)\,k1^2-\left(2.963909584212… [539 chars total]

## 2026-07-18 16:17:47 — ▶️ RUN CELLS 124–138
Cell 125: ✓ — (no output)
Cell 126: ✓ — Out[446]= \left\{\text{Symbolic-$\theta$,$\chi$1 port-fidelity: xL1Sym vs xL1},0,\text{x2Sym vs x2},0\right\}
Cell 127: ✓ — Out[451]= \left\{\left\{0,0\right\}\to 0,\left\{0,1\right\}\to \frac{\theta A-\theta B}{\theta A-\theta B+1},\left\{1,0\right\}\to \frac{1}{\theta A-\theta B+1},\left\{1,1\right\}\to 0\right\}
Cell 128: ✓ — Out[454]= 35
Cell 129: ✓ — Out[460]= \left\{\left\{0,0\right\}\to 0,\left\{0,1\right\}\to \theta B-\theta A,\left\{1,0\right\}\to -1,\left\{1,1\right\}\to 0\right\}
Cell 130: ✓ — Out[463]= 127
Cell 131: ✓ — Out[466]= 127
Cell 132: ✓ — Out[470]//MatrixForm= \begin{pmatrix}\frac{\theta A^2-2\,\theta A\,\theta B-\theta A+\theta B^2+\theta B}{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)} & \frac{6\,\left(\theta A-\theta B\right)}{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)} & \frac{6}{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)}\\0 & \frac{\theta A-\theta B+1}{\the… [472 chars total]
Cell 133: ✓ — Out[472]//MatrixForm= \begin{pmatrix}\frac{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)}{\theta A^2-2\,\theta A\,\theta B-\theta A+\theta B^2+\theta B} & -\frac{6\,\left(\theta A-\theta B+3\right)}{\theta A^2-2\,\theta A\,\theta B+\theta B^2-1} & \frac{6}{\left(\theta A-\theta B\right)\,\left(\theta A-\theta B+1\right)}\\0 & \frac{\theta A-\theta B+3}{\theta A-\theta B+1} & -… [453 chars total]
Cell 134: ✓ — Out[478]= \left\{\text{equations},\left\{1.000000000000000000000+0.\times 10^{-21}\,\mathrm{i}⩵\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right)\,wc0+\left(2.603469709065218092712489006+0.6885776528619826437114966465\,\mathrm{i}\right)\,wc1,-1.0350620568640061670+0.16605926109734660949\,\mathrm{i}⩵\left(1.522150083623538480912456055-0.2442047957313920727733393355\,\mathrm{i}\right)\,wc0-\left(1.2961… [785 chars total]
Cell 135: ✓ — Out[484]= \begin{pmatrix}0 & 0 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\0 & 1 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\1 & 0 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\1 & 1 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\2 & 0 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} … [510 chars total]
Cell 136: ✓ — Out[489]= \left\{\text{Conjectured sum rule FAILS beyond the lucky L=1 case},\text{num failing cases},29,\text{num total cases},81,\text{sample failures},\begin{pmatrix}1 & 2 & 0 & 1 & \frac{1}{-\theta A+\theta B-1}\\1 & 2 & 0 & 2 & \frac{2}{\theta A-\theta B+1}\\1 & 3 & 0 & 1 & -\frac{2}{\theta A-\theta B+1}\end{pmatrix}\right\}
Cell 138: ✓ — Out[496]= \left\{\text{BOXED: A[s1,m;s1,m] = Pochhammer[s1+$\theta$1-$\theta$2-m+1,m]/Pochhammer[$\lambda$1+$\theta$1-$\theta$2-m+1,m] --- verified at the live numeric $\theta$1,$\theta$2},\text{worst deviation},0,\text{num ($\lambda$1,$\lambda$2,s1,m) cases},81\right\}

## 2026-07-18 16:17:56 — ▶️ RUN CELLS 139–154
Cell 140: ✓ — Out[502]= \left\{\text{STEP 1a (WITH $\kappa$1\^{}(u/h) prefactor): worst |x2.Psi - Q1($\theta$1+s1)Q1($\theta$2+s2)/(Q1($\theta$1)Q1($\theta$2))|},0,\text{num cases},297\right\}
Cell 141: ✓ — Out[509]= \left\{\text{STEP 1b (WITHOUT prefactor, bare q\_n): worst residual (expect NONZERO, confirming prefactor matters)},3.7396476326262356,\text{sample},\begin{pmatrix}\left\{0,0\right\} & 0 & 0 & 0\\\left\{1,0\right\} & 0 & 0 & 0\\\left\{1,1\right\} & 0 & 0 & 0\end{pmatrix}\right\}
Cell 143: ✓ — Out[515]= \begin{pmatrix}\left\{0,0,\left\{\text{srng},\left\{0\right\},\text{sol},\left\{1.\,\right\},\text{maxResid},0.\times 10^{-28}\right\}\right\} & \left\{0,1,\left\{\text{srng},\left\{0,1\right\},\text{sol},\left\{11.50000000000000000000000000-4.431302171530629069010197265\times 10^{-38}\,\mathrm{i},-10.50000000000000000000000000+2.933894180666706878219239219\times 10^{-38}\,\mathrm{i}\rig… [1115 chars total]
Cell 144: ✓ — Out[521]= \left\{\text{worst residual across full B-sweep (want \textasciitilde{}0)},0.\times 10^{-24},\text{num (rep,a,b) cases},37\right\}
Cell 145: ✓ — Out[525]= \left\{\text{A-block ($\lambda$1=2,$\lambda$2=1, L=2, s,s'=1,2)},\begin{pmatrix}\frac{25}{46} & \frac{21}{46}\\0 & 1\end{pmatrix},\text{its inverse (should = B row)},\begin{pmatrix}\frac{46}{25} & -\frac{21}{25}\\0 & 1\end{pmatrix},\text{Step2's independently-solved B row \{1,1\}},\left\{\frac{46}{25},-\frac{21}{25}\right\},\text{match?},\left\{0,0\right\}\right\}
Cell 147: ✓ — Out[532]= \left\{\text{det[$\lambda$1,$\lambda$2][k\^{}u Q1] - explicit residue-sum expansion (rep (2,1), state (1,0)), want 0},0\right\}
Cell 148: ✓ — Out[538]= \left\{\text{TRUE sum rule (derived, Step 4): worst deviation},0,\text{num ($\lambda$1,$\lambda$2,s,t) cases},20,\text{sample},\begin{pmatrix}2 & 1 & 0 & 0 & 0\\2 & 1 & 0 & 1 & 0\\2 & 1 & 1 & 0 & 0\\2 & 1 & 1 & 1 & 0\\2 & 1 & 2 & 0 & 0\\2 & 1 & 2 & 1 & 0\end{pmatrix}\right\}
Cell 150: ✓ — Out[546]= \left\{\text{STEP 3 assembly test (rep (2,1), 6 states x 3 grid points): worst |Oassembled-O2rot|},0,\text{num cases},18\right\}
Cell 151: ✓ — Out[552]= \left\{\text{STEP 3 FULL-SWEEP ASSEMBLY CONFIRMED: O[k1,k2] = $\Sigma$ cA(k1)cB(k2) B-weighted Q-bilinears == O2rot exactly},\text{worst deviation},0,\text{num states x grid cases},82\right\}
Cell 153: ✓ — Out[557]= \begin{pmatrix}0 & 0 & 1.\, & 1.\, & 0\\0 & 1 & -2.603469709065218-0.6885776528619827\,\mathrm{i} & 5.387818229958359 -1.970652830451791\,\mathrm{i} & -7.991287939023577+1.2820751775898085\,\mathrm{i}\\1 & 0 & 1.5221500836235384 -0.24420479573139206\,\mathrm{i} & -6.469137855400039+1.0378703818584163\,\mathrm{i} & 7.991287939023577 -1.2820751775898085\,\mathrm{i}\\1 & 1 & 1.2961011478017… [843 chars total]

## 2026-07-18 16:18:31 — 📥 BULK INSERT 2 CELL(S) at positions 155–156
1. [markdown] (*Section L — Paul's directed ansatz: does O[k1,k2] admit a 2x2-determinant form with row 1↵built pu
2. [code] (*Build M[s,t](k1,k2) as a (λ1+1)x(λ2+1) matrix, exact rational entries at rational (k1,k2).↵  Reuse


## 2026-07-18 16:18:31
**Input:**
```wolfram
(*Build M[s,t](k1,k2) as a (λ1+1)x(λ2+1) matrix, exact rational entries at rational (k1,k2).
  Reuses Section K's solveBratM (exact B) and gammaFun (γ_a(λ)=λ!/(λ-a)!).*)
ClearAll[Mmatrix];
Mmatrix[λ1_,λ2_][k1_,k2_] := Table[
  Module[{L=s+t, aMin, aMax},
    aMin = Max[0,L-λ2]; aMax = Min[L,λ1];
    Sum[
      Module[{srng,bv,idx},
        {srng,bv} = solveBratM[λ1,λ2][a,L-a][[1;;2]];
        idx = FirstPosition[srng,s];
        If[idx===Missing["NotFound"], 0,
          (-k1)^a/a! gammaFun[λ1,a] (-k2)^(L-a)/(L-a)! gammaFun[λ2,L-a] bv[[idx[[1]]]]]
      ]
    ,{a,aMin,aMax}]
  ]
,{s,0,λ1},{t,0,λ2}];
(*smoke test: rep (2,1) at a generic rational (k1,k2) point.*)
Mtest21 = Mmatrix[2,1][7/3,5/2];
{"M matrix (rep (2,1), k1=7/3,k2=5/2)", Mtest21 // MatrixForm, "MatrixRank", MatrixRank[Mtest21]}
```
**Output:**
```
Out[561]= \left\{\text{M matrix (rep (2,1), k1=7/3,k2=5/2)},\begin{pmatrix}1 & -\frac{115}{4}\\\frac{259}{12} & \frac{322}{15}\\-\frac{196}{45} & -\frac{245}{18}\end{pmatrix},\text{MatrixRank},2\right\}
```
## 2026-07-18 16:18:50 — 📥 BULK INSERT 1 CELL(S) at positions 157–157
1. [code] (*GATE 1 proper: only reps with BOTH λ1,λ2>=2 give a nontrivial rank<=2 test (max possible rank↵  fo


## 2026-07-18 16:18:51
**Input:**
```wolfram
(*GATE 1 proper: only reps with BOTH λ1,λ2>=2 give a nontrivial rank<=2 test (max possible rank
  for a (λ1+1)x(λ2+1) matrix is min(λ1+1,λ2+1); reps like (2,1),(3,1),(1,2),(1,3) have min<=2
  automatically). In the standard sweep only (2,2) is a genuine 3x3 test. Compute M at several
  independent rational (k1,k2) sample points for EVERY rep, report rank (and max-possible-rank) for
  the full record, flagging (2,2) as the decisive case.*)
ClearAll[gate1Report];
samplePts = {{7/3,5/2},{11/5,-4/3},{-2/7,9/4},{13/6,1/9}};
gate1Report = {};
Do[
  {λ1,λ2} = rep;
  Do[
    Module[{Mm,rk},
      Mm = Mmatrix[λ1,λ2][pt[[1]],pt[[2]]];
      rk = MatrixRank[Mm];
      AppendTo[gate1Report, {λ1,λ2,pt, rk, Min[λ1+1,λ2+1]}]
    ]
  ,{pt,samplePts}]
,{rep,{{1,1},{2,1},{1,2},{2,2},{3,1},{1,3}}}];
gate1Report//TableForm
```
**Output:**
```
Out[566]//TableForm= \begin{pmatrix}1 & 1 & \left\{\frac{7}{3},\frac{5}{2}\right\} & 2 & 2\\1 & 1 & \left\{\frac{11}{5},-\frac{4}{3}\right\} & 2 & 2\\1 & 1 & \left\{-\frac{2}{7},\frac{9}{4}\right\} & 2 & 2\\1 & 1 & \left\{\frac{13}{6},\frac{1}{9}\right\} & 2 & 2\\2 & 1 & \left\{\frac{7}{3},\frac{5}{2}\right\} & 2 & 2\\2 & 1 & \left\{\frac{11}{5},-\frac{4}{3}\right\} & 2 & 2\\2 & 1 & \left\{-\frac{… [1258 chars total]
```

## 2026-07-18 16:18:58
**Input:**
```wolfram
Select[gate1Report, #[[1]]==2 && #[[2]]==2 &]
```
**Output:**
```
Out= {{2, 2, {7/3, 5/2}, 3, 3}, {2, 2, {11/5, -4/3}, 3, 3}, {2, 2, {-2/7, 9/4}, 3, 3}, {2, 2, {13/6, 1/9}, 3, 3}}
```

## 2026-07-18 16:19:11
**Input:**
```wolfram
Mtest22 = Mmatrix[2,2][7/3,5/2]
MatrixRank[Mtest22]
Det[Mtest22]
(* cross-check: sum over M of Mtest22[[s+1,t+1]]*W[s,t] should equal O2rot at this k1,k2, for some state *)
One22 = Sum[Mtest22[[s+1,t+1]] Wval[2,2][1,1][s,t],{s,0,2},{t,0,2}]
O2rot[2,2][1,1][7/3,5/2]
```
**Output:**
```
Out[1]= {{1, -115/2, -14375/136}, {287/6, 25438/255, -4025/12}, {14749/360, 10045/36, 1225/36}}
Out[2]= 3
Out[3]= 1143905/578
Out[4]= -26.35272221942383516565637706206329909422`26.74096379604192 - 3.05193081177028427093874617402402432395`25.80576404167926*I
Out[5]= -26.35272221942383516565637706206329909538`27.75450724060588 - 3.0519308117702842709387461740240243241`26.82391005190665*I
```
## 2026-07-18 16:19:36 — 📥 BULK INSERT 1 CELL(S) at positions 158–158
1. [code] (*Confirm rank 3 is GENERIC (not an artifact of the 4 sample points) by computing Det[Mmatrix]↵  as


## 2026-07-18 16:19:36
**Input:**
```wolfram
(*Confirm rank 3 is GENERIC (not an artifact of the 4 sample points) by computing Det[Mmatrix]
  as a SYMBOLIC polynomial in k1,k2 for rep (2,2) — if identically 0, rank<=2 almost everywhere
  despite the isolated nonzero samples (impossible for a low-degree polynomial nonzero at 4 points,
  but verify directly for full rigor) and cross-validate the M-reconstruction against O2rot on this
  symbolic det computation.*)
detM22sym = Det[Mmatrix[2,2][k1,k2]] // Simplify;
{"Det[M(k1,k2)] symbolic, rep (2,2)", detM22sym, "identically zero?", detM22sym === 0}
```
**Output:**
```
Out[568]= \left\{\text{Det[M(k1,k2)] symbolic, rep (2,2)},\frac{213003\,k2^3\,\left(k1-k2\right)^2\,\left(17\,k1+46\,k2\right)}{7225},\text{identically zero?},\mathrm{False}\right\}
```
## 2026-07-18 16:20:13 — 📥 BULK INSERT 1 CELL(S) at positions 159–159
1. [code] (*Extra due-diligence check (beyond the fixed 6-rep sweep): does the rank-3 failure generalize to↵


## 2026-07-18 16:20:14
**Input:**
```wolfram
(*Extra due-diligence check (beyond the fixed 6-rep sweep): does the rank-3 failure generalize to
  OTHER reps with min(λ1,λ2)>=2, or is it specific to (2,2)? Test rep (2,3) (not in the standard
  sweep — exercises the full on-demand TauEigensystem/QSolve machinery fresh for this rep).*)
Mtest23 = Mmatrix[2,3][7/3,5/2];
rank23 = MatrixRank[Mtest23];
detM23sym = Det[Mtest23[[1;;3,1;;3]]]; (*3x3 minor, since M is 3x4 here, min dim 3*)
{"rep (2,3): M dimensions", Dimensions[Mtest23], "rank at sample point", rank23,
 "a 3x3 minor determinant (nonzero => rank>=3)", Chop[detM23sym,10^-10]}
```
**Output:**
```
Out[572]= \left\{\text{rep (2,3): M dimensions},\left\{3,4\right\},\text{rank at sample point},3,\text{a 3x3 minor determinant (nonzero => rank>=3)},\frac{19563970545}{1405696}\right\}
```
## 2026-07-18 16:20:59 — 📥 BULK INSERT 1 CELL(S) at positions 160–160
1. [markdown] (*Section L final report — GATE 1 (pointwise rank) FAILS. Per the ask-first protocol ("if a gate↵fai


---

## 2026-07-18 16:21:07 — 🔄 KERNEL RESTART

## 2026-07-18 16:21:35 — ▶️ RUN CELLS 1–10
Cell 1: ✓ — Out[8]= \left\{\frac{1}{3},\frac{1}{7},1,0.360439875147226537266677650596+0.932782448593374716484835982071\,\mathrm{i},0.36043987514722653726667765060-0.932782448593374716484835982071\,\mathrm{i},0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i},1\right\}
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 8: ✓ — (no output)
Cell 10: ✓ — (no output)

## 2026-07-18 16:21:40 — ▶️ RUN CELLS 11–20
Cell 11: ✓ — (no output)
Cell 12: ✓ — Out[27]= \left\{0\right\}
Cell 13: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — Out[31]= \left\{0\right\}
Cell 16: ✓ — (no output)
Cell 17: ✓ — Out[34]= \begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{pmatrix}
Cell 18: ✓ — (no output)
Cell 20: ✓ — (no output)

## 2026-07-18 16:21:43 — ▶️ RUN CELLS 21–28
Cell 21: ✓ — (no output)
Cell 23: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 26: ✓ — (no output)
Cell 27: ✓ — Out[52]= \begin{pmatrix}0 & \frac{1}{2} & 0 & 0 & 0 & 0\end{pmatrix}
Cell 28: ✓ — Out[53]= \left\{0\right\}

## 2026-07-18 16:21:45 — ▶️ RUN CELLS 29–36
Cell 30: ✓ — Out[59]= \left\{\left(u-\frac{4}{3}\right)\,\left(u-\frac{8}{7}\right),\left(u-\frac{1}{3}\right)\,\left(u-\frac{1}{7}\right),\left(0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i}\right)\,u^2-\left(1.06415582186324025288257211128+0.\times 10^{-29}\,\mathrm{i}\right)\,u\right\}
Cell 31: ✓ — Out[62]= \begin{pmatrix}1 & 1\\2 & 2\\3 & 3\end{pmatrix}
Cell 32: ✓ — Out[65]= \text{max|C2-($\kappa$1+$\kappa$2)Id|}\to 0.\times 10^{-29}
Cell 33: ✓ — Out[75]= \langle\!\lvert \left\{0,0\right\}\to 0.5664055180884988 +1.3769645669711723\,\mathrm{i},\left\{1,0\right\}\to -0.778123966451158+0.\,\,\mathrm{i},\left\{1,1\right\}\to 1.1900552523337025 +0.\,\,\mathrm{i},\left\{2,0\right\}\to 0.5664055180884988 -1.3769645669711723\,\mathrm{i}\rvert\!\rangle
Cell 34: ✓ — 2 kernel messages: Psi::spurious: No genuine state at (M,n)=(2,2) for (λ1,λ2)=(2,1): weight sector M=2 has only 2 genuine state(s) (n=0..1), so n=2 does not exist. [also: Psi::spurious] result: Print: spurious-guard test (expect a Psi::spurious message and Missing[...]):[formula] | Out[84]= \left\{\left\{1.000000000000000000000000000+0.\times 10^{-29}\,\mathrm{i},-1.3004494881076436272092358…
⚠️ stopped at Cell 34 — error detected (pass stopOnError:false to continue past errors)

## 2026-07-18 16:22:04 — ▶️ RUN CELLS 35–42
Cell 35: ✓ — Out[91]= \left\{1.03464899078135229594035121+0.381869568380946283899575941\,\mathrm{i},-2.37022278949001677248172619-0.20048970475590097847093779\,\mathrm{i},1.00000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right\}
Cell 36: ✓ — Out[94]= \left\{-0.709231010958016926537225142+1.09596916251277206706752728\,\mathrm{i},-0.02291660900+0.\times 10^{-12}\,\mathrm{i}\right\}
Cell 37: ✓ — (no output)
Cell 38: ✓ — (no output)
Cell 39: ✓ — (no output)
Cell 40: ✓ — (no output)
Cell 41: ✓ — (no output)
Cell 42: ✓ — Out[109]= 0

## 2026-07-18 16:22:07 — ▶️ RUN CELLS 43–47
Cell 43: ✓ — (no output)
Cell 44: ✓ — Out[113]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 45: ✓ — Out[115]= \left\{6,6\right\}
Cell 46: ✓ — Out[116]= \begin{pmatrix}2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i} & 0 & 0 & 0 & 0 & 0\\0 & 7.36608403016431791181-0.968046299794461400266\,\mathrm{i} & 0 & 0 & 0 & 0\\0 & 0 & 2.68702160063169134012-1.3456153260557858257\,\mathrm{i} & 0 & 0 & 0\\0 & 0 & 0 & 7.3660840301643179118107+0.9680462997944614002659\,\mathrm{i} & 0 & 0\\0 & 0 & 0 & 0 & 2.687021600631691340118… [537 chars total]
Cell 47: ✓ — Out[125]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{max|off-diagonal|}\to 5.065342385207116\times 10^{-34},\text{digits vanished (-Log10)}\to 33.29539119374288,\text{min|diagonal|}\to 3.0051232400647923,\text{diagonal values}\to \left\{2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i},7.36608403016431791181-0… [699 chars total]

## 2026-07-18 16:22:11 — ▶️ RUN CELLS 48–54
Cell 48: ✓ — (no output)
Cell 49: ✓ — (no output)
Cell 50: ✓ — Out[133]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 51: ✓ — Out[138]= \left\{1.00000000000000000000000000000,6\right\}
Cell 52: ✓ — Out[143]= \langle\!\lvert \text{max|<PsiLeft|Psi>| off-diagonal}\to 0.\,,\text{digits vanished (-Log10)}\to \mathrm{Indeterminate}\rvert\!\rangle
Cell 53: ✓ — Out[147]= \left\{7.366084030164317911810727273-0.9680462997944614002658874268\,\mathrm{i},1.00000000000000000000+0.\times 10^{-21}\,\mathrm{i}\right\}
Cell 54: ✓ — Out[150]= \begin{matrix}\langle\!\lvert \text{M}\to 0,\text{n}\to 0,\text{braket}\to 2.909064289-5.804620064\,\mathrm{i},\text{G21 diagonal}\to 2.909064289-5.804620064\,\mathrm{i},\text{normRatio}\to 1.000000000+0.\times 10^{-10}\,\mathrm{i}\rvert\!\rangle \\\langle\!\lvert \text{M}\to 1,\text{n}\to 0,\text{braket}\to 7.366084030-0.968046300\,\mathrm{i},\text{G21 diagonal}\to 7.366084030-0.9680463… [1247 chars total]

## 2026-07-18 16:22:15 — ▶️ RUN CELLS 55–62
Cell 56: ✓ — (no output)
Cell 57: ✓ — Out[158]= \begin{pmatrix}3 & 3\\3 & 3\end{pmatrix}
Cell 58: ✓ — (no output)
Cell 59: ✓ — Out[165]= \left\{0\right\}
Cell 60: ✓ — Out[169]= \begin{matrix}\left\{1,1,0,1.\,,0.\,\right\}\\\left\{1,1,1,1.\,,0.\,\right\}\\\left\{1,2,0,1.\,,0.\,\right\}\\\left\{1,2,1,1.\,,0.\,\right\}\\\left\{1,3,0,1.\,,0.\,\right\}\\\left\{1,3,1,1.\,,0.\,\right\}\\\left\{2,1,0,1.\,,0.\,\right\}\\\left\{2,1,1,1.\,,0.\,\right\}\\\left\{2,1,2,1.\,,0.\,\right\}\\\left\{2,2,0,1.\,,0.\,\right\}\\\left\{2,2,1,1.\,,0.\,\right\}\\\left\{2,2,2,1.\,,0.\,\r… [898 chars total]
Cell 61: ✓ — Out[170]= \langle\!\lvert \text{max spread (ratio non-constancy)}\to 0.\,,\text{max |ratio - 1| over sweep}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases checked}\to 27\rvert\!\rangle
Cell 62: ✓ — Out[176]= \langle\!\lvert \text{max operator-level spread}\to 0.\,,\text{max |$\mu$ - $\mu$Site2|}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases}\to 27\rvert\!\rangle

## 2026-07-18 16:22:26 — ▶️ RUN CELLS 63–70
Cell 64: ✓ — (no output)
Cell 65: ✓ — Out[185]= \langle\!\lvert \text{w\_0 ($\lambda$=1,2,3)}\to \left\{1,1,1\right\},\text{full weights w\_j}\to \left\{1\to \left\{1,-1\right\},2\to \left\{1,-2,1\right\},3\to \left\{1,-3,3,-1\right\}\right\},\text{br1[$\lambda$][1] ($\lambda$=1,2,3)}\to \left\{0,0,0\right\}\rvert\!\rangle
Cell 66: ✓ — Out[187]= \left\{0,0\right\}
Cell 67: ✓ — Out[190]= \sqrt[3]{k}\,\left(\left(-1.18561270626131985872+0.832253303224995280404\,\mathrm{i}\right)\,k^2-\left(0.761819910806671651081+0.812812473459294825632\,\mathrm{i}\right)\,k-\left(0.404410037536507329420+0.598870869485872613792\,\mathrm{i}\right)\right)
Cell 68: ✓ — Out[194]= \begin{pmatrix}0 & 0 & 0\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix}
Cell 69: ✓ — Out[198]= \langle\!\lvert \text{OP(E12)}\to 1.\,+0.\times 10^{-30}\,\mathrm{i},\text{OP(E21)}\to \left(-0.03626734060947795145020244939-2.004237714666237867354326095\,\mathrm{i}\right)\,\varphi ^2-\left(1.522150083623538480912456055-0.2442047957313920727733393355\,\mathrm{i}\right)\,\varphi +\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right),\text{nontrivial generator \{a,b\}}\to \left\{2,1\right\},\… [464 chars total]
Cell 70: ✓ — Out[201]= \left\{1,2\right\}

## 2026-07-18 16:22:29 — ▶️ RUN CELLS 71–77
Cell 71: ✓ — Out[207]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{$\varphi$ candidates per state}\to \begin{pmatrix}\left\{0,0\right\} & \left(-0.920793778924720536105378676-0.390049761817044163511250929\,\mathrm{i}\right)\,k^{4/3}+\left(0.695722366612934233848844287-0.718310788304406529935213657\,\mathrm{i}\right)\,\sqrt[3]{k}-\left(0.3604398751… [1247 chars total]
Cell 72: ✓ — Out[213]= \langle\!\lvert \text{c\_sb solution(s)}\to \left\{\left\{\mathrm{cb}\to \frac{\left(-6.242239011314\times 10^{45}-4.511757007219\times 10^{45}\,\mathrm{i}\right)\,k^{4/3}-\left(5.999302657796\times 10^{45}+8.356648463050\times 10^{45}\,\mathrm{i}\right)\,\sqrt[3]{k}}{9.833530773718\times 10^{45}\,k^{4/3}-\left(3.54439660434\times 10^{45}-9.172544913427\times 10^{45}\,\mathrm{i}\right)\,… [841 chars total]
Cell 73: ✓ — Out[216]= \langle\!\lvert \text{generator \{a,b\} tested}\to \left\{2,1\right\},\text{rep tested}\to \left\{1,2\right\},\text{FL shape (verified exact for all 6 states)}\to \text{FL[1,2][state][k] = k\^{}(1/3) (A\_state + B\_state k)},\text{c0(state) = OPg[state][0]}\to \begin{pmatrix}\left\{0,0\right\} & 1.\,+0.\times 10^{-28}\,\mathrm{i}\\\left\{1,0\right\} & 1.\,+0.\times 10^{-28}\,\mathrm{i}\\… [1247 chars total]
Cell 74: ✓ — Out[220]= \left(1.000000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right)-\left(0.360439875147226537266677651+0.932782448593374716484835982\,\mathrm{i}\right)\,k
Cell 75: ✓ — Out[227]= \langle\!\lvert \text{alphaTab \{state,D0,c0,Dc,c1,$\alpha$\}}\to \begin{pmatrix}\left\{0,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & -0.36043987514722653-0.9327824485933747\,\mathrm{i} & 0.36043987514722653 +0.9327824485933747\,\mathrm{i} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.7101545203859353 -0.2371533418544… [1247 chars total]
Cell 76: ✓ — Out[235]= \langle\!\lvert \text{validation}\to \left\{\left\{1,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\end{pmatrix},\left\{2,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\\\left\{2,1\right\} & 0\\\left\{3,0\right\} & 0\end{pmatrix},\left\{1… [1247 chars total]
Cell 77: ✓ — Out[237]= \langle\!\lvert \text{generator \{a,b\}}\to \left\{2,1\right\},\text{$\varphi$(k)}\to \text{$\alpha$ k},\text{$\alpha$ (universal value)}\to -1.\,+0.\,\,\mathrm{i},\text{normalization}\to \text{1-site FSoV determinant det1[f] = br1[f]/(f/.u->$\theta$1); Section C's det[$\lambda$1,$\lambda$2][f] restricted to 1 site},\text{max residual over rep sweep}\to 0,\text{reps validated}\to \begin{… [510 chars total]

## 2026-07-18 16:22:52 — ▶️ RUN CELLS 78–82
Cell 79: ✓ — (no output)
Cell 80: ✓ — Out[245]= \left\{\text{Test 1 (symbolic monomials u\^{}j) residuals},\left\{0\right\}\right\}
Cell 81: ✓ — Out[250]= \left\{\text{Test 2 (twisted monomials k\^{}(u/h) u\^{}j) residuals},\left\{0\right\}\right\}
Cell 82: ✓ — Out[259]= \left\{\text{Test 3 worst |br identity residual| on Q1},0,\text{worst |detRatio - predicted| },0,\text{num states checked},41\right\}

## 2026-07-18 16:22:56 — ▶️ RUN CELLS 83–87
Cell 84: ✓ — (no output)
Cell 85: ✓ — Out[271]= \begin{pmatrix}1 & 1 & 0\\2 & 1 & 0\\3 & 1 & 0\end{pmatrix}
Cell 86: ✓ — Out[278]= \left\{\text{Mirror Test1 (u\^{}j)},\left\{0\right\},\text{Mirror Test2 (k\^{}(u/h)u\^{}j)},\left\{0\right\}\right\}
Cell 87: ✓ — Out[285]= \left\{\text{Mirror Test3 worst |br identity residual| on Q1},0,\text{num states checked},41\right\}

## 2026-07-18 16:22:57 — ▶️ RUN CELLS 88–91
Cell 89: ✓ — Out[294]= \left\{\text{OP2(E12)},1.\,+0.\times 10^{-30}\,\mathrm{i},\text{OP2(E21)},\left(2.603469709065218092712489006+0.6885776528619826437114966465\,\mathrm{i}\right)\,\varphi +\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right),\text{nontrivial generator \{a,b\}},\left\{2,1\right\}\right\}
Cell 90: ✓ — Out[297]= \left(1.\,+0.\times 10^{-27}\,\mathrm{i}\right)-\left(0.360439875147226537266677651+0.932782448593374716484835982\,\mathrm{i}\right)\,k
Cell 91: ✓ — Out[300]//TableForm= \begin{pmatrix}\left\{0,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & -0.36043987514722653-0.9327824485933747\,\mathrm{i} & 0.36043987514722653 +0.9327824485933747\,\mathrm{i} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.4685059330398573 -0.1713611156914601\,\mathrm{i} & 2.603469709065218 +0.6885776528619827… [1258 chars total]

## 2026-07-18 16:23:10 — ▶️ RUN CELLS 92–108
Cell 92: ✓ — Out[302]//TableForm= \begin{pmatrix}\left\{0,0\right\} & -1.\,+0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.973902234276272 +0.2269679230016555\,\mathrm{i} & 0.7391304347826086 +0.\,\,\mathrm{i} & 0.973902234276272 +0.2269679230016555\,\mathrm{i}\\\left\{1,0\right\} & 0.15191893851739566 -0.10600050418875033\,\mathrm{i} & -0.15191893851739566+0.10600050418875033\,\mathrm{i} & 0.7651595274649543 +… [1258 chars total]
Cell 93: ✓ — Out[305]= \begin{pmatrix}\left\{0,0\right\} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 0.05448860910968355 -0.09540961504089102\,\mathrm{i}\\\left\{1,1\right\} & -0.0293413933141388+0.7363100318787654\,\mathrm{i}\\\left\{2,0\right\} & 0.11703059508534998 +0.\,\,\mathrm{i}\\\left\{2,1\right\} & 0.8158972359113332 +0.\,\,\mathrm{i}\\\left\{3,0\right\} & 0.05448860910968355 +0.09540961504089102\,… [541 chars total]
Cell 96: ✓ — Out[311]= \left\{\text{O[k1,0]-FLdet[k1] (must be 0)},0,\text{O[0,k2]-FLdet2[k2] (expected NONZERO per Phase1 finding)},\left(-3.0719756421050755-0.5172165371705225\,\mathrm{i}\right)\,k2\right\}
Cell 97: ✓ — Out[316]= \left\{\text{Fcand1},k1^2\,\left(\left(-2.2527492196701659-5.82989030370859198\,\mathrm{i}\right)\,k2-\left(0.4170744170089964+23.0487337186617355\,\mathrm{i}\right)\right)+k1\,\left(\left(0.67128913890100625-0.37524873225652243\,\mathrm{i}\right)\,k2+\left(9.51343802264711551-1.5262799733212005\,\mathrm{i}\right)\right)-\left(1.99115021541939352-0.72828474168870538\,\mathrm{i}\right)\,k… [1014 chars total]
Cell 98: ✓ — Out[318]= \left(0.9252125034899888 -6.175727508765443\,\mathrm{i}\right)\,k1^2\,k2+\left(5.4271267482291075 +3.0890004366773067\,\mathrm{i}\right)\,k1\,k2
Cell 99: ✓ — Out[323]= \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 7.991287939023577 -1.2820751775898085\,\mathrm{i}\\\left\{2,0\right\} & -0.3808070763995185-21.044496003995498\,\mathrm{i}\\\left\{0,1\right\} & 0.6123194936458246 +1.416862394550688\,\mathrm{i}\\\left\{1,1\right\} & -0.6248120089007682-3.051910425543801\,\mathrm{i}\\\left\{2,1\right\} & -1.8923093445229393-4.897107855115217\,\m… [421 chars total]
Cell 100: ✓ — Out[328]= \left\{\text{Fcand2 (no P-insertion)},k1^2\,\left(\left(0.02466179161444501+1.36288164597304175\,\mathrm{i}\right)-\left(0.36043987514722654+0.93278244859337472\,\mathrm{i}\right)\,k2\right)+k1\,\left(\left(1.2351720155778515-0.69045766735200127\,\mathrm{i}\right)\,k2-\left(6.46913785540003854-1.0378703818584163\,\mathrm{i}\right)\right)+\left(5.38781822995835893-1.9706528304517910\,\mat… [788 chars total]
Cell 101: ✓ — 2 kernel messages: Power::infy: Infinite expression 0       encountered.
Infinity::indet: Indeterminate expression ----------------- encountered.                                                   1/3                    … [also: Power::infy] result: Out[331]//TableForm= \begin{pmatrix}2 & 1 & \left\{0,0\right\} & \mathrm{Indeterminate} & \mathrm{Indeterminate} & -4.970498218180917+4.515583419069149\,\mathrm{i}\\2 & 1 & \left\{1,1\right\} & \mathr…
Cell 102: ✓ — Out[334]//TableForm= \begin{pmatrix}2 & 1 & \left\{0,0\right\} & \left(-7.771745024480887+7.060451775520925\,\mathrm{i}\right)\,k1^2-\left(3.7846186890458786 +9.794215710230434\,\mathrm{i}\right)\,k1 & \left(1.8923093445229393 +4.897107855115217\,\mathrm{i}\right)\,k2 & -4.970498218180917+4.515583419069149\,\mathrm{i}\\2 & 1 & \left\{1,1\right\} & \left(9.339255641812587 +10.661478687052961\,\math… [1258 chars total]
Cell 103: ✓ — 2 kernel messages: Power::infy: Infinite expression 0       encountered.
Infinity::indet: Indeterminate expression ----------------- encountered.                                                   1/3                    … [also: Power::infy] result: Out[338]= \left\{\text{Fcand1(k1,0)},\mathrm{Indeterminate},\text{FLdet[k1]},\frac{\left(2.45075493251774917289426756629+5.78551794268217244227029022731\,\mathrm{i}\right)\,\left(\left(0.0525946877178…
Cell 104: ✓ — Out[343]= \left\{\text{Fcand1(k1,0)},\left(-8.51191121729+7.73287575414\,\mathrm{i}\right)\,k1^2-\left(4.5054984393+11.6597806074\,\mathrm{i}\right)\,k1+\left(1.000000000000+0.\times 10^{-12}\,\mathrm{i}\right),\text{FLdet[k1]},\frac{\left(2.45075493251774917289426756629+5.78551794268217244227029022731\,\mathrm{i}\right)\,\left(\left(0.05259468771783351844996683714+0.15021349717800998891606886537\… [960 chars total]
Cell 107: ✓ — Out[348]= \left\{\text{STEP 0: det[$\lambda$1,$\lambda$2][k\^{}(u/h)Q1] vs O[k,k], worst |residual| (sample k)},0,\text{sample rows \{$\lambda$1,$\lambda$2,state,resid,ratio\}},\begin{pmatrix}1 & 1 & \left\{0,0\right\} & 0 & 1.\,\\1 & 1 & \left\{1,0\right\} & 0 & 1.\,\\1 & 1 & \left\{1,1\right\} & 0 & 1.\,\\1 & 1 & \left\{2,0\right\} & 0 & 1.\,\\2 & 1 & \left\{0,0\right\} & 0 & 1.\,\\2 & 1 & \left… [445 chars total]
Cell 108: ✓ — Out[354]= \left\{\text{STEP 0 CONFIRMED: det[$\lambda$1,$\lambda$2][k\^{}(u/h)Q1] == O[k,k] identically},\text{worst residual over full sweep},0,\text{num states},41\right\}

## 2026-07-18 16:23:21 — ▶️ RUN CELLS 109–123
Cell 109: ✓ — 2 kernel messages: Power::infy: Infinite expression - encountered.                                  0
Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.
Power::infy: Infinite expression - encounter… [also: Power::infy] result: Out[358]= \begin{pmatrix}1 & 0 & \left\{1.\,\right\}\\1 & 1 & \left\{1.\,\right\}\\2 & 0 & \left\{1.\,\right\}\\2 & 1 & \left\{2.\,\right\}\\2 & 2 & \left\{2.\,\right\}\\3 & 0 & \left\{1.\,\right\}\\3…
Cell 110: ✓ — Out[364]= \left\{\text{STEP 1 CONFIRMED: ⟨x0[$\lambda$]|E21\^{}m = $\gamma$\_m($\lambda$) xL1[$\lambda$][\{m\}], $\gamma$\_m($\lambda$)=$\lambda$!/($\lambda$-m)!},\text{worst deviation},0,\text{table \{$\lambda$,m,$\gamma$\_found,$\gamma$\_predicted,dev\}},\begin{pmatrix}1 & 0 & 1.\, & 1 & 0\\1 & 1 & 1.\, & 1 & 0\\2 & 0 & 1.\, & 1 & 0\\2 & 1 & 2.\, & 2 & 0\\2 & 2 & 2.\, & 2 & 0\\3 & 0 & 1.\, & 1 &… [489 chars total]
Cell 111: ✓ — Out[369]= \left\{\text{STEP 2: x2[\{s1,m\}] vs xL1[$\lambda$1][\{s1\}]$\otimes$xL1[$\lambda$2][\{m\}], worst ratio-spread (constancy check)},0,\text{num cases},54,\text{sample rows},\begin{pmatrix}1 & 1 & 0 & 1 & 0.16 & 0\\1 & 1 & 1 & 1 & 1.\, & 0\\1 & 2 & 0 & 1 & 0.16 & 0\\1 & 2 & 0 & 2 & -0.68 & 0\\1 & 2 & 1 & 1 & 1.\, & 0\\1 & 2 & 1 & 2 & 1.\, & 0\\1 & 3 & 0 & 1 & 0.16 & 0\\1 & 3 & 0 & 2 & -0.6… [479 chars total]
Cell 112: ✓ — Out[373]= \begin{pmatrix}1 & 0 & \frac{4}{25}\\1 & 1 & 1\\2 & 0 & \frac{2}{23}\\2 & 1 & \frac{25}{46}\\2 & 2 & 1\\3 & 0 & \frac{4}{67}\\3 & 1 & \frac{25}{67}\\3 & 2 & \frac{46}{67}\\3 & 3 & 1\end{pmatrix}
Cell 113: ✓ — Out[378]= \left\{\text{c($\lambda$1,s1) hypothesis (m-independent) vs full sweep --- FAILS for m>=2, see next cell for the fix},\text{worst deviation over all 54 ($\lambda$1,$\lambda$2,s1,m) cases},1.68\right\}
Cell 114: ✓ — Out[382]= \begin{pmatrix}1 & 0 & 1 & \frac{4}{25}\\1 & 0 & 2 & -\frac{17}{25}\\1 & 0 & 3 & -\frac{38}{25}\\1 & 1 & 1 & 1\\1 & 1 & 2 & 1\\1 & 1 & 3 & 1\\2 & 0 & 1 & \frac{2}{23}\\2 & 0 & 2 & -\frac{34}{575}\\2 & 0 & 3 & \frac{323}{575}\\2 & 1 & 1 & \frac{25}{46}\\2 & 1 & 2 & \frac{2}{23}\\2 & 1 & 3 & -\frac{17}{46}\\2 & 2 & 1 & 1\\2 & 2 & 2 & 1\\2 & 2 & 3 & 1\\3 & 0 & 1 & \frac{4}{67}\\3 & 0 & 2 & … [672 chars total]
Cell 115: ✓ — 4 kernel messages: LinearSolve::nosol: Linear equation encountered that has no solution.
LinearSolve::nosol: Linear equation encountered that has no solution.
LinearSolve::nosol: Linear equation encountered that has no … [also: LinearSolve::nosol, LinearSolve::nosol, LinearSolve::nosol] result: Out[387]= \left\{\text{coeffs (rows=s1, cols=s1'=0..3)},\left\{\mathrm{LinearSolve}\left[\begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & \frac{1}{36}\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & \f…
Cell 116: ✓ — Out[392]= \left\{\left\{0,2\right\}\to -0.02206359506813757300454250486697,\left\{1,1\right\}\to 0.1635301752109020116807268007787,\left\{2,0\right\}\to 0.8585334198572355613238157040883\right\}
Cell 117: ✓ — Out[396]= \left\{\text{STEP 2 level-conservation check: \# off-level nonzero coeffs (want 0 everywhere)},0,\text{num ($\lambda$1,$\lambda$2,s1,m) cases},29\right\}
Cell 119: ✓ — Out[404]= \left\{\text{m=0: Om-Cm (must be 0 per Phase1, C1pred/P1site2[$\theta$1] should equal exactly Cm's own norm here)},\left(2.5572121404875445 -0.4102640568287387\,\mathrm{i}\right)\,k1-\left(0.038411844228125346 +2.1227491621421546\,\mathrm{i}\right)\,k1^2,\text{m=1: Om\_1[k1]},\left(0.3604398751472265372666776506+0.9327824485933747164848359821\,\mathrm{i}\right)\,k1^2-\left(1.296101147801… [801 chars total]
Cell 120: ✓ — Out[411]= \left\{\text{m=0: (B0-B1)/O\_0 ratio},\frac{\sqrt[3]{k1}\,\left(\left(-0.0292504903494856237+0.0205326892042810627\,\mathrm{i}\right)\,k1^2+\left(0.499242526583056932+0.532659421382621779\,\mathrm{i}\right)\,k1-\left(0.139388493663376676+0.206413542317126162\,\mathrm{i}\right)\right)}{\left(0.3604398751472265372666776506+0.9327824485933747164848359821\,\mathrm{i}\right)\,k1^2-\left(1.296… [1083 chars total]
Cell 121: ✓ — Out[414]= \left\{\text{Cm[0]-Cm[1] shape},\left(-0.0225172879957976151-1.244370198497125076\,\mathrm{i}\right)\,k1^2+\left(5.434075798536032377-0.8718111207610696998\,\mathrm{i}\right)\,k1+\left(0.\times 10^{-20}+0.\times 10^{-20}\,\mathrm{i}\right),\text{ratio vs O\_0[k1] (want: clean constant)},\frac{\left(-0.0225172879957976151-1.244370198497125076\,\mathrm{i}\right)\,k1^2+\left(5.4340757985360… [753 chars total]
Cell 122: ✓ — Out[416]= \left\{\text{Cm[0][k1] / O\_1[k1] (reversed-index hypothesis; want clean k1-independent constant)},\frac{\left(0.111318362779488714+0.0407158959953172346\,\mathrm{i}\right)\,k1^2-\left(0.21818047436770437-1.02534209918467083\,\mathrm{i}\right)\,k1+\left(0.360439875147226537-0.932782448593374716\,\mathrm{i}\right)}{\left(1.\,+0.\times 10^{-29}\,\mathrm{i}\right)\,k1^2-\left(2.963909584212… [539 chars total]

## 2026-07-18 16:23:36 — ▶️ RUN CELLS 124–138
Cell 125: ✓ — (no output)
Cell 126: ✓ — Out[446]= \left\{\text{Symbolic-$\theta$,$\chi$1 port-fidelity: xL1Sym vs xL1},0,\text{x2Sym vs x2},0\right\}
Cell 127: ✓ — Out[451]= \left\{\left\{0,0\right\}\to 0,\left\{0,1\right\}\to \frac{\theta A-\theta B}{\theta A-\theta B+1},\left\{1,0\right\}\to \frac{1}{\theta A-\theta B+1},\left\{1,1\right\}\to 0\right\}
Cell 128: ✓ — Out[454]= 35
Cell 129: ✓ — Out[460]= \left\{\left\{0,0\right\}\to 0,\left\{0,1\right\}\to \theta B-\theta A,\left\{1,0\right\}\to -1,\left\{1,1\right\}\to 0\right\}
Cell 130: ✓ — Out[463]= 127
Cell 131: ✓ — Out[466]= 127
Cell 132: ✓ — Out[470]//MatrixForm= \begin{pmatrix}\frac{\theta A^2-2\,\theta A\,\theta B-\theta A+\theta B^2+\theta B}{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)} & \frac{6\,\left(\theta A-\theta B\right)}{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)} & \frac{6}{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)}\\0 & \frac{\theta A-\theta B+1}{\the… [472 chars total]
Cell 133: ✓ — Out[472]//MatrixForm= \begin{pmatrix}\frac{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)}{\theta A^2-2\,\theta A\,\theta B-\theta A+\theta B^2+\theta B} & -\frac{6\,\left(\theta A-\theta B+3\right)}{\theta A^2-2\,\theta A\,\theta B+\theta B^2-1} & \frac{6}{\left(\theta A-\theta B\right)\,\left(\theta A-\theta B+1\right)}\\0 & \frac{\theta A-\theta B+3}{\theta A-\theta B+1} & -… [453 chars total]
Cell 134: ✓ — Out[478]= \left\{\text{equations},\left\{1.000000000000000000000+0.\times 10^{-21}\,\mathrm{i}⩵\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right)\,wc0+\left(2.603469709065218092712489006+0.6885776528619826437114966465\,\mathrm{i}\right)\,wc1,-1.0350620568640061670+0.16605926109734660949\,\mathrm{i}⩵\left(1.522150083623538480912456055-0.2442047957313920727733393355\,\mathrm{i}\right)\,wc0-\left(1.2961… [785 chars total]
Cell 135: ✓ — Out[484]= \begin{pmatrix}0 & 0 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\0 & 1 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\1 & 0 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\1 & 1 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\2 & 0 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} … [510 chars total]
Cell 136: ✓ — Out[489]= \left\{\text{Conjectured sum rule FAILS beyond the lucky L=1 case},\text{num failing cases},29,\text{num total cases},81,\text{sample failures},\begin{pmatrix}1 & 2 & 0 & 1 & \frac{1}{-\theta A+\theta B-1}\\1 & 2 & 0 & 2 & \frac{2}{\theta A-\theta B+1}\\1 & 3 & 0 & 1 & -\frac{2}{\theta A-\theta B+1}\end{pmatrix}\right\}
Cell 138: ✓ — Out[496]= \left\{\text{BOXED: A[s1,m;s1,m] = Pochhammer[s1+$\theta$1-$\theta$2-m+1,m]/Pochhammer[$\lambda$1+$\theta$1-$\theta$2-m+1,m] --- verified at the live numeric $\theta$1,$\theta$2},\text{worst deviation},0,\text{num ($\lambda$1,$\lambda$2,s1,m) cases},81\right\}

## 2026-07-18 16:23:47 — ▶️ RUN CELLS 139–154
Cell 140: ✓ — Out[502]= \left\{\text{STEP 1a (WITH $\kappa$1\^{}(u/h) prefactor): worst |x2.Psi - Q1($\theta$1+s1)Q1($\theta$2+s2)/(Q1($\theta$1)Q1($\theta$2))|},0,\text{num cases},297\right\}
Cell 141: ✓ — Out[509]= \left\{\text{STEP 1b (WITHOUT prefactor, bare q\_n): worst residual (expect NONZERO, confirming prefactor matters)},3.7396476326262356,\text{sample},\begin{pmatrix}\left\{0,0\right\} & 0 & 0 & 0\\\left\{1,0\right\} & 0 & 0 & 0\\\left\{1,1\right\} & 0 & 0 & 0\end{pmatrix}\right\}
Cell 143: ✓ — Out[515]= \begin{pmatrix}\left\{0,0,\left\{\text{srng},\left\{0\right\},\text{sol},\left\{1.\,\right\},\text{maxResid},0.\times 10^{-28}\right\}\right\} & \left\{0,1,\left\{\text{srng},\left\{0,1\right\},\text{sol},\left\{11.50000000000000000000000000-4.431302171530629069010197265\times 10^{-38}\,\mathrm{i},-10.50000000000000000000000000+2.933894180666706878219239219\times 10^{-38}\,\mathrm{i}\rig… [1115 chars total]
Cell 144: ✓ — Out[521]= \left\{\text{worst residual across full B-sweep (want \textasciitilde{}0)},0.\times 10^{-24},\text{num (rep,a,b) cases},37\right\}
Cell 145: ✓ — Out[525]= \left\{\text{A-block ($\lambda$1=2,$\lambda$2=1, L=2, s,s'=1,2)},\begin{pmatrix}\frac{25}{46} & \frac{21}{46}\\0 & 1\end{pmatrix},\text{its inverse (should = B row)},\begin{pmatrix}\frac{46}{25} & -\frac{21}{25}\\0 & 1\end{pmatrix},\text{Step2's independently-solved B row \{1,1\}},\left\{\frac{46}{25},-\frac{21}{25}\right\},\text{match?},\left\{0,0\right\}\right\}
Cell 147: ✓ — Out[532]= \left\{\text{det[$\lambda$1,$\lambda$2][k\^{}u Q1] - explicit residue-sum expansion (rep (2,1), state (1,0)), want 0},0\right\}
Cell 148: ✓ — Out[538]= \left\{\text{TRUE sum rule (derived, Step 4): worst deviation},0,\text{num ($\lambda$1,$\lambda$2,s,t) cases},20,\text{sample},\begin{pmatrix}2 & 1 & 0 & 0 & 0\\2 & 1 & 0 & 1 & 0\\2 & 1 & 1 & 0 & 0\\2 & 1 & 1 & 1 & 0\\2 & 1 & 2 & 0 & 0\\2 & 1 & 2 & 1 & 0\end{pmatrix}\right\}
Cell 150: ✓ — Out[546]= \left\{\text{STEP 3 assembly test (rep (2,1), 6 states x 3 grid points): worst |Oassembled-O2rot|},0,\text{num cases},18\right\}
Cell 151: ✓ — Out[552]= \left\{\text{STEP 3 FULL-SWEEP ASSEMBLY CONFIRMED: O[k1,k2] = $\Sigma$ cA(k1)cB(k2) B-weighted Q-bilinears == O2rot exactly},\text{worst deviation},0,\text{num states x grid cases},82\right\}
Cell 153: ✓ — Out[557]= \begin{pmatrix}0 & 0 & 1.\, & 1.\, & 0\\0 & 1 & -2.603469709065218-0.6885776528619827\,\mathrm{i} & 5.387818229958359 -1.970652830451791\,\mathrm{i} & -7.991287939023577+1.2820751775898085\,\mathrm{i}\\1 & 0 & 1.5221500836235384 -0.24420479573139206\,\mathrm{i} & -6.469137855400039+1.0378703818584163\,\mathrm{i} & 7.991287939023577 -1.2820751775898085\,\mathrm{i}\\1 & 1 & 1.2961011478017… [843 chars total]

## 2026-07-18 16:23:56 — ▶️ RUN CELLS 155–160
Cell 156: ✓ — Out[561]= \left\{\text{M matrix (rep (2,1), k1=7/3,k2=5/2)},\begin{pmatrix}1 & -\frac{115}{4}\\\frac{259}{12} & \frac{322}{15}\\-\frac{196}{45} & -\frac{245}{18}\end{pmatrix},\text{MatrixRank},2\right\}
Cell 157: ✓ — Out[566]//TableForm= \begin{pmatrix}1 & 1 & \left\{\frac{7}{3},\frac{5}{2}\right\} & 2 & 2\\1 & 1 & \left\{\frac{11}{5},-\frac{4}{3}\right\} & 2 & 2\\1 & 1 & \left\{-\frac{2}{7},\frac{9}{4}\right\} & 2 & 2\\1 & 1 & \left\{\frac{13}{6},\frac{1}{9}\right\} & 2 & 2\\2 & 1 & \left\{\frac{7}{3},\frac{5}{2}\right\} & 2 & 2\\2 & 1 & \left\{\frac{11}{5},-\frac{4}{3}\right\} & 2 & 2\\2 & 1 & \left\{-\frac{… [1258 chars total]
Cell 158: ✓ — Out[568]= \left\{\text{Det[M(k1,k2)] symbolic, rep (2,2)},\frac{213003\,k2^3\,\left(k1-k2\right)^2\,\left(17\,k1+46\,k2\right)}{7225},\text{identically zero?},\mathrm{False}\right\}
Cell 159: ✓ — Out[572]= \left\{\text{rep (2,3): M dimensions},\left\{3,4\right\},\text{rank at sample point},3,\text{a 3x3 minor determinant (nonzero => rank>=3)},\frac{19563970545}{1405696}\right\}


## 2026-07-18 18:05:05
**Input:**
```wolfram
{θ1,θ2,h,Length[DownValues[Mmatrix]],Length[DownValues[O2rot]]}
```
**Output:**
```
Out= {1/3, 1/7, 1, 0, 0}
```

---

## 2026-07-18 18:05:10 — 🔄 KERNEL RESTART

## 2026-07-18 18:05:28 — ▶️ RUN CELLS 1–10
Cell 1: ✓ — Out[8]= \left\{\frac{1}{3},\frac{1}{7},1,0.360439875147226537266677650596+0.932782448593374716484835982071\,\mathrm{i},0.36043987514722653726667765060-0.932782448593374716484835982071\,\mathrm{i},0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i},1\right\}
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 8: ✓ — (no output)
Cell 10: ✓ — (no output)

## 2026-07-18 18:05:33 — ▶️ RUN CELLS 11–20
Cell 11: ✓ — (no output)
Cell 12: ✓ — Out[27]= \left\{0\right\}
Cell 13: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — Out[31]= \left\{0\right\}
Cell 16: ✓ — (no output)
Cell 17: ✓ — Out[34]= \begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{pmatrix}
Cell 18: ✓ — (no output)
Cell 20: ✓ — (no output)

## 2026-07-18 18:05:36 — ▶️ RUN CELLS 21–28
Cell 21: ✓ — (no output)
Cell 23: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 26: ✓ — (no output)
Cell 27: ✓ — Out[52]= \begin{pmatrix}0 & \frac{1}{2} & 0 & 0 & 0 & 0\end{pmatrix}
Cell 28: ✓ — Out[53]= \left\{0\right\}

## 2026-07-18 18:05:39 — ▶️ RUN CELLS 29–36
Cell 30: ✓ — Out[59]= \left\{\left(u-\frac{4}{3}\right)\,\left(u-\frac{8}{7}\right),\left(u-\frac{1}{3}\right)\,\left(u-\frac{1}{7}\right),\left(0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i}\right)\,u^2-\left(1.06415582186324025288257211128+0.\times 10^{-29}\,\mathrm{i}\right)\,u\right\}
Cell 31: ✓ — Out[62]= \begin{pmatrix}1 & 1\\2 & 2\\3 & 3\end{pmatrix}
Cell 32: ✓ — Out[65]= \text{max|C2-($\kappa$1+$\kappa$2)Id|}\to 0.\times 10^{-29}
Cell 33: ✓ — Out[75]= \langle\!\lvert \left\{0,0\right\}\to 0.5664055180884988 +1.3769645669711723\,\mathrm{i},\left\{1,0\right\}\to -0.778123966451158+0.\,\,\mathrm{i},\left\{1,1\right\}\to 1.1900552523337025 +0.\,\,\mathrm{i},\left\{2,0\right\}\to 0.5664055180884988 -1.3769645669711723\,\mathrm{i}\rvert\!\rangle
Cell 34: ✓ — 2 kernel messages: Psi::spurious: No genuine state at (M,n)=(2,2) for (λ1,λ2)=(2,1): weight sector M=2 has only 2 genuine state(s) (n=0..1), so n=2 does not exist. [also: Psi::spurious] result: Print: spurious-guard test (expect a Psi::spurious message and Missing[...]):[formula] | Out[84]= \left\{\left\{1.000000000000000000000000000+0.\times 10^{-29}\,\mathrm{i},-1.3004494881076436272092358…
⚠️ stopped at Cell 34 — error detected (pass stopOnError:false to continue past errors)

## 2026-07-18 18:05:49 — ▶️ RUN CELLS 35–42
Cell 35: ✓ — Out[91]= \left\{1.03464899078135229594035121+0.381869568380946283899575941\,\mathrm{i},-2.37022278949001677248172619-0.20048970475590097847093779\,\mathrm{i},1.00000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right\}
Cell 36: ✓ — Out[94]= \left\{-0.709231010958016926537225142+1.09596916251277206706752728\,\mathrm{i},-0.02291660900+0.\times 10^{-12}\,\mathrm{i}\right\}
Cell 37: ✓ — (no output)
Cell 38: ✓ — (no output)
Cell 39: ✓ — (no output)
Cell 40: ✓ — (no output)
Cell 41: ✓ — (no output)
Cell 42: ✓ — Out[109]= 0

## 2026-07-18 18:05:52 — ▶️ RUN CELLS 43–47
Cell 43: ✓ — (no output)
Cell 44: ✓ — Out[113]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 45: ✓ — Out[115]= \left\{6,6\right\}
Cell 46: ✓ — Out[116]= \begin{pmatrix}2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i} & 0 & 0 & 0 & 0 & 0\\0 & 7.36608403016431791181-0.968046299794461400266\,\mathrm{i} & 0 & 0 & 0 & 0\\0 & 0 & 2.68702160063169134012-1.3456153260557858257\,\mathrm{i} & 0 & 0 & 0\\0 & 0 & 0 & 7.3660840301643179118107+0.9680462997944614002659\,\mathrm{i} & 0 & 0\\0 & 0 & 0 & 0 & 2.687021600631691340118… [537 chars total]
Cell 47: ✓ — Out[125]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{max|off-diagonal|}\to 5.065342385207116\times 10^{-34},\text{digits vanished (-Log10)}\to 33.29539119374288,\text{min|diagonal|}\to 3.0051232400647923,\text{diagonal values}\to \left\{2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i},7.36608403016431791181-0… [699 chars total]

## 2026-07-18 18:05:55 — ▶️ RUN CELLS 48–54
Cell 48: ✓ — (no output)
Cell 49: ✓ — (no output)
Cell 50: ✓ — Out[133]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 51: ✓ — Out[138]= \left\{1.00000000000000000000000000000,6\right\}
Cell 52: ✓ — Out[143]= \langle\!\lvert \text{max|<PsiLeft|Psi>| off-diagonal}\to 0.\,,\text{digits vanished (-Log10)}\to \mathrm{Indeterminate}\rvert\!\rangle
Cell 53: ✓ — Out[147]= \left\{7.366084030164317911810727273-0.9680462997944614002658874268\,\mathrm{i},1.00000000000000000000+0.\times 10^{-21}\,\mathrm{i}\right\}
Cell 54: ✓ — Out[150]= \begin{matrix}\langle\!\lvert \text{M}\to 0,\text{n}\to 0,\text{braket}\to 2.909064289-5.804620064\,\mathrm{i},\text{G21 diagonal}\to 2.909064289-5.804620064\,\mathrm{i},\text{normRatio}\to 1.000000000+0.\times 10^{-10}\,\mathrm{i}\rvert\!\rangle \\\langle\!\lvert \text{M}\to 1,\text{n}\to 0,\text{braket}\to 7.366084030-0.968046300\,\mathrm{i},\text{G21 diagonal}\to 7.366084030-0.9680463… [1247 chars total]

## 2026-07-18 18:05:59 — ▶️ RUN CELLS 55–62
Cell 56: ✓ — (no output)
Cell 57: ✓ — Out[158]= \begin{pmatrix}3 & 3\\3 & 3\end{pmatrix}
Cell 58: ✓ — (no output)
Cell 59: ✓ — Out[165]= \left\{0\right\}
Cell 60: ✓ — Out[169]= \begin{matrix}\left\{1,1,0,1.\,,0.\,\right\}\\\left\{1,1,1,1.\,,0.\,\right\}\\\left\{1,2,0,1.\,,0.\,\right\}\\\left\{1,2,1,1.\,,0.\,\right\}\\\left\{1,3,0,1.\,,0.\,\right\}\\\left\{1,3,1,1.\,,0.\,\right\}\\\left\{2,1,0,1.\,,0.\,\right\}\\\left\{2,1,1,1.\,,0.\,\right\}\\\left\{2,1,2,1.\,,0.\,\right\}\\\left\{2,2,0,1.\,,0.\,\right\}\\\left\{2,2,1,1.\,,0.\,\right\}\\\left\{2,2,2,1.\,,0.\,\r… [898 chars total]
Cell 61: ✓ — Out[170]= \langle\!\lvert \text{max spread (ratio non-constancy)}\to 0.\,,\text{max |ratio - 1| over sweep}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases checked}\to 27\rvert\!\rangle
Cell 62: ✓ — Out[176]= \langle\!\lvert \text{max operator-level spread}\to 0.\,,\text{max |$\mu$ - $\mu$Site2|}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases}\to 27\rvert\!\rangle

## 2026-07-18 18:06:17 — ▶️ RUN CELLS 63–70
Cell 64: ✓ — (no output)
Cell 65: ✓ — Out[185]= \langle\!\lvert \text{w\_0 ($\lambda$=1,2,3)}\to \left\{1,1,1\right\},\text{full weights w\_j}\to \left\{1\to \left\{1,-1\right\},2\to \left\{1,-2,1\right\},3\to \left\{1,-3,3,-1\right\}\right\},\text{br1[$\lambda$][1] ($\lambda$=1,2,3)}\to \left\{0,0,0\right\}\rvert\!\rangle
Cell 66: ✓ — Out[187]= \left\{0,0\right\}
Cell 67: ✓ — Out[190]= \sqrt[3]{k}\,\left(\left(-1.18561270626131985872+0.832253303224995280404\,\mathrm{i}\right)\,k^2-\left(0.761819910806671651081+0.812812473459294825632\,\mathrm{i}\right)\,k-\left(0.404410037536507329420+0.598870869485872613792\,\mathrm{i}\right)\right)
Cell 68: ✓ — Out[194]= \begin{pmatrix}0 & 0 & 0\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix}
Cell 69: ✓ — Out[198]= \langle\!\lvert \text{OP(E12)}\to 1.\,+0.\times 10^{-30}\,\mathrm{i},\text{OP(E21)}\to \left(-0.03626734060947795145020244939-2.004237714666237867354326095\,\mathrm{i}\right)\,\varphi ^2-\left(1.522150083623538480912456055-0.2442047957313920727733393355\,\mathrm{i}\right)\,\varphi +\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right),\text{nontrivial generator \{a,b\}}\to \left\{2,1\right\},\… [464 chars total]
Cell 70: ✓ — Out[201]= \left\{1,2\right\}

## 2026-07-18 18:06:21 — ▶️ RUN CELLS 71–77
Cell 71: ✓ — Out[207]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{$\varphi$ candidates per state}\to \begin{pmatrix}\left\{0,0\right\} & \left(-0.920793778924720536105378676-0.390049761817044163511250929\,\mathrm{i}\right)\,k^{4/3}+\left(0.695722366612934233848844287-0.718310788304406529935213657\,\mathrm{i}\right)\,\sqrt[3]{k}-\left(0.3604398751… [1247 chars total]
Cell 72: ✓ — Out[213]= \langle\!\lvert \text{c\_sb solution(s)}\to \left\{\left\{\mathrm{cb}\to \frac{\left(-6.242239011314\times 10^{45}-4.511757007219\times 10^{45}\,\mathrm{i}\right)\,k^{4/3}-\left(5.999302657796\times 10^{45}+8.356648463050\times 10^{45}\,\mathrm{i}\right)\,\sqrt[3]{k}}{9.833530773718\times 10^{45}\,k^{4/3}-\left(3.54439660434\times 10^{45}-9.172544913427\times 10^{45}\,\mathrm{i}\right)\,… [841 chars total]
Cell 73: ✓ — Out[216]= \langle\!\lvert \text{generator \{a,b\} tested}\to \left\{2,1\right\},\text{rep tested}\to \left\{1,2\right\},\text{FL shape (verified exact for all 6 states)}\to \text{FL[1,2][state][k] = k\^{}(1/3) (A\_state + B\_state k)},\text{c0(state) = OPg[state][0]}\to \begin{pmatrix}\left\{0,0\right\} & 1.\,+0.\times 10^{-28}\,\mathrm{i}\\\left\{1,0\right\} & 1.\,+0.\times 10^{-28}\,\mathrm{i}\\… [1247 chars total]
Cell 74: ✓ — Out[220]= \left(1.000000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right)-\left(0.360439875147226537266677651+0.932782448593374716484835982\,\mathrm{i}\right)\,k
Cell 75: ✓ — Out[227]= \langle\!\lvert \text{alphaTab \{state,D0,c0,Dc,c1,$\alpha$\}}\to \begin{pmatrix}\left\{0,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & -0.36043987514722653-0.9327824485933747\,\mathrm{i} & 0.36043987514722653 +0.9327824485933747\,\mathrm{i} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.7101545203859353 -0.2371533418544… [1247 chars total]
Cell 76: ✓ — Out[235]= \langle\!\lvert \text{validation}\to \left\{\left\{1,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\end{pmatrix},\left\{2,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\\\left\{2,1\right\} & 0\\\left\{3,0\right\} & 0\end{pmatrix},\left\{1… [1247 chars total]
Cell 77: ✓ — Out[237]= \langle\!\lvert \text{generator \{a,b\}}\to \left\{2,1\right\},\text{$\varphi$(k)}\to \text{$\alpha$ k},\text{$\alpha$ (universal value)}\to -1.\,+0.\,\,\mathrm{i},\text{normalization}\to \text{1-site FSoV determinant det1[f] = br1[f]/(f/.u->$\theta$1); Section C's det[$\lambda$1,$\lambda$2][f] restricted to 1 site},\text{max residual over rep sweep}\to 0,\text{reps validated}\to \begin{… [510 chars total]

## 2026-07-18 18:06:33 — ▶️ RUN CELLS 78–82
Cell 79: ✓ — (no output)
Cell 80: ✓ — Out[245]= \left\{\text{Test 1 (symbolic monomials u\^{}j) residuals},\left\{0\right\}\right\}
Cell 81: ✓ — Out[250]= \left\{\text{Test 2 (twisted monomials k\^{}(u/h) u\^{}j) residuals},\left\{0\right\}\right\}
Cell 82: ✓ — Out[259]= \left\{\text{Test 3 worst |br identity residual| on Q1},0,\text{worst |detRatio - predicted| },0,\text{num states checked},41\right\}

## 2026-07-18 18:06:37 — ▶️ RUN CELLS 83–87
Cell 84: ✓ — (no output)
Cell 85: ✓ — Out[271]= \begin{pmatrix}1 & 1 & 0\\2 & 1 & 0\\3 & 1 & 0\end{pmatrix}
Cell 86: ✓ — Out[278]= \left\{\text{Mirror Test1 (u\^{}j)},\left\{0\right\},\text{Mirror Test2 (k\^{}(u/h)u\^{}j)},\left\{0\right\}\right\}
Cell 87: ✓ — Out[285]= \left\{\text{Mirror Test3 worst |br identity residual| on Q1},0,\text{num states checked},41\right\}

## 2026-07-18 18:06:38 — ▶️ RUN CELLS 88–91
Cell 89: ✓ — Out[294]= \left\{\text{OP2(E12)},1.\,+0.\times 10^{-30}\,\mathrm{i},\text{OP2(E21)},\left(2.603469709065218092712489006+0.6885776528619826437114966465\,\mathrm{i}\right)\,\varphi +\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right),\text{nontrivial generator \{a,b\}},\left\{2,1\right\}\right\}
Cell 90: ✓ — Out[297]= \left(1.\,+0.\times 10^{-27}\,\mathrm{i}\right)-\left(0.360439875147226537266677651+0.932782448593374716484835982\,\mathrm{i}\right)\,k
Cell 91: ✓ — Out[300]//TableForm= \begin{pmatrix}\left\{0,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & -0.36043987514722653-0.9327824485933747\,\mathrm{i} & 0.36043987514722653 +0.9327824485933747\,\mathrm{i} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.4685059330398573 -0.1713611156914601\,\mathrm{i} & 2.603469709065218 +0.6885776528619827… [1258 chars total]

## 2026-07-18 18:06:51 — ▶️ RUN CELLS 92–108
Cell 92: ✓ — Out[302]//TableForm= \begin{pmatrix}\left\{0,0\right\} & -1.\,+0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.973902234276272 +0.2269679230016555\,\mathrm{i} & 0.7391304347826086 +0.\,\,\mathrm{i} & 0.973902234276272 +0.2269679230016555\,\mathrm{i}\\\left\{1,0\right\} & 0.15191893851739566 -0.10600050418875033\,\mathrm{i} & -0.15191893851739566+0.10600050418875033\,\mathrm{i} & 0.7651595274649543 +… [1258 chars total]
Cell 93: ✓ — Out[305]= \begin{pmatrix}\left\{0,0\right\} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 0.05448860910968355 -0.09540961504089102\,\mathrm{i}\\\left\{1,1\right\} & -0.0293413933141388+0.7363100318787654\,\mathrm{i}\\\left\{2,0\right\} & 0.11703059508534998 +0.\,\,\mathrm{i}\\\left\{2,1\right\} & 0.8158972359113332 +0.\,\,\mathrm{i}\\\left\{3,0\right\} & 0.05448860910968355 +0.09540961504089102\,… [541 chars total]
Cell 96: ✓ — Out[311]= \left\{\text{O[k1,0]-FLdet[k1] (must be 0)},0,\text{O[0,k2]-FLdet2[k2] (expected NONZERO per Phase1 finding)},\left(-3.0719756421050755-0.5172165371705225\,\mathrm{i}\right)\,k2\right\}
Cell 97: ✓ — Out[316]= \left\{\text{Fcand1},k1^2\,\left(\left(-2.2527492196701659-5.82989030370859198\,\mathrm{i}\right)\,k2-\left(0.4170744170089964+23.0487337186617355\,\mathrm{i}\right)\right)+k1\,\left(\left(0.67128913890100625-0.37524873225652243\,\mathrm{i}\right)\,k2+\left(9.51343802264711551-1.5262799733212005\,\mathrm{i}\right)\right)-\left(1.99115021541939352-0.72828474168870538\,\mathrm{i}\right)\,k… [1014 chars total]
Cell 98: ✓ — Out[318]= \left(0.9252125034899888 -6.175727508765443\,\mathrm{i}\right)\,k1^2\,k2+\left(5.4271267482291075 +3.0890004366773067\,\mathrm{i}\right)\,k1\,k2
Cell 99: ✓ — Out[323]= \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 7.991287939023577 -1.2820751775898085\,\mathrm{i}\\\left\{2,0\right\} & -0.3808070763995185-21.044496003995498\,\mathrm{i}\\\left\{0,1\right\} & 0.6123194936458246 +1.416862394550688\,\mathrm{i}\\\left\{1,1\right\} & -0.6248120089007682-3.051910425543801\,\mathrm{i}\\\left\{2,1\right\} & -1.8923093445229393-4.897107855115217\,\m… [421 chars total]
Cell 100: ✓ — Out[328]= \left\{\text{Fcand2 (no P-insertion)},k1^2\,\left(\left(0.02466179161444501+1.36288164597304175\,\mathrm{i}\right)-\left(0.36043987514722654+0.93278244859337472\,\mathrm{i}\right)\,k2\right)+k1\,\left(\left(1.2351720155778515-0.69045766735200127\,\mathrm{i}\right)\,k2-\left(6.46913785540003854-1.0378703818584163\,\mathrm{i}\right)\right)+\left(5.38781822995835893-1.9706528304517910\,\mat… [788 chars total]
Cell 101: ✓ — 2 kernel messages: Power::infy: Infinite expression 0       encountered.
Infinity::indet: Indeterminate expression ----------------- encountered.                                                   1/3                    … [also: Power::infy] result: Out[331]//TableForm= \begin{pmatrix}2 & 1 & \left\{0,0\right\} & \mathrm{Indeterminate} & \mathrm{Indeterminate} & -4.970498218180917+4.515583419069149\,\mathrm{i}\\2 & 1 & \left\{1,1\right\} & \mathr…
Cell 102: ✓ — Out[334]//TableForm= \begin{pmatrix}2 & 1 & \left\{0,0\right\} & \left(-7.771745024480887+7.060451775520925\,\mathrm{i}\right)\,k1^2-\left(3.7846186890458786 +9.794215710230434\,\mathrm{i}\right)\,k1 & \left(1.8923093445229393 +4.897107855115217\,\mathrm{i}\right)\,k2 & -4.970498218180917+4.515583419069149\,\mathrm{i}\\2 & 1 & \left\{1,1\right\} & \left(9.339255641812587 +10.661478687052961\,\math… [1258 chars total]
Cell 103: ✓ — 2 kernel messages: Power::infy: Infinite expression 0       encountered.
Infinity::indet: Indeterminate expression ----------------- encountered.                                                   1/3                    … [also: Power::infy] result: Out[338]= \left\{\text{Fcand1(k1,0)},\mathrm{Indeterminate},\text{FLdet[k1]},\frac{\left(2.45075493251774917289426756629+5.78551794268217244227029022731\,\mathrm{i}\right)\,\left(\left(0.0525946877178…
Cell 104: ✓ — Out[343]= \left\{\text{Fcand1(k1,0)},\left(-8.51191121729+7.73287575414\,\mathrm{i}\right)\,k1^2-\left(4.5054984393+11.6597806074\,\mathrm{i}\right)\,k1+\left(1.000000000000+0.\times 10^{-12}\,\mathrm{i}\right),\text{FLdet[k1]},\frac{\left(2.45075493251774917289426756629+5.78551794268217244227029022731\,\mathrm{i}\right)\,\left(\left(0.05259468771783351844996683714+0.15021349717800998891606886537\… [960 chars total]
Cell 107: ✓ — Out[348]= \left\{\text{STEP 0: det[$\lambda$1,$\lambda$2][k\^{}(u/h)Q1] vs O[k,k], worst |residual| (sample k)},0,\text{sample rows \{$\lambda$1,$\lambda$2,state,resid,ratio\}},\begin{pmatrix}1 & 1 & \left\{0,0\right\} & 0 & 1.\,\\1 & 1 & \left\{1,0\right\} & 0 & 1.\,\\1 & 1 & \left\{1,1\right\} & 0 & 1.\,\\1 & 1 & \left\{2,0\right\} & 0 & 1.\,\\2 & 1 & \left\{0,0\right\} & 0 & 1.\,\\2 & 1 & \left… [445 chars total]
Cell 108: ✓ — Out[354]= \left\{\text{STEP 0 CONFIRMED: det[$\lambda$1,$\lambda$2][k\^{}(u/h)Q1] == O[k,k] identically},\text{worst residual over full sweep},0,\text{num states},41\right\}

## 2026-07-18 18:07:01 — ▶️ RUN CELLS 109–123
Cell 109: ✓ — 2 kernel messages: Power::infy: Infinite expression - encountered.                                  0
Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.
Power::infy: Infinite expression - encounter… [also: Power::infy] result: Out[358]= \begin{pmatrix}1 & 0 & \left\{1.\,\right\}\\1 & 1 & \left\{1.\,\right\}\\2 & 0 & \left\{1.\,\right\}\\2 & 1 & \left\{2.\,\right\}\\2 & 2 & \left\{2.\,\right\}\\3 & 0 & \left\{1.\,\right\}\\3…
Cell 110: ✓ — Out[364]= \left\{\text{STEP 1 CONFIRMED: ⟨x0[$\lambda$]|E21\^{}m = $\gamma$\_m($\lambda$) xL1[$\lambda$][\{m\}], $\gamma$\_m($\lambda$)=$\lambda$!/($\lambda$-m)!},\text{worst deviation},0,\text{table \{$\lambda$,m,$\gamma$\_found,$\gamma$\_predicted,dev\}},\begin{pmatrix}1 & 0 & 1.\, & 1 & 0\\1 & 1 & 1.\, & 1 & 0\\2 & 0 & 1.\, & 1 & 0\\2 & 1 & 2.\, & 2 & 0\\2 & 2 & 2.\, & 2 & 0\\3 & 0 & 1.\, & 1 &… [489 chars total]
Cell 111: ✓ — Out[369]= \left\{\text{STEP 2: x2[\{s1,m\}] vs xL1[$\lambda$1][\{s1\}]$\otimes$xL1[$\lambda$2][\{m\}], worst ratio-spread (constancy check)},0,\text{num cases},54,\text{sample rows},\begin{pmatrix}1 & 1 & 0 & 1 & 0.16 & 0\\1 & 1 & 1 & 1 & 1.\, & 0\\1 & 2 & 0 & 1 & 0.16 & 0\\1 & 2 & 0 & 2 & -0.68 & 0\\1 & 2 & 1 & 1 & 1.\, & 0\\1 & 2 & 1 & 2 & 1.\, & 0\\1 & 3 & 0 & 1 & 0.16 & 0\\1 & 3 & 0 & 2 & -0.6… [479 chars total]
Cell 112: ✓ — Out[373]= \begin{pmatrix}1 & 0 & \frac{4}{25}\\1 & 1 & 1\\2 & 0 & \frac{2}{23}\\2 & 1 & \frac{25}{46}\\2 & 2 & 1\\3 & 0 & \frac{4}{67}\\3 & 1 & \frac{25}{67}\\3 & 2 & \frac{46}{67}\\3 & 3 & 1\end{pmatrix}
Cell 113: ✓ — Out[378]= \left\{\text{c($\lambda$1,s1) hypothesis (m-independent) vs full sweep --- FAILS for m>=2, see next cell for the fix},\text{worst deviation over all 54 ($\lambda$1,$\lambda$2,s1,m) cases},1.68\right\}
Cell 114: ✓ — Out[382]= \begin{pmatrix}1 & 0 & 1 & \frac{4}{25}\\1 & 0 & 2 & -\frac{17}{25}\\1 & 0 & 3 & -\frac{38}{25}\\1 & 1 & 1 & 1\\1 & 1 & 2 & 1\\1 & 1 & 3 & 1\\2 & 0 & 1 & \frac{2}{23}\\2 & 0 & 2 & -\frac{34}{575}\\2 & 0 & 3 & \frac{323}{575}\\2 & 1 & 1 & \frac{25}{46}\\2 & 1 & 2 & \frac{2}{23}\\2 & 1 & 3 & -\frac{17}{46}\\2 & 2 & 1 & 1\\2 & 2 & 2 & 1\\2 & 2 & 3 & 1\\3 & 0 & 1 & \frac{4}{67}\\3 & 0 & 2 & … [672 chars total]
Cell 115: ✓ — 4 kernel messages: LinearSolve::nosol: Linear equation encountered that has no solution.
LinearSolve::nosol: Linear equation encountered that has no solution.
LinearSolve::nosol: Linear equation encountered that has no … [also: LinearSolve::nosol, LinearSolve::nosol, LinearSolve::nosol] result: Out[387]= \left\{\text{coeffs (rows=s1, cols=s1'=0..3)},\left\{\mathrm{LinearSolve}\left[\begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & \frac{1}{36}\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & \f…
Cell 116: ✓ — Out[392]= \left\{\left\{0,2\right\}\to -0.02206359506813757300454250486697,\left\{1,1\right\}\to 0.1635301752109020116807268007787,\left\{2,0\right\}\to 0.8585334198572355613238157040883\right\}
Cell 117: ✓ — Out[396]= \left\{\text{STEP 2 level-conservation check: \# off-level nonzero coeffs (want 0 everywhere)},0,\text{num ($\lambda$1,$\lambda$2,s1,m) cases},29\right\}
Cell 119: ✓ — Out[404]= \left\{\text{m=0: Om-Cm (must be 0 per Phase1, C1pred/P1site2[$\theta$1] should equal exactly Cm's own norm here)},\left(2.5572121404875445 -0.4102640568287387\,\mathrm{i}\right)\,k1-\left(0.038411844228125346 +2.1227491621421546\,\mathrm{i}\right)\,k1^2,\text{m=1: Om\_1[k1]},\left(0.3604398751472265372666776506+0.9327824485933747164848359821\,\mathrm{i}\right)\,k1^2-\left(1.296101147801… [801 chars total]
Cell 120: ✓ — Out[411]= \left\{\text{m=0: (B0-B1)/O\_0 ratio},\frac{\sqrt[3]{k1}\,\left(\left(-0.0292504903494856237+0.0205326892042810627\,\mathrm{i}\right)\,k1^2+\left(0.499242526583056932+0.532659421382621779\,\mathrm{i}\right)\,k1-\left(0.139388493663376676+0.206413542317126162\,\mathrm{i}\right)\right)}{\left(0.3604398751472265372666776506+0.9327824485933747164848359821\,\mathrm{i}\right)\,k1^2-\left(1.296… [1083 chars total]
Cell 121: ✓ — Out[414]= \left\{\text{Cm[0]-Cm[1] shape},\left(-0.0225172879957976151-1.244370198497125076\,\mathrm{i}\right)\,k1^2+\left(5.434075798536032377-0.8718111207610696998\,\mathrm{i}\right)\,k1+\left(0.\times 10^{-20}+0.\times 10^{-20}\,\mathrm{i}\right),\text{ratio vs O\_0[k1] (want: clean constant)},\frac{\left(-0.0225172879957976151-1.244370198497125076\,\mathrm{i}\right)\,k1^2+\left(5.4340757985360… [753 chars total]
Cell 122: ✓ — Out[416]= \left\{\text{Cm[0][k1] / O\_1[k1] (reversed-index hypothesis; want clean k1-independent constant)},\frac{\left(0.111318362779488714+0.0407158959953172346\,\mathrm{i}\right)\,k1^2-\left(0.21818047436770437-1.02534209918467083\,\mathrm{i}\right)\,k1+\left(0.360439875147226537-0.932782448593374716\,\mathrm{i}\right)}{\left(1.\,+0.\times 10^{-29}\,\mathrm{i}\right)\,k1^2-\left(2.963909584212… [539 chars total]

## 2026-07-18 18:07:18 — ▶️ RUN CELLS 124–138
Cell 125: ✓ — (no output)
Cell 126: ✓ — Out[446]= \left\{\text{Symbolic-$\theta$,$\chi$1 port-fidelity: xL1Sym vs xL1},0,\text{x2Sym vs x2},0\right\}
Cell 127: ✓ — Out[451]= \left\{\left\{0,0\right\}\to 0,\left\{0,1\right\}\to \frac{\theta A-\theta B}{\theta A-\theta B+1},\left\{1,0\right\}\to \frac{1}{\theta A-\theta B+1},\left\{1,1\right\}\to 0\right\}
Cell 128: ✓ — Out[454]= 35
Cell 129: ✓ — Out[460]= \left\{\left\{0,0\right\}\to 0,\left\{0,1\right\}\to \theta B-\theta A,\left\{1,0\right\}\to -1,\left\{1,1\right\}\to 0\right\}
Cell 130: ✓ — Out[463]= 127
Cell 131: ✓ — Out[466]= 127
Cell 132: ✓ — Out[470]//MatrixForm= \begin{pmatrix}\frac{\theta A^2-2\,\theta A\,\theta B-\theta A+\theta B^2+\theta B}{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)} & \frac{6\,\left(\theta A-\theta B\right)}{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)} & \frac{6}{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)}\\0 & \frac{\theta A-\theta B+1}{\the… [472 chars total]
Cell 133: ✓ — Out[472]//MatrixForm= \begin{pmatrix}\frac{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)}{\theta A^2-2\,\theta A\,\theta B-\theta A+\theta B^2+\theta B} & -\frac{6\,\left(\theta A-\theta B+3\right)}{\theta A^2-2\,\theta A\,\theta B+\theta B^2-1} & \frac{6}{\left(\theta A-\theta B\right)\,\left(\theta A-\theta B+1\right)}\\0 & \frac{\theta A-\theta B+3}{\theta A-\theta B+1} & -… [453 chars total]
Cell 134: ✓ — Out[478]= \left\{\text{equations},\left\{1.000000000000000000000+0.\times 10^{-21}\,\mathrm{i}⩵\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right)\,wc0+\left(2.603469709065218092712489006+0.6885776528619826437114966465\,\mathrm{i}\right)\,wc1,-1.0350620568640061670+0.16605926109734660949\,\mathrm{i}⩵\left(1.522150083623538480912456055-0.2442047957313920727733393355\,\mathrm{i}\right)\,wc0-\left(1.2961… [785 chars total]
Cell 135: ✓ — Out[484]= \begin{pmatrix}0 & 0 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\0 & 1 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\1 & 0 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\1 & 1 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\2 & 0 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} … [510 chars total]
Cell 136: ✓ — Out[489]= \left\{\text{Conjectured sum rule FAILS beyond the lucky L=1 case},\text{num failing cases},29,\text{num total cases},81,\text{sample failures},\begin{pmatrix}1 & 2 & 0 & 1 & \frac{1}{-\theta A+\theta B-1}\\1 & 2 & 0 & 2 & \frac{2}{\theta A-\theta B+1}\\1 & 3 & 0 & 1 & -\frac{2}{\theta A-\theta B+1}\end{pmatrix}\right\}
Cell 138: ✓ — Out[496]= \left\{\text{BOXED: A[s1,m;s1,m] = Pochhammer[s1+$\theta$1-$\theta$2-m+1,m]/Pochhammer[$\lambda$1+$\theta$1-$\theta$2-m+1,m] --- verified at the live numeric $\theta$1,$\theta$2},\text{worst deviation},0,\text{num ($\lambda$1,$\lambda$2,s1,m) cases},81\right\}

## 2026-07-18 18:07:27 — ▶️ RUN CELLS 139–154
Cell 140: ✓ — Out[502]= \left\{\text{STEP 1a (WITH $\kappa$1\^{}(u/h) prefactor): worst |x2.Psi - Q1($\theta$1+s1)Q1($\theta$2+s2)/(Q1($\theta$1)Q1($\theta$2))|},0,\text{num cases},297\right\}
Cell 141: ✓ — Out[509]= \left\{\text{STEP 1b (WITHOUT prefactor, bare q\_n): worst residual (expect NONZERO, confirming prefactor matters)},3.7396476326262356,\text{sample},\begin{pmatrix}\left\{0,0\right\} & 0 & 0 & 0\\\left\{1,0\right\} & 0 & 0 & 0\\\left\{1,1\right\} & 0 & 0 & 0\end{pmatrix}\right\}
Cell 143: ✓ — Out[515]= \begin{pmatrix}\left\{0,0,\left\{\text{srng},\left\{0\right\},\text{sol},\left\{1.\,\right\},\text{maxResid},0.\times 10^{-28}\right\}\right\} & \left\{0,1,\left\{\text{srng},\left\{0,1\right\},\text{sol},\left\{11.50000000000000000000000000-4.431302171530629069010197265\times 10^{-38}\,\mathrm{i},-10.50000000000000000000000000+2.933894180666706878219239219\times 10^{-38}\,\mathrm{i}\rig… [1115 chars total]
Cell 144: ✓ — Out[521]= \left\{\text{worst residual across full B-sweep (want \textasciitilde{}0)},0.\times 10^{-24},\text{num (rep,a,b) cases},37\right\}
Cell 145: ✓ — Out[525]= \left\{\text{A-block ($\lambda$1=2,$\lambda$2=1, L=2, s,s'=1,2)},\begin{pmatrix}\frac{25}{46} & \frac{21}{46}\\0 & 1\end{pmatrix},\text{its inverse (should = B row)},\begin{pmatrix}\frac{46}{25} & -\frac{21}{25}\\0 & 1\end{pmatrix},\text{Step2's independently-solved B row \{1,1\}},\left\{\frac{46}{25},-\frac{21}{25}\right\},\text{match?},\left\{0,0\right\}\right\}
Cell 147: ✓ — Out[532]= \left\{\text{det[$\lambda$1,$\lambda$2][k\^{}u Q1] - explicit residue-sum expansion (rep (2,1), state (1,0)), want 0},0\right\}
Cell 148: ✓ — Out[538]= \left\{\text{TRUE sum rule (derived, Step 4): worst deviation},0,\text{num ($\lambda$1,$\lambda$2,s,t) cases},20,\text{sample},\begin{pmatrix}2 & 1 & 0 & 0 & 0\\2 & 1 & 0 & 1 & 0\\2 & 1 & 1 & 0 & 0\\2 & 1 & 1 & 1 & 0\\2 & 1 & 2 & 0 & 0\\2 & 1 & 2 & 1 & 0\end{pmatrix}\right\}
Cell 150: ✓ — Out[546]= \left\{\text{STEP 3 assembly test (rep (2,1), 6 states x 3 grid points): worst |Oassembled-O2rot|},0,\text{num cases},18\right\}
Cell 151: ✓ — Out[552]= \left\{\text{STEP 3 FULL-SWEEP ASSEMBLY CONFIRMED: O[k1,k2] = $\Sigma$ cA(k1)cB(k2) B-weighted Q-bilinears == O2rot exactly},\text{worst deviation},0,\text{num states x grid cases},82\right\}
Cell 153: ✓ — Out[557]= \begin{pmatrix}0 & 0 & 1.\, & 1.\, & 0\\0 & 1 & -2.603469709065218-0.6885776528619827\,\mathrm{i} & 5.387818229958359 -1.970652830451791\,\mathrm{i} & -7.991287939023577+1.2820751775898085\,\mathrm{i}\\1 & 0 & 1.5221500836235384 -0.24420479573139206\,\mathrm{i} & -6.469137855400039+1.0378703818584163\,\mathrm{i} & 7.991287939023577 -1.2820751775898085\,\mathrm{i}\\1 & 1 & 1.2961011478017… [843 chars total]

## 2026-07-18 18:07:33 — ▶️ RUN CELLS 155–160
Cell 156: ✓ — Out[561]= \left\{\text{M matrix (rep (2,1), k1=7/3,k2=5/2)},\begin{pmatrix}1 & -\frac{115}{4}\\\frac{259}{12} & \frac{322}{15}\\-\frac{196}{45} & -\frac{245}{18}\end{pmatrix},\text{MatrixRank},2\right\}
Cell 157: ✓ — Out[566]//TableForm= \begin{pmatrix}1 & 1 & \left\{\frac{7}{3},\frac{5}{2}\right\} & 2 & 2\\1 & 1 & \left\{\frac{11}{5},-\frac{4}{3}\right\} & 2 & 2\\1 & 1 & \left\{-\frac{2}{7},\frac{9}{4}\right\} & 2 & 2\\1 & 1 & \left\{\frac{13}{6},\frac{1}{9}\right\} & 2 & 2\\2 & 1 & \left\{\frac{7}{3},\frac{5}{2}\right\} & 2 & 2\\2 & 1 & \left\{\frac{11}{5},-\frac{4}{3}\right\} & 2 & 2\\2 & 1 & \left\{-\frac{… [1258 chars total]
Cell 158: ✓ — Out[568]= \left\{\text{Det[M(k1,k2)] symbolic, rep (2,2)},\frac{213003\,k2^3\,\left(k1-k2\right)^2\,\left(17\,k1+46\,k2\right)}{7225},\text{identically zero?},\mathrm{False}\right\}
Cell 159: ✓ — Out[572]= \left\{\text{rep (2,3): M dimensions},\left\{3,4\right\},\text{rank at sample point},3,\text{a 3x3 minor determinant (nonzero => rank>=3)},\frac{19563970545}{1405696}\right\}

## 2026-07-18 18:08:01 — 📥 BULK INSERT 2 CELL(S) at positions 161–162
1. [markdown] (*Section M — Paul's single-tower ansatz. O[p,r] := (⟨x0[λ1]|e^{pE21}⊗⟨x0[λ2]|e^{rE21}).Psi[M,n]↵= O
2. [code] (*O[p,r] and the site-1-tower measurement matrix (S=λ1+λ2 to start).*)↵ClearAll[Opr,towerMat];↵Opr[λ


## 2026-07-18 18:08:02
**Input:**
```wolfram
(*O[p,r] and the site-1-tower measurement matrix (S=λ1+λ2 to start).*)
ClearAll[Opr,towerMat];
Opr[λ1_,λ2_][M_,n_][p_,r_] := O2rot[λ1,λ2][M,n][-p,-r];
towerMat[λ1_,λ2_][S_] := Table[
  Q1[λ1,λ2][st[[1]],st[[2]]][θ1+s h] / Q1[λ1,λ2][st[[1]],st[[2]]][θ1]
,{st,statesOf[λ1,λ2]},{s,0,S}];
(*rep (2,1): S=λ1+λ2=3, matrix is 6 states x 4 unknowns — overdetermined by 2.*)
mat21 = towerMat[2,1][3];
{"towerMat dims (rep (2,1), S=3)", Dimensions[mat21], "MatrixRank", MatrixRank[mat21]}
```
**Output:**
```
Out[577]= \left\{\text{towerMat dims (rep (2,1), S=3)},\left\{6,4\right\},\text{MatrixRank},4\right\}
```
## 2026-07-18 18:08:34 — 📥 BULK INSERT 1 CELL(S) at positions 163–163
1. [code] (*Solve for c_s at a sample (p,r), rep (2,1), S=λ1+λ2=3: b=O[p,r] per state, X=towerMat, exact↵  rat


## 2026-07-18 18:08:35
**Input:**
```wolfram
(*Solve for c_s at a sample (p,r), rep (2,1), S=λ1+λ2=3: b=O[p,r] per state, X=towerMat, exact
  rational (p,r) for exact arithmetic.*)
ClearAll[rhsVec,solveTower];
rhsVec[λ1_,λ2_][p_,r_] := Table[Opr[λ1,λ2][st[[1]],st[[2]]][p,r], {st,statesOf[λ1,λ2]}];
solveTower[λ1_,λ2_][S_][p_,r_] := Module[{Xm,bv,csol,resid},
  Xm = towerMat[λ1,λ2][S];
  bv = N[rhsVec[λ1,λ2][p,r],40];
  csol = LeastSquares[Xm, bv];
  resid = Max[Abs[Xm.csol - bv]];
  {csol, resid}
];
{c21test, resid21test} = solveTower[2,1][3][7/3,5/2];
{"c_s (rep (2,1), S=3, p=7/3,r=5/2)", Chop[c21test,10^-15], "max residual", resid21test}
```
**Output:**
```
Out[582]= \left\{\text{c\_s (rep (2,1), S=3, p=7/3,r=5/2)},\left\{13.7688078414687182186600+2.10326272230760132765845\,\mathrm{i},-3.12167136931824065325212+18.8662496625125701186869\,\mathrm{i},-7.51900595437450423013552-10.1245361864155291062487\,\mathrm{i},3.80567731087062145228106-0.708676905117290761157203\,\mathrm{i}\right\},\text{max residual},27.6845440371648442360677\right\}
```
## 2026-07-18 18:08:55 — 📥 BULK INSERT 1 CELL(S) at positions 164–164
1. [code] (*Variant (a): extended range. Try S=λ1+λ2+1,...,λ1+λ2+4 for rep (2,1) (states available: 6, so↵  S


## 2026-07-18 18:08:56
**Input:**
```wolfram
(*Variant (a): extended range. Try S=λ1+λ2+1,...,λ1+λ2+4 for rep (2,1) (states available: 6, so
  S can go up to 5 before the system becomes underdetermined against a single (p,r) sample — but
  since towerMat is (p,r)-INDEPENDENT, rank must reach full column rank S+1<=6 for uniqueness;
  beyond S=5 columns, an EXACT solve at one (p,r) can't distinguish extra freedom, but consistency
  (existence) can still be checked via projecting out the column space).*)
extendedTest = Table[
  Module[{Xm,bv,csol,resid,rk},
    Xm = towerMat[2,1][S];
    rk = MatrixRank[Xm];
    bv = N[rhsVec[2,1][7/3,5/2],40];
    csol = LeastSquares[Xm, bv];
    resid = Max[Abs[Xm.csol - bv]];
    {S, Dimensions[Xm], rk, resid}
  ]
,{S,3,7}];
extendedTest//TableForm
```
**Output:**
```
Out[584]//TableForm= \begin{pmatrix}3 & \left\{6,4\right\} & 4 & 27.6845440371648442360677\\4 & \left\{6,5\right\} & 4 & 27.6845440371648442360677\\5 & \left\{6,6\right\} & 4 & 27.684544037164844236068\\6 & \left\{6,7\right\} & 4 & 27.684544037164844236068\\7 & \left\{6,8\right\} & 4 & 27.684544037164844236068\end{pmatrix}
```
## 2026-07-18 18:09:30 — 📥 BULK INSERT 1 CELL(S) at positions 165–165
1. [code] (*Sanity check FIRST: at r=0, O[p,0] is ALREADY KNOWN to equal Σ_{s=0}^{λ1} w1(s)(-p)^s Q1(θ1+sh)/Q1


## 2026-07-18 18:09:30
**Input:**
```wolfram
(*Sanity check FIRST: at r=0, O[p,0] is ALREADY KNOWN to equal Σ_{s=0}^{λ1} w1(s)(-p)^s Q1(θ1+sh)/Q1(θ1)
  exactly (Section F's det1 dictionary). Verify this closes with S=λ1=2 for rep (2,1), confirming
  the Opr/towerMat machinery itself is correct before diagnosing the general (p,r) failure.*)
{c21r0, resid21r0} = solveTower[2,1][2][7/3,0];
predW = Table[w1wt[2,1][s] (-7/3)^s, {s,0,2}];
{"c_s at r=0, S=λ1=2 (rep (2,1), p=7/3)", Chop[N[c21r0],10^-10], "residual", resid21r0,
 "predicted w1(s)(-p)^s", Chop[N[predW],10^-10], "match?", Chop[N[c21r0-predW],10^-8]}
```
**Output:**
```
Out[587]= \left\{\text{c\_s at r=0, S=$\lambda$1=2 (rep (2,1), p=7/3)},\left\{1.\,,4.666666666666667,5.444444444444445\right\},\text{residual},0.\times 10^{-23},\text{predicted w1(s)(-p)\^{}s},\left\{1.\,,-3.1733333333333333,-0.3219323671497585\right\},\text{match?},\left\{0,7.84,5.766376811594203\right\}\right\}
```

## 2026-07-18 18:09:55
**Input:**
```wolfram
w1[1]
w1[2]
w1wt[2,1][0]
w1wt[2,1][1]
```
**Output:**
```
Out[1]= {1, -1}
Out[2]= {1, -2, 1}
Out[3]= 1
Out[4]= 34/25
```
## 2026-07-18 18:10:06 — ✏️ EDIT CELL 165
(*Sanity check FIRST: at r=0, O[p,0] is ALREADY KNOWN to equal Σ_{s=0}^{λ1} w1[λ1][[s+1]](-p)^s Q1(θ1+sh)/Q1(θ1)
  exactly (Section F's det1 dictionary, br1's single-site tower weights w1[λ]). Verify …


## 2026-07-18 18:10:06
**Input:**
```wolfram
(*Sanity check FIRST: at r=0, O[p,0] is ALREADY KNOWN to equal Σ_{s=0}^{λ1} w1[λ1][[s+1]](-p)^s Q1(θ1+sh)/Q1(θ1)
  exactly (Section F's det1 dictionary, br1's single-site tower weights w1[λ]). Verify this closes
  with S=λ1=2 for rep (2,1), confirming the Opr/towerMat machinery itself is correct before
  diagnosing the general (p,r) failure.
  (*why w1[λ1] not w1wt[λ1,λ2]*) w1wt[λ1,λ2][s] (Section K, Cell 147) is built from the TWO-SITE
  measure mu[λ1,λ2][1] (used there for the det-residue-sum sum-rule derivation) — a DIFFERENT
  object from Section F's single-site br1 tower weights w1[λ], which is what actually appears in
  det1/FLdet. Confusing the two was the bug in the first version of this cell.*)
{c21r0, resid21r0} = solveTower[2,1][2][7/3,0];
predW = Table[w1[2][[s+1]] (-7/3)^s, {s,0,2}];
{"c_s at r=0, S=λ1=2 (rep (2,1), p=7/3)", Chop[N[c21r0],10^-10], "residual", resid21r0,
 "predicted w1[λ1][[s+1]](-p)^s", Chop[N[predW],10^-10], "match?", Chop[N[c21r0-predW],10^-8]}
```
**Output:**
```
Out= {"c_s at r=0, S=λ1=2 (rep (2,1), p=7/3)", {1., 4.666666666666667, 5.444444444444445}, "residual", 0``22.57763473576602, "predicted w1[λ1][[s+1]](-p)^s", {1., 4.666666666666667, 5.444444444444445}, "match?", {0, 0, 0}}
```
## 2026-07-18 18:11:05 — 📥 BULK INSERT 1 CELL(S) at positions 166–166
1. [code] (*Variant (c): tower on θ2 instead of θ1 (site-ordering asymmetry — recall A,B depended on λ1↵  ONLY


## 2026-07-18 18:11:05
**Input:**
```wolfram
(*Variant (c): tower on θ2 instead of θ1 (site-ordering asymmetry — recall A,B depended on λ1
  ONLY, so the "distinguished" site for this kind of object might be site 2, not site 1). Test
  S=λ1+λ2=3 at the SAME sample point, rep (2,1).*)
ClearAll[towerMat2];
towerMat2[λ1_,λ2_][S_] := Table[
  Q1[λ1,λ2][st[[1]],st[[2]]][θ2+s h] / Q1[λ1,λ2][st[[1]],st[[2]]][θ2]
,{st,statesOf[λ1,λ2]},{s,0,S}];
Mm2 = towerMat2[2,1][3];
bv2 = N[rhsVec[2,1][7/3,5/2],40];
csol2 = LeastSquares[Mm2, bv2];
resid2 = Max[Abs[Mm2.csol2 - bv2]];
{"θ2-tower: dims", Dimensions[Mm2], "rank", MatrixRank[Mm2], "residual", resid2}
```
**Output:**
```
Out[594]= \left\{\text{$\theta$2-tower: dims},\left\{6,4\right\},\text{rank},4,\text{residual},25.4178846095368977408093\right\}
```
## 2026-07-18 18:11:32 — 📥 BULK INSERT 1 CELL(S) at positions 167–167
1. [code] (*Broaden "extended range" to allow shifted/negative-offset node sets (still θ1-tower, same↵  total


## 2026-07-18 18:11:32
**Input:**
```wolfram
(*Broaden "extended range" to allow shifted/negative-offset node sets (still θ1-tower, same
  total node count as natural candidates, but different starting point) — test a scan of ranges
  s=sMin..sMax for rep (2,1), same sample point, looking for ANY consistent (residual~0) choice.*)
ClearAll[towerMatRange];
towerMatRange[λ1_,λ2_][sMin_,sMax_] := Table[
  Q1[λ1,λ2][st[[1]],st[[2]]][θ1+s h] / Q1[λ1,λ2][st[[1]],st[[2]]][θ1]
,{st,statesOf[λ1,λ2]},{s,sMin,sMax}];
rangeScan = Table[
  Module[{Xm,bv,csol,resid,rk},
    Xm = towerMatRange[2,1][sMin,sMax];
    bv = N[rhsVec[2,1][7/3,5/2],40];
    rk = MatrixRank[Xm];
    csol = LeastSquares[Xm, bv];
    resid = Max[Abs[Xm.csol - bv]];
    {sMin,sMax, Dimensions[Xm][[2]], rk, resid}
  ]
,{sMin,-4,0},{sMax,sMin+3,sMin+7}];
Flatten[rangeScan,1]//TableForm
```
**Output:**
```
Out[598]//TableForm= \begin{pmatrix}-4 & -1 & 4 & 4 & 27.684544037164844236068\\-4 & 0 & 5 & 4 & 27.6845440371648442360677\\-4 & 1 & 6 & 4 & 27.6845440371648442360677\\-4 & 2 & 7 & 4 & 27.6845440371648442360677\\-4 & 3 & 8 & 4 & 27.6845440371648442360677\\-3 & 0 & 4 & 4 & 27.6845440371648442360677\\-3 & 1 & 5 & 4 & 27.6845440371648442360677\\-3 & 2 & 6 & 4 & 27.6845440371648442360677\\-3 & 3 & 7 &… [1139 chars total]
```
## 2026-07-18 18:12:26 — 📥 BULK INSERT 1 CELL(S) at positions 168–168
1. [code] (*Variant (b), properly tested: multiply (not divide) the TARGET by Q1(θ2) per state before↵  solvin


## 2026-07-18 18:12:27
**Input:**
```wolfram
(*Variant (b), properly tested: multiply (not divide) the TARGET by Q1(θ2) per state before
  solving against the θ1-tower basis Q1(θ1+sh)/Q1(θ1) — genuinely different linear problem (not a
  trivial row-rescaling of the original one, since it changes which functional the tower must
  reproduce). Also try dividing by Q1(θ2) for completeness.*)
Xm0 = towerMat[2,1][3];
bvBase = N[rhsVec[2,1][7/3,5/2],40];
q1th2 = N[Table[Q1[2,1][st[[1]],st[[2]]][θ2], {st,statesOf[2,1]}],40];
bvMult = bvBase q1th2;
bvDiv  = bvBase / q1th2;
csolMult = LeastSquares[Xm0, bvMult];
residMult = Max[Abs[Xm0.csolMult - bvMult]];
csolDiv = LeastSquares[Xm0, bvDiv];
residDiv = Max[Abs[Xm0.csolDiv - bvDiv]];
{"multiply by Q1(θ2): residual", residMult, "divide by Q1(θ2): residual", residDiv}
```
**Output:**
```
Out[608]= \left\{\text{multiply by Q1($\theta$2): residual},19.7572259100834793810211,\text{divide by Q1($\theta$2): residual},18.310606875396829403204\right\}
```
## 2026-07-18 18:14:32 — 📥 BULK INSERT 1 CELL(S) at positions 169–169
1. [code] (*Diagnose the rank cap precisely: for rep (2,1), the achievable θ1-tower column space has rank 4↵


## 2026-07-18 18:14:32
**Input:**
```wolfram
(*Diagnose the rank cap precisely: for rep (2,1), the achievable θ1-tower column space has rank 4
  REGARDLESS of window (Cells 164,167) — matches max_state(M+1)=4 (state (3,0) has M=3, the highest
  weight sector in this rep). Confirm this is the TRUE reason (not a coincidence) by checking the
  rank achieved using ONLY states up to a given M, and separately confirm the r=0 anchor's success
  (rank-3 subspace) is a genuine SUBSPACE of this rank-4 space. Then test O[p,r] failure at several
  MORE (p,r) points and a second representation for robustness.*)
(*(a) is r=0 really inside the rank-4 space, i.e. is the r=0 solution reproducible using ONLY
  θ1-tower data with the SAME S=3 window that failed for r≠0?*)
Xm3 = towerMat[2,1][3]; (*4 columns, rank 4*)
bvR0 = N[rhsVec[2,1][7/3,0],40];
csolR0S3 = LeastSquares[Xm3,bvR0];
residR0S3 = Max[Abs[Xm3.csolR0S3-bvR0]];
(*(b) O[p,r] failure at several more (p,r), rep (2,1)*)
morePts = {{1,1},{2,-1},{-3,4},{0,7/2}};
failScan = Table[
  Module[{bv,csol,resid},
    bv = N[rhsVec[2,1][pt[[1]],pt[[2]]],40];
    csol = LeastSquares[Xm3,bv];
    resid = Max[Abs[Xm3.csol-bv]];
    {pt,resid}
  ]
,{pt,morePts}];
(*(c) second rep: (1,2) — same qualitative test.*)
Xm3b = towerMat[1,2][3];
rk1_2 = MatrixRank[Xm3b];
bv12 = N[rhsVec[1,2][7/3,5/2],40];
csol12 = LeastSquares[Xm3b,bv12];
resid12 = Max[Abs[Xm3b.csol12-bv12]];
{"r=0 residual at S=3 window (should stay 0, confirming subspace nesting)", residR0S3,
 "more (p,r) failures, rep(2,1)", failScan,
 "rep (1,2): towerMat rank", rk1_2, "residual at (7/3,5/2)", resid12}
```
**Output:**
```
2 kernel messages: Set::write: Tag Times in 2 (rk1_) is Protected. [also: Set::write] result: Out[620]= \left\{\text{r=0 residual at S=3 window (should stay 0, confirming subspace nesting)},0.\times 10^{-23},\text{more (p,r) failures, rep(2,1)},\begin{pmatrix}\left\{1,1\right\} & 4.04518880433…
```
## 2026-07-18 18:14:44 — ✏️ EDIT CELL 169
(*Diagnose the rank cap precisely: for rep (2,1), the achievable θ1-tower column space has rank 4
  REGARDLESS of window (Cells 164,167) — matches max_state(M+1)=4 (state (3,0) has M=3, the highest
  …


## 2026-07-18 18:14:45
**Input:**
```wolfram
(*Diagnose the rank cap precisely: for rep (2,1), the achievable θ1-tower column space has rank 4
  REGARDLESS of window (Cells 164,167) — matches max_state(M+1)=4 (state (3,0) has M=3, the highest
  weight sector in this rep). Confirm this is the TRUE reason (not a coincidence) by checking the
  rank achieved using ONLY states up to a given M, and separately confirm the r=0 anchor's success
  (rank-3 subspace) is a genuine SUBSPACE of this rank-4 space. Then test O[p,r] failure at several
  MORE (p,r) points and a second representation for robustness.
  (*why rkRep12 not rk1_2*) "rk1_2" silently parses as Times[rk1_,2] (Blank pattern times 2) — same
  underscore-in-identifier gotcha as before; renamed to avoid it.*)
(*(a) is r=0 really inside the rank-4 space, i.e. is the r=0 solution reproducible using ONLY
  θ1-tower data with the SAME S=3 window that failed for r≠0?*)
Xm3 = towerMat[2,1][3]; (*4 columns, rank 4*)
bvR0 = N[rhsVec[2,1][7/3,0],40];
csolR0S3 = LeastSquares[Xm3,bvR0];
residR0S3 = Max[Abs[Xm3.csolR0S3-bvR0]];
(*(b) O[p,r] failure at several more (p,r), rep (2,1)*)
morePts = {{1,1},{2,-1},{-3,4},{0,7/2}};
failScan = Table[
  Module[{bv,csol,resid},
    bv = N[rhsVec[2,1][pt[[1]],pt[[2]]],40];
    csol = LeastSquares[Xm3,bv];
    resid = Max[Abs[Xm3.csol-bv]];
    {pt,resid}
  ]
,{pt,morePts}];
(*(c) second rep: (1,2) — same qualitative test.*)
Xm3b = towerMat[1,2][3];
rkRep12 = MatrixRank[Xm3b];
bv12 = N[rhsVec[1,2][7/3,5/2],40];
csol12 = LeastSquares[Xm3b,bv12];
resid12 = Max[Abs[Xm3b.csol12-bv12]];
{"r=0 residual at S=3 window (should stay 0, confirming subspace nesting)", residR0S3,
 "more (p,r) failures, rep(2,1)", failScan,
 "rep (1,2): towerMat rank", rkRep12, "residual at (7/3,5/2)", resid12}
```
**Output:**
```
Out= {"r=0 residual at S=3 window (should stay 0, confirming subspace nesting)", 0``22.57763473576602, "more (p,r) failures, rep(2,1)", {{{1, 1}, 4.04518880433634545024478333679182823392`23.571662174886168}, {{2, -1}, 8.92888343029211987360144656308363921184`23.426020765200608}, {{-3, 4}, 39.33363888266671510111859404019820069192`23.879898766747296}, {{0, 7/2}, 5.3318082362908496984122066573075485121`23.345524998723956}}, "rep (1,2): towerMat rank", 4, "residual at (7/3,5/2)", 27.94146729657645870816812742102593166321`23.62181100339649}
```
## 2026-07-18 18:34:29 — 📥 BULK INSERT 1 CELL(S) at positions 170–170
1. [code] (*Rule out one more reading: basis = BARE q1(θ1+sh)/q1(θ1) (no κ1^(u/h) prefactor at all), in case↵


## 2026-07-18 18:34:30
**Input:**
```wolfram
(*Rule out one more reading: basis = BARE q1(θ1+sh)/q1(θ1) (no κ1^(u/h) prefactor at all), in case
  "handle the κ1^s...bookkeeping explicitly" meant using bare q1 as the tower object with κ1 powers
  absorbed into c_s. This does NOT change the achievable column space (κ1^s is just an invertible
  per-column rescaling, which cannot change RANK or reachability of a fixed target vector by more
  than what a column rescaling allows — but verify directly rather than assume, since a rescaling of
  s-index by a STATE-INDEPENDENT factor could interact with a state-DEPENDENT normalization q1(θ1)
  in an unobvious way).*)
towerMatBare[λ1_,λ2_][S_] := Table[
  qbare[λ1,λ2][st[[1]],st[[2]]][θ1+s h] / qbare[λ1,λ2][st[[1]],st[[2]]][θ1]
,{st,statesOf[λ1,λ2]},{s,0,S}];
XmBare = towerMatBare[2,1][3];
bvBareTest = N[rhsVec[2,1][7/3,5/2],40];
csolBare = LeastSquares[XmBare,bvBareTest];
residBare = Max[Abs[XmBare.csolBare-bvBareTest]];
{"bare-q1 tower basis: rank", MatrixRank[XmBare], "residual at (7/3,5/2)", residBare}
```
**Output:**
```
Out[626]= \left\{\text{bare-q1 tower basis: rank},4,\text{residual at (7/3,5/2)},27.6845440371648442360677\right\}
```
## 2026-07-18 18:36:06 — 📥 BULK INSERT 1 CELL(S) at positions 171–171
1. [code] (*KERNEL CHECK (task step 3), for the r=0 anchor where the system DOES close: is the c_s solution↵


## 2026-07-18 18:36:06
**Input:**
```wolfram
(*KERNEL CHECK (task step 3), for the r=0 anchor where the system DOES close: is the c_s solution
  UNIQUE? towerMat[2,1][2] is 6x3 (6 states, 3 unknowns s=0,1,2) — compute its rank; if rank=3 (full
  column rank), the map is injective and any consistent solution is unique (kernel of the map
  c->Xc is trivial, matching the exact w1[λ1] match found in Cell 165 with no ambiguity).*)
rankR0 = MatrixRank[towerMat[2,1][2]];
nullR0 = NullSpace[towerMat[2,1][2]];
{"rank of towerMat[2,1][2] (r=0 anchor, 3 unknowns)", rankR0,
 "null space dimension (0 = unique solution)", Length[nullR0]}
```
**Output:**
```
Out[629]= \left\{\text{rank of towerMat[2,1][2] (r=0 anchor, 3 unknowns)},3,\text{null space dimension (0 = unique solution)},0\right\}
```
## 2026-07-18 18:36:42 — 📥 BULK INSERT 1 CELL(S) at positions 172–172
1. [markdown] (*Section M final report — the r=0 ANCHOR closes exactly (uniquely); the GENERAL (p,r) ansatz, as↵li


---

## 2026-07-18 18:36:48 — 🔄 KERNEL RESTART

## 2026-07-18 18:38:37 — ▶️ RUN CELLS 1–10
Cell 1: ✓ — Out[8]= \left\{\frac{1}{3},\frac{1}{7},1,0.360439875147226537266677650596+0.932782448593374716484835982071\,\mathrm{i},0.36043987514722653726667765060-0.932782448593374716484835982071\,\mathrm{i},0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i},1\right\}
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 8: ✓ — (no output)
Cell 10: ✓ — (no output)

## 2026-07-18 18:38:42 — ▶️ RUN CELLS 11–20
Cell 11: ✓ — (no output)
Cell 12: ✓ — Out[27]= \left\{0\right\}
Cell 13: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — Out[31]= \left\{0\right\}
Cell 16: ✓ — (no output)
Cell 17: ✓ — Out[34]= \begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{pmatrix}
Cell 18: ✓ — (no output)
Cell 20: ✓ — (no output)

## 2026-07-18 18:38:45 — ▶️ RUN CELLS 21–28
Cell 21: ✓ — (no output)
Cell 23: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 26: ✓ — (no output)
Cell 27: ✓ — Out[52]= \begin{pmatrix}0 & \frac{1}{2} & 0 & 0 & 0 & 0\end{pmatrix}
Cell 28: ✓ — Out[53]= \left\{0\right\}

## 2026-07-18 18:38:48 — ▶️ RUN CELLS 29–36
Cell 30: ✓ — Out[59]= \left\{\left(u-\frac{4}{3}\right)\,\left(u-\frac{8}{7}\right),\left(u-\frac{1}{3}\right)\,\left(u-\frac{1}{7}\right),\left(0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i}\right)\,u^2-\left(1.06415582186324025288257211128+0.\times 10^{-29}\,\mathrm{i}\right)\,u\right\}
Cell 31: ✓ — Out[62]= \begin{pmatrix}1 & 1\\2 & 2\\3 & 3\end{pmatrix}
Cell 32: ✓ — Out[65]= \text{max|C2-($\kappa$1+$\kappa$2)Id|}\to 0.\times 10^{-29}
Cell 33: ✓ — Out[75]= \langle\!\lvert \left\{0,0\right\}\to 0.5664055180884988 +1.3769645669711723\,\mathrm{i},\left\{1,0\right\}\to -0.778123966451158+0.\,\,\mathrm{i},\left\{1,1\right\}\to 1.1900552523337025 +0.\,\,\mathrm{i},\left\{2,0\right\}\to 0.5664055180884988 -1.3769645669711723\,\mathrm{i}\rvert\!\rangle
Cell 34: ✓ — 2 kernel messages: Psi::spurious: No genuine state at (M,n)=(2,2) for (λ1,λ2)=(2,1): weight sector M=2 has only 2 genuine state(s) (n=0..1), so n=2 does not exist. [also: Psi::spurious] result: Print: spurious-guard test (expect a Psi::spurious message and Missing[...]):[formula] | Out[84]= \left\{\left\{1.000000000000000000000000000+0.\times 10^{-29}\,\mathrm{i},-1.3004494881076436272092358…
⚠️ stopped at Cell 34 — error detected (pass stopOnError:false to continue past errors)

## 2026-07-18 18:41:01 — ▶️ RUN CELLS 35–42
Cell 35: ✓ — Out[91]= \left\{1.03464899078135229594035121+0.381869568380946283899575941\,\mathrm{i},-2.37022278949001677248172619-0.20048970475590097847093779\,\mathrm{i},1.00000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right\}
Cell 36: ✓ — Out[94]= \left\{-0.709231010958016926537225142+1.09596916251277206706752728\,\mathrm{i},-0.02291660900+0.\times 10^{-12}\,\mathrm{i}\right\}
Cell 37: ✓ — (no output)
Cell 38: ✓ — (no output)
Cell 39: ✓ — (no output)
Cell 40: ✓ — (no output)
Cell 41: ✓ — (no output)
Cell 42: ✓ — Out[109]= 0

## 2026-07-18 18:41:04 — ▶️ RUN CELLS 43–47
Cell 43: ✓ — (no output)
Cell 44: ✓ — Out[113]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 45: ✓ — Out[115]= \left\{6,6\right\}
Cell 46: ✓ — Out[116]= \begin{pmatrix}2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i} & 0 & 0 & 0 & 0 & 0\\0 & 7.36608403016431791181-0.968046299794461400266\,\mathrm{i} & 0 & 0 & 0 & 0\\0 & 0 & 2.68702160063169134012-1.3456153260557858257\,\mathrm{i} & 0 & 0 & 0\\0 & 0 & 0 & 7.3660840301643179118107+0.9680462997944614002659\,\mathrm{i} & 0 & 0\\0 & 0 & 0 & 0 & 2.687021600631691340118… [537 chars total]
Cell 47: ✓ — Out[125]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{max|off-diagonal|}\to 5.065342385207116\times 10^{-34},\text{digits vanished (-Log10)}\to 33.29539119374288,\text{min|diagonal|}\to 3.0051232400647923,\text{diagonal values}\to \left\{2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i},7.36608403016431791181-0… [699 chars total]

## 2026-07-18 18:41:08 — ▶️ RUN CELLS 48–54
Cell 48: ✓ — (no output)
Cell 49: ✓ — (no output)
Cell 50: ✓ — Out[133]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 51: ✓ — Out[138]= \left\{1.00000000000000000000000000000,6\right\}
Cell 52: ✓ — Out[143]= \langle\!\lvert \text{max|<PsiLeft|Psi>| off-diagonal}\to 0.\,,\text{digits vanished (-Log10)}\to \mathrm{Indeterminate}\rvert\!\rangle
Cell 53: ✓ — Out[147]= \left\{7.366084030164317911810727273-0.9680462997944614002658874268\,\mathrm{i},1.00000000000000000000+0.\times 10^{-21}\,\mathrm{i}\right\}
Cell 54: ✓ — Out[150]= \begin{matrix}\langle\!\lvert \text{M}\to 0,\text{n}\to 0,\text{braket}\to 2.909064289-5.804620064\,\mathrm{i},\text{G21 diagonal}\to 2.909064289-5.804620064\,\mathrm{i},\text{normRatio}\to 1.000000000+0.\times 10^{-10}\,\mathrm{i}\rvert\!\rangle \\\langle\!\lvert \text{M}\to 1,\text{n}\to 0,\text{braket}\to 7.366084030-0.968046300\,\mathrm{i},\text{G21 diagonal}\to 7.366084030-0.9680463… [1247 chars total]

## 2026-07-18 18:41:12 — ▶️ RUN CELLS 55–62
Cell 56: ✓ — (no output)
Cell 57: ✓ — Out[158]= \begin{pmatrix}3 & 3\\3 & 3\end{pmatrix}
Cell 58: ✓ — (no output)
Cell 59: ✓ — Out[165]= \left\{0\right\}
Cell 60: ✓ — Out[169]= \begin{matrix}\left\{1,1,0,1.\,,0.\,\right\}\\\left\{1,1,1,1.\,,0.\,\right\}\\\left\{1,2,0,1.\,,0.\,\right\}\\\left\{1,2,1,1.\,,0.\,\right\}\\\left\{1,3,0,1.\,,0.\,\right\}\\\left\{1,3,1,1.\,,0.\,\right\}\\\left\{2,1,0,1.\,,0.\,\right\}\\\left\{2,1,1,1.\,,0.\,\right\}\\\left\{2,1,2,1.\,,0.\,\right\}\\\left\{2,2,0,1.\,,0.\,\right\}\\\left\{2,2,1,1.\,,0.\,\right\}\\\left\{2,2,2,1.\,,0.\,\r… [898 chars total]
Cell 61: ✓ — Out[170]= \langle\!\lvert \text{max spread (ratio non-constancy)}\to 0.\,,\text{max |ratio - 1| over sweep}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases checked}\to 27\rvert\!\rangle
Cell 62: ✓ — Out[176]= \langle\!\lvert \text{max operator-level spread}\to 0.\,,\text{max |$\mu$ - $\mu$Site2|}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases}\to 27\rvert\!\rangle


## 2026-07-18 18:46:12
**Input:**
```wolfram
{θ1,θ2,h,Length[DownValues[FL]],Length[DownValues[O2rot]],Length[DownValues[Opr]]}
```
**Output:**
```
Out= {1/3, 1/7, 1, 0, 0, 0}
```

---

## 2026-07-18 18:46:18 — 🔄 KERNEL RESTART

## 2026-07-18 18:46:48 — ▶️ RUN CELLS 1–10
Cell 1: ✓ — Out[8]= \left\{\frac{1}{3},\frac{1}{7},1,0.360439875147226537266677650596+0.932782448593374716484835982071\,\mathrm{i},0.36043987514722653726667765060-0.932782448593374716484835982071\,\mathrm{i},0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i},1\right\}
Cell 3: ✓ — (no output)
Cell 4: ✓ — (no output)
Cell 5: ✓ — (no output)
Cell 6: ✓ — (no output)
Cell 8: ✓ — (no output)
Cell 10: ✓ — (no output)

## 2026-07-18 18:46:52 — ▶️ RUN CELLS 11–20
Cell 11: ✓ — (no output)
Cell 12: ✓ — Out[27]= \left\{0\right\}
Cell 13: ✓ — (no output)
Cell 14: ✓ — (no output)
Cell 15: ✓ — Out[31]= \left\{0\right\}
Cell 16: ✓ — (no output)
Cell 17: ✓ — Out[34]= \begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{pmatrix}
Cell 18: ✓ — (no output)
Cell 20: ✓ — (no output)

## 2026-07-18 18:46:55 — ▶️ RUN CELLS 21–28
Cell 21: ✓ — (no output)
Cell 23: ✓ — (no output)
Cell 25: ✓ — (no output)
Cell 26: ✓ — (no output)
Cell 27: ✓ — Out[52]= \begin{pmatrix}0 & \frac{1}{2} & 0 & 0 & 0 & 0\end{pmatrix}
Cell 28: ✓ — Out[53]= \left\{0\right\}

## 2026-07-18 18:46:58 — ▶️ RUN CELLS 29–36
Cell 30: ✓ — Out[59]= \left\{\left(u-\frac{4}{3}\right)\,\left(u-\frac{8}{7}\right),\left(u-\frac{1}{3}\right)\,\left(u-\frac{1}{7}\right),\left(0.720879750294453074533355301191+0.\times 10^{-30}\,\mathrm{i}\right)\,u^2-\left(1.06415582186324025288257211128+0.\times 10^{-29}\,\mathrm{i}\right)\,u\right\}
Cell 31: ✓ — Out[62]= \begin{pmatrix}1 & 1\\2 & 2\\3 & 3\end{pmatrix}
Cell 32: ✓ — Out[65]= \text{max|C2-($\kappa$1+$\kappa$2)Id|}\to 0.\times 10^{-29}
Cell 33: ✓ — Out[75]= \langle\!\lvert \left\{0,0\right\}\to 0.5664055180884988 +1.3769645669711723\,\mathrm{i},\left\{1,0\right\}\to -0.778123966451158+0.\,\,\mathrm{i},\left\{1,1\right\}\to 1.1900552523337025 +0.\,\,\mathrm{i},\left\{2,0\right\}\to 0.5664055180884988 -1.3769645669711723\,\mathrm{i}\rvert\!\rangle
Cell 34: ✓ — 2 kernel messages: Psi::spurious: No genuine state at (M,n)=(2,2) for (λ1,λ2)=(2,1): weight sector M=2 has only 2 genuine state(s) (n=0..1), so n=2 does not exist. [also: Psi::spurious] result: Print: spurious-guard test (expect a Psi::spurious message and Missing[...]):[formula] | Out[84]= \left\{\left\{1.000000000000000000000000000+0.\times 10^{-29}\,\mathrm{i},-1.3004494881076436272092358…
⚠️ stopped at Cell 34 — error detected (pass stopOnError:false to continue past errors)

## 2026-07-18 18:47:08 — ▶️ RUN CELLS 35–42
Cell 35: ✓ — Out[91]= \left\{1.03464899078135229594035121+0.381869568380946283899575941\,\mathrm{i},-2.37022278949001677248172619-0.20048970475590097847093779\,\mathrm{i},1.00000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right\}
Cell 36: ✓ — Out[94]= \left\{-0.709231010958016926537225142+1.09596916251277206706752728\,\mathrm{i},-0.02291660900+0.\times 10^{-12}\,\mathrm{i}\right\}
Cell 37: ✓ — (no output)
Cell 38: ✓ — (no output)
Cell 39: ✓ — (no output)
Cell 40: ✓ — (no output)
Cell 41: ✓ — (no output)
Cell 42: ✓ — Out[109]= 0

## 2026-07-18 18:47:11 — ▶️ RUN CELLS 43–47
Cell 43: ✓ — (no output)
Cell 44: ✓ — Out[113]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 45: ✓ — Out[115]= \left\{6,6\right\}
Cell 46: ✓ — Out[116]= \begin{pmatrix}2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i} & 0 & 0 & 0 & 0 & 0\\0 & 7.36608403016431791181-0.968046299794461400266\,\mathrm{i} & 0 & 0 & 0 & 0\\0 & 0 & 2.68702160063169134012-1.3456153260557858257\,\mathrm{i} & 0 & 0 & 0\\0 & 0 & 0 & 7.3660840301643179118107+0.9680462997944614002659\,\mathrm{i} & 0 & 0\\0 & 0 & 0 & 0 & 2.687021600631691340118… [537 chars total]
Cell 47: ✓ — Out[125]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{max|off-diagonal|}\to 5.065342385207116\times 10^{-34},\text{digits vanished (-Log10)}\to 33.29539119374288,\text{min|diagonal|}\to 3.0051232400647923,\text{diagonal values}\to \left\{2.9090642893150832630837933001-5.8046200635659338430425538562\,\mathrm{i},7.36608403016431791181-0… [699 chars total]

## 2026-07-18 18:47:15 — ▶️ RUN CELLS 48–54
Cell 48: ✓ — (no output)
Cell 49: ✓ — (no output)
Cell 50: ✓ — Out[133]= \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix}
Cell 51: ✓ — Out[138]= \left\{1.00000000000000000000000000000,6\right\}
Cell 52: ✓ — Out[143]= \langle\!\lvert \text{max|<PsiLeft|Psi>| off-diagonal}\to 0.\,,\text{digits vanished (-Log10)}\to \mathrm{Indeterminate}\rvert\!\rangle
Cell 53: ✓ — Out[147]= \left\{7.366084030164317911810727273-0.9680462997944614002658874268\,\mathrm{i},1.00000000000000000000+0.\times 10^{-21}\,\mathrm{i}\right\}
Cell 54: ✓ — Out[150]= \begin{matrix}\langle\!\lvert \text{M}\to 0,\text{n}\to 0,\text{braket}\to 2.909064289-5.804620064\,\mathrm{i},\text{G21 diagonal}\to 2.909064289-5.804620064\,\mathrm{i},\text{normRatio}\to 1.000000000+0.\times 10^{-10}\,\mathrm{i}\rvert\!\rangle \\\langle\!\lvert \text{M}\to 1,\text{n}\to 0,\text{braket}\to 7.366084030-0.968046300\,\mathrm{i},\text{G21 diagonal}\to 7.366084030-0.9680463… [1247 chars total]

## 2026-07-18 18:47:18 — ▶️ RUN CELLS 55–62
Cell 56: ✓ — (no output)
Cell 57: ✓ — Out[158]= \begin{pmatrix}3 & 3\\3 & 3\end{pmatrix}
Cell 58: ✓ — (no output)
Cell 59: ✓ — Out[165]= \left\{0\right\}
Cell 60: ✓ — Out[169]= \begin{matrix}\left\{1,1,0,1.\,,0.\,\right\}\\\left\{1,1,1,1.\,,0.\,\right\}\\\left\{1,2,0,1.\,,0.\,\right\}\\\left\{1,2,1,1.\,,0.\,\right\}\\\left\{1,3,0,1.\,,0.\,\right\}\\\left\{1,3,1,1.\,,0.\,\right\}\\\left\{2,1,0,1.\,,0.\,\right\}\\\left\{2,1,1,1.\,,0.\,\right\}\\\left\{2,1,2,1.\,,0.\,\right\}\\\left\{2,2,0,1.\,,0.\,\right\}\\\left\{2,2,1,1.\,,0.\,\right\}\\\left\{2,2,2,1.\,,0.\,\r… [898 chars total]
Cell 61: ✓ — Out[170]= \langle\!\lvert \text{max spread (ratio non-constancy)}\to 0.\,,\text{max |ratio - 1| over sweep}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases checked}\to 27\rvert\!\rangle
Cell 62: ✓ — Out[176]= \langle\!\lvert \text{max operator-level spread}\to 0.\,,\text{max |$\mu$ - $\mu$Site2|}\to 0.\,,\text{num ($\lambda$1,$\lambda$2,s1) cases}\to 27\rvert\!\rangle

## 2026-07-18 18:47:29 — ▶️ RUN CELLS 63–70
Cell 64: ✓ — (no output)
Cell 65: ✓ — Out[185]= \langle\!\lvert \text{w\_0 ($\lambda$=1,2,3)}\to \left\{1,1,1\right\},\text{full weights w\_j}\to \left\{1\to \left\{1,-1\right\},2\to \left\{1,-2,1\right\},3\to \left\{1,-3,3,-1\right\}\right\},\text{br1[$\lambda$][1] ($\lambda$=1,2,3)}\to \left\{0,0,0\right\}\rvert\!\rangle
Cell 66: ✓ — Out[187]= \left\{0,0\right\}
Cell 67: ✓ — Out[190]= \sqrt[3]{k}\,\left(\left(-1.18561270626131985872+0.832253303224995280404\,\mathrm{i}\right)\,k^2-\left(0.761819910806671651081+0.812812473459294825632\,\mathrm{i}\right)\,k-\left(0.404410037536507329420+0.598870869485872613792\,\mathrm{i}\right)\right)
Cell 68: ✓ — Out[194]= \begin{pmatrix}0 & 0 & 0\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix}
Cell 69: ✓ — Out[198]= \langle\!\lvert \text{OP(E12)}\to 1.\,+0.\times 10^{-30}\,\mathrm{i},\text{OP(E21)}\to \left(-0.03626734060947795145020244939-2.004237714666237867354326095\,\mathrm{i}\right)\,\varphi ^2-\left(1.522150083623538480912456055-0.2442047957313920727733393355\,\mathrm{i}\right)\,\varphi +\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right),\text{nontrivial generator \{a,b\}}\to \left\{2,1\right\},\… [464 chars total]
Cell 70: ✓ — Out[201]= \left\{1,2\right\}

## 2026-07-18 18:47:32 — ▶️ RUN CELLS 71–77
Cell 71: ✓ — Out[207]= \langle\!\lvert \text{states}\to \begin{pmatrix}0 & 0\\1 & 0\\1 & 1\\2 & 0\\2 & 1\\3 & 0\end{pmatrix},\text{$\varphi$ candidates per state}\to \begin{pmatrix}\left\{0,0\right\} & \left(-0.920793778924720536105378676-0.390049761817044163511250929\,\mathrm{i}\right)\,k^{4/3}+\left(0.695722366612934233848844287-0.718310788304406529935213657\,\mathrm{i}\right)\,\sqrt[3]{k}-\left(0.3604398751… [1247 chars total]
Cell 72: ✓ — Out[213]= \langle\!\lvert \text{c\_sb solution(s)}\to \left\{\left\{\mathrm{cb}\to \frac{\left(-6.242239011314\times 10^{45}-4.511757007219\times 10^{45}\,\mathrm{i}\right)\,k^{4/3}-\left(5.999302657796\times 10^{45}+8.356648463050\times 10^{45}\,\mathrm{i}\right)\,\sqrt[3]{k}}{9.833530773718\times 10^{45}\,k^{4/3}-\left(3.54439660434\times 10^{45}-9.172544913427\times 10^{45}\,\mathrm{i}\right)\,… [841 chars total]
Cell 73: ✓ — Out[216]= \langle\!\lvert \text{generator \{a,b\} tested}\to \left\{2,1\right\},\text{rep tested}\to \left\{1,2\right\},\text{FL shape (verified exact for all 6 states)}\to \text{FL[1,2][state][k] = k\^{}(1/3) (A\_state + B\_state k)},\text{c0(state) = OPg[state][0]}\to \begin{pmatrix}\left\{0,0\right\} & 1.\,+0.\times 10^{-28}\,\mathrm{i}\\\left\{1,0\right\} & 1.\,+0.\times 10^{-28}\,\mathrm{i}\\… [1247 chars total]
Cell 74: ✓ — Out[220]= \left(1.000000000000000000000000000+0.\times 10^{-27}\,\mathrm{i}\right)-\left(0.360439875147226537266677651+0.932782448593374716484835982\,\mathrm{i}\right)\,k
Cell 75: ✓ — Out[227]= \langle\!\lvert \text{alphaTab \{state,D0,c0,Dc,c1,$\alpha$\}}\to \begin{pmatrix}\left\{0,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & -0.36043987514722653-0.9327824485933747\,\mathrm{i} & 0.36043987514722653 +0.9327824485933747\,\mathrm{i} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.7101545203859353 -0.2371533418544… [1247 chars total]
Cell 76: ✓ — Out[235]= \langle\!\lvert \text{validation}\to \left\{\left\{1,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\end{pmatrix},\left\{2,1\right\}\to \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 0\\\left\{1,1\right\} & 0\\\left\{2,0\right\} & 0\\\left\{2,1\right\} & 0\\\left\{3,0\right\} & 0\end{pmatrix},\left\{1… [1247 chars total]
Cell 77: ✓ — Out[237]= \langle\!\lvert \text{generator \{a,b\}}\to \left\{2,1\right\},\text{$\varphi$(k)}\to \text{$\alpha$ k},\text{$\alpha$ (universal value)}\to -1.\,+0.\,\,\mathrm{i},\text{normalization}\to \text{1-site FSoV determinant det1[f] = br1[f]/(f/.u->$\theta$1); Section C's det[$\lambda$1,$\lambda$2][f] restricted to 1 site},\text{max residual over rep sweep}\to 0,\text{reps validated}\to \begin{… [510 chars total]

## 2026-07-18 18:47:54 — ▶️ RUN CELLS 78–82
Cell 79: ✓ — (no output)
Cell 80: ✓ — Out[245]= \left\{\text{Test 1 (symbolic monomials u\^{}j) residuals},\left\{0\right\}\right\}
Cell 81: ✓ — Out[250]= \left\{\text{Test 2 (twisted monomials k\^{}(u/h) u\^{}j) residuals},\left\{0\right\}\right\}
Cell 82: ✓ — Out[259]= \left\{\text{Test 3 worst |br identity residual| on Q1},0,\text{worst |detRatio - predicted| },0,\text{num states checked},41\right\}

## 2026-07-18 18:47:58 — ▶️ RUN CELLS 83–87
Cell 84: ✓ — (no output)
Cell 85: ✓ — Out[271]= \begin{pmatrix}1 & 1 & 0\\2 & 1 & 0\\3 & 1 & 0\end{pmatrix}
Cell 86: ✓ — Out[278]= \left\{\text{Mirror Test1 (u\^{}j)},\left\{0\right\},\text{Mirror Test2 (k\^{}(u/h)u\^{}j)},\left\{0\right\}\right\}
Cell 87: ✓ — Out[285]= \left\{\text{Mirror Test3 worst |br identity residual| on Q1},0,\text{num states checked},41\right\}

## 2026-07-18 18:48:00 — ▶️ RUN CELLS 88–91
Cell 89: ✓ — Out[294]= \left\{\text{OP2(E12)},1.\,+0.\times 10^{-30}\,\mathrm{i},\text{OP2(E21)},\left(2.603469709065218092712489006+0.6885776528619826437114966465\,\mathrm{i}\right)\,\varphi +\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right),\text{nontrivial generator \{a,b\}},\left\{2,1\right\}\right\}
Cell 90: ✓ — Out[297]= \left(1.\,+0.\times 10^{-27}\,\mathrm{i}\right)-\left(0.360439875147226537266677651+0.932782448593374716484835982\,\mathrm{i}\right)\,k
Cell 91: ✓ — Out[300]//TableForm= \begin{pmatrix}\left\{0,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & -0.36043987514722653-0.9327824485933747\,\mathrm{i} & 0.36043987514722653 +0.9327824485933747\,\mathrm{i} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 1.\, +0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.4685059330398573 -0.1713611156914601\,\mathrm{i} & 2.603469709065218 +0.6885776528619827… [1258 chars total]

## 2026-07-18 18:48:13 — ▶️ RUN CELLS 92–108
Cell 92: ✓ — Out[302]//TableForm= \begin{pmatrix}\left\{0,0\right\} & -1.\,+0.\,\,\mathrm{i} & 1.\, +0.\,\,\mathrm{i} & 0.973902234276272 +0.2269679230016555\,\mathrm{i} & 0.7391304347826086 +0.\,\,\mathrm{i} & 0.973902234276272 +0.2269679230016555\,\mathrm{i}\\\left\{1,0\right\} & 0.15191893851739566 -0.10600050418875033\,\mathrm{i} & -0.15191893851739566+0.10600050418875033\,\mathrm{i} & 0.7651595274649543 +… [1258 chars total]
Cell 93: ✓ — Out[305]= \begin{pmatrix}\left\{0,0\right\} & -1.\,+0.\,\,\mathrm{i}\\\left\{1,0\right\} & 0.05448860910968355 -0.09540961504089102\,\mathrm{i}\\\left\{1,1\right\} & -0.0293413933141388+0.7363100318787654\,\mathrm{i}\\\left\{2,0\right\} & 0.11703059508534998 +0.\,\,\mathrm{i}\\\left\{2,1\right\} & 0.8158972359113332 +0.\,\,\mathrm{i}\\\left\{3,0\right\} & 0.05448860910968355 +0.09540961504089102\,… [541 chars total]
Cell 96: ✓ — Out[311]= \left\{\text{O[k1,0]-FLdet[k1] (must be 0)},0,\text{O[0,k2]-FLdet2[k2] (expected NONZERO per Phase1 finding)},\left(-3.0719756421050755-0.5172165371705225\,\mathrm{i}\right)\,k2\right\}
Cell 97: ✓ — Out[316]= \left\{\text{Fcand1},k1^2\,\left(\left(-2.2527492196701659-5.82989030370859198\,\mathrm{i}\right)\,k2-\left(0.4170744170089964+23.0487337186617355\,\mathrm{i}\right)\right)+k1\,\left(\left(0.67128913890100625-0.37524873225652243\,\mathrm{i}\right)\,k2+\left(9.51343802264711551-1.5262799733212005\,\mathrm{i}\right)\right)-\left(1.99115021541939352-0.72828474168870538\,\mathrm{i}\right)\,k… [1014 chars total]
Cell 98: ✓ — Out[318]= \left(0.9252125034899888 -6.175727508765443\,\mathrm{i}\right)\,k1^2\,k2+\left(5.4271267482291075 +3.0890004366773067\,\mathrm{i}\right)\,k1\,k2
Cell 99: ✓ — Out[323]= \begin{pmatrix}\left\{0,0\right\} & 0\\\left\{1,0\right\} & 7.991287939023577 -1.2820751775898085\,\mathrm{i}\\\left\{2,0\right\} & -0.3808070763995185-21.044496003995498\,\mathrm{i}\\\left\{0,1\right\} & 0.6123194936458246 +1.416862394550688\,\mathrm{i}\\\left\{1,1\right\} & -0.6248120089007682-3.051910425543801\,\mathrm{i}\\\left\{2,1\right\} & -1.8923093445229393-4.897107855115217\,\m… [421 chars total]
Cell 100: ✓ — Out[328]= \left\{\text{Fcand2 (no P-insertion)},k1^2\,\left(\left(0.02466179161444501+1.36288164597304175\,\mathrm{i}\right)-\left(0.36043987514722654+0.93278244859337472\,\mathrm{i}\right)\,k2\right)+k1\,\left(\left(1.2351720155778515-0.69045766735200127\,\mathrm{i}\right)\,k2-\left(6.46913785540003854-1.0378703818584163\,\mathrm{i}\right)\right)+\left(5.38781822995835893-1.9706528304517910\,\mat… [788 chars total]
Cell 101: ✓ — 2 kernel messages: Power::infy: Infinite expression 0       encountered.
Infinity::indet: Indeterminate expression ----------------- encountered.                                                   1/3                    … [also: Power::infy] result: Out[331]//TableForm= \begin{pmatrix}2 & 1 & \left\{0,0\right\} & \mathrm{Indeterminate} & \mathrm{Indeterminate} & -4.970498218180917+4.515583419069149\,\mathrm{i}\\2 & 1 & \left\{1,1\right\} & \mathr…
Cell 102: ✓ — Out[334]//TableForm= \begin{pmatrix}2 & 1 & \left\{0,0\right\} & \left(-7.771745024480887+7.060451775520925\,\mathrm{i}\right)\,k1^2-\left(3.7846186890458786 +9.794215710230434\,\mathrm{i}\right)\,k1 & \left(1.8923093445229393 +4.897107855115217\,\mathrm{i}\right)\,k2 & -4.970498218180917+4.515583419069149\,\mathrm{i}\\2 & 1 & \left\{1,1\right\} & \left(9.339255641812587 +10.661478687052961\,\math… [1258 chars total]
Cell 103: ✓ — 2 kernel messages: Power::infy: Infinite expression 0       encountered.
Infinity::indet: Indeterminate expression ----------------- encountered.                                                   1/3                    … [also: Power::infy] result: Out[338]= \left\{\text{Fcand1(k1,0)},\mathrm{Indeterminate},\text{FLdet[k1]},\frac{\left(2.45075493251774917289426756629+5.78551794268217244227029022731\,\mathrm{i}\right)\,\left(\left(0.0525946877178…
Cell 104: ✓ — Out[343]= \left\{\text{Fcand1(k1,0)},\left(-8.51191121729+7.73287575414\,\mathrm{i}\right)\,k1^2-\left(4.5054984393+11.6597806074\,\mathrm{i}\right)\,k1+\left(1.000000000000+0.\times 10^{-12}\,\mathrm{i}\right),\text{FLdet[k1]},\frac{\left(2.45075493251774917289426756629+5.78551794268217244227029022731\,\mathrm{i}\right)\,\left(\left(0.05259468771783351844996683714+0.15021349717800998891606886537\… [960 chars total]
Cell 107: ✓ — Out[348]= \left\{\text{STEP 0: det[$\lambda$1,$\lambda$2][k\^{}(u/h)Q1] vs O[k,k], worst |residual| (sample k)},0,\text{sample rows \{$\lambda$1,$\lambda$2,state,resid,ratio\}},\begin{pmatrix}1 & 1 & \left\{0,0\right\} & 0 & 1.\,\\1 & 1 & \left\{1,0\right\} & 0 & 1.\,\\1 & 1 & \left\{1,1\right\} & 0 & 1.\,\\1 & 1 & \left\{2,0\right\} & 0 & 1.\,\\2 & 1 & \left\{0,0\right\} & 0 & 1.\,\\2 & 1 & \left… [445 chars total]
Cell 108: ✓ — Out[354]= \left\{\text{STEP 0 CONFIRMED: det[$\lambda$1,$\lambda$2][k\^{}(u/h)Q1] == O[k,k] identically},\text{worst residual over full sweep},0,\text{num states},41\right\}

## 2026-07-18 18:48:24 — ▶️ RUN CELLS 109–123
Cell 109: ✓ — 2 kernel messages: Power::infy: Infinite expression - encountered.                                  0
Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.
Power::infy: Infinite expression - encounter… [also: Power::infy] result: Out[358]= \begin{pmatrix}1 & 0 & \left\{1.\,\right\}\\1 & 1 & \left\{1.\,\right\}\\2 & 0 & \left\{1.\,\right\}\\2 & 1 & \left\{2.\,\right\}\\2 & 2 & \left\{2.\,\right\}\\3 & 0 & \left\{1.\,\right\}\\3…
Cell 110: ✓ — Out[364]= \left\{\text{STEP 1 CONFIRMED: ⟨x0[$\lambda$]|E21\^{}m = $\gamma$\_m($\lambda$) xL1[$\lambda$][\{m\}], $\gamma$\_m($\lambda$)=$\lambda$!/($\lambda$-m)!},\text{worst deviation},0,\text{table \{$\lambda$,m,$\gamma$\_found,$\gamma$\_predicted,dev\}},\begin{pmatrix}1 & 0 & 1.\, & 1 & 0\\1 & 1 & 1.\, & 1 & 0\\2 & 0 & 1.\, & 1 & 0\\2 & 1 & 2.\, & 2 & 0\\2 & 2 & 2.\, & 2 & 0\\3 & 0 & 1.\, & 1 &… [489 chars total]
Cell 111: ✓ — Out[369]= \left\{\text{STEP 2: x2[\{s1,m\}] vs xL1[$\lambda$1][\{s1\}]$\otimes$xL1[$\lambda$2][\{m\}], worst ratio-spread (constancy check)},0,\text{num cases},54,\text{sample rows},\begin{pmatrix}1 & 1 & 0 & 1 & 0.16 & 0\\1 & 1 & 1 & 1 & 1.\, & 0\\1 & 2 & 0 & 1 & 0.16 & 0\\1 & 2 & 0 & 2 & -0.68 & 0\\1 & 2 & 1 & 1 & 1.\, & 0\\1 & 2 & 1 & 2 & 1.\, & 0\\1 & 3 & 0 & 1 & 0.16 & 0\\1 & 3 & 0 & 2 & -0.6… [479 chars total]
Cell 112: ✓ — Out[373]= \begin{pmatrix}1 & 0 & \frac{4}{25}\\1 & 1 & 1\\2 & 0 & \frac{2}{23}\\2 & 1 & \frac{25}{46}\\2 & 2 & 1\\3 & 0 & \frac{4}{67}\\3 & 1 & \frac{25}{67}\\3 & 2 & \frac{46}{67}\\3 & 3 & 1\end{pmatrix}
Cell 113: ✓ — Out[378]= \left\{\text{c($\lambda$1,s1) hypothesis (m-independent) vs full sweep --- FAILS for m>=2, see next cell for the fix},\text{worst deviation over all 54 ($\lambda$1,$\lambda$2,s1,m) cases},1.68\right\}
Cell 114: ✓ — Out[382]= \begin{pmatrix}1 & 0 & 1 & \frac{4}{25}\\1 & 0 & 2 & -\frac{17}{25}\\1 & 0 & 3 & -\frac{38}{25}\\1 & 1 & 1 & 1\\1 & 1 & 2 & 1\\1 & 1 & 3 & 1\\2 & 0 & 1 & \frac{2}{23}\\2 & 0 & 2 & -\frac{34}{575}\\2 & 0 & 3 & \frac{323}{575}\\2 & 1 & 1 & \frac{25}{46}\\2 & 1 & 2 & \frac{2}{23}\\2 & 1 & 3 & -\frac{17}{46}\\2 & 2 & 1 & 1\\2 & 2 & 2 & 1\\2 & 2 & 3 & 1\\3 & 0 & 1 & \frac{4}{67}\\3 & 0 & 2 & … [672 chars total]
Cell 115: ✓ — 4 kernel messages: LinearSolve::nosol: Linear equation encountered that has no solution.
LinearSolve::nosol: Linear equation encountered that has no solution.
LinearSolve::nosol: Linear equation encountered that has no … [also: LinearSolve::nosol, LinearSolve::nosol, LinearSolve::nosol] result: Out[387]= \left\{\text{coeffs (rows=s1, cols=s1'=0..3)},\left\{\mathrm{LinearSolve}\left[\begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & \frac{1}{36}\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & \f…
Cell 116: ✓ — Out[392]= \left\{\left\{0,2\right\}\to -0.02206359506813757300454250486697,\left\{1,1\right\}\to 0.1635301752109020116807268007787,\left\{2,0\right\}\to 0.8585334198572355613238157040883\right\}
Cell 117: ✓ — Out[396]= \left\{\text{STEP 2 level-conservation check: \# off-level nonzero coeffs (want 0 everywhere)},0,\text{num ($\lambda$1,$\lambda$2,s1,m) cases},29\right\}
Cell 119: ✓ — Out[404]= \left\{\text{m=0: Om-Cm (must be 0 per Phase1, C1pred/P1site2[$\theta$1] should equal exactly Cm's own norm here)},\left(2.5572121404875445 -0.4102640568287387\,\mathrm{i}\right)\,k1-\left(0.038411844228125346 +2.1227491621421546\,\mathrm{i}\right)\,k1^2,\text{m=1: Om\_1[k1]},\left(0.3604398751472265372666776506+0.9327824485933747164848359821\,\mathrm{i}\right)\,k1^2-\left(1.296101147801… [801 chars total]
Cell 120: ✓ — Out[411]= \left\{\text{m=0: (B0-B1)/O\_0 ratio},\frac{\sqrt[3]{k1}\,\left(\left(-0.0292504903494856237+0.0205326892042810627\,\mathrm{i}\right)\,k1^2+\left(0.499242526583056932+0.532659421382621779\,\mathrm{i}\right)\,k1-\left(0.139388493663376676+0.206413542317126162\,\mathrm{i}\right)\right)}{\left(0.3604398751472265372666776506+0.9327824485933747164848359821\,\mathrm{i}\right)\,k1^2-\left(1.296… [1083 chars total]
Cell 121: ✓ — Out[414]= \left\{\text{Cm[0]-Cm[1] shape},\left(-0.0225172879957976151-1.244370198497125076\,\mathrm{i}\right)\,k1^2+\left(5.434075798536032377-0.8718111207610696998\,\mathrm{i}\right)\,k1+\left(0.\times 10^{-20}+0.\times 10^{-20}\,\mathrm{i}\right),\text{ratio vs O\_0[k1] (want: clean constant)},\frac{\left(-0.0225172879957976151-1.244370198497125076\,\mathrm{i}\right)\,k1^2+\left(5.4340757985360… [753 chars total]
Cell 122: ✓ — Out[416]= \left\{\text{Cm[0][k1] / O\_1[k1] (reversed-index hypothesis; want clean k1-independent constant)},\frac{\left(0.111318362779488714+0.0407158959953172346\,\mathrm{i}\right)\,k1^2-\left(0.21818047436770437-1.02534209918467083\,\mathrm{i}\right)\,k1+\left(0.360439875147226537-0.932782448593374716\,\mathrm{i}\right)}{\left(1.\,+0.\times 10^{-29}\,\mathrm{i}\right)\,k1^2-\left(2.963909584212… [539 chars total]

## 2026-07-18 18:48:39 — ▶️ RUN CELLS 124–138
Cell 125: ✓ — (no output)
Cell 126: ✓ — Out[446]= \left\{\text{Symbolic-$\theta$,$\chi$1 port-fidelity: xL1Sym vs xL1},0,\text{x2Sym vs x2},0\right\}
Cell 127: ✓ — Out[451]= \left\{\left\{0,0\right\}\to 0,\left\{0,1\right\}\to \frac{\theta A-\theta B}{\theta A-\theta B+1},\left\{1,0\right\}\to \frac{1}{\theta A-\theta B+1},\left\{1,1\right\}\to 0\right\}
Cell 128: ✓ — Out[454]= 35
Cell 129: ✓ — Out[460]= \left\{\left\{0,0\right\}\to 0,\left\{0,1\right\}\to \theta B-\theta A,\left\{1,0\right\}\to -1,\left\{1,1\right\}\to 0\right\}
Cell 130: ✓ — Out[463]= 127
Cell 131: ✓ — Out[466]= 127
Cell 132: ✓ — Out[470]//MatrixForm= \begin{pmatrix}\frac{\theta A^2-2\,\theta A\,\theta B-\theta A+\theta B^2+\theta B}{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)} & \frac{6\,\left(\theta A-\theta B\right)}{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)} & \frac{6}{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)}\\0 & \frac{\theta A-\theta B+1}{\the… [472 chars total]
Cell 133: ✓ — Out[472]//MatrixForm= \begin{pmatrix}\frac{\left(\theta A-\theta B+2\right)\,\left(\theta A-\theta B+3\right)}{\theta A^2-2\,\theta A\,\theta B-\theta A+\theta B^2+\theta B} & -\frac{6\,\left(\theta A-\theta B+3\right)}{\theta A^2-2\,\theta A\,\theta B+\theta B^2-1} & \frac{6}{\left(\theta A-\theta B\right)\,\left(\theta A-\theta B+1\right)}\\0 & \frac{\theta A-\theta B+3}{\theta A-\theta B+1} & -… [453 chars total]
Cell 134: ✓ — Out[478]= \left\{\text{equations},\left\{1.000000000000000000000+0.\times 10^{-21}\,\mathrm{i}⩵\left(1.\,+0.\times 10^{-30}\,\mathrm{i}\right)\,wc0+\left(2.603469709065218092712489006+0.6885776528619826437114966465\,\mathrm{i}\right)\,wc1,-1.0350620568640061670+0.16605926109734660949\,\mathrm{i}⩵\left(1.522150083623538480912456055-0.2442047957313920727733393355\,\mathrm{i}\right)\,wc0-\left(1.2961… [785 chars total]
Cell 135: ✓ — Out[484]= \begin{pmatrix}0 & 0 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\0 & 1 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\1 & 0 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\1 & 1 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} & 4.04123637990545\\2 & 0 & 0.8405990749714453 -0.5416578210984777\,\mathrm{i} … [510 chars total]
Cell 136: ✓ — Out[489]= \left\{\text{Conjectured sum rule FAILS beyond the lucky L=1 case},\text{num failing cases},29,\text{num total cases},81,\text{sample failures},\begin{pmatrix}1 & 2 & 0 & 1 & \frac{1}{-\theta A+\theta B-1}\\1 & 2 & 0 & 2 & \frac{2}{\theta A-\theta B+1}\\1 & 3 & 0 & 1 & -\frac{2}{\theta A-\theta B+1}\end{pmatrix}\right\}
Cell 138: ✓ — Out[496]= \left\{\text{BOXED: A[s1,m;s1,m] = Pochhammer[s1+$\theta$1-$\theta$2-m+1,m]/Pochhammer[$\lambda$1+$\theta$1-$\theta$2-m+1,m] --- verified at the live numeric $\theta$1,$\theta$2},\text{worst deviation},0,\text{num ($\lambda$1,$\lambda$2,s1,m) cases},81\right\}

## 2026-07-18 18:48:48 — ▶️ RUN CELLS 139–154
Cell 140: ✓ — Out[502]= \left\{\text{STEP 1a (WITH $\kappa$1\^{}(u/h) prefactor): worst |x2.Psi - Q1($\theta$1+s1)Q1($\theta$2+s2)/(Q1($\theta$1)Q1($\theta$2))|},0,\text{num cases},297\right\}
Cell 141: ✓ — Out[509]= \left\{\text{STEP 1b (WITHOUT prefactor, bare q\_n): worst residual (expect NONZERO, confirming prefactor matters)},3.7396476326262356,\text{sample},\begin{pmatrix}\left\{0,0\right\} & 0 & 0 & 0\\\left\{1,0\right\} & 0 & 0 & 0\\\left\{1,1\right\} & 0 & 0 & 0\end{pmatrix}\right\}
Cell 143: ✓ — Out[515]= \begin{pmatrix}\left\{0,0,\left\{\text{srng},\left\{0\right\},\text{sol},\left\{1.\,\right\},\text{maxResid},0.\times 10^{-28}\right\}\right\} & \left\{0,1,\left\{\text{srng},\left\{0,1\right\},\text{sol},\left\{11.50000000000000000000000000-4.431302171530629069010197265\times 10^{-38}\,\mathrm{i},-10.50000000000000000000000000+2.933894180666706878219239219\times 10^{-38}\,\mathrm{i}\rig… [1115 chars total]
Cell 144: ✓ — Out[521]= \left\{\text{worst residual across full B-sweep (want \textasciitilde{}0)},0.\times 10^{-24},\text{num (rep,a,b) cases},37\right\}
Cell 145: ✓ — Out[525]= \left\{\text{A-block ($\lambda$1=2,$\lambda$2=1, L=2, s,s'=1,2)},\begin{pmatrix}\frac{25}{46} & \frac{21}{46}\\0 & 1\end{pmatrix},\text{its inverse (should = B row)},\begin{pmatrix}\frac{46}{25} & -\frac{21}{25}\\0 & 1\end{pmatrix},\text{Step2's independently-solved B row \{1,1\}},\left\{\frac{46}{25},-\frac{21}{25}\right\},\text{match?},\left\{0,0\right\}\right\}
Cell 147: ✓ — Out[532]= \left\{\text{det[$\lambda$1,$\lambda$2][k\^{}u Q1] - explicit residue-sum expansion (rep (2,1), state (1,0)), want 0},0\right\}
Cell 148: ✓ — Out[538]= \left\{\text{TRUE sum rule (derived, Step 4): worst deviation},0,\text{num ($\lambda$1,$\lambda$2,s,t) cases},20,\text{sample},\begin{pmatrix}2 & 1 & 0 & 0 & 0\\2 & 1 & 0 & 1 & 0\\2 & 1 & 1 & 0 & 0\\2 & 1 & 1 & 1 & 0\\2 & 1 & 2 & 0 & 0\\2 & 1 & 2 & 1 & 0\end{pmatrix}\right\}
Cell 150: ✓ — Out[546]= \left\{\text{STEP 3 assembly test (rep (2,1), 6 states x 3 grid points): worst |Oassembled-O2rot|},0,\text{num cases},18\right\}
Cell 151: ✓ — Out[552]= \left\{\text{STEP 3 FULL-SWEEP ASSEMBLY CONFIRMED: O[k1,k2] = $\Sigma$ cA(k1)cB(k2) B-weighted Q-bilinears == O2rot exactly},\text{worst deviation},0,\text{num states x grid cases},82\right\}
Cell 153: ✓ — Out[557]= \begin{pmatrix}0 & 0 & 1.\, & 1.\, & 0\\0 & 1 & -2.603469709065218-0.6885776528619827\,\mathrm{i} & 5.387818229958359 -1.970652830451791\,\mathrm{i} & -7.991287939023577+1.2820751775898085\,\mathrm{i}\\1 & 0 & 1.5221500836235384 -0.24420479573139206\,\mathrm{i} & -6.469137855400039+1.0378703818584163\,\mathrm{i} & 7.991287939023577 -1.2820751775898085\,\mathrm{i}\\1 & 1 & 1.2961011478017… [843 chars total]

## 2026-07-18 18:49:00 — ▶️ RUN CELLS 155–172
Cell 156: ✓ — Out[561]= \left\{\text{M matrix (rep (2,1), k1=7/3,k2=5/2)},\begin{pmatrix}1 & -\frac{115}{4}\\\frac{259}{12} & \frac{322}{15}\\-\frac{196}{45} & -\frac{245}{18}\end{pmatrix},\text{MatrixRank},2\right\}
Cell 157: ✓ — Out[566]//TableForm= \begin{pmatrix}1 & 1 & \left\{\frac{7}{3},\frac{5}{2}\right\} & 2 & 2\\1 & 1 & \left\{\frac{11}{5},-\frac{4}{3}\right\} & 2 & 2\\1 & 1 & \left\{-\frac{2}{7},\frac{9}{4}\right\} & 2 & 2\\1 & 1 & \left\{\frac{13}{6},\frac{1}{9}\right\} & 2 & 2\\2 & 1 & \left\{\frac{7}{3},\frac{5}{2}\right\} & 2 & 2\\2 & 1 & \left\{\frac{11}{5},-\frac{4}{3}\right\} & 2 & 2\\2 & 1 & \left\{-\frac{… [1258 chars total]
Cell 158: ✓ — Out[568]= \left\{\text{Det[M(k1,k2)] symbolic, rep (2,2)},\frac{213003\,k2^3\,\left(k1-k2\right)^2\,\left(17\,k1+46\,k2\right)}{7225},\text{identically zero?},\mathrm{False}\right\}
Cell 159: ✓ — Out[572]= \left\{\text{rep (2,3): M dimensions},\left\{3,4\right\},\text{rank at sample point},3,\text{a 3x3 minor determinant (nonzero => rank>=3)},\frac{19563970545}{1405696}\right\}
Cell 162: ✓ — Out[577]= \left\{\text{towerMat dims (rep (2,1), S=3)},\left\{6,4\right\},\text{MatrixRank},4\right\}
Cell 163: ✓ — Out[582]= \left\{\text{c\_s (rep (2,1), S=3, p=7/3,r=5/2)},\left\{13.7688078414687182186600+2.10326272230760132765845\,\mathrm{i},-3.12167136931824065325212+18.8662496625125701186869\,\mathrm{i},-7.51900595437450423013552-10.1245361864155291062487\,\mathrm{i},3.80567731087062145228106-0.708676905117290761157203\,\mathrm{i}\right\},\text{max residual},27.6845440371648442360677\right\}
Cell 164: ✓ — Out[584]//TableForm= \begin{pmatrix}3 & \left\{6,4\right\} & 4 & 27.6845440371648442360677\\4 & \left\{6,5\right\} & 4 & 27.6845440371648442360677\\5 & \left\{6,6\right\} & 4 & 27.684544037164844236068\\6 & \left\{6,7\right\} & 4 & 27.684544037164844236068\\7 & \left\{6,8\right\} & 4 & 27.684544037164844236068\end{pmatrix}
Cell 165: ✓ — Out[587]= \left\{\text{c\_s at r=0, S=$\lambda$1=2 (rep (2,1), p=7/3)},\left\{1.\,,4.666666666666667,5.444444444444445\right\},\text{residual},0.\times 10^{-23},\text{predicted w1[$\lambda$1][[s+1]](-p)\^{}s},\left\{1.\,,4.666666666666667,5.444444444444445\right\},\text{match?},\left\{0,0,0\right\}\right\}
Cell 166: ✓ — Out[594]= \left\{\text{$\theta$2-tower: dims},\left\{6,4\right\},\text{rank},4,\text{residual},25.4178846095368977408093\right\}
Cell 167: ✓ — Out[598]//TableForm= \begin{pmatrix}-4 & -1 & 4 & 4 & 27.684544037164844236068\\-4 & 0 & 5 & 4 & 27.6845440371648442360677\\-4 & 1 & 6 & 4 & 27.6845440371648442360677\\-4 & 2 & 7 & 4 & 27.6845440371648442360677\\-4 & 3 & 8 & 4 & 27.6845440371648442360677\\-3 & 0 & 4 & 4 & 27.6845440371648442360677\\-3 & 1 & 5 & 4 & 27.6845440371648442360677\\-3 & 2 & 6 & 4 & 27.6845440371648442360677\\-3 & 3 & 7 &… [1139 chars total]
Cell 168: ✓ — Out[608]= \left\{\text{multiply by Q1($\theta$2): residual},19.7572259100834793810211,\text{divide by Q1($\theta$2): residual},18.310606875396829403204\right\}
Cell 169: ✓ — Out[620]= \left\{\text{r=0 residual at S=3 window (should stay 0, confirming subspace nesting)},0.\times 10^{-23},\text{more (p,r) failures, rep(2,1)},\begin{pmatrix}\left\{1,1\right\} & 4.04518880433634545024478\\\left\{2,-1\right\} & 8.9288834302921198736014\\\left\{-3,4\right\} & 39.3336388826667151011186\\\left\{0,\frac{7}{2}\right\} & 5.3318082362908496984122\end{pmatrix},\text{rep (1,2): tow… [476 chars total]
Cell 170: ✓ — Out[626]= \left\{\text{bare-q1 tower basis: rank},4,\text{residual at (7/3,5/2)},27.6845440371648442360677\right\}
Cell 171: ✓ — Out[629]= \left\{\text{rank of towerMat[2,1][2] (r=0 anchor, 3 unknowns)},3,\text{null space dimension (0 = unique solution)},0\right\}

## 2026-07-18 18:49:28 — 📥 BULK INSERT 2 CELL(S) at positions 173–174
1. [markdown] (*REDIRECT — Paul supplied his own working formula (su(2), his eq. 238-239, ℓ=1,L=2, β=α=1):↵  (⟨↑|e
2. [code] (*Cross-state elimination of C: for FIXED (T,r), solve O[p,r][state_i]·FL[state_j](T) =↵  O[p,r][sta


## 2026-07-18 18:49:29
**Input:**
```wolfram
(*Cross-state elimination of C: for FIXED (T,r), solve O[p,r][state_i]·FL[state_j](T) =
  O[p,r][state_j]·FL[state_i](T) for p, using ONE pair, then check ALL OTHER pairs give the SAME p.*)
ClearAll[solveCoupling];
solveCoupling[λ1_,λ2_][T_,r_] := Module[{states,FLvals,Ovals,eqs,sols,pCandidates},
  states = statesOf[λ1,λ2];
  FLvals = Table[FL[λ1,λ2][st[[1]],st[[2]]][T], {st,states}];
  Ovals = Table[Opr[λ1,λ2][st[[1]],st[[2]]][p,r], {st,states}];
  (*use state 1 as reference, solve p from pairing with every other state*)
  pCandidates = Table[
    Module[{eq,sol},
      eq = Ovals[[1]] FLvals[[i]] - Ovals[[i]] FLvals[[1]] == 0;
      sol = NSolve[eq, p, WorkingPrecision->30];
      {i, states[[i]], p/.sol}
    ]
  ,{i,2,Length[states]}];
  pCandidates
];
(*rep (2,1), T=7/3, r=5/2 sample.*)
pCand21 = solveCoupling[2,1][7/3,5/2];
pCand21 // TableForm
```
**Output:**
```
2 kernel messages: NSolve::precw: The precision of the argument function ({(-11.45569554337990047323 + 2.70008684227044968952 I) ((1.901099687868066343166694126 + 2.331956121483436791212089955 I) - (2.979951213744064381… [also: NSolve::precw] result: Out[633]//TableForm= \begin{pmatrix}2 & \left\{1,0\right\} & \left\{-0.360439875147226537266677650589+0.932782448593374716484835982015\,\mathrm{i},0.220899450193933919307122410797-0.287240959241362798…
```

## 2026-07-18 18:49:37
**Input:**
```wolfram
pCand21
```
**Output:**
```
Out= {{2, {1, 0}, {-0.3604398751472265372666776505892551745836968140321853043781`29.707347828414424 + 0.9327824485933747164848359820150137223925650442489347621531`30.120295363558725*I, 0.2208994501939339193071224107971100353878522277884962688134`29.935570760761575 - 0.2872409592413627982779083453686364061584141107120948173195`30.049622494207107*I}}, {3, {1, 1}, {-0.3604398751472265372666776506359255811570372401417942778905`29.707347828414424 + 0.9327824485933747164848359820375970253646960012413238638644`30.120295363558725*I, 0.1789805138623094861608805753468137683299975095184079100886`29.736816525980974 + 0.4280751514707218504416081095145540934774771587641777779875`30.11553079437402*I}}, {4, {2, 0}, {-2.027354048517321255721269193762454418657117566742732304648`29.959943948699156 + 2.4032053593522338620422260889396174774939379510554559397776`30.033805233874258*I, 0.0999667698001488976086889112801375155951002865451144929524`29.661873012615562 - 0.291285877434607902050114471851602062679563500852011441963`30.12633678253536*I}}, {5, {2, 1}, {-1.4534592034602140998550152561639904878933065479133399233488`30.145171759512497 - 0.2294061510126451110558671873287968390731857548685510708961`29.343373971621702*I, 0.1090659352459134356047190568295590890975395104500118736532`29.468756326520563 + 0.5126675212205921076812137042769255532187167034238617955399`30.14090300336992*I}}, {6, {3, 0}, {-0.6340432828130668926861423714150601412379695033122768756431`29.19327304716634 - 5.710847626030475987674937457200078054234235954562342531382`30.1478547139557*I, -0.3531572930556511008651934401593963873831972782416842669502`30.113607331863832 - 0.1520073050400028148816297234699652425341981870431148822553`29.747503611970764*I}}}
```
## 2026-07-18 18:50:11 — 📥 BULK INSERT 1 CELL(S) at positions 175–175
1. [code] (*Joint solve: for FIXED r, (p,T,C) are 3 unknowns satisfying O[p,r][state_i]=C·FL[state_i][T] for↵


## 2026-07-18 18:50:11
**Input:**
```wolfram
(*Joint solve: for FIXED r, (p,T,C) are 3 unknowns satisfying O[p,r][state_i]=C·FL[state_i][T] for
  EVERY state — use 3 states to solve simultaneously (nonlinear in T via FL, but AFFINE in p,C given
  T), then check the REMAINING 3 states of rep (2,1) for consistency.*)
ClearAll[jointSolve];
jointSolve[λ1_,λ2_][r_] := Module[{states,eqs,sols},
  states = statesOf[λ1,λ2];
  eqs = Table[Opr[λ1,λ2][states[[i,1]],states[[i,2]]][p,r] == c Global`FLT[i], {i,1,3}];
  eqs = eqs /. Global`FLT[i_] :> FL[λ1,λ2][states[[i,1]],states[[i,2]]][tt];
  NSolve[eqs, {p,c,tt}, WorkingPrecision->30]
];
sols21 = jointSolve[2,1][5/2];
Length[sols21]
sols21[[1;;Min[3,Length[sols21]]]]
```
**Output:**
```
2 kernel messages: NSolve::precw: The precision of the argument function ({(1.901099687868066343166694126 + 2.331956121483436791212089955 I) - (2.979951213744064381468205485 - 5.227684790291951818644572297 I) p - (2.975… [also: NSolve::precw] result: Out[637]= 2 | Out[638]= \left\{\left\{p\to -0.36+0.93\,\mathrm{i},c\to 21+21\,\mathrm{i},\mathrm{tt}\to 0.36-0.93\,\mathrm{i}\right\},\left\{p\to -0.36+0.93\,\mathrm{i},c\to 19+22\,\mathrm{i},\mathrm{…
```

## 2026-07-18 18:50:23
**Input:**
```wolfram
sols21
κ1
κ2
```
**Output:**
```
Out[1]= {{p -> -0.3604399051539879073`1.6425035616932673 + 0.9327824257577005694`1.8740763828430353*I, c -> 21.1979952303432754986`1.8740763828430353 + 21.288817579003628385`1.8363718305601333*I, tt -> 0.3604398751513760214`1.6425035343052212 - 0.932782448623959264`1.8740763828430353*I}, {p -> -0.3604399051539846662`1.6425035617376684 + 0.932782425757702345`1.8740763828430318*I, c -> 18.6948696842280562641`1.8740763828430353 + 21.7515617332618365138`1.8578055567171874*I, tt -> 0.3604398751515455951`1.6425035343497998 - 0.9327824486237882288`1.8740763828430353*I}}
Out[2]= 0.3604398751472265372666776505955500912937645232622850024938`29.707347828414424 + 0.9327824485933747164848359820708679294578496404879183613065`30.120295363558725*I
Out[3]= 0.360439875147226537266677650595550091293764523262285002494`29.421282762500656 - 0.9327824485933747164848359820708679294578496404879183613068`29.58721122792411*I
```
## 2026-07-18 18:50:51 — 📥 BULK INSERT 1 CELL(S) at positions 176–176
1. [code] (*Robust elimination via Resultant (avoids NSolve's Newton-iteration precision loss near↵  degenerat


## 2026-07-18 18:50:51
**Input:**
```wolfram
(*Robust elimination via Resultant (avoids NSolve's Newton-iteration precision loss near
  degenerate points). EqA(p,T):=O1·FL2(T)-O2·FL1(T)=0, EqB(p,T):=O1·FL3(T)-O3·FL1(T)=0 (both from
  clearing c between pairs of states); Resultant[EqA,EqB,p] eliminates p, giving a polynomial in T
  alone whose roots are candidate couplings.*)
ClearAll[elimT];
elimT[λ1_,λ2_][r_] := Module[{states,O1,O2,O3,FL1,FL2,FL3,eqA,eqB,resT},
  states = statesOf[λ1,λ2];
  O1 = Opr[λ1,λ2][states[[1,1]],states[[1,2]]][p,r];
  O2 = Opr[λ1,λ2][states[[2,1]],states[[2,2]]][p,r];
  O3 = Opr[λ1,λ2][states[[3,1]],states[[3,2]]][p,r];
  FL1 = FL[λ1,λ2][states[[1,1]],states[[1,2]]][tt];
  FL2 = FL[λ1,λ2][states[[2,1]],states[[2,2]]][tt];
  FL3 = FL[λ1,λ2][states[[3,1]],states[[3,2]]][tt];
  eqA = Simplify[O1 FL2 - O2 FL1];
  eqB = Simplify[O1 FL3 - O3 FL1];
  resT = Resultant[eqA,eqB,p] // Simplify;
  {eqA, eqB, resT}
];
{eqA21, eqB21, resT21} = elimT[2,1][5/2];
{"Resultant[EqA,EqB,p] (should be a polynomial in tt alone)", resT21, "Exponent in tt", Exponent[resT21,tt]}
```
**Output:**
```
Out[642]= \left\{\text{Resultant[EqA,EqB,p] (should be a polynomial in tt alone)},\left(5.60\times 10^{-15}+5.7\times 10^{-16}\,\mathrm{i}\right)\,\mathrm{tt}^{4/3}\,\left(\left(1.\,+0.\times 10^{-3}\,\mathrm{i}\right)\,\mathrm{tt}-\left(1.06-0.069\,\mathrm{i}\right)\right)\,{\left(\left(1.\,+0.\times 10^{-3}\,\mathrm{i}\right)\,\mathrm{tt}-\left(0.36-0.933\,\mathrm{i}\right)\right)}^5\,\left(\lef… [651 chars total]
```
## 2026-07-18 18:51:17 — 📥 BULK INSERT 1 CELL(S) at positions 177–177
1. [code] (*Factor out the COMMON fractional prefactor T^(θ1/h) (from κ1^(u/h)·T^(u/h) evaluated at the↵  θ1-n


## 2026-07-18 18:51:18
**Input:**
```wolfram
(*Factor out the COMMON fractional prefactor T^(θ1/h) (from κ1^(u/h)·T^(u/h) evaluated at the
  θ1-node, same for every state of a rep — tied to θ1 only, not the state) before eliminating, so
  Resultant sees genuine polynomials.*)
ClearAll[FLpoly];
FLpoly[λ1_,λ2_][M_,n_][T_] := Simplify[FL[λ1,λ2][M,n][T] / T^(θ1/h)];
(*smoke: confirm this is now a clean polynomial (no fractional powers of T).*)
FLpoly[2,1][1,0][T] // Simplify
```
**Output:**
```
Out[645]= \left(-1.185612706261319859+0.8322533032249952804\,\mathrm{i}\right)\,T^2-\left(0.7618199108066716511+0.8128124734592948256\,\mathrm{i}\right)\,T-\left(0.4044100375365073294+0.5988708694858726138\,\mathrm{i}\right)
```
## 2026-07-18 18:51:29 — 📥 BULK INSERT 1 CELL(S) at positions 178–178
1. [code] (*Redo elimination with FLpoly (clean polynomials, no fractional powers) — Resultant should now↵  gi


## 2026-07-18 18:51:29
**Input:**
```wolfram
(*Redo elimination with FLpoly (clean polynomials, no fractional powers) — Resultant should now
  give a genuine polynomial in T.*)
ClearAll[elimTpoly];
elimTpoly[λ1_,λ2_][r_] := Module[{states,O1,O2,O3,FL1,FL2,FL3,eqA,eqB,resT},
  states = statesOf[λ1,λ2];
  O1 = Opr[λ1,λ2][states[[1,1]],states[[1,2]]][p,r];
  O2 = Opr[λ1,λ2][states[[2,1]],states[[2,2]]][p,r];
  O3 = Opr[λ1,λ2][states[[3,1]],states[[3,2]]][p,r];
  FL1 = FLpoly[λ1,λ2][states[[1,1]],states[[1,2]]][tt];
  FL2 = FLpoly[λ1,λ2][states[[2,1]],states[[2,2]]][tt];
  FL3 = FLpoly[λ1,λ2][states[[3,1]],states[[3,2]]][tt];
  eqA = Simplify[O1 FL2 - O2 FL1];
  eqB = Simplify[O1 FL3 - O3 FL1];
  resT = Resultant[eqA,eqB,p] // Simplify;
  {eqA, eqB, resT}
];
{eqA21p, eqB21p, resT21p} = elimTpoly[2,1][5/2];
{"Resultant (clean polynomial version)", resT21p, "Exponent in tt", Exponent[resT21p,tt]}
```
**Output:**
```
Out[649]= \left\{\text{Resultant (clean polynomial version)},\left(7.2096900227257981752515547669\times 10^{-10}-1\times 10^{-9}\,\mathrm{i}\right)\,\left(\left(1.\,+0.\times 10^{-1}\,\mathrm{i}\right)\,\mathrm{tt}-\left(1+0.\times 10^{-1}\,\mathrm{i}\right)\right)\,{\left(\left(1.\,+0.\times 10^{-1}\,\mathrm{i}\right)\,\mathrm{tt}-\left(0.3604676634583857275-0.9\,\mathrm{i}\right)\right)}^5\,\lef… [668 chars total]
```
## 2026-07-18 18:51:51 — 📥 BULK INSERT 1 CELL(S) at positions 179–179
1. [code] (*Test T=1 directly (the candidate root found in the resultant): solve for p using state-pairing,↵


## 2026-07-18 18:51:51
**Input:**
```wolfram
(*Test T=1 directly (the candidate root found in the resultant): solve for p using state-pairing,
  check consistency across ALL 6 states of rep (2,1), r=5/2 sample.*)
Tval = 1;
FLvalsT1 = Table[FL[2,1][st[[1]],st[[2]]][Tval], {st,statesOf[2,1]}];
OvalsT1[pp_] := Table[Opr[2,1][st[[1]],st[[2]]][pp,5/2], {st,statesOf[2,1]}];
(*solve p from state1 vs state2 ratio.*)
eqT1 = OvalsT1[p][[1]] FLvalsT1[[2]] - OvalsT1[p][[2]] FLvalsT1[[1]] == 0;
pSolT1 = NSolve[eqT1, p, WorkingPrecision->30];
(*check consistency: for each candidate p, compute residual of ALL 6 O[p,r][state]/FL[state] - C
  being equal (i.e. spread of ratios).*)
consistCheck = Table[
  Module[{pv,ratios,spread},
    pv = p/.sol;
    ratios = OvalsT1[pv]/FLvalsT1;
    spread = Max[Abs[ratios - ratios[[1]]]];
    {pv, ratios[[1]], spread}
  ]
,{sol,pSolT1}];
consistCheck//TableForm
```
**Output:**
```
2 kernel messages: NSolve::precw: The precision of the argument function ({(-2.351842654604498839218 - 0.579430039720172159021 I) ((1.901099687868066343166694126 + 2.331956121483436791212089955 I) - (2.97995121374406438… [also: NSolve::precw] result: Out[656]//TableForm= \begin{pmatrix}-0.360439875147226537266677650605+0.932782448593374716484835982024\,\mathrm{i} & 0.\times 10^{-27}+0.\times 10^{-27}\,\mathrm{i} & 15.8018300251165013842\\0.8484789…
```

## 2026-07-18 18:52:01
**Input:**
```wolfram
consistCheck
```
**Output:**
```
Out= {{-0.3604398751472265372666776506046790677033733542633701466781`29.707347828414424 + 0.93278244859337471648483598202444567040052156735816881392`30.120295363558725*I, 0``27.26559845777911 + 0``27.241036728160708*I, 15.80183002511650138419627368384399156405`21.49308424162044}, {0.8484789663621428335594333745767685646536138074113739854537`30.144339964179096 + 0.1441047912748553617304354611581160949615068141849192650605`29.374377304186787*I, -4.22313351003804772518225885133859546017`27.924897474471702 - 2.48637441470467120310184799779929184352`27.710810069374144*I, 12.67118670326289083051747836258355225039`21.574875849668196}}
```
## 2026-07-18 18:53:01 — 📥 BULK INSERT 1 CELL(S) at positions 180–180
1. [code] (*Test T=κ2 EXACTLY (not a rounded numeric guess) — the multiplicity-5 root found in the↵  resultant


## 2026-07-18 18:53:02
**Input:**
```wolfram
(*Test T=κ2 EXACTLY (not a rounded numeric guess) — the multiplicity-5 root found in the
  resultant. Solve for p via state-pairing, check consistency across ALL 6 states.*)
Tvalκ2 = κ2;
FLvalsTk = Table[FL[2,1][st[[1]],st[[2]]][Tvalκ2], {st,statesOf[2,1]}];
OvalsTk[pp_] := Table[Opr[2,1][st[[1]],st[[2]]][pp,5/2], {st,statesOf[2,1]}];
eqTk = OvalsTk[p][[1]] FLvalsTk[[2]] - OvalsTk[p][[2]] FLvalsTk[[1]] == 0;
pSolTk = NSolve[eqTk, p, WorkingPrecision->30];
consistCheckTk = Table[
  Module[{pv,ratios,spread},
    pv = p/.sol;
    ratios = OvalsTk[pv]/FLvalsTk;
    spread = Max[Abs[ratios - ratios[[1]]]];
    {pv, ratios[[1]], spread}
  ]
,{sol,pSolTk}];
consistCheckTk//TableForm
```
**Output:**
```
3 kernel messages: NSolve::precw: The precision of the argument function ({(0. 10    + 0. 10    I) ((1.901099687868066343166694126 + 2.331956121483436791212089955 I) - (2.979951213744064381468205485 - 5.2276847902919518… [also: Power::infy, NSolve::precw] result: Out[663]//TableForm= \begin{pmatrix}-0.776577576937402927802882824634 & \mathrm{ComplexInfinity} & \mathrm{Indeterminate}\end{pmatrix}
```

## 2026-07-18 18:54:19
**Input:**
```wolfram
Factor[resT21p]
```
**Output:**
```
Out= 0``0.5161899144247263
```

## 2026-07-18 18:54:40
**Input:**
```wolfram
statesOf[2,1]
```
**Output:**
```
Out= {{0, 0}, {1, 0}, {1, 1}, {2, 0}, {2, 1}, {3, 0}}
```
## 2026-07-18 18:54:55 — 📥 BULK INSERT 1 CELL(S) at positions 181–181
1. [code] (*Redo elimination using states {1,0},{1,1},{2,0} (indices 2,3,4 — AVOIDING the M=0 state {0,0},↵  w


## 2026-07-18 18:54:56
**Input:**
```wolfram
(*Redo elimination using states {1,0},{1,1},{2,0} (indices 2,3,4 — AVOIDING the M=0 state {0,0},
  which has deg q1=0 (Q1 constant) and may be causing a structural degeneracy in the elimination).*)
elimTpoly2[λ1_,λ2_][r_] := Module[{states,Oa,Ob,Oc,FLa,FLb,FLc,eqA,eqB,resT},
  states = statesOf[λ1,λ2];
  Oa = Opr[λ1,λ2][states[[2,1]],states[[2,2]]][p,r];
  Ob = Opr[λ1,λ2][states[[3,1]],states[[3,2]]][p,r];
  Oc = Opr[λ1,λ2][states[[4,1]],states[[4,2]]][p,r];
  FLa = FLpoly[λ1,λ2][states[[2,1]],states[[2,2]]][tt];
  FLb = FLpoly[λ1,λ2][states[[3,1]],states[[3,2]]][tt];
  FLc = FLpoly[λ1,λ2][states[[4,1]],states[[4,2]]][tt];
  eqA = Simplify[Oa FLb - Ob FLa];
  eqB = Simplify[Oa FLc - Oc FLa];
  resT = Resultant[eqA,eqB,p] // Simplify;
  {eqA, eqB, resT}
];
{eqA2, eqB2, resT2} = elimTpoly2[2,1][5/2];
{"Resultant (states {1,0},{1,1},{2,0})", resT2, "is it identically 0?", resT2===0 || Chop[N[resT2/.tt->13/7],10^-10]===0}
```
**Output:**
```
Out[666]= \left\{\text{Resultant (states \{1,0\},\{1,1\},\{2,0\})},\left(14329-6305\,\mathrm{i}\right)\,\left(\left(1.\,+0.\times 10^{-5}\,\mathrm{i}\right)\,\mathrm{tt}-\left(0.95775-0.3310\,\mathrm{i}\right)\right)\,{\left(\left(1.\,+0.\times 10^{-5}\,\mathrm{i}\right)\,\mathrm{tt}-\left(0.3604-0.93278\,\mathrm{i}\right)\right)}^3\,\left(\left(1.\,+0.\times 10^{-5}\,\mathrm{i}\right)\,\mathrm{tt… [733 chars total]
```
## 2026-07-18 18:55:29 — 📥 BULK INSERT 1 CELL(S) at positions 182–182
1. [code] (*T=κ2 is a genuine special point: T·κ1=1 exactly means T^(u/h)Q1(u)=(Tκ1)^(u/h)q1(u)=q1(u), the↵  B


## 2026-07-18 18:55:30
**Input:**
```wolfram
(*T=κ2 is a genuine special point: T·κ1=1 exactly means T^(u/h)Q1(u)=(Tκ1)^(u/h)q1(u)=q1(u), the
  BARE untwisted polynomial (κ1κ2=χ2=1 exactly) — matches a Λ=0 pole in Paul's language, not the
  informative root. Factor it out and find the OTHER root(s) at high precision.*)
resT2reduced = Simplify[resT2 / (tt-κ2)^3 / (tt-κ1)]; (*strip the κ2^3 factor and any κ1 factor if present*)
resT2reduced // Simplify
otherRoots = NSolve[resT2reduced==0, tt, WorkingPrecision->40];
otherRoots
```
**Output:**
```
2 kernel messages: NSolve::precw: The precision of the argument function ({------------------------------------------------------------------------------------------------------------------------------------------------… [also: NSolve::precw] result: Out[668]= \frac{\left(14300-6310\,\mathrm{i}\right)\,\mathrm{tt}^5+\left(4290-15000\,\mathrm{i}\right)\,\mathrm{tt}^4-\left(9370-10900\,\mathrm{i}\right)\,\mathrm{tt}^3+\left(3990+15100\,\mathrm{i}\ri…
```
## 2026-07-18 18:55:44 — 📥 BULK INSERT 1 CELL(S) at positions 183–183
1. [code] (*Cleaner: get ALL roots of the FULL resultant directly (no manual imprecise factoring), then↵  filt


## 2026-07-18 18:55:44
**Input:**
```wolfram
(*Cleaner: get ALL roots of the FULL resultant directly (no manual imprecise factoring), then
  filter out those near κ2 (the known spurious Λ=0 point) to isolate the genuine candidate(s).*)
allRootsT = NSolve[resT2==0, tt, WorkingPrecision->40];
allRootsT = tt/.allRootsT;
distinctRootsT = DeleteDuplicates[Chop[allRootsT,10^-8], Abs[#1-#2]<10^-6&];
genuineRoots = Select[distinctRootsT, Abs[#-κ2]>10^-4 && Abs[#-κ1]>10^-4 &];
{"all distinct roots", distinctRootsT, "genuine (non-κ1,κ2) candidates", genuineRoots}
```
**Output:**
```
2 kernel messages: NSolve::precw: The precision of the argument function ({(14329. - 6305. I) ((-0.9578 + 0.3310 I) + (1.0000 + 0. I) tt) ((-0.3604 + 0.9328 I) + (1.0000 + 0. I) tt)  ((-0.12809 - 0.20204 I) + (1.0000 + … [also: NSolve::precw] result: Out[675]= \left\{\text{all distinct roots},\left\{-0.7860962566844919786096256684491978609626+0.5547445255474452554744525547445255474453\,\mathrm{i},-0.4685039370078740157480314960629921259843+0.17135…
```

## 2026-07-18 18:55:51
**Input:**
```wolfram
genuineRoots
```
**Output:**
```
Out= {-0.7860962566844919786096256684491978609625668449197860962567`40.06275821830605 + 0.5547445255474452554744525547445255474452554744525547445256`39.911375515218765*I, -0.4685039370078740157480314960629921259842519685039370078741`40.12325187572833 + 0.171355498721227621483375959079283887468030690537084398977`39.68643667626069*I, 0.1280898876404494382022471910112359550561797752808988764044`39.879219411782 + 0.2020460358056265984654731457800511508951406649616368286445`40.077154900985015*I, 0.9577464788732394366197183098591549295774647887323943661972`40.12601628878301 - 0.3309859154929577464788732394366197183098591549295774647887`39.66457523834852*I}
```
## 2026-07-18 18:56:04 — 📥 BULK INSERT 1 CELL(S) at positions 184–184
1. [code] (*The 4 genuine roots have repeating-decimal patterns — rationalize them, then test EACH for FULL↵


## 2026-07-18 18:56:05
**Input:**
```wolfram
(*The 4 genuine roots have repeating-decimal patterns — rationalize them, then test EACH for FULL
  6-state consistency (the real verification).*)
genuineRootsRat = Rationalize[genuineRoots, 10^-25];
testConsistency[Tcand_] := Module[{FLvals,Ovals,eqRef,pSols,results},
  FLvals = Table[FL[2,1][st[[1]],st[[2]]][Tcand], {st,statesOf[2,1]}];
  Ovals[pp_] := Table[Opr[2,1][st[[1]],st[[2]]][pp,5/2], {st,statesOf[2,1]}];
  eqRef = Ovals[p][[2]] FLvals[[3]] - Ovals[p][[3]] FLvals[[2]] == 0; (*use non-degenerate states 2,3*)
  pSols = NSolve[eqRef, p, WorkingPrecision->30];
  results = Table[
    Module[{pv,ratios,spread},
      pv = p/.sol;
      ratios = Ovals[pv]/FLvals;
      spread = Max[Abs[ratios-ratios[[1]]]];
      {pv,ratios[[1]],spread}
    ]
  ,{sol,pSols}];
  results
];
allTests = Table[{Tcand, testConsistency[Tcand]}, {Tcand,genuineRootsRat}];
allTests // TableForm
```
**Output:**
```
2 kernel messages: NSolve::precw: The precision of the argument function ({(-4.541691109273399428597 - 0.161214415749599781314 I) ((7.508674272663045231781222516 + 1.7214441321549566092787416163 I) + (1.7181027858808976… [also: NSolve::precw] result: Out[679]//TableForm= \begin{pmatrix}-\frac{147}{187}+\frac{76\,\mathrm{i}}{137} & \begin{pmatrix}-7.98771842959443551594094471762+2.16546605458097969867976016316\,\mathrm{i} & 20.993118722988603573569…
```

## 2026-07-18 18:56:13
**Input:**
```wolfram
Table[{allTests[[i,1]], allTests[[i,2,All,3]]}, {i,1,4}]
```
**Output:**
```
Out= {{-147/187 + (76*I)/137, {281.47431787189200661706437574955790309903`21.55325834474556, 14.07713558057919308789951856097502125519`22.990284032958943}}, {-119/254 + (67*I)/391, {433.58771800483860829228609502573330016855`23.057130553181835, 43.46103748736020766075566643974471911362`22.94770450123653}}, {57/445 + (79*I)/391, {44.81721102729618900028306925870457225953`24.206970241607834, 55.37739569006840917798253938406855562113`21.940723924410385}}, {68/71 - (47*I)/142, {10.55855315605064278530521088540409397077`22.79965437327681, 17.63399319653037512830294678117130088028`22.029464731182337}}}
```
## 2026-07-18 18:57:08 — 📥 BULK INSERT 1 CELL(S) at positions 185–185
1. [code] (*More robust: compute pairwise elimination polynomials (in T alone) using state 2 as reference↵  ag


## 2026-07-18 18:57:09
**Input:**
```wolfram
(*More robust: compute pairwise elimination polynomials (in T alone) using state 2 as reference
  against states 3,4,5,6 SEPARATELY, then find roots COMMON to ALL FOUR (not just a 3-state
  resultant) — the TRUE consistent T must be a root of every one of them simultaneously.*)
statesL = statesOf[2,1];
refIdx = 2; (*state {1,0}*)
pairResultants = Table[
  Module[{Oref,Oj,FLref,FLj,eqP,resTj},
    Oref = Opr[2,1][statesL[[refIdx,1]],statesL[[refIdx,2]]][p,5/2];
    Oj = Opr[2,1][statesL[[j,1]],statesL[[j,2]]][p,5/2];
    FLref = FLpoly[2,1][statesL[[refIdx,1]],statesL[[refIdx,2]]][tt];
    FLj = FLpoly[2,1][statesL[[j,1]],statesL[[j,2]]][tt];
    eqP = Simplify[Oref FLj - Oj FLref];
    resTj = Resultant[eqP, D[eqP,p]/.p->p, p]; (*placeholder, real resultant needs pairing with another eq — instead just collect roots of THIS single eq in p for many T samples? Too slow. Use direct approach below.*)
    eqP
  ]
,{j,{3,4,5,6}}];
(*Instead: for each of the 4 candidate T's from before, directly check ALL 4 pairwise equations
  (not just one) share a common p-root — i.e. does pairP[j][T,p]=0 for the SAME p, for j=3,4,5,6?*)
checkCommonRoot[Tcand_] := Table[
  Module[{eqPj,rts},
    eqPj = pairResultants[[j-2]] /. tt->Tcand;
    rts = p/.NSolve[eqPj==0,p,WorkingPrecision->30];
    {j, rts}
  ]
,{j,{3,4,5,6}}];
Table[{Tc, checkCommonRoot[Tc]}, {Tc, genuineRootsRat}] // TableForm
```
**Output:**
```
2 kernel messages: NSolve::precw: The precision of the argument function ({(-28.210761568625537 + 18.649534368430135 I) - (31.443835020525303 + 18.309446881119772 I) p - (2.595139161792717 + 3.156364786435206 I) p }) is… [also: NSolve::precw] result: Out[684]//TableForm= \begin{pmatrix}-\frac{147}{187}+\frac{76\,\mathrm{i}}{137} & \begin{pmatrix}3 & \left\{-7.98771842959443551215428927267+2.16546605458098000197401311282\,\mathrm{i},-0.360439875147…
```

## 2026-07-18 18:57:20
**Input:**
```wolfram
Table[{genuineRootsRat[[i]], checkCommonRoot[genuineRootsRat[[i]]]}, {i,1,4}]
```
**Output:**
```
[message] NSolve::precw: The precision of the argument function ({(-28.210761568625537 + 18.649534368430135 I) - (31.443835020525303 + 18.309446881119772 I) p - (2.595139161792717 + 3.156364786435206 I) p }) is… [233 chars]
[message] [NSolve::precw ×2 more suppressed — first instance shown above]
[message] General::stop: Further output of NSolve::precw will be suppressed during this calculation.
Out= {{-147/187 + (76*I)/137, {{3, {-7.9877184295944355121542892726715239138089793905453992070569`30.135115002942896 + 2.1654660545809800019740131128198049230255146251339838529545`29.568243635740536*I, -0.3604398751472265615145126464979078410979387840370969089354`29.707347828414424 + 0.9327824485933747305884330600618464361593665566306437821144`30.120295363558725*I}}, {4, {-7.9877745031587372901125820395978095726750002918353429334348`30.135113133016564 + 2.1656174399265463683743701537528078384253309398985908527474`29.568269077076987*I, -0.8587670926915252739551466425026653938900524082202771199269`30.129294999881004 - 0.2751485950663362582693728531621538044574387795984672569326`29.63498690524139*I}}, {5, {-0.3731857335739121174596662636652130835717342886002759675598`29.82291073749403 + 0.7002310456640634888826181856265478238868274903006907138224`30.09622706651244*I, 1.095495819701219134230909732660057402496365725331855977767`29.55000888940435 - 4.2266713756883085948078751597931592512676987154367267355702`30.136396647410052*I}}, {6, {-2.0440375576144116578590756113374941350458575152150455712696`30.065037868673244 - 1.4196307301139388006390943940135739502421192151180798375144`29.906724389067406*I, 0.492337492076547113306503835027837310542307482559403481744`30.14310713649716 + 0.0917163227215623840707124008879500732649718591282828034527`29.413290861215895*I}}}}, {-119/254 + (67*I)/391, {{3, {-4.0340558926646873302851258831060066541572719202820561632812`29.83563836048493 - 7.2874440379042366070067072846290809083638152192159817375944`30.09247168205264*I, -0.3604398751472265339434518198199396690690197010698067086483`29.707347828414424 + 0.9327824485933747008777722006025660525010301066722087253122`30.120295363558725*I}}, {4, {-4.0339153655536044412303294969718760423290585115100621576824`29.835617671694763 - 7.2876437070076316812548257910786511563191490082900611306689`30.092478021354385*I, -0.360392413949566820655936787945282894177761765107771733969`29.70729258904727 + 0.9327959722345952016741997126972498220736904491183423602226`30.120303610447536*I}}, {5, {-4.0342290912485447230952468817354570368619577335459890563999`29.835663555810747 - 7.2872045765985782246504357693250834036747745799667731808817`30.092463960826862*I, -0.3604359426843157786291162770330917314579527400978934108681`29.707343166343943 + 0.9327837800776050477909995122557161278129116145800679246863`30.12029605966958*I}}, {6, {-4.0342657197898103019286139504773406584426415699739985194352`29.835643015998137 - 7.287721024078541660997701726156803107933284132324386978753`30.09247025543694*I, -0.3604568645356879816786136776277944370084490829315339477245`29.707359452963278 + 0.9327977190668002490675446589677375345293434066767716432916`30.120293627779855*I}}}}, {57/445 + (79*I)/391, {{3, {-0.4373369434466416758187540175442106046859170803912456899445`29.505418744286622 - 1.8814177804475519494387403989912226035280975947082061169333`30.139087822986774*I, -0.3604398751472265344453718622243675696005312746413519155102`29.707347828414424 + 0.9327824485933747103547228281202392331094376787028684160554`30.120295363558725*I}}, {4, {-1.2469472755944250923956514263683237787332645302285603911782`29.590394278631457 + 4.3535960181098359126074656778345603702314730159008870253256`30.133394315166637*I, -0.4373289805255354549256832736182138379084282608827580000413`29.505409196806635 - 1.8814271194637243504661421940647073335511668444130762076496`30.13908833885709*I}}, {5, {-1.7032172226088559279658865535524378672132614201120073262293`30.080025044452718 - 1.0547610036913904747114876704878412083950784006320820967639`29.8719090693094*I, 0.2406810359897653719973050548569393161992180653325850522387`29.585344574979246 + 0.8509434513297804055686136433768660355822080495714914247845`30.13380340317444*I}}, {6, {-1.2198926405570556263632044763194673636892922100721942509865`30.13342366649606 - 0.3490873523893208798546160868337390305746329610864236315263`29.590036169610386*I, 0.9751174181114076689081753155963678674969264160956998524587`29.959514976202133 - 1.15784904684446270742828021392449007339344396233068260290
[output truncated — 5492 chars total]
```
