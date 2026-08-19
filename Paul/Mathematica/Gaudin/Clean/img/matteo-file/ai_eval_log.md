# AI Action Log — matteo-file

## 2026-07-23 12:59:53 — ▶️ RUN CELL 1 (0.01 s)
**In [1]:** `ClearAll[fund, antifund, adj];↵(*Rep of gl(4)*)↵fund = {1, 0, 0, 0};↵antifund = {1, 1, 1, 0};↵adj = {2, 1, 1, 0};`
**Out:** (no output)

## 2026-07-23 12:59:53 — ▶️ RUN CELL 17
**In [17]:** `ClearAll[h];↵h[rep1_, rep2_][1, 1] := ↵  h[rep1, rep2][1, 1] = ↵   rep1[[1]] + rep1[[2]] + rep1[[3]] + rep1[[4]];↵h[rep1_, rep2_][2, 1] := ↵  h[rep1, rep2][2, 1] = ↵   rep1[[3]] (-1 + rep1[[4]]) + ↵  `
**Out:** (no output)


## 2026-07-23 12:59:56
**Input:**
```wolfram
oldVals = Table[{reps, k, h[reps[[1]], reps[[2]]][k, 1], h[reps[[1]], reps[[2]]][k, 2]}, {reps, {{fund, antifund}, {antifund, adj}, {adj, fund}, {fund, fund}, {antifund, antifund}, {adj, adj}}}, {k, 1, 4}]
```
**Output:**
```
Out= {{{{{1, 0, 0, 0}, {1, 1, 1, 0}}, 1, 1, 3}, {{{1, 0, 0, 0}, {1, 1, 1, 0}}, 2, -3, -3}, {{{1, 0, 0, 0}, {1, 1, 1, 0}}, 3, 6, 0}, {{{1, 0, 0, 0}, {1, 1, 1, 0}}, 4, -6, 0}}, {{{{1, 1, 1, 0}, {2, 1, 1, 0}}, 1, 3, 4}, {{{1, 1, 1, 0}, {2, 1, 1, 0}}, 2, -3, -4}, {{{1, 1, 1, 0}, {2, 1, 1, 0}}, 3, 0, 0}, {{{1, 1, 1, 0}, {2, 1, 1, 0}}, 4, 0, 0}}, {{{{2, 1, 1, 0}, {1, 0, 0, 0}}, 1, 4, 1}, {{{2, 1, 1, 0}, {1, 0, 0, 0}}, 2, -4, -3}, {{{2, 1, 1, 0}, {1, 0, 0, 0}}, 3, 0, 6}, {{{2, 1, 1, 0}, {1, 0, 0, 0}}, 4, 0, -6}}, {{{{1, 0, 0, 0}, {1, 0, 0, 0}}, 1, 1, 1}, {{{1, 0, 0, 0}, {1, 0, 0, 0}}, 2, -3, -3}, {{{1, 0, 0, 0}, {1, 0, 0, 0}}, 3, 6, 6}, {{{1, 0, 0, 0}, {1, 0, 0, 0}}, 4, -6, -6}}, {{{{1, 1, 1, 0}, {1, 1, 1, 0}}, 1, 3, 3}, {{{1, 1, 1, 0}, {1, 1, 1, 0}}, 2, -3, -3}, {{{1, 1, 1, 0}, {1, 1, 1, 0}}, 3, 0, 0}, {{{1, 1, 1, 0}, {1, 1, 1, 0}}, 4, 0, 0}}, {{{{2, 1, 1, 0}, {2, 1, 1, 0}}, 1, 4, 4}, {{{2, 1, 1, 0}, {2, 1, 1, 0}}, 2, -4, -4}, {{{2, 1, 1, 0}, {2, 1, 1, 0}}, 3, 0, 0}, {{{2, 1, 1, 0}, {2, 1, 1, 0}}, 4, 0, 0}}}
```

## 2026-07-23 13:00:08
**Input:**
```wolfram
ClearAll[htest];
htest[rep_, a_] := Sum[Product[rep[[subset[[m]]]] - (m - 1), {m, 1, a}], {subset, Subsets[Range[Length[rep]], {a}]}];
newVals = Table[{reps, k, htest[reps[[1]], k], htest[reps[[2]], k]}, {reps, {{fund, antifund}, {antifund, adj}, {adj, fund}, {fund, fund}, {antifund, antifund}, {adj, adj}}}, {k, 1, 4}];
newVals === oldVals
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= True
```
## 2026-07-23 13:00:36 — ✏️ EDIT CELL 17 [batch]
ClearAll[h];
(* h[rep1,rep2][a,a,0] and h[rep1,rep2][a,0,a]: order-a elementary-symmetric-type invariants of a single
   site's weight/Young-diagram vector, closed form for general a (subset size). Fo…


## 2026-07-23 13:00:36
**Input:**
```wolfram
ClearAll[h];
(* h[rep1,rep2][a,a,0] and h[rep1,rep2][a,0,a]: order-a elementary-symmetric-type invariants of a single
   site's weight/Young-diagram vector, closed form for general a (subset size). For a subset
   {i1<...<ia} of Range[Length[rep]], the m-th chosen entry is shifted by its RANK within the
   subset (m-1), not its absolute index -- this single formula reproduces all of k=1..4 for both
   sites (verified against the old hand-expanded k=1,2,3,4 formulas below before replacing them). *)
h[rep1_, rep2_][a_, a_, 0] :=
  h[rep1, rep2][a, a, 0] =
   Sum[Product[rep1[[subset[[m]]]] - (m - 1), {m, 1, a}], {subset, Subsets[Range[Length[rep1]], {a}]}];
h[rep1_, rep2_][a_, 0, a_] :=
  h[rep1, rep2][a, 0, a] =
   Sum[Product[rep2[[subset[[m]]]] - (m - 1), {m, 1, a}], {subset, Subsets[Range[Length[rep2]], {a}]}];
```
**Output:**
```
Out= (no output)
```
## 2026-07-23 13:00:36 — ✏️ EDIT CELL 19 [batch]
ClearAll[ϵ1, ϵ2, ϵ3];
\[Epsilon]1[rep1_, rep2_, \[Phi]1_, \[Phi]2_, \[Phi]3_, \[Phi]4_][
   z_] := \[Epsilon]1[rep1, rep2, \[Phi]1, \[Phi]2, \[Phi]3, \[Phi]4][
    z] = h[rep1, rep2][1, 1, 0]/(z - w1)…


## 2026-07-23 13:00:37
**Input:**
```wolfram
ClearAll[ϵ1, ϵ2, ϵ3];
\[Epsilon]1[rep1_, rep2_, \[Phi]1_, \[Phi]2_, \[Phi]3_, \[Phi]4_][
   z_] := \[Epsilon]1[rep1, rep2, \[Phi]1, \[Phi]2, \[Phi]3, \[Phi]4][
    z] = h[rep1, rep2][1, 1, 0]/(z - w1) + h[rep1, rep2][1, 0, 1]/(z - w2) +
     h[\[Phi]1, \[Phi]2, \[Phi]3, \[Phi]4][1, 3];
\[Epsilon]2[rep1_, rep2_, 
    repAll_, \[Phi]1_, \[Phi]2_, \[Phi]3_, \[Phi]4_][
   z_] := \[Epsilon]2[rep1, rep2, 
     repAll, \[Phi]1, \[Phi]2, \[Phi]3, \[Phi]4][z] = 
   h[rep1, rep2][2, 2, 0]/(z - w1)^2 + h[rep1, rep2][2, 0, 2]/(z - w2)^2 + 
    h[rep1, rep2, repAll, \[Phi]1, \[Phi]2, \[Phi]3, \[Phi]4][2, 
      3]/(z - w1) + h[2, 4]/(z - w2) + 
    h[\[Phi]1, \[Phi]2, \[Phi]3, \[Phi]4][2, 5];
\[Epsilon]3[rep1_, rep2_, 
    repAll_, \[Phi]1_, \[Phi]2_, \[Phi]3_, \[Phi]4_][
   z_] := \[Epsilon]3[rep1, rep2, 
     repAll, \[Phi]1, \[Phi]2, \[Phi]3, \[Phi]4][z] = 
   h[rep1, rep2][3, 3, 0]/(z - w1)^3 + h[rep1, rep2][3, 0, 3]/(z - w2)^3 + 
    h[3, 3]/(z - w1)^2 + h[3, 4]/(z - w2)^2 + 
    h[rep1, rep2, repAll, \[Phi]1, \[Phi]2, \[Phi]3, \[Phi]4][3, 
      5]/(z - w1) + h[3, 6]/(z - w2) + 
    h[\[Phi]1, \[Phi]2, \[Phi]3, \[Phi]4][3, 7];
\[Epsilon]4[rep1_, rep2_, 
    repAll_, \[Phi]1_, \[Phi]2_, \[Phi]3_, \[Phi]4_][
   z_] := \[Epsilon]4[rep1, rep2, 
     repAll, \[Phi]1, \[Phi]2, \[Phi]3, \[Phi]4][z] = 
   h[rep1, rep2][4, 4, 0]/(z - w1)^4 + h[rep1, rep2][4, 0, 4]/(z - w2)^4 + 
    h[4, 3]/(z - w1)^3 + h[4, 4]/(z - w2)^3 + h[4, 5]/(z - w1)^2 + 
    h[4, 6]/(z - w2)^2 + 
    h[rep1, rep2, repAll, \[Phi]1, \[Phi]2, \[Phi]3, \[Phi]4][4, 
      7]/(z - w1) + h[4, 8]/(z - w2) + 
    h[\[Phi]1, \[Phi]2, \[Phi]3, \[Phi]4][4, 9];
```
**Output:**
```
Out= (no output)
```
## 2026-07-23 13:00:37 — ✏️ EDIT CELL 20 [batch]
(* H1 *)
h[φ1_, φ2_, φ3_, φ4_][1, 3] :=
    h[φ1, φ2, φ3, φ4][1, 3] = φ1 + φ2 + φ3 + φ4;

(* H2 *)
h[φ1_, φ2_, φ3_, φ4_][2, 5] :=
    h[φ1, φ2, φ3, φ4][2, 5] =
     φ1 φ2 + φ1 φ3 + φ2 φ3 + φ1 φ4 + φ2 …


## 2026-07-23 13:00:37
**Input:**
```wolfram
(* H1 *)
h[φ1_, φ2_, φ3_, φ4_][1, 3] :=
    h[φ1, φ2, φ3, φ4][1, 3] = φ1 + φ2 + φ3 + φ4;

(* H2 *)
h[φ1_, φ2_, φ3_, φ4_][2, 5] :=
    h[φ1, φ2, φ3, φ4][2, 5] =
     φ1 φ2 + φ1 φ3 + φ2 φ3 + φ1 φ4 + φ2 φ4 + φ3 φ4;

(* H3 *)
h[φ1_, φ2_, φ3_, φ4_][3, 7] :=
    h[φ1, φ2, φ3, φ4][3, 7] =
     φ1 φ2 φ3 + φ1 φ2 φ4 + φ1 φ3 φ4 + φ2 φ3 φ4;

(* H4 *)
h[φ1_, φ2_, φ3_, φ4_][4, 9] :=
    h[φ1, φ2, φ3, φ4][4, 9] = φ1*φ2*φ3*φ4;

(* H2 *)
h[rep1_, rep2_, repAll_, φ1_, φ2_, φ3_, φ4_][2, 3] :=
    h[rep1, rep2, repAll, φ1, φ2, φ3, φ4][2, 3] =
     -h[2, 4] + (φ2 + φ3 + φ4) (repAll[[1]]) +
         (φ1 + φ3 + φ4) (repAll[[2]]) +
         (φ1 + φ2 + φ4) (repAll[[3]]) +
         (φ1 + φ2 + φ3) (repAll[[4]]);

(* H3 *)
h[rep1_, rep2_, repAll_, φ1_, φ2_, φ3_, φ4_][3, 5] :=
    h[rep1, rep2, repAll, φ1, φ2, φ3, φ4][3, 5] =
     -h[3, 6] + φ3 φ4 (repAll[[1]] + repAll[[2]]) +
         φ2 φ4 (repAll[[1]] + repAll[[3]]) +
         φ1 φ4 (repAll[[2]] + repAll[[3]]) +
         φ2 φ3 (repAll[[1]] + repAll[[4]]) +
         φ1 φ3 (repAll[[2]] + repAll[[4]]) +
         φ1 φ2 (repAll[[3]] + repAll[[4]]);

(* H4 *)
h[rep1_, rep2_, repAll_, φ1_, φ2_, φ3_, φ4_][4, 7] :=
    h[rep1, rep2, repAll, φ1, φ2, φ3, φ4][4, 7] =
     -h[4, 8] + φ2 φ3 φ4 (repAll[[1]]) +
         φ1 φ3 φ4 (repAll[[2]]) +
         φ1 φ2 φ4 (repAll[[3]]) +
         φ1 φ2 φ3 (repAll[[4]]);

ClearAll[ϵ1, ϵ2, ϵ3];
ϵ1[rep1_, rep2_, φ1_, φ2_, φ3_, φ4_][z_] :=
    ϵ1[rep1, rep2, φ1, φ2, φ3, φ4][z] =
     h[rep1, rep2][1, 1, 0]/(z - w1) + h[rep1, rep2][1, 0, 1]/(z - w2) +
        h[φ1, φ2, φ3, φ4][1, 3];

ϵ2[rep1_, rep2_, repAll_, φ1_, φ2_, φ3_, φ4_][z_] :=
    ϵ2[rep1, rep2, repAll, φ1, φ2, φ3, φ4][z] =
     h[rep1, rep2][2, 2, 0]/(z - w1)^2 + h[rep1, rep2][2, 0, 2]/(z - w2)^2 +
        h[rep1, rep2, repAll, φ1, φ2, φ3, φ4][2, 3]/(z - w1) +
        h[2, 4]/(z - w2) + h[φ1, φ2, φ3, φ4][2, 5];

ϵ3[rep1_, rep2_, repAll_, φ1_, φ2_, φ3_, φ4_][z_] :=
    ϵ3[rep1, rep2, repAll, φ1, φ2, φ3, φ4][z] =
     h[rep1, rep2][3, 3, 0]/(z - w1)^3 + h[rep1, rep2][3, 0, 3]/(z - w2)^3 +
        h[3, 3]/(z - w1)^2 + h[3, 4]/(z - w2)^2 +
        h[rep1, rep2, repAll, φ1, φ2, φ3, φ4][3, 5]/(z - w1) +
        h[3, 6]/(z - w2) + h[φ1, φ2, φ3, φ4][3, 7];

ϵ4[rep1_, rep2_, repAll_, φ1_, φ2_, φ3_, φ4_][z_] :=
    ϵ4[rep1, rep2, repAll, φ1, φ2, φ3, φ4][z] =
     h[rep1, rep2][4, 4, 0]/(z - w1)^4 + h[rep1, rep2][4, 0, 4]/(z - w2)^4 +
        h[4, 3]/(z - w1)^3 + h[4, 4]/(z - w2)^3 + h[4, 5]/(z - w1)^2 +
        h[4, 6]/(z - w2)^2 +
        h[rep1, rep2, repAll, φ1, φ2, φ3, φ4][4, 7]/(z - w1) +
        h[4, 8]/(z - w2) + h[φ1, φ2, φ3, φ4][4, 9];
```
**Output:**
```
Out= (no output)
```

## 2026-07-23 13:01:02
**Input:**
```wolfram
newValsFromNotebook = Table[{reps, k, h[reps[[1]], reps[[2]]][k, k, 0], h[reps[[1]], reps[[2]]][k, 0, k]}, {reps, {{fund, antifund}, {antifund, adj}, {adj, fund}, {fund, fund}, {antifund, antifund}, {adj, adj}}}, {k, 1, 4}];
(* compare against the coordinates recorded before the refactor (old [k,1]/[k,2] site values) *)
oldValsCheck = {{{{{1, 0, 0, 0}, {1, 1, 1, 0}}, 1, 1, 3}, {{{1, 0, 0, 0}, {1, 1, 1, 0}}, 2, -3, -3}, {{{1, 0, 0, 0}, {1, 1, 1, 0}}, 3, 6, 0}, {{{1, 0, 0, 0}, {1, 1, 1, 0}}, 4, -6, 0}}, {{{{1, 1, 1, 0}, {2, 1, 1, 0}}, 1, 3, 4}, {{{1, 1, 1, 0}, {2, 1, 1, 0}}, 2, -3, -4}, {{{1, 1, 1, 0}, {2, 1, 1, 0}}, 3, 0, 0}, {{{1, 1, 1, 0}, {2, 1, 1, 0}}, 4, 0, 0}}, {{{{2, 1, 1, 0}, {1, 0, 0, 0}}, 1, 4, 1}, {{{2, 1, 1, 0}, {1, 0, 0, 0}}, 2, -4, -3}, {{{2, 1, 1, 0}, {1, 0, 0, 0}}, 3, 0, 6}, {{{2, 1, 1, 0}, {1, 0, 0, 0}}, 4, 0, -6}}, {{{{1, 0, 0, 0}, {1, 0, 0, 0}}, 1, 1, 1}, {{{1, 0, 0, 0}, {1, 0, 0, 0}}, 2, -3, -3}, {{{1, 0, 0, 0}, {1, 0, 0, 0}}, 3, 6, 6}, {{{1, 0, 0, 0}, {1, 0, 0, 0}}, 4, -6, -6}}, {{{{1, 1, 1, 0}, {1, 1, 1, 0}}, 1, 3, 3}, {{{1, 1, 1, 0}, {1, 1, 1, 0}}, 2, -3, -3}, {{{1, 1, 1, 0}, {1, 1, 1, 0}}, 3, 0, 0}, {{{1, 1, 1, 0}, {1, 1, 1, 0}}, 4, 0, 0}}, {{{{2, 1, 1, 0}, {2, 1, 1, 0}}, 1, 4, 4}, {{{2, 1, 1, 0}, {2, 1, 1, 0}}, 2, -4, -4}, {{{2, 1, 1, 0}, {2, 1, 1, 0}}, 3, 0, 0}, {{{2, 1, 1, 0}, {2, 1, 1, 0}}, 4, 0, 0}}};
newValsFromNotebook === oldValsCheck
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= True
```
## 2026-07-23 13:40:39 — 📥 BULK INSERT 1 CELL(S) at positions 1–1
1. [markdown] ## Setup: gl(4) representations

## 2026-07-23 13:40:40 — 📥 BULK INSERT 1 CELL(S) at positions 4–4
1. [markdown] ## Gelfand–Tsetlin basis & generators

## 2026-07-23 13:40:41 — 📥 BULK INSERT 1 CELL(S) at positions 6–6
1. [markdown] ## Matrix generators

## 2026-07-23 13:40:42 — 📥 BULK INSERT 1 CELL(S) at positions 11–11
1. [markdown] ## Weight-finding utilities

## 2026-07-23 13:40:43 — 📥 BULK INSERT 1 CELL(S) at positions 16–16
1. [markdown] ## Twist parameters & fusion-hierarchy Hamiltonians

## 2026-07-23 13:40:44 — 📥 BULK INSERT 1 CELL(S) at positions 18–18
1. [markdown] ## Site Hamiltonians h[rep1,rep2][...] and ε_a(z)


---

## 2026-07-23 13:53:37 — 🔄 KERNEL RESTART

## 2026-07-23 17:26:54 — ▶️ RUN CELL 30
**In [30]:** `ClearAll[F];↵F[rep1_List, rep2_List][z_]:=↵    1/((z - w1)^(rep1[[1]]) * (z - w2)^(rep2[[1]]));`
**Out:** (no output)

## 2026-07-23 17:26:56 — ▶️ RUN CELL 33 (0.01 s)
**In [33]:** `(*asymptotics of Q's for a given triple of reps*)↵M[j_][rep1_List,rep2_List,repall_List]:=rep1[[1]]+rep2[[1]]-repall[[j]];`
**Out:** (no output)

## 2026-07-23 17:27:18 — 📥 BULK INSERT 2 CELL(S) at positions 36–37
1. [markdown] # Baxter solver↵↵Self-contained solver section (needs the Hamiltonian rule cells, `F` and `M` above 
2. [code] (*Shifted Hamiltonian coefficients: acting with \[PartialD]^4+\[CapitalSigma] H[a]\[PartialD]^(4-a)


## 2026-07-23 17:27:18
**Input:**
```wolfram
(*Shifted Hamiltonian coefficients: acting with \[PartialD]^4+\[CapitalSigma] H[a]\[PartialD]^(4-a) on F(z)q(z), the coefficient of q^(4-k)/F is HSc[k] = \[CapitalSigma]_a Binomial[4-a,k-a] H[a] F^(k-a)/F. Reproduces the \[Epsilon]kS of Matteo's file, including the k=4 term never extracted above.*)
ClearAll[HSc, hUnknowns];
HSc[k_][rep1_List, rep2_List, repall_List][z_] := HSc[k][rep1, rep2, repall][z] =
  Module[{zz},
   Sum[Binomial[4 - aa, k - aa]*If[aa == 0, 1, H[aa][zz]]*
      D[F[rep1, rep2][zz], {zz, k - aa}]/F[rep1, rep2][zz], {aa, 0, k}] //. 
     expandH /. applyleadingH /. applycentralH /. applyrep[rep1, rep2, repall] /. 
    zz -> z // Together];
(*the 9 moments not fixed by centrality or the large-z sum rule*)
hUnknowns = {h[2,1,0], h[3,1,0], h[3,2,0], h[3,0,2], h[4,1,0], h[4,2,0], h[4,3,0], h[4,0,2], h[4,0,3]};
```
**Output:**
```
(no output)
```

## 2026-07-23 17:27:22
**Input:**
```wolfram
HSc[1][fund, antifund, adj][z]
```
**Output:**
```
Out= (3 - 4*z - z*ϕ[1] + z^2*ϕ[1] - z*ϕ[2] + z^2*ϕ[2] - z*ϕ[3] + z^2*ϕ[3] - z*ϕ[4] + z^2*ϕ[4])/((-1 + z)*z)
```

## 2026-07-23 17:27:30
**Input:**
```wolfram
{Complement[Cases[HSc[4][fund, antifund, adj][z], h[__], Infinity] // DeleteDuplicates, hUnknowns], Exponent[Denominator[HSc[4][fund, antifund, adj][z]], z]}
```
**Output:**
```
Out= {{}, 6}
```
## 2026-07-23 17:27:48 — 📥 BULK INSERT 1 CELL(S) at positions 38–38
1. [code] (*monic polynomial part q_a of the a-th Q-function; the exp/F factors are handled by the operator*)↵


## 2026-07-23 17:27:49
**Input:**
```wolfram
(*monic polynomial part q_a of the a-th Q-function; the exp/F factors are handled by the operator*)
ClearAll[qAns, dShift, BaxEqn];
qAns[a_][rep1_List, rep2_List, repall_List][z_] :=
  Module[{M = M[a][rep1, rep2, repall]},
   z^M + Sum[c[a, n] z^n, {n, 0, M - 1}]];

(*(\[PartialD]-\[Phi])^m applied to a polynomial in z: shift accounts for the stripped e^{-\[Phi] z}*)
dShift[expr_, phi_, m_, z_] := Nest[(D[#, z] - phi #) &, expr, m];

(*a-th Baxter equation, polynomial in z (F and e^{-\[Phi]_a z} already divided out)*)
BaxEqn[a_][rep1_List, rep2_List, repall_List, phis_List][z_] :=
  Module[{q, hsc},
   q = qAns[a][rep1, rep2, repall][z];
   Sum[(HSc[k][rep1, rep2, repall][z] /. HSc0 -> 1 /. Thread[(ϕ /@ Range[4]) -> phis]) *
       dShift[q, phis[[a]], 4 - k, z], {k, 0, 4}] /.
     HSc[0][__][_] -> 1 // Together // Numerator // Expand];
(*note HSc[0]=1: the leading \[PartialD]^4 term*)
```
**Output:**
```
(no output)
```

## 2026-07-23 17:27:58
**Input:**
```wolfram
{HSc[0][fund,antifund,adj][z], BaxEqn[1][fund, antifund, adj, {1, 1/3, 1/5, 1/7}][z]}
```
**Output:**
```
Out= {1, -834 + 1946*z - 279*z^2 - 1182*z^3 - 461*z^4 - 210*h[2, 1, 0] + 420*z*h[2, 1, 0] - 105*z^2*h[2, 1, 0] - 210*z^3*h[2, 1, 0] - 105*z^4*h[2, 1, 0] + 105*z^2*h[3, 0, 2] - 105*z^3*h[3, 0, 2] - 105*z^4*h[3, 0, 2] + 105*z*h[3, 1, 0] - 210*z^2*h[3, 1, 0] + 105*z^4*h[3, 1, 0] + 105*h[3, 2, 0] - 315*z*h[3, 2, 0] + 210*z^2*h[3, 2, 0] + 105*z^3*h[3, 2, 0] - 105*z^4*h[3, 2, 0] - 105*z^3*h[4, 0, 2] + 105*z^4*h[4, 0, 2] + 105*z^3*h[4, 0, 3] - 105*z^2*h[4, 1, 0] + 210*z^3*h[4, 1, 0] - 105*z^4*h[4, 1, 0] - 105*z*h[4, 2, 0] + 315*z^2*h[4, 2, 0] - 315*z^3*h[4, 2, 0] + 105*z^4*h[4, 2, 0] - 105*h[4, 3, 0] + 315*z*h[4, 3, 0] - 315*z^2*h[4, 3, 0] + 105*z^3*h[4, 3, 0]}
```
## 2026-07-23 17:28:24 — 📥 BULK INSERT 2 CELL(S) at positions 39–40
1. [code] (*Joint solver. All four Baxter equations are polynomial in z; Solve handles the↵  full (mildly nonl
2. [code] (*Independent multiplicity: number of times weight repall appears in rep1⊗rep2,↵  using the notebook


## 2026-07-23 17:28:24
**Input:**
```wolfram
(*Joint solver. All four Baxter equations are polynomial in z; Solve handles the
  full (mildly nonlinear, via c*h cross terms) system. Returns one rule-list per
  Gaudin eigenstate carrying global weight repall.*)
ClearAll[SolveBax, cVars];
cVars[rep1_, rep2_, repall_] :=
  Flatten[Table[Table[c[a, n], {n, 0, M[a][rep1, rep2, repall] - 1}], {a, 1, 4}]];
SolveBax[rep1_List, rep2_List, repall_List, phis_List] :=
  Module[{Ms, eqs, vars},
   Ms = Table[M[a][rep1, rep2, repall], {a, 1, 4}];
   If[Min[Ms] < 0, Message[SolveBax::badrep, Ms]; Return[$Failed]];
   eqs = Flatten[Table[
       Thread[CoefficientList[BaxEqn[a][rep1, rep2, repall, phis][z], z] == 0], {a, 1, 4}]];
   vars = Join[cVars[rep1, rep2, repall], hUnknowns];
   Solve[eqs, vars]];
SolveBax::badrep = "repall not in tensor product: M = `1`";
```
**Output:**
```
(no output)
```

## 2026-07-23 17:28:25
**Input:**
```wolfram
(*Independent multiplicity: number of times weight repall appears in rep1⊗rep2,
  using the notebook's own weight enumerator. = expected number of Baxter branches.*)
ClearAll[weightMult];
weightMult[rep1_, rep2_, wt_] := Count[FindWeightAll[rep1, rep2], wt];
```
**Output:**
```
(no output)
```

## 2026-07-23 17:28:35
**Input:**
```wolfram
phiTest = {1, 1/3, 1/5, 1/7};
solAdj = SolveBax[fund, antifund, adj, phiTest];
{Length[solAdj], weightMult[fund, antifund, adj], M[#][fund,antifund,adj] & /@ Range[4]}
```
**Output:**
```
Out= {1, 1, {0, 1, 1, 2}}
```

## 2026-07-23 17:28:41
**Input:**
```wolfram
Table[BaxEqn[a][fund, antifund, adj, phiTest][z] /. solAdj[[1]] // Expand, {a, 1, 4}]
```
**Output:**
```
Out= {0, 0, 0, 0}
```

## 2026-07-23 17:28:59
**Input:**
```wolfram
wts = FindWeightDistinct[fund, antifund];
sweep = Table[Module[{s = SolveBax[fund, antifund, wt, phiTest]},
    {wt, weightMult[fund, antifund, wt], If[MatchQ[s, _List], Length[s], s],
     If[MatchQ[s, _List] && s =!= {},
       Max[Table[Max@Abs@CoefficientList[BaxEqn[a][fund,antifund,wt,phiTest][z] /. s[[b]], z],
         {a,1,4}, {b,1,Length[s]}]], 0]}], {wt, wts}];
Column[{
  "columns: {weight, weightMult, nBranches, maxResidual}",
  Grid[sweep, Frame -> All],
  "sum nBranches (should be 16): " <> ToString[Total[sweep[[All, 3]] /. Except[_Integer] -> 0]],
  "all match & residual 0: " <> ToString[
    And @@ Table[sweep[[i,2]] === sweep[[i,3]] && sweep[[i,4]] === 0, {i, Length[sweep]}]]}]
```
**Output:**
```
[message] N::meprec: Internal precision limit $MaxExtraPrecision = 50. reached while evaluating 630 - ------------------------------------------------------------------------------------------------------------… [1922 chars]
[message] [N::meprec ×2 more suppressed — first instance shown above]
[message] General::stop: Further output of N::meprec will be suppressed during this calculation.
[message] Max::meprec: Internal precision limit $MaxExtraPrecision = 50. reached while evaluating Abs[630 - ------------------------------------------------------------------------------------------------------… [1935 chars]
[message] [Max::meprec ×2 more suppressed — first instance shown above]
[message] [General::stop ×1 more suppressed — first instance shown above]
Out= Column[{"columns: {weight, weightMult, nBranches, maxResidual}", Grid[{{{2, 1, 1, 0}, 1, 1, -Infinity}, {{2, 1, 0, 1}, 1, 1, -Infinity}, {{2, 0, 1, 1}, 1, 1, -Infinity}, {{1, 1, 1, 1}, 4, 4, Max[Abs[5670 - 3726*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0] - 4*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]*(-4725 - 792*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0] + 344*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^2 - 16*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^3) + 2*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]*(4725 + 792*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0] - 344*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^2 + 16*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^3)], Abs[15750 - (1375*(63 + 444*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0] + 112*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^2 - 8*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^3))/27 - (50*(-4725 - 792*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0] + 344*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^2 - 16*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^3)*(63 + 444*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0] + 112*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^2 - 8*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^3))/729 + (25*(63 + 444*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0] + 112*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^2 - 8*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^3)*(4725 + 792*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0] - 344*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^2 + 16*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^3))/729], Abs[630 - (139*(-315 - 672*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0] + 184*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^2 - 8*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^3))/81 - (2*(-4725 - 792*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0] + 344*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^2 - 16*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^3)*(-315 - 672*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0] + 184*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^2 - 8*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^3))/2187 + ((-315 - 672*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0] + 184*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^2 - 8*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^3)*(4725 + 792*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0] - 344*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^2 + 16*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^3))/2187], Abs[30870 - (4802*(-180 - 573*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0] - 260*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^2 + 16*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^3))/81 - (196*(-4725 - 792*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0] + 344*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^2 - 16*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^3)*(-180 - 573*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0] - 260*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^2 + 16*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^3))/2187 + (98*(4725 + 792*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0] - 344*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^2 + 16*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^3)*(-180 - 573*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0] - 260*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^4 & , 1, 0]^2 + 16*Root[945 + 4104*#1 + 792*#1^2 - 344*#1^3 + 16*#1^
[output truncated — 249106 chars total]
```

## 2026-07-23 17:29:14
**Input:**
```wolfram
sweep2 = Table[Module[{s = SolveBax[fund, antifund, wt, phiTest], res},
    res = If[MatchQ[s, _List] && s =!= {},
      Max[Table[Max[Abs[N[CoefficientList[BaxEqn[a][fund,antifund,wt,phiTest][z] /. s[[b]], z], 30]]],
         {a,1,4}, {b,1,Length[s]}]], 0];
    {wt, weightMult[fund, antifund, wt], Length[s], Chop[res, 10^-20]}], {wt, wts}];
{Grid[Prepend[sweep2, {"weight","mult","nBranch","maxResidual"}], Frame -> All],
 "sum nBranches = " <> ToString[Total[sweep2[[All,3]]]],
 "mult==nBranch all: " <> ToString[And @@ (#[[2]]===#[[3]] & /@ sweep2)],
 "all residuals 0: " <> ToString[And @@ (#[[4]]===0 & /@ sweep2)]}
```
**Output:**
```
[message] N::meprec: Internal precision limit $MaxExtraPrecision = 50. reached while evaluating {630 - -----------------------------------------------------------------------------------------------------------… [27678 chars]
[message] [N::meprec ×2 more suppressed — first instance shown above]
[message] General::stop: Further output of N::meprec will be suppressed during this calculation.
Out= {Grid[{{"weight", "mult", "nBranch", "maxResidual"}, {{2, 1, 1, 0}, 1, 1, -Infinity}, {{2, 1, 0, 1}, 1, 1, -Infinity}, {{2, 0, 1, 1}, 1, 1, -Infinity}, {{1, 1, 1, 1}, 4, 4, 0}, {{1, 2, 1, 0}, 1, 1, -Infinity}, {{1, 2, 0, 1}, 1, 1, -Infinity}, {{0, 2, 1, 1}, 1, 1, -Infinity}, {{1, 1, 2, 0}, 1, 1, -Infinity}, {{1, 0, 2, 1}, 1, 1, -Infinity}, {{0, 1, 2, 1}, 1, 1, -Infinity}, {{1, 1, 0, 2}, 1, 1, -Infinity}, {{1, 0, 1, 2}, 1, 1, -Infinity}, {{0, 1, 1, 2}, 1, 1, -Infinity}}, Frame -> All], "sum nBranches = 16", "mult==nBranch all: True", "all residuals 0: False"}
```

## 2026-07-23 17:29:29
**Input:**
```wolfram
residClean = sweep2[[All,4]] /. -Infinity -> 0;
{"weights tested: " <> ToString[Length[wts]],
 "sum branches = dim(4x4): " <> ToString[Total[sweep2[[All,3]]] == 16],
 "mult==nBranch everywhere: " <> ToString[And @@ (#[[2]]===#[[3]] & /@ sweep2)],
 "all residuals vanish: " <> ToString[And @@ (Chop[#,10^-20]===0 & /@ residClean)]}
```
**Output:**
```
Out= {"weights tested: 13", "sum branches = dim(4x4): True", "mult==nBranch everywhere: True", "all residuals vanish: True"}
```
## 2026-07-23 17:29:50 — 📥 BULK INSERT 2 CELL(S) at positions 41–42
1. [markdown] ## Validation↵↵`fund ⊗ antifund`: for every distinct global weight, the number of Baxter solution br
2. [code] phiTest = {1, 1/3, 1/5, 1/7};↵wts = FindWeightDistinct[fund, antifund];↵validation = Table[Module[{s


## 2026-07-23 17:29:52
**Input:**
```wolfram
phiTest = {1, 1/3, 1/5, 1/7};
wts = FindWeightDistinct[fund, antifund];
validation = Table[Module[{s = SolveBax[fund, antifund, wt, phiTest], res},
    res = If[MatchQ[s, _List] && s =!= {},
       Max[Table[Max[Abs[N[CoefficientList[BaxEqn[a][fund, antifund, wt, phiTest][z] /. s[[b]], z], 30]]],
          {a, 1, 4}, {b, 1, Length[s]}]], 0] /. -Infinity -> 0;  (* -Infinity = identically-zero residual polynomial *)
    {wt, weightMult[fund, antifund, wt], Length[s], Chop[res, 10^-20]}], {wt, wts}];
Column[{
  Grid[Prepend[validation, {"weight", "mult", "nBranch", "maxResidual"}], Frame -> All],
  "Σ branches = dim(4×4)=16 : " <> ToString[Total[validation[[All, 3]]] == 16],
  "mult == nBranch ∀ : " <> ToString[And @@ (#[[2]] === #[[3]] & /@ validation)],
  "all residuals vanish : " <> ToString[And @@ (Chop[#[[4]], 10^-20] === 0 & /@ validation)]}]
```
**Output:**
```
2 kernel messages: N::meprec: Internal precision limit $MaxExtraPrecision = 50. reached while evaluating {630 - -----------------------------------------------------------------------------------------------------------… [also: N::meprec] result: Out[184]= \begin{matrix}\begin{array}{|c|c|c|c|}\hline\text{weight} & \text{mult} & \text{nBranch} & \text{maxResidual}\\\hline\left\{2,1,1,0\right\} & 1 & 1 & 0\\\hline\left\{2,1,0,1\right\} & 1 & 1 …
```

## 2026-07-23 17:30:14
**Input:**
```wolfram
(* build plain unshifted oper coeffs H[a][z] as explicit functions of z, with solved h's *)
checkFull[rep1_, rep2_, repall_, phis_, sol_] := Module[{Hc, Qfull, oper, resid},
  Hc = Table[(H[a][zz] //. expandH /. applyleadingH /. applycentralH /.
       applyrep[rep1, rep2, repall] /. Thread[(\[Phi] /@ Range[4]) -> phis]) /. sol, {a, 1, 4}];
  Table[
   Qfull = F[rep1, rep2][zz] Exp[-phis[[a]] zz] qAns[a][rep1, rep2, repall][zz] /. sol[[br]] /. z -> zz;
   oper = D[Qfull, {zz, 4}] + Sum[(Hc[[k]][[br]]) D[Qfull, {zz, 4 - k}], {k, 1, 4}];
   (oper/(F[rep1, rep2][zz] Exp[-phis[[a]] zz]) // Together // Numerator // Expand),
   {a, 1, 4}, {br, 1, Length[sol]}]];
(* adjoint: single branch, exact *)
resAdj = checkFull[fund, antifund, adj, phiTest, solAdj];
Max[Abs[Flatten[resAdj /. z -> 17/13]]]
```
**Output:**
```
Out= 0
```

## 2026-07-23 17:30:21
**Input:**
```wolfram
solM4 = SolveBax[fund, antifund, {1,1,1,1}, phiTest];
resM4 = checkFull[fund, antifund, {1,1,1,1}, phiTest, solM4];
{Length[solM4], Chop[Max[Abs[N[Flatten[resM4 /. z -> 17/13], 40]]], 10^-25]}
```
**Output:**
```
[message] N::meprec: Internal precision limit $MaxExtraPrecision = 50. reached while evaluating {-17146080 - 9185400 zz - 1825740 zz  - 69768216 Root[945 + 4104 #1 + 792 #1  - 344 #1  + 16 #1  & , 1] - 37736280… [16275 chars]
Out= {4, 0}
```
## 2026-07-23 17:30:39 — 📥 BULK INSERT 2 CELL(S) at positions 43–44
1. [markdown] ### Independent cross-check↵↵`SolveBax` works with the conjugated/shifted coefficients `HSc`. As an 
2. [code] ClearAll[checkFull];↵checkFull[rep1_, rep2_, repall_, phis_, sol_] := Module[{Hc, Qfull, oper},↵  Hc


## 2026-07-23 17:30:40
**Input:**
```wolfram
ClearAll[checkFull];
checkFull[rep1_, rep2_, repall_, phis_, sol_] := Module[{Hc, Qfull, oper},
  Hc = Table[(H[a][zz] //. expandH /. applyleadingH /. applycentralH /.
        applyrep[rep1, rep2, repall] /. Thread[(ϕ /@ Range[4]) -> phis]) /. sol, {a, 1, 4}];
  Table[
   Qfull = F[rep1, rep2][zz] Exp[-phis[[a]] zz] qAns[a][rep1, rep2, repall][zz] /. sol[[br]] /. z -> zz;
   oper = D[Qfull, {zz, 4}] + Sum[Hc[[k]][[br]] D[Qfull, {zz, 4 - k}], {k, 1, 4}];
   oper/(F[rep1, rep2][zz] Exp[-phis[[a]] zz]) // Together // Numerator // Expand,
   {a, 1, 4}, {br, 1, Length[sol]}]];

(* run it on both a mult-1 (adjoint) and a mult-4 weight; residuals at a random rational z *)
With[{r1 = checkFull[fund, antifund, adj, phiTest, SolveBax[fund, antifund, adj, phiTest]],
      r4 = checkFull[fund, antifund, {1,1,1,1}, phiTest, SolveBax[fund, antifund, {1,1,1,1}, phiTest]]},
 {"adjoint (mult 1) residual : " <> ToString[Chop[Max[Abs[N[Flatten[r1 /. z -> 17/13], 40]]], 10^-25]],
  "{1,1,1,1} (mult 4) residual : " <> ToString[Chop[Max[Abs[N[Flatten[r4 /. z -> 17/13], 40]]], 10^-25]]}]
```
**Output:**
```
2 kernel messages: N::meprec: Internal precision limit $MaxExtraPrecision = 50. reached while evaluating {-17146080 - 9185400 zz - 1825740 zz  - 69768216 Root[945 + 4104 #1 + 792 #1  - 344 #1  + 16 #1  & , 1] - 37736280… [also: N::meprec] result: Out[187]= \left\{\text{adjoint (mult 1) residual : 0},\text{\{1,1,1,1\} (mult 4) residual : 0}\right\}
```
## 2026-07-23 19:07:09 — ✏️ EDIT CELL 33
(*monic polynomial part q_a of the a-th Q-function; the F and e^{-\[Phi]_a z} factors
  are handled by the operator (HSc absorbs F; the shift \[PartialD]->\[PartialD]-\[Phi]_a absorbs the exp).
  HSc[…


## 2026-07-23 19:07:09
**Input:**
```wolfram
(*monic polynomial part q_a of the a-th Q-function; the F and e^{-\[Phi]_a z} factors
  are handled by the operator (HSc absorbs F; the shift \[PartialD]->\[PartialD]-\[Phi]_a absorbs the exp).
  HSc[0] evaluates to 1 automatically, so the k=0 term is the leading \[PartialD]^4.*)
ClearAll[qAns, BaxEqn];
qAns[a_][rep1_List, rep2_List, repall_List][z_] :=
  z^M[a][rep1, rep2, repall] + Sum[c[a, n] z^n, {n, 0, M[a][rep1, rep2, repall] - 1}];

BaxEqn[a_][rep1_List, rep2_List, repall_List, phis_List][z_] :=
  Module[{q = qAns[a][rep1, rep2, repall][z]},
   Sum[(HSc[k][rep1, rep2, repall][z] /. Thread[(ϕ /@ Range[4]) -> phis]) *
       Nest[D[#, z] - phis[[a]] # &, q, 4 - k], {k, 0, 4}] // Together // Numerator // Expand];
```
**Output:**
```
Out= Null
```
## 2026-07-23 19:07:14 — ✏️ EDIT CELL 34
(*Joint solver: all four Baxter equations are polynomial in z; Solve handles the
  full system. Returns one rule-list per Gaudin eigenstate with global weight repall.*)
ClearAll[SolveBax];
SolveBax[re…


## 2026-07-23 19:07:14
**Input:**
```wolfram
(*Joint solver: all four Baxter equations are polynomial in z; Solve handles the
  full system. Returns one rule-list per Gaudin eigenstate with global weight repall.*)
ClearAll[SolveBax];
SolveBax[rep1_List, rep2_List, repall_List, phis_List] :=
  Solve[
   Flatten[Table[
     Thread[CoefficientList[BaxEqn[a][rep1, rep2, repall, phis][z], z] == 0], {a, 1, 4}]],
   Join[Flatten[Table[c[a, n], {a, 1, 4}, {n, 0, M[a][rep1, rep2, repall] - 1}]], hUnknowns]];
```
**Output:**
```
Out= Null
```

## 2026-07-23 19:07:20
**Input:**
```wolfram
phiTest = {1, 1/3, 1/5, 1/7};
sA = SolveBax[fund, antifund, adj, phiTest];
s4 = SolveBax[fund, antifund, {1,1,1,1}, phiTest];
{Length[sA], Max[Abs[Table[BaxEqn[a][fund,antifund,adj,phiTest][z] /. sA[[1]], {a,4}]]],
 Length[s4]}
```
**Output:**
```
Out= {1, 0, 4}
```
## 2026-07-23 19:14:49 — ✏️ EDIT CELL 34
(*monic polynomial part q_a of the a-th Q-function; the F and e^{-\[Phi]_a z} factors
  are handled by the operator (HSc absorbs F; the shift \[PartialD]->\[PartialD]-\[Phi]_a absorbs the exp).
  HSc[…

## 2026-07-23 19:14:55 — ✏️ EDIT CELL 35
(*Joint solver: all four Baxter equations are polynomial in z; Solve handles the
  full system. Returns one rule-list per Gaudin eigenstate with global weight repall.
  Twists are fixed numerically vi…

## 2026-07-23 19:14:59 — ✏️ EDIT CELL 37
SolveBax[{1,1,0,0},{1,1,0,0},{2,2,0,0}]//Simplify

## 2026-07-23 19:15:04 — ▶️ RUN CELLS 33–37
Cell 33: ✓ — (no output)
Cell 34: ✓ — (no output)
Cell 35: ✓ — (no output)
Cell 36: ✓ — (no output)
Cell 37: ✓ — Out[209]= \left\{\left\{c\left(3,0\right)\to \frac{1145}{8},c\left(3,1\right)\to -\frac{37}{2},c\left(4,0\right)\to \frac{5509}{72},c\left(4,1\right)\to -\frac{83}{6},h\left(2,1,0\right)\to \frac{2}{105},h\left(3,1,0\right)\to -\frac{6}{35},h\left(3,2,0\right)\to \frac{68}{105},h\left(3,0,2\right)\to -\frac{772}{105},h\left(4,1,0\right)\to -\frac{2}{105},h\left(4,2,0\right)\to \frac{8}{35},h\left(… [516 chars total]


## 2026-07-23 19:15:11
**Input:**
```wolfram
{Length[SolveBax[fund, antifund, adj]],
 Max[Abs[Table[BaxEqn[a][fund,antifund,adj][z] /. SolveBax[fund,antifund,adj][[1]], {a,4}]]],
 Length[SolveBax[fund, antifund, {1,1,1,1}]]}
```
**Output:**
```
Out= {1, 0, 4}
```

---

## 2026-07-23 19:17:10 — 🔄 KERNEL RESTART

## 2026-07-23 19:21:36 — 📥 BULK INSERT 1 CELL(S) at positions 36–36
1. [code] (*Full Q-functions from the solved coefficients. Per Gaudin eigenstate (branch of↵  SolveBax) return

## 2026-07-23 19:21:54 — ▶️ RUN CELLS 1–38
Cell 12: ✓ — Syntax::sntxf: "the weights on each site*" cannot be followed by ")". | Syntax::sntxf: "the weights on each site*" cannot be followed by ")".
  ⚠️ SYNTAX MESSAGE DETECTED — definitions loaded before this error may be stale.


## 2026-07-23 19:22:01
**Input:**
```wolfram
BuildQs[fund, antifund, adj]
```
**Output:**
```
Out= {{E^(-z), (-3/2 + z)/E^(z/3), (-5/4 + z)/E^(z/5), (2191/72 - (313*z)/12 + z^2)/E^(z/7)}}
```

## 2026-07-23 19:22:20
**Input:**
```wolfram
sol = SolveBax[fund, antifund, adj][[1]];
Qs = BuildQs[fund, antifund, adj][[1]];
Table[
  (D[Qs[[a]], {z, 4}] + Sum[(HSc[k][fund, antifund, adj][z] /. sol /. numphis) D[Qs[[a]], {z, 4 - k}], {k, 1, 4}]) //
   Together // Simplify,
  {a, 4}]
```
**Output:**
```
Out= {0, 0, 0, 0}
```

## 2026-07-23 19:25:39
**Input:**
```wolfram
rep = {1,1,0,0};
sols = SolveBax[rep, rep, {1,1,1,1}];
{"weightMult = " <> ToString[weightMult[rep, rep, {1,1,1,1}]],
 "nBranches = " <> ToString[Length[sols]],
 "M[a] = " <> ToString[Table[M[a][rep,rep,{1,1,1,1}], {a,4}]],
 "distinct h[2,1,0] across branches: " <> ToString[Sort[h[2,1,0] /. sols]],
 "all residuals 0: " <> ToString[
   And @@ Flatten[Table[BaxEqn[a][rep,rep,{1,1,1,1}][z] /. sols[[b]] // Expand,
      {a,4},{b,Length[sols]}] /. 0 -> True]]}
```
**Output:**
```
Out= {"weightMult = 6", "nBranches = 6", "M[a] = {1, 1, 1, 1}", "distinct h[2,1,0] across branches:                                                   2           3          4         5        6                                                 2           3          4         5        6       2                                          2           3          4         5        6       3                                         2           3          4         5        6       4                                       2           3          4         5        6       5                                                   2           3          4         5        6                                                 2           3          4         5        6       2                                          2           3          4         5        6       3                                         2           3          4         5        6       4                                       2           3          4         5        6       5                                                   2           3          4         5        6                                                 2           3          4         5        6       2                                          2           3          4         5        6       3                                         2           3          4         5        6       4                                       2           3          4         5        6       5                                                   2           3          4         5        6                                                 2           3          4         5        6       2                                          2           3          4         5        6       3                                         2           3          4         5        6       4                                       2           3          4         5        6       5                                                   2           3          4         5        6                                                 2           3          4         5        6       2                                          2           3          4         5        6       3                                         2           3          4         5        6       4                                       2           3          4         5        6       5                                                   2           3          4         5        6                                                 2           3          4         5        6       2                                          2           3          4         5        6       3                                         2           3          4         5        6       4                                       2           3          4         5        6       5\n -89127 + 36711 Root[102087 + 176823 #1 - 73908 #1  - 21296 #1  + 8292 #1  - 672 #1  + 16 #1  & , 1] + 9446 Root[102087 + 176823 #1 - 73908 #1  - 21296 #1  + 8292 #1  - 672 #1  + 16 #1  & , 1]  - 4156 Root[102087 + 176823 #1 - 73908 #1  - 21296 #1  + 8292 #1  - 672 #1  + 16 #1  & , 1]  + 340 Root[102087 + 176823 #1 - 73908 #1  - 21296 #1  + 8292 #1  - 672 #1  + 16 #1  & , 1]  - 8 Root[102087 + 176823 #1 - 73908 #1  - 21296 #1  + 8292 #1  - 672 #1  + 16 #1  & , 1]   -89127 + 36711 Root[102087 + 176823 #1 - 73908 #1  - 21296 #1  + 8292 #1  - 672 #1  + 16 #1  & , 2] + 9446 Root[102087 + 176823 #1 - 73908 #1  - 21296 #1  + 8292 #1  - 672 #1  + 16 #1  & , 2]  - 4156 Root[102087 + 176823 #1 - 73908 #1  - 21296 #1  + 8292 #1  - 672 #1  + 16 #1  & , 2]  + 340 Root[102087 + 176823 #1 - 73908 #1  - 21296 #1  + 8292 #1  - 672 #1  + 16 #1  & , 2]  - 8 Root[102087 + 176823 #1 - 73908 #1  - 21296 #1  + 8292 #1  - 672 #1  + 16 #1  & , 2]   -89127 + 36711 Root[102087 + 176823 #1 - 73908 #1  - 21296 #1  + 8292 #1  - 672 #1  + 16 #1  & , 3] + 9446 Root[102087 + 176823 #1 - 73908 #1  - 21296 #1  + 8292 #1  - 672 #1  + 16 #1  & , 3]  - 4156 Root[102087 + 176823 #1 - 73908 #1  - 21296 #1 
[output truncated — 328130 chars total]
```

## 2026-07-23 19:25:47
**Input:**
```wolfram
N[Sort[h[2,1,0] /. sols], 8]
```
**Output:**
```
Out= {-1.4510512344447966228`8., -5.5244749720648766174`8., 0.5681501314906029865`8., -1.5387865434569617331`8., -1.5789701652464906992`8., 0.6108470694368084002`8.}
```

---

## 2026-07-24 07:24:56 — 🔄 KERNEL RESTART


---

## 2026-07-24 07:40:53 — 🔄 KERNEL RESTART


## 2026-07-24 07:47:36
**Input:**
```wolfram
Coefficient[detOp, f''[z]] // Short
```
**Output:**
```
Out= L[1, 1][z] . L[2, 2][z]   L[1, 1][z] . L[3, 3][z]   L[1, 1][z] . L[2, 2][z] + L[1, 1][z] . L[3, 3][z]   L[1, 1][z] . L[4, 4][z]   L[1, 1][z] . L[2, 2][z] + L[1, 1][z] . L[4, 4][z]   L[1, 1][z] . L[3, 3][z] + L[1, 1][z] . L[4, 4][z]   L[1, 1][z] . L[2, 2][z] + L[1, 1][z] . L[3, 3][z] + L[1, 1][z] . L[4, 4][z]   L[1, 2][z] . L[2, 1][z]   L[1, 3][z] . L[3, 1][z]   L[1, 4][z] . L[4, 1][z]   L[2, 1][z] . L[1, 2][z]   L[2, 2][z] . L[1, 1][z]   L[2, 2][z] . L[3, 3][z]   L[2, 2][z] . L[1, 1][z] + L[2, 2][z] . L[3, 3][z]   L[2, 2][z] . L[4, 4][z]   L[2, 2][z] . L[1, 1][z] + L[2, 2][z] . L[4, 4][z]   L[2, 2][z] . L[3, 3][z] + L[2, 2][z] . L[4, 4][z]   L[2, 2][z] . L[1, 1][z] + L[2, 2][z] . L[3, 3][z] + L[2, 2][z] . L[4, 4][z]   L[2, 3][z] . L[3, 2][z]   L[2, 4][z] . L[4, 2][z]   L[3, 1][z] . L[1, 3][z]   L[3, 2][z] . L[2, 3][z]   L[3, 3][z] . L[1, 1][z]   L[3, 3][z] . L[2, 2][z]   L[3, 3][z] . L[1, 1][z] + L[3, 3][z] . L[2, 2][z]   L[3, 3][z] . L[4, 4][z]   L[3, 3][z] . L[1, 1][z] + L[3, 3][z] . L[4, 4][z]   L[3, 3][z] . L[2, 2][z] + L[3, 3][z] . L[4, 4][z]   L[3, 3][z] . L[1, 1][z] + L[3, 3][z] . L[2, 2][z] + L[3, 3][z] . L[4, 4][z]   L[3, 4][z] . L[4, 3][z]   L[4, 1][z] . L[1, 4][z]   L[4, 2][z] . L[2, 4][z]   L[4, 3][z] . L[3, 4][z]   L[4, 4][z] . L[1, 1][z]   L[4, 4][z] . L[2, 2][z]   L[4, 4][z] . L[1, 1][z] + L[4, 4][z] . L[2, 2][z]   L[4, 4][z] . L[3, 3][z]   L[4, 4][z] . L[1, 1][z] + L[4, 4][z] . L[3, 3][z]   L[4, 4][z] . L[2, 2][z] + L[4, 4][z] . L[3, 3][z]   L[4, 4][z] . L[1, 1][z] + L[4, 4][z] . L[2, 2][z] + L[4, 4][z] . L[3, 3][z]   3 (L[1, 1])'[z]   3 (L[2, 2])'[z]   3 (L[3, 3])'[z]   3 (L[4, 4])'[z]
----------------------- + ----------------------- + ------------------------------------------------- + ----------------------- + ------------------------------------------------- + ------------------------------------------------- + --------------------------------------------------------------------------- - ----------------------- - ----------------------- - ----------------------- - ----------------------- + ----------------------- + ----------------------- + ------------------------------------------------- + ----------------------- + ------------------------------------------------- + ------------------------------------------------- + --------------------------------------------------------------------------- - ----------------------- - ----------------------- - ----------------------- - ----------------------- + ----------------------- + ----------------------- + ------------------------------------------------- + ----------------------- + ------------------------------------------------- + ------------------------------------------------- + --------------------------------------------------------------------------- - ----------------------- - ----------------------- - ----------------------- - ----------------------- + ----------------------- + ----------------------- + ------------------------------------------------- + ----------------------- + ------------------------------------------------- + ------------------------------------------------- + --------------------------------------------------------------------------- + --------------- + --------------- + --------------- + ---------------
          12                        12                                     12                                     12                                     12                                                  12                                                                4                                                   2                         2                         2                         2                        12                        12                                     12                                     12                                     12                                                  12                                                                4                                                   2                         2                         2                         2                  
[output truncated — 4882 chars total]
```

## 2026-07-24 07:47:43
**Input:**
```wolfram
Derivative[1][L[1,1]][z] /. L -> LL[{1,1,0,0},{1,1,0,0}] // Normal // Short
```
**Output:**
```
Out=                                        -2                                -2                                                                                                                                                     -2                                -2                                                                                                                                                     -2                                -2                                                                                                                                              -2                                                                                                                                              -2                                                                                                                                              -2                                                                                                                                                     -2                                -2                                                                                                                                                     -2                                -2                                                                                                                                                     -2                                -2                                                                                                                                              -2                                                                                                                                              -2                                                                                                                                              -2                                                                                                                                                     -2                                -2                                                                                                                                                     -2                                -2                                                                                                                                                     -2                                -2                                                                                                                                              -2                                                                                                                                              -2                                                                                                                                              -2                                                                                                                                                     -2                                                                                                                                                     -2                                                                                                                                                     -2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        -2                                                                                                                                                     -2                                                                                                                                                     -2
[output truncated — 10118 chars total]
```

## 2026-07-24 07:47:49
**Input:**
```wolfram
flawed = Derivative[1][L[1,1]][z] /. L -> LL[{1,1,0,0},{1,1,0,0}];
robust = D[LL[{1,1,0,0},{1,1,0,0}][1,1][z], z];
{leakInFlawed, leakInRobust} = {!FreeQ[flawed, System`Private`DerivativeX], !FreeQ[robust, System`Private`DerivativeX]}
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= {False, False}
```

## 2026-07-24 07:47:57
**Input:**
```wolfram
res = Hop[2][{1,1,0,0},{1,1,0,0}][z];
{Head[res], !FreeQ[res, System`Private`DerivativeX], Count[res, System`Private`DerivativeX[_], Infinity, Heads->True]}
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= {SparseArray, False, 0}
```

---

## 2026-07-24 07:49:21 — 🔄 KERNEL RESTART

## 2026-07-24 07:49:43 — ▶️ RUN CELLS 1–32
## 2026-07-24 07:49:48 — ▶️ RUN CELLS 40–41
Cell 40: ✓ — Out[53]= \left\{\frac{2}{z}+\frac{2}{z-1}+\frac{176}{105},-\frac{4}{z^2}+\frac{2}{105\,z}+\frac{422}{105\,\left(z-1\right)}-\frac{4}{\left(z-1\right)^2}+\frac{86}{105},\frac{4}{z^3}+\frac{68}{105\,z^2}-\frac{6}{35\,z}+\frac{6}{5\,\left(z-1\right)}-\frac{772}{105\,\left(z-1\right)^2}+\frac{4}{\left(z-1\right)^3}+\frac{16}{105},-\frac{4}{3\,z^3}+\frac{8}{35\,z^2}-\frac{2}{105\,z}+\frac{2}{21\,\left(… [498 chars total]
Cell 41: ✓ — Out[54]= \left(\mathrm{Hop}\left(2\right)\left(\left\{1,1,0,0\right\},\left\{1,1,0,0\right\}\right)\right)\left(z\right)

## 2026-07-24 07:50:02 — ▶️ RUN CELLS 36–41
Cell 36: ✓ — (no output)
Cell 37: ✓ — (no output)
Cell 38: ✓ — (no output)
Cell 39: ✓ — (no output)
Cell 40: ✓ — Out[65]= \left\{\frac{2}{z}+\frac{2}{z-1}+\frac{176}{105},-\frac{4}{z^2}+\frac{2}{105\,z}+\frac{422}{105\,\left(z-1\right)}-\frac{4}{\left(z-1\right)^2}+\frac{86}{105},\frac{4}{z^3}+\frac{68}{105\,z^2}-\frac{6}{35\,z}+\frac{6}{5\,\left(z-1\right)}-\frac{772}{105\,\left(z-1\right)^2}+\frac{4}{\left(z-1\right)^3}+\frac{16}{105},-\frac{4}{3\,z^3}+\frac{8}{35\,z^2}-\frac{2}{105\,z}+\frac{2}{21\,\left(… [498 chars total]
Cell 41: ✓ — Out[66]= \begin{pmatrix}2\,\left(-\frac{1}{z^2}-\frac{1}{\left(z-1\right)^2}\right)+3\,\left(-\frac{1}{\text{System$\grave{ }$Private$\grave{ }$DerivativeX}\left(1\right)^2}-\frac{1}{\left(\text{System$\grave{ }$Private$\grave{ }$DerivativeX}\left(1\right)-1\right)^2}\right)+\frac{1}{6}\,\left(\phi \left(1\right)+\frac{1}{z-1}+\frac{1}{z}\right)\,\left(\phi \left(2\right)+\frac{1}{z-1}+\frac{1}{z}… [1248 chars total]


## 2026-07-24 07:50:14
**Input:**
```wolfram
t1 = Derivative[1][L[1,1]][z] /. L -> LL[{1,1,0,0},{1,1,0,0}];
leak1 = !FreeQ[t1, System`Private`DerivativeX];
t2 = D[LL[{1,1,0,0},{1,1,0,0}][1,1][z], z];
leak2 = !FreeQ[t2, System`Private`DerivativeX];
{leak1, leak2}
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= (suppressed)
Out[5]= {False, False}
```

## 2026-07-24 07:51:24
**Input:**
```wolfram
applyLL[rep1_,rep2_][expr_] := expr /. {
   Derivative[n_][L[a_,b_]][zz_] :> D[LL[rep1,rep2][a,b][zz], {zz, n}],
   L[a_,b_][zz_] :> LL[rep1,rep2][a,b][zz]};
fixed = applyLL[{1,1,0,0},{1,1,0,0}][Coefficient[detOp, f''[z]]] // ExpandAll;
{Head[fixed], !FreeQ[fixed, System`Private`DerivativeX], !FreeQ[fixed, L]}
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= {SparseArray, False, False}
```

## 2026-07-24 07:52:26
**Input:**
```wolfram
(* fixed Hop for all 4, using explicit differentiation *)
derivOf = {f'''[z], f''[z], f'[z], f[z]};
HopFixed[k_][r1_,r2_][zz_] := applyLL[r1,r2][Coefficient[detOp, derivOf[[k]]]] /. z->zz // ExpandAll;
(* energies from Baxter/analytic side *)
ens = ComputeEnergy[{1,1,0,0},{1,1,0,0},{2,2,0,0}];
zval = 2;
enNum = (ens /. numphis /. z->zval);
(* operator matrices at z=2 with numeric twists *)
mats = Table[Normal[HopFixed[k][{1,1,0,0},{1,1,0,0}][zval]] /. numphis, {k,1,4}];
(* check each analytic energy is an eigenvalue of the corresponding operator matrix *)
Table[MemberQ[Chop[N[Eigenvalues[mats[[k]]]]], _?(Abs[#-N[enNum[[k]]]]<10^-8 &)], {k,1,4}]
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= (suppressed)
Out[5]= (suppressed)
Out[6]= (suppressed)
Out[7]= {False, False, False, False}
```

## 2026-07-24 07:52:52
**Input:**
```wolfram
Column[{
  "enNum (analytic H_k at z=2, numeric twist):",
  N[enNum],
  "Hop[1] eigenvalues at z=2:",
  Sort[Chop[N[Eigenvalues[mats[[1]]]]]],
  "Hop[2] eigenvalues at z=2:",
  Sort[Chop[N[Eigenvalues[mats[[2]]]]]]
}]
```
**Output:**
```
Out= Column[{"enNum (analytic H_k at z=2, numeric twist):", {4.6761904761904765, -0.1523809523809524, -1.4238095238095239, 5.5095238095238095}, "Hop[1] eigenvalues at z=2:", {(0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z), (0.01904761904761905*(-105. + 122.*z + 88.*z^2))/((-1. + z)*z)}, "Hop[2] eigenvalues at z=2:", {(0.01904761904761905*(-210. + 421.*z - 590.*z^2 + 126.*z^3 + 43.*z^4))/((-1. + z)^2*z^2), (0.01904761904761905*(-210. + 428.*z - 611.*z^2 + 140.*z^3 + 43.*z^4))/((-1. + z)^2*z^2), (0.01904761904761905*(-210. + 431.*z - 620.*z^2 + 146.*z^3 + 43.*z^4))/((-1. + z)^2*z^2), (0.01904761904761905*(-210. + 463.*z - 716.*z^2 + 210.*z^3 + 43.*z^4))/((-1. + z)^2*z^2), (0.01904761904761905*(-210. + 466.*z - 725.*z^2 + 216.*z^3 + 43.*z^4))/((-1. + z)^2*z^2), (0.01904761904761905*(-210. + 473.*z - 746.*z^2 + 230.*z^3 + 43.*z^4))/((-1. + z)^2*z^2), (0.009523809523809525*(-420. + 782.*z - 1210.*z^2 + 342.*z^3 + 86.*z^4 - 114.23659658795863*Sqrt[z^2 - 2.*z^3 + z^4]))/((-1. + z)^2*z^2), (0.009523809523809525*(-420. + 796.*z - 1252.*z^2 + 370.*z^3 + 86.*z^4 - 114.23659658795863*Sqrt[z^2 - 2.*z^3 + z^4]))/((-1. + z)^2*z^2), (0.009523809523809525*(-420. + 779.*z - 1201.*z^2 + 336.*z^3 + 86.*z^4 - 113.08846094982458*Sqrt[z^2 - 2.*z^3 + z^4]))/((-1. + z)^2*z^2), (0.009523809523809525*(-420. + 799.*z - 1261.*z^2 + 376.*z^3 + 86.*z^4 - 113.08846094982458*Sqrt[z^2 - 2.*z^3 + z^4]))/((-1. + z)^2*z^2), (0.009523809523809525*(-420. + 786.*z - 1222.*z^2 + 350.*z^3 + 86.*z^4 - 110.67971810589329*Sqrt[z^2 - 2.*z^3 + z^4]))/((-1. + z)^2*z^2), (0.009523809523809525*(-420. + 792.*z - 1240.*z^2 + 362.*z^3 + 86.*z^4 - 110.67971810589329*Sqrt[z^2 - 2.*z^3 + z^4]))/((-1. + z)^2*z^2), (0.009523809523809525*(-420. + 747.*z - 1105.*z^2 + 272.*z^3 + 86.*z^4 - 105.47511554864494*Sqrt[z^2 - 2.*z^3 + z^4]))/((-1. + z)^2*z^2), (0.009523809523809525*(-420. + 831.*z - 1357.*z^2 + 440.*z^3 + 86.*z^4 - 105.4751155
[output truncated — 18605 chars total]
```

## 2026-07-24 07:53:39
**Input:**
```wolfram
matsNum = Table[ Normal[HopFixed[k][{1,1,0,0},{1,1,0,0}][z]] /. numphis /. z->2 // N, {k,1,4}];
enN = N[enNum];
Table[{k, enN[[k]], AnyTrue[Eigenvalues[matsNum[[k]]], Abs[#-enN[[k]]]<10^-7&]}, {k,1,4}]
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= {{1, 4.6761904761904765, True}, {2, -0.1523809523809524, True}, {3, -1.4238095238095239, True}, {4, 5.5095238095238095, True}}
```

## 2026-07-24 07:53:50
**Input:**
```wolfram
checkCommon[zv_] := Module[{ms, es, evec, resid},
  ms = Table[Normal[HopFixed[k][{1,1,0,0},{1,1,0,0}][z]] /. numphis /. z->zv // N, {k,1,4}];
  es = N[ens /. numphis /. z->zv];
  (* find a vector that is simultaneously an eigenvector of all 4 with eigenvalue es[[k]] *)
  evec = NullSpace[ms[[2]] - es[[2]] IdentityMatrix[36]];
  resid = If[evec==={}, "no eigvec for H2",
    Table[Max@Abs[(ms[[k]] - es[[k]] IdentityMatrix[36]).First[evec]], {k,1,4}]];
  {Length[evec], resid}];
{checkCommon[2], checkCommon[7/3]}
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= {{1, {0., 0., 0., 0.}}, {1, {0., 0., 0., 0.}}}
```
## 2026-07-24 07:54:25 — ✏️ EDIT CELL 39
ClearAll[Hop, applyLL];
(* Substitute the symbolic Lax matrix L -> LL[rep1,rep2].  Derivatives of L
   are handled EXPLICITLY: (L[a,b])^(n)[z] -> D[LL[..][a,b][z],{z,n}], i.e. we
   differentiate the …


## 2026-07-24 07:54:25
**Input:**
```wolfram
ClearAll[Hop, applyLL];
(* Substitute the symbolic Lax matrix L -> LL[rep1,rep2].  Derivatives of L
   are handled EXPLICITLY: (L[a,b])^(n)[z] -> D[LL[..][a,b][z],{z,n}], i.e. we
   differentiate the concrete matrix in the real variable z.  The old blanket
   rule  /. L -> LL[rep1,rep2]  turned (L[a,b])'[z] into
   Derivative[1][LL[rep1,rep2][a,b]][z] and let Mathematica differentiate the
   SubValue-defined LL through its internal machinery, which (order-dependently)
   leaked the formal variable System`Private`DerivativeX into the output. *)
applyLL[rep1_, rep2_][expr_] := expr /. {
    Derivative[n_][L[a_, b_]][zz_] :> D[LL[rep1, rep2][a, b][zz], {zz, n}],
    L[a_, b_][zz_] :> LL[rep1, rep2][a, b][zz]};
Hop[1][rep1_, rep2_][z_] := applyLL[rep1, rep2][Coefficient[detOp, f'''[z]]] // ExpandAll;
Hop[2][rep1_, rep2_][z_] := applyLL[rep1, rep2][Coefficient[detOp, f''[z]]]  // ExpandAll;
Hop[3][rep1_, rep2_][z_] := applyLL[rep1, rep2][Coefficient[detOp, f'[z]]]   // ExpandAll;
Hop[4][rep1_, rep2_][z_] := applyLL[rep1, rep2][Coefficient[detOp, f[z]]]    // ExpandAll;
```
**Output:**
```
Out= Null
```

---

## 2026-07-24 07:54:29 — 🔄 KERNEL RESTART

## 2026-07-24 07:54:48 — ▶️ RUN CELLS 1–39
## 2026-07-24 07:54:53 — ▶️ RUN CELLS 40–41
Cell 40: ✓ — Out[66]= \left\{\frac{2}{z}+\frac{2}{z-1}+\frac{176}{105},-\frac{4}{z^2}+\frac{2}{105\,z}+\frac{422}{105\,\left(z-1\right)}-\frac{4}{\left(z-1\right)^2}+\frac{86}{105},\frac{4}{z^3}+\frac{68}{105\,z^2}-\frac{6}{35\,z}+\frac{6}{5\,\left(z-1\right)}-\frac{772}{105\,\left(z-1\right)^2}+\frac{4}{\left(z-1\right)^3}+\frac{16}{105},-\frac{4}{3\,z^3}+\frac{8}{35\,z^2}-\frac{2}{105\,z}+\frac{2}{21\,\left(… [498 chars total]
Cell 41: ✓ — Out[67]= \begin{pmatrix}5\,\left(-\frac{1}{z^2}-\frac{1}{\left(z-1\right)^2}\right)+\frac{1}{6}\,\left(\phi \left(1\right)+\frac{1}{z-1}+\frac{1}{z}\right)\,\left(\phi \left(2\right)+\frac{1}{z-1}+\frac{1}{z}\right)+\frac{1}{6}\,\left(\phi \left(1\right)+\frac{1}{z-1}+\frac{1}{z}\right)\,\phi \left(3\right)+\frac{1}{6}\,\left(\phi \left(2\right)+\frac{1}{z-1}+\frac{1}{z}\right)\,\phi \left(3\right… [1248 chars total]


## 2026-07-24 07:55:01
**Input:**
```wolfram
r = Hop[2][{1,1,0,0},{1,1,0,0}][z];
{Head[r], !FreeQ[r, System`Private`DerivativeX], !FreeQ[r, L], Dimensions[r]}
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= {SparseArray, False, False, {36, 36}}
```

---

## 2026-07-24 08:00:05 — 🔄 KERNEL RESTART


## 2026-07-24 08:16:40
**Input:**
```wolfram
SolveBax[{1,1,0,0},{1,1,0,0},{2,2,0,0}]
```
**Output:**
```
Out=             1145               37              5509               83                  2                   6                  68                   772                    2                  8                   4                   8                 20
{c[3, 0] -> ----, c[3, 1] -> -(--), c[4, 0] -> ----, c[4, 1] -> -(--), h[2, 1, 0] -> ---, h[3, 1, 0] -> -(--), h[3, 2, 0] -> ---, h[3, 0, 2] -> -(---), h[4, 1, 0] -> -(---), h[4, 2, 0] -> --, h[4, 3, 0] -> -(-), h[4, 0, 2] -> -(-), h[4, 0, 3] -> --}
             8                 2                72                6                  105                  35                 105                  105                   105                 35                  3                   7                 3
```

## 2026-07-24 08:16:40
**Input:**
```wolfram
{Head[ComputeEnergy[{1,1,0,0},{1,1,0,0},{2,2,0,0}]], Length[ComputeEnergy[{1,1,0,0},{1,1,0,0},{2,2,0,0}]], Depth[ComputeEnergy[{1,1,0,0},{1,1,0,0},{2,2,0,0}]]}
```
**Output:**
```
Out= {List, 4, 6}
```

## 2026-07-24 08:17:55
**Input:**
```wolfram
Short[Coefficient[detOp, f'''[z]], 10]
```
**Output:**
```
Out= Short[L[1, 1][z] + L[2, 2][z] + L[3, 3][z] + L[4, 4][z], 10]
```

## 2026-07-24 08:17:55
**Input:**
```wolfram
{Length[Coefficient[detOp, f[z]]//If[Head[#]===Plus,List@@#,{#}]&], Length[Coefficient[detOp, f'[z]]//If[Head[#]===Plus,List@@#,{#}]&], Length[Coefficient[detOp, f''[z]]//If[Head[#]===Plus,List@@#,{#}]&], Length[Coefficient[detOp, f'''[z]]//If[Head[#]===Plus,List@@#,{#}]&]}
```
**Output:**
```
Out= {1036, 280, 44, 4}
```

## 2026-07-24 08:18:09
**Input:**
```wolfram
Timing[hop1big = Hop[1][{3,3,0,0},{5,5,0,0}][z];]
```
**Output:**
```
Out= {0.092392, Null}
```

## 2026-07-24 08:18:22
**Input:**
```wolfram
Timing[hop2big = Hop[2][{3,3,0,0},{5,5,0,0}][z];]
```
**Output:**
```
Out= {9.778583, Null}
```

## 2026-07-24 08:21:22
**Input:**
```wolfram
Timing[hop3big = Hop[3][{3,3,0,0},{5,5,0,0}][z];]
```
**Output:**
```
Out= {177.115344, Null}
```

## 2026-07-24 08:23:50
**Input:**
```wolfram
ClearAll[splitDotTerm, HopVec];
splitDotTerm[term_] := Module[{d, s}, Which[Head[term] === Dot, {1, List @@ term}, Head[term] === Times, d = Cases[term, _Dot, {1}]; If[d === {}, {term, {1}}, s = term/First[d]; {s, List @@ First[d]}], True, {1, {term}}]];
HopVec[k_][rep1_, rep2_][z_][v_] := Module[{coeffExpr, terms, subst}, coeffExpr = Coefficient[detOp, D[f[z], {z, 4 - k}]]; terms = If[Head[coeffExpr] === Plus, List @@ coeffExpr, {coeffExpr}]; subst = {L[a_, b_][zz_] :> LL[rep1, rep2][a, b][zz], Derivative[n_][L[a_, b_]][zz_] :> D[LL[rep1, rep2][a, b][zz], {zz, n}]}; Total[(With[{parts = splitDotTerm[#]}, parts[[1]]*Fold[(#2 . #1) &, v, Reverse[parts[[2]] /. subst]]]) & /@ terms]];
```
**Output:**
```
Out= Null
```

## 2026-07-24 08:24:57
**Input:**
```wolfram
N36 = dim[{1,1,0,0}]*dim[{1,1,0,0}];
testv = RandomInteger[{-3,3}, N36];
Table[
  full = Hop[k][{1,1,0,0},{1,1,0,0}][z] . testv // Together;
  fast = HopVec[k][{1,1,0,0},{1,1,0,0}][z][testv] // Together;
  Max[Abs[(full - fast) /. z -> RandomReal[{2,5}]]],
{k,1,4}]
```
**Output:**
```
Timed out after 59s — kernel is still alive.
Simplify the expression, wrap it in TimeConstrained manually, or increase timeoutSeconds.
```

## 2026-07-24 08:25:02
**Input:**
```wolfram
N36 = dim[{1,1,0,0}]*dim[{1,1,0,0}]; testv = RandomInteger[{-3,3}, N36]; hop1s = Hop[1][{1,1,0,0},{1,1,0,0}][z]; hop2s = Hop[2][{1,1,0,0},{1,1,0,0}][z];
```
**Output:**
```
Out= Null
```

## 2026-07-24 08:25:06
**Input:**
```wolfram
hop3s = Hop[3][{1,1,0,0},{1,1,0,0}][z];
```
**Output:**
```
Out= Null
```

## 2026-07-24 08:25:06
**Input:**
```wolfram
hop4s = Hop[4][{1,1,0,0},{1,1,0,0}][z];
```
**Output:**
```
Out= Null
```

## 2026-07-24 08:25:12
**Input:**
```wolfram
fast1 = HopVec[1][{1,1,0,0},{1,1,0,0}][z][testv]; fast2 = HopVec[2][{1,1,0,0},{1,1,0,0}][z][testv]; fast3 = HopVec[3][{1,1,0,0},{1,1,0,0}][z][testv]; fast4 = HopVec[4][{1,1,0,0},{1,1,0,0}][z][testv];
```
**Output:**
```
Out= Null
```

## 2026-07-24 08:25:14
**Input:**
```wolfram
zval = 17/5;
{Max[Abs[N[(hop1s . testv - fast1) /. z->zval]]],
 Max[Abs[N[(hop2s . testv - fast2) /. z->zval]]],
 Max[Abs[N[(hop3s . testv - fast3) /. z->zval]]],
 Max[Abs[N[(hop4s . testv - fast4) /. z->zval]]]}
```
**Output:**
```
Out= {Max[0., Abs[3.*(0.29411764705882354 + ϕ[1.]) + 3.*(0.4166666666666667 + ϕ[2.]) + 3.*(0.7107843137254902 + ϕ[3.]) + 3.*ϕ[4.] - 3.*(1.4215686274509804 + ϕ[1.] + ϕ[2.] + ϕ[3.] + ϕ[4.])], Abs[3.*ϕ[1.] + 3.*(0.7107843137254902 + ϕ[2.]) + 3.*(0.7107843137254902 + ϕ[3.]) + 3.*ϕ[4.] - 3.*(1.4215686274509804 + ϕ[1.] + ϕ[2.] + ϕ[3.] + ϕ[4.])], Abs[3.*(0.7107843137254902 + ϕ[1.]) + 3.*(0.4166666666666667 + ϕ[2.]) + 3.*ϕ[3.] + 3.*(0.29411764705882354 + ϕ[4.]) - 3.*(1.4215686274509804 + ϕ[1.] + ϕ[2.] + ϕ[3.] + ϕ[4.])], Abs[3.*(0.7107843137254902 + ϕ[1.]) + 3.*(0.29411764705882354 + ϕ[2.]) + 3.*ϕ[3.] + 3.*(0.4166666666666667 + ϕ[4.]) - 3.*(1.4215686274509804 + ϕ[1.] + ϕ[2.] + ϕ[3.] + ϕ[4.])], Abs[3.*(0.29411764705882354 + ϕ[1.]) + 3.*(0.4166666666666667 + ϕ[2.]) + 3.*(0.29411764705882354 + ϕ[3.]) + 3.*(0.4166666666666667 + ϕ[4.]) - 3.*(1.4215686274509804 + ϕ[1.] + ϕ[2.] + ϕ[3.] + ϕ[4.])], Abs[3.*(0.4166666666666667 + ϕ[1.]) + 3.*(0.29411764705882354 + ϕ[2.]) + 3.*ϕ[3.] + 3.*(0.7107843137254902 + ϕ[4.]) - 3.*(1.4215686274509804 + ϕ[1.] + ϕ[2.] + ϕ[3.] + ϕ[4.])], Abs[2.*(0.7107843137254902 + ϕ[1.]) + 2.*(0.7107843137254902 + ϕ[2.]) + 2.*ϕ[3.] + 2.*ϕ[4.] - 2.*(1.4215686274509804 + ϕ[1.] + ϕ[2.] + ϕ[3.] + ϕ[4.])], Abs[2.*(0.4166666666666667 + ϕ[1.]) + 2.*(0.29411764705882354 + ϕ[2.]) + 2.*(0.7107843137254902 + ϕ[3.]) + 2.*ϕ[4.] - 2.*(1.4215686274509804 + ϕ[1.] + ϕ[2.] + ϕ[3.] + ϕ[4.])], Abs[2.*(0.29411764705882354 + ϕ[1.]) + 2.*(0.4166666666666667 + ϕ[2.]) + 2.*(0.4166666666666667 + ϕ[3.]) + 2.*(0.29411764705882354 + ϕ[4.]) - 2.*(1.4215686274509804 + ϕ[1.] + ϕ[2.] + ϕ[3.] + ϕ[4.])], Abs[2.*ϕ[1.] + 2.*(0.7107843137254902 + ϕ[2.]) + 2.*(0.29411764705882354 + ϕ[3.]) + 2.*(0.4166666666666667 + ϕ[4.]) - 2.*(1.4215686274509804 + ϕ[1.] + ϕ[2.] + ϕ[3.] + ϕ[4.])], Abs[-2.*(0.29411764705882354 + ϕ[1.]) - 2.*(0.7107843137254902 + ϕ[2.]) - 2.*(0.4166666666666667 + ϕ[3.]) - 2.*ϕ[4.] + 2.*(1.4215686274509804 + ϕ[1.] + ϕ[2.] + ϕ[3.] + ϕ[4.])], Abs[-2.*(0.4166666666666667 + ϕ[1.]) - 2.*(0.7107843137254902 + ϕ[2.]) - 2.*ϕ[3.] - 2.*(0.29411764705882354 + ϕ[4.]) + 2.*(1.4215686274509804 + ϕ[1.] + ϕ[2.] + ϕ[3.] + ϕ[4.])], Abs[-2.*ϕ[1.] - 2.*(0.7107843137254902 + ϕ[2.]) - 2.*(0.4166666666666667 + ϕ[3.]) - 2.*(0.29411764705882354 + ϕ[4.]) + 2.*(1.4215686274509804 + ϕ[1.] + ϕ[2.] + ϕ[3.] + ϕ[4.])], Abs[-2.*(0.4166666666666667 + ϕ[1.]) - 2.*ϕ[2.] - 2.*(0.7107843137254902 + ϕ[3.]) - 2.*(0.29411764705882354 + ϕ[4.]) + 2.*(1.4215686274509804 + ϕ[1.] + ϕ[2.] + ϕ[3.] + ϕ[4.])], Abs[-2.*ϕ[1.] - 2.*(0.4166666666666667 + ϕ[2.]) - 2.*(0.7107843137254902 + ϕ[3.]) - 2.*(0.29411764705882354 + ϕ[4.]) + 2.*(1.4215686274509804 + ϕ[1.] + ϕ[2.] + ϕ[3.] + ϕ[4.])], Abs[-3.*(0.7107843137254902 + ϕ[1.]) - 3.*(0.4166666666666667 + ϕ[2.]) - 3.*(0.29411764705882354 + ϕ[3.]) - 3.*ϕ[4.] + 3.*(1.4215686274509804 + ϕ[1.] + ϕ[2.] + ϕ[3.] + ϕ[4.])], Abs[-3.*(0.7107843137254902 + ϕ[1.]) - 3.*ϕ[2.] - 3.*(0.7107843137254902 + ϕ[3.]) - 3.*ϕ[4.] + 3.*(1.4215686274509804 + ϕ[1.] + ϕ[2.] + ϕ[3.] + ϕ[4.])], Abs[-3.*(0.4166666666666667 + ϕ[1.]) - 3.*(0.29411764705882354 + ϕ[2.]) - 3.*(0.4166666666666667 + ϕ[3.]) - 3.*(0.29411764705882354 + ϕ[4.]) + 3.*(1.4215686274509804 + ϕ[1.] + ϕ[2.] + ϕ[3.] + ϕ[4.])], Abs[-3.*(0.4166666666666667 + ϕ[1.]) - 3.*(0.29411764705882354 + ϕ[2.]) - 3.*(0.29411764705882354 + ϕ[3.]) - 3.*(0.4166666666666667 + ϕ[4.]) + 3.*(1.4215686274509804 + ϕ[1.] + ϕ[2.] + ϕ[3.] + ϕ[4.])], Abs[-3.*(0.29411764705882354 + ϕ[1.]) - 3.*(0.29411764705882354 + ϕ[2.]) - 3.*(0.4166666666666667 + ϕ[3.]) - 3.*(0.4166666666666667 + ϕ[4.]) + 3.*(1.4215686274509804 + ϕ[1.] + ϕ[2.] + ϕ[3.] + ϕ[4.])], Abs[-3.*ϕ[1.] - 3.*(0.7107843137254902 + ϕ[2.]) - 3.*ϕ[3.] - 3.*(0.7107843137254902 + ϕ[4.]) + 3.*(1.4215686274509804 + ϕ[1.] + ϕ[2.] + ϕ[3.] + ϕ[4.])], Abs[-3.*ϕ[1.] - 3.*(0.29411764705882354 + ϕ[2.]) - 3.*(0.4166666666666667 + ϕ[3.]) - 3.*(0.7107843137254902 + ϕ[4.]) + 3.*(1.4215686274509804 + ϕ[1.] + ϕ[2.] + ϕ[3.] + ϕ[4.])]], Max[Abs[-0.08333333333333333*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1
[output truncated — 8005079 chars total]
```

## 2026-07-24 08:25:22
**Input:**
```wolfram
zval = 17/5;
diffs = Table[
   (Hop[k][{1,1,0,0},{1,1,0,0}][z] . testv - HopVec[k][{1,1,0,0},{1,1,0,0}][z][testv]),
   {k, 1, 2}];
Max[Abs[N[diffs /. numphis /. z -> zval]]]
```
**Output:**
```
Out= Max[0., Abs[-6.16101808370407 - 0.08333333333333333*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*(L[1., 1.][3.4] . L[2., 2.][3.4] + L[1., 1.][3.4] . L[3., 3.][3.4]) - 0.08333333333333333*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*(L[1., 1.][3.4] . L[2., 2.][3.4] + L[1., 1.][3.4] . L[4., 4.][3.4]) - 0.08333333333333333*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*(L[1., 1.][3.4] . L[3., 3.][3.4] + L[1., 1.][3.4] . L[4., 4.][3.4]) - 0.25*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*(L[1., 1.][3.4] . L[2., 2.][3.4] + L[1., 1.][3.4] . L[3., 3.][3.4] + L[1., 1.][3.4] . L[4., 4.][3.4]) - 0.08333333333333333*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*(L[2., 2.][3.4] . L[1., 1.][3.4] + L[2., 2.][3.4] . L[3., 3.][3.4]) - 0.08333333333333333*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*(L[2., 2.][3.4] . L[1., 1.][3.4] + L[2., 2.][3.4] . L[4., 4.][3.4]) - 0.08333333333333333*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*(L[2., 2.][3.4] . L[3., 3.][3.4] + L[2., 2.][3.4] . L[4., 4.][3.4]) - 0.25*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*(L[2., 2.][3.4] . L[1., 1.][3.4] + L[2., 2.][3.4] . L[3., 3.][3.4] + L[2., 2.][3.4] . L[4., 4.][3.4]) - 0.08333333333333333*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*(L[3., 3.][3.4] . L[1., 1.][3.4] + L[3., 3.][3.4] . L[2., 2.][3.4]) - 0.08333333333333333*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*(L[3., 3.][3.4] . L[1., 1.][3.4] + L[3., 3.][3.4] . L[4., 4.][3.4]) - 0.08333333333333333*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*(L[3., 3.][3.4] . L[2., 2.][3.4] + L[3., 3.][3.4] . L[4., 4.][3.4]) - 0.25*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*(L[3., 3.][3.4] . L[1., 1.][3.4] + L[3., 3.][3.4] . L[2., 2.][3.4] + L[3., 3.][3.4] . L[4., 4.][3.4]) - 0.08333333333333333*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*(L[4., 4.][3.4] . L[1., 1.][3.4] + L[4., 4.][3.4] . L[2., 2.][3.4]) - 0.08333333333333333*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*(L[4., 4.][3.4] . L[1., 1.][3.4] + L[4., 4.][3.4] . L[3., 3.][3.4]) - 0.08333333333333333*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*(L[4., 4.][3.4] . L[2., 2.][3.4] + L[4., 4.][3.4] . L[3., 3.][3.4]) - 0.25*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*(L[4., 4.][3.4] . L[1., 1.][
[output truncated — 178717 chars total]
```

## 2026-07-24 08:25:39
**Input:**
```wolfram
c2 = Coefficient[detOp, f''[z]]; terms2 = List @@ c2; Take[terms2, 5]
```
**Output:**
```
Out= {L[1, 1][z] . L[2, 2][z]/12, L[1, 1][z] . L[3, 3][z]/12, (L[1, 1][z] . L[2, 2][z] + L[1, 1][z] . L[3, 3][z])/12, L[1, 1][z] . L[4, 4][z]/12, (L[1, 1][z] . L[2, 2][z] + L[1, 1][z] . L[4, 4][z])/12}
```

## 2026-07-24 08:25:51
**Input:**
```wolfram
c2exp = Expand[c2]; termsE = List @@ c2exp; Length[termsE]; Take[termsE,5]
```
**Output:**
```
Out= {L[1, 1][z] . L[2, 2][z]/2, L[1, 1][z] . L[3, 3][z]/2, L[1, 1][z] . L[4, 4][z]/2, -1/2*L[1, 2][z] . L[2, 1][z], -1/2*L[1, 3][z] . L[3, 1][z]}
```

## 2026-07-24 08:26:01
**Input:**
```wolfram
ClearAll[splitDotTerm, HopVec];
splitDotTerm[term_] := Module[{d, s}, Which[Head[term] === Dot, {1, List @@ term}, Head[term] === Times, d = Cases[term, _Dot, {1}]; If[Length[d] != 1, Print["PARSE FAIL: ", term]; {term, {1}}, s = term/First[d]; {s, List @@ First[d]}], True, {1, {term}}]];
HopVec[k_][rep1_, rep2_][z_][v_] := Module[{coeffExpr, terms, subst}, coeffExpr = Expand[Coefficient[detOp, D[f[z], {z, 4 - k}]]]; terms = If[Head[coeffExpr] === Plus, List @@ coeffExpr, {coeffExpr}]; subst = {L[a_, b_][zz_] :> LL[rep1, rep2][a, b][zz], Derivative[n_][L[a_, b_]][zz_] :> D[LL[rep1, rep2][a, b][zz], {zz, n}]}; Total[(With[{parts = splitDotTerm[#]}, parts[[1]]*Fold[(#2 . #1) &, v, Reverse[parts[[2]] /. subst]]]) & /@ terms]];
```
**Output:**
```
Out= Null
```

## 2026-07-24 08:26:07
**Input:**
```wolfram
zval = 17/5;
diffs = Table[
   (Hop[k][{1,1,0,0},{1,1,0,0}][z] . testv - HopVec[k][{1,1,0,0},{1,1,0,0}][z][testv]),
   {k, 1, 4}];
Max[Abs[N[diffs /. numphis /. z -> zval]]]
```
**Output:**
```
Print:
*STR*PARSE FAIL: \034BoxData[RowBox[{FractionBox["3", "2"], " ", RowBox[{SuperscriptBox[RowBox[{"L", "(", RowBox[{"1", ",", "1"}], ")"}], "\342\200\262", MultilineFunction -> None], "\357\235\255", RowBox[{"(", "z", ")"}]}]}]]*STR*PARSE FAIL: \034BoxData[RowBox[{FractionBox["3", "2"], " ", RowBox[{SuperscriptBox[RowBox[{"L", "(", RowBox[{"2", ",", "2"}], ")"}], "\342\200\262", MultilineFunction -> None], "\357\235\255", RowBox[{"(", "z", ")"}]}]}]]*STR*PARSE FAIL: \034BoxData[RowBox[{FractionBox["3", "2"], " ", RowBox[{SuperscriptBox[RowBox[{"L", "(", RowBox[{"3", ",", "3"}], ")"}], "\342\200\262", MultilineFunction -> None], "\357\235\255", RowBox[{"(", "z", ")"}]}]}]]*STR*PARSE FAIL: \034BoxData[RowBox[{FractionBox["3", "2"], " ", RowBox[{SuperscriptBox[RowBox[{"L", "(", RowBox[{"4", ",", "4"}], ")"}], "\342\200\262", MultilineFunction -> None], "\357\235\255", RowBox[{"(", "z", ")"}]}]}]]*STR*PARSE FAIL: \034BoxData[RowBox[{FractionBox["1", "4"], " ", RowBox[{SuperscriptBox[RowBox[{"L", "(", RowBox[{"1", ",", "1"}], ")"}], TagBox[RowBox[{"(", "3", ")"}], Derivative], MultilineFunction -> None], "\357\235\255", RowBox[{"(", "z", ")"}]}]}]]*STR*PARSE FAIL: \034BoxData[RowBox[{FractionBox["1", "4"], " ", RowBox[{SuperscriptBox[RowBox[{"L", "(", RowBox[{"2", ",", "2"}], ")"}], TagBox[RowBox[{"(", "3", ")"}], Derivative], MultilineFunction -> None], "\357\235\255", RowBox[{"(", "z", ")"}]}]}]]*STR*PARSE FAIL: \034BoxData[RowBox[{FractionBox["1", "4"], " ", RowBox[{SuperscriptBox[RowBox[{"L", "(", RowBox[{"3", ",", "3"}], ")"}], TagBox[RowBox[{"(", "3", ")"}], Derivative], MultilineFunction -> None], "\357\235\255", RowBox[{"(", "z", ")"}]}]}]]*STR*PARSE FAIL: \034BoxData[RowBox[{FractionBox["1", "4"], " ", RowBox[{SuperscriptBox[RowBox[{"L", "(", RowBox[{"4", ",", "4"}], ")"}], TagBox[RowBox[{"(", "3", ")"}], Derivative], MultilineFunction -> None], "\357\235\255", RowBox[{"(", "z", ")"}]}]}]]
Out= Max[0., Abs[-2.3410467128027683 - 1.5*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*Derivative[1][L[1., 1.]][3.4] - 1.5*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*Derivative[1][L[2., 2.]][3.4] - 1.5*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*Derivative[1][L[3., 3.]][3.4] - 1.5*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*Derivative[1][L[4., 4.]][3.4]], Abs[-1.5606978085351788 - 1.5*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*Derivative[1][L[1., 1.]][3.4] - 1.5*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*Derivative[1][L[2., 2.]][3.4] - 1.5*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*Derivative[1][L[3., 3.]][3.4] - 1.5*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*Derivative[1][L[4., 4.]][3.4]], Abs[-0.7803489042675894 - 1.5*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*Derivative[1][L[1., 1.]][3.4] - 1.5*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*Derivative[1][L[2., 2.]][3.4] - 1.5*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*Derivative[1][L[3., 3.]][3.4] - 1.5*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*Derivative[1][L[4., 4.]][3.4]], Abs[-1.5*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*Derivative[1][L[1., 1.]][3.4] - 1.5*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*Derivative[1][L[2., 2.]][3.4] - 1.5*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*Derivative[1][L[3., 3.]][3.4] - 1.5*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*Derivative[1][L[4., 4.]][3.4]], Abs[0.7803489042675894 - 1.5*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*Derivative[1][L[1., 1.]][3.4] - 1.5*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*Derivative[1][L[2., 2.]][3.4] - 1.5*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*Derivative[1][L[3., 3.]][3.4] - 1.5*1. . {-2., -1., -3., 2., -1., 3., 3., 3., 1., -3., -3., 1., -3., 1., 0., -2., -1., 1., -1., -2., 3., -3., -2., 1., 2., 3., -3., 2., 3., 3., 1., 2., 1., 2., 1., 1.}*Derivative[1
[output truncated — 11545 chars total]
```

## 2026-07-24 08:26:26
**Input:**
```wolfram
ClearAll[splitDotTerm, HopVec];
splitDotTerm[term_] := Module[{fs, scalarFs, matFs, mexpr},
  Which[
   Head[term] === Dot, {1, List @@ term},
   Head[term] === Times,
    fs = List @@ term;
    scalarFs = Select[fs, FreeQ[#, L] &];
    matFs = Select[fs, !FreeQ[#, L] &];
    If[Length[matFs] != 1, Print["PARSE FAIL: ", term]; {term, {1}},
     mexpr = First[matFs];
     {Times @@ scalarFs, If[Head[mexpr] === Dot, List @@ mexpr, {mexpr}]}],
   True, {1, {term}}]];
HopVec[k_][rep1_, rep2_][z_][v_] := Module[{coeffExpr, terms, subst}, coeffExpr = Expand[Coefficient[detOp, D[f[z], {z, 4 - k}]]]; terms = If[Head[coeffExpr] === Plus, List @@ coeffExpr, {coeffExpr}]; subst = {L[a_, b_][zz_] :> LL[rep1, rep2][a, b][zz], Derivative[n_][L[a_, b_]][zz_] :> D[LL[rep1, rep2][a, b][zz], {zz, n}]}; Total[(With[{parts = splitDotTerm[#]}, parts[[1]]*Fold[(#2 . #1) &, v, Reverse[parts[[2]] /. subst]]]) & /@ terms]];
```
**Output:**
```
Out= Null
```

## 2026-07-24 08:26:33
**Input:**
```wolfram
zval = 17/5;
diffs = Table[
   (Hop[k][{1,1,0,0},{1,1,0,0}][z] . testv - HopVec[k][{1,1,0,0},{1,1,0,0}][z][testv]),
   {k, 1, 4}];
Max[Abs[N[diffs /. numphis /. z -> zval]]]
```
**Output:**
```
Out= 0.
```

## 2026-07-24 08:26:45
**Input:**
```wolfram
zvals = {17/5, -3/7, 11};
Max[Table[
  Max[Abs[N[(Hop[k][{1,1,0,0},{1,1,0,0}][z] . testv - HopVec[k][{1,1,0,0},{1,1,0,0}][z][testv]) /. numphis /. z -> zv]]],
  {k,1,4},{zv,zvals}]]
```
**Output:**
```
Out= 0.
```

## 2026-07-24 08:26:52
**Input:**
```wolfram
Nbig = dim[{3,3,0,0}]*dim[{5,5,0,0}];
weightsBig = Transpose[Table[Diagonal[EMAll[i,i,{3,3,0,0},{5,5,0,0}]], {i,1,4}]];
idxExtreme = First[Flatten[Position[weightsBig, {8,8,0,0}]]];
idxOther = First[Flatten[Position[weightsBig, {8,0,8,0}]]];
{Nbig, idxExtreme, idxOther}
```
**Output:**
```
[message] First::nofirst: {} has zero length and no first element.
[message] [First::nofirst ×1 more suppressed — first instance shown above]
Out= {9800, First[{}], First[{}]}
```

## 2026-07-24 08:26:58
**Input:**
```wolfram
Take[weightsBig, 3]
Union[Total /@ weightsBig]
```
**Output:**
```
[message] Thread::tdlen: Objects of unequal length in {{8, 8, 0, 0}, {8, 7, 1, 0}, {8, 7, 0, 1}} {16} cannot be combined.
[message] [Thread::tdlen ×2 more suppressed — first instance shown above]
[message] General::stop: Further output of Thread::tdlen will be suppressed during this calculation.
Out= {16}*{{8, 8, 0, 0}, {8, 7, 1, 0}, {8, 7, 0, 1}}
```

## 2026-07-24 08:27:10
**Input:**
```wolfram
idxExtreme = First[Flatten[Position[weightsBig, {8,8,0,0}]]];
idxOther = First[Flatten[Position[weightsBig, {8,0,8,0}]]];
{idxExtreme, idxOther}
```
**Output:**
```
[message] First::nofirst: {} has zero length and no first element.
Out[1]= (suppressed)
[message] First::nofirst: {} has zero length and no first element.
Out[2]= (suppressed)
Out[3]= {First[{}], First[{}]}
```

## 2026-07-24 08:27:18
**Input:**
```wolfram
e1 = weightsBig[[1]]; {e1, e1 === {8,8,0,0}, Head/@e1, ArrayDepth[weightsBig]}
```
**Output:**
```
Out= {SparseArray[Automatic, {4}, 0, {1, {{0, 2}, {{1}, {2}}}, {8, 8}}], False, SparseArray[Automatic, {4}, Integer, {1, {{0, 2}, {{1}, {2}}}, {Integer, Integer}}], 2}
```

## 2026-07-24 08:27:28
**Input:**
```wolfram
weightsBig = Transpose[Table[Normal[Diagonal[EMAll[i,i,{3,3,0,0},{5,5,0,0}]]], {i,1,4}]];
idxExtreme = First[Flatten[Position[weightsBig, {8,8,0,0}]]];
idxOther = First[Flatten[Position[weightsBig, {8,0,8,0}]]];
{idxExtreme, idxOther, Length[Position[weightsBig,{8,8,0,0}]], Length[Position[weightsBig,{8,0,8,0}]]}
```
**Output:**
```
Out= {1, 1192, 1, 1}
```

## 2026-07-24 08:29:32
**Input:**
```wolfram
Timing[colTest = HopVec[4][{3,3,0,0},{5,5,0,0}][z][SparseArray[idxExtreme -> 1, Nbig]];]
```
**Output:**
```
Timed out after 119s — kernel is still alive.
Simplify the expression, wrap it in TimeConstrained manually, or increase timeoutSeconds.
```

## 2026-07-24 08:30:03
**Input:**
```wolfram
{Head[dim[{1,1,0,0}]], ValueQ[detOp], Names["Global`Psi*"], Names["Global`Hop"]}
```
**Output:**
```
Out= {Integer, True, {}, {"Hop"}}
```

## 2026-07-24 08:30:16
**Input:**
```wolfram
dim[{1,1,0,0}]
Dimensions[Hop[1][{1,1,0,0},{1,1,0,0}][z] /. numphis]
ComputeEnergy[{1,1,0,0},{1,1,0,0},{2,2,0,0}]
```
**Output:**
```
Out[1]= 6
Out[2]= {36, 36}
Out[3]= {176/105 + 2/(-1 + z) + 2/z, 86/105 - 4/(-1 + z)^2 + 422/(105*(-1 + z)) - 4/z^2 + 2/(105*z), 16/105 + 4/(-1 + z)^3 - 772/(105*(-1 + z)^2) + 6/(5*(-1 + z)) + 4/z^3 + 68/(105*z^2) - 6/(35*z), 1/105 + 20/(3*(-1 + z)^3) - 8/(7*(-1 + z)^2) + 2/(21*(-1 + z)) - 4/(3*z^3) + 8/(35*z^2) - 2/(105*z)}
```

## 2026-07-24 08:30:47
**Input:**
```wolfram
r1={1,1,0,0};r2={1,1,0,0};ra={2,2,0,0};
nn=dim[r1] dim[r2];
psiV=Array[ψ,nn];
En=ComputeEnergy[r1,r2,ra];
Hk=Table[Hop[k][r1,r2][z]/.numphis,{k,4}];
res=Flatten[Table[Hk[[k]].psiV-En[[k]] psiV,{k,4}]];
lin=Flatten[CoefficientList[Numerator[Together[#]],z]&/@res];
mat=Normal[CoefficientArrays[lin,psiV][[2]]];
ns=NullSpace[mat];
{Length[lin],Dimensions[mat],Length[ns]}
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= (suppressed)
Out[5]= (suppressed)
Out[6]= (suppressed)
Out[7]= (suppressed)
Out[8]= (suppressed)
Out[9]= (suppressed)
Out[10]= {576, {576, 36}, 0}
```

## 2026-07-24 08:31:27
**Input:**
```wolfram
wts=Transpose[Table[Diagonal[Normal[EMAll[i,i,r1,r2]]],{i,4}]];
eidx=Position[wts,{2,2,0,0}];
{eidx,Length[eidx]}
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= {{{1}}, 1}
```

## 2026-07-24 08:31:34
**Input:**
```wolfram
ee=UnitVector[nn,1];
he=Table[(Hk[[k]].ee),{k,4}];
(*is the image proportional to ee? check entry 1 and whether off-entries vanish*)
proj=Table[he[[k,1]],{k,4}];
offnorm=Table[Simplify[Total[Abs[Delete[he[[k]],1]]]],{k,4}];
{Simplify[proj-En],offnorm}
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= (suppressed)
Out[5]= {{-176/105 + ϕ[1] + ϕ[2] + ϕ[3] + ϕ[4], -1/105*(-212 + 105*ϕ[1] + 105*ϕ[2] + 210*ϕ[3] + 210*ϕ[4] + z*(338 - 420*ϕ[3] - 420*ϕ[4] + 105*ϕ[3]*ϕ[4] + 105*ϕ[2]*(-2 + ϕ[3] + ϕ[4]) + 105*ϕ[1]*(-2 + ϕ[2] + ϕ[3] + ϕ[4])) - z^2*(-86 + 105*ϕ[3]*ϕ[4] + 105*ϕ[2]*(ϕ[3] + ϕ[4]) + 105*ϕ[1]*(ϕ[2] + ϕ[3] + ϕ[4])))/((-1 + z)*z), (-2*(-176 + 105*ϕ[1] + 105*ϕ[2] + 105*ϕ[3] + 105*ϕ[4]) + 7*z*(-98 + 30*ϕ[3] + 30*ϕ[4] + 30*ϕ[3]*ϕ[4] + 15*ϕ[1]*(4 + ϕ[3] + ϕ[4]) + 15*ϕ[2]*(4 + ϕ[3] + ϕ[4])) + z^2*(778 - 210*ϕ[3] + 105*ϕ[2]*(-4 + ϕ[3]*(-3 + ϕ[4]) - 3*ϕ[4]) - 210*ϕ[4] - 630*ϕ[3]*ϕ[4] + 105*ϕ[1]*(-4 + (-3 + ϕ[2])*ϕ[4] + ϕ[3]*(-3 + ϕ[2] + ϕ[4]))) + z^3*(-76 - 210*ϕ[2]*(ϕ[3]*(-1 + ϕ[4]) - ϕ[4]) + 420*ϕ[3]*ϕ[4] - 210*ϕ[1]*((-1 + ϕ[2])*ϕ[4] + ϕ[3]*(-1 + ϕ[2] + ϕ[4]))) + z^4*(-16 + 105*ϕ[2]*ϕ[3]*ϕ[4] + 105*ϕ[1]*(ϕ[3]*ϕ[4] + ϕ[2]*(ϕ[3] + ϕ[4]))))/(105*(-1 + z)^2*z^2), (-70*(-4 + 3*ϕ[1] + 3*ϕ[2]) + z^6*(-1 + 105*ϕ[1]*ϕ[2]*ϕ[3]*ϕ[4]) - 5*z^5*(1 - 42*ϕ[2]*ϕ[3]*ϕ[4] + 21*ϕ[1]*(-2 + 3*ϕ[2])*ϕ[3]*ϕ[4]) + 3*z*(35*ϕ[1]*(6 + ϕ[3] + ϕ[4]) + 35*ϕ[2]*(6 + ϕ[3] + ϕ[4]) - 2*(136 + 35*ϕ[3] + 35*ϕ[4])) - z^2*(105*ϕ[1]*(3*(2 + ϕ[4]) + ϕ[3]*(3 + ϕ[4])) + 105*ϕ[2]*(3*(2 + ϕ[4]) + ϕ[3]*(3 + ϕ[4])) - 2*(383 + 315*ϕ[4] + 105*ϕ[3]*(3 + ϕ[4]))) + z^3*(-611 - 420*ϕ[3] - 420*ϕ[4] - 420*ϕ[3]*ϕ[4] + 420*ϕ[2]*(1 + ϕ[3])*(1 + ϕ[4]) - 105*ϕ[1]*(-4*(1 + ϕ[4]) + ϕ[3]*(-4 + (-4 + ϕ[2])*ϕ[4]))) + z^4*(107 + 210*ϕ[3]*ϕ[4] - 105*ϕ[2]*(2*ϕ[4] + ϕ[3]*(2 + 5*ϕ[4])) + 105*ϕ[1]*(-2*ϕ[4] + ϕ[3]*(-2 + (-5 + 3*ϕ[2])*ϕ[4]))))/(105*(-1 + z)^3*z^3)}, {0, 0, 0, 0}}
```

## 2026-07-24 08:32:30
**Input:**
```wolfram
FreeQ[Hk[[1]],ϕ]
Simplify[(proj-En)/.numphis]
FreeQ[Hop[1][r1,r2][z],ϕ]
```
**Output:**
```
Out[1]= True
Out[2]= {0, 0, 0, 0}
Out[3]= True
```

## 2026-07-24 08:33:20
**Input:**
```wolfram
{DownValues[ϕ], OwnValues[ϕ]}
{ϕ[1],ϕ[2],ϕ[3],ϕ[4]}
FreeQ[En,ϕ]
```
**Output:**
```
Out[1]= {{}, {}}
Out[2]= {ϕ[1], ϕ[2], ϕ[3], ϕ[4]}
Out[3]= True
```

## 2026-07-24 08:34:21
**Input:**
```wolfram
hopSyms=DeleteDuplicates[Cases[Hop[1][r1,r2][z]/.numphis, h_Symbol[i_Integer]:>h, Infinity]];
nh=numphis[[1,1,0]];
{"Hop-heads"->({SymbolName[#],ToCharacterCode[SymbolName[#]]}&/@hopSyms),"numphis-head"->{SymbolName[nh],ToCharacterCode[SymbolName[nh]]}}
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= {"Hop-heads" -> {}, "numphis-head" -> {"ϕ", {981}}}
```

## 2026-07-24 08:34:52
**Input:**
```wolfram
r1={1,1,0,0};r2={1,1,0,0};ra={2,2,0,0};
nn=dim[r1] dim[r2];
psiV=Array[ψ,nn];
Hk=Table[Together[Hop[k][r1,r2][z]/.numphis],{k,4}];
En=Together/@ComputeEnergy[r1,r2,ra];
{FreeQ[Hk,nh],FreeQ[En,nh]}
res=Flatten[Table[Hk[[k]].psiV-En[[k]] psiV,{k,4}]];
lin=Flatten[CoefficientList[Numerator[Together[#]],z]&/@res];
mat=Normal[CoefficientArrays[lin,psiV][[2]]];
ns=NullSpace[mat];
{Dimensions[mat],Length[ns],If[ns=!={},{Position[ns[[1]],x_/;x=!=0],ns[[1]]//Normal//Short},None]}
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= (suppressed)
Out[5]= (suppressed)
Out[6]= {True, True}
Out[7]= (suppressed)
Out[8]= (suppressed)
Out[9]= (suppressed)
Out[10]= (suppressed)
Out[11]= {{576, 36}, 0, None}
```

## 2026-07-24 08:35:25
**Input:**
```wolfram
res0=res/.Table[ψ[i]->KroneckerDelta[i,1],{i,nn}];
{"res at e1 all zero?"->(Simplify[res0]//Union),
 "mat.e1 zero?"->(Union[mat . UnitVector[nn,1]]),
 "mat rank"->MatrixRank[mat],
 "any pure-const lin entries?"->(Union[Cases[lin,_?(FreeQ[#,ψ]&)]])}
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= {"res at e1 all zero?" -> {0, -176/105 + ϕ[1] + ϕ[2] + ϕ[3] + ϕ[4], -1/105*(-212 + 105*ϕ[1] + 105*ϕ[2] + 210*ϕ[3] + 210*ϕ[4] + z*(338 - 420*ϕ[3] - 420*ϕ[4] + 105*ϕ[3]*ϕ[4] + 105*ϕ[2]*(-2 + ϕ[3] + ϕ[4]) + 105*ϕ[1]*(-2 + ϕ[2] + ϕ[3] + ϕ[4])) - z^2*(-86 + 105*ϕ[3]*ϕ[4] + 105*ϕ[2]*(ϕ[3] + ϕ[4]) + 105*ϕ[1]*(ϕ[2] + ϕ[3] + ϕ[4])))/((-1 + z)*z), (-2*(-176 + 105*ϕ[1] + 105*ϕ[2] + 105*ϕ[3] + 105*ϕ[4]) + 7*z*(-98 + 30*ϕ[3] + 30*ϕ[4] + 30*ϕ[3]*ϕ[4] + 15*ϕ[1]*(4 + ϕ[3] + ϕ[4]) + 15*ϕ[2]*(4 + ϕ[3] + ϕ[4])) + z^2*(778 - 210*ϕ[3] + 105*ϕ[2]*(-4 + ϕ[3]*(-3 + ϕ[4]) - 3*ϕ[4]) - 210*ϕ[4] - 630*ϕ[3]*ϕ[4] + 105*ϕ[1]*(-4 + (-3 + ϕ[2])*ϕ[4] + ϕ[3]*(-3 + ϕ[2] + ϕ[4]))) + z^3*(-76 - 210*ϕ[2]*(ϕ[3]*(-1 + ϕ[4]) - ϕ[4]) + 420*ϕ[3]*ϕ[4] - 210*ϕ[1]*((-1 + ϕ[2])*ϕ[4] + ϕ[3]*(-1 + ϕ[2] + ϕ[4]))) + z^4*(-16 + 105*ϕ[2]*ϕ[3]*ϕ[4] + 105*ϕ[1]*(ϕ[3]*ϕ[4] + ϕ[2]*(ϕ[3] + ϕ[4]))))/(105*(-1 + z)^2*z^2), (-70*(-4 + 3*ϕ[1] + 3*ϕ[2]) + z^6*(-1 + 105*ϕ[1]*ϕ[2]*ϕ[3]*ϕ[4]) - 5*z^5*(1 - 42*ϕ[2]*ϕ[3]*ϕ[4] + 21*ϕ[1]*(-2 + 3*ϕ[2])*ϕ[3]*ϕ[4]) + 3*z*(35*ϕ[1]*(6 + ϕ[3] + ϕ[4]) + 35*ϕ[2]*(6 + ϕ[3] + ϕ[4]) - 2*(136 + 35*ϕ[3] + 35*ϕ[4])) - z^2*(105*ϕ[1]*(3*(2 + ϕ[4]) + ϕ[3]*(3 + ϕ[4])) + 105*ϕ[2]*(3*(2 + ϕ[4]) + ϕ[3]*(3 + ϕ[4])) - 2*(383 + 315*ϕ[4] + 105*ϕ[3]*(3 + ϕ[4]))) + z^3*(-611 - 420*ϕ[3] - 420*ϕ[4] - 420*ϕ[3]*ϕ[4] + 420*ϕ[2]*(1 + ϕ[3])*(1 + ϕ[4]) - 105*ϕ[1]*(-4*(1 + ϕ[4]) + ϕ[3]*(-4 + (-4 + ϕ[2])*ϕ[4]))) + z^4*(107 + 210*ϕ[3]*ϕ[4] - 105*ϕ[2]*(2*ϕ[4] + ϕ[3]*(2 + 5*ϕ[4])) + 105*ϕ[1]*(-2*ϕ[4] + ϕ[3]*(-2 + (-5 + 3*ϕ[2])*ϕ[4]))))/(105*(-1 + z)^3*z^3)}, "mat.e1 zero?" -> {0, 280 - 210*ϕ[1] - 210*ϕ[2], 352 - 210*ϕ[1] - 210*ϕ[2] - 210*ϕ[3] - 210*ϕ[4], 212 - 105*ϕ[1] - 105*ϕ[2] - 210*ϕ[3] - 210*ϕ[4], -176 + 105*ϕ[1] + 105*ϕ[2] + 105*ϕ[3] + 105*ϕ[4], -816 + 630*ϕ[1] + 630*ϕ[2] - 210*ϕ[3] + 105*ϕ[1]*ϕ[3] + 105*ϕ[2]*ϕ[3] - 210*ϕ[4] + 105*ϕ[1]*ϕ[4] + 105*ϕ[2]*ϕ[4], -338 + 210*ϕ[1] + 210*ϕ[2] - 105*ϕ[1]*ϕ[2] + 420*ϕ[3] - 105*ϕ[1]*ϕ[3] - 105*ϕ[2]*ϕ[3] + 420*ϕ[4] - 105*ϕ[1]*ϕ[4] - 105*ϕ[2]*ϕ[4] - 105*ϕ[3]*ϕ[4], -86 + 105*ϕ[1]*ϕ[2] + 105*ϕ[1]*ϕ[3] + 105*ϕ[2]*ϕ[3] + 105*ϕ[1]*ϕ[4] + 105*ϕ[2]*ϕ[4] + 105*ϕ[3]*ϕ[4], -686 + 420*ϕ[1] + 420*ϕ[2] + 210*ϕ[3] + 105*ϕ[1]*ϕ[3] + 105*ϕ[2]*ϕ[3] + 210*ϕ[4] + 105*ϕ[1]*ϕ[4] + 105*ϕ[2]*ϕ[4] + 210*ϕ[3]*ϕ[4], -76 + 210*ϕ[1]*ϕ[3] + 210*ϕ[2]*ϕ[3] - 210*ϕ[1]*ϕ[2]*ϕ[3] + 210*ϕ[1]*ϕ[4] + 210*ϕ[2]*ϕ[4] - 210*ϕ[1]*ϕ[2]*ϕ[4] + 420*ϕ[3]*ϕ[4] - 210*ϕ[1]*ϕ[3]*ϕ[4] - 210*ϕ[2]*ϕ[3]*ϕ[4], 766 - 630*ϕ[1] - 630*ϕ[2] + 630*ϕ[3] - 315*ϕ[1]*ϕ[3] - 315*ϕ[2]*ϕ[3] + 630*ϕ[4] - 315*ϕ[1]*ϕ[4] - 315*ϕ[2]*ϕ[4] + 210*ϕ[3]*ϕ[4] - 105*ϕ[1]*ϕ[3]*ϕ[4] - 105*ϕ[2]*ϕ[3]*ϕ[4], -16 + 105*ϕ[1]*ϕ[2]*ϕ[3] + 105*ϕ[1]*ϕ[2]*ϕ[4] + 105*ϕ[1]*ϕ[3]*ϕ[4] + 105*ϕ[2]*ϕ[3]*ϕ[4], 778 - 420*ϕ[1] - 420*ϕ[2] - 210*ϕ[3] - 315*ϕ[1]*ϕ[3] - 315*ϕ[2]*ϕ[3] + 105*ϕ[1]*ϕ[2]*ϕ[3] - 210*ϕ[4] - 315*ϕ[1]*ϕ[4] - 315*ϕ[2]*ϕ[4] + 105*ϕ[1]*ϕ[2]*ϕ[4] - 630*ϕ[3]*ϕ[4] + 105*ϕ[1]*ϕ[3]*ϕ[4] + 105*ϕ[2]*ϕ[3]*ϕ[4], -5 + 210*ϕ[1]*ϕ[3]*ϕ[4] + 210*ϕ[2]*ϕ[3]*ϕ[4] - 315*ϕ[1]*ϕ[2]*ϕ[3]*ϕ[4], -611 + 420*ϕ[1] + 420*ϕ[2] - 420*ϕ[3] + 420*ϕ[1]*ϕ[3] + 420*ϕ[2]*ϕ[3] - 420*ϕ[4] + 420*ϕ[1]*ϕ[4] + 420*ϕ[2]*ϕ[4] - 420*ϕ[3]*ϕ[4] + 420*ϕ[1]*ϕ[3]*ϕ[4] + 420*ϕ[2]*ϕ[3]*ϕ[4] - 105*ϕ[1]*ϕ[2]*ϕ[3]*ϕ[4], -1 + 105*ϕ[1]*ϕ[2]*ϕ[3]*ϕ[4], 107 - 210*ϕ[1]*ϕ[3] - 210*ϕ[2]*ϕ[3] - 210*ϕ[1]*ϕ[4] - 210*ϕ[2]*ϕ[4] + 210*ϕ[3]*ϕ[4] - 525*ϕ[1]*ϕ[3]*ϕ[4] - 525*ϕ[2]*ϕ[3]*ϕ[4] + 315*ϕ[1]*ϕ[2]*ϕ[3]*ϕ[4]}, "mat rank" -> 36, "any pure-const lin entries?" -> {}}
```

## 2026-07-24 08:36:10
**Input:**
```wolfram
resSyms=DeleteDuplicates[Cases[res0,s_Symbol[_Integer]:>s,Infinity]];
{"res-ϕ"->({SymbolName[#],Context[#],ToCharacterCode[SymbolName[#]]}&/@resSyms),
 "numphis-ϕ"->{SymbolName[nh],Context[nh],ToCharacterCode[SymbolName[nh]]},
 "same?"->(MemberQ[resSyms,nh])}
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= {"res-ϕ" -> {{"ϕ", "Global`", {981}}}, "numphis-ϕ" -> {"ϕ", "Global`", {981}}, "same?" -> True}
```

## 2026-07-24 08:36:49
**Input:**
```wolfram
h1=Hop[1][r1,r2][z];
{"Head"->Head[h1],
 "FreeQ raw sparse /.numphis"->FreeQ[h1/.numphis,nh],
 "FreeQ Normal-then-subst"->FreeQ[Normal[h1]/.numphis,nh],
 "FreeQ subst-then-Normal"->FreeQ[Normal[h1/.numphis],nh]}
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= {"Head" -> SparseArray, "FreeQ raw sparse /.numphis" -> True, "FreeQ Normal-then-subst" -> True, "FreeQ subst-then-Normal" -> False}
```

## 2026-07-24 08:37:14
**Input:**
```wolfram
Hk=Table[Normal[Hop[k][r1,r2][z]]/.numphis,{k,4}];
En=ComputeEnergy[r1,r2,ra];
{"Hk clean"->FreeQ[Hk,nh],"En clean"->FreeQ[En,nh]}
res=Flatten[Table[Hk[[k]].psiV-En[[k]] psiV,{k,4}]];
lin=Flatten[CoefficientList[Numerator[Together[#]],z]&/@res];
mat=Normal[CoefficientArrays[lin,psiV][[2]]];
ns=NullSpace[mat];
{"nullspace dim"->Length[ns],"support"->If[ns=!={},Flatten[Position[ns[[1]],x_/;x=!=0]],None],"vec"->If[ns=!={},ns[[1]],None]}
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= {"Hk clean" -> True, "En clean" -> True}
Out[4]= (suppressed)
Out[5]= (suppressed)
Out[6]= (suppressed)
Out[7]= (suppressed)
Out[8]= {"nullspace dim" -> 1, "support" -> {0, 1}, "vec" -> {-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}
```

## 2026-07-24 08:38:06
**Input:**
```wolfram
ClearAll[PsiSolve];
PsiSolve[rep1_,rep2_,repall_]:=Module[{nn,psiV,Hk,En,res,lin,mat},
  nn=dim[rep1] dim[rep2];
  psiV=Array[ψ,nn];
  Hk=Table[Normal[Hop[k][rep1,rep2][z]]/.numphis,{k,4}];
  En=ComputeEnergy[rep1,rep2,repall];
  res=Flatten[Table[Hk[[k]].psiV-En[[k]] psiV,{k,4}]];
  lin=Flatten[CoefficientList[Numerator[Together[#]],z]&/@res];
  mat=Normal[CoefficientArrays[lin,psiV][[2]]];
  NullSpace[mat]];
(* validate over ALL good weights of the small case: each should give a 1-dim nullspace,
   a single basis vector, whose EMAll-weight equals repall *)
gw=FindGoodWeight2[{1,1,0,0},{1,1,0,0}];
wtmat=Transpose[Table[Diagonal[Normal[EMAll[i,i,{1,1,0,0},{1,1,0,0}]]],{i,4}]];
check=Table[Module[{ns=PsiSolve[{1,1,0,0},{1,1,0,0},w],supp,wt},
   supp=If[ns=!={},Flatten[Position[ns[[1]],x_/;x=!=0]],{}];
   wt=If[Length[supp]==1,wtmat[[supp[[1]]]],"multi"];
   <|"repall"->w,"nsdim"->Length[ns],"support"->supp,"stateWeight"->wt,"weightOK"->(wt===w)|>],{w,gw}];
Column[check]
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= (suppressed)
Out[5]= (suppressed)
Out[6]= Column[{<|"repall" -> {2, 2, 0, 0}, "nsdim" -> 1, "support" -> {0, 1}, "stateWeight" -> "multi", "weightOK" -> False|>, <|"repall" -> {2, 0, 2, 0}, "nsdim" -> 1, "support" -> {0, 8}, "stateWeight" -> "multi", "weightOK" -> False|>, <|"repall" -> {2, 0, 0, 2}, "nsdim" -> 1, "support" -> {0, 15}, "stateWeight" -> "multi", "weightOK" -> False|>, <|"repall" -> {0, 2, 2, 0}, "nsdim" -> 1, "support" -> {0, 22}, "stateWeight" -> "multi", "weightOK" -> False|>, <|"repall" -> {0, 2, 0, 2}, "nsdim" -> 1, "support" -> {0, 29}, "stateWeight" -> "multi", "weightOK" -> False|>, <|"repall" -> {0, 0, 2, 2}, "nsdim" -> 1, "support" -> {0, 36}, "stateWeight" -> "multi", "weightOK" -> False|>}]
```

## 2026-07-24 08:38:58
**Input:**
```wolfram
check2=Table[Module[{ns=PsiSolve[{1,1,0,0},{1,1,0,0},w],pos,wt},
   pos=Flatten[Position[ns[[1]],_?(#=!=0&),{1}]];
   wt=If[Length[pos]==1,wtmat[[pos[[1]]]],"multi:"<>ToString[pos]];
   <|"repall"->w,"nsdim"->Length[ns],"stateIndex"->pos,"stateWeight"->wt,"weightOK"->(wt===w)|>],{w,gw}];
{Column[check2],"ALL PASS"->AllTrue[check2,#["nsdim"]==1&&#["weightOK"]&]}
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= {Column[{<|"repall" -> {2, 2, 0, 0}, "nsdim" -> 1, "stateIndex" -> {0, 1}, "stateWeight" -> "multi:{0, 1}", "weightOK" -> False|>, <|"repall" -> {2, 0, 2, 0}, "nsdim" -> 1, "stateIndex" -> {0, 8}, "stateWeight" -> "multi:{0, 8}", "weightOK" -> False|>, <|"repall" -> {2, 0, 0, 2}, "nsdim" -> 1, "stateIndex" -> {0, 15}, "stateWeight" -> "multi:{0, 15}", "weightOK" -> False|>, <|"repall" -> {0, 2, 2, 0}, "nsdim" -> 1, "stateIndex" -> {0, 22}, "stateWeight" -> "multi:{0, 22}", "weightOK" -> False|>, <|"repall" -> {0, 2, 0, 2}, "nsdim" -> 1, "stateIndex" -> {0, 29}, "stateWeight" -> "multi:{0, 29}", "weightOK" -> False|>, <|"repall" -> {0, 0, 2, 2}, "nsdim" -> 1, "stateIndex" -> {0, 36}, "stateWeight" -> "multi:{0, 36}", "weightOK" -> False|>}], "ALL PASS" -> False}
```

## 2026-07-24 08:39:35
**Input:**
```wolfram
nzIdx[v_]:=Select[Range[Length[v]],v[[#]]=!=0&];
check3=Table[Module[{ns=PsiSolve[{1,1,0,0},{1,1,0,0},w],pos,wt},
   pos=nzIdx[ns[[1]]];
   wt=If[Length[pos]==1,wtmat[[pos[[1]]]],"multi"];
   <|"repall"->w,"nsdim"->Length[ns],"stateIndex"->pos,"stateWeight"->wt,"weightOK"->(wt===w)|>],{w,gw}];
{Column[check3],"ALL PASS"->AllTrue[check3,#["nsdim"]==1&&#["weightOK"]&]}
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= {Column[{<|"repall" -> {2, 2, 0, 0}, "nsdim" -> 1, "stateIndex" -> {1}, "stateWeight" -> {2, 2, 0, 0}, "weightOK" -> True|>, <|"repall" -> {2, 0, 2, 0}, "nsdim" -> 1, "stateIndex" -> {8}, "stateWeight" -> {2, 0, 2, 0}, "weightOK" -> True|>, <|"repall" -> {2, 0, 0, 2}, "nsdim" -> 1, "stateIndex" -> {15}, "stateWeight" -> {2, 0, 0, 2}, "weightOK" -> True|>, <|"repall" -> {0, 2, 2, 0}, "nsdim" -> 1, "stateIndex" -> {22}, "stateWeight" -> {0, 2, 2, 0}, "weightOK" -> True|>, <|"repall" -> {0, 2, 0, 2}, "nsdim" -> 1, "stateIndex" -> {29}, "stateWeight" -> {0, 2, 0, 2}, "weightOK" -> True|>, <|"repall" -> {0, 0, 2, 2}, "nsdim" -> 1, "stateIndex" -> {36}, "stateWeight" -> {0, 0, 2, 2}, "weightOK" -> True|>}], "ALL PASS" -> True}
```

## 2026-07-24 08:40:03
**Input:**
```wolfram
{"dim{3,3,0,0}"->dim[{3,3,0,0}],"dim{5,5,0,0}"->dim[{5,5,0,0}],"tensor"->dim[{3,3,0,0}]dim[{5,5,0,0}],
 "dim{2,2,0,0}"->dim[{2,2,0,0}],"2x2 tensor"->dim[{2,2,0,0}]^2}
```
**Output:**
```
Out[1]= {"dim{3,3,0,0}" -> 50, "dim{5,5,0,0}" -> 196, "tensor" -> 9800, "dim{2,2,0,0}" -> 20, "2x2 tensor" -> 400}
```

## 2026-07-24 08:41:30
**Input:**
```wolfram
wtmat22=Transpose[Table[Diagonal[Normal[EMAll[i,i,{2,2,0,0},{2,2,0,0}]]],{i,4}]];
{tt,ns22}=AbsoluteTiming[PsiSolve[{2,2,0,0},{2,2,0,0},{4,4,0,0}]];
pos22=nzIdx[ns22[[1]]];
{"time"->tt,"nsdim"->Length[ns22],"stateIndex"->pos22,"weight"->If[Length[pos22]==1,wtmat22[[pos22[[1]]]],"multi"],"weightOK"->(Length[pos22]==1&&wtmat22[[pos22[[1]]]]==={4,4,0,0})}
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= {"time" -> 63.997281, "nsdim" -> 1, "stateIndex" -> {1}, "weight" -> {4, 4, 0, 0}, "weightOK" -> True}
```

## 2026-07-24 08:41:57
**Input:**
```wolfram
TimeConstrained[
  First[AbsoluteTiming[hbig1=Hop[1][{3,3,0,0},{5,5,0,0}][z];]],
  150, "TIMED OUT (>150s)"]
```
**Output:**
```
Out= 0.113062
```

## 2026-07-24 08:46:23
**Input:**
```wolfram
r1b={3,3,0,0};r2b={5,5,0,0};
wtB=Transpose[Table[Diagonal[Normal[EMAll[i,i,r1b,r2b]]],{i,4}]];
HkB=Table[Hop[k][r1b,r2b][z],{k,4}];
nB=dim[r1b] dim[r2b];
directCheck[w_]:=Module[{p,En,resid},
  p=Position[wtB,w][[1,1]];
  En=ComputeEnergy[r1b,r2b,w];
  resid=Table[Together[(Normal[HkB[[k]][[All,p]]]/.numphis)-En[[k]] UnitVector[nB,p]],{k,4}];
  <|"w"->w,"stateIndex"->p,"eigenOK"->(Simplify[Total[Abs[Flatten[resid]]]]===0)|>];
directCheck[{8,8,0,0}]
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]: timed out — kernel interrupted. Retry with a higher timeoutSeconds.
```

## 2026-07-24 08:46:47
**Input:**
```wolfram
pex=Position[wtB,{8,8,0,0}][[1,1]];
col1=HkB[[1]][[All,pex]];
{"alive"->True,"stateIndex"->pex,"nnz col(H1)"->Length[Most[ArrayRules[col1]]],
 "nnz cols H1..H4"->Table[Length[Most[ArrayRules[HkB[[k]][[All,pex]]]]],{k,4}]}
```
**Output:**
```
Out[1]= (suppressed)
[message] Part::partd: Part specification HkB[[1]] is longer than depth of object.
[message] [Part::partd ×1 more suppressed — first instance shown above]
Out[2]= (suppressed)
[message] ArrayRules::rect: Nonrectangular array encountered.
[message] ArrayRules::argt: ArrayRules called with 0 arguments; 1 or 2 arguments are expected.
[message] Part::partd: Part specification HkB[[1]] is longer than depth of object.
[message] [Part::partd ×2 more suppressed — first instance shown above]
[message] [ArrayRules::rect ×2 more suppressed — first instance shown above]
[message] [ArrayRules::argt ×2 more suppressed — first instance shown above]
[message] General::stop: Further output of Part::partd will be suppressed during this calculation.
[message] [General::stop ×2 more suppressed — first instance shown above]
Out[3]= {"alive" -> True, "stateIndex" -> 1, "nnz col(H1)" -> 0, "nnz cols H1..H4" -> {0, 0, 0, 0}}
```

## 2026-07-24 08:48:59
**Input:**
```wolfram
ClearAll[HkB];
Do[HkB[k]=Hop[k][r1b,r2b][z],{k,4}];
{"heads"->Table[Head[HkB[k]],{k,4}],"dims"->Table[Dimensions[HkB[k]],{k,4}]}
```
**Output:**
```
Out[1]= (suppressed)
Out[2]: timed out — kernel interrupted. Retry with a higher timeoutSeconds.
```

## 2026-07-24 08:49:26
**Input:**
```wolfram
{"alive"->(1+1),"PsiSolve def?"->(Head[PsiSolve]===Symbol&&DownValues[PsiSolve]=!={}),"small case recheck"->AllTrue[check3,#["weightOK"]&]}
```
**Output:**
```
Out[1]= {"alive" -> 2, "PsiSolve def?" -> True, "small case recheck" -> True}
```

## 2026-07-24 08:52:57
**Input:**
```wolfram
TimeConstrained[First[AbsoluteTiming[h4=Hop[4][{3,3,0,0},{3,3,0,0}][z];]],120,"TIMEOUT>120s"]
```
**Output:**
```
Out= TIMEOUT>120s
```

## 2026-07-24 08:56:04
**Input:**
```wolfram
dc[1]=Coefficient[detOp,f'''[z]];
dc[2]=Coefficient[detOp,f''[z]];
dc[3]=Coefficient[detOp,f'[z]];
dc[4]=Coefficient[detOp,f[z]];
ClearAll[applyLLnum,HopNum];
applyLLnum[rep1_,rep2_,zval_][expr_]:=expr/.{
  Derivative[n_][L[a_,b_]][zz_]:>(D[LL[rep1,rep2][a,b][ζ],{ζ,n}]/.ζ->zval),
  L[a_,b_][zz_]:>(LL[rep1,rep2][a,b][zz]/.zz->zval)};
HopNum[k_][rep1_,rep2_,zval_]:=applyLLnum[rep1,rep2,zval/.numphis][dc[k]/.numphis];
TimeConstrained[First[AbsoluteTiming[hn4=HopNum[4][{3,3,0,0},{3,3,0,0},2];]],120,"TIMEOUT"]
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= (suppressed)
Out[5]= (suppressed)
Out[6]= (suppressed)
Out[7]= (suppressed)
Out[8]= TIMEOUT
```

## 2026-07-24 08:58:30
**Input:**
```wolfram
ClearAll[LLn,MvecN,HopActEp];
LLn[rep1_,rep2_][i_,j_]:=LLn[rep1,rep2][i,j]=(ϕ[j]/.numphis) KroneckerDelta[i,j] SparseArray[IdentityMatrix[dim[rep1] dim[rep2]]]+KroneckerProduct[EM[i,j,rep1],Id[rep2]]/z+KroneckerProduct[Id[rep1],EM[i,j,rep2]]/(z-1);
MvecN[rep1_,rep2_][{i_,j_}][g_]:=D[g,z] KroneckerDelta[i,j]+LLn[rep1,rep2][i,j].g;
HopActEp[rep1_,rep2_,p_]:=Module[{ep,g0,res,pm=Permutations[Range[4]]},
  ep=SparseArray[{p->1},{dim[rep1] dim[rep2]}];
  g0=f[z] ep;
  res=(1/24)Sum[Signature[ip]Signature[jp] MvecN[rep1,rep2][{ip[[1]],jp[[1]]}][MvecN[rep1,rep2][{ip[[2]],jp[[2]]}][MvecN[rep1,rep2][{ip[[3]],jp[[3]]}][MvecN[rep1,rep2][{ip[[4]],jp[[4]]}][g0]]]],{ip,pm},{jp,pm}];
  Table[Together[Coefficient[res,Derivative[4-k][f][z]]],{k,4}]];
(* validate against known small case: Hop.e1 should equal En_k * e1 *)
En11=ComputeEnergy[{1,1,0,0},{1,1,0,0},{2,2,0,0}];
act11=HopActEp[{1,1,0,0},{1,1,0,0},1];
Table[Simplify[Total[Abs[Normal[act11[[k]]]-En11[[k]] UnitVector[36,1]]]],{k,4}]
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= (suppressed)
Out[5]= (suppressed)
Out[6]= (suppressed)
Out[7]= {0, 0, 0, 0}
```

## 2026-07-24 08:59:36
**Input:**
```wolfram
ClearAll[Psi,weightList];
(* joint spectrum of the global Cartans EMAll[j,j] (already diagonal in the GT tensor basis) *)
weightList[rep1_,rep2_]:=weightList[rep1,rep2]=Transpose[Table[Diagonal[Normal[EMAll[j,j,rep1,rep2]]],{j,4}]];
(* good weight => unique basis state; return it as the (unit) eigenvector *)
Psi[rep1_,rep2_,repall_]:=Module[{pos=Position[weightList[rep1,rep2],repall]},
  If[Length[pos]!=1,Return[$Failed[ "multiplicity="<>ToString[Length[pos]] ]]];
  SparseArray[{pos[[1,1]]->1},{dim[rep1] dim[rep2]}]];
(* validate: small case matches the earlier PsiSolve supports {1,8,15,22,29,36} *)
{"small-case indices"->(FirstPosition[Normal[Psi[{1,1,0,0},{1,1,0,0},#]],1][[1]]&/@FindGoodWeight2[{1,1,0,0},{1,1,0,0}])}
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= (suppressed)
Out[4]= {"small-case indices" -> {1, 8, 15, 22, 29, 36}}
```

## 2026-07-24 08:59:45
**Input:**
```wolfram
gw33=FindGoodWeight2[{3,3,0,0},{3,3,0,0}];
tt33=First[AbsoluteTiming[res33=Table[Module[{p=Psi[{3,3,0,0},{3,3,0,0},w]},
   <|"w"->w,"ok"->(Head[p]===SparseArray),"index"->If[Head[p]===SparseArray,FirstPosition[Normal[p],1][[1]],p]|>],{w,gw33}];]];
{"tensor dim"->dim[{3,3,0,0}]^2,"time(s)"->tt33,"num good weights"->Length[gw33],"all unique?"->AllTrue[res33,#["ok"]&],Column[res33]}
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= (suppressed)
Out[3]= {"tensor dim" -> 2500, "time(s)" -> 0.067774, "num good weights" -> 6, "all unique?" -> True, Column[{<|"w" -> {6, 6, 0, 0}, "ok" -> True, "index" -> 1|>, <|"w" -> {6, 0, 6, 0}, "ok" -> True, "index" -> 307|>, <|"w" -> {6, 0, 0, 6}, "ok" -> True, "index" -> 460|>, <|"w" -> {0, 6, 6, 0}, "ok" -> True, "index" -> 2041|>, <|"w" -> {0, 6, 0, 6}, "ok" -> True, "index" -> 2194|>, <|"w" -> {0, 0, 6, 6}, "ok" -> True, "index" -> 2500|>}]}
```
## 2026-07-24 09:00:12 — 📥 BULK INSERT 3 CELL(S) at positions 42–44
1. [markdown] # Eigenvectors for good (multiplicity-1) weights↵↵A "good" weight (`FindGoodWeight2`) appears with m
2. [code] (* joint spectrum of the global Cartans EMAll[j,j] (diagonal in the GT tensor basis):↵   row i = wei
3. [code] (* example: all good weights of {3,3,0,0}x{3,3,0,0} (2500-dim) -> unique basis states, instant.↵   '


## 2026-07-24 09:00:12
**Input:**
```wolfram
(* joint spectrum of the global Cartans EMAll[j,j] (diagonal in the GT tensor basis):
   row i = weight of the i-th tensor basis state *)
ClearAll[weightList, Psi];
weightList[rep1_, rep2_] := weightList[rep1, rep2] =
   Transpose[Table[Diagonal[Normal[EMAll[j, j, rep1, rep2]]], {j, 4}]];

(* good weight => unique basis state; return it as the (unit) eigenvector.
   $Failed[mult] if the weight is absent or not multiplicity-1 (i.e. not good). *)
Psi[rep1_, rep2_, repall_] := Module[{pos = Position[weightList[rep1, rep2], repall]},
   If[Length[pos] != 1, Return[$Failed["multiplicity=" <> ToString[Length[pos]]]]];
   SparseArray[{pos[[1, 1]] -> 1}, {dim[rep1] dim[rep2]}]];
```
**Output:**
```
(no output)
```

## 2026-07-24 09:00:13
**Input:**
```wolfram
(* example: all good weights of {3,3,0,0}x{3,3,0,0} (2500-dim) -> unique basis states, instant.
   'index' is the position of the nonzero component in the GT tensor basis. *)
Table[{w, FirstPosition[Normal[Psi[{3,3,0,0},{3,3,0,0},w]], 1][[1]]}, {w, FindGoodWeight2[{3,3,0,0},{3,3,0,0}]}] // Column
```
**Output:**
```
Out[151]= \begin{matrix}\left\{\left\{6,6,0,0\right\},1\right\}\\\left\{\left\{6,0,6,0\right\},307\right\}\\\left\{\left\{6,0,0,6\right\},460\right\}\\\left\{\left\{0,6,6,0\right\},2041\right\}\\\left\{\left\{0,6,0,6\right\},2194\right\}\\\left\{\left\{0,0,6,6\right\},2500\right\}\end{matrix}
```

## 2026-07-24 09:09:07
**Input:**
```wolfram
{"FindWeight2 defined?"->(DownValues[FindWeight2]=!={}),
 "FindGoodWeight2 def"->(DownValues[FindGoodWeight2]/.HoldPattern[a_:>b_]:>b),
 "weightList defined?"->(DownValues[weightList]=!={})}
```
**Output:**
```
Out[1]= {"FindWeight2 defined?" -> False, "FindGoodWeight2 def" -> {{rep1 + rep2, rep1 + rep2}}, "weightList defined?" -> True}
```
## 2026-07-24 09:09:59 — ✏️ EDIT CELL 43
(* Shared weight-space primitive: joint spectrum of the global Cartans EMAll[j,j]
   (diagonal in the GT tensor basis). Row i = weight of the i-th tensor basis state.
   This is the pre-Tally body of …


## 2026-07-24 09:09:59
**Input:**
```wolfram
(* Shared weight-space primitive: joint spectrum of the global Cartans EMAll[j,j]
   (diagonal in the GT tensor basis). Row i = weight of the i-th tensor basis state.
   This is the pre-Tally body of the earlier FindWeight2 — kept as one named atom so the
   weight-multiplicity tool and Psi both build on it instead of re-deriving it. *)
ClearAll[weightList, FindWeight2, Psi];
weightList[rep1_, rep2_] := weightList[rep1, rep2] =
   Transpose[Table[Diagonal[Normal[EMAll[j, j, rep1, rep2]]], {j, 4}]];

(* general weight-multiplicity table (restored, now on top of weightList): {weight, mult}.
   FindGoodWeight2 stays the fast Permutations[rep1+rep2] rectangular shortcut from above. *)
FindWeight2[rep1_, rep2_] := Tally[weightList[rep1, rep2]];

(* good weight => unique basis state; return it as the (unit) eigenvector.
   $Failed[mult] if the weight is absent or not multiplicity-1 (i.e. not good). *)
Psi[rep1_, rep2_, repall_] := Module[{pos = Position[weightList[rep1, rep2], repall]},
   If[Length[pos] != 1, Return[$Failed["multiplicity=" <> ToString[Length[pos]]]]];
   SparseArray[{pos[[1, 1]] -> 1}, {dim[rep1] dim[rep2]}]];
```
**Output:**
```
Out= Null
```

## 2026-07-24 09:10:07
**Input:**
```wolfram
(* general good-weights from restored FindWeight2 vs the Permutations shortcut — should match as sets *)
genGood=Cases[FindWeight2[{1,1,0,0},{1,1,0,0}],{w_,1}:>w];
{"general==shortcut"->(Sort[genGood]===Sort[FindGoodWeight2[{1,1,0,0},{1,1,0,0}]]),
 "Psi still works @2500"->(Head[Psi[{3,3,0,0},{3,3,0,0},{6,0,6,0}]]===SparseArray),
 "one weightList def"->(Length[DownValues[weightList]])}
```
**Output:**
```
Out[1]= (suppressed)
Out[2]= {"general==shortcut" -> True, "Psi still works @2500" -> True, "one weightList def" -> 3}
```
