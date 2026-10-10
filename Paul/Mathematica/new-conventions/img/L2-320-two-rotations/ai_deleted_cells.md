## Deleted Cell 42 — 2026-10-09 16:21:52
**Kind:** markdown
```markdown
### Result (2026-10-09): lowest-weight frame and point A for $[3,1,0]$ (is $[4,2,2]$ fully determined?)

**Task 1, LW frame of the point-B family.** $e_{\rm low}$ = GT index 15 (weight $(0,1,3)$); $e_{\rm low}^T E_{ij}=0$ ($i<j$) and $E_{ij}e_{\rm low}=0$ ($i>j$). $v_2=e_1^T\rho(n)=s\,e_{\rm low}^T\rho(m)$ exactly for all six $w_2$ ($m\in N_-$ from $w_0^{-1}n=b_+m$, $\rho(m)=e^{m_{21}E_{21}}e^{m_{31}E_{31}}e^{m_{32}E_{32}}$). 10 LW kets (weights = reverses of the HW weights), $(\rho(m)\otimes\rho(m))LW=LW$, $\langle f_1\otimes v_2|LW\rangle=s\langle\varphi'\otimes e_{\rm low}^T|LW\rangle$ exact. Per-weight ranks of $\varphi'=f_1\rho(m)^{-1}$ (six $w_2$): full except $(2,1,1)$ rank 1 of 2 and $(3,0,1),(3,1,0)$ identically 0 — the exact mirror of the HW frame (which lacks one direction at $(1,1,2)$ and $(0,1,3),(1,0,3)$). Hence the LW frame determines $[5,2,1]^2$ (needs $(1,1,2)$), $[4,3,1]^2$, $[3,3,2]$, $[5,3,0]$, $[4,4,0]$, $[6,2,0]$, $[6,1,1]$, but **not** $[4,2,2]$ (needs $(2,1,1)$, rank 1). LW count torus-completed: 10/10. Combined HW-or-LW verdict: all irreps determined except $[4,2,2]$ (the missing direction is one in each frame, and the two frames do not help each other).

**Task 2, point A** ($y_A=w_2(\theta_{12}-h)/(\theta_{12}-2h)$, rotRanks (6,4)). $v_{2,A}\parallel e_1^T\rho(n_A)$ with $n_A$ the $\mathrm{rOm}$ element at $y_A$ directly (no torus-rotation needed), undo exact on all HW kets. Per-weight ranks of $\varphi_A$: full at all 12 weights (including $(1,1,2)$ rank 2), nothing vanishing; so $A$ alone, and $A\cup B$, give rank 2 at $(1,1,2)$. HW count torus-completed: $A$ 10, $A\cup B$ 10.

**Conclusion.** $[4,2,2]$ is NOT fully determined by point B in either frame (HW frame: rank 1 of 2 at $(1,1,2)$; LW frame: rank 1 of 2 at $(2,1,1)$). It becomes fully determined only through point A (rank 2 at $(1,1,2)$), i.e. conditionally: point A's overlap is not yet known as a Q-expression (it is not on the ladder $\det_2$ curve), so this fix needs the A overlap.
```

## Deleted Cell 41 — 2026-10-09 16:21:52
**Kind:** code
```wolfram
(* result tables for section 14: LW frame (point B), point A alone and A union B, counts, per-irrep verdict *)
u14wtRows = Table[{u14rankTab[[k, 1]], u14rankTab[[k, 2]], u13ranks["a"][[k, 3]], u14rankTab[[k, 3]], u14rkA[[k, 3]], u14rkU[[k, 3]]}, {k, Length[u13wtList]}];
u14rkAt[rk_, μ_] := Cases[rk, {μ, dm_, r_, _} :> {dm, r}][[1]];
u14verdTab = Table[Module[{ν = v["irrep"], μH, μL, A, U}, μH = v["HW site-1 wt"]; μL = v["LW site-1 wt"]; {ν, v["mult"], μH, v["HW dim,rank"], μL, v["LW dim,rank"], If[v["HW determined"], "yes", "no"], If[v["LW determined"], "yes", "no"], If[v["either"], "yes", "NO"], u14rkAt[u14rkA, μH][[2]], u14rkAt[u14rkU, μH][[2]]}], {v, u14verd}];
Grid[{{"(a) per-weight ranks over six w2: {μ, dim V_μ, HW-frame B, LW-frame B, HW-frame A, HW-frame A∪B}", Column[u14wtRows]}, {"(a') per irrep {ν, mult, HW wt, HW (dim,rank), LW wt, LW (dim,rank), HW det., LW det., HW or LW det., rank A at HW wt, rank A∪B at HW wt}", Column[u14verdTab]}, {"(b) point A rotRank (site1, site2); v2_A ∥ e1.rho(n_A) with n_A from rOm at y_A", {Union[checks["u14 point A rotRank (site1, site2) per w2"]], And @@ checks["u14 v2_A parallel to e1.rho(n_A) with n_A = rOm[z,{w1,y_A,w3}] ?"]}}, {"(c) LW count torus-completed (12 random (w2,D)), fixed six w2", {checks["u14 LW count torus-completed (expect 10)"], checks["u14 LW count fixed torus six w2"]}}, {"(c) HW count fixed torus: A (6 rows), A∪B (12 rows)", checks["u14 HW rank fixed torus: A, A+B"]}, {"(c) HW count torus-completed: A, A∪B (B alone: 10)", checks["u14 HW rank torus-completed: A alone, A+B"]}}, Frame -> All, Alignment -> Left]
```

## Deleted Cell 33 — 2026-10-09 16:21:52
**Kind:** markdown
```markdown
### Result (2026-10-09): undoing the site-2 rotation for $[3,1,0]$ (spec `site1-Qform/7-undo-site2-rotation-span.md` §5/§7)

**Mechanism.** With $y^*=w_2(\theta_{12}-h)/(\theta_{12}-3h)$, $v_2=e_1^T\rho(n)$ with $n\in N_+$ the Gauss element of the vacuum, $v_1=e_1^T(\alpha+\beta E_{12})\rho(n)\rho(\Delta)$, $\Delta=\mathrm{diag}(1/(w_2-z_k))$. Since $\rho(n)\otimes\rho(n)$ fixes every highest-weight ket (checked exactly), $\langle v_1\otimes v_2|HW\rangle=s\langle\varphi\otimes e_1^T|HW\rangle$ with $\varphi=v_1\rho(n)^{-1}=e_1^T(\alpha+\beta E_{12})\rho(n\Delta n^{-1})$. The undone 3x3 element $n'=\Delta^{-1}n\Delta n^{-1}$ has $b'\neq0$ exactly because $w_2\neq y^*$, and $b'/(a'c')=(y^*-w_2)/(y^*-z_2)$ (verified symbolically, also with symbolic $z$): the same mechanism as the adjoint. With the resolvent point set to $y^*$ ($b'\equiv0$) and $\beta=0$ the control reproduces the old blind spot.

**Outcome.** (i) Point-B $\beta/\alpha$: $392/2335,\ 64/731,\ -32/2409,\ 8/77,\ 1/23,\ 4312/2693$ at $w_2=7/3,-8/9,4/27,3,-1/2,11/7$. (ii) HW kets of $[3,1,0]^{\otimes2}$: 10 (the $E_{ij}$ with $i<j$ annihilate the ket $e_1$). (iii) $[3,1,0]$ has 12 distinct weights (not 10). Per-weight spans over six $w_2$ (table above): the rotation part alone ($\beta=0$) spans all 12 weight spaces. In the actual family (with the dressing) the components at weights (0,1,3), (1,0,3) vanish identically and (1,1,2) has rank 1 < 2; these do not obstruct the HW count (the weights relevant to the HW kets are (1,1,2) -> $[4,2,2]$ and (0,2,2) -> $[3,3,2]$, still reached through the torus). (iv) HW count: with the torus (12 random $(w_2,D)$) rank **10/10**, all irreps incl. $[4,2,2]$ and $[3,3,2]$ reached, versus 8/10 for the 2026-09-23 families; at fixed torus 7/10 with 12 random $w_2$ (cross-irrep linear relations only: each irrep separately has full rank; the six-$w_2$ rank 6 just saturates at the number of rows). The $\beta=0$ part alone also gives 10 with the torus; the control (resolvent $=y^*$, $\beta=0$) gives 6 and misses $[4,2,2]$, $[3,3,2]$ (zero) and loses one direction in each octet pair — the old blind spot; with the dressing $\beta\neq0$ the control also reaches 10 (the dressing alone rescues it).
```

## Deleted Cell 23 — 2026-10-09 16:21:52
**Kind:** code
```wolfram
(* Review check (main session, independent of the fitted curve): the grading claim c_λ(y) = a_λ y^(6-|λ|) tested at two FRESH y values by direct least-squares fits of f_y (not via the symbolic curve); and a second wrong-y control for the end-to-end identity. *)
revFit[yy_] := LeastSquares[basisT, fL1[yy]];
revC1 = revFit[2]; revC2 = revFit[-5/3];
checks["REVIEW: f_y in span at y = 2, -5/3 (rel residuals)"] = {N[nrmQ[basisT . revC1 - fL1[2]]/nrmQ[fL1[2]], 3], N[nrmQ[basisT . revC2 - fL1[-5/3]]/nrmQ[fL1[-5/3]], 3]};
checks["REVIEW: c_λ(2)/c_λ(-5/3) vs (2/(-5/3))^(6-|λ|) (max rel dev)"] = N[Max[Abs[(revC1/revC2)/((2/(-5/3))^(6 - box)) - 1]], 3];
checks["REVIEW: a_λ from y = 2 vs the rationalized a_λ (max rel dev)"] = N[Max[Abs[revC1/(aLam 2^(6 - box)) - 1]], 3];
checks["REVIEW: end-to-end control with y = 2 (should fail)"] = With[{fW = fL1[2]}, N[constDev[Table[(ovlB[[B]]/x0ovl[[B]])/((det6List[[B]]/psiQ[B]) fW[[B]]), {B, nSt}]], 3]];
Grid[{{"f_y in span at fresh y = 2, -5/3", checks["REVIEW: f_y in span at y = 2, -5/3 (rel residuals)"]}, {"grading ratio test", checks["REVIEW: c_λ(2)/c_λ(-5/3) vs (2/(-5/3))^(6-|λ|) (max rel dev)"]}, {"a_λ from fresh y vs rationalized", checks["REVIEW: a_λ from y = 2 vs the rationalized a_λ (max rel dev)"]}, {"end-to-end control y = 2 (should fail)", checks["REVIEW: end-to-end control with y = 2 (should fail)"]}}, Frame -> All, Alignment -> Left]
```

## Deleted Cell 22 — 2026-10-09 16:21:52
**Kind:** markdown
```markdown
### Q-side result (2026-10-09)

All four tasks of `L2-310-Q-side-spec.md` succeed at P1 ($h=1/2$, $\theta=(2/7,-3/5)$, $z=(3/2,5/4,2/3)$, $w=(5/2,7/3,3/5)$), 60-digit numerics, 225 states. Whole section reruns in about 10 s.

1. **t's, eigenstates, Q's** (cell 16). $t_a(u)$ of degree $(2,4,6)$ from exact 7-point interpolation (0.3 s). The quantum determinant $t_3=\chi_3\nu_1(u)\nu_2(u-h)\nu_3(u-2h)$ holds symbolically. The spectrum of $t_1+\pi t_2$ at $3/7$ is simple: the minimum gap within a sector is $6.4\cdot10^{-3}$. Eigen-residual $7\cdot10^{-54}$ (7 nodes plus the off-node point 9/11). Q's: all six found for **225/225** states. Degrees are exactly $\deg q^j=\mu_j$ and $\deg q_i=6-\mu_i$. Max Baxter residual $1.5\cdot10^{-43}$. Timing: 0.02 s for one state, 3.4 s for all.
2. **Vacuum Q's and det$_6$** (cell 17). $\Omega\otimes\Omega$ is an exact eigen-covector of $t_a(u,G')$ with the interpolated $\tau_a$. The vacuum direct $q$'s have degrees **(0, 4, 6)**. With all 8 poles ($n=0..3$), insertions $1, e^{2\pi ix/h}$, and the **first-guess reference factor** $\Psi_Q=\prod_s[Q^2(\theta_s-3h)Q^3(\theta_s-4h)-Q^3(\theta_s-3h)Q^2(\theta_s-4h)]\,Q^3(\theta_s-h)$: $\langle\Omega|\psi\rangle/\langle x_0|\psi\rangle\propto\det_6[B,F_{vac}]/\Psi_Q$ with max $|\text{ratio}/\text{ratio}_1-1| = 5.1\cdot10^{-45}$. No scan or per-ladder variant was needed. det$_6$ for all states takes 0.3 s.
3. **Ladder-1 object** (cells 18–20). $f_y=\det_{j\in\{3,2\}}[-\langle Q^jE\rangle_1,\langle Q^{j[-]}E\rangle_1]/\text{minor}_1$, with one-site weights on $n=0..3$ and $\text{minor}_1=Q^2(\theta_1-3h)Q^3(\theta_1-4h)-Q^3(\theta_1-3h)Q^2(\theta_1-4h)$, **lies in span$\{1,\tau_\lambda(\theta_1-3h)\}$**. Relative residuals: $8\cdot10^{-48}$, $6\cdot10^{-47}$, $2\cdot10^{-46}$, $8\cdot10^{-48}$ at $y=-6/5, 3, 7/4, -63/43$. Controls: a random combination gives $10^{-54}$, $\langle x_0|\psi\rangle/\langle\Omega|\psi\rangle$ gives 0.91. Basis conditioning: min/max singular value 0.024. $T_{(3,3)}(\theta_1-3h)\equiv0$ as an operator. With symbolic $y$, every $y$-power fits ($5\cdot10^{-46}$) and the curve rationalizes exactly (to $10^{-53}$) to **monomials** $c_\lambda(y)=a_\lambda\,y^{6-|\lambda|}$:
   $c(y)=\big(\tfrac{16}{9}y^6,\ -\tfrac{2240}{387}y^5,\ \tfrac{78400}{15093}y^4,\ \tfrac{19600}{387}y^4,\ -\tfrac{2744000}{45279}y^3,\ \tfrac{96040000}{1946997}y^2,\ \tfrac{2744000}{31347}y^3,\ -\tfrac{48020000}{407511}y^2,\ \tfrac{6722800000}{52568919}y\big)$ on $(1,T_{10},T_{11},T_{20},T_{21},T_{22},T_{30},T_{31},T_{32})$.
   The only nonzero $y$ with $c(y)\parallel$ `cP3` is **$y=-63/43=y^*$** (exact `Solve`), and there $c(y^*)/c(y^*)_{T_{20}}=$ `cP3` exactly, with $\kappa=c(y^*)_{T_{20}}=34306448400/147008443$. No $y$ gives $c(y)\parallel$ `cP2`, so point A is not on the curve: its $c$ vanishes on $1, T_{10}, T_{20}, T_{30}$, whereas every entry of $c(y)$ is nonzero.
4. **End-to-end, point B** (cell 21). $\Phi_B=\sum_k$ `cP3`$_k\langle\Omega|O_k$ is a product covector (rank 1, exact), and
   $$\frac{\langle v_1\otimes v_2|\Psi_B\rangle}{\langle x_0|\Psi_B\rangle}\ \propto\ \frac{\det_6[Q_B,Q^{vac}_{(w_1,w_2)}]}{\Psi_Q(x_0)}\times\frac{\det_2^{(1)}[Q_B;\,y^*]}{\text{minor}_1},\qquad y^*=w_2\frac{\theta_{12}-h}{\theta_{12}-3h}=-\tfrac{63}{43},$$
   with max $|\text{ratio}/\text{ratio}_1-1|=5.1\cdot10^{-45}$ over the 225 states. ratio4/ratio6 $=1/\kappa$ to $5\cdot10^{-47}$. Control: $y=-63/8$ gives 77.

**Structure.** The adjoint's two-rotation formula carries over to $[3,1,0]$ with only the rep-generic shifts: poles $n=0..\lambda_1$, vacuum triple $(\theta-4h,\theta-3h,\theta-h)$, $\text{minor}_1$ at $(\theta_1-3h,\theta_1-4h)$, insertion point $\theta_1-\lambda_1h$, and $y^*$ with $\theta_{12}-\lambda_1h$. This holds even though site 1 is a *dressed* vacuum and not a rotation. New in the Q-language: the ladder-1 det$_2$ curve is diagonal in the box grading, $c_\lambda\propto y^{-|\lambda|}$. Single numeric point only; the other two points of cell 12 are not redone on the Q side.
```

## Deleted Cell 21 — 2026-10-09 16:21:52
**Kind:** code
```wolfram
(* Task 4. End-to-end for point B: Φ_B = Σ_k cP3_k <Ω|O_k> (O = 1, T_λ(θ1-3h); allMats of cell 6), a product covector v1⊗v2 (rank 1 as a 15x15 matrix, exact). Claim: <Φ_B|ψ>/<x0|ψ> ∝ det6[B, vac (w1,w2)]/Ψ_Q(x0) × f_{y*}(B), y* = w2(θ12-h)/(θ12-3h) = -63/43, state-independent over the 225 states *)
ΦB = Flatten[Sum[cP3[[k]] allMats[[k]], {k, Length[allMats]}]];
checks["Q4: Φ_B is a product covector (rank of 15x15 matrix, exact)"] = MatrixRank[Partition[ΦB, d1]];
yStarB = wN[[2]] (θ12 - hN)/(θ12 - 3 hN);
fB = fL1[yStarB];
ovlB = Table[ΦB . vecSt[[B]], {B, nSt}];
ratio4 = Table[(ovlB[[B]]/x0ovl[[B]])/((det6List[[B]]/psiQ[B]) fB[[B]]), {B, nSt}];
checks["Q4: y* = w2(θ12-h)/(θ12-3h)"] = yStarB;
checks["Q4: MAIN <Φ_B|ψ>/<x0|ψ> ∝ det6(w)/Ψ_Q × f_{y*} (max |ratio/ratio[1] - 1|, 225 states)"] = sizeQ[constDev[ratio4]];
checks["Q4: the constant ratio4[1]; ratio4[1]/ratio6[1]"] = N[{ratio4[[1]], ratio4[[1]]/ratio6[[1]]}, 12];
(* consistency with Task 3: ratio4/ratio6 should be 1/κ with κ = c(y★)[T20] (cP3 has T20 entry 1) *)
checks["Q4: ratio4/ratio6 = 1/κ (κ from the exact curve)"] = sizeQ[Max[Abs[ratio4/ratio6 κ0 - 1]] /. κ0 -> checks["Q3x: c(-63/43) = κ cP3 exactly; κ"][[2]]];
(* control: a wrong y in f breaks state independence *)
checks["Q4: control with y = -63/8 (should fail)"] = With[{fA = fL1[-63/8]}, N[constDev[Table[(ovlB[[B]]/x0ovl[[B]])/((det6List[[B]]/psiQ[B]) fA[[B]]), {B, nSt}]], 3]];
Grid[{{"Φ_B rank (1 = product)", checks["Q4: Φ_B is a product covector (rank of 15x15 matrix, exact)"]}, {"y*", yStarB}, {"MAIN: max |ratio4/ratio4[1] - 1|", checks["Q4: MAIN <Φ_B|ψ>/<x0|ψ> ∝ det6(w)/Ψ_Q × f_{y*} (max |ratio/ratio[1] - 1|, 225 states)"]}, {"ratio4[1], ratio4[1]/ratio6[1]", checks["Q4: the constant ratio4[1]; ratio4[1]/ratio6[1]"]}, {"ratio4/ratio6 = 1/κ", checks["Q4: ratio4/ratio6 = 1/κ (κ from the exact curve)"]}, {"control y = -63/8", checks["Q4: control with y = -63/8 (should fail)"]}}, Frame -> All, Alignment -> Left]
```

## Deleted Cell 20 — 2026-10-09 16:21:52
**Kind:** code
```wolfram
(* Task 3, exact form of the curve. The fitted y-power coefficients rationalize (to 1e-53); every entry of c(y) is a single monomial y^(6-|λ|) (|λ| = λ1 + λ2, |1| = 0), so the curve is c_λ(y) = a_λ y^(6-|λ|): the insertion y^n acts as a grading by the box number. The roots y = 0 of the previous cell are the trivial point c(0) = 0; here the parallelism with cP3 / cP2 is solved exactly on the curve *)
cPowEx = Rationalize[Chop[cPow, 10^-40], 10^-35];
cYex = Factor /@ (ySym^Range[0, degY] . cPowEx);
box = Join[{0}, Total /@ Most[shapes]];
aLam = Table[Coefficient[cYex[[k]], ySym, 6 - box[[k]]], {k, 9}];
checks["Q3x: rationalized curve reproduces the fit (max |Δ| of y-power coefficients)"] = sizeQ[Max[Abs[cPowEx - cPow]]];
checks["Q3x: c_λ(y) = a_λ y^(6-|λ|) exactly"] = (cYex === aLam ySym^(6 - box));
(* exact parallelism with the product points on the curve *)
parY[ref_] := If[# === {}, "none", ySym /. #] &@Solve[Join[Thread[Flatten[Minors[{cYex/ySym, ref}, 2]] == 0], {ySym != 0}], ySym];
checks["Q3x: y with c(y) ∥ cP3 (exact Solve, y ≠ 0)"] = parY[cP3[[1 ;; 9]]];
checks["Q3x: y with c(y) ∥ cP2 (exact Solve, y ≠ 0)"] = parY[cP2[[1 ;; 9]]];
checks["Q3x: c(-63/43) = κ cP3 exactly; κ"] = With[{cv = cYex /. ySym -> -63/43}, {cv/cv[[4]] === cP3[[1 ;; 9]], cv[[4]]}];
(* a_λ relative to cP3 at y*: cP3_λ / (a_λ y*^(6-|λ|)) should be one constant *)
checks["Q3x: cP3_λ/(a_λ y*^(6-|λ|)) (one constant)"] = Union[cP3[[1 ;; 9]]/(aLam (-63/43)^(6 - box))];
Grid[{{"fit vs rationalized", checks["Q3x: rationalized curve reproduces the fit (max |Δ| of y-power coefficients)"]}, {"c_λ(y) = a_λ y^(6-|λ|)", checks["Q3x: c_λ(y) = a_λ y^(6-|λ|) exactly"]}, {"c(y) (1, T10, T11, T20, T21, T22, T30, T31, T32)", cYex}, {"a_λ factored", FactorInteger /@ aLam}, {"nonzero y with c(y) ∥ cP3", checks["Q3x: y with c(y) ∥ cP3 (exact Solve, y ≠ 0)"]}, {"nonzero y with c(y) ∥ cP2", checks["Q3x: y with c(y) ∥ cP2 (exact Solve, y ≠ 0)"]}, {"c(-63/43)/c[T20] == cP3 exactly; c(-63/43)[T20]", checks["Q3x: c(-63/43) = κ cP3 exactly; κ"]}, {"cP3/(a_λ y*^(6-|λ|))", checks["Q3x: cP3_λ/(a_λ y*^(6-|λ|)) (one constant)"]}}, Frame -> All, Alignment -> Left]
```

## Deleted Cell 19 — 2026-10-09 16:21:52
**Kind:** code
```wolfram
(* Task 3, symbolic y (adjoint cell 94): E = y^n at x = θ1 - n h, so f_y(B) is a polynomial in y; fit each y-power coefficient vector to span{1, τ_λ} and assemble the coefficient curve c(y) (order (1, T10, T11, T20, T21, T22, T30, T31, T32)). Then: where is c(y) parallel to cP3 (point B) resp. cP2 (point A)? ySym is a fresh symbol (y, yS, yA are taken) *)
su2det1S[B_] := Det[Table[{-br1[Function[x, QU[B, j, x] ySym^(-(x - th1)/hN)]], br1[Function[x, QU[B, j, x - hN] ySym^(-(x - th1)/hN)]]}, {j, {3, 2}}]];
{fySTime, fyS} = AbsoluteTiming[Table[Expand[su2det1S[B]/minor1[B]], {B, nSt}]];
degY = Max[Exponent[fyS, ySym]]; minDegY = Min[Table[Exponent[f, ySym, Min], {f, fyS}]];
coefY[p_] := Coefficient[fyS, ySym, p];
cPow = Table[LeastSquares[basisT, coefY[p]], {p, 0, degY}];
checks["Q3s: every y-power coefficient of f_y in span{1, τ} (max rel residual; zero powers skipped)"] = sizeQ[Max[Table[If[nrmQ[coefY[p]] < 10^-40, 0, nrmQ[basisT . cPow[[p + 1]] - coefY[p]]/nrmQ[coefY[p]]], {p, 0, degY}]]];
cY = Chop[Sum[ySym^p cPow[[p + 1]], {p, 0, degY}], 10^-45];
(* parallelism with the product points: all 2x2 minors of {c(y), c_point} *)
yB = -63/43; yAp = -63/8;
minorsAt[cv_, ref_, yy_] := Max[Abs[Flatten[Minors[{cv /. ySym -> yy, ref}, 2]]]]/(nrmQ[cv /. ySym -> yy] nrmQ[ref]);
checks["Q3s: c(y) ∥ cP3 at y = -63/43 (max |2x2 minor|, normalized)"] = sizeQ[minorsAt[cY, cP3[[1 ;; 9]], yB]];
checks["Q3s: c(y) ∥ cP2 at y = -63/8 (max |2x2 minor|, normalized)"] = sizeQ[minorsAt[cY, cP2[[1 ;; 9]], yAp]];
(* the c-vector at y*, normalized like cP3 (T20 entry = 1), recognized as rationals *)
cAtB = (cY /. ySym -> yB)/((cY /. ySym -> yB)[[4]]);
checks["Q3s: c(y*)/c(y*)[T20] rationalized == cP3[[1;;9]]"] = (Rationalize[Chop[cAtB, 10^-40], 10^-40] === cP3[[1 ;; 9]]);
(* all roots y where c(y) is parallel to cP3 (resp. cP2): common roots of the minors *)
commonRoots[ref_] := Module[{mins = DeleteCases[Chop[Expand[Flatten[Minors[{cY, ref}, 2]]], 10^-40], 0], r}, r = ySym /. NSolve[mins[[1]] == 0, ySym, WorkingPrecision -> 40]; Select[r, Max[Abs[mins /. ySym -> #]]/Max[Abs[CoefficientList[mins[[1]], ySym]]] < 10^-25 &]];
checks["Q3s: common roots of c(y) ∥ cP3"] = N[commonRoots[cP3[[1 ;; 9]]], 10];
checks["Q3s: common roots of c(y) ∥ cP2"] = N[commonRoots[cP2[[1 ;; 9]]], 10];
Grid[{{"f_y symbolic: time (s), y-degree range", {fySTime, {minDegY, degY}}}, {"max rel residual of the y-power fits", checks["Q3s: every y-power coefficient of f_y in span{1, τ} (max rel residual; zero powers skipped)"]}, {"c(y) (order 1, T10, T11, T20, T21, T22, T30, T31, T32), 8 digits", Column[N[Expand[cY], 8]]}, {"c(y) ∥ cP3 at y* = -63/43", checks["Q3s: c(y) ∥ cP3 at y = -63/43 (max |2x2 minor|, normalized)"]}, {"c(y*) normalized and rationalized = cP3?", checks["Q3s: c(y*)/c(y*)[T20] rationalized == cP3[[1;;9]]"]}, {"c(y) ∥ cP2 at y = -63/8", checks["Q3s: c(y) ∥ cP2 at y = -63/8 (max |2x2 minor|, normalized)"]}, {"all y with c(y) ∥ cP3; with c(y) ∥ cP2", {checks["Q3s: common roots of c(y) ∥ cP3"], checks["Q3s: common roots of c(y) ∥ cP2"]}}}, Frame -> All, Alignment -> Left]
```

## Deleted Cell 18 — 2026-10-09 16:21:52
**Kind:** code
```wolfram
(* Task 3. The ladder-1 object: one-site weights wgt1(n) = 1/Π_{m≠n}(m-n)h on n = 0..3 (adjoint cell 68), E(x) = y^(-(x-θ1)/h), f_y(B) = det_{j∈{3,2}}[-<Q^j E>1, <Q^j[-] E>1]/minor1(B) with minor1 = Q2(θ1-3h)Q3(θ1-4h) - Q3(θ1-3h)Q2(θ1-4h) (the ladder-1 part of the reference factor that worked in Task 2). Span test against {1, τ_λ(θ1-3h)} over the 225 states, then symbolic y and the coefficient curve c(y) *)
KeyDropFrom[checks, "Q3: max |τ_(3,3)| over the states"];
poles1 = Range[Λ[[3]], Λ[[1]]]; wgt1[n_] := 1/Product[(m - n) hN, {m, DeleteCases[poles1, n]}];
th1 = thN[[1]];
br1[f_] := Sum[wgt1[n] f[th1 - n hN], {n, poles1}];
su2det1[B_, yy_] := Module[{E1f = Function[x, yy^(-(x - th1)/hN)]}, Det[Table[{-br1[Function[x, QU[B, j, x] E1f[x]]], br1[Function[x, QU[B, j, x - hN] E1f[x]]]}, {j, {3, 2}}]]];
minor1[B_] := QU[B, 2, th1 - 3 hN] QU[B, 3, th1 - 4 hN] - QU[B, 3, th1 - 3 hN] QU[B, 2, th1 - 4 hN];
fL1[yy_] := Table[su2det1[B, N[yy, precQ]]/minor1[B], {B, nSt}];
(* T-eigenvalues at v0 = θ1-3h on the eigenstates (largest component) and their residual *)
Tv0N = Association[Table[λ -> N[Tv0[λ], precQ], {λ, shapes}]];
τSt = Association[Table[λ -> Table[With[{k = Ordering[Abs[vecSt[[B]]], -1][[1]]}, (Tv0N[λ] . vecSt[[B]])[[k]]/vecSt[[B, k]]], {B, nSt}], {λ, shapes}]];
checks["Q3: T_λ(θ1-3h) eigen-residual on the 225 states"] = sizeQ[Max[Table[Max[Abs[Tv0N[λ] . vecSt[[B]] - τSt[λ][[B]] vecSt[[B]]]], {λ, shapes}, {B, nSt}]]];
checks["Q3: T_(3,3)(θ1-3h) is the zero operator (exact)"] = SubsetQ[{0}, Union[Tv0[{3, 3}]["NonzeroValues"]]];
ones = ConstantArray[1, nSt];
basisT = Transpose[Join[{ones}, Table[τSt[λ], {λ, Most[shapes]}]]];   (* order (1, T10, T11, T20, T21, T22, T30, T31, T32) *)
checks["Q3: basis {1, τ_λ} (9 functions) conditioning: min sv/max sv"] = With[{sv = SingularValueList[basisT, Tolerance -> 0]}, N[{Length[sv], Last[sv]/First[sv]}, 3]];
(* Norm of an all-zero 60-digit vector returns an exact 0, so use the explicit 2-norm *)
nrmQ[v_] := Sqrt[Total[Abs[Flatten[v]]^2]];
spanRes[f_, A_] := Module[{cf = LeastSquares[A, f]}, nrmQ[A . cf - f]/nrmQ[f]];
yTest = {-6/5, 3, 7/4, -63/43};
{fTime, spanTab} = AbsoluteTiming[Table[{yv0, sizeQ[spanRes[fL1[yv0], basisT]]}, {yv0, yTest}]];
checks["Q3: f_y in span{1, τ_λ}: rel residual at y = -6/5, 3, 7/4, -63/43"] = spanTab;
(* control: a random combination of the basis is fitted to 1e-50, a generic function (here <x0|ψ>) is not *)
checks["Q3: control rel residuals (random combination, <x0|ψ>/<Ω|ψ>)"] = {sizeQ[spanRes[basisT . RandomReal[{-1, 1}, 9, WorkingPrecision -> precQ], basisT]], sizeQ[spanRes[x0ovl/ovl, basisT]]};
Grid[Join[{{"T eigen-residual", checks["Q3: T_λ(θ1-3h) eigen-residual on the 225 states"]}, {"T_(3,3)(θ1-3h) = 0 identically (exact)", checks["Q3: T_(3,3)(θ1-3h) is the zero operator (exact)"]}, {"basis {1, τ_λ}: # sv, min/max sv", checks["Q3: basis {1, τ_λ} (9 functions) conditioning: min sv/max sv"]}, {"controls (random combo, <x0|ψ>/<Ω|ψ>)", checks["Q3: control rel residuals (random combination, <x0|ψ>/<Ω|ψ>)"]}, {"span test: {y, rel residual}; time (s)", fTime}}, spanTab], Frame -> All, Alignment -> Left]
```

## Deleted Cell 17 — 2026-10-09 16:21:52
**Kind:** code
```wolfram
(* Task 2. Vacuum Q's of the bra twist w (direct Baxter equation with w_i in place of z_i, eigenvalues of t_a(G') on Ω⊗Ω interpolated exactly) and the vacuum overlap determinant: det6 with all 8 poles {s, n}, n = λ3..λ1 = 0..3, insertions 1 and e^{2πix/h} (adjoint cells 7-8); reference factor = first guess with the vacuum SoV triple (θ-4h, θ-3h, θ-h) *)
tWp = Table[Sum[lagQ[m] tL2[a, uPts[[m]], GP], {m, 7}], {a, 3}];
kOm = Ordering[Abs[Ω2], -1][[1]];
tauVac = Table[Expand[Sum[lagQ[m] (Ω2 . tL2[a, uPts[[m]], GP])[[kOm]]/Ω2[[kOm]], {m, 7}]], {a, 3}];
checks["Q2: Ω⊗Ω eigen-covector of t_a(u, G') with the interpolated τ_a (exact, nodes and 2/9)"] = And @@ Flatten[Table[SubsetQ[{0}, Union[Normal[Ω2 . tL2[a, x, GP] - (tauVac[[a]] /. uQ -> x) Ω2]]], {a, 3}, {x, Append[uPts, 2/9]}]];
qdVac = Table[solveQ[opDir, N[tauVac, precQ], wN[[i]], -1], {i, 3}];
checks["Q2: vacuum direct q's degrees (expect 0,4,6)"] = If[FreeQ[qdVac, {}], Exponent[qdVac, xQ], "not found"];
checks["Q2: vacuum Baxter residual"] = sizeQ[Max[Table[qResid[opDir, N[tauVac, precQ], wN[[i]], -1, qdVac[[i]]], {i, 3}]]];
(* full Q's at rational points, memoized; A = state index, "vac" = the w-vacuum *)
QU[A_, j_, x_] := QU[A, j, x] = N[zN[[j]]^(x/hN), precQ] (quS[[A, j]] /. xQ -> x);
QD[A_Integer, i_, x_] := QD[A, i, x] = N[zN[[i]]^(-x/hN), precQ] (qdS[[A, i]] /. xQ -> x);
QD["vac", i_, x_] := QD["vac", i, x] = N[wN[[i]]^(-x/hN), precQ] (qdVac[[i]] /. xQ -> x);
poles = Flatten[Table[{s, n}, {s, nL}, {n, Λ[[3]], Λ[[1]]}], 1];
xP[{s_, n_}] := thN[[s]] - n hN;
wgt[p_] := wgt[p] = 1/Product[xP[p] - xP[q], {q, DeleteCases[poles, p]}];
ins[1][x_] := 1; ins[2][x_] := Exp[2 Pi I x/hN];
prsL1 = {{1, 3}, {1, 2}, {2, 3}};
rowG[B_, F_, {i_, j_}, k_] := Sum[N[wgt[p] ins[k][xP[p]], precQ] With[{x = xP[p]}, {-QU[B, j, x] F[i, x], -x QU[B, j, x] F[i, x], QU[B, j, x] F[i, x + hN], x QU[B, j, x] F[i, x + hN], QU[B, j, x - hN] F[i, x], x QU[B, j, x - hN] F[i, x]}], {p, poles}];
det6[B_, F_] := Det[Flatten[Table[rowG[B, F, ij, k], {k, 2}, {ij, prsL1}], 1]];
Fvac = Function[{i, x}, QD["vac", i, x]];
(* reference factor Ψ_Q(x0) with a general triple (x1, x21, x22) = θ_s - (n1, n21, n22) h; first guess (4, 3, 1) *)
psiQg[B_, {n1_, n21_, n22_}] := Product[With[{t = thN[[s]]}, (QU[B, 2, t - n21 hN] QU[B, 3, t - n1 hN] - QU[B, 3, t - n21 hN] QU[B, 2, t - n1 hN]) QU[B, 3, t - n22 hN]], {s, nL}];
psiQ[B_] := psiQg[B, {4, 3, 1}];
ovl = Table[Ω2 . vecSt[[B]], {B, nSt}]; x0ovl = Table[x02 . vecSt[[B]], {B, nSt}];
constDev[l_List] := Max[Abs[l/l[[1]] - 1]];
{d6Time1, d6one} = AbsoluteTiming[det6[1, Fvac]];
{d6Time, det6List} = AbsoluteTiming[Table[det6[B, Fvac], {B, nSt}]];
ratio6 = Table[(ovl[[B]]/x0ovl[[B]])/(det6List[[B]]/psiQ[B]), {B, nSt}];
checks["Q2: min |<x0|ψ>|, min |<Ω|ψ>|, min |det6|, min |Ψ_Q| (nonzero, not zero tests)"] = N[{Min[Abs[x0ovl]], Min[Abs[ovl]], Min[Abs[det6List]], Min[Abs[psiQ /@ Range[nSt]]]}, 3];
checks["Q2: 6x6: <Ω|ψ>/<x0|ψ> ∝ det6[B, vac]/Ψ_Q(x0) (max |ratio/ratio[1] - 1|)"] = sizeQ[constDev[ratio6]];
checks["Q2: the constant ratio6[1]"] = N[ratio6[[1]], 12];
Grid[{{"Ω⊗Ω eigen-covector (exact)", checks["Q2: Ω⊗Ω eigen-covector of t_a(u, G') with the interpolated τ_a (exact, nodes and 2/9)"]}, {"vacuum q degrees", checks["Q2: vacuum direct q's degrees (expect 0,4,6)"]}, {"vacuum Baxter residual", checks["Q2: vacuum Baxter residual"]}, {"det6 time: one state, 225 states (s)", {d6Time1, d6Time}}, {"min |<x0|ψ>|, |<Ω|ψ>|, |det6|, |Ψ_Q|", checks["Q2: min |<x0|ψ>|, min |<Ω|ψ>|, min |det6|, min |Ψ_Q| (nonzero, not zero tests)"]}, {"max |ratio6/ratio6[1] - 1| (225 states)", checks["Q2: 6x6: <Ω|ψ>/<x0|ψ> ∝ det6[B, vac]/Ψ_Q(x0) (max |ratio/ratio[1] - 1|)"]}, {"ratio6[1]", checks["Q2: the constant ratio6[1]"]}}, Frame -> All, Alignment -> Left]
```

## Deleted Cell 16 — 2026-10-09 16:21:52
**Kind:** code
```wolfram
(* Task 1. Polynomial t_a(u) (degree 2a) by exact Lagrange interpolation on 7 points; eigenstates of t1 + π t2 at u = 3/7 per weight sector at 60 digits; eigenvalue polynomials from the 7 point values (largest component); then minimal-degree monic q's of the direct and dual Baxter equations (operators verbatim from the adjoint notebook, cell 5) *)
precQ = 60; uPts = 1/7 + Range[0, 6]/3;
(* size of a 60-digit residual: a numerically exact zero is reported by its accuracy bound 10^-Accuracy *)
sizeQ[x_] := If[x == 0, 10.^-Accuracy[x], N[Abs[x], 3]];
lagQ[k_] := Product[(uQ - uPts[[l]])/(uPts[[k]] - uPts[[l]]), {l, DeleteCases[Range[7], k]}];
{tpTime, tZp} = AbsoluteTiming[Table[Sum[lagQ[m] tL2[a, uPts[[m]], Zd], {m, 7}], {a, 3}]];
checks["Q1: qdet t3(u) = χ3 ν1(u)ν2(u-h)ν3(u-2h) (symbolic, exact)"] = SubsetQ[{0}, Union[Expand[(tZp[[3]] - χ3z nuP[1][uQ] nuP[2][uQ - hN] nuP[3][uQ - 2 hN] id2)["NonzeroValues"]]]];
checks["Q1: degrees of t_a(u) in u"] = Table[Max[Exponent[Expand[tZp[[a]]["NonzeroValues"]], uQ]], {a, 3}];
(* eigenstates *)
mixQ = N[tL2[1, 3/7, Zd] + Pi tL2[2, 3/7, Zd], precQ];
sectQ = GatherBy[Range[d2], wtBasis[[#]] &];
{eigTime, stQ} = AbsoluteTiming[Flatten[Table[With[{es = Eigensystem[Normal[mixQ[[sec, sec]]]]}, Table[<|"mu" -> wtBasis[[sec[[1]]]], "gap" -> If[Length[sec] > 1, Min[Abs[Delete[es[[1]], r] - es[[1, r]]]], Infinity], "vec" -> Normal[SparseArray[Thread[sec -> es[[2, r]]], {d2}]]|>, {r, Length[sec]}]], {sec, sectQ}], 1]];
vecSt = stQ[[All, "vec"]]; muSt = stQ[[All, "mu"]]; nSt = Length[stQ];
tNum[a_, x_] := tNum[a, x] = N[tL2[a, x, Zd], precQ];
{valTime, valsSt} = AbsoluteTiming[Table[With[{k = Ordering[Abs[ψ], -1][[1]]}, Table[Expand[Sum[lagQ[m] (tNum[a, uPts[[m]]] . ψ)[[k]]/ψ[[k]], {m, 7}]], {a, 3}]], {ψ, vecSt}]];
(* residual at the 7 nodes and at the off-node point 9/11 (tests the polynomial, not just the node values) *)
checks["Q1: eigen-residual max |t_a(x)ψ - t_a^B(x)ψ| (x = nodes and 9/11)"] = sizeQ[Max[Table[Max[Abs[tNum[a, x] . vecSt[[A]] - (valsSt[[A, a]] /. uQ -> x) vecSt[[A]]]], {A, nSt}, {a, 3}, {x, Append[uPts, 9/11]}]]];
checks["Q1: simple spectrum, min eigenvalue gap of t1+πt2 within a sector"] = N[Min[stQ[[All, "gap"]]], 3];
(* Q-functions *)
opDir[tt_, Q_, x_] := nuP[1][x] nuP[1][x + hN] nuP[1][x + 2 hN] Q[x - hN] - (tt[[1]] /. uQ -> x) nuP[1][x + hN] nuP[1][x + 2 hN] Q[x] + (tt[[2]] /. uQ -> x + hN) nuP[1][x + 2 hN] Q[x + hN] - (tt[[3]] /. uQ -> x + 2 hN) Q[x + 2 hN];
opDual[tt_, Q_, x_] := nuP[3][x] nuP[3][x - hN] nuP[3][x - 2 hN] Q[x + hN] - (tt[[1]] /. uQ -> x) nuP[3][x - hN] nuP[3][x - 2 hN] Q[x] + (tt[[2]] /. uQ -> x) nuP[3][x - 2 hN] Q[x - hN] - (tt[[3]] /. uQ -> x) Q[x - 2 hN];
(* minimal-degree monic q with Q = zz^(s x/h) q solving op; {} if none up to degree 8 *)
solveQ[op_, tt_, zz_, s_] := Module[{cq, eqs, ns, res = {}}, Do[If[res === {}, eqs = CoefficientList[Chop[Expand[op[tt, Function[yy, zz^(s (yy - xQ)/hN) Sum[cq[j] yy^j, {j, 0, d}]], xQ]], 10^-40], xQ]; ns = If[eqs === {}, {UnitVector[d + 1, d + 1]}, NullSpace[N[Table[Coefficient[e, cq[j]], {e, eqs}, {j, 0, d}], precQ], Tolerance -> 10^-30]]; If[Length[ns] > 0, res = Expand[#/Coefficient[#, xQ, d]] &[Sum[ns[[1, j + 1]] xQ^j, {j, 0, d}]]]], {d, 0, 8}]; res];
qResid[op_, tt_, zz_, s_, q_] := Max[Abs[CoefficientList[Expand[op[tt, Function[yy, zz^(s (yy - xQ)/hN) (q /. xQ -> yy)], xQ]], xQ]]];
{q1Time, q1Deg} = AbsoluteTiming[{Exponent[Table[solveQ[opDir, valsSt[[1]], zN[[i]], -1], {i, 3}], xQ], Exponent[Table[solveQ[opDual, valsSt[[1]], zN[[j]], 1], {j, 3}], xQ]}];
{qTime, {qdS, quS}} = AbsoluteTiming[{Table[solveQ[opDir, valsSt[[A]], zN[[i]], -1], {A, nSt}, {i, 3}], Table[solveQ[opDual, valsSt[[A]], zN[[j]], 1], {A, nSt}, {j, 3}]}];
notFound = Select[Range[nSt], ! FreeQ[{qdS[[#]], quS[[#]]}, {}] &];
checks["Q1: # states with all six Q's found"] = nSt - Length[notFound];
checks["Q1: deg q^j = μ_j for all states"] = And @@ Flatten[Table[Exponent[quS[[A, j]], xQ] == muSt[[A, j]], {A, nSt}, {j, 3}]];
checks["Q1: deg q_i = 6 - μ_i for all states"] = And @@ Flatten[Table[Exponent[qdS[[A, i]], xQ] == 6 - muSt[[A, i]], {A, nSt}, {i, 3}]];
checks["Q1: max Baxter residual (coefficients, monic q)"] = sizeQ[Max[Table[{qResid[opDir, valsSt[[A]], zN[[i]], -1, qdS[[A, i]]], qResid[opDual, valsSt[[A]], zN[[i]], 1, quS[[A, i]]]}, {A, nSt}, {i, 3}]]];
Grid[{{"t_a(u) interpolation time (s)", tpTime}, {"qdet symbolic exact", checks["Q1: qdet t3(u) = χ3 ν1(u)ν2(u-h)ν3(u-2h) (symbolic, exact)"]}, {"deg_u t_a", checks["Q1: degrees of t_a(u) in u"]}, {"# states, eigensystem time (s), eigenvalue-polynomial time (s)", {nSt, eigTime, valTime}}, {"eigen-residual (60 digits)", checks["Q1: eigen-residual max |t_a(x)ψ - t_a^B(x)ψ| (x = nodes and 9/11)"]}, {"min sector gap (simple spectrum)", checks["Q1: simple spectrum, min eigenvalue gap of t1+πt2 within a sector"]}, {"state 1 (μ = " <> ToString[muSt[[1]]] <> "): deg q_i, deg q^j; time (s)", {q1Deg, q1Time}}, {"all 225 states: Q time (s)", qTime}, {"# states with all six Q's", checks["Q1: # states with all six Q's found"]}, {"deg q^j = μ_j / deg q_i = 6 - μ_i", {checks["Q1: deg q^j = μ_j for all states"], checks["Q1: deg q_i = 6 - μ_i for all states"]}}, {"max Baxter residual", checks["Q1: max Baxter residual (coefficients, monic q)"]}}, Frame -> All, Alignment -> Left]
```

## Deleted Cell 15 — 2026-10-09 16:21:52
**Kind:** markdown
```markdown
## Q-side

The Q-function dressing of point B (spec: `L2-310-Q-side-spec.md`). Patterns copied from `L2-rotation-sov.wb` (cells 3–8, 27, 68, 87, 94, 102) with the rep-generic pole set $n=\lambda_3..\lambda_1=0..3$. New global names: `uQ`, `xQ` (formal variables), `uPts`, `precQ`, `lagQ`, `tNum`, `tZp`, `mixQ`, `sectQ`, `stQ`, `vecSt`, `muSt`, `nSt`, `valsSt`, `opDir`, `opDual`, `solveQ`, `qdS`, `quS`, `QU`, `QD`, … (none used in cells 1–14; `u`, `yS`, `yA` hold values from the earlier cells and are avoided).
```

## Deleted Cell 14 — 2026-10-09 16:21:52
**Kind:** markdown
```markdown
### Result (2026-10-09): the §10 recipe for $[3,1,0]$ — same mechanism, two isolated product points, but the site-1 factor is a *dressed* vacuum (vacuum + first descendant), not a rotation of the vacuum

1. **Localization at $\theta_1-3h$** (cell 5): $(1,1),(2,2),(3,2)$ are site-1 local (as in `L2-310-localization.md`); every other shape with $\lambda_1\le3$ has rank 2 with $\Omega_1$ inside its site-2 span; $(3,3)$ vanishes on $\langle\Omega|$. The joint site-2 span of all shapes is **2-dimensional**, $\{\Omega_1,\omega_2\}$, with $\omega_2$ a rotation and the pencil $\Omega_1+\kappa\omega_2$ inside the rotation orbit (cell 6): exactly the adjoint's structure. Shapes with $\lambda_1=4$ add nothing (cell 8, S1).
2. **Product locus** in the 10-dim span $\langle\Omega|(c_0+\sum_\lambda c_\lambda T_\lambda(\theta_1-3h))$ (105 quadrics, exact `Solve` 0.5 s): the localizing span (5 free parameters) and two isolated points.
   - **Point A**: $c$ = the adjoint's $c^*=(-2187/1600,-243/80,-27/43,1,320/1161)$ carried verbatim to the shapes shifted by one column, $(1,1),(2,1),(2,2),(3,1),(3,2)$. Site 2 is a torus-rotated vacuum with invariant $y=w_2(\theta_{12}-h)/(\theta_{12}-2h)$ (the adjoint's $y^*$) at all three numeric points (cell 12); site 1 has `rotRank` 6 (not a rotation; not in the vacuum/descendant spans tested, cells 9, 13).
   - **Point B** (all nine shapes and $1$; at P1 $c=(81/1075,36/215,4/39,1,860/1053,4300/9477,-860/729,-92450/85293,-1849000/2302911)$ on $(1,T_{10},T_{11},T_{20},T_{21},T_{22},T_{30},T_{31},T_{32})$): site 2 is a torus-rotated vacuum with invariant
   $$y^*=w_2\,\frac{\theta_{12}-h}{\theta_{12}-3h}=w_2\,\frac{\theta_{12}-h}{\theta_{12}-\lambda_1h}$$
   at all three numeric points (cell 12), and site 1 is
   $$v_1\ \propto\ e_1^T(\alpha+\beta E_{12})\,R\big(z,(w_1,y^*)\big)\,\rho\big((w_2-G^T)^{-1}\big),$$
   the adjoint's $v_1=\Omega_{(w_1,y^*)}\rho(D_{w_2})$ with the vacuum covector $e_1^T$ replaced by vacuum + first descendant (cell 13: in the 2-dim span, exactly; $\beta/\alpha=392/2335$ at P1, $294/2255$ at $h=1/3$, $5400/337687$ at the third point). For the adjoint $e_1^T(\alpha+\beta E_{12})=\alpha\,e_1^T e^{(\beta/\alpha)E_{12}}$ is itself a rotation ($E_{12}^2$ kills $e_1^T$ there), which is why site 1 was a rotation; for $[3,1,0]$ it is not (`rotRank` 5 = the orbit type of $e_1^TE_{12}$; exact orbit test, 2 solutions, cell 9).
3. **Mixed-point span** ($\theta_1-3h$ shapes plus $\lambda_1\le2$ shapes at $\theta_1-2h$; site-2 span 9-dim, 3549 quadrics, 29 s; cells 8, 10): 18 components, no new point with a site-2 factor other than $\Omega_1$; A and B reappear; the only site-1 rotations come with site 2 $=\Omega_1$ (trivial).
4. **Overlap.** $\langle v_1\otimes v_2|\Psi_B\rangle=c\cdot(1,\tau_\lambda(B))(\theta_1-3h)\,\langle\Omega|\Psi_B\rangle$ automatically (the $T$'s are diagonal). The $Q$-form of this $T$-polynomial (the $[3,1,0]$ analogue of the ladder-1 contour $\det_2$) is **not done**: it needs the $[3,1,0]$ Q-functions/Baxter machinery and the general-rep vacuum determinant.
5. **Torus.** As for the adjoint (`L2-rotation-sov.wb` §11), $y$ is a torus invariant, so point B gives for every site-2 rotation $g$ exactly one partner $e_1^T(\alpha+\beta E_{12})RD\,\rho(\text{torus})$.

**Verdict.** "Two distinct rotations of the vacuum with a factorized overlap" does **not** carry over to $[3,1,0]$ through the $\theta_1-\lambda_1h$ insertions; what carries over is "rotated vacuum on site 2 $\times$ rotated *dressed* vacuum on site 1", with the same $y^*$ pattern. **Open:** closed form of $\beta/\alpha(h,\theta,z,w)$; the Q-form of $c_B\cdot\tau$; the structure of point A's site-1 factor; a meaning for the dressing $\alpha+\beta E_{12}$ (e.g. a rotation with its $E_{12}^2$ component removed, i.e. a sum of two rotations $\alpha e_1^Te^{(\beta/\alpha)E_{12}}-\tfrac{\beta^2}{2\alpha}e_1^TE_{12}^2$).
```

## Deleted Cell 13 — 2026-10-09 16:21:52
**Kind:** code
```wolfram
(* site 1 of point B as a combination of the vacuum and its first descendant with ONE rotation: the orbit solutions have a1 ≠ 0, i.e. site 1 ∝ e1.(α + E12).D.U+ (for the adjoint e1.(α+E12) = α e1.exp(E12/α) is itself a rotation, since E12^2 kills e1 there; for [3,1,0] it is not). Test: is f1 in the 2-dim span {e1.R.D, e1.E12.R.D} (descendant before the rotation), {e1.R.D, e1.R.D.E12} (after), or {e1.R.D, e1.R.E12.D}, for R = Rrot(z,(w1,y*)), D = ρdiag[Dy[w2]], and the (w2,y*)/Dy[w1] and (w1,w2)/no-D variants? *)
inSpan2[f_, u_, v_] := MatrixRank[Normal /@ {u, v, f}] == 2;
spanTests = Flatten[Table[With[{R = Rrot[zN, {pq[[1]], pq[[2]], 1}], D = ρdiag[dd]}, {ToString[pq] <> ", D=" <> dn, "E12 before R" -> inSpan2[f1P3, e1 . R . D, e1 . Ee[1, 2] . R . D], "E12 after D" -> inSpan2[f1P3, e1 . R . D, e1 . R . D . Ee[1, 2]], "E12 between" -> inSpan2[f1P3, e1 . R . D, e1 . R . Ee[1, 2] . D]}], {pq, {{wN[[1]], yB}, {wN[[2]], yB}, {wN[[1]], wN[[2]]}, {wN[[1]], yA}}}, {dd, {Dy[wN[[2]]], Dy[wN[[1]]], {1, 1, 1}}}, {dn, {"Dy[w2]", "Dy[w1]", "1"}}] /. {a_String, b_, c_, d_, e_String} :> Nothing, 1];
spanTests = Flatten[Table[Module[{R = Rrot[zN, {pq[[1]], pq[[2]], 1}], D = ρdiag[dd[[1]]]}, {ToString[pq] <> ", D = " <> dd[[2]], inSpan2[f1P3, e1 . R . D, e1 . Ee[1, 2] . R . D], inSpan2[f1P3, e1 . R . D, e1 . R . D . Ee[1, 2]], inSpan2[f1P3, e1 . R . D, e1 . R . Ee[1, 2] . D]}], {pq, {{wN[[1]], yB}, {wN[[2]], yB}, {wN[[1]], wN[[2]]}, {wN[[1]], yA}}}, {dd, {{Dy[wN[[2]]], "Dy[w2]"}, {Dy[wN[[1]]], "Dy[w1]"}, {{1, 1, 1}, "1"}}}], 1];
checks["site 1 of B in span{e1 R D, e1 E12 R D} etc."] = spanTests;
Grid[Join[{{"(p,q), D", "E12 before R", "E12 after D", "E12 between R and D"}}, spanTests, {{"r(z,(w1,y*)), r(z,w)", rStar, rW, ""}}], Frame -> All, Alignment -> Left]
```

## Deleted Cell 11 — 2026-10-09 16:21:52
**Kind:** code
```wolfram
(* the site-1 factor of the genuine [3,1,0] point (component 3 of cell 6, = #17 of S2): it is (e1.E12).ρ(g). Read off g from the exact solutions of the orbit test (Gauss form U- D U+): a1 is fixed by the seed, a2, a3 and the torus lie in the stabilizer of [e1.E12] (so they stay free); the informative data are the upper-unipotent b's. Compare with the one-site vacuum Gauss parameters r_k(z, (w1, y)) for y = y* and the adjoint's y, and with the resolvent torus Dy[w2]. Also test the direct hypotheses H1: (e1.E12).Rrot(z,(w1,y*)).ρdiag[Dy[w2]], H2: e1.Rrot(z,(w1,y*)).ρdiag[Dy[w2]].E12 (lowering after the rotation), with y* the site-2 invariant and w2 free. *)
yB = -63/43; yA = -63/8;
solOrbB = Select[Solve[Thread[Normal[Normal[e1 . Ee[1, 2]] . gGauss] == sv f1P3], {a1, a2, a3, d1v, d2v, b1, b2, b3, sv}], (sv /. #) =!= 0 &];
checks["orbit solutions for site 1 of point B (Gauss parameters)"] = solOrbB;
rStar = rOm[zN, {wN[[1]], yB, 1}]; rAdj = rOm[zN, {wN[[1]], yA, 1}]; rW = rOm[zN, wN];
sameLine[u_, v_] := MatrixRank[{Normal[u], Normal[v]}] == 1;
hyp = {"H1 (e1.E12).Rrot(w1,y*).ρ(Dy[w2])" -> sameLine[e1 . Ee[1, 2] . Rrot[zN, {wN[[1]], yB, 1}] . ρdiag[Dy[wN[[2]]]], f1P3], "H2 e1.Rrot(w1,y*).ρ(Dy[w2]).E12" -> sameLine[e1 . Rrot[zN, {wN[[1]], yB, 1}] . ρdiag[Dy[wN[[2]]]] . Ee[1, 2], f1P3], "H3 (e1.E12).Rrot(w2,y*).ρ(Dy[w1])" -> sameLine[e1 . Ee[1, 2] . Rrot[zN, {wN[[2]], yB, 1}] . ρdiag[Dy[wN[[1]]]], f1P3], "H4 e1.Rrot(w2,y*).ρ(Dy[w1]).E12" -> sameLine[e1 . Rrot[zN, {wN[[2]], yB, 1}] . ρdiag[Dy[wN[[1]]]] . Ee[1, 2], f1P3], "H5 (e1.E12).Rrot(w1,w2)" -> sameLine[e1 . Ee[1, 2] . Rrot[zN, wN], f1P3], "H6 e1.Rrot(w1,w2).E12" -> sameLine[e1 . Rrot[zN, wN] . Ee[1, 2], f1P3]};
(* H7: (e1.E12).Rrot(z,(p,q)).ρdiag[d] with (p, q, d) free: exact Solve (this is the 'torus-rotated vacuum descendant' form) *)
solH7 = Select[Solve[Thread[Normal[e1 . Ee[1, 2] . Rrot[zN, {pv, qv, 1}] . ρdiag[{d1v, d2v, 1}]] == sv f1P3], {pv, qv, d1v, d2v, sv}], (sv /. #) =!= 0 &];
checks["site 1 of point B = (e1.E12).Ω-rotation(p,q).torus: solutions"] = solH7;
Grid[{{"orbit solutions (Gauss form), site 1 of point B", Column[solOrbB]}, {"r_k(z,(w1,y*)), r_k(z,(w1,y_adj)), r_k(z,w)", {rStar, rAdj, rW}}, {"Dy[w2], Dy[w1]", {Dy[wN[[2]]], Dy[wN[[1]]]}}, {"hypotheses", Column[hyp]}, {"H7 (e1.E12).Rrot(p,q).ρ(d) = site 1: solutions (p, q, d)", Column[solH7]}}, Frame -> All, Alignment -> Left]
```

## Deleted Cell 10 — 2026-10-09 16:21:52
**Kind:** code
```wolfram
(* S2 in detail (shapes λ1 ≤ 3 at θ1-3h plus λ1 ≤ 2 at θ1-2h): keep the solutions; for every component report {free params (scale and the vanishing T(3,3)@-3 are always free), rotRank site 1, rotRank site 2, rank{site1,Ω1}, rank{site2,Ω1}, support}. For product points with a ROTATION on site 1 distinct from Ω1: torus invariants y of both factors and whether site 1 = site 2 · ρ(torus) (same invariant). *)
namesS2 = Join[{"1"}, ("T" <> ToString[#] <> "@-3" &) /@ shapes, ("T" <> ToString[#] <> "@-2" &) /@ shapes2];
βS2 = Module[{b = {Ω1}}, Do[If[MatrixRank[Append[b, row]] > Length[b], AppendTo[b, row/row[[Ordering[Abs[row], -1][[1]]]]]], {row, Flatten[matsS2, 1]}]; b];
AcS2 = Table[Transpose[Table[LinearSolve[Transpose[βS2], M[[i]]], {i, d1}]], {M, matsS2}];
cvS2 = Array[cc, Length[matsS2]]; eqsS2 = DeleteCases[Union[Flatten[Minors[Sum[cvS2[[k]] AcS2[[k]], {k, Length[matsS2]}], 2]]], 0];
{tSolS2, solS2} = AbsoluteTiming[Solve[eqsS2 == 0, cvS2]];
SeedRandom[312];
facS2[s_] := Module[{c1 = cvS2 /. s, free, M, rowi, colj}, free = Variables[c1]; c1 = c1 /. Thread[free -> RandomInteger[{2, 9}, Length[free]]]; M = Sum[c1[[k]] matsS2[[k]], {k, Length[matsS2]}];
   If[Union[Flatten[M]] === {0}, {Length[free], "zero", c1, 0, 0}, rowi = First[FirstPosition[Total[Abs[M], {2}], _?(# != 0 &)]]; colj = First[FirstPosition[Total[Abs[M], {1}], _?(# != 0 &)]]; {Length[free], "ok", c1, M[[All, colj]]/M[[rowi, colj]], M[[rowi, All]]/M[[rowi, colj]]}]];
facsS2 = facS2 /@ solS2;
supportOf[c1_] := Pick[namesS2, Map[# =!= 0 &, c1]];
classS2 = Table[With[{f = facsS2[[k]]}, If[f[[2]] === "zero", {k, f[[1]], "zero", "", "", "", ""}, {k, f[[1]], rotRank[f[[4]]], rotRank[f[[5]]], MatrixRank[{f[[4]], Ω1}], MatrixRank[{f[[5]], Ω1}], supportOf[f[[3]]]}]], {k, Length[solS2]}];
(* product points with a rotation on site 1 that is not Ω1 *)
rotPts = Select[classS2, #[[3]] === 4 && #[[5]] == 2 &];
yInv[g_] := Module[{s = solveVac[g]}, If[s === {}, "none", yv /. s]];
rotDetail = Table[With[{f = facsS2[[p[[1]]]]}, {p[[1]], "y(site1)" -> yInv[f[[4]]], "y(site2)" -> yInv[f[[5]]], "site2 = Ω1?" -> (MatrixRank[{f[[5]], Ω1}] == 1)}], {p, rotPts}];
checks["S2 classification"] = classS2; checks["S2 rotation-on-site-1 points"] = rotDetail;
Grid[Join[{{"k", "free", "rotRank s1", "rotRank s2", "rank{s1,Ω1}", "rank{s2,Ω1}", "support"}}, classS2, {{"Solve time (s)", tSolS2, "", "", "", "", ""}}, {{"rotation-on-site-1 points: torus invariants", Column[rotDetail], "", "", "", "", ""}}], Frame -> All, Alignment -> Left]
```

## Deleted Cell 8 — 2026-10-09 16:21:52
**Kind:** code
```wolfram
(* enlarged operator spans: does a rotation on site 1 appear? S1: all shapes with λ1 ≤ 4 at θ1-3h; S2: shapes λ1 ≤ 3 at θ1-3h plus shapes λ1 ≤ 2 at θ1-2h (the adjoint's point). Same machinery as cell 6, packaged; Solve is time-limited. *)
locusFor[mats_, names_] := Module[{β, rb, Ac, cv, aC, eqs, sol, tm, cls},
  β = Module[{b = {Ω1}}, Do[If[MatrixRank[Append[b, row]] > Length[b], AppendTo[b, row/row[[Ordering[Abs[row], -1][[1]]]]]], {row, Flatten[mats, 1]}]; b]; rb = Length[β];
  Ac = Table[Transpose[Table[LinearSolve[Transpose[β], M[[i]]], {i, d1}]], {M, mats}];
  cv = Array[cc, Length[mats]]; aC = Sum[cv[[k]] Ac[[k]], {k, Length[mats]}];
  eqs = DeleteCases[Union[Flatten[Minors[aC, 2]]], 0];
  {tm, sol} = AbsoluteTiming[TimeConstrained[Solve[eqs == 0, cv], 200, $Failed]];
  If[sol === $Failed, {"site-2 span dim" -> rb, "# quadrics" -> Length[eqs], "Solve" -> "timed out (200 s)"},
   cls = Table[Module[{c1 = cv /. s, free, M, rowi, colj, f1, f2}, free = Variables[c1]; c1 = c1 /. Thread[free -> RandomInteger[{2, 9}, Length[free]]]; M = Sum[c1[[k]] mats[[k]], {k, Length[mats]}];
      If[Union[Flatten[M]] === {0}, {Length[free], "zero"}, rowi = First[FirstPosition[Total[Abs[M], {2}], _?(# != 0 &)]]; colj = First[FirstPosition[Total[Abs[M], {1}], _?(# != 0 &)]]; f1 = M[[All, colj]]; f2 = M[[rowi, All]];
       {Length[free], rotRank[f1], rotRank[f2], MatrixRank[{f2, Ω1}], Pick[names, Thread[(cv /. s) =!= 0] /. Thread[free -> 1]]}]], {s, sol}];
   {"site-2 span dim" -> rb, "# quadrics" -> Length[eqs], "Solve time (s)" -> N[tm, 3], "# components" -> Length[sol], "components {free, rotRank site1, rotRank site2, rank{site2,Ω1}, support}" -> cls}]];
shapes4 = Join[shapes, {{4, 0}, {4, 1}, {4, 2}, {4, 3}, {4, 4}}];
{tS1, matsS1} = AbsoluteTiming[Join[{Partition[Ω2, d1]}, Table[Partition[Normal[Ω2 . TL2[λ, v0, Zd]], d1], {λ, shapes4}]]];
SeedRandom[311]; resS1 = locusFor[matsS1, Join[{"1"}, ("T" <> ToString[#] <> "@-3" &) /@ shapes4]];
shapes2 = {{1, 0}, {1, 1}, {2, 0}, {2, 1}};
{tS2, matsS2} = AbsoluteTiming[Join[{Partition[Ω2, d1]}, Table[Partition[Normal[Ω2 . TL2[λ, v0, Zd]], d1], {λ, shapes}], Table[Partition[Normal[Ω2 . TL2[λ, v0 + hN, Zd]], d1], {λ, shapes2}]]];
resS2 = locusFor[matsS2, Join[{"1"}, ("T" <> ToString[#] <> "@-3" &) /@ shapes, ("T" <> ToString[#] <> "@-2" &) /@ shapes2]];
Grid[{{"S1: λ1 ≤ 4 at θ1-3h (build time s)", tS1}, {"S1 result", Column[resS1]}, {"S2: λ1 ≤ 3 at θ1-3h + λ1 ≤ 2 at θ1-2h (build time s)", tS2}, {"S2 result", Column[resS2]}}, Frame -> All, Alignment -> Left]
```

## Deleted Cell 26 — 2026-10-09 16:35:23
**Kind:** code
```wolfram
Grid[{{"test", 1 + 1, $KernelID}}]
```

## Deleted Cell 25 — 2026-10-09 16:35:23
**Kind:** code
```wolfram
(* torus-completed HW count for point A alone and A union B: the 12 random (w2, D) pairs of cell 29 (fresh pipeline runs with f1, f2 recorded), rows <f1_A rho(D) x f2_A rho(D)|HW>. *)
{u14tRnd, u14rndRuns} = AbsoluteTiming[Association[Table[c -> u14pipe[zN, {u13w1, c, u13w3}, hN, thN], {c, u13rndW2}]]];
u14rndA = Association[Table[c -> Select["nontrivial points" /. u14rndRuns[c], ListQ["y(site2)" /. #] && ("y(site2)" /. #)[[1]] == u14yA[c] &], {c, u13rndW2}]];
checks["u14 torus: point A found for all 12 random w2"] = Values[Length /@ u14rndA];
u14torusRowA[i_] := With[{c = u13rndW2[[i]], D = u13rndD[[i]]}, With[{f1 = Normal["f1" /. u14rndA[c][[1]]], f2 = Normal["f2" /. u14rndA[c][[1]]], rD = ρdiag[D]}, kron[f1 . rD, f2 . rD] . u13HW]];
u14MTA = Table[u14torusRowA[i], {i, 12}];
u14MTU = Join[u14MTA, u13MT["a"]];
checks["u14 HW rank torus-completed: A alone, A+B"] = {MatrixRank[u14MTA], MatrixRank[u14MTU]};
u14perIrr[M_] := Table[{ν, Length[u13colsOf[ν]], MatrixRank[M[[All, u13colsOf[ν]]]]}, {ν, u13irreps}];
checks["u14 per-irrep torus-completed rank A"] = u14perIrr[u14MTA];
checks["u14 per-irrep torus-completed rank A+B"] = u14perIrr[u14MTU];
Grid[{{"HW rank torus-completed: A, A+B", checks["u14 HW rank torus-completed: A alone, A+B"]}, {"per irrep {irrep, mult, rank} A", Column[checks["u14 per-irrep torus-completed rank A"]]}, {"per irrep A+B", Column[checks["u14 per-irrep torus-completed rank A+B"]]}, {"time (s)", u14tRnd}}, Frame -> All, Alignment -> Left]
```

## Deleted Cell 24 — 2026-10-09 16:35:23
**Kind:** code
```wolfram
(* point A undo: n_A = Gauss element at y_A (rOm), f2 = sA e1.rho(n_A), phi_A = f1_A.rho(n_A)^-1, <f1_A x f2|HW> = sA <phi_A x e1|HW>. Per-weight ranks of phi_A over the six w2, and the UNION with the point-B phi's (variant a, alpha = 1; cell 28). *)
u14rhoA[w2_] := u14rhoNr[u14rA0[w2]];
u14sA = Association[Table[w2 -> u14Af2[w2][[1]]/Normal[e1 . u14rhoA[w2]][[1]], {w2, u13w2s}]];
checks["u14 f2_A == sA e1.rho(n_A) exactly"] = Table[u14Af2[w2] == u14sA[w2] Normal[e1 . u14rhoA[w2]], {w2, u13w2s}];
u14phiA[w2_] := u14Af1[w2] . Inverse[u14rhoA[w2]];
checks["u14 (rho(n_A) x rho(n_A)).HW == HW (two w2)"] = Table[With[{R = u14rhoA[w2]}, KroneckerProduct[R, R] . u13HW == u13HW], {w2, {7/3, -1/2}}];
checks["u14 <f1_A x f2_A|HW> == sA <phi_A x e1|HW> (all w2, all 10 HW)"] = Table[kron[u14Af1[w2], u14Af2[w2]] . u13HW == u14sA[w2] kron[u14phiA[w2], e1] . u13HW, {w2, u13w2s}];
u14PhiA = Table[u14phiA[w2], {w2, u13w2s}];
u14PhiB = u13Phi["a"];
u14PhiU = Join[u14PhiA, u14PhiB];
u14rk[Phi_] := Table[With[{M = Phi[[All, u13idx[μ]]]}, {μ, Length[u13idx[μ]], MatrixRank[M], Flatten[Position[Transpose[M], _?(Union[#] === {0} &)]] // Map[u13idx[μ][[#]] &]}], {μ, u13wtList}];
u14rkA = u14rk[u14PhiA]; u14rkB = u14rk[u14PhiB]; u14rkU = u14rk[u14PhiU];
checks["u14 per-weight ranks A, B, A+B"] = {u14rkA[[All, 3]], u14rkB[[All, 3]], u14rkU[[All, 3]]};
checks["u14 total span rank of phi (15): A, B, A+B"] = {MatrixRank[u14PhiA], MatrixRank[u14PhiB], MatrixRank[u14PhiU]};
(* fixed-torus HW counts (six w2 each) *)
u14MfixA = Table[kron[u14PhiA[[i]], e1] . u13HW, {i, 6}];
u14MfixU = Join[u14MfixA, u13MF["a"]];
checks["u14 HW rank fixed torus: A, A+B"] = {MatrixRank[u14MfixA], MatrixRank[u14MfixU]};
Grid[Join[{{"weight μ", "dim V_μ", "rank A", "rank B", "rank A+B", "vanishing idx A"}}, Table[{u14rkA[[k, 1]], u14rkA[[k, 2]], u14rkA[[k, 3]], u14rkB[[k, 3]], u14rkU[[k, 3]], u14rkA[[k, 4]]}, {k, Length[u13wtList]}], {{"total rank (15)", "", MatrixRank[u14PhiA], MatrixRank[u14PhiB], MatrixRank[u14PhiU], ""}}], Frame -> All, Alignment -> Left]
```

## Deleted Cell 23 — 2026-10-09 16:35:23
**Kind:** code
```wolfram
(* v2_A versus the vacuum rotation e1.rho(n_A), n_A the N+ element of rOm[zN,{w1,y_A,w3}]; if only a torus-rotated vacuum, extract the N+ element n_A' by an exact Solve of e1.rho(exp(r1 E12) exp(r2 E13) exp(r3 E23)) == s f2. *)
u14rhoNr[r_] := MatrixExp[r[[1]] Normal[Ee[1, 2]]] . MatrixExp[r[[2]] Normal[Ee[1, 3]]] . MatrixExp[r[[3]] Normal[Ee[2, 3]]];
u14rA0[w2_] := rOm[zN, {u13w1, u14yA[w2], u13w3}];
checks["u14 v2_A parallel to e1.rho(n_A) with n_A = rOm[z,{w1,y_A,w3}] ?"] = Table[MatrixRank[{u14Af2[w2], Normal[e1 . u14rhoNr[u14rA0[w2]]]}] == 1, {w2, u13w2s}];
u14solA[w2_] := Module[{R = u14rhoNr[{q1, q2, q3}], sol}, sol = Solve[Thread[Normal[e1 . R] == ss u14Af2[w2]], {q1, q2, q3, ss}]; sol];
u14rAs = Association[Table[w2 -> u14solA[w2], {w2, u13w2s}]];
checks["u14 v2_A = s e1.rho(n_A') solvable, # solutions"] = Length /@ Values[u14rAs];
Grid[{{"v2_A ∥ e1.rho(n_A) (rOm at y_A)?", checks["u14 v2_A parallel to e1.rho(n_A) with n_A = rOm[z,{w1,y_A,w3}] ?"]}, {"# solutions of the exact Solve", checks["u14 v2_A = s e1.rho(n_A') solvable, # solutions"]}, {"r(w2=7/3) from rOm at y_A", u14rA0[7/3]}, {"Solve r(w2=7/3)", u14rAs[7/3]}}, Frame -> All, Alignment -> Left]
```

## Deleted Cell 22 — 2026-10-09 16:35:23
**Kind:** code
```wolfram
(* Task 2, point A. Re-run the pipeline with BOTH site covectors recorded (u14pipe = cell 12 source with f1 and f2 added to each point record), six w2 (same as point B), then pick the point whose site-2 invariant is y_A = w2 (th12-h)/(th12-2h). *)
u14def = StringReplace[u13def, {"u13pipe" -> "u14pipe", "{\"f1\" -> f1, " -> "{\"f1\" -> f1, \"f2\" -> f2, "}];
ToExpression[u14def];
u14yA[w2_] := w2 (θ12 - hN)/(θ12 - 2 hN);
{u14tA, u14runs} = AbsoluteTiming[Association[Table[w2 -> u14pipe[zN, {u13w1, w2, u13w3}, hN, thN], {w2, u13w2s}]]];
u14pickA[w2_] := Select["nontrivial points" /. u14runs[w2], ListQ["y(site2)" /. #] && ("y(site2)" /. #)[[1]] == u14yA[w2] &];
u14A = Association[Table[w2 -> u14pickA[w2], {w2, u13w2s}]];
checks["u14 point A found for each w2 (exactly one)"] = Values[Length /@ u14A];
u14Af1 = Association[Table[w2 -> Normal["f1" /. u14A[w2][[1]]], {w2, u13w2s}]];
u14Af2 = Association[Table[w2 -> Normal["f2" /. u14A[w2][[1]]], {w2, u13w2s}]];
checks["u14 point A rotRank (site1, site2) per w2"] = Table["rotRank" /. u14A[w2][[1]], {w2, u13w2s}];
(* consistency: the same runs give point B with the same f1 as before (parallel) *)
checks["u14 point B f1 of the new runs parallel to u13Bf1"] = Table[With[{pb = Select["nontrivial points" /. u14runs[w2], ListQ["y(site2)" /. #] && ("y(site2)" /. #)[[1]] == u13ystar[w2] &]}, MatrixRank[{Normal["f1" /. pb[[1]]], Normal[u13Bf1[w2]]}] == 1], {w2, u13w2s}];
Grid[Join[{{"w2", "y_A", "rotRank (site1, site2)"}}, Table[{w2, u14yA[w2], "rotRank" /. u14A[w2][[1]]}, {w2, u13w2s}], {{"time (s)", u14tA, ""}, {"B f1 consistent", checks["u14 point B f1 of the new runs parallel to u13Bf1"], ""}}], Frame -> All, Alignment -> Left]
```

## Deleted Cell 21 — 2026-10-09 16:42:42
**Kind:** code
```wolfram
(* LW frame: per-weight ranks of phi' over the six w2, identically vanishing components, per-irrep determination (rank at the site-1 weight where the LW slice of nu sits: nu_low - (0,1,3)), torus-completed LW count, and the HW-frame determination (site-1 weight nu - (3,1,0), variant (a) of cell 28) for the combined verdict. *)
u14Phi = Table[u14phi[w2], {w2, u13w2s}];
u14rankTab = Table[With[{M = u14Phi[[All, u13idx[μ]]]}, {μ, Length[u13idx[μ]], MatrixRank[M], Flatten[Position[Transpose[M], _?(Union[#] === {0} &)]] // Map[u13idx[μ][[#]] &]}], {μ, u13wtList}];
checks["u14 LW-frame per-weight ranks {mu, dim, rank}"] = u14rankTab[[All, 1 ;; 3]];
checks["u14 LW-frame total rank of phi' span (of 15)"] = MatrixRank[u14Phi];
u14rkOfLW[μ_] := Cases[u14rankTab, {μ, dm_, rk_, _} :> {dm, rk}][[1]];
u14rkOfHW[μ_] := Cases[u13ranks["a"], {μ, dm_, rk_, _} :> {dm, rk}][[1]];
u14irr = Union[u13HWwt];
u14verd = Table[Module[{lw = Reverse[ν], μL, μH, rl, rh}, μL = lw - {0, 1, 3}; μH = ν - Λ; rl = u14rkOfLW[μL]; rh = u14rkOfHW[μH]; <|"irrep" -> ν, "mult" -> Count[u13HWwt, ν], "HW site-1 wt" -> μH, "HW dim,rank" -> rh, "HW determined" -> (rh[[1]] == rh[[2]]), "LW site-1 wt" -> μL, "LW dim,rank" -> rl, "LW determined" -> (rl[[1]] == rl[[2]]), "either" -> (rh[[1]] == rh[[2]] || rl[[1]] == rl[[2]])|>], {ν, u14irr}];
(* torus-completed LW count with the 12 random (w2, D) pairs of cell 29: rows <v1 rho(D) x v2 rho(D)|LW> *)
u14torusRow[c_, D_] := Module[{v2 = u13v2[c], v1 = u13v1[c, u13rndBeta[c], c], rD = ρdiag[D]}, kron[v1 . rD, v2 . rD] . u14LW];
u14MT = Table[u14torusRow[u13rndW2[[i]], u13rndD[[i]]], {i, 12}];
u14colsOf[ν_] := Flatten[Position[u14LWwt, Reverse[ν]]];
u14lwCount = {MatrixRank[u14MT], Table[{ν, Length[u14colsOf[ν]], MatrixRank[u14MT[[All, u14colsOf[ν]]]]}, {ν, u14irr}]};
checks["u14 LW count torus-completed (expect 10)"] = u14lwCount[[1]];
u14MFix = Table[kron[u14Phi[[i]], u14elow] . u14LW, {i, 6}];
checks["u14 LW count fixed torus six w2"] = MatrixRank[u14MFix];
Grid[Join[{{"weight μ", "dim V_μ", "rank φ' (LW frame)", "identically vanishing indices", "rank φ (HW frame, point B)"}}, Table[{u14rankTab[[k, 1]], u14rankTab[[k, 2]], u14rankTab[[k, 3]], u14rankTab[[k, 4]], u13ranks["a"][[k, 3]]}, {k, Length[u13wtList]}], {{"total rank (15)", "", MatrixRank[u14Phi], "", MatrixRank[u13Phi["a"]]}}], Frame -> All, Alignment -> Left]
```

