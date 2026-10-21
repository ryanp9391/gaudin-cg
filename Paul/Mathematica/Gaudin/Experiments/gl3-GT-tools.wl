(* ::Package:: *)

(* gl3-GT-tools.wl -- the gl(3) Gelfand-Tsetlin machinery of Minimal_complete_working.nb (repo root), extracted verbatim
   (cells: GT patterns/generators/vertex operator; computeHWS; the vertex operator V) so that it can be loaded without the .nb:
       Get["<repo>/Paul/Mathematica/Gaudin/Experiments/gl3-GT-tools.wl"]
   plus computeV with two changes vs the notebook: Subscript[h,i,j,k] -> hh[i,j,k] (a Block-local symbol; the notebook's subscripts
   escape the Block and make the second call reuse the previous lambda, mu) and the final NSolve replaced by pencilSolve
   (linear-pencil solver of Minimal_complete_working-fast.wl, which must be loaded AFTER this file). computeV[lam,mu][nu] returns
   {V, q1, q2, q3} per solution; the sites sit at u = 1 (site 1, rep lam) and u = 0 (site 2, rep mu), theta12 = 1.
   Conventions: ps[{lam}] = GT patterns, Reverse@Sort order, index 1 = highest weight; tp[pat1, pat2] = product state;
   ee[i,j][n] = E_ij on site n (rational, non-unitary GT normalisation); V = the vertex operator (leading-order Hamiltonian). *)

(* ---- GT patterns, local/global Cartan, raising/lowering, generators ee[i,j][n] (notebook cell 'GT patterns, generators and vertex operator') ---- *)
ClearAll[ps, H, J0, JJp, JJP, ee];
ps[rep_] := ps[rep] = Block[{res = (Prepend[rep, #1] & ) /@ Tuples[Table[Range[rep[[1,i + 1]], rep[[1,i]]], {i, 1, Length[rep[[1]]] - 1}]]}, If[Length[res[[1,1]]] == 1, Reverse[Sort[res]], Reverse[Sort[Flatten[ps /@ res, 1]]]]];
J0[k_][pat_] := J0[k][pat] = Sum[pat[[k,i]], {i, 1, k}] - Sum[pat[[k - 1,i]], {i, 1, k - 1}];
H[k_][pats_] := H[k][pats] = J0[k][pats[[1]]] + J0[k][pats[[2]]];
JJp[k_][n_][(a_) + (b_)] := JJp[k][n][a] + JJp[k][n][b];
JJp[k_][n_][(a_)*(b_)] := a*JJp[k][n][b] /; FreeQ[a, tp];
JJp[k_][n_][0] := JJp[k][n][0] = 0;
JJp[k_][1][tp[pat1_, pat2_]] := JJp[k][1][tp[pat1, pat2]] = -Sum[(Product[pat1[[k,i]] - i - pat1[[k + 1,j]] + j, {j, Range[k + 1]}]/Product[pat1[[k,i]] - i - pat1[[k,j]] + j, {j, Complement[Range[k], {i}]}])*tp[ReplacePart[pat1, {k, i} -> pat1[[k,i]] + 1], pat2], {i, 1, k}];
JJp[k_][2][tp[pat1_, pat2_]] := JJp[k][2][tp[pat1, pat2]] = -Sum[(Product[pat2[[k,i]] - i - pat2[[k + 1,j]] + j, {j, Range[k + 1]}]/Product[pat2[[k,i]] - i - pat2[[k,j]] + j, {j, Complement[Range[k], {i}]}])*tp[pat1, ReplacePart[pat2, {k, i} -> pat2[[k,i]] + 1]], {i, 1, k}];
JJP[k_][a_] := JJP[k][a] = JJp[k][1][a] + JJp[k][2][a];
JJm[k_][n_][(a_) + (b_)] := JJm[k][n][a] + JJm[k][n][b];
JJm[k_][n_][(a_)*(b_)] := a*JJm[k][n][b] /; FreeQ[a, tp];
JJm[k_][n_][0] := JJm[k][n][0] = 0;
JJm[k_][1][tp[pat1_, pat2_]] := JJm[k][1][tp[pat1, pat2]] = Sum[(Product[pat1[[k,i]] - i - pat1[[k - 1,j]] + j, {j, Range[k - 1]}]/Product[pat1[[k,i]] - i - pat1[[k,j]] + j, {j, Complement[Range[k], {i}]}])*tp[ReplacePart[pat1, {k, i} -> pat1[[k,i]] - 1], pat2], {i, 1, k}];
JJm[k_][2][tp[pat1_, pat2_]] := JJm[k][2][tp[pat1, pat2]] = Sum[(Product[pat2[[k,i]] - i - pat2[[k - 1,j]] + j, {j, Range[k - 1]}]/Product[pat2[[k,i]] - i - pat2[[k,j]] + j, {j, Complement[Range[k], {i}]}])*tp[pat1, ReplacePart[pat2, {k, i} -> pat2[[k,i]] - 1]], {i, 1, k}];
ClearAll[testpat];
testpat[pat_] := testpat[pat] = pat[[1,1]] <= pat[[2,1]] && pat[[2,2]] <= pat[[2,1]] && pat[[2,2]] <= pat[[1,1]] && pat[[2,2]] <= pat[[3,2]] && pat[[2,2]] >= pat[[3,3]] && pat[[2,1]] <= pat[[3,1]] && pat[[2,1]] >= pat[[3,2]]
removebadpats = {tp[pat1_, pat2_] :> If[testpat[pat1] && testpat[pat2], tp[pat1, pat2], 0]};
ee[i_, j_][n_][(a_) + (b_)] := ee[i, j][n][a] + ee[i, j][n][b];
ee[i_, j_][n_][(a_)*(b_)] := a*ee[i, j][n][b] /; FreeQ[a, tp];
ee[1, 2][n_][a_] := ee[1, 2][n][a] = JJp[1][n][a] /. removebadpats;
ee[2, 3][n_][a_] := ee[2, 3][n][a] = JJp[2][n][a] /. removebadpats;
ee[1, 3][n_][a_] := ee[1, 3][n][a] = ee[1, 2][n][ee[2, 3][n][a]] - ee[2, 3][n][ee[1, 2][n][a]];
ee[2, 1][n_][a_] := ee[2, 1][n][a] = JJm[1][n][a] /. removebadpats;
ee[3, 2][n_][a_] := ee[3, 2][n][a] = JJm[2][n][a] /. removebadpats;
ee[3, 1][n_][a_] := ee[3, 1][n][a] = ee[3, 2][n][ee[2, 1][n][a]] - ee[2, 1][n][ee[3, 2][n][a]];
ee[i_, i_][1][tp[pat1_, pat2_]] := ee[i, i][1][pat1, pat2] = J0[i][pat1]*tp[pat1, pat2];
ee[i_, i_][2][tp[pat1_, pat2_]] := ee[i, i][2][pat1, pat2] = J0[i][pat2]*tp[pat1, pat2];

(* ---- computeHWS: the HW space of weight nu in lam x mu as a sum of cc[n] tp[GT1,GT2] (symbolic; slow, reference only) ---- *)
ClearAll[computeHWS];
computeHWS[λ_, μ_][ν_] := computeHWS[λ, μ][ν] = Block[{states, weights, weightstates, vec, cc, eqns, sp, wp, tab, rules}, states = Tuples[{ps[{λ}], ps[{μ}]}]; weights = findweight[λ, μ][ν]; weightstates = Select[states, H[1][#1] == weights[[1]] && H[2][#1] == weights[[2]] && H[3][#1] == weights[[3]] & ]; vec = Table[cc[n], {n, 1, Length[weightstates]}]*tp @@@ weightstates; sp = ExpandAll[{Plus @@ (JJP[1][#1] & ) /@ vec, Plus @@ (JJP[2][#1] & ) /@ vec} /. removebadpats]; eqns = Union[Flatten[List @@@ (Collect[#1, tp[__]] & ) /@ sp /. tp[a__] -> 1]]; (Collect[#1, tp[__], Expand] & )[Plus @@ vec /. Last[Solve[eqns == 0]]]];

(* ---- the vertex operator V (36 terms; = 1/2 Sum_{s,t in S3} sgn(s) sgn(t) e1_{s1 t1} e1_{s2 t2} e2_{s3 t3}) ---- *)
ClearAll[V];
V[x_] := V[x] = (1/2)*ee[1, 1][1][ee[2, 2][1][ee[3, 3][2][x]]] - (1/2)*ee[1, 1][1][ee[2, 3][1][ee[3, 2][2][x]]] - (1/2)*ee[1, 1][1][ee[3, 2][1][ee[2, 3][2][x]]] + (1/2)*ee[1, 1][1][ee[3, 3][1][ee[2, 2][2][x]]] - (1/2)*ee[1, 2][1][ee[2, 1][1][ee[3, 3][2][x]]] + (1/2)*ee[1, 2][1][ee[2, 3][1][ee[3, 1][2][x]]] + (1/2)*ee[1, 2][1][ee[3, 1][1][ee[2, 3][2][x]]] - (1/2)*ee[1, 2][1][ee[3, 3][1][ee[2, 1][2][x]]] + (1/2)*ee[1, 3][1][ee[2, 1][1][ee[3, 2][2][x]]] - (1/2)*ee[1, 3][1][ee[2, 2][1][ee[3, 1][2][x]]] - (1/2)*ee[1, 3][1][ee[3, 1][1][ee[2, 2][2][x]]] + (1/2)*ee[1, 3][1][ee[3, 2][1][ee[2, 1][2][x]]] - (1/2)*ee[2, 1][1][ee[1, 2][1][ee[3, 3][2][x]]] + (1/2)*ee[2, 1][1][ee[1, 3][1][ee[3, 2][2][x]]] + (1/2)*ee[2, 1][1][ee[3, 2][1][ee[1, 3][2][x]]] - (1/2)*ee[2, 1][1][ee[3, 3][1][ee[1, 2][2][x]]] + (1/2)*ee[2, 2][1][ee[1, 1][1][ee[3, 3][2][x]]] - (1/2)*ee[2, 2][1][ee[1, 3][1][ee[3, 1][2][x]]] - (1/2)*ee[2, 2][1][ee[3, 1][1][ee[1, 3][2][x]]] + (1/2)*ee[2, 2][1][ee[3, 3][1][ee[1, 1][2][x]]] - (1/2)*ee[2, 3][1][ee[1, 1][1][ee[3, 2][2][x]]] + (1/2)*ee[2, 3][1][ee[1, 2][1][ee[3, 1][2][x]]] + (1/2)*ee[2, 3][1][ee[3, 1][1][ee[1, 2][2][x]]] - (1/2)*ee[2, 3][1][ee[3, 2][1][ee[1, 1][2][x]]] + (1/2)*ee[3, 1][1][ee[1, 2][1][ee[2, 3][2][x]]] - (1/2)*ee[3, 1][1][ee[1, 3][1][ee[2, 2][2][x]]] - (1/2)*ee[3, 1][1][ee[2, 2][1][ee[1, 3][2][x]]] + (1/2)*ee[3, 1][1][ee[2, 3][1][ee[1, 2][2][x]]] - (1/2)*ee[3, 2][1][ee[1, 1][1][ee[2, 3][2][x]]] + (1/2)*ee[3, 2][1][ee[1, 3][1][ee[2, 1][2][x]]] + (1/2)*ee[3, 2][1][ee[2, 1][1][ee[1, 3][2][x]]] - (1/2)*ee[3, 2][1][ee[2, 3][1][ee[1, 1][2][x]]] + (1/2)*ee[3, 3][1][ee[1, 1][1][ee[2, 2][2][x]]] - (1/2)*ee[3, 3][1][ee[1, 2][1][ee[2, 1][2][x]]] - (1/2)*ee[3, 3][1][ee[2, 1][1][ee[1, 2][2][x]]] + (1/2)*ee[3, 3][1][ee[2, 2][1][ee[1, 1][2][x]]];

(* ---- computeV: the untwisted (leading-order) Q-system; needs pencilSolve from Minimal_complete_working-fast.wl ---- *)
ClearAll[findweight, computeV];
findweight[λ_, μ_][ν_] := findweight[λ, μ][ν] = {(2*ν[[1]])/3 - ν[[2]]/3 - ν[[3]]/3, -(ν[[1]]/3) + (2*ν[[2]])/3 - ν[[3]]/3, -(ν[[1]]/3) - ν[[2]]/3 + (2*ν[[3]])/3} + (1/3)*Plus @@ (λ + μ);
computeV[λ_, μ_][ν_] := computeV[λ, μ][ν] = Block[{h, ℰ, M, l, Q, q, Bax, V, diffs, min, minpos, lowest, hh}, h[1][z_] = Plus @@ λ/(z - 1) + Plus @@ μ/z; h[2][z_] = Sum[(hh[2, r, s]/((z - 1)^r*z^s))*Boole[r + s == 2], {r, 0, 2}, {s, 0, 2}]; h[3][z_] = Sum[(hh[3, r, s]/((z - 1)^r*z^s))*Boole[r + s == 3], {r, 0, 3}, {s, 0, 3}]; hh[2, 2, 0] = -λ[[1]] + λ[[1]]*λ[[2]] + λ[[3]] + λ[[1]]*λ[[3]] + λ[[2]]*λ[[3]]; hh[2, 0, 2] = -μ[[1]] + μ[[1]]*μ[[2]] + μ[[3]] + μ[[1]]*μ[[3]] + μ[[2]]*μ[[3]]; hh[3, 3, 0] = λ[[1]]/3 - (2*λ[[2]])/3 - λ[[1]]*λ[[2]] + λ[[3]]/3 + λ[[2]]*λ[[3]] + λ[[1]]*λ[[2]]*λ[[3]]; hh[3, 0, 3] = μ[[1]]/3 - (2*μ[[2]])/3 - μ[[1]]*μ[[2]] + μ[[3]]/3 + μ[[2]]*μ[[3]] + μ[[1]]*μ[[2]]*μ[[3]]; hh[3, 2, 1] = V; hh[3, 1, 2] = -V - (1 + ℰ[1])*ℰ[2]^2 + λ[[2]] + λ[[1]]*λ[[2]] - λ[[2]]*λ[[3]] - λ[[1]]*λ[[2]]*λ[[3]] + μ[[2]] + μ[[1]]*μ[[2]] - (1 + ℰ[1])*ℰ[2]*(1 + ℰ[1] - λ[[1]] - λ[[2]] - λ[[3]] - μ[[1]] - μ[[2]] - μ[[3]]) - μ[[2]]*μ[[3]] - μ[[1]]*μ[[2]]*μ[[3]]; hh[2, 1, 1] = -ℰ[1]^2 - ℰ[2]^2 + 2*λ[[1]] + λ[[2]] - λ[[1]]*λ[[2]] - λ[[1]]*λ[[3]] - λ[[2]]*λ[[3]] + 2*μ[[1]] + μ[[2]] - μ[[1]]*μ[[2]] - μ[[1]]*μ[[3]] - μ[[2]]*μ[[3]] + ℰ[2]*(-1 + λ[[1]] + λ[[2]] + λ[[3]] + μ[[1]] + μ[[2]] + μ[[3]]) + ℰ[1]*(-2 - ℰ[2] + λ[[1]] + λ[[2]] + λ[[3]] + μ[[1]] + μ[[2]] + μ[[3]]); ℰ[1] = findweight[λ, μ][ν][[1]]; ℰ[2] = findweight[λ, μ][ν][[2]]; ℰ[3] = findweight[λ, μ][ν][[3]]; M[1] = λ[[1]] + μ[[1]] - ℰ[1]; M[2] = λ[[1]] + μ[[1]] - ℰ[2] + 1; M[3] = λ[[1]] + μ[[1]] - ℰ[3] + 2; l[u_] = 1/((-1 + u)^λ[[1]]*u^μ[[1]]); Q[k_][u_] = l[u]*q[k][u]; q[1][u_] = u^M[1] + Sum[c[1, n]*u^n, {n, 0, M[1] - 1}]; q[2][u_] = u^M[2] + Sum[c[2, n]*u^n*Boole[n != M[1]], {n, 0, M[2] - 1}]; q[3][u_] = u^M[3] + Sum[c[3, n]*u^n*Boole[n != M[2]]*Boole[n != M[1]], {n, 0, M[3] - 1}]; Bax[k_] = Numerator[Together[Derivative[1][Q[k]][u]*(h[2][u] + Derivative[1][h[1]][u]) + (h[3][u] + (1/2)*Derivative[1][h[2]][u] + (1/3)*Derivative[2][h[1]][u])*Q[k][u] + h[1][u]*Derivative[2][Q[k]][u] + Derivative[3][Q[k]][u]]]; pencilSolve[{Bax[1], Bax[2], Bax[3]}, {q[1][u], q[2][u], q[3][u]}, V]];
