## Deleted Cell 7 — 2026-07-21 11:59:33
**Kind:** code
```wolfram
ClearAll[occ, wt, evec, La, polyFac];
(* GT basis, weights, and transfer-matrix eigenvalues La[a,oc][u], L=1 rep [1,1,0,0] *)
basis    = ps[{{1,1,0,0}}];
wt[s_]   := Table[Diagonal[Ee[i,i]][[s]], {i,4}];     (* weight (E11,E22,E33,E44) *)
occ[s_]  := Flatten@Position[wt[s], 1];               (* occupied pair {i,j}      *)
evec[s_] := UnitVector[6, s];                         (* eigenvector             *)
polyFac[Ii_, O_][u_] := Module[{e = Length@Complement[Ii,O]},
   Product[u-θ1-s h, {s,0,e-1}] Product[u-θ1-s h, {s,e+1,Length[Ii]}]];
La[a_, O_][u_] := Sum[(Times@@(z/@Ii)) polyFac[Ii,O][u], {Ii, Subsets[{1,2,3,4},{a}]}];
```

## Deleted Cell 10 — 2026-07-21 11:59:45
**Kind:** code
```wolfram
ClearAll[qpart, Qfun];
(* single-index Baxter Q on state oc={i,j}: Qfun[oc,k][u] = z[k]^(u/h)*qpart, monic *)
qpart[oc_, k_][uu_] := If[MemberQ[oc, k], 1,
   (uu - θ1) + h (z[oc[[1]]] z[oc[[2]]] - z[k]^2)/((z[k]-z[oc[[1]]]) (z[k]-z[oc[[2]]]))];
Qfun[oc_, k_][uu_] := z[k]^(uu/h) qpart[oc, k][uu];
```

