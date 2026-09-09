(* TARGET: chi (x) coh[p] (x) coh[a] with p != a, starting from the global state (a,a,a).
   Route: A-tower at w1 (site 1 free) then A-tower at w2 (site 3 provably untouched). *)

qca = qcoh[qra];
qglob = qst[qra, qra, qra];

(* two-insertion family *)
qV = Flatten[Table[qglob.qA[mu][qw[[1]]].qA[s][qw[[2]]], {mu, 0, 3}, {s, 0, 3}], 1];

(* is site 3 always still coh[a]? *)
qsite3ok = Table[With[{m = ArrayReshape[v, {qd^2, qd}]},
    MatrixRank[m] == 1 && MatrixRank[{m[[First[Flatten[Position[Norm /@ m, x_ /; x =!= 0, {1}, 1]]]]], qca}] == 1],
   {v, qV}];

(* strip site 3 -> 16-vectors on sites (1,2) *)
qstrip = Table[ArrayReshape[v, {qd^2, qd}][[All, qd]]/qca[[qd]], {v, qV}];
qdimV = MatrixRank[qstrip];

(* the target subspace: chi (x) coh[p], chi free  ->  4-dim, for given p *)
qW[p_] := Table[Flatten[KroneckerProduct[UnitVector[qd, i], qcoh[p]]], {i, qd}];

(* is chi (x) coh[p] reachable for a generic p != a ? *)
qptest = 5/3;
qdimInt = MatrixRank[qstrip] + MatrixRank[qW[qptest]] - MatrixRank[Join[qstrip, qW[qptest]]];

{"site3 still coh[a] for all 16 products?" -> Union[qsite3ok],
 "dim span on sites (1,2)" -> qdimV,
 "dim( span  cap  C^4 (x) coh[5/3] )" -> qdimInt}
