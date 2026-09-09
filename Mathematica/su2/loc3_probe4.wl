(* scan: which single ops / pairs of ops leave which sites untouched, starting from the GLOBAL state *)
qmid[v_] := MatrixRank[Flatten /@ Transpose[ArrayReshape[v, {qd, qd, qd}], {2, 1, 3}]];
qr1[v_] := MatrixRank[ArrayReshape[v, {qd, qd^2}]];
qr3[v_] := MatrixRank[ArrayReshape[v, {qd^2, qd}]];
qsig[v_] := {qr1[v], qmid[v], qr3[v]};   (* rank in split site1|rest , site2|rest , site3|rest *)

qglob = qst[qra, qra, qra];
qops = Join[
   Table[{"A@w" <> ToString[k] <> " s" <> ToString[s], qA[s][qw[[k]]]}, {k, 3}, {s, 1, 3}],
   Table[{"B@w" <> ToString[k] <> "+hS s" <> ToString[s], qB[s][qw[[k]] + qh qS]}, {k, 3}, {s, 1, 3}]];
qops = Flatten[qops, 2];
qops = Partition[qops, 2];

qsingle = Table[{qops[[i, 1]], qsig[qglob.qops[[i, 2]]]}, {i, Length[qops]}];

qpairs = Flatten[Table[
    {qops[[i, 1]] <> " * " <> qops[[j, 1]], qsig[qglob.qops[[i, 2]].qops[[j, 2]]]},
    {i, Length[qops]}, {j, Length[qops]}], 1];
(* keep only pairs where at least two of the three sites stay factorized *)
qgood = Select[qpairs, Count[#[[2]], 1] >= 2 &];

{"legend" -> "{rank site1|rest, rank site2|rest, rank site3|rest}; 1 = that site factorizes off",
 "singles" -> qsingle,
 "pairs with >=2 factorized sites" -> qgood}
