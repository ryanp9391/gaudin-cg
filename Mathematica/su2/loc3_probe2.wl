(* assumes loc3_probe.wl already loaded (q-symbols defined) *)
qcb = qcoh[qrb]; qcc = qcoh[qrc]; qca = qcoh[qra];

(* strip a factorized site-3 (resp site-1) and return the remaining 4x4 block *)
qstrip3[v_, cb_] := ArrayReshape[ArrayReshape[v, {qd^2, qd}][[All, qd]]/cb[[qd]], {qd, qd}];
qstrip1[v_, ca_] := ArrayReshape[ArrayReshape[v, {qd, qd^2}][[qd]]/ca[[qd]], {qd, qd}];

(* --- E1: composite of the two ends on the global state --- *)
qE1 = Table[With[{v = qst[qra, qra, qra].qA[s][qw[[1]]].qB[sp][qw[[3]] + qh qS]},
    {s, sp, qrk1[v], qrk3[v]}], {s, 0, 3}, {sp, 0, 3}];

(* --- E2: (a,b,b), A@w2 : span of 4 states, site3 stripped, look for rank-1 combos --- *)
qM2 = Table[qstrip3[qst[qra, qrb, qrb].qA[s][qw[[2]]], qcb], {s, 0, 3}];
qcv = {c0, c1, c2, c3};
qN2 = Sum[qcv[[s + 1]] qM2[[s + 1]], {s, 0, 3}];
qmin2 = Flatten[Table[Det[qN2[[{i1, i2}, {j1, j2}]]], {i1, 4}, {i2, i1 + 1, 4}, {j1, 4}, {j2, j1 + 1, 4}]];
qsol2 = Solve[Join[Thread[qmin2 == 0], {c0 + c1 + c2 + c3 == 1}], qcv];

(* --- E3: (a,a,c), B@w2+hS : span of 4 states, site1 stripped, look for rank-1 combos --- *)
qM3 = Table[qstrip1[qst[qra, qra, qrc].qB[s][qw[[2]] + qh qS], qca], {s, 0, 3}];
qN3 = Sum[qcv[[s + 1]] qM3[[s + 1]], {s, 0, 3}];
qmin3 = Flatten[Table[Det[qN3[[{i1, i2}, {j1, j2}]]], {i1, 4}, {i2, i1 + 1, 4}, {j1, 4}, {j2, j1 + 1, 4}]];
qsol3 = Solve[Join[Thread[qmin3 == 0], {c0 + c1 + c2 + c3 == 1}], qcv];

{"E1 composite ends {s,sp,rk(1|23),rk(12|3)}" -> Flatten[qE1, 1],
 "E2 rank-1 combos in span A@w2 on (a,b,b)" -> qsol2,
 "E3 rank-1 combos in span B@w2+hS on (a,a,c)" -> qsol3,
 "E2 span dim" -> MatrixRank[Flatten /@ qM2], "E3 span dim" -> MatrixRank[Flatten /@ qM3]}
