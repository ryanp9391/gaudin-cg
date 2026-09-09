(* private-namespace probe of L=3 localization structure; reuses ps/dim/JJ from the kernel *)
qS = 3; qLam = {qS, 0}; qd = 4; qDD = 64;
qh = 1; qw = {0, 1/5, 7/11}; qz = {2, 5};

qE[{i_, j_}] := qE[{i, j}] = JJ[i, j, qLam];
qL[{i_, j_}][u_] := u IdentityMatrix[qd] KroneckerDelta[i, j] - qh qE[{j, i}];
qT[{i_, j_}][u_] := qT[{i, j}][u] =
   Sum[KroneckerProduct[qL[{i, k}][u - qw[[1]]], qL[{k, l}][u - qw[[2]]], qL[{l, j}][u - qw[[3]]]], {k, 2}, {l, 2}];
qnu1[u_] := Product[u - qw[[n]] - qh qLam[[1]], {n, 3}];
qnu2[u_] := Product[u - qw[[n]] - qh qLam[[2]], {n, 3}];

qTtA[{i_, j_}][u_] := qTtA[{i, j}][u] = qT[{i, j}][u]/qnu1[u];
qA1[u_] := qA1[u] = (qz[[1]] + qz[[2]]) qTtA[{1, 1}][u] - qz[[1]] qz[[2]] qTtA[{2, 1}][u] + qTtA[{1, 2}][u];
qA2[u_] := qA2[u] = qz[[1]] qz[[2]] (qTtA[{1, 1}][u].qTtA[{2, 2}][u - qh] - qTtA[{2, 1}][u].qTtA[{1, 2}][u - qh]);
qA[0][u_] := IdentityMatrix[qDD];
qA[1][u_] := qA1[u];
qA[s_][u_] := qA[s][u] = qA[s - 1][u].qA1[u + (s - 1) qh] - qA2[u + (s - 1) qh].qA[s - 2][u];

qTtB[{i_, j_}][u_] := qTtB[{i, j}][u] = qT[{i, j}][u]/qnu2[u];
qB1[u_] := qB1[u] = (qz[[1]] + qz[[2]]) qTtB[{1, 1}][u] - qz[[1]] qz[[2]] qTtB[{2, 1}][u] + qTtB[{1, 2}][u];
qB2[u_] := qB2[u] = qz[[1]] qz[[2]] (qTtB[{1, 1}][u].qTtB[{2, 2}][u - qh] - qTtB[{2, 1}][u].qTtB[{1, 2}][u - qh]);
qB[0][u_] := IdentityMatrix[qDD];
qB[1][u_] := qB1[u];
qB[s_][u_] := qB[s][u] = qB[s - 1][u].qB1[u - (s - 1) qh] - qB2[u - (s - 2) qh].qB[s - 2][u];

qvac = Table[KroneckerDelta[i, qd], {i, qd}];
qcoh[r_] := qvac.MatrixExp[r qE[{2, 1}]];
qst[r1_, r2_, r3_] := Flatten[KroneckerProduct[qcoh[r1], qcoh[r2], qcoh[r3]]];

qrk1[v_] := MatrixRank[ArrayReshape[v, {qd, qd^2}]];
qrk3[v_] := MatrixRank[ArrayReshape[v, {qd^2, qd}]];
qf1[v_] := ArrayReshape[v, {qd, qd^2}][[All, First[Flatten[Position[Total[Abs[ArrayReshape[v, {qd, qd^2}]]], x_ /; x =!= 0, {1}, 1]]]]];
qf3[v_] := ArrayReshape[v, {qd^2, qd}][[First[Flatten[Position[Total[Abs[Transpose[ArrayReshape[v, {qd^2, qd}]]]], x_ /; x =!= 0, {1}, 1]]]]];
qprop[a_, b_] := MatrixRank[{a, b}] == 1;

qra = 3/7; qrb = 11/5; qrc = 2/9;

qrep[label_, v_, s3_, s1_] := {label, qrk1[v], qrk3[v],
   If[qrk3[v] == 1, qprop[qf3[v], qcoh[s3]], "-"],
   If[qrk1[v] == 1, qprop[qf1[v], qcoh[s1]], "-"]};

qres = {{"case", "rk(1|23)", "rk(12|3)", "site3=coh", "site1=coh"}};
Do[AppendTo[qres, qrep["T1 (a,a,a) A@w1 s=" <> ToString[s], qst[qra, qra, qra].qA[s][qw[[1]]], qra, qra]], {s, 0, 3}];
Do[AppendTo[qres, qrep["T2 (a,b,b) A@w2 s=" <> ToString[s], qst[qra, qrb, qrb].qA[s][qw[[2]]], qrb, qra]], {s, 0, 3}];
Do[AppendTo[qres, qrep["T2b (a,b,c) A@w2 s=" <> ToString[s], qst[qra, qrb, qrc].qA[s][qw[[2]]], qrc, qra]], {s, 0, 3}];
Do[AppendTo[qres, qrep["T3 (a,a,c) B@w2+hS s=" <> ToString[s], qst[qra, qra, qrc].qB[s][qw[[2]] + qh qS], qrc, qra]], {s, 0, 3}];
Do[AppendTo[qres, qrep["T3b (a,a,c) A@w1 s=" <> ToString[s], qst[qra, qra, qrc].qA[s][qw[[1]]], qrc, qra]], {s, 0, 3}];
Do[AppendTo[qres, qrep["T4 (a,a,a) B@w3+hS s=" <> ToString[s], qst[qra, qra, qra].qB[s][qw[[3]] + qh qS], qra, qra]], {s, 0, 3}];
qres
