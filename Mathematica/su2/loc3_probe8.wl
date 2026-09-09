(* a[mu,s] symbolically in the two target parameters *)
qM = Transpose[qstrip];                                  (* rows indexed by (mu,s), stripped of site 3 *)
qasym = LinearSolve[qM, Flatten[KroneckerProduct[qcoh[qX], qcoh[qY]]]];
qaArr = ArrayReshape[Expand /@ qasym, {4, 4}];

qdeg = Map[{Exponent[#, qX], Exponent[#, qY]} &, qaArr, {2}];
qfac = Map[Factor, qaArr, {2}];

(* one-site tower: eigenvalues of t1[{1,s}][w1], to identify the polynomials p_s *)
qTt1[{i_, j_}][u_] := qL[{i, j}][u - qw[[1]]]/(u - qw[[1]] - qh qLam[[1]]);
qa1[u_] := (qz[[1]] + qz[[2]]) qTt1[{1, 1}][u] - qz[[1]] qz[[2]] qTt1[{2, 1}][u] + qTt1[{1, 2}][u];
qa2[u_] := qz[[1]] qz[[2]] (qTt1[{1, 1}][u].qTt1[{2, 2}][u - qh] - qTt1[{2, 1}][u].qTt1[{1, 2}][u - qh]);
qt1[0][u_] := IdentityMatrix[qd]; qt1[1][u_] := qa1[u];
qt1[s_][u_] := qt1[s - 1][u].qa1[u + (s - 1) qh] - qa2[u + (s - 1) qh].qt1[s - 2][u];
qrhoG = qz[[1]] qz[[2]] 0 + Sum[{{qz[[1]] + qz[[2]], -qz[[1]] qz[[2]]}, {1, 0}}[[i, j]] qE[{i, j}], {i, 2}, {j, 2}];
qev = Sort[Eigenvalues[qrhoG]];
qspec = Table[Sort[Eigenvalues[qt1[s][qw[[1]]]]], {s, 0, 3}];

{"a[mu,s] bidegree in (q,p)" -> qdeg,
 "a[mu,s] factored" -> qfac,
 "eigenvalues of rho(G) in the spin-S rep" -> qev,
 "eigenvalues of t1[{1,s}][w1], s=0..3" -> qspec}
