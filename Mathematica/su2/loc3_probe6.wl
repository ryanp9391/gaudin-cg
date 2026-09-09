(* Explicit construction of the fully independent 3-parameter state from the global one *)
qca = qcoh[qra]; qglob = qst[qra, qra, qra];
qVv = Flatten[Table[qglob.qA[mu][qw[[1]]].qA[s][qw[[2]]], {mu, 0, 3}, {s, 0, 3}], 1];
qidx = Flatten[Table[{mu, s}, {mu, 0, 3}, {s, 0, 3}], 1];

(* (i) every one of the 16 products leaves site 3 exactly at coh[a] *)
qnz[m_] := m[[First[Flatten[Position[Map[Norm, m], x_ /; x =!= 0, {1}, 1]]]]];
qsite3 = Table[With[{m = ArrayReshape[v, {qd^2, qd}]},
    {MatrixRank[m], MatrixRank[{qnz[m], qca}]}], {v, qVv}];

(* (ii) strip site 3, solve for the coefficients hitting a generic coherent pair *)
qstrip = Table[ArrayReshape[v, {qd^2, qd}][[All, qd]]/qca[[qd]], {v, qVv}];
qq = 13/6; qp = 5/3;            (* two NEW parameters, both different from a = 3/7 *)
qtgt = Flatten[KroneckerProduct[qcoh[qq], qcoh[qp]]];
qsolA = LinearSolve[Transpose[qstrip], qtgt];

(* (iii) exact check in the full 64-dim space *)
qbuilt = Sum[qsolA[[i]] qVv[[i]], {i, 16}];
qcheck = Simplify[qbuilt - qst[qq, qp, qra]];

{"ranks {12|3, is-coh[a]} for the 16 products" -> Union[qsite3],
 "target (q,p,r)" -> {qq, qp, qra},
 "coefficient matrix a[mu,s]" -> ArrayReshape[qsolA, {4, 4}],
 "rank of coefficient matrix" -> MatrixRank[ArrayReshape[qsolA, {4, 4}]],
 "built state - coh[q](x)coh[p](x)coh[r]" -> Union[qcheck]}
