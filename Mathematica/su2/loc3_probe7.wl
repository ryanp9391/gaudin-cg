(* end-to-end: <coh(q),coh(p),coh(r)|Psi>  ==  (Sum a[mu,s] tau_mu(w1) tau_s(w2)) <Omega(r)|Psi> *)
qu0 = 17/13;
qevec = Transpose[Eigenvectors[N[qA1[qu0], 40]]];      (* columns = right eigenvectors *)
qnev = Length[qevec[[1]]];

(* joint-eigenvalue check + the identity, for every eigenvector *)
qtest = Table[
   Module[{psi = qevec[[All, k]], tau1, tau2, lhs, rhs, ok},
    tau1 = Table[(qA[mu][qw[[1]]].psi)[[1]]/psi[[1]], {mu, 0, 3}];
    tau2 = Table[(qA[s][qw[[2]]].psi)[[1]]/psi[[1]], {s, 0, 3}];
    ok = Max[Abs[Flatten[{Table[qA[mu][qw[[1]]].psi - tau1[[mu + 1]] psi, {mu, 0, 3}],
         Table[qA[s][qw[[2]]].psi - tau2[[s + 1]] psi, {s, 0, 3}]}]]];
    lhs = qst[qq, qp, qra].psi;
    rhs = Sum[qsolA[[4 mu + s + 1]] tau1[[mu + 1]] tau2[[s + 1]], {mu, 0, 3}, {s, 0, 3}] (qst[qra, qra, qra].psi);
    {ok, Abs[lhs - rhs]/Max[Abs[lhs], 10^-30]}], {k, qnev}];

{"max joint-eigenvector residual over all 64 states" -> Max[qtest[[All, 1]]],
 "max relative error of the overlap identity" -> Max[qtest[[All, 2]]],
 "dim span for L=2-style single insertion at w1" ->
   MatrixRank[Table[ArrayReshape[qst[qra, qra, qra].qA[mu][qw[[1]]], {qd, qd^2}][[All, qd^2]], {mu, 0, 3}]]}
