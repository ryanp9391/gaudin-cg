(* ::Package:: *)

ClearAll[schur, tensorDecomp, dimGL4, irreps];

(*schur polynomial in x[1],x[2],x[3],x[4]*)
schur[lam_, x_List] := Cancel[Det[Table[x[[i]]^(lam[[j]] + Length[x] - j), {i, Length[x]}, {j, Length[x]}]]/Det[Table[x[[i]]^(Length[x] - j), {i, Length[x]}, {j, Length[x]}]]];

(*dimension of the representation lam = {\[Lambda]1,\[Lambda]2,\[Lambda]3,\[Lambda]4}*)
dimGL4[lam_] := Product[(lam[[i]] - lam[[j]] + j - i)/(j - i), {i, 4}, {j, i + 1, 4}];

(*decompose the tensor product \[Lambda] x \[Mu] into irreps*)
tensorDecomp[lam_, mu_] := Module[{x = Array[xx, 4], prod, rules, nu, mult, out = {}},
   prod = Expand[schur[lam, x] schur[mu, x]];
   While[prod =!= 0, rules = Association[CoefficientRules[prod, x]]; nu = Last[Sort[Keys[rules]]];
     mult = rules[nu]; AppendTo[out, nu -> mult]; prod = Expand[prod - mult schur[nu, x]]];
   out];

(* the irreps mu (all multiplicity 1 for this family) appearing in [S,S,0,0] x [P,P,0,0] *)
irreps[S_, P_] := tensorDecomp[{S, S, 0, 0}, {P, P, 0, 0}][[All, 1]];


(*example - all irreps appearing in the tensor product of [2,2,0,0] x [3,3,0,0]*)
irreps[2,3]


(* --- Casimir eigenvalues on the output irrep mu, and the Q-polynomial degrees ---
   Column-ordered generating function (Kazakov et al conventions): G4[u] = Prod_i (u - h(mu_i+4-i)),
   whose coefficients are the Casimirs. Eigenvalues are symmetric functions of the centered
   shifted weights xw_i = mu_i + 5/2 - i (= mu_i + rho_i); Pw_k = Sum_i xw_i^k. *)
ClearAll[muhat, xw, Pw, Cas2, Cas4, degQ];
muhat[mu_] := Table[mu[[i]] + 4 - i, {i, 4}];
xw[mu_] := Table[mu[[i]] + 5/2 - i, {i, 4}];
Pw[mu_, k_] := Total[xw[mu]^k];
Cas2[mu_] := Pw[mu, 2] - 5;
Cas4[mu_] := Pw[mu, 4] + 4 Pw[mu, 3] - Pw[mu, 1] Pw[mu, 2] + Pw[mu, 2]/2 - 2 Pw[mu, 1]^2 - 10 Pw[mu, 1] - 51/4;
(* the SAME shifted weights fix the Q-degrees: deg q_k = (S+P+3) - muhat_k = (S+P) - mu_k + (k-1)
   -- convention: arguments always ordered (S,P,mu) *)
degQ[SS_, PP_, mu_] := Table[(SS + PP + 3) - muhat[mu][[k]], {k, 4}]


(* --- Gaudin ODE --- 
   Q'''' coefficient is
   8(z-w1)^2(z-w2)^2 -- double singular points at the two Gaudin sites, as expected. C2,C4 are
   the abstract Casimirs of the output irrep mu; S,P label the two input reps [S,S,0,0],[P,P,0,0]*)
ClearAll[gaudinODE];
gaudinODE = (2*C2 + C2^2 - 2*C4 + 24*P + 2*C2*P + 20*P^2 + 4*P^3 + 24*S + 2*C2*S + 40*P*S + 12*P^2*S + 20*S^2 + 12*P*S^2 + 4*S^3)*Q[z]
   + 4*(C2 - 4*P - 2*P^2 - 4*S - 4*P*S - 2*S^2)*(2*z + P*z + S*z - w1 - P*w1 - w2 - S*w2)*Derivative[1][Q][z]
   - 4*(C2*z^2 - 6*P*z^2 - 4*P^2*z^2 - 6*S*z^2 - 8*P*S*z^2 - 4*S^2*z^2 - C2*z*w1 + 8*P*z*w1 + 6*P^2*z*w1 + 4*S*z*w1 + 8*P*S*z*w1
      + 2*S^2*z*w1 - 2*P*w1^2 - 2*P^2*w1^2 - C2*z*w2 + 4*P*z*w2 + 2*P^2*z*w2 + 8*S*z*w2 + 8*P*S*z*w2 + 6*S^2*z*w2
      + C2*w1*w2 - 4*P*w1*w2 - 2*P^2*w1*w2 - 4*S*w1*w2 - 8*P*S*w1*w2 - 2*S^2*w1*w2 - 2*S*w2^2 - 2*S^2*w2^2)*Derivative[2][Q][z]
   - 16*(z - w1)*(z - w2)*(P*z + S*z - P*w1 - S*w2)*Derivative[3][Q][z]
   + 8*(z - w1)^2*(z - w2)^2*Derivative[4][Q][z];
(* punctures fixed at w1=0, w2=1 for this solver *)
wsub = {w1 -> 0, w2 -> 1};



