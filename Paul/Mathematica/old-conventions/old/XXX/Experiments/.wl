(* ---- 2. Numerical simultaneous eigenstates and product reference vectors. ---- *)
ClearAll[precisionL2, weightList2, sectors2, rho2, stateRecords2,
  rightMatrix2, leftMatrix2, vacuumState2, alphaConstraints2,
  betaRules2, alphaValues2, betaValues2, bra01, ket01, bra02, ket02,
  rightReferencePhysical2, leftReferencePhysical2, physicalScalar2];

precisionL2 = 90;
weightList2 = Transpose[Table[Diagonal[Etot2[i, i]], {i, 3}]];
sectors2 = Union[weightList2];
rho2[a_, k_] := Prime[7 a + k + 1]/(11 + 2 a + k);

stateRecords2 = Flatten[Table[
   Module[{idx, K, vals, vecs, ord},
    idx = Flatten@Position[weightList2, wt];
    K = Sum[rho2[a, k] N[Tcoeff2[a, k][[idx, idx]], precisionL2],
      {a, 2}, {k, 0, 2 a}];
    {vals, vecs} = Eigensystem[K];
    ord = Ordering[({Re[#], Im[#]} &) /@ vals];
    Table[
     Module[{vv = vecs[[ord[[r]]]], p, t1, t2, full, eigResidual},
      p = First@Ordering[-Abs[vv], 1];
      vv = vv/vv[[p]];
      t1 = Sum[
        ((N[Tcoeff2[1, k][[idx, idx]], precisionL2].vv)[[p]]/
           vv[[p]]) u^k, {k, 0, 2}];
      t2 = Sum[
        ((N[Tcoeff2[2, k][[idx, idx]], precisionL2].vv)[[p]]/
           vv[[p]]) u^k, {k, 0, 4}];
      full = ReplacePart[ConstantArray[SetPrecision[0, precisionL2], dchain],
        Thread[idx -> vv]];
      eigResidual = Max[Flatten@Table[
         Abs[N[Tcoeff2[a, k], precisionL2].full
          - Coefficient[{t1, t2}[[a]], u, k] full],
         {a, 2}, {k, 0, 2 a}]];
      <|"mu" -> wt, "vec" -> full, "T1" -> Expand[t1],
        "T2" -> Expand[t2], "eigResidual" -> eigResidual|>],
     {r, Length[ord]}]],
   {wt, sectors2}], 1];

stateRecords2 = SortBy[stateRecords2,
   {(-#["mu"][[1]]) &, (-#["mu"][[2]]) &, (-#["mu"][[3]]) &}];
rightMatrix2 = Transpose[stateRecords2[[All, "vec"]]];
leftMatrix2 = Inverse[rightMatrix2];
vacuumState2 = First@FirstPosition[stateRecords2[[All, "mu"]], {4, 2, 0}];

alphaConstraints2 = {
  -zTwist[[2]] (zTwist[[1]] - zTwist[[3]])/
    (zTwist[[1]] (zTwist[[2]] - zTwist[[3]])),
  -(zTwist[[1]] - zTwist[[3]])/(zTwist[[2]] - zTwist[[3]])};
betaRules2 = Table[{
    b21 -> -zTwist[[1]] (zTwist[[1]] - zTwist[[3]])/
      (a12 zTwist[[2]] (zTwist[[2]] - zTwist[[3]])),
    b31 -> {-zTwist[[1]]^2/(a13 zTwist[[3]]^2),
       -zTwist[[1]]/(a13 zTwist[[3]])}[[q]],
    b32 -> -zTwist[[2]]^2 (zTwist[[1]] - zTwist[[3]])/
      (a23 (zTwist[[1]] - zTwist[[2]]) zTwist[[3]]^2)},
   {q, 2}];
alphaValues2[q_] := {a13 -> 3/5, a23 -> 4/7,
   a12 -> alphaConstraints2[[q]] (3/5)/(4/7)};
betaValues2[q_] := betaRules2[[q]] /. alphaValues2[q];

bra01[q_] := N[Dot[
   UnitVector[dsite, 1],
   MatrixExp[(a12 /. alphaValues2[q]) EElocal[1, 2]],
   MatrixExp[(a13 /. alphaValues2[q]) EElocal[1, 3]],
   MatrixExp[(a23 /. alphaValues2[q]) EElocal[2, 3]]],
   precisionL2];
ket01[q_] := N[Dot[
   MatrixExp[(b21 /. betaValues2[q]) EElocal[2, 1]],
   MatrixExp[(b31 /. betaValues2[q]) EElocal[3, 1]],
   MatrixExp[(b32 /. betaValues2[q]) EElocal[3, 2]],
   UnitVector[dsite, 1]], precisionL2];
bra02[q_] := Flatten[KroneckerProduct[bra01[q], bra01[q]]];
ket02[q_] := Flatten[KroneckerProduct[ket01[q], ket01[q]]];
rightReferencePhysical2[q_] := bra02[q].rightMatrix2;
leftReferencePhysical2[q_] := leftMatrix2.ket02[q];
physicalScalar2[q_] := IdentityMatrix[dchain]/
  Outer[Times, leftReferencePhysical2[q], rightReferencePhysical2[q]];

Grid[{
  {"number of states", Length[stateRecords2]},
  {"weight multiplicities", Counts[stateRecords2[[All, "mu"]]]},
  {"max simultaneous-eigenvector residual",
   N[Max[stateRecords2[[All, "eigResidual"]]], 6]},
  {"biorthogonal inverse residual",
   N[Max[Abs[leftMatrix2.rightMatrix2 - IdentityMatrix[dchain]]], 6]},
  {"vacuum state index and weight",
   {vacuumState2, stateRecords2[[vacuumState2, "mu"]]}},
  {"branch I: a12 a23/a13, beta rules",
   {alphaConstraints2[[1]], betaRules2[[1]]}},
  {"branch II: a12 a23/a13, beta rules",
   {alphaConstraints2[[2]], betaRules2[[2]]}}
  }, Frame -> All, Alignment -> Left]