qcohQ[v_] := Module[{u = v/v[[qd]]}, {u[[3]], Simplify[u[[2]] - u[[3]]^2/2], Simplify[u[[1]] - u[[3]]^3/6]}];
qfactor[mm_] := Module[{i, row, j, col},
   i = First[Flatten[Position[Map[Total[Abs[#]] &, mm], x_ /; x =!= 0, {1}, 1]]];
   row = mm[[i]];
   j = First[Flatten[Position[Abs[row], x_ /; x =!= 0, {1}, 1]]];
   col = mm[[All, j]];
   {col/col[[qd]], row/row[[qd]]}];

qc3 = {136/89, -47/89, 0, 0};
qN3s = Sum[qc3[[s + 1]] qM3[[s + 1]], {s, 0, 3}];
qf3s = qfactor[qN3s];
qc2 = {226/443, 217/443, 0, 0};
qN2s = Sum[qc2[[s + 1]] qM2[[s + 1]], {s, 0, 3}];
qf2s = qfactor[qN2s];

{"a,b,c" -> {qra, qrb, qrc},
 "E3: state (a,a,c), op = c0 + c1 tB[1](w2+hS), site1=coh[a]" -> "",
 "  site2 factor" -> qf3s[[1]], "  site2 coherent?" -> qcohQ[qf3s[[1]]],
 "  site3 factor" -> qf3s[[2]], "  site3 coherent?" -> qcohQ[qf3s[[2]]],
 "E2: state (a,b,b), op = c0 + c1 tA[1](w2), site3=coh[b]" -> "",
 "  site1 factor" -> qf2s[[1]], "  site1 coherent?" -> qcohQ[qf2s[[1]]],
 "  site2 factor" -> qf2s[[2]], "  site2 coherent?" -> qcohQ[qf2s[[2]]]}
