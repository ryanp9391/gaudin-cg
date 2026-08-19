(* ::Package:: *)

(* ::Title::Closed:: *)
(*RUN: Introduction*)


(* ::Text:: *)
(*The gl(3) representation are labeled by the eigenvalues of the Cartans of the highest weight state:*)
(*Rep1= {\[Lambda][1,1],\[Lambda][1,2],\[Lambda][1,3]}*)
(*Rep2={\[Lambda][2,1],\[Lambda][2,2],\[Lambda][2,3]}*)
(**)
(*We have that:*)
(*\[Lambda][*,1]>=\[Lambda][*,2]>=\[Lambda][*,3]*)


(* ::Input:: *)
(*ClearAll[w1,w2];*)
(*w1=0;*)
(*w2=1;*)


(* ::Input:: *)
(*ClearAll[FUNDAMENTAL, ANTIFUNDAMENTAL, ADJOINT];*)
(*(*Rep of gl(4)*)*)
(*FUNDAMENTAL={1,0,0,0};*)
(*ANTIFUNDAMENTAL={1,1,1,0};*)
(*ADJOINT={2,1,1,0};*)


(* ::Input:: *)
(*Ng=4;*)


(* ::Title::Closed:: *)
(*RUN: GL(4) generators*)


(* ::Subsection::Closed:: *)
(*Gl(4) generators*)


(* ::Input:: *)
(*ps[rep_]:=ps[rep]=Block[{res=(Prepend[rep,#]&/@(Tuples[Table[Range[rep[[1,i+1]],rep[[1,i]]],{i,1,Length[rep[[1]]]-1}]]))},If[Length[res[[1,1]]]==1,Reverse@Sort@res,Reverse@Sort@Flatten[ps/@res,1]]];*)
(*dim[rep_]:=dim[rep]=Length[ps[{rep}]];*)
(*J0[k_][pat_]:=J0[k][pat]=Sum[pat[[k,i]],{i,1,k}]-Sum[pat[[k-1,i]],{i,1,k-1}];*)
(*Jp[k_][pat_]:=Jp[k][pat]=-Sum[Product[pat[[k,i]]-i-pat[[k+1,j]]+j,{j,Range[k+1]}]/Product[pat[[k,i]]-i-pat[[k,j]]+j,{j,Complement[Range[k],{i}]}] v[ReplacePart[pat,{k,i}->pat[[k,i]]+1]],{i,1,k}];*)
(*Jm[k_][pat_]:=Jm[k][pat]=Sum[Product[pat[[k,i]]-i-pat[[k-1,j]]+j,{j,Range[k-1]}]/Product[pat[[k,i]]-i-pat[[k,j]]+j,{j,Complement[Range[k],{i}]}] v[ReplacePart[pat,{k,i}->pat[[k,i]]-1]],{i,1,k}];*)
(*JJ[k_,l_,rep_]:=JJ[k,l,rep]=Which[k==l,DiagonalMatrix[J0[k]/@ps[{rep}]],l==k+1,Transpose[Table[Coefficient[Jp[k][ps[{rep}][[i]]],v[#]]&/@ps[{rep}],{i,1,Length[ps[{rep}]]}]],l==k-1,Transpose[Table[Coefficient[Jm[l][ps[{rep}][[i]]],v[#]]&/@ps[{rep}],{i,1,Length[ps[{rep}]]}]],l>k+1,JJ[k,k+1,rep] . JJ[k+1,l,rep]-JJ[k+1,l,rep] . JJ[k,k+1,rep],l<k-1,JJ[k,k-1,rep] . JJ[k-1,l,rep]-JJ[k-1,l,rep] . JJ[k,k-1,rep]];*)
(*CC[k_,set_,rep_]:=CC[k,set,rep]=Sum[Dot[Sequence@@Table[JJ[j[i+1],j[i],rep],{i,1,k}]]/.j[k+1]->j[1],Evaluate[Sequence@@Array[{j[#],set}&,k]]];*)
(*CCDiag[k_,set_,rep_]:=CCDiag[k,set,rep]=Diagonal[Inverse[Transpose[Eigenvectors[CC[Length[set],set,rep]]]] . CC[k,set,rep] . Transpose[Eigenvectors[CC[Length[set],set,rep]]]]*)


(* ::Subsection::Closed:: *)
(*Hermitian generators*)


(* ::Input:: *)
(*ClearAll[JJH,PM,SolvePH,PH,PElement];*)
(*PM[rep_]:=PM[rep]=DiagonalMatrix[Table[PElement[i],{i,1,dim[rep]}]];*)
(*SolvePH[rep_]:=SolvePH[rep]=FullSimplify[(Solve[(DeleteDuplicates[Table[PM[rep] . JJ[k,l,rep] . Inverse[PM[rep]]-Inverse[PM[rep]] . ConjugateTranspose[JJ[l,k,rep]] . PM[rep],{k,1,Ng},{l,1,Ng}]//Flatten])==0,Table[PElement[i],{i,1,dim[rep]}]])];(*If PM is real*)*)
(*PH[rep_,Optional[M_,1]]:=PH[rep,Optional[M,1]]=(PM[rep]/.SolvePH[rep])[[M]]//Simplify;*)
(*(*If PM is complex*)(*Solve[(DeleteDuplicates[Table[PM[rep].JJ[k,l,rep].Inverse[PM[rep]]-ConjugateTranspose[PM[rep].JJ[l,k,rep].Inverse[PM[rep]]],{k,1,3},{l,1,3}]//Flatten])==0,Table[PElement[i],{i,1,dim[rep]}]];*)*)
(*JJH[k_,l_,rep_,Optional[M_,1]]:=JJH[k,l,rep,Optional[M,1]]=(PH[rep,M] . JJ[k,l,rep] . Inverse[PH[rep,M]])//Simplify;*)


(* ::Subsection::Closed:: *)
(*Definition matrices for Lax connection*)


(* ::Input:: *)
(*ClearAll[e,EM,Id];*)
(*e[a_,b_,rep1_]:=e[a,b,rep1]=SparseArray[{a,b}->1,{dim[rep1],dim[rep1]}];*)
(*EM[a_,b_,rep1_]:=EM[a,b,rep1]=JJ[a,b,rep1];(*Non Hermitian generators*)*)
(*(*JJH[a,b,rep1];*)*)
(*Id[rep1_]:=Id[rep1]=SparseArray[{a_,a_}->1,{dim[rep1],dim[rep1]}];*)


(* ::Text:: *)
(*Generator of overall rotation (E^All)*)


(* ::Input:: *)
(*ClearAll[EMAll];*)
(*EMAll[a_,b_,rep1_,rep2_]:=EMAll[a,b,rep1,rep2]=KroneckerProduct[EM[a,b,rep1],Id[rep2]]+KroneckerProduct[Id[rep1],EM[a,b,rep2]];*)


(* ::Subsection::Closed:: *)
(*Module to find good weights for the rep*)


(* ::Input:: *)
(*ClearAll[FindWeight,FindGoodWeight];*)
(*FindWeight[rep_]:=Table[{EM[1,1,rep][[n,n]],EM[2,2,rep][[n,n]],EM[3,3,rep][[n,n]],EM[4,4,rep][[n,n]]},{n,1,dim[rep]}];*)
(*FindGoodWeight[rep_]:=Pick[FindWeight[rep],Lookup[Counts[FindWeight[rep]],FindWeight[rep]],1]*)


(* ::Input:: *)
(*ClearAll[FindWeight2All,FindWeight2,FindGoodWeight2]*)
(*FindWeight2All[rep1_,rep2_]:=FindWeight2All[rep1,rep2]=Flatten[Table[{EM[1,1,rep1][[n,n]]+EM[1,1,rep2][[m,m]],*)
(*EM[2,2,rep1][[n,n]]+EM[2,2,rep2][[m,m]],*)
(*EM[3,3,rep1][[n,n]]+EM[3,3,rep2][[m,m]],EM[4,4,rep1][[n,n]]+EM[4,4,rep2][[m,m]]},{n,1,dim[rep1]},{m,1,dim[rep2]}],1];*)
(*FindWeight2[rep1_,rep2_]:=FindWeight2[rep1,rep2]=DeleteDuplicates[FindWeight2All[rep1,rep2]];*)
(*FindGoodWeight2[rep1_,rep2_]:=FindGoodWeight2[rep1,rep2]=Pick[FindWeight2All[rep1,rep2],Lookup[Counts[FindWeight2All[rep1,rep2]],FindWeight2All[rep1,rep2]],1]*)


(* ::Subsection::Closed:: *)
(*From Dynkin (\[CapitalDelta],l) label of sl(2) to weights of gl(3)*)


(* ::Input:: *)
(*ClearAll[WeightDynkin,\[Lambda]1,\[Lambda]2,\[Lambda]3];*)
(*WeightDynkin[\[CapitalDelta]_,l_]:=Solve[\[Lambda]1-\[Lambda]2==(\[CapitalDelta])&&\[Lambda]2-\[Lambda]3==l&&\[Lambda]1+\[Lambda]2+\[Lambda]3==0,{\[Lambda]1,\[Lambda]2,\[Lambda]3}]*)
(*(*Fundamental*)*)
(*WeightDynkin[1,0]*)
(*(*AntiFundamental*)*)
(*WeightDynkin[0,1]*)
(*(*Adjoint*)*)
(*WeightDynkin[1,1]*)


(* ::Title::Closed:: *)
(*RUN: Q space*)


(* ::Chapter::Closed:: *)
(*With regulators*)


(* ::Section::Closed:: *)
(*hs*)


(* ::Subsection::Closed:: *)
(*Cartan*)


(* ::Input:: *)
(*ClearAll[h];*)
(*h[rep1_,rep2_][1,1]:=h[rep1,rep2][1,1]=rep1[[1]]+rep1[[2]]+rep1[[3]]+rep1[[4]];*)
(*h[rep1_,rep2_][2,1]:=h[rep1,rep2][2,1]=rep1[[3]] (-1+rep1[[4]])+rep1[[2]] (-2+rep1[[3]]+rep1[[4]])+rep1[[1]] (-3+rep1[[2]]+rep1[[3]]+rep1[[4]]);*)
(*h[rep1_,rep2_][3,1]:=h[rep1,rep2][3,1]=rep1[[2]] (-1+rep1[[3]]) (-2+rep1[[4]])+rep1[[1]] ((-2+rep1[[3]]) (-3+rep1[[4]])+rep1[[2]] (-4+rep1[[3]]+rep1[[4]]));*)
(*h[rep1_,rep2_][4,1]:=h[rep1,rep2][4,1]=rep1[[1]] (-1+rep1[[2]]) (-2+rep1[[3]]) (-3+rep1[[4]]);*)
(**)
(*h[rep1_,rep2_][1,2]:=h[rep1,rep2][1,2]=rep2[[1]]+rep2[[2]]+rep2[[3]]+rep2[[4]];*)
(*h[rep1_,rep2_][2,2]:=h[rep1,rep2][2,2]=rep2[[3]] (-1+rep2[[4]])+rep2[[2]] (-2+rep2[[3]]+rep2[[4]])+rep2[[1]] (-3+rep2[[2]]+rep2[[3]]+rep2[[4]]);*)
(*h[rep1_,rep2_][3,2]:=h[rep1,rep2][3,2]=rep2[[2]] (-1+rep2[[3]]) (-2+rep2[[4]])+rep2[[1]] ((-2+rep2[[3]]) (-3+rep2[[4]])+rep2[[2]] (-4+rep2[[3]]+rep2[[4]]));*)
(*h[rep1_,rep2_][4,2]:=h[rep1,rep2][4,2]=rep2[[1]] (-1+rep2[[2]]) (-2+rep2[[3]]) (-3+rep2[[4]]);*)


(* ::Subsubsection::Closed:: *)
(*z^0 (dependence only on regulators)*)


(* ::Input:: *)
(*(*H1*)*)
(*h[\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][1,3]:=h[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][1,3]=\[Phi]1+\[Phi]2+\[Phi]3+\[Phi]4;*)
(**)
(*(*H2*)*)
(*h[\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][2,5]:=h[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][2,5]=\[Phi]1 \[Phi]2+\[Phi]1 \[Phi]3+\[Phi]2 \[Phi]3+\[Phi]1 \[Phi]4+\[Phi]2 \[Phi]4+\[Phi]3 \[Phi]4;*)
(**)
(*(*H3*)*)
(*h[\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][3,7]:=h[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][3,7]=\[Phi]1 \[Phi]2 \[Phi]3+\[Phi]1 \[Phi]2 \[Phi]4+\[Phi]1 \[Phi]3 \[Phi]4+\[Phi]2 \[Phi]3 \[Phi]4;*)
(**)
(*(*H4*)*)
(*h[\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][4,9]:=h[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][4,9]=\[Phi]1*\[Phi]2*\[Phi]3*\[Phi]4;*)


(* ::Subsection::Closed:: *)
(*Relation (dependence on repAll)*)


(* ::Text:: *)
(*By expanding at z=\[Infinity] and taking the coefficient of 1/z*)


(* ::Input:: *)
(*(*H2*)*)
(*h[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][2,3]:=h[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][2,3]=-h[2,4]+(\[Phi]2+\[Phi]3+\[Phi]4) (repAll[[1]])+(\[Phi]1+\[Phi]3+\[Phi]4) (repAll[[2]])+(\[Phi]1+\[Phi]2+\[Phi]4) (repAll[[3]])+(\[Phi]1+\[Phi]2+\[Phi]3) (repAll[[4]]);*)
(**)
(*(*H3*)*)
(*h[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][3,5]:=h[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][3,5]=-h[3,6]+\[Phi]3 \[Phi]4 (repAll[[1]]+repAll[[2]])+\[Phi]2 \[Phi]4 (repAll[[1]]+repAll[[3]])+\[Phi]1 \[Phi]4 (repAll[[2]]+repAll[[3]])+\[Phi]2 \[Phi]3 (repAll[[1]]+repAll[[4]])+\[Phi]1 \[Phi]3 (repAll[[2]]+repAll[[4]])+\[Phi]1 \[Phi]2 (repAll[[3]]+repAll[[4]]);*)
(**)
(*(*H4*)*)
(*h[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][4,7]:=h[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][4,7]=-h[4,8]+\[Phi]2 \[Phi]3 \[Phi]4 (repAll[[1]])+\[Phi]1 \[Phi]3 \[Phi]4 (repAll[[2]])+\[Phi]1 \[Phi]2 \[Phi]4 (repAll[[3]])+\[Phi]1 \[Phi]2 \[Phi]3 (repAll[[4]]);*)


(* ::Subsection::Closed:: *)
(*Hamiltonian in Q-space*)


(* ::Input:: *)
(*ClearAll[\[Epsilon]1,\[Epsilon]2,\[Epsilon]3];*)
(*\[Epsilon]1[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]1[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=h[rep1,rep2][1,1]/(z-w1)+h[rep1,rep2][1,2]/(z-w2)+h[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][1,3];*)
(*\[Epsilon]2[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]2[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=h[rep1,rep2][2,1]/(z-w1)^2+h[rep1,rep2][2,2]/(z-w2)^2+h[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][2,3]/(z-w1)+h[2,4]/(z-w2)+h[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][2,5];*)
(*\[Epsilon]3[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]3[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=h[rep1,rep2][3,1]/(z-w1)^3+h[rep1,rep2][3,2]/(z-w2)^3+h[3,3]/(z-w1)^2+h[3,4]/(z-w2)^2+h[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][3,5]/(z-w1)+h[3,6]/(z-w2)+h[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][3,7];*)
(*\[Epsilon]4[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]4[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=h[rep1,rep2][4,1]/(z-w1)^4+h[rep1,rep2][4,2]/(z-w2)^4+h[4,3]/(z-w1)^3+h[4,4]/(z-w2)^3+h[4,5]/(z-w1)^2+h[4,6]/(z-w2)^2+h[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][4,7]/(z-w1)+h[4,8]/(z-w2)+h[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][4,9];*)


(* ::Section::Closed:: *)
(*Baxter equation with F*)


(* ::Input:: *)
(*ClearAll[F];*)
(*F[rep1_,rep2_,C1_,*)
(*C2_][z_]:=F[rep1,rep2,C1,C2][z]=Cases[{{C1,C2}},{a_,b_}:>1/((z-w1)^(rep1[[a]]-a+1)*(z-w2)^(rep2[[b]]-b+1))][[1]];*)


(* ::Text:: *)
(*Definition of shifted Hamiltonian*)


(* ::Input:: *)
(*ClearAll[\[Epsilon]1S,\[Epsilon]2S,\[Epsilon]3S]*)
(*\[Epsilon]1S[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_,C2_][z_]:=\[Epsilon]1S[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2][z]=\[Epsilon]1[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+(4 Derivative[1][F[rep1,rep2,C1,C2]][z])/F[rep1,rep2,C1,C2][z];*)
(*\[Epsilon]2S[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_,C2_][z_]:=\[Epsilon]2S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2][z]=\[Epsilon]2[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+(3 \[Epsilon]1[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] Derivative[1][F[rep1,rep2,C1,C2]][z])/F[rep1,rep2,C1,C2][z]+(6 (F[rep1,rep2,C1,C2]^\[Prime]\[Prime])[z])/F[rep1,rep2,C1,C2][z];*)
(*\[Epsilon]3S[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_,C2_][z_]:=\[Epsilon]3S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2][z]=\[Epsilon]3[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+2*\[Epsilon]2[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*Derivative[1][F[rep1,rep2,C1,C2]][z]/F[rep1,rep2,C1,C2][z]+3*(\[Epsilon]1[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] (F[rep1,rep2,C1,C2]^\[Prime]\[Prime])[z])/F[rep1,rep2,C1,C2][z]+4*\!\(\*SuperscriptBox[\(F[rep1, rep2, C1, C2]\), *)
(*TagBox[*)
(*RowBox[{"(", "3", ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[z]/F[rep1,rep2,C1,C2][z];*)
(*\[Epsilon]4S[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_,C2_][z_]:=\[Epsilon]4S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2][z]=\[Epsilon]4[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+\[Epsilon]3[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*Derivative[1][F[rep1,rep2,C1,C2]][z]/F[rep1,rep2,C1,C2][z]+(\[Epsilon]2[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] (F[rep1,rep2,C1,C2]^\[Prime]\[Prime])[z])/F[rep1,rep2,C1,C2][z]+\[Epsilon]1[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*\!\(\*SuperscriptBox[\(F[rep1, rep2, C1, C2]\), *)
(*TagBox[*)
(*RowBox[{"(", "3", ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[z]/F[rep1,rep2,C1,C2][z]+\!\(\*SuperscriptBox[\(F[rep1, rep2, C1, C2]\), *)
(*TagBox[*)
(*RowBox[{"(", "4", ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[z]/F[rep1,rep2,C1,C2][z];*)


(* ::Input:: *)
(*ClearAll[BaxOp];*)
(*BaxOp[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_,C2_]:=D[#,{z,4}]+\[Epsilon]1S[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2][z]*D[#,{z,3}]+\[Epsilon]2S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2][z]*D[#,{z,2}]+\[Epsilon]3S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2][z]*D[#,z]+\[Epsilon]4S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2][z]*#&;*)


(* ::Section::Closed:: *)
(*Q functions*)


(* ::Input:: *)
(*ClearAll[QPoly,Q,c1,c2,c3,c4];*)
(*QPoly[c_,M_,r_][z_]:=QPoly[c,M,r][z]=Exp[r z]*(z^M+Sum[c[k]*z^k,{k,0,M-1}]);*)
(*Q[1][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q[1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly[c1,rep1[[1]]+rep2[[1]]-repAll[[1]],-\[Phi]1][z];*)
(*Q[2][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q[2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly[c2,rep1[[1]]+rep2[[1]]-repAll[[2]],-\[Phi]2][z];*)
(*Q[3][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q[3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly[c3,rep1[[1]]+rep2[[1]]-repAll[[3]],-\[Phi]3][z];*)
(*Q[4][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q[4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly[c4,rep1[[1]]+rep2[[1]]-repAll[[4]],-\[Phi]4][z];*)


(* ::Section::Closed:: *)
(*Baxter equations*)


(* ::Text:: *)
(*Baxter equations*)


(* ::Input:: *)
(*ClearAll[Bax];*)
(*Bax[1][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_,C2_]:=Bax[1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2]=((BaxOp[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2][Q[1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[\[Phi]1 z])//Simplify;*)
(*Bax[2][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_,C2_]:=Bax[2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2]=((BaxOp[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2][Q[2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[\[Phi]2 z])//Simplify;*)
(*Bax[3][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_,C2_]:=Bax[3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2]=((BaxOp[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2][Q[3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[\[Phi]3 z])//Simplify;*)
(*Bax[4][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_,C2_]:=Bax[4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2]=((BaxOp[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2][Q[4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[\[Phi]4 z])//Simplify;*)


(* ::Text:: *)
(*Unknown coefficients*)


(* ::Input:: *)
(*ClearAll[Vars1,Vars2,Vars3,Vars4];*)
(*Vars1[rep1_,rep2_,repAll_]:=Vars1[rep1,rep2,repAll]=Table[c1[n],{n,0,rep1[[1]]+rep2[[1]]-repAll[[1]]-1}];*)
(*Vars2[rep1_,rep2_,repAll_]:=Vars2[rep1,rep2,repAll]=Table[c2[n],{n,0,rep1[[1]]+rep2[[1]]-repAll[[2]]-1}];*)
(*Vars3[rep1_,rep2_,repAll_]:=Vars3[rep1,rep2,repAll]=Table[c3[n],{n,0,rep1[[1]]+rep2[[1]]-repAll[[3]]-1}];*)
(*Vars4[rep1_,rep2_,repAll_]:=Vars4[rep1,rep2,repAll]=Table[c4[n],{n,0,rep1[[1]]+rep2[[1]]-repAll[[4]]-1}];*)


(* ::Text:: *)
(*Solve Baxter equation*)


(* ::Input:: *)
(*ClearAll[SolveBaxAll];*)
(*SolveBaxAll[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,c1_,c2_]:=SolveBaxAll[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,c1,c2]=Solve[CoefficientList[Bax[1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,c1,c2]//Together//Numerator,z]==0&&*)
(*CoefficientList[Bax[2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,c1,c2]//Together//Numerator,z]==0&&*)
(*CoefficientList[Bax[3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,c1,c2]//Together//Numerator,z]==0&&CoefficientList[Bax[4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,c1,c2]//Together//Numerator,z]==0,Join[Vars1[rep1,rep2,repAll],Vars2[rep1,rep2,repAll],Vars3[rep1,rep2,repAll],Vars4[rep1,rep2,repAll],{h[2,4],h[3,3],h[3,4],h[3,6],h[4,3],h[4,4],h[4,5],h[4,6],h[4,8]}]]//Simplify;*)


(* ::Text:: *)
(*Spectrum of Hamiltonian*)


(* ::Input:: *)
(*ClearAll[HSpectrum];*)
(*HSpectrum[1][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=HSpectrum[1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]1[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]/.SolveBaxAll[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1];*)
(**)
(*HSpectrum[2][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=HSpectrum[2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]2[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]/.SolveBaxAll[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1];*)
(**)
(*HSpectrum[3][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=HSpectrum[3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]3[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]/.SolveBaxAll[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1];*)
(**)
(*HSpectrum[4][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=HSpectrum[4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]4[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]/.SolveBaxAll[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1];*)


(* ::Subsubsection:: *)
(*Example of use*)


(* ::Input:: *)
(*Sol=SolveBaxAll[FUNDAMENTAL,FUNDAMENTAL,{2,0,0,0},1,1/3,1/5,1/7,1,1]*)


(* ::Input:: *)
(*H1ST=HSpectrum[1][FUNDAMENTAL,FUNDAMENTAL,{2,0,0,0},1,1/3,1/5,1/7][z];*)
(*H2ST=HSpectrum[2][FUNDAMENTAL,FUNDAMENTAL,{2,0,0,0},1,1/3,1/5,1/7][z];*)
(*H3ST=HSpectrum[3][FUNDAMENTAL,FUNDAMENTAL,{2,0,0,0},1,1/3,1/5,1/7][z];*)
(*H4ST=HSpectrum[4][FUNDAMENTAL,FUNDAMENTAL,{2,0,0,0},1,1/3,1/5,1/7][z];*)


(* ::Input:: *)
(*SolveBaxAll[ADJOINT,ADJOINT,ADJOINT+ADJOINT,1,1/3,1/5,1/7,1,1]*)


(* ::Chapter::Closed:: *)
(*Mean value*)


(* ::Section::Closed:: *)
(*Measure*)


(* ::Subsection::Closed:: *)
(*2-sites*)


(* ::Text:: *)
(*Measure such that <f Subscript[O, F]g> \!\(\*OverscriptBox[\(=\), \(!\)]\) <g SuperDagger[Subscript[O, F]] f>*)


(* ::Input:: *)
(*ClearAll[\[Mu]];*)
(*\[Mu][rep1_,rep2_][z_]:=\[Mu][rep1,rep2][z]=(-1+z)^(-rep2[[1]]+rep2[[4]]) z^(-rep1[[1]]+rep1[[4]]);*)


(* ::Text:: *)
(*Mean value*)


(* ::Input:: *)
(*ClearAll[MeanValue];*)
(*MeanValue[f_,\[Tau]_,rep1_,rep2_]:=MeanValue[f,\[Tau],rep1,rep2]=*)
(*2 \[Pi] I Residue[(f*\[Mu][rep1,rep2][z]),{z,\[Tau]}];*)


(* ::Subsection::Closed:: *)
(*1-site*)


(* ::Text:: *)
(*Measure such that <f Subscript[O, F]g> \!\(\*OverscriptBox[\(=\), \(!\)]\) <g SuperDagger[Subscript[O, F]] f>*)


(* ::Input:: *)
(*ClearAll[\[Mu]1,\[Mu]2];*)
(*\[Mu]1[rep1_][z_]:=\[Mu]1[rep1][z]=z^(-rep1[[1]]+rep1[[4]]);*)
(*\[Mu]2[rep2_][z_]:=\[Mu]2[rep2][z]=(-1+z)^(-rep2[[1]]+rep2[[4]]);*)


(* ::Text:: *)
(*Mean value*)


(* ::Input:: *)
(*ClearAll[MeanValue1,MeanValue2];*)
(*MeanValue1[f_,\[Tau]_,rep1_]:=MeanValue1[f,\[Tau],rep1]=*)
(*2 \[Pi] I Residue[(f*\[Mu]1[rep1][z]),{z,\[Tau]}];*)
(*MeanValue2[f_,\[Tau]_,rep2_]:=MeanValue2[f,\[Tau],rep2]=*)
(*2 \[Pi] I Residue[(f*\[Mu]2[rep2][z]),{z,\[Tau]}];*)


(* ::Subsubsection:: *)
(*Test*)


(* ::Input:: *)
(*MeanValue[f[z]*1/(z*(z-1)),0,{1,1,0,0},{1,1,0,0}]*)
(*MeanValue[f[z]*1/(z*(z-1)),1,{1,1,0,0},{1,1,0,0}]*)
(*MeanValue1[f[z]*1/z,0,{1,1,0,0}]*)


(* ::Chapter::Closed:: *)
(*ADJOINT BAXTER*)


(* ::Section::Closed:: *)
(*Adjoint Baxter equation with SuperDagger[F]*)


(* ::Text:: *)
(*HERE we consider only the case C1=1,C2=1 in the definition of M*)


(* ::Text:: *)
(*Shifted Hamiltonians (already with SuperDagger[F])*)


(* ::Input:: *)
(*ClearAll[\[Epsilon]1S\[Dagger],\[Epsilon]2S\[Dagger],\[Epsilon]3S\[Dagger]];*)
(*\[Epsilon]1S\[Dagger][rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]1S\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=-\[Epsilon]1S[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z]+4*D[\[Mu][rep1,rep2][z],z]/\[Mu][rep1,rep2][z];*)
(*\[Epsilon]2S\[Dagger][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]2S\[Dagger][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]2S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z]-3*D[\[Epsilon]1S[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z],z]-3*\[Epsilon]1S[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z]*D[\[Mu][rep1,rep2][z],z]/\[Mu][rep1,rep2][z]+6*D[\[Mu][rep1,rep2][z],{z,2}]/\[Mu][rep1,rep2][z];*)
(*\[Epsilon]3S\[Dagger][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]3S\[Dagger][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=-\[Epsilon]3S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z]+2*D[\[Epsilon]2S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z],z]-3*D[\[Epsilon]1S[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z],{z,2}]+D[\[Mu][rep1,rep2][z],z]/\[Mu][rep1,rep2][z]*(2*\[Epsilon]2S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z]-6*D[\[Epsilon]1S[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z],z])-3*D[\[Mu][rep1,rep2][z],{z,2}]/\[Mu][rep1,rep2][z]*\[Epsilon]1S[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z]+4*D[\[Mu][rep1,rep2][z],{z,3}]/\[Mu][rep1,rep2][z];*)
(*\[Epsilon]4S\[Dagger][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]4S\[Dagger][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]4S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z]-D[\[Epsilon]3S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z],z]+D[\[Epsilon]2S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z],{z,2}]-D[\[Epsilon]1S[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z],{z,3}]+D[\[Mu][rep1,rep2][z],z]/\[Mu][rep1,rep2][z]*(-\[Epsilon]3S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z]+2*D[\[Epsilon]2S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z],z]-3*D[\[Epsilon]1S[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z],{z,2}])+D[\[Mu][rep1,rep2][z],{z,2}]/\[Mu][rep1,rep2][z]*(\[Epsilon]2S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z]-3D[\[Epsilon]1S[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z],z])-\[Epsilon]1S[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z]*D[\[Mu][rep1,rep2][z],{z,3}]/\[Mu][rep1,rep2][z]+D[\[Mu][rep1,rep2][z],{z,4}]/\[Mu][rep1,rep2][z];*)


(* ::Input:: *)
(*ClearAll[BaxOp\[Dagger]];*)
(*BaxOp\[Dagger][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=D[#,{z,4}]+\[Epsilon]1S\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*D[#,{z,3}]+\[Epsilon]2S\[Dagger][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*D[#,{z,2}]+\[Epsilon]3S\[Dagger][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*D[#,z]+\[Epsilon]4S\[Dagger][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*#&;*)


(* ::Section::Closed:: *)
(*Q functions*)


(* ::Input:: *)
(*ClearAll[QPoly\[Dagger],Q\[Dagger],c1\[Dagger],c2\[Dagger],c3\[Dagger],c4\[Dagger]];*)
(*QPoly\[Dagger][c_,M_,r_][z_]:=QPoly\[Dagger][c,M,r][z]=Exp[r z]*(z^M+Sum[c[k]*z^k,{k,0,M-1}]);*)
(*Q\[Dagger][1][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q\[Dagger][1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly\[Dagger][c1\[Dagger],-rep1[[4]]-rep2[[4]]+repAll[[1]],\[Phi]1][z];*)
(*Q\[Dagger][2][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q\[Dagger][2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly\[Dagger][c2\[Dagger],-rep1[[4]]-rep2[[4]]+repAll[[2]],\[Phi]2][z];*)
(*Q\[Dagger][3][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q\[Dagger][3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly\[Dagger][c3\[Dagger],-rep1[[4]]-rep2[[4]]+repAll[[3]],\[Phi]3][z];*)
(*Q\[Dagger][4][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q\[Dagger][4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly\[Dagger][c4\[Dagger],-rep1[[4]]-rep2[[4]]+repAll[[4]],\[Phi]4][z];*)


(* ::Section::Closed:: *)
(*Baxter equations*)


(* ::Input:: *)
(*ClearAll[Bax\[Dagger]];*)
(*Bax\[Dagger][1][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=Bax\[Dagger][1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=((BaxOp\[Dagger][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][Q\[Dagger][1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[-\[Phi]1 z])//Simplify;*)
(*Bax\[Dagger][2][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=Bax\[Dagger][2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=((BaxOp\[Dagger][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][Q\[Dagger][2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[-\[Phi]2 z])//Simplify;*)
(*Bax\[Dagger][3][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=Bax\[Dagger][3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=((BaxOp\[Dagger][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][Q\[Dagger][3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[-\[Phi]3 z])//Simplify;*)
(*Bax\[Dagger][4][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=Bax\[Dagger][4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=((BaxOp\[Dagger][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][Q\[Dagger][4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[-\[Phi]4 z])//Simplify;*)


(* ::Text:: *)
(*Unknown coefficients*)


(* ::Input:: *)
(*ClearAll[Vars1\[Dagger],Vars2\[Dagger],Vars3\[Dagger],Vars4\[Dagger]];*)
(*Vars1\[Dagger][rep1_,rep2_,repAll_]:=Vars1\[Dagger][rep1,rep2,repAll]=Table[c1\[Dagger][n],{n,0,-rep1[[4]]-rep2[[4]]+repAll[[1]]-1}];*)
(*Vars2\[Dagger][rep1_,rep2_,repAll_]:=Vars2\[Dagger][rep1,rep2,repAll]=Table[c2\[Dagger][n],{n,0,-rep1[[4]]-rep2[[4]]+repAll[[2]]-1}];*)
(*Vars3\[Dagger][rep1_,rep2_,repAll_]:=Vars3\[Dagger][rep1,rep2,repAll]=Table[c3\[Dagger][n],{n,0,-rep1[[4]]-rep2[[4]]+repAll[[3]]-1}];*)
(*Vars4\[Dagger][rep1_,rep2_,repAll_]:=Vars4\[Dagger][rep1,rep2,repAll]=Table[c4\[Dagger][n],{n,0,-rep1[[4]]-rep2[[4]]+repAll[[4]]-1}];*)


(* ::Subsection::Closed:: *)
(*Solving Baxter equations giving hs*)


(* ::Input:: *)
(*ClearAll[SolveBaxh\[Dagger],SolveBaxAllh\[Dagger]];*)
(*SolveBaxh\[Dagger][1][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,h24_,h33_,h34_,h36_,h43_,h44_,h45_,h46_,h48_]:=SolveBaxh\[Dagger][1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h24,h33,h34,h36,h43,h44,h45,h46,h48]=Solve[CoefficientList[Bax\[Dagger][1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]/.{h[2,4]->h24,h[3,3]->h33,h[3,4]->h34,h[3,6]->h36,*)
(*h[4,3]->h43,h[4,4]->h44,h[4,5]->h45,h[4,6]->h46,h[4,8]->h48}//Together//Numerator,z]==0,Vars1\[Dagger][rep1,rep2,repAll]]//Simplify//Flatten;*)
(**)
(*SolveBaxh\[Dagger][2][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,h24_,h33_,h34_,h36_,h43_,h44_,h45_,h46_,h48_]:=SolveBaxh\[Dagger][2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h24,h33,h34,h36,h43,h44,h45,h46,h48]=Solve[CoefficientList[Bax\[Dagger][2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]/.{h[2,4]->h24,h[3,3]->h33,h[3,4]->h34,h[3,6]->h36,*)
(*h[4,3]->h43,h[4,4]->h44,h[4,5]->h45,h[4,6]->h46,h[4,8]->h48}//Together//Numerator,z]==0,Vars2\[Dagger][rep1,rep2,repAll]]//Simplify//Flatten;*)
(**)
(*SolveBaxh\[Dagger][3][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,h24_,h33_,h34_,h36_,h43_,h44_,h45_,h46_,h48_]:=SolveBaxh\[Dagger][3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h24,h33,h34,h36,h43,h44,h45,h46,h48]=Solve[CoefficientList[Bax\[Dagger][3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]/.{h[2,4]->h24,h[3,3]->h33,h[3,4]->h34,h[3,6]->h36,*)
(*h[4,3]->h43,h[4,4]->h44,h[4,5]->h45,h[4,6]->h46,h[4,8]->h48}//Together//Numerator,z]==0,Vars3\[Dagger][rep1,rep2,repAll]]//Simplify//Flatten;*)
(**)
(*SolveBaxh\[Dagger][4][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,h24_,h33_,h34_,h36_,h43_,h44_,h45_,h46_,h48_]:=SolveBaxh\[Dagger][4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h24,h33,h34,h36,h43,h44,h45,h46,h48]=Solve[CoefficientList[Bax\[Dagger][4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]/.{h[2,4]->h24,h[3,3]->h33,h[3,4]->h34,h[3,6]->h36,*)
(*h[4,3]->h43,h[4,4]->h44,h[4,5]->h45,h[4,6]->h46,h[4,8]->h48}//Together//Numerator,z]==0,Vars4\[Dagger][rep1,rep2,repAll]]//Simplify//Flatten;*)
(**)
(*SolveBaxAllh\[Dagger][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,h24_,h33_,h34_,h36_,h43_,h44_,h45_,h46_,h48_]:=SolveBaxAllh\[Dagger][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h24,h33,h34,h36,h43,h44,h45,h46,h48]=Join[SolveBaxh\[Dagger][1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h24,h33,h34,h36,h43,h44,h45,h46,h48],SolveBaxh\[Dagger][2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h24,h33,h34,h36,h43,h44,h45,h46,h48],SolveBaxh\[Dagger][3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h24,h33,h34,h36,h43,h44,h45,h46,h48],SolveBaxh\[Dagger][4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h24,h33,h34,h36,h43,h44,h45,h46,h48]]//Simplify;*)


(* ::Subsubsection::Closed:: *)
(*Example of use*)


(* ::Input:: *)
(*Solh=(SolveBaxAll[FUNDAMENTAL,FUNDAMENTAL,FUNDAMENTAL+FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1]//Simplify)[[1]]*)
(**)
(*SolveBaxh\[Dagger][1][FUNDAMENTAL,FUNDAMENTAL,FUNDAMENTAL+FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h[2,4]/.Solh,h[3,3]/.Solh,h[3,4]/.Solh,h[3,6]/.Solh,*)
(*h[4,3]/.Solh,h[4,4]/.Solh,h[4,5]/.Solh,h[4,6]/.Solh,h[4,8]/.Solh]*)
(*SolveBaxh\[Dagger][2][FUNDAMENTAL,FUNDAMENTAL,FUNDAMENTAL+FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h[2,4]/.Solh,h[3,3]/.Solh,h[3,4]/.Solh,h[3,6]/.Solh,*)
(*h[4,3]/.Solh,h[4,4]/.Solh,h[4,5]/.Solh,h[4,6]/.Solh,h[4,8]/.Solh]*)
(*SolveBaxh\[Dagger][3][FUNDAMENTAL,FUNDAMENTAL,FUNDAMENTAL+FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h[2,4]/.Solh,h[3,3]/.Solh,h[3,4]/.Solh,h[3,6]/.Solh,*)
(*h[4,3]/.Solh,h[4,4]/.Solh,h[4,5]/.Solh,h[4,6]/.Solh,h[4,8]/.Solh]*)
(*SolveBaxh\[Dagger][4][FUNDAMENTAL,FUNDAMENTAL,FUNDAMENTAL+FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h[2,4]/.Solh,h[3,3]/.Solh,h[3,4]/.Solh,h[3,6]/.Solh,*)
(*h[4,3]/.Solh,h[4,4]/.Solh,h[4,5]/.Solh,h[4,6]/.Solh,h[4,8]/.Solh]*)
(**)
(**)
(*Print["Hi"]*)
(*SolveBaxAllh\[Dagger][FUNDAMENTAL,FUNDAMENTAL,FUNDAMENTAL+FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h[2,4]/.Solh,h[3,3]/.Solh,h[3,4]/.Solh,h[3,6]/.Solh,*)
(*h[4,3]/.Solh,h[4,4]/.Solh,h[4,5]/.Solh,h[4,6]/.Solh,h[4,8]/.Solh]*)


(* ::Section::Closed:: *)
(*Q Hodge*)


(* ::Input:: *)
(*ClearAll[Ng,SetN];*)
(*Ng=4;*)
(*SetN=Range[Ng];*)


(* ::Input:: *)
(*ClearAll[Hodge];*)
(*Hodge[Qfour_,rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,z_,i_]:=Hodge[Qfour,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,z,i]=*)
(*(Exp[z*(\[Phi]1+\[Phi]2+\[Phi]3+\[Phi]4)]*Signature[Join[Complement[SetN,{i}],{i}]]*Wronskian[{Qfour[[Complement[SetN,{i}]]][[1]],Qfour[[Complement[SetN,{i}]]][[2]],Qfour[[Complement[SetN,{i}]]][[3]]},z]*(((z)^-(2*rep1[[1]]-rep1[[2]]-rep1[[3]])*(z-1)^-(2*rep2[[1]]-rep2[[2]]-rep2[[3]])))/((({\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[1]]]]-{\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[2]]]]))*({\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[1]]]]-{\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[3]]]])*({\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[2]]]]-{\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[3]]]])))//Simplify;*)
(*(*The Q argument should be Q1[z], not Q1!*)*)
(**)
(*(**(((z)^(-2(rep1[[1]]-rep1[[2]]))*(z-1)^(-2(rep2[[1]]-rep2[[2]]))))*)*)
(*(*/({\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[1]]]]-{\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[2]]]])*)*)


(* ::Section::Closed:: *)
(*Check Q Hodge*)


(* ::Subsubsection::Closed:: *)
(*Check Table*)


(* ::Input:: *)
(*ClearAll[rep1T,rep2T,repAllT]*)
(*rep1T={3,0,0,0};*)
(*rep2T={4,3,2,1};*)
(*repAllT[i_]:=FindGoodWeight2[rep1T,rep2T][[i]];*)
(*\[Phi]1T=1;*)
(*\[Phi]2T=1/2;*)
(*\[Phi]3T=1/3;*)
(*\[Phi]4T=1/5;*)
(**)
(*(*rep1T=ADJOINT;*)
(*rep2T=ADJOINT;*)
(*repAllT[1]={4,2,0};*)
(*repAllT[2]={3,2,1};*)
(*repAllT[3]={2,2,2};*)
(*\[Phi]1T=3;*)
(*\[Phi]2T=5;*)
(*\[Phi]3T=7;*)*)
(**)
(*(*rep1T={3,2,1};*)
(*rep2T={2,0,0};*)
(*repAllT[1]={5,2,1};*)
(*repAllT[2]={4,3,1};*)
(*repAllT[3]={4,2,2};*)
(*\[Phi]1T=1;*)
(*\[Phi]2T=1/2;*)
(*\[Phi]3T=1/7;*)*)
(**)
(*(*rep1T=FUNDAMENTAL;*)
(*rep2T=FUNDAMENTAL;*)
(*repAllT[1]={1,1,0};*)
(*repAllT[2]={2,0,0};*)
(*repAllT[3]={1,1,0};*)*)
(**)
(*Table[*)
(*Solh=(SolveBaxAll[rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,1,1]//Simplify)[[1]];*)
(**)
(*Qh[1]=Q[1][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh;*)
(*Qh[2]=Q[2][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh;*)
(*Qh[3]=Q[3][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh;*)
(*Qh[4]=Q[4][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh;*)
(**)
(*(*SolveBaxh\[Dagger][1][rep1T,rep2T,repAllT,\[Phi]1,\[Phi]2,\[Phi]3,h[2,5]/.Solh,h[3,5]/.Solh,h[3,7]/.Solh,h[3,9]/.Solh];*)
(*SolveBaxh\[Dagger][2][rep1T,rep2T,repAllT,\[Phi]1,\[Phi]2,\[Phi]3,h[2,5]/.Solh,h[3,5]/.Solh,h[3,7]/.Solh,h[3,9]/.Solh];*)
(*SolveBaxh\[Dagger][3][rep1T,rep2T,repAllT,\[Phi]1,\[Phi]2,\[Phi]3,h[2,5]/.Solh,h[3,5]/.Solh,h[3,7]/.Solh,h[3,9]/.Solh];*)*)
(**)
(**)
(*Print["Hi"];*)
(*Sol\[Dagger]=SolveBaxAllh\[Dagger][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,h[2,4]/.Solh,h[3,3]/.Solh,h[3,4]/.Solh,h[3,6]/.Solh,*)
(*h[4,3]/.Solh,h[4,4]/.Solh,h[4,5]/.Solh,h[4,6]/.Solh,h[4,8]/.Solh];*)
(**)
(*QH[1]=Q\[Dagger][1][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Sol\[Dagger];*)
(*QH[2]=Q\[Dagger][2][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Sol\[Dagger];*)
(*QH[3]=Q\[Dagger][3][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Sol\[Dagger];*)
(*QH[4]=Q\[Dagger][4][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Sol\[Dagger];*)
(**)
(*(*Print["Hodge"]*)
(*Hodge[{Qh[1],Qh[2],Qh[3]},rep1T,rep2T,\[Phi]1,\[Phi]2,\[Phi]3,z,2,1,1]//Simplify*)
(**)
(*Print["Q^\[Dagger]"]*)
(*QH[2]//Simplify*)*)
(**)
(*Print["M=",M,", ",Table[(1/QH[t] Hodge[{Qh[1],Qh[2],Qh[3],Qh[4]},rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,z,t])//Simplify,{t,1,4}]],{M,1,4}];*)


(* ::Subsubsection::Closed:: *)
(*Testing*)


(* ::Input:: *)
(*ClearAll[rep1T,rep2T,repAllT]*)
(*rep1T=ADJOINT;*)
(*rep2T=FUNDAMENTAL;*)
(*repAllT[i_]:=FindGoodWeight2[rep1T,rep2T][[i]];*)
(*\[Phi]1T=1;*)
(*\[Phi]2T=1/2;*)
(*\[Phi]3T=1/5;*)
(*\[Phi]4T=1/7;*)
(*M=3;*)
(*M2=2;*)
(*Solh=(SolveBaxAll[rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,1,1]//Simplify)[[1]];*)
(*Solh2=(SolveBaxAll[rep1T,rep2T,repAllT[M2],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,1,1]//Simplify)[[1]];*)
(**)
(*Qh[1]=Q[1][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh;*)
(*Qh[2]=Q[2][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh;*)
(*Qh[3]=Q[3][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh;*)
(*Qh[4]=Q[4][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh;*)
(**)
(*Qh2[1]=Q[1][rep1T,rep2T,repAllT[M2],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh2;*)
(*Qh2[2]=Q[2][rep1T,rep2T,repAllT[M2],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh2;*)
(*Qh2[3]=Q[3][rep1T,rep2T,repAllT[M2],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh2;*)
(*Qh2[4]=Q[4][rep1T,rep2T,repAllT[M2],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh2;*)
(**)
(*(*SolveBaxh\[Dagger][1][rep1T,rep2T,repAllT,\[Phi]1,\[Phi]2,\[Phi]3,h[2,5]/.Solh,h[3,5]/.Solh,h[3,7]/.Solh,h[3,9]/.Solh];*)
(*SolveBaxh\[Dagger][2][rep1T,rep2T,repAllT,\[Phi]1,\[Phi]2,\[Phi]3,h[2,5]/.Solh,h[3,5]/.Solh,h[3,7]/.Solh,h[3,9]/.Solh];*)
(*SolveBaxh\[Dagger][3][rep1T,rep2T,repAllT,\[Phi]1,\[Phi]2,\[Phi]3,h[2,5]/.Solh,h[3,5]/.Solh,h[3,7]/.Solh,h[3,9]/.Solh];*)*)
(**)
(**)
(*Print["Hi"];*)
(*Sol\[Dagger]=SolveBaxAllh\[Dagger][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,h[2,4]/.Solh,h[3,3]/.Solh,h[3,4]/.Solh,h[3,6]/.Solh,*)
(*h[4,3]/.Solh,h[4,4]/.Solh,h[4,5]/.Solh,h[4,6]/.Solh,h[4,8]/.Solh];*)
(**)
(*QH[1]=Q\[Dagger][1][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Sol\[Dagger];*)
(*QH[2]=Q\[Dagger][2][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Sol\[Dagger];*)
(*QH[3]=Q\[Dagger][3][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Sol\[Dagger];*)
(*QH[4]=Q\[Dagger][4][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Sol\[Dagger];*)
(**)
(*(*Print["Hodge"]*)
(*Hodge[{Qh[1],Qh[2],Qh[3]},rep1T,rep2T,\[Phi]1,\[Phi]2,\[Phi]3,z,2,1,1]//Simplify*)
(**)
(*Print["Q^\[Dagger]"]*)
(*QH[2]//Simplify*)*)
(**)
(*Table[(1/QH[t] Hodge[{Qh[1],Qh[2],Qh[3],Qh[4]},rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,z,t])//Simplify,{t,1,4}]*)
(**)
(**)
(**)
(*Print["Hello!"];*)
(*Q3[{Qh[1],Qh[2],Qh[3],Qh[4]},z,1,2,3];*)
(*Hodge[{Qh[1],Qh[2],Qh[3],Qh[4]},rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,z,4]*)
(*QH[4]//Simplify*)
(**)
(*QH[1]//Simplify*)
(*Wronskian[{Qh[2],Qh[3],Qh[4]},z]*Exp[(\[Phi]1T+\[Phi]2T+\[Phi]3T+\[Phi]4T)*z]//Simplify*)
(**)
(*QHT={QH[1],QH[2],QH[3],QH[4]};*)
(*QhT={Qh[1],Qh[2],Qh[3],Qh[4]};*)
(*Qh2T={Qh2[1],Qh2[2],Qh2[3],Qh2[4]};*)
(*OrthogRen[QHT,QhT,QHT,QhT,rep1T,rep2T]*)
(*OrthogRen[QHT,Qh2T,QHT,Qh2T,rep1T,rep2T]*)


(* ::Chapter::Closed:: *)
(*2-sites Scalar product*)


(* ::Section::Closed:: *)
(*Renormalised scalar product*)


(* ::Text:: *)
(*Input:*)
(*Q={Subscript[Q, 1],Subscript[Q, 2],Subscript[Q, 3]}*)
(*QH={Subscript[QH, 1],Subscript[QH, 2],Subscript[QH, 3]}*)


(* ::Input:: *)
(*ClearAll[Normalisation];*)
(*Normalisation[QHL_,QL_,QHR_,QR_,rep1_,rep2_]:=Normalisation[QHL,QL,QHR,QR,rep1,rep2]=(Limit[QHL[[1]],z->0,Direction->"FromAbove"])*(Limit[Wronskian[{QHL[[1]],QHL[[2]]},z]*1/(z^(rep1[[3]]-rep1[[4]])*(z-1)^(rep2[[3]]-rep2[[4]])),z->0,Direction->"FromAbove"])*(Limit[QL[[4]],z->0,Direction->"FromAbove"])*(Limit[QR[[1]],z->0,Direction->"FromAbove"])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z]*1/(z^(rep1[[1]]-rep1[[2]])*(z-1)^(rep2[[1]]-rep2[[2]])),z->0,Direction->"FromAbove"])*(Limit[QHR[[4]],z->0,Direction->"FromAbove"])**)
(*(Limit[QHL[[1]],z->1,Direction->"FromAbove"])*(Limit[Wronskian[{QHL[[1]],QHL[[2]]},z]*1/(z^(rep1[[3]]-rep1[[4]])*(z-1)^(rep2[[3]]-rep2[[4]])),z->1,Direction->"FromAbove"])*(Limit[QL[[4]],z->1,Direction->"FromAbove"])*(Limit[QR[[1]],z->1,Direction->"FromAbove"])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z]*1/(z^(rep1[[1]]-rep1[[2]])*(z-1)^(rep2[[1]]-rep2[[2]])),z->1,Direction->"FromAbove"])*(Limit[QHR[[4]],z->1,Direction->"FromAbove"]);*)


(* ::Input:: *)
(*ClearAll[OrthogRen,RawOrthogT,RawOrthog12];*)
(*RawOrthogT[QH_,Q_,rep1_,rep2_,w_,\[Alpha]_,i_,j_]:=RawOrthogT[QH,Q,rep1,rep2,w,\[Alpha],i,j]={MeanValue[QH[[i]]*D[Q[[j]],{z,2}]*1/(z-w),\[Alpha],rep1,rep2],MeanValue[QH[[i]]*D[Q[[j]],{z,1}]*1/(z-w)^2,\[Alpha],rep1,rep2],MeanValue[QH[[i]]*D[Q[[j]],{z,1}]*1/(z-w),\[Alpha],rep1,rep2],MeanValue[QH[[i]]*Q[[j]]*1/(z-w)^3,\[Alpha],rep1,rep2],MeanValue[QH[[i]]*Q[[j]]*1/(z-w)^2,\[Alpha],rep1,rep2],MeanValue[QH[[i]]*Q[[j]]*1/(z-w),\[Alpha],rep1,rep2]};*)
(**)
(*RawOrthog12[QH_,Q_,rep1_,rep2_,\[Alpha]_,i_,j_]:=RawOrthog12[QH,Q,rep1,rep2,\[Alpha],i,j]=Join[RawOrthogT[QH,Q,rep1,rep2,0,\[Alpha],i,j],RawOrthogT[QH,Q,rep1,rep2,1,\[Alpha],i,j]];*)
(**)
(*OrthogRen[QHL_,QL_,QHR_,QR_,rep1_,rep2_]:=OrthogRen[QHL,QL,QHR,QR,rep1,rep2]=(1/Normalisation[QHL,QL,QHR,QR,rep1,rep2])*Det[{RawOrthog12[QHL,QR,rep1,rep2,0,1,1],RawOrthog12[QHL,QR,rep1,rep2,0,1,2],RawOrthog12[QHL,QR,rep1,rep2,0,2,1],RawOrthog12[QHL,QR,rep1,rep2,0,2,2],RawOrthog12[QHL,QR,rep1,rep2,0,1,3],RawOrthog12[QHL,QR,rep1,rep2,0,3,1],RawOrthog12[QHL,QR,rep1,rep2,1,1,1],RawOrthog12[QHL,QR,rep1,rep2,1,1,2],RawOrthog12[QHL,QR,rep1,rep2,1,2,1],RawOrthog12[QHL,QR,rep1,rep2,1,2,2],RawOrthog12[QHL,QR,rep1,rep2,1,1,3],RawOrthog12[QHL,QR,rep1,rep2,1,3,1]}];*)


(* ::Subsubsection::Closed:: *)
(*TEST*)


(* ::Input:: *)
(*ClearAll[rep1T,rep2T,repAllT,M,M2,Solh,Solh2,Qh,Qh2,QH]*)
(*rep1T=ADJOINT;*)
(*rep2T=ADJOINT;*)
(*repAllT[i_]:=FindGoodWeight2[rep1T,rep2T][[i]];*)
(*\[Phi]1T=1;*)
(*\[Phi]2T=1/2;*)
(*\[Phi]3T=1/5;*)
(*\[Phi]4T=1/7;*)
(*M=1;*)
(*M2=3;*)
(*Solh=(SolveBaxAll[rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,1,1]//Simplify)[[1]];*)
(*Solh2=(SolveBaxAll[rep1T,rep2T,repAllT[M2],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,1,1]//Simplify)[[1]];*)
(**)
(*Qh[1]=Q[1][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh;*)
(*Qh[2]=Q[2][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh;*)
(*Qh[3]=Q[3][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh;*)
(*Qh[4]=Q[4][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh;*)
(**)
(*Qh2[1]=Q[1][rep1T,rep2T,repAllT[M2],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh2;*)
(*Qh2[2]=Q[2][rep1T,rep2T,repAllT[M2],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh2;*)
(*Qh2[3]=Q[3][rep1T,rep2T,repAllT[M2],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh2;*)
(*Qh2[4]=Q[4][rep1T,rep2T,repAllT[M2],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh2;*)
(**)
(*(*SolveBaxh\[Dagger][1][rep1T,rep2T,repAllT,\[Phi]1,\[Phi]2,\[Phi]3,h[2,5]/.Solh,h[3,5]/.Solh,h[3,7]/.Solh,h[3,9]/.Solh];*)
(*SolveBaxh\[Dagger][2][rep1T,rep2T,repAllT,\[Phi]1,\[Phi]2,\[Phi]3,h[2,5]/.Solh,h[3,5]/.Solh,h[3,7]/.Solh,h[3,9]/.Solh];*)
(*SolveBaxh\[Dagger][3][rep1T,rep2T,repAllT,\[Phi]1,\[Phi]2,\[Phi]3,h[2,5]/.Solh,h[3,5]/.Solh,h[3,7]/.Solh,h[3,9]/.Solh];*)*)
(**)
(**)
(*Print["Hi"];*)
(*Sol\[Dagger]=SolveBaxAllh\[Dagger][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,h[2,4]/.Solh,h[3,3]/.Solh,h[3,4]/.Solh,h[3,6]/.Solh,*)
(*h[4,3]/.Solh,h[4,4]/.Solh,h[4,5]/.Solh,h[4,6]/.Solh,h[4,8]/.Solh];*)
(**)
(*QH[1]=Q\[Dagger][1][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Sol\[Dagger];*)
(*QH[2]=Q\[Dagger][2][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Sol\[Dagger];*)
(*QH[3]=Q\[Dagger][3][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Sol\[Dagger];*)
(*QH[4]=Q\[Dagger][4][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Sol\[Dagger];*)
(**)
(*(*Print["Hodge"]*)
(*Hodge[{Qh[1],Qh[2],Qh[3]},rep1T,rep2T,\[Phi]1,\[Phi]2,\[Phi]3,z,2,1,1]//Simplify*)
(**)
(*Print["Q^\[Dagger]"]*)
(*QH[2]//Simplify*)*)
(**)
(*Table[(1/QH[t] Hodge[{Qh[1],Qh[2],Qh[3],Qh[4]},rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,z,t])//Simplify,{t,1,4}];*)
(**)
(**)
(*OrthogRen[{QH[1],QH[2],QH[3],QH[4]},{Qh[1],Qh[2],Qh[3],Qh[4]},{QH[1],QH[2],QH[3],QH[4]},{Qh[1],Qh[2],Qh[3],Qh[4]},rep1T,rep2T]*)
(*OrthogRen[{QH[1],QH[2],QH[3],QH[4]},{Qh2[1],Qh2[2],Qh2[3],Qh2[4]},{QH[1],QH[2],QH[3],QH[4]},{Qh2[1],Qh2[2],Qh2[3],Qh2[4]},rep1T,rep2T]*)
(**)
(*(**)
(*Print["Normalisation"];*)
(*Normalisation[{QH[1],QH[2],QH[3],QH[4]},{Qh[1],Qh[2],Qh[3],Qh[4]},{QH[1],QH[2],QH[3],QH[4]},{Qh[1],Qh[2],Qh[3],Qh[4]},rep1T,rep2T]*)
(*Normalisation[{QH[1],QH[2],QH[3],QH[4]},{Qh2[1],Qh2[2],Qh2[3],Qh2[4]},{QH[1],QH[2],QH[3],QH[4]},{Qh2[1],Qh2[2],Qh2[3],Qh2[4]},rep1T,rep2T]*)
(**)*)
(**)
(*(**)
(*Print["Norm0"]*)
(*(Limit[QH[1],z->0,Direction->"FromAbove"])*)
(*(Limit[Wronskian[{QH[1],QH[2]},z]*1/(z^(rep1T[[2]]-rep1T[[3]])*(z-1)^(rep2T[[2]]-rep2T[[3]])),z->0,Direction->"FromAbove"])*)
(*(Limit[Qh[1],z->0,Direction->"FromAbove"])*)
(*(Limit[Qh[1],z->0,Direction->"FromAbove"])*)
(*(Limit[Wronskian[{Qh[1],Qh[2]},z]*1/(z^(rep1T[[1]]-rep1T[[2]])*(z-1)^(rep2T[[1]]-rep2T[[2]])),z->0,Direction->"FromAbove"])*)
(*(Limit[QH[4],z->0,Direction->"FromAbove"])*)
(**)
(*Print["Norm1"]*)
(*(Limit[QH[1],z->1,Direction->"FromAbove"])*)
(*(Limit[Wronskian[{QH[1],QH[2]},z]*1/(z^(rep1T[[2]]-rep1T[[3]])*(z-1)^(rep2T[[2]]-rep2T[[3]])),z->1,Direction->"FromAbove"])*)
(*(Limit[Qh[1],z->1,Direction->"FromAbove"])*)
(*(Limit[Qh[1],z->1,Direction->"FromAbove"])*)
(*(Limit[Wronskian[{Qh[1],Qh[2]},z]*1/(z^(rep1T[[1]]-rep1T[[2]])*(z-1)^(rep2T[[1]]-rep2T[[2]])),z->1,Direction->"FromAbove"])*)
(*(Limit[QH[4],z->1,Direction->"FromAbove"])*)
(**)*)
(**)
(*(**)
(*Print["Reg Wr Q"];*)
(*Wronskian[{Qh[1],Qh[2]},z]*)
(**)
(*Print["Reg Wr QH"];*)
(*Wronskian[{QH[1],QH[2]},z]*)
(**)*)


(* ::Chapter::Closed:: *)
(*1-site Scalar product*)


(* ::Section::Closed:: *)
(*Renormalised scalar product*)


(* ::Text:: *)
(*Input:*)
(*Q={Subscript[Q, 1],Subscript[Q, 2],Subscript[Q, 3]}*)
(*QH={Subscript[QH, 1],Subscript[QH, 2],Subscript[QH, 3]}*)


(* ::Subsubsection::Closed:: *)
(*rep1*)


(* ::Input:: *)
(*ClearAll[Normalisation1];*)
(*Normalisation1[QHL_,QL_,QHR_,QR_,rep1_]:=Normalisation1[QHL,QL,QHR,QR,rep1]=(Limit[QHL[[1]],z->0,Direction->"FromAbove"])*(Limit[Wronskian[{QHL[[1]],QHL[[2]]},z]*1/z^(rep1[[3]]-rep1[[4]]),z->0,Direction->"FromAbove"])*(Limit[QL[[4]],z->0,Direction->"FromAbove"])*(Limit[QR[[1]],z->0,Direction->"FromAbove"])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),z->0,Direction->"FromAbove"])*(Limit[QHR[[4]],z->0,Direction->"FromAbove"]);*)


(* ::Input:: *)
(*ClearAll[RawOrthog1];*)
(*RawOrthog1[QH_,Q_,rep1_,w_,\[Alpha]_,i_,j_]:=RawOrthog1[QH,Q,rep1,w,\[Alpha],i,j]={MeanValue1[QH[[i]]*D[Q[[j]],{z,2}]*1/(z-w),\[Alpha],rep1],MeanValue1[QH[[i]]*D[Q[[j]],{z,1}]*1/(z-w)^2,\[Alpha],rep1],MeanValue1[QH[[i]]*D[Q[[j]],{z,1}]*1/(z-w),\[Alpha],rep1],MeanValue1[QH[[i]]*Q[[j]]*1/(z-w)^3,\[Alpha],rep1],MeanValue1[QH[[i]]*Q[[j]]*1/(z-w)^2,\[Alpha],rep1],MeanValue1[QH[[i]]*Q[[j]]*1/(z-w),\[Alpha],rep1]};*)
(*ClearAll[OrthogRen1];*)
(*OrthogRen1[QHL_,QL_,QHR_,QR_,rep1_]:=OrthogRen1[QHL,QL,QHR,QR,rep1]=(1/Normalisation1[QHL,QL,QHR,QR,rep1])*Det[{RawOrthog1[QHL,QR,rep1,0,0,1,1],RawOrthog1[QHL,QR,rep1,0,0,1,2],RawOrthog1[QHL,QR,rep1,0,0,2,1],RawOrthog1[QHL,QR,rep1,0,0,2,2],RawOrthog1[QHL,QR,rep1,0,0,1,3],RawOrthog1[QHL,QR,rep1,0,0,3,1]}];*)


(* ::Subsubsection::Closed:: *)
(*rep2*)


(* ::Input:: *)
(*ClearAll[Normalisation2];*)
(*Normalisation2[QHL_,QL_,QHR_,QR_,rep2_]:=Normalisation2[QHL,QL,QHR,QR,rep2]=*)
(*(Limit[QHL[[1]],z->1,Direction->"FromAbove"])*(Limit[Wronskian[{QHL[[1]],QHL[[2]]},z]*1/(z-1)^(rep2[[3]]-rep2[[4]]),z->1,Direction->"FromAbove"])*(Limit[QL[[4]],z->1,Direction->"FromAbove"])*(Limit[QR[[1]],z->1,Direction->"FromAbove"])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z]*1/(z-1)^(rep2[[1]]-rep2[[2]]),z->1,Direction->"FromAbove"])*(Limit[QHR[[4]],z->1,Direction->"FromAbove"]);*)


(* ::Input:: *)
(*ClearAll[RawOrthog2];*)
(*RawOrthog2[QH_,Q_,rep2_,w_,\[Alpha]_,i_,j_]:=RawOrthog2[QH,Q,rep2,w,\[Alpha],i,j]={MeanValue2[QH[[i]]*D[Q[[j]],{z,2}]*1/(z-w),\[Alpha],rep2],MeanValue2[QH[[i]]*D[Q[[j]],{z,1}]*1/(z-w)^2,\[Alpha],rep2],MeanValue2[QH[[i]]*D[Q[[j]],{z,1}]*1/(z-w),\[Alpha],rep2],MeanValue2[QH[[i]]*Q[[j]]*1/(z-w)^3,\[Alpha],rep2],MeanValue2[QH[[i]]*Q[[j]]*1/(z-w)^2,\[Alpha],rep2],MeanValue2[QH[[i]]*Q[[j]]*1/(z-w),\[Alpha],rep2]};*)
(*ClearAll[OrthogRen2];*)
(*OrthogRen2[QHL_,QL_,QHR_,QR_,rep2_]:=OrthogRen2[QHL,QL,QHR,QR,rep2]=(1/Normalisation2[QHL,QL,QHR,QR,rep2])*Det[{RawOrthog2[QHL,QR,rep2,0,0,1,1],RawOrthog2[QHL,QR,rep2,0,0,1,2],RawOrthog2[QHL,QR,rep2,0,0,2,1],RawOrthog2[QHL,QR,rep2,0,0,2,2],RawOrthog2[QHL,QR,rep2,0,0,1,3],RawOrthog2[QHL,QR,rep2,0,0,3,1]}];*)


(* ::Title::Closed:: *)
(*RUN: \[CapitalPsi] space*)


(* ::Chapter::Closed:: *)
(*Definition Lax connection*)


(* ::Subsection::Closed:: *)
(*Lax connection*)


(* ::TextIndent:: *)
(*Lax connection Subscript[L, ij] with regulator \[Phi]*)


(* ::Input:: *)
(*ClearAll[LIndex];*)
(*LIndex[a_,b_,rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=1/(z-w1)*KroneckerProduct[EM[a,b,rep1],Id[rep2]]+1/(z-w2)*KroneckerProduct[Id[rep1],EM[a,b,rep2]]+KroneckerDelta[a,1]*KroneckerDelta[b,1]*\[Phi]1*KroneckerProduct[Id[rep1],Id[rep2]]+KroneckerDelta[a,2]*KroneckerDelta[b,2]*\[Phi]2*KroneckerProduct[Id[rep1],Id[rep2]]+KroneckerDelta[a,3]*KroneckerDelta[b,3]*\[Phi]3*KroneckerProduct[Id[rep1],Id[rep2]]+KroneckerDelta[a,4]*KroneckerDelta[b,4]*\[Phi]4*KroneckerProduct[Id[rep1],Id[rep2]];*)


(* ::Subsubsection::Closed:: *)
(*Example of use*)


(* ::Input:: *)
(*ClearAll[rep1,rep2];*)
(*rep1=FUNDAMENTAL;*)
(*rep2=FUNDAMENTAL;*)
(*LIndex[1,1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal*)


(* ::Chapter::Closed:: *)
(*Hamiltonians*)


(* ::Section::Closed:: *)
(*H1*)


(* ::Subsection::Closed:: *)
(*Regulated Hamiltonian tr(L)*)


(* ::Input:: *)
(*ClearAll[H1];*)
(*H1[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=H1[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{a,1,Ng}]//Simplify;*)


(* ::Subsection::Closed:: *)
(*Definition Residue of H*)


(* ::Input:: *)
(*Clear[ResH1];*)
(*ResH1[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_]:=ResH1[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w]=SeriesCoefficient[H1[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,-1}]*)


(* ::Subsection:: *)
(*Example of use*)


(* ::Input:: *)
(*H1[FUNDAMENTAL,FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal;*)


(* ::Section::Closed:: *)
(*H2*)


(* ::Subsection:: *)
(*Regulated Hamiltonian 3/2*tr(L')+trA(L^2)*)


(* ::Subsubsection::Closed:: *)
(*Antisymm. def*)


(* ::Input:: *)
(*ClearAll[trAL2];*)
(*trAL2[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trAL2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=1/2*(Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[b,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[b,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{a,1,Ng},{b,1,Ng}]);*)


(* ::Subsubsection::Closed:: *)
(*Def. H2*)


(* ::Input:: *)
(*ClearAll[H2];*)
(*H2[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=H2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=(trAL2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+3/2*Sum[D[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z],{a,1,Ng}])//Simplify*)


(* ::Subsection::Closed:: *)
(*Definition Residue of H*)


(* ::Input:: *)
(*Clear[ResH2];*)
(*ResH2[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_]:=ResH2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w]=SeriesCoefficient[H2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,-1}]*)


(* ::Input:: *)
(*Clear[ResH2O];*)
(*ResH2O[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_,o_]:=ResH2O[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w,o]=SeriesCoefficient[H2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,o}]*)


(* ::Subsection::Closed:: *)
(*Example of use*)


(* ::Input:: *)
(*H2[FUNDAMENTAL,FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal;*)


(* ::Section::Closed:: *)
(*H3*)


(* ::Subsection:: *)
(*H3 Unsimplified*)


(* ::TextIndent:: *)
(*Regulated Hamiltonian trA(L^3)+4/3 trA(L\[CenterDot]L')+1/2 trA(L'\[CenterDot]L)+tr(L'')*)


(* ::Subsubsection::Closed:: *)
(*Antisymm. def*)


(* ::Input:: *)
(*ClearAll[trAL3];*)
(*trAL3[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trAL3[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=1/6*Sum[LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{i1,1,Ng},{i2,1,Ng},{i3,1,Ng}];*)


(* ::Input:: *)
(*ClearAll[trALpL,trALLp];*)
(*trALLp[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trALLp[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=1/2*(Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex[b,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]-LIndex[b,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z],{a,1,Ng},{b,1,Ng}]);*)
(*trALpL[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trALpL[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=1/2*(Sum[D[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex[b,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-D[LIndex[b,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{a,1,Ng},{b,1,Ng}]);*)


(* ::Subsubsection::Closed:: *)
(*Def H3*)


(* ::Input:: *)
(*ClearAll[H3U];*)
(*H3U[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=H3U[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=trAL3[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+4/3*trALLp[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+2/3*trALpL[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+Sum[D[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{z,2}],{a,1,Ng}]*)


(* ::Subsection::Closed:: *)
(*H3 Simplified*)


(* ::Input:: *)
(*ClearAll[H3];*)
(*H3[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=H3[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=H3U[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*)


(* ::Subsection::Closed:: *)
(*Definition Residue of H*)


(* ::Input:: *)
(*Clear[ResH3];*)
(*ResH3[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_]:=ResH3[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w]=SeriesCoefficient[H3[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,-1}]*)


(* ::Input:: *)
(*Clear[ResH3O];*)
(*ResH3O[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_,o_]:=ResH3O[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w,o]=SeriesCoefficient[H3[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,o}]*)


(* ::Subsection::Closed:: *)
(*Example of use*)


(* ::Input:: *)
(*H3[FUNDAMENTAL,FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal//Simplify*)


(* ::Section::Closed:: *)
(*H4*)


(* ::Subsection:: *)
(*H4 Unsimplified*)


(* ::TextIndent:: *)
(*Regulated Hamiltonian trA(L^3)+4/3 trA(L\[CenterDot]L')+1/2 trA(L'\[CenterDot]L)+tr(L'')*)


(* ::Subsubsection::Closed:: *)
(*Antisymm. def*)


(* ::Input:: *)
(*ClearAll[trAL4];*)
(*trAL4[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trAL4[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=1/4!*Sum[LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i4,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i4,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i3,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i4,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i4,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i2,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i4,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i3,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i4,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i2,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i4,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i4,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i3,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i4,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i4,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i1,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i4,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i3,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i4,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i1,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i4,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i4,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i2,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i4,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i4,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i1,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i4,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i2,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i4,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i1,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i4,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i3,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i4,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i2,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i4,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i3,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i4,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i1,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i4,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i2,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i4,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i1,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{i1,1,Ng},{i2,1,Ng},{i3,1,Ng},{i4,1,Ng}];*)


(* ::Input:: *)
(*(*L^3*)*)
(*ClearAll[trALLLp,trALLpL,trALpLL];*)
(*trALLLp[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trALLLp[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=1/6*Sum[LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]-LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]-LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]+LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]+LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]-LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z],{i1,1,Ng},{i2,1,Ng},{i3,1,Ng}];*)
(**)
(*trALLpL[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trALLpL[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=1/6*Sum[LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{i1,1,Ng},{i2,1,Ng},{i3,1,Ng}];*)
(**)
(*trALpLL[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trALpLL[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=1/6*Sum[D[LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-D[LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-D[LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+D[LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+D[LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-D[LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{i1,1,Ng},{i2,1,Ng},{i3,1,Ng}];*)


(* ::Input:: *)
(*(*L^2*)*)
(*ClearAll[trALLpp,trALpLp,trALppL];*)
(*trALLpp[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trALLpp[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=1/2*(Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex[b,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{z,2}]-LIndex[b,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{z,2}],{a,1,Ng},{b,1,Ng}]);*)
(**)
(*trALpLp[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trALpLp[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=1/2*(Sum[D[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . D[LIndex[b,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]-D[LIndex[b,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . D[LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z],{a,1,Ng},{b,1,Ng}]);*)
(**)
(*trALppL[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trALppL[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=1/2*(Sum[D[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{z,2}] . LIndex[b,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-D[LIndex[b,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{z,2}] . LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{a,1,Ng},{b,1,Ng}]);*)


(* ::Subsubsection::Closed:: *)
(*Def H4*)


(* ::Input:: *)
(*ClearAll[H4U];*)
(*H4U[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=H4U[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=trAL4[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+3/4*trALLLp[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+1/2*trALLpL[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+1/4*trALpLL[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+*)
(*1/2*trALLpp[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+1/2*trALpLp[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+1/6*trALppL[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+*)
(*1/4*Sum[D[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{z,3}],{a,1,Ng}]*)


(* ::Subsection::Closed:: *)
(*H4 Simplified*)


(* ::Input:: *)
(*ClearAll[H4];*)
(*H4[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=H4[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=H4U[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*)


(* ::Subsection::Closed:: *)
(*Definition Residue of H*)


(* ::Input:: *)
(*Clear[ResH4];*)
(*ResH4[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_]:=ResH4[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w]=SeriesCoefficient[H4[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,-1}]*)


(* ::Input:: *)
(*Clear[ResH4O];*)
(*ResH4O[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_,o_]:=ResH4O[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w,o]=SeriesCoefficient[H4[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,o}]*)


(* ::Subsection:: *)
(*Example of use*)


(* ::Input:: *)
(*H4[FUNDAMENTAL,FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal//Simplify;*)


(* ::Section::Closed:: *)
(*Check eigenvalues*)


(* ::Input:: *)
(*UPR[rep1_,rep2_]:=UPR[rep1,rep2]=Table[{KroneckerDelta[j,1]},{j,1,dim[rep1]*dim[rep2]}];*)


(* ::Input:: *)
(*rep1=FUNDAMENTAL;*)
(*rep2=ADJOINT;*)


(* ::Subsection::Closed:: *)
(*H2*)


(* ::Input:: *)
(*H2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . UPR[rep1,rep2]*)


(* ::Input:: *)
(*SeriesCoefficient[((H2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . UPR[rep1,rep2])//Flatten)[[1]],{z,0,-2}]//Normal*)
(*SeriesCoefficient[((H2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . UPR[rep1,rep2])//Flatten)[[1]],{z,1,-2}]//Normal*)


(* ::Input:: *)
(*rep1[[3]] (-1+rep1[[4]])+rep1[[2]] (-2+rep1[[3]]+rep1[[4]])+rep1[[1]] (-3+rep1[[2]]+rep1[[3]]+rep1[[4]])*)


(* ::Input:: *)
(*rep2[[3]] (-1+rep2[[4]])+rep2[[2]] (-2+rep2[[3]]+rep2[[4]])+rep2[[1]] (-3+rep2[[2]]+rep2[[3]]+rep2[[4]])*)


(* ::Subsection::Closed:: *)
(*H3*)


(* ::Input:: *)
(*SeriesCoefficient[((H3U[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . UPR[rep1,rep2])//Flatten)[[1]],{z,0,-3}]//Normal*)
(*SeriesCoefficient[((H3U[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . UPR[rep1,rep2])//Flatten)[[1]],{z,1,-3}]//Normal*)


(* ::Input:: *)
(*rep1[[2]] (-1+rep1[[3]]) (-2+rep1[[4]])+rep1[[1]] ((-2+rep1[[3]]) (-3+rep1[[4]])+rep1[[2]] (-4+rep1[[3]]+rep1[[4]]))*)


(* ::Input:: *)
(*rep2[[2]] (-1+rep2[[3]]) (-2+rep2[[4]])+rep2[[1]] ((-2+rep2[[3]]) (-3+rep2[[4]])+rep2[[2]] (-4+rep2[[3]]+rep2[[4]]))*)


(* ::Subsection::Closed:: *)
(*H4*)


(* ::Input:: *)
(*SeriesCoefficient[((H4U[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . UPR[rep1,rep2])//Flatten)[[1]],{z,0,-4}]//Normal*)
(*SeriesCoefficient[((H4U[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . UPR[rep1,rep2])//Flatten)[[1]],{z,1,-4}]//Normal*)


(* ::Input:: *)
(*rep1[[1]] (-1+rep1[[2]]) (-2+rep1[[3]]) (-3+rep1[[4]])*)
(*rep2[[1]] (-1+rep2[[2]]) (-2+rep2[[3]]) (-3+rep2[[4]])*)


(* ::Chapter::Closed:: *)
(*Adjoint Hamiltonians TO DO*)


(* ::Section::Closed:: *)
(*H1\[Dagger]*)


(* ::Subsection::Closed:: *)
(*Regulated Hamiltonian tr(L)*)


(* ::Input:: *)
(*ClearAll[H1\[Dagger]];*)
(*H1\[Dagger][rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_][z_]:=H1\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]=H1[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z](*Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}]//Simplify;*)*)


(* ::Subsection::Closed:: *)
(*Definition Residue of H*)


(* ::Input:: *)
(*Clear[ResH1\[Dagger]];*)
(*ResH1\[Dagger][rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,w_]:=ResH1\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,w]=SeriesCoefficient[H1\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]//Normal,{z,w,-1}]*)


(* ::Subsubsection::Closed:: *)
(*Example of use*)


(* ::Input:: *)
(*ClearAll[rep1,rep2];*)
(*rep1=FUNDAMENTAL;*)
(*rep2=FUNDAMENTAL;*)
(*H1\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]//Normal;*)


(* ::Section::Closed:: *)
(*H2\[Dagger]*)


(* ::Subsection::Closed:: *)
(*Regulated Hamiltonian -tr(L')+1/2 (tr(L))^2-1/2 tr(L^2)*)


(* ::Input:: *)
(*ClearAll[H2\[Dagger]];*)
(*H2\[Dagger][rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_][z_]:=H2\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]=(-Sum[D[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],z],{a,1,3}]+1/2*(Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}]) . (Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}])-1/2*Sum[LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z] . LIndex[b,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3},{b,1,3}])//Simplify*)


(* ::Subsection::Closed:: *)
(*Definition Residue of H*)


(* ::Input:: *)
(*Clear[ResH2\[Dagger]];*)
(*ResH2\[Dagger][rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,w_]:=ResH2\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,w]=SeriesCoefficient[H2\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]//Normal,{z,w,-1}]*)


(* ::Section::Closed:: *)
(*H3\[Dagger]*)


(* ::Subsection::Closed:: *)
(*H3\[Dagger] simplified 1/3 tr(L'')+tr(L' . L)-2/3 tr(L') . tr(L)-1/2 tr(L) . tr(L^2)+1/6 tr(L)^3+1/3 tr(L^3)*)


(* ::Input:: *)
(*ClearAll[H3\[Dagger]];*)
(*H3\[Dagger][rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_][z_]:=H3\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]=(1/3*Sum[D[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{z,2}],{a,1,3}]+Sum[D[LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],z] . LIndex[b,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3},{b,1,3}]-2/3*Sum[D[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],z],{a,1,3}] . Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}]-1/2*Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}] . Sum[LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z] . LIndex[b,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3},{b,1,3}]+1/6*(Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}]) . (Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}]) . (Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}])+1/3*Sum[LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z] . LIndex[b,c,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z] . LIndex[c,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3},{b,1,3},{c,1,3}])//Simplify;*)


(* ::Subsection::Closed:: *)
(*Definition Residue of H*)


(* ::Input:: *)
(*Clear[ResH3\[Dagger]];*)
(*ResH3\[Dagger][rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,w_]:=ResH3\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,w]=SeriesCoefficient[H3[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]//Normal,{z,w,-1}]*)


(* ::Subsection::Closed:: *)
(*Example of use*)


(* ::Input:: *)
(*Series[H3\[Dagger][ADJOINT,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]//Normal,{z,0,-3}]//Simplify//Normal;*)


(* ::Chapter::Closed:: *)
(*Eigenvectors (FIRST RUN Q SPACE)*)


(* ::Section::Closed:: *)
(*Right eigenvectors (column)*)


(* ::Input:: *)
(*ClearAll[VecR,VecR1];*)
(*VecR[rep1_,rep2_]:=VecR[rep1,rep2]=Table[{v[j]},{j,1,dim[rep1]*dim[rep2]}];*)
(*VecR1[rep1_]:=VecR1[rep1]=Table[{v[j]},{j,1,dim[rep1]}];*)


(* ::Input:: *)
(*ClearAll[EigenR];*)
(*EigenR[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=EigenR[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=Table[Quiet[Flatten[VecR[rep1,rep2]/.Solve[*)
(*CoefficientList[((H2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . VecR[rep1,rep2]-HSpectrum[2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z][[idx]]* VecR[rep1,rep2]))//Together//Numerator,z]==0&&*)
(*CoefficientList[((H3[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . VecR[rep1,rep2]-HSpectrum[3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z][[idx]]* VecR[rep1,rep2]))//Together//Numerator,z]==0&&*)
(*CoefficientList[((H4[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . VecR[rep1,rep2]-HSpectrum[4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z][[idx]]* VecR[rep1,rep2]))//Together//Numerator,z]==0,Flatten[VecR[rep1,rep2]]],1]],{idx,1,Dimensions[HSpectrum[1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]][[1]]}]*)


(* ::Subsubsection::Closed:: *)
(*Example*)


(* ::Input:: *)
(*EigenR[FUNDAMENTAL,FUNDAMENTAL,{1,1,0,0},3,5,7,9]*)


(* ::Section::Closed:: *)
(*Left eigenvectors (row)*)


(* ::Input:: *)
(*ClearAll[VecL,VecL1];*)
(*VecL[rep1_,rep2_]:=VecL[rep1,rep2]=Transpose[VecR[rep1,rep2]];*)
(*VecL1[rep1_]:=VecL1[rep1]=Transpose[VecR1[rep1]];*)


(* ::Input:: *)
(*ClearAll[EigenL];*)
(*EigenL[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=EigenL[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=Table[Quiet[Flatten[VecL[rep1,rep2]/.Solve[*)
(*CoefficientList[((VecL[rep1,rep2] . H2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-HSpectrum[2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z][[idx]]* VecL[rep1,rep2]))//Together//Numerator,z]==0&&*)
(*CoefficientList[((VecL[rep1,rep2] . H3[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-HSpectrum[3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z][[idx]]* VecL[rep1,rep2]))//Together//Numerator,z]==0&&*)
(*CoefficientList[((VecL[rep1,rep2] . H4[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-HSpectrum[4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z][[idx]]* VecL[rep1,rep2]))//Together//Numerator,z]==0,Flatten[VecL[rep1,rep2]]],1]],{idx,1,Dimensions[HSpectrum[1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]][[1]]}]*)


(* ::Subsubsection::Closed:: *)
(*Example*)


(* ::Input:: *)
(*EigenL[FUNDAMENTAL,FUNDAMENTAL,{2,0,0,0},3,5,7,9]*)


(* ::Section::Closed:: *)
(*Normalisation*)


(* ::Subsection::Closed:: *)
(*Highest weight state*)


(* ::Text:: *)
(*h.w.s. 2-sites*)


(* ::Input:: *)
(*ClearAll[UPR,UPL];*)
(*UPR[rep1_,rep2_]:=UPR[rep1,rep2]=Table[{KroneckerDelta[j,1]},{j,1,dim[rep1]*dim[rep2]}];*)
(*UPL[rep1_,rep2_]:=UPL[rep1,rep2]=Transpose[UPR[rep1,rep2]];*)


(* ::Text:: *)
(*h.w.s. 1-site*)


(* ::Input:: *)
(*ClearAll[UPR1,UPL1];*)
(*UPR1[rep1_]:=UPR1[rep1]=Table[{KroneckerDelta[j,1]},{j,1,dim[rep1]}];*)
(*UPL1[rep1_]:=UPL1[rep1]=Transpose[UPR1[rep1]];*)


(* ::Subsection::Closed:: *)
(*Vacuum vector*)


(* ::Text:: *)
(*Vacuum vector 2-sites*)


(* ::Input:: *)
(*ClearAll[VacR,VacL];*)
(*VacR[rep1_,rep2_,\[Beta]12_,\[Beta]13_,\[Beta]14_,\[Beta]23_,\[Beta]24_,\[Beta]34_]:=VacR[rep1,rep2,\[Beta]12,\[Beta]13,\[Beta]14,\[Beta]23,\[Beta]24,\[Beta]34]=*)
(*MatrixExp[\[Beta]24(EMAll[4,2,rep1,rep2])] . MatrixExp[\[Beta]23(EMAll[3,2,rep1,rep2])] . MatrixExp[\[Beta]14(EMAll[4,1,rep1,rep2])] . MatrixExp[\[Beta]13(EMAll[3,1,rep1,rep2])] . MatrixExp[\[Beta]34(EMAll[4,3,rep1,rep2])] . MatrixExp[\[Beta]12(EMAll[2,1,rep1,rep2])] . UPR[rep1,rep2];*)
(*VacL[rep1_,rep2_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=VacL[rep1,rep2,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=UPL[rep1,rep2] . MatrixExp[\[Alpha]12(EMAll[1,2,rep1,rep2])] . MatrixExp[\[Alpha]34(EMAll[3,4,rep1,rep2])] . MatrixExp[\[Alpha]13(EMAll[1,3,rep1,rep2])] . MatrixExp[\[Alpha]14(EMAll[1,4,rep1,rep2])] . MatrixExp[\[Alpha]23(EMAll[2,3,rep1,rep2])] . MatrixExp[\[Alpha]24(EMAll[2,4,rep1,rep2])];*)


(* ::Text:: *)
(*Vacuum vector 1-site*)


(* ::Input:: *)
(*ClearAll[VacR1,VacL1];*)
(*VacR1[rep1_,\[Beta]12_,\[Beta]13_,\[Beta]14_,\[Beta]23_,\[Beta]24_,\[Beta]34_]:=VacR1[rep1,\[Beta]23,\[Beta]12,\[Beta]13,\[Beta]14,\[Beta]23,\[Beta]24,\[Beta]34]=*)
(*MatrixExp[\[Beta]24(EM[4,2,rep1])] . MatrixExp[\[Beta]23(EM[3,2,rep1])] . MatrixExp[\[Beta]14(EM[4,1,rep1])] . MatrixExp[\[Beta]13(EM[3,1,rep1])] . MatrixExp[\[Beta]34(EM[4,3,rep1])] . MatrixExp[\[Beta]12(EM[2,1,rep1])] . UPR1[rep1];*)
(*VacL1[rep1_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=VacL1[rep1,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=UPL1[rep1] . MatrixExp[\[Alpha]12(EM[1,2,rep1])] . MatrixExp[\[Alpha]34(EM[3,4,rep1])] . MatrixExp[\[Alpha]13(EM[1,3,rep1])] . MatrixExp[\[Alpha]14(EM[1,4,rep1])] . MatrixExp[\[Alpha]23(EM[2,3,rep1])] . MatrixExp[\[Alpha]24(EM[2,4,rep1])];*)


(* ::Subsection::Closed:: *)
(*Renormalisation*)


(* ::Input:: *)
(*ClearAll[ReNormR,ReNormL,\[Alpha]23,\[Alpha]13,\[Alpha]12,\[Alpha]24,\[Alpha]14,\[Alpha]34,\[Beta]23,\[Beta]13,\[Beta]12,\[Beta]24,\[Beta]14,\[Beta]34];*)
(*ReNormR[rep1_,rep2_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_,\[Psi]_]:=ReNormR[rep1,rep2,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34,\[Psi]]=(\[Psi]*(1/(((VacL[rep1,rep2,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34] . \[Psi])//Flatten)[[1]])))//Simplify;*)
(*ReNormL[rep1_,rep2_,\[Beta]12_,\[Beta]13_,\[Beta]14_,\[Beta]23_,\[Beta]24_,\[Beta]34_,\[Psi]_]:=ReNormL[rep1,rep2,\[Beta]12,\[Beta]13,\[Beta]14,\[Beta]23,\[Beta]24,\[Beta]34,\[Psi]]=(\[Psi]*(1/((\[Psi] . VacR[rep1,rep2,\[Beta]12,\[Beta]13,\[Beta]14,\[Beta]23,\[Beta]24,\[Beta]34]//Flatten)[[1]])))//Simplify;*)


(* ::Subsubsection:: *)
(*Example of use*)


(* ::Input:: *)
(*ReNormR[FUNDAMENTAL,FUNDAMENTAL,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34,EigenR[FUNDAMENTAL,FUNDAMENTAL,{1,1,0,0},3,5,7,9][[1]]]*)


(* ::Input:: *)
(*ReNormL[FUNDAMENTAL,FUNDAMENTAL,\[Beta]12,\[Beta]13,\[Beta]14,\[Beta]23,\[Beta]24,\[Beta]34,EigenL[FUNDAMENTAL,FUNDAMENTAL,{1,1,0,0},3,5,7,9][[1]]]*)


(* ::Input:: *)
(*((ReNormL[FUNDAMENTAL,FUNDAMENTAL,\[Beta]12,\[Beta]13,\[Beta]14,\[Beta]23,\[Beta]24,\[Beta]34,EigenL[FUNDAMENTAL,FUNDAMENTAL,{1,1,0,0},3,5,7,9][[1]]] . ReNormR[FUNDAMENTAL,FUNDAMENTAL,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34,EigenR[FUNDAMENTAL,FUNDAMENTAL,{1,1,0,0},3,5,7,9][[1]]])//Flatten)[[1]]//FullSimplify*)


(* ::Subsection::Closed:: *)
(*Correct normalisation*)


(* ::Subsubsection::Closed:: *)
(*\[Beta] without fixing \[Alpha]*)


(* ::Input:: *)
(*ClearAll[\[Beta]rT,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];*)
(*\[Beta]rT[12][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=\[Beta]rT[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=(\[Alpha]23 (\[Phi]1-\[Phi]4))/(\[Alpha]13 (\[Phi]2-\[Phi]4));*)
(*\[Beta]rT[13][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=\[Beta]rT[13][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=((\[Phi]2-\[Phi]3) (\[Phi]1-\[Phi]4))/(\[Alpha]13 (\[Phi]1-\[Phi]2) (\[Phi]3-\[Phi]4));*)
(*\[Beta]rT[14][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=\[Beta]rT[14][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=-(((\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))/(\[Alpha]14 (\[Phi]1-\[Phi]2) (\[Phi]3-\[Phi]4)));*)
(*\[Beta]rT[23][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=\[Beta]rT[23][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=-(((\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))/(\[Alpha]23 (\[Phi]1-\[Phi]2) (\[Phi]3-\[Phi]4)));*)
(*\[Beta]rT[24][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=\[Beta]rT[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=(\[Alpha]13 (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))/(\[Alpha]14 \[Alpha]23 (\[Phi]1-\[Phi]2) (\[Phi]3-\[Phi]4));*)
(*\[Beta]rT[34][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=\[Beta]rT[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=(\[Alpha]13 (\[Phi]2-\[Phi]4))/(\[Alpha]14 (-\[Phi]2+\[Phi]3));*)


(* ::Subsubsection::Closed:: *)
(*fixing \[Alpha]*)


(* ::Input:: *)
(*ClearAll[\[Alpha]r,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];*)
(*\[Alpha]r[12][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=(\[Alpha]13 (-\[Phi]1+\[Phi]3))/(\[Alpha]23 (\[Phi]2-\[Phi]3));*)
(*\[Alpha]r[24][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=(\[Alpha]14 \[Alpha]23 (\[Phi]2-\[Phi]3) (\[Phi]1-\[Phi]4))/(\[Alpha]13 (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4));*)
(*\[Alpha]r[34][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=(\[Alpha]14 (\[Phi]1-\[Phi]4))/(\[Alpha]13 (\[Phi]1-\[Phi]3));*)


(* ::Subsubsection::Closed:: *)
(*\[Beta] with fixed \[Alpha]*)


(* ::Input:: *)
(*ClearAll[\[Beta]r,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];*)
(*\[Beta]r[12][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Beta]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=\[Beta]rT[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]];*)
(*\[Beta]r[13][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Beta]r[13][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=\[Beta]rT[13][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]];*)
(*\[Beta]r[14][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Beta]r[14][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=\[Beta]rT[14][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]];*)
(*\[Beta]r[23][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Beta]r[23][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=\[Beta]rT[23][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]];*)
(*\[Beta]r[24][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Beta]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=\[Beta]rT[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]];*)
(*\[Beta]r[34][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Beta]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=\[Beta]rT[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]];*)


(* ::Subsection::Closed:: *)
(*Renormalised eigenvectors*)


(* ::Input:: *)
(*ClearAll[ReEigenR,ReEigenL];*)
(*ReEigenR[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=ReEigenR[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=ReNormR[rep1,rep2,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],#]&/@EigenR[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];*)
(*ReEigenL[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=ReEigenL[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=ReNormL[rep1,rep2,\[Beta]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Beta]r[13][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Beta]r[14][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Beta]r[23][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Beta]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Beta]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],#]&/@EigenL[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];*)
(**)
(**)
(*(**)
(*ReEigenR[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]14_,\[Alpha]34_]:=ReEigenR[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34]=ReNormR[rep1,rep2,\[Alpha]r[23][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Alpha]r[13][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Alpha]12,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Alpha]14,\[Alpha]34,#]&/@EigenR[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];*)
(*ReEigenL[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]14_,\[Alpha]34_]:=ReEigenL[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34]=ReNormL[rep1,rep2,\[Beta]r[23][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Beta]r[13][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Beta]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Beta]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Beta]r[14][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Beta]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],#]&/@EigenL[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];*)
(**)*)
(*(**)
(*ReEigenR[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]14_,\[Alpha]34_]:=ReEigenR[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34]=ReNormR[rep1,rep2,\[Alpha]r[23][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Alpha]r[13][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Alpha]12,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Alpha]14,\[Alpha]34,#]&/@EigenR[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];*)
(*ReEigenL[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]14_,\[Alpha]34_]:=ReEigenL[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34]=ReNormL[rep1,rep2,\[Beta]r[23][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Beta]r[13][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Beta]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Beta]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Beta]r[14][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Beta]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],#]&/@EigenL[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];*)
(**)*)


(* ::Subsubsection::Closed:: *)
(*Example of use*)


(* ::Input:: *)
(*ReEigenL[FUNDAMENTAL,FUNDAMENTAL,{2,0,0,0},\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34]*)


(* ::Title::Closed:: *)
(*RUN: CREATE BASIS*)


(* ::Section::Closed:: *)
(*Q Base*)


(* ::Subsection::Closed:: *)
(*Q Base given repAll*)


(* ::Input:: *)
(*ClearAll[QBase,T];*)
(*QBase[T_][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=QBase[T][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=*)
(*Q[T][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]/.SolveBaxAll[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1];*)


(* ::Section::Closed:: *)
(*Q Hodge Base*)


(* ::Subsection::Closed:: *)
(*Q Hodge Base given repAll*)


(* ::Input:: *)
(*ClearAll[QHBase,T,U];*)
(*QHBase[T_][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=QHBase[T][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=Table[Hodge[{QBase[1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]],QBase[2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]],QBase[3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]],QBase[4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]]},rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,z,T],{U,1,Length[QBase[1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]]}];*)
(**)
(**)
(**)
(*(*Q\[Dagger][T][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]/.SolveBaxAllh\[Dagger][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1];*)*)
(**)
(*(*With Hodge*)*)
(*(*Table[Hodge[{QBase[1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]],QBase[2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]],QBase[3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]],QBase[4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]]},rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,z,T],{U,1,Length[QBase[1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]]}];*)*)


(* ::Section::Closed:: *)
(*Q Vacuum*)


(* ::Input:: *)
(*ClearAll[QVac2,QHVac2,QvQv2,QvQv2\[Omega]\[Phi]];*)
(*(*QvQv*)*)
(*(*VACUUM STATES*)*)
(*QVac2[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=QVac2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]={QBase[1][rep1,rep2,rep1+rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]],QBase[2][rep1,rep2,rep1+rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]],QBase[3][rep1,rep2,rep1+rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]],QBase[4][rep1,rep2,rep1+rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]]};*)
(*QHVac2[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=QHVac2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]={QHBase[1][rep1,rep2,rep1+rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]],QHBase[2][rep1,rep2,rep1+rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]],QHBase[3][rep1,rep2,rep1+rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]],QHBase[4][rep1,rep2,rep1+rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]]};*)
(**)
(*QvQv2[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=QvQv2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=OrthogRen[QHVac2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4],QVac2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4],QHVac2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4],QVac2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4],rep1,rep2];*)
(**)
(*QvQv2\[Omega]\[Phi][rep1_,rep2_,\[Omega]1_,\[Omega]2_,\[Omega]3_,\[Omega]4_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=QvQv2\[Omega]\[Phi][rep1,rep2,\[Omega]1,\[Omega]2,\[Omega]3,\[Omega]4,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=OrthogRen[QHVac2[rep1,rep2,\[Omega]1,\[Omega]2,\[Omega]3,\[Omega]4],QVac2[rep1,rep2,\[Omega]1,\[Omega]2,\[Omega]3,\[Omega]4],QHVac2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4],QVac2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4],rep1,rep2];*)


(* ::Title::Closed:: *)
(*RUN: Q space ONE SITE*)


(* ::Chapter::Closed:: *)
(*With regulators*)


(* ::Section::Closed:: *)
(*hs*)


(* ::Subsection::Closed:: *)
(*Cartan*)


(* ::Input:: *)
(*ClearAll[h1];*)
(*h1[rep1_][1,1]:=h1[rep1][1,1]=rep1[[1]]+rep1[[2]]+rep1[[3]]+rep1[[4]];*)
(*h1[rep1_][2,1]:=h1[rep1][2,1]=rep1[[3]] (-1+rep1[[4]])+rep1[[2]] (-2+rep1[[3]]+rep1[[4]])+rep1[[1]] (-3+rep1[[2]]+rep1[[3]]+rep1[[4]]);*)
(*h1[rep1_][3,1]:=h1[rep1][3,1]=rep1[[2]] (-1+rep1[[3]]) (-2+rep1[[4]])+rep1[[1]] ((-2+rep1[[3]]) (-3+rep1[[4]])+rep1[[2]] (-4+rep1[[3]]+rep1[[4]]));*)
(*h1[rep1_][4,1]:=h1[rep1][4,1]=rep1[[1]] (-1+rep1[[2]]) (-2+rep1[[3]]) (-3+rep1[[4]]);*)


(* ::Subsubsection::Closed:: *)
(*z^0 (dependence only on regulators)*)


(* ::Input:: *)
(*(*H1*)*)
(*h1[\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][1,3]:=h1[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][1,3]=\[Phi]1+\[Phi]2+\[Phi]3+\[Phi]4;*)
(**)
(*(*H2*)*)
(*h1[\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][2,5]:=h1[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][2,5]=\[Phi]1 \[Phi]2+\[Phi]1 \[Phi]3+\[Phi]2 \[Phi]3+\[Phi]1 \[Phi]4+\[Phi]2 \[Phi]4+\[Phi]3 \[Phi]4;*)
(**)
(*(*H3*)*)
(*h1[\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][3,7]:=h1[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][3,7]=\[Phi]1 \[Phi]2 \[Phi]3+\[Phi]1 \[Phi]2 \[Phi]4+\[Phi]1 \[Phi]3 \[Phi]4+\[Phi]2 \[Phi]3 \[Phi]4;*)
(**)
(*(*H4*)*)
(*h1[\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][4,9]:=h1[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][4,9]=\[Phi]1*\[Phi]2*\[Phi]3*\[Phi]4;*)


(* ::Subsection::Closed:: *)
(*Relation (dependence on weight1)*)


(* ::Text:: *)
(*By expanding at z=\[Infinity] and taking the coefficient of 1/z*)


(* ::Input:: *)
(*(*H2*)*)
(*h1[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][2,3]:=h1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][2,3]=(\[Phi]2+\[Phi]3+\[Phi]4) (weight1[[1]])+(\[Phi]1+\[Phi]3+\[Phi]4) (weight1[[2]])+(\[Phi]1+\[Phi]2+\[Phi]4) (weight1[[3]])+(\[Phi]1+\[Phi]2+\[Phi]3) (weight1[[4]]);*)
(**)
(*(*H3*)*)
(*h1[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][3,5]:=h1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][3,5]=\[Phi]3 \[Phi]4 (weight1[[1]]+weight1[[2]])+\[Phi]2 \[Phi]4 (weight1[[1]]+weight1[[3]])+\[Phi]1 \[Phi]4 (weight1[[2]]+weight1[[3]])+\[Phi]2 \[Phi]3 (weight1[[1]]+weight1[[4]])+\[Phi]1 \[Phi]3 (weight1[[2]]+weight1[[4]])+\[Phi]1 \[Phi]2 (weight1[[3]]+weight1[[4]]);*)
(**)
(*(*H4*)*)
(*h1[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][4,7]:=h1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][4,7]=\[Phi]2 \[Phi]3 \[Phi]4 (weight1[[1]])+\[Phi]1 \[Phi]3 \[Phi]4 (weight1[[2]])+\[Phi]1 \[Phi]2 \[Phi]4 (weight1[[3]])+\[Phi]1 \[Phi]2 \[Phi]3 (weight1[[4]]);*)


(* ::Subsection::Closed:: *)
(*Hamiltonian in Q-space*)


(* ::Input:: *)
(*ClearAll[\[Epsilon]11,\[Epsilon]12,\[Epsilon]13];*)
(*\[Epsilon]11[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]11[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=h1[rep1][1,1]/(z-w1)+h1[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][1,3];*)
(*\[Epsilon]12[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]12[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=h1[rep1][2,1]/(z-w1)^2+h1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][2,3]/(z-w1)+h1[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][2,5];*)
(*\[Epsilon]13[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]13[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=h1[rep1][3,1]/(z-w1)^3+h1[3,3]/(z-w1)^2+h1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][3,5]/(z-w1)+h1[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][3,7];*)
(*\[Epsilon]14[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]14[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=h1[rep1][4,1]/(z-w1)^4+h1[4,3]/(z-w1)^3+h1[4,5]/(z-w1)^2+h1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][4,7]/(z-w1)+h1[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][4,9];*)


(* ::Section::Closed:: *)
(*Baxter equation with F*)


(* ::Input:: *)
(*ClearAll[F1];*)
(*F1[rep1_,C1_][z_]:=F1[rep1,C1][z]=1/(z-w1)^(rep1[[C1]]-C1+1);*)


(* ::Text:: *)
(*Definition of shifted Hamiltonian*)


(* ::Input:: *)
(*ClearAll[\[Epsilon]11S,\[Epsilon]12S,\[Epsilon]13S]*)
(*\[Epsilon]11S[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_][z_]:=\[Epsilon]11S[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1][z]=\[Epsilon]11[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+(4 Derivative[1][F1[rep1,C1]][z])/F1[rep1,C1][z];*)
(*\[Epsilon]12S[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_][z_]:=\[Epsilon]12S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1][z]=\[Epsilon]12[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+(3 \[Epsilon]11[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] Derivative[1][F1[rep1,C1]][z])/F1[rep1,C1][z]+(6 (F1[rep1,C1]^\[Prime]\[Prime])[z])/F1[rep1,C1][z];*)
(*\[Epsilon]13S[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_][z_]:=\[Epsilon]13S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1][z]=\[Epsilon]13[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+2*\[Epsilon]12[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*Derivative[1][F1[rep1,C1]][z]/F1[rep1,C1][z]+3*(\[Epsilon]11[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] (F1[rep1,C1]^\[Prime]\[Prime])[z])/F1[rep1,C1][z]+4*\!\(\*SuperscriptBox[\(F1[rep1, C1]\), *)
(*TagBox[*)
(*RowBox[{"(", "3", ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[z]/F1[rep1,C1][z];*)
(*\[Epsilon]14S[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_][z_]:=\[Epsilon]14S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1][z]=\[Epsilon]14[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+\[Epsilon]13[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*Derivative[1][F1[rep1,C1]][z]/F1[rep1,C1][z]+(\[Epsilon]12[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] (F1[rep1,C1]^\[Prime]\[Prime])[z])/F1[rep1,C1][z]+\[Epsilon]11[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*\!\(\*SuperscriptBox[\(F1[rep1, C1]\), *)
(*TagBox[*)
(*RowBox[{"(", "3", ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[z]/F1[rep1,C1][z]+\!\(\*SuperscriptBox[\(F1[rep1, C1]\), *)
(*TagBox[*)
(*RowBox[{"(", "4", ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[z]/F1[rep1,C1][z];*)


(* ::Input:: *)
(*ClearAll[BaxOp1];*)
(*BaxOp1[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_]:=D[#,{z,4}]+\[Epsilon]11S[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1][z]*D[#,{z,3}]+\[Epsilon]12S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1][z]*D[#,{z,2}]+\[Epsilon]13S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1][z]*D[#,z]+\[Epsilon]14S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1][z]*#&;*)


(* ::Section::Closed:: *)
(*Q functions*)


(* ::Input:: *)
(*ClearAll[QPoly,Q1,c11,c12,c13,c14];*)
(*QPoly[c_,M_,r_][z_]:=QPoly[c,M,r][z]=Exp[r z]*(z^M+Sum[c[k]*z^k,{k,0,M-1}]);*)
(*Q1[1][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q1[1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly[c11,rep1[[1]]-weight1[[1]],-\[Phi]1][z];*)
(*Q1[2][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q1[2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly[c12,rep1[[1]]-weight1[[2]],-\[Phi]2][z];*)
(*Q1[3][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q1[3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly[c13,rep1[[1]]-weight1[[3]],-\[Phi]3][z];*)
(*Q1[4][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q1[4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly[c14,rep1[[1]]-weight1[[4]],-\[Phi]4][z];*)


(* ::Section::Closed:: *)
(*Baxter equations*)


(* ::Text:: *)
(*Baxter equations*)


(* ::Input:: *)
(*ClearAll[Bax1];*)
(*Bax1[1][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_]:=Bax1[1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1]=((BaxOp1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1][Q1[1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[\[Phi]1 z])//Simplify;*)
(*Bax1[2][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_]:=Bax1[2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1]=((BaxOp1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1][Q1[2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[\[Phi]2 z])//Simplify;*)
(*Bax1[3][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_]:=Bax1[3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1]=((BaxOp1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1][Q1[3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[\[Phi]3 z])//Simplify;*)
(*Bax1[4][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_]:=Bax1[4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1]=((BaxOp1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1][Q1[4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[\[Phi]3 z])//Simplify;*)


(* ::Text:: *)
(*Unknown coefficients*)


(* ::Input:: *)
(*ClearAll[Vars11,Vars12,Vars13,Vars14];*)
(*Vars11[rep1_,weight1_]:=Vars11[rep1,weight1]=Table[c11[n],{n,0,rep1[[1]]-weight1[[1]]-1}];*)
(*Vars12[rep1_,weight1_]:=Vars12[rep1,weight1]=Table[c12[n],{n,0,rep1[[1]]-weight1[[2]]-1}];*)
(*Vars13[rep1_,weight1_]:=Vars13[rep1,weight1]=Table[c13[n],{n,0,rep1[[1]]-weight1[[3]]-1}];*)
(*Vars14[rep1_,weight1_]:=Vars14[rep1,weight1]=Table[c14[n],{n,0,rep1[[1]]-weight1[[4]]-1}];*)


(* ::Text:: *)
(*Solve Baxter equation*)


(* ::Input:: *)
(*ClearAll[SolveBaxAll1];*)
(*SolveBaxAll1[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,c1_]:=SolveBaxAll1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,c1]=Solve[CoefficientList[Bax1[1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,c1]//Together//Numerator,z]==0&&*)
(*CoefficientList[Bax1[2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,c1]//Together//Numerator,z]==0&&CoefficientList[Bax1[3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,c1]//Together//Numerator,z]==0&&CoefficientList[Bax1[4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,c1]//Together//Numerator,z]==0,Join[Vars11[rep1,weight1],Vars12[rep1,weight1],Vars13[rep1,weight1],Vars14[rep1,weight1],{h1[3,3],h1[4,3],h1[4,5]}]]//Simplify;*)


(* ::Text:: *)
(*Spectrum of Hamiltonian*)


(* ::Input:: *)
(*ClearAll[HSpectrum1];*)
(*HSpectrum1[1][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=HSpectrum1[1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]11[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]/.SolveBaxAll1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1];*)
(**)
(*HSpectrum1[2][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=HSpectrum1[2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]12[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]/.SolveBaxAll1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1];*)
(**)
(*HSpectrum1[3][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=HSpectrum1[3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]13[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]/.SolveBaxAll1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1];*)
(**)
(*HSpectrum1[4][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=HSpectrum1[4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]14[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]/.SolveBaxAll1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1];*)


(* ::Text:: *)
(*Shifted Spectrum of Hamiltonian*)


(* ::Input:: *)
(*ClearAll[HSpectrum1S];*)
(*HSpectrum1S[1][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=HSpectrum1S[1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]11S[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]/.SolveBaxAll1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1];*)
(**)
(*HSpectrum1S[2][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=HSpectrum1S[2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]12S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]/.SolveBaxAll1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1];*)
(**)
(*HSpectrum1S[3][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=HSpectrum1S[3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]13S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]/.SolveBaxAll1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1];*)
(**)
(*HSpectrum1S[4][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=HSpectrum1S[4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]14S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]/.SolveBaxAll1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1];*)


(* ::Subsubsection:: *)
(*Example of use*)


(* ::Input:: *)
(*Sol=SolveBaxAll1[FUNDAMENTAL,FUNDAMENTAL,1,1/3,1/5,1/7,1]*)


(* ::Input:: *)
(*H1ST=HSpectrum1[1][FUNDAMENTAL,FUNDAMENTAL,1,1/3,1/5,1/7][z]*)
(*H2ST=HSpectrum1[2][FUNDAMENTAL,FUNDAMENTAL,1,1/3,1/5,1/7][z]*)
(*H3ST=HSpectrum1[3][FUNDAMENTAL,FUNDAMENTAL,1,1/3,1/5,1/7][z]*)
(*H4ST=HSpectrum1[4][FUNDAMENTAL,FUNDAMENTAL,1,1/3,1/5,1/7][z]*)


(* ::Chapter::Closed:: *)
(*Mean value*)


(* ::Section::Closed:: *)
(*Measure*)


(* ::Subsection::Closed:: *)
(*1-site*)


(* ::Text:: *)
(*Measure such that <f Subscript[O, F]g> \!\(\*OverscriptBox[\(=\), \(!\)]\) <g SuperDagger[Subscript[O, F]] f>*)


(* ::Input:: *)
(*ClearAll[\[Mu]1,\[Mu]2];*)
(*\[Mu]1[rep1_][z_]:=\[Mu]1[rep1][z]=z^(-rep1[[1]]+rep1[[4]]);*)
(*\[Mu]2[rep2_][z_]:=\[Mu]2[rep2][z]=(-1+z)^(-rep2[[1]]+rep2[[4]]);*)


(* ::Text:: *)
(*Mean value*)


(* ::Input:: *)
(*ClearAll[MeanValue1,MeanValue2];*)
(*MeanValue1[f_,\[Tau]_,rep1_]:=MeanValue1[f,\[Tau],rep1]=*)
(*2 \[Pi] I Residue[(f*\[Mu]1[rep1][z]),{z,\[Tau]}];*)
(*MeanValue2[f_,\[Tau]_,rep2_]:=MeanValue2[f,\[Tau],rep2]=*)
(*2 \[Pi] I Residue[(f*\[Mu]2[rep2][z]),{z,\[Tau]}];*)


(* ::Chapter::Closed:: *)
(*ADJOINT BAXTER*)


(* ::Section::Closed:: *)
(*Adjoint Baxter equation with SuperDagger[F]*)


(* ::Text:: *)
(*HERE we consider only the case C1=1,C2=1 in the definition of M*)


(* ::Text:: *)
(*Shifted Hamiltonians (already with SuperDagger[F])*)


(* ::Input:: *)
(*ClearAll[\[Epsilon]11S\[Dagger],\[Epsilon]12S\[Dagger],\[Epsilon]13S\[Dagger]];*)
(*\[Epsilon]11S\[Dagger][rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]11S\[Dagger][rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=-\[Epsilon]11S[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]+4*D[\[Mu]1[rep1][z],z]/\[Mu]1[rep1][z];*)
(*\[Epsilon]12S\[Dagger][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]12S\[Dagger][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]12S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]-3*D[\[Epsilon]11S[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z],z]-3*\[Epsilon]11S[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]*D[\[Mu]1[rep1][z],z]/\[Mu]1[rep1][z]+6*D[\[Mu]1[rep1][z],{z,2}]/\[Mu]1[rep1][z];*)
(*\[Epsilon]13S\[Dagger][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]13S\[Dagger][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=-\[Epsilon]13S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]+2*D[\[Epsilon]12S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z],z]-3*D[\[Epsilon]11S[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z],{z,2}]+D[\[Mu]1[rep1][z],z]/\[Mu]1[rep1][z]*(2*\[Epsilon]12S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]-6*D[\[Epsilon]11S[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z],z])-3*D[\[Mu]1[rep1][z],{z,2}]/\[Mu]1[rep1][z]*\[Epsilon]11S[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]+4*D[\[Mu]1[rep1][z],{z,3}]/\[Mu]1[rep1][z];*)
(*\[Epsilon]14S\[Dagger][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]14S\[Dagger][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]14S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]-D[\[Epsilon]13S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z],z]+D[\[Epsilon]12S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z],{z,2}]-D[\[Epsilon]11S[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z],{z,3}]+D[\[Mu]1[rep1][z],z]/\[Mu]1[rep1][z]*(-\[Epsilon]13S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]+2*D[\[Epsilon]12S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z],z]-3*D[\[Epsilon]11S[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z],{z,2}])+D[\[Mu]1[rep1][z],{z,2}]/\[Mu]1[rep1][z]*(\[Epsilon]12S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]-3D[\[Epsilon]11S[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z],z])-\[Epsilon]11S[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]*D[\[Mu]1[rep1][z],{z,3}]/\[Mu]1[rep1][z]+D[\[Mu]1[rep1][z],{z,4}]/\[Mu]1[rep1][z];*)


(* ::Input:: *)
(*ClearAll[BaxOp1\[Dagger]];*)
(*BaxOp1\[Dagger][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=D[#,{z,4}]+\[Epsilon]11S\[Dagger][rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*D[#,{z,3}]+\[Epsilon]12S\[Dagger][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*D[#,{z,2}]+\[Epsilon]13S\[Dagger][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*D[#,z]+\[Epsilon]14S\[Dagger][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*#&;*)


(* ::Section::Closed:: *)
(*Q functions*)


(* ::Input:: *)
(*ClearAll[QPoly\[Dagger],Q1\[Dagger],c11\[Dagger],c12\[Dagger],c13\[Dagger],c14\[Dagger]];*)
(*QPoly\[Dagger][c_,M_,r_][z_]:=QPoly\[Dagger][c,M,r][z]=Exp[r z]*(z^M+Sum[c[k]*z^k,{k,0,M-1}]);*)
(*Q1\[Dagger][1][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q1\[Dagger][1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly\[Dagger][c11\[Dagger],-rep1[[4]]+weight1[[1]],\[Phi]1][z];*)
(*Q1\[Dagger][2][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q1\[Dagger][2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly\[Dagger][c12\[Dagger],-rep1[[4]]+weight1[[2]],\[Phi]2][z];*)
(*Q1\[Dagger][3][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q1\[Dagger][3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly\[Dagger][c13\[Dagger],-rep1[[4]]+weight1[[3]],\[Phi]3][z];*)
(*Q1\[Dagger][4][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q1\[Dagger][4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly\[Dagger][c14\[Dagger],-rep1[[4]]+weight1[[4]],\[Phi]4][z];*)


(* ::Section::Closed:: *)
(*Baxter equations*)


(* ::Input:: *)
(*ClearAll[Bax1\[Dagger]];*)
(*Bax1\[Dagger][1][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=Bax1\[Dagger][1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=((BaxOp1\[Dagger][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][Q1\[Dagger][1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[-\[Phi]1 z])//Simplify;*)
(*Bax1\[Dagger][2][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=Bax1\[Dagger][2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=((BaxOp1\[Dagger][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][Q1\[Dagger][2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[-\[Phi]2 z])//Simplify;*)
(*Bax1\[Dagger][3][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=Bax1\[Dagger][3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=((BaxOp1\[Dagger][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][Q1\[Dagger][3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[-\[Phi]3 z])//Simplify;*)
(*Bax1\[Dagger][4][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=Bax1\[Dagger][4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=((BaxOp1\[Dagger][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][Q1\[Dagger][4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[-\[Phi]4 z])//Simplify;*)


(* ::Text:: *)
(*Unknown coefficients*)


(* ::Input:: *)
(*ClearAll[Vars11\[Dagger],Vars12\[Dagger],Vars13\[Dagger],Vars14\[Dagger]];*)
(*Vars11\[Dagger][rep1_,weight1_]:=Vars11\[Dagger][rep1,weight1]=Table[c11\[Dagger][n],{n,0,-rep1[[4]]+weight1[[1]]-1}];*)
(*Vars12\[Dagger][rep1_,weight1_]:=Vars12\[Dagger][rep1,weight1]=Table[c12\[Dagger][n],{n,0,-rep1[[4]]+weight1[[2]]-1}];*)
(*Vars13\[Dagger][rep1_,weight1_]:=Vars13\[Dagger][rep1,weight1]=Table[c13\[Dagger][n],{n,0,-rep1[[4]]+weight1[[3]]-1}];*)
(*Vars14\[Dagger][rep1_,weight1_]:=Vars14\[Dagger][rep1,weight1]=Table[c14\[Dagger][n],{n,0,-rep1[[4]]+weight1[[4]]-1}];*)


(* ::Subsection::Closed:: *)
(*Solving Baxter equations giving hs*)


(* ::Input:: *)
(*ClearAll[SolveBaxh1\[Dagger],SolveBaxAllh1\[Dagger]];*)
(*SolveBaxh1\[Dagger][1][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,h33_,h43_,h45_]:=SolveBaxh1\[Dagger][1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h33,h43,h45]=Solve[CoefficientList[Bax1\[Dagger][1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]/.{h1[3,3]->h33,*)
(*h1[4,3]->h43,h1[4,5]->h45}//Together//Numerator,z]==0,Vars11\[Dagger][rep1,weight1]]//Simplify//Flatten;*)
(**)
(*SolveBaxh1\[Dagger][2][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,h33_,h43_,h45_]:=SolveBaxh1\[Dagger][2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h33,h43,h45]=Solve[CoefficientList[Bax1\[Dagger][2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]/.{h1[3,3]->h33,*)
(*h1[4,3]->h43,h1[4,5]->h45}//Together//Numerator,z]==0,Vars12\[Dagger][rep1,weight1]]//Simplify//Flatten;*)
(**)
(*SolveBaxh1\[Dagger][3][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,h33_,h43_,h45_]:=SolveBaxh1\[Dagger][3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h33,h43,h45]=Solve[CoefficientList[Bax1\[Dagger][3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]/.{h1[3,3]->h33,*)
(*h1[4,3]->h43,h1[4,5]->h45}//Together//Numerator,z]==0,Vars13\[Dagger][rep1,weight1]]//Simplify//Flatten;*)
(**)
(*SolveBaxh1\[Dagger][4][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,h33_,h43_,h45_]:=SolveBaxh1\[Dagger][4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h33,h43,h45]=Solve[CoefficientList[Bax1\[Dagger][4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]/.{h1[3,3]->h33,*)
(*h1[4,3]->h43,h1[4,5]->h45}//Together//Numerator,z]==0,Vars14\[Dagger][rep1,weight1]]//Simplify//Flatten;*)
(**)
(*SolveBaxAllh1\[Dagger][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,h33_,h43_,h45_]:=SolveBaxAllh1\[Dagger][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h33,h43,h45]=Join[SolveBaxh1\[Dagger][1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h33,h43,h45],SolveBaxh1\[Dagger][2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h33,h43,h45],SolveBaxh1\[Dagger][3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h33,h43,h45],SolveBaxh1\[Dagger][4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h33,h43,h45]]//Simplify;*)


(* ::Subsubsection:: *)
(*Example of use*)


(* ::Input:: *)
(*Solh1=(SolveBaxAll1[FUNDAMENTAL,FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1]//Simplify)[[1]]*)
(**)
(*SolveBaxh1\[Dagger][1][FUNDAMENTAL,FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h1[3,3]/.Solh1,h1[4,3]/.Solh1,h1[4,5]/.Solh1]*)
(*SolveBaxh1\[Dagger][2][FUNDAMENTAL,FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h1[3,3]/.Solh1,h1[4,3]/.Solh1,h1[4,5]/.Solh1]*)
(*SolveBaxh1\[Dagger][3][FUNDAMENTAL,FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h1[3,3]/.Solh1,h1[4,3]/.Solh1,h1[4,5]/.Solh1]*)
(*SolveBaxh1\[Dagger][4][FUNDAMENTAL,FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h1[3,3]/.Solh1,h1[4,3]/.Solh1,h1[4,5]/.Solh1]*)
(**)
(*Print["Hi"]*)
(*SolveBaxAllh1\[Dagger][FUNDAMENTAL,FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h1[3,3]/.Solh1,h1[4,3]/.Solh1,h1[4,5]/.Solh1]*)


(* ::Section::Closed:: *)
(*Q Hodge*)


(* ::Input:: *)
(*ClearAll[Ng,SetN];*)
(*Ng=4;*)
(*SetN=Range[Ng];*)


(* ::Input:: *)
(*ClearAll[Hodge1];*)
(*Hodge1[Qfour_,rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,z_,i_]:=*)
(*(Exp[z*(\[Phi]1+\[Phi]2+\[Phi]3+\[Phi]4)]*Signature[Join[Complement[SetN,{i}],{i}]]*Wronskian[{Qfour[[Complement[SetN,{i}]]][[1]],Qfour[[Complement[SetN,{i}]]][[2]],Qfour[[Complement[SetN,{i}]]][[3]]},z]*(((z)^-(2*rep1[[1]]-rep1[[2]]-rep1[[3]])))/((({\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[1]]]]-{\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[2]]]]))*({\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[1]]]]-{\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[3]]]])*({\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[2]]]]-{\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[3]]]])))//Simplify;*)
(*(*The Q argument should be Q1[z], not Q1!*)*)
(**)
(*(**(((z)^(-2(rep1[[1]]-rep1[[2]]))*(z-1)^(-2(rep2[[1]]-rep2[[2]]))))*)*)
(*(**)*)


(* ::Chapter::Closed:: *)
(*1-site Scalar product*)


(* ::Section::Closed:: *)
(*Renormalised scalar product*)


(* ::Text:: *)
(*Input:*)
(*Q={Subscript[Q, 1],Subscript[Q, 2],Subscript[Q, 3]}*)
(*QH={Subscript[QH, 1],Subscript[QH, 2],Subscript[QH, 3]}*)


(* ::Subsubsection::Closed:: *)
(*rep1*)


(* ::Input:: *)
(*ClearAll[Normalisation1];*)
(*Normalisation1[QHL_,QL_,QHR_,QR_,rep1_]:=Normalisation1[QHL,QL,QHR,QR,rep1]=(Limit[QHL[[1]],z->0,Direction->"FromAbove"])*(Limit[Wronskian[{QHL[[1]],QHL[[2]]},z]*1/z^(rep1[[3]]-rep1[[4]]),z->0,Direction->"FromAbove"])*(Limit[QL[[4]],z->0,Direction->"FromAbove"])*(Limit[QR[[1]],z->0,Direction->"FromAbove"])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),z->0,Direction->"FromAbove"])*(Limit[QHR[[4]],z->0,Direction->"FromAbove"]);*)


(* ::Input:: *)
(*ClearAll[RawOrthog1];*)
(*RawOrthog1[QH_,Q_,rep1_,w_,\[Alpha]_,i_,j_]:=RawOrthog1[QH,Q,rep1,w,\[Alpha],i,j]={MeanValue1[QH[[i]]*D[Q[[j]],{z,2}]*1/(z-w),\[Alpha],rep1],MeanValue1[QH[[i]]*D[Q[[j]],{z,1}]*1/(z-w)^2,\[Alpha],rep1],MeanValue1[QH[[i]]*D[Q[[j]],{z,1}]*1/(z-w),\[Alpha],rep1],MeanValue1[QH[[i]]*Q[[j]]*1/(z-w)^3,\[Alpha],rep1],MeanValue1[QH[[i]]*Q[[j]]*1/(z-w)^2,\[Alpha],rep1],MeanValue1[QH[[i]]*Q[[j]]*1/(z-w),\[Alpha],rep1]};*)
(*ClearAll[OrthogRen1];*)
(*OrthogRen1[QHL_,QL_,QHR_,QR_,rep1_]:=OrthogRen1[QHL,QL,QHR,QR,rep1]=(1/Normalisation1[QHL,QL,QHR,QR,rep1])*Det[{RawOrthog1[QHL,QR,rep1,0,0,1,1],RawOrthog1[QHL,QR,rep1,0,0,1,2],RawOrthog1[QHL,QR,rep1,0,0,2,1],RawOrthog1[QHL,QR,rep1,0,0,2,2],RawOrthog1[QHL,QR,rep1,0,0,1,3],RawOrthog1[QHL,QR,rep1,0,0,3,1]}];*)


(* ::Subsubsection::Closed:: *)
(*rep2*)


(* ::Input:: *)
(*ClearAll[Normalisation2];*)
(*Normalisation2[QHL_,QL_,QHR_,QR_,rep2_]:=Normalisation2[QHL,QL,QHR,QR,rep2]=*)
(*(Limit[QHL[[1]],z->1,Direction->"FromAbove"])*(Limit[Wronskian[{QHL[[1]],QHL[[2]]},z]*1/(z-1)^(rep2[[3]]-rep2[[4]]),z->1,Direction->"FromAbove"])*(Limit[QL[[4]],z->1,Direction->"FromAbove"])*(Limit[QR[[1]],z->1,Direction->"FromAbove"])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z]*1/(z-1)^(rep2[[1]]-rep2[[2]]),z->1,Direction->"FromAbove"])*(Limit[QHR[[4]],z->1,Direction->"FromAbove"]);*)


(* ::Input:: *)
(*ClearAll[RawOrthog2];*)
(*RawOrthog2[QH_,Q_,rep2_,w_,\[Alpha]_,i_,j_]:=RawOrthog2[QH,Q,rep2,w,\[Alpha],i,j]={MeanValue2[QH[[i]]*D[Q[[j]],{z,2}]*1/(z-w),\[Alpha],rep2],MeanValue2[QH[[i]]*D[Q[[j]],{z,1}]*1/(z-w)^2,\[Alpha],rep2],MeanValue2[QH[[i]]*D[Q[[j]],{z,1}]*1/(z-w),\[Alpha],rep2],MeanValue2[QH[[i]]*Q[[j]]*1/(z-w)^3,\[Alpha],rep2],MeanValue2[QH[[i]]*Q[[j]]*1/(z-w)^2,\[Alpha],rep2],MeanValue2[QH[[i]]*Q[[j]]*1/(z-w),\[Alpha],rep2]};*)
(*ClearAll[OrthogRen2];*)
(*OrthogRen2[QHL_,QL_,QHR_,QR_,rep2_]:=OrthogRen2[QHL,QL,QHR,QR,rep2]=(1/Normalisation2[QHL,QL,QHR,QR,rep2])*Det[{RawOrthog2[QHL,QR,rep2,0,0,1,1],RawOrthog2[QHL,QR,rep2,0,0,1,2],RawOrthog2[QHL,QR,rep2,0,0,2,1],RawOrthog2[QHL,QR,rep2,0,0,2,2],RawOrthog2[QHL,QR,rep2,0,0,1,3],RawOrthog2[QHL,QR,rep2,0,0,3,1]}];*)


(* ::Title::Closed:: *)
(*RUN: \[CapitalPsi] space ONE SITE*)


(* ::Chapter::Closed:: *)
(*Definition Lax connection*)


(* ::Subsection::Closed:: *)
(*Lax connection*)


(* ::TextIndent:: *)
(*Lax connection Subscript[L, ij] with regulator \[Phi]*)


(* ::Input:: *)
(*ClearAll[LIndex1];*)
(*LIndex1[a_,b_,rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=LIndex1[a,b,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=1/(z-w1)*EM[a,b,rep1]+KroneckerDelta[a,1]*KroneckerDelta[b,1]*\[Phi]1*Id[rep1]+KroneckerDelta[a,2]*KroneckerDelta[b,2]*\[Phi]2*Id[rep1]+KroneckerDelta[a,3]*KroneckerDelta[b,3]*\[Phi]3*Id[rep1]+KroneckerDelta[a,4]*KroneckerDelta[b,4]*\[Phi]4*Id[rep1];*)


(* ::Subsubsection::Closed:: *)
(*Example of use*)


(* ::Input:: *)
(*ClearAll[rep1T];*)
(*rep1T=FUNDAMENTAL;*)
(*LIndex[1,1,rep1T,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal*)


(* ::Chapter::Closed:: *)
(*Hamiltonians*)


(* ::Section::Closed:: *)
(*H1*)


(* ::Subsection::Closed:: *)
(*Regulated Hamiltonian tr(L)*)


(* ::Input:: *)
(*ClearAll[H11];*)
(*H11[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=H11[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=Sum[LIndex1[a,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{a,1,Ng}]//Simplify;*)


(* ::Subsection::Closed:: *)
(*Definition Residue of H*)


(* ::Input:: *)
(*Clear[ResH11];*)
(*ResH11[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_]:=ResH11[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w]=SeriesCoefficient[H11[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,-1}]*)


(* ::Subsection::Closed:: *)
(*Example of use*)


(* ::Input:: *)
(*H11[FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal;*)


(* ::Section::Closed:: *)
(*H2*)


(* ::Subsection::Closed:: *)
(*Regulated Hamiltonian 3/2*tr(L')+trA(L^2)*)


(* ::Subsubsection::Closed:: *)
(*Antisymm. def*)


(* ::Input:: *)
(*ClearAll[trA1L2];*)
(*trA1L2[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trA1L2[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=1/2*(Sum[LIndex1[a,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[b,b,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[b,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[a,b,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{a,1,Ng},{b,1,Ng}]);*)


(* ::Subsubsection::Closed:: *)
(*Def. H2*)


(* ::Input:: *)
(*ClearAll[H12];*)
(*H12[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=H12[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=(trA1L2[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+3/2*Sum[D[LIndex1[a,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z],{a,1,Ng}])//Simplify*)


(* ::Subsection::Closed:: *)
(*Definition Residue of H*)


(* ::Input:: *)
(*Clear[ResH12];*)
(*ResH12[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_]:=ResH12[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w]=SeriesCoefficient[H12[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,-1}]*)


(* ::Input:: *)
(*Clear[ResH12O];*)
(*ResH12O[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_,o_]:=ResH12O[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w,o]=SeriesCoefficient[H12[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,o}]*)


(* ::Subsection::Closed:: *)
(*Example of use*)


(* ::Input:: *)
(*H12[FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal;*)


(* ::Section::Closed:: *)
(*H3*)


(* ::Subsection::Closed:: *)
(*H3 Unsimplified*)


(* ::TextIndent:: *)
(*Regulated Hamiltonian trA(L^3)+4/3 trA(L\[CenterDot]L')+1/2 trA(L'\[CenterDot]L)+tr(L'')*)


(* ::Subsubsection::Closed:: *)
(*Antisymm. def*)


(* ::Input:: *)
(*ClearAll[trA1L3];*)
(*trA1L3[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trA1L3[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=1/6*Sum[LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{i1,1,Ng},{i2,1,Ng},{i3,1,Ng}];*)


(* ::Input:: *)
(*ClearAll[trA1LpL,trA1LLp];*)
(*trA1LLp[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trA1LLp[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=1/2*(Sum[LIndex1[a,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex1[b,b,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]-LIndex1[b,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex1[a,b,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z],{a,1,Ng},{b,1,Ng}]);*)
(*trA1LpL[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trA1LpL[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=1/2*(Sum[D[LIndex1[a,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex1[b,b,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-D[LIndex1[b,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex1[a,b,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{a,1,Ng},{b,1,Ng}]);*)


(* ::Subsubsection::Closed:: *)
(*Def H3*)


(* ::Input:: *)
(*ClearAll[H13U];*)
(*H13U[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=H3U[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=trA1L3[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+4/3*trA1LLp[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+2/3*trA1LpL[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+Sum[D[LIndex1[a,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{z,2}],{a,1,Ng}]*)


(* ::Subsection::Closed:: *)
(*H3 Simplified*)


(* ::Input:: *)
(*ClearAll[H13];*)
(*H13[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=H13[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=H13U[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*)


(* ::Subsection::Closed:: *)
(*Definition Residue of H*)


(* ::Input:: *)
(*Clear[ResH13];*)
(*ResH13[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_]:=ResH13[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w]=SeriesCoefficient[H13[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,-1}]*)


(* ::Input:: *)
(*Clear[ResH13O];*)
(*ResH13O[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_,o_]:=ResH13O[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w,o]=SeriesCoefficient[H13[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,o}]*)


(* ::Subsection::Closed:: *)
(*Example of use*)


(* ::Input:: *)
(*H13[FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal//Simplify*)


(* ::Section::Closed:: *)
(*H4*)


(* ::Subsection::Closed:: *)
(*H4 Unsimplified*)


(* ::TextIndent:: *)
(*Regulated Hamiltonian trA(L^3)+4/3 trA(L\[CenterDot]L')+1/2 trA(L'\[CenterDot]L)+tr(L'')*)


(* ::Subsubsection::Closed:: *)
(*Antisymm. def*)


(* ::Input:: *)
(*ClearAll[trA1L4];*)
(*trA1L4[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trA1L4[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=1/4!*Sum[LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i4,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i4,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i3,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i4,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i4,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i2,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i4,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i3,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i4,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i2,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i4,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i4,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i3,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i4,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i4,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i1,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i4,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i3,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i4,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i1,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i4,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i4,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i2,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i4,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i4,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i1,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i4,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i2,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i4,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i1,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i4,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i3,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i4,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i2,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i4,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i3,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i4,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i1,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i4,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i2,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i4,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i1,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{i1,1,Ng},{i2,1,Ng},{i3,1,Ng},{i4,1,Ng}];*)


(* ::Input:: *)
(*(*L^3*)*)
(*ClearAll[trA1LLLp,trA1LLpL,trA1LpLL];*)
(*trA1LLLp[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trA1LLLp[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=1/6*Sum[LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]-LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]-LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]+LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]+LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]-LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z],{i1,1,Ng},{i2,1,Ng},{i3,1,Ng}];*)
(**)
(*trA1LLpL[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trA1LLpL[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=1/6*Sum[LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{i1,1,Ng},{i2,1,Ng},{i3,1,Ng}];*)
(**)
(*trA1LpLL[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trA1LpLL[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=1/6*Sum[D[LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-D[LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-D[LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+D[LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+D[LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-D[LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{i1,1,Ng},{i2,1,Ng},{i3,1,Ng}];*)


(* ::Input:: *)
(*(*L^2*)*)
(*ClearAll[trA1LLpp,trA1LpLp,trA1LppL];*)
(*trA1LLpp[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trA1LLpp[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=1/2*(Sum[LIndex1[a,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex1[b,b,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{z,2}]-LIndex1[b,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . D[LIndex1[a,b,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{z,2}],{a,1,Ng},{b,1,Ng}]);*)
(**)
(*trA1LpLp[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trA1LpLp[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=1/2*(Sum[D[LIndex1[a,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . D[LIndex1[b,b,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]-D[LIndex1[b,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z] . D[LIndex1[a,b,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z],{a,1,Ng},{b,1,Ng}]);*)
(**)
(*trA1LppL[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trA1LppL[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=1/2*(Sum[D[LIndex1[a,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{z,2}] . LIndex1[b,b,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-D[LIndex1[b,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{z,2}] . LIndex1[a,b,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{a,1,Ng},{b,1,Ng}]);*)


(* ::Subsubsection::Closed:: *)
(*Def H4*)


(* ::Input:: *)
(*ClearAll[H14U];*)
(*H14U[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=H14U[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=trA1L4[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+3/4*trA1LLLp[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+1/2*trA1LLpL[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+1/4*trA1LpLL[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+*)
(*1/2*trA1LLpp[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+1/2*trA1LpLp[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+1/6*trA1LppL[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+*)
(*1/4*Sum[D[LIndex1[a,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{z,3}],{a,1,Ng}]*)


(* ::Subsection::Closed:: *)
(*H4 Simplified*)


(* ::Input:: *)
(*ClearAll[H14];*)
(*H14[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=H14[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=H14U[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*)


(* ::Subsection::Closed:: *)
(*Definition Residue of H*)


(* ::Input:: *)
(*Clear[ResH14];*)
(*ResH14[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_]:=ResH14[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w]=SeriesCoefficient[H14[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,-1}]*)


(* ::Input:: *)
(*Clear[ResH14O];*)
(*ResH14O[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_,o_]:=ResH14O[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w,o]=SeriesCoefficient[H14[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,o}]*)


(* ::Subsection::Closed:: *)
(*Example of use*)


(* ::Input:: *)
(*H14[FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal//Simplify;*)


(* ::Chapter::Closed:: *)
(*Adjoint Hamiltonians TO DO*)


(* ::Section::Closed:: *)
(*H1\[Dagger]*)


(* ::Subsection::Closed:: *)
(*Regulated Hamiltonian tr(L)*)


(* ::Input:: *)
(*ClearAll[H1\[Dagger]];*)
(*H1\[Dagger][rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_][z_]:=H1\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]=H1[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z](*Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}]//Simplify;*)*)


(* ::Subsection::Closed:: *)
(*Definition Residue of H*)


(* ::Input:: *)
(*Clear[ResH1\[Dagger]];*)
(*ResH1\[Dagger][rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,w_]:=ResH1\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,w]=SeriesCoefficient[H1\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]//Normal,{z,w,-1}]*)


(* ::Subsubsection::Closed:: *)
(*Example of use*)


(* ::Input:: *)
(*ClearAll[rep1,rep2];*)
(*rep1=FUNDAMENTAL;*)
(*rep2=FUNDAMENTAL;*)
(*H1\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]//Normal;*)


(* ::Section::Closed:: *)
(*H2\[Dagger]*)


(* ::Subsection::Closed:: *)
(*Regulated Hamiltonian -tr(L')+1/2 (tr(L))^2-1/2 tr(L^2)*)


(* ::Input:: *)
(*ClearAll[H2\[Dagger]];*)
(*H2\[Dagger][rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_][z_]:=H2\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]=(-Sum[D[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],z],{a,1,3}]+1/2*(Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}]) . (Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}])-1/2*Sum[LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z] . LIndex[b,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3},{b,1,3}])//Simplify*)


(* ::Subsection::Closed:: *)
(*Definition Residue of H*)


(* ::Input:: *)
(*Clear[ResH2\[Dagger]];*)
(*ResH2\[Dagger][rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,w_]:=ResH2\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,w]=SeriesCoefficient[H2\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]//Normal,{z,w,-1}]*)


(* ::Section::Closed:: *)
(*H3\[Dagger]*)


(* ::Subsection::Closed:: *)
(*H3\[Dagger] simplified 1/3 tr(L'')+tr(L' . L)-2/3 tr(L') . tr(L)-1/2 tr(L) . tr(L^2)+1/6 tr(L)^3+1/3 tr(L^3)*)


(* ::Input:: *)
(*ClearAll[H3\[Dagger]];*)
(*H3\[Dagger][rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_][z_]:=H3\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]=(1/3*Sum[D[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{z,2}],{a,1,3}]+Sum[D[LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],z] . LIndex[b,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3},{b,1,3}]-2/3*Sum[D[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],z],{a,1,3}] . Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}]-1/2*Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}] . Sum[LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z] . LIndex[b,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3},{b,1,3}]+1/6*(Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}]) . (Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}]) . (Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}])+1/3*Sum[LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z] . LIndex[b,c,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z] . LIndex[c,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3},{b,1,3},{c,1,3}])//Simplify;*)


(* ::Subsection::Closed:: *)
(*Definition Residue of H*)


(* ::Input:: *)
(*Clear[ResH3\[Dagger]];*)
(*ResH3\[Dagger][rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,w_]:=ResH3\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,w]=SeriesCoefficient[H3[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]//Normal,{z,w,-1}]*)


(* ::Subsection::Closed:: *)
(*Example of use*)


(* ::Input:: *)
(*Series[H3\[Dagger][ADJOINT,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]//Normal,{z,0,-3}]//Simplify//Normal;*)


(* ::Chapter::Closed:: *)
(*Eigenvectors (FIRST RUN Q SPACE)*)


(* ::Section::Closed:: *)
(*Right eigenvectors (column)*)


(* ::Input:: *)
(*ClearAll[VecR1];*)
(*VecR1[rep1_]:=VecR1[rep1]=Table[{v[j]},{j,1,dim[rep1]}];*)


(* ::Input:: *)
(*ClearAll[EigenR1];*)
(*EigenR1[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=EigenR1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=Table[Quiet[Flatten[VecR1[rep1]/.Solve[*)
(*CoefficientList[((H12[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . VecR1[rep1]-HSpectrum1[2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z][[idx]]* VecR1[rep1]))//Together//Numerator,z]==0&&*)
(*CoefficientList[((H13[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . VecR1[rep1]-HSpectrum1[3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z][[idx]]* VecR1[rep1]))//Together//Numerator,z]==0&&*)
(*CoefficientList[((H14[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] . VecR1[rep1]-HSpectrum1[4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z][[idx]]* VecR1[rep1]))//Together//Numerator,z]==0,Flatten[VecR1[rep1]]],1]],{idx,1,Dimensions[HSpectrum1[1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]][[1]]}]*)


(* ::Subsubsection::Closed:: *)
(*Example*)


(* ::Input:: *)
(*EigenR1[FUNDAMENTAL,FUNDAMENTAL,3,5,7,9]*)


(* ::Section::Closed:: *)
(*Left eigenvectors (row)*)


(* ::Input:: *)
(*ClearAll[VecL1];*)
(*VecL1[rep1_]:=VecL1[rep1]=Transpose[VecR1[rep1]];*)


(* ::Input:: *)
(*ClearAll[EigenL1];*)
(*EigenL1[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=EigenL1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=Table[Quiet[Flatten[VecL1[rep1]/.Solve[*)
(*CoefficientList[((VecL1[rep1] . H12[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-HSpectrum1[2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z][[idx]]* VecL1[rep1]))//Together//Numerator,z]==0&&*)
(*CoefficientList[((VecL1[rep1] . H13[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-HSpectrum1[3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z][[idx]]* VecL1[rep1]))//Together//Numerator,z]==0&&*)
(*CoefficientList[((VecL1[rep1] . H14[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-HSpectrum1[4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z][[idx]]* VecL1[rep1]))//Together//Numerator,z]==0,Flatten[VecL1[rep1]]],1]],{idx,1,Dimensions[HSpectrum1[1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]][[1]]}]*)


(* ::Subsubsection::Closed:: *)
(*Example*)


(* ::Input:: *)
(*EigenL1[FUNDAMENTAL,FUNDAMENTAL,3,5,7,9]*)


(* ::Section::Closed:: *)
(*Normalisation*)


(* ::Subsection::Closed:: *)
(*Highest weight state*)


(* ::Text:: *)
(*h.w.s. 1-site*)


(* ::Input:: *)
(*ClearAll[UPR1,UPL1];*)
(*UPR1[rep1_]:=UPR1[rep1]=Table[{KroneckerDelta[j,1]},{j,1,dim[rep1]}];*)
(*UPL1[rep1_]:=UPL1[rep1]=Transpose[UPR1[rep1]];*)


(* ::Subsection::Closed:: *)
(*Vacuum vector*)


(* ::Text:: *)
(*Vacuum vector 1-site*)


(* ::Input:: *)
(*ClearAll[VacR1,VacL1];*)
(*VacR1[rep1_,\[Beta]12_,\[Beta]13_,\[Beta]14_,\[Beta]23_,\[Beta]24_,\[Beta]34_]:=VacR1[rep1,\[Beta]23,\[Beta]12,\[Beta]13,\[Beta]14,\[Beta]23,\[Beta]24,\[Beta]34]=*)
(*MatrixExp[\[Beta]24(EM[4,2,rep1])] . MatrixExp[\[Beta]23(EM[3,2,rep1])] . MatrixExp[\[Beta]14(EM[4,1,rep1])] . MatrixExp[\[Beta]13(EM[3,1,rep1])] . MatrixExp[\[Beta]34(EM[4,3,rep1])] . MatrixExp[\[Beta]12(EM[2,1,rep1])] . UPR1[rep1];*)
(*VacL1[rep1_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=VacL1[rep1,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=UPL1[rep1] . MatrixExp[\[Alpha]12(EM[1,2,rep1])] . MatrixExp[\[Alpha]34(EM[3,4,rep1])] . MatrixExp[\[Alpha]13(EM[1,3,rep1])] . MatrixExp[\[Alpha]14(EM[1,4,rep1])] . MatrixExp[\[Alpha]23(EM[2,3,rep1])] . MatrixExp[\[Alpha]24(EM[2,4,rep1])];*)


(* ::Subsection::Closed:: *)
(*Renormalisation*)


(* ::Input:: *)
(*ClearAll[ReNormR1,ReNormL1,\[Alpha]23,\[Alpha]13,\[Alpha]12,\[Alpha]24,\[Alpha]14,\[Alpha]34,\[Beta]23,\[Beta]13,\[Beta]12,\[Beta]24,\[Beta]14,\[Beta]34];*)
(*ReNormR1[rep1_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_,\[Psi]_]:=ReNormR1[rep1,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34,\[Psi]]=(\[Psi]*(1/(((VacL1[rep1,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34] . \[Psi])//Flatten)[[1]])))//Simplify;*)
(*ReNormL1[rep1_,\[Beta]12_,\[Beta]13_,\[Beta]14_,\[Beta]23_,\[Beta]24_,\[Beta]34_,\[Psi]_]:=ReNormL1[rep1,\[Beta]12,\[Beta]13,\[Beta]14,\[Beta]23,\[Beta]24,\[Beta]34,\[Psi]]=(\[Psi]*(1/((\[Psi] . VacR1[rep1,\[Beta]12,\[Beta]13,\[Beta]14,\[Beta]23,\[Beta]24,\[Beta]34]//Flatten)[[1]])))//Simplify;*)


(* ::Subsubsection::Closed:: *)
(*Example of use*)


(* ::Input:: *)
(*ReNormR1[FUNDAMENTAL,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34,EigenR1[FUNDAMENTAL,FUNDAMENTAL,3,5,7,9][[1]]]*)


(* ::Input:: *)
(*ReNormL1[FUNDAMENTAL,\[Beta]12,\[Beta]13,\[Beta]14,\[Beta]23,\[Beta]24,\[Beta]34,EigenL1[FUNDAMENTAL,FUNDAMENTAL,3,5,7,9][[1]]]*)


(* ::Input:: *)
(*((ReNormL1[FUNDAMENTAL,\[Beta]12,\[Beta]13,\[Beta]14,\[Beta]23,\[Beta]24,\[Beta]34,EigenL1[FUNDAMENTAL,FUNDAMENTAL,3,5,7,9][[1]]] . ReNormR1[FUNDAMENTAL,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34,EigenR1[FUNDAMENTAL,FUNDAMENTAL,3,5,7,9][[1]]])//Flatten)[[1]]//FullSimplify*)


(* ::Subsection::Closed:: *)
(*Correct normalisation*)


(* ::Subsubsection::Closed:: *)
(*\[Beta] without fixing \[Alpha]*)


(* ::Input:: *)
(*ClearAll[\[Beta]rT,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];*)
(*\[Beta]rT[12][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=\[Beta]rT[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=(\[Alpha]23 (\[Phi]1-\[Phi]4))/(\[Alpha]13 (\[Phi]2-\[Phi]4));*)
(*\[Beta]rT[13][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=\[Beta]rT[13][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=((\[Phi]2-\[Phi]3) (\[Phi]1-\[Phi]4))/(\[Alpha]13 (\[Phi]1-\[Phi]2) (\[Phi]3-\[Phi]4));*)
(*\[Beta]rT[14][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=\[Beta]rT[14][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=-(((\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))/(\[Alpha]14 (\[Phi]1-\[Phi]2) (\[Phi]3-\[Phi]4)));*)
(*\[Beta]rT[23][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=\[Beta]rT[23][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=-(((\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))/(\[Alpha]23 (\[Phi]1-\[Phi]2) (\[Phi]3-\[Phi]4)));*)
(*\[Beta]rT[24][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=\[Beta]rT[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=(\[Alpha]13 (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))/(\[Alpha]14 \[Alpha]23 (\[Phi]1-\[Phi]2) (\[Phi]3-\[Phi]4));*)
(*\[Beta]rT[34][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=\[Beta]rT[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=(\[Alpha]13 (\[Phi]2-\[Phi]4))/(\[Alpha]14 (-\[Phi]2+\[Phi]3));*)


(* ::Subsubsection::Closed:: *)
(*fixing \[Alpha]*)


(* ::Input:: *)
(*ClearAll[\[Alpha]r,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];*)
(*\[Alpha]r[12][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=(\[Alpha]13 (-\[Phi]1+\[Phi]3))/(\[Alpha]23 (\[Phi]2-\[Phi]3));*)
(*\[Alpha]r[24][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=(\[Alpha]14 \[Alpha]23 (\[Phi]2-\[Phi]3) (\[Phi]1-\[Phi]4))/(\[Alpha]13 (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4));*)
(*\[Alpha]r[34][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=(\[Alpha]14 (\[Phi]1-\[Phi]4))/(\[Alpha]13 (\[Phi]1-\[Phi]3));*)


(* ::Subsubsection::Closed:: *)
(*\[Beta] with fixed \[Alpha]*)


(* ::Input:: *)
(*ClearAll[\[Beta]r,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];*)
(*\[Beta]r[12][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Beta]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=\[Beta]rT[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]];*)
(*\[Beta]r[13][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Beta]r[13][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=\[Beta]rT[13][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]];*)
(*\[Beta]r[14][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Beta]r[14][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=\[Beta]rT[14][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]];*)
(*\[Beta]r[23][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Beta]r[23][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=\[Beta]rT[23][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]];*)
(*\[Beta]r[24][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Beta]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=\[Beta]rT[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]];*)
(*\[Beta]r[34][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Beta]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=\[Beta]rT[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]];*)


(* ::Subsection::Closed:: *)
(*Renormalised eigenvectors*)


(* ::Input:: *)
(*ClearAll[ReEigenR1,ReEigenL1];*)
(*ReEigenR1[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=ReEigenR1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=ReNormR1[rep1,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],#]&/@EigenR1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];*)
(*ReEigenL1[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=ReEigenL1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=ReNormL1[rep1,\[Beta]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Beta]r[13][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Beta]r[14][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Beta]r[23][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Beta]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Beta]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],#]&/@EigenL1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];*)


(* ::Title::Closed:: *)
(*RUN: CREATE BASIS ONE SITE*)


(* ::Section::Closed:: *)
(*Q Base*)


(* ::Subsection::Closed:: *)
(*Q Base given repAll*)


(* ::Input:: *)
(*ClearAll[QBase1,T];*)
(*QBase1[T_][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=QBase1[T][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=*)
(*Q1[T][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]/.SolveBaxAll1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1];*)


(* ::Section::Closed:: *)
(*Q Hodge Base*)


(* ::Subsection::Closed:: *)
(*Q Hodge Base given repAll*)


(* ::Input:: *)
(*ClearAll[QHBase1,T,U];*)
(*QHBase1[T_][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=QHBase1[T][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=Table[Hodge1[{QBase1[1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]],QBase1[2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]],QBase1[3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]],QBase1[4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]]},rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,z,T],{U,1,Length[QBase1[1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]]}];*)
(**)
(*(*Q1\[Dagger][T][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]/.SolveBaxAllh1\[Dagger][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h1[3,3]/.SolveBaxAll1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1],h1[4,3]/.SolveBaxAll1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1],h1[4,5]/.SolveBaxAll1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1]];*)*)
(**)
(*(*With Hodge*)*)
(*(*Table[Hodge[{QBase1[1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]],QBase1[2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]],QBase1[3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]],QBase1[4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]]},rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,z,T],{U,1,Length[QBase1[1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]]}];*)*)


(* ::Section::Closed:: *)
(*Q Vacuum*)


(* ::Input:: *)
(*ClearAll[QVac1,QHVac1,QvQv1,QvQv1\[Omega]\[Phi]];*)
(*(*QvQv*)*)
(*(*VACUUM STATES*)*)
(*QVac1[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=QVac1[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]={QBase1[1][rep1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]],QBase1[2][rep1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]],QBase1[3][rep1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]],QBase1[4][rep1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]]};*)
(*QHVac1[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=QHVac1[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]={QHBase1[1][rep1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]],QHBase1[2][rep1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]],QHBase1[3][rep1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]],QHBase1[4][rep1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]]};*)
(**)
(*QvQv1[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=QvQv1[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=OrthogRen1[QHVac1[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4],QVac1[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4],QHVac1[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4],QVac1[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4],rep1];*)
(**)
(*QvQv1\[Omega]\[Phi][rep1_,\[Omega]1_,\[Omega]2_,\[Omega]3_,\[Omega]4_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=QvQv1\[Omega]\[Phi][rep1,\[Omega]1,\[Omega]2,\[Omega]3,\[Omega]4,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=OrthogRen1[QHVac1[rep1,\[Omega]1,\[Omega]2,\[Omega]3,\[Omega]4],QVac1[rep1,\[Omega]1,\[Omega]2,\[Omega]3,\[Omega]4],QHVac1[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4],QVac1[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4],rep1];*)


(* ::Title:: *)
(*K vector 2-sites*)


(* ::Section::Closed:: *)
(*RUN: Definition K vector*)


(* ::Subsection::Closed:: *)
(*<K|*)


(* ::Subsubsection::Closed:: *)
(*For {p,p,0,0}, ONLY SU(2) scalar product*)


(* ::Input:: *)
(*ClearAll[K1];*)
(*K1[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]\[Phi]13_,\[Alpha]\[Phi]14_,\[Alpha]\[Phi]23_,k1_,k2_,k3_,k4_]:=K1[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4]=UPL1[rep1] . MatrixExp[Log[p11]*EM[1,1,rep1]] . MatrixExp[Log[1/p11]*EM[2,2,rep1]] . MatrixExp[Log[p33]*EM[3,3,rep1]] . MatrixExp[Log[p44]*EM[4,4,rep1]] . MatrixExp[p12(EM[1,2,rep1])] . MatrixExp[p34(EM[3,4,rep1])] . MatrixExp[(\[Alpha]\[Phi]13 (k4-\[Phi]2-\[Phi]3))/(k4-\[Phi]1-\[Phi]2) (EM[1,3,rep1])] . MatrixExp[(\[Alpha]\[Phi]14 (k4-\[Phi]2-\[Phi]4))/(k4-\[Phi]1-\[Phi]2) (EM[1,4,rep1])] . MatrixExp[(\[Alpha]\[Phi]23 (k4-\[Phi]1-\[Phi]3))/(k4-\[Phi]1-\[Phi]2) (EM[2,3,rep1])] . MatrixExp[(\[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (\[Phi]2-\[Phi]3) (k4-\[Phi]1-\[Phi]4) (\[Phi]1-\[Phi]4))/(\[Alpha]\[Phi]13 (-k4+\[Phi]1+\[Phi]2) (-\[Phi]1+\[Phi]3) (\[Phi]2-\[Phi]4)) (EM[2,4,rep1])];*)
(**)
(*K22[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]\[Phi]13_,\[Alpha]\[Phi]14_,\[Alpha]\[Phi]23_,k1_,k2_,k3_,k4_]:=K22[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4]=KroneckerProduct[K1[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4],UPL1[rep2] . MatrixExp[Log[P11]*EM[1,1,rep2]] . MatrixExp[Log[P22]*EM[2,2,rep2]] . MatrixExp[Log[P33]*EM[3,3,rep2]] . MatrixExp[Log[P44]*EM[4,4,rep2]] . MatrixExp[P12(EM[1,2,rep1])] . MatrixExp[P34(EM[3,4,rep1])] . MatrixExp[P13(EM[1,3,rep1])] . MatrixExp[P14(EM[1,4,rep1])] . MatrixExp[P23(EM[2,3,rep1])] . MatrixExp[P24(EM[2,4,rep1])]];*)


(* ::Subsubsection::Closed:: *)
(*For {p,p,0,0}, SU(2) and SU(3) scalar product*)


(* ::Text:: *)
(*Using SU(2) scalar product and SU(3) scalar product*)


(* ::Input:: *)
(**)
(*ClearAll[K1B];*)
(*K1B[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]\[Phi]13_,\[Alpha]\[Phi]14_,\[Alpha]\[Phi]23_,k1_,k2_,k3_,k4_]:=K1B[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4]=UPL1[rep1] . MatrixExp[Log[p11]*EM[1,1,rep1]] . MatrixExp[Log[1/p11]*EM[2,2,rep1]] . MatrixExp[Log[p33]*EM[3,3,rep1]] . MatrixExp[Log[p44]*EM[4,4,rep1]] . MatrixExp[p12(EM[1,2,rep1])] . MatrixExp[p34(EM[3,4,rep1])] . MatrixExp[(\[Alpha]\[Phi]13 (k3-\[Phi]3) (-k4+\[Phi]2+\[Phi]3))/((k3-\[Phi]1) (-k4+\[Phi]1+\[Phi]2)) (EM[1,3,rep1])] . MatrixExp[(\[Alpha]\[Phi]14 (k3-\[Phi]4) (-k4+\[Phi]2+\[Phi]4))/((k3-\[Phi]1) (-k4+\[Phi]1+\[Phi]2)) (EM[1,4,rep1])] . MatrixExp[(\[Alpha]\[Phi]23 (k3-\[Phi]3) (-k4+\[Phi]1+\[Phi]3))/((k3-\[Phi]2) (-k4+\[Phi]1+\[Phi]2)) (EM[2,3,rep1])] . MatrixExp[-((\[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (\[Phi]2-\[Phi]3) (k3-\[Phi]4) (-\[Phi]1+\[Phi]4) (-k4+\[Phi]1+\[Phi]4))/(\[Alpha]\[Phi]13 (k3-\[Phi]2) (-k4+\[Phi]1+\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4)))(EM[2,4,rep1])];*)
(**)
(**)
(*K22B[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]\[Phi]13_,\[Alpha]\[Phi]14_,\[Alpha]\[Phi]23_,k1_,k2_,k3_,k4_]:=K22B[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4]=KroneckerProduct[K1B[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4],UPL1[rep2] . MatrixExp[Log[P11]*EM[1,1,rep2]] . MatrixExp[Log[P22]*EM[2,2,rep2]] . MatrixExp[Log[P33]*EM[3,3,rep2]] . MatrixExp[Log[P44]*EM[4,4,rep2]] . MatrixExp[P12(EM[1,2,rep1])] . MatrixExp[P34(EM[3,4,rep1])] . MatrixExp[P13(EM[1,3,rep1])] . MatrixExp[P14(EM[1,4,rep1])] . MatrixExp[P23(EM[2,3,rep1])] . MatrixExp[P24(EM[2,4,rep1])]];*)
(**)


(* ::Section:: *)
(*RUN (+experiment): Modules for testing*)


(* ::Subsection::Closed:: *)
(*RUN: Definition SU(2) scalar product*)


(* ::Input:: *)
(*ClearAll[SU2MeanValue];*)
(*SU2MeanValue[f_,\[Tau]_,l1_,l2_]:=*)
(*2 \[Pi] I Residue[(f*1/((z-w1)^(2 l1 +1) (z-w2)^(2 l2 +1))),{z,\[Tau]}];*)


(* ::Input:: *)
(*ClearAll[SU2MeanValue1];*)
(*SU2MeanValue1[f_,\[Tau]_,l1_]:=*)
(*2 \[Pi] I Residue[(f*1/(z-w1)^(2 l1+1 ) ),{z,\[Tau]}];*)


(* ::Input:: *)
(*(*2-sites scalar product*)*)
(*ClearAll[SU2OrthogRen];*)
(*SU2OrthogRen[Q1HODGE_,Q2_,l1_,l2_]:=(SU2MeanValue[Q1HODGE*Q2,w1,l1,l2]*SU2MeanValue[Q1HODGE*Q2*z,w2,l1,l2]-SU2MeanValue[Q1HODGE*Q2,w2,l1,l2]*SU2MeanValue[Q1HODGE*Q2*z,w1,l1,l2])/((Limit[Q1HODGE,z->w1])*(Limit[Q1HODGE,z->w2])*(Limit[Q2,z->w1])*(Limit[Q2,z->w2]))//Simplify;*)


(* ::Input:: *)
(*(*1-site scalar product*)*)
(*ClearAll[SU2OrthogRen1];*)
(*SU2OrthogRen1[Q1HODGE_,Q2_,l1_]:=(SU2MeanValue1[Q1HODGE*Q2,w1,l1])/((Limit[Q1HODGE,z->w1])*(Limit[Q2,z->w1]))//Simplify;*)


(* ::Subsection::Closed:: *)
(*RUN: Definition SU(3) scalar product*)


(* ::Input:: *)
(*ClearAll[SU3MeanValue];*)
(*SU3MeanValue[f_,\[Tau]_,rep1_,rep2_]:=*)
(*2 \[Pi] I Residue[(f*1/((z-w1)^(rep1[[1]]-rep1[[3]]) (z-w2)^(rep2[[1]]-rep2[[3]]))),{z,\[Tau]}];*)


(* ::Input:: *)
(*ClearAll[SU3MeanValue1];*)
(*SU3MeanValue1[f_,\[Tau]_,rep1_]:=*)
(*2 \[Pi] I Residue[(f*1/(z-w1)^(rep1[[1]]-rep1[[3]])),{z,\[Tau]}];;*)


(* ::Subsubsection::Closed:: *)
(*2-sites scalar product*)


(* ::Text:: *)
(*Input:*)
(*Q={Subscript[Q, 1],Subscript[Q, 2],Subscript[Q, 3]}*)
(*QH={Subscript[QH, 1],Subscript[QH, 2],Subscript[QH, 3]}*)


(* ::Input:: *)
(*ClearAll[SU3Normalisation];*)
(*SU3Normalisation[QHL_,QL_,QHR_,QR_]:=SU3Normalisation[QHL,QL,QHR,QR]=(Limit[QHL[[1]],z->0,Direction->"FromAbove"])*(Limit[QL[[3]],z->0,Direction->"FromAbove"])*(Limit[QHR[[3]],z->0,Direction->"FromAbove"])*(Limit[QR[[1]],z->0,Direction->"FromAbove"])**)
(*(Limit[QHL[[1]],z->1,Direction->"FromBelow"])*(Limit[QL[[3]],z->1,Direction->"FromBelow"])*(Limit[QHR[[3]],z->1,Direction->"FromBelow"])*(Limit[QR[[1]],z->1,Direction->"FromBelow"]);*)


(* ::Input:: *)
(*ClearAll[SU3OrthogRen];*)
(*SU3OrthogRen[QHL_,QL_,QHR_,QR_,rep1_,rep2_]:=SU3OrthogRen[QHL,QL,QHR,QR,rep1,rep2]=(1/SU3Normalisation[QHL,QL,QHR,QR])*Det[{{SU3MeanValue[QHL[[1]]*D[QR[[1]],z]*1/z,0,rep1,rep2],SU3MeanValue[QHL[[1]]*D[QR[[1]],z]*1/(z-1),0,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[1]]*1/z^2,0,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[1]]*1/(z-1)^2,0,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[1]]*1/z,0,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[1]]*1/(z-1),0,rep1,rep2]},*)
(**)
(*{SU3MeanValue[QHL[[1]]*D[QR[[2]],z]*1/z,0,rep1,rep2],SU3MeanValue[QHL[[1]]*D[QR[[2]],z]*1/(z-1),0,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[2]]*1/z^2,0,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[2]]*1/(z-1)^2,0,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[2]]*1/z,0,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[2]]*1/(z-1),0,rep1,rep2]},*)
(**)
(*{SU3MeanValue[QHL[[2]]*D[QR[[1]],z]*1/z,0,rep1,rep2],SU3MeanValue[QHL[[2]]*D[QR[[1]],z]*1/(z-1),0,rep1,rep2],SU3MeanValue[QHL[[2]]*QR[[1]]*1/z^2,0,rep1,rep2],SU3MeanValue[QHL[[2]]*QR[[1]]*1/(z-1)^2,0,rep1,rep2],SU3MeanValue[QHL[[2]]*QR[[1]]*1/z,0,rep1,rep2],SU3MeanValue[QHL[[2]]*QR[[1]]*1/(z-1),0,rep1,rep2]},*)
(**)
(**)
(**)
(**)
(**)
(*{SU3MeanValue[QHL[[1]]*D[QR[[1]],z]*1/z,1,rep1,rep2],SU3MeanValue[QHL[[1]]*D[QR[[1]],z]*1/(z-1),1,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[1]]*1/z^2,1,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[1]]*1/(z-1)^2,1,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[1]]*1/z,1,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[1]]*1/(z-1),1,rep1,rep2]},*)
(**)
(*{SU3MeanValue[QHL[[1]]*D[QR[[2]],z]*1/z,1,rep1,rep2],SU3MeanValue[QHL[[1]]*D[QR[[2]],z]*1/(z-1),1,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[2]]*1/z^2,1,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[2]]*1/(z-1)^2,1,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[2]]*1/z,1,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[2]]*1/(z-1),1,rep1,rep2]},*)
(**)
(*{SU3MeanValue[QHL[[2]]*D[QR[[1]],z]*1/z,1,rep1,rep2],SU3MeanValue[QHL[[2]]*D[QR[[1]],z]*1/(z-1),1,rep1,rep2],SU3MeanValue[QHL[[2]]*QR[[1]]*1/z^2,1,rep1,rep2],SU3MeanValue[QHL[[2]]*QR[[1]]*1/(z-1)^2,1,rep1,rep2],SU3MeanValue[QHL[[2]]*QR[[1]]*1/z,1,rep1,rep2],SU3MeanValue[QHL[[2]]*QR[[1]]*1/(z-1),1,rep1,rep2]}}]//Simplify;*)


(* ::Subsubsection::Closed:: *)
(*1-site scalar product*)


(* ::Input:: *)
(*ClearAll[SU3Normalisation1];*)
(*SU3Normalisation1[QHL_,QL_,QHR_,QR_]:=SU3Normalisation1[QHL,QL,QHR,QR]=(Limit[QHL[[1]],z->0,Direction->"FromAbove"])*(Limit[QL[[3]],z->0,Direction->"FromAbove"])*(Limit[QHR[[3]],z->0,Direction->"FromAbove"])*(Limit[QR[[1]],z->0,Direction->"FromAbove"]);*)


(* ::Input:: *)
(*ClearAll[SU3OrthogRen1];*)
(*SU3OrthogRen1[QHL_,QL_,QHR_,QR_,rep1_]:=SU3OrthogRen1[QHL,QL,QHR,QR,rep1]=(1/SU3Normalisation1[QHL,QL,QHR,QR])*Det[{{SU3MeanValue1[QHL[[1]]*D[QR[[1]],z]*1/z,0,rep1],SU3MeanValue1[QHL[[1]]*QR[[1]]*1/z^2,0,rep1],SU3MeanValue1[QHL[[1]]*QR[[1]]*1/z,0,rep1]},*)
(**)
(*{SU3MeanValue1[QHL[[1]]*D[QR[[2]],z]*1/z,0,rep1],SU3MeanValue1[QHL[[1]]*QR[[2]]*1/z^2,0,rep1],SU3MeanValue1[QHL[[1]]*QR[[2]]*1/z,0,rep1]},*)
(**)
(*{SU3MeanValue1[QHL[[2]]*D[QR[[1]],z]*1/z,0,rep1],SU3MeanValue1[QHL[[2]]*QR[[1]]*1/z^2,0,rep1],SU3MeanValue1[QHL[[2]]*QR[[1]]*1/z,0,rep1]}}]//Simplify;*)


(* ::Subsection:: *)
(*IF YOU WANT TO EXPERIMENT MODIFY THIS: Definition scalar product for {p,p,0,0}x{q,q,0,0}*)


(* ::Subsubsection::Closed:: *)
(*Definition OrthogRen22*)


(* ::Input:: *)
(*ClearAll[SU3Normalisation1Sub];*)
(*SU3Normalisation1Sub[QHL_,QL_,QHR_,QR_,rep1_]:=SU3Normalisation1Sub[QHL,QL,QHR,QR,rep1]=(Limit[QHL[[3]],z->0,Direction->"FromAbove"])*(Limit[D[QHL[[3]],z],z->0,Direction->"FromAbove"])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),z->0,Direction->"FromAbove"]);*)


(* ::Input:: *)
(*ClearAll[SU3OrthogRen1Sub];*)
(*SU3OrthogRen1Sub[QHL_,QL_,QHR_,QR_,rep1_]:=SU3OrthogRen1Sub[QHL,QL,QHR,QR,rep1]=(1/SU3Normalisation1Sub[QHL,QL,QHR,QR,rep1])*Det[{{SU3MeanValue1[QHL[[3]]*D[QR[[1]],z]*1/z,0,rep1],SU3MeanValue1[QHL[[3]]*QR[[1]]*1/z^2,0,rep1]},*)
(**)
(*{SU3MeanValue1[QHL[[3]]*D[QR[[2]],z]*1/z,0,rep1],SU3MeanValue1[QHL[[3]]*QR[[2]]*1/z^2,0,rep1]}}]//Simplify;*)


(* ::Input:: *)
(*ClearAll[OrthogRen22];*)
(*OrthogRen22[QHL_,QL_,QHR_,QR_,rep1_,rep2_]:=OrthogRen22[QHL,QL,QHR,QR,rep1,rep2]=SU2OrthogRen[QHL[[4]],Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),(rep1[[1]]-rep1[[4]])/2,(rep2[[1]]-rep2[[4]])/2];*)
(**)
(*(*22 July Guess*)
(*SU2OrthogRen[QHL[[4]],Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),(rep1[[1]]-rep1[[4]])/2,(rep2[[1]]-rep2[[4]])/2];*)
(*My attempt using SU(3) 2-sites*)
(*(*Another guess*)*)
(*SU2OrthogRen1[QHL[[4]],Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),(rep1[[1]]-rep1[[4]])/2];*)
(**)*)
(**)
(**)
(**)
(*(*Main Guess*)*)
(*(**)
(*SU3OrthogRen1Sub[QHL,QL,QHR,QR,rep1]*SU2OrthogRen[QHL[[4]],Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),(rep1[[1]]-rep1[[4]])/2,(rep2[[1]]-rep2[[4]])/2];*)
(**)*)
(**)
(**)
(**)
(*(*My guess*)*)
(*(**)
(*SU2OrthogRen1[QHL[[4]],Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),(rep1[[1]]-rep1[[4]])/2]*SU2OrthogRen[QHL[[3]],Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),(rep1[[1]]-rep1[[4]])/2,(rep2[[1]]-rep2[[4]])/2];*)
(**)*)
(**)
(**)
(**)
(**)
(**)
(*(*NEW GUESSES*)*)
(*(*First NEW GUESS*)*)
(*(**)
(*SU3OrthogRen[{QHL[[2]],QHL[[3]],QHL[[4]]},{QL[[2]],QL[[3]],QL[[4]]},{QHR[[1]],QHR[[2]],QHR[[3]]},{QR[[1]],QR[[2]],QR[[3]]},rep1,rep2]*SU2OrthogRen1[QHL[[4]],Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),(rep1[[1]]-rep1[[4]])/2];*)
(**)*)
(**)
(*(*Second NEW GUESS*)*)
(*(**)
(*SU3OrthogRen1[{QHL[[2]],QHL[[3]],QHL[[4]]},{QL[[2]],QL[[3]],QL[[4]]},{QHR[[1]],QHR[[2]],QHR[[3]]},{QR[[1]],QR[[2]],QR[[3]]},rep1]*SU2OrthogRen[QHL[[4]],Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),(rep1[[1]]-rep1[[4]])/2,(rep2[[1]]-rep2[[4]])/2]*)
(**)*)
(**)
(**)
(**)
(**)
(**)
(**)
(**)
(**)
(**)
(**)
(**)
(*(*OLD GUESS*)*)
(*(*SU3OrthogRen[{QHL[[2]],QHL[[3]],QHL[[4]]},{QL[[2]],QL[[3]],QL[[4]]},{QHR[[2]],QHR[[3]],QHR[[4]]},{QR[[2]],QR[[3]],QR[[4]]},rep1,rep2]*SU2OrthogRen1[QHL[[4]],QR[[2]],(rep1[[1]]-rep1[[4]])/2];*)
(**)*)
(**)
(*(*SU3OrthogRen[{QHL[[2]],QHL[[3]],QHL[[4]]},{QL[[2]],QL[[3]],QL[[4]]},{QHR[[2]],QHR[[3]],QHR[[4]]},{QR[[2]],QR[[3]],QR[[4]]},rep1,rep2]*SU3OrthogRen1[{QHL[[3]],QHL[[4]],QHL[[1]]},{QL[[3]],QL[[4]],QL[[1]]},{QHR[[2]],QHR[[3]],QHR[[4]]},{QR[[2]],QR[[3]],QR[[4]]},rep1];*)
(**)*)
(**)
(*(*SU2OrthogRen1[QHL[[2]],QR[[1]],(rep1[[1]]-rep1[[4]])/2]*SU2OrthogRen[QHL[[3]],QR[[1]],(rep1[[1]]-rep1[[4]])/2,(rep2[[1]]-rep2[[4]])/2];*)*)


(* ::Subsubsection::Closed:: *)
(*Definition OrthogRen22A*)


(* ::Input:: *)
(*ClearAll[OrthogRen22A];*)
(*OrthogRen22A[QHL_,QL_,QHR_,QR_,rep1_,rep2_]:=OrthogRen22A[QHL,QL,QHR,QR,rep1,rep2]=SU2OrthogRen1[QHL[[4]],Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),(rep1[[1]]-rep1[[4]])/2];*)
(**)
(*(*22 July Guess*)
(*SU2OrthogRen[QHL[[4]],Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),(rep1[[1]]-rep1[[4]])/2,(rep2[[1]]-rep2[[4]])/2];*)
(*My attempt using SU(3) 2-sites*)
(*(*Another guess*)*)
(*SU2OrthogRen1[QHL[[4]],Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),(rep1[[1]]-rep1[[4]])/2];*)
(**)*)
(**)
(**)
(**)
(*(*Main Guess*)*)
(*(**)
(*SU3OrthogRen1Sub[QHL,QL,QHR,QR,rep1]*SU2OrthogRen[QHL[[4]],Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),(rep1[[1]]-rep1[[4]])/2,(rep2[[1]]-rep2[[4]])/2];*)
(**)*)
(**)
(**)
(**)
(*(*My guess*)*)
(*(**)
(*SU2OrthogRen1[QHL[[4]],Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),(rep1[[1]]-rep1[[4]])/2]*SU2OrthogRen[QHL[[3]],Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),(rep1[[1]]-rep1[[4]])/2,(rep2[[1]]-rep2[[4]])/2];*)
(**)*)
(**)
(**)
(**)
(**)
(**)
(*(*NEW GUESSES*)*)
(*(*First NEW GUESS*)*)
(*(**)
(*SU3OrthogRen[{QHL[[2]],QHL[[3]],QHL[[4]]},{QL[[2]],QL[[3]],QL[[4]]},{QHR[[1]],QHR[[2]],QHR[[3]]},{QR[[1]],QR[[2]],QR[[3]]},rep1,rep2]*SU2OrthogRen1[QHL[[4]],Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),(rep1[[1]]-rep1[[4]])/2];*)
(**)*)
(**)
(*(*Second NEW GUESS*)*)
(*(**)
(*SU3OrthogRen1[{QHL[[2]],QHL[[3]],QHL[[4]]},{QL[[2]],QL[[3]],QL[[4]]},{QHR[[1]],QHR[[2]],QHR[[3]]},{QR[[1]],QR[[2]],QR[[3]]},rep1]*SU2OrthogRen[QHL[[4]],Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),(rep1[[1]]-rep1[[4]])/2,(rep2[[1]]-rep2[[4]])/2]*)
(**)*)
(**)
(**)
(**)
(**)
(**)
(**)
(**)
(**)
(**)
(**)
(**)
(*(*OLD GUESS*)*)
(*(*SU3OrthogRen[{QHL[[2]],QHL[[3]],QHL[[4]]},{QL[[2]],QL[[3]],QL[[4]]},{QHR[[2]],QHR[[3]],QHR[[4]]},{QR[[2]],QR[[3]],QR[[4]]},rep1,rep2]*SU2OrthogRen1[QHL[[4]],QR[[2]],(rep1[[1]]-rep1[[4]])/2];*)
(**)*)
(**)
(*(*SU3OrthogRen[{QHL[[2]],QHL[[3]],QHL[[4]]},{QL[[2]],QL[[3]],QL[[4]]},{QHR[[2]],QHR[[3]],QHR[[4]]},{QR[[2]],QR[[3]],QR[[4]]},rep1,rep2]*SU3OrthogRen1[{QHL[[3]],QHL[[4]],QHL[[1]]},{QL[[3]],QL[[4]],QL[[1]]},{QHR[[2]],QHR[[3]],QHR[[4]]},{QR[[2]],QR[[3]],QR[[4]]},rep1];*)
(**)*)
(**)
(*(*SU2OrthogRen1[QHL[[2]],QR[[1]],(rep1[[1]]-rep1[[4]])/2]*SU2OrthogRen[QHL[[3]],QR[[1]],(rep1[[1]]-rep1[[4]])/2,(rep2[[1]]-rep2[[4]])/2];*)*)


(* ::Subsubsection::Closed:: *)
(*Definition OrthogRen22B*)


(* ::Input:: *)
(**)
(*ClearAll[SU3Normalisation2Sub];*)
(*SU3Normalisation2Sub[QHL_,QL_,QHR_,QR_,rep1_,rep2_]:=SU3Normalisation2Sub[QHL,QL,QHR,QR,rep1,rep2]=(Limit[QHL[[3]],z->0,Direction->"FromAbove"])*(Limit[D[QHL[[3]],z],z->0,Direction->"FromAbove"])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),z->0,Direction->"FromAbove"])**)
(*(Limit[QHL[[3]],z->1,Direction->"FromAbove"])*(Limit[D[QHL[[3]],z],z->1,Direction->"FromAbove"])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep2[[1]]-rep2[[2]]),z->1,Direction->"FromAbove"]);*)
(*ClearAll[SU3OrthogRen2Sub];*)
(*SU32SubRow[QHL_,QL_,QHR_,QR_,rep1_,rep2_][i_,j_,\[Alpha]_,w_]:=SU32SubRow[QHL,QL,QHR,QR,rep1,rep2][i,j,\[Alpha],w]={SU3MeanValue[QHL[[i]]*D[QR[[j]],z]*1/(z-w),\[Alpha],rep1,rep2],SU3MeanValue[QHL[[i]]*QR[[j]]*1/(z-w)^2,\[Alpha],rep1,rep2]};*)
(*SU3OrthogRen2Sub[QHL_,QL_,QHR_,QR_,rep1_,rep2_]:=SU3OrthogRen2Sub[QHL,QL,QHR,QR,rep1,rep2]=(1/SU3Normalisation2Sub[QHL,QL,QHR,QR,rep1,rep2])*Det[{Join[SU32SubRow[QHL,QL,QHR,QR,rep1,rep2][3,1,0,0],SU32SubRow[QHL,QL,QHR,QR,rep1,rep2][3,1,0,1]],*)
(*Join[SU32SubRow[QHL,QL,QHR,QR,rep1,rep2][3,2,0,0],SU32SubRow[QHL,QL,QHR,QR,rep1,rep2][3,2,0,1]],*)
(*Join[SU32SubRow[QHL,QL,QHR,QR,rep1,rep2][3,1,1,0],SU32SubRow[QHL,QL,QHR,QR,rep1,rep2][3,1,1,1]],*)
(*Join[SU32SubRow[QHL,QL,QHR,QR,rep1,rep2][3,2,1,0],SU32SubRow[QHL,QL,QHR,QR,rep1,rep2][3,2,1,1]]}]//Simplify;*)
(**)
(*ClearAll[OrthogRen22B];*)
(*OrthogRen22B[QHL_,QL_,QHR_,QR_,rep1_,rep2_]:=OrthogRen22B[QHL,QL,QHR,QR,rep1,rep2]=SU3OrthogRen2Sub[QHL,QL,QHR,QR,rep1,rep2]*SU2OrthogRen1[QHL[[4]],Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),(rep1[[1]]-rep1[[4]])/2];*)
(**)


(* ::Subsection::Closed:: *)
(*RUN: Modules*)


(* ::Input:: *)
(*ClearAll[Q2Base,QH2Base,Q1Base,QH1Base];*)
(*ClearAll[rep1T,rep2T,repAllT,Weight,QT\[Phi],QHT\[Phi],QT1\[Phi],QHT1\[Phi],M,W,MaxM,MaxW,KT,T,REig\[Phi],LEig\[Phi]];*)
(*ClearAll[\[Phi]1T,\[Phi]1,\[Phi]2T,\[Phi]2,\[Phi]3T,\[Phi]3,\[Phi]4T,\[Phi]4,k1,k2,k3,k4];*)
(**)
(*(*Create Q base*)*)
(*Q2Base[rep1T_,rep2T_,repAllT_,Weight_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_][M_]:=Q2Base[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][M]=Table[{QBase[1][rep1T,rep2T,repAllT[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],QBase[2][rep1T,rep2T,repAllT[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],*)
(*QBase[3][rep1T,rep2T,repAllT[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],QBase[4][rep1T,rep2T,repAllT[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]]},{T,1,Length[QBase[1][rep1T,rep2T,repAllT[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T]]}];*)
(**)
(*QH2Base[rep1T_,rep2T_,repAllT_,Weight_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_][M_]:=QH2Base[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][M]=Table[{QHBase[1][rep1T,rep2T,repAllT[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],QHBase[2][rep1T,rep2T,repAllT[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],*)
(*QHBase[3][rep1T,rep2T,repAllT[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],QHBase[4][rep1T,rep2T,repAllT[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]]},{T,1,Length[QHBase[1][rep1T,rep2T,repAllT[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T]]}];*)
(**)
(*Q1Base[rep1T_,rep2T_,repAllT_,Weight_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_][W_]:=Q1Base[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][W]=Table[{QBase1[1][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],QBase1[2][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],*)
(*QBase1[3][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],QBase1[4][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]]},{T,1,Length[QBase1[1][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T]]}];*)
(**)
(*QH1Base[rep1T_,rep2T_,repAllT_,Weight_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_][W_]:=QH1Base[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][W]=Table[{QHBase1[1][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],QHBase1[2][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],*)
(*QHBase1[3][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],QHBase1[4][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]]},{T,1,Length[QHBase1[1][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T]]}];*)
(**)
(*(*Create \[CapitalPsi] base*)*)
(*RBase[rep1T_,rep2T_,repAllT_,Weight_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_,\[Alpha]\[Phi]13_,\[Alpha]\[Phi]14_,\[Alpha]\[Phi]23_][M_]:=RBase[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23][M]=ReEigenR[rep1T,rep2T,repAllT[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23]//Simplify;*)
(**)
(*LBase[rep1T_,rep2T_,repAllT_,Weight_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_,\[Alpha]\[Phi]13_,\[Alpha]\[Phi]14_,\[Alpha]\[Phi]23_][M_]:=LBase[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23][M]=ReEigenL[rep1T,rep2T,repAllT[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23]//Simplify;*)
(**)
(*(*Create \[CapitalPsi] base*)*)
(*RBase1P[rep1T_,rep2T_,repAllT_,Weight_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_,\[Alpha]\[Phi]13_,\[Alpha]\[Phi]14_,\[Alpha]\[Phi]23_][M_]:=RBase1P[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23][M]=ReEigenR1[rep1T,Weight[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23]//Simplify;*)
(**)
(*LBase1P[rep1T_,rep2T_,repAllT_,Weight_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_,\[Alpha]\[Phi]13_,\[Alpha]\[Phi]14_,\[Alpha]\[Phi]23_][M_]:=LBase1P[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23][M]=ReEigenL1[rep1T,Weight[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23]//Simplify;*)
(**)
(**)


(* ::Text:: *)
(*Type of rep:*)
(*{1,0,0,0}: 1*)
(*{1,1,0,0}: 2*)
(*{1,1,1,0}:3*)


(* ::Input:: *)
(*ClearAll[TESTK22,TESTK22B];*)
(*(*Using K22*)*)
(*TESTK22[rep1T_,rep2T_,repAllT_,Weight_,MaxM_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_,\[Alpha]\[Phi]13_,\[Alpha]\[Phi]14_,\[Alpha]\[Phi]23_,k1_,k2_,k3_,k4_][L_]:=TESTK22[rep1T,rep2T,repAllT,Weight,MaxM,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4][L]=*)
(*Table[*)
(*((K22[rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4] . RBase[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23][L][[U]]//Flatten)[[1]])-\[Mu]R*((OrthogRen22[QHVac1[rep1T,k1,k2,k3,k4],QVac1[rep1T,k1,k2,k3,k4],QH2Base[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][L][[U]],Q2Base[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][L][[U]],rep1T,rep2T])/(OrthogRen22[QHVac1[rep1T,k1,k2,k3,k4],QVac1[rep1T,k1,k2,k3,k4],QHVac2[rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T],QVac2[rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T],rep1T,rep2T])),{U,1,Dimensions[RBase[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23][L]][[1]]}];*)
(**)
(*TESTK22A[rep1T_,rep2T_,repAllT_,Weight_,MaxM_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_,\[Alpha]\[Phi]13_,\[Alpha]\[Phi]14_,\[Alpha]\[Phi]23_,k1_,k2_,k3_,k4_][L_]:=TESTK22A[rep1T,rep2T,repAllT,Weight,MaxM,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4][L]=*)
(*Table[*)
(*((K22[rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4] . RBase[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23][L][[U]]//Flatten)[[1]])-\[Mu]R*((OrthogRen22A[QHVac1[rep1T,k1,k2,k3,k4],QVac1[rep1T,k1,k2,k3,k4],QH2Base[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][L][[U]],Q2Base[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][L][[U]],rep1T,rep2T])/(OrthogRen22A[QHVac1[rep1T,k1,k2,k3,k4],QVac1[rep1T,k1,k2,k3,k4],QHVac2[rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T],QVac2[rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T],rep1T,rep2T])),{U,1,Dimensions[RBase[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23][L]][[1]]}];*)
(**)
(*TESTK22B[rep1T_,rep2T_,repAllT_,Weight_,MaxM_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_,\[Alpha]\[Phi]13_,\[Alpha]\[Phi]14_,\[Alpha]\[Phi]23_,k1_,k2_,k3_,k4_][L_]:=TESTK22B[rep1T,rep2T,repAllT,Weight,MaxM,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4][L]=*)
(*Table[*)
(*((K22B[rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4] . RBase[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23][L][[U]]//Flatten)[[1]])-\[Mu]R*((OrthogRen22B[QHVac1[rep1T,k1,k2,k3,k4],QVac1[rep1T,k1,k2,k3,k4],QH2Base[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][L][[U]],Q2Base[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][L][[U]],rep1T,rep2T])/(OrthogRen22B[QHVac1[rep1T,k1,k2,k3,k4],QVac1[rep1T,k1,k2,k3,k4],QHVac2[rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T],QVac2[rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T],rep1T,rep2T])),{U,1,Dimensions[RBase[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23][L]][[1]]}];*)
(**)


(* ::Section:: *)
(*Find K vector*)


(* ::Subsection:: *)
(*Test {1,1,0,0}x{1,1,0,0}, ONLY SU(2) scalar product*)


(* ::Subsubsection::Closed:: *)
(*{1,1,0,0}x{1,1,0,0}, SU(2) 2-sites*)


(* ::Input:: *)
(*rep1T={1,1,0,0};*)
(*rep2T={1,1,0,0};*)
(*(*FindGoodWeight2[rep1T,rep2T]*)*)
(*repAllT=FindWeight2[rep1T,rep2T];*)
(*(*Permutations[rep1T+rep2T]*)*)
(*MaxM=Length[repAllT];*)
(*Weight=FindGoodWeight[rep1T];*)
(*MaxW=Length[FindGoodWeight[rep1T]];*)
(**)
(*\[Phi]1T=1;*)
(*\[Phi]2T=1/2;*)
(*\[Phi]3T=1/3;*)
(*\[Phi]4T=1/5;*)
(**)
(*(**)
(*\[Phi]1T=\[Phi]1;*)
(*\[Phi]2T=\[Phi]2;*)
(*\[Phi]3T=\[Phi]3;*)
(*\[Phi]4T=\[Phi]4;*)
(**)*)


(* ::Text:: *)
(*TEST K22 crates the equations:*)
(*<K|Subscript[Psi, A]>-((<(e^(k z)),Subscript[Q, A]>)/(<(e^(k z)),Subscript[Q, Vac]>)),*)
(*where <\[CenterDot],\[CenterDot]> is <(e^kz),\!\( *)
(*\*SubscriptBox[\(Q\), \(12, A\)]*)
(*\*SubscriptBox[\(>\), \(SU(2), \ 2 - sites\)]\)*)


(* ::Input:: *)
(*ClearAll[EQS];*)
(*EQS=Flatten[Table[TESTK22[rep1T,rep2T,repAllT,Weight,MaxM,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4][i],{i,1,MaxM}]]//Simplify*)


(* ::Input:: *)
(*Length[EQS]*)


(* ::Text:: *)
(*We get only the trivial solution:*)


(* ::Input:: *)
(*Solve[EQS==0,{p11,p33,p44,p12,p34,P11,P22,P33,P44,P12,P13,P14,P23,P24,P34,\[Mu]R}]//Simplify*)


(* ::Subsubsection::Closed:: *)
(*{1,1,0,0}x{1,1,0,0}, SU(2) 1-site*)


(* ::Input:: *)
(*rep1T={1,1,0,0};*)
(*rep2T={1,1,0,0};*)
(*(*FindGoodWeight2[rep1T,rep2T]*)*)
(*repAllT=FindWeight2[rep1T,rep2T];*)
(*(*Permutations[rep1T+rep2T]*)*)
(*MaxM=Length[repAllT];*)
(*Weight=FindGoodWeight[rep1T];*)
(*MaxW=Length[FindGoodWeight[rep1T]];*)
(**)
(*\[Phi]1T=1;*)
(*\[Phi]2T=1/2;*)
(*\[Phi]3T=1/3;*)
(*\[Phi]4T=1/5;*)
(**)
(*(**)
(*\[Phi]1T=\[Phi]1;*)
(*\[Phi]2T=\[Phi]2;*)
(*\[Phi]3T=\[Phi]3;*)
(*\[Phi]4T=\[Phi]4;*)
(**)*)


(* ::Text:: *)
(*TEST K22 crates the equations:*)
(*<K|Subscript[Psi, A]>-((<(e^(k z)),Subscript[Q, A]>)/(<(e^(k z)),Subscript[Q, Vac]>)),*)
(*where <\[CenterDot],\[CenterDot]> is <(e^kz),\!\( *)
(*\*SubscriptBox[\(Q\), \(12, A\)]*)
(*\*SubscriptBox[\(>\), \(SU(2), \ 1 - sites\)]\)*)


(* ::Input:: *)
(*ClearAll[EQS];*)
(*EQS=Flatten[Table[TESTK22A[rep1T,rep2T,repAllT,Weight,MaxM,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4][i],{i,1,MaxM}]]//Simplify*)


(* ::Input:: *)
(*Length[EQS]*)


(* ::Text:: *)
(*We get only the trivial solution:*)


(* ::Input:: *)
(*Solve[EQS==0,{p11,p33,p44,p12,p34,P11,P22,P33,P44,P12,P13,P14,P23,P24,P34,\[Mu]R}]//Simplify*)


(* ::Subsection::Closed:: *)
(*Test {1,1,0,0}x{1,1,0,0}, SU(3) and SU(2) scalar product*)


(* ::Subsubsection:: *)
(*{1,1,0,0}x{1,1,0,0}*)


(* ::Input:: *)
(*rep1T={1,1,0,0};*)
(*rep2T={1,1,0,0};*)
(*(*FindGoodWeight2[rep1T,rep2T]*)*)
(*repAllT=FindWeight2[rep1T,rep2T];*)
(*(*Permutations[rep1T+rep2T]*)*)
(*MaxM=Length[repAllT];*)
(*Weight=FindGoodWeight[rep1T];*)
(*MaxW=Length[FindGoodWeight[rep1T]];*)
(**)
(*\[Phi]1T=1;*)
(*\[Phi]2T=1/2;*)
(*\[Phi]3T=1/3;*)
(*\[Phi]4T=1/5;*)
(**)
(*(**)
(*\[Phi]1T=\[Phi]1;*)
(*\[Phi]2T=\[Phi]2;*)
(*\[Phi]3T=\[Phi]3;*)
(*\[Phi]4T=\[Phi]4;*)
(**)*)


(* ::Text:: *)
(*TEST K22 crates the equations:*)
(*<K|Subscript[Psi, A]>-((<(e^(k z)),Subscript[Q, A]>)/(<(e^(k z)),Subscript[Q, Vac]>)),*)
(*where <\[CenterDot],\[CenterDot]> is <(e^kz),\!\( *)
(*\*SubscriptBox[\(Q\), \(12, A\)]*)
(*\*SubscriptBox[\(>\), \(SU(2), \ 1 - sites\)]\) <(e^(Subscript[k, 2] z)),\!\(Q*)
(*\*SubscriptBox[\(>\), \(\(SU(3)\)\ 2 - site\)]\)*)


(* ::Input:: *)
(*ClearAll[EQS];*)
(*EQS=Flatten[Table[TESTK22B[rep1T,rep2T,repAllT,Weight,MaxM,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4][i],{i,1,MaxM}]]//Simplify*)


(* ::Input:: *)
(*Length[EQS]*)


(* ::Text:: *)
(*We get only the trivial solution:*)


(* ::Input:: *)
(*Solve[EQS==0,{p11,p33,p44,p12,p34,P11,P22,P33,P44,P12,P13,P14,P23,P24,P34,\[Mu]R}]//Simplify*)


(* ::Title:: *)
(*K vector 1-site*)


(* ::Section::Closed:: *)
(*RUN: Definition K vector*)


(* ::Subsection::Closed:: *)
(*<K|*)


(* ::Input:: *)
(*ClearAll[K1];*)
(**)
(*K1[rep1_,p12_,p13_,p14_,p23_,p24_,p34_,p11_,p22_,p33_,p44_]:=K1[rep1,p12,p13,p14,p23,p24,p34,p11,p22,p33,p44]=UPL1[rep1] . MatrixExp[Log[p11]*EM[1,1,rep1]] . MatrixExp[Log[p22]*EM[2,2,rep1]] . MatrixExp[Log[p33]*EM[3,3,rep1]] . MatrixExp[Log[p44]*EM[4,4,rep1]] . MatrixExp[p12(EM[1,2,rep1])] . MatrixExp[p34(EM[3,4,rep1])] . MatrixExp[p13(EM[1,3,rep1])] . MatrixExp[p14(EM[1,4,rep1])] . MatrixExp[p23(EM[2,3,rep1])] . MatrixExp[p24(EM[2,4,rep1])];*)


(* ::Section::Closed:: *)
(*RUN: Modules for testing*)


(* ::Subsection::Closed:: *)
(*Definition SU(2) scalar product*)


(* ::Input:: *)
(*ClearAll[SU2MeanValue];*)
(*SU2MeanValue[f_,\[Tau]_,l1_,l2_]:=*)
(*2 \[Pi] I Residue[(f*1/((z-w1)^(2 l1 +1) (z-w2)^(2 l2 +1))),{z,\[Tau]}];*)


(* ::Input:: *)
(*ClearAll[SU2MeanValue1];*)
(*SU2MeanValue1[f_,\[Tau]_,l1_]:=*)
(*2 \[Pi] I Residue[(f*1/(z-w1)^(2 l1+1 ) ),{z,\[Tau]}];*)


(* ::Input:: *)
(*(*2-sites scalar product*)*)
(*ClearAll[SU2OrthogRen];*)
(*SU2OrthogRen[Q1HODGE_,Q2_,l1_,l2_]:=(SU2MeanValue[Q1HODGE*Q2,w1,l1,l2]*SU2MeanValue[Q1HODGE*Q2*z,w2,l1,l2]-SU2MeanValue[Q1HODGE*Q2,w2,l1,l2]*SU2MeanValue[Q1HODGE*Q2*z,w1,l1,l2])/((Limit[Q1HODGE,z->w1])*(Limit[Q1HODGE,z->w2])*(Limit[Q2,z->w1])*(Limit[Q2,z->w2]))//Simplify;*)


(* ::Input:: *)
(*(*1-site scalar product*)*)
(*ClearAll[SU2OrthogRen1];*)
(*SU2OrthogRen1[Q1HODGE_,Q2_,l1_]:=(SU2MeanValue1[Q1HODGE*Q2,w1,l1])/((Limit[Q1HODGE,z->w1])*(Limit[Q2,z->w1]))//Simplify;*)


(* ::Subsection::Closed:: *)
(*Definition SU(3) scalar product*)


(* ::Input:: *)
(*ClearAll[SU3MeanValue];*)
(*SU3MeanValue[f_,\[Tau]_,rep1_,rep2_]:=*)
(*2 \[Pi] I Residue[(f*1/((z-w1)^(rep1[[1]]-rep1[[3]]) (z-w2)^(rep2[[1]]-rep2[[3]]))),{z,\[Tau]}];*)


(* ::Input:: *)
(*ClearAll[SU3MeanValue1];*)
(*SU3MeanValue1[f_,\[Tau]_,rep1_]:=*)
(*2 \[Pi] I Residue[(f*1/(z-w1)^(rep1[[1]]-rep1[[3]])),{z,\[Tau]}];;*)


(* ::Subsubsection::Closed:: *)
(*2-sites scalar product*)


(* ::Text:: *)
(*Input:*)
(*Q={Subscript[Q, 1],Subscript[Q, 2],Subscript[Q, 3]}*)
(*QH={Subscript[QH, 1],Subscript[QH, 2],Subscript[QH, 3]}*)


(* ::Input:: *)
(*ClearAll[SU3Normalisation];*)
(*SU3Normalisation[QHL_,QL_,QHR_,QR_]:=SU3Normalisation[QHL,QL,QHR,QR]=(Limit[QHL[[1]],z->0,Direction->"FromAbove"])*(Limit[QL[[3]],z->0,Direction->"FromAbove"])*(Limit[QHR[[3]],z->0,Direction->"FromAbove"])*(Limit[QR[[1]],z->0,Direction->"FromAbove"])**)
(*(Limit[QHL[[1]],z->1,Direction->"FromBelow"])*(Limit[QL[[3]],z->1,Direction->"FromBelow"])*(Limit[QHR[[3]],z->1,Direction->"FromBelow"])*(Limit[QR[[1]],z->1,Direction->"FromBelow"]);*)


(* ::Input:: *)
(*ClearAll[SU3OrthogRen];*)
(*SU3OrthogRen[QHL_,QL_,QHR_,QR_,rep1_,rep2_]:=SU3OrthogRen[QHL,QL,QHR,QR,rep1,rep2]=(1/SU3Normalisation[QHL,QL,QHR,QR])*Det[{{SU3MeanValue[QHL[[1]]*D[QR[[1]],z]*1/z,0,rep1,rep2],SU3MeanValue[QHL[[1]]*D[QR[[1]],z]*1/(z-1),0,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[1]]*1/z^2,0,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[1]]*1/(z-1)^2,0,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[1]]*1/z,0,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[1]]*1/(z-1),0,rep1,rep2]},*)
(**)
(*{SU3MeanValue[QHL[[1]]*D[QR[[2]],z]*1/z,0,rep1,rep2],SU3MeanValue[QHL[[1]]*D[QR[[2]],z]*1/(z-1),0,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[2]]*1/z^2,0,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[2]]*1/(z-1)^2,0,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[2]]*1/z,0,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[2]]*1/(z-1),0,rep1,rep2]},*)
(**)
(*{SU3MeanValue[QHL[[2]]*D[QR[[1]],z]*1/z,0,rep1,rep2],SU3MeanValue[QHL[[2]]*D[QR[[1]],z]*1/(z-1),0,rep1,rep2],SU3MeanValue[QHL[[2]]*QR[[1]]*1/z^2,0,rep1,rep2],SU3MeanValue[QHL[[2]]*QR[[1]]*1/(z-1)^2,0,rep1,rep2],SU3MeanValue[QHL[[2]]*QR[[1]]*1/z,0,rep1,rep2],SU3MeanValue[QHL[[2]]*QR[[1]]*1/(z-1),0,rep1,rep2]},*)
(**)
(**)
(**)
(**)
(**)
(*{SU3MeanValue[QHL[[1]]*D[QR[[1]],z]*1/z,1,rep1,rep2],SU3MeanValue[QHL[[1]]*D[QR[[1]],z]*1/(z-1),1,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[1]]*1/z^2,1,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[1]]*1/(z-1)^2,1,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[1]]*1/z,1,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[1]]*1/(z-1),1,rep1,rep2]},*)
(**)
(*{SU3MeanValue[QHL[[1]]*D[QR[[2]],z]*1/z,1,rep1,rep2],SU3MeanValue[QHL[[1]]*D[QR[[2]],z]*1/(z-1),1,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[2]]*1/z^2,1,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[2]]*1/(z-1)^2,1,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[2]]*1/z,1,rep1,rep2],SU3MeanValue[QHL[[1]]*QR[[2]]*1/(z-1),1,rep1,rep2]},*)
(**)
(*{SU3MeanValue[QHL[[2]]*D[QR[[1]],z]*1/z,1,rep1,rep2],SU3MeanValue[QHL[[2]]*D[QR[[1]],z]*1/(z-1),1,rep1,rep2],SU3MeanValue[QHL[[2]]*QR[[1]]*1/z^2,1,rep1,rep2],SU3MeanValue[QHL[[2]]*QR[[1]]*1/(z-1)^2,1,rep1,rep2],SU3MeanValue[QHL[[2]]*QR[[1]]*1/z,1,rep1,rep2],SU3MeanValue[QHL[[2]]*QR[[1]]*1/(z-1),1,rep1,rep2]}}]//Simplify;*)


(* ::Subsubsection::Closed:: *)
(*1-site scalar product*)


(* ::Input:: *)
(*ClearAll[SU3Normalisation1];*)
(*SU3Normalisation1[QHL_,QL_,QHR_,QR_]:=SU3Normalisation1[QHL,QL,QHR,QR]=(Limit[QHL[[1]],z->0,Direction->"FromAbove"])*(Limit[QL[[3]],z->0,Direction->"FromAbove"])*(Limit[QHR[[3]],z->0,Direction->"FromAbove"])*(Limit[QR[[1]],z->0,Direction->"FromAbove"]);*)


(* ::Input:: *)
(*ClearAll[SU3OrthogRen1];*)
(*SU3OrthogRen1[QHL_,QL_,QHR_,QR_,rep1_]:=SU3OrthogRen1[QHL,QL,QHR,QR,rep1]=(1/SU3Normalisation1[QHL,QL,QHR,QR])*Det[{{SU3MeanValue1[QHL[[1]]*D[QR[[1]],z]*1/z,0,rep1],SU3MeanValue1[QHL[[1]]*QR[[1]]*1/z^2,0,rep1],SU3MeanValue1[QHL[[1]]*QR[[1]]*1/z,0,rep1]},*)
(**)
(*{SU3MeanValue1[QHL[[1]]*D[QR[[2]],z]*1/z,0,rep1],SU3MeanValue1[QHL[[1]]*QR[[2]]*1/z^2,0,rep1],SU3MeanValue1[QHL[[1]]*QR[[2]]*1/z,0,rep1]},*)
(**)
(*{SU3MeanValue1[QHL[[2]]*D[QR[[1]],z]*1/z,0,rep1],SU3MeanValue1[QHL[[2]]*QR[[1]]*1/z^2,0,rep1],SU3MeanValue1[QHL[[2]]*QR[[1]]*1/z,0,rep1]}}]//Simplify;*)


(* ::Subsection::Closed:: *)
(*Definition scalar product for {p,p,0,0}x{q,q,0,0}*)


(* ::Subsubsection::Closed:: *)
(*Using su(3) scalar product*)


(* ::Input:: *)
(*ClearAll[SU3Normalisation1Sub];*)
(*SU3Normalisation1Sub[QHL_,QL_,QHR_,QR_,rep1_]:=SU3Normalisation1Sub[QHL,QL,QHR,QR,rep1]=(Limit[QHL[[3]],z->0,Direction->"FromAbove"])*(Limit[D[QHL[[3]],z],z->0,Direction->"FromAbove"])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),z->0,Direction->"FromAbove"]);*)
(*(*(Limit[QHL[[1]],z->0,Direction->"FromAbove"])*(Limit[QL[[3]],z->0,Direction->"FromAbove"])*(Limit[QHR[[3]],z->0,Direction->"FromAbove"])*(Limit[QR[[1]],z->0,Direction->"FromAbove"]);*)*)


(* ::Input:: *)
(*ClearAll[SU3OrthogRen1Sub];*)
(*SU3OrthogRen1Sub[QHL_,QL_,QHR_,QR_,rep1_]:=SU3OrthogRen1Sub[QHL,QL,QHR,QR,rep1]=(1/SU3Normalisation1Sub[QHL,QL,QHR,QR,rep1])*Det[{{SU3MeanValue1[QHL[[3]]*D[QR[[1]],z]*1/z,0,rep1],SU3MeanValue1[QHL[[3]]*QR[[1]]*1/z^2,0,rep1]},*)
(**)
(*{SU3MeanValue1[QHL[[3]]*D[QR[[2]],z]*1/z,0,rep1],SU3MeanValue1[QHL[[3]]*QR[[2]]*1/z^2,0,rep1]}}]//Simplify;*)
(**)
(*(*Det[{{SU3MeanValue1[QHL[[4]]*QR[[1]]*1/z^2,0,rep1],SU3MeanValue1[QHL[[4]]*QR[[1]]*1/z,0,rep1]},*)
(**)
(*{SU3MeanValue1[QHL[[4]]*QR[[2]]*1/z^2,0,rep1],SU3MeanValue1[QHL[[4]]*QR[[2]]*1/z,0,rep1]}}]*)*)
(**)
(*(*Full SU(3) column*)*)
(*(*Det[{{SU3MeanValue1[QHL[[1]]*D[QR[[1]],z]*1/z,0,rep1],SU3MeanValue1[QHL[[1]]*QR[[1]]*1/z^2,0,rep1],SU3MeanValue1[QHL[[1]]*QR[[1]]*1/z,0,rep1]},*)
(**)
(*{SU3MeanValue1[QHL[[1]]*D[QR[[2]],z]*1/z,0,rep1],SU3MeanValue1[QHL[[1]]*QR[[2]]*1/z^2,0,rep1],SU3MeanValue1[QHL[[1]]*QR[[2]]*1/z,0,rep1]}}]//Simplify;*)*)


(* ::Input:: *)
(*ClearAll[OrthogRen1Site22];*)
(*OrthogRen1Site22[QHL_,QL_,QHR_,QR_,rep1_]:=OrthogRen1Site22[QHL,QL,QHR,QR,rep1]=SU2OrthogRen1[QHL[[4]],Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),(rep1[[1]]-rep1[[4]])/2];*)
(**)
(*ClearAll[OrthogRen1Site22B]*)
(*OrthogRen1Site22B[QHL_,QL_,QHR_,QR_,rep1_]:=OrthogRen1Site22B[QHL,QL,QHR,QR,rep1]=SU3OrthogRen1Sub[QHL,QL,QHR,QR,rep1]*SU2OrthogRen1[QHL[[4]],Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),(rep1[[1]]-rep1[[4]])/2];*)
(**)
(*(*What WORKS for 1-site*)*)
(*(**)
(*SU2OrthogRen1[QHL[[4]],Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),(rep1[[1]]-rep1[[4]])/2]*)
(**)*)
(**)
(**)
(*(*Main Guess*)*)
(*(**)
(*SU3OrthogRen1Sub[QHL,QL,QHR,QR,rep1]*SU2OrthogRen[QHL[[4]],Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),(rep1[[1]]-rep1[[4]])/2,(rep2[[1]]-rep2[[4]])/2];*)
(**)*)
(**)
(**)
(**)
(*(*My guess*)*)
(*(**)
(*SU2OrthogRen1[QHL[[4]],Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),(rep1[[1]]-rep1[[4]])/2]*SU2OrthogRen[QHL[[3]],Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),(rep1[[1]]-rep1[[4]])/2,(rep2[[1]]-rep2[[4]])/2];*)
(**)*)
(**)
(**)
(**)
(**)
(**)
(*(*NEW GUESSES*)*)
(*(*First NEW GUESS*)*)
(*(**)
(*SU3OrthogRen[{QHL[[2]],QHL[[3]],QHL[[4]]},{QL[[2]],QL[[3]],QL[[4]]},{QHR[[1]],QHR[[2]],QHR[[3]]},{QR[[1]],QR[[2]],QR[[3]]},rep1,rep2]*SU2OrthogRen1[QHL[[4]],Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),(rep1[[1]]-rep1[[4]])/2];*)
(**)*)
(**)
(*(*Second NEW GUESS*)*)
(*(**)
(*SU3OrthogRen1[{QHL[[2]],QHL[[3]],QHL[[4]]},{QL[[2]],QL[[3]],QL[[4]]},{QHR[[1]],QHR[[2]],QHR[[3]]},{QR[[1]],QR[[2]],QR[[3]]},rep1]*SU2OrthogRen[QHL[[4]],Wronskian[{QR[[1]],QR[[2]]},z]*1/z^(rep1[[1]]-rep1[[2]]),(rep1[[1]]-rep1[[4]])/2,(rep2[[1]]-rep2[[4]])/2]*)
(**)*)
(**)
(**)
(**)
(**)
(**)
(**)
(**)
(**)
(**)
(**)
(**)
(*(*OLD GUESS*)*)
(*(*SU3OrthogRen[{QHL[[2]],QHL[[3]],QHL[[4]]},{QL[[2]],QL[[3]],QL[[4]]},{QHR[[2]],QHR[[3]],QHR[[4]]},{QR[[2]],QR[[3]],QR[[4]]},rep1,rep2]*SU2OrthogRen1[QHL[[4]],QR[[2]],(rep1[[1]]-rep1[[4]])/2];*)
(**)*)
(**)
(*(*SU3OrthogRen[{QHL[[2]],QHL[[3]],QHL[[4]]},{QL[[2]],QL[[3]],QL[[4]]},{QHR[[2]],QHR[[3]],QHR[[4]]},{QR[[2]],QR[[3]],QR[[4]]},rep1,rep2]*SU3OrthogRen1[{QHL[[3]],QHL[[4]],QHL[[1]]},{QL[[3]],QL[[4]],QL[[1]]},{QHR[[2]],QHR[[3]],QHR[[4]]},{QR[[2]],QR[[3]],QR[[4]]},rep1];*)
(**)*)
(**)
(*(*SU2OrthogRen1[QHL[[2]],QR[[1]],(rep1[[1]]-rep1[[4]])/2]*SU2OrthogRen[QHL[[3]],QR[[1]],(rep1[[1]]-rep1[[4]])/2,(rep2[[1]]-rep2[[4]])/2];*)*)


(* ::Subsection::Closed:: *)
(*Modules*)


(* ::Input:: *)
(*ClearAll[Q1Base,QH1Base];*)
(*ClearAll[rep1T,Weight,QT\[Phi],QHT\[Phi],QT1\[Phi],QHT1\[Phi],M,W,MaxM,MaxW,KT,T,REig\[Phi],LEig\[Phi]];*)
(*ClearAll[\[Phi]1T,\[Phi]1,\[Phi]2T,\[Phi]2,\[Phi]3T,\[Phi]3,\[Phi]4T,\[Phi]4,k1,k2,k3,k4];*)
(**)
(*Q1Base[rep1T_,Weight_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_][W_]:=Q1Base[rep1T,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][W]=Table[{QBase1[1][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],QBase1[2][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],*)
(*QBase1[3][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],QBase1[4][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]]},{T,1,Length[QBase1[1][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T]]}];*)
(**)
(*QH1Base[rep1T_,Weight_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_][W_]:=QH1Base[rep1T,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][W]=Table[{QHBase1[1][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],QHBase1[2][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],*)
(*QHBase1[3][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],QHBase1[4][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]]},{T,1,Length[QHBase1[1][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T]]}];*)
(**)
(*(*Create \[CapitalPsi] base*)*)
(*RBase1P[rep1T_,Weight_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_,\[Alpha]\[Phi]13_,\[Alpha]\[Phi]14_,\[Alpha]\[Phi]23_][M_]:=RBase1P[rep1T,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23][M]=ReEigenR1[rep1T,Weight[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23]//Simplify;*)
(**)
(*LBase1P[rep1T_,Weight_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_,\[Alpha]\[Phi]13_,\[Alpha]\[Phi]14_,\[Alpha]\[Phi]23_][M_]:=LBase1P[rep1T,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23][M]=ReEigenL1[rep1T,Weight[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23]//Simplify;*)
(**)
(**)


(* ::Text:: *)
(*Type of rep:*)
(*{1,0,0,0}: 1*)
(*{1,1,0,0}: 2*)
(*{1,1,1,0}:3*)


(* ::Input:: *)
(*ClearAll[TESTK1Site22];*)
(*(*With 1-site*)*)
(*TESTK1Site22[rep1T_,Weight_,MaxM_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_,\[Alpha]\[Phi]13_,\[Alpha]\[Phi]14_,\[Alpha]\[Phi]23_,k1_,k2_,k3_,k4_][L_]:=TESTK1Site22[rep1T,Weight,MaxM,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4][L]=*)
(*Table[*)
(*((K1[rep1T,p12,p13,p14,p23,p24,p34,p11,p22,p33,p44] . RBase1P[rep1T,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23][L][[U]]//Flatten)[[1]])-1*((OrthogRen1Site22[{Exp[k1 z],Exp[k2 z],Exp[k3 z],Exp[k4 z]},{Exp[-k1 z],Exp[-k2 z],Exp[-k3 z],Exp[-k4 z]},QH1Base[rep1T,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][L][[U]],Q1Base[rep1T,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][L][[U]],rep1T])/(OrthogRen1Site22[{Exp[k1 z],Exp[k2 z],Exp[k3 z],Exp[k4 z]},{Exp[-k1 z],Exp[-k2 z],Exp[-k3 z],Exp[-k4 z]},QH1Base[rep1T,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][1][[1]],Q1Base[rep1T,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][1][[1]],rep1T])),{U,1,Dimensions[RBase1P[rep1T,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23][L]][[1]]}];*)
(**)
(*ClearAll[TESTK1Site22B]*)
(*TESTK1Site22B[rep1T_,Weight_,MaxM_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_,\[Alpha]\[Phi]13_,\[Alpha]\[Phi]14_,\[Alpha]\[Phi]23_,k1_,k2_,k3_,k4_][L_]:=TESTK1Site22B[rep1T,Weight,MaxM,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4][L]=*)
(*Table[*)
(*((K1[rep1T,p12,p13,p14,p23,p24,p34,p11,p22,p33,p44] . RBase1P[rep1T,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23][L][[U]]//Flatten)[[1]])-1*((OrthogRen1Site22B[{Exp[k1 z],Exp[k2 z],Exp[k3 z],Exp[k4 z]},{Exp[-k1 z],Exp[-k2 z],Exp[-k3 z],Exp[-k4 z]},QH1Base[rep1T,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][L][[U]],Q1Base[rep1T,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][L][[U]],rep1T])/(OrthogRen1Site22B[{Exp[k1 z],Exp[k2 z],Exp[k3 z],Exp[k4 z]},{Exp[-k1 z],Exp[-k2 z],Exp[-k3 z],Exp[-k4 z]},QH1Base[rep1T,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][1][[1]],Q1Base[rep1T,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][1][[1]],rep1T])),{U,1,Dimensions[RBase1P[rep1T,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23][L]][[1]]}];*)
(**)
(**)
(*(*With QHVac it takes too much time*)*)
(**)
(*(**)
(*TESTK1Site22[rep1T_,Weight_,MaxM_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_,\[Alpha]\[Phi]12_,\[Alpha]\[Phi]14_,\[Alpha]\[Phi]34_,k1_,k2_,k3_,k4_][L_]:=TESTK1Site22[rep1T,Weight,MaxM,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]12,\[Alpha]\[Phi]14,\[Alpha]\[Phi]34,k1,k2,k3,k4][L]=*)
(*Table[*)
(*((K1[rep1T,p23,p13,p12,p24,p14,p34,p11,p22,p33,p44].RBase1P[rep1T,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]12,\[Alpha]\[Phi]14,\[Alpha]\[Phi]34][L][[U]]//Flatten)[[1]])-1*((OrthogRen1Site22[QHVac1[rep1T,k1,k2,k3,k4],QVac1[rep1T,k1,k2,k3,k4],QH1Base[rep1T,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][L][[U]],Q1Base[rep1T,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][L][[U]],rep1T])/(OrthogRen1Site22[QHVac1[rep1T,k1,k2,k3,k4],QVac1[rep1T,k1,k2,k3,k4],QH1Base[rep1T,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][1][[1]],Q1Base[rep1T,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][1][[1]],rep1T])),{U,1,Dimensions[RBase1P[rep1T,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]12,\[Alpha]\[Phi]14,\[Alpha]\[Phi]34][L]][[1]]}];*)
(**)*)


(* ::Section:: *)
(*Find K vector*)


(* ::Subsection:: *)
(*Test {1,1,0,0}x{1,1,0,0}, ONLY SU(2) scalar product*)


(* ::Text:: *)
(*Check the equations:*)
(*<K|Subscript[Psi, A]>-((<(e^(k z)),Subscript[Q, A]>)/(<(e^(k z)),Subscript[Q, Vac]>)),*)
(*where <\[CenterDot],\[CenterDot]> is <(e^kz),\!\( *)
(*\*SubscriptBox[\(Q\), \(12, A\)]*)
(*\*SubscriptBox[\(>\), \(SU(2), \ 1 - sites\)]\)*)


(* ::Subsubsection::Closed:: *)
(*{1,1,0,0}*)


(* ::Input:: *)
(*rep1T={1,1,0,0};*)
(*Weight=FindGoodWeight[rep1T];*)
(*MaxW=Length[FindGoodWeight[rep1T]];*)
(*WeightAll=FindWeight [rep1T];*)
(*MaxWAll=Length[FindWeight[rep1T]];*)
(**)
(*(**)
(*\[Phi]1T=1;*)
(*\[Phi]2T=1/2;*)
(*\[Phi]3T=1/3;*)
(*\[Phi]4T=1/5;*)
(**)*)
(**)
(**)
(**)
(*\[Phi]1T=\[Phi]1;*)
(*\[Phi]2T=\[Phi]2;*)
(*\[Phi]3T=\[Phi]3;*)
(*\[Phi]4T=\[Phi]4;*)


(* ::Input:: *)
(*ClearAll[EQS];*)
(*EQS=Flatten[Table[TESTK1Site22[rep1T,WeightAll,MaxWAll,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4][i],{i,1,MaxWAll}]]//Simplify*)


(* ::Input:: *)
(*Solve[EQS==0,{p11,p22,p33,p44,p12,p13,p14,p23,p24,p34}]//Simplify*)


(* ::Input:: *)
(*UPL1[rep1] . MatrixExp[Log[p11]*EM[1,1,rep1]] . MatrixExp[Log[1/p11]*EM[2,2,rep1]] . MatrixExp[Log[p33]*EM[3,3,rep1]] . MatrixExp[Log[p44]*EM[4,4,rep1]] . MatrixExp[p12(EM[1,2,rep1])] . MatrixExp[p34(EM[3,4,rep1])] . MatrixExp[(\[Alpha]\[Phi]13 (k4-\[Phi]2-\[Phi]3))/(k4-\[Phi]1-\[Phi]2) (EM[1,3,rep1])] . MatrixExp[(\[Alpha]\[Phi]14 (k4-\[Phi]2-\[Phi]4))/(k4-\[Phi]1-\[Phi]2) (EM[1,4,rep1])] . MatrixExp[(\[Alpha]\[Phi]23 (k4-\[Phi]1-\[Phi]3))/(k4-\[Phi]1-\[Phi]2) (EM[2,3,rep1])] . MatrixExp[(\[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (\[Phi]2-\[Phi]3) (k4-\[Phi]1-\[Phi]4) (\[Phi]1-\[Phi]4))/(\[Alpha]\[Phi]13 (-k4+\[Phi]1+\[Phi]2) (-\[Phi]1+\[Phi]3) (\[Phi]2-\[Phi]4)) (EM[2,4,rep1])];*)


(* ::Subsubsection::Closed:: *)
(*{2,2,0,0}*)


(* ::Input:: *)
(*rep1T={2,2,0,0};*)
(*Weight=FindGoodWeight[rep1T];*)
(*MaxW=Length[FindGoodWeight[rep1T]];*)
(*WeightAll=FindWeight [rep1T];*)
(*MaxWAll=Length[FindWeight[rep1T]];*)
(**)
(**)
(*\[Phi]1T=1;*)
(*\[Phi]2T=1/2;*)
(*\[Phi]3T=1/3;*)
(*\[Phi]4T=1/5;*)
(**)
(**)
(**)
(*(**)
(*\[Phi]1T=\[Phi]1;*)
(*\[Phi]2T=\[Phi]2;*)
(*\[Phi]3T=\[Phi]3;*)
(*\[Phi]4T=\[Phi]4;*)
(**)*)


(* ::Input:: *)
(*ClearAll[EQS];*)
(*EQS=Flatten[Table[TESTK1Site22[rep1T,WeightAll,Max,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4][i],{i,1,MaxWAll}]]//Simplify*)


(* ::Input:: *)
(*EQS/.{p22->1/p11,p13->(\[Alpha]\[Phi]13 (k4-\[Phi]2T-\[Phi]3T))/(k4-\[Phi]1T-\[Phi]2T),p14->(\[Alpha]\[Phi]14 (k4-\[Phi]2T-\[Phi]4T))/(k4-\[Phi]1T-\[Phi]2T),p23->(\[Alpha]\[Phi]23 (k4-\[Phi]1T-\[Phi]3T))/(k4-\[Phi]1T-\[Phi]2T),p24->(\[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (\[Phi]2T-\[Phi]3T) (k4-\[Phi]1T-\[Phi]4T) (\[Phi]1T-\[Phi]4T))/(\[Alpha]\[Phi]13 (-k4+\[Phi]1T+\[Phi]2T) (-\[Phi]1T+\[Phi]3T) (\[Phi]2T-\[Phi]4T))}//PowerExpand//Simplify*)


(* ::Subsubsection::Closed:: *)
(*{3,3,0,0}*)


(* ::Input:: *)
(*rep1T={3,3,0,0};*)
(*Weight=FindGoodWeight[rep1T];*)
(*MaxW=Length[FindGoodWeight[rep1T]];*)
(*WeightAll=FindWeight [rep1T];*)
(*MaxWAll=Length[FindWeight[rep1T]];*)
(**)
(**)
(*\[Phi]1T=1;*)
(*\[Phi]2T=1/2;*)
(*\[Phi]3T=1/3;*)
(*\[Phi]4T=1/5;*)
(**)
(**)
(**)
(*(**)
(*\[Phi]1T=\[Phi]1;*)
(*\[Phi]2T=\[Phi]2;*)
(*\[Phi]3T=\[Phi]3;*)
(*\[Phi]4T=\[Phi]4;*)
(**)*)


(* ::Input:: *)
(*ClearAll[EQS];*)
(*EQS=Flatten[Table[TESTK1Site22[rep1T,WeightAll,Max,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4][i],{i,1,MaxWAll}]]//Simplify*)


(* ::Input:: *)
(*EQS/.{p22->1/p11,p13->(\[Alpha]\[Phi]13 (k4-\[Phi]2T-\[Phi]3T))/(k4-\[Phi]1T-\[Phi]2T),p14->(\[Alpha]\[Phi]14 (k4-\[Phi]2T-\[Phi]4T))/(k4-\[Phi]1T-\[Phi]2T),p23->(\[Alpha]\[Phi]23 (k4-\[Phi]1T-\[Phi]3T))/(k4-\[Phi]1T-\[Phi]2T),p24->(\[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (\[Phi]2T-\[Phi]3T) (k4-\[Phi]1T-\[Phi]4T) (\[Phi]1T-\[Phi]4T))/(\[Alpha]\[Phi]13 (-k4+\[Phi]1T+\[Phi]2T) (-\[Phi]1T+\[Phi]3T) (\[Phi]2T-\[Phi]4T))}//PowerExpand//Simplify*)


(* ::Subsection:: *)
(*Test {1,1,0,0}x{1,1,0,0} SU(3)xSU(2)*)


(* ::Text:: *)
(*TEST K22 crates the equations:*)
(*<K|Subscript[Psi, A]>-((<(e^(k z)),Subscript[Q, A]>)/(<(e^(k z)),Subscript[Q, Vac]>)),*)
(*where <\[CenterDot],\[CenterDot]> is <(e^kz),\!\( *)
(*\*SubscriptBox[\(Q\), \(12, A\)]*)
(*\*SubscriptBox[\(>\), \(SU(2), \ 1 - sites\)]\) <(e^(Subscript[k, 2] z)),\!\(Q*)
(*\*SubscriptBox[\(>\), \(\(SU(3)\)\ 1 - site\)]\)*)


(* ::Subsubsection::Closed:: *)
(*{1,1,0,0}*)


(* ::Input:: *)
(*rep1T={1,1,0,0};*)
(*Weight=FindGoodWeight[rep1T];*)
(*MaxW=Length[FindGoodWeight[rep1T]];*)
(*WeightAll=FindWeight [rep1T];*)
(*MaxWAll=Length[FindWeight[rep1T]];*)
(**)
(*(**)
(*\[Phi]1T=1;*)
(*\[Phi]2T=1/2;*)
(*\[Phi]3T=1/3;*)
(*\[Phi]4T=1/5;*)
(**)*)
(**)
(**)
(**)
(*\[Phi]1T=\[Phi]1;*)
(*\[Phi]2T=\[Phi]2;*)
(*\[Phi]3T=\[Phi]3;*)
(*\[Phi]4T=\[Phi]4;*)
(**)


(* ::Input:: *)
(*ClearAll[EQS];*)
(*EQS=Flatten[Table[TESTK1Site22B[rep1T,WeightAll,MaxWAll,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4][i],{i,1,MaxWAll}]]//Simplify*)


(* ::Input:: *)
(*Solve[EQS==0,{p11,p22,p33,p44,p12,p13,p14,p23,p24,p34}]//Simplify*)


(* ::Subsubsection::Closed:: *)
(*{2,2,0,0}*)


(* ::Input:: *)
(*rep1T={2,2,0,0};*)
(*Weight=FindGoodWeight[rep1T];*)
(*MaxW=Length[FindGoodWeight[rep1T]];*)
(*WeightAll=FindWeight [rep1T];*)
(*MaxWAll=Length[FindWeight[rep1T]];*)
(**)
(**)
(*\[Phi]1T=1;*)
(*\[Phi]2T=1/2;*)
(*\[Phi]3T=1/3;*)
(*\[Phi]4T=1/5;*)
(**)
(**)
(**)
(*(**)
(*\[Phi]1T=\[Phi]1;*)
(*\[Phi]2T=\[Phi]2;*)
(*\[Phi]3T=\[Phi]3;*)
(*\[Phi]4T=\[Phi]4;*)
(**)*)


(* ::Input:: *)
(*ClearAll[EQS];*)
(*EQS=Flatten[Table[TESTK1Site22B[rep1T,WeightAll,Max,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4][i],{i,1,MaxWAll}]]//Simplify*)


(* ::Input:: *)
(*EQS/.{p22->1/p11,p13->(\[Alpha]\[Phi]13 (k3-\[Phi]3T) (-k4+\[Phi]2T+\[Phi]3T))/((k3-\[Phi]1T) (-k4+\[Phi]1T+\[Phi]2T)),p14->(\[Alpha]\[Phi]14 (k3-\[Phi]4T) (-k4+\[Phi]2T+\[Phi]4T))/((k3-\[Phi]1T) (-k4+\[Phi]1T+\[Phi]2T)),p23->(\[Alpha]\[Phi]23 (k3-\[Phi]3T) (-k4+\[Phi]1T+\[Phi]3T))/((k3-\[Phi]2T) (-k4+\[Phi]1T+\[Phi]2T)),p24->-((\[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (\[Phi]2T-\[Phi]3T) (k3-\[Phi]4T) (-\[Phi]1T+\[Phi]4T) (-k4+\[Phi]1T+\[Phi]4T))/(\[Alpha]\[Phi]13 (k3-\[Phi]2T) (-k4+\[Phi]1T+\[Phi]2T) (\[Phi]1T-\[Phi]3T) (\[Phi]2T-\[Phi]4T)))}//Simplify*)


(* ::Subsubsection::Closed:: *)
(*{3,3,0,0}*)


(* ::Input:: *)
(*rep1T={3,3,0,0};*)
(*Weight=FindGoodWeight[rep1T];*)
(*MaxW=Length[FindGoodWeight[rep1T]];*)
(*WeightAll=FindWeight [rep1T];*)
(*MaxWAll=Length[FindWeight[rep1T]];*)
(**)
(**)
(*\[Phi]1T=1;*)
(*\[Phi]2T=1/2;*)
(*\[Phi]3T=1/3;*)
(*\[Phi]4T=1/5;*)
(**)
(**)
(**)
(*(**)
(*\[Phi]1T=\[Phi]1;*)
(*\[Phi]2T=\[Phi]2;*)
(*\[Phi]3T=\[Phi]3;*)
(*\[Phi]4T=\[Phi]4;*)
(**)*)


(* ::Input:: *)
(*ClearAll[EQS];*)
(*EQS=Flatten[Table[TESTK1Site22B[rep1T,WeightAll,Max,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4][i],{i,1,MaxWAll}]]//Simplify*)


(* ::Input:: *)
(*EQS/.{p22->1/p11,p13->(\[Alpha]\[Phi]13 (k3-\[Phi]3T) (-k4+\[Phi]2T+\[Phi]3T))/((k3-\[Phi]1T) (-k4+\[Phi]1T+\[Phi]2T)),p14->(\[Alpha]\[Phi]14 (k3-\[Phi]4T) (-k4+\[Phi]2T+\[Phi]4T))/((k3-\[Phi]1T) (-k4+\[Phi]1T+\[Phi]2T)),p23->(\[Alpha]\[Phi]23 (k3-\[Phi]3T) (-k4+\[Phi]1T+\[Phi]3T))/((k3-\[Phi]2T) (-k4+\[Phi]1T+\[Phi]2T)),p24->-((\[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (\[Phi]2T-\[Phi]3T) (k3-\[Phi]4T) (-\[Phi]1T+\[Phi]4T) (-k4+\[Phi]1T+\[Phi]4T))/(\[Alpha]\[Phi]13 (k3-\[Phi]2T) (-k4+\[Phi]1T+\[Phi]2T) (\[Phi]1T-\[Phi]3T) (\[Phi]2T-\[Phi]4T)))}//Simplify*)


(* ::Title::Closed:: *)
(*TEST matching with Article Paul*)


(* ::Subsubsection:: *)
(*SU(2)*)


(* ::Input:: *)
(*SU2OrthogRen[g[z],f[z],1/2,1/2]*)


(* ::Input:: *)
(*SU2OrthogRen1[g[z],f[z],1/2]*)


(* ::Subsubsection:: *)
(*SU(3)*)


(* ::Input:: *)
(*SU3OrthogRen1[{QHL[1][z],QHL[2][z],QHL[3][z]},{QL[1][z],QL[2][z],QL[3][z]},{QHR[1][z],QHR[2][z],QHR[3][z]},{QR[1][z],QR[2][z],QR[3][z]},{1,1,0}]*)


(* ::Subsubsection:: *)
(*SU(3) scalar product matrix*)


(* ::Input:: *)
(*ClearAll[SU3OrthogRen1M];*)
(*SU3OrthogRen1M[QHL_,QL_,QHR_,QR_,rep1_]:=SU3OrthogRen1M[QHL,QL,QHR,QR,rep1]={{SU3MeanValue1[QHL[[1]]*D[QR[[1]],z]*1/z,0,rep1],SU3MeanValue1[QHL[[1]]*QR[[1]]*1/z^2,0,rep1],SU3MeanValue1[QHL[[1]]*QR[[1]]*1/z,0,rep1]},*)
(**)
(*{SU3MeanValue1[QHL[[1]]*D[QR[[2]],z]*1/z,0,rep1],SU3MeanValue1[QHL[[1]]*QR[[2]]*1/z^2,0,rep1],SU3MeanValue1[QHL[[1]]*QR[[2]]*1/z,0,rep1]},*)
(**)
(*{SU3MeanValue1[QHL[[2]]*D[QR[[1]],z]*1/z,0,rep1],SU3MeanValue1[QHL[[2]]*QR[[1]]*1/z^2,0,rep1],SU3MeanValue1[QHL[[2]]*QR[[1]]*1/z,0,rep1]}};*)


(* ::Input:: *)
(*SU3OrthogRen1M[{QHL[1][z],QHL[2][z],QHL[3][z]},{QL[1][z],QL[2][z],QL[3][z]},{QHR[1][z],QHR[2][z],QHR[3][z]},{QR[1][z],QR[2][z],QR[3][z]},{1,1,0}]*)
