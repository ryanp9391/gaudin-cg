(* Converted from SU(4)_Paul_Normalisation_EXPLICIT.nb; notebook presentation/output cells removed. *)


(* === RUN: Introduction === *)

(* \<\ The gl(3) representation are labeled by the eigenvalues of the Cartans of the \ highest weight state: Rep1= {\[Lambda][1,1],\[Lambda][1,2],\[Lambda][1,3]} Rep2={\[Lambda][2,1],\[Lambda][2,2],\[Lambda][2,3]}  We have that: \[Lambda][*,1]\[GreaterEqual]\[Lambda][*,2]\[GreaterEqual]\[Lambda][*,3]\ \> *)

ClearAll[w1,w2];
w1=0;
w2=1;

ClearAll[FUNDAMENTAL, ANTIFUNDAMENTAL, ADJOINT];
(*Rep of gl(4)*)
FUNDAMENTAL={1,0,0,0};
ANTIFUNDAMENTAL={1,1,1,0};
ADJOINT={2,1,1,0};

Ng=4;


(* === RUN: GL(4) generators === *)


(* === Gl(4) generators === *)

ps[rep_]:=ps[rep]=Block[{res=(Prepend[rep,#]&/@(Tuples[Table[Range[rep[[1,i+1]],rep[[1,i]]],{i,1,Length[rep[[1]]]-1}]]))},If[Length[res[[1,1]]]==1,Reverse@Sort@res,Reverse@Sort@Flatten[ps/@res,1]]];
dim[rep_]:=dim[rep]=Length[ps[{rep}]];
J0[k_][pat_]:=J0[k][pat]=Sum[pat[[k,i]],{i,1,k}]-Sum[pat[[k-1,i]],{i,1,k-1}];
Jp[k_][pat_]:=Jp[k][pat]=-Sum[Product[pat[[k,i]]-i-pat[[k+1,j]]+j,{j,Range[k+1]}]/Product[pat[[k,i]]-i-pat[[k,j]]+j,{j,Complement[Range[k],{i}]}] v[ReplacePart[pat,{k,i}->pat[[k,i]]+1]],{i,1,k}];
Jm[k_][pat_]:=Jm[k][pat]=Sum[Product[pat[[k,i]]-i-pat[[k-1,j]]+j,{j,Range[k-1]}]/Product[pat[[k,i]]-i-pat[[k,j]]+j,{j,Complement[Range[k],{i}]}] v[ReplacePart[pat,{k,i}->pat[[k,i]]-1]],{i,1,k}];
JJ[k_,l_,rep_]:=JJ[k,l,rep]=Which[k==l,DiagonalMatrix[J0[k]/@ps[{rep}]],l==k+1,Transpose[Table[Coefficient[Jp[k][ps[{rep}][[i]]],v[#]]&/@ps[{rep}],{i,1,Length[ps[{rep}]]}]],l==k-1,Transpose[Table[Coefficient[Jm[l][ps[{rep}][[i]]],v[#]]&/@ps[{rep}],{i,1,Length[ps[{rep}]]}]],l>k+1,JJ[k,k+1,rep].JJ[k+1,l,rep]-JJ[k+1,l,rep].JJ[k,k+1,rep],l<k-1,JJ[k,k-1,rep].JJ[k-1,l,rep]-JJ[k-1,l,rep].JJ[k,k-1,rep]];
CC[k_,set_,rep_]:=CC[k,set,rep]=Sum[Dot[Sequence@@Table[JJ[j[i+1],j[i],rep],{i,1,k}]]/.j[k+1]->j[1],Evaluate[Sequence@@Array[{j[#],set}&,k]]];
CCDiag[k_,set_,rep_]:=CCDiag[k,set,rep]=Diagonal[Inverse[Transpose[Eigenvectors[CC[Length[set],set,rep]]]].CC[k,set,rep].Transpose[Eigenvectors[CC[Length[set],set,rep]]]]


(* === Hermitian generators === *)

ClearAll[JJH,PM,SolvePH,PH,PElement];
PM[rep_]:=PM[rep]=DiagonalMatrix[Table[PElement[i],{i,1,dim[rep]}]];
SolvePH[rep_]:=SolvePH[rep]=FullSimplify[(Solve[(DeleteDuplicates[Table[PM[rep].JJ[k,l,rep].Inverse[PM[rep]]-Inverse[PM[rep]].ConjugateTranspose[JJ[l,k,rep]].PM[rep],{k,1,Ng},{l,1,Ng}]//Flatten])==0,Table[PElement[i],{i,1,dim[rep]}]])];(*If PM is real*)
PH[rep_,Optional[M_,1]]:=PH[rep,Optional[M,1]]=(PM[rep]/.SolvePH[rep])[[M]]//Simplify;
(*If PM is complex*)(*Solve[(DeleteDuplicates[Table[PM[rep].JJ[k,l,rep].Inverse[PM[rep]]-ConjugateTranspose[PM[rep].JJ[l,k,rep].Inverse[PM[rep]]],{k,1,3},{l,1,3}]//Flatten])==0,Table[PElement[i],{i,1,dim[rep]}]];*)
JJH[k_,l_,rep_,Optional[M_,1]]:=JJH[k,l,rep,Optional[M,1]]=(PH[rep,M].JJ[k,l,rep].Inverse[PH[rep,M]])//Simplify;


(* === Definition matrices for Lax connection === *)

ClearAll[e,EM,Id];
e[a_,b_,rep1_]:=e[a,b,rep1]=SparseArray[{a,b}->1,{dim[rep1],dim[rep1]}];
EM[a_,b_,rep1_]:=EM[a,b,rep1]=JJ[a,b,rep1];(*Non Hermitian generators*)
(*JJH[a,b,rep1];*)
Id[rep1_]:=Id[rep1]=SparseArray[{a_,a_}->1,{dim[rep1],dim[rep1]}];

(* TextData[{Generator of overall rotation , Cell[(((E)^(All)), ExpressionUUID->   "97ad50d3-58db-461b-9869-c6ebd068c312"], )}] *)

ClearAll[EMAll];
EMAll[a_,b_,rep1_,rep2_]:=EMAll[a,b,rep1,rep2]=KroneckerProduct[EM[a,b,rep1],Id[rep2]]+KroneckerProduct[Id[rep1],EM[a,b,rep2]];


(* === Module to find good weights for the rep === *)

ClearAll[FindWeight,FindGoodWeight];
FindWeight[rep_]:=Table[{EM[1,1,rep][[n,n]],EM[2,2,rep][[n,n]],EM[3,3,rep][[n,n]],EM[4,4,rep][[n,n]]},{n,1,dim[rep]}];
FindGoodWeight[rep_]:=Pick[FindWeight[rep],Lookup[Counts[FindWeight[rep]],FindWeight[rep]],1]

ClearAll[FindWeight2All,FindWeight2,FindGoodWeight2]
FindWeight2All[rep1_,rep2_]:=FindWeight2All[rep1,rep2]=Flatten[Table[{EM[1,1,rep1][[n,n]]+EM[1,1,rep2][[m,m]],
EM[2,2,rep1][[n,n]]+EM[2,2,rep2][[m,m]],
EM[3,3,rep1][[n,n]]+EM[3,3,rep2][[m,m]],EM[4,4,rep1][[n,n]]+EM[4,4,rep2][[m,m]]},{n,1,dim[rep1]},{m,1,dim[rep2]}],1];
FindWeight2[rep1_,rep2_]:=FindWeight2[rep1,rep2]=DeleteDuplicates[FindWeight2All[rep1,rep2]];
FindGoodWeight2[rep1_,rep2_]:=FindGoodWeight2[rep1,rep2]=Pick[FindWeight2All[rep1,rep2],Lookup[Counts[FindWeight2All[rep1,rep2]],FindWeight2All[rep1,rep2]],1]

ClearAll[WeightDynkin,\[Lambda]1,\[Lambda]2,\[Lambda]3];
WeightDynkin[\[CapitalDelta]_,l_]:=Solve[\[Lambda]1-\[Lambda]2==(\[CapitalDelta])&&\[Lambda]2-\[Lambda]3==l&&\[Lambda]1+\[Lambda]2+\[Lambda]3==0,{\[Lambda]1,\[Lambda]2,\[Lambda]3}]
(*Fundamental*)
WeightDynkin[1,0]
(*AntiFundamental*)
WeightDynkin[0,1]
(*Adjoint*)
WeightDynkin[1,1]


(* === RUN: Q space === *)


(* === hs === *)


(* === Cartan === *)

ClearAll[h];
h[rep1_,rep2_][1,1]:=h[rep1,rep2][1,1]=rep1[[1]]+rep1[[2]]+rep1[[3]]+rep1[[4]];
h[rep1_,rep2_][2,1]:=h[rep1,rep2][2,1]=rep1\[LeftDoubleBracket]3\[RightDoubleBracket] (-1+rep1\[LeftDoubleBracket]4\[RightDoubleBracket])+rep1\[LeftDoubleBracket]2\[RightDoubleBracket] (-2+rep1\[LeftDoubleBracket]3\[RightDoubleBracket]+rep1\[LeftDoubleBracket]4\[RightDoubleBracket])+rep1\[LeftDoubleBracket]1\[RightDoubleBracket] (-3+rep1\[LeftDoubleBracket]2\[RightDoubleBracket]+rep1\[LeftDoubleBracket]3\[RightDoubleBracket]+rep1\[LeftDoubleBracket]4\[RightDoubleBracket]);
h[rep1_,rep2_][3,1]:=h[rep1,rep2][3,1]=rep1\[LeftDoubleBracket]2\[RightDoubleBracket] (-1+rep1\[LeftDoubleBracket]3\[RightDoubleBracket]) (-2+rep1\[LeftDoubleBracket]4\[RightDoubleBracket])+rep1\[LeftDoubleBracket]1\[RightDoubleBracket] ((-2+rep1\[LeftDoubleBracket]3\[RightDoubleBracket]) (-3+rep1\[LeftDoubleBracket]4\[RightDoubleBracket])+rep1\[LeftDoubleBracket]2\[RightDoubleBracket] (-4+rep1\[LeftDoubleBracket]3\[RightDoubleBracket]+rep1\[LeftDoubleBracket]4\[RightDoubleBracket]));
h[rep1_,rep2_][4,1]:=h[rep1,rep2][4,1]=rep1\[LeftDoubleBracket]1\[RightDoubleBracket] (-1+rep1\[LeftDoubleBracket]2\[RightDoubleBracket]) (-2+rep1\[LeftDoubleBracket]3\[RightDoubleBracket]) (-3+rep1\[LeftDoubleBracket]4\[RightDoubleBracket]);

h[rep1_,rep2_][1,2]:=h[rep1,rep2][1,2]=rep2[[1]]+rep2[[2]]+rep2[[3]]+rep2[[4]];
h[rep1_,rep2_][2,2]:=h[rep1,rep2][2,2]=rep2\[LeftDoubleBracket]3\[RightDoubleBracket] (-1+rep2\[LeftDoubleBracket]4\[RightDoubleBracket])+rep2\[LeftDoubleBracket]2\[RightDoubleBracket] (-2+rep2\[LeftDoubleBracket]3\[RightDoubleBracket]+rep2\[LeftDoubleBracket]4\[RightDoubleBracket])+rep2\[LeftDoubleBracket]1\[RightDoubleBracket] (-3+rep2\[LeftDoubleBracket]2\[RightDoubleBracket]+rep2\[LeftDoubleBracket]3\[RightDoubleBracket]+rep2\[LeftDoubleBracket]4\[RightDoubleBracket]);
h[rep1_,rep2_][3,2]:=h[rep1,rep2][3,2]=rep2\[LeftDoubleBracket]2\[RightDoubleBracket] (-1+rep2\[LeftDoubleBracket]3\[RightDoubleBracket]) (-2+rep2\[LeftDoubleBracket]4\[RightDoubleBracket])+rep2\[LeftDoubleBracket]1\[RightDoubleBracket] ((-2+rep2\[LeftDoubleBracket]3\[RightDoubleBracket]) (-3+rep2\[LeftDoubleBracket]4\[RightDoubleBracket])+rep2\[LeftDoubleBracket]2\[RightDoubleBracket] (-4+rep2\[LeftDoubleBracket]3\[RightDoubleBracket]+rep2\[LeftDoubleBracket]4\[RightDoubleBracket]));
h[rep1_,rep2_][4,2]:=h[rep1,rep2][4,2]=rep2\[LeftDoubleBracket]1\[RightDoubleBracket] (-1+rep2\[LeftDoubleBracket]2\[RightDoubleBracket]) (-2+rep2\[LeftDoubleBracket]3\[RightDoubleBracket]) (-3+rep2\[LeftDoubleBracket]4\[RightDoubleBracket]);


(* === TextData[{Cell[((z)^(0)), ExpressionUUID->   "b172c059-a326-4b3d-a1b1-9de1d980270e"],  (dependence only on regulators)}] === *)

(*H1*)
h[\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][1,3]:=h[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][1,3]=\[Phi]1+\[Phi]2+\[Phi]3+\[Phi]4;

(*H2*)
h[\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][2,5]:=h[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][2,5]=\[Phi]1 \[Phi]2+\[Phi]1 \[Phi]3+\[Phi]2 \[Phi]3+\[Phi]1 \[Phi]4+\[Phi]2 \[Phi]4+\[Phi]3 \[Phi]4;

(*H3*)
h[\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][3,7]:=h[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][3,7]=\[Phi]1 \[Phi]2 \[Phi]3+\[Phi]1 \[Phi]2 \[Phi]4+\[Phi]1 \[Phi]3 \[Phi]4+\[Phi]2 \[Phi]3 \[Phi]4;

(*H4*)
h[\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][4,9]:=h[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][4,9]=\[Phi]1*\[Phi]2*\[Phi]3*\[Phi]4;


(* === Relation (dependence on repAll) === *)

(* TextData[{By expanding at z=\[Infinity] and taking the coefficient of , Cell[((1)/(z)), ExpressionUUID->   "197dd7a4-3488-44a5-be31-62f5b45469f9"]}] *)

(*H2*)
h[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][2,3]:=h[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][2,3]=-h[2,4]+(\[Phi]2+\[Phi]3+\[Phi]4) (repAll[[1]])+(\[Phi]1+\[Phi]3+\[Phi]4) (repAll[[2]])+(\[Phi]1+\[Phi]2+\[Phi]4) (repAll[[3]])+(\[Phi]1+\[Phi]2+\[Phi]3) (repAll[[4]]);

(*H3*)
h[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][3,5]:=h[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][3,5]=-h[3,6]+\[Phi]3 \[Phi]4 (repAll[[1]]+repAll[[2]])+\[Phi]2 \[Phi]4 (repAll[[1]]+repAll[[3]])+\[Phi]1 \[Phi]4 (repAll[[2]]+repAll[[3]])+\[Phi]2 \[Phi]3 (repAll[[1]]+repAll[[4]])+\[Phi]1 \[Phi]3 (repAll[[2]]+repAll[[4]])+\[Phi]1 \[Phi]2 (repAll[[3]]+repAll[[4]]);

(*H4*)
h[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][4,7]:=h[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][4,7]=-h[4,8]+\[Phi]2 \[Phi]3 \[Phi]4 (repAll[[1]])+\[Phi]1 \[Phi]3 \[Phi]4 (repAll[[2]])+\[Phi]1 \[Phi]2 \[Phi]4 (repAll[[3]])+\[Phi]1 \[Phi]2 \[Phi]3 (repAll[[4]]);


(* === Hamiltonian in Q-space === *)

ClearAll[\[Epsilon]1,\[Epsilon]2,\[Epsilon]3];
\[Epsilon]1[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]1[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((h[rep1,rep2][1,1])/(z-w1))+((h[rep1,rep2][1,2])/(z-w2))+h[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][1,3];
\[Epsilon]2[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]2[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((h[rep1,rep2][2,1])/((((z-w1))^(2))))+((h[rep1,rep2][2,2])/((((z-w2))^(2))))+((h[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][2,3])/((z-w1)))+((h[2,4])/((z-w2)))+h[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][2,5];
\[Epsilon]3[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]3[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((h[rep1,rep2][3,1])/((((z-w1))^(3))))+((h[rep1,rep2][3,2])/((((z-w2))^(3))))+((h[3,3])/((((z-w1))^(2))))+((h[3,4])/((((z-w2))^(2))))+((h[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][3,5])/(z-w1))+((h[3,6])/(z-w2))+h[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][3,7];
\[Epsilon]4[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]4[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((h[rep1,rep2][4,1])/((((z-w1))^(4))))+((h[rep1,rep2][4,2])/((((z-w2))^(4))))+((h[4,3])/((((z-w1))^(3))))+((h[4,4])/((((z-w2))^(3))))+((h[4,5])/((((z-w1))^(2))))+((h[4,6])/((((z-w2))^(2))))+((h[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][4,7])/(z-w1))+((h[4,8])/(z-w2))+h[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][4,9];


(* === Baxter equation with F === *)

ClearAll[F];
F[rep1_,rep2_,C1_,
C2_][z_]:=F[rep1,rep2,C1,C2][z]=Cases[{{C1,C2}},{a_,b_}:>((1)/((((z-w1))^(rep1[[a]]-a+1))*(((z-w2))^(rep2[[b]]-b+1))))][[1]];

(* Definition of shifted Hamiltonian *)

ClearAll[\[Epsilon]1S,\[Epsilon]2S,\[Epsilon]3S]
\[Epsilon]1S[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_,C2_][z_]:=\[Epsilon]1S[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2][z]=\[Epsilon]1[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+((4 ((F[rep1,rep2,C1,C2])^(\[Prime]))[z])/(F[rep1,rep2,C1,C2][z]));
\[Epsilon]2S[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_,C2_][z_]:=\[Epsilon]2S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2][z]=\[Epsilon]2[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+((3 \[Epsilon]1[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] ((F[rep1,rep2,C1,C2])^(\[Prime]))[z])/(F[rep1,rep2,C1,C2][z]))+((6 ((F[rep1,rep2,C1,C2])^(\[Prime]\[Prime]))[z])/(F[rep1,rep2,C1,C2][z]));
\[Epsilon]3S[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_,C2_][z_]:=\[Epsilon]3S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2][z]=\[Epsilon]3[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+2*\[Epsilon]2[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*((((F[rep1,rep2,C1,C2])^(\[Prime]))[z])/(F[rep1,rep2,C1,C2][z]))+3*((\[Epsilon]1[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] ((F[rep1,rep2,C1,C2])^(\[Prime]\[Prime]))[z])/(F[rep1,rep2,C1,C2][z]))+4*((((F[rep1,rep2,C1,C2])^((3)))[z])/(F[rep1,rep2,C1,C2][z]));
\[Epsilon]4S[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_,C2_][z_]:=\[Epsilon]4S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2][z]=\[Epsilon]4[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+\[Epsilon]3[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*((((F[rep1,rep2,C1,C2])^(\[Prime]))[z])/(F[rep1,rep2,C1,C2][z]))+((\[Epsilon]2[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] ((F[rep1,rep2,C1,C2])^(\[Prime]\[Prime]))[z])/(F[rep1,rep2,C1,C2][z]))+\[Epsilon]1[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*((((F[rep1,rep2,C1,C2])^((3)))[z])/(F[rep1,rep2,C1,C2][z]))+((((F[rep1,rep2,C1,C2])^((4)))[z])/(F[rep1,rep2,C1,C2][z]));

ClearAll[BaxOp];
BaxOp[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_,C2_]:=D[#,{z,4}]+\[Epsilon]1S[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2][z]*D[#,{z,3}]+\[Epsilon]2S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2][z]*D[#,{z,2}]+\[Epsilon]3S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2][z]*D[#,z]+\[Epsilon]4S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2][z]*#&;


(* === Q functions === *)

ClearAll[QPoly,Q,c1,c2,c3,c4];
QPoly[c_,M_,r_][z_]:=QPoly[c,M,r][z]=Exp[r z]*(((z)^(M))+Sum[c[k]*((z)^(k)),{k,0,M-1}]);
Q[1][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q[1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly[c1,rep1\[LeftDoubleBracket]1\[RightDoubleBracket]+rep2\[LeftDoubleBracket]1\[RightDoubleBracket]-repAll\[LeftDoubleBracket]1\[RightDoubleBracket],-\[Phi]1][z];
Q[2][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q[2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly[c2,rep1\[LeftDoubleBracket]1\[RightDoubleBracket]+rep2\[LeftDoubleBracket]1\[RightDoubleBracket]-repAll\[LeftDoubleBracket]2\[RightDoubleBracket],-\[Phi]2][z];
Q[3][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q[3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly[c3,rep1\[LeftDoubleBracket]1\[RightDoubleBracket]+rep2\[LeftDoubleBracket]1\[RightDoubleBracket]-repAll\[LeftDoubleBracket]3\[RightDoubleBracket],-\[Phi]3][z];
Q[4][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q[4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly[c4,rep1\[LeftDoubleBracket]1\[RightDoubleBracket]+rep2\[LeftDoubleBracket]1\[RightDoubleBracket]-repAll\[LeftDoubleBracket]4\[RightDoubleBracket],-\[Phi]4][z];


(* === Baxter equations === *)

(* Baxter equations *)

ClearAll[Bax];
Bax[1][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_,C2_]:=Bax[1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2]=((BaxOp[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2][Q[1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[\[Phi]1 z])//Simplify;
Bax[2][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_,C2_]:=Bax[2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2]=((BaxOp[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2][Q[2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[\[Phi]2 z])//Simplify;
Bax[3][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_,C2_]:=Bax[3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2]=((BaxOp[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2][Q[3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[\[Phi]3 z])//Simplify;
Bax[4][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_,C2_]:=Bax[4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2]=((BaxOp[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1,C2][Q[4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[\[Phi]4 z])//Simplify;

(* Unknown coefficients *)

ClearAll[Vars1,Vars2,Vars3,Vars4];
Vars1[rep1_,rep2_,repAll_]:=Vars1[rep1,rep2,repAll]=Table[c1[n],{n,0,rep1\[LeftDoubleBracket]1\[RightDoubleBracket]+rep2\[LeftDoubleBracket]1\[RightDoubleBracket]-repAll\[LeftDoubleBracket]1\[RightDoubleBracket]-1}];
Vars2[rep1_,rep2_,repAll_]:=Vars2[rep1,rep2,repAll]=Table[c2[n],{n,0,rep1\[LeftDoubleBracket]1\[RightDoubleBracket]+rep2\[LeftDoubleBracket]1\[RightDoubleBracket]-repAll\[LeftDoubleBracket]2\[RightDoubleBracket]-1}];
Vars3[rep1_,rep2_,repAll_]:=Vars3[rep1,rep2,repAll]=Table[c3[n],{n,0,rep1\[LeftDoubleBracket]1\[RightDoubleBracket]+rep2\[LeftDoubleBracket]1\[RightDoubleBracket]-repAll\[LeftDoubleBracket]3\[RightDoubleBracket]-1}];
Vars4[rep1_,rep2_,repAll_]:=Vars4[rep1,rep2,repAll]=Table[c4[n],{n,0,rep1\[LeftDoubleBracket]1\[RightDoubleBracket]+rep2\[LeftDoubleBracket]1\[RightDoubleBracket]-repAll\[LeftDoubleBracket]4\[RightDoubleBracket]-1}];

(* Solve Baxter equation *)

ClearAll[SolveBaxAll];
SolveBaxAll[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,c1_,c2_]:=SolveBaxAll[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,c1,c2]=Solve[CoefficientList[Bax[1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,c1,c2]//Together//Numerator,z]==0&&
CoefficientList[Bax[2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,c1,c2]//Together//Numerator,z]==0&&
CoefficientList[Bax[3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,c1,c2]//Together//Numerator,z]==0&&CoefficientList[Bax[4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,c1,c2]//Together//Numerator,z]==0,Join[Vars1[rep1,rep2,repAll],Vars2[rep1,rep2,repAll],Vars3[rep1,rep2,repAll],Vars4[rep1,rep2,repAll],{h[2,4],h[3,3],h[3,4],h[3,6],h[4,3],h[4,4],h[4,5],h[4,6],h[4,8]}]]//Simplify;

(* Spectrum of Hamiltonian *)

ClearAll[HSpectrum];
HSpectrum[1][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=HSpectrum[1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]1[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]/.SolveBaxAll[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1];

HSpectrum[2][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=HSpectrum[2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]2[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]/.SolveBaxAll[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1];

HSpectrum[3][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=HSpectrum[3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]3[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]/.SolveBaxAll[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1];

HSpectrum[4][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=HSpectrum[4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]4[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]/.SolveBaxAll[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1];


(* === Example of use === *)

Sol=SolveBaxAll[FUNDAMENTAL,FUNDAMENTAL,{2,0,0,0},1,((1)/(3)),((1)/(5)),((1)/(7)),1,1]

H1ST=HSpectrum[1][FUNDAMENTAL,FUNDAMENTAL,{2,0,0,0},1,((1)/(3)),((1)/(5)),((1)/(7))][z];
H2ST=HSpectrum[2][FUNDAMENTAL,FUNDAMENTAL,{2,0,0,0},1,((1)/(3)),((1)/(5)),((1)/(7))][z];
H3ST=HSpectrum[3][FUNDAMENTAL,FUNDAMENTAL,{2,0,0,0},1,((1)/(3)),((1)/(5)),((1)/(7))][z];
H4ST=HSpectrum[4][FUNDAMENTAL,FUNDAMENTAL,{2,0,0,0},1,((1)/(3)),((1)/(5)),((1)/(7))][z];

SolveBaxAll[ADJOINT,ADJOINT,ADJOINT+ADJOINT,1,((1)/(3)),((1)/(5)),((1)/(7)),1,1]


(* === Measure === *)


(* === 2-sites === *)

(* TextData[{Measure such that <f , Cell[Subscript[O, F]g> Overscript[=, !], ExpressionUUID->   "2832d562-1c3b-4ec8-bb01-a613c658cb0e"],  <g , Cell[((Subscript[O, F])^(\[Dagger])) f>, ExpressionUUID->"0015f24a-3c6d-4dca-aa2a-020cecfd5118"]}] *)

ClearAll[\[Mu]];
\[Mu][rep1_,rep2_][z_]:=\[Mu][rep1,rep2][z]=(((-1+z))^(-rep2[[1]]+rep2[[4]])) ((z)^(-rep1[[1]]+rep1[[4]]));

(* Mean value *)

ClearAll[MeanValue];
MeanValue[f_,\[Tau]_,rep1_,rep2_]:=MeanValue[f,\[Tau],rep1,rep2]=
2 \[Pi] I Residue[(f*\[Mu][rep1,rep2][z]),{z,\[Tau]}];


(* === 1-site === *)

(* TextData[{Measure such that <f , Cell[Subscript[O, F]g> Overscript[=, !], ExpressionUUID->   "4a94f0b5-fdcb-45ff-9339-a06cb39f04e4"],  <g , Cell[((Subscript[O, F])^(\[Dagger])) f>, ExpressionUUID->"51ba7f18-7e8e-4cd8-a7d2-c8804fb1f9f5"]}] *)

ClearAll[\[Mu]1,\[Mu]2];
\[Mu]1[rep1_][z_]:=\[Mu]1[rep1][z]=((z)^(-rep1[[1]]+rep1[[4]]));
\[Mu]2[rep2_][z_]:=\[Mu]2[rep2][z]=(((-1+z))^(-rep2[[1]]+rep2[[4]]));

(* Mean value *)

ClearAll[MeanValue1,MeanValue2];
MeanValue1[f_,\[Tau]_,rep1_]:=MeanValue1[f,\[Tau],rep1]=
2 \[Pi] I Residue[(f*\[Mu]1[rep1][z]),{z,\[Tau]}];
MeanValue2[f_,\[Tau]_,rep2_]:=MeanValue2[f,\[Tau],rep2]=
2 \[Pi] I Residue[(f*\[Mu]2[rep2][z]),{z,\[Tau]}];


(* === Test === *)

MeanValue[f[z]*((1)/(z*(z-1))),0,{1,1,0,0},{1,1,0,0}]
MeanValue[f[z]*((1)/(z*(z-1))),1,{1,1,0,0},{1,1,0,0}]
MeanValue1[f[z]*((1)/(z)),0,{1,1,0,0}]


(* === TextData[{Adjoint Baxter equation with , Cell[((F)^(\[Dagger])), ExpressionUUID->   "f17c81db-e1df-4939-ab15-43e078aad370"]}] === *)

(* \<\ HERE we consider only the case C1=1,C2=1 in the definition of M\ \> *)

(* TextData[{Shifted Hamiltonians (already with , Cell[((F)^(\[Dagger])), ExpressionUUID->   "8fb23737-07e7-487e-98f5-aa203054a6ae"], )}] *)

ClearAll[\[Epsilon]1S\[Dagger],\[Epsilon]2S\[Dagger],\[Epsilon]3S\[Dagger]];
\[Epsilon]1S\[Dagger][rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]1S\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=-\[Epsilon]1S[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z]+4*((D[\[Mu][rep1,rep2][z],z])/(\[Mu][rep1,rep2][z]));
\[Epsilon]2S\[Dagger][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]2S\[Dagger][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]2S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z]-3*D[\[Epsilon]1S[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z],z]-3*\[Epsilon]1S[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z]*((D[\[Mu][rep1,rep2][z],z])/(\[Mu][rep1,rep2][z]))+6*((D[\[Mu][rep1,rep2][z],{z,2}])/(\[Mu][rep1,rep2][z]));
\[Epsilon]3S\[Dagger][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]3S\[Dagger][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=-\[Epsilon]3S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z]+2*D[\[Epsilon]2S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z],z]-3*D[\[Epsilon]1S[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z],{z,2}]+((D[\[Mu][rep1,rep2][z],z])/(\[Mu][rep1,rep2][z]))*(2*\[Epsilon]2S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z]-6*D[\[Epsilon]1S[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z],z])-3*((D[\[Mu][rep1,rep2][z],{z,2}])/(\[Mu][rep1,rep2][z]))*\[Epsilon]1S[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z]+4*((D[\[Mu][rep1,rep2][z],{z,3}])/(\[Mu][rep1,rep2][z]));
\[Epsilon]4S\[Dagger][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]4S\[Dagger][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]4S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z]-D[\[Epsilon]3S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z],z]+D[\[Epsilon]2S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z],{z,2}]-D[\[Epsilon]1S[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z],{z,3}]+((D[\[Mu][rep1,rep2][z],z])/(\[Mu][rep1,rep2][z]))*(-\[Epsilon]3S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z]+2*D[\[Epsilon]2S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z],z]-3*D[\[Epsilon]1S[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z],{z,2}])+((D[\[Mu][rep1,rep2][z],{z,2}])/(\[Mu][rep1,rep2][z]))*(\[Epsilon]2S[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z]-3D[\[Epsilon]1S[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z],z])-\[Epsilon]1S[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1][z]*((D[\[Mu][rep1,rep2][z],{z,3}])/(\[Mu][rep1,rep2][z]))+((D[\[Mu][rep1,rep2][z],{z,4}])/(\[Mu][rep1,rep2][z]));

ClearAll[BaxOp\[Dagger]];
BaxOp\[Dagger][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=D[#,{z,4}]+\[Epsilon]1S\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*D[#,{z,3}]+\[Epsilon]2S\[Dagger][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*D[#,{z,2}]+\[Epsilon]3S\[Dagger][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*D[#,z]+\[Epsilon]4S\[Dagger][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*#&;


(* === Q functions === *)

ClearAll[QPoly\[Dagger],Q\[Dagger],c1\[Dagger],c2\[Dagger],c3\[Dagger],c4\[Dagger]];
QPoly\[Dagger][c_,M_,r_][z_]:=QPoly\[Dagger][c,M,r][z]=Exp[r z]*(((z)^(M))+Sum[c[k]*((z)^(k)),{k,0,M-1}]);
Q\[Dagger][1][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q\[Dagger][1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly\[Dagger][c1\[Dagger],-rep1\[LeftDoubleBracket]4\[RightDoubleBracket]-rep2\[LeftDoubleBracket]4\[RightDoubleBracket]+repAll\[LeftDoubleBracket]1\[RightDoubleBracket],\[Phi]1][z];
Q\[Dagger][2][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q\[Dagger][2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly\[Dagger][c2\[Dagger],-rep1\[LeftDoubleBracket]4\[RightDoubleBracket]-rep2\[LeftDoubleBracket]4\[RightDoubleBracket]+repAll\[LeftDoubleBracket]2\[RightDoubleBracket],\[Phi]2][z];
Q\[Dagger][3][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q\[Dagger][3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly\[Dagger][c3\[Dagger],-rep1\[LeftDoubleBracket]4\[RightDoubleBracket]-rep2\[LeftDoubleBracket]4\[RightDoubleBracket]+repAll\[LeftDoubleBracket]3\[RightDoubleBracket],\[Phi]3][z];
Q\[Dagger][4][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q\[Dagger][4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly\[Dagger][c4\[Dagger],-rep1\[LeftDoubleBracket]4\[RightDoubleBracket]-rep2\[LeftDoubleBracket]4\[RightDoubleBracket]+repAll\[LeftDoubleBracket]4\[RightDoubleBracket],\[Phi]4][z];


(* === Baxter equations === *)

ClearAll[Bax\[Dagger]];
Bax\[Dagger][1][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=Bax\[Dagger][1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=((BaxOp\[Dagger][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][Q\[Dagger][1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[-\[Phi]1 z])//Simplify;
Bax\[Dagger][2][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=Bax\[Dagger][2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=((BaxOp\[Dagger][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][Q\[Dagger][2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[-\[Phi]2 z])//Simplify;
Bax\[Dagger][3][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=Bax\[Dagger][3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=((BaxOp\[Dagger][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][Q\[Dagger][3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[-\[Phi]3 z])//Simplify;
Bax\[Dagger][4][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=Bax\[Dagger][4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=((BaxOp\[Dagger][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][Q\[Dagger][4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[-\[Phi]4 z])//Simplify;

(* Unknown coefficients *)

ClearAll[Vars1\[Dagger],Vars2\[Dagger],Vars3\[Dagger],Vars4\[Dagger]];
Vars1\[Dagger][rep1_,rep2_,repAll_]:=Vars1\[Dagger][rep1,rep2,repAll]=Table[c1\[Dagger][n],{n,0,-rep1\[LeftDoubleBracket]4\[RightDoubleBracket]-rep2\[LeftDoubleBracket]4\[RightDoubleBracket]+repAll\[LeftDoubleBracket]1\[RightDoubleBracket]-1}];
Vars2\[Dagger][rep1_,rep2_,repAll_]:=Vars2\[Dagger][rep1,rep2,repAll]=Table[c2\[Dagger][n],{n,0,-rep1\[LeftDoubleBracket]4\[RightDoubleBracket]-rep2\[LeftDoubleBracket]4\[RightDoubleBracket]+repAll\[LeftDoubleBracket]2\[RightDoubleBracket]-1}];
Vars3\[Dagger][rep1_,rep2_,repAll_]:=Vars3\[Dagger][rep1,rep2,repAll]=Table[c3\[Dagger][n],{n,0,-rep1\[LeftDoubleBracket]4\[RightDoubleBracket]-rep2\[LeftDoubleBracket]4\[RightDoubleBracket]+repAll\[LeftDoubleBracket]3\[RightDoubleBracket]-1}];
Vars4\[Dagger][rep1_,rep2_,repAll_]:=Vars4\[Dagger][rep1,rep2,repAll]=Table[c4\[Dagger][n],{n,0,-rep1\[LeftDoubleBracket]4\[RightDoubleBracket]-rep2\[LeftDoubleBracket]4\[RightDoubleBracket]+repAll\[LeftDoubleBracket]4\[RightDoubleBracket]-1}];


(* === Solving Baxter equations giving hs === *)

ClearAll[SolveBaxh\[Dagger],SolveBaxAllh\[Dagger]];
SolveBaxh\[Dagger][1][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,h24_,h33_,h34_,h36_,h43_,h44_,h45_,h46_,h48_]:=SolveBaxh\[Dagger][1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h24,h33,h34,h36,h43,h44,h45,h46,h48]=Solve[CoefficientList[Bax\[Dagger][1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]/.{h[2,4]->h24,h[3,3]->h33,h[3,4]->h34,h[3,6]->h36,
h[4,3]->h43,h[4,4]->h44,h[4,5]->h45,h[4,6]->h46,h[4,8]->h48}//Together//Numerator,z]==0,Vars1\[Dagger][rep1,rep2,repAll]]//Simplify//Flatten;

SolveBaxh\[Dagger][2][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,h24_,h33_,h34_,h36_,h43_,h44_,h45_,h46_,h48_]:=SolveBaxh\[Dagger][2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h24,h33,h34,h36,h43,h44,h45,h46,h48]=Solve[CoefficientList[Bax\[Dagger][2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]/.{h[2,4]->h24,h[3,3]->h33,h[3,4]->h34,h[3,6]->h36,
h[4,3]->h43,h[4,4]->h44,h[4,5]->h45,h[4,6]->h46,h[4,8]->h48}//Together//Numerator,z]==0,Vars2\[Dagger][rep1,rep2,repAll]]//Simplify//Flatten;

SolveBaxh\[Dagger][3][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,h24_,h33_,h34_,h36_,h43_,h44_,h45_,h46_,h48_]:=SolveBaxh\[Dagger][3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h24,h33,h34,h36,h43,h44,h45,h46,h48]=Solve[CoefficientList[Bax\[Dagger][3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]/.{h[2,4]->h24,h[3,3]->h33,h[3,4]->h34,h[3,6]->h36,
h[4,3]->h43,h[4,4]->h44,h[4,5]->h45,h[4,6]->h46,h[4,8]->h48}//Together//Numerator,z]==0,Vars3\[Dagger][rep1,rep2,repAll]]//Simplify//Flatten;

SolveBaxh\[Dagger][4][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,h24_,h33_,h34_,h36_,h43_,h44_,h45_,h46_,h48_]:=SolveBaxh\[Dagger][4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h24,h33,h34,h36,h43,h44,h45,h46,h48]=Solve[CoefficientList[Bax\[Dagger][4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]/.{h[2,4]->h24,h[3,3]->h33,h[3,4]->h34,h[3,6]->h36,
h[4,3]->h43,h[4,4]->h44,h[4,5]->h45,h[4,6]->h46,h[4,8]->h48}//Together//Numerator,z]==0,Vars4\[Dagger][rep1,rep2,repAll]]//Simplify//Flatten;

SolveBaxAllh\[Dagger][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,h24_,h33_,h34_,h36_,h43_,h44_,h45_,h46_,h48_]:=SolveBaxAllh\[Dagger][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h24,h33,h34,h36,h43,h44,h45,h46,h48]=Join[SolveBaxh\[Dagger][1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h24,h33,h34,h36,h43,h44,h45,h46,h48],SolveBaxh\[Dagger][2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h24,h33,h34,h36,h43,h44,h45,h46,h48],SolveBaxh\[Dagger][3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h24,h33,h34,h36,h43,h44,h45,h46,h48],SolveBaxh\[Dagger][4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h24,h33,h34,h36,h43,h44,h45,h46,h48]]//Simplify;


(* === Example of use === *)

Solh=(SolveBaxAll[FUNDAMENTAL,FUNDAMENTAL,FUNDAMENTAL+FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1]//Simplify)[[1]]

SolveBaxh\[Dagger][1][FUNDAMENTAL,FUNDAMENTAL,FUNDAMENTAL+FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h[2,4]/.Solh,h[3,3]/.Solh,h[3,4]/.Solh,h[3,6]/.Solh,
h[4,3]/.Solh,h[4,4]/.Solh,h[4,5]/.Solh,h[4,6]/.Solh,h[4,8]/.Solh]
SolveBaxh\[Dagger][2][FUNDAMENTAL,FUNDAMENTAL,FUNDAMENTAL+FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h[2,4]/.Solh,h[3,3]/.Solh,h[3,4]/.Solh,h[3,6]/.Solh,
h[4,3]/.Solh,h[4,4]/.Solh,h[4,5]/.Solh,h[4,6]/.Solh,h[4,8]/.Solh]
SolveBaxh\[Dagger][3][FUNDAMENTAL,FUNDAMENTAL,FUNDAMENTAL+FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h[2,4]/.Solh,h[3,3]/.Solh,h[3,4]/.Solh,h[3,6]/.Solh,
h[4,3]/.Solh,h[4,4]/.Solh,h[4,5]/.Solh,h[4,6]/.Solh,h[4,8]/.Solh]
SolveBaxh\[Dagger][4][FUNDAMENTAL,FUNDAMENTAL,FUNDAMENTAL+FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h[2,4]/.Solh,h[3,3]/.Solh,h[3,4]/.Solh,h[3,6]/.Solh,
h[4,3]/.Solh,h[4,4]/.Solh,h[4,5]/.Solh,h[4,6]/.Solh,h[4,8]/.Solh]


Print["\<Hi\>"]
SolveBaxAllh\[Dagger][FUNDAMENTAL,FUNDAMENTAL,FUNDAMENTAL+FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h[2,4]/.Solh,h[3,3]/.Solh,h[3,4]/.Solh,h[3,6]/.Solh,
h[4,3]/.Solh,h[4,4]/.Solh,h[4,5]/.Solh,h[4,6]/.Solh,h[4,8]/.Solh]


(* === Q Hodge === *)

ClearAll[Ng,SetN];
Ng=4;
SetN=Range[Ng];

ClearAll[Hodge];
Hodge[Qfour_,rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,z_,i_]:=Hodge[Qfour,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,z,i]=
(Exp[z*(\[Phi]1+\[Phi]2+\[Phi]3+\[Phi]4)]*Signature[Join[Complement[SetN,{i}],{i}]]*Wronskian[{Qfour[[Complement[SetN,{i}]]][[1]],Qfour[[Complement[SetN,{i}]]][[2]],Qfour[[Complement[SetN,{i}]]][[3]]},z]*(((((z))^(-(2*rep1[[1]]-rep1[[2]]-rep1[[3]])))*(((z-1))^(-(2*rep2[[1]]-rep2[[2]]-rep2[[3]])))))/((({\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[1]]]]-{\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[2]]]]))*({\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[1]]]]-{\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[3]]]])*({\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[2]]]]-{\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[3]]]])))//Simplify;
(*The Q argument should be Q1[z], not Q1!*)

(**(((((z))^(-2(rep1[[1]]-rep1[[2]])))*(((z-1))^(-2(rep2[[1]]-rep2[[2]])))))*)
(*/({\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[1]]]]-{\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[2]]]])*)


(* === Check Q Hodge === *)


(* === Check Table === *)

ClearAll[rep1T,rep2T,repAllT]
rep1T={3,0,0,0};
rep2T={4,3,2,1};
repAllT[i_]:=FindGoodWeight2[rep1T,rep2T][[i]];
\[Phi]1T=1;
\[Phi]2T=((1)/(2));
\[Phi]3T=((1)/(3));
\[Phi]4T=((1)/(5));

(*rep1T=ADJOINT;
rep2T=ADJOINT;
repAllT[1]={4,2,0};
repAllT[2]={3,2,1};
repAllT[3]={2,2,2};
\[Phi]1T=3;
\[Phi]2T=5;
\[Phi]3T=7;*)

(*rep1T={3,2,1};
rep2T={2,0,0};
repAllT[1]={5,2,1};
repAllT[2]={4,3,1};
repAllT[3]={4,2,2};
\[Phi]1T=1;
\[Phi]2T=((1)/(2));
\[Phi]3T=((1)/(7));*)

(*rep1T=FUNDAMENTAL;
rep2T=FUNDAMENTAL;
repAllT[1]={1,1,0};
repAllT[2]={2,0,0};
repAllT[3]={1,1,0};*)

Table[
Solh=(SolveBaxAll[rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,1,1]//Simplify)[[1]];

Qh[1]=Q[1][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh;
Qh[2]=Q[2][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh;
Qh[3]=Q[3][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh;
Qh[4]=Q[4][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh;

(*SolveBaxh\[Dagger][1][rep1T,rep2T,repAllT,\[Phi]1,\[Phi]2,\[Phi]3,h[2,5]/.Solh,h[3,5]/.Solh,h[3,7]/.Solh,h[3,9]/.Solh];
SolveBaxh\[Dagger][2][rep1T,rep2T,repAllT,\[Phi]1,\[Phi]2,\[Phi]3,h[2,5]/.Solh,h[3,5]/.Solh,h[3,7]/.Solh,h[3,9]/.Solh];
SolveBaxh\[Dagger][3][rep1T,rep2T,repAllT,\[Phi]1,\[Phi]2,\[Phi]3,h[2,5]/.Solh,h[3,5]/.Solh,h[3,7]/.Solh,h[3,9]/.Solh];*)


Print["\<Hi\>"];
Sol\[Dagger]=SolveBaxAllh\[Dagger][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,h[2,4]/.Solh,h[3,3]/.Solh,h[3,4]/.Solh,h[3,6]/.Solh,
h[4,3]/.Solh,h[4,4]/.Solh,h[4,5]/.Solh,h[4,6]/.Solh,h[4,8]/.Solh];

QH[1]=Q\[Dagger][1][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Sol\[Dagger];
QH[2]=Q\[Dagger][2][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Sol\[Dagger];
QH[3]=Q\[Dagger][3][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Sol\[Dagger];
QH[4]=Q\[Dagger][4][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Sol\[Dagger];

(*Print["\<Hodge\>"]
Hodge[{Qh[1],Qh[2],Qh[3]},rep1T,rep2T,\[Phi]1,\[Phi]2,\[Phi]3,z,2,1,1]//Simplify

Print["\<\!\(\*SuperscriptBox[\(Q\), \(\[Dagger]\)]\)\>"]
QH[2]//Simplify*)

Print["\<M=\>",M,"\<, \>",Table[(((1)/(QH[t]))Hodge[{Qh[1],Qh[2],Qh[3],Qh[4]},rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,z,t])//Simplify,{t,1,4}]],{M,1,4}];


(* === Testing === *)

ClearAll[rep1T,rep2T,repAllT]
rep1T=ADJOINT;
rep2T=FUNDAMENTAL;
repAllT[i_]:=FindGoodWeight2[rep1T,rep2T][[i]];
\[Phi]1T=1;
\[Phi]2T=((1)/(2));
\[Phi]3T=((1)/(5));
\[Phi]4T=((1)/(7));
M=3;
M2=2;
Solh=(SolveBaxAll[rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,1,1]//Simplify)[[1]];
Solh2=(SolveBaxAll[rep1T,rep2T,repAllT[M2],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,1,1]//Simplify)[[1]];

Qh[1]=Q[1][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh;
Qh[2]=Q[2][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh;
Qh[3]=Q[3][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh;
Qh[4]=Q[4][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh;

Qh2[1]=Q[1][rep1T,rep2T,repAllT[M2],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh2;
Qh2[2]=Q[2][rep1T,rep2T,repAllT[M2],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh2;
Qh2[3]=Q[3][rep1T,rep2T,repAllT[M2],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh2;
Qh2[4]=Q[4][rep1T,rep2T,repAllT[M2],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh2;

(*SolveBaxh\[Dagger][1][rep1T,rep2T,repAllT,\[Phi]1,\[Phi]2,\[Phi]3,h[2,5]/.Solh,h[3,5]/.Solh,h[3,7]/.Solh,h[3,9]/.Solh];
SolveBaxh\[Dagger][2][rep1T,rep2T,repAllT,\[Phi]1,\[Phi]2,\[Phi]3,h[2,5]/.Solh,h[3,5]/.Solh,h[3,7]/.Solh,h[3,9]/.Solh];
SolveBaxh\[Dagger][3][rep1T,rep2T,repAllT,\[Phi]1,\[Phi]2,\[Phi]3,h[2,5]/.Solh,h[3,5]/.Solh,h[3,7]/.Solh,h[3,9]/.Solh];*)


Print["\<Hi\>"];
Sol\[Dagger]=SolveBaxAllh\[Dagger][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,h[2,4]/.Solh,h[3,3]/.Solh,h[3,4]/.Solh,h[3,6]/.Solh,
h[4,3]/.Solh,h[4,4]/.Solh,h[4,5]/.Solh,h[4,6]/.Solh,h[4,8]/.Solh];

QH[1]=Q\[Dagger][1][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Sol\[Dagger];
QH[2]=Q\[Dagger][2][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Sol\[Dagger];
QH[3]=Q\[Dagger][3][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Sol\[Dagger];
QH[4]=Q\[Dagger][4][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Sol\[Dagger];

(*Print["\<Hodge\>"]
Hodge[{Qh[1],Qh[2],Qh[3]},rep1T,rep2T,\[Phi]1,\[Phi]2,\[Phi]3,z,2,1,1]//Simplify

Print["\<\!\(\*SuperscriptBox[\(Q\), \(\[Dagger]\)]\)\>"]
QH[2]//Simplify*)

Table[(((1)/(QH[t]))Hodge[{Qh[1],Qh[2],Qh[3],Qh[4]},rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,z,t])//Simplify,{t,1,4}]



Print["\<Hello!\>"];
Q3[{Qh[1],Qh[2],Qh[3],Qh[4]},z,1,2,3];
Hodge[{Qh[1],Qh[2],Qh[3],Qh[4]},rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,z,4]
QH[4]//Simplify

QH[1]//Simplify
Wronskian[{Qh[2],Qh[3],Qh[4]},z]*Exp[(\[Phi]1T+\[Phi]2T+\[Phi]3T+\[Phi]4T)*z]//Simplify

QHT={QH[1],QH[2],QH[3],QH[4]};
QhT={Qh[1],Qh[2],Qh[3],Qh[4]};
Qh2T={Qh2[1],Qh2[2],Qh2[3],Qh2[4]};
OrthogRen[QHT,QhT,QHT,QhT,rep1T,rep2T]
OrthogRen[QHT,Qh2T,QHT,Qh2T,rep1T,rep2T]


(* === Renormalised scalar product === *)

(* TextData[{Input: , Cell[Q={Subscript[Q, 1],Subscript[Q, 2],Subscript[Q, 3]}, ExpressionUUID->   "3d028b11-383d-4a70-a660-547d034fcdc1"],  , Cell[QH={Subscript[QH, 1],Subscript[QH, 2],Subscript[QH, 3]}, ExpressionUUID->   "11d676a9-beb9-4c5b-a3c6-194f29502c4c"]}] *)

ClearAll[Normalisation];
Normalisation[QHL_,QL_,QHR_,QR_,rep1_,rep2_]:=Normalisation[QHL,QL,QHR,QR,rep1,rep2]=(Limit[QHL[[1]],z->0,Direction->"\<FromAbove\>"])*(Limit[Wronskian[{QHL[[1]],QHL[[2]]},z]*((1)/(((z)^(rep1[[3]]-rep1[[4]]))*(((z-1))^(rep2[[3]]-rep2[[4]])))),z->0,Direction->"\<FromAbove\>"])*(Limit[QL[[4]],z->0,Direction->"\<FromAbove\>"])*(Limit[QR[[1]],z->0,Direction->"\<FromAbove\>"])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z]*((1)/(((z)^(rep1[[1]]-rep1[[2]]))*(((z-1))^(rep2[[1]]-rep2[[2]])))),z->0,Direction->"\<FromAbove\>"])*(Limit[QHR[[4]],z->0,Direction->"\<FromAbove\>"])*
(Limit[QHL[[1]],z->1,Direction->"\<FromAbove\>"])*(Limit[Wronskian[{QHL[[1]],QHL[[2]]},z]*((1)/(((z)^(rep1[[3]]-rep1[[4]]))*(((z-1))^(rep2[[3]]-rep2[[4]])))),z->1,Direction->"\<FromAbove\>"])*(Limit[QL[[4]],z->1,Direction->"\<FromAbove\>"])*(Limit[QR[[1]],z->1,Direction->"\<FromAbove\>"])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z]*((1)/(((z)^(rep1[[1]]-rep1[[2]]))*(((z-1))^(rep2[[1]]-rep2[[2]])))),z->1,Direction->"\<FromAbove\>"])*(Limit[QHR[[4]],z->1,Direction->"\<FromAbove\>"]);

ClearAll[OrthogRen,RawOrthogT,RawOrthog12];
RawOrthogT[QH_,Q_,rep1_,rep2_,w_,\[Alpha]_,i_,j_]:=RawOrthogT[QH,Q,rep1,rep2,w,\[Alpha],i,j]={MeanValue[QH[[i]]*D[Q[[j]],{z,2}]*((1)/(z-w)),\[Alpha],rep1,rep2],MeanValue[QH[[i]]*D[Q[[j]],{z,1}]*((1)/((((z-w))^(2)))),\[Alpha],rep1,rep2],MeanValue[QH[[i]]*D[Q[[j]],{z,1}]*((1)/(z-w)),\[Alpha],rep1,rep2],MeanValue[QH[[i]]*Q[[j]]*((1)/((((z-w))^(3)))),\[Alpha],rep1,rep2],MeanValue[QH[[i]]*Q[[j]]*((1)/((((z-w))^(2)))),\[Alpha],rep1,rep2],MeanValue[QH[[i]]*Q[[j]]*((1)/(z-w)),\[Alpha],rep1,rep2]};

RawOrthog12[QH_,Q_,rep1_,rep2_,\[Alpha]_,i_,j_]:=RawOrthog12[QH,Q,rep1,rep2,\[Alpha],i,j]=Join[RawOrthogT[QH,Q,rep1,rep2,0,\[Alpha],i,j],RawOrthogT[QH,Q,rep1,rep2,1,\[Alpha],i,j]];

OrthogRen[QHL_,QL_,QHR_,QR_,rep1_,rep2_]:=OrthogRen[QHL,QL,QHR,QR,rep1,rep2]=(1/Normalisation[QHL,QL,QHR,QR,rep1,rep2])*Det[{RawOrthog12[QHL,QR,rep1,rep2,0,1,1],RawOrthog12[QHL,QR,rep1,rep2,0,1,2],RawOrthog12[QHL,QR,rep1,rep2,0,2,1],RawOrthog12[QHL,QR,rep1,rep2,0,2,2],RawOrthog12[QHL,QR,rep1,rep2,0,1,3],RawOrthog12[QHL,QR,rep1,rep2,0,3,1],RawOrthog12[QHL,QR,rep1,rep2,1,1,1],RawOrthog12[QHL,QR,rep1,rep2,1,1,2],RawOrthog12[QHL,QR,rep1,rep2,1,2,1],RawOrthog12[QHL,QR,rep1,rep2,1,2,2],RawOrthog12[QHL,QR,rep1,rep2,1,1,3],RawOrthog12[QHL,QR,rep1,rep2,1,3,1]}];


(* === TEST === *)

ClearAll[rep1T,rep2T,repAllT,M,M2,Solh,Solh2,Qh,Qh2,QH]
rep1T=ADJOINT;
rep2T=ADJOINT;
repAllT[i_]:=FindGoodWeight2[rep1T,rep2T][[i]];
\[Phi]1T=1;
\[Phi]2T=((1)/(2));
\[Phi]3T=((1)/(5));
\[Phi]4T=((1)/(7));
M=1;
M2=3;
Solh=(SolveBaxAll[rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,1,1]//Simplify)[[1]];
Solh2=(SolveBaxAll[rep1T,rep2T,repAllT[M2],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,1,1]//Simplify)[[1]];

Qh[1]=Q[1][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh;
Qh[2]=Q[2][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh;
Qh[3]=Q[3][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh;
Qh[4]=Q[4][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh;

Qh2[1]=Q[1][rep1T,rep2T,repAllT[M2],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh2;
Qh2[2]=Q[2][rep1T,rep2T,repAllT[M2],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh2;
Qh2[3]=Q[3][rep1T,rep2T,repAllT[M2],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh2;
Qh2[4]=Q[4][rep1T,rep2T,repAllT[M2],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Solh2;

(*SolveBaxh\[Dagger][1][rep1T,rep2T,repAllT,\[Phi]1,\[Phi]2,\[Phi]3,h[2,5]/.Solh,h[3,5]/.Solh,h[3,7]/.Solh,h[3,9]/.Solh];
SolveBaxh\[Dagger][2][rep1T,rep2T,repAllT,\[Phi]1,\[Phi]2,\[Phi]3,h[2,5]/.Solh,h[3,5]/.Solh,h[3,7]/.Solh,h[3,9]/.Solh];
SolveBaxh\[Dagger][3][rep1T,rep2T,repAllT,\[Phi]1,\[Phi]2,\[Phi]3,h[2,5]/.Solh,h[3,5]/.Solh,h[3,7]/.Solh,h[3,9]/.Solh];*)


Print["\<Hi\>"];
Sol\[Dagger]=SolveBaxAllh\[Dagger][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,h[2,4]/.Solh,h[3,3]/.Solh,h[3,4]/.Solh,h[3,6]/.Solh,
h[4,3]/.Solh,h[4,4]/.Solh,h[4,5]/.Solh,h[4,6]/.Solh,h[4,8]/.Solh];

QH[1]=Q\[Dagger][1][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Sol\[Dagger];
QH[2]=Q\[Dagger][2][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Sol\[Dagger];
QH[3]=Q\[Dagger][3][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Sol\[Dagger];
QH[4]=Q\[Dagger][4][rep1T,rep2T,repAllT[M],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][z]/.Sol\[Dagger];

(*Print["\<Hodge\>"]
Hodge[{Qh[1],Qh[2],Qh[3]},rep1T,rep2T,\[Phi]1,\[Phi]2,\[Phi]3,z,2,1,1]//Simplify

Print["\<\!\(\*SuperscriptBox[\(Q\), \(\[Dagger]\)]\)\>"]
QH[2]//Simplify*)

Table[(((1)/(QH[t]))Hodge[{Qh[1],Qh[2],Qh[3],Qh[4]},rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,z,t])//Simplify,{t,1,4}];


OrthogRen[{QH[1],QH[2],QH[3],QH[4]},{Qh[1],Qh[2],Qh[3],Qh[4]},{QH[1],QH[2],QH[3],QH[4]},{Qh[1],Qh[2],Qh[3],Qh[4]},rep1T,rep2T]
OrthogRen[{QH[1],QH[2],QH[3],QH[4]},{Qh2[1],Qh2[2],Qh2[3],Qh2[4]},{QH[1],QH[2],QH[3],QH[4]},{Qh2[1],Qh2[2],Qh2[3],Qh2[4]},rep1T,rep2T]

(*
Print["\<Normalisation\>"];
Normalisation[{QH[1],QH[2],QH[3],QH[4]},{Qh[1],Qh[2],Qh[3],Qh[4]},{QH[1],QH[2],QH[3],QH[4]},{Qh[1],Qh[2],Qh[3],Qh[4]},rep1T,rep2T]
Normalisation[{QH[1],QH[2],QH[3],QH[4]},{Qh2[1],Qh2[2],Qh2[3],Qh2[4]},{QH[1],QH[2],QH[3],QH[4]},{Qh2[1],Qh2[2],Qh2[3],Qh2[4]},rep1T,rep2T]
*)

(*
Print["\<Norm0\>"]
(Limit[QH[1],z->0,Direction->"\<FromAbove\>"])
(Limit[Wronskian[{QH[1],QH[2]},z]*((1)/(((z)^(rep1T[[2]]-rep1T[[3]]))*(((z-1))^(rep2T[[2]]-rep2T[[3]])))),z->0,Direction->"\<FromAbove\>"])
(Limit[Qh[1],z->0,Direction->"\<FromAbove\>"])
(Limit[Qh[1],z->0,Direction->"\<FromAbove\>"])
(Limit[Wronskian[{Qh[1],Qh[2]},z]*((1)/(((z)^(rep1T[[1]]-rep1T[[2]]))*(((z-1))^(rep2T[[1]]-rep2T[[2]])))),z->0,Direction->"\<FromAbove\>"])
(Limit[QH[4],z->0,Direction->"\<FromAbove\>"])

Print["\<Norm1\>"]
(Limit[QH[1],z->1,Direction->"\<FromAbove\>"])
(Limit[Wronskian[{QH[1],QH[2]},z]*((1)/(((z)^(rep1T[[2]]-rep1T[[3]]))*(((z-1))^(rep2T[[2]]-rep2T[[3]])))),z->1,Direction->"\<FromAbove\>"])
(Limit[Qh[1],z->1,Direction->"\<FromAbove\>"])
(Limit[Qh[1],z->1,Direction->"\<FromAbove\>"])
(Limit[Wronskian[{Qh[1],Qh[2]},z]*((1)/(((z)^(rep1T[[1]]-rep1T[[2]]))*(((z-1))^(rep2T[[1]]-rep2T[[2]])))),z->1,Direction->"\<FromAbove\>"])
(Limit[QH[4],z->1,Direction->"\<FromAbove\>"])
*)

(*
Print["\<Reg Wr Q\>"];
Wronskian[{Qh[1],Qh[2]},z]

Print["\<Reg Wr QH\>"];
Wronskian[{QH[1],QH[2]},z]
*)


(* === Renormalised scalar product === *)

(* TextData[{Input: , Cell[Q={Subscript[Q, 1],Subscript[Q, 2],Subscript[Q, 3]}, ExpressionUUID->   "c3fd177f-bbcf-4744-86d0-9eb57821b3b0"],  , Cell[QH={Subscript[QH, 1],Subscript[QH, 2],Subscript[QH, 3]}, ExpressionUUID->   "f0210aab-2898-432d-8d93-c1b6dab927a3"]}] *)


(* === rep1 === *)

ClearAll[Normalisation1];
Normalisation1[QHL_,QL_,QHR_,QR_,rep1_]:=Normalisation1[QHL,QL,QHR,QR,rep1]=(Limit[QHL[[1]],z->0,Direction->"\<FromAbove\>"])*(Limit[Wronskian[{QHL[[1]],QHL[[2]]},z]*((1)/(((z)^(rep1[[3]]-rep1[[4]])))),z->0,Direction->"\<FromAbove\>"])*(Limit[QL[[4]],z->0,Direction->"\<FromAbove\>"])*(Limit[QR[[1]],z->0,Direction->"\<FromAbove\>"])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z]*((1)/(((z)^(rep1[[1]]-rep1[[2]])))),z->0,Direction->"\<FromAbove\>"])*(Limit[QHR[[4]],z->0,Direction->"\<FromAbove\>"]);

ClearAll[RawOrthog1];
RawOrthog1[QH_,Q_,rep1_,w_,\[Alpha]_,i_,j_]:=RawOrthog1[QH,Q,rep1,w,\[Alpha],i,j]={MeanValue1[QH[[i]]*D[Q[[j]],{z,2}]*((1)/(z-w)),\[Alpha],rep1],MeanValue1[QH[[i]]*D[Q[[j]],{z,1}]*((1)/((((z-w))^(2)))),\[Alpha],rep1],MeanValue1[QH[[i]]*D[Q[[j]],{z,1}]*((1)/(z-w)),\[Alpha],rep1],MeanValue1[QH[[i]]*Q[[j]]*((1)/((((z-w))^(3)))),\[Alpha],rep1],MeanValue1[QH[[i]]*Q[[j]]*((1)/((((z-w))^(2)))),\[Alpha],rep1],MeanValue1[QH[[i]]*Q[[j]]*((1)/(z-w)),\[Alpha],rep1]};
ClearAll[OrthogRen1];
OrthogRen1[QHL_,QL_,QHR_,QR_,rep1_]:=OrthogRen1[QHL,QL,QHR,QR,rep1]=(1/Normalisation1[QHL,QL,QHR,QR,rep1])*Det[{RawOrthog1[QHL,QR,rep1,0,0,1,1],RawOrthog1[QHL,QR,rep1,0,0,1,2],RawOrthog1[QHL,QR,rep1,0,0,2,1],RawOrthog1[QHL,QR,rep1,0,0,2,2],RawOrthog1[QHL,QR,rep1,0,0,1,3],RawOrthog1[QHL,QR,rep1,0,0,3,1]}];


(* === rep2 === *)

ClearAll[Normalisation2];
Normalisation2[QHL_,QL_,QHR_,QR_,rep2_]:=Normalisation2[QHL,QL,QHR,QR,rep2]=
(Limit[QHL[[1]],z->1,Direction->"\<FromAbove\>"])*(Limit[Wronskian[{QHL[[1]],QHL[[2]]},z]*((1)/((((z-1))^(rep2[[3]]-rep2[[4]])))),z->1,Direction->"\<FromAbove\>"])*(Limit[QL[[4]],z->1,Direction->"\<FromAbove\>"])*(Limit[QR[[1]],z->1,Direction->"\<FromAbove\>"])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z]*((1)/((((z-1))^(rep2[[1]]-rep2[[2]])))),z->1,Direction->"\<FromAbove\>"])*(Limit[QHR[[4]],z->1,Direction->"\<FromAbove\>"]);

ClearAll[RawOrthog2];
RawOrthog2[QH_,Q_,rep2_,w_,\[Alpha]_,i_,j_]:=RawOrthog2[QH,Q,rep2,w,\[Alpha],i,j]={MeanValue2[QH[[i]]*D[Q[[j]],{z,2}]*((1)/(z-w)),\[Alpha],rep2],MeanValue2[QH[[i]]*D[Q[[j]],{z,1}]*((1)/((((z-w))^(2)))),\[Alpha],rep2],MeanValue2[QH[[i]]*D[Q[[j]],{z,1}]*((1)/(z-w)),\[Alpha],rep2],MeanValue2[QH[[i]]*Q[[j]]*((1)/((((z-w))^(3)))),\[Alpha],rep2],MeanValue2[QH[[i]]*Q[[j]]*((1)/((((z-w))^(2)))),\[Alpha],rep2],MeanValue2[QH[[i]]*Q[[j]]*((1)/(z-w)),\[Alpha],rep2]};
ClearAll[OrthogRen2];
OrthogRen2[QHL_,QL_,QHR_,QR_,rep2_]:=OrthogRen2[QHL,QL,QHR,QR,rep2]=(1/Normalisation2[QHL,QL,QHR,QR,rep2])*Det[{RawOrthog2[QHL,QR,rep2,0,0,1,1],RawOrthog2[QHL,QR,rep2,0,0,1,2],RawOrthog2[QHL,QR,rep2,0,0,2,1],RawOrthog2[QHL,QR,rep2,0,0,2,2],RawOrthog2[QHL,QR,rep2,0,0,1,3],RawOrthog2[QHL,QR,rep2,0,0,3,1]}];


(* === RUN: \[CapitalPsi] space === *)


(* === Lax connection === *)

ClearAll[LIndex];
LIndex[a_,b_,rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((1)/(z-w1))*KroneckerProduct[EM[a,b,rep1],Id[rep2]]+((1)/(z-w2))*KroneckerProduct[Id[rep1],EM[a,b,rep2]]+KroneckerDelta[a,1]*KroneckerDelta[b,1]*\[Phi]1*KroneckerProduct[Id[rep1],Id[rep2]]+KroneckerDelta[a,2]*KroneckerDelta[b,2]*\[Phi]2*KroneckerProduct[Id[rep1],Id[rep2]]+KroneckerDelta[a,3]*KroneckerDelta[b,3]*\[Phi]3*KroneckerProduct[Id[rep1],Id[rep2]]+KroneckerDelta[a,4]*KroneckerDelta[b,4]*\[Phi]4*KroneckerProduct[Id[rep1],Id[rep2]];


(* === Example of use === *)

ClearAll[rep1,rep2];
rep1=FUNDAMENTAL;
rep2=FUNDAMENTAL;
LIndex[1,1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal


(* === H1 === *)


(* === TextData[{Regulated Hamiltonian , Cell[tr(L), ExpressionUUID->   "a5b2117a-86bf-42e4-926d-43fae8c35d94"]}] === *)

ClearAll[H1];
H1[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=H1[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{a,1,Ng}]//Simplify;


(* === Definition Residue of H === *)

Clear[ResH1];
ResH1[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_]:=ResH1[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w]=SeriesCoefficient[H1[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,-1}]


(* === Example of use === *)

H1[FUNDAMENTAL,FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal;


(* === H2 === *)


(* === TextData[{Regulated Hamiltonian , Cell[((3)/(2))*tr(L')+trA(((L)^(2))), ExpressionUUID->   "df3421e7-3a88-4156-b6d2-09ca81d354ee"]}] === *)


(* === Antisymm. def === *)

ClearAll[trAL2];
trAL2[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trAL2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((1)/(2))*(Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[b,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[b,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{a,1,Ng},{b,1,Ng}]);


(* === Def. H2 === *)

ClearAll[H2];
H2[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=H2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=(trAL2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+((3)/(2))*Sum[D[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z],{a,1,Ng}])//Simplify


(* === Definition Residue of H === *)

Clear[ResH2];
ResH2[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_]:=ResH2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w]=SeriesCoefficient[H2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,-1}]

Clear[ResH2O];
ResH2O[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_,o_]:=ResH2O[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w,o]=SeriesCoefficient[H2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,o}]


(* === Example of use === *)

H2[FUNDAMENTAL,FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal;


(* === H3 === *)


(* === H3 Unsimplified === *)


(* === Antisymm. def === *)

ClearAll[trAL3];
trAL3[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trAL3[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((1)/(6))*Sum[LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{i1,1,Ng},{i2,1,Ng},{i3,1,Ng}];

ClearAll[trALpL,trALLp];
trALLp[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trALLp[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((1)/(2))*(Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex[b,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]-LIndex[b,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z],{a,1,Ng},{b,1,Ng}]);
trALpL[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trALpL[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((1)/(2))*(Sum[D[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex[b,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-D[LIndex[b,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{a,1,Ng},{b,1,Ng}]);


(* === Def H3 === *)

ClearAll[H3U];
H3U[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=H3U[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=trAL3[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+((4)/(3))*trALLp[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+((2)/(3))*trALpL[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+Sum[D[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{z,2}],{a,1,Ng}]


(* === H3 Simplified === *)

ClearAll[H3];
H3[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=H3[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=H3U[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]


(* === Definition Residue of H === *)

Clear[ResH3];
ResH3[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_]:=ResH3[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w]=SeriesCoefficient[H3[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,-1}]

Clear[ResH3O];
ResH3O[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_,o_]:=ResH3O[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w,o]=SeriesCoefficient[H3[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,o}]


(* === Example of use === *)

H3[FUNDAMENTAL,FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal//Simplify


(* === H4 === *)


(* === H4 Unsimplified === *)


(* === Antisymm. def === *)

ClearAll[trAL4];
trAL4[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trAL4[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((1)/(4!))*Sum[LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i4,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i4,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i3,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i4,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i4,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i2,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i4,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i3,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i4,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i2,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i4,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i4,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i3,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i4,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i4,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i1,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i4,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i3,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i4,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i1,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i4,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i4,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i2,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i4,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i4,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i1,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i4,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i2,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i4,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i1,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i4,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i3,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i4,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i2,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i4,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i3,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i4,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i1,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i4,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i2,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i4,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i1,i4,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{i1,1,Ng},{i2,1,Ng},{i3,1,Ng},{i4,1,Ng}];

(*((L)^(3))*)
ClearAll[trALLLp,trALLpL,trALpLL];
trALLLp[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trALLLp[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((1)/(6))*Sum[LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]-LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]-LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]+LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]+LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]-LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z],{i1,1,Ng},{i2,1,Ng},{i3,1,Ng}];

trALLpL[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trALLpL[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((1)/(6))*Sum[LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{i1,1,Ng},{i2,1,Ng},{i3,1,Ng}];

trALpLL[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trALpLL[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((1)/(6))*Sum[D[LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-D[LIndex[i1,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-D[LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i3,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+D[LIndex[i2,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex[i3,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+D[LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex[i1,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i2,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-D[LIndex[i3,i1,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex[i2,i2,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex[i1,i3,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{i1,1,Ng},{i2,1,Ng},{i3,1,Ng}];

(*((L)^(2))*)
ClearAll[trALLpp,trALpLp,trALppL];
trALLpp[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trALLpp[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((1)/(2))*(Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex[b,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{z,2}]-LIndex[b,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{z,2}],{a,1,Ng},{b,1,Ng}]);

trALpLp[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trALpLp[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((1)/(2))*(Sum[D[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].D[LIndex[b,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]-D[LIndex[b,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].D[LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z],{a,1,Ng},{b,1,Ng}]);

trALppL[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trALppL[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((1)/(2))*(Sum[D[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{z,2}].LIndex[b,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-D[LIndex[b,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{z,2}].LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{a,1,Ng},{b,1,Ng}]);


(* === Def H4 === *)

ClearAll[H4U];
H4U[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=H4U[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=trAL4[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+((3)/(4))*trALLLp[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+((1)/(2))*trALLpL[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+((1)/(4))*trALpLL[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+
((1)/(2))*trALLpp[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+((1)/(2))*trALpLp[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+((1)/(6))*trALppL[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+
((1)/(4))*Sum[D[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{z,3}],{a,1,Ng}]


(* === H4 Simplified === *)

ClearAll[H4];
H4[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=H4[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=H4U[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]


(* === Definition Residue of H === *)

Clear[ResH4];
ResH4[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_]:=ResH4[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w]=SeriesCoefficient[H4[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,-1}]

Clear[ResH4O];
ResH4O[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_,o_]:=ResH4O[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w,o]=SeriesCoefficient[H4[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,o}]


(* === Example of use === *)

H4[FUNDAMENTAL,FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal//Simplify;


(* === Check eigenvalues === *)

UPR[rep1_,rep2_]:=UPR[rep1,rep2]=Table[{KroneckerDelta[j,1]},{j,1,dim[rep1]*dim[rep2]}];

rep1=FUNDAMENTAL;
rep2=ADJOINT;


(* === H2 === *)

H2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].UPR[rep1,rep2]

SeriesCoefficient[((H2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].UPR[rep1,rep2])//Flatten)[[1]],{z,0,-2}]//Normal
SeriesCoefficient[((H2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].UPR[rep1,rep2])//Flatten)[[1]],{z,1,-2}]//Normal

rep1\[LeftDoubleBracket]3\[RightDoubleBracket] (-1+rep1\[LeftDoubleBracket]4\[RightDoubleBracket])+rep1\[LeftDoubleBracket]2\[RightDoubleBracket] (-2+rep1\[LeftDoubleBracket]3\[RightDoubleBracket]+rep1\[LeftDoubleBracket]4\[RightDoubleBracket])+rep1\[LeftDoubleBracket]1\[RightDoubleBracket] (-3+rep1\[LeftDoubleBracket]2\[RightDoubleBracket]+rep1\[LeftDoubleBracket]3\[RightDoubleBracket]+rep1\[LeftDoubleBracket]4\[RightDoubleBracket])

rep2\[LeftDoubleBracket]3\[RightDoubleBracket] (-1+rep2\[LeftDoubleBracket]4\[RightDoubleBracket])+rep2\[LeftDoubleBracket]2\[RightDoubleBracket] (-2+rep2\[LeftDoubleBracket]3\[RightDoubleBracket]+rep2\[LeftDoubleBracket]4\[RightDoubleBracket])+rep2\[LeftDoubleBracket]1\[RightDoubleBracket] (-3+rep2\[LeftDoubleBracket]2\[RightDoubleBracket]+rep2\[LeftDoubleBracket]3\[RightDoubleBracket]+rep2\[LeftDoubleBracket]4\[RightDoubleBracket])


(* === H3 === *)

SeriesCoefficient[((H3U[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].UPR[rep1,rep2])//Flatten)[[1]],{z,0,-3}]//Normal
SeriesCoefficient[((H3U[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].UPR[rep1,rep2])//Flatten)[[1]],{z,1,-3}]//Normal

rep1\[LeftDoubleBracket]2\[RightDoubleBracket] (-1+rep1\[LeftDoubleBracket]3\[RightDoubleBracket]) (-2+rep1\[LeftDoubleBracket]4\[RightDoubleBracket])+rep1\[LeftDoubleBracket]1\[RightDoubleBracket] ((-2+rep1\[LeftDoubleBracket]3\[RightDoubleBracket]) (-3+rep1\[LeftDoubleBracket]4\[RightDoubleBracket])+rep1\[LeftDoubleBracket]2\[RightDoubleBracket] (-4+rep1\[LeftDoubleBracket]3\[RightDoubleBracket]+rep1\[LeftDoubleBracket]4\[RightDoubleBracket]))

rep2\[LeftDoubleBracket]2\[RightDoubleBracket] (-1+rep2\[LeftDoubleBracket]3\[RightDoubleBracket]) (-2+rep2\[LeftDoubleBracket]4\[RightDoubleBracket])+rep2\[LeftDoubleBracket]1\[RightDoubleBracket] ((-2+rep2\[LeftDoubleBracket]3\[RightDoubleBracket]) (-3+rep2\[LeftDoubleBracket]4\[RightDoubleBracket])+rep2\[LeftDoubleBracket]2\[RightDoubleBracket] (-4+rep2\[LeftDoubleBracket]3\[RightDoubleBracket]+rep2\[LeftDoubleBracket]4\[RightDoubleBracket]))


(* === H4 === *)

SeriesCoefficient[((H4U[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].UPR[rep1,rep2])//Flatten)[[1]],{z,0,-4}]//Normal
SeriesCoefficient[((H4U[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].UPR[rep1,rep2])//Flatten)[[1]],{z,1,-4}]//Normal

rep1\[LeftDoubleBracket]1\[RightDoubleBracket] (-1+rep1\[LeftDoubleBracket]2\[RightDoubleBracket]) (-2+rep1\[LeftDoubleBracket]3\[RightDoubleBracket]) (-3+rep1\[LeftDoubleBracket]4\[RightDoubleBracket])
rep2\[LeftDoubleBracket]1\[RightDoubleBracket] (-1+rep2\[LeftDoubleBracket]2\[RightDoubleBracket]) (-2+rep2\[LeftDoubleBracket]3\[RightDoubleBracket]) (-3+rep2\[LeftDoubleBracket]4\[RightDoubleBracket])


(* === H1\[Dagger] === *)


(* === TextData[{Regulated Hamiltonian , Cell[tr(L), ExpressionUUID->   "20d533a6-6470-4bb3-97c3-d9035e757156"]}] === *)

ClearAll[H1\[Dagger]];
H1\[Dagger][rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_][z_]:=H1\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]=H1[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z](*Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}]//Simplify;*)


(* === Definition Residue of H === *)

Clear[ResH1\[Dagger]];
ResH1\[Dagger][rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,w_]:=ResH1\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,w]=SeriesCoefficient[H1\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]//Normal,{z,w,-1}]


(* === Example of use === *)

ClearAll[rep1,rep2];
rep1=FUNDAMENTAL;
rep2=FUNDAMENTAL;
H1\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]//Normal;


(* === H2\[Dagger] === *)


(* === TextData[{Regulated Hamiltonian -, Cell[tr(L')+((1)/(2))(((tr(L)))^(2))-((1)/(2))tr(((L)^(2))), ExpressionUUID->   "4460f527-dedc-454e-b278-5827f745c6a5"]}] === *)

ClearAll[H2\[Dagger]];
H2\[Dagger][rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_][z_]:=H2\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]=(-Sum[D[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],z],{a,1,3}]+((1)/(2))*(Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}]).(Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}])-((1)/(2))*Sum[LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z].LIndex[b,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3},{b,1,3}])//Simplify


(* === Definition Residue of H === *)

Clear[ResH2\[Dagger]];
ResH2\[Dagger][rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,w_]:=ResH2\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,w]=SeriesCoefficient[H2\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]//Normal,{z,w,-1}]


(* === H3\[Dagger] === *)


(* === TextData[{H3\[Dagger] simplified , Cell[((1)/(3))tr(L'')+tr(L'.L)-((2)/(3))tr(L').tr(L)-((1)/(2))tr(L).tr(((L)^(2)))+((1)/(6))((tr(L))^(3))+((1)/(3))tr(((L)^(3))), ExpressionUUID->   "7937eed1-0007-4610-93e3-e1dbcf97a25c"]}] === *)

ClearAll[H3\[Dagger]];
H3\[Dagger][rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_][z_]:=H3\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]=(((1)/(3))*Sum[D[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{z,2}],{a,1,3}]+Sum[D[LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],z].LIndex[b,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3},{b,1,3}]-((2)/(3))*Sum[D[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],z],{a,1,3}].Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}]-((1)/(2))*Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}].Sum[LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z].LIndex[b,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3},{b,1,3}]+((1)/(6))*(Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}]).(Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}]).(Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}])+((1)/(3))*Sum[LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z].LIndex[b,c,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z].LIndex[c,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3},{b,1,3},{c,1,3}])//Simplify;


(* === Definition Residue of H === *)

Clear[ResH3\[Dagger]];
ResH3\[Dagger][rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,w_]:=ResH3\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,w]=SeriesCoefficient[H3[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]//Normal,{z,w,-1}]


(* === Example of use === *)

Series[H3\[Dagger][ADJOINT,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]//Normal,{z,0,-3}]//Simplify//Normal;


(* === Right eigenvectors (column) === *)

ClearAll[VecR,VecR1];
VecR[rep1_,rep2_]:=VecR[rep1,rep2]=Table[{v[j]},{j,1,dim[rep1]*dim[rep2]}];
VecR1[rep1_]:=VecR1[rep1]=Table[{v[j]},{j,1,dim[rep1]}];

ClearAll[EigenR];
EigenR[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=EigenR[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=Table[Quiet[Flatten[VecR[rep1,rep2]/.Solve[
CoefficientList[((H2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].VecR[rep1,rep2]-HSpectrum[2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z][[idx]]* VecR[rep1,rep2]))//Together//Numerator,z]==0&&
CoefficientList[((H3[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].VecR[rep1,rep2]-HSpectrum[3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z][[idx]]* VecR[rep1,rep2]))//Together//Numerator,z]==0&&
CoefficientList[((H4[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].VecR[rep1,rep2]-HSpectrum[4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z][[idx]]* VecR[rep1,rep2]))//Together//Numerator,z]==0,Flatten[VecR[rep1,rep2]]],1]],{idx,1,Dimensions[HSpectrum[1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]][[1]]}]


(* === Example === *)

EigenR[FUNDAMENTAL,FUNDAMENTAL,{1,1,0,0},3,5,7,9]


(* === Left eigenvectors (row) === *)

ClearAll[VecL,VecL1];
VecL[rep1_,rep2_]:=VecL[rep1,rep2]=Transpose[VecR[rep1,rep2]];
VecL1[rep1_]:=VecL1[rep1]=Transpose[VecR1[rep1]];

ClearAll[EigenL];
EigenL[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=EigenL[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=Table[Quiet[Flatten[VecL[rep1,rep2]/.Solve[
CoefficientList[((VecL[rep1,rep2].H2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-HSpectrum[2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z][[idx]]* VecL[rep1,rep2]))//Together//Numerator,z]==0&&
CoefficientList[((VecL[rep1,rep2].H3[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-HSpectrum[3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z][[idx]]* VecL[rep1,rep2]))//Together//Numerator,z]==0&&
CoefficientList[((VecL[rep1,rep2].H4[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-HSpectrum[4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z][[idx]]* VecL[rep1,rep2]))//Together//Numerator,z]==0,Flatten[VecL[rep1,rep2]]],1]],{idx,1,Dimensions[HSpectrum[1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]][[1]]}]


(* === Example === *)

EigenL[FUNDAMENTAL,FUNDAMENTAL,{2,0,0,0},3,5,7,9]


(* === Normalisation === *)


(* === Highest weight state === *)

(* h.w.s. 2-sites *)

ClearAll[UPR,UPL];
UPR[rep1_,rep2_]:=UPR[rep1,rep2]=Table[{KroneckerDelta[j,1]},{j,1,dim[rep1]*dim[rep2]}];
UPL[rep1_,rep2_]:=UPL[rep1,rep2]=Transpose[UPR[rep1,rep2]];

(* h.w.s. 1-site *)

ClearAll[UPR1,UPL1];
UPR1[rep1_]:=UPR1[rep1]=Table[{KroneckerDelta[j,1]},{j,1,dim[rep1]}];
UPL1[rep1_]:=UPL1[rep1]=Transpose[UPR1[rep1]];


(* === Vacuum vector === *)

(* Vacuum vector 2-sites *)

ClearAll[VacR,VacL];
VacR[rep1_,rep2_,\[Beta]12_,\[Beta]13_,\[Beta]14_,\[Beta]23_,\[Beta]24_,\[Beta]34_]:=VacR[rep1,rep2,\[Beta]12,\[Beta]13,\[Beta]14,\[Beta]23,\[Beta]24,\[Beta]34]=
MatrixExp[\[Beta]24(EMAll[4,2,rep1,rep2])].MatrixExp[\[Beta]23(EMAll[3,2,rep1,rep2])].MatrixExp[\[Beta]14(EMAll[4,1,rep1,rep2])].MatrixExp[\[Beta]13(EMAll[3,1,rep1,rep2])].MatrixExp[\[Beta]34(EMAll[4,3,rep1,rep2])].MatrixExp[\[Beta]12(EMAll[2,1,rep1,rep2])].UPR[rep1,rep2];
VacL[rep1_,rep2_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=VacL[rep1,rep2,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=UPL[rep1,rep2].MatrixExp[\[Alpha]12(EMAll[1,2,rep1,rep2])].MatrixExp[\[Alpha]34(EMAll[3,4,rep1,rep2])].MatrixExp[\[Alpha]13(EMAll[1,3,rep1,rep2])].MatrixExp[\[Alpha]14(EMAll[1,4,rep1,rep2])].MatrixExp[\[Alpha]23(EMAll[2,3,rep1,rep2])].MatrixExp[\[Alpha]24(EMAll[2,4,rep1,rep2])];

(* Vacuum vector 1-site *)

ClearAll[VacR1,VacL1];
VacR1[rep1_,\[Beta]12_,\[Beta]13_,\[Beta]14_,\[Beta]23_,\[Beta]24_,\[Beta]34_]:=VacR1[rep1,\[Beta]23,\[Beta]12,\[Beta]13,\[Beta]14,\[Beta]23,\[Beta]24,\[Beta]34]=
MatrixExp[\[Beta]24(EM[4,2,rep1])].MatrixExp[\[Beta]23(EM[3,2,rep1])].MatrixExp[\[Beta]14(EM[4,1,rep1])].MatrixExp[\[Beta]13(EM[3,1,rep1])].MatrixExp[\[Beta]34(EM[4,3,rep1])].MatrixExp[\[Beta]12(EM[2,1,rep1])].UPR1[rep1];
VacL1[rep1_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=VacL1[rep1,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=UPL1[rep1].MatrixExp[\[Alpha]12(EM[1,2,rep1])].MatrixExp[\[Alpha]34(EM[3,4,rep1])].MatrixExp[\[Alpha]13(EM[1,3,rep1])].MatrixExp[\[Alpha]14(EM[1,4,rep1])].MatrixExp[\[Alpha]23(EM[2,3,rep1])].MatrixExp[\[Alpha]24(EM[2,4,rep1])];


(* === Renormalisation === *)

ClearAll[ReNormR,ReNormL,\[Alpha]23,\[Alpha]13,\[Alpha]12,\[Alpha]24,\[Alpha]14,\[Alpha]34,\[Beta]23,\[Beta]13,\[Beta]12,\[Beta]24,\[Beta]14,\[Beta]34];
ReNormR[rep1_,rep2_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_,\[Psi]_]:=ReNormR[rep1,rep2,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34,\[Psi]]=(\[Psi]*(1/(((VacL[rep1,rep2,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34].\[Psi])//Flatten)[[1]])))//Simplify;
ReNormL[rep1_,rep2_,\[Beta]12_,\[Beta]13_,\[Beta]14_,\[Beta]23_,\[Beta]24_,\[Beta]34_,\[Psi]_]:=ReNormL[rep1,rep2,\[Beta]12,\[Beta]13,\[Beta]14,\[Beta]23,\[Beta]24,\[Beta]34,\[Psi]]=(\[Psi]*(1/((\[Psi].VacR[rep1,rep2,\[Beta]12,\[Beta]13,\[Beta]14,\[Beta]23,\[Beta]24,\[Beta]34]//Flatten)[[1]])))//Simplify;


(* === Example of use === *)

ReNormR[FUNDAMENTAL,FUNDAMENTAL,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34,EigenR[FUNDAMENTAL,FUNDAMENTAL,{1,1,0,0},3,5,7,9][[1]]]

ReNormL[FUNDAMENTAL,FUNDAMENTAL,\[Beta]12,\[Beta]13,\[Beta]14,\[Beta]23,\[Beta]24,\[Beta]34,EigenL[FUNDAMENTAL,FUNDAMENTAL,{1,1,0,0},3,5,7,9][[1]]]

((ReNormL[FUNDAMENTAL,FUNDAMENTAL,\[Beta]12,\[Beta]13,\[Beta]14,\[Beta]23,\[Beta]24,\[Beta]34,EigenL[FUNDAMENTAL,FUNDAMENTAL,{1,1,0,0},3,5,7,9][[1]]].ReNormR[FUNDAMENTAL,FUNDAMENTAL,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34,EigenR[FUNDAMENTAL,FUNDAMENTAL,{1,1,0,0},3,5,7,9][[1]]])//Flatten)[[1]]//FullSimplify


(* === Correct normalisation === *)


(* === \[Beta] without fixing \[Alpha] === *)

ClearAll[\[Beta]rT,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];
\[Beta]rT[12][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=\[Beta]rT[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=((\[Alpha]23 (\[Phi]1-\[Phi]4))/(\[Alpha]13 (\[Phi]2-\[Phi]4)));
\[Beta]rT[13][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=\[Beta]rT[13][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=(((\[Phi]2-\[Phi]3) (\[Phi]1-\[Phi]4))/(\[Alpha]13 (\[Phi]1-\[Phi]2) (\[Phi]3-\[Phi]4)));
\[Beta]rT[14][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=\[Beta]rT[14][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=-(((\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))/(\[Alpha]14 (\[Phi]1-\[Phi]2) (\[Phi]3-\[Phi]4)));
\[Beta]rT[23][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=\[Beta]rT[23][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=-(((\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))/(\[Alpha]23 (\[Phi]1-\[Phi]2) (\[Phi]3-\[Phi]4)));
\[Beta]rT[24][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=\[Beta]rT[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=((\[Alpha]13 (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))/(\[Alpha]14 \[Alpha]23 (\[Phi]1-\[Phi]2) (\[Phi]3-\[Phi]4)));
\[Beta]rT[34][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=\[Beta]rT[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=((\[Alpha]13 (\[Phi]2-\[Phi]4))/(\[Alpha]14 (-\[Phi]2+\[Phi]3)));


(* === fixing \[Alpha] === *)

ClearAll[\[Alpha]r,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];
\[Alpha]r[12][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=((\[Alpha]13 (-\[Phi]1+\[Phi]3))/(\[Alpha]23 (\[Phi]2-\[Phi]3)));
\[Alpha]r[24][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=((\[Alpha]14 \[Alpha]23 (\[Phi]2-\[Phi]3) (\[Phi]1-\[Phi]4))/(\[Alpha]13 (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4)));
\[Alpha]r[34][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=((\[Alpha]14 (\[Phi]1-\[Phi]4))/(\[Alpha]13 (\[Phi]1-\[Phi]3)));


(* === \[Beta] with fixed \[Alpha] === *)

ClearAll[\[Beta]r,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];
\[Beta]r[12][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Beta]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=\[Beta]rT[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]];
\[Beta]r[13][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Beta]r[13][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=\[Beta]rT[13][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]];
\[Beta]r[14][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Beta]r[14][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=\[Beta]rT[14][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]];
\[Beta]r[23][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Beta]r[23][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=\[Beta]rT[23][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]];
\[Beta]r[24][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Beta]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=\[Beta]rT[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]];
\[Beta]r[34][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Beta]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=\[Beta]rT[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]];


(* === Renormalised eigenvectors === *)

ClearAll[ReEigenR,ReEigenL];
ReEigenR[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=ReEigenR[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=ReNormR[rep1,rep2,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],#]&/@EigenR[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];
ReEigenL[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=ReEigenL[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=ReNormL[rep1,rep2,\[Beta]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Beta]r[13][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Beta]r[14][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Beta]r[23][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Beta]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Beta]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],#]&/@EigenL[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];


(*
ReEigenR[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]14_,\[Alpha]34_]:=ReEigenR[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34]=ReNormR[rep1,rep2,\[Alpha]r[23][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Alpha]r[13][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Alpha]12,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Alpha]14,\[Alpha]34,#]&/@EigenR[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];
ReEigenL[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]14_,\[Alpha]34_]:=ReEigenL[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34]=ReNormL[rep1,rep2,\[Beta]r[23][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Beta]r[13][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Beta]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Beta]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Beta]r[14][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Beta]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],#]&/@EigenL[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];
*)
(*
ReEigenR[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]14_,\[Alpha]34_]:=ReEigenR[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34]=ReNormR[rep1,rep2,\[Alpha]r[23][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Alpha]r[13][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Alpha]12,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Alpha]14,\[Alpha]34,#]&/@EigenR[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];
ReEigenL[rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]14_,\[Alpha]34_]:=ReEigenL[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34]=ReNormL[rep1,rep2,\[Beta]r[23][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Beta]r[13][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Beta]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Beta]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Beta]r[14][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],\[Beta]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34],#]&/@EigenL[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];
*)


(* === Example of use === *)

ReEigenL[FUNDAMENTAL,FUNDAMENTAL,{2,0,0,0},\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]14,\[Alpha]34]


(* === RUN: CREATE BASIS === *)


(* === Q Base === *)


(* === Q Base given repAll === *)

ClearAll[QBase,T];
QBase[T_][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=QBase[T][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=
Q[T][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]/.SolveBaxAll[rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1];


(* === Q Hodge Base === *)


(* === Q Hodge Base given repAll === *)

ClearAll[QHBase,T,U];
QHBase[T_][rep1_,rep2_,repAll_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=QHBase[T][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=Table[Hodge[{QBase[1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]],QBase[2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]],QBase[3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]],QBase[4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]]},rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,z,T],{U,1,Length[QBase[1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]]}];



(*Q\[Dagger][T][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]/.SolveBaxAllh\[Dagger][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1,1];*)

(*With Hodge*)
(*Table[Hodge[{QBase[1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]],QBase[2][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]],QBase[3][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]],QBase[4][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]]},rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,z,T],{U,1,Length[QBase[1][rep1,rep2,repAll,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]]}];*)


(* === Q Vacuum === *)

ClearAll[QVac2,QHVac2,QvQv2,QvQv2\[Omega]\[Phi]];
(*QvQv*)
(*VACUUM STATES*)
QVac2[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=QVac2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]={QBase[1][rep1,rep2,rep1+rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]],QBase[2][rep1,rep2,rep1+rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]],QBase[3][rep1,rep2,rep1+rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]],QBase[4][rep1,rep2,rep1+rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]]};
QHVac2[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=QHVac2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]={QHBase[1][rep1,rep2,rep1+rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]],QHBase[2][rep1,rep2,rep1+rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]],QHBase[3][rep1,rep2,rep1+rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]],QHBase[4][rep1,rep2,rep1+rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]]};

QvQv2[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=QvQv2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=OrthogRen[QHVac2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4],QVac2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4],QHVac2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4],QVac2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4],rep1,rep2];

QvQv2\[Omega]\[Phi][rep1_,rep2_,\[Omega]1_,\[Omega]2_,\[Omega]3_,\[Omega]4_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=QvQv2\[Omega]\[Phi][rep1,rep2,\[Omega]1,\[Omega]2,\[Omega]3,\[Omega]4,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=OrthogRen[QHVac2[rep1,rep2,\[Omega]1,\[Omega]2,\[Omega]3,\[Omega]4],QVac2[rep1,rep2,\[Omega]1,\[Omega]2,\[Omega]3,\[Omega]4],QHVac2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4],QVac2[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4],rep1,rep2];


(* === RUN: Q space ONE SITE === *)


(* === hs === *)


(* === Cartan === *)

ClearAll[h1];
h1[rep1_][1,1]:=h1[rep1][1,1]=rep1[[1]]+rep1[[2]]+rep1[[3]]+rep1[[4]];
h1[rep1_][2,1]:=h1[rep1][2,1]=rep1\[LeftDoubleBracket]3\[RightDoubleBracket] (-1+rep1\[LeftDoubleBracket]4\[RightDoubleBracket])+rep1\[LeftDoubleBracket]2\[RightDoubleBracket] (-2+rep1\[LeftDoubleBracket]3\[RightDoubleBracket]+rep1\[LeftDoubleBracket]4\[RightDoubleBracket])+rep1\[LeftDoubleBracket]1\[RightDoubleBracket] (-3+rep1\[LeftDoubleBracket]2\[RightDoubleBracket]+rep1\[LeftDoubleBracket]3\[RightDoubleBracket]+rep1\[LeftDoubleBracket]4\[RightDoubleBracket]);
h1[rep1_][3,1]:=h1[rep1][3,1]=rep1\[LeftDoubleBracket]2\[RightDoubleBracket] (-1+rep1\[LeftDoubleBracket]3\[RightDoubleBracket]) (-2+rep1\[LeftDoubleBracket]4\[RightDoubleBracket])+rep1\[LeftDoubleBracket]1\[RightDoubleBracket] ((-2+rep1\[LeftDoubleBracket]3\[RightDoubleBracket]) (-3+rep1\[LeftDoubleBracket]4\[RightDoubleBracket])+rep1\[LeftDoubleBracket]2\[RightDoubleBracket] (-4+rep1\[LeftDoubleBracket]3\[RightDoubleBracket]+rep1\[LeftDoubleBracket]4\[RightDoubleBracket]));
h1[rep1_][4,1]:=h1[rep1][4,1]=rep1\[LeftDoubleBracket]1\[RightDoubleBracket] (-1+rep1\[LeftDoubleBracket]2\[RightDoubleBracket]) (-2+rep1\[LeftDoubleBracket]3\[RightDoubleBracket]) (-3+rep1\[LeftDoubleBracket]4\[RightDoubleBracket]);


(* === TextData[{Cell[((z)^(0)), ExpressionUUID->   "b137a0c3-3c28-43aa-8700-0da3f629ba0b"],  (dependence only on regulators)}] === *)

(*H1*)
h1[\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][1,3]:=h1[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][1,3]=\[Phi]1+\[Phi]2+\[Phi]3+\[Phi]4;

(*H2*)
h1[\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][2,5]:=h1[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][2,5]=\[Phi]1 \[Phi]2+\[Phi]1 \[Phi]3+\[Phi]2 \[Phi]3+\[Phi]1 \[Phi]4+\[Phi]2 \[Phi]4+\[Phi]3 \[Phi]4;

(*H3*)
h1[\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][3,7]:=h1[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][3,7]=\[Phi]1 \[Phi]2 \[Phi]3+\[Phi]1 \[Phi]2 \[Phi]4+\[Phi]1 \[Phi]3 \[Phi]4+\[Phi]2 \[Phi]3 \[Phi]4;

(*H4*)
h1[\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][4,9]:=h1[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][4,9]=\[Phi]1*\[Phi]2*\[Phi]3*\[Phi]4;


(* === Relation (dependence on weight1) === *)

(* TextData[{By expanding at z=\[Infinity] and taking the coefficient of , Cell[((1)/(z)), ExpressionUUID->   "9880c854-af1e-4c9d-b9ac-ac74245e5014"]}] *)

(*H2*)
h1[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][2,3]:=h1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][2,3]=(\[Phi]2+\[Phi]3+\[Phi]4) (weight1[[1]])+(\[Phi]1+\[Phi]3+\[Phi]4) (weight1[[2]])+(\[Phi]1+\[Phi]2+\[Phi]4) (weight1[[3]])+(\[Phi]1+\[Phi]2+\[Phi]3) (weight1[[4]]);

(*H3*)
h1[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][3,5]:=h1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][3,5]=\[Phi]3 \[Phi]4 (weight1[[1]]+weight1[[2]])+\[Phi]2 \[Phi]4 (weight1[[1]]+weight1[[3]])+\[Phi]1 \[Phi]4 (weight1[[2]]+weight1[[3]])+\[Phi]2 \[Phi]3 (weight1[[1]]+weight1[[4]])+\[Phi]1 \[Phi]3 (weight1[[2]]+weight1[[4]])+\[Phi]1 \[Phi]2 (weight1[[3]]+weight1[[4]]);

(*H4*)
h1[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][4,7]:=h1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][4,7]=\[Phi]2 \[Phi]3 \[Phi]4 (weight1[[1]])+\[Phi]1 \[Phi]3 \[Phi]4 (weight1[[2]])+\[Phi]1 \[Phi]2 \[Phi]4 (weight1[[3]])+\[Phi]1 \[Phi]2 \[Phi]3 (weight1[[4]]);


(* === Hamiltonian in Q-space === *)

ClearAll[\[Epsilon]11,\[Epsilon]12,\[Epsilon]13];
\[Epsilon]11[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]11[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((h1[rep1][1,1])/(z-w1))+h1[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][1,3];
\[Epsilon]12[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]12[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((h1[rep1][2,1])/((((z-w1))^(2))))+((h1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][2,3])/((z-w1)))+h1[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][2,5];
\[Epsilon]13[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]13[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((h1[rep1][3,1])/((((z-w1))^(3))))+((h1[3,3])/((((z-w1))^(2))))+((h1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][3,5])/(z-w1))+h1[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][3,7];
\[Epsilon]14[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]14[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((h1[rep1][4,1])/((((z-w1))^(4))))+((h1[4,3])/((((z-w1))^(3))))+((h1[4,5])/((((z-w1))^(2))))+((h1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][4,7])/(z-w1))+h1[\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][4,9];


(* === Baxter equation with F === *)

ClearAll[F1];
F1[rep1_,C1_][z_]:=F1[rep1,C1][z]=((1)/((((z-w1))^(rep1[[C1]]-C1+1))));

(* Definition of shifted Hamiltonian *)

ClearAll[\[Epsilon]11S,\[Epsilon]12S,\[Epsilon]13S]
\[Epsilon]11S[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_][z_]:=\[Epsilon]11S[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1][z]=\[Epsilon]11[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+((4 ((F1[rep1,C1])^(\[Prime]))[z])/(F1[rep1,C1][z]));
\[Epsilon]12S[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_][z_]:=\[Epsilon]12S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1][z]=\[Epsilon]12[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+((3 \[Epsilon]11[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] ((F1[rep1,C1])^(\[Prime]))[z])/(F1[rep1,C1][z]))+((6 ((F1[rep1,C1])^(\[Prime]\[Prime]))[z])/(F1[rep1,C1][z]));
\[Epsilon]13S[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_][z_]:=\[Epsilon]13S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1][z]=\[Epsilon]13[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+2*\[Epsilon]12[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*((((F1[rep1,C1])^(\[Prime]))[z])/(F1[rep1,C1][z]))+3*((\[Epsilon]11[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] ((F1[rep1,C1])^(\[Prime]\[Prime]))[z])/(F1[rep1,C1][z]))+4*((((F1[rep1,C1])^((3)))[z])/(F1[rep1,C1][z]));
\[Epsilon]14S[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_][z_]:=\[Epsilon]14S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1][z]=\[Epsilon]14[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+\[Epsilon]13[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*((((F1[rep1,C1])^(\[Prime]))[z])/(F1[rep1,C1][z]))+((\[Epsilon]12[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z] ((F1[rep1,C1])^(\[Prime]\[Prime]))[z])/(F1[rep1,C1][z]))+\[Epsilon]11[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*((((F1[rep1,C1])^((3)))[z])/(F1[rep1,C1][z]))+((((F1[rep1,C1])^((4)))[z])/(F1[rep1,C1][z]));

ClearAll[BaxOp1];
BaxOp1[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_]:=D[#,{z,4}]+\[Epsilon]11S[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1][z]*D[#,{z,3}]+\[Epsilon]12S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1][z]*D[#,{z,2}]+\[Epsilon]13S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1][z]*D[#,z]+\[Epsilon]14S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1][z]*#&;


(* === Q functions === *)

ClearAll[QPoly,Q1,c11,c12,c13,c14];
QPoly[c_,M_,r_][z_]:=QPoly[c,M,r][z]=Exp[r z]*(((z)^(M))+Sum[c[k]*((z)^(k)),{k,0,M-1}]);
Q1[1][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q1[1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly[c11,rep1\[LeftDoubleBracket]1\[RightDoubleBracket]-weight1\[LeftDoubleBracket]1\[RightDoubleBracket],-\[Phi]1][z];
Q1[2][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q1[2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly[c12,rep1\[LeftDoubleBracket]1\[RightDoubleBracket]-weight1\[LeftDoubleBracket]2\[RightDoubleBracket],-\[Phi]2][z];
Q1[3][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q1[3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly[c13,rep1\[LeftDoubleBracket]1\[RightDoubleBracket]-weight1\[LeftDoubleBracket]3\[RightDoubleBracket],-\[Phi]3][z];
Q1[4][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q1[4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly[c14,rep1\[LeftDoubleBracket]1\[RightDoubleBracket]-weight1\[LeftDoubleBracket]4\[RightDoubleBracket],-\[Phi]4][z];


(* === Baxter equations === *)

(* Baxter equations *)

ClearAll[Bax1];
Bax1[1][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_]:=Bax1[1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1]=((BaxOp1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1][Q1[1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[\[Phi]1 z])//Simplify;
Bax1[2][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_]:=Bax1[2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1]=((BaxOp1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1][Q1[2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[\[Phi]2 z])//Simplify;
Bax1[3][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_]:=Bax1[3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1]=((BaxOp1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1][Q1[3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[\[Phi]3 z])//Simplify;
Bax1[4][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,C1_]:=Bax1[4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1]=((BaxOp1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,C1][Q1[4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[\[Phi]3 z])//Simplify;

(* Unknown coefficients *)

ClearAll[Vars11,Vars12,Vars13,Vars14];
Vars11[rep1_,weight1_]:=Vars11[rep1,weight1]=Table[c11[n],{n,0,rep1\[LeftDoubleBracket]1\[RightDoubleBracket]-weight1\[LeftDoubleBracket]1\[RightDoubleBracket]-1}];
Vars12[rep1_,weight1_]:=Vars12[rep1,weight1]=Table[c12[n],{n,0,rep1\[LeftDoubleBracket]1\[RightDoubleBracket]-weight1\[LeftDoubleBracket]2\[RightDoubleBracket]-1}];
Vars13[rep1_,weight1_]:=Vars13[rep1,weight1]=Table[c13[n],{n,0,rep1\[LeftDoubleBracket]1\[RightDoubleBracket]-weight1\[LeftDoubleBracket]3\[RightDoubleBracket]-1}];
Vars14[rep1_,weight1_]:=Vars14[rep1,weight1]=Table[c14[n],{n,0,rep1\[LeftDoubleBracket]1\[RightDoubleBracket]-weight1\[LeftDoubleBracket]4\[RightDoubleBracket]-1}];

(* Solve Baxter equation *)

ClearAll[SolveBaxAll1];
SolveBaxAll1[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,c1_]:=SolveBaxAll1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,c1]=Solve[CoefficientList[Bax1[1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,c1]//Together//Numerator,z]==0&&
CoefficientList[Bax1[2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,c1]//Together//Numerator,z]==0&&CoefficientList[Bax1[3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,c1]//Together//Numerator,z]==0&&CoefficientList[Bax1[4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,c1]//Together//Numerator,z]==0,Join[Vars11[rep1,weight1],Vars12[rep1,weight1],Vars13[rep1,weight1],Vars14[rep1,weight1],{h1[3,3],h1[4,3],h1[4,5]}]]//Simplify;

(* Spectrum of Hamiltonian *)

ClearAll[HSpectrum1];
HSpectrum1[1][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=HSpectrum1[1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]11[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]/.SolveBaxAll1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1];

HSpectrum1[2][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=HSpectrum1[2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]12[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]/.SolveBaxAll1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1];

HSpectrum1[3][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=HSpectrum1[3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]13[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]/.SolveBaxAll1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1];

HSpectrum1[4][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=HSpectrum1[4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]14[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]/.SolveBaxAll1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1];

(* Shifted Spectrum of Hamiltonian *)

ClearAll[HSpectrum1S];
HSpectrum1S[1][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=HSpectrum1S[1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]11S[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]/.SolveBaxAll1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1];

HSpectrum1S[2][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=HSpectrum1S[2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]12S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]/.SolveBaxAll1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1];

HSpectrum1S[3][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=HSpectrum1S[3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]13S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]/.SolveBaxAll1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1];

HSpectrum1S[4][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=HSpectrum1S[4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]14S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]/.SolveBaxAll1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1];


(* === Example of use === *)

Sol=SolveBaxAll1[FUNDAMENTAL,FUNDAMENTAL,1,((1)/(3)),((1)/(5)),((1)/(7)),1]

H1ST=HSpectrum1[1][FUNDAMENTAL,FUNDAMENTAL,1,((1)/(3)),((1)/(5)),((1)/(7))][z]
H2ST=HSpectrum1[2][FUNDAMENTAL,FUNDAMENTAL,1,((1)/(3)),((1)/(5)),((1)/(7))][z]
H3ST=HSpectrum1[3][FUNDAMENTAL,FUNDAMENTAL,1,((1)/(3)),((1)/(5)),((1)/(7))][z]
H4ST=HSpectrum1[4][FUNDAMENTAL,FUNDAMENTAL,1,((1)/(3)),((1)/(5)),((1)/(7))][z]


(* === Measure === *)


(* === 1-site === *)

(* TextData[{Measure such that <f , Cell[Subscript[O, F]g> Overscript[=, !], ExpressionUUID->   "9a620987-a867-4819-8e05-0087d06d66c5"],  <g , Cell[((Subscript[O, F])^(\[Dagger])) f>, ExpressionUUID->"05e7e2de-6033-4e8b-8273-63cfde6d96be"]}] *)

ClearAll[\[Mu]1,\[Mu]2];
\[Mu]1[rep1_][z_]:=\[Mu]1[rep1][z]=((z)^(-rep1[[1]]+rep1[[4]]));
\[Mu]2[rep2_][z_]:=\[Mu]2[rep2][z]=(((-1+z))^(-rep2[[1]]+rep2[[4]]));

(* Mean value *)

ClearAll[MeanValue1,MeanValue2];
MeanValue1[f_,\[Tau]_,rep1_]:=MeanValue1[f,\[Tau],rep1]=
2 \[Pi] I Residue[(f*\[Mu]1[rep1][z]),{z,\[Tau]}];
MeanValue2[f_,\[Tau]_,rep2_]:=MeanValue2[f,\[Tau],rep2]=
2 \[Pi] I Residue[(f*\[Mu]2[rep2][z]),{z,\[Tau]}];


(* === TextData[{Adjoint Baxter equation with , Cell[((F)^(\[Dagger])), ExpressionUUID->   "69c4bc11-ec2a-49b2-9957-9c1ff9ff946a"]}] === *)

(* \<\ HERE we consider only the case C1=1,C2=1 in the definition of M\ \> *)

(* TextData[{Shifted Hamiltonians (already with , Cell[((F)^(\[Dagger])), ExpressionUUID->   "55fe020f-a3ff-4fd7-98d4-22a059b4a84a"], )}] *)

ClearAll[\[Epsilon]11S\[Dagger],\[Epsilon]12S\[Dagger],\[Epsilon]13S\[Dagger]];
\[Epsilon]11S\[Dagger][rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]11S\[Dagger][rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=-\[Epsilon]11S[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]+4*((D[\[Mu]1[rep1][z],z])/(\[Mu]1[rep1][z]));
\[Epsilon]12S\[Dagger][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]12S\[Dagger][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]12S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]-3*D[\[Epsilon]11S[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z],z]-3*\[Epsilon]11S[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]*((D[\[Mu]1[rep1][z],z])/(\[Mu]1[rep1][z]))+6*((D[\[Mu]1[rep1][z],{z,2}])/(\[Mu]1[rep1][z]));
\[Epsilon]13S\[Dagger][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]13S\[Dagger][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=-\[Epsilon]13S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]+2*D[\[Epsilon]12S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z],z]-3*D[\[Epsilon]11S[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z],{z,2}]+((D[\[Mu]1[rep1][z],z])/(\[Mu]1[rep1][z]))*(2*\[Epsilon]12S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]-6*D[\[Epsilon]11S[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z],z])-3*((D[\[Mu]1[rep1][z],{z,2}])/(\[Mu]1[rep1][z]))*\[Epsilon]11S[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]+4*((D[\[Mu]1[rep1][z],{z,3}])/(\[Mu]1[rep1][z]));
\[Epsilon]14S\[Dagger][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=\[Epsilon]14S\[Dagger][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=\[Epsilon]14S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]-D[\[Epsilon]13S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z],z]+D[\[Epsilon]12S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z],{z,2}]-D[\[Epsilon]11S[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z],{z,3}]+((D[\[Mu]1[rep1][z],z])/(\[Mu]1[rep1][z]))*(-\[Epsilon]13S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]+2*D[\[Epsilon]12S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z],z]-3*D[\[Epsilon]11S[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z],{z,2}])+((D[\[Mu]1[rep1][z],{z,2}])/(\[Mu]1[rep1][z]))*(\[Epsilon]12S[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]-3D[\[Epsilon]11S[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z],z])-\[Epsilon]11S[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1][z]*((D[\[Mu]1[rep1][z],{z,3}])/(\[Mu]1[rep1][z]))+((D[\[Mu]1[rep1][z],{z,4}])/(\[Mu]1[rep1][z]));

ClearAll[BaxOp1\[Dagger]];
BaxOp1\[Dagger][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=D[#,{z,4}]+\[Epsilon]11S\[Dagger][rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*D[#,{z,3}]+\[Epsilon]12S\[Dagger][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*D[#,{z,2}]+\[Epsilon]13S\[Dagger][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*D[#,z]+\[Epsilon]14S\[Dagger][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]*#&;


(* === Q functions === *)

ClearAll[QPoly\[Dagger],Q1\[Dagger],c11\[Dagger],c12\[Dagger],c13\[Dagger],c14\[Dagger]];
QPoly\[Dagger][c_,M_,r_][z_]:=QPoly\[Dagger][c,M,r][z]=Exp[r z]*(((z)^(M))+Sum[c[k]*((z)^(k)),{k,0,M-1}]);
Q1\[Dagger][1][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q1\[Dagger][1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly\[Dagger][c11\[Dagger],-rep1\[LeftDoubleBracket]4\[RightDoubleBracket]+weight1\[LeftDoubleBracket]1\[RightDoubleBracket],\[Phi]1][z];
Q1\[Dagger][2][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q1\[Dagger][2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly\[Dagger][c12\[Dagger],-rep1\[LeftDoubleBracket]4\[RightDoubleBracket]+weight1\[LeftDoubleBracket]2\[RightDoubleBracket],\[Phi]2][z];
Q1\[Dagger][3][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q1\[Dagger][3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly\[Dagger][c13\[Dagger],-rep1\[LeftDoubleBracket]4\[RightDoubleBracket]+weight1\[LeftDoubleBracket]3\[RightDoubleBracket],\[Phi]3][z];
Q1\[Dagger][4][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z]:=Q1\[Dagger][4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=QPoly\[Dagger][c14\[Dagger],-rep1\[LeftDoubleBracket]4\[RightDoubleBracket]+weight1\[LeftDoubleBracket]4\[RightDoubleBracket],\[Phi]4][z];


(* === Baxter equations === *)

ClearAll[Bax1\[Dagger]];
Bax1\[Dagger][1][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=Bax1\[Dagger][1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=((BaxOp1\[Dagger][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][Q1\[Dagger][1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[-\[Phi]1 z])//Simplify;
Bax1\[Dagger][2][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=Bax1\[Dagger][2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=((BaxOp1\[Dagger][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][Q1\[Dagger][2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[-\[Phi]2 z])//Simplify;
Bax1\[Dagger][3][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=Bax1\[Dagger][3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=((BaxOp1\[Dagger][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][Q1\[Dagger][3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[-\[Phi]3 z])//Simplify;
Bax1\[Dagger][4][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=Bax1\[Dagger][4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=((BaxOp1\[Dagger][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][Q1\[Dagger][4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]])*Exp[-\[Phi]4 z])//Simplify;

(* Unknown coefficients *)

ClearAll[Vars11\[Dagger],Vars12\[Dagger],Vars13\[Dagger],Vars14\[Dagger]];
Vars11\[Dagger][rep1_,weight1_]:=Vars11\[Dagger][rep1,weight1]=Table[c11\[Dagger][n],{n,0,-rep1\[LeftDoubleBracket]4\[RightDoubleBracket]+weight1\[LeftDoubleBracket]1\[RightDoubleBracket]-1}];
Vars12\[Dagger][rep1_,weight1_]:=Vars12\[Dagger][rep1,weight1]=Table[c12\[Dagger][n],{n,0,-rep1\[LeftDoubleBracket]4\[RightDoubleBracket]+weight1\[LeftDoubleBracket]2\[RightDoubleBracket]-1}];
Vars13\[Dagger][rep1_,weight1_]:=Vars13\[Dagger][rep1,weight1]=Table[c13\[Dagger][n],{n,0,-rep1\[LeftDoubleBracket]4\[RightDoubleBracket]+weight1\[LeftDoubleBracket]3\[RightDoubleBracket]-1}];
Vars14\[Dagger][rep1_,weight1_]:=Vars14\[Dagger][rep1,weight1]=Table[c14\[Dagger][n],{n,0,-rep1\[LeftDoubleBracket]4\[RightDoubleBracket]+weight1\[LeftDoubleBracket]4\[RightDoubleBracket]-1}];


(* === Solving Baxter equations giving hs === *)

ClearAll[SolveBaxh1\[Dagger],SolveBaxAllh1\[Dagger]];
SolveBaxh1\[Dagger][1][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,h33_,h43_,h45_]:=SolveBaxh1\[Dagger][1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h33,h43,h45]=Solve[CoefficientList[Bax1\[Dagger][1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]/.{h1[3,3]->h33,
h1[4,3]->h43,h1[4,5]->h45}//Together//Numerator,z]==0,Vars11\[Dagger][rep1,weight1]]//Simplify//Flatten;

SolveBaxh1\[Dagger][2][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,h33_,h43_,h45_]:=SolveBaxh1\[Dagger][2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h33,h43,h45]=Solve[CoefficientList[Bax1\[Dagger][2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]/.{h1[3,3]->h33,
h1[4,3]->h43,h1[4,5]->h45}//Together//Numerator,z]==0,Vars12\[Dagger][rep1,weight1]]//Simplify//Flatten;

SolveBaxh1\[Dagger][3][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,h33_,h43_,h45_]:=SolveBaxh1\[Dagger][3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h33,h43,h45]=Solve[CoefficientList[Bax1\[Dagger][3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]/.{h1[3,3]->h33,
h1[4,3]->h43,h1[4,5]->h45}//Together//Numerator,z]==0,Vars13\[Dagger][rep1,weight1]]//Simplify//Flatten;

SolveBaxh1\[Dagger][4][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,h33_,h43_,h45_]:=SolveBaxh1\[Dagger][4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h33,h43,h45]=Solve[CoefficientList[Bax1\[Dagger][4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]/.{h1[3,3]->h33,
h1[4,3]->h43,h1[4,5]->h45}//Together//Numerator,z]==0,Vars14\[Dagger][rep1,weight1]]//Simplify//Flatten;

SolveBaxAllh1\[Dagger][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,h33_,h43_,h45_]:=SolveBaxAllh1\[Dagger][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h33,h43,h45]=Join[SolveBaxh1\[Dagger][1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h33,h43,h45],SolveBaxh1\[Dagger][2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h33,h43,h45],SolveBaxh1\[Dagger][3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h33,h43,h45],SolveBaxh1\[Dagger][4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h33,h43,h45]]//Simplify;


(* === Example of use === *)

Solh1=(SolveBaxAll1[FUNDAMENTAL,FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1]//Simplify)[[1]]

SolveBaxh1\[Dagger][1][FUNDAMENTAL,FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h1[3,3]/.Solh1,h1[4,3]/.Solh1,h1[4,5]/.Solh1]
SolveBaxh1\[Dagger][2][FUNDAMENTAL,FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h1[3,3]/.Solh1,h1[4,3]/.Solh1,h1[4,5]/.Solh1]
SolveBaxh1\[Dagger][3][FUNDAMENTAL,FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h1[3,3]/.Solh1,h1[4,3]/.Solh1,h1[4,5]/.Solh1]
SolveBaxh1\[Dagger][4][FUNDAMENTAL,FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h1[3,3]/.Solh1,h1[4,3]/.Solh1,h1[4,5]/.Solh1]

Print["\<Hi\>"]
SolveBaxAllh1\[Dagger][FUNDAMENTAL,FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h1[3,3]/.Solh1,h1[4,3]/.Solh1,h1[4,5]/.Solh1]


(* === Q Hodge === *)

ClearAll[Ng,SetN];
Ng=4;
SetN=Range[Ng];

ClearAll[Hodge1];
Hodge1[Qfour_,rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,z_,i_]:=
(Exp[z*(\[Phi]1+\[Phi]2+\[Phi]3+\[Phi]4)]*Signature[Join[Complement[SetN,{i}],{i}]]*Wronskian[{Qfour[[Complement[SetN,{i}]]][[1]],Qfour[[Complement[SetN,{i}]]][[2]],Qfour[[Complement[SetN,{i}]]][[3]]},z]*(((((z))^(-(2*rep1[[1]]-rep1[[2]]-rep1[[3]])))))/((({\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[1]]]]-{\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[2]]]]))*({\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[1]]]]-{\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[3]]]])*({\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[2]]]]-{\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4}[[Complement[SetN,{i}][[3]]]])))//Simplify;
(*The Q argument should be Q1[z], not Q1!*)

(**(((((z))^(-2(rep1[[1]]-rep1[[2]])))*(((z-1))^(-2(rep2[[1]]-rep2[[2]])))))*)
(**)


(* === Renormalised scalar product === *)

(* TextData[{Input: , Cell[Q={Subscript[Q, 1],Subscript[Q, 2],Subscript[Q, 3]}, ExpressionUUID->   "65876c6a-eee5-4ee6-ba76-d24c880daa37"],  , Cell[QH={Subscript[QH, 1],Subscript[QH, 2],Subscript[QH, 3]}, ExpressionUUID->   "b781f233-2228-4500-9d63-dc6c7d646ea7"]}] *)


(* === rep1 === *)

ClearAll[Normalisation1];
Normalisation1[QHL_,QL_,QHR_,QR_,rep1_]:=Normalisation1[QHL,QL,QHR,QR,rep1]=(Limit[QHL[[1]],z->0,Direction->"\<FromAbove\>"])*(Limit[Wronskian[{QHL[[1]],QHL[[2]]},z]*((1)/(((z)^(rep1[[3]]-rep1[[4]])))),z->0,Direction->"\<FromAbove\>"])*(Limit[QL[[4]],z->0,Direction->"\<FromAbove\>"])*(Limit[QR[[1]],z->0,Direction->"\<FromAbove\>"])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z]*((1)/(((z)^(rep1[[1]]-rep1[[2]])))),z->0,Direction->"\<FromAbove\>"])*(Limit[QHR[[4]],z->0,Direction->"\<FromAbove\>"]);

ClearAll[RawOrthog1];
RawOrthog1[QH_,Q_,rep1_,w_,\[Alpha]_,i_,j_]:=RawOrthog1[QH,Q,rep1,w,\[Alpha],i,j]={MeanValue1[QH[[i]]*D[Q[[j]],{z,2}]*((1)/(z-w)),\[Alpha],rep1],MeanValue1[QH[[i]]*D[Q[[j]],{z,1}]*((1)/((((z-w))^(2)))),\[Alpha],rep1],MeanValue1[QH[[i]]*D[Q[[j]],{z,1}]*((1)/(z-w)),\[Alpha],rep1],MeanValue1[QH[[i]]*Q[[j]]*((1)/((((z-w))^(3)))),\[Alpha],rep1],MeanValue1[QH[[i]]*Q[[j]]*((1)/((((z-w))^(2)))),\[Alpha],rep1],MeanValue1[QH[[i]]*Q[[j]]*((1)/(z-w)),\[Alpha],rep1]};
ClearAll[OrthogRen1];
OrthogRen1[QHL_,QL_,QHR_,QR_,rep1_]:=OrthogRen1[QHL,QL,QHR,QR,rep1]=(1/Normalisation1[QHL,QL,QHR,QR,rep1])*Det[{RawOrthog1[QHL,QR,rep1,0,0,1,1],RawOrthog1[QHL,QR,rep1,0,0,1,2],RawOrthog1[QHL,QR,rep1,0,0,2,1],RawOrthog1[QHL,QR,rep1,0,0,2,2],RawOrthog1[QHL,QR,rep1,0,0,1,3],RawOrthog1[QHL,QR,rep1,0,0,3,1]}];


(* === rep2 === *)

ClearAll[Normalisation2];
Normalisation2[QHL_,QL_,QHR_,QR_,rep2_]:=Normalisation2[QHL,QL,QHR,QR,rep2]=
(Limit[QHL[[1]],z->1,Direction->"\<FromAbove\>"])*(Limit[Wronskian[{QHL[[1]],QHL[[2]]},z]*((1)/((((z-1))^(rep2[[3]]-rep2[[4]])))),z->1,Direction->"\<FromAbove\>"])*(Limit[QL[[4]],z->1,Direction->"\<FromAbove\>"])*(Limit[QR[[1]],z->1,Direction->"\<FromAbove\>"])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z]*((1)/((((z-1))^(rep2[[1]]-rep2[[2]])))),z->1,Direction->"\<FromAbove\>"])*(Limit[QHR[[4]],z->1,Direction->"\<FromAbove\>"]);

ClearAll[RawOrthog2];
RawOrthog2[QH_,Q_,rep2_,w_,\[Alpha]_,i_,j_]:=RawOrthog2[QH,Q,rep2,w,\[Alpha],i,j]={MeanValue2[QH[[i]]*D[Q[[j]],{z,2}]*((1)/(z-w)),\[Alpha],rep2],MeanValue2[QH[[i]]*D[Q[[j]],{z,1}]*((1)/((((z-w))^(2)))),\[Alpha],rep2],MeanValue2[QH[[i]]*D[Q[[j]],{z,1}]*((1)/(z-w)),\[Alpha],rep2],MeanValue2[QH[[i]]*Q[[j]]*((1)/((((z-w))^(3)))),\[Alpha],rep2],MeanValue2[QH[[i]]*Q[[j]]*((1)/((((z-w))^(2)))),\[Alpha],rep2],MeanValue2[QH[[i]]*Q[[j]]*((1)/(z-w)),\[Alpha],rep2]};
ClearAll[OrthogRen2];
OrthogRen2[QHL_,QL_,QHR_,QR_,rep2_]:=OrthogRen2[QHL,QL,QHR,QR,rep2]=(1/Normalisation2[QHL,QL,QHR,QR,rep2])*Det[{RawOrthog2[QHL,QR,rep2,0,0,1,1],RawOrthog2[QHL,QR,rep2,0,0,1,2],RawOrthog2[QHL,QR,rep2,0,0,2,1],RawOrthog2[QHL,QR,rep2,0,0,2,2],RawOrthog2[QHL,QR,rep2,0,0,1,3],RawOrthog2[QHL,QR,rep2,0,0,3,1]}];


(* === RUN: \[CapitalPsi] space ONE SITE === *)


(* === Lax connection === *)

ClearAll[LIndex1];
LIndex1[a_,b_,rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=LIndex1[a,b,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((1)/(z-w1))*EM[a,b,rep1]+KroneckerDelta[a,1]*KroneckerDelta[b,1]*\[Phi]1*Id[rep1]+KroneckerDelta[a,2]*KroneckerDelta[b,2]*\[Phi]2*Id[rep1]+KroneckerDelta[a,3]*KroneckerDelta[b,3]*\[Phi]3*Id[rep1]+KroneckerDelta[a,4]*KroneckerDelta[b,4]*\[Phi]4*Id[rep1];


(* === Example of use === *)

ClearAll[rep1T];
rep1T=FUNDAMENTAL;
LIndex[1,1,rep1T,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal


(* === H1 === *)


(* === TextData[{Regulated Hamiltonian , Cell[tr(L), ExpressionUUID->   "df38a54a-22fa-47df-a00d-bcefa85f1cdb"]}] === *)

ClearAll[H11];
H11[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=H11[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=Sum[LIndex1[a,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{a,1,Ng}]//Simplify;


(* === Definition Residue of H === *)

Clear[ResH11];
ResH11[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_]:=ResH11[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w]=SeriesCoefficient[H11[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,-1}]


(* === Example of use === *)

H11[FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal;


(* === H2 === *)


(* === TextData[{Regulated Hamiltonian , Cell[((3)/(2))*tr(L')+trA(((L)^(2))), ExpressionUUID->   "52d33b0c-9727-4065-b62e-cbe4d9675f17"]}] === *)


(* === Antisymm. def === *)

ClearAll[trA1L2];
trA1L2[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trA1L2[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((1)/(2))*(Sum[LIndex1[a,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[b,b,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[b,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[a,b,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{a,1,Ng},{b,1,Ng}]);


(* === Def. H2 === *)

ClearAll[H12];
H12[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=H12[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=(trA1L2[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+((3)/(2))*Sum[D[LIndex1[a,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z],{a,1,Ng}])//Simplify


(* === Definition Residue of H === *)

Clear[ResH12];
ResH12[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_]:=ResH12[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w]=SeriesCoefficient[H12[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,-1}]

Clear[ResH12O];
ResH12O[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_,o_]:=ResH12O[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w,o]=SeriesCoefficient[H12[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,o}]


(* === Example of use === *)

H12[FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal;


(* === H3 === *)


(* === H3 Unsimplified === *)


(* === Antisymm. def === *)

ClearAll[trA1L3];
trA1L3[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trA1L3[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((1)/(6))*Sum[LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{i1,1,Ng},{i2,1,Ng},{i3,1,Ng}];

ClearAll[trA1LpL,trA1LLp];
trA1LLp[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trA1LLp[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((1)/(2))*(Sum[LIndex1[a,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex1[b,b,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]-LIndex1[b,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex1[a,b,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z],{a,1,Ng},{b,1,Ng}]);
trA1LpL[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trA1LpL[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((1)/(2))*(Sum[D[LIndex1[a,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex1[b,b,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-D[LIndex1[b,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex1[a,b,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{a,1,Ng},{b,1,Ng}]);


(* === Def H3 === *)

ClearAll[H13U];
H13U[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=H3U[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=trA1L3[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+((4)/(3))*trA1LLp[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+((2)/(3))*trA1LpL[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+Sum[D[LIndex1[a,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{z,2}],{a,1,Ng}]


(* === H3 Simplified === *)

ClearAll[H13];
H13[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=H13[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=H13U[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]


(* === Definition Residue of H === *)

Clear[ResH13];
ResH13[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_]:=ResH13[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w]=SeriesCoefficient[H13[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,-1}]

Clear[ResH13O];
ResH13O[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_,o_]:=ResH13O[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w,o]=SeriesCoefficient[H13[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,o}]


(* === Example of use === *)

H13[FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal//Simplify


(* === H4 === *)


(* === H4 Unsimplified === *)


(* === Antisymm. def === *)

ClearAll[trA1L4];
trA1L4[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trA1L4[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((1)/(4!))*Sum[LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i4,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i4,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i3,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i4,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i4,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i2,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i4,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i3,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i4,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i2,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i4,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i4,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i3,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i4,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i4,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i1,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i4,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i3,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i4,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i1,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i4,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i4,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i2,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i4,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i4,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i1,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i4,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i2,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i4,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i1,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i4,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i3,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i4,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i2,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i4,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i3,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i4,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i1,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i4,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i2,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i4,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i1,i4,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{i1,1,Ng},{i2,1,Ng},{i3,1,Ng},{i4,1,Ng}];

(*((L)^(3))*)
ClearAll[trA1LLLp,trA1LLpL,trA1LpLL];
trA1LLLp[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trA1LLLp[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((1)/(6))*Sum[LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]-LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]-LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]+LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]+LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]-LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z],{i1,1,Ng},{i2,1,Ng},{i3,1,Ng}];

trA1LLpL[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trA1LLpL[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((1)/(6))*Sum[LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{i1,1,Ng},{i2,1,Ng},{i3,1,Ng}];

trA1LpLL[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trA1LpLL[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((1)/(6))*Sum[D[LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-D[LIndex1[i1,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-D[LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i3,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+D[LIndex1[i2,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex1[i3,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+D[LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex1[i1,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i2,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-D[LIndex1[i3,i1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].LIndex1[i2,i2,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].LIndex1[i1,i3,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{i1,1,Ng},{i2,1,Ng},{i3,1,Ng}];

(*((L)^(2))*)
ClearAll[trA1LLpp,trA1LpLp,trA1LppL];
trA1LLpp[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trA1LLpp[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((1)/(2))*(Sum[LIndex1[a,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex1[b,b,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{z,2}]-LIndex1[b,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].D[LIndex1[a,b,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{z,2}],{a,1,Ng},{b,1,Ng}]);

trA1LpLp[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trA1LpLp[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((1)/(2))*(Sum[D[LIndex1[a,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].D[LIndex1[b,b,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z]-D[LIndex1[b,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z].D[LIndex1[a,b,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],z],{a,1,Ng},{b,1,Ng}]);

trA1LppL[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=trA1LppL[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=((1)/(2))*(Sum[D[LIndex1[a,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{z,2}].LIndex1[b,b,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-D[LIndex1[b,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{z,2}].LIndex1[a,b,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{a,1,Ng},{b,1,Ng}]);


(* === Def H4 === *)

ClearAll[H14U];
H14U[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=H14U[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=trA1L4[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+((3)/(4))*trA1LLLp[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+((1)/(2))*trA1LLpL[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+((1)/(4))*trA1LpLL[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+
((1)/(2))*trA1LLpp[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+((1)/(2))*trA1LpLp[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+((1)/(6))*trA1LppL[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]+
((1)/(4))*Sum[D[LIndex1[a,a,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z],{z,3}],{a,1,Ng}]


(* === H4 Simplified === *)

ClearAll[H14];
H14[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_][z_]:=H14[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]=H14U[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]


(* === Definition Residue of H === *)

Clear[ResH14];
ResH14[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_]:=ResH14[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w]=SeriesCoefficient[H14[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,-1}]

Clear[ResH14O];
ResH14O[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,w_,o_]:=ResH14O[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,w,o]=SeriesCoefficient[H14[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal,{z,w,o}]


(* === Example of use === *)

H14[FUNDAMENTAL,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]//Normal//Simplify;


(* === H1\[Dagger] === *)


(* === TextData[{Regulated Hamiltonian , Cell[tr(L), ExpressionUUID->   "770a6f65-79c6-42b1-9632-6373dc9f7d1b"]}] === *)

ClearAll[H1\[Dagger]];
H1\[Dagger][rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_][z_]:=H1\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]=H1[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z](*Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}]//Simplify;*)


(* === Definition Residue of H === *)

Clear[ResH1\[Dagger]];
ResH1\[Dagger][rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,w_]:=ResH1\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,w]=SeriesCoefficient[H1\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]//Normal,{z,w,-1}]


(* === Example of use === *)

ClearAll[rep1,rep2];
rep1=FUNDAMENTAL;
rep2=FUNDAMENTAL;
H1\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]//Normal;


(* === H2\[Dagger] === *)


(* === TextData[{Regulated Hamiltonian -, Cell[tr(L')+((1)/(2))(((tr(L)))^(2))-((1)/(2))tr(((L)^(2))), ExpressionUUID->   "c407857c-548c-4cb2-8d4b-e721db1c1733"]}] === *)

ClearAll[H2\[Dagger]];
H2\[Dagger][rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_][z_]:=H2\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]=(-Sum[D[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],z],{a,1,3}]+((1)/(2))*(Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}]).(Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}])-((1)/(2))*Sum[LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z].LIndex[b,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3},{b,1,3}])//Simplify


(* === Definition Residue of H === *)

Clear[ResH2\[Dagger]];
ResH2\[Dagger][rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,w_]:=ResH2\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,w]=SeriesCoefficient[H2\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]//Normal,{z,w,-1}]


(* === H3\[Dagger] === *)


(* === TextData[{H3\[Dagger] simplified , Cell[((1)/(3))tr(L'')+tr(L'.L)-((2)/(3))tr(L').tr(L)-((1)/(2))tr(L).tr(((L)^(2)))+((1)/(6))((tr(L))^(3))+((1)/(3))tr(((L)^(3))), ExpressionUUID->   "e54d8888-4526-4975-ad4d-10addffadc7c"]}] === *)

ClearAll[H3\[Dagger]];
H3\[Dagger][rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_][z_]:=H3\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]=(((1)/(3))*Sum[D[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{z,2}],{a,1,3}]+Sum[D[LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],z].LIndex[b,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3},{b,1,3}]-((2)/(3))*Sum[D[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],z],{a,1,3}].Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}]-((1)/(2))*Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}].Sum[LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z].LIndex[b,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3},{b,1,3}]+((1)/(6))*(Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}]).(Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}]).(Sum[LIndex[a,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3}])+((1)/(3))*Sum[LIndex[a,b,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z].LIndex[b,c,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z].LIndex[c,a,rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z],{a,1,3},{b,1,3},{c,1,3}])//Simplify;


(* === Definition Residue of H === *)

Clear[ResH3\[Dagger]];
ResH3\[Dagger][rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,w_]:=ResH3\[Dagger][rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,w]=SeriesCoefficient[H3[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]//Normal,{z,w,-1}]


(* === Example of use === *)

Series[H3\[Dagger][ADJOINT,rep2,\[Phi]1,\[Phi]2,\[Phi]3][z]//Normal,{z,0,-3}]//Simplify//Normal;


(* === Right eigenvectors (column) === *)

ClearAll[VecR1];
VecR1[rep1_]:=VecR1[rep1]=Table[{v[j]},{j,1,dim[rep1]}];

ClearAll[EigenR1];
EigenR1[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=EigenR1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=Table[Quiet[Flatten[VecR1[rep1]/.Solve[
CoefficientList[((H12[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].VecR1[rep1]-HSpectrum1[2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z][[idx]]* VecR1[rep1]))//Together//Numerator,z]==0&&
CoefficientList[((H13[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].VecR1[rep1]-HSpectrum1[3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z][[idx]]* VecR1[rep1]))//Together//Numerator,z]==0&&
CoefficientList[((H14[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z].VecR1[rep1]-HSpectrum1[4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z][[idx]]* VecR1[rep1]))//Together//Numerator,z]==0,Flatten[VecR1[rep1]]],1]],{idx,1,Dimensions[HSpectrum1[1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]][[1]]}]


(* === Example === *)

EigenR1[FUNDAMENTAL,FUNDAMENTAL,3,5,7,9]


(* === Left eigenvectors (row) === *)

ClearAll[VecL1];
VecL1[rep1_]:=VecL1[rep1]=Transpose[VecR1[rep1]];

ClearAll[EigenL1];
EigenL1[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=EigenL1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=Table[Quiet[Flatten[VecL1[rep1]/.Solve[
CoefficientList[((VecL1[rep1].H12[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-HSpectrum1[2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z][[idx]]* VecL1[rep1]))//Together//Numerator,z]==0&&
CoefficientList[((VecL1[rep1].H13[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-HSpectrum1[3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z][[idx]]* VecL1[rep1]))//Together//Numerator,z]==0&&
CoefficientList[((VecL1[rep1].H14[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]-HSpectrum1[4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z][[idx]]* VecL1[rep1]))//Together//Numerator,z]==0,Flatten[VecL1[rep1]]],1]],{idx,1,Dimensions[HSpectrum1[1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]][[1]]}]


(* === Example === *)

EigenL1[FUNDAMENTAL,FUNDAMENTAL,3,5,7,9]


(* === Normalisation === *)


(* === Highest weight state === *)

(* h.w.s. 1-site *)

ClearAll[UPR1,UPL1];
UPR1[rep1_]:=UPR1[rep1]=Table[{KroneckerDelta[j,1]},{j,1,dim[rep1]}];
UPL1[rep1_]:=UPL1[rep1]=Transpose[UPR1[rep1]];


(* === Vacuum vector === *)

(* Vacuum vector 1-site *)

ClearAll[VacR1,VacL1];
VacR1[rep1_,\[Beta]12_,\[Beta]13_,\[Beta]14_,\[Beta]23_,\[Beta]24_,\[Beta]34_]:=VacR1[rep1,\[Beta]23,\[Beta]12,\[Beta]13,\[Beta]14,\[Beta]23,\[Beta]24,\[Beta]34]=
MatrixExp[\[Beta]24(EM[4,2,rep1])].MatrixExp[\[Beta]23(EM[3,2,rep1])].MatrixExp[\[Beta]14(EM[4,1,rep1])].MatrixExp[\[Beta]13(EM[3,1,rep1])].MatrixExp[\[Beta]34(EM[4,3,rep1])].MatrixExp[\[Beta]12(EM[2,1,rep1])].UPR1[rep1];
VacL1[rep1_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=VacL1[rep1,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=UPL1[rep1].MatrixExp[\[Alpha]12(EM[1,2,rep1])].MatrixExp[\[Alpha]34(EM[3,4,rep1])].MatrixExp[\[Alpha]13(EM[1,3,rep1])].MatrixExp[\[Alpha]14(EM[1,4,rep1])].MatrixExp[\[Alpha]23(EM[2,3,rep1])].MatrixExp[\[Alpha]24(EM[2,4,rep1])];


(* === Renormalisation === *)

ClearAll[ReNormR1,ReNormL1,\[Alpha]23,\[Alpha]13,\[Alpha]12,\[Alpha]24,\[Alpha]14,\[Alpha]34,\[Beta]23,\[Beta]13,\[Beta]12,\[Beta]24,\[Beta]14,\[Beta]34];
ReNormR1[rep1_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_,\[Psi]_]:=ReNormR1[rep1,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34,\[Psi]]=(\[Psi]*(1/(((VacL1[rep1,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34].\[Psi])//Flatten)[[1]])))//Simplify;
ReNormL1[rep1_,\[Beta]12_,\[Beta]13_,\[Beta]14_,\[Beta]23_,\[Beta]24_,\[Beta]34_,\[Psi]_]:=ReNormL1[rep1,\[Beta]12,\[Beta]13,\[Beta]14,\[Beta]23,\[Beta]24,\[Beta]34,\[Psi]]=(\[Psi]*(1/((\[Psi].VacR1[rep1,\[Beta]12,\[Beta]13,\[Beta]14,\[Beta]23,\[Beta]24,\[Beta]34]//Flatten)[[1]])))//Simplify;


(* === Example of use === *)

ReNormR1[FUNDAMENTAL,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34,EigenR1[FUNDAMENTAL,FUNDAMENTAL,3,5,7,9][[1]]]

ReNormL1[FUNDAMENTAL,\[Beta]12,\[Beta]13,\[Beta]14,\[Beta]23,\[Beta]24,\[Beta]34,EigenL1[FUNDAMENTAL,FUNDAMENTAL,3,5,7,9][[1]]]

((ReNormL1[FUNDAMENTAL,\[Beta]12,\[Beta]13,\[Beta]14,\[Beta]23,\[Beta]24,\[Beta]34,EigenL1[FUNDAMENTAL,FUNDAMENTAL,3,5,7,9][[1]]].ReNormR1[FUNDAMENTAL,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34,EigenR1[FUNDAMENTAL,FUNDAMENTAL,3,5,7,9][[1]]])//Flatten)[[1]]//FullSimplify


(* === Correct normalisation === *)


(* === \[Beta] without fixing \[Alpha] === *)

ClearAll[\[Beta]rT,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];
\[Beta]rT[12][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=\[Beta]rT[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=((\[Alpha]23 (\[Phi]1-\[Phi]4))/(\[Alpha]13 (\[Phi]2-\[Phi]4)));
\[Beta]rT[13][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=\[Beta]rT[13][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=(((\[Phi]2-\[Phi]3) (\[Phi]1-\[Phi]4))/(\[Alpha]13 (\[Phi]1-\[Phi]2) (\[Phi]3-\[Phi]4)));
\[Beta]rT[14][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=\[Beta]rT[14][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=-(((\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))/(\[Alpha]14 (\[Phi]1-\[Phi]2) (\[Phi]3-\[Phi]4)));
\[Beta]rT[23][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=\[Beta]rT[23][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=-(((\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))/(\[Alpha]23 (\[Phi]1-\[Phi]2) (\[Phi]3-\[Phi]4)));
\[Beta]rT[24][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=\[Beta]rT[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=((\[Alpha]13 (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))/(\[Alpha]14 \[Alpha]23 (\[Phi]1-\[Phi]2) (\[Phi]3-\[Phi]4)));
\[Beta]rT[34][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]12_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_,\[Alpha]24_,\[Alpha]34_]:=\[Beta]rT[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34]=((\[Alpha]13 (\[Phi]2-\[Phi]4))/(\[Alpha]14 (-\[Phi]2+\[Phi]3)));


(* === fixing \[Alpha] === *)

ClearAll[\[Alpha]r,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];
\[Alpha]r[12][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=((\[Alpha]13 (-\[Phi]1+\[Phi]3))/(\[Alpha]23 (\[Phi]2-\[Phi]3)));
\[Alpha]r[24][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=((\[Alpha]14 \[Alpha]23 (\[Phi]2-\[Phi]3) (\[Phi]1-\[Phi]4))/(\[Alpha]13 (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4)));
\[Alpha]r[34][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=((\[Alpha]14 (\[Phi]1-\[Phi]4))/(\[Alpha]13 (\[Phi]1-\[Phi]3)));


(* === \[Beta] with fixed \[Alpha] === *)

ClearAll[\[Beta]r,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];
\[Beta]r[12][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Beta]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=\[Beta]rT[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]];
\[Beta]r[13][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Beta]r[13][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=\[Beta]rT[13][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]];
\[Beta]r[14][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Beta]r[14][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=\[Beta]rT[14][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]];
\[Beta]r[23][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Beta]r[23][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=\[Beta]rT[23][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]];
\[Beta]r[24][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Beta]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=\[Beta]rT[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]];
\[Beta]r[34][\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=\[Beta]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=\[Beta]rT[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]];


(* === Renormalised eigenvectors === *)

ClearAll[ReEigenR1,ReEigenL1];
ReEigenR1[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=ReEigenR1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=ReNormR1[rep1,\[Alpha]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],#]&/@EigenR1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];
ReEigenL1[rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]13_,\[Alpha]14_,\[Alpha]23_]:=ReEigenL1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23]=ReNormL1[rep1,\[Beta]r[12][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Beta]r[13][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Beta]r[14][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Beta]r[23][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Beta]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Beta]r[34][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]13,\[Alpha]14,\[Alpha]23],#]&/@EigenL1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4];


(* === RUN: CREATE BASIS ONE SITE === *)


(* === Q Base === *)


(* === Q Base given repAll === *)

ClearAll[QBase1,T];
QBase1[T_][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=QBase1[T][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=
Q1[T][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]/.SolveBaxAll1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1];


(* === Q Hodge Base === *)


(* === Q Hodge Base given repAll === *)

ClearAll[QHBase1,T,U];
QHBase1[T_][rep1_,weight1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=QHBase1[T][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=Table[Hodge1[{QBase1[1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]],QBase1[2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]],QBase1[3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]],QBase1[4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]]},rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,z,T],{U,1,Length[QBase1[1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]]}];

(*Q1\[Dagger][T][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][z]/.SolveBaxAllh1\[Dagger][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,h1[3,3]/.SolveBaxAll1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1],h1[4,3]/.SolveBaxAll1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1],h1[4,5]/.SolveBaxAll1[rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,1]];*)

(*With Hodge*)
(*Table[Hodge[{QBase1[1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]],QBase1[2][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]],QBase1[3][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]],QBase1[4][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[U]]},rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,z,T],{U,1,Length[QBase1[1][rep1,weight1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]]}];*)


(* === Q Vacuum === *)

ClearAll[QVac1,QHVac1,QvQv1,QvQv1\[Omega]\[Phi]];
(*QvQv*)
(*VACUUM STATES*)
QVac1[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=QVac1[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]={QBase1[1][rep1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]],QBase1[2][rep1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]],QBase1[3][rep1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]],QBase1[4][rep1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]]};
QHVac1[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=QHVac1[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]={QHBase1[1][rep1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]],QHBase1[2][rep1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]],QHBase1[3][rep1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]],QHBase1[4][rep1,rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4][[1]]};

QvQv1[rep1_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=QvQv1[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=OrthogRen1[QHVac1[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4],QVac1[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4],QHVac1[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4],QVac1[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4],rep1];

QvQv1\[Omega]\[Phi][rep1_,\[Omega]1_,\[Omega]2_,\[Omega]3_,\[Omega]4_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:=QvQv1\[Omega]\[Phi][rep1,\[Omega]1,\[Omega]2,\[Omega]3,\[Omega]4,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4]=OrthogRen1[QHVac1[rep1,\[Omega]1,\[Omega]2,\[Omega]3,\[Omega]4],QVac1[rep1,\[Omega]1,\[Omega]2,\[Omega]3,\[Omega]4],QHVac1[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4],QVac1[rep1,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4],rep1];


(* === K vector 2-sites EXPLICIT === *)

(* \<\ The first part contains the definition for general {p,p,0,0}x{q,q,0,0} reps. The second part contains explicit definition of K vector for \ {1,1,0,0}x{1,1,0,0} rep.\ \> *)


(* === K vector === *)


(* === (SU(2) 1-site)\[CenterDot](SU(3) 2-sites) === *)

ClearAll[K22B];
K22B[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]\[Phi]13_,\[Alpha]\[Phi]14_,\[Alpha]\[Phi]23_,k1_,k2_,k3_,k4_]:=K22B[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4]=KroneckerProduct[UPL1[rep1].MatrixExp[((\[Alpha]\[Phi]13 (k3-\[Phi]3) (-k4+\[Phi]2+\[Phi]3))/((k3-\[Phi]1) (-k4+\[Phi]1+\[Phi]2)))(EM[1,3,rep1])].MatrixExp[((\[Alpha]\[Phi]14 (k3-\[Phi]4) (-k4+\[Phi]2+\[Phi]4))/((k3-\[Phi]1) (-k4+\[Phi]1+\[Phi]2)))(EM[1,4,rep1])].MatrixExp[((\[Alpha]\[Phi]23 (k3-\[Phi]3) (-k4+\[Phi]1+\[Phi]3))/((k3-\[Phi]2) (-k4+\[Phi]1+\[Phi]2)))(EM[2,3,rep1])].MatrixExp[-((\[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (\[Phi]2-\[Phi]3) (k3-\[Phi]4) (-\[Phi]1+\[Phi]4) (-k4+\[Phi]1+\[Phi]4))/(\[Alpha]\[Phi]13 (k3-\[Phi]2) (-k4+\[Phi]1+\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4)))(EM[2,4,rep1])],UPL1[rep2].MatrixExp[((\[Alpha]\[Phi]13*(k3-\[Phi]3))/((k3-\[Phi]1)))(EM[1,3,rep2])].MatrixExp[((\[Alpha]\[Phi]14*(k3-\[Phi]4))/((k3-\[Phi]1)))(EM[1,4,rep2])].MatrixExp[\[Alpha]\[Phi]23*(((k3-\[Phi]3))/((k3-\[Phi]2)))(EM[2,3,rep2])].MatrixExp[\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23]*(((k3-\[Phi]4))/((k3-\[Phi]2)))(EM[2,4,rep2])]];


(* === (SU(2) 2-sites) \[CenterDot] (SU(3) 1-site) === *)

ClearAll[K22C];
K22C[rep1_,rep2_,\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_,\[Alpha]\[Phi]13_,\[Alpha]\[Phi]14_,\[Alpha]\[Phi]23_,k1_,k2_,k3_,k4_]:=K22C[rep1,rep2,\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4]=KroneckerProduct[UPL1[rep1].MatrixExp[\[Alpha]\[Phi]13*(((k4-\[Phi]2-\[Phi]3)*(k3-\[Phi]3))/((k4-\[Phi]1-\[Phi]2)*(k3-\[Phi]1)))(EM[1,3,rep1])].MatrixExp[\[Alpha]\[Phi]14*(((k3-\[Phi]4)*(k4-\[Phi]2-\[Phi]4))/((k3-\[Phi]1)*(k4-\[Phi]1-\[Phi]2)))(EM[1,4,rep2])].MatrixExp[\[Alpha]\[Phi]23*(((k3-\[Phi]3)*(k4-\[Phi]3-\[Phi]1))/((k3-\[Phi]2)*(k4-\[Phi]1-\[Phi]2)))(EM[2,3,rep2])].MatrixExp[\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23]*(((k4-\[Phi]1-\[Phi]4)*(k3-\[Phi]4))/((k4-\[Phi]1-\[Phi]2)*(k3-\[Phi]2)))(EM[2,4,rep2])],UPL1[rep2].MatrixExp[\[Alpha]\[Phi]13*(((k4-\[Phi]2-\[Phi]3))/((k4-\[Phi]1-\[Phi]2)))(EM[1,3,rep2])].MatrixExp[\[Alpha]\[Phi]14*(((k4-\[Phi]2-\[Phi]4))/((k4-\[Phi]1-\[Phi]2)))(EM[1,4,rep2])].MatrixExp[\[Alpha]\[Phi]23*((k4-\[Phi]3-\[Phi]1)/(k4-\[Phi]1-\[Phi]2))(EM[2,3,rep2])].MatrixExp[\[Alpha]r[24][\[Phi]1,\[Phi]2,\[Phi]3,\[Phi]4,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23]*(((k4-\[Phi]1-\[Phi]4))/((k4-\[Phi]1-\[Phi]2)))(EM[2,4,rep2])]];


(* === Q scalar product === *)


(* === (SU(2) 1-site)\[CenterDot](SU(3) 2-sites) === *)

ClearAll[OrthogRen22B];
OrthogRen22B[QHL_,QL_,QHR_,QR_,rep1_,rep2_]:=OrthogRen22B[QHL,QL,QHR,QR,rep1,rep2]=
2 \[Pi] I Residue[(QHL[[4]]*Wronskian[{QR[[1]],QR[[2]]},z]*(((1)/((((z-w1))^((rep1[[1]]-rep1[[4]])+1 )) )))),{z,0}]/((Limit[QHL[[4]],z->w1])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z],z->w1]))*(1/((Limit[QHL[[3]],z->0,Direction->"\<FromAbove\>"])*(Limit[D[QHL[[3]],z],z->0,Direction->"\<FromAbove\>"])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z],z->0,Direction->"\<FromAbove\>"])*
(Limit[QHL[[3]],z->1,Direction->"\<FromAbove\>"])*(Limit[D[QHL[[3]],z],z->1,Direction->"\<FromAbove\>"])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z],z->1,Direction->"\<FromAbove\>"])))*
Det[{{2 \[Pi] I Residue[((QHL[[3]] D[QR[[1]],z])/(z))*(((1)/((((z-w1))^(rep1[[1]]-rep1[[3]])) (((z-w2))^(rep2[[1]]-rep2[[3]]))))),{z,0}],
2 \[Pi] I Residue[((QHL[[3]] QR[[1]])/(((z)^(2))))*(((1)/((((z-w1))^(rep1[[1]]-rep1[[3]])) (((z-w2))^(rep2[[1]]-rep2[[3]]))))),{z,0}],
2 \[Pi] I Residue[((QHL[[3]] D[QR[[1]],z])/(-1+z))*(((1)/((((z-w1))^(rep1[[1]]-rep1[[3]])) (((z-w2))^(rep2[[1]]-rep2[[3]]))))),{z,0}],
2 \[Pi] I Residue[((QHL[[3]] QR[[1]])/((((-1+z))^(2))))*(((1)/((((z-w1))^(rep1[[1]]-rep1[[3]])) (((z-w2))^(rep2[[1]]-rep2[[3]]))))),{z,0}]},

{2 \[Pi] I Residue[((QHL[[3]] D[QR[[2]],z])/(z))*(((1)/((((z-w1))^(rep1[[1]]-rep1[[3]])) (((z-w2))^(rep2[[1]]-rep2[[3]]))))),{z,0}],
2 \[Pi] I Residue[((QHL[[3]] QR[[2]])/(((z)^(2))))*(((1)/((((z-w1))^(rep1[[1]]-rep1[[3]])) (((z-w2))^(rep2[[1]]-rep2[[3]]))))),{z,0}],
2 \[Pi] I Residue[((QHL[[3]] D[QR[[2]],z])/(-1+z))*(((1)/((((z-w1))^(rep1[[1]]-rep1[[3]])) (((z-w2))^(rep2[[1]]-rep2[[3]]))))),{z,0}],
2 \[Pi] I Residue[((QHL[[3]] QR[[2]])/((((-1+z))^(2))))*(((1)/((((z-w1))^(rep1[[1]]-rep1[[3]])) (((z-w2))^(rep2[[1]]-rep2[[3]]))))),{z,0}]},


{2 \[Pi] I Residue[((QHL[[3]] D[QR[[1]],z])/(z))*(((1)/((((z-w1))^(rep1[[1]]-rep1[[3]])) (((z-w2))^(rep2[[1]]-rep2[[3]]))))),{z,1}],
2 \[Pi] I Residue[((QHL[[3]] QR[[1]])/(((z)^(2))))*(((1)/((((z-w1))^(rep1[[1]]-rep1[[3]])) (((z-w2))^(rep2[[1]]-rep2[[3]]))))),{z,1}],
2 \[Pi] I Residue[((QHL[[3]] D[QR[[1]],z])/(-1+z))*(((1)/((((z-w1))^(rep1[[1]]-rep1[[3]])) (((z-w2))^(rep2[[1]]-rep2[[3]]))))),{z,1}],
2 \[Pi] I Residue[((QHL[[3]] QR[[1]])/((((-1+z))^(2))))*(((1)/((((z-w1))^(rep1[[1]]-rep1[[3]])) (((z-w2))^(rep2[[1]]-rep2[[3]]))))),{z,1}]},

{2 \[Pi] I Residue[((QHL[[3]] D[QR[[2]],z])/(z))*(((1)/((((z-w1))^(rep1[[1]]-rep1[[3]])) (((z-w2))^(rep2[[1]]-rep2[[3]]))))),{z,1}],
2 \[Pi] I Residue[((QHL[[3]] QR[[2]])/(((z)^(2))))*(((1)/((((z-w1))^(rep1[[1]]-rep1[[3]])) (((z-w2))^(rep2[[1]]-rep2[[3]]))))),{z,1}],
2 \[Pi] I Residue[((QHL[[3]] D[QR[[2]],z])/(-1+z))*(((1)/((((z-w1))^(rep1[[1]]-rep1[[3]])) (((z-w2))^(rep2[[1]]-rep2[[3]]))))),{z,1}],
2 \[Pi] I Residue[((QHL[[3]] QR[[2]])/((((-1+z))^(2))))*(((1)/((((z-w1))^(rep1[[1]]-rep1[[3]])) (((z-w2))^(rep2[[1]]-rep2[[3]]))))),{z,1}]}}];


(* === (SU(2) 2-sites) \[CenterDot] (SU(3) 1-site) === *)

ClearAll[OrthogRen22C];
OrthogRen22C[QHL_,QL_,QHR_,QR_,rep1_,rep2_]:=OrthogRen22C[QHL,QL,QHR,QR,rep1,rep2]=(2 \[Pi] I Residue[(QHL[[4]]*Wronskian[{QR[[1]],QR[[2]]},z]*(((1)/((((z-w1))^((rep1[[1]]-rep1[[4]]) +1)) (((z-w2))^((rep2[[1]]-rep2[[4]]) +1)))))),{z,w1}]*2 \[Pi] I Residue[(QHL[[4]]*Wronskian[{QR[[1]],QR[[2]]},z]*z*(((1)/((((z-w1))^((rep1[[1]]-rep1[[4]]) +1)) (((z-w2))^((rep2[[1]]-rep2[[4]]) +1)))))),{z,w2}]-2 \[Pi] I Residue[(QHL[[4]]*Wronskian[{QR[[1]],QR[[2]]},z]*(((1)/((((z-w1))^((rep1[[1]]-rep1[[4]]) +1)) (((z-w2))^((rep2[[1]]-rep2[[4]]) +1)))))),{z,w2}]*2 \[Pi] I Residue[(QHL[[4]]*Wronskian[{QR[[1]],QR[[2]]},z]*z*(((1)/((((z-w1))^((rep1[[1]]-rep1[[4]]) +1)) (((z-w2))^((rep2[[1]]-rep2[[4]]) +1)))))),{z,w1}])/((Limit[QHL[[4]],z->w1])*(Limit[QHL[[4]],z->w2])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z],z->w1])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z],z->w2]))*
Det[{{2 \[Pi] I Residue[QHL[[3]]*D[QR[[1]],z]*((1)/(z))*(((1)/((((z-w1))^(rep1[[1]]-rep1[[3]]))))),{z,0}],2 \[Pi] I Residue[QHL[[3]]*QR[[1]]*((1)/(((z)^(2))))*(((1)/((((z-w1))^(rep1[[1]]-rep1[[3]]))))),{z,0}]},

{2 \[Pi] I Residue[QHL[[3]]*D[QR[[2]],z]*((1)/(z))*(((1)/((((z-w1))^(rep1[[1]]-rep1[[3]]))))),{z,0}],2 \[Pi] I Residue[QHL[[3]]*QR[[2]]*((1)/(((z)^(2))))*(((1)/((((z-w1))^(rep1[[1]]-rep1[[3]]))))),{z,0}]}}]/((Limit[QHL[[3]],z->0,Direction->"\<FromAbove\>"])*(Limit[D[QHL[[3]],z],z->0,Direction->"\<FromAbove\>"])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z],z->0,Direction->"\<FromAbove\>"]));


(* === Modules for testing === *)


(* === RUN: Modules === *)

ClearAll[Q2Base,QH2Base,Q1Base,QH1Base];
ClearAll[rep1T,rep2T,repAllT,Weight,QT\[Phi],QHT\[Phi],QT1\[Phi],QHT1\[Phi],M,W,MaxM,MaxW,KT,T,REig\[Phi],LEig\[Phi]];
ClearAll[\[Phi]1T,\[Phi]1,\[Phi]2T,\[Phi]2,\[Phi]3T,\[Phi]3,\[Phi]4T,\[Phi]4,k1,k2,k3,k4];

(*Create Q base*)
Q2Base[rep1T_,rep2T_,repAllT_,Weight_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_][M_]:=Q2Base[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][M]=Table[{QBase[1][rep1T,rep2T,repAllT[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],QBase[2][rep1T,rep2T,repAllT[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],
QBase[3][rep1T,rep2T,repAllT[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],QBase[4][rep1T,rep2T,repAllT[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]]},{T,1,Length[QBase[1][rep1T,rep2T,repAllT[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T]]}];

QH2Base[rep1T_,rep2T_,repAllT_,Weight_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_][M_]:=QH2Base[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][M]=Table[{QHBase[1][rep1T,rep2T,repAllT[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],QHBase[2][rep1T,rep2T,repAllT[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],
QHBase[3][rep1T,rep2T,repAllT[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],QHBase[4][rep1T,rep2T,repAllT[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]]},{T,1,Length[QHBase[1][rep1T,rep2T,repAllT[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T]]}];

Q1Base[rep1T_,rep2T_,repAllT_,Weight_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_][W_]:=Q1Base[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][W]=Table[{QBase1[1][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],QBase1[2][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],
QBase1[3][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],QBase1[4][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]]},{T,1,Length[QBase1[1][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T]]}];

QH1Base[rep1T_,rep2T_,repAllT_,Weight_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_][W_]:=QH1Base[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][W]=Table[{QHBase1[1][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],QHBase1[2][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],
QHBase1[3][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]],QHBase1[4][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][[T]]},{T,1,Length[QHBase1[1][rep1T,Weight[[W]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T]]}];

(*Create \[CapitalPsi] base*)
RBase[rep1T_,rep2T_,repAllT_,Weight_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_,\[Alpha]\[Phi]13_,\[Alpha]\[Phi]14_,\[Alpha]\[Phi]23_][M_]:=RBase[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23][M]=ReEigenR[rep1T,rep2T,repAllT[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23]//Simplify;

LBase[rep1T_,rep2T_,repAllT_,Weight_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_,\[Alpha]\[Phi]13_,\[Alpha]\[Phi]14_,\[Alpha]\[Phi]23_][M_]:=LBase[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23][M]=ReEigenL[rep1T,rep2T,repAllT[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23]//Simplify;

(*Create \[CapitalPsi] base*)
RBase1P[rep1T_,rep2T_,repAllT_,Weight_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_,\[Alpha]\[Phi]13_,\[Alpha]\[Phi]14_,\[Alpha]\[Phi]23_][M_]:=RBase1P[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23][M]=ReEigenR1[rep1T,Weight[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23]//Simplify;

LBase1P[rep1T_,rep2T_,repAllT_,Weight_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_,\[Alpha]\[Phi]13_,\[Alpha]\[Phi]14_,\[Alpha]\[Phi]23_][M_]:=LBase1P[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23][M]=ReEigenL1[rep1T,Weight[[M]],\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23]//Simplify;

(* \<\ Type of rep: {1,0,0,0}: 1 {1,1,0,0}: 2 {1,1,1,0}:3\ \> *)

ClearAll[TESTK22B,TESTK22C];
(*SU(2) 1-site \[CenterDot] SU(3) 2-sites*)
TESTK22B[rep1T_,rep2T_,repAllT_,Weight_,MaxM_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_,\[Alpha]\[Phi]13_,\[Alpha]\[Phi]14_,\[Alpha]\[Phi]23_,k1_,k2_,k3_,k4_][L_]:=TESTK22B[rep1T,rep2T,repAllT,Weight,MaxM,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4][L]=
Table[
((K22B[rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4].RBase[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23][L][[U]]//Flatten)[[1]])-\[Mu]R*((OrthogRen22B[QHVac1[rep1T,k1,k2,k3,k4],QVac1[rep1T,k1,k2,k3,k4],QH2Base[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][L][[U]],Q2Base[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][L][[U]],rep1T,rep2T])/(OrthogRen22B[QHVac1[rep1T,k1,k2,k3,k4],QVac1[rep1T,k1,k2,k3,k4],QHVac2[rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T],QVac2[rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T],rep1T,rep2T])),{U,1,Dimensions[RBase[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23][L]][[1]]}];

(*SU(2) 2-site \[CenterDot] SU(3) 1-sites*)
TESTK22C[rep1T_,rep2T_,repAllT_,Weight_,MaxM_,\[Phi]1T_,\[Phi]2T_,\[Phi]3T_,\[Phi]4T_,\[Alpha]\[Phi]13_,\[Alpha]\[Phi]14_,\[Alpha]\[Phi]23_,k1_,k2_,k3_,k4_][L_]:=TESTK22C[rep1T,rep2T,repAllT,Weight,MaxM,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4][L]=
Table[
((K22C[rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23,k1,k2,k3,k4].RBase[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23][L][[U]]//Flatten)[[1]])-\[Mu]R*((OrthogRen22C[QHVac1[rep1T,k1,k2,k3,k4],QVac1[rep1T,k1,k2,k3,k4],QH2Base[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][L][[U]],Q2Base[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][L][[U]],rep1T,rep2T])/(OrthogRen22C[QHVac1[rep1T,k1,k2,k3,k4],QVac1[rep1T,k1,k2,k3,k4],QHVac2[rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T],QVac2[rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T],rep1T,rep2T])),{U,1,Dimensions[RBase[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]\[Phi]13,\[Alpha]\[Phi]14,\[Alpha]\[Phi]23][L]][[1]]}];


(* === Test reconstruction === *)


(* === (SU(2) 1-site) \[CenterDot] (SU(3) 2-sites) === *)

rep1T={1,1,0,0};
rep2T={1,1,0,0};
(*FindGoodWeight2[rep1T,rep2T]*)
repAllT=FindWeight2[rep1T,rep2T]
(*Permutations[rep1T+rep2T]*)
MaxM=Length[repAllT];
Weight=FindGoodWeight[rep1T];
MaxW=Length[FindGoodWeight[rep1T]];
\[Phi]1T=1;
\[Phi]2T=((1)/(2));
\[Phi]3T=((1)/(3));
\[Phi]4T=((1)/(5));

ClearAll[vecT,M,Mult,SolCartNorm,vecT2,SolF,vecT3,ExpectedRes];
vecT=VecR[rep1T,rep2T];
M=3;
Mult=1;
Print["\<Rep\>"];
repAllT[[M]]
Print["\<Imposing Cartan eigenvalues\>"];
VacL[rep1T,rep2T,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34];
SolCartNorm=Solve[EMAll[1,1,rep1T,rep2T].vecT==repAllT[[M]][[1]]*vecT&&EMAll[2,2,rep1T,rep2T].vecT==repAllT[[M]][[2]]*vecT&&EMAll[3,3,rep1T,rep2T].vecT==repAllT[[M]][[3]]*vecT&&EMAll[4,4,rep1T,rep2T].vecT==repAllT[[M]][[4]]*vecT&&(((VacL[rep1T,rep2T,\[Alpha]r[12][\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]13,\[Alpha]14,\[Alpha]23]].vecT)//Flatten)[[1]])==1,Flatten[vecT]]
Print["\<Substituting solution\>"];
vecT2=vecT/.SolCartNorm[[1]]

ClearAll[k1,k2,k3,k4];
Print["\<Solution from expansion\>"];
SolF=Solve[(CoefficientList[((((K22B[rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]13,\[Alpha]14,\[Alpha]23,k1,k2,k3,k4].vecT2)//Flatten)[[1]])//Simplify)-((OrthogRen22B[QHVac1[rep1T,k1,k2,k3,k4],QVac1[rep1T,k1,k2,k3,k4],QH2Base[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][M][[Mult]],Q2Base[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][M][[Mult]],rep1T,rep2T])/(OrthogRen22B[QHVac1[rep1T,k1,k2,k3,k4],QVac1[rep1T,k1,k2,k3,k4],QHVac2[rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T],QVac2[rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T],rep1T,rep2T]))//Together//Numerator,{k2,k3,k4}]//Flatten)==0,Flatten[vecT]]
vecT3=vecT2/.SolF[[1]]//Simplify
ExpectedRes=RBase[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]13,\[Alpha]14,\[Alpha]23][M][[Mult]]
vecT3==ExpectedRes//FullSimplify


(* === (SU(2) 2-sites) \[CenterDot] (SU(3) 1-site) === *)

rep1T={1,1,0,0};
rep2T={1,1,0,0};
(*FindGoodWeight2[rep1T,rep2T]*)
repAllT=FindWeight2[rep1T,rep2T]
(*Permutations[rep1T+rep2T]*)
MaxM=Length[repAllT];
Weight=FindGoodWeight[rep1T];
MaxW=Length[FindGoodWeight[rep1T]];
\[Phi]1T=1;
\[Phi]2T=((1)/(2));
\[Phi]3T=((1)/(3));
\[Phi]4T=((1)/(5));

ClearAll[vecT,M,Mult,SolCartNorm,vecT2,SolF,vecT3,ExpectedRes];
vecT=VecR[rep1T,rep2T];
M=3;
Mult=1;
Print["\<Rep\>"];
repAllT[[M]]
Print["\<Imposing Cartan eigenvalues\>"];
VacL[rep1T,rep2T,\[Alpha]12,\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]24,\[Alpha]34];
SolCartNorm=Solve[EMAll[1,1,rep1T,rep2T].vecT==repAllT[[M]][[1]]*vecT&&EMAll[2,2,rep1T,rep2T].vecT==repAllT[[M]][[2]]*vecT&&EMAll[3,3,rep1T,rep2T].vecT==repAllT[[M]][[3]]*vecT&&EMAll[4,4,rep1T,rep2T].vecT==repAllT[[M]][[4]]*vecT&&(((VacL[rep1T,rep2T,\[Alpha]r[12][\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]13,\[Alpha]14,\[Alpha]23,\[Alpha]r[24][\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]13,\[Alpha]14,\[Alpha]23],\[Alpha]r[34][\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]13,\[Alpha]14,\[Alpha]23]].vecT)//Flatten)[[1]])==1,Flatten[vecT]]
Print["\<Substituting solution\>"];
vecT2=vecT/.SolCartNorm[[1]]

ClearAll[k1,k2,k3,k4];
K22B[rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]13,\[Alpha]14,\[Alpha]23,k1,k2,k3,k4].vecT2

SolF=Solve[(CoefficientList[((((K22C[rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]13,\[Alpha]14,\[Alpha]23,k1,k2,k3,k4].vecT2)//Flatten)[[1]])//Simplify)-((OrthogRen22C[QHVac1[rep1T,k1,k2,k3,k4],QVac1[rep1T,k1,k2,k3,k4],QH2Base[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][M][[Mult]],Q2Base[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T][M][[Mult]],rep1T,rep2T])/(OrthogRen22C[QHVac1[rep1T,k1,k2,k3,k4],QVac1[rep1T,k1,k2,k3,k4],QHVac2[rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T],QVac2[rep1T,rep2T,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T],rep1T,rep2T]))//Together//Numerator,{k2,k3,k4}]//Flatten)==0,Flatten[vecT]]
vecT3=vecT2/.SolF[[1]]//Simplify
ExpectedRes=RBase[rep1T,rep2T,repAllT,Weight,\[Phi]1T,\[Phi]2T,\[Phi]3T,\[Phi]4T,\[Alpha]13,\[Alpha]14,\[Alpha]23][M][[Mult]]
vecT3==ExpectedRes//FullSimplify


(* === Explicitly for {1,1,0,0}x{1,1,0,0} === *)


(* === K vectors === *)


(* === (SU(2) 1-site)\[CenterDot](SU(3) 2-sites) === *)

KBVec={{1,((\[Alpha]\[Phi]23 (k3-\[Phi]3))/(k3-\[Phi]2)),-((2 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (\[Phi]2-\[Phi]3) (k3-\[Phi]4) (-\[Phi]1+\[Phi]4))/(\[Alpha]\[Phi]13 (k3-\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))),-((\[Alpha]\[Phi]13 (k3-\[Phi]3))/(k3-\[Phi]1)),-((2 \[Alpha]\[Phi]14 (k3-\[Phi]4))/(k3-\[Phi]1)),-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (k3-\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2)))-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (\[Phi]2-\[Phi]3) (k3-\[Phi]4) (-\[Phi]1+\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))),((\[Alpha]\[Phi]23 (k3-\[Phi]3) (-k4+\[Phi]1+\[Phi]3))/((k3-\[Phi]2) (-k4+\[Phi]1+\[Phi]2))),((((\[Alpha]\[Phi]23)^(2)) (((k3-\[Phi]3))^(2)) (-k4+\[Phi]1+\[Phi]3))/((((k3-\[Phi]2))^(2)) (-k4+\[Phi]1+\[Phi]2))),-((2 \[Alpha]\[Phi]14 ((\[Alpha]\[Phi]23)^(2)) (k3-\[Phi]3) (\[Phi]2-\[Phi]3) (-k4+\[Phi]1+\[Phi]3) (k3-\[Phi]4) (-\[Phi]1+\[Phi]4))/(\[Alpha]\[Phi]13 (((k3-\[Phi]2))^(2)) (-k4+\[Phi]1+\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))),-((\[Alpha]\[Phi]13 \[Alpha]\[Phi]23 (((k3-\[Phi]3))^(2)) (-k4+\[Phi]1+\[Phi]3))/((k3-\[Phi]1) (k3-\[Phi]2) (-k4+\[Phi]1+\[Phi]2))),-((2 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (-k4+\[Phi]1+\[Phi]3) (k3-\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (-k4+\[Phi]1+\[Phi]2))),((\[Alpha]\[Phi]23 (k3-\[Phi]3) (-k4+\[Phi]1+\[Phi]3) (-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (k3-\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2)))-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (\[Phi]2-\[Phi]3) (k3-\[Phi]4) (-\[Phi]1+\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4)))))/((k3-\[Phi]2) (-k4+\[Phi]1+\[Phi]2))),-((2 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (\[Phi]2-\[Phi]3) (k3-\[Phi]4) (-\[Phi]1+\[Phi]4) (-k4+\[Phi]1+\[Phi]4))/(\[Alpha]\[Phi]13 (k3-\[Phi]2) (-k4+\[Phi]1+\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))),-((2 \[Alpha]\[Phi]14 ((\[Alpha]\[Phi]23)^(2)) (k3-\[Phi]3) (\[Phi]2-\[Phi]3) (k3-\[Phi]4) (-\[Phi]1+\[Phi]4) (-k4+\[Phi]1+\[Phi]4))/(\[Alpha]\[Phi]13 (((k3-\[Phi]2))^(2)) (-k4+\[Phi]1+\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))),((4 ((\[Alpha]\[Phi]14)^(2)) ((\[Alpha]\[Phi]23)^(2)) (((\[Phi]2-\[Phi]3))^(2)) (((k3-\[Phi]4))^(2)) (((-\[Phi]1+\[Phi]4))^(2)) (-k4+\[Phi]1+\[Phi]4))/(((\[Alpha]\[Phi]13)^(2)) (((k3-\[Phi]2))^(2)) (-k4+\[Phi]1+\[Phi]2) (((\[Phi]1-\[Phi]3))^(2)) (((\[Phi]2-\[Phi]4))^(2)))),((2 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (\[Phi]2-\[Phi]3) (k3-\[Phi]4) (-\[Phi]1+\[Phi]4) (-k4+\[Phi]1+\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (-k4+\[Phi]1+\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))),((4 ((\[Alpha]\[Phi]14)^(2)) \[Alpha]\[Phi]23 (\[Phi]2-\[Phi]3) (((k3-\[Phi]4))^(2)) (-\[Phi]1+\[Phi]4) (-k4+\[Phi]1+\[Phi]4))/(\[Alpha]\[Phi]13 (k3-\[Phi]1) (k3-\[Phi]2) (-k4+\[Phi]1+\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))),-((2 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (\[Phi]2-\[Phi]3) (k3-\[Phi]4) (-\[Phi]1+\[Phi]4) (-k4+\[Phi]1+\[Phi]4) (-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (k3-\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2)))-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (\[Phi]2-\[Phi]3) (k3-\[Phi]4) (-\[Phi]1+\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4)))))/(\[Alpha]\[Phi]13 (k3-\[Phi]2) (-k4+\[Phi]1+\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))),-((\[Alpha]\[Phi]13 (k3-\[Phi]3) (-k4+\[Phi]2+\[Phi]3))/((k3-\[Phi]1) (-k4+\[Phi]1+\[Phi]2))),-((\[Alpha]\[Phi]13 \[Alpha]\[Phi]23 (((k3-\[Phi]3))^(2)) (-k4+\[Phi]2+\[Phi]3))/((k3-\[Phi]1) (k3-\[Phi]2) (-k4+\[Phi]1+\[Phi]2))),((2 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (\[Phi]2-\[Phi]3) (-k4+\[Phi]2+\[Phi]3) (k3-\[Phi]4) (-\[Phi]1+\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (-k4+\[Phi]1+\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))),((((\[Alpha]\[Phi]13)^(2)) (((k3-\[Phi]3))^(2)) (-k4+\[Phi]2+\[Phi]3))/((((k3-\[Phi]1))^(2)) (-k4+\[Phi]1+\[Phi]2))),((2 \[Alpha]\[Phi]13 \[Alpha]\[Phi]14 (k3-\[Phi]3) (-k4+\[Phi]2+\[Phi]3) (k3-\[Phi]4))/((((k3-\[Phi]1))^(2)) (-k4+\[Phi]1+\[Phi]2))),-((\[Alpha]\[Phi]13 (k3-\[Phi]3) (-k4+\[Phi]2+\[Phi]3) (-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (k3-\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2)))-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (\[Phi]2-\[Phi]3) (k3-\[Phi]4) (-\[Phi]1+\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4)))))/((k3-\[Phi]1) (-k4+\[Phi]1+\[Phi]2))),-((2 \[Alpha]\[Phi]14 (k3-\[Phi]4) (-k4+\[Phi]2+\[Phi]4))/((k3-\[Phi]1) (-k4+\[Phi]1+\[Phi]2))),-((2 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (k3-\[Phi]4) (-k4+\[Phi]2+\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (-k4+\[Phi]1+\[Phi]2))),((4 ((\[Alpha]\[Phi]14)^(2)) \[Alpha]\[Phi]23 (\[Phi]2-\[Phi]3) (((k3-\[Phi]4))^(2)) (-\[Phi]1+\[Phi]4) (-k4+\[Phi]2+\[Phi]4))/(\[Alpha]\[Phi]13 (k3-\[Phi]1) (k3-\[Phi]2) (-k4+\[Phi]1+\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))),((2 \[Alpha]\[Phi]13 \[Alpha]\[Phi]14 (k3-\[Phi]3) (k3-\[Phi]4) (-k4+\[Phi]2+\[Phi]4))/((((k3-\[Phi]1))^(2)) (-k4+\[Phi]1+\[Phi]2))),((4 ((\[Alpha]\[Phi]14)^(2)) (((k3-\[Phi]4))^(2)) (-k4+\[Phi]2+\[Phi]4))/((((k3-\[Phi]1))^(2)) (-k4+\[Phi]1+\[Phi]2))),-((2 \[Alpha]\[Phi]14 (k3-\[Phi]4) (-k4+\[Phi]2+\[Phi]4) (-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (k3-\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2)))-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (\[Phi]2-\[Phi]3) (k3-\[Phi]4) (-\[Phi]1+\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4)))))/((k3-\[Phi]1) (-k4+\[Phi]1+\[Phi]2))),-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (\[Phi]2-\[Phi]3) (-k4+\[Phi]2+\[Phi]3) (k3-\[Phi]4) (-\[Phi]1+\[Phi]4) (-k4+\[Phi]1+\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (((-k4+\[Phi]1+\[Phi]2))^(2)) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4)))-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (-k4+\[Phi]1+\[Phi]3) (k3-\[Phi]4) (-k4+\[Phi]2+\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (((-k4+\[Phi]1+\[Phi]2))^(2)))),((\[Alpha]\[Phi]23 (k3-\[Phi]3) (-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (\[Phi]2-\[Phi]3) (-k4+\[Phi]2+\[Phi]3) (k3-\[Phi]4) (-\[Phi]1+\[Phi]4) (-k4+\[Phi]1+\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (((-k4+\[Phi]1+\[Phi]2))^(2)) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4)))-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (-k4+\[Phi]1+\[Phi]3) (k3-\[Phi]4) (-k4+\[Phi]2+\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (((-k4+\[Phi]1+\[Phi]2))^(2))))))/(k3-\[Phi]2)),-((2 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (\[Phi]2-\[Phi]3) (k3-\[Phi]4) (-\[Phi]1+\[Phi]4) (-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (\[Phi]2-\[Phi]3) (-k4+\[Phi]2+\[Phi]3) (k3-\[Phi]4) (-\[Phi]1+\[Phi]4) (-k4+\[Phi]1+\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (((-k4+\[Phi]1+\[Phi]2))^(2)) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4)))-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (-k4+\[Phi]1+\[Phi]3) (k3-\[Phi]4) (-k4+\[Phi]2+\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (((-k4+\[Phi]1+\[Phi]2))^(2))))))/(\[Alpha]\[Phi]13 (k3-\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))),-((\[Alpha]\[Phi]13 (k3-\[Phi]3) (-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (\[Phi]2-\[Phi]3) (-k4+\[Phi]2+\[Phi]3) (k3-\[Phi]4) (-\[Phi]1+\[Phi]4) (-k4+\[Phi]1+\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (((-k4+\[Phi]1+\[Phi]2))^(2)) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4)))-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (-k4+\[Phi]1+\[Phi]3) (k3-\[Phi]4) (-k4+\[Phi]2+\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (((-k4+\[Phi]1+\[Phi]2))^(2))))))/(k3-\[Phi]1)),-((2 \[Alpha]\[Phi]14 (k3-\[Phi]4) (-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (\[Phi]2-\[Phi]3) (-k4+\[Phi]2+\[Phi]3) (k3-\[Phi]4) (-\[Phi]1+\[Phi]4) (-k4+\[Phi]1+\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (((-k4+\[Phi]1+\[Phi]2))^(2)) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4)))-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (-k4+\[Phi]1+\[Phi]3) (k3-\[Phi]4) (-k4+\[Phi]2+\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (((-k4+\[Phi]1+\[Phi]2))^(2))))))/(k3-\[Phi]1)),(-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (k3-\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2)))-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (\[Phi]2-\[Phi]3) (k3-\[Phi]4) (-\[Phi]1+\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4)))) (-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (\[Phi]2-\[Phi]3) (-k4+\[Phi]2+\[Phi]3) (k3-\[Phi]4) (-\[Phi]1+\[Phi]4) (-k4+\[Phi]1+\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (((-k4+\[Phi]1+\[Phi]2))^(2)) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4)))-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (-k4+\[Phi]1+\[Phi]3) (k3-\[Phi]4) (-k4+\[Phi]2+\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (((-k4+\[Phi]1+\[Phi]2))^(2)))))}};


(* === (SU(2) 2-sites) \[CenterDot] (SU(3) 1-site) === *)

KCVec={{1,((\[Alpha]\[Phi]23 (k4-\[Phi]1-\[Phi]3))/(k4-\[Phi]1-\[Phi]2)),((2 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (\[Phi]2-\[Phi]3) (k4-\[Phi]1-\[Phi]4) (\[Phi]1-\[Phi]4))/(\[Alpha]\[Phi]13 (k4-\[Phi]1-\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))),-((\[Alpha]\[Phi]13 (k4-\[Phi]2-\[Phi]3))/(k4-\[Phi]1-\[Phi]2)),-((2 \[Alpha]\[Phi]14 (k4-\[Phi]2-\[Phi]4))/(k4-\[Phi]1-\[Phi]2)),-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k4-\[Phi]1-\[Phi]3) (k4-\[Phi]2-\[Phi]4))/((((k4-\[Phi]1-\[Phi]2))^(2))))+((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k4-\[Phi]2-\[Phi]3) (\[Phi]2-\[Phi]3) (k4-\[Phi]1-\[Phi]4) (\[Phi]1-\[Phi]4))/((((k4-\[Phi]1-\[Phi]2))^(2)) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))),((\[Alpha]\[Phi]23 (k3-\[Phi]3) (-k4+\[Phi]1+\[Phi]3))/((k3-\[Phi]2) (-k4+\[Phi]1+\[Phi]2))),((((\[Alpha]\[Phi]23)^(2)) (k3-\[Phi]3) (k4-\[Phi]1-\[Phi]3) (-k4+\[Phi]1+\[Phi]3))/((k3-\[Phi]2) (k4-\[Phi]1-\[Phi]2) (-k4+\[Phi]1+\[Phi]2))),((2 \[Alpha]\[Phi]14 ((\[Alpha]\[Phi]23)^(2)) (k3-\[Phi]3) (\[Phi]2-\[Phi]3) (-k4+\[Phi]1+\[Phi]3) (k4-\[Phi]1-\[Phi]4) (\[Phi]1-\[Phi]4))/(\[Alpha]\[Phi]13 (k3-\[Phi]2) (k4-\[Phi]1-\[Phi]2) (-k4+\[Phi]1+\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))),-((\[Alpha]\[Phi]13 \[Alpha]\[Phi]23 (k3-\[Phi]3) (k4-\[Phi]2-\[Phi]3) (-k4+\[Phi]1+\[Phi]3))/((k3-\[Phi]2) (k4-\[Phi]1-\[Phi]2) (-k4+\[Phi]1+\[Phi]2))),-((2 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (-k4+\[Phi]1+\[Phi]3) (k4-\[Phi]2-\[Phi]4))/((k3-\[Phi]2) (k4-\[Phi]1-\[Phi]2) (-k4+\[Phi]1+\[Phi]2))),((\[Alpha]\[Phi]23 (k3-\[Phi]3) (-k4+\[Phi]1+\[Phi]3) (-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k4-\[Phi]1-\[Phi]3) (k4-\[Phi]2-\[Phi]4))/((((k4-\[Phi]1-\[Phi]2))^(2))))+((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k4-\[Phi]2-\[Phi]3) (\[Phi]2-\[Phi]3) (k4-\[Phi]1-\[Phi]4) (\[Phi]1-\[Phi]4))/((((k4-\[Phi]1-\[Phi]2))^(2)) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4)))))/((k3-\[Phi]2) (-k4+\[Phi]1+\[Phi]2))),((2 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (\[Phi]2-\[Phi]3) (k3-\[Phi]4) (k4-\[Phi]1-\[Phi]4) (\[Phi]1-\[Phi]4))/(\[Alpha]\[Phi]13 (k3-\[Phi]2) (k4-\[Phi]1-\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))),((2 \[Alpha]\[Phi]14 ((\[Alpha]\[Phi]23)^(2)) (k4-\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]3) (k3-\[Phi]4) (k4-\[Phi]1-\[Phi]4) (\[Phi]1-\[Phi]4))/(\[Alpha]\[Phi]13 (k3-\[Phi]2) (((k4-\[Phi]1-\[Phi]2))^(2)) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))),((4 ((\[Alpha]\[Phi]14)^(2)) ((\[Alpha]\[Phi]23)^(2)) (((\[Phi]2-\[Phi]3))^(2)) (k3-\[Phi]4) (((k4-\[Phi]1-\[Phi]4))^(2)) (((\[Phi]1-\[Phi]4))^(2)))/(((\[Alpha]\[Phi]13)^(2)) (k3-\[Phi]2) (((k4-\[Phi]1-\[Phi]2))^(2)) (((\[Phi]1-\[Phi]3))^(2)) (((\[Phi]2-\[Phi]4))^(2)))),-((2 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k4-\[Phi]2-\[Phi]3) (\[Phi]2-\[Phi]3) (k3-\[Phi]4) (k4-\[Phi]1-\[Phi]4) (\[Phi]1-\[Phi]4))/((k3-\[Phi]2) (((k4-\[Phi]1-\[Phi]2))^(2)) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))),-((4 ((\[Alpha]\[Phi]14)^(2)) \[Alpha]\[Phi]23 (\[Phi]2-\[Phi]3) (k3-\[Phi]4) (k4-\[Phi]1-\[Phi]4) (\[Phi]1-\[Phi]4) (k4-\[Phi]2-\[Phi]4))/(\[Alpha]\[Phi]13 (k3-\[Phi]2) (((k4-\[Phi]1-\[Phi]2))^(2)) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))),((2 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (\[Phi]2-\[Phi]3) (-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k4-\[Phi]1-\[Phi]3) (k4-\[Phi]2-\[Phi]4))/((((k4-\[Phi]1-\[Phi]2))^(2))))+((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k4-\[Phi]2-\[Phi]3) (\[Phi]2-\[Phi]3) (k4-\[Phi]1-\[Phi]4) (\[Phi]1-\[Phi]4))/((((k4-\[Phi]1-\[Phi]2))^(2)) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4)))) (k3-\[Phi]4) (k4-\[Phi]1-\[Phi]4) (\[Phi]1-\[Phi]4))/(\[Alpha]\[Phi]13 (k3-\[Phi]2) (k4-\[Phi]1-\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))),-((\[Alpha]\[Phi]13 (k3-\[Phi]3) (-k4+\[Phi]2+\[Phi]3))/((k3-\[Phi]1) (-k4+\[Phi]1+\[Phi]2))),-((\[Alpha]\[Phi]13 \[Alpha]\[Phi]23 (k3-\[Phi]3) (k4-\[Phi]1-\[Phi]3) (-k4+\[Phi]2+\[Phi]3))/((k3-\[Phi]1) (k4-\[Phi]1-\[Phi]2) (-k4+\[Phi]1+\[Phi]2))),-((2 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (\[Phi]2-\[Phi]3) (-k4+\[Phi]2+\[Phi]3) (k4-\[Phi]1-\[Phi]4) (\[Phi]1-\[Phi]4))/((k3-\[Phi]1) (k4-\[Phi]1-\[Phi]2) (-k4+\[Phi]1+\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))),((((\[Alpha]\[Phi]13)^(2)) (k3-\[Phi]3) (k4-\[Phi]2-\[Phi]3) (-k4+\[Phi]2+\[Phi]3))/((k3-\[Phi]1) (k4-\[Phi]1-\[Phi]2) (-k4+\[Phi]1+\[Phi]2))),((2 \[Alpha]\[Phi]13 \[Alpha]\[Phi]14 (k3-\[Phi]3) (-k4+\[Phi]2+\[Phi]3) (k4-\[Phi]2-\[Phi]4))/((k3-\[Phi]1) (k4-\[Phi]1-\[Phi]2) (-k4+\[Phi]1+\[Phi]2))),-((\[Alpha]\[Phi]13 (k3-\[Phi]3) (-k4+\[Phi]2+\[Phi]3) (-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k4-\[Phi]1-\[Phi]3) (k4-\[Phi]2-\[Phi]4))/((((k4-\[Phi]1-\[Phi]2))^(2))))+((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k4-\[Phi]2-\[Phi]3) (\[Phi]2-\[Phi]3) (k4-\[Phi]1-\[Phi]4) (\[Phi]1-\[Phi]4))/((((k4-\[Phi]1-\[Phi]2))^(2)) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4)))))/((k3-\[Phi]1) (-k4+\[Phi]1+\[Phi]2))),-((2 \[Alpha]\[Phi]14 (k3-\[Phi]4) (-k4+\[Phi]2+\[Phi]4))/((k3-\[Phi]1) (-k4+\[Phi]1+\[Phi]2))),-((2 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k4-\[Phi]1-\[Phi]3) (k3-\[Phi]4) (-k4+\[Phi]2+\[Phi]4))/((k3-\[Phi]1) (k4-\[Phi]1-\[Phi]2) (-k4+\[Phi]1+\[Phi]2))),-((4 ((\[Alpha]\[Phi]14)^(2)) \[Alpha]\[Phi]23 (\[Phi]2-\[Phi]3) (k3-\[Phi]4) (k4-\[Phi]1-\[Phi]4) (\[Phi]1-\[Phi]4) (-k4+\[Phi]2+\[Phi]4))/(\[Alpha]\[Phi]13 (k3-\[Phi]1) (k4-\[Phi]1-\[Phi]2) (-k4+\[Phi]1+\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))),((2 \[Alpha]\[Phi]13 \[Alpha]\[Phi]14 (k4-\[Phi]2-\[Phi]3) (k3-\[Phi]4) (-k4+\[Phi]2+\[Phi]4))/((k3-\[Phi]1) (k4-\[Phi]1-\[Phi]2) (-k4+\[Phi]1+\[Phi]2))),((4 ((\[Alpha]\[Phi]14)^(2)) (k3-\[Phi]4) (k4-\[Phi]2-\[Phi]4) (-k4+\[Phi]2+\[Phi]4))/((k3-\[Phi]1) (k4-\[Phi]1-\[Phi]2) (-k4+\[Phi]1+\[Phi]2))),-((2 \[Alpha]\[Phi]14 (-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k4-\[Phi]1-\[Phi]3) (k4-\[Phi]2-\[Phi]4))/((((k4-\[Phi]1-\[Phi]2))^(2))))+((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k4-\[Phi]2-\[Phi]3) (\[Phi]2-\[Phi]3) (k4-\[Phi]1-\[Phi]4) (\[Phi]1-\[Phi]4))/((((k4-\[Phi]1-\[Phi]2))^(2)) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4)))) (k3-\[Phi]4) (-k4+\[Phi]2+\[Phi]4))/((k3-\[Phi]1) (-k4+\[Phi]1+\[Phi]2))),((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (\[Phi]2-\[Phi]3) (-k4+\[Phi]2+\[Phi]3) (k3-\[Phi]4) (k4-\[Phi]1-\[Phi]4) (\[Phi]1-\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (k4-\[Phi]1-\[Phi]2) (-k4+\[Phi]1+\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4)))-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (-k4+\[Phi]1+\[Phi]3) (k3-\[Phi]4) (-k4+\[Phi]2+\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (((-k4+\[Phi]1+\[Phi]2))^(2)))),((\[Alpha]\[Phi]23 (k4-\[Phi]1-\[Phi]3) (((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (\[Phi]2-\[Phi]3) (-k4+\[Phi]2+\[Phi]3) (k3-\[Phi]4) (k4-\[Phi]1-\[Phi]4) (\[Phi]1-\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (k4-\[Phi]1-\[Phi]2) (-k4+\[Phi]1+\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4)))-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (-k4+\[Phi]1+\[Phi]3) (k3-\[Phi]4) (-k4+\[Phi]2+\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (((-k4+\[Phi]1+\[Phi]2))^(2))))))/(k4-\[Phi]1-\[Phi]2)),((2 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (\[Phi]2-\[Phi]3) (k4-\[Phi]1-\[Phi]4) (\[Phi]1-\[Phi]4) (((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (\[Phi]2-\[Phi]3) (-k4+\[Phi]2+\[Phi]3) (k3-\[Phi]4) (k4-\[Phi]1-\[Phi]4) (\[Phi]1-\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (k4-\[Phi]1-\[Phi]2) (-k4+\[Phi]1+\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4)))-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (-k4+\[Phi]1+\[Phi]3) (k3-\[Phi]4) (-k4+\[Phi]2+\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (((-k4+\[Phi]1+\[Phi]2))^(2))))))/(\[Alpha]\[Phi]13 (k4-\[Phi]1-\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4))),-((\[Alpha]\[Phi]13 (k4-\[Phi]2-\[Phi]3) (((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (\[Phi]2-\[Phi]3) (-k4+\[Phi]2+\[Phi]3) (k3-\[Phi]4) (k4-\[Phi]1-\[Phi]4) (\[Phi]1-\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (k4-\[Phi]1-\[Phi]2) (-k4+\[Phi]1+\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4)))-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (-k4+\[Phi]1+\[Phi]3) (k3-\[Phi]4) (-k4+\[Phi]2+\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (((-k4+\[Phi]1+\[Phi]2))^(2))))))/(k4-\[Phi]1-\[Phi]2)),-((2 \[Alpha]\[Phi]14 (k4-\[Phi]2-\[Phi]4) (((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (\[Phi]2-\[Phi]3) (-k4+\[Phi]2+\[Phi]3) (k3-\[Phi]4) (k4-\[Phi]1-\[Phi]4) (\[Phi]1-\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (k4-\[Phi]1-\[Phi]2) (-k4+\[Phi]1+\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4)))-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (-k4+\[Phi]1+\[Phi]3) (k3-\[Phi]4) (-k4+\[Phi]2+\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (((-k4+\[Phi]1+\[Phi]2))^(2))))))/(k4-\[Phi]1-\[Phi]2)),(-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k4-\[Phi]1-\[Phi]3) (k4-\[Phi]2-\[Phi]4))/((((k4-\[Phi]1-\[Phi]2))^(2))))+((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k4-\[Phi]2-\[Phi]3) (\[Phi]2-\[Phi]3) (k4-\[Phi]1-\[Phi]4) (\[Phi]1-\[Phi]4))/((((k4-\[Phi]1-\[Phi]2))^(2)) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4)))) (((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (\[Phi]2-\[Phi]3) (-k4+\[Phi]2+\[Phi]3) (k3-\[Phi]4) (k4-\[Phi]1-\[Phi]4) (\[Phi]1-\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (k4-\[Phi]1-\[Phi]2) (-k4+\[Phi]1+\[Phi]2) (\[Phi]1-\[Phi]3) (\[Phi]2-\[Phi]4)))-((4 \[Alpha]\[Phi]14 \[Alpha]\[Phi]23 (k3-\[Phi]3) (-k4+\[Phi]1+\[Phi]3) (k3-\[Phi]4) (-k4+\[Phi]2+\[Phi]4))/((k3-\[Phi]1) (k3-\[Phi]2) (((-k4+\[Phi]1+\[Phi]2))^(2)))))}};


(* === Q scalar product === *)


(* === (SU(2) 1-site)\[CenterDot](SU(3) 2-sites) === *)

ClearAll[OrthogRen22B];
OrthogRen22B[QHL_,QL_,QHR_,QR_,{1,1,0,0},{1,1,0,0}]:=2 \[Pi] I Residue[(QHL[[4]]*Wronskian[{QR[[1]],QR[[2]]},z]*(((1)/((((z-w1))^(2 )) )))),{z,0}]/((Limit[QHL[[4]],z->w1])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z],z->w1]))*(1/((Limit[QHL[[3]],z->0,Direction->"\<FromAbove\>"])*(Limit[D[QHL[[3]],z],z->0,Direction->"\<FromAbove\>"])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z],z->0,Direction->"\<FromAbove\>"])*
(Limit[QHL[[3]],z->1,Direction->"\<FromAbove\>"])*(Limit[D[QHL[[3]],z],z->1,Direction->"\<FromAbove\>"])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z],z->1,Direction->"\<FromAbove\>"])))*
Det[{{2 \[Pi] I Residue[((QHL[[3]] D[QR[[1]],z])/(z))*(((1)/((((z-w1))^(1)) (((z-w2))^(1))))),{z,0}],
2 \[Pi] I Residue[((QHL[[3]] QR[[1]])/(((z)^(2))))*(((1)/((((z-w1))^(1)) (((z-w2))^(1))))),{z,0}],
2 \[Pi] I Residue[((QHL[[3]] D[QR[[1]],z])/(-1+z))*(((1)/((((z-w1))^(1)) (((z-w2))^(1))))),{z,0}],
2 \[Pi] I Residue[((QHL[[3]] QR[[1]])/((((-1+z))^(2))))*(((1)/((((z-w1))^(1)) (((z-w2))^(1))))),{z,0}]},

{2 \[Pi] I Residue[((QHL[[3]] D[QR[[2]],z])/(z))*(((1)/((((z-w1))^(1)) (((z-w2))^(1))))),{z,0}],
2 \[Pi] I Residue[((QHL[[3]] QR[[2]])/(((z)^(2))))*(((1)/((((z-w1))^(1)) (((z-w2))^(1))))),{z,0}],
2 \[Pi] I Residue[((QHL[[3]] D[QR[[2]],z])/(-1+z))*(((1)/((((z-w1))^(1)) (((z-w2))^(1))))),{z,0}],
2 \[Pi] I Residue[((QHL[[3]] QR[[2]])/((((-1+z))^(2))))*(((1)/((((z-w1))^(1)) (((z-w2))^(1))))),{z,0}]},


{2 \[Pi] I Residue[((QHL[[3]] D[QR[[1]],z])/(z))*(((1)/((((z-w1))^(1)) (((z-w2))^(1))))),{z,1}],
2 \[Pi] I Residue[((QHL[[3]] QR[[1]])/(((z)^(2))))*(((1)/((((z-w1))^(1)) (((z-w2))^(1))))),{z,1}],
2 \[Pi] I Residue[((QHL[[3]] D[QR[[1]],z])/(-1+z))*(((1)/((((z-w1))^(1)) (((z-w2))^(1))))),{z,1}],
2 \[Pi] I Residue[((QHL[[3]] QR[[1]])/((((-1+z))^(2))))*(((1)/((((z-w1))^(1)) (((z-w2))^(1))))),{z,1}]},

{2 \[Pi] I Residue[((QHL[[3]] D[QR[[2]],z])/(z))*(((1)/((((z-w1))^(1)) (((z-w2))^(1))))),{z,1}],
2 \[Pi] I Residue[((QHL[[3]] QR[[2]])/(((z)^(2))))*(((1)/((((z-w1))^(1)) (((z-w2))^(1))))),{z,1}],
2 \[Pi] I Residue[((QHL[[3]] D[QR[[2]],z])/(-1+z))*(((1)/((((z-w1))^(1)) (((z-w2))^(1))))),{z,1}],
2 \[Pi] I Residue[((QHL[[3]] QR[[2]])/((((-1+z))^(2))))*(((1)/((((z-w1))^(1)) (((z-w2))^(1))))),{z,1}]}}];


(* === (SU(2) 2-sites) \[CenterDot] (SU(3) 1-site) === *)

ClearAll[OrthogRen22C];
OrthogRen22C[QHL_,QL_,QHR_,QR_,{1,1,0,0},{1,1,0,0}]:=(2 \[Pi] I Residue[(QHL[[4]]*Wronskian[{QR[[1]],QR[[2]]},z]*(((1)/((((z-w1))^(2)) (((z-w2))^(2)))))),{z,w1}]*2 \[Pi] I Residue[(QHL[[4]]*Wronskian[{QR[[1]],QR[[2]]},z]*z*(((1)/((((z-w1))^(2)) (((z-w2))^(2)))))),{z,w2}]-2 \[Pi] I Residue[(QHL[[4]]*Wronskian[{QR[[1]],QR[[2]]},z]*(((1)/((((z-w1))^(2)) (((z-w2))^(2)))))),{z,w2}]*2 \[Pi] I Residue[(QHL[[4]]*Wronskian[{QR[[1]],QR[[2]]},z]*z*(((1)/((((z-w1))^(2)) (((z-w2))^(2)))))),{z,w1}])/((Limit[QHL[[4]],z->w1])*(Limit[QHL[[4]],z->w2])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z],z->w1])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z],z->w2]))*
Det[{{2 \[Pi] I Residue[QHL[[3]]*D[QR[[1]],z]*((1)/(z))*(((1)/((((z-w1))^(1))))),{z,0}],2 \[Pi] I Residue[QHL[[3]]*QR[[1]]*((1)/(((z)^(2))))*(((1)/((((z-w1))^(1))))),{z,0}]},

{2 \[Pi] I Residue[QHL[[3]]*D[QR[[2]],z]*((1)/(z))*(((1)/((((z-w1))^(1))))),{z,0}],2 \[Pi] I Residue[QHL[[3]]*QR[[2]]*((1)/(((z)^(2))))*(((1)/((((z-w1))^(1))))),{z,0}]}}]/((Limit[QHL[[3]],z->0,Direction->"\<FromAbove\>"])*(Limit[D[QHL[[3]],z],z->0,Direction->"\<FromAbove\>"])*(Limit[Wronskian[{QR[[1]],QR[[2]]},z],z->0,Direction->"\<FromAbove\>"]));
