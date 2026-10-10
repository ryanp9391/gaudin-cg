ClearAll[GTPatternsRec,GTPatterns,GTWeight,GTIndex,GTTransitionMatrix,GTGenerator];

GTPatternsRec[rows_]:=Module[{next},
 next=Prepend[rows,#]&/@Tuples[
   Table[Range[rows[[1,i+1]],rows[[1,i]]],{i,1,Length[rows[[1]]]-1}]];
 If[Length[next[[1,1]]]==1,Reverse@Sort[next],Reverse@Sort@Flatten[GTPatternsRec/@next,1]]
];
GTPatterns[rep_List]:=GTPatterns[rep]=GTPatternsRec[{rep}];
GTIndex[rep_List]:=GTIndex[rep]=AssociationThread[GTPatterns[rep],Range[Length[GTPatterns[rep]]]];
GTWeight[pat_]:=Table[
 Sum[pat[[k,i]],{i,1,k}]-If[k==1,0,Sum[pat[[k-1,i]],{i,1,k-1}]],{k,1,4}];

GTTransitionMatrix[k_Integer,direction_String,rep_List]:=GTTransitionMatrix[k,direction,rep]=Module[
 {patterns=GTPatterns[rep],index=GTIndex[rep],dimension,rules},
 dimension=Length[patterns];
 rules=Reap[Do[
   Module[{pat=patterns[[source]],denominator,numerator,coefficient,target},
    Do[
     denominator=Product[pat[[k,i]]-i-pat[[k,j]]+j,
       {j,Complement[Range[k],{i}]}];
     If[direction=="Raise",
      numerator=Product[pat[[k,i]]-i-pat[[k+1,j]]+j,{j,1,k+1}];
      coefficient=-numerator/denominator;
      target=ReplacePart[pat,{k,i}->pat[[k,i]]+1],
      numerator=Product[pat[[k,i]]-i-pat[[k-1,j]]+j,{j,1,k-1}];
      coefficient=numerator/denominator;
      target=ReplacePart[pat,{k,i}->pat[[k,i]]-1]
     ];
     If[!PossibleZeroQ[coefficient]&&KeyExistsQ[index,target],
      Sow[{index[target],source}->coefficient]],
     {i,1,k}]
   ],{source,1,dimension}]][[2]];
 SparseArray[If[rules==={}, {}, First[rules]],{dimension,dimension}]
];

GTGenerator[a_Integer,b_Integer,rep_List]:=GTGenerator[a,b,rep]=Which[
 a==b,SparseArray@DiagonalMatrix[GTWeight[#][[a]]&/@GTPatterns[rep]],
 b==a+1,GTTransitionMatrix[a,"Raise",rep],
 a==b+1,GTTransitionMatrix[b,"Lower",rep],
 b>a+1,GTGenerator[a,a+1,rep].GTGenerator[a+1,b,rep]-
   GTGenerator[a+1,b,rep].GTGenerator[a,a+1,rep],
 a>b+1,GTGenerator[a,a-1,rep].GTGenerator[a-1,b,rep]-
   GTGenerator[a-1,b,rep].GTGenerator[a,a-1,rep]
];

ClearAll[RectangularIrreps,ChannelParameters,TensorSectorSupport,TotalRaisingBlock];
RectangularIrreps[s_Integer?NonNegative,p_Integer?NonNegative]:=Module[
 {total=s+p,depth=Min[s,p]},
 Flatten[Table[{total-a,total-b,b,a},{a,0,depth},{b,a,depth}],1]
];
ChannelParameters[s_Integer?NonNegative,p_Integer?NonNegative,mu_List]:=Module[
 {total=s+p,a,b},
 If[Length[mu]=!=4,Return[$Failed]];
 a=total-mu[[1]];b=total-mu[[2]];
 If[TrueQ[IntegerQ[a]&&IntegerQ[b]&&0<=a<=b<=Min[s,p]&&
   mu=={total-a,total-b,b,a}],{a,b},$Failed]
];

TensorSectorSupport[rep1_List,rep2_List,weight_List]:=TensorSectorSupport[rep1,rep2,weight]=Module[
 {weights1=GTWeight/@GTPatterns[rep1],weights2=GTWeight/@GTPatterns[rep2],d2},
 d2=Length[weights2];
 Flatten@Reap[Do[If[TrueQ[weights1[[i]]+weights2[[j]]==weight],
   Sow[(i-1)d2+j]],{i,Length[weights1]},{j,Length[weights2]}]][[2]]
];

TotalRaisingBlock[rep1_List,rep2_List,mu_List,k_Integer]:=Module[
 {d1=Length[GTPatterns[rep1]],d2=Length[GTPatterns[rep2]],source,target,root,operator},
 source=TensorSectorSupport[rep1,rep2,mu];
 root=UnitVector[4,k]-UnitVector[4,k+1];
 target=TensorSectorSupport[rep1,rep2,mu+root];
 If[target==={},Return[SparseArray[{}, {0,Length[source]}]]];
 operator=KroneckerProduct[GTGenerator[k,k+1,rep1],IdentityMatrix[d2,SparseArray]]+
   KroneckerProduct[IdentityMatrix[d1,SparseArray],GTGenerator[k,k+1,rep2]];
 operator[[target,source]]
];

ClearAll[SparseMatrixPower,FirstSiteTuple,FirstSiteBra,FirstSiteOverlap,ClosedFormOverlap,ClosedFormPolynomial];
SparseMatrixPower[m_,0]:=IdentityMatrix[Length[m],SparseArray];
SparseMatrixPower[m_,n_Integer?Positive]:=MatrixPower[m,n];
FirstSiteTuple[a_Integer,b_Integer,r_Integer]:={r,a-r,b-r,r};
FirstSiteBra[rep_List,{n13_,n14_,n23_,n24_}]:=Module[
 {dimension=Length[GTPatterns[rep]],bra},
 bra=UnitVector[dimension,1];
 bra . SparseMatrixPower[GTGenerator[1,3,rep],n13] .
   SparseMatrixPower[GTGenerator[1,4,rep],n14] .
   SparseMatrixPower[GTGenerator[2,3,rep],n23] .
   SparseMatrixPower[GTGenerator[2,4,rep],n24]
];
FirstSiteOverlap[state_,rep1_List,rep2_List,n_List]:=Module[
 {bra1=FirstSiteBra[rep1,n],top2},
 top2=UnitVector[Length[GTPatterns[rep2]],1];
 Flatten[KroneckerProduct[bra1,top2]].state
];
ClosedFormOverlap[b_Integer?NonNegative,r_Integer?NonNegative]:=(-1)^r/Binomial[b,r];
ClosedFormPolynomial[a_Integer?NonNegative,b_Integer?NonNegative]:=
 Expand[u23^(b-a)(u14 u23-u13 u24)^a/(a! b!)];

ClearAll[DirectHWSData,VerifyCGFormula];
Options[DirectHWSData]={"NumericalDiagonalization"->True};
DirectHWSData[s_Integer?Positive,p_Integer?Positive,mu_List,OptionsPattern[]]:=Module[
 {rep1={s,s,0,0},rep2={p,p,0,0},parameters,a,b,d1,d2,support,blocks,constraint,
  nullspace,sectorVector,state,tuples,normalization,overlaps,expected,residuals,
  raisingResidual,gram,eigenvalues,eigenvectors,eigenPosition,numericalVector,exactUnit,diagonalizationError},
 parameters=ChannelParameters[s,p,mu];
 If[parameters===$Failed,Return[Failure["IrrepNotInProduct",<|
   "RequestedRepresentation"->mu,"AvailableRepresentations"->RectangularIrreps[s,p]|>]]];
 {a,b}=parameters;
 d1=Length[GTPatterns[rep1]];d2=Length[GTPatterns[rep2]];
 support=TensorSectorSupport[rep1,rep2,mu];
 blocks=Table[TotalRaisingBlock[rep1,rep2,mu,k],{k,1,3}];
 constraint=If[Total[First/@(Dimensions/@blocks)]==0,
   SparseArray[{}, {0,Length[support]}],Join@@blocks];
 nullspace=If[First[Dimensions[constraint]]==0,
   {UnitVector[Length[support],1]},NullSpace[Normal[constraint]]];
 If[Length[nullspace]=!=1,Return[Failure["HWSNotUnique",<|
   "SectorDimension"->Length[support],"Nullity"->Length[nullspace]|>]]];
 sectorVector=First[nullspace];
 state=SparseArray[Thread[support->sectorVector],d1 d2];
 tuples=Table[FirstSiteTuple[a,b,r],{r,0,a}];
 normalization=FirstSiteOverlap[state,rep1,rep2,First[tuples]];
 sectorVector=FullSimplify[sectorVector/normalization];
 state=SparseArray[Thread[support->sectorVector],d1 d2];
 overlaps=AssociationThread[tuples,
   FullSimplify[FirstSiteOverlap[state,rep1,rep2,#]]&/@tuples];
 expected=AssociationThread[tuples,Table[ClosedFormOverlap[b,r],{r,0,a}]];
 residuals=FullSimplify[Values[overlaps]-Values[expected]];
 raisingResidual=If[First[Dimensions[constraint]]==0,{},constraint.sectorVector];
 diagonalizationError=Missing["NotComputed"];
 If[TrueQ[OptionValue["NumericalDiagonalization"]],
  gram=Transpose[constraint].constraint;
  {eigenvalues,eigenvectors}=Eigensystem[N[Normal[gram]]];
  eigenPosition=First@Ordering[Abs[eigenvalues],1];
  numericalVector=Normalize[eigenvectors[[eigenPosition]]];
  exactUnit=Normalize[N[sectorVector]];
  diagonalizationError=Chop[1-Abs[Conjugate[numericalVector].exactUnit]]
 ];
 <|"InputRepresentations"->{rep1,rep2},"Representation"->mu,
   "ChannelParameters"->{a,b},"TensorDimension"->d1 d2,
   "WeightSectorDimension"->Length[support],"RaisingConstraintDimensions"->Dimensions[constraint],
   "HighestWeightNullity"->Length[nullspace],"SectorSupport"->support,
   "NormalizedSectorVector"->sectorVector,"FullSparseState"->state,
   "Normalization"->HoldForm[bra[0,a,b,0,0,0,0,0].Psi==1],
   "DirectOverlaps"->overlaps,"ClosedFormOverlaps"->expected,
   "OverlapResiduals"->residuals,
   "GeneratingPolynomial"->ClosedFormPolynomial[a,b],
   "RaisingResidual"->raisingResidual,
   "SmallestGramEigenvalue"->If[ValueQ[eigenvalues],Min[Abs[eigenvalues]],Missing["NotComputed"]],
   "DiagonalizedVectorProjectiveError"->diagonalizationError,
   "HighestWeightExact"->TrueQ[And@@(PossibleZeroQ/@Flatten[raisingResidual])],
   "FormulaExact"->TrueQ[And@@(PossibleZeroQ/@residuals)]|>
];

VerifyCGFormula[s_Integer?Positive,p_Integer?Positive,mu_List]:=Module[
 {timing,result},
 {timing,result}=AbsoluteTiming[DirectHWSData[s,p,mu]];
 If[FailureQ[result],Return[result]];
 Join[result,<|"EvaluationSeconds"->timing|>]
];

