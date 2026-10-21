(* ::Package:: *)

(* Fast drop-ins for Minimal_complete_working.nb.
   Load AFTER the notebook's "GT patterns, generators and vertex operator" cells
   (needs ps, J0, ee, findweight, V) and after the computeV cell.

   1. pencilSolve: replaces the final NSolve line of computeV. Also replace every
      Subscript[h,i,j,k] by hh[i,j,k] and add hh to the Block variables -- the
      subscripts live on Subscript, escape the Block, and h[2]/h[3] (immediate Set,
      defined before the coefficients) silently pick up the PREVIOUS call's values.
      Use the last line
        pencilSolve[{Bax[1],Bax[2],Bax[3]}, {q[1][u],q[2][u],q[3][u]}, V]
      Returns one {V, q1, q2, q3} per solution (the original returned {q1, q2}).

   2. fastEig[lam,mu][nu]: replaces computeHWS + V[hws] + the per-sigma Solve of
      eigenvec. Returns {eigenvalues of V, eigenvectors (coefficients on W), W},
      W = list of {i,j} with tp[ps[{lam}][[i]], ps[{mu}][[j]]] of weight nu.
      The eigenvalues coincide with the V's of pencilSolve. *)

(* V admissible iff [A|b] has a null vector with last entry 1: two random square
   projections, exact Det interpolated in V, then GCD removes projection artefacts *)
charPoly[A_, b_, V_] := Module[{aug = Transpose[Append[Transpose[A], b]], n, m, proj},
  {m, n} = Dimensions[aug];
  proj[] := Module[{RA = RandomInteger[{-9, 9}, {n, m}] . aug, pts = Range[0, n]},
    Expand@InterpolatingPolynomial[Transpose[{pts, Det[RA /. V -> #] & /@ pts}], V]];
  PolynomialGCD[proj[], proj[]]];

(* Q1's pencil alone has spurious roots (e.g. 2a+1 for the singlet); GCD with Q2's removes them *)
pencilSolve[baxs_, qs_, V_] := Module[{sys, P, roots},
  sys = Table[Module[{eqs = DeleteCases[Flatten[CoefficientList[baxs[[k]], u]], 0], vars},
      vars = Union@Cases[eqs, c[k, _], Infinity];
      {vars, Normal@CoefficientArrays[eqs, vars]}], {k, Length[baxs]}];
  P = PolynomialGCD @@ (charPoly[#[[2, 2]], #[[2, 1]], V] & /@ sys[[;; 2]]);
  roots = V /. NSolve[P == 0, V, WorkingPrecision -> 50];
  Table[Prepend[Table[qs[[k]] /. Thread[sys[[k, 1]] ->
        LeastSquares[sys[[k, 2, 2]] /. V -> v, -(sys[[k, 2, 1]] /. V -> v)]], {k, Length[qs]}], v],
   {v, roots}]];

(* the 36 terms {coeff, {a,b}, {c,d}, {e,f}} of the notebook's V = sum coeff ee[a,b][1] ee[c,d][1] ee[e,f][2]:
   V = 1/2 sum_{s,t in S3} sgn(s) sgn(t) e1_{s1 t1} e1_{s2 t2} e2_{s3 t3}  (checked term-by-term against the V cell) *)
vTriples = Flatten[Table[{Signature[s] Signature[t]/2, {s[[1]], t[[1]]}, {s[[2]], t[[2]]}, {s[[3]], t[[3]]}},
    {s, Permutations[{1, 2, 3}]}, {t, Permutations[{1, 2, 3}]}], 1];

(* single-site matrices E[[new,old]] read off the notebook's ee[i,j][1] acting on tp[p, ref] *)
siteMats[rep_] := siteMats[rep] = Module[{pats = ps[{rep}], pos, ref, D},
   D = Length[pats]; ref = pats[[1]]; pos = AssociationThread[pats -> Range[D]];
   Association@Table[{i, j} -> SparseArray[Flatten@Table[
         With[{img = Expand[ee[i, j][1][tp[pats[[k]], ref]]]},
          Cases[If[Head[img] === Plus, List @@ img, {img}], cf_. tp[p_, ref] :> ({pos[p], k} -> cf)]], {k, D}], {D, D}],
     {i, 3}, {j, 3}]];

(* per-rep caches shared by all irreps nu of a product: site matrices at the working precision, and the
   site-1 products E_ab . E_cd appearing in V *)
siteMatsP[rep_, Infinity] := siteMats[rep];
siteMatsP[rep_, prec_] := siteMatsP[rep, prec] = N[#, prec] & /@ siteMats[rep];
sitePairs[rep_, prec_] := sitePairs[rep, prec] = With[{E = siteMatsP[rep, prec]},
   Association@Table[{t[[2]], t[[3]]} -> E[t[[2]]] . E[t[[3]]], {t, vTriples}]];

(* outer multiplicity of the gl(3) weight wt in lam x mu (Racah-Speiser, GT weight multiplicities of mu) *)
multRS[lam_, mu_, wt_] := Module[{wm = Counts[Table[J0[k][p], {p, ps[{mu}]}, {k, 3}]], rho = {2, 1, 0}},
  Sum[Signature[s] Lookup[wm, Key[(wt + rho)[[s]] - rho - lam], 0], {s, Permutations[{1, 2, 3}]}]];

(* a two-site vector is a D1 x D2 matrix X, X[[i,j]] = coefficient of tp[p1_i, p2_j]:
   site-1 operators act from the left, site-2 operators as X . E^T.
   prec = Infinity: exact NullSpace.  prec = MachinePrecision: the k = multRS lowest eigenvectors of the
   Gram matrix A^T A by sparse shift-invert Arnoldi (~20x faster at a=7), orthonormal, accurate to ~1e-13;
   aborts if the (k+1)-th eigenvalue doesn't separate (null eigenvalues ~1e-13, the next O(10) in tests) *)
fastEig[lam_, mu_][nu_, prec_ : Infinity] := Module[{p1 = ps[{lam}], p2 = ps[{mu}], E1, E2, wt, w1, w2, byW, W, nW, D1, D2, cons, rows, N0, Xs, P, Y, Vm, es, vecs, k, G, gs},
  E1 = siteMatsP[lam, prec]; E2 = siteMatsP[mu, prec]; D1 = Length[p1]; D2 = Length[p2];
  wt = findweight[lam, mu][nu];
  w1 = Table[J0[k][p], {p, p1}, {k, 3}]; w2 = Table[J0[k][p], {p, p2}, {k, 3}];
  byW = GroupBy[Range[D2], w2[[#]] &];
  W = Join @@ Table[{i, #} & /@ Lookup[byW, Key[wt - w1[[i]]], {}], {i, D1}];
  nW = Length[W];
  (* column w = image of basis pair W[[w]] = (p,q) under total E_ij: E1[[r,p]] at (r,q), E2[[s,q]] at (p,s);
     assembled from the column nonzeros directly (building each D1 x D2 image and flattening was the bottleneck) *)
  cons = Join @@ Table[Module[{col1, col2},
      (* no position collisions: (r,q) = (p,s) would need a diagonal entry of an off-diagonal E_ij *)
      col1 = GroupBy[Most@ArrayRules[E1[ij]], #[[1, 2]] & -> ({#[[1, 1]], #[[2]]} &)];
      col2 = GroupBy[Most@ArrayRules[E2[ij]], #[[1, 2]] & -> ({#[[1, 1]], #[[2]]} &)];
      SparseArray[Flatten[Table[With[{p = W[[w, 1]], q = W[[w, 2]]},
          Join[({(#[[1]] - 1) D2 + q, w} -> #[[2]]) & /@ Lookup[col1, p, {}],
               ({(p - 1) D2 + #[[1]], w} -> #[[2]]) & /@ Lookup[col2, q, {}]]], {w, nW}], 1], {D1 D2, nW}]],
     {ij, {{1, 2}, {2, 3}}}];
  (* dropping the all-zero rows and using OneStepRowReduction: ~10x faster exact NullSpace than the default *)
  rows = Union[First /@ Most[ArrayRules[cons]][[All, 1]]];
  (* no constraints at all (e.g. nu = lam + mu): the whole weight space is highest weight *)
  N0 = Which[
    rows === {}, IdentityMatrix[nW],
    prec === Infinity, NullSpace[cons[[rows]], Method -> "OneStepRowReduction"],
    True,
     k = multRS[lam, mu, wt]; G = Transpose[cons[[rows]]] . cons[[rows]];
     If[k == 0, {},
      (* Arnoldi can skip members of the degenerate null cluster (k+1 or k+5 requested eigenpairs miss
         some at a=7), so ask for k+20 and double until exactly k eigenvalues are ~0; dense once that's most of G *)
      Module[{extra = 20, m, z},
       While[True,
        m = k + extra;
        {gs, N0} = If[2 m >= nW, Eigensystem[Normal[G]],
          Quiet@Eigensystem[G, -m, Method -> {"Arnoldi", "Shift" -> -1.}]];
        gs = Abs[gs]; z = Count[gs, x_ /; x < 10^-8 Max[1, Max[gs]]];
        If[z == k || 2 m >= nW, Break[]]; extra *= 2];
       If[z != k, Print["fastEig: found ", z, " null vectors, expected ", k, " for ", wt]; Abort[]];
       N0[[Ordering[gs][[;; k]]]]]]];
  If[N0 === {}, Return[{{}, {}, W}, Module]];
  Xs = SparseArray[Thread[W -> #], {D1, D2}] & /@ N0;
  P = sitePairs[lam, prec];
  Y = Table[Sum[t[[1]] P[{t[[2]], t[[3]]}] . X . Transpose[E2[t[[4]]]], {t, vTriples}], {X, Xs}];
  (* V preserves the highest-weight space: its mult x mult matrix in the NullSpace basis *)
  Vm = LeastSquares[Transpose[N0], Transpose[Table[Extract[y, W], {y, Y}]]];
  {es, vecs} = Eigensystem[N[Vm, 50]];
  {es, vecs . N0, W}];
