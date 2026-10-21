(* ::Package:: *)

(* Gaudin-Qsystem-scale.wl -- WP-Q0 of round 6 (agent Q0, 2026-10-21): the untwisted (eps = 0) Q-system of the two-site sl(3) Gaudin
   model for EVERY irrep of Lam x Lam, Lam = {2a,a,0}, at scale (a = 7, 10, 14 and the test cases a = 3, 4).

   LOAD ORDER (all three with Get, in Global context, on the shared kernel):
       Get[".../Paul/Mathematica/Gaudin/Experiments/gl3-GT-tools.wl"];
       Get[".../Minimal_complete_working-fast.wl"];
       Get[".../Paul/Mathematica/Gaudin/Experiments/Gaudin-Qsystem-scale.wl"];
       runAll[7]        (* reproduces the a = 7 table of Gaudin-Qsystem-scale.md (shipped baseline, corrected solver, fastEig); ~20 s *)
       runAll[14, "FastEig" -> False]   (* 421 irreps, ~4 min including the (slow, wrong) shipped baseline; "Shipped" -> False skips it *)
       qsWriteRef[{14,7,0}, ".../ref/qsystem-spectrum-1470.txt", 100]   (* the reference files *)
   This file never edits the two inputs; every variant it needs is defined here with the prefix qs.

   WHAT IT DOES
   (1) qsIrreps / qsMult / qsDim / qsNW: Racah-Speiser enumeration of the irreps mu of Lam x Lam (and the weight-space dimension n_W).
   (2) qsPencil: the three Baxter equations (k = 1,2,3) in exact polynomial-coefficient form.  With Q = l q, l = (u-1)^-lam1 u^-mu1,
         (u^3 (u-1)^3) B[l q] / l  =  p3 q''' + p2 q'' + p1 q' + (p0 + V u (u-1)) q ,   p0..p3 polynomials in u with rational coefficients;
       V (the vertex-operator eigenvalue) enters ONLY through the multiplication operator V u(u-1).  The systems
       (A0 + V A1) c + (b0 + V b1) = 0 for the coefficients c of the monic echelon q_k (computeV's gauge) are assembled directly on the
       monomials u^n (banded, exact integers) -- no symbolic Together/CoefficientArrays.
   (3) qsConsPoly / qsSolve: the consistency polynomials P_k(V) (V admissible for k iff [A|b] has a null vector with last entry 1) by
       exact Det interpolation of random square projections + GCD (the charPoly idea of Minimal_complete_working-fast.wl), and
       the spectrum = roots of  GCD(P1, P2, P3)  -- ALL THREE (P3 is computed adaptively: only when deg GCD(P1,P2) differs from the
       Racah-Speiser multiplicity; option "ExactThird" -> True forces it).  The shipped pencilSolve uses only P1, P2; that is WRONG on
       C(a,2) + a + 1 irreps of Lam x Lam (29/113 at a = 7, 106/421 at a = 14; see Gaudin-Qsystem-scale.md), because (i) when M1 = 0
       (mu1 = 4a, a+1 irreps) system 1 has no unknowns and CoefficientArrays[eqs, {}] silently takes V itself as the unknown
       (P1 == 0, GCD = P2 alone) and (ii) for the C(a,2) irreps with mu3 >= a+1, mu2 <= 2a-1, P1 and P2 have common roots at which the
       degree-M3 equation is inconsistent.  qsSolve then solves the three linear systems at every root (QR + back substitution at
       working precision) and records the relative residuals of all three.  qsPencilSolvePatched is the minimal two-change patch of the
       shipped function (verified at a <= 7).
   (4) qsLaurent: the untwisted Laurent data (H_s, a_s, b_s) from the Baxter coefficients (map below), and qsEta / qsAlpha: W4's closed forms.
   (5) cross-checks: qsCompareFastEig (a <= 7), qsCompareW4 (W4's exact BRANCH records of ref/untwisted-exact-*.txt), degree rule.
   (6) qsShippedTimed: step-by-step timing of the SHIPPED pipeline (computeV's Baxter assembly + pencilSolve), for the "five slowest" study.
   (7) qsWriteRef: the reference files ref/qsystem-spectrum-<l1><l2>0.txt.

   CONVENTIONS (identical to gl3-GT-tools.wl): site 1 (rep lam) at u = theta1 = 1, site 2 (rep mu = lam here) at u = theta2 = 0, theta12 = 1;
   findweight[lam,mu][nu] = nu - |nu|/3 + |lam+mu|/3, so nu = mu_sector - mu_sector[[3]] {1,1,1} gives E = mu_sector (sum E = |lam|+|mu|);
   Baxter:  Q''' + h1 Q'' + (h2 + h1') Q' + (h3 + h2'/2 + h1''/3) Q = 0,
     h1 = |lam|/(u-1) + |mu|/u,  h2 = hh220/(u-1)^2 + hh211/((u-1)u) + hh020/u^2,
     h3 = hh330/(u-1)^3 + hh321/((u-1)^2 u) + hh312/((u-1)u^2) + hh303/u^3,  hh321 = V (the pencil parameter),
   q_k monic of degree M_k = lam1 + mu1 - E_k + (k-1) (echelon: q2 has no u^M1 term, q3 no u^M1, u^M2 term).

   MAP  Baxter coefficients -> Laurent data (derived from  Baxter operator = -(formal adjoint of  d^3 - S1 d^2 + S2 d - S3) ):
       h1 = S1,   h2 = S2 + S1',   h3 = S3 + S2'/2 + S1''/6      <=>   S2 = h2 - h1',   S3 = h3 - h2'/2 + h1''/3 ,
       S2 = Sum_s [H_s/(u-th_s) + c2/(u-th_s)^2],   S3 = Sum_s [b_s/(u-th_s) + a_s/(u-th_s)^2]   (no constant, no triple pole; checked),
   which gives, for theta = (1, 0) and any (lam, mu) with lam3 = mu3 = 0:
       H1 = -H2 = hh211,   a1 = hh321 + hh211/2 = V + H1/2,   a2 = -hh312 - hh211/2,   b1 = -b2 = hh312 - hh321.
   (Hence a1 + a2 + theta1 b1 + theta2 b2 = 0 and b1 + b2 = 0 hold identically.)  CLOSED FORMS (equal sites Lam = {l1,l2,0}, sector weight E; proved as
   polynomial identities in l1, l2, E1, E2 with Simplify, and checked on every branch): hh211 = eta = |Lam|^2 - (C2(E) - 2 C2(Lam))/2 and
   hh312 + V = 2 alpha - eta with alpha = (E1+2)(E2+1)E3/2.  Hence  H1 = eta,  a1 = V + eta/2,  b1 = 2 (alpha - a1),  a2 = a1 - 2 alpha,  b2 = -b1:
   the Laurent data of a branch are fixed by the single number V; multiplicity-1: b = 0, a1 = alpha (W4); sum_k b1 = 0 over the copies of every sector. *)

ClearAll["Global`qs*"]; ClearAll[runAll];

$qsDir = If[StringQ[$InputFileName] && $InputFileName =!= "", DirectoryName[$InputFileName],
   "C:\\Users\\pryan\\git-projects\\gaudin-cg\\Paul\\Mathematica\\Gaudin\\Experiments\\"];
$qsRefDir = ExpandFileName[FileNameJoin[{$qsDir, "..", "ref"}]];

(* ================= 1. representation theory bookkeeping ================= *)
qsDim[m_List] := (m[[1]] - m[[2]] + 1) (m[[2]] - m[[3]] + 1) (m[[1]] - m[[3]] + 2)/2;

(* weight multiplicities of V_lam from the GT patterns (k | p1 p2 | lam): weight {k, p1+p2-k, |lam|-p1-p2} *)
qsWtMults[lam_List] := qsWtMults[lam] = Module[{cnt = <||>, w},
   Do[w = {k, p1 + p2 - k, Total[lam] - p1 - p2}; cnt[w] = Lookup[cnt, Key[w], 0] + 1,
    {p1, lam[[2]], lam[[1]]}, {p2, lam[[3]], lam[[2]]}, {k, p2, p1}]; cnt];

(* Racah-Speiser: multiplicity of V_mu in V_lam x V_lam2 *)
qsMult[lam_List, lam2_List, mu_List] := Module[{rho = {2, 1, 0}, wm = qsWtMults[lam2]},
   Sum[Signature[s] Lookup[wm, Key[(mu + rho)[[s]] - rho - lam], 0], {s, Permutations[{1, 2, 3}]}]];

(* all irreps {mu, mult} of V_lam x V_lam2, in decreasing lexicographic order of mu *)
qsIrreps[lam_List, lam2_List] := Module[{tot = Total[lam] + Total[lam2], top = lam[[1]] + lam2[[1]], res = {}, m},
   Do[Do[With[{m3 = tot - m1 - m2}, If[m3 >= 0 && m3 <= m2, m = qsMult[lam, lam2, {m1, m2, m3}]; If[m > 0, AppendTo[res, {{m1, m2, m3}, m}]]]],
     {m2, 0, Min[m1, tot - m1]}], {m1, Ceiling[tot/3], top}];
   Reverse[SortBy[res, First]]];

(* n_W = dim of the weight-mu subspace of V_lam x V_lam2 (what the chain-side null space has to work in) *)
qsNW[lam_List, lam2_List, mu_List] := Module[{w1 = qsWtMults[lam], w2 = qsWtMults[lam2]},
   Total[KeyValueMap[#2 Lookup[w2, Key[mu - #1], 0] &, w1]]];

(* ================= 2. the Baxter equation: coefficients, polynomial operator, exact pencils ================= *)
(* the seven constants hh of computeV (copied from gl3-GT-tools.wl; validated by the Baxter residual of computeV's own q's) *)
qsHH[lam_List, mu_List, nu_List, Vv_] := Module[{en = findweight[lam, mu][nu], l1, l2, l3, m1, m2, m3, T},
   {l1, l2, l3} = lam; {m1, m2, m3} = mu; T = Total[lam] + Total[mu];
   <|"220" -> -l1 + l1 l2 + l3 + l1 l3 + l2 l3, "020" -> -m1 + m1 m2 + m3 + m1 m3 + m2 m3,
    "330" -> l1/3 - (2 l2)/3 - l1 l2 + l3/3 + l2 l3 + l1 l2 l3, "303" -> m1/3 - (2 m2)/3 - m1 m2 + m3/3 + m2 m3 + m1 m2 m3,
    "321" -> Vv,
    "312" -> -Vv - (1 + en[[1]]) en[[2]]^2 + l2 + l1 l2 - l2 l3 - l1 l2 l3 + m2 + m1 m2 - (1 + en[[1]]) en[[2]] (1 + en[[1]] - T) - m2 m3 - m1 m2 m3,
    "211" -> -en[[1]]^2 - en[[2]]^2 + 2 l1 + l2 - l1 l2 - l1 l3 - l2 l3 + 2 m1 + m2 - m1 m2 - m1 m3 - m2 m3 + en[[2]] (-1 + T) + en[[1]] (-2 - en[[2]] + T)|>];

qsH[lam_List, mu_List, nu_List, Vv_] := Module[{hh = qsHH[lam, mu, nu, Vv], z},
   {Total[lam]/(z - 1) + Total[mu]/z,
     hh["220"]/(z - 1)^2 + hh["211"]/((z - 1) z) + hh["020"]/z^2,
     hh["330"]/(z - 1)^3 + hh["321"]/((z - 1)^2 z) + hh["312"]/((z - 1) z^2) + hh["303"]/z^3} /. z -> u];

(* coefficient functions of  u^3 (u-1)^3 B[l q]/l = p3 q''' + p2 q'' + p1 q' + p0 q  (all polynomials; V only in p0, as + V u (u-1)) *)
qsOpPolys[lam_List, mu_List, nu_List] := Module[{Vs, h1, h2, h3, g, c, pl, den = u^3 (u - 1)^3, ok},
   {h1, h2, h3} = qsH[lam, mu, nu, Vs];
   g = -lam[[1]]/(u - 1) - mu[[1]]/u;   (* l'/l *)
   c = {D[g, {u, 2}] + 3 g D[g, u] + g^3 + h1 (D[g, u] + g^2) + (h2 + D[h1, u]) g + h3 + D[h2, u]/2 + D[h1, {u, 2}]/3,
     3 (D[g, u] + g^2) + 2 g h1 + h2 + D[h1, u], 3 g + h1, 1};
   pl = Cancel[Together[den #1]] & /@ c;
   ok = And @@ (PolynomialQ[#1, u] & /@ pl) && Expand[pl[[1]] - (pl[[1]] /. Vs -> 0) - Vs u (u - 1)] === 0;
   If[! ok, Return[$Failed]];
   {CoefficientList[pl[[1]] /. Vs -> 0, u], CoefficientList[pl[[2]], u], CoefficientList[pl[[3]], u], CoefficientList[pl[[4]], u]}];

(* exact pencils: for k = 1,2,3  (A0 + V A1) c + (b0 + V b1) = 0, stored as augmented matrices aug0 = [A0|b0], aug1 = [A1|b1]
   (rows = powers u^0 .. u^(M_k+2); the top row u^(M_k+3) vanishes identically).  S[[k]] = exponents of the unknown c[k,n]. *)
qsPencil[lam_List, mu_List, nu_List] := Module[{P, en, M, Nmax, falling, colOp, colV, S},
   P = qsOpPolys[lam, mu, nu]; If[P === $Failed, Return[$Failed]];
   en = findweight[lam, mu][nu];
   M = {lam[[1]] + mu[[1]] - en[[1]], lam[[1]] + mu[[1]] - en[[2]] + 1, lam[[1]] + mu[[1]] - en[[3]] + 2};
   Nmax = M[[3]] + 3;
   falling[n_, j_] := Product[n - i, {i, 0, j - 1}];
   colOp[n_] := Module[{v = ConstantArray[0, Nmax + 1]},
     Do[With[{ff = falling[n, j], Pj = P[[j + 1]]},
       If[ff != 0, Do[v[[n - j + s + 1]] += ff Pj[[s + 1]], {s, 0, Length[Pj] - 1}]]], {j, 0, 3}]; v];
   colV[n_] := Module[{v = ConstantArray[0, Nmax + 1]}, v[[n + 3]] += 1; v[[n + 2]] -= 1; v];
   S = {Range[0, M[[1]] - 1], DeleteCases[Range[0, M[[2]] - 1], M[[1]]], DeleteCases[Range[0, M[[3]] - 1], M[[1]] | M[[2]]]};
   <|"M" -> M, "E" -> en, "S" -> S,
    "aug0" -> Table[Transpose[Append[colOp /@ S[[k]], colOp[M[[k]]]]][[1 ;; M[[k]] + 3]], {k, 3}],
    "aug1" -> Table[Transpose[Append[colV /@ S[[k]], colV[M[[k]]]]][[1 ;; M[[k]] + 3]], {k, 3}]|>];

(* ================= 3. consistency polynomials and the spectrum ================= *)
(* P(V) = GCD of nproj random square projections of det([A|b](V)), by exact interpolation at V = 0..n; returns {P, tDet, tGcd}.
   nproj = 2 is the charPoly of Minimal_complete_working-fast.wl; nproj = 1 keeps the projection artefacts (random extra roots), which the
   GCD with the polynomial of another k removes (the artefacts of different k and different projections do not coincide). *)
qsConsPoly[aug0_, aug1_, Vs_, nproj_ : 2] := Module[{m, n, proj, ps, p, tD = 0., tG = 0., t},
   {m, n} = Dimensions[aug0];
   proj[] := Module[{R = RandomInteger[{-9, 9}, {n, m}], RA0, RA1, pts = Range[0, n], dets},
     RA0 = R . aug0; RA1 = R . aug1;
     {t, dets} = AbsoluteTiming[Det[RA0 + #1 RA1] & /@ pts]; tD += t;
     Expand@InterpolatingPolynomial[Transpose[{pts, dets}], Vs]];
   ps = Table[proj[], {nproj}]; p = First[ps];
   If[nproj > 1, {tG, p} = AbsoluteTiming[PolynomialGCD @@ ps]];
   {p, tD, tG}];

(* full coefficient list (ascending powers, length M_k + 1) of q_k from the unknowns c *)
qsFullCoeffs[Sk_List, Mk_Integer, c_List] := Module[{v = ConstantArray[0, Mk + 1]}, v[[Mk + 1]] = 1; Do[v[[Sk[[i]] + 1]] = c[[i]], {i, Length[Sk]}]; v];

(* solve the three systems at one root v (working precision wp); returns {qs, relative residuals} *)
qsSolveRoot[pp_, v_, solver_ : "QR"] := Module[{rows},
   rows = Table[Module[{aug = pp["aug0"][[k]] + v pp["aug1"][[k]], A, b, c, res, scale, qq, rr},
      A = aug[[All, 1 ;; -2]]; b = aug[[All, -1]];
      (* QR + back substitution is ~25x faster than LeastSquares on these high-precision systems and equal to 1e-60 *)
      c = If[pp["S"][[k]] === {}, {}, Quiet@If[solver === "QR", {qq, rr} = QRDecomposition[A]; LinearSolve[rr, qq . (-b)], LeastSquares[A, -b]]];
      res = Max[Abs[If[c === {}, b, A . c + b]]];
      (* matrix-norm scale of the system (no division by a numerical zero): |aug0| + |v| |aug1|, times max(1, |c|) *)
      scale = (Max[Abs[pp["aug0"][[k]]]] + Abs[v] Max[Abs[pp["aug1"][[k]]]]) Max[1, If[c === {}, 1, Max[Abs[c]]]];
      {qsFullCoeffs[pp["S"][[k]], pp["M"][[k]], c], N[res/scale, 6]}], {k, 3}];
   {rows[[All, 1]], rows[[All, 2]]}];

(* The corrected solver.  Options:
     "ExactThird" -> True:  spectrum = roots of GCD(P1,P2,P3) always (exact; this is what the first full sweeps used);
                      False: roots of GCD(P1,P2) only, spurious ones dropped afterwards by the residual of the third system (not exact);
                      Automatic (default): adaptive -- P3 is computed only when deg GCD(P1,P2) differs from the Racah-Speiser multiplicity
                      (the genuine roots are always among the roots of GCD(P1,P2), so deg = mult leaves nothing to remove); retried with new
                      random projections up to 3 times if deg GCD != mult.
     "Projections" -> number of random square projections per consistency polynomial (2 = charPoly of the fast file; 1 default).
     "Solver" -> "QR" (default) | "LeastSquares" for the linear systems at the roots; "WorkingPrecision" -> digits; "Seed".
   Returns an Association: M, mexp, degP (degrees of the P_k computed), degG12, degG (of the GCD used), usedP3, attempts, roots (sorted by Re),
   sols ({V, {q1,q2,q3} coefficient lists, rel. residuals}), t (step times: build, det, gcd, nsolve, solve), pp (the pencils). *)
Options[qsSolve] = {"WorkingPrecision" -> 60, "ExactThird" -> Automatic, "Seed" -> 20261021, "Solver" -> "QR", "Projections" -> 1};
qsSolve[lam_List, mu_List, nu_List, OptionsPattern[]] := Module[{wp = OptionValue["WorkingPrecision"], ex3 = OptionValue["ExactThird"],
    np = OptionValue["Projections"], solver = OptionValue["Solver"], seed = OptionValue["Seed"], pp, en, mexp, tB, tDet = 0., tGcd = 0.,
    Vs = Unique["Vq"], P = <||>, getP, G12, G, deg12, usedP3 = False, tR, roots, tS, sols, t, attempt = 0},
   en = findweight[lam, mu][nu]; mexp = qsMult[lam, mu, en];
   {tB, pp} = AbsoluteTiming[qsPencil[lam, mu, nu]];
   getP[k_] := Module[{r = qsConsPoly[pp["aug0"][[k]], pp["aug1"][[k]], Vs, np]}, tDet += r[[2]]; tGcd += r[[3]]; r[[1]]];
   While[True,
    attempt++;
    BlockRandom[SeedRandom[seed + attempt - 1];
     P[1] = getP[1]; P[2] = getP[2];
     {t, G12} = AbsoluteTiming[PolynomialGCD[P[1], P[2]]]; tGcd += t;
     deg12 = Exponent[G12, Vs]; G = G12;
     If[ex3 === True || (ex3 === Automatic && deg12 != mexp), usedP3 = True; P[3] = getP[3]; {t, G} = AbsoluteTiming[PolynomialGCD[G12, P[3]]]; tGcd += t]];
    If[ex3 =!= Automatic || Exponent[G, Vs] == mexp || attempt >= 3, Break[]]];
   {tR, roots} = AbsoluteTiming[If[Exponent[G, Vs] < 1, {}, Sort[Vs /. NSolve[G == 0, Vs, WorkingPrecision -> wp], Re[#1] < Re[#2] &]]];
   {tS, sols} = AbsoluteTiming[Table[Module[{qr = qsSolveRoot[pp, N[v, wp], solver]}, <|"V" -> N[v, wp], "q" -> qr[[1]], "resid" -> qr[[2]]|>], {v, roots}]];
   <|"M" -> pp["M"], "mexp" -> mexp, "degP" -> Table[If[KeyExistsQ[P, k], Exponent[P[k], Vs], Missing[]], {k, 3}], "degG12" -> deg12, "degG" -> Exponent[G, Vs],
    "usedP3" -> usedP3, "attempts" -> attempt, "roots" -> roots, "sols" -> sols, "G" -> (G /. Vs -> qsVar),
    "t" -> <|"build" -> tB, "det" -> tDet, "gcd" -> tGcd, "nsolve" -> tR, "solve" -> tS|>, "pp" -> pp|>];

(* ================= 4. Laurent data and W4's closed forms ================= *)
qsC2[l_List] := Sum[l[[i]] (l[[i]] + 4 - 2 i), {i, 3}];
qsEta[Lam_List, mu_List] := Total[Lam]^2 - (qsC2[mu] - 2 qsC2[Lam])/2;     (* equal site representations Lam *)
qsAlpha[mu_List] := (mu[[1]] + 2) (mu[[2]] + 1) mu[[3]]/2;                 (* W4's closed form, multiplicity-1 sectors *)

(* principal parts of S2, S3 at theta = (1, 0), symbolic in qsVar; also the checks that S2, S3 equal their principal parts *)
qsLaurentSym[lam_List, mu_List, nu_List] := qsLaurentSym[lam, mu, nu] = Module[{h1, h2, h3, S2, S3, th = {1, 0}, pr},
   {h1, h2, h3} = qsH[lam, mu, nu, qsVar];
   S2 = Together[h2 - D[h1, u]]; S3 = Together[h3 - D[h2, u]/2 + D[h1, {u, 2}]/3];
   pr = Table[<|"H" -> SeriesCoefficient[S2, {u, th[[s]], -1}], "c2" -> SeriesCoefficient[S2, {u, th[[s]], -2}], "c3" -> SeriesCoefficient[S2, {u, th[[s]], -3}],
      "b" -> SeriesCoefficient[S3, {u, th[[s]], -1}], "a" -> SeriesCoefficient[S3, {u, th[[s]], -2}], "t3" -> SeriesCoefficient[S3, {u, th[[s]], -3}]|>, {s, 2}];
   <|"H" -> pr[[All, "H"]], "a" -> pr[[All, "a"]], "b" -> pr[[All, "b"]], "c2" -> pr[[All, "c2"]], "c3" -> pr[[All, "c3"]], "t3" -> pr[[All, "t3"]],
    "S2rem" -> Together[S2 - Sum[pr[[s]]["H"]/(u - th[[s]]) + pr[[s]]["c2"]/(u - th[[s]])^2, {s, 2}]],
    "S3rem" -> Together[S3 - Sum[pr[[s]]["b"]/(u - th[[s]]) + pr[[s]]["a"]/(u - th[[s]])^2, {s, 2}]]|>];

(* Laurent data {H1, a1, b1, H2, a2, b2} at a root value V *)
qsLaurent[lam_List, mu_List, nu_List, Vv_] := Module[{ls = qsLaurentSym[lam, mu, nu]},
   {ls["H"][[1]], ls["a"][[1]], ls["b"][[1]], ls["H"][[2]], ls["a"][[2]], ls["b"][[2]]} /. qsVar -> Vv];

(* ================= 5. cross-checks ================= *)
qsClearComputeVMemo[] := (SubValues[computeV] = Select[SubValues[computeV], ! FreeQ[First[#1], _Pattern] &]; Null);

(* W4's exact BRANCH records: list of <|mu, copy, D0, eta, alpha, beta|> (pairs evaluated; theta12-scaled) *)
qsReadW4[label_String] := Module[{file = FileNameJoin[{$qsRefDir, "untwisted-exact-" <> label <> ".txt"}], lines},
   lines = Select[Import[file, "Lines"], StringStartsQ[#1, "BRANCH\t"] &];
   Map[Function[ln, Module[{f = StringSplit[ln, "\t"], D0, pv},
      D0 = ToExpression[f[[4]]]; pv[p_] := p[[1]] + p[[2]] Sqrt[D0];
      <|"mu" -> ToExpression[f[[2]]], "copy" -> ToExpression[f[[3]]], "D0" -> D0, "eta" -> (pv /@ ToExpression[f[[8]]]),
       "alpha" -> (pv /@ ToExpression[f[[9]]]), "beta" -> (pv /@ ToExpression[f[[10]]])|>]], lines]];

(* compare qsSolve's Laurent data with W4's exact branches for the representation Lam (equal sites); returns per-branch max deviations *)
qsCompareW4[Lam_List] := Module[{label = StringJoin[ToString /@ Lam], w4, mus},
   w4 = qsReadW4[label]; mus = DeleteDuplicates[#1["mu"] & /@ w4];
   Flatten[Table[Module[{nu = mu - mu[[3]] {1, 1, 1}, sol, mine, theirs},
       sol = qsSolve[Lam, Lam, nu]; theirs = SortBy[Select[w4, #1["mu"] === mu &], N[#1["alpha"][[1]]] &];
       mine = Table[qsLaurent[Lam, Lam, nu, s["V"]], {s, sol["sols"]}];
       If[Length[mine] != Length[theirs], {<|"mu" -> mu, "mismatch" -> {Length[mine], Length[theirs]}|>},
        Table[Module[{x = mine[[k]], y = theirs[[k]]},
          <|"mu" -> mu, "copy" -> k, "dH" -> Max[Abs[N[{x[[1]] - y["eta"][[1]], x[[4]] - y["eta"][[2]]}, 20]]],
           "da" -> Max[Abs[N[{x[[2]] - y["alpha"][[1]], x[[5]] - y["alpha"][[2]]}, 20]]],
           "db" -> Max[Abs[N[{x[[3]] - y["beta"][[1]], x[[6]] - y["beta"][[2]]}, 20]]]|>], {k, Length[mine]}]]], {mu, mus}], 1]];

(* roots versus fastEig's V spectrum (sorted); returns max |difference| *)
qsCompareFastEig[Lam_List, mu_List, roots_List] := Module[{nu = mu - mu[[3]] {1, 1, 1}, es},
   es = Sort[Re[fastEig[Lam, Lam][nu, MachinePrecision][[1]]]];
   If[Length[es] != Length[roots], Infinity, If[roots === {}, 0, Max[Abs[es - Sort[N[roots]]]]]]];

(* ================= 6. the shipped pipeline, step by step (verbatim algorithm of computeV + pencilSolve, with timers) ================= *)
(* returns <|"roots", "t" -> <|setup, bax, arrays, det, gcdproj, gcd12, nsolve, ls|>|> *)
qsShippedTimed[lam_List, mu_List, nu_List] := Block[{h, ev, M, l, Q, q, Bax, V, hh, tm = <||>, baxs, qs, sys, P, roots, Ps, projT = 0., dT = 0., gT = 0., lsT = 0., sol, wp = 50, vv},
   tm["setup"] = First@AbsoluteTiming[
     h[1][z_] = Plus @@ lam/(z - 1) + Plus @@ mu/z;
     h[2][z_] = Sum[(hh[2, r, s]/((z - 1)^r*z^s))*Boole[r + s == 2], {r, 0, 2}, {s, 0, 2}];
     h[3][z_] = Sum[(hh[3, r, s]/((z - 1)^r*z^s))*Boole[r + s == 3], {r, 0, 3}, {s, 0, 3}];
     With[{hq = qsHH[lam, mu, nu, V]}, hh[2, 2, 0] = hq["220"]; hh[2, 0, 2] = hq["020"]; hh[3, 3, 0] = hq["330"]; hh[3, 0, 3] = hq["303"]; hh[3, 2, 1] = V; hh[3, 1, 2] = hq["312"]; hh[2, 1, 1] = hq["211"]];
     With[{en = findweight[lam, mu][nu]}, ev[1] = en[[1]]; ev[2] = en[[2]]; ev[3] = en[[3]]];
     M[1] = lam[[1]] + mu[[1]] - ev[1]; M[2] = lam[[1]] + mu[[1]] - ev[2] + 1; M[3] = lam[[1]] + mu[[1]] - ev[3] + 2;
     l[u_] = 1/((-1 + u)^lam[[1]]*u^mu[[1]]); Q[k_][u_] = l[u]*q[k][u];
     q[1][u_] = u^M[1] + Sum[c[1, n]*u^n, {n, 0, M[1] - 1}];
     q[2][u_] = u^M[2] + Sum[c[2, n]*u^n*Boole[n != M[1]], {n, 0, M[2] - 1}];
     q[3][u_] = u^M[3] + Sum[c[3, n]*u^n*Boole[n != M[2]]*Boole[n != M[1]], {n, 0, M[3] - 1}]];
   tm["bax"] = First@AbsoluteTiming[
     Bax[k_] = Numerator[Together[Derivative[1][Q[k]][u]*(h[2][u] + Derivative[1][h[1]][u]) + (h[3][u] + (1/2)*Derivative[1][h[2]][u] + (1/3)*Derivative[2][h[1]][u])*Q[k][u] + h[1][u]*Derivative[2][Q[k]][u] + Derivative[3][Q[k]][u]]];
     baxs = {Bax[1], Bax[2], Bax[3]}; qs = {q[1][u], q[2][u], q[3][u]}];
   (* pencilSolve, verbatim, with timers *)
   tm["arrays"] = First@AbsoluteTiming[
     sys = Table[Module[{eqs = DeleteCases[Flatten[CoefficientList[baxs[[k]], u]], 0], vars},
        vars = Union@Cases[eqs, c[k, _], Infinity]; {vars, Normal@CoefficientArrays[eqs, vars]}], {k, 3}]];
   tm["charPoly"] = First@AbsoluteTiming[Ps = Table[Module[{A = sys[[k, 2, 2]], b = sys[[k, 2, 1]], aug, n, m, p},
        aug = Transpose[Append[Transpose[A], b]]; {m, n} = Dimensions[aug];
        p = Table[Module[{RA = RandomInteger[{-9, 9}, {n, m}] . aug, pts = Range[0, n], dets, tt},
           {tt, dets} = AbsoluteTiming[Det[RA /. V -> #1] & /@ pts]; dT += tt;
           Expand@InterpolatingPolynomial[Transpose[{pts, dets}], V]], {2}];
        gT += First@AbsoluteTiming[p = PolynomialGCD[p[[1]], p[[2]]]]; p], {k, 2}]];
   tm["det"] = dT; tm["gcdproj"] = gT;
   tm["gcd12"] = First@AbsoluteTiming[P = PolynomialGCD @@ Ps];
   tm["nsolve"] = First@AbsoluteTiming[roots = V /. NSolve[P == 0, V, WorkingPrecision -> 50]];
   tm["ls"] = First@AbsoluteTiming[
     sol = Table[Prepend[Table[qs[[k]] /. Quiet@Thread[sys[[k, 1]] -> LeastSquares[sys[[k, 2, 2]] /. V -> vv, -(sys[[k, 2, 1]] /. V -> vv)]], {k, 3}], vv], {vv, roots}]];
   <|"roots" -> roots, "t" -> tm, "nroots" -> Length[roots]|>];

(* The MINIMAL patch of the shipped pencilSolve (two changes, everything else verbatim):
     (1) GCD over ALL THREE consistency polynomials (shipped: sys[[;; 2]] only),
     (2) a guard for a system without unknowns (M1 = 0 makes q1 == 1: CoefficientArrays[eqs, {}] would take V itself as the unknown).
   Use:  Block[{pencilSolve = qsPencilSolvePatched}, computeV[lam, mu][nu]]   (clear the memo first: qsClearComputeVMemo[]). *)
qsPencilSolvePatched[baxs_, qs_, V_] := Module[{sys, P, roots},
   sys = Table[Module[{eqs = DeleteCases[Flatten[CoefficientList[baxs[[k]], u]], 0], vars, ca},
      vars = Union@Cases[eqs, c[k, _], Infinity];
      ca = If[vars === {}, {eqs, ConstantArray[0, {Length[eqs], 0}]}, Normal@CoefficientArrays[eqs, vars]];
      {vars, ca}], {k, Length[baxs]}];
   P = PolynomialGCD @@ (charPoly[#1[[2, 2]], #1[[2, 1]], V] & /@ sys);
   roots = V /. NSolve[P == 0, V, WorkingPrecision -> 50];
   Table[Prepend[Table[qs[[k]] /. If[sys[[k, 1]] === {}, {}, Thread[sys[[k, 1]] -> LeastSquares[sys[[k, 2, 2]] /. V -> v, -(sys[[k, 2, 1]] /. V -> v)]]], {k, Length[qs]}], v], {v, roots}]];

(* counts of roots for every irrep of Lam x Lam with the shipped computeV (shipped = True) or with the patched pencilSolve; returns {mu, mult, #roots, time} rows *)
qsCountRoots[a_Integer, patched_ : False] := Module[{Lam = {2 a, a, 0}},
   qsClearComputeVMemo[];
   Table[Module[{mu = r[[1]], nu, t, sol}, nu = mu - mu[[3]] {1, 1, 1};
      {t, sol} = AbsoluteTiming[Quiet[If[patched, Block[{pencilSolve = qsPencilSolvePatched}, computeV[Lam, Lam][nu]], computeV[Lam, Lam][nu]]]];
      {mu, r[[2]], Length[sol], t}], {r, qsIrreps[Lam, Lam]}]];

(* ================= 7. number formatting and reference files ================= *)
(* nd significant digits, properly rounded; plain decimal for 10^-4 <= |x| < 10^10, else d.ddd*^e (Wolfram InputForm exponent) *)
qsNumStr[x_, nd_ : 30] := Module[{x0, e10, n, s, sgn, ds},
   Which[
    x === 0 || (NumericQ[x] && PossibleZeroQ[x] && Head[x] =!= Real && Head[x] =!= Complex), "0",
    IntegerQ[x], ToString[x],
    Head[x] === Rational, ToString[x, InputForm],
    True,
    x0 = N[x, nd + 10];
    If[Head[x0] === Complex, Return["COMPLEX[" <> ToString[x0, InputForm] <> "]"]];
    If[x0 == 0, Return["0"]];
    sgn = If[x0 < 0, "-", ""]; x0 = Abs[x0];
    e10 = Floor[Log10[x0]]; n = Round[x0 10^(nd - 1 - e10)];
    If[n >= 10^nd, n = Round[n/10]; e10++];
    s = StringJoin[ToString /@ IntegerDigits[n, 10, nd]];
    Which[0 <= e10 <= 9, sgn <> StringTake[s, e10 + 1] <> "." <> StringDrop[s, e10 + 1],
     -4 <= e10 < 0, sgn <> "0." <> StringJoin[ConstantArray["0", -e10 - 1]] <> s,
     True, sgn <> StringTake[s, 1] <> "." <> StringDrop[s, 1] <> "*^" <> ToString[e10]]]];
qsListStr[l_List] := "{" <> StringRiffle[qsNumStr /@ l, ","] <> "}";
qsMuStr[mu_List] := "{" <> StringRiffle[ToString /@ mu, ","] <> "}";

(* chop coefficients that are zero to working precision relative to the largest (structural zeros: gauge or valuation) *)
qsChopList[l_List, wp_] := With[{mx = Max[Abs[l]]}, If[Abs[#1] < 10^-(wp - 12) mx, 0, #1] & /@ l];

(* all irreps: <|mu, mult, nu, sol, laurent|>; solves with qsSolve at working precision wp *)
qsSolveAll[Lam_List, wp_ : 100, progress_ : False] := Table[Module[{mu = r[[1]], nu, sol},
     nu = mu - mu[[3]] {1, 1, 1}; sol = qsSolve[Lam, Lam, nu, "WorkingPrecision" -> wp];
     <|"mu" -> mu, "mult" -> r[[2]], "nu" -> nu, "sol" -> sol|>], {r, qsIrreps[Lam, Lam]}];

qsWriteRef[Lam_List, file_String, wp_ : 100, data_ : Automatic] := Module[{res, stream, T = Total[Lam], c2 = Lam[[2]] (Lam[[1]] + 1), nb, checks, tol},
   res = If[data === Automatic, qsSolveAll[Lam, wp], data];
   nb = Total[#1["mult"] & /@ res];
   (* the checks announced in the header, enforced: nothing is written if one fails *)
   tol = 10^-(wp/2);
   checks = Table[Module[{sols = r["sol"]["sols"], lau},
       lau = Table[qsLaurent[Lam, Lam, r["nu"], s["V"]], {s, sols}];
       <|"mu" -> r["mu"], "count" -> (Length[sols] == r["mult"]), "resid" -> (sols =!= {} && Max[Flatten[sols[[All, "resid"]]]] < tol),
        "eta" -> (Union[lau[[All, 1]]] === {qsEta[Lam, r["mu"]]}), "bsum" -> (sols =!= {} && Abs[Total[lau[[All, 3]]]] < tol),
        "alpha1" -> (r["mult"] != 1 || Abs[lau[[1, 2]] - qsAlpha[r["mu"]]] < tol)|>], {r, res}];
   If[! (And @@ Flatten[Values /@ (KeyDrop[#1, "mu"] & /@ checks)]), Print["qsWriteRef: CHECK FAILED on ", Select[checks, ! (And @@ Values[KeyDrop[#1, "mu"]]) &]]; Return[$Failed]];
   stream = OpenWrite[file, CharacterEncoding -> "UTF-8", PageWidth -> Infinity];
   WriteString[stream, StringRiffle[{
      "# " <> FileNameTake[file] <> "   (agent Q0, Paul/Mathematica/Gaudin/Experiments/Gaudin-Qsystem-scale.wl, 2026-10-21)",
      "# Untwisted (epsilon = 0) Q-system of the two-site gl(3) Gaudin chain, Lambda = " <> qsMuStr[Lam] <> " on both sites, ALL irreps mu of Lambda x Lambda, all multiplicity copies.",
      "# CONVENTIONS (those of Paul/Mathematica/Gaudin/Experiments/gl3-GT-tools.wl = computeV): site 1 (rep Lambda) at u = theta1 = 1, site 2 (rep Lambda) at u = theta2 = 0, theta12 = theta1 - theta2 = 1, UNTWISTED.",
      "#   V = the pencil parameter = eigenvalue of the vertex operator V = 1/2 Sum_{s,t in S3} sgn(s) sgn(t) e1_{s1 t1} e1_{s2 t2} e2_{s3 t3} of the notebook (GT matrices ee[i,j][n], rational non-unitary normalisation) on the HW kets of weight mu;",
      "#   the branch copies k = 1..mult of a sector are the roots of the exact polynomial GCD of the Baxter consistency polynomials P1, P2 (and P3 whenever GCD(P1,P2) has more roots than the multiplicity), degree = mult on every sector,",
      "#   ordered by INCREASING V (real simple roots; this is also the increasing a1 order, copy 1 = smallest a1 as in W4's untwisted-exact-*.txt).  The shipped pencilSolve (GCD of P1,P2 only, and P1 == 0 when M1 = 0) is NOT used: it returns spurious roots.",
      "# Baxter equation: Q''' + h1 Q'' + (h2 + h1') Q' + (h3 + h2'/2 + h1''/3) Q = 0, Q_k = q_k(u)/((u-1)^L1 u^L1) with L1 = Lambda_1 (the gauge l(u) of computeV), k = 1,2,3,",
      "#   h1 = |Lambda|/(u-1) + |Lambda|/u,  h2 = hh220/(u-1)^2 + hh211/((u-1)u) + hh020/u^2,  h3 = hh330/(u-1)^3 + hh321/((u-1)^2 u) + hh312/((u-1)u^2) + hh303/u^3, hh321 = V (formulas: qsHH in the .wl).",
      "#   E = mu (the sector weight; findweight with the central shift). q_k are MONIC polynomials of degree M_k = 2 L1 - mu_k + (k-1), echelon gauge (computeV): q2 has no u^M1 term, q3 has no u^M1 and no u^M2 term.",
      "# Laurent data of the oper (W4/W5 convention): S2 = Sum_s [H_s/(u-theta_s) + c2/(u-theta_s)^2], S3 = Sum_s [b_s/(u-theta_s) + a_s/(u-theta_s)^2], c2 = Lambda2 (Lambda1+1) = " <> ToString[c2] <> ", |Lambda| = " <> ToString[T] <> ".",
      "#   dictionary from the Baxter coefficients: S2 = h2 - h1', S3 = h3 - h2'/2 + h1''/3  (Baxter operator = -(formal adjoint of d^3 - S1 d^2 + S2 d - S3), S1 = h1); for theta = (1,0):",
      "#   H1 = -H2 = hh211 (= W4's eta = |Lambda|^2 - (C2(mu) - 2 C2(Lambda))/2, exact, integer),  a1 = V + H1/2,  a2 = -hh312 - H1/2,  b1 = -b2 = hh312 - V;  a1 + a2 + theta1 b1 + theta2 b2 = 0.",
      "#   CLOSED FORMS (exact polynomial identities in the sector weight, verified symbolically and on every branch written here): H1 = eta;  K := hh312 + V = 2 alpha - eta with alpha = (mu1+2)(mu2+1)mu3/2 (W4's product);",
      "#   hence  b1 = 2 (alpha - a1),  a2 = a1 - 2 alpha,  a1 = V + eta/2: the Laurent data of a branch are fixed by V alone.  Multiplicity-1 sectors: b1 = b2 = 0, a1 = -a2 = alpha, V = alpha - eta/2.",
      "#   The sum over the copies of b1 vanishes on every sector (verified to 1e-54 for all irreps at a = 1..14; expected from the site-exchange symmetry of equal sites), i.e. the MEAN of a1 over the copies of any sector is alpha (mean V = alpha - eta/2).",
      "# NUMBERS: 30 significant digits (rounded), computed at working precision " <> ToString[wp] <> "; plain decimal, or d.ddd*^e for |x| < 1e-4 or >= 1e10 (Python: replace '*^' by 'e'); integers/rationals exact.",
      "# FORMAT: TAB-separated records, first field = tag; lists in braces without spaces; lines starting with # are comments.",
      "# RECORDS: REP {l1,l2,l3} | THETA {theta1,theta2} | NSECTORS n | NBRANCHES n | SECTOR mu mult | BRANCH mu k V H1 a1 b1 H2 a2 b2 | QPOLY mu k {coeffs q1} {coeffs q2} {coeffs q3}",
      "#   SECTOR mu mult : mu = {mu1,mu2,mu3} HW of the irrep (|mu| = 2|Lambda|), mult = Racah-Speiser multiplicity (= number of roots, asserted).",
      "#   BRANCH mu k V H1 a1 b1 H2 a2 b2 : copy k = 1..mult (increasing V); V the pencil root; (H_s, a_s, b_s) the untwisted Laurent data at site s (s = 1: theta = 1, s = 2: theta = 0) at theta12 = 1.",
      "#     (theta12 scaling, W4: H_s = eta_s/theta12, a_s = alpha_s/theta12, b_s = beta_s/theta12^2.)",
      "#   QPOLY mu k {q1} {q2} {q3} : coefficients in ASCENDING powers of u (index 0 = constant term; the last entry of each list is the leading coefficient 1; exact zeros = gauge or valuation) of the monic echelon q_1,q_2,q_3 of branch (mu,k).",
      "# CHECKS performed before writing (see Gaudin-Qsystem-scale.md): #roots = multiplicity on every sector; all three Baxter systems satisfied at every root (relative residual < 1e-" <> ToString[Floor[wp/2]] <> ");",
      "#   H1 = eta exactly; sum over copies of b1 = 0; multiplicity-1: V = alpha - eta/2.",
      "REP\t" <> qsMuStr[Lam], "THETA\t{1,0}", "NSECTORS\t" <> ToString[Length[res]], "NBRANCHES\t" <> ToString[nb]}, "\n"] <> "\n"];
   Do[Module[{mu = r["mu"], nu = r["nu"], sol = r["sol"], sols},
     sols = sol["sols"];
     WriteString[stream, "SECTOR\t" <> qsMuStr[mu] <> "\t" <> ToString[r["mult"]] <> "\n"];
     Do[Module[{s = sols[[k]], lau = qsLaurent[Lam, Lam, nu, sols[[k]]["V"]], qc},
       WriteString[stream, StringRiffle[Join[{"BRANCH", qsMuStr[mu], ToString[k], qsNumStr[s["V"]]}, qsNumStr /@ lau], "\t"] <> "\n"];
       qc = qsChopList[#1, wp] & /@ s["q"];
       WriteString[stream, StringRiffle[Join[{"QPOLY", qsMuStr[mu], ToString[k]}, qsListStr /@ qc], "\t"] <> "\n"]], {k, Length[sols]}]], {r, res}];
   Close[stream]; file];

(* per-irrep table of a list of runAll results (e.g. {runAll[7], runAll[10], runAll[14]}): the roots to 30 digits, counts and times *)
qsWriteTable[results_List, file_String] := Module[{stream},
   stream = OpenWrite[file, CharacterEncoding -> "UTF-8", PageWidth -> Infinity];
   WriteString[stream, StringRiffle[{
      "# " <> FileNameTake[file] <> "   (agent Q0, Gaudin-Qsystem-scale.wl runAll[a], 2026-10-21)",
      "# Pencil roots V (the V-spectrum) of EVERY irrep mu of Lambda x Lambda, Lambda = {2a,a,0}, from the corrected solver qsSolve (conventions of qsystem-spectrum-*.txt: theta = (1,0), theta12 = 1, V = hh321 of computeV).",
      "# RECORDS: LAMBDA a {2a,a,0} nirreps ncopies  |  ROOTS mu mult nroots nShipped M tShippedWall tCorrectedWall nW {V1,...,Vm}",
      "#   mult = Racah-Speiser multiplicity; nroots = number of roots of the corrected solver (= mult on every irrep, asserted); nShipped = number of roots the SHIPPED computeV/pencilSolve returned",
      "#   (wrong where nShipped != mult: see Gaudin-Qsystem-scale.md); M = {M1,M2,M3} degrees of q1,q2,q3; tShippedWall/tCorrectedWall = wall seconds on the shared i7 box (noisy, +-50%);",
      "#   nW = dimension of the weight-mu subspace of V_Lambda x V_Lambda (the chain-side problem size); V_k in increasing order, 30 significant digits ('*^' exponent as in Wolfram InputForm; Python: replace by 'e').",
      "# FORMAT: TAB-separated, lists in braces without spaces; lines starting with # are comments."}, "\n"] <> "\n"];
   Do[Module[{a = r["a"], Lam = r["Lam"]},
      WriteString[stream, "LAMBDA\t" <> ToString[a] <> "\t" <> qsMuStr[Lam] <> "\t" <> ToString[Length[r["table"]]] <> "\t" <> ToString[Total[r["table"][[All, "mult"]]]] <> "\n"];
      Do[WriteString[stream, StringRiffle[{"ROOTS", qsMuStr[t["mu"]], ToString[t["mult"]], ToString[t["nRoots"]], ToString[t["nShipped"]], qsMuStr[t["M"]],
          ToString[If[NumericQ[t["tShipped"]], Round[1000 t["tShipped"]]/1000., "NA"]], ToString[Round[1000 t["tFixed"]]/1000.], ToString[t["nW"]], qsListStr[t["roots"]]}, "\t"] <> "\n"], {t, r["table"]}]], {r, results}];
   Close[stream]; file];

(* markdown rows of a runAll result (the a = 7 table of the report): | mu | m | M | n_W | #roots shipped | #roots corrected | t_ship ms | t_corr ms | V roots | *)
qsMarkdownTable[res_Association] := StringRiffle[Table[StringJoin["| ", StringRiffle[{StringReplace[ToString[t["mu"]], {"{" -> "(", "}" -> ")", " " -> ""}], ToString[t["mult"]],
      StringRiffle[ToString /@ t["M"], ","], ToString[t["nW"]], ToString[t["nShipped"]] <> If[t["nShipped"] =!= t["mult"], " **x**", ""], ToString[t["nRoots"]],
      ToString[If[NumericQ[t["tShipped"]], Round[1000 t["tShipped"]], "NA"]], ToString[Round[1000 t["tFixed"]]],
      StringRiffle[If[Abs[#1 - Round[#1]] < 10^-30, ToString[Round[#1]], qsNumStr[#1, 9]] & /@ t["roots"], ", "]}, " | "], " |"], {t, res["table"]}], "\n"];

(* ================= 8. the driver ================= *)
(* one-line number formatting for the summary text: 4 digits, '*^' exponent; an inexact zero prints as its accuracy bound *)
qsFmt[x_] := Which[x === 0, "0", Head[x] === Real && x == 0 && Accuracy[x] < Infinity, "<1e-" <> ToString[Floor[Accuracy[x]]], True, qsNumStr[x, 4]];
Options[runAll] = {"WorkingPrecision" -> 60, "FastEig" -> Automatic, "Shipped" -> True, "Verbose" -> True, "SolverOptions" -> {}};
(* runAll[a]: Lam = {2a,a,0}; every irrep: shipped computeV count/time, corrected qsSolve (exact GCD of P1,P2,P3), Laurent checks,
   fastEig cross-check (a <= 7).  Returns <|"table" -> list of per-irrep Associations, "summary" -> text, ...|>. *)
runAll[a_Integer, OptionsPattern[]] := Module[{Lam = {2 a, a, 0}, wp = OptionValue["WorkingPrecision"], irr, fe, table, summary, nIrr, nCop, tTot},
   irr = qsIrreps[Lam, Lam]; nIrr = Length[irr]; nCop = Total[Last /@ irr];
   fe = If[OptionValue["FastEig"] === Automatic, a <= 7, TrueQ[OptionValue["FastEig"]]];
   If[OptionValue["Shipped"], qsClearComputeVMemo[]];
   table = Table[Module[{mu = r[[1]], m = r[[2]], nu, ship, tShip, sol, tFix, lau, etaOK, alphaOK, traceOK, bsum, d, rec, fdev},
       nu = mu - mu[[3]] {1, 1, 1};
       If[OptionValue["Shipped"], {tShip, ship} = AbsoluteTiming[Quiet[computeV[Lam, Lam][nu]]], tShip = Missing[]; ship = Missing[]];
       {tFix, sol} = AbsoluteTiming[qsSolve[Lam, Lam, nu, "WorkingPrecision" -> wp, Sequence @@ OptionValue["SolverOptions"]]];
       lau = Table[qsLaurent[Lam, Lam, nu, s["V"]], {s, sol["sols"]}];
       etaOK = (Length[lau] == 0) || (Union[lau[[All, 1]]] === {qsEta[Lam, mu]});
       alphaOK = If[m == 1 && Length[lau] == 1, Abs[N[lau[[1, 2]] - qsAlpha[mu], 20]], Missing[]];
       bsum = If[Length[lau] > 0, Abs[Total[lau[[All, 3]]]], 0];
       fdev = If[fe, qsCompareFastEig[Lam, mu, sol["roots"]], Missing[]];
       <|"mu" -> mu, "mult" -> m, "M" -> sol["M"], "nW" -> qsNW[Lam, Lam, mu], "nShipped" -> If[ship === Missing[], Missing[], Length[ship]], "tShipped" -> tShip,
        "nRoots" -> Length[sol["roots"]], "degG12" -> sol["degG12"], "degG" -> sol["degG"], "usedP3" -> sol["usedP3"], "attempts" -> sol["attempts"], "tFixed" -> tFix, "t" -> sol["t"],
        "roots" -> sol["roots"], "maxRelRes" -> If[sol["sols"] === {}, 0, Max[sol["sols"][[All, "resid"]]]], "etaOK" -> etaOK, "alphaDev" -> alphaOK,
        "bSum" -> bsum, "fastEigDev" -> fdev, "realRoots" -> (Length[sol["roots"]] == 0 || Max[Abs[Im[sol["roots"]]]] < 10^-(wp/2))|>], {r, irr}];
   tTot = <|"shipped" -> If[OptionValue["Shipped"], Total[table[[All, "tShipped"]]], Missing[]], "fixed" -> Total[table[[All, "tFixed"]]]|>;
   summary = StringRiffle[{
      "a = " <> ToString[a] <> ", Lambda = " <> ToString[Lam] <> ": " <> ToString[nIrr] <> " irreps, " <> ToString[nCop] <> " copies, sum m dim = dim^2: " <> ToString[Total[#1[[2]] qsDim[#1[[1]]] & /@ irr] == qsDim[Lam]^2],
      If[OptionValue["Shipped"], "shipped computeV: total " <> qsFmt[tTot["shipped"]] <> " s; #roots != mult on " <> ToString[Count[table, t_ /; t["nShipped"] =!= t["mult"]]] <> "/" <> ToString[nIrr] <> " irreps (M1 = 0: " <> ToString[Count[table, t_ /; t["nShipped"] =!= t["mult"] && t["M"][[1]] == 0]] <> ")", "shipped skipped"],
      "corrected qsSolve (" <> ToString[OptionValue["SolverOptions"]] <> "): total " <> qsFmt[tTot["fixed"]] <> " s; #roots = mult on " <> ToString[Count[table, t_ /; t["nRoots"] === t["mult"]]] <> "/" <> ToString[nIrr] <> "; deg GCD(P1,P2) = mult on " <> ToString[Count[table, t_ /; t["degG12"] === t["mult"]]] <> "/" <> ToString[nIrr] <> " (P3 needed/used on " <> ToString[Count[table, t_ /; TrueQ[t["usedP3"]]]] <> "); deg of the final GCD = mult on " <> ToString[Count[table, t_ /; t["degG"] === t["mult"]]] <> "/" <> ToString[nIrr] <> "; max rel. residual " <> qsFmt[Max[table[[All, "maxRelRes"]]]] <> "; real roots on all: " <> ToString[And @@ table[[All, "realRoots"]]],
      "H1 = eta exactly on " <> ToString[Count[table[[All, "etaOK"]], True]] <> "/" <> ToString[nIrr] <> "; mult-1 sectors: " <> ToString[Count[table, t_ /; t["mult"] == 1]] <> ", max |a1 - alpha| = " <> qsFmt[Max[Cases[table[[All, "alphaDev"]], _?NumericQ]]] <> "; max |sum_k b1| = " <> qsFmt[Max[table[[All, "bSum"]]]],
      If[fe, "fastEig: max |V - eigenvalue| = " <> qsFmt[Max[table[[All, "fastEigDev"]]]], "fastEig skipped"]}, "\n"];
   If[OptionValue["Verbose"], Print[summary]];
   <|"a" -> a, "Lam" -> Lam, "table" -> table, "summary" -> summary, "totals" -> tTot|>];
