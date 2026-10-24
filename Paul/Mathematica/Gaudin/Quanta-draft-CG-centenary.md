# A Century-Old Bookkeeping Problem Gets an Elegant Solution From Spinning Magnets

*Draft for a Quanta-style feature. Working title alternatives: "The Hidden Music Behind Quantum Addition"; "How a Solvable Magnet Learned to Add Symmetries".*

---

**Since 1927, physicists have used tables of "Clebsch–Gordan coefficients" to combine quantum states. For the larger symmetries of particle physics, such as SU(4), these tables are slow to compute, and no shortcut has been available. A new method computes them hundreds of times faster by combining ideas from integrable systems with an old theorem of representation theory.**

---

When two quantum objects meet, their properties do not simply add. Combine two electrons, each with spin one-half, and you do not get "spin one" and stop there. You get a spin-one triplet and a spin-zero singlet, and each of these is a specific blend of the original states. The numbers that specify the blend are called Clebsch–Gordan coefficients. They appear wherever quantum mechanics meets symmetry: in atomic spectra, nuclear shell models, the quark model, lattice gauge theory and quantum chemistry.

The name honours two 19th-century German mathematicians, Alfred Clebsch and Paul Gordan, who met the same algebra in the 1860s and 1870s while studying invariants of polynomial forms. The coefficients entered physics almost a century ago. In 1927 the young Eugene Wigner, in a paper in *Zeitschrift für Physik* on the term structure of atomic spectra, used group theory to explain how angular momenta combine. [CHECK: Wigner 1927, Z. Phys. 43, 624; confirm which coefficients appear there explicitly versus in his 1931 book.] Hermann Weyl had just finished building the representation theory of continuous groups (1925–26). Within a few years "group pest", as some physicists grumbled, had become standard equipment. In 2027 the coefficients reach their centenary as working tools of quantum physics.

For ordinary rotations, the symmetry group SU(2), the problem has long been solved. Closed formulas date from Wigner and Giulio Racah in the 1940s, and every physics student meets the tables. The difficulty begins with larger symmetries.

## When the labels run out

Particle physics needs SU(3), the symmetry of quark colour and of Gell-Mann's "eightfold way", and also SU(4) and beyond. For these groups a new problem appears. When two representations are combined, the same outcome can occur several times. Two octets of SU(3) contain *two* copies of an octet. The standard quantum numbers cannot tell those copies apart, so the coefficients are not even defined until someone chooses how to separate the copies.

This is the "missing label problem", and it has a fifty-year history. For SU(3) one extra label is missing. Jerry Draayer and Yoshimi Akiyama dealt with it in 1973 with a careful prescription that remains the basis of today's SU(3) software. For SU(4), three labels are missing, and nothing comparable exists: no specialised method has ever been developed for it. It is known what the missing labels *should* be. Mathematicians proved in the 1990s that suitable operators exist. Physicists and mathematicians studying integrable systems found a natural choice of such operators, discussed below. What nobody had was a fast way to *use* them.

So the general codes take the brute-force route. The most widely used one, written by a group in Munich in 2011, builds a large system of linear equations whose size grows with the full dimension of the combined space, and then extracts a null space. It is reliable, but the cost rises steeply. On a laptop, the highest-weight SU(4) coefficients for one standard family of products take about 85 seconds at moderate size. One step larger, the program needs tens of gigabytes of memory.

## A magnet that can be solved exactly

The new work, by a team at DESY in Hamburg led by Volker Schomerus, with Till Bargheer, who heads an Emmy Noether research group there, the postdoctoral researcher Paul Ryan and the doctoral student Matteo Baroni, sets out to make the computation fast for any SU(N), and especially for SU(4), where no specialised method exists. It combines two ingredients, one taken from integrable systems and one from classical representation theory.

In 1976 the French physicist Michel Gaudin wrote down a model of spins interacting through a central magnetic field, a toy version of the electron spins in a quantum dot that couple to the nuclear spins of the surrounding lattice. Gaudin's magnet belongs to the same family as the spin chain that Hans Bethe solved in 1931. Like that chain, it is *integrable*: it has as many independent conserved quantities as it has degrees of freedom, and those quantities all commute with one another. This makes such models unusually rigid and unusually solvable. Their conserved quantities can be packed into a single generating object, from which Bethe equations, "Baxter Q-functions" and differential operators called opers follow in a uniform way.

The first ingredient is a known connection. Take a Gaudin magnet with just *two* sites and put one of the representations to be combined on each site. The total symmetry is then the very product that the Clebsch–Gordan problem asks about. Gaudin's conserved quantities commute with that total symmetry, so they act inside each set of repeated copies, and they commute with each other, so they can be diagonalised simultaneously. Their joint eigenvalues are the missing labels. Mathematicians, among them Evgeny Mukhin, Vitaly Tarasov, Alexander Varchenko and Leonid Rybnikov, worked out this structure, the "Bethe algebra", in depth. As a computational tool, though, it was no better than the brute-force approach: writing down the operators on the full space and diagonalising them is just as expensive as the null-space method.

"The integrable model tells you what the right labels are. The challenge was to make that knowledge cheap," Schomerus said. [QUOTE PLACEHOLDER — to be approved by Volker Schomerus.]

For SU(3) the two-site model supplies exactly one extra operator, a cubic one, which matches the one missing label. For SU(4) the model supplies exactly three operators, one cubic and two quartic, which matches the three missing labels. They come out as the coefficients of a single expansion, a "column determinant" built from the model's Lax matrix, so no one has to choose them.

## A surprise in the symmetry

Nature still had one subtlety in reserve. The natural guess is that the cubic operator alone would separate the SU(4) copies, as its counterpart does for SU(3). It does not. Swapping the two sites of the magnet flips the sign of the cubic operator, up to a constant, so whenever a representation is combined with itself some copies are forced to share an eigenvalue. In the tests the team ran, one of the quartic operators breaks every such tie, and the three operators together always give each copy its own distinct label. [Status: verified on every case computed so far; not yet proven in general.]

The eigenvalues also connect to a classical object. Each copy corresponds to a fourth-order differential equation, an "oper", whose polynomial solutions have degrees fixed by the representations. In the language of integrability, every Clebsch–Gordan copy is a solution of Bethe equations. The team checked this correspondence exactly, by an independent algebraic solve, on hundreds of cases.

## Only a slice

The second ingredient, and the heart of the new work, is what makes the method fast. A 1960s theorem of Kalyanapuram Parthasarathy, R. Ranga Rao and Veeravalli Varadarajan (PRV), refined by Bertram Kostant, says that the highest-weight vector of each combined representation is determined by a small "slice" of the full state, namely its components along one weight space of one factor. That slice must satisfy a few simple "string conditions". The team's "slice engine" works only with this small piece and rebuilds the full vector by a short recursion. The key step is to reduce the Gaudin operators so that they act on the slice directly. The missing labels are then computed in the same small space, and the full space is never built. The reduced operators turn out to be universal polynomials in the labels of one site, written once and reused for every product.

The large linear system disappears. On the same laptop, in the benchmark family where the Munich code needs about 85 seconds, the slice engine finishes in under a tenth of a second, roughly 900 times faster. Against an improved version of the Munich program the gap is about 40 times, and it widens to several hundred times at the next sizes up. Cases that no version of the earlier code could run, with hundreds of thousands of basis states, take seconds. The two approaches agree to machine precision on every comparison made: 42 multiplicity-free cases vector by vector, and 100 sectors with repeated copies.

For SU(3), where the method was first developed, a further simplification emerged: the reduced operator becomes a *tridiagonal* (Jacobi) matrix, the same structure that appears in classical orthogonal polynomials. The missing label is then a root of a polynomial, and each coefficient is that polynomial evaluated at the root.

## Old mathematics, new pattern

Some of the most suggestive results are still conjectures. When the team followed the Gaudin operators around the lattice of states that mathematicians use to label SU(4) bases, a signed sum over the 24 elements of the Weyl group, the symmetry group of the root system, exactly counted the physical copies in all 832 sectors tested. This looks like a shadow of the Bernstein–Gelfand–Gelfand resolution, a cornerstone of 1970s representation theory. A proof would explain *why* the integrable model knows which eigenvectors are physical. [Status: exact on every tested case; proof sketched, one step missing.]

"A hundred years ago, Wigner showed physicists that symmetry organises quantum mechanics," said [EXTERNAL EXPERT — e.g. a representation theorist or integrability expert for an independent comment]. "What's striking here is that integrability, which is a stronger kind of symmetry, organises the symmetry itself."

The practical result is that SU(4), which for fifty years has had no method of its own, now has one that is faster than the general-purpose codes by orders of magnitude. The same construction works for any SU(N). The remaining open questions are proofs for all representations, the extension to larger N, and the full tables of "descendant" coefficients needed by lattice and nuclear physicists.

A century after Wigner, the coefficients that tell quantum states how to combine can be read off from a solvable magnet, and now, quickly.

---

### Fact-check / editorial notes (not for publication)

- **Anniversary framing.** Clebsch and Gordan's invariant-theory work dates to about 1868–75 (Gordan's finiteness memoir, Crelle 69, 1869; Clebsch, *Theorie der binären algebraischen Formen*, 1872), so roughly 150 years. The about-100-year mark is the arrival in quantum mechanics: Weyl's representation theory (1925–26), Wigner, *Z. Phys.* 43, 624 (1927), Weyl's *Gruppentheorie und Quantenmechanik* (1928) and Wigner's book (1931). **2027 = centenary of Wigner 1927** is the cleanest hook. Confirm before using it as a headline claim.
- Gaudin model: 1976 (Gaudin, *J. Physique* 37). Bethe ansatz: 1931.
- Benchmark numbers come from `Paul/Python/gaudin_gl4/RESULTS-G7-munich-comparison.md` (same machine, single thread, S = 8: 84.6 s vs 0.093 s; repaired Munich 3.74 s). Quote them as "on one laptop, for highest-weight vectors", not as full CG tables.
- Prior art the article must not overclaim against: Klink–Ton-That (labelling by invariants is a theorem for U(N)); MTV/Rybnikov (Bethe algebra, opers, completeness). See `Literature-check-SU4-CG.md`.
- Claims marked [Status] are numerical/exact-on-tested-cases, not theorems.
- Munich code: Alex, Kalus, Huckleberry, von Delft, *J. Math. Phys.* 52 (2011).
