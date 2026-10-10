# L=2 diagonal-frame overlap formula — design (2026-10-01)

Goal: generalize row 19 of `Companion-B-C-diagonalize-clean.md` (Ω_zw·Ψ_diag = det N / normalization)
to two sites, adjoint⊗adjoint [2,1,0]², working purely in the diagonal-twist frame.
Notebook: `Paul/Mathematica/new-conventions/L2-diagonal-frame-overlap.wb`.

Approach: derive det N from Baxter/dual adjointness with a two-site bracket (6 varying integrals →
6×6 det N), check same-twist orthogonality, and only then fit normalization (c₀, reference factors,
⟨x₀|y₀⟩). No ansatz-search loops: stop and report after orthogonality, and again if the
site-product reference factors fail.

Sections: (1) setup, ℒ=u+hE_ij, sites θ1,θ2, T^z=T·diag(z); (2) t_a, qdet, numeric diagonalization
per weight sector; (3) Q_i, Q^j from Baxter eqs with ν→∏_s ν, degrees, QQ with F=∏_s∏_n(u−θ_s+nh);
(4) w-vacuum bra = product of rotated HW covectors, eigen of t(G'), G'=P_aux G(w) P_aux⁻¹;
(5) two-site bracket, 6×6 det N, orthogonality for all 64 states — STOP & report;
(6) overlap formula with site-product reference factors, c₀ symbolic; (7) pass/fail table.

Numeric point: h=1/2, θ1=2/7, θ2=−3/5, z=(3/2,5/4,2/3), w=(5/2,7/3,3/5), 60 digits.
