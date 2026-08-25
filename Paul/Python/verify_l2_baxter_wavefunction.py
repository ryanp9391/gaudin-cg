import numpy as np
from numpy.polynomial import Polynomial as P

h = 0.73
theta1 = -0.41
theta2 = 0.86
lam = 1.67
chi = lam + 1 / lam

I2 = np.eye(2, dtype=complex)
E11 = np.array([[1, 0], [0, 0]], dtype=complex)
E22 = np.array([[0, 0], [0, 1]], dtype=complex)
E12 = np.array([[0, -1], [0, 0]], dtype=complex)  # -d_z
E21 = np.array([[0, 0], [-1, 0]], dtype=complex)  # z^2 d_z-z
E = [[E11, E12], [E21, E22]]


def embed(mat, site):
    return np.kron(mat, I2) if site == 0 else np.kron(I2, mat)


def transfer(u, lax_sign=-1):
    sites = []
    for site, theta in enumerate((theta1, theta2)):
        L = [[None, None], [None, None]]
        for i in range(2):
            for j in range(2):
                L[i][j] = ((u - theta) * np.eye(4, dtype=complex) if i == j else np.zeros((4, 4), dtype=complex))
                L[i][j] += lax_sign * h * embed(E[j][i], site)
        sites.append(L)
    # Monodromy L1 L2. Reversing the order is also tested below.
    T = [[sum((sites[0][i][k] @ sites[1][k][j] for k in range(2)),
              np.zeros((4, 4), dtype=complex)) for j in range(2)] for i in range(2)]
    return chi * T[0][0] + T[0][1] - T[1][0]


def transfer_coefficients(lax_sign=-1):
    tm = transfer(-1.0, lax_sign)
    t0 = transfer(0.0, lax_sign)
    tp = transfer(1.0, lax_sign)
    c0 = t0
    c1 = (tp - tm) / 2
    c2 = (tp + tm) / 2 - t0
    return c0, c1, c2


def eigensystem(lax_sign=-1):
    c0, c1, c2 = transfer_coefficients(lax_sign)
    vals, vecs = np.linalg.eig(c0 + 0.371 * c1)
    states = []
    for k in range(4):
        v = vecs[:, k]
        v /= np.linalg.norm(v)
        coeff = np.array([np.vdot(v, c @ v) for c in (c0, c1, c2)])
        residual = max(np.linalg.norm(c @ v - coeff[j] * v) for j, c in enumerate((c0, c1, c2)))
        states.append((coeff, v, residual))
    return states


qtheta = P([theta1 * theta2, -(theta1 + theta2), 1])
qtheta_plus = P([(h - theta1) * (h - theta2), 2 * h - theta1 - theta2, 1])


def shift_poly(p, shift):
    return p(P([shift, 1]))


def baxter_nullvector(tcoeff, degree):
    tpoly = P(tcoeff)
    cols = []
    maxdeg = degree + 2
    for k in range(degree + 1):
        basis = P([0] * k + [1])
        out = (1 / lam) * qtheta_plus * shift_poly(basis, -h)
        out -= tpoly * basis
        out += lam * qtheta * shift_poly(basis, h)
        cols.append(np.pad(out.coef, (0, maxdeg + 1 - len(out.coef))))
    mat = np.column_stack(cols)
    _, singular, vh = np.linalg.svd(mat)
    vector = vh.conj().T[:, -1]
    residual = np.linalg.norm(mat @ vector)
    if abs(vector[-1]) > 1e-12:
        vector /= vector[-1]
    return vector, residual, singular[-1]


def localized_weights(site):
    tj, tk = ((theta1, theta2) if site == 0 else (theta2, theta1))
    rho = 1 - np.exp(2j * np.pi * (tj - tk) / h)
    w0 = rho / (-h * (tj - tk) * (tj - tk - h))
    w1 = rho / (h * (tj - tk) * (tj - tk + h))
    return w0, w1


def local_data(qcoef, shift=-h):
    q = P(qcoef)
    data = []
    for site, tj in enumerate((theta1, theta2)):
        w0, w1 = localized_weights(site)
        # Shifted Q evaluated at tj and tj+h.
        a = w0 * lam ** ((tj + shift) / h) * q(tj + shift)
        b = w1 * lam ** ((tj + h + shift) / h) * q(tj + h + shift)
        data.append((a, b))
    return data


def wavefunction(qcoef, shift=-h):
    data = local_data(qcoef, shift)
    (a1, b1), (a2, b2) = data
    delta = theta2 - theta1
    # Kronecker basis is {1,z2,z1,z1 z2}.
    return np.array([
        delta * a1 * a2,
        (delta + h) * a1 * b2,
        (delta - h) * b1 * a2,
        delta * b1 * b2,
    ], dtype=complex)


def product_wavefunction(qcoef, shift=-h):
    (a1, b1), (a2, b2) = local_data(qcoef, shift)
    return np.array([a1 * a2, a1 * b2, b1 * a2, b1 * b2], dtype=complex)


def normalized_residual(mat, vec, eig):
    return np.linalg.norm(mat @ vec - eig * vec) / np.linalg.norm(vec)


def normalized_left_residual(mat, vec, eig):
    return np.linalg.norm(mat.T @ vec - eig * vec) / np.linalg.norm(vec)


for lax_sign in (-1, 1):
    print(f"\nLAX SIGN {lax_sign:+d}")
    states = eigensystem(lax_sign)
    direct_columns = []
    psi_columns = []
    for idx, (tcoeff, direct, common_residual) in enumerate(states):
        matches = []
        for degree in range(3):
            qcoef, residual, smallest = baxter_nullvector(tcoeff, degree)
            matches.append((degree, residual, qcoef))
        good = [item for item in matches if item[1] < 1e-9]
        best = good[0] if good else min(matches, key=lambda item: item[1])
        degree, tq_residual, qcoef = best
        psi = wavefunction(qcoef)
        product = product_wavefunction(qcoef)
        # Compare against the same transfer polynomial at two generic u values.
        eig0 = sum(tcoeff[k] * 0.23 ** k for k in range(3))
        eig1 = sum(tcoeff[k] * (-0.67) ** k for k in range(3))
        psi_residual = max(
            normalized_residual(transfer(0.23, lax_sign), psi, eig0),
            normalized_residual(transfer(-0.67, lax_sign), psi, eig1),
        )
        psi_left_residual = max(
            normalized_left_residual(transfer(0.23, lax_sign), psi, eig0),
            normalized_left_residual(transfer(-0.67, lax_sign), psi, eig1),
        )
        product_residual = max(
            normalized_residual(transfer(0.23, lax_sign), product, eig0),
            normalized_residual(transfer(-0.67, lax_sign), product, eig1),
        )
        overlap = abs(np.vdot(direct, psi)) / np.linalg.norm(psi)
        direct_columns.append(direct)
        psi_columns.append(psi / np.linalg.norm(psi))
        print({
            "state": idx,
            "tcoeff(I0,I1,chi)": np.real_if_close(tcoeff),
            "common_eigenvector_residual": common_residual,
            "best_Q_degree": degree,
            "TQ_residual": tq_residual,
            "q_coefficients_ascending": np.real_if_close(qcoef),
            "Psi_transfer_residual": psi_residual,
            "Psi_left_transfer_residual": psi_left_residual,
            "product_transfer_residual": product_residual,
            "direct_overlap": overlap,
        })
    if lax_sign == 1:
        direct_matrix = np.column_stack(direct_columns)
        psi_matrix = np.column_stack(psi_columns)
        intertwiner = direct_matrix @ np.linalg.inv(psi_matrix)
        print("universal map from determinant Psi columns to direct eigenvectors:")
        print(np.real_if_close(intertwiner, tol=1000))
        print("off-diagonal norm:", np.linalg.norm(intertwiner - np.diag(np.diag(intertwiner))))


# A state-independent inverse transform from the failed jet coordinates to the
# physical polynomial coordinates.  Direct eigenvectors are fixed by setting
# their first component to one; jet columns retain their literal residue scale.
plus_states = eigensystem(1)
physical_columns = []
jet_columns = []
for state in plus_states:
    solution = qSolve = None
    candidates = [baxter_nullvector(state[0], degree) for degree in range(3)]
    for degree, (qcoef, residual, _) in enumerate(candidates):
        if residual < 1e-9:
            solution = (degree, qcoef)
            break
    degree, qcoef = solution
    direct = state[1].copy()
    direct /= direct[0]
    physical_columns.append(direct)
    jet_columns.append(wavefunction(qcoef))

physical_matrix = np.column_stack(physical_columns)
jet_matrix = np.column_stack(jet_columns)
inverse_sov = physical_matrix @ np.linalg.inv(jet_matrix)
corrected = inverse_sov @ jet_matrix
print("\nCALIBRATED INVERSE-SOV MAP")
print(np.real_if_close(inverse_sov, tol=1000))
print("column reconstruction residual:", np.linalg.norm(corrected - physical_matrix))

# Kronecker-rank diagnostic: rank one would mean a product of independent
# one-site frame changes.  The reshuffling below tests S=A tensor B.
reshuffled = np.empty((4, 4), dtype=complex)
for i1 in range(2):
    for i2 in range(2):
        for j1 in range(2):
            for j2 in range(2):
                reshuffled[2 * i1 + j1, 2 * i2 + j2] = inverse_sov[2 * i1 + i2, 2 * j1 + j2]
print("inverse-SoV Kronecker singular values:", np.linalg.svd(reshuffled, compute_uv=False))


# Canonical SoV reconstruction test using B(u)=T_11(u).  Diagonalize the
# commuting coefficient matrices of B, identify each separated point from the
# two roots of its eigenvalue, and test X Psi = row_weight * Product Q(x_a).
def monodromy(u, lax_sign=1):
    sites = []
    for site, theta in enumerate((theta1, theta2)):
        L = [[None, None], [None, None]]
        for i in range(2):
            for j in range(2):
                L[i][j] = ((u - theta) * np.eye(4, dtype=complex) if i == j else np.zeros((4, 4), dtype=complex))
                L[i][j] += lax_sign * h * embed(E[j][i], site)
        sites.append(L)
    return [[sum((sites[0][i][k] @ sites[1][k][j] for k in range(2)),
                 np.zeros((4, 4), dtype=complex)) for j in range(2)] for i in range(2)]


def b_coefficients():
    bm = monodromy(-1.0)[0][0]
    b0 = monodromy(0.0)[0][0]
    bp = monodromy(1.0)[0][0]
    return b0, (bp - bm) / 2, (bp + bm) / 2 - b0


bcs = b_coefficients()
_, left_columns = np.linalg.eig((bcs[0] + 0.319 * bcs[1]).T)
Xrows = []
root_rows = []
for k in range(4):
    row = left_columns[:, k]
    coeff = np.array([(row @ (c @ np.ones(4))) / (row @ np.ones(4)) for c in bcs])
    # More stable eigenvalue extraction using any component of row B = b row.
    coeff = []
    for c in bcs:
        rc = row @ c
        j = np.argmax(np.abs(row))
        coeff.append(rc[j] / row[j])
    roots = np.roots([coeff[2], coeff[1], coeff[0]])
    Xrows.append(row)
    root_rows.append(roots)
X = np.vstack(Xrows)

q_solutions = []
for state in plus_states:
    for degree in range(3):
        qcoef, residual, _ = baxter_nullvector(state[0], degree)
        if residual < 1e-9:
            q_solutions.append(qcoef)
            break

print("\nSOV B-EIGENVALUE ROOTS")
print(np.real_if_close(np.array(root_rows), tol=1000))
for shift in (-h, 0.0, h):
    ratios = np.zeros((4, 4), dtype=complex)
    for row in range(4):
        roots = root_rows[row]
        for state_index, (state, qcoef) in enumerate(zip(plus_states, q_solutions)):
            q = P(qcoef)
            qprod = np.prod([lam ** ((root + shift) / h) * q(root + shift) for root in roots])
            ratios[row, state_index] = (X[row] @ state[1]) / qprod
    spreads = [np.std(ratios[row] / ratios[row, 0]) for row in range(4)]
    ratio_singular = np.linalg.svd(ratios, compute_uv=False)
    print("SoV factorization shift", shift, "row spreads", spreads,
          "ratio singular values", ratio_singular)


# Normalize the four B-eigenvector rows on the degree-zero reference state,
# then invert the SoV transform.  This uses no fitted physical vector for the
# other three states.
reference_index = next(i for i, qcoef in enumerate(q_solutions) if len(qcoef) == 1)
reference_state = plus_states[reference_index]
reference_q = P(q_solutions[reference_index])
reference_W = np.array([
    np.prod([lam ** (root / h) * reference_q(root) for root in roots])
    for roots in root_rows
])
row_scales = reference_W / (X @ reference_state[1])
X_normalized = np.diag(row_scales) @ X

inverse_sov_residuals = []
inverse_sov_overlaps = []
for state, qcoef in zip(plus_states, q_solutions):
    q = P(qcoef)
    W = np.array([
        np.prod([lam ** (root / h) * q(root) for root in roots])
        for roots in root_rows
    ])
    reconstructed = np.linalg.solve(X_normalized, W)
    eig0 = sum(state[0][k] * 0.23 ** k for k in range(3))
    eig1 = sum(state[0][k] * (-0.67) ** k for k in range(3))
    residual = max(
        normalized_residual(transfer(0.23, 1), reconstructed, eig0),
        normalized_residual(transfer(-0.67, 1), reconstructed, eig1),
    )
    overlap = abs(np.vdot(state[1], reconstructed)) / np.linalg.norm(reconstructed)
    inverse_sov_residuals.append(residual)
    inverse_sov_overlaps.append(overlap)
print("inverse-SoV physical residuals:", inverse_sov_residuals)
print("inverse-SoV direct overlaps:", inverse_sov_overlaps)


print("\nALTERNATIVE ANSATZ RESIDUALS (sign-matched Lax)")
for name, builder in [
    ("jet Q(u-h)", lambda q: wavefunction(q, -h)),
    ("jet Q(u)", lambda q: wavefunction(q, 0.0)),
    ("jet Q(u+h)", lambda q: wavefunction(q, h)),
    ("product Q(u-h)", lambda q: product_wavefunction(q, -h)),
    ("product Q(u)", lambda q: product_wavefunction(q, 0.0)),
    ("product Q(u+h)", lambda q: product_wavefunction(q, h)),
    ("site-swapped jet Q(u-h)", lambda q: wavefunction(q, -h)[[0, 2, 1, 3]]),
]:
    residuals = []
    for state, qcoef in zip(plus_states, q_solutions):
        candidate = builder(qcoef)
        eig0 = sum(state[0][k] * 0.23 ** k for k in range(3))
        eig1 = sum(state[0][k] * (-0.67) ** k for k in range(3))
        residuals.append(max(
            normalized_residual(transfer(0.23, 1), candidate, eig0),
            normalized_residual(transfer(-0.67, 1), candidate, eig1),
        ))
    print(name, residuals)


# Purely functional alternatives: apply the fundamental rational R-operator
# (delta Identity +/- h Permutation), or its inverse, to the product of the two
# one-site contour transforms.  P swaps z1 and z2, hence components 1 and 2 in
# the basis {1,z2,z1,z1 z2}.
Pswap = np.array([[1, 0, 0, 0],
                  [0, 0, 1, 0],
                  [0, 1, 0, 0],
                  [0, 0, 0, 1]], dtype=complex)
delta12 = theta1 - theta2
functional_operators = {
    "R(delta)=delta I+hP": delta12 * np.eye(4) + h * Pswap,
    "R(delta)=delta I-hP": delta12 * np.eye(4) - h * Pswap,
    "R(-delta)=-delta I+hP": -delta12 * np.eye(4) + h * Pswap,
    "R(-delta)=-delta I-hP": -delta12 * np.eye(4) - h * Pswap,
}
functional_operators |= {
    name + " inverse": np.linalg.inv(op) for name, op in list(functional_operators.items())
}
print("\nR-OPERATOR FUNCTIONAL ANSATZ RESIDUALS")
for shift in (-h, 0.0, h):
    for name, op in functional_operators.items():
        residuals = []
        for state, qcoef in zip(plus_states, q_solutions):
            candidate = op @ product_wavefunction(qcoef, shift)
            eig0 = sum(state[0][k] * 0.23 ** k for k in range(3))
            eig1 = sum(state[0][k] * (-0.67) ** k for k in range(3))
            residuals.append(max(
                normalized_residual(transfer(0.23, 1), candidate, eig0),
                normalized_residual(transfer(-0.67, 1), candidate, eig1),
            ))
        if max(residuals) < 1e-8:
            print("WORKS", shift, name, residuals)
        elif max(residuals) < 0.25:
            print("near", shift, name, residuals)


# Closed functional factorizing twist.  No numerical SoV eigenvector matrix:
# B(u)=T11(u) is triangular in {1,z2,z1,z1z2}, and its sole nontrivial left
# eigenvector is <10|-(h/(theta1-theta2))<01|.
delta = theta1 - theta2
X_closed = np.array([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, -h / delta, 1, 0],
    [0, 0, 0, 1],
], dtype=complex)
root_closed = np.array([
    [theta1 - h, theta2 - h],
    [theta1 - h, theta2],
    [theta1, theta2 - h],
    [theta1, theta2],
], dtype=complex)

# The reference is the product of the known one-site functions
# phi_0(z)=z-1/lambda.  Overall scalar factors cancel from the construction.
vacuum_polynomial = np.array([1 / lam**2, -1 / lam, -1 / lam, 1], dtype=complex)
vacuum_W = np.array([np.prod([lam ** (x / h) for x in roots]) for roots in root_closed])
closed_row_weights = (X_closed @ vacuum_polynomial) / vacuum_W


def closed_functional_wavefunction(qcoef):
    q = P(qcoef)
    W = np.array([
        np.prod([lam ** (x / h) * q(x) for x in roots])
        for roots in root_closed
    ])
    separated = closed_row_weights * W
    return np.array([
        separated[0],
        separated[1],
        separated[2] + (h / delta) * separated[1],
        separated[3],
    ])


closed_residuals = []
for state, qcoef in zip(plus_states, q_solutions):
    candidate = closed_functional_wavefunction(qcoef)
    eig0 = sum(state[0][k] * 0.23 ** k for k in range(3))
    eig1 = sum(state[0][k] * (-0.67) ** k for k in range(3))
    closed_residuals.append(max(
        normalized_residual(transfer(0.23, 1), candidate, eig0),
        normalized_residual(transfer(-0.67, 1), candidate, eig1),
    ))
print("closed functional factorizing-twist residuals:", closed_residuals)
print("closed functional row weights:", closed_row_weights)


# User-supplied diagonal datum: the 2x2 scalar-product determinant with
# QA=z^(u/h).  For two spin-1/2 sites it fixes A, B1+B2, C.  Solve the transfer
# equation for the one missing antisymmetric coefficient.
def lift_from_diagonal(qcoef, state):
    polarized = wavefunction(qcoef, 0.0)  # raw Q(u), as in QA*QB
    A = polarized[0]
    Bsum = polarized[1] + polarized[2]
    C = polarized[3]
    # Basis {1,z2,z1,z1z2}: c = base + x direction, with
    # c_z2=x and c_z1=Bsum-x.
    base = np.array([A, 0, Bsum, C], dtype=complex)
    direction = np.array([0, 1, -1, 0], dtype=complex)
    lhs_blocks = []
    rhs_blocks = []
    for spectral in (0.23, -0.67):
        eig = sum(state[0][k] * spectral ** k for k in range(3))
        operator = transfer(spectral, 1) - eig * np.eye(4)
        lhs_blocks.append(operator @ direction)
        rhs_blocks.append(-(operator @ base))
    lhs = np.concatenate(lhs_blocks)
    rhs = np.concatenate(rhs_blocks)
    x = np.vdot(lhs, rhs) / np.vdot(lhs, lhs)
    lifted = base + x * direction
    residual = max(
        normalized_residual(transfer(0.23, 1), lifted,
                            sum(state[0][k] * 0.23 ** k for k in range(3))),
        normalized_residual(transfer(-0.67, 1), lifted,
                            sum(state[0][k] * (-0.67) ** k for k in range(3))),
    )
    return lifted, x, residual, (A, Bsum, C)


diagonal_lift_results = [
    lift_from_diagonal(qcoef, state)
    for state, qcoef in zip(plus_states, q_solutions)
]
print("diagonal-determinant lift residuals:", [item[2] for item in diagonal_lift_results])
print("solved z2 linear coefficients:", [item[1] for item in diagonal_lift_results])


# Verify the symbolic lift formula itself using the exact diagonal restriction
# of each direct eigenvector.  This separates the PDE lift from the unresolved
# convention mismatch in our implementation of the bracket determinant.
def differential_diagonal_lift(diagonal_coefficients, state, lax_sign=1):
    A, Bsum, C = diagonal_coefficients
    I0_state = state[0][0]
    k = lax_sign * h
    x_z2 = (
        A * (-I0_state + chi * (k**2 - k * theta1 - k * theta2 + theta1 * theta2))
        + Bsum * k * (k - theta2)
    ) / (k * (k + theta1 - theta2))
    return np.array([A, x_z2, Bsum - x_z2, C], dtype=complex)


direct_diagonal_lift_residuals = []
for state in plus_states:
    v = state[1]
    lifted = differential_diagonal_lift((v[0], v[1] + v[2], v[3]), state)
    eig0 = sum(state[0][k] * 0.23 ** k for k in range(3))
    eig1 = sum(state[0][k] * (-0.67) ** k for k in range(3))
    direct_diagonal_lift_residuals.append(max(
        normalized_residual(transfer(0.23, 1), lifted, eig0),
        normalized_residual(transfer(-0.67, 1), lifted, eig1),
    ))
print("differential lift from exact diagonal residuals:", direct_diagonal_lift_residuals)

literal_diagonal_lift_residuals = []
for state in eigensystem(-1):
    v = state[1]
    lifted = differential_diagonal_lift(
        (v[0], v[1] + v[2], v[3]), state, lax_sign=-1
    )
    eig0 = sum(state[0][k] * 0.23 ** k for k in range(3))
    eig1 = sum(state[0][k] * (-0.67) ** k for k in range(3))
    literal_diagonal_lift_residuals.append(max(
        normalized_residual(transfer(0.23, -1), lifted, eig0),
        normalized_residual(transfer(-0.67, -1), lifted, eig1),
    ))
print("literal L=u-hE^t diagonal lift residuals:", literal_diagonal_lift_residuals)


# The closed factorizing-twist wavefunction is itself one 2x2 determinant.
# Put a=Qtw(theta1-h), b=Qtw(theta1), c=Qtw(theta2-h), d=Qtw(theta2)
# and r=h/(theta1-theta2).  Its determinant is
# det [[a,b],[(1-r)z1 c-z1 z2 d, c-(z2+r z1)d]].
def full_shifted_q_determinant(qcoef):
    q = P(qcoef)
    qtw = lambda x: lam ** (x / h) * q(x)
    a = qtw(theta1 - h)
    b = qtw(theta1)
    c = qtw(theta2 - h)
    d = qtw(theta2)
    r = h / (theta1 - theta2)
    common_normalization = closed_row_weights[0]
    return common_normalization * np.array([
        a * c,
        -a * d,
        (r - 1) * b * c - r * a * d,
        b * d,
    ], dtype=complex)


full_determinant_differences = []
full_determinant_residuals = []
for state, qcoef in zip(plus_states, q_solutions):
    determinant_vector = full_shifted_q_determinant(qcoef)
    closed_vector = closed_functional_wavefunction(qcoef)
    full_determinant_differences.append(
        np.linalg.norm(determinant_vector - closed_vector)
        / np.linalg.norm(closed_vector)
    )
    eig0 = sum(state[0][k] * 0.23 ** k for k in range(3))
    eig1 = sum(state[0][k] * (-0.67) ** k for k in range(3))
    full_determinant_residuals.append(max(
        normalized_residual(transfer(0.23, 1), determinant_vector, eig0),
        normalized_residual(transfer(-0.67, 1), determinant_vector, eig1),
    ))
print("full shifted-Q determinant vs closed formula:", full_determinant_differences)
print("full shifted-Q determinant transfer residuals:", full_determinant_residuals)


# Operator-level proof.  On each SoV corner x_i in {theta_i-h,theta_i},
# interpolate the quadratic t(u) from t(x1),t(x2), while the endpoint Baxter
# operator supplies t(x_i).  The triangular factorizing twist M intertwines
# this finite-difference operator with the physical transfer matrix.
def sov_transfer(spectral):
    corners = [
        (theta1 - h, theta2 - h),
        (theta1 - h, theta2),
        (theta1, theta2 - h),
        (theta1, theta2),
    ]
    result = np.zeros((4, 4), dtype=complex)
    index = {(round(x1, 12), round(x2, 12)): j
             for j, (x1, x2) in enumerate(corners)}
    for source, (x1, x2) in enumerate(corners):
        result[source, source] += chi * (spectral - x1) * (spectral - x2)

        if abs(x1 - (theta1 - h)) < 1e-10:
            target_x1 = theta1
            baxter1 = qtheta(x1)
        else:
            target_x1 = theta1 - h
            baxter1 = qtheta(x1 + h)
        target = index[(round(target_x1, 12), round(x2, 12))]
        result[source, target] += (spectral - x2) / (x1 - x2) * baxter1

        if abs(x2 - (theta2 - h)) < 1e-10:
            target_x2 = theta2
            baxter2 = qtheta(x2)
        else:
            target_x2 = theta2 - h
            baxter2 = qtheta(x2 + h)
        target = index[(round(x1, 12), round(target_x2, 12))]
        result[source, target] += (x1 - spectral) / (x1 - x2) * baxter2
    return result


r_twist = h / (theta1 - theta2)
factorizing_twist = np.array([
    [1, 0, 0, 0],
    [0, -1, 0, 0],
    [0, -r_twist, r_twist - 1, 0],
    [0, 0, 0, 1],
], dtype=complex)
intertwining_residuals = []
for spectral in (0.23, -0.67, 1.11):
    lhs = transfer(spectral, 1) @ factorizing_twist
    rhs = factorizing_twist @ sov_transfer(spectral)
    intertwining_residuals.append(
        np.linalg.norm(lhs - rhs) / max(np.linalg.norm(lhs), 1e-30)
    )
print("factorizing-twist operator intertwining residuals:", intertwining_residuals)


# Diagonal-preserving differential twist acting on the shifted jet determinant.
def diagonal_preserving_twist(qcoef):
    jet = wavefunction(qcoef, shift=-h)
    delta = theta2 - theta1
    twist = h / (delta - h)
    z2_linear = jet[1]
    return jet + twist * z2_linear * np.array([0, 1, -1, 0])


diagonal_twist_residuals = []
for state, qcoef in zip(plus_states, q_solutions):
    candidate = diagonal_preserving_twist(qcoef)
    eig0 = sum(state[0][k] * 0.23 ** k for k in range(3))
    eig1 = sum(state[0][k] * (-0.67) ** k for k in range(3))
    diagonal_twist_residuals.append(max(
        normalized_residual(transfer(0.23, 1), candidate, eig0),
        normalized_residual(transfer(-0.67, 1), candidate, eig1),
    ))
print("diagonal-preserving shifted-jet twist residuals:", diagonal_twist_residuals)

shifted_vs_unshifted_diagonal = []
for qcoef in q_solutions:
    shifted = wavefunction(qcoef, shift=-h)
    unshifted = wavefunction(qcoef, shift=0.0)
    gs = np.array([shifted[0], shifted[1] + shifted[2], shifted[3]])
    gu = np.array([unshifted[0], unshifted[1] + unshifted[2], unshifted[3]])
    scale = np.vdot(gu, gs) / np.vdot(gu, gu)
    shifted_vs_unshifted_diagonal.append(
        np.linalg.norm(gs - scale * gu) / np.linalg.norm(gs)
    )
print("shifted vs unshifted diagonal proportionality residuals:",
      shifted_vs_unshifted_diagonal)
