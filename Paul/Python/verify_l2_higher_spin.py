"""Numerical experiments for two-site symmetric gl2 representations."""

import math

import numpy as np
from numpy.polynomial import Polynomial as P


H = 0.73
THETA1 = -0.41
THETA2 = 0.86
LAM = 1.67
CHI = LAM + 1 / LAM


def qtheta(x):
    return (x - THETA1) * (x - THETA2)


def generators(boxes):
    dim = boxes + 1
    e11 = np.diag([boxes - k for k in range(dim)]).astype(complex)
    e22 = np.diag([k for k in range(dim)]).astype(complex)
    e12 = np.zeros((dim, dim), dtype=complex)
    e21 = np.zeros((dim, dim), dtype=complex)
    for k in range(dim):
        if k > 0:
            e12[k - 1, k] = -k
        if k < boxes:
            e21[k + 1, k] = k - boxes
    return [[e11, e12], [e21, e22]]


def transfer(boxes, spectral):
    dim = boxes + 1
    ident = np.eye(dim, dtype=complex)
    ident2 = np.eye(dim * dim, dtype=complex)
    egrid = generators(boxes)

    def embed(matrix, site):
        return np.kron(matrix, ident) if site == 0 else np.kron(ident, matrix)

    def lax(theta, site):
        return [[
            (spectral - theta) * (ident2 if i == j else np.zeros_like(ident2))
            + H * embed(egrid[j][i], site)
            for j in range(2)
        ] for i in range(2)]

    left = lax(THETA1, 0)
    right = lax(THETA2, 1)
    monodromy = [[
        sum((left[i][k] @ right[k][j] for k in range(2)), np.zeros_like(ident2))
        for j in range(2)
    ] for i in range(2)]
    return CHI * monodromy[0][0] + monodromy[0][1] - monodromy[1][0]


def transfer_coefficients(boxes):
    samples = (-0.83, 0.17, 1.21)
    vandermonde = np.array([[1, x, x * x] for x in samples], dtype=complex)
    values = np.stack([transfer(boxes, x) for x in samples])
    coefficients = np.tensordot(np.linalg.inv(vandermonde), values, axes=(1, 0))
    return coefficients


def common_states(boxes):
    coefficients = transfer_coefficients(boxes)
    _, vectors = np.linalg.eig(coefficients[0] + 0.371 * coefficients[1])
    states = []
    for index in range(vectors.shape[1]):
        vector = vectors[:, index]
        vector /= np.linalg.norm(vector)
        tcoeff = np.array([
            np.vdot(vector, matrix @ vector) / np.vdot(vector, vector)
            for matrix in coefficients
        ])
        residual = max(np.linalg.norm(matrix @ vector - eig * vector)
                       for matrix, eig in zip(coefficients, tcoeff))
        states.append((tcoeff, vector, residual))
    return states


def baxter_matrix(tcoeff, degree, boxes):
    tpoly = P(tcoeff)
    columns = []
    qtheta_poly = P([THETA1 * THETA2, -(THETA1 + THETA2), 1])
    for power in range(degree + 1):
        basis = P([0] * power + [1])
        first = qtheta_poly(P([boxes * H, 1])) * basis(P([-H, 1])) / LAM
        second = qtheta_poly * basis(P([H, 1])) * LAM
        residual = first - tpoly * basis + second
        columns.append(np.pad(residual.coef, (0, max(0, degree + 3 - len(residual.coef))))
                       [:degree + 3])
    return np.stack(columns, axis=1)


def solve_q(tcoeff, boxes):
    candidates = []
    for degree in range(2 * boxes + 1):
        matrix = baxter_matrix(tcoeff, degree, boxes)
        _, singular, vh = np.linalg.svd(matrix)
        vector = vh.conj().T[:, -1]
        residual = np.linalg.norm(matrix @ vector)
        if abs(vector[-1]) > 1e-12:
            vector /= vector[-1]
        candidates.append((degree, vector, residual, singular[-1]))
    good = [item for item in candidates if item[2] < 1e-8]
    return good[0] if good else min(candidates, key=lambda item: item[2])


def localized_residue_weight(boxes, site, k):
    theta_i, theta_j = ((THETA1, THETA2) if site == 0 else (THETA2, THETA1))
    pole = theta_i + k * H
    rho = 1 - np.exp(2j * np.pi * (pole - theta_j) / H)
    derivative_factor = theta_i - theta_j
    denominator = derivative_factor
    for m in range(boxes + 1):
        if m == k:
            continue
        denominator *= qtheta(pole - m * H)
    return rho / denominator


def qfull(qcoeff, x):
    return LAM ** (x / H) * P(qcoeff)(x)


def moment_vector(qcoeff, boxes):
    dim = boxes + 1
    data = []
    for site, theta_i in enumerate((THETA1, THETA2)):
        row = []
        for k in range(dim):
            weight = localized_residue_weight(boxes, site, k)
            qvalue = qfull(qcoeff, theta_i + (k - boxes) * H)
            row.append(weight * qvalue)
        data.append(row)
    vector = np.zeros(dim * dim, dtype=complex)
    for k in range(dim):
        for ell in range(dim):
            x1 = THETA1 + k * H
            x2 = THETA2 + ell * H
            vector[k * dim + ell] = (x2 - x1) * data[0][k] * data[1][ell]
    return vector


def diagonal_coefficients(vector, boxes):
    dim = boxes + 1
    diagonal = np.zeros(2 * boxes + 1, dtype=complex)
    for k in range(dim):
        for ell in range(dim):
            diagonal[k + ell] += vector[k * dim + ell]
    return diagonal


def normalized_residual(matrix, vector, eigenvalue):
    return np.linalg.norm(matrix @ vector - eigenvalue * vector) / np.linalg.norm(vector)


def rising(x, count):
    value = 1.0 + 0j
    for offset in range(count):
        value *= x + offset
    return value


def universal_twist(boxes):
    """Diagonal-preserving finite factorizing twist in the monomial basis."""
    dim = boxes + 1
    d = (THETA1 - THETA2) / H
    matrix = np.zeros((dim * dim, dim * dim), dtype=complex)
    for a in range(dim):
        for b in range(dim):
            column = a * dim + b
            diagonal = rising(d + a - b + 1, b) / rising(d + boxes - b + 1, b)
            for j in range(min(boxes - a, b) + 1):
                falling_a = math.factorial(boxes - a) / math.factorial(boxes - a - j)
                falling_b = math.factorial(b) / math.factorial(b - j)
                coefficient = (diagonal * falling_a * falling_b
                               / (math.factorial(j) * rising(d + a - b + 1, j)))
                row = (a + j) * dim + (b - j)
                matrix[row, column] = coefficient
    return matrix


def experiment(boxes):
    states = common_states(boxes)
    solutions = [solve_q(state[0], boxes) for state in states]
    physical_columns = []
    moment_columns = []
    diagonal_residuals = []
    tq_residuals = []
    transfer_residuals = []

    for state, solution in zip(states, solutions):
        degree, qcoeff, tq_residual, _ = solution
        tq_residuals.append((degree, tq_residual))
        moment = moment_vector(qcoeff, boxes)
        direct = state[1]
        direct_diagonal = diagonal_coefficients(direct, boxes)
        moment_diagonal = diagonal_coefficients(moment, boxes)
        scale = np.vdot(direct_diagonal, moment_diagonal) / np.vdot(direct_diagonal, direct_diagonal)
        scaled_direct = scale * direct
        diagonal_residuals.append(
            np.linalg.norm(diagonal_coefficients(scaled_direct, boxes) - moment_diagonal)
            / max(np.linalg.norm(moment_diagonal), 1e-30)
        )
        physical_columns.append(scaled_direct)
        moment_columns.append(moment)
        eig0 = sum(state[0][k] * 0.23 ** k for k in range(3))
        transfer_residuals.append(normalized_residual(transfer(boxes, 0.23), scaled_direct, eig0))

    physical = np.stack(physical_columns, axis=1)
    moments = np.stack(moment_columns, axis=1)
    twist = physical @ np.linalg.inv(moments)
    twist_formula = universal_twist(boxes)
    intertwining = []
    for spectral in (0.23, -0.67, 1.11):
        # The conjugated operator is the finite-difference/SoV transfer in moment coordinates.
        separated = np.linalg.solve(twist, transfer(boxes, spectral) @ twist)
        intertwining.append(np.linalg.norm(transfer(boxes, spectral) @ twist - twist @ separated))

    return {
        "boxes": boxes,
        "state_common_residuals": [state[2] for state in states],
        "tq": tq_residuals,
        "diagonal_residuals": diagonal_residuals,
        "transfer_residuals": transfer_residuals,
        "twist": twist,
        "twist_formula_residual": np.linalg.norm(twist - twist_formula),
        "twist_singular_values": np.linalg.svd(twist, compute_uv=False),
        "intertwining_residuals": intertwining,
    }


if __name__ == "__main__":
    for boxes in (1, 2, 3):
        result = experiment(boxes)
        print("\nS =", boxes)
        print("TQ degree/residual:", result["tq"])
        print("diagonal proportionality residuals:", result["diagonal_residuals"])
        print("physical transfer residuals:", result["transfer_residuals"])
        print("twist singular values:", result["twist_singular_values"])
        print("universal twist residual:", result["twist_formula_residual"])
        print("twist matrix:\n", np.round(result["twist"], 8))
