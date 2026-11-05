from __future__ import annotations

import numpy as np

from .hamiltonian import pauli_matrix


def basis_state(bits: str) -> np.ndarray:
    v = np.zeros(2 ** len(bits), dtype=complex)
    v[int(bits, 2)] = 1.0
    return v


def ansatz_state(theta: float, hf: str = "1100", doubly: str = "0011") -> np.ndarray:
    """cos(theta)|HF> + sin(theta)|double excitation>: a one-parameter UCC-like state (real)."""
    return np.cos(theta) * basis_state(hf) + np.sin(theta) * basis_state(doubly)


def term_expectations(psi: np.ndarray, labels: list[str]) -> dict[str, float]:
    return {lab: float(np.real(psi.conj() @ pauli_matrix(lab) @ psi)) for lab in labels}


def noisy_energy(
    psi: np.ndarray, coeffs: dict[str, float], shots: int, rng: np.random.Generator
) -> float:
    """Energy estimate with ``shots`` measurements per non-identity Pauli term (no grouping)."""
    e = coeffs.get("I" * int(np.log2(len(psi))), 0.0)
    ex = term_expectations(psi, [k for k in coeffs if set(k) != {"I"}])
    for lab, v in ex.items():
        p = (1 + v) / 2
        k = rng.binomial(shots, min(max(p, 0.0), 1.0))
        e += coeffs[lab] * (2 * k / shots - 1)
    return float(e)


def rotosolve(
    coeffs: dict[str, float], shots: int, rng: np.random.Generator, hf="1100", doubly="0011"
):
    """Minimise E(theta) = a + b cos(2 theta) + c sin(2 theta) from three noisy evaluations.

    Returns (theta, estimated energy at theta, exact energy of the returned state is computed by the caller).
    """
    th = np.array([0.0, np.pi / 3, 2 * np.pi / 3])
    E = np.array([noisy_energy(ansatz_state(t, hf, doubly), coeffs, shots, rng) for t in th])
    A = np.column_stack([np.ones(3), np.cos(2 * th), np.sin(2 * th)])
    a, b, c = np.linalg.solve(A, E)
    theta = 0.5 * np.arctan2(-c, -b)
    return float(theta), float(a - np.hypot(b, c))
