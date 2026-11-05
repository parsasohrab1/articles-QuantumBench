from __future__ import annotations

import itertools

import numpy as np

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]])
Z = np.diag([1.0, -1.0]).astype(complex)
LOWER = np.array([[0, 1], [0, 0]], dtype=complex)  # |0><1| : annihilates an occupied (|1>) orbital


def kron(ops):
    out = np.array([[1.0 + 0j]])
    for o in ops:
        out = np.kron(out, o)
    return out


def annihilation(p: int, n: int) -> np.ndarray:
    """Jordan-Wigner annihilation operator of spin-orbital p on n qubits (qubit 0 is the leftmost factor)."""
    return kron([Z] * p + [LOWER] + [I2] * (n - p - 1))


def molecular_hamiltonian(h1: np.ndarray, eri: np.ndarray, e_nuc: float) -> np.ndarray:
    """Dense qubit Hamiltonian from spatial-orbital integrals (chemist notation eri[p,q,r,s] = (pq|rs))."""
    norb = h1.shape[0]
    n = 2 * norb
    a = [annihilation(p, n) for p in range(n)]
    ad = [x.conj().T for x in a]
    H = np.zeros((2**n, 2**n), dtype=complex)
    for p, q in itertools.product(range(n), repeat=2):
        if p % 2 == q % 2:
            H += h1[p // 2, q // 2] * ad[p] @ a[q]
    for p, q, r, s in itertools.product(range(n), repeat=4):
        if p % 2 == s % 2 and q % 2 == r % 2:
            H += 0.5 * eri[p // 2, s // 2, q // 2, r // 2] * ad[p] @ ad[q] @ a[r] @ a[s]
    return H + e_nuc * np.eye(2**n)


def pauli_decomposition(H: np.ndarray, tol: float = 1e-10) -> dict[str, float]:
    n = int(np.log2(H.shape[0]))
    mats = {"I": I2, "X": X, "Y": Y, "Z": Z}
    out = {}
    for labels in itertools.product("IXYZ", repeat=n):
        P = kron([mats[c] for c in labels])
        c = np.trace(P.conj().T @ H).real / H.shape[0]
        if abs(c) > tol:
            out["".join(labels)] = float(c)
    return out


def pauli_matrix(label: str) -> np.ndarray:
    mats = {"I": I2, "X": X, "Y": Y, "Z": Z}
    return kron([mats[c] for c in label])
