import numpy as np

from quantumbench.hamiltonian import annihilation, molecular_hamiltonian, pauli_decomposition
from quantumbench.vqe import ansatz_state, noisy_energy, rotosolve


def test_fermionic_anticommutation():
    n = 4
    a = [annihilation(p, n) for p in range(n)]
    for p in range(n):
        for q in range(n):
            ac = a[p] @ a[q].conj().T + a[q].conj().T @ a[p]
            assert np.allclose(ac, np.eye(16) if p == q else 0)


def _toy():
    h1 = np.array([[-1.25, 0.0], [0.0, -0.47]])
    eri = np.zeros((2, 2, 2, 2))
    eri[0, 0, 0, 0] = 0.67
    eri[1, 1, 1, 1] = 0.70
    eri[0, 0, 1, 1] = eri[1, 1, 0, 0] = 0.66
    eri[0, 1, 0, 1] = eri[1, 0, 1, 0] = eri[0, 1, 1, 0] = eri[1, 0, 0, 1] = 0.18
    return molecular_hamiltonian(h1, eri, 0.7)


def test_hamiltonian_hermitian_and_pauli_roundtrip():
    H = _toy()
    assert np.allclose(H, H.conj().T)
    coeffs = pauli_decomposition(H)
    from quantumbench.hamiltonian import pauli_matrix

    R = sum(c * pauli_matrix(k) for k, c in coeffs.items())
    assert np.allclose(R, H, atol=1e-8)


def test_noiseless_limit_and_rotosolve_finds_minimum():
    H = _toy()
    coeffs = pauli_decomposition(H)
    rng = np.random.default_rng(0)
    theta, _ = rotosolve(coeffs, 10**7, rng)
    psi = ansatz_state(theta)
    e = float(np.real(psi.conj() @ H @ psi))
    grid = [
        float(np.real(ansatz_state(t).conj() @ H @ ansatz_state(t)))
        for t in np.linspace(0, np.pi, 400)
    ]
    assert e <= min(grid) + 1e-3
    assert abs(noisy_energy(psi, coeffs, 10**7, rng) - e) < 5e-3
