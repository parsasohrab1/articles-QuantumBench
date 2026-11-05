# Software Requirements Specification – QuantumBench

- SR-1 `hamiltonian.molecular_hamiltonian`: 16x16 qubit Hamiltonian from PySCF one- and two-electron integrals by Jordan–Wigner mapping (spin orbitals interleaved); `pauli_decomposition`: Pauli coefficients. Unit tests check fermionic anticommutation, Hermiticity and the Pauli round trip.
- SR-2 `vqe.ansatz_state`: cos(theta)|HF> + sin(theta)|double excitation> (exact ground-state manifold of H2/STO-3G).
- SR-3 `vqe.noisy_energy`: energy from binomial sampling of every non-identity Pauli term separately with ``shots`` measurements each (no grouping).
- SR-4 `vqe.rotosolve`: three noisy energy evaluations determine a + b cos 2theta + c sin 2theta and its minimum.
- SR-5 `research/run_study.py`: bond lengths 0.5–2.5 Å, shots per term 1e2–1e5, 200 repetitions; reports the root-mean-square error of the *reported* energy (a final measurement at the returned parameter) and the share of reports within 1.6 mHa.
