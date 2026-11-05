# TRL assessment – QuantumBench

**Claimed level: TRL 4** (simulation-based validation). Not TRL 5: no hardware, no noise model, one 4-qubit molecule.

H2/STO-3G, 14 non-identity Pauli terms. Total shots for one VQE run (3 optimisation evaluations + 1 final measurement) by shots per term — 100: 5,600, 1000: 56,000, 10000: 560,000, 100000: 5,600,000.

Reported-energy RMSE in mHa (share of reports within 1.6 mHa) by shots per term; classical times in ms (HF / FCI):

| Bond (Å) | HF error (mHa) | HF / FCI time (ms) | 100 | 1 000 | 10 000 | 100 000 |
|---|---|---|---|---|---|---|
| 0.5 | 12.2 | 105 / 1 | 12.01 (0.12) | 3.48 (0.32) | 1.25 (0.79) | 0.37 (1.00) |
| 0.735 | 20.3 | 16 / 1 | 12.27 (0.08) | 3.64 (0.33) | 1.29 (0.80) | 0.44 (1.00) |
| 1.0 | 35.0 | 16 / 1 | 12.28 (0.09) | 4.00 (0.34) | 1.30 (0.81) | 0.41 (1.00) |
| 1.5 | 87.3 | 16 / 1 | 12.49 (0.08) | 3.96 (0.27) | 1.27 (0.82) | 0.37 (1.00) |
| 2.0 | 164.8 | 15 / 1 | 10.86 (0.10) | 3.34 (0.36) | 0.98 (0.90) | 0.31 (1.00) |
| 2.5 | 233.1 | 12 / 1 | 8.28 (0.15) | 2.84 (0.44) | 0.92 (0.90) | 0.26 (1.00) |

## Findings and caveats
- The one-parameter ansatz contains the exact ground state, so the returned *state* is within 1 mHa of FCI already at 1 000 shots per term. The *reported energy*, however, is limited by measurement noise: a report within chemical accuracy needs about 10^4 shots per term for roughly 80–90 % of runs (about 5.6 x 10^5 shots in total) and 10^5 per term (about 5.6 x 10^6 shots) for all runs.
- Classical Hartree–Fock and exact FCI take milliseconds on this system and FCI is exact; **the quantum algorithm gives no advantage here** and the comparison only quantifies its measurement cost.
- Hartree–Fock misses 12 mHa at 0.5 Å and 233 mHa at 2.5 Å (static correlation), so the variational result is much better than Hartree–Fock — but any exact classical method for two electrons is also.
- No grouping of commuting Pauli terms and no error mitigation: grouping would reduce the shot count by a modest factor. No gate or readout noise was simulated.
- The fully simulated 'quantum' side costs more classical computation than FCI itself.
