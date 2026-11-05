# articles-QuantumBench

Honest benchmark of a variational quantum eigensolver (VQE) against classical methods (Hartree–Fock, exact FCI) for the
hydrogen molecule H2 (STO-3G, 4 qubits), along a bond-length scan, with shot noise. Everything is a **noiseless statevector
simulation plus binomial shot noise** — no hardware, no device noise. TRL 4 (simulation-based).

```
pip install -e ".[dev]" pyscf && python -m pytest
python research/run_study.py   # writes research/results.json
```
