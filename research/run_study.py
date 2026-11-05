"""H2 (STO-3G): classical (Hartree-Fock, FCI) versus VQE with shot noise, bond-length scan.

The molecular integrals come from PySCF; the qubit Hamiltonian is built here (Jordan-Wigner, 4 qubits).
The 'quantum' side is a noiseless statevector simulation plus binomial shot noise: no hardware noise, no device.
usage: python research/run_study.py
"""

import json
import time
from pathlib import Path

import numpy as np
from pyscf import fci, gto, scf

from quantumbench.hamiltonian import molecular_hamiltonian, pauli_decomposition
from quantumbench.vqe import ansatz_state, noisy_energy, rotosolve

BONDS = (0.5, 0.735, 1.0, 1.5, 2.0, 2.5)
SHOTS = (100, 1000, 10000, 100000)
REPS = 200
CHEM_ACC = 1.6e-3  # Hartree


def one_bond(r: float, rng):
    mol = gto.M(atom=f"H 0 0 0; H 0 0 {r}", basis="sto-3g", verbose=0)
    t0 = time.perf_counter()
    mf = scf.RHF(mol).run()
    t_hf = time.perf_counter() - t0
    t0 = time.perf_counter()
    cis = fci.FCI(mf).run()
    t_fci = time.perf_counter() - t0
    c = mf.mo_coeff
    h1 = c.T @ mf.get_hcore() @ c
    from pyscf import ao2mo

    eri = ao2mo.restore(1, ao2mo.kernel(mol, c), c.shape[1])
    H = molecular_hamiltonian(h1, eri, mol.energy_nuc())
    coeffs = pauli_decomposition(H)
    w = np.linalg.eigvalsh(H)
    # exact ground energy of the 2-electron sector equals FCI; the lowest eigenvalue of the full Fock space may be lower
    res = {
        "r_angstrom": r,
        "e_hf": float(mf.e_tot),
        "e_fci": float(cis.e_tot),
        "hf_error_mha": 1e3 * (mf.e_tot - cis.e_tot),
        "n_pauli_terms": len(coeffs) - 1,
        "t_hf_s": t_hf,
        "t_fci_s": t_fci,
        "fock_space_min_minus_fci": float(w[0] - cis.e_tot),
        "shots": {},
    }
    for s in SHOTS:
        errs, cnt, est_errs, est_cnt = [], 0, [], 0
        for _ in range(REPS):
            theta, _est = rotosolve(coeffs, s, rng)
            psi = ansatz_state(theta)
            e = float(
                np.real(psi.conj() @ H @ psi)
            )  # exact energy of the state the noisy optimiser returned
            errs.append(e - cis.e_tot)
            cnt += abs(e - cis.e_tot) < CHEM_ACC
            est = noisy_energy(psi, coeffs, s, rng)  # what the experimenter would actually report
            est_errs.append(est - cis.e_tot)
            est_cnt += abs(est - cis.e_tot) < CHEM_ACC
        res["shots"][str(s)] = {
            "mean_energy_error_mha": 1e3 * float(np.mean(errs)),
            "p95_energy_error_mha": 1e3 * float(np.quantile(np.abs(errs), 0.95)),
            "share_chemical_accuracy": cnt / REPS,
            "reported_energy_rmse_mha": 1e3 * float(np.sqrt(np.mean(np.square(est_errs)))),
            "reported_share_chemical_accuracy": est_cnt / REPS,
            "total_shots": 4 * s * (len(coeffs) - 1),
        }
    return res


def main():
    rng = np.random.default_rng(0)
    out = {
        "note": "H2/STO-3G; classical = PySCF HF and FCI; quantum = noiseless statevector + shot noise, one-parameter "
        "UCC-like ansatz, Rotosolve (3 evaluations); 200 repetitions per setting; chemical accuracy 1.6 mHa",
        "bonds": [one_bond(r, rng) for r in BONDS],
    }
    Path(__file__).with_name("results.json").write_text(json.dumps(out, indent=1))
    for b in out["bonds"]:
        print(
            b["r_angstrom"],
            round(b["hf_error_mha"], 2),
            b["n_pauli_terms"],
            round(b["fock_space_min_minus_fci"], 4),
            {
                s: (
                    round(v["mean_energy_error_mha"], 3),
                    round(v["reported_energy_rmse_mha"], 2),
                    float(v["reported_share_chemical_accuracy"]),
                )
                for s, v in b["shots"].items()
            },
        )


if __name__ == "__main__":
    main()
