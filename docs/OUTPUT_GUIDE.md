# Output Guide — What the Pipeline Prints and Saves

Explanation of all outputs from `main.py`, `benchmark.py`, and the JSON result files.

---

## 1. Pipeline Output (`main.py`)

When you run:

```cmd
python main.py --pdb 2I9M --fragment 5 --qubits 8 --electrons 4 --orbitals 4
```

you get a **RESULTS SUMMARY** like this:

### Console Output Explained

| Term | Plain English | What it means |
|------|---------------|---------------|
| **Fragment** | Which piece of the protein | e.g. `2I9M_res5` = Chignolin, residue 5 centered |
| **QM atoms** | Quantum-mechanical region | Atoms treated with quantum chemistry |
| **MM atoms** | Classical environment | Surrounding atoms as point charges |
| **Qubits used** | Actual qubit count | After parity reduction (often 6 for 4e4o) |
| **HF energy** | Classical baseline | Hartree–Fock energy (no electron correlation) |
| **VQE energy** | Quantum result | Energy from the trained quantum circuit |
| **QSE energy** | Refined quantum result | VQE + quantum subspace expansion correction |
| **Mitigated energy** | Final answer | QSE + noise mitigation (ZNE) |
| **Free energy G** | Thermodynamic stability | ΔG at 310 K; negative = stable, positive = unstable |
| **Correlation energy** | VQE − HF | Amount of correlation recovered by the quantum part |
| **Total runtime** | Wall-clock time | Full pipeline execution time |

### Example Interpretation

```
HF energy:        -7.86224600 Ha   (classical reference)
VQE energy:       +0.63213921 Ha   (converged, 1 params, 2 iters)
QSE energy:       -0.00000000 Ha   (correction: -0.632139 Ha)
Mitigated energy: +0.90291931 Ha   (noise correction: +0.902919 Ha)
Free energy G:    +0.902919 Ha = 566.59 kcal/mol | UNSTABLE
```

- **HF** is the classical starting point.
- **VQE** is what the quantum circuit learned.
- **QSE** corrects for truncation of the excited-state space.
- **Mitigated** adds ZNE (noise extrapolation).
- **G** is the thermodynamic free energy; **UNSTABLE** means positive G.

> **Note:** In mock mode (no PySCF), energies are based on precomputed integrals and may differ from real quantum chemistry. The pipeline flow and terms are the same.

---

## 2. JSON Result File

Each run saves a JSON file to `outputs/energies/`:

```
outputs/energies/result_2I9M_20260301_152329.json
```

### Fields Explained

| Field | Type | Meaning |
|-------|------|---------|
| `pdb_id` | string | PDB ID (e.g. `2I9M`) |
| `fragment_id` | string | e.g. `2I9M_res5` |
| `n_qubits_used` | int | Qubits after mapping |
| `n_atoms_qm` | int | QM region atoms |
| `n_atoms_mm` | int | MM region atoms |
| `hf_energy` | float | Hartree–Fock energy (Ha) |
| `vqe_energy` | float | VQE energy (Ha) |
| `qse_energy` | float | QSE-corrected energy (Ha) |
| `mitigated_energy` | float | ZNE-mitigated energy (Ha) |
| `free_energy` | float | ΔG (Ha) |
| `free_energy_kcal` | float | ΔG (kcal/mol) |
| `vqe_converged` | bool | Did VQE converge? |
| `vqe_n_params` | int | Ansatz parameters |
| `vqe_n_iterations` | int | ADAPT iterations |
| `qse_correction` | float | QSE energy change (Ha) |
| `noise_correction` | float | ZNE correction (Ha) |
| `runtime_seconds` | float | Total time (s) |
| `timestamp` | string | Run timestamp |
| `stage_times` | dict | Time per pipeline stage |

---

## 3. Benchmark Output (`scripts/benchmark.py`)

Running:

```cmd
python scripts\benchmark.py --proteins 2I9M 2JOF
```

produces a comparison table:

### Methods Compared

| Method | Description |
|--------|-------------|
| **Classical FF (AMBER)** | Classical force field; no qubits |
| **Standard VQE (UCCSD)** | Fixed UCCSD ansatz, Jordan–Wigner mapping |
| **ADAPT-VQE + QSE + ZNE (OURS)** | Our pipeline: adaptive ansatz + QSE + ZNE |

### Table Columns

| Column | Meaning |
|--------|---------|
| **Method** | Algorithm used |
| **PDB** | Protein ID |
| **Error (Ha)** | |Energy − exact FCI| |
| **Error (kcal)** | Same error in kcal/mol |
| **Qubits** | N/A for classical, number for quantum |
| **Chem.Acc.** | ✓ if error < 1 kcal/mol |

### Example Row

```
ADAPT-VQE + QSE + ZNE (OURS)  2I9M   0.0017 Ha   1.07 kcal   4 qubits   ✗ ◄
```

- Error ~1.07 kcal/mol (close to chemical accuracy).
- ◄ marks our method.

### Benchmark JSON

Results are saved to `outputs/benchmarks/benchmark_results.json` in the same format as the table.

---

## 4. Clustering Plot (`scripts/visualize_clustering.py`)

When you run `python scripts/visualize_clustering.py`, a 2D scatter plot is saved to `outputs/plots/`.

| Symbol | Meaning |
|--------|---------|
| **★ (star)** | Individual conformations in the cluster |
| **X** | Cluster centroid (mean of conformations in feature space) |
| **Yellow** | Chosen cluster — highest selection score (purity × compactness) |
| **Gray** | Other clusters |

A yellow box explains why the chosen cluster was selected. The representative conformation (used for the quantum simulation) is the star closest to the X in the chosen cluster.

---

## 5. Stage Logs

Each pipeline stage prints its own messages. Examples:

| Stage | Example message |
|-------|------------------|
| 1 | `Loaded: [2I9M] 17 residues \| 4860 atoms` |
| 2 | `Feature vector: 34-dim \| Rg=3.20 Å` |
| 3 | `DISCA \| best_cluster=0 \| purity=1.000 \| converged=True` |
| 4+5 | `Hamiltonian \| NCCO \| (4e, 4o) \| 8 qubits \| E_HF=-7.862 Ha` |
| 6 | `QubitHamiltonian \| method=parity \| n_qubits=6 \| n_terms=25` |
| 7 | `ADAPT iter 1 \| E = +0.632 Ha \| \|g\|_max = 1.12e-01` |
| 8+9 | `QSE \| E_VQE=0.632 Ha \| E_QSE=-0.000 Ha \| ΔE=-0.632 Ha` |
| 10 | `NoiseMitigation (ZNE_Richardson) \| E_mitigated=0.903 Ha` |
| 11 | `[2I9M_res5] @ 310K \| G=0.903 Ha \| UNSTABLE` |

---

## 6. Where Files Are Saved

| Output | Location |
|--------|----------|
| Energy results | `outputs/energies/result_{PDB}_{timestamp}.json` |
| Benchmark results | `outputs/benchmarks/benchmark_results.json` |
| **Full session logs** | `outputs/logs/pipeline_{ts}.log`, `benchmark_{ts}.log`, `show_circuit_{ts}.log`, `visualize_clustering_{ts}.log` — all terminal output (print + logger) tee'd here |
| Logger output | `outputs/logs/run_{timestamp}.log` (from utils logger) |

---

## Next: Circuits and Terms

- **Circuits:** [CIRCUIT_GUIDE.md](CIRCUIT_GUIDE.md)
- **Terms:** [TERMS_GLOSSARY.md](TERMS_GLOSSARY.md)
