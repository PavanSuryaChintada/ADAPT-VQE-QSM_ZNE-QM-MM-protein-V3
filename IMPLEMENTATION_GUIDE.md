# Implementation Guide — Step-by-Step Walkthrough

Step-by-step guide to set up, run, and understand the hybrid quantum-classical protein fragment solver.

---

## Prerequisites

Follow [PREREQUISITES.md](PREREQUISITES.md) first. Ensure:

- Python 3.10+
- Virtual environment created and activated
- All packages installed
- `python scripts\verify_install.py` passes

---

## Step 1: Open Terminal and Go to Project

```cmd
cd C:\quantum_protein
```

*(Replace with your actual project path.)*

---

## Step 2: Activate Virtual Environment

**Do this every time you open a new terminal.**

```cmd
venv\Scripts\activate
```

Your prompt should show:

```
(venv) C:\quantum_protein>
```

The `(venv)` at the start means it's active. ✅

---

## Step 3: Download the Protein

The pipeline downloads PDB files automatically. To pre-download Chignolin (best starting protein):

```cmd
python -c "from core.data_acquisition import PDBDownloader; PDBDownloader('data/raw').fetch('2I9M')"
```

You should see:

```
Downloading 2I9M from RCSB...
Downloaded 2I9M (45 KB) → data\raw\2I9M.pdb
Loaded: [2I9M] 10 residues | 138 atoms
```

### Verify download

```cmd
dir data\raw\
```

You should see `2I9M.pdb` listed.

---

## Step 4: Run the Pipeline

```cmd
python main.py --pdb 2I9M --fragment 5 --qubits 8 --electrons 4 --orbitals 4
```

Typical run time: **3–5 minutes**.

The pipeline will:

1. Load PDB structure
2. Extract QM/MM fragment
3. Build YOPO features
4. Run DISCA clustering
5. Construct Hamiltonian (QM/MM + PySCF or mock)
6. Map to qubits (parity)
7. Run ADAPT-VQE
8. Apply QSE refinement
9. Apply noise mitigation (ZNE)
10. Compute free energy

---

## Step 5: Run the Benchmark (for Paper / Presentation)

```cmd
python scripts\benchmark.py --proteins 2I9M 2JOF
```

This generates a comparison table: Classical FF vs Standard VQE vs ADAPT-VQE + QSE + ZNE.

---

## Step 6: Inspect Results

```cmd
dir outputs\energies\
```

You will see JSON files such as `result_2I9M_20260301_152329.json`.

---

## Step 7: Download Second Protein (2JOF — Trp-cage)

```cmd
python -c "from core.data_acquisition import PDBDownloader; PDBDownloader('data/raw').fetch('2JOF')"
```

Then run the pipeline:

```cmd
python main.py --pdb 2JOF --fragment 5 --qubits 8 --electrons 4 --orbitals 4
```

---

## Manual PDB Download

If automatic download fails (no internet, firewall, etc.):

1. Open: [https://files.rcsb.org/download/2I9M.pdb](https://files.rcsb.org/download/2I9M.pdb)
2. Save as `data\raw\2I9M.pdb` in your project folder
3. Run the pipeline again

---

## Understanding the Output

### What each term means (for presentations)

| Term | Plain English |
|------|---------------|
| **HF energy** | Classical baseline — what a physics textbook gives you |
| **VQE energy** | Your quantum circuit’s answer after training |
| **QSE energy** | VQE answer improved by quantum subspace expansion |
| **Mitigated energy** | Final answer after correcting for hardware noise |
| **Error (Ha)** | Distance from exact — lower is better |
| **Chemical accuracy** | Gold standard: error < 0.0016 Ha (≈1 kcal/mol) |
| **Free energy G** | Is this protein conformation stable? **Negative = stable** |

### Units

- **Ha (Hartree):** Energy unit in atomic physics. 1 Ha ≈ 627.5 kcal/mol
- **Chemical accuracy:** ~1 kcal/mol ≈ 0.0016 Ha
- **kcal/mol:** Common unit in biochemistry and protein folding

---

## Quick Reference — Daily Usage

### Open CMD, run pipeline

```cmd
# 1. Go to project
cd C:\quantum_protein

# 2. Activate environment (every time)
venv\Scripts\activate

# 3. Run pipeline on Chignolin
python main.py --pdb 2I9M --fragment 5 --qubits 8 --electrons 4 --orbitals 4

# 4. Run benchmark
python scripts\benchmark.py --proteins 2I9M 2JOF

# 5. Verify install
python scripts\verify_install.py
```

### Quantum circuit inspection

```cmd
python scripts\show_circuit.py
```

With full Hamiltonian details:

```cmd
python scripts\show_circuit.py --full
```

For simpler molecules (e.g. H₂):

```cmd
python scripts\show_circuit.py --molecule H2 --electrons 2 --orbitals 2
```

---

## Command Cheat Sheet

| Task | Command |
|------|---------|
| Verify install | `python scripts\verify_install.py` |
| Download Chignolin | `python -c "from core.data_acquisition import PDBDownloader; PDBDownloader('data/raw').fetch('2I9M')"` |
| Download Trp-cage | `python -c "from core.data_acquisition import PDBDownloader; PDBDownloader('data/raw').fetch('2JOF')"` |
| Run pipeline (Chignolin) | `python main.py --pdb 2I9M --fragment 5 --qubits 8 --electrons 4 --orbitals 4` |
| Run pipeline (Trp-cage) | `python main.py --pdb 2JOF --fragment 5 --qubits 8 --electrons 4 --orbitals 4` |
| Run benchmark | `python scripts\benchmark.py --proteins 2I9M 2JOF` |
| Show circuit | `python scripts\show_circuit.py` |
| Show circuit (full) | `python scripts\show_circuit.py --full` |
| Simpler molecule | `python scripts\show_circuit.py --molecule H2 --electrons 2 --orbitals 2` |
| **Visualize clustering** | `python scripts\visualize_clustering.py` |
| Clustering (2JOF) | `python scripts\visualize_clustering.py --pdb 2JOF --save` |

---

## Clustering Visualization — Understanding the Plot

Run `python scripts\visualize_clustering.py` to generate a 2D plot of DISCA clusters.

### Symbols

| Symbol | Meaning |
|--------|---------|
| **★ (star)** | Individual conformations — each point is one conformation from the ensemble |
| **X** | Cluster centroid — the mean position of all conformations in that cluster |
| **Yellow** | Chosen cluster — highlighted because it has the highest selection score |

### Selection Criterion

- **Selection score** = purity × compactness (both in [0,1])
- **Purity** = mean responsibility (confidence that members belong to the cluster)
- **Compactness** = 1 / (1 + within-cluster variance); higher = tighter cluster
- The **chosen cluster** (yellow) is the one with the highest score
- The **representative conformation** used for quantum simulation is the star closest to the X in the chosen cluster

### Output

Plot saved to `outputs/plots/clustering_{PDB}_res{residue}.png`.

---

## Common Errors

See [TROUBLESHOOTING.md](TROUBLESHOOTING.md) for detailed fixes.

| Error | Quick fix |
|-------|-----------|
| `python` not recognized | Reinstall Python; check "Add to PATH" |
| No module named X | Run `venv\Scripts\activate`, then `pip install X` |
| PySCF not found | Normal on Windows — mock mode works |
| ConnectionError | Manual PDB download (see above) |
| MemoryError | Use `--electrons 4 --orbitals 4 --qubits 8` |
| `activate` fails | Run `Set-ExecutionPolicy RemoteSigned` in PowerShell (Admin) |

---

## Where to Look Next

- **Prerequisites:** [PREREQUISITES.md](PREREQUISITES.md)
- **Troubleshooting:** [TROUBLESHOOTING.md](TROUBLESHOOTING.md)
- **Setup details:** [SETUP.md](SETUP.md)
- **Project overview:** [README.md](README.md)
