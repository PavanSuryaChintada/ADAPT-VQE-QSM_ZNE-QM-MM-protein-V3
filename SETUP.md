# SETUP.md — Complete Installation Guide
# Plain Python on Windows 11 (No Anaconda Required)

---

## Step 0 — Check Your Python Version

Open Command Prompt (cmd) or PowerShell:

```cmd
python --version
```

You need Python 3.10 or higher.
If not installed: https://www.python.org/downloads/
During install: ✅ CHECK "Add Python to PATH"

---

## Step 1 — Create a Virtual Environment

This keeps your project isolated from system Python.

```cmd
cd C:\path\to\quantum_protein
python -m venv venv
```

Activate it (do this EVERY time you work on the project):

```cmd
venv\Scripts\activate
```

You should see `(venv)` in your prompt.

---

## Step 2 — Upgrade pip

```cmd
python -m pip install --upgrade pip
```

---

## Step 3 — Install Core Packages First

```cmd
pip install numpy scipy matplotlib scikit-learn pandas
```

---

## Step 4 — Install Quantum Packages

```cmd
pip install qiskit==0.45.3
pip install qiskit-aer==0.13.3
pip install qiskit-algorithms==0.3.0
pip install qiskit-nature==0.7.2
pip install openfermion==1.6.1
pip install mitiq==0.36.0
```

---

## Step 5 — Install PySCF (Quantum Chemistry Engine)

### Option A: Direct (may work on Windows with build tools):
```cmd
pip install pyscf
```

### Option B: If Option A fails (recommended for Windows):
Install WSL2 (Windows Subsystem for Linux):
1. Open PowerShell as Administrator
2. Run: `wsl --install`
3. Restart computer
4. Open Ubuntu terminal from Start Menu
5. Run all the above pip installs inside Ubuntu
6. Run the project from Ubuntu terminal

### Option C: Use mock mode (no PySCF):
The project runs in MOCK MODE automatically if PySCF is missing.
Mock mode uses pre-computed integrals. 
Good for: testing pipeline, debugging, small runs.
NOT for: final paper results (need real PySCF).

---

## Step 6 — Install Remaining Packages

```cmd
pip install biopython loguru rich tqdm click pyyaml h5py
```

---

## Step 7 — Verify Installation

```cmd
python scripts/verify_install.py
```

Expected output:
```
[OK] numpy 1.26.4
[OK] scipy 1.13.0
[OK] sklearn 1.5.0
[OK] qiskit
[OK] openfermion
[OK] biopython
[OK] loguru
[WARN] pyscf — not found (mock mode active)
[OK] mitiq

Installation: READY (mock mode)
```

---

## Step 8 — Run Your First Test

```cmd
python main.py --pdb 2I9M --fragment 5 --qubits 8
```

This will:
1. Download Chignolin (2I9M) from PDB automatically
2. Run all 12 pipeline stages
3. Print energy results
4. Save output to outputs/energies/

---

## Troubleshooting

### "python is not recognized"
→ Reinstall Python and check "Add to PATH"

### "No module named X"
→ Make sure venv is activated: `venv\Scripts\activate`

### PySCF install fails
→ Use WSL2 (Option B above) or mock mode (Option C)

### Memory error at 14+ qubits
→ Reduce: `python main.py --electrons 4 --orbitals 4`

### SSL certificate error downloading PDB
→ Run: `pip install certifi` then retry

---

## Quick Reference — Daily Usage

```cmd
# 1. Open cmd, navigate to project
cd C:\path\to\quantum_protein

# 2. Activate environment
venv\Scripts\activate

# 3. Run pipeline
python main.py --pdb 2I9M --fragment 5 --qubits 8

# 4. Run benchmark (for paper)
python scripts/benchmark.py --proteins 2I9M 2JOF

# 5. Check results
type outputs\energies\result_2I9M_latest.json
```
