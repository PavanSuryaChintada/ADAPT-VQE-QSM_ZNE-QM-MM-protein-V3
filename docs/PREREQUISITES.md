# Prerequisites — Before You Start

Everything you need before running the hybrid quantum-classical protein fragment solver.

---

## System Requirements

| Requirement | Minimum |
|-------------|---------|
| **OS** | Windows 10/11, Linux, or macOS |
| **Python** | 3.10 or higher |
| **RAM** | 8 GB (16 GB recommended for 8+ qubits) |
| **Disk** | ~2 GB (Python + packages + virtual env) |
| **Internet** | Required for first-time PDB downloads |

---

## 1. Check Python Version

```cmd
python --version
```

You need **Python 3.10 or higher**.

- If not installed: [python.org/downloads](https://www.python.org/downloads/)
- During install: ✅ **Check "Add Python to PATH"**

---

## 2. Create Virtual Environment

Navigate to your project folder, then:

```cmd
cd C:\quantum_protein
python -m venv venv
```

Wait ~10 seconds. A `venv` folder will be created.

---

## 3. Activate the Environment

**Do this every time you open a new terminal to work on the project:**

```cmd
venv\Scripts\activate
```

Your prompt should change to:

```
(venv) C:\quantum_protein>
```

The `(venv)` at the start means it's working. ✅

---

## 4. Upgrade pip

```cmd
python -m pip install --upgrade pip
```

---

## 5. Install Packages

### Core science packages (~2 min)

```cmd
pip install numpy scipy matplotlib scikit-learn pandas
```

### Quantum packages (~5 min total)

```cmd
pip install qiskit==0.45.3
```
Wait ~1 min

```cmd
pip install qiskit-aer
```
Wait ~2 min

```cmd
pip install qiskit-algorithms
pip install qiskit-nature
```

### Chemistry + noise mitigation

```cmd
pip install openfermion
pip install mitiq
```

### Biology + utilities

```cmd
pip install biopython loguru rich pyyaml tqdm
```

### PySCF (optional, Windows may fail)

```cmd
pip install pyscf
```

- **Success:** Full quantum chemistry; real Hamiltonians from PDB atoms.
- **Failure:** Pipeline runs in **mock mode** (pre-computed integrals). Fine for testing; not for final paper results.

---

## 6. Verify Installation

```cmd
python scripts\verify_install.py
```

Expected output:

```
[OK] numpy 1.26.4
[OK] scipy 1.13.0
[OK] sklearn 1.0.x
[OK] qiskit
[OK] openfermion
[OK] biopython
[OK] loguru
[WARN] pyscf — not found (mock mode active)
[OK] mitiq

Installation: READY (mock mode)
```

---

## Package Summary

| Category | Packages |
|----------|----------|
| Core | numpy, scipy, matplotlib, scikit-learn, pandas |
| Quantum | qiskit, qiskit-aer, qiskit-algorithms, qiskit-nature, openfermion, mitiq |
| Chemistry | pyscf (optional) |
| Biology | biopython |
| Utilities | loguru, rich, pyyaml, tqdm |

---

## Next Steps

After prerequisites are met, see [IMPLEMENTATION_GUIDE.md](IMPLEMENTATION_GUIDE.md) for the full walkthrough.
