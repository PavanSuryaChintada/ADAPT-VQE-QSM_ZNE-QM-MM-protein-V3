# Troubleshooting — Common Errors and Fixes

Quick reference for common issues and their solutions.

---

## Common Errors

| Error message | What it means | Fix |
|---------------|---------------|-----|
| `'python' is not recognized` | Python not in PATH | Reinstall Python; during setup, check **"Add Python to PATH"** |
| `No module named 'qiskit'` | venv not activated or package not installed | Run `venv\Scripts\activate` first, then `pip install qiskit` |
| `No module named 'X'` | Package missing | Activate venv, then `pip install X` |
| `ModuleNotFoundError: pyscf` | PySCF not installed | Normal on Windows. Pipeline runs in **mock mode** — fine for testing |
| `ConnectionError: Could not download` | No internet or RCSB unreachable | Download PDB manually (see [IMPLEMENTATION_GUIDE.md](IMPLEMENTATION_GUIDE.md#manual-pdb-download)) |
| `MemoryError` | Too many qubits | Use `--qubits 8 --electrons 4 --orbitals 4` |
| `venv\Scripts\activate` gives error | PowerShell execution policy | Run `Set-ExecutionPolicy RemoteSigned -Scope CurrentUser` in PowerShell |
| `SSL certificate error` | Certificate / network issue | Run `pip install certifi` then retry |
| `ImportError: cannot import name 'X'` | Code / version mismatch | Ensure all files from the project are present; reinstall packages |
| `pip is not recognized` | pip or PATH issue | Run `python -m pip` instead of `pip` |
| `Permission denied` | Writing to protected folder | Run cmd/PowerShell as Administrator or use a non-system path |

---

## Detailed Fixes

### 1. `'python' is not recognized`

**Cause:** Python is not in your system PATH.

**Fix:**
1. Reinstall Python from [python.org](https://www.python.org/downloads/)
2. During installation, **check** "Add Python to PATH"
3. Restart the terminal
4. Or use the full path, e.g. `C:\Python310\python.exe`

---

### 2. `No module named 'X'` (qiskit, numpy, etc.)

**Cause:** Virtual environment not activated, or package not installed.

**Fix:**
1. Ensure venv is active: `venv\Scripts\activate` — you should see `(venv)` in the prompt
2. Install the missing package: `pip install qiskit` (or whatever is missing)
3. Verify: `python scripts\verify_install.py`

---

### 3. PySCF install fails (Windows)

**Cause:** PySCF has C/Fortran extensions; build tools may be missing on Windows.

**Fix:**
- **Option A:** Use **mock mode** (default). Pipeline runs with pre-computed integrals; good for testing.
- **Option B:** Use WSL2 (Ubuntu) and install packages inside Linux.
- **Option C:** Install Visual Studio Build Tools, then retry `pip install pyscf`.

---

### 4. ConnectionError: Could not download PDB

**Cause:** No internet, firewall, or RCSB server issue.

**Fix:**
1. Check internet connection
2. Manually download from [files.rcsb.org/download/2I9M.pdb](https://files.rcsb.org/download/2I9M.pdb)
3. Save to `data\raw\2I9M.pdb` in your project folder
4. Run the pipeline again

---

### 5. MemoryError or slow / stuck run

**Cause:** Active space (orbitals/electrons) too large; Hamiltonian / state vector grows exponentially.

**Fix:**
- Reduce: `python main.py --pdb 2I9M --electrons 4 --orbitals 4 --qubits 8`
- For 16 GB RAM, stay at 8 qubits or below (4e, 4o).

---

### 6. PowerShell: `venv\Scripts\activate` fails

**Cause:** Execution policy blocks running scripts.

**Fix (PowerShell as Administrator):**
```powershell
Set-ExecutionPolicy RemoteSigned -Scope CurrentUser
```

Or use **Command Prompt (cmd)** instead of PowerShell for `activate`.

---

### 7. Wrong path / "can't find project"

**Fix:** Always `cd` into your project folder first:
```cmd
cd C:\quantum_protein
venv\Scripts\activate
python main.py --pdb 2I9M --fragment 5 --qubits 8 --electrons 4 --orbitals 4
```

---

### 8. VQE not converging / NaN energies

**Cause:** Hamiltonian or optimizer settings; sometimes noisy mock data.

**Fix:**
- Check `configs/config.yaml` for `gradient_threshold`, `convergence_threshold`
- Try `--seed 42` for reproducibility
- If using mock mode, exact energies are approximate; focus on pipeline behavior

---

## Quick Diagnostics

```cmd
# 1. Python version
python --version

# 2. venv active?
# Look for (venv) in prompt

# 3. Packages OK
python scripts\verify_install.py

# 4. PDB present
dir data\raw\
```

---

## Where to Look Next

- Full setup: [PREREQUISITES.md](PREREQUISITES.md)
- Step-by-step run: [IMPLEMENTATION_GUIDE.md](IMPLEMENTATION_GUIDE.md)
- Installation details: [SETUP.md](SETUP.md)
