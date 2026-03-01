# WSL Setup Guide — Quantum Protein Pipeline

Follow these steps to set up Windows Subsystem for Linux (WSL) and run the pipeline.

---

## 1. Install WSL

### Option A: PowerShell (Admin)
1. Open **PowerShell as Administrator**
2. Run:
   ```powershell
   wsl --install
   ```
3. Restart your PC if prompted.
4. On first launch, create a **username** and **password** for your Linux user.

### Option B: Manual (Windows Features)
1. Press `Win + R`, type `optionalfeatures`, Enter.
2. Enable **Windows Subsystem for Linux** and **Virtual Machine Platform**.
3. Restart.
4. Open **Microsoft Store** → search **Ubuntu** → Install.

---

## 2. Install WSL 2 (recommended)

In PowerShell:
```powershell
wsl --set-default-version 2
```

Check status:
```powershell
wsl -l -v
```

---

## 3. Launch WSL

- **Start Menu** → Ubuntu (or your distro).
- Or in PowerShell/CMD: `wsl`

---

## 4. Set Up Python and Dependencies

Copy-paste these commands in your WSL terminal (Ubuntu):

```bash
sudo apt update
sudo apt install python3 python3-pip python3-venv -y
cd /mnt/c/Users/praveeeee/Desktop/qq\ protein/quantum_protein_v3-3
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install numpy scipy matplotlib scikit-learn pandas
pip install qiskit qiskit-aer qiskit-algorithms qiskit-nature openfermion mitiq
pip install pyscf biopython loguru rich pyyaml tqdm
```

---

## 5. Run the Pipeline

```bash
source venv/bin/activate   # if not already active
python main.py --pdb 2I9M --fragment 5 --qubits 8 --electrons 4 --orbitals 4
```

---

## Tips

| Issue | Fix |
|-------|-----|
| `venv` not found | Make sure you ran `cd` to the project folder first |
| Permission denied | Use `chmod +x scripts/wsl_setup.sh` or run with `bash scripts/wsl_setup.sh` |
| PySCF build fails | Run `sudo apt install build-essential` first |
| Out of memory | Close other apps; pipeline uses first PDB model only now |

---

## Path Reference

- Windows path: `C:\Users\praveeeee\Desktop\qq protein\quantum_protein_v3-3`
- WSL path: `/mnt/c/Users/praveeeee/Desktop/qq protein/quantum_protein_v3-3`
