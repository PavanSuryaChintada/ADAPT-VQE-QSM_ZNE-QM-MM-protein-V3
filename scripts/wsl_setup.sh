#!/bin/bash
# WSL setup script for quantum protein pipeline
# Run from WSL: bash scripts/wsl_setup.sh

set -e

echo "=== Updating system ==="
sudo apt update

echo "=== Installing Python ==="
sudo apt install -y python3 python3-pip python3-venv build-essential

echo "=== Navigating to project ==="
cd /mnt/c/Users/praveeeee/Desktop/qq\ protein/quantum_protein_v3-3

echo "=== Creating venv ==="
python3 -m venv venv
source venv/bin/activate

echo "=== Upgrading pip ==="
pip install --upgrade pip

echo "=== Installing core deps ==="
pip install numpy scipy matplotlib scikit-learn pandas

echo "=== Installing quantum deps ==="
pip install qiskit qiskit-aer qiskit-algorithms qiskit-nature openfermion mitiq

echo "=== Installing chemistry & utils ==="
pip install pyscf biopython loguru rich pyyaml tqdm

echo "=== Running pipeline ==="
python main.py --pdb 2I9M --fragment 5 --qubits 8 --electrons 4 --orbitals 4

echo "=== Done ==="
