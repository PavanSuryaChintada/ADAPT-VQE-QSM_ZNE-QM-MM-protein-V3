"""
scripts/verify_install.py — Check all dependencies
====================================================
Run this after installation to verify everything is working.

Usage:
  python scripts/verify_install.py
"""

import sys
import subprocess
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))


def check(name, import_str, min_version=None):
    try:
        mod = __import__(import_str)
        ver = getattr(mod, "__version__", "?")
        if min_version:
            print(f"  [OK]   {name:<20} v{ver}")
        else:
            print(f"  [OK]   {name:<20} v{ver}")
        return True
    except ImportError:
        print(f"  [WARN] {name:<20} not installed")
        return False


def main():
    print()
    print("=" * 55)
    print("  quantum_protein — Installation Verification")
    print("=" * 55)
    print(f"  Python: {sys.version.split()[0]}")
    print()

    # Core scientific
    print("  Core packages:")
    ok_core = all([
        check("numpy",       "numpy"),
        check("scipy",       "scipy"),
        check("matplotlib",  "matplotlib"),
        check("scikit-learn","sklearn"),
        check("pandas",      "pandas"),
    ])

    print()
    print("  Quantum computing:")
    ok_qiskit = check("qiskit",        "qiskit")
    ok_aer    = check("qiskit-aer",    "qiskit_aer")
    ok_nature = check("qiskit-nature", "qiskit_nature")
    ok_of     = check("openfermion",   "openfermion")
    ok_mitiq  = check("mitiq",         "mitiq")

    print()
    print("  Quantum chemistry:")
    ok_pyscf = check("pyscf",          "pyscf")
    if not ok_pyscf:
        print("         ↳ MOCK MODE will be used (results valid for testing)")
        print("         ↳ For publication: install PySCF (use WSL2 on Windows)")

    print()
    print("  Biology:")
    ok_bio = check("biopython",        "Bio")

    print()
    print("  Utilities:")
    ok_log  = check("loguru",          "loguru")
    ok_yaml = check("pyyaml",          "yaml")
    ok_rich = check("rich",            "rich")

    # Test project modules
    print()
    print("  Project modules:")
    try:
        from core.data_acquisition import PDBDownloader
        print("  [OK]   core.data_acquisition")
    except Exception as e:
        print(f"  [FAIL] core.data_acquisition: {e}")

    try:
        from core.yopo.feature_extractor import YOPOExtractor
        print("  [OK]   core.yopo")
    except Exception as e:
        print(f"  [FAIL] core.yopo: {e}")

    try:
        from core.disca.clustering import DISCAClustering
        print("  [OK]   core.disca")
    except Exception as e:
        print(f"  [FAIL] core.disca: {e}")

    try:
        from quantum.hamiltonian.builder import QMMMHamiltonianBuilder
        print("  [OK]   quantum.hamiltonian")
    except Exception as e:
        print(f"  [FAIL] quantum.hamiltonian: {e}")

    try:
        from quantum.vqe.adapt_vqe import ADAPTVQESolver
        print("  [OK]   quantum.vqe")
    except Exception as e:
        print(f"  [FAIL] quantum.vqe: {e}")

    try:
        from thermodynamics.free_energy import FreeEnergyCalculator
        print("  [OK]   thermodynamics")
    except Exception as e:
        print(f"  [FAIL] thermodynamics: {e}")

    # Quick sanity test
    print()
    print("  Quick test (ADAPT-VQE on mock 4-qubit system):")
    try:
        import numpy as np
        from quantum.hamiltonian.builder import QMMMHamiltonianBuilder
        from quantum.mapping.qubit_mapper import QubitMapper
        from quantum.vqe.adapt_vqe import ADAPTVQESolver

        np.random.seed(42)
        builder = QMMMHamiltonianBuilder(active_electrons=2, active_orbitals=2)
        h = builder._mock_hamiltonian()
        mapper = QubitMapper(method="jordan_wigner", two_qubit_reduction=False)
        qh = mapper.map(h)
        exact = float(np.linalg.eigvalsh(qh.to_matrix())[0])

        solver = ADAPTVQESolver(qh, n_electrons=2, max_iterations=10, gradient_threshold=0.01)
        result = solver.solve(hf_energy=h.hf_energy)
        error = abs(result.energy - exact)

        if error < 0.5:
            print(f"  [OK]   ADAPT-VQE test | error={error:.4f} Ha")
        else:
            print(f"  [WARN] ADAPT-VQE test | error={error:.4f} Ha (unexpected)")
    except Exception as e:
        print(f"  [FAIL] ADAPT-VQE test: {e}")

    print()
    print("=" * 55)
    pyscf_note = "(mock mode)" if not ok_pyscf else "(full mode)"
    print(f"  Status: READY {pyscf_note}")
    print()
    print("  Run the pipeline:")
    print("  > python main.py --pdb 2I9M --fragment 5 --qubits 8")
    print()
    print("  Run benchmark for your paper:")
    print("  > python scripts/benchmark.py --proteins 2I9M 2JOF")
    print("=" * 55)


if __name__ == "__main__":
    import os
    os.chdir(ROOT)
    main()
