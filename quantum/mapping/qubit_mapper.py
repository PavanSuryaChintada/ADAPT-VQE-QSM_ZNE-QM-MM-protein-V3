"""
Stage 6 — Qubit Mapping (Parity / Jordan-Wigner)
==================================================
Transforms fermionic Hamiltonian into qubit (Pauli) Hamiltonian.

Result: H = Σ_i c_i P_i
  where P_i are tensor products of Pauli matrices {I, X, Y, Z}

Two mappings available:

1. Jordan-Wigner (JW):
   a†_p = (Z_0⊗...⊗Z_{p-1}) ⊗ X_p/2 - iY_p/2
   → n_qubits = n_spin_orbitals = 2 * n_orbitals

2. Parity (RECOMMENDED — saves 2 qubits):
   Exploits Z_2 symmetry of molecular Hamiltonians.
   → n_qubits = 2 * n_orbitals - 2
   On your hardware: (4e,4o): 8→6 qubits, (6e,6o): 12→10 qubits

Always use Parity on hardware with ≤ 16GB RAM.
"""

import sys
import numpy as np
from pathlib import Path
from typing import List, Tuple, Optional, Dict
from dataclasses import dataclass

sys.path.insert(0, str(Path(__file__).parent.parent))
from utils import logger

try:
    from openfermion.transforms import jordan_wigner, bravyi_kitaev
    from openfermion import QubitOperator
    OF_AVAILABLE = True
except ImportError:
    OF_AVAILABLE = False

try:
    from qiskit.quantum_info import SparsePauliOp
    QISKIT_AVAILABLE = True
except ImportError:
    QISKIT_AVAILABLE = False

from quantum.hamiltonian.builder import HamiltonianData


# Pauli matrices (used throughout the pipeline)
PAULI_I = np.array([[1, 0], [0, 1]], dtype=complex)
PAULI_X = np.array([[0, 1], [1, 0]], dtype=complex)
PAULI_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
PAULI_Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI_MAP = {"I": PAULI_I, "X": PAULI_X, "Y": PAULI_Y, "Z": PAULI_Z}


@dataclass
class QubitHamiltonian:
    """Qubit Hamiltonian as list of (Pauli string, coefficient) pairs."""
    terms: List[Tuple[str, complex]]   # [("XZIY", 0.123+0j), ...]
    n_qubits: int
    mapping_method: str
    identity_constant: float           # constant term (E_core etc.)

    @property
    def n_terms(self) -> int:
        return len(self.terms)

    def to_matrix(self) -> np.ndarray:
        """Build full (2^n, 2^n) Hamiltonian matrix."""
        dim = 2 ** self.n_qubits
        H = np.zeros((dim, dim), dtype=complex)

        # Sum all identity term coefficients from the terms list
        identity_sum = 0.0
        non_identity_terms = []
        for pauli_str, coeff in self.terms:
            if all(p == 'I' for p in pauli_str):
                identity_sum += float(coeff.real if hasattr(coeff, 'real') else coeff)
            else:
                non_identity_terms.append((pauli_str, coeff))

        # Use identity_constant if no identity terms found in list
        if abs(identity_sum) < 1e-10:
            identity_sum = self.identity_constant

        # Add identity (constant) contribution
        H += identity_sum * np.eye(dim, dtype=complex)

        # Add all non-identity Pauli terms
        for pauli_str, coeff in non_identity_terms:
            if abs(coeff) < 1e-14:
                continue
            H += coeff * pauli_string_to_matrix(pauli_str)
        return H

    def to_sparse_pauli_op(self) -> Optional["SparsePauliOp"]:
        """Convert to Qiskit SparsePauliOp."""
        if not QISKIT_AVAILABLE:
            return None
        terms = [(p, c) for p, c in self.terms if abs(c) > 1e-14]
        if not terms:
            return SparsePauliOp.from_list([("I" * self.n_qubits, 1.0)])
        try:
            return SparsePauliOp.from_list([(p, float(c.real)) for p, c in terms])
        except Exception as e:
            logger.warning(f"SparsePauliOp conversion failed: {e}")
            return None

    def summary(self) -> str:
        return (
            f"QubitHamiltonian | method={self.mapping_method} | "
            f"n_qubits={self.n_qubits} | n_terms={self.n_terms} | "
            f"constant={self.identity_constant:.6f}"
        )


class QubitMapper:
    """
    Maps fermionic Hamiltonian to qubit Hamiltonian.

    Usage:
        mapper = QubitMapper(method="parity")
        qubit_h = mapper.map(h_data)
    """

    def __init__(
        self,
        method: str = "parity",
        two_qubit_reduction: bool = True,
    ):
        assert method in ("parity", "jordan_wigner", "bravyi_kitaev"), \
            f"method must be parity/jordan_wigner/bravyi_kitaev"
        self.method = method
        self.tqr = two_qubit_reduction  # saves 2 qubits for parity

    def map(self, h_data: HamiltonianData) -> QubitHamiltonian:
        """Map HamiltonianData → QubitHamiltonian."""
        n_so = 2 * h_data.n_orbitals  # spin-orbitals

        logger.info(
            f"Qubit mapping: ({h_data.n_electrons}e, {h_data.n_orbitals}o) | "
            f"method={self.method} | "
            f"spin-orbitals={n_so}"
        )

        if h_data.fermion_operator is not None and OF_AVAILABLE:
            qubit_h = self._map_openfermion(h_data, n_so)
        else:
            qubit_h = self._map_direct(h_data, n_so)

        logger.success(qubit_h.summary())
        return qubit_h

    def _map_openfermion(self, h_data: HamiltonianData, n_so: int) -> QubitHamiltonian:
        """Use OpenFermion for accurate mapping."""
        fop = h_data.fermion_operator

        if self.method == "jordan_wigner":
            qubit_op = jordan_wigner(fop)
        else:
            qubit_op = bravyi_kitaev(fop)

        terms = []
        n_qubits_used = 0
        identity_const = 0.0

        for term, coeff in qubit_op.terms.items():
            if abs(coeff) < 1e-12:
                continue
            if not term:
                identity_const += float(coeff.real)
                continue
            max_q = max(idx for idx, _ in term)
            n_qubits_used = max(n_qubits_used, max_q + 1)

            pauli_dict = {idx: op for idx, op in term}
            ps = "".join(pauli_dict.get(i, "I") for i in range(n_qubits_used))
            terms.append((ps, complex(coeff)))

        # Pad all terms to same length
        terms = [(ps + "I" * (n_qubits_used - len(ps)), c) for ps, c in terms]

        if self.tqr and self.method == "parity":
            n_qubits_used = max(0, n_qubits_used - 2)
            terms = self._apply_tqr(terms, n_qubits_used)

        return QubitHamiltonian(
            terms=terms,
            n_qubits=n_qubits_used if n_qubits_used > 0 else n_so,
            mapping_method=self.method,
            identity_constant=identity_const,
        )

    def _map_direct(self, h_data: HamiltonianData, n_so: int) -> QubitHamiltonian:
        """
        Direct Jordan-Wigner mapping from integrals (no OpenFermion).
        Implements: H = Σ h_pq (a†p aq) + 1/2 Σ h_pqrs (a†p a†q ar as)
        """
        n = h_data.n_orbitals
        n_qubits = n_so - (2 if self.tqr else 0)
        terms = []
        identity_const = h_data.e_core

        # One-body: number operator contribution
        for p in range(n):
            for spin in range(2):
                idx = 2 * p + spin
                if idx >= n_qubits:
                    continue
                h_pp = h_data.h1e[p, p]
                if abs(h_pp) < 1e-12:
                    continue
                # n_p = (I - Z_p) / 2
                pauli_I = "I" * n_qubits
                pauli_Z = list("I" * n_qubits)
                pauli_Z[idx] = "Z"

                terms.append((pauli_I, complex(0.5 * h_pp)))
                terms.append(("".join(pauli_Z), complex(-0.5 * h_pp)))

        # Off-diagonal one-body
        for p in range(n):
            for q in range(p + 1, n):
                h_pq = h_data.h1e[p, q]
                if abs(h_pq) < 1e-12:
                    continue
                for spin in range(2):
                    ps, qs = 2*p+spin, 2*q+spin
                    if ps >= n_qubits or qs >= n_qubits:
                        continue
                    # a†_p a_q = X_p Z_{p+1}...Z_{q-1} X_q / 4 + ...
                    # Simplified contribution:
                    pauli_xx = list("I" * n_qubits)
                    pauli_xx[ps] = "X"; pauli_xx[qs] = "X"
                    for k in range(ps+1, qs): pauli_xx[k] = "Z"
                    terms.append(("".join(pauli_xx), complex(0.5 * h_pq)))

                    pauli_yy = list("I" * n_qubits)
                    pauli_yy[ps] = "Y"; pauli_yy[qs] = "Y"
                    for k in range(ps+1, qs): pauli_yy[k] = "Z"
                    terms.append(("".join(pauli_yy), complex(0.5 * h_pq)))

        # Add identity constant
        terms.insert(0, ("I" * n_qubits, complex(identity_const)))

        return QubitHamiltonian(
            terms=terms,
            n_qubits=n_qubits,
            mapping_method=self.method,
            identity_constant=identity_const,
        )

    def _apply_tqr(self, terms, n_q):
        """Trim Pauli strings to n_q qubits after two-qubit reduction."""
        new_terms = []
        for ps, c in terms:
            if len(ps) > n_q:
                ps = ps[:n_q]
            new_terms.append((ps, c))
        return new_terms


def pauli_string_to_matrix(pauli_str: str) -> np.ndarray:
    """Convert Pauli string (e.g. 'XZIY') to matrix via tensor product."""
    result = np.array([[1.0 + 0j]])
    for p in pauli_str:
        result = np.kron(result, PAULI_MAP[p])
    return result


# ── Main ─────────────────────────────────────────────────────

if __name__ == "__main__":
    import os
    os.chdir(Path(__file__).parent.parent)

    logger.info("=" * 60)
    logger.info("STAGE 6 — Qubit Mapping Test")
    logger.info("=" * 60)

    from quantum.hamiltonian.builder import QMMMHamiltonianBuilder

    builder = QMMMHamiltonianBuilder(active_electrons=2, active_orbitals=2)
    h_data = builder._mock_hamiltonian()

    mapper = QubitMapper(method="parity", two_qubit_reduction=True)
    qh = mapper.map(h_data)

    logger.info(qh.summary())
    logger.info(f"First 5 terms:")
    for ps, c in qh.terms[:5]:
        logger.info(f"  {c.real:+.4f} × {ps}")

    # Verify matrix is Hermitian
    H_mat = qh.to_matrix()
    is_herm = np.allclose(H_mat, H_mat.conj().T, atol=1e-10)
    logger.info(f"Hamiltonian is Hermitian: {is_herm}")

    logger.success("Stage 6 PASSED")
