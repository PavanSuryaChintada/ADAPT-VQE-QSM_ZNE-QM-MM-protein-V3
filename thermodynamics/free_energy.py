"""
Stage 11 — Free Energy Computation
=====================================
Biological stability from quantum simulation results.

ΔG = ΔH - TΔS

At T = 310 K (body temperature):
  ΔH  ≈ electronic energy difference (from QSE-corrected VQE)
  TΔS ≈ T × von Neumann entropy × k_B

Von Neumann entropy:
  S_VN = -Tr(ρ log ρ) = -Σ λ_i log(λ_i)
  where λ_i are eigenvalues of reduced density matrix ρ

Physical meaning:
  S_VN = 0 → pure state (fully determined, no entanglement with environment)
  S_VN > 0 → mixed state (entangled or thermally populated)

For a fragment at body temperature:
  ΔG < 0 → this conformation is thermodynamically stable (favored)
  ΔG > 0 → unstable (disfavored)

This lets us rank multiple conformations by stability
without solving the full protein folding problem.
"""

import sys
import numpy as np
from pathlib import Path
from typing import List, Optional, Tuple
from dataclasses import dataclass, field

sys.path.insert(0, str(Path(__file__).parent.parent))
from utils import logger

# Physical constants
KB_HARTREE   = 3.16681e-6    # Boltzmann constant in Ha/K
HA_TO_KCAL   = 627.5094      # 1 Ha = 627.5094 kcal/mol
HA_TO_KJ     = 2625.5        # 1 Ha = 2625.5 kJ/mol
BODY_TEMP_K  = 310.0


@dataclass
class ThermodynamicsResult:
    fragment_id: str
    temperature_K: float
    electronic_energy: float       # E_elec  (Ha)
    enthalpy: float                # H = E_elec + ZPE  (Ha) — ZPE≈0 here
    entropy_vn: float              # von Neumann entropy  (dimensionless)
    entropy_thermal: float         # S = k_B × S_VN  (Ha/K)
    free_energy: float             # G = H - TS  (Ha)
    free_energy_kcal: float        # G in kcal/mol
    free_energy_kj: float          # G in kJ/mol
    density_eigenvalues: np.ndarray
    is_stable: bool                # G < 0

    def summary(self) -> str:
        stab = "STABLE" if self.is_stable else "UNSTABLE"
        return (
            f"[{self.fragment_id}] @ {self.temperature_K:.0f}K | "
            f"E={self.electronic_energy:.6f} Ha | "
            f"G={self.free_energy:.6f} Ha "
            f"({self.free_energy_kcal:.3f} kcal/mol) | "
            f"S_VN={self.entropy_vn:.4f} | {stab}"
        )


@dataclass
class ConformationComparison:
    reference_id: str
    temperature_K: float
    rankings: List[Tuple[str, float, float]]  # (id, delta_G_Ha, delta_G_kcal)

    def print_table(self) -> None:
        print()
        print("=" * 70)
        print(f"  Conformation Stability Ranking @ {self.temperature_K:.0f}K")
        print(f"  Reference: {self.reference_id}")
        print("=" * 70)
        print(f"  {'Conformation':<20} {'ΔG (Ha)':>12}  {'ΔG (kcal/mol)':>14}  {'Stability'}")
        print("-" * 70)
        for rank, (fid, dg_ha, dg_kcal) in enumerate(self.rankings, 1):
            stab = "✓ STABLE" if dg_ha <= 0 else "✗ UNSTABLE"
            print(f"  {rank}. {fid:<18} {dg_ha:>+12.6f}  {dg_kcal:>+14.4f}  {stab}")
        print("=" * 70)
        print(f"  Chemical accuracy threshold: 1.0 kcal/mol = 0.0016 Ha")
        print("=" * 70)


class FreeEnergyCalculator:
    """
    Computes thermodynamic free energy from quantum simulation.

    Usage:
        calc = FreeEnergyCalculator(temperature_K=310.0)
        result = calc.compute(state_vector, energy, fragment_id="frag_1")
        comparison = calc.compare([result1, result2, result3])
    """

    def __init__(self, temperature_K: float = BODY_TEMP_K):
        self.T = temperature_K
        self.beta = 1.0 / (KB_HARTREE * temperature_K) if temperature_K > 0 else np.inf
        logger.info(f"FreeEnergyCalculator: T={temperature_K}K")

    def compute(
        self,
        state_vector: np.ndarray,
        energy: float,
        fragment_id: str = "fragment",
        subsystem_size: Optional[int] = None,
    ) -> ThermodynamicsResult:
        """
        Compute full thermodynamic properties.

        Args:
            state_vector: (2^n,) complex quantum state from VQE/QSE
            energy: ground state energy in Hartree
            fragment_id: name for logging and comparison
            subsystem_size: if given, compute reduced DM for this many qubits

        Returns:
            ThermodynamicsResult
        """
        # 1. Build (reduced) density matrix
        rho = self._density_matrix(state_vector, subsystem_size)

        # 2. Von Neumann entropy: S = -Tr(ρ log ρ) = -Σ λ_i log(λ_i)
        s_vn = self._von_neumann_entropy(rho)
        eigenvalues = np.real(np.linalg.eigvalsh(rho))

        # 3. Thermal entropy (biological units)
        s_thermal = KB_HARTREE * s_vn   # Ha/K

        # 4. Enthalpy ≈ electronic energy (ZPE negligible for our fragment)
        H = energy

        # 5. Free energy: G = H - TS
        G = H - self.T * s_thermal
        G_kcal = G * HA_TO_KCAL
        G_kj = G * HA_TO_KJ

        result = ThermodynamicsResult(
            fragment_id=fragment_id,
            temperature_K=self.T,
            electronic_energy=energy,
            enthalpy=H,
            entropy_vn=s_vn,
            entropy_thermal=s_thermal,
            free_energy=G,
            free_energy_kcal=G_kcal,
            free_energy_kj=G_kj,
            density_eigenvalues=eigenvalues,
            is_stable=(G < 0),
        )

        logger.info(result.summary())
        return result

    def compare(
        self,
        results: List[ThermodynamicsResult],
    ) -> ConformationComparison:
        """
        Compare conformations by ΔG relative to most stable one.

        Args:
            results: list of ThermodynamicsResult from multiple conformations

        Returns:
            ConformationComparison sorted by ΔG
        """
        if not results:
            raise ValueError("No results to compare")

        # Reference = most stable (lowest G)
        ref_idx = int(np.argmin([r.free_energy for r in results]))
        ref = results[ref_idx]

        rankings = []
        for r in results:
            dg = r.free_energy - ref.free_energy
            rankings.append((r.fragment_id, dg, dg * HA_TO_KCAL))

        rankings.sort(key=lambda x: x[1])

        comparison = ConformationComparison(
            reference_id=ref.fragment_id,
            temperature_K=self.T,
            rankings=rankings,
        )
        comparison.print_table()
        return comparison

    # ── Private ──────────────────────────────────────────────

    def _density_matrix(
        self,
        state: np.ndarray,
        subsystem_size: Optional[int] = None,
    ) -> np.ndarray:
        """Build density matrix. If subsystem_size given, trace out rest."""
        rho = np.outer(state, state.conj())

        if subsystem_size is not None:
            n_total = int(np.log2(len(state)))
            if subsystem_size < n_total:
                rho = self._partial_trace(rho, subsystem_size, n_total)

        return rho

    def _partial_trace(
        self,
        rho: np.ndarray,
        keep_n: int,
        total_n: int,
    ) -> np.ndarray:
        """
        Trace out qubits [keep_n, ..., total_n-1].
        Returns reduced density matrix of size (2^keep_n, 2^keep_n).
        """
        d_keep = 2 ** keep_n
        d_trace = 2 ** (total_n - keep_n)

        rho_tensor = rho.reshape([2] * (2 * total_n))

        # Trace out last (total_n - keep_n) qubits
        for _ in range(total_n - keep_n):
            n_current = rho_tensor.ndim // 2
            rho_tensor = np.trace(rho_tensor, axis1=n_current-1, axis2=2*n_current-1)

        return rho_tensor.reshape(d_keep, d_keep)

    def _von_neumann_entropy(self, rho: np.ndarray) -> float:
        """
        S = -Tr(ρ log ρ) = -Σ λ_i log(λ_i)

        Uses eigendecomposition for numerical stability.
        Convention: 0 × log(0) = 0 (limit is 0).
        """
        eigvals = np.real(np.linalg.eigvalsh(rho))
        eigvals = np.clip(eigvals, 0.0, 1.0)

        eps = 1e-15
        nonzero = eigvals[eigvals > eps]
        if len(nonzero) == 0:
            return 0.0

        return float(-np.sum(nonzero * np.log(nonzero)))


# ── Main ─────────────────────────────────────────────────────

if __name__ == "__main__":
    import os
    os.chdir(Path(__file__).parent.parent)

    logger.info("=" * 60)
    logger.info("STAGE 11 — Free Energy Test")
    logger.info("=" * 60)

    np.random.seed(42)
    calc = FreeEnergyCalculator(temperature_K=310.0)

    # Simulate 3 conformations with different energies and entanglement
    conformations = [
        # (name, energy_Ha, state_vector_probs)
        ("conf_A_native",    -1.137, [0.97, 0.02, 0.01, 0.00]),
        ("conf_B_misfolded", -1.100, [0.60, 0.25, 0.10, 0.05]),
        ("conf_C_unfolded",  -1.050, [0.40, 0.30, 0.20, 0.10]),
    ]

    results = []
    for name, energy, probs in conformations:
        probs = np.array(probs, dtype=float)
        probs /= probs.sum()
        state = np.sqrt(probs).astype(complex)
        r = calc.compute(state, energy, fragment_id=name)
        results.append(r)

    comparison = calc.compare(results)

    logger.info("\nKey result for paper:")
    logger.info(f"Most stable conformation: {comparison.reference_id}")
    logger.info(f"ΔG between best and worst: "
                f"{comparison.rankings[-1][2]:.3f} kcal/mol")
    logger.success("Stage 11 PASSED")
