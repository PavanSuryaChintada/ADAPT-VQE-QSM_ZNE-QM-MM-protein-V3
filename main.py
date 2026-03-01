"""
main.py — Hybrid Quantum-Classical Protein Fragment Solver
===========================================================
Runs all 12 pipeline stages end-to-end.

Usage:
  python main.py --pdb 2I9M --fragment 5 --qubits 8
  python main.py --pdb 2JOF --electrons 4 --orbitals 4
  python main.py --benchmark --proteins 2I9M 2JOF

Recommended starting point:
  python main.py --pdb 2I9M --fragment 5 --qubits 8 --electrons 4 --orbitals 4
"""

import sys
import os
import time
import json
import argparse
import numpy as np
from pathlib import Path
from datetime import datetime
from dataclasses import dataclass, asdict, field
from typing import Optional, List, Dict

# Project root on path
ROOT = Path(__file__).parent
sys.path.insert(0, str(ROOT))

from utils import logger


# ─────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────

@dataclass
class Config:
    # Protein
    pdb_id: str = "2I9M"
    center_residue: int = 5
    n_neighbors: int = 2

    # Quantum
    max_qubits: int = 8
    active_electrons: int = 4
    active_orbitals: int = 4
    basis: str = "sto-3g"
    shots: int = 1024

    # ADAPT-VQE
    max_adapt_iterations: int = 30
    gradient_threshold: float = 1e-3
    convergence_threshold: float = 1e-6

    # QSE
    n_qse_operators: int = 4

    # Thermodynamics
    temperature_K: float = 310.0

    # Output
    save_results: bool = True
    results_dir: str = "outputs/energies"
    random_seed: int = 42

    def validate(self):
        n_q = 2 * self.active_orbitals
        if n_q > self.max_qubits:
            raise ValueError(
                f"Active space ({self.active_electrons}e, {self.active_orbitals}o) "
                f"requires {n_q} qubits > max_qubits={self.max_qubits}.\n"
                f"Either increase --qubits or decrease --orbitals."
            )
        return self


# ─────────────────────────────────────────────────────────────
# Result container
# ─────────────────────────────────────────────────────────────

@dataclass
class PipelineResult:
    pdb_id: str
    fragment_id: str
    n_qubits_used: int
    n_atoms_qm: int
    n_atoms_mm: int

    hf_energy: float = 0.0
    vqe_energy: float = 0.0
    qse_energy: float = 0.0
    mitigated_energy: float = 0.0
    free_energy: float = 0.0
    free_energy_kcal: float = 0.0

    vqe_converged: bool = False
    vqe_n_params: int = 0
    vqe_n_iterations: int = 0
    qse_correction: float = 0.0
    noise_correction: float = 0.0

    runtime_seconds: float = 0.0
    timestamp: str = ""
    stage_times: Dict[str, float] = field(default_factory=dict)


# ─────────────────────────────────────────────────────────────
# Pipeline
# ─────────────────────────────────────────────────────────────

class QuantumProteinPipeline:
    """
    Runs all 12 stages and produces a complete PipelineResult.

    Each stage is a separate method for clarity and testability.
    """

    def __init__(self, config: Config):
        self.config = config
        np.random.seed(config.random_seed)
        self._stage_times = {}

    def run(self) -> PipelineResult:
        t0 = time.time()
        cfg = self.config
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")

        self._header()

        # ─── Stage 1: Data ──────────────────────────────────
        structure, qm_frag, mm_atoms, mm_coords, mm_charges = self._stage1()

        # ─── Stage 2: YOPO ──────────────────────────────────
        features = self._stage2(qm_frag)

        # ─── Stage 3: DISCA ─────────────────────────────────
        cluster_result = self._stage3(features)

        # ─── Stage 4+5: Hamiltonian ─────────────────────────
        h_data = self._stage45(qm_frag, mm_coords, mm_charges)

        # ─── Stage 6: Mapping ───────────────────────────────
        qubit_h = self._stage6(h_data)

        # ─── Stage 7: ADAPT-VQE ─────────────────────────────
        vqe_result = self._stage7(qubit_h, h_data)

        # ─── Stage 8+9: NatGrad + QSE ───────────────────────
        qse_result = self._stage89(vqe_result, qubit_h)

        # ─── Stage 10: Noise ────────────────────────────────
        noise_result = self._stage10(qse_result)

        # ─── Stage 11: Thermodynamics ───────────────────────
        thermo = self._stage11(noise_result, vqe_result)

        # ─── Stage 12: Reassembly ───────────────────────────
        final_energy = self._stage12(qse_result, noise_result)

        # ─── Compile result ─────────────────────────────────
        runtime = time.time() - t0

        result = PipelineResult(
            pdb_id=cfg.pdb_id,
            fragment_id=f"{cfg.pdb_id}_res{cfg.center_residue}",
            n_qubits_used=qubit_h.n_qubits,
            n_atoms_qm=qm_frag.n_atoms,
            n_atoms_mm=len(mm_atoms),
            hf_energy=h_data.hf_energy,
            vqe_energy=vqe_result.energy,
            qse_energy=qse_result.gs_energy,
            mitigated_energy=noise_result.mitigated_energy,
            free_energy=thermo.free_energy,
            free_energy_kcal=thermo.free_energy_kcal,
            vqe_converged=vqe_result.converged,
            vqe_n_params=vqe_result.n_parameters,
            vqe_n_iterations=vqe_result.n_iterations,
            qse_correction=qse_result.correction,
            noise_correction=noise_result.noise_correction,
            runtime_seconds=runtime,
            timestamp=ts,
            stage_times=self._stage_times,
        )

        self._print_results(result)
        if cfg.save_results:
            self._save(result)

        return result

    # ─── Stage implementations ───────────────────────────────

    def _stage1(self):
        logger.info("\n[STAGE 1] Data Acquisition")
        t = time.time()
        from core.data_acquisition import PDBDownloader, FragmentExtractor

        dl = PDBDownloader(data_dir="data/raw")
        structure = dl.fetch(self.config.pdb_id)

        extractor = FragmentExtractor(max_qm_atoms=15)
        qm_frag, mm_atoms, mm_coords, mm_charges = extractor.extract(
            structure,
            center_residue=self.config.center_residue,
            n_neighbors=self.config.n_neighbors,
        )
        self._stage_times["stage1"] = time.time() - t
        return structure, qm_frag, mm_atoms, mm_coords, mm_charges

    def _stage2(self, qm_frag):
        logger.info("\n[STAGE 2] YOPO Feature Extraction")
        t = time.time()
        from core.yopo.feature_extractor import YOPOExtractor

        coords = qm_frag.get_coords()
        ex = YOPOExtractor(n_eigenvalues=10)
        features = ex.extract(coords, pdb_id=self.config.pdb_id)
        logger.info(f"Feature vector: {len(features.feature_vector)}-dim | "
                    f"Rg={features.radius_of_gyration:.2f} Å")
        self._stage_times["stage2"] = time.time() - t
        return features

    def _stage3(self, features):
        logger.info("\n[STAGE 3] DISCA Clustering")
        t = time.time()
        from core.disca.clustering import DISCAClustering

        # For single structure: simulate small ensemble with noise
        rng = np.random.default_rng(self.config.random_seed)
        fv = features.feature_vector
        ensemble = np.vstack([
            fv + rng.normal(0, 0.05, len(fv)) for _ in range(20)
        ])

        clusterer = DISCAClustering(n_components=2, purity_threshold=0.70)
        result = clusterer.fit(ensemble)
        self._stage_times["stage3"] = time.time() - t
        return result

    def _stage45(self, qm_frag, mm_coords, mm_charges):
        logger.info("\n[STAGE 4+5] QM/MM + Hamiltonian Construction")
        t = time.time()
        from quantum.hamiltonian.builder import QMMMHamiltonianBuilder, FragmentExtractor as FE

        builder = QMMMHamiltonianBuilder(
            basis=self.config.basis,
            active_electrons=self.config.active_electrons,
            active_orbitals=self.config.active_orbitals,
            max_qubits=self.config.max_qubits,
        )

        from core.data_acquisition import FragmentExtractor
        fex = FragmentExtractor()
        atoms = fex.atoms_to_pyscf_list(qm_frag.atoms)

        h_data = builder.build(atoms, mm_coords, mm_charges)
        logger.info(h_data.summary())
        self._stage_times["stage45"] = time.time() - t
        return h_data

    def _stage6(self, h_data):
        logger.info("\n[STAGE 6] Qubit Mapping (Parity)")
        t = time.time()
        from quantum.mapping.qubit_mapper import QubitMapper

        mapper = QubitMapper(method="parity", two_qubit_reduction=True)
        qh = mapper.map(h_data)
        logger.info(qh.summary())
        self._stage_times["stage6"] = time.time() - t
        return qh

    def _stage7(self, qubit_h, h_data):
        logger.info("\n[STAGE 7] ADAPT-VQE")
        t = time.time()
        from quantum.vqe.adapt_vqe import ADAPTVQESolver

        solver = ADAPTVQESolver(
            qubit_hamiltonian=qubit_h,
            n_electrons=self.config.active_electrons,
            max_iterations=self.config.max_adapt_iterations,
            gradient_threshold=self.config.gradient_threshold,
            convergence_threshold=self.config.convergence_threshold,
        )
        result = solver.solve(hf_energy=h_data.hf_energy)
        logger.info(result.summary())
        self._stage_times["stage7"] = time.time() - t
        return result

    def _stage89(self, vqe_result, qubit_h):
        logger.info("\n[STAGE 8+9] Natural Gradient + QSE")
        t = time.time()
        from quantum.optimization.natural_gradient import QSERefinement

        H_mat = qubit_h.to_matrix()
        qse = QSERefinement(n_operators=self.config.n_qse_operators)
        result = qse.refine(
            vqe_result, H_mat,
            n_qubits=qubit_h.n_qubits,
            n_electrons=self.config.active_electrons,
        )
        logger.info(result.summary())
        self._stage_times["stage89"] = time.time() - t
        return result

    def _stage10(self, qse_result):
        logger.info("\n[STAGE 10] Noise Mitigation (ZNE + Pauli Twirling)")
        t = time.time()
        from quantum.noise.mitigation import NoiseMitigationPipeline

        pipeline = NoiseMitigationPipeline(
            n_electrons=self.config.active_electrons,
            scale_factors=[1.0, 2.0, 3.0],
        )
        result = pipeline.mitigate(
            raw_energy=qse_result.gs_energy,
            noise_strength=0.003,   # realistic NISQ noise
            n_gates=2 * self.config.active_orbitals ** 2,
            seed=self.config.random_seed,
        )
        self._stage_times["stage10"] = time.time() - t
        return result

    def _stage11(self, noise_result, vqe_result):
        logger.info("\n[STAGE 11] Free Energy Computation")
        t = time.time()
        from thermodynamics.free_energy import FreeEnergyCalculator

        calc = FreeEnergyCalculator(temperature_K=self.config.temperature_K)

        # Use VQE state vector
        state = vqe_result.state_vector
        result = calc.compute(
            state_vector=state,
            energy=noise_result.mitigated_energy,
            fragment_id=f"{self.config.pdb_id}_res{self.config.center_residue}",
        )
        self._stage_times["stage11"] = time.time() - t
        return result

    def _stage12(self, qse_result, noise_result):
        logger.info("\n[STAGE 12] Global Reassembly")
        # Combine corrections additively
        final = noise_result.mitigated_energy
        logger.info(f"Final assembled energy: {final:.8f} Ha")
        return final

    # ─── Utilities ───────────────────────────────────────────

    def _header(self):
        cfg = self.config
        logger.info("\n" + "=" * 70)
        logger.info("  HYBRID QUANTUM-CLASSICAL PROTEIN FRAGMENT SOLVER")
        logger.info("=" * 70)
        logger.info(f"  PDB:              {cfg.pdb_id}")
        logger.info(f"  Center residue:   {cfg.center_residue}")
        logger.info(f"  Active space:     ({cfg.active_electrons}e, {cfg.active_orbitals}o)")
        logger.info(f"  Max qubits:       {cfg.max_qubits}")
        logger.info(f"  Temperature:      {cfg.temperature_K}K")
        logger.info(f"  Basis set:        {cfg.basis}")
        logger.info("=" * 70)

    def _print_results(self, r: PipelineResult):
        logger.info("\n" + "=" * 70)
        logger.info("  RESULTS SUMMARY")
        logger.info("=" * 70)
        logger.info(f"  Fragment:         {r.fragment_id}")
        logger.info(f"  QM atoms:         {r.n_atoms_qm}  |  MM atoms: {r.n_atoms_mm}")
        logger.info(f"  Qubits used:      {r.n_qubits_used}")
        logger.info("")
        logger.info(f"  HF energy:        {r.hf_energy:+.8f} Ha   (classical reference)")
        logger.info(f"  VQE energy:       {r.vqe_energy:+.8f} Ha   "
                    f"({'converged' if r.vqe_converged else 'NOT converged'}, "
                    f"{r.vqe_n_params} params, {r.vqe_n_iterations} iters)")
        logger.info(f"  QSE energy:       {r.qse_energy:+.8f} Ha   "
                    f"(correction: {r.qse_correction:+.6f} Ha)")
        logger.info(f"  Mitigated energy: {r.mitigated_energy:+.8f} Ha   "
                    f"(noise correction: {r.noise_correction:+.6f} Ha)")
        logger.info(f"  Free energy G:    {r.free_energy:+.8f} Ha   "
                    f"= {r.free_energy_kcal:+.4f} kcal/mol")
        logger.info("")
        logger.info(f"  Correlation energy: {r.vqe_energy - r.hf_energy:+.8f} Ha")
        logger.info(f"  Total runtime:    {r.runtime_seconds:.1f} s")
        logger.info("=" * 70)

    def _save(self, r: PipelineResult):
        outdir = Path(self.config.results_dir)
        outdir.mkdir(parents=True, exist_ok=True)
        path = outdir / f"result_{r.pdb_id}_{r.timestamp}.json"
        with open(path, "w") as f:
            json.dump(asdict(r), f, indent=2, default=str)
        logger.success(f"Results saved: {path}")


# ─────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Hybrid Quantum-Classical Protein Fragment Solver",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--pdb",        default="2I9M",  help="PDB ID  (2I9M, 2JOF, 1VII, 1CRN, 1UBQ)")
    p.add_argument("--fragment",   type=int, default=5, help="Center residue")
    p.add_argument("--qubits",     type=int, default=8,  help="Max qubits (≤14 safe for 16GB)")
    p.add_argument("--electrons",  type=int, default=4,  help="Active electrons")
    p.add_argument("--orbitals",   type=int, default=4,  help="Active orbitals")
    p.add_argument("--shots",      type=int, default=1024)
    p.add_argument("--temp",       type=float, default=310.0, help="Temperature (K)")
    p.add_argument("--basis",      default="sto-3g")
    p.add_argument("--seed",       type=int, default=42)
    p.add_argument("--no-save",    action="store_true")

    # Benchmark mode
    p.add_argument("--benchmark",  action="store_true", help="Run benchmark vs classical")
    p.add_argument("--proteins",   nargs="+", default=["2I9M", "2JOF"],
                   help="Proteins for benchmark mode")
    return p.parse_args()


def run_benchmark(proteins: List[str], config: Config):
    """Run pipeline on multiple proteins and compare methods."""
    from scripts.benchmark import run_full_benchmark
    run_full_benchmark(proteins, config)


if __name__ == "__main__":
    args = parse_args()

    cfg = Config(
        pdb_id=args.pdb,
        center_residue=args.fragment,
        max_qubits=args.qubits,
        active_electrons=args.electrons,
        active_orbitals=args.orbitals,
        shots=args.shots,
        temperature_K=args.temp,
        basis=args.basis,
        random_seed=args.seed,
        save_results=not args.no_save,
    )

    try:
        cfg.validate()
    except ValueError as e:
        logger.error(str(e))
        sys.exit(1)

    if args.benchmark:
        run_benchmark(args.proteins, cfg)
    else:
        pipeline = QuantumProteinPipeline(cfg)
        result = pipeline.run()
