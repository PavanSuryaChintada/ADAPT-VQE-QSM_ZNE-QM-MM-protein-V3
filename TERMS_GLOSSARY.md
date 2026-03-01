# Terms Glossary — Definitions and Explanations

Alphabetical reference for terms used in the pipeline, circuits, and outputs.

---

## A

**Active space** — A subset of orbitals and electrons used for quantum simulation (e.g. 4 electrons, 4 orbitals). Larger active spaces require more qubits and gates.

**ADAPT-VQE** — *Adaptive* Variational Quantum Eigensolver. Builds the ansatz iteratively by adding excitation operators with the largest gradients instead of a fixed UCCSD ansatz.

**Ansatz** — The variational quantum circuit **U(θ)** that prepares the trial state |ψ(θ)⟩. Parameters θ are optimized to minimize ⟨ψ|H|ψ⟩.

**AO2MO** — Atomic orbitals to molecular orbitals. Transforms one- and two-electron integrals to the MO basis.

---

## B

**Barren plateau** — A region where the energy gradient becomes very small, making optimization hard. ADAPT-VQE helps avoid this by using gradient-based operator selection.

**Basis set** — Set of functions used to represent molecular orbitals (e.g. STO-3G). Smaller basis sets mean fewer orbitals and qubits.

**BIC** — Bayesian Information Criterion. Used in DISCA clustering to choose the number of components.

**Compactness** — 1 / (1 + within-cluster variance). Higher = tighter cluster = more representative conformation. Used with purity to select the best cluster.

---

## C

**Chemical accuracy** — Error threshold of ~1 kcal/mol (≈ 0.0016 Ha). Widely used standard for molecular energies.

**CNOT** — Controlled-NOT gate. Flips the target qubit if the control qubit is |1⟩. Used in ADAPT-VQE for excitation operators.

**Correlation energy** — E_corr = E_exact − E_HF. Electron correlation beyond Hartree–Fock.

**Circuit depth** — Number of gate layers. Lower depth is better for noisy hardware.

---

## D

**DISCA** — Clustering method (GMM) for conformational states. Used for preprocessing structural ensembles.

**Double excitation** — Operator that moves two electrons between orbitals (e.g. occupied → virtual).

---

## E

**Error (Ha)** — |E_method − E_exact|. Distance from the exact (FCI) energy.

**Exact diagonalization** — Direct computation of the ground state by diagonalizing the full Hamiltonian. Exact but exponentially expensive.

**Excitation** — Moving one or more electrons from occupied to virtual orbitals. Basis for building correlation in VQE.

**Expectation value** — ⟨ψ|H|ψ⟩. Energy of the state |ψ⟩ under Hamiltonian H.

---

## F

**FCI** — Full Configuration Interaction. Exact solution in the given basis set. Used as reference in benchmarks.

**Feature vector** — Fixed-size vector (e.g. 34-dim from YOPO) used to represent structure or conformation.

**Fragment** — Region of the protein treated as the QM region. E.g. center residue ± neighbors.

**Free energy G** — ΔG = ΔH − TΔS at 310 K. Negative G = thermodynamically stable conformation.

---

## G

**Gradient** — ∂E/∂θ. Used to select ADAPT operators and update parameters.

**GMM** — Gaussian Mixture Model. Used for DISCA clustering.

**Ground state** — Lowest energy eigenstate of the Hamiltonian.

---

## H

**Ha (Hartree)** — Atomic unit of energy. 1 Ha ≈ 627.5 kcal/mol. Energy scale used in quantum chemistry.

**Hamiltonian (H)** — Operator for total energy. In second quantization: H = Σ h_pq a†_p a_q + (1/2) Σ h_pqrs a†_p a†_q a_r a_s + E_core.

**Hartree–Fock (HF)** — Mean-field method. No electron correlation. Used as starting state and reference energy.

**HF state |HF⟩** — Reference state with electrons in lowest orbitals. |HF⟩ = |111...000⟩ in qubit notation.

---

## J

**Jordan–Wigner (JW)** — Fermion-to-qubit mapping. Uses n qubits for n spin-orbitals, with long-range Z strings.

---

## M

**Mapping** — Fermion → qubit encoding (e.g. Jordan–Wigner, parity).

**Mitigated energy** — Energy after Zero Noise Extrapolation (ZNE) to reduce hardware noise.

**MM (molecular mechanics)** — Classical region (e.g. point charges) surrounding the QM region.

---

## N

**Natural gradient** — Optimization using the quantum Fisher information matrix (QFIM) instead of plain gradients.

**NISQ** — Noisy Intermediate-Scale Quantum. Current quantum hardware regime.

**Noise mitigation** — Techniques (e.g. ZNE, Pauli twirling) to reduce effects of hardware noise.

**Noise correction** — Energy shift from ZNE. Positive: raw energy was too low; negative: raw energy was too high.

---

## O

**Operator pool** — Set of excitation operators (single, double) from which ADAPT-VQE picks.

**Orbital** — Molecular orbital. In the active space, each orbital has 2 spin-orbitals (α and β).

---

## P

**Parity mapping** — Fermion-to-qubit encoding that exploits Z₂ symmetry to remove 2 qubits vs Jordan–Wigner.

**Parameter-shift rule** — ∂E/∂θ = ½[E(θ + π/2) − E(θ − π/2)]. Used for exact gradients on quantum hardware.

**Pauli twirling** — Noise mitigation by averaging over Pauli rotations to symmetrize noise.

**PDB** — Protein Data Bank. Source of protein structures (e.g. 2I9M, 2JOF).

**Purity** — Mean responsibility (confidence) that members belong to a DISCA cluster. Higher = more coherent.

**Selection score** — purity × compactness. Best cluster = argmax(selection score). When purity ties, compactness breaks the tie (tighter cluster wins).

**PySCF** — Python library for quantum chemistry. Builds Hamiltonians from molecular geometry.

---

## Q

**QFIM** — Quantum Fisher Information Matrix. Used in natural gradient optimization.

**QM (quantum mechanics)** — Region treated with quantum chemistry (e.g. fragment atoms).

**QM/MM** — Hybrid model: QM region + MM environment.

**QSE** — Quantum Subspace Expansion. Expands the variational space after VQE to improve the energy estimate.

**Qubit Hamiltonian** — H expressed as a sum of Pauli strings: H = Σ c_i P_i.

---

## R

**Rg (radius of gyration)** — Measure of protein compactness in Å.

**Richardson extrapolation** — Extrapolates to zero noise from runs at different noise scales (ZNE).

**RY gate** — Rotation around the Y-axis. RY(2θ) implements single excitations in ADAPT-VQE.

---

## S

**Second quantization** — Formulation using creation (a†) and annihilation (a) operators.

**Shots** — Number of measurements per Pauli term. More shots → lower variance.

**Single excitation** — Operator that moves one electron between orbitals.

**Spin-orbital** — Orbital with spin (α or β). n orbitals → 2n spin-orbitals.

**STO-3G** — Minimal basis set. Good for small systems and prototyping.

**STABLE / UNSTABLE** — Sign of ΔG: negative = STABLE, positive = UNSTABLE.

---

## T

**Thermodynamics** — Computation of ΔG at 310 K for stability.

**Two-qubit reduction** — Parity mapping removes 2 qubits using particle-number symmetry.

---

## U

**UCCSD** — Unitary Coupled Cluster Singles and Doubles. Fixed ansatz with all single and double excitations.

---

## V

**Variational** — Optimize parameters to minimize energy. No need to prepare the exact eigenstate.

**VQE** — Variational Quantum Eigensolver. Uses a quantum computer to prepare |ψ(θ)⟩ and a classical optimizer to minimize ⟨ψ|H|ψ⟩.

**VQE energy** — ⟨ψ(θ)|H|ψ(θ)⟩ at optimized θ.

---

## Y

**YOPO** — Rotation-invariant structural fingerprint. Used for conformational analysis.

---

## Z

**ZNE** — Zero Noise Extrapolation. Runs at different noise scales and extrapolates to zero noise.

**ZNE correction** — Energy shift from extrapolation: E_mitigated − E_raw.

---

## Units Quick Reference

| Unit | Conversion |
|------|------------|
| 1 Ha | ≈ 627.5 kcal/mol |
| Chemical accuracy | 1 kcal/mol ≈ 0.0016 Ha |
| Temperature | 310 K ≈ 37 °C (physiological) |

---

## Related

- **Outputs:** [OUTPUT_GUIDE.md](OUTPUT_GUIDE.md)
- **Circuits:** [CIRCUIT_GUIDE.md](CIRCUIT_GUIDE.md)
