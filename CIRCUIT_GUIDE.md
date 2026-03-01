# Circuit Guide — Understanding the Quantum Circuit

Explanation of the quantum circuit produced by `show_circuit.py` and how it relates to the pipeline output.

---

## How to Generate the Circuit Output

```cmd
python scripts\show_circuit.py           # Default (Chignolin fragment)
python scripts\show_circuit.py --full    # Add Hamiltonian terms
python scripts\show_circuit.py --molecule H2 --electrons 2 --orbitals 2   # Simpler molecule
```

---

## 1. Circuit Overview

The ADAPT-VQE circuit has three main stages:

```
|0...0⟩  →  [HF preparation]  →  [ADAPT ansatz]  →  [Measurement]  →  ⟨ψ|H|ψ⟩
```

| Stage | Role | What happens |
|-------|------|--------------|
| **HF preparation** | Initialize | Put electrons into the Hartree–Fock reference state |
| **ADAPT ansatz** | Variational | Apply excitation operators chosen by gradient |
| **Measurement** | Energy | Measure Pauli terms to get ⟨ψ\|H\|ψ⟩ |

---

## 2. Full Circuit Diagram

Example for Chignolin (4e, 4o) after parity mapping → 6 qubits:

```
  q5 |0⟩ ───────────────[M]
  q4 |0⟩ ───────────────[M]
  q3 |0⟩ ──────[RY]───[M]
  q2 |0⟩ ──[X]──[Z]─────[M]
  q1 |0⟩ ──[X]──[●]────[M]
  q0 |0⟩ ──[X]───────────[M]
```

### Gate Legend

| Symbol | Gate | Meaning |
|--------|------|---------|
| **[X]** | Pauli-X | Flip \|0⟩ → \|1⟩ (occupy orbital) |
| **[●]** | Control (CNOT) | Part of excitation operator |
| **[Z]** | Pauli-Z | Jordan–Wigner phase string |
| **[RY]** | RY(2θ) | Rotation on target qubit (excitation angle) |
| **[M]** | Measurement | Read-out in Z basis |

---

## 3. Stage 1: Hartree–Fock Preparation

**Purpose:** Start from the HF reference state (electrons in lowest orbitals).

**For 3 electrons in 6 orbitals (after parity):**

```
|HF⟩ = |111000⟩  (q0,q1,q2 = |1⟩, q3,q4,q5 = |0⟩)
```

Each occupied orbital gets a Pauli-X on its qubit:

```
   q0  |0⟩ ───┤ X ├── |1⟩  ← occupied
   q1  |0⟩ ───┤ X ├── |1⟩  ← occupied
   q2  |0⟩ ───┤ X ├── |1⟩  ← occupied
   q3  |0⟩ ─────────────── |0⟩  ← virtual
   q4  |0⟩ ─────────────── |0⟩
   q5  |0⟩ ─────────────── |0⟩
```

---

## 4. Stage 2: ADAPT-VQE Ansatz

**Purpose:** Apply excitations chosen by the gradient.

Form: **U(θ) = ∏ exp(θₖ Aₖ)**, where each **Aₖ** is a single or double excitation operator.

### Single Excitation Example

For excitation from orbital 1 → 3:

```
   q0  ─────────────────────
   q1  ──────── ●  ────────   ← control (occupied)
   q2  ───────── Z ─────────   ← Jordan–Wigner phase
   q3  ────────┤RY├────────   ← rotation angle θ
   q4  ─────────────────────
   q5  ─────────────────────
```

Steps:

1. CNOT ladder for the Z string.
2. **RY(2θ)** on the virtual orbital.
3. CNOT ladder (uncompute).

### Why ADAPT-VQE Is Compact

- Only operators with large gradient are added.
- Result: typically 1–5 operators (vs 15+ in full UCCSD).
- Helps avoid barren plateaus in fixed ansätze.

---

## 5. Stage 3: Pauli Measurement

**Purpose:** Estimate **E(θ) = ⟨ψ(θ)|H|ψ(θ)⟩** by measuring Pauli strings.

```
E(θ) = Σᵢ cᵢ ⟨ψ(θ)|Pᵢ|ψ(θ)⟩
```

- **H** is written as a sum of Pauli operators **Pᵢ** with coefficients **cᵢ**.
- Each **Pᵢ** is estimated from measurement outcomes.
- Energy is the weighted sum.

### Parameter-Shift Rule (for gradients)

```
∂E/∂θₖ = ½ [ E(θₖ + π/2) − E(θₖ − π/2) ]
```

ADAPT-VQE uses these gradients to pick new operators and update angles.

---

## 6. Qubit Mapping: Parity vs Jordan–Wigner

| Mapping | Qubits for (4e, 4o) | Advantage |
|---------|----------------------|-----------|
| Jordan–Wigner | 8 | Simple, local terms |
| **Parity** | **6** | Uses Z₂ symmetry → 2 fewer qubits |

Parity mapping lets us drop 2 qubits while preserving physics.

---

## 7. Circuit Statistics (from `show_circuit.py`)

Example table:

```
  Gate type                   ADAPT-VQE (ours)    Standard UCCSD
  ------------------------- ------------------  ----------------
  Pauli-X (init)                             6                 6
  CNOT (2-qubit)                              4                 90
  RY rotations                                1                 24
  Pauli-Z (JW)                                1               N/A
  ------------------------- ------------------  ----------------
  Total gates                                12                180
  Circuit depth                               7                360
  Parameters (θ)                              1                 15
```

- **ADAPT-VQE** is much shorter and shallower than UCCSD.
- Fewer gates → less noise on real hardware.
- Fewer parameters → easier optimization.

---

## 8. ADAPT-VQE Convergence Table

```
  Iter     Energy (Ha)            ΔE      |∇|max  Status
  ----  --------------  ------------  ----------  --------------------
    HF     -7.86224600             —           —  Hartree–Fock start
     1      0.63213921             —    0.111725  optimizing...
     2      0.63213921   +0.00000000    0.111725  optimizing...
```

| Column | Meaning |
|--------|---------|
| **Iter** | ADAPT iteration (HF = initial) |
| **Energy** | Current VQE energy (Ha) |
| **ΔE** | Change in energy from previous iteration |
| **\|\|∇\|\|max** | Largest gradient component |
| **Status** | Convergence / optimization status |

Stopping rules:

- Gradient convergence: \|\|∇\|\| < threshold
- Energy convergence: ΔE < threshold

---

## 9. Hamiltonian Terms (with `--full`)

`show_circuit.py --full` prints the qubit Hamiltonian:

```
  H = Σᵢ cᵢ Pᵢ
```

Example:

```
  Identity terms:  7  (constant energy offset)
    +1.0507 × IIIIII
    +0.0745 × IIIIII
    ...

  Single-qubit terms:  6  (local fields)
    -0.0745 × ZIIIII
    -0.0745 × IZIIII
    ...

  Two-qubit terms:  0
  Higher-body terms:  12
    -0.0279 × XZXIII
    ...
```

These are the Pauli strings and coefficients used for measurement.

---

## 10. H₂ Example (simplest)

```cmd
python scripts\show_circuit.py --molecule H2 --electrons 2 --orbitals 2
```

- **2 electrons, 2 orbitals** → parity → **2 qubits**
- Circuit:

```
  q1 |0⟩ ──────[●]────[M]
  q0 |0⟩ ──[X]───────────[M]
```

- HF: \|10⟩ (one electron in lowest orbital)
- One excitation operator (angle θ = 0 at convergence).
- Good minimal example to inspect the circuit and pipeline.

---

## Summary

| Concept | Summary |
|---------|---------|
| **Stages** | HF prep → ADAPT ansatz → measurement |
| **HF prep** | Pauli-X on occupied orbitals |
| **ADAPT** | Gradient-selected excitations: CNOT + RY + CNOT |
| **Measurement** | Estimate Pauli expectations → energy |
| **Parity** | 2 fewer qubits than Jordan–Wigner |
| **Stats** | ADAPT-VQE is much shallower than UCCSD |

---

## Related

- **Outputs:** [OUTPUT_GUIDE.md](OUTPUT_GUIDE.md)
- **Terms:** [TERMS_GLOSSARY.md](TERMS_GLOSSARY.md)
