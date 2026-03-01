"""
Stage 1 — Data Acquisition
===========================
Downloads protein structures directly from RCSB Protein Data Bank.
No external biology library required — pure Python urllib.

Best proteins for your project (priority order):
  1. 2I9M — Chignolin        (10 res) — NISQ gold standard, start here
  2. 2JOF — Trp-cage         (20 res) — IBM quantum paper benchmark
  3. 1VII — Villin HP35      (35 res) — fast-folding, well studied
  4. 1CRN — Crambin          (46 res) — highest resolution X-ray
  5. 1UBQ — Ubiquitin        (76 res) — use fragments only

Why these proteins:
  - Small enough to extract 10-20 atom QM fragments
  - Well-studied: exact energies available for validation
  - Used in published quantum chemistry / NISQ papers
  - Clean PDB structures (no missing atoms, no solvent issues)
"""

import os
import sys
import urllib.request
import urllib.error
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, field

sys.path.insert(0, str(Path(__file__).parent.parent))
from utils import logger


# ─────────────────────────────────────────────────────────────
# Data Structures
# ─────────────────────────────────────────────────────────────

@dataclass
class Atom:
    index: int
    name: str
    element: str
    residue_name: str
    residue_id: int
    chain: str
    x: float
    y: float
    z: float

    @property
    def coords(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z], dtype=float)

    @property
    def is_backbone(self) -> bool:
        return self.name.strip() in {"N", "CA", "C", "O"}

    @property
    def is_heavy(self) -> bool:
        return self.element.strip().upper() != "H"


@dataclass
class ProteinStructure:
    pdb_id: str
    atoms: List[Atom] = field(default_factory=list)
    residues: Dict[int, List[Atom]] = field(default_factory=dict)
    filepath: str = ""

    @property
    def n_atoms(self) -> int:
        return len(self.atoms)

    @property
    def n_residues(self) -> int:
        return len(self.residues)

    def get_coords(self) -> np.ndarray:
        if not self.atoms:
            return np.zeros((0, 3))
        return np.array([a.coords for a in self.atoms])

    def get_backbone_atoms(self) -> List[Atom]:
        return [a for a in self.atoms if a.is_backbone]

    def get_heavy_atoms(self) -> List[Atom]:
        return [a for a in self.atoms if a.is_heavy]

    def get_fragment(self, start_res: int, end_res: int) -> "ProteinStructure":
        frag = ProteinStructure(pdb_id=f"{self.pdb_id}_frag{start_res}-{end_res}")
        for rid in range(start_res, end_res + 1):
            if rid in self.residues:
                for atom in self.residues[rid]:
                    frag.atoms.append(atom)
                frag.residues[rid] = self.residues[rid]
        return frag

    def summary(self) -> str:
        return (
            f"[{self.pdb_id}] "
            f"{self.n_residues} residues | "
            f"{self.n_atoms} atoms | "
            f"{len(self.get_heavy_atoms())} heavy atoms"
        )


# ─────────────────────────────────────────────────────────────
# PDB Downloader
# ─────────────────────────────────────────────────────────────

class PDBDownloader:
    """
    Downloads and parses PDB files from RCSB.
    No biopython required — pure Python.

    Usage:
        dl = PDBDownloader(data_dir="data/raw")
        structure = dl.fetch("2I9M")
    """

    RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
    RCSB_FALLBACK = "https://www.rcsb.org/structure/{pdb_id}"

    # Partial charges for common protein atoms (simple estimates)
    # For publication accuracy: use AMBER/CHARMM force field charges
    PARTIAL_CHARGES = {
        "N":  -0.470, "CA":  0.070, "C":   0.510,
        "O":  -0.510, "CB":  0.070, "S":  -0.090,
        "H":   0.310, "HA":  0.090, "HN":  0.340,
    }

    # Standard element-to-mass mapping
    ATOMIC_MASSES = {
        "H": 1.008, "C": 12.011, "N": 14.007,
        "O": 15.999, "S": 32.06, "P": 30.974,
        "FE": 55.845, "ZN": 65.38,
    }

    def __init__(self, data_dir: str = "data/raw"):
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)

    def fetch(self, pdb_id: str, force_download: bool = False) -> ProteinStructure:
        """Download (if needed) and parse a PDB structure."""
        pdb_id = pdb_id.upper().strip()
        filepath = self.data_dir / f"{pdb_id}.pdb"

        if not filepath.exists() or force_download:
            self._download(pdb_id, filepath)
        else:
            logger.info(f"Using cached: {filepath}")

        structure = self._parse_pdb(pdb_id, filepath)
        logger.info(f"Loaded: {structure.summary()}")
        return structure

    def fetch_benchmark_set(self) -> Dict[str, ProteinStructure]:
        """Download all 5 recommended benchmark proteins."""
        ids = ["2I9M", "2JOF", "1VII", "1CRN", "1UBQ"]
        results = {}
        logger.info("Downloading benchmark protein set...")
        for pid in ids:
            try:
                results[pid] = self.fetch(pid)
            except Exception as e:
                logger.error(f"Failed {pid}: {e}")
        logger.success(f"Fetched {len(results)}/{len(ids)} proteins")
        return results

    def _download(self, pdb_id: str, filepath: Path) -> None:
        url = self.RCSB_URL.format(pdb_id=pdb_id)
        logger.info(f"Downloading {pdb_id} from RCSB...")
        try:
            urllib.request.urlretrieve(url, filepath)
            size_kb = filepath.stat().st_size // 1024
            logger.success(f"Downloaded {pdb_id} ({size_kb} KB) → {filepath}")
        except urllib.error.URLError as e:
            raise ConnectionError(
                f"Could not download {pdb_id}.\n"
                f"Check internet connection or manually download from:\n"
                f"  https://www.rcsb.org/structure/{pdb_id}\n"
                f"  Save as: {filepath}\n"
                f"Error: {e}"
            )

    def _parse_pdb(self, pdb_id: str, filepath: Path) -> ProteinStructure:
        """Parse ATOM records from PDB file.
        For NMR multi-model files, uses only the first MODEL to avoid OOM.
        """
        structure = ProteinStructure(pdb_id=pdb_id, filepath=str(filepath))

        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                record = line[:6].strip()
                # Stop after first model in NMR structures (avoids 20× atoms)
                if record == "ENDMDL":
                    break
                if record not in ("ATOM", "HETATM"):
                    continue

                try:
                    atom_idx   = int(line[6:11].strip())
                    atom_name  = line[12:16].strip()
                    res_name   = line[17:20].strip()
                    chain      = line[21].strip() or "A"
                    res_id     = int(line[22:26].strip())
                    x          = float(line[30:38].strip())
                    y          = float(line[38:46].strip())
                    z          = float(line[46:54].strip())

                    # Element: column 77-78 if present, else infer from atom name
                    raw_elem = line[76:78].strip() if len(line) > 76 else ""
                    if not raw_elem:
                        raw_elem = "".join(c for c in atom_name if c.isalpha())[:1]
                    element = raw_elem.upper()

                    # Skip solvent HETATM (HOH = water)
                    if record == "HETATM" and res_name in ("HOH", "WAT", "SOL"):
                        continue

                    atom = Atom(
                        index=atom_idx,
                        name=atom_name,
                        element=element,
                        residue_name=res_name,
                        residue_id=res_id,
                        chain=chain,
                        x=x, y=y, z=z,
                    )
                    structure.atoms.append(atom)
                    structure.residues.setdefault(res_id, []).append(atom)

                except (ValueError, IndexError):
                    continue  # skip malformed lines

        if structure.n_atoms == 0:
            raise ValueError(f"No ATOM records found in {filepath}")

        return structure


# ─────────────────────────────────────────────────────────────
# Fragment Extractor
# ─────────────────────────────────────────────────────────────

class FragmentExtractor:
    """
    Extracts QM/MM fragment from a protein structure.

    Strategy:
      1. Select center residue + N neighbors on each side
      2. Keep only backbone atoms (N, CA, C, O) to stay within qubit limit
      3. Cap dangling bonds (simplified: just track boundary atoms)
      4. Everything outside QM region = MM environment (electrostatics)

    Rule: Stay at max_qm_atoms ≤ 15 for safe 8-12 qubit simulation.
    """

    # Residue one-letter codes for logging
    AA_CODES = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }

    def __init__(self, max_qm_atoms: int = 15):
        self.max_qm_atoms = max_qm_atoms

    def extract(
        self,
        structure: ProteinStructure,
        center_residue: int,
        n_neighbors: int = 2,
    ) -> Tuple[ProteinStructure, List[Atom], np.ndarray, np.ndarray]:
        """
        Extract QM region and MM environment.

        Args:
            structure: full protein
            center_residue: residue ID at center of QM region
            n_neighbors: residues on each side

        Returns:
            (qm_fragment, mm_atoms, mm_coords, mm_charges)
        """
        all_res = sorted(structure.residues.keys())
        if center_residue not in all_res:
            center_residue = all_res[len(all_res) // 2]
            logger.warning(f"Residue not found. Using center: {center_residue}")

        start = max(all_res[0], center_residue - n_neighbors)
        end   = min(all_res[-1], center_residue + n_neighbors)

        qm_fragment = structure.get_fragment(start, end)

        # Trim to backbone only if too many atoms
        heavy = qm_fragment.get_heavy_atoms()
        if len(heavy) > self.max_qm_atoms:
            backbone = qm_fragment.get_backbone_atoms()
            qm_fragment.atoms = backbone
            qm_fragment.residues = {}
            for a in backbone:
                qm_fragment.residues.setdefault(a.residue_id, []).append(a)
            logger.warning(
                f"QM region trimmed to backbone only: "
                f"{len(backbone)} atoms (was {len(heavy)} heavy atoms)"
            )

        # MM region = all atoms NOT in QM fragment
        qm_indices = {a.index for a in qm_fragment.atoms}
        mm_atoms = [a for a in structure.atoms if a.index not in qm_indices]

        # MM coordinates and charges
        if mm_atoms:
            mm_coords  = np.array([a.coords for a in mm_atoms])
            mm_charges = np.array([
                PDBDownloader.PARTIAL_CHARGES.get(a.name, 0.0)
                for a in mm_atoms
            ])
        else:
            mm_coords  = np.zeros((0, 3))
            mm_charges = np.zeros(0)

        # Log the sequence of the QM fragment
        seq = "".join(
            self.AA_CODES.get(structure.residues[r][0].residue_name, "?")
            for r in range(start, end + 1)
            if r in structure.residues
        )

        logger.info(
            f"QM fragment: residues {start}-{end} | "
            f"sequence: {seq} | "
            f"atoms: {qm_fragment.n_atoms} | "
            f"MM atoms: {len(mm_atoms)}"
        )

        return qm_fragment, mm_atoms, mm_coords, mm_charges

    def atoms_to_pyscf_list(
        self,
        atoms: List[Atom],
        coords_angstrom: bool = True,
    ) -> List[Tuple[str, Tuple[float, float, float]]]:
        """
        Convert atoms to PySCF mol.atom format.

        Returns:
            [("C", (x, y, z)), ("N", (x, y, z)), ...]
            Coordinates in Angstrom (PySCF default) or Bohr.
        """
        ANGSTROM_TO_BOHR = 1.8897259886
        factor = 1.0 if coords_angstrom else ANGSTROM_TO_BOHR

        result = []
        for a in atoms:
            if not a.is_heavy:
                continue  # skip hydrogens for minimal basis
            elem = a.element if a.element else "C"
            result.append((elem, (a.x * factor, a.y * factor, a.z * factor)))
        return result


# ─────────────────────────────────────────────────────────────
# Main — test this stage
# ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import os
    os.chdir(Path(__file__).parent.parent)

    logger.info("=" * 60)
    logger.info("STAGE 1 — Data Acquisition Test")
    logger.info("=" * 60)

    dl = PDBDownloader(data_dir="data/raw")
    structure = dl.fetch("2I9M")     # Chignolin — best for your hardware
    logger.info(structure.summary())

    extractor = FragmentExtractor(max_qm_atoms=15)
    qm_frag, mm_atoms, mm_coords, mm_charges = extractor.extract(
        structure, center_residue=5, n_neighbors=2
    )

    logger.info(f"QM region atoms:   {qm_frag.n_atoms}")
    logger.info(f"MM environment:    {len(mm_atoms)} atoms")
    logger.info(f"MM coords shape:   {mm_coords.shape}")
    logger.success("Stage 1 PASSED")
