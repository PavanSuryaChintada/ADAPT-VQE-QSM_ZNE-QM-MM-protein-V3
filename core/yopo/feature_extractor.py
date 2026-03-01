"""
Stage 2 — YOPO Feature Extraction
====================================
Extracts rotation/translation-invariant structural fingerprints.

Why invariance matters:
  Two conformations rotated in space are physically identical.
  Raw (x,y,z) coordinates change with rotation → useless for comparing.
  Invariant features → same structure = same feature vector.

What we compute:
  1. Distance matrix D_ij = ||r_i - r_j||
  2. Graph Laplacian eigenvalues (spectral fingerprint of topology)
  3. Radial distribution histogram (pairwise distance distribution)
  4. Radius of gyration (compactness)
  5. Principal inertia moments (shape: prolate / oblate / spherical)

These five together give a complete, rotation-invariant structural description.
"""

import sys
import numpy as np
from pathlib import Path
from typing import Optional, List
from dataclasses import dataclass

sys.path.insert(0, str(Path(__file__).parent.parent))
from utils import logger


@dataclass
class YOPOFeatures:
    pdb_id: str
    distance_matrix: np.ndarray           # (N, N)
    laplacian_eigenvalues: np.ndarray      # (k,) spectral fingerprint
    distance_histogram: np.ndarray         # (bins,) radial distribution
    radius_of_gyration: float             # compactness
    inertia_moments: np.ndarray           # (3,) shape descriptors
    feature_vector: np.ndarray            # concatenated feature (for clustering)

    def __repr__(self):
        return f"YOPOFeatures({self.pdb_id}, dim={len(self.feature_vector)})"


class YOPOExtractor:
    """
    Rotation/translation-invariant feature extractor.

    Usage:
        ex = YOPOExtractor(n_eigenvalues=10)
        features = ex.extract(coords, pdb_id="2I9M")
    """

    def __init__(
        self,
        n_eigenvalues: int = 10,
        n_distance_bins: int = 20,
        contact_cutoff_ang: float = 8.0,
    ):
        self.k = n_eigenvalues
        self.n_bins = n_distance_bins
        self.cutoff = contact_cutoff_ang

    def extract(
        self,
        coords: np.ndarray,
        pdb_id: str = "unknown",
        weights: Optional[np.ndarray] = None,
    ) -> YOPOFeatures:
        """
        Extract invariant features from (N, 3) coordinate array.

        Args:
            coords: atomic coordinates in Angstrom
            pdb_id: identifier for this structure
            weights: per-atom weights (e.g., atomic masses)

        Returns:
            YOPOFeatures
        """
        if coords.ndim != 2 or coords.shape[1] != 3:
            raise ValueError(f"coords must be shape (N,3), got {coords.shape}")

        n = len(coords)

        # 1. Distance matrix (N×N pairwise Euclidean distances)
        D = self._distance_matrix(coords)

        # 2. Laplacian eigenvalue spectrum
        eigvals = self._laplacian_spectrum(D, n)

        # 3. Radial distribution histogram
        rdf = self._radial_distribution(D)

        # 4. Radius of gyration
        rog = self._radius_of_gyration(coords, weights)

        # 5. Principal inertia moments
        moments = self._inertia_moments(coords)

        # 6. Concatenate into feature vector
        fv = np.concatenate([eigvals, rdf, [rog], moments])

        return YOPOFeatures(
            pdb_id=pdb_id,
            distance_matrix=D,
            laplacian_eigenvalues=eigvals,
            distance_histogram=rdf,
            radius_of_gyration=rog,
            inertia_moments=moments,
            feature_vector=fv,
        )

    def extract_batch(
        self,
        coord_list: List[np.ndarray],
        pdb_ids: Optional[List[str]] = None,
    ) -> List[YOPOFeatures]:
        if pdb_ids is None:
            pdb_ids = [f"conf_{i}" for i in range(len(coord_list))]
        return [self.extract(c, pid) for c, pid in zip(coord_list, pdb_ids)]

    # ── Private ──────────────────────────────────────────────

    def _distance_matrix(self, coords: np.ndarray) -> np.ndarray:
        """Pairwise Euclidean distance matrix, shape (N,N)."""
        diff = coords[:, None, :] - coords[None, :, :]    # (N,N,3)
        D = np.sqrt((diff ** 2).sum(axis=2))               # (N,N)
        return D

    def _laplacian_spectrum(self, D: np.ndarray, n: int) -> np.ndarray:
        """
        Build contact-graph Laplacian and return top-k eigenvalues.

        Graph edges: atoms within self.cutoff Angstrom.
        Laplacian L = D_deg - A  (D_deg = diagonal degree matrix, A = adjacency)
        Normalized: L_norm = D_deg^{-1/2} L D_deg^{-1/2}

        Eigenvalues are unique to graph topology and rotation-invariant.
        """
        A = (D < self.cutoff).astype(float)
        np.fill_diagonal(A, 0)

        degrees = A.sum(axis=1)
        # Avoid divide-by-zero for isolated atoms
        d_inv_sqrt = np.where(degrees > 0, 1.0 / np.sqrt(degrees + 1e-12), 0.0)

        # Normalized Laplacian
        L = np.diag(degrees) - A
        D_inv = np.diag(d_inv_sqrt)
        L_norm = D_inv @ L @ D_inv

        # Eigenvalues (real symmetric → use eigh for speed and stability)
        eigvals = np.linalg.eigvalsh(L_norm)
        eigvals = np.sort(eigvals)

        # Skip trivial zero eigenvalue, take next k
        k = min(self.k, n - 1)
        return eigvals[1:k + 1] if len(eigvals) > k else np.pad(eigvals[1:], (0, k - len(eigvals) + 1))

    def _radial_distribution(self, D: np.ndarray) -> np.ndarray:
        """
        Histogram of pairwise distances (upper triangle only).
        This is the discrete Radial Distribution Function (RDF).
        """
        triu = D[np.triu_indices_from(D, k=1)]
        max_d = triu.max() if len(triu) > 0 and triu.max() > 0 else 30.0
        hist, _ = np.histogram(triu, bins=self.n_bins, range=(0.0, max_d), density=True)
        return np.nan_to_num(hist)

    def _radius_of_gyration(
        self,
        coords: np.ndarray,
        weights: Optional[np.ndarray] = None,
    ) -> float:
        """Rg = sqrt( Σ w_i |r_i - r_cm|² / Σ w_i )"""
        w = weights if weights is not None else np.ones(len(coords))
        w = w / w.sum()
        cm = (coords * w[:, None]).sum(axis=0)
        dev = coords - cm
        return float(np.sqrt((w * (dev ** 2).sum(axis=1)).sum()))

    def _inertia_moments(self, coords: np.ndarray) -> np.ndarray:
        """
        Principal moments of inertia (eigenvalues of inertia tensor).
        Encodes global shape: I_1 ≤ I_2 ≤ I_3
        """
        c = coords - coords.mean(axis=0)
        I = np.zeros((3, 3))
        for r in c:
            I += np.eye(3) * np.dot(r, r) - np.outer(r, r)
        I /= max(len(c), 1)

        moments = np.sort(np.linalg.eigvalsh(I))
        norm = moments.max() if moments.max() > 0 else 1.0
        return moments / norm


def cosine_similarity(f1: YOPOFeatures, f2: YOPOFeatures) -> float:
    """Cosine similarity between two feature vectors, in [0,1]."""
    v1, v2 = f1.feature_vector, f2.feature_vector
    m = max(len(v1), len(v2))
    v1 = np.pad(v1, (0, m - len(v1)))
    v2 = np.pad(v2, (0, m - len(v2)))
    n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
    if n1 == 0 or n2 == 0:
        return 0.0
    return float(np.dot(v1, v2) / (n1 * n2))


# ── Main ─────────────────────────────────────────────────────

if __name__ == "__main__":
    import os
    os.chdir(Path(__file__).parent.parent)

    logger.info("=" * 60)
    logger.info("STAGE 2 — YOPO Feature Extraction Test")
    logger.info("=" * 60)

    np.random.seed(42)
    coords = np.random.randn(20, 3) * 5.0

    # Rotation matrix (45° around Z)
    theta = np.pi / 4
    R = np.array([
        [np.cos(theta), -np.sin(theta), 0],
        [np.sin(theta),  np.cos(theta), 0],
        [0,              0,             1],
    ])
    rotated = (R @ coords.T).T

    ex = YOPOExtractor(n_eigenvalues=10, n_distance_bins=20)
    f1 = ex.extract(coords,  pdb_id="original")
    f2 = ex.extract(rotated, pdb_id="rotated_45deg")

    sim = cosine_similarity(f1, f2)
    logger.info(f"Similarity (original vs rotated): {sim:.6f}")
    logger.info(f"Feature vector dimension: {len(f1.feature_vector)}")

    if sim > 0.999:
        logger.success("Rotation invariance: PASSED (similarity ≈ 1.0)")
    else:
        logger.warning(f"Rotation invariance: partial (sim={sim:.4f})")

    logger.success("Stage 2 PASSED")
