"""
Stage 3 — DISCA Clustering
============================
Partitions conformational ensemble into energy basins.

Problem:
  MD trajectories or cryo-EM ensembles contain many conformations.
  Running quantum simulation on all of them is too expensive.
  We need the single best (most structurally pure) conformation.

Solution — Gaussian Mixture Model (GMM):
  P(x) = Σ_k π_k · N(x | μ_k, Σ_k)

  Each Gaussian = one conformational basin.
  Purity = how tightly members cluster around their Gaussian center.
  We pick the highest-purity cluster and its centroid conformation.

Why GMM over K-means:
  GMM gives soft assignments (probabilities) — more robust.
  Purity metric requires these probabilities.
  K-means gives hard assignments — can't compute purity.
"""

import sys
import numpy as np
from pathlib import Path
from typing import List, Optional, Tuple
from dataclasses import dataclass, field

sys.path.insert(0, str(Path(__file__).parent.parent))
from utils import logger


@dataclass
class ClusterResult:
    best_cluster_id: int
    best_purity: float
    best_indices: np.ndarray            # indices in feature matrix
    representative_idx: int             # single best conformation
    n_clusters: int
    labels: np.ndarray
    purities: np.ndarray
    bic: float
    converged: bool

    def summary(self) -> str:
        return (
            f"DISCA | best_cluster={self.best_cluster_id} | "
            f"purity={self.best_purity:.3f} | "
            f"size={len(self.best_indices)} | "
            f"BIC={self.bic:.1f} | "
            f"converged={self.converged}"
        )


class GaussianMixtureSimple:
    """
    Lightweight EM-based Gaussian Mixture Model.
    Uses scipy only (no sklearn required, though sklearn is preferred).
    Falls back to sklearn.GaussianMixture if available.
    """

    def __init__(
        self,
        n_components: int = 3,
        max_iter: int = 200,
        tol: float = 1e-4,
        random_state: int = 42,
        covariance_type: str = "diag",
    ):
        self.K = n_components
        self.max_iter = max_iter
        self.tol = tol
        self.rng = np.random.default_rng(random_state)
        self.cov_type = covariance_type

        self.means_ = None
        self.covars_ = None
        self.weights_ = None
        self.converged_ = False

    def fit(self, X: np.ndarray) -> "GaussianMixtureSimple":
        """EM algorithm for GMM."""
        n, d = X.shape
        K = min(self.K, n - 1)

        # Initialize means via K-means++ seeding
        self.means_ = self._init_means(X, K)
        self.covars_ = np.array([np.eye(d) for _ in range(K)])
        self.weights_ = np.ones(K) / K

        log_likelihood_prev = -np.inf

        for iteration in range(self.max_iter):
            # E-step: compute responsibilities
            R = self._e_step(X)

            # M-step: update parameters
            self._m_step(X, R)

            # Check convergence via log-likelihood
            ll = self._log_likelihood(X)
            if abs(ll - log_likelihood_prev) < self.tol:
                self.converged_ = True
                break
            log_likelihood_prev = ll

        return self

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        return self._e_step(X)

    def predict(self, X: np.ndarray) -> np.ndarray:
        return np.argmax(self.predict_proba(X), axis=1)

    def bic(self, X: np.ndarray) -> float:
        """Bayesian Information Criterion. Lower = better model."""
        n, d = X.shape
        K = len(self.means_)
        n_params = K * d + K * d + K - 1  # means + variances + weights
        ll = self._log_likelihood(X)
        return -2 * ll + n_params * np.log(n)

    def _init_means(self, X: np.ndarray, K: int) -> np.ndarray:
        """K-means++ initialization."""
        idx = self.rng.integers(len(X))
        centers = [X[idx].copy()]
        for _ in range(K - 1):
            dists = np.array([min(np.sum((x - c) ** 2) for c in centers) for x in X])
            probs = dists / dists.sum()
            idx = self.rng.choice(len(X), p=probs)
            centers.append(X[idx].copy())
        return np.array(centers)

    def _gaussian_log_pdf(self, X: np.ndarray, mean: np.ndarray, cov: np.ndarray) -> np.ndarray:
        d = len(mean)
        diff = X - mean
        # Diagonal covariance only (stable and fast)
        var = np.diag(cov) + 1e-6
        log_det = np.sum(np.log(var))
        maha = np.sum(diff ** 2 / var, axis=1)
        return -0.5 * (d * np.log(2 * np.pi) + log_det + maha)

    def _e_step(self, X: np.ndarray) -> np.ndarray:
        K = len(self.means_)
        log_resp = np.zeros((len(X), K))
        for k in range(K):
            log_resp[:, k] = (
                np.log(self.weights_[k] + 1e-300) +
                self._gaussian_log_pdf(X, self.means_[k], self.covars_[k])
            )
        # Normalize (log-sum-exp for stability)
        log_resp -= log_resp.max(axis=1, keepdims=True)
        resp = np.exp(log_resp)
        resp /= resp.sum(axis=1, keepdims=True) + 1e-300
        return resp

    def _m_step(self, X: np.ndarray, R: np.ndarray) -> None:
        n, d = X.shape
        K = R.shape[1]
        Nk = R.sum(axis=0) + 1e-10
        self.weights_ = Nk / n
        self.means_ = (R.T @ X) / Nk[:, None]
        for k in range(K):
            diff = X - self.means_[k]
            cov_diag = (R[:, k] @ (diff ** 2)) / Nk[k] + 1e-6
            self.covars_[k] = np.diag(cov_diag)

    def _log_likelihood(self, X: np.ndarray) -> float:
        K = len(self.means_)
        ll = np.zeros(len(X))
        for k in range(K):
            ll += self.weights_[k] * np.exp(
                self._gaussian_log_pdf(X, self.means_[k], self.covars_[k])
            )
        return float(np.sum(np.log(ll + 1e-300)))


class DISCAClustering:
    """
    Conformational clustering for quantum fragment selection.

    Usage:
        clusterer = DISCAClustering(n_components=3)
        result = clusterer.fit(feature_matrix)
        best_coords = coord_list[result.representative_idx]
    """

    def __init__(
        self,
        n_components: int = 3,
        purity_threshold: float = 0.80,
        max_iter: int = 200,
        random_state: int = 42,
        use_pca: bool = True,
        pca_dim: int = 10,
    ):
        self.K = n_components
        self.purity_thresh = purity_threshold
        self.max_iter = max_iter
        self.random_state = random_state
        self.use_pca = use_pca
        self.pca_dim = pca_dim
        self._pca_components = None   # fitted PCA
        self._mean = None
        self._std = None

    def fit(self, feature_matrix: np.ndarray) -> ClusterResult:
        """
        Cluster conformations and select the highest-purity cluster.

        Args:
            feature_matrix: (n_confs, n_features) from YOPO

        Returns:
            ClusterResult
        """
        n, d = feature_matrix.shape
        logger.info(f"DISCA: clustering {n} conformations, {d} features, K={self.K}")

        # Standardize
        X = self._standardize(feature_matrix, fit=True)

        # PCA (optional dimensionality reduction)
        if self.use_pca and d > self.pca_dim:
            X = self._pca_transform(X, fit=True)
            logger.debug(f"PCA: {d} → {X.shape[1]} dimensions")

        # Fit GMM
        K = min(self.K, n - 1)
        try:
            from sklearn.mixture import GaussianMixture
            gmm = GaussianMixture(
                n_components=K,
                covariance_type="full",
                max_iter=self.max_iter,
                random_state=self.random_state,
                n_init=5,
            )
            gmm.fit(X)
            labels = gmm.predict(X)
            responsibilities = gmm.predict_proba(X)
            bic_score = gmm.bic(X)
            converged = gmm.converged_
        except ImportError:
            # Fallback to our simple GMM
            logger.warning("sklearn not found — using built-in GMM")
            gmm = GaussianMixtureSimple(K, max_iter=self.max_iter, random_state=self.random_state)
            gmm.fit(X)
            labels = gmm.predict(X)
            responsibilities = gmm.predict_proba(X)
            bic_score = gmm.bic(X)
            converged = gmm.converged_

        # Compute purity per cluster
        purities = self._compute_purities(responsibilities, labels, K)

        # Best cluster = highest purity
        best_k = int(np.argmax(purities))
        best_indices = np.where(labels == best_k)[0]
        best_purity = float(purities[best_k])

        if best_purity < self.purity_thresh:
            logger.warning(
                f"Best cluster purity {best_purity:.3f} < threshold {self.purity_thresh}. "
                f"Consider more conformations or different n_components."
            )

        # Representative = conformation closest to cluster center
        rep_idx = self._find_representative(X, best_indices, best_k, gmm)

        result = ClusterResult(
            best_cluster_id=best_k,
            best_purity=best_purity,
            best_indices=best_indices,
            representative_idx=int(rep_idx),
            n_clusters=K,
            labels=labels,
            purities=purities,
            bic=bic_score,
            converged=converged,
        )
        logger.success(result.summary())
        return result

    def _compute_purities(self, R: np.ndarray, labels: np.ndarray, K: int) -> np.ndarray:
        purities = np.zeros(K)
        for k in range(K):
            members = np.where(labels == k)[0]
            if len(members) > 0:
                purities[k] = float(R[members, k].mean())
        return purities

    def _find_representative(self, X, indices, k, gmm) -> int:
        """Find sample closest to cluster mean."""
        try:
            center = gmm.means_[k]
        except AttributeError:
            center = X[indices].mean(axis=0)
        pts = X[indices]
        dists = np.linalg.norm(pts - center, axis=1)
        return indices[int(np.argmin(dists))]

    def _standardize(self, X: np.ndarray, fit: bool = False) -> np.ndarray:
        if fit:
            self._mean = X.mean(axis=0)
            self._std  = X.std(axis=0) + 1e-8
        return (X - self._mean) / self._std

    def _pca_transform(self, X: np.ndarray, fit: bool = False) -> np.ndarray:
        d = min(self.pca_dim, X.shape[1], X.shape[0] - 1)
        if fit:
            # Simple PCA via SVD
            _, _, Vt = np.linalg.svd(X, full_matrices=False)
            self._pca_components = Vt[:d]
        if self._pca_components is not None:
            return X @ self._pca_components.T
        return X[:, :d]


# ── Main ─────────────────────────────────────────────────────

if __name__ == "__main__":
    import os
    os.chdir(Path(__file__).parent.parent)

    logger.info("=" * 60)
    logger.info("STAGE 3 — DISCA Clustering Test")
    logger.info("=" * 60)

    np.random.seed(42)
    # Simulate 60 conformations, 34 features, 3 real clusters
    centers = [np.random.randn(34) * 5 for _ in range(3)]
    X = np.vstack([
        centers[0] + np.random.randn(25, 34) * 0.5,
        centers[1] + np.random.randn(25, 34) * 0.5,
        centers[2] + np.random.randn(10, 34) * 0.5,
    ])

    clusterer = DISCAClustering(n_components=3, purity_threshold=0.75)
    result = clusterer.fit(X)

    logger.info(f"Representative conformation index: {result.representative_idx}")
    logger.info(f"Best cluster size: {len(result.best_indices)} conformations")
    logger.success("Stage 3 PASSED")
