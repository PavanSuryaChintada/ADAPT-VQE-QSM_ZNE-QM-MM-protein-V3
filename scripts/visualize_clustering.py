"""
scripts/visualize_clustering.py — DISCA Clustering Visualization
==================================================================
Runs stages 1–3 (Data → YOPO → DISCA) and plots the clustering result
for easier understanding.

Usage:
  python scripts/visualize_clustering.py
  python scripts/visualize_clustering.py --pdb 2JOF --fragment 5 --save
"""

import sys
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))


def run_pipeline_and_visualize(
    pdb_id: str = "2I9M",
    center_residue: int = 5,
    n_ensemble: int = 20,
    n_components: int = 2,
    save: bool = True,
    show: bool = True,
):
    """Run stages 1–3 and plot DISCA clustering."""
    from core.data_acquisition import PDBDownloader, FragmentExtractor
    from core.yopo.feature_extractor import YOPOExtractor
    from core.disca.clustering import DISCAClustering

    from utils.tee_stdout import start_tee_session
    log_path = start_tee_session("visualize_clustering")

    print(f"\n{'='*60}")
    print("  DISCA Clustering Visualization")
    print("="*60)
    print(f"  PDB: {pdb_id}  |  Center residue: {center_residue}")
    print("="*60)

    # Stage 1: Data
    dl = PDBDownloader(data_dir="data/raw")
    structure = dl.fetch(pdb_id)
    extractor = FragmentExtractor(max_qm_atoms=15)
    qm_frag, _, _, _ = extractor.extract(
        structure, center_residue=center_residue, n_neighbors=2
    )

    # Stage 2: YOPO features
    coords = qm_frag.get_coords()
    ex = YOPOExtractor(n_eigenvalues=10)
    features = ex.extract(coords, pdb_id=pdb_id)
    fv = features.feature_vector

    # Create ensemble (same as main.py)
    rng = np.random.default_rng(42)
    ensemble = np.vstack([fv + rng.normal(0, 0.05, len(fv)) for _ in range(n_ensemble)])

    # Standardize
    X_mean = ensemble.mean(axis=0)
    X_std = ensemble.std(axis=0) + 1e-8
    X_stdized = (ensemble - X_mean) / X_std

    # PCA to 2D for visualization
    d = min(2, X_stdized.shape[1], X_stdized.shape[0] - 1)
    _, _, Vt = np.linalg.svd(X_stdized, full_matrices=False)
    X_2d = X_stdized @ Vt[:d].T

    # Stage 3: DISCA
    clusterer = DISCAClustering(n_components=n_components, purity_threshold=0.70)
    result = clusterer.fit(ensemble)

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))

    labels = result.labels
    purities = result.purities
    compactness = result.compactness
    sizes = result.sizes
    best_k = result.best_cluster_id
    selection_score = purities * compactness
    # Chosen cluster = yellow; others = muted gray/blue
    YELLOW = "#FFD700"
    OTHER_COLORS = plt.cm.Greys(np.linspace(0.4, 0.7, result.n_clusters))

    for k in range(result.n_clusters):
        mask = labels == k
        pts = X_2d[mask]
        purity = purities[k]
        comp = compactness[k]
        score = selection_score[k]
        n_members = int(sizes[k])
        is_chosen = k == best_k
        color = YELLOW if is_chosen else OTHER_COLORS[k]
        ax.scatter(
            pts[:, 0],
            pts[:, 1],
            c=[color],
            marker="*" if is_chosen else "o",
            label=(
                f"Cluster {k}: purity={purity:.2f}, compactness={comp:.2f}, "
                f"score={score:.2f}, n={n_members}"
                + (" ← CHOSEN" if is_chosen else "")
            ),
            s=140 if is_chosen else 60,
            alpha=0.9 if is_chosen else 0.6,
            edgecolors="black",
            linewidths=1.5 if is_chosen else 0.5,
            zorder=10 if is_chosen else 1,
        )
        # Cluster center
        center = X_2d[labels == k].mean(axis=0)
        ax.scatter(
            center[0],
            center[1],
            c=[color],
            marker="X",
            s=250 if is_chosen else 150,
            edgecolors="black",
            linewidths=2.5 if is_chosen else 2,
            zorder=11 if is_chosen else 5,
        )

    # Box explaining why this cluster was chosen
    best_score = result.best_purity * result.best_compactness
    why_text = (
        f"CHOSEN: Cluster {best_k} (yellow)\n"
        f"Reason: Highest selection score = purity × compactness\n"
        f"        = {result.best_purity:.3f} × {result.best_compactness:.3f} = {best_score:.3f}\n"
        f"        (higher score = more representative conformation)"
    )
    props = dict(boxstyle="round,pad=0.5", facecolor="yellow", alpha=0.3, edgecolor="gold", linewidth=2)
    ax.text(0.02, 0.98, why_text, transform=ax.transAxes, fontsize=9, verticalalignment="top",
            bbox=props, family="monospace")

    ax.set_xlabel("PC1", fontsize=12)
    ax.set_ylabel("PC2", fontsize=12)
    ax.set_title(
        f"DISCA Clustering — {pdb_id} (res {center_residue})\n"
        f"Chosen cluster: {best_k} (yellow) | BIC: {result.bic:.1f}",
        fontsize=11,
    )
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_aspect("equal", adjustable="datalim")
    plt.tight_layout(rect=[0, 0.04, 1, 1])
    fig.text(0.5, 0.01, "Selection: best = argmax(purity × compactness). "
             "Compactness = 1/(1+variance). Chosen cluster highlighted in yellow.", ha="center", fontsize=8, style="italic")

    if save:
        out_dir = ROOT / "outputs" / "plots"
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / f"clustering_{pdb_id}_res{center_residue}.png"
        plt.savefig(path, dpi=150)
        print(f"\n  Plot saved: {path}")
    print(f"  Full session log: {log_path}")

    if show:
        plt.show()
    else:
        plt.close()

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Visualize DISCA clustering for a protein fragment"
    )
    parser.add_argument("--pdb", default="2I9M", help="PDB ID")
    parser.add_argument("--fragment", type=int, default=5, help="Center residue")
    parser.add_argument("--ensemble", type=int, default=20, help="Number of conformations")
    parser.add_argument("--clusters", type=int, default=2, help="Number of clusters")
    parser.add_argument("--save", action="store_true", help="Save plot to outputs/plots/")
    parser.add_argument("--no-show", action="store_true", help="Don't display plot")
    args = parser.parse_args()

    os.chdir(ROOT)

    run_pipeline_and_visualize(
        pdb_id=args.pdb,
        center_residue=args.fragment,
        n_ensemble=args.ensemble,
        n_components=args.clusters,
        save=args.save or True,
        show=not args.no_show,
    )


if __name__ == "__main__":
    main()
