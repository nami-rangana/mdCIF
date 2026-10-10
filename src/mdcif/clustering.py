"""Clusters from per-frame labels, with one representative (medoid) per cluster.

mdcif is not tied to one algorithm: any clustering that gives one label per
analysed frame (-1 = noise) can be written, from a scikit-learn-style estimator
(``fit_predict``) or from labels computed with any other tool. Built in is
HDBSCAN on the projections onto the top global PCA modes (the default of
``mdcif build``); since those modes are stored in the file, anyone can redraw
the clusters from it.

Clusters are sorted by population (cluster 0 = most populated). The
representative is the medoid (frame with the lowest mean RMSD to the rest of
its cluster), which maps onto ``ihm_model_representative``.
"""

from dataclasses import dataclass, field

import numpy as np

from .trajectory import rmsd_to

DENSITY = "Density based threshold-clustering"
# ihm_ensemble_info.ensemble_clustering_method for well-known algorithms (by class name);
# anything else is "Other"
IHM_METHODS = {
    "HDBSCAN": DENSITY, "DBSCAN": DENSITY, "OPTICS": DENSITY,
    "KMeans": "Partitioning (k-means)", "MiniBatchKMeans": "Partitioning (k-means)",
    "BisectingKMeans": "Partitioning (k-means)",
    "AgglomerativeClustering": "Hierarchical", "Birch": "Hierarchical",
}


@dataclass
class ClusterResult:
    labels: np.ndarray            # (n_frames,) cluster index per frame, -1 = noise
    representatives: np.ndarray   # (n_clusters,) frame index of each medoid
    populations: np.ndarray       # (n_clusters,) number of frames per cluster
    precision: np.ndarray         # (n_clusters,) mean RMSD of members to the medoid (angstrom)
    method: str                   # algorithm name, e.g. "HDBSCAN"
    feature: str                  # what was clustered, e.g. "PCA projection (PCs 1-3)"
    ihm_method: str = "Other"     # ihm_ensemble_info.ensemble_clustering_method
    ihm_feature: str = "other"    # ihm_ensemble_info.ensemble_clustering_feature
    software: tuple | None = None  # (name, version, url) of the clustering software
    details: str = ""
    params: dict = field(default_factory=dict)  # algorithm parameters actually used

    @property
    def n_clusters(self):
        return len(self.representatives)

    @property
    def noise_fraction(self):
        return float(np.mean(self.labels < 0))

    def members(self, cluster):
        return np.flatnonzero(self.labels == cluster)


def default_min_cluster_size(n_frames):
    """HDBSCAN default used by mdcif: 1% of the frames, at least 25."""
    return max(25, int(round(0.01 * n_frames)))


def hdbscan_labels(features, min_cluster_size=None, min_samples=None, cluster_selection="eom"):
    """Built-in clustering: scikit-learn HDBSCAN on ``features`` (PCs 1-3).

    Returns ``(labels, params)``. ``cluster_selection`` "eom" gives few, large
    clusters; "leaf" splits them into substates and leaves more frames as noise.
    """
    from sklearn.cluster import HDBSCAN

    features = np.asarray(features, dtype=np.float64)
    min_cluster_size = min_cluster_size or default_min_cluster_size(len(features))
    labels = HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples,
                     cluster_selection_method=cluster_selection, copy=True).fit_predict(features)
    if np.all(labels < 0):
        raise ValueError("HDBSCAN put every frame in noise; lower min_cluster_size "
                         f"(was {min_cluster_size}) or min_samples")
    params = {"min_cluster_size": min_cluster_size,
              "min_samples": min_samples if min_samples is not None else min_cluster_size,
              "cluster_selection": cluster_selection}
    return labels, params


def medoid(X, n_atoms, exact_max=2000, n_candidates=256):
    """Index of the row of ``X`` (frames x 3N) with the lowest summed RMSD to the others.

    Exact all-pairs for up to ``exact_max`` frames. Above that, the medoid is
    searched among the ``n_candidates`` frames closest to the geometric median,
    each scored exactly against all frames.
    """
    X = np.asarray(X, dtype=np.float64)
    n = len(X)
    sq = (X ** 2).sum(axis=1)
    if n <= exact_max:
        cand = np.arange(n)
    else:
        g = X.mean(axis=0)
        for _ in range(200):
            w = 1.0 / np.maximum(np.linalg.norm(X - g, axis=1), 1e-12)
            g_new = w @ X / w.sum()
            done = np.linalg.norm(g_new - g) <= 1e-7 * (1.0 + np.linalg.norm(g))
            g = g_new
            if done:
                break
        cand = np.argpartition(np.linalg.norm(X - g, axis=1), n_candidates - 1)[:n_candidates]
    total = np.empty(len(cand))
    for start in range(0, len(cand), 1000):
        c = cand[start:start + 1000]
        d2 = sq[c, None] + sq[None, :] - 2.0 * X[c] @ X.T
        total[start:start + 1000] = np.sqrt(np.maximum(d2, 0.0) / n_atoms).sum(axis=1)
    return int(cand[np.argmin(total)])


def clusters_from_labels(coords, labels, *, method, feature, ihm_method=None,
                         ihm_feature="other", software=None, details="",
                         params=None) -> ClusterResult:
    """Clusters, medoids and precision from per-frame ``labels`` (-1 = noise).

    ``coords`` are the aligned analysis-atom coordinates ``(n_frames, n_atoms, 3)``
    used for medoids and precision whatever was clustered. Clusters are
    renumbered by decreasing population.
    """
    coords = np.asarray(coords, dtype=np.float64)
    labels = np.asarray(labels).astype(int)
    n_frames, n_atoms, _ = coords.shape
    if labels.shape != (n_frames,):
        raise ValueError(f"need one label per analysed frame ({n_frames}), got {labels.shape}")
    found = [c for c in np.unique(labels) if c >= 0]
    if not found:
        raise ValueError("no clusters: every frame is labelled as noise (-1)")
    order = sorted(found, key=lambda c: -np.sum(labels == c))
    remap = {old: new for new, old in enumerate(order)}
    labels = np.array([remap.get(c, -1) for c in labels])

    X = coords.reshape(n_frames, -1)
    reps, pops, prec = [], [], []
    for c in range(len(order)):
        idx = np.flatnonzero(labels == c)
        rep = idx[medoid(X[idx], n_atoms)]
        reps.append(rep)
        pops.append(len(idx))
        prec.append(float(rmsd_to(coords[idx], coords[rep]).mean()))
    return ClusterResult(labels=labels, representatives=np.array(reps),
                         populations=np.array(pops), precision=np.array(prec),
                         method=method, feature=feature,
                         ihm_method=ihm_method or IHM_METHODS.get(method, "Other"),
                         ihm_feature=ihm_feature, software=software, details=details,
                         params=dict(params or {}))
