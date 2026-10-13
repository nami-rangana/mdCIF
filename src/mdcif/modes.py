"""Essential modes (PCA / quasi-harmonic): global, per cluster, and the ENM baseline.

The global set is the PCA of all analysed frames (aligned on their average);
its top modes are the axes the frames are clustered on. For each cluster,
frames are aligned to the cluster representative, the covariance matrix of the
selected atoms is diagonalized, and the top ``k`` eigenvectors are kept.
Vectors are stored as unit-length Cartesian displacements in the frame of the
deposited models (see docs/spec/mdcif_spec.md). ``cartesian_variance`` is the
variance of the motion along each of them in square angstroms: the eigenvalue
for PCA; for mass-weighted (quasi-harmonic) modes, whose eigenvalues are in
amu * square angstroms, it is not stored but follows from the eigenvalue and
the standard masses of the atoms' elements.
"""

from dataclasses import dataclass

import numpy as np
from scipy.linalg import eigh

from .trajectory import superpose


@dataclass
class ModeSet:
    eigenvalues: np.ndarray        # (k,) angstrom^2 (amu*angstrom^2 if mass_weighted)
    eigenvectors: np.ndarray       # (k, n_atoms, 3) unit 3N vectors
    variance_fraction: np.ndarray  # (k,)
    method: str                    # "PCA" | "quasi-harmonic" | "ANM"
    mass_weighted: bool
    selection: str
    n_frames: int
    origin: np.ndarray | None = None      # (n_atoms, 3) mean structure the modes are about
    total_variance: float | None = None   # trace of the covariance matrix
    # (k,) angstrom^2: variance of the Cartesian motion along each unit vector, the display
    # amplitude (= eigenvalues unless mass_weighted)
    cartesian_variance: np.ndarray | None = None

    def __post_init__(self):
        if self.cartesian_variance is None and not self.mass_weighted:
            self.cartesian_variance = np.asarray(self.eigenvalues, dtype=np.float64)

    @property
    def n_modes(self):
        return len(self.eigenvalues)

    @property
    def cumulative_variance_fraction(self):
        return np.cumsum(self.variance_fraction)


def cluster_modes(coords, n_modes=3, masses=None, reference=None, selection="name CA") -> ModeSet:
    """PCA (or mass-weighted quasi-harmonic analysis if ``masses`` is given) of one cluster.

    ``coords`` is ``(n_frames, n_atoms, 3)``. If ``reference`` (``(n_atoms, 3)``,
    normally the cluster representative) is given, frames are first superposed
    on it so the modes are expressed in the representative's frame. The sign of
    each vector is fixed so that its largest component is positive.

    Mass-weighted eigenvectors ``u`` are returned as Cartesian unit vectors
    ``v ~ M^-1/2 u``; the variance of the Cartesian motion along ``v`` is then
    ``lambda / sum_j m_j |v_j|^2`` (``cartesian_variance``, angstrom^2).
    """
    coords = np.asarray(coords, dtype=np.float64)
    if reference is not None:
        coords = superpose(coords, reference)[0]
    n_frames, n_atoms, _ = coords.shape
    if n_frames < 2:
        raise ValueError(
            f"Need at least 2 frames to compute modes. "
            f"Given trajectory contains only {n_frames}"
        )
    X = coords.reshape(n_frames, -1)
    w = None
    if masses is not None:
        w = np.sqrt(np.repeat(np.asarray(masses, dtype=np.float64), 3))
        X = X * w
    mean = X.mean(axis=0)
    X = X - mean
    total = float(np.vdot(X, X)) / (n_frames - 1)  # trace of the covariance
    evals, vecs = _top_eigenpairs(X, n_modes)
    evals = np.maximum(evals / (n_frames - 1), 0.0)
    if w is not None:
        vecs = vecs / w
    vecs /= np.linalg.norm(vecs, axis=1, keepdims=True)
    cart = evals if w is None else evals / (vecs ** 2 @ w ** 2)
    # sign is arbitrary (spec); fixed so the output does not depend on solver or LAPACK build
    vecs *= np.sign(vecs[np.arange(len(vecs)), np.abs(vecs).argmax(axis=1)])[:, None]
    return ModeSet(
        eigenvalues=evals,
        eigenvectors=vecs.reshape(len(evals), n_atoms, 3),
        variance_fraction=evals / total,
        method="PCA" if masses is None else "quasi-harmonic",
        mass_weighted=masses is not None,
        selection=selection,
        n_frames=n_frames,
        origin=mean.reshape(n_atoms, 3) / (w.reshape(n_atoms, 3) if w is not None else 1.0),
        total_variance=total,
        cartesian_variance=cart,
    )


def _top_eigenpairs(X, k):
    """Largest ``k`` eigenpairs of ``X^T X`` for centred ``X`` ``(n_frames, 3N)``, largest first.

    Exact, and only the top ``k`` pairs are computed. With fewer frames than
    coordinates the small ``X X^T`` (frames x frames) is diagonalized instead:
    same non-zero eigenvalues, eigenvectors ``X^T u`` (returned unnormalized),
    and the 3N x 3N covariance is never formed.
    """
    n, d = X.shape
    k = min(k, d)
    if k < n < d:
        lam, u = eigh(X @ X.T, subset_by_index=[n - k, n - 1])
        if lam[0] > 1e-12 * lam[-1]:  # X^T u is undefined for a zero eigenvalue
            return lam[::-1], (X.T @ u[:, ::-1]).T
    lam, v = eigh(X.T @ X, subset_by_index=[d - k, d - 1])
    return lam[::-1], v[:, ::-1].T


def project(coords, ms: ModeSet) -> np.ndarray:
    """Coordinates of already superposed frames ``(n_frames, n_atoms, 3)`` along each mode.

    ``(x - origin) . v_i`` in angstroms, for modes of an unweighted PCA (whose
    vectors are orthonormal); returns ``(n_frames, n_modes)``.
    """
    if ms.mass_weighted:
        raise ValueError("projections are only defined for unweighted PCA modes")
    d = np.asarray(coords, dtype=np.float64) - ms.origin
    return d.reshape(len(d), -1) @ ms.eigenvectors.reshape(ms.n_modes, -1).T


def anm_hessian(coords, cutoff=15.0, gamma=1.0):
    """Hessian ``(3N, 3N)`` of the anisotropic network model (Atilgan et al. 2001).

    Off-diagonal blocks ``-gamma * d d^T / |d|^2`` for pairs closer than
    ``cutoff``; diagonal blocks make every block row sum to zero (as ProDy's ANM).
    """
    x = np.asarray(coords, dtype=np.float64)
    n = len(x)
    d = x[None, :, :] - x[:, None, :]
    r2 = (d ** 2).sum(axis=2)
    k = np.zeros((n, n))
    contact = (r2 < cutoff ** 2) & ~np.eye(n, dtype=bool)
    k[contact] = -gamma / r2[contact]
    h = k[:, :, None, None] * d[:, :, :, None] * d[:, :, None, :]
    h[np.arange(n), np.arange(n)] = -h.sum(axis=1)
    return h.transpose(0, 2, 1, 3).reshape(3 * n, 3 * n)


def anm_modes(coords, n_modes=3, cutoff=15.0, selection="name CA") -> ModeSet:
    """Anisotropic network model modes of a single structure, as a baseline.

    Eigenvalues are reported as ``1 / lambda`` (relative variances, arbitrary
    units) so they rank like PCA eigenvalues; ``variance_fraction`` is relative
    to the sum over all non-trivial ANM modes.

    The rigid-body modes (6; more if the network falls apart at ``cutoff``) are
    those below ``1e-10`` of the largest eigenvalue: they sit at round-off
    (~1e-16), while the softest bending modes of long extended chains reach
    ~1e-9 (a relative ``1e-6`` cut dropped some of these).
    """
    h = anm_hessian(coords, cutoff)
    evals = eigh(h, eigvals_only=True)  # all of them, for the zero count and the total
    n_zero = int(np.sum(evals < 1e-10 * evals[-1]))
    var = 1.0 / evals[n_zero:]
    k = min(n_modes, len(var))
    vecs = eigh(h, subset_by_index=[n_zero, n_zero + k - 1])[1].T
    return ModeSet(
        eigenvalues=var[:k],
        eigenvectors=vecs.reshape(k, -1, 3),
        variance_fraction=var[:k] / var.sum(),
        method="ANM",
        mass_weighted=False,
        selection=selection,
        n_frames=1,
    )


def overlap_matrix(a: ModeSet, b: ModeSet) -> np.ndarray:
    """``|v_i . w_j|`` between two mode sets, e.g. MD modes vs ANM modes."""
    va = a.eigenvectors.reshape(a.n_modes, -1)
    vb = b.eigenvectors.reshape(b.n_modes, -1)
    return np.abs(va @ vb.T)


def cumulative_overlap(a: ModeSet, b: ModeSet) -> np.ndarray:
    """Cumulative overlap of each mode of ``a`` with the space of ``b`` (Yang et al. 2008).

    ``sqrt(sum_j (v_i . w_j)^2)``: 1 if ``v_i`` is a combination of the ``b``
    modes, 0 if orthogonal to all of them; ``b`` must be orthonormal (PCA or
    ANM, not mass-weighted). Computed by ``MDAnalysis.analysis.pca``.
    """
    from MDAnalysis.analysis.pca import cumulative_overlap as overlap

    va = a.eigenvectors.reshape(a.n_modes, -1)
    vb = b.eigenvectors.reshape(b.n_modes, -1)
    return np.array([overlap(va, vb, i) for i in range(a.n_modes)])
