"""Validate mdcif files against PDBx + IHM + the mdcif extension dictionary.

Loads mmcif_pdbx_v50.dic (from mmcif.wwpdb.org) and mmcif_ihm.dic (from
IHMCIF on GitHub), downloaded once per python-ihm version and cached, plus
docs/spec/mdcif_ext.dic, with ``ihm.dictionary.read``; combines them and runs
the combined dictionary's validator on the file. Also checks mdcif-specific
rules the dictionary can't express (unit-length vectors, one vector row per
selected atom, eigenvalues within the total variance, a standard mass for
every mode atom of a quasi-harmonic set, orthonormal global modes with an
origin, projections pointing at existing models).
"""

import os
import urllib.request
from pathlib import Path

import numpy as np

# IHM comes from IHMCIF's dist/, the dictionary python-ihm writes against: the copy on
# mmcif.wwpdb.org can lag a release behind (python-ihm 2.12 writes
# _pdbx_poly_seq_scheme.ihm_model_id_list, new in IHM 1.29, while wwpdb still had 1.28).
DICTIONARY_URLS = {
    "mmcif_pdbx_v50.dic": "https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic",
    "mmcif_ihm.dic": "https://raw.githubusercontent.com/ihmwg/IHMCIF/master/dist/mmcif_ihm.dic",
}
EXTENSION_DICTIONARY = Path(__file__).resolve().parents[2] / "docs" / "spec" / "mdcif_ext.dic"


class MdcifValidationError(Exception):
    """All problems found in a file, one per line."""

    def __init__(self, problems):
        self.problems = list(problems)
        super().__init__("\n".join(self.problems))


def cache_dir():
    base = os.environ.get("MDCIF_CACHE") or Path.home() / ".cache" / "mdcif"
    path = Path(base) / "dictionaries"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _dictionary(name):
    import ihm

    # one copy per python-ihm version: a newer python-ihm may write items that only a newer
    # dictionary defines, so a cached older one must not be reused after an upgrade
    path = cache_dir() / f"python-ihm-{ihm.__version__}" / name
    if not path.exists():
        path.parent.mkdir(exist_ok=True)
        with urllib.request.urlopen(DICTIONARY_URLS[name], timeout=300) as r:
            path.write_bytes(r.read())
    return path


def load_dictionary(extension_dictionary=None):
    """PDBx + IHM + mdcif extension, combined into one ``ihm.dictionary.Dictionary``."""
    import ihm.dictionary

    combined = None
    paths = [_dictionary(n) for n in DICTIONARY_URLS]
    paths.append(Path(extension_dictionary or EXTENSION_DICTIONARY))
    for p in paths:
        with open(p, encoding="utf-8", errors="replace") as fh:
            d = ihm.dictionary.read(fh)
        combined = d if combined is None else combined + d
    return combined


def _check_mode_set(label, ms, tol):
    problems = []
    norms = np.linalg.norm(ms.eigenvectors.reshape(ms.n_modes, -1), axis=1)
    for i, nrm in enumerate(norms):
        if abs(nrm - 1) > tol:
            problems.append(f"{label} mode {i + 1}: |v| = {nrm:.4f}, expected 1")
    vf = ms.variance_fraction  # eigenvalue / total_variance
    if np.any((vf < 0) | (vf > 1)) or vf.sum() > 1 + tol:
        problems.append(f"{label}: eigenvalues negative or above the total variance")
    if len(np.unique(ms.atom_indices)) != len(ms.atom_indices):
        problems.append(f"{label}: duplicate atoms in mode vectors")
    if ms.mass_weighted and not np.all(np.isfinite(ms.cartesian_variance)):
        problems.append(f"{label}: quasi-harmonic modes on atoms without a standard mass "
                        "(type_symbol), so their size in angstroms is unknown")
    return problems


def check_modes(path, tol=1e-3):
    """mdcif-specific rules; returns a list of problems (empty if none)."""
    from .reader import read_mdcif

    problems = []
    try:
        data = read_mdcif(path)
    except ValueError as exc:
        return [f"cannot map modes onto models: {exc}"]
    model_ids = {int(m._id) for _, m in data.system._all_models()}
    for model_id, ms in data.mode_sets.items():
        if model_id not in model_ids:
            problems.append(f"modes refer to unknown model {model_id}")
        problems.extend(_check_mode_set(f"model {model_id}", ms, tol))
    g = data.global_modes
    if g is not None:
        problems.extend(_check_mode_set("global", g, tol))
        if g.origin is None:
            problems.append("global mode set has no _mdcif_mode_origin rows")
        vecs = g.eigenvectors.reshape(g.n_modes, -1)
        if not np.allclose(vecs @ vecs.T, np.eye(g.n_modes), atol=10 * tol):
            problems.append("global modes are not orthonormal")
    p = data.projection
    if p is not None:
        if g is None:
            problems.append("_mdcif_projection without a global mode set")
        unknown = set(np.unique(p["model_id"])) - model_ids - {0}
        if unknown:
            problems.append(f"projections refer to unknown models {sorted(unknown)}")
    return problems


def validate_mdcif(path, extension_dictionary=None):
    """Raise an error describing every problem found in ``path``."""
    import ihm.dictionary

    path = Path(path)
    problems = []
    dictionary = load_dictionary(extension_dictionary)
    binary = path.suffix.lower() == ".bcif"
    with open(path, "rb" if binary else "r", **({} if binary else {"encoding": "utf-8"})) as fh:
        try:
            dictionary.validate(fh, format="BCIF" if binary else "mmCIF")
        except ihm.dictionary.ValidatorError as exc:
            problems.extend(str(exc).strip().splitlines())
    problems.extend(check_modes(path))
    if problems:
        raise MdcifValidationError(problems)
