# mdcif — a shareable mmCIF object for MD trajectories

> Status: **pre-alpha**. The core pipeline (`mdcif build | info | validate | landscape`) works
> end to end (format draft v0.2); the PyMOL plugin works and is tested; the ChimeraX bundle is
> implemented but has not yet been run inside ChimeraX.

The PDB works because everyone shares one object (PDB / mmCIF) without needing to know how it
was produced. MD trajectories have no such object: they are large, engine-specific and hard to
use outside the MD community. Cluster centroids are easy to share but lose the dynamics.

**mdcif** packages an MD simulation into a single, PDB-IHM-compatible mmCIF file containing:

1. **Models** — the medoids of the main clusters (multi-model `atom_site`). Frames are
   clustered with HDBSCAN on their projections onto the top 3 global PCs.
2. **Provenance** — simulation engine, force field and conditions, plus how clustering and PCA
   were done (existing PDBx/IHM categories wherever possible, a small `_mdcif_` extension otherwise).
   Given a PDB entry of the same molecule (`--reference 3GB1`), the file reuses its entity,
   chains, numbering, source and citations, and the models are superposed on it.
3. **Essential modes** — the global top-3 PCs (the clustering axes, with their mean structure)
   and, per model, the top *k* (default 3) PCA / quasi-harmonic eigenvectors of its own cluster,
   with eigenvalues and % variance explained. This replaces the elastic-network approximation
   with modes taken from the actual simulated ensemble.
4. **The landscape** — every frame's (or a subsample's) coordinates on PCs 1–3 and its cluster,
   so anyone can redraw the clusters from the file alone (`mdcif landscape`).

…and viewer plugins for **ChimeraX** and **PyMOL** to morph between models, animate along each
model's own modes (or the global PCs), and draw porcupine plots.

Worked example: GB1 folded with MELD × Amber, written as an edited copy of PDB 3GB1 —
[examples/3gb1/](examples/3gb1/README.md).

See [docs/proposal.md](docs/proposal.md) for goals and milestones and
[docs/spec/mdcif_spec.md](docs/spec/mdcif_spec.md) for the draft format.

## Repository layout

| Path | Contents |
|------|----------|
| `src/mdcif/` | Core Python package: trajectory loading, clustering, modes, provenance, mmCIF writer/reader/validator, CLI |
| `viewers/chimerax/` | ChimeraX bundle `ChimeraX-MDCIF` (`mdcif ...` commands) |
| `viewers/pymol/mdcif_pymol/` | PyMOL plugin (`mdcif_*` commands) |
| `docs/spec/` | Format spec, IHM/PDBx category mapping, DDL2 extension dictionary |
| `examples/` | Provenance template, 3GB1 MELD example, notebooks; `examples/data/` holds local data (git-ignored) |
| `tests/` | pytest suite; small fixtures in `tests/data/` |
| `scripts/` | One-off utilities |

## Setup

### 1. Install the tools (once)

- **Miniforge** (conda + mamba, conda-forge by default): <https://github.com/conda-forge/miniforge>
- **Git for Windows**: <https://git-scm.com/download/win>
- **UCSF ChimeraX** ≥ 1.6 (needed for pyproject-based bundles): <https://www.cgl.ucsf.edu/chimerax/download.html>

### 2. Create the environment

From the repository root, in a Miniforge Prompt:

```bash
mamba env create -f environment.yml
```

```bash
conda activate mdcif
```

```bash
pip install --no-deps -r requirements-pip.txt
```

```bash
pip install -e . --no-deps
```

```bash
pytest
```

`requirements-pip.txt` adds MDANCE (optional k-means NANI clustering) and its PyTorch
dependency, which conda-forge can't provide here. The default HDBSCAN clustering only needs
scikit-learn. `mdcif` is installed in editable mode, so code changes take effect without
reinstalling. Always use `--no-deps`: conda provides all other dependencies, and letting pip
re-resolve them makes it try to compile ProDy from source on Windows. After editing
`environment.yml`, run `mamba env update -f environment.yml --prune`.

On Windows with **Smart App Control** on, conda-forge's scipy and scikit-learn DLLs are blocked;
replace them with the PyPI wheels as described at the top of `requirements-pip.txt`. The env sets
`MKL_THREADING_LAYER=TBB` because conda-forge MKL's default threading layer crashes on win-64.

On Linux (e.g. HiPerGator) the same steps work, with two differences: install torch from the
CPU index (`pip install --no-deps torch --index-url https://download.pytorch.org/whl/cpu`)
before the rest of `requirements-pip.txt`, and clear `PYTHONPATH` if a module (such as Amber's)
puts another Python's site-packages on it. Run trajectory analyses on a compute node
(`srun ...`), not a login node.

### Quick run

```bash
python scripts/fetch_mdposit.py MD-A000GI --out examples/data/mdposit
mdcif build examples/data/mdposit/MD-A000GI/structure.pdb examples/data/mdposit/MD-A000GI/replica_1.xtc --provenance examples/data/mdposit/MD-A000GI/provenance.yaml --reference 2GB1 --output gb1.cif
mdcif info gb1.cif
mdcif validate gb1.cif
mdcif landscape gb1.cif
```

Useful `mdcif build` options: `--start/--stop/--step` (frame range, e.g. skip equilibration),
`--min-cluster-size`, `--cluster-selection eom|leaf` (HDBSCAN), `--max-models`,
`--method nani|kmeans --n-clusters k` (k-means instead of HDBSCAN), `--reference <PDB ID or
mmCIF>`, `--max-projection-points`, `--format` (if the trajectory format can't be guessed;
Amber `name.nc.001` is recognised).

### 3. Initialize git

```bash
git init -b main
```

```bash
pre-commit install
```

```bash
git add .
```

```bash
git commit -m "Scaffold mdcif repository"
```

Then create an empty repository on GitHub and follow its instructions to add the remote and push.

### 4. Load the ChimeraX bundle (development mode)

In the ChimeraX command line (use forward slashes):

```
devel install C:/path/to/mdcif/viewers/chimerax
mdcif open gb1.cif
mdcif modes #1
mdcif morph #1
mdcif animate #1.1 mode 1 amplitude 3
mdcif animate #1.2 mode 1 scope global
mdcif porcupine #1.1 mode 1 scale 3 color orange
```

Re-run `devel install` after changing the bundle. `amplitude`/`scale` are in standard deviations
of the mode.

### 5. Load the PyMOL plugin (development mode)

With the `mdcif` environment active, start `pymol` and type this in its command line (the leading
`/` runs a line of Python):

```
/import sys; sys.path.insert(0, "C:/path/to/mdcif/viewers/pymol"); import mdcif_pymol; mdcif_pymol.__init_plugin__()
mdcif_load gb1.cif, gb1
mdcif_morph gb1
mdcif_animate gb1, model=1, mode=1, amplitude=3
mdcif_animate gb1, model=2, mode=1, scope=global
mdcif_porcupine gb1, model=1, mode=1, scale=3
```

Or install it permanently via *Plugin → Plugin Manager → Install New Plugin* and choose
`viewers/pymol/mdcif_pymol/__init__.py`.

## Data

Trajectories are never committed (see `.gitignore`). Tests and notebooks use the small public
adenylate kinase (AdK) trajectories shipped with `MDAnalysisTests`. Put any other local data in
`examples/data/`.

## License

MIT — see [LICENSE](LICENSE).
