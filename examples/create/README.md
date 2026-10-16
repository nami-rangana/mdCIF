# Creating mdCIF files: tutorials

Four notebooks, all on the same simulation: the protein G B1 domain (GB1, 56 residues) folded
with MELD × Amber, written as an edited copy of the PDB entry 3GB1. Start with `1_basic`;
the others are short and assume you have run it.

| Folder | What you learn | Runs in |
|---|---|---|
| [`1_basic/`](1_basic/1_basic.ipynb) | the whole method, step by step in plain language: the folding as a movie, PCA of the trajectory from 50 ns on, HDBSCAN states on PCs 1–3, one medoid per state, each state's three largest modes, a tour of the six new `_mdcif_` categories read back from the file, 3D views, and a glossary | < 1 min |
| [`2_custom_clustering/`](2_custom_clustering/2_custom_clustering.ipynb) | your own clustering instead of HDBSCAN: a scikit-learn-style estimator (a Gaussian mixture) or labels made anywhere (here an RMSD rule) | < 1 min |
| [`3_cli/`](3_cli/3_cli.ipynb) | the same from a terminal: `mdcif build`, `info`, `validate`, `landscape`, and the main options | ~1 min |
| [`4_massweighted_heavy_atoms/`](4_massweighted_heavy_atoms/4_massweighted_heavy_atoms.ipynb) | all C, N and O atoms instead of Cα, and mass-weighted (quasi-harmonic) modes for each state; what changes in the file; the same with `mdcif build --selection ... --mass-weighted` | < 1 min |

Each folder holds its notebook, its `provenance.yaml` (the run conditions that go into the file)
and everything the notebook writes:

```
1_basic/               gb1_meld_300K.cif / .bcif, figures/ (Figs. 2-4: vector PDF + 600-dpi PNG)
2_custom_clustering/   gb1_gmm.cif, gb1_rmsd.cif, labels_rmsd.txt, figures/ (Fig. 1)
3_cli/                 gb1.cif, gb1_landscape.png, gb1_leaf.bcif, gb1_labels.cif
4_massweighted_heavy_atoms/   gb1_heavy_mw.cif
```

The `--labels` example in `3_cli` reads `../2_custom_clustering/labels_rmsd.txt`, so run 2
before 3.

## The simulation

GB1 folded from an extended chain with MELD (30-replica temperature + Hamiltonian exchange,
ff19SB, GB-Neck2 implicit solvent, Amber 24). The mdCIF files summarize the 300 K replica and
are written as an edited copy of 3GB1 (solution NMR): same entity, chain, numbering, source and
citations; the 32 NMR models are replaced by the MD cluster medoids, superposed on NMR model 1.
`provenance.yaml` records the run conditions (from the run's `leap.in`, `meld.mdin.001`,
`meld.groupfile`, `restraints.indxf`).

## Data (not committed)

```
examples/data/3gb1/meld/   3gb1_HMR.prmtop, meld.nc.001 (replica 1 = 300 K, 2 GB), run inputs
```

3GB1 itself is downloaded from RCSB on first use and cached in `~/.cache/mdcif/pdb`. In Amber
MELD runs each `meld.nc.NNN` holds one condition (constant `temp0` / `remd_indices`), so replica
1 needs no demultiplexing; `name.nc.001` is read as NetCDF.

## Running

On HiPerGator, run the notebooks on a compute node, not a login node:

```bash
srun --account=alberto.perezant --qos=alberto.perezant --cpus-per-task=8 --mem=32gb \
     --time=01:00:00 --pty bash
conda activate mdcif          # module load conda first
export PYTHONPATH=            # the Amber module's python3.12 site-packages leak in otherwise
cd examples/create/1_basic
jupyter nbconvert --to notebook --execute --inplace 1_basic.ipynb
```

The whole of `1_basic` as one command (from `examples/create/1_basic/`):

```bash
mdcif build ../../data/3gb1/meld/3gb1_HMR.prmtop ../../data/3gb1/meld/meld.nc.001 \
    --provenance provenance.yaml --start 25000 --step 10 --output gb1_meld_300K.cif
```

`--start 25000 --step 10` = 50–403.65 ns every 20 ps (17,683 frames; ~6 s).

## Result (run stopped at 403.65 of 500 ns)

| | HDBSCAN `eom` (default) | HDBSCAN `leaf` | Gaussian mixture, 3 states |
|---|---|---|---|
| global PCA, % variance PC1/2/3 | 60 / 11 / 8 | same | same |
| clusters (noise) | 2 (0%) | 4 (51%) | 3 (none) |
| model 1 | 95.3% of frames, CA RMSD to 3GB1 2.61 Å, modes 26/20/8% | 26.2%, 2.97 Å | 53.0%, 2.54 Å |
| model 2 | 4.7%, 9.04 Å (misfolded) | 16.3%, 2.78 Å | 42.2%, 3.06 Å |
| models 3–4 | – | 4.7% (misfolded, 9.04 Å), 1.5% (2.63 Å) | 4.7% (misfolded, 9.04 Å) |
| file size | 280 KB (.cif), 148 KB (.bcif) | 180 KB (.bcif) | 315 KB (.cif) |

PC1 separates the misfolded state from the folded basin; PCs 2–3 separate substates inside the
folded basin, which `eom` keeps as one state whose own modes describe the motion between them.
The NMR models project at the edge of the folded basin.

## View

Python, no plugin: [`../view/python/view_mdcif.ipynb`](../view/python/view_mdcif.ipynb) reads
`1_basic/gb1_meld_300K.cif` with `mdcif.reader.read_mdcif` and, with `mdcif.view.MdcifViewer`,
shows the models, draws each state's three modes as porcupine plots and plays them as movies
(py3Dmol). ChimeraX: [`../view/chimerax/`](../view/chimerax/README.md).

PyMOL (`viewers/pymol`), from the repository root:

```
mdcif_load examples/create/1_basic/gb1_meld_300K.cif, gb1
mdcif_morph gb1
mdcif_animate gb1, model=1, mode=1, amplitude=3
mdcif_porcupine gb1, model=1, mode=1, scale=3
mdcif_animate gb1, model=2, mode=1, scope=global
```

ChimeraX (`viewers/chimerax`, after `devel install`):

```
mdcif open examples/create/1_basic/gb1_meld_300K.cif
mdcif morph #1
mdcif animate #1.1 mode 1 amplitude 3
mdcif porcupine #1.1 mode 1 scale 3
open 3gb1          # already superposed: the models sit on NMR model 1
```
