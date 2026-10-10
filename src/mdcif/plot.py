"""Conformational landscape: frames on the global PCs, coloured by cluster.

The frames stored in ``_mdcif_projection`` are plotted on global PCs 1-3,
coloured by the deposited model whose cluster they belong to. Each model's
own coordinates are projected with the global origin and vectors (so the
medoids land inside their clusters), and optionally the models of a reference
PDB entry too.
"""

import numpy as np

from .modes import project
from .reader import read_mdcif
from .trajectory import superpose

# Reference palette slots 1-3 (validated all-pairs for scatter); clusters beyond the third
# fold into a neutral "other", noise is a light neutral.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]
OTHER, NOISE = "#8c8b86", "#d9d8d3"
SURFACE, TEXT, TEXT_2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e0"
# medoid label offsets (points), alternated so nearby medoids don't print on top of each other
LABEL_OFFSETS = [(6, 6), (-22, -14), (6, -14), (-22, 6)]


def _atoms(model):
    return list(model.get_atoms())


def model_projections(data):
    """``{model_id: (n_modes,)}`` projection of every deposited model on the global modes."""
    g = data.global_modes
    out = {}
    for _, model in data.system._all_models():
        xyz = np.array([[a.x, a.y, a.z] for a in _atoms(model)])[g.atom_indices]
        out[int(model._id)] = project(superpose(xyz[None], g.origin)[0], g)[0]
    return out


def reference_projections(data, reference):
    """Projection of every model of a reference entry onto the global modes.

    Chains are matched by sequence (as in ``mdcif build --reference``), atoms by residue
    and name.
    """
    from .reference import entity_sequence, load_reference, match_chains

    g = data.global_modes
    first = next(m for _, m in data.system._all_models() if int(m._id) == g.model_id)
    atoms = _atoms(first)
    ref = load_reference(reference)
    cmap = match_chains(ref, {a.asym_unit._id: entity_sequence(a.asym_unit.entity)
                              for a in atoms})
    keys = []
    for i in g.atom_indices:
        asym, offset = cmap[atoms[i].asym_unit._id]
        keys.append((id(asym), offset + atoms[i].seq_id, atoms[i].atom_id))
    out = []
    for _, model in ref.system._all_models():
        xyz = {}
        for a in model.get_atoms():
            xyz.setdefault((id(a.asym_unit), a.seq_id, a.atom_id), (a.x, a.y, a.z))
        if not all(k in xyz for k in keys):
            raise ValueError(f"{ref.label} lacks some of the {len(keys)} atoms of the modes")
        coords = np.array([xyz[k] for k in keys])
        out.append(project(superpose(coords[None], g.origin)[0], g)[0])
    return ref.label, np.array(out)


def landscape_figure(coords, model_ids, medoids, variance_fraction, title="", subtitle="",
                     ref_label=None, ref_proj=None):
    """PC1-PC2 / PC1-PC3 / PC2-PC3 scatter of frames coloured by model; returns the Figure.

    ``coords`` ``(n, >=2)`` frame coordinates on the PCs; ``model_ids`` ``(n,)`` the
    model (1-based) of each frame's cluster, 0 for noise; ``medoids`` ``{model: (k,)}``.
    Works in notebooks and headless (no pyplot state, no backend switch).
    """
    import matplotlib
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D

    coords, model_ids = np.asarray(coords), np.asarray(model_ids)
    models = sorted(medoids)
    n_pc = min(3, coords.shape[1])
    pairs = [(0, 1), (0, 2), (1, 2)][: {1: 0, 2: 1, 3: 3}[n_pc]]
    if not pairs:
        raise ValueError("need at least 2 PCs to plot")
    # drawing order, back to front: noise, models beyond the third, models 1-3
    groups = []
    if np.any(model_ids == 0):
        groups.append((0, NOISE, f"noise / undeposited ({100 * np.mean(model_ids == 0):.0f}%)"))
    rest = models[len(SERIES):]
    if rest:
        share = 100 * np.isin(model_ids, rest).mean()
        which = f"model {rest[0]}" if len(rest) == 1 else f"models {rest[0]}-{rest[-1]}"
        groups.append((rest, OTHER, f"{which} ({share:.0f}%)"))
    for i, m in enumerate(models[: len(SERIES)]):
        groups.append((m, SERIES[i], f"model {m} ({100 * np.mean(model_ids == m):.0f}%)"))
    var = 100 * np.asarray(variance_fraction)
    style = {"font.size": 9, "axes.edgecolor": TEXT_2, "axes.labelcolor": TEXT_2,
             "xtick.color": TEXT_2, "ytick.color": TEXT_2}
    with matplotlib.rc_context(style):
        fig = Figure(figsize=(4.2 * len(pairs), 4.4), facecolor=SURFACE)
        axes = fig.subplots(1, len(pairs), squeeze=False)[0]
        for ax, (a, b) in zip(axes, pairs, strict=True):
            ax.set_facecolor(SURFACE)
            ax.grid(color=GRID, linewidth=0.6)
            ax.set_axisbelow(True)
            for key, color, label in groups:
                sel = np.isin(model_ids, key) if isinstance(key, list) else model_ids == key
                if sel.any():
                    ax.scatter(coords[sel, a], coords[sel, b], s=4, color=color,
                               alpha=0.55 if key else 0.8, linewidths=0, label=label,
                               rasterized=True)
            if ref_proj is not None:
                ax.scatter(ref_proj[:, a], ref_proj[:, b], s=22, marker="x", color=TEXT,
                           linewidths=1.0, label=f"{ref_label} models")
            for i, m in enumerate(models):
                x, y = medoids[m][a], medoids[m][b]
                color = SERIES[i] if i < len(SERIES) else OTHER
                ax.scatter([x], [y], s=70, color=color, edgecolors=SURFACE, linewidths=2,
                           zorder=5)
                ax.annotate(f"M{m}", (x, y), xytext=LABEL_OFFSETS[i % len(LABEL_OFFSETS)],
                            textcoords="offset points", color=TEXT, fontsize=9,
                            fontweight="bold", zorder=6)
            ax.set_xlabel(f"PC{a + 1} ({var[a]:.0f}% var) / Å")
            ax.set_ylabel(f"PC{b + 1} ({var[b]:.0f}% var) / Å")
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)

        def dot(color, label, marker="o", size=7):
            return Line2D([], [], marker=marker, linestyle="", color=color, markersize=size,
                          label=label)

        # legend: models 1-3, then the other models, then noise
        n_back = len(groups) - min(len(models), len(SERIES))
        legend = [dot(c, label) for _, c, label in groups[n_back:] + groups[:n_back][::-1]]
        if ref_proj is not None:
            legend.append(dot(TEXT, f"{ref_label} models", marker="x", size=6))
        fig.legend(handles=legend, loc="lower center", ncol=len(legend), frameon=False,
                   labelcolor=TEXT, bbox_to_anchor=(0.5, -0.01))
        fig.suptitle(title or "Conformational landscape on the global PCs", x=0.01,
                     ha="left", color=TEXT, fontsize=12, fontweight="bold")
        if subtitle:
            fig.text(0.01, 0.905, subtitle, color=TEXT_2, fontsize=9)
        fig.tight_layout(rect=(0, 0.07, 1, 0.9))
    return fig


def plot_landscape(path, output=None, reference=None):
    """Landscape of an mdcif file from the file alone (``reference``: also its models).

    Saves to ``output`` if given (and returns the path), else returns the Figure.
    """
    data = read_mdcif(path)
    p, g = data.projection, data.global_modes
    if p is None or g is None:
        raise ValueError(f"{path} has no global modes / frame projections")
    ref_label, ref_proj = (reference_projections(data, reference) if reference
                           else (None, None))
    fig = landscape_figure(
        p["coords"], p["model_id"], model_projections(data), g.variance_fraction,
        subtitle=f"{data.system.title} — {len(p['frame'])} frames stored, {g.n_frames} "
                 "analysed; M = cluster medoid (deposited model)",
        ref_label=ref_label, ref_proj=ref_proj)
    if output is None:
        return fig
    fig.savefig(output, dpi=150, facecolor=SURFACE, bbox_inches="tight")
    return output
