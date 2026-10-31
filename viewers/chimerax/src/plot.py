"""Distance / angle / dihedral plots along a mode animation ('mdcif plot', the MDCIF Plot tool).

A plot is one animation - an 'mdcif animate' copy, one coordinate set per frame - and one kind
of measurement, with up to eight series of that kind (one per set of atoms). Atoms are kept as
(chain, residue number, atom name) keys, so a plot can move to another model of the same file.
The session's plots are in ``plots(session)``: the command and the window (``tool.py``) add to
the same ones, and a window listens to its plot through ``Plot.listeners``.

No Qt or ChimeraX imports here (matplotlib only, imported when drawing): the command saves
images and CSV files without a window, and tests load this module without ChimeraX.
"""

import csv
import os

import numpy as np

from . import __version__
from .mdcif_data import MEASURES

UNITS = {"distance": "A", "angle": "deg", "dihedral": "deg"}          # log and CSV
A_KIND = {"distance": "a distance", "angle": "an angle", "dihedral": "a dihedral"}
AXIS_UNITS = {"distance": "Å", "angle": "°", "dihedral": "°"}
IMAGE_SUFFIXES = (".png", ".svg", ".pdf")
IMAGE_SIZE, IMAGE_DPI = (6.0, 3.5), 300                                 # inches, saved images
TICK_STEPS = (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)                 # frame axis
# Categorical slots of the reference palette in their fixed order (slots 1-3 = SERIES in
# src/mdcif/plot.py, which the viewers may not import); all eight pass the adjacent-pair
# colour-vision checks for lines. A series keeps its colour when others are removed.
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6",
          "#e34948"]
SURFACE, TEXT, TEXT_2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"


PLOT_MADE = "_mdcif_plot_made"          # attribute of copies animated for a plot (not saved)


def atom_label(key):
    chain, number, name = key
    return f"/{chain}:{number}@{name}"


class Series:
    """One measurement: atom keys in order, its colour and its value at every frame."""

    def __init__(self, keys, values, color):
        self.keys = [tuple(str(v) for v in k) for k in keys]
        self.values = np.asarray(values, float)
        self.color = color

    @property
    def label(self):
        return " - ".join(atom_label(k) for k in self.keys)


class Plot:
    """One kind of measurement along the frames of one animate copy.

    ``source``: the model the copy animates (None if unknown), which a new animation of the
    same model starts from. ``own``: the copy was made for a plot (``PLOT_MADE`` set on it), so
    a plot may close it when it animates again; a copy the user made is only stopped and
    hidden. The animation settings are kept when the copy is closed, so the data can still be
    saved.
    """

    def __init__(self, kind, structure, source=None):
        if kind not in MEASURES:
            raise ValueError(f"unknown measurement {kind!r}")
        self.kind, self.series, self.listeners = kind, [], []
        self.set_animation(structure, source)

    def set_animation(self, structure, source=None):
        self.structure, self.source = structure, source
        self.data = structure.mdcif_data
        self.path = self.data.path
        self.model = structure.mdcif_model
        self.scope = structure.mdcif_scope
        self.mode = structure.mdcif_mode
        self.amplitude = structure.mdcif_amplitude
        self.n_frames = structure.num_coordsets
        self.spec = f"#{structure.id_string}"

    @property
    def own(self):
        return bool(getattr(self.structure, PLOT_MADE, False))

    @property
    def alive(self):
        """Whether the animate copy is still open."""
        s = self.structure
        return s is not None and not getattr(s, "deleted", False)

    @property
    def unit(self):
        return UNITS[self.kind]

    def describe(self):
        return (f"{self.scope} mode {self.mode} of model {self.model}, "
                f"+/-{abs(self.amplitude):g} sd, {self.n_frames} frames")

    def add(self, keys, values):
        """A new series in the first free colour; ValueError if all eight are used."""
        if len(keys) != MEASURES[self.kind]:
            raise ValueError(f"{A_KIND[self.kind]} needs {MEASURES[self.kind]} atoms, got "
                             f"{len(keys)}")
        used = {s.color for s in self.series}
        free = [c for c in COLORS if c not in used]
        if not free:
            raise ValueError(f"a plot holds at most {len(COLORS)} series (start another plot "
                             "with 'new true')")
        series = Series(keys, values, free[0])
        self.series.append(series)
        return series

    def changed(self):
        """Tell the windows showing this plot to draw it again."""
        for listener in list(self.listeners):
            listener(self)


def plots(session):
    """The session's plots, oldest first; those whose copy is closed and that no window
    shows any more are dropped."""
    found = getattr(session, "_mdcif_plots", None)
    if found is None:
        found = session._mdcif_plots = []
    found[:] = [p for p in found if p.alive or p.listeners]
    return found


# --------------------------------------------------------------------------- drawing

class FrameMarker:
    """The current-frame marks of a drawn plot: a thin vertical line and a larger dot (series
    colour, dark edge) on each series; ``set(frame)`` moves them.

    The marks are animated artists: a full draw of the figure leaves them out, and the window
    draws them over a saved copy of the rest (blitting), so a new frame costs about 1 ms
    instead of a whole redraw. ``artists`` are the ones to draw.
    """

    def __init__(self, ax, plot, frame):
        self.plot = plot
        self.line = ax.axvline(frame, color=MUTED, lw=0.8, zorder=2, animated=True)
        xy = self._points(frame)                # placed on the data: the y range stays as is
        self.dots = ax.scatter(xy[:, 0], xy[:, 1], s=64, c=[s.color for s in plot.series],
                               edgecolors=TEXT, linewidths=1.2, zorder=4, animated=True)
        self.artists = (self.line, self.dots)

    def _points(self, frame):
        ys = [s.values[frame - 1] if 1 <= frame <= len(s.values) else np.nan
              for s in self.plot.series]
        return np.column_stack([np.full(len(ys), float(frame)), ys])

    def set(self, frame):
        frame = int(frame)
        self.line.set_xdata([frame, frame])
        self.dots.set_offsets(self._points(frame))


def frame_ticks(n_frames):
    """Frame-axis ticks: 1, then every few frames (at most about ten)."""
    step = next((s for s in TICK_STEPS if n_frames / s <= 10), TICK_STEPS[-1])
    return sorted({1, *range(step, n_frames + 1, step)})


def draw(figure, plot, current=None, hint=None):
    """Draw ``plot`` (or None) on a matplotlib figure, which is cleared first.

    Points at every frame joined by a thin line, one colour per series; with two or more
    series a legend above the axes, with one its atoms as the title. ``current``: the frame to
    mark (the window's playing frame; saved images have none; an empty plot has none either);
    ``hint``: text for an empty plot. Returns the ``FrameMarker`` or None.
    """
    figure.clear()
    figure.set_layout_engine("constrained")
    figure.set_facecolor(SURFACE)
    ax = figure.add_subplot()
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=TEXT_2, labelsize=8)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    if plot is None or not plot.series:
        if hint:
            ax.text(0.5, 0.5, hint, transform=ax.transAxes, ha="center", va="center",
                    color=MUTED, fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
    if plot is None:
        return None
    frames = np.arange(1, plot.n_frames + 1)
    for s in plot.series:
        ax.plot(frames[:len(s.values)], s.values, color=s.color, lw=1.2, marker="o", ms=4.5,
                mfc=s.color, mec=SURFACE, mew=0.8, label=s.label, zorder=3)
    ax.set_xlim(0.5, plot.n_frames + 0.5)
    if plot.series:
        ax.set_xticks(frame_ticks(plot.n_frames))
    ax.set_xlabel(f"Frame ({plot.scope} mode {plot.mode} of model {plot.model}, "
                  f"±{abs(plot.amplitude):g} sd, {os.path.basename(plot.path)})",
                  color=TEXT_2, fontsize=8)
    ax.set_ylabel(f"{plot.kind.capitalize()} ({AXIS_UNITS[plot.kind]})", color=TEXT_2,
                  fontsize=9)
    if len(plot.series) == 1:
        ax.set_title(plot.series[0].label, loc="left", color=TEXT, fontsize=9)
    elif plot.series:                           # legend rows read in series order: matplotlib
        n = len(plot.series)                    # fills columns, so hand it column by column
        ncols = min(n, 3 if plot.kind == "distance" else 2)
        order = [k for j in range(ncols) for k in range(j, n, ncols)]
        handles = [ax.get_lines()[k] for k in order]
        figure.legend(handles, [h.get_label() for h in handles], loc="outside upper left",
                      frameon=False, fontsize=8, labelcolor=TEXT, ncols=ncols)
    if current is None or not plot.series:
        return None
    return FrameMarker(ax, plot, current)


def save_image(plot, path):
    """PNG, SVG or PDF (by suffix) of the whole plot, 6 x 3.5 in at 300 dpi, no frame marks."""
    from matplotlib.figure import Figure

    if os.path.splitext(str(path))[1].lower() not in IMAGE_SUFFIXES:
        raise ValueError(f"{path}: the image must be {', '.join(IMAGE_SUFFIXES)}")
    figure = Figure(figsize=IMAGE_SIZE)
    draw(figure, plot)
    figure.savefig(path, dpi=IMAGE_DPI, facecolor=SURFACE)


def header(plot):
    """The '#' lines at the top of a CSV file, without the '# '."""
    return [f"mdcif plot (ChimeraX-MDCIF {__version__})",
            f"file: {os.path.abspath(plot.path)}",
            f"model: {plot.model} ({plot.data.model_label(plot.model)})",
            f"scope: {plot.scope}",
            f"mode: {plot.mode}",
            f"amplitude: {plot.amplitude:g} sd",
            f"frames: {plot.n_frames}",
            f"kind: {plot.kind} ({plot.unit})"] + (
        ["dihedrals unwrapped from frame 1: values can pass +/-180 (subtract or add 360 for "
         "ChimeraX's torsion, -180 to 180)"] if plot.kind == "dihedral" else [])


def write_csv(plot, path):
    """``#`` header lines, then ``frame`` and one column per series named by its atoms
    (pandas: ``read_csv(path, comment="#")``)."""
    with open(path, "w", newline="", encoding="utf-8") as fh:
        for line in header(plot):
            fh.write(f"# {line}\n")
        out = csv.writer(fh)
        out.writerow(["frame"] + [s.label for s in plot.series])
        for k in range(plot.n_frames):
            out.writerow([k + 1] + [f"{s.values[k]:.6g}" for s in plot.series])
