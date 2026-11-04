"""Distance / angle / dihedral plots along a mode animation ('mdcif plot', the MDCIF Plot tool).

A plot is one animation - an 'mdcif animate' copy, one coordinate set per frame - and the
series measured on it, one per set of atoms; the number of atoms gives a series' kind (2 =
distance, 3 = angle, 4 = dihedral). Each kind present gets a panel, stacked in that order on
one frame axis, with up to eight series. Atoms are kept as (chain, residue number, atom name)
keys, so a plot can move to another model of the same file; ``parse_atom`` / ``match_key``
turn typed atoms (``38CA``, ``A:38@CA``) into keys. The session's plots are in
``plots(session)``: the command and the window (``tool.py``) add to the same ones, and a
window listens to its plot through ``Plot.listeners``.

No Qt or ChimeraX imports here (matplotlib only, imported when drawing): the command saves
images and CSV files without a window, and tests load this module without ChimeraX.
"""

import csv
import os
import re

import numpy as np

from . import __version__
from .mdcif_data import MEASURES

UNITS = {"distance": "A", "angle": "deg", "dihedral": "deg"}          # log and CSV
A_KIND = {"distance": "a distance", "angle": "an angle", "dihedral": "a dihedral"}
AXIS_UNITS = {"distance": "Å", "angle": "°", "dihedral": "°"}
KIND_OF = {n: k for k, n in MEASURES.items()}                           # number of atoms -> kind
IMAGE_SUFFIXES = (".png", ".svg", ".pdf")
IMAGE_SIZE, IMAGE_DPI = (6.0, 3.5), 300                 # inches (one panel), saved images
PANEL_HEIGHT = 2.2                                      # inches added by each further panel
PAD = 0.12                                              # inches of white around the plot
LEGEND_SIZE, HANDLE, HANDLE_PAD, COLUMN_GAP = 8, 1.8, 0.8, 1.4      # points; widths in em
MAX_COLUMNS, LEFT_ROOM = 4, 1.0                 # legend columns; inches beside a panel's box
TICK_STEPS = (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)                 # frame axis
MINOR_STEPS = dict(zip(TICK_STEPS, (None, 1, 1, 2, 5, 10, 20, 50, 100, 200), strict=True))
# Categorical slots of the reference palette in their fixed order (slots 1-3 = SERIES in
# src/mdcif/plot.py, which the viewers may not import); all eight pass the adjacent-pair
# colour-vision checks for lines, also on FACE. A series keeps its colour when others are
# removed; each panel has its own eight.
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6",
          "#e34948"]
PAPER, FACE, GRID = "#ffffff", "#efefed", "#ffffff"     # margins, plot area, grid lines
TEXT, TEXT_2, MUTED, FRAME = "#0b0b0b", "#52514e", "#898781", "#52514e"


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
    def kind(self):
        return KIND_OF[len(self.keys)]

    @property
    def label(self):
        return " - ".join(atom_label(k) for k in self.keys)


class Plot:
    """Distances, angles and dihedrals along the frames of one animate copy.

    ``source``: the model the copy animates (None if unknown), which a new animation of the
    same model starts from. ``own``: the copy was made for a plot (``PLOT_MADE`` set on it), so
    a plot may close it when it animates again; a copy the user made is only stopped and
    hidden. The animation settings are kept when the copy is closed, so the data can still be
    saved.
    """

    def __init__(self, structure, source=None):
        self.series, self.listeners = [], []
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
    def kinds(self):
        """The kinds of the series, in panel order (distance, angle, dihedral)."""
        present = {s.kind for s in self.series}
        return [k for k in MEASURES if k in present]

    def series_of(self, kind):
        return [s for s in self.series if s.kind == kind]

    def describe(self):
        return (f"{self.scope} mode {self.mode} of model {self.model}, "
                f"+/-{abs(self.amplitude):g} sd, {self.n_frames} frames")

    def add(self, keys, values):
        """A new series (its kind from the number of atoms) in the first colour free in its
        panel; ValueError if all eight are used there."""
        kind = KIND_OF.get(len(keys))
        if kind is None:
            raise ValueError(f"a distance, angle or dihedral needs 2, 3 or 4 atoms, got "
                             f"{len(keys)}")
        used = {s.color for s in self.series_of(kind)}
        free = [c for c in COLORS if c not in used]
        if not free:
            raise ValueError(f"a plot holds at most {len(COLORS)} {kind}s (start another plot "
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


# --------------------------------------------------------------------------- typed atoms

_ATOM_TEXT = re.compile(r"""\s*(?:/?(?P<chain>[^\s/:@]*):)?\s*     # /A: or A: or :
                            (?P<number>-?\d+)\s*                    # residue number
                            (?:@\s*(?P<at>[^\s/:@]+)                # @name: any name
                            |(?P<name>[A-Za-z][^\s/:@]*))\s*        # name: starts with a letter
                        """, re.VERBOSE)


def parse_atom(text):
    """``(chain or None, residue number, ATOM NAME)`` from typed text: ``38CA``, ``38 CA``,
    ``38@CA``, ``:38@CA``, ``A:38CA``, ``/A:38@CA``. The name is upper-cased; one that starts
    with a digit needs the ``@`` (``38@1HB``). ValueError if the text is not an atom."""
    found = _ATOM_TEXT.fullmatch(text)
    if found is None:
        what = f"{text.strip()!r} is not an atom" if text.strip() else "no atom typed"
        raise ValueError(f"{what}: type a residue number and atom name, e.g. 38CA or 38@CA "
                         "(A:38CA for chain A)")
    chain = found["chain"] or None
    return chain, int(found["number"]), (found["at"] or found["name"]).upper()


def match_key(parsed, keys):
    """The key among a structure's (chain, number, name) ``keys`` that ``parsed`` (from
    ``parse_atom``) names. Without a chain, the atom must be in one chain only; a chain or
    name that matches only when ignoring case is taken too. ValueError naming what is
    missing or ambiguous."""
    chain, number, name = parsed
    number = str(number)
    residue = [tuple(str(v) for v in k) for k in keys if str(k[1]) == number]
    if chain is not None:
        residue = ([k for k in residue if k[0] == chain]
                   or [k for k in residue if k[0].upper() == chain.upper()])
    if not residue:
        raise ValueError(f"no residue {'/' + chain if chain else ''}:{number}")
    hits = list(dict.fromkeys([k for k in residue if k[2] == name]
                              or [k for k in residue if k[2].upper() == name.upper()]))
    if not hits:
        chains = list(dict.fromkeys(k[0] for k in residue))
        where = f"/{chains[0]}:{number}" if len(chains) == 1 else f":{number}"
        names = list(dict.fromkeys(k[2] for k in residue))
        more = ", ..." if len(names) > 12 else ""
        raise ValueError(f"no atom {where}@{name} (residue {number} has "
                         f"{', '.join(names[:12])}{more})")
    found = list(dict.fromkeys(k[0] for k in hits))
    if len(found) > 1:
        at = "@" if name[0].isdigit() else ""
        raise ValueError(f"{number}{at}{name} is in chains {', '.join(found)}: write "
                         f"{found[0]}:{number}{at}{name} for chain {found[0]}")
    return hits[0]


def digit_hint(text, parsed, keys):
    """What to add to the error for typed ``text`` that names no atom, when it has no ``@``
    and its residue number may end with the first digit of the atom name ('381HB' read as
    residue 381, atom HB): that reading if it is an atom, else a reminder if the residue is
    missing too; otherwise '' (also when no atom name in ``keys`` starts with a digit)."""
    chain, number, name = parsed
    digits = str(number)
    if ("@" in text or len(digits.lstrip("-")) < 2
            or not any(k[2][:1].isdigit() for k in keys)):
        return ""
    other = (chain, int(digits[:-1]), digits[-1] + name)
    try:
        match_key(other, keys)
    except ValueError:
        if any(str(k[1]) == digits for k in keys):      # the residue is there: another mistake
            return ""
        return "; a name that starts with a digit needs @, e.g. 38@1HB"
    return (f"; did you mean {chain + ':' if chain else ''}{other[1]}@{other[2]}? (a name that "
            "starts with a digit needs @)")


# --------------------------------------------------------------------------- drawing

class FrameMarker:
    """The current-frame marks of a drawn plot: in each panel a thin vertical line and a
    larger dot (series colour, dark edge) on each series; ``set(frame)`` moves them all.

    The marks are animated artists: a full draw of the figure leaves them out, and the window
    draws them over a saved copy of the rest (blitting), so a new frame costs about 1 ms
    instead of a whole redraw. Each dot is a one-point Line2D: a scatter with one point takes
    about 1.8 ms to draw (matplotlib 3.10, Agg), and an angle or dihedral is often alone in its
    panel. ``artists`` are the ones to draw; ``panels`` maps each kind to its (line, [dot of
    each series]), ``line`` / ``dots`` are the top panel's, and ``points(kind)`` gives where
    a panel's dots are.
    """

    def __init__(self, axes, plot, frame):            # axes: {kind: matplotlib Axes}
        self.plot, self.panels, self._series = plot, {}, {}
        for kind, ax in axes.items():
            self._series[kind] = series = plot.series_of(kind)
            line = ax.axvline(frame, color=MUTED, lw=0.8, zorder=2, animated=True)
            dots = [ax.plot([frame], [self._y(s, frame)], "o", ms=8, mfc=s.color, mec=TEXT,
                            mew=1.2, zorder=4, animated=True)[0]      # on the data: the y
                    for s in series]                                  # range stays as is
            self.panels[kind] = (line, dots)
        self.line, self.dots = next(iter(self.panels.values()))
        self.artists = tuple(a for line, dots in self.panels.values() for a in (line, *dots))

    @staticmethod
    def _y(series, frame):
        return series.values[frame - 1] if 1 <= frame <= len(series.values) else np.nan

    def set(self, frame):
        frame = int(frame)
        for kind, (line, dots) in self.panels.items():
            line.set_xdata([frame, frame])
            for dot, s in zip(dots, self._series[kind], strict=True):
                dot.set_data([frame], [self._y(s, frame)])

    def points(self, kind):
        """``(n, 2)`` frame and value of each dot in the ``kind`` panel."""
        dots = self.panels[kind][1]
        return np.array([[d.get_xdata()[0], d.get_ydata()[0]] for d in dots], float)


def _frame_step(n_frames):
    return next((s for s in TICK_STEPS if n_frames / s <= 10), TICK_STEPS[-1])


def frame_ticks(n_frames):
    """Frame-axis ticks: 1, then every few frames (at most about ten)."""
    step = _frame_step(n_frames)
    return sorted({1, *range(step, n_frames + 1, step)})


def frame_minor_step(n_frames):
    """Frames between minor ticks (None when every frame has a major tick). Not matplotlib's
    AutoMinorLocator: the major ticks start at 1, so their first gap is uneven."""
    return MINOR_STEPS[_frame_step(n_frames)]


def figure_size(n_panels):
    """Image size in inches: 6 x 3.5 for one panel, 2.2 in taller for each further one."""
    width, height = IMAGE_SIZE
    return width, height + PANEL_HEIGHT * (max(n_panels, 1) - 1)


def _style(ax):
    """A panel: light gray plot area with white grid lines, closed dark frame, ticks inward
    on all four sides."""
    ax.set_facecolor(FACE)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color(FRAME)
        spine.set_linewidth(0.8)
    ax.tick_params(which="both", direction="in", top=True, right=True, color=FRAME,
                   labelcolor=TEXT_2, labelsize=8)
    ax.tick_params(which="major", length=4, width=0.8)
    ax.tick_params(which="minor", length=2.2, width=0.6)
    ax.grid(which="major", color=GRID, lw=0.8)
    ax.set_axisbelow(True)


def legend_layout(figure, plot):
    """Legend columns of each panel with two or more series, ``{kind: n}``: as many as fit
    across the panel at the figure's width (the widest label sets the column width), at most
    four. A legend wider than its panel would make constrained layout narrow every panel."""
    from matplotlib.font_manager import FontProperties
    from matplotlib.textpath import text_to_path

    prop, em = FontProperties(size=LEGEND_SIZE), LEGEND_SIZE / 72      # em in inches
    room = figure.get_figwidth() - LEFT_ROOM - PAD
    layout = {}
    for kind in plot.kinds:
        labels = [s.label for s in plot.series_of(kind)]
        if len(labels) < 2:
            continue
        widest = max(text_to_path.get_text_width_height_descent(t, prop, ismath=False)[0]
                     for t in labels) / 72
        column, gap = widest + (HANDLE + HANDLE_PAD) * em, COLUMN_GAP * em
        layout[kind] = max(1, min(MAX_COLUMNS, len(labels), int((room + gap) // (column + gap))))
    return layout


def _legend(ax, lines, ncols):
    """Several series of a panel: their atoms in a legend above it, read row by row in series
    order (matplotlib fills columns, so it gets them column by column)."""
    n = len(lines)
    order = [k for j in range(ncols) for k in range(j, n, ncols)]
    handles = [lines[k] for k in order]
    ax.legend(handles, [h.get_label() for h in handles], loc="lower left",
              bbox_to_anchor=(0.0, 1.0), borderaxespad=0.2, borderpad=0.0, frameon=False,
              fontsize=LEGEND_SIZE, labelcolor=TEXT, ncols=ncols, handlelength=HANDLE,
              handletextpad=HANDLE_PAD, columnspacing=COLUMN_GAP)


def draw(figure, plot, current=None, hint=None):
    """Draw ``plot`` (or None) on a matplotlib figure, which is cleared first.

    One panel per kind present (distances, then angles, then dihedrals) on a shared frame
    axis; points at every frame joined by a thin line, one colour per series; in a panel with
    two or more series a legend above it, with one its atoms as the panel title. The
    animation is named in a small line at the top left. ``current``: the frame to mark in
    every panel (the window's playing frame; saved images have none; an empty plot has none
    either); ``hint``: text for an empty plot. Returns the ``FrameMarker`` or None.
    """
    from matplotlib.ticker import AutoMinorLocator, MultipleLocator
    from matplotlib.transforms import ScaledTranslation

    figure.clear()
    figure.set_layout_engine("constrained", w_pad=PAD, h_pad=PAD, hspace=0.03)
    figure.set_facecolor(PAPER)
    kinds = plot.kinds if plot is not None else []
    axes = list(figure.subplots(max(len(kinds), 1), 1, sharex=True, squeeze=False)[:, 0])
    for ax in axes:
        _style(ax)
    bottom = axes[-1]
    if not kinds:
        if hint:
            bottom.text(0.5, 0.5, hint, transform=bottom.transAxes, ha="center", va="center",
                        color=MUTED, fontsize=9)
        bottom.set_xticks([])
        bottom.set_yticks([])
    if plot is None:
        return None
    frames = np.arange(1, plot.n_frames + 1)
    panels = dict(zip(kinds, axes, strict=False))     # (empty for an empty plot)
    columns = legend_layout(figure, plot)
    for kind, ax in panels.items():
        series = plot.series_of(kind)
        lines = [ax.plot(frames[:len(s.values)], s.values, color=s.color, lw=1.2, marker="o",
                         ms=4.5, mfc=s.color, mec=FACE, mew=0.8, label=s.label, zorder=3)[0]
                 for s in series]
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.set_ylabel(f"{kind.capitalize()} ({AXIS_UNITS[kind]})", color=TEXT, fontsize=9,
                      fontweight="bold")
        if len(series) == 1:
            ax.set_title(series[0].label, loc="left", color=TEXT, fontsize=9)
        else:
            _legend(ax, lines, columns[kind])
    bottom.set_xlim(0.5, plot.n_frames + 0.5)       # shared by every panel
    if kinds:
        bottom.set_xticks(frame_ticks(plot.n_frames))
        minor = frame_minor_step(plot.n_frames)
        if minor:
            bottom.xaxis.set_minor_locator(MultipleLocator(minor))
    bottom.set_xlabel("Frame", color=TEXT, fontsize=9, fontweight="bold")
    if len(axes) > 1:
        figure.align_ylabels(axes)
    left = figure.transFigure + ScaledTranslation(PAD, 0, figure.dpi_scale_trans)
    figure.suptitle(f"{plot.scope} mode {plot.mode} of model {plot.model}, "
                    f"±{abs(plot.amplitude):g} sd, {os.path.basename(plot.path)}", x=0.0,
                    ha="left", transform=left, color=TEXT_2, fontsize=8)
    if current is None or not kinds:
        return None
    return FrameMarker(panels, plot, current)


def save_image(plot, path):
    """PNG, SVG or PDF (by suffix) of the whole plot, every panel, at 300 dpi, no frame
    marks: 6 in wide, 3.5 in high for one panel and 2.2 in more for each further one."""
    from matplotlib.figure import Figure

    if os.path.splitext(str(path))[1].lower() not in IMAGE_SUFFIXES:
        raise ValueError(f"{path}: the image must be {', '.join(IMAGE_SUFFIXES)}")
    figure = Figure(figsize=figure_size(len(plot.kinds)))
    draw(figure, plot)
    figure.savefig(path, dpi=IMAGE_DPI, facecolor=PAPER)


def header(plot):
    """The '#' lines at the top of a CSV file, without the '# '."""
    kinds = ", ".join(f"{k} ({UNITS[k]})" for k in plot.kinds) or "none"
    return [f"mdcif plot (ChimeraX-MDCIF {__version__})",
            f"file: {os.path.abspath(plot.path)}",
            f"model: {plot.model} ({plot.data.model_label(plot.model)})",
            f"scope: {plot.scope}",
            f"mode: {plot.mode}",
            f"amplitude: {plot.amplitude:g} sd",
            f"frames: {plot.n_frames}",
            f"kind: {kinds}"] + (
        ["dihedrals unwrapped from frame 1: values can pass +/-180 (subtract or add 360 for "
         "ChimeraX's torsion, -180 to 180)"] if "dihedral" in plot.kinds else [])


def write_csv(plot, path):
    """``#`` header lines, then ``frame`` and one column per series named by its atoms (2
    atoms = a distance, 3 = an angle, 4 = a dihedral; pandas: ``read_csv(path,
    comment="#")``)."""
    with open(path, "w", newline="", encoding="utf-8") as fh:
        for line in header(plot):
            fh.write(f"# {line}\n")
        out = csv.writer(fh)
        out.writerow(["frame"] + [s.label for s in plot.series])
        for k in range(plot.n_frames):
            out.writerow([k + 1] + [f"{s.values[k]:.6g}" for s in plot.series])
