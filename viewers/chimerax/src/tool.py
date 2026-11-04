"""The MDCIF Plot window: distances, angles and dihedrals along a mode animation (GUI only).

Opened from Tools > Structure Analysis > MDCIF Plot (empty) or by 'mdcif plot' with the GUI;
imported only then (Qt). The plot itself (state, drawing, files) is in plot.py and the
ChimeraX work in cmd.py, so the window and the command share their plots:

    Model [1: cluster 1] Scope [cluster] Mode [1] Amplitude [3.00 sd] Frames [30]
    Measure [Distance] Atoms [10CA] [40CA] [Add] [Add selected atoms]   [Remove] [Remove all]
    (series: colour, kind and atoms of each)
    (the plot: a panel per kind - distances, then angles, then dihedrals - on one frame axis,
     points at every frame, the current frame marked in each; click or drag = stop there)
    [Play] [Stop]  frame 12 / 30: 8.43, 10.12 Å · 141.3°     [Save plot...] [Save data...]

Measure sets what Add and Add selected atoms measure, and shows that many atom boxes (2 for a
distance, 3 for an angle, 4 for a dihedral); each kind gets its own panel, so a plot can show
all three. Atoms are typed as residue number and atom name (38CA, 38 CA, 38@CA; A:38CA when
the number is in several chains) and looked up on the plot's animation (on its model if the
animation was closed, which is then animated again), or in an empty window on the first shown
mdcif structure; the boxes keep their text. Add selected atoms takes the selection in picking
order.

Changing Model, Scope, Mode, Amplitude or Frames animates again (``cmd.reanimate``) and
measures every series anew. An empty window takes the first atoms' animation: atoms on an
animate copy use that copy (the controls then show its settings), atoms on a model are
animated with the controls' settings. Play plays from frame 1 (a 'coordset' loop can't
start mid-range). Plot windows are not saved in sessions.

While the motion plays only the frame marks are redrawn: they are animated matplotlib
artists (one-point lines, not scatters), drawn over a copy of the rest of the canvas taken
after each full draw, and only the panels are repainted; frame changes that arrive before Qt
runs are shown once. The values line under the plot is cut to the window's width (the whole
of it is its tooltip), so it never widens the window.

For scripts (``runscript``, ``mdcif.chimerax.ChimeraX``): ``windows(session)`` lists the open
windows. Each has ``plot`` (a ``plot.Plot``: ``series``, ``kinds``, ``structure``, ...),
``figure`` / ``canvas`` (matplotlib; ``figure.axes`` = the panels, top to bottom),
``go_to_frame(k)`` (what a click does), ``play()``, ``stop()``, ``add_atoms(atoms)``,
``add_text(*texts)`` (typed atoms: ``add_text("10CA", "40CA")``), ``set_settings(model=,
scope=, mode=, amplitude=, frames=)``, ``set_kind(kind)``, ``save_plot(path)``,
``save_data(path)``, ``shown_frame`` and ``canvas_position(frame, series=0)`` (where a click
on that frame of ``plot.series[series]`` lands, canvas pixels from top left).
"""

import os
from contextlib import contextmanager

from chimerax.core.commands import run
from chimerax.core.errors import UserError
from chimerax.core.tools import ToolInstance
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from Qt.QtCore import Qt, QTimer
from Qt.QtGui import QColor, QIcon, QPixmap
from Qt.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QSizePolicy,
    QSpinBox,
    QVBoxLayout,
)

from . import attach, cmd
from .mdcif_data import MEASURES
from .plot import (
    A_KIND,
    IMAGE_SIZE,
    KIND_OF,
    Plot,
    digit_hint,
    draw,
    legend_layout,
    match_key,
    parse_atom,
    plots,
)

TOOL_NAME = "MDCIF Plot"
KINDS = list(MEASURES)
EXAMPLES = {"distance": "10CA 40CA", "angle": "10CA 11CA 12CA",
            "dihedral": "10CA 11CA 12CA 13CA"}
INFO_UNITS = {"distance": " Å", "angle": "°", "dihedral": "°"}       # 8.43 Å, 141.3°
ATOM_TIP = ("Residue number and atom name: 38CA, 38 CA or 38@CA (38@1HB when the name starts "
            "with a digit); A:38CA for chain A when the number is in several chains. Return "
            "adds.")
IMAGE_FILTERS = ["PNG image (*.png)", "SVG image (*.svg)", "PDF document (*.pdf)"]
DATA_FILTERS = ["CSV file (*.csv)"]
ERROR_COLOR, NOTE_COLOR = "#d03b3b", "#52514e"
CANVAS_HEIGHT, PANEL_PIXELS = 200, 150          # minimum canvas height, + per further panel


def windows(session):
    """The open MDCIF Plot windows."""
    return session.tools.find_by_class(MdcifPlot)


def show_plot(session, plot):
    """Show ``plot`` in its window: the one showing it, else an empty one, else a new one."""
    found = windows(session)
    tool = (next((t for t in found if t.plot is plot), None)
            or next((t for t in found if t.plot is None), None))
    if tool is None:
        tool = MdcifPlot(session, plot=plot)
    elif tool.plot is not plot:
        tool.set_plot(plot)
    tool.display(True)
    return tool


class _ElidedLabel(QLabel):
    """A one-line label that never sets the window's width: text too long for it is cut
    with an ellipsis, and the whole text is its tooltip (and what ``text()`` returns)."""

    def __init__(self):
        super().__init__()
        self._full = ""
        self.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)

    def setText(self, text):                            # Qt's names: overridden
        self._full = text
        self.setToolTip(text)
        self._elide()

    def text(self):
        return self._full

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._elide()

    def _elide(self):
        super().setText(self.fontMetrics().elidedText(self._full, Qt.TextElideMode.ElideRight,
                                                      max(self.width(), 0)))


class MdcifPlot(ToolInstance):
    SESSION_SAVE = False                # the animate copy is saved in sessions, the window not

    def __init__(self, session, tool_name=TOOL_NAME, plot=None):
        from chimerax.core.models import END_DELETE_MODELS
        from chimerax.ui import MainToolWindow

        super().__init__(session, tool_name)
        self.plot = None
        self.kind = KINDS[0]                            # what Add / Add selected atoms measure
        # settings for the first animation of an empty window: mode, scope, amplitude, frames
        self.settings = dict(zip(("mode", "scope", "amplitude", "frames"),
                                 cmd.ANIMATE_DEFAULTS, strict=True))
        self.shown_frame = None                         # frame the marks are on
        self._marker = self._watched = self._watch_handler = None
        self._background = None                         # the canvas without the frame marks
        self._columns = {}                              # legend columns of the last full draw
        self._boxes = []                                # the panels' boxes: what blits repaint
        self._frame_due = False                         # a frame change waits for Qt
        self._error = self._notice = ""                 # red note; grey one, till next action
        self.tool_window = MainToolWindow(self)
        self._build(self.tool_window.ui_area)
        self._delete_handler = session.triggers.add_handler(END_DELETE_MODELS,
                                                            self._models_deleted)
        self.tool_window.manage(placement=None)         # floating: a plot needs the width
        self.set_plot(plot)

    # ------------------------------------------------------------------ layout
    def _build(self, area):
        layout = QVBoxLayout(area)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(4)

        row = QHBoxLayout()
        self.model_box = QComboBox()                    # cluster names are long: cut them short
        self.model_box.setSizeAdjustPolicy(
            QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.model_box.setMinimumContentsLength(14)
        self.scope_box = QComboBox()
        self.mode_box = QSpinBox()
        self.amplitude_box = QDoubleSpinBox()
        self.amplitude_box.setRange(-20.0, 20.0)
        self.amplitude_box.setSingleStep(0.5)
        self.amplitude_box.setDecimals(2)
        self.amplitude_box.setSuffix(" sd")
        self.frames_box = QSpinBox()
        self.frames_box.setRange(cmd.MIN_ANIMATE_FRAMES, 1000)
        for label, widget in (("Model", self.model_box), ("Scope", self.scope_box),
                              ("Mode", self.mode_box), ("Amplitude", self.amplitude_box)):
            row.addWidget(QLabel(label))
            row.addWidget(widget)
        row.addSpacing(6)
        row.addWidget(QLabel("Frames"))
        row.addWidget(self.frames_box)
        row.addStretch(1)
        layout.addLayout(row)
        for box in (self.mode_box, self.amplitude_box, self.frames_box):
            box.setKeyboardTracking(False)              # typed values count on Enter / leaving
            box.valueChanged.connect(self._setting_changed)
        for box in (self.model_box, self.scope_box):
            box.currentIndexChanged.connect(self._setting_changed)

        row = QHBoxLayout()
        self.kind_box = QComboBox()
        self.kind_box.addItems([k.capitalize() for k in KINDS])
        self.kind_box.setToolTip("What Add and Add selected atoms measure; each kind gets its "
                                 "own panel")
        self.kind_box.currentIndexChanged.connect(
            lambda i: self._guard(self.set_kind, KINDS[i]))
        self.atom_boxes = []                            # as many shown as Measure needs
        for _ in range(max(MEASURES.values())):
            box = QLineEdit()
            box.setPlaceholderText("e.g. 38CA")
            box.setToolTip(ATOM_TIP)
            box.setFixedWidth(box.fontMetrics().horizontalAdvance("e.g. 38CA") + 16)
            box.returnPressed.connect(lambda: self._guard(self.add_typed))
            self.atom_boxes.append(box)
        self.add_typed_button = QPushButton("Add")
        self.add_typed_button.setToolTip("Adds the atoms typed in the boxes")
        self.add_typed_button.clicked.connect(lambda: self._guard(self.add_typed))
        self.add_button = QPushButton("Add selected atoms")
        self.add_button.setToolTip("Adds the selected atoms, in picking order. Show the atoms "
                                   "first: a click on the cartoon selects whole residues.")
        self.add_button.clicked.connect(lambda: self._guard(self.add_selected))
        self.remove_button = QPushButton("Remove")
        self.remove_button.clicked.connect(lambda: self._guard(self.remove_selected))
        self.clear_button = QPushButton("Remove all")
        self.clear_button.clicked.connect(lambda: self._guard(self.remove_all))
        row.addWidget(QLabel("Measure"))
        row.addWidget(self.kind_box)
        row.addSpacing(6)
        row.addWidget(QLabel("Atoms"))
        for box in self.atom_boxes:
            row.addWidget(box)
        row.addWidget(self.add_typed_button)
        row.addWidget(self.add_button)
        row.addStretch(1)
        row.addWidget(self.remove_button)
        row.addWidget(self.clear_button)
        layout.addLayout(row)

        self.series_list = QListWidget()
        self.series_list.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.series_list.setMaximumHeight(4 * self.series_list.fontMetrics().height() + 12)
        layout.addWidget(self.series_list)
        self.note = QLabel()
        self.note.setWordWrap(True)
        self.note.hide()
        layout.addWidget(self.note)

        self.figure = Figure(figsize=IMAGE_SIZE, dpi=100)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setMinimumSize(320, CANVAS_HEIGHT)
        layout.addWidget(self.canvas, stretch=1)
        self._mouse = [self.canvas.mpl_connect("button_press_event", self._mouse_event),
                       self.canvas.mpl_connect("motion_notify_event", self._mouse_event),
                       self.canvas.mpl_connect("draw_event", self._drawn),
                       self.canvas.mpl_connect("resize_event", self._resized)]

        row = QHBoxLayout()
        self.play_button = QPushButton("Play")
        self.play_button.setToolTip("Play the animation in a loop, from frame 1")
        self.play_button.clicked.connect(lambda: self._guard(self.play))
        self.stop_button = QPushButton("Stop")
        self.stop_button.setToolTip("Stop the animation (a click on the plot also stops it, "
                                    "at that frame)")
        self.stop_button.clicked.connect(lambda: self._guard(self.stop))
        self.info = _ElidedLabel()                      # values at the frame: may be long
        self.save_plot_button = QPushButton("Save plot...")
        self.save_plot_button.clicked.connect(lambda: self._guard(self.save_plot))
        self.save_data_button = QPushButton("Save data...")
        self.save_data_button.clicked.connect(lambda: self._guard(self.save_data))
        row.addWidget(self.play_button)
        row.addWidget(self.stop_button)
        row.addWidget(self.info, stretch=1)
        row.addWidget(self.save_plot_button)
        row.addWidget(self.save_data_button)
        layout.addLayout(row)
        self._show_boxes()

    @contextmanager
    def _quiet(self):
        """Set the controls without their signals (which would animate again)."""
        boxes = (self.model_box, self.scope_box, self.mode_box, self.amplitude_box,
                 self.frames_box, self.kind_box)
        was = [b.blockSignals(True) for b in boxes]
        try:
            yield
        finally:
            for b, w in zip(boxes, was, strict=True):
                b.blockSignals(w)

    def _guard(self, func, *args, **kw):
        """Run a button's action; a UserError goes to the window (and the Log as a warning)
        instead of a dialog, and the controls are set back to the plot. A notice (grey) is
        shown until the next action."""
        shown, self._notice = bool(self._error or self._notice), ""
        try:
            func(*args, **kw)
        except UserError as exc:
            self._error = str(exc)
            self.session.logger.warning(self._error)
            self.refresh()
        else:
            self._error = ""
            if shown or self._notice:
                self._update_note()

    def _show_boxes(self):
        """Show as many atom boxes as the measurement needs."""
        for k, box in enumerate(self.atom_boxes):
            box.setVisible(k < MEASURES[self.kind])

    # ------------------------------------------------------------------ state
    def set_plot(self, plot):
        """Show ``plot`` (a plot.Plot, or None for an empty window)."""
        if self.plot is not None and self._plot_changed in self.plot.listeners:
            self.plot.listeners.remove(self._plot_changed)
        self.plot = plot
        if plot is not None:
            plot.listeners.append(self._plot_changed)
            if plot not in plots(self.session):
                plots(self.session).append(plot)
            if plot.series:                             # Measure: the newest series' kind
                self.kind = plot.series[-1].kind
        self.refresh()

    def _plot_changed(self, plot):
        self._notice = ""                               # the plot changed: the grey note is old
        self.refresh()

    def _alive(self):
        return self.plot is not None and self.plot.alive

    def refresh(self):
        """Controls, series list, title and plot from the plot's current state."""
        plot = self.plot
        self._watch(plot.structure if self._alive() else None)
        with self._quiet():
            self.model_box.clear()
            self.scope_box.clear()
            if plot is None:
                mode, scope, amplitude, frames = (self.settings[k] for k in
                                                  ("mode", "scope", "amplitude", "frames"))
                self.scope_box.addItems(["cluster", "global"])
                n_modes = max(mode, 20)
            else:
                mode, scope, amplitude, frames = (plot.mode, plot.scope, plot.amplitude,
                                                  plot.n_frames)
                for m in plot.data.models:
                    self.model_box.addItem(f"{m}: {plot.data.model_label(m)}", m)
                    self.model_box.setItemData(self.model_box.count() - 1,
                                               plot.data.model_label(m),
                                               Qt.ItemDataRole.ToolTipRole)
                self.model_box.setCurrentIndex(max(plot.data.models.index(plot.model), 0))
                self.scope_box.addItems(["cluster"] + (["global"] if plot.data.global_set
                                                       else []))
                n_modes = self._n_modes(plot.model, plot.scope)
            self.scope_box.setCurrentIndex(max(self.scope_box.findText(scope), 0))
            self.mode_box.setRange(1, max(n_modes, 1))
            self.mode_box.setValue(mode)
            self.amplitude_box.setValue(amplitude)
            self.frames_box.setValue(frames)
            self.kind_box.setCurrentIndex(KINDS.index(self.kind))
        self._show_boxes()
        self.model_box.setEnabled(plot is not None)
        self.series_list.clear()
        for s in plot.series if plot is not None else []:      # row k = plot.series[k]
            pixmap = QPixmap(12, 12)
            pixmap.fill(QColor(s.color))
            self.series_list.addItem(QListWidgetItem(QIcon(pixmap),
                                                     f"{s.kind.capitalize()}: {s.label}"))
        has_series = plot is not None and bool(plot.series)
        for button in (self.remove_button, self.clear_button):
            button.setEnabled(has_series)
        for button in (self.play_button, self.stop_button):
            button.setEnabled(self._alive())
        for button in (self.save_plot_button, self.save_data_button):
            button.setEnabled(has_series)
        kinds = plot.kinds if plot is not None else []
        if plot is None:
            self.tool_window.title = TOOL_NAME
        else:
            what = f"{', '.join(kinds)} on " if kinds else ""
            self.tool_window.title = f"{TOOL_NAME}: {what}{plot.spec}"
        self.canvas.setMinimumHeight(CANVAS_HEIGHT + PANEL_PIXELS * max(len(kinds) - 1, 0))
        self._update_note()
        self._redraw()

    def _n_modes(self, model, scope):
        try:
            return self.plot.data.mode_set(model, scope == "global").n_modes
        except ValueError:
            return 0

    def _update_note(self):
        text, color = self._error, ERROR_COLOR
        if not text:
            text, color = self._notice, NOTE_COLOR
        if not text and self.plot is not None and not self.plot.alive:
            text = (f"{self.plot.spec}, the animation of this plot, was closed: Save plot and "
                    "Save data still work; change a setting or add atoms to animate again.")
        self.note.setText(text)
        self.note.setStyleSheet(f"color: {color}")
        self.note.setVisible(bool(text))

    def _shown_structure(self):
        """The first shown structure with mdcif data (a model or an animate copy), or None."""
        from chimerax.atomic import AtomicStructure

        return next((s for s in self.session.models.list(type=AtomicStructure)
                     if getattr(s, "mdcif_data", None) is not None and s.visible), None)

    def _show_atoms(self, structure=None):
        """The commands that make atoms pickable (a click on the cartoon selects whole
        residues), for ``structure``, else the plot's copy, else a shown mdcif structure."""
        if structure is None and self._alive():
            structure = self.plot.structure
        if structure is None:
            structure = self._shown_structure()
        spec = f"#{structure.id_string}" if structure is not None else "#1"
        return f"cartoon hide {spec}; show {spec}@CA atoms"

    def _redraw(self):
        """Draw the whole plot again (new series, settings, ...); the marks follow by blits."""
        frame = self.plot.structure.active_coordset_id if self._alive() else None
        n = MEASURES[self.kind]
        hint = (f"Type {n} atoms (e.g. {EXAMPLES[self.kind]}) and click Add, or\n"
                f"show the atoms ({self._show_atoms()}),\n"
                f"select {n} in order (Ctrl-click, then Ctrl-Shift-click)\n"
                "and click Add selected atoms")
        self._background = None                         # until the next full draw
        self._marker = draw(self.figure, self.plot, current=frame, hint=hint)
        self._columns = legend_layout(self.figure, self.plot) if self.plot is not None else {}
        self.canvas.draw_idle()
        self._show_frame(frame)

    def _drawn(self, _event):
        """After every full draw of the canvas (a new plot, a resize, ...): keep the picture
        without the frame marks, then draw the marks over it. The marks stay inside their
        panels (clipped to the axes), so blits repaint only the panels' boxes, a pixel wider
        for rounding: about a third of the canvas is legends, titles, labels and margins."""
        from matplotlib.transforms import Bbox

        self._background = self.canvas.copy_from_bbox(self.figure.bbox)
        self._boxes = [Bbox.from_extents(b.x0 - 1, b.y0 - 1, b.x1 + 1, b.y1 + 1)
                       for b in (ax.bbox for ax in self.figure.axes)]
        self._draw_marks()

    def _resized(self, _event):
        """The figure has the new size; a full draw follows. Legends that no longer fit
        across their panels (or could have more columns) are made again first."""
        self._background = None
        if self.plot is not None and legend_layout(self.figure, self.plot) != self._columns:
            self._redraw()

    def _draw_marks(self):
        if self._marker is not None:
            for artist in self._marker.artists:
                self.figure.draw_artist(artist)

    def _blit_marks(self):
        """Put moved marks on screen: the saved background, the marks over it, each panel
        repainted (about 1 ms, against tens of ms for a full draw). Without a background yet
        the pending full draw shows them."""
        if self._background is None:
            self.canvas.draw_idle()
            return
        self.canvas.restore_region(self._background)
        self._draw_marks()
        for box in self._boxes:
            self.canvas.blit(box)

    def _show_frame(self, frame):
        """Move the current-frame marks to ``frame`` and print the values there, by kind."""
        self.shown_frame = frame
        plot = self.plot
        if plot is None:
            self.info.setText("")
            return
        if frame is None or not self._alive():
            self.info.setText(f"{len(plot.series)} series, {plot.n_frames} frames")
            return
        if self._marker is not None:
            self._marker.set(frame)
            self._blit_marks()
        parts = []
        for kind in plot.kinds:
            digits = 2 if kind == "distance" else 1
            values = ", ".join(f"{s.values[frame - 1]:.{digits}f}" for s in plot.series_of(kind)
                               if 1 <= frame <= len(s.values))
            if values:
                parts.append(values + INFO_UNITS[kind])
        self.info.setText(f"frame {frame} / {plot.n_frames}"
                          + (f": {' · '.join(parts)}" if parts else ""))

    # ------------------------------------------------------------------ ChimeraX triggers
    def _watch(self, structure):
        """Follow the current frame of ``structure`` (None: stop following)."""
        if structure is self._watched:
            return
        if self._watch_handler is not None:
            self._watch_handler.remove()
        self._watched, self._watch_handler = structure, None
        if structure is not None:
            self._watch_handler = structure.triggers.add_handler("changes",
                                                                 self._structure_changed)

    def _structure_changed(self, _trigger, data):
        """A new frame on screen: update the marks once Qt runs, once for any number of
        changes before then."""
        structure, changes = data
        if ("active_coordset changed" in changes.structure_reasons()
                and structure is self._watched and not self._frame_due):
            self._frame_due = True
            QTimer.singleShot(0, self._frame_changed)

    def _frame_changed(self):
        self._frame_due = False
        if self._watched is not None and self._alive():      # not closed in the meantime
            self._show_frame(self._watched.active_coordset_id)

    def _models_deleted(self, _trigger, _models):
        if self._watched is not None and getattr(self._watched, "deleted", False):
            self._error = ""                            # the copy was closed: data kept
            self.refresh()

    # ------------------------------------------------------------------ actions
    def _spec(self):
        if not self._alive():
            raise UserError(f"{self.plot.spec if self.plot else 'this plot'}: no animation is "
                            "open (change a setting or add atoms)")
        return f"#{self.plot.structure.id_string}"

    def go_to_frame(self, frame):
        """What a click on the plot does: stop the motion and show frame ``frame``."""
        spec = self._spec()
        frame = min(max(int(frame), 1), self.plot.n_frames)
        run(self.session, f"coordset stop {spec}", log=False)    # 'coordset #N k' alone
        run(self.session, f"coordset {spec} {frame}", log=False)  # doesn't stop a loop
        self._show_frame(frame)

    def play(self):
        spec = self._spec()
        run(self.session, f"coordset stop {spec}", log=False)
        run(self.session, f"coordset {spec} 1,{self.plot.n_frames} loop 1000", log=False)

    def stop(self):
        run(self.session, f"coordset stop {self._spec()}", log=False)

    def _mouse_event(self, event):
        """A click or left-button drag on any panel (they share the frame axis)."""
        from matplotlib.backend_bases import MouseButton

        if event.name == "button_press_event":
            if event.button != MouseButton.LEFT:
                return
        elif MouseButton.LEFT not in event.buttons:     # a drag with the left button
            return
        if event.inaxes is None or event.xdata is None or not self._alive():
            return
        frame = min(max(round(event.xdata), 1), self.plot.n_frames)
        if event.name == "motion_notify_event" and frame == self.shown_frame:
            return
        self._guard(self.go_to_frame, frame)

    def canvas_position(self, frame, series=0):
        """Where frame ``frame`` of ``plot.series[series]`` is drawn (in its kind's panel), in
        canvas widget pixels from the top left (Qt logical pixels): where a mouse click on
        that point lands."""
        self.canvas.draw()                              # lay the axes out
        plot = self.plot
        if plot is not None and plot.series:
            s = plot.series[series]
            ax = self.figure.axes[plot.kinds.index(s.kind)]     # panels top to bottom
            y = s.values[frame - 1]
        else:
            ax = self.figure.axes[0]
            y = sum(ax.get_ylim()) / 2
        x_px, y_px = ax.transData.transform((frame, y))
        ratio = self.canvas.device_pixel_ratio
        return x_px / ratio, (self.figure.bbox.height - y_px) / ratio

    def set_kind(self, kind):
        """What Add and Add selected atoms measure (the Measure chooser): shows that many atom
        boxes. Each kind has its own panel, so the plot keeps the series of other kinds."""
        if kind not in MEASURES:
            raise UserError(f"unknown measurement {kind!r}")
        self.kind = kind
        with self._quiet():
            self.kind_box.setCurrentIndex(KINDS.index(kind))
        self._show_boxes()
        if self.plot is None or not self.plot.series:
            self._redraw()                              # the hint names the number of atoms

    def _control_values(self):
        return {"model": self.model_box.currentData(), "scope": self.scope_box.currentText(),
                "mode": self.mode_box.value(), "amplitude": self.amplitude_box.value(),
                "frames": self.frames_box.value()}

    def _setting_changed(self, *_):
        values = self._control_values()
        if self.plot is None:
            values.pop("model")
            self.settings.update(values)
            return
        self._guard(self.set_settings, **values)

    def set_settings(self, model=None, scope=None, mode=None, amplitude=None, frames=None):
        """Animate again with other settings (None = keep; what the controls do) and measure
        every series anew."""
        plot = self.plot
        if plot is None:
            raise UserError("the plot is empty: add atoms first")
        model = plot.model if model is None else model
        scope = plot.scope if scope is None else scope
        n_modes = self._n_modes(model, scope)
        if n_modes == 0:
            raise UserError(f"no {scope} modes for model {model}")
        mode = min(plot.mode if mode is None else mode, n_modes)
        amplitude = plot.amplitude if amplitude is None else amplitude
        frames = plot.n_frames if frames is None else frames
        if plot.alive and (model, scope, mode, amplitude, frames) == (
                plot.model, plot.scope, plot.mode, plot.amplitude, plot.n_frames):
            return
        cmd.reanimate(self.session, plot, model=model, mode=mode, scope=scope,
                      amplitude=amplitude, frames=frames)

    def add_selected(self):
        """Add selected atoms: the selection, in picking order."""
        from chimerax.atomic import selected_atoms

        self.add_atoms(selected_atoms(self.session))

    def add_atoms(self, atoms):
        """Add a series from ``atoms`` (as many as Measure needs, in order, of one
        structure)."""
        structures = list(atoms.unique_structures)
        if len(structures) > 1:
            raise UserError("the selected atoms must be in one structure")
        n = MEASURES[self.kind]
        if len(atoms) != n:
            tip = ""
            if len(atoms) > n and len(atoms.unique_residues.atoms) == len(atoms):
                tip = (" - a click on the cartoon selects whole residues: show the atoms "
                       f"first ({self._show_atoms(structures[0])})")
            raise UserError(f"{A_KIND[self.kind]} needs {n} selected atoms, {len(atoms)} "
                            f"{'is' if len(atoms) == 1 else 'are'} selected (used in picking "
                            f"order){tip}")
        return self._add_keys(structures[0], cmd._atom_keys(atoms))

    def add_typed(self):
        """Add (or Return in an atom box): the atoms typed in the boxes Measure shows."""
        n = MEASURES[self.kind]
        texts = [box.text() for box in self.atom_boxes[:n]]
        if not all(t.strip() for t in texts):
            raise UserError(f"{A_KIND[self.kind]} needs {n} atoms: type one in each box "
                            f"(e.g. {EXAMPLES[self.kind]})")
        return self.add_text(*texts)

    def add_text(self, *texts):
        """Add a series from typed atoms, residue number and atom name (``add_text("10CA",
        "40CA")``; 2, 3 or 4 of them: a distance, angle or dihedral), looked up where
        ``_target_structure`` says."""
        if len(texts) not in KIND_OF:
            raise UserError(f"type 2, 3 or 4 atoms (a distance, angle or dihedral), not "
                            f"{len(texts)}")
        structure = self._target_structure()
        names, keys = cmd._keys(structure), []
        for k, text in enumerate(texts, start=1):
            try:
                parsed = parse_atom(text)
            except ValueError as exc:
                raise UserError(f"atom {k}: {exc}") from None
            try:
                keys.append(match_key(parsed, names))
            except ValueError as exc:
                raise UserError(f"atom {k}, on #{structure.id_string}: {exc}"
                                f"{digit_hint(text, parsed, names)}") from None
        return self._add_keys(structure, keys)

    def _target_structure(self):
        """Where typed atoms are looked up: the plot's animate copy; with that closed, the
        model it animated (animated again when the atoms are added); in an empty window, the
        first shown mdcif structure."""
        plot = self.plot
        if self._alive():
            return plot.structure
        if plot is not None:
            found = (plot.source if cmd._alive(plot.source)
                     else cmd._model_structure(self.session, plot.path, plot.model))
            if found is None:
                raise UserError(f"no structure of model {plot.model} of "
                                f"{os.path.basename(plot.path)} is open")
            return found
        found = self._shown_structure()
        if found is None:
            raise UserError("no mdcif model is shown: open an mdcif file (or show its model) "
                            "first")
        return found

    def _add_keys(self, structure, keys):
        """Add the atoms ``keys`` of ``structure`` (a model or an animate copy) to the plot,
        starting one in an empty window or animating again a plot whose copy was closed."""
        data, _ = cmd._data(self.session, structure)
        plot = self.plot
        if plot is None:
            if getattr(structure, "mdcif_copy", False):
                copy, source = cmd.copy_for(self.session, structure)
            else:
                s = self.settings
                copy, source = cmd.copy_for(self.session, structure, s["mode"], s["scope"],
                                            s["amplitude"], s["frames"])
            plot = Plot(copy, source)
        elif attach.path_key(data.path) != attach.path_key(plot.path):
            raise UserError(f"the atoms must be on a model of {os.path.basename(plot.path)} "
                            "or one of its animations")
        elif not plot.alive:
            cmd.reanimate(self.session, plot)
        n = len(plot.series)
        series = cmd.add_series(self.session, plot, keys)
        if len(plot.series) == n:                       # also in the Log
            self._notice = f"{plot.spec} {series.kind} {series.label} is already plotted"
            self._update_note()                         # shown now, also from scripts
        if self.plot is not plot:
            self.set_plot(plot)
        return series

    def remove_selected(self):
        rows = sorted({self.series_list.row(i) for i in self.series_list.selectedItems()})
        if not rows:
            raise UserError("choose the series to remove in the list")
        for k in reversed(rows):
            del self.plot.series[k]
        self.plot.changed()

    def remove_all(self):
        if self.plot is not None:
            self.plot.series.clear()
            self.plot.changed()

    def _ask(self, title, filters, suffix):
        """A file name from a save dialog (Qt's own, as ChimeraX's), or None."""
        plot = self.plot
        base = os.path.splitext(os.path.basename(plot.path))[0]
        default = (f"{base}_model{plot.model}_{plot.scope}_mode{plot.mode}_"
                   f"{'-'.join(plot.kinds)}{suffix}")
        dialog = QFileDialog(self.tool_window.ui_area, title)
        dialog.setAcceptMode(QFileDialog.AcceptMode.AcceptSave)
        dialog.setFileMode(QFileDialog.FileMode.AnyFile)
        dialog.setOption(QFileDialog.Option.DontUseNativeDialog)
        dialog.setNameFilters(filters)
        dialog.selectFile(os.path.join(os.getcwd(), default))
        if not dialog.exec() or not dialog.selectedFiles():
            return None
        path = dialog.selectedFiles()[0]
        if not os.path.splitext(path)[1]:               # suffix from the chosen file type
            path += dialog.selectedNameFilter().split("(*")[-1].rstrip(")")
        return path

    def save_plot(self, path=None):
        """The plot as PNG, SVG or PDF (by suffix; every panel, 6 in wide, 3.5 in high for one
        panel and 2.2 in more for each further one, 300 dpi, no frame marks). Without
        ``path`` a save dialog asks. Returns the path or None."""
        if self.plot is None or not self.plot.series:
            raise UserError("nothing to save: the plot is empty")
        path = path or self._ask("Save Plot", IMAGE_FILTERS, ".png")
        if path:
            cmd._save(self.session, self.plot, path, "image")
        return path

    def save_data(self, path=None):
        """The values as CSV ('#' header lines, then frame + one column per series)."""
        if self.plot is None or not self.plot.series:
            raise UserError("nothing to save: the plot is empty")
        path = path or self._ask("Save Data", DATA_FILTERS, ".csv")
        if path:
            cmd._save(self.session, self.plot, path, "data")
        return path

    # ------------------------------------------------------------------ cleanup
    def delete(self):
        """Closing the window: its handlers go, and its plot with it (the copy stays)."""
        self._watch(None)                               # also stops a pending frame update
        self.session.triggers.remove_handler(self._delete_handler)
        for cid in self._mouse:
            self.canvas.mpl_disconnect(cid)
        plot, self.plot = self.plot, None
        if plot is not None:
            if self._plot_changed in plot.listeners:
                plot.listeners.remove(self._plot_changed)
            found = plots(self.session)
            if not plot.listeners and plot in found:
                found.remove(plot)
        super().delete()
