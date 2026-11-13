"""The Geometry tab of the MDCIF Analyze window: distances, angles and dihedrals along a mode
animation (GUI only; ported from the ChimeraX bundle).

The window around it (``analyze_tool.MdcifAnalyze``, opened from MDCIF > Analyze..., or by
'mdcif_plot' with the GUI) hands it a page to build into and gives it its title. The plot itself
(state, drawing, files) is in plot.py and the PyMOL work in commands.py, so the tab and the
command share their plots:

    Model [1: cluster 1] Scope [cluster] Mode [1] Amplitude [3.00 sd] Frames [30]
    Measure [Distance] Atoms [10CA] [40CA] [Add] [Pick atoms] [Add picked]  [Remove] [Remove all]
    (series: colour, kind and atoms of each)
    (the plot: a panel per kind - distances, then angles, then dihedrals - on one frame axis,
     points at every frame, the current frame marked in each; click or drag = stop there)
    [Play] [Stop]  frame 12 / 30: 8.43, 10.12 Å · 141.3°     [Save plot...] [Save data...]

Measure sets what Add, Pick atoms and Add picked measure, and shows that many atom boxes (2 for
a distance, 3 for an angle, 4 for a dihedral); each kind gets its own panel, so a plot can show
all three. Atoms are typed as residue number and atom name (38CA, 38 CA, 38@CA; A:38CA when the
number is in several chains) and looked up on the plot's animation (on its model if the
animation was deleted, which is then animated again), or in an empty tab on its target (the
animation of the Animate window's Plot geometry... link) or else the first shown mdcif object;
the boxes keep their text. Pick atoms starts a small PyMOL wizard: click the atoms in the viewer
in order (a click on the cartoon picks the residue's CA), and the series is added when there are
enough. Add picked takes PyMOL's pk1, pk2, ... (Ctrl-middle-click), in that order.

Changing Model, Scope, Mode, Amplitude or Frames animates again (``commands.reanimate``) and
measures every series anew. An empty tab takes the first atoms' animation: atoms on an animate
copy use that copy (the controls then show its settings), atoms on a model are animated with the
controls' settings. Play plays from frame 1 ('frame 1; mplay'); a click on the plot is 'mstop'
and 'frame k'. Analyze windows are not saved in sessions.

While the motion plays only the frame marks are redrawn: they are animated matplotlib artists,
drawn over a copy of the rest of the canvas taken after each full draw (blitting). Without
matplotlib in PyMOL's Python the tab says so (the plugin installs it when it loads, deps.py:
a tab opened after that has its plot); the values, the list and Save data still work.

For scripts: ``analyze_tool.windows()`` lists the Geometry tabs of the open Analyze windows,
oldest first. Each has ``host`` (its window), ``plot`` (a ``plot.Plot``), ``figure`` /
``canvas`` (matplotlib; ``figure.axes`` = the panels, top to bottom), ``go_to_frame(k)`` (what a
click does), ``play()``, ``stop()``, ``add_text(*texts)`` (``add_text("10CA", "40CA")``),
``add_selection(*selections)``, ``add_picked()``, ``pick_atoms()``, ``set_settings(model=,
scope=, mode=, amplitude=, frames=)``, ``set_kind(kind)``, ``save_plot(path)``,
``save_data(path)``, ``shown_frame`` and ``canvas_position(frame, series=0)``.
"""

import os

from pymol import CmdException, cmd
from pymol.Qt import QtCore, QtGui, QtWidgets
from pymol.wizard import Wizard

from . import attach, deps
from . import commands as cx
from .mdcif_data import MEASURES
from .menu import ANALYZE
from .plot import (
    A_KIND,
    IMAGE_SIZE,
    KIND_OF,
    Plot,
    atom_label,
    digit_hint,
    draw,
    legend_layout,
    match_key,
    parse_atom,
    plots,
)
from .widgets import (
    ElidedLabel,
    FrameWatcher,
    Guarded,
    Note,
    ask_save_path,
    guarded_call,
    quiet,
    shown_frame,
)

Qt = QtCore.Qt
KINDS = list(MEASURES)
EXAMPLES = {"distance": "10CA 40CA", "angle": "10CA 11CA 12CA",
            "dihedral": "10CA 11CA 12CA 13CA"}
INFO_UNITS = {"distance": " Å", "angle": "°", "dihedral": "°"}       # 8.43 Å, 141.3°
ATOM_TIP = ("Residue number and atom name: 38CA, 38 CA or 38@CA (38@1HB when the name starts "
            "with a digit); A:38CA for chain A when the number is in several chains. Return "
            "adds.")
IMAGE_FILTERS = ["PNG image (*.png)", "SVG image (*.svg)", "PDF document (*.pdf)"]
DATA_FILTERS = ["CSV file (*.csv)"]
CANVAS_HEIGHT, PANEL_PIXELS = 200, 150          # minimum canvas height, + per further panel
PICKED = "mdcif_picked"                 # the selection showing the atoms picked so far


def matplotlib_canvas():
    """matplotlib's Figure and Qt canvas classes, or (None, None) when PyMOL's Python lacks
    matplotlib. Looked up for each new tab: the plugin may have installed it since (deps.py)."""
    try:
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
        from matplotlib.figure import Figure
    except ImportError:
        return None, None
    return Figure, FigureCanvasQTAgg


class PickAtoms(Wizard):
    """A PyMOL wizard that collects ``n`` atoms clicked in the viewer, in order (as the
    Measurement wizard: the mouse selects atoms, each click becomes a pick), then calls
    ``done([(object, key), ...])`` once Qt runs. Cancel or another wizard ends it."""

    def __init__(self, n, kind, done, _self=cmd):
        super().__init__(_self)
        self.n, self.kind, self.done = n, kind, done
        self.picked = []
        self.selection_mode = self.cmd.get_setting_int("mouse_selection_mode")
        self.cmd.set("mouse_selection_mode", 0)        # atoms
        self.cmd.deselect()

    def get_event_mask(self):
        return Wizard.event_mask_pick + Wizard.event_mask_select

    def get_prompt(self):
        got = ", ".join(atom_label(k) for _, k in self.picked)
        return [f"MDCIF: click atom {len(self.picked) + 1} of {self.n} for {A_KIND[self.kind]}"
                + (f" (so far {got})" if got else "")]

    def get_panel(self):
        return [[1, "MDCIF: pick atoms", ""],
                [2, "Undo last atom", "cmd.get_wizard().undo()"],
                [2, "Cancel", "cmd.set_wizard()"]]

    def undo(self):
        if self.picked:
            self.picked.pop()
            self._show()
        self.cmd.refresh_wizard()

    def do_select(self, name):                          # map selects into picks
        self.cmd.unpick()
        try:
            self.cmd.select("pk1", name)
            self.cmd.delete(name)
            self.do_pick(0)
        except CmdException:
            pass

    def do_pick(self, bondFlag):
        if bondFlag:
            print(" mdcif: please pick an atom, not a bond")
            return
        found = cx._atoms_of("pk1")
        self.cmd.unpick()
        if len(found) != 1:
            return
        if found[0] in self.picked:
            print(f" mdcif: {atom_label(found[0][1])} is picked already")
            return
        self.picked.append(found[0])
        self._show()
        if len(self.picked) >= self.n:
            picked = list(self.picked)
            QtCore.QTimer.singleShot(0, lambda: self._finish(picked))
        self.cmd.refresh_wizard()

    def _finish(self, picked):
        if self.cmd.get_wizard() is self:
            self.cmd.set_wizard()
        guarded_call(self.done, picked)

    def _show(self):
        self.cmd.delete(PICKED)
        if self.picked:
            sele = " or ".join(f"/{o}//{c}/{r}/{a}" for o, (c, r, a) in self.picked)
            self.cmd.select(PICKED, sele, enable=1)

    def cleanup(self):
        self.cmd.set("mouse_selection_mode", self.selection_mode)
        self.cmd.delete(PICKED)
        super().cleanup()


class GeometryTab(Guarded):
    """The Geometry tab, built into ``page`` (a QWidget of the ``host`` window, which has
    ``set_title(text)``, ``display(show)`` and ``delete()``). The host calls ``set_plot`` once
    the window is built, and ``cleanup`` when it closes."""

    def __init__(self, host, page):
        self.host = host
        self.plot = None
        self.kind = KINDS[0]                            # what Add / Pick atoms measure
        # settings for the first animation of an empty tab: mode, scope, amplitude, frames
        self.settings = dict(zip(("mode", "scope", "amplitude", "frames"),
                                 cx.ANIMATE_DEFAULTS, strict=True))
        self.target = self._target_source = None        # where an empty tab takes typed atoms
        self.shown_frame = None                         # frame the marks are on
        self._marker = None
        self._watcher = FrameWatcher(self._frame_changed)
        self._background = None                         # the canvas without the frame marks
        self._columns = {}                              # legend columns of the last full draw
        self._boxes = []                                # the panels' boxes: what blits repaint
        self._mouse = []
        self._wizard = None
        self._error = self._notice = ""
        self.page = page
        self._build(page)
        self._poll_timer = QtCore.QTimer(page)
        self._poll_timer.setInterval(500)
        self._poll_timer.timeout.connect(lambda: guarded_call(self._poll))
        self._poll_timer.start()
        self._watched_alive = None

    def display(self, show):
        self.host.display(show)

    def delete(self):
        self.host.delete()

    # ------------------------------------------------------------------ layout
    def _build(self, area):
        layout = QtWidgets.QVBoxLayout(area)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(4)

        row = QtWidgets.QHBoxLayout()
        self.model_box = QtWidgets.QComboBox()          # cluster names are long: cut them short
        self.model_box.setSizeAdjustPolicy(
            QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.model_box.setMinimumContentsLength(14)
        self.scope_box = QtWidgets.QComboBox()
        self.mode_box = QtWidgets.QSpinBox()
        self.amplitude_box = QtWidgets.QDoubleSpinBox()
        self.amplitude_box.setRange(-20.0, 20.0)
        self.amplitude_box.setSingleStep(0.5)
        self.amplitude_box.setDecimals(2)
        self.amplitude_box.setSuffix(" sd")
        self.frames_box = QtWidgets.QSpinBox()
        self.frames_box.setRange(cx.MIN_ANIMATE_FRAMES, 1000)
        for label, widget in (("Model", self.model_box), ("Scope", self.scope_box),
                              ("Mode", self.mode_box), ("Amplitude", self.amplitude_box)):
            row.addWidget(QtWidgets.QLabel(label))
            row.addWidget(widget)
        row.addSpacing(6)
        row.addWidget(QtWidgets.QLabel("Frames"))
        row.addWidget(self.frames_box)
        row.addStretch(1)
        layout.addLayout(row)
        for box in (self.mode_box, self.amplitude_box, self.frames_box):
            box.setKeyboardTracking(False)              # typed values count on Enter / leaving
            box.valueChanged.connect(lambda *_: guarded_call(self._setting_changed))
        for box in (self.model_box, self.scope_box):
            box.currentIndexChanged.connect(lambda *_: guarded_call(self._setting_changed))

        row = QtWidgets.QHBoxLayout()
        self.kind_box = QtWidgets.QComboBox()
        self.kind_box.addItems([k.capitalize() for k in KINDS])
        self.kind_box.setToolTip("What Add, Pick atoms and Add picked measure; each kind gets "
                                 "its own panel")
        self.kind_box.currentIndexChanged.connect(
            lambda i: self._guard(self.set_kind, KINDS[i]))
        self.atom_boxes = []                            # as many shown as Measure needs
        for _ in range(max(MEASURES.values())):
            box = QtWidgets.QLineEdit()
            box.setPlaceholderText("e.g. 38CA")
            box.setToolTip(ATOM_TIP)
            box.setFixedWidth(box.fontMetrics().horizontalAdvance("e.g. 38CA") + 16)
            box.returnPressed.connect(lambda: self._guard(self.add_typed))
            self.atom_boxes.append(box)
        self.add_typed_button = QtWidgets.QPushButton("Add")
        self.add_typed_button.setToolTip("Adds the atoms typed in the boxes")
        self.add_typed_button.clicked.connect(lambda: self._guard(self.add_typed))
        self.pick_button = QtWidgets.QPushButton("Pick atoms")
        self.pick_button.setToolTip("Then click the atoms in the viewer, in order (a click on "
                                    "the cartoon picks the residue's CA); added when there are "
                                    "enough")
        self.pick_button.clicked.connect(lambda: self._guard(self.pick_atoms))
        self.add_button = QtWidgets.QPushButton("Add picked")
        self.add_button.setToolTip("Adds PyMOL's picked atoms pk1, pk2, ... in that order "
                                   "(Ctrl-middle-click picks atoms)")
        self.add_button.clicked.connect(lambda: self._guard(self.add_picked))
        self.remove_button = QtWidgets.QPushButton("Remove")
        self.remove_button.clicked.connect(lambda: self._guard(self.remove_selected))
        self.clear_button = QtWidgets.QPushButton("Remove all")
        self.clear_button.clicked.connect(lambda: self._guard(self.remove_all))
        row.addWidget(QtWidgets.QLabel("Measure"))
        row.addWidget(self.kind_box)
        row.addSpacing(6)
        row.addWidget(QtWidgets.QLabel("Atoms"))
        for box in self.atom_boxes:
            row.addWidget(box)
        row.addWidget(self.add_typed_button)
        row.addWidget(self.pick_button)
        row.addWidget(self.add_button)
        row.addStretch(1)
        row.addWidget(self.remove_button)
        row.addWidget(self.clear_button)
        layout.addLayout(row)

        self.series_list = QtWidgets.QListWidget()
        self.series_list.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection)
        self.series_list.setMaximumHeight(4 * self.series_list.fontMetrics().height() + 12)
        layout.addWidget(self.series_list)
        self.note = Note()
        layout.addWidget(self.note)

        Figure, FigureCanvasQTAgg = matplotlib_canvas()
        if Figure is not None:
            self.figure = Figure(figsize=IMAGE_SIZE, dpi=100)
            self.canvas = FigureCanvasQTAgg(self.figure)
            self.canvas.setMinimumSize(320, CANVAS_HEIGHT)
            layout.addWidget(self.canvas, stretch=1)
            self._mouse = [self.canvas.mpl_connect("button_press_event", self._mouse_event),
                           self.canvas.mpl_connect("motion_notify_event", self._mouse_event),
                           self.canvas.mpl_connect("draw_event", self._drawn),
                           self.canvas.mpl_connect("resize_event", self._resized)]
        else:
            self.figure = self.canvas = None
            missing = QtWidgets.QLabel(deps.missing_text())
            missing.setWordWrap(True)
            missing.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
            missing.setStyleSheet("color: #d03b3b")
            layout.addWidget(missing, stretch=1)

        row = QtWidgets.QHBoxLayout()
        self.play_button = QtWidgets.QPushButton("Play")
        self.play_button.setToolTip("Play the animation in a loop, from frame 1")
        self.play_button.clicked.connect(lambda: self._guard(self.play))
        self.stop_button = QtWidgets.QPushButton("Stop")
        self.stop_button.setToolTip("Stop the animation (a click on the plot also stops it, "
                                    "at that frame)")
        self.stop_button.clicked.connect(lambda: self._guard(self.stop))
        self.info = ElidedLabel()                       # values at the frame: may be long
        self.save_plot_button = QtWidgets.QPushButton("Save plot...")
        self.save_plot_button.clicked.connect(lambda: self._guard(self.save_plot))
        self.save_data_button = QtWidgets.QPushButton("Save data...")
        self.save_data_button.clicked.connect(lambda: self._guard(self.save_data))
        row.addWidget(self.play_button)
        row.addWidget(self.stop_button)
        row.addWidget(self.info, stretch=1)
        row.addWidget(self.save_plot_button)
        row.addWidget(self.save_data_button)
        layout.addLayout(row)
        self._show_boxes()

    def _quiet(self):
        """Set the controls without their signals (which would animate again)."""
        return quiet(self.model_box, self.scope_box, self.mode_box, self.amplitude_box,
                     self.frames_box, self.kind_box)

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
            if plot not in plots():
                plots().append(plot)
            if plot.series:                             # Measure: the newest series' kind
                self.kind = plot.series[-1].kind
        self.refresh()

    def set_target(self, e, source=None):
        """For an empty tab: take typed atoms on entry ``e`` (the Animate window's copy, from
        its Plot geometry... button; ``source``: the model it animates), whose settings the
        controls then show."""
        self.target, self._target_source = e, source
        if e is not None and e.copy:
            have = {"mode": e.mode, "scope": e.scope, "amplitude": e.amplitude,
                    "frames": e.n_frames}
            self.settings.update({k: v for k, v in have.items() if v is not None})
        self.refresh()

    def _plot_changed(self, plot):
        self._notice = ""                               # the plot changed: the grey note is old
        guarded_call(self.refresh)

    def _alive(self):
        return self.plot is not None and self.plot.alive

    def refresh(self):
        """Controls, series list, title and plot from the plot's current state."""
        plot = self.plot
        self._watcher.watch(plot.entry if self._alive() else None)
        self._watched_alive = self._alive()
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
            pixmap = QtGui.QPixmap(12, 12)
            pixmap.fill(QtGui.QColor(s.color))
            self.series_list.addItem(QtWidgets.QListWidgetItem(
                QtGui.QIcon(pixmap), f"{s.kind.capitalize()}: {s.label}"))
        has_series = plot is not None and bool(plot.series)
        for button in (self.remove_button, self.clear_button):
            button.setEnabled(has_series)
        for button in (self.play_button, self.stop_button):
            button.setEnabled(self._alive())
        self.save_plot_button.setEnabled(has_series and self.figure is not None)
        self.save_data_button.setEnabled(has_series)
        kinds = plot.kinds if plot is not None else []
        if plot is None:
            self.host.set_title(ANALYZE)
        else:
            what = f"{', '.join(kinds)} on " if kinds else ""
            self.host.set_title(f"{ANALYZE}: {what}{plot.spec}")
        if self.canvas is not None:
            self.canvas.setMinimumHeight(CANVAS_HEIGHT + PANEL_PIXELS * max(len(kinds) - 1, 0))
        self._update_note()
        self._redraw()

    def _n_modes(self, model, scope):
        try:
            return self.plot.data.mode_set(model, scope == "global").n_modes
        except ValueError:
            return 0

    def _update_note(self):
        text = self._error
        if not text:
            text = self._notice
        if not text and self.plot is not None and not self.plot.alive:
            text = (f"{self.plot.spec}, the animation of this plot, was deleted: Save plot and "
                    "Save data still work; change a setting or add atoms to animate again.")
        self.note.show_text(text, error=bool(self._error))

    def _shown_entry(self):
        """The first shown object with mdcif data (a model or an animate copy), or None."""
        enabled = set(cmd.get_names("objects", enabled_only=1))
        return next((e for e in attach.entries() if e.data is not None and e.name in enabled),
                    None)

    def _hint(self):
        n = MEASURES[self.kind]
        return (f"Type {n} atoms (e.g. {EXAMPLES[self.kind]}) and click Add, or\n"
                f"click Pick atoms, then click {n} atoms in the viewer, in order\n"
                "(a click on the cartoon picks the residue's CA;\n"
                "'show sticks' to pick other atoms)")

    def _redraw(self):
        """Draw the whole plot again (new series, settings, ...); the marks follow by blits."""
        frame = shown_frame(self.plot.entry) if self._alive() else None
        if self.figure is not None:
            self._background = None                     # until the next full draw
            self._marker = draw(self.figure, self.plot, current=frame, hint=self._hint())
            self._columns = (legend_layout(self.figure, self.plot) if self.plot is not None
                             else {})
            self.canvas.draw_idle()
        self._show_frame(frame)

    def _drawn(self, _event):
        """After every full draw of the canvas: keep the picture without the frame marks, then
        draw the marks over it (blits repaint only the panels' boxes)."""
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
            guarded_call(self._redraw)

    def _draw_marks(self):
        if self._marker is not None:
            for artist in self._marker.artists:
                self.figure.draw_artist(artist)

    def _blit_marks(self):
        """Put moved marks on screen: the saved background, the marks over it, each panel
        repainted. Without a background yet the pending full draw shows them."""
        if self.canvas is None:
            return
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

    # ------------------------------------------------------------------ PyMOL changes
    def _frame_changed(self, e):
        """A new frame on screen (``FrameWatcher``): move the marks."""
        if self._alive():
            self._show_frame(shown_frame(e))

    def _poll(self):
        """The plot's copy deleted (the data stay), or an empty tab's target deleted."""
        if self.plot is not None and self._alive() != self._watched_alive:
            self._error = ""
            self.refresh()
        elif self.plot is None and self.target is not None and not self.target.alive:
            self.target = self._target_source = None
            self.refresh()

    # ------------------------------------------------------------------ actions
    def _name(self):
        if not self._alive():
            cx._error(f"{self.plot.spec if self.plot else 'this plot'}: no animation is "
                      "there (change a setting or add atoms)")
        return self.plot.entry.name

    def go_to_frame(self, frame):
        """What a click on the plot does: stop the motion and show frame ``frame``."""
        self._name()
        frame = min(max(int(frame), 1), self.plot.n_frames)
        cmd.mstop()
        cmd.frame(frame)
        self._show_frame(frame)

    def play(self):
        name = self._name()
        if name not in cmd.get_names("objects", enabled_only=1):
            cmd.enable(name)
        cx._play()

    def stop(self):
        self._name()
        cmd.mstop()

    def _mouse_event(self, event):
        """A click or left-button drag on any panel (they share the frame axis)."""
        from matplotlib.backend_bases import MouseButton

        if event.name == "button_press_event":
            if event.button != MouseButton.LEFT:
                return
        else:                                           # a drag with the left button
            buttons = getattr(event, "buttons", None)   # (matplotlib >= 3.9)
            if (event.button != MouseButton.LEFT if buttons is None
                    else MouseButton.LEFT not in buttons):
                return
        if event.inaxes is None or event.xdata is None or not self._alive():
            return
        frame = min(max(round(event.xdata), 1), self.plot.n_frames)
        if event.name == "motion_notify_event" and frame == self.shown_frame:
            return
        self._guard(self.go_to_frame, frame)

    def canvas_position(self, frame, series=0):
        """Where frame ``frame`` of ``plot.series[series]`` is drawn (in its kind's panel), in
        canvas widget pixels from the top left: where a mouse click on that point lands."""
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
        """What Add, Pick atoms and Add picked measure (the Measure chooser): shows that many
        atom boxes. Each kind has its own panel, so the plot keeps the series of other
        kinds."""
        if kind not in MEASURES:
            cx._error(f"unknown measurement {kind!r}")
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
            cx._error("the plot is empty: add atoms first")
        model = plot.model if model is None else model
        scope = plot.scope if scope is None else scope
        n_modes = self._n_modes(model, scope)
        if n_modes == 0:
            cx._error(f"no {scope} modes for model {model}")
        mode = min(plot.mode if mode is None else mode, n_modes)
        amplitude = plot.amplitude if amplitude is None else amplitude
        frames = plot.n_frames if frames is None else frames
        if plot.alive and (model, scope, mode, amplitude, frames) == (
                plot.model, plot.scope, plot.mode, plot.amplitude, plot.n_frames):
            return
        cx.reanimate(plot, model=model, mode=mode, scope=scope, amplitude=amplitude,
                     frames=frames)

    def pick_atoms(self):
        """Pick atoms: a PyMOL wizard collects as many clicked atoms as Measure needs, in
        order, then adds them."""
        n, kind = MEASURES[self.kind], self.kind
        self._wizard = PickAtoms(n, kind, lambda picked: self._guard(self._add_picked_list,
                                                                     kind, picked))
        cmd.set_wizard(self._wizard)
        self._notice = (f"Click {n} atoms in the viewer, in order (Cancel in the wizard "
                        "panel stops).")
        self._update_note()
        return self._wizard

    def _add_picked_list(self, kind, picked):
        objects = list(dict.fromkeys(o for o, _ in picked))
        if len(objects) > 1:
            cx._error(f"the atoms must be in one object, not {', '.join(objects)}")
        if kind != self.kind:
            self.set_kind(kind)
        return self._add_object_keys(objects[0], [k for _, k in picked])

    def add_picked(self):
        """Add picked: PyMOL's pk1, pk2, ... in that order, as many as Measure needs."""
        n = MEASURES[self.kind]
        names = set(cmd.get_names("selections"))
        have = [f"pk{k}" for k in range(1, n + 1) if f"pk{k}" in names
                and cmd.count_atoms(f"pk{k}") == 1]
        if len(have) < n:
            cx._error(f"{A_KIND[self.kind]} needs {n} picked atoms (pk1 to pk{n}), "
                      f"{len(have)} are picked: click Pick atoms, or Ctrl-middle-click the "
                      "atoms in order")
        return self.add_selection(*have)

    def add_selection(self, *selections):
        """Add a series from selections: one atom each, in order (or one selection of exactly
        as many atoms as Measure needs, in atom order)."""
        obj, keys = cx.plot_atoms(self.kind, selections)
        return self._add_object_keys(obj, keys)

    def _add_object_keys(self, obj, keys):
        e = attach.entry(obj)
        if e is None:
            cx._error(f"{obj} has no mdcif data (pick atoms on an mdcif model or animation)")
        return self._add_keys(e, keys)

    def add_typed(self):
        """Add (or Return in an atom box): the atoms typed in the boxes Measure shows."""
        n = MEASURES[self.kind]
        texts = [box.text() for box in self.atom_boxes[:n]]
        if not all(t.strip() for t in texts):
            cx._error(f"{A_KIND[self.kind]} needs {n} atoms: type one in each box "
                      f"(e.g. {EXAMPLES[self.kind]})")
        return self.add_text(*texts)

    def add_text(self, *texts):
        """Add a series from typed atoms, residue number and atom name (``add_text("10CA",
        "40CA")``; 2, 3 or 4 of them: a distance, angle or dihedral), looked up where
        ``_target_entry`` says."""
        if len(texts) not in KIND_OF:
            cx._error(f"type 2, 3 or 4 atoms (a distance, angle or dihedral), not "
                      f"{len(texts)}")
        e = self._target_entry()
        names, keys = cx._keys(e.name), []
        for k, text in enumerate(texts, start=1):
            try:
                parsed = parse_atom(text)
            except ValueError as exc:
                cx._error(f"atom {k}: {exc}")
            try:
                keys.append(match_key(parsed, names))
            except ValueError as exc:
                cx._error(f"atom {k}, on {e.name}: {exc}{digit_hint(text, parsed, names)}")
        return self._add_keys(e, keys)

    def _target_entry(self):
        """Where typed atoms are looked up: the plot's animate copy; with that deleted, the
        model it animated (animated again when the atoms are added); in an empty tab, its
        target (``set_target``), else the first shown mdcif object."""
        plot = self.plot
        if self._alive():
            return plot.entry
        if plot is not None:
            found = (plot.source if plot.source is not None and plot.source.alive
                     else cx.model_entry(plot.path, plot.model))
            if found is None:
                cx._error(f"no object of model {plot.model} of "
                          f"{os.path.basename(plot.path)} is loaded")
            return found
        if self.target is not None and self.target.alive:
            return self.target
        found = self._shown_entry()
        if found is None:
            cx._error("no mdcif model is shown: load an mdcif file (or show its model) first")
        return found

    def _add_keys(self, e, keys):
        """Add the atoms ``keys`` of entry ``e`` (a model or an animate copy) to the plot,
        starting one in an empty tab or animating again a plot whose copy was deleted."""
        cx._check_data(e)
        plot = self.plot
        if plot is None:
            if e.copy:
                copy, source = cx.copy_for(e)
                if source is None and copy is self.target:      # the Animate window's model
                    source = self._target_source
            else:
                s = self.settings
                copy, source = cx.copy_for(e, s["mode"], s["scope"], s["amplitude"],
                                           s["frames"])
            plot = Plot(copy, source)
        elif attach.path_key(e.path) != attach.path_key(plot.path):
            cx._error(f"the atoms must be on a model of {os.path.basename(plot.path)} or one "
                      "of its animations")
        elif not plot.alive:
            cx.reanimate(plot)
        n = len(plot.series)
        series = cx.add_series(plot, keys)
        if len(plot.series) == n:                       # also in PyMOL's output
            self._notice = f"{plot.spec} {series.kind} {series.label} is already plotted"
            self._update_note()
        if self.plot is not plot:
            self.set_plot(plot)
        return series

    def remove_selected(self):
        rows = sorted({self.series_list.row(i) for i in self.series_list.selectedItems()})
        if not rows:
            cx._error("choose the series to remove in the list")
        for k in reversed(rows):
            del self.plot.series[k]
        self.plot.changed()

    def remove_all(self):
        if self.plot is not None:
            self.plot.series.clear()
            self.plot.changed()

    def _ask(self, title, filters, suffix):
        """A file name from a save dialog, or None."""
        plot = self.plot
        base = os.path.splitext(os.path.basename(plot.path))[0]
        default = (f"{base}_model{plot.model}_{plot.scope}_mode{plot.mode}_"
                   f"{'-'.join(plot.kinds)}{suffix}")
        return ask_save_path(self.page, title, filters, default)

    def save_plot(self, path=None):
        """The plot as PNG, SVG or PDF (by suffix; every panel, 6 in wide, 3.5 in high for one
        panel and 2.2 in more for each further one, 300 dpi, no frame marks). Without
        ``path`` a save dialog asks. Returns the path or None."""
        if self.plot is None or not self.plot.series:
            cx._error("nothing to save: the plot is empty")
        path = path or self._ask("Save Plot", IMAGE_FILTERS, ".png")
        if path:
            cx.save_plot_file(self.plot, path, "image")
        return path

    def save_data(self, path=None):
        """The values as CSV ('#' header lines, then frame + one column per series)."""
        if self.plot is None or not self.plot.series:
            cx._error("nothing to save: the plot is empty")
        path = path or self._ask("Save Data", DATA_FILTERS, ".csv")
        if path:
            cx.save_plot_file(self.plot, path, "data")
        return path

    # ------------------------------------------------------------------ cleanup
    def cleanup(self):
        """The window closes: the tab's timers go, and its plot with it (the copy stays)."""
        self._watcher.watch(None)
        self._poll_timer.stop()
        if self._wizard is not None and cmd.get_wizard() is self._wizard:
            cmd.set_wizard()
        if self.canvas is not None:
            for cid in self._mouse:
                self.canvas.mpl_disconnect(cid)
        plot, self.plot = self.plot, None
        if plot is not None:
            if self._plot_changed in plot.listeners:
                plot.listeners.remove(self._plot_changed)
            found = plots()
            if not plot.listeners and plot in found:
                found.remove(plot)
