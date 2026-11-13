"""The MDCIF Porcupine window: porcupine arrows for one or more modes, with sliders (GUI only).

    ┌ MDCIF Porcupine ───────────────────────────────────────────────┐
    │ Object      [gb1_m1  Cluster 1 medoid (frame 142700, ...)   ▾] │
    │ Scope       (•) cluster  ( ) global                            │
    │ Modes       [x] 1   26% variance   sd 6.12 Å   [colour]        │
    │             [ ] 2   20% variance   sd 5.34 Å   [colour]        │
    │ Length      ├─────●──────┤ [3.0 sd]     0.5 to 10 sd           │
    │ Thickness   ├──●─────────┤ [0.15 Å]     0.05 to 0.6 Å          │
    │ Min length  ├─●──────────┤ [0.5 Å]      0 to 3 Å               │
    │ [x] Live  [Draw]  [Keep]  [Clear]           [Animate this mode] │
    │ gb1_m1_mode1_porcupine: 51 arrows (5 shorter than 0.5 Å ...)    │
    └────────────────────────────────────────────────────────────────┘

One per session (MDCIF > Porcupine..., Plugin > MDCIF, 'mdcif_window porcupine'). Ported from
the ChimeraX bundle. The arrows are those of 'mdcif_porcupine' (``commands.porcupine_cgo``):
Length is its ``scale`` (sd), Thickness its ``radius`` and Min length its ``min_length``. The
objects are the mdcif models and their 'mdcif_animate' copies (arrows on a copy start from its
state 1, the deposited model, as the command's do).

Every ticked mode is drawn, one CGO object per mode (``<object>_mode<N>_porcupine``, or with
``_2``, ``_3``... when that name is taken), in its colour button's colour: mode 1 the command's
orange, modes 2 and up the plot palette. With Live on, every slider move, tick or colour change
redraws the window's arrows in place (the same objects); a typed number counts on Enter. With
Live off, Draw applies the settings. Keep leaves the arrows as they are and starts new ones on
the next change (to compare settings side by side); Clear deletes the window's arrows that were
not kept. Another object or scope deletes the window's arrows (not the kept ones) and draws
anew. Closing the window keeps its arrows.

PyMOL's output gets one 'mdcif_porcupine ...' line per ticked mode once a change has settled
(slider released, Enter, or 0.4 s without a change), not one per slider step. When the chosen
model is hidden because its animation plays, a note offers to draw on the animation. "Animate
this mode" opens the Animate window on the object, scope and lowest ticked mode.
``preset(entry, scope, modes)`` is what the Animate window's "Arrows for this mode" calls.
"""

from pymol import cmd
from pymol.Qt import QtCore, QtWidgets

from . import commands as cx
from .menu import PORCUPINE
from .plot import COLORS
from .widgets import (
    NOTE_COLOR,
    POLL_MS,
    SETTLE_MS,
    ColorButton,
    Guarded,
    Lever,
    Note,
    ObjectBox,
    ToolWindow,
    guarded_call,
    quiet,
)

LENGTH = (0.5, 10.0, 0.1, 3.0)          # min, max, step, start: 'scale', sd
THICKNESS = (0.05, 0.6, 0.01, 0.15)     # 'radius', A
MIN_LENGTH = (0.0, 3.0, 0.1, 0.5)       # 'min_length', A
MODE_ROWS = 6                           # mode rows shown before the list scrolls


def _rgb255(hex_color):
    return tuple(int(hex_color[k:k + 2], 16) for k in (1, 3, 5))


# mode 1: the command's orange; modes 2 and up: the plot palette without its orange
MODE_COLORS = [tuple(cx.PORCUPINE_COLOR)] + [_rgb255(c) for c in COLORS if c != "#eb6834"]


def _exists(name):
    return name is not None and name in cmd.get_names("objects")


class MdcifPorcupine(Guarded, ToolWindow):

    def __init__(self):
        super().__init__(PORCUPINE, width=460)
        self.entry, self.scope = None, "cluster"
        self.ticked = {1}                               # modes drawn
        self.colors = {}                                # mode -> RGB chosen (else MODE_COLORS)
        self.live = {}                                  # mode -> CGO object the sliders move
        self.drawn = {}                                 # mode -> (arrows, too short): last draw
        self.failed = {}                                # mode -> why it has no arrows
        self._rows = {}                                 # mode -> (tick box, colour button)
        self._redraw_due = self._log_due = False
        self._watch = None                              # (enabled, live objects) last seen
        self._build(self.ui_area)
        self._log_timer = QtCore.QTimer(self)
        self._log_timer.setSingleShot(True)
        self._log_timer.setInterval(SETTLE_MS)
        self._log_timer.timeout.connect(lambda: guarded_call(self._flush_log))
        self._poll_timer = QtCore.QTimer(self)
        self._poll_timer.setInterval(POLL_MS)
        self._poll_timer.timeout.connect(lambda: guarded_call(self._poll))
        self._poll_timer.start()
        self._structure_picked(self.box.current())

    # ------------------------------------------------------------------ layout
    def _build(self, area):
        layout = QtWidgets.QVBoxLayout(area)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(4)
        grid = QtWidgets.QGridLayout()
        grid.setColumnStretch(1, 1)
        row = 0

        self.box = ObjectBox(copies=True)
        self.box.setToolTip("The mdcif model (or animation) to draw the arrows on")
        self.box.picked.connect(lambda e: self._guard(self._structure_picked, e))
        grid.addWidget(QtWidgets.QLabel("Object"), row, 0)
        grid.addWidget(self.box, row, 1)
        row += 1

        scopes = QtWidgets.QHBoxLayout()
        self.scope_buttons = {}
        self._scope_group = QtWidgets.QButtonGroup(area)
        for name, tip in (("cluster", "The model's own modes (its cluster's frames)"),
                          ("global", "The modes of all frames (the clustering axes)")):
            button = QtWidgets.QRadioButton(name)
            button.setToolTip(tip)
            button.toggled.connect(lambda on, name=name: on and self._guard(self.set_scope,
                                                                            name))
            self._scope_group.addButton(button)
            self.scope_buttons[name] = button
            scopes.addWidget(button)
        scopes.addStretch(1)
        grid.addWidget(QtWidgets.QLabel("Scope"), row, 0)
        grid.addLayout(scopes, row, 1)
        row += 1

        self.modes_widget = QtWidgets.QWidget()
        self.modes_grid = QtWidgets.QGridLayout(self.modes_widget)
        self.modes_grid.setContentsMargins(0, 0, 0, 0)
        self.modes_grid.setColumnStretch(4, 1)
        self.modes_area = QtWidgets.QScrollArea()
        self.modes_area.setWidget(self.modes_widget)
        self.modes_area.setWidgetResizable(True)
        self.modes_area.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        modes_label = QtWidgets.QLabel("Modes")
        modes_label.setToolTip("Tick the modes to draw; the button sets a mode's colour")
        grid.addWidget(modes_label, row, 0)
        grid.addWidget(self.modes_area, row, 1)
        row += 1

        self.length = Lever(*LENGTH[:3], LENGTH[3], decimals=1, suffix=" sd",
                            tip="Arrow length: how far an atom moves at this many standard "
                                "deviations of the mode ('scale')")
        self.thickness = Lever(*THICKNESS[:3], THICKNESS[3], decimals=2, suffix=" Å",
                               tip="Radius of the arrow shafts ('radius'); the heads are up to "
                                   "2.7 times wider")
        self.min_length = Lever(*MIN_LENGTH[:3], MIN_LENGTH[3], decimals=1, suffix=" Å",
                                tip="Arrows shorter than this are not drawn ('min_length')")
        for label, lever in (("Length", self.length), ("Thickness", self.thickness),
                             ("Min length", self.min_length)):
            lever.moved.connect(lambda *_: guarded_call(self._changed))
            lever.settled.connect(lambda *_: guarded_call(self._flush_log))
            grid.addWidget(QtWidgets.QLabel(label), row, 0)
            grid.addWidget(lever, row, 1)
            row += 1
        layout.addLayout(grid)

        buttons = QtWidgets.QHBoxLayout()
        self.live_box = QtWidgets.QCheckBox("Live")
        self.live_box.setChecked(True)
        self.live_box.setToolTip("Redraw at once on every change (off: Draw applies the "
                                 "settings)")
        self.live_box.toggled.connect(lambda on: self._guard(self._live_toggled, on))
        self.draw_button = QtWidgets.QPushButton("Draw")
        self.draw_button.setToolTip("Draw the ticked modes with these settings now")
        self.draw_button.clicked.connect(lambda: self._guard(self.draw_now))
        self.keep_button = QtWidgets.QPushButton("Keep")
        self.keep_button.setToolTip("Leave these arrows in the scene; the next change draws "
                                    "new ones (to compare settings)")
        self.keep_button.clicked.connect(lambda: self._guard(self.keep))
        self.clear_button = QtWidgets.QPushButton("Clear")
        self.clear_button.setToolTip("Delete the arrows of this window (not the kept ones)")
        self.clear_button.clicked.connect(lambda: self._guard(self.clear))
        self.animate_button = QtWidgets.QPushButton("Animate this mode")
        self.animate_button.setToolTip("Play the lowest ticked mode in the MDCIF Animate window")
        self.animate_button.clicked.connect(lambda: self._guard(self.animate_mode))
        for w in (self.live_box, self.draw_button, self.keep_button, self.clear_button):
            buttons.addWidget(w)
        buttons.addStretch(1)
        buttons.addWidget(self.animate_button)
        layout.addLayout(buttons)

        hidden = QtWidgets.QHBoxLayout()
        self.hidden_label = QtWidgets.QLabel()
        self.hidden_label.setWordWrap(True)
        self.hidden_label.setStyleSheet(f"color: {NOTE_COLOR}")
        self.use_shown_button = QtWidgets.QPushButton()
        self.use_shown_button.clicked.connect(lambda: self._guard(self._use_shown))
        hidden.addWidget(self.hidden_label, stretch=1)
        hidden.addWidget(self.use_shown_button)
        self.hidden_row = QtWidgets.QWidget()
        self.hidden_row.setLayout(hidden)
        hidden.setContentsMargins(0, 0, 0, 0)
        self.hidden_row.hide()
        layout.addWidget(self.hidden_row)
        self.note = Note()
        layout.addWidget(self.note)
        layout.addStretch(1)

    def _color(self, mode):
        found = self.colors.get(mode)
        return found if found is not None else MODE_COLORS[(mode - 1) % len(MODE_COLORS)]

    def _mode_set(self):
        e = self.entry
        if e is None or not e.alive:
            return None
        try:
            return cx.mode_set(e, self.scope)[0]
        except Exception:
            return None

    def _fill_modes(self):
        """One row per mode of the object's set: tick box, variance, sd, colour."""
        while self.modes_grid.count():
            w = self.modes_grid.takeAt(0).widget()
            if w is not None:
                w.hide()
                w.deleteLater()
        self._rows = {}
        ms = self._mode_set()
        if ms is None or ms.n_modes == 0:
            what = "no mdcif object" if self.entry is None else f"no {self.scope} modes"
            self.modes_grid.addWidget(QtWidgets.QLabel(what), 0, 0)
            self._fit_modes(1)
            return
        self.ticked = {m for m in self.ticked if 1 <= m <= ms.n_modes} or {1}
        for k in range(1, ms.n_modes + 1):
            tick = QtWidgets.QCheckBox(str(k))
            tick.setChecked(k in self.ticked)
            tick.toggled.connect(lambda on, k=k: self._guard(self.tick, k, on))
            try:
                sd = f"sd {ms.sigma(k):.2f} Å"
            except ValueError:
                sd = "sd ?"
            color = ColorButton(self._color(k))
            color.setToolTip(f"Colour of the mode {k} arrows")
            color.color_changed.connect(lambda rgb, k=k: self._guard(self.set_color, k, rgb))
            for col, w in enumerate((tick, QtWidgets.QLabel(
                    f"{100 * ms.variance_fraction[k - 1]:.0f}% variance"),
                    QtWidgets.QLabel(sd), color)):
                self.modes_grid.addWidget(w, k - 1, col)
            self._rows[k] = (tick, color)
        self._fit_modes(ms.n_modes)

    def _fit_modes(self, rows):
        """The mode list as tall as its rows, up to MODE_ROWS (more scroll)."""
        for k in range(self.modes_grid.count()):
            w = self.modes_grid.itemAt(k).widget()
            if w is not None:
                w.show()
        self.modes_grid.invalidate()
        height = self.modes_widget.sizeHint().height()
        self.modes_area.setFixedHeight(height * min(rows, MODE_ROWS) // rows)

    # ------------------------------------------------------------------ settings
    def _structure_picked(self, e):
        """Another object chosen (or the chosen one deleted): its modes, new arrows."""
        self._close_live()
        self.entry = e if e is not None and e.alive else None
        self.refresh()
        self._changed()

    def set_scope(self, scope):
        if scope == self.scope:
            return
        self._close_live()
        self.scope = scope
        self.refresh()
        self._changed()

    def tick(self, mode, on=True):
        """Draw mode ``mode`` too (``on``), or no longer."""
        (self.ticked.add if on else self.ticked.discard)(mode)
        self._changed()

    def set_color(self, mode, rgb):
        """The colour of the mode's arrows (its colour button): redrawn at once with Live."""
        self.colors[mode] = tuple(int(c) for c in list(rgb)[:3])
        if mode in self._rows:
            with quiet(self._rows[mode][1]):
                self._rows[mode][1].set_color(self.colors[mode])
        if self.live_box.isChecked() and _exists(self.live.get(mode)):
            self._changed()

    def settings(self):
        """``(scale, radius, min_length)`` from the sliders."""
        return self.length.value(), self.thickness.value(), self.min_length.value()

    def set_settings(self, scale=None, radius=None, min_length=None):
        """Set the sliders (None = keep) and redraw as a slider move would (for scripts)."""
        for lever, value in ((self.length, scale), (self.thickness, radius),
                             (self.min_length, min_length)):
            if value is not None:
                lever.set_value(value)
        self._changed()

    def preset(self, e, scope="cluster", modes=(1,)):
        """Show entry ``e`` with ``scope`` and just ``modes`` ticked (the Animate window's
        Arrows for this mode)."""
        if not self.box.set_current(e):
            cx._error(f"{e.name} has no mdcif data")
        self.ticked = set(modes)
        self._close_live()
        self.entry, self.scope = e, scope
        self.refresh()
        self._changed()

    def _changed(self, *_):
        """Something that changes the arrows: redraw once Qt runs (Live)."""
        if self.live_box.isChecked() and not self._redraw_due:
            self._redraw_due = True
            QtCore.QTimer.singleShot(0, self._redraw_now)
        self._update_note()

    def _redraw_now(self):
        if self._redraw_due and not self.deleted:
            self._redraw_due = False
            self._guard(self.draw)

    def _live_toggled(self, on):
        if on:
            self._changed()
        self._update_note()

    # ------------------------------------------------------------------ drawing
    def _new_name(self, mode):
        """A name for new arrows of ``mode``: the command's, unless that object is there."""
        base = cx.porcupine_name(self.entry, self.scope, mode)
        taken = set(cmd.get_names("all"))
        k = 1
        while (name := base if k == 1 else f"{base}_{k}") in taken:
            k += 1
        return name

    def draw(self):
        """Draw the ticked modes' arrows with the sliders' settings: the window's arrows are
        changed in place, a mode without them gets a new object, an unticked mode's are
        deleted. A mode that has no arrows at these settings is disabled and named in the
        note."""
        self._redraw_due = False
        e = self.entry
        self.drawn, self.failed = {}, {}
        if e is None or not e.alive:
            self._update_note()
            return
        scale, radius, min_length = self.settings()
        for mode in sorted(self.ticked):
            name = self.live.get(mode)
            rgb = tuple(c / 255 for c in self._color(mode))
            try:
                cgo, n, n_short, _ = cx.porcupine_cgo(e, mode, scale, self.scope, min_length,
                                                      radius, rgb)
            except Exception as exc:
                self.failed[mode] = cx.message(exc)
                if _exists(name):
                    cmd.disable(name)
                continue
            if _exists(name):
                cmd.load_cgo(cgo, name, 1, zoom=0)      # state 1 replaced: same object
                cmd.enable(name)
            else:
                name = self.live[mode] = self._new_name(mode)
                cmd.load_cgo(cgo, name, 1, zoom=0)
            self.drawn[mode] = (n, n_short)
        for mode in [m for m in self.live if m not in self.ticked]:
            self._close([self.live.pop(mode)])
        self._log_due = bool(self.drawn)
        if self._log_due:                               # logged when no change follows soon
            self._log_timer.start()
        self._update_note()

    def draw_now(self):
        """The Draw button: draw now and log it."""
        self.draw()
        self._flush_log()

    def keep(self):
        """Leave the window's arrows in the scene; the next change draws new ones."""
        self._flush_log()
        kept = [n for n in self.live.values() if _exists(n)]
        self.live, self.drawn, self.failed = {}, {}, {}
        if kept:
            self._notice = f"Kept {', '.join(kept)}; the next change draws new arrows."
        self._update_note()

    def clear(self):
        """Delete the window's arrows (not the kept ones)."""
        self._log_timer.stop()
        self._log_due = False
        self._close_live()
        self.drawn = {}
        self._update_note()

    def _close_live(self):
        live, self.live = list(self.live.values()), {}
        self._close(live)

    def _close(self, names):
        for name in names:
            if _exists(name):
                cmd.delete(name)

    def _flush_log(self, *_):
        """Show the drawn arrows as 'mdcif_porcupine' commands, one per mode, once a change
        has settled (a pending redraw is done first)."""
        self._log_timer.stop()
        if self._redraw_due:
            self._redraw_now()
        e = self.entry
        if not self._log_due or e is None or not e.alive:
            return
        self._log_due = False
        scale, radius, min_length = self.settings()
        for mode in sorted(m for m in self.drawn if _exists(self.live.get(m))):
            name = self.live[mode]
            default = cx.porcupine_name(e, self.scope, mode)
            text = cx.porcupine_command(e.name, mode, scale, self._color(mode), self.scope,
                                        min_length, radius)
            cx.log_command(text + (f", name={name}" if name != default else ""))

    # ------------------------------------------------------------------ links
    def animate_mode(self):
        """Animate this mode: the Animate window on this object, scope and the lowest ticked
        mode."""
        from .menu import singleton

        if self.entry is None or not self.entry.alive:
            cx._error("no mdcif object is chosen")
        mode = min(self.ticked) if self.ticked else 1
        window = singleton("animate")
        window.preset(self.entry, self.scope, mode)
        return window

    def _shown_copy(self):
        """The enabled animation of the chosen (disabled) model, if any."""
        e = self.entry
        if e is None or not e.alive or e.copy:
            return None
        found = cx.shown_copies(e)
        return found[-1] if found else None

    def _use_shown(self):
        copy = self._shown_copy()
        if copy is not None:
            self.preset(copy, self.scope, sorted(self.ticked))

    # ------------------------------------------------------------------ state
    def refresh(self):
        """Controls from the state: scope buttons, mode rows, buttons, note."""
        e = self.entry
        data = e.data if e is not None else None
        has_global = data is not None and data.global_set is not None
        if self.scope == "global" and not has_global:
            self.scope = "cluster"
        with quiet(*self.scope_buttons.values()):
            self.scope_buttons[self.scope].setChecked(True)
        self.scope_buttons["global"].setEnabled(has_global)
        self._fill_modes()
        self._update_note()

    def _update_note(self):
        """Buttons and the lines under them: the hidden-model note, then what was drawn,
        errors (red), or a notice."""
        has_live = any(_exists(v) for v in self.live.values())
        self.keep_button.setEnabled(has_live)
        self.clear_button.setEnabled(has_live)
        e = self.entry
        alive = e is not None and e.alive
        self.draw_button.setEnabled(alive)
        self.animate_button.setEnabled(alive and bool(self.ticked))
        enabled = set(cmd.get_names("objects", enabled_only=1))
        copy = self._shown_copy() if alive and e.name not in enabled else None
        if copy is not None:
            self.hidden_label.setText(f"{e.name} is hidden; {copy.name} is its animation on "
                                      "screen.")
            self.use_shown_button.setText(f"Draw on {copy.name}")
        self.hidden_row.setVisible(copy is not None)
        if self._error:
            self.note.show_text(self._error, error=True)
            return
        if not alive:
            self.note.show_text("Load an mdcif file (File > Open) to draw its modes' arrows.")
            return
        parts = []
        for mode in sorted(self.drawn):
            n, short = self.drawn[mode]
            name = self.live.get(mode)
            where = f"{name}: " if _exists(name) else f"mode {mode}: "
            parts.append(f"{where}{n} arrow{'' if n == 1 else 's'}" + (
                f" ({short} shorter than {self.min_length.value():g} Å not drawn)"
                if short else ""))
        failed = [f"mode {m}: {why}" for m, why in sorted(self.failed.items())]
        text = " · ".join(parts + failed)
        if not self.live_box.isChecked():
            text = " ".join(t for t in ("Live is off: Draw applies the settings.", text) if t)
        if self._notice:
            text = f"{self._notice} {text}" if text else self._notice
        self.note.show_text(text, error=bool(failed))

    def _poll(self):
        """Arrows deleted by hand, or the model hidden or shown: the note follows."""
        dead = [m for m, n in self.live.items() if not _exists(n)]
        for m in dead:
            del self.live[m]
            self.drawn.pop(m, None)
        enabled = set(cmd.get_names("objects", enabled_only=1))
        seen = (self.entry is not None and self.entry.name in enabled,
                tuple(sorted(self.live.values())), len(enabled))
        if dead or seen != self._watch:
            self._watch = seen
            self._update_note()

    # ------------------------------------------------------------------ cleanup
    def discard_pending(self):
        """Forget a redraw or log line not done yet."""
        self._log_timer.stop()
        self._log_due = self._redraw_due = False

    def cleanup(self):
        """Closing the window keeps its arrows (a change not yet logged is logged now)."""
        self._flush_log()
        self._log_timer.stop()
        self._poll_timer.stop()
        self.box.delete()
