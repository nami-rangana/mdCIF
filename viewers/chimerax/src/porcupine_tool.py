"""The MDCIF Porcupine window: porcupine arrows for one or more modes, with sliders (GUI only).

    ┌ MDCIF Porcupine ───────────────────────────────────────────────┐
    │ Structure   [#1.1 Cluster 1 medoid (frame 142700, ...)      ▾] │
    │ Scope       (•) cluster  ( ) global                            │
    │ Modes       [x] 1   26% variance   sd 1.21 Å   [colour]        │
    │             [ ] 2   14% variance   sd 0.89 Å   [colour]        │
    │ Length      ├─────●──────┤ [3.0 sd]     0.5 to 10 sd           │
    │ Thickness   ├──●─────────┤ [0.15 Å]     0.05 to 0.6 Å          │
    │ Min length  ├─●──────────┤ [0.5 Å]      0 to 3 Å               │
    │ [x] Live  [Draw]  [Keep]  [Clear]           [Animate this mode] │
    │ #1.1.1 mode 1: 51 arrows (5 shorter than 0.5 Å not drawn)      │
    └────────────────────────────────────────────────────────────────┘

One per session (MDCIF > Porcupine..., Tools > Structure Analysis, 'ui tool show "MDCIF
Porcupine"'). The arrows are those of 'mdcif porcupine' (``cmd.porcupine_geometry``): Length is
its ``scale`` (sd), Thickness its ``radius`` and Min length its ``minLength``. The structures
are the open mdcif models and their 'mdcif animate' copies (arrows on a copy start from its
first frame, the deposited model, as the command's do).

Every ticked mode is drawn, one arrow model (a Surface) per mode under the structure, in its
colour button's colour: mode 1 the command's orange, modes 2 and up the plot palette. With
Live on, every slider move, tick or colour change redraws the window's arrows in place (the same
models, so their ids stay); a typed number counts on Enter. With Live off, Draw applies the
settings. Keep leaves the arrows as they are and starts new ones on the next change (to compare
settings side by side); Clear closes the window's arrows that were not kept. Another structure
or scope closes the window's arrows (not the kept ones) and draws anew. Closing the window
keeps its arrows.

The Log gets one 'mdcif porcupine ...' line per ticked mode once a change has settled (slider
released, Enter, or 0.4 s without a change), not one per slider step. A hidden structure (the
model an animation was made from) hides its arrows: a note says so and offers the shown
animation of that model. "Animate this mode" opens the Animate window on the structure, scope
and lowest ticked mode. ``preset(structure, scope, modes)`` is what the Animate window's
"Arrows for this mode" calls. Not saved in sessions (the arrows are).
"""

import numpy as np
from chimerax.core.commands import log_equivalent_command, run
from chimerax.core.errors import UserError
from chimerax.core.tools import ToolInstance
from Qt.QtCore import QTimer
from Qt.QtWidgets import (
    QButtonGroup,
    QCheckBox,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QRadioButton,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from . import cmd
from .menu import ANIMATE, PORCUPINE
from .plot import COLORS
from .widgets import NOTE_COLOR, SETTLE_MS, Guarded, Lever, Note, StructureBox, quiet

LENGTH = (0.5, 10.0, 0.1, 3.0)          # min, max, step, start: 'scale', sd
THICKNESS = (0.05, 0.6, 0.01, 0.15)     # 'radius', A
MIN_LENGTH = (0.0, 3.0, 0.1, 0.5)       # 'minLength', A
MODE_ROWS = 6                           # mode rows shown before the list scrolls


def _rgba(hex_color):
    return np.array([int(hex_color[k:k + 2], 16) for k in (1, 3, 5)] + [255], dtype=np.uint8)


# mode 1: the command's orange; modes 2 and up: the plot palette without its orange
MODE_COLORS = [np.array(cmd.PORCUPINE_COLOR, dtype=np.uint8)] + [
    _rgba(c) for c in COLORS if c != "#eb6834"]


class MdcifPorcupine(Guarded, ToolInstance):
    SESSION_SAVE = False                # the arrows are saved in sessions, the window not

    def __init__(self, session, tool_name=PORCUPINE):
        from chimerax.core.models import END_DELETE_MODELS, MODEL_DISPLAY_CHANGED
        from chimerax.ui import MainToolWindow

        super().__init__(session, tool_name)
        self.structure, self.scope = None, "cluster"
        self.ticked = {1}                               # modes drawn
        self.colors = {}                                # mode -> RGBA chosen (else MODE_COLORS)
        self.live = {}                                  # mode -> Surface the sliders move
        self.drawn = {}                                 # mode -> (arrows, too short): last draw
        self.failed = {}                                # mode -> why it has no arrows
        self._rows = {}                                 # mode -> (tick box, colour button)
        self._redraw_due = self._log_due = self._note_due = self._deleted = False
        self.tool_window = MainToolWindow(self)
        self._build(self.tool_window.ui_area)
        self._log_timer = QTimer()
        self._log_timer.setSingleShot(True)
        self._log_timer.setInterval(SETTLE_MS)
        self._log_timer.timeout.connect(self._flush_log)
        triggers = session.triggers
        self._handlers = [triggers.add_handler(END_DELETE_MODELS, self._models_deleted),
                          triggers.add_handler(MODEL_DISPLAY_CHANGED, self._display_changed)]
        self.tool_window.manage(placement=None)
        self._structure_picked(self.box.current())

    # ------------------------------------------------------------------ layout
    def _build(self, area):
        from chimerax.ui.widgets import ColorButton

        self._ColorButton = ColorButton
        layout = QVBoxLayout(area)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(4)
        grid = QGridLayout()
        grid.setColumnStretch(1, 1)
        row = 0

        self.box = StructureBox(self.session, copies=True)
        self.box.setToolTip("The mdcif model (or animation) to draw the arrows on")
        self.box.picked.connect(lambda s: self._guard(self._structure_picked, s))
        grid.addWidget(QLabel("Structure"), row, 0)
        grid.addWidget(self.box, row, 1)
        row += 1

        scopes = QHBoxLayout()
        self.scope_buttons = {}
        self._scope_group = QButtonGroup(area)
        for name, tip in (("cluster", "The model's own modes (its cluster's frames)"),
                          ("global", "The modes of all frames (the clustering axes)")):
            button = QRadioButton(name)
            button.setToolTip(tip)
            button.toggled.connect(lambda on, name=name: on and self._guard(self.set_scope,
                                                                            name))
            self._scope_group.addButton(button)
            self.scope_buttons[name] = button
            scopes.addWidget(button)
        scopes.addStretch(1)
        grid.addWidget(QLabel("Scope"), row, 0)
        grid.addLayout(scopes, row, 1)
        row += 1

        self.modes_widget = QWidget()
        self.modes_grid = QGridLayout(self.modes_widget)
        self.modes_grid.setContentsMargins(0, 0, 0, 0)
        self.modes_grid.setColumnStretch(4, 1)
        self.modes_area = QScrollArea()
        self.modes_area.setWidget(self.modes_widget)
        self.modes_area.setWidgetResizable(True)
        self.modes_area.setFrameShape(QFrame.Shape.NoFrame)
        modes_label = QLabel("Modes")
        modes_label.setToolTip("Tick the modes to draw; the button sets a mode's colour")
        grid.addWidget(modes_label, row, 0)
        grid.addWidget(self.modes_area, row, 1)
        row += 1

        self.length = Lever(*LENGTH[:3], LENGTH[3], decimals=1, suffix=" sd",
                            tip="Arrow length: how far an atom moves at this many standard "
                                "deviations of the mode ('scale')")
        self.thickness = Lever(*THICKNESS[:3], THICKNESS[3], decimals=2, suffix=" Å",
                               tip="Radius of the arrow shafts ('radius'); the heads are 2.7 "
                                   "times wider")
        self.min_length = Lever(*MIN_LENGTH[:3], MIN_LENGTH[3], decimals=1, suffix=" Å",
                                tip="Arrows shorter than this are not drawn ('minLength')")
        for label, lever in (("Length", self.length), ("Thickness", self.thickness),
                             ("Min length", self.min_length)):
            lever.moved.connect(self._changed)
            lever.settled.connect(self._flush_log)
            grid.addWidget(QLabel(label), row, 0)
            grid.addWidget(lever, row, 1)
            row += 1
        layout.addLayout(grid)

        buttons = QHBoxLayout()
        self.live_box = QCheckBox("Live")
        self.live_box.setChecked(True)
        self.live_box.setToolTip("Redraw at once on every change (off: Draw applies the "
                                 "settings)")
        self.live_box.toggled.connect(lambda on: self._guard(self._live_toggled, on))
        self.draw_button = QPushButton("Draw")
        self.draw_button.setToolTip("Draw the ticked modes with these settings now")
        self.draw_button.clicked.connect(lambda: self._guard(self.draw_now))
        self.keep_button = QPushButton("Keep")
        self.keep_button.setToolTip("Leave these arrows in the scene; the next change draws "
                                    "new ones (to compare settings)")
        self.keep_button.clicked.connect(lambda: self._guard(self.keep))
        self.clear_button = QPushButton("Clear")
        self.clear_button.setToolTip("Close the arrows of this window (not the kept ones)")
        self.clear_button.clicked.connect(lambda: self._guard(self.clear))
        self.animate_button = QPushButton("Animate this mode")
        self.animate_button.setToolTip("Play the lowest ticked mode in the MDCIF Animate window")
        self.animate_button.clicked.connect(lambda: self._guard(self.animate_mode))
        for w in (self.live_box, self.draw_button, self.keep_button, self.clear_button):
            buttons.addWidget(w)
        buttons.addStretch(1)
        buttons.addWidget(self.animate_button)
        layout.addLayout(buttons)

        hidden = QHBoxLayout()
        self.hidden_label = QLabel()
        self.hidden_label.setWordWrap(True)
        self.hidden_label.setStyleSheet(f"color: {NOTE_COLOR}")
        self.use_shown_button = QPushButton()
        self.use_shown_button.clicked.connect(lambda: self._guard(self._use_shown))
        hidden.addWidget(self.hidden_label, stretch=1)
        hidden.addWidget(self.use_shown_button)
        self.hidden_row = QWidget()
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
        if not cmd._alive(self.structure):
            return None
        try:
            return cmd._mode_set(self.session, self.structure, self.scope)[0]
        except UserError:
            return None

    def _fill_modes(self):
        """One row per mode of the structure's set: tick box, variance, sd, colour."""
        while self.modes_grid.count():
            w = self.modes_grid.takeAt(0).widget()
            if w is not None:
                w.hide()                                # gone now, deleted when Qt runs
                w.deleteLater()
        self._rows = {}
        ms = self._mode_set()
        if ms is None or ms.n_modes == 0:
            what = "no mdcif structure" if self.structure is None else f"no {self.scope} modes"
            self.modes_grid.addWidget(QLabel(what), 0, 0)
            self._fit_modes(1)
            return
        self.ticked = {m for m in self.ticked if 1 <= m <= ms.n_modes} or {1}
        for k in range(1, ms.n_modes + 1):
            tick = QCheckBox(str(k))
            tick.setChecked(k in self.ticked)
            tick.toggled.connect(lambda on, k=k: self._guard(self.tick, k, on))
            try:
                sd = f"sd {ms.sigma(k):.2f} Å"
            except ValueError:
                sd = "sd ?"
            color = self._ColorButton(max_size=(24, 18))
            color.set_color(self._color(k))
            color.setToolTip(f"Colour of the mode {k} arrows")
            color.color_changed.connect(lambda rgba, k=k: self._guard(self.set_color, k, rgba))
            for col, w in enumerate((tick, QLabel(f"{100 * ms.variance_fraction[k - 1]:.0f}% "
                                                  "variance"), QLabel(sd), color)):
                self.modes_grid.addWidget(w, k - 1, col)
            self._rows[k] = (tick, color)
        self._fit_modes(ms.n_modes)

    def _fit_modes(self, rows):
        """The mode list as tall as its rows, up to MODE_ROWS (more scroll). Rows added to a
        window that is on screen stay hidden until Qt runs, and hidden rows don't count in the
        size hint (which was 0, and the list blank): show them before measuring."""
        for k in range(self.modes_grid.count()):
            w = self.modes_grid.itemAt(k).widget()
            if w is not None:
                w.show()
        self.modes_grid.invalidate()
        height = self.modes_widget.sizeHint().height()
        self.modes_area.setFixedHeight(height * min(rows, MODE_ROWS) // rows)

    # ------------------------------------------------------------------ settings
    def _structure_picked(self, structure):
        """Another structure chosen (or the chosen one closed): its modes, new arrows."""
        self._close_live()
        self.structure = structure if cmd._alive(structure) else None
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

    def set_color(self, mode, rgba):
        """The colour of the mode's arrows (its colour button): changed at once with Live."""
        self.colors[mode] = np.asarray(rgba, dtype=np.uint8)
        surf = self.live.get(mode)
        if self.live_box.isChecked() and cmd._alive(surf):
            surf.color = self.colors[mode]
            self._log_due = True
            self._log_timer.start()
        if mode in self._rows:
            with quiet(self._rows[mode][1]):
                self._rows[mode][1].set_color(self.colors[mode])

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

    def preset(self, structure, scope="cluster", modes=(1,)):
        """Show ``structure`` with ``scope`` and just ``modes`` ticked (the Animate window's
        Arrows for this mode)."""
        if not self.box.set_current(structure):
            raise UserError(f"#{structure.id_string} has no mdcif data")
        self.ticked = set(modes)
        self._close_live()
        self.structure, self.scope = structure, scope
        self.refresh()
        self._changed()

    def _changed(self, *_):
        """Something that changes the arrows: redraw once Qt runs (Live)."""
        if self.live_box.isChecked() and not self._redraw_due:
            self._redraw_due = True
            QTimer.singleShot(0, self._redraw_now)
        self._update_note()

    def _redraw_now(self):
        if self._redraw_due and not self._deleted:
            self._redraw_due = False
            self._guard(self.draw)

    def _live_toggled(self, on):
        if on:
            self._changed()
        self._update_note()

    # ------------------------------------------------------------------ drawing
    def draw(self):
        """Draw the ticked modes' arrows with the sliders' settings: the window's arrows are
        changed in place, a mode without them gets a new arrow model, an unticked mode's are
        closed. A mode that has no arrows at these settings is hidden and named in the note."""
        self._redraw_due = False
        s = self.structure
        self.drawn, self.failed = {}, {}
        if not cmd._alive(s):
            self._update_note()
            return
        scale, radius, min_length = self.settings()
        for mode in sorted(self.ticked):
            surf = self.live.get(mode)
            try:
                va, na, ta, n, n_short, _, _ = cmd.porcupine_geometry(
                    self.session, s, mode, scale, self.scope, min_length, radius)
            except UserError as exc:
                self.failed[mode] = str(exc)
                if cmd._alive(surf):
                    surf.display = False
                continue
            if cmd._alive(surf):
                surf.set_geometry(va, na, ta)
                surf.display = True
            else:
                surf = cmd.arrow_surface(self.session, self.scope, mode)
                surf.set_geometry(va, na, ta)
                s.add([surf])
                self.live[mode] = surf
            surf.color = self._color(mode)
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
        kept = [s for s in self.live.values() if cmd._alive(s)]
        self.live, self.drawn, self.failed = {}, {}, {}
        if kept:
            ids = ", ".join(f"#{s.id_string}" for s in kept)
            self._notice = f"Kept {ids}; the next change draws new arrows."

    def clear(self):
        """Close the window's arrows (not the kept ones)."""
        self._log_timer.stop()
        self._log_due = False
        self._close_live()
        self.drawn = {}
        self._update_note()

    def _close_live(self):
        live, self.live = list(self.live.values()), {}
        self._close(live)

    def _close(self, surfaces):
        alive = [s for s in surfaces if cmd._alive(s)]
        if alive:
            self.session.models.close(alive)

    def _flush_log(self, *_):
        """Log the drawn arrows as 'mdcif porcupine' commands, one per mode, once a change
        has settled (a pending redraw is done first)."""
        self._log_timer.stop()
        if self._redraw_due:
            self._redraw_now()
        if not self._log_due or not cmd._alive(self.structure):
            return
        self._log_due = False
        scale, radius, min_length = self.settings()
        for mode in sorted(m for m in self.drawn if cmd._alive(self.live.get(m))):
            log_equivalent_command(self.session, cmd.porcupine_command(
                self.structure, mode, scale, self._color(mode), self.scope, min_length, radius))

    # ------------------------------------------------------------------ links
    def animate_mode(self):
        """Animate this mode: the Animate window on this structure, scope and the lowest
        ticked mode."""
        from chimerax.core.commands import StringArg
        from chimerax.core.tools import get_singleton

        from .animate_tool import MdcifAnimate

        if not cmd._alive(self.structure):
            raise UserError("no mdcif structure is chosen")
        mode = min(self.ticked) if self.ticked else 1
        run(self.session, f"ui tool show {StringArg.unparse(ANIMATE)}")
        tool = get_singleton(self.session, MdcifAnimate, ANIMATE)
        tool.preset(self.structure, self.scope, mode)
        return tool

    def _shown_copy(self):
        """The shown animation of the chosen (hidden) model, if any."""
        s = self.structure
        if not cmd._alive(s) or getattr(s, "mdcif_copy", False):
            return None
        found = cmd._shown_copies(self.session, s)
        return found[-1] if found else None

    def _use_shown(self):
        copy = self._shown_copy()
        if copy is not None:
            self.preset(copy, self.scope, sorted(self.ticked))

    # ------------------------------------------------------------------ state
    def refresh(self):
        """Controls from the state: scope buttons, mode rows, buttons, note."""
        s = self.structure
        data = getattr(s, "mdcif_data", None) if s is not None else None
        has_global = data is not None and data.global_set is not None
        if self.scope == "global" and not has_global:
            self.scope = "cluster"
        with quiet(*self.scope_buttons.values()):
            self.scope_buttons[self.scope].setChecked(True)
        self.scope_buttons["global"].setEnabled(has_global)
        self._fill_modes()
        self._update_note()

    def _update_note(self):
        """Buttons and the lines under them: the hidden-structure note, then what was drawn,
        errors (red), or a notice."""
        self._note_due = False
        has_live = any(cmd._alive(v) for v in self.live.values())
        self.keep_button.setEnabled(has_live)
        self.clear_button.setEnabled(has_live)
        s = self.structure
        alive = cmd._alive(s)
        self.draw_button.setEnabled(alive)
        self.animate_button.setEnabled(alive and bool(self.ticked))
        hidden = alive and not s.visible
        if hidden:
            copy = self._shown_copy()
            self.hidden_label.setText(f"#{s.id_string} is hidden, and so are its arrows"
                                      + (f"; #{copy.id_string} is its animation on screen."
                                         if copy is not None else
                                         f" ('show #{s.id_string} models' shows them)."))
            self.use_shown_button.setText(f"Draw on #{copy.id_string}" if copy else "")
            self.use_shown_button.setVisible(copy is not None)
        self.hidden_row.setVisible(hidden)
        if self._error:
            self.note.show_text(self._error, error=True)
            return
        if not alive:
            self.note.show_text("Open an mdcif file (File > Open) to draw its modes' arrows.")
            return
        parts = []
        for mode in sorted(self.drawn):
            n, short = self.drawn[mode]
            surf = self.live.get(mode)
            where = f"#{surf.id_string} " if cmd._alive(surf) else ""
            parts.append(f"{where}mode {mode}: {n} arrow{'' if n == 1 else 's'}" + (
                f" ({short} shorter than {self.min_length.value():g} Å not drawn)"
                if short else ""))
        failed = [f"mode {m}: {why}" for m, why in sorted(self.failed.items())]
        text = " · ".join(parts + failed)
        if not self.live_box.isChecked():
            text = " ".join(t for t in ("Live is off: Draw applies the settings.", text) if t)
        if self._notice:
            text = f"{self._notice} {text}" if text else self._notice
        self.note.show_text(text, error=bool(failed))

    def _display_changed(self, _trigger, _models):
        if not self._note_due and not self._deleted:
            self._note_due = True
            QTimer.singleShot(0, self._note_later)

    def _note_later(self):
        if not self._deleted:
            self._update_note()

    def _models_deleted(self, _trigger, _models):
        dead = [m for m, s in self.live.items() if not cmd._alive(s)]
        for m in dead:                                  # closed by hand
            del self.live[m]
            self.drawn.pop(m, None)
        if dead:
            self._update_note()

    # ------------------------------------------------------------------ cleanup
    def discard_pending(self):
        """Forget a redraw or Log line not done yet (the bundle is being uninstalled)."""
        self._log_timer.stop()
        self._log_due = self._redraw_due = False

    def delete(self):
        """Closing the window keeps its arrows (a change not yet logged is logged now)."""
        self._flush_log()
        self._deleted = True
        self._log_timer.stop()
        for handler in self._handlers:
            self.session.triggers.remove_handler(handler)
        self.box.delete()
        super().delete()
