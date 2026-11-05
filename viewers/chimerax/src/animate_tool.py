"""The MDCIF Animate window: play a model along one of its modes, with sliders (GUI only).

    ┌ MDCIF Animate ─────────────────────────────────────────────────┐
    │ Structure  [#1.1 Cluster 1 medoid (frame 142700, ...)        ▾] │
    │ Scope      (•) cluster  ( ) global     Mode [1 ▾]  26% variance │
    │ Amplitude  ├────────●───────┤ [3.0 sd]       -10 to 10 sd       │
    │ Speed      ├●───────────────┤ [1/1]          1/1 to 1/10        │
    │ Frames     [30]                                                 │
    │ [Play] [Stop]  Frame ├──────●───────────┤ 12 / 30               │
    │ [Arrows for this mode]  [Plot geometry...]  [Save movie...]     │
    │ #2: cluster mode 1 of model 1, +/-3 sd, 26% variance (playing)  │
    └─────────────────────────────────────────────────────────────────┘

One per session (MDCIF > Animate..., Tools > Structure Analysis, 'ui tool show "MDCIF
Animate"'). The first Play, or any change of a setting, makes an animate copy of the chosen
model exactly as 'mdcif animate' does (that command is logged) and marks it as the window's
(``cmd.ANIMATE_MADE``). After that a new mode, scope, amplitude or number of frames rewrites
that copy's frames in place (``cmd.retune_animation``: same model id, and it keeps playing);
once the change has settled the Log says what the copy now plays and the 'mdcif animate'
command that makes a copy like it. A copy that a plot measures on is not rewritten (the plot's
values would no longer match): a new copy is made instead. Another structure also gets a new
copy. The copy the window moves off is closed if the window made it, else stopped and hidden
(``cmd.retire_copy``, the rule plots follow). The copy closed elsewhere: the next Play makes a
new one.

Speed is ChimeraX's 'coordset ... pauseFrames N' (each frame shown N times: 1/N of full speed).
Play is 'coordset #N 1,F loop 1000 pauseFrames P' (logged; a loop always starts at frame 1),
Stop 'coordset stop #N'. The frame slider follows the animation; dragging it stops the
animation on that frame. Save movie... records one cycle at the current speed ('movie record',
'coordset', 'wait', 'movie encode', all logged; mp4, webm, mov or avi by suffix); the window is
busy for that cycle, and the animation plays on afterwards if it was playing.

Arrows for this mode opens the Porcupine window on the copy (or the model) with this mode
ticked; Plot geometry... opens the Analyze window on the copy's Geometry tab (its plot, or an
empty tab that takes atoms on the copy). ``preset(structure, scope, mode)`` is what the
Porcupine window's "Animate this mode" calls. Not saved in sessions (the copy is).
"""

import os

from chimerax.core.commands import FileNameArg, log_equivalent_command, run
from chimerax.core.errors import UserError
from chimerax.core.tools import ToolInstance
from Qt.QtCore import Qt, QTimer
from Qt.QtWidgets import (
    QButtonGroup,
    QComboBox,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QRadioButton,
    QSlider,
    QSpinBox,
    QVBoxLayout,
)

from . import cmd
from .menu import ANIMATE, PORCUPINE
from .widgets import (
    SETTLE_MS,
    ElidedLabel,
    FrameWatcher,
    Guarded,
    Lever,
    Note,
    StructureBox,
    ask_save_path,
    quiet,
)

AMPLITUDE = (-10.0, 10.0, 0.1)          # sd
SPEED = (1, 10, 1)                      # pauseFrames: 1/1 to 1/10 of full speed
MOVIE_FILTERS = ["MP4 video (*.mp4)", "WebM video (*.webm)", "QuickTime movie (*.mov)",
                 "AVI video (*.avi)"]
MOVIE_SUFFIXES = (".mp4", ".webm", ".mov", ".avi")
MOVIE_FRAMERATE = 25


class MdcifAnimate(Guarded, ToolInstance):
    SESSION_SAVE = False                # the animate copy is saved in sessions, the window not

    def __init__(self, session, tool_name=ANIMATE):
        from chimerax.core.models import END_DELETE_MODELS
        from chimerax.ui import MainToolWindow

        super().__init__(session, tool_name)
        self.structure = None                           # chosen model
        self.copy = self.source = None                  # the window's copy, the model it plays
        self.scope = cmd.ANIMATE_DEFAULTS[1]
        self._apply_due = self._retuned = self._speed_changed = self._deleted = False
        self._status = ""
        self._watcher = FrameWatcher(self._frame_changed)
        self.tool_window = MainToolWindow(self)
        self._build(self.tool_window.ui_area)
        self._log_timer = QTimer()
        self._log_timer.setSingleShot(True)
        self._log_timer.setInterval(SETTLE_MS)
        self._log_timer.timeout.connect(self._flush_log)
        self._delete_handler = session.triggers.add_handler(END_DELETE_MODELS,
                                                            self._models_deleted)
        self.tool_window.manage(placement=None)
        self._structure_picked(self.box.current())

    # ------------------------------------------------------------------ layout
    def _build(self, area):
        mode, scope, amplitude, frames = cmd.ANIMATE_DEFAULTS
        layout = QVBoxLayout(area)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(4)
        grid = QGridLayout()
        grid.setColumnStretch(1, 1)

        self.box = StructureBox(self.session)
        self.box.setToolTip("The mdcif model to animate")
        self.box.picked.connect(lambda s: self._guard(self._structure_picked, s))
        grid.addWidget(QLabel("Structure"), 0, 0)
        grid.addWidget(self.box, 0, 1)

        row = QHBoxLayout()
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
            row.addWidget(button)
        row.addSpacing(12)
        row.addWidget(QLabel("Mode"))
        self.mode_box = QComboBox()
        self.mode_box.currentIndexChanged.connect(lambda _: self._mode_changed())
        row.addWidget(self.mode_box)
        self.variance_label = QLabel()
        row.addWidget(self.variance_label)
        row.addStretch(1)
        grid.addWidget(QLabel("Scope"), 1, 0)
        grid.addLayout(row, 1, 1)

        self.amplitude = Lever(*AMPLITUDE, amplitude, decimals=1, suffix=" sd",
                               tip="How far the model swings each way, in standard deviations "
                                   "of the mode (negative: the other way first)")
        self.amplitude.moved.connect(lambda _: self._changed())
        self.amplitude.settled.connect(self._flush_log)
        self.speed = Lever(*SPEED, 1, decimals=0, prefix="1/",
                           range_text="1/1 (full speed) to 1/10",
                           tip="Each frame is shown N times (ChimeraX 'coordset ... pauseFrames "
                               "N')")
        self.speed.moved.connect(lambda _: self._guard(self._speed_moved))
        self.speed.settled.connect(self._flush_log)
        self.frames_box = QSpinBox()
        self.frames_box.setRange(cmd.MIN_ANIMATE_FRAMES, 1000)
        self.frames_box.setValue(frames)
        self.frames_box.setKeyboardTracking(False)      # a typed number counts on Enter
        self.frames_box.setToolTip("Frames in one back-and-forth cycle")
        self.frames_box.valueChanged.connect(lambda _: self._setting_done())
        frames_row = QHBoxLayout()
        frames_row.addWidget(self.frames_box)
        frames_row.addStretch(1)
        for k, (label, w) in enumerate((("Amplitude", self.amplitude), ("Speed", self.speed),
                                        ("Frames", frames_row)), start=2):
            grid.addWidget(QLabel(label), k, 0)
            (grid.addLayout if isinstance(w, QHBoxLayout) else grid.addWidget)(w, k, 1)
        layout.addLayout(grid)

        row = QHBoxLayout()
        self.play_button = QPushButton("Play")
        self.play_button.setToolTip("Play the animation in a loop, from frame 1 (the first "
                                    "time: make it, as 'mdcif animate' does)")
        self.play_button.clicked.connect(lambda: self._guard(self.play))
        self.stop_button = QPushButton("Stop")
        self.stop_button.clicked.connect(lambda: self._guard(self.stop))
        self.frame_slider = QSlider(Qt.Orientation.Horizontal)
        self.frame_slider.setMinimum(1)
        self.frame_slider.setToolTip("The frame shown; dragging stops the animation there")
        self.frame_slider.valueChanged.connect(lambda k: self._guard(self.go_to_frame, k))
        self.frame_label = QLabel()
        row.addWidget(self.play_button)
        row.addWidget(self.stop_button)
        row.addSpacing(8)
        row.addWidget(QLabel("Frame"))
        row.addWidget(self.frame_slider, stretch=1)
        row.addWidget(self.frame_label)
        layout.addLayout(row)

        row = QHBoxLayout()
        self.arrows_button = QPushButton("Arrows for this mode")
        self.arrows_button.setToolTip("Porcupine arrows of this mode, in the MDCIF Porcupine "
                                      "window")
        self.arrows_button.clicked.connect(lambda: self._guard(self.arrows_for_mode))
        self.plot_button = QPushButton("Plot geometry...")
        self.plot_button.setToolTip("Distances, angles or dihedrals along this animation, in "
                                    "the MDCIF Analyze window")
        self.plot_button.clicked.connect(lambda: self._guard(self.plot_geometry))
        self.movie_button = QPushButton("Save movie...")
        self.movie_button.setToolTip("Record one cycle at this speed as a movie (mp4, webm, "
                                     "mov or avi)")
        self.movie_button.clicked.connect(lambda: self._guard(self.save_movie))
        for w in (self.arrows_button, self.plot_button, self.movie_button):
            row.addWidget(w)
        row.addStretch(1)
        layout.addLayout(row)
        self.status = ElidedLabel()
        layout.addWidget(self.status)
        self.note = Note()
        layout.addWidget(self.note)
        layout.addStretch(1)
        with quiet(self.scope_buttons[scope]):
            self.scope_buttons[scope].setChecked(True)

    # ------------------------------------------------------------------ settings
    def settings(self):
        """``(mode, scope, amplitude, frames)`` from the controls."""
        return (self.mode_box.currentIndex() + 1 if self.mode_box.count() else 1, self.scope,
                round(self.amplitude.value(), 3), self.frames_box.value())

    def pause_frames(self):
        return int(round(self.speed.value()))

    def set_settings(self, mode=None, scope=None, amplitude=None, frames=None, speed=None):
        """Set the controls (None = keep) and apply them now, as the controls do (scripts)."""
        if scope is not None:
            self.set_scope(scope)
        with quiet(self.mode_box, self.frames_box):
            if mode is not None:
                if not 1 <= mode <= self.mode_box.count():
                    raise UserError(f"mode {mode} not in 1..{self.mode_box.count()}")
                self.mode_box.setCurrentIndex(mode - 1)
            if frames is not None:
                self.frames_box.setValue(frames)
        if amplitude is not None:
            self.amplitude.set_value(amplitude)
        if speed is not None:
            self.speed.set_value(speed)
            self._speed_moved()
        self.apply()
        self._flush_log()

    def _structure_picked(self, structure):
        self.structure = structure if cmd._alive(structure) else None
        self._fill_modes()
        if cmd._alive(self.copy) and self.structure is not None:
            self._changed()                             # a new copy of the new model
        self.refresh()

    def set_scope(self, scope):
        if scope == self.scope:
            return
        self.scope = scope
        self._fill_modes()
        self._setting_done()

    def _mode_changed(self):
        self._show_variance()
        self._setting_done()

    def _setting_done(self):
        """A choice (mode, scope, frames): apply it, log it soon."""
        self._changed()
        self._log_timer.start()

    def _mode_set(self):
        if self.structure is None:
            return None
        try:
            return cmd._mode_set(self.session, self.structure, self.scope)[0]
        except UserError:
            return None

    def _fill_modes(self):
        """The Mode list of the chosen structure and scope (global off without a global set)."""
        data = getattr(self.structure, "mdcif_data", None)
        has_global = data is not None and data.global_set is not None
        if self.scope == "global" and not has_global:
            self.scope = "cluster"
        with quiet(*self.scope_buttons.values(), self.mode_box):
            self.scope_buttons[self.scope].setChecked(True)
            self.scope_buttons["global"].setEnabled(has_global)
            keep = self.mode_box.currentIndex()
            self.mode_box.clear()
            ms = self._mode_set()
            for k in range(1, (ms.n_modes if ms is not None else 0) + 1):
                self.mode_box.addItem(str(k))
                self.mode_box.setItemData(k - 1, f"{100 * ms.variance_fraction[k - 1]:.0f}% "
                                          "of the variance", Qt.ItemDataRole.ToolTipRole)
            self.mode_box.setCurrentIndex(min(max(keep, 0), self.mode_box.count() - 1))
        self._show_variance()

    def _show_variance(self):
        ms = self._mode_set()
        k = self.mode_box.currentIndex()
        self.variance_label.setText(f"{100 * ms.variance_fraction[k]:.0f}% variance"
                                    if ms is not None and 0 <= k < ms.n_modes else "")

    # ------------------------------------------------------------------ the animation
    def _changed(self):
        """A setting moved: apply it once Qt runs (once for any number of moves)."""
        if not self._apply_due:
            self._apply_due = True
            QTimer.singleShot(0, self._apply_now)

    def _apply_now(self):
        if self._apply_due and not self._deleted:
            self._guard(self.apply)

    def apply(self):
        """Make the animation match the controls: the window's copy rewritten in place, or a
        new copy (the first time, for another model, or when a plot measures on the copy).
        Returns the copy."""
        self._apply_due = False
        s = self.structure
        if s is None:
            raise UserError("no mdcif model is open: open an mdcif file first")
        mode, scope, amplitude, frames = self.settings()
        copy = self.copy
        if cmd._alive(copy) and self.source is s:
            if cmd._matches(copy, mode, scope, amplitude, frames):
                return copy
            if cmd.retune_animation(self.session, copy, mode, scope, amplitude,
                                    frames) is not None:
                self._retuned = True
                self._log_timer.start()
                self.refresh()
                return copy
        return self._new_copy(s, mode, scope, amplitude, frames)

    def _new_copy(self, structure, mode, scope, amplitude, frames):
        """'mdcif animate' (logged) for the window; the copy it moves off is put away."""
        old = self.copy
        copy = cmd._animate(self.session, structure, mode, scope, amplitude, frames,
                            mark=cmd.ANIMATE_MADE)
        self.copy, self.source, self._retuned = copy, structure, False
        if cmd._alive(old) and old is not copy:
            cmd.retire_copy(self.session, old, bool(getattr(old, cmd.ANIMATE_MADE, False)))
        if self.pause_frames() > 1:                     # 'mdcif animate' plays at full speed
            self._loop(copy, log=True)
        self._watcher.watch(copy)
        self.refresh()
        return copy

    def _loop(self, copy, log):
        run(self.session, f"coordset stop #{copy.id_string}", log=False)
        run(self.session, cmd.loop_command(copy, self.pause_frames()), log=log)

    def _speed_moved(self):
        """A new speed: a playing animation starts again at it (logged once settled)."""
        if cmd._alive(self.copy) and cmd.player(self.session, self.copy) is not None:
            self._loop(self.copy, log=False)
            self._speed_changed = True
        self.refresh()

    def play(self):
        """Play the animation in a loop from frame 1 (made first if needed). Returns the copy."""
        old = self.copy if cmd._alive(self.copy) else None
        copy = self.apply()
        if not copy.display:
            run(self.session, f"show #{copy.id_string} models")
        if copy is old:                                 # a new copy plays already, at this
            self._loop(copy, log=True)                  # speed ('mdcif animate', _new_copy)
        self._speed_changed = False
        self.refresh()
        return copy

    def stop(self):
        copy = self._live_copy()
        run(self.session, f"coordset stop #{copy.id_string}")
        self.refresh()

    def go_to_frame(self, frame):
        """The frame slider: stop the animation on frame ``frame``."""
        copy = self._live_copy()
        frame = min(max(int(frame), 1), copy.num_coordsets)
        run(self.session, f"coordset stop #{copy.id_string}", log=False)   # 'coordset #N k'
        run(self.session, f"coordset #{copy.id_string} {frame}", log=False)  # alone doesn't
        self._show_frame(frame)                                             # stop a loop
        self.refresh()

    def _live_copy(self):
        if not cmd._alive(self.copy):
            raise UserError("nothing is animated yet: click Play")
        return self.copy

    def _flush_log(self, *_):
        """Once a change has settled: say what the copy now plays (rewritten in place), and
        log a new speed of a playing animation as its 'coordset' command."""
        self._log_timer.stop()
        if self._apply_due:
            self._apply_now()
        copy = self.copy
        if not cmd._alive(copy):
            self._retuned = self._speed_changed = False
            return
        if self._retuned:
            self._retuned = False
            mode, scope, amplitude = copy.mdcif_mode, copy.mdcif_scope, copy.mdcif_amplitude
            frames = copy.num_coordsets
            ms, model = cmd._mode_set(self.session, copy, scope)
            again = (cmd.animate_command(self.source, mode, amplitude, frames, scope)
                     if cmd._alive(self.source) else "")
            self.session.logger.info(
                f"#{copy.id_string} now plays {scope} mode {mode} of model {model}, "
                f"{100 * ms.variance_fraction[mode - 1]:.0f}% variance, +/-{abs(amplitude):g} "
                f"sd, {frames} frames (changed in place" + (
                    f"; a new copy like it: {again})" if again else ")"))
        if self._speed_changed:
            self._speed_changed = False
            if cmd.player(self.session, copy) is not None:
                log_equivalent_command(self.session,
                                       cmd.loop_command(copy, self.pause_frames()))

    # ------------------------------------------------------------------ movie
    def save_movie(self, path=None):
        """Record one cycle of the animation at the current speed as a movie (mp4, webm, mov
        or avi, by suffix); without ``path`` a save dialog asks. Returns the path or None."""
        copy = self.apply()
        if path is None:
            base = os.path.splitext(os.path.basename(copy.mdcif_data.path))[0]
            path = ask_save_path(self.tool_window.ui_area, "Save Movie", MOVIE_FILTERS,
                                 f"{base}_model{copy.mdcif_model}_{copy.mdcif_scope}_mode"
                                 f"{copy.mdcif_mode}.mp4")
            if not path:
                return None
        if os.path.splitext(path)[1].lower() not in MOVIE_SUFFIXES:
            raise UserError(f"{path}: a movie must be {', '.join(MOVIE_SUFFIXES)}")
        folder = os.path.dirname(os.path.abspath(os.path.expanduser(path)))
        if not os.path.isdir(folder) or not os.access(folder, os.W_OK):
            raise UserError(f"{path}: can't write in {folder}")
        spec, frames, pause = f"#{copy.id_string}", copy.num_coordsets, self.pause_frames()
        playing = cmd.player(self.session, copy) is not None
        shown = copy.active_coordset_id
        if not copy.display:
            run(self.session, f"show {spec} models")
        self._busy(f"Recording one cycle of {spec} ({frames * pause} frames)...")
        try:
            run(self.session, f"coordset stop {spec}", log=False)
            run(self.session, "movie record")
            try:
                run(self.session, f"coordset {spec} 1,{frames}"
                    + (f" pauseFrames {pause}" if pause > 1 else ""))
                run(self.session, f"wait {frames * pause}")
            except BaseException:
                run(self.session, "movie abort", log=False)
                raise
            run(self.session, f"movie encode {FileNameArg.unparse(path)} framerate "
                              f"{MOVIE_FRAMERATE} quality good wait true")
        finally:
            self._busy(None)
            if cmd._alive(copy):
                if playing:
                    self._loop(copy, log=False)
                else:
                    copy.active_coordset_id = shown
        self._notice = f"Movie saved to {path}"
        return path

    def _busy(self, text):
        """While a movie is recorded: the controls off and a note (drawn at once)."""
        self.tool_window.ui_area.setEnabled(text is None)
        if text is not None:
            self.note.show_text(text)
            self.note.repaint()

    # ------------------------------------------------------------------ links
    def arrows_for_mode(self):
        """Arrows for this mode: the Porcupine window on the copy (or the model, before
        Play), with this mode ticked."""
        from chimerax.core.commands import StringArg
        from chimerax.core.tools import get_singleton

        from .porcupine_tool import MdcifPorcupine

        target = self.copy if cmd._alive(self.copy) else self.structure
        if target is None:
            raise UserError("no mdcif model is open: open an mdcif file first")
        mode, scope, _, _ = self.settings()
        if target is self.copy:
            mode, scope = self.copy.mdcif_mode, self.copy.mdcif_scope
        run(self.session, f"ui tool show {StringArg.unparse(PORCUPINE)}")
        tool = get_singleton(self.session, MdcifPorcupine, PORCUPINE)
        tool.preset(target, scope, [mode])
        return tool

    def plot_geometry(self):
        """Plot geometry...: the Analyze window on this copy (made first if needed)."""
        from .analyze_tool import plot_geometry

        copy = self.apply()
        return plot_geometry(self.session, copy, self.source)

    def preset(self, structure, scope="cluster", mode=1):
        """Animate ``structure`` (a model; an animate copy stands for its model) along
        ``scope`` mode ``mode`` (the Porcupine window's Animate this mode)."""
        if getattr(structure, "mdcif_copy", False):
            if structure is self.copy and cmd._alive(self.source):
                structure = self.source
            else:
                model = cmd._model_structure(self.session, structure.mdcif_data.path,
                                             structure.mdcif_model)
                if model is None:
                    raise UserError(f"#{structure.id_string} is an animation of model "
                                    f"{structure.mdcif_model}, which is not open")
                structure = model
        if not self.box.set_current(structure):
            raise UserError(f"#{structure.id_string} has no mdcif data")
        self.structure, self.scope = structure, scope
        self._fill_modes()
        if not 1 <= mode <= self.mode_box.count():
            raise UserError(f"mode {mode} not in 1..{self.mode_box.count()}")
        with quiet(self.mode_box):
            self.mode_box.setCurrentIndex(mode - 1)
        self._show_variance()
        copy = self.apply()
        self._log_timer.start()
        return copy

    # ------------------------------------------------------------------ state
    def _frame_changed(self, structure):
        if structure is self.copy:
            self._show_frame(structure.active_coordset_id)
            self._update_status()

    def _show_frame(self, frame):
        with quiet(self.frame_slider):
            self.frame_slider.setValue(frame)
        n = self.copy.num_coordsets if cmd._alive(self.copy) else self.frames_box.value()
        self.frame_label.setText(f"{frame} / {n}")

    def refresh(self):
        """Buttons, frame slider and status line from the state."""
        copy = self.copy if cmd._alive(self.copy) else None
        has_structure = self.structure is not None
        for w in (self.play_button, self.arrows_button, self.plot_button, self.movie_button):
            w.setEnabled(has_structure)
        self.stop_button.setEnabled(copy is not None)
        self.frame_slider.setEnabled(copy is not None)
        with quiet(self.frame_slider):
            self.frame_slider.setMaximum(copy.num_coordsets if copy is not None
                                         else self.frames_box.value())
        if copy is not None:
            self._show_frame(copy.active_coordset_id)
        else:
            self.frame_label.setText(f"- / {self.frames_box.value()}")
        self._update_status()
        self._update_note()

    def _update_status(self):
        copy = self.copy if cmd._alive(self.copy) else None
        if copy is not None:                            # (also called on every frame)
            try:
                ms, model = cmd._mode_set(self.session, copy, copy.mdcif_scope)
                variance = f", {100 * ms.variance_fraction[copy.mdcif_mode - 1]:.0f}% variance"
            except (UserError, IndexError):
                model, variance = copy.mdcif_model, ""
            playing = cmd.player(self.session, copy) is not None
            state = ("playing" if playing else
                     f"stopped at frame {copy.active_coordset_id}")
            text = (f"#{copy.id_string}: {copy.mdcif_scope} mode {copy.mdcif_mode} of model "
                    f"{model}, ±{abs(copy.mdcif_amplitude):g} sd{variance} "
                    f"({state}{'' if copy.display else ', hidden: Play shows it'})")
        elif self.structure is not None:
            text = (f"Play (or a change) animates a copy of #{self.structure.id_string}, as "
                    "'mdcif animate' does.")
        else:
            text = "Open an mdcif file (File > Open) to animate its modes."
        if text != self._status:
            self._status = text
            self.status.setText(text)

    def _update_note(self):
        self.note.show_text(self._error or self._notice, error=bool(self._error))

    def _models_deleted(self, _trigger, _models):
        if self.copy is not None and not cmd._alive(self.copy):
            self.copy = self.source = None              # closed elsewhere: Play makes a new one
            self._watcher.watch(None)
            self.refresh()

    # ------------------------------------------------------------------ cleanup
    def discard_pending(self):
        """Forget a change or Log line not done yet (the bundle is being uninstalled)."""
        self._log_timer.stop()
        self._apply_due = self._retuned = self._speed_changed = False

    def delete(self):
        """Closing the window leaves its animation playing."""
        self._flush_log()
        self._deleted = True
        self._log_timer.stop()
        self._watcher.watch(None)
        self.session.triggers.remove_handler(self._delete_handler)
        self.box.delete()
        super().delete()
