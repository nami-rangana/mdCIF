"""The MDCIF Animate window: play a model along one of its modes, with sliders (GUI only).

    ┌ MDCIF Animate ─────────────────────────────────────────────────┐
    │ Object     [gb1_m1  Cluster 1 medoid (frame 142700, ...)     ▾] │
    │ Scope      (•) cluster  ( ) global     Mode [1 ▾]  26% variance │
    │ Amplitude  ├────────●───────┤ [3.0 sd]       -10 to 10 sd       │
    │ Speed      ├──────●─────────┤ [30 fps]       1 to 60 fps        │
    │ Frames     [30]                                                 │
    │ [Play] [Stop]  Frame ├──────●───────────┤ 12 / 30               │
    │ [Arrows for this mode]  [Plot geometry...]  [Save movie...]     │
    │ gb1_m1_mode1: cluster mode 1 of model 1, ±3 sd, 26% ... playing │
    └─────────────────────────────────────────────────────────────────┘

One per session (MDCIF > Animate..., Plugin > MDCIF, 'mdcif_window animate'). Ported from the
ChimeraX bundle. The first Play, or any change of a setting, makes an animate copy of the chosen
model exactly as 'mdcif_animate' does (that command is shown) and marks it as the window's
(``commands.ANIMATE_MADE``). After that a new mode, scope, amplitude or number of frames rewrites
that copy's states in place (``commands.retune_animation``: the same object, renamed to follow
the mode when that name is free, and it keeps playing); once the change has settled PyMOL's
output says what the copy now plays and the 'mdcif_animate' command that makes a copy like it.
A copy that a plot measures on is not rewritten (the plot's values would no longer match): a new
copy is made instead. Another model also gets a new copy. The copy the window moves off is
deleted if the window made it, else disabled (``commands.retire_copy``, the rule plots follow).
The copy deleted elsewhere: the next Play makes a new one.

The animation plays on PyMOL's movie player, which every object with several states follows:
Play is 'frame 1; mplay', Stop 'mstop', Speed PyMOL's ``movie_fps`` ('set movie_fps, N'). The
frame slider follows the frame shown; dragging it stops the animation on that frame. Save
movie... records one cycle at this speed ('mdcif_movie': GIF always, MP4 / WebM / MOV with
ffmpeg); the window is busy meanwhile, and the animation plays on afterwards if it was playing.

Arrows for this mode opens the Porcupine window on the copy (or the model) with this mode
ticked; Plot geometry... opens the Analyze window on the copy's Geometry tab (its plot, or an
empty tab that takes atoms on the copy). ``preset(entry, scope, mode)`` is what the Porcupine
window's "Animate this mode" calls.
"""

import os

from pymol import cmd
from pymol.Qt import QtCore, QtWidgets

from . import commands as cx
from .menu import ANIMATE
from .widgets import (
    POLL_MS,
    SETTLE_MS,
    ElidedLabel,
    FrameWatcher,
    Guarded,
    Lever,
    Note,
    ObjectBox,
    ToolWindow,
    ask_save_path,
    guarded_call,
    quiet,
    shown_frame,
)

Qt = QtCore.Qt
AMPLITUDE = (-10.0, 10.0, 0.1)          # sd
SPEED = (1, 60, 1)                      # movie_fps
GIF_FILTER = "GIF animation (*.gif)"
VIDEO_FILTERS = ["MP4 video (*.mp4)", "WebM video (*.webm)", "QuickTime movie (*.mov)"]


def _alive(e):
    return e is not None and e.alive


class MdcifAnimate(Guarded, ToolWindow):

    def __init__(self):
        super().__init__(ANIMATE, width=480)
        self.entry = None                               # chosen model
        self.copy = self.source = None                  # the window's copy, the model it plays
        self.scope = cx.ANIMATE_DEFAULTS[1]
        self._apply_due = self._retuned = self._speed_changed = False
        self._status = ""
        self._playing = None
        self._watcher = FrameWatcher(self._frame_changed)
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
        mode, scope, amplitude, frames = cx.ANIMATE_DEFAULTS
        layout = QtWidgets.QVBoxLayout(area)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(4)
        grid = QtWidgets.QGridLayout()
        grid.setColumnStretch(1, 1)

        self.box = ObjectBox()
        self.box.setToolTip("The mdcif model to animate")
        self.box.picked.connect(lambda e: self._guard(self._structure_picked, e))
        grid.addWidget(QtWidgets.QLabel("Object"), 0, 0)
        grid.addWidget(self.box, 0, 1)

        row = QtWidgets.QHBoxLayout()
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
            row.addWidget(button)
        row.addSpacing(12)
        row.addWidget(QtWidgets.QLabel("Mode"))
        self.mode_box = QtWidgets.QComboBox()
        self.mode_box.currentIndexChanged.connect(lambda _: guarded_call(self._mode_changed))
        row.addWidget(self.mode_box)
        self.variance_label = QtWidgets.QLabel()
        row.addWidget(self.variance_label)
        row.addStretch(1)
        grid.addWidget(QtWidgets.QLabel("Scope"), 1, 0)
        grid.addLayout(row, 1, 1)

        self.amplitude = Lever(*AMPLITUDE, amplitude, decimals=1, suffix=" sd",
                               tip="How far the model swings each way, in standard deviations "
                                   "of the mode (negative: the other way first)")
        self.amplitude.moved.connect(lambda _: guarded_call(self._changed))
        self.amplitude.settled.connect(lambda _: guarded_call(self._flush_log))
        fps = int(round(min(max(cmd.get_setting_float("movie_fps"), SPEED[0]), SPEED[1])))
        self.speed = Lever(*SPEED, fps, decimals=0, suffix=" fps",
                           tip="Frames per second (PyMOL's movie_fps); every animation plays "
                               "at this speed")
        self.speed.moved.connect(lambda _: self._guard(self._speed_moved))
        self.speed.settled.connect(lambda _: guarded_call(self._flush_log))
        self.frames_box = QtWidgets.QSpinBox()
        self.frames_box.setRange(cx.MIN_ANIMATE_FRAMES, 1000)
        self.frames_box.setValue(frames)
        self.frames_box.setKeyboardTracking(False)      # a typed number counts on Enter
        self.frames_box.setToolTip("Frames in one back-and-forth cycle")
        self.frames_box.valueChanged.connect(lambda _: guarded_call(self._setting_done))
        frames_row = QtWidgets.QHBoxLayout()
        frames_row.addWidget(self.frames_box)
        frames_row.addStretch(1)
        for k, (label, w) in enumerate((("Amplitude", self.amplitude), ("Speed", self.speed),
                                        ("Frames", frames_row)), start=2):
            grid.addWidget(QtWidgets.QLabel(label), k, 0)
            if isinstance(w, QtWidgets.QHBoxLayout):
                grid.addLayout(w, k, 1)
            else:
                grid.addWidget(w, k, 1)
        layout.addLayout(grid)

        row = QtWidgets.QHBoxLayout()
        self.play_button = QtWidgets.QPushButton("Play")
        self.play_button.setToolTip("Play the animation in a loop, from frame 1 (the first "
                                    "time: make it, as 'mdcif_animate' does)")
        self.play_button.clicked.connect(lambda: self._guard(self.play))
        self.stop_button = QtWidgets.QPushButton("Stop")
        self.stop_button.clicked.connect(lambda: self._guard(self.stop))
        self.frame_slider = QtWidgets.QSlider(Qt.Orientation.Horizontal)
        self.frame_slider.setMinimum(1)
        self.frame_slider.setToolTip("The frame shown; dragging stops the animation there")
        self.frame_slider.valueChanged.connect(lambda k: self._guard(self.go_to_frame, k))
        self.frame_label = QtWidgets.QLabel()
        row.addWidget(self.play_button)
        row.addWidget(self.stop_button)
        row.addSpacing(8)
        row.addWidget(QtWidgets.QLabel("Frame"))
        row.addWidget(self.frame_slider, stretch=1)
        row.addWidget(self.frame_label)
        layout.addLayout(row)

        row = QtWidgets.QHBoxLayout()
        self.arrows_button = QtWidgets.QPushButton("Arrows for this mode")
        self.arrows_button.setToolTip("Porcupine arrows of this mode, in the MDCIF Porcupine "
                                      "window")
        self.arrows_button.clicked.connect(lambda: self._guard(self.arrows_for_mode))
        self.plot_button = QtWidgets.QPushButton("Plot geometry...")
        self.plot_button.setToolTip("Distances, angles or dihedrals along this animation, in "
                                    "the MDCIF Analyze window")
        self.plot_button.clicked.connect(lambda: self._guard(self.plot_geometry))
        self.movie_button = QtWidgets.QPushButton("Save movie...")
        self.movie_button.setToolTip("Record one cycle at this speed as a movie (GIF; MP4, WebM "
                                     "or MOV with ffmpeg)")
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

    def fps(self):
        return int(round(self.speed.value()))

    def set_settings(self, mode=None, scope=None, amplitude=None, frames=None, speed=None):
        """Set the controls (None = keep) and apply them now, as the controls do (scripts);
        ``speed`` in frames per second."""
        if scope is not None:
            self.set_scope(scope)
        with quiet(self.mode_box, self.frames_box):
            if mode is not None:
                if not 1 <= mode <= self.mode_box.count():
                    cx._error(f"mode {mode} not in 1..{self.mode_box.count()}")
                self.mode_box.setCurrentIndex(mode - 1)
                self._show_variance()
            if frames is not None:
                self.frames_box.setValue(frames)
        if amplitude is not None:
            self.amplitude.set_value(amplitude)
        if speed is not None:
            self.speed.set_value(speed)
            self._speed_moved()
        self.apply()
        self._flush_log()

    def _structure_picked(self, e):
        self.entry = e if _alive(e) else None
        self._fill_modes()
        if _alive(self.copy) and self.entry is not None:
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
        if self.entry is None:
            return None
        try:
            return cx.mode_set(self.entry, self.scope)[0]
        except Exception:
            return None

    def _fill_modes(self):
        """The Mode list of the chosen model and scope (global off without a global set)."""
        data = self.entry.data if self.entry is not None else None
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
            QtCore.QTimer.singleShot(0, self._apply_now)

    def _apply_now(self):
        if self._apply_due and not self.deleted:
            self._guard(self.apply)

    def apply(self):
        """Make the animation match the controls: the window's copy rewritten in place, or a
        new copy (the first time, for another model, or when a plot measures on the copy).
        Returns the copy's entry."""
        self._apply_due = False
        e = self.entry
        if e is None or not e.alive:
            cx._error("no mdcif model is loaded: load an mdcif file first")
        mode, scope, amplitude, frames = self.settings()
        copy = self.copy
        if _alive(copy) and self.source is e:
            if cx._matches(copy, mode, scope, amplitude, frames):
                return copy
            if cx.retune_animation(copy, mode, scope, amplitude, frames) is not None:
                self._retuned = True
                self._log_timer.start()
                self.refresh()
                return copy
        return self._new_copy(e, mode, scope, amplitude, frames)

    def _new_copy(self, e, mode, scope, amplitude, frames):
        """'mdcif_animate' (shown) for the window; the copy it moves off is put away."""
        old = self.copy
        copy = cx._animate(e, mode, scope, amplitude, frames, made=cx.ANIMATE_MADE)
        self.copy, self.source, self._retuned = copy, e, False
        if _alive(old) and old is not copy:
            cx.retire_copy(old, old.made == cx.ANIMATE_MADE)
        self._watcher.watch(copy)
        self.refresh()
        return copy

    def _speed_moved(self):
        """A new speed: PyMOL's movie_fps at once (logged once settled)."""
        cmd.set("movie_fps", self.fps())
        self._speed_changed = True
        self.refresh()

    def play(self):
        """Play the animation in a loop from frame 1 (made first if needed). Returns the
        copy."""
        old = self.copy if _alive(self.copy) else None
        copy = self.apply()
        if copy.name not in cmd.get_names("objects", enabled_only=1):
            cmd.enable(copy.name)
            cx.log_command(f"enable {copy.name}")
        if copy is old:                                 # a new copy plays already
            cx._play()
            cx.log_command("frame 1; mplay")
        self.refresh()
        return copy

    def stop(self):
        self._live_copy()
        cmd.mstop()
        cx.log_command("mstop")
        self.refresh()

    def go_to_frame(self, frame):
        """The frame slider: stop the animation on frame ``frame``."""
        copy = self._live_copy()
        frame = min(max(int(frame), 1), copy.n_frames)
        cmd.mstop()
        cmd.frame(frame)
        self._show_frame(frame)
        self.refresh()

    def _live_copy(self):
        if not _alive(self.copy):
            cx._error("nothing is animated yet: click Play")
        return self.copy

    def _flush_log(self, *_):
        """Once a change has settled: say what the copy now plays (rewritten in place), and
        show a new speed as its 'set movie_fps' command."""
        self._log_timer.stop()
        if self._apply_due:
            self._apply_now()
        if self._speed_changed:
            self._speed_changed = False
            cx.log_command(f"set movie_fps, {self.fps()}")
        copy = self.copy
        if not _alive(copy):
            self._retuned = False
            return
        if self._retuned:
            self._retuned = False
            mode, scope, amplitude = copy.mode, copy.scope, copy.amplitude
            frames = copy.n_frames
            ms, model = cx.mode_set(copy, scope)
            again = (cx.animate_command(self.source.name, mode, amplitude, frames, scope)
                     if _alive(self.source) else "")
            print(f" {copy.name} now plays {scope} mode {mode} of model {model}, "
                  f"{100 * ms.variance_fraction[mode - 1]:.0f}% variance, +/-{abs(amplitude):g} "
                  f"sd, {frames} frames (changed in place" + (
                      f"; a new copy like it: {again})" if again else ")"))

    # ------------------------------------------------------------------ movie
    def save_movie(self, path=None):
        """Record one cycle of the animation at the current speed as a movie (GIF; MP4, WebM
        or MOV with ffmpeg; by suffix); without ``path`` a save dialog asks. Returns the path
        or None."""
        copy = self.apply()
        video = cx.ffmpeg_exe() is not None
        if path is None:
            base = os.path.splitext(os.path.basename(copy.path))[0]
            suffix = ".mp4" if video else ".gif"
            filters = (VIDEO_FILTERS + [GIF_FILTER]) if video else [GIF_FILTER]
            path = ask_save_path(self, "Save Movie", filters,
                                 f"{base}_model{copy.model}_{copy.scope}_mode{copy.mode}"
                                 f"{suffix}")
            if not path:
                return None
        if copy.name not in cmd.get_names("objects", enabled_only=1):
            cmd.enable(copy.name)
        n = copy.n_frames
        self._busy(f"Recording one cycle of {copy.name} ({n} frames)...")
        try:
            cx.log_command(f"mdcif_movie {path}, {copy.name}")
            path = cx.mdcif_movie(path, copy.name)
        finally:
            self._busy(None)
        self._notice = f"Movie saved to {path}"
        return path

    def _busy(self, text):
        """While a movie is recorded: the controls off and a note (drawn at once)."""
        self.ui_area.setEnabled(text is None)
        if text is not None:
            self.note.show_text(text)
            self.note.repaint()
            QtWidgets.QApplication.processEvents(
                QtCore.QEventLoop.ProcessEventsFlag.ExcludeUserInputEvents)

    # ------------------------------------------------------------------ links
    def arrows_for_mode(self):
        """Arrows for this mode: the Porcupine window on the copy (or the model, before
        Play), with this mode ticked."""
        from .menu import singleton

        target = self.copy if _alive(self.copy) else self.entry
        if target is None:
            cx._error("no mdcif model is loaded: load an mdcif file first")
        mode, scope, _, _ = self.settings()
        if target is self.copy:
            mode, scope = self.copy.mode, self.copy.scope
        cx.log_command("mdcif_window porcupine")
        window = singleton("porcupine")
        window.preset(target, scope, [mode])
        return window

    def plot_geometry(self):
        """Plot geometry...: the Analyze window on this copy (made first if needed)."""
        from .analyze_tool import plot_geometry

        copy = self.apply()
        return plot_geometry(copy, self.source)

    def preset(self, e, scope="cluster", mode=1):
        """Animate entry ``e`` (a model; an animate copy stands for its model) along
        ``scope`` mode ``mode`` (the Porcupine window's Animate this mode)."""
        if e.copy:
            if e is self.copy and _alive(self.source):
                e = self.source
            else:
                model = cx.model_entry(e.path, e.model)
                if model is None:
                    cx._error(f"{e.name} is an animation of model {e.model}, which is not "
                              "loaded")
                e = model
        if not self.box.set_current(e):
            cx._error(f"{e.name} has no mdcif data")
        self.entry, self.scope = e, scope
        self._fill_modes()
        if not 1 <= mode <= self.mode_box.count():
            cx._error(f"mode {mode} not in 1..{self.mode_box.count()}")
        with quiet(self.mode_box):
            self.mode_box.setCurrentIndex(mode - 1)
        self._show_variance()
        copy = self.apply()
        self._log_timer.start()
        return copy

    # ------------------------------------------------------------------ state
    def _frame_changed(self, e):
        if e is self.copy:
            self._show_frame(shown_frame(e))
            self._update_status()

    def _show_frame(self, frame):
        frame = frame or 1
        with quiet(self.frame_slider):
            self.frame_slider.setValue(frame)
        n = self.copy.n_frames if _alive(self.copy) else self.frames_box.value()
        self.frame_label.setText(f"{frame} / {n}")

    def refresh(self):
        """Buttons, frame slider and status line from the state."""
        copy = self.copy if _alive(self.copy) else None
        has_model = self.entry is not None
        for w in (self.play_button, self.arrows_button, self.plot_button, self.movie_button):
            w.setEnabled(has_model)
        self.stop_button.setEnabled(copy is not None)
        self.frame_slider.setEnabled(copy is not None)
        with quiet(self.frame_slider):
            self.frame_slider.setMaximum(copy.n_frames if copy is not None
                                         else self.frames_box.value())
        if copy is not None:
            self._show_frame(shown_frame(copy))
        else:
            self.frame_label.setText(f"- / {self.frames_box.value()}")
        self._update_status()
        self._update_note()

    def _update_status(self):
        copy = self.copy if _alive(self.copy) else None
        if copy is not None:                            # (also called on every frame)
            try:
                ms, model = cx.mode_set(copy, copy.scope)
                variance = f", {100 * ms.variance_fraction[copy.mode - 1]:.0f}% variance"
            except Exception:
                model, variance = copy.model, ""
            self._playing = cx.playing()
            state = ("playing" if self._playing else f"stopped at frame {shown_frame(copy)}")
            shown = copy.name in cmd.get_names("objects", enabled_only=1)
            text = (f"{copy.name}: {copy.scope} mode {copy.mode} of model {model}, "
                    f"±{abs(copy.amplitude):g} sd{variance} "
                    f"({state}{'' if shown else ', hidden: Play shows it'})")
        elif self.entry is not None:
            text = (f"Play (or a change) animates a copy of {self.entry.name}, as "
                    "'mdcif_animate' does.")
        else:
            text = "Load an mdcif file (File > Open) to animate its modes."
        if text != self._status:
            self._status = text
            self.status.setText(text)

    def _update_note(self):
        self.note.show_text(self._error or self._notice, error=bool(self._error))

    def _poll(self):
        """The copy deleted elsewhere (the next Play makes a new one), or the player started or
        stopped by a typed command: the status follows."""
        if self.copy is not None and not self.copy.alive:
            self.copy = self.source = None
            self._watcher.watch(None)
            self.refresh()
        elif self.copy is not None and cx.playing() != self._playing:
            self._update_status()

    # ------------------------------------------------------------------ cleanup
    def discard_pending(self):
        """Forget a change or log line not done yet."""
        self._log_timer.stop()
        self._apply_due = self._retuned = self._speed_changed = False

    def cleanup(self):
        """Closing the window leaves its animation playing."""
        self._flush_log()
        self._log_timer.stop()
        self._poll_timer.stop()
        self._watcher.watch(None)
        self.box.delete()
