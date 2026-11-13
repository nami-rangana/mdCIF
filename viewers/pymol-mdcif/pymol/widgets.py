"""Qt pieces shared by the MDCIF windows (GUI only; imported only with a window).

    ToolWindow      a floating window docked to PyMOL's main window (QDockWidget); its
                    ``ui_area`` holds the controls, ``display(show)`` / ``delete()``
    Lever           a slider with a number box beside it, kept in step: ``moved`` on every
                    change (live redraws), ``settled`` once a change is done (logging)
    ObjectBox       chooser of the objects with mdcif data, kept up to date (polling)
    ColorButton     a colour swatch that opens Qt's colour dialog; ``color_changed(rgb)``
    FrameWatcher    follows the frame PyMOL shows (polling ``cmd.get_state``)
    ElidedLabel     one-line label cut to its width, the whole text as its tooltip
    Note            the grey / red line under a window's controls
    Guarded         mixin: a button's CmdException goes to the note line, not a traceback
    quiet           sets widgets without their signals
    ask_save_path   Qt's save dialog, with a file type list

Ported from the ChimeraX bundle's widgets.py. PyMOL has no "models added / closed" or "frame
changed" events, so ``ObjectBox`` and ``FrameWatcher`` poll on a QTimer (cheap: object names
and ``cmd.get_state``, which takes no lock).
"""

import os
import traceback
from contextlib import contextmanager

from pymol import CmdException, cmd
from pymol.Qt import QtCore, QtGui, QtWidgets

from . import attach
from .commands import message

Qt = QtCore.Qt
ERROR_COLOR, NOTE_COLOR = "#d03b3b", "#52514e"
SETTLE_MS = 400                         # a change has settled after this long without another
POLL_MS = 500                           # objects added, deleted or renamed
FRAME_MS = 80                           # the frame shown


@contextmanager
def quiet(*widgets):
    """Set the widgets without their signals (which would redraw, animate again, ...)."""
    was = [w.blockSignals(True) for w in widgets]
    try:
        yield
    finally:
        for w, b in zip(widgets, was, strict=True):
            w.blockSignals(b)


def guarded_call(func, *args, **kw):
    """Run ``func`` from a Qt signal or timer: a Python error is printed with its traceback
    instead of reaching Qt (which may abort PyMOL). Returns the result or None."""
    try:
        return func(*args, **kw)
    except CmdException as exc:
        print(f" mdcif: {message(exc)}")
    except Exception:
        traceback.print_exc()
    return None


class ToolWindow(QtWidgets.QDockWidget):
    """A window of the MDCIF menu: floating at first, dockable at the side of PyMOL's main
    window. ``closed`` is emitted once when it is closed (×) or ``delete()``d."""

    closed = QtCore.Signal()

    def __init__(self, title, width=None):
        from .commands import gui_window

        main = gui_window()
        super().__init__(title, main)
        self.setObjectName(title)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
        self.setAllowedAreas(Qt.DockWidgetArea.LeftDockWidgetArea
                             | Qt.DockWidgetArea.RightDockWidgetArea)
        self.ui_area = QtWidgets.QWidget()
        self.setWidget(self.ui_area)
        self._closed = False
        self._width = width
        if main is not None:
            main.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self)
        self.setFloating(True)

    @property
    def title(self):
        return self.windowTitle()

    @title.setter
    def title(self, text):
        self.setWindowTitle(text)

    def place(self):
        """First show: a fitting size, beside the right edge of the main window."""
        self.adjustSize()
        if self._width:
            self.resize(max(self.width(), self._width), self.height())
        main = self.parentWidget()
        if main is not None and main.isVisible():
            geo = main.frameGeometry()
            x = max(geo.right() - self.width() - 40, geo.left())
            self.move(x, geo.top() + 80)

    def display(self, show=True):
        if show:
            if not self.isVisible():
                self.show()
            self.raise_()
            self.activateWindow()
        else:
            self.hide()

    def closeEvent(self, event):                        # Qt's name: overridden
        self._finish()
        super().closeEvent(event)

    def delete(self):
        self._finish()
        self.close()

    def _finish(self):
        if not self._closed:
            self._closed = True
            guarded_call(self.cleanup)
            self.closed.emit()

    def cleanup(self):
        """Stop timers, forget the scene (subclasses)."""

    @property
    def deleted(self):
        return self._closed


class ElidedLabel(QtWidgets.QLabel):
    """A one-line label that never sets the window's width: text too long for it is cut
    with an ellipsis, and the whole text is its tooltip (and what ``text()`` returns)."""

    def __init__(self):
        super().__init__()
        self._full = ""
        self.setSizePolicy(QtWidgets.QSizePolicy.Policy.Ignored,
                           QtWidgets.QSizePolicy.Policy.Preferred)

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


class Note(QtWidgets.QLabel):
    """The line under a window's controls: red for an error, grey for a note, hidden when
    empty. Wraps, so it may make the window taller, never wider (its minimum height follows
    its wrapped height, so the rows above are never squeezed)."""

    def __init__(self):
        super().__init__()
        self.setWordWrap(True)
        self.setSizePolicy(QtWidgets.QSizePolicy.Policy.Ignored,
                           QtWidgets.QSizePolicy.Policy.Preferred)
        self.setAlignment(Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft)
        self.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.hide()

    def show_text(self, text, error=False):
        self.setText(text)
        self.setStyleSheet(f"color: {ERROR_COLOR if error else NOTE_COLOR}")
        self._fit(self._layout_width())
        self.setVisible(bool(text))

    def resizeEvent(self, event):                       # Qt's name: overridden
        super().resizeEvent(event)
        if event.size().width() != event.oldSize().width():
            self._fit(event.size().width())

    def _layout_width(self):
        parent = self.parentWidget()
        if parent is None or parent.layout() is None:
            return self.width()
        margins = parent.layout().contentsMargins()
        return parent.contentsRect().width() - margins.left() - margins.right()

    def _fit(self, width):
        fits = self.heightForWidth(width) if self.text() and width > 0 else 0
        self.setMinimumHeight(max(fits, 0))


class Guarded:
    """Mixin for the windows: ``_guard(func, ...)`` runs a button's action. A CmdException
    goes to the window (``_error``, and the console) instead of a traceback, and
    ``refresh()`` sets the controls back; a notice (``_notice``, grey) is shown until the
    next action. Any other Python error is printed with its traceback (a bug), and shown in
    red. The class provides ``refresh()`` and ``_update_note()``."""

    _error = _notice = ""

    def _guard(self, func, *args, **kw):
        shown, self._notice = bool(self._error or self._notice), ""
        try:
            result = func(*args, **kw)
        except CmdException as exc:
            self._error = message(exc)
            print(f" mdcif: {self._error}")
            guarded_call(self.refresh)
            return None
        except Exception as exc:
            traceback.print_exc()
            self._error = f"internal error: {exc!r} (details in PyMOL's output)"
            guarded_call(self.refresh)
            return None
        self._error = ""
        if shown or self._notice:
            guarded_call(self._update_note)
        return result


def ask_save_path(parent, title, filters, default_name):
    """A file name from a save dialog (starting in the current folder), with the suffix of
    the chosen file type added when none was typed; or None."""
    dialog = QtWidgets.QFileDialog(parent, title)
    dialog.setAcceptMode(QtWidgets.QFileDialog.AcceptMode.AcceptSave)
    dialog.setFileMode(QtWidgets.QFileDialog.FileMode.AnyFile)
    dialog.setNameFilters(filters)
    dialog.selectFile(os.path.join(os.getcwd(), default_name))
    if not dialog.exec_() or not dialog.selectedFiles():
        return None
    path = dialog.selectedFiles()[0]
    if not os.path.splitext(path)[1]:
        path += dialog.selectedNameFilter().split("(*")[-1].rstrip(")")
    return path


class Lever(QtWidgets.QWidget):
    """A slider and a number box kept in step, with the range written beside them.

    ``value()`` / ``set_value(v)`` (no signals). ``moved(value)`` on every change: a slider
    drag, the mouse wheel or arrow keys on it, a number typed and Enter, the box's arrows.
    ``settled(value)`` once a change is done: the slider released, Enter in the box, or
    ``SETTLE_MS`` with no further change (wheel, keys, the box's arrows). Typed numbers count on
    Enter or leaving the box; the box takes any value in range, the slider moves in ``step``s.
    """

    moved = QtCore.Signal(float)
    settled = QtCore.Signal(float)

    def __init__(self, minimum, maximum, step, value, *, decimals=2, suffix="", prefix="",
                 tip="", range_text=None):
        super().__init__()
        self._min, self._step = float(minimum), float(step)
        self._dirty = False
        self.slider = QtWidgets.QSlider(Qt.Orientation.Horizontal)
        n = round((maximum - minimum) / step)
        self.slider.setRange(0, n)
        self.slider.setPageStep(max(1, n // 10))
        self.slider.setMinimumWidth(140)
        self.box = QtWidgets.QDoubleSpinBox()
        self.box.setRange(minimum, maximum)
        self.box.setDecimals(decimals)
        self.box.setSingleStep(step)
        self.box.setSuffix(suffix)
        self.box.setPrefix(prefix)
        self.box.setKeyboardTracking(False)
        self.box.setAlignment(Qt.AlignmentFlag.AlignRight)
        unit = suffix.strip()
        self.range_label = QtWidgets.QLabel(
            range_text or f"{prefix}{minimum:g} to {prefix}{maximum:g}" + (f" {unit}" if unit
                                                                          else ""))
        self.range_label.setStyleSheet(f"color: {NOTE_COLOR}")
        for w in (self.slider, self.box):
            w.setToolTip(tip)
        layout = QtWidgets.QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.slider, stretch=1)
        layout.addWidget(self.box)
        layout.addWidget(self.range_label)
        self._timer = QtCore.QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(SETTLE_MS)
        self._timer.timeout.connect(self._settle)
        self.slider.valueChanged.connect(self._slid)
        self.slider.sliderReleased.connect(self._settle)
        self.box.valueChanged.connect(self._typed)
        self.box.editingFinished.connect(self._settle)
        self.set_value(value)

    def value(self):
        return self.box.value()

    def set_value(self, value):
        with quiet(self.slider, self.box):
            self.box.setValue(value)
            self.slider.setValue(self._position(self.box.value()))

    def _position(self, value):
        return round((value - self._min) / self._step)

    def _slid(self, position):
        with quiet(self.box):
            self.box.setValue(self._min + position * self._step)
        self._changed()

    def _typed(self, value):
        with quiet(self.slider):
            self.slider.setValue(self._position(value))
        self._changed()

    def _changed(self):
        self._dirty = True
        self.moved.emit(self.value())
        if not self.slider.isSliderDown():          # a drag settles on release
            self._timer.start()

    def _settle(self):
        self._timer.stop()
        if self._dirty:
            self._dirty = False
            self.settled.emit(self.value())


class ColorButton(QtWidgets.QPushButton):
    """A small button showing a colour; a click opens Qt's colour dialog. ``color()`` is RGB
    0-255; ``set_color(rgb)`` (no signal); ``color_changed(rgb)`` when the user picks one."""

    color_changed = QtCore.Signal(object)

    def __init__(self, rgb=(255, 128, 0), size=(24, 18)):
        super().__init__()
        self.setFixedSize(*size)
        self._rgb = tuple(int(c) for c in rgb[:3])
        self.clicked.connect(lambda: guarded_call(self._ask))
        self._paint()

    def color(self):
        return self._rgb

    def set_color(self, rgb):
        self._rgb = tuple(int(round(float(c))) for c in list(rgb)[:3])
        self._paint()

    def _paint(self):
        r, g, b = self._rgb
        self.setStyleSheet(f"background-color: rgb({r},{g},{b}); border: 1px solid #52514e;")

    def _ask(self):
        chosen = QtWidgets.QColorDialog.getColor(QtGui.QColor(*self._rgb), self.window(),
                                                 self.toolTip() or "Colour")
        if chosen.isValid():
            self.set_color((chosen.red(), chosen.green(), chosen.blue()))
            self.color_changed.emit(self._rgb)


def entry_label(e):
    """How the windows name an mdcif object: ``gb1_m1  Cluster 1 medoid (...)``, or for an
    animate copy ``gb1_m1_mode1  (cluster mode 1 of model 1)``."""
    if e.copy:
        return f"{e.name}  ({e.scope} mode {e.mode} of model {e.model})"
    label = e.data.model_label(e.model) if e.data is not None else f"model {e.model}"
    return f"{e.name}  {label}"


class ObjectBox(QtWidgets.QComboBox):
    """The objects with mdcif data, oldest first: ``gb1_m1  Cluster 1 medoid (...)``; with
    ``copies``, also the 'mdcif_animate' copies.

    Kept up to date by polling (``POLL_MS``): objects loaded, deleted or renamed. ``current()``
    is the chosen entry (None when there is none); ``set_current(e)`` chooses one without
    signals. ``picked(entry or None)`` when the choice changes: by the user, or because the
    chosen object was deleted (the first one is chosen then). Call ``delete()`` when the window
    closes.
    """

    picked = QtCore.Signal(object)

    def __init__(self, copies=False):
        super().__init__()
        self.copies = copies
        self._entries, self._last, self._snapshot = [], None, None
        self.setSizeAdjustPolicy(
            QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.setMinimumContentsLength(26)
        self.currentIndexChanged.connect(lambda *_: guarded_call(self._index_changed))
        self._timer = QtCore.QTimer(self)
        self._timer.setInterval(POLL_MS)
        self._timer.timeout.connect(lambda: guarded_call(self._poll))
        self._timer.start()
        self.refresh()

    def entries(self):
        """The entries to choose from, oldest first."""
        return [e for e in attach.entries() if e.data is not None
                and (self.copies or not e.copy)]

    def current(self):
        k = self.currentIndex()
        e = self._entries[k] if 0 <= k < len(self._entries) else None
        return e if e is not None and e.alive else None

    def set_current(self, e):
        """Choose entry ``e`` (no ``picked`` signal); False if it is not in the list."""
        if e not in self._entries:
            self.refresh()
        if e not in self._entries:
            return False
        with quiet(self):
            self.setCurrentIndex(self._entries.index(e))
        self._last = e
        return True

    def _snap(self, found):
        return tuple((id(e), e.name, e.mode, e.scope) for e in found)

    def refresh(self):
        found = self.entries()
        self._snapshot = self._snap(found)
        keep = self._last if self._last in found else None
        with quiet(self):
            self.clear()
            self._entries = found
            for e in found:
                label = entry_label(e)
                self.addItem(label)
                self.setItemData(self.count() - 1, label, Qt.ItemDataRole.ToolTipRole)
            self.setCurrentIndex(found.index(keep) if keep is not None else 0 if found else -1)
        self.setEnabled(bool(found))
        if self.current() is not self._last:
            self._index_changed()

    def _poll(self):
        if self._snap(self.entries()) != self._snapshot:
            self.refresh()

    def _index_changed(self, *_):
        self._last = self.current()
        self.picked.emit(self._last)

    def delete(self):
        self._timer.stop()


class FrameWatcher:
    """Follows the frame PyMOL shows: ``callback(entry)`` when ``cmd.get_state()`` changed
    while an entry is watched (polled every ``FRAME_MS``). ``watch(entry)`` (None: stop)."""

    def __init__(self, callback):
        self.callback, self.entry, self._state = callback, None, None
        self._timer = QtCore.QTimer()
        self._timer.setInterval(FRAME_MS)
        self._timer.timeout.connect(lambda: guarded_call(self._poll))

    def watch(self, e):
        self.entry, self._state = e, None
        if e is None:
            self._timer.stop()
        else:
            self._timer.start()

    def _poll(self):
        e = self.entry
        if e is None:
            return
        if not e.alive:
            self.watch(None)
            return
        state = cmd.get_state()
        if state != self._state:
            self._state = state
            self.callback(e)


def shown_frame(e):
    """The frame (state) an animate copy shows: PyMOL's current state, within its states."""
    n = e.n_frames if e is not None else 0
    return min(max(cmd.get_state(), 1), n) if n else None
