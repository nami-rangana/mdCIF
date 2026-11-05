"""Qt pieces shared by the MDCIF windows (GUI only; imported only with the window).

    Lever           a slider with a number box beside it, kept in step: ``moved`` on every
                    change (live redraws), ``settled`` once a change is done (logging)
    StructureBox    chooser of the open structures with mdcif data, kept up to date
    FrameWatcher    follows the frame (active coordinate set) a structure shows
    ElidedLabel     one-line label cut to its width, the whole text as its tooltip
    Note            the grey / red line under a window's controls
    Guarded         mixin: a button's UserError goes to the note line, not a dialog
    quiet           sets widgets without their signals
    ask_save_path   Qt's own save dialog (as ChimeraX's), with a file type list

The ChimeraX colour button is ``chimerax.ui.widgets.ColorButton``; its ``FloatSlider`` has no
number box, hence ``Lever``.
"""

import os
from contextlib import contextmanager

from chimerax.core.errors import UserError
from Qt.QtCore import Qt, QTimer, Signal
from Qt.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QSizePolicy,
    QSlider,
    QWidget,
)

ERROR_COLOR, NOTE_COLOR = "#d03b3b", "#52514e"
SETTLE_MS = 400                         # a change has settled after this long without another


@contextmanager
def quiet(*widgets):
    """Set the widgets without their signals (which would redraw, animate again, ...)."""
    was = [w.blockSignals(True) for w in widgets]
    try:
        yield
    finally:
        for w, b in zip(widgets, was, strict=True):
            w.blockSignals(b)


class ElidedLabel(QLabel):
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


class Note(QLabel):
    """The line under a window's controls: red for an error, grey for a note, hidden when
    empty. Wraps, so it may make the window taller, never wider.

    A floating tool window's minimum size ignores a wrapped label's height for its width, so
    the window would not grow and the layout would squeeze the rows above the note (they then
    overlapped). The note therefore sets its minimum height to its wrapped height: when its
    text is set, at the width of its parent's layout (never narrower than the note, so never
    too tall), and again whenever its width changes."""

    def __init__(self):
        super().__init__()
        self.setWordWrap(True)
        self.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        self.setAlignment(Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft)   # under the
        self.hide()                                     # buttons when the window is taller

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
    """Mixin for the windows: ``_guard(func, ...)`` runs a button's action. A UserError goes to
    the window (``_error``, and the Log as a warning) instead of a dialog, and ``refresh()``
    sets the controls back; a notice (``_notice``, grey) is shown until the next action. The
    class provides ``session``, ``refresh()`` and ``_update_note()``."""

    _error = _notice = ""

    def _guard(self, func, *args, **kw):
        shown, self._notice = bool(self._error or self._notice), ""
        try:
            result = func(*args, **kw)
        except UserError as exc:
            self._error = str(exc)
            self.session.logger.warning(self._error)
            self.refresh()
            return None
        self._error = ""
        if shown or self._notice:
            self._update_note()
        return result


def ask_save_path(parent, title, filters, default_name):
    """A file name from a save dialog (Qt's own, as ChimeraX's, starting in the current
    folder), with the suffix of the chosen file type added when none was typed; or None."""
    dialog = QFileDialog(parent, title)
    dialog.setAcceptMode(QFileDialog.AcceptMode.AcceptSave)
    dialog.setFileMode(QFileDialog.FileMode.AnyFile)
    dialog.setOption(QFileDialog.Option.DontUseNativeDialog)
    dialog.setNameFilters(filters)
    dialog.selectFile(os.path.join(os.getcwd(), default_name))
    if not dialog.exec() or not dialog.selectedFiles():
        return None
    path = dialog.selectedFiles()[0]
    if not os.path.splitext(path)[1]:
        path += dialog.selectedNameFilter().split("(*")[-1].rstrip(")")
    return path


class Lever(QWidget):
    """A slider and a number box kept in step, with the range written beside them.

    ``value()`` / ``set_value(v)`` (no signals). ``moved(value)`` on every change: a slider
    drag, the mouse wheel or arrow keys on it, a number typed and Enter, the box's arrows.
    ``settled(value)`` once a change is done: the slider released, Enter in the box, or
    ``SETTLE_MS`` with no further change (wheel, keys, the box's arrows). Typed numbers count on
    Enter or leaving the box; the box takes any value in range, the slider moves in ``step``s.
    """

    moved = Signal(float)
    settled = Signal(float)

    def __init__(self, minimum, maximum, step, value, *, decimals=2, suffix="", prefix="",
                 tip="", range_text=None):
        super().__init__()
        self._min, self._step = float(minimum), float(step)
        self._dirty = False
        self.slider = QSlider(Qt.Orientation.Horizontal)
        n = round((maximum - minimum) / step)
        self.slider.setRange(0, n)
        self.slider.setPageStep(max(1, n // 10))
        self.slider.setMinimumWidth(140)
        self.box = QDoubleSpinBox()
        self.box.setRange(minimum, maximum)
        self.box.setDecimals(decimals)
        self.box.setSingleStep(step)
        self.box.setSuffix(suffix)
        self.box.setPrefix(prefix)
        self.box.setKeyboardTracking(False)
        self.box.setAlignment(Qt.AlignmentFlag.AlignRight)
        unit = suffix.strip()
        self.range_label = QLabel(range_text or f"{prefix}{minimum:g} to {prefix}{maximum:g}"
                                  + (f" {unit}" if unit else ""))
        self.range_label.setStyleSheet(f"color: {NOTE_COLOR}")
        for w in (self.slider, self.box):
            w.setToolTip(tip)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.slider, stretch=1)
        layout.addWidget(self.box)
        layout.addWidget(self.range_label)
        self._timer = QTimer(self)
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


class StructureBox(QComboBox):
    """The open structures with mdcif data, by id: ``#1.1 Cluster 1 medoid (...)``; with
    ``copies``, also the 'mdcif animate' copies.

    Kept up to date one Qt tick after models are added, closed or renamed (so the "add models"
    handler of ``attach`` has run first). ``current()`` is the chosen structure (None when there
    is none); ``set_current(s)`` chooses one without signals. ``picked(structure or None)`` when
    the choice changes: by the user, or because the chosen structure was closed (the first one
    is chosen then). Call ``delete()`` when the window closes.
    """

    picked = Signal(object)

    def __init__(self, session, copies=False):
        from chimerax.core.models import ADD_MODELS, MODEL_NAME_CHANGED, REMOVE_MODELS

        super().__init__()
        self.session, self.copies = session, copies
        self._structures, self._last, self._due = [], None, False
        self.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.setMinimumContentsLength(26)
        self.currentIndexChanged.connect(self._index_changed)
        triggers = session.triggers
        self._handlers = [triggers.add_handler(name, self._models_changed)
                          for name in (ADD_MODELS, REMOVE_MODELS, MODEL_NAME_CHANGED)]
        self.refresh()

    def structures(self):
        """The structures to choose from, by id."""
        from chimerax.atomic import AtomicStructure

        return sorted((s for s in self.session.models.list(type=AtomicStructure)
                       if getattr(s, "mdcif_data", None) is not None
                       and not getattr(s, "deleted", False)
                       and (self.copies or not getattr(s, "mdcif_copy", False))),
                      key=lambda s: s.id)

    def current(self):
        k = self.currentIndex()
        s = self._structures[k] if 0 <= k < len(self._structures) else None
        return None if s is None or getattr(s, "deleted", False) else s

    def set_current(self, structure):
        """Choose ``structure`` (no ``picked`` signal); False if it is not in the list."""
        if structure not in self._structures:
            self.refresh()
        if structure not in self._structures:
            return False
        with quiet(self):
            self.setCurrentIndex(self._structures.index(structure))
        self._last = structure
        return True

    def refresh(self):
        self._due = False
        found = self.structures()
        keep = self._last if self._last in found else None
        with quiet(self):
            self.clear()
            self._structures = found
            for s in found:
                self.addItem(f"#{s.id_string} {s.name}")
                self.setItemData(self.count() - 1, f"#{s.id_string} {s.name}",
                                 Qt.ItemDataRole.ToolTipRole)
            self.setCurrentIndex(found.index(keep) if keep is not None else 0 if found else -1)
        self.setEnabled(bool(found))
        if self.current() is not self._last:
            self._index_changed()

    def _index_changed(self, *_):
        self._last = self.current()
        self.picked.emit(self._last)

    def _models_changed(self, *_):
        if not self._due:
            self._due = True
            QTimer.singleShot(0, self._refresh_due)

    def _refresh_due(self):
        if self._handlers is not None:          # not closed in the meantime
            self.refresh()

    def delete(self):
        for handler in self._handlers or []:
            self.session.triggers.remove_handler(handler)
        self._handlers = None


class FrameWatcher:
    """Follows the frame a structure shows: ``callback(structure)`` once Qt runs after its
    active coordinate set changed (once for any number of changes before then). ``watch(s)``
    (None: stop following; also drops a pending call)."""

    def __init__(self, callback):
        self.callback, self.structure = callback, None
        self._handler, self._due = None, False

    def watch(self, structure):
        if structure is self.structure:
            return
        if self._handler is not None:
            self._handler.remove()
        self.structure, self._handler = structure, None
        if structure is not None:
            self._handler = structure.triggers.add_handler("changes", self._changed)

    def _changed(self, _trigger, data):
        structure, changes = data
        if ("active_coordset changed" in changes.structure_reasons()
                and structure is self.structure and not self._due):
            self._due = True
            QTimer.singleShot(0, self._fire)

    def _fire(self):
        self._due = False
        s = self.structure
        if s is not None and not getattr(s, "deleted", False):
            self.callback(s)
