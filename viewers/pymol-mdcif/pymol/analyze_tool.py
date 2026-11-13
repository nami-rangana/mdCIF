"""The MDCIF Analyze window: analyses of an mdcif file, one tab each (GUI only).

    ┌ MDCIF Analyze: distance, angle on gb1_m1_mode1 ───────┐
    │ [Geometry]                                            │
    │ (the Geometry tab: geometry_tab.GeometryTab)          │
    └───────────────────────────────────────────────────────┘

Ported from the ChimeraX bundle. Geometry, the only tab so far, is the plot of distances, angles
and dihedrals along a mode animation; its title is the window's. One window per plot:
'mdcif_plot' (with the GUI) shows its plot in the window showing it, else in an empty one, else
a new one, and switches to the Geometry tab (``show_plot``). MDCIF > Analyze... (also
Plugin > MDCIF, 'mdcif_window analyze') opens a new window when none is open, else brings an
empty one to the front, else the newest (``menu_window``). The Animate window's Plot
geometry... button calls ``plot_geometry``. Windows are not saved in sessions.

For scripts: ``windows()`` lists the Geometry tabs, oldest window first, with the API described
in geometry_tab.py (``plot``, ``figure``, ``canvas``, ``go_to_frame``, ``add_text``,
``set_settings``, ...); ``analyze_windows()`` lists the windows themselves.
"""

from itertools import count

from pymol.Qt import QtWidgets

from .geometry_tab import GeometryTab
from .menu import ANALYZE
from .plot import plots
from .widgets import ToolWindow

_serial = count(1)                      # windows in the order they were opened
_WINDOWS = []


def analyze_windows():
    """The open Analyze windows, oldest first."""
    _WINDOWS[:] = [w for w in _WINDOWS if not w.deleted]
    return list(_WINDOWS)


def windows():
    """The Geometry tabs of the open Analyze windows, oldest first."""
    return [w.geometry for w in analyze_windows()]


def show_plot(plot):
    """Show ``plot`` on a Geometry tab: the one showing it, else an empty one, else a new
    window's."""
    found = windows()
    tab = (next((t for t in found if t.plot is plot), None)
           or next((t for t in found if t.plot is None), None))
    if tab is None:
        return MdcifAnalyze(plot=plot).show_geometry()
    if tab.plot is not plot:
        tab.set_plot(plot)
    return tab.host.show_geometry()


def menu_window(create=True):
    """The window MDCIF > Analyze... shows: the newest empty one, else the newest one, else
    (``create``) a new one."""
    found = analyze_windows()
    empty = [w for w in found if w.geometry.plot is None]
    if empty or found:
        return (empty or found)[-1]
    return MdcifAnalyze() if create else None


def plot_geometry(copy, source=None):
    """The Animate window's Plot geometry... button: the newest plot measuring on the animate
    copy ``copy`` (an entry), else an empty Geometry tab ready to take atoms on it
    (``source``: the model it animates). Returns the window."""
    found = [p for p in plots() if p.alive and p.entry is copy]
    if found:
        return show_plot(found[-1]).host
    tab = next((t for t in windows() if t.plot is None), None)
    host = tab.host if tab is not None else MdcifAnalyze()
    host.geometry.set_target(copy, source)
    host.show_geometry()
    return host


class MdcifAnalyze(ToolWindow):

    def __init__(self, plot=None):
        super().__init__(ANALYZE, width=760)
        self.serial = next(_serial)
        layout = QtWidgets.QVBoxLayout(self.ui_area)
        layout.setContentsMargins(0, 0, 0, 0)
        self.tabs = QtWidgets.QTabWidget()
        layout.addWidget(self.tabs)
        page = QtWidgets.QWidget()
        self.geometry = GeometryTab(self, page)
        self.tabs.addTab(page, "Geometry")
        self.tabs.setTabToolTip(0, "Distances, angles and dihedrals at every frame of a mode "
                                   "animation")
        _WINDOWS.append(self)
        self.geometry.set_plot(plot)
        self.place()
        k = len(analyze_windows()) - 1                  # several: not all on one spot
        if k:
            self.move(self.x() - 30 * k, self.y() + 30 * k)

    def set_title(self, title):
        self.title = title

    def show_geometry(self):
        """Bring the window to the front on its Geometry tab; returns the tab."""
        self.tabs.setCurrentIndex(self.tabs.indexOf(self.geometry.page))
        self.display(True)
        return self.geometry

    def cleanup(self):
        """Closing the window: its tabs' timers go, and the plot with them (the copy stays)."""
        self.geometry.cleanup()
