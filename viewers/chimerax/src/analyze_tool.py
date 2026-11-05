"""The MDCIF Analyze window: analyses of an mdcif file, one tab each (GUI only).

    ┌ MDCIF Analyze: distance, angle on #2 ─────────────────┐
    │ [Geometry]                                            │
    │ (the Geometry tab: geometry_tab.GeometryTab)          │
    └───────────────────────────────────────────────────────┘

Geometry, the only tab so far, is the plot of distances, angles and dihedrals along a mode
animation (the "MDCIF Plot" window until 0.4.0); its title is the window's. One window per
plot, as before: 'mdcif plot' (with the GUI) shows its plot in the window showing it, else in
an empty one, else a new one, and switches to the Geometry tab (``show_plot``). MDCIF >
Analyze... (also Tools > Structure Analysis, 'ui tool show "MDCIF Analyze"') opens a new
window when none is open, else brings an empty one to the front, else the newest
(``menu_window``). The Animate window's Plot geometry... button calls ``plot_geometry``.
Windows are not saved in sessions.

For scripts: ``windows(session)`` lists the Geometry tabs, oldest window first, with the API
described in geometry_tab.py (``plot``, ``figure``, ``canvas``, ``go_to_frame``, ``add_text``,
``set_settings``, ...); ``analyze_windows(session)`` lists the windows themselves.
"""

from itertools import count

from chimerax.core.tools import ToolInstance
from Qt.QtWidgets import QTabWidget, QVBoxLayout, QWidget

from .geometry_tab import GeometryTab
from .menu import ANALYZE
from .plot import plots

_serial = count(1)                      # windows in the order they were opened


def analyze_windows(session):
    """The open Analyze windows, oldest first (the session keeps its tools in a set)."""
    return sorted(session.tools.find_by_class(MdcifAnalyze), key=lambda t: t.serial)


def windows(session):
    """The Geometry tabs of the open Analyze windows, oldest first."""
    return [t.geometry for t in analyze_windows(session)]


def show_plot(session, plot):
    """Show ``plot`` on a Geometry tab: the one showing it, else an empty one, else a new
    window's."""
    found = windows(session)
    tab = (next((t for t in found if t.plot is plot), None)
           or next((t for t in found if t.plot is None), None))
    if tab is None:
        return MdcifAnalyze(session, plot=plot).show_geometry()
    if tab.plot is not plot:
        tab.set_plot(plot)
    return tab.host.show_geometry()


def menu_window(session, create=True):
    """The window MDCIF > Analyze... shows: the newest empty one, else the newest one, else
    (``create``) a new one."""
    found = analyze_windows(session)
    empty = [t for t in found if t.geometry.plot is None]
    if empty or found:
        return (empty or found)[-1]
    return MdcifAnalyze(session) if create else None


def plot_geometry(session, copy, source=None):
    """The Animate window's Plot geometry... button: the newest plot measuring on the animate
    copy ``copy``, else an empty Geometry tab ready to take atoms on it (``source``: the model
    it animates). Returns the window."""
    found = [p for p in plots(session) if p.alive and p.structure is copy]
    if found:
        return show_plot(session, found[-1]).host
    tab = next((t for t in windows(session) if t.plot is None), None)
    host = tab.host if tab is not None else MdcifAnalyze(session)
    host.geometry.set_target(copy, source)
    host.show_geometry()
    return host


class MdcifAnalyze(ToolInstance):
    SESSION_SAVE = False                # the animate copy is saved in sessions, the window not

    def __init__(self, session, tool_name=ANALYZE, plot=None):
        from chimerax.ui import MainToolWindow

        super().__init__(session, tool_name)
        self.serial = next(_serial)
        self.tool_window = MainToolWindow(self)
        area = self.tool_window.ui_area
        layout = QVBoxLayout(area)
        layout.setContentsMargins(0, 0, 0, 0)
        self.tabs = QTabWidget()
        layout.addWidget(self.tabs)
        page = QWidget()
        self.geometry = GeometryTab(self, page)
        self.tabs.addTab(page, "Geometry")
        self.tabs.setTabToolTip(0, "Distances, angles and dihedrals at every frame of a mode "
                                   "animation")
        self.tool_window.manage(placement=None)         # floating: a plot needs the width
        self.geometry.set_plot(plot)

    def set_title(self, title):
        self.tool_window.title = title

    def show_geometry(self):
        """Bring the window to the front on its Geometry tab; returns the tab."""
        self.tabs.setCurrentIndex(self.tabs.indexOf(self.geometry.page))
        self.display(True)
        return self.geometry

    def delete(self):
        """Closing the window: its tabs' handlers go, and the plot with them (the copy
        stays)."""
        self.geometry.cleanup()
        super().delete()
