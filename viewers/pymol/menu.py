"""The MDCIF menu in PyMOL's menu bar, just before Help (GUI only):

    MDCIF > Porcupine...    the window "MDCIF Porcupine" (porcupine_tool.py)
            Animate...      the window "MDCIF Animate"   (animate_tool.py)
            Analyze...      the window "MDCIF Analyze"   (analyze_tool.py)

The same entries are under Plugin > MDCIF. Each entry does what 'mdcif_window porcupine'
(animate, analyze) does, and shows that command as typed: it starts the window or brings it to
the front. Porcupine and Animate are one per session; Analyze is one per plot, and the menu
brings an empty one to the front, else the newest, else opens one.

``register`` (the plugin's init) adds the menu; called again (the plugin loaded twice) it
replaces it, so the menu bar never has two. The Plugin menu entries are added once per PyMOL.
"""

MENU = "MDCIF"
PORCUPINE, ANIMATE, ANALYZE = "MDCIF Porcupine", "MDCIF Animate", "MDCIF Analyze"
WINDOWS = {"porcupine": PORCUPINE, "animate": ANIMATE, "analyze": ANALYZE}
ENTRIES = (
    ("Porcupine...", "porcupine", "Arrows for one or more modes of an mdcif model, with "
                                  "sliders"),
    ("Animate...", "animate", "Play a model along one of its modes, with sliders; save a "
                              "movie"),
    ("Analyze...", "analyze", "Plot distances, angles and dihedrals along a mode animation"),
)

_singletons = {}                        # "porcupine" / "animate" -> the open window


def register():
    """Add the MDCIF menu before Help (replacing an older one) and the Plugin > MDCIF
    entries (once)."""
    from pymol.Qt import QtWidgets

    from .commands import gui_window
    from .widgets import guarded_call

    main = gui_window()
    if main is None:
        return
    bar = main.menuBar()
    _remove_menus(bar)
    menu = QtWidgets.QMenu("&" + MENU, bar)
    menu.setObjectName(MENU)
    menu.setToolTipsVisible(True)
    for label, key, tip in ENTRIES:
        action = menu.addAction(label)
        action.setToolTip(tip)
        action.triggered.connect(lambda *_, key=key: guarded_call(menu_show, key))
    before = next((a for a in bar.actions()
                   if a.text().replace("&", "").strip().lower() == "help"), None)
    if before is not None:
        bar.insertMenu(before, menu)
    else:
        bar.addMenu(menu)
    if not getattr(main, "_mdcif_plugin_entries", False):
        try:
            from pymol import plugins

            for label, key, _ in ENTRIES:
                plugins.addmenuitemqt(f"{MENU}|{label}",
                                      lambda key=key: guarded_call(menu_show, key))
            main._mdcif_plugin_entries = True
        except Exception:                           # an old PyMOL: the menu bar is enough
            pass


def _remove_menus(bar):
    for action in list(bar.actions()):
        menu = action.menu()
        if menu is not None and menu.objectName() == MENU:
            bar.removeAction(action)
            menu.deleteLater()


def menu_show(key):
    """A menu entry: 'mdcif_window KEY', shown as typed."""
    from .commands import log_command

    log_command(f"mdcif_window {key}")
    return show_window(key)


def show_window(key):
    """Start the window ``key`` ("porcupine", "animate", "analyze") or bring it to the front;
    returns it."""
    if key == "analyze":
        from .analyze_tool import menu_window

        window = menu_window()
        window.display(True)
        return window
    window = _singletons.get(key)
    if window is None or window.deleted:
        if key == "porcupine":
            from .porcupine_tool import MdcifPorcupine as Window
        else:
            from .animate_tool import MdcifAnimate as Window
        window = _singletons[key] = Window()
        window.place()
    window.display(True)
    return window


def singleton(key, create=True):
    """The open Porcupine or Animate window (started if ``create``), or None."""
    window = _singletons.get(key)
    if window is not None and not window.deleted:
        return window
    return show_window(key) if create else None


def close_windows():
    """Close every MDCIF window (arrows and copies stay)."""
    for window in list(_singletons.values()):
        if not window.deleted:
            window.delete()
    _singletons.clear()
    try:
        from .analyze_tool import analyze_windows

        for window in analyze_windows():
            window.delete()
    except ImportError:                             # no matplotlib, or no Qt
        pass
