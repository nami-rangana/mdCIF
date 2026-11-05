"""The MDCIF menu in ChimeraX's menu bar, just before Help (GUI only):

    MDCIF > Porcupine...    the window "MDCIF Porcupine" (porcupine_tool.py)
            Animate...      the window "MDCIF Animate"   (animate_tool.py)
            Analyze...      the window "MDCIF Analyze"   (analyze_tool.py)

Each entry runs 'ui tool show "MDCIF ..."', logged, which starts the window or brings it to the
front (the windows are also in Tools > Structure Analysis). Porcupine and Animate are one per
session. Analyze is one per plot: 'ui tool show' raises all of them, and the menu then puts an
empty one on top, else the newest.

``register`` (bundle custom init) adds the entries with ``MainWindow.add_menu_entry``, which
makes a missing top-level menu before Help (chimerax/ui/gui.py). At startup ChimeraX
initializes the bundle before its main window exists, so the entries are added on the GUI's
'ready' trigger, as ChimeraX's own open_command adds its File menu entries; a 'devel install'
during a session finds the window and adds them at once.

``unregister`` takes the entries and the menu out of the menu bar at once. ChimeraX 1.12 does
not call a bundle's ``finish`` on 'toolshed uninstall', so ``__init__.py`` also calls it (and
``close_windows``) from the toolshed's "bundle uninstalled" trigger. ChimeraX's
``remove_menu`` only schedules the menu's deletion (``deleteLater``), and a 'devel install' in
the same session would find that dying menu by name and add its entries to it. ``register``
also removes an MDCIF menu that is already there, so the menu bar never has two.

No Qt imports: the menu bar is reached through the main window (Qt objects' own methods), so the
tests can use stand-ins.
"""

MENU = "MDCIF"
PORCUPINE, ANIMATE, ANALYZE = "MDCIF Porcupine", "MDCIF Animate", "MDCIF Analyze"
ENTRIES = (
    ("Porcupine...", PORCUPINE, "Arrows for one or more modes of an mdcif model, with sliders"),
    ("Animate...", ANIMATE, "Play a model along one of its modes, with sliders; save a movie"),
    ("Analyze...", ANALYZE, "Plot distances, angles and dihedrals along a mode animation"),
)

_registered = {}                        # id(session) -> {"actions": [...], "ready": handler}


def register(session):
    """Add the MDCIF menu (now if the main window exists, else when the GUI is ready)."""
    ui = session.ui
    if not getattr(ui, "is_gui", False):
        return
    unregister(session)
    state = _registered[id(session)] = {"actions": [], "ready": None}
    if getattr(ui, "main_window", None) is not None:
        _add_entries(session)
    else:
        state["ready"] = ui.triggers.add_handler("ready", lambda *_: _add_entries(session))


def unregister(session):
    """Take the MDCIF menu out of the menu bar, at once."""
    state = _registered.pop(id(session), None)
    ui = session.ui
    if state is not None and state["ready"] is not None:
        ui.triggers.remove_handler(state["ready"])
    main_window = getattr(ui, "main_window", None)
    if getattr(ui, "is_gui", False) and main_window is not None:
        _remove_menus(main_window)


def is_registered(session):
    return id(session) in _registered


def close_windows(session):
    """Close the open MDCIF windows (the bundle is uninstalled: their code is gone, and a
    window would go on logging commands that no longer exist). Arrows and copies stay."""
    names = {name for _, name, _ in ENTRIES}
    for tool in list(session.tools.list()):
        if getattr(tool, "tool_name", None) in names:
            discard = getattr(tool, "discard_pending", None)
            if discard is not None:                 # no Log line for a command that is gone
                discard()
            tool.delete()


def _add_entries(session):
    state = _registered.get(id(session))
    if state is None:                           # unregistered before the GUI was ready
        return
    main_window = session.ui.main_window
    _remove_menus(main_window)                  # never two MDCIF menus, never entries twice
    state["actions"] = [
        main_window.add_menu_entry([MENU], label, lambda *_, name=name: show_tool(session, name),
                                   tool_tip=tip)
        for label, name, tip in ENTRIES]


def _menus(main_window):
    """The MDCIF menus in the menu bar: found by object name, as ``add_menu_entry`` does."""
    return [c for c in main_window.menuBar().children()
            if hasattr(c, "menuAction") and c.objectName() == MENU]


def _remove_menus(main_window):
    """Remove the MDCIF menus now: their entries (actions belong to the main window, so they
    are deleted too), their place in the menu bar and their name; Qt deletes the menus later."""
    bar = main_window.menuBar()
    for menu in _menus(main_window):
        for action in menu.actions():
            menu.removeAction(action)
            action.deleteLater()
        bar.removeAction(menu.menuAction())
        menu.setObjectName(f"{MENU} (removed)")     # add_menu_entry finds menus by this name
        menu.deleteLater()


def show_tool(session, name):
    """A menu entry: 'ui tool show NAME' (logged); for Analyze, then an empty window on top,
    else the newest. Returns the window."""
    from chimerax.core.commands import StringArg, run

    tool = run(session, f"ui tool show {StringArg.unparse(name)}")
    if name == ANALYZE:
        from .analyze_tool import menu_window

        tool = menu_window(session, create=False) or tool
        if tool is not None:
            tool.display(True)
    return tool
