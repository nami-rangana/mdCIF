"""ChimeraX-MDCIF bundle: viewer commands and windows for mdcif files.

Self-contained on purpose: depends only on ChimeraX and numpy, never on the
``mdcif`` core package (ChimeraX has its own Python environment).

Custom init (``custom-init = true`` in pyproject.toml): ``initialize`` runs at startup
and puts the "add models" handler of ``attach`` in place, so mdcif files opened any way
get their data, and adds the MDCIF menu to the menu bar (``menu``, with the window);
``finish`` removes both when the bundle is unloaded (reinstall). ChimeraX 1.12 does not call
``finish`` on 'toolshed uninstall' (it deregisters without a session), so ``initialize`` also
watches the toolshed's "bundle uninstalled" trigger and then removes the menu and the handler
and closes the MDCIF windows; otherwise the next 'open' raised ModuleNotFoundError in the
handler and the stale menu entries gave error dialogs.

The windows (Qt, imported only with the GUI), in the MDCIF menu and in Tools > Structure
Analysis: "MDCIF Porcupine" (``porcupine_tool.py``), "MDCIF Animate" (``animate_tool.py``)
and "MDCIF Analyze" (``analyze_tool.py``, whose Geometry tab, ``geometry_tab.py``, is also
opened by ``mdcif plot``). Shared Qt pieces are in ``widgets.py``.
"""

from chimerax.core.toolshed import BundleAPI

__version__ = "0.5.0"

# Command name (as declared in pyproject.toml) -> function name in cmd.py
_COMMANDS = {
    "mdcif": "mdcif",
    "mdcif open": "mdcif_open",
    "mdcif modes": "mdcif_modes",
    "mdcif morph": "mdcif_morph",
    "mdcif animate": "mdcif_animate",
    "mdcif porcupine": "mdcif_porcupine",
    "mdcif plot": "mdcif_plot",
}


_uninstall_handlers = {}                # id(session) -> toolshed "bundle uninstalled" handler


def _watch_uninstall(session, bundle_info):
    """On 'toolshed uninstall' of this bundle: what ``finish`` would do, plus closing the
    windows. The modules are bound here: once uninstalled, ``chimerax.mdcif`` is out of
    ``sys.modules`` and its files are gone, so a relative import would fail."""
    from chimerax.core.toolshed import TOOLSHED_BUNDLE_UNINSTALLED, get_toolshed
    from chimerax.core.triggerset import DEREGISTER

    from . import attach, menu

    def uninstalled(_trigger, bi):
        if getattr(bi, "name", None) != bundle_info.name:
            return None
        menu.close_windows(session)
        attach.unregister(session)
        menu.unregister(session)
        _uninstall_handlers.pop(id(session), None)
        return DEREGISTER

    _unwatch_uninstall(session)
    _uninstall_handlers[id(session)] = get_toolshed().triggers.add_handler(
        TOOLSHED_BUNDLE_UNINSTALLED, uninstalled)


def _unwatch_uninstall(session):
    from chimerax.core.toolshed import get_toolshed

    handler = _uninstall_handlers.pop(id(session), None)
    if handler is not None:
        get_toolshed().triggers.remove_handler(handler)


class _MDCIFAPI(BundleAPI):
    api_version = 1

    @staticmethod
    def initialize(session, bundle_info):
        from . import attach, menu

        attach.register(session)
        menu.register(session)
        _watch_uninstall(session, bundle_info)

    @staticmethod
    def finish(session, bundle_info):
        from . import attach, menu

        _unwatch_uninstall(session)
        attach.unregister(session)
        menu.unregister(session)

    @staticmethod
    def register_command(bi, ci, logger):
        from chimerax.core.commands import register

        from . import cmd

        name = _COMMANDS[ci.name]
        register(ci.name, getattr(cmd, name + "_desc"), getattr(cmd, name), logger=logger)

    @staticmethod
    def start_tool(session, bi, ti):
        """The windows of the MDCIF menu (Qt: GUI only). Porcupine and Animate are one per
        session; Analyze: an empty window, else the newest, else a new one."""
        from chimerax.core.tools import get_singleton

        from .menu import ANIMATE, PORCUPINE

        if ti.name == PORCUPINE:
            from .porcupine_tool import MdcifPorcupine

            return get_singleton(session, MdcifPorcupine, PORCUPINE)
        if ti.name == ANIMATE:
            from .animate_tool import MdcifAnimate

            return get_singleton(session, MdcifAnimate, ANIMATE)
        from .analyze_tool import menu_window

        return menu_window(session)


bundle_api = _MDCIFAPI()
