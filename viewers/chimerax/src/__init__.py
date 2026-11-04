"""ChimeraX-MDCIF bundle: viewer commands for mdcif files.

Self-contained on purpose: depends only on ChimeraX and numpy, never on the
``mdcif`` core package (ChimeraX has its own Python environment).

Custom init (``custom-init = true`` in pyproject.toml): ``initialize`` runs at startup
and puts the "add models" handler of ``attach`` in place, so mdcif files opened any way
get their data; ``finish`` removes it when the bundle is unloaded (reinstall).

The tool "MDCIF Plot" (Tools > Structure Analysis, also opened by ``mdcif plot``) is in
``tool.py``, imported only with the GUI.
"""

from chimerax.core.toolshed import BundleAPI

__version__ = "0.4.0"

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


class _MDCIFAPI(BundleAPI):
    api_version = 1

    @staticmethod
    def initialize(session, bundle_info):
        from . import attach

        attach.register(session)

    @staticmethod
    def finish(session, bundle_info):
        from . import attach

        attach.unregister(session)

    @staticmethod
    def register_command(bi, ci, logger):
        from chimerax.core.commands import register

        from . import cmd

        name = _COMMANDS[ci.name]
        register(ci.name, getattr(cmd, name + "_desc"), getattr(cmd, name), logger=logger)

    @staticmethod
    def start_tool(session, bi, ti):
        """Tools > Structure Analysis > MDCIF Plot: an empty plot window (Qt: GUI only)."""
        from .tool import MdcifPlot

        return MdcifPlot(session, ti.name)


bundle_api = _MDCIFAPI()
