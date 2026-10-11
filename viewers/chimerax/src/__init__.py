"""ChimeraX-MDCIF bundle: viewer commands for mdcif files.

Self-contained on purpose: depends only on ChimeraX and numpy, never on the
``mdcif`` core package (ChimeraX has its own Python environment).
"""

from chimerax.core.toolshed import BundleAPI

__version__ = "0.1.0"

# Command name (as declared in pyproject.toml) -> function name in cmd.py
_COMMANDS = {
    "mdcif": "mdcif",
    "mdcif open": "mdcif_open",
    "mdcif modes": "mdcif_modes",
    "mdcif morph": "mdcif_morph",
    "mdcif animate": "mdcif_animate",
    "mdcif porcupine": "mdcif_porcupine",
}


class _MDCIFAPI(BundleAPI):
    api_version = 1

    @staticmethod
    def register_command(bi, ci, logger):
        from chimerax.core.commands import register

        from . import cmd

        name = _COMMANDS[ci.name]
        register(ci.name, getattr(cmd, name + "_desc"), getattr(cmd, name), logger=logger)


bundle_api = _MDCIFAPI()
