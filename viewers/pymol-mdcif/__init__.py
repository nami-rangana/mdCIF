"""PyMOL plugin for mdcif files: commands and windows (the ChimeraX bundle's, for PyMOL).

Usage (each model of the file becomes an object, in a group named after the file)::

    load gb1.cif                                group gb1: gb1_m1, gb1_m2, ... (or mdcif_load)
    mdcif_modes gb1                             list models, populations and modes again
    mdcif_porcupine gb1_m1, mode=1, scale=3     CGO arrows gb1_m1_mode1_porcupine
    mdcif_animate gb1_m1, mode=1                new object gb1_m1_mode1 swinging along mode 1
    mdcif_animate gb1_m2, mode=1, scope=global  ... along global PC1 (the clustering axis)
    mdcif_morph gb1, frames=20                  new object gb1_morph through the models
    mdcif_plot distance, gb1_m1///10/CA, gb1_m1///40/CA     a distance at every frame of a mode
    mdcif_movie gb1_m1_mode1.gif, gb1_m1_mode1  one cycle of the animation as a movie
    mdcif_window porcupine                      a window of the MDCIF menu (also animate, analyze)
    mdcif                                       the plugin's status

``amplitude`` / ``scale`` are in standard deviations of the mode (k in k * sigma * v, sigma in
angstroms: sqrt(lambda) for PCA, with the elements' standard masses taken out for quasi-harmonic
modes). Self-contained: needs only PyMOL and numpy, and matplotlib for the plots of the Analyze
window and plot images, which the plugin installs into PyMOL's Python with pip when it loads
with PyMOL's window and PyMOL lacks it (``deps``).

Modules: ``attach`` (the load hook, the registry of mdcif objects, sessions), ``commands``
(the commands and their Qt-free helpers), ``plot`` (plots: data, drawing, files),
``mdcif_data`` (the file, shared with the ChimeraX bundle), ``deps`` (matplotlib). The windows
(Qt, imported only with the GUI) are in the MDCIF menu (``menu.py``) and under Plugin > MDCIF:
"MDCIF Porcupine"
(``porcupine_tool.py``), "MDCIF Animate" (``animate_tool.py``) and "MDCIF Analyze"
(``analyze_tool.py``, whose Geometry tab, ``geometry_tab.py``, is also opened by
``mdcif_plot``); shared Qt pieces in ``widgets.py``.
"""

__version__ = "0.2.0"

from .commands import (  # noqa: E402  (the modules import __version__ from here)
    mdcif,
    mdcif_animate,
    mdcif_load,
    mdcif_modes,
    mdcif_morph,
    mdcif_movie,
    mdcif_plot,
    mdcif_porcupine,
    mdcif_window,
)

_OBJECT_ARG = ("mdcif_modes", "mdcif_morph", "mdcif_animate", "mdcif_porcupine")


def __init_plugin__(app=None):
    """The plugin's init (PyMOL's plugin manager, or by hand): the commands, the load hook,
    the session tasks, and with PyMOL's window the MDCIF menu and, if PyMOL lacks it, an
    install of matplotlib in the background. Safe to call again."""
    from pymol import cmd

    from . import attach, commands, deps

    attach.register()
    for f in (mdcif, mdcif_load, mdcif_modes, mdcif_morph, mdcif_animate, mdcif_porcupine,
              mdcif_plot, mdcif_movie, mdcif_window):
        cmd.extend(f.__name__, f)
    cmd.extend("mdcif_open", mdcif_load)
    try:                                            # Tab completes object names
        for name in _OBJECT_ARG:
            cmd.auto_arg[0][name] = cmd.auto_arg[0]["enable"]
    except Exception:
        pass
    gui = commands.gui_window() is not None
    if gui:
        from . import menu

        menu.register()
    print(f" mdcif_pymol {__version__} loaded: mdcif_load, mdcif_modes, mdcif_morph, "
          "mdcif_animate, mdcif_porcupine, mdcif_plot, mdcif_movie, mdcif_window (type 'mdcif' "
          "for more)")
    if gui:
        deps.ensure_matplotlib()
