"""matplotlib, the plugin's one optional dependency: found, or installed into PyMOL's Python.

The plots of the MDCIF Analyze window and ``mdcif_plot ... save_image=`` need matplotlib, which
many PyMOLs lack (incentive PyMOL; open-source PyMOL from conda-forge). When the plugin loads
with PyMOL's window (right after Plugin Manager > Install New Plugin, and at each start while
matplotlib is missing) ``ensure_matplotlib`` runs, in the background, PyMOL's own Python::

    <PyMOL's python> -m pip install --only-binary=:all: matplotlib numpy==<PyMOL's numpy>

so PyMOL stays usable; pip's main lines and the result go to PyMOL's output. numpy is pinned to
the version PyMOL has, so pip picks a matplotlib that works with it instead of upgrading the
numpy PyMOL was built against; ``--only-binary`` keeps pip from compiling anything. When PyMOL's
folder isn't writable (an all-users install) pip installs into the user's site-packages by
itself, which is then put on ``sys.path``. A Python owned by the system's package manager
(PEP 668: Ubuntu's or Fedora's own PyMOL) is left alone: the output gives the package manager's
command. ``hint()`` is the command to type by hand for this PyMOL; ``install()`` does the same
work as the background install, but blocks (scripts, tests).

PyMOL's interpreter is looked up from ``sys.prefix``, as ``sys.executable`` can be PyMOL's own
program (``PyMOLWin.exe``, ``PyMOL.app/Contents/MacOS/PyMOL``), and checked by running it. No
PyMOL or Qt imports at module level (Qt only for the background install's timer).
"""

import collections
import importlib
import importlib.metadata
import importlib.util
import os
import queue
import re
import shutil
import site
import subprocess
import sys
import sysconfig
import threading

PACKAGE = "matplotlib"
PIP_TIMEOUT = 900                       # s; a slow connection, a first-time font cache
POLL_MS = 250                           # background install: output and end checked this often
SHOWN = ("Collecting", "Downloading", "Installing collected", "Successfully",
         "Defaulting to user")          # pip's lines worth printing (errors go in the message)
SYSTEM_PACKAGES = (("apt", "sudo apt install python3-matplotlib"),
                   ("dnf", "sudo dnf install python3-matplotlib"),
                   ("zypper", "sudo zypper install python3-matplotlib"),
                   ("pacman", "sudo pacman -S python-matplotlib"))

_state = {"status": None, "message": ""}    # status: None (not tried), "installing", "done",
_jobs = []                                  # "failed", "system"; [timer] of a running install


def version():
    """matplotlib's version if PyMOL's Python has it, else None."""
    if not available():
        return None
    try:
        return importlib.metadata.version(PACKAGE)
    except importlib.metadata.PackageNotFoundError:
        return getattr(sys.modules.get(PACKAGE), "__version__", "?")


def available():
    """True if matplotlib can be imported (also right after an install in this session)."""
    if PACKAGE in sys.modules:
        return True
    importlib.invalidate_caches()
    try:
        return importlib.util.find_spec(PACKAGE) is not None
    except (ImportError, ValueError):
        return False


def installing():
    return _state["status"] == "installing"


def externally_managed():
    """True when the system's package manager owns this Python (PEP 668), so pip may not
    install into it (a virtual environment never is)."""
    if sys.prefix != sys.base_prefix:
        return False
    stdlib = sysconfig.get_paths().get("stdlib")
    return bool(stdlib) and os.path.isfile(os.path.join(stdlib, "EXTERNALLY-MANAGED"))


def python_exe():
    """The path of PyMOL's Python interpreter, or None. Not checked here (``_probe`` is)."""
    exe = sys.executable or ""
    candidates = []
    base = os.path.basename(exe).lower()
    if re.fullmatch(r"python(w)?[0-9.]*(\.exe)?", base):
        if base.startswith("pythonw"):              # no console; python.exe is the same Python
            candidates.append(os.path.join(os.path.dirname(exe), "python.exe"))
        candidates.append(exe)
    major, minor = sys.version_info[:2]
    for prefix in dict.fromkeys((sys.prefix, sys.exec_prefix)):
        if os.name == "nt":
            candidates.append(os.path.join(prefix, "python.exe"))
        else:
            candidates += [os.path.join(prefix, "bin", n)
                           for n in ("python3", "python", f"python{major}.{minor}")]
    return next((os.path.normpath(c) for c in candidates if os.path.isfile(c)), None)


def _quoted(path):
    return f'"{path}"' if " " in path else path


def hint():
    """The command that installs matplotlib for this PyMOL, to type in a terminal."""
    if externally_managed():
        for tool, command in SYSTEM_PACKAGES:
            if shutil.which(tool):
                return command
        return "your system's package manager: its package python3-matplotlib"
    exe = python_exe()
    return f"{_quoted(exe) if exe else '<PyMOL folder>/python'} -m pip install matplotlib"


def pip_command(exe):
    """The pip command line: matplotlib, with numpy held at the version PyMOL has."""
    args = [exe, "-m", "pip", "install", "--disable-pip-version-check", "--no-input",
            "--no-warn-script-location", "--only-binary=:all:", PACKAGE]
    try:
        args.append(f"numpy=={importlib.metadata.version('numpy')}")
    except importlib.metadata.PackageNotFoundError:
        pass
    return args


def _no_window():
    """subprocess options: no console window popping up from PyMOL on Windows."""
    return {"creationflags": subprocess.CREATE_NO_WINDOW} if os.name == "nt" else {}


def _same_path(a, b):
    return os.path.normcase(os.path.realpath(a)) == os.path.normcase(os.path.realpath(b))


def _probe(exe):
    """None if ``exe`` is PyMOL's own Python with pip, else what is wrong."""
    code = ("import sys, importlib.util as u; print(sys.prefix); "
            "print(u.find_spec('pip') is not None)")
    try:
        out = subprocess.run([exe, "-c", code], capture_output=True, text=True, timeout=120,
                             stdin=subprocess.DEVNULL, **_no_window())
    except (OSError, subprocess.SubprocessError) as exc:
        return f"could not run {exe}: {exc}"
    lines = out.stdout.split()
    if out.returncode or len(lines) < 2:
        return f"could not run {exe}: {(out.stderr.strip().splitlines() or ['?'])[-1]}"
    if not _same_path(lines[0], sys.prefix):
        return f"{exe} is another Python than PyMOL's ({lines[0]}, not {sys.prefix})"
    if lines[1] != "True":
        return f"PyMOL's Python ({exe}) has no pip"
    return None


def _stream(args, log):
    """Run ``args``, giving ``log`` pip's main lines. Returns (exit code, last lines)."""
    tail = collections.deque(maxlen=8)
    try:
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                stdin=subprocess.DEVNULL, text=True, errors="replace",
                                **_no_window())
    except OSError as exc:
        return 1, [str(exc)]
    killed = []
    watchdog = threading.Timer(PIP_TIMEOUT, lambda: (killed.append(1), proc.kill()))
    watchdog.start()
    try:
        for line in proc.stdout:
            line = line.strip()
            if line:
                tail.append(line)
                if line.startswith(SHOWN):
                    log(f"   {line}")
        code = proc.wait()
    finally:
        watchdog.cancel()
    if killed:
        tail.append(f"ERROR: stopped after {PIP_TIMEOUT} s")
    return code, list(tail)


def _pip(log):
    """The work of an install up to pip's end (safe in a thread). Returns (ok, message)."""
    exe = python_exe()
    if exe is None:
        return False, f"PyMOL's Python was not found in {sys.prefix}"
    problem = _probe(exe)
    if problem:
        return False, problem
    code, tail = _stream(pip_command(exe), log)
    if code:
        errors = [t for t in tail if t.startswith("ERROR")] or tail[-3:]
        return False, "pip failed: " + " / ".join(errors)
    return True, ""


def _finish(ok, message):
    """After pip (in PyMOL's main thread): make the new package importable and check it."""
    if ok:
        user_site = site.getusersitepackages() if site.ENABLE_USER_SITE else None
        if user_site and os.path.isdir(user_site) and user_site not in sys.path:
            site.addsitedir(user_site)          # created by pip's fallback to a user install
        importlib.invalidate_caches()
        try:
            importlib.import_module(f"{PACKAGE}.figure")
        except Exception as exc:                # e.g. built for another numpy
            ok, message = False, f"matplotlib was installed but does not load: {exc}"
        else:
            message = (f"matplotlib {version()} installed: the Analyze window draws its plots "
                       "(a window opened before now: close it and open it again)")
    _state.update(status="done" if ok else "failed", message=message)
    return ok, message


def install(log=print):
    """Install matplotlib into PyMOL's Python now (blocks: seconds to minutes). ``log`` gets
    pip's main lines. Returns (ok, message)."""
    if available():
        return True, f"matplotlib {version()} is already installed"
    if externally_managed():
        message = f"this PyMOL's Python belongs to the system: install matplotlib with {hint()}"
        _state.update(status="system", message=message)
        return False, message
    _state.update(status="installing", message="")
    return _finish(*_pip(log))


def failure_text(reason):
    return (f" mdcif: matplotlib could not be installed ({reason}). Without it the Analyze "
            "window has no plot (its values and Save data still work). Restart PyMOL to try "
            "again, or close PyMOL and type in a terminal (Command Prompt on Windows):\n"
            f"    {hint()}")


def ensure_matplotlib():
    """With PyMOL's window, when the plugin loads: install matplotlib in the background if
    PyMOL's Python lacks it (once per PyMOL session). PyMOL's output follows the install."""
    if available() or _state["status"] is not None:
        return
    if externally_managed():
        _state.update(status="system", message=hint())
        print(" mdcif: the Analyze window's plots need matplotlib, which this PyMOL lacks. Its "
              f"Python belongs to the system; install it with:\n    {hint()}")
        return
    from pymol.Qt import QtCore

    lines, result = queue.SimpleQueue(), []
    _state.update(status="installing", message="")
    print(" mdcif: the Analyze window's plots need matplotlib, which this PyMOL lacks: "
          "installing it into PyMOL's Python with pip, in the background (PyMOL can be used "
          "meanwhile; it takes a minute or so) ...")

    def work():
        try:
            result.append(_pip(lines.put))
        except Exception as exc:                # never leave the timer waiting
            result.append((False, f"{type(exc).__name__}: {exc}"))

    def poll():                                 # Qt's main thread: print, then finish
        while not lines.empty():
            print(lines.get())
        if not result:
            return
        timer.stop()
        _jobs.clear()
        ok, message = _finish(*result[0])
        print(f" mdcif: {message}" if ok else failure_text(message))

    timer = QtCore.QTimer()
    timer.timeout.connect(poll)
    timer.start(POLL_MS)
    _jobs.append(timer)                         # keeps the timer alive
    threading.Thread(target=work, name="mdcif-install-matplotlib", daemon=True).start()


def missing_text():
    """What the Analyze window says in place of its plot when matplotlib is missing."""
    if installing():
        return ("The plot needs matplotlib, which the mdcif plugin is installing into PyMOL's "
                "Python now (PyMOL's output shows how far it is). Close this window and open "
                "it again once the output says it is installed. The values, the list and Save "
                "data work meanwhile.")
    return ("The plot needs matplotlib in PyMOL's Python, which the mdcif plugin could not "
            "install (PyMOL's output says why). Restart PyMOL to try again, or close PyMOL and "
            "type in a terminal (Command Prompt on Windows):\n"
            f"    {hint()}\nThe values, the list and Save data work without it.")


def status():
    """One line for the 'mdcif' command."""
    v = version()
    if v:
        return f"matplotlib {v}"
    if installing():
        return "matplotlib being installed"
    return f"no matplotlib (to install it, in a terminal: {hint()})"
