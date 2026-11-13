"""Make the figures of tutorial.md in PyMOL's window, by the tutorial's own commands.

Run it in PyMOL with its window (the windows need the GUI), with the script's full path (PyMOL
starts in a folder of its own)::

    pymol /path/to/mdcif/examples/view/pymol/make_figures.py

It loads the plugin from ``viewers/pymol-mdcif`` (no install needed), types the commands of each
section of the tutorial as you would (``cmd.do``: they are shown in PyMOL's output), uses the
MDCIF windows through their controls, and saves every picture and movie in
``examples/view/pymol/figures/``; then PyMOL quits. PyMOL's output is also written to
``figures/make_figures.log`` (not kept in the repository), which is where the output quoted in
the tutorial comes from. The figures of the tutorial were made with open-source PyMOL 3.1
(incentive PyMOL without a license file writes "For Evaluation Only" on every picture).

Scenes are ray-traced at 700 x 606 pixels, movies are GIFs at 480 x 416, windows are grabbed
from the screen buffer (Qt's own grab, plus the 3D view's frame buffer for the main window).
The steps run on Qt timers, so the movie player, the windows and their timers work as they do
for a user; a full run takes about two minutes.
"""

import os
import sys
import time
import traceback

from pymol import cmd
from pymol.Qt import QtCore, QtGui, QtWidgets

HERE = os.path.dirname(os.path.abspath(sys._getframe().f_code.co_filename))  # PyMOL's 'run'
                                                                            # sets no __file__
REPO = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
FIG = os.path.join(HERE, "figures")
GB1 = os.path.join(REPO, "examples", "create", "1_basic", "gb1_meld_300K.cif")
SCENE = (700, 606)                      # ray-traced pictures
MOVIE = (480, 416)                      # GIF movies
os.makedirs(FIG, exist_ok=True)
_log = open(os.path.join(FIG, "make_figures.log"), "w", encoding="utf-8")


class _Tee:
    """PyMOL's output, also into the log file."""

    def __init__(self, other):
        self.other = other

    def write(self, text):
        _log.write(text)
        _log.flush()
        try:
            return self.other.write(text)
        except Exception:
            return len(text)

    def flush(self):
        try:
            self.other.flush()
        except Exception:
            pass

    def __getattr__(self, name):
        return getattr(self.other, name)


sys.stdout = _Tee(sys.stdout)


# --------------------------------------------------------------------------- helpers

def do(text):
    """Type ``text`` at PyMOL's command line (shown, logged, run now)."""
    for line in text.strip().splitlines():
        cmd.do(line.strip(), echo=1, flush=1)


def section(title):
    """A marker in the log file (not in PyMOL's output, which the pictures show)."""
    _log.write(f"\n######## {title}\n")


def spin(ms):
    """Let PyMOL's window run for ``ms`` milliseconds (frames play, timers fire)."""
    app = QtWidgets.QApplication.instance()
    end = time.time() + ms / 1000
    while time.time() < end:
        app.processEvents(QtCore.QEventLoop.ProcessEventsFlag.AllEvents, 20)
        time.sleep(0.005)


def scene(name, size=SCENE):
    """A ray-traced picture of the scene."""
    cmd.png(os.path.join(FIG, name), size[0], size[1], ray=1, quiet=1)


def _main():
    from pymol import gui

    return gui.get_qtwindow()


def grab_main():
    """PyMOL's main window as a QPixmap: Qt's grab of the widgets, with the 3D view (and
    PyMOL's object panel, drawn in it) from its frame buffer."""
    main = _main()
    spin(300)
    pixmap = main.grab()
    gl = main.pymolwidget
    frame = (gl.grabFramebuffer() if hasattr(gl, "grabFramebuffer")
             else gl.grabFrameBuffer())
    painter = QtGui.QPainter(pixmap)
    painter.drawImage(gl.mapTo(main, QtCore.QPoint(0, 0)), frame)
    painter.end()
    return pixmap


def save_widget(widget, name):
    spin(300)
    widget.grab().save(os.path.join(FIG, name))


def side_by_side(window, name, height=None):
    """A window grab next to a draw of the scene at the window's height."""
    spin(300)
    shot = window.grab()
    h = height or max(shot.height(), 420)
    w = int(h * 1.1)
    data = cmd.png(None, w, h, ray=1, quiet=1)
    scene_img = QtGui.QImage.fromData(data, "PNG")
    out = QtGui.QPixmap(w + 12 + shot.width(), max(h, shot.height()))
    out.fill(QtGui.QColor("white"))
    painter = QtGui.QPainter(out)
    painter.drawImage(0, 0, scene_img)
    painter.drawPixmap(w + 12, 0, shot)
    painter.end()
    out.save(os.path.join(FIG, name))


def menu_picture(name):
    """The MDCIF menu open: the top of the main window with the menu drawn under its title."""
    main = _main()
    bar = main.menuBar()
    action = next(a for a in bar.actions() if a.text().replace("&", "") == "MDCIF")
    menu = action.menu()
    menu.popup(bar.mapToGlobal(bar.actionGeometry(action).bottomLeft()))
    menu.setActiveAction(menu.actions()[0])
    spin(400)
    shot_menu = menu.grab()
    menu.hide()
    top = grab_main()
    x = bar.actionGeometry(action).left()
    y = bar.mapTo(main, bar.actionGeometry(action).bottomLeft()).y()
    width = min(top.width(), x + shot_menu.width() + 160)
    out = top.copy(0, 0, width, y + shot_menu.height() + 40)
    painter = QtGui.QPainter(out)
    painter.drawPixmap(x, y, shot_menu)
    painter.end()
    out.save(os.path.join(FIG, name))


def window(key):
    from mdcif_pymol import menu

    return menu.singleton(key, create=False)


# --------------------------------------------------------------------------- the tutorial

S = {}


def s2_plugin():
    section("2. Add the plugin")
    import importlib.util

    spec = importlib.util.spec_from_file_location(  # the plugin as the package mdcif_pymol
        "mdcif_pymol", os.path.join(REPO, "viewers", "pymol-mdcif", "__init__.py"))
    mdcif_pymol = importlib.util.module_from_spec(spec)
    sys.modules["mdcif_pymol"] = mdcif_pymol
    spec.loader.exec_module(mdcif_pymol)
    mdcif_pymol.__init_plugin__()
    main = _main()
    main.resize(1100, 800)
    do("mdcif")


def s2_menu():
    menu_picture("mdcif_menu.png")


def s3_open():
    section("3. Open the file")
    cmd.cd(REPO)
    do("load examples/create/1_basic/gb1_meld_300K.cif, gb1")
    do("""color 0x2a78d6, gb1_m1
          color 0xeb6834, gb1_m2
          bg_color white
          orient gb1""")
    cmd.set("ray_opaque_background", 1)
    S["view"] = cmd.get_view()
    cmd.group("gb1", action="open")                 # the Objects panel lists gb1_m1, gb1_m2
    spin(500)
    grab_main().save(os.path.join(FIG, "window_open.png"))


def s4_porcupine():
    section("4.1 porcupine")
    do("""disable gb1_m2
          mdcif_porcupine gb1_m1, mode=1""")
    scene("porcupine_1_mode_1.png")
    section("4.2 options")
    do("""delete gb1_m1_mode1_porcupine
          mdcif_porcupine gb1_m1, mode=2, scale=2, color=forest""")
    scene("porcupine_1_mode_2.png")
    do("""delete gb1_m1_mode2_porcupine
          disable gb1_m1
          enable gb1_m2
          mdcif_porcupine gb1_m2, mode=1, scope=global, color=0x6a3d9a""")
    scene("porcupine_2_global_1.png")


def s4_window():
    section("4.3 window")
    do("""delete gb1_m2_global1_porcupine
          disable gb1_m2
          enable gb1_m1
          color 0xc8c8c8, gb1_m1""")
    from mdcif_pymol import menu
    menu.menu_show("porcupine")
    spin(800)


def s4_window_modes():
    p = window("porcupine")
    p._rows[2][0].setChecked(True)
    p._rows[3][0].setChecked(True)
    spin(1200)
    side_by_side(p, "window_porcupine.png")


def s4_window_mode2():
    p = window("porcupine")
    p._rows[1][0].setChecked(False)
    p._rows[3][0].setChecked(False)
    p.length.box.setValue(2.0)
    p.length.box.editingFinished.emit()
    spin(300)
    p.set_color(2, (0x22, 0x8b, 0x22))
    spin(1200)
    print(f"[note] {p.note.text()}")
    p.min_length.box.setValue(3.0)
    p.length.box.setValue(0.5)
    p.length.box.editingFinished.emit()
    spin(1200)
    print(f"[note] {p.note.text()}")
    p.min_length.box.setValue(0.5)
    p.length.box.setValue(3.0)
    p.length.box.editingFinished.emit()
    spin(1200)
    p.clear()
    p.delete()
    do("color 0x2a78d6, gb1_m1")


def s5_animate():
    section("5.1 animate")
    do("mdcif_animate gb1_m1, mode=1")
    spin(500)
    do("mdcif_movie examples/view/pymol/figures/model_1_mode_1.gif, gb1_m1_mode1, "
       "width=480, height=416, ray=1")
    do("""mstop
          frame 8""")
    do("""delete gb1_m1_mode1
          disable gb1_m1
          enable gb1_m2
          mdcif_animate gb1_m2, mode=1""")
    spin(500)
    do("mdcif_movie examples/view/pymol/figures/model_2_mode_1.gif, gb1_m2_mode1, "
       "width=480, height=416, ray=1")


def s5_window():
    section("5.2 window")
    do("""delete gb1_m2_mode1
          disable gb1_m2
          enable gb1_m1""")
    from mdcif_pymol import menu
    menu.menu_show("animate")
    spin(800)
    a = window("animate")
    print(f"[status] {a.status.text()}")
    a.play_button.click()
    spin(1500)
    print(f"[status] {a.status.text()}")
    side_by_side(a, "window_animate.png", height=360)


def s5_window_steps():
    a = window("animate")
    a.speed.box.setValue(15)
    a.speed.box.editingFinished.emit()
    spin(600)
    a.frame_slider.setValue(8)
    spin(500)
    print(f"[status] {a.status.text()}")
    a.play_button.click()
    spin(500)
    a.arrows_button.click()
    spin(1000)
    p = window("porcupine")
    print(f"[porcupine note] {p.note.text()}")
    p.clear()
    p.delete()
    a.plot_button.click()
    spin(800)
    from mdcif_pymol import analyze_tool
    tab = analyze_tool.windows()[-1]
    print(f"[analyze title] {tab.host.title}")
    tab.host.delete()
    a.mode_box.setCurrentIndex(1)
    spin(1500)
    print(f"[status] {a.status.text()}")


def s5_movie():
    section("5.3 movie")
    a = window("animate")
    path = os.path.join(REPO, "gb1_meld_300K_model1_cluster_mode2.gif")
    a.save_movie(path)
    print(f"[note] {a.note.text()}")
    os.remove(path)
    a.mode_box.setCurrentIndex(0)
    spin(800)
    a.delete()


def s6_morph():
    section("6. morph")
    do("""delete gb1_m1_mode*
          enable gb1_m1
          enable gb1_m2
          mdcif_morph gb1""")
    spin(500)
    do("mdcif_movie examples/view/pymol/figures/morph.gif, gb1_morph, width=480, "
       "height=416, ray=1")


def s7_plot():
    section("7.1 plot")
    do("""delete gb1_morph
          disable gb1_m2
          enable gb1_m1
          mdcif_plot distance, /gb1_m1//A/10/CA, /gb1_m1//A/40/CA""")
    spin(1000)
    from mdcif_pymol import analyze_tool
    tab = analyze_tool.windows()[-1]
    tab.atom_boxes[0].setText("12CA")
    tab.atom_boxes[1].setText("38CA")
    tab.add_typed_button.click()
    spin(800)
    tab.host.resize(820, 560)
    spin(800)
    tab.go_to_frame(8)
    spin(500)
    print(f"[title] {tab.host.title}")
    print(f"[info] {tab.info.text()}")
    save_widget(tab.host, "plot_window.png")
    S["tab"] = tab


def s7_structure():
    section("7.3 distances on the structure")
    do("""distance d1, /gb1_m1_mode1//A/10/CA, /gb1_m1_mode1//A/40/CA
          distance d2, /gb1_m1_mode1//A/12/CA, /gb1_m1_mode1//A/38/CA
          color black, d1
          color black, d2
          set label_size, 24""")
    cmd.frame(8)
    spin(300)
    scene("plot_structure.png")


def s7_angles():
    section("7.4 angles and dihedrals")
    do("mdcif_plot angle, /gb1_m1_mode1//A/10/CA, /gb1_m1_mode1//A/11/CA, "
       "/gb1_m1_mode1//A/12/CA")
    do("mdcif_plot dihedral, /gb1_m1_mode1//A/10/N, /gb1_m1_mode1//A/10/CA, "
       "/gb1_m1_mode1//A/10/C, /gb1_m1_mode1//A/11/N")
    spin(1000)
    tab = S["tab"]
    tab.host.resize(820, 760)
    spin(800)
    tab.go_to_frame(8)
    spin(500)
    print(f"[title] {tab.host.title}")
    print(f"[info] {tab.info.text()}")
    save_widget(tab.host, "plot_stacked.png")


def s7_save():
    section("7.5 save")
    do("mdcif_plot distance, /gb1_m1_mode1//A/10/CA, /gb1_m1_mode1//A/40/CA, "
       "save_image=examples/view/pymol/figures/plot_distance.png, "
       "save_data=examples/view/pymol/figures/gb1_mode1.csv")
    do("mdcif_plot dihedral, /gb1_m1_mode1//A/10/N, /gb1_m1_mode1//A/10/CA, "
       "/gb1_m1_mode1//A/10/C, /gb1_m1_mode1//A/11/N, new=1, "
       "save_image=examples/view/pymol/figures/plot_dihedral.png")
    spin(800)
    from mdcif_pymol import analyze_tool
    analyze_tool.windows()[-1].host.delete()
    with open(os.path.join(FIG, "gb1_mode1.csv")) as fh:
        print("[csv]\n" + "".join(fh.readlines()[:12]))
    os.remove(os.path.join(FIG, "gb1_mode1.csv"))


def s7_pick():
    section("7.6 pick")
    tab = S["tab"]
    tab.set_kind("distance")
    w = tab.pick_atoms()
    print(f"[prompt] {w.get_prompt()}")
    for sele in ("/gb1_m1_mode1//A/40/CA", "/gb1_m1_mode1//A/12/CA"):
        cmd.select("sele", sele)
        w.do_select("sele")
        spin(300)
    spin(500)
    print(f"[series] {[s.label for s in tab.plot.series]}")
    tab._guard(tab.add_picked)
    print(f"[note] {tab.note.text()}")


def s7_tidy():
    section("7.7 tidy")
    do("""delete d1
          delete d2""")
    S["tab"].host.delete()
    from mdcif_pymol import menu
    menu.menu_show("analyze")
    spin(800)
    from mdcif_pymol import analyze_tool
    analyze_tool.windows()[-1].host.delete()


def s8_more():
    section("8.1 residues")
    do("""delete gb1_m1_mode1
          enable gb1_m1
          disable gb1_m2
          show sticks, gb1_m1 and resi 9-13 and not name N+C+O
          color atomic, gb1_m1 and resi 9-13 and not elem C
          label gb1_m1 and resi 9-13 and name CA, "%s%s" % (resn, resi)
          set label_size, 20
          orient gb1_m1 and resi 9-13
          zoom gb1_m1 and resi 9-13, 4""")
    scene("turn.png")
    section("8.2 some residues")
    do("""label
          hide sticks
          enable gb1_m2
          hide cartoon
          show cartoon, (gb1_m1 or gb1_m2) and resi 1-20
          orient (gb1_m1 or gb1_m2) and resi 1-20""")
    scene("hairpin.png")
    section("8.3 porcupine and residues")
    do("""show cartoon
          disable gb1_m2
          orient gb1
          mdcif_porcupine gb1_m1, mode=1
          color 0xb8336a, gb1_m1 and resi 9-13""")
    scene("porcupine_turn.png")
    section("8.4 effects")
    do("""delete gb1_m1_mode1_porcupine
          spectrum resi, rainbow, gb1_m1
          show surface, gb1_m1
          set transparency, 0.7
          set ray_trace_mode, 1
          bg_color 0xf4f4f2""")
    scene("effects.png")


def s9_sessions():
    section("9. sessions")
    do("""set ray_trace_mode, 0
          hide surface
          mdcif_animate gb1_m1, mode=2
          save examples/view/pymol/figures/gb1.pse
          reinitialize
          load examples/view/pymol/figures/gb1.pse
          mdcif_porcupine gb1_m1, mode=3""")
    os.remove(os.path.join(FIG, "gb1.pse"))


STEPS = [s2_plugin, s2_menu, s3_open, s4_porcupine, s4_window, s4_window_modes,
         s4_window_mode2, s5_animate, s5_window, s5_window_steps, s5_movie, s6_morph,
         s7_plot, s7_structure, s7_angles, s7_save, s7_pick, s7_tidy, s8_more, s9_sessions]


def _run(k=0):
    if k >= len(STEPS):
        _log.write("\n######## done\n")
        _log.flush()
        QtCore.QTimer.singleShot(500, cmd.quit)
        return
    try:
        STEPS[k]()
    except Exception:
        _log.write(f"\n######## FAILED in {STEPS[k].__name__}\n")
        traceback.print_exc(file=sys.stdout)
    QtCore.QTimer.singleShot(300, lambda: _run(k + 1))


QtCore.QTimer.singleShot(1500, _run)
