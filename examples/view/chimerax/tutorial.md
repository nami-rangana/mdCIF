# Viewing an mdCIF file in ChimeraX: a step-by-step tutorial

An mdCIF file captures a simulation's main states as model structures, plus the three largest
**modes** of each state. A mode shows how the atoms move together. This tutorial looks at such
a file in [UCSF ChimeraX](https://www.cgl.ucsf.edu/chimerax/) with the mdcif plugin for
ChimeraX (the ChimeraX-MDCIF bundle in [`viewers/chimerax`](../../../viewers/chimerax)). The
plugin adds a few `mdcif` commands, and a menu, **MDCIF**, whose windows do the main ones with
sliders and buttons instead of typing:

| Command | Or, in the MDCIF menu | What it does |
|---|---|---|
| `mdcif porcupine` | Porcupine… | arrows showing how each Cα moves in one mode |
| `mdcif animate` | Animate… | a copy of a model that swings back and forth along one mode |
| `mdcif morph` | | a movie that morphs from one model to the next |
| `mdcif plot` | Analyze… | a distance, angle or dihedral at every frame of a mode, in a window you can click |
| `mdcif modes` | | lists the models and their modes in the Log |
| `mdcif open` | | another way to open the file (the same as File → Open) |

We use the GB1 example, `examples/create/1_basic/gb1_meld_300K.cif`, made by
[`1_basic.ipynb`](../../create/1_basic/1_basic.ipynb): a 300 K MELD simulation of the protein G
B1 domain, summarised as two states.

**How to follow along.** Open ChimeraX and type each block of commands at its command line (the
box at the bottom of the window), in order: each block continues from the one before. Where it
helps, the text that ChimeraX's **Log** shows comes after the commands. Every picture and movie
here was made by ChimeraX 1.12 with these same commands. The notebook
[`view_mdcif.ipynb`](view_mdcif.ipynb) runs the same steps for you, if you would rather not type
them.

**Commands or windows.** Sections 4, 5 and 7 show each step twice: first as commands, then with
the window from the MDCIF menu that does the same (sections 4.3, 5.2 and 7.2). Use whichever
you like, or both. Whatever you do in a window, it writes the matching command to the Log, so
the windows are also a good way to learn the commands.

**You'll need:** ChimeraX 1.7 or newer (we used 1.12), a copy of this repository (for the
plugin), and the example file above.

**Contents:** [1. Install ChimeraX](#1-install-chimerax) ·
[2. Add the plugin](#2-add-the-mdcif-plugin) · [3. Open the file](#3-open-the-file) ·
[4. Porcupine plots](#4-porcupine-plots) · [5. Movies of the modes](#5-bring-the-modes-to-life) ·
[6. Morph](#6-morph-between-the-models) ·
[7. Plot a distance, angle or dihedral](#7-plot-a-distance-angle-or-dihedral-along-a-mode) ·
[8. More with ChimeraX](#8-more-with-chimerax-pick-residues-change-the-look) ·
[9. Sessions](#9-sessions) · [10. Cheat sheet](#10-cheat-sheet-commands-and-windows)

## 1. Install ChimeraX

1. Go to the ChimeraX website, <https://www.cgl.ucsf.edu/chimerax/>.
2. Click **Download** and pick the installer for your computer (Windows, macOS or Linux).
   ChimeraX is free for academic and other non-commercial use; the download page asks you to
   accept its license.
3. Run the installer as you would for any other program, then start ChimeraX.

The plugin needs ChimeraX 1.7 or newer (we used 1.12). To check which one you have, type
`version` at the command line:

```
UCSF ChimeraX version: 1.12 (2026-06-12)
```

## 2. Add the mdcif plugin

Type these two commands, with the path to your copy of this repository:

```
devel clean /path/to/mdcif/viewers/chimerax
devel install /path/to/mdcif/viewers/chimerax
```

`devel clean` removes files left by an earlier build: without it, `devel install` can install
an old copy it finds there. The Log ends with

```
Installed ChimeraX-MDCIF (0.5.0)
```

and the plugin stays installed when you restart ChimeraX. On Windows the commands are the same,
with forward slashes: `devel install C:/Users/you/code/mdcif/viewers/chimerax` (not tested
yet). To check the plugin, type `mdcif`:

```
ChimeraX-MDCIF 0.5.0. Subcommands: mdcif open | modes | morph | animate | porcupine | plot. Attaching mdcif data to opened .cif files: on. Windows: menu MDCIF > Porcupine, Animate, Analyze.
```

"Attaching mdcif data to opened .cif files: on" matters: the plugin attaches mdcif data to every
`.cif` file you open, so you open an mdCIF file the normal way and the `mdcif` commands just
work.

The plugin also adds the menu **MDCIF** to ChimeraX's menu bar, just before Help. Its three
entries open the windows of sections 4.3, 5.2 and 7. Hold the pointer over an entry to see what
it does:

![The MDCIF menu, open, with the pointer on Porcupine...](figures/mdcif_menu.png)

The same windows are also in **Tools → Structure Analysis** (MDCIF Porcupine, MDCIF Animate and
MDCIF Analyze). A window opened before any mdCIF file says so ("Open an mdcif file (File >
Open) …"), with its buttons greyed out, and fills in as soon as you open one. The windows float
on top of ChimeraX's main window: drag one aside by its title bar if it covers the protein.

**To update** the plugin after its code has changed, remove it first and install it again, in
the same ChimeraX:

```
toolshed uninstall ChimeraX-MDCIF
devel clean /path/to/mdcif/viewers/chimerax
devel install /path/to/mdcif/viewers/chimerax
```

The new version works at once. `toolshed uninstall` takes the MDCIF menu away and closes its
windows, and `devel install` brings back a new menu. Don't just run `devel install` again: the
plugin is in use from the moment ChimeraX starts (it watches for mdCIF files being opened), so
ChimeraX asks to restart, and after the restart it keeps the old code if the version number
hasn't changed. If an install ends with "custom initialization failed", restart ChimeraX.

## 3. Open the file

Open the mdCIF file the way you open any structure in ChimeraX: **File → Open…** and pick the
file, drag the file onto the window, or type

```
open /path/to/mdcif/examples/create/1_basic/gb1_meld_300K.cif
```

![ChimeraX's File → Open dialog with gb1_meld_300K.cif picked](figures/file_open.png)

The Log shows ChimeraX's usual report on the file, then the plugin's summary of what the mdCIF
file holds:

```
gb1_meld_300K.cif: 2 models
  global set 1: 3 PCA modes on 56 atoms (name CA), 17683 frames, eigenvalues 203/35.5/25.6 A^2, sd 14/6/5.1 A, variance 60/11/8%
  5000 frames projected on global modes 1-3, 0% in no deposited cluster
  model 1: Cluster 1 medoid (frame 142700, t = 285.402 ns) (95% of frames)
    cluster set 2: 3 PCA modes on 56 atoms (name CA), 16847 frames, eigenvalues 37.4/28.5/10.8 A^2, sd 6.1/5.3/3.3 A, variance 26/20/8%
  model 2: Cluster 2 medoid (frame 127120, t = 254.242 ns) (5% of frames)
    cluster set 3: 3 PCA modes on 56 atoms (name CA), 836 frames, eigenvalues 7.59/5.57/3.35 A^2, sd 2.8/2.4/1.8 A, variance 20/14/9%
```

For each set of modes it gives the number of frames it was computed from, the size of each mode
(`sd`, one standard deviation of the motion, in Å) and each mode's share of the motion
(`variance`). There are two kinds of modes:

- **cluster** modes (the default in every command) are a model's own: how the protein moves
  while it stays in that state;
- **global** modes are computed over the whole simulation. The first three are the directions
  along which the frames were sorted into states.

ChimeraX also says that the file "contains IHM (integrative hybrid model) data" and suggests
opening it with `format ihm`. Ignore that: the normal open gives the `mdcif` commands everything
they need, and they don't work on a file opened as an Integrative Hybrid Model. Use the text
`.cif` file: the plugin doesn't read BinaryCIF (`.bcif`) or gzipped files.

Each model is its own structure. Colour them, model 1 blue and model 2 orange, on a white
background:

```
color #1.1 #2a78d6; color #1.2 #eb6834
set bgColor white; lighting soft; view
```

![ChimeraX's window after opening the file](figures/window_open.png)

**Figure 1.** ChimeraX's window after opening the file. The **Models** panel (bottom right)
lists the two models under the file's name: **#1.1 is model 1**, the state seen most often, and
**#1.2 is model 2**. Each is named after its state's most typical frame (its medoid) and when it
occurred in the simulation. The **Log** (top right) keeps everything the commands report, and
the command line at the bottom is where you type.

Both models are aligned to the NMR structure 3GB1, so they overlap and you can compare their
shapes: the folded model 1 and the misfolded model 2.

`mdcif open` with the file's path does the same as `open`; it also reads a file whose name
doesn't end in `.cif`. `mdcif modes` prints the summary again, which is handy after the Log has
scrolled away.

The pictures below were saved with `save picture.png width 700 supersample 3`.

## 4. Porcupine plots

### 4.1. Draw the arrows

A porcupine plot turns a mode into a snapshot: an arrow on each Cα shows which way it moves and
how far it goes at 3σ (option `scale`). Tiny movements under 0.5 Å (option `minLength`) are left
out. Modes move in both directions, so the arrows could point the other way too.

```
hide #1.2 models
mdcif porcupine #1.1 mode 1
```

```
#1.1.1: 51 arrows, cluster mode 1 of model 1, 3 sd
```

![Porcupine plot of mode 1 of model 1](figures/porcupine_1_mode_1.png)

The arrows are a model of their own inside the structure (#1.1.1 here; the Log gives the
number), so they can be hidden or closed like any model: `close #1.1.1`. The longest arrows sit
on the β1–β2 turn (residues 9–13) and on the loop after the helix.

### 4.2. Options

All options are optional and go in any order:

| Option | Meaning | Default |
|---|---|---|
| `mode 2` | which mode (1 = largest) | 1 |
| `scale 2` | arrow length, in standard deviations (σ) | 3 |
| `color forestgreen` | arrow colour: any ChimeraX colour, also `#ff000080` with transparency | orange |
| `scope global` | the global modes instead of the model's own | cluster |
| `minLength 1` | leave out arrows shorter than this (Å; 0 keeps them all) | 0.5 |
| `radius 0.3` | arrow thickness (Å) | 0.15 |

Mode 2 of model 1, with shorter (2σ) green arrows:

```
close #1.1.1
mdcif porcupine #1.1 mode 2 scale 2 color forestgreen
```

```
#1.1.1: 34 arrows, cluster mode 2 of model 1, 2 sd
```

![Porcupine plot of mode 2 of model 1, green arrows at 2 sd](figures/porcupine_1_mode_2.png)

And global mode 1, the main direction of the whole simulation (the one along which the frames
were sorted into states), drawn in purple on the misfolded model 2:

```
close #1.1.1; hide #1.1 models; show #1.2 models
mdcif porcupine #1.2 mode 1 scope global color #6a3d9a
```

```
#1.2.1: 55 arrows, global mode 1 of model 2, 3 sd
```

![Porcupine plot of global mode 1 on model 2, purple arrows](figures/porcupine_2_global_1.png)

### 4.3. With the window: MDCIF → Porcupine…

The **MDCIF Porcupine** window draws the same arrows with sliders, and can show several modes
at once. Put model 1 back on screen, in light grey this time so that the blue arrows of mode 2
stand out:

```
close #1.2.1; hide #1.2 models; show #1.1 models
color #1.1 #c8c8c8
```

Then choose **MDCIF → Porcupine…** in the menu bar. The window opens on model 1 (#1.1) with
mode 1 ticked, and draws its arrows at once: the same arrows as `mdcif porcupine #1.1 mode 1`,
as model #1.1.1. The Log shows what happened, as commands:

```
> ui tool show "MDCIF Porcupine"
> mdcif porcupine #1.1 mode 1 scale 3 color #ff8c00 radius 0.15 minLength 0.5 scope cluster
```

Now tick modes **2** and **3** too. Each mode gets its own arrows, in the colour of the button
at the end of its row:

![The MDCIF Porcupine window with modes 1-3 ticked, and their arrows on model 1](figures/window_porcupine.png)

**Figure 2.** Modes 1, 2 and 3 of model 1 at once: orange, blue and green arrows, each set its
own model (#1.1.1, #1.1.2 and #1.1.3). Each row gives the mode's share of the motion and its
size (one standard deviation, sd). The line at the bottom says how many arrows each mode has.

Each control of the window is an option of `mdcif porcupine` (section 4.2):

| In the window | Option of `mdcif porcupine` |
|---|---|
| **Structure** | the model to draw on (#1.1, #1.2, or a playing copy from section 5) |
| **Scope**: cluster or global | `scope` ("global" is greyed out for a file without global modes) |
| **Modes**: the tick boxes | `mode`: each ticked mode gets its own arrows |
| the colour button of a mode | `color`: click it for ChimeraX's colour chooser |
| **Length** (0.5 to 10 sd) | `scale` |
| **Thickness** (0.05 to 0.6 Å) | `radius` |
| **Min length** (0 to 3 Å) | `minLength` |

To draw what section 4.2 drew with a command, mode 2 alone with green arrows at 2σ:

1. Untick modes **1** and **3**.
2. Type `2` in the box beside **Length** and press Enter (or drag the slider to 2.0).
3. Click the colour button of mode 2. In the colour chooser, pick a green, or type `#228b22`
   in its **HTML** box and press Return.

The arrows change as you go, and they stay the same model (#1.1.2), so nothing piles up. The
line at the bottom then reads `#1.1.2 mode 2: 34 arrows (22 shorter than 0.5 Å not drawn)`: the
34 arrows of section 4.2. Once you stop changing things (let go of the slider, press Enter),
the Log gets one command per ticked mode, not one for every step of the slider. Here the last
one is

```
> mdcif porcupine #1.1 mode 2 scale 2 color #228b22 radius 0.15 minLength 0.5 scope cluster
```

What else the window does:

- **Live** is on at the start: the arrows follow the sliders while you drag. With Live off,
  nothing changes until you click **Draw**.
- **Keep** leaves the arrows you see in the scene; the next change draws new ones next to them
  (#1.1.4, …), for example to compare two thicknesses. **Clear** removes the window's arrows
  that were not kept.
- A mode with no arrow left at your settings (try Min length 3 with Length 0.5) is named in red
  at the bottom, for example `mode 2: no arrows to draw for cluster mode 2 at 0.5 sd with
  minLength 3 A`.
- **Animate this mode** opens the Animate window of section 5.2 on the same model and mode.
- Closing the window (the × in its corner) leaves its arrows in the scene.

Before the next section, click **Clear**, close the window, and colour model 1 blue again:

```
color #1.1 #2a78d6
```

## 5. Bring the modes to life

Each mode gives every Cα a direction and a typical displacement, measured in Å. `mdcif animate`
turns that into a movie: it makes a copy of the model (a new model, #2 here) that swings back
and forth up to 3σ (option `amplitude`) and keeps playing. The original is hidden. Other atoms
follow their residue's Cα.

| Option | Meaning | Default |
|---|---|---|
| `mode 2` | which mode | 1 |
| `amplitude 2` | how far to swing, in standard deviations (σ) | 3 |
| `frames 10` | frames per back-and-forth cycle (at least 3) | 30 |
| `scope global` | the global modes | cluster |

These modes follow straight paths, so large motions can stretch bonds or loops. That's a
limitation of the approximation, not of the file: real protein motions follow more complex
paths.

### 5.1. With commands

Mode 1 of the folded model 1:

```
close #1.2.1; hide #1.2 models; show #1.1 models
mdcif animate #1.1 mode 1
```

```
#2: cluster mode 1 of model 1, 26% variance, +/-3 sd ('coordset stop #2' to stop)
```

![Model 1 swinging along its mode 1](figures/model_1_mode_1.gif)

**Movie 1.** Model 1 swinging along its mode 1, ±3σ (two cycles; also as
[mp4](figures/model_1_mode_1.mp4)). The β1–β2 turn and the loop after the helix move the most,
as the porcupine plot showed.

The copy keeps playing until you stop it; the Log says how (`coordset stop #2`). To hold one
frame for a picture, stop first and then pick the frame: `coordset #2 8` is the +3σ turning
point with the default 30 frames. The copy keeps its model's modes, so `mdcif porcupine #2`
works on it too.

The misfolded model 2 along its own mode 1:

```
coordset stop #2; close #2; hide #1.1 models; show #1.2 models
mdcif animate #1.2 mode 1
```

```
#2: cluster mode 1 of model 2, 20% variance, +/-3 sd ('coordset stop #2' to stop)
```

![Model 2 swinging along its mode 1](figures/model_2_mode_1.gif)

**Movie 2.** Model 2 along its mode 1, ±3σ ([mp4](figures/model_2_mode_1.mp4)). This state was
seen in only 5% of the frames, and its modes are smaller (σ = 2.8 Å for mode 1, against 6.1 Å
in model 1).

### 5.2. With the window: MDCIF → Animate…

The **MDCIF Animate** window plays a mode with sliders, and changes it while it plays. Stop the
copy of model 2 and put model 1 back:

```
coordset stop #2; close #2; hide #1.2 models; show #1.1 models
```

Choose **MDCIF → Animate…**. The window opens on model 1 (#1.1), mode 1, and nothing moves
yet: the line at the bottom says `Play (or a change) animates a copy of #1.1, as 'mdcif
animate' does.` Click **Play**. The window makes a playing copy (#2), exactly as `mdcif animate`
does, and the Log shows that command:

```
> mdcif animate #1.1 mode 1 amplitude 3 frames 30 scope cluster
#2: cluster mode 1 of model 1, 26% variance, +/-3 sd ('coordset stop #2' to stop)
```

![The MDCIF Animate window while #2 plays](figures/window_animate.png)

**Figure 3.** The MDCIF Animate window while copy #2 plays. The line at the bottom says what
plays, and the **Frame** slider follows the motion.

The controls are the options of `mdcif animate` (section 5.1), plus a speed:

| In the window | Option of `mdcif animate` |
|---|---|
| **Structure** | the model to play |
| **Scope**: cluster or global, and **Mode** | `scope`, `mode` (the mode's share of the motion is shown next to it) |
| **Amplitude** (−10 to 10 sd) | `amplitude` (a negative number swings the other way first) |
| **Frames** (type a number, press Enter) | `frames` |
| **Speed** (1/1 = full speed, down to 1/10) | none: ChimeraX's `coordset … pauseFrames`, each frame shown 2, 3, … times |

Try these, in order:

1. Type `2` in the box beside **Speed** and press Enter: the copy now plays at half speed, and
   the Log shows the command for that, `coordset #2 1,30 loop 1000 pauseFrames 2`.
2. Drag the **Frame** slider to 8. The motion stops on that frame, the +3σ turning point, and
   the line at the bottom ends `(stopped at frame 8)`.
3. Click **Play** to start again. Play always starts at frame 1; **Stop** stops the motion.
4. Click **Save movie…** to record one back-and-forth cycle as a movie (section 5.3).
5. Click **Arrows for this mode**. The Porcupine window of section 4.3 opens on the copy with
   this mode ticked, and draws its arrows (#2.1). Click **Clear** there and close that window.
6. Click **Plot geometry…**. The Analyze window of section 7 opens for this copy, ready for the
   atoms to measure. Close it for now.
7. Choose **2** in the **Mode** list. The copy changes in place, and #2 stays #2. Once you stop
   changing things, the Log says what the copy now plays, and which command makes a copy like it:

   ```
   #2 now plays cluster mode 2 of model 1, 20% variance, +/-3 sd, 30 frames (changed in place; a new copy like it: mdcif animate #1.1 mode 2 amplitude 3 frames 30 scope cluster)
   ```

   Amplitude, Frames and Scope change the copy in place the same way. Choosing another
   Structure makes a new copy and closes the old one.

Close the Animate window with the × in its corner: the copy keeps playing.

### 5.3. Record a movie

In the Animate window, **Save movie…** asks for a file name. The suggested name says what
plays, for example `gb1_meld_300K_model1_cluster_mode1.mp4` in ChimeraX's current folder, and
you can choose MP4, WebM, QuickTime (`.mov`) or AVI. It records one cycle at the speed you set:
the 30 frames at half speed give a movie of 60 frames (2.4 s at 25 frames per second). The window
is greyed out while it records. Then the line at the bottom says `Movie saved to …`, and the
motion plays on. The Log shows the commands it ran:

```
> movie record
> coordset #2 1,30 pauseFrames 2
> wait 60
> movie encode /path/to/gb1_meld_300K_model1_cluster_mode1.mp4 framerate 25 quality good wait true
```

With commands, close the window's copy first (or `mdcif animate` makes a second one). Then start
recording, start the motion, let it play for as many frames as you want (60 = two cycles at
full speed), and write the file:

```
close #2
movie record
mdcif animate #1.1 mode 1
wait 60
movie encode ~/model_1_mode_1.mp4
```

## 6. Morph between the models

`mdcif morph` uses ChimeraX's own `morph` command to make a new model (#2 here) that changes
smoothly from model 1 to model 2 and back to model 1: 20 frames per step, 41 in all. It plays
once.

```
coordset stop #2; close #2; show #1.1,2 models
mdcif morph #1
```

```
> morph #1.1 #1.2 frames 20 wrap true play true
Computed 41 frame morph #2
> coordset #2 1,41
```

![Morph from model 1 to model 2 and back](figures/morph.gif)

**Movie 3.** The morph from the folded model 1 to the misfolded model 2 and back
([mp4](figures/morph.mp4)). A morph is a straight interpolation between the two shapes: it shows
how they differ, not the path the simulation took.

To play it again, type `coordset #2`; to repeat it a few times, `coordset #2 1,41 loop 10`.
Option `frames 10` sets the frames per step, and `wrap false` stops at the last model. ChimeraX
also shows a slider to step through the frames (bottom right):

![ChimeraX's window after the morph, with the frame slider](figures/window_morph.png)

## 7. Plot a distance, angle or dihedral along a mode

A movie shows a mode, and a plot puts numbers on it. `mdcif plot` plays a mode, the same way
`mdcif animate` does, and measures a distance, an angle or a dihedral at every frame. The result
opens in the **MDCIF Analyze** window, on its **Geometry** tab, with one point per frame.
Clicking a point stops the motion at that frame.

This section starts with the command, because it opens the window already filled in. Without
typing, the same window comes from **MDCIF → Analyze…**, or from **Plot geometry…** in the
Animate window (section 5.2); you then type or pick the atoms in it (sections 7.2 and 7.6).

### 7.1. A distance

Show model 1 on its own and measure the distance between the Cα atoms of Lys 10 (in the β1–β2
turn) and Asp 40 (in the loop after the helix). These are the two parts that move most in
mode 1.

```
coordset stop #2; close #2; hide #1.2 models; show #1.1 models
mdcif plot distance #1.1/A:10@CA #1.1/A:40@CA
```

```
> mdcif animate #1.1 mode 1 amplitude 3 frames 30 scope cluster
#2: cluster mode 1 of model 1, 26% variance, +/-3 sd ('coordset stop #2' to stop)
#2 distance /A:10@CA - /A:40@CA along cluster mode 1 of model 1, +/-3 sd, 30 frames: min 6.64, max 13.99, mean 9.97 A
```

The command made the same copy as `mdcif animate #1.1 mode 1` (#2, playing). It then measured
the distance in each of the copy's 30 frames. As the mode swings, the turn and the loop move
between 6.6 and 14.0 Å apart.

A second distance along the same mode goes into the same window, in a new colour, and uses the
same copy. Type it as a command:

```
mdcif plot distance #1.1/A:12@CA #1.1/A:38@CA
```

or type the two atoms into the window instead: `12CA` in the first **Atoms** box, `38CA` in the
second, then click **Add** (or press Return). Either way the Log shows

```
#2 distance /A:12@CA - /A:38@CA along cluster mode 1 of model 1, +/-3 sd, 30 frames: min 10.77, max 16.54, mean 12.93 A
```

### 7.2. The MDCIF Analyze window

![The MDCIF Analyze window with two distances, stopped at frame 8](figures/plot_window.png)

**Figure 4.** The MDCIF Analyze window with the two distances, stopped at frame 8 by a click on
the plot. Its title says what it plots, `MDCIF Analyze: distance on #2`. The plot is on the
**Geometry** tab, the only tab so far. Across the bottom are the frame numbers, and up the side
the distance. The second row holds the atoms typed for the second distance (`12CA`, `38CA`).

- **Type atoms** in the second row. Choose what to measure under **Measure** (Distance, Angle
  or Dihedral). The window then shows 2, 3 or 4 **Atoms** boxes. Type a residue number and an
  atom name in each box, for example `38CA`. `38 CA`, `38@CA` and `/A:38@CA` work too, and lower
  case is fine. Then click **Add**, or press Return in any box. The atoms are looked up on the
  copy the window shows. The boxes keep what you typed, so for a neighbouring residue you only
  change a number. If an atom isn't there, a red message under the list says why and nothing
  is added, for example `atom 1, on #2: no atom /A:38@CB (residue 38 has N, CA, C, O)`
  (Gly 38 has no side chain). The message lists the atoms the residue does have.
- While the copy plays, the large dots and the thin vertical line move to the frame on screen.
  The line at the bottom gives the values there, for example `frame 8 / 30: 13.99, 16.54 Å`.
  The motion plays as fast with the window open as without it.
- **Click a point** to stop the motion at that frame. Frame 8 is the +3σ turning point, where
  both distances are largest. You can also hold the button down and drag along the plot to step
  through the frames. **Play** always starts the motion again from frame 1 (not from the frame
  shown), and **Stop** stops it.
- The top row changes the motion: **Model**, **Scope** (cluster or global), **Mode**,
  **Amplitude** (in σ) and **Frames**. After a change the window plays the new motion and
  measures every line again. The Log shows the `mdcif animate` command, and the copy the window
  made before is closed. A typed number counts when you press Enter.
- **Remove** removes the lines you choose in the list, and **Remove all** removes all of them.
  A panel goes away with its last line.

### 7.3. See the distances on the structure

ChimeraX's own `distance` command draws the two distances on the copy. Their labels follow the
motion while it plays.

```
distance #2/A:10@CA #2/A:40@CA; distance #2/A:12@CA #2/A:38@CA
distance style color black decimalPlaces 1; label height 1.5
```

![The copy at frame 8 with the two distances drawn](figures/plot_structure.png)

At frame 8 the labels read 14.0 and 16.5 Å, the same values the Analyze window gives. We checked
all 30 frames of the plot against ChimeraX's `distance`, `angle` and `torsion` commands, and
they agree.

### 7.4. Angles and dihedrals

An angle takes 3 atoms and a dihedral takes 4, in order. Atoms on the copy (#2) use its frames,
so no new copy is made. They go into the same window as the distances: the angle gets a panel
under the distances, and the dihedral a panel under that. All panels share the frame axis, so
you can read every value at the same frame.

```
mdcif plot angle #2/A:10@CA #2/A:11@CA #2/A:12@CA
mdcif plot dihedral #2/A:10@N #2/A:10@CA #2/A:10@C #2/A:11@N
```

```
#2 angle /A:10@CA - /A:11@CA - /A:12@CA along cluster mode 1 of model 1, +/-3 sd, 30 frames: min 119.3, max 164.2, mean 138.7 deg
#2 dihedral /A:10@N - /A:10@CA - /A:10@C - /A:11@N along cluster mode 1 of model 1, +/-3 sd, 30 frames: min 28.1, max 147.0, mean 94.3 deg
```

![The MDCIF Analyze window with the distances, the angle and the dihedral, stopped at frame 8](figures/plot_stacked.png)

**Figure 5.** The same window after the two commands, stopped at frame 8: the distances at the
top, the angle in the middle, the dihedral at the bottom. The window is now called `MDCIF
Analyze: distance, angle, dihedral on #2`. The list names the kind of each line, and the line at the
bottom gives the values at the shown frame, kind by kind: `frame 8 / 30: 13.99, 16.54 Å ·
119.3° · 28.1°`. A click or drag on any panel stops the motion at that frame. The large dots
and the vertical line then move in every panel.

In the window, you can do the same by choosing **Angle** under **Measure** and typing `10CA`,
`11CA` and `12CA` in the three boxes, or **Dihedral** and `10N`, `10CA`, `10C` and `11N` in four.
You can switch Measure at any time, and the lines already plotted stay. Up to 8 lines fit in
each panel. To put an angle or a dihedral in a separate plot instead, add `new true` to the
command (it opens a second window, see section 7.5).

The angle is measured at the middle atom, Thr 11. The dihedral is ψ of Lys 10, with the same
sign as ChimeraX's `torsion` command gives. A dihedral starts between −180° and 180° at frame 1,
and after that it changes smoothly, so it can go past ±180°. For example, ψ of Thr 11 reaches
256.4°, which is the same angle as −103.6°. A dihedral's CSV file (section 7.5) says so in an
extra line: `# dihedrals unwrapped from frame 1: values can pass +/-180 (subtract or add 360 for
ChimeraX's torsion, -180 to 180)`.

The same atoms in reverse order are the same measurement, so
`mdcif plot distance #2/A:40@CA #2/A:10@CA` only logs `... is already plotted`.

### 7.5. Save the plot and the numbers

In a window, **Save plot…** saves the picture (PNG, SVG or PDF) and **Save data…** saves the
numbers (a CSV file). Both save the whole plot, every panel and every line. The suggested names
say what is in the file, here `gb1_meld_300K_model1_cluster_mode1_distance-angle-dihedral.png`
and `.csv`. You can also save from the command line, which works without the window too:

```
mdcif plot distance #2/A:10@CA #2/A:40@CA saveImage gb1_mode1.png saveData gb1_mode1.csv
```

```
#2 distance /A:10@CA - /A:40@CA is already plotted
#2 distance + angle + dihedral plot: data saved to gb1_mode1.csv
#2 distance + angle + dihedral plot: image saved to gb1_mode1.png
```

This distance is already plotted, so nothing is added: the command saves its whole plot, here
the two distances, the angle and the dihedral. The files go to ChimeraX's current folder, which
`pwd` shows (`cd` changes it), or give a full path. The picture is 6 inches wide and 300 dpi,
without the current-frame marks. It is 3.5 inches high for one panel and 2.2 inches more for
each further panel, so 7.9 inches here (1800 × 2370 pixels):

![The saved plot: the two distances, the angle and the dihedral, one panel each](figures/plot_distance.png)

For a picture of one measurement on its own, start a separate plot with `new true`. It opens a
second window, which you can close afterwards:

```
mdcif plot dihedral #2/A:10@N #2/A:10@CA #2/A:10@C #2/A:11@N new true saveImage gb1_psi10.png
```

```
#2 dihedral /A:10@N - /A:10@CA - /A:10@C - /A:11@N along cluster mode 1 of model 1, +/-3 sd, 30 frames: min 28.1, max 147.0, mean 94.3 deg
#2 dihedral plot: image saved to gb1_psi10.png
```

![The saved plot of the dihedral psi of Lys 10 alone](figures/plot_dihedral.png)

The CSV file starts with a few lines beginning with `#` that say what was measured. The `kind`
line lists the kinds in the plot with their units, and the extra dihedral line explains the
values past ±180° (section 7.4). Then comes one column per line of the plot, named by its atoms,
in the order the lines were added. Two atoms in a name mean a distance, three an angle and four
a dihedral. In Python, `pandas.read_csv("gb1_mode1.csv", comment="#")` reads it.

```
# mdcif plot (ChimeraX-MDCIF 0.5.0)
# file: /path/to/mdcif/examples/create/1_basic/gb1_meld_300K.cif
# model: 1 (Cluster 1 medoid (frame 142700, t = 285.402 ns))
# scope: cluster
# mode: 1
# amplitude: 3 sd
# frames: 30
# kind: distance (A), angle (deg), dihedral (deg)
# dihedrals unwrapped from frame 1: values can pass +/-180 (subtract or add 360 for ChimeraX's torsion, -180 to 180)
frame,/A:10@CA - /A:40@CA,/A:12@CA - /A:38@CA,/A:10@CA - /A:11@CA - /A:12@CA,/A:10@N - /A:10@CA - /A:10@C - /A:11@N
1,9.64629,11.1256,141.271,106.734
2,10.4915,11.8198,134.21,87.1693
...
```

### 7.6. Pick the atoms with the mouse

Atoms are easier to pick when they hold still, so if the copy is playing, click a point in the
plot to stop it. A click on the cartoon selects a whole residue, so also hide the cartoon and
show the Cα atoms:

```
cartoon hide #2; show #2@CA atoms
```

Then pick the atoms in order. **Ctrl-click** the Cα of Asp 40, then **Ctrl-Shift-click** the Cα
of Leu 12. In the Analyze window, with **Distance** chosen under **Measure**, click **Add selected
atoms**: the line `/A:40@CA - /A:12@CA` is added to the distance panel. The window uses the
order in which you picked the atoms. Typing `mdcif plot distance sel` does the same, with the
atoms in picking order.

If the wrong number of atoms is selected, a red message appears in the window (and in the Log),
and nothing is added. If you clicked the cartoon instead of the atoms, the message says so and
gives the commands above:

```
a distance needs 2 selected atoms, 12 are selected (used in picking order) - a click on the cartoon selects whole residues: show the atoms first (cartoon hide #2; show #2@CA atoms)
```

**Start from an empty window.** When no Analyze window is open, **MDCIF → Analyze…** (also
Tools → Structure Analysis → MDCIF Analyze) opens an empty one. When Analyze windows are open,
it brings them to the front instead, an empty one on top, else the newest. **Plot geometry…** in
the Animate window (section 5.2) also opens an empty one, for the copy it plays. The empty
window says what to do, with the number of the model on screen, here the copy #2:

```
Type 2 atoms (e.g. 10CA 40CA) and click Add, or
show the atoms (cartoon hide #2; show #2@CA atoms),
select 2 in order (Ctrl-click, then Ctrl-Shift-click)
and click Add selected atoms
```

Choose Distance, Angle or Dihedral under **Measure**, then type the atoms or add selected atoms.
Typed atoms are looked up on the first mdcif model on screen. If no mdcif model is shown, the
window says so in red. Atoms on a model (#1.1) are measured on a copy of it that is already on
screen, if there is one, else played with the settings in the window's top row. Atoms on a copy
use that copy.

### 7.7. Tidy up

```
select clear; cartoon #2; hide #2 atoms; ~distance
```

Close the Analyze windows with the × in their corner. Closing a window leaves its copy playing, and
a later `mdcif plot` on model 1 uses that copy again instead of making a second one. If you
close the copy first, the window stays: it says the motion was closed, and Save plot and Save
data still work.

## 8. More with ChimeraX: pick residues, change the look

Everything ChimeraX can do works on these models. Two ideas cover most of it:

- An **atom specification** picks atoms: `#1.1` is model 1, `#1.1:9-13` its residues 9 to 13,
  `#1.1,2:12,39` residues 12 and 39 of both models, `:leu` every leucine, `@ca` the Cα atoms.
  Combine them with `&` (and), `|` (or) and `~` (not): `#1.1:9-13 & sidechain`. Residue
  numbers are the ones in the file.
- **Commands** say how to draw them: `cartoon` and `cartoon hide`, `show … atoms` and
  `hide … atoms`, `style … stick` (or `sphere`, `ball`), `color`, `label`, `surface`,
  `transparency`, `lighting`, and `view` to point the camera at them.

### 8.1. Highlight residues

The β1–β2 turn (residues 9–13) moves the most in all three modes of the folded state. Here it is
in model 1, with its side chains as sticks and every residue labelled (`pad` leaves some of the
protein around it in view):

```
close #2; show #1.1 models; hide #1.2 models
show #1.1:9-13 atoms; style #1.1:9-13 stick; color #1.1:9-13 byhetero
label #1.1:9-13 residues height 0.8; view #1.1:9-13 pad 0.4
```

![The β1–β2 turn of model 1 as labelled sticks](figures/turn.png)

### 8.2. Show only some residues

To draw only part of the protein, hide all cartoons and show the residues you want. Here are
residues 1–20 (β1, the turn and β2): the hairpin of the folded model 1 (blue) next to the same
residues in the misfolded model 2 (orange).

```
label delete; hide #1.1:9-13 atoms; show #1.2 models
cartoon hide #1; cartoon #1.1,2:1-20; view #1.1,2:1-20
```

![Residues 1-20 of both models](figures/hairpin.png)

### 8.3. Residues in porcupine plots and movies

The arrows and the animated copies are ordinary models, so colours and selections work on them
the usual way. Here the turn is coloured pink in the porcupine plot of mode 1, which shows that
it carries the longest arrows. Hiding residues works less well in porcupine plots: the arrows
are not atoms, so they stay where they are (`close` them, or draw them again with `minLength`).

```
cartoon #1; hide #1.2 models; view
mdcif porcupine #1.1 mode 1
color #1.1:9-13 #b8336a target c
```

![Porcupine plot of mode 1 with the β1–β2 turn in pink](figures/porcupine_turn.png)

### 8.4. Visual effects

A few more of ChimeraX's options: a rainbow along the chain (`rainbow`, from the N-terminus in
blue to the C-terminus in red), a see-through molecular surface, outlines
(`graphics silhouettes`), shadows (`lighting full`) and a background colour:

```
close #1.1.1; rainbow #1.1
surface #1.1; transparency #1.1 70 target s
graphics silhouettes true; lighting full; set bgColor #f4f4f2; view
```

![Model 1 in rainbow colours with a see-through surface, outlines and shadows](figures/effects.png)

`turn y 1 360` spins the model once around the vertical axis.

## 9. Sessions

Save your work as a ChimeraX session, with `save` or File → Save…. When you open the session
again, the modes come back too, also for animated copies, as long as the mdCIF `.cif` file is
still where it was:

```
save ~/gb1.cxs
close session
open ~/gb1.cxs
mdcif porcupine #1.1 mode 3
```

```
mdcif: data of gb1_meld_300K.cif attached to #1.1 #1.2
opened ChimeraX session
#1.1.2: 39 arrows, cluster mode 3 of model 1, 3 sd
```

(The arrows are #1.1.2 this time because the surface from 8.4 is #1.1.1.) If the `.cif` file
has moved, the `mdcif` commands say "its file … is gone": open the `.cif` file again.

Porcupine arrows come back as they were, in their colours, whether you drew them with
`mdcif porcupine` or with the Porcupine window.

**The MDCIF windows (Porcupine, Animate and Analyze, with its plots) are not saved in a
session.** `close session` closes them, and opening the session brings none back: open them
again from the MDCIF menu. An animated copy is saved, and it remembers its mode. After you open
the session again, type the plot command again, for example
`mdcif plot distance #2/A:10@CA #2/A:40@CA` with the copy's atoms, or the same atoms on the
model (`#1.1/…`). Either way it uses the restored copy's frames and makes no new copy. (The
Animate window doesn't take over a restored copy: its Play makes a new one.)

## 10. Cheat sheet: commands and windows

| To… | type | or, in the MDCIF menu |
|---|---|---|
| open an mdCIF file | `open file.cif` (or File → Open…, or drag and drop); `mdcif open file.cif` | |
| list the models and their modes | `mdcif modes` | (each window lists the modes of its model) |
| draw the arrows of mode 2 of model 1 | `mdcif porcupine #1.1 mode 2` (options `scale`, `color`, `scope global`, `minLength`, `radius`) | Porcupine…: Structure #1.1, tick 2; sliders Length, Thickness, Min length |
| draw several modes at once | one `mdcif porcupine` per mode | Porcupine…: tick them all |
| remove those arrows | `close #1.1.1` (the Log gives the number) | Porcupine…: Clear |
| play mode 1 of model 2 | `mdcif animate #1.2 mode 1` (options `amplitude`, `frames`, `scope global`) | Animate…: Structure #1.2, Mode 1, Play |
| change what plays | a new `mdcif animate` (a new copy) | Animate…: Mode, Amplitude, Frames (the same copy, changed in place) |
| stop it, or hold one frame | `coordset stop #2`, then `coordset #2 8` | Animate…: Stop; drag the Frame slider |
| play it slower | `coordset #2 1,30 loop 1000 pauseFrames 2` | Animate…: Speed |
| morph through the models | `mdcif morph #1` (options `frames`, `wrap false`); again: `coordset #2` | |
| plot a distance along mode 1 of model 1 | `mdcif plot distance #1.1/A:10@CA #1.1/A:40@CA` (options `mode`, `scope global`, `amplitude`, `frames`, `new true`) | Analyze… (or Plot geometry… in Animate): type `10CA` `40CA`, click Add |
| plot an angle or a dihedral | `mdcif plot angle` with 3 atoms, `mdcif plot dihedral` with 4, in order (a panel under the distances; `new true` for a plot of its own); `sel` = the atoms you picked, in picking order | Analyze…: choose Measure, type 3 or 4 atoms, Add (or Add selected atoms) |
| stop the motion at one frame of the plot | | click that point in any panel of the Analyze window |
| save a plot and its numbers | `mdcif plot … saveImage plot.png saveData plot.csv` | Analyze…: Save plot… / Save data… (every panel) |
| open an empty Analyze window | | Analyze… (or Tools → Structure Analysis → MDCIF Analyze) |
| every option of a command | `usage mdcif porcupine` | |
| install or update the plugin | `toolshed uninstall ChimeraX-MDCIF` (if installed), `devel clean PATH`, `devel install PATH` | |
| colour a model (or only its cartoon) | `color #1.1 cornflowerblue` (`… target c`) | |
| show some residues as sticks | `show #1.1:9-13 atoms; style #1.1:9-13 stick` | |
| save a picture | `save picture.png width 1200 supersample 3` | |
| record a movie | `movie record`, then `mdcif animate …; wait 60`, then `movie encode movie.mp4` | Animate…: Save movie… (one cycle) |
| save everything | `save mywork.cxs` (the MDCIF windows are not saved) | |

**See more.** The ChimeraX User Guide explains every command, with all its options:

- [open](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/open.html),
  [save](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/save.html) and
  [close](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/close.html)
- [atom specifications](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/atomspec.html),
  [select](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/select.html) and
  [selecting with the mouse](https://www.cgl.ucsf.edu/chimerax/docs/user/selection.html)
  (Ctrl-click, Ctrl-Shift-click)
- [distance](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/distance.html),
  [angle](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/angle.html) and
  [torsion](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/torsion.html), ChimeraX's own
  measurements; [pwd](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/pwd.html) and
  [cd](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/cd.html) for the folder files are
  saved in
- [color](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/color.html) (also `rainbow`),
  [cartoon](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/cartoon.html),
  [show and hide](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/show.html),
  [style](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/style.html),
  [label](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/label.html)
- [surface](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/surface.html),
  [transparency](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/transparency.html),
  [lighting](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/lighting.html),
  [graphics](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/graphics.html),
  [view](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/view.html),
  [turn](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/turn.html)
- [coordset](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/coordset.html),
  [morph](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/morph.html),
  [movie](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/movie.html) and
  [wait](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/wait.html)
- [devel](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/devel.html) and
  [toolshed](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/toolshed.html), for
  installing the plugin
- the [Log](https://www.cgl.ucsf.edu/chimerax/docs/user/tools/log.html) and the
  [Models panel](https://www.cgl.ucsf.edu/chimerax/docs/user/tools/modelpanel.html)

