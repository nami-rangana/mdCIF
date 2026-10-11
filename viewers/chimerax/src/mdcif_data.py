"""mdcif file contents for viewers: pure Python + numpy, no mdcif / gemmi / host imports.

Identical copies live in viewers/pymol/mdcif_pymol/ and viewers/chimerax/src/ (the
viewers run in their host's own Python and must stay self-contained); keep them in
sync - tests/test_viewers.py checks that they are identical.

Reads text mmCIF only (not BinaryCIF). Atoms are matched between the file and the
viewer by (auth chain ID, auth residue number, atom name), as strings.

Display conventions (docs/spec/mdcif_spec.md): displacement = k * sqrt(lambda) * v,
k in standard deviations. Atoms without their own vector component (side chains
when the modes are on CA atoms) move with the vector of their residue's mode atom.
"""

import numpy as np

CATEGORIES = ("_atom_site", "_ihm_model_list", "_mdcif_mode_set", "_mdcif_mode",
              "_mdcif_mode_vector", "_mdcif_mode_origin", "_mdcif_projection")


# --------------------------------------------------------------------------- CIF reading

def _tokens(text):
    """Yield ``(value, quoted)`` for every CIF token (semicolon text fields included)."""
    lines = text.splitlines()
    i, n_lines = 0, len(lines)
    while i < n_lines:
        line = lines[i]
        if line.startswith(";"):
            buf = [line[1:]]
            i += 1
            while i < n_lines and not lines[i].startswith(";"):
                buf.append(lines[i])
                i += 1
            i += 1
            yield "\n".join(buf).strip(), True
            continue
        pos, n = 0, len(line)
        while pos < n:
            c = line[pos]
            if c in " \t":
                pos += 1
            elif c == "#":
                break
            elif c in "'\"":
                end = pos + 1
                while True:
                    end = line.find(c, end)
                    if end < 0 or end + 1 >= n or line[end + 1] in " \t":
                        break
                    end += 1
                end = n if end < 0 else end
                yield line[pos + 1:end], True
                pos = end + 1
            else:
                end = pos
                while end < n and line[end] not in " \t":
                    end += 1
                yield line[pos:end], False
                pos = end
        i += 1


def read_tables(path, categories=CATEGORIES):
    """``{category: {item: [values]}}`` for the wanted categories of a text mmCIF file.

    Unquoted ``.`` and ``?`` become ``None``.
    """
    with open(path, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    if text[:1] not in ("d", "#", " ", "\n", "\r", "\t"):
        raise ValueError(f"{path}: not a text mmCIF file (BinaryCIF is not supported here)")
    wanted = set(categories)
    tables = {}
    toks = _tokens(text)

    def value(tok):
        v, quoted = tok
        return None if not quoted and v in (".", "?") else v

    pending = next(toks, None)
    while pending is not None:
        tok, quoted = pending
        pending = None
        if not quoted and tok.lower() == "loop_":
            names = []
            pending = next(toks, None)
            while pending is not None and not pending[1] and pending[0].startswith("_"):
                names.append(pending[0])
                pending = next(toks, None)
            cat = names[0].split(".")[0] if names else ""
            keep = cat in wanted
            items = [nm.split(".", 1)[1] for nm in names]
            cols = [[] for _ in items]
            k = 0
            while pending is not None and (pending[1] or not (
                    pending[0].startswith("_") or pending[0].lower() in ("loop_",)
                    or pending[0].lower().startswith(("data_", "save_")))):
                if keep:
                    cols[k % len(items)].append(value(pending))
                k += 1
                pending = next(toks, None)
            if keep:
                tables.setdefault(cat, {}).update(zip(items, cols, strict=True))
        elif not quoted and tok.startswith("_"):
            val = next(toks, None)
            cat, _, item = tok.partition(".")
            if cat in wanted and val is not None:
                tables.setdefault(cat, {})[item] = [value(val)]
        pending = pending or next(toks, None)
    return tables


# --------------------------------------------------------------------------- file model

def _floats(col):
    return np.array([np.nan if v is None else float(v) for v in col])


class ModeSet:
    """One mode set: eigenvalues, variance fractions and per-atom vectors."""

    def __init__(self, info):
        self.id = int(info["id"])
        self.scope = info.get("scope") or "cluster"
        self.model = int(info["model_id"])
        self.method = info.get("method") or "PCA"
        self.selection = info.get("atom_selection") or ""
        self.n_frames = int(info["num_frames"]) if info.get("num_frames") else None
        self.eigenvalues = np.zeros(0)
        self.variance_fraction = np.zeros(0)
        self.keys = []                          # atom key of every vector component
        self.vectors = np.zeros((0, 0, 3))      # (n_modes, n_atoms, 3)
        self.origin = None                      # (n_atoms, 3) for the global set

    @property
    def n_modes(self):
        return len(self.eigenvalues)

    def summary(self):
        var = "/".join(f"{100 * v:.0f}" for v in self.variance_fraction)
        ev = "/".join(f"{v:.3g}" for v in self.eigenvalues)
        return (f"{self.scope} set {self.id}: {self.n_modes} {self.method} modes on "
                f"{len(self.keys)} atoms ({self.selection}), {self.n_frames} frames, "
                f"eigenvalues {ev} A^2, variance {var}%")


class MdcifFile:
    """Models, mode sets and frame projections of one mdcif file."""

    def __init__(self, path):
        self.path = str(path)
        t = read_tables(path)
        site = t.get("_atom_site")
        if not site:
            raise ValueError(f"{path}: no _atom_site")
        model_col = site.get("pdbx_PDB_model_num") or ["1"] * len(site["id"])
        chain = site.get("auth_asym_id") or site.get("label_asym_id")
        resnum = site.get("auth_seq_id") or site.get("label_seq_id")
        name = site.get("label_atom_id") or site.get("auth_atom_id")
        self.atom_key = {}                      # _atom_site.id -> (chain, resnum, name)
        self.atom_model = {}                    # _atom_site.id -> model number
        self.model_keys = {}                    # model number -> keys in file order
        for i, aid in enumerate(site["id"]):
            m = int(model_col[i])
            key = (str(chain[i]), str(resnum[i]), str(name[i]))
            self.atom_key[aid] = key
            self.atom_model[aid] = m
            self.model_keys.setdefault(m, []).append(key)
        self.models = sorted(self.model_keys)
        ml = t.get("_ihm_model_list", {})
        self.model_names = {int(i): n for i, n in zip(ml.get("model_id", []),
                                                       ml.get("model_name", []), strict=True)}
        self.mode_sets = self._mode_sets(t)
        self.projection = self._projection(t)

    def _mode_sets(self, t):
        info = t.get("_mdcif_mode_set")
        if not info:
            if t.get("_mdcif_mode"):
                raise ValueError(f"{self.path}: mdcif v0.1 file (no _mdcif_mode_set); rebuild it")
            return {}
        rows = [dict(zip(info, vals, strict=True)) for vals in zip(*info.values(), strict=True)]
        sets = {int(r["id"]): ModeSet(r) for r in rows}
        modes = t.get("_mdcif_mode", {})
        mode_of = {}                            # mode id -> (set id, mode index)
        per_set = {}
        for mid, sid, idx, ev, vf in zip(modes.get("id", []), modes.get("mode_set_id", []),
                                         modes.get("mode_index", []),
                                         modes.get("eigenvalue", []),
                                         modes.get("variance_fraction", []), strict=True):
            per_set.setdefault(int(sid), []).append((int(idx), mid, float(ev),
                                                     np.nan if vf is None else float(vf)))
        for sid, rows in per_set.items():
            rows.sort()
            ms = sets[sid]
            ms.eigenvalues = np.array([r[2] for r in rows])
            ms.variance_fraction = np.array([r[3] for r in rows])
            for k, r in enumerate(rows):
                mode_of[r[1]] = (sid, k)
        vec = t.get("_mdcif_mode_vector", {})
        comps = {}
        for mid, aid, dx, dy, dz in zip(vec.get("mode_id", []), vec.get("atom_id", []),
                                        _floats(vec.get("dx", [])), _floats(vec.get("dy", [])),
                                        _floats(vec.get("dz", [])), strict=True):
            sid, k = mode_of[mid]
            comps.setdefault(sid, {}).setdefault(aid, {})[k] = (dx, dy, dz)
        for sid, by_atom in comps.items():
            ms = sets[sid]
            aids = sorted(by_atom, key=int)
            ms.keys = [self.atom_key[a] for a in aids]
            ms.vectors = np.array([[by_atom[a][k] for a in aids] for k in range(ms.n_modes)])
        org = t.get("_mdcif_mode_origin", {})
        origins = {}
        for sid, aid, x, y, z in zip(org.get("mode_set_id", []), org.get("atom_id", []),
                                     _floats(org.get("x", [])), _floats(org.get("y", [])),
                                     _floats(org.get("z", [])), strict=True):
            origins.setdefault(int(sid), {})[self.atom_key[aid]] = (x, y, z)
        for sid, xyz in origins.items():
            ms = sets[sid]
            ms.origin = np.array([xyz[k] for k in ms.keys])
        return sets

    def _projection(self, t):
        p = t.get("_mdcif_projection")
        if not p:
            return None
        return {"frame": np.array([int(v) for v in p["frame"]]),
                "time": _floats(p.get("time", [None] * len(p["frame"]))),
                "model": np.array([int(v) if v is not None else 0 for v in p["model_id"]]),
                "coords": np.column_stack([_floats(p.get(f"proj_{i}", [None] * len(p["frame"])))
                                           for i in (1, 2, 3)])}

    # ------------------------------------------------------------------ lookups
    @property
    def global_set(self):
        return next((ms for ms in self.mode_sets.values() if ms.scope == "global"), None)

    def cluster_set(self, model):
        return next((ms for ms in self.mode_sets.values()
                     if ms.scope == "cluster" and ms.model == model), None)

    def mode_set(self, model, use_global=False):
        """The cluster set of ``model`` or the global set; raises if missing."""
        ms = self.global_set if use_global else self.cluster_set(model)
        if ms is None:
            what = "global modes" if use_global else f"modes for model {model}"
            raise ValueError(f"{self.path}: no {what}")
        return ms

    def model_label(self, model):
        return self.model_names.get(model) or f"model {model}"

    def summary(self):
        lines = [f"{self.path}: {len(self.models)} models"]
        g = self.global_set
        if g is not None:
            lines.append("  " + g.summary())
        if self.projection is not None:
            p = self.projection
            lines.append(f"  {len(p['frame'])} frames projected on global modes 1-3, "
                         f"{100 * np.mean(p['model'] == 0):.0f}% in no deposited cluster")
        for m in self.models:
            ms = self.cluster_set(m)
            share = ""
            if self.projection is not None:
                share = f" ({100 * np.mean(self.projection['model'] == m):.0f}% of frames)"
            lines.append(f"  model {m}: {self.model_label(m)}{share}")
            if ms is not None:
                lines.append("    " + ms.summary())
        return "\n".join(lines)


# --------------------------------------------------------------------------- geometry

def displacements(keys, ms, mode_index):
    """Unit-vector components for viewer atoms ``keys`` (list of (chain, resnum, name)).

    Returns ``(d, own)``: ``(n, 3)`` components (atoms without their own component
    take their residue's) and a boolean mask of atoms with their own component.
    """
    if not 1 <= mode_index <= ms.n_modes:
        raise ValueError(f"mode {mode_index} not in 1..{ms.n_modes}")
    v = ms.vectors[mode_index - 1]
    exact = {k: v[i] for i, k in enumerate(ms.keys)}
    by_res = {}
    for i, k in enumerate(ms.keys):
        by_res.setdefault(k[:2], v[i])
    d = np.zeros((len(keys), 3))
    own = np.zeros(len(keys), dtype=bool)
    for i, k in enumerate(keys):
        if k in exact:
            d[i], own[i] = exact[k], True
        elif k[:2] in by_res:
            d[i] = by_res[k[:2]]
    if not own.any():
        raise ValueError("no atom of the structure matches the mode atoms")
    return d, own


def mode_frames(coords, keys, ms, mode_index, amplitude=3.0, n_frames=30):
    """``(n_frames, n, 3)`` coordinates oscillating along a mode: x0 + k sqrt(l) sin(phase) d."""
    d, _ = displacements(keys, ms, mode_index)
    sigma = np.sqrt(max(ms.eigenvalues[mode_index - 1], 0.0))
    phase = np.sin(2 * np.pi * np.arange(n_frames) / n_frames)
    return np.asarray(coords, float)[None] + amplitude * sigma * phase[:, None, None] * d[None]


def morph_frames(model_coords, steps=10, wrap=True):
    """Linear interpolation through a list of ``(n, 3)`` model coordinates (same atoms)."""
    seq = list(model_coords) + ([model_coords[0]] if wrap and len(model_coords) > 1 else [])
    out = []
    for a, b in zip(seq[:-1], seq[1:], strict=True):
        for s in range(steps):
            out.append(a + (b - a) * s / steps)
    if not wrap or len(model_coords) == 1:
        out.append(seq[-1])
    return np.array(out)


def arrows(coords, keys, ms, mode_index, scale=3.0, min_length=0.5):
    """Porcupine arrows ``(starts, ends)`` for atoms with their own vector component.

    Arrow = scale * sqrt(lambda) * v_i (``scale`` in standard deviations); arrows
    shorter than ``min_length`` angstroms are dropped.
    """
    d, own = displacements(keys, ms, mode_index)
    vec = scale * np.sqrt(max(ms.eigenvalues[mode_index - 1], 0.0)) * d[own]
    start = np.asarray(coords, float)[own]
    keep = np.linalg.norm(vec, axis=1) >= min_length
    return start[keep], start[keep] + vec[keep]


def arrow_mesh(starts, ends, radius=0.15, head_radius=0.4, head_length=1.0, sides=12):
    """Triangle mesh ``(vertices, normals, triangles)`` of arrows (open cylinder + cone)."""
    verts, norms, tris = [], [], []
    ang = 2 * np.pi * np.arange(sides) / sides
    circle = np.column_stack([np.cos(ang), np.sin(ang), np.zeros(sides)])
    for s, e in zip(np.asarray(starts, float), np.asarray(ends, float), strict=True):
        axis = e - s
        length = np.linalg.norm(axis)
        if length < 1e-6:
            continue
        z = axis / length
        x = np.cross(z, [1.0, 0, 0] if abs(z[0]) < 0.9 else [0, 1.0, 0])
        x /= np.linalg.norm(x)
        y = np.cross(z, x)
        rot = np.array([x, y, z])               # rows: local axes in world coordinates
        head = min(head_length, 0.5 * length)
        neck = s + axis * (1 - head / length)
        ring = circle @ rot                     # unit circle in the plane normal to z
        base = len(verts)
        # shaft: two rings
        for p in (s, neck):
            verts.extend(p + radius * ring)
            norms.extend(ring)
        for k in range(sides):
            a, b = base + k, base + (k + 1) % sides
            tris.extend([(a, b, a + sides), (b, b + sides, a + sides)])
        # head: base ring + tip, slanted normals; plus a cap disc under the cone
        base = len(verts)
        slant = ring * (head / np.hypot(head, head_radius)) + z * (
            head_radius / np.hypot(head, head_radius))
        verts.extend(neck + head_radius * ring)
        norms.extend(slant)
        verts.extend([e] * sides)
        norms.extend(slant)
        for k in range(sides):
            a, b = base + k, base + (k + 1) % sides
            tris.append((a, b, base + sides + k))
        base = len(verts)
        verts.extend(neck + head_radius * ring)
        norms.extend([-z] * sides)
        verts.append(neck)
        norms.append(-z)
        for k in range(sides):
            tris.append((base + (k + 1) % sides, base + k, base + sides))
    return (np.array(verts, dtype=np.float32).reshape(-1, 3),
            np.array(norms, dtype=np.float32).reshape(-1, 3),
            np.array(tris, dtype=np.int32).reshape(-1, 3))
