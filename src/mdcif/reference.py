"""Experimental reference structure: a PDB entry in PDBx/mmCIF.

``mdcif build --reference 3GB1`` (a PDB ID, downloaded once from RCSB and
cached, or a local ``.cif`` / ``.bcif`` file) makes the mdcif file an edited
copy of the PDB entry's description of the molecule:

- its entities, asym units (chain IDs, author numbering), source organism,
  sequence database references and citations are reused, so the MD models are
  described exactly as the PDB describes the experimental structure;
- its experimental models are not copied; the entry is cited as a dataset and
  the deposited models are superposed onto its first model;
- the CA RMSD of each deposited model to that first model is recorded in an
  IHM validation step.
"""

import os
import re
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

RCSB_URL = "https://files.rcsb.org/download/{}.cif"
PDB_ID = re.compile(r"[0-9][A-Za-z0-9]{3}")


def cache_dir():
    base = os.environ.get("MDCIF_CACHE") or Path.home() / ".cache" / "mdcif"
    path = Path(base) / "pdb"
    path.mkdir(parents=True, exist_ok=True)
    return path


def fetch(pdb_id):
    """Path of ``<pdb_id>.cif`` in the cache, downloading it from RCSB the first time."""
    path = cache_dir() / f"{pdb_id.lower()}.cif"
    if not path.exists():
        with urllib.request.urlopen(RCSB_URL.format(pdb_id.upper()), timeout=120) as r:
            path.write_bytes(r.read())
    return path


@dataclass
class Reference:
    system: object # ihm.System read from the entry
    path: Path
    pdb_id: str | None  # entry ID if it looks like a PDB ID
    method: str | None  # _exptl.method, e.g. "SOLUTION NMR"
    chain_map: dict = field(default_factory=dict)  # MD chain -> (ihm.AsymUnit, seq_id offset)

    @property
    def label(self): # print
        return f"PDB {self.pdb_id}" if self.pdb_id else self.path.name

    def first_model(self): # return first model
        for _, model in self.system._all_models():
            return model
        raise ValueError(f"{self.path}: no models in the reference entry")


def load_reference(ref) -> Reference:
    """Read a reference entry from a file path or a PDB ID."""
    import gemmi
    import ihm.reader

    path = Path(ref)
    if not path.exists():
        if not PDB_ID.fullmatch(str(ref)):
            raise FileNotFoundError(f"reference {ref!r} is neither a file nor a PDB ID")
        path = fetch(str(ref))
    binary = path.suffix.lower() == ".bcif"
    with open(path, "rb" if binary else "r", **({} if binary else {"encoding": "utf-8"})) as fh:
        system = ihm.reader.read(fh, format="BCIF" if binary else "mmCIF")[0]
    method = None
    if not binary:
        method = gemmi.cif.read(str(path)).sole_block().find_value("_exptl.method")
        method = gemmi.cif.as_string(method) if method else None
    pdb_id = system.id if system.id and PDB_ID.fullmatch(system.id) else None
    for entity in system.entities:
        _repair_references(entity)
    return Reference(system=system, path=path, pdb_id=pdb_id.upper() if pdb_id else None,
                     method=method)


def _aligned(entity, ref):
    """Whether ``ref`` (an ihm.reference.Sequence) matches the entity over its alignments."""
    seq, refseq = entity_sequence(entity), ref.sequence
    alignments = list(ref.alignments) or [None]
    for a in alignments:
        if a is not None and a.seq_dif:
            continue  # mutations: leave the check to python-ihm
        db_b, en_b = (a.db_begin, a.entity_begin) if a is not None else (1, 1)
        db_e = a.db_end if a is not None and a.db_end else len(refseq)
        en_e = a.entity_end if a is not None and a.entity_end else len(seq)
        if refseq[db_b - 1:db_e] != seq[en_b - 1:en_e]:
            return False
    return True


def _repair_references(entity):
    """Keep only sequence-database references (struct_ref) consistent with the entity.

    python-ihm pads ``pdbx_seq_one_letter_code`` with ``db_align_beg - 1`` gaps, which
    misaligns entries that store the full database sequence (e.g. 2GB1); such
    references are repaired by removing the padding, others that still do not match
    are dropped with a warning (python-ihm refuses to write them).
    """
    import warnings

    kept = []
    for ref in entity.references:
        if not isinstance(getattr(ref, "sequence", None), str) or _aligned(entity, ref):
            kept.append(ref)
            continue
        padded = ref.sequence
        ref.sequence = padded.lstrip("-")
        if _aligned(entity, ref):
            kept.append(ref)
            continue
        ref.sequence = padded
        warnings.warn(f"dropped sequence reference {ref} of entity {entity.description!r}: "
                      "its alignment does not match the entity sequence", stacklevel=3)
    entity.references[:] = kept


def entity_sequence(entity):
    """One-letter sequence of an ihm.Entity (canonical codes)."""
    return "".join(getattr(c, "code_canonical", None) or "X" for c in entity.sequence)


def match_chains(reference: Reference, chain_sequences):
    """Map each MD chain onto a reference asym unit whose entity sequence contains it.

    ``chain_sequences`` is ``{chain_id: one-letter sequence}``. A reference
    chain with the same ID is preferred. The MD residues get the entity's
    ``seq_id`` numbering: ``seq_id = offset + index`` (1-based index).
    Fills and returns ``reference.chain_map``.
    """
    asyms = [a for a in reference.system.asym_units if a.entity.is_polymeric()]
    used = set()
    for cid, seq in chain_sequences.items():
        candidates = sorted((a for a in asyms if id(a) not in used),
                            key=lambda a: (a.strand_id or a.id) != cid)
        for asym in candidates:
            offset = entity_sequence(asym.entity).find(seq)
            if offset >= 0:
                reference.chain_map[cid] = (asym, offset)
                used.add(id(asym))
                break
        else:
            raise ValueError(f"sequence of chain {cid} ({len(seq)} residues) is not found in "
                             f"any polymer entity of {reference.label}")
    return reference.chain_map


def reference_coordinates(reference: Reference, keys, atom_name="CA"):
    """Coordinates of ``atom_name`` in the reference's first model for each ``(asym, seq_id)``.

    Returns ``(found, xyz)``: a boolean mask over ``keys`` and the coordinates
    of the atoms found (first alternate location only).
    """
    xyz = {}
    for atom in reference.first_model().get_atoms():
        if atom.atom_id != atom_name:
            continue
        key = (id(atom.asym_unit), atom.seq_id)
        xyz.setdefault(key, (atom.x, atom.y, atom.z))
    found = np.array([(id(a), s) in xyz for a, s in keys])
    coords = np.array([xyz[(id(a), s)] for (a, s), f in zip(keys, found, strict=True) if f],
                      dtype=np.float64).reshape(-1, 3)
    return found, coords
