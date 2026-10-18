"""
Invariants of pyEQL.utils.standardize_formula, checked against every aqueous species defined in the PHREEQC
databases supported by pyEQL and against every entry in the pyEQL property database.
"""

import json
import re
from collections import defaultdict
from importlib import resources
from pathlib import Path

import pytest
from pymatgen.core.ion import Ion

import pyEQL
from pyEQL.utils import _VALENCE_LABEL_RE, standardize_formula

try:
    import phreeqpython
except ImportError:  # phreeqpython is not available on all platforms
    phreeqpython = None

_PHREEQC_DIR = resources.files("pyEQL.phreeqc.database")
_PYEQL_DIR = Path(pyEQL.__file__).parent / "database"
_PHREEQPYTHON_DIR = Path(phreeqpython.__file__).parent / "database" if phreeqpython else None
PHREEQC_DATABASES = {
    "phreeqc.dat": _PHREEQC_DIR / "phreeqc.dat",
    "pitzer.dat": _PHREEQC_DIR / "pitzer.dat",
    "llnl.dat": _PYEQL_DIR / "llnl.dat",
    "geothermal.dat": _PYEQL_DIR / "geothermal.dat",
    # these databases are distributed only with phreeqpython
    "vitens.dat": _PHREEQPYTHON_DIR / "vitens.dat" if _PHREEQPYTHON_DIR else None,
    "wateq4f_PWN.dat": _PHREEQPYTHON_DIR / "wateq4f_PWN.dat" if _PHREEQPYTHON_DIR else None,
}


def _solution_species(path) -> set[str]:
    """Return the names of all species that appear in the SOLUTION_SPECIES reactions of a PHREEQC database."""
    species = set()
    in_block = False
    for raw in path.read_bytes().decode(errors="ignore").splitlines():
        line = raw.split("#")[0].rstrip()
        if not line.strip():
            continue
        if re.match(r"^[A-Z_]+(\s|$)", line):
            in_block = line.split()[0] == "SOLUTION_SPECIES"
            continue
        if not in_block or "=" not in line or line.strip().startswith("-"):
            continue
        for side in line.split("="):
            # terms are separated by ' + ' (charges such as Ca+2 contain '+' without spaces)
            for term in re.split(r"\s\+\s", side):
                if m := re.match(r"^\s*[\d.]*\s*(\S+)\s*$", term):
                    species.add(m.group(1))
    return species - {"e-", "H2O"}


def _identity(ion: Ion) -> tuple:
    return frozenset((str(el), round(amt, 6)) for el, amt in ion.composition.items()), round(ion.charge, 6)


def _parse(name: str):
    """Composition and charge of a PHREEQC species name, or None for names that are not formulas (e.g. 'Humate-2')."""
    try:
        return _identity(Ion.from_formula(_VALENCE_LABEL_RE.sub(r"\1", name)))
    except ValueError:
        return None


@pytest.mark.parametrize("db", PHREEQC_DATABASES)
def test_phreeqc_species(db):
    """Standardized formulas preserve composition and charge, are idempotent, and never collide."""
    if PHREEQC_DATABASES[db] is None:
        pytest.skip(f"{db} is distributed with phreeqpython, which is not installed")
    by_key = defaultdict(set)
    errors = []
    for name in sorted(_solution_species(PHREEQC_DATABASES[db])):
        identity = _parse(name)
        if identity is None:
            continue
        key = standardize_formula(name)
        if _identity(Ion.from_formula(key)) != identity:
            errors.append(f"{name} -> {key} changes the composition or charge")
        if standardize_formula(key) != key:
            errors.append(f"{name} -> {key} is not idempotent")
        by_key[key].add(identity)
    errors.extend(f"{key} is shared by distinct species" for key, ids in by_key.items() if len(ids) > 1)
    assert not errors, "\n".join(errors)


def test_database_keys():
    """Every database key is already standardized, and its informatics fields match its formula."""
    with open(Path(pyEQL.__file__).parent / "database" / "pyeql_db.json", encoding="utf-8") as f:
        docs = json.load(f)
    errors = []
    for doc in docs:
        key = doc["formula"]
        if standardize_formula(key) != key:
            errors.append(f"{key} standardizes to {standardize_formula(key)}")
            continue
        ion = Ion.from_formula(key)
        if doc["charge"] != ion.charge:
            errors.append(f"{key}: charge {doc['charge']} != {ion.charge}")
        if doc["n_atoms"] != ion.num_atoms:
            errors.append(f"{key}: n_atoms {doc['n_atoms']} != {ion.num_atoms}")
        stored = {k: v for k, v in doc["pmg_ion"].items() if not k.startswith("@")}
        if stored != ion.as_dict():
            errors.append(f"{key}: pmg_ion {stored} != {ion.as_dict()}")
        mw = float(doc["molecular_weight"].split()[0])
        if mw != pytest.approx(float(ion.weight), rel=1e-4):
            errors.append(f"{key}: molecular_weight {mw} != {float(ion.weight)}")
    assert not errors, "\n".join(errors)
