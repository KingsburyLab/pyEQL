"""
pyEQL utilities

:copyright: 2013-2024 by Ryan S. Kingsbury
:license: LGPL, see LICENSE for more details.

"""

import logging
import math
import re
from collections import UserDict
from fractions import Fraction
from functools import lru_cache
from typing import Any

from iapws import IAPWS95, IAPWS97
from pymatgen.core import Composition, Element
from pymatgen.core.ion import Ion
from pymatgen.util.string import charge_string, formula_double_format

from pyEQL import ureg

logger = logging.getLogger(__name__)


def translate_units(unit: str) -> str:
    """
    Translate commonly used environmental units such as 'ppm' into strings that `pint` can understand.

    Args:
        unit: string representing the unit to translate

    Returns: a unit that pint can understand
    """
    if unit == "m":  # molal
        return "mol/kg"
    if unit == "ppm":
        return "mg/L"
    if unit == "ppb":
        return "ug/L"
    if unit == "ppt":
        return "ng/L"
    # if all else fails, return the unit we were provided
    return unit


def _translate_pint_quantity(amount: str):
    """
    Helper method to split a pint quantity string into magnitude and units.
    """
    import re  # noqa: PLC0415

    from pint import Quantity  # noqa: PLC0415

    # skip if already a pint Quantity
    if isinstance(amount, Quantity):
        return amount.magnitude, str(amount.units)

    match = re.match(r"^\s*([0-9eE+\-*/().]+)\s*(.*)$", amount)

    if match is None:
        return amount

    _value, _unit = match.groups()
    # handle python ** expression in Pint quantity
    _value = eval(_value) if "**" in _value else float(_value)
    unit = translate_units(_unit)
    return (float(_value), unit)


# Formulas of polymers such as "(CO2)2", which PHREEQC reports for some solutions, are left intact.
_POLYMER_RE = re.compile(r"^\([A-Za-z0-9+-]+\)\d+$")

# Element labels that PHREEQC databases (e.g. vitens.dat) use to decouple redox states, e.g. "[N-3]H4+" or "H3[As+3]O3".
# Charges in brackets (e.g. "Na[+1]") never start with a capital letter, so they are not matched.
_VALENCE_LABEL_RE = re.compile(r"\[([A-Z][a-z]?)[+-]\d*\]")

# A formula ending in a number followed by a sign, e.g. "NH4+", "Ca2+", or "SO42-". The final digit can be either a
# subscript (NH4+) or a charge number (Ca2+).
_AMBIGUOUS_CHARGE_RE = re.compile(r"^(.*[A-Za-z)])(\d+)([+-])$")
_ELEMENT_RE = re.compile(r"[A-Z][a-z]?")

# Keys produced by older versions of pyEQL whose integer multiple cannot be inferred from the charge alone.
_LEGACY_FORMULAS = {
    "S[-0.4]": "S5[-2]",
    "S[-0.5]": "S4[-2]",
    "S[-0.66666667]": "S3[-2]",
}

# Preferred (condensed or conventional) names for specific species. Entries are matched by composition and charge,
# so they apply regardless of how the formula was written (e.g. "C2H3O2-" -> "CH3COO[-1]").
_PREFERRED_FORMULAS = [
    "NH3(aq)",
    "NH4NO3(aq)",
    "(NH4)2SO4(aq)",
    "NH4SO4[-1]",
    "H3PO4(aq)",
    "CH3COOH(aq)",
    "CH3COO[-1]",
    "CH3OH(aq)",
    "C2H5OH(aq)",
    "C3H7OH(aq)",
    "C4H9OH(aq)",
    "HCOOH(aq)",
    "HCOO[-1]",
    "CF3SO3[-1]",
    "CF3COO[-1]",
    "CCl3COO[-1]",
    "CBr3COO[-1]",
    "CI3COO[-1]",
    "CFCl2COO[-1]",
    "CF2ClCOO[-1]",
    "CBrCl2COO[-1]",
    "CBr2ClCOO[-1]",
    "CICl2COO[-1]",
    "CI2ClCOO[-1]",
]

# Rewrites applied to electronegativity-ordered (inorganic) formulas so that common groups appear in their
# conventional order, including inside complexes (e.g. "CaPHO4" -> "CaHPO4"). The lookahead prevents a match from
# running into a following subscript or lowercase letter (e.g. "H4N2" is not rewritten).
_FRAGMENT_REWRITES = [
    # hydroxide groups, e.g. "Al(HO)4" -> "Al(OH)4" (but not "U2(HO3)2")
    (re.compile(r"HO(?![0-9a-z])"), "OH"),
    (re.compile(r"P\(HO2\)2(?![0-9a-z])"), "H2PO4"),
    (re.compile(r"PHO4(?![0-9a-z])"), "HPO4"),
    (re.compile(r"CSN(?![0-9a-z])"), "SCN"),
    (re.compile(r"H4N(?![0-9a-z])"), "NH4"),
]


def _normalize_formula_text(formula: str) -> str:
    """Replace unicode superscripts, subscripts, and dashes with their ASCII equivalents."""
    # fix permuted sign and charge number (e.g. Co²⁺)
    for s, rep in zip(["²⁺", "³⁺", "⁴⁺", "²⁻", "³⁻", "⁴⁻"], ["+2", "+3", "+4", "-2", "-3", "-4"], strict=False):
        formula = formula.replace(s, rep)

    # replace superscripts with non superscripts
    for char, rep in zip("⁻⁺⁰¹²³⁴⁵⁶⁷⁸⁹", "-+0123456789", strict=False):
        formula = formula.replace(char, rep)

    # replace subscripts with non subscripts
    for char, rep in zip("₀₁₂₃₄₅₆₇₈₉", "0123456789", strict=False):
        formula = formula.replace(char, rep)

    # replace different types of dashes with a minus sign
    for char in [r"‑", r"‐", r"‒", r"–", r"—", r"−"]:  # noqa: RUF001
        formula = formula.replace(char, "-")

    return formula


def _is_charge_balanced(formula: str, charge: int) -> bool:
    """Whether pymatgen can find a set of oxidation states that sums to `charge` for `formula`."""
    try:
        return bool(Composition(formula).oxi_state_guesses(target_charge=charge))
    except Exception:
        return False


def _resolve_ambiguous_charge(formula: str) -> str:
    """
    Disambiguate a formula that ends in a number followed by a sign, e.g. "Ca2+" or "SO42-".

    `pymatgen` reads the number as a subscript, so "Ca2+" would mean Ca2 with a charge of +1. That reading is kept
    in most cases (e.g. "NH4+", "NO3-", "ClO3-", "I3-"). The final digit is instead read as the charge number for
    1) monatomic metal ions for which that charge is a known oxidation state (e.g. "Ca2+", "Fe3+"), 2) monatomic
    chalcogenide ions (e.g. "S2-"), and 3) polyatomic ions in which the subscript reading would give the last element
    a subscript of 10 or more and the charge reading can be charge balanced (e.g. "SO42-", "Cr2O72-").

    A warning is logged whenever the charge reading is used, and also when the subscript reading is kept for a metal
    oxo-cation whose charge reading is also plausible, e.g. "VO2+" could mean VO2[+1] or vanadyl, VO[+2].
    """
    match = _AMBIGUOUS_CHARGE_RE.match(formula)
    if match is None:
        return formula

    head, digits, sign = match.groups()
    charge = int(f"{sign}{digits[-1]}")
    if digits[:-1] == "" and _ELEMENT_RE.fullmatch(head):
        try:
            el = Element(head)
        except ValueError:
            return formula
        reinterpret = (el.is_metal and charge in el.oxidation_states) or (head in {"S", "Se", "Te"} and charge == -2)
    else:
        reinterpret = int(digits) >= 10 and _is_charge_balanced(f"{head}{digits[:-1]}", charge)

    as_subscript = f"{head}{digits}[{sign}1]"
    as_charge = f"{head}{digits[:-1]}[{sign}{digits[-1]}]"
    if not reinterpret:
        if len(digits) == 1 and _is_metal_oxo_cation(head, charge):
            logger.warning(
                f"Interpreting ambiguous formula {formula} as {as_subscript}, not {as_charge}. Write the charge "
                f"after the sign (e.g. 'VO+2') or in brackets (e.g. 'VO[+2]') to avoid ambiguity."
            )
        return formula

    logger.warning(
        f"Interpreting ambiguous formula {formula} as {as_charge}, not {as_subscript}. Write the charge after the "
        f"sign (e.g. 'Ca+2') or in brackets (e.g. 'Ca[+2]') to avoid ambiguity."
    )
    return as_charge


def _is_metal_oxo_cation(formula: str, charge: int) -> bool:
    """
    Whether `formula` with `charge` is a plausible oxo-cation of a metal, e.g. VO[+2] (vanadyl) or TiO[+2].

    Actinides are excluded because their mono-oxo cations (e.g. UO[+2]) are not aqueous species, whereas their
    dioxo cations (e.g. UO2[+1], PuO2[+1]) are.
    """
    if charge <= 0:
        return False
    try:
        symbols = {el.symbol for el in Composition(formula).elements}
    except Exception:
        return False
    if len(symbols) != 2 or "O" not in symbols:
        return False
    metal = Element((symbols - {"O"}).pop())
    return metal.is_metal and not metal.is_actinoid and _is_charge_balanced(formula, charge)


def _is_organic(el_amt: dict[str, float]) -> bool:
    """
    Whether a species should be written in Hill order, based only on its composition.

    Carbon-containing species are treated as organic, except for inorganic carbon: carbonates, bicarbonates, and
    carbonic acid (at least 3 O per C, not counting 3 O per S), cyanides, cyanates, and thiocyanates (at least 1 N
    and at most 1 H per C), and CO, CO2, CS2, and COS. Complexes of these species with metals are also inorganic.
    """
    n_c = el_amt.get("C", 0)
    if n_c == 0:
        return False
    n_h = el_amt.get("H", 0)
    if el_amt.get("O", 0) - 3 * el_amt.get("S", 0) >= 3 * n_c:
        return False
    if el_amt.get("N", 0) >= n_c and n_h <= n_c:
        return False
    return not (n_c == 1 and n_h == 0 and set(el_amt) <= {"C", "O", "S"})


def _electronegativity(symbol: str) -> float:
    try:
        x = Element(symbol).X
    except ValueError:  # pymatgen DummySpecies, e.g. "Ethylene" in llnl.dat
        return float("inf")
    return float("inf") if math.isnan(x) else x


def _format_formate_salt(el_amt: dict[str, float]) -> str | None:
    """
    Write formate salts and complexes with the metal first, e.g. "NaHCOO" or "Ca(HCOO)2".

    Returns None unless the species consists of one metal and n formate (HCOO) groups.
    """
    n = el_amt.get("C", 0)
    metals = set(el_amt) - {"C", "H", "O"}
    if len(metals) != 1 or el_amt.get("H") != n or el_amt.get("O") != 2 * n:
        return None
    (metal,) = metals
    try:
        if not Element(metal).is_metal:
            return None
    except ValueError:  # pymatgen DummySpecies
        return None
    groups = "HCOO" if n == 1 else f"(HCOO){formula_double_format(n)}"
    return f"{metal}{formula_double_format(el_amt[metal])}{groups}"


def _format_formula(el_amt: dict[str, float]) -> str:
    """
    Write a formula (without charge) from an element: amount dict, without reducing it.

    Organic species are written in Hill order (C, then H, then other elements alphabetically), except for formate
    salts, which are written with the metal first (e.g. "NaHCOO"). Other species are written in order of increasing
    electronegativity, with the two most electronegative elements grouped in parentheses when they share a common
    factor (e.g. "Al(OH)4"), consistent with `pymatgen`.
    """
    if _is_organic(el_amt):
        if (formate := _format_formate_salt(el_amt)) is not None:
            return formate
        syms = ["C"] + (["H"] if "H" in el_amt else []) + sorted(set(el_amt) - {"C", "H"})
        return "".join(f"{s}{formula_double_format(el_amt[s])}" for s in syms)

    syms = sorted(el_amt, key=lambda s: (_electronegativity(s), s))
    group = ""
    if len(syms) >= 3 and _electronegativity(syms[-1]) - _electronegativity(syms[-2]) < 1.65:
        a, b = syms[-2:]
        if el_amt[a] == int(el_amt[a]) and el_amt[b] == int(el_amt[b]):
            factor = math.gcd(int(el_amt[a]), int(el_amt[b]))
            if factor > 1:
                n_a, n_b = formula_double_format(el_amt[a] / factor), formula_double_format(el_amt[b] / factor)
                group = f"({a}{n_a}{b}{n_b}){factor}"
                syms = syms[:-2]
    formula = "".join(f"{s}{formula_double_format(el_amt[s])}" for s in syms) + group
    for pattern, rep in _FRAGMENT_REWRITES:
        formula = pattern.sub(rep, formula)
    return formula


def _composition_key(el_amt: dict[str, float], charge: float) -> tuple[frozenset, float]:
    return frozenset(el_amt.items()), charge


@lru_cache(maxsize=1)
def _preferred_formulas() -> dict[tuple[frozenset, float], str]:
    """Map (composition, charge) to the preferred formula (without charge) for each entry in _PREFERRED_FORMULAS."""
    preferred = {}
    for f in _PREFERRED_FORMULAS:
        ion = Ion.from_formula(f)
        base = f[: f.index("[")] if "[" in f else f.removesuffix("(aq)")
        preferred[_composition_key(ion.get_el_amt_dict(), ion.charge)] = base
    return preferred


@lru_cache
def standardize_formula(formula: str):
    """
    Convert a chemical formula into standard form.

    Args:
        formula: the chemical formula to standardize.

    Returns:
        A standardized chemical formula

    Raises:
        ValueError if `formula` cannot be processed or is invalid.

    Notes:
        The formula is parsed by `pymatgen.core.ion.Ion` into element amounts and a charge, and the standardized
        formula is written from those amounts *without reducing them*, so that multimers such as `Hg2[+2]`,
        `S2O8[-2]`, or `Fe2(OH)2[+4]` keep their identity. In the standardized formula:

        1. Charged species always list the charge number explicitly, enclosed in square brackets, e.g. 'Na+',
           'Na+1', and 'Na[+]' all standardize to "Na[+1]". Uncharged species end in "(aq)".
        2. Organic species are written in Hill order (e.g. "C2H6(aq)"); other species are written in order of
           increasing electronegativity (e.g. "NaCl(aq)", "Al(OH)4[-1]").
        3. A small number of species use a preferred conventional name (e.g. "NH4[+1]", "CH3COO[-1]").

        Formulas that end in a number followed by a sign without brackets (e.g. "Ca2+" or "SO42-") are ambiguous.
        They are interpreted as described in `_resolve_ambiguous_charge` and a warning is logged.

    Examples:
        >>> from pyEQL.utils import standardize_formula
        >>> standardize_formula("Na+")
        'Na[+1]'
        >>> standardize_formula("Hg2+2")
        'Hg2[+2]'
        >>> standardize_formula("Fe(OH)2+")
        'Fe(OH)2[+1]'
        >>> standardize_formula("SO42-")
        'SO4[-2]'
    """
    formula = _normalize_formula_text(formula)

    # Do not modify any dimers etc (Phreeqc reports a small amount of
    # "(CO2)2" in a water solution with C(4), for example.
    if _POLYMER_RE.match(formula):
        return formula

    if formula in _LEGACY_FORMULAS:
        return _LEGACY_FORMULAS[formula]

    formula = _VALENCE_LABEL_RE.sub(r"\1", formula)
    formula = _resolve_ambiguous_charge(formula)

    ion = Ion.from_formula(formula)
    el_amt = ion.get_el_amt_dict()
    charge = ion.charge

    # Fractional charges (e.g. "I[-0.33333333]") are produced only by reducing the formula of a multimer, as was done
    # by older versions of pyEQL. Restore the smallest multiple that has an integer charge (e.g. "I3[-1]").
    if not math.isclose(charge, round(charge), abs_tol=1e-6):
        factor = Fraction(charge).limit_denominator(12).denominator
        el_amt = {el: amt * factor for el, amt in el_amt.items()}
        charge = round(charge * factor)

    base = _preferred_formulas().get(_composition_key(el_amt, charge))
    if base is None:
        base = _format_formula(el_amt)
    return base + charge_string(charge)


def format_solutes_dict(solute_dict: dict, units: str):
    """
    Formats a dictionary of solutes by converting the amount to a string with the provided units suitable for passing to
    use with the Solution class. Note that all solutes must be given in the same units.

    Args:
        solute_dict: The dictionary to format. This must be of the form dict{str: Number}
            e.g. {"Na+": 0.5, "Cl-": 0.9}
        units: The units to use for the solute. e.g. "mol/kg"

    Returns:
        A formatted solute dictionary.

    Raises:
        TypeError if `solute_dict` is not a dictionary.
    """
    if not isinstance(solute_dict, dict):
        raise TypeError("solute_dict must be a dictionary. Refer to the doc for proper formatting.")

    return {key: f"{value!s} {units}" for key, value in solute_dict.items()}


@lru_cache
@ureg.wraps(ret=None, args=["K", "MPa"], strict=False)
def create_water_substance(temperature: float, pressure: float):
    """
    Instantiate a water substance model from IAPWS.

    Args:
        temperature: the desired temperature in K
        pressure: the desired pressure in MPa

    Notes:
        The IAPWS97 model is much faster than IAPWS95, but the latter can do temp
        below zero. See https://github.com/jjgomera/iapws/issues/14. Hence,
        IAPWS97 will be used except when `temperature` is less than 0 degC.

    Returns:
        A IAPWS97 or IAPWS95 instance
    """
    if temperature >= 273.15:
        return IAPWS97(T=temperature, P=pressure)
    return IAPWS95(T=temperature, P=pressure)


class FormulaDict(UserDict):
    """
    Automatically converts keys on get/set using pyEQL.utils.standardize_formula(key).

    This allows getting/setting/updating of Solution.components using flexible
    formula notation (e.g., "Na+", "Na+1", "Na[+]" all have the same effect)
    """

    def __init__(self, *args, **kwargs) -> None:
        # Whether self.data is currently in sorted (descending-by-value) order. Set before
        # calling super().__init__ because that may trigger __setitem__ during construction.
        self._sorted = True
        super().__init__(*args, **kwargs)

    def _ensure_sorted(self) -> None:
        """Sort the underlying data by value (descending) if a mutation invalidated the order."""
        if not self._sorted:
            self.data = dict(sorted(self.data.items(), key=lambda x: x[1], reverse=True))
            self._sorted = True

    def __getitem__(self, key) -> Any:
        return super().__getitem__(standardize_formula(key))

    def __setitem__(self, key, value) -> None:
        # ensure that all values are stored as python floats, not numpy types
        # see https://numpy.org/doc/stable/release/2.0.0-notes.html#representation-of-numpy-scalars-changed
        super().__setitem__(standardize_formula(key), float(value))
        # defer sorting until the contents are next iterated (lazy sort via _ensure_sorted),
        # which avoids an O(n log n) re-sort on every assignment
        self._sorted = False

    def __iter__(self):
        # Iteration is the read barrier: keys(), values(), items(), list(), dict(), and
        # `for` loops all funnel through here, so sorting now guarantees callers observe
        # the same descending-by-amount order they did under eager sorting.
        self._ensure_sorted()
        return iter(self.data)

    def __repr__(self) -> str:
        # UserDict.__repr__ reads self.data directly, bypassing __iter__; sort first so the
        # printed representation matches iteration order.
        self._ensure_sorted()
        return super().__repr__()

    # Necessary to define this so that .get() works properly in python 3.12+
    # see https://github.com/python/cpython/issues/105524
    def __contains__(self, key) -> bool:
        return standardize_formula(key) in self.data

    def __delitem__(self, key) -> None:
        super().__delitem__(standardize_formula(key))
        # deleting a key preserves the relative order of the remaining items, so the
        # sorted/dirty state is unchanged.
