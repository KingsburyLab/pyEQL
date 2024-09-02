"""
pyEQL utilities

:copyright: 2013-2024 by Ryan S. Kingsbury
:license: LGPL, see LICENSE for more details.

"""

import logging
from collections import UserDict
from functools import lru_cache
from typing import Any

from iapws import IAPWS95, IAPWS97
from pymatgen.core.ion import Ion

from pyEQL import ureg

logger = logging.getLogger(__name__)


def interpret_units(unit: str) -> str:
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
        Currently this method standardizes formulae by passing them through `pymatgen.core.ion.Ion.reduced_formula()`.
        For ions, this means that 1) the charge number will always be listed explicitly and 2) the charge number will
        be enclosed in square brackets to remove any ambiguity in the meaning of the formula. For example, 'Na+',
        'Na+1', and 'Na[+]' will all standardize to "Na[+1]"
    """
    # fix permuted sign and charge number (e.g. Co2+)
    for str, rep in zip(["²⁺", "³⁺", "⁴⁺", "²⁻", "³⁻", "⁴⁻"], ["+2", "+3", "+4", "-2", "-3", "-4"]):
        formula = formula.replace(str, rep)

    # replace superscripts with non superscripts
    for char, rep in zip("⁻⁺⁰¹²³⁴⁵⁶⁷⁸⁹", "-+0123456789"):
        formula = formula.replace(char, rep)

    # replace subscripts with non subscripts
    for char, rep in zip("₀₁₂₃₄₅₆₇₈₉", "0123456789"):
        formula = formula.replace(char, rep)

    sform = Ion.from_formula(formula).reduced_formula

    # TODO - manual formula adjustments. May be implemented upstream in pymatgen in the future
    # thanks to @xiaoxiaozhu123 for pointing out these issues in
    # https://github.com/KingsburyLab/pyEQL/issues/136

    # ammonia
    if sform == "H4N[+1]":
        sform = "NH4[+1]"
    elif sform == "H3N(aq)":
        sform = "NH3(aq)"
    # phosphoric acid system
    elif sform == "PH3O4(aq)":
        sform = "H3PO4(aq)"
    elif sform == "PHO4[-2]":
        sform = "HPO4[-2]"
    elif sform == "P(HO2)2[-1]":
        sform = "H2PO4[-1]"
    # thiocyanate
    elif sform == "CSN[-1]":
        sform = "SCN[-1]"
    # triiodide, nitride, an phosphide
    elif sform == "I[-0.33333333]":
        sform = "I3[-1]"
    elif sform == "N[-0.33333333]":
        sform = "N3[-1]"
    elif sform == "P[-0.33333333]":
        sform = "P3[-1]"
    # formate
    elif sform == "HCOO[-1]":
        sform = "HCO2[-1]"
    # oxalate
    elif sform == "CO2[-1]":
        sform = "C2O4[-2]"
    # triflate
    elif sform == "CS(OF)3[-1]":
        sform = "CF3SO3[-1]"
    # haloacetic acids of F, Cl, Br, I
    elif sform == "C2Cl3O2[-1]":
        sform = "CCl3COO[-1]"
    elif sform == "C2O2F3[-1]":
        sform = "CF3COO[-1]"
    elif sform == "C2I3O2[-1]":
        sform = "CI3COO[-1]"
    elif sform == "C2Br3O2[-1]":
        sform = "CBr3COO[-1]"

    # Cl+F
    elif sform == "C2Cl2O2F[-1]":
        sform = "CFCl2COO[-1]"
    elif sform == "C2Cl(OF)2[-1]":
        sform = "CF2ClCOO[-1]"

    # Cl+Br
    elif sform == "C2Br(ClO)2[-1]":
        sform = "CBrCl2COO[-1]"
    elif sform == "C2Br2ClO2[-1]":
        sform = "CBr2ClCOO[-1]"

    # Cl+I
    elif sform == "C2I(ClO)2[-1]":
        sform = "CICl2COO[-1]"
    elif sform == "C2I2ClO2[-1]":
        sform = "CI2ClCOO[-1]"

    # ammonium sulfate salts
    elif sform == "H8S(NO2)2(aq)":
        sform = "(NH4)2SO4(aq)"
    elif sform == "H4SNO4[-1]":
        sform = "NH4SO4[-1]"

    # TODO - consider adding recognition of special formulas like MeOH for methanol or Cit for citrate
    return sform


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
    Automatically converts keys on get/set using pymatgen.core.Ion.from_formula(key).reduced_formula.

    This allows getting/setting/updating of Solution.components using flexible
    formula notation (e.g., "Na+", "Na+1", "Na[+]" all have the same effect)
    """

    def __getitem__(self, key) -> Any:
        return super().__getitem__(standardize_formula(key))

    def __setitem__(self, key, value) -> None:
        super().__setitem__(standardize_formula(key), value)
        # sort contents anytime an item is set
        self.data = dict(sorted(self.items(), key=lambda x: x[1], reverse=True))

    # Necessary to define this so that .get() works properly in python 3.12+
    # see https://github.com/python/cpython/issues/105524
    def __contains__(self, key) -> bool:
        return standardize_formula(key) in self.data

    def __delitem__(self, key) -> None:
        super().__delitem__(standardize_formula(key))
