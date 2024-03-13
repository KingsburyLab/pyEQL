"""
pyEQL utilities

:copyright: 2013-2024 by Ryan S. Kingsbury
:license: LGPL, see LICENSE for more details.

"""

import logging
from collections import UserDict
from functools import lru_cache

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
    return Ion.from_formula(formula).reduced_formula


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

    def __getitem__(self, key):
        return super().__getitem__(standardize_formula(key))

    def __setitem__(self, key, value):
        super().__setitem__(standardize_formula(key), value)
        # sort contents anytime an item is set
        self.data = dict(sorted(self.items(), key=lambda x: x[1], reverse=True))

    # Necessary to define this so that .get() works properly in python 3.12+
    # see https://github.com/python/cpython/issues/105524
    def __contains__(self, key):
        return standardize_formula(key) in self.data

    def __delitem__(self, key):
        super().__delitem__(standardize_formula(key))
