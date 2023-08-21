"""
pyEQL utilities

:copyright: 2023 by Ryan S. Kingsbury
:license: LGPL, see LICENSE for more details.

"""

from collections import UserDict

from pymatgen.core.ion import Ion


class FormulaDict(UserDict):
    """
    Automatically converts keys on get/set using pymatgen.core.Ion.from_formula(key).reduced_formula.

    This allows getting/setting/updating of Solution.components using flexible
    formula notation (e.g., "Na+", "Na+1", "Na[+]" all have the same effect)
    """

    def __getitem__(self, key):
        return super().__getitem__(Ion.from_formula(key).reduced_formula)

    def __setitem__(self, key, value):
        super().__setitem__(Ion.from_formula(key).reduced_formula, value)

    def __delitem__(self, key):
        super().__delitem__(Ion.from_formula(key).reduced_formula)
