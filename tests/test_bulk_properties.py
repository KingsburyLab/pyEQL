"""
pyEQL test suite for bulk property calculations
===============================================

This file contains tests for the bulk properties calculated by
Solution class methods. Currently included methods are:

- hardness

"""


import pyEQL


# an empty solution should have zero hardness
def test_empty_solution():
    s1 = pyEQL.Solution()
    result = s1.hardness.magnitude
    expected = 0

    assert result == expected


# a solution with only monovalent ions should have zero hardness
def test_hardness_1():
    s1 = pyEQL.Solution([["Na+", "0.2 mol/L"], ["Cl-", "0.2 mol/L"]])
    result = s1.hardness.magnitude
    expected = 0

    assert result == expected


# a solution with only multivalent anions should have zero hardness
def test_hardness_2():
    s1 = pyEQL.Solution([["Na+", "0.4 mol/L"], ["SO4-2", "0.2 mol/L"]])
    result = s1.hardness.magnitude
    expected = 0

    assert result == expected


# the hardness should return the equivalent concentration, not just the
# molar concentration (e.g. multiply by the charge)
def test_hardness_3():
    s1 = pyEQL.Solution([["Fe+3", "0.1 mol/L"], ["Cl-", "0.3 mol/L"]])
    result = s1.hardness.magnitude
    expected = 15013.5

    assert round(abs(result - expected), 7) == 0


# the hardness should account for multiple cations but count only those
# that are multivalent
def test_hardness_4():
    s1 = pyEQL.Solution(
        [
            ["Na+", "0.1 mol/L"],
            ["K+", "0.1 mol/L"],
            ["Mg+2", "0.1 mol/L"],
            ["Ca+2", "0.1 mol/L"],
            ["Fe+3", "0.1 mol/L"],
            ["Cl-", "0.1 mol/L"],
            ["F-", "0.1 mol/L"],
            ["SO4-2", "0.2 mol/L"],
            ["PO4-3", "0.1 mol/L"],
        ]
    )
    result = s1.hardness.magnitude
    expected = 35031.5

    assert round(abs(result - expected), 7) == 0


# the hardness should return g/L units
def test_hardness_5():
    s1 = pyEQL.Solution([["Fe+3", "0.1 mol/L"], ["Cl-", "0.3 mol/L"]])
    result = str(s1.hardness.dimensionality)
    expected = "[mass] / [length] ** 3"

    assert result == expected
