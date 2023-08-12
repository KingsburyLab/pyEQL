"""
pyEQL volume and concentration methods test suite
=================================================

This file contains tests for the volume and concentration-related methods
used by pyEQL's Solution class
"""

import numpy as np
import pytest
from pyEQL import Solution


@pytest.fixture()
def s1():
    return Solution(volume="2 L")


@pytest.fixture()
def s2():
    return Solution([["Na+", "4 mol/L"], ["Cl-", "4 mol/L"]], volume="2 L")


@pytest.fixture()
def s3():
    return Solution([["Na+", "4 mol/kg"], ["Cl-", "4 mol/kg"]], volume="2 L")


@pytest.fixture()
def s4():
    return Solution([["Na+", "8 mol"], ["Cl-", "8 mol"]], volume="2 L")


def test_empty_solution_3():
    # create an empty solution
    s1 = Solution(database=None)
    # It should return type Solution
    assert isinstance(s1, Solution)
    # It should have exactly 1L volume
    assert s1.get_volume().to("L").magnitude == 1.0
    #  the solvent should be water
    assert s1.solvent == "H2O"
    # It should have 0.997 kg water mass
    assert np.isclose(s1.get_solvent_mass().to("kg").magnitude, 0.9970415)
    # the temperature should be 25 degC
    assert s1.temperature.to("degC").magnitude == 25
    # the pressure should be 1 atm
    assert s1.pressure.to("atm").magnitude == 1
    # the pH should be 7.0
    assert np.isclose(s1.get_activity("H+"), 1e-7, atol=1e-9)
    assert np.isclose(s1.pH, 7.0, atol=0.01)
    assert np.isclose(s1.pE, 8.5)
    # it should contain H2O, H+, and OH- species
    assert set(s1.list_solutes()) == {"H2O", "OH-", "H+"}


# create an empty and test solutions with the same volume using substance / volume,
# substance/mass, and substance units
def test_solute_addition(s2, s3, s4):
    # if solutes are added at creation-time with substance / volume units,
    # then the total volume of the solution should not change (should remain at 2 L)
    assert s2.get_volume().to("L").magnitude == 2

    # if solutes are added at creation-time with substance / volume units,
    # then the resulting mol/L concentrations should be exactly what was specified
    assert s2.get_amount("Na+", "mol/L").magnitude == 4

    # if solutes are added at creation-time with substance / mass units,
    # then the resulting mol/kg concentrations should be exactly what was specified
    assert s3.get_amount("Na+", "mol/kg").magnitude == 4

    # the water mass of solution s2 should be less than that of s3, because
    # of the volume recalculation
    result_molL = s2.get_solvent_mass().to("kg").magnitude
    result_molkg = s3.get_solvent_mass().to("kg").magnitude
    assert result_molL < result_molkg

    # if solutes are added at creation-time with substance units,
    # then the resulting mol amounts should be exactly what was specified
    assert s4.get_amount("Na+", "mol").magnitude == 8

    # the water mass of solution s2 should be less than that of s4, because
    # of the volume recalculation
    result_molL = s2.get_solvent_mass().to("kg").magnitude
    result_mol = s4.get_solvent_mass().to("kg").magnitude
    assert result_molL < result_mol


# def test_serialization(s1, s2, s3):
#     assert isinstance(s1.as_dict(), dict)
#     s1_new = Solution.from_dict(s1.as_dict())
#     assert s1_new.volume.magnitude == 2
