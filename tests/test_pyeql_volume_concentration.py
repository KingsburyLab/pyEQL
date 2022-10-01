"""
pyEQL volume and concentration methods test suite
=================================================

This file contains tests for the volume and concentration-related methods
used by pyEQL's Solution class
"""

import numpy as np
import pytest

import pyEQL


@pytest.fixture
def s1():
    return pyEQL.Solution(volume="2 L")


@pytest.fixture
def s2():
    return pyEQL.Solution([["Na+", "4 mol/L"], ["Cl-", "4 mol/L"]], volume="2 L")


@pytest.fixture
def s3():
    return pyEQL.Solution([["Na+", "4 mol/kg"], ["Cl-", "4 mol/kg"]], volume="2 L")


@pytest.fixture
def s4():
    return pyEQL.Solution([["Na+", "8 mol"], ["Cl-", "8 mol"]], volume="2 L")


class Test_empty_solution:
    """
    test behavior when creating an empty solution
    ------------------------------------------------

    """

    def test_empty_solution_3(self):
        # create an empty solution
        s1 = pyEQL.Solution()
        # It should return type Solution
        assert isinstance(s1, pyEQL.solution.Solution)
        # It should have exactly 1L volume
        assert s1.get_volume().to("L").magnitude == 1.0
        #  the solvent should be water
        assert s1.get_solvent().get_name() == "H2O"
        # It should have 0.997 kg water mass
        assert np.isclose(s1.get_solvent_mass().to("kg").magnitude, 0.9970415)
        # the temperature should be 25 degC
        assert s1.temperature.to("degC").magnitude == 25
        # the pressure should be 1 atm
        assert s1.pressure.to("atm").magnitude == 1
        # the pH should be 7.0
        assert np.isclose(s1.get_activity("H+"), 1e-7, atol=1e-9)
        # it should contain H2O, H+, and OH- species
        assert set(s1.list_solutes()) == set(["H2O", "OH-", "H+"])


class Test_solute_addition:
    """
    test behavior of various methods for adding solutes to a solution
    -----------------------------------------------------------------

    """

    # create an empty and test solutions with the same volume using substance / volume,
    # substance/mass, and substance units
    def test_solute_addition(self, s2, s3, s4):

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

    def test_set_amount_1(self, s2):
        """
        Tests for set_amount() method
        """
        # If the concentration of a solute is directly set with a substance / volume
        # unit, the volume should not change
        s2.set_amount("Na+", "5 mol/L")
        s2.set_amount("Cl-", "5 mol/L")
        assert s2.get_volume().to("L").magnitude == 2

    def test_set_amount_2(self, s2):
        # If the concentration of a solute is directly set with a substance / volume
        # unit, the water mass should be reduced
        original = s2.get_solvent_mass().to("kg").magnitude
        s2.set_amount("Na+", "5 mol/L")
        s2.set_amount("Cl-", "5 mol/L")
        assert s2.get_solvent_mass().to("kg").magnitude < original

    def test_set_amount_3(self, s2):
        # If the concentration of a solute is directly set with a substance / volume
        # unit, the resulting concentration should be exactly what was specified
        s2.set_amount("Na+", "5 mol/L")
        s2.set_amount("Cl-", "5 mol/L")
        assert np.allclose(s2.get_amount("Na+", "mol/L").magnitude, 5)

    def test_set_amount_4(self, s2):
        # If the concentration of a solute is directly set with a substance / mass
        # unit, the volume should increase
        original = s2.get_volume().to("L").magnitude
        s2.set_amount("Na+", "5 mol/kg")
        s2.set_amount("Cl-", "5 mol/kg")
        assert s2.get_volume().to("L").magnitude > original

    def test_set_amount_5(self, s2):
        # If the concentration of a solute is directly set with a substance / mass
        # unit, the water mass should not change
        original = s2.get_solvent_mass().to("kg").magnitude
        s2.set_amount("Na+", "5 mol/kg")
        s2.set_amount("Cl-", "5 mol/kg")
        assert np.allclose(s2.get_solvent_mass().to("kg").magnitude, original)

    def test_set_amount_6(self, s2):
        # If the concentration of a solute is directly set with a substance / mass
        # unit, the resulting concentration should be exactly what was specified
        s2.set_amount("Na+", "5 mol/kg")
        s2.set_amount("Cl-", "5 mol/kg")
        assert np.allclose(s2.get_amount("Na+", "mol/kg").magnitude, 5)

    def test_set_amount_7(self, s2):
        # If the concentration of a solute is directly set with a substance
        # unit, the volume should increase
        original = s2.get_volume().to("L").magnitude
        s2.set_amount("Na+", "10 mol")
        s2.set_amount("Cl-", "10 mol")
        assert s2.get_volume().to("L").magnitude > original

    def test_set_amount_8(self, s2):
        # If the concentration of a solute is directly set with a substance
        # unit, the water mass should not change
        original = s2.get_solvent_mass().to("kg").magnitude
        s2.set_amount("Na+", "10 mol")
        s2.set_amount("Cl-", "10 mol")
        assert np.allclose(s2.get_solvent_mass().to("kg").magnitude, original)

    def test_set_amount_9(self, s2):
        # If the concentration of a solute is directly set with a substance / mass
        # unit, the resulting concentration should be exactly what was specified
        s2.set_amount("Na+", "10 mol")
        s2.set_amount("Cl-", "10 mol")
        assert np.allclose(s2.get_amount("Na+", "mol").magnitude, 10)

    def test_add_amount_1(self, s2):
        """
        Tests for add_amount() method
        """
        # substance / volume units
        # If the concentration of a solute is directly increased with a substance / volume
        # unit, the volume should not change
        s2.add_amount("Na+", "1 mol/L")
        s2.add_amount("Cl-", "1 mol/L")
        assert np.allclose(s2.get_volume().to("L").magnitude, 2)

    def test_add_amount_2(self, s2):
        # If the concentration of a solute is directly increased with a substance / volume
        # unit, the water mass should be reduced
        original = s2.get_solvent_mass().to("kg").magnitude
        s2.add_amount("Na+", "1 mol/L")
        s2.add_amount("Cl-", "1 mol/L")
        assert s2.get_solvent_mass().to("kg").magnitude < original

    def test_add_amount_3(self, s2):
        # If the concentration of a solute is directly increased with a substance / volume
        # unit, the resulting concentration should be exactly what was specified
        s2.add_amount("Na+", "1 mol/L")
        s2.add_amount("Cl-", "1 mol/L")
        assert np.allclose(s2.get_amount("Na+", "mol/L").magnitude, 5)

    def test_add_amount_4(self, s3):
        # substance / mass units

        # If the concentration of a solute is directly increased with a substance / mass
        # unit, the volume should increase
        original = s3.get_volume().to("L").magnitude
        s3.add_amount("Na+", "1 mol/kg")
        s3.add_amount("Cl-", "1 mol/kg")
        assert s3.get_volume().to("L").magnitude > original

    def test_add_amount_5(self, s3):
        # If the concentration of a solute is directly increased with a substance / mass
        # unit, the water mass should not change
        original = s3.get_solvent_mass().to("kg").magnitude
        s3.add_amount("Na+", "1 mol/kg")
        s3.add_amount("Cl-", "1 mol/kg")
        assert np.allclose(s3.get_solvent_mass().to("kg").magnitude, original)

    def test_add_amount_12(self, s3):
        # If the concentration of a solute is directly increased with a substance / mass
        # unit, the resulting concentration should be exactly what was specified
        s3.add_amount("Na+", "1 mol/kg")
        s3.add_amount("Cl-", "1 mol/kg")
        assert np.allclose(s3.get_amount("Na+", "mol/kg").magnitude, 5)

    def test_add_amount_6(self, s2):
        # substance units

        # If the concentration of a solute is directly increased with a substance
        # unit, the volume should increase
        original = s2.get_volume().to("L").magnitude
        s2.add_amount("Na+", "2 mol")
        s2.add_amount("Cl-", "2 mol")
        assert s2.get_volume().to("L").magnitude > original

    def test_add_amount_7(self, s2):
        # If the concentration of a solute is directly increased with a substance
        # unit, the water mass should not change
        original = s2.get_solvent_mass().to("kg").magnitude
        s2.add_amount("Na+", "2 mol")
        s2.add_amount("Cl-", "2 mol")
        assert np.allclose(s2.get_solvent_mass().to("kg").magnitude, original)

    def test_add_amount_8(self, s2):
        # If the concentration of a solute is directly increased with a substance
        # unit, the resulting concentration should be exactly what was specified
        s2.add_amount("Na+", "2 mol")
        s2.add_amount("Cl-", "2 mol")
        assert np.allclose(s2.get_amount("Na+", "mol").magnitude, 10)

    def test_add_amount_9(self, s2):
        # negative substance units
        # If the concentration of a solute is directly decreased with a substance
        # unit, the volume should decrease
        original = s2.get_volume().to("L").magnitude
        s2.add_amount("Na+", "-2 mol")
        s2.add_amount("Cl-", "-2 mol")
        assert s2.get_volume().to("L").magnitude < original

    def test_add_amount_10(self, s2):
        # If the concentration of a solute is directly changed with a substance
        # unit, the water mass should not change
        original = s2.get_solvent_mass().to("kg").magnitude
        s2.add_amount("Na+", "-2 mol")
        s2.add_amount("Cl-", "-2 mol")
        assert np.allclose(s2.get_solvent_mass().to("kg").magnitude, original)

    def test_add_amount_11(self, s2):
        # If the concentration of a solute is directly changed with a substance
        # unit, the resulting concentration should be exactly what was specified
        s2.add_amount("Na+", "-2 mol")
        s2.add_amount("Cl-", "-2 mol")
        assert np.allclose(s2.get_amount("Na+", "mol").magnitude, 6)


class Test_get_amount:
    """
    test the get_amount() method on a 1 mol/L NaCl solution
    ----------------------------
    1 mol NaCl / L = 58.44 g/L
    Na+ = 22.98977 g/mol

    """

    def test_get_amount_nacl(self):
        # create the 1 M NaCl solution
        s1 = pyEQL.Solution([["Na+", "1 mol/L"], ["Cl-", "1 mol/L"]])

        # get_amount() - mol/L
        assert np.isclose(s1.get_amount("Na+", "mol/L").magnitude, 1)

        # get_amount() - mol/kg
        assert np.isclose(s1.get_amount("Na+", "mol/kg").magnitude, 1.02181221888)

        # get_amount() - g/L
        assert np.isclose(s1.get_amount("Na+", "g/L").magnitude, 22.98977)

        # get_amount() - mg
        assert np.isclose(s1.get_amount("Na+", "mg").magnitude, 22989.77)

        # get_amount() - mol
        assert np.isclose(s1.get_amount("Na+", "mol").magnitude, 1)

        # get_amount() - fraction
        assert np.isclose(s1.get_amount("Na+", "fraction"), 0.01775457254)
