"""
pyEQL volume and concentration methods test suite
=================================================

This file contains tests for the volume and concentration-related methods
used by pyEQL's Solution class
"""

import numpy as np
import pytest

import pyEQL


@pytest.fixture()
def s1():
    return pyEQL.Solution(volume="2 L")


@pytest.fixture()
def s2():
    return pyEQL.Solution([["Na+", "4 mol/L"], ["Cl-", "4 mol/L"]], volume="2 L")


@pytest.fixture()
def s3():
    return pyEQL.Solution([["Na+", "4 mol/kg"], ["Cl-", "4 mol/kg"]], volume="2 L")


@pytest.fixture()
def s4():
    return pyEQL.Solution([["Na+", "8 mol"], ["Cl-", "8 mol"]], volume="2 L")


class Test_solute_addition:
    """
    test behavior of various methods for adding solutes to a solution
    -----------------------------------------------------------------

    """

    def test_set_amount_1(self, s2):
        """
        Tests for set_amount() method
        """
        # If the concentration of a solute is directly set with a substance / volume
        # unit, the volume should not change
        s2.set_amount("Na+", "5 mol/L")
        s2.set_amount("Cl-", "5 mol/L")
        assert s2.volume.to("L").magnitude == 2

    def test_set_amount_2(self, s2):
        # If the concentration of a solute is directly set with a substance / volume
        # unit, the water mass should be reduced
        original = s2.solvent_mass.to("kg").magnitude
        s2.set_amount("Na+", "5 mol/L")
        s2.set_amount("Cl-", "5 mol/L")
        assert s2.solvent_mass.to("kg").magnitude < original

    def test_set_amount_3(self, s2):
        # If the concentration of a solute is directly set with a substance / volume
        # unit, the resulting concentration should be exactly what was specified
        s2.set_amount("Na+", "5 mol/L")
        s2.set_amount("Cl-", "5 mol/L")
        assert np.allclose(s2.get_amount("Na+", "mol/L").magnitude, 5)

    def test_set_amount_4(self, s2):
        # If the concentration of a solute is directly set with a substance / mass
        # unit, the volume should increase
        original = s2.volume.to("L").magnitude
        s2.set_amount("Na+", "5 mol/kg")
        s2.set_amount("Cl-", "5 mol/kg")
        assert s2.volume.to("L").magnitude > original

    def test_set_amount_5(self, s2):
        # If the concentration of a solute is directly set with a substance / mass
        # unit, the water mass should not change
        original = s2.solvent_mass.to("kg").magnitude
        s2.set_amount("Na+", "5 mol/kg")
        s2.set_amount("Cl-", "5 mol/kg")
        assert np.allclose(s2.solvent_mass.to("kg").magnitude, original)

    def test_set_amount_6(self, s2):
        # If the concentration of a solute is directly set with a substance / mass
        # unit, the resulting concentration should be exactly what was specified
        s2.set_amount("Na+", "5 mol/kg")
        s2.set_amount("Cl-", "5 mol/kg")
        assert np.allclose(s2.get_amount("Na+", "mol/kg").magnitude, 5)

    def test_set_amount_7(self, s2):
        # If the concentration of a solute is directly set with a substance
        # unit, the volume should increase
        original = s2.volume.to("L").magnitude
        s2.set_amount("Na+", "10 mol")
        s2.set_amount("Cl-", "10 mol")
        assert s2.volume.to("L").magnitude > original

    def test_set_amount_8(self, s2):
        # If the concentration of a solute is directly set with a substance
        # unit, the water mass should not change
        original = s2.solvent_mass.to("kg").magnitude
        s2.set_amount("Na+", "10 mol")
        s2.set_amount("Cl-", "10 mol")
        assert np.allclose(s2.solvent_mass.to("kg").magnitude, original)

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
        assert np.allclose(s2.volume.to("L").magnitude, 2)

    def test_add_amount_2(self, s2):
        # If the concentration of a solute is directly increased with a substance / volume
        # unit, the water mass should be reduced
        original = s2.solvent_mass.to("kg").magnitude
        s2.add_amount("Na+", "1 mol/L")
        s2.add_amount("Cl-", "1 mol/L")
        assert s2.solvent_mass.to("kg").magnitude < original

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
        original = s3.volume.to("L").magnitude
        s3.add_amount("Na+", "1 mol/kg")
        s3.add_amount("Cl-", "1 mol/kg")
        assert s3.volume.to("L").magnitude > original

    def test_add_amount_5(self, s3):
        # If the concentration of a solute is directly increased with a substance / mass
        # unit, the water mass should not change
        original = s3.solvent_mass.to("kg").magnitude
        V_orig = s3.volume.to("L").magnitude
        s3.add_amount("Na+", "1 mol/kg")
        s3.add_amount("Cl-", "1 mol/kg")
        assert np.allclose(s3.solvent_mass.to("kg").magnitude, original)
        assert s3.volume.to("L").magnitude > V_orig
        # pH will be slightly higher than 7 b/c the addition of solute caused the
        # solution volume to increase, so lower mol/L = higher pH
        assert np.isclose(s3.pH, 7.0, atol=0.05)
        assert np.isclose(s3.pE, 8.5)

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
        original = s2.volume.to("L").magnitude
        s2.add_amount("Na+", "2 mol")
        s2.add_amount("Cl-", "2 mol")
        assert s2.volume.to("L").magnitude > original

    def test_add_amount_7(self, s2):
        # If the concentration of a solute is directly increased with a substance
        # unit, the water mass should not change
        original = s2.solvent_mass.to("kg").magnitude
        s2.add_amount("Na+", "2 mol")
        s2.add_amount("Cl-", "2 mol")
        assert np.allclose(s2.solvent_mass.to("kg").magnitude, original)

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
        original = s2.volume.to("L").magnitude
        s2.add_amount("Na+", "-2 mol")
        s2.add_amount("Cl-", "-2 mol")
        assert s2.volume.to("L").magnitude < original

    def test_add_amount_10(self, s2):
        # If the concentration of a solute is directly changed with a substance
        # unit, the water mass should not change
        original = s2.solvent_mass.to("kg").magnitude
        s2.add_amount("Na+", "-2 mol")
        s2.add_amount("Cl-", "-2 mol")
        assert np.allclose(s2.solvent_mass.to("kg").magnitude, original)

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

        # get_amount() - count
        assert np.isclose(s1.get_amount("Na+", "count").magnitude, 6.02214e23)
