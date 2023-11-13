"""
pyEQL solute properties test suite
============================================

This file contains tests of the Solution class methods that retrieve or
calculate properties of individual solutes. Methods currently included
in the testing are:

- get_transport_number()
- get_molar_conductivity()
- get_mobility()


"""

import numpy as np

from pyEQL import Solution, ureg

# relative tolerance between experimental and computed properties for this test file
RTOL = 0.05


class Test_transport_number:
    """
    test get_transport_number
    ------------------------------------------------
    """

    def test_transport_number(self):
        # create a test solution
        s1 = Solution([["K+", "0.1 mol/L"], ["Cl-", "0.1 mol/L"], ["FeO", "0.2 mol/L"]])

        # The transport number of any ion should be between 0 and 1
        assert s1.get_transport_number("K+") >= 0
        assert s1.get_transport_number("K+") <= 1

        # The transport number of water should be 0
        assert s1.get_transport_number("H2O") == 0

        # The transport number of any uncharged species should be 0
        assert s1.get_transport_number("FeO") == 0

        # The transport numbers should add up to 1
        total_t = 0
        for item in s1.components:
            total_t += s1.get_transport_number(item)

        assert round(abs(total_t.magnitude - 1), 7) == 0


class Test_molar_conductivity:
    """
    test get_molar_conductivity

    Published molar conductivity values at 25 degC and infinite dilution
    are found in CRC Handbook of Chemistry and Physics, 92nd ed., 2011.

    This reference is also the source of diffusion coefficient data
    in the default database, so the values molar conductivity (which
    is calculated from the diffusion coefficient) should match.
    ----------------------------------------------------------------
    """

    def test_molar_conductivity_potassium(self):
        # K+ - 73.48 x 10 ** -4 m ** 2 S / mol
        s1 = Solution([["K+", "0.001 mol/L"], ["Cl-", "0.001 mol/L"]], temperature="25 degC")
        result = s1.get_molar_conductivity("K+").to("m**2*S/mol").magnitude
        expected = ureg.Quantity("73.48e-4 m**2 * S / mol").magnitude

        assert np.isclose(result, expected, rtol=RTOL)

    def test_molar_conductivity_sodium(self):
        # Na+ - 50.08 x 10 ** -4 m ** 2 S / mol
        s1 = Solution([["Na+", "0.001 mol/L"], ["Cl-", "0.001 mol/L"]], temperature="25 degC")
        result = s1.get_molar_conductivity("Na+").to("m**2*S/mol").magnitude
        expected = ureg.Quantity("50.08e-4 m**2 * S / mol").magnitude

        assert np.isclose(result, expected, rtol=RTOL)

    def test_molar_conductivity_magnesium(self):
        # Mg+2 - 106 x 10 ** -4 m ** 2 S / mol
        s1 = Solution([["Mg+2", "0.001 mol/L"], ["Cl-", "0.002 mol/L"]], temperature="25 degC")
        result = s1.get_molar_conductivity("Mg+2").to("m**2*S/mol").magnitude
        expected = ureg.Quantity("106e-4 m**2 * S / mol").magnitude

        assert np.isclose(result, expected, atol=0.005)

    def test_molar_conductivity_chloride(self):
        # Cl- - 76.31 x 10 ** -4 m ** 2 S / mol
        s1 = Solution([["Na+", "0.001 mol/L"], ["Cl-", "0.001 mol/L"]], temperature="25 degC")
        result = s1.get_molar_conductivity("Cl-").to("m**2*S/mol").magnitude
        expected = ureg.Quantity("76.31e-4 m**2 * S / mol").magnitude

        assert np.isclose(result, expected, rtol=RTOL)

    def test_molar_conductivity_fluoride(self):
        # F- - 55.4 x 10 ** -4 m ** 2 S / mol
        s1 = Solution([["Na+", "0.001 mol/L"], ["F-", "0.001 mol/L"]], temperature="25 degC")
        result = s1.get_molar_conductivity("F-").to("m**2*S/mol").magnitude
        expected = ureg.Quantity("55.4e-4 m**2 * S / mol").magnitude

        assert np.isclose(result, expected, rtol=RTOL)

    def test_molar_conductivity_sulfate(self):
        # SO4-2 - 160 x 10 ** -4 m ** 2 S / mol
        s1 = Solution([["Na+", "0.002 mol/L"], ["SO4-2", "0.001 mol/L"]], temperature="25 degC")
        result = s1.get_molar_conductivity("SO4-2").to("m**2*S/mol").magnitude
        expected = ureg.Quantity("160.0e-4 m**2 * S / mol").magnitude

        assert np.isclose(result, expected, atol=0.002)

    def test_molar_conductivity_hydroxide(self):
        # OH- - 198 x 10 ** -4 m ** 2 S / mol
        s1 = Solution(temperature="25 degC")
        result = s1.get_molar_conductivity("OH-").to("m**2*S/mol").magnitude
        expected = ureg.Quantity("198e-4 m**2 * S / mol").magnitude

        assert np.isclose(result, expected, rtol=RTOL)

    def test_molar_conductivity_hydrogen(self):
        # H+ - 349.65 x 10 ** -4 m ** 2 S / mol
        s1 = Solution(temperature="25 degC")
        result = s1.get_molar_conductivity("H+").to("m**2*S/mol").magnitude
        expected = ureg.Quantity("349.65e-4 m**2 * S / mol").magnitude

        assert np.isclose(result, expected, rtol=RTOL)

    # molar conductivity of a neutral solute should be zero
    def test_molar_conductivity_neutral(self):
        s1 = Solution([["FeCl3", "0.001 mol/L"]], temperature="25 degC")
        result = s1.get_molar_conductivity("FeCl3").to("m**2*S/mol").magnitude
        expected = ureg.Quantity(0, "m**2 * S / mol").magnitude

        assert round(abs(result - expected), 5) == 0

    # molar conductivity of water should be zero
    def test_molar_conductivity_water(self):
        s1 = Solution(temperature="25 degC")
        result = s1.get_molar_conductivity("H2O").to("m**2*S/mol").magnitude
        expected = ureg.Quantity(0, "m**2 * S / mol").magnitude

        assert round(abs(result - expected), 5) == 0


class Test_mobility:
    """
    test get_mobility

    Since the mobility is simply the molar conductivity divided by
    Faraday's constant and the ion charge, and we have already
    tested get_molar_conductivity, here we only do a couple of tests
    to make sure get_mobility returns a value that is correctly
    scaled related to get_molar_conductivity
    ----------------------------------------------------------------
    """

    # relative error tolerance for assertWithinExperimentalError

    def test_mobility_potassium(self):
        s1 = Solution([["K+", "0.001 mol/L"], ["Cl-", "0.001 mol/L"]], temperature="25 degC")
        molar_conductivity = s1.get_molar_conductivity("K+").to("m**2*S/mol")
        expected = s1.get_mobility("K+").to("m**2/s/V").magnitude
        charge = s1.get_property("K+", "charge")
        # calculate the mobility from get_molar_conductivity, then compare with get_mobility
        result = (molar_conductivity / (ureg.N_A * ureg.e * abs(charge))).to("m**2/s/V").magnitude

        assert np.isclose(result, expected, rtol=RTOL)

    def test_mobility_chloride(self):
        s1 = Solution([["K+", "0.001 mol/L"], ["Cl-", "0.001 mol/L"]], temperature="25 degC")
        molar_conductivity = s1.get_molar_conductivity("Cl-").to("m**2*S/mol")
        expected = s1.get_mobility("Cl-").to("m**2/s/V").magnitude
        charge = s1.get_property("Cl-", "charge")
        # calculate the mobility from get_molar_conductivity, then compare with get_mobility
        result = (molar_conductivity / (ureg.N_A * ureg.e * abs(charge))).to("m**2/s/V").magnitude

        assert np.isclose(result, expected, rtol=RTOL)

    def test_mobility_magnesium(self):
        s1 = Solution([["Mg+2", "0.001 mol/L"], ["Cl-", "0.002 mol/L"]], temperature="25 degC")
        molar_conductivity = s1.get_molar_conductivity("Mg+2").to("m**2*S/mol")
        expected = s1.get_mobility("Mg+2").to("m**2/s/V").magnitude
        charge = s1.get_property("Mg+2", "charge")
        # calculate the mobility from get_molar_conductivity, then compare with get_mobility
        result = (molar_conductivity / (ureg.N_A * ureg.e * abs(charge))).to("m**2/s/V").magnitude

        assert np.isclose(result, expected, rtol=RTOL)

    def test_mobility_sulfate(self):
        s1 = Solution([["K+", "0.002 mol/L"], ["SO4-2", "0.001 mol/L"]], temperature="25 degC")
        molar_conductivity = s1.get_molar_conductivity("SO4-2").to("m**2*S/mol")
        expected = s1.get_mobility("SO4-2").to("m**2/s/V").magnitude
        charge = s1.get_property("SO4-2", "charge")
        # calculate the mobility from get_molar_conductivity, then compare with get_mobility
        result = (molar_conductivity / (ureg.N_A * ureg.e * abs(charge))).to("m**2/s/V").magnitude

        assert np.isclose(result, expected, rtol=RTOL)
