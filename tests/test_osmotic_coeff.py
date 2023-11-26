"""
pyEQL osmotic coefficient test suite
============================================

This file contains tests for the osmotic coefficient method
employed by pyEQL.

NOTE: generally, these tests check the module output against experimental
data rather than the theoretical result of the respective functions.
"""

import numpy as np

from pyEQL import Solution


def test_osmotic_pressure():
    """
    The osmotic pressure of seawater is approximately 27 atm
    """
    # TODO - at present this test is inaccurate because in the complex matrix
    # of seawater, pyEQL falls back to using an ideal solution model with
    # unit osmotic coefficient.
    empty = Solution()
    assert np.isclose(empty.osmotic_pressure.to("atm").magnitude, 0, atol=1e-5)
    sea = Solution.from_preset("seawater")
    assert np.isclose(sea.osmotic_pressure.to("atm").magnitude, 27, rtol=0.15)


class Test_osmotic_pitzer:
    """
    test osmotic coefficient based on the Pitzer model
    ------------------------------------------------

    """

    def test_dimensionality(self):
        s1 = Solution([["Na+", "0.1 mol/L"], ["Cl-", "0.1 mol/L"]])
        assert s1.get_osmotic_coefficient().dimensionality == ""
        assert s1.get_osmotic_coefficient() >= 0

    def test_osmotic_pitzer_ammoniumnitrate(self):
        """
        calculate the osmotic coefficient at each concentration and compare
        to experimental data for ammonium nitrate

        References
        ----------
        May, P. M., Rowland, D., Hefter, G., & Königsberger, E. (2011).
        A Generic and Updatable Pitzer Characterization of Aqueous Binary Electrolyte Solutions at 1 bar and 25 °C.
        Journal of Chemical & Engineering Data, 56(12), 5066-5077. doi:10.1021/je2009329

        """
        # list of concentrations to test, mol/kg
        conc_list = [0.25, 0.5, 0.75, 1, 1.5, 2]

        # list of published experimental osmotic coefficients
        pub_osmotic_coeff = [0.86, 0.855, 0.83, 0.825, 0.80, 0.78]

        for i, conc in enumerate(conc_list):
            conc = str(conc) + "mol/kg"
            sol = Solution()
            sol.add_solute("NH4+", conc)
            sol.add_solute("NO3-", conc)
            result = sol.get_osmotic_coefficient()
            expected = pub_osmotic_coeff[i]

            assert np.isclose(result, expected, rtol=0.05)

    def test_osmotic_pitzer_coppersulfate(self):
        """
        calculate the osmotic coefficient at each concentration and compare
        to experimental data for copper sulate

        References
        ----------
        May, P. M., Rowland, D., Hefter, G., & Königsberger, E. (2011).
        A Generic and Updatable Pitzer Characterization of Aqueous Binary Electrolyte Solutions at 1 bar and 25 °C.
        Journal of Chemical & Engineering Data, 56(12), 5066-5077. doi:10.1021/je2009329

        """
        # list of concentrations to test, mol/kg
        conc_list = [0.25, 0.5, 0.75, 1]

        # list of published experimental osmotic coefficients
        pub_osmotic_coeff = [0.5, 0.485, 0.48, 0.485, 0.5]

        for i, conc in enumerate(conc_list):
            conc = str(conc) + "mol/kg"
            sol = Solution()
            sol.add_solute("Cu+2", conc)
            sol.add_solute("SO4-2", conc)
            result = sol.get_osmotic_coefficient()
            expected = pub_osmotic_coeff[i]

            assert np.isclose(result, expected, rtol=0.05)
