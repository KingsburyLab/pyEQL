"""
pyEQL density test suite
============================================

This file contains tests that check the density of electrolyte solutions
that are computed using the Pitzer model to get partial molar volumes.

NOTE: generally, these tests check the module output against experimental
data rather than the theoretical result of the respective functions. In some
cases, the output is also tested against a well-established model published
by USGS(PHREEQC)
"""

import numpy as np
import pyEQL


class Test_density_nacl:
    """
    test Pitzer model for density of NaCl
    ------------------------------------------------
    """

    def test_density_pitzer_nacl_1(self):
        """
        calculate the density at each concentration and compare
        to experimental data

        Published density values at 25 degC over a range of NaCl concentrations up
        to saturation are found in J. Phys. Chem. Reference Data Vol 13 (1), 1984, p.84

        """
        # list of concentrations to test, mol/kg
        conc_list = [0.1, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6]

        # list of published experimental densities
        pub_density = [
            1.00117,
            1.00722,
            1.01710,
            1.02676,
            1.03623,
            1.07228,
            1.10577,
            1.13705,
            1.16644,
            1.1942,
        ]

        for i, conc in enumerate(conc_list):
            conc = str(conc) + "mol/kg"
            sol = pyEQL.Solution()
            sol.add_solute("Na+", conc)
            sol.add_solute("Cl-", conc)
            result = sol.density.to("g/mL").magnitude
            expected = pub_density[i]

            assert np.isclose(result, expected, rtol=0.01)

    def test_density_pitzer_nacl_phreeqc_1(self):
        """
        calculate the density at each concentration and compare
        to the output of the PHREEQC model

        PHREEQC version 3.1.4 was used to calculate density, conductivity, water
        activity, and NaCl activity coefficient for NaCl solutions up to 6m.
        The Pitzer model (pitzer.dat) database was used.
        <http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/>

        """
        # list of concentrations to test, mol/kg
        conc_list = [0.1, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6]

        # list of modeled densities
        phreeqc_pitzer_density = [
            1.00115,
            1.00718,
            1.01702,
            1.02664,
            1.03606,
            1.07204,
            1.10562,
            1.13705,
            1.16651,
            1.19417,
        ]

        for i, conc in enumerate(conc_list):
            conc = str(conc) + "mol/kg"
            sol = pyEQL.Solution()
            sol.add_solute("Na+", conc)
            sol.add_solute("Cl-", conc)
            result = sol.density.to("g/mL").magnitude
            expected = phreeqc_pitzer_density[i]

            assert np.isclose(result, expected, rtol=0.01)
