"""
pyEQL activity correction methods test suite
============================================

This file contains tests for some of the activity correction methods
employed by pyEQL.

NOTE: generally, these tests check the module output against experimental
data rather than the theoretical result of the respective functions. In some
cases, the output is also tested against a well-established model published
by USGS(PHREEQC)
"""
import platform

import numpy as np
import pytest

from pyEQL.activity_correction import _debye_parameter_activity, _debye_parameter_B
from pyEQL.solution import Solution

## Tests of the pitzer model


def test_units_and_equality():
    s1 = Solution([["Na+", "0.1 mol/L"], ["Cl-", "0.1 mol/L"]])

    # TODO - add a quantitative test for molal, molar, and rational conversion
    # the activity coefficient should be dimensionless
    assert s1.get_activity_coefficient("Na+").dimensionality == ""
    assert s1.get_activity_coefficient("Na+", scale="molar").dimensionality == ""
    assert s1.get_activity_coefficient("Na+", scale="rational").dimensionality == ""

    assert s1.get_activity("Na+").dimensionality == ""
    assert s1.get_activity("Na+", scale="molar").dimensionality == ""
    assert s1.get_activity("Na+", scale="rational").dimensionality == ""

    assert s1.get_osmotic_coefficient(scale="molal").dimensionality == ""
    assert s1.get_osmotic_coefficient(scale="rational").dimensionality == ""
    assert s1.get_osmotic_coefficient(scale="fugacity").dimensionality == ""
    with pytest.raises(ValueError, match="Invalid scale argument"):
        s1.get_activity_coefficient("Na+", scale="random")
    with pytest.raises(ValueError, match="Invalid scale argument"):
        s1.get_activity("Na+", scale="random")
    with pytest.raises(ValueError, match="Invalid scale argument"):
        s1.get_osmotic_coefficient(scale="random")

    # the activity should be dimensionless
    assert s1.get_activity("Na+").dimensionality == ""
    assert s1.get_activity("H2O").dimensionality == ""

    # the activity coefficient of both the Na+ and Cl- should be the same
    a1 = s1.get_activity_coefficient("Na+")
    a2 = s1.get_activity_coefficient("Cl-")
    assert a1 == a2


def test_debye_params():
    # tests of the various Debye Huckel parameters
    # A should be equal to 0.509 at 25 C, for log base 10
    assert np.isclose(_debye_parameter_activity().magnitude / 2.303, 0.509, atol=1e-3)
    assert np.isclose(_debye_parameter_B().to("nm**-1 * kg**0.5/mol**0.5").magnitude, 3.29, atol=1e-2)


def test_activity_crc_HCl():
    """
    calculate the activity coefficient of HCl at each concentration and compare
    to experimental data

    Experimental activity coefficient values at 25 degC are found in
    *CRC Handbook of Chemistry and Physics*, Mean Activity Coefficients of Electrolytes as a Function of Concentration,
    in: W.M. Haynes (Ed.), 92nd ed., 2011.

    """
    cation = "H+"
    nu_cation = 1
    anion = "Cl-"
    nu_anion = 1

    # list of concentrations to test, mol/kg
    conc_list = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5]

    # list of published experimental activity coefficients
    pub_activity_coeff = [
        0.965,
        0.952,
        0.929,
        0.905,
        0.876,
        0.832,
        0.797,
        0.768,
        0.759,
        0.811,
        1.009,
        2.380,
    ]

    for i, conc in enumerate(conc_list):
        conc_c = str(conc * nu_cation) + "mol/kg"
        conc_a = str(conc * nu_anion) + "mol/kg"
        sol = Solution()
        sol.add_solute(cation, conc_c)
        sol.add_solute(anion, conc_a)
        act_cat = sol.get_activity_coefficient(cation)
        act_an = sol.get_activity_coefficient(anion)
        result = (act_cat**nu_cation * act_an**nu_anion) ** (1 / (nu_cation + nu_anion))
        expected = pub_activity_coeff[i]

        assert np.isclose(result, expected, rtol=0.05)


def test_activity_crc_CsI():
    """
    calculate the activity coefficient of CsI at each concentration and compare
    to experimental data

    Experimental activity coefficient values at 25 degC are found in
    *CRC Handbook of Chemistry and Physics*, Mean Activity Coefficients of Electrolytes as a Function of Concentration,
    in: W.M. Haynes (Ed.), 92nd ed., 2011.

    """
    cation = "Cs+"
    nu_cation = 1
    anion = "I-"
    nu_anion = 1

    # list of concentrations to test, mol/kg
    conc_list = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2]

    # list of published experimental activity coefficients
    pub_activity_coeff = [
        0.965,
        0.951,
        0.925,
        0.898,
        0.863,
        0.804,
        0.749,
        0.688,
        0.601,
        0.534,
        0.470,
    ]

    for i, conc in enumerate(conc_list):
        conc_c = str(conc * nu_cation) + "mol/kg"
        conc_a = str(conc * nu_anion) + "mol/kg"
        sol = Solution()
        sol.add_solute(cation, conc_c)
        sol.add_solute(anion, conc_a)
        act_cat = sol.get_activity_coefficient(cation)
        act_an = sol.get_activity_coefficient(anion)
        result = (act_cat**nu_cation * act_an**nu_anion) ** (1 / (nu_cation + nu_anion))
        expected = pub_activity_coeff[i]

        assert np.isclose(result, expected, rtol=0.05)


def test_activity_crc_bacl2():
    """
    calculate the activity coefficient of BaCl2 at each concentration and compare
    to experimental data

    Experimental activity coefficient values at 25 degC are found in
    *CRC Handbook of Chemistry and Physics*, Mean Activity Coefficients of Electrolytes as a Function of Concentration,
    in: W.M. Haynes (Ed.), 92nd ed., 2011.

    """
    # list of concentrations to test, mol/kg
    conc_list = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1]

    # list of published experimental activity coefficients
    pub_activity_coeff = [
        0.887,
        0.849,
        0.782,
        0.721,
        0.653,
        0.559,
        0.492,
        0.436,
        0.391,
        0.393,
    ]

    for i, conc in enumerate(conc_list):
        conc_c = str(conc) + "mol/kg"
        conc_a = str(conc * 2) + "mol/kg"
        sol = Solution()
        sol.add_solute("Ba+2", conc_c)
        sol.add_solute("Cl-", conc_a)
        act_cat = sol.get_activity_coefficient("Ba+2")
        act_an = sol.get_activity_coefficient("Cl-")
        result = (act_cat**1 * act_an**2) ** (1 / 3)
        expected = pub_activity_coeff[i]

        assert np.isclose(result, expected, rtol=0.05)


def test_activity_crc_licl():
    """
    calculate the activity coefficient of LiCl at each concentration and compare
    to experimental data

    Experimental activity coefficient values at 25 degC are found in
    *CRC Handbook of Chemistry and Physics*, Mean Activity Coefficients of Electrolytes as a Function of Concentration,
    in: W.M. Haynes (Ed.), 92nd ed., 2011.

    """
    # list of concentrations to test, mol/kg
    conc_list = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5]

    # list of published experimental activity coefficients
    pub_activity_coeff = [
        0.965,
        0.952,
        0.928,
        0.904,
        0.874,
        0.827,
        0.789,
        0.756,
        0.739,
        0.775,
        0.924,
        2.0,
    ]

    for i, conc in enumerate(conc_list):
        conc_c = str(conc) + "mol/kg"
        conc_a = str(conc) + "mol/kg"
        sol = Solution()
        sol.add_solute("Li+", conc_c)
        sol.add_solute("Cl-", conc_a)
        act_cat = sol.get_activity_coefficient("Li+")
        act_an = sol.get_activity_coefficient("Cl-")
        result = (act_cat**1 * act_an**1) ** (1 / 2)
        expected = pub_activity_coeff[i]

        assert np.isclose(result, expected, rtol=0.05)


def test_activity_crc_rbcl():
    """
    calculate the activity coefficient of RbCl at each concentration and compare
    to experimental data

    Experimental activity coefficient values at 25 degC are found in
    *CRC Handbook of Chemistry and Physics*, Mean Activity Coefficients of Electrolytes as a Function of Concentration,
    in: W.M. Haynes (Ed.), 92nd ed., 2011.

    """
    # list of concentrations to test, mol/kg
    conc_list = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5]

    # list of published experimental activity coefficients
    pub_activity_coeff = [
        0.965,
        0.951,
        0.926,
        0.900,
        0.867,
        0.811,
        0.761,
        0.707,
        0.633,
        0.583,
        0.546,
        0.544,
    ]

    for i, conc in enumerate(conc_list):
        conc_c = str(conc) + "mol/kg"
        conc_a = str(conc) + "mol/kg"
        sol = Solution()
        sol.add_solute("Rb+", conc_c)
        sol.add_solute("Cl-", conc_a)
        result = sol.get_activity_coefficient("Rb+")
        expected = pub_activity_coeff[i]

        assert np.isclose(result, expected, rtol=0.05)


@pytest.mark.skipif(platform.machine() == "arm64" and platform.system() == "Darwin", reason="arm64 not supported")
def test_activity_crc_MgCl2():
    """
    calculate the activity coefficient of MgCl2 at each concentration and compare
    to experimental data

    Experimental activity coefficient values at 25 degC are found in
    *CRC Handbook of Chemistry and Physics*, Mean Activity Coefficients of Electrolytes as a Function of Concentration,
    in: W.M. Haynes (Ed.), 92nd ed., 2011.

    """
    cation = "Mg+2"
    nu_cation = 1
    anion = "Cl-"
    nu_anion = 2

    # list of concentrations to test, mol/kg
    conc_list = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5]

    # list of published experimental activity coefficients
    pub_activity_coeff = [
        0.889,
        0.852,
        0.790,
        0.734,
        0.672,
        0.590,
        0.535,
        0.493,
        0.485,
        0.577,
        1.065,
        14.40,
    ]

    for i, conc in enumerate(conc_list):
        conc_c = str(conc * nu_cation) + "mol/kg"
        conc_a = str(conc * nu_anion) + "mol/kg"
        sol = Solution()
        sol.add_solute(cation, conc_c)
        sol.add_solute(anion, conc_a)
        act_cat = sol.get_activity_coefficient(cation)
        act_an = sol.get_activity_coefficient(anion)
        result = (act_cat**nu_cation * act_an**nu_anion) ** (1 / (nu_cation + nu_anion))
        expected = pub_activity_coeff[i]
        print(sol.ionic_strength)

        assert np.isclose(result, expected, rtol=0.05)

        # ignore the highest concentration; I=8.35 when speciated vs. 15 without speciation
        if conc <= 2:
            sol.equilibrate()
            # NOTE: speciating the solution results in a decrease in the overall ionic strength, because some of the
            # Mg+2 is converted to monovalent complexes like MgOH+. Hence, the activity coefficients deviate a bit from
            # the published values.
            print(sol.ionic_strength)
            act_cat = sol.get_activity_coefficient(cation)
            act_an = sol.get_activity_coefficient(anion)
            result = (act_cat**nu_cation * act_an**nu_anion) ** (1 / (nu_cation + nu_anion))
            assert np.isclose(result, expected, rtol=0.05 * (1 + (4 * conc)))  # max rtol = 0.45


def test_activity_crc_KBr():
    """
    calculate the activity coefficient of KBr at each concentration and compare
    to experimental data

    Experimental activity coefficient values at 25 degC are found in
    *CRC Handbook of Chemistry and Physics*, Mean Activity Coefficients of Electrolytes as a Function of Concentration,
    in: W.M. Haynes (Ed.), 92nd ed., 2011.

    """
    cation = "K+"
    nu_cation = 1
    anion = "Br-"
    nu_anion = 1

    # list of concentrations to test, mol/kg
    conc_list = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5]

    # list of published experimental activity coefficients
    pub_activity_coeff = [
        0.965,
        0.952,
        0.927,
        0.902,
        0.870,
        0.817,
        0.771,
        0.722,
        0.658,
        0.617,
        0.593,
        0.626,
    ]

    for i, conc in enumerate(conc_list):
        conc_c = str(conc * nu_cation) + "mol/kg"
        conc_a = str(conc * nu_anion) + "mol/kg"
        sol = Solution()
        sol.add_solute(cation, conc_c)
        sol.add_solute(anion, conc_a)
        act_cat = sol.get_activity_coefficient(cation)
        act_an = sol.get_activity_coefficient(anion)
        result = (act_cat**nu_cation * act_an**nu_anion) ** (1 / (nu_cation + nu_anion))
        expected = pub_activity_coeff[i]

        assert np.isclose(result, expected, rtol=0.05)


def test_activity_crc_k2so4():
    """
    calculate the activity coefficient of K2SO4 at each concentration and compare
    to experimental data

    Experimental activity coefficient values at 25 degC are found in
    *CRC Handbook of Chemistry and Physics*, Mean Activity Coefficients of Electrolytes as a Function of Concentration,
    in: W.M. Haynes (Ed.), 92nd ed., 2011.

    """
    # list of concentrations to test, mol/kg
    conc_list = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]

    # list of published experimental activity coefficients
    pub_activity_coeff = [
        0.885,
        0.844,
        0.772,
        0.704,
        0.625,
        0.511,
        0.424,
        0.343,
        0.251,
    ]

    for i, conc in enumerate(conc_list):
        conc_c = str(conc * 2) + "mol/kg"
        conc_a = str(conc) + "mol/kg"
        sol = Solution()
        sol.add_solute("K+", conc_c)
        sol.add_solute("SO4-2", conc_a)
        act_cat = sol.get_activity_coefficient("K+")
        act_an = sol.get_activity_coefficient("SO4-2")
        result = (act_cat**2 * act_an**1) ** (1 / 3)
        expected = pub_activity_coeff[i]

        assert np.isclose(result, expected, rtol=0.05)


def test_activity_pitzer_nacl_1():
    """
    calculate the activity coefficient at each concentration and compare
    to experimental data

    Experimental activity coefficient values at 25 degC are found in
    *J. Phys. Chem. Reference Data* Vol 13 (1), 1984, p.53.

    """
    # list of concentrations to test, mol/kg
    conc_list = [0.1, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6]

    # list of published experimental activity coefficients
    pub_activity_coeff = [
        0.778,
        0.72,
        0.681,
        0.665,
        0.657,
        0.669,
        0.714,
        0.782,
        0.873,
        0.987,
    ]

    for i, conc in enumerate(conc_list):
        conc = str(conc) + "mol/kg"
        sol = Solution()
        sol.add_solute("Na+", conc)
        sol.add_solute("Cl-", conc)
        result = sol.get_activity_coefficient("Na+")
        expected = pub_activity_coeff[i]

        assert np.isclose(result, expected, rtol=0.05)


# The pitzer model diverges a bit from experimental data at high concentration
@pytest.mark.xfail()
def test_water_activity_pitzer_nacl_1():
    r"""
    calculate the water activity at each concentration and compare
    to experimental data

    Experimental osmotic coefficients for NaCl are found in:
    Pitzer and Pelper, 1984. "Thermodyamic Properties of Aqueous Sodium Chloride Solutions"
    *J. Phys. Chem. Ref. Data* 13(1).

    Osmotic coefficients were converted into water activity according to the equation
    found in
    Blandamer, Mike J., Engberts, Jan B. F. N., Gleeson, Peter T., Reis,
    Joao Carlos R., 2005. "Activity of water in aqueous systems: A frequently
    neglected property." *Chemical Society Review* 34, 440-458.

    .. math:: ln a_w = - \Phi M_w \sum_i m_i

    Where :math:`M_w` is the molar mass of water (0.018015 kg/mol) and :math:`m_i` is the molal concentration
    of each species.

    """
    # list of concentrations to test, mol/kg
    conc_list = [0.1, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6]

    # list of published experimental water activity
    pub_water_activity = [
        0.9964,
        0.9910,
        0.9821,
        0.9733,
        0.9646,
        0.9304,
        0.8975,
        0.8657,
        0.8351,
        0.8055,
    ]

    for i, conc in enumerate(conc_list):
        conc = str(conc) + "mol/kg"
        sol = Solution()
        sol.add_solute("Na+", conc)
        sol.add_solute("Cl-", conc)
        result = sol.get_water_activity()
        expected = pub_water_activity[i]

        assert np.isclose(result, expected, rtol=0.05)


def test_activity_pitzer_phreeqc_nacl_2():
    """
    calculate the activity coefficient at each concentration and compare
    to the output of the PHREEQC model

    PHREEQC version 3.1.4 was used to calculate density, conductivity, water
    activity, and NaCl activity coefficient for NaCl solutions up to 6m.
    The Pitzer model (pitzer.dat) database was used.
    <http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/>

    """
    # list of concentrations to test, mol/kg
    conc_list = [0.1, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6]

    # list of modeled activity coefficients
    phreeqc_pitzer_activity_coeff = [
        0.7771,
        0.7186,
        0.6799,
        0.6629,
        0.6561,
        0.6675,
        0.7128,
        0.7825,
        0.8729,
        0.9874,
    ]

    for i, conc in enumerate(conc_list):
        conc = str(conc) + "mol/kg"
        sol = Solution()
        sol.add_solute("Na+", conc)
        sol.add_solute("Cl-", conc)
        result = sol.get_activity_coefficient("Na+")
        expected = phreeqc_pitzer_activity_coeff[i]

        assert np.isclose(result, expected, rtol=0.05)


def test_water_activity_phreeqc_pitzer_nacl_2():
    """
    calculate the water activity at each concentration and compare
    to the output of the PHREEQC model

    PHREEQC version 3.1.4 was used to calculate density, conductivity, water
    activity, and NaCl activity coefficient for NaCl solutions up to 6m.
    The Pitzer model (pitzer.dat) database was used.
    <http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/>

    """
    # list of concentrations to test, mol/kg
    conc_list = [0.1, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6]

    # list of modeled water activities
    phreeqc_pitzer_water_activity = [
        0.997,
        0.992,
        0.984,
        0.975,
        0.967,
        0.932,
        0.893,
        0.851,
        0.807,
        0.759,
    ]

    for i, conc in enumerate(conc_list):
        sol = Solution(
            {
                "Na+": f"{conc} mol/kg",
                "Cl-": f"{conc} mol/kg",
            }
        )
        result = sol.get_water_activity()
        expected = phreeqc_pitzer_water_activity[i]

        assert np.isclose(result, expected, rtol=0.05)
        # to get pi in Pa, need V_w in m3/mol
        osmotic_pressure = -8.314 * 298.15 / 0.000018015 * np.log(result)
        assert np.isclose(sol.osmotic_pressure.to("Pa").magnitude, osmotic_pressure, rtol=0.05), f"{osmotic_pressure}"
