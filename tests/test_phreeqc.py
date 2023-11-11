"""
pyEQL volume and concentration methods test suite
=================================================

This file contains tests for the volume and concentration-related methods
used by pyEQL's Solution class
"""

import numpy as np
import pytest

from pyEQL import Solution
from pyEQL.engines import PhreeqcEOS


@pytest.fixture()
def s1():
    return Solution(volume="2 L", engine="phreeqc")


@pytest.fixture()
def s2():
    return Solution([["Na+", "4 mol/L"], ["Cl-", "4 mol/L"]], volume="2 L", engine="phreeqc")


@pytest.fixture()
def s3():
    return Solution([["Na+", "4 mol/kg"], ["Cl-", "4 mol/kg"]], volume="2 L", engine="phreeqc")


@pytest.fixture()
def s5():
    # 100 mg/L as CaCO3 ~ 1 mM
    return Solution([["Ca+2", "40.078 mg/L"], ["CO3-2", "60.0089 mg/L"]], volume="1 L", engine="phreeqc")


@pytest.fixture()
def s5_pH():
    # 100 mg/L as CaCO3 ~ 1 mM
    return Solution([["Ca+2", "40.078 mg/L"], ["CO3-2", "60.0089 mg/L"]], volume="1 L", balance_charge="pH")


@pytest.fixture()
def s6():
    # non-electroneutral solution with lots of hardness
    # alk = -118 meq/L * 50 = -5900 mg/L, hardness = 12*50 = 600 mg/L as CaCO3
    # charge balance = 2+10+10+10-120-20-12 = -120 meq/L
    return Solution(
        [
            ["Ca+2", "1 mM"],  # 2 meq/L
            ["Mg+2", "5 mM"],  # 10 meq/L
            ["Na+1", "10 mM"],  # 10 meq/L
            ["Ag+1", "10 mM"],  # no contribution to alk or hardness
            ["CO3-2", "6 mM"],  # no contribution to alk or hardness
            ["SO4-2", "60 mM"],  # -120 meq/L
            ["Br-", "20 mM"],
        ],  # -20 meq/L
        volume="1 L",
    )


@pytest.fixture()
def s6_Ca():
    # non-electroneutral solution with lots of hardness
    # alk = -118 meq/L * 50 = -5900 mg/L, hardness = 12*50 = 600 mg/L as CaCO3
    # charge balance = 2+10+10+10-120-20-12 = -120 meq/L
    return Solution(
        [
            ["Ca+2", "1 mM"],  # 2 meq/L
            ["Mg+2", "5 mM"],  # 10 meq/L
            ["Na+1", "10 mM"],  # 10 meq/L
            ["Ag+1", "10 mM"],  # no contribution to alk or hardness
            ["CO3-2", "6 mM"],  # no contribution to alk or hardness
            ["SO4-2", "60 mM"],  # -120 meq/L
            ["Br-", "20 mM"],
        ],  # -20 meq/L
        volume="1 L",
        balance_charge="Ca+2",
    )


def test_empty_solution_3():
    # create an empty solution
    s1 = Solution(database=None, engine="phreeqc")
    # It should return type Solution
    assert isinstance(s1, Solution)
    # It should have exactly 1L volume
    assert s1.volume.to("L").magnitude == 1.0
    #  the solvent should be water
    assert s1.solvent == "H2O(aq)"
    # It should have 0.997 kg water mass
    assert np.isclose(s1.solvent_mass.to("kg").magnitude, 0.997, atol=1e-3)
    # the temperature should be 25 degC
    assert s1.temperature.to("degC").magnitude == 25
    # the pressure should be 1 atm
    assert s1.pressure.to("atm").magnitude == 1
    # the pH should be 7.0
    assert np.isclose(s1.get_activity("H+"), 1e-7, atol=1e-9)
    # assert np.isclose(s1.pH, 7.0, atol=0.01)
    assert np.isclose(s1.pE, 8.5)
    # it should contain H2O, H+, and OH- species
    assert set(s1.components.keys()) == {"H2O(aq)", "OH[-1]", "H[+1]"}


def test_init_engines():
    """
    Test passing an EOS instance as well as the ideal and native EOS
    """
    s = Solution([["Na+", "4 mol/L"], ["Cl-", "4 mol/L"]], engine="phreeqc")
    assert isinstance(s.engine, PhreeqcEOS)
    assert s.get_activity_coefficient("Na+").magnitude * s.get_activity_coefficient("Cl-").magnitude < 1
    assert s.get_osmotic_coefficient().magnitude == 1
    # with pytest.warns(match="Solute Mg+2 not found"):
    assert s.get_activity_coefficient("Mg+2").magnitude == 1
    assert s.get_activity("Mg+2").magnitude == 0


def test_conductivity(s1, s2):
    # even an empty solution should have some conductivity
    assert s1.conductivity > 0
    # per CRC handbook "standard Kcl solutions for calibratinG conductiVity cells", 0.1m KCl has a conductivity of 12.824 mS/cm at 25 C
    s_kcl = Solution({"K+": "0.1 mol/kg", "Cl-": "0.1 mol/kg"})
    assert np.isclose(s_kcl.conductivity.magnitude, 1.2824, atol=0.02)  # conductivity is in S/m

    # TODO - expected failures due to limited temp adjustment of diffusion coeff
    s_kcl.temperature = "5 degC"
    assert np.isclose(s_kcl.conductivity.magnitude, 0.81837, atol=0.02)

    s_kcl.temperature = "50 degC"
    assert np.isclose(s_kcl.conductivity.magnitude, 1.91809, atol=0.05)

    # TODO - conductivity model not very accurate at high conc.
    s_kcl = Solution({"K+": "1 mol/kg", "Cl-": "1 mol/kg"})
    assert np.isclose(s_kcl.conductivity.magnitude, 10.862, rtol=0.2)


def test_equilibrate(s1, s2, s5_pH):
    assert "H2(aq)" not in s1.components
    orig_pH = s1.pH
    orig_pE = s1.pE
    s1.equilibrate()
    assert "H2(aq)" in s1.components
    assert np.isclose(s1.charge_balance, 0, atol=1e-7)
    assert np.isclose(s1.pH, orig_pH, atol=0.01)
    assert np.isclose(s1.pE, orig_pE)

    assert "NaOH(aq)" not in s2.components
    s2.equilibrate()
    orig_pH = s2.pH
    orig_pE = s2.pE
    orig_density = s2.density.magnitude
    orig_solv_mass = s2.solvent_mass.magnitude
    assert "NaOH(aq)" in s2.components

    # total element concentrations should be conserved after equilibrating
    assert np.isclose(s2.get_total_amount("Na", "mol").magnitude, 8)
    assert np.isclose(s2.get_total_amount("Cl", "mol").magnitude, 8)
    assert np.isclose(s2.solvent_mass.magnitude, orig_solv_mass)
    assert np.isclose(s2.density.magnitude, orig_density)
    assert np.isclose(s2.charge_balance, 0, atol=1e-7)
    assert np.isclose(s2.pH, orig_pH, atol=0.01)
    assert np.isclose(s2.pE, orig_pE)

    # this solution is the only one in the test that contains alkalinity
    # and equilibrating it results in a shift in the pH
    # the CO3-2 initially present reacts with the water to consume H+ and
    # increase the pH by approximately 0.0006 M (b/c at pH 7 virtually all
    # carbonate is present as HCO3-) -log10(0.001) =
    assert "HCO3[-1]" not in s5_pH.components
    assert np.isclose(s5_pH.charge_balance, 0)
    orig_pH = s5_pH.pH
    orig_pE = s5_pH.pE
    orig_density = s5_pH.density.magnitude
    orig_solv_mass = s5_pH.solvent_mass.magnitude
    set(s5_pH.components.keys())
    s5_pH.equilibrate()
    assert np.isclose(s5_pH.get_total_amount("Ca", "mol").magnitude, 0.001)
    assert np.isclose(s5_pH.get_total_amount("C(4)", "mol").magnitude, 0.001)
    # due to the large pH shift, water mass and density need not be perfectly conserved
    assert np.isclose(s5_pH.solvent_mass.magnitude, orig_solv_mass, atol=1e-3)
    assert np.isclose(s5_pH.density.magnitude, orig_density, atol=1e-3)
    assert np.isclose(s5_pH.charge_balance, 0)
    assert "CaOH[+1]" in s5_pH.components
    assert "HCO3[-1]" in s5_pH.components
    assert s5_pH.pH > orig_pH
    assert np.isclose(s5_pH.pE, orig_pE)
