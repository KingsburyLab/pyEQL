"""
pyEQL volume and concentration methods test suite
=================================================

This file contains tests for the volume and concentration-related methods
used by pyEQL's Solution class
"""

import logging
import platform

import numpy as np
import pytest

from pyEQL import Solution
from pyEQL.engines import PhreeqcEOS

if platform.machine().startswith("arm64") and platform.system().startswith("Darwin"):
    pytest.skip("skipping PHREEQC tests because arm64 is not supported", allow_module_level=True)


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
    return Solution(
        [["Ca+2", "40.078 mg/L"], ["CO3-2", "60.0089 mg/L"]], volume="1 L", balance_charge="pH", engine="phreeqc"
    )


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
        engine="phreeqc",
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
        engine="phreeqc",
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
    s.engine._destroy_ppsol()
    assert s.engine.ppsol is None


def test_conductivity(s1):
    # even an empty solution should have some conductivity
    assert s1.conductivity > 0

    for conc, cond in zip([0.001, 0.05, 0.1], [123.68, 111.01, 106.69]):
        s1 = Solution({"Na+": f"{conc} mol/L", "Cl-": f"{conc} mol/L"})
        assert np.isclose(
            s1.conductivity.to("S/m").magnitude, conc * cond / 10, atol=0.5
        ), f"Conductivity test failed for NaCl at {conc} mol/L. Result = {s1.conductivity.to('S/m').magnitude}"

    # higher concentration data points from Appelo, 2017 Figure 4.
    s1 = Solution({"Na+": "2 mol/kg", "Cl-": "2 mol/kg"})
    assert np.isclose(s1.conductivity.to("mS/cm").magnitude, 145, atol=10)

    # MgCl2
    for conc, cond in zip([0.001, 0.05, 0.1], [124.15, 114.49, 97.05]):
        s1 = Solution({"Mg+2": f"{conc} mol/L", "Cl-": f"{2*conc} mol/L"})
        assert np.isclose(
            s1.conductivity.to("S/m").magnitude, 2 * conc * cond / 10, atol=1
        ), f"Conductivity test failed for MgCl2 at {conc} mol/L. Result = {s1.conductivity.to('S/m').magnitude}"

    # per CRC handbook "standard KCl solutions for calibrating conductiVity cells",
    # 0.1m KCl has a conductivity of 12.824 mS/cm at 25 C
    s_kcl = Solution({"K+": "0.1 mol/kg", "Cl-": "0.1 mol/kg"})
    assert np.isclose(s_kcl.conductivity.magnitude, 1.2824, atol=0.02)  # conductivity is in S/m

    s_kcl.temperature = "5 degC"
    assert np.isclose(s_kcl.conductivity.magnitude, 0.81837, atol=0.06)

    s_kcl.temperature = "50 degC"
    assert np.isclose(s_kcl.conductivity.magnitude, 1.91809, atol=0.18)

    # TODO - conductivity model not very accurate at high conc.
    s_kcl = Solution({"K+": "1 mol/kg", "Cl-": "1 mol/kg"})
    assert np.isclose(s_kcl.conductivity.magnitude, 10.862, rtol=0.05)


def test_equilibrate(s1, s2, s5_pH, s6_Ca, caplog):
    assert "H2(aq)" not in s1.components
    orig_pH = s1.pH
    orig_pE = s1.pE
    orig_mass = s1.mass
    orig_density = s2.density.magnitude
    orig_solv_mass = s2.solvent_mass.magnitude
    s1.equilibrate()
    assert "H2(aq)" in s1.components
    assert np.isclose(s1.charge_balance, 0, atol=1e-8)
    assert np.isclose(s1.pH, orig_pH, atol=0.01)
    assert np.isclose(s1.pE, orig_pE)
    # assert np.isclose(s1.mass, orig_mass)
    # assert np.isclose(s1.density.magnitude, orig_density)
    # assert np.isclose(s1.solvent_mass.magnitude, orig_solv_mass)

    assert "NaOH(aq)" not in s2.components
    orig_density = s2.density.magnitude
    orig_solv_mass = s2.solvent_mass.magnitude
    orig_pH = s2.pH
    orig_pE = s2.pE
    orig_mass = s2.mass
    s2.equilibrate()
    assert "NaOH(aq)" in s2.components

    # total element concentrations should be conserved after equilibrating
    assert np.isclose(s2.get_total_amount("Na", "mol").magnitude, 8)
    assert np.isclose(s2.get_total_amount("Cl", "mol").magnitude, 8)
    assert np.isclose(s2.solvent_mass.magnitude, orig_solv_mass)
    assert np.isclose(s2.density.magnitude, orig_density)
    # the pH should drop due to hydrolysis of Ca+2
    assert s2.pH < orig_pH
    assert np.isclose(s2.pE, orig_pE)
    assert np.isclose(s2.mass, orig_mass)
    # this solution has balance_charge=None, therefore, the charge balance
    # may be off after equilibration
    assert not np.isclose(s2.charge_balance, 0, atol=1e-8)
    eq_Hplus = s2.components["H+"]
    s2.balance_charge = "pH"
    s2.equilibrate()
    assert np.isclose(s2.charge_balance, 0, atol=1e-8)
    assert s2.components["H+"] > eq_Hplus

    # test log message if there is a species not present in the phreeqc database
    s_zr = Solution(
        {"Zr+4": "0.05 mol/kg", "Na+": "0.05 mol/kg", "Cl-": "0.1 mol/kg"}, engine="phreeqc", log_level="WARNING"
    )
    totzr = s_zr.get_total_amount("Zr", "mol")
    with caplog.at_level(logging.WARNING, "pyEQL.engines"):
        s_zr.equilibrate()
        assert "likely absent from its database" in caplog.text
        assert "Zr[+4]" in s_zr.components
        assert s_zr.get_total_amount("Zr", "mol") == totzr

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

    # test equilibrate() with a non-pH balancing species
    assert np.isclose(s6_Ca.charge_balance, 0, atol=1e-8)
    initial_Ca = s6_Ca.get_total_amount("Ca", "mol").magnitude
    assert s6_Ca.balance_charge == "Ca[+2]"
    s6_Ca.equilibrate()
    assert s6_Ca.get_total_amount("Ca", "mol").magnitude != initial_Ca
    assert np.isclose(s6_Ca.charge_balance, 0, atol=1e-8)
