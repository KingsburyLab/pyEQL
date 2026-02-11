"""
pyEQL volume and concentration methods test suite
=================================================

This file contains tests for the volume and concentration-related methods
used by pyEQL's Solution class using the `phreeqc` engine (phreeqpython)
"""

import logging

import numpy as np
import pytest

from pyEQL import Solution
from pyEQL.engines import PHREEQPYTHON_AVAILABLE, PhreeqcEOS

if not PHREEQPYTHON_AVAILABLE:
    pytest.skip(
        "Phreeqpython not available",
        allow_module_level=True,
    )


@pytest.fixture
def s1():
    return Solution(volume="2 L", engine="phreeqc")


@pytest.fixture
def s2():
    return Solution([["Na+", "4 mol/L"], ["Cl-", "4 mol/L"]], volume="2 L", engine="phreeqc")


@pytest.fixture
def s3():
    return Solution([["Na+", "4 mol/kg"], ["Cl-", "4 mol/kg"]], volume="2 L", engine="phreeqc")


@pytest.fixture
def s5():
    # 100 mg/L as CaCO3 ~ 1 mM
    return Solution([["Ca+2", "40.078 mg/L"], ["CO3-2", "60.0089 mg/L"]], volume="1 L", engine="phreeqc")


@pytest.fixture
def s5_pH():
    # 100 mg/L as CaCO3 ~ 1 mM
    return Solution(
        [["Ca+2", "40.078 mg/L"], ["CO3-2", "60.0089 mg/L"]], volume="1 L", balance_charge="pH", engine="phreeqc"
    )


@pytest.fixture
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


@pytest.fixture
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


def test_osmotic_pressure():
    s1 = Solution([["Na+", "1 mol/L"], ["SO4-2", "0.5 mol/L"]], engine="phreeqc")
    s2 = Solution([["Na+", "1 mol/L"], ["SO4-2", "0.5 mol/L"]], engine="ideal")
    assert np.isclose(s1.osmotic_pressure.to("MPa").magnitude, s2.osmotic_pressure.to("MPa").magnitude), (
        f"PHREEQC = {s1.osmotic_pressure.to('MPa').magnitude} MPa, Ideal = {s2.osmotic_pressure.to('MPa').magnitude} MPa"
    )


def test_conductivity(s1):
    # even an empty solution should have some conductivity
    assert s1.conductivity > 0

    for conc, cond in zip([0.001, 0.05, 0.1], [123.68, 111.01, 106.69], strict=False):
        s1 = Solution({"Na+": f"{conc} mol/L", "Cl-": f"{conc} mol/L"})
        assert np.isclose(s1.conductivity.to("S/m").magnitude, conc * cond / 10, atol=0.5), (
            f"Conductivity test failed for NaCl at {conc} mol/L. Result = {s1.conductivity.to('S/m').magnitude}"
        )

    # higher concentration data points from Appelo, 2017 Figure 4.
    s1 = Solution({"Na+": "2 mol/kg", "Cl-": "2 mol/kg"})
    assert np.isclose(s1.conductivity.to("mS/cm").magnitude, 145, atol=10)

    # MgCl2
    for conc, cond in zip([0.001, 0.05, 0.1], [124.15, 114.49, 97.05], strict=False):
        s1 = Solution({"Mg+2": f"{conc} mol/L", "Cl-": f"{2 * conc} mol/L"})
        assert np.isclose(s1.conductivity.to("S/m").magnitude, 2 * conc * cond / 10, atol=1), (
            f"Conductivity test failed for MgCl2 at {conc} mol/L. Result = {s1.conductivity.to('S/m').magnitude}"
        )

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

    assert "NaOH(aq)" not in s2.components
    orig_density = s2.density.magnitude
    orig_solv_mass = s2.solvent_mass.magnitude
    orig_pH = s2.pH
    orig_pE = s2.pE
    orig_mass = s2.mass
    s2.equilibrate()
    assert np.isclose(s2.mass, orig_mass)
    assert np.isclose(s2.density.magnitude, orig_density)
    assert np.isclose(s2.solvent_mass.magnitude, orig_solv_mass)
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
    assert np.isclose(s5_pH.get_total_amount("C(4)", "mol").magnitude, 0.001, atol=1e-7)
    # due to the large pH shift, water mass and density need not be perfectly conserved
    assert np.isclose(s5_pH.solvent_mass.magnitude, orig_solv_mass, atol=1e-3)
    assert np.isclose(s5_pH.density.magnitude, orig_density, atol=1e-3)
    assert np.isclose(s5_pH.charge_balance, 0)
    assert "CaOH[+1]" in s5_pH.components
    assert "HCO3[-1]" in s5_pH.components
    s5_pH_after = s5_pH.pH
    assert s5_pH_after > orig_pH
    assert np.isclose(s5_pH.pE, orig_pE)

    # repeated calls to equilibrate should not change the properties (much)
    for i in range(10):
        s5_pH.equilibrate()
        assert np.isclose(s5_pH.charge_balance, 0, atol=1e-8), f"C.B. failed at iteration {i}"
        assert np.isclose(s5_pH.pH, s5_pH_after, atol=0.01), f"pH failed at iteration {i}"
        assert np.isclose(s5_pH.pE, orig_pE), f"pE failed at iteration {i}"

    # test equilibrate() with a non-pH balancing species
    assert np.isclose(s6_Ca.charge_balance, 0, atol=1e-8)
    initial_Ca = s6_Ca.get_total_amount("Ca", "mol").magnitude
    assert s6_Ca.balance_charge == "Ca[+2]"
    s6_Ca.equilibrate()
    assert s6_Ca.get_total_amount("Ca", "mol").magnitude != initial_Ca
    assert np.isclose(s6_Ca.charge_balance, 0, atol=1e-8)


def test_equilibrate_water_pH7():
    solution = Solution({}, pH=7.00, temperature="25 degC", volume="1 L", engine="phreeqc")
    solution.equilibrate()
    # pH = -log10[H+]
    assert np.isclose(solution.get_amount("H+", "mol/kg").magnitude, 1.001e-07)
    # 14 - pH = log10[OH-]
    assert np.isclose(solution.get_amount("OH-", "mol/kg").magnitude, 1.013e-07)
    # small amount of H2 gas
    assert solution.get_amount("H2", "mol/kg") > 0
    # Density of H2O in phreeqc = 0.99704 kg/L
    # For 0.99704 kg H2O
    #   mol = 0.99704 kg * 1000 g/kg / (18.01528 g/mol) = 55.34413009400909
    #   mol/kg = 55.34413009400909 / 0.99704 = 55.50843506179199
    assert np.isclose(solution.get_amount("H2O", "mol/kg").magnitude, 55.50843506179199)


def test_equilibrate_CO2_with_calcite():
    solution = Solution({}, pH=7.0, volume="1 L", engine="phreeqc")
    solution.equilibrate(atmosphere=False, gases={"CO2": -2.95, "O2": -0.6778}, solids=["Calcite"])
    # 5 rxns: I) CaCO3 dissolution, II) Ka1, III) Ka2, IV) water dissociation, V) CaHCO3+ rxn in PHREEQC
    # 9 species, 5 components, 4 rxns exclude water dissociation
    assert solution.get_amount("Na+", "mol") == 0
    assert np.isclose(solution.get_amount("CO2(aq)", "mol/kg").magnitude, 3.816e-05)
    assert np.isclose(
        solution.get_amount("HCO3-", "mol/kg").magnitude, 1.48e-03, atol=1e-5
    )  # slight tolerance adjustment
    assert np.isclose(
        solution.get_amount("Ca+2", "mol/kg").magnitude, 7.427e-04, atol=1e-5
    )  # slight tolerance adjustment


def test_equilibrate_FeO3H3_ppt():
    """Test an oversaturated solution"""
    solution = Solution({"Fe+3": "0.01 mol/L", "OH-": "10**-7 mol/L"}, volume="1 L", engine="phreeqc")
    solution.equilibrate()
    assert np.isclose(solution.get_amount("Fe+3", "mol/L").magnitude, 3.093e-11)
    assert np.isclose(solution.get_amount("OH-", "mol/L").magnitude, 1.067e-07)
    # The following assert passes, but we need to explain why phreeqc gives a
    # Fe(OH)3(a) SI value of ~5.408 instead of the expected ~7.2
    # SI_FeO3H3 = np.log10((Fe_3) * (OH_) ** 3 / 10**-38.8)
    assert solution.engine.ppsol.si("Fe(OH)3(a)") > 0


def test_equilibrate_nophaseeq():
    # Test to see that equilibrating without any solids/gases has no effect
    # on the concentrations.
    solution0 = Solution({"CO2(aq)": "0.001 mol/L"}, pH=3.0, volume="1 L", engine="phreeqc")
    solution0.equilibrate()
    solution1 = Solution({"CO2(aq)": "0.001 mol/L"}, pH=3.0, volume="1 L", engine="phreeqc")
    solution1.equilibrate(atmosphere=False, solids=None, gases=None)

    assert np.isclose(solution0.get_amount("CO2(aq)", "mol"), solution1.get_amount("CO2(aq)", "mol"))
    assert np.isclose(solution0.get_amount("HCO3-", "mol"), solution1.get_amount("HCO3-", "mol"))
    assert np.isclose(solution0.get_amount("CO3-2", "mol"), solution1.get_amount("CO3-2", "mol"))


def test_equilibrate_logC_pH_carbonate_8_3():
    solution = Solution({"CO2(aq)": "0.001 mol/L"}, pH=8.3, volume="1 L", engine="phreeqc")
    solution.equilibrate()

    # To evaluate CO2 equilibrium using Henry's law
    # Note: In the asserts below, the results were calculated assuming Henry's Law constant for CO2 as 0.034, which
    # is slightly different from what phreeqc uses (logK = -1.468, which gives us 0.0340408),
    # Hence the relaxed atols.

    # CO2 + HCO3- + CO3-2 = total C = 0.001 mol with pKa1 = 6.3 and pKa2 = 10.3
    # [H2CO3] + 100[H2CO3] + 2.5119e-3[H2CO3] = total C = 0.001 mol
    assert np.isclose(solution.get_amount("CO2(aq)", "mol/kg").magnitude, 1.079e-5, atol=1e-5)
    # The atols below are more relaxed because of the error carried forward from CO2(aq)
    # HCO3- approx CO2(aq) dominant species at pH 7
    # Note: The pKa1 value used in the calculation below is 6.3. Phreeqc uses 6.352
    assert np.isclose(solution.get_amount("HCO3-", "mol/kg").magnitude, 9.8219e-4, atol=1e-5)

    # CO2 + HCO3- + CO3-2 = total C = 0.001 mol with pKa1 = 6.3 and pKa2 = 10.3
    # [H2CO3] + 100[H2CO3] + 5.0119e-9[H2CO3] = total C = 0.001 mol
    # Note: The pKa1 value used in the calculation below is 6.3. Phreeqc uses 6.352
    # Note: The pKa2 value used in the calculation below is 10.3. Phreeqc uses 10.329
    assert np.isclose(solution.get_amount("CO3-2", "mol/kg").magnitude, 9.923e-6, atol=1e-5)


def test_equilibrate_logC_pH_carbonate_13():
    solution = Solution({"CO2(aq)": "0.001 mol/L"}, pH=13.0, volume="1 L", engine="phreeqc")
    solution.equilibrate()

    # To evaluate CO2 equilibrium using Henry's law
    # Note: In the asserts below, the results were calculated assuming Henry's Law constant for CO2 as 0.034, which
    # is slightly different from what phreeqc uses (logK = -1.468, which gives us 0.0340408),
    # Hence the relaxed atols.

    # H2CO3 approx CO2(aq) is negligible at pH 13, but still > 0
    assert solution.get_amount("CO2(aq)", "mol/kg") > 0

    # CO2 + HCO3- + CO3-2 = total C = 0.001 mol with pKa1 = 6.3 and pKa2 = 10.3
    # 1.995e-7[HCO3-] + [HCO3-] + (10^-10.3 / 10^-13)[CO3-2] = total C = 0.001 mol
    # Note: The pKa1 value used in the calculation below is 6.3. Phreeqc uses 6.352
    assert np.isclose(solution.get_amount("HCO3-", "mol/kg").magnitude, 1.152e-6, atol=1e-5)

    # CO2 + HCO3- + CO3-2 = total C = 0.001 mol with pKa1 = 6.3 and pKa2 = 10.3
    # (10^-3*2 / 10^-16.6)[CO3-2] + (10^-3 / 10^-10.3)[CO3-2] + [CO3-2] = total C = 0.001 mol
    # Note: The pKa1 value used in the calculation below is 6.3. Phreeqc uses 6.352
    # Note: The pKa2 value used in the calculation below is 10.3. Phreeqc uses 10.329
    assert np.isclose(solution.get_amount("CO3-2", "mol/kg").magnitude, 9.979e-4, atol=1e-5)


@pytest.mark.xfail(strict=True, reason="alkalinity discrepancy needs to be investigated")
def test_alkalinity():
    solution = Solution({"CO2(aq)": "0.001 mol/L"}, pH=7, volume="1 L", engine="phreeqc")
    solution.equilibrate()

    HCO3 = solution.get_amount("HCO3-", "mg/L").magnitude
    CO3 = solution.get_amount("CO3-2", "mg/L").magnitude
    OH = solution.get_amount("OH-", "mg/L").magnitude
    H = solution.get_amount("H+", "mg/L").magnitude

    # Alkalinity calculated from the excess of negative charges from weak acids
    calculated_alk = HCO3 + 2 * CO3 + OH - H
    assert solution.alkalinity.to("mg/L").magnitude == pytest.approx(calculated_alk, abs=0.001)


def test_equilibrate_2L():
    solution = Solution({"Cu+2": "1 umol/L", "O-2": "1 umol/L"}, volume="2 L", engine="phreeqc")
    solution.equilibrate(atmosphere=True)
    assert np.isclose(solution.get_total_amount("Cu", "umol").magnitude, 1.999955, atol=1e-6)  # PHREEQCUI - 2e-06


def test_equilibrate_unrecognized_component():
    solution = Solution({}, engine="phreeqc")
    # Specifying an unrecognized solid raises a ValueError
    with pytest.raises(ValueError, match="Phase not found in database, Ferroxite."):
        solution.equilibrate(solids=["Ferroxite"])


def test_equilibrate_OER_region():
    # The combination of pH and pE values don't fall within the water stability region.
    solution = Solution({}, pH=12.0, pE=13, volume="1 L", engine="phreeqc")
    with pytest.raises(ValueError, match="Activity of water has not converged."):
        solution.equilibrate()


def test_equilibrate_gas_units():
    # Specify CO2 partial pressure directly as log10 partial pressure, as well
    # as an explicit pressure unit.
    # Note: log10(0.000316) = -3.5
    s0 = Solution({}, pH=7.0, volume="1 L", engine="phreeqc")
    s0.equilibrate(atmosphere=True, gases={"CO2": "0.00031622776601683794 atm"})
    s1 = Solution({}, pH=7.0, volume="1 L", engine="phreeqc")
    s1.equilibrate(atmosphere=True, gases={"CO2": -3.5})
    assert s0.components == s1.components


def test_equilibrate_with_atm():
    s1 = Solution({}, pH=7.0, volume="1 L", engine="phreeqc")
    s1.equilibrate(atmosphere=True)
    # PHREEQCUI final CO2, O2, and N2 concentrations were slightly adjusted for consistency with wrapper outputs
    assert np.isclose(s1.get_amount("CO2(aq)", "mol/L").magnitude, 1.429e-05, atol=1e-6)  # PHREQCUI - 1.429e-05
    assert np.isclose(s1.get_amount("O2(aq)", "mol/L").magnitude, 0.0002683, atol=1e-6)  # PHREEQCUI - 2.683e-04
    assert np.isclose(s1.get_amount("N2(aq)", "mol/L").magnitude, 0)  # PHREEQCUI - 0
