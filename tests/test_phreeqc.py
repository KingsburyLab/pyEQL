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
    solution = Solution([], pH=7.00, temperature="25 degC", volume="1 L", engine="phreeqc")
    solution.equilibrate()
    # pH = -log10[H+]
    assert solution.get_amount("H+", "mol").magnitude == pytest.approx(1.0e-07 * 0.997, rel=0.01)
    # 14 - pH = log10[OH-]
    assert solution.get_amount("OH-", "mol").magnitude == pytest.approx(1.01e-07 * 0.997, rel=0.01)
    # small amount of H2 gas
    assert solution.get_amount("H2", "mol").magnitude == pytest.approx(0 * 0.997, rel=0.01)
    # 55.5 mol/kg * 0.997 kg = total moles of water
    assert solution.get_amount("H2O", "mol").magnitude == pytest.approx(55.5 * 0.997, rel=0.01)


def test_equilibrate_CO2_with_calcite():
    solution = Solution([], pH=7.0, volume="1 L", engine="phreeqc")
    solution.equilibrate(atmosphere=True, gases={"CO2": -2.95}, solids=["Calcite"])
    # 5 reactions: I) CaCO3 dissolution, II) Ka1, III) Ka2, IV) water dissociation, V) CaHCO3+ rxn in PHREEQC
    # 9 species, 5 components, 4 rxns exclude water dissociation
    assert solution.get_amount("Na+", "mol").magnitude == pytest.approx(0, rel=0.01)
    assert solution.get_amount("CO2(aq)", "mol").magnitude == pytest.approx(3.8e-05, rel=0.01)
    assert solution.get_amount("HCO3-", "mol").magnitude == pytest.approx(1.48e-03, rel=0.01)


def test_equilibrate_FeO3H3_ppt():
    # pint.errors.UndefinedUnitError: 'Fe' is not defined in the unit registry, hence using depreciated syntax
    solution = Solution([["Fe+3", "0.01 mol/L"], ["OH-", "10**-7 mol/L"]], volume="1 L", engine="phreeqc")
    solution.equilibrate()
    # Under supersaturation, Fe(OH)3 precipitates out with Ksp = 10**-38.8
    Fe3_act_coef = solution.get_activity_coefficient("Fe+3")
    OH_act_coef = solution.get_activity_coefficient("OH-")
    Fe3_conc = solution.get_amount("Fe+3", "mol/L").magnitude
    OH_conc = solution.get_amount("OH-", "mol/L").magnitude
    SI_FeO3H3 = np.log10((Fe3_act_coef.magnitude * Fe3_conc) * (OH_act_coef.magnitude * OH_conc) ** 3 / 10**-38.8)
    assert solution.engine.ppsol.si("Fe(OH)3(a)") > 0
    assert solution.engine.ppsol.si("Fe(OH)3(a)") == pytest.approx(SI_FeO3H3, rel=0.5)


def test_equilibrate_logC_pH_carbonate_3():
    # pint.errors.UndefinedUnitError: 'CO2' is not defined in the unit registry, hence using depreciated syntax
    solution = Solution([["CO2(aq)", "0.001 mol/L"]], pH=3.0, volume="1 L", engine="phreeqc")
    solution.equilibrate()
    # H2CO3 approx CO2(aq) dominant species at pH 3
    assert solution.get_amount("CO2(aq)", "mol").magnitude == pytest.approx(1.0e-3, rel=0.01)
    # CO2 + HCO3- + CO3-2 = total C = 0.001 mol with Ka1 = 10^6.3 and Ka2 = 10^-10.3
    # (10^-3 / 10^-6.3)[HCO3-] + [HCO3-] = total C = 0.001 mol
    assert solution.get_amount("HCO3-", "mol").magnitude == pytest.approx(5e-07, rel=0.1)
    # CO3-2 negligible at pH 3
    assert solution.get_amount("CO3-2", "mol").magnitude == pytest.approx(0, rel=0.01)


def test_equilibrate_logC_pH_carbonate_8_3():
    solution = Solution([["CO2(aq)", "0.001 mol/L"]], pH=8.3, volume="1 L", engine="phreeqc")
    solution.equilibrate()
    # CO2 + HCO3- + CO3-2 = total C = 0.001 mol with Ka1 = 10^6.3 and Ka2 = 10^-10.3
    # [H2CO3] + 100[H2CO3] + 2.5119e-3[H2CO3] = total C = 0.001 mol
    assert solution.get_amount("CO2(aq)", "mol").magnitude == pytest.approx(9.9e-6, rel=0.1)
    # HCO3- approx CO2(aq) dominant species at pH 7
    assert solution.get_amount("HCO3-", "mol").magnitude == pytest.approx(1.0e-3, rel=0.1)
    # CO2 + HCO3- + CO3-2 = total C = 0.001 mol with Ka1 = 10^6.3 and Ka2 = 10^-10.3
    # [H2CO3] + 100[H2CO3] + 5.0119e-9[H2CO3] = total C = 0.001 mol
    assert solution.get_amount("CO3-2", "mol").magnitude == pytest.approx(9.9e-6, rel=0.1)


def test_equilibrate_logC_pH_carbonate_13():
    solution = Solution([["CO2(aq)", "0.001 mol/L"]], pH=13.0, volume="1 L", engine="phreeqc")
    solution.equilibrate()
    # H2CO3 approx CO2(aq) is negligible at pH 13
    assert solution.get_amount("CO2(aq)", "mol").magnitude == pytest.approx(0, rel=0.01)
    # CO2 + HCO3- + CO3-2 = total C = 0.001 mol with Ka1 = 10^6.3 and Ka2 = 10^-10.3
    # 1.995e-7[HCO3-] + [HCO3-] + (10^-10.3 / 10^-13)[CO3-2] = total C = 0.001 mol
    assert solution.get_amount("HCO3-", "mol").magnitude == pytest.approx(1.991e-06, rel=0.5)  # flag loose tolerance*
    # CO2 + HCO3- + CO3-2 = total C = 0.001 mol with Ka1 = 10^6.3 and Ka2 = 10^-10.3
    # (10^-3*2 / 10^-16.6)[CO3-2] + (10^-3 / 10^-10.3)[CO3-2] + [CO3-2] = total C = 0.001 mol
    assert solution.get_amount("CO3-2", "mol").magnitude == pytest.approx(9.98e-4, rel=0.01)


def test_alkalinity():
    solution = Solution([["CO2(aq)", "0.001 mol/L"]], pH=7, volume="1 L", engine="phreeqc")
    solution.equilibrate()
    alk = solution.alkalinity
    # Total alkalinity
    HCO3 = solution.get_amount("HCO3-", "mol/L").magnitude
    CO3 = solution.get_amount("CO3-2", "mol/L").magnitude
    OH = solution.get_amount("OH-", "mol/L").magnitude
    H = solution.get_amount("H+", "mol/L").magnitude
    # Alkalinity calculated from the excess of negative charges from weak acids
    total_alk = HCO3 + 2 * CO3 + OH - H
    assert alk.to("mg/L").magnitude == pytest.approx(total_alk, abs=0.01)


def test_equilibrate_liquid():
    solution = Solution([["Cu+2", "4 mol/L"], ["O-2", "4 mol/L"]], volume="2 L", engine="phreeqc")
    solution.equilibrate()
    assert solution.get_total_amount("Cu", "mol").magnitude == pytest.approx(1.654340558452752)


def test_equilibrate_with_atm():
    solution = Solution([["Cu+2", "4 mol/L"], ["O-2", "4 mol/L"]], volume="2 L", engine="phreeqc")
    solution.equilibrate(atmosphere=True, solids=["Calcite"])
    # assert solution.get_total_amount("Cu", "mol").magnitude == pytest.approx(1.0422861990051868)
    assert solution.get_total_amount("N(0)", "mol").magnitude == pytest.approx(0)


#
# def test_equilibrate_with_n2_pp():
#     solution = Solution([["Cu+2", "4 mol/L"], ["O-2", "4 mol/L"]], volume="2 L", engine="phreeqc")
#     solution.equilibrate(gases={"N2": -0.1079}, solids=["Calcite"])
#     assert solution.get_total_amount("N(0)", "mol").magnitude == pytest.approx(0.0007084184487814338)
#
# def test_equilibrate_with_co2_pp_and_atm():
#     solution = Solution([["Cu+2", "4 mol/L"], ["O-2", "4 mol/L"]], volume="2 L", engine="phreeqc")
#     solution.equilibrate(atmosphere=True, gases={"CO2": -2})
#     assert solution.get_total_amount("Cu", "mol").magnitude == pytest.approx(2.0966131693539927)
#
# def test_equilibrate_with_co2_pp_atm():
#     solution = Solution([["Cu+2", "4 mol/L"], ["O-2", "4 mol/L"]], volume="2 L", engine="phreeqc")
#     solution.equilibrate(gases={"CO2": "0.01 atm"})
#     assert solution.get_total_amount("Cu", "mol").magnitude == pytest.approx(2.0966131693539927)
#
# def test_equilibrate_with_calcite():
#     solution = Solution([["Cu+2", "4 mol/L"], ["O-2", "4 mol/L"]], volume="2 L", engine="phreeqc")
#     solution.equilibrate(solids=["Calcite"])
#     assert solution.get_total_amount("Cu", "mol").magnitude == pytest.approx(1.6597326160055588)
#
# def test_equilibrate_with_calcite_and_atm():
#     solution = Solution([["Cu+2", "4 mol/L"], ["O-2", "4 mol/L"]], volume="2 L", engine="phreeqc")
#     solution.equilibrate(atmosphere=True, solids=["Calcite"])
#     assert solution.get_total_amount("Cu", "mol").magnitude == pytest.approx(1.0422861990051895)
#
# def test_equilibrate_unrecognized_component():
#     solution = Solution([["Cu+2", "4 mol/L"], ["O-2", "4 mol/L"]], volume="2 L", engine="phreeqc")
#     with pytest.raises(Exception):
#         solution.equilibrate(solids=["Ferroxite"])
