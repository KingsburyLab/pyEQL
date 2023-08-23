"""
pyEQL volume and concentration methods test suite
=================================================

This file contains tests for the volume and concentration-related methods
used by pyEQL's Solution class
"""

import numpy as np
import pytest
from pyEQL import Solution, ureg
from pyEQL.engines import IdealEOS, NativeEOS


@pytest.fixture()
def s1():
    return Solution(volume="2 L")


@pytest.fixture()
def s2():
    return Solution([["Na+", "4 mol/L"], ["Cl-", "4 mol/L"]], volume="2 L")


@pytest.fixture()
def s3():
    return Solution([["Na+", "4 mol/kg"], ["Cl-", "4 mol/kg"]], volume="2 L")


@pytest.fixture()
def s4():
    return Solution([["Na+", "8 mol"], ["Cl-", "8 mol"]], volume="2 L")


@pytest.fixture()
def s5():
    # 100 mg/L as CaCO3
    return Solution([["Ca+2", "40 mg/L"], ["CO3-2", "60 mg/L"]], volume="1 L")


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


def test_empty_solution_3():
    # create an empty solution
    s1 = Solution(database=None)
    # It should return type Solution
    assert isinstance(s1, Solution)
    # It should have exactly 1L volume
    assert s1.volume.to("L").magnitude == 1.0
    #  the solvent should be water
    assert s1.solvent == "H2O(aq)"
    # It should have 0.997 kg water mass
    assert np.isclose(s1.solvent_mass.to("kg").magnitude, 0.9970415)
    # the temperature should be 25 degC
    assert s1.temperature.to("degC").magnitude == 25
    # the pressure should be 1 atm
    assert s1.pressure.to("atm").magnitude == 1
    # the pH should be 7.0
    assert np.isclose(s1.get_activity("H+"), 1e-7, atol=1e-9)
    assert np.isclose(s1.pH, 7.0, atol=0.01)
    assert np.isclose(s1.pE, 8.5)
    # it should contain H2O, H+, and OH- species
    assert set(s1.list_solutes()) == {"H2O(aq)", "OH[-1]", "H[+1]"}


def test_init_raises():
    with pytest.raises(ValueError, match="random is not a valid value"):
        Solution(engine="random")
    with pytest.raises(ValueError, match="Non-aqueous solvent detected"):
        Solution(solvent="D2O")
    with pytest.raises(ValueError, match="Multiple solvents"):
        Solution(solvent=["D2O", "MeOH"])


def test_init_engines():
    """
    Test passing an EOS instance as well as the ideal and native EOS
    """
    ideal = IdealEOS()
    s = Solution([["Na+", "4 mol/L"], ["Cl-", "4 mol/L"]], engine=ideal)
    assert s.engine == ideal
    s = Solution([["Na+", "4 mol/L"], ["Cl-", "4 mol/L"]], engine="ideal")
    assert isinstance(s.engine, IdealEOS)
    assert s.get_activity_coefficient("Na+").magnitude == 1
    assert s.get_osmotic_coefficient().magnitude == 1
    s = Solution([["Na+", "4 mol/L"], ["Cl-", "4 mol/L"]], engine="native")
    assert isinstance(s.engine, NativeEOS)
    assert s.get_activity_coefficient("Na+").magnitude < 1
    assert s.get_osmotic_coefficient().magnitude != 1


# create an empty and test solutions with the same volume using substance / volume,
# substance/mass, and substance units
def test_solute_addition(s2, s3, s4):
    # if solutes are added at creation-time with substance / volume units,
    # then the total volume of the solution should not change (should remain at 2 L)
    assert s2.volume.to("L").magnitude == 2

    # if solutes are added at creation-time with substance / volume units,
    # then the resulting mol/L concentrations should be exactly what was specified
    assert s2.get_amount("Na+", "mol/L").magnitude == 4

    # if solutes are added at creation-time with substance / mass units,
    # then the resulting mol/kg concentrations should be exactly what was specified
    assert s3.get_amount("Na+", "mol/kg").magnitude == 4

    # the water mass of solution s2 should be less than that of s3, because
    # of the volume recalculation
    result_molL = s2.solvent_mass.to("kg").magnitude
    result_molkg = s3.solvent_mass.to("kg").magnitude
    assert result_molL < result_molkg

    # if solutes are added at creation-time with substance units,
    # then the resulting mol amounts should be exactly what was specified
    assert s4.get_amount("Na+", "mol").magnitude == 8

    # the water mass of solution s2 should be less than that of s4, because
    # of the volume recalculation
    result_molL = s2.solvent_mass.to("kg").magnitude
    result_mol = s4.solvent_mass.to("kg").magnitude
    assert result_molL < result_mol


def test_chempot_energy(s1, s2):
    # TODO - double check calculation formula and update with quantitative tests
    # in an empty solution, chempot should be zero
    # assert np.isclose(s1.get_chemical_potential_energy().magnitude, 0)
    # assert np.isclose(s1.get_chemical_potential_energy(activity_correction = False).magnitude, 0)
    pass


def test_alkalinity_hardness_chargebalance(s3, s5, s6):
    assert np.isclose(s3.charge_balance, 0)
    assert np.isclose(s3.hardness, 0)
    assert np.isclose(s3.alkalinity, 0)

    assert np.isclose(s5.alkalinity.magnitude, 100, rtol=0.005)
    assert np.isclose(s5.hardness.magnitude, 100, rtol=0.005)
    assert np.isclose(s5.charge_balance, 0, atol=1e-5)

    assert np.isclose(s6.alkalinity.magnitude, -5900, rtol=0.005)
    assert np.isclose(s6.hardness.magnitude, 600, rtol=0.005)
    assert np.isclose(s6.charge_balance, -0.12)


def test_pressure_temperature(s5):
    orig_V = s5.volume
    s5.temperature = "50 degC"
    assert s5.temperature == ureg.Quantity("50 degC")
    assert s5.volume > orig_V
    intermediate_V = s5.volume
    s5.pressure = "2 atm"
    assert s5.pressure == ureg.Quantity("2 atm")
    assert s5.volume < intermediate_V


def test_get_el_amt_dict(s6):
    """ """
    water_mol = s6.components["H2O(aq)"]
    # scale volume to 8L
    s6 *= 8
    d = s6.get_el_amt_dict()
    for el, amt in zip(
        ["H(1)", "O(-2)", "Ca(2)", "Mg(2)", "Na(1)", "Ag(1)", "C(4)", "S(6)", "Br(-1)"],
        [water_mol * 2 * 8, (water_mol + 0.018 + 0.24) * 8, 0.008, 0.040, 0.08, 0.08, 0.048, 0.48, 0.16],
    ):
        assert np.isclose(d[el], amt, atol=1e-3)

    s = Solution({"Fe+2": "1 mM", "Fe+3": "5 mM", "FeCl2": "1 mM", "FeCl3": "5 mM"})
    d = s.get_el_amt_dict()
    for el, amt in zip(["Fe(2)", "Fe(3)", "Cl(-1)"], [0.002, 0.01, 0.002 + 0.015]):
        assert np.isclose(d[el], amt, atol=1e-3)


def test_p(s2):
    assert np.isclose(s2.p("Na+"), -1 * np.log10(s2.get_activity("Na+")))
    assert np.isclose(s2.p("Na+", activity=False), -1 * np.log10(s2.get_amount("Na+", "M").magnitude))
    assert np.isclose(s2.p("Mg++"), 0)


def test_get_amount(s3, s5):
    TEST_UNITS = [
        "mol/L",
        "mmol/L",
        "umol/L",
        "M",
        "eq/L",
        "meq/L",
        "ueq/L",
        "eq",
        "meq",
        "mol/kg",
        "mmol/kg",
        "umol/kg",
        "m",
        "fraction",
        "count",
        "%",
        "g/L",
        "mg/L",
        "ug/L",
        "ng/L",
        "g",
        "ng",
        "mg",
        "ug",
        "kg",
        "ppm",
        "ppb",
        "ppt",
        "mol",
        "eq",
    ]
    # TODO - make this test more precise i.e. test numerical values
    for u in TEST_UNITS:
        qty = s3.get_amount("Na+", u)
        assert isinstance(qty, ureg.Quantity), f"get_amount() failed for unit {u}"
        assert qty.magnitude > 0
    assert s3.get_amount("Na+", "ppm") == s3.get_amount("Na+", "mg/L")
    assert s3.get_amount("Na+", "ppb") == s3.get_amount("Na+", "ug/L")
    assert s3.get_amount("Na+", "ppt") == s3.get_amount("Na+", "ng/L")
    assert s3.get_amount("Na+", "eq/L") == s3.get_amount("Na+", "M")
    assert s3.get_amount("Na+", "meq/L") == s3.get_amount("Na+", "mmol/L")
    assert s5.get_amount("CO3-2", "eq/L") == -2 * s5.get_amount("CO3-2", "M")
    assert s5.get_amount("CO3-2", "eq") == -2 * s5.get_amount("CO3-2", "mol")
    # TODO - pint does not consider "mM" and "mmol/L" equivalent. Consider filing bug report? Or perhaps an issue with
    # my unit definition file
    # assert s3.get_amount('Na+', "mmol/L") == s3.get_amount('Na+', "mM")


def test_tds(s1, s2, s5):
    assert s1.total_dissolved_solids.magnitude == 0
    assert np.isclose(s2.total_dissolved_solids.magnitude, 4 * 58442.769)
    assert s2.total_dissolved_solids == s2.TDS
    assert np.isclose(s5.TDS.magnitude, 100)


def test_conductivity(s1, s2):
    # even an empty solution should have some conductivity
    assert s1.conductivity > 0
    # per CRC handbook "standard Kcl solutions for calibratinG conductiVity cells", 0.1m KCl has a conductivity of 12.824 mS/cm at 25 C
    s_kcl = Solution({"K+": "0.1 mol/kg", "Cl-": "0.1 mol/kg"})
    assert np.isclose(s_kcl.conductivity.magnitude, 1.2824, atol=0.02)  # conductivity is in S/m

    # TODO - expected failures due to limited temp adjustment of diffusion coeff
    # s_kcl.temperature = '5 degC'
    # assert np.isclose(s_kcl.conductivity.magnitude, 0.81837, atol=0.02)

    # s_kcl.temperature = '50 degC'
    # assert np.isclose(s_kcl.conductivity.magnitude, 1.91809, atol=0.02)

    # TODO - conductivity model not very accurate at high conc.
    s_kcl = Solution({"K+": "1 mol/kg", "Cl-": "1 mol/kg"})
    assert np.isclose(s_kcl.conductivity.magnitude, 10.862, rtol=0.2)


def test_arithmetic_and_copy(s2, s6):
    s6_scale = s6.copy()
    s6_scale *= 1.5
    assert s6_scale.volume == 1.5 * s6.volume
    assert s6_scale.pressure == s6.pressure
    for s, amt in s6_scale.components.items():
        assert amt == 1.5 * s6.components[s]
    s6_scale /= 2
    assert s6_scale.volume == 1.5 / 2 * s6.volume
    assert s6_scale.pressure == s6.pressure
    for s, amt in s6_scale.components.items():
        assert amt == 1.5 / 2 * s6.components[s]

    with pytest.raises(NotImplementedError):
        s6 - s6_scale

    # TODO - test pH and pE
    s2.temperature = "35 degC"
    s2.pressure = "1.1 atm"
    initial_mix_vol = s2.volume.to("L").magnitude + s6.volume.to("L").magnitude
    mix = s2 + s6
    assert isinstance(mix, Solution)

    assert mix.get_amount("Na+", "mol").magnitude == 8.01  # 4 M x 2 L + 10 mM x 1 L
    assert mix.get_amount("Na+", "mol").magnitude == 8.01
    assert mix.get_amount("Na+1", "mol").magnitude == 8.01
    assert mix.get_amount("Cl-", "mol").magnitude == 8.0
    assert mix.get_amount("Br-", "mol").magnitude == 0.02
    assert np.isclose(
        mix.volume.to("L").magnitude, initial_mix_vol, atol=0.15
    )  # 0.15 L tolerance; deviation is due to non-idealities
    assert np.isclose(
        mix.temperature.to("K").magnitude, (np.sum([(273.15 + 35) * 2, (273.15 + 25) * 1]) / initial_mix_vol), atol=1
    )  # 1 K tolerance
    assert np.isclose(
        mix.pressure.to("atm").magnitude, np.sum([1.1 * 2, 1 * 1]) / initial_mix_vol, atol=0.01
    )  # 0.01 atm tolerance

    s_bad = Solution()
    # workaround necessary b/c it's not currently possible to init a solution with a non-water solvent
    s_bad.solvent = "D2O"
    with pytest.raises(ValueError, match="Cannot add Solution with different solvents!"):
        s2 + s_bad

    s_bad = Solution(engine="ideal")
    with pytest.raises(ValueError, match="Cannot add Solution with different engines!"):
        s2 + s_bad

    s_bad = Solution()
    # bad workaround
    s_bad.database = "random_database.json"
    with pytest.raises(ValueError, match="Cannot add Solution with different databases!"):
        s2 + s_bad


def test_serialization(s1, s2):
    assert isinstance(s1.as_dict(), dict)
    s1_new = Solution.from_dict(s1.as_dict())
    assert s1_new.volume.magnitude == 2
    assert s1_new.components == s1.components
    assert np.isclose(s1_new.pH, s1.pH)
    assert np.isclose(s1_new.pE, s1.pE)
    assert s1_new.temperature == s1.temperature
    assert s1_new.pressure == s1.pressure
    assert s1_new.solvent == s1.solvent
    assert s1_new._engine == s1._engine
    # the solutions should point to different EOS instances
    assert s1_new.engine != s1.engine
    # also should point to different Store instances
    # TODO currently this test will fail due to a bug in maggma's __eq__
    # assert s1_new.database != s1.database

    s2_new = Solution.from_dict(s2.as_dict())
    assert s2_new.volume == s2.volume
    # components concentrations should be the same
    assert s2_new.components == s2.components
    # but not point to the same instances
    assert s2_new.components is not s2.components
    assert np.isclose(s2_new.pH, s2.pH)
    assert np.isclose(s2_new.pE, s2.pE)
    assert s2_new.temperature == s2.temperature
    assert s2_new.pressure == s2.pressure
    assert s2_new.solvent == s2.solvent
    assert s2_new._engine == s2._engine
    # the solutions should point to different EOS instances
    assert s2_new.engine != s2.engine
    # also should point to different Store instances
    # TODO currently this test will fail due to a bug in maggma's __eq__
    # assert s2_new.database != s2.database
