"""
pyEQL volume and concentration methods test suite
=================================================

This file contains tests for the volume and concentration-related methods
used by pyEQL's Solution class
"""

import copy
import logging
import platform
from importlib.resources import files
from itertools import zip_longest

import numpy as np
import pytest
import yaml
from monty.serialization import dumpfn, loadfn
from pint import Quantity

import pyEQL
import pyEQL.activity_correction as ac
from pyEQL import Solution, engines, ureg
from pyEQL.engines import IdealEOS, NativeEOS
from pyEQL.salt_ion_match import Salt
from pyEQL.solution import UNKNOWN_OXI_STATE


@pytest.fixture
def s1():
    return Solution(volume="2 L")


@pytest.fixture
def s2():
    return Solution([["Na+", "4 mol/L"], ["Cl-", "4 mol/L"]], volume="2 L")


@pytest.fixture
def s3():
    return Solution([["Na+", "4 mol/kg"], ["Cl-", "4 mol/kg"]], volume="2 L")


@pytest.fixture
def s4():
    return Solution([["Na+", "8 mol"], ["Cl-", "8 mol"]], volume="2 L")


@pytest.fixture
def s5():
    # 100 mg/L as CaCO3 ~ 1 mM
    return Solution([["Ca+2", "40.078 mg/L"], ["CO3-2", "60.0089 mg/L"]])


@pytest.fixture
def s5_pH():
    # 100 mg/L as CaCO3 ~ 1 mM
    return Solution([["Ca+2", "40.078 mg/L"], ["CO3-2", "60.0089 mg/L"]], balance_charge="pH")


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
    )


def test_empty_solution():
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
    assert set(s1.components.keys()) == {"H2O(aq)", "OH[-1]", "H[+1]"}
    assert np.isclose(s1.density.to("kg/m**3").magnitude, 997.0479, atol=0.1)
    assert np.isclose(s1.viscosity_kinematic.to("mm**2/s").magnitude, 0.8917, atol=1e-3)  # 1 cSt = 1 mm2/s
    assert np.isclose(s1.viscosity_dynamic, s1.viscosity_kinematic * s1.density, atol=1e-8)


@pytest.mark.skipif(platform.machine() == "arm64" and platform.system() == "Darwin", reason="arm64 not supported")
def test_oxi_state_handling():
    # see https://github.com/KingsburyLab/pyEQL/issues/116
    # and https://github.com/materialsproject/pymatgen/issues/3687
    s = Solution({"Na+": "0.5 mol/kg", "Br-": "0.5 mol/kg"}, pH=7, balance_charge="Br-")
    s.equilibrate()
    assert f"Br({UNKNOWN_OXI_STATE})" in s.get_components_by_element()


def test_diffusion_transport(s1, s2):
    # test ionic strength adjustment
    assert s1.get_diffusion_coefficient("H+") > s2.get_diffusion_coefficient("H+")

    # for Na+, d=122, a1=1.52, a2=3.7, A=1.173802/2.303 at 25 DegC, B = 3.2843078+10
    factor = np.exp(
        -1.52
        * 1.173802
        / 2.303
        * 1
        * np.sqrt(s2.ionic_strength.magnitude)
        / (1 + 3.2843078e10 * np.sqrt(s2.ionic_strength.magnitude) * 3.7 / (1 + s2.ionic_strength.magnitude**0.75))
    )
    assert np.isclose(
        factor * s2.get_diffusion_coefficient("Na+").magnitude,
        s2.get_diffusion_coefficient("Na+").magnitude,
        atol=5e-11,
    )
    s_dilute = Solution({"Na+": "1 mmol/L", "Cl-": "1 mmol/L"})
    assert np.isclose(
        s_dilute.get_diffusion_coefficient("Na+", activity_correction=False).magnitude, 1.334e-9, atol=1e-12
    )
    assert np.isclose(s_dilute.get_transport_number("Na+"), 0.396, atol=1e-3)
    assert np.isclose(s_dilute.get_transport_number("Cl-"), 0.604, atol=1e-3)

    # test setting a default value
    s2.default_diffusion_coeff = 0
    assert s2.get_diffusion_coefficient("Cs+").magnitude == 0
    s2.default_diffusion_coeff = 1e-9
    assert s2.get_diffusion_coefficient("Cs+", activity_correction=False).magnitude == 1e-9
    s2.default_diffusion_coeff = 0
    assert s2.get_diffusion_coefficient("Cs+", activity_correction=True).magnitude < 1e-9
    d25 = s2.get_diffusion_coefficient("Na+", activity_correction=False).magnitude
    nu25 = s2.water_substance.nu
    s2.temperature = "40 degC"
    d40 = s2.get_diffusion_coefficient("Na+", activity_correction=False).magnitude
    nu40 = s2.water_substance.nu
    assert np.isclose(
        d40,
        d25 * np.exp(122 / (273.15 + 40) - 122 / 298.15) * (nu25 / nu40),
        atol=5e-11,
    )

    # test correction factors for concentration, as per Appelo 2017 Fig 5
    D1 = Solution({"Na+": "1 umol/L", "Cl-": "1 umol/L"}).get_diffusion_coefficient("Na+").magnitude
    D2 = Solution({"Na+": "1.7 mol/kg", "Cl-": "1.7 mol/kg"}).get_diffusion_coefficient("Na+").magnitude
    assert np.isclose(D2 / D1, 0.54, atol=1e-2)

    D1 = Solution({"K+": "1 umol/L", "Cl-": "1 umol/L"}).get_diffusion_coefficient("K+").magnitude
    D2 = Solution({"K+": "0.5 mol/kg", "Cl-": "0.5 mol/kg"}).get_diffusion_coefficient("K+").magnitude
    assert np.isclose(D2 / D1, 0.80, atol=1e-2)


def test_init_raises(caplog):
    with pytest.raises(ValueError, match="random is not a valid value"):
        Solution(engine="random")
    with pytest.raises(ValueError, match="Non-aqueous solvent detected"):
        Solution(solvent="D2O")
    with pytest.raises(ValueError, match="Multiple solvents"):
        Solution(solvent=["D2O", "MeOH"])
    with pytest.raises(ValueError, match="Cannot specify both a non-default pH and H+"):
        Solution({"HCO3-": "1 mM", "CO3--": "1 mM", "H+": "1 mM"}, pH=10)
    module_log = logging.getLogger("pyEQL")
    with caplog.at_level(logging.WARNING, "pyEQL"):
        Solution({"HCO3-": "1 mM", "CO3--": "1 mM", "H+": "1 mM"}, log_level="warning")
        assert module_log.level == logging.WARNING
        assert "WARNING" in caplog.text


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


def test_component_subsets(s2):
    assert list(s2.cations.keys()) == ["Na[+1]", "H[+1]"]
    assert list(s2.anions.keys()) == ["Cl[-1]", "OH[-1]"]
    assert list(s2.neutrals.keys()) == ["H2O(aq)"]


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


@pytest.mark.skipif(platform.machine() == "arm64" and platform.system() == "Darwin", reason="arm64 not supported")
def test_charge_balance(s3, s5, s5_pH, s6, s6_Ca):
    assert np.isclose(s3.charge_balance, 0)
    assert np.isclose(s5.charge_balance, 0, atol=1e-5)
    assert np.isclose(s5_pH.charge_balance, 0, atol=1e-8)
    assert np.isclose(s6.charge_balance, -0.12)
    assert np.isclose(s6_Ca.charge_balance, 0, atol=1e-8)

    # test auto charge balance
    s = Solution(
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
        balance_charge="auto",
    )
    assert s.balance_charge == "auto"
    assert s._cb_species == "Na[+1]"
    assert np.isclose(s.charge_balance, 0, atol=1e-8)
    s.equilibrate()
    assert s.balance_charge == "auto"
    assert s._cb_species == "Na[+1]"
    assert np.isclose(s.charge_balance, 0, atol=1e-8)

    s = Solution({"Na+": "2 mM", "Cl-": "1 mM"}, balance_charge="auto")
    assert s.balance_charge == "auto"
    assert s._cb_species == "Cl[-1]"
    assert np.isclose(s.charge_balance, 0, atol=1e-8)
    s.equilibrate()
    assert s.balance_charge == "auto"
    assert s._cb_species == "Cl[-1]"
    assert np.isclose(s.charge_balance, 0, atol=1e-8)

    # check 'pH' when the solution needs to be made more POSITIVE
    s = Solution({"Na+": "2 mM", "Cl-": "1 mM"}, balance_charge="pH", pH=4)
    assert s.balance_charge == "pH"
    assert s._cb_species == "H[+1]"
    assert np.isclose(s.charge_balance, 0, atol=1e-8)
    assert s.pH > 4
    s.equilibrate()
    assert s.balance_charge == "pH"
    assert s._cb_species == "H[+1]"
    assert np.isclose(s.charge_balance, 0, atol=1e-8)

    # check 'pH' when the imbalance is extreme
    s = Solution({"Na+": "2 mM", "Cl-": "1 M"}, balance_charge="pH", pH=4)
    assert s.balance_charge == "pH"
    assert s._cb_species == "H[+1]"
    assert np.isclose(s.charge_balance, 0, atol=1e-8)
    assert np.isclose(s.pH, 0, atol=0.1)
    s.equilibrate()
    assert s.balance_charge == "pH"
    assert s._cb_species == "H[+1]"
    assert np.isclose(s.charge_balance, 0, atol=1e-8)

    # check warning when there isn't enough to balance
    s = Solution({"Na+": "1 M", "K+": "2 mM", "Cl-": "2 mM"}, balance_charge="K+")
    assert s.get_amount("K+", "mol/L") == 0

    # check "auto" with an electroneutral solution
    s = Solution({"Na+": "2 mM", "Cl-": "2 mM"}, balance_charge="auto")
    assert s.balance_charge == "auto"
    assert s._cb_species == "Na[+1]"
    assert np.isclose(s.charge_balance, 0, atol=1e-8)
    s.equilibrate()
    assert s.balance_charge == "auto"
    assert s._cb_species == "Na[+1]"
    assert np.isclose(s.charge_balance, 0, atol=1e-8)

    with pytest.raises(ValueError, match=r"Charge balancing species Zr\[\+4\] was not found"):
        s = Solution({"Na+": "2 mM", "Cl-": "2 mM"}, balance_charge="Zr[+4]")


def test_alkalinity_hardness(s3, s5, s6):
    assert np.isclose(s3.hardness, 0)
    assert np.isclose(s3.alkalinity, 0)

    assert np.isclose(s5.alkalinity.magnitude, 100, rtol=0.005)
    assert np.isclose(s5.hardness.magnitude, 100, rtol=0.005)

    assert np.isclose(s6.alkalinity.magnitude, -5900, rtol=0.005)
    assert np.isclose(s6.hardness.magnitude, 600, rtol=0.005)


def test_pressure_temperature(s5):
    orig_V = s5.volume
    s5.temperature = "50 degC"
    assert s5.temperature == ureg.Quantity("50 degC")
    assert s5.volume > orig_V
    intermediate_V = s5.volume
    s5.pressure = "2 atm"
    assert s5.pressure == ureg.Quantity("2 atm")
    assert s5.volume < intermediate_V


def test_elements(s5, s6):
    assert s6.elements == sorted({"Ag", "Br", "C", "Ca", "H", "Mg", "Na", "O", "S"})
    assert s6.chemical_system == "-".join(s6.elements)
    assert s5.chemical_system == "C-Ca-H-O"


def test_get_el_amt_dict(s6):
    """ """
    water_mol = s6.components["H2O(aq)"]
    # scale volume to 8L
    s6 *= 8
    d = s6.get_el_amt_dict()
    d_nested = s6.get_el_amt_dict(nested=True)
    for el, amt in zip(
        ["H(1.0)", "O(-2.0)", "Ca(2.0)", "Mg(2.0)", "Na(1.0)", "Ag(1.0)", "C(4.0)", "S(6.0)", "Br(-1.0)"],
        [water_mol * 2 * 8, (water_mol + 0.018 + 0.24) * 8, 0.008, 0.040, 0.08, 0.08, 0.048, 0.48, 0.16],
        strict=False,
    ):
        assert np.isclose(d[el], amt, atol=1e-3)

        el_no_valence = el.split("(")[0]
        valence = float(el.split("(")[1].split(")")[0])
        assert np.isclose(d_nested[el_no_valence][valence], amt, atol=1e-3)

    s = Solution({"Fe+2": "1 mM", "Fe+3": "5 mM", "FeCl2": "1 mM", "FeCl3": "5 mM"})
    d = s.get_el_amt_dict()
    d_nested = s.get_el_amt_dict(nested=True)
    for el, amt in zip(["Fe(2.0)", "Fe(3.0)", "Cl(-1.0)"], [0.002, 0.01, 0.002 + 0.015], strict=False):
        assert np.isclose(d[el], amt, atol=1e-3)

        el_no_valence = el.split("(")[0]
        valence = float(el.split("(")[1].split(")")[0])
        assert np.isclose(d_nested[el_no_valence][valence], amt, atol=1e-3)


def test_p(s2):
    assert np.isclose(s2.p("Na+"), -1 * np.log10(s2.get_activity("Na+")))
    assert np.isclose(s2.p("Na+", activity=False), -1 * np.log10(s2.get_amount("Na+", "M").magnitude))
    assert np.isnan(s2.p("Mg++"))


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


def test_components_by_element(s1, s2):
    assert s1.get_components_by_element() == {
        "H(1.0)": [
            "H2O(aq)",
            "OH[-1]",
            "H[+1]",
        ],
        "O(-2.0)": ["H2O(aq)", "OH[-1]"],
    }
    assert s2.get_components_by_element() == {
        "H(1.0)": [
            "H2O(aq)",
            "OH[-1]",
            "H[+1]",
        ],
        "O(-2.0)": ["H2O(aq)", "OH[-1]"],
        "Na(1.0)": ["Na[+1]"],
        "Cl(-1.0)": ["Cl[-1]"],
    }
    if platform.machine() == "arm64" and platform.system() == "Darwin":
        pytest.skip(reason="arm64 not supported")
    s2.equilibrate()
    assert s2.get_components_by_element() == {
        "H(1.0)": ["H2O(aq)", "OH[-1]", "H[+1]", "HCl(aq)", "NaOH(aq)", "HClO(aq)", "HClO2(aq)"],
        "H(0.0)": ["H2(aq)"],
        "O(-2.0)": [
            "H2O(aq)",
            "OH[-1]",
            "NaOH(aq)",
            "HClO(aq)",
            "ClO[-1]",
            "ClO2[-1]",
            "ClO3[-1]",
            "ClO4[-1]",
            "HClO2(aq)",
        ],
        "O(0.0)": ["O2(aq)"],
        "Na(1.0)": ["Na[+1]", "NaCl(aq)", "NaOH(aq)"],
        "Cl(-1.0)": ["Cl[-1]", "NaCl(aq)", "HCl(aq)"],
        "Cl(1.0)": ["HClO(aq)", "ClO[-1]"],
        "Cl(3.0)": ["ClO2[-1]", "HClO2(aq)"],
        "Cl(5.0)": ["ClO3[-1]"],
        "Cl(7.0)": ["ClO4[-1]"],
    }


def test_components_by_element_nested(s1, s2):
    assert s1.get_components_by_element(nested=True) == {
        "H": {
            1.0: ["H2O(aq)", "OH[-1]", "H[+1]"],
        },
        "O": {
            -2.0: ["H2O(aq)", "OH[-1]"],
        },
    }

    assert s2.get_components_by_element(nested=True) == {
        "H": {
            1.0: ["H2O(aq)", "OH[-1]", "H[+1]"],
        },
        "O": {
            -2.0: ["H2O(aq)", "OH[-1]"],
        },
        "Na": {
            1.0: ["Na[+1]"],
        },
        "Cl": {
            -1.0: ["Cl[-1]"],
        },
    }

    if platform.machine() == "arm64" and platform.system() == "Darwin":
        pytest.skip(reason="arm64 not supported")

    s2.equilibrate()

    assert s2.get_components_by_element(nested=True) == {
        "H": {
            1.0: [
                "H2O(aq)",
                "OH[-1]",
                "H[+1]",
                "HCl(aq)",
                "NaOH(aq)",
                "HClO(aq)",
                "HClO2(aq)",
            ],
            0.0: ["H2(aq)"],
        },
        "O": {
            -2.0: [
                "H2O(aq)",
                "OH[-1]",
                "NaOH(aq)",
                "HClO(aq)",
                "ClO[-1]",
                "ClO2[-1]",
                "ClO3[-1]",
                "ClO4[-1]",
                "HClO2(aq)",
            ],
            0.0: ["O2(aq)"],
        },
        "Na": {
            1.0: ["Na[+1]", "NaCl(aq)", "NaOH(aq)"],
        },
        "Cl": {
            -1.0: ["Cl[-1]", "NaCl(aq)", "HCl(aq)"],
            1.0: ["HClO(aq)", "ClO[-1]"],
            3.0: ["ClO2[-1]", "HClO2(aq)"],
            5.0: ["ClO3[-1]"],
            7.0: ["ClO4[-1]"],
        },
    }


def test_get_total_amount(s2):
    assert np.isclose(s2.get_total_amount("Na(1)", "g").magnitude, 8 * 58, 44)
    assert np.isclose(s2.get_total_amount("Na", "mol").magnitude, 8)
    assert np.isclose(s2.get_total_amount("N", "mol").magnitude, 0)
    assert np.isclose(s2.get_total_amount("Na", "ppm").magnitude, 4 * 23300, rtol=0.02)
    sox = Solution({"Fe+2": "10 mM", "Fe+3": "40 mM", "Cl-": "50 mM"}, pH=3)
    assert np.isclose(sox.get_total_amount("Fe(2)", "mol/L").magnitude, 0.01)
    assert np.isclose(sox.get_total_amount("Fe(3)", "mol/L").magnitude, 0.04)
    assert np.isclose(sox.get_total_amount("Fe", "mol").magnitude, 0.05)


@pytest.mark.skipif(platform.machine() == "arm64" and platform.system() == "Darwin", reason="arm64 not supported")
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
    assert np.isclose(s5_pH.get_total_amount("Ca", "mol").magnitude, 0.001, atol=1e-7)
    assert np.isclose(s5_pH.get_total_amount("C(4)", "mol").magnitude, 0.001, atol=1e-7)
    # due to the large pH shift, water mass and density need not be perfectly conserved
    assert np.isclose(s5_pH.solvent_mass.magnitude, orig_solv_mass, atol=1e-3)
    assert np.isclose(s5_pH.density.magnitude, orig_density, atol=1e-3)
    assert np.isclose(s5_pH.charge_balance, 0)
    assert "CaOH[+1]" in s5_pH.components
    assert "HCO3[-1]" in s5_pH.components
    assert s5_pH.pH > orig_pH
    assert np.isclose(s5_pH.pE, orig_pE)


def test_tds(s1, s2, s5):
    assert s1.total_dissolved_solids.magnitude == 0
    assert np.isclose(s2.total_dissolved_solids.magnitude, 4 * 58442.769)
    assert s2.total_dissolved_solids == s2.TDS
    assert np.isclose(s5.TDS.magnitude, 40.078 + 60.0089)


def test_conductivity(s1, s2):
    # per CRC handbook - "electrical conductiVity of Water" , conductivity of pure water
    # at 25 and 100 C is 0.0550 and 0.765 uS/cm
    assert np.isclose(s1.conductivity.to("uS/cm").magnitude, 0.055, atol=1e-3)

    # TODO - seems to be a possible bug related to setting temperature here
    # s1.temperature = "100 degC"
    # s2 = Solution(temperature='100 degC')
    # assert np.isclose(s1.conductivity.to('uS/cm').magnitude, 0.765, atol=1e-3)

    # CRC handbook table - "equivalent conductivity of electrolytes in aqueous solution"
    # nacl
    for conc, cond in zip([0.001, 0.05, 0.1], [123.68, 111.01, 106.69], strict=False):
        s1 = Solution({"Na+": f"{conc} mol/L", "Cl-": f"{conc} mol/L"})
        fail_msg = f"Conductivity test failed for NaCl at {conc} mol/L. Result = {s1.conductivity.to('S/m').magnitude}"
        assert np.isclose(s1.conductivity.to("S/m").magnitude, conc * cond / 10, atol=0.5), fail_msg

    # higher concentration data points from Appelo, 2017 Figure 4.
    s1 = Solution({"Na+": "2 mol/kg", "Cl-": "2 mol/kg"})
    assert np.isclose(s1.conductivity.to("mS/cm").magnitude, 145, atol=10)

    # MgCl2
    for conc, cond in zip([0.001, 0.05, 0.1], [124.15, 114.49, 97.05], strict=False):
        s1 = Solution({"Mg+2": f"{conc} mol/L", "Cl-": f"{2 * conc} mol/L"})
        fail_msg = f"Conductivity test failed for MgCl2 at {conc} mol/L. Result = {s1.conductivity.to('S/m').magnitude}"
        assert np.isclose(s1.conductivity.to("S/m").magnitude, 2 * conc * cond / 10, atol=1), fail_msg

    # per CRC handbook "standard KCl solutions for calibrating conductivity cells", 0.1m KCl has a conductivity of 12.824 mS/cm at 25 C
    s_kcl = Solution({"K+": "0.1 mol/kg", "Cl-": "0.1 mol/kg"})
    assert np.isclose(s_kcl.conductivity.magnitude, 1.2824, atol=0.25)  # conductivity is in S/m

    # TODO - expected failures due to limited temp adjustment of diffusion coeff
    s_kcl.temperature = "5 degC"
    assert np.isclose(s_kcl.conductivity.magnitude, 0.81837, atol=0.2)

    s_kcl.temperature = "50 degC"
    assert np.isclose(s_kcl.conductivity.magnitude, 1.91809, atol=0.2)

    # TODO - conductivity model not very accurate at high conc.
    s_kcl = Solution({"K+": "1 mol/kg", "Cl-": "1 mol/kg"})
    assert np.isclose(s_kcl.conductivity.magnitude, 10.862, atol=0.45)


@pytest.mark.skipif(platform.machine() == "arm64" and platform.system() == "Darwin", reason="arm64 not supported")
def test_arithmetic_and_copy(s2, s6):
    s6_scale = copy.deepcopy(s6)
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

    # confirm that a copied solution can be equilibrated
    s6_scale.equilibrate()

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


def test_as_from_dict(s1, s2):
    assert isinstance(s1.as_dict(), dict)
    s1_new = Solution.from_dict(s1.as_dict())
    assert s1_new.volume.magnitude == 2
    assert s1_new._solutes["H[+1]"] == "2e-07 mol"
    assert s1_new.get_total_moles_solute() == s1.get_total_moles_solute()
    assert s1_new.components == s1.components
    assert np.isclose(s1_new.pH, s1.pH)
    assert np.isclose(s1_new._pH, s1._pH)
    assert np.isclose(s1_new.pE, s1.pE)
    assert np.isclose(s1_new._pE, s1._pE)
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
    assert s2_new.get_total_moles_solute() == s2.get_total_moles_solute()
    assert np.isclose(s2_new.pH, s2.pH)
    assert np.isclose(s2_new._pH, s2._pH)
    assert np.isclose(s2_new.pE, s2.pE)
    assert np.isclose(s2_new._pE, s2._pE)
    assert s2_new.temperature == s2.temperature
    assert s2_new.pressure == s2.pressure
    assert s2_new.solvent == s2.solvent
    assert s2_new._engine == s2._engine
    # the solutions should point to different EOS instances
    assert s2_new.engine != s2.engine
    # also should point to different Store instances
    # TODO currently this test will fail due to a bug in maggma's __eq__
    # assert s2_new.database != s2.database


def test_serialization(s1, s2, tmp_path):
    dumpfn(s1, str(tmp_path / "s1.json"))
    s1_new = loadfn(str(tmp_path / "s1.json"))
    assert s1_new.volume.magnitude == 2
    assert s1_new._solutes["H[+1]"] == "2e-07 mol"
    assert s1_new.get_total_moles_solute() == s1.get_total_moles_solute()
    assert s1_new.components == s1.components
    assert np.isclose(s1_new.pH, s1.pH)
    assert np.isclose(s1_new._pH, s1._pH)
    assert np.isclose(s1_new.pE, s1.pE)
    assert np.isclose(s1_new._pE, s1._pE)
    assert s1_new.temperature == s1.temperature
    assert s1_new.pressure == s1.pressure
    assert s1_new.solvent == s1.solvent
    assert s1_new._engine == s1._engine
    # the solutions should point to different EOS instances
    assert s1_new.engine != s1.engine
    # also should point to different Store instances
    # TODO currently this test will fail due to a bug in maggma's __eq__
    # assert s1_new.database != s1.database

    dumpfn(s2, str(tmp_path / "s2.json"))
    s2_new = loadfn(str(tmp_path / "s2.json"))
    assert s2_new.volume == s2.volume
    # components concentrations should be the same
    assert s2_new.components == s2.components
    # but not point to the same instances
    assert s2_new.components is not s2.components
    assert s2_new.get_total_moles_solute() == s2.get_total_moles_solute()
    assert np.isclose(s2_new.pH, s2.pH)
    assert np.isclose(s2_new._pH, s2._pH)
    assert np.isclose(s2_new.pE, s2.pE)
    assert np.isclose(s2_new._pE, s2._pE)
    assert s2_new.temperature == s2.temperature
    assert s2_new.pressure == s2.pressure
    assert s2_new.solvent == s2.solvent
    assert s2_new._engine == s2._engine
    # the solutions should point to different EOS instances
    assert s2_new.engine != s2.engine
    # also should point to different Store instances
    # TODO currently this test will fail due to a bug in maggma's __eq__
    # assert s2_new.database != s2.database


def test_from_preset(tmp_path):
    preset_name = "seawater"
    solution = Solution.from_preset(preset_name)
    preset_path = files("pyEQL") / "presets" / "seawater.yaml"

    with open(str(preset_path)) as file:
        data = yaml.load(file, Loader=yaml.FullLoader)
    assert isinstance(solution, Solution)
    assert solution.temperature.to("degC") == ureg.Quantity(data["temperature"])
    assert solution.pressure == ureg.Quantity(data["pressure"])
    assert np.isclose(solution.pH, data["pH"], atol=0.01)
    for solute in solution._solutes:
        assert solute in data["solutes"]
    # test invalid preset
    with pytest.raises(FileNotFoundError):
        Solution.from_preset("nonexistent_preset")
    # test json as preset
    json_preset = tmp_path / "test.json"
    dumpfn(solution, json_preset)
    solution_json = Solution.from_preset(tmp_path / "test")
    assert isinstance(solution_json, Solution)
    assert solution_json.temperature.to("degC") == ureg.Quantity(data["temperature"])
    assert solution_json.pressure == ureg.Quantity(data["pressure"])
    assert np.isclose(solution_json.pH, data["pH"], atol=0.01)


def test_to_from_file(tmp_path):
    s1 = Solution(volume="2 L", pH=5)
    for f in ["test.json", "test.yaml"]:
        filename = tmp_path / f
        s1.to_file(filename)
        assert filename.exists()
        loaded_s1 = Solution.from_file(filename)
        assert isinstance(loaded_s1, Solution)
        assert pytest.approx(loaded_s1.volume.to("L").magnitude) == s1.volume.to("L").magnitude
    # test invalid extension raises error
    filename = tmp_path / "test_solution.txt"
    with pytest.raises(ValueError, match=r"File extension must be .json or .yaml"):
        s1.to_file(filename)
    with pytest.raises(FileNotFoundError, match=r"File .* not found!"):
        Solution.from_file(filename)


"""
The section below generates values to be used for test parametrization.
"""

_CATIONS = []
_ANIONS = []

for doc in pyEQL.IonDB.query(criteria={"size.molar_volume": {"$ne": None}}):
    if doc["charge"] > 0:
        _CATIONS.append(doc["formula"])
    elif doc["charge"] < 0:
        _ANIONS.append(doc["formula"])

_fill_value = _CATIONS[0] if len(_CATIONS) < len(_ANIONS) else _ANIONS[0]
# These cation-anion pairs include all ions for which molar volumes are available
_SOLUTES = list(zip_longest(_CATIONS, _ANIONS, fillvalue=_fill_value))
_FORMULAS_TO_SALTS = {f"{Salt(anion, cation).formula}(aq)": Salt(anion, cation) for anion, cation in _SOLUTES}
_criteria = {
    "model_parameters.molar_volume_pitzer.Beta0": {"$ne": None},
    "charge": 0.0,
    "formula": {"$in": list(_FORMULAS_TO_SALTS)},
}
# These salts include all salts for which Pitzer molar volume parameters are available
_SALTS = [_FORMULAS_TO_SALTS[doc["formula"]] for doc in pyEQL.IonDB.query(criteria=_criteria)]


@pytest.fixture(name="salt", params=_SOLUTES)
def fixture_salt(request: pytest.FixtureRequest) -> Salt:
    cation, anion = request.param
    return Salt(cation=cation, anion=anion)


def _get_solute_volume(
    ionic_strength: float,
    conc: Quantity,
    alphas: tuple[float, float],
    param: dict[str, dict[str, Quantity]],
    salt: Salt,
    temp: str,
) -> float:
    return ac.get_apparent_volume_pitzer(
        ionic_strength,
        conc,
        alphas[0],
        alphas[1],
        ureg.Quantity(param["Beta0"]["value"]).magnitude,
        ureg.Quantity(param["Beta1"]["value"]).magnitude,
        ureg.Quantity(param["Beta2"]["value"]).magnitude,
        ureg.Quantity(param["Cphi"]["value"]).magnitude,
        ureg.Quantity(param["V_o"]["value"]).magnitude,
        salt.z_cation,
        salt.z_anion,
        salt.nu_cation,
        salt.nu_anion,
        temp,
    )


class TestSolutionAdd:
    @staticmethod
    @pytest.fixture(name="salt_conc", params=[0.0, 1.0, 2.0])
    def fixture_salt_conc(request: pytest.FixtureRequest) -> float:
        return float(request.param)

    @staticmethod
    @pytest.fixture(name="engine", params=["ideal", "native", "phreeqc"])
    def fixture_engine(request: pytest.FixtureRequest) -> str:
        return str(request.param)

    @staticmethod
    @pytest.fixture(name="solution_sum")
    def fixture_solution_sum(solution: Solution) -> Solution:
        return solution + solution

    @staticmethod
    @pytest.mark.parametrize("engine", ["ideal"])
    def test_should_conserve_volume_with_ideal_engine(solution: Solution, solution_sum: Solution) -> None:
        assert np.isclose(solution_sum.volume.m, 2 * solution.volume.m)

    @staticmethod
    def test_should_preserve_engine_when_adding_solutions(solution: Solution, solution_sum: Solution) -> None:
        assert solution._engine == solution_sum._engine

    @staticmethod
    def test_should_preserve_the_number_of_moles_when_adding_solutions(
        solution: Solution, solution_sum: Solution
    ) -> None:
        moles_conserved = []
        for component, moles in solution.components.items():
            moles_conserved.append(solution_sum.components[component] == 2 * moles)
        assert all(moles_conserved)

    @staticmethod
    def test_should_add_all_components_to_new_solution(solution: Solution, solution_sum: Solution) -> None:
        assert sorted(solution.components) == sorted(solution_sum.components)

    @staticmethod
    def test_should_preserve_solution_solvent(solution: Solution, solution_sum: Solution) -> None:
        assert solution.solvent == solution_sum.solvent

    @staticmethod
    def test_should_preserve_solution_database(solution: Solution, solution_sum: Solution) -> None:
        assert solution.database == solution_sum.database

    @staticmethod
    @pytest.mark.parametrize("engine", ["native"])
    @pytest.mark.skipif(platform.machine() == "arm64" and platform.system() == "Darwin", reason="arm64 not supported")
    def test_should_replace_monatomic_species_from_engine(engine, caplog) -> None:
        # When initializing a solution without specifying the charge on the ion,
        # `.equilibrate()` should replace the ion with the ion with the charge
        # defined in the phreeqc database.
        solution = Solution({"Na": "1 mg/L"}, balance_charge="auto", engine=engine)
        assert "Na(aq)" in solution.components
        assert "Na[+1]" not in solution.components
        orig_el_amount = solution.get_total_amount("Na", "mol")

        with caplog.at_level(logging.INFO, "pyEQL"):
            solution.equilibrate()

        assert "amounts of species ['Na(aq)'] were not modified by PHREEQC" in caplog.text
        assert "Na[+1]" in solution.components  # correct charge assignment
        assert "Na(aq)" not in solution.components
        new_el_amount = solution.get_total_amount("Na", "mol")

        assert np.isclose(new_el_amount, orig_el_amount)

    @staticmethod
    @pytest.mark.parametrize("engine", ["native"])
    @pytest.mark.skipif(platform.machine() == "arm64" and platform.system() == "Darwin", reason="arm64 not supported")
    def test_should_replace_diatomic_species_from_engine(engine, caplog) -> None:
        # When initializing a solution by specifying the charge on the ion
        # that is different from the one determined by phreeqc,
        # `.equilibrate()` should replace the ion with the ion with the charge
        # determined by phreeqc.
        solution = Solution({"ReO4-2": "0.001 mg/L"}, balance_charge="auto", engine=engine)
        assert "ReO4[-2]" in solution.components
        assert "ReO4[-1]" not in solution.components
        orig_el_amount = solution.get_total_amount("Re", "mol")

        with caplog.at_level(logging.INFO, "pyEQL"):
            solution.equilibrate()

        # [ReO4-2] is not in phreeqc, but the element Re[+7] is, so it comes
        # up with ReO4[-1] as the species and replaces our incorrect ReO4[-2].
        assert "amounts of species ['ReO4[-2]'] were not modified by PHREEQC" in caplog.text
        assert "ReO4[-1]" in solution.components
        assert "ReO4[-2]" not in solution.components
        new_el_amount = solution.get_total_amount("Re", "mol")

        assert np.isclose(new_el_amount, orig_el_amount)

    @staticmethod
    @pytest.mark.parametrize("engine", ["native"])
    @pytest.mark.skipif(platform.machine() == "arm64" and platform.system() == "Darwin", reason="arm64 not supported")
    def test_should_not_discard_missing_species_from_engine(engine, caplog) -> None:
        # When initializing a solution by specifying a species with an element
        # that is not found in phreeqc, the species should not be discarded.
        solution = Solution({"Rh+3": "0.001 mg/L", "Rh2O3": "0.001 mg/L"}, balance_charge="auto", engine=engine)
        assert "Rh[+3]" in solution.components
        assert "Rh2O3(aq)" in solution.components
        orig_el_amount = solution.get_total_amount("Rh", "mol")

        with caplog.at_level(logging.INFO, "pyEQL"):
            solution.equilibrate()

        assert "amounts of species ['Rh2O3(aq)', 'Rh[+3]'] were not modified by PHREEQC" in caplog.text
        assert (
            "PHREEQC discarded element Rh during equilibration. Adding all components for this element." in caplog.text
        )
        assert "Rh[+3]" in solution.components  # still there
        assert "Rh2O3(aq)" in solution.components  # still there
        new_el_amount = solution.get_total_amount("Rh", "mol")

        assert np.isclose(new_el_amount, orig_el_amount)


class TestZeroSoluteVolume:
    @staticmethod
    @pytest.mark.parametrize("engine", ["ideal"])
    def test_should_return_zero_solute_volume_for_ideal_engine(solution: Solution) -> None:
        assert solution._get_solute_volume() == 0.0


class TestLinearCombinationSoluteVolume:
    @staticmethod
    @pytest.mark.parametrize(("salt_conc", "salt_conc_units"), [(1e-7, "mol/kg")])
    def test_should_return_solute_volume_equal_to_linear_combination_of_molar_solute_volumes_for_dilute_solutions(
        solution: Solution,
    ) -> None:
        sum_of_molar_volumes = ureg.Quantity(0.0, "L")

        for solute, component in solution.components.items():
            if solute != solution.solvent:
                molar_volume = solution.get_property(solute, "size.molar_volume")
                sum_of_molar_volumes += ureg.Quantity(component, "mol") * molar_volume

        assert solution._get_solute_volume().m == sum_of_molar_volumes.m

    @staticmethod
    @pytest.mark.parametrize("salt", _SALTS)
    def test_should_log_debug_message_when_using_pitzer_model(
        solution: Solution, salt: Salt, caplog: pytest.LogCaptureFixture
    ) -> None:
        caplog.set_level(logging.DEBUG, logger=solution.logger.name)
        _ = solution._get_solute_volume()
        expected_record = (
            engines.logger.name,
            logging.DEBUG,
            f"Updated solution volume using Pitzer model for solute {salt.formula}",
        )
        assert expected_record in caplog.record_tuples

    @staticmethod
    @pytest.mark.parametrize("salt", _SALTS)
    @pytest.mark.parametrize("salt_conc_units", ["mol/kg"])
    def test_should_use_major_salt_molar_volume_to_calculate_solute_volume_when_parameters_exist(
        solution: Solution, salt: Salt, salt_conc: float, salt_conc_units: str, alphas: tuple[float, float], volume: str
    ) -> None:
        param = solution.get_property(salt.formula, "model_parameters.molar_volume_pitzer")
        conc = ureg.Quantity(salt_conc, salt_conc_units)
        molality = (1 / 2) * (salt.nu_cation + salt.nu_anion) * conc
        expected_solute_volume = (
            _get_solute_volume(
                solution.ionic_strength,
                molality,
                alphas,
                param,
                salt,
                str(solution.temperature),
            )
            * conc
            * solution.solvent_mass
        ).to("L")
        solute_volume_without_protons_and_hydroxide = solution._get_solute_volume().to("L")

        if salt.cation != "H[+1]":
            solute_volume_without_protons_and_hydroxide -= solution.get_amount("H+", "mol") * solution.get_property(
                "H+", "size.molar_volume"
            )
        if salt.anion != "OH[-1]":
            solute_volume_without_protons_and_hydroxide -= solution.get_amount("OH-", "mol") * solution.get_property(
                "OH-", "size.molar_volume"
            )
        assert solute_volume_without_protons_and_hydroxide.m == expected_solute_volume.m
