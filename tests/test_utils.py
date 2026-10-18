"""
Tests of pyEQL.utils module

"""

import numpy as np
import pytest
from iapws import IAPWS95, IAPWS97
from pytest import raises

from pyEQL import ureg
from pyEQL.utils import (
    FormulaDict,
    _translate_pint_quantity,
    create_water_substance,
    format_solutes_dict,
    standardize_formula,
)


def test_standardize_formula():
    """
    Test formula standardization
    """
    # test weird hyphens, super/subscripts
    assert standardize_formula("Mg⁺²") == "Mg[+2]"
    assert standardize_formula("VO²⁺") == "VO[+2]"
    assert standardize_formula("SO₄²⁻") == "SO4[-2]"
    for dash in ["‑", "‐", "‒", "–", "—", "−"]:  # noqa: RUF001
        assert standardize_formula(f"Cl{dash}") == "Cl[-1]"

    assert standardize_formula("H2O") == "H2O(aq)"
    assert standardize_formula("Na+") == "Na[+1]"
    assert standardize_formula("Na[+]") == "Na[+1]"
    assert standardize_formula("SO4--") == "SO4[-2]"
    assert standardize_formula("Mg+2") == "Mg[+2]"
    assert standardize_formula("O2") == "O2(aq)"
    assert standardize_formula("NH4+") == "NH4[+1]"
    assert standardize_formula("NH3") == "NH3(aq)"
    assert standardize_formula("HPO4--") == "HPO4[-2]"
    assert standardize_formula("H2PO4-") == "H2PO4[-1]"
    assert standardize_formula("SCN-") == "SCN[-1]"
    assert standardize_formula("I3-") == "I3[-1]"
    assert standardize_formula("N3-") == "N3[-1]"
    assert standardize_formula("P3-") == "P3[-1]"
    assert standardize_formula("Br3-") == "Br3[-1]"
    assert standardize_formula("HCOO-") == "HCO2[-1]"
    assert standardize_formula("CO2-1") == "CO2[-1]"
    assert standardize_formula("C2O4--") == "C2O4[-2]"
    assert standardize_formula("H3PO4") == "H3PO4(aq)"
    assert standardize_formula("H2SO4") == "H2SO4(aq)"
    assert standardize_formula("HClO4") == "HClO4(aq)"
    assert standardize_formula("CF3SO3-") == "CF3SO3[-1]"
    assert standardize_formula("S5-2") == "S5[-2]"
    assert standardize_formula("S5--") == "S5[-2]"
    assert standardize_formula("S[-0.4]") == "S5[-2]"
    assert standardize_formula("S[-0.5]") == "S4[-2]"
    assert standardize_formula("S[-0.66666667]") == "S3[-2]"
    # superscripts, subscripts, and permuted sign/charge number
    assert standardize_formula("PO₄³⁻") == "PO4[-3]"
    assert standardize_formula("Co²⁺") == "Co[+2]"
    # haloacetic acids
    assert standardize_formula("CCl3COO-") == "CCl3COO[-1]"
    assert standardize_formula("CF3COO-") == "CF3COO[-1]"
    assert standardize_formula("CI3COO-") == "CI3COO[-1]"
    assert standardize_formula("CBr3COO-") == "CBr3COO[-1]"
    # Cl+F
    assert standardize_formula("CCl2FCOO-") == "CFCl2COO[-1]"
    assert standardize_formula("CClF2COO-") == "CF2ClCOO[-1]"
    # Cl+I
    assert standardize_formula("CCl2ICOO-") == "CICl2COO[-1]"
    assert standardize_formula("CClI2COO-") == "CI2ClCOO[-1]"
    # Cl+Br
    assert standardize_formula("CBrCl2COO-") == "CBrCl2COO[-1]"
    assert standardize_formula("CBr2ClCOO-") == "CBr2ClCOO[-1]"
    assert standardize_formula("(NH4)2SO4") == "(NH4)2SO4(aq)"
    assert standardize_formula("NH4SO4-") == "NH4SO4[-1]"

    # Polymers should be left intact
    assert standardize_formula("(CO2)2") == "(CO2)2"
    assert standardize_formula("(H2O)3") == "(H2O)3"


@pytest.mark.parametrize(
    ("formula", "expected"),
    [
        # multimers keep their molecular formula (not reduced to an empirical formula)
        ("Hg2+2", "Hg2[+2]"),
        ("Hg2[+2]", "Hg2[+2]"),
        ("Hg2Cl2", "Hg2Cl2(aq)"),
        ("Fe2(OH)2+4", "Fe2(OH)2[+4]"),
        ("Cu2(OH)2+2", "Cu2(OH)2[+2]"),
        ("S2O4-2", "S2O4[-2]"),
        ("S4O6-2", "S4O6[-2]"),
        ("S2O8-2", "S2O8[-2]"),
        ("S2O6-2", "S2O6[-2]"),
        ("V10O28-6", "V10O28[-6]"),
        ("N2O4", "N2O4(aq)"),
        ("H2F2", "H2F2(aq)"),
        ("Br2", "Br2(aq)"),
        ("I2", "I2(aq)"),
        # legacy keys with fractional charges are restored to the smallest integer multiple
        ("I[-0.33333333]", "I3[-1]"),
        ("Br[-0.33333333]", "Br3[-1]"),
        ("S[-0.33333333]", "S3[-1]"),
        # ambiguous charge notation
        ("Ca2+", "Ca[+2]"),
        ("Fe3+", "Fe[+3]"),
        ("Hg2+", "Hg[+2]"),
        ("S2-", "S[-2]"),
        ("SO42-", "SO4[-2]"),
        ("CO32-", "CO3[-2]"),
        ("Cr2O72-", "Cr2O7[-2]"),
        ("UO22+", "UO2[+2]"),
        ("NH4+", "NH4[+1]"),
        ("NO3-", "NO3[-1]"),
        ("ClO3-", "ClO3[-1]"),
        ("PuO2+", "PuO2[+1]"),
        # PHREEQC valence-state labels
        ("[N-3]H4+", "NH4[+1]"),
        ("H2[S-2]", "H2S(aq)"),
        ("H3[As+3]O3", "As(OH)3(aq)"),
        ("[Fe+2]+2", "Fe[+2]"),
        # organic species are written in Hill order
        ("C2H6", "C2H6(aq)"),
        ("C2H4", "C2H4(aq)"),
        ("C6H12O6", "C6H12O6(aq)"),
        ("NaCH3COO", "C2H3NaO2(aq)"),
        ("Ba(CH3COO)2", "C4H6BaO4(aq)"),
        ("B(C6H5)4-", "C24H20B[-1]"),
        ("C6H5O7-3", "C6H5O7[-3]"),
        ("C8F15O2-", "C8F15O2[-1]"),
        ("CO(NH2)2", "CH4N2O(aq)"),
        # ... except for species with preferred names
        ("C2H3O2-", "CH3COO[-1]"),
        ("C2H4O2", "CH3COOH(aq)"),
        ("CH4O", "CH3OH(aq)"),
        ("HCOO-", "HCO2[-1]"),
        # inorganic carbon species are not written in Hill order
        ("HCO3-", "HCO3[-1]"),
        ("H2CO3", "H2CO3(aq)"),
        ("CaHCO3+", "CaHCO3[+1]"),
        ("Na2CO3", "Na2CO3(aq)"),
        ("CO2", "CO2(aq)"),
        ("Fe(CN)6-3", "Fe(CN)6[-3]"),
        ("HCN", "HCN(aq)"),
        ("KSCN", "KSCN(aq)"),
        ("OCN-", "CNO[-1]"),
        # fragments are written in conventional order inside complexes
        ("CaHPO4", "CaHPO4(aq)"),
        ("MgH2PO4+", "MgH2PO4[+1]"),
        ("NH4ClO4", "NH4ClO4(aq)"),
    ],
)
def test_standardize_formula_molecular(formula, expected):
    assert standardize_formula(formula) == expected
    # standardization is idempotent
    assert standardize_formula(expected) == expected


def test_formula_dict():
    """
    Make sure flexible key getting/setting works
    """
    d = FormulaDict({"Na+": 0.5, "Cl-1": 1.5, "Mg++": 0.5})
    # assert d == {"Na[+1]": 0.5, "Cl[-1]": 1.5, "Mg[+2]": 0.5}
    # sanitizing of keys should happen automagically
    # assert "Na+" in d
    del d["Mg+2"]
    assert "Mg[+2]" not in d
    assert not d.get("Mg[++]")
    d.update({"Cl-": np.float64(0.5)})
    print(d, d.data)
    assert d["Cl[-1]"] == 0.5
    assert d.get("Na+") == 0.5
    d.update({"Br-": 2})
    assert d["Br[-]"] == 2


def test_format_solute():
    """
    Test formatting solute dictionaries
    """
    test = {"Na+": 0.5, "Cl-": 0.5}
    base = {"Na+": "0.5 mol/kg", "Cl-": "0.5 mol/kg"}
    assert format_solutes_dict(test, units="mol/kg") == base

    bad_test = [["Na+", 0.5], ["Cl-", 0.5]]
    error_msg = "solute_dict must be a dictionary. Refer to the doc for proper formatting."
    with raises(TypeError, match=error_msg):
        format_solutes_dict(bad_test, units="mol/kg")


def test_create_water_substance():
    assert isinstance(create_water_substance(300, 0.1), IAPWS97)
    assert isinstance(create_water_substance(ureg.Quantity(-15, "degC"), 0.1), IAPWS95)


def test_translate_pint_quantity():
    import pytest  # noqa: PLC0415

    assert _translate_pint_quantity("5mol/kg") == (5, "mol/kg")
    assert _translate_pint_quantity("0.1g/L") == (0.1, "g/L")
    assert _translate_pint_quantity("1ppb") == (1, "ug/L")
    assert _translate_pint_quantity("0.01 ppm") == (0.01, "mg/L")
    assert _translate_pint_quantity("10**(-10.3) mol/L") == pytest.approx((5.011872336e-11, "mol/L"))
