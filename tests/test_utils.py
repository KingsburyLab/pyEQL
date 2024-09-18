"""
Tests of pyEQL.utils module

"""

from iapws import IAPWS95, IAPWS97
from pytest import raises

from pyEQL import ureg
from pyEQL.utils import FormulaDict, create_water_substance, format_solutes_dict, standardize_formula


def test_standardize_formula():
    """
    Test formula standardization
    """
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
    assert standardize_formula("HCOO-") == "HCO2[-1]"
    assert standardize_formula("CO2-1") == "C2O4[-2]"
    assert standardize_formula("C2O4--") == "C2O4[-2]"
    assert standardize_formula("H3PO4") == "H3PO4(aq)"
    assert standardize_formula("H2SO4") == "H2SO4(aq)"
    assert standardize_formula("HClO4") == "HClO4(aq)"
    assert standardize_formula("CF3SO3-") == "CF3SO3[-1]"
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
    d.update({"Cl-": 0.5})
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
