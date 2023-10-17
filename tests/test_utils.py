"""
Tests of pyEQL.utils module

"""

from pyEQL.utils import FormulaDict, standardize_formula


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
