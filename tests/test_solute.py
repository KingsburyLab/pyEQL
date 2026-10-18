"""
Tests for the solute.py module
"""

from pyEQL.solute import Solute


def test_from_formula():
    s = Solute.from_formula("Mg+2")
    assert s.formula == "Mg[+2]"
    assert s.formula_pretty == "Mg^+2"
    assert s.formula_html == "Mg<sup>+2</sup>"
    assert s.formula_hill == "Mg"
    assert s.formula_latex == "Mg$^{+2}$"
    assert s.chemsys == "Mg"
    assert s.charge == 2
    assert s.n_atoms == 1
    assert s.n_elements == 1
    assert s.oxi_state_guesses == {"Mg": 2.0}
    # test behavior when oxi_state_guesses fails
    assert Solute.from_formula("Na[+2]").oxi_state_guesses == {}
    assert s.molecular_weight == "24.305 g/mol"
    s2 = Solute.from_formula("O3")
    assert s2.formula == "O3(aq)"
    assert s2.molecular_weight == "47.9982 g/mol"
    assert s2.oxi_state_guesses == {"O": 0.0}


def test_from_formula_multimer():
    """All fields describe the molecular (not empirical) formula, e.g. Hg2[+2] is not reduced to Hg[+1]."""
    s = Solute.from_formula("Hg2+2")
    assert s.formula == "Hg2[+2]"
    assert s.charge == 2
    assert s.molecular_weight == "401.18 g/mol"
    assert s.n_atoms == 2
    assert s.formula_pretty == "Hg2^+2"
    assert s.formula_html == "Hg<sub>2</sub><sup>+2</sup>"
    assert s.formula_latex == "Hg$_{2}$$^{+2}$"
    assert s.formula_hill == "Hg2"

    s = Solute.from_formula("S2O8-2")
    assert s.formula == "S2O8[-2]"
    assert s.molecular_weight == "192.1252 g/mol"
    assert s.formula_pretty == "S2O8^-2"

    # legacy keys with fractional charges are restored to the molecular formula
    s = Solute.from_formula("Br[-0.33333333]")
    assert s.formula == "Br3[-1]"
    assert s.charge == -1
    assert s.molecular_weight == "239.712 g/mol"

    # O6 is not reduced to O3
    s = Solute.from_formula("O6")
    assert s.formula == "O6(aq)"
    assert s.molecular_weight == "95.9964 g/mol"
