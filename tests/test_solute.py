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
    assert Solute.from_formula("Br[-0.33333333]").oxi_state_guesses == {}
    assert s.molecular_weight == "24.305 g/mol"
    s2 = Solute.from_formula("O6")
    assert s2.formula == "O3(aq)"
    assert s2.molecular_weight == "47.9982 g/mol"
    assert s2.oxi_state_guesses == {"O": 0.0}
