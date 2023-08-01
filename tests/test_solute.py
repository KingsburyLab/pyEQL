"""
Tests for the solute.py module
"""
import numpy as np
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
    assert s.oxi_state_guesses == ({"Mg": 2.0},)
    assert np.isclose(s.molecular_weight, 24.305)
