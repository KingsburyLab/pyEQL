"""
Tests of pyEQL.functions module

"""
import platform

import numpy as np
import pytest

from pyEQL import Solution
from pyEQL.functions import entropy_mix, gibbs_mix


@pytest.fixture()
def s1():
    return Solution(volume="2 L")


@pytest.fixture()
def s2():
    return Solution({"Na+": "1 mol/L", "Cl-": "1 mol/L"}, volume="10 L")


@pytest.fixture()
def s1_p():
    return Solution(volume="2 L", engine="phreeqc")


@pytest.fixture()
def s2_p():
    return Solution({"Na+": "1 mol/L", "Cl-": "1 mol/L"}, volume="10 L", engine="phreeqc")


@pytest.fixture()
def s1_i():
    return Solution(volume="2 L", engine="ideal")


@pytest.fixture()
def s2_i():
    return Solution({"Na+": "1 mol/L", "Cl-": "1 mol/L"}, volume="10 L", engine="ideal")


@pytest.mark.skipif(platform.machine() == "arm64" and platform.system() == "Darwin", reason="arm64 not supported")
def test_mixing_functions(s1, s2, s1_p, s2_p, s1_i, s2_i):
    # mixing energy and entropy of any solution with itself should be zero
    assert np.isclose(gibbs_mix(s1, s1).magnitude, 0)
    assert np.isclose(entropy_mix(s2, s2).magnitude, 0)
    assert np.isclose(gibbs_mix(s1_p, s1_p).magnitude, 0)
    assert np.isclose(entropy_mix(s2_p, s2_p).magnitude, 0)
    assert np.isclose(gibbs_mix(s1_i, s1_i).magnitude, 0, atol=1e-6)
    assert np.isclose(entropy_mix(s2_i, s2_i).magnitude, 0)

    # TODO - I have not tested how equilibrate() affects the results
    for dil, conc in zip([s1, s1_p, s1_i], [s2, s2_p, s2_i]):
        # for mixing 1 and 2, we should have
        # H20: 55.5 * 2 mol + 55.5 * 10 mol, x1 = 0.9999 x2 = 0.9645, mixture = 0.9703 = approximately -9043 J
        s_theoretical = (
            8.314
            * 298.15
            * (
                (dil + conc).get_amount("H2O", "mol").magnitude
                * np.log((dil + conc).get_amount("H2O", "fraction").magnitude)
                + (dil + conc).get_amount("Na+", "mol").magnitude
                * np.log((dil + conc).get_amount("Na+", "fraction").magnitude)
                + (dil + conc).get_amount("Cl-", "mol").magnitude
                * np.log((dil + conc).get_amount("Cl-", "fraction").magnitude)
                - dil.get_amount("H2O", "mol").magnitude * np.log(dil.get_amount("H2O", "fraction").magnitude)
                - conc.get_amount("H2O", "mol").magnitude * np.log(conc.get_amount("H2O", "fraction").magnitude)
                - conc.get_amount("Na+", "mol").magnitude * np.log(conc.get_amount("Na+", "fraction").magnitude)
                - conc.get_amount("Cl-", "mol").magnitude * np.log(conc.get_amount("Cl-", "fraction").magnitude)
            )
        )
        assert np.isclose(entropy_mix(dil, conc).magnitude, s_theoretical, rtol=0.005)
        g_theoretical = (
            8.314
            * 298.15
            * (
                (dil + conc).get_amount("H2O", "mol").magnitude * np.log((dil + conc).get_activity("H2O").magnitude)
                + (dil + conc).get_amount("Na+", "mol").magnitude * np.log((dil + conc).get_activity("Na+").magnitude)
                + (dil + conc).get_amount("Cl-", "mol").magnitude * np.log((dil + conc).get_activity("Cl-").magnitude)
                - dil.get_amount("H2O", "mol").magnitude * np.log(dil.get_activity("H2O").magnitude)
                - conc.get_amount("H2O", "mol").magnitude * np.log(conc.get_activity("H2O").magnitude)
                - conc.get_amount("Na+", "mol").magnitude * np.log(conc.get_activity("Na+").magnitude)
                - conc.get_amount("Cl-", "mol").magnitude * np.log(conc.get_activity("Cl-").magnitude)
            )
        )
        assert np.isclose(gibbs_mix(dil, conc).magnitude, g_theoretical, rtol=0.005)
