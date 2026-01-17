"""
pyEQL volume and concentration methods test suite
=================================================

This file contains tests for the volume and concentration-related methods
used by pyEQL's Solution class
"""

import numpy as np
import pytest

from pyEQL import Solution


@pytest.fixture
def s1():
    return Solution(volume="2 L", engine="pyeql")


@pytest.fixture
def s2():
    return Solution([["Na+", "4 mol/L"], ["Cl-", "4 mol/L"]], volume="2 L", engine="pyeql")


def test_diffusion_transport(s1, s2):
    # test ionic strength adjustment
    # assert s1.get_diffusion_coefficient("H+", use_engine=True) > s2.get_diffusion_coefficient("H+", use_engine=True)
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
        factor * s2.get_diffusion_coefficient("Na+", use_engine=True).magnitude,
        s2.get_diffusion_coefficient("Na+", use_engine=True).magnitude,
        atol=5e-11,
    )
    s_dilute = Solution({"Na+": "1 mmol/L", "Cl-": "1 mmol/L"}, engine="pyeql")
    assert np.isclose(
        s_dilute.get_diffusion_coefficient("Na+", activity_correction=False, use_engine=True).magnitude,
        1.334e-9,
        atol=1e-11,
    )
    assert np.isclose(s_dilute.get_transport_number("Na+", use_engine=True), 0.396, atol=1e-3)
    assert np.isclose(s_dilute.get_transport_number("Cl-", use_engine=True), 0.604, atol=1e-3)

    # test setting a default value
    s2.default_diffusion_coeff = 0
    assert s2.get_diffusion_coefficient("Cs+", use_engine=True).magnitude == 0
    s2.default_diffusion_coeff = 1e-9
    assert s2.get_diffusion_coefficient("Cs+", activity_correction=False, use_engine=True).magnitude == 1e-9
    s2.default_diffusion_coeff = 0
    assert s2.get_diffusion_coefficient("Cs+", activity_correction=True, use_engine=True).magnitude < 1e-9
    d25 = s2.get_diffusion_coefficient("Na+", activity_correction=False, use_engine=True).magnitude
    nu25 = s2.water_substance.nu
    s2.temperature = "40 degC"
    d40 = s2.get_diffusion_coefficient("Na+", activity_correction=False, use_engine=True).magnitude
    nu40 = s2.water_substance.nu
    assert np.isclose(
        d40,
        d25 * np.exp(122 / (273.15 + 40) - 122 / 298.15) * (nu25 / nu40),
        atol=5e-10,
    )

    # test correction factors for concentration, as per Appelo 2017 Fig 5
    D1 = (
        Solution({"Na+": "1 umol/L", "Cl-": "1 umol/L"}, engine="pyeql")
        .get_diffusion_coefficient("Na+", use_engine=True)
        .magnitude
    )
    D2 = (
        Solution({"Na+": "1.7 mol/kg", "Cl-": "1.7 mol/kg"}, engine="pyeql")
        .get_diffusion_coefficient("Na+", use_engine=True)
        .magnitude
    )
    assert np.isclose(D2 / D1, 1, atol=1e-2)

    D1 = (
        Solution({"K+": "1 umol/L", "Cl-": "1 umol/L"}, engine="pyeql")
        .get_diffusion_coefficient("K+", use_engine=True)
        .magnitude
    )
    D2 = (
        Solution({"K+": "0.5 mol/kg", "Cl-": "0.5 mol/kg"}, engine="pyeql")
        .get_diffusion_coefficient("K+", use_engine=True)
        .magnitude
    )
    assert np.isclose(D2 / D1, 1, atol=1e-2)
