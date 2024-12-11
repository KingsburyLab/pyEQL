"""
pyEQL debye length test suite
============================================

This file contains tests that check the Debye Length
computations of pyEQL

"""

import numpy as np

import pyEQL

# relative tolerance between experimental and computed properties for this test file
RTOL = 0.03


class Test_debye_length:
    """
    test the Debye Length calculations of various solutions
    ------------------------------------------------

    Reference: [1] M. Hu, B. Mi, Enabling graphene oxide nanosheets as water separation membranes,
    Environ. Sci. Technol. 47 (2013) 3715-3723. doi:10.1021/es400571g.

    0.1 mM NaCl: 31nm
    10 mM NaCl: 3.1 nm
    0.1 mM Na2SO4: 18nm
    10 mM Na2SO4: 1.8nm

    """

    def test_debye_length_1(self):
        """ """
        s1 = pyEQL.Solution([["Na+", "0.1 mmol/L"], ["Cl-", "0.1 mmol/L"]])

        result = s1.debye_length.magnitude
        expected = 31

        assert np.isclose(result, expected, rtol=RTOL)

    def test_debye_length_2(self):
        """ """
        s1 = pyEQL.Solution([["Na+", "10 mmol/L"], ["Cl-", "10 mmol/L"]])

        result = s1.debye_length.magnitude
        expected = 3.1

        assert np.isclose(result, expected, rtol=RTOL)

    def test_debye_length_3(self):
        """ """
        s1 = pyEQL.Solution([["Na+", "0.2 mmol/L"], ["SO4-2", "0.1 mmol/L"]])

        result = s1.debye_length.magnitude
        expected = 18

        assert np.isclose(result, expected, rtol=RTOL)

    def test_debye_length_4(self):
        """ """
        s1 = pyEQL.Solution([["Na+", "20 mmol/L"], ["SO4-2", "10 mmol/L"]])

        result = s1.debye_length.magnitude
        expected = 1.8

        assert np.isclose(result, expected, rtol=RTOL)
