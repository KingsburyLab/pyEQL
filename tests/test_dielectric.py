"""
pyEQL dielectric constant test suite
============================================

This file contains tests that check the dielectric constant
computations of pyEQL

"""

import numpy as np
import pyEQL
import pytest

# relative tolerance between experimental and computed properties for this test file
RTOL = 0.04


class Test_dielectric:
    """
    test the Dielectric Constant calculations of various solutions
    ------------------------------------------------

    Reference: A. Zuber, L. Cardozo-Filho, V.F. Cabral, R.F. Checoni, M. Castier, An empirical equation for the dielectric constant in aqueous and
    nonaqueous electrolyte mixtures, Fluid Phase Equilib. 376 (2014) 116-123. doi:10.1016/j.fluid.2014.05.037.

    """

    def test_dielectric_constant(self):
        """
        4.4 mol/kg NaCl = 46
        """
        s1 = pyEQL.Solution([["Na+", "4.4 mol/kg"], ["Cl-", "4.4 mol/kg"]])

        result = s1.dielectric_constant.magnitude
        expected = 46

        assert np.isclose(result, expected, rtol=RTOL)

    def test_dielectric_constant2(self):
        """
        2 mol/kg NaCl = 58
        """
        s1 = pyEQL.Solution([["Na+", "2 mol/kg"], ["Cl-", "2 mol/kg"]])

        result = s1.dielectric_constant.magnitude
        expected = 58

        assert np.isclose(result, expected, rtol=RTOL)

    def test_dielectric_constant3(self):
        """
        1 mol/kg NaCl = 66
        """
        s1 = pyEQL.Solution([["Na+", "1 mol/kg"], ["Cl-", "1 mol/kg"]])

        result = s1.dielectric_constant.magnitude
        expected = 66

        assert np.isclose(result, expected, rtol=RTOL)

    def test_dielectric_constant4(self):
        """
        1 mol/kg KBr = 67
        """
        s1 = pyEQL.Solution([["K+", "1 mol/kg"], ["Br-", "1 mol/kg"]])

        result = s1.dielectric_constant.magnitude
        expected = 67

        assert np.isclose(result, expected, rtol=RTOL)

    def test_dielectric_constant5(self):
        """
        3.4 mol/kg KBr = 51
        """
        s1 = pyEQL.Solution([["K+", "3.4 mol/kg"], ["Br-", "3.4 mol/kg"]])

        result = s1.dielectric_constant.magnitude
        expected = 51

        assert np.isclose(result, expected, rtol=RTOL)

    def test_dielectric_constant6(self):
        """
        5 mol/kg LiCl = 39
        """
        s1 = pyEQL.Solution([["Li+", "5 mol/kg"], ["Cl-", "5 mol/kg"]])

        result = s1.dielectric_constant.magnitude
        expected = 39

        assert np.isclose(result, expected, rtol=RTOL)

    def test_dielectric_constant7(self):
        """
        1 mol/kg LiCl = 64
        """
        s1 = pyEQL.Solution([["Li+", "1 mol/kg"], ["Cl-", "1 mol/kg"]])

        result = s1.dielectric_constant.magnitude
        expected = 64

        assert np.isclose(result, expected, rtol=RTOL)

    @pytest.mark.xfail()
    def test_dielectric_constant8(self):
        """
        12 mol/kg LiCl = 24
        """
        s1 = pyEQL.Solution([["Li+", "12 mol/kg"], ["Cl-", "12 mol/kg"]])

        result = s1.dielectric_constant.magnitude
        expected = 24

        assert np.isclose(result, expected, rtol=RTOL)

    def test_dielectric_constant9(self):
        """
        6.5 mol/kg RbCl = 43
        """
        s1 = pyEQL.Solution([["Rb+", "6.5 mol/kg"], ["Cl-", "6.5 mol/kg"]])

        result = s1.dielectric_constant.magnitude
        expected = 43

        assert np.isclose(result, expected, rtol=RTOL)

    def test_dielectric_constant10(self):
        """
        2.1 mol/kg RbCl = 59
        """
        s1 = pyEQL.Solution([["Rb+", "2.1 mol/kg"], ["Cl-", "2.1 mol/kg"]])

        result = s1.dielectric_constant.magnitude
        expected = 59

        assert np.isclose(result, expected, rtol=RTOL)

    def test_dielectric_constant11(self):
        """
        0.5 mol/kg RbCl = 73
        """
        s1 = pyEQL.Solution([["Rb+", "0.5 mol/kg"], ["Cl-", "0.5 mol/kg"]])

        result = s1.dielectric_constant.magnitude
        expected = 73

        assert np.isclose(result, expected, rtol=RTOL)
