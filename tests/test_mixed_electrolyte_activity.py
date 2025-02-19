"""
pyEQL test suite for mixed electrolyte activity
===============================================

This file contains tests for the Effective Pitzer Model
implemented in pyEQL.

This test suite creates several mixed salt solutions and compares
pyEQL's calculations with experimental data for activity coefficients
and/or osmotic coefficients.

The Effective Pitzer Model is described in Mistry et al.
DOI: 10.1016/j.desal.2013.03.015

Experimental data for mixed electrolytes is obtained from Rodil et al.
DOI: 10.1021/je9004432

"""

import numpy as np

import pyEQL

# relative tolerance between experimental and computed properties for this test file
RTOL = 0.45


class Test_nano3_kno3_activity:
    """
    test mean activity coefficients in a NaNO3 + KNO3 mixture
    ---------------------------------------------------------
        Note that the values given in Table 3 of Rodil et al. are single-ion
        activity coefficients that must be averaged (geometrically) to obtained
        the mean values below

    """

    def test_activity_Na_XNa_75(self):
        # test the activity coefficient of Na+ in mixed NaNO3 and KNO3 when the
        # mole fraction of Na+ is 0.75
        # molality refers to the total nitrate ion molality
        molality = [0.3997, 0.992, 1.584, 2.373]
        expected = [0.6288, 0.5323, 0.4812, 0.4361]

        for item in range(len(molality)):
            s1 = pyEQL.Solution(
                [
                    ["Na+", str(0.75 * molality[item]) + "mol/kg"],
                    ["K+", str(0.25 * molality[item]) + "mol/kg"],
                    ["NO3-", str(molality[item]) + "mol/kg"],
                ]
            )
            result = s1.get_activity_coefficient("Na+").magnitude

            assert np.isclose(result, expected[item], RTOL)

    def test_activity_K_XNa_75(self):
        # test the activity coefficient of Na+ in mixed NaNO3 and KNO3 when the
        # mole fraction of Na+ is 0.75
        # molality refers to the total nitrate ion molality
        molality = [0.3997, 0.992, 1.584, 2.373]
        expected = [0.5939, 0.4727, 0.4059, 0.3370]

        for item in range(len(molality)):
            s1 = pyEQL.Solution(
                [
                    ["Na+", str(0.75 * molality[item]) + "mol/kg"],
                    ["K+", str(0.25 * molality[item]) + "mol/kg"],
                    ["NO3-", str(molality[item]) + "mol/kg"],
                ]
            )
            result = s1.get_activity_coefficient("K+").magnitude

            assert np.isclose(result, expected[item], RTOL)

    def test_activity_Na_XNa_50(self):
        # test the activity coefficient of Na+ in mixed NaNO3 and KNO3 when the
        # mole fraction of Na+ is 0.50
        # molality refers to the total nitrate ion molality
        molality = [0.4005, 0.9926, 1.787, 2.384]
        expected = [0.588, 0.463, 0.3754, 0.3245]

        for item in range(len(molality)):
            s1 = pyEQL.Solution(
                [
                    ["Na+", str(0.5 * molality[item]) + "mol/kg"],
                    ["K+", str(0.5 * molality[item]) + "mol/kg"],
                    ["NO3-", str(molality[item]) + "mol/kg"],
                ]
            )
            result = s1.get_activity_coefficient("Na+").magnitude

            assert np.isclose(result, expected[item], RTOL)

    def test_activity_K_XNa_50(self):
        # test the activity coefficient of Na+ in mixed NaNO3 and KNO3 when the
        # mole fraction of Na+ is 0.50
        # molality refers to the total nitrate ion molality
        molality = [0.4005, 0.9926, 1.787, 2.384]
        expected = [0.5879, 0.4630, 0.3754, 0.3245]

        for item in range(len(molality)):
            s1 = pyEQL.Solution(
                [
                    ["Na+", str(0.5 * molality[item]) + "mol/kg"],
                    ["K+", str(0.5 * molality[item]) + "mol/kg"],
                    ["NO3-", str(molality[item]) + "mol/kg"],
                ]
            )
            result = s1.get_activity_coefficient("K+").magnitude

            assert np.isclose(result, expected[item], RTOL)

    def test_activity_Na_XNa_25(self):
        # test the activity coefficient of Na+ in mixed NaNO3 and KNO3 when the
        # mole fraction of Na+ is 0.25
        # molality refers to the total nitrate ion molality
        molality = [0.4021, 0.9976, 1.794, 2.393]
        expected = [0.6134, 0.5042, 0.4295, 0.3931]
        for item in range(len(molality)):
            s1 = pyEQL.Solution(
                [
                    ["Na+", str(0.25 * molality[item]) + "mol/kg"],
                    ["K+", str(0.75 * molality[item]) + "mol/kg"],
                    ["NO3-", str(molality[item]) + "mol/kg"],
                ]
            )
            result = s1.get_activity_coefficient("Na+").magnitude

            assert np.isclose(result, expected[item], RTOL)

    def test_activity_K_XNa_25(self):
        # test the activity coefficient of Na+ in mixed NaNO3 and KNO3 when the
        # mole fraction of Na+ is 0.25
        # molality refers to the total nitrate ion molality
        molality = [0.4021, 0.9976, 1.794, 2.393]
        expected = [0.582, 0.4533, 0.3642, 0.3137]

        for item in range(len(molality)):
            s1 = pyEQL.Solution(
                [
                    ["Na+", str(0.25 * molality[item]) + "mol/kg"],
                    ["K+", str(0.75 * molality[item]) + "mol/kg"],
                    ["NO3-", str(molality[item]) + "mol/kg"],
                ]
            )
            result = s1.get_activity_coefficient("K+").magnitude

            assert np.isclose(result, expected[item], RTOL)
