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
import pytest

# relative tolerance between experimental and computed properties for this test file
RTOL = 0.25


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
        expected = [0.613, 0.5039, 0.4448, 0.3928]

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
        expected = [0.582, 0.4523, 0.3827, 0.3138]

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
        expected = [0.6211, 0.5181, 0.4481, 0.4132]

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
        expected = [0.582, 0.4529, 0.3635, 0.3133]

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
        expected = [0.6293, 0.5327, 0.4680, 0.4364]

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
        expected = [0.5942, 0.4731, 0.3878, 0.3370]

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
