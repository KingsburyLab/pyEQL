"""
pyEQL salt matching test suite
==============================

This file contains tests for the salt-matching algorithm used by pyEQL in
salt_ion_match.py
"""
import platform

import numpy as np
import pytest

import pyEQL
from pyEQL.salt_ion_match import Salt


def test_salt_init():
    s = Salt("Na[+1]", "Cl[-1]")
    assert s.formula == "NaCl"
    assert s.cation == "Na[+1]"
    assert s.anion == "Cl[-1]"
    assert s.nu_anion == 1
    assert s.nu_cation == 1
    assert s.z_cation == 1
    assert s.z_anion == -1

    s = Salt("Fe+3", "OH-1")
    assert s.formula == "Fe(OH)3"
    assert s.cation == "Fe[+3]"
    assert s.anion == "OH[-1]"
    assert s.nu_anion == 3
    assert s.nu_cation == 1
    assert s.z_cation == 3
    assert s.z_anion == -1


def test_empty_solution():
    """
    test matching a solution that contains no solutes other than water
    """
    s1 = pyEQL.Solution()
    assert isinstance(s1.get_salt(), pyEQL.salt_ion_match.Salt)
    assert s1.get_salt().formula == "HOH"
    assert s1.get_salt().cation == "H[+1]"
    assert s1.get_salt().anion == "OH[-1]"
    assert s1.get_salt().nu_cation == 1
    assert s1.get_salt().nu_anion == 1


def test_single_salt_mono():
    """
    test matching a solution with a single monovalent salt
    """
    s1 = pyEQL.Solution([["Na+", "2 mol/L"], ["Cl-", "2 mol/L"]])
    assert isinstance(s1.get_salt(), pyEQL.salt_ion_match.Salt)
    assert s1.get_salt().formula == "NaCl"
    assert s1.get_salt().cation == "Na[+1]"
    assert s1.get_salt().anion == "Cl[-1]"
    assert s1.get_salt().nu_cation == 1
    assert s1.get_salt().nu_anion == 1


def test_single_salt_di():
    """
    test matching a solution with a single divalent salt
    """
    s1 = pyEQL.Solution([["Na+", "4 mol/L"], ["SO4-2", "2 mol/L"]])
    assert isinstance(s1.get_salt(), pyEQL.salt_ion_match.Salt)
    assert s1.get_salt().formula == "Na2SO4"
    assert s1.get_salt().cation == "Na[+1]"
    assert s1.get_salt().anion == "SO4[-2]"
    assert s1.get_salt().nu_cation == 2
    assert s1.get_salt().nu_anion == 1


def test_single_salt_di2():
    """
    test matching a solution with a single divalent salt
    """
    s1 = pyEQL.Solution([["Fe+3", "1 mol/L"], ["Cl-", "3 mol/L"]])
    assert isinstance(s1.get_salt(), pyEQL.salt_ion_match.Salt)
    assert s1.get_salt().formula == "FeCl3"
    assert s1.get_salt().cation == "Fe[+3]"
    assert s1.get_salt().anion == "Cl[-1]"
    assert s1.get_salt().nu_cation == 1
    assert s1.get_salt().nu_anion == 3


def test_single_ion():
    """
    test matching a solution containing only a single ion
    """
    s1 = pyEQL.Solution([["Fe+3", "1 mol/L"]])
    assert isinstance(s1.get_salt(), pyEQL.salt_ion_match.Salt)
    assert s1.get_salt().formula == "Fe(OH)3"
    assert s1.get_salt().cation == "Fe[+3]"
    assert s1.get_salt().anion == "OH[-1]"
    assert s1.get_salt().nu_cation == 1
    assert s1.get_salt().nu_anion == 3


@pytest.mark.skipif(platform.machine() == "arm64" and platform.system() == "Darwin", reason="arm64 not supported")
def test_salt_with_equilibration():
    """
    test matching a solution containing a salt, before and after equilibration.
    Due to speciation changes, the concentration of the salt will decrease unless
    get_salt_dict() uses total concentrations
    """
    s1 = pyEQL.Solution({"Mg+2": "1 mol/L", "Cl-": "2 mol/L"})
    assert isinstance(s1.get_salt(), pyEQL.salt_ion_match.Salt)
    assert s1.get_salt().formula == "MgCl2"
    assert s1.get_salt().cation == "Mg[+2]"
    assert s1.get_salt().anion == "Cl[-1]"
    assert s1.get_salt().nu_cation == 1
    assert s1.get_salt().nu_anion == 2
    assert np.isclose(s1.get_salt_dict()["MgCl2"]["mol"], 1)
    s1.equilibrate()
    assert s1.get_salt().formula == "MgCl2"
    assert np.isclose(s1.get_salt_dict()["MgCl2"]["mol"], 1)


def test_salt_asymmetric():
    """
    test matching a solution where the cation and anion concentrations
    are not equal
    """
    s1 = pyEQL.Solution([["Na+", "1 mol/kg"], ["Cl-", "4 mol/kg"]])
    assert isinstance(s1.get_salt(), pyEQL.salt_ion_match.Salt)
    assert s1.get_salt().formula == "NaCl"
    assert s1.get_salt().cation == "Na[+1]"
    assert s1.get_salt().anion == "Cl[-1]"
    assert s1.get_salt().nu_cation == 1
    assert s1.get_salt().nu_anion == 1
