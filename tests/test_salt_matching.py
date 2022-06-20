"""
pyEQL salt matching test suite
==============================

This file contains tests for the salt-matching algorithm used by pyEQL in
salt_ion_match.py
"""

import unittest

import numpy as np
import pytest

import pyEQL


class Test_empty_solution(unittest.TestCase):
    """
    test matching a solution that contains no solutes other than water

    """

    # simple case - single-salt
    def setUp(self):
        self.s1 = pyEQL.Solution()

    # The return type should be a salt object
    def test_salt_type(self):
        assert isinstance(self.s1.get_salt(), pyEQL.salt_ion_match.Salt)

    # The salt should be 'HOH'
    def test_salt_formula(self):
        assert self.s1.get_salt().formula == "HOH"

    # The cation should be 'Na+'
    def test_salt_cation(self):
        assert self.s1.get_salt().cation == "H+"

    # The anion should be 'Cl-'
    def test_salt_anion(self):
        assert self.s1.get_salt().anion == "OH-"

    # The cation coefficient should be 1
    def test_salt_nu_cation(self):
        assert self.s1.get_salt().nu_cation == 1

    # The anion coefficient should be 1
    def test_salt_nu_anion(self):
        assert self.s1.get_salt().nu_anion == 1


class Test_single_salt_mono(unittest.TestCase):
    """
    test matching a solution with a single monovalent salt
    ------------------------------------------------

    """

    # simple case - single-salt
    def setUp(self):
        self.s1 = pyEQL.Solution([["Na+", "2 mol/L"], ["Cl-", "2 mol/L"]])

    # The return type should be a salt object
    def test_salt_type(self):
        assert isinstance(self.s1.get_salt(), pyEQL.salt_ion_match.Salt)

    # The salt should be 'NaCl'
    def test_salt_formula(self):
        assert self.s1.get_salt().formula == "NaCl"

    # The cation should be 'Na+'
    def test_salt_cation(self):
        assert self.s1.get_salt().cation == "Na+"

    # The anion should be 'Cl-'
    def test_salt_anion(self):
        assert self.s1.get_salt().anion == "Cl-"

    # The cation coefficient should be 1
    def test_salt_nu_cation(self):
        assert self.s1.get_salt().nu_cation == 1

    # The anion coefficient should be 1
    def test_salt_nu_anion(self):
        assert self.s1.get_salt().nu_anion == 1


class Test_single_salt_di(unittest.TestCase):
    """
    test matching a solution with a single divalent salt
    ------------------------------------------------

    """

    # simple case - single-salt
    def setUp(self):
        self.s1 = pyEQL.Solution([["Na+", "4 mol/L"], ["SO4-2", "2 mol/L"]])

    # The return type should be a salt object
    def test_salt_type(self):
        assert isinstance(self.s1.get_salt(), pyEQL.salt_ion_match.Salt)

    # The salt should be 'Na2SO4'
    def test_salt_formula(self):
        assert self.s1.get_salt().formula == "Na2SO4"

    # The cation should be 'Na+'
    def test_salt_cation(self):
        assert self.s1.get_salt().cation == "Na+"

    # The anion should be 'SO4-2'
    def test_salt_anion(self):
        assert self.s1.get_salt().anion == "SO4-2"

    # The cation coefficient should be 2
    def test_salt_nu_cation(self):
        assert self.s1.get_salt().nu_cation == 2

    # The anion coefficient should be 1
    def test_salt_nu_anion(self):
        assert self.s1.get_salt().nu_anion == 1


class Test_single_salt_di2(unittest.TestCase):
    """
    test matching a solution with a single divalent salt


    """

    # simple case - single-salt
    def setUp(self):
        self.s1 = pyEQL.Solution([["Fe+3", "1 mol/L"], ["Cl-", "3 mol/L"]])

    # The return type should be a salt object
    def test_salt_type(self):
        assert isinstance(self.s1.get_salt(), pyEQL.salt_ion_match.Salt)

    # The salt should be 'FeCl3'
    def test_salt_formula(self):
        assert self.s1.get_salt().formula == "FeCl3"

    # The cation should be 'Fe+3+'
    def test_salt_cation(self):
        assert self.s1.get_salt().cation == "Fe+3"

    # The anion should be 'Cl-'
    def test_salt_anion(self):
        assert self.s1.get_salt().anion == "Cl-"

    # The cation coefficient should be 1
    def test_salt_nu_cation(self):
        assert self.s1.get_salt().nu_cation == 1

    # The anion coefficient should be 3
    def test_salt_nu_anion(self):
        assert self.s1.get_salt().nu_anion == 3


class Test_single_ion(unittest.TestCase):
    """
    test matching a solution containing only a single ion


    """

    # simple case - single-salt
    def setUp(self):
        self.s1 = pyEQL.Solution([["Fe+3", "1 mol/L"]])

    # The return type should be a salt object
    def test_salt_type(self):
        assert isinstance(self.s1.get_salt(), pyEQL.salt_ion_match.Salt)

    # The salt should be 'Fe(OH)3'
    def test_salt_formula(self):
        assert self.s1.get_salt().formula == "Fe(OH)3"

    # The cation should be 'Fe+3'
    def test_salt_cation(self):
        assert self.s1.get_salt().cation == "Fe+3"

    # The anion should be 'OH-'
    def test_salt_anion(self):
        assert self.s1.get_salt().anion == "OH-"

    # The cation coefficient should be 1
    def test_salt_nu_cation(self):
        assert self.s1.get_salt().nu_cation == 1

    # The anion coefficient should be 3
    def test_salt_nu_anion(self):
        assert self.s1.get_salt().nu_anion == 3


class Test_salt_asymmetric(unittest.TestCase):
    """
    test matching a solution where the cation and anion concentrations
    are not equal

    """

    # simple case - single-salt
    def setUp(self):
        self.s1 = pyEQL.Solution([["Na+", "1 mol/kg"], ["Cl-", "4 mol/kg"]])

    # The return type should be a salt object
    def test_salt_type(self):
        assert isinstance(self.s1.get_salt(), pyEQL.salt_ion_match.Salt)

    # The salt should be 'NaCl'
    def test_salt_formula(self):
        assert self.s1.get_salt().formula == "NaCl"

    # The cation should be 'Na+'
    def test_salt_cation(self):
        assert self.s1.get_salt().cation == "Na+"

    # The anion should be 'Cl-'
    def test_salt_anion(self):
        assert self.s1.get_salt().anion == "Cl-"

    # The cation coefficient should be 1
    def test_salt_nu_cation(self):
        assert self.s1.get_salt().nu_cation == 1

    # The anion coefficient should be 1
    def test_salt_nu_anion(self):
        assert self.s1.get_salt().nu_anion == 1
