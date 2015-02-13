'''
pyEQL salt matching test suite
==============================

This file contains tests for the salt-matching algorithm used by pyEQL in
salt_ion_match.py
'''

import pyEQL
import unittest


class Test_empty_solution(unittest.TestCase):
    '''
    test matching a solution that contains no solutes other than water
    
    '''
    # simple case - single-salt
    def setUp(self):
        self.s1 = pyEQL.Solution()
    
    # The return type shoudl be a salt object
    def test_salt_type(self):
        result = self.s1.get_salt()
        expected = pyEQL.salt_ion_match.Salt
        
        self.assertIsInstance(result,expected)
        
    # The salt should be 'HOH'
    def test_salt_formula(self):
        result = self.s1.get_salt().formula
        expected = 'HOH'
        
        self.assertEqual(result,expected)
    
    # The cation should be 'Na+'
    def test_salt_cation(self):
        result = self.s1.get_salt().cation
        expected = 'H+'
        
        self.assertEqual(result,expected)
    
    # The anion should be 'Cl-'
    def test_salt_anion(self):
        result = self.s1.get_salt().anion
        expected = 'OH-'
        
        self.assertEqual(result,expected)
        
    # The cation coefficient should be 1
    def test_salt_nu_cation(self):
        result = self.s1.get_salt().nu_cation
        expected = 1
        
        self.assertEqual(result,expected)
        
    # The anion coefficient should be 1
    def test_salt_nu_anion(self):
        result = self.s1.get_salt().nu_anion
        expected = 1
        
        self.assertEqual(result,expected)

class Test_single_salt_mono(unittest.TestCase):
    '''
    test matching a solution with a single monovalent salt
    ------------------------------------------------
    
    '''
    # simple case - single-salt
    def setUp(self):
        self.s1 = pyEQL.Solution([['Na+','2 mol/L'],['Cl-','2 mol/L']])
    
    # The return type shoudl be a salt object
    def test_salt_type(self):
        result = self.s1.get_salt()
        expected = pyEQL.salt_ion_match.Salt
        
        self.assertIsInstance(result,expected)
        
    # The salt should be 'NaCl'
    def test_salt_formula(self):
        result = self.s1.get_salt().formula
        expected = 'NaCl'
        
        self.assertEqual(result,expected)
    
    # The cation should be 'Na+'
    def test_salt_cation(self):
        result = self.s1.get_salt().cation
        expected = 'Na+'
        
        self.assertEqual(result,expected)
    
    # The anion should be 'Cl-'
    def test_salt_anion(self):
        result = self.s1.get_salt().anion
        expected = 'Cl-'
        
        self.assertEqual(result,expected)
        
    # The cation coefficient should be 1
    def test_salt_nu_cation(self):
        result = self.s1.get_salt().nu_cation
        expected = 1
        
        self.assertEqual(result,expected)
        
    # The anion coefficient should be 1
    def test_salt_nu_anion(self):
        result = self.s1.get_salt().nu_anion
        expected = 1
        
        self.assertEqual(result,expected)
        
class Test_single_salt_di(unittest.TestCase):
    '''
    test matching a solution with a single divalent salt
    ------------------------------------------------
    
    '''
    # simple case - single-salt
    def setUp(self):
        self.s1 = pyEQL.Solution([['Na+','4 mol/L'],['SO4-2','2 mol/L']])
    
    # The return type shoudl be a salt object
    def test_salt_type(self):
        result = self.s1.get_salt()
        expected = pyEQL.salt_ion_match.Salt
        
        self.assertIsInstance(result,expected)
        
    # The salt should be 'Na2SO4'
    def test_salt_formula(self):
        result = self.s1.get_salt().formula
        expected = 'Na2SO4'
        
        self.assertEqual(result,expected)
    
    # The cation should be 'Na+'
    def test_salt_cation(self):
        result = self.s1.get_salt().cation
        expected = 'Na+'
        
        self.assertEqual(result,expected)
    
    # The anion should be 'SO4-2'
    def test_salt_anion(self):
        result = self.s1.get_salt().anion
        expected = 'SO4-2'
        
        self.assertEqual(result,expected)
        
    # The cation coefficient should be 2
    def test_salt_nu_cation(self):
        result = self.s1.get_salt().nu_cation
        expected = 2
        
        self.assertEqual(result,expected)
        
    # The anion coefficient should be 1
    def test_salt_nu_anion(self):
        result = self.s1.get_salt().nu_anion
        expected = 1
        
        self.assertEqual(result,expected)

class Test_single_salt_di2(unittest.TestCase):
    '''
    test matching a solution with a single divalent salt
    
    
    '''
    # simple case - single-salt
    def setUp(self):
        self.s1 = pyEQL.Solution([['Fe+3','1 mol/L'],['Cl-','3 mol/L']])
    
    # The return type shoudl be a salt object
    def test_salt_type(self):
        result = self.s1.get_salt()
        expected = pyEQL.salt_ion_match.Salt
        
        self.assertIsInstance(result,expected)
        
    # The salt should be 'FeCl3'
    def test_salt_formula(self):
        result = self.s1.get_salt().formula
        expected = 'FeCl3'
        
        self.assertEqual(result,expected)
    
    # The cation should be 'Fe+3+'
    def test_salt_cation(self):
        result = self.s1.get_salt().cation
        expected = 'Fe+3'
        
        self.assertEqual(result,expected)
    
    # The anion should be 'Cl-'
    def test_salt_anion(self):
        result = self.s1.get_salt().anion
        expected = 'Cl-'
        
        self.assertEqual(result,expected)
        
    # The cation coefficient should be 1
    def test_salt_nu_cation(self):
        result = self.s1.get_salt().nu_cation
        expected = 1
        
        self.assertEqual(result,expected)
        
    # The anion coefficient should be 3
    def test_salt_nu_anion(self):
        result = self.s1.get_salt().nu_anion
        expected = 3
        
        self.assertEqual(result,expected)

class Test_single_ion(unittest.TestCase):
    '''
    test matching a solution containing only a single ion
    
    
    '''
    # simple case - single-salt
    def setUp(self):
        self.s1 = pyEQL.Solution([['Fe+3','1 mol/L']])
    
    # The return type shoudl be a salt object
    def test_salt_type(self):
        result = self.s1.get_salt()
        expected = pyEQL.salt_ion_match.Salt
        
        self.assertIsInstance(result,expected)
        
    # The salt should be 'Fe(OH)3'
    def test_salt_formula(self):
        result = self.s1.get_salt().formula
        expected = 'Fe(OH)3'
        
        self.assertEqual(result,expected)
    
    # The cation should be 'Fe+3'
    def test_salt_cation(self):
        result = self.s1.get_salt().cation
        expected = 'Fe+3'
        
        self.assertEqual(result,expected)
    
    # The anion should be 'OH-'
    def test_salt_anion(self):
        result = self.s1.get_salt().anion
        expected = 'OH-'
        
        self.assertEqual(result,expected)
        
    # The cation coefficient should be 1
    def test_salt_nu_cation(self):
        result = self.s1.get_salt().nu_cation
        expected = 1
        
        self.assertEqual(result,expected)
        
    # The anion coefficient should be 3
    def test_salt_nu_anion(self):
        result = self.s1.get_salt().nu_anion
        expected = 3
        
        self.assertEqual(result,expected)

class Test_salt_asymmetric(unittest.TestCase):
    '''
    test matching a solution where the cation and anion concentrations
    are not equal    
    
    '''
    # simple case - single-salt
    def setUp(self):
        self.s1 = pyEQL.Solution([['Na+','1 mol/kg'],['Cl-','4 mol/kg']])
    
    # The return type shoudl be a salt object
    def test_salt_type(self):
        result = self.s1.get_salt()
        expected = pyEQL.salt_ion_match.Salt
        
        self.assertIsInstance(result,expected)
        
    # The salt should be 'NaCl'
    def test_salt_formula(self):
        result = self.s1.get_salt().formula
        expected = 'NaCl'
        
        self.assertEqual(result,expected)
    
    # The cation should be 'Na+'
    def test_salt_cation(self):
        result = self.s1.get_salt().cation
        expected = 'Na+'
        
        self.assertEqual(result,expected)
    
    # The anion should be 'Cl-'
    def test_salt_anion(self):
        result = self.s1.get_salt().anion
        expected = 'Cl-'
        
        self.assertEqual(result,expected)
        
    # The cation coefficient should be 1
    def test_salt_nu_cation(self):
        result = self.s1.get_salt().nu_cation
        expected = 1
        
        self.assertEqual(result,expected)
        
    # The anion coefficient should be 1
    def test_salt_nu_anion(self):
        result = self.s1.get_salt().nu_anion
        expected = 1
        
        self.assertEqual(result,expected)

if __name__ == '__main__':
    unittest.main()
