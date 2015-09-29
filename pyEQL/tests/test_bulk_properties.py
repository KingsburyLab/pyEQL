'''
pyEQL test suite for bulk property calculations
===============================================

This file contains tests for the bulk properties calculated by
Solution class methods. Currently included methods are:

- get_hardness()

'''

import pyEQL
import unittest

class test_hardness(unittest.TestCase):
    '''
    test the get_hardness() method
    ------------------------------   
    '''
    # an empty solution should have zero hardness
    def test_empty_solution(self):
        s1 = pyEQL.Solution()
        result = s1.get_hardness().magnitude
        expected = 0
        
        self.assertEqual(result,expected)
    
    # a solution with only monovalent ions should have zero hardness
    def test_hardness_1(self):
        s1 = pyEQL.Solution([['Na+','0.2 mol/L'],['Cl-','0.2 mol/L']])
        result = s1.get_hardness().magnitude
        expected = 0
        
        self.assertEqual(result,expected)
        
    # a solution with only multivalent anions should have zero hardness
    def test_hardness_2(self):
        s1 = pyEQL.Solution([['Na+','0.4 mol/L'],['SO4-2','0.2 mol/L']])
        result = s1.get_hardness().magnitude
        expected = 0
        
        self.assertEqual(result,expected)
    
    # the hardness should return the equivalent concentration, not just the
    # molar concentration (e.g. multiply by the charge)
    def test_hardness_3(self):
        s1 = pyEQL.Solution([['Fe+3','0.1 mol/L'],['Cl-','0.3 mol/L']])
        result = s1.get_hardness().magnitude
        expected = 15013.5
        
        self.assertAlmostEqual(result,expected)
    
    # the hardness should account for multiple cations but count only those
    # that are multivalent
    def test_hardness_4(self):
        s1 = pyEQL.Solution([['Na+','0.1 mol/L'],['K+','0.1 mol/L'],['Mg+2','0.1 mol/L'],['Ca+2','0.1 mol/L'],['Fe+3','0.1 mol/L'],['Cl-','0.1 mol/L'],['F-','0.1 mol/L'],['SO4-2','0.2 mol/L'],['PO4-3','0.1 mol/L']])
        result = s1.get_hardness().magnitude
        expected = 35031.5
        
        self.assertAlmostEqual(result,expected)
    
    # the hardness should return g/L units
    def test_hardness_5(self):
        s1 = pyEQL.Solution([['Fe+3','0.1 mol/L'],['Cl-','0.3 mol/L']])
        result = str(s1.get_hardness().dimensionality)
        expected = '[mass] / [length] ** 3'
        
        self.assertEqual(result,expected)
        
if __name__ == '__main__':
    unittest.main()
