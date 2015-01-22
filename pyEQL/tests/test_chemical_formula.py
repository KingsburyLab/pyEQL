'''
chemical_formula.py test suite
==============================

This file contains tests for the chemical formula interpreter module of pyEQL.
'''

from pyEQL import chemical_formula as cf
import unittest

class Test_check_formula(unittest.TestCase):
    '''
    tests for the _check_formula() function
    ---------------------------------------
    This function is the core of the library and enforces the rules for 
    chemical formula formatting.
    '''
      
    def test_check_formula_1(self):
        input = 'Fe2(SO4)3'
        result = cf._check_formula(input)        
        expected = ['Fe', '2', '(', 'S', 'O', '4', ')', '3']
        
        self.assertEqual(result,expected)
    
    def test_check_formula_2(self):
        input = 'C7H16'
        result = cf._check_formula(input)        
        expected = ['C', '7', 'H', '16']
        
        self.assertEqual(result,expected)
    
    def test_check_formula_3(self):
        input = '(NH3)2SO4'
        result = cf._check_formula(input)        
        expected = ['(', 'N', 'H', '3', ')', '2', 'S', 'O', '4']
        
        self.assertEqual(result,expected)
    
    def test_check_formula_4(self):
        input = 'MgCl2'
        result = cf._check_formula(input)        
        expected = ['Mg', 'Cl', '2']
        
        self.assertEqual(result,expected)
    
    def test_check_formula_5(self):
        input = 'C100H202'
        result = cf._check_formula(input)        
        expected = ['C', '100', 'H', '202']
        
        self.assertEqual(result,expected)
    
    def test_check_formula_6(self):
        input = 'Fe+++'
        result = cf._check_formula(input)        
        expected = ['Fe', '+++']
        
        self.assertEqual(result,expected)
    
    def test_check_formula_7(self):
        input = 'V+4'
        result = cf._check_formula(input)        
        expected = ['V', '+', '4']
        
        self.assertEqual(result,expected)
#    def test_choice(self):
#        element = random.choice(self.seq)
#        self.assertTrue(element in self.seq)
#
#    def test_sample(self):
#        with self.assertRaises(ValueError):
#            random.sample(self.seq, 20)
#        for element in random.sample(self.seq, 5):
#            self.assertTrue(element in self.seq)

class Test_is_valid_formula(unittest.TestCase):
    '''
    tests for is_valid_formula()
    ----------------------------
    '''
    # A formula must start with a letter or an open parenthesis
    def test_is_valid_formula_1(self):
        input = '(NH3)2'
        result = cf.is_valid_formula(input)        
        expected = True
        
        self.assertEqual(result,expected)
    
    # A formula must start with a letter or an open parenthesis
    def test_is_valid_formula_2(self):
        input = '3CO3-'
        result = cf.is_valid_formula(input)        
        expected = False
        
        self.assertEqual(result,expected)
        
    # A formula cannot contain any non-alphanumeric characters beside '(', ')', '+', and '-'
    def test_is_valid_formula_3(self):
        input = 'Na^+'
        result = cf.is_valid_formula(input)        
        expected = False
        
        self.assertEqual(result,expected)

    # A formula cannot contain both '+' and '-'     
    def test_is_valid_formula_4(self):
        input = 'Na+-+'
        result = cf.is_valid_formula(input)        
        expected = False
        
        self.assertEqual(result,expected)
        
    # An ionic formula must end with either a number, a '+', or a '-'  
    def test_is_valid_formula_5(self):
        input = 'HCO3-'
        result = cf.is_valid_formula(input)        
        expected = True
        
        self.assertEqual(result,expected)
    
    # An ionic formula must end with either a number, a '+', or a '-'  
    def test_is_valid_formula_6(self):
        input = 'Fe++'
        result = cf.is_valid_formula(input)        
        expected = True
        
        self.assertEqual(result,expected)
    
    # An ionic formula must end with either a number, a '+', or a '-'  
    def test_is_valid_formula_7(self):
        input = 'Mg+2'
        result = cf.is_valid_formula(input)        
        expected = True
        
        self.assertEqual(result,expected)
    
    # An ionic formula must end with either a number, a '+', or a '-'  
    def test_is_valid_formula_8(self):
        input = 'V+5+'
        result = cf.is_valid_formula(input)        
        expected = False
        
        self.assertEqual(result,expected)
    
    # Formulas must contain only valid atomic symbols that start with capital letters
    def test_is_valid_formula_9(self):
        input = 'NaOH'
        result = cf.is_valid_formula(input)        
        expected = True
        
        self.assertEqual(result,expected)
    
    # Formulas must contain only valid atomic symbols that start with capital letters
    def test_is_valid_formula_9(self):
        input = 'naOH'
        result = cf.is_valid_formula(input)        
        expected = False
        
        self.assertEqual(result,expected)
    
    # Formulas must contain only valid atomic symbols that start with capital letters
    def test_is_valid_formula_10(self):
        input = 'HzCl'
        result = cf.is_valid_formula(input)        
        expected = False
        
        self.assertEqual(result,expected)
    
    # A formula with parentheses must have the same number of '(' and ')'
    def test_is_valid_formula_11(self):
        input = '(NH3)2(NO3)2'
        result = cf.is_valid_formula(input)        
        expected = True
        
        self.assertEqual(result,expected)
    
    # A formula with parentheses must have the same number of '(' and ')'
    def test_is_valid_formula_12(self):
        input = 'Mg)(OH)2'
        result = cf.is_valid_formula(input)        
        expected = False
        
        self.assertEqual(result,expected)
    
    # A formula cannot end with an open parenthesis
    def test_is_valid_formula_13(self):
        input = 'Na+('
        result = cf.is_valid_formula(input)        
        expected = False
        
        self.assertEqual(result,expected)
    
    # An open parenthesis must always be followed by an atomic symbol
    def test_is_valid_formula_14(self):
        input = '(3)Na+'
        result = cf.is_valid_formula(input)        
        expected = False
        
        self.assertEqual(result,expected)
    
    # An open parenthesis must always be followed by an atomic symbol
    def test_is_valid_formula_15(self):
        input = '()'
        result = cf.is_valid_formula(input)        
        expected = False
        
        self.assertEqual(result,expected)
    
    # A closed parenthesis may be followed by a number or used to designate a group
    def test_is_valid_formula_16(self):
        input = 'CH3C(O)CH3'
        result = cf.is_valid_formula(input)        
        expected = True
        
        self.assertEqual(result,expected)
    
    # A closed parenthesis may be followed by a number or used to designate a group
    def test_is_valid_formula_17(self):
        input = 'Mg(OH)2'
        result = cf.is_valid_formula(input)        
        expected = True
        
        self.assertEqual(result,expected)
    
    # An open parenthesis must always precede the nearest closed parenthesis
    def test_is_valid_formula_17(self):
        input = 'CH3(CH)2(CH)2'
        result = cf.is_valid_formula(input)        
        expected = True
        
        self.assertEqual(result,expected)
    
    # An open parenthesis must always precede the nearest closed parenthesis
    def test_is_valid_formula_18(self):
        input = ')Na+('
        result = cf.is_valid_formula(input)        
        expected = False
        
        self.assertEqual(result,expected)

class Test_consolidate_formula(unittest.TestCase):
    '''
    tests for _consolidate_formula()
    --------------------------------
    This function is important for calculating molecular weights correctly
    '''
    
    def test_is_valid_formula_1(self):
        input = 'Fe2(SO4)4'
        result = cf._consolidate_formula(input)        
        expected = ['Fe', 2, 'S', 4, 'O', 16]
        
        self.assertEqual(result,expected)
    
    def test_is_valid_formula_2(self):
        input = '(NH4)3PO4'
        result = cf._consolidate_formula(input)        
        expected = ['N', 3, 'H', 12, 'P', 1, 'O', 4]
        
        self.assertEqual(result,expected)
    
    def test_is_valid_formula_3(self):
        input = 'CH3(CH2)6CH3'
        result = cf._consolidate_formula(input)        
        expected = ['C', 8, 'H', 18]
        
        self.assertEqual(result,expected)


if __name__ == '__main__':
    unittest.main()
