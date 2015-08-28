'''
pyEQL volume and concentration methods test suite
=================================================

This file contains tests for the volume and concentration-related methods
used by pyEQL's Solution class
'''

import pyEQL
import unittest

class Test_empty_solution(unittest.TestCase):
    '''
    test behavior when creating an empty solution
    ------------------------------------------------
    
    '''
    # create an empty solution
    def setUp(self):
        self.s1 = pyEQL.Solution()
        
    # It should return type Solution
    def test_empty_solution_1(self):
        expected = pyEQL.solution.Solution
        
        self.assertIsInstance(self.s1,expected)
    
    # It should have exactly 1L volume
    def test_empty_solution_2(self):
        result = self.s1.get_volume().to('L').magnitude
        expected = 1.0
        
        self.assertEqual(result,expected)
    
    #  the solvent should be water
    def test_empty_solution_3(self):
        result = self.s1.get_solvent().get_name()
        expected = 'H2O'
        
        self.assertEqual(result,expected)
    
    # It should have 0.997 kg water mass
    def test_empty_solution_4(self):
        result = self.s1.get_solvent_mass().to('kg').magnitude
        expected = 0.9970415
        
        self.assertAlmostEqual(result,expected)
    
    # the temperature should be 25 degC
    def test_empty_solution_5(self):
        result = self.s1.get_temperature().to('degC').magnitude
        expected = 25
        
        self.assertEqual(result,expected)
    
    # the pressure should be 1 atm
    def test_empty_solution_6(self):
        result = self.s1.get_pressure().to('atm').magnitude
        expected = 1
        
        self.assertEqual(result,expected)
    
    # the pH should be 7.0
    def test_empty_solution_7(self):
        result = self.s1.get_activity('H+')
        expected = 1e-7
        
        self.assertAlmostEqual(result,expected,9)
    
    # it should contain H2O, H+, and OH- species
    def test_empty_solution_8(self):
        result = self.s1.list_solutes()
        expected = ['H2O', 'OH-', 'H+']
        
        self.assertCountEqual(result,expected)
    
class Test_solute_addition(unittest.TestCase):
    '''
    test behavior of various methods for adding solutes to a solution
    -----------------------------------------------------------------
    
    '''
    # create an empty and test solutions with the same volume using substance / volume,
    # substance/mass, and substance units
    def setUp(self):
        self.s1 = pyEQL.Solution(volume='2 L')
        self.s2 = pyEQL.Solution([['Na+','4 mol/L'],['Cl-','4 mol/L']],volume='2 L')
        self.s3 = pyEQL.Solution([['Na+','4 mol/kg'],['Cl-','4 mol/kg']],volume='2 L')
        self.s4 = pyEQL.Solution([['Na+','8 mol'],['Cl-','8 mol']],volume='2 L')
        
    # if solutes are added at creation-time with substance / volume units, 
    # then the total volume of the solution should not change (should remain at 2 L)
    def test_solute_addition_1(self):
        result = self.s2.get_volume().to('L').magnitude
        expected = 2
        
        self.assertEqual(result,expected)
    
    # if solutes are added at creation-time with substance / volume units, 
    # then the resulting mol/L concentrations should be exactly what was specified
    def test_solute_addition_2(self):
        result = self.s2.get_amount('Na+','mol/L').magnitude
        expected = 4
        
        self.assertEqual(result,expected)
        
    # if solutes are added at creation-time with substance / mass units, 
    # then the resulting mol/kg concentrations should be exactly what was specified
    def test_solute_addition_3(self):
        result = self.s3.get_amount('Na+','mol/kg').magnitude
        expected = 4
        
        self.assertEqual(result,expected)
        
    # the water mass of solution s2 should be less than that of s3, because
    # of the volume recalculation
    def test_solute_addition_4(self):
        result_molL = self.s2.get_solvent_mass().to('kg').magnitude
        result_molkg = self.s3.get_solvent_mass().to('kg').magnitude
        
        self.assertLess(result_molL,result_molkg)
    
    # if solutes are added at creation-time with substance units, 
    # then the resulting mol amounts should be exactly what was specified
    def test_solute_addition_4a(self):
        result = self.s4.get_amount('Na+','mol').magnitude
        expected = 8
        
        self.assertEqual(result,expected)
        
    # the water mass of solution s2 should be less than that of s4, because
    # of the volume recalculation
    def test_solute_addition_4b(self):
        result_molL = self.s2.get_solvent_mass().to('kg').magnitude
        result_mol = self.s4.get_solvent_mass().to('kg').magnitude
        
        self.assertLess(result_molL,result_mol)
    
    '''
    Tests for set_amount() method
    '''
    
    # If the concentration of a solute is directly set with a substance / volume
    # unit, the volume should not change
    def test_solute_addition_5(self):
        self.s2.set_amount('Na+','5 mol/L')
        self.s2.set_amount('Cl-','5 mol/L')
        result = self.s2.get_volume().to('L').magnitude
        expected = 2
        
        self.assertEqual(result,expected)
    
    # If the concentration of a solute is directly set with a substance / volume
    # unit, the water mass should be reduced
    def test_solute_addition_6(self):
        original = self.s2.get_solvent_mass().to('kg').magnitude
        self.s2.set_amount('Na+','5 mol/L')
        self.s2.set_amount('Cl-','5 mol/L')
        result = self.s2.get_solvent_mass().to('kg').magnitude
        
        self.assertLess(result,original)
    
    # If the concentration of a solute is directly set with a substance / volume
    # unit, the resulting concentration should be exactly what was specified
    def test_solute_addition_7(self):
        self.s2.set_amount('Na+','5 mol/L')
        self.s2.set_amount('Cl-','5 mol/L')
        result = self.s2.get_amount('Na+','mol/L').magnitude
        expected = 5
        
        self.assertEqual(result,expected)
    
    # If the concentration of a solute is directly set with a substance / mass
    # unit, the volume should increase
    def test_solute_addition_8(self):
        original = self.s2.get_volume().to('L').magnitude        
        self.s2.set_amount('Na+','5 mol/kg')
        self.s2.set_amount('Cl-','5 mol/kg')
        result = self.s2.get_volume().to('L').magnitude
        
        self.assertGreater(result,original)
    
    # If the concentration of a solute is directly set with a substance / mass
    # unit, the water mass should not change
    def test_solute_addition_9(self):
        original = self.s2.get_solvent_mass().to('kg').magnitude
        self.s2.set_amount('Na+','5 mol/kg')
        self.s2.set_amount('Cl-','5 mol/kg')
        result = self.s2.get_solvent_mass().to('kg').magnitude
        
        self.assertEqual(result,original)
    
    # If the concentration of a solute is directly set with a substance / mass
    # unit, the resulting concentration should be exactly what was specified
    def test_solute_addition_10(self):
        self.s2.set_amount('Na+','5 mol/kg')
        self.s2.set_amount('Cl-','5 mol/kg')
        result = self.s2.get_amount('Na+','mol/kg').magnitude
        expected = 5
        
        self.assertEqual(result,expected)
    
    # If the concentration of a solute is directly set with a substance
    # unit, the volume should increase
    def test_solute_addition_8a(self):
        original = self.s2.get_volume().to('L').magnitude        
        self.s2.set_amount('Na+','10 mol')
        self.s2.set_amount('Cl-','10 mol')
        result = self.s2.get_volume().to('L').magnitude
        
        self.assertGreater(result,original)
    
    # If the concentration of a solute is directly set with a substance
    # unit, the water mass should not change
    def test_solute_addition_9a(self):
        original = self.s2.get_solvent_mass().to('kg').magnitude
        self.s2.set_amount('Na+','10 mol')
        self.s2.set_amount('Cl-','10 mol')
        result = self.s2.get_solvent_mass().to('kg').magnitude
        
        self.assertEqual(result,original)
    
    # If the concentration of a solute is directly set with a substance / mass
    # unit, the resulting concentration should be exactly what was specified
    def test_solute_addition_10a(self):
        self.s2.set_amount('Na+','10 mol')
        self.s2.set_amount('Cl-','10 mol')
        result = self.s2.get_amount('Na+','mol').magnitude
        expected = 10
        
        self.assertEqual(result,expected)
    
    '''
    Tests for add_amount() method
    '''
    
    ## substance / volume units
    
    # If the concentration of a solute is directly increased with a substance / volume
    # unit, the volume should not change
    def test_solute_addition_11(self):
        self.s2.add_amount('Na+','1 mol/L')
        self.s2.add_amount('Cl-','1 mol/L')
        result = self.s2.get_volume().to('L').magnitude
        expected = 2
        
        self.assertEqual(result,expected)
    
    # If the concentration of a solute is directly increased with a substance / volume
    # unit, the water mass should be reduced
    def test_solute_addition_12(self):
        original = self.s2.get_solvent_mass().to('kg').magnitude
        self.s2.add_amount('Na+','1 mol/L')
        self.s2.add_amount('Cl-','1 mol/L')
        result = self.s2.get_solvent_mass().to('kg').magnitude
        
        self.assertLess(result,original)
    
    # If the concentration of a solute is directly increased with a substance / volume
    # unit, the resulting concentration should be exactly what was specified
    def test_solute_addition_13(self):
        self.s2.add_amount('Na+','1 mol/L')
        self.s2.add_amount('Cl-','1 mol/L')
        result = self.s2.get_amount('Na+','mol/L').magnitude
        expected = 5
        
        self.assertEqual(result,expected)
    
    ## substance / mass units    
    
    # If the concentration of a solute is directly increased with a substance / mass
    # unit, the volume should increase
    def test_solute_addition_14(self):
        original = self.s3.get_volume().to('L').magnitude        
        self.s3.add_amount('Na+','1 mol/kg')
        self.s3.add_amount('Cl-','1 mol/kg')
        result = self.s3.get_volume().to('L').magnitude
        
        self.assertGreater(result,original)
    
    # If the concentration of a solute is directly increased with a substance / mass
    # unit, the water mass should not change
    def test_solute_addition_15(self):
        original = self.s3.get_solvent_mass().to('kg').magnitude
        self.s3.add_amount('Na+','1 mol/kg')
        self.s3.add_amount('Cl-','1 mol/kg')
        result = self.s3.get_solvent_mass().to('kg').magnitude
        
        self.assertEqual(result,original)
    
    # If the concentration of a solute is directly increased with a substance / mass
    # unit, the resulting concentration should be exactly what was specified
    def test_solute_addition_16(self):
        self.s3.add_amount('Na+','1 mol/kg')
        self.s3.add_amount('Cl-','1 mol/kg')
        result = self.s3.get_amount('Na+','mol/kg').magnitude
        expected = 5
        
        self.assertEqual(result,expected)
    
    ## substance units
    
    # If the concentration of a solute is directly increased with a substance
    # unit, the volume should increase
    def test_solute_addition_14a(self):
        original = self.s2.get_volume().to('L').magnitude        
        self.s2.add_amount('Na+','2 mol')
        self.s2.add_amount('Cl-','2 mol')
        result = self.s2.get_volume().to('L').magnitude
        
        self.assertGreater(result,original)
    
    # If the concentration of a solute is directly increased with a substance
    # unit, the water mass should not change
    def test_solute_addition_15a(self):
        original = self.s2.get_solvent_mass().to('kg').magnitude
        self.s2.add_amount('Na+','2 mol')
        self.s2.add_amount('Cl-','2 mol')
        result = self.s2.get_solvent_mass().to('kg').magnitude
        
        self.assertEqual(result,original)
    
    # If the concentration of a solute is directly increased with a substance
    # unit, the resulting concentration should be exactly what was specified
    def test_solute_addition_16a(self):
        self.s2.add_amount('Na+','2 mol')
        self.s2.add_amount('Cl-','2 mol')
        result = self.s2.get_amount('Na+','mol').magnitude
        expected = 10
        
        self.assertEqual(result,expected)
    
    ## negative substance units
    # If the concentration of a solute is directly decreased with a substance
    # unit, the volume should decrease
    def test_solute_addition_14b(self):
        original = self.s2.get_volume().to('L').magnitude        
        self.s2.add_amount('Na+','-2 mol')
        self.s2.add_amount('Cl-','-2 mol')
        result = self.s2.get_volume().to('L').magnitude
        
        self.assertLess(result,original)
    
    # If the concentration of a solute is directly changed with a substance
    # unit, the water mass should not change
    def test_solute_addition_15b(self):
        original = self.s2.get_solvent_mass().to('kg').magnitude
        self.s2.add_amount('Na+','-2 mol')
        self.s2.add_amount('Cl-','-2 mol')
        result = self.s2.get_solvent_mass().to('kg').magnitude
        
        self.assertEqual(result,original)
    
    # If the concentration of a solute is directly changed with a substance
    # unit, the resulting concentration should be exactly what was specified
    def test_solute_addition_16b(self):
        self.s2.add_amount('Na+','-2 mol')
        self.s2.add_amount('Cl-','-2 mol')
        result = self.s2.get_amount('Na+','mol').magnitude
        expected = 6
        
        self.assertEqual(result,expected)    
    
class Test_get_amount(unittest.TestCase):
    '''
    test the get_amount() method on a 1 mol/L NaCl solution
    ----------------------------
    1 mol NaCl / L = 58.44 g/L
    Na+ = 22.98977 g/mol
    
    '''
    # create the 1 M NaCl solution
    def setUp(self):
        self.s1 = pyEQL.Solution([['Na+','1 mol/L'],['Cl-','1 mol/L']])
            
    # get_amount() - mol/L
    def test_get_amount_molL(self):
        result = self.s1.get_amount('Na+','mol/L').magnitude
        expected = 1
        
        self.assertAlmostEqual(result,expected,9)
    
    # get_amount() - mol/kg
    def test_get_amount_molkg(self):
        result = self.s1.get_amount('Na+','mol/kg').magnitude
        expected = 1.02181221888
        
        self.assertAlmostEqual(result,expected,9)
        
    # get_amount() - g/L
    def test_get_amount_gL(self):
        result = self.s1.get_amount('Na+','g/L').magnitude
        expected = 22.98977
        
        self.assertAlmostEqual(result,expected,9)
    
    # get_amount() - mg
    def test_get_amount_mg(self):
        result = self.s1.get_amount('Na+','mg').magnitude
        expected = 22989.77
        
        self.assertAlmostEqual(result,expected,9)
    
    # get_amount() - mol
    def test_get_amount_mol(self):
        result = self.s1.get_amount('Na+','mol').magnitude
        expected = 1
        
        self.assertAlmostEqual(result,expected,9)
    
    # get_amount() - fraction
    def test_get_amount_fraction(self):
        result = self.s1.get_amount('Na+','fraction')
        expected = 0.01775457254
        
        self.assertAlmostEqual(result,expected,9)    
    
    
    
    
if __name__ == '__main__':
    unittest.main()
