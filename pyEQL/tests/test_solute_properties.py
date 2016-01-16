'''
pyEQL solute properties test suite
============================================

This file contains tests of the Solution class methods that retrieve or
calculate properties of individual solutes. Methods currently included
in the testing are:

- get_transport_number()
- get_molar_conductivity()
- get_mobility()


'''

import pyEQL
import unittest

class Test_transport_number(unittest.TestCase):
    '''
    test get_transport_number
    ------------------------------------------------
    '''
    # create a test solution
    def setUp(self):
        self.s1 = pyEQL.Solution([['K+','0.1 mol/L'],['Cl-','0.1 mol/L'],['FeO','0.2 mol/L']])
        
    # The transport number of any ion should be between 0 and 1
    def test_transport_number_1(self):
        actual = self.s1.get_transport_number('K+')
        
        self.assertTrue((0 <= actual) and (actual <= 1))
    
    # The transport number of water should be 0
    def test_transport_number_2(self):
        actual = self.s1.get_transport_number('H2O')
        
        self.assertEqual(actual.magnitude,0)
        
    # The transport number of any uncharged species should be 0
    def test_transport_number_3(self):
        actual = self.s1.get_transport_number('FeO')
        
        self.assertEqual(actual.magnitude,0)
    
    # The transport numbers should add up to 1
    def test_transport_number_4(self):
        total_t = 0
        for item in self.s1.components:
            total_t += self.s1.get_transport_number(item)
        
        self.assertAlmostEqual(total_t.magnitude,1)
    

class Test_molar_conductivity(unittest.TestCase):
    '''
    test get_molar_conductivity
    
    Published molar conductivity values at 25 degC and infinite dilution
    are found in CRC Handbook of Chemistry and Physics, 92nd ed., 2011.
    
    This reference is also the source of diffusion coefficient data
    in the default database, so the values molar conductivity (which
    is calculated from the diffusion coefficient) should match.
    ----------------------------------------------------------------
    '''
    # decimal places for AlmostEqual
    def setUp(self):
        self.places = 5
    
    def test_molar_conductivity_potassium(self):
        # K+ - 73.48 x 10 ** -4 m ** 2 S / mol
        self.s1 = pyEQL.Solution([['K+','0.001 mol/L'],['Cl-','0.001 mol/L']],temperature='25 degC')
        actual = self.s1.get_molar_conductivity('K+').to('m**2*S/mol')
        expected = pyEQL.unit('73.48e-4 m**2 * S / mol')
        result = (actual-expected).magnitude
        
        self.assertAlmostEqual(result,0,self.places)
        
    def test_molar_conductivity_sodium(self):
        # Na+ - 50.08 x 10 ** -4 m ** 2 S / mol
        self.s1 = pyEQL.Solution([['Na+','0.001 mol/L'],['Cl-','0.001 mol/L']],temperature='25 degC')
        actual = self.s1.get_molar_conductivity('Na+').to('m**2*S/mol')
        expected = pyEQL.unit('50.08e-4 m**2 * S / mol')
        result = (actual-expected).magnitude
        
        self.assertAlmostEqual(result,0,self.places)
    
    def test_molar_conductivity_magnesium(self):
        # Mg+2 - 106 x 10 ** -4 m ** 2 S / mol
        self.s1 = pyEQL.Solution([['Mg+2','0.001 mol/L'],['Cl-','0.002 mol/L']],temperature='25 degC')
        actual = self.s1.get_molar_conductivity('Mg+2').to('m**2*S/mol')
        expected = pyEQL.unit('106e-4 m**2 * S / mol')
        result = (actual-expected).magnitude
        
        self.assertAlmostEqual(result,0,self.places)
        
    def test_molar_conductivity_chloride(self):
        # Cl- - 76.31 x 10 ** -4 m ** 2 S / mol
        self.s1 = pyEQL.Solution([['Na+','0.001 mol/L'],['Cl-','0.001 mol/L']],temperature='25 degC')
        actual = self.s1.get_molar_conductivity('Cl-').to('m**2*S/mol')
        expected = pyEQL.unit('76.31e-4 m**2 * S / mol')
        result = (actual-expected).magnitude
        
        self.assertAlmostEqual(result,0,self.places)

    def test_molar_conductivity_fluoride(self):
        # F- - 55.4 x 10 ** -4 m ** 2 S / mol
        self.s1 = pyEQL.Solution([['Na+','0.001 mol/L'],['F-','0.001 mol/L']],temperature='25 degC')
        actual = self.s1.get_molar_conductivity('F-').to('m**2*S/mol')
        expected = pyEQL.unit('55.4e-4 m**2 * S / mol')
        result = (actual-expected).magnitude
        
        self.assertAlmostEqual(result,0,self.places)
        
    def test_molar_conductivity_sulfate(self):
        # SO4-2 - 160 x 10 ** -4 m ** 2 S / mol
        self.s1 = pyEQL.Solution([['Na+','0.002 mol/L'],['SO4-2','0.001 mol/L']],temperature='25 degC')
        actual = self.s1.get_molar_conductivity('SO4-2').to('m**2*S/mol')
        expected = pyEQL.unit('160.0e-4 m**2 * S / mol')
        result = (actual-expected).magnitude
        
        self.assertAlmostEqual(result,0,self.places)
        
    def test_molar_conductivity_hydroxide(self):
        # OH- - 198 x 10 ** -4 m ** 2 S / mol
        self.s1 = pyEQL.Solution(temperature='25 degC')
        actual = self.s1.get_molar_conductivity('OH-').to('m**2*S/mol')
        expected = pyEQL.unit('198e-4 m**2 * S / mol')
        result = (actual-expected).magnitude
        
        self.assertAlmostEqual(result,0,self.places)
    
    def test_molar_conductivity_hydrogen(self):
        # H+ - 349.65 x 10 ** -4 m ** 2 S / mol
        self.s1 = pyEQL.Solution(temperature='25 degC')
        actual = self.s1.get_molar_conductivity('H+').to('m**2*S/mol')
        expected = pyEQL.unit('349.65e-4 m**2 * S / mol')
        result = (actual-expected).magnitude
        
        self.assertAlmostEqual(result,0,self.places)
    
    # molar conductivity of a neutral solute should be zero
    def test_molar_conductivity_neutral(self):
        self.s1 = pyEQL.Solution([['FeCl3','0.001 mol/L']],temperature='25 degC')
        actual = self.s1.get_molar_conductivity('FeCl3').to('m**2*S/mol')
        expected = pyEQL.unit('0 m**2 * S / mol')
        result = (actual-expected).magnitude
        
        self.assertAlmostEqual(result,0,self.places)
        
    # molar conductivity of water should be zero
    def test_molar_conductivity_water(self):
        self.s1 = pyEQL.Solution(temperature='25 degC')
        actual = self.s1.get_molar_conductivity('H2O').to('m**2*S/mol')
        expected = pyEQL.unit('0 m**2 * S / mol')
        result = (actual-expected).magnitude
        
        self.assertAlmostEqual(result,0,self.places)

class Test_mobility(unittest.TestCase):
    '''
    test get_mobility
    
    Since the mobility is simply the molar conductivity divided by
    Faraday's constant and the ion charge, and we have already
    tested get_molar_conductivity, here we only do a couple of tests
    to make sure get_mobility returns a value that is correctly
    scaled related to get_molar_conductivity
    ----------------------------------------------------------------
    '''
    # decimal places for AlmostEqual
    def setUp(self):
        self.places = 5
    
    def test_mobility_potassium(self):
        self.s1 = pyEQL.Solution([['K+','0.001 mol/L'],['Cl-','0.001 mol/L']],temperature='25 degC')
        molar_conductivity = self.s1.get_molar_conductivity('K+').to('m**2*S/mol')
        mobility = self.s1.get_mobility('K+')
        charge = self.s1.get_solute('K+').get_formal_charge()
        # calculate the mobility from get_molar_conductivity, then compare with get_mobility
        result = (molar_conductivity / (pyEQL.unit.N_A * pyEQL.unit.e * abs(charge)) - mobility).magnitude
        
        self.assertAlmostEqual(result,0,self.places)
        
    def test_mobility_chloride(self):
        self.s1 = pyEQL.Solution([['K+','0.001 mol/L'],['Cl-','0.001 mol/L']],temperature='25 degC')
        molar_conductivity = self.s1.get_molar_conductivity('Cl-').to('m**2*S/mol')
        mobility = self.s1.get_mobility('Cl-')
        charge = self.s1.get_solute('Cl-').get_formal_charge()
        # calculate the mobility from get_molar_conductivity, then compare with get_mobility
        result = (molar_conductivity / (pyEQL.unit.N_A * pyEQL.unit.e * abs(charge)) - mobility).magnitude
        
        self.assertAlmostEqual(result,0,self.places)
    
    def test_mobility_magnesium(self):
        self.s1 = pyEQL.Solution([['Mg+2','0.001 mol/L'],['Cl-','0.002 mol/L']],temperature='25 degC')
        molar_conductivity = self.s1.get_molar_conductivity('Mg+2').to('m**2*S/mol')
        mobility = self.s1.get_mobility('Mg+2')
        charge = self.s1.get_solute('Mg+2').get_formal_charge()
        # calculate the mobility from get_molar_conductivity, then compare with get_mobility
        result = (molar_conductivity / (pyEQL.unit.N_A * pyEQL.unit.e * abs(charge)) - mobility).magnitude
        
        self.assertAlmostEqual(result,0,self.places)
        
    def test_mobility_sulfate(self):
        self.s1 = pyEQL.Solution([['K+','0.002 mol/L'],['SO4-2','0.001 mol/L']],temperature='25 degC')
        molar_conductivity = self.s1.get_molar_conductivity('SO4-2').to('m**2*S/mol')
        mobility = self.s1.get_mobility('SO4-2')
        charge = self.s1.get_solute('SO4-2').get_formal_charge()
        # calculate the mobility from get_molar_conductivity, then compare with get_mobility
        result = (molar_conductivity / (pyEQL.unit.N_A * pyEQL.unit.e * abs(charge)) - mobility).magnitude
        
        self.assertAlmostEqual(result,0,self.places)
        
if __name__ == '__main__':
    unittest.main()
