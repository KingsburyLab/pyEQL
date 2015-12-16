'''
pyEQL osmotic coefficient test suite
============================================

This file contains tests for the osmotic coefficient method
employed by pyEQL.

NOTE: generally, these tests check the module output against experimental
data rather than the theoretical result of the respective functions.
'''

import pyEQL
import unittest

class Test_osmotic_pitzer(unittest.TestCase):
    '''
    test osmotic coefficient based on the Pitzer model
    ------------------------------------------------

    '''
    def setUp(self):
        self.s1 = pyEQL.Solution([['Na+','0.1 mol/L'],['Cl-','0.1 mol/L']])
        
    def test_osmotic_pitzer_coeff_units(self):
        # the osmotic coefficient should be dimensionless
        result = self.s1.get_osmotic_coefficient().dimensionality
        self.assertEqual(result,'')

    def test_activity_pitzer_magnitude(self):
        # the osmotic coefficient should be greater than zero
        result = self.s1.get_osmotic_coefficient()
        self.assertGreaterEqual(result,0)

    def test_osmotic_pitzer_ammoniumnitrate(self):
        '''        
        calculate the osmotic coefficient at each concentration and compare
        to experimental data for ammonium nitrate
        
        References
        ----------
        May, P. M., Rowland, D., Hefter, G., & Königsberger, E. (2011). 
        A Generic and Updatable Pitzer Characterization of Aqueous Binary Electrolyte Solutions at 1 bar and 25 °C. 
        Journal of Chemical & Engineering Data, 56(12), 5066–5077. doi:10.1021/je2009329
    
        '''
        # list of concentrations to test, mol/kg
        conc_list = [0.25,0.5,0.75,1,1.5,2]
    
        # list of published experimental osmotic coefficients
        pub_osmotic_coeff = [0.86,0.855,0.83,0.825,0.80,0.78]
    
        for i in range(len(conc_list)):
            with self.subTest(conc=conc_list[i]):
                conc = str(conc_list[i]) + 'mol/kg'
                sol = pyEQL.Solution()
                sol.add_solute('NH4+',conc)
                sol.add_solute('NO3-',conc)
                result=sol.get_osmotic_coefficient()
                expected = pub_osmotic_coeff[i]
                
                self.assertAlmostEqual(result,expected,1)
    
    def test_osmotic_pitzer_coppersulfate(self):
        '''        
        calculate the osmotic coefficient at each concentration and compare
        to experimental data for copper sulate
        
        References
        ----------
        May, P. M., Rowland, D., Hefter, G., & Königsberger, E. (2011). 
        A Generic and Updatable Pitzer Characterization of Aqueous Binary Electrolyte Solutions at 1 bar and 25 °C. 
        Journal of Chemical & Engineering Data, 56(12), 5066–5077. doi:10.1021/je2009329
    
        '''
        # list of concentrations to test, mol/kg
        conc_list = [0.25,0.5,0.75,1]
    
        # list of published experimental osmotic coefficients
        pub_osmotic_coeff = [0.5,0.485,0.48,0.485,0.5]
    
        for i in range(len(conc_list)):
            with self.subTest(conc=conc_list[i]):
                conc = str(conc_list[i]) + 'mol/kg'
                sol = pyEQL.Solution()
                sol.add_solute('Cu+2',conc)
                sol.add_solute('SO4-2',conc)
                result=sol.get_osmotic_coefficient()
                expected = pub_osmotic_coeff[i]
                
                self.assertAlmostEqual(result,expected,1)        

if __name__ == '__main__':
    unittest.main()
