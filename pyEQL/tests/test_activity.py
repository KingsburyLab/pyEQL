'''
pyEQL activity correction methods test suite
============================================

This file contains tests for some of the activity correction methods
employed by pyEQL.

NOTE: generally, these tests check the module output against experimental
data rather than the theoretical result of the respective functions. In some
cases, the output is also tested against a well-established model published
by USGS(PHREEQC)
'''

import pyEQL
import unittest

class Test_activity_pitzer_nacl(unittest.TestCase):
    '''
    test Pitzer model for activity of NaCl
    ------------------------------------------------
    '''
    
    def setUp(self):
        self.s1 = pyEQL.Solution([['Na+','0.1 mol/L'],['Cl-','0.1 mol/L']])
        # list of molal concentrations for published activity coefficients

    def test_activity_pitzer_coeff_units(self):
        # the activity coefficient should be a float
        result = self.s1.get_activity_coefficient('Na+')
        self.assertIsInstance(result,float)

    def test_activity_pitzer_units(self):
        # the activity should be a float
        result = self.s1.get_activity('Na+')
        self.assertIsInstance(result,float)          
        
    def test_activity_pitzer_equality(self):
        # the activity coefficient of both the Na+ and Cl- should be the same
        a1 = self.s1.get_activity_coefficient('Na+')
        a2 = self.s1.get_activity_coefficient('Cl-')
        self.assertEqual(a1,a2)
        
    def test_activity_pitzer_nacl_1(self):
        '''        
        calculate the activity coefficient at each concentration and compare
        to experimental data
                    
        Experimental activity coefficient values at 25 degC are found in 
        J. Phys. Chem. Reference Data Vol 13 (1), 1984, p.53.
    
        '''
        # list of concentrations to test, mol/kg
        conc_list = [0.1,0.25,0.5,0.75,1,2,3,4,5,6]
    
        # list of published experimental activity coefficients
        pub_activity_coeff = [0.778,0.72,0.681,0.665,0.657,0.669,0.714,0.782,0.873,0.987]
    
        for i in range(len(conc_list)):
            with self.subTest(conc=conc_list[i]):
                conc = str(conc_list[i]) + 'mol/kg'
                sol = pyEQL.Solution()
                sol.add_solute('Na+',conc)
                sol.add_solute('Cl-',conc)
                result=sol.get_activity_coefficient('Na+')
                expected = pub_activity_coeff[i]
                
                self.assertAlmostEqual(result,expected,2)
                
    # The pitzer model diverges a bit from experimental data at high concentration
    @unittest.expectedFailure
    def test_water_activity_pitzer_nacl_1(self):
        '''        
        calculate the water activity at each concentration and compare
        to experimental data
                    
        Experimental osmotic coefficients for NaCl are found in:
        Pitzer and Pelper, 1984. "Thermodyamic Properties of Aqueous Sodium Chloride Solutions"
        J. Phys. Chem. Ref. Data 13(1).
        
        Osmotic coefficients were converted into water activity according to the equation
        found in 
        Blandamer, Mike J., Engberts, Jan B. F. N., Gleeson, Peter T., Reis, 
        Joao Carlos R., 2005. "Activity of water in aqueous systems: A frequently 
        neglected property." //Chemical Society Review// 34, 440-458.
        
        $$ ln a_w = - \Phi M_w \sum_i m_i    $$
                
        Where M_w is the molar mass of water (0.018015 kg/mol) and m_i is the molal concentration
        of each species.
    
        '''
        # list of concentrations to test, mol/kg
        conc_list = [0.1,0.25,0.5,0.75,1,2,3,4,5,6]
    
        # list of published experimental water activity
        pub_water_activity = [0.9964,0.9910,0.9821,0.9733,0.9646,0.9304,0.8975,0.8657,0.8351,0.8055]
    
        for i in range(len(conc_list)):
            with self.subTest(conc=conc_list[i]):
                conc = str(conc_list[i]) + 'mol/kg'
                sol = pyEQL.Solution()
                sol.add_solute('Na+',conc)
                sol.add_solute('Cl-',conc)
                result=sol.get_water_activity()
                expected = pub_water_activity[i]
                
                self.assertAlmostEqual(result,expected,2)
        
    def test_activity_pitzer_phreeqc_nacl_2(self):
        '''        
        calculate the activity coefficient at each concentration and compare
        to the output of the PHREEQC model
                    
        PHREEQC version 3.1.4 was used to calculate density, conductivity, water
        activity, and NaCl activity coefficient for NaCl solutions up to 6m. 
        The Pitzer model (pitzer.dat) database was used.
        <http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/>
    
        '''
        # list of concentrations to test, mol/kg
        conc_list = [0.1,0.25,0.5,0.75,1,2,3,4,5,6]
    
        # list of modeled activity coefficients
        phreeqc_pitzer_activity_coeff = [0.7771,0.7186,0.6799,0.6629,0.6561,0.6675,0.7128,0.7825,0.8729,0.9874]
    
        for i in range(len(conc_list)):
            with self.subTest(conc=conc_list[i]):
                conc = str(conc_list[i]) + 'mol/kg'
                sol = pyEQL.Solution()
                sol.add_solute('Na+',conc)
                sol.add_solute('Cl-',conc)
                result=sol.get_activity_coefficient('Na+')
                expected = phreeqc_pitzer_activity_coeff[i]
                
                self.assertAlmostEqual(result,expected,2)
    
            
    def test_water_activity_phreeqc_pitzer_nacl_2(self):
        '''        
        calculate the water activity at each concentration and compare
        to the output of the PHREEQC model
                    
        PHREEQC version 3.1.4 was used to calculate density, conductivity, water
        activity, and NaCl activity coefficient for NaCl solutions up to 6m. 
        The Pitzer model (pitzer.dat) database was used.
        <http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/>
    
        '''
        # list of concentrations to test, mol/kg
        conc_list = [0.1,0.25,0.5,0.75,1,2,3,4,5,6]
    
        # list of modeled water activities
        phreeqc_pitzer_water_activity = [0.997,0.992,0.984,0.975,0.967,0.932,0.893,0.851,0.807,0.759]
    
        for i in range(len(conc_list)):
            with self.subTest(conc=conc_list[i]):
                conc = str(conc_list[i]) + 'mol/kg'
                sol = pyEQL.Solution()
                sol.add_solute('Na+',conc)
                sol.add_solute('Cl-',conc)
                result=sol.get_water_activity()
                expected = phreeqc_pitzer_water_activity[i]
                
                self.assertAlmostEqual(result,expected,2)


if __name__ == '__main__':
    unittest.main()
