'''
pyEQL dielectric constant test suite
============================================

This file contains tests that check the dielectric constant
computations of pyEQL

'''

import pyEQL
import unittest

class Test_dielectric(unittest.TestCase,pyEQL.CustomAssertions):
    '''
    test the Dielectric Constant calculations of various solutions
    ------------------------------------------------
    
    Reference: A. Zuber, L. Cardozo-Filho, V.F. Cabral, R.F. Checoni, M. Castier, An empirical equation for the dielectric constant in aqueous and 
    nonaqueous electrolyte mixtures, Fluid Phase Equilib. 376 (2014) 116–123. doi:10.1016/j.fluid.2014.05.037.
    
    '''
    def setUp(self):
        # relative error tolerance for assertWithinExperimentalError
        self.tol = 0.01
        
    def test_dielectric_constant(self):
        '''        
        4.4 mol/kg NaCl = 46
        '''
        s1=pyEQL.Solution([['Na+','4.4 mol/kg'],['Cl-','4.4 mol/kg']])

        result=s1.get_dielectric_constant().magnitude
        expected = 46
                
        self.assertWithinExperimentalError(result,expected,self.tol)
        
    
    def test_dielectric_constant2(self):
        '''        
        2 mol/kg NaCl = 58
        '''
        s1=pyEQL.Solution([['Na+','2 mol/kg'],['Cl-','2 mol/kg']])

        result=s1.get_dielectric_constant().magnitude
        expected = 58
                
        self.assertWithinExperimentalError(result,expected,self.tol)
    
    
    def test_dielectric_constant3(self):
        '''        
        1 mol/kg NaCl = 66
        '''
        s1=pyEQL.Solution([['Na+','1 mol/kg'],['Cl-','1 mol/kg']])

        result=s1.get_dielectric_constant().magnitude
        expected = 66
                
        self.assertWithinExperimentalError(result,expected,self.tol)
    
    def test_dielectric_constant4(self):
        '''        
        1 mol/kg KBr = 67
        '''
        s1=pyEQL.Solution([['K+','1 mol/kg'],['Br-','1 mol/kg']])

        result=s1.get_dielectric_constant().magnitude
        expected = 67
                
        self.assertWithinExperimentalError(result,expected,self.tol)

    def test_dielectric_constant5(self):
        '''        
        3.4 mol/kg KBr = 51
        '''
        s1=pyEQL.Solution([['K+','3.4 mol/kg'],['Br-','3.4 mol/kg']])

        result=s1.get_dielectric_constant().magnitude
        expected = 51
                
        self.assertWithinExperimentalError(result,expected,self.tol)
    
    def test_dielectric_constant6(self):
        '''        
        5 mol/kg LiCl = 39
        '''
        s1=pyEQL.Solution([['Li+','5 mol/kg'],['Cl-','5 mol/kg']])

        result=s1.get_dielectric_constant().magnitude
        expected = 39
                
        self.assertWithinExperimentalError(result,expected,self.tol)
    
    def test_dielectric_constant7(self):
        '''        
        1 mol/kg LiCl = 64
        '''
        s1=pyEQL.Solution([['Li+','1 mol/kg'],['Cl-','1 mol/kg']])

        result=s1.get_dielectric_constant().magnitude
        expected = 64
                
        self.assertWithinExperimentalError(result,expected,self.tol)
    
    @unittest.expectedFailure
    def test_dielectric_constant8(self):
        '''        
        12 mol/kg LiCl = 24
        '''
        s1=pyEQL.Solution([['Li+','12 mol/kg'],['Cl-','12 mol/kg']])

        result=s1.get_dielectric_constant().magnitude
        expected = 24
                
        self.assertWithinExperimentalError(result,expected,self.tol)
        
    def test_dielectric_constant9(self):
        '''        
        6.5 mol/kg RbCl = 43
        '''
        s1=pyEQL.Solution([['Rb+','6.5 mol/kg'],['Cl-','6.5 mol/kg']])

        result=s1.get_dielectric_constant().magnitude
        expected = 43
                
        self.assertWithinExperimentalError(result,expected,self.tol)

    def test_dielectric_constant9(self):
        '''        
        2.1 mol/kg RbCl = 59
        '''
        s1=pyEQL.Solution([['Rb+','2.1 mol/kg'],['Cl-','2.1 mol/kg']])

        result=s1.get_dielectric_constant().magnitude
        expected = 59
                
        self.assertWithinExperimentalError(result,expected,self.tol)   

    def test_dielectric_constant9(self):
        '''        
        0.5 mol/kg RbCl = 73
        '''
        s1=pyEQL.Solution([['Rb+','0.5 mol/kg'],['Cl-','0.5 mol/kg']])

        result=s1.get_dielectric_constant().magnitude
        expected = 73
                
        self.assertWithinExperimentalError(result,expected,self.tol)
         
if __name__ == '__main__':
    unittest.main()
