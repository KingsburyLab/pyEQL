'''
pyEQL Volume and Concentration methods test
-------------------------------------------
This test runs through the different methods of adding or changing solute
amounts in pyEQL solutions. 

It has been tested primarily with mol/kg and mol/L units. In principle other units
e.g. g/L or g/kg should work, but use with caution

'''

from pyEQL import *
'''
Tests for add_solute() / solution init method
'''
print('\n Test 1')
# if solutes are added with substance / volume units, then the total volume
# of the solution should not change (should remain at 1 L)
soln1 = Solution([['Na+','4 mol/L'],['Cl-','4 mol/L']],volume='2 L')
# should read 4.0 mol/L
soln1.list_concentrations('mol/L')
soln1.list_concentrations('mol/kg')
# should be 2 L
print('Volume: '+str(soln1.get_volume()))
print('Water Mass: '+str(soln1.get_solvent_mass()))

print('\n Test 2')
# if solutes are added with substance / mass units, then the total volume
# of the solution will be recalculated based on the composition
soln1 = Solution([['Na+','4 mol/kg'],['Cl-','4 mol/kg']],volume='2 L')
# 
soln1.list_concentrations('mol/L')
soln1.list_concentrations('mol/kg')
# should be 2 L
print('Volume: '+str(soln1.get_volume()))
print('Water Mass: '+str(soln1.get_solvent_mass()))

'''
Tests for set_amount()
'''
print('\n Test 3')
# If the concentration of a solute is directly changed with a substance / volume
# unit, the volume should not change, but the water mass will be reduced
soln1 = Solution([['Na+','4 mol/L'],['Cl-','4 mol/L']],volume='2 L')
soln1.set_amount('Na+','5 mol/L')
soln1.set_amount('Cl-','5 mol/L')
# Na+ should be 5.0 mol/L
soln1.list_concentrations('mol/L')
soln1.list_concentrations('mol/kg')
# Volume should still be 2L
print('Volume: '+str(soln1.get_volume()))
print('Water Mass: '+str(soln1.get_solvent_mass()))

print('\n Test 4')
# If the concentration of a solute is directly set with a substance / mass
# unit, the volume should be updated, but the amount of water should not change
soln1 = Solution([['Na+','4 mol/kg'],['Cl-','4 mol/kg']],volume='2 L')
soln1.set_amount('Na+','5 mol/kg')
soln1.set_amount('Cl-','5 mol/kg')
# Na+ should be 5.0 mol/L
soln1.list_concentrations('mol/L')
soln1.list_concentrations('mol/kg')
# Volume should still be 2L
print('Volume: '+str(soln1.get_volume()))
print('Water Mass: '+str(soln1.get_solvent_mass()))

'''
Tests for add_amount()
'''
print('\n Test 5')
# If the concentration of a solute is directly changed with a substance / volume
# unit, the volume should not change, but the water mass will be reduced
soln1 = Solution([['Na+','4 mol/L'],['Cl-','4 mol/L']],volume='2 L')
soln1.add_amount('Na+','1 mol/L')
soln1.add_amount('Cl-','1 mol/L')
# Na+ should be 5.0 mol/L
soln1.list_concentrations('mol/L')
soln1.list_concentrations('mol/kg')
# Volume should still be 2L
print('Volume: '+str(soln1.get_volume()))
print('Water Mass: '+str(soln1.get_solvent_mass()))

print('\n Test 6')
# If the concentration of a solute is directly set with a substance / mass
# unit, the volume should be updated, but the amount of water should not change
soln1 = Solution([['Na+','4 mol/kg'],['Cl-','4 mol/kg']],volume='2 L')
soln1.add_amount('Na+','1 mol/kg')
soln1.add_amount('Cl-','1 mol/kg')
# Na+ should be 5.0 mol/L
soln1.list_concentrations('mol/L')
soln1.list_concentrations('mol/kg')
# Volume should still be 2L
print('Volume: '+str(soln1.get_volume()))
print('Water Mass: '+str(soln1.get_solvent_mass()))
