chemical_formula.py test suite
==============================

This file contains tests for the chemical formula interpreter module of pyEQL.

>>> import chemical_formula as cf


tests for the _check_formula() function
---------------------------------------
This function is the core of the library and enforces the rules for 
chemical formula formatting.

>>> cf._check_formula('Fe2(SO4)3')
['Fe', '2', '(', 'S', 'O', '4', ')', '3']

>>> cf._check_formula('C7H16')
['C', '7', 'H', '16']

>>> cf._check_formula('(NH3)2SO4')
['(', 'N', 'H', '3', ')', '2', 'S', 'O', '4']

>>> cf._check_formula('MgCl2')
['Mg', 'Cl', '2']

>>> cf._check_formula('C100H202')
['C', '100', 'H', '202']

>>> cf._check_formula('Fe+++')
['Fe', '+++']

>>> cf._check_formula('V+4')
['V', '+', '4']


tests for is_valid_formula()
---------------------------

A formula must start with a letter or an open parenthesis
>>> cf.is_valid_formula('(NH3)2')
True

>>> cf.is_valid_formula('3CO3-')
False

A formula cannot contain any non-alphanumeric characters beside '(', ')', '+', and '-'
>>> cf.is_valid_formula('Na^+')
False

A formula cannot contain both '+' and '-'
>>> cf.is_valid_formula('Na+-+')
False

An ionic formula must end with either a number, a '+', or a '-'
>>> cf.is_valid_formula('HCO3-')
True

>>> cf.is_valid_formula('Fe++')
True

>>> cf.is_valid_formula('Mg+2')
True

>>> cf.is_valid_formula('V+5+')
False

Formulas must contain only valid atomic symbols that start with capital letters
>>> cf.is_valid_formula('NaOH')
True

>>> cf.is_valid_formula('naOH')
False

>>> cf.is_valid_formula('HzCl')
False

A formula with parentheses must have the same number of '(' and ')'
>>> cf.is_valid_formula('(NH3)2(NO3)2')
True

>>> cf.is_valid_formula('Mg)(OH)2')
False

A formula cannot end with an open parenthesis
>>> cf.is_valid_formula('Na+(')
False

An open parenthesis must always be followed by an atomic symbol
>>> cf.is_valid_formula('(3)Na+')
False

>>> cf.is_valid_formula('()')
False

A closed parenthesis may be followed by a number or used to designate a group
>>> cf.is_valid_formula('CH3C(O)CH3')
True

>>> cf.is_valid_formula('Mg(OH)2')
True

An open parenthesis must always precede the nearest closed parenthesis
>>> cf.is_valid_formula('CH3(CH)2(CH)2')
True

>>> cf.is_valid_formula(')Na+(')
False

tests for _consolidate_formula()
----------------------------
This function is important for calculating molecular weights correctly

>>> cf._consolidate_formula('Fe2(SO4)4')
['Fe', 2, 'S', 4, 'O', 16]

>>> cf._consolidate_formula('(NH4)3PO4')
['N', 3, 'H', 12, 'P', 1, 'O', 4]

>>> cf._consolidate_formula('CH3(CH2)6CH3')
['C', 8, 'H', 18]
