.. _chemistry:

Chemical Formulas
*****************

Representing Chemical Substances in pyEQL
=========================================

pyEQL interprets the chemical formula of a substance to calculate its molecular
weight and formal charge. The formula is also used as a key to search the 
database for parameters (e.g. diffusion coefficient) that are used in subsequent
calculations.

How to Enter Valid Chemical Formulas
------------------------------------

Generally speaking, type the chemical formula of your solute the "normal" way
and pyEQL should be able to inerpret it. Here are some examples:

- Sodium Chloride - NaCl
- Sodium Sulfate - Na(SO4)2
- Methanol - CH4OH or CH5O
- Magnesium Ion - Mg+2
- Chloride Ion - Cl-

Formula Rules:

1. Are composed of valid atomic symbols that start with capital letters
2. Contain no non-alphanumeric characters other than '(', ')', '+', or '-'
3. If a '+' or '-' is present, the formula must contain ONLY '+' or '-' 
   (e.g. 'Na+-' is invalid) and the formula must end with either a series of 
   charges (e.g. 'Fe+++') or a numeric charge (e.g. 'Fe+3')
4. Formula must contain matching numbers of '(' and ')'
5. Open parentheses must precede closed parentheses  

Alternate Formulas and Isomers
------------------------------

Many complex molecules can be written in multiple ways. pyEQL cares only about
the number and identity of the elements and the formal charge on the molecule,
so you can use any form you choose. The hill_order() method takes a formula
and reduces it to its simplest form, like so::

    >>> pyEQL.chemical_formula.hill_order('CH2(CH3)4COOH')
    'C6H15O2'

When searching the parameters database, pyEQL uses this method to reduce
both user-entered formulas AND keys in the database. So even if you created
a solution containing 'ClNa', pyEQL would still match it with parameters for
'NaCl'.

Currently pyEQL **does not distinguish between isomers**. 

API Documentation (chemical_formula.py)
=======================================

.. automodule:: pyEQL.chemical_formula
   :members:
