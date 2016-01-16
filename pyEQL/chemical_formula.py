# -*- coding: utf-8 -*-
"""
This module contains classes, functions, and methods to facilitate the
input, output, and parsing of chemical formulas for pyEQL.

The correct case must be used when specifying elements.

:copyright: 2013-2016 by Ryan S. Kingsbury
:license: LGPL, see LICENSE for more details.

"""

# logging system
import logging
logger = logging.getLogger(__name__)

# add a filter to emit only unique log messages to the handler
import pyEQL.logging_system
unique = pyEQL.logging_system.Unique()
logger.addFilter(unique)

# add a handler for console output, since pyEQL is meant to be used interactively
ch = logging.StreamHandler()

# create formatter for the log
formatter = logging.Formatter('(%(name)s) - %(levelname)s - %(message)s')

# add formatter to the handler
ch.setFormatter(formatter)
logger.addHandler(ch)

## Formula validation and processing functions. These internal routines
## parse chemical formulas into a format that can be easily processed
## by user-facing functions.
def _invalid_formula(reason):
    raise ValueError('Invalid chemical formula specified - %s' % reason)
    return None
    
def _check_formula(formula):
    '''
    Parse a chemical formula into a list that separates atomic symbols, 
    numbers, and parentheses, and check the formula for compliance with 
    formatting rules. 
    
    Similar to Python's default list() function for strings.
    
    Parameters
    ----------
    formula: str
        String representing a molecular formula. e.g. 'H2O' or 'FeOH+'
        Valid molecular formulas must meet the following criteria:
        
        #. Are composed of valid atomic symbols that start with capital letters
        #. Contain no non-alphanumeric characters other than '(', ')',
           '+', or '-'
        #. If a '+' or '-' is present, the formula must contain ONLY '+' or
           '-' (e.g. 'Na+-' is invalid) and the formula must end with either
           a series of charges (e.g. 'Fe+++') or a numeric charge (e.g. 'Fe+3')
        #. Formula must contain matching numbers of '(' and ')'
        #. Open parentheses must precede closed parentheses
        
        
    Examples
    --------
    >>> _check_formula('Fe2(SO4)3')
    ['Fe', '2', '(', 'S', 'O', '4', ')', '3']
    >>> _check_formula('C10H12')
    ['C', '10', 'H', '12']
    '''
    
    # check that formula starts with a letter or open parenthesis
    if not (formula[0].isalpha() or formula[0] == '('):
        _invalid_formula('formula must begin with an element or open parenthesis')
    
    # check for mismatched charges 
    if '+' in formula and '-' in formula:
        _invalid_formula('ionic formulas cannot contain mismatched charge symbols')
    
    # check that ionic formulas end with either a number of a charge symbol
    if '+' in formula:
        if formula.count('+') == 1 and not (formula[-1] == '+' or formula[-1].isnumeric()):
            _invalid_formula('ionic formulas must end with one or more charge symbols or a single charge symbol and a number')    
        elif formula.count('+') > 1:
            start = formula.find('+')
            for char in formula[start:]:
                if char != '+':
                    _invalid_formula('ionic formulas must end with one or more charge symbols or a single charge symbol and a number')    
    elif '-' in formula:
        if formula.count('-') == 1 and not (formula[-1] == '-' or formula[-1].isnumeric()):
            _invalid_formula('ionic formulas must end with one or more charge symbols or a single charge symbol and a number')    
        elif formula.count('-') > 1:           
            start = formula.find('-')
            for char in formula[start:]:
                if char != '-':
                    _invalid_formula('ionic formulas must end with one or more charge symbols or a single charge symbol and a number')    

    # check for equal parentheses
    if formula.count('(') != formula.count(')'):
        _invalid_formula('parentheses mismatch')
    
    # make sure open parenthesis doesn't end the formula
    if formula.endswith('('):
        _invalid_formula('formula cannot end with open parenthesis')

    # split the formula string into a list of characters
    input_list = list(formula)    

    
    for i in range(len(input_list)):
        try:
            # check for invalid characters
            parentheses = ['(',')']            
            charge_symbols = ['+','-']
            
            if not (input_list[i].isalnum() or input_list[i] in \
            parentheses or input_list[i] in charge_symbols):
                _invalid_formula('contains invalid character')
            
            elif input_list[i] == '(':
                # check that open parentheses are followed by an atomic symbol
                if not input_list[i+1].isalpha():
                    _invalid_formula('parentheses must contain elements')
                # make sure that the open parenthesis precedes the nearest closed
                # parenthesis
                try:
                    if not i < input_list.index(')',i):
                        _invalid_formula('open parenthesis must precede closed parenthesis')
                # add exception for ValueError, in case there is no closed parenthesis
                # after index i
                except ValueError:    
                    _invalid_formula('open parenthesis must precede closed parenthesis')
             
            # removed this rule to allow for organic structural formulas
            # check that closed parenthesis are followed by a number
#            elif input_list[i] == ')':
#                try:
#                    if not input_list[i+1].isnumeric():
#                        _invalid_formula('parnetheses must be followed by numbers')
#                except IndexError:
#                    _invalid_formula('parentheses must be followed by numbers')
            
            # concatenate any uppercase letters with up to two subsequent lowercase letters
            elif input_list[i].isupper():
                try:                
                    if input_list[i+1].islower():
                        try:
                            if input_list[i+2].islower():
                                char = input_list.pop(i+2)
                                input_list[i+1]+= char
                        except IndexError:
                            pass
                        char = input_list.pop(i+1)
                        input_list[i] += char
                except IndexError:
                    pass
                    
            # concatenate any adjacent numbers together
            elif input_list[i].isnumeric():
                try:
                    j = i+1
                    while input_list[j].isnumeric():
                        char = input_list.pop(j)
                        input_list[i] += char
                except IndexError:
                    pass
            
            # concatenate adjacent + or -
            elif input_list[i] == '+' or input_list[i] == '-':
                try:
                    if input_list[i+1] == '+' or input_list[i+1] == '-':
                        j = i+1
                        while input_list[j] == '+' or input_list[j] == '-':
                            char = input_list.pop(j)
                            input_list[i] += char
                except IndexError:
                    pass
                
            else:
                pass
            
        except IndexError:
            pass
    
    # check that all elements are valid
    for item in input_list:
        if item.isalpha():
            if not is_valid_element(item):
                _invalid_formula('invalid element symbol')
    
    return input_list   
   
    
def _remove_parentheses(formula):
    '''
    Remove parentheses from a formula and distribute the associated numbers
    as appropriate. 
    
    NOTE: does not support nested parentheses as these violate
    the formatting rules for chemical formulas
    
    >>> _remove_parentheses('(Fe2)(SO4)3')
    ['Fe', '2', 'S', '3', 'O', '12']

    
    See Also
    --------
    _check_formula()
    '''
 
    # perform validity check and return a list of the chemical formula's components
    input_list = _check_formula(formula)
    output_list = []
    
    # remove all parentheses from the formula and distribute numbers accordingly   
    i = 0    
    while i < len(input_list):
        if input_list[i] == '(':
            # locate the beginning and end indices of the parenthetical group
            start = i
            stop = input_list.index(')',i)

            # locate the number after the group, if any
            if input_list[stop+1].isnumeric():
                num = int(input_list[stop+1])
                # skip past the parenthetical group once the loop is done
                i = stop+2
            else:
                num = 1
                # skip past the parenthetical group once the loop is done
                i = stop+1
                
            # loop through the elements / numbers contained in the group
            for j in range(start+1,stop):            
                if input_list[j].isalpha() and input_list[j+1].isnumeric():
                    output_list.append(input_list[j])
                    output_list.append(str(int(input_list[j+1])*num))
                elif input_list[j].isalpha():
                    output_list.append(input_list[j])
                    if num > 1:
                        output_list.append(str(num))

        else:
            output_list.append(input_list[i])
            # advance to the next list element
            i = i+1
        
    return output_list

def _consolidate_formula(formula):
    '''
    Consolidate a formula into its simplest form, containing only one
    instance of each element and no parentheses
    
    Examples
    --------
    >>> _consolidate_formula('CH3(CH2)6CH3')
    ['C', 8, 'H', 18]
    >>> _consolidate_formula('(Fe)2(SO4)4')
    ['Fe', 2, 'S', 4, 'O', 16]
    >>> _consolidate_formula('Fe(OH)2+')
    ['Fe', 1, 'O', 2, 'H', 2, '+1']
    
    '''
    # perform validity check and return a list of the chemical formula's components
    input_list = _remove_parentheses(formula)
    output_list = []
    
    for i in range (0,len(input_list)):
        # is the item an element, a number, or a charge?
        if input_list[i].isalpha():
            # is it followed by a number?
            try:
                if input_list[i+1].isnumeric():
                    quantity = input_list[i+1]
                else:
                    quantity = 1
            except IndexError:
                quantity = 1
                
            # have we seen it before?
            if input_list[i] in output_list:
                # yes, so find it in the output list
                index = output_list.index(input_list[i])
                # add the quantity
                output_list[index+1] += int(quantity)
            else:
                # no, so add it and its quantity
                output_list.append(input_list[i])
                output_list.append(int(quantity))
          
    # include any charge symbols
    charge = get_formal_charge(formula)
    if charge > 0:
        output_list.append('+' + str(charge))
    elif charge < 0:
        output_list.append(str(charge))

    return output_list

## Truth Functions
def is_valid_element(formula):
    '''
    Check whether a string is a valid atomic symbol
    
    Parameters
    ----------
    :formula: str
            String representing an atomic symbol. First letter must be 
            uppercase, second letter must be lowercase.
    
    Returns
    -------
    bool
            True if the string is a valid atomic symbol. False otherwise.     
    
    Examples
    --------
    >>> is_valid_element('Cu')
    True
    >>> is_valid_element('Na+')
    False
    '''
    if formula in atomic_numbers:
        return True
    else:
        _invalid_formula('invalid element symbol')
        return False

def is_valid_formula(formula):
    '''
    Check that a molecular formula is formatted correctly
    
    Parameters
    ----------
    formula: str
        String representing a molecular formula. e.g. 'H2O' or 'FeOH+'
        Valid molecular formulas must meet the following criteria:
        
        #. Are composed of valid atomic symbols that start with capital letters
        #. Contain no non-alphanumeric characters other than '(', ')',
           '+', or '-'
        #. If a '+' or '-' is present, the formula must contain ONLY '+' or
           '-' (e.g. 'Na+-' is invalid) and the formula must end with either
           a series of charges (e.g. 'Fe+++') or a numeric charge (e.g. 'Fe+3')
        #. Formula must contain matching numbers of '(' and ')'
        #. Open parentheses must precede closed parentheses
    
    Returns
    -------
    bool
            True if the formula is valid. False otherwise.
            
    Examples
    --------
    >>> is_valid_formula('Fe2(SO4)3')
    True
    >>> is_valid_formula('2Na+')
    False
    >>> is_valid_formula('HCO3-')
    True
    >>> is_valid_formula('Na+-')
    False
    >>> is_valid_formula('C10h12')
    False 
    '''
    
    try:
        _check_formula(formula)
        return True
    except:
        return False
        
def contains(formula,element):
    '''
    Check whether a formula contains a given element.
    
    Parameters
    ----------
    formula: str
        String representing a molecular formula. e.g. 'H2O' or 'FeOH+'
        Valid molecular formulas must meet the following criteria:
        
        #. Are composed of valid atomic symbols that start with capital letters
        #. Contain no non-alphanumeric characters other than '(', ')',
           '+', or '-'
        #. If a '+' or '-' is present, the formula must contain ONLY '+' or
           '-' (e.g. 'Na+-' is invalid) and the formula must end with either
           a series of charges (e.g. 'Fe+++') or a numeric charge (e.g. 'Fe+3')
        #. Formula must contain matching numbers of '(' and ')'
        #. Open parentheses must precede closed parentheses
    element: str
        String representing the element to check for. Must be a valid element 
        name.
    
    Returns
    -------
    bool
            True if the formula contains the element. False otherwise.
            
    Examples
    --------
    >>> contains('Fe2(SO4)3','Fe')
    True
    >>> contains('NaCOOH','S')
    False
    '''
    if is_valid_element(element):
        if element in get_elements(formula):
            return True
        else:
            return False  
    
## Information Retrieval Functions    
def get_element_numbers(formula):
    '''
    Return the atomic numbers of the elements in a chemical formula
    
    Parameters
    ----------
    formula: str
            String representing a chemical formula
    
    Examples
    --------
    >>> get_element_numbers('FeSO4')
    [26, 16, 8]
    

    '''
    # perform validity check and return a list of the chemical formula's components
    input_list = get_elements(formula)
    output_list = []
    
    for item in input_list:
        output_list.append(atomic_numbers[item][0])        
    
    return output_list
    
def get_element_names(formula):
    '''
    Return the names of the elements in a chemical formula
    
    Parameters
    ----------
    formula: str
            String representing a chemical formula
    
    Examples
    --------
    >>> get_element_names('FeSO4')
    ['Iron', 'Sulfur', 'Oxygen']
    

    '''
    # perform validity check and return a list of the chemical formula's components
    input_list = get_elements(formula)
    output_list = []
    
    for item in input_list:
        output_list.append(atomic_numbers[item][1])        
    
    return output_list
    
def hill_order(formula):
    '''
    Return a string representing the simplest form of 'formula'
    in the Hill order (Carbon, Hydrgen, then other elements
    in alphabetical order). If no Carbon is present, then
    all elements are listed in alphabetical order.
    
    NOTE: this function does NOT (yet) honor exceptions to the Hill Order
    for acids, hydroxides, oxides, and ionic compounds. It follows the
    rule above no matter what.
    
    Examples
    --------
    >>> hill_order('CH2(CH3)4COOH')
    'C6H15O2'

    >>> hill_order('NaCl')
    'ClNa'
    
    >>> hill_order('NaHCO2') == hill_order('HCOONa')
    True
    
    >>> hill_order('Fe+2') == hill_order('Fe+3')
    False

    '''
    ### TODO - add exceptions for oxides (end in O2), acids (start with H),  
    ## ions (cation first), and hydroxides (ends in OH)
    temp_list = _consolidate_formula(formula)
    hill = ''
    
    # start the formula with C and H if Carbon is present
    if 'C' in temp_list:
        for item in ['C','H']:
            if item in temp_list:
                index=temp_list.index(item)
                hill += item
                # copy the number only if greater than 1        
                if temp_list[index+1] > 1:
                    hill += str(temp_list.pop(index+1))
                elif temp_list[index+1] == 1:
                    temp_list.pop(index+1)
                temp_list.remove(item)
        
    # convert any remaining list entries into tuples of (element,number)
    # so they can be sorted
    tuple_list=[]
    for item in temp_list:
        index=temp_list.index(item)
        try:
            if item.isalpha():
                tuple_list.append((item,temp_list[index+1]))
        except AttributeError:
            continue
    
    # put the remaining elements in alphabetical order
    tuple_list.sort()
    
    # add them to the formula
    for item in tuple_list:
        if item[1] ==1:
            hill += str(item[0])
        else:
            hill += str(item[0]) + str(item[1])
    
    # append the formal charge to the end of the formula, if not zero
    charge = get_formal_charge(formula)
    if charge != 0:
        hill += str(charge)
    
    return hill  
        
def get_elements(formula):
    '''
    Return a list of strings representing the elements in a 
    molecular formula, with no duplicates.
    
    Examples
    --------
    >>> get_elements('FeSO4')
    ['Fe', 'S', 'O']
    >>> get_elements('CH3(CH2)4(CO)3')
    ['C', 'H', 'O']
    
    See Also
    --------
    _check_formula()
    '''
    # perform validity check and return a parsed list of the chemical formula
    input_list = _consolidate_formula(formula)
    output_list = []

    for item in input_list:
        if item in atomic_numbers:
            output_list.append(item)

    return output_list

def get_formal_charge(formula):
    '''
    Return the formal charge on a molecule based on its formula
    
    Examples
    --------
    >>> get_formal_charge('Na+')
    1
    >>> get_formal_charge('PO4-3')
    -3
    >>> get_formal_charge('Fe+++')
    3
    
    See Also
    --------
    _check_formula()
    
    '''
    # perform validity check and return a parsed list of the chemical formula
    input_list = _check_formula(formula)
    
    if '+' in input_list:
        index= input_list.index('+')
        try:
            formal_charge = 1 * int(input_list[index+1])
        except:
            formal_charge = 1
    elif '-' in input_list:
        index= input_list.index('-')
        try:
            formal_charge = -1 * int(input_list[index+1])
        except:
            formal_charge = -1
    elif '+' in input_list[-1] :
        formal_charge = int(1 * input_list[-1].count('+'))
    elif '-' in input_list[-1] :
        formal_charge = int(-1 * input_list[-1].count('-'))
    else:
        formal_charge = 0
        
    return formal_charge

def get_element_mole_ratio(formula,element):
    '''
    compute the  moles of a specific element per mole of formula
    
    Parameters
    ----------
    formula: str
        String representing a molecular formula. e.g. 'H2O' or 'FeOH+'
        Valid molecular formulas must meet the following criteria:
        
        #. Are composed of valid atomic symbols that start with capital letters
        #. Contain no non-alphanumeric characters other than '(', ')',
           '+', or '-'
        #. If a '+' or '-' is present, the formula must contain ONLY '+' or
           '-' (e.g. 'Na+-' is invalid) and the formula must end with either
           a series of charges (e.g. 'Fe+++') or a numeric charge (e.g. 'Fe+3')
        #. Formula must contain matching numbers of '(' and ')'
        #. Open parentheses must precede closed parentheses
    element: str
        String representing the element to check for. Must be a valid element 
        name.
    
    Returns
    -------
    number
            The number of moles of element per mole of formula, mol/mol.
    
    >>> get_element_mole_ratio('NaCl','Na')
    1
    >>> get_element_mole_ratio('H2O','H')
    2
    >>> get_element_mole_ratio('H2O','Br')
    0
    >>> get_element_mole_ratio('CH3CH2CH3','C')
    3
    
    See Also
    --------
    contains
    consolidate_formula
    get_element_weight
    get_element_weight_fraction
    
    '''
    # perform validity check and return a parsed list of the chemical formula
    if contains(formula,element):
        input_list = _consolidate_formula(formula)
        index = input_list.index(element)
        moles = input_list[index+1]              
    # return 0 weight if the element isn't present in the formula
    else:
        moles = 0

    return moles
    
def get_element_weight(formula,element):
    '''
    compute the  weight of a specific element in a formula
    
    Parameters
    ----------
    formula: str
        String representing a molecular formula. e.g. 'H2O' or 'FeOH+'
        Valid molecular formulas must meet the following criteria:
        
        #. Are composed of valid atomic symbols that start with capital letters
        #. Contain no non-alphanumeric characters other than '(', ')',
           '+', or '-'
        #. If a '+' or '-' is present, the formula must contain ONLY '+' or
           '-' (e.g. 'Na+-' is invalid) and the formula must end with either
           a series of charges (e.g. 'Fe+++') or a numeric charge (e.g. 'Fe+3')
        #. Formula must contain matching numbers of '(' and ')'
        #. Open parentheses must precede closed parentheses
    element: str
        String representing the element to check for. Must be a valid element 
        name.
    
    Returns
    -------
    number
            The weight of the specified element within the formula, g/mol.
    
    >>> get_element_weight('NaCl','Na')
    22.98977
    >>> get_element_weight('H2O','H')
    2.01588
    >>> get_element_weight('H2O','Br')
    0.0
    >>> get_element_weight('CH3CH2CH3','C')
    36.0321
    
    See Also
    --------
    contains()
    _consolidate_formula()
    elements
    get_element_mole_ratio
    
    '''
    # find the number of moles of element per mole of formula
    moles = get_element_mole_ratio(formula,element)
    
    if moles != 0:
        # import elements.py - used to retreive various molecular data
        from pyEQL.elements import ELEMENTS    
            
        # look up the molecular weight for the element
        mass = ELEMENTS[element].mass
        
        wt = mass * moles
    else:
        wt = 0.0

    return wt

def get_element_weight_fraction(formula,element):
    '''
    compute the  weight fraction of a specific element in a formula
    
    Parameters
    ----------
    formula: str
        String representing a molecular formula. e.g. 'H2O' or 'FeOH+'
        Valid molecular formulas must meet the following criteria:
        
        #. Are composed of valid atomic symbols that start with capital letters
        #. Contain no non-alphanumeric characters other than '(', ')',
           '+', or '-'
        #. If a '+' or '-' is present, the formula must contain ONLY '+' or
           '-' (e.g. 'Na+-' is invalid) and the formula must end with either
           a series of charges (e.g. 'Fe+++') or a numeric charge (e.g. 'Fe+3')
        #. Formula must contain matching numbers of '(' and ')'
        #. Open parentheses must precede closed parentheses
    element: str
        String representing the element to check for. Must be a valid element 
        name.
    
    Returns
    -------
    number
            The weight fraction of the specified element within the formula.
    
    >>> get_element_weight_fraction('NaCl','Na')
    0.39337...
    >>> get_element_weight_fraction('H2O','H')
    0.111898...
    >>> get_element_weight_fraction('H2O','Br')
    0.0
    >>> get_element_weight_fraction('CH3CH2CH3','C')
    0.8171355...
    
    See Also
    --------
    get_element_weight
    contains
    _consolidate_formula
    elements
    
    '''
    # calculate the element weight in the formula
    wt = get_element_weight(formula,element)
    
    # calculate the fraction
    frac = wt  / get_molecular_weight(formula)
    
    return frac

def get_molecular_weight(formula):
    '''
    compute the molecular weight of a formula
    
    >>> get_molecular_weight('Na+')
    22.98977
    >>> get_molecular_weight('H2O')
    18.01528
    >>> get_molecular_weight('CH3CH2CH3')
    44.09562
    
    See Also
    --------
    _consolidate_formula()
    elements
    
    '''
    # import elements.py - used to retreive various molecular data
    from pyEQL.elements import ELEMENTS    
    
    # perform validity check and return a parsed list of the chemical formula
    input_list = _consolidate_formula(formula)
    mw = 0
    
    for item in input_list:
        try:
            if item.isalpha():
                index = input_list.index(item)
                quantity = input_list[index+1]              
                # look up the molecular weight for the element
                mass = ELEMENTS[item].mass
                mw += mass * quantity
                
        # if the list item is a number or a charge, move on
        except AttributeError:
            pass
    
    return mw
         
## Output functions

def print_latex(formula):
    '''
    Print a LaTeX - formatted version of the formula
    
    Examples
    ---------
    >>> print_latex('Fe2SO4')
    Fe_2SO_4
    >>> print_latex('CH3CH2CH3')
    CH_3CH_2CH_3
    >>> print_latex('Fe2(OH)2+2')
    Fe_2(OH)_2^+^2
    
    '''
    output = ''
    for i in range(len(formula)):
        # insert LaTeX subscripts (underscore) before all numbers, unless they start
        # the formula or follow a charge symbol
        if formula[i].isnumeric() and i > 0:
            # insert superscript before a charge
            if formula[i-1] == '+' or formula[i-1] == '-':
                output += '^'
            else:
                output += '_'
        
        # insert superscripts before charges
        if formula[i] == '+' or formula[i] == '-':
            output += '^'
        
        # pass the rest of the formula to output
        output += formula[i]
    
    print(output)
            

# Dictionary to correlate atomic symbols, names, and numbers
atomic_numbers={
    'H':(1,'Hydrogen'),
    'He':(2,'Helium'),
    'Li':(3,'Lithium'),
    'Be':(4,'Beryllium'),
    'B':(5,'Boron'),
    'C':(6,'Carbon'),
    'N':(7,'Nitrogen'),
    'O':(8,'Oxygen'),
    'F':(9,'Fluorine'),
    'Ne':(10,'Neon'),
    'Na':(11,'Sodium'),
    'Mg':(12,'Magnesium'),
    'Al':(13,'Aluminum'),
    'Si':(14,'Silicon'),
    'P':(15,'Phosphorus'),
    'S':(16,'Sulfur'),
    'Cl':(17,'Chlorine'),
    'Ar':(18,'Argon'),
    'K':(19,'Potassium'),
    'Ca':(20,'Calcium'),
    'Sc':(21,'Scandium'),
    'Ti':(22,'Titanium'),
    'V':(23,'Vanadium'),
    'Cr':(24,'Chromium'),
    'Mn':(25,'Manganese'),
    'Fe':(26,'Iron'),
    'Co':(27,'Cobalt'),
    'Ni':(28,'Nickel'),
    'Cu':(29,'Copper'),
    'Zn':(30,'Zinc'),
    'Ga':(31,'Gallium'),
    'Ge':(32,'Germanium'),
    'As':(33,'Arsenic'),
    'Se':(34,'Selenium'),
    'Br':(35,'Bromine'),
    'Kr':(36,'Krypton'),
    'Rb':(37,'Rubidium'),
    'Sr':(38,'Strontium'),
    'Y':(39,'Yttrium'),
    'Zr':(40,'Zirconium'),
    'Nb':(41,'Niobium'),
    'Mo':(42,'Molybdenum'),
    'Tc':(43,'Technetium'),
    'Ru':(44,'Ruthenium'),
    'Rh':(45,'Rhodium'),
    'Pd':(46,'Palladium'),
    'Ag':(47,'Silver'),
    'Cd':(48,'Cadmium'),
    'In':(49,'Indium'),
    'Sn':(50,'Tin'),
    'Sb':(51,'Antimony'),
    'Te':(52,'Tellurium'),
    'I':(53,'Iodine'),
    'Xe':(54,'Xenon'),
    'Cs':(55,'Cesium'),
    'Ba':(56,'Barium'),
    'La':(57,'Lanthanum'),
    'Ce':(58,'Cerium'),
    'Pr':(59,'Praseodymium'),
    'Nd':(60,'Neodymium'),
    'Pm':(61,'Promethium'),
    'Sm':(62,'Samarium'),
    'Eu':(63,'Europium'),
    'Gd':(64,'Gadolinium'),
    'Tb':(65,'Terbium'),
    'Dy':(66,'Dysprosium'),
    'Ho':(67,'Holmium'),
    'Er':(68,'Erbium'),
    'Tm':(69,'Thulium'),
    'Yb':(70,'Ytterbium'),
    'Lu':(71,'Lutetium'),
    'Hf':(72,'Hafnium'),
    'Ta':(73,'Tantalum'),
    'W':(74,'Tungsten'),
    'Re':(75,'Rhenium'),
    'Os':(76,'Osmium'),
    'Ir':(77,'Iridium'),
    'Pt':(78,'Platinum'),
    'Au':(79,'Gold'),
    'Hg':(80,'Mercury'),
    'Tl':(81,'Thallium'),
    'Pb':(82,'Lead'),
    'Bi':(83,'Bismuth'),
    'Po':(84,'Polonium'),
    'At':(85,'Astatine'),
    'Rn':(86,'Radon'),
    'Fr':(87,'Francium'),
    'Ra':(88,'Radium'),
    'Ac':(89,'Actinium'),
    'Th':(90,'Thorium'),
    'Pa':(91,'Protactinium'),
    'U':(92,'Uranium'),
    'Np':(93,'Neptunium'),
    'Pu':(94,'Plutonium'),
    'Am':(95,'Americium'),
    'Cm':(96,'Curium'),
    'Bk':(97,'Berkelium'),
    'Cf':(98,'Californium'),
    'Es':(99,'Einsteinium'),
    'Fm':(100,'Fermium'),
    'Md':(101,'Mendelevium'),
    'No':(102,'Nobelium'),
    'Lr':(103,'Lawrencium'),
    'Rf':(104,'Rutherfordium'),
    'Db':(105,'Dubnium'),
    'Sg':(106,'Seaborgium'),
    'Bh':(107,'Bohrium'),
    'Hs':(108,'Hassium'),
    'Mt':(109,'Meitnerium'),
    'Ds':(110,'Darmstadtium'),
    'Rg':(111,'Roentgenium'),
    'Cn':(112,'Copernicium'),
    'Uut':(113,'Ununtrium'),
    'Fl':(114,'Flerovium'),
    'Uup':(115,'Ununpentium'),
    'Lv':(116,'Livermorium'),
    'Uus':(117,'Ununseptium'),
    'Uuo':(118,'Ununoctium')
    }

# TODO - turn doctest back on when the nosigint error is gone    
#if __name__ == "__main__":
 #   import doctest
  #  doctest.testfile('tests/test_chemical_formula.rst')
   # doctest.testmod()
