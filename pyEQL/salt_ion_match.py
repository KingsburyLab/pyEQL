'''
pyEQL salt matching library

This file contains functions that allow a pyEQL Solution object composed of
individual species (usually ions) to be mapped to a solution of one or more
salts. This mapping is necessary because some parameters (such as activity
coefficient data) can only be determined for salts (e.g. NaCl) and not individual
species (e.g. Na+)

:copyright: 2013-2016 by Ryan S. Kingsbury
:license: LGPL, see LICENSE for more details.

'''
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

import pyEQL.chemical_formula as chem

class Salt:
    '''
    Class to represent a salt.
    '''
    def __init__(self,cation,anion):
        self.cation=cation
        self.anion=anion
        
        '''
        Create a salt object based on its component ions
        
        Parameters:
        ----------
        cation, anion : str
                Chemical formula of the cation and anion, respectively
        
        Returns:
        -------
        Salt : An object representing the properties of the salt
        
        Examples:
        --------
        >>> Salt('Na+','Cl-').formula
        'NaCl'
        
        >>> Salt('Mg++','Cl-').formula
        'MgCl2'
        
        '''
    
        # get the charges on cation and anion
        self.z_cation = chem.get_formal_charge(cation)
        self.z_anion = chem.get_formal_charge(anion)
        
        # assign stoichiometric coefficients by finding a common multiple
        self.nu_cation = abs(self.z_anion)
        self.nu_anion = abs(self.z_cation)
        
        # if both coefficients are the same, set each to one
        if self.nu_cation == self.nu_anion:
            self.nu_cation = 1
            self.nu_anion = 1
            
        # start building the formula, cation first
        salt_formula=''
        if self.nu_cation > 1:
            # add parentheses if the cation is a polyatomic ion
            if len(chem.get_elements(cation)) > 1:
                salt_formula+='('
                salt_formula+= _trim_formal_charge(cation)
                salt_formula+=')'
            else:
                salt_formula+= _trim_formal_charge(cation)
            salt_formula+=str(self.nu_cation)
        else:
            salt_formula+= _trim_formal_charge(cation)
        
        if self.nu_anion > 1:
            # add parentheses if the anion is a polyatomic ion
            if len(chem.get_elements(anion)) > 1:
                salt_formula+='('
                salt_formula+= _trim_formal_charge(anion)
                salt_formula+=')'
            else:
                salt_formula+= _trim_formal_charge(anion)
            salt_formula+=str(self.nu_anion)
        else:
            salt_formula+= _trim_formal_charge(anion)
        
        self.formula = salt_formula

def _sort_components(Solution):
    '''
    Sort the components of a solution in descending order (by mole fraction).
    
    Parameters:
    ----------
    Solution : Solution object
    
    Returns:
    -------
    A list whose keys are the component names (formulas) and whose
    values are the component objects themselves
    
    
    '''
    formula_list = []
        
    # populate a list with component names
    for item in Solution.components:
        formula_list.append(item)
        
    # populate a dictionary with formula:concentration pairs
    mol_list={}
    for item in Solution.components.keys():
        mol_list.update({item:Solution.get_amount(item,'mol')})
    
    return sorted(formula_list,key=mol_list.__getitem__,reverse=True)
    
def identify_salt(Solution):
    '''
    Analyze the components of a solution and identify the salt that most closely
    approximates it.
    (e.g., if a solution contains 0.5 mol/kg of Na+ and Cl-, plus traces of H+
    and OH-, the matched salt is 0.5 mol/kg NaCl)
    
    Create a Salt object for this salt.
    
    Returns:
    -------
    A Salt object.
    '''
    # sort the components by moles
    sort_list = _sort_components(Solution)
   
    # default to returning water as the salt
    cation = 'H+'
    anion = 'OH-'    
        
    # return water if there are no solutes        
    if len(sort_list) < 3 and sort_list[0] == 'H2O':
        logger.info('Salt matching aborted because there are not enough solutes.')
        return Salt(cation,anion)
    
    # warn if something other than water is the predominant component    
    if sort_list[0] != 'H2O':
        logger.warning('H2O is not the most prominent component')
    
    # take the dominant cation and anion and assemble a salt from them
    for item in sort_list:
        if chem.get_formal_charge(item) > 0 and cation == 'H+':
            cation = item
        elif chem.get_formal_charge(item) < 0 and anion == 'OH-':
            anion = item
        else:
            pass
            
    # assemble the salt
    return Salt(cation,anion)
    
def _trim_formal_charge(formula):
    '''
    remove the formal charge from a chemical formula
    
    Examples:
    --------
    >>> _trim_formal_charge('Fe+++')
    'Fe'
    >>> _trim_formal_charge('SO4-2')
    'SO4'
    >>> _trim_formal_charge('Na+')
    'Na'
    
    '''
    charge = chem.get_formal_charge(formula)
    output = ''
    if charge > 0:
        output = formula.split('+')[0]
    elif charge < 0:
        output = formula.split('-')[0]

    return output

# TODO - turn doctest back on when the nosigint error is gone        
#if __name__ == "__main__":
 #   import doctest
  #  doctest.testmod()
