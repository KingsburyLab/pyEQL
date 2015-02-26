'''
pyEQL solute class

This file contains functions and methods for managing properties of 
individual solutes

:copyright: 2013-2015 by Ryan S. Kingsbury
:license: LGPL, see LICENSE for more details.

'''

# the pint unit registry
from pyEQL import unit

# logging system
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# add a filter to emit only unique log messages to the handler
import pyEQL.logging_system
unique = pyEQL.logging_system.Unique()
logger.addFilter(unique)

# import the parameters database
from pyEQL import paramsDB as db

class Solute:
    '''represent each chemical species as an object containing its formal charge, 
    transport numbers, concentration, activity, etc. 
    
    '''
    
    def __init__(self,formula,amount,volume,solvent_mass,parameters={}):
        '''
        Parameters
        ----------
        formula : str
                    Chemical formula for the solute. 
                    Charged species must contain a + or - and (for polyvalent solutes) a number representing the net charge (e.g. 'SO4-2').
        amount : str
                    The amount of substance in the specified unit system. The string should contain both a quantity and
                    a pint-compatible representation of a unit. e.g. '5 mol/kg' or '0.1 g/L'
        volume : pint Quantity
                    The volume of the solution
        solvent_mass : pint Quantity
                    The mass of solvent in the parent solution.
        parameters : dictionary, optional
                    Dictionary of custom parameters, such as diffusion coefficients, transport numbers, etc. Specify parameters as key:value pairs separated by commas within curly braces, e.g. {diffusion_coeff:5e-10,transport_number:0.8}. The 'key' is the name that will be used to access the parameter, the value is its value.
        '''             
        # import the chemical formula interpreter module
        import pyEQL.chemical_formula as chem
        
        # check that 'formula' is a valid chemical formula
        if not chem.is_valid_formula:
            logger.error('Invalid chemical formula specified.')
            return None
        else:
            self.formula = formula

            # set molecular weight 
            self.mw = chem.get_molecular_weight(formula) * unit('g/mol')
            
            # set formal charge
            self.charge = chem.get_formal_charge(formula)
            
            # translate the 'amount' string into a pint Quantity
            quantity = unit(amount)
            
            self.moles = quantity.to('moles','chem',mw=self.mw,volume=volume,solvent_mass=solvent_mass)                

            # trigger the function that checks whether parameters already exist for this species, and if not,
            # searches the database files and creates them
            db.search_parameters(self.formula)

    def get_parameter(self,parameter,temperature=None,pressure=None,ionic_strength=None):
        '''
        Return the value of the parameter named 'parameter'
        
        '''
        # Search for this solute in the database of parameters
        
        # TODO - handle the case where multiple parameters with same name exist. Need function
        # to compare specified conditions and choose the most appropriate parameter
        param = db.get_parameter(self.formula,parameter)
        
        return param.get_value(temperature,pressure,ionic_strength)
                
    def add_parameter(self,name,magnitude,units='',**kwargs):
        '''manually insert a parameter into the parameters database for a solute
        
        See pyEQL.parameters documentation for a description of the arguments
        
        '''
        import pyEQL.parameter as pm
        newparam = pm.Parameter(name,magnitude,units,**kwargs)
        db.add_parameter(self.get_name(),newparam)
        
    def get_name(self):
        return self.formula
        
    def get_formal_charge(self):
        return self.charge
        
    def get_molecular_weight(self):
        return self.mw
    
    def get_moles(self):
        return self.moles
    
    def add_moles(self,amount,volume,solvent_mass):
        '''Increase or decrease the amount of a substance present in the solution
        
        Parameters
        ----------
        amount: str quantity
                Amount of substance to add. Must be in mass or substance units.
                Negative values indicate subtraction of material.
        
        '''
        quantity = unit(amount)
        self.moles += quantity.to('moles','chem',mw=self.mw,volume=volume,solvent_mass=solvent_mass)
    
    def set_moles(self,amount,volume,solvent_mass):
        quantity = unit(amount)
        self.moles = quantity.to('moles','chem',mw=self.mw,volume=volume,solvent_mass=solvent_mass)  
  
    def get_molar_conductivity(self,temperature=25*unit('degC')):
        # TODO - requires diffusion coefficient which may not be present
        '''(float,int,number) -> float
        Calculate the molar (equivalent) conductivity at infinte dilution for a species
        
        Parameters
        ----------
        temperature : Quantity, optional
                    The temperature of the parent solution. Defaults to 25 degC if omitted.
            
        Returns
        -------
        float
                The molar or equivalent conductivity of the species at infinte dilution.
        
        Notes
        -----
        Molar conductivity is calculated from the Nernst-Einstein relation:[1]_
            
        .. math:: \\kappa_i = {z_i^2 D_i F^2 \\over RT}
        
        Note that the diffusion coefficient is strongly variable with temperature.
        
        References
        ----------
        
        .. [1] Smedley, Stuart. The Interpretation of Ionic Conductivity in Liquids, pp 1-9. Plenum Press, 1980.
        
        TODO Examples
#             --------
#             
        '''
        diffusion_coefficient = self.get_parameter('diffusion_coefficient')
        
        molar_cond = diffusion_coefficient * (unit.e * unit.N_A) ** 2 * self.get_formal_charge() ** 2 / (unit.R * temperature)
        
        return molar_cond.to('mS / cm / (mol/L)')
            
    def get_mobility(self,temperature=25*unit('degC')):
        '''Return the ionic mobility of a species
        
        Parameters
        ----------
        temperature : Quantity, optional
                    The temperature of the parent solution. Defaults to 25 degC if omitted.
        
        Returns
        -------
        float : the ionic mobility
        
        
        Notes
        -----
        This function uses the Einstein relation to convert a diffusion coefficient
        into an ionic mobility [1]_
        
        .. math:: \mu_i = {F |z_i| D_i \over RT}
        
        References
        ----------
        .. [1] Smedley, Stuart I. The Interpretation of Ionic Conductivity in Liquids. Plenum Press, 1980.
        
        '''
        mobility = unit.N_A * unit.e * abs(self.get_formal_charge()) * self.get_parameter('diffusion_coefficient') / (unit.R * temperature)
        
        logger.info('Computed ionic mobility as %s from D = %s at T=%S' % mobility,self.get_diffusion_coefficient(),temperature)
        
        return mobility
        
    #set output of the print() statement
    def __str__(self):
        return 'Species ' + str(self.get_name()) + ' MW=' + str(self.get_molecular_weight()) +' Formal Charge='+str(self.get_formal_charge()) + ' Amount= ' + str(self.get_moles())