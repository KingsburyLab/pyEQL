# -*- coding: utf-8 -*-
"""
This module implements the Parameter() class, which is used to store
values, units, uncertainties, and reference data for various quantities
used throughout pyEQL.

Created on Wed Jun 11 09:27:01 2014

@author: ryan
"""

## Logging System
''' Create a logging system using Python's built-in module. 
Add the null handler to avoid errors in case the calling application doesn't configure any handlers.

NOTE: make sure to set the disable_existing_loggers option in the log configuration
options of the calling application in order to avoid disabling the pyEQL module's log
 
The default logging levels are mapped to pyEQL events as follows:
 
DEBUG       -   detailed messages about function execution including methods used, data sources,
                temperature adjustments, etc.
INFO        -   Messages indicating calculation steps, function calls, etc.
WARNING     -   assumptions or limitations of module output
ERROR       -   Module could not complete a task due to invalid input or other problem
CRITICAL    -   not used

'''
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# per the pint documentation, it's important that pint and its associated Unit
# Registry are only imported once.
from pint import UnitRegistry
# here we assign the identifier 'unit' to the UnitRegistry
unit = UnitRegistry()
# append custom unit definitions and contexts
unit.load_definitions('pint_custom_units.txt') 


class Parameter:
    '''
    Class for storing and retrieving measured parameter values together with their
    units, context, and reference information.

    
    '''
    def __init__(self,name,magnitude,units,**kwargs):
        '''
        Parameters:
        ----------
        name : str
                    A short name (akin to a variable name) for the parameter
        magnitude : float, int, or tuple of floats or ints
                    The numerical value of the parameter. In most cases this will only be a single value.
                    However, in some cases it may be desirable to store a group of parameters (such as coefficients
                    for an equation) together in a tuple.
        units : str
                    A string representing the units of measure for the parameter value
                    given in 'magnitude.' See the pint documentation for proper syntax. In general
                    common abbreviations or unit names can be used, but pythonic math syntax
                    must be followed. e.g. 'meters ** 2 / second' and 'm **2 /s' are valid but
                    'm2/s' is not.
                    
                    If a parameter has no untis, leave the second element blank or enter 'dimensionless'
                    
                    Note that if a parameter DOES have units but they are not specified, all
                    calculations involving this parameter will return incorrect units.
                    
        Optional Keyword Arguments:
        --------------------------
        reference : str, optional
                    A string containing reference information documenting the source of
                    the parameter.
        uncertainty : tuple of floats or ints, optional
                    The uncertainty of the parameter value as a percentage (0-100). If only one value is supplied in the tuple,
                    the same uncertainty will be used for positive and negative. If two values are supplied, the first
                    is the positive uncertainty and the second is the negative.
        uncertainty_type: str, optional
                    A string specifying different ways of measuring uncertainty.
        temperature : tuple, optional
                    The temperature at which 'magnitude' was measured in degrees Celsius.
                    The first element of the tuple is a float or int representing the
                    temperature. The second element is a string representing the unit.
                    Valid temperature units are 'degC', 'degF', 'kelvin', and 'degR'
                    Defaults to 25 degrees C if omitted.
        pressure : tuple, optional
                    The pressure at which 'magnitude' was measured in Pascals
                    The first element of the tuple is a float or int representing the 
                    pressure. The second element is a string representing the unit.
                    Typical valid units are 'Pa', 'atm', or 'torr'.
                    Defaults to 1 atm (101325) Pa if omitted
       ionic_strength : float or int, optional
                    The ionic strength of the solution in which 'magnitude' was measured. 
                    Defaults to 0 (infinite dilution) if omitted.
        description : str, optional
                    A string contiaining a longer name describing the paramter. For example
                    'Diffusion Coefficient' or 'Hydrated Ionic Radius'
        comments : str, optional
                    A string containing additional notes pertaining to the context,
                    conditions, or assumptions that may restrict the use of 'value'
                    
        Notes:
        -----
        In general, parameter values are assumed to be entered in fundamental 
        SI units (m, kg, s, etc.). The 'units' field is required to call attention 
        to this fact and provide a levelof error-checking in calculations involving
        the parameter.
        
        Examples:
        --------
        >>> sodium_diffusion = Parameter('diffusion coefficient',(1.334e-9,),'m2/s','CRC Handbook of Chemistry and Physics, 92nd Ed., pp. 5-77 to 5-79',(),25,101325,0)
        >>> print(sodium_diffusion)
        Parameter diffusion coefficient
        <BLANKLINE>
        -------------------------------------------
        Value: (1.334e-09,) m2/s at 25 degrees C.
        Notes: 
        Reference: CRC Handbook of Chemistry and Physics, 92nd Ed., pp. 5-77 to 5-79
        <BLANKLINE>
        '''
        self.name = name
        self.value = magnitude * unit(units)
        
        # process optional keyword arguments - reference conditions
        self.base_temperature = 'Not specified'
        self.base_pressure = 'Not specified'
        self.base_ionic_strength = 'Not specified'
        
        if 'temperature' in kwargs:
            self.temperature_set=True
            self.base_temperature = kwargs["temperature"][0]*unit(kwargs["temperature"][1])
        if 'pressure' in kwargs:
            self.pressure_set=True
            self.base_pressure = kwargs["pressure"][0]*unit(kwargs["pressure"][1])
        if 'ionic_strength' in kwargs:
            self.ionic_strength_set=True
            self.base_ionic_strength = kwargs["ionic_strength"][0]*unit(kwargs["ionic_strength"][1])
            
        # process optional keyword arguments - descriptive and reference
        self.reference = 'Not specified'
        self.description = ''
        self.comment = ''         
            
        if 'reference' in kwargs:
            self.reference=kwargs["reference"]
        if 'description' in kwargs:
            self.description = kwargs['description']
        if 'comment' in kwargs:
            self.comment = kwargs['comment']
            
        # process optional keyword arguments - uncertanty
            # TODO - uncertainty functions / parameters
        
        # process optional keyword arguments  - temperature adjustment
        # TODO - temperature adjustment functions / parameters
     
       
    def get_value(self,temperature=None,pressure=None,ionic_strength=None):
        # TODO- implement temperature correction and make conditions mandatory
        '''return a temperature-adjusted paramter value and log any qualifying
        assumptions
        '''
        # use the base values for temperature, pressure, ionic strength if none
        # are specified
        if temperature is None: temperature = self.base_temperature
        if pressure is None: pressure = self.base_pressure
        if ionic_strength is None: ionic_strength = self.base_ionic_strength
    
        return self.value
        
    def get_magnitude(self,temperature):
        '''return the temperature-adjusted magnitude of the parameter
        '''
        return self.get_value(temperature).magnitude
        
    def get_units(self,temperature):
        '''return the temperature-adjusted magnitude of the parameter
        '''
        return self.get_value(temperature).units
    
    def get_dimensions(self,temperature):
        '''return the temperature-adjusted magnitude of the parameter
        '''
        return self.get_value(temperature).dimensionality
    
    def __str__(self):
        '''
        Set the output of the print() statement for a parameter value
        '''
        return 'Parameter '+str(self.name)+'\n'+str(self.description)+'\n' \
        +'-------------------------------------------'+'\n'+ \
        'Value: '+str(self.get_value())+'\n'+ \
        'Conditions (Temperature, Pressure, Ionic Strength): '+str(self.base_temperature)+', '+str(self.base_pressure)+', '+ \
        str(self.base_ionic_strength)+'\n'+ \
        'Notes: '+str(self.comment)+'\n'+ \
        'Reference: '+str(self.reference)+'\n'
