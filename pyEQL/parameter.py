# -*- coding: utf-8 -*-
"""
This module implements the Parameter() class, which is used to store
values, units, uncertainties, and reference data for various quantities
used throughout pyEQL.

:copyright: 2013-2018 by Ryan S. Kingsbury
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

## Units handling
# per the pint documentation, it's important that pint and its associated Unit
# Registry are only imported once.
from pint import UnitRegistry
# here we assign the identifier 'unit' to the UnitRegistry
unit = UnitRegistry()

#use this to enable legacy handling of offset units
# TODO fix this to handle offsets the way pint wants us to since 0.7
unit.autoconvert_offset_to_baseunit = True

# append custom unit definitions and contexts
import os
directory = os.path.dirname(__file__)
unit.load_definitions(directory+'/pint_custom_units.txt') 
# activate the "chemistry" context globally
unit.enable_contexts('chem')
# set the default string formatting for pint quantities
unit.default_format = 'P~'

def testfunc(val):
    list = []
    try:
        list.append(float(val) * unit(''))
        return list
    except ValueError:
        print('Value Error')
        return None


class Parameter:
    '''
    Class for storing and retrieving measured parameter values together with their
    units, context, and reference information.
    
    Some pyEQL functions search for specific parameter names, such as:
    diffusion_coefficient
    

    
    '''
    def __init__(self,name,magnitude,units='',**kwargs):
        '''
        Parameters
        ----------
        name : str
                    A short name (akin to a variable name) for the parameter
                    
        magnitude : float, int, str, tuple or list of floats, ints, or strs
                    The value of the parameter. In most cases this will only be a single numerical value.
                    However, in some cases it may be desirable to store a group of parameters (such as coefficients
                    for an equation) together in a tuple or list.
                    
                    Numeric values can be input as strings (or lists of strings) and they will be converted to 
                    floats. 
                    
                    Non-numeric values are permissible as well. When specifying non-numeric values, the units
                    argument must either be 'dimensionless' or left blank.
                    
                    Lists of non-numeric strings are not permitted.
                    
        units : str, optional
                    A string representing the units of measure for the parameter value
                    given in 'magnitude.' See the pint documentation for proper syntax. In general
                    common abbreviations or unit names can be used, but pythonic math syntax
                    must be followed. e.g. 'meters ** 2 / second' and 'm **2 /s' are valid but
                    'm2/s' is not.
                    
                    If a parameter has no units, leave blank or enter 'dimensionless'
                    
                    Note that if a parameter DOES have units but they are not specified, all
                    calculations involving this parameter will return incorrect units.
                    
        Optional Keyword Arguments
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
        temperature : str, optional
                    The temperature at which 'magnitude' was measured in degrees Celsius.
                    Specify the temperature as a string containing the magnitude and
                    a unit, e.g. '25 degC', '32 degF', '298 kelvin', and '500 degR'                    
        pressure : str, optional
                    The pressure at which 'magnitude' was measured in Pascals
                    Specify the pressure as a string containing the magnitude and a
                    unit. e.g. '101 kPa'.
                    Typical valid units are 'Pa', 'atm', or 'torr'.                   
        ionic_strength : str, optional
                    The ionic strength of the solution in which 'magnitude' was measured. Specify
                    the ionic strength as a string containing the magnitude and a unit. e.g. '2 mol/kg'
        description : str, optional
                    A string contiaining a longer name describing the paramter. For example
                    'Diffusion Coefficient' or 'Hydrated Ionic Radius'
        comments : str, optional
                    A string containing additional notes pertaining to the context,
                    conditions, or assumptions that may restrict the use of 'value'
                    
        Notes
        -----
        In general, parameter values are assumed to be entered in fundamental 
        SI units (m, kg, s, etc.). The 'units' field is required to call attention 
        to this fact and provide a levelof error-checking in calculations involving
        the parameter.
        
        Examples
        --------
        # TODO fix this example
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
        use_units = ''
        # turn numeric parameters into quantities with associated units
        # if units were specified as 'None', convert into something pint will understand
        if units == 'None' or units == 'none' or units == '' or units == 'dimensionless':
            use_units = 'dimensionless'
        else:
            use_units = units
        
        # see if the input value is a list or tuple. If so, create a list of
        # quantities (including units), and convert the list to a tuple
        if isinstance(magnitude,(tuple,list)):
            # check whether each element is a number
            temp_list=[]
            for item in magnitude:
                try:
                    temp_list.append(float(item) * unit(use_units))
                except ValueError:
                    print('Value Error on %s' % item)
                    # Throw an error if units are assigned to a non-numeric parameter
                    if not (use_units == 'dimensionless'):
                        logger.error('A non-numeric parameter cannot have units, but units of %s were specified' % units)
                    temp_list.append(item)
                    
            # convert the resulting list into a tuple
            self.value = tuple(temp_list)
            
        # if the input is a single item, try to convert it to a number. If that
        # doesn't work, it must be a str and will be passed on as-is
        else:
            try:
                self.value=float(magnitude) * unit(use_units)
            except ValueError:
                # Throw an error if units are assigned to a non-numeric parameter
                if not (use_units == 'dimensionless'):
                    logger.error('A non-numeric parameter cannot have units, but units of %s were specified' % units)
                
                self.value = magnitude
            
        # process optional keyword arguments - reference conditions
        self.base_temperature = 'Not specified'
        self.base_pressure = 'Not specified'
        self.base_ionic_strength = 'Not specified'
        
        if 'temperature' in kwargs and kwargs['temperature'] != '':
            self.temperature_set=True
            self.base_temperature = unit(kwargs['temperature'])
        if 'pressure' in kwargs and kwargs['pressure'] != '':
            self.pressure_set=True
            self.base_pressure = unit(kwargs['pressure'])
        if 'ionic_strength' in kwargs and kwargs['ionic_strength'] != '':
            self.ionic_strength_set=True
            self.base_ionic_strength = unit(kwargs['ionic_strength'])
            
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
       
    def get_name(self):
        '''
        Return the name of the parameter.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        str
            The name of the parameter
        '''
        return self.name
            
    def get_value(self,temperature=None,pressure=None,ionic_strength=None):
        '''
        Return the value of a parameter at the specified conditions.
        
        Parameters
        ----------        
        temperature : str, optional
                    The temperature at which 'magnitude' was measured in degrees Celsius.
                    Specify the temperature as a string containing the magnitude and
                    a unit, e.g. '25 degC', '32 degF', '298 kelvin', and '500 degR'                    
        pressure : str, optional
                    The pressure at which 'magnitude' was measured in Pascals
                    Specify the pressure as a string containing the magnitude and a
                    unit. e.g. '101 kPa'.
                    Typical valid units are 'Pa', 'atm', or 'torr'.                   
        ionic_strength : str, optional
                    The ionic strength of the solution in which 'magnitude' was measured. Specify
                    the ionic strength as a string containing the magnitude and a unit. e.g. '2 mol/kg' 
        
        Returns
        -------
        Quantity
            The value of the parameter at the specified conditions.
            
        '''
        # if the user does not specify conditions, return the value at base_temperature,
        # base_pressure, and/or base_ionic_strength
        if temperature is None: 
            temperature = self.base_temperature
            logger.info('Temperature not specified for '+str(self.name)+'. Returning value at '+str(temperature)+'.')
        else:
            temperature = unit(temperature)
        if pressure is None: 
            pressure = self.base_pressure
            logger.info('Pressure not specified for '+str(self.name)+'. Returning value at '+str(pressure)+'.')
        else:
            pressure = unit(pressure)
        if ionic_strength is None: 
            ionic_strength = self.base_ionic_strength
            logger.info('Ionic Strength not specified for '+str(self.name)+'. Returning value at '+str(ionic_strength)+'.')
        else:
            ionic_strength = unit(ionic_strength)
        
        # compare requested conditions with base conditions
        if temperature != self.base_temperature:        
            # TODO- implement temperature correction
            logger.warning('Requested temperature for '+str(self.name)+ \
            ' ('+str(temperature)+') differs from measurement conditions.' \
            +'Returning value at '+str(self.base_temperature))
            
        if pressure != self.base_pressure:        
            # TODO - implement pressure correction
            logger.warning('Requested pressure for '+str(self.name)+ \
            ' ('+str(pressure)+') differs from measurement conditions.' \
            +'Returning value at '+str(self.base_pressure))
        
        
        if ionic_strength != self.base_ionic_strength:        
            logger.warning('Requested ionic strength for '+str(self.name)+ \
            ' ('+str(ionic_strength)+') differs from measurement conditions.' \
            +'Returning value at '+str(self.base_ionic_strength))
        
        return self.value
        
    def get_magnitude(self,temperature=None,pressure=None,ionic_strength=None):
        '''
        Return the magnitude of a parameter at the specified conditions.
        
        Parameters
        ----------        
        temperature : str, optional
                    The temperature at which 'magnitude' was measured in degrees Celsius.
                    Specify the temperature as a string containing the magnitude and
                    a unit, e.g. '25 degC', '32 degF', '298 kelvin', and '500 degR'                    
        pressure : str, optional
                    The pressure at which 'magnitude' was measured in Pascals
                    Specify the pressure as a string containing the magnitude and a
                    unit. e.g. '101 kPa'.
                    Typical valid units are 'Pa', 'atm', or 'torr'.                   
        ionic_strength : str, optional
                    The ionic strength of the solution in which 'magnitude' was measured. Specify
                    the ionic strength as a string containing the magnitude and a unit. e.g. '2 mol/kg' 
        
        Returns
        -------
        Number
            The magnitude of the parameter at the specified conditions.
            
        '''
        return self.get_value(temperature,pressure,ionic_strength).magnitude
        
    def get_units(self):
        '''
        Return the units of a parameter
        '''
        return self.get_value().units
    
    def get_dimensions(self):
        '''
        Return the dimensions of the parameter.
        '''
        return self.get_value().dimensionality
    
    def __str__(self):
        '''
        Set the output of the print() statement for a parameter value
        '''
        return '\n'+'----------------------------------------------------------------------------'+'\n'+ \
        'Parameter '+str(self.name)+'\n\n'+str(self.description)+'\n' \
        'Value: '+str(self.get_value())+'\n'+ \
        'Conditions (T,P,Ionic Strength): '+str(self.base_temperature)+', '+str(self.base_pressure)+', ' \
        +str(self.base_ionic_strength)+'\n'+ \
        'Notes: '+str(self.comment)+'\n'+ \
        'Reference: '+str(self.reference)+'\n' + \
        '--------------------------------------------------------------------------------------'+'\n'

# TODO - turn doctest back on when the nosigint error is gone
## Tests
#if __name__ == "__main__":
 #   import doctest
  #  doctest.testmod()
