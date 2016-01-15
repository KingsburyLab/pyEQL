'''
pyEQL methods for chemical equilibrium calculations (e.g. acid/base, reactions,
redox, complexation, etc.)

NOTE: these methods are not currently used but are here for the future.

:copyright: 2013-2016 by Ryan S. Kingsbury
:license: LGPL, see LICENSE for more details.

'''

## Dependencies
# import libraries for scientific functions
import math

# the pint unit registry
from pyEQL import unit

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

def adjust_temp_pitzer(c1,c2,c3,c4,c5,temp,temp_ref=unit('298.15 K')):
    '''
    Calculate a parameter for th e Pitzer model based on temperature-dependent
    coefficients c1,c2,c3,c4,and c5.
    
    Parameters
    ----------
    c1, c2, c3, c4, c5: float
                Temperature-dependent coefficients for the pitzer parameter of 
                interest.
    temp: Quantity
                The temperature at which the Pitzer parameter is to be calculated
    temp_ref: Quantity, optional
                The reference temperature on which the parameters are based.
                298.15 K if omitted.
    
    As described in the PHREEQC documentation
    
    
    '''
    pitzer_param = c1 + c2 * (1/temp + 1/temp_ref) + c2 * math.log(temp/temp_ref) \
    + c3 * (temp - temp_ref) + c4 * (temp ** 2 - temp_ref ** 2) + c5 * (temp ** -2 - temp_ref ** -2)
    
    return pitzer_param

def adjust_temp_vanthoff(equilibrium_constant,enthalpy,temperature,reference_temperature = 25*unit('degC')):
    '''(float,float,number, optional number) -> float
    
    Adjust a reaction equilibrium constant from one temperature to another.
    
    Parameters
    ----------
    equilibrium_constant : float
                           The reaction equilibrium constant for the reaction
    enthalpy : Quantity
               The enthalpy change (delta H) for the reaction in kJ/mol. Assumed
               independent of temperature (see Notes).
    temperature : Quantity
                  the desired reaction temperature in degrees Celsius
    reference_temperature : Quantity, optional
                      the temperature at which equilibrium_constant is valid. (25 degrees C if omitted).
   
    Returns
    -------
    float
        adjusted reaction equilibrium constant
    
    Notes
    -----
    This function implements the Van't Hoff equation to adjust measured 
    equilibrium constants to other temperatures. 
    
    .. math::
        ln(K2 / K1) = {\delta H \over R} ( {1 \over T_1} - {1 \over T_2} )
    
    This implementation assumes that the enthalpy is independent of temperature
    over the range of interest.[1]
    
    .. [1] Stumm, Werner and Morgan, James J. Aquatic Chemistry, 3rd ed, pp 53. 
        Wiley Interscience, 1996.
    
    Examples
    --------
    >>> adjust_temp_vanthoff(0.15,-197.6*unit('kJ/mol'),42*unit('degC'),25*unit('degC')) #doctest: +ELLIPSIS
    0.00203566...
    
    If the 'ref_temperature' parameter is omitted, a default of 25 C is used.
    
    >>> adjust_temp_vanthoff(0.15,-197.6*unit('kJ/mol'),42*unit('degC')) #doctest: +ELLIPSIS
    0.00203566...
    
    '''
    output = equilibrium_constant * math.exp( enthalpy / unit.R * ( 1 / reference_temperature.to('K') - 1 / temperature.to('K')))
    
    logger.info('Adjusted equilibrium constant K=%s from %s to %s degrees Celsius with Delta H = %s. Adjusted K = %s % equilibrium_constant,reference_temperature,temperature,enthalpy,output')
    
    logger.warning("Van't Hoff equation assumes enthalpy is independent of temperature over the range of interest")
    return output

def adjust_temp_arrhenius(rate_constant,activation_energy,temperature,reference_temperature = 25*unit('degC')):
    '''(float,float,number, optional number) -> float
    
    Adjust a reaction equilibrium constant from one temperature to another.
    
    Parameters
    ----------
    rate_constant : Quantity
                The parameter value (usually a rate constant) being adjusted
    activation_energy : Quantity
               The activation energy of the process, in kJ/mol
    temperature : Quantity
                  the desired reaction temperature.
    reference_temperature : Quantity, optional
                      the temperature at which equilibrium_constant is valid
                      Defaults to 25 degrees C if omitted.
   
    Returns
    -------
    Quantity
        The adjusted reaction equilibrium constant

    See Also
    --------
    kelvin
    
    Notes
    -----
    This function implements the Arrhenius equation to adjust measured rate
    constants to other temperatures. [1]
    
    .. math::
        ln(K2 / K1) = {E_a \over R} ( {1 \over T_1} - {1 \over T_2} )
    
    .. [1] http://chemwiki.ucdavis.edu/Physical_Chemistry/Kinetics/Reaction_Rates/Temperature_Dependence_of_Reaction_Rates/Arrhenius_Equation
    TODO - add better reference
    
    Examples
    --------
    >>> adjust_temp_arrhenius(7,900*unit('kJ/mol'),37*unit('degC'),97*unit('degC')) #doctest: +ELLIPSIS
    1.8867225...e-24
    
    '''
    output = rate_constant * math.exp( activation_energy / unit.R * ( 1 / reference_temperature.to('K') - 1 / temperature.to('K')))
    
    logger.info('Adjusted parameter %s from %s to %s degrees Celsius with Activation Energy = %s kJ/mol. Adjusted value = %s % rate_constant,reference_temperature,temperature,activation_energy,output')
    
    return output

def alpha(n,pH,pKa_list):
    '''(int,number,list of numbers)
    Returns the acid-base distribution coefficient (alpha) of an acid in the 
    n-deprotonated form at a given pH.
    
    Parameters
    ----------
    n : int
        The number of protons that have been lost by the desired form of the
        acid. Also the subscript on the alpha value. E.g. for bicarbonate
        (HCO3-), n=1 because 1 proton has been lost from the fully-protonated
        carbonic acid (H2CO3) form.
    pH : float or int
         The pH of the solution.
    pKa_list : list of floats or ints
               The pKa values (negative log of equilibrium constants) for the acid
               of interest. There must be a minimum of n pKa values in the list.
    
    Returns
    -------
    float
        The fraction of total acid present in the specified form.
    
    Notes
    -----
    The acid-base distribution coefficient is calculated as follows:[1]
        
    .. math::
        \alpha_n = {term_n \over [H+]^n + k_{a1}[H+]^n-1 + k_{a1}k_{a2}[H+]^n-2 ... k_{a1}k_{a2}...k_{an} }
    
    Where :math: '\term_n' refers to the nth term in the denominator, starting from 0
    
    .. [1] Stumm, Werner and Morgan, James J. Aquatic Chemistry, 3rd ed, 
        pp 127-130. Wiley Interscience, 1996.
    
    Examples
    --------
    >>> alpha(1,8,[4.7]) #doctest: +ELLIPSIS
    0.999...
    
    The sum of all alpha values should equal 1
    
    >>> alpha(0,8,[6.35,10.33]) #doctest: +ELLIPSIS
    0.021...
    >>> alpha(1,8,[6.35,10.33]) #doctest: +ELLIPSIS
    0.979...
    >>> alpha(2,8,[6.35,10.33]) #doctest: +ELLIPSIS
    2.043...e-09
    
    If pH is equal to one of the pKa values the function should return 0.5.
    
    >>> alpha(1,6.35,[6.35,10.33])
    0.5
    
#     The function will return an error if the number ofpKa's is less than n.
#     
#     >>> alpha(2,8,[])
#     ERROR: insufficient number of pKa values given
#     0.5   

    '''
    #generate an error if no pKa values are specified
    if len(pKa_list) == 0:
        logger.error('No pKa values given. Cannot calculate distribution coeffiicent.')
        return None
    
    #generate an error if n > number of pKa values
    if len(pKa_list) < n:
        logger.error('Insufficient number of pKa values given. Cannot calculate distribution coeffiicent.')
        return None
        
    #convert pH to hydrogen ion concentration
    Hplus = 10 ** -pH
    
    #determine how many protons the acid has
    num_protons = len(pKa_list)
    
    #build a list of terms where the term subscript corresponds to the list index
    terms_list = []
    k_term = 1
    
    #the 'item' index counts from 0 to the number of protons, inclusive
    for item in range(0,num_protons+1):
        #multiply the preceding k values together
        for i in range(len(pKa_list[:item])):
            k_term *= 10 ** -pKa_list[i]
        
        #add the term to the list
        terms_list.append(k_term * Hplus ** (num_protons - item))
    
    #build the expression
    numerator = terms_list[n]
    denominator = 0
    for item in terms_list:
        denominator += item
        
    #return the desired distribution factor
    alpha = numerator / denominator
    logger.info('Calculated %s-deprotonated acid distribution coefficient of %s for pKa=%s at pH %s % n,alpha,pKa_list,pH')
    return alpha
