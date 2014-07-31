'''pyEQL - Python Library for Aquatic Chemistry Calculations
Developed by: Ryan Kingsbury
RyanSKingsbury@alumni.unc.edu

TODO - LICENSE INFO
'''

## Dependencies
# import libraries for scientific functions
import math
import numpy as np
# used for fundamental constants
from scipy import constants as spc

# FUTURE - to be used for plotting
#import matplotlib as mpl
#import matplotlib.pyplot as plt
# FUTURE - to be used for automatic calculation of molecular weights
#from elements import ELEMENTS as pte

# internal pyEQL imports
# the parameter handling module
import parameter as pm
# the pint unit registry
from parameter import unit
# functions to manage importing paramters from database files and making them accessible to pyEQL
import database
from database import parameters_database as db

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

## Automated Tests
''' This section will only run when the pyEQL module is executed on its own (i.e.
not called from another program). It creates a couple of standard test solutions
that are used for automated tests of the functions.
'''

## Fundamental Constants
# Values for fundamental constants are provided by the Scipy package
# float values (commented out) are from Wikipedia

# Avogadro's number, #/mole
#CONST_Na = 6.02214129e23
CONST_Na = spc.N_A

# Universal gas constant, Joule/mole-Kelvin
#CONST_R = 8.3144621 
CONST_R = spc.R

# Fundamental charge, coulombs
#CONST_e = 1.602176565e-19
CONST_e = spc.e

# Permittivity of free space, Farad/meter
#CONST_Eo = 8.854187817e-12
CONST_Eo = spc.epsilon_0

# Faraday constant, coulombs/mole (derived)
#9.64853399e4 
CONST_F = CONST_e * CONST_Na

# Boltzmann constant, joule/Kelvin (derived)
#1.3806488e-23
#CONST_kb = CONST_R / CONST_Na
CONST_kb = spc.k

# Acceleration due to gravity, meters/second^2
#CONST_g = 9.80665
CONST_g = spc.g

## Temperature Functions

def kelvin(temp_celsius):
    '''(number) -> float
    Convert a temperature from degrees Celsius to Kelvin
    
    Parameters:
    ----------
    temp_celsius : float or int
                   temperature in degrees Celsius
    
    Returns:
    -------
    float
        temperature in Kelvin
    
    Examples:
    --------
    >>> kelvin(25)
    298.15
    '''
    output = temp_celsius + 273.15
    logger.debug('Converted %s degrees Celsius to %s kelvin' % (temp_celsius, output))
    return output
    
    
def celsius(temp_kelvin):
    '''(number) -> float
    Convert a temperature in Kelvin to degrees Celsius
    
    Parameters:
    ----------
    temp_kelvin : float or int
                  temperature in Kelvin
    
    Returns:
    -------
    float
        temperature in degrees Celsius
    
    Examples:
    --------
    >>> celsius(298) #doctest: +ELLIPSIS
    24.85...
    '''
    output = temp_kelvin - 273.15
    logger.debug('Converted %s kelvin to %s degrees Celsius' % (temp_kelvin,output))
    return output
    
    
def adjust_temp_vanthoff(equilibrium_constant,enthalpy,temperature,reference_temperature = 25):
    '''(float,float,number, optional number) -> float
    
    Adjust a reaction equilibrium constant from one temperature to another.
    
    Parameters:
    ----------
    equilibrium_constant : float
                           The reaction equilibrium constant for the reaction
    enthalpy : float
               The enthalpy change (delta H) for the reaction in kJ/mol. Assumed
               independent of temperature (see Notes).
    temperature : float or int
                  the desired reaction temperature in degrees Celsius
    reference_temperature : float or int, optional
                      the temperature at which equilibrium_constant is valid
                      in degrees Celsius. (25 degrees C if omitted).
   
    Returns:
    -------
    float
        adjusted reaction equilibrium constant

    See Also:
    --------
    kelvin
    
    Notes:
    -----
    This function implements the Van't Hoff equation to adjust measured 
    equilibrium constants to other temperatures. 
    
    .. math::
        ln(K2 / K1) = {\delta H \over R} ( {1 \over T_1} - {1 \over T_2} )
    
    This implementation assumes that the enthalpy is independent of temperature
    over the range of interest.[1]
    
    .. [1] Stumm, Werner and Morgan, James J. Aquatic Chemistry, 3rd ed, pp 53. 
        Wiley Interscience, 1996.
    
    Examples:
    --------
    >>> adjust_temp_vanthoff(0.15,-197.6,42,25) #doctest: +ELLIPSIS
    0.00203566...
    
    If the 'ref_temperature' parameter is omitted, a default of 25 C is used.
    
    >>> adjust_temp_vanthoff(0.15,-197.6,42) #doctest: +ELLIPSIS
    0.00203566...
    
    '''
    output = equilibrium_constant * math.exp( enthalpy * 1000 / CONST_R * ( 1 / kelvin(reference_temperature) - 1 / kelvin(temperature)))
    
    logger.info('Adjusted equilibrium constant K=%s from %s to %s degrees Celsius with Delta H = %s. Adjusted K = %s % equilibrium_constant,reference_temperature,temperature,enthalpy,output')
    
    logger.warning("Van't Hoff equation assumes enthalpy is independent of temperature over the range of interest")
    return output

def adjust_temp_arrhenius(rate_constant,activation_energy,temperature,reference_temperature = 25):
    '''(float,float,number, optional number) -> float
    
    Adjust a reaction equilibrium constant from one temperature to another.
    
    Parameters:
    ----------
    rate_constant : float
                The parameter value (usually a rate constant) being adjusted
    activation_energy : float
               The activation energy of the process, in kJ/mol
    temperature : float or int
                  the desired reaction temperature in degrees Celsius
    reference_temperature : float or int, optional
                      the temperature at which equilibrium_constant is valid
                      in degrees Celsius. (25 degrees C if omitted).
   
    Returns:
    -------
    float
        adjusted reaction equilibrium constant

    See Also:
    --------
    kelvin
    
    Notes:
    -----
    This function implements the Arrhenius equation to adjust measured rate
    constants to other temperatures. [1]
    
    .. math::
        ln(K2 / K1) = {E_a \over R} ( {1 \over T_1} - {1 \over T_2} )
    
    .. [1] http://chemwiki.ucdavis.edu/Physical_Chemistry/Kinetics/Reaction_Rates/Temperature_Dependence_of_Reaction_Rates/Arrhenius_Equation
    TODO - add better reference
    
    Examples:
    --------
    >>> adjust_temp_arrhenius(7,900,37,97) #doctest: +ELLIPSIS
    1.8867225...e-24
    
    '''
    output = rate_constant * math.exp( activation_energy * 1000 / CONST_R * ( 1 / kelvin(reference_temperature) - 1 / kelvin(temperature)))
    
    logger.info('Adjusted parameter %s from %s to %s degrees Celsius with Activation Energy = %s kJ/mol. Adjusted value = %s % rate_constant,reference_temperature,temperature,activation_energy,output')
    
    return output
    
def adjust_temp_diffusion(diffusion_coefficient,temperature,ref_temperature=25):
    '''Adjust a diffusion coefficient for a different temperature
    The diffusion coefficients are corrected to temperature T (Kelvin) of the solution with:
        (Dw)T = (Dw)298 × (T / 298) × (η298 / ηT), where η is the viscosity of water. 
    TODO - FIND A LEGIT REFERENCE FOR THAT EQUATION
    
    .. [1] http://www.hydrochemistry.eu/exmpls/sc.html
    '''
    ## TODO - check this function (does it need dynamic or kinematic viscosity or does it matter?)
    return diffusion_coefficient * temperature / kelvin(ref_temperature) * water_viscosity_dynamic(ref_temperature)/water_viscosity_dynamic(temperature)

## Properties of Water

def water_density(temperature=25,pressure=101325):
    # TODO add pressure??
    # TODO more up to date equation??
    '''(number) -> float
    
    Return the density of water in kg/m3 at the specified temperature and pressure.
    
    Parameters:
    ----------
    temperature : float or int, optional
                  The temperature in Celsius. Defaults to 25 degrees if not specified.
    pressure    : float or int, optional
                  The ambient pressure of the solution in Pascals (N/m2). 
                  Defaults to atmospheric pressure (101325 Pa) if not specified.
    
    Returns:
    -------
    float
            The density of water in kg/m3.
    
    Notes:
    -----
    Based on the following empirical equation reported in [1]
    
    $$ \rho_W = 999.65 + 0.20438 T - 6.1744e-2 T ^ 1.5
    
    Where T is the temperature in Celsius.
    
    
    ..[1] Sohnel, O and Novotny, P. //Densities of Aqueous Solutions of Inorganic Substances.// Elsevier Science, Amsterdam, 1985.
    
    Examples:
    --------
    >>> water_density(25) #doctest: +ELLIPSIS
    997.04...
    
    '''
    density = 999.65 + 0.20438 * temperature - 6.1744e-2 * temperature ** 1.5
    logger.info('Computed density of water as %s kg/m3 at T= %s degrees C and P = %s Pa' % (density,temperature,pressure))
    logger.debug('Computed density of water using empirical relation in Sohnel and Novotny, "Densities of Aqueous Solutions of Inorganic Substances," 1985' )
    return density
    
    
def water_specific_weight(temperature,pressure=101325):
    '''(number) -> float
    
    Return the specific weight of water in N/m3 at the specified temperature and pressure.
    
    Parameters:
    ----------
    temperature : float or int, optional
                  The temperature in Celsius. Defaults to 25 degrees if not specified.
    pressure    : float or int, optional
                  The ambient pressure of the solution in Pascals (N/m2). 
                  Defaults to atmospheric pressure (101325 Pa) if not specified.
                  
    Returns:
    -------
    float
            The specific weight of water in N/m3.  
            
    See Also:
    --------
    water_density
    
    '''
    spweight = water_density(temperature,pressure) * CONST_g
    logger.info('Computed specific weight of water as %s N/m3 at T=%S degrees C and P = %s Pa' % (spweight,temperature,pressure))
    return spweight


def water_viscosity_dynamic(temperature=25,pressure=101325):
    '''
    Return the dynamic (absolute) viscosity of water in N-s/m2 = Pa-s = kg/m-s
    at the specified temperature.
    
    Parameters:
    ----------
    temperature : float or int, optional
                  The temperature in Celsius. Defaults to 25 degrees if not specified.
    pressure    : float or int, optional
                  The ambient pressure of the solution in Pascals (N/m2). 
                  Defaults to atmospheric pressure (101325 Pa) if not specified.
    
    Returns:
    -------
    float 
                The dynamic (absolute) viscosity of water in N-s/m2 = Pa-s = kg/m-s
                  
    Notes:
    -----
    Implements the international equation for viscosity of water as specified by NIST[1]
    
    Valid for 273 < temperature < 1073 K and 0 < pressure < 100,000,000 Pa
    
    .. [1] Sengers, J.V. "Representative Equations for the Viscosity of Water Substance." 
        J. Phys. Chem. Ref. Data 13(1), 1984.http://www.nist.gov/data/PDFfiles/jpcrd243.pdf
    
    Examples:
    --------
    >>> water_viscosity_dynamic(20) #doctest: +ELLIPSIS
    8.934...e-0.7
    >>> water_viscosity_dynamic(100,25000000) #doctest: +ELLIPSIS
    2.979...e-0.7
    >>> water_viscosity_dynamic(300,100000000) #doctest: +ELLIPSIS
    1.329...e-0.7
    #TODO - check these again after I implement pressure-dependent density function
    
    '''
    # generate warnings if temp or pressure are outside valid range of equation
    if kelvin(temperature) < 273 or kelvin(temperature)>1073:
        logger.error('Specified temperature (%s K) exceeds valid range of NIST equation for viscosity of water. Cannot extrapolate. % kelvin(temperature)')
        return None
        
    if pressure < 0 or pressure > 100000000:
        logger.error('Specified pressure (%s Pa) exceeds valid range of NIST equation for viscosity of water. Cannot extrapolate. % pressure')
        return None
    
    # calculate dimensionless temperature and pressure
    T_star = 647.27 #K
    P_star = 22115000 #Pa
    rho_star = 317.763 #kg/m3
    
    T_bar = kelvin(temperature) / T_star
    P_bar = pressure / P_star
    rho_bar = water_density(temperature,pressure) / rho_star
    
    # calculate the first function, mu_o
    mu_star = 1e-6 #Pa-s
    a = [0.0181583,0.0177624,0.0105287,-0.0036477]
    sum_o = 0
    mu_temp = 0
    for index in range(len(a)):
        sum_o += a[index] * T_bar ** -index
    
    mu_o = mu_star * math.sqrt(T_bar) / sum_o
    
    # calculate the second fucntion, mu_1
    b=[[0.501938,0.235622,-0.274637,0.145831,-0.0270448],[0.162888,0.789393,-0.743539,0.263129,-0.0253093],[-0.130356,0.673665,-0.959456,0.347247,-0.0267758],[0.907919,1.207552,-0.687343,0.213486,-0.0822904],[-0.551119,0.0670665,-0.497089,0.100754,0.0602253],[0.146543,-0.0843370,0.195286,-0.032932,-0.0202595]]
    mu_1 = 0
    
    for i in range(len(b)):
        for j in range(len(b[i])):
            mu_temp += rho_bar * b[i][j] * (1/T_bar -1 ) ** i * (rho_bar -1) ** j
    
    mu_1 = math.exp(mu_temp)
    # multiply the functions to return the viscosity
    viscosity = mu_o * mu_1
    
    logger.info('Computed dynamic (absolute) viscosity of water as %s kg/m-s at T=%S degrees C and P = %s Pa % viscosity,temperature,pressure') 
    
    logger.debug('Computed dynamic (absolute) viscosity of water using empirical NIST equation described in Sengers, J.V. "Representative Equations for the Viscosity of Water Substance." J. Phys. Chem. Ref. Data 13(1), 1984.')
    
    return viscosity


def water_viscosity_kinematic(temperature=25,pressure=101325):
    '''
    Return the kinematic viscosity of water in m2/s = Stokes
    at the specified temperature.
    
    Parameters:
    ----------
    temperature : float or int, optional
                  The temperature in Celsius. Defaults to 25 degrees if not specified.
    pressure    : float or int, optional
                  The ambient pressure of the solution in Pascals (N/m2). 
                  Defaults to atmospheric pressure (101325 Pa) if not specified.
                  
    Returns:
    -------
    float
            The kinematic viscosity of water in Stokes (m2/s)
            
    See Also:
    --------
    water_viscosity_dynamic
    water_density
    
    '''
    kviscosity = water_viscosity_dynamic(temperature,pressure) / water_density(temperature,pressure)
    logger.info('Computed kinematic viscosity of water as %s m2/s at T=%S degrees C and P = %s Pa' % (kviscosity,temperature,pressure)) 
    return kviscosity
    

def water_dielectric_constant(temperature=25):
    '''(number) -> float
    
    Return the dielectric constant of water at the specified temperature.
    
    Parameters:
    ----------
    temperature : float or int, optional
                  The temperature in Celsius. Must be between 0 and 74 C. See 
                  notes. Defaults to 25 degrees if not specified.
                  
    Returns:
    -------
    float
            The dielectric constant (or permittivity) of water relative to the
            permittivity of a vacuum. Dimensionless.
    
    Notes:
    -----
    This function implements a quadratic fit of measured permittivity data as
    reported in the CRC Handbook[1]. The parameters given are valid over the
    range 273 K to 372 K. Permittivity should not be extrapolated beyond this
    range.
    
    :math: \epsilon(T) = a + b T + c T^2
    
    .. [1] "Permittivity (Dielectric Constant) of Liquids." CRC Handbook of 
            Chemistry and Physics, 92nd ed, pp 6-187 - 6-208.
    
    Examples:
    --------
    >>> water_dielectric_constant(20) #doctest: +ELLIPSIS
    80.15060...
    
    Display an error if 'temperature' is outside the valid range
    
#     TODO >>> water_dielectric_constant(-5)
#     ERROR: Temperature specified exceeds range of data. Cannot extrapolate dielectric constant.
#     
    '''
    # do not return anything if 'temperature' is outside the range for which
    # this fit applies
    if kelvin(temperature) < 273 or kelvin(temperature) > 372:
        logger.error('Specified temperature (%s K) exceeds valid range of data. Cannot extrapolate. % kelvin(temperature)')
        return None
    
    # otherwise, calculate the dielectric constant using the quadratic fit    
    a = 0.24921e3
    b = -0.79069e0
    c = 0.72997e-3
    dielectric = a + b * kelvin(temperature) + c * kelvin(temperature) ** 2
    
    logger.info('Computed dielectric constant of water as %s at %s degrees Celsius' % (dielectric,temperature))
    
    logger.debug('Computed dielectric constant of water using empirical equation given in "Permittivity (Dielectric Constant) of Liquids." CRC Handbook of Chemistry and Physics, 92nd ed, pp 6-187 - 6-208.')
    
    return dielectric
    
    
def water_conductivity(temperature):
    pass

def water_activity():
    pass
    
def water_debye_parameter_activity(temperature=25):
    '''(number) -> float
    return the constant A for use in the Debye-Huckel limiting law (base 10)
    
    Parameters:
    ----------
    temperature : float or int, optional
                  The temperature in Celsius. Defaults to 25 degrees if not specified.
    
    Returns:
    -------
    float          The parameter A for use in the Debye-Huckel limiting law (base 10)
    
    Notes:
    -----
    ### TODO - FIX THIS TO INCLUDE DENSITY in kg/m3
     The parameter A is equal to:[1]
         
     $$ A = {e^3 \over 8 \pi} sqrt{ 2 N_a \over (\epsilon_r \epsilon_o k_B T) ^3}

      This should not be confused with the Debye-Huckel constant for osmotic coefficients.
     
     TODO - FIND MORE CREDIBLE REFERENCE
     .. [1] http://en.wikipedia.org/wiki/Debye%E2%80%93H%C3%BCckel_equation
     
     
     
     The parameter A is equal to:[1]
     
     $$ A = 1.82e6 (\epsilon_r T) ^ -3/2 $$
    
    Note that when used in conjunction with the Debye-Huckel limiting law or related equations,
    this parameter is valid only when ionic strength is calculated from molar (mol/L) scale concentrations.
    
    .. [1] Stumm, Werner and Morgan, James J. Aquatic Chemistry, 3rd ed, 
        pp 103. Wiley Interscience, 1996.
        
    Examples:
    --------
    >>> water_debye_parameter_activity()
    0.509...
    
    See Also:
    --------
    water_debye_parameter_osmotic
    
    '''
#       this works, but I can't figure out how to reconcile the units with density included
#     return CONST_e ** 3 / (8 * math.pi) * math.sqrt(2 * CONST_Na / ( 
#     (water_dielectric_constant(temperature) * CONST_Eo * CONST_kb * 
#     kelvin(temperature)) ** 3)) * math.sqrt(water_density(temperature)) * math.log10(math.e)
    
    # use this from Stumm and Morgan Instead
    debyeparam = 1.8246e6 * (water_dielectric_constant(temperature) * kelvin(temperature)) ** -1.5
    logger.info('Computed Debye-Huckel Limiting Law Constant A = %s at %s degrees Celsius' % (debyeparam,temperature))
    return debyeparam
    
    # or this from MWH treatment book, page 302
    # TODO - document reference
    #return 1.29e6 * math.sqrt(2) * (water_dielectric_constant(temperature) * kelvin(temperature)) ** -1.5

def water_debye_parameter_osmotic(temperature=25):
    '''(number) -> float
    return the constant A_phi for use in calculating the osmotic coefficient according to Debye-Huckel theory
    
    Parameters:
    ----------
    temperature : float or int, optional
                  The temperature in Celsius. Defaults to 25 degrees if not specified.
    
    Notes:
    -----
    Not to be confused with the Debye-Huckel constant used for activity coefficients in the limiting law. Takes the value 0.392 at 25 C.
    This constant is calculated according to:[1]

     $$ A_{phi} = {1 \over 3} ({ 2 \pi N_A \rho_w \over 1000})^0.5 * ({e^2 \over \epsilon_o \epsilon_r k T})^1.5 $$
    
    
    .. [1] Kim, Hee-Talk and Frederick, William Jr, 1988. "Evaluation of Pitzer Ion Interaction Parameters of Aqueous Electrolytes at 25 C. 1. Single Salt Parameters,"
    //J. Chemical Engineering Data// 33, pp.177-184.
    
    Examples:
    --------
    >>> water_debye_parameter_osmotic() #doctest: +ELLIPSIS
    0.3920009...
    
    See Also:
    --------
    water_debye_parameter
    
    '''
    # TODO - the factor 0.710 is a mystery number needed to make this equation return the correct value at 25C. I don't know why
    return 0.71049 * 1/3 * (2 * math.pi * CONST_Na * water_density(temperature) / 1000 ) ** 0.5 * ( CONST_e ** 2 / (water_dielectric_constant(temperature) * CONST_Eo * CONST_kb * kelvin(temperature) ) ) ** 1.5


## Acid - Base Functions


def p(x):
    ''' (number) -> float
    Negative log of x. Generally used for expressing concentration of hydrogen
    ions (pH)
    
    Parameters:
    ----------
    x : float or int
        Any number (usually a species concentration)
    
    Returns:
    -------
    float
        The negative log10 of the input value.
        
    Examples:
    --------
    >>> p(1e-7)
    7.0
    >>> p(1.568e-9)
    8.80465394165158
    
    '''
    return -1 * math.log10(x)
    

def alpha(n,pH,pKa_list):
    '''(int,number,list of numbers)
    Returns the acid-base distribution coefficient (alpha) of an acid in the 
    n-deprotonated form at a given pH.
    
    Parameters:
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
    
    Returns:
    -------
    float
        The fraction of total acid present in the specified form.
    
    Notes:
    -----
    The acid-base distribution coefficient is calculated as follows:[1]
        
    .. math::
        \alpha_n = {term_n \over [H+]^n + k_{a1}[H+]^n-1 + k_{a1}k_{a2}[H+]^n-2 ... k_{a1}k_{a2}...k_{an} }
    
    Where :math: '\term_n' refers to the nth term in the denominator, starting from 0
    
    .. [1] Stumm, Werner and Morgan, James J. Aquatic Chemistry, 3rd ed, 
        pp 127-130. Wiley Interscience, 1996.
    
    Examples:
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

### Functions that operate on Solution Objects

def gibbsmix(Solution1, Solution2,temperature=25):
    '''(Solution, Solution) -> float
    Return the Gibbs energychange associated with mixing two solutions

    Parameters:
    ----------
    Solution1, Solution2 : Solution objects
        The two solutions to be mixed.
    temperature : float or int, optional
                  The temperature in Celsius. Defaults to 25 degrees if not specified.
        
    Returns:
    -------
    float
        The change in Gibbs eneryg associated with complete mixing of the
        Solutions, in Joules.
    
    Notes:
    -----
    
    The Gibbs energy of mixing is calculated as follows:[1]
        
    .. math::
        \Delta_{mix} G = \sum_i (n_c + n_d) R T \ln a_b - \sum_i n_c R T \ln a_c - \sum_i n_d R T \ln a_d
    
    Where n is the number of moles of substance, T is the temperature in kelvin,
    and  subscripts b, c, and refer to the concentrated, dilute, and blended
    Solutions, respectively. 
    
    Note that dissociated ions must be counted as separate components,
    so a simple salt dissolved in water is a three component solution (cation,
    anion, and water).
    
    .. [1] Koga, Yoshikata, 2007. //Solution Thermodynamics and its Application to Aqueous Solutions: 
    A differential approach.// Elsevier, 2007, pp. 23-37.
    
    Examples:
    --------
 
    '''
    concentrate = Solution1
    dilute = Solution2
    blend = mix(Solution1,Solution2)
    term_list = {concentrate:0, dilute:0, blend:0}

    # calculte the entropy change and number of moles solute for each solution
    for solution in term_list:
        for solute in solution.components:
            #print(solution.list_concentrations())
            if not solution.get_amount(solute,'fraction') == 0:
                term_list[solution] += solution.get_amount(solute,'mol') * math.log(solution.get_activity(solute))

    return CONST_R * kelvin(temperature) * (term_list[blend] - term_list[concentrate] - term_list[dilute])

def entropy_mix(Solution1, Solution2,temperature=25):
    '''(Solution, Solution) -> float
    Return the ideal mixing entropy associated with mixing two solutions

    Parameters:
    ----------
    Solution1, Solution2 : Solution objects
        The two solutions to be mixed.
    temperature : float or int, optional
                  The temperature in Celsius. Defaults to 25 degrees if not specified.
        
    Returns:
    -------
    float
        The ideal mixing entropy associated with complete mixing of the
        Solutions, in Joules.
    
    Notes:
    -----
    
    The ideal entropy of mixing is calculated as follows:[1]
        
    .. math::
        \Delta_{mix} S = \sum_i (n_c + n_d) R T \ln x_b - \sum_i n_c R T \ln x_c - \sum_i n_d R T \ln x_d
    
    Where n is the number of moles of substance, T is the temperature in kelvin,
    and  subscripts b, c, and refer to the concentrated, dilute, and blended
    Solutions, respectively. 
    
    Note that dissociated ions must be counted as separate components,
    so a simple salt dissolved in water is a three component solution (cation,
    anion, and water).
    
    .. [1] Koga, Yoshikata, 2007. //Solution Thermodynamics and its Application to Aqueous Solutions: 
    A differential approach.// Elsevier, 2007, pp. 23-37.
    
    Examples:
    --------
 
    '''
    concentrate = Solution1
    dilute = Solution2
    blend = mix(Solution1,Solution2)
    term_list = {concentrate:0, dilute:0, blend:0}

    # calculte the entropy change and number of moles solute for each solution
    for solution in term_list:
        #mole_fraction_water = solution.get_moles_water() / (solution.get_total_moles_solute() + solution.get_moles_water())
        for solute in solution.components:
            #print(solution.list_concentrations())
            if not solution.get_amount(solute,'fraction') == 0:
                term_list[solution] += solution.get_amount(solute,'mol') * math.log(solution.get_amount(solute,'fraction'))

    return CONST_R * kelvin(temperature) * (term_list[blend] - term_list[concentrate] - term_list[dilute])


def mix(Solution1, Solution2):
    '''(Solution, Solution) -> Solution
    Returns a new Solution object that results from the mixing of Solution1
    and Solution2
    
    '''
    # check to see if the two solutions have the same solvent
    if not Solution1.solvent_name == Solution2.solvent_name:
        logger.error('mix() function does not support solutions with different solvents. Aborting.')
    
    if not Solution1.solvent_name == 'H2O' or Solution1.solvent_name == 'water':
        logger.error('mix() function does not support non-water solvents. Aborting.')
    
    # set solvent as water and retrieve total mass from solutions
    solvent_mass = Solution1.get_solvent_mass() + Solution2.get_solvent_mass()
    solvent = [Solution1.solvent_name,18,solvent_mass]
    
    # determine the volume for the new solution (assume additive)
    mix_vol = Solution1.get_volume() + Solution2.get_volume()
    
    # determine the density for the new solution. Assume volume is additive.
    # Convert from kg/L to kg/m3
    mix_dense = (Solution1.get_mass() + Solution2.get_mass()) / mix_vol *1000
    
    #conductivity will be initialized to zero in the new solution TODO
    
    #temperature and pressure will be left at defaults TODO

    # Loop through all the solutes in Solution1 and Solution2 and make a list
    mix_solute_list=[]
    for item in Solution1.components.keys():
            if not item in mix_solute_list:
                mix_solute_list.append(item)
    for item in Solution2.components.keys():
        if not item in mix_solute_list:
            mix_solute_list.append(item)
            
    # add each solute to the new solution
    component_list=[]
    for item in mix_solute_list:
        if item in Solution1.components and item in Solution2.components:
            mix_moles = Solution1.components[item].get_moles() + Solution2.components[item].get_moles()
            
            component_list.append([Solution1.get_solute(item).get_name(),Solution1.get_solute(item).get_molecular_weight(),mix_moles,'mol',Solution1.get_solute(item).parameters])
        
        elif item in Solution1.components:
            mix_moles = Solution1.components[item].get_moles()
        
            component_list.append([Solution1.get_solute(item).get_name(),Solution1.get_solute(item).get_molecular_weight(),mix_moles,'mol',Solution1.get_solute(item).parameters])
                        
        elif item in Solution2.components:
            mix_moles = Solution2.components[item].get_moles()
            
            component_list.append([Solution2.get_solute(item).get_name(),Solution2.get_solute(item).get_molecular_weight(),mix_moles,'mol',Solution2.get_solute(item).parameters])
    
    #create a new Solution object for the blend
    Blend = Solution(component_list,solvent,mix_dense)
        
    #copy over any activity coefficient parameters and others TODO
    for item in Blend.components:
        if item in Solution1.components:
            if Solution1.components[item].parameters_TCPC:
                Blend.components[item].set_parameters_TCPC(Solution1.components[item].get_parameters_TCPC('S'),Solution1.components[item].get_parameters_TCPC('b'),Solution1.components[item].get_parameters_TCPC('n'),Solution1.components[item].get_parameters_TCPC('z_plus'),Solution1.components[item].get_parameters_TCPC('z_minus'),Solution1.components[item].get_parameters_TCPC('nu_plus'),Solution1.components[item].get_parameters_TCPC('nu_minus'))

        elif item in Solution2.componenents:
            if Solution2.components[item].parameters_TCPC:
                Blend.components[item].set_parameters_TCPC(Solution2.components[item].get_parameters_TCPC('S'),Solution2.components[item].get_parameters_TCPC('b'),Solution2.components[item].get_parameters_TCPC('n'),Solution2.components[item].get_parameters_TCPC('z_plus'),Solution2.components[item].get_parameters_TCPC('z_minus'),Solution2.components[item].get_parameters_TCPC('nu_plus'),Solution2.components[item].get_parameters_TCPC('nu_minus'))

    return Blend
    
### Activity Functions
# Individual functions for activity coefficients are defined here so that they can be used independently of a 
# pyEQL solution object. Normally, these functions are called from within the get_activity_coefficient method of
# the Solution class.

def get_activity_coefficient_debyehuckel(ionic_strength,valence=1,temperature=25):
    '''Return the activity coefficient of solute in the parent solution according to the Debye-Huckel limiting law.
    
    Parameters:
    ----------
    valence : int, optional      
                        The charge on the solute, including sign. Defaults to +1 if not specified.
    ionic_strength : number
                        The ionic strength of the parent solution, mol/kg
    temperature : float or int, optional
                        The solution temperature in degrees Celsius. 
                        Defaults to 25 degrees if omitted.
    Returns:
    -------
    float
         The mean molal (mol/kg) scale ionic activity coefficient of solute

    See Also:
    --------
    water_debye_parameter_activity
    get_ionic_strength
    
    Notes:
    ------
    Valid only for I < 0.005
    
    .. [1] Stumm, Werner and Morgan, James J. Aquatic Chemistry, 3rd ed, 
    pp 103. Wiley Interscience, 1996.
    '''
    # check if this method is valid for the given ionic strength
    if not ionic_strength < 0.005:
        logger.warning('Ionic strength exceeds valid range of the Debye-Huckel limiting law')
    
    log_f = - water_debye_parameter_activity(temperature) *valence ** 2 * math.sqrt(ionic_strength)
    
    return 10 ** log_f

def get_activity_coefficient_guntelberg(ionic_strength,valence=1,temperature=25):
    '''Return the activity coefficient of solute in the parent solution according to the Guntelberg approximation.
    
    Parameters:
    ----------
    valence : int, optional          
                        The charge on the solute, including sign. Defaults to +1 if not specified.
    ionic_strength : number
                        The ionic strength of the parent solution, mol/kg
    temperature : float or int, optional
                        The solution temperature in degrees Celsius. 
                        Defaults to 25 degrees if omitted.
    Returns:
    -------
    float
         The mean molal (mol/kg) scale ionic activity coefficient of solute
         
    See Also:
    --------
    water_debye_parameter_activity
    get_ionic_strength
    
    Notes:
    ------
    Valid for I < 0.1
    
    .. [1] Stumm, Werner and Morgan, James J. Aquatic Chemistry, 3rd ed, 
    pp 103. Wiley Interscience, 1996.
    '''
    # check if this method is valid for the given ionic strength
    if not ionic_strength < 0.1:
        logger.warning('Ionic strength exceeds valid range of the Guntelberg approximation')
    
    log_f = - water_debye_parameter_activity(temperature) * valence ** 2 * math.sqrt(ionic_strength) / (1+math.sqrt(ionic_strength))

    return 10 ** log_f
    
def get_activity_coefficient_davies(ionic_strength,valence=1,temperature=25):
    '''Return the activity coefficient of solute in the parent solution according to the Davies equation.
    
    Parameters:
    ----------
    valence : int, optional           
                        The charge on the solute, including sign. Defaults to +1 if not specified.
    ionic_strength : number
                        The ionic strength of the parent solution, mol/kg
    temperature : float or int, optional
                        The solution temperature in degrees Celsius. 
                        Defaults to 25 degrees if omitted.
    Returns:
    -------
    float
         The mean molal (mol/kg) scale ionic activity coefficient of solute

    See Also:
    --------
    water_debye_parameter_activity
    get_ionic_strength
    
    Notes:
    ------
    Valid for 0.1 < I < 0.5
    
    .. [1] Stumm, Werner and Morgan, James J. Aquatic Chemistry, 3rd ed, 
    pp 103. Wiley Interscience, 1996.
    '''
    # check if this method is valid for the given ionic strength
    if not ionic_strength < 0.5 and ionic_strength > 0.1:
        logger.warning('Ionic strength exceeds valid range of the Davies equation')
    
    log_f = - water_debye_parameter_activity(temperature) * valence ** 2 * ( math.sqrt(ionic_strength) / (1+math.sqrt(ionic_strength)) - 0.2 * ionic_strength)
    
    return 10 ** log_f
    
def get_activity_coefficient_TCPC(ionic_strength,S,b,n,valence=1,counter_valence=-1,stoich_coeff=1,counter_stoich_coeff=1,temperature=25):
    '''Return the activity coefficient of solute in the parent solution according to the modified TCPC model.
    
    Parameters:
    ----------
    ionic_strength : number
                        The ionic strength of the parent solution, mol/kg
    S : float
                        The solvation parameter for the parent salt. See Reference.
    b : float
                        The approaching parameter for the parent salt. See Reference.
    n : float       
                        The n parameter for the parent salt. See Reference.
    valence : int, optional           
                        The charge on the solute, including sign. Defaults to +1 if not specified.
    counter_valence : int, optional           
                        The charge on the solute's complementary ion, including sign. Defaults to -1 if not specified.
                        E.g. if the solute is Na+ and the salt is NaCl, counter_valence = -1
    stoich_coeff : int, optional
                        The stoichiometric coefficient of the solute in its parent salt. Defaults to1 if not specified.
                        E.g. for Zn+2 in ZnCl2, stoich_coeff = 1
    counter_stoich_coeff : int, optional
                        The stoichiometric coefficient of the solute's complentary ion in its parent salt. Defaults to 1 if not specified.
                        E.g. for Cl- in ZnCl2, stoich_coeff = 2
    temperature : float or int, optional
                        The solution temperature in degrees Celsius. 
                        Defaults to 25 degrees if omitted.
    Returns:
    -------
    float
        The mean molal (mol/kgL) scale ionic activity coefficient of solute

    See Also:
    --------
    water_debye_parameter_osmotic
    get_ionic_strength
    
    Notes:
    ------
    Valid for concentrated solutions up to saturation. Accuracy compares well with the Pitzer approach. See Reference [1] for a compilation of the appropriate parameters for a variety of commonly-encountered electrolytes.
    
    .. [1] Ge, Xinlei, Wang, Xidong, Zhang, Mei, and Seetharaman, Seshadri. "Correlation and Prediction of Activity and Osmotic Coefficients of Aqueous Electrolytes at 298.15 K by the Modified TCPC Model." J. Chemical Engineering Data 52, pp.538-547, 2007.
    '''
    # compute the PDF parameter
    PDH = - math.fabs(valence * counter_valence) * water_debye_parameter_osmotic(temperature) * ( ionic_strength ** 0.5 / (1 + b * ionic_strength ** 0.5) + 2/b * math.log(1 + b * ionic_strength ** 0.5))
    # compute the SV parameter
    SV = S / kelvin(temperature) * ionic_strength  ** (2*n) / (stoich_coeff + counter_stoich_coeff)
    # add and exponentiate to eliminate the log
    return math.exp(PDH + SV)
    
def get_activity_coefficient_pitzer():
    '''Return the activity coefficient of solute in the parent solution according to the Pitzer model.
    
    Returns:
    -------
    float
        The mean molal (mol/kg) scale ionic activity coefficient of solute
    
    See also:
    --------
    water_debye_parameter_activity
    
    '''
    pass

def get_osmotic_coefficient_TCPC(ionic_strength,S,b,n,valence=1,counter_valence=-1,stoich_coeff=1,counter_stoich_coeff=1,temperature=25):
    '''Return the osmotic coefficient of solute in the parent solution according to the modified TCPC model.
    
    Parameters:
    ----------
    ionic_strength : number
                        The ionic strength of the parent solution, mol/kg
    S : float
                        The solvation parameter for the parent salt. See Reference.
    b : float
                        The approaching parameter for the parent salt. See Reference.
    n : float       
                        The n parameter for the parent salt. See Reference.
    valence : int, optional           
                        The charge on the solute, including sign. Defaults to +1 if not specified.
    counter_valence : int, optional           
                        The charge on the solute's complementary ion, including sign. Defaults to -1 if not specified.
                        E.g. if the solute is Na+ and the salt is NaCl, counter_valence = -1
    stoich_coeff : int, optional
                        The stoichiometric coefficient of the solute in its parent salt. Defaults to1 if not specified.
                        E.g. for Zn+2 in ZnCl2, stoich_coeff = 1
    counter_stoich_coeff : int, optional
                        The stoichiometric coefficient of the solute's complentary ion in its parent salt. Defaults to 1 if not specified.
                        E.g. for Cl- in ZnCl2, stoich_coeff = 2
    temperature : float or int, optional
                        The solution temperature in degrees Celsius. 
                        Defaults to 25 degrees if omitted.
    Returns:
    -------
    float
        The osmotic coefficient of the solute

    See Also:
    --------
    water_debye_parameter_osmotic
    get_ionic_strength
    
    Notes:
    ------
    Valid for concentrated solutions up to saturation. Accuracy compares well with the Pitzer approach. See Reference [1] for a compilation of the appropriate parameters for a variety of commonly-encountered electrolytes.
    
    .. [1] Ge, Xinlei, Wang, Xidong, Zhang, Mei, and Seetharaman, Seshadri. "Correlation and Prediction of Activity and Osmotic Coefficients of Aqueous Electrolytes at 298.15 K by the Modified TCPC Model." J. Chemical Engineering Data 52, pp.538-547, 2007.
    '''
    # compute the 2nd term
    term2 = - math.fabs(valence * counter_valence) * water_debye_parameter_osmotic(temperature) * ionic_strength ** 0.5 / (1 + b * ionic_strength ** 0.5)
    # compute the 3rd term
    term3 = S / (kelvin(temperature) * ( stoich_coeff + counter_stoich_coeff)) * 2 * n / (2 * n + 1) * ionic_strength  ** (2 * n)
    # add and return the osmotic coefficient
    return 1 - term2 + term3

### Other Stuff - TODO

def mobility(diffusion_coefficient,valence,temperature=25):
    '''Return the ionic mobility of a species
    
    Parameters:
    ----------
    diffusion_coefficient : float
                The diffusion coefficient of the species in m2/s
    valence : float or int
                The charge on the species, including sign
    temperature : float or int, optional
                The solution temperature in degrees Celsius. 
                Defaults to 25 degrees if omitted.
    
    Returns:
    -------
    float : the ionic mobility in m2/V-s
    
    
    Notes:
    -----
    This function uses the Einstein relation to convert a diffusion coefficient
    into an ionic mobility[1]
    
    .. math::
        \mu_i = {F |z_i| D_i \over RT}
    
    .. [1] Smedley, Stuart I. The Interpretation of Ionic Conductivity in Liquids. Plenum Press, 1980.
    
    '''
    mobility = CONST_F * math.fabs(valence) * diffusion_coefficient / (CONST_R * kelvin(temperature))
    
    logger.info('Computed ionic mobility as %s m/s-V from D = %s m2/s at T=%S degrees C % mobility,diffusion_coeffiicent,temperature')
    
    return mobility

def debye_length(dielectric_constant,ionic_strength,temp):
    '''(number,number,number) -> float
    Return the Debye length of a solution in meters
    
    dielectric_constant is the dielectric constant of the solution
    ionic_strength is the ionic strength in moles per cubic meter
    temp is the temperature in degrees celsius
    
    '''
    return math.sqrt( dielectric_constant * CONST_Eo * CONST_R * kelvin(temp) / (2 * CONST_F ** 2 * ionic_strength) )

#####CLASS DEFINITIONS
'''

Idiom for Class Methods
-----------------------

All class methods begin with a verb such as "get," "set," or "calc." The Solution
class provides a few methods with "add" and "remove" methods as well. 
The following meaning is assigned to these verbs:

get_ - Method returns an attribute from the object. Object data is not modified.

list_ - Method returns a human-readable representation of some object data. Not
        intended for use with functions.

set_ - Method modifies an object's attribute. These methods are often used to 
       override default parameter values, or replace calculated values with 
       more precise data. Does not return a value.
        
calc_ - Method calculates and updates an object's attribute using data internal
        to the object and possibly additional arguments. Often used for triggering
        updates to dynamic attributes, such as species concentrations before and
        after a chemical reaction.
    
add_ - Method adds data to an existing object. Typically used to add chemical
       species to a Solution object.
       
remove_ - Method removes data from an existing object. Typically used to remove
          a chemical species from a Solution object.

'''


class Solution:
    '''Class representing the properties of a solution. Instances of this class contain information about the solutes, solvent, and bulk properties.
    
    Parameters:
    ----------
    solutes : list of lists, optional
                See add_solute() documentation for formatting of this list
                Defaults to empty (pure solvent) if omitted
    volume : pint quantity, optional
                Volume of the solution. Defaults to 1L if omitted.
    solvent : list of lists, optional
                String representing the chemical formula of the solvent. 
                Defaults to 1 kg water if omitted.
    temperature : pint quantity, optional
                The solution temperature. Defaults to 25 degrees C if omitted.
    pressure : pint quantity
                The ambient pressure of the solution. 
                Defaults to 1 atm if omitted.
    
    Returns:
    -------
    A Solution object.
    
    Examples:
    --------
    # Defining a 0.5M NaCl solution
    >>> solutes = [['Na+',23,0.0115],['Cl-',35,0.0175]]
    >>> solvent = ['H2O',18,1]
    >>> my_solution = Solution(solutes,solvent)
    >>> print(my_solution)
    Components: ['Na+', 'H2O', 'Cl-']
    Volume: 1.0290000000000001L
    Density: 1000 kg/m3
    
    
    See Also:
    --------
    add_solute
    
    '''
    
    '''THE PLAN FOR SOLUTION INIT
    A) Solvent: specify mass and bulk density Solutes: amount per mass units
    Calculate total solution volume
    
    OR 
    
    B) Solvent: total solution volume and density. Solutes: amount per volume units
    Calculate solvent mass
    
    conductivity is calculated from solutes/database and can be directly set with a method
    When solvent=H2O, pH is calculated through speciation / reaction and can be directly set with a method
    
    temperature is always set as a bulk property
    
    
    
    '''
    
    def __init__(self,solutes=[],solvent=['H2O','0.998 kg'],volume='1 L',temperature='25 degC',pressure='1 atm'):
        self.volume = unit(volume)
        self.temperature = unit(temperature)
        self.pressure = unit(pressure)
        self.components={}
        self.solvent_name=solvent[0]
        
        # warn if the solvent is anything besides water
        if not solvent[0] == 'H2O' or solvent[0] == 'water' :
            logger.error('Non-aqueous solvent detected. These are not yet supported!')
        # define the solvent
        self.add_solvent(*solvent)
        
        # populate the solutes
        for item in solutes:
            self.add_solute(*item)
        
        # TODO - move the code below to a separate method to calculate conductivity
#        # calculate the conductivity based on the molar conductivity of the solutes
#        '''
#        Conductivity is calculated by summing the molar conductivities of the respective
#        solutes, but they are activity-corrected and adjusted using an empricial exponent.
#        This approach is used in PHREEQC and Aqion as described at these URLs:        
#        <http://www.aqion.de/site/77>        
#        <http://www.hydrochemistry.eu/exmpls/sc.html>
#          
#        
#        '''
#        EC = 0
#        
#        for item in solutes:
#            z = math.abs(self.get_formal_charge)
#            # determine the value of the exponent alpha 
#            if self.get_ionic_strength < 0.36 * z:
#                alpha = 0.6 / z ** 0.5
#            else:
#                alpha = self.get_ionic_strength ** 0.5 / z
#                
#            EC += item.get_molar_conductivity * self.get_activity_coefficient(item) ** alpha * self.get_amount(item,'mol/L')
#             
#             
#        self.conductivity = EC
        

    def add_solute(self,formula,amount,parameters={}):
        '''Primary method for adding substances to a pyEQL solution
        
        Parameters:
        ----------
        formula : str
                    Chemical formula for the solute. 
                    Charged species must contain a + or - and (for polyvalent solutes) a number representing the net charge (e.g. 'SO4-2').
        amount : str
                    The amount of substance in the specified unit system. The string should contain both a quantity and
                    a pint-compatible representation of a unit. e.g. '5 mol/kg' or '0.1 g/L'
        parameters : dictionary, optional
                    Dictionary of custom parameters, such as diffusion coefficients, transport numbers, etc. Specify parameters as key:value pairs separated by commas within curly braces, e.g. {diffusion_coeff:5e-10,transport_number:0.8}. The 'key' is the name that will be used to access the parameter, the value is its value.
                        
        '''        
        new_solute = self.Solute(formula,amount,self.get_volume(),self.get_solvent_mass(),parameters)
        self.components.update({new_solute.get_name():new_solute})
    
    def add_solvent(self,formula,amount,parameters={}):
        '''Same as add_solute but omits the need to pass solvent mass to pint
        '''
        new_solvent = self.Solute(formula,amount,self.get_volume(),amount,parameters)
        self.components.update({new_solvent.get_name():new_solvent})
                        
    def get_solute(self,i):
        return self.components[i]
    
    def get_solvent(self):
        return self.components[self.solvent_name]
    
    class Solute:
        '''represent each chemical species as an object containing its formal charge, transport numbers, concentration, activity, etc. 
        
        '''
    
        def __init__(self,formula,amount,volume,solvent_mass,parameters={}):
            '''
            Parameters:
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
            import chemical_formula as chem
            
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
                            
                # TODO - deprecate in favor of the parameter module
                self.parameters_TCPC={}
                
                # trigger the function that checks whether parameters already exist for this species, and if not,
                # searches the database files and creates them
                database.search_parameters(self.formula)

        def get_parameter(self,parameter,temperature=None,pressure=None,ionic_strength=None):
            '''
            Return the value of the parameter named 'parameter'
            
            '''
            # Search for this solute in the database of parameters
            if self.formula in db:
                # loop over the parameters in the python set associated with this solute species to find the one 
                #we're looking for
            
                # TODO - handle the case where multiple parameters with same name exist. Need function
                # to compare specified conditions and choose the most appropriate parameter
            
                for item in db[self.formula]:
                    if item.get_name() == parameter:
                        return item.get_value(temperature,pressure,ionic_strength)
                # TODO - fix this bug - throws a log message at every entry
                    else:
                        logger.error('Required parameter %s not found for species %s.' % (parameter,self.formula))
            else:
                logger.error('No entry for species %s in parameters database' % self.formula)
                    
        
        # TODO - deprecate in favor of the parameter module
        def set_parameters_TCPC(self,S,b,n,valence=1,counter_valence=-1,stoich_coeff=1,counter_stoich_coeff=1):
            '''Use this function to store parameters for the TCPC activity model
            
            Parameters:
            ----------
            S : float
                            The solvation parameter for the parent salt. See Reference.
            b : float
                                The approaching parameter for the parent salt. See Reference.
            n : float       
                                The n parameter for the parent salt. See Reference.
            valence : int, optional           
                                The charge on the solute, including sign. Defaults to +1 if not specified.
            counter_valence : int, optional           
                                The charge on the solute's complementary ion, including sign. Defaults to -1 if not specified.
                                E.g. if the solute is Na+ and the salt is NaCl, counter_valence = -1
            stoich_coeff : int, optional
                                The stoichiometric coefficient of the solute in its parent salt. Defaults to1 if not specified.
                                E.g. for Zn+2 in ZnCl2, stoich_coeff = 1
            counter_stoich_coeff : int, optional
                                The stoichiometric coefficient of the solute's complentary ion in its parent salt. Defaults to 1 if not specified.
                                E.g. for Cl- in ZnCl2, stoich_coeff = 2
            
            Returns:
            -------
            No return value. Parameter values are stored in a dictionary.
            
            See Also:
            get_parameters_TCPC
            
            '''
            self.parameters_TCPC={'S':S,'b':b,'n':n,'z_plus':valence,'z_minus':counter_valence,'nu_plus':stoich_coeff,'nu_minus':counter_stoich_coeff}
        
        # TODO - deprecate in favor of the parameter module
        def get_parameters_TCPC(self,name):
            '''Retrieve modeling parameters used for activity coefficient modeling
            
            Parameters:
            ----------
            name : str
                        String identifying the specific parameter to be retrieved. Must correspond to a key in the 'name' dictionary
            
            Returns:
            -------
            The parameter stored at key 'name' in the 'model' dictionary
            
            '''
            return self.parameters_TCPC[name]
        
        
        def get_name(self):
            return self.formula
            
        def get_formal_charge(self):
            return self.charge
            
        def get_molecular_weight(self):
            return self.mw
        
        def get_moles(self):
            return self.moles
        
        def set_moles(self,moles):
            self.moles = moles
            
        def get_activity(self):
            return self.activity
        
        def set_activity(self,activity):
            self.activity = activity

        def get_molar_conductivity(self,temperature=25):
            # TODO - requires diffusion coefficient which may not be present
            '''(float,int,number) -> float
            Calculate the molar (equivalent) conductivity for a species in 
            Siemens-meters^2/mole
            
            Parameters:
            ----------
            temperature : float or int, optional
                                    The solution temperature in degrees Celsius. 
                                    Defaults to 25 degrees if omitted.
                
            Returns:
            -------
            float
                    The molar or equivalent conductivity of the species, 
                    in Siemens-meters^2/mole
            
            Notes:
            -----
            Molar conductivity is calculated from the Nernst-Einstein relation:[1]
                
            .. math::
                \DELTA_i = {z_i^2 D_i F^2 \over RT}
            
            Note that the diffusion coefficient is strongly variable with temperature.
            
            .. [1] Smedley, Stuart. The Interpretation of Ionic Conductivity in Liquids, pp 1-9. Plenum Press, 1980.
            
            TODO Examples:
#             --------
#             For sodium ion:
#                 
#             >>> molar_conductivity(1.334e-9,1) #doctest: +ELLIPSIS
#             0.0050096...
#             
#             For sulfate ion at 30 C:
#                 
#             >>> molar_conductivity(1.065e-9,2,30) #doctest: +ELLIPSIS
#             0.0157340...
            '''
            return diffusion_coefficient * CONST_F ** 2 * self.get_formal_charge() ** 2 / (CONST_R * kelvin(temperature))
        
            
        #set output of the print() statement
        def __str__(self):
            return 'Species ' + str(self.get_name()) + ' MW=' + str(self.get_molecular_weight()) +' Valence='+str(self.get_formal_charge()) + ' Amount= ' + str(self.get_moles()) + 'moles  Activity= ' + str(self.get_activity())
    
    class Solvent:
        '''subclass of Solute. Adds density'''
        pass 
    
    def get_solvent_mass(self):
        # return the total mass (kg) of the solvent
        solvent = self.get_solvent()
        mw = solvent.get_molecular_weight()
    
        return solvent.get_moles().to('kg','chem',mw=mw)
            
    def get_volume(self,temperature=25):
        return self.volume
    
    def get_mass(self):
        '''returns the total solution mass in kg'''
        total_mass = 0
        for item in self.components:
            total_mass+= self.get_amount(item,'kg')
        return total_mass
        
    def get_density(self):
        return self.get_mass() / self.get_volume()
        
    def get_viscosity_dynamic(self,temperature=25):
        '''
        Return the dynamic (absolute) viscosity of water in N-s/m2 = Pa-s = kg/m-s
        at the specified temperature.
    
        Parameters:
        ----------
        temperature : float or int, optional
                    The temperature in Celsius. Defaults to 25 degrees if not specified.
        
        Notes:
        ------
        Calculated using the Einstein relation, strictly valid only for dilute solutions[1]:
        
        .. math::
        \eta = \eta_o (1 + 2.5 \sum_i^j \theta_i
        
        \theta_i = {4 \pi r_i ^3 \over 3} {Na C_i \over 1000}
        
        Where $\C_i$ is the molar concentration and $r_i$ is the hydrodynamic radius in meters.
    
        .. [1] Smedley, Stuart. The Interpretation of Ionic Conductivity in Liquids, pp 13-14. Plenum Press, 1980.
        
        '''
        return water_viscosity_dynamic(temperature)
        #TODO
        
    
    def get_viscosity_kinematic(self,temperature=25):
        '''
        Return the kinematic viscosity of water in m2/s = Stokes
        at the specified temperature.
        
        Parameters:
        ----------
        temperature : float or int, optional
                    The temperature in Celsius. Defaults to 25 degrees if not specified.
        
        Returns:
        -------
        float
                The kinematic viscosity of water in Stokes (m2/s)
                
        See Also:
        --------
        get_density_dynamic
        '''
        return self.get_viscosity_dynamic(temperature) / self.get_density(temperature)
        
    def get_osmotic_pressure(self):
        ''' 
        Return the osmotic pressure of the solution relative to pure water
        
        Parameters:
        ----------
        none
        
        Returns:
        -------
        float
                The osmotic pressure of the solution relative to pure water in Pa
                
        References:
        ----------
        Sata, Toshikatsu. /Ion Exchange Membranes: Preparation, Characterization, and Modification./ Royal Society of Chemistry, 2004, p. 10.
        
        http://en.wikipedia.org/wiki/Osmotic_pressure#Derivation_of_osmotic_pressure
        
        Examples:
        --------
        TODO - add examples for osmotic pressure
        
        '''
        # TODO - tie this into parameter() and solvent() objects
        partial_molar_volume_water = 18e-6
        
        osmotic_pressure = const.R * kelvin(self.temperature) / partial_molar_volume_water * math.log (self.get_water_activity())
        logger.info('Computed osmotic pressure of solution as %s Pa at T= %s degrees C' % (osmotic_pressure,self.temperature))
        return osmotic_pressure
    
    def get_conductivity(self):
        return self.conductivity
        
    # to be deprecated TODO
    def get_unit_cost(self):
        print('DEPRECATE!')
        return self.unit_storage_cost


## Concentration  Methods        
    
    def get_amount(self,solute,unit,temperature=25):
        '''returns the amount of 'solute' in the parent solution
       
        Parameters:
        ----------
        solute : str 
                    String representing the name of the solute of interest
        unit : str
                    Units desired for the output. Valid units are 'mol/L','mol/kg','mol','fraction', 'kg', and 'g/L'
        temperature : float or int, optional
                    The temperature in Celsius. Defaults to 25 degrees if not specified.
        
        Returns:
        -------
        The amount of the solute in question, in the specified units
        
        Notes:
        -----
        Amount-per-volume units such as mol/L and g/L are dependent on temperature. Others are not.
        
        See Also:
        --------
        '''
        moles = self.get_solute(solute).get_moles()
        mw = self.get_solute(solute).get_molecular_weight()
        
        # with pint unit conversions enable, we just pass the unit to pint
        return moles.to(unit,'chem',mw=mw,volume=self.get_volume(temperature),solvent_mass=self.get_solvent_mass())
        
#        if unit == 'mol':
#            return moles
#        elif unit == 'mmol':
#            return moles * 1000
#        elif unit == 'mol/kg':
#            return moles / self.get_solvent_mass()
#        elif unit == 'fraction':
#            return moles / (self.get_total_moles_solute() + self.get_moles_water())
#        # mass per volume units
#        elif unit == 'mol/L':
#            return moles / self.get_volume(temperature)
#        elif unit == 'mmol/L':
#            return moles * 1000 / self.get_volume(temperature)
#        elif unit == 'ng/L':
#            return moles * self.components[solute].get_molecular_weight() * 1e9 / self.get_volume(temperature)
#        elif unit == 'ug/L':
#            return moles * self.components[solute].get_molecular_weight() * 1e6 / self.get_volume(temperature)
#        elif unit == 'mg/L':
#            return moles * self.components[solute].get_molecular_weight() * 1e3 / self.get_volume(temperature)
#        elif unit == 'g/L':
#            return moles * self.components[solute].get_molecular_weight() / self.get_volume(temperature)
#        elif unit == 'kg/L':
#            return moles * self.components[solute].get_molecular_weight() / 1e3 / self.get_volume(temperature)
#        # mass units
#        elif unit == 'ng':
#            return moles * self.components[solute].get_molecular_weight() * 1e9
#        elif unit == 'ug':
#            return moles * self.components[solute].get_molecular_weight() * 1e6
#        elif unit == 'mg':
#            return moles * self.components[solute].get_molecular_weight() * 1e3
#        elif unit == 'g':
#            return moles * self.components[solute].get_molecular_weight()
#        elif unit == 'kg':
#            return moles * self.components[solute].get_molecular_weight() / 1e3
#        # mass fraction units
#        elif unit == 'ppt':
#            return moles * self.components[solute].get_molecular_weight() / (self.get_mass * 1e3) * 1e8
#        elif unit == 'ppb':
#            return moles * self.components[solute].get_molecular_weight() / (self.get_mass * 1e3) * 1e7
#        elif unit == 'ppm' or unit == 'mg/kg':
#            return moles * self.components[solute].get_molecular_weight() / (self.get_mass * 1e3) * 1e6
#        elif unit == '%':
#            return moles * self.components[solute].get_molecular_weight() / (self.get_mass * 1e3) * 100
#        elif unit == 'g/g' or unit == 'mg/mg' or unit == 'kg/kg':
#            return moles * self.components[solute].get_molecular_weight() / (self.get_mass * 1e3)
#        else:
#            print('Invalid unit %s specified for amount % unit')
#            return None
    
    def set_amount(self,solute,amount,unit,temperature=25):
        '''returns the amount of 'solute' in the parent solution
       
        Parameters:
        ----------
        solute : str 
                    String representing the name of the solute of interest
        unit : str
                    Units desired for the output. Valid units are 'mol/L','mol/kg','mol','fraction', 'kg', and 'g/L'
        temperature : float or int, optional
                    The temperature in Celsius. Defaults to 25 degrees if not specified.
        
        Returns:
        -------
        The amount of the solute in question, in the specified units
        
        Notes:
        -----
        Amount-per-volume units such as mol/L and g/L are dependent on temperature. Others are not.
        
        See Also:
        --------
        '''
        
        if unit == 'mol':
            self.components[solute].set_moles(amount)
        elif unit == 'mol/L':
            self.components[solute].set_moles(amount * self.get_volume(temperature))
        elif unit == 'mol/kg':
            self.components[solute].set_moles(amount * self.get_solvent_mass())
        elif unit == 'g/L':
            self.components[solute].set_moles(amount / molecular_weight * self.get_volume(temperature))
        elif unit == 'fraction':
            self.components[solute].set_moles(amount * (self.get_total_moles_solute() + self.get_moles_water()))
        elif unit == 'kg':
            self.components[solute].set_moles(amount / self.components[solute].get_molecular_weight() * 1000)
        else:
            print('Invalid unit %s specified for amount' % unit)
            return None
        # to be deprecated TODO
    
    def get_moles(self,solute):
        '''(str) -> float
        Return the total moles of 'solute' in the parent Solution
        
        Parameters:
        ----------
        solute : str 
                 String representing the name of the solute of interest
    
        Returns:
        -------
        float
            The total moles of 'solute' in the parent Solution
    
        See Also:
        --------
        get_solvent_mass()
        
        Examples:
        --------
        TBD
        
        '''
        print('DEPRECATE!')
        return self.get_amount(solute,'mol')
        
    def get_total_moles_solute(self):
        '''Return the total moles of all solute in the solution'''
        tot_mol = 0
        for item in self.components:
            if item != self.solvent_name:
                tot_mol += self.components[item].get_moles()
        return tot_mol
    
    #to be deprecated TODO
    def get_mole_fraction(self,solute):
        '''(Solute) -> float
        Return the mole fraction of 'solute' in the solution
        
        
        Parameters:
        ----------
        solute : str 
                 String representing the name of the solute of interest
    
        Returns:
        -------
        float
            The mole fraction of 'solute' in the parent Solution object
    
        See Also:
        --------
        get_solvent_mass()
        
        Notes:
        -----
        This function assumes water is the solvent with MW = 18
 
        Examples:
        --------
        TBD
        
        '''
        print('DEPRECATE!')
        return self.get_amount(solute,'fraction')
    
    def get_moles_water(self):
        return self.get_amount(self.solvent_name,'mol',self.temperature)
    
    # to be deprecated TODO
    def get_molar_concentration(self,solute):
        '''(Solute) -> float
        Return the molar concentration of 'solute' in the solution
        
        
        Parameters:
        ----------
        solute : str 
                 String representing the name of the solute of interest
    
        Returns:
        -------
        float
            The mole fraction of 'solute' in the parent Solution object
    
        See Also:
        --------
        get_solvent_mass()
        
        Notes:
        -----
        This function assumes water is the solvent with MW = 18
 
        Examples:
        --------
        TBD
        '''
        print('DEPRECATE!')
        return self.get_amount(solute,'mol/L')
    
    
## Activity-related methods
    def get_activity_coefficient(self,solute,temperature=25):
        '''Routine to determine the activity coefficient of a solute in solution. The correct function is chosen based on the ionic strength of the parent solution.
        
        Parameters:
        ----------
        solute : str 
                    String representing the name of the solute of interest
        temperature : float or int, optional
                    The temperature in Celsius. Defaults to 25 degrees if not specified.
        
        Returns:
        -------
        The molar (mol/L) scale mean ion activity coefficient of the solute in question
        
        See Also:
        --------
        get_activity_coefficient_debyehuckel
        get_activity_coefficient_guntelberg
        get_activity_coefficient_davies
        get_activity_coefficient_pitzer
        get_activity_coefficient_TCPC
        '''
        ion = self.components[solute]
        # for very low ionic strength, use the Debye-Huckel limiting law
        
        if self.get_ionic_strength() <= 0.005:
            logger.info('Ionic strength = %s. Using Debye-Huckel to calculate activity coefficient.' % self.get_ionic_strength())
            return get_activity_coefficient_debyehuckel(self.get_ionic_strength(),temperature)
            
        # use the Guntelberg approximation for 0.005 < I < 0.1
        elif self.get_ionic_strength() <= 0.1:
            logger.info('Ionic strength = %s. Using Guntelberg to calculate activity coefficient.' % self.get_ionic_strength())
            return get_activity_coefficient_guntelberg(self.get_ionic_strength(),ion.get_formal_charge(),temperature)
            
        # use the Davies equation for 0.1 < I < 0.5
        elif self.get_ionic_strength() <= 0.5:
            logger.info('Ionic strength = %s. Using Davies equation to calculate activity coefficient.' % self.get_ionic_strength())
            return get_activity_coefficient_davies(self.get_ionic_strength(),ion.get_formal_charge(),temperature)
            
        # use the TCPC model for higher ionic strengths, if the parameters have been set
        elif self.components[solute].parameters_TCPC:
            logger.info('Ionic strength = %s. Using TCPC model to calculate activity coefficient.' % self.get_ionic_strength())
            return get_activity_coefficient_TCPC(self.get_ionic_strength(),ion.get_parameters_TCPC('S'),ion.get_parameters_TCPC('b'),ion.get_parameters_TCPC('n'),ion.get_parameters_TCPC('z_plus'),ion.get_parameters_TCPC('z_minus'),ion.get_parameters_TCPC('nu_plus'),ion.get_parameters_TCPC('nu_minus'),temperature)
            
        else:
            print('WARNING: Ionic strength too high to estimate activity. Specify parameters for Pitzer or TCPC methods. Returning unit activity coefficient')
            return 1.0
        # TODO - NEED TO TEST THIS FUNCTION
    
    def get_activity(self,solute,temperature=25):
        '''returns the thermodynamic activity of the solute in solution
       
        Parameters:
        ----------
        solute : str 
                    String representing the name of the solute of interest
        temperature : float or int, optional
                    The temperature in Celsius. Defaults to 25 degrees if not specified.
        
        Returns:
        -------
        The thermodynamic activity of the solute in question (dimensionless)
        
        Notes:
        -----
        The thermodynamic activity is independent of the concentration scale used. However,
        the concentration and the activity coefficient must use corresponding scales.[1][2]
        In this module, ionic strength, activity coefficients, and activities are all
        calculated based on the molal (mol/kg) concentration scale.
        
        References:
        ----------
        ..[1] http://adsorption.org/awm/utils/Activity.htm
        ..[2] http://en.wikipedia.org/wiki/Thermodynamic_activity#Activity_coefficient
        
        See Also:
        --------
        get_activity_coefficient
        get_ionic_strength
        
        '''
        # switch to the water activity function if the species is H2O
        if solute == 'H2O' or solute == 'water':
            # find the predominant non-solvent solute
            most = 0
            predominant_solute = ''
            for item in self.components:
                mass = self.get_amount(item,'mol') 
                if item != self.solvent_name and mass > most:
                    most = mass
                    predominant_solute = item
            if most > 0:       
                activity = self.get_water_activity(predominant_solute,temperature)
                logger.info('Calculated activity of solvent (water) as %s using osmotic coefficient based on solute %s.' % (activity,predominant_solute))
                
            # return 1.0 water activity if there are no solutes
            else: 
                activity = 1
                logger.info('Calculated activity of solvent (water) as 1.0 because no solutes were found.')
# TODO fix this for multivalent salts e.g. MgCl2
        else:
            activity = self.get_activity_coefficient(solute,temperature) * self.get_amount(solute,'mol/kg').magnitude
            logger.info('Calculated activity of solute %s as %s' % (solute,activity))
        
        return activity

    def get_osmotic_coefficient(self,solute,temperature=25):
        '''calculate the osmotic coefficient for a given solute
        
        Parameters:
        ----------
        solute : str 
                    String representing the name of the solute of interest
        temperature : float or int, optional
                    The temperature in Celsius. Defaults to 25 degrees if not specified.
        
        Returns:
        -------
        float : 
            The practical osmotic coefficient, based on 'solute'
        '''
        
        ion = self.components[solute]
        
        if self.components[solute].parameters_TCPC:
            osmotic_coefficient= get_osmotic_coefficient_TCPC(self.get_ionic_strength(),ion.get_parameters_TCPC('S'),ion.get_parameters_TCPC('b'),ion.get_parameters_TCPC('n'),ion.get_parameters_TCPC('z_plus'),ion.get_parameters_TCPC('z_minus'),ion.get_parameters_TCPC('nu_plus'),ion.get_parameters_TCPC('nu_minus'),temperature)
            logger.info('Calculated osmotic coefficient of water as %s based on solute %s using TCPC model' % (osmotic_coefficient,solute))
            return osmotic_coefficient
        else:
            logger.error('Cannot calculate water activity because TCPC parameters for solute are not specified. Returning unit osmotic coefficient')
            return 1
        
    def get_water_activity(self,solute,temperature=25):
        '''return the water activity based on a given solute
        
        Parameters:
        ----------
        solute : str 
                    String representing the name of the solute of interest
        temperature : float or int, optional
                    The temperature in Celsius. Defaults to 25 degrees if not specified.
        
        Returns:
        -------
        float : 
            The thermodynamic activity of water in the solution.
        
        Notes:
        -----
        Water activity is related to the osmotic coefficient in a solution containing i solutes by:[1]
        
        ## ln a_w = - \Phi M_w \sum_i m_i
        
        Where M_w is the molar mass of water (0.018015 kg/mol)
        
        References:
        ----------
        Blandamer, Mike J., Engberts, Jan B. F. N., Gleeson, Peter T., Reis, Joao Carlos R., 2005. "Activity of water in aqueous systems: A frequently neglected property."
        //Chemical Society Review// 34, 440-458.
        
        '''
        return math.exp(- self.get_osmotic_coefficient(solute,temperature) * 0.018015 * self.get_total_moles_solute())
        
    def get_ionic_strength(self):
        '''() -> float
        
        Return the ionic strength of the solution, calculated as 1/2 * sum ( molality * charge ^2) over all the ions.
        Molal (mol/kg) scale concentrations are used for compatibility with the activity correction formulas, but
        the returned value does not carry units.
        
        
        Returns:
        -------
        float : 
            The molal scale ionic strength of the parent solution.
        
        Examples:
        --------
        TODO
#         >>> conc_soln.list_concentrations()
#         {'Na+': 5.999375074924214, 'Cl-': 5.999904143046362, 'HCO3-': 0, 'NaCO3-': 0, 'NaHCO3': 0}
#         >>> conc_soln.get_ionic_strength()
#         5.999639608985288
        '''
        self.ionic_strength=0
        for solute in self.components.keys():
            self.ionic_strength += 0.5 * self.get_amount(solute,'mol/kg') * self.components[solute].get_formal_charge() ** 2
       
        return self.ionic_strength.magnitude
            
    ## informational methods
    def list_solutes(self):
        return list(self.components.keys())
    
    # TODO - deprecate this
    def list_concentrations(self,unit):
        print('DEPRECATE!')
        list_components(self,unit)
    
    def list_components(self,unit):
        '''() -> dict
        
        Return a dictionary containing a list of the species in solution paired with their amount in the specified units
        '''
        self.mol_list={}
        for i in self.components.keys():
            self.mol_list.update({i:self.get_amount(i,unit)})
        print('Component amounts (%s):' % unit,self.mol_list )
        
    def list_activities(self):
        '''() -> dict
        
        Return a dictionary containing a list of the species in solution paired with their molal activity
        '''
        self.act_list={}
        for i in self.components.keys():
            self.act_list.update({i:self.components[i].get_activity()})
        print('Component activities:',self.act_list )
    
    # TODO deprecate this
    def list_mole_fractions(self):
        '''() -> dict
        
        Return a dictionary containing a list of the species in solution paired with their mole fraction
        '''
        self.fraction_list={}
        for i in self.components.keys():
            self.fraction_list.update({i:self.get_mole_fraction(i)})
        # add mole fraction for water
        self.fraction_list.update({'H2O':self.get_moles_water()/(self.get_moles_water() + self.get_total_moles_solute())})
        print('DEPRECATE!')
        return self.fraction_list
        
   
    def __str__(self):
        #set output of the print() statement for the solution     
        return 'Components: '+str(self.list_solutes()) + '\n' + 'Volume: '+str(self.get_volume()) + '\n' + 'Density: '+str(self.get_density())

class Membrane:
    '''Class representing the properties of various kinds of water treatment membranes'''
    
    
    def __init__(self,name,type,permselectivity,area_resistance,cost,thickness,fixed_charge_density):
        '''(str,str,float,float,float,number,number) -> Membrane
        
        name is a str describing the membrane
        type indicates the kind of membrane. Valid types are 'aem' 'cem' 'bpem' 'mf' 'uf' and 'ro'
        perm is a number representing the membrane permselectivity (0 < perm < 1)
        resist is the areal resistance of the membrane, in ohm-m2
        cost is the unit cost of the membrane, in $/m2
        thickness is the thickness of the membrane in m
        fixed_charge_density is the concentration of fixed charges (for IX membranes), eq/m3
        
        '''
        #warn if an invalid membrane type is given
        types = ['aem','cem','bpem','mf','uf','ro','fo']
        if type in types:
            self.mem_type = type
        else:
            self.mem_type = 'Invalid'
            print("ERROR: Invalid membrane type specified.")
       
        self.title = name
        self.permselectivity = permselectivity
        self.resistance = area_resistance
        self.unit_cost = cost
        self.thickness = thickness
        self.fixed_charge_density = fixed_charge_density
    
    #simple methods to access the main properties
    def get_mem_type(self):
        return self.mem_type
    def get_permselectivity(self):
        return self.permselectivity
    def get_resistance(self):
        return self.resistance
    def get_unit_cost(self):
        return self.unit_cost
    def get_thickness(self):
        return self.thickness
    def get_fixed_charge_density(self):
        return self.fixed_charge_density
        
    #set output of the print() statement for the solution     
    def __str__(self):
        return self.title + ' -- Type: ' + self.mem_type + '  Permselectivity: ' +str(round(self.permselectivity,3)) + ' Resistance: ' + str(self.resistance) + ' ohm-m2  Cost: ' +  str(self.unit_cost)+' $/m2'
        
# TODO - turn doctest back on when the nosigint error is gone        
#if __name__ == "__main__":
 #   import doctest
  #  doctest.testmod()