'''
pyEQL water properties library

This file contains functions for retrieving various physical properties
of water substance

'''
import math

from parameter import unit

import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def water_density(temperature=25*unit('degC'),pressure=1*unit('atm')):
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
    997.0415 kilogram/meter3 
    
    '''
    # calculate the magnitude
    density = 999.65 + 0.20438 * temperature.to('degC').magnitude - 6.1744e-2 * temperature.to('degC').magnitude ** 1.5
    # assign the proper units
    density = density  * unit('kg/m**3')
    logger.info('Computed density of water as %s at T= %s and P = %s' % (density,temperature,pressure))
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
    

def water_dielectric_constant(temperature=25*unit('degC')):
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
    if temperature.to('K').magnitude < 273 or temperature.to('K').magnitude > 372:
        logger.error('Specified temperature (%s K) exceeds valid range of data. Cannot extrapolate. % kelvin(temperature)')
        return None
    
    # otherwise, calculate the dielectric constant using the quadratic fit    
    a = 0.24921e3
    b = -0.79069e0
    c = 0.72997e-3
    dielectric = a + b * temperature.to('K').magnitude + c * temperature.to('K').magnitude ** 2
    
    logger.info('Computed dielectric constant of water as %s at %s degrees Celsius' % (dielectric,temperature))
    
    logger.debug('Computed dielectric constant of water using empirical equation given in "Permittivity (Dielectric Constant) of Liquids." CRC Handbook of Chemistry and Physics, 92nd ed, pp 6-187 - 6-208.')
    
    return dielectric
    
    
def water_conductivity(temperature):
    pass

