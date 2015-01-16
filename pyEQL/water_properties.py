'''
pyEQL water properties library

This file contains functions for retrieving various physical properties
of water substance

'''
import math

from pyEQL.parameter import unit
# TODO fix this to handle offsets the way pint wants us to since 0.7
unit.autoconvert_offset_to_baseunit = True

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
    
    
    .. math:: \\rho_W = 999.65 + 0.20438 T - 6.1744e-2 T ^ {1.5}
    
    Where T is the temperature in Celsius.
    
    
    ..[1] Sohnel, O and Novotny, P. //Densities of Aqueous Solutions of Inorganic Substances.// Elsevier Science, Amsterdam, 1985.
    
    Examples:
    --------
    >>> water_density(25*unit('degC')) #doctest: +ELLIPSIS
    <Quantity(997.0415, 'kilogram / meter ** 3')>
    
    '''
    # calculate the magnitude
    density = 999.65 + 0.20438 * temperature.to('degC').magnitude - 6.1744e-2 * temperature.to('degC').magnitude ** 1.5
    # assign the proper units
    density = density  * unit('kg/m**3')
    logger.info('Computed density of water as %s at T= %s and P = %s' % (density,temperature,pressure))
    logger.debug('Computed density of water using empirical relation in Sohnel and Novotny, "Densities of Aqueous Solutions of Inorganic Substances," 1985' )
    return density.to('kg/m**3')
    
def water_specific_weight(temperature=25*unit('degC'),pressure=1*unit('atm')):
    '''(number) -> float
    
    Return the specific weight of water in N/m3 at the specified temperature and pressure.
    
    Parameters:
    ----------
    temperature : Quantity, optional
                  The temperature. Defaults to 25 degC if omitted.
    pressure    : Quantity, optional
                  The ambient pressure of the solution. 
                  Defaults to atmospheric pressure (1 atm) if omitted.
                  
    Returns:
    -------
    Quantity
            The specific weight of water in N/m3.  
    
    Examples:
    --------
    >>> water_specific_weight() #doctest: +ELLIPSIS
    <Quantity(9777.637025975, 'newton / meter ** 3')>
            
    See Also:
    --------
    water_density
    
    '''
    spweight = water_density(temperature,pressure) * unit.g_n
    logger.info('Computed specific weight of water as %s at T=%s and P = %s' % (spweight,temperature,pressure))
    return spweight.to('N/m ** 3')

def water_viscosity_dynamic(temperature=25*unit('degC'),pressure=1*unit('atm')):
    '''
    Return the dynamic (absolute) viscosity of water in N-s/m2 = Pa-s = kg/m-s
    at the specified temperature.
    
    Parameters:
    ----------
    temperature : Quantity, optional
                  The temperature. Defaults to 25 degC if omitted.
    pressure    : Quantity, optional
                  The ambient pressure of the solution. 
                  Defaults to atmospheric pressure (1 atm) if omitted.
    
    Returns:
    -------
    Quantity 
                The dynamic (absolute) viscosity of water in N-s/m2 = Pa-s = kg/m-s
                  
    Notes:
    -----
    Implements the international equation for viscosity of water as specified by NIST[1]
    
    Valid for 273 < temperature < 1073 K and 0 < pressure < 100,000,000 Pa
    
    .. [1] Sengers, J.V. "Representative Equations for the Viscosity of Water Substance." 
        J. Phys. Chem. Ref. Data 13(1), 1984.http://www.nist.gov/data/PDFfiles/jpcrd243.pdf
    
    Examples:
    --------
    >>> water_viscosity_dynamic(20*unit('degC')) #doctest: +ELLIPSIS
    <Quantity(0.000998588610804179, 'kilogram / meter / second')>
    >>> water_viscosity_dynamic(unit('100 degC'),unit('25 MPa')) #doctest: +ELLIPSIS
    <Quantity(0.00028165034364318573, 'kilogram / meter / second')>
    >>> water_viscosity_dynamic(25*unit('degC'),0.1*unit('MPa')) #doctest: +ELLIPSIS
    <Quantity(0.0008872817880143659, 'kilogram / meter / second')>
    
    #TODO - check these again after I implement pressure-dependent density function
    
    '''
    # generate warnings if temp or pressure are outside valid range of equation
    if temperature < 273 * unit('K') or temperature > 1073 * unit('K'):
        logger.error('Specified temperature (%s) exceeds valid range of NIST equation for viscosity of water. Cannot extrapolate.' % temperature)
        return None
        
    if pressure < 0 * unit('Pa') or pressure > 100000000 * unit ('Pa'):
        logger.error('Specified pressure (%s) exceeds valid range of NIST equation for viscosity of water. Cannot extrapolate.' % pressure)
        return None
    
    # calculate dimensionless temperature and pressure
    T_star = 647.27*unit('K') #K
    P_star = 22115000*unit('Pa') #Pa
    rho_star = 317.763*unit('kg/m**3') #kg/m3
    
    T_bar = temperature.to('K') / T_star
    P_bar = pressure.to('Pa') / P_star
    rho_bar = water_density(temperature,pressure) / rho_star
    
    # calculate the first function, mu_o
    mu_star = 1e-6*unit('Pa * s') #Pa-s
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
    
    logger.info('Computed dynamic (absolute) viscosity of water as %s at T=%s and P = %s'  % (viscosity,temperature,pressure)) 
    
    logger.debug('Computed dynamic (absolute) viscosity of water using empirical NIST equation described in Sengers, J.V. "Representative Equations for the Viscosity of Water Substance." J. Phys. Chem. Ref. Data 13(1), 1984.')
    
    return viscosity.to('kg/m/s')


def water_viscosity_kinematic(temperature=25*unit('degC'),pressure=1*unit('atm')):
    '''
    Return the kinematic viscosity of water in m2/s = Stokes
    at the specified temperature.
    
    Parameters:
    ----------
    temperature : Quantity, optional
                  The temperature. Defaults to 25 degC if omitted.
    pressure    : Quantity, optional
                  The ambient pressure of the solution. 
                  Defaults to atmospheric pressure (1 atm) if omitted.
                  
    Returns:
    -------
    Quantity
            The kinematic viscosity of water in Stokes (m2/s)
    
    Examples:
    --------
    >>> water_viscosity_kinematic()  #doctest: +ELLIPSIS
    <Quantity(8.899146003595295e-07, 'meter ** 2 / second')>
            
    See Also:
    --------
    water_viscosity_dynamic
    water_density
    
    '''
    kviscosity = water_viscosity_dynamic(temperature,pressure) / water_density(temperature,pressure)
    logger.info('Computed kinematic viscosity of water as %s at T=%s and P = %s ' % (kviscosity,temperature,pressure)) 
    return kviscosity.to('m**2 / s')
    

def water_dielectric_constant(temperature=25*unit('degC')):
    '''(number) -> float
    
    Return the dielectric constant of water at the specified temperature.
    
    Parameters:
    ----------
    temperature : Quantity, optional
                  The temperature. Defaults to 25 degC if omitted.
                  
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
    
    .. math:: \\epsilon(T) = a + b T + c T^2
    
    .. [1] "Permittivity (Dielectric Constant) of Liquids." CRC Handbook of 
            Chemistry and Physics, 92nd ed, pp 6-187 - 6-208.
    
    Examples:
    --------
    >>> water_dielectric_constant(unit('20 degC')) #doctest: +ELLIPSIS
    80.15060...
    
    Display an error if 'temperature' is outside the valid range
    
    >>> water_dielectric_constant(-5*unit('degC'))
    
     
    '''
    # do not return anything if 'temperature' is outside the range for which
    # this fit applies
    if temperature < 273 * unit('K') or temperature > 372 * unit('K'):
        logger.error('Specified temperature (%s) exceeds valid range of data. Cannot extrapolate.' % temperature.to('K'))
        return None
    
    # otherwise, calculate the dielectric constant using the quadratic fit    
    a = 0.24921e3
    b = -0.79069e0
    c = 0.72997e-3
    dielectric = a + b * temperature.to('K').magnitude + c * temperature.to('K').magnitude ** 2
    
    logger.info('Computed dielectric constant of water as %s at %s' % (dielectric,temperature))
    
    logger.debug('Computed dielectric constant of water using empirical equation given in "Permittivity (Dielectric Constant) of Liquids." CRC Handbook of Chemistry and Physics, 92nd ed, pp 6-187 - 6-208.')
    
    return dielectric
    
def water_conductivity(temperature):
    pass

# TODO - turn doctest back on when the nosigint error is gone        
if __name__ == "__main__":
    import doctest
    doctest.testmod()