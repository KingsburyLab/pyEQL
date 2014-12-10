'''
pyEQL activity correction library

This file contains functions for computing molal-scale activity coefficients 
of ions and salts in aqueous solution.

Individual functions for activity coefficients are defined here so that they 
can be used independently of a pyEQL solution object. Normally, these functions
are called from within the get_activity_coefficient method of the Solution class.

'''
import math

# functions for properties of water
import pyEQL.water_properties as h2o

# the pint unit registry
from pyEQL.parameter import unit
# TODO fix this to handle offsets the way pint wants us to since 0.7
unit.autoconvert_offset_to_baseunit = True

# logging system
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def _debye_parameter_B(temperature='25 degC'):
    '''
    return the constant B used in the extended Debye-Huckel equation
    
    Parameters:
    ----------
    temperature : str Quantity, optional
                  String representing the temperature of the solution. Defaults to '25 degC' if not specified.
    
    Notes:
    -----
    The parameter B is equal to:[1]

    .. math:: B = ( {8 \\pi N_A e^2 \\over 1000 \\epsilon k T} ) ^ {1 \\over 2}
    
    .. [1] Bockris and Reddy. /Modern Electrochemistry/, vol 1. Plenum/Rosetta, 1977, p.210.    
    
    Examples:
    --------
    >>> _debye_parameter_B() #doctest: +ELLIPSIS
    0.3291...
    
    '''
    # TODO - fix this and resolve units
    param_B = ( 8 * math.pi * unit.avogadro_number * unit.elementary_charge ** 2 
    / (h2o.water_density(unit(temperature)) * unit.epsilon_0 * h2o.water_dielectric_constant(unit(temperature)) * unit.boltzmann_constant * unit(temperature)) )** 0.5
    return param_B.to_base_units()
    
def _debye_parameter_activity(temperature='25 degC'):
    '''(number) -> float
    return the constant A for use in the Debye-Huckel limiting law (base 10)
    
    Parameters:
    ----------
    temperature : str Quantity, optional
                  String representing the temperature of the solution. Defaults to '25 degC' if not specified.
    
    Returns:
    -------
    Quantity          The parameter A for use in the Debye-Huckel limiting law (base e)
    
    Notes:
    -----
     
    The parameter A is equal to:[1]
     
    ..  math::    
        A^{\\gamma} = {e^3 ( 2 \\pi N_A {\\rho})^{0.5} \\over (4 \\pi \\epsilon_o \\epsilon_r k T)^{1.5}}
    
    Note that this equation returns the parameter value that can be used to calculate
    the natural logarithm of the activity coefficient. For base 10, divide the
    value returned by 2.303. The value is often given in base 10 terms (0.509 at
    25 degC) in older textbooks.
    
    .. [1] Archer, Donald G. and Wang, Peiming. "The Dielectric Constant of Water \
    and Debye-Huckel Limiting Law Slopes." /J. Phys. Chem. Ref. Data/ 19(2), 1990.
        
    Examples:
    --------
    >>> _debye_parameter_activity() #doctest: +ELLIPSIS
    1.17499...
    
    See Also:
    --------
    _debye_parameter_osmotic
    
    '''
    
    debyeparam = unit.elementary_charge ** 3 * ( 2 * math.pi * unit.avogadro_number * h2o.water_density(unit(temperature))) ** 0.5 \
    / ( 4 * math.pi * unit.epsilon_0 * h2o.water_dielectric_constant(unit(temperature)) * unit.boltzmann_constant * unit(temperature)) ** 1.5
    
    logger.info('Computed Debye-Huckel Limiting Law Constant A = %s at %s' % (debyeparam,temperature))
    return debyeparam.to('kg ** 0.5 / mol ** 0.5')

def _debye_parameter_osmotic(temperature='25 degC'):
    '''(number) -> float
    return the constant A_phi for use in calculating the osmotic coefficient according to Debye-Huckel theory
    
    Parameters:
    ----------
    temperature : str Quantity, optional
                  String representing the temperature of the solution. Defaults to '25 degC' if not specified.
    
    Notes:
    -----
    Not to be confused with the Debye-Huckel constant used for activity coefficients in the limiting law. 
    Takes the value 0.392 at 25 C.
    This constant is calculated according to:[1][2]

     .. math:: A^{\\phi} = {1 \\over 3} A^{\\gamma}
    
    .. [1] Kim, Hee-Talk and Frederick, William Jr, 1988. "Evaluation of Pitzer Ion Interaction Parameters of Aqueous Electrolytes at 25 C. 1. Single Salt Parameters,"
    //J. Chemical Engineering Data// 33, pp.177-184.
    
    .. [2] Archer, Donald G. and Wang, Peiming. "The Dielectric Constant of Water \
    and Debye-Huckel Limiting Law Slopes." /J. Phys. Chem. Ref. Data/ 19(2), 1990.
        
    
    Examples:
    --------
    >>> _debye_parameter_osmotic() #doctest: +ELLIPSIS
    0.3916... 
  
    
    See Also:
    --------
    _debye_parameter_activity
    
    '''
    
    output = 1/3 * _debye_parameter_activity(temperature)
    return output.to('kg ** 0.5 /mol ** 0.5')

def _debye_parameter_volume(temperature='25 degC'):
    '''
    return the constant A_V, the Debye-Huckel limiting slope for apparent
    molar volume.
    
    Parameters:
    ----------
    temperature : str Quantity, optional
                  String representing the temperature of the solution. Defaults to '25 degC' if not specified.
    
    Notes:
    -----
    Takes the value 1.898 cm **3 * kg ** 0.5 /  mol ** 1.5 at 25 C.
    This constant is calculated according to:[1]

     .. math:: A_V = -2 A_{\\phi} R T [ {3 \\over \\epsilon} {{\\partial \\epsilon \\over \\partial p} \
     } - {{1 \\over \\rho}{\\partial \\rho \\over \\partial p} }]
    
    .. [1] Archer, Donald G. and Wang, Peiming. "The Dielectric Constant of Water \
    and Debye-Huckel Limiting Law Slopes." /J. Phys. Chem. Ref. Data/ 19(2), 1990.
        
    
    Examples:
    --------
    TODO
  
    
    See Also:
    --------
    _debye_parameter_osmotic
    '''
    
    # TODO - add partial derivatives to calculation
    #result = -2 * _debye_parameter_osmotic(temperature) * unit.R * temperature * (unit('1 * Pa ** -1'))
    result = unit('1.898 cm ** 3 * kg ** 0.5 /  mol ** 1.5')
    
    if unit(temperature) != unit('25 degC'):
        logger.warning('Debye-Huckel limiting slope for volume is valid only at 25 degC')
        
    return result.to('cm ** 3 * kg ** 0.5 /  mol ** 1.5')


def get_activity_coefficient_debyehuckel(ionic_strength,formal_charge=1,temperature='25 degC'):
    '''Return the activity coefficient of solute in the parent solution according to the Debye-Huckel limiting law.
    
    Parameters
    ----------
    formal_charge : int, optional      
                    The charge on the solute, including sign. Defaults to +1 if not specified.
    ionic_strength : Quantity
                     The ionic strength of the parent solution, mol/kg
    temperature :    str Quantity, optional
                     String representing the temperature of the solution. Defaults to '25 degC' if not specified.
                  
    Returns
    -------
    float
         The mean molal (mol/kg) scale ionic activity coefficient of solute

    See Also:
    --------
    _debye_parameter_activity
    
    Notes
    -----
    Activity coefficient is calculated according to: 

    .. math:: \\ln \\gamma = A^{\\gamma} z_i^2 \sqrt I
    
    Valid only for I < 0.005
    
    .. [1] Stumm, Werner and Morgan, James J. Aquatic Chemistry, 3rd ed, 
    pp 103. Wiley Interscience, 1996.
    '''
    # check if this method is valid for the given ionic strength
    if not ionic_strength.magnitude <= 0.005:
        logger.warning('Ionic strength exceeds valid range of the Debye-Huckel limiting law')
    
    log_f = - _debye_parameter_activity(temperature) * formal_charge ** 2 * ionic_strength ** 0.5

    return math.exp(log_f)

def get_activity_coefficient_guntelberg(ionic_strength,formal_charge=1,temperature='25 degC'):
    '''Return the activity coefficient of solute in the parent solution according to the Guntelberg approximation.
    
    Parameters
    ----------
    formal_charge : int, optional      
                    The charge on the solute, including sign. Defaults to +1 if not specified.
    ionic_strength : Quantity
                     The ionic strength of the parent solution, mol/kg
    temperature :    str Quantity, optional
                     String representing the temperature of the solution. Defaults to '25 degC' if not specified.
                  
    Returns
    -------
    float
         The mean molal (mol/kg) scale ionic activity coefficient of solute
         
    See Also
    --------
    _debye_parameter_activity
    
    Notes
    ------
    Activity coefficient is calculated according to: 

    .. math:: \\ln \\gamma = A^{\\gamma} z_i^2 {\sqrt I \\over (1 + \sqrt I)}
    
    Valid for I < 0.1
    
    .. [1] Stumm, Werner and Morgan, James J. Aquatic Chemistry, 3rd ed, 
    pp 103. Wiley Interscience, 1996.
    '''
    # check if this method is valid for the given ionic strength
    if not ionic_strength.magnitude <= 0.1:
        logger.warning('Ionic strength exceeds valid range of the Guntelberg approximation')
    
    log_f = - _debye_parameter_activity(temperature) * formal_charge ** 2 * ionic_strength ** 0.5 / (1+ionic_strength.magnitude ** 0.5)

    return math.exp(log_f)
    
def get_activity_coefficient_davies(ionic_strength,formal_charge=1,temperature='25 degC'):
    '''Return the activity coefficient of solute in the parent solution according to the Davies equation.
    
    Parameters
    ----------
    formal_charge : int, optional      
                    The charge on the solute, including sign. Defaults to +1 if not specified.
    ionic_strength : Quantity
                     The ionic strength of the parent solution, mol/kg
    temperature :    str Quantity, optional
                     String representing the temperature of the solution. Defaults to '25 degC' if not specified.
                  
    Returns
    -------
    float
         The mean molal (mol/kg) scale ionic activity coefficient of solute

    See Also
    --------
    _debye_parameter_activity
    
    Notes
    -----
    Activity coefficient is calculated according to: 

    .. math:: \\ln \\gamma = A^{\\gamma} z_i^2 ({\sqrt I \\over (1 + \sqrt I)} + 0.2 I)
    
    Valid for 0.1 < I < 0.5
    
    .. [1] Stumm, Werner and Morgan, James J. Aquatic Chemistry, 3rd ed, 
    pp 103. Wiley Interscience, 1996.
    '''
    # check if this method is valid for the given ionic strength
    if not ionic_strength.magnitude <= 0.5 and ionic_strength.magnitude >= 0.1:
        logger.warning('Ionic strength exceeds valid range of the Davies equation')
    
    # the units in this empirical equation don't work out, so we must use magnitudes
    log_f = - _debye_parameter_activity(temperature).magnitude * formal_charge ** 2 * (ionic_strength.magnitude ** 0.5 / (1+ionic_strength.magnitude ** 0.5) - 0.2 * ionic_strength.magnitude)
    
    return math.exp(log_f)
    
def get_activity_coefficient_TCPC(ionic_strength,S,b,n,formal_charge=1,counter_formal_charge=-1,stoich_coeff=1,counter_stoich_coeff=1,temperature=25*unit('degC')):
    '''Return the activity coefficient of solute in the parent solution according to the modified TCPC model.
    
    Parameters:
    ----------
    ionic_strength : Quantity
                        The ionic strength of the parent solution, mol/kg
    S : float
                        The solvation parameter for the parent salt. See Reference.
    b : float
                        The approaching parameter for the parent salt. See Reference.
    n : float       
                        The n parameter for the parent salt. See Reference.
    formal_charge : int, optional           
                        The charge on the solute, including sign. Defaults to +1 if not specified.
    counter_formal_charge : int, optional           
                        The charge on the solute's complementary ion, including sign. Defaults to -1 if not specified.
                        E.g. if the solute is Na+ and the salt is NaCl, counter_formal_charge = -1
    stoich_coeff : int, optional
                        The stoichiometric coefficient of the solute in its parent salt. Defaults to1 if not specified.
                        E.g. for Zn+2 in ZnCl2, stoich_coeff = 1
    counter_stoich_coeff : int, optional
                        The stoichiometric coefficient of the solute's complentary ion in its parent salt. Defaults to 1 if not specified.
                        E.g. for Cl- in ZnCl2, stoich_coeff = 2
    temperature :    Quantity, optional
                     The temperature of the solution. Defaults to 25 degrees C if omitted
    
    Returns:
    -------
    float
        The mean molal (mol/kgL) scale ionic activity coefficient of solute

    See Also:
    --------
    _debye_parameter_osmotic
    
    Notes:
    ------
    Valid for concentrated solutions up to saturation. Accuracy compares well with the Pitzer approach. See Reference [1] for a compilation of the appropriate parameters for a variety of commonly-encountered electrolytes.
    
    .. [1] Ge, Xinlei, Wang, Xidong, Zhang, Mei, and Seetharaman, Seshadri. "Correlation and Prediction of Activity and Osmotic Coefficients of Aqueous Electrolytes at 298.15 K by the Modified TCPC Model." J. Chemical Engineering Data 52, pp.538-547, 2007.
    '''
    # compute the PDF parameter
    PDH = - math.fabs(formal_charge * counter_formal_charge) * _debye_parameter_osmotic(temperature) * ( ionic_strength ** 0.5 / (1 + b * ionic_strength ** 0.5) + 2/b * math.log(1 + b * ionic_strength ** 0.5))
    # compute the SV parameter
    SV = S / temperature.to('K') * ionic_strength  ** (2*n) / (stoich_coeff + counter_stoich_coeff)
    # add and exponentiate to eliminate the log
    return math.exp(PDH + SV)
    
def get_activity_coefficient_pitzer(ionic_strength,molality,alpha1,alpha2,beta0,beta1,beta2,C_phi,z_cation,z_anion,nu_cation,nu_anion,temperature='25 degC',b=1.2):
    '''Return the activity coefficient of solute in the parent solution according to the Pitzer model.
    
    Parameters:
    ----------
    ionic_strength: Quantity
                    The ionic strength of the parent solution, mol/kg
    molality:       Quantity
                    The molal concentration of the parent salt, mol/kg
    alpha1, alpha2: number
                    Coefficients for the Pitzer model. This function assigns the coefficients
                    proper units of kg ** 0.5 / mol ** 0.5 after they are entered.
    beta0, beta1, beta2, C_phi: number
                    Coefficients for the Pitzer model. These ion-interaction parameters are
                    specific to each salt system.
    z_cation, z_anion: int
                    The formal charge on the cation and anion, respectively
    nu_cation, nu_anion: int
                    The stoichiometric coefficient of the cation and anion in the salt
    temperature:    str Quantity
                    String representing the temperature of the solution. Defaults to '25 degC' if not specified.
    b:              number, optional
                    Coefficient. Usually set equal to 1.2 and 
                    considered independent of temperature and pressure. If provided, this
                    coefficient is assigned proper units of kg ** 0.5 / mol ** 0.5  after
                    entry.
    
    Returns:
    -------
    float
        The mean molal (mol/kg) scale ionic activity coefficient of solute, dimensionless
    
    Examples:
    --------  
    >>> get_activity_coefficient_pitzer(0.5*unit('mol/kg'),0.5*unit('mol/kg'),1,0.5,-.0181191983,-.4625822071,.4682,.000246063,1,-1,1,1,b=1.2)
    ￼0.61915...   
        
    >>> get_activity_coefficient_pitzer(5.6153*unit('mol/kg'),5.6153*unit('mol/kg'),3,0.5,0.0369993,0.354664,0.0997513,-0.00171868,1,-1,1,1,b=1.2)
    ￼0.76331...
    
    NOTE: the examples below are for comparison with experimental and modeling data presented in
    the May et al reference below. 
    
    10 mol/kg ammonium nitrate. Estimated result (from graph) = 0.2725    
    >>> get_activity_coefficient_pitzer(10*unit('mol/kg'),10*unit('mol/kg'),2,0,-0.01709,0.09198,0,0.000419,1,-1,1,1,b=1.2)
    0.22595 ...
    
    5 mol/kg ammonium nitrate. Estimated result (from graph) = 0.3011
    >>> get_activity_coefficient_pitzer(5*unit('mol/kg'),5*unit('mol/kg'),2,0,-0.01709,0.09198,0,0.000419,1,-1,1,1,b=1.2)
    0.30249 ...
    
    18 mol/kg ammonium nitrate. Estimated result (from graph) = 0.1653
    >>> get_activity_coefficient_pitzer(18*unit('mol/kg'),18*unit('mol/kg'),2,0,-0.01709,0.09198,0,0.000419,1,-1,1,1,b=1.2)
    0.16241 ...
    
    
    References:
    ----------
    Scharge, T., Munoz, A.G., and Moog, H.C. (2012). Activity Coefficients of Fission Products in Highly
    Salinary Solutions of Na+, K+, Mg2+, Ca2+, Cl-, and SO42- : Cs+.
    /Journal of Chemical& Engineering Data (57), p. 1637-1647.
    
    Kim, H., & Jr, W. F. (1988). Evaluation of Pitzer ion interaction parameters of aqueous electrolytes at 25 degree C. 1. Single salt parameters. 
    Journal of Chemical and Engineering Data, (2), 177–184. 
    
    May, P. M., Rowland, D., Hefter, G., & Königsberger, E. (2011). 
    A Generic and Updatable Pitzer Characterization of Aqueous Binary Electrolyte Solutions at 1 bar and 25 °C. 
    Journal of Chemical & Engineering Data, 56(12), 5066–5077. doi:10.1021/je2009329
    
    Beyer, R., & Steiger, M. (2010). Vapor Pressure Measurements of NaHCOO + H 2 O and KHCOO + H 2 O from 278 to 308 K 
    and Representation with an Ion Interaction (Pitzer) Model. 
    Journal of Chemical & Engineering Data, 55(2), 830–838. doi:10.1021/je900487a
    
    See also:
    --------
    _debye_parameter_activity
    _pitzer_B_MX
    _pitzer_B_gamma
    _pitzer_B_phi
    _pitzer_log_gamma

    
    '''
    # assign proper units to alpha1, alpha2, and b
    alpha1 = alpha1* unit('kg ** 0.5 / mol ** 0.5')
    alpha2 = alpha2* unit('kg ** 0.5 / mol ** 0.5')
    b = b * unit('kg ** 0.5 / mol ** 0.5')
    C_phi = C_phi * unit('kg ** 2 /mol ** 2')
    
    # assign units appropriate for the activity parameters
    BMX = _pitzer_B_MX(ionic_strength,alpha1,alpha2,beta0,beta1,beta2) * unit('kg/mol')
    Bphi = _pitzer_B_phi(ionic_strength,alpha1,alpha2,beta0,beta1,beta2) * unit('kg/mol')
    
    loggamma = _pitzer_log_gamma(ionic_strength,molality,BMX,Bphi,C_phi,z_cation,z_anion,nu_cation,nu_anion,temperature,b)
    
    return math.exp(loggamma) 
    
def get_apparent_volume_pitzer(ionic_strength,molality,alpha1,alpha2,beta0,beta1,beta2,C_phi,V_o,z_cation,z_anion,nu_cation,nu_anion,temperature='25 degC',b=1.2):
    '''Return the apparent molar volume of solute in the parent solution according to the Pitzer model.
    
    Parameters:
    ----------
    ionic_strength: Quantity
                    The ionic strength of the parent solution, mol/kg
    molality:       Quantity
                    The molal concentration of the parent salt, mol/kg
    alpha1, alpha2: number
                    Coefficients for the Pitzer model. This function assigns the coefficients
                    proper units of kg ** 0.5 / mol ** 0.5 after they are entered.
    beta0, beta1, beta2, C_phi: number
                    Pitzer coefficients for the apparent molar volume. 
                    These ion-interaction parameters are specific to each salt system.
    V_o:  number
                    The V^o Pitzer coefficient for the apparent molar volume.
    z_cation, z_anion: int
                    The formal charge on the cation and anion, respectively
    nu_cation, nu_anion: int
                    The stoichiometric coefficient of the cation and anion in the salt
    temperature:    str Quantity
                    String representing the temperature of the solution. Defaults to '25 degC' if not specified.
    b:              number, optional
                    Coefficient. Usually set equal to 1.2 and 
                    considered independent of temperature and pressure. If provided, this
                    coefficient is assigned proper units of kg ** 0.5 / mol ** 0.5  after
                    entry.
    
    Returns:
    -------
    float
        The apparent molar volume of the solute, cm **3 / mol
    
    Examples:
    --------  
    NOTE: the example below is for comparison with experimental and modeling data presented in
    the Krumgalz et al reference below. 
    
    0.25 mol/kg CuSO4. Expected result (from graph) = 0.5 cm ** 3 / mol 
    >>> get_apparent_volume_pitzer(1.0*unit('mol/kg'),0.25*unit('mol/kg'),1.4,12,0.001499,-0.008124,0.2203,-0.0002589,-6,2,-2,1,1,b=1.2)
    0.404...    
    
    1.0 mol/kg CuSO4. Expected result (from graph) = 4 cm ** 3 / mol 
    >>> get_apparent_volume_pitzer(4.0*unit('mol/kg'),1.0*unit('mol/kg'),1.4,12,0.001499,-0.008124,0.2203,-0.0002589,-6,2,-2,1,1,b=1.2)
    4.424...
    
    10.0 mol/kg ammonium nitrate. Expected result (from graph) = 50.3 cm ** 3 / mol 
    >>> get_apparent_volume_pitzer(10.0*unit('mol/kg'),10.0*unit('mol/kg'),2,0,0.000001742,0.0002926,0,0.000000424,46.9,1,-1,1,1,b=1.2)
    50.286...
    
    20.0 mol/kg ammonium nitrate. Expected result (from graph) = 51.2 cm ** 3 / mol 
    >>> get_apparent_volume_pitzer(20.0*unit('mol/kg'),20.0*unit('mol/kg'),2,0,0.000001742,0.0002926,0,0.000000424,46.9,1,-1,1,1,b=1.2)
    51.145...
    
    NOTE: the examples below are for comparison with experimental and modeling data presented in
    the Krumgalz et al reference below. 
    
    0.8 mol/kg NaF. Expected result = 0.03  
    >>> get_apparent_volume_pitzer(0.8*unit('mol/kg'),0.8*unit('mol/kg'),2,0,0.000024693,0.00003169,0,-0.000004068,-2.426,1,-1,1,1,b=1.2)
    0.22595 ...
    
    
    References:
    ----------  
    May, P. M., Rowland, D., Hefter, G., & Königsberger, E. (2011). 
    A Generic and Updatable Pitzer Characterization of Aqueous Binary Electrolyte Solutions at 1 bar and 25 °C. 
    Journal of Chemical & Engineering Data, 56(12), 5066–5077. doi:10.1021/je2009329

    Krumgalz, Boris S., Pogorelsky, Rita (1996).
    Volumetric Properties of Single Aqueous Electrolytes from Zero to Saturation Concentration at 298.15 K 
    Represented by Pitzer's Ion-Interaction Equations.
    Journal of Physical Chemical Reference Data, 25(2), 663-689.
    
    See also:
    --------
    _debye_parameter_volume
    _pitzer_B_MX

    
    '''
    # TODO - find a cleaner way to make sure coefficients are assigned the proper units
    # if they aren't, the calculation gives very wrong results
    alpha1 = alpha1* unit('kg ** 0.5 / mol ** 0.5')
    alpha2 = alpha2* unit('kg ** 0.5 / mol ** 0.5')
    b = b * unit('kg ** 0.5 / mol ** 0.5')
    C_phi = C_phi * unit('kg ** 2 /mol ** 2 / dabar')
    V_o = V_o * unit('cm ** 3 / mol')
    
    # assign units appropriate for the volume parameter
    BMX = _pitzer_B_MX(ionic_strength,alpha1,alpha2,beta0,beta1,beta2) * unit('kg /mol/dabar')
    
    second_term = (nu_cation + nu_anion) * abs(z_cation * z_anion) * (_debye_parameter_volume(temperature) / 2 / b) * math.log((1+b*ionic_strength ** 0.5))
    
    third_term = nu_cation * nu_anion * unit.R * unit(temperature) * \
    (2 * molality * BMX + molality ** 2 * C_phi * (nu_cation * nu_anion) ** 0.5)
    
    volume = V_o + second_term + third_term
    
    return volume.to('cm ** 3 / mol')
    
def _pitzer_f1(x):
    '''
    The function of ionic strength used to calculate \beta_MX in the Pitzer ion intercation model.
    
    .. math:: f(x) = 2 [ 1- (1+x) \exp(-x)] / x ^ 2
    
    References:
    ----------
    Scharge, T., Munoz, A.G., and Moog, H.C. (2012). Activity Coefficients of Fission Products in Highly
    Salinary Solutions of Na+, K+, Mg2+, Ca2+, Cl-, and SO42- : Cs+.
    /Journal of Chemical& Engineering Data (57), p. 1637-1647.
    
    Kim, H., & Jr, W. F. (1988). Evaluation of Pitzer ion interaction parameters of aqueous electrolytes at 25 degree C. 1. Single salt parameters. 
    Journal of Chemical and Engineering Data, (2), 177–184. 
    
    '''
    # return 0 if the input is 0
    if x == 0:
        return 0
    else:
        return 2 * ( 1 - (1+x) * math.exp(-x)) / x ** 2

def _pitzer_f2(x):
    '''
    The function of ionic strength used to calculate \beta_\gamma in the Pitzer ion intercation model.
    
    .. math:: f(x) = -{2 \\over x ^ 2} [ 1 - ({1+x+ x^2 \\over 2}) \\exp(-x)] 
    
    References:
    ----------
    Scharge, T., Munoz, A.G., and Moog, H.C. (2012). Activity Coefficients of Fission Products in Highly
    Salinary Solutions of Na+, K+, Mg2+, Ca2+, Cl-, and SO42- : Cs+.
    /Journal of Chemical& Engineering Data (57), p. 1637-1647.
    
    Kim, H., & Jr, W. F. (1988). Evaluation of Pitzer ion interaction parameters of aqueous electrolytes at 25 degree C. 1. Single salt parameters. 
    Journal of Chemical and Engineering Data, (2), 177–184. 
    
    '''
    # return 0 if the input is 0
    if x == 0:
        return 0
    else:
        return -2 * ( 1 - (1 + x + x ** 2 / 2) * math.exp(-x)) / x ** 2

def _pitzer_B_MX(ionic_strength,alpha1,alpha2,beta0,beta1,beta2):
    '''
    Return the B_MX coefficient for the Pitzer ion interaction model.
    
    .. math:: B_MX = \\beta_0 + \\beta_1 f1(\\alpha_1 I ^ {0.5}) + \\beta_2 f2(\\alpha_2 I ^ {0.5})
    
    Parameters:
    ----------
    ionic_strength: number
                    The ionic strength of the parent solution, mol/kg
    alpha1, alpha2: number
                    Coefficients for the Pitzer model, kg ** 0.5 / mol ** 0.5
    beta0, beta1, beta2: number
                    Coefficients for the Pitzer model. These ion-interaction parameters are
                    specific to each salt system.
                    
    Returns:
    -------
    float
            The B_MX parameter for the Pitzer ion interaction model.
    
    References:
    ----------
    Scharge, T., Munoz, A.G., and Moog, H.C. (2012). Activity Coefficients of Fission Products in Highly
    Salinary Solutions of Na+, K+, Mg2+, Ca2+, Cl-, and SO42- : Cs+.
    /Journal of Chemical& Engineering Data (57), p. 1637-1647.
    
    Kim, H., & Jr, W. F. (1988). Evaluation of Pitzer ion interaction parameters of aqueous electrolytes at 25 degree C. 1. Single salt parameters. 
    Journal of Chemical and Engineering Data, (2), 177–184.
    
    See Also:
    --------
    _pitzer_f1
    
    '''
    coeff = beta0 + beta1 * _pitzer_f1(alpha1 * ionic_strength ** 0.5) + beta2 * _pitzer_f1(alpha2 * ionic_strength ** 0.5)
    return coeff.magnitude

#def _pitzer_B_gamma(ionic_strength,alpha1,alpha2,beta1,beta2):
#    '''
#    Return the B^\gamma coefficient for the Pitzer ion interaction model.
#    
#    .. math:: B_\gamma = [ \beta_1 f2(\alpha_1 I ^ 0.5) + beta_2 f2(\alpha_2 I^0.5) ] / I
#    
#    Parameters:
#    ----------
#    ionic_strength: number
#                    The ionic strength of the parent solution, mol/kg
#    alpha1, alpha2: number
#                    Coefficients for the Pitzer model, kg ** 0.5 / mol ** 0.5.
#    beta1, beta2: number
#                    Coefficients for the Pitzer model. These ion-interaction parameters are
#                    specific to each salt system.
#                    
#    Returns:
#    -------
#    float
#            The B^gamma parameter for the Pitzer ion interaction model.
#    
#    References:
#    ----------
#    Scharge, T., Munoz, A.G., and Moog, H.C. (2012). Activity Coefficients of Fission Products in Highly
#    Salinary Solutions of Na+, K+, Mg2+, Ca2+, Cl-, and SO42- : Cs+.
#    /Journal of Chemical& Engineering Data (57), p. 1637-1647.
#    
#    Kim, H., & Jr, W. F. (1988). Evaluation of Pitzer ion interaction parameters of aqueous electrolytes at 25 degree C. 1. Single salt parameters. 
#    Journal of Chemical and Engineering Data, (2), 177–184.
#    
#    See Also:
#    --------
#    _pitzer_f2
#    
#    '''
#    coeff = (beta1 * _pitzer_f2(alpha1 * ionic_strength ** 0.5) + beta2 * _pitzer_f2(alpha2 * ionic_strength ** 0.5)) / ionic_strength
#    return coeff * unit('kg/mol')

def _pitzer_B_phi(ionic_strength,alpha1,alpha2,beta0,beta1,beta2):
    '''
    Return the B^\Phi coefficient for the Pitzer ion interaction model.
    
    .. math:: B^\Phi = \\beta_0 + \\beta1 \exp(-\\alpha_1 I ^{0.5}) + \\beta_2 \exp(-\\alpha_2 I ^ {0.5})
    
    or 
    
    .. math:: B^\\Phi = B^\\gamma - B_{MX}
    
    Parameters:
    ----------
    ionic_strength: number
                    The ionic strength of the parent solution, mol/kg
    alpha1, alpha2: number
                    Coefficients for the Pitzer model, kg ** 0.5 / mol ** 0.5
    beta0, beta1, beta2: number
                    Coefficients for the Pitzer model. These ion-interaction parameters are
                    specific to each salt system.
                    
    Returns:
    -------
    float
            The B^Phi parameter for the Pitzer ion interaction model.
    
    References:
    ----------
    Scharge, T., Munoz, A.G., and Moog, H.C. (2012). Activity Coefficients of Fission Products in Highly
    Salinary Solutions of Na+, K+, Mg2+, Ca2+, Cl-, and SO42- : Cs+.
    /Journal of Chemical& Engineering Data (57), p. 1637-1647.
    
    Kim, H., & Jr, W. F. (1988). Evaluation of Pitzer ion interaction parameters of aqueous electrolytes at 25 degree C. 1. Single salt parameters. 
    Journal of Chemical and Engineering Data, (2), 177–184.
    
    Beyer, R., & Steiger, M. (2010). Vapor Pressure Measurements of NaHCOO + H 2 O and KHCOO + H 2 O from 278 to 308 K 
    and Representation with an Ion Interaction (Pitzer) Model. 
    Journal of Chemical & Engineering Data, 55(2), 830–838. doi:10.1021/je900487a
    
    '''
    coeff = beta0 + beta1 * math.exp(-alpha1 * ionic_strength ** 0.5) + beta2 * math.exp(-alpha2* ionic_strength ** 0.5)
    return coeff

#def _pitzer_C_MX(C_phi,z_cation,z_anion):
#    '''
#    Return the C^\Phi coefficient for the Pitzer ion interaction model.
#    
#    .. math:: C_MX = C^\Phi / 2 \sqrt( \abs(z_+ z_-))
#    
#    Parameters:
#    ----------
#    C_phi: number
#                    The C_phi paramter for the Pitzer ion interaction model.
#    z_cation, z_anion: int
#                    The formal charge on the cation and anion, respectively
#                    
#    Returns:
#    -------
#    float
#            The C_MX parameter for the Pitzer ion interaction model.
#    
#    References:
#    ----------       
#    Kim, H., & Jr, W. F. (1988). 
#    Evaluation of Pitzer ion interaction parameters of aqueous electrolytes at 25 degree C. 1. Single salt parameters. 
#    Journal of Chemical and Engineering Data, (2), 177–184.
#    
#    May, P. M., Rowland, D., Hefter, G., & Königsberger, E. (2011). 
#    A Generic and Updatable Pitzer Characterization of Aqueous Binary Electrolyte Solutions at 1 bar and 25 °C. 
#    Journal of Chemical & Engineering Data, 56(12), 5066–5077. doi:10.1021/je2009329
#    '''
#    
#    coeff = C_phi / ( 2 * abs(z_cation * z_anion) ** 0.5 ) 
#    return coeff * unit('kg ** 2 /mol ** 2')

def _pitzer_log_gamma(ionic_strength,molality,B_MX,B_phi,C_phi,z_cation,z_anion,nu_cation,nu_anion,temperature='25 degC',b=1.2):
    '''
    Return the natural logarithm of the binary activity coefficient calculated by the Pitzer
    ion interaction model.
    
    .. math:: \\ln \\gamma_{MX} = -{|z_+ z_-| A^{Phi} ( I ^ {0.5} \\over (1 + b I ^ {0.5})} + {2 \\over b }\\ln (1 + b I ^ {0.5}) )+\
    + {m (2 \\nu_+ \\nu_-) \\over (\\nu_+ + \\nu_-)} (B_{MX} + B_{MX}^\\Phi) + {m^2(3 (\\nu_+ \\nu_-)^{1.5} \\over (\\nu_+ + \\nu_-))} C_{MX}^\\Phi    
    
    
    Parameters:
    ----------
    ionic_strength: number
                    The ionic strength of the parent solution, mol/kg
    molality:       number
                    The concentration of the salt, mol/kg
    B_MX,B_phi,C_phi: number
                    Calculated paramters for the Pitzer ion interaction model.
    z_cation, z_anion: int
                    The formal charge on the cation and anion, respectively
    nu_cation, nu_anion: int
                    The stoichiometric coefficient of the cation and anion in the salt
    temperature:    str Quantity
                    String representing the temperature of the solution. Defaults to '25 degC' if not specified.
    b:              number, optional
                    Coefficient. Usually set equal to 1.2 kg ** 0.5 / mol ** 0.5 and considered independent of temperature and pressure
                    
    Returns:
    -------
    float
            The natural logarithm of the binary activity coefficient calculated by the Pitzer ion interaction model.
    
    References:
    ----------       
    Kim, H., & Jr, W. F. (1988). 
    Evaluation of Pitzer ion interaction parameters of aqueous electrolytes at 25 degree C. 1. Single salt parameters. 
    Journal of Chemical and Engineering Data, (2), 177–184.
    
    May, P. M., Rowland, D., Hefter, G., & Königsberger, E. (2011). 
    A Generic and Updatable Pitzer Characterization of Aqueous Binary Electrolyte Solutions at 1 bar and 25 °C. 
    Journal of Chemical & Engineering Data, 56(12), 5066–5077. doi:10.1021/je2009329
    '''
    first_term = -1 * abs(z_cation * z_anion) * _debye_parameter_osmotic(temperature) * (ionic_strength ** 0.5 / (1+ b*ionic_strength ** 0.5) + 2/b * math.log(1+b*ionic_strength**0.5))
    second_term = 2 * molality * nu_cation * nu_anion / (nu_cation + nu_anion) * (B_MX + B_phi)
    third_term = 3 * molality ** 2 * (nu_cation * nu_anion) ** 1.5 / (nu_cation + nu_anion) * C_phi
    
    ln_gamma = first_term + second_term + third_term
    
    return ln_gamma

def get_osmotic_coefficient_TCPC(ionic_strength,S,b,n,formal_charge=1,counter_formal_charge=-1,stoich_coeff=1,counter_stoich_coeff=1,temperature=25):
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
    formal_charge : int, optional           
                        The charge on the solute, including sign. Defaults to +1 if not specified.
    counter_formal_charge : int, optional           
                        The charge on the solute's complementary ion, including sign. Defaults to -1 if not specified.
                        E.g. if the solute is Na+ and the salt is NaCl, counter_formal_charge = -1
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
    _debye_parameter_osmotic
    get_ionic_strength
    
    Notes:
    ------
    Valid for concentrated solutions up to saturation. Accuracy compares well with the Pitzer approach. See Reference [1] for a compilation of the appropriate parameters for a variety of commonly-encountered electrolytes.
    
    .. [1] Ge, Xinlei, Wang, Xidong, Zhang, Mei, and Seetharaman, Seshadri. "Correlation and Prediction of Activity and Osmotic Coefficients of Aqueous Electrolytes at 298.15 K by the Modified TCPC Model." J. Chemical Engineering Data 52, pp.538-547, 2007.
    '''
    # compute the 2nd term
    term2 = - math.fabs(formal_charge * counter_formal_charge) * _debye_parameter_osmotic(temperature) * ionic_strength ** 0.5 / (1 + b * ionic_strength ** 0.5)
    # compute the 3rd term
    term3 = S / (kelvin(temperature) * ( stoich_coeff + counter_stoich_coeff)) * 2 * n / (2 * n + 1) * ionic_strength  ** (2 * n)
    # add and return the osmotic coefficient
    return 1 - term2 + term3
    
def get_osmotic_coefficient_pitzer(ionic_strength,molality,alpha1,alpha2,beta0,beta1,beta2,C_phi,z_cation,z_anion,nu_cation,nu_anion,temperature='25 degC',b=1.2):
    '''Return the osmotic coefficient of water in an electrolyte solution according to the Pitzer model.
    
    Parameters:
    ----------
    ionic_strength: Quantity
                    The ionic strength of the parent solution, mol/kg
    molality:       Quantity
                    The molal concentration of the parent salt, mol/kg
    alpha1, alpha2: number
                    Coefficients for the Pitzer model. This function assigns the coefficients
                    proper units of kg ** 0.5 / mol ** 0.5 after they are entered.
    beta0, beta1, beta2, C_phi
                    Coefficients for the Pitzer model. These ion-interaction parameters are
                    specific to each salt system.
    z_cation, z_anion: int
                    The formal charge on the cation and anion, respectively
    nu_cation, nu_anion: int
                    The stoichiometric coefficient of the cation and anion in the salt
    temperature:    str Quantity
                    String representing the temperature of the solution. Defaults to '25 degC' if not specified.
    b:              number, optional
                    Coefficient. Usually set equal to 1.2 and 
                    considered independent of temperature and pressure. If provided, this
                    coefficient is assigned proper units of kg ** 0.5 / mol ** 0.5  after
                    entry.
    
    Returns:
    -------
    Quantity
        The osmotic coefficient of water, dimensionless
    
    Examples:
    --------
    Experimental value according to Beyer and Stieger reference is 1.3550
    
    >>> get_osmotic_coefficient_pitzer(10.175*unit('mol/kg'),10.175*unit('mol/kg'),1,0.5,-.0181191983,-.4625822071,.4682,.000246063,1,-1,1,1,b=1.2)
    1.3552 ...
    
    Experimental value according to Beyer and Stieger reference is 1.084
    
    >>> get_osmotic_coefficient_pitzer(5.6153*unit('mol/kg'),5.6153*unit('mol/kg'),3,0.5,0.0369993,0.354664,0.0997513,-0.00171868,1,-1,1,1,b=1.2)
    1.0850 ...  
    
    NOTE: the examples below are for comparison with experimental and modeling data presented in
    the May et al reference below. 
    
    10 mol/kg ammonium nitrate. Estimated result (from graph) = 0.62    
    >>> get_osmotic_coefficient_pitzer(10*unit('mol/kg'),10*unit('mol/kg'),2,0,-0.01709,0.09198,0,0.000419,1,-1,1,1,b=1.2)
    0.6143 ...
    
    5 mol/kg ammonium nitrate. Estimated result (from graph) = 0.7
    >>> get_osmotic_coefficient_pitzer(5*unit('mol/kg'),5*unit('mol/kg'),2,0,-0.01709,0.09198,0,0.000419,1,-1,1,1,b=1.2)
    0.6925 ...
    
    18 mol/kg ammonium nitrate. Estimated result (from graph) = 0.555
    >>> get_osmotic_coefficient_pitzer(18*unit('mol/kg'),18*unit('mol/kg'),2,0,-0.01709,0.09198,0,0.000419,1,-1,1,1,b=1.2)
    0.5556 ...
    
    References:
    ----------
    Scharge, T., Munoz, A.G., and Moog, H.C. (2012). Activity Coefficients of Fission Products in Highly
    Salinary Solutions of Na+, K+, Mg2+, Ca2+, Cl-, and SO42- : Cs+.
    /Journal of Chemical& Engineering Data (57), p. 1637-1647.
    
    Kim, H., & Jr, W. F. (1988). Evaluation of Pitzer ion interaction parameters of aqueous electrolytes at 25 degree C. 1. Single salt parameters. 
    Journal of Chemical and Engineering Data, (2), 177–184. 
    
    May, P. M., Rowland, D., Hefter, G., & Königsberger, E. (2011). 
    A Generic and Updatable Pitzer Characterization of Aqueous Binary Electrolyte Solutions at 1 bar and 25 °C. 
    Journal of Chemical & Engineering Data, 56(12), 5066–5077. doi:10.1021/je2009329
    
    Beyer, R., & Steiger, M. (2010). Vapor Pressure Measurements of NaHCOO + H 2 O and KHCOO + H 2 O from 278 to 308 K 
    and Representation with an Ion Interaction (Pitzer) Model. 
    Journal of Chemical & Engineering Data, 55(2), 830–838. doi:10.1021/je900487a
    
    See also:
    --------
    _debye_parameter_activity
    _pitzer_B_MX
    _pitzer_B_gamma
    _pitzer_B_phi
    _pitzer_log_gamma
       
    
    '''
    # assign proper units to alpha1, alpha2, and b
    alpha1 = alpha1* unit('kg ** 0.5 / mol ** 0.5')
    alpha2 = alpha2* unit('kg ** 0.5 / mol ** 0.5')
    b = b * unit('kg ** 0.5 / mol ** 0.5')
    C_phi = C_phi * unit('kg ** 2 /mol ** 2')
    B_phi = _pitzer_B_phi(ionic_strength,alpha1,alpha2,beta0,beta1,beta2) * unit('kg/mol')
    
    first_term = 1 - _debye_parameter_osmotic(temperature) * abs(z_cation * z_anion) * ionic_strength ** 0.5 / (1 + b * ionic_strength ** 0.5)
    second_term = molality * 2 * nu_cation * nu_anion / (nu_cation + nu_anion) * B_phi
    third_term = molality ** 2 * ( 2 * (nu_cation * nu_anion) ** 1.5 / (nu_cation + nu_anion)) * C_phi
 
    osmotic_coefficient = first_term + second_term + third_term
    
    return osmotic_coefficient