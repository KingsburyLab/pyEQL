'''pyEQL - Python Library for Aquatic Chemistry Calculations
Developed by: Ryan Kingsbury
RyanSKingsbury@alumni.unc.edu

WIP - LICENSE INFO
'''

import math

## Fundamental Constants
# Reference for all values: Wikipedia

# Avogadro's number, #/mole
CONST_Na = 6.02214129e23

# Universal gas constant, Joule/mole-Kelvin
CONST_R = 8.3144621 

# Fundamental charge, coulombs
CONST_e = 1.602176565e-19

# Permittivity of free space, Farad/meter
CONST_Eo = 8.854187817e-12

# Faraday constant, coulombs/mole (derived)
#9.64853399e4 
CONST_F = CONST_e * CONST_Na

# Boltzmann constant, joule/Kelvin (derived)
#1.3806488e-23
CONST_kb = CONST_R / CONST_Na

# Acceleration due to gravity, meters/second^2
CONST_g = 9.80665

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
    >>>kelvin(25)
    298.15
    '''
    return temp_celsius + 273.15
    
    
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
    >>>celsius(298)
    24.85
    '''
    return temp_kelvin - 273.15
    
    
def adjust_equilibrium_constant(equilibrium_constant,enthalpy,temperature,ref_temperature = 25):
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
    ref_temperature : float or int, optional
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
    >>>vanthoff(0.15,-197.6,42,25)
    0.0020356642823557407
    
    If the 'ref_temperature' parameter is omitted, a default of 25 C is used.
    
    >>>vanthoff(0.15,-197.6,42)
    0.0020356642823557407
    
    '''
    return equilibrium_constant * math.exp( enthalpy * 1000 / CONST_R * ( 1 / kelvin(ref_temperature) - 1 / kelvin(temperature)))
    
    
def adjust_diffusion_coefficient(diffusion_coefficient,temperature,ref_temperature=25):
    '''Adjust a diffusion coefficient for a different temperature
    The diffusion coefficients are corrected to temperature T (Kelvin) of the solution with:
        (Dw)T = (Dw)298 × (T / 298) × (η298 / ηT), where η is the viscosity of water. 
    WIP - FIND A LEGIT REFERENCE FOR THAT EQUATION
    
    .. [1] http://www.hydrochemistry.eu/exmpls/sc.html
    '''
    ## WIP - check this function (does it need dynamic or kinematic viscosity or does it matter?)
    return diffusion_coefficient * temperature / kelvin(ref_temperature) * water_viscosity_dynamic(ref_temperature)/water_viscosity_dynamic(temperature)

## Properties of Water

def water_density(temperature=25,pressure=101325):
    # WIP add pressure??
    # WIP more up to date equation??
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
    >>>water_density(25)
    997.04....
    
    '''
    return 999.65 + 0.20438 * temperature - 6.1744e-2 * temperature ** 1.5
    
    
def water_specific_weight(temperature):
    '''(number) -> float
    
    Return the specific weight of water in N/m3 at the specified temperature and pressure.
    
    Parameters:
    ----------
    temperature : float or int, optional
                  The temperature in Celsius. Defaults to 25 degrees if not specified.
    
    Returns:
    -------
    float
            The specific weight of water in N/m3.  
            
    See Also:
    --------
    water_density
    
    '''
    return water_density(temperature) * CONST_g


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
                  
    Notes:
    -----
    Implements the international equation for viscosity of water as specified by NIST[1]
    
    Valid for 273 < temperature < 1073 K and 0 < pressure < 100,000,000 Pa
    #WIP - add error checking for press and temp range
    
    .. [1] Sengers, J.V. "Representative Equations for the Viscosity of Water Substance." 
        J. Phys. Chem. Ref. Data 13(1), 1984.http://www.nist.gov/data/PDFfiles/jpcrd243.pdf
    
    Examples:
    --------
    >>>water_density_dynamic(25,100000)
    8.934...e-0.7
    >>>water_density_dynamic(100,25000000)
    2.979...e-0.7
    >>>water_density_dynamic(300,100000000)
    1.329...e-0.7
    #WIP - check these again after I implement density function
    
    '''
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
    
    return mu_o * mu_1


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
    return water_viscosity_dynamic(temperature,pressure) / water_density(temperature)
    

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
    >>>water_dielectric_constant(20)
    80.15...
    
    Display an error if 'temperature' is outside the valid range
    
    >>>water_dielectric_constant(-5)
    ERROR: Temperature specified exceeds range of data. Cannot extrapolate dielectric constant.
    
    '''
    # do not return anything if 'temperature' is outside the range for which
    # this fit applies
    if kelvin(temperature) < 273 or kelvin(temperature) > 372:
        print('ERROR: Temperature specified exceeds range of data. Cannot extrapolate dielectric constant.')
        return None
    
    # otherwise, calculate the dielectric constant using the quadratic fit    
    a = 0.24921e3
    b = -0.79069e0
    c = 0.72997e-3
    return a + b * kelvin(temperature) + c * kelvin(temperature) ** 2
    
    
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
    
    Notes:
    -----
    ### WIP - FIX THIS TO INCLUDE DENSITY in kg/m3
     The parameter A is equal to:[1]
         
     $$ A = {e^3 \over 8 \pi} sqrt{ 2 N_a \over (\epsilon_r \epsilon_o k_B T) ^3}

      This should not be confused with the Debye-Huckel constant for osmotic coefficients.
     
     WIP - FIND MORE CREDIBLE REFERENCE
     .. [1] http://en.wikipedia.org/wiki/Debye%E2%80%93H%C3%BCckel_equation
     
     
     
     The parameter A is equal to:[1]
     
     $$ A = 1.82e6 (\epsilon_r T) ^ -3/2 $$
    
    Note that when used in conjunction with the Debye-Huckel limiting law or related equations,
    this parameter is valid only when ionic strength is calculated from molar (mol/L) scale concentrations.
    
    .. [1] Stumm, Werner and Morgan, James J. Aquatic Chemistry, 3rd ed, 
        pp 103. Wiley Interscience, 1996.
        
    Examples:
    --------
    >>>water_debye_parameter_activity()
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
    return 1.8246e6 * (water_dielectric_constant(temperature) * kelvin(temperature)) ** -1.5
    
    # or this from MWH treatment book, page 302
    # WIP - document reference
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
    >>>water_debye_parameter_osmotic()
    0.392...
    
    See Also:
    --------
    water_debye_parameter
    
    '''
    # WIP - the factor 0.710 is a mystery number needed to make this equation return the correct value at 25C. I don't know why
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
    >>>p(1e-7)
    7.0
    >>>p(1.568e-9)
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
    >>>alpha(1,8,[4.7])
    1.00
    
    The sum of all alpha values should equal 1
    
    >>>alpha(0,8,[6.35,10.33])
    0.021..
    >>>alpha(1,8,[6.35,10.33])
    0.979..
    >>>alpha(2,8,[6.35,10.33])
    2.0435468121248873e-09
    
    If pH is equal to one of the pKa values the function should return 0.5.
    
    >>>alpha(1,6.35,[6.35,10.33])
    0.5
    
    The function will return an error if the number ofpKa's is less than n.
    
    >>> alpha(2,8,[])
    ERROR: insufficient number of pKa values given
    0.5   

    '''
    #generate an error if no pKa values are specified
    if len(pKa_list) == 0:
        print('ERROR: No pKa values given')
        return None
    
    #generate an error if n > number of pKa values
    if len(pKa_list) < n:
        print('ERROR: insufficient number of pKa values given')
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
    return numerator / denominator

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
        for solute in solution.ion_species:
            #print(solution.list_concentrations())
            my_solute = solute
            if not solution.get_mole_fraction(solute) == 0:
                term_list[solution] += solution.get_moles(solute) * math.log(solution.get_activity(solute))
        term_list[solution] += solution.get_moles_water() * math.log(solution.get_water_activity(my_solute))

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
        mole_fraction_water = solution.get_moles_water() / (solution.get_total_moles_solute() + solution.get_moles_water())
        for solute in solution.ion_species:
            #print(solution.list_concentrations())
            if not solution.get_mole_fraction(solute) == 0:
                term_list[solution] += solution.get_moles(solute) * math.log(solution.get_mole_fraction(solute))
        term_list[solution] += solution.get_moles_water() * math.log(mole_fraction_water)

    return CONST_R * kelvin(temperature) * (term_list[blend] - term_list[concentrate] - term_list[dilute])


def mix(Solution1, Solution2):
    '''(Solution, Solution) -> Solution
    Returns a new Solution object that results from the mixing of Solution1
    and Solution2
    
    '''
    # determine the volume for the new solution (assume additive)
    mix_vol = Solution1.get_volume() + Solution2.get_volume()
    
    # determine the density for the new solution. Assume volume is additive.
    # Convert from kg/L to kg/m3
    mix_dense = (Solution1.get_mass() + Solution2.get_mass()) / mix_vol *1000
    
    # determine the total mass of water in the blend
    mix_water_mass = Solution1.get_water_mass() + Solution2.get_water_mass()
    
    #conductivity will be initialized to zero in the new solution WIP
    
    #create a new Solution object for the blend
    Blend = Solution(mix_vol,mix_dense)
    #set the water mass, since we don't have any solutes yet
    Blend.water_mass = mix_water_mass
    
    # Loop through all the solutes in Solution1 and Solution2
    # determine the total moles, then divide by the new water mass to get concentration
    for item in Solution1.ion_species:
        Blend.add_solute(Solution1.ion_species[item])
        if item in Solution2.ion_species:
            Blend.ion_species[item].set_concentration((Solution1.get_moles(item) + Solution2.get_moles(item)) / mix_water_mass )
        else:
            Blend.ion_species[item].set_concentration(Solution1.get_moles(item) / mix_water_mass )
        if Solution1.ion_species[item].parameters_TCPC:
            Blend.ion_species[item].set_parameters_TCPC(Solution1.ion_species[item].get_parameters_TCPC('S'),Solution1.ion_species[item].get_parameters_TCPC('b'),Solution1.ion_species[item].get_parameters_TCPC('n'),Solution1.ion_species[item].get_parameters_TCPC('z_plus'),Solution1.ion_species[item].get_parameters_TCPC('z_minus'),Solution1.ion_species[item].get_parameters_TCPC('nu_plus'),Solution1.ion_species[item].get_parameters_TCPC('nu_minus'))
    
    for item in Solution2.ion_species:
        if not item in Solution1.ion_species:
            Blend.add_solute(Solution2.ion_species[item])
            Blend.ion_species[item].set_concentration(Solution2.get_moles(item) / mix_water_mass )
        if Solution2.ion_species[item].parameters_TCPC:
            Blend.ion_species[item].set_parameters_TCPC(Solution2.ion_species[item].get_parameters_TCPC('S'),Solution2.ion_species[item].get_parameters_TCPC('b'),Solution2.ion_species[item].get_parameters_TCPC('n'),Solution2.ion_species[item].get_parameters_TCPC('z_plus'),Solution2.ion_species[item].get_parameters_TCPC('z_minus'),Solution2.ion_species[item].get_parameters_TCPC('nu_plus'),Solution2.ion_species[item].get_parameters_TCPC('nu_minus'))
    
    return Blend
    
    
### Activity Functions
def get_activity_coefficient_debyehuckel(ionic_strength,valence=1,temperature=25):
    '''Return the activity coefficient of solute in the parent solution according to the Debye-Huckel limiting law.
    
    Parameters:
    ----------
    valence : int, optional      
                        The charge on the solute, including sign. Defaults to +1 if not specified.
    ionic_strength : number
                        The ionic strength of the parent solution, dimensionless
    temperature : float or int, optional
                        The solution temperature in degrees Celsius. 
                        Defaults to 25 degrees if omitted.
    Returns:
    -------
    float
         The mean molar (mol/L) scale ionic activity coefficient of solute

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
        print('WARNING: Ionic strength exceeds valid range of the Debye-Huckel limiting law')
    
    return - water_debye_parameter_activity(temperature) *valence ** 2 * math.sqrt(ionic_strength)

def get_activity_coefficient_guntelberg(ionic_strength,valence=1,temperature=25):
    '''Return the activity coefficient of solute in the parent solution according to the Guntelberg approximation.
    
    Parameters:
    ----------
    valence : int, optional          
                        The charge on the solute, including sign. Defaults to +1 if not specified.
    ionic_strength : number
                        The ionic strength of the parent solution, dimensionless
    temperature : float or int, optional
                        The solution temperature in degrees Celsius. 
                        Defaults to 25 degrees if omitted.
    Returns:
    -------
    float
         The mean molar (mol/L) scale ionic activity coefficient of solute
         
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
        print('WARNING: Ionic strength exceeds valid range of the Guntelberg approximation')
    
        return - water_debye_parameter_activity(temperature) * valence ** 2 * math.sqrt(ionic_strength) / (1+math.sqrt(ionic_strength))

def get_activity_coefficient_davies(ionic_strength,valence=1,temperature=25):
    '''Return the activity coefficient of solute in the parent solution according to the Davies equation.
    
    Parameters:
    ----------
    valence : int, optional           
                        The charge on the solute, including sign. Defaults to +1 if not specified.
    ionic_strength : number
                        The ionic strength of the parent solution, dimensionless
    temperature : float or int, optional
                        The solution temperature in degrees Celsius. 
                        Defaults to 25 degrees if omitted.
    Returns:
    -------
    float
         The mean molar (mol/L) scale ionic activity coefficient of solute

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
        print('WARNING: Ionic strength exceeds valid range of the Davies equation')
    
        return - water_debye_parameter_activity(temperature) * valence ** 2 * ( math.sqrt(ionic_strength) / (1+math.sqrt(ionic_strength)) - 0.2 * ionic_strength)
    
def get_activity_coefficient_TCPC(ionic_strength,S,b,n,valence=1,counter_valence=-1,stoich_coeff=1,counter_stoich_coeff=1,temperature=25):
    '''Return the activity coefficient of solute in the parent solution according to the modified TCPC model.
    
    Parameters:
    ----------
    ionic_strength : number
                        The ionic strength of the parent solution, dimensionless
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
        The mean molar (mol/L) scale ionic activity coefficient of solute

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
        The mean molar (mol/L) scale ionic activity coefficient of solute
    
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
                        The ionic strength of the parent solution, dimensionless
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

### Other Stuff - WIP
    
def debye_length(dielectric_constant,ionic_strength,temp):
    '''(number,number,number) -> float
    Return the Debye length of a solution in meters
    
    dielectric_constant is the dielectric constant of the solution
    ionic_strenght is the ionic strength in moles per cubic meter
    temp is the temperature in degrees celsius
    
    '''
    return math.sqrt( dielectric_constant * CONST_Eo * CONST_R * kelvin(temp) / (2 * CONST_F ** 2 * ionic_strength) )


def molar_conductivity(diffusion_coefficient,valence,temperature=25):
    '''(float,int,number) -> float
    Calculate the molar (equivalent) conductivity for a species in 
    Siemens-meters^2/mole
    
    Parameters:
    ----------
    diffusion_coefficient : float
                            The diffusion coefficient for the species 
                            in meters^2/second
    valence : int           
                            The charge on the species, including sign
                            
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
    
    Examples:
    --------
    For sodium ion:
        
    >>>molar_conductivity(1.334e-9,1)
    0.0050096...
    
    For sulfate ion at 30 C:
        
    >>>molar_conductivity(1.065e-9,2,30)
    0.0157340...
    '''
    return diffusion_coefficient * CONST_F ** 2 * valence ** 2 / (CONST_R * kelvin(temperature))



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

class Solute:
    '''represent each chemical species as an object containing its valence, transport numbers, concentration, activity, etc. '''
    def __init__(self,formula,molecular_weight=0,transport_num_cem = 0,transport_num_aem = 0,diffusion_coefficient=0,molality=0,activity=None,cost=0):
        self.formula = formula
        self.mw = molecular_weight
        self.t_cem = transport_num_cem
        self.t_aem = transport_num_aem
        self.molal = molality
        self.act = activity
        if activity == None:
            self.act = molality
        self.D = diffusion_coefficient
        self.unit_cost = cost
        self.parameters_TCPC={}
  
   #compute the ion's valence from the formula (only -3 to +3 supported right now)
        if self.formula.endswith('-'):
            self.charge = -1
        elif self.formula.endswith('+'):
            self.charge = 1
        elif self.formula.endswith('-2'):
            self.charge = -2
        elif self.formula.endswith('+2'):
            self.charge = 2
        elif self.formula.endswith('-3'):
            self.charge = -3
        elif self.formula.endswith('+3'):
            self.charge = 3
        elif self.formula.find('-')>=0 or self.formula.find('+')>=0:
            self.charge=0
            print('Warning: valence of ion %s out of bounds' % self.formula)
        else:
            self.charge=0
    
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
        
    def get_valence(self):
        return self.charge
        
    def get_transport_num_cem(self):
        return self.t_cem
        
    def get_transport_num_aem(self):
        return self.t_aem
        
    def get_diffusion_coefficient(self):
        return self.D
        
    def get_molecular_weight(self):
        return self.mw
    
    def get_concentration(self):
        return self.molal
    
    def set_concentration(self,conc):
        self.molal = conc
        
    def get_activity(self):
        return self.act
    
    def set_activity(self,act):
        self.act = act
        
    def get_unit_cost(self):
        return self.unit_cost
    
    def set_unit_cost(self,cost):
        self.unit_cost = cost
        
    def molar_conductivity(self):
        '''
        Return the molar conductivity (S/m /mol/m3) of the species
        Precondition: diffusion coefficient of species > 0
        '''
        #Calculate the molar conductivity at infinite dilution
        delta_o = self.charge ** 2 * self.D * const_F ** 2 / (const_R * temp_kelvin)
        return delta_o
        
    #set output of the print() statement
    def __str__(self):
        return 'Species ' + str(self.get_name()) + ' MW=' + str(self.get_molecular_weight()) +' Valence='+str(self.get_valence()) + ' t_cem='+str(self.get_transport_num_cem())+' t_aem='+str(self.get_transport_num_aem())+ ' Conc= ' + str(self.get_concentration()) + 'm  Activity= ' + str(self.get_activity())
        
class Solution:
    '''represent each solution with its own set of ions and bulk properties like volume(L), density(kg/m3) conductivity(S/m), etc.'''
    def __init__(self,volume=1,density=1000,conductivity=0):
        self.volume = volume
        self.density = density
        self.cond = conductivity
        self.ion_species={}
        self.ix_species=[]
        self.water_mass=0
        #parameter for the storage cost, $/L
        self.unit_storage_cost = 0
        self.calc_water_mass()
        
    def add_solute(self,solute):
        '''(solute object) -> None
        Creates a new copy of the solute object and associates it with Solution
        '''
        new_solute = Solute(solute.get_name(),solute.get_molecular_weight(),solute.get_transport_num_cem(),solute.get_transport_num_aem(),
        solute.get_diffusion_coefficient(),solute.get_concentration(),solute.get_activity(),solute.get_unit_cost())
        
        self.ion_species.update({new_solute.get_name():new_solute})
        
    def calc_solute_mass(self,ion):
        '''(str) -> float
        
        Return the total mass (kg) of an ion present in the solution
        '''
        solute_mass = self.ion_species[ion].get_concentration() * self.ion_species[ion].get_molecular_weight() * self.get_water_mass()/ 1000
        return solute_mass
        
    def calc_water_mass(self):
        '''()->float
        Used during __init__ to calculate the solvent (water) mass given the solution volume, density, and species concentrations
        '''

        #calculate total solution mass, kg
        solution_mass = self.density * self.volume / 1000
        
        #calculate total mass of solutes per kg water
        mass_per_kgw = 1
        for i in self.ion_species:
            mass_per_kgw += self.ion_species[i].get_concentration() * self.ion_species[i].get_molecular_weight() / 1000
        
        #calculate water mass by dividing
        self.water_mass = solution_mass / mass_per_kgw
        return self.water_mass
        
    def get_water_mass(self):
        return self.water_mass   
    
    ## WIP - deprecate this
    def set_storage_cost(self,cost):
        '''set the unit storage cost, $/L'''
        self.unit_storage_cost = cost
            
    def get_volume(self):
        return self.volume
    
    def get_mass(self):
        return self.volume * self.density /1000
        
    def get_density(self):
        return self.density
        
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
        #WIP
        
    
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
        
        
    def get_conductivity(self):
        return self.cond
        
    def get_unit_cost(self):
        return self.unit_storage_cost
        
    def get_solute(self,i):
        return self.ion_species[i]
        
    def list_solutes(self):
        return list(self.ion_species.keys())
    
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
        get_water_mass()
        
        Examples:
        --------
        TBD
        
        '''
        return self.ion_species[solute].get_concentration() * self.get_water_mass()
        
    def get_total_moles_solute(self):
        '''Return the total moles of all solute in the solution'''
        tot_mol = 0
        for item in self.ion_species:
            tot_mol += self.get_moles(item)
        return tot_mol
    
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
        get_water_mass()
        
        Notes:
        -----
        This function assumes water is the solvent with MW = 18
 
        Examples:
        --------
        TBD
        
        '''
        # calculate total moles of solvent and solutes
        tot_mol = self.get_moles_water() + self.get_total_moles_solute()
        # compute the fraction
        return self.get_moles(solute) / tot_mol
    
    def get_moles_water(self):
        return self.get_water_mass() * 1000 / 18
    
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
        get_water_mass()
        
        Notes:
        -----
        This function assumes water is the solvent with MW = 18
 
        Examples:
        --------
        TBD
        '''
        return self.get_moles(solute) / self.get_volume()
        
    def list_concentrations(self):
        '''() -> dict
        
        Return a dictionary containing a list of the species in solution paired with their molal concentration
        '''
        self.mol_list={}
        for i in self.ion_species.keys():
            self.mol_list.update({i:self.ion_species[i].get_concentration()})
        return self.mol_list
        
    def list_activity(self):
        '''() -> dict
        
        Return a dictionary containing a list of the species in solution paired with their molal activity
        '''
        self.act_list={}
        for i in self.ion_species.keys():
            self.act_list.update({i:self.ion_species[i].get_activity()})
        return self.act_list
    
    def list_mole_fractions(self):
        '''() -> dict
        
        Return a dictionary containing a list of the species in solution paired with their mole fraction
        '''
        self.fraction_list={}
        for i in self.ion_species.keys():
            self.fraction_list.update({i:self.get_mole_fraction(i)})
        # add mole fraction for water
        self.fraction_list.update({H2O:self.get_moles_water/self.get_
        return self.fraction_list
    
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
        ion = self.ion_species[solute]
        # for very low ionic strength, use the Debye-Huckel limiting law
        if self.get_ionic_strength() <= 0.005:
            return get_activity_coefficient_debyehuckel(self.get_ionic_strength(),temperature)
        # use the Guntelberg approximation for 0.005 < I < 0.1
        elif self.get_ionic_strength() <= 0.1:
            return get_activity_coefficient_guntelberg(self.get_ionic_strength(),ion.get_valence(),temperature)
        # use the Davies equation for 0.1 < I < 0.5
        elif self.get_ionic_strength() <= 0.5:
            return get_activity_coefficient_davies(self.get_ionic_strength(),ion.get_valence(),temperature)
        # use the TCPC model for higher ionic strengths, if the parameters have been set
        elif self.ion_species[solute].parameters_TCPC:
            return get_activity_coefficient_TCPC(self.get_ionic_strength(),ion.get_parameters_TCPC('S'),ion.get_parameters_TCPC('b'),ion.get_parameters_TCPC('n'),ion.get_parameters_TCPC('z_plus'),ion.get_parameters_TCPC('z_minus'),ion.get_parameters_TCPC('nu_plus'),ion.get_parameters_TCPC('nu_minus'),temperature)
            
        else:
            print('WARNING: Ionic strength too high to estimate activity. Specify parameters for Pitzer or TCPC methods. Returning unit activity coefficient')
            return 1.0
        # WIP - NEED TO TEST THIS FUNCTION
    
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
        The thermodynamic activity of the solute in question
        
        Notes:
        -----
        The thermodynamic activity is independent of the concentration scale used. However,
        the concentration and the activity coefficient must use corresponding scales.[1][2]
        In this module, ionic strenght, activity coefficients, and activities are all
        based on the molar (mol/L) concentration scale.
        
        References:
        ----------
        ..[1] http://adsorption.org/awm/utils/Activity.htm
        ..[2] http://en.wikipedia.org/wiki/Thermodynamic_activity#Activity_coefficient
        
        See Also:
        --------
        get_activity_coefficient
        get_ionic_strength
        
        '''
        return self.get_activity_coefficient(solute,temperature) * self.get_molar_concentration(solute)
        
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
        
        ion = self.ion_species[solute]
        
        if self.ion_species[solute].parameters_TCPC:
            return get_osmotic_coefficient_TCPC(self.get_ionic_strength(),ion.get_parameters_TCPC('S'),ion.get_parameters_TCPC('b'),ion.get_parameters_TCPC('n'),ion.get_parameters_TCPC('z_plus'),ion.get_parameters_TCPC('z_minus'),ion.get_parameters_TCPC('nu_plus'),ion.get_parameters_TCPC('nu_minus'),temperature)
        else:
            print('Cannot calculate water activity because TCPC parameters for solute are not specified. Returning unit osmotic coefficient')
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
        
        Return the ionic strength of the solution, calculated as 1/2 * sum ( molarity * valence ^2) over all the ions.
        Molar (mol/L) scale concentrations are used for compatibility with the activity correction formulas.
        
        
        Returns:
        -------
        float : 
            The molar scale ionic strength of the parent solution.
        
        >>> conc_soln.list_concentrations()
        {'Na+': 5.999375074924214, 'Cl-': 5.999904143046362, 'HCO3-': 0, 'NaCO3-': 0, 'NaHCO3': 0}
        >>> conc_soln.get_ionic_strength()
        5.999639608985288
        '''
        self.ionic_strength=0
        for solute in self.ion_species.keys():
            self.ionic_strength += 0.5 * self.get_molar_concentration(solute) * self.ion_species[solute].get_valence() ** 2
        return self.ionic_strength
            
    
   #set output of the print() statement for the solution     
    def __str__(self):
        return 'Ionic Species: '+str(self.list_solutes())+' Volume: '+str(self.get_volume())+'L  Density: '+str(self.get_density())+' kg/m3  Conductivity: '+str(self.get_conductivity())+'S/m'

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