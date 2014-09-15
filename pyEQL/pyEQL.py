'''pyEQL - Python Library for Aquatic Chemistry Calculations
Developed by: Ryan Kingsbury
RyanSKingsbury@alumni.unc.edu

TODO - LICENSE INFO
'''

## Dependencies
# import libraries for scientific functions
import math

# internal pyEQL imports
import pyEQL.activity_correction as ac
import pyEQL.water_properties as h2o

# the pint unit registry
from pyEQL.parameter import unit
# functions to manage importing paramters from database files and making them accessible to pyEQL
import pyEQL.database as database
from pyEQL.database import parameters_database as db


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

def adjust_temp_pitzer(c1,c2,c3,c4,c5,temp,temp_ref=unit('298.15 K')):
    '''
    Calculate a parameter for th e Pitzer model based on temperature-dependent
    coefficients c1,c2,c3,c4,and c5.
    
    Parameters:
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
    
    Parameters:
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
   
    Returns:
    -------
    float
        adjusted reaction equilibrium constant
    
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
    
    Parameters:
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
   
    Returns:
    -------
    Quantity
        The adjusted reaction equilibrium constant

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
    >>> adjust_temp_arrhenius(7,900*unit('kJ/mol'),37*unit('degC'),97*unit('degC')) #doctest: +ELLIPSIS
    1.8867225...e-24
    
    '''
    output = rate_constant * math.exp( activation_energy / unit.R * ( 1 / reference_temperature.to('K') - 1 / temperature.to('K')))
    
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
    return diffusion_coefficient * temperature / ref_temperature * water_viscosity_dynamic(ref_temperature)/water_viscosity_dynamic(temperature)

## Properties of Water


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

def gibbsmix(Solution1, Solution2,temperature=25*unit('degC')):
    '''(Solution, Solution) -> float
    Return the Gibbs energychange associated with mixing two solutions

    Parameters:
    ----------
    Solution1, Solution2 : Solution objects
        The two solutions to be mixed.
    temperature : Quantity, optional
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

    return unit.R * temperature.to('K') * (term_list[blend] - term_list[concentrate] - term_list[dilute])

def entropy_mix(Solution1, Solution2,temperature=25*unit('degC')):
    '''(Solution, Solution) -> float
    Return the ideal mixing entropy associated with mixing two solutions

    Parameters:
    ----------
    Solution1, Solution2 : Solution objects
        The two solutions to be mixed.
    temperature : Quantity, optional
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

    return unit.R * temperature.to('K') * (term_list[blend] - term_list[concentrate] - term_list[dilute])

def has_parameter(formula,name):
    '''
    Boolean test to determine whether a parameter exists in the database for a given species
    '''
    found = False
    try:
        for item in db[formula]:
            if item.get_name() == name:
                found = True
            else:
                continue
        
        return found
    except KeyError:
        logger.error('Species %s not found in database' % formula)
        return None    

def donnan_eql(solution,fixed_charge):
    ''' Return a solution object in equilibrium with fixed_charge
    
    Parameters:
    ----------
    Solution : Solution object
        The external solution to be brought into equilibrium with the fixed
        charges
    fixed_charge : str
        String representing the concentration of fixed charges, including sign. 
        May be specified in mol/L or mol/kg units. e.g. '1 mol/kg'
        
    Returns:
    -------
    Solution
        A solution that has established Donnan equilibrium with the external
        (input) Solution
    
    Notes:
    -----
    
    The general equation representing the equilibrium between an external 
    electrolyte solution and an ion-exchange medium containing fixed charges
    is:[1]
    
    $$ {a_-^s \over a_-^m}^{1 \over z_-} {a_+^s \over a_+^m}^{1 \over z_+} = exp({\DELTA \pi V \over {RT z_+ \nu_+}}) $$
    
    Where subscripts + and - indicate the cation and anion, respectively, 
    superscripts s and m indicate the membrane phase and the solution phase,
    a represents activity, z represents charge, nu represents the stoichiometric
    coefficient, V represents the partial molar volume of the salt, and 
    delta pi is the difference in osmotic pressure between the membrane and the
    solution phase.
    
    In addition, electroneutrality must prevail within the membrane phase:
    
    $$ C_+^m z_+ + X + C_-^m z_- = 0 $$
    
    Where C represents concentration and X is the fixed charge concentration
    in the membrane or ion exchange phase.
    
    This function solves these two equations simultaneously to arrive at the 
    concentrations of the cation and anion in the membrane phase. It returns
    a solution equal to the input solution except that the concentrations of
    the predominant cation and anion have been adjusted according to this 
    equilibrium.
    
    NOTE that this treatment is only capable of equilibrating a single salt.
    This salt is identified by the get_salt() method.
    
    .. [1] Strathmann, Heiner, ed. //Membrane Science and Technology// vol. 9, 2004. 
    Chapter 2, p. 51. 
    http://dx.doi.org/10.1016/S0927-5193(04)80033-0

    
    Examples:
    --------
    TODO
    
    See Also:
    --------
    get_salt()
    
    '''
    # identify the salt
    salt = solution.get_salt()
    
    # convert fixed_charge in to a quantity
    fixed_charge = unit(fixed_charge)
    
    # identify variables from the external solution
    conc_cation_soln = solution.get_amount(salt.cation,str(fixed_charge.units))
    conc_anion_soln = solution.get_amount(salt.anion,str(fixed_charge.units))
    act_cation_soln = solution.get_activity(salt.cation)
    act_anion_soln = solution.get_activity(salt.anion)
    z_cation= salt.z_cation
    z_anion = salt.z_anion
    nu_cation = salt.nu_cation
    
    # get the partial molar volume for the salt
    # TODO - consider how to incorporate pitzer parameters
    cation_vol = solution.get_solute(salt.cation).get_parameter('partial_molar_volume')
    anion_vol = solution.get_solute(salt.anion).get_parameter('partial_molar_volume')
    molar_volume = cation_vol + anion_vol
    
    # initialize the equilibrated solution - start with a direct copy of the 
    # input / external solution
    donnan_soln = solution.copy()
    
    # do nothing if either of the ion concentrations is zero
    if conc_cation_soln.magnitude == 0 or conc_anion_soln.magnitude == 0:
        return donnan_soln
    
    # define a function representing the donnan equilibrium as a function
    # of the two unknown actvities to feed to the nonlinear solver
    def donnan_solve(x):
        '''Where x is a numpy array of [co-ion concentration]
        '''
        # solve for the counter-ion concentration by enforcing electroneutrality
        # match the units given for fixed_charge
        if fixed_charge.magnitude >= 0:
            # counter-ion is the anion
            conc_cation_mem = unit(str(x[0]) + str(fixed_charge.units)) / abs(z_cation)
            conc_anion_mem = -(conc_cation_mem * z_cation + fixed_charge) / z_anion
        elif fixed_charge.magnitude < 0:
            # counter-ion is the cation
            conc_anion_mem = unit(str(x[0]) + str(fixed_charge.units)) / abs(z_anion)
            conc_cation_mem = -(conc_anion_mem * z_anion + fixed_charge) / z_cation
          
        # set the cation and anion concentrations in the membrane phase equal
        # to the current guess
        donnan_soln.set_amount(salt.cation,str(conc_cation_mem))
        donnan_soln.set_amount(salt.anion,str(conc_anion_mem))
        
        # get the new concentrations and activities
        act_cation_mem = donnan_soln.get_activity(salt.cation)
        act_anion_mem = donnan_soln.get_activity(salt.anion)
        
        # compute the difference in osmotic pressure
        delta_pi = donnan_soln.get_osmotic_pressure() - solution.get_osmotic_pressure()
        
        return (act_cation_mem/act_cation_soln) ** (1/z_cation) * (act_anion_soln/act_anion_mem)**(1/z_anion) - math.exp(delta_pi * molar_volume / (unit.R * solution.get_temperature() * z_cation * nu_cation))

    # solve the function above using one of scipy's nonlinear solvers
    # the initial guess is to set the cation concentration in the membrane
    # equal to that in the solution
    from scipy.optimize import broyden1
    
    # determine which ion concentration represents the co-ion
    # call a nonlinear solver to adjust the concentrations per the donnan
    # equilibrium, unless the membrane is uncharged
    if fixed_charge.magnitude >0:
        x = [conc_cation_soln.magnitude]
        broyden1(donnan_solve,x)
    elif fixed_charge.magnitude <0:
        x = [conc_anion_soln.magnitude]
        broyden1(donnan_solve,x)
    else:
        pass

    # return the equilibrated solution
    return donnan_soln

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
    volume : str, optional
                Volume of the solvent, including the unit. Defaults to '1 L' if omitted.
                Note that the total solution volume will be computed using partial molar
                volumes of the respective solutes as they are added to the solution.
    temperature : str, optional
                The solution temperature, including the unit. Defaults to '25 degC' if omitted.
    pressure : Quantity, optional
                The ambient pressure of the solution, including the unit. 
                Defaults to '1 atm' if omitted.
    
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
    
    def __init__(self,solutes=[],**kwargs):
        
        # initialize the voluem
        if 'volume' in kwargs:
            volume_set = True
            self.volume = unit(kwargs['volume'])
        else:
            volume_set = False
            self.volume = unit('1 L')
            
        # set the temperature
        if 'temperature' in kwargs:
            self.temperature = unit(kwargs['temperature'])
        else:
            self.temperature = unit('25 degC')
            
        # set the pressure
        if 'pressure' in kwargs:
            self.pressure = unit(kwargs['pressure'])
        else:
            self.pressure = unit('1 atm')
                        
        # create an empty dictionary of components
        self.components={}        
        
        # define the solvent
        if 'solvent' in kwargs:
            self.solvent_name=kwargs['solvent'][0]
            # warn if the solvent is anything besides water
            if not kwargs['solvent'][0] == 'H2O' or kwargs['solvent'][0] == 'water' :
                logger.error('Non-aqueous solvent detected. These are not yet supported!')
            
            # raise an error if the solvent volume has also been given
            if volume_set is True:
                logger.error('Solvent volume and mass cannot both be specified. Calculating volume based on solvent mass.')
            
            # add the solvent and the mass
            self.add_solvent(self.solvent_name,kwargs['solvent'][1])
        else:
            self.solvent_name = 'H2O'
            
            # calculate the solvent (water) mass based on the density and the solution volume
            self.add_solvent(self.solvent_name,str(self.volume * h2o.water_density(self.temperature)))

        
        # populate the solutes
        for item in solutes:
            self.add_solute(*item)        

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
        
        # add the volume occupied by the solute.
        self._update_volume()
        
    def add_solvent(self,formula,amount):
        '''Same as add_solute but omits the need to pass solvent mass to pint
        '''
        new_solvent = self.Solute(formula,amount,self.get_volume(),amount)
        self.components.update({new_solvent.get_name():new_solvent})
                        
    def get_solute(self,i):
        return self.components[i]
    
    def get_solvent(self):
        return self.components[self.solvent_name]
    
    def get_temperature(self):
        return self.temperature.to('degC')
    
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
                found = False
                
                for item in db[self.formula]:
                    if item.get_name() == parameter:
                        found = True
                        return item.get_value(temperature,pressure,ionic_strength)
                
                # Log an error if the parameter was not found
                if found == False:
                    logger.error('Required parameter %s not found for species %s.' % (parameter,self.formula))
            else:
                logger.error('No entry for species %s in parameters database' % self.formula)
                return None
                    
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
        
        def set_moles(self,amount,volume,solvent_mass):
            quantity = unit(amount)
            self.moles = quantity.to('moles','chem',mw=self.mw,volume=volume,solvent_mass=solvent_mass)  
  
        def get_activity(self):
            return self.activity
        
        def set_activity(self,activity):
            self.activity = activity

        def get_molar_conductivity(self,temperature=25*unit('degC')):
            # TODO - requires diffusion coefficient which may not be present
            '''(float,int,number) -> float
            Calculate the molar (equivalent) conductivity at infinte dilution for a species
            
            Parameters:
            ----------
            temperature : Quantity, optional
                        The temperature of the parent solution. Defaults to 25 degC if omitted.
                
            Returns:
            -------
            float
                    The molar or equivalent conductivity of the species at infinte dilution.
            
            Notes:
            -----
            Molar conductivity is calculated from the Nernst-Einstein relation:[1]
                
            .. math::
                \DELTA_i = {z_i^2 D_i F^2 \over RT}
            
            Note that the diffusion coefficient is strongly variable with temperature.
            
            .. [1] Smedley, Stuart. The Interpretation of Ionic Conductivity in Liquids, pp 1-9. Plenum Press, 1980.
            
            TODO Examples:
#             --------
#             
            '''
            diffusion_coefficient = self.get_parameter('diffusion_coefficient')
            
            molar_cond = diffusion_coefficient * (unit.e * unit.N_A) ** 2 * self.get_formal_charge() ** 2 / (unit.R * temperature)
            
            return molar_cond.to('mS / cm / (mol/L)')
                
        def mobility(self,temperature=25*unit('degC')):
            '''Return the ionic mobility of a species
            
            Parameters:
            ----------
            temperature : Quantity, optional
                        The temperature of the parent solution. Defaults to 25 degC if omitted.
            
            Returns:
            -------
            float : the ionic mobility
            
            
            Notes:
            -----
            This function uses the Einstein relation to convert a diffusion coefficient
            into an ionic mobility[1]
            
            .. math::
                \mu_i = {F |z_i| D_i \over RT}
            
            .. [1] Smedley, Stuart I. The Interpretation of Ionic Conductivity in Liquids. Plenum Press, 1980.
            
            '''
            mobility = unit.N_A * unit.e * abs(self.get_formal_charge()) * self.get_parameter('diffusion_coefficient') / (unit.R * temperature)
            
            logger.info('Computed ionic mobility as %s from D = %s at T=%S' % mobility,self.get_diffusion_coefficient(),temperature)
            
            return mobility
            
        #set output of the print() statement
        def __str__(self):
            return 'Species ' + str(self.get_name()) + ' MW=' + str(self.get_molecular_weight()) +' Valence='+str(self.get_formal_charge()) + ' Amount= ' + str(self.get_moles()) + 'moles  Activity= ' + str(self.get_activity())
    
    def get_solvent_mass(self):
        # return the total mass (kg) of the solvent
        solvent = self.get_solvent()
        mw = solvent.get_molecular_weight()
    
        return solvent.get_moles().to('kg','chem',mw=mw)
            
    def get_volume(self):
        return self.volume.to('L')
        
    def set_volume(self,volume):
        '''Change the total solution volume to volume, while preserving
        all component concentrations
        
        Parameters:
        ----------
        volume : str quantity
                Total volume of the solution, including the unit, e.g. '1 L' 
        
        Examples:
        ---------
        >>> mysol = Solution([['Na+','2 mol/L'],['Cl-','0.01 mol/L']],volume='500 mL')
        >>> print(mysol.get_volume())
        0.5000883925072983 l
        >>> mysol.list_concentrations()
        {'H2O': '55.508435061791985 mol/kg', 'Cl-': '0.00992937605907076 mol/kg', 'Na+': '2.0059345573880325 mol/kg'}
        >>> mysol.set_volume('200 mL')
        >>> print(mysol.get_volume())
        0.2 l
        >>> mysol.list_concentrations()
        {'H2O': '55.50843506179199 mol/kg', 'Cl-': '0.00992937605907076 mol/kg', 'Na+': '2.0059345573880325 mol/kg'}

        '''
        # figure out the factor to multiply the old concentrations by 
        scale_factor = unit(volume)/self.get_volume()
        
        # scale down the amount of all the solutes according to the factor        
        for item in self.components:
            self.set_amount(item,str(self.get_solute(item).get_moles()*scale_factor))
        
        # update the solution volume
        self.volume = unit(volume)
        
    def get_mass(self):
        '''returns the total solution mass in kg'''
        total_mass = 0
        for item in self.components:
            total_mass+= self.get_amount(item,'kg')
        return total_mass.to('kg')
        
    def get_density(self):
        return self.get_mass() / self.get_volume()
        
    def get_viscosity_dynamic(self):
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
        return h2o.water_viscosity_dynamic(self.get_temperature())
        #TODO
        
    
    def get_viscosity_kinematic(self):
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
        return self.get_viscosity_dynamic(self.get_temperature()) / self.get_density(self.get_temperature())
        
    
    def get_conductivity(self):
        '''
        Conductivity is calculated by summing the molar conductivities of the respective
        solutes, but they are activity-corrected and adjusted using an empricial exponent.
        This approach is used in PHREEQC and Aqion as described at these URLs:        
        <http://www.aqion.de/site/77>        
        <http://www.hydrochemistry.eu/exmpls/sc.html>        
        
        See Also:
        --------
        get_ionic_strength()
        get_molar_conductivity()
        get_activity_coefficient()
        
        '''
        EC = 0 * unit('S/m')
        temperature = self.get_temperature()
        
        for item in self.components:
            z = abs(self.get_solute(item).get_formal_charge())
            # ignore uncharged species            
            if not z == 0:
                # determine the value of the exponent alpha 
                if self.get_ionic_strength().magnitude < 0.36 * z:
                    alpha = 0.6 / z ** 0.5
                else:
                    alpha = self.get_ionic_strength().magnitude ** 0.5 / z
                
                EC += self.get_solute(item).get_molar_conductivity(temperature) * self.get_activity_coefficient(item) ** alpha * self.get_amount(item,'mol/L')
                             
        return EC.to('S/m')
    
    def get_osmotic_pressure(self):
        ''' 
        Return the osmotic pressure of the solution relative to pure water
        
        Parameters:
        ----------
        none
        
        Returns:
        -------
        Quantity
                The osmotic pressure of the solution relative to pure water in Pa
                
        References:
        ----------
        Sata, Toshikatsu. /Ion Exchange Membranes: Preparation, Characterization, and Modification./ Royal Society of Chemistry, 2004, p. 10.
        
        http://en.wikipedia.org/wiki/Osmotic_pressure#Derivation_of_osmotic_pressure
        
        Examples:
        --------
        If 'soln' is pure water, return 0
        >>> soln.get_osmotic_pressure()
        0.0
        
        If 'soln' is 0.5 mol/kg NaCl
        >>> soln.get_osmotic_pressure()
        2262808... pascal
        
        
        '''
        # TODO - tie this into parameter() and solvent() objects
        partial_molar_volume_water = 1.82e-5 *unit('m ** 3/mol')
        
        osmotic_pressure = - unit.R * self.get_temperature() / partial_molar_volume_water * math.log (self.get_water_activity())
        logger.info('Computed osmotic pressure of solution as %s Pa at T= %s degrees C' % (osmotic_pressure,self.get_temperature()))
        return osmotic_pressure.to('Pa')

## Concentration  Methods        
    
    def get_amount(self,solute,unit,temperature=25):
        '''returns the amount of 'solute' in the parent solution
       
        Parameters:
        ----------
        solute : str 
                    String representing the name of the solute of interest
        unit : str
                    Units desired for the output. Valid units are 'mol/L','mol/kg','mol', 'kg', and 'g/L'
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
        
        # with pint unit conversions enabled, we just pass the unit to pint
        return moles.to(unit,'chem',mw=mw,volume=self.get_volume(),solvent_mass=self.get_solvent_mass())
        
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
#            return moles / self.get_volume()
#        elif unit == 'mmol/L':
#            return moles * 1000 / self.get_volume()
#        elif unit == 'ng/L':
#            return moles * self.components[solute].get_molecular_weight() * 1e9 / self.get_volume()
#        elif unit == 'ug/L':
#            return moles * self.components[solute].get_molecular_weight() * 1e6 / self.get_volume()
#        elif unit == 'mg/L':
#            return moles * self.components[solute].get_molecular_weight() * 1e3 / self.get_volume()
#        elif unit == 'g/L':
#            return moles * self.components[solute].get_molecular_weight() / self.get_volume()
#        elif unit == 'kg/L':
#            return moles * self.components[solute].get_molecular_weight() / 1e3 / self.get_volume()
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
    
    def set_amount(self,solute,amount):
        '''Sets the amount of 'solute' in the parent solution.
       
        Parameters:
        ----------
        solute : str 
                    String representing the name of the solute of interest
        amount : str
                    String representing the concentration desired, e.g. '1 mol/kg'

        Returns:
        -------
        Nothing. The concentration of solute is modified.
        

        See Also:
        --------
        Solute.set_moles()
        '''
        # change the amount of the solute present
        self.get_solute(solute).set_moles(amount,self.get_volume(),self.get_solvent_mass())
        
        # update the solution volume
        self._update_volume()
        
    def get_total_moles_solute(self):
        '''Return the total moles of all solute in the solution'''
        tot_mol = 0
        for item in self.components:
            if item != self.solvent_name:
                tot_mol += self.components[item].get_moles()
        return tot_mol
    
    #TODO - figure out how best to integrate with pint / units
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
        return self.get_amount(solute,'moles') / (self.get_moles_water() + self.get_total_moles_solute())
    
    def get_moles_water(self):
        return self.get_amount(self.solvent_name,'mol',self.temperature)
    
    def get_salt(self):
        # identify the predominant salt in the solution
        import pyEQL.salt_ion_match as salt
        return salt.identify_salt(self)
      
## Activity-related methods
    def get_activity_coefficient(self,solute):
        '''Routine to determine the activity coefficient of a solute in solution. The correct function is chosen based on the ionic strength of the parent solution.
        
        Parameters:
        ----------
        solute : str 
                    String representing the name of the solute of interest
        
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
        temperature = self.get_temperature()
        # for very low ionic strength, use the Debye-Huckel limiting law
        
        if self.get_ionic_strength().magnitude <= 0.005:
            logger.info('Ionic strength = %s. Using Debye-Huckel to calculate activity coefficient.' % self.get_ionic_strength())
            return ac.get_activity_coefficient_debyehuckel(self.get_ionic_strength(),ion.get_formal_charge(),temperature)
            
        # use the Guntelberg approximation for 0.005 < I < 0.1
        elif self.get_ionic_strength().magnitude <= 0.1:
            logger.info('Ionic strength = %s. Using Guntelberg to calculate activity coefficient.' % self.get_ionic_strength())
            return ac.get_activity_coefficient_guntelberg(self.get_ionic_strength(),ion.get_formal_charge(),temperature)
            
        # use the Davies equation for 0.1 < I < 0.5
        elif self.get_ionic_strength().magnitude <= 0.5:
            logger.info('Ionic strength = %s. Using Davies equation to calculate activity coefficient.' % self.get_ionic_strength())
            return ac.get_activity_coefficient_davies(self.get_ionic_strength(),ion.get_formal_charge(),temperature)
              
        else:
            # identify the predominant salt in the solution
            Salt = self.get_salt()
            
            # search the database for pitzer parameters for 'salt'
            database.search_parameters(Salt.formula)
        
            found = False
            # use the Pitzer model for higher ionic strenght, if the parameters are available
            for item in db[Salt.formula]:
                if item.get_name() == 'pitzer_parameters_activity':
                    found == True
                    # TODO - fix inputs for alpha1 and alpha2
                    activity_coefficient=ac.get_activity_coefficient_pitzer(self.get_ionic_strength(), \
                    self.get_amount(solute,'mol/kg'),2,0,item.get_value()[0],item.get_value()[1],item.get_value()[2],item.get_value()[3], \
                    Salt.z_cation,Salt.z_anion,Salt.nu_cation,Salt.nu_anion,temperature)
                    
                    logger.info('Calculated activity coefficient of species %s as %s based on salt %s using Pitzer model' % (solute,activity_coefficient,Salt))
                    return activity_coefficient            
                    
            # TODO - fix TCPC implementation
            # use the TCPC model for higher ionic strengths, if the parameters have been set
            if self.components[solute].parameters_TCPC:
                logger.info('Ionic strength = %s. Using TCPC model to calculate activity coefficient.' % self.get_ionic_strength())
                return ac.get_activity_coefficient_TCPC(self.get_ionic_strength(),ion.get_parameters_TCPC('S'),ion.get_parameters_TCPC('b'),ion.get_parameters_TCPC('n'),ion.get_parameters_TCPC('z_plus'),ion.get_parameters_TCPC('z_minus'),ion.get_parameters_TCPC('nu_plus'),ion.get_parameters_TCPC('nu_minus'),temperature)
            
            if found == False:
                print('WARNING: Ionic strength too high to estimate activity. Specify parameters for Pitzer or TCPC models. Returning unit activity coefficient')
                return unit('1 dimensionless')
        # TODO - NEED TO TEST THIS FUNCTION
    
    def get_activity(self,solute):
        '''returns the thermodynamic activity of the solute in solution
       
        Parameters:
        ----------
        solute : str 
                    String representing the name of the solute of interest
        temperature :    Quantity, optional
                     The temperature of the solution. Defaults to 25 degrees C if omitted
        
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
            activity = self.get_water_activity()
        # TODO fix this for multivalent salts e.g. MgCl2
        else:
            activity = self.get_activity_coefficient(solute) * self.get_amount(solute,'mol/kg').magnitude
            logger.info('Calculated activity of solute %s as %s' % (solute,activity))
        
        return activity

    def get_osmotic_coefficient(self):
        '''calculate the osmotic coefficient

        Returns:
        -------
        Quantity : 
            The osmotic coefficient
            
        Notes:
        -----
        For ionic strengths below 0.5 mol/kg, the osmotic coefficient is assumed to equal 1.0.
        1.0 will also be returned at higher ionic strengths if appropriate Pitzer
        parameters are not supplied.
        '''
        temperature = self.get_temperature()
        ionic_strength = self.get_ionic_strength()
        
        import pyEQL.salt_ion_match as salt
        
        # identify the predominant salt in the solution
        Salt = self.get_salt()
        
        # set the concentration as the average concentration of the cation and
        # anion in the salt, accounting for stoichiometry
        concentration = (self.get_amount(Salt.cation,'mol/kg')/Salt.nu_cation + \
        self.get_amount(Salt.anion,'mol/kg')/Salt.nu_anion)/2
        
        # search the database for pitzer parameters for 'salt'
        database.search_parameters(Salt.formula)
        
        found = False        
        
        for item in db[Salt.formula]:
            if item.get_name() == 'pitzer_parameters_activity':
                found == True
                # TODO - fix inputs for alpha1 and alpha2
                osmotic_coefficient=ac.get_osmotic_coefficient_pitzer(ionic_strength, \
                concentration,2,0,item.get_value()[0],item.get_value()[1],item.get_value()[2],item.get_value()[3], \
                Salt.z_cation,Salt.z_anion,Salt.nu_cation,Salt.nu_anion,temperature)
                
                logger.info('Calculated osmotic coefficient of water as %s based on salt %s using Pitzer model' % (osmotic_coefficient,salt))
                return osmotic_coefficient
            
            # TODO - either deprecate or update to parameter framework
#            elif self.components[solute].parameters_TCPC:
#                ion = self.components[solute]
#                osmotic_coefficient= ac.get_osmotic_coefficient_TCPC(self.get_ionic_strength(),ion.get_parameters_TCPC('S'),ion.get_parameters_TCPC('b'),ion.get_parameters_TCPC('n'),ion.get_parameters_TCPC('z_plus'),ion.get_parameters_TCPC('z_minus'),ion.get_parameters_TCPC('nu_plus'),ion.get_parameters_TCPC('nu_minus'),temperature)
#                logger.info('Calculated osmotic coefficient of water as %s based on solute %s using TCPC model' % (osmotic_coefficient,solute))
#                return osmotic_coefficient
                
            else:
                continue
            
        if found == False:
            logger.warning('Cannot calculate osmotic coefficient because Pitzer parameters for solute are not specified. Returning unit osmotic coefficient')
            return unit('1 dimensionless')
    
    def get_water_activity(self):
        '''return the water activity
        
        Returns:
        -------
        float : 
            The thermodynamic activity of water in the solution.
            
        Examples:
        --------
        If 'soln' is a 0.5 mol/kg NaCl solution at 25 degC:
        >>> soln.get_water_activity()
        0.9835...
        
        If 'soln' is a 5.11 mol/kg NaHCO2 (sodium formate) solution at 25 degC:
        (literature value from Cabot specialty fluids is 0.82)
        >>> soln.get_water_activity()
        0.8631...
        
        Notes:
        -----
        Water activity is related to the osmotic coefficient in a solution containing i solutes by:[1]
        
        $$ ln a_w = - \Phi M_w \sum_i m_i    $$
        
        Where M_w is the molar mass of water (0.018015 kg/mol) and m_i is the molal concentration
        of each species.
        
        If appropriate Pitzer or TCPC model parameters are not available, the
        water activity is assumed equal to the mole fraction of water.
        
        References:
        ----------
        Blandamer, Mike J., Engberts, Jan B. F. N., Gleeson, Peter T., Reis, Joao Carlos R., 2005. "Activity of water in aqueous systems: A frequently neglected property."
        //Chemical Society Review// 34, 440-458.
        
        '''
        '''
        pseudo code
        
        identify predominant salt for coefficients
        check if coefficients exist for that salt
        if so => calc osmotic coefficient and log an info message
        
        if not = > return mole fraction and log a warning message
        
        '''
        osmotic_coefficient = self.get_osmotic_coefficient()
        
        if osmotic_coefficient == 1:
            logger.warning('Pitzer parameters not found. Water activity set equal to mole fraction')
            return self.get_mole_fraction('H2O')
        else:
            concentration_sum = unit('0 mol/kg')
            for item in self.components:                
                if item == 'H2O':
                    pass
                else:
                    concentration_sum += self.get_amount(item,'mol/kg')
                    
            logger.info('Calculated water activity using osmotic coefficient')  
            
            return math.exp(- self.get_osmotic_coefficient() * 0.018015*unit('kg/mol') * concentration_sum)
    
    def get_ionic_strength(self):
        '''() -> float
        
        Return the ionic strength of the solution, calculated as 1/2 * sum ( molality * charge ^2) over all the ions.
        Molal (mol/kg) scale concentrations are used for compatibility with the activity correction formulas, but
        the returned value does not carry units.
        
        
        Returns:
        -------
        float : 
            The ionic strength of the parent solution, mol/kg.
        
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
       
        return self.ionic_strength
    

    def get_debye_length(self):
        '''(number,number,number) -> float
        Return the Debye length of a solution in meters
        
        dielectric_constant is the dielectric constant of the solution
        ionic_strength is the ionic strength in moles per cubic meter
        temp is the temperature in degrees celsius
        
        '''
        # TODO - make dielectric constant dependent on ionic strength
        temperature = self.get_temperature()
        ionic_strength= self.get_ionic_strength()
        dielectric_constant = h2o.water_dielectric_constant()
        
        logger.warning('Debye length is being calculated using the dielectric constant for pure water. The influence \
        of ionic strength is not yet accounted for')
        
        return math.sqrt(dielectric_constant * unit.epsilon_0 * unit.R * temperature / (2 * unit.N_A * unit.e ** 2 * ionic_strength) )

    def get_transport_number(self,solute):
        '''Calculate the transport number of the solute in the solution
        
        $$ t_i = {D_i z_i^2 C_i \over \sum D_i z_i^2 C_i} $$
        
        Where C is the concentration in mol/L
        TODO - add more references / explanation here        

        '''
        denominator= 0
        numerator = 0
        
        for item in self.components:
            
            if item == 'H2O':
                continue
            
            term = self.get_solute(item).get_parameter('diffusion_coefficient') * \
            self.get_solute(item).get_formal_charge() ** 2 * \
            self.get_amount(item,'mol/L')
            
            if item == solute:
                numerator = term
            
            denominator += term
        
        return numerator / denominator
    
    def get_lattice_distance(self,solute):
        '''Calculate the average distance between molecules of the given solute,
        assuming that the molecules are uniformly distributed throughout the
        solution.
        
        Parameters:
        ----------
        solute : str 
                    String representing the name of the solute of interest
        
        Returns:
        -------
        Quantity : The average distance between solute molecules
        
        Examples:
        --------
        >>> soln = Solution([['Na+','0.5 mol/kg'],['Cl-','0.5 mol/kg']])
        >>> soln.get_lattice_distance('Na+')
        1.492964.... nanometer
        
        '''
        # calculate the volume per particle as the reciprocal of the molar concentration
        # (times avogadro's number). Take the cube root of the volume to get 
        # the average distance between molecules
        distance = (self.get_amount(solute,'mol/L') * unit.N_A) ** (-1/3)     
        
        return distance.to('nm')
    
    def _update_volume(self):
        '''
        Recalculate the solution volume based on composition
        
        '''    
        
        temperature = self.get_temperature()
        
        # identify the predominant salt in the solution
        Salt = self.get_salt()
        
        # search the database for pitzer parameters for 'salt'
        database.search_parameters(Salt.formula)
        
        # calculate the volume of the pure solvent
        self.volume = self.get_solvent_mass() / h2o.water_density(temperature)
        
        for item in self.components:
            
            solute = self.get_solute(item)            
            
            if item in ['H2O','H+','OH-','HOH']:
                continue
            
            # use the pitzer approach if parameters are available
            elif has_parameter(Salt.formula,'pitzer_parameters_volume'):
                
                for params in db[Salt.formula]:
                    if params.get_name() == 'pitzer_parameters_volume':
                
                        apparent_vol = ac.get_apparent_volume_pitzer(self.get_ionic_strength(), \
                        self.get_amount(item,'mol/kg'),2,0,params.get_value()[0],params.get_value()[1],params.get_value()[2],params.get_value()[3], \
                        params.get_value()[4],Salt.z_cation,Salt.z_anion,Salt.nu_cation,Salt.nu_anion,temperature)
                
                # TODO - fix so that this works for multivalent salts                
                self.volume += apparent_vol * solute.get_moles()
                logger.info('Updated solution volume using Pitzer model for solute %s' % item)
                
                # since the salt accounts for all solutes, stop                
                break
            
            elif has_parameter(item,'partial_molar_volume'):
                self.volume += solute.get_parameter('partial_molar_volume') * solute.get_moles()
                logger.info('Updated solution volume using direct partial molar volume for solute %s' % item)
                
            else:
                logger.warning('Partial molar volume data not available for solute %s. Solution volume will not be corrected.' % item)
        
    def copy(self):
        '''Return a copy of the solution
        
        TODO - clarify whether this is a deep or shallow copy
        '''
        # prepare to copy the bulk properties        
        new_temperature = str(self.get_temperature())
        new_pressure = str(self.pressure)
        new_solvent = self.solvent_name
        new_solvent_mass = str(self.get_solvent_mass())
        
        # create a list of solutes
        new_solutes = []        
        for item in self.components:
            # ignore the solvent
            if item == self.solvent_name:
                pass
            else:
                new_solutes.append([item,str(self.get_amount(item,'mol'))])
        
        # create the new solution
        return Solution(new_solutes,solvent=[new_solvent,new_solvent_mass],temperature=new_temperature,pressure=new_pressure)

    ## informational methods
        
    def list_solutes(self):
        return list(self.components.keys())
    
    def list_concentrations(self,unit='mol/kg'):
        '''() -> dict
        
        Return a dictionary containing a list of the species in solution paired with their amount in the specified units
        '''
        self.mol_list={}
        for i in self.components.keys():
            self.mol_list.update({i:str(self.get_amount(i,unit))})
        print('Component amounts (%s):\n' % unit,self.mol_list )
        
    def list_activities(self):
        '''() -> dict
        
        Return a dictionary containing a list of the species in solution paired with their molal activity
        '''
        self.act_list={}
        for i in self.components.keys():
            self.act_list.update({i:self.components[i].get_activity()})
        print('Component activities:\n',self.act_list )   

     
    def __str__(self):
        #set output of the print() statement for the solution     
        return 'Components: \n'+str(self.list_solutes()) + '\n' + 'Volume: '+str(self.get_volume()) + '\n' + 'Density: '+str(self.get_density())
      
        
# TODO - turn doctest back on when the nosigint error is gone        
#if __name__ == "__main__":
 #   import doctest
  #  doctest.testmod()