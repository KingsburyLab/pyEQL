'''pyEQL - Python Library for Aquatic Chemistry Calculations
Developed by: Ryan Kingsbury
RyanSKingsbury@alumni.unc.edu

TODO - LICENSE INFO
'''

## Dependencies
# import libraries for scientific functions
import math

# internal pyEQL imports
from pyEQL import solution

# the pint unit registry
from pyEQL.parameter import unit
# TODO fix this to handle offsets the way pint wants us to since 0.7
unit.autoconvert_offset_to_baseunit = True

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

### Functions that operate on Solution Objects

def gibbs_mix(Solution1, Solution2):
    '''(Solution, Solution) -> float
    Return the Gibbs energy change associated with mixing two solutions

    Parameters
    ----------
    Solution1, Solution2 : Solution objects
        The two solutions to be mixed.
        
    Returns
    -------
    float
        The change in Gibbs eneryg associated with complete mixing of the
        Solutions, in Joules.
    
    Notes
    -----
    
    The Gibbs energy of mixing is calculated as follows: [1]_
        
    .. math::
        \Delta_{mix} G = \sum_i (n_c + n_d) R T \ln a_b - \sum_i n_c R T \ln a_c - \sum_i n_d R T \ln a_d
    
    Where n is the number of moles of substance, T is the temperature in kelvin,
    and  subscripts b, c, and refer to the concentrated, dilute, and blended
    Solutions, respectively. 
    
    Note that dissociated ions must be counted as separate components,
    so a simple salt dissolved in water is a three component solution (cation,
    anion, and water).
    
    References
    ----------
    
    .. [1] Koga, Yoshikata, 2007. //Solution Thermodynamics and its Application to Aqueous Solutions: 
    A differential approach.// Elsevier, 2007, pp. 23-37.
    
    Examples
    --------
 
    '''
    concentrate = Solution1
    dilute = Solution2
    blend = mix(Solution1,Solution2)
    term_list = {concentrate:0, dilute:0, blend:0}
    temperature = blend.get_temperature()

    # calculte the entropy change and number of moles solute for each solution
    for solution in term_list:
        for solute in solution.components:
            #print(solution.list_concentrations())
            if not solution.get_amount(solute,'fraction') == 0:
                term_list[solution] += solution.get_amount(solute,'mol') * math.log(solution.get_activity(solute))

    return (unit.R * temperature.to('K') * (term_list[blend] - term_list[concentrate] - term_list[dilute])).to('J')

def entropy_mix(Solution1, Solution2):
    '''(Solution, Solution) -> float
    Return the ideal mixing entropy associated with mixing two solutions

    Parameters
    ----------
    Solution1, Solution2 : Solution objects
        The two solutions to be mixed.
        
    Returns
    -------
    float
        The ideal mixing entropy associated with complete mixing of the
        Solutions, in Joules.
    
    Notes
    -----
    
    The ideal entropy of mixing is calculated as follows:[1]_
        
    .. math::
        \Delta_{mix} S = \sum_i (n_c + n_d) R T \ln x_b - \sum_i n_c R T \ln x_c - \sum_i n_d R T \ln x_d
    
    Where n is the number of moles of substance, T is the temperature in kelvin,
    and  subscripts b, c, and refer to the concentrated, dilute, and blended
    Solutions, respectively. 
    
    Note that dissociated ions must be counted as separate components,
    so a simple salt dissolved in water is a three component solution (cation,
    anion, and water).
    
    References
    ----------
    
    .. [1] Koga, Yoshikata, 2007. //Solution Thermodynamics and its Application to Aqueous Solutions: 
    A differential approach.// Elsevier, 2007, pp. 23-37.
    
    Examples
    --------
 
    '''
    concentrate = Solution1
    dilute = Solution2
    blend = mix(Solution1,Solution2)
    term_list = {concentrate:0, dilute:0, blend:0}
    temperature = blend.get_temperature()

    # calculate the entropy change and number of moles solute for each solution
    for solution in term_list:
        
        for solute in solution.components:
            #print(solution.list_concentrations())
            if not solution.get_amount(solute,'fraction') == 0:
                term_list[solution] += solution.get_amount(solute,'mol') * math.log(solution.get_amount(solute,'fraction'))

    return (unit.R * temperature.to('K') * (term_list[blend] - term_list[concentrate] - term_list[dilute])).to('J')

def donnan_eql(solution,fixed_charge):
    ''' Return a solution object in equilibrium with fixed_charge
    
    Parameters
    ----------
    Solution : Solution object
        The external solution to be brought into equilibrium with the fixed
        charges
    fixed_charge : str quantity
        String representing the concentration of fixed charges, including sign. 
        May be specified in mol/L or mol/kg units. e.g. '1 mol/kg'
        
    Returns
    -------
    Solution
        A solution that has established Donnan equilibrium with the external
        (input) Solution
    
    Notes
    -----
    
    The general equation representing the equilibrium between an external 
    electrolyte solution and an ion-exchange medium containing fixed charges
    is:[1]_
    
    .. math:: {a_- \\over \\bar a_-}^{1 \\over z_-} {\\bar a_+ \\over a_+}^{1 \\over z_+} \
    = exp({\\Delta \\pi \\bar V \\over {RT z_+ \\nu_+}})
    
    Where subscripts + and - indicate the cation and anion, respectively, 
    the overbar indicates the membrane phase,
    a represents activity, z represents charge, nu represents the stoichiometric
    coefficient, V represents the partial molar volume of the salt, and 
    delta pi is the difference in osmotic pressure between the membrane and the
    solution phase.
    
    In addition, electroneutrality must prevail within the membrane phase:
    
    .. math:: \\bar C_+ z_+ + \\bar X + \\bar C_- z_- = 0
    
    Where C represents concentration and X is the fixed charge concentration
    in the membrane or ion exchange phase.
    
    This function solves these two equations simultaneously to arrive at the 
    concentrations of the cation and anion in the membrane phase. It returns
    a solution equal to the input solution except that the concentrations of
    the predominant cation and anion have been adjusted according to this 
    equilibrium.
    
    NOTE that this treatment is only capable of equilibrating a single salt.
    This salt is identified by the get_salt() method.
    
    References
    ----------
    
    .. [1] Strathmann, Heiner, ed. //Membrane Science and Technology// vol. 9, 2004. \
    Chapter 2, p. 51. http://dx.doi.org/10.1016/S0927-5193(04)80033-0

    
    Examples
    --------
    TODO
    
    See Also
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
    
    # the stuff in the term below doesn't change on iteration, so calculate it up-front
    # assign it the correct units and extract the magnitude for a performance gain
    exp_term =  (molar_volume / (unit.R * solution.get_temperature() * z_cation * nu_cation)).to('1/Pa').magnitude
    
    def donnan_solve(x):
        '''Where x is the magnitude of co-ion concentration
        '''
        # solve for the counter-ion concentration by enforcing electroneutrality
        # using only floats / ints here instead of quantities helps performance
        if fixed_charge.magnitude >= 0:
            # counter-ion is the anion
            conc_cation_mem = x / abs(z_cation)
            conc_anion_mem = -(conc_cation_mem * z_cation + fixed_charge.magnitude) / z_anion
        elif fixed_charge.magnitude < 0:
            # counter-ion is the cation
            conc_anion_mem = x / abs(z_anion) 
            conc_cation_mem = -(conc_anion_mem * z_anion + fixed_charge.magnitude) / z_cation
        
        # match the units given for fixed_charge
        units = str(fixed_charge.units)
        
        # set the cation and anion concentrations in the membrane phase equal
        # to the current guess
        donnan_soln.set_amount(salt.cation,str(conc_cation_mem)+units)
        donnan_soln.set_amount(salt.anion,str(conc_anion_mem)+units)

        # get the new concentrations and activities
        act_cation_mem = donnan_soln.get_activity(salt.cation)
        act_anion_mem = donnan_soln.get_activity(salt.anion)
        
        # compute the difference in osmotic pressure
        # using the magnitudes here helps performance
        delta_pi = donnan_soln.get_osmotic_pressure().magnitude - solution.get_osmotic_pressure().magnitude
        
        return (act_cation_mem/act_cation_soln) ** (1/z_cation) * (act_anion_soln/act_anion_mem)**(1/z_anion) - math.exp(delta_pi * exp_term)

    # solve the function above using one of scipy's nonlinear solvers

    from scipy.optimize import brentq
    
    # determine which ion concentration represents the co-ion
    # call a nonlinear solver to adjust the concentrations per the donnan
    # equilibrium, unless the membrane is uncharged
    # the initial guess is to set the co-ion concentration in the membrane
    # equal to that in the solution
    if fixed_charge.magnitude >0:
        x = conc_cation_soln.magnitude
        brentq(donnan_solve,1e-10,x,xtol=0.001)
    elif fixed_charge.magnitude <0:
        x = conc_anion_soln.magnitude
        brentq(donnan_solve,1e-10,x,xtol=0.001)
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

    # set the pressure for the new solution
    p1 = Solution1.get_pressure()
    t1 = Solution1.get_temperature()
    v1 = Solution1.get_volume()
    p2 = Solution2.get_pressure()
    t2 = Solution2.get_temperature()
    v2 = Solution2.get_volume()

    # check to see if the solutions have the same temperature and pressure
    if not p1 == p2:
        logger.info('mix() function called between two solutions of different pressure. Pressures will be averaged (weighted by volume)')
          
    blend_pressure = str((p1 * v1 + p2 * v2) / (v1 + v2))
    
    if not t1 == t2:
        logger.info('mix() function called between two solutions of different temperature. Temperatures will be averaged (weighted by volume)')
    
    blend_temperature = str((t1 * v1 + t2 * v2) / (v1 + v2))
          
    # retrieve the amount of each component in the parent solution and 
    # store in a list.
    mix_species={}
    for item in Solution1.components:
        mix_species.update({item:str(Solution1.get_amount(item,'mol'))})
    for item in Solution2.components:
        if item in mix_species:
            new_amt = str(unit(mix_species[item]) + Solution2.get_amount(item,'mol'))
            mix_species.update({item:new_amt})
        else:
            mix_species.update({item:Solution2.get_amount(item,'mol')}) 
            
    
    # create an empty solution for the mixture
    Blend = solution.Solution(temperature = blend_temperature,pressure= blend_pressure)
    
    # set or add the appropriate amount of all the components
    for item in mix_species.keys():
        if item in Blend.components:
            # if already present (e.g. H2O, H+), modify the amount
            Blend.set_amount(item,mix_species[item])
        else:
            # if not already present, add the component
            Blend.add_solute(item,mix_species[item])
            
    return Blend