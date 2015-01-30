'''
pyEQL functions that take Solution objects as inputs or return Solution objects

'''

## Dependencies
# import libraries for scientific functions
import math

# internal pyEQL imports
import pyEQL

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
    Blend = pyEQL.Solution(temperature = blend_temperature,pressure= blend_pressure)
    
    # set or add the appropriate amount of all the components
    for item in mix_species.keys():
        if item in Blend.components:
            # if already present (e.g. H2O, H+), modify the amount
            Blend.set_amount(item,mix_species[item])
        else:
            # if not already present, add the component
            Blend.add_solute(item,mix_species[item])
            
    return Blend