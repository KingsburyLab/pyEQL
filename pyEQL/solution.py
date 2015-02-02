'''
pyEQL Solution Class


'''

## Dependencies
# import libraries for scientific functions
import math

# internal pyEQL imports
import pyEQL.activity_correction as ac
import pyEQL.water_properties as h2o
import pyEQL.solute as sol

# the pint unit registry
from pyEQL import unit

# functions to manage importing paramters from database files and making them accessible to pyEQL
import pyEQL.database as database
from pyEQL.database import parameters_database as db

# logging system
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# add a filter to emit only unique log messages to the handler
import pyEQL.logging_system
unique = pyEQL.logging_system.Unique()
logger.addFilter(unique)


class Solution:
    '''Class representing the properties of a solution. Instances of this class 
    contain information about the solutes, solvent, and bulk properties.
    
    Parameters
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
    pH : number, optional
                Negative log of H+ activity. If omitted, the solution will be 
                initialized to pH 7 (neutral) with appropriate quantities of 
                H+ and OH- ions
    
    Returns
    -------
    A Solution object.
    
    Examples
    --------
    # Defining a 0.5M NaCl solution
    >>> solutes = [['Na+',23,0.0115],['Cl-',35,0.0175]]
    >>> solvent = ['H2O',18,1]
    >>> my_solution = Solution(solutes,solvent)
    >>> print(my_solution)
    Components: ['Na+', 'H2O', 'Cl-']
    Volume: 1.0290000000000001L
    Density: 1000 kg/m3
    
    
    See Also
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
        
        # initialize the volume
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

        # initialize the volume recalculation flag
        self.volume_update_required = False        
        
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

        # set the pH with H+ and OH-
        if 'pH' in kwargs:
            pH = kwargs['pH']
        else:
            pH = 7

        self.add_solute('H+',str(10 ** (-1 * pH))+'mol/L')
        self.add_solute('OH-',str(10 ** (-1 * (14-pH)))+'mol/L')        
        
        # populate the other solutes
        for item in solutes:
            self.add_solute(*item)        

    def add_solute(self,formula,amount,parameters={}):
        '''Primary method for adding substances to a pyEQL solution
        
        Parameters
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
        
        # if units are given on a per-volume basis, 
        # iteratively solve for the amount of solute that will preserve the
        # original volume and result in the desired concentration  
        if unit(amount).dimensionality == ('[substance]/[length]**3' or '[mass]/[length]**3'):
            
            # store the original volume for later
            orig_volume = self.get_volume()
            
            # add the new solute
            new_solute = sol.Solute(formula,amount,self.get_volume(),self.get_solvent_mass(),parameters)
            self.components.update({new_solute.get_name():new_solute})            
            
            # calculate the volume occupied by all the solutes
            solute_vol = self._get_solute_volume()
            
            # determine the volume of solvent that will preserve the original volume
            target_vol = orig_volume - solute_vol
            
            # adjust the amount of solvent
            target_mass = target_vol * h2o.water_density(self.get_temperature())
            mw = self.get_solvent().get_molecular_weight()
            target_mol = target_mass / mw
            self.get_solvent().moles = target_mol
            
        else:
            
            # add the new solute
            new_solute = sol.Solute(formula,amount,self.get_volume(),self.get_solvent_mass(),parameters)
            self.components.update({new_solute.get_name():new_solute})            
            
            # update the volume to account for the space occupied by all the solutes
            # make sure that there is still solvent present in the first place
            if self.get_solvent_mass() <= unit('0 kg'):
                logger.error('All solvent has been depleted from the solution')
                return None
            else:
                # set the volume recalculation flag
                self.volume_update_required = True
        
    def add_solvent(self,formula,amount):
        '''Same as add_solute but omits the need to pass solvent mass to pint
        '''
        new_solvent = sol.Solute(formula,amount,self.get_volume(),amount)
        self.components.update({new_solvent.get_name():new_solvent})
                        
    def get_solute(self,i):
        return self.components[i]
    
    def get_solvent(self):
        return self.components[self.solvent_name]
    
    def get_temperature(self):
        return self.temperature.to('K')
    
    def set_temperature(self,temperature):
        '''
        Set the solution temperature.
        
        Parameters
        ----------
        temperature : str
            String representing the temperature, e.g. '25 degC'
        '''
        self.temperature = unit(temperature)
    
    def get_pressure(self):
        return self.pressure.to('atm')
        
    def set_pressure(self,pressure):
        '''
        Set the hydrostatic pressure of the solution.
        
        Parameters
        ----------
        pressure : str
            String representing the temperature, e.g. '25 degC'
        '''
        self.pressure = unit(pressure)
    
    def get_solvent_mass(self):
        # return the total mass (kg) of the solvent
        solvent = self.get_solvent()
        mw = solvent.get_molecular_weight()
    
        return solvent.get_moles().to('kg','chem',mw=mw)
            
    def get_volume(self):
        # if the composition has changed, recalculate the volume first
        if self.volume_update_required is True:
            self._update_volume()
            self.volume_update_required = False
            
        return self.volume.to('L')
        
    def set_volume(self,volume):
        '''Change the total solution volume to volume, while preserving
        all component concentrations
        
        Parameters
        ----------
        volume : str quantity
                Total volume of the solution, including the unit, e.g. '1 L' 
        
        Examples
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
            self.get_solute(item).moles = self.get_solute(item).moles * scale_factor
        
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
    
    def get_viscosity_relative(self):
        '''
        Return the viscosity of the solution relative to that of water
        
        This is calculated using a simplified form of the Jones-Dole equation:
        
        .. math:: \\eta_{rel} = 1 + \\sum_i B_i m_i
        
        Where m is the molal concentration and B is an empirical parameter.
        
        See 
        <http://downloads.olisystems.com/ResourceCD/TransportProperties/Viscosity-Aqueous.pdf>
        <http://www.nrcresearchpress.com/doi/pdf/10.1139/v77-148>
        <http://apple.csgi.unifi.it/~fratini/chen/pdf/14.pdf>
        '''
        #if self.get_ionic_strength().magnitude > 0.2:
         #   logger.warning('Viscosity calculation has limited accuracy above 0.2m')
        
#        viscosity_rel = 1
#        for item in self.components:
#            # ignore water
#            if item != 'H2O':
#                # skip over solutes that don't have parameters
#                try:
#                    conc = self.get_amount(item,'mol/kg').magnitude
#                    coefficients= self.get_solute(item).get_parameter('jones_dole_viscosity')
#                    viscosity_rel += coefficients[0] * conc ** 0.5 + coefficients[1] * conc + \
#                    coefficients[2] * conc ** 2
#                except TypeError:
#                    continue
        viscosity_rel = self.get_viscosity_dynamic() / h2o.water_viscosity_dynamic(self.get_temperature(),self.get_pressure())    
        
        return viscosity_rel
        
    def get_viscosity_dynamic(self):
        '''
        Return the dynamic (absolute) viscosity of the solution.
        
        Calculated from the kinematic viscosity
    
        See Also
        --------
        get_viscosity_kinematic
        get_viscosity_relative
        '''
        return self.get_viscosity_kinematic() * self.get_density()    
    
    def get_viscosity_kinematic(self):
        '''
        Return the kinematic viscosity of the solution.
        
        Notes
        -----
        The calculation is based on a model derived from the Eyring equation
        and presented in [1]_
        
        .. math:: \\ln \\nu = \\ln {\\nu_w MW_w \over \sum_i x_i MW_i } + \
        15 x_+^2 + x_+^3  \delta G^*_{123} + 3 x_+ \delta G^*_{23} (1-0.05x_+)
        
        Where:
        
        .. math:: \delta G^*_{123} = a_o + a_1 (T)^{0.75}
        .. math:: \delta G^*_{23} = b_o + b_1 (T)^{0.5}
        
        In which `\\nu` is the kinematic viscosity, MW is the molecular weight,
        `x_+` is the mole fraction of cations, and T is the temperature in degrees C.
        
        The a and b fitting parameters for a variety of common salts are included in the
        database.        
        
        References
        ----------  
        .. [1] Vásquez-Castillo, G.; Iglesias-Silva, G. a.; Hall, K. R. An extension
        of the McAllister model to correlate kinematic viscosity of electrolyte solutions.
        Fluid Phase Equilib. 2013, 358, 44–49.
                
        See Also
        --------
        get_density_dynamic
        get_viscosity_relative
        '''
        # identify the main salt in the solution
        salt = self.get_salt()
        cation = salt.cation
        
        # search the database for parameters for 'salt'
        database.search_parameters(salt.formula)
        
        a0=a1=b0=b1 = 0
        
        found = False                
        
        # retrieve the parameters for the delta G equations
        for item in db[salt.formula]:

            if item.get_name() == 'erying_viscosity_coefficients':
                found = True
    
                a0 = item.get_value()[0]
                a1 = item.get_value()[1]
                b0 = item.get_value()[2]
                b1 = item.get_value()[3]

        # compute the delta G parameters
        temperature = self.get_temperature().to('degC')
        G_123 = a0 + a1 * (temperature.magnitude) ** 0.75
        G_23 = b0 + b1 * (temperature.magnitude) ** 0.5
        
        # get the kinematic viscosity of water
        nu_w = h2o.water_viscosity_kinematic(temperature,self.get_pressure()).to('m**2 / s').magnitude
        
        # compute the effective molar mass of the solution
        MW= self.get_mass() / (self.get_moles_solvent() + self.get_total_moles_solute())
        
        # get the MW of water
        MW_w = self.get_solvent().get_molecular_weight()
        
        # calculate the cation mole fraction
        x_cat = self.get_amount(cation,'fraction')
        
        # calculate the kinematic viscosity
        nu = math.log(nu_w * MW_w / MW) + 15 * x_cat ** 2 + x_cat ** 3 * G_123 + \
        3 * x_cat * G_23 * (1-0.05*x_cat)

        return math.exp(nu) * unit('m**2 / s')
        
    def get_conductivity(self):
        '''
        Compute the electrical conductivity of the solution.
        
        Parameters
        ----------
        None        
        
        Returns
        -------
        Quantity
            The electrical conductivity of the solution in Siemens / meter.
            
        Notes
        -----
        Conductivity is calculated by summing the molar conductivities of the respective
        solutes, but they are activity-corrected and adjusted using an empricial exponent.
        This approach is used in PHREEQC and Aqion models [1]_ [2]_
        
        .. math:: EC = {F^2 \\over R T} \\sum_i D_i z_i ^ 2 \\gamma_i ^ {\\alpha} m_i
        
        Where:
        
        .. math:: \\alpha = \\begin{cases} {0.6 \\over \\sqrt{|z_i|}} & {I < 0.36|z_i|} \\\ {\\sqrt{I} \\over |z_i|} & otherwise \\end{cases}
        
        Note: PHREEQC uses the molal rather than molar concentration according to
        <http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/phreeqc3-html/phreeqc3-43.htm>

        References
        ----------
        .. [1] <http://www.aqion.de/site/77>        
        .. [2] <http://www.hydrochemistry.eu/exmpls/sc.html>
                
        
        See Also
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
                
                diffusion_coefficient = self.get_property(item,'diffusion_coefficient')
        
                molar_cond = diffusion_coefficient * (unit.e * unit.N_A) ** 2 * self.get_solute(item).get_formal_charge() ** 2 / (unit.R * temperature)
                           
                EC += molar_cond * self.get_activity_coefficient(item) ** alpha * self.get_amount(item,'mol/L')
                             
        return EC.to('S/m')
    
    def get_osmotic_pressure(self):
        ''' 
        Return the osmotic pressure of the solution relative to pure water
        
        Parameters
        ----------
        None
        
        Returns
        -------
        Quantity
                The osmotic pressure of the solution relative to pure water in Pa
                
        Notes
        -----
        Osmotic pressure is calculated based on the water activity [1]_ [2]_ :
        
        .. math:: \\Pi = {RT \\over V_w} \ln a_w
        
        Where :math:`\\Pi` is the osmotic pressure, :math:`V_w` is the partial
        molar volume of water (18.2 cm**3/mol), and :math:`a_w` is the water
        activity.
                
        References
        ----------
        .. [1] Sata, Toshikatsu. Ion Exchange Membranes: Preparation, Characterization, and Modification. Royal Society of Chemistry, 2004, p. 10.
        
        .. [2] http://en.wikipedia.org/wiki/Osmotic_pressure#Derivation_of_osmotic_pressure
        
        Examples
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
    
    def p(self,solute,activity=True):
        ''' (number) -> float
        Return the negative log of the activity of solute.
        Generally used for expressing concentration of hydrogen ions (pH)
        
        Parameters
        ----------
        solute : str
            String representing the formula of the solute
        activity: bool, optional
            If False, the function will use the molar concentration rather 
            than the activity to calculate p. Defaults to True.
        
        Returns
        -------
        Quantity
            The negative log10 of the activity (or molar concentration if
            activity = False) of the solute.
            
        Examples
        --------

        
        '''
        if activity is True:
            return -1 * math.log10(self.get_activity(solute))
        elif activity is False:
            return -1 * math.log10(self.get_amount(solute,'mol/L').magnitude)

    def get_amount(self,solute,units):
        '''returns the amount of 'solute' in the parent solution
       
        Parameters
        ----------
        solute : str 
                    String representing the name of the solute of interest
        units : str
                    Units desired for the output. Examples of valid units are 
                    'mol/L','mol/kg','mol', 'kg', and 'g/L'
                    Use 'fraction' to return the mole fraction.

        Returns
        -------
        The amount of the solute in question, in the specified units
        
        
        See Also
        --------
        '''
        if units == 'fraction':
            return self.get_mole_fraction(solute)
        else:
            moles = self.get_solute(solute).get_moles()
            mw = self.get_solute(solute).get_molecular_weight()
        
        # with pint unit conversions enabled, we just pass the unit to pint
        # the logic tests here ensure that only the required arguments are 
        # passed to pint for the unit conversion. This avoids unecessary 
        # function calls.
        if unit(units).dimensionality == ('[substance]/[length]**3' or '[mass]/[length]**3'):
            return moles.to(units,'chem',mw=mw,volume=self.get_volume())
        elif unit(units).dimensionality == ('[substance]/[mass]' or '[mass]/[mass]'):
            return moles.to(units,'chem',mw=mw,solvent_mass=self.get_solvent_mass())
        elif unit(units).dimensionality == ('[mass]'):
            return moles.to(units,'chem',mw=mw)
        elif unit(units).dimensionality == ('[substance]'):
            return moles.to(units)
        else:
            logger.error('Unsupported unit specified for get_amount')
            return None

    def add_amount(self,solute,amount):
        '''Adds the amount of 'solute' to the parent solution.
       
        Parameters
        ----------
        solute : str 
                    String representing the name of the solute of interest
        amount : str quantity
                    String representing the concentration desired, e.g. '1 mol/kg'
                    If the units are given on a per-volume basis, the solution 
                    volume is not recalculated
                    If the units are given on a mass, substance, per-mass, or 
                    per-substance basis, then the solution volume is recalculated
                    based on the new composition

        Returns
        -------
        Nothing. The concentration of solute is modified.
        

        See Also
        --------
        Solute.add_moles()
        '''
        
        # if units are given on a per-volume basis, 
        # iteratively solve for the amount of solute that will preserve the
        # original volume and result in the desired concentration  
        if unit(amount).dimensionality == ('[substance]/[length]**3' or '[mass]/[length]**3'):
            
            # store the original volume for later
            orig_volume = self.get_volume()
            
            # change the amount of the solute present to match the desired amount
            self.get_solute(solute).add_moles(amount,self.get_volume(),self.get_solvent_mass())            
            
            # set the amount to zero and log a warning if the desired amount
            # change would result in a negative concentration
            if self.get_amount(solute,'mol').magnitude < 0:
                logger.warning('Attempted to set a negative concentration for solute %s. Concentration set to 0' % solute)
                self.set_amount(solute,'0 mol')
            
            # calculate the volume occupied by all the solutes
            solute_vol = self._get_solute_volume()
            
            # determine the volume of solvent that will preserve the original volume
            target_vol = orig_volume - solute_vol
            
            # adjust the amount of solvent
            target_mass = target_vol * h2o.water_density(self.get_temperature())
            mw = self.get_solvent().get_molecular_weight()
            target_mol = target_mass / mw
            self.get_solvent().moles = target_mol
            
        else:

            # change the amount of the solute present
            self.get_solute(solute).add_moles(amount,self.get_volume(),self.get_solvent_mass())
        
            # set the amount to zero and log a warning if the desired amount
            # change would result in a negative concentration
            if self.get_amount(solute,'mol').magnitude < 0:
                logger.warning('Attempted to set a negative concentration for solute %s. Concentration set to 0' % solute)
                self.set_amount(solute,'0 mol')
            
            # update the volume to account for the space occupied by all the solutes
            # make sure that there is still solvent present in the first place
            if self.get_solvent_mass() <= unit('0 kg'):
                logger.error('All solvent has been depleted from the solution')
                return None
            else:
                # set the volume recalculation flag
                self.volume_update_required = True

    def set_amount(self,solute,amount):
        '''Sets the amount of 'solute' in the parent solution.
       
        Parameters
        ----------
        solute : str 
                    String representing the name of the solute of interest
        amount : str Quantity
                    String representing the concentration desired, e.g. '1 mol/kg'
                    If the units are given on a per-volume basis, the solution 
                    volume is not recalculated and the molar concentrations of 
                    other components in the solution are not altered, while the
                    molal concentrations are modified.
                    
                    If the units are given on a mass, substance, per-mass, or 
                    per-substance basis, then the solution volume is recalculated
                    based on the new composition and the molal concentrations of
                    other components are not altered, while the molar concentrations
                    are modified.

        Returns
        -------
        Nothing. The concentration of solute is modified.
        

        See Also
        --------
        Solute.set_moles()
        '''
        # raise an error if a negative amount is specified
        if unit(amount).magnitude < 0:
            logger.error('Negative amount specified for solute %s. Concentration not changed.' % solute)
        
        # if units are given on a per-volume basis, 
        # iteratively solve for the amount of solute that will preserve the
        # original volume and result in the desired concentration  
        elif unit(amount).dimensionality == ('[substance]/[length]**3' or '[mass]/[length]**3'):
            
            # store the original volume for later
            orig_volume = self.get_volume()
            
            # change the amount of the solute present to match the desired amount
            self.get_solute(solute).set_moles(amount,self.get_volume(),self.get_solvent_mass())            

            # calculate the volume occupied by all the solutes
            solute_vol = self._get_solute_volume()
            
            # determine the volume of solvent that will preserve the original volume
            target_vol = orig_volume - solute_vol
            
            # adjust the amount of solvent
            target_mass = target_vol * h2o.water_density(self.get_temperature())
            mw = self.get_solvent().get_molecular_weight()
            target_mol = target_mass / mw
            self.get_solvent().moles = target_mol
            
        else:
            
            # change the amount of the solute present
            self.get_solute(solute).set_moles(amount,self.get_volume(),self.get_solvent_mass())

            # update the volume to account for the space occupied by all the solutes
            # make sure that there is still solvent present in the first place
            if self.get_solvent_mass() <= unit('0 kg'):
                logger.error('All solvent has been depleted from the solution')
                return None
            else:
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
        
        
        Parameters
        ----------
        solute : str 
                 String representing the name of the solute of interest
    
        Returns
        -------
        float
            The mole fraction of 'solute' in the parent Solution object
    
        See Also
        --------
        get_solvent_mass()
        
        Notes
        -----
        This function assumes water is the solvent with MW = 18
 
        Examples
        --------
        TBD
        
        '''
        return (self.get_amount(solute,'moles') / (self.get_moles_solvent() + self.get_total_moles_solute())).magnitude
    
    def get_moles_solvent(self):
        return self.get_amount(self.solvent_name,'mol')
    
    def get_salt(self):
        # identify the predominant salt in the solution
        import pyEQL.salt_ion_match as salt
        return salt.identify_salt(self)
      
## Activity-related methods
    def get_activity_coefficient(self,solute):
        '''Routine to determine the activity coefficient of a solute in solution. The correct function is chosen based on the ionic strength of the parent solution.
        
        Parameters
        ----------
        solute : str 
                    String representing the name of the solute of interest
        
        Returns
        -------
        The molal (mol/kg) scale mean ion activity coefficient of the solute in question
        
        See Also
        --------
        get_activity_coefficient_debyehuckel
        get_activity_coefficient_guntelberg
        get_activity_coefficient_davies
        get_activity_coefficient_pitzer
        get_activity_coefficient_TCPC
        '''
        ion = self.components[solute]
        temperature = str(self.get_temperature())
        
        # identify the predominant salt in the solution
        Salt = self.get_salt()
        
        # search the database for pitzer parameters for 'salt'
        database.search_parameters(Salt.formula)
    

        # use the Pitzer model for higher ionic strenght, if the parameters are available
        
        # search for Pitzer parameters    
        found = False
        for item in db[Salt.formula]:
            if item.get_name() == 'pitzer_parameters_activity':
                found = True
                break       
                
        if found == True:
            # determine alpha1 and alpha2 based on the type of salt
            # see the May reference for the rules used to determine
            # alpha1 and alpha2 based on charge
            if Salt.nu_cation >= 2 and Salt.nu_anion >=2:
                if Salt.nu_cation >=3 or Salt.nu_anion >=3:
                    alpha1 = 2
                    alpha2 = 50
                else:
                    alpha1 = 1.4
                    alpha2 = 12
            else:
                alpha1 = 2
                alpha2 = 0
            
            # determine the average molality of the salt
            # this is necessary for solutions inside e.g. an ion exchange
            # membrane, where the cation and anion concentrations may be
            # unequal
            molality = (self.get_amount(Salt.cation,'mol/kg')+self.get_amount(Salt.anion,'mol/kg'))/2
            
            activity_coefficient=ac.get_activity_coefficient_pitzer(self.get_ionic_strength(), \
            molality,alpha1,alpha2,item.get_value()[0],item.get_value()[1],item.get_value()[2],item.get_value()[3], \
            Salt.z_cation,Salt.z_anion,Salt.nu_cation,Salt.nu_anion,temperature)
            
            logger.info('Calculated activity coefficient of species %s as %s based on salt %s using Pitzer model' % (solute,activity_coefficient,Salt))
            return activity_coefficient            

        # for very low ionic strength, use the Debye-Huckel limiting law
        elif self.get_ionic_strength().magnitude <= 0.005:
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
            logger.warning('Ionic strength too high to estimate activity for species %s. Specify parameters for Pitzer or TCPC models. Returning unit activity coefficient' % solute)
            
            return unit('1 dimensionless')
                
            # TODO - fix TCPC implementation
            # use the TCPC model for higher ionic strengths, if the parameters have been set
#            if self.components[solute].parameters_TCPC:
#                logger.info('Ionic strength = %s. Using TCPC model to calculate activity coefficient.' % self.get_ionic_strength())
#                return ac.get_activity_coefficient_TCPC(self.get_ionic_strength(),ion.get_parameters_TCPC('S'),ion.get_parameters_TCPC('b'),ion.get_parameters_TCPC('n'),ion.get_parameters_TCPC('z_plus'),ion.get_parameters_TCPC('z_minus'),ion.get_parameters_TCPC('nu_plus'),ion.get_parameters_TCPC('nu_minus'),temperature)
#    
    def get_activity(self,solute):
        '''returns the thermodynamic activity of the solute in solution
       
        Parameters
        ----------
        solute : str 
                    String representing the name of the solute of interest
        temperature :    Quantity, optional
                     The temperature of the solution. Defaults to 25 degrees C if omitted
        
        Returns
        -------
        The thermodynamic activity of the solute in question (dimensionless)
        
        Notes
        -----
        The thermodynamic activity is independent of the concentration scale used. However,
        the concentration and the activity coefficient must use corresponding scales.[1]_ [2]_
        In this module, ionic strength, activity coefficients, and activities are all
        calculated based on the molal (mol/kg) concentration scale.
        
        References
        ----------
        .. [1] http://adsorption.org/awm/utils/Activity.htm
        .. [2] http://en.wikipedia.org/wiki/Thermodynamic_activity#Activity_coefficient
        
        See Also
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

        Returns
        -------
        Quantity : 
            The osmotic coefficient
            
        Notes
        -----
        For ionic strengths below 0.5 mol/kg, the osmotic coefficient is assumed to equal 1.0.
        1.0 will also be returned at higher ionic strengths if appropriate Pitzer
        parameters are not supplied.
        '''
        temperature = str(self.get_temperature())
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
                found = True
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
        
        Returns
        -------
        float : 
            The thermodynamic activity of water in the solution.
            
        Examples
        --------
        If 'soln' is a 0.5 mol/kg NaCl solution at 25 degC:
        >>> soln.get_water_activity()
        0.9835...
        
        If 'soln' is a 5.11 mol/kg NaHCO2 (sodium formate) solution at 25 degC:
        (literature value from Cabot specialty fluids is 0.82)
        >>> soln.get_water_activity()
        0.8631...
        
        Notes
        -----
        Water activity is related to the osmotic coefficient in a solution containing i solutes by:[1]_
        
        .. math:: \ln a_w = - \\Phi M_w \\sum_i m_i
        
        Where M_w is the molar mass of water (0.018015 kg/mol) and m_i is the molal concentration
        of each species.
        
        If appropriate Pitzer or TCPC model parameters are not available, the
        water activity is assumed equal to the mole fraction of water.
        
        References
        ----------
        .. [1] Blandamer, Mike J., Engberts, Jan B. F. N., Gleeson, Peter T., Reis, Joao Carlos R., 2005. "Activity of water in aqueous systems: A frequently neglected property."
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
            
            return math.exp(- osmotic_coefficient * 0.018015*unit('kg/mol') * concentration_sum)
    
    def get_ionic_strength(self):
        '''() -> float
        
        Return the ionic strength of the solution, calculated as 1/2 * sum ( molality * charge ^2) over all the ions.
        Molal (mol/kg) scale concentrations are used for compatibility with the activity correction formulas, but
        the returned value does not carry units.
        
        
        Returns
        -------
        float : 
            The ionic strength of the parent solution, mol/kg.
        
        Examples
        --------
        TODO
#         >>> conc_soln.list_concentrations()
#         {'Na+': 5.999375074924214, 'Cl-': 5.999904143046362, 'HCO3-': 0, 'NaCO3-': 0, 'NaHCO3': 0}
#         >>> conc_soln.get_ionic_strength()
#         5.999639608985288
        
        Notes
        -----
        The ionic strength is calculated according to:
        
        .. math:: I = \sum_i m_i z_i^2
        
        Where m_i is the molal concentration and z_i is the charge on species i.
        
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
    
    def get_transport_number(self,solute,activity_correction = False):
        '''Calculate the transport number of the solute in the solution
        
        Parameters
        ----------
        solute : str
            String identifying the solute for which the transport number is
            to be calculated.
        activity_correction: bool
            If True, the transport number will be corrected for activity following
            the same method used for solution conductivity. Defaults to False
            if omitted.
        
        Returns
        -------
        float
            The transport number of `solute`
        
        Notes
        -----
        Transport number is calculated according to:[1]_
        
        .. math:: t_i = {D_i z_i^2 C_i \\over \sum D_i z_i^2 C_i}
        
        Where C is the concentration in mol/L.
        
        If `activity_correction` is True, the contribution of each ion to the
        transport number is corrected with an activity factor. See the documentation
        for get_conductivity() for an explanation of this correction.
        
        References
        ----------
        .. [1] Geise, G. M.; Cassady, H. J.; Paul, D. R.; Logan, E.; Hickner, M. A. Specific ion effects on membrane potential and the permselectivity of ion exchange membranes. Phys. Chem. Chem. Phys. 2014, 16, 21673–21681.

        See Also
        --------
        get_conductivity()
        '''
        denominator= 0
        numerator = 0
        
        for item in self.components:
            
            if item == 'H2O':
                continue
            
            z = self.get_solute(item).get_formal_charge()
            term = self.get_property(item,'diffusion_coefficient') * \
            z ** 2 * self.get_amount(item,'mol/L')
            
            if activity_correction == True:
                gamma = self.get_activity_coefficient(item)
                
                if self.get_ionic_strength().magnitude < 0.36 * z:
                    alpha = 0.6 / z ** 0.5
                else:
                    alpha = self.get_ionic_strength().magnitude ** 0.5 / z   
                
                if item == solute:
                    numerator = term * gamma ** alpha
            
                denominator += term  * gamma ** alpha
                        
            else:
                if item == solute:
                    numerator = term
            
                denominator += term
        
        return numerator / denominator
    
    def get_property(self,solute,name):
        '''Retrieve a thermodynamic property (such as diffusion coefficient)
        for solute, and adjust it from the reference conditions to the conditions
        of the solution
        
        Parameters
        ----------
        solute: str
            String representing the chemical formula of the solute species
        name: str
            The name of the property needed, e.g.
            'diffusion coefficient'
        
        Returns
        -------
        Quantity: The desired parameter
        
        '''
        # retrieve the base value and the conditions of measurement from the
        # database
        
        if solute != 'H2O':
            if database.has_parameter(solute,name):
                base_value = self.get_solute(solute).get_parameter(name)
            else:
                base_value = None                
                
        base_temperature = unit('25 degC')
        base_pressure = unit ('1 atm')
        
        # perform temperature-corrections or other adjustments for certain
        # parameter types        
        if name == 'diffusion_coefficient':
            if base_value is not None:
                # correct for temperature and viscosity
                # $$ D_1 \over D_2 = T_1 \over T_2 * \mu_2 \over \mu_1 $$
                # where $\mu$ is the dynamic viscosity
                # assume that the base viscosity is that of pure water
                return base_value * self.get_temperature() / base_temperature * h2o.water_viscosity_dynamic(base_temperature,base_pressure) / self.get_viscosity_dynamic()
            else:
                logger.warning('Diffusion coefficient not found for species %s. Assuming zero.' % (solute))
                return unit('0 m**2/s')
            
        # just return the base-value molar volume for now; find a way to adjust for
        # concentration later
        if name == 'partial_molar_volume':
            # calculate the partial molar volume for water since it isn't in the database            
            if solute == 'H2O':
                vol = self.get_solute('H2O').get_molecular_weight() / h2o.water_density(self.get_temperature())
                return vol.to('cm **3 / mol')
            else:
                if base_value is not None:
                    return base_value
                else:
                    logger.warning('Partial molar volume not found for species %s. Assuming zero.' % (solute))
                    return unit ('0 cm **3 / mol')
                
    def get_chemical_potential_energy(self,activity_correction=True):
        '''
        Return the total chemical potential energy of a solution (not including
        pressure or electric effects)
    
        Parameters
        ----------
        activity_correction : bool, optional
            If True, activities will be used to calculate the true chemical 
            potential. If False, mole fraction will be used, resulting in 
            a calculation of the ideal chemical potential.
            
        Returns
        -------
        Quantity
            The actual or ideal chemical potential energy of the solution, in Joules.
        
        Notes
        -----
        
        The chemical potential energy (related to the Gibbs mixing energy) is
        calculated as follows: [1]_
            
        .. math::
            E = R T \sum_i n_i  \ln a_i
            
            or 
            
            E = R T \sum_i n_i \ln x_i
        
        Where n is the number of moles of substance, T is the temperature in kelvin,
        R the ideal gas constant, x the mole fraction, and a the activity of
        each component.
        
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
        temperature = self.get_temperature()
        
        E = unit('0 J')
        
        # loop through all the components and add their potential energy
        for item in self.components:
            if activity_correction is True:
                E += unit.R * temperature.to('K') * self.get_amount(item,'mol') * math.log(self.get_activity(item))
            else:
                E += unit.R * temperature.to('K') * self.get_amount(item,'mol') * math.log(self.get_amount(item,'fraction'))
    
        return E.to('J')

    def get_lattice_distance(self,solute):
        '''Calculate the average distance between molecules of the given solute,
        assuming that the molecules are uniformly distributed throughout the
        solution.
        
        Parameters
        ----------
        solute : str 
                    String representing the name of the solute of interest
        
        Returns
        -------
        Quantity : The average distance between solute molecules
        
        Examples
        --------
        >>> soln = Solution([['Na+','0.5 mol/kg'],['Cl-','0.5 mol/kg']])
        >>> soln.get_lattice_distance('Na+')
        1.492964.... nanometer
        
        Notes
        -----
        The lattice distance is related to the molar concentration as follows:
        
        .. math:: d = ( C_i N_A ) ^ {-{1\over3}}
        
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
        self.volume = self._get_solvent_volume() + self._get_solute_volume()
        
    def _get_solvent_volume(self):
        '''
        Return the volume of the pure solvent
        
        '''        
        # calculate the volume of the pure solvent
        solvent_vol = self.get_solvent_mass() / h2o.water_density(self.get_temperature())
        
        return solvent_vol.to('L')
        
    def _get_solute_volume(self):
        '''
        Return the volume of only the solutes
        
        '''   
        temperature = str(self.get_temperature())
        
        # identify the predominant salt in the solution
        Salt = self.get_salt()
        
        # search the database for pitzer parameters for 'salt'
        database.search_parameters(Salt.formula)
         
        solute_vol = 0 * unit('L')

        # use the pitzer approach if parameters are available
        
        pitzer_calc = False

        if database.has_parameter(Salt.formula,'pitzer_parameters_volume'):
            
            for params in db[Salt.formula]:
                if params.get_name() == 'pitzer_parameters_volume':
                    
                    # determine the average molality of the salt
                    # this is necessary for solutions inside e.g. an ion exchange
                    # membrane, where the cation and anion concentrations may be
                    # unequal
                    molality = (self.get_amount(Salt.cation,'mol/kg')+self.get_amount(Salt.anion,'mol/kg'))/2
                    
                    # determine alpha1 and alpha2 based on the type of salt
                    # see the May reference for the rules used to determine
                    # alpha1 and alpha2 based on charge
                    if Salt.nu_cation >= 2 and Salt.nu_anion >=2:
                        if Salt.nu_cation >=3 or Salt.nu_anion >=3:
                            alpha1 = 2
                            alpha2 = 50
                        else:
                            alpha1 = 1.4
                            alpha2 = 12
                    else:
                        alpha1 = 2
                        alpha2 = 0
                        
                    apparent_vol = ac.get_apparent_volume_pitzer(self.get_ionic_strength(), \
                    molality,alpha1,alpha2,params.get_value()[0],params.get_value()[1],params.get_value()[2],params.get_value()[3], \
                    params.get_value()[4],Salt.z_cation,Salt.z_anion,Salt.nu_cation,Salt.nu_anion,temperature)
            
            solute_vol += apparent_vol * (self.get_amount(Salt.cation,'mol')/Salt.nu_cation \
            +self.get_amount(Salt.anion,'mol')/Salt.nu_anion)/2
            
            pitzer_calc = True
            
            logger.info('Updated solution volume using Pitzer model for solute %s' % Salt.formula)
            
            # add the partial molar volume of any other solutes, except for water
            # which is already accounted for by the Pitzer parameters
            for item in self.components:
                
                solute = self.get_solute(item)            
                
                # ignore water
                if item in ['H2O','HOH']:
                    continue
                
                # ignore the salt cation and anion, if already accounted for by Pitzer
                if pitzer_calc is True and item in [Salt.anion,Salt.cation]:
                    continue                
                
                if database.has_parameter(item,'partial_molar_volume'):
                    solute_vol += solute.get_parameter('partial_molar_volume') * solute.get_moles()
                    logger.info('Updated solution volume using direct partial molar volume for solute %s' % item)
                    
                else:
                    logger.warning('Partial molar volume data not available for solute %s. Solution volume will not be corrected.' % item)
                    
        return solute_vol.to('L')
            
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
            self.act_list.update({i:str(self.get_activity(i))})
        print('Component activities:\n',self.act_list )
     
    def __str__(self):
        #set output of the print() statement for the solution     
        return 'Components: \n'+str(self.list_solutes()) + '\n' + 'Volume: '+str(self.get_volume()) + '\n' + 'Density: '+str(self.get_density())
