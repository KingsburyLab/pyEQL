"""
pyEQL Solution Class

:copyright: 2013-2022 by Ryan S. Kingsbury
:license: LGPL, see LICENSE for more details.

"""

# logging system
import logging

# import libraries for scientific functions
import math
from typing import Literal, Optional, Union, List, Dict

from iapws import IAPWS95
from monty.dev import deprecated
from pint import DimensionalityError

# internal pyEQL imports
import pyEQL.solute as sol

# import the parameters database
# the pint unit registry
from pyEQL import paramsDB as db
from pyEQL import unit
from pyEQL.engines import EOS, IdealEOS, NativeEOS

# add a filter to emit only unique log messages to the handler
from pyEQL.logging_system import Unique
from pyEQL.salt_ion_match import generate_salt_list, identify_salt

logger = logging.getLogger(__name__)
unique = Unique()
logger.addFilter(unique)

# add a handler for console output, since pyEQL is meant to be used interactively
ch = logging.StreamHandler()

# create formatter for the log
formatter = logging.Formatter("(%(name)s) - %(levelname)s - %(message)s")

# add formatter to the handler
ch.setFormatter(formatter)
logger.addHandler(ch)


class Solution:
    """
    Class representing the properties of a solution. Instances of this class
    contain information about the solutes, solvent, and bulk properties.
    """

    def __init__(
        self,
        solutes: Optional[Union[List[List[str]], Dict[str, str]]] = None,
        volume: Optional[str] = None,
        temperature: str = "298.15 K",
        pressure: str = "1 atm",
        pH: float = 7,
        engine: Literal["native", "ideal"] = "native",
        **kwargs,
    ):
        """

        Args:
            solutes : dict, optional. Keys must be the chemical formula, while values must be
                        str Quantity representing the amount. For example:

                        {"Na+": "0.1 mol/L", "Cl-": "0.1 mol/L"}

                        Note that an older "list of lists" syntax is also supported; however this
                        will be deprecated in the future and is no longer recommended. The equivalent
                        list syntax for the above example is

                        [["Na+", "0.1 mol/L"], ["Cl-", "0.1 mol/L"]]

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

        Examples:
            >>> s1 = pyEQL.Solution([['Na+','1 mol/L'],['Cl-','1 mol/L']],temperature='20 degC',volume='500 mL')
            >>> print(s1)
            Components:
            ['H2O', 'Cl-', 'H+', 'OH-', 'Na+']
            Volume: 0.5 l
            Density: 1.0383030844030992 kg/l

        See Also:
            add_solute
        """

        # initialize the volume with a flag to distinguish user-specified volume
        if volume is not None:
            volume_set = True
            self.volume = unit(volume).to("L")
        else:
            volume_set = False
            self.volume = unit("1 L")
        self._temperature = unit(temperature)
        self._pressure = unit(pressure)

        # instantiate a water substance for property retrieval
        self.water_substance = IAPWS95(
            T=self.temperature.magnitude,
            P=self.pressure.to("MPa").magnitude,
        )

        # create an empty dictionary of components
        self.components: dict = {}

        # initialize the volume recalculation flag
        self.volume_update_required = False

        # set the equation of state engine
        # self.engine: Optional[EOS] = None
        if engine == "ideal":
            self.engine: EOS = IdealEOS()
        elif engine == "native":
            self.engine = NativeEOS()
        else:
            raise ValueError(f'{engine} is not a valid value for the "engine" kwarg!')

        # define the solvent
        if "solvent" in kwargs:
            self.solvent_name = kwargs["solvent"][0]
            # warn if the solvent is anything besides water
            if not kwargs["solvent"][0] == "H2O" or kwargs["solvent"][0] == "water":
                logger.error(
                    "Non-aqueous solvent detected. These are not yet supported!"
                )

            # raise an error if the solvent volume has also been given
            if volume_set is True:
                logger.error(
                    "Solvent volume and mass cannot both be specified. Calculating volume based on solvent mass."
                )

            # add the solvent and the mass
            self.add_solvent(self.solvent_name, kwargs["solvent"][1])
        else:
            self.solvent_name = "H2O"

            # calculate the solvent (water) mass based on the density and the solution volume
            self.add_solvent(
                self.solvent_name,
                str(
                    self.volume.magnitude
                    / 1000
                    * self.water_substance.rho
                    * unit.Quantity("1 kg")
                ),
            )

        # set the pH with H+ and OH-
        self.add_solute("H+", str(10 ** (-1 * pH)) + "mol/L")
        self.add_solute("OH-", str(10 ** (-1 * (14 - pH))) + "mol/L")

        # populate the other solutes
        if isinstance(solutes, dict):
            for k, v in solutes.items():
                self.add_solute(k, v)
        elif isinstance(solutes, list):
            for item in solutes:
                self.add_solute(*item)
        elif solutes is not None:
            raise ValueError("Solutes must be given as a list or dict!")

    def add_solute(self, formula, amount, parameters={}):
        """Primary method for adding substances to a pyEQL solution

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

        """

        # if units are given on a per-volume basis,
        # iteratively solve for the amount of solute that will preserve the
        # original volume and result in the desired concentration
        if unit(amount).dimensionality in (
            "[substance]/[length]**3",
            "[mass]/[length]**3",
        ):

            # store the original volume for later
            orig_volume = self.get_volume()

            # add the new solute
            new_solute = sol.Solute(
                formula, amount, self.get_volume(), self.get_solvent_mass(), parameters
            )
            self.components.update({new_solute.get_name(): new_solute})

            # calculate the volume occupied by all the solutes
            solute_vol = self._get_solute_volume()

            # determine the volume of solvent that will preserve the original volume
            target_vol = orig_volume - solute_vol

            # adjust the amount of solvent
            # density is returned in kg/m3 = g/L
            target_mass = (
                target_vol.magnitude * self.water_substance.rho * unit.Quantity("1 g")
            )
            mw = self.get_solvent().get_molecular_weight()
            target_mol = target_mass / mw
            self.get_solvent().moles = target_mol

        else:

            # add the new solute
            new_solute = sol.Solute(
                formula, amount, self.get_volume(), self.get_solvent_mass(), parameters
            )
            self.components.update({new_solute.get_name(): new_solute})

            # update the volume to account for the space occupied by all the solutes
            # make sure that there is still solvent present in the first place
            if self.get_solvent_mass() <= unit("0 kg"):
                logger.error("All solvent has been depleted from the solution")
                return None
            else:
                # set the volume recalculation flag
                self.volume_update_required = True

    def add_solvent(self, formula, amount):
        """
        Same as add_solute but omits the need to pass solvent mass to pint
        """
        new_solvent = sol.Solute(formula, amount, self.get_volume(), amount)
        self.components.update({new_solvent.get_name(): new_solvent})

    def get_solute(self, i):
        """
        Return the specified solute object.

        """
        return self.components[i]

    def get_solvent(self):
        """
        Return the solvent object.

        """
        return self.components[self.solvent_name]

    @property
    def temperature(self):
        """
        Return the temperature of the solution in Kelvin.
        """
        return self._temperature.to("K")

    @temperature.setter
    def temperature(self, temperature: str):
        """
        Set the solution temperature.

        Args:
            temperature: pint-compatible string, e.g. '25 degC'
        """
        self._temperature = unit(temperature)
        # recalculate the volume
        self._update_volume()

    @deprecated(
        message="get_temperature() will be removed in the next release. Access the temperature directly via the attribute Solution.temperature"
    )
    def get_temperature(self):
        """
        Return the temperature of the solution.

        Parameters
        ----------
        None

        Returns
        -------
        Quantity: The temperature of the solution, in Kelvin.
        """
        return self.temperature

    @deprecated(
        message="set_temperature() will be removed in the next release. Set the temperature directly via the attribute Solution.temperature"
    )
    def set_temperature(self, temperature):
        """
        Set the solution temperature.

        Parameters
        ----------
        temperature : str
            String representing the temperature, e.g. '25 degC'
        """
        self.temperature = unit(temperature)

        # recalculate the volume
        self._update_volume()

    @property
    def pressure(self):
        """
        Return the hydrostatic pressure of the solution in atm.
        """
        return self._pressure.to("atm")

    @pressure.setter
    def pressure(self, pressure: str):
        """
        Set the solution pressure.

        Args:
            pressure: pint-compatible string, e.g. '1.2 atmC'
        """
        self._pressure = unit(pressure)
        # recalculate the volume
        self._update_volume()

    @deprecated(
        message="get_pressure() will be removed in the next release. Access the pressure directly via the attribute Solution.pressure"
    )
    def get_pressure(self):
        """
        Return the hydrostatic pressure of the solution.

        Returns
        -------
        Quantity: The hydrostatic pressure of the solution, in atm.
        """
        return self.pressure

    @deprecated(
        message="set_pressure() will be removed in the next release. Set the pressure directly via Solution.pressure"
    )
    def set_pressure(self, pressure):
        """
        Set the hydrostatic pressure of the solution.

        Parameters
        ----------
        pressure : str
            String representing the temperature, e.g. '25 degC'
        """
        self._pressure = unit(pressure)

        # recalculate the volume
        self._update_volume()

    def get_solvent_mass(self):
        """
        Return the mass of the solvent.

        This method is used whenever mol/kg (or similar) concentrations
        are requested by get_amount()

        Parameters
        ----------
        None

        Returns
        -------
        Quantity: the mass of the solvent, in kg

        See Also
        --------
        get_amount()
        """
        # return the total mass (kg) of the solvent
        solvent = self.get_solvent()
        mw = solvent.get_molecular_weight()

        return solvent.get_moles().to("kg", "chem", mw=mw)

    def get_volume(self):
        """
        Return the volume of the solution.

        Parameters
        ----------
        None

        Returns
        -------
        Quantity: the volume of the solution, in L
        """

        # if the composition has changed, recalculate the volume first
        if self.volume_update_required is True:
            self._update_volume()
            self.volume_update_required = False

        return self.volume.to("L")

    def set_volume(self, volume):
        """Change the total solution volume to volume, while preserving
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

        """
        # figure out the factor to multiply the old concentrations by
        scale_factor = unit(volume) / self.get_volume()

        # scale down the amount of all the solutes according to the factor
        for item in self.components:
            self.get_solute(item).moles = self.get_solute(item).moles * scale_factor

        # update the solution volume
        self.volume = unit(volume)

    def get_mass(self):
        """
        Return the total mass of the solution.

        The mass is calculated each time this method is called.
        Parameters
        ----------
        None

        Returns
        -------
        Quantity: the mass of the solution, in kg

        """
        total_mass = 0
        for item in self.components:
            total_mass += self.get_amount(item, "kg")
        return total_mass.to("kg")

    def get_density(self):
        """
        Return the density of the solution.

        Density is calculated from the mass and volume each time this method is called.

        Returns
        -------
        Quantity: The density of the solution.
        """
        return self.get_mass() / self.get_volume()

    def get_dielectric_constant(self):
        """
        Returns the dielectric constant of the solution.

        Parameters
        ----------
        None

        Returns
        -------
        Quantity: the dielectric constant of the solution, dimensionless.

        Notes
        -----
        Implements the following equation as given by [#]_

        .. math:: \\epsilon = \\epsilon_solvent \\over 1 + \\sum_i \\alpha_i x_i

        where :math:`\\alpha_i` is a coefficient specific to the solvent and ion, and :math:`x_i`
        is the mole fraction of the ion in solution.


        References
        ----------
        .. [#] [1] A. Zuber, L. Cardozo-Filho, V.F. Cabral, R.F. Checoni, M. Castier,
        An empirical equation for the dielectric constant in aqueous and nonaqueous
        electrolyte mixtures, Fluid Phase Equilib. 376 (2014) 116–123.
        doi:10.1016/j.fluid.2014.05.037.
        """
        di_water = self.water_substance.epsilon

        denominator = 1
        for item in self.components:
            # ignore water
            if item != "H2O":
                # skip over solutes that don't have parameters
                try:
                    fraction = self.get_amount(item, "fraction")
                    coefficient = self.get_solute(item).get_parameter(
                        "dielectric_parameter_water"
                    )
                    denominator += coefficient * fraction
                except TypeError:
                    logger.warning(
                        "No dielectric parameters found for species %s." % item
                    )
                    continue

        dielectric_constant = di_water / denominator

        return dielectric_constant

    def get_viscosity_relative(self):
        """
        Return the viscosity of the solution relative to that of water

        This is calculated using a simplified form of the Jones-Dole equation:

        .. math:: \\eta_{rel} = 1 + \\sum_i B_i m_i

        Where :math:`m` is the molal concentration and :math:`B` is an empirical parameter.

        See
        <http://downloads.olisystems.com/ResourceCD/TransportProperties/Viscosity-Aqueous.pdf>
        <http://www.nrcresearchpress.com/doi/pdf/10.1139/v77-148>
        <http://apple.csgi.unifi.it/~fratini/chen/pdf/14.pdf>
        """
        # if self.get_ionic_strength().magnitude > 0.2:
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
        viscosity_rel = (
            self.get_viscosity_dynamic()
            / self.water_substance.mu
            * unit.Quantity("1 Pa*s")
        )

        return viscosity_rel

    def get_viscosity_dynamic(self):
        """
        Return the dynamic (absolute) viscosity of the solution.

        Calculated from the kinematic viscosity

        See Also
        --------
        get_viscosity_kinematic
        get_viscosity_relative
        """
        return self.get_viscosity_kinematic() * self.get_density()

    def get_viscosity_kinematic(self):
        """
        Return the kinematic viscosity of the solution.

        Notes
        -----
        The calculation is based on a model derived from the Eyring equation
        and presented in [#]_

        .. math::

            \\ln \\nu = \\ln {\\nu_w MW_w \\over \\sum_i x_i MW_i } +
            15 x_+^2 + x_+^3  \\delta G^*_{123} + 3 x_+ \\delta G^*_{23} (1-0.05x_+)

        Where:

        .. math:: \\delta G^*_{123} = a_o + a_1 (T)^{0.75}
        .. math:: \\delta G^*_{23} = b_o + b_1 (T)^{0.5}

        In which :math: `\\nu` is the kinematic viscosity, MW is the molecular weight,
        `x_+` is the mole fraction of cations, and T is the temperature in degrees C.

        The a and b fitting parameters for a variety of common salts are included in the
        database.

        References
        ----------
        .. [#] Vásquez-Castillo, G.; Iglesias-Silva, G. a.; Hall, K. R. An extension
               of the McAllister model to correlate kinematic viscosity of electrolyte solutions.
               Fluid Phase Equilib. 2013, 358, 44–49.

        See Also
        --------
        get_viscosity_dynamic
        get_viscosity_relative
        """
        # identify the main salt in the solution
        salt = self.get_salt()
        cation = salt.cation

        # search the database for parameters for 'salt'
        db.search_parameters(salt.formula)

        a0 = a1 = b0 = b1 = 0

        # retrieve the parameters for the delta G equations
        if db.has_parameter(salt.formula, "erying_viscosity_coefficients"):
            params = db.get_parameter(salt.formula, "erying_viscosity_coefficients")

            a0 = params.get_value()[0]
            a1 = params.get_value()[1]
            b0 = params.get_value()[2]
            b1 = params.get_value()[3]
        else:
            # proceed with the coefficients equal to zero and log a warning
            logger.warning(
                "Viscosity coefficients for %s not found. Viscosity will be approximate."
                % salt.formula
            )

        # compute the delta G parameters
        temperature = self.temperature.to("degC")
        G_123 = a0 + a1 * (temperature.magnitude) ** 0.75
        G_23 = b0 + b1 * (temperature.magnitude) ** 0.5

        # get the kinematic viscosity of water, returned by IAPWS in m2/s
        nu_w = self.water_substance.nu

        # compute the effective molar mass of the solution
        MW = self.get_mass() / (
            self.get_moles_solvent() + self.get_total_moles_solute()
        )

        # get the MW of water
        MW_w = self.get_solvent().get_molecular_weight()

        # calculate the cation mole fraction
        x_cat = self.get_amount(cation, "fraction")

        # calculate the kinematic viscosity
        nu = (
            math.log(nu_w * MW_w / MW)
            + 15 * x_cat**2
            + x_cat**3 * G_123
            + 3 * x_cat * G_23 * (1 - 0.05 * x_cat)
        )

        return math.exp(nu) * unit("m**2 / s")

    def get_conductivity(self):
        """
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
        This approach is used in PHREEQC and Aqion models [#]_ [#]_

        .. math::

            EC = {F^2 \\over R T} \\sum_i D_i z_i ^ 2 \\gamma_i ^ {\\alpha} m_i

        Where:

        .. math::

            \\alpha = \\begin{cases} {0.6 \\over \\sqrt{|z_i|}} & {I < 0.36|z_i|} \\ {\\sqrt{I} \\over |z_i|} & otherwise \\end{cases}

        Note: PHREEQC uses the molal rather than molar concentration according to
        http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/phreeqc3-html/phreeqc3-43.htm

        References
        ----------
        .. [#] http://www.aqion.de/site/77
        .. [#] http://www.hydrochemistry.eu/exmpls/sc.html

        See Also
        --------
        get_ionic_strength()
        get_molar_conductivity()
        get_activity_coefficient()

        """
        EC = 0 * unit("S/m")

        for item in self.components:
            z = abs(self.get_solute(item).get_formal_charge())
            # ignore uncharged species
            if not z == 0:
                # determine the value of the exponent alpha
                if self.get_ionic_strength().magnitude < 0.36 * z:
                    alpha = 0.6 / z**0.5
                else:
                    alpha = self.get_ionic_strength().magnitude ** 0.5 / z

                diffusion_coefficient = self.get_property(item, "diffusion_coefficient")

                molar_cond = (
                    diffusion_coefficient
                    * (unit.e * unit.N_A) ** 2
                    * self.get_solute(item).get_formal_charge() ** 2
                    / (unit.R * self.temperature)
                )

                EC += (
                    molar_cond
                    * self.get_activity_coefficient(item) ** alpha
                    * self.get_amount(item, "mol/L")
                )

        return EC.to("S/m")

    def get_osmotic_pressure(self):
        """
        Return the osmotic pressure of the solution relative to pure water.

        Parameters
        ----------
        None

        Returns
        -------
        Quantity
                The osmotic pressure of the solution relative to pure water in Pa

        See Also
        --------
        get_water_activity
        get_osmotic_coefficient
        get_salt

        Notes
        -----
        Osmotic pressure is calculated based on the water activity [#]_ [#]_ :

        .. math:: \\Pi = {RT \\over V_w} \\ln a_w

        Where :math:`\\Pi` is the osmotic pressure, :math:`V_w` is the partial
        molar volume of water (18.2 cm**3/mol), and :math:`a_w` is the water
        activity.


        References
        ----------
        .. [#] Sata, Toshikatsu. Ion Exchange Membranes: Preparation, Characterization, and Modification. Royal Society of Chemistry, 2004, p. 10.

        .. [#] http://en.wikipedia.org/wiki/Osmotic_pressure#Derivation_of_osmotic_pressure

        Examples
        --------
        >>> s1=pyEQL.Solution()
        >>> s1.get_osmotic_pressure()
        0.0

        >>> s1 = pyEQL.Solution([['Na+','0.2 mol/kg'],['Cl-','0.2 mol/kg']])
        >>> soln.get_osmotic_pressure()
        <Quantity(906516.7318131207, 'pascal')>
        """
        # TODO - tie this into parameter() and solvent() objects
        partial_molar_volume_water = 1.82e-5 * unit("m ** 3/mol")

        osmotic_pressure = (
            -1
            * unit.R
            * self.temperature
            / partial_molar_volume_water
            * math.log(self.get_water_activity())
        )
        logger.info(
            "Computed osmotic pressure of solution as %s Pa at T= %s degrees C"
            % (osmotic_pressure, self.temperature)
        )
        return osmotic_pressure.to("Pa")

    # Concentration  Methods

    def p(self, solute, activity=True):
        """
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
        TODO

        """
        try:
            if activity is True:
                return -1 * math.log10(self.get_activity(solute))
            elif activity is False:
                return -1 * math.log10(self.get_amount(solute, "mol/L").magnitude)
        # if the solute has zero concentration, the log will generate a ValueError
        except ValueError:
            return 0

    def get_amount(self, solute, units):
        """
        Return the amount of 'solute' in the parent solution.

        The amount of a solute can be given in a variety of unit types.
        1. substance per volume (e.g., 'mol/L')
        1. substance per mass of solvent (e.g., 'mol/kg')
        1. mass of substance (e.g., 'kg')
        1. moles of substance ('mol')
        1. mole fraction ('fraction')
        1. percent by weight (%)

        Parameters
        ----------
        solute : str
                    String representing the name of the solute of interest
        units : str
                    Units desired for the output. Examples of valid units are
                    'mol/L','mol/kg','mol', 'kg', and 'g/L'
                    Use 'fraction' to return the mole fraction.
                    Use '%' to return the mass percent

        Returns
        -------
        The amount of the solute in question, in the specified units


        See Also
        --------
        add_amount
        set_amount
        get_total_amount
        get_osmolarity
        get_osmolality
        get_solvent_mass
        get_mass
        get_total_moles_solute
        """
        # retrieve the number of moles of solute and its molecular weight
        try:
            moles = self.get_solute(solute).get_moles()
            mw = self.get_solute(solute).get_molecular_weight()
        # if the solute is not present in the solution, we'll get a KeyError
        # In that case, the amount is zero
        except KeyError:
            try:
                return 0 * unit(units)
            except DimensionalityError:
                logger.error("Unsupported unit specified for get_amount")
                return 0

        # with pint unit conversions enabled, we just pass the unit to pint
        # the logic tests here ensure that only the required arguments are
        # passed to pint for the unit conversion. This avoids unecessary
        # function calls.
        if units == "fraction":
            return moles / (self.get_moles_solvent() + self.get_total_moles_solute())
        elif units == "%":
            return moles.to("kg", "chem", mw=mw) / self.get_mass().to("kg") * 100
        elif unit(units).dimensionality in (
            "[substance]/[length]**3",
            "[mass]/[length]**3",
        ):
            return moles.to(units, "chem", mw=mw, volume=self.get_volume())
        elif unit(units).dimensionality in ("[substance]/[mass]", "[mass]/[mass]"):
            return moles.to(units, "chem", mw=mw, solvent_mass=self.get_solvent_mass())
        elif unit(units).dimensionality == "[mass]":
            return moles.to(units, "chem", mw=mw)
        elif unit(units).dimensionality == "[substance]":
            return moles.to(units)
        else:
            logger.error("Unsupported unit specified for get_amount")
            return None

    def get_total_amount(self, element, units):
        """
        Return the total amount of 'element' (across all solutes) in the solution.

        Parameters
        ----------
        element : str
                    String representing the name of the element of interest
        units : str
                    Units desired for the output. Examples of valid units are
                    'mol/L','mol/kg','mol', 'kg', and 'g/L'

        Returns
        -------
        The total amount of the element in the solution, in the specified units

        Notes
        -----
        There is currently no way to distinguish between different oxidation
        states of the same element (e.g. TOTFe(II) vs. TOTFe(III)). This
        is planned for a future release. (TODO)

        See Also
        --------
        get_amount
        """
        import pyEQL.chemical_formula as ch

        TOT = 0 * unit(units)

        # loop through all the solutes, process each one containing element
        for item in self.components:
            # check whether the solute contains the element
            if ch.contains(item, element):
                # start with the amount of the solute in the desired units
                amt = self.get_amount(item, units)

                # convert the solute amount into the amount of element by
                # either the mole / mole or weight ratio
                if unit(units).dimensionality in (
                    "[substance]",
                    "[substance]/[length]**3",
                    "[substance]/[mass]",
                ):
                    TOT += amt * ch.get_element_mole_ratio(item, element)

                elif unit(units).dimensionality in (
                    "[mass]",
                    "[mass]/[length]**3",
                    "[mass]/[mass]",
                ):
                    TOT += amt * ch.get_element_weight_fraction(item, element)

        return TOT

    def add_amount(self, solute, amount):
        """
        Add the amount of 'solute' to the parent solution.

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
        """

        # if units are given on a per-volume basis,
        # iteratively solve for the amount of solute that will preserve the
        # original volume and result in the desired concentration
        if unit(amount).dimensionality in (
            "[substance]/[length]**3",
            "[mass]/[length]**3",
        ):

            # store the original volume for later
            orig_volume = self.get_volume()

            # change the amount of the solute present to match the desired amount
            self.get_solute(solute).add_moles(
                amount, self.get_volume(), self.get_solvent_mass()
            )

            # set the amount to zero and log a warning if the desired amount
            # change would result in a negative concentration
            if self.get_amount(solute, "mol").magnitude < 0:
                logger.warning(
                    "Attempted to set a negative concentration for solute %s. Concentration set to 0"
                    % solute
                )
                self.set_amount(solute, "0 mol")

            # calculate the volume occupied by all the solutes
            solute_vol = self._get_solute_volume()

            # determine the volume of solvent that will preserve the original volume
            target_vol = orig_volume - solute_vol

            # adjust the amount of solvent
            # volume in L, density in kg/m3 = g/L
            target_mass = (
                target_vol.magnitude * self.water_substance.rho * unit.Quantity("1 g")
            )

            mw = self.get_solvent().get_molecular_weight()
            target_mol = target_mass / mw
            self.get_solvent().moles = target_mol

        else:

            # change the amount of the solute present
            self.get_solute(solute).add_moles(
                amount, self.get_volume(), self.get_solvent_mass()
            )

            # set the amount to zero and log a warning if the desired amount
            # change would result in a negative concentration
            if self.get_amount(solute, "mol").magnitude < 0:
                logger.warning(
                    "Attempted to set a negative concentration for solute %s. Concentration set to 0"
                    % solute
                )
                self.set_amount(solute, "0 mol")

            # update the volume to account for the space occupied by all the solutes
            # make sure that there is still solvent present in the first place
            if self.get_solvent_mass() <= unit("0 kg"):
                logger.error("All solvent has been depleted from the solution")
                return None
            else:
                # set the volume recalculation flag
                self.volume_update_required = True

    def set_amount(self, solute, amount):
        """
        Set the amount of 'solute' in the parent solution.

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
        """
        # raise an error if a negative amount is specified
        if unit(amount).magnitude < 0:
            logger.error(
                "Negative amount specified for solute %s. Concentration not changed."
                % solute
            )

        # if units are given on a per-volume basis,
        # iteratively solve for the amount of solute that will preserve the
        # original volume and result in the desired concentration
        elif unit(amount).dimensionality in (
            "[substance]/[length]**3",
            "[mass]/[length]**3",
        ):

            # store the original volume for later
            orig_volume = self.get_volume()

            # change the amount of the solute present to match the desired amount
            self.get_solute(solute).set_moles(
                amount, self.get_volume(), self.get_solvent_mass()
            )

            # calculate the volume occupied by all the solutes
            solute_vol = self._get_solute_volume()

            # determine the volume of solvent that will preserve the original volume
            target_vol = orig_volume - solute_vol

            # adjust the amount of solvent
            target_mass = (
                target_vol.magnitude
                / 1000
                * self.water_substance.rho
                * unit.Quantity("1 kg")
            )
            mw = self.get_solvent().get_molecular_weight()
            target_mol = target_mass / mw
            self.get_solvent().moles = target_mol

        else:

            # change the amount of the solute present
            self.get_solute(solute).set_moles(
                amount, self.get_volume(), self.get_solvent_mass()
            )

            # update the volume to account for the space occupied by all the solutes
            # make sure that there is still solvent present in the first place
            if self.get_solvent_mass() <= unit("0 kg"):
                logger.error("All solvent has been depleted from the solution")
                return None
            else:
                self._update_volume()

    def get_osmolarity(self, activity_correction=False):
        """Return the osmolarity of the solution in Osm/L

        Parameters
        ----------
        activity_correction : bool
                If TRUE, the osmotic coefficient is used to calculate the
                osmolarity. This correction is appropriate when trying to predict
                the osmolarity that would be measured from e.g. freezing point
                depression. Defaults to FALSE if omitted.
        """
        if activity_correction is True:
            factor = self.get_osmotic_coefficient()
        else:
            factor = 1
        return factor * self.get_total_moles_solute() / self.get_volume().to("L")

    def get_osmolality(self, activity_correction=False):
        """Return the osmolality of the solution in Osm/kg

        Parameters
        ----------
        activity_correction : bool
                If TRUE, the osmotic coefficient is used to calculate the
                osmolarity. This correction is appropriate when trying to predict
                the osmolarity that would be measured from e.g. freezing point
                depression. Defaults to FALSE if omitted.
        """
        if activity_correction is True:
            factor = self.get_osmotic_coefficient()
        else:
            factor = 1
        return factor * self.get_total_moles_solute() / self.get_solvent_mass().to("kg")

    def get_total_moles_solute(self):
        """Return the total moles of all solute in the solution"""
        tot_mol = 0
        for item in self.components:
            if item != self.solvent_name:
                tot_mol += self.components[item].get_moles()
        return tot_mol

    def get_mole_fraction(self, solute):
        """
        Return the mole fraction of 'solute' in the solution

        Notes
        -----
        This function is DEPRECATED and will raise a warning when called.
        Use get_amount() instead and specify 'fraction' as the unit type.
        """
        logger.warning("get_mole_fraction is DEPRECATED! Use get_amount() instead.")
        return self.get_amount(solute, "fraction")

    def get_moles_solvent(self):
        """
        Return the moles of solvent present in the solution

        Parameters
        ----------
        None

        Returns
        -------
        Quantity
            The moles of solvent in the solution.

        """

        return self.get_amount(self.solvent_name, "mol")

    def get_salt(self):
        """
        Determine the predominant salt in a solution of ions.

        Many empirical equations for solution properties such as activity coefficient,
        partial molar volume, or viscosity are based on the concentration of
        single salts (e.g., NaCl). When multiple ions are present (e.g., a solution
        containing Na+, Cl-, and Mg+2), it is generally not possible to direclty model
        these quantities. pyEQL works around this problem by treating such solutions
        as single salt solutions.

        The get_salt() method examines the ionic composition of a solution and returns
        an object that identifies the single most predominant salt in the solution, defined
        by the cation and anion with the highest mole fraction. The Salt object contains
        information about the stoichiometry of the salt to enable its effective concentration
        to be calculated (e.g., 1 M MgCl2 yields 1 M Mg+2 and 2 M Cl-).

        Parameters
        ----------
        None

        Returns
        -------
        Salt
            Salt object containing information about the parent salt.

        See Also
        --------
        get_activity
        get_activity_coefficient
        get_water_activity
        get_osmotic_coefficient
        get_osmotic_pressure
        get_viscosity_kinematic

        Examples
        --------
        >>> s1 = Solution([['Na+','0.5 mol/kg'],['Cl-','0.5 mol/kg']])
        >>> s1.get_salt()
        <pyEQL.salt_ion_match.Salt object at 0x7fe6d3542048>
        >>> s1.get_salt().formula
        'NaCl'
        >>> s1.get_salt().nu_cation
        1
        >>> s1.get_salt().z_anion
        -1

        >>> s2 = pyEQL.Solution([['Na+','0.1 mol/kg'],['Mg+2','0.2 mol/kg'],['Cl-','0.5 mol/kg']])
        >>> s2.get_salt().formula
        'MgCl2'
        >>> s2.get_salt().nu_anion
        2
        >>> s2.get_salt().z_cation
        2
        """
        # identify the predominant salt in the solution
        return identify_salt(self)

    def get_salt_list(self):
        """
        Determine the predominant salt in a solution of ions.

        Many empirical equations for solution properties such as activity coefficient,
        partial molar volume, or viscosity are based on the concentration of
        single salts (e.g., NaCl). When multiple ions are present (e.g., a solution
        containing Na+, Cl-, and Mg+2), it is generally not possible to direclty model
        these quantities.

        The get_salt_list() method examines the ionic composition of a solution and
        simplifies it into a list of salts. The method retuns a dictionary of
        Salt objects where the keys are the salt formulas (e.g., 'NaCl'). The
        Salt object contains information about the stoichiometry of the salt to
        enable its effective concentration to be calculated
        (e.g., 1 M MgCl2 yields 1 M Mg+2 and 2 M Cl-).

        Parameters
        ----------
        None

        Returns
        -------
        dict
            A dictionary of Salt objects, keyed to the salt formula

        See Also
        --------
        get_activity
        get_activity_coefficient
        get_water_activity
        get_osmotic_coefficient
        get_osmotic_pressure
        get_viscosity_kinematic

        """
        # identify the predominant salt in the solution
        return generate_salt_list(self, unit="mol/kg")

    # Activity-related methods
    def get_activity_coefficient(
        self,
        solute: str,
        scale: Literal["molal", "molar", "fugacity", "rational"] = "molal",
        verbose: bool = False,
    ):
        """
        Return the activity coefficient of a solute in solution.

        The model used to calculte the activity coefficient is determined by the Solution's equation of state
        engine.

        Args:
            solute: The solute for which to retrieve the activity coefficient
            scale:  The activity coefficient concentration scale
            verbose: If True, pyEQL will print a message indicating the parent salt
                     that is being used for activity calculations. This option is
                     useful when modeling multicomponent solutions. False by default.

        Returns:
            Quantity: the activity coefficient as a dimensionless pint Quantity
        """
        # return unit activity coefficient if the concentration of the solute is zero
        if self.get_amount(solute, "mol").magnitude == 0:
            return unit("1 dimensionless")
        else:
            try:
                # get the molal-scale activity coefficient from the EOS engine
                molal = self.engine.get_activity_coefficient(
                    solution=self, solute=solute
                )
            except ValueError:
                logger.warning(
                    "Calculation unsuccessful. Returning unit activity coefficient."
                )
                return unit("1 dimensionless")

        # if necessary, convert the activity coefficient to another scale, and return the result
        if scale == "molal":
            return molal
        elif scale == "molar":
            total_molality = self.get_total_moles_solute() / self.get_solvent_mass()
            total_molarity = self.get_total_moles_solute() / self.get_volume()
            return (
                molal
                * self.water_substance.rho
                * unit.Quantity("1 g/L")
                * total_molality
                / total_molarity
            ).to("dimensionless")
        elif scale == "rational":
            return molal * (
                1
                + unit("0.018 kg/mol")
                * self.get_total_moles_solute()
                / self.get_solvent_mass()
            )
        else:
            logger.warning(
                "Invalid scale argument. Returning molal-scale activity coefficient"
            )
            return molal

    def get_activity(
        self,
        solute: str,
        scale: Literal["molal", "molar", "rational"] = "molal",
        verbose: bool = False,
    ):
        """
        Return the thermodynamic activity of the solute in solution on the chosen concentration scale.

        Args:
            solute : str
                    String representing the name of the solute of interest
            scale : str, optional
                    The concentration scale for the returned activity.
                    Valid options are "molal", "molar", and "rational" (i.e., mole fraction).
                    By default, the molal scale activity is returned.
            verbose : bool, optional
                    If True, pyEQL will print a message indicating the parent salt
                    that is being used for activity calculations. This option is
                    useful when modeling multicomponent solutions. False by default.

        Returns
        -------
        The thermodynamic activity of the solute in question (dimensionless)

        See Also
        --------
        get_activity_coefficient
        get_ionic_strength
        get_salt

        Notes
        -----
        The thermodynamic activity depends on the concentration scale used [#].
        By default, the ionic strength, activity coefficients, and activities are all
        calculated based on the molal (mol/kg) concentration scale.

        References
        ----------
        .. [#] Robinson, R. A.; Stokes, R. H. Electrolyte Solutions: Second Revised
               Edition; Butterworths: London, 1968, p.32.

        """
        # switch to the water activity function if the species is H2O
        if solute == "H2O" or solute == "water":
            activity = self.get_water_activity()
        else:
            # determine the concentration units to use based on the desired scale
            if scale == "molal":
                unit = "mol/kg"
            elif scale == "molar":
                unit = "mol/L"
            elif scale == "rational":
                unit = "fraction"
            else:
                logger.error("Invalid scale argument. Returning molal-scale activity.")
                unit = "mol/kg"
                scale = "molal"

            activity = (
                self.get_activity_coefficient(solute, scale=scale, verbose=verbose)
                * self.get_amount(solute, unit).magnitude
            )
            logger.info(
                "Calculated %s scale activity of solute %s as %s"
                % (scale, solute, activity)
            )

        return activity

    # TODO - engine method
    def get_osmotic_coefficient(
        self, scale: Literal["molal", "molar", "rational"] = "molal"
    ):
        """
        Return the osmotic coefficient of an aqueous solution.

        The method used depends on the Solution object's equation of state engine.

        """
        molal_phi = self.engine.get_osmotic_coefficient(self)

        if scale == "molal":
            return molal_phi
        elif scale == "rational":
            solvent = self.get_solvent().formula
            return (
                -molal_phi
                * unit("0.018 kg/mol")
                * self.get_total_moles_solute()
                / self.get_solvent_mass()
                / math.log(self.get_amount(solvent, "fraction"))
            )
        elif scale == "fugacity":
            solvent = self.get_solvent().formula
            return math.exp(
                -molal_phi
                * unit("0.018 kg/mol")
                * self.get_total_moles_solute()
                / self.get_solvent_mass()
                - math.log(self.get_amount(solvent, "fraction"))
            )
        else:
            logger.warning(
                "Invalid scale argument. Returning molal-scale osmotic coefficient"
            )
            return molal_phi

    def get_water_activity(self):
        """
        Return the water activity.

        Returns
        -------
        Quantity :
            The thermodynamic activity of water in the solution.

        See Also
        --------
        get_osmotic_coefficient
        get_ionic_strength
        get_salt

        Notes
        -----
        Water activity is related to the osmotic coefficient in a solution containing i solutes by: [#]_

        .. math:: \\ln a_w = - \\Phi M_w \\sum_i m_i

        Where :math:`M_w` is the molar mass of water (0.018015 kg/mol) and :math:`m_i` is the molal concentration
        of each species.

        If appropriate Pitzer model parameters are not available, the
        water activity is assumed equal to the mole fraction of water.

        References
        ----------
        .. [#] Blandamer, Mike J., Engberts, Jan B. F. N., Gleeson, Peter T., Reis, Joao Carlos R., 2005. "Activity of
        water in aqueous systems: A frequently neglected property." *Chemical Society Review* 34, 440-458.

        Examples
        --------
        >>> s1 = pyEQL.Solution([['Na+','0.3 mol/kg'],['Cl-','0.3 mol/kg']])
        >>> s1.get_water_activity()
        <Quantity(0.9900944932888518, 'dimensionless')>
        """
        """
        pseudo code

        identify predominant salt for coefficients
        check if coefficients exist for that salt
        if so => calc osmotic coefficient and log an info message

        if not = > return mole fraction and log a warning message

        """
        osmotic_coefficient = self.get_osmotic_coefficient()

        if osmotic_coefficient == 1:
            logger.warning(
                "Pitzer parameters not found. Water activity set equal to mole fraction"
            )
            return self.get_amount("H2O", "fraction")
        else:
            concentration_sum = unit("0 mol/kg")
            for item in self.components:
                if item == "H2O":
                    pass
                else:
                    concentration_sum += self.get_amount(item, "mol/kg")

            logger.info("Calculated water activity using osmotic coefficient")

            return math.exp(
                -osmotic_coefficient * 0.018015 * unit("kg/mol") * concentration_sum
            ) * unit("1 dimensionless")

    def get_ionic_strength(self):
        """
        Return the ionic strength of the solution.

        Return the ionic strength of the solution, calculated as 1/2 * sum ( molality * charge ^2) over all the ions.
        Molal (mol/kg) scale concentrations are used for compatibility with the activity correction formulas.

        Returns
        -------
        Quantity :
            The ionic strength of the parent solution, mol/kg.

        See Also
        --------
        get_activity
        get_water_activity

        Notes
        -----
        The ionic strength is calculated according to:

        .. math:: I = \\sum_i m_i z_i^2

        Where :math:`m_i` is the molal concentration and :math:`z_i` is the charge on species i.

        Examples
        --------
        >>> s1 = pyEQL.Solution([['Na+','0.2 mol/kg'],['Cl-','0.2 mol/kg']])
        >>> s1.get_ionic_strength()
        <Quantity(0.20000010029672785, 'mole / kilogram')>

        >>> s1 = pyEQL.Solution([['Mg+2','0.3 mol/kg'],['Na+','0.1 mol/kg'],['Cl-','0.7 mol/kg']],temperature='30 degC')
        >>> s1.get_ionic_strength()
        <Quantity(1.0000001004383303, 'mole / kilogram')>
        """
        self.ionic_strength = 0
        for solute in self.components.keys():
            self.ionic_strength += (
                0.5
                * self.get_amount(solute, "mol/kg")
                * self.components[solute].get_formal_charge() ** 2
            )

        return self.ionic_strength

    def get_charge_balance(self):
        """
        Return the charge balance of the solution.

        Return the charge balance of the solution. The charge balance represents the net electric charge
        on the solution and SHOULD equal zero at all times, but due to numerical errors will usually
        have a small nonzero value.

        Returns
        -------
        float :
            The charge balance of the solution, in equivalents.

        Notes
        -----
        The charge balance is calculated according to:

        .. math:: CB = F \\sum_i n_i z_i

        Where :math:`n_i` is the number of moles, :math:`z_i` is the charge on species i, and :math:`F` is the Faraday constant.

        """
        self.charge_balance = 0
        for solute in self.components.keys():
            self.charge_balance += (
                self.get_amount(solute, "mol")
                * self.components[solute].get_formal_charge()
                * unit.e
                * unit.N_A
            )

        return self.charge_balance.magnitude

    def get_alkalinity(self):
        """
        Return the alkalinity or acid neutralizing capacity of a solution

        Returns
        -------
        Quantity :
            The alkalinity of the solution in mg/L as CaCO3

        Notes
        -----
        The alkalinity is calculated according to: [#]_

        .. math:: Alk = F \\sum_i z_i C_B - \\sum_i z_i C_A

        Where :math:`C_B` and :math:`C_A` are conservative cations and anions, respectively
        (i.e. ions that do not participate in acid-base reactions), and :math:`z_i` is their charge.
        In this method, the set of conservative cations is all Group I and Group II cations, and the conservative anions
        are all the anions of strong acids.

        References
        ----------
        .. [#] Stumm, Werner and Morgan, James J. Aquatic Chemistry, 3rd ed,
               pp 165. Wiley Interscience, 1996.
        """
        alkalinity = 0 * unit("mol/L")
        equiv_wt_CaCO3 = 100.09 / 2 * unit("g/mol")

        base_cations = [
            "Li+",
            "Na+",
            "K+",
            "Rb+",
            "Cs+",
            "Fr+",
            "Be+2",
            "Mg+2",
            "Ca+2",
            "Sr+2",
            "Ba+2",
            "Ra+2",
        ]
        acid_anions = ["Cl-", "Br-", "I-", "SO4-2", "NO3-", "ClO4-", "ClO3-"]

        for item in self.components:
            if item in base_cations:
                z = self.get_solute(item).get_formal_charge()
                alkalinity += self.get_amount(item, "mol/L") * z
            if item in acid_anions:
                z = self.get_solute(item).get_formal_charge()
                alkalinity -= self.get_amount(item, "mol/L") * z

        # convert the alkalinity to mg/L as CaCO3
        return (alkalinity * equiv_wt_CaCO3).to("mg/L")

    def get_hardness(self):
        """
        Return the hardness of a solution.

        Hardness is defined as the sum of the equivalent concentrations
        of multivalent cations as calcium carbonate.

        NOTE: at present pyEQL cannot distinguish between mg/L as CaCO3
        and mg/L units. Use with caution.

        Parameters
        ----------
        None

        Returns
        -------
        Quantity
            The hardness of the solution in mg/L as CaCO3

        """
        hardness = 0 * unit("mol/L")
        equiv_wt_CaCO3 = 100.09 / 2 * unit("g/mol")

        for item in self.components:
            z = self.get_solute(item).get_formal_charge()
            if z > 1:
                hardness += z * self.get_amount(item, "mol/L")

        # convert the hardness to mg/L as CaCO3
        return (hardness * equiv_wt_CaCO3).to("mg/L")

    def get_debye_length(self):
        """
        Return the Debye length of a solution

        Debye length is calculated as [#]_

        .. math::

            \\kappa^{-1} = \\sqrt({\\epsilon_r \\epsilon_o k_B T \\over (2 N_A e^2 I)})

        where :math:`I` is the ionic strength, :math:`epsilon_r` and :math:`epsilon_r`
        are the relative permittivity and vacuum permittivity, :math:`k_B` is the
        Boltzmann constant, and :math:`T` is the temperature, :math:`e` is the
        elementary charge, and :math:`N_A` is Avogadro's number.

        Parameters
        ----------
        None

        Returns
        -------
        Quantity
            The Debye length, in nanometers.

        References
        ----------
        .. [#] https://en.wikipedia.org/wiki/Debye_length#Debye_length_in_an_electrolyte

        See Also
        --------
        get_ionic_strength()
        get_dielectric_constant()

        """
        # to preserve dimensionality, convert the ionic strength into mol/L units
        ionic_strength = self.get_ionic_strength().magnitude * unit("mol/L")
        dielectric_constant = self.get_dielectric_constant()

        debye_length = (
            dielectric_constant
            * unit.epsilon_0
            * unit.k
            * self.temperature
            / (2 * unit.N_A * unit.e**2 * ionic_strength)
        ) ** 0.5

        return debye_length.to("nm")

    def get_bjerrum_length(self):
        """
        Return the Bjerrum length of a solution

        Bjerrum length representes the distance at which electrostatic
        interactions between particles become comparable in magnitude
        to the thermal energy.:math:`\\lambda_B` is calculated as [#]_

        .. math::

            \\lambda_B = {e^2 \\over (4 \\pi \\epsilon_r \\epsilon_o k_B T)}

        where :math:`e` is the fundamental charge, :math:`epsilon_r` and :math:`epsilon_r`
        are the relative permittivity and vacuum permittivity, :math:`k_B` is the
        Boltzmann constant, and :math:`T` is the temperature.

        Parameters
        ----------
        None

        Returns
        -------
        Quantity
            The Bjerrum length, in nanometers.

        References
        ----------
        .. [#] https://en.wikipedia.org/wiki/Bjerrum_length

        Examples
        --------
        >>> s1 = pyEQL.Solution()
        >>> s1.get_bjerrum_length()
        <Quantity(0.7152793009386953, 'nanometer')>

        See Also
        --------
        get_dielectric_constant()

        """
        dielectric_constant = self.get_dielectric_constant()

        bjerrum_length = unit.e**2 / (
            4
            * math.pi
            * dielectric_constant
            * unit.epsilon_0
            * unit.k
            * self.temperature
        )
        return bjerrum_length.to("nm")

    def get_transport_number(self, solute, activity_correction=False):
        """Calculate the transport number of the solute in the solution

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
        Transport number is calculated according to [#]_ :

        .. math::

            t_i = {D_i z_i^2 C_i \\over \\sum D_i z_i^2 C_i}

        Where :math:`C_i` is the concentration in mol/L, :math:`D_i` is the diffusion
        coefficient, and :math:`z_i` is the charge, and the summation extends
        over all species in the solution.

        If `activity_correction` is True, the contribution of each ion to the
        transport number is corrected with an activity factor. See the documentation
        for get_conductivity() for an explanation of this correction.

        References
        ----------
        .. [#] Geise, G. M.; Cassady, H. J.; Paul, D. R.; Logan, E.; Hickner, M. A. "Specific
        ion effects on membrane potential and the permselectivity of ion exchange membranes.""
        *Phys. Chem. Chem. Phys.* 2014, 16, 21673–21681.

        """
        denominator = 0
        numerator = 0

        for item in self.components:

            z = self.get_solute(item).get_formal_charge()
            term = (
                self.get_property(item, "diffusion_coefficient")
                * z**2
                * self.get_amount(item, "mol/L")
            )

            if activity_correction is True:
                gamma = self.get_activity_coefficient(item)

                if self.get_ionic_strength().magnitude < 0.36 * z:
                    alpha = 0.6 / z**0.5
                else:
                    alpha = self.get_ionic_strength().magnitude ** 0.5 / z

                if item == solute:
                    numerator = term * gamma**alpha

                denominator += term * gamma**alpha

            else:
                if item == solute:
                    numerator = term

                denominator += term

        return (numerator / denominator).to("dimensionless")

    def get_molar_conductivity(self, solute):

        """
        Calculate the molar (equivalent) conductivity for a solute

        Parameters
        ----------
        solute : str
            String identifying the solute for which the molar conductivity is
            to be calculated.

        Returns
        -------
        float
                The molar or equivalent conductivity of the species in the solution.
                Zero if the solute is not charged.

        Notes
        -----
        Molar conductivity is calculated from the Nernst-Einstein relation [#]_

        .. math::

            \\kappa_i = {z_i^2 D_i F^2 \\over RT}

        Note that the diffusion coefficient is strongly variable with temperature.

        References
        ----------

        .. [#] Smedley, Stuart. The Interpretation of Ionic Conductivity in Liquids, pp 1-9. Plenum Press, 1980.

        Examples
        --------
        TODO

        """
        D = self.get_property(solute, "diffusion_coefficient")

        molar_cond = (
            D
            * (unit.e * unit.N_A) ** 2
            * self.get_solute(solute).get_formal_charge() ** 2
            / (unit.R * self.temperature)
        )

        logger.info(
            "Computed molar conductivity as %s from D = %s at T=%s"
            % (molar_cond, str(D), self.temperature)
        )

        return molar_cond.to("mS / cm / (mol/L)")

    def get_mobility(self, solute):
        """
        Calculate the ionic mobility of the solute

        Parameters
        ----------
        solute : str
            String identifying the solute for which the mobility is
            to be calculated.

        Returns
        -------
        float : The ionic mobility. Zero if the solute is not charged.


        Notes
        -----
        This function uses the Einstein relation to convert a diffusion coefficient
        into an ionic mobility [#]_

        .. math::

            \\mu_i = {F |z_i| D_i \\over RT}

        References
        ----------
        .. [#] Smedley, Stuart I. The Interpretation of Ionic Conductivity in Liquids. Plenum Press, 1980.

        """
        D = self.get_property(solute, "diffusion_coefficient")

        mobility = (
            unit.N_A
            * unit.e
            * abs(self.get_solute(solute).get_formal_charge())
            * D
            / (unit.R * self.temperature)
        )

        logger.info(
            "Computed ionic mobility as %s from D = %s at T=%s"
            % (mobility, str(D), self.temperature)
        )

        return mobility.to("m**2/V/s")

    def get_property(self, solute, name):
        """Retrieve a thermodynamic property (such as diffusion coefficient)
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

        """
        # retrieve the base value and the conditions of measurement from the
        # database

        if db.has_parameter(solute, name):
            base_value = self.get_solute(solute).get_parameter(name)
        else:
            base_value = None

        base_temperature = unit("25 degC")
        # base_pressure = unit("1 atm")

        # perform temperature-corrections or other adjustments for certain
        # parameter types
        if name == "diffusion_coefficient":
            if base_value is not None:
                # correct for temperature and viscosity
                # .. math:: D_1 \over D_2 = T_1 \over T_2 * \mu_2 \over \mu_1
                # where :math:`\mu` is the dynamic viscosity
                # assume that the base viscosity is that of pure water
                return (
                    base_value
                    * self.temperature
                    / base_temperature
                    * self.water_substance.mu
                    * unit.Quantity("1 Pa*s")
                    / self.get_viscosity_dynamic()
                )
            else:
                logger.warning(
                    "Diffusion coefficient not found for species %s. Assuming zero."
                    % (solute)
                )
                return unit("0 m**2/s")

        # just return the base-value molar volume for now; find a way to adjust for
        # concentration later
        if name == "partial_molar_volume":
            # calculate the partial molar volume for water since it isn't in the database
            if solute == "H2O":
                vol = (
                    self.get_solute("H2O").get_molecular_weight()
                    / self.water_substance.rho
                    * unit.Quantity("1 g/L")
                )

                return vol.to("cm **3 / mol")
            else:
                if base_value is not None:
                    return base_value
                    if self.temperature != base_temperature:
                        logger.warning(
                            "Partial molar volume for species %s not corrected for temperature"
                            % solute
                        )
                else:
                    logger.warning(
                        "Partial molar volume not found for species %s. Assuming zero."
                        % solute
                    )
                    return unit("0 cm **3 / mol")

        # for parameters not named above, just return the base value
        else:
            logger.warning("%s has not been corrected for solution conditions" % name)
            return base_value

    def get_chemical_potential_energy(self, activity_correction=True):
        """
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
        calculated as follows: [#]_

        .. math::      E = R T \\sum_i n_i  \\ln a_i

        or

        .. math::      E = R T \\sum_i n_i \\ln x_i

        Where :math:`n` is the number of moles of substance, :math:`T` is the temperature in kelvin,
        :math:`R` the ideal gas constant, :math:`x` the mole fraction, and :math:`a` the activity of
        each component.

        Note that dissociated ions must be counted as separate components,
        so a simple salt dissolved in water is a three component solution (cation,
        anion, and water).

        References
        ----------
        .. [#] Koga, Yoshikata, 2007. *Solution Thermodynamics and its Application to Aqueous Solutions: A differential approach.* Elsevier, 2007, pp. 23-37.

        Examples
        --------

        """
        E = unit("0 J")

        # loop through all the components and add their potential energy
        for item in self.components:
            try:
                if activity_correction is True:
                    E += (
                        unit.R
                        * self.temperature.to("K")
                        * self.get_amount(item, "mol")
                        * math.log(self.get_activity(item))
                    )
                else:
                    E += (
                        unit.R
                        * self.temperature.to("K")
                        * self.get_amount(item, "mol")
                        * math.log(self.get_amount(item, "fraction"))
                    )
            # If we have a solute with zero concentration, we will get a ValueError
            except ValueError:
                continue

        return E.to("J")

    def get_lattice_distance(self, solute):
        """
        Calculate the average distance between molecules

        Calculate the average distance between molecules of the given solute,
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

        .. math:: d = ( C_i N_A ) ^ {-{1 \\over 3}}

        """
        # calculate the volume per particle as the reciprocal of the molar concentration
        # (times avogadro's number). Take the cube root of the volume to get
        # the average distance between molecules
        distance = (self.get_amount(solute, "mol/L") * unit.N_A) ** (-1 / 3)

        return distance.to("nm")

    def _update_volume(self):
        """
        Recalculate the solution volume based on composition

        """
        self.volume = self._get_solvent_volume() + self._get_solute_volume()

    def _get_solvent_volume(self):
        """
        Return the volume of the pure solvent

        """
        # calculate the volume of the pure solvent
        solvent_vol = self.get_solvent_mass() / (
            self.water_substance.rho * unit.Quantity("1 g/L")
        )

        return solvent_vol.to("L")

    def _get_solute_volume(self):
        """
        Return the volume of only the solutes

        """
        return self.engine.get_solute_volume(self)

    def copy(self):
        """Return a copy of the solution

        TODO - clarify whether this is a deep or shallow copy
        """
        # prepare to copy the bulk properties
        new_temperature = str(self.temperature)
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
                new_solutes.append([item, str(self.get_amount(item, "mol"))])

        # create the new solution
        return Solution(
            new_solutes,
            solvent=[new_solvent, new_solvent_mass],
            temperature=new_temperature,
            pressure=new_pressure,
        )

    # informational methods
    def list_solutes(self):
        """
        List all the solutes in the solution.

        """
        return list(self.components.keys())

    def list_concentrations(self, unit="mol/kg", decimals=4, type="all"):
        """
        List the concentration of each species in a solution.

        Parameters
        ----------
        unit: str
            String representing the desired concentration unit.
        decimals: int
            The number of decimal places to display. Defaults to 4.
        type     : str
            The type of component to be sorted. Defaults to 'all' for all
            solutes. Other valid arguments are 'cations' and 'anions' which
            return lists of cations and anions, respectively.

        Returns
        -------
        dict
            Dictionary containing a list of the species in solution paired with their amount in the specified units

        """
        result_list = []
        # populate a list with component names

        if type == "all":
            print("Component Concentrations:\n")
            print("========================\n")
            for item in self.components:
                amount = self.get_amount(item, unit)
                result_list.append([item, amount])
                print(
                    item
                    + ":"
                    + "\t {0:0.{decimals}f~}".format(amount, decimals=decimals)
                )
        elif type == "cations":
            print("Cation Concentrations:\n")
            print("========================\n")
            for item in self.components:
                if self.components[item].get_formal_charge() > 0:
                    amount = self.get_amount(item, unit)
                    result_list.append([item, amount])
                    print(
                        item
                        + ":"
                        + "\t {0:0.{decimals}f~}".format(amount, decimals=decimals)
                    )
        elif type == "anions":
            print("Anion Concentrations:\n")
            print("========================\n")
            for item in self.components:
                if self.components[item].get_formal_charge() < 0:
                    amount = self.get_amount(item, unit)
                    result_list.append([item, amount])
                    print(
                        item
                        + ":"
                        + "\t {0:0.{decimals}f~}".format(amount, decimals=decimals)
                    )

        return result_list

    def list_salts(self, unit="mol/kg", decimals=4):
        list = generate_salt_list(self, unit)
        for item in list:
            print(
                item.formula
                + "\t {:0.{decimals}f}".format(list[item], decimals=decimals)
            )

    def list_activities(self, decimals=4):
        """
        List the activity of each species in a solution.

        Parameters
        ----------
        decimals: int
            The number of decimal places to display. Defaults to 4.

        Returns
        -------
        dict
            Dictionary containing a list of the species in solution paired with their activity

        """
        print("Component Activities:\n")
        print("=====================\n")
        for i in self.components.keys():
            print(
                i
                + ":"
                + "\t {0.magnitude:0.{decimals}f}".format(
                    self.get_activity(i), decimals=decimals
                )
            )

    def __str__(self):
        # set output of the print() statement for the solution
        str1 = f"Volume: {self.get_volume():.3f~}\n"
        str2 = f"Pressure: {self.pressure:.3f~}\n"
        str3 = f"Temperature: {self.temperature:.3f~}\n"
        str4 = f"Components: {self.list_solutes():}\n"
        return str1 + str2 + str3 + str4
