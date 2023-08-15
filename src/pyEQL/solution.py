"""
pyEQL Solution Class.

:copyright: 2013-2023 by Ryan S. Kingsbury
:license: LGPL, see LICENSE for more details.

"""

import math
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Literal, Optional, Union

from iapws import IAPWS95
from maggma.stores import JSONStore, Store
from monty.dev import deprecated
from monty.json import MSONable
from pint import DimensionalityError, Quantity
from pymatgen.core.ion import Ion

from pyEQL import unit
from pyEQL.engines import EOS, IdealEOS, NativeEOS

# logging system
from pyEQL.logging_system import logger
from pyEQL.salt_ion_match import generate_salt_list, identify_salt


class Solution(MSONable):
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
        pE: float = 8.5,
        solvent: Union[str, list] = "H2O",
        engine: Literal["native", "ideal"] = "native",
        database: Optional[Union[str, Path, Store]] = None,
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
            pe: the pE value (redox potential) of the solution.     Lower values = more reducing,
                higher values = more oxidizing. At pH 7, water is stable between approximately
                -7 to +14. The default value corresponds to a pE value typical of natural
                waters in equilibrium with the atmosphere.
            solvent: Formula of the solvent. Solvents other than water are not supported at
                this time.
            engine:
            database: path to a .json file (str or Path) or maggma Store instance that
                contains serialized SoluteDocs. `None` (default) will use the built-in pyEQL database.

        Examples:
            >>> s1 = pyEQL.Solution([['Na+','1 mol/L'],['Cl-','1 mol/L']],temperature='20 degC',volume='500 mL')
            >>> print(s1)
            Components:
            ['H2O', 'Cl-', 'H+', 'OH-', 'Na+']
            Volume: 0.5 l
            Density: 1.0383030844030992 kg/l
        """
        # create a logger attached to this class
        # self.logger = logging.getLogger(type(self).__name__)

        # per-instance cache of get_property calls
        self.get_property = lru_cache(maxsize=None)(self._get_property)

        # initialize the volume with a flag to distinguish user-specified volume
        if volume is not None:
            # volume_set = True
            self.volume = unit.Quantity(volume).to("L")
        else:
            # volume_set = False
            self.volume = unit.Quantity("1 L")
        # store the initial conditions as private variables in case they are
        # changed later
        self._temperature = unit.Quantity(temperature)
        self._pressure = unit.Quantity(pressure)
        self._pE = pE
        self._pH = pH
        self.pE = self._pE

        # instantiate a water substance for property retrieval
        self.water_substance = IAPWS95(
            T=self.temperature.magnitude,
            P=self.pressure.to("MPa").magnitude,
        )

        # create an empty dictionary of components. This dict comprises {formula: moles}
        #  where moles is the number of moles in the solution.
        self.components: dict = {}

        # initialize the volume recalculation flag
        self.volume_update_required = False

        # connect to the desired property database
        if not isinstance(database, Store):
            if database is None:
                from pkg_resources import resource_filename

                database_dir = resource_filename("pyEQL", "database")
                json = Path(database_dir) / "pyeql_db.json"
            else:
                json = database if isinstance(database, str) else str(database)
            db_store = JSONStore(json, key="formula")
            logger.info(f"Created maggma JSONStore from .json file {database}")
        else:
            db_store = database
        self.database = db_store
        self.database.connect()
        logger.info(f"Connected to property database {self.database!s}")

        # set the equation of state engine
        self._engine = engine
        # self.engine: Optional[EOS] = None
        if self._engine == "ideal":
            self.engine: EOS = IdealEOS()
        elif self._engine == "native":
            self.engine = NativeEOS()
        else:
            raise ValueError(f'{engine} is not a valid value for the "engine" kwarg!')

        # define the solvent. Allow for list input to support future use of mixed solvents
        if not isinstance(solvent, list):
            solvent = [solvent]
        if len(solvent) > 1:
            raise ValueError("Multiple solvents are not yet supported!")
        if solvent[0] not in ["H2O", "H2O(aq)", "water", "Water", "HOH"]:
            raise ValueError("Non-aqueous solvent detected. These are not yet supported!")
        self.solvent = solvent[0]

        # TODO - do I need the ability to specify the solvent mass?
        # # raise an error if the solvent volume has also been given
        # if volume_set is True:
        #     logger.error(
        #         "Solvent volume and mass cannot both be specified. Calculating volume based on solvent mass."
        #     )
        # # add the solvent and the mass
        # self.add_solvent(self.solvent, kwargs["solvent"][1])

        # calculate the moles of solvent (water) on the density and the solution volume
        moles = self.volume / unit.Quantity("55.55 mol/L")
        self.components["H2O"] = moles.magnitude

        # set the pH with H+ and OH-
        self.add_solute("H+", str(10 ** (-1 * pH)) + "mol/L")
        self.add_solute("OH-", str(10 ** (-1 * (14 - pH))) + "mol/L")

        # populate the other solutes
        self._solutes = solutes
        if self._solutes is None:
            self._solutes = {}
        if isinstance(self._solutes, dict):
            for k, v in self._solutes.items():
                self.add_solute(k, v)
        elif isinstance(self._solutes, list):
            logger.warning(
                'List input of solutes (e.g., [["Na+", "0.5 mol/L]]) is deprecated! Use dictionary formatted input (e.g., {"Na+":"0.5 mol/L"} instead.)'
            )
            for item in self._solutes:
                self.add_solute(*item)
        elif self._solutes is not None:
            raise ValueError("Solutes must be given as a list or dict!")

    def add_solute(self, formula, amount):
        """Primary method for adding substances to a pyEQL solution.

        Parameters
        ----------
        formula : str
                    Chemical formula for the solute.
                    Charged species must contain a + or - and (for polyvalent solutes) a number representing the net charge (e.g. 'SO4-2').
        amount : str
                    The amount of substance in the specified unit system. The string should contain both a quantity and
                    a pint-compatible representation of a unit. e.g. '5 mol/kg' or '0.1 g/L'
        """
        # if units are given on a per-volume basis,
        # iteratively solve for the amount of solute that will preserve the
        # original volume and result in the desired concentration
        if unit.Quantity(amount).dimensionality in (
            "[substance]/[length]**3",
            "[mass]/[length]**3",
        ):
            # store the original volume for later
            orig_volume = self.get_volume()

            # add the new solute
            quantity = unit.Quantity(amount)
            mw = self.get_property(formula, "molecular_weight")  # returns a quantity
            target_mol = quantity.to("moles", "chem", mw=mw, volume=self.volume, solvent_mass=self.get_solvent_mass())
            self.components[formula] = target_mol.to("moles").magnitude

            # calculate the volume occupied by all the solutes
            solute_vol = self._get_solute_volume()

            # determine the volume of solvent that will preserve the original volume
            target_vol = orig_volume - solute_vol

            # adjust the amount of solvent
            # density is returned in kg/m3 = g/L
            target_mass = target_vol.to("L").magnitude * self.water_substance.rho * unit.Quantity("1 g")
            # mw = unit.Quantity(self.get_property(self.solvent_name, "molecular_weight"))
            mw = self.get_property(self.solvent, "molecular_weight")
            if mw is None:
                raise ValueError(
                    f"Molecular weight for solvent {self.solvent} not found in database. This is required to proceed."
                )
            target_mol = target_mass.to("g") / mw.to("g/mol")
            self.components[self.solvent] = target_mol.magnitude

        else:
            # add the new solute
            quantity = unit.Quantity(amount)
            mw = unit.Quantity(self.get_property(formula, "molecular_weight"))
            target_mol = quantity.to("moles", "chem", mw=mw, volume=self.volume, solvent_mass=self.get_solvent_mass())
            self.components[formula] = target_mol.to("moles").magnitude

            # update the volume to account for the space occupied by all the solutes
            # make sure that there is still solvent present in the first place
            if self.get_solvent_mass() <= unit.Quantity("0 kg"):
                logger.error("All solvent has been depleted from the solution")
                return
            # set the volume recalculation flag
            self.volume_update_required = True

    # TODO - deprecate this method. Solvent should be added to the dict like anything else
    # and solvent_name will track which component it is.
    def add_solvent(self, formula, amount):
        """Same as add_solute but omits the need to pass solvent mass to pint."""
        quantity = unit.Quantity(amount)
        mw = unit.Quantity(self.get_property(formula, "molecular_weight"))
        target_mol = quantity.to("moles", "chem", mw=mw, volume=self.volume, solvent_mass=self.get_solvent_mass())
        self.components[formula] = target_mol.to("moles").magnitude

    @property
    def temperature(self) -> Quantity:
        """Return the temperature of the solution in Kelvin."""
        return self._temperature.to("K")

    @temperature.setter
    def temperature(self, temperature: str):
        """
        Set the solution temperature.

        Args:
            temperature: pint-compatible string, e.g. '25 degC'
        """
        self._temperature = unit.Quantity(temperature)
        # recalculate the volume
        self.volume_update_required = True

    @property
    def pH(self) -> Quantity:
        """Return the pH of the solution."""
        return self.p("H+", activity=True)

    @property
    def pressure(self) -> Quantity:
        """Return the hydrostatic pressure of the solution in atm."""
        return self._pressure.to("atm")

    @pressure.setter
    def pressure(self, pressure: str):
        """
        Set the solution pressure.

        Args:
            pressure: pint-compatible string, e.g. '1.2 atmC'
        """
        self._pressure = unit.Quantity(pressure)
        # recalculate the volume
        self.volume_update_required = True

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
        :py:meth:`get_amount()`
        """
        # return the total mass (kg) of the solvent
        mw = self.get_property(self.solvent, "molecular_weight").to("kg/mol").magnitude

        return self.components[self.solvent] * mw * unit.Quantity("1 kg")

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
        all component concentrations.

        Parameters
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

        """
        # figure out the factor to multiply the old concentrations by
        scale_factor = unit.Quantity(volume) / self.get_volume()

        # scale down the amount of all the solutes according to the factor
        for solute in self.components:
            self.components[solute] *= scale_factor.magnitude

        # update the solution volume
        self.volume = unit.Quantity(volume)

    @property
    def mass(self) -> Quantity:
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

    @property
    def density(self) -> Quantity:
        """
        Return the density of the solution.

        Density is calculated from the mass and volume each time this method is called.

        Returns
        -------
        Quantity: The density of the solution.
        """
        return self.mass / self.get_volume()

    @property
    def dielectric_constant(self) -> Quantity:
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
        Implements the following equation as given by Zuber et al.

        .. math:: \\epsilon = \\epsilon_{solvent} \\over 1 + \\sum_i \\alpha_i x_i

        where :math:`\\alpha_i` is a coefficient specific to the solvent and ion, and :math:`x_i`
        is the mole fraction of the ion in solution.


        References
        ----------
        .A. Zuber, L. Cardozo-Filho, V.F. Cabral, R.F. Checoni, M. Castier,
        An empirical equation for the dielectric constant in aqueous and nonaqueous
        electrolyte mixtures, Fluid Phase Equilib. 376 (2014) 116-123.
        doi:10.1016/j.fluid.2014.05.037.
        """
        di_water = self.water_substance.epsilon

        denominator = 1
        for item in self.components:
            # ignore water
            if item != "H2O":
                # skip over solutes that don't have parameters
                # try:
                fraction = self.get_amount(item, "fraction")
                coefficient = self.get_property(item, "model_parameters.dielectric_zuber")
                if coefficient is not None:
                    denominator += coefficient * fraction
                # except TypeError:
                #     logger.warning("No dielectric parameters found for species %s." % item)
                # continue

        return unit.Quantity(di_water / denominator, "dimensionless")

    # TODO - need tests for viscosity
    @property
    def viscosity_dynamic(self) -> Quantity:
        """
        Return the dynamic (absolute) viscosity of the solution.

        Calculated from the kinematic viscosity

        See Also:
        --------
        viscosity_kinematic
        """
        return self.viscosity_kinematic * self.density

    # TODO - before deprecating get_viscosity_relative, consider whether the Jones-Dole
    # model should be integrated here as a fallback, in case salt parameters for the
    # other model are not available.
    # if self.ionic_strength.magnitude > 0.2:
    #   logger.warning('Viscosity calculation has limited accuracy above 0.2m')

    #        viscosity_rel = 1
    #        for item in self.components:
    #            # ignore water
    #            if item != 'H2O':
    #                # skip over solutes that don't have parameters
    #                try:
    #                    conc = self.get_amount(item,'mol/kg').magnitude
    #                    coefficients= self.get_property(item, 'jones_dole_viscosity')
    #                    viscosity_rel += coefficients[0] * conc ** 0.5 + coefficients[1] * conc + \
    #                    coefficients[2] * conc ** 2
    #                except TypeError:
    #                    continue
    # return (
    #     self.viscosity_dynamic / self.water_substance.mu * unit.Quantity("1 Pa*s")
    # )
    @property
    def viscosity_kinematic(self):
        """
        Return the kinematic viscosity of the solution.

        Notes
        -----
        The calculation is based on a model derived from the Eyring equation
        and presented in

        .. math::

            \\ln \\nu = \\ln {\\nu_w MW_w \\over \\sum_i x_i MW_i } +
            15 x_+^2 + x_+^3  \\delta G^*_{123} + 3 x_+ \\delta G^*_{23} (1-0.05x_+)

        Where:

        .. math:: \\delta G^*_{123} = a_o + a_1 (T)^{0.75}
        .. math:: \\delta G^*_{23} = b_o + b_1 (T)^{0.5}

        In which :math:`\\nu` is the kinematic viscosity, MW is the molecular weight,
        :math:`x_{+}` is the mole fraction of cations, and :math:`T` is the temperature in degrees C.

        The a and b fitting parameters for a variety of common salts are included in the
        database.

        References
        ----------
        Vásquez-Castillo, G.; Iglesias-Silva, G. a.; Hall, K. R. An extension of the McAllister model to correlate kinematic viscosity of electrolyte solutions. Fluid Phase Equilib. 2013, 358, 44-49.

        See Also:
        --------
        :py:meth:`viscosity_dynamic`
        """
        # identify the main salt in the solution
        salt = self.get_salt()
        # reverse-convert the sanitized formula back to whatever was in self.components
        for i in self.components:
            if Ion.from_formula(i).reduced_formula == salt.cation:
                cation = i

        a0 = a1 = b0 = b1 = 0

        # retrieve the parameters for the delta G equations
        params = self.get_property(salt.formula, "model_parameters.viscosity_eyring")
        if params is not None:
            a0 = unit.Quantity(params["a0"]["value"]).magnitude
            a1 = unit.Quantity(params["a1"]["value"]).magnitude
            b0 = unit.Quantity(params["b0"]["value"]).magnitude
            b1 = unit.Quantity(params["b1"]["value"]).magnitude
        else:
            # proceed with the coefficients equal to zero and log a warning
            logger.warning("Viscosity coefficients for %s not found. Viscosity will be approximate." % salt.formula)

        # compute the delta G parameters
        temperature = self.temperature.to("degC").magnitude
        G_123 = a0 + a1 * (temperature) ** 0.75
        G_23 = b0 + b1 * (temperature) ** 0.5

        # get the kinematic viscosity of water, returned by IAPWS in m2/s
        nu_w = self.water_substance.nu

        # compute the effective molar mass of the solution
        MW = self.mass / (self.get_moles_solvent() + self.get_total_moles_solute())

        # get the MW of water
        MW_w = unit.Quantity(self.get_property(self.solvent, "molecular_weight"))

        # calculate the cation mole fraction
        x_cat = self.get_amount(cation, "fraction")

        # calculate the kinematic viscosity
        nu = math.log(nu_w * MW_w / MW) + 15 * x_cat**2 + x_cat**3 * G_123 + 3 * x_cat * G_23 * (1 - 0.05 * x_cat)

        return math.exp(nu) * unit.Quantity("m**2 / s")

    # TODO - need tests of conductivity
    @property
    def conductivity(self):
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
        This approach is used in PHREEQC and Aqion models [aq]_ [hc]_

        .. math::

            EC = {F^2 \\over R T} \\sum_i D_i z_i ^ 2 \\gamma_i ^ {\\alpha} m_i

        Where:

        .. math::

            \\alpha =
            \\begin{cases}
                {\\frac{0.6}{\\sqrt{| z_{i} | }}} & {I < 0.36 | z_{i} | }
                {\\frac{\\sqrt{I}}{| z_i |}} & otherwise
            \\end{cases}

        Note: PHREEQC uses the molal rather than molar concentration according to
        http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/phreeqc3-html/phreeqc3-43.htm

        References
        ----------
        .. [aq] https://www.aqion.de/site/electrical-conductivity
        .. [hc] http://www.hydrochemistry.eu/exmpls/sc.html

        See Also
        --------
        :py:attr:`ionic_strength`
        :py:meth:`get_molar_conductivity()`
        :py:meth:`get_activity_coefficient()`

        """
        EC = 0 * unit.Quantity("S/m")

        for item in self.components:
            z = abs(z=self.get_property(item, "charge"))
            # ignore uncharged species
            if z != 0:
                # determine the value of the exponent alpha
                if self.ionic_strength.magnitude < 0.36 * z:
                    alpha = 0.6 / z**0.5
                else:
                    alpha = self.ionic_strength.magnitude**0.5 / z

                diffusion_coefficient = self.get_property(item, "transport.diffusion_coefficient")

                molar_cond = (
                    diffusion_coefficient
                    * (unit.e * unit.N_A) ** 2
                    * self.get_property(item, "charge") ** 2
                    / (unit.R * self.temperature)
                )

                EC += molar_cond * self.get_activity_coefficient(item) ** alpha * self.get_amount(item, "mol/L")

        return EC.to("S/m")

    @property
    def ionic_strength(self) -> Quantity:
        """
        Return the ionic strength of the solution.

        Return the ionic strength of the solution, calculated as 1/2 * sum ( molality * charge ^2) over all the ions.

        Molal (mol/kg) scale concentrations are used for compatibility with the activity correction formulas.

        Returns
        -------
        Quantity :
            The ionic strength of the parent solution, mol/kg.

        See Also:
        --------
        :py:meth:`get_activity`
        :py:meth:`get_water_activity`

        Notes
        -----
        The ionic strength is calculated according to:

        .. math:: I = \\sum_i m_i z_i^2

        Where :math:`m_i` is the molal concentration and :math:`z_i` is the charge on species i.

        Examples:
        --------
        >>> s1 = pyEQL.Solution([['Na+','0.2 mol/kg'],['Cl-','0.2 mol/kg']])
        >>> s1.ionic_strength
        <Quantity(0.20000010029672785, 'mole / kilogram')>

        >>> s1 = pyEQL.Solution([['Mg+2','0.3 mol/kg'],['Na+','0.1 mol/kg'],['Cl-','0.7 mol/kg']],temperature='30 degC')
        >>> s1.ionic_strength
        <Quantity(1.0000001004383303, 'mole / kilogram')>
        """
        ionic_strength = 0
        for solute in self.components:
            ionic_strength += 0.5 * self.get_amount(solute, "mol/kg") * self.get_property(solute, "charge") ** 2

        return ionic_strength

    @property
    def charge_balance(self) -> float:
        """
        Return the charge balance of the solution.

        Return the charge balance of the solution. The charge balance represents the net electric charge
        on the solution and SHOULD equal zero at all times, but due to numerical errors will usually
        have a small nonzero value. It is calculated according to:

        .. math:: CB = F \\sum_i n_i z_i

        where :math:`n_i` is the number of moles, :math:`z_i` is the charge on species i, and :math:`F` is the Faraday constant.

        Returns
        -------
        float :
            The charge balance of the solution, in equivalents.

        """
        charge_balance = 0
        F = (unit.e * unit.N_A).magnitude
        for solute in self.components:
            charge_balance += self.get_amount(solute, "mol").magnitude * self.get_property(solute, "charge") * F

        return charge_balance

    # TODO - need tests for alkalinity
    @property
    def alkalinity(self):
        """
        Return the alkalinity or acid neutralizing capacity of a solution.

        Returns
        -------
        Quantity :
            The alkalinity of the solution in mg/L as CaCO3

        Notes
        -----
        The alkalinity is calculated according to [stm]_

        .. math::   Alk = F \\sum_{i} z_{i} C_{B} - \\sum_{i} z_{i} C_{A}

        Where :math:`C_{B}` and :math:`C_{A}` are conservative cations and anions, respectively
        (i.e. ions that do not participate in acid-base reactions), and :math:`z_{i}` is their charge.
        In this method, the set of conservative cations is all Group I and Group II cations, and the
        conservative anions are all the anions of strong acids.

        References
        ----------

        .. [stm] Stumm, Werner and Morgan, James J. Aquatic Chemistry, 3rd ed, pp 165. Wiley Interscience, 1996.

        """
        alkalinity = 0 * unit.Quantity("mol/L")
        equiv_wt_CaCO3 = 100.09 / 2 * unit.Quantity("g/mol")

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
                z = self.get_property(item, "charge")
                alkalinity += self.get_amount(item, "mol/L") * z
            if item in acid_anions:
                z = self.get_property(item, "charge")
                alkalinity -= self.get_amount(item, "mol/L") * z

        # convert the alkalinity to mg/L as CaCO3
        return (alkalinity * equiv_wt_CaCO3).to("mg/L")

    @property
    def hardness(self):
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
        hardness = 0 * unit.Quantity("mol/L")
        equiv_wt_CaCO3 = 100.09 / 2 * unit.Quantity("g/mol")

        for item in self.components:
            z = self.get_property(item, "charge")
            if z > 1:
                hardness += z * self.get_amount(item, "mol/L")

        # convert the hardness to mg/L as CaCO3
        return (hardness * equiv_wt_CaCO3).to("mg/L")

    @property
    def debye_length(self) -> Quantity:
        """
        Return the Debye length of a solution.

        Debye length is calculated as [wk3]_

        .. math::

            \\kappa^{-1} = \\sqrt({\\epsilon_r \\epsilon_o k_B T \\over (2 N_A e^2 I)})

        where :math:`I` is the ionic strength, :math:`\\epsilon_r` and :math:`\\epsilon_r`
        are the relative permittivity and vacuum permittivity, :math:`k_B` is the
        Boltzmann constant, and :math:`T` is the temperature, :math:`e` is the
        elementary charge, and :math:`N_A` is Avogadro's number.

        Returns The Debye length, in nanometers.

        References
        .. [wk3] https://en.wikipedia.org/wiki/Debye_length#Debye_length_in_an_electrolyte

        See Also:
            :attr:`ionic_strength`
            :attr:`dielectric_constant`

        """
        # to preserve dimensionality, convert the ionic strength into mol/L units
        ionic_strength = self.ionic_strength.magnitude * unit.Quantity("mol/L")
        dielectric_constant = self.dielectric_constant

        debye_length = (
            dielectric_constant
            * unit.epsilon_0
            * unit.k
            * self.temperature
            / (2 * unit.N_A * unit.e**2 * ionic_strength)
        ) ** 0.5

        return debye_length.to("nm")

    @property
    def bjerrum_length(self) -> Quantity:
        """
        Return the Bjerrum length of a solution.

        Bjerrum length represents the distance at which electrostatic
        interactions between particles become comparable in magnitude
        to the thermal energy.:math:`\\lambda_B` is calculated as

        .. math::

            \\lambda_B = {e^2 \\over (4 \\pi \\epsilon_r \\epsilon_o k_B T)}

        where :math:`e` is the fundamental charge, :math:`\\epsilon_r` and :math:`\\epsilon_r`
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
        https://en.wikipedia.org/wiki/Bjerrum_length

        Examples
        --------
        >>> s1 = pyEQL.Solution()
        >>> s1.bjerrum_length
        <Quantity(0.7152793009386953, 'nanometer')>

        See Also
        --------
        :attr:`dielectric_constant`

        """
        bjerrum_length = unit.e**2 / (
            4 * math.pi * self.dielectric_constant * unit.epsilon_0 * unit.k * self.temperature
        )
        return bjerrum_length.to("nm")

    def get_osmotic_pressure(self):
        """
        Return the osmotic pressure of the solution relative to pure water.

        Returns
            The osmotic pressure of the solution relative to pure water in Pa

        See Also:
            get_water_activity
            get_osmotic_coefficient
            get_salt

        Notes:
            Osmotic pressure is calculated based on the water activity [sata]_ [wk]_

            .. math:: \\Pi = \\frac{RT}{V_{w}} \\ln a_{w}

            Where :math:`\\Pi` is the osmotic pressure, :math:`V_{w}` is the partial
            molar volume of water (18.2 cm**3/mol), and :math:`a_{w}` is the water
            activity.

        References
            .. [sata] Sata, Toshikatsu. Ion Exchange Membranes: Preparation, Characterization, and Modification.
                Royal Society of Chemistry, 2004, p. 10.

            .. [wk] http://en.wikipedia.org/wiki/Osmotic_pressure#Derivation_of_osmotic_pressure

        Examples:
            >>> s1=pyEQL.Solution()
            >>> s1.get_osmotic_pressure()
            0.0

            >>> s1 = pyEQL.Solution([['Na+','0.2 mol/kg'],['Cl-','0.2 mol/kg']])
            >>> soln.get_osmotic_pressure()
            <Quantity(906516.7318131207, 'pascal')>
        """
        # TODO - tie this into parameter() and solvent() objects
        partial_molar_volume_water = 1.82e-5 * unit.Quantity("m ** 3/mol")

        osmotic_pressure = (
            -1 * unit.R * self.temperature / partial_molar_volume_water * math.log(self.get_water_activity())
        )
        logger.info(
            f"Computed osmotic pressure of solution as {osmotic_pressure} Pa at T= {self.temperature} degrees C"
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

        Examples:
        --------

        Todo:

        """
        try:
            if activity is True:
                return -1 * math.log10(self.get_activity(solute))
            if activity is False:
                return -1 * math.log10(self.get_amount(solute, "mol/L").magnitude)
        # if the solute has zero concentration, the log will generate a ValueError
        except ValueError:
            return 0

    def get_amount(self, solute, units):
        """
        Return the amount of 'solute' in the parent solution.

        The amount of a solute can be given in a variety of unit types.
        1. substance per volume (e.g., 'mol/L')
        2. substance per mass of solvent (e.g., 'mol/kg')
        3. mass of substance (e.g., 'kg')
        4. moles of substance ('mol')
        5. mole fraction ('fraction')
        6. percent by weight (%)
        7. number of molecules ('count')

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
            moles = unit.Quantity(self.components[solute], "mol")
        # if the solute is not present in the solution, we'll get a KeyError
        # In that case, the amount is zero
        except KeyError:
            try:
                return 0 * unit.Quantity(units)
            except DimensionalityError:
                logger.error("Unsupported unit specified for get_amount")
                return 0

        # with pint unit conversions enabled, we just pass the unit to pint
        # the logic tests here ensure that only the required arguments are
        # passed to pint for the unit conversion. This avoids unnecessary
        # function calls.
        if units == "count":
            return round((moles * unit.N_A).to("dimensionless"), 0)
        if units == "fraction":
            return moles / (self.get_moles_solvent() + self.get_total_moles_solute())
        mw = unit.Quantity(self.get_property(solute, "molecular_weight")).to("g/mol")
        if units == "%":
            return moles.to("kg", "chem", mw=mw) / self.mass.to("kg") * 100
        if unit.Quantity(units).check("[substance]"):
            return moles.to(units)
        qty = unit.Quantity(units)
        if qty.check("[substance]/[length]**3") or qty.check("[mass]/[length]**3"):
            return moles.to(units, "chem", mw=mw, volume=self.volume)
        if qty.check("[substance]/[mass]") or qty.check("[mass]/[mass]"):
            return moles.to(units, "chem", mw=mw, solvent_mass=self.get_solvent_mass())
        if qty.check("[mass]"):
            return moles.to(units, "chem", mw=mw)

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
        from pymatgen.core import Element

        el = str(Element(element))

        TOT = 0 * unit.Quantity(units)

        # loop through all the solutes, process each one containing element
        for item in self.components:
            # check whether the solute contains the element
            # if ch.contains(item, element):
            if el in self.get_property(item, "elements"):
                # start with the amount of the solute in the desired units
                amt = self.get_amount(item, units)
                ion = Ion.from_formula(item)

                # convert the solute amount into the amount of element by
                # either the mole / mole or weight ratio
                if unit.Quantity(units).dimensionality in (
                    "[substance]",
                    "[substance]/[length]**3",
                    "[substance]/[mass]",
                ):
                    TOT += amt * ion.get_el_amt_dict[el]  # returns {el: mol per formula unit}

                elif unit.Quantity(units).dimensionality in (
                    "[mass]",
                    "[mass]/[length]**3",
                    "[mass]/[mass]",
                ):
                    TOT += amt * ion.to_weight_dict["el"]  # returns {el: wt fraction}

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
        """
        # if units are given on a per-volume basis,
        # iteratively solve for the amount of solute that will preserve the
        # original volume and result in the desired concentration
        if unit.Quantity(amount).dimensionality in (
            "[substance]/[length]**3",
            "[mass]/[length]**3",
        ):
            # store the original volume for later
            orig_volume = self.get_volume()

            # change the amount of the solute present to match the desired amount
            self.components[solute] += (
                unit.Quantity(amount)
                .to(
                    "moles",
                    "chem",
                    mw=unit.Quantity(self.get_property(solute, "molecular_weight")),
                    volume=self.volume,
                    solvent_mass=self.get_solvent_mass(),
                )
                .magnitude
            )

            # set the amount to zero and log a warning if the desired amount
            # change would result in a negative concentration
            if self.get_amount(solute, "mol").magnitude < 0:
                logger.warning(
                    "Attempted to set a negative concentration for solute %s. Concentration set to 0" % solute
                )
                self.set_amount(solute, "0 mol")

            # calculate the volume occupied by all the solutes
            solute_vol = self._get_solute_volume()

            # determine the volume of solvent that will preserve the original volume
            target_vol = orig_volume - solute_vol

            # adjust the amount of solvent
            # volume in L, density in kg/m3 = g/L
            target_mass = target_vol.magnitude * self.water_substance.rho * unit.Quantity("1 g")

            mw = unit.Quantity(self.get_property(self.solvent, "molecular_weight"))
            target_mol = target_mass / mw
            self.components[self.solvent] = target_mol.magnitude

        else:
            # change the amount of the solute present
            self.components[solute] += (
                unit.Quantity(amount)
                .to(
                    "moles",
                    "chem",
                    mw=unit.Quantity(self.get_property(solute, "molecular_weight")),
                    volume=self.volume,
                    solvent_mass=self.get_solvent_mass(),
                )
                .magnitude
            )

            # set the amount to zero and log a warning if the desired amount
            # change would result in a negative concentration
            if self.get_amount(solute, "mol").magnitude < 0:
                logger.warning(
                    "Attempted to set a negative concentration for solute %s. Concentration set to 0" % solute
                )
                self.set_amount(solute, "0 mol")

            # update the volume to account for the space occupied by all the solutes
            # make sure that there is still solvent present in the first place
            if self.get_solvent_mass() <= unit.Quantity("0 kg"):
                logger.error("All solvent has been depleted from the solution")
                return

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

        """
        # raise an error if a negative amount is specified
        if unit.Quantity(amount).magnitude < 0:
            logger.error("Negative amount specified for solute %s. Concentration not changed." % solute)

        # if units are given on a per-volume basis,
        # iteratively solve for the amount of solute that will preserve the
        # original volume and result in the desired concentration
        elif unit.Quantity(amount).dimensionality in (
            "[substance]/[length]**3",
            "[mass]/[length]**3",
        ):
            # store the original volume for later
            orig_volume = self.get_volume()

            # change the amount of the solute present to match the desired amount
            self.components[solute] = (
                unit.Quantity(amount)
                .to(
                    "moles",
                    "chem",
                    mw=unit.Quantity(self.get_property(solute, "molecular_weight")),
                    volume=self.volume,
                    solvent_mass=self.get_solvent_mass(),
                )
                .magnitude
            )

            # calculate the volume occupied by all the solutes
            solute_vol = self._get_solute_volume()

            # determine the volume of solvent that will preserve the original volume
            target_vol = orig_volume - solute_vol

            # adjust the amount of solvent
            target_mass = target_vol.magnitude / 1000 * self.water_substance.rho * unit.Quantity("1 kg")
            mw = self.get_property(self.solvent, "molecular_weight")
            target_mol = target_mass / mw
            self.components[self.solvent] = target_mol.to("mol").magnitude

        else:
            # change the amount of the solute present
            self.components[solute] = (
                unit.Quantity(amount)
                .to(
                    "moles",
                    "chem",
                    mw=unit.Quantity(self.get_property(solute, "molecular_weight")),
                    volume=self.volume,
                    solvent_mass=self.get_solvent_mass(),
                )
                .magnitude
            )

            # update the volume to account for the space occupied by all the solutes
            # make sure that there is still solvent present in the first place
            if self.get_solvent_mass() <= unit.Quantity("0 kg"):
                logger.error("All solvent has been depleted from the solution")
                return

            self._update_volume()

    def get_osmolarity(self, activity_correction=False):
        """Return the osmolarity of the solution in Osm/L.

        Parameters
        ----------
        activity_correction : bool
            If TRUE, the osmotic coefficient is used to calculate the
            osmolarity. This correction is appropriate when trying to predict
            the osmolarity that would be measured from e.g. freezing point
            depression. Defaults to FALSE if omitted.
        """
        factor = self.get_osmotic_coefficient() if activity_correction is True else 1
        return factor * self.get_total_moles_solute() / self.get_volume().to("L")

    def get_osmolality(self, activity_correction=False):
        """Return the osmolality of the solution in Osm/kg.

        Parameters
        ----------
        activity_correction : bool
            If TRUE, the osmotic coefficient is used to calculate the
            osmolarity. This correction is appropriate when trying to predict
            the osmolarity that would be measured from e.g. freezing point
            depression. Defaults to FALSE if omitted.
        """
        factor = self.get_osmotic_coefficient() if activity_correction is True else 1
        return factor * self.get_total_moles_solute() / self.get_solvent_mass().to("kg")

    def get_total_moles_solute(self) -> Quantity:
        """Return the total moles of all solute in the solution."""
        tot_mol = 0
        for item in self.components:
            if item != self.solvent:
                tot_mol += self.components[item]
        return unit.Quantity(tot_mol, "mol")

    def get_moles_solvent(self) -> Quantity:
        """
        Return the moles of solvent present in the solution.

        Returns
            The moles of solvent in the solution.

        """
        return self.get_amount(self.solvent, "mol")

    def get_salt(self):
        """
        Determine the predominant salt in a solution of ions.

        Many empirical equations for solution properties such as activity coefficient,
        partial molar volume, or viscosity are based on the concentration of
        single salts (e.g., NaCl). When multiple ions are present (e.g., a solution
        containing Na+, Cl-, and Mg+2), it is generally not possible to directly model
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
        :py:meth:`get_activity`
        :py:meth:`get_activity_coefficient`
        :py:meth:`get_water_activity`
        :py:meth:`get_osmotic_coefficient`
        :py:meth:`get_osmotic_pressure`
        :py:meth:`get_viscosity_kinematic`

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
        containing Na+, Cl-, and Mg+2), it is generally not possible to directly model
        these quantities.

        The get_salt_list() method examines the ionic composition of a solution and
        simplifies it into a list of salts. The method returns a dictionary of
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

        See Also:
        --------
        :py:meth:`get_activity`
        :py:meth:`get_activity_coefficient`
        :py:meth:`get_water_activity`
        :py:meth:`get_osmotic_coefficient`
        :py:meth:`get_osmotic_pressure`
        :py:meth:`get_viscosity_kinematic`

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

        Returns
            Quantity: the activity coefficient as a dimensionless pint Quantity
        """
        # return unit activity coefficient if the concentration of the solute is zero
        if self.get_amount(solute, "mol").magnitude == 0:
            return unit.Quantity("1 dimensionless")

        try:
            # get the molal-scale activity coefficient from the EOS engine
            molal = self.engine.get_activity_coefficient(solution=self, solute=solute)
        except ValueError:
            logger.warning("Calculation unsuccessful. Returning unit activity coefficient.")
            return unit.Quantity("1 dimensionless")

        # if necessary, convert the activity coefficient to another scale, and return the result
        if scale == "molal":
            return molal
        if scale == "molar":
            total_molality = self.get_total_moles_solute() / self.get_solvent_mass()
            total_molarity = self.get_total_moles_solute() / self.get_volume()
            return (molal * self.water_substance.rho * unit.Quantity("1 g/L") * total_molality / total_molarity).to(
                "dimensionless"
            )
        if scale == "rational":
            return molal * (1 + unit.Quantity("0.018 kg/mol") * self.get_total_moles_solute() / self.get_solvent_mass())

        logger.warning("Invalid scale argument. Returning molal-scale activity coefficient")
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
            solute:
                String representing the name of the solute of interest
            scale:
                The concentration scale for the returned activity.
                Valid options are "molal", "molar", and "rational" (i.e., mole fraction).
                By default, the molal scale activity is returned.
            verbose:
                If True, pyEQL will print a message indicating the parent salt
                that is being used for activity calculations. This option is
                useful when modeling multicomponent solutions. False by default.

        Returns
            The thermodynamic activity of the solute in question (dimensionless)

        Notes:
            The thermodynamic activity depends on the concentration scale used [rs]_ .
            By default, the ionic strength, activity coefficients, and activities are all
            calculated based on the molal (mol/kg) concentration scale.

        References:
            .. [rs] Robinson, R. A.; Stokes, R. H. Electrolyte Solutions: Second Revised
                Edition; Butterworths: London, 1968, p.32.

        See Also:
            :py:meth:`get_activity_coefficient`
            :attr:`ionic_strength`
            :py:meth:`get_salt`

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
            logger.info(f"Calculated {scale} scale activity of solute {solute} as {activity}")

        return activity

    # TODO - engine method
    def get_osmotic_coefficient(self, scale: Literal["molal", "molar", "rational"] = "molal"):
        """
        Return the osmotic coefficient of an aqueous solution.

        The method used depends on the Solution object's equation of state engine.

        """
        molal_phi = self.engine.get_osmotic_coefficient(self)

        if scale == "molal":
            return molal_phi
        if scale == "rational":
            return (
                -molal_phi
                * unit.Quantity("0.018 kg/mol")
                * self.get_total_moles_solute()
                / self.get_solvent_mass()
                / math.log(self.get_amount(self.solvent, "fraction"))
            )
        if scale == "fugacity":
            return math.exp(
                -molal_phi * unit.Quantity("0.018 kg/mol") * self.get_total_moles_solute() / self.get_solvent_mass()
                - math.log(self.get_amount(self.solvent, "fraction"))
            )

        logger.warning("Invalid scale argument. Returning molal-scale osmotic coefficient")
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
        :py:meth:`get_activity_coefficient`
        :attr:`ionic_strength`
        :py:meth:`get_salt`

        Notes
        -----
        Water activity is related to the osmotic coefficient in a solution containing i solutes by:

        .. math:: \\ln a_{w} = - \\Phi M_{w} \\sum_{i} m_{i}

        Where :math:`M_{w}` is the molar mass of water (0.018015 kg/mol) and :math:`m_{i}` is the molal concentration
        of each species.

        If appropriate Pitzer model parameters are not available, the
        water activity is assumed equal to the mole fraction of water.

        References
        ----------
        Blandamer, Mike J., Engberts, Jan B. F. N., Gleeson, Peter T., Reis, Joao Carlos R., 2005. "Activity of
        water in aqueous systems: A frequently neglected property." *Chemical Society Review* 34, 440-458.

        Examples:
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
            logger.warning("Pitzer parameters not found. Water activity set equal to mole fraction")
            return self.get_amount("H2O", "fraction")

        concentration_sum = unit.Quantity("0 mol/kg")
        for item in self.components:
            if item == "H2O":
                pass
            else:
                # TODO - use magnitude instead of quantity; add unit at end of loop
                concentration_sum += self.get_amount(item, "mol/kg")

        logger.info("Calculated water activity using osmotic coefficient")

        return math.exp(-osmotic_coefficient * 0.018015 * unit.Quantity("kg/mol") * concentration_sum) * unit.Quantity(
            "1 dimensionless"
        )

    def get_transport_number(self, solute, activity_correction=False):
        """Calculate the transport number of the solute in the solution.

        Args:
            solute : String identifying the solute for which the transport number is
                to be calculated.

            activity_correction: If True, the transport number will be corrected for activity following
                the same method used for solution conductivity. Defaults to False if omitted.

            Returns
                The transport number of `solute`

            Notes:
                Transport number is calculated according to :

                .. math::

                    t_i = {D_i z_i^2 C_i \\over \\sum D_i z_i^2 C_i}

                Where :math:`C_i` is the concentration in mol/L, :math:`D_i` is the diffusion
                coefficient, and :math:`z_i` is the charge, and the summation extends
                over all species in the solution.

                If `activity_correction` is True, the contribution of each ion to the
                transport number is corrected with an activity factor. See the documentation
                for Solution.conductivity for an explanation of this correction.

            References:
                Geise, G. M.; Cassady, H. J.; Paul, D. R.; Logan, E.; Hickner, M. A. "Specific
                ion effects on membrane potential and the permselectivity of ion exchange membranes.""
                *Phys. Chem. Chem. Phys.* 2014, 16, 21673-21681.

        """
        denominator = unit.Quantity("0  mol / m / s")
        numerator = unit.Quantity("0  mol / m / s")

        for item in self.components:
            z = self.get_property(item, "charge")
            # neutral solutes do not contribute to transport number
            if z == 0:
                continue

            term = self.get_property(item, "transport.diffusion_coefficient") * z**2 * self.get_amount(item, "mol/L")

            if activity_correction is True:
                gamma = self.get_activity_coefficient(item)

                if self.ionic_strength.magnitude < 0.36 * z:
                    alpha = 0.6 / z**0.5
                else:
                    alpha = self.ionic_strength.magnitude**0.5 / z

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
        Calculate the molar (equivalent) conductivity for a solute.

        Args:
            solute: String identifying the solute for which the molar conductivity is
                to be calculated.

        Returns
            The molar or equivalent conductivity of the species in the solution.
            Zero if the solute is not charged.

        Notes:
            Molar conductivity is calculated from the Nernst-Einstein relation [smed]_

            .. math::

                \\kappa_i = {z_i^2 D_i F^2 \\over RT}

            Note that the diffusion coefficient is strongly variable with temperature.

        References:
            .. [smed] Smedley, Stuart. The Interpretation of Ionic Conductivity in Liquids, pp 1-9. Plenum Press, 1980.
        """
        D = self.get_property(solute, "transport.diffusion_coefficient")

        if D is not None:
            molar_cond = (
                D * (unit.e * unit.N_A) ** 2 * self.get_property(solute, "charge") ** 2 / (unit.R * self.temperature)
            )
        else:
            molar_cond = unit.Quantity("0 mS / cm / (mol/L)")

        logger.info(f"Computed molar conductivity as {molar_cond} from D = {D!s} at T={self.temperature}")

        return molar_cond.to("mS / cm / (mol/L)")

    def get_mobility(self, solute):
        """
        Calculate the ionic mobility of the solute.

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
        into an ionic mobility [smed]_

        .. math::

            \\mu_i = {F |z_i| D_i \\over RT}

        References
        ----------
        .. [smed] Smedley, Stuart I. The Interpretation of Ionic Conductivity in Liquids. Plenum Press, 1980.

        """
        D = self.get_property(solute, "transport.diffusion_coefficient")

        mobility = unit.N_A * unit.e * abs(self.get_property(solute, "charge")) * D / (unit.R * self.temperature)

        logger.info(f"Computed ionic mobility as {mobility} from D = {D!s} at T={self.temperature}")

        return mobility.to("m**2/V/s")

    def _get_property(self, solute: str, name: str) -> Optional[Quantity]:
        """Retrieve a thermodynamic property (such as diffusion coefficient)
        for solute, and adjust it from the reference conditions to the conditions
        of the solution.

        Parameters
        ----------
        solute: str
            String representing the chemical formula of the solute species
        name: str
            The name of the property needed, e.g.
            'diffusion coefficient'

        Returns
        -------
        Quantity: The desired parameter or None if not found

        """
        base_temperature = unit.Quantity("25 degC")
        # base_pressure = unit.Quantity("1 atm")

        # query the database using the sanitized formula
        rform = Ion.from_formula(solute).reduced_formula
        # TODO - there seems to be a bug in mongomock / JSONStore wherein properties does
        # not properly return dot-notation fields, e.g. size.molar_volume will not be returned.
        # also $exists:True does not properly return dot notated fields.
        # for now, just set properties=[] to return everything
        # data = list(self.database.query({"formula": rform, name: {"$ne": None}}, properties=["formula", name]))
        data = list(self.database.query({"formula": rform, name: {"$ne": None}}))
        # formulas should always be unique in the database. len==0 indicates no
        # data. len>1 indicates duplicate data.
        if len(data) == 0:
            # try to determine basic properties using pymatgen
            if name == "charge":
                return Ion.from_formula(solute).charge
            if name == "molecular_weight":
                return f"{float(Ion.from_formula(solute).weight)} g/mol"  # weight is a FloatWithUnit

            logger.warning(f"Property {name} for solute {solute} not found in database. Returning None.")
            return None
        if len(data) > 1:
            logger.warning(f"Duplicate database entries for solute {solute} found!")

        data = data[0]

        # perform temperature-corrections or other adjustments for certain
        # parameter types
        if name == "transport.diffusion_coefficient":
            base_value = data["transport"]["diffusion_coefficient"]["value"]

            # correct for temperature and viscosity
            # .. math:: D_1 \over D_2 = T_1 \over T_2 * \mu_2 \over \mu_1
            # where :math:`\mu` is the dynamic viscosity
            # assume that the base viscosity is that of pure water
            return (
                unit.Quantity(base_value)
                * self.temperature
                / base_temperature
                * self.water_substance.mu
                * unit.Quantity("1 Pa*s")
                / self.get_viscosity_dynamic()
            )

        # logger.warning("Diffusion coefficient not found for species %s. Assuming zero." % (solute))
        # return unit.Quantity("0 m**2/s")

        # just return the base-value molar volume for now; find a way to adjust for
        # concentration later
        if name == "size.molar_volume":
            # calculate the partial molar volume for water since it isn't in the database
            if rform == "H2O(aq)":
                vol = (
                    unit.Quantity(self.get_property("H2O", "molecular_weight"))
                    / self.water_substance.rho
                    * unit.Quantity("1 g/L")
                )

                return vol.to("cm **3 / mol")

            base_value = unit.Quantity(data["size"]["molar_volume"]["value"])
            if self.temperature != base_temperature:
                logger.warning("Partial molar volume for species %s not corrected for temperature" % solute)
            return base_value

        if name == "model_parameters.dielectric_zuber":
            return unit.Quantity(data["model_parameters"]["dielectric_zuber"]["value"])

        if name == "model_parameters.activity_pitzer":
            # return a dict
            if data["model_parameters"]["activity_pitzer"].get("Beta0") is not None:
                return data["model_parameters"]["activity_pitzer"]
            return None

        if name == "model_parameters.molar_volume_pitzer":
            # return a dict
            if data["model_parameters"]["molar_volume_pitzer"].get("Beta0") is not None:
                return data["model_parameters"]["molar_volume_pitzer"]
            return None

        # for parameters not named above, just return the base value
        val = data.get(name) if not isinstance(data.get(name), dict) else data[name].get("value")
        # logger.warning("%s has not been corrected for solution conditions" % name)
        if val is not None:
            return unit.Quantity(val)
        return None

    def get_chemical_potential_energy(self, activity_correction=True):
        """
        Return the total chemical potential energy of a solution (not including
        pressure or electric effects).

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
        calculated as follows: [koga]_

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
        .. [koga] Koga, Yoshikata, 2007. *Solution Thermodynamics and its Application to Aqueous Solutions: A differential approach.* Elsevier, 2007, pp. 23-37.

        """
        E = unit.Quantity("0 J")

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
        Calculate the average distance between molecules.

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
        """Recalculate the solution volume based on composition."""
        self.volume = self._get_solvent_volume() + self._get_solute_volume()

    def _get_solvent_volume(self):
        """Return the volume of the pure solvent."""
        # calculate the volume of the pure solvent
        solvent_vol = self.get_solvent_mass() / (self.water_substance.rho * unit.Quantity("1 g/L"))

        return solvent_vol.to("L")

    def _get_solute_volume(self):
        """Return the volume of only the solutes."""
        return self.engine.get_solute_volume(self)

    def copy(self):
        """Return a copy of the solution.

        TODO - clarify whether this is a deep or shallow copy
        """
        # prepare to copy the bulk properties
        new_temperature = str(self.temperature)
        new_pressure = str(self.pressure)
        new_solvent = self.solvent
        new_solvent_mass = str(self.get_solvent_mass())

        # create a list of solutes
        new_solutes = []
        for item in self.components:
            # ignore the solvent
            if item == self.solvent:
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

    def as_dict(self) -> dict:
        """
        Convert the Solution into a dict representation that can be serialized to .json or other format.

        This method is mostly inherited from MSONable
        """
        d = super().as_dict()
        # replace solutes with the current composition
        d["solutes"] = {k: v * unit.Quantity("1 mol") for k, v in self.components.items()}
        # replace the engine with the associated str
        d["engine"] = self._engine
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "Solution":
        """
        Instantiate a Solution from a dictionary generated by as_dict().

        This method is inherited from MSONable
        """
        return super().from_dict(d)

    # informational methods
    def list_solutes(self):
        """List all the solutes in the solution."""
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
                print(item + ":" + "\t {0:0.{decimals}f~}".format(amount, decimals=decimals))
        elif type == "cations":
            print("Cation Concentrations:\n")
            print("========================\n")
            for item in self.components:
                if self.components[item].charge > 0:
                    amount = self.get_amount(item, unit)
                    result_list.append([item, amount])
                    print(item + ":" + "\t {0:0.{decimals}f~}".format(amount, decimals=decimals))
        elif type == "anions":
            print("Anion Concentrations:\n")
            print("========================\n")
            for item in self.components:
                if self.components[item].charge < 0:
                    amount = self.get_amount(item, unit)
                    result_list.append([item, amount])
                    print(item + ":" + "\t {0:0.{decimals}f~}".format(amount, decimals=decimals))

        return result_list

    def list_salts(self, unit="mol/kg", decimals=4):
        list = generate_salt_list(self, unit)
        for item in list:
            print(item.formula + "\t {:0.{decimals}f}".format(list[item], decimals=decimals))

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
        for i in self.components:
            print(i + ":" + "\t {0.magnitude:0.{decimals}f}".format(self.get_activity(i), decimals=decimals))

    def __str__(self):
        # set output of the print() statement for the solution
        str1 = f"Volume: {self.get_volume():.3f~}\n"
        str2 = f"Pressure: {self.pressure:.3f~}\n"
        str3 = f"Temperature: {self.temperature:.3f~}\n"
        str4 = f"Components: {self.list_solutes():}\n"
        return str1 + str2 + str3 + str4

    """
    Legacy methods to be deprecated in a future release.
    """

    @deprecated(
        message="get_solute() is deprecated and will be removed in the next release! Access solutes via the Solution.components attribute and their properties via Solution.get_property(solute, ...)"
    )
    def get_solute(self, i):
        """Return the specified solute object.

        :meta private:
        """
        return self.components[i]

    @deprecated(
        message="get_solvent is deprecated and will be removed in the next release! Use Solution.solvent instead."
    )
    def get_solvent(self):
        """Return the solvent object.

        :meta private:
        """
        return self.components[self.solvent]

    @deprecated(
        message="get_temperature() will be removed in the next release. Access the temperature directly via the property Solution.temperature"
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

        :meta private:
        """
        return self.temperature

    @deprecated(
        message="set_temperature() will be removed in the next release. Set the temperature directly via the property Solution.temperature"
    )
    def set_temperature(self, temperature):
        """
        Set the solution temperature.

        Parameters
        ----------
        temperature : str
            String representing the temperature, e.g. '25 degC'

        :meta private:
        """
        self.temperature = unit.Quantity(temperature)

        # recalculate the volume
        self._update_volume()

    @deprecated(
        message="get_pressure() will be removed in the next release. Access the pressure directly via the property Solution.pressure"
    )
    def get_pressure(self):
        """
        Return the hydrostatic pressure of the solution.

        Returns
        -------
        Quantity: The hydrostatic pressure of the solution, in atm.

        :meta private:
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

        :meta private:
        """
        self._pressure = unit.Quantity(pressure)

    @deprecated(message="get_mass() will be removed in the next release. Use the Solution.mass property instead.")
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

        :meta private:

        """
        return self.mass

    @deprecated(message="get_density() will be removed in the next release. Use the Solution.density property instead.")
    def get_density(self):
        """
        Return the density of the solution.

        Density is calculated from the mass and volume each time this method is called.

        Returns
        -------
        Quantity: The density of the solution.

        :meta private:
        """
        return self.density

    @deprecated(message="get_viscosity_relative() will be removed in the next release.")
    def get_viscosity_relative(self):
        """
        Return the viscosity of the solution relative to that of water.

        This is calculated using a simplified form of the Jones-Dole equation:

        .. math:: \\eta_{rel} = 1 + \\sum_i B_i m_i

        Where :math:`m` is the molal concentration and :math:`B` is an empirical parameter.

        See
        <http://www.nrcresearchpress.com/doi/pdf/10.1139/v77-148>

        :meta private:

        """
        # if self.ionic_strength.magnitude > 0.2:
        #   logger.warning('Viscosity calculation has limited accuracy above 0.2m')

        #        viscosity_rel = 1
        #        for item in self.components:
        #            # ignore water
        #            if item != 'H2O':
        #                # skip over solutes that don't have parameters
        #                try:
        #                    conc = self.get_amount(item,'mol/kg').magnitude
        #                    coefficients= self.get_property(item, 'jones_dole_viscosity')
        #                    viscosity_rel += coefficients[0] * conc ** 0.5 + coefficients[1] * conc + \
        #                    coefficients[2] * conc ** 2
        #                except TypeError:
        #                    continue
        return self.viscosity_dynamic / self.water_substance.mu * unit.Quantity("1 Pa*s")

    @deprecated(
        message="get_viscosity_dynamic() will be removed in the next release. Access directly via the property Solution.viscosity_dynamic."
    )
    def get_viscosity_dynamic(self):
        """
        Return the dynamic (absolute) viscosity of the solution.

        Calculated from the kinematic viscosity

        See Also:
        --------
        get_viscosity_kinematic

        :meta private:
        """
        return self.viscosity_dynamic

    @deprecated(
        message="get_viscosity_kinematic() will be removed in the next release. Access directly via the property Solution.viscosity_kinematic."
    )
    def get_viscosity_kinematic(self):
        """
        Return the kinematic viscosity of the solution.

        Notes
        -----
        The calculation is based on a model derived from the Eyring equation
        and presented by Vásquez-Castillo et al.

        .. math::

            \\ln \\nu = \\ln {\\nu_w MW_w \\over \\sum_i x_i MW_i } +
            15 x_+^2 + x_+^3  \\delta G^*_{123} + 3 x_+ \\delta G^*_{23} (1-0.05x_+)

        Where:

        .. math:: \\delta G^*_{123} = a_o + a_1 (T)^{0.75}
        .. math:: \\delta G^*_{23} = b_o + b_1 (T)^{0.5}

        In which :math:`\\nu` is the kinematic viscosity, MW is the molecular weight,
        `x_+` is the mole fraction of cations, and T is the temperature in degrees C.

        The a and b fitting parameters for a variety of common salts are included in the
        database.

        References
        ----------
        Vásquez-Castillo, G.; Iglesias-Silva, G. a.; Hall, K. R. An extension
        of the McAllister model to correlate kinematic viscosity of electrolyte solutions.
        Fluid Phase Equilib. 2013, 358, 44-49.

        See Also:
        --------
        viscosity_dynamic

        :meta private:

        """
        return self.viscosity_kinematic

    @deprecated(
        message="get_conductivity() will be removed in the next release. Access directly via the property Solution.conductivity."
    )
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
        .. [#] https://www.aqion.de/site/electrical-conductivity
        .. [#] http://www.hydrochemistry.eu/exmpls/sc.html

        See Also:
        --------
        ionic_strength
        get_molar_conductivity()
        get_activity_coefficient()

        :meta private:

        """
        return self.conductivity

    @deprecated(
        replacement=get_amount,
        message="get_mole_fraction() will be removed in the next release. Use get_amount() with units='fraction' instead.",
    )
    def get_mole_fraction(self, solute):
        """
        Return the mole fraction of 'solute' in the solution.

        Notes
        -----
        This function is DEPRECATED.
        Use get_amount() instead and specify 'fraction' as the unit type.

        :meta private:
        """

    @deprecated(
        message="get_ionic_strength() will be removed in the next release. Access directly via the property Solution.ionic_strength"
    )
    def get_ionic_strength(self):
        """
        Return the ionic strength of the solution.

        Return the ionic strength of the solution, calculated as 1/2 * sum ( molality * charge ^2) over all the ions.
        Molal (mol/kg) scale concentrations are used for compatibility with the activity correction formulas.

        Returns
        -------
        Quantity :
            The ionic strength of the parent solution, mol/kg.

        See Also:
        --------
        get_activity
        get_water_activity

        Notes
        -----
        The ionic strength is calculated according to:

        .. math:: I = \\sum_i m_i z_i^2

        Where :math:`m_i` is the molal concentration and :math:`z_i` is the charge on species i.

        Examples:
        --------
        >>> s1 = pyEQL.Solution([['Na+','0.2 mol/kg'],['Cl-','0.2 mol/kg']])
        >>> s1.ionic_strength
        <Quantity(0.20000010029672785, 'mole / kilogram')>

        >>> s1 = pyEQL.Solution([['Mg+2','0.3 mol/kg'],['Na+','0.1 mol/kg'],['Cl-','0.7 mol/kg']],temperature='30 degC')
        >>> s1.ionic_strength
        <Quantity(1.0000001004383303, 'mole / kilogram')>

        :meta private:
        """
        return self.ionic_strength

    @deprecated(
        message="get_charge_balance() will be removed in the next release. Access directly via the property Solution.charge_balance"
    )
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

        :meta private:

        """
        return self.charge_balance

    @deprecated(
        message="get_alkalinity() will be removed in the next release. Access directly via the property Solution.alkalinity"
    )
    def get_alkalinity(self):
        """
        Return the alkalinity or acid neutralizing capacity of a solution.

        Returns
        -------
        Quantity :
            The alkalinity of the solution in mg/L as CaCO3

        Notes
        -----
        The alkalinity is calculated according to:

        .. math:: Alk = F \\sum_i z_i C_B - \\sum_i z_i C_A

        Where :math:`C_B` and :math:`C_A` are conservative cations and anions, respectively
        (i.e. ions that do not participate in acid-base reactions), and :math:`z_i` is their charge.
        In this method, the set of conservative cations is all Group I and Group II cations, and the conservative anions are all the anions of strong acids.

        References
        ----------
        Stumm, Werner and Morgan, James J. Aquatic Chemistry, 3rd ed,
               pp 165. Wiley Interscience, 1996.

        :meta private:

        """
        return self.alkalinity

    @deprecated(
        message="get_hardness() will be removed in the next release. Access directly via the property Solution.hardness"
    )
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

        :meta private:

        """
        return self.hardness

    @deprecated(
        message="get_debye_length() will be removed in the next release. Access directly via the property Solution.debye_length"
    )
    def get_debye_length(self):
        """
        Return the Debye length of a solution.

        Debye length is calculated as

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
        https://en.wikipedia.org/wiki/Debye_length#Debye_length_in_an_electrolyte

        See Also:
        --------
        ionic_strength
        get_dielectric_constant()

        :meta private:

        """
        return self.debye_length

    @deprecated(
        message="get_bjerrum_length() will be removed in the next release. Access directly via the property Solution.bjerrum_length"
    )
    def get_bjerrum_length(self):
        """
        Return the Bjerrum length of a solution.

        Bjerrum length represents the distance at which electrostatic
        interactions between particles become comparable in magnitude
        to the thermal energy.:math:`\\lambda_B` is calculated as

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
        https://en.wikipedia.org/wiki/Bjerrum_length

        Examples:
        --------
        >>> s1 = pyEQL.Solution()
        >>> s1.get_bjerrum_length()
        <Quantity(0.7152793009386953, 'nanometer')>

        See Also:
        --------
        get_dielectric_constant()

        :meta private:

        """
        return self.bjerrum_length

    @deprecated(
        message="get_dielectric_constant() will be removed in the next release. Access directly via the property Solution.dielectric_constant"
    )
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
        Implements the following equation as given by [zub]_

        .. math:: \\epsilon = \\epsilon_solvent \\over 1 + \\sum_i \\alpha_i x_i

        where :math:`\\alpha_i` is a coefficient specific to the solvent and ion, and :math:`x_i`
        is the mole fraction of the ion in solution.


        References
        ----------
        .. [zub] A. Zuber, L. Cardozo-Filho, V.F. Cabral, R.F. Checoni, M. Castier,
        An empirical equation for the dielectric constant in aqueous and nonaqueous
        electrolyte mixtures, Fluid Phase Equilib. 376 (2014) 116-123.
        doi:10.1016/j.fluid.2014.05.037.

        :meta private:
        """
        return self.dielectric_constant
