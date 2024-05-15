"""
pyEQL Solution Class.

:copyright: 2013-2024 by Ryan S. Kingsbury
:license: LGPL, see LICENSE for more details.

"""

from __future__ import annotations

import logging
import math
import os
import warnings
from functools import lru_cache
from importlib.resources import files
from pathlib import Path
from typing import Any, Literal

import numpy as np
from maggma.stores import JSONStore, Store
from monty.dev import deprecated
from monty.json import MontyDecoder, MSONable
from monty.serialization import dumpfn, loadfn
from pint import DimensionalityError, Quantity
from pymatgen.core import Element
from pymatgen.core.ion import Ion

from pyEQL import IonDB, ureg
from pyEQL.activity_correction import _debye_parameter_activity, _debye_parameter_B
from pyEQL.engines import EOS, IdealEOS, NativeEOS, PhreeqcEOS
from pyEQL.salt_ion_match import Salt
from pyEQL.solute import Solute
from pyEQL.utils import FormulaDict, create_water_substance, interpret_units, standardize_formula

EQUIV_WT_CACO3 = ureg.Quantity(100.09 / 2, "g/mol")
# string to denote unknown oxidation states
UNKNOWN_OXI_STATE = "unk"


class Solution(MSONable):
    """
    Class representing the properties of a solution. Instances of this class
    contain information about the solutes, solvent, and bulk properties.
    """

    def __init__(
        self,
        solutes: list[list[str]] | dict[str, str] | None = None,
        volume: str | None = None,
        temperature: str = "298.15 K",
        pressure: str = "1 atm",
        pH: float = 7,
        pE: float = 8.5,
        balance_charge: str | None = None,
        solvent: str | list = "H2O",
        engine: EOS | Literal["native", "ideal", "phreeqc"] = "native",
        database: str | Path | Store | None = None,
        default_diffusion_coeff: float = 1.6106e-9,
        log_level: Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"] | None = "ERROR",
    ):
        """
        Instantiate a Solution from a composition.

        Args:
            solutes: dict, optional. Keys must be the chemical formula, while values must be
                str Quantity representing the amount. For example:

                {"Na+": "0.1 mol/L", "Cl-": "0.1 mol/L"}

                Note that an older "list of lists" syntax is also supported; however this
                will be deprecated in the future and is no longer recommended. The equivalent
                list syntax for the above example is

                [["Na+", "0.1 mol/L"], ["Cl-", "0.1 mol/L"]]

                Defaults to empty (pure solvent) if omitted
            volume: str, optional
                Volume of the solvent, including the unit. Defaults to '1 L' if omitted.
                Note that the total solution volume will be computed using partial molar
                volumes of the respective solutes as they are added to the solution.
            temperature: str, optional
                The solution temperature, including the ureg. Defaults to '25 degC' if omitted.
            pressure: Quantity, optional
                The ambient pressure of the solution, including the unit.
                Defaults to '1 atm' if omitted.
            pH: number, optional
                Negative log of H+ activity. If omitted, the solution will be
                initialized to pH 7 (neutral) with appropriate quantities of
                H+ and OH- ions
            pE: the pE value (redox potential) of the solution.     Lower values = more reducing,
                higher values = more oxidizing. At pH 7, water is stable between approximately
                -7 to +14. The default value corresponds to a pE value typical of natural
                waters in equilibrium with the atmosphere.
            balance_charge: The strategy for balancing charge during init and equilibrium calculations. Valid options
                are 'pH', which will adjust the solution pH to balance charge, 'pE' which will adjust the
                redox equilibrium to balance charge, or the name of a dissolved species e.g. 'Ca+2' or 'Cl-'
                that will be added/subtracted to balance charge. If set to None, no charge balancing will be
                performed either on init or when equilibrate() is called. Note that in this case, equilibrate()
                can distort the charge balance!
            solvent: Formula of the solvent. Solvents other than water are not supported at this time.
            engine: Electrolyte modeling engine to use. See documentation for details on the available engines.
            database: path to a .json file (str or Path) or maggma Store instance that
                contains serialized SoluteDocs. `None` (default) will use the built-in pyEQL database.
            log_level: Log messages of this or higher severity will be printed to stdout. Defaults to 'ERROR', meaning
                that ERROR and CRITICAL messages will be shown, while WARNING, INFO, and DEBUG messages are not. If set to None, nothing will be printed.
            default_diffusion_coeff: Diffusion coefficient value in m^2/s to use in
                calculations when there is no diffusion coefficient for a species in the database. This affects several
                important property calculations including conductivity and transport number, which are related to the
                weighted sums of diffusion coefficients of all species. Setting this argument to zero will exclude any
                species that does not have a tabulated diffusion coefficient from these calculations, possibly resulting
                in underestimation of the conductivity and/or inaccurate transport numbers.

                Missing diffusion coefficients are especially likely in complex electrolytes containing, for example,
                complexes or paired species such as NaSO4[-1]. In such cases, setting default_diffusion_coeff  to zero
                is likely to result in the above errors.

                By default, this argument is set to the diffusion coefficient of NaCl salt, 1.61x10^-9 m2/s.

        Examples:
            >>> s1 = pyEQL.Solution({'Na+': '1 mol/L','Cl-': '1 mol/L'},temperature='20 degC',volume='500 mL')
            >>> print(s1)
            Components:
            Volume: 0.500 l
            Pressure: 1.000 atm
            Temperature: 293.150 K
            Components: ['H2O(aq)', 'H[+1]', 'OH[-1]', 'Na[+1]', 'Cl[-1]']
        """
        # create a logger and attach it to this class
        self.log_level = log_level.upper()
        self.logger = logging.getLogger("pyEQL")
        if self.log_level is not None:
            # set the level of the module logger
            self.logger.setLevel(self.log_level)
            # clear handlers and add a StreamHandler
            self.logger.handlers.clear()
            # use rich for pretty log formatting, if installed
            try:
                from rich.logging import RichHandler

                sh = RichHandler(rich_tracebacks=True)
            except ImportError:
                sh = logging.StreamHandler()
            # the formatter determines what our logs will look like
            formatter = logging.Formatter("[%(asctime)s] [%(levelname)8s] --- %(message)s (%(filename)s:%(lineno)d)")
            sh.setFormatter(formatter)
            self.logger.addHandler(sh)

        # per-instance cache of get_property and other calls that do not depend
        # on composition
        # see https://rednafi.com/python/lru_cache_on_methods/
        self.get_property = lru_cache()(self._get_property)
        self.get_molar_conductivity = lru_cache()(self._get_molar_conductivity)
        self.get_mobility = lru_cache()(self._get_mobility)
        self.default_diffusion_coeff = default_diffusion_coeff
        self.get_diffusion_coefficient = lru_cache()(self._get_diffusion_coefficient)

        # initialize the volume recalculation flag
        self.volume_update_required = False

        # initialize the volume with a flag to distinguish user-specified volume
        if volume is not None:
            # volume_set = True
            self._volume = ureg.Quantity(volume).to("L")
        else:
            # volume_set = False
            self._volume = ureg.Quantity(1, "L")
        # store the initial conditions as private variables in case they are
        # changed later
        self._temperature = ureg.Quantity(temperature)
        self._pressure = ureg.Quantity(pressure)
        self._pE = pE
        self._pH = pH
        self.pE = self._pE
        if isinstance(balance_charge, str) and balance_charge not in ["pH", "pE"]:
            self.balance_charge = standardize_formula(balance_charge)
        else:
            self.balance_charge = balance_charge  #: Standardized formula of the species used for charge balancing.

        # instantiate a water substance for property retrieval
        self.water_substance = create_water_substance(self.temperature, self.pressure)
        """IAPWS instance describing water properties."""

        # create an empty dictionary of components. This dict comprises {formula: moles}
        #  where moles is the number of moles in the solution.
        self.components = FormulaDict({})
        """Special dictionary where keys are standardized formula and values are the moles present in Solution."""

        # connect to the desired property database
        if database is None:
            # load the default database, which is a JSONStore
            db_store = IonDB
        elif isinstance(database, (str, Path)):
            db_store = JSONStore(str(database), key="formula")
            self.logger.debug(f"Created maggma JSONStore from .json file {database}")
        else:
            db_store = database
        self.database = db_store
        """`Store` instance containing the solute property database."""
        self.database.connect()
        self.logger.debug(f"Connected to property database {self.database!s}")

        # set the equation of state engine
        self._engine = engine
        # self.engine: Optional[EOS] = None
        if isinstance(self._engine, EOS):
            self.engine: EOS = self._engine
        elif self._engine == "ideal":
            self.engine = IdealEOS()
        elif self._engine == "native":
            self.engine = NativeEOS()
        elif self._engine == "phreeqc":
            self.engine = PhreeqcEOS()
        else:
            raise ValueError(f'{engine} is not a valid value for the "engine" kwarg!')

        # define the solvent. Allow for list input to support future use of mixed solvents
        if not isinstance(solvent, list):
            solvent = [solvent]
        if len(solvent) > 1:
            raise ValueError("Multiple solvents are not yet supported!")
        if solvent[0] not in ["H2O", "H2O(aq)", "water", "Water", "HOH"]:
            raise ValueError("Non-aqueous solvent detected. These are not yet supported!")
        self.solvent = standardize_formula(solvent[0])
        """Formula of the component that is set as the solvent (currently only H2O(aq) is supported)."""

        # TODO - do I need the ability to specify the solvent mass?
        # # raise an error if the solvent volume has also been given
        # if volume_set is True:
        #     self.logger.error(
        #         "Solvent volume and mass cannot both be specified. Calculating volume based on solvent mass."
        #     )
        # # add the solvent and the mass
        # self.add_solvent(self.solvent, kwargs["solvent"][1])

        # calculate the moles of solvent (water) on the density and the solution volume
        moles = self.volume.magnitude / 55.55  # molarity of pure water
        self.components["H2O"] = moles

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
            msg = (
                'List input of solutes (e.g., [["Na+", "0.5 mol/L]]) is deprecated! Use dictionary formatted input '
                '(e.g., {"Na+":"0.5 mol/L"} instead.)'
            )
            self.logger.warning(msg)
            warnings.warn(msg, DeprecationWarning)
            for item in self._solutes:
                self.add_solute(*item)

        # adjust the charge balance, if necessary
        cb = self.charge_balance
        if not np.isclose(cb, 0, atol=1e-8) and self.balance_charge is not None:
            balanced = False
            self.logger.info(
                f"Solution is not electroneutral (C.B. = {cb} eq/L). Adding {balance_charge} to compensate."
            )
            if self.balance_charge == "pH":
                self.components["H+"] += (
                    -1 * cb * self.volume.to("L").magnitude
                )  # if C.B. is negative, we need to add cations. H+ is 1 eq/mol
                balanced = True
            elif self.balance_charge == "pE":
                raise NotImplementedError("Balancing charge via redox (pE) is not yet implemented!")
            else:
                ions = set().union(*[self.cations, self.anions])  # all ions
                if self.balance_charge not in ions:
                    raise ValueError(
                        f"Charge balancing species {self.balance_charge} was not found in the solution!. "
                        f"Species {ions} were found."
                    )
                z = self.get_property(balance_charge, "charge")
                self.components[balance_charge] += -1 * cb / z * self.volume.to("L").magnitude
                balanced = True

            if not balanced:
                warnings.warn(f"Unable to balance charge using species {self.balance_charge}")

    @property
    def mass(self) -> Quantity:
        """
        Return the total mass of the solution.

        The mass is calculated each time this method is called.

        Returns: The mass of the solution, in kg

        """
        mass = np.sum([self.get_amount(item, "kg").magnitude for item in self.components])
        return ureg.Quantity(mass, "kg")

    @property
    def solvent_mass(self) -> Quantity:
        """
        Return the mass of the solvent.

        This property is used whenever mol/kg (or similar) concentrations
        are requested by get_amount()

        Returns:
            The mass of the solvent, in kg

        See Also:
            :py:meth:`get_amount()`
        """
        return self.get_amount(self.solvent, "kg")

    @property
    def volume(self) -> Quantity:
        """
        Return the volume of the solution.

        Returns:
            Quantity: the volume of the solution, in L
        """
        # if the composition has changed, recalculate the volume first
        if self.volume_update_required is True:
            self._update_volume()
            self.volume_update_required = False

        return self._volume.to("L")

    @volume.setter
    def volume(self, volume: str):
        """Change the total solution volume to volume, while preserving
        all component concentrations.

        Args:
            volume : Total volume of the solution, including the unit, e.g. '1 L'

        Examples:
            >>> mysol = Solution([['Na+','2 mol/L'],['Cl-','0.01 mol/L']],volume='500 mL')
            >>> print(mysol.volume)
            0.5000883925072983 l
            >>> mysol.list_concentrations()
            {'H2O': '55.508435061791985 mol/kg', 'Cl-': '0.00992937605907076 mol/kg', 'Na+': '2.0059345573880325 mol/kg'}
            >>> mysol.volume = '200 mL')
            >>> print(mysol.volume)
            0.2 l
            >>> mysol.list_concentrations()
            {'H2O': '55.50843506179199 mol/kg', 'Cl-': '0.00992937605907076 mol/kg', 'Na+': '2.0059345573880325 mol/kg'}

        """
        # figure out the factor to multiply the old concentrations by
        scale_factor = ureg.Quantity(volume) / self.volume

        # scale down the amount of all the solutes according to the factor
        for solute in self.components:
            self.components[solute] *= scale_factor.magnitude

        # update the solution volume
        self._volume *= scale_factor.magnitude

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
        self._temperature = ureg.Quantity(temperature)

        # update the water substance
        self.water_substance = create_water_substance(self.temperature, self.pressure)

        # recalculate the volume
        self.volume_update_required = True

        # clear any cached solute properties that may depend on temperature
        self.get_property.cache_clear()
        self.get_molar_conductivity.cache_clear()
        self.get_mobility.cache_clear()
        self.get_diffusion_coefficient.cache_clear()

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
        self._pressure = ureg.Quantity(pressure)

        # update the water substance
        self.water_substance = create_water_substance(self.temperature, self.pressure)

        # recalculate the volume
        self.volume_update_required = True

    @property
    def pH(self) -> float | None:
        """Return the pH of the solution."""
        return self.p("H+", activity=False)

    def p(self, solute: str, activity=True) -> float | None:
        """
        Return the negative log of the activity of solute.

        Generally used for expressing concentration of hydrogen ions (pH)

        Args:
            solute : str
                String representing the formula of the solute
            activity: bool, optional
                If False, the function will use the molar concentration rather
                than the activity to calculate p. Defaults to True.

        Returns:
            Quantity
                The negative log10 of the activity (or molar concentration if
                activity = False) of the solute.
        """
        try:
            if activity is True:
                return -1 * math.log10(self.get_activity(solute))
            return -1 * math.log10(self.get_amount(solute, "mol/L").magnitude)
        # if the solute has zero concentration, the log will generate a ValueError
        except ValueError:
            return 0

    @property
    def density(self) -> Quantity:
        """
        Return the density of the solution.

        Density is calculated from the mass and volume each time this method is called.

        Returns:
            Quantity: The density of the solution.
        """
        return self.mass / self.volume

    @property
    def dielectric_constant(self) -> Quantity:
        r"""
        Returns the dielectric constant of the solution.

        Args:
            None

        Returns:
            Quantity: the dielectric constant of the solution, dimensionless.

        Notes:
            Implements the following equation as given by Zuber et al.

            .. math:: \epsilon = \epsilon_{solvent} \over 1 + \sum_i \alpha_i x_i

            where :math:`\alpha_i` is a coefficient specific to the solvent and ion, and :math:`x_i`
            is the mole fraction of the ion in solution.


        References:
            A. Zuber, L. Cardozo-Filho, V.F. Cabral, R.F. Checoni, M. Castier,
            An empirical equation for the dielectric constant in aqueous and nonaqueous
            electrolyte mixtures, Fluid Phase Equilib. 376 (2014) 116-123.
            doi:10.1016/j.fluid.2014.05.037.
        """
        di_water = self.water_substance.epsilon

        denominator = 1
        for item in self.components:
            # ignore water
            if item != "H2O(aq)":
                # skip over solutes that don't have parameters
                # try:
                fraction = self.get_amount(item, "fraction")
                coefficient = self.get_property(item, "model_parameters.dielectric_zuber")
                if coefficient is not None:
                    denominator += coefficient * fraction
                # except TypeError:
                #     self.logger.warning("No dielectric parameters found for species %s." % item)
                # continue

        return ureg.Quantity(di_water / denominator, "dimensionless")

    @property
    def chemical_system(self) -> str:
        """
        Return the chemical system of the Solution as a "-" separated list of elements, sorted alphabetically. For
        example, a solution containing CaCO3 would have a chemical system of "C-Ca-H-O".
        """
        return "-".join(self.elements)

    @property
    def elements(self) -> list:
        """
        Return a list of elements that are present in the solution.

        For example, a solution containing CaCO3 would return ["C", "Ca", "H", "O"]
        """
        els = []
        for s in self.components:
            els.extend(self.get_property(s, "elements"))
        return sorted(set(els))

    @property
    def cations(self) -> dict[str, float]:
        """
        Returns the subset of `components` that are cations.

        The returned dict is sorted by amount in descending order.
        """
        return {k: v for k, v in self.components.items() if self.get_property(k, "charge") > 0}

    @property
    def anions(self) -> dict[str, float]:
        """
        Returns the subset of `components` that are anions.

        The returned dict is sorted by amount in descending order.
        """
        return {k: v for k, v in self.components.items() if self.get_property(k, "charge") < 0}

    @property
    def neutrals(self) -> dict[str, float]:
        """
        Returns the subset of `components` that are neutral (not charged).

        The returned dict is sorted by amount in descending order.
        """
        return {k: v for k, v in self.components.items() if self.get_property(k, "charge") == 0}

    # TODO - need tests for viscosity
    @property
    def viscosity_dynamic(self) -> Quantity:
        """
        Return the dynamic (absolute) viscosity of the solution.

        Calculated from the kinematic viscosity

        See Also:
            :attr:`viscosity_kinematic`
        """
        return self.viscosity_kinematic * self.density

    # TODO - before deprecating get_viscosity_relative, consider whether the Jones-Dole
    # model should be integrated here as a fallback, in case salt parameters for the
    # other model are not available.
    # if self.ionic_strength.magnitude > 0.2:
    #   self.logger.warning('Viscosity calculation has limited accuracy above 0.2m')

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
    #     self.viscosity_dynamic / self.water_substance.mu * ureg.Quantity("1 Pa*s")
    # )
    @property
    def viscosity_kinematic(self) -> Quantity:
        r"""
        Return the kinematic viscosity of the solution.

        Notes:
            The calculation is based on a model derived from the Eyring equation
            and presented in

            .. math::

                \ln \nu = \ln {\nu_w MW_w \over \sum_i x_i MW_i } +
                15 x_+^2 + x_+^3  \delta G^*_{123} + 3 x_+ \delta G^*_{23} (1-0.05x_+)

            Where:

            .. math:: \delta G^*_{123} = a_o + a_1 (T)^{0.75}
            .. math:: \delta G^*_{23} = b_o + b_1 (T)^{0.5}

            In which :math:`\nu` is the kinematic viscosity, MW is the molecular weight,
            :math:`x_{+}` is the mole fraction of cations, and :math:`T` is the temperature in degrees C.

            The a and b fitting parameters for a variety of common salts are included in the
            database.

        References:
            Vásquez-Castillo, G.; Iglesias-Silva, G. a.; Hall, K. R. An extension of the McAllister model to correlate
            kinematic viscosity of electrolyte solutions. Fluid Phase Equilib. 2013, 358, 44-49.

        See Also:
            :py:meth:`viscosity_dynamic`
        """
        # identify the main salt in the solution
        salt = self.get_salt()

        a0 = a1 = b0 = b1 = 0

        # retrieve the parameters for the delta G equations
        params = self.get_property(salt.formula, "model_parameters.viscosity_eyring")
        if params is not None:
            a0 = ureg.Quantity(params["a0"]["value"]).magnitude
            a1 = ureg.Quantity(params["a1"]["value"]).magnitude
            b0 = ureg.Quantity(params["b0"]["value"]).magnitude
            b1 = ureg.Quantity(params["b1"]["value"]).magnitude

            # compute the delta G parameters
            temperature = self.temperature.to("degC").magnitude
            G_123 = a0 + a1 * (temperature) ** 0.75
            G_23 = b0 + b1 * (temperature) ** 0.5
        else:
            # TODO - fall back to the Jones-Dole model! There are currently no eyring parameters in the database!
            # proceed with the coefficients equal to zero and log a warning
            self.logger.warning(f"Viscosity coefficients for {salt.formula} not found. Viscosity will be approximate.")
            G_123 = G_23 = 0

        # get the kinematic viscosity of water, returned by IAPWS in m2/s
        nu_w = self.water_substance.nu

        # compute the effective molar mass of the solution
        total_moles = np.sum([v for k, v in self.components.items()])
        MW = self.mass.to("g").magnitude / total_moles

        # get the MW of water
        MW_w = self.get_property(self.solvent, "molecular_weight").magnitude

        # calculate the cation mole fraction
        # x_cat = self.get_amount(cation, "fraction")
        x_cat = self.get_amount(salt.cation, "fraction").magnitude

        # calculate the kinematic viscosity
        nu = math.log(nu_w * MW_w / MW) + 15 * x_cat**2 + x_cat**3 * G_123 + 3 * x_cat * G_23 * (1 - 0.05 * x_cat)

        return ureg.Quantity(np.exp(nu), "m**2 / s")

    @property
    def conductivity(self) -> Quantity:
        r"""
        Compute the electrical conductivity of the solution.

        Returns:
            The electrical conductivity of the solution in Siemens / meter.

        Notes:
            Conductivity is calculated by summing the molar conductivities of the respective
            solutes.

            .. math::

                EC = {F^2 \over R T} \sum_i D_i z_i ^ 2 m_i = \sum_i \lambda_i m_i

            Where :math:`D_i` is the diffusion coefficient, :math:`m_i` is the molal concentration,
            :math:`z_i` is the charge, and the summation extends over all species in the solution.
            Alternatively, :math:`\lambda_i` is the molar conductivity of solute i.

            Diffusion coefficients :math:`D_i` (and molar conductivities :math:`\lambda_i`) are
            adjusted for the effects of temperature and ionic strength using the method implemented
            in PHREEQC >= 3.4. [aq]_ [hc]_  See `get_diffusion_coefficient for` further details.

        References:
            .. [aq] https://www.aqion.de/site/electrical-conductivity
            .. [hc] http://www.hydrochemistry.eu/exmpls/sc.html

        See Also:
            :py:attr:`ionic_strength`
            :py:meth:`get_diffusion_coefficient`
            :py:meth:`get_molar_conductivity`
        """
        EC = ureg.Quantity(
            np.asarray(
                [
                    self.get_molar_conductivity(i).to("S*L/mol/m").magnitude * self.get_amount(i, "mol/L").magnitude
                    for i in self.components
                ]
            ),
            "S/m",
        )
        return np.sum(EC)

    @property
    def ionic_strength(self) -> Quantity:
        r"""
        Return the ionic strength of the solution.

        Return the ionic strength of the solution, calculated as 1/2 * sum ( molality * charge ^2) over all the ions.

        Molal (mol/kg) scale concentrations are used for compatibility with the activity correction formulas.

        Returns:
            Quantity:
                The ionic strength of the parent solution, mol/kg.

        See Also:
            :py:meth:`get_activity`
            :py:meth:`get_water_activity`

        Notes:
            The ionic strength is calculated according to:

            .. math:: I = \sum_i m_i z_i^2

            Where :math:`m_i` is the molal concentration and :math:`z_i` is the charge on species i.

        Examples:
            >>> s1 = pyEQL.Solution([['Na+','0.2 mol/kg'],['Cl-','0.2 mol/kg']])
            >>> s1.ionic_strength
            <Quantity(0.20000010029672785, 'mole / kilogram')>

            >>> s1 = pyEQL.Solution([['Mg+2','0.3 mol/kg'],['Na+','0.1 mol/kg'],['Cl-','0.7 mol/kg']],temperature='30 degC')
            >>> s1.ionic_strength
            <Quantity(1.0000001004383303, 'mole / kilogram')>
        """
        # compute using magnitudes only, for performance reasons
        ionic_strength = np.sum(
            [mol * self.get_property(solute, "charge") ** 2 for solute, mol in self.components.items()]
        )
        ionic_strength /= self.solvent_mass.to("kg").magnitude  # convert to mol/kg
        ionic_strength *= 0.5
        return ureg.Quantity(ionic_strength, "mol/kg")

    @property
    def charge_balance(self) -> float:
        r"""
        Return the charge balance of the solution.

        Return the charge balance of the solution. The charge balance represents the net electric charge
        on the solution and SHOULD equal zero at all times, but due to numerical errors will usually
        have a small nonzero value. It is calculated according to:

        .. math:: CB = \sum_i C_i z_i

        where :math:`C_i` is the molar concentration, and :math:`z_i` is the charge on species i.

        Returns:
            float :
                The charge balance of the solution, in equivalents (mol of charge) per L.

        """
        charge_balance = 0
        for solute in self.components:
            charge_balance += self.get_amount(solute, "eq/L").magnitude

        return charge_balance

    # TODO - consider adding guard statements to prevent alkalinity from being negative
    @property
    def alkalinity(self) -> Quantity:
        r"""
        Return the alkalinity or acid neutralizing capacity of a solution.

        Returns:
            Quantity: The alkalinity of the solution in mg/L as CaCO3

        Notes:
            The alkalinity is calculated according to [stm]_

            .. math::   Alk = \sum_{i} z_{i} C_{B} + \sum_{i} z_{i} C_{A}

            Where :math:`C_{B}` and :math:`C_{A}` are conservative cations and anions, respectively
            (i.e. ions that do not participate in acid-base reactions), and :math:`z_{i}` is their signed charge.
            In this method, the set of conservative cations is all Group I and Group II cations, and the
            conservative anions are all the anions of strong acids.

        References:
            .. [stm] Stumm, Werner and Morgan, James J. Aquatic Chemistry, 3rd ed, pp 165. Wiley Interscience, 1996.

        """
        alkalinity = ureg.Quantity(0, "mol/L")

        base_cations = {
            "Li[+1]",
            "Na[+1]",
            "K[+1]",
            "Rb[+1]",
            "Cs[+1]",
            "Fr[+1]",
            "Be[+2]",
            "Mg[+2]",
            "Ca[+2]",
            "Sr[+2]",
            "Ba[+2]",
            "Ra[+2]",
        }
        acid_anions = {"Cl[-1]", "Br[-1]", "I[-1]", "SO4[-2]", "NO3[-1]", "ClO4[-1]", "ClO3[-1]"}

        for item in self.components:
            if item in base_cations.union(acid_anions):
                z = self.get_property(item, "charge")
                alkalinity += self.get_amount(item, "mol/L") * z

        # convert the alkalinity to mg/L as CaCO3
        return (alkalinity * EQUIV_WT_CACO3).to("mg/L")

    @property
    def hardness(self) -> Quantity:
        """
        Return the hardness of a solution.

        Hardness is defined as the sum of the equivalent concentrations
        of multivalent cations as calcium carbonate.

        NOTE: at present pyEQL cannot distinguish between mg/L as CaCO3
        and mg/L units. Use with caution.

        Returns:
            Quantity:
                The hardness of the solution in mg/L as CaCO3

        """
        hardness = ureg.Quantity(0, "mol/L")

        for item in self.components:
            z = self.get_property(item, "charge")
            if z > 1:
                hardness += z * self.get_amount(item, "mol/L")

        # convert the hardness to mg/L as CaCO3
        return (hardness * EQUIV_WT_CACO3).to("mg/L")

    @property
    def total_dissolved_solids(self) -> Quantity:
        """
        Total dissolved solids in mg/L (equivalent to ppm) including both charged and uncharged species.

        The TDS is defined as the sum of the concentrations of all aqueous solutes (not including the solvent),
        except for H[+1] and OH[-1]].
        """
        tds = ureg.Quantity(0, "mg/L")
        for s in self.components:
            # ignore pure water and dissolved gases, but not CO2
            if s in ["H2O(aq)", "H[+1]", "OH[-1]"]:
                continue
            tds += self.get_amount(s, "mg/L")

        return tds

    @property
    def TDS(self) -> Quantity:
        """Alias of :py:meth:`total_dissolved_solids`."""
        return self.total_dissolved_solids

    @property
    def debye_length(self) -> Quantity:
        r"""
        Return the Debye length of a solution.

        Debye length is calculated as [wk3]_

        .. math::

            \kappa^{-1} = \sqrt({\epsilon_r \epsilon_o k_B T \over (2 N_A e^2 I)})

        where :math:`I` is the ionic strength, :math:`\epsilon_r` and :math:`\epsilon_r`
        are the relative permittivity and vacuum permittivity, :math:`k_B` is the
        Boltzmann constant, and :math:`T` is the temperature, :math:`e` is the
        elementary charge, and :math:`N_A` is Avogadro's number.

        Returns The Debye length, in nanometers.

        References:
            .. [wk3] https://en.wikipedia.org/wiki/Debye_length#Debye_length_in_an_electrolyte

        See Also:
            :attr:`ionic_strength`
            :attr:`dielectric_constant`

        """
        # to preserve dimensionality, convert the ionic strength into mol/L units
        ionic_strength = ureg.Quantity(self.ionic_strength.magnitude, "mol/L")
        dielectric_constant = self.dielectric_constant

        debye_length = (
            dielectric_constant
            * ureg.epsilon_0
            * ureg.k
            * self.temperature
            / (2 * ureg.N_A * ureg.e**2 * ionic_strength)
        ) ** 0.5

        return debye_length.to("nm")

    @property
    def bjerrum_length(self) -> Quantity:
        r"""
        Return the Bjerrum length of a solution.

        Bjerrum length represents the distance at which electrostatic
        interactions between particles become comparable in magnitude
        to the thermal energy.:math:`\lambda_B` is calculated as

        .. math::

            \lambda_B = {e^2 \over (4 \pi \epsilon_r \epsilon_o k_B T)}

        where :math:`e` is the fundamental charge, :math:`\epsilon_r` and :math:`\epsilon_r`
        are the relative permittivity and vacuum permittivity, :math:`k_B` is the
        Boltzmann constant, and :math:`T` is the temperature.

        Returns:
            Quantity:
                The Bjerrum length, in nanometers.

        References:
            https://en.wikipedia.org/wiki/Bjerrum_length

        Examples:
            >>> s1 = pyEQL.Solution()
            >>> s1.bjerrum_length
            <Quantity(0.7152793009386953, 'nanometer')>

        See Also:
            :attr:`dielectric_constant`

        """
        bjerrum_length = ureg.e**2 / (
            4 * math.pi * self.dielectric_constant * ureg.epsilon_0 * ureg.k * self.temperature
        )
        return bjerrum_length.to("nm")

    @property
    def osmotic_pressure(self) -> Quantity:
        r"""
        Return the osmotic pressure of the solution relative to pure water.

        Returns:
            The osmotic pressure of the solution relative to pure water in Pa

        See Also:
            :attr:`get_water_activity`
            :attr:`get_osmotic_coefficient`
            :attr:`get_salt`

        Notes:
            Osmotic pressure is calculated based on the water activity [sata]_ [wk]_

            .. math:: \Pi = -\frac{RT}{V_{w}} \ln a_{w}

            Where :math:`\Pi` is the osmotic pressure, :math:`V_{w}` is the partial
            molar volume of water (18.2 cm**3/mol), and :math:`a_{w}` is the water
            activity.

        References:
            .. [sata] Sata, Toshikatsu. Ion Exchange Membranes: Preparation, Characterization, and Modification.
                Royal Society of Chemistry, 2004, p. 10.

            .. [wk] http://en.wikipedia.org/wiki/Osmotic_pressure#Derivation_of_osmotic_pressure

        Examples:
            >>> s1=pyEQL.Solution()
            >>> s1.osmotic_pressure
            <Quantity(0.495791416, 'pascal')>

            >>> s1 = pyEQL.Solution([['Na+','0.2 mol/kg'],['Cl-','0.2 mol/kg']])
            >>> soln.osmotic_pressure
            <Quantity(906516.7318131207, 'pascal')>
        """
        partial_molar_volume_water = self.get_property(self.solvent, "size.molar_volume")

        osmotic_pressure = (
            -1 * ureg.R * self.temperature / partial_molar_volume_water * math.log(self.get_water_activity())
        )
        self.logger.debug(
            f"Calculated osmotic pressure of solution as {osmotic_pressure} Pa at T= {self.temperature} degrees C"
        )
        return osmotic_pressure.to("Pa")

    # Concentration  Methods

    def get_amount(self, solute: str, units: str = "mol/L") -> Quantity:
        """
        Return the amount of 'solute' in the parent solution.

        The amount of a solute can be given in a variety of unit types.
        1. substance per volume (e.g., 'mol/L', 'M')
        2. equivalents (i.e., moles of charge) per volume (e.g., 'eq/L', 'meq/L')
        3. substance per mass of solvent (e.g., 'mol/kg', 'm')
        4. mass of substance (e.g., 'kg')
        5. moles of substance ('mol')
        6. mole fraction ('fraction')
        7. percent by weight (%)
        8. number of molecules ('count')
        9. "parts-per-x" units, where ppm = mg/L, ppb = ug/L ppt = ng/L

        Args:
            solute : str
                        String representing the name of the solute of interest
            units : str
                        Units desired for the output. Examples of valid units are
                        'mol/L','mol/kg','mol', 'kg', and 'g/L'
                        Use 'fraction' to return the mole fraction.
                        Use '%' to return the mass percent

        Returns:
            The amount of the solute in question, in the specified units

        See Also:
            :attr:`mass`
            :meth:`add_amount`
            :meth:`set_amount`
            :meth:`get_total_amount`
            :meth:`get_osmolarity`
            :meth:`get_osmolality`
            :meth:`get_total_moles_solute`
            :func:`pyEQL.utils.interpret_units`
        """
        z = 1
        # sanitized unit to be passed to pint
        if "eq" in units:
            _units = units.replace("eq", "mol")
            z = self.get_property(solute, "charge")
            if z == 0:  # uncharged solutes have zero equiv concentration
                return ureg.Quantity(0, _units)
        else:
            _units = interpret_units(units)

        # retrieve the number of moles of solute and its molecular weight
        try:
            moles = ureg.Quantity(self.components[solute], "mol")
        # if the solute is not present in the solution, we'll get a KeyError
        # In that case, the amount is zero
        except KeyError:
            try:
                return ureg.Quantity(0, _units)
            except DimensionalityError:
                self.logger.error(
                    f"Unsupported unit {units} specified for zero-concentration solute {solute}. Returned 0."
                )
                return ureg.Quantity(0, "dimensionless")

        # with pint unit conversions enabled, we just pass the unit to pint
        # the logic tests here ensure that only the required arguments are
        # passed to pint for the unit conversion. This avoids unnecessary
        # function calls.
        if units == "count":
            return round((moles * ureg.N_A).to("dimensionless"), 0)
        if units == "fraction":
            return moles / (self.get_moles_solvent() + self.get_total_moles_solute())
        mw = self.get_property(solute, "molecular_weight").to("g/mol")
        if units == "%":
            return moles.to("kg", "chem", mw=mw) / self.mass.to("kg") * 100
        qty = ureg.Quantity(_units)
        if _units in ["eq", "mol", "moles"] or qty.check("[substance]"):
            return z * moles.to(_units)
        if (
            _units in ["mol/L", "eq/L", "g/L", "mg/L", "ug/L"]
            or qty.check("[substance]/[length]**3")
            or qty.check("[mass]/[length]**3")
        ):
            return z * moles.to(_units, "chem", mw=mw, volume=self.volume)
        if _units in ["mol/kg"] or qty.check("[substance]/[mass]") or qty.check("[mass]/[mass]"):
            return z * moles.to(_units, "chem", mw=mw, solvent_mass=self.solvent_mass)
        if _units in ["kg", "g"] or qty.check("[mass]"):
            return moles.to(_units, "chem", mw=mw)

        raise ValueError(f"Unsupported unit {units} specified for get_amount")

    def get_components_by_element(self) -> dict[str, list]:
        """
        Return a list of all species associated with a given element.

        Elements (keys) are suffixed with their oxidation state in parentheses, e.g.,

        {"Na(1.0)":["Na[+1]", "NaOH(aq)"]}

        Species associated with each element are sorted in descending order of the amount
        present (i.e., the first species listed is the most abundant).
        """
        d = {}
        # by sorting the components according to amount, we ensure that the species
        # are sorted in descending order of concentration in the resulting dict
        for s in self.components:
            # determine the element and oxidation state
            elements = self.get_property(s, "elements")

            for el in elements:
                try:
                    oxi_states = self.get_property(s, "oxi_state_guesses")
                    oxi_state = oxi_states.get(el, UNKNOWN_OXI_STATE)
                except (TypeError, IndexError):
                    self.logger.error(f"No oxidation state found for element {el}. Assigning '{UNKNOWN_OXI_STATE}'")
                    oxi_state = UNKNOWN_OXI_STATE
                key = f"{el}({oxi_state})"
                if d.get(key):
                    d[key].append(s)
                else:
                    d[key] = [s]

        return d

    def get_el_amt_dict(self):
        """
        Return a dict of Element: amount in mol.

        Elements (keys) are suffixed with their oxidation state in parentheses,
        e.g. "Fe(2.0)", "Cl(-1.0)".
        """
        d = {}
        for s, mol in self.components.items():
            elements = self.get_property(s, "elements")
            pmg_ion_dict = self.get_property(s, "pmg_ion")
            oxi_states = self.get_property(s, "oxi_state_guesses")

            for el in elements:
                # stoichiometric coefficient, mol element per mol solute
                stoich = pmg_ion_dict.get(el)
                try:
                    oxi_states = self.get_property(s, "oxi_state_guesses")
                    oxi_state = oxi_states.get(el, UNKNOWN_OXI_STATE)
                except (TypeError, IndexError):
                    self.logger.error(f"No oxidation state found for element {el}. Assigning '{UNKNOWN_OXI_STATE}'")
                    oxi_state = UNKNOWN_OXI_STATE
                key = f"{el}({oxi_state})"
                if d.get(key):
                    d[key] += stoich * mol
                else:
                    d[key] = stoich * mol

        return d

    def get_total_amount(self, element: str, units: str) -> Quantity:
        """
        Return the total amount of 'element' (across all solutes) in the solution.

        Args:
            element: The symbol of the element of interest. The symbol can optionally be followed by the
                oxidation state in parentheses, e.g., "Na(1.0)", "Fe(2.0)", or "O(0.0)". If no oxidation state
                is given, the total concentration of the element (over all oxidation states) is returned.
            units : str
                Units desired for the output. Any unit understood by `get_amount` can be used. Examples of valid
                units are 'mol/L','mol/kg','mol', 'kg', and 'g/L'.

        Returns:
            The total amount of the element in the solution, in the specified units

        See Also:
            :meth:`get_amount`
            :func:`pyEQL.utils.interpret_units`
        """
        TOT: Quantity = 0

        # standardize the element formula and units
        el = str(Element(element.split("(")[0]))
        units = interpret_units(units)

        # enumerate the species whose concentrations we need
        comp_by_element = self.get_components_by_element()

        # compile list of species in different ways depending whether there is an oxidation state
        if "(" in element and UNKNOWN_OXI_STATE not in element:
            ox = float(element.split("(")[-1].split(")")[0])
            key = f"{el}({ox})"
            species = comp_by_element.get(key)
        else:
            species = []
            for k, v in comp_by_element.items():
                if el in k:
                    species.extend(v)

        # loop through the species of interest, adding moles of element
        for item, amt in self.components.items():
            if item in species:
                amt = self.get_amount(item, units)
                ion = Ion.from_formula(item)

                # convert the solute amount into the amount of element by
                # either the mole / mole or weight ratio
                if ureg.Quantity(units).dimensionality in (
                    "[substance]",
                    "[substance]/[length]**3",
                    "[substance]/[mass]",
                ):
                    TOT += amt * ion.get_el_amt_dict()[el]  # returns {el: mol per formula unit}

                elif ureg.Quantity(units).dimensionality in (
                    "[mass]",
                    "[mass]/[length]**3",
                    "[mass]/[mass]",
                ):
                    TOT += amt * ion.to_weight_dict[el]  # returns {el: wt fraction}

        return TOT

    def add_solute(self, formula: str, amount: str):
        """Primary method for adding substances to a pyEQL solution.

        Args:
            formula (str): Chemical formula for the solute. Charged species must contain a + or - and
            (for polyvalent solutes) a number representing the net charge (e.g. 'SO4-2').
            amount (str): The amount of substance in the specified unit system. The string should contain
            both a quantity and a pint-compatible representation of a ureg. e.g. '5 mol/kg' or '0.1 g/L'.
        """
        # if units are given on a per-volume basis,
        # iteratively solve for the amount of solute that will preserve the
        # original volume and result in the desired concentration
        if ureg.Quantity(amount).dimensionality in (
            "[substance]/[length]**3",
            "[mass]/[length]**3",
        ):
            # store the original volume for later
            orig_volume = self.volume

            # add the new solute
            quantity = ureg.Quantity(amount)
            mw = self.get_property(formula, "molecular_weight")  # returns a quantity
            target_mol = quantity.to("moles", "chem", mw=mw, volume=self.volume, solvent_mass=self.solvent_mass)
            self.components[formula] = target_mol.to("moles").magnitude

            # calculate the volume occupied by all the solutes
            solute_vol = self._get_solute_volume()

            # determine the volume of solvent that will preserve the original volume
            target_vol = orig_volume - solute_vol

            # adjust the amount of solvent
            # density is returned in kg/m3 = g/L
            target_mass = target_vol * ureg.Quantity(self.water_substance.rho, "g/L")
            # mw = ureg.Quantity(self.get_property(self.solvent_name, "molecular_weight"))
            mw = self.get_property(self.solvent, "molecular_weight")
            if mw is None:
                raise ValueError(f"Molecular weight for solvent {self.solvent} not found in database. Cannot proceed.")
            target_mol = target_mass.to("g") / mw.to("g/mol")
            self.components[self.solvent] = target_mol.magnitude

        else:
            # add the new solute
            quantity = ureg.Quantity(amount)
            mw = ureg.Quantity(self.get_property(formula, "molecular_weight"))
            target_mol = quantity.to("moles", "chem", mw=mw, volume=self.volume, solvent_mass=self.solvent_mass)
            self.components[formula] = target_mol.to("moles").magnitude

            # update the volume to account for the space occupied by all the solutes
            # make sure that there is still solvent present in the first place
            if self.solvent_mass <= ureg.Quantity(0, "kg"):
                self.logger.error("All solvent has been depleted from the solution")
                return
            # set the volume recalculation flag
            self.volume_update_required = True

    # TODO - deprecate this method. Solvent should be added to the dict like anything else
    # and solvent_name will track which component it is.
    def add_solvent(self, formula: str, amount: str):
        """Same as add_solute but omits the need to pass solvent mass to pint."""
        quantity = ureg.Quantity(amount)
        mw = self.get_property(formula, "molecular_weight")
        target_mol = quantity.to("moles", "chem", mw=mw, volume=self.volume, solvent_mass=self.solvent_mass)
        self.components[formula] = target_mol.to("moles").magnitude

    def add_amount(self, solute: str, amount: str):
        """
        Add the amount of 'solute' to the parent solution.

        Args:
            solute : str
                String representing the name of the solute of interest
            amount : str quantity
                String representing the concentration desired, e.g. '1 mol/kg'
                If the units are given on a per-volume basis, the solution
                volume is not recalculated
                If the units are given on a mass, substance, per-mass, or
                per-substance basis, then the solution volume is recalculated
                based on the new composition

        Returns:
            Nothing. The concentration of solute is modified.
        """
        # if units are given on a per-volume basis,
        # iteratively solve for the amount of solute that will preserve the
        # original volume and result in the desired concentration
        if ureg.Quantity(amount).dimensionality in (
            "[substance]/[length]**3",
            "[mass]/[length]**3",
        ):
            # store the original volume for later
            orig_volume = self.volume

            # change the amount of the solute present to match the desired amount
            self.components[solute] += (
                ureg.Quantity(amount)
                .to(
                    "moles",
                    "chem",
                    mw=self.get_property(solute, "molecular_weight"),
                    volume=self.volume,
                    solvent_mass=self.solvent_mass,
                )
                .magnitude
            )

            # set the amount to zero and log a warning if the desired amount
            # change would result in a negative concentration
            if self.get_amount(solute, "mol").magnitude < 0:
                self.logger.error(
                    "Attempted to set a negative concentration for solute %s. Concentration set to 0" % solute
                )
                self.set_amount(solute, "0 mol")

            # calculate the volume occupied by all the solutes
            solute_vol = self._get_solute_volume()

            # determine the volume of solvent that will preserve the original volume
            target_vol = orig_volume - solute_vol

            # adjust the amount of solvent
            # volume in L, density in kg/m3 = g/L
            target_mass = target_vol * ureg.Quantity(self.water_substance.rho, "g/L")

            mw = self.get_property(self.solvent, "molecular_weight")
            target_mol = target_mass / mw
            self.components[self.solvent] = target_mol.magnitude

        else:
            # change the amount of the solute present
            self.components[solute] += (
                ureg.Quantity(amount)
                .to(
                    "moles",
                    "chem",
                    mw=self.get_property(solute, "molecular_weight"),
                    volume=self.volume,
                    solvent_mass=self.solvent_mass,
                )
                .magnitude
            )

            # set the amount to zero and log a warning if the desired amount
            # change would result in a negative concentration
            if self.get_amount(solute, "mol").magnitude < 0:
                self.logger.error(
                    "Attempted to set a negative concentration for solute %s. Concentration set to 0" % solute
                )
                self.set_amount(solute, "0 mol")

            # update the volume to account for the space occupied by all the solutes
            # make sure that there is still solvent present in the first place
            if self.solvent_mass <= ureg.Quantity(0, "kg"):
                self.logger.error("All solvent has been depleted from the solution")
                return

            # set the volume recalculation flag
            self.volume_update_required = True

    def set_amount(self, solute: str, amount: str):
        """
        Set the amount of 'solute' in the parent solution.

        Args:
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

        Returns:
            Nothing. The concentration of solute is modified.

        """
        # raise an error if a negative amount is specified
        if ureg.Quantity(amount).magnitude < 0:
            raise ValueError(f"Negative amount specified for solute {solute}. Concentration not changed.")

        # if units are given on a per-volume basis,
        # iteratively solve for the amount of solute that will preserve the
        # original volume and result in the desired concentration
        if ureg.Quantity(amount).dimensionality in (
            "[substance]/[length]**3",
            "[mass]/[length]**3",
        ):
            # store the original volume for later
            orig_volume = self.volume

            # change the amount of the solute present to match the desired amount
            self.components[solute] = (
                ureg.Quantity(amount)
                .to(
                    "moles",
                    "chem",
                    mw=ureg.Quantity(self.get_property(solute, "molecular_weight")),
                    volume=self.volume,
                    solvent_mass=self.solvent_mass,
                )
                .magnitude
            )

            # calculate the volume occupied by all the solutes
            solute_vol = self._get_solute_volume()

            # determine the volume of solvent that will preserve the original volume
            target_vol = orig_volume - solute_vol

            # adjust the amount of solvent
            target_mass = target_vol * ureg.Quantity(self.water_substance.rho, "g/L")
            mw = self.get_property(self.solvent, "molecular_weight")
            target_mol = target_mass / mw
            self.components[self.solvent] = target_mol.to("mol").magnitude

        else:
            # change the amount of the solute present
            self.components[solute] = (
                ureg.Quantity(amount)
                .to(
                    "moles",
                    "chem",
                    mw=ureg.Quantity(self.get_property(solute, "molecular_weight")),
                    volume=self.volume,
                    solvent_mass=self.solvent_mass,
                )
                .magnitude
            )

            # update the volume to account for the space occupied by all the solutes
            # make sure that there is still solvent present in the first place
            if self.solvent_mass <= ureg.Quantity(0, "kg"):
                self.logger.critical("All solvent has been depleted from the solution")
                return

            self._update_volume()

    def get_total_moles_solute(self) -> Quantity:
        """Return the total moles of all solute in the solution."""
        tot_mol = 0
        for item in self.components:
            if item != self.solvent:
                tot_mol += self.components[item]
        return ureg.Quantity(tot_mol, "mol")

    def get_moles_solvent(self) -> Quantity:
        """
        Return the moles of solvent present in the solution.

        Returns:
            The moles of solvent in the solution.

        """
        return self.get_amount(self.solvent, "mol")

    def get_osmolarity(self, activity_correction=False) -> Quantity:
        """Return the osmolarity of the solution in Osm/L.

        Args:
            activity_correction : bool
                If TRUE, the osmotic coefficient is used to calculate the
                osmolarity. This correction is appropriate when trying to predict
                the osmolarity that would be measured from e.g. freezing point
                depression. Defaults to FALSE if omitted.
        """
        factor = self.get_osmotic_coefficient() if activity_correction is True else 1
        return factor * self.get_total_moles_solute() / self.volume.to("L")

    def get_osmolality(self, activity_correction=False) -> Quantity:
        """Return the osmolality of the solution in Osm/kg.

        Args:
            activity_correction : bool
                If TRUE, the osmotic coefficient is used to calculate the
                osmolarity. This correction is appropriate when trying to predict
                the osmolarity that would be measured from e.g. freezing point
                depression. Defaults to FALSE if omitted.
        """
        factor = self.get_osmotic_coefficient() if activity_correction is True else 1
        return factor * self.get_total_moles_solute() / self.solvent_mass.to("kg")

    def get_salt(self) -> Salt:
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
        to be calculated (e.g., if a solution contains 0.5 mol/kg of Na+ and Cl-, plus traces
        of H+ and OH-, the matched salt is 0.5 mol/kg NaCl).

        Returns:
            Salt object containing information about the parent salt.

        See Also:
            :py:meth:`get_activity`
            :py:meth:`get_activity_coefficient`
            :py:meth:`get_water_activity`
            :py:meth:`get_osmotic_coefficient`
            :py:attr:`osmotic_pressure`
            :py:attr:`viscosity_kinematic`

        Examples:
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
        d = self.get_salt_dict()
        first_key = next(iter(d.keys()))
        return Salt(d[first_key]["cation"], d[first_key]["anion"])

    # TODO - modify? deprecate? make a salts property?
    def get_salt_dict(self, cutoff: float = 0.01, use_totals: bool = True) -> dict[str, dict]:
        """
        Returns a dict of salts that approximates the composition of the Solution. Like `components`, the dict is
        keyed by formula and the values are the total moles present in the solution, e.g., {"NaCl(aq)": 1}. If the
        Solution is pure water, the returned dict contains only 'HOH'.

        Args:
            cutoff: Lowest salt concentration to consider. Analysis will stop once the concentrations of Salts being
                analyzed goes below this value. Useful for excluding analysis of trace anions.
            use_totals: Whether to base the analysis on total element concentrations or individual species
                concentrations.

        Notes:
            Salts are identified by pairing the predominant cations and anions in the solution, in descending order
            of their respective equivalent amounts.

        Many empirical equations for solution properties such as activity coefficient,
        partial molar volume, or viscosity are based on the concentration of
        single salts (e.g., NaCl). When multiple ions are present (e.g., a solution
        containing Na+, Cl-, and Mg+2), it is generally not possible to directly model
        these quantities.

        The get_salt_dict() method examines the ionic composition of a solution and
        simplifies it into a list of salts. The method returns a dictionary of
        Salt objects where the keys are the salt formulas (e.g., 'NaCl'). The
        Salt object contains information about the stoichiometry of the salt to
        enable its effective concentration to be calculated
        (e.g., 1 M MgCl2 yields 1 M Mg+2 and 2 M Cl-).

        Returns:
            dict
                A dictionary of Salt objects, keyed to the salt formula

        See Also:
            :py:attr:`osmotic_pressure`
            :py:attr:`viscosity_kinematic`
            :py:meth:`get_activity`
            :py:meth:`get_activity_coefficient`
            :py:meth:`get_water_activity`
            :py:meth:`get_osmotic_coefficient`
        """
        """
        Returns a dict of salts that approximates the composition of the Solution. Like `components`, the dict is
        keyed by formula and the values are the total moles of salt present in the solution, e.g., {"NaCl(aq)": 1}

        Notes:
            Salts are identified by pairing the predominant cations and anions in the solution, in descending order
            of their respective equivalent amounts.

        See Also:
            :attr:`components`
            :attr:`cations`
            :attr:`anions`
        """
        salt_dict: dict[str, float] = {}

        if use_totals:
            # # use only the predominant species for each element
            components = {}
            for el, lst in self.get_components_by_element().items():
                components[lst[0]] = self.get_total_amount(el, "mol").magnitude
            # add H+ and OH-, which would otherwise be excluded
            for k in ["H[+1]", "OH[-1]"]:
                if self.components.get(k):
                    components[k] = self.components[k]
        else:
            components = self.components
        components = dict(sorted(components.items(), key=lambda x: x[1], reverse=True))

        # warn if something other than water is the predominant component
        if next(iter(components)) != "H2O(aq)":
            self.logger.warning("H2O(aq) is not the most prominent component in this Solution!")

        # equivalents (charge-weighted moles) of cations and anions
        cations = set(self.cations.keys()).intersection(components.keys())
        anions = set(self.anions.keys()).intersection(components.keys())

        # calculate the charge-weighted (equivalent) concentration of each ion
        cation_equiv = {k: self.get_property(k, "charge") * components[k] for k in cations}
        anion_equiv = {
            k: -1 * self.get_property(k, "charge") * components[k] for k in anions
        }  # make sure amounts are positive

        # sort in descending order of equivalent concentration
        cation_equiv = dict(sorted(cation_equiv.items(), key=lambda x: x[1], reverse=True))
        anion_equiv = dict(sorted(anion_equiv.items(), key=lambda x: x[1], reverse=True))

        len_cat = len(cation_equiv)
        len_an = len(anion_equiv)

        # Only ions are H+ and OH-; return a Salt represnting water (with no amount)
        if len_cat <= 1 and len_an <= 1 and self.solvent == "H2O(aq)":
            x = Salt("H[+1]", "OH[-1]")
            salt_dict.update({x.formula: x.as_dict()})
            salt_dict[x.formula]["mol"] = self.get_amount("H2O", "mol")
            return salt_dict

        # start with the first cation and anion
        index_cat = 0
        index_an = 0

        # list(dict) returns a list of [(key, value), ]
        cation_list = list(cation_equiv.items())
        anion_list = list(anion_equiv.items())

        # calculate the equivalent concentrations of each ion
        c1 = cation_list[index_cat][-1]
        a1 = anion_list[index_an][-1]

        while index_cat < len_cat and index_an < len_an:
            # if the cation concentration is greater, there will be leftover cations
            if c1 > a1:
                # create the salt
                x = Salt(cation_list[index_cat][0], anion_list[index_an][0])
                # there will be leftover cation, so use the anion amount
                salt_dict.update({x.formula: x.as_dict()})
                salt_dict[x.formula]["mol"] = a1 / abs(x.z_anion * x.nu_anion)
                # adjust the amounts of the respective ions
                c1 = c1 - a1
                # move to the next anion
                index_an += 1
                try:
                    a1 = anion_list[index_an][-1]
                    if a1 < cutoff:
                        continue
                except IndexError:
                    continue
            # if the anion concentration is greater, there will be leftover anions
            if c1 < a1:
                # create the salt
                x = Salt(cation_list[index_cat][0], anion_list[index_an][0])
                # there will be leftover anion, so use the cation amount
                salt_dict.update({x.formula: x.as_dict()})
                salt_dict[x.formula]["mol"] = c1 / x.z_cation * x.nu_cation
                # calculate the leftover cation amount
                a1 = a1 - c1
                # move to the next cation
                index_cat += 1
                try:
                    a1 = cation_list[index_cat][-1]
                    if a1 < cutoff:
                        continue
                except IndexError:
                    continue
            if np.isclose(c1, a1):
                # create the salt
                x = Salt(cation_list[index_cat][0], anion_list[index_an][0])
                # there will be nothing leftover, so it doesn't matter which ion you use
                salt_dict.update({x.formula: x.as_dict()})
                salt_dict[x.formula]["mol"] = c1 / x.z_cation * x.nu_cation
                # move to the next cation and anion
                index_an += 1
                index_cat += 1
                try:
                    c1 = cation_list[index_cat][-1]
                    a1 = anion_list[index_an][-1]
                    if (c1 < cutoff) or (a1 < cutoff):
                        continue
                except IndexError:
                    continue

        return salt_dict

    def equilibrate(self, **kwargs) -> None:
        """
        Update the composition of the Solution using the thermodynamic engine.

        Any kwargs specified are passed through to self.engine.equilibrate()

        Returns:
            Nothing. The .components attribute of the Solution is updated.
        """
        self.engine.equilibrate(self, **kwargs)

    # Activity-related methods
    def get_activity_coefficient(
        self,
        solute: str,
        scale: Literal["molal", "molar", "fugacity", "rational"] = "molal",
    ) -> Quantity:
        """
        Return the activity coefficient of a solute in solution.

        The model used to calculate the activity coefficient is determined by the Solution's equation of state
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
            return ureg.Quantity(1, "dimensionless")

        try:
            # get the molal-scale activity coefficient from the EOS engine
            molal = self.engine.get_activity_coefficient(solution=self, solute=solute)
        except (ValueError, ZeroDivisionError):
            self.logger.error("Calculation unsuccessful. Returning unit activity coefficient.", exc_info=True)
            return ureg.Quantity(1, "dimensionless")

        # if necessary, convert the activity coefficient to another scale, and return the result
        if scale == "molal":
            return molal
        if scale == "molar":
            total_molality = self.get_total_moles_solute() / self.solvent_mass
            total_molarity = self.get_total_moles_solute() / self.volume
            return (molal * ureg.Quantity(self.water_substance.rho, "g/L") * total_molality / total_molarity).to(
                "dimensionless"
            )
        if scale == "rational":
            return molal * (1 + ureg.Quantity(0.018015, "kg/mol") * self.get_total_moles_solute() / self.solvent_mass)

        raise ValueError("Invalid scale argument. Pass 'molal', 'molar', or 'rational'.")

    def get_activity(
        self,
        solute: str,
        scale: Literal["molal", "molar", "rational"] = "molal",
    ) -> Quantity:
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

        Returns:
            The thermodynamic activity of the solute in question (dimensionless Quantity)

        Notes:
            The thermodynamic activity depends on the concentration scale used [rs]_ .
            By default, the ionic strength, activity coefficients, and activities are all
            calculated based on the molal (mol/kg) concentration scale.

        References:
            .. [rs] Robinson, R. A.; Stokes, R. H. Electrolyte Solutions: Second Revised
                Edition; Butterworths: London, 1968, p.32.

        See Also:
            :attr:`ionic_strength`
            :py:meth:`get_activity_coefficient`
            :py:meth:`get_salt`

        """
        # switch to the water activity function if the species is H2O
        if solute in ["H2O(aq)", "water", "H2O", "HOH"]:
            activity = self.get_water_activity()
        else:
            # determine the concentration units to use based on the desired scale
            if scale == "molal":
                units = "mol/kg"
            elif scale == "molar":
                units = "mol/L"
            elif scale == "rational":
                units = "fraction"
            else:
                raise ValueError("Invalid scale argument. Pass 'molal', 'molar', or 'rational'.")

            activity = (self.get_activity_coefficient(solute, scale=scale) * self.get_amount(solute, units)).magnitude
            self.logger.debug(f"Calculated {scale} scale activity of solute {solute} as {activity}")

        return ureg.Quantity(activity, "dimensionless")

    # TODO - engine method
    def get_osmotic_coefficient(self, scale: Literal["molal", "molar", "rational"] = "molal") -> Quantity:
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
                * ureg.Quantity(0.018015, "kg/mol")
                * self.get_total_moles_solute()
                / self.solvent_mass
                / math.log(self.get_amount(self.solvent, "fraction"))
            )
        if scale == "fugacity":
            return np.exp(
                -molal_phi * ureg.Quantity(0.018015, "kg/mol") * self.get_total_moles_solute() / self.solvent_mass
                - math.log(self.get_amount(self.solvent, "fraction"))
            ) * ureg.Quantity(1, "dimensionless")

        raise ValueError("Invalid scale argument. Pass 'molal', 'rational', or 'fugacity'.")

    def get_water_activity(self) -> Quantity:
        r"""
        Return the water activity.

        Returns:
            Quantity:
                The thermodynamic activity of water in the solution.

        See Also:
            :attr:`ionic_strength`
            :py:meth:`get_activity_coefficient`
            :py:meth:`get_salt`

        Notes:
            Water activity is related to the osmotic coefficient in a solution containing i solutes by:

            .. math:: \ln a_{w} = - \Phi M_{w} \sum_{i} m_{i}

            Where :math:`M_{w}` is the molar mass of water (0.018015 kg/mol) and :math:`m_{i}` is the molal
            concentration of each species.

            If appropriate Pitzer model parameters are not available, the
            water activity is assumed equal to the mole fraction of water.

        References:
            Blandamer, Mike J., Engberts, Jan B. F. N., Gleeson, Peter T., Reis, Joao Carlos R., 2005. "Activity of
            water in aqueous systems: A frequently neglected property." *Chemical Society Review* 34, 440-458.

        Examples:
            >>> s1 = pyEQL.Solution([['Na+','0.3 mol/kg'],['Cl-','0.3 mol/kg']])
            >>> s1.get_water_activity()
            <Quantity(0.9900944932888518, 'dimensionless')>
        """
        osmotic_coefficient = self.get_osmotic_coefficient()

        if osmotic_coefficient == 1:
            self.logger.warning("Pitzer parameters not found. Water activity set equal to mole fraction")
            return self.get_amount("H2O", "fraction")

        concentration_sum = np.sum([mol for item, mol in self.components.items() if item != "H2O(aq)"])
        concentration_sum /= self.solvent_mass.to("kg").magnitude  # converts to mol/kg

        self.logger.debug("Calculated water activity using osmotic coefficient")

        return ureg.Quantity(np.exp(-osmotic_coefficient * 0.018015 * concentration_sum), "dimensionless")

    def get_chemical_potential_energy(self, activity_correction: bool = True) -> Quantity:
        r"""
        Return the total chemical potential energy of a solution (not including
        pressure or electric effects).

        Args:
            activity_correction : bool, optional
                If True, activities will be used to calculate the true chemical
                potential. If False, mole fraction will be used, resulting in
                a calculation of the ideal chemical potential.

        Returns:
            Quantity
                The actual or ideal chemical potential energy of the solution, in Joules.

        Notes:
            The chemical potential energy (related to the Gibbs mixing energy) is
            calculated as follows: [koga]_

            .. math::      E = R T \sum_i n_i  \ln a_i

            or

            .. math::      E = R T \sum_i n_i \ln x_i

            Where :math:`n` is the number of moles of substance, :math:`T` is the temperature in kelvin,
            :math:`R` the ideal gas constant, :math:`x` the mole fraction, and :math:`a` the activity of
            each component.

            Note that dissociated ions must be counted as separate components,
            so a simple salt dissolved in water is a three component solution (cation,
            anion, and water).

        References:
            .. [koga] Koga, Yoshikata, 2007. *Solution Thermodynamics and its Application to Aqueous Solutions:
            A differential approach.* Elsevier, 2007, pp. 23-37.

        """
        E = ureg.Quantity(0, "J")

        # loop through all the components and add their potential energy
        for item in self.components:
            try:
                if activity_correction is True:
                    E += (
                        ureg.R
                        * self.temperature.to("K")
                        * self.get_amount(item, "mol")
                        * math.log(self.get_activity(item))
                    )
                else:
                    E += (
                        ureg.R
                        * self.temperature.to("K")
                        * self.get_amount(item, "mol")
                        * math.log(self.get_amount(item, "fraction"))
                    )
            # If we have a solute with zero concentration, we will get a ValueError
            except ValueError:
                continue

        return E.to("J")

    def _get_property(self, solute: str, name: str) -> Any | None:
        """Retrieve a thermodynamic property (such as diffusion coefficient)
        for solute, and adjust it from the reference conditions to the conditions
        of the solution.

        Args:
            solute: str
                String representing the chemical formula of the solute species
            name: str
                The name of the property needed, e.g.
                'diffusion coefficient'

        Returns:
            Quantity: The desired parameter or None if not found

        """
        base_temperature = ureg.Quantity(25, "degC")
        # base_pressure = ureg.Quantity("1 atm")

        # query the database using the standardized formula
        rform = standardize_formula(solute)
        # TODO - there seems to be a bug in mongomock / JSONStore wherein properties does
        # not properly return dot-notation fields, e.g. size.molar_volume will not be returned.
        # also $exists:True does not properly return dot notated fields.
        # for now, just set properties=[] to return everything
        # data = list(self.database.query({"formula": rform, name: {"$ne": None}}, properties=["formula", name]))
        data = list(self.database.query({"formula": rform, name: {"$ne": None}}))
        # formulas should always be unique in the database. len==0 indicates no
        # data. len>1 indicates duplicate data.
        if len(data) > 1:
            self.logger.warning(f"Duplicate database entries for solute {solute} found!")
        if len(data) == 0:
            # TODO - add molar volume of water to database?
            if name == "size.molar_volume" and rform == "H2O(aq)":
                # calculate the partial molar volume for water since it isn't in the database
                vol = ureg.Quantity(self.get_property("H2O", "molecular_weight")) / (
                    ureg.Quantity(self.water_substance.rho, "g/L")
                )

                return vol.to("cm **3 / mol")

            # try to determine basic properties using pymatgen
            doc = Solute.from_formula(rform).as_dict()
            data = [doc]

        doc: dict = data[0]

        try:
            # perform temperature-corrections or other adjustments for certain
            # parameter types
            if name == "transport.diffusion_coefficient":
                data = doc["transport"]["diffusion_coefficient"]
                if data is not None:
                    return ureg.Quantity(data["value"]).to("m**2/s")

            # just return the base-value molar volume for now; find a way to adjust for concentration later
            if name == "size.molar_volume":
                data = doc["size"]["molar_volume"]
                if data is not None:
                    base_value = ureg.Quantity(doc["size"]["molar_volume"].get("value"))
                    if self.temperature != base_temperature:
                        self.logger.warning(f"Partial molar volume for species {solute} not corrected for temperature")
                    return base_value
                return data

            if name == "model_parameters.dielectric_zuber":
                return ureg.Quantity(doc["model_parameters"]["dielectric_zuber"]["value"])

            if name == "model_parameters.activity_pitzer":
                # return a dict
                if doc["model_parameters"]["activity_pitzer"].get("Beta0") is not None:
                    return doc["model_parameters"]["activity_pitzer"]
                return None

            if name == "model_parameters.molar_volume_pitzer":
                # return a dict
                if doc["model_parameters"]["molar_volume_pitzer"].get("Beta0") is not None:
                    return doc["model_parameters"]["molar_volume_pitzer"]
                return None

            if name == "molecular_weight":
                return ureg.Quantity(doc.get(name))

            if name == "elements":
                return doc.get(name)

            if name == "oxi_state_guesses":
                # ensure that all oxi states are returned as floats
                return {k: float(v) for k, v in doc.get(name).items()}

            # for parameters not named above, just return the base value
            if name == "pmg_ion" or not isinstance(doc.get(name), dict):
                # if the queried value is not a dict, it is a root level key and should be returned as is
                return doc.get(name)

            val = doc[name].get("value")
            # self.logger.warning("%s has not been corrected for solution conditions" % name)
            if val is not None:
                return ureg.Quantity(val)

        except KeyError:
            self.logger.error(f"Property {name} for solute {solute} not found in database. Returning None.")
            return None

        if name == "model_parameters.molar_volume_pitzer":
            # return a dict
            if doc["model_parameters"]["molar_volume_pitzer"].get("Beta0") is not None:
                return doc["model_parameters"]["molar_volume_pitzer"]
            return None

        if name == "molecular_weight":
            return ureg.Quantity(doc.get(name))

        if name == "oxi_state_guesses":
            return doc.get(name)

        # for parameters not named above, just return the base value
        if name == "pmg_ion" or not isinstance(doc.get(name), dict):
            # if the queried value is not a dict, it is a root level key and should be returned as is
            return doc.get(name)

        val = doc[name].get("value")
        # self.logger.warning("%s has not been corrected for solution conditions" % name)
        if val is not None:
            return ureg.Quantity(val)
        return None

    def get_transport_number(self, solute: str) -> Quantity:
        r"""Calculate the transport number of the solute in the solution.

        Args:
            solute: Formula of the solute for which the transport number is
                to be calculated.

        Returns:
                The transport number of `solute`, as a dimensionless Quantity.

        Notes:
            Transport number is calculated according to :

                .. math::

                    t_i = {D_i z_i^2 C_i \over \sum D_i z_i^2 C_i}

                Where :math:`C_i` is the concentration in mol/L, :math:`D_i` is the diffusion
                coefficient, and :math:`z_i` is the charge, and the summation extends
                over all species in the solution.

                Diffusion coefficients :math:`D_i` are adjusted for the effects of temperature
                and ionic strength using the method implemented in PHREEQC >= 3.4.
                See `get_diffusion_coefficient for` further details.


        References:
                Geise, G. M.; Cassady, H. J.; Paul, D. R.; Logan, E.; Hickner, M. A. "Specific
                ion effects on membrane potential and the permselectivity of ion exchange membranes.""
                *Phys. Chem. Chem. Phys.* 2014, 16, 21673-21681.

        See Also:
            :py:meth:`get_diffusion_coefficient`
            :py:meth:`get_molar_conductivity`
        """
        solute = standardize_formula(solute)
        denominator = numerator = 0

        for item, mol in self.components.items():
            # the molar conductivity of each species is F/RT D * z^2, and the F/RT factor
            # cancels out
            # using species amounts in mol is equivalent to using concentrations in mol/L
            # since there is only one solution volume, and it's much faster.
            term = self.get_molar_conductivity(item).magnitude * mol

            if item == solute:
                numerator = term

            denominator += term

        return ureg.Quantity(numerator / denominator, "dimensionless")

    def _get_molar_conductivity(self, solute: str) -> Quantity:
        r"""
        Calculate the molar (equivalent) conductivity for a solute.

        Args:
            solute: String identifying the solute for which the molar conductivity is
                to be calculated.

        Returns:
            The molar or equivalent conductivity of the species in the solution.
            Zero if the solute is not charged.

        Notes:
            Molar conductivity is calculated from the Nernst-Einstein relation [smed]_

            .. math::

                \lambda_i = \frac{F^2}{RT} D_i z_i^2

            Diffusion coefficients :math:`D_i` are adjusted for the effects of temperature
            and ionic strength using the method implemented in PHREEQC >= 3.4.  See `get_diffusion_coefficient`
            for further details.

        References:
            1. .. [smed] Smedley, Stuart. The Interpretation of Ionic Conductivity in Liquids, pp 1-9. Plenum Press, 1980.

            2. https://www.hydrochemistry.eu/exmpls/sc.html

            3. Appelo, C.A.J. Solute transport solved with the Nernst-Planck equation for concrete pores with `free'
               water and a double layer. Cement and Concrete Research 101, 2017.
               https://dx.doi.org/10.1016/j.cemconres.2017.08.030

            4. CRC Handbook of Chemistry and Physics

        See Also:
            :py:meth:`get_diffusion_coefficient`
        """
        D = self.get_diffusion_coefficient(solute)

        if D != 0:
            molar_cond = (
                D * (ureg.e * ureg.N_A) ** 2 * self.get_property(solute, "charge") ** 2 / (ureg.R * self.temperature)
            )
        else:
            molar_cond = ureg.Quantity(0, "mS / cm / (mol/L)")

        self.logger.debug(f"Calculated molar conductivity as {molar_cond} from D = {D!s} at T={self.temperature}")

        return molar_cond.to("mS / cm / (mol/L)")

    def _get_diffusion_coefficient(self, solute: str, activity_correction: bool = True) -> Quantity:
        r"""
        Get the **temperature-adjusted** diffusion coefficient of a solute.

        Args:
            solute: the solute for which to retrieve the diffusion coefficient.
            activity_correction: If True (default), adjusts the diffusion coefficient for the effects of ionic
                strength using a model from Ref 2.

        Notes:
            This method is equivalent to self.get_property(solute, "transport.diffusion_coefficient")
            ONLY when the Solution temperature is the same as the reference temperature for the diffusion coefficient
            in the database (usually 25 C).

            Otherwise, the reference D value is adjusted based on the Solution temperature and (optionally),
            ionic strength. The adjustments are

            .. math::

                D_T = D_{298} \exp(\frac{d}{T} - \frac{d}{298}) \frac{\nu_{298}}{\nu_T}

            .. math::

                D_{\gamma} = D^0 \exp(\frac{-a1 A |z_i| \sqrt{I}}{1+\kappa a}

            .. math::

                 \kappa a = B \sqrt{I} \frac{a2}{1+I^{0.75}}

            where a1, a2, and d are parameters from Ref. 2, A and B are the parameters used in the Debye Huckel
            equation, and I is the ionic strength. If the model parameters for a particular solute are not available,
            default values of d=0, a1=1.6, and a2=4.73 (as recommended in Ref. 2) are used instead.

        References:
            1. https://www.hydrochemistry.eu/exmpls/sc.html
            2. Appelo, C.A.J. Solute transport solved with the Nernst-Planck equation for concrete pores with `free'
               water and a double layer. Cement and Concrete Research 101, 2017.
               https://dx.doi.org/10.1016/j.cemconres.2017.08.030
            3. CRC Handbook of Chemistry and Physics

        See Also:
            pyEQL.activity_correction._debye_parameter_B
            pyEQL.activity_correction._debye_parameter_activity

        """
        D = self.get_property(solute, "transport.diffusion_coefficient")
        rform = standardize_formula(solute)
        if D is None or D.magnitude == 0:
            self.logger.warning(
                f"Diffusion coefficient not found for species {rform}. Using default value of "
                f"{self.default_diffusion_coeff} m**2/s."
            )
            D = ureg.Quantity(self.default_diffusion_coeff, "m**2/s")

        # assume reference temperature is 298.15 K (this is the case for all current DB entries)
        T_ref = 298.15
        mu_ref = 0.0008900225512925807  # water viscosity from IAPWS97 at 298.15 K
        T_sol = self.temperature.to("K").magnitude
        mu = self.water_substance.mu

        # skip temperature correction if within 1 degree
        if abs(T_sol - T_ref) > 1 or activity_correction is True:
            # get the a1, a2, and d parameters required by the PHREEQC model
            try:
                doc = self.database.query_one({"formula": rform})
                d = doc["model_parameters"]["diffusion_temp_smolyakov"]["d"]["value"]
                a1 = doc["model_parameters"]["diffusion_temp_smolyakov"]["a1"]["value"]
                a2 = doc["model_parameters"]["diffusion_temp_smolyakov"]["a2"]["value"]
                # values will be a str, e.g. "1 dimensionless"
                d = float(d.split(" ")[0])
                a1 = float(a1.split(" ")[0])
                a2 = float(a2.split(" ")[0])
            except TypeError:
                # this means the database doesn't contain a d value.
                # according to Ref 2, the following are recommended default parameters
                self.logger.warning(
                    f"Temperature and ionic strength correction parameters for solute {rform} diffusion "
                    "coefficient not in database. Using recommended default values of a1=1.6, a2=4.73, and d=0."
                )
                d = 0
                a1 = 1.6
                a2 = 4.73

            # use the PHREEQC model from Ref 2 to correct for temperature
            D_final = D * np.exp(d / T_sol - d / T_ref) * mu_ref / mu

            if activity_correction:
                A = _debye_parameter_activity(str(self.temperature)).to("kg**0.5/mol**0.5").magnitude / 2.303
                B = _debye_parameter_B(str(self.temperature)).to("1/angstrom * kg**0.5/mol**0.5").magnitude
                z = self.get_property(solute, "charge")
                IS = self.ionic_strength.magnitude
                kappaa = B * IS**0.5 * a2 / (1 + IS**0.75)
                # correct for ionic strength
                D_final *= np.exp(-a1 * A * abs(z) * IS**0.5 / (1 + kappaa))
            # else:
            #     # per CRC handbook, D increases by 2-3% per degree above 25 C
            #     return D * (1 + 0.025 * (T_sol - T_ref))
        else:
            D_final = D

        return D_final

    def _get_mobility(self, solute: str) -> Quantity:
        r"""
        Calculate the ionic mobility of the solute.

        Args:
            solute (str): String identifying the solute for which the mobility is to be calculated.

        Returns:
            float: The ionic mobility. Zero if the solute is not charged.

        Note:
            This function uses the Einstein relation to convert a diffusion coefficient into an ionic mobility [smed]_

            .. math::

                \mu_i = {F |z_i| D_i \over RT}

        References:
            Smedley, Stuart I. The Interpretation of Ionic Conductivity in Liquids. Plenum Press, 1980.
        """
        D = self.get_diffusion_coefficient(solute)

        mobility = ureg.N_A * ureg.e * abs(self.get_property(solute, "charge")) * D / (ureg.R * self.temperature)

        self.logger.debug(f"Calculated ionic mobility as {mobility} from D = {D!s} at T={self.temperature}")

        return mobility.to("m**2/V/s")

    def get_lattice_distance(self, solute: str) -> Quantity:
        r"""
        Calculate the average distance between molecules.

        Calculate the average distance between molecules of the given solute,
        assuming that the molecules are uniformly distributed throughout the
        solution.

        Args:
            solute : str
                String representing the name of the solute of interest

        Returns:
            Quantity: The average distance between solute molecules

        Examples:
            >>> soln = Solution([['Na+','0.5 mol/kg'],['Cl-','0.5 mol/kg']])
            >>> soln.get_lattice_distance('Na+')
            1.492964.... nanometer

        Notes:
            The lattice distance is related to the molar concentration as follows:

            .. math:: d = ( C_i N_A ) ^ {-{1 \over 3}}

        """
        # calculate the volume per particle as the reciprocal of the molar concentration
        # (times avogadro's number). Take the cube root of the volume to get
        # the average distance between molecules
        distance = (self.get_amount(solute, "mol/L") * ureg.N_A) ** (-1 / 3)

        return distance.to("nm")

    def _update_volume(self):
        """Recalculate the solution volume based on composition."""
        self._volume = self._get_solvent_volume() + self._get_solute_volume()

    def _get_solvent_volume(self):
        """Return the volume of the pure solvent."""
        # calculate the volume of the pure solvent
        solvent_vol = self.solvent_mass / ureg.Quantity(self.water_substance.rho, "g/L")

        return solvent_vol.to("L")

    def _get_solute_volume(self):
        """Return the volume of only the solutes."""
        return self.engine.get_solute_volume(self)

    def as_dict(self) -> dict:
        """Convert the Solution into a dict representation that can be serialized to .json or other format."""
        # clear the volume update flag, if required
        if self.volume_update_required:
            self._update_volume()
        d = super().as_dict()
        for k, v in d.items():
            # convert all Quantity to str
            if isinstance(v, Quantity):
                d[k] = str(v)
        # replace solutes with the current composition
        d["solutes"] = {k: f"{v} mol" for k, v in self.components.items()}
        # replace the engine with the associated str
        d["engine"] = self._engine
        # d["logger"] = self.logger.__dict__
        return d

    @classmethod
    def from_dict(cls, d: dict) -> Solution:
        """Instantiate a Solution from a dictionary generated by as_dict()."""
        # because of the automatic volume updating that takes place during the __init__ process,
        # care must be taken here to recover the exact quantities of solute and volume
        # first we store the volume of the serialized solution
        orig_volume = ureg.Quantity(d["volume"])
        # then instantiate a new one
        decoded = {k: MontyDecoder().process_decoded(v) for k, v in d.items() if not k.startswith("@")}
        new_sol = cls(**decoded)
        # now determine how different the new solution volume is from the original
        scale_factor = (orig_volume / new_sol.volume).magnitude
        # reset the new solution volume to that of the original. In the process of
        # doing this, all the solute amounts are scaled by new_sol.volume / volume
        new_sol.volume = str(orig_volume)
        # undo the scaling by diving by that scale factor
        for sol in new_sol.components:
            new_sol.components[sol] /= scale_factor
        # ensure that another volume update won't be triggered by these changes
        # (this line should in principle be unnecessary, but it doesn't hurt anything)
        new_sol.volume_update_required = False
        return new_sol

    @classmethod
    def from_preset(
        cls, preset: Literal["seawater", "rainwater", "wastewater", "urine", "normal saline", "Ringers lactate"]
    ) -> Solution:
        """Instantiate a solution from a preset composition.

        Args:
            preset (str): String representing the desired solution.
              Valid entries are 'seawater', 'rainwater', 'wastewater',
              'urine', 'normal saline' and 'Ringers lactate'.

        Returns:
            A pyEQL Solution object.

        Raises:
            FileNotFoundError: If the given preset file doesn't exist on the file system.

        Notes:
            The following sections explain the different solution options:

            - 'rainwater' - pure water in equilibrium with atmospheric CO2 at pH 6
            - 'seawater' or 'SW'- Standard Seawater. See Table 4 of the Reference for Composition [1]_
            - 'wastewater' or 'WW' - medium strength domestic wastewater. See Table 3-18 of [2]_
            - 'urine' - typical human urine. See Table 3-15 of [2]_
            - 'normal saline' or 'NS' - normal saline solution used in medicine [3]_
            - 'Ringers lacatate' or 'RL' - Ringer's lactate solution used in medicine [4]_

        References:
            .. [1] Millero, Frank J. "The composition of Standard Seawater and the definition of
                   the Reference-Composition Salinity Scale." *Deep-sea Research. Part I* 55(1), 2008, 50-72.

            .. [2] Metcalf & Eddy, Inc. et al. *Wastewater Engineering: Treatment and Resource Recovery*, 5th Ed.
                   McGraw-Hill, 2013.

            .. [3] https://en.wikipedia.org/wiki/Saline_(medicine)

            .. [4] https://en.wikipedia.org/wiki/Ringer%27s_lactate_solution
        """
        # preset_dir = files("pyEQL") / "presets"
        # Path to the YAML and JSON files corresponding to the preset
        yaml_path = files("pyEQL") / "presets" / f"{preset}.yaml"
        json_path = files("pyEQL") / "presets" / f"{preset}.json"

        # Check if the file exists
        if yaml_path.exists():
            preset_path = yaml_path
        elif json_path.exists():
            preset_path = json_path
        else:
            raise FileNotFoundError(f"Invalid preset! File '{yaml_path}' or '{json_path} not found!")

        # Create and return a Solution object
        return cls().from_file(preset_path)

    def to_file(self, filename: str | Path) -> None:
        """Saving to a .yaml or .json file.

        Args:
            filename (str | Path): The path to the file to save Solution.
              Valid extensions are .json or .yaml.
        """
        str_filename = str(filename)
        if not ("yaml" in str_filename.lower() or "json" in str_filename.lower()):
            self.logger.error("Invalid file extension entered - %s" % str_filename)
            raise ValueError("File extension must be .json or .yaml")
        if "yaml" in str_filename.lower():
            solution_dict = self.as_dict()
            solution_dict.pop("database")
            dumpfn(solution_dict, filename)
        else:
            dumpfn(self, filename)

    @classmethod
    def from_file(self, filename: str | Path) -> Solution:
        """Loading from a .yaml or .json file.

        Args:
            filename (str | Path): Path to the .json or .yaml file (including extension) to load the Solution from.
              Valid extensions are .json or .yaml.

        Returns:
            A pyEQL Solution object.

        Raises:
            FileNotFoundError: If the given filename doesn't exist on the file system.
        """
        if not os.path.exists(filename):
            raise FileNotFoundError(f"File '{filename}' not found!")
        str_filename = str(filename)
        if "yaml" in str_filename.lower():
            true_keys = [
                "solutes",
                "volume",
                "temperature",
                "pressure",
                "pH",
                "pE",
                "balance_charge",
                "solvent",
                "engine",
                # "database",
            ]
            solution_dict = loadfn(filename)
            keys_to_delete = [key for key in solution_dict if key not in true_keys]
            for key in keys_to_delete:
                solution_dict.pop(key)
            return Solution(**solution_dict)
        return loadfn(filename)

    # arithmetic operations
    def __add__(self, other: Solution):
        """
        Solution addition: mix two solutions together.

        Args:
            other: The Solutions to be mixed with this solution.

        Returns:
            A Solution object that represents the result of mixing this solution and other.

        Notes:
            The initial volume of the mixed solution is set as the sum of the volumes of this solution and other.
            The pressure and temperature are volume-weighted averages. The pH and pE values are currently APPROXIMATE
            because they are calculated assuming H+ and e- mix conservatively (i.e., the mixing process does not
            incorporate any equilibration reactions or buffering). Such support is planned in a future release.
        """
        # check to see if the two solutions have the same solvent
        if self.solvent != other.solvent:
            raise ValueError("Cannot add Solution with different solvents!")

        if self._engine != other._engine:
            raise ValueError("Cannot add Solution with different engines!")

        if self.database != other.database:
            raise ValueError("Cannot add Solution with different databases!")

        # set the pressure for the new solution
        p1 = self.pressure
        t1 = self.temperature
        v1 = self.volume
        p2 = other.pressure
        t2 = other.temperature
        v2 = other.volume

        # set the initial volume as the sum of the volumes
        mix_vol = v1 + v2

        # check to see if the solutions have the same temperature and pressure
        if p1 != p2:
            self.logger.info(
                "Adding two solutions of different pressure. Pressures will be averaged (weighted by volume)"
            )

        mix_pressure = (p1 * v1 + p2 * v2) / (mix_vol)

        if t1 != t2:
            self.logger.info(
                "Adding two solutions of different temperature. Temperatures will be averaged (weighted by volume)"
            )

        # do all temperature conversions in Kelvin to avoid ambiguity associated with "offset units". See pint docs.
        mix_temperature = (t1.to("K") * v1 + t2.to("K") * v2) / (mix_vol)

        # retrieve the amount of each component in the parent solution and
        # store in a list.
        mix_species = FormulaDict({})
        for sol, amt in self.components.items():
            mix_species.update({sol: f"{amt} mol"})
        for sol2, amt2 in other.components.items():
            if mix_species.get(sol2):
                orig_amt = float(mix_species[sol2].split(" ")[0])
                mix_species[sol2] = f"{orig_amt+amt2} mol"
            else:
                mix_species.update({sol2: f"{amt2} mol"})

        # TODO - call equilibrate() here once the method is functional to get new pH and pE, instead of the below
        warnings.warn(
            "The pH and pE value of the mixed solution is approximate! More accurate addition (mixing) of"
            "this property is planned for a future release."
        )
        # calculate the new pH and pE (before reactions) by mixing
        mix_pH = -math.log10(float(mix_species["H+"].split(" ")[0]) / mix_vol.to("L").magnitude)

        # pE = -log[e-], so calculate the moles of e- in each solution and mix them
        mol_e_self = 10 ** (-1 * self.pE) * self.volume.to("L").magnitude
        mol_e_other = 10 ** (-1 * other.pE) * other.volume.to("L").magnitude
        mix_pE = -math.log10((mol_e_self + mol_e_other) / mix_vol.to("L").magnitude)

        # create a new solution
        return Solution(
            mix_species.data,  # pass a regular dict instead of the FormulaDict
            volume=str(mix_vol),
            pressure=str(mix_pressure),
            temperature=str(mix_temperature.to("K")),
            pH=mix_pH,
            pE=mix_pE,
        )

    def __sub__(self, other: Solution):
        raise NotImplementedError("Subtraction of solutions is not implemented.")

    def __mul__(self, factor: float):
        """
        Solution multiplication: scale all components by a factor. For example, Solution * 2 will double the moles of
        every component (including solvent). No other properties will change.
        """
        self.volume *= factor
        return self

    def __truediv__(self, factor: float):
        """
        Solution division: scale all components by a factor. For example, Solution / 2 will remove half of the moles
        of every compoonents (including solvent). No other properties will change.
        """
        self.volume /= factor
        return self

    # informational methods

    def print(
        self,
        mode: Literal["all", "ions", "cations", "anions", "neutrals"] = "all",
        units: Literal["ppm", "mol", "mol/kg", "mol/L", "%", "activity"] = "mol",
        places=4,
    ):
        """
        Print details about the Solution.

        Args:
            mode: Whether to list the amounts of all solutes, or only anions, cations, any ion, or any neutral solute.
            units: The units to list solute amounts in. "activity" will list dimensionless activities instead of
                concentrations.
            places: The number of decimal places to round the solute amounts.
        """
        print(self)
        str1 = "Activities" if units == "activity" else "Amounts"
        str2 = f" ({units})" if units != "activity" else ""
        header = f"\nComponent {str1}{str2}:"
        print(header)
        print("=" * (len(header) - 1))
        for i in self.components:
            if mode != "all":
                z = self.get_property(i, "charge")
                if (
                    (z != 0 and mode == "neutrals")
                    or (z >= 0 and mode == "anions")
                    or (z <= 0 and mode == "cations")
                    or (z == 0 and mode == "ions")
                ):
                    continue

            amt = self.get_activity(i).magnitude if units == "activity" else self.get_amount(i, units).magnitude

            print(f"{i}:\t {amt:0.{places}f}")

    def __str__(self):
        # set output of the print() statement for the solution
        l1 = f"Volume: {self.volume:.3f~}"
        l2 = f"Temperature: {self.temperature:.3f~}"
        l3 = f"Pressure: {self.pressure:.3f~}"
        l4 = f"pH: {self.pH:.1f}"
        l5 = f"pE: {self.pE:.1f}"
        l6 = f"Solvent: {self.solvent}"
        l7 = f"Components: {self.list_solutes():}"
        return f"{l1}\n{l2}\n{l3}\n{l4}\n{l5}\n{l6}\n{l7}"

    """
    Legacy methods to be deprecated in a future release.
    """

    @deprecated(
        message="list_salts() is deprecated and will be removed in the next release! Use Solution.get_salt_dict() instead.)"
    )
    def list_salts(self, unit="mol/kg", decimals=4):  # pragma: no cover
        for k, v in self.get_salt_dict().items():
            print(k + "\t {:0.{decimals}f}".format(v, decimals=decimals))

    @deprecated(
        message="list_solutes() is deprecated and will be removed in the next release! Use Solution.components.keys() instead.)"
    )
    def list_solutes(self):  # pragma: no cover
        """List all the solutes in the solution."""
        return list(self.components.keys())

    @deprecated(
        message="list_concentrations() is deprecated and will be removed in the next release! Use Solution.print() instead.)"
    )
    def list_concentrations(self, unit="mol/kg", decimals=4, type="all"):  # pragma: no cover
        """
        List the concentration of each species in a solution.

        Parameters
        ----------
        unit: str
            String representing the desired concentration ureg.
        decimals: int
            The number of decimal places to display. Defaults to 4.
        type     : str
            The type of component to be sorted. Defaults to 'all' for all
            solutes. Other valid arguments are 'cations' and 'anions' which
            return lists of cations and anions, respectively.

        Returns:
        -------
        dict
            Dictionary containing a list of the species in solution paired with their amount in the specified units
        :meta private:
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

    @deprecated(
        message="list_activities() is deprecated and will be removed in the next release! Use Solution.print() instead.)"
    )
    def list_activities(self, decimals=4):  # pragma: no cover
        """
        List the activity of each species in a solution.

        Parameters
        ----------
        decimals: int
            The number of decimal places to display. Defaults to 4.

        Returns:
        -------
        dict
            Dictionary containing a list of the species in solution paired with their activity

        :meta private:
        """
        print("Component Activities:\n")
        print("=====================\n")
        for i in self.components:
            print(i + ":" + "\t {0.magnitude:0.{decimals}f}".format(self.get_activity(i), decimals=decimals))
