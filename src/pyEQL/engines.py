"""
pyEQL engines for computing aqueous equilibria (e.g., speciation, redox, etc.).

:copyright: 2013-2024 by Ryan S. Kingsbury
:license: LGPL, see LICENSE for more details.

"""

import logging
import os
import warnings
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Literal

from phreeqpython import PhreeqPython

import pyEQL.activity_correction as ac
from pyEQL import ureg
from pyEQL.salt_ion_match import Salt
from pyEQL.utils import standardize_formula

# These are the only elements that are allowed to have parenthetical oxidation states
# PHREEQC will ignore others (e.g., 'Na(1)')
SPECIAL_ELEMENTS = ["S", "C", "N", "Cu", "Fe", "Mn"]

logger = logging.getLogger(__name__)


class EOS(ABC):
    """
    Abstract base class for pyEQL equation of state classes.

    The intent is that concrete implementations of this class make use of the
    standalone functions available in pyEQL.activity_correction and pyEQL.equilibrium
    as much as possible. This facilitates robust unit testing while allowing users
    to "mix and match" or customize the various models as needed.
    """

    @abstractmethod
    def get_activity_coefficient(self, solution, solute):
        """
        Return the *molal scale* activity coefficient of solute, given a Solution
        object.

        Args:
            solution: pyEQL Solution object
            solute: str identifying the solute of interest

        Returns
            Quantity: dimensionless quantity object

        Raises:
            ValueError if the calculation cannot be completed, e.g. due to insufficient number of parameters.
        """

    @abstractmethod
    def get_osmotic_coefficient(self, solution):
        """
        Return the *molal scale* osmotic coefficient of a Solution.

        Args:
            solution: pyEQL Solution object

        Returns
            Quantity: dimensionless molal scale osmotic coefficient

        Raises:
            ValueError if the calculation cannot be completed, e.g. due to insufficient number of parameters.
        """

    @abstractmethod
    def get_solute_volume(self):
        """
        Return the volume of only the solutes.

        Args:
            solution: pyEQL Solution object

        Returns
            Quantity: solute volume in L

        Raises:
            ValueError if the calculation cannot be completed, e.g. due to insufficient number of parameters.
        """

    @abstractmethod
    def equilibrate(self, solution):
        """
        Adjust the speciation and pH of a Solution object to achieve chemical equilibrium.

        The Solution should be modified in-place, likely using add_moles / set_moles, etc.

        Args:
            solution: pyEQL Solution object

        Returns
            Nothing. The speciation of the Solution is modified in-place.

        Raises:
            ValueError if the calculation cannot be completed, e.g. due to insufficient number of parameters or lack of convergence.
        """


class IdealEOS(EOS):
    """Ideal solution equation of state engine."""

    def get_activity_coefficient(self, solution, solute):
        """
        Return the *molal scale* activity coefficient of solute, given a Solution
        object.
        """
        return ureg.Quantity(1, "dimensionless")

    def get_osmotic_coefficient(self, solution):
        """
        Return the *molal scale* osmotic coefficient of solute, given a Solution
        object.
        """
        return ureg.Quantity(1, "dimensionless")

    def get_solute_volume(self, solution):
        """Return the volume of the solutes."""
        return ureg.Quantity(0, "L")

    def equilibrate(self, solution):
        """Adjust the speciation of a Solution object to achieve chemical equilibrium."""
        warnings.warn("equilibrate() has no effect in IdealEOS!")
        return


class NativeEOS(EOS):
    """
    pyEQL's native EOS. Uses the Pitzer model when possible, falls
    back to other models (e.g. Debye-Huckel) based on ionic strength
    if sufficient parameters are not available.
    """

    def __init__(
        self,
        phreeqc_db: Literal["vitens.dat", "wateq4f_PWN.dat", "pitzer.dat", "llnl.dat", "geothermal.dat"] = "llnl.dat",
    ):
        """
        Args:
            phreeqc_db: Name of the PHREEQC database file to use for solution thermodynamics
                and speciation calculations. Generally speaking, `llnl.dat` is recommended
                for moderate salinity water and prediction of mineral solubilities,
                `wateq4f_PWN.dat` is recommended for low to moderate salinity waters. It is
                similar to vitens.dat but has many more species. `pitzer.dat` is recommended
                when accurate activity coefficients in solutions above 1 M TDS are desired, but
                it has fewer species than the other databases. `llnl.dat` and `geothermal.dat`
                may offer improved prediction of LSI but currently these databases are not
                usable because they do not allow for conductivity calculations.
        """
        self.phreeqc_db = phreeqc_db
        # database files in this list are not distributed with phreeqpython
        self.db_path = (
            Path(os.path.dirname(__file__)) / "database" if self.phreeqc_db in ["llnl.dat", "geothermal.dat"] else None
        )
        # create the PhreeqcPython instance
        # try/except added to catch unsupported architectures, such as Apple Silicon
        try:
            self.pp = PhreeqPython(database=self.phreeqc_db, database_directory=self.db_path)
        except OSError:
            logger.error(
                "OSError encountered when trying to instantiate phreeqpython. Most likely this means you"
                " are running on an architecture that is not supported by PHREEQC, such as Apple M1/M2 chips."
                " pyEQL will work, but equilibrate() will have no effect."
            )
        # attributes to hold the PhreeqPython solution.
        self.ppsol = None
        # store the solution composition to see whether we need to re-instantiate the solution
        self._stored_comp = None

    def _setup_ppsol(self, solution):
        """
        Helper method to set up a PhreeqPython solution for subsequent analysis.
        """
        self._stored_comp = solution.components.copy()
        solv_mass = solution.solvent_mass.to("kg").magnitude
        # inherit bulk solution properties
        d = {
            "temp": solution.temperature.to("degC").magnitude,
            "units": "mol/kgw",  # to avoid confusion about volume, use mol/kgw which seems more robust in PHREEQC
            "pH": solution.pH,
            "pe": solution.pE,
            "redox": "pe",  # hard-coded to use the pe
            # PHREEQC will assume 1 kg if not specified, there is also no direct way to specify volume, so we
            # really have to specify the solvent mass in 1 liter of solution
            "water": solv_mass,
        }
        if solution.balance_charge == "pH":
            d["pH"] = str(d["pH"]) + " charge"
        if solution.balance_charge == "pE":
            d["pe"] = str(d["pe"]) + " charge"

        # add the composition to the dict
        # also, skip H and O
        for el, mol in solution.get_el_amt_dict().items():
            # strip off the oxi state
            bare_el = el.split("(")[0]
            if bare_el in SPECIAL_ELEMENTS:
                # PHREEQC will ignore float-formatted oxi states. Need to make sure we are
                # passing, e.g. 'C(4)' and not 'C(4.0)'
                key = f'{bare_el}({int(float(el.split("(")[-1].split(")")[0]))})'
            elif bare_el in ["H", "O"]:
                continue
            else:
                key = bare_el

            d[key] = str(mol / solv_mass)

            # tell PHREEQC which species to use for charge balance
            if (
                solution.balance_charge is not None
                and solution.balance_charge in solution.get_components_by_element()[el]
            ):
                d[key] += " charge"

        # create the PHREEQC solution object
        try:
            ppsol = self.pp.add_solution(d)
        except Exception as e:
            print(d)
            # catch problems with the input to phreeqc
            raise ValueError(
                "There is a problem with your input. The error message received from "
                f" phreeqpython is:\n\n {e}\n Check your input arguments, especially "
                "the composition dictionary, and try again."
            )

        self.ppsol = ppsol

    def _destroy_ppsol(self):
        """Remove the PhreeqPython solution from memory"""
        if self.ppsol is not None:
            self.ppsol.forget()
            self.ppsol = None

    def get_activity_coefficient(self, solution, solute):
        r"""
        Whenever the appropriate parameters are available, the Pitzer model [may]_ is used.
        If no Pitzer parameters are available, then the appropriate equations are selected
        according to the following logic: [stumm]_.

        I <= 0.0005: Debye-Huckel equation
        0.005 < I <= 0.1:  Guntelberg approximation
        0.1 < I <= 0.5: Davies equation
        I > 0.5: Raises a warning and returns activity coefficient = 1

        The ionic strength, activity coefficients, and activities are all
        calculated based on the molal (mol/kg) concentration scale. If a different
        scale is given as input, then the molal-scale activity coefficient :math:`\gamma_\pm` is
        converted according to [rbs]_

        .. math:: f_\pm = \gamma_\pm * (1 + M_w \sum_i \nu_i m_i)

        .. math:: y_\pm = \frac{m \rho_w}{C \gamma_\pm}

        where :math:`f_\pm` is the rational activity coefficient, :math:`M_w` is
        the molecular weight of water, the summation represents the total molality of
        all solute  species, :math:`y_\pm` is the molar activity coefficient,
        :math:`\rho_w` is the density of pure water, :math:`m` and :math:`C` are
        the molal and molar concentrations of the chosen salt (not individual solute), respectively.

        Args:
            solute: String representing the name of the solute of interest
            scale: The concentration scale for the returned activity coefficient.
                Valid options are "molal", "molar", and "rational" (i.e., mole fraction).
                By default, the molal scale activity coefficient is returned.

        Returns:
            The mean ion activity coefficient of the solute in question on  the selected scale.


        Notes:
            For multicomponent mixtures, pyEQL implements the "effective Pitzer model"
            presented by Mistry et al. [mistry]_. In this model, the activity coefficient
            of a salt in a multicomponent mixture is calculated using an "effective
            molality," which is the molality that would result in a single-salt
            mixture with the same total ionic strength as the multicomponent solution.

            .. math:: m_{effective} = \frac{2 I}{(\nu_{+} z_{+}^2 + \nu_{-}- z_{-}^2)}

        References:
            .. [may] May, P. M., Rowland, D., Hefter, G., & Königsberger, E. (2011).
                A Generic and Updatable Pitzer Characterization of Aqueous Binary Electrolyte Solutions at 1 bar and 25 °C.
               *Journal of Chemical & Engineering Data*, 56(12), 5066-5077. doi:10.1021/je2009329

        .. [stumm] Stumm, Werner and Morgan, James J. *Aquatic Chemistry*, 3rd ed,
               pp 165. Wiley Interscience, 1996.

        .. [rbs] Robinson, R. A.; Stokes, R. H. Electrolyte Solutions: Second Revised
               Edition; Butterworths: London, 1968, p.32.

        .. [mistry] Mistry, K. H.; Hunter, H. a.; Lienhard V, J. H. Effect of composition and nonideal solution behavior on
               desalination calculations for mixed electrolyte solutions with comparison to seawater. Desalination 2013, 318, 34-47.

        See Also:
            :attr:`pyEQL.solution.Solution.ionic_strength`
            :func:`pyEQL.activity_correction.get_activity_coefficient_debyehuckel`
            :func:`pyEQL.activity_correction.get_activity_coefficient_guntelberg`
            :func:`pyEQL.activity_correction.get_activity_coefficient_davies`
            :func:`pyEQL.activity_correction.get_activity_coefficient_pitzer`
        """
        # identify the predominant salt that this ion is a member of
        salt = None
        rform = standardize_formula(solute)
        for v in solution.get_salt_dict().values():
            if v == "HOH":
                continue
            if rform == v["cation"] or rform == v["anion"]:
                del v["mol"]
                salt = Salt.from_dict(v)
                break

        # show an error if no salt can be found that contains the solute
        if salt is None:
            logger.error("No salts found that contain solute %s. Returning unit activity coefficient." % solute)
            return ureg.Quantity(1, "dimensionless")

        # use the Pitzer model for higher ionic strength, if the parameters are available
        # search for Pitzer parameters
        param = solution.get_property(salt.formula, "model_parameters.activity_pitzer")
        if param is not None:
            # TODO - consider re-enabling a log message recording what salt(s) are used as basis for activity calculation
            logger.info(f"Calculating activity coefficient based on parent salt {salt.formula}")

            # determine alpha1 and alpha2 based on the type of salt
            # see the May reference for the rules used to determine
            # alpha1 and alpha2 based on charge
            if salt.nu_cation >= 2 and salt.nu_anion <= -2:
                if salt.nu_cation >= 3 or salt.nu_anion <= -3:
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
            # molality = (solution.get_amount(salt.cation,'mol/kg')/salt.nu_cation+solution.get_amount(salt.anion,'mol/kg')/salt.nu_anion)/2

            # determine the effective molality of the salt in the solution
            molality = salt.get_effective_molality(solution.ionic_strength)

            activity_coefficient = ac.get_activity_coefficient_pitzer(
                solution.ionic_strength,
                molality,
                alpha1,
                alpha2,
                ureg.Quantity(param["Beta0"]["value"]).magnitude,
                ureg.Quantity(param["Beta1"]["value"]).magnitude,
                ureg.Quantity(param["Beta2"]["value"]).magnitude,
                ureg.Quantity(param["Cphi"]["value"]).magnitude,
                salt.z_cation,
                salt.z_anion,
                salt.nu_cation,
                salt.nu_anion,
                str(solution.temperature),
            )

            logger.debug(
                f"Calculated activity coefficient of species {solute} as {activity_coefficient} based on salt"
                f" {salt} using Pitzer model"
            )
            molal = activity_coefficient

        # for very low ionic strength, use the Debye-Huckel limiting law
        elif solution.ionic_strength.magnitude <= 0.005:
            logger.debug(
                f"Ionic strength = {solution.ionic_strength}. Using Debye-Huckel to calculate activity coefficient."
            )
            molal = ac.get_activity_coefficient_debyehuckel(
                solution.ionic_strength,
                solution.get_property(solute, "charge"),
                str(solution.temperature),
            )

        # use the Guntelberg approximation for 0.005 < I < 0.1
        elif solution.ionic_strength.magnitude <= 0.1:
            logger.debug(
                f"Ionic strength = {solution.ionic_strength}. Using Guntelberg to calculate activity coefficient."
            )
            molal = ac.get_activity_coefficient_guntelberg(
                solution.ionic_strength,
                solution.get_property(solute, "charge"),
                str(solution.temperature),
            )

        # use the Davies equation for 0.1 < I < 0.5
        elif solution.ionic_strength.magnitude <= 0.5:
            logger.debug(
                f"Ionic strength = {solution.ionic_strength}. Using Davies equation to calculate activity coefficient."
            )
            molal = ac.get_activity_coefficient_davies(
                solution.ionic_strength,
                solution.get_property(solute, "charge"),
                str(solution.temperature),
            )

        else:
            logger.error(
                f"Ionic strength too high to estimate activity for species {solute}. Specify parameters for Pitzer "
                "model. Returning unit activity coefficient"
            )

            molal = ureg.Quantity(1, "dimensionless")

        return molal

    def get_osmotic_coefficient(self, solution):
        r"""
        Return the *molal scale* osmotic coefficient of solute, given a Solution
        object.

        Osmotic coefficient is calculated using the Pitzer model. [may]_ If appropriate parameters for
        the model are not available, then pyEQL raises a WARNING and returns an osmotic
        coefficient of 1.

        If the 'rational' scale is given as input, then the molal-scale osmotic
        coefficient :math:`\phi` is converted according to [rbs]_

        .. math:: g = - \phi M_{w} \frac{\sum_{i} \nu_{i} m_{i}}{\ln x_{w}}

        where :math:`g` is the rational osmotic coefficient, :math:`M_{w}` is
        the molecular weight of water, the summation represents the total molality of
        all solute  species, and :math:`x_{w}` is the mole fraction of water.

        Args:
            scale: The concentration scale for the returned osmotic coefficient. Valid options are "molal",
                "rational" (i.e., mole fraction), and "fugacity".  By default, the molal scale osmotic
                coefficient is returned.

        Returns:
            Quantity:
                The osmotic coefficient

        See Also:
            :meth:`pyEQL.solution.Solution.get_water_activity`
            :meth:`pyEQL.solution.Solution.get_salt`
            :attr:`pyEQL.solution.Solution.ionic_strength`

        Notes:
            For multicomponent mixtures, pyEQL adopts the "effective Pitzer model"
            presented by Mistry et al. [mstry]_. In this approach, the osmotic coefficient of
            each individual salt is calculated using the normal Pitzer model based
            on its respective concentration. Then, an effective osmotic coefficient
            is calculated as the concentration-weighted average of the individual
            osmotic coefficients.

            For example, in a mixture of 0.5 M NaCl and 0.5 M KBr, one would calculate
            the osmotic coefficient for each salt using a concentration of 0.5 M and
            an ionic strength of 1 M. Then, one would average the two resulting
            osmotic coefficients to obtain an effective osmotic coefficient for the
            mixture.

            (Note: in the paper referenced below, the effective osmotic coefficient is determined by weighting
            using the "effective molality" rather than the true molality. Subsequent checking and correspondence with
            the author confirmed that the weight factor should be the true molality, and that is what is implemented
            in pyEQL.)

        References:
            .. [may] May, P. M., Rowland, D., Hefter, G., & Königsberger, E. (2011).
                A Generic and Updatable Pitzer Characterization of Aqueous Binary Electrolyte Solutions at 1 bar and
                25 °C. Journal of Chemical & Engineering Data, 56(12), 5066-5077. doi:10.1021/je2009329

            .. [rbs] Robinson, R. A.; Stokes, R. H. Electrolyte Solutions: Second Revised
                Edition; Butterworths: London, 1968, p.32.

            .. [mstry] Mistry, K. H.; Hunter, H. a.; Lienhard V, J. H. Effect of composition and nonideal solution
               behavior on desalination calculations for mixed electrolyte solutions with comparison to seawater. Desalination 2013, 318, 34-47.

        Examples:
            >>> s1 = pyEQL.Solution({'Na+': '0.2 mol/kg', 'Cl-': '0.2 mol/kg'})
            >>> s1.get_osmotic_coefficient()
            <Quantity(0.923715281, 'dimensionless')>

            >>> s1 = pyEQL.Solution({'Mg+2': '0.3 mol/kg', 'Cl-': '0.6 mol/kg'},temperature='30 degC')
            >>> s1.get_osmotic_coefficient()
            <Quantity(0.891409618, 'dimensionless')>

        """
        ionic_strength = solution.ionic_strength

        effective_osmotic_sum = 0
        molality_sum = 0

        # loop through all the salts in the solution, calculate the osmotic
        # coefficint for each, and average them into an effective osmotic
        # coefficient
        for d in solution.get_salt_dict().values():
            item = Salt(d["cation"], d["anion"])
            # ignore HOH in the salt list
            if item.formula == "HOH":
                continue

            # determine alpha1 and alpha2 based on the type of salt
            # see the May reference for the rules used to determine
            # alpha1 and alpha2 based on charge
            if item.z_cation >= 2 and item.z_anion <= -2:
                if item.z_cation >= 3 or item.z_anion <= -3:
                    alpha1 = 2.0
                    alpha2 = 50.0
                else:
                    alpha1 = 1.4
                    alpha2 = 12.0
            else:
                alpha1 = 2.0
                alpha2 = 0

            # set the concentration as the average concentration of the cation and
            # anion in the salt, accounting for stoichiometry
            # concentration = (solution.get_amount(Salt.cation,'mol/kg')/Salt.nu_cation + \
            # solution.get_amount(Salt.anion,'mol/kg')/Salt.nu_anion)/2

            # get the effective molality of the salt
            concentration = ureg.Quantity(d["mol"], "mol") / solution.solvent_mass

            molality_sum += concentration

            param = solution.get_property(item.formula, "model_parameters.activity_pitzer")
            if param is not None:
                osmotic_coefficient = ac.get_osmotic_coefficient_pitzer(
                    ionic_strength,
                    concentration,
                    alpha1,
                    alpha2,
                    ureg.Quantity(param["Beta0"]["value"]).magnitude,
                    ureg.Quantity(param["Beta1"]["value"]).magnitude,
                    ureg.Quantity(param["Beta2"]["value"]).magnitude,
                    ureg.Quantity(param["Cphi"]["value"]).magnitude,
                    item.z_cation,
                    item.z_anion,
                    item.nu_cation,
                    item.nu_anion,
                    str(solution.temperature),
                )

                logger.debug(
                    f"Calculated osmotic coefficient of water as {osmotic_coefficient} based on salt "
                    f"{item.formula} using Pitzer model"
                )
                effective_osmotic_sum += concentration * osmotic_coefficient

            else:
                logger.debug(
                    f"Returning unit osmotic coefficient for salt {item.formula} because Pitzer parameters are not"
                    "available in database."
                )
                effective_osmotic_sum += concentration * 1

        try:
            return effective_osmotic_sum / molality_sum
        except ZeroDivisionError:
            # this means the solution is empty
            return 1

    def get_solute_volume(self, solution):
        """Return the volume of the solutes."""
        # identify the predominant salt in the solution
        salt = solution.get_salt()
        solute_vol = ureg.Quantity(0, "L")

        # use the pitzer approach if parameters are available
        pitzer_calc = False

        param = solution.get_property(salt.formula, "model_parameters.molar_volume_pitzer")
        if param is not None:
            # determine the average molality of the salt
            # this is necessary for solutions inside e.g. an ion exchange
            # membrane, where the cation and anion concentrations may be
            # unequal
            molality = (solution.get_amount(salt.cation, "mol/kg") + solution.get_amount(salt.anion, "mol/kg")) / 2

            # determine alpha1 and alpha2 based on the type of salt
            # see the May reference for the rules used to determine
            # alpha1 and alpha2 based on charge
            if salt.nu_cation >= 2 and salt.nu_anion >= 2:
                if salt.nu_cation >= 3 or salt.nu_anion >= 3:
                    alpha1 = 2
                    alpha2 = 50
                else:
                    alpha1 = 1.4
                    alpha2 = 12
            else:
                alpha1 = 2
                alpha2 = 0

            apparent_vol = ac.get_apparent_volume_pitzer(
                solution.ionic_strength,
                molality,
                alpha1,
                alpha2,
                ureg.Quantity(param["Beta0"]["value"]).magnitude,
                ureg.Quantity(param["Beta1"]["value"]).magnitude,
                ureg.Quantity(param["Beta2"]["value"]).magnitude,
                ureg.Quantity(param["Cphi"]["value"]).magnitude,
                ureg.Quantity(param["V_o"]["value"]).magnitude,
                salt.z_cation,
                salt.z_anion,
                salt.nu_cation,
                salt.nu_anion,
                str(solution.temperature),
            )

            solute_vol += (
                apparent_vol
                * (
                    solution.get_amount(salt.cation, "mol") / salt.nu_cation
                    + solution.get_amount(salt.anion, "mol") / salt.nu_anion
                )
                / 2
            )

            pitzer_calc = True

            logger.debug("Updated solution volume using Pitzer model for solute %s" % salt.formula)

        # add the partial molar volume of any other solutes, except for water
        # or the parent salt, which is already accounted for by the Pitzer parameters
        for solute, mol in solution.components.items():
            # ignore water
            if solute in ["H2O", "HOH", "H2O(aq)"]:
                continue

            # ignore the salt cation and anion, if already accounted for by Pitzer
            if pitzer_calc is True and solute in [salt.anion, salt.cation]:
                continue

            part_vol = solution.get_property(solute, "size.molar_volume")
            if part_vol is not None:
                solute_vol += part_vol * ureg.Quantity(mol, "mol")
                logger.debug("Updated solution volume using direct partial molar volume for solute %s" % solute)

            else:
                logger.warning(
                    f"Volume of solute {solute} will be ignored because partial molar volume data are not available."
                )

        return solute_vol.to("L")

    def equilibrate(self, solution):
        """Adjust the speciation of a Solution object to achieve chemical equilibrium."""
        if self.ppsol is not None:
            self.ppsol.forget()
        self._setup_ppsol(solution)

        # store the original solvent mass
        orig_solvent_moles = solution.components[solution.solvent]

        # use the output from PHREEQC to update the Solution composition
        # the .species_moles attribute should return MOLES (not moles per ___)
        for s, mol in self.ppsol.species_moles.items():
            solution.components[s] = mol

        # make sure all species are accounted for
        charge_adjust = 0
        assert set(self._stored_comp.keys()) - set(solution.components.keys()) == set()

        # log a message if any components were not touched by PHREEQC
        # if that was the case, re-adjust the charge balance to account for those species (since PHREEQC did not)
        missing_species = set(self._stored_comp.keys()) - {standardize_formula(s) for s in self.ppsol.species}
        if len(missing_species) > 0:
            logger.warning(
                f"After equilibration, the amounts of species {missing_species} were not modified "
                "by PHREEQC. These species are likely absent from its database."
            )
        for s in missing_species:
            charge_adjust += -1 * solution.get_amount(s, "eq").magnitude
        if charge_adjust != 0:
            logger.warning(
                "After equilibration, the charge balance of the solution was not electroneutral."
                f" {charge_adjust} eq of charge were added via {solution.balance_charge}"
            )

            # re-adjust charge balance
            if solution.balance_charge is None:
                pass
            elif solution.balance_charge == "pH":
                solution.components["H+"] += charge_adjust
            elif solution.balance_charge == "pE":
                raise NotImplementedError
            else:
                z = solution.get_property(solution.balance_charge, "charge")
                solution.add_amount(solution.balance_charge, f"{charge_adjust/z} mol")

        # rescale the solvent mass to ensure the total mass of solution does not change
        # this is important because PHREEQC and the pyEQL database may use slightly different molecular
        # weights for water. Since water amount is passed to PHREEQC in kg but returned in moles, each
        # call to equilibrate can thus result in a slight change in the Solution mass.
        solution.components[solution.solvent] = orig_solvent_moles

    def __deepcopy__(self, memo):
        # custom deepcopy required because the PhreeqPython instance used by the Native and Phreeqc engines
        # is not pickle-able.
        import copy

        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            if k == "pp":
                result.pp = PhreeqPython(database=self.phreeqc_db, database_directory=self.db_path)
                continue
            setattr(result, k, copy.deepcopy(v, memo))
        return result


class PhreeqcEOS(NativeEOS):
    """Engine based on the PhreeqC model, as implemented via the phreeqpython package."""

    def __init__(
        self,
        phreeqc_db: Literal[
            "vitens.dat", "wateq4f_PWN.dat", "pitzer.dat", "llnl.dat", "geothermal.dat"
        ] = "phreeqc.dat",
    ):
        """
        Args:
        phreeqc_db: Name of the PHREEQC database file to use for solution thermodynamics
                and speciation calculations. Generally speaking, `llnl.dat` is recommended
                for moderate salinity water and prediction of mineral solubilities,
                `wateq4f_PWN.dat` is recommended for low to moderate salinity waters. It is
                similar to vitens.dat but has many more species. `pitzer.dat` is recommended
                when accurate activity coefficients in solutions above 1 M TDS are desired, but
                it has fewer species than the other databases. `llnl.dat` and `geothermal.dat`
                may offer improved prediction of LSI but currently these databases are not
                usable because they do not allow for conductivity calculations.
        """
        super().__init__(phreeqc_db=phreeqc_db)

    def get_activity_coefficient(self, solution, solute):
        """
        Return the *molal scale* activity coefficient of solute, given a Solution
        object.
        """
        if self.ppsol is None or solution.components != self._stored_comp:
            self._destroy_ppsol()
            self._setup_ppsol(solution)

        # translate the species into keys that phreeqc will understand
        k = standardize_formula(solute)
        spl = k.split("[")
        el = spl[0]
        chg = spl[1].split("]")[0]
        if chg[-1] == "1":
            chg = chg[0]  # just pass + or -, not +1 / -1
        k = el + chg

        # calculate the molal scale activity coefficient
        # act = self.ppsol.activity(k, "mol") / self.ppsol.molality(k, "mol")
        act = self.ppsol.pp.ip.get_activity(self.ppsol.number, k) / self.ppsol.pp.ip.get_molality(self.ppsol.number, k)

        return ureg.Quantity(act, "dimensionless")

    def get_osmotic_coefficient(self, solution):
        """
        Return the *molal scale* osmotic coefficient of solute, given a Solution
        object.

        PHREEQC appears to assume a unit osmotic coefficient unless the pitzer database
        is used. Unfortunately, there is no easy way to access the osmotic coefficient
        via phreeqcpython
        """
        # TODO - find a way to access or calculate osmotic coefficient
        return ureg.Quantity(1, "dimensionless")

    def get_solute_volume(self, solution):
        """Return the volume of the solutes."""
        # TODO - phreeqc seems to have no concept of volume, but it does calculate density
        return ureg.Quantity(0, "L")
