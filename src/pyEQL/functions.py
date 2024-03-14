"""
pyEQL functions that take Solution objects as inputs or return Solution objects.

:copyright: 2013-2024 by Ryan S. Kingsbury
:license: LGPL, see LICENSE for more details.

"""

import logging
import math
from typing import Literal

from monty.dev import deprecated

from pyEQL import Solution, ureg

logger = logging.getLogger(__name__)


def gibbs_mix(solution1: Solution, solution2: Solution):
    r"""
    Return the Gibbs energy change associated with mixing two solutions.

    Args:
        solution1: a solution to be mixed.
        solution2: a solution to be mixed.

    Returns:
        The change in Gibbs energy associated with complete mixing of the
        Solutions, in Joules.

    Notes:
        The Gibbs energy of mixing is calculated as follows

        .. math::

            \Delta_{mix} G = \sum_i {(n_c + n_d) R T \ln a_b} - \sum_i {n_c R T \ln a_c} - \sum_i {n_d R T \ln a_d}

        Where :math:`n` is the number of moles of substance, :math:`T` is the temperature in kelvin,
        and  subscripts :math:`b`, :math:`c`, and :math:`d` refer to the concentrated, dilute, and blended
        Solutions, respectively.

        Note that dissociated ions must be counted as separate components,
        so a simple salt dissolved in water is a three component solution (cation,
        anion, and water).

    References:
        Koga, Yoshikata, 2007. Solution Thermodynamics and its Application to Aqueous Solutions:
            A differential approach. Elsevier, 2007, pp. 23-37.

    """
    concentrate = solution1
    dilute = solution2
    blend = solution1 + solution2
    term_list = {concentrate: 0, dilute: 0, blend: 0}

    # calculate the entropy change and number of moles solute for each solution
    for solution in term_list:
        for solute in solution.components:
            if solution.get_amount(solute, "fraction") != 0:
                term_list[solution] += solution.get_amount(solute, "mol") * math.log(solution.get_activity(solute))

    return (ureg.R * blend.temperature.to("K") * (term_list[blend] - term_list[concentrate] - term_list[dilute])).to(
        "J"
    )


def entropy_mix(solution1: Solution, solution2: Solution):
    r"""
    Return the ideal mixing entropy associated with mixing two solutions.

    Parameters:
        solution1, solution2: The two solutions to be mixed.

    Returns:
        The ideal mixing entropy associated with complete mixing of the
        Solutions, in Joules.

    Notes:
        The ideal entropy of mixing is calculated as follows

        .. math::

            \Delta_{mix} S = \sum_i {(n_c + n_d) R T \ln x_b} - \sum_i {n_c R T \ln x_c} - \sum_i {n_d R T \ln x_d}

        Where :math:`n` is the number of moles of substance, :math:`T` is the temperature in kelvin,
        and  subscripts :math:`b`, :math:`c`, and :math:`d` refer to the concentrated, dilute, and blended
        Solutions, respectively.

        Note that dissociated ions must be counted as separate components,
        so a simple salt dissolved in water is a three component solution (cation,
        anion, and water).

    References:
        Koga, Yoshikata, 2007. *Solution Thermodynamics and its Application to Aqueous Solutions:
            A differential approach.* Elsevier, 2007, pp. 23-37.

    """
    concentrate = solution1
    dilute = solution2
    blend = solution1 + solution2
    term_list = {concentrate: 0, dilute: 0, blend: 0}

    # calculate the entropy change and number of moles solute for each solution
    for solution in term_list:
        for solute in solution.components:
            if solution.get_amount(solute, "fraction") != 0:
                term_list[solution] += solution.get_amount(solute, "mol") * math.log(
                    solution.get_amount(solute, "fraction")
                )

    return (ureg.R * blend.temperature.to("K") * (term_list[blend] - term_list[concentrate] - term_list[dilute])).to(
        "J"
    )


def donnan_eql(solution: Solution, fixed_charge: str):
    r"""
    Return a solution object in equilibrium with fixed_charge.

    Args:
        solution : Solution object
            The external solution to be brought into equilibrium with the fixed
            charges
        fixed_charge : str quantity
            String representing the concentration of fixed charges, including sign.
            May be specified in mol/L or mol/kg units. e.g. '1 mol/kg'

    Returns:
        A Solution that has established Donnan equilibrium with the external
        (input) Solution

    Notes:
        The general equation representing the equilibrium between an external
        electrolyte solution and an ion-exchange medium containing fixed charges
        is

        .. math::

            \big(\frac{a_{-}}{\bar a_{-}} \big)^(\frac{1}{z_{-})
            \big(\frac{\bar a_{+}}{a_{+}}\big)^(\frac{1}{z_{+})
            \exp(\frac{\Delta \pi \bar V}{RT z_{+} \nu_{+}})

        Where subscripts :math:`+` and :math:`-` indicate the cation and anion, respectively,
        the overbar indicates the membrane phase,
        :math:`a` represents activity, :math:`z` represents charge, :math:`\nu` represents the stoichiometric
        coefficient, :math:`V` represents the partial molar volume of the salt, and
        :math:`\Delta \pi` is the difference in osmotic pressure between the membrane and the
        solution phase.

        In addition, electroneutrality must prevail within the membrane phase:

        .. math:: \bar C_{+} z_{+} + \bar X + \bar C_{-} z_{-} = 0

        Where :math:`C` represents concentration and :math:`X` is the fixed charge concentration
        in the membrane or ion exchange phase.

        This function solves these two equations simultaneously to arrive at the
        concentrations of the cation and anion in the membrane phase. It returns
        a solution equal to the input solution except that the concentrations of
        the predominant cation and anion have been adjusted according to this
        equilibrium.

        NOTE that this treatment is only capable of equilibrating a single salt.
        This salt is identified by the get_salt() method.

    References:
        Strathmann, Heiner, ed. *Membrane Science and Technology* vol. 9, 2004. Chapter 2, p. 51.
           http://dx.doi.org/10.1016/S0927-5193(04)80033-0

    See Also:
        get_salt()

    """
    # identify the salt
    salt = solution.get_salt()

    # convert fixed_charge in to a quantity
    fixed_charge = ureg.Quantity(fixed_charge)

    # identify variables from the external solution
    conc_cation_soln = solution.get_amount(salt.cation, str(fixed_charge.units))
    conc_anion_soln = solution.get_amount(salt.anion, str(fixed_charge.units))
    act_cation_soln = solution.get_activity(salt.cation)
    act_anion_soln = solution.get_activity(salt.anion)
    z_cation = salt.z_cation
    z_anion = salt.z_anion
    nu_cation = salt.nu_cation

    # get the partial molar volume for the salt, or calculate it from the ions
    # TODO - consider how to incorporate pitzer parameters
    molar_volume = solution.get_property(salt.formula, "size.molar_volume")
    if molar_volume is None:
        cation_vol = solution.get_property(salt.cation, "size.molar_volume")
        anion_vol = solution.get_property(salt.anion, "size.molar_volume")
        if cation_vol is not None and anion_vol is not None:
            molar_volume = cation_vol + anion_vol
        else:
            logger.critical("Required partial molar volume information not available. Aborting.")
            return None

    # initialize the equilibrated solution - start with a direct copy of the
    # input / external solution
    donnan_soln = solution.copy()

    # do nothing if either of the ion concentrations is zero
    if conc_cation_soln.magnitude == 0 or conc_anion_soln.magnitude == 0:
        return donnan_soln

    # define a function representing the donnan equilibrium as a function
    # of the two unknown activities to feed to the nonlinear solver

    # the stuff in the term below doesn't change on iteration, so calculate it up-front
    # assign it the correct units and extract the magnitude for a performance gain
    exp_term = (molar_volume / (ureg.R * solution.temperature * z_cation * nu_cation)).to("1/Pa").magnitude

    def donnan_solve(x):
        """Where x is the magnitude of co-ion concentration."""
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
        donnan_soln.set_amount(salt.cation, str(conc_cation_mem) + units)
        donnan_soln.set_amount(salt.anion, str(conc_anion_mem) + units)

        # get the new concentrations and activities
        act_cation_mem = donnan_soln.get_activity(salt.cation)
        act_anion_mem = donnan_soln.get_activity(salt.anion)

        # compute the difference in osmotic pressure
        # using the magnitudes here helps performance
        delta_pi = donnan_soln.osmotic_pressure.magnitude - solution.osmotic_pressure.magnitude

        return (act_cation_mem / act_cation_soln) ** (1 / z_cation) * (act_anion_soln / act_anion_mem) ** (
            1 / z_anion
        ) - math.exp(delta_pi * exp_term)

    # solve the function above using one of scipy's nonlinear solvers

    from scipy.optimize import brentq

    # determine which ion concentration represents the co-ion
    # call a nonlinear solver to adjust the concentrations per the donnan
    # equilibrium, unless the membrane is uncharged
    # the initial guess is to set the co-ion concentration in the membrane
    # equal to that in the solution
    if fixed_charge.magnitude > 0:
        x = conc_cation_soln.magnitude
        brentq(donnan_solve, 1e-10, x, xtol=0.001)
    elif fixed_charge.magnitude < 0:
        x = conc_anion_soln.magnitude
        brentq(donnan_solve, 1e-10, x, xtol=0.001)
    else:
        pass

    # return the equilibrated solution
    return donnan_soln


@deprecated(
    message="mix() is deprecated and will be removed in the next release! You can now mix solutions using the addition operator, e.g. s_mix = s1 + s2."
)
def mix(s1, s2):  # pragma: no cover
    """
    Mix two solutions together.

    Args:
        s1, s2: The two solutions to be mixed.

    Returns:
        A Solution object that represents the result of mixing s1 and s2.

    Notes:
        The initial volume of the mixed solution is set as the sum of the volumes of s1 and s2. The pressure and
        temperature are volume-weighted averages. The pH and pE values are currently APPROXIMATE because they are
        calculated assuming H+ and e- mix conservatively (i.e., the mixing process does not incorporate any
        equilibration reactions or buffering). Such support is planned in a future release.
    """
    return s1 + s2


@deprecated(
    message="autogenerate() is deprecated and will be removed in the next release! Use Solution.from_preset() instead.)"
)
def autogenerate(
    solution: Literal["seawater", "rainwater", "wastewater", "urine", "Ringers lactate", "normal saline"]
):  # pragma: no cover
    """
    This method provides a quick way to create Solution objects representing
    commonly-encountered solutions, such as seawater, rainwater, and wastewater.

    Parameters
    ----------
    solution : str
                String representing the desired solution
                Valid entries are 'seawater', 'rainwater',
                'wastewater',and 'urine'

    Returns
    -------
    Solution
        A pyEQL Solution object.

    Notes
    -----
    The following sections explain the different solution options:

    - '' - empty solution, equivalent to pyEQL.Solution()
    - 'rainwater' - pure water in equilibrium with atmospheric CO2 at pH 6
    - 'seawater' or 'SW'- Standard Seawater. See Table 4 of the Reference for Composition [1]_
    - 'wastewater' or 'WW' - medium strength domestic wastewater. See Table 3-18 of [2]_
    - 'urine' - typical human urine. See Table 3-15 of [2]_
    - 'normal saline' or 'NS' - normal saline solution used in medicine [3]_
    - 'Ringers lacatate' or 'RL' - Ringer's lactate solution used in medicine [4]_

    References:
    ----------
        .. [1] Millero, Frank J. "The composition of Standard Seawater and the definition of
            the Reference-Composition Salinity Scale." *Deep-sea Research. Part I* 55(1), 2008, 50-72.

        .. [2] Metcalf & Eddy, Inc. et al. *Wastewater Engineering: Treatment and Resource Recovery*, 5th Ed.
                McGraw-Hill, 2013.

        .. [3] https://en.wikipedia.org/wiki/Saline_(medicine)

        .. [4] https://en.wikipedia.org/wiki/Ringer%27s_lactate_solution

    """
    if solution == "":
        temperature = "25 degC"
        pressure = "1 atm"
        pH = 7
        solutes = []
    elif solution == "seawater" or solution == "SW":
        temperature = "25 degC"
        pressure = "1 atm"
        pH = 8.1
        solutes = [
            ["Na+", "10.78145 g/kg"],
            ["Mg+2", "1.28372 g/kg"],
            ["Ca+2", "0.41208 g/kg"],
            ["K+", "0.39910 g/kg"],
            ["Sr+2", "0.00795 g/kg"],
            ["Cl-", "19.35271 g/kg"],
            ["SO4-2", "2.71235 g/kg"],
            ["HCO3-", "0.10481 g/kg"],
            ["Br-", "0.06728 g/kg"],
            ["CO3-2", "0.01434 g/kg"],
            ["B(OH)4", "0.00795 g/kg"],
            ["F-", "0.00130 g/kg"],
            ["OH-", "0.00014 g/kg"],
            ["B(OH)3", "0.01944 g/kg"],
            ["CO2", "0.00042 g/kg"],
        ]
    elif solution == "rainwater":
        temperature = "25 degC"
        pressure = "1 atm"
        pH = 6
        solutes = [["HCO3-", "10^-5.5 mol/L"], ["CO3-2", "10^-9 mol/L"]]
    elif solution == "wastewater" or solution == "WW":
        temperature = "25 degC"
        pressure = "1 atm"
        pH = 7
        solutes = [
            ["NH3", "24.3 mg/L"],
            ["PO4-3", "7.6 mg/L"],
            ["C6H12O6", "410 mg/L"],
            ["K+", "16 mg/L"],
            ["Cl-", "59 mg/L"],
            ["SO4-2", "26 mg/L"],
        ]
        logger.warning("Total organic carbon in wastewater is approximated as glucose")
    elif solution == "urine":
        temperature = "25 degC"
        pressure = "1 atm"
        pH = 7
        solutes = [
            ["CON2H4", "20,000 mg/L"],
            ["C4H7N3O", "1,000 mg/L"],
            ["C5H4N4O3", "300 mg/L"],
            ["NH4+", "500 mg/L"],
            ["HCO3-", "300 mg/L"],
            ["NH4+", "500 mg/L"],
            ["Mg+2", "100 mg/L"],
            ["PO4-3", "1200 mg/L"],
            ["Na+", "6000 mg/L"],
            ["K+", "1500 mg/L"],
            ["Cl-", "1900 mg/L"],
            ["SO4-2", "1800 mg/L"],
        ]
    elif solution == "normal saline" or solution == "NS":
        temperature = "25 degC"
        pressure = "1 atm"
        pH = 7
        solutes = [["Na+", "154 mmol/L"], ["Cl-", "154 mmol/L"]]
    elif solution == "Ringers lactate" or solution == "RL":
        temperature = "25 degC"
        pressure = "1 atm"
        pH = 6.5
        solutes = [
            ["Na+", "130 mmol/L"],
            ["Cl-", "109 mmol/L"],
            ["K+", "4 mmol/L"],
            ["Ca+2", "1.5 mmol/L"],
            ["C3H5O3-", "28 mmol/L"],
        ]
    else:
        logger.error("Invalid solution entered - %s" % solution)
        return None

    return Solution(solutes, temperature=temperature, pressure=pressure, pH=pH)
