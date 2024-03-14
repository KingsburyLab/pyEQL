"""
pyEQL methods for chemical equilibrium calculations (e.g. acid/base, reactions,
redox, complexation, etc.).

NOTE: these methods are not currently used but are here for the future.

:copyright: 2013-2024 by Ryan S. Kingsbury
:license: LGPL, see LICENSE for more details.

"""

# import libraries for scientific functions
import logging
import math

from pyEQL import ureg

logger = logging.getLogger(__name__)

# TODO - not used. Remove?
SPECIES_ALIAISES = {
    "Sodium": "Na+",
    "Potassium": "K+",
    "Calcium": "Ca+2",
    "Barium": "Ba+2",
    "Strontium": "Sr+2",
    "Magnesium": "Mg+2",
    "Chloride": "Cl-",
    "Fluoride": "F-",
    "Nitrate": "NO3-",
    "Ammonium": "NH4+",
    "Sulfate": "SO4-2",
    "Phosphate": "PO4-3",
    "Carbonate": "CO3-2",
    "Bicarbonate": "HCO3-",
}


def adjust_temp_pitzer(c1, c2, c3, c4, c5, temp, temp_ref=ureg.Quantity("298.15 K")):
    """
    Calculate a parameter for the Pitzer model based on temperature-dependent coefficients c1,c2,c3,c4,and c5.

    Args:
        c1, c2, c3, c4, c5 (float): Temperature-dependent coefficients for the pitzer parameter of interest.
        temp (Quantity): The temperature at which the Pitzer parameter is to be calculated.
        temp_ref (Quantity, optional): The reference temperature on which the parameters are based. Defaults to 298.15 K if omitted.

    Note:
        As described in the PHREEQC documentation.
    """
    return (
        c1
        + c2 * (1 / temp + 1 / temp_ref)
        + c2 * math.log(temp / temp_ref)
        + c3 * (temp - temp_ref)
        + c4 * (temp**2 - temp_ref**2)
        + c5 * (temp**-2 - temp_ref**-2)
    )


def adjust_temp_vanthoff(equilibrium_constant, enthalpy, temperature, reference_temperature=ureg.Quantity(25, "degC")):
    r"""
    Adjust a reaction equilibrium constant from one temperature to another.

    Args:
        equilibrium_constant (float): The reaction equilibrium constant for the reaction.
        enthalpy (Quantity): The enthalpy change (delta H) for the reaction in kJ/mol. Assumed independent of temperature (see Notes).
        temperature (Quantity): The desired reaction temperature in degrees Celsius.
        reference_temperature (Quantity, optional): The temperature at which equilibrium_constant is valid. Defaults to 25 degrees C if omitted.

    Returns:
        float: The adjusted reaction equilibrium constant.

    Note:
        This function implements the Van't Hoff equation to adjust measured equilibrium constants to other temperatures.

        .. math::
            ln(K2 / K1) = {\delta H \over R} ( {1 \over T_1} - {1 \over T_2} )

        This implementation assumes that the enthalpy is independent of temperature over the range of interest.

    References:
        Stumm, Werner and Morgan, James J. Aquatic Chemistry, 3rd ed, pp 53. Wiley Interscience, 1996.

    Examples:
        >>> adjust_temp_vanthoff(0.15,ureg.Quantity('-197.6 kJ/mol'),ureg.Quantity('42 degC'),ureg.Quantity(' 25degC')) #doctest: +ELLIPSIS
        0.00203566...

        If the 'ref_temperature' parameter is omitted, a default of 25 C is used.

        >>> adjust_temp_vanthoff(0.15,ureg.Quantity('-197.6 kJ/mol'),ureg.Quantity('42 degC')) #doctest: +ELLIPSIS
        0.00203566...
    """
    output = equilibrium_constant * math.exp(
        enthalpy / ureg.R * (1 / reference_temperature.to("K") - 1 / temperature.to("K"))
    )

    logger.debug(
        "Adjusted equilibrium constant K=%s from %s to %s degrees Celsius with Delta H = %s. Adjusted K = %s % equilibrium_constant,reference_temperature,temperature,enthalpy,output"
    )

    logger.info(
        "Note that the Van't Hoff equation assumes enthalpy is independent of temperature over the range of interest"
    )
    return output


def adjust_temp_arrhenius(
    rate_constant,
    activation_energy,
    temperature,
    reference_temperature=ureg.Quantity(25, "degC"),
):
    r"""
    Adjust a reaction equilibrium constant from one temperature to another.

    Args:
        rate_constant (Quantity): The parameter value (usually a rate constant) being adjusted.
        activation_energy (Quantity): The activation energy of the process, in kJ/mol.
        temperature (Quantity): The desired reaction temperature.
        reference_temperature (Quantity, optional): The temperature at which equilibrium_constant is valid. Defaults to 25 degrees C if omitted.

    Returns:
        Quantity: The adjusted reaction equilibrium constant.

    Notes:
        This function implements the Arrhenius equation to adjust measured rate constants to other temperatures.
        TODO - add better reference

        .. math::
            ln(\frac{K2}{K1} = \frac{E_a}{R} ( \frac{1}{T_{1}} - {\frac{1}{T_2}} )

    References:
        http://chemwiki.ucdavis.edu/Physical_Chemistry/Kinetics/Reaction_Rates/Temperature_Dependence_of_Reaction_Rates/Arrhenius_Equation

    Examples:
        >>> adjust_temp_arrhenius(7,900*ureg.Quantity('kJ/mol'),37*ureg.Quantity('degC'),97*ureg.Quantity('degC')) #doctest: +ELLIPSIS
        1.8867225...e-24
    """
    output = rate_constant * math.exp(
        activation_energy / ureg.R * (1 / reference_temperature.to("K") - 1 / temperature.to("K"))
    )

    logger.debug(
        f"Adjusted parameter {rate_constant} from {reference_temperature} to {temperature} degrees Celsius with"
        f"Activation Energy = {activation_energy}s kJ/mol. Adjusted value = {output}"
    )

    return output


def alpha(n, pH, pKa_list):
    r"""
    Returns the acid-base distribution coefficient (alpha) of an acid in the n-deprotonated form at a given pH.

    Args:
        n (int): The number of protons that have been lost by the desired form of the acid. Also the subscript on the
            alpha value. E.g. for bicarbonate (HCO3-), n=1 because 1 proton has been lost from the fully-protonated
            carbonic acid (H2CO3) form.
        pH (float or int): The pH of the solution.
        pKa_list (list of floats or ints): The pKa values (negative log of equilibrium constants) for the acid of
            interest. There must be a minimum of n pKa values in the list.

    Returns:
        float: The fraction of total acid present in the specified form.

    Notes:
        The acid-base coefficient is calculated as follows: [stm]_

        .. math::

            \alpha_n = \frac{term_n}{[H+]^n + k_{a1}[H+]^{n-1} + k_{a1}k_{a2}[H+]^{n-2} ... k_{a1}k_{a2}...k_{an} }

        Where :math: '\term_n' refers to the nth term in the denominator, starting from 0

    References:
        Stumm, Werner and Morgan, James J. Aquatic Chemistry, 3rd ed, pp 127-130. Wiley Interscience, 1996.

    Examples:
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
    """
    # generate an error if no pKa values are specified
    if len(pKa_list) == 0:
        raise ValueError("No pKa values given. Cannot calculate distribution coefficient.")

    # generate an error if n > number of pKa values
    if len(pKa_list) < n:
        raise ValueError("Insufficient number of pKa values given. Cannot calculate distribution coefficient.")

    # sort the pKa list
    pKa_list = sorted(pKa_list)

    # determine how many protons the acid has
    num_protons = len(pKa_list)

    # build a list of terms where the term subscript corresponds to the list index
    terms_list = [10 ** (-pH * num_protons)]
    for i in range(len(pKa_list)):
        # multiply together all K values with
        k_term = 1
        for pK in pKa_list[0 : i + 1]:
            k_term *= 10**-pK
        terms_list.append(k_term * 10 ** (-pH * (num_protons - (i + 1))))
    print(terms_list)

    # return the desired distribution factor
    alpha = terms_list[n] / sum(terms_list)
    logger.debug(
        "Calculated %s-deprotonated acid distribution coefficient of %s for pKa=%s at pH %s % n,alpha,pKa_list,pH"
    )
    return alpha
