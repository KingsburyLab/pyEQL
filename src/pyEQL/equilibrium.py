"""
pyEQL methods for chemical equilibrium calculations (e.g. acid/base, reactions,
redox, complexation, etc.).

NOTE: these methods are not currently used but are here for the future.

:copyright: 2013-2023 by Ryan S. Kingsbury
:license: LGPL, see LICENSE for more details.

"""
# import libraries for scientific functions
import math
import os
from pathlib import Path
from typing import Literal

from phreeqpython import PhreeqPython

# the pint unit registry
from pyEQL import ureg
from pyEQL.logging_system import logger
from pyEQL.utils import standardize_formula

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

# These are the only elements that are allowed to have parenthetical oxidation states
# PHREEQC will ignore others (e.g., 'Na(1)')
SPECIAL_ELEMENTS = ["S", "C", "N", "Cu", "Fe", "Mn"]


def equilibrate_phreeqc(
    solution,
    phreeqc_db: Literal["vitens.dat", "wateq4f_PWN.dat", "pitzer.dat", "llnl.dat", "geothermal.dat"] = "vitens.dat",
):
    """Adjust the speciation of a Solution object to achieve chemical equilibrium.

    Args:
        phreeqc_db: Name of the PHREEQC database file to use for solution thermodynamics
                and speciation calculations. Generally speaking, `llnl.dat` is recommended
                for moderate salinity water and prediction of mineral solubilities,
                `wateq4f_PWN.dat` is recommended for low to moderate salinity waters. It is
                similar to vitens.dat but has many more species. `pitzer.dat` is recommended
                when accurate activity coefficients in solutions above 1 M TDS are desired, but
                it has fewer species than the other databases. `llnl.dat` and `geothermal.dat`
                may offer improved prediction of LSI.
    """
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
        "density": solution.density.to("g/mL").magnitude,
    }
    balance_charge = solution.balance_charge
    if balance_charge == "pH":
        d["pH"] = str(d["pH"]) + " charge"
    if balance_charge == "pE":
        d["pe"] = str(d["pe"]) + " charge"
    initial_comp = solution.components.copy()

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

        # tell PHREEQC which species to use for charge balance
        if el == balance_charge:
            key += " charge"
        d[key] = mol / solv_mass

    # database files in this list are not distributed with phreeqpython
    db_path = Path(os.path.dirname(__file__)) / "database" if phreeqc_db in ["llnl.dat", "geothermal.dat"] else None
    # create the PhreeqcPython instance
    pp = PhreeqPython(database=phreeqc_db, database_directory=db_path)

    # # equalize with atmospheric air (optional)
    # if EQUALIZE:
    #     phases = [("CO2(g)", -3.5), ("O2(g)", -0.67)]
    #     self.ppsol.equalize([t[0] for t in phases], [t[1] for t in phases])

    # create the PHREEQC solution object
    try:
        ppsol = pp.add_solution(d)
    except Exception as e:
        print(d)
        # catch problems with the input to phreeqc
        raise ValueError(
            "There is a problem with your input. The error message received from "
            f" phreeqpython is:\n\n {e}\n Check your input arguments, especially "
            "the composition dictionary, and try again."
        )

    # use the output from PHREEQC to update the Solution composition
    # the .species attribute should return MOLES (not moles per ___)
    for s, mol in ppsol.species.items():
        solution.components[s] = mol

    # make sure all species are accounted for
    assert set(initial_comp.keys()) - set(solution.components.keys()) == set()

    # log a message if any components were not touched by PHREEQC
    missing_species = set(initial_comp.keys()) - {standardize_formula(s) for s in ppsol.species}
    if len(missing_species) > 0:
        logger.info(
            f"After equilibration, the amounts of species {missing_species} were not modified "
            "by PHREEQC. These species are likely absent from its database."
        )

    # remove the PPSol from the phreeqcpython instance
    pp.remove_solutions([0])


def adjust_temp_pitzer(c1, c2, c3, c4, c5, temp, temp_ref=ureg.Quantity("298.15 K")):
    """
    Calculate a parameter for the Pitzer model based on temperature-dependent
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
    r"""(float,float,number, optional number) -> float.

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
    over the range of interest.

    References
    ----------
    Stumm, Werner and Morgan, James J. Aquatic Chemistry, 3rd ed, pp 53.
        Wiley Interscience, 1996.

    Examples:
    --------
    >>> adjust_temp_vanthoff(0.15,ureg.Quantity('-197.6 kJ/mol'),ureg.Quantity('42 degC'),ureg.Quantity(' 25degC')) #doctest: +ELLIPSIS
    0.00203566...

    If the 'ref_temperature' parameter is omitted, a default of 25 C is used.

    >>> adjust_temp_vanthoff(0.15,ureg.Quantity('-197.6 kJ/mol'),ureg.Quantity('42 degC')) #doctest: +ELLIPSIS
    0.00203566...

    """
    output = equilibrium_constant * math.exp(
        enthalpy / ureg.R * (1 / reference_temperature.to("K") - 1 / temperature.to("K"))
    )

    logger.info(
        "Adjusted equilibrium constant K=%s from %s to %s degrees Celsius with Delta H = %s. Adjusted K = %s % equilibrium_constant,reference_temperature,temperature,enthalpy,output"
    )

    logger.warning("Van't Hoff equation assumes enthalpy is independent of temperature over the range of interest")
    return output


def adjust_temp_arrhenius(
    rate_constant,
    activation_energy,
    temperature,
    reference_temperature=ureg.Quantity(25, "degC"),
):
    r"""(float,float,number, optional number) -> float.

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

    See Also:
    --------
    kelvin

    Notes
    -----
    This function implements the Arrhenius equation to adjust measured rate
    constants to other temperatures.
    TODO - add better reference

    .. math::
        ln(\\frac{K2}{K1} = \\frac{E_a}{R} ( \\frac{1}{T_{1}} - {\\frac{1}{T_2}} )

    References
    ----------

    http://chemwiki.ucdavis.edu/Physical_Chemistry/Kinetics/Reaction_Rates/Temperature_Dependence_of_Reaction_Rates/Arrhenius_Equation


    Examples:
    --------
    >>> adjust_temp_arrhenius(7,900*ureg.Quantity('kJ/mol'),37*ureg.Quantity('degC'),97*ureg.Quantity('degC')) #doctest: +ELLIPSIS
    1.8867225...e-24

    """
    output = rate_constant * math.exp(
        activation_energy / ureg.R * (1 / reference_temperature.to("K") - 1 / temperature.to("K"))
    )

    logger.info(
        "Adjusted parameter %s from %s to %s degrees Celsius with Activation Energy = %s kJ/mol. Adjusted value = %s % rate_constant,reference_temperature,temperature,activation_energy,output"
    )

    return output


def alpha(n, pH, pKa_list):
    """Returns the acid-base distribution coefficient (alpha) of an acid in the
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
        The acid-base cient is calculated as follows: [stm]_

        .. math::

            \\alpha_n = \\frac{term_n}{[H+]^n + k_{a1}[H+]^{n-1} + k_{a1}k_{a2}[H+]^{n-2} ... k_{a1}k_{a2}...k_{an} }

        Where :math: '\term_n' refers to the nth term in the denominator, starting from 0

    References
    ----------
        .. [stm] Stumm, Werner and Morgan, James J. Aquatic Chemistry, 3rd ed, pp 127-130. Wiley Interscience, 1996.

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

    """
    # generate an error if no pKa values are specified
    if len(pKa_list) == 0:
        logger.error("No pKa values given. Cannot calculate distribution coeffiicent.")
        return None

    # generate an error if n > number of pKa values
    if len(pKa_list) < n:
        logger.error("Insufficient number of pKa values given. Cannot calculate distribution coeffiicent.")
        return None

    # convert pH to hydrogen ion concentration
    Hplus = 10**-pH

    # determine how many protons the acid has
    num_protons = len(pKa_list)

    # build a list of terms where the term subscript corresponds to the list index
    terms_list = []
    k_term = 1

    # the 'item' index counts from 0 to the number of protons, inclusive
    for item in range(num_protons + 1):
        # multiply the preceding k values together
        for i in range(len(pKa_list[:item])):
            k_term *= 10 ** -pKa_list[i]

        # add the term to the list
        terms_list.append(k_term * Hplus ** (num_protons - item))

    # build the expression
    numerator = terms_list[n]
    denominator = 0
    for item in terms_list:
        denominator += item

    # return the desired distribution factor
    alpha = numerator / denominator
    logger.info(
        "Calculated %s-deprotonated acid distribution coefficient of %s for pKa=%s at pH %s % n,alpha,pKa_list,pH"
    )
    return alpha
