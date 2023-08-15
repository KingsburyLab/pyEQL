"""
pyEQL salt matching library.

This file contains functions that allow a pyEQL Solution object composed of
individual species (usually ions) to be mapped to a solution of one or more
salts. This mapping is necessary because some parameters (such as activity
coefficient data) can only be determined for salts (e.g. NaCl) and not individual
species (e.g. Na+)

:copyright: 2013-2023 by Ryan S. Kingsbury
:license: LGPL, see LICENSE for more details.

"""
from pymatgen.core.ion import Ion

from pyEQL.logging_system import logger


class Salt:
    """Class to represent a salt."""

    def __init__(self, cation, anion):
        """
        Create a Salt object based on its component ions.

        Parameters:
        ----------
        cation, anion : str
                Chemical formula of the cation and anion, respectively

        Returns:
        -------
        Salt : An object representing the properties of the salt

        Examples:
        --------
        >>> Salt('Na+','Cl-').formula
        'NaCl'

        >>> Salt('Mg++','Cl-').formula
        'MgCl2'

        """
        # create pymatgen Ion objects
        pmg_cat = Ion.from_formula(cation)
        pmg_an = Ion.from_formula(anion)
        # sanitize the cation and anion formulas
        self.cation = pmg_cat.reduced_formula
        self.anion = pmg_an.reduced_formula

        # get the charges on cation and anion
        self.z_cation = pmg_cat.charge
        self.z_anion = pmg_an.charge

        # assign stoichiometric coefficients by finding a common multiple
        self.nu_cation = int(abs(self.z_anion))
        self.nu_anion = int(abs(self.z_cation))

        # if both coefficients are the same, set each to one
        if self.nu_cation == self.nu_anion:
            self.nu_cation = 1
            self.nu_anion = 1

        # start building the formula, cation first
        salt_formula = ""
        if self.nu_cation > 1:
            # add parentheses if the cation is a polyatomic ion
            if len(pmg_cat.elements) > 1:
                salt_formula += "("
                salt_formula += self.cation.split("[")[0]
                salt_formula += ")"
            else:
                salt_formula += self.cation.split("[")[0]
            salt_formula += str(self.nu_cation)
        else:
            salt_formula += self.cation.split("[")[0]

        if self.nu_anion > 1:
            # add parentheses if the anion is a polyatomic ion
            if len(pmg_an.elements) > 1:
                salt_formula += "("
                salt_formula += self.anion.split("[")[0]
                salt_formula += ")"
            else:
                salt_formula += self.anion.split("[")[0]
            salt_formula += str(self.nu_anion)
        else:
            salt_formula += self.anion.split("[")[0]

        self.formula = salt_formula

    def get_effective_molality(self, ionic_strength):
        """Calculate the effective molality according to [#]_.

        .. math:: 2 I \\over (\\nu_+ z_+^2 + \\nu_- z_- ^2)

        Parameters
        ----------
        ionic_strength: Quantity
                        The ionic strength of the parent solution, mol/kg

        Returns:
        -------
        Quantity: the effective molality of the salt in the parent solution

        References:
        ----------
        .. [#] Mistry, K. H.; Hunter, H. a.; Lienhard V, J. H. Effect of
        composition and nonideal solution behavior on desalination calculations
        for mixed electrolyte solutions with comparison to seawater.
        Desalination 2013, 318, 34-47.
        """
        m_effective = 2 * ionic_strength / (self.nu_cation * self.z_cation**2 + self.nu_anion * self.z_anion**2)

        return m_effective.to("mol/kg")


def _sort_components(Solution, type="all"):
    """
    Sort the components of a solution in descending order (by mol).

    Parameters:
    ----------
    Solution : Solution object
    type     : The type of component to be sorted. Defaults to 'all' for all
                solutes. Other valid arguments are 'cations' and 'anions' which
                return sorted lists of cations and anions, respectively.

    Returns:
    -------
    A list whose keys are the component names (formulas) and whose
    values are the component objects themselves


    """
    formula_list = []

    # populate a list with component names
    for item in Solution.components:
        z = Solution.get_property(item, "charge")
        if type == "all" or (type == "cations" and z > 0) or (type == "anions" and z < 0):
            formula_list.append(item)

    # populate a dictionary with formula:concentration pairs
    mol_list = {item: Solution.get_amount(item, "mol") for item in formula_list}

    return sorted(formula_list, key=mol_list.__getitem__, reverse=True)


def identify_salt(sol):
    """
    Analyze the components of a solution and identify the salt that most closely
    approximates it.
    (e.g., if a solution contains 0.5 mol/kg of Na+ and Cl-, plus traces of H+
    and OH-, the matched salt is 0.5 mol/kg NaCl).

    Create a Salt object for this salt.

    Returns:
    -------
    A Salt object.
    """
    # sort the components by moles
    sort_list = _sort_components(sol)

    # default to returning water as the salt
    cation = "H+"
    anion = "OH-"

    # return water if there are no solutes
    if len(sort_list) < 3 and sort_list[0] == "H2O":
        logger.info("Salt matching aborted because there are not enough solutes.")
        return Salt(cation, anion)

    # warn if something other than water is the predominant component
    if sort_list[0] != "H2O":
        logger.warning("H2O is not the most prominent component")

    # take the dominant cation and anion and assemble a salt from them
    for item in sort_list:
        pmg_ion = Ion.from_formula(item)
        if pmg_ion.charge > 0 and cation == "H+":
            cation = item
        elif pmg_ion.charge < 0 and anion == "OH-":
            anion = item
        else:
            pass

    # assemble the salt
    return Salt(cation, anion)


def generate_salt_list(sol, unit="mol/kg"):
    """
    Generate a list of salts that represents the ionic composition of a
    solution.

    Returns:
    -------
    dict
        A dictionary of Salt objects, where Salt objects are the keys and
        the amounts are the values.

    """
    salt_list = {}

    # sort the cations and anions by moles
    cation_list = _sort_components(sol, type="cations")
    anion_list = _sort_components(sol, type="anions")

    # iterate through the lists of ions
    # create salts by matching the equivalent concentrations of cations
    # and anions along the way
    len_cat = len(cation_list)
    len_an = len(anion_list)

    # start with the first cation and anion
    index_cat = 0
    index_an = 0

    # TODO - add an equivalent concnetration method to get_amount
    # calculate the equivalent concentrations of each ion
    c1 = sol.get_amount(cation_list[index_cat], unit) * sol.get_property(cation_list[index_cat], "charge")
    a1 = sol.get_amount(anion_list[index_an], unit) * abs(sol.get_property(anion_list[index_an], "charge"))

    while index_cat < len_cat and index_an < len_an:
        # if the cation concentration is greater, there will be leftover cations
        if c1 > a1:
            # create the salt
            x = Salt(cation_list[index_cat], anion_list[index_an])
            # there will be leftover cation, so use the anion amount
            amount = a1 / abs(x.z_anion)
            # add it to the list
            salt_list.update({x: amount})
            # adjust the amounts of the respective ions
            c1 = c1 - a1
            # move to the next anion
            index_an += 1
            try:
                a1 = sol.get_amount(anion_list[index_an], unit) * abs(sol.get_property(anion_list[index_an], "charge"))
            except IndexError:
                continue
        # if the anion concentration is greater, there will be leftover anions
        if c1 < a1:
            # create the salt
            x = Salt(cation_list[index_cat], anion_list[index_an])
            # there will be leftover anion, so use the cation amount
            amount = c1 / x.z_cation
            # add it to the list
            salt_list.update({x: amount})
            # calculate the leftover cation amount
            a1 = a1 - c1
            # move to the next cation
            index_cat += 1
            try:
                c1 = sol.get_amount(cation_list[index_cat], unit) * sol.get_property(cation_list[index_cat], "charge")
            except IndexError:
                continue
        if c1 == a1:
            # create the salt
            x = Salt(cation_list[index_cat], anion_list[index_an])
            # there will be nothing leftover, so it doesn't matter which ion you use
            amount = c1 / x.z_cation
            # add it to the list
            salt_list.update({x: amount})
            # move to the next cation and anion
            index_an += 1
            index_cat += 1
            try:
                c1 = sol.get_amount(cation_list[index_cat], unit) * sol.get_property(cation_list[index_cat], "charge")
                a1 = sol.get_amount(anion_list[index_an], unit) * abs(sol.get_property(anion_list[index_an], "charge"))
            except IndexError:
                continue

    return salt_list
