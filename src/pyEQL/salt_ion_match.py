"""
pyEQL salt matching library.

This file contains functions that allow a pyEQL Solution object composed of
individual species (usually ions) to be mapped to a solution of one or more
salts. This mapping is necessary because some parameters (such as activity
coefficient data) can only be determined for salts (e.g. NaCl) and not individual
species (e.g. Na+)

:copyright: 2013-2024 by Ryan S. Kingsbury
:license: LGPL, see LICENSE for more details.

"""

from monty.json import MSONable
from pymatgen.core.ion import Ion

from pyEQL.utils import standardize_formula


class Salt(MSONable):
    """Class to represent a salt."""

    def __init__(self, cation, anion):
        """
        Create a Salt object based on its component ions.

        Parameters:
            cation, anion: (str) Chemical formula of the cation and anion, respectively.

        Returns:
            Salt : An object representing the properties of the salt

        Examples:
            >>> Salt('Na+','Cl-').formula
            'NaCl'

            >>> Salt('Mg++','Cl-').formula
            'MgCl2'
        """
        # create pymatgen Ion objects
        pmg_cat = Ion.from_formula(cation)
        pmg_an = Ion.from_formula(anion)
        # standardize the cation and anion formulas
        self.cation = standardize_formula(cation)
        self.anion = standardize_formula(anion)

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

    # TODO - consider whether this should be adjusted to be based on total concentrations or not
    # NOTE: speciating the solution results in a decrease in the overall ionic strength, because some of the
    # Mg+2 is converted to monovalent complexes like MgOH+. Hence, the activity coefficients deviate a bit from
    # the published values.
    def get_effective_molality(self, ionic_strength):
        r"""Calculate the effective molality according to [mistry]_.

        .. math:: 2 I \over (\nu_+ z_+^2 + \nu_- z_- ^2)

        Args:
            ionic_strength: Quantity
                The ionic strength of the parent solution, mol/kg

        Returns:
            Quantity: the effective molality of the salt in the parent solution

        References:
            .. [mistry] Mistry, K. H.; Hunter, H. a.; Lienhard V, J. H. Effect of composition and nonideal solution behavior
                on desalination calculations for mixed electrolyte solutions with comparison to seawater. Desalination
                2013, 318, 34-47.
        """
        m_effective = 2 * ionic_strength / (self.nu_cation * self.z_cation**2 + self.nu_anion * self.z_anion**2)

        return m_effective.to("mol/kg")
