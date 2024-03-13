"""
pyEQL Solute class.

This file contains functions and methods for managing properties of
individual solutes. The Solute class contains methods for accessing
ONLY those properties that DO NOT depend on solution composition.
Solute properties such as activity coefficient or concentration
that do depend on compsition are accessed via Solution class methods.

:copyright: 2013-2024 by Ryan S. Kingsbury
:license: LGPL, see LICENSE for more details.

"""

from __future__ import annotations

import logging
import warnings
from dataclasses import asdict, dataclass, field
from typing import Literal

import numpy as np
from pymatgen.core.ion import Ion

from pyEQL.utils import standardize_formula

logger = logging.getLogger(__name__)


@dataclass
class Datum:
    """Document containing data for a single computed or experimental property."""

    value: str
    reference: str | None = None
    data_type: Literal["computed", "experimental", "fitted", "unknown"] = "unknown"

    @property
    def magnitude(self):
        """Return the numerical value of a Datum."""
        return float(self.value.split(" ")[0])

    @property
    def unit(self):
        """Return the unit of a Datum."""
        return self.value.split(" ")[-1]

    @property
    def uncertainty(self):
        """Return the uncertainty of a Datum."""
        if len(self.value.split(" ")) > 3:
            return float(self.value.split(" ")[2])
        return np.nan

    def as_dict(self):
        """Return a dictionary representation of the Datum."""
        return dict(asdict(self).items())


@dataclass
class Solute:
    """
    represent each chemical species as an object containing its formal charge,
    transport numbers, concentration, activity, etc.

    Args:
        formula: Chemical formula for the solute. Charged species must contain a + or - and (for polyvalent solutes)
            a number representing the net charge (e.g. 'SO4-2').
    """

    formula: str
    charge: int
    molecular_weight: str
    elements: list
    chemsys: str
    pmg_ion: Ion
    formula_html: str
    formula_latex: str
    formula_hill: str
    formula_pretty: str
    oxi_state_guesses: dict[str, float]
    n_atoms: int
    n_elements: int
    size: dict = field(
        default_factory=lambda: {
            "radius_ionic": None,
            "radius_hydrated": None,
            "radius_vdw": None,
            "molar_volume": None,
        }
    )
    thermo: dict = field(default_factory=lambda: {"ΔG_hydration": None, "ΔG_formation": None})
    transport: dict = field(default_factory=lambda: {"diffusion_coefficient": None})
    model_parameters: dict = field(
        default_factory=lambda: {
            "activity_pitzer": {"Beta0": None, "Beta1": None, "Beta2": None, "Cphi": None, "Max_C": None},
            "molar_volume_pitzer": {
                "Beta0": None,
                "Beta1": None,
                "Beta2": None,
                "Cphi": None,
                "V_o": None,
                "Max_C": None,
            },
            "viscosity_jones_dole": {"B": None},
            "diffusion_temp_smolyakov": {"a1": None, "a2": None, "d": None},
        }
    )

    @classmethod
    def from_formula(cls, formula: str):
        """
        Create an Ion document from a chemical formula. The formula is passed to
        pymatgen.core.Ion.from_formula() and used to populate the basic chemical
        informatics fields (e.g., formula, charge, molecular weight, elements, etc.)
        of the IonDoc.
        """
        pmg_ion = Ion.from_formula(formula)
        f, factor = pmg_ion.get_reduced_formula_and_factor()
        rform = standardize_formula(formula)
        charge = int(pmg_ion.charge)
        els = [str(el) for el in pmg_ion.elements]
        mw = f"{float(pmg_ion.weight / factor)} g/mol"  # weight is a FloatWithUnit
        chemsys = pmg_ion.chemical_system
        # store only the most likely oxi_state guesses
        try:
            oxi_states = pmg_ion.oxi_state_guesses(all_oxi_states=True)[0]
        except (IndexError, ValueError):
            warnings.warn(f"Guessing oxi states failed for {formula}")
            oxi_states = {}

        return cls(
            rform,
            charge=charge,
            molecular_weight=mw,
            elements=els,
            chemsys=chemsys,
            pmg_ion=pmg_ion,
            formula_html=pmg_ion.to_html_string(),
            formula_latex=pmg_ion.to_latex_string(),
            formula_hill=pmg_ion.hill_formula,
            formula_pretty=pmg_ion.to_pretty_string(),
            oxi_state_guesses=oxi_states,
            n_atoms=int(pmg_ion.num_atoms),
            n_elements=len(els),
        )

    def as_dict(self):
        """Return a dictionary representation of the Solute."""
        return dict(asdict(self).items())

    # set output of the print() statement
    def __str__(self):
        return (
            "Species "
            + str(self.formula)
            + " MW="
            + str(self.mw)
            + " Formal Charge="
            + str(self.charge)
            + " Amount= "
            + str(self.moles)
        )
