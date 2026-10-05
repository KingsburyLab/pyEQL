"""This module implements Compatibility corrections for mixing runs of different
functionals.
"""

from __future__ import annotations

import copy
import os
import warnings
from importlib.resources import files
from typing import TYPE_CHECKING, TypeAlias

import numpy as np
from monty.design_patterns import cached_class
from monty.serialization import loadfn
from pymatgen.entries.compatibility import Compatibility, CompatibilityError, MaterialsProject2020Compatibility
from pymatgen.entries.computed_entries import (
    CompositionEnergyAdjustment,
    ComputedEntry,
    ComputedStructureEntry,
    ConstantEnergyAdjustment,
    EnergyAdjustment,
    TemperatureEnergyAdjustment,
)
from pymatgen.util.due import Doi, due

if TYPE_CHECKING:
    from typing import Literal


__author__ = "Amanda Wang, Ryan Kingsbury, Shyue Ping Ong, Anubhav Jain, Stephen Dacek, Sai Jayaraman"
__copyright__ = "Copyright 2012-2020, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "April 2020"

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
MU_H2O = -2.4583  # Free energy of formation of water, eV/H2O, used by MaterialsProjectAqueousCompatibility
MU_CO2 = -4.0004  # Free energy of formation of CO2, eV/CO2
MU_C = MU_CO2 - 2 * MU_H2O
ENTROPY_DATABASE = loadfn(files("pyEQL") / "pourbaix" / "phonon_database.json")

AnyComputedEntry: TypeAlias = ComputedEntry | ComputedStructureEntry


@cached_class
@due.dcite(
    Doi(
        "10.1103/PhysRevB.85.235438",
        "Pourbaix scheme to combine calculated and experimental data",
    )
)
class MaterialsProjectAqueousCompatibility(Compatibility):
    """This class implements the Aqueous energy referencing scheme for constructing
    Pourbaix diagrams from DFT energies, as described in Persson et al.

    This scheme applies various energy adjustments to convert DFT energies into
    Gibbs free energies of formation at 298 K and to guarantee that the experimental
    formation free energy of H2O is reproduced. Briefly, the steps are:

        1. Beginning with the DFT energy of O2, adjust the energy of H2 so that
           the experimental reaction energy of -2.458 eV/H2O is reproduced.
        2. Add entropy to the DFT energy of any compounds that are liquid or
           gaseous at room temperature
        3. Adjust the DFT energies of solid hydrate compounds (compounds that
           contain water, e.g. FeO.nH2O) such that the energies of the embedded
           H2O molecules are equal to the experimental free energy

    The above energy adjustments are computed dynamically based on the input
    Entries.

    References:
        K.A. Persson, B. Waldwick, P. Lazic, G. Ceder, Prediction of solid-aqueous
        equilibria: Scheme to combine first-principles calculations of solids with
        experimental aqueous states, Phys. Rev. B - Condens. Matter Mater. Phys. # codespell:ignore=Mater
        85 (2012) 1-12. doi:10.1103/PhysRevB.85.235438.
    """

    def __init__(
        self,
        solid_compat: Compatibility | type[Compatibility] | None = MaterialsProject2020Compatibility,
        o2_energy: float | None = None,
        h2o_energy: float | None = None,
        h2o_adjustments: float | None = None,
        universal_solid_shift_eV_per_atom: float = -0.055,
    ) -> None:
        """Initialize the MaterialsProjectAqueousCompatibility class.

        Note that this class requires as inputs the ground-state DFT energies of O2 and H2O, plus the value of any
        energy adjustments applied to an H2O molecule. If these parameters are not provided in __init__, they can
        be automatically populated by including ComputedEntry for the ground state of O2 and H2O in a list of
        entries passed to process_entries. process_entries will fail if one or the other is not provided.

        Args:
            solid_compat: Compatibility scheme used to pre-process solid DFT energies prior to applying aqueous
                energy adjustments. May be passed as a class (e.g. MaterialsProject2020Compatibility) or an instance
                (e.g., MaterialsProject2020Compatibility()). If None, solid DFT energies are used as-is.
                Default: MaterialsProject2020Compatibility
            o2_energy: The ground-state DFT energy of oxygen gas, including any adjustments or corrections, in eV/atom.
                If not set, this value will be determined from any O2 entries passed to process_entries.
                Default: None
            h2o_energy: The ground-state DFT energy of water, including any adjustments or corrections, in eV/atom.
                If not set, this value will be determined from any H2O entries passed to process_entries.
                Default: None
            h2o_adjustments: Total energy adjustments applied to one water molecule, in eV/atom.
                If not set, this value will be determined from any H2O entries passed to process_entries.
                Default: None
            universal_solid_shift_eV_per_atom: Uniform energy correction applied to all solid compounds, in eV/atom, to account for systematic errors in DFT formation energies relative to the Wagman-NBS tabulated free energies of formation. This correction is specific to pyEQL and is not part of pymatgen's MaterialsProjectAqueousCompatibility. To disable, set to 0. Default: 0.055.
        """
        self.solid_compat = None
        # check whether solid_compat has been instantiated
        if solid_compat is None:
            self.solid_compat = None
        elif isinstance(solid_compat, type) and issubclass(solid_compat, Compatibility):
            self.solid_compat = solid_compat()
        elif issubclass(type(solid_compat), Compatibility):
            self.solid_compat = solid_compat
        else:
            raise ValueError("Expected a Compatibility class, instance of a Compatibility or None")

        self.o2_energy = o2_energy
        self.h2o_energy = h2o_energy
        self.h2_energy = None
        self.h2o_adjustments = h2o_adjustments

        if not all([self.o2_energy, self.h2o_energy, self.h2o_adjustments]):
            warnings.warn(
                f"You did not provide the required O2 and H2O energies. {type(self).__name__} "
                "needs these energies in order to compute the appropriate energy adjustments. It will try "
                "to determine the values from ComputedEntry for O2 and H2O passed to process_entries, but "
                "will fail if these entries are not provided.",
                stacklevel=2,
            )

        # Standard state entropy of pure elements, molecular, gases, and reference solids at 298K (-T delta S)
        # from Wagman-NBS tables and UMA phonon calculations (eV/atom)
        self.cpd_entropies = {
            # exp anion entropy at 300 K
            "O2": 0.311731,
            "N2": 0.300729,
            "F2": 0.315562,
            "Cl2": 0.339373,
            "Br": 0.235039,
            "Hg": 0.234421,
            "H2O": 0.076963,  # 0.23079 eV/H2O
        }

        # pymatgen's reference entropies that must not be overritten by either the UMA phonon calculations or the NIST-NBS experimental data:
        mp_entropies = {
            "H2",
            "O2",
            "N2",
            "F2",
            "Cl2",
            "Br",
            "Hg",
            "H2O",
        }

        entropy_lib = [
            ("uMLIP", "uma-s-1p1"),
            ("Experiment", "NIST-NBS"),
        ]

        for etr in ENTROPY_DATABASE:
            formula = etr["formula"]

            if formula in mp_entropies:
                continue

            entropy_data = etr.get("entropy_per_atom", {})

            entropy = next(
                (
                    entropy_data.get(source, {}).get(method, {}).get("entropy_eV_per_atom")
                    for source, method in entropy_lib
                    if (entropy_data.get(source, {}).get(method, {}).get("entropy_eV_per_atom")) is not None
                ),
                None,
            )

            if entropy is not None:
                self.cpd_entropies[formula] = entropy

        self.name = "MP Aqueous free energy adjustment"
        super().__init__()

        self.universal_solid_shift_eV_per_atom = universal_solid_shift_eV_per_atom

    def get_adjustments(self, entry: ComputedEntry) -> list[EnergyAdjustment]:
        """Get the Aqueous corrections applied to DFT entries.

        Args:
            entry: A ComputedEntry object.

        Returns:
            list[EnergyAdjustment]: Energy adjustments to be applied to entry.

        Raises:
            CompatibilityError if the required O2 and H2O energies have not been provided to
            MaterialsProjectAqueousCompatibility during init or in the list of entries passed to process_entries.
        """
        adjustments = []
        if self.o2_energy is None or self.h2o_energy is None or self.h2o_adjustments is None:
            raise CompatibilityError(
                "You did not provide the required O2 and H2O energies. "
                f"{type(self).__name__} needs these energies in order to compute "
                "the appropriate energy adjustments. Either specify the energies as arguments "
                f"to {type(self).__name__}.__init__ or run process_entries on a list that includes ComputedEntry for "
                "the ground state of O2 and H2O."
            )

        # compute the free energies of H2 and H2O (eV/atom) to guarantee that the
        # formation-free energy of H2O is equal to -2.4583 eV/H2O from experiments
        # (MU_H2O from Pourbaix module)

        # Free energy of H2 in eV/atom, fitted using Eq. 40 of Persson et al. PRB 2012 85(23)
        # https://journals.aps.org/prb/abstract/10.1103/PhysRevB.85.235438
        # for this calculation ONLY, we need the (corrected) DFT energy of water
        self.fit_h2_energy = round(
            0.5
            * (
                3 * (self.h2o_energy - self.cpd_entropies["H2O"]) - (self.o2_energy - self.cpd_entropies["O2"]) - MU_H2O
            ),
            6,
        )

        comp = entry.composition
        rform = comp.reduced_formula

        # use fit_h2_energy to adjust the energy of all H2 polymorphs such that
        # the lowest energy polymorph has the correct experimental value
        # if H2O and O2 energies have been set explicitly via kwargs, then
        # all H2 polymorphs will get the same energy.
        if rform == "H2":
            if self.h2_energy is None:
                raise ValueError("H2 energy not set")
            adjustments.append(
                ConstantEnergyAdjustment(
                    (self.fit_h2_energy - self.h2_energy) * comp.num_atoms,
                    uncertainty=np.nan,
                    name="MP Aqueous H2 / H2O referencing",
                    cls=self.as_dict(),
                    description="Adjusts the H2 energy to reproduce the experimental "
                    "Gibbs formation free energy of H2O, based on the DFT energy "
                    "of Oxygen and H2O",
                )
            )

        # add minus T delta S to the DFT energy (enthalpy) of compounds that are
        # molecular-like at room temperature
        if rform in self.cpd_entropies:
            adjustments.append(
                TemperatureEnergyAdjustment(
                    -self.cpd_entropies[rform] / 300,
                    300,
                    comp.num_atoms,
                    uncertainty_per_deg=np.nan,
                    name="Compound entropy at room temperature",
                    cls=self.as_dict(),
                    description="Adds the entropy (T delta S) to energies of compounds that "
                    "are gaseous or liquid at standard state",
                )
            )

        MU_N_CORRECTION = 0.26
        # For nitrogen compounds, we apply a correction to the DFT energy
        if rform != "N2" and "N" in comp:
            n_N = comp["N"]
            n_correction = MU_N_CORRECTION

            adjustments.append(
                CompositionEnergyAdjustment(
                    n_correction,
                    n_N,
                    uncertainty_per_atom=np.nan,
                    name="MP Aqueous Nitrogen correction",
                    cls=self.as_dict(),
                    description="Adjust the energy of solid nitrogen compounds so that the"
                    "free energies match the experimental"
                    " value enforced by the MP Aqueous energy referencing scheme.",
                )
            )

        # Universal correction for solid/compounds
        if self.universal_solid_shift_eV_per_atom:
            is_element = comp.is_element

            molecular_like_rforms = {"O2", "N2", "F2", "Cl2", "Br", "Hg"}
            is_molecular_standard_state = rform in molecular_like_rforms
            is_special_ref = rform in {"H2", "H2O", "O2"}

            # universal solid shift is applied to compounds, except for reference elements and molecular standard states
            apply_shift = (not is_element) and (not is_molecular_standard_state) and (not is_special_ref)

            if apply_shift:
                total_shift = self.universal_solid_shift_eV_per_atom * comp.num_atoms

                adjustments.append(
                    ConstantEnergyAdjustment(
                        total_shift,
                        uncertainty=np.nan,
                        name="User universal solid shift (eV/atom)",
                        cls=self.as_dict(),
                        description=(
                            f"Applies a user-defined shift of {self.universal_solid_shift_eV_per_atom:+.4f} eV/atom "
                            f"({total_shift:+.4f} eV total) to selected entries."
                        ),
                    )
                )

        return adjustments

    def process_entries(
        self,
        entries: list[AnyComputedEntry],
        clean: bool = False,
        verbose: bool = False,
        inplace: bool = True,
        n_workers: int = 1,
        on_error: Literal["ignore", "warn", "raise"] = "ignore",
    ) -> list[AnyComputedEntry]:
        """Process a sequence of entries with the chosen Compatibility scheme.

        Args:
            entries (list[ComputedEntry | ComputedStructureEntry]): Entries to be processed.
            clean (bool): Whether to remove any previously-applied energy adjustments.
                If True, all EnergyAdjustment are removed prior to processing the Entry.
                Default is False.
            verbose (bool): Whether to display progress bar for processing multiple entries.
                Default is False.
            inplace (bool): Whether to modify the entries in place. If False, a copy of the
                entries is made and processed. Default is True.
            n_workers (int): Number of workers to use for parallel processing. Default is 1.
            on_error ('ignore' | 'warn' | 'raise'): What to do when get_adjustments(entry)
                raises CompatibilityError. Defaults to 'ignore'.

        Returns:
            list[AnyComputedEntry]: Adjusted entries. Entries in the original list incompatible with
                chosen correction scheme are excluded from the returned list.
        """
        # Convert input arg to a list if not already
        if isinstance(entries, ComputedEntry):
            entries = [entries]

        # If not inplace, process entries on a copy
        if not inplace:
            entries = copy.deepcopy(entries)

        # Pre-process entries with the given solid compatibility class
        if self.solid_compat:
            entries = self.solid_compat.process_entries(entries, clean=True, inplace=inplace, n_workers=n_workers)

        # when processing single entries, all H2 polymorphs will get assigned the
        # same energy
        if len(entries) == 1 and entries[0].reduced_formula == "H2":
            warnings.warn(
                "Processing single H2 entries will result in the all polymorphs "
                "being assigned the same energy. This should not cause problems "
                "with Pourbaix diagram construction, but may be confusing. "
                "Pass all entries to process_entries() at once in if you want to "
                "preserve H2 polymorph energy differences.",
                stacklevel=2,
            )

        # extract the DFT energies of oxygen and water from the list of entries, if present
        # do not do this when processing a single entry, as it might lead to unintended
        # results
        if len(entries) > 1:
            if not self.o2_energy and (o2_entries := [e for e in entries if e.reduced_formula == "O2"]):
                self.o2_energy = min(e.energy_per_atom for e in o2_entries)

            if not self.h2o_energy and not self.h2o_adjustments:  # noqa: SIM102
                if h2o_entries := [e for e in entries if e.reduced_formula == "H2O"]:
                    h2o_entries = sorted(h2o_entries, key=lambda e: e.energy_per_atom)
                    self.h2o_energy = h2o_entries[0].energy_per_atom
                    self.h2o_adjustments = h2o_entries[0].correction / h2o_entries[0].composition.num_atoms

        if h2_entries := [e for e in entries if e.reduced_formula == "H2"]:
            h2_entries = sorted(h2_entries, key=lambda e: e.energy_per_atom)
            self.h2_energy = h2_entries[0].energy_per_atom  # type: ignore[assignment]

        return super().process_entries(
            entries,
            clean=clean,
            verbose=verbose,
            inplace=inplace,
            n_workers=n_workers,
            on_error=on_error,
        )
