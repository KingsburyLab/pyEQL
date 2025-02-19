from __future__ import annotations

import copy
import os

import pymatgen.entries
import pytest
from pymatgen.core.composition import Composition
from pymatgen.entries.computed_entries import ComputedEntry
from pytest import approx

from pyEQL.pourbaix.compatibility import MU_H2O, CompatibilityError, MaterialsProjectAqueousCompatibility

PMG_ENTRIES_DIR = os.path.dirname(os.path.abspath(pymatgen.entries.__file__))

### Not sure to keep this line for the warning !!!!
# @pytest.mark.filterwarnings("ignore:MaterialsProjectCompatibility is deprecated")


class TestMaterialsProjectAqueousCompatibility:
    """
    Test MaterialsProjectAqueousCompatibility.

    -x- formation energy of H2O should always be -2.458 eV/H2O
    -x- H2 energy should always be the same value
    -x- H2O energy should always be the same value
    -x- Should get warnings if you init without all energy args
    -x- Should get CompatibilityError if you get_entry without all energy args
    -x- energy args should auto-populate from entries passed to process_entries
    -x- check compound entropies appropriately added
    -x- check hydrate adjustment appropriately applied

    Notes:
        Argument values from MaterialsProjectCompatibility as of April 2020:
            corrected DFT energy of H2O = -15.5875 eV/H2O (mp-697111) or -5.195 eV/atom
            corrected DFT energy of O2 = -4.9276 eV/atom (mp-12957)
            total energy corrections applied to H2O (eV/H2O) -0.70229 eV/H2O or -0.234 eV/atom
    """

    def test_h_h2o_energy_with_args_single(self):
        compat = MaterialsProjectAqueousCompatibility(
            o2_energy=-4.9276,
            h2o_energy=-5,
            h2o_adjustments=-0.234,
            solid_compat=None,
        )

        h2o_entry_1 = ComputedEntry(Composition("H2O"), -15)  # -5 eV/atom
        h2o_entry_2 = ComputedEntry(Composition("H4O2"), -6)  # -1 eV/atom
        h2_entry_1 = ComputedEntry(Composition("H8"), -100)  # -12.5 eV/atom
        h2_entry_2 = ComputedEntry(Composition("H2"), -16)  # -8 eV/atom

        for entry in [h2o_entry_1, h2o_entry_2]:
            compat.process_entries(entry)

        for entry in [h2_entry_1, h2_entry_2]:
            with pytest.warns(UserWarning, match="Processing single H2 entries"):
                compat.process_entries(entry)

        # the corrections should set the energy of any H2 polymorph the same, because
        # we have only processed one entry at time. Energy differences of H2O
        # polymorphs should be preserved.
        assert h2o_entry_2.energy_per_atom == approx(h2o_entry_1.energy_per_atom + 4)
        assert h2_entry_2.energy_per_atom == approx(h2_entry_1.energy_per_atom)

        o2_entry_1 = ComputedEntry(Composition("O2"), -4.9276 * 2)
        o2_entry_1 = compat.process_entries(o2_entry_1)[0]

        h2o_form_e = 3 * h2o_entry_1.energy_per_atom - 2 * h2_entry_2.energy_per_atom - o2_entry_1.energy_per_atom
        assert h2o_form_e == approx(MU_H2O)

    def test_h_h2o_energy_with_args_multi(self):
        compat = MaterialsProjectAqueousCompatibility(
            o2_energy=-4.9276,
            h2o_energy=-5,
            h2o_adjustments=-0.234,
            solid_compat=None,
        )

        h2o_entry_1 = ComputedEntry(Composition("H2O"), -15)  # -5 eV/atom
        h2o_entry_2 = ComputedEntry(Composition("H4O2"), -6)  # -1 eV/atom
        h2_entry_1 = ComputedEntry(Composition("H8"), -100)  # -12.5 eV/atom
        h2_entry_2 = ComputedEntry(Composition("H2"), -16)  # -8 eV/atom

        compat.process_entries([h2o_entry_1, h2o_entry_2, h2_entry_1, h2_entry_2])

        # Energy differences of H2O and H2 polymorphs should be preserved.
        assert h2o_entry_2.energy_per_atom == approx(h2o_entry_1.energy_per_atom + 4)
        assert h2_entry_2.energy_per_atom == approx(h2_entry_1.energy_per_atom + 4.5)

        o2_entry_1 = ComputedEntry(Composition("O2"), -4.9276 * 2)
        o2_entry_1 = compat.process_entries(o2_entry_1)[0]

        h2o_form_e = 3 * h2o_entry_1.energy_per_atom - 2 * h2_entry_1.energy_per_atom - o2_entry_1.energy_per_atom
        assert h2o_form_e == approx(MU_H2O)

    def test_h_h2o_energy_no_args(self):
        with pytest.warns(UserWarning, match="You did not provide the required O2 and H2O energies."):
            compat = MaterialsProjectAqueousCompatibility(solid_compat=None)

        h2o_entry_1 = ComputedEntry(Composition("H2O"), (-5.195 + 0.234) * 3, correction=-0.234 * 3)  # -5.195 eV/atom
        h2o_entry_2 = ComputedEntry(Composition("H4O2"), -6)  # -1 eV/atom
        h2_entry_1 = ComputedEntry(Composition("H8"), -100)  # -12.5 eV/atom``
        h2_entry_2 = ComputedEntry(Composition("H2"), -16)  # -8 eV/atom
        o2_entry_1 = ComputedEntry(Composition("O2"), -4.9276 * 2)

        with pytest.raises(CompatibilityError, match="Either specify the energies as arguments to "):
            compat.get_adjustments(h2_entry_1)

        compat.process_entries([h2o_entry_1, h2o_entry_2, h2_entry_1, h2_entry_2, o2_entry_1])

        assert compat.o2_energy == approx(-4.9276)
        assert compat.h2o_energy == approx(-5.195)
        assert compat.h2o_adjustments == approx(-0.234)

        # the corrections should preserve the difference in energy among H2O and H2 polymorphs
        assert h2o_entry_2.energy_per_atom == approx(h2o_entry_1.energy_per_atom + 4.195)
        assert h2_entry_2.energy_per_atom == approx(h2_entry_1.energy_per_atom + 4.5)

        # the water formation energy, calculated from the lowest energy polymorphs,
        # should equal the experimental value
        h2o_form_e = 3 * h2o_entry_1.energy_per_atom - 2 * h2_entry_1.energy_per_atom - o2_entry_1.energy_per_atom
        assert h2o_form_e == approx(MU_H2O)

    def test_compound_entropy(self):
        compat = MaterialsProjectAqueousCompatibility(
            o2_energy=-10, h2o_energy=-20, h2o_adjustments=-0.5, solid_compat=None
        )

        o2_entry_1 = ComputedEntry(Composition("O2"), -4.9276 * 2)

        initial_energy = o2_entry_1.energy_per_atom
        o2_entry_1 = compat.process_entries(o2_entry_1)[0]
        processed_energy = o2_entry_1.energy_per_atom

        assert initial_energy - processed_energy == approx(compat.cpd_entropies["O2"])

    def test_hydrate_adjustment(self):
        compat = MaterialsProjectAqueousCompatibility(
            o2_energy=-10, h2o_energy=-20, h2o_adjustments=-0.5, solid_compat=None
        )

        hydrate_entry = ComputedEntry(Composition("FeH4O2"), -10)  # nH2O = 2
        hydrate_entry2 = ComputedEntry(Composition("Li2O2H2"), -10)  # nH2O = 0

        compat.process_entries([hydrate_entry, hydrate_entry2])

        assert hydrate_entry.uncorrected_energy - hydrate_entry.energy == approx(
            2 * (compat.h2o_adjustments * 3 + MU_H2O)
        )
        assert hydrate_entry2.uncorrected_energy - hydrate_entry2.energy == 0

    def test_processing_entries_inplace(self):
        h2o_entry = ComputedEntry(Composition("H2O"), (-5.195 + 0.234) * 3, correction=-0.234 * 3)  # -5.195 eV/atom
        o2_entry = ComputedEntry(Composition("O2"), -4.9276 * 2)
        # check that compatibility scheme does not change input entries themselves
        entries = [h2o_entry, o2_entry]
        entries_copy = copy.deepcopy(entries)
        MaterialsProjectAqueousCompatibility().process_entries(entries, inplace=False)
        assert all(e.correction == e_copy.correction for e, e_copy in zip(entries, entries_copy, strict=True))

    def test_parallel_process_entries(self):
        hydrate_entry = ComputedEntry(Composition("FeH4O2"), -10)  # nH2O = 2
        hydrate_entry2 = ComputedEntry(Composition("Li2O2H2"), -10)  # nH2O = 0

        entry_list = [hydrate_entry, hydrate_entry2]

        compat = MaterialsProjectAqueousCompatibility(
            o2_energy=-10, h2o_energy=-20, h2o_adjustments=-0.5, solid_compat=None
        )

        with pytest.raises(
            ValueError,
            match="Parallel processing is not possible with for 'inplace=True'",
        ):
            entries = compat.process_entries(entry_list, inplace=True, n_workers=2)

        entries = compat.process_entries(entry_list, inplace=False, n_workers=2, on_error="raise")
        assert len(entries) == 2
