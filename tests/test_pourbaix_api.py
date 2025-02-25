import itertools
import os
import random
import importlib

import numpy as np
import pytest

from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.analysis.pourbaix_diagram import IonEntry, PourbaixDiagram, PourbaixEntry
from pymatgen.core.ion import Ion
from pymatgen.entries.compatibility import MaterialsProjectAqueousCompatibility

from mp_api.client import MPRester

from pyEQL.pourbaix.pourbaix_api import Pourbaix_api

@pytest.fixture()
def mpr():
    rester = MPRester()
    yield rester
    rester.session.close()


@pytest.mark.skipif(os.getenv("MP_API_KEY", None) is None, reason="No API key found.")
class TestMPRester:
    fake_mp_api_key = "12345678901234567890123456789012"  
    default_endpoint = "https://api.materialsproject.org/"

    @pytest.mark.skip(reason="SSL issues")
    def test_get_ion_entries(self, mpr):
        entries = mpr.get_entries_in_chemsys("Ti-O-H")
        pd = PhaseDiagram(entries)
        pourbaix_api = Pourbaix_api(mpr) # instantiated in test
        ion_entry_data = pourbaix_api.get_ion_reference_data_for_chemsys("Ti-O-H")
        ion_entries = pourbaix_api.get_ion_entries(pd, ion_entry_data)
        assert len(ion_entries) == 5
        assert all([isinstance(i, IonEntry) for i in ion_entries])
        bi_v_entry_data = pourbaix_api.get_ion_reference_data_for_chemsys("Bi-V")
        bi_data = pourbaix_api.get_ion_reference_data_for_chemsys("Bi")
        v_data = pourbaix_api.get_ion_reference_data_for_chemsys("V")
        assert len(bi_v_entry_data) == len(bi_data) + len(v_data)

        # test an incomplete phase diagram
        entries = mpr.get_entries_in_chemsys("Ti-O")
        pd = PhaseDiagram(entries)
        with pytest.raises(ValueError, match="The phase diagram chemical system"):
            pourbaix_api.get_ion_entries(pd)

        # test ion energy calculation
        ion_data = pourbaix_api.get_ion_reference_data_for_chemsys("S")
        ion_ref_comps = [
            Ion.from_formula(d["data"]["RefSolid"]).composition for d in ion_data
        ]
        ion_ref_elts = set(
            itertools.chain.from_iterable(i.elements for i in ion_ref_comps)
        )
        ion_ref_entries = mpr.get_entries_in_chemsys(
            [*map(str, ion_ref_elts), "O", "H"]
        )
        mpc = MaterialsProjectAqueousCompatibility()
        ion_ref_entries = mpc.process_entries(ion_ref_entries)
        ion_ref_pd = PhaseDiagram(ion_ref_entries)
        ion_entries = pourbaix_api.get_ion_entries(ion_ref_pd, ion_ref_data=ion_data)

        # In ion ref data, SO4-2 is -744.27 kJ/mol; ref solid is -1,279.0 kJ/mol
        # so the ion entry should have an energy (-744.27 +1279) = 534.73 kJ/mol
        # or 5.542 eV/f.u. above the energy of Na2SO4
        so4_two_minus = [e for e in ion_entries if e.ion.reduced_formula == "SO4[-2]"][
            0
        ]

        # the ref solid is Na2SO4, ground state mp-4770
        # the rf factor correction is necessary to make sure the composition
        # of the reference solid is normalized to a single formula unit
        ref_solid_entry = [e for e in ion_ref_entries if e.entry_id == "mp-4770-GGA"][0]
        rf = ref_solid_entry.composition.get_reduced_composition_and_factor()[1]
        solid_energy = ion_ref_pd.get_form_energy(ref_solid_entry) / rf

        assert np.allclose(so4_two_minus.energy, solid_energy + 5.542, atol=1e-3)
    
    @pytest.mark.skip(reason="SSL issues")
    def test_get_pourbaix_entries(self, mpr):
        # test input chemsys as a list of elements
        pourbaix_api = Pourbaix_api(mpr) # instantiated in test
        pbx_entries = pourbaix_api.get_pourbaix_entries(["Fe", "Cr"])
        for pbx_entry in pbx_entries:
            assert isinstance(pbx_entry, PourbaixEntry)

        # test input chemsys as a string
        pbx_entries = pourbaix_api.get_pourbaix_entries("Fe-Cr")
        for pbx_entry in pbx_entries:
            assert isinstance(pbx_entry, PourbaixEntry)

        # test use_gibbs kwarg
        pbx_entries = pourbaix_api.get_pourbaix_entries("Li-O", use_gibbs=300)
        for pbx_entry in pbx_entries:
            assert isinstance(pbx_entry, PourbaixEntry)

        # test solid_compat kwarg
        with pytest.raises(ValueError, match="Solid compatibility can only be"):
            pourbaix_api.get_pourbaix_entries("Ti-O", solid_compat=None)

        # test removal of extra elements from reference solids
        # Li-Zn-S has Na in reference solids
        pbx_entries = pourbaix_api.get_pourbaix_entries("Li-Zn-S")
        assert not any(e for e in pbx_entries if "Na" in e.composition)

        # Ensure entries are pourbaix compatible
        PourbaixDiagram(pbx_entries)

        # TODO - old tests copied from pymatgen with specific energy values. Update or delete
        # fe_two_plus = [e for e in pbx_entries if e.entry_id == "ion-0"][0]
        # self.assertAlmostEqual(fe_two_plus.energy, -1.12369, places=3)
        #
        # feo2 = [e for e in pbx_entries if e.entry_id == "mp-25332"][0]
        # self.assertAlmostEqual(feo2.energy, 3.56356, places=3)
        #
        # # Test S, which has Na in reference solids
        # pbx_entries = self.rester.get_pourbaix_entries(["S"])
        # so4_two_minus = pbx_entries[9]
        # self.assertAlmostEqual(so4_two_minus.energy, 0.301511, places=3)