import itertools
import os
from importlib.resources import files

import numpy as np
import pytest

pytest.importorskip("mp_api", reason="mp_api not installed or incompatible with this Python version")
from mp_api.client import MPRester
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.analysis.pourbaix_diagram import IonEntry, PourbaixDiagram, PourbaixEntry
from pymatgen.core.composition import Composition
from pymatgen.core.ion import Ion
from pymatgen.entries.compatibility import MaterialsProjectAqueousCompatibility
from pymatgen.entries.computed_entries import ComputedEntry

from pyEQL.pourbaix.pourbaix_api import Pourbaix_api


@pytest.fixture
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
        entries = mpr.get_entries_in_chemsys("Ti-O-H", additional_criteria={"thermo_types": ["GGA_GGA+U"]})
        pd = PhaseDiagram(entries)
        pourbaix_api = Pourbaix_api(mpr)  # instantiated in test
        ion_entry_data = pourbaix_api.get_ion_reference_data_for_chemsys("Ti-O-H")
        ion_entries = pourbaix_api.get_ion_entries(pd, ion_entry_data)
        assert len(ion_entries) == 5
        assert all(isinstance(i, IonEntry) for i in ion_entries)
        bi_v_entry_data = pourbaix_api.get_ion_reference_data_for_chemsys("Bi-V")
        bi_data = pourbaix_api.get_ion_reference_data_for_chemsys("Bi")
        v_data = pourbaix_api.get_ion_reference_data_for_chemsys("V")
        assert len(bi_v_entry_data) == len(bi_data) + len(v_data)

        # test an incomplete phase diagram
        entries = mpr.get_entries_in_chemsys("Ti-O", additional_criteria={"thermo_types": ["GGA_GGA+U"]})
        pd = PhaseDiagram(entries)
        with pytest.raises(ValueError, match="The phase diagram chemical system"):
            pourbaix_api.get_ion_entries(pd)

        # test ion energy calculation
        ion_data = pourbaix_api.get_ion_reference_data_for_chemsys("S")
        ion_ref_comps = [Ion.from_formula(d["data"]["RefSolid"]).composition for d in ion_data]
        ion_ref_elts = set(itertools.chain.from_iterable(i.elements for i in ion_ref_comps))
        ion_ref_entries = mpr.get_entries_in_chemsys(
            [*map(str, ion_ref_elts), "O", "H"], additional_criteria={"thermo_types": ["GGA_GGA+U"]}
        )
        mpc = MaterialsProjectAqueousCompatibility()
        ion_ref_entries = mpc.process_entries(ion_ref_entries)
        ion_ref_pd = PhaseDiagram(ion_ref_entries)
        ion_entries = pourbaix_api.get_ion_entries(ion_ref_pd, ion_ref_data=ion_data)

        # In ion ref data, SO4-2 is -744.27 kJ/mol; ref solid is -1,279.0 kJ/mol
        # so the ion entry should have an energy (-744.27 +1279) = 534.73 kJ/mol
        # or 5.542 eV/f.u. above the energy of Na2SO4
        so4_two_minus = next(e for e in ion_entries if e.ion.reduced_formula == "SO4[-2]")

        # the ref solid is Na2SO4, ground state mp-4770
        # the rf factor correction is necessary to make sure the composition
        # of the reference solid is normalized to a single formula unit
        ref_solid_entry = next(e for e in ion_ref_entries if e.entry_id == "mp-4770-GGA")
        rf = ref_solid_entry.composition.get_reduced_composition_and_factor()[1]
        solid_energy = ion_ref_pd.get_form_energy(ref_solid_entry) / rf

        assert np.allclose(so4_two_minus.energy, solid_energy + 5.542, atol=1e-3)

    @pytest.mark.skip(reason="SSL issues")
    def test_get_pourbaix_entries(self, mpr):
        # test input chemsys as a list of elements
        pourbaix_api = Pourbaix_api(mpr)  # instantiated in test
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


def test_get_ion_reference_data_for_chemsys(tmp_path):
    mp_ref_dir = files("pyEQL") / "pourbaix" / "mpr_reference_ion_database.json"

    pourbaix_api = object.__new__(Pourbaix_api)
    pourbaix_api.json_path = mp_ref_dir

    entries = pourbaix_api.get_ion_reference_data_for_chemsys("B-V")
    assert {entry["identifier"] for entry in entries} == {
        "V2O7[4-]",
        "HVO4(aq)",
        "VOH2O2[3+]",
        "H2V10O28[4-]",
        "VO2[+]",
        "VO3[-]",
        "VO4[3-]",
        "HVO4[2-]",
        "HV2O7[3-]",
        "VO[2+]",
        "VO4[-]",
        "H2VO4[+]",
        "H3V2O7[-]",
        "HV10O28[5-]",
        "B4O7[2-]",
        "H2B4O7(aq)",
        "H5(BO3)2(H2O2)2[-]",
        "BH4[-]",
        "BO2[-]",
        "H2BO3[-]",
        "H3BO3(aq)",
        "H2BO3(H2O2)[-]",
        "HB4O7[-]",
        "B(OH)4[-]",
    }
    assert {entry["data"]["MajElements"] for entry in entries} == {"B", "V"}
    assert {entry["data"]["RefSolid"] for entry in entries} == {"B2O3", "VO2"}

    entries = pourbaix_api.get_ion_reference_data_for_chemsys(["Na"])
    assert {entry["identifier"] for entry in entries} == {"Na[+]"}
    assert {entry["formula"] for entry in entries} == {"Na[+1]"}
    assert {entry["data"]["MajElements"] for entry in entries} == {"Na"}
    assert {entry["data"]["RefSolid"] for entry in entries} == {"Na2O"}


def test_get_ion_entries_from_phase_diagram():
    entries = [
        ComputedEntry("H2", 0.0),
        ComputedEntry("O2", 0.0),
        ComputedEntry("Na", 0.0),
        ComputedEntry("S", 0.0),
        ComputedEntry("Na2SO4", -10.0),
    ]
    pd = PhaseDiagram(entries)

    ion_ref_data = [
        {
            "formula": "SO4[-2]",
            "data": {
                "MajElements": "S",
                "RefSolid": "Na2SO4",
                "ΔGᶠRefSolid": {"value": -0.1, "unit": "MJ/mol"},
                "ΔGᶠ": {"value": -0.05, "unit": "MJ/mol"},
            },
        }
    ]

    pourbaix_api = object.__new__(Pourbaix_api)
    ion_entries = pourbaix_api.get_ion_entries(pd, ion_ref_data)

    ref_solid = entries[-1]
    ref_solid_energy = -0.1 / 96485
    ion_free_energy = -0.05 / 96485
    expected_energy = ion_free_energy + (pd.get_form_energy(ref_solid) - ref_solid_energy)

    assert len(ion_entries) == 1
    assert isinstance(ion_entries[0], IonEntry)
    assert ion_entries[0].ion.reduced_formula == "SO4[-2]"
    assert ion_entries[0].energy == pytest.approx(expected_energy)

    ion_ref_data[0]["data"]["ΔGᶠRefSolid"] = {
        "value": -100.0,
        "unit": "kJ/mol",
    }
    ion_ref_data[0]["data"]["ΔGᶠ"] = {
        "value": -50.0,
        "unit": "kJ/mol",
    }

    ion_entries = pourbaix_api.get_ion_entries(pd, ion_ref_data)
    ref_solid_energy = -100.0 / 96.485
    ion_free_energy = -50.0 / 96.485
    expected_energy = ion_free_energy + (pd.get_form_energy(ref_solid) - ref_solid_energy)

    assert ion_entries[0].energy == pytest.approx(expected_energy)


def test_get_pourbaix_entries(monkeypatch):
    from pymatgen.analysis.compatibility import MaterialsProjectAqueousCompatibility  # noqa: PLC0415

    class DummyMPR:
        def get_entries_in_chemsys(self, *args, **kwargs):
            return [
                ComputedEntry(Composition("Fe"), 0),
                ComputedEntry(Composition("O2"), 0),
                ComputedEntry(Composition("H2"), 0),
                ComputedEntry(Composition("Fe2O3"), -100),
            ]

    class DummyGibbsComputedStructureEntry:
        temp = None

        @classmethod
        def from_entries(cls, entries, temp):
            cls.temp = temp
            return entries

    pourbaix_api = Pourbaix_api(DummyMPR())

    monkeypatch.setattr(
        pourbaix_api,
        "get_ion_reference_data_for_chemsys",
        lambda chemsys: [
            {
                "identifier": "Fe[2+]",
                "formula": "Fe[2+]",
                "data": {
                    "MajElements": "Fe",
                    "RefSolid": "Fe",
                    "ΔGᶠ": {"value": 0, "unit": "kJ/mol"},
                    "ΔGᶠRefSolid": {"value": 0, "unit": "kJ/mol"},
                },
            }
        ],
    )

    monkeypatch.setattr(
        pourbaix_api,
        "get_ion_entries",
        lambda pd, ion_ref_data=None: [IonEntry(Ion.from_formula("Fe[2+]"), 0)],
    )

    monkeypatch.setattr(
        MaterialsProjectAqueousCompatibility,
        "process_entries",
        lambda self, entries: entries,
    )

    monkeypatch.setattr(
        "pymatgen.entries.computed_entries.GibbsComputedStructureEntry",
        DummyGibbsComputedStructureEntry,
    )

    pbx_entries = pourbaix_api.get_pourbaix_entries("Fe")

    assert pbx_entries
    assert all(isinstance(e, PourbaixEntry) for e in pbx_entries)

    pbx_entries_with_gibbs = pourbaix_api.get_pourbaix_entries("Fe", use_gibbs=300)
    assert pbx_entries_with_gibbs
    assert all(isinstance(e, PourbaixEntry) for e in pbx_entries_with_gibbs)


def test_nbs_table_ion_data():
    class DummyMPR:
        pass

    pourbaix_api = Pourbaix_api(DummyMPR())
    nbs_db = pourbaix_api.NBS_table_ion_data()

    assert "NO3[-1]" in nbs_db
    assert "SO4[-2]" in nbs_db
    assert "Al[+3]" in nbs_db

    entry = nbs_db["SO4[-2]"]

    assert "exp_form_E" in entry
    assert "exp_entropy" in entry

    assert entry["exp_form_E"]["units"] == "kJ/mol"
    assert entry["exp_entropy"]["units"] == "J/(mol*K)"
