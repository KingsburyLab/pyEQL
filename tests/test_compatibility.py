from __future__ import annotations

import copy
import math
import os
import sys
from collections import defaultdict
from importlib.resources import files
from typing import TYPE_CHECKING

import pytest
from monty.json import MontyDecoder
from pymatgen.core import Species
from pymatgen.core.composition import Composition
from pymatgen.core.entries import ComputedEntry, ComputedStructureEntry, ConstantEnergyAdjustment
from pymatgen.core.structure import Structure
from pytest import approx

from pyEQL.pourbaix.compatibility import (
    MP2020_COMPAT_CONFIG,
    MP_COMPAT_CONFIG,
    MU_H2O,
    AnionCorrection,
    AqueousCorrection,
    Compatibility,
    CompatibilityError,
    MaterialsProject2020Compatibility,
    MaterialsProjectAqueousCompatibility,
    MaterialsProjectCompatibility,
    needs_u_correction,
)

# from tests.testing import TEST_FILES_DIR

if TYPE_CHECKING:
    from pymatgen.util.typing import CompositionLike

import pymatgen

PMG_ENTRIES_DIR = os.path.dirname(os.path.abspath(pymatgen.analysis.compatibility.__file__))

mp_compat_dir = files("pyEQL") / "pourbaix" / "MITCompatibility.yaml"


@pytest.mark.filterwarnings("ignore:MaterialsProjectCompatibility is deprecated")
# abstract Compatibility tests
class DummyCompatibility(Compatibility):
    """Dummy class to test abstract Compatibility interface."""

    def get_adjustments(self, entry):
        return [ConstantEnergyAdjustment(-10, name="Dummy adjustment")]


def test_process_entries_return_type():
    """process_entries should accept single entries or a list, and always return a list."""
    entry = ComputedEntry("Fe2O3", -2)
    compat = DummyCompatibility()

    assert isinstance(compat.process_entries(entry), list)
    assert isinstance(compat.process_entries([entry]), list)


def test_no_duplicate_corrections():
    """Compatibility should never apply the same correction twice."""
    entry = ComputedEntry("Fe2O3", -2)
    compat = DummyCompatibility()

    assert entry.correction == 0
    compat.process_entries(entry)
    assert entry.correction == -10
    compat.process_entries(entry)
    assert entry.correction == -10
    compat.process_entries(entry, clean=True)
    assert entry.correction == -10


def test_clean_arg():
    """
    clean=False should preserve existing corrections, clean=True should delete
    them before processing.
    """
    entry = ComputedEntry("Fe2O3", -2, correction=-4)
    compat = DummyCompatibility()

    assert entry.correction == -4
    compat.process_entries(entry, clean=False)
    assert entry.correction == -14
    compat.process_entries(entry)
    assert entry.correction == -10


def test_energy_adjustment_normalize():
    """
    Both manual and automatically generated energy adjustments should be scaled
    by the normalize method.
    """
    entry = ComputedEntry("Fe4O6", -2, correction=-4)
    entry = entry.normalize()
    for ea in entry.energy_adjustments:
        if "Manual" in ea.name:
            assert ea.value == -2

    compat = DummyCompatibility()
    entry = ComputedEntry("Fe4O6", -2, correction=-4)
    entry = compat.process_entries(entry)[0]
    entry = entry.normalize()
    for ea in entry.energy_adjustments:
        if "Dummy" in ea.name:
            assert ea.value == -5


def test_overlapping_adjustments():
    """
    Compatibility should raise a CompatibilityError if there is already a
    correction with the same name, but a different value, and process_entries
    should skip that entry.
    """
    ea = ConstantEnergyAdjustment(-5, name="Dummy adjustment")
    entry = ComputedEntry("Fe2O3", -2, energy_adjustments=[ea])
    compat = DummyCompatibility()

    assert entry.correction == -5

    # in case of a collision between EnergyAdjustment, check for a UserWarning
    with pytest.warns(UserWarning, match="already has an energy adjustment called Dummy"):
        processed = compat.process_entries(entry, clean=False)

    assert len(processed) == 0


@pytest.mark.filterwarnings("ignore:MaterialsProjectCompatibility is deprecated")
class TestMaterialsProjectCompatibility:
    def setup_method(self):
        self.entry1 = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",  # codespell:ignore=titel
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",  # codespell:ignore=titel
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )

        self.entry_sulfide = ComputedEntry(
            "FeS",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": False,
                "run_type": "GGA",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",  # codespell:ignore=titel
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE S 08Apr2002",  # codespell:ignore=titel
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )

        self.entry4 = ComputedEntry(
            "H8",
            -27.1,
            correction=0.0,
            parameters={
                "run_type": "LDA",
                "is_hubbard": False,
                "pseudo_potential": {
                    "functional": "PBE",
                    "labels": ["H"],
                    "pot_type": "paw",
                },
                "hubbards": {},
                "potcar_symbols": ["PBE H"],
                "oxide_type": "None",
            },
        )

        self.entry2 = ComputedEntry(
            "Fe3O4",
            -2,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",  # codespell:ignore=titel
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",  # codespell:ignore=titel
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        self.entry3 = ComputedEntry(
            "FeO",
            -2,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 4.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",  # codespell:ignore=titel
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",  # codespell:ignore=titel
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )

        self.compat = MaterialsProjectCompatibility(check_potcar_hash=False)
        self.gga_compat = MaterialsProjectCompatibility("GGA", check_potcar_hash=False)

    def test_process_entry(self):
        # Correct parameters
        assert self.compat.process_entry(self.entry1) is not None
        assert self.gga_compat.process_entry(self.entry1) is None

        # Correct parameters
        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": False,
                "hubbards": {},
                "run_type": "GGA",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",  # codespell:ignore=titel
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",  # codespell:ignore=titel
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is None
        assert self.gga_compat.process_entry(entry) is not None

        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",  # codespell:ignore=titel
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",  # codespell:ignore=titel
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is not None

    def test_correction_values(self):
        # test_corrections
        assert self.compat.process_entry(self.entry1).correction == approx(-2.733 * 2 - 0.70229 * 3)

        entry = ComputedEntry(
            "FeF3",
            -2,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "F": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",  # codespell:ignore=titel
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE F 08Apr2002",  # codespell:ignore=titel
                        "hash": "180141c33d032bfbfff30b3bea9d23dd",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is not None

        # Check actual correction
        assert self.compat.process_entry(entry).correction == approx(-2.733)

        assert self.compat.process_entry(self.entry_sulfide).correction == approx(-0.66346)

    def test_u_values(self):
        # Wrong U value
        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.2, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",  # codespell:ignore=titel
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",  # codespell:ignore=titel
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is None

        # GGA run of U
        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": False,
                "hubbards": None,
                "run_type": "GGA",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",  # codespell:ignore=titel
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",  # codespell:ignore=titel
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is None

        # GGA+U run of non-U
        entry = ComputedEntry(
            "Al2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Al": 5.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Al 06Sep2000",  # codespell:ignore=titel
                        "hash": "805c888bbd2793e462311f6a20d873d9",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",  # codespell:ignore=titel
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is None

        # Materials project should not have a U for sulfides
        entry = ComputedEntry(
            "FeS2",
            -2,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "S": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",  # codespell:ignore=titel
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE S 08Apr2002",  # codespell:ignore=titel
                        "hash": "f7f8e4a74a6cbb8d63e41f4373b54df2",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is None

    def test_wrong_psp(self):
        # Wrong psp
        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe 06Sep2000",  # codespell:ignore=titel
                        "hash": "9530da8244e4dac17580869b4adab115",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",  # codespell:ignore=titel
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is None

    def test_element_processing(self):
        entry = ComputedEntry(
            "O",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": False,
                "hubbards": {},
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE O 08Apr2002",  # codespell:ignore=titel
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    }
                ],
                "run_type": "GGA",
            },
        )
        entry = self.compat.process_entry(entry)
        # assert entry.entry_id == -8
        assert entry.energy == approx(-1)
        assert self.gga_compat.process_entry(entry).energy == approx(-1)

    def test_get_explanation_dict(self):
        compat = MaterialsProjectCompatibility(check_potcar_hash=False)
        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",  # codespell:ignore=titel
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",  # codespell:ignore=titel
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        dct = compat.get_explanation_dict(entry)
        assert dct["corrections"][0]["name"] == "MPRelaxSet Potcar Correction"

        compat.explain(entry)

    def test_get_corrections_dict(self):
        compat = MaterialsProjectCompatibility(check_potcar_hash=False)
        gga_compat = MaterialsProjectCompatibility("GGA", check_potcar_hash=False)

        # Correct parameters
        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",  # codespell:ignore=titel
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",  # codespell:ignore=titel
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        comp = compat.get_corrections_dict(entry)[0]
        assert comp["MP Anion Correction"] == approx(-2.10687)
        assert comp["MP Advanced Correction"] == approx(-5.466)

        entry.parameters["is_hubbard"] = False
        del entry.parameters["hubbards"]
        comp = gga_compat.get_corrections_dict(entry)[0]
        assert "MP Advanced Correction" not in comp

    def test_process_entries(self):
        entries = self.compat.process_entries([self.entry1, self.entry2, self.entry3, self.entry4])
        assert len(entries) == 2

    @pytest.mark.skipif(sys.platform.startswith("win"), reason="Windows broken permissions.")
    def test_parallel_process_entries(self):
        # TODO: get DeprecationWarning: This process (pid=xxxx) is multi-threaded,
        # use of fork() may lead to deadlocks in the child.
        # pid = os.fork()
        with pytest.raises(
            ValueError,
            match="Parallel processing is not possible with for 'inplace=True'",
        ):
            entries = self.compat.process_entries(
                [self.entry1, self.entry2, self.entry3, self.entry4],
                inplace=True,
                n_workers=2,
            )

        entries = self.compat.process_entries(
            [self.entry1, self.entry2, self.entry3, self.entry4],
            inplace=False,
            n_workers=2,
        )
        assert len(entries) == 2

    def test_msonable(self):
        compat_dict = self.compat.as_dict()
        decoder = MontyDecoder()
        temp_compat = decoder.process_decoded(compat_dict)
        assert isinstance(temp_compat, MaterialsProjectCompatibility)


class TestMaterialsProjectCompatibility2020:
    def setup_method(self):
        self.entry1 = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",  # codespell:ignore=titel
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",  # codespell:ignore=titel
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )

        self.entry_sulfide = ComputedEntry(
            "FeS",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": False,
                "run_type": "GGA",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",  # codespell:ignore=titel
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE S 08Apr2002",  # codespell:ignore=titel
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )

        self.entry2 = ComputedEntry(
            "Fe3O4",
            -2,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",  # codespell:ignore=titel
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",  # codespell:ignore=titel
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        self.entry3 = ComputedEntry(
            "FeO",
            -2,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 4.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",  # codespell:ignore=titel
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",  # codespell:ignore=titel
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )

        self.compat = MaterialsProject2020Compatibility(check_potcar_hash=False)
        self.gga_compat = MaterialsProject2020Compatibility("GGA", check_potcar_hash=False)

        self.entry_many_anions = ComputedEntry(
            "CuI9",
            -1,
            0.0,
            parameters={
                "is_hubbard": False,
                "run_type": "GGA",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Cu_pv 06Sep2000",  # codespell:ignore=titel
                        "hash": "2d718b6be91068094207c9e861e11a89",
                    },
                    {
                        "titel": "PAW_PBE I 08Apr2002",  # codespell:ignore=titel
                        "hash": "f4ff16a495dd361ff5824ee61b418bb0",
                    },
                ],
            },
        )

    def test_process_entry(self):
        # Correct parameters
        assert self.compat.process_entry(self.entry1) is not None
        assert self.gga_compat.process_entry(self.entry1) is None

        # Correct parameters
        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": False,
                "hubbards": {},
                "run_type": "GGA",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",  # codespell:ignore=titel
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",  # codespell:ignore=titel
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is None
        assert self.gga_compat.process_entry(entry) is not None

        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",  # codespell:ignore=titel
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",  # codespell:ignore=titel
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is not None

    def test_oxi_state_guess(self):
        # An entry where Composition.oxi_state_guesses will return an empty list
        entry_blank = ComputedEntry(
            "Ga3Te",
            -12.1900,
            correction=0.0,
            parameters={
                "run_type": "GGA",
                "is_hubbard": False,
                "pseudo_potential": {
                    "functional": "PBE",
                    "labels": ["Ga_d", "Te"],
                    "pot_type": "paw",
                },
                "hubbards": {},
                "potcar_symbols": ["PBE Ga_d", "PBE Te"],
                "oxide_type": "None",
            },
        )

        # An entry where one anion will only be corrected if oxidation_states is populated
        entry_oxi = ComputedEntry(
            "Mo2Cl8O",
            -173.0655,
            correction=0.0,
            parameters={
                "run_type": "GGA+U",
                "is_hubbard": True,
                "pseudo_potential": {
                    "functional": "PBE",
                    "labels": ["Mo_pv", "Cl", "O"],
                    "pot_type": "paw",
                },
                "hubbards": {"Mo": 4.38, "Cl": 0.0, "O": 0.0},
                "potcar_symbols": ["PBE Mo_pv", "PBE Cl", "PBE O"],
                "oxide_type": "oxide",
            },
        )

        # An entry that should receive multiple anion corrections if oxidation
        # states are populated
        entry_multi_anion = ComputedEntry(
            "C8N4Cl4",
            -87.69656726,
            correction=0.0,
            parameters={
                "run_type": "GGA",
                "is_hubbard": False,
                "pseudo_potential": {
                    "functional": "PBE",
                    "labels": ["C", "N", "Cl"],
                    "pot_type": "paw",
                },
                "hubbards": {},
                "potcar_symbols": ["PBE C", "PBE N", "PBE Cl"],
                "oxide_type": "None",
            },
        )

        with pytest.warns(UserWarning, match="Failed to guess oxidation state"):
            e1 = self.compat.process_entry(entry_blank)
        assert e1.correction == approx(-0.422)

        e2 = self.compat.process_entry(entry_oxi)
        assert e2.correction == approx(-0.687 + -3.202 * 2 + -0.614 * 8)

        e3 = self.compat.process_entry(entry_multi_anion)
        assert e3.correction == approx(-0.361 * 4 + -0.614 * 4)

    def test_correction_values(self):
        # test_corrections
        assert self.compat.process_entry(self.entry1).correction == approx(-2.256 * 2 - 0.687 * 3)

        entry = ComputedEntry(
            "FeF3",
            -2,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "F": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",  # codespell:ignore=titel
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE F 08Apr2002",  # codespell:ignore=titel
                        "hash": "180141c33d032bfbfff30b3bea9d23dd",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is not None

        # Check actual correction
        assert self.compat.process_entry(entry).correction == approx(-0.462 * 3 + -2.256)

        assert self.compat.process_entry(self.entry_sulfide).correction == approx(-0.503)

    def test_oxidation_by_electronegativity(self):
        # make sure anion corrections are only applied when the element has
        # a negative oxidation state (e.g., correct CaSi but not SiO2 for Si)
        # as determined by electronegativity (i.e., the data.oxidation_states key is absent)

        entry1 = ComputedEntry.from_dict(
            {
                "@module": "pymatgen.entries.computed_entries",
                "@class": "ComputedEntry",
                "energy": -17.01015622,
                "composition": defaultdict(float, {"Si": 2.0, "Ca": 2.0}),
                "energy_adjustments": [],
                "parameters": {
                    "run_type": "GGA",
                    "is_hubbard": False,
                    "pseudo_potential": {
                        "functional": "PBE",
                        "labels": ["Ca_sv", "Si"],
                        "pot_type": "paw",
                    },
                    "hubbards": {},
                    "potcar_symbols": ["PBE Ca_sv", "PBE Si"],
                    "oxide_type": "None",
                },
                "data": {"oxide_type": "None"},
                "entry_id": "mp-1563",
                "correction": 0.0,
            }
        )

        entry2 = ComputedEntry.from_dict(
            {
                "@module": "pymatgen.entries.computed_entries",
                "@class": "ComputedEntry",
                "energy": -47.49120119,
                "composition": defaultdict(float, {"Si": 2.0, "O": 4.0}),
                "energy_adjustments": [],
                "parameters": {
                    "run_type": "GGA",
                    "is_hubbard": False,
                    "pseudo_potential": {
                        "functional": "PBE",
                        "labels": ["Si", "O"],
                        "pot_type": "paw",
                    },
                    "hubbards": {},
                    "potcar_symbols": ["PBE Si", "PBE O"],
                    "oxide_type": "oxide",
                },
                "data": {"oxide_type": "oxide"},
                "entry_id": "mp-546794",
                "correction": 0.0,
            }
        )

        # CaSi; only correction should be Si
        assert self.compat.process_entry(entry1).correction == approx(0.071 * 2)

        # SiO2; only corrections should be oxide
        assert self.compat.process_entry(entry2).correction == approx(-0.687 * 4)

    def test_oxidation(self):
        # make sure anion corrections are only applied when the element has
        # a negative oxidation state (e.g., correct CaSi but not SiO2 for Si)
        # as determined by the data.oxidation_states key

        entry1 = ComputedEntry.from_dict(
            {
                "@module": "pymatgen.entries.computed_entries",
                "@class": "ComputedEntry",
                "energy": -17.01015622,
                "composition": defaultdict(float, {"Si": 2.0, "Ca": 2.0}),
                "energy_adjustments": [],
                "parameters": {
                    "run_type": "GGA",
                    "is_hubbard": False,
                    "pseudo_potential": {
                        "functional": "PBE",
                        "labels": ["Ca_sv", "Si"],
                        "pot_type": "paw",
                    },
                    "hubbards": {},
                    "potcar_symbols": ["PBE Ca_sv", "PBE Si"],
                    "oxide_type": "None",
                },
                "data": {
                    "oxide_type": "None",
                    "oxidation_states": {"Ca": 2.0, "Si": -2.0},
                },
                "entry_id": "mp-1563",
                "correction": 0.0,
            }
        )

        entry2 = ComputedEntry.from_dict(
            {
                "@module": "pymatgen.entries.computed_entries",
                "@class": "ComputedEntry",
                "energy": -47.49120119,
                "composition": defaultdict(float, {"Si": 2.0, "O": 4.0}),
                "energy_adjustments": [],
                "parameters": {
                    "run_type": "GGA",
                    "is_hubbard": False,
                    "pseudo_potential": {
                        "functional": "PBE",
                        "labels": ["Si", "O"],
                        "pot_type": "paw",
                    },
                    "hubbards": {},
                    "potcar_symbols": ["PBE Si", "PBE O"],
                    "oxide_type": "oxide",
                },
                "data": {
                    "oxide_type": "oxide",
                    "oxidation_states": {"Si": 4.0, "O": -2.0},
                },
                "entry_id": "mp-546794",
                "correction": 0.0,
            }
        )

        # CaSi; only correction should be Si
        assert self.compat.process_entry(entry1).correction == approx(0.071 * 2)

        # SiO2; only corrections should be oxide
        assert self.compat.process_entry(entry2).correction == approx(-0.687 * 4)

    def test_u_values(self):
        # Wrong U value
        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.2, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",  # codespell:ignore=titel
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",  # codespell:ignore=titel
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is None

        # GGA run of U
        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": False,
                "hubbards": None,
                "run_type": "GGA",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",  # codespell:ignore=titel
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",  # codespell:ignore=titel
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is None

        # GGA+U run of non-U
        entry = ComputedEntry(
            "Al2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Al": 5.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Al 06Sep2000",  # codespell:ignore=titel
                        "hash": "805c888bbd2793e462311f6a20d873d9",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",  # codespell:ignore=titel
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is None

        # Materials project should not have a U for sulfides
        entry = ComputedEntry(
            "FeS2",
            -2,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "S": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",  # codespell:ignore=titel
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE S 08Apr2002",  # codespell:ignore=titel
                        "hash": "f7f8e4a74a6cbb8d63e41f4373b54df2",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is None

    def test_wrong_psp(self):
        # Wrong psp
        params = {
            "is_hubbard": True,
            "hubbards": {"Fe": 5.3, "O": 0},
            "run_type": "GGA+U",
            "potcar_spec": [
                {
                    "titel": "PAW_PBE Fe 06Sep2000",  # codespell:ignore=titel
                    "hash": "9530da8244e4dac17580869b4adab115",
                },
                {
                    "titel": "PAW_PBE O 08Apr2002",  # codespell:ignore=titel
                    "hash": "7a25bc5b9a5393f46600a4939d357982",
                },
            ],
        }
        entry = ComputedEntry("Fe2O3", -1, correction=0.0, parameters=params)
        assert self.compat.process_entry(entry) is None

    def test_element_processing(self):
        params = {
            "is_hubbard": False,
            "hubbards": {},
            "potcar_spec": [
                {
                    "titel": "PAW_PBE O 08Apr2002",  # codespell:ignore=titel
                    "hash": "7a25bc5b9a5393f46600a4939d357982",
                }
            ],
            "run_type": "GGA",
        }
        entry = ComputedEntry("O", -1, correction=0.0, parameters=params)
        entry = self.compat.process_entry(entry)
        assert entry.energy == approx(-1)
        assert self.gga_compat.process_entry(entry).energy == approx(-1)

    def test_energy_adjustments(self):
        compat = MaterialsProject2020Compatibility(check_potcar_hash=False)
        gga_compat = MaterialsProject2020Compatibility("GGA", check_potcar_hash=False)

        # Fe 4 Co 2 O 8 (Fe2CoO4)
        entry = {
            "@module": "pymatgen.entries.computed_entries",
            "@class": "ComputedEntry",
            "energy": -91.94962744,
            "composition": defaultdict(float, {"Fe": 4.0, "Co": 2.0, "O": 8.0}),
            "energy_adjustments": [],
            "parameters": {
                "run_type": "GGA+U",
                "is_hubbard": True,
                "pseudo_potential": {
                    "functional": "PBE",
                    "labels": ["Fe_pv", "Co", "O"],
                    "pot_type": "paw",
                },
                "hubbards": {"Fe": 5.3, "Co": 3.32, "O": 0.0},
                "potcar_symbols": ["PBE Fe_pv", "PBE Co", "PBE O"],
                "oxide_type": "oxide",
            },
            "data": {"oxide_type": "oxide"},
            "entry_id": "mp-753222",
            "correction": 0,
        }
        entry = ComputedEntry.from_dict(entry)

        processed = compat.process_entry(entry)
        assert "MP2020 anion correction (oxide)" in [ea.name for ea in processed.energy_adjustments]
        assert "MP2020 GGA/GGA+U mixing correction (Fe)" in [ea.name for ea in processed.energy_adjustments]
        assert "MP2020 GGA/GGA+U mixing correction (Co)" in [ea.name for ea in processed.energy_adjustments]

        for ea in processed.energy_adjustments:
            if ea.name == "MP2020 GGA/GGA+U mixing correction (Fe)":
                assert ea.value == approx(-2.256 * 4)
                assert ea.uncertainty == approx(0.0101 * 4)
            elif ea.name == "MP2020 GGA/GGA+U mixing correction (Co)":
                assert ea.value == approx(-1.638 * 2)
                assert ea.uncertainty == approx(0.006 * 2)
            elif ea.name == "MP2020 anion correction (oxide)":
                assert ea.value == approx(-0.687 * 8)
                assert ea.uncertainty == approx(0.002 * 8)

        entry.parameters["is_hubbard"] = False
        del entry.parameters["hubbards"]
        processed = gga_compat.process_entry(entry)
        assert "MP2020 GGA/GGA+U mixing correction" not in [ea.name for ea in processed.energy_adjustments]

        # Fe 4 H 4 O 8 (FeHO2)
        entry = {
            "@module": "pymatgen.entries.computed_entries",
            "@class": "ComputedEntry",
            "energy": -96.71,
            "composition": defaultdict(float, {"Fe": 4.0, "H": 4.0, "O": 8.0}),
            "energy_adjustments": [],
            "parameters": {
                "run_type": "GGA+U",
                "is_hubbard": True,
                "pseudo_potential": {
                    "functional": "PBE",
                    "labels": ["Fe_pv", "H", "O"],
                    "pot_type": "paw",
                },
                "hubbards": {"Fe": 5.3, "H": 0.0, "O": 0.0},
                "potcar_symbols": ["PBE Fe_pv", "PBE H", "PBE O"],
                "oxide_type": "hydroxide",
            },
            "data": {"oxide_type": "hydroxide"},
            "entry_id": "mp-605437",
            "correction": 0,
        }
        entry = ComputedEntry.from_dict(entry)
        adjustments = compat.get_adjustments(entry)

        assert any(ea.name == "MP2020 anion correction (oxide)" for ea in adjustments)
        assert not any(ea.name == "MP2020 anion correction (hydroxide)" for ea in adjustments)

        # Li 4 O 4 (Li2O2)
        entry = {
            "@module": "pymatgen.entries.computed_entries",
            "@class": "ComputedEntry",
            "energy": -38.71,
            "composition": defaultdict(float, {"Li": 4.0, "O": 4.0}),
            "energy_adjustments": [],
            "parameters": {
                "run_type": "GGA",
                "is_hubbard": True,
                "pseudo_potential": {
                    "functional": "PBE",
                    "labels": ["Li", "O"],
                    "pot_type": "paw",
                },
                "hubbards": {"Li": 0.0, "O": 0.0},
                "potcar_symbols": ["PBE Li_sv", "PBE O"],
                "oxide_type": "peroxide",
            },
            "data": {"oxide_type": "peroxide"},
            "entry_id": "mp-841",
            "correction": 0,
        }
        entry = ComputedEntry.from_dict(entry)
        adjustments = compat.get_adjustments(entry)

        assert any(ea.name == "MP2020 anion correction (peroxide)" for ea in adjustments)

        # K 1 O 2 (KO2)
        entry = {
            "@module": "pymatgen.entries.computed_entries",
            "@class": "ComputedEntry",
            "energy": -13.62,
            "composition": defaultdict(float, {"K": 1.0, "O": 2.0}),
            "energy_adjustments": [],
            "parameters": {
                "run_type": "GGA",
                "is_hubbard": True,
                "pseudo_potential": {
                    "functional": "PBE",
                    "labels": ["K", "O"],
                    "pot_type": "paw",
                },
                "hubbards": {"K": 0.0, "O": 0.0},
                "potcar_symbols": ["PBE K_sv", "PBE O"],
                "oxide_type": "superoxide",
            },
            "data": {"oxide_type": "superoxide"},
            "entry_id": "mp-1866",
            "correction": 0,
        }
        entry = ComputedEntry.from_dict(entry)
        adjustments = compat.get_adjustments(entry)

        assert any(ea.name == "MP2020 anion correction (superoxide)" for ea in adjustments)

        # K 4 O 12 (KO3)
        entry = {
            "@module": "pymatgen.entries.computed_entries",
            "@class": "ComputedEntry",
            "energy": -74.14,
            "composition": defaultdict(float, {"K": 4.0, "O": 12.0}),
            "energy_adjustments": [],
            "parameters": {
                "run_type": "GGA",
                "is_hubbard": True,
                "pseudo_potential": {
                    "functional": "PBE",
                    "labels": ["K", "O"],
                    "pot_type": "paw",
                },
                "hubbards": {"K": 0.0, "O": 0.0},
                "potcar_symbols": ["PBE K_sv", "PBE O"],
                "oxide_type": "ozonide",
            },
            "data": {"oxide_type": "ozonide"},
            "entry_id": "mp-1726",
            "correction": 0,
        }
        entry = ComputedEntry.from_dict(entry)
        adjustments = compat.get_adjustments(entry)

        assert any(ea.name == "MP2020 anion correction (ozonide)" for ea in adjustments)

    def test_process_entries(self):
        entries = self.compat.process_entries([self.entry1, self.entry2, self.entry3])
        assert len(entries) == 2

    def test_msonable(self):
        compat_dict = self.compat.as_dict()
        decoder = MontyDecoder()
        temp_compat = decoder.process_decoded(compat_dict)
        assert isinstance(temp_compat, MaterialsProject2020Compatibility)

    def test_check_potcar(self):
        MaterialsProject2020Compatibility(check_potcar=False).process_entries(self.entry1)
        entry = self.entry1.copy()
        del entry.parameters["potcar_spec"]

        # default behavior of checking POTCAR symbols should fail if entry params are missing
        with pytest.raises(KeyError, match="potcar_symbols"):
            MaterialsProject2020Compatibility(check_potcar=True).process_entries(entry)

        # but should work fine if we disable POTCAR checking
        MaterialsProject2020Compatibility(check_potcar=False).process_entries(entry)

    def test_process_entry_with_oxidation_state(self):
        params = {
            "is_hubbard": True,
            "hubbards": {"Fe": 5.3, "O": 0},
            "run_type": "GGA+U",
        }
        entry = ComputedEntry({Species("Fe2+"): 2, Species("O2-"): 3}, -1, parameters=params)

        # Test that MaterialsProject2020Compatibility can process entries with oxidation states
        # https://github.com/materialsproject/pymatgen/issues/3154
        compat = MaterialsProject2020Compatibility(check_potcar=False)
        processed_entry = compat.process_entry(entry, clean=True, inplace=False)

        assert len(processed_entry.energy_adjustments) == 2
        assert processed_entry.energy_adjustments[0].name == "MP2020 anion correction (oxide)"
        assert processed_entry.energy_adjustments[1].name == "MP2020 GGA/GGA+U mixing correction (Fe)"
        assert processed_entry.correction == approx(-6.572999)
        assert processed_entry.energy == approx(-1 + -6.572999)

        # for https://github.com/materialsproject/pymatgen/issues/3425
        frac_coords = [
            [0.5, 0.5, 0.3797505],
            [0.0, 0.0, 0.6202495],
            [0.5, 0.5, 0.8632525],
            [0.0, 0.0, 0.1367475],
            [0.5, 0.0, 0.3608245],
            [0.0, 0.5, 0.0985135],
            [0.5, 0.0, 0.9014865],
            [0.0, 0.5, 0.6391755],
        ]
        lattice = [
            [2.86877900, 0.00000000e00, 1.75662051e-16],
            [-2.83779749e-16, 4.63447500e00, 2.83779749e-16],
            [0.00000000e00, 0.00000000e00, 5.83250700e00],
        ]
        species = ["Li+", "Li+", "Mn3+", "Mn3+", "O2-", "O2-", "O2-", "O2-"]
        li_mn_o = Structure(lattice, species, frac_coords)

        params = {"hubbards": {"Mn": 3.9, "O": 0, "Li": 0}, "run_type": "GGA+U"}
        cse = ComputedStructureEntry(li_mn_o, -58.97, parameters=params)
        processed_entry = compat.process_entry(cse, clean=True, inplace=False)

        assert len(processed_entry.energy_adjustments) == 2
        assert processed_entry.energy_adjustments[0].name == "MP2020 anion correction (oxide)"
        assert processed_entry.energy_adjustments[1].name == "MP2020 GGA/GGA+U mixing correction (Mn)"
        assert processed_entry.correction == approx(-6.084)
        assert processed_entry.energy == approx(-58.97 + -6.084)

    def test_many_anions(self):
        compat = MaterialsProject2020Compatibility(strict_anions="require_bound")
        processed_entry = compat.process_entry(self.entry_many_anions)
        assert processed_entry.energy == -1

        compat = MaterialsProject2020Compatibility(strict_anions="no_check")
        processed_entry = compat.process_entry(self.entry_many_anions)
        assert processed_entry.energy == approx(-4.411)

        compat = MaterialsProject2020Compatibility(strict_anions="require_exact")
        processed_entry = compat.process_entry(self.entry_many_anions)
        assert processed_entry.energy == -1

    def test_custom_config_file(self):
        config_file = str(files("pyEQL") / "pourbaix" / "MP2020Compatibility.yaml")

        compat = MaterialsProject2020Compatibility(
            check_potcar_hash=False,
            config_file=config_file,
        )

        assert compat.config_file == config_file


class TestSulfideTypeCorrection2020:
    def setup_method(self):
        self.compat = MaterialsProject2020Compatibility(check_potcar_hash=False)

    def test_struct_no_struct(self):
        # Processing an Entry should produce the same correction whether or not
        # that entry has a Structure attached to it.

        # Na2S2, entry mp-2400, with and without structure

        entry_struct_as_dict = {
            "@module": "pymatgen.entries.computed_entries",
            "@class": "ComputedStructureEntry",
            "energy": -28.42580746,
            "composition": defaultdict(float, {"Na": 4.0, "S": 4.0}),
            "correction": 0,
            "parameters": {
                "run_type": "GGA",
                "is_hubbard": False,
                "pseudo_potential": {
                    "functional": "PBE",
                    "labels": ["Na_pv", "S"],
                    "pot_type": "paw",
                },
                "hubbards": {},
                "potcar_symbols": ["PBE Na_pv", "PBE S"],
                "oxide_type": "None",
            },
            "data": {"oxide_type": "None"},
            "entry_id": "mp-2400",
            "structure": {
                "@module": "pymatgen.core.structure",
                "@class": "Structure",
                "charge": None,
                "lattice": {
                    "matrix": [
                        [4.5143094, 0.0, 0.0],
                        [-2.2571547, 3.90950662, 0.0],
                        [0.0, 0.0, 10.28414905],
                    ],
                    "a": 4.5143094,
                    "b": 4.514309399183436,
                    "c": 10.28414905,
                    "alpha": 90.0,
                    "beta": 90.0,
                    "gamma": 120.00000000598358,
                    "volume": 181.50209256783256,
                },
                "sites": [
                    {
                        "species": [{"element": "Na", "occu": 1}],
                        "abc": [0.0, 0.0, 0.0],
                        "xyz": [0.0, 0.0, 0.0],
                        "label": "Na",
                        "properties": {"magmom": 0.0},
                    },
                    {
                        "species": [{"element": "Na", "occu": 1}],
                        "abc": [0.0, 0.0, 0.5],
                        "xyz": [0.0, 0.0, 5.142074525],
                        "label": "Na",
                        "properties": {"magmom": 0.0},
                    },
                    {
                        "species": [{"element": "Na", "occu": 1}],
                        "abc": [0.33333333, 0.66666667, 0.25],
                        "xyz": [
                            -2.2571547075855847e-08,
                            2.6063377596983557,
                            2.5710372625,
                        ],
                        "label": "Na",
                        "properties": {"magmom": 0.0},
                    },
                    {
                        "species": [{"element": "Na", "occu": 1}],
                        "abc": [0.66666667, 0.33333333, 0.75],
                        "xyz": [2.2571547225715474, 1.3031688603016447, 7.7131117875],
                        "label": "Na",
                        "properties": {"magmom": 0.0},
                    },
                    {
                        "species": [{"element": "S", "occu": 1}],
                        "abc": [0.33333333, 0.66666667, 0.644551],
                        "xyz": [
                            -2.2571547075855847e-08,
                            2.6063377596983557,
                            6.62865855432655,
                        ],
                        "label": "S",
                        "properties": {"magmom": 0.0},
                    },
                    {
                        "species": [{"element": "S", "occu": 1}],
                        "abc": [0.66666667, 0.33333333, 0.144551],
                        "xyz": [
                            2.2571547225715474,
                            1.3031688603016447,
                            1.4865840293265502,
                        ],
                        "label": "S",
                        "properties": {"magmom": 0.0},
                    },
                    {
                        "species": [{"element": "S", "occu": 1}],
                        "abc": [0.66666667, 0.33333333, 0.355449],
                        "xyz": [
                            2.2571547225715474,
                            1.3031688603016447,
                            3.65549049567345,
                        ],
                        "label": "S",
                        "properties": {"magmom": 0.0},
                    },
                    {
                        "species": [{"element": "S", "occu": 1}],
                        "abc": [0.33333333, 0.66666667, 0.855449],
                        "xyz": [
                            -2.2571547075855847e-08,
                            2.6063377596983557,
                            8.79756502067345,
                        ],
                        "label": "S",
                        "properties": {"magmom": 0.0},
                    },
                ],
            },
        }

        entry_no_struct_as_dict = {
            "@module": "pymatgen.entries.computed_entries",
            "@class": "ComputedEntry",
            "energy": -28.42580746,
            "composition": defaultdict(float, {"Na": 4.0, "S": 4.0}),
            "correction": 0,
            "parameters": {
                "run_type": "GGA",
                "is_hubbard": False,
                "pseudo_potential": {
                    "functional": "PBE",
                    "labels": ["Na_pv", "S"],
                    "pot_type": "paw",
                },
                "hubbards": {},
                "potcar_symbols": ["PBE Na_pv", "PBE S"],
                "oxide_type": "None",
            },
            "data": {"oxide_type": "None"},
            "entry_id": "mp-2400",
        }

        na2s2_entry_struct = ComputedStructureEntry.from_dict(entry_struct_as_dict)
        na2s2_entry_nostruct = ComputedEntry.from_dict(entry_no_struct_as_dict)

        struct_corrected = self.compat.process_entry(na2s2_entry_struct)
        nostruct_corrected = self.compat.process_entry(na2s2_entry_nostruct)

        assert struct_corrected.correction == approx(nostruct_corrected.correction)


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
        hydrate_entry2 = ComputedEntry(Composition("Ca2O2H2"), -10)  # nH2O = 0

        compat.process_entries([hydrate_entry, hydrate_entry2])

        assert hydrate_entry.uncorrected_energy - hydrate_entry.energy == 0
        assert hydrate_entry2.uncorrected_energy - hydrate_entry2.energy == 0

    def test_processing_entries_inplace(self):
        h2o_entry = ComputedEntry(Composition("H2O"), (-5.195 + 0.234) * 3, correction=-0.234 * 3)  # -5.195 eV/atom
        o2_entry = ComputedEntry(Composition("O2"), -4.9276 * 2)
        # check that compatibility scheme does not change input entries themselves
        entries = [h2o_entry, o2_entry]
        entries_copy = copy.deepcopy(entries)
        MaterialsProjectAqueousCompatibility().process_entries(entries, inplace=False)
        assert all(e.correction == e_copy.correction for e, e_copy in zip(entries, entries_copy, strict=True))

    def test_nitrogen_correction(self):
        import numpy as np  # noqa: PLC0415

        compat = MaterialsProjectAqueousCompatibility(
            o2_energy=-10,
            h2o_energy=-20,
            h2o_adjustments=-0.5,
            solid_compat=None,
        )

        entry = ComputedEntry(Composition("Li3N"), -10)
        correction = next(
            adj.value for adj in compat.get_adjustments(entry) if adj.name == "MP Aqueous Nitrogen correction"
        )
        assert np.isclose(correction, 0.26, atol=1e-8)

        entry = ComputedEntry(Composition("N2"), -10)
        assert not any(adj.name == "MP Aqueous Nitrogen correction" for adj in compat.get_adjustments(entry))

    def test_universal_solid_shift_adjustment(self):
        import numpy as np  # noqa: PLC0415

        compat = MaterialsProjectAqueousCompatibility(
            o2_energy=-10,
            h2o_energy=-20,
            h2o_adjustments=-0.5,
            solid_compat=None,
            universal_solid_shift_eV_per_atom=0.1,
        )

        li2o_entry = ComputedEntry(Composition("Li2O"), -10)
        correction = next(
            adj.value
            for adj in compat.get_adjustments(li2o_entry)
            if adj.name == "User universal solid shift (eV/atom)"
        )
        assert np.isclose(correction, 0.3, atol=1e-8)

    # def test_solid_compat_args_propagation(self): # Currently commented out
    #     hydrate_entry = ComputedEntry(Composition("FeH4O2"), -10)

    #     compat = MaterialsProjectAqueousCompatibility(
    #         o2_energy=-10, h2o_energy=-20, h2o_adjustments=-0.5, solid_compat=MaterialsProject2020Compatibility()
    #     )

    #     # the solid compatibility object raises the error
    #     with pytest.raises(CompatibilityError, match="invalid run type"):
    #         entries = compat.process_entries([hydrate_entry], on_error="raise")

    #     entries = compat.process_entries([hydrate_entry], on_error="ignore")
    #     assert len(entries) == 0

    @pytest.mark.skipif(sys.platform.startswith("win"), reason="Windows broken permissions.")
    def test_parallel_process_entries(self):
        hydrate_entry = ComputedEntry(Composition("FeH4O2"), -10)  # nH2O = 2
        hydrate_entry2 = ComputedEntry(Composition("Ca2O2H2"), -10)  # nH2O = 0

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


class TestAqueousCorrection:
    def setup_method(self):
        fp = mp_compat_dir
        self.corr = AqueousCorrection(fp)

    def test_compound_energy(self):
        O2_entry = self.corr.correct_entry(
            ComputedEntry(Composition("O2"), -4.9355 * 2, parameters={"run_type": "GGA"})
        )
        H2_entry = self.corr.correct_entry(ComputedEntry(Composition("H2"), 3, parameters={"run_type": "GGA"}))
        H2O_entry = self.corr.correct_entry(ComputedEntry(Composition("H2O"), 3, parameters={"run_type": "GGA"}))
        H2O_formation_energy = H2O_entry.energy - (H2_entry.energy + O2_entry.energy / 2.0)
        assert H2O_formation_energy == approx(-2.46)

        entry = ComputedEntry(Composition("H2O"), -16, parameters={"run_type": "GGA"})
        entry = self.corr.correct_entry(entry)
        assert entry.energy == approx(-14.916)

        entry = ComputedEntry(Composition("H2O"), -24, parameters={"run_type": "GGA"})
        entry = self.corr.correct_entry(entry)
        assert entry.energy == approx(-14.916)

        entry = ComputedEntry(Composition("Cl"), -24, parameters={"run_type": "GGA"})
        entry = self.corr.correct_entry(entry)
        assert entry.energy == approx(-24.344373)

    def test_correction_str(self):
        assert str(self.corr) == f"{self.corr.name} Aqueous Correction"

    def test_cpd_error_File(self):
        assert isinstance(self.corr.cpd_errors, defaultdict)


class TestCorrectionErrors2020Compatibility:
    def setup_method(self):
        self.compat = MaterialsProject2020Compatibility()

        params = {
            "is_hubbard": True,
            "hubbards": {"Fe": 5.3, "O": 0},
            "run_type": "GGA+U",
            "potcar_spec": [
                {
                    "titel": "PAW_PBE Fe_pv 06Sep2000",  # codespell:ignore=titel
                    "hash": "994537de5c4122b7f1b77fb604476db4",
                },
                {
                    "titel": "PAW_PBE O 08Apr2002",  # codespell:ignore=titel
                    "hash": "7a25bc5b9a5393f46600a4939d357982",
                },
            ],
        }
        self.entry1 = ComputedEntry("Fe2O3", -1, correction=0.0, parameters=params)

        params = {
            "is_hubbard": False,
            "run_type": "GGA",
            "potcar_spec": [
                {
                    "titel": "PAW_PBE Fe_pv 06Sep2000",  # codespell:ignore=titel
                    "hash": "994537de5c4122b7f1b77fb604476db4",
                },
                {
                    "titel": "PAW_PBE S 08Apr2002",  # codespell:ignore=titel
                    "hash": "7a25bc5b9a5393f46600a4939d357982",
                },
            ],
        }
        self.entry_sulfide = ComputedEntry("FeS", -1, 0.0, parameters=params)

        params = {
            "is_hubbard": True,
            "hubbards": {"Fe": 5.3, "O": 0},
            "run_type": "GGA+U",
            "potcar_spec": [
                {
                    "titel": "PAW_PBE Fe_pv 06Sep2000",  # codespell:ignore=titel
                    "hash": "994537de5c4122b7f1b77fb604476db4",
                },
                {
                    "titel": "PAW_PBE O 08Apr2002",  # codespell:ignore=titel
                    "hash": "7a25bc5b9a5393f46600a4939d357982",
                },
            ],
        }
        self.entry2 = ComputedEntry("Fe3O4", -2, correction=0.0, parameters=params)

        params = {
            "is_hubbard": True,
            "hubbards": {"Fe": 5.3, "F": 0},
            "run_type": "GGA+U",
            "potcar_spec": [
                {
                    "titel": "PAW_PBE Fe_pv 06Sep2000",  # codespell:ignore=titel
                    "hash": "994537de5c4122b7f1b77fb604476db4",
                },
                {
                    "titel": "PAW_PBE F 08Apr2002",  # codespell:ignore=titel
                    "hash": "180141c33d032bfbfff30b3bea9d23dd",
                },
            ],
        }
        self.entry_fluoride = ComputedEntry("FeF3", -2, correction=0.0, parameters=params)
        potcar_spec = [
            {
                "titel": "PAW_PBE Li_sv 10Sep2004",  # codespell:ignore=titel
                "hash": "8245d7383d7556214082aa40a887cd96",
            },
            {
                "titel": "PAW_PBE H 15Jun2001",  # codespell:ignore=titel
                "hash": "bb43c666e3d36577264afe07669e9582",
            },
        ]
        self.entry_hydride = ComputedEntry(
            "LiH",
            -2,
            correction=0.0,
            parameters={
                "is_hubbard": False,
                "run_type": "GGA",
                "potcar_spec": potcar_spec,
            },
        )

    def test_errors(self):
        for entry, expected in (
            (self.entry1, math.sqrt((2 * 0.0101) ** 2 + (3 * 0.002) ** 2)),
            (self.entry2, math.sqrt((3 * 0.0101) ** 2 + (4 * 0.002) ** 2)),
            (self.entry_sulfide, 0.0093),
            (self.entry_fluoride, math.sqrt((3 * 0.0026) ** 2 + 0.0101**2)),
            (self.entry_hydride, 0.0013),
        ):
            corrected_entry = self.compat.process_entry(entry)
            assert corrected_entry.correction_uncertainty == approx(expected)


@pytest.mark.parametrize(
    "u_config",
    [
        MP2020_COMPAT_CONFIG["Corrections"]["GGAUMixingCorrections"],
        MP_COMPAT_CONFIG["Advanced"]["UCorrections"],
    ],
)
@pytest.mark.parametrize(
    ("comp", "expected"),
    [
        ("Fe2O3", {"Fe", "O"}),
        ("Fe3O4", {"Fe", "O"}),
        ("FeS", set()),
        ("FeF3", {"Fe", "F"}),
        ("LiH", set()),
        ("H", set()),
        (Composition("MnO"), {"Mn", "O"}),
        (Composition("MnO2"), {"Mn", "O"}),
        (Composition("LiFePO4"), {"Fe", "O"}),
        (Composition("LiFePS4"), set()),
    ],
)
def test_needs_u_correction(comp: CompositionLike, expected: set[str], u_config: dict):
    assert needs_u_correction(comp, u_config=u_config) == expected


def test_explain_with_adjustments(capsys):
    class DummyCompatibility(Compatibility):
        def get_adjustments(self, entry):
            return [ConstantEnergyAdjustment(-10, name="Dummy adjustment")]

    compat = DummyCompatibility()
    entry = compat.process_entry(ComputedEntry("Fe2O3", -10))

    compat.explain(entry)
    out = capsys.readouterr().out

    assert "The uncorrected energy of Fe2 O3" in out
    # since -10 + (-10) = -20
    assert "The final energy after adjustments is -20.000 eV" in out

    entry = ComputedEntry("Fe2O3", -10)

    compat.explain(entry)
    out = capsys.readouterr().out

    assert "No energy adjustments have been applied to this entry." in out


def test_anion_correction():
    import numpy as np  # noqa: PLC0415
    from pymatgen.core import Lattice, Structure  # noqa: PLC0415
    from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry  # noqa: PLC0415

    assert AnionCorrection.__module__ == "pyEQL.pourbaix.compatibility"

    corr_dir = files("pyEQL") / "pourbaix" / "MPCompatibility.yaml"
    correction = AnionCorrection(corr_dir, correct_peroxide=True)

    entry = ComputedEntry("FeS2", 0, data={"sulfide_type": "polysulfide"})
    value = correction.get_correction(entry)
    expected = correction.sulfide_correction["sulfide"] * entry.composition["S"]
    assert np.isclose(value.nominal_value, expected, atol=1e-8)

    entry = ComputedEntry("Li2O2", 0, data={"oxide_type": "peroxide"})
    value = correction.get_correction(entry)
    expected = correction.oxide_correction["peroxide"] * entry.composition["O"]
    assert np.isclose(value.nominal_value, expected, atol=1e-8)

    entry = ComputedEntry("FeHO2", 0, data={"oxide_type": "hydroxide"})
    value = correction.get_correction(entry)
    expected = correction.oxide_correction["oxide"] * entry.composition["O"]
    assert np.isclose(value.nominal_value, expected, atol=1e-8)

    entry = ComputedEntry("Fe2O3", 0)

    with pytest.warns(UserWarning, match="No structure or oxide_type"):
        value = correction.get_correction(entry)

    expected = correction.oxide_correction["oxide"] * entry.composition["O"]
    assert np.isclose(value.nominal_value, expected, atol=1e-8)

    for formula, correction_key in [
        ("Li2O2", "peroxide"),
        ("KO2", "superoxide"),
        ("KO3", "ozonide"),
    ]:
        entry = ComputedEntry(formula, 0)

        value = correction.get_correction(entry)

        expected = correction.oxide_correction[correction_key] * entry.composition["O"]

        assert np.isclose(value.nominal_value, expected, atol=1e-8)

    # Generate a structure to avoid mp-api
    entry = ComputedStructureEntry(
        Structure(
            Lattice.cubic(10),
            ["Mg", "O"],
            [[0, 0, 0], [0.5, 0.5, 0.5]],
        ),
        0,
    )

    value = correction.get_correction(entry)
    expected = correction.oxide_correction["oxide"] * entry.composition["O"]
    assert np.isclose(value.nominal_value, expected, atol=1e-8)
