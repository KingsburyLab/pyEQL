"""
pyEQL salt matching test suite
==============================

This file contains tests for the salt-matching algorithm used by pyEQL in
salt_ion_match.py
"""

import logging
import platform
from io import StringIO
from itertools import combinations, product

import numpy as np
import pytest

import pyEQL
from pyEQL.salt_ion_match import Salt


def test_salt_init() -> None:
    s = Salt("Na[+1]", "Cl[-1]")
    assert s.formula == "NaCl"
    assert s.cation == "Na[+1]"
    assert s.anion == "Cl[-1]"
    assert s.nu_anion == 1
    assert s.nu_cation == 1
    assert s.z_cation == 1
    assert s.z_anion == -1

    s = Salt("Fe+3", "OH-1")
    assert s.formula == "Fe(OH)3"
    assert s.cation == "Fe[+3]"
    assert s.anion == "OH[-1]"
    assert s.nu_anion == 3
    assert s.nu_cation == 1
    assert s.z_cation == 3
    assert s.z_anion == -1


def test_empty_solution() -> None:
    """
    test matching a solution that contains no solutes other than water
    """
    s1 = pyEQL.Solution()
    assert isinstance(s1.get_salt(), pyEQL.salt_ion_match.Salt)
    assert s1.get_salt().formula == "HOH"
    assert s1.get_salt().cation == "H[+1]"
    assert s1.get_salt().anion == "OH[-1]"
    assert s1.get_salt().nu_cation == 1
    assert s1.get_salt().nu_anion == 1


def test_single_salt_mono() -> None:
    """
    test matching a solution with a single monovalent salt
    """
    s1 = pyEQL.Solution([["Na+", "2 mol/L"], ["Cl-", "2 mol/L"]])
    assert isinstance(s1.get_salt(), pyEQL.salt_ion_match.Salt)
    assert s1.get_salt().formula == "NaCl"
    assert s1.get_salt().cation == "Na[+1]"
    assert s1.get_salt().anion == "Cl[-1]"
    assert s1.get_salt().nu_cation == 1
    assert s1.get_salt().nu_anion == 1


def test_single_salt_di() -> None:
    """
    test matching a solution with a single divalent salt
    """
    s1 = pyEQL.Solution([["Na+", "4 mol/L"], ["SO4-2", "2 mol/L"]])
    assert isinstance(s1.get_salt(), pyEQL.salt_ion_match.Salt)
    assert s1.get_salt().formula == "Na2SO4"
    assert s1.get_salt().cation == "Na[+1]"
    assert s1.get_salt().anion == "SO4[-2]"
    assert s1.get_salt().nu_cation == 2
    assert s1.get_salt().nu_anion == 1


def test_single_salt_di2() -> None:
    """
    test matching a solution with a single divalent salt
    """
    s1 = pyEQL.Solution([["Fe+3", "1 mol/L"], ["Cl-", "3 mol/L"]])
    assert isinstance(s1.get_salt(), pyEQL.salt_ion_match.Salt)
    assert s1.get_salt().formula == "FeCl3"
    assert s1.get_salt().cation == "Fe[+3]"
    assert s1.get_salt().anion == "Cl[-1]"
    assert s1.get_salt().nu_cation == 1
    assert s1.get_salt().nu_anion == 3


def test_single_ion() -> None:
    """
    test matching a solution containing only a single ion
    """
    s1 = pyEQL.Solution([["Fe+3", "1 mol/L"]])
    assert isinstance(s1.get_salt(), pyEQL.salt_ion_match.Salt)
    assert s1.get_salt().formula == "Fe(OH)3"
    assert s1.get_salt().cation == "Fe[+3]"
    assert s1.get_salt().anion == "OH[-1]"
    assert s1.get_salt().nu_cation == 1
    assert s1.get_salt().nu_anion == 3


@pytest.mark.skipif(platform.machine() == "arm64" and platform.system() == "Darwin", reason="arm64 not supported")
def test_salt_with_equilibration() -> None:
    """
    test matching a solution containing a salt, before and after equilibration.
    Due to speciation changes, the concentration of the salt will decrease unless
    get_salt_dict() uses total concentrations
    """
    s1 = pyEQL.Solution({"Mg+2": "1 mol/L", "Cl-": "2 mol/L"})
    assert isinstance(s1.get_salt(), pyEQL.salt_ion_match.Salt)
    assert s1.get_salt().formula == "MgCl2"
    assert s1.get_salt().cation == "Mg[+2]"
    assert s1.get_salt().anion == "Cl[-1]"
    assert s1.get_salt().nu_cation == 1
    assert s1.get_salt().nu_anion == 2
    assert np.isclose(s1.get_salt_dict()["MgCl2"]["mol"], 1)
    s1.equilibrate()
    assert s1.get_salt().formula == "MgCl2"
    assert np.isclose(s1.get_salt_dict()["MgCl2"]["mol"], 1)


def test_salt_asymmetric() -> None:
    """
    test matching a solution where the cation and anion concentrations
    are not equal
    """
    s1 = pyEQL.Solution([["Na+", "1 mol/kg"], ["Cl-", "4 mol/kg"]])
    assert isinstance(s1.get_salt(), pyEQL.salt_ion_match.Salt)
    assert s1.get_salt().formula == "NaCl"
    assert s1.get_salt().cation == "Na[+1]"
    assert s1.get_salt().anion == "Cl[-1]"
    assert s1.get_salt().nu_cation == 1
    assert s1.get_salt().nu_anion == 1


"""
Solution.get_salt_dict pytest.fixture reference
===============================================

This section summarizes the organization of pytest fixtures for the Solution.get_salt_dict unit tests.

Generally, unit tests request and inspect the `salt_dict` and `salts` fixtures. The `salt_dict`
fixture returns the result of Solution.get_salt_dict called with values, `cutoff` and `use_totals`. The Solution object
on which Solution.get_salt_dict is called along with the values of `cutoff` and `use_totals` are all set by the
`solution`, `cutoff` and `use_totals` fixtures, respectively. In particular, the `solution` fixture is configured
from fixtures corresponding to keyword arguments to the Solution constructor. (Currently, only "solutes" and "volume"
are needed.) The composition of `solution` is ultimately determined by the `salts` fixture (a list of Salt objects).
Each Salt in `salts` is added to `solution` in amounts determined by the floats corresponding to the `salt_conc`,
`cation_scale`, `anion_scale`, and `salt_ratio` fixtures.

- `salt_conc`: determines the *base* concentration of salts; if `cation_scale = anion_scale = salt_ratio = 1.0`, then
  the concentration of each ion will be that resulting from dissolving a concentration of the salt equal to `salt_conc`.
- `cation_scale`/`anion_scale`: a factor with which to scale the concentration of the cation/anion in a salt relative
  to `salt_conc`. This is applied to each pair of ions in `salts`. If `cation_scale = 0.0`/`anion_scale = 0.0`, then
  only the anions/cations are added to `solution`.
- `salt_ratio`: This is the ratio of the concentration of any given Salt in `salts` relative to the previous Salt in
  the list. Note that if this value is 0.0, then no more than one Salt is added to `solution`.

The Salt objects in `salts` are constructed from salt-ion pairs in `_SOLUTES` or the Cartesian product of `_CATIONS` and
`_CONJUGATE_BASES`.

Under different parametrizations of the compositions fixtures, several Solution test cases are covered including:
- an empty solution
- stochiometric/unequal ratios of single/double salt-ion pairs in a solution, including the case where the same ion
  belongs to multiple salts (e.g., Solution(solutes={'Na[+1]': '3 mol/L', 'Cl[-1]': '1 mol/L', 'SO4[-2]': '1 mol/L'}))
- single/double anion/cation only solutions (e.g., a solution of Na+ and/or Ca[+2] ions)
- salts composed of mono-/di-/trivalent ions
"""


_CATIONS = ["Na[+1]", "Ca[+2]", "Fe[+3]"]
_ANIONS = ["Cl[-1]", "SO4[-2]", "PO4[-3]"]
_SALTS = list(product(_CATIONS, _ANIONS))
_ACIDS = [("H[+1]", anion) for anion in _ANIONS]
_BASES = [(cation, "OH[-1]") for cation in _CATIONS]
_CONJUGATE_BASES = [
    ("HCO3[-1]", "CO3[-2]"),
    ("H2PO4[-1]", "PO4[-3]"),
]
_SOLUTES = [*_SALTS, *_ACIDS, *_BASES]


@pytest.fixture(name="salt", params=_SOLUTES)
def fixture_salt(request: pytest.FixtureRequest) -> Salt:
    cation, anion = request.param
    return Salt(cation=cation, anion=anion)


@pytest.fixture(name="salts")
def fixture_salts(salt: Salt) -> list[Salt]:
    return [salt]


# When cation_scale = anion_scale = 1.0, this is equal to the concentration of the first salt in salts
@pytest.fixture(name="salt_conc", params=[1.0])
def fixture_salt_conc(request: pytest.FixtureRequest) -> float:
    return float(request.param)


# Used to scale cation concentration relative to anion concentration
@pytest.fixture(name="cation_scale", params=[1.0])
def fixture_cation_scale(request: pytest.FixtureRequest) -> float:
    return float(request.param)


# Used to scale anion concentration relative to cation concentration
@pytest.fixture(name="anion_scale", params=[1.0])
def fixture_anion_scale(request: pytest.FixtureRequest) -> float:
    return float(request.param)


# Ratio of the concentration of each Salt in salts to that of the previous Salt in the list
# This must be low enough such that when ordered according to concentration, anionic and cationic solutes of the same
# index correspond to a Salt in `salts` (e.g., consider a solution of 1 M NaCl and 0.75 M K2SO4, corresponding to
# salt_conc = cation_scale = anion_scale = 1.0 and salt_ratio = 0.75. Calling .get_salt_dict on such a solution should
# return entries for 1 M KCl, 0.25 M K2SO4, and 0.5 M Na2SO4, which would cause
# test_should_calculate_correct_concentration_for_salts to fail. To fix this, salt_ratio must be less than 0.5.)
@pytest.fixture(name="salt_ratio", params=[0.1])
def fixture_salt_ratio(request: pytest.FixtureRequest) -> float:
    return float(request.param)


@pytest.fixture(name="solutes")
def fixture_solutes(
    salts: list[Salt], salt_conc: float, cation_scale: float, anion_scale: float, salt_ratio: float
) -> dict[str, str]:
    solute_values = {salt.anion: 0 for salt in salts}
    solute_values.update({salt.cation: 0 for salt in salts})

    for i, salt in enumerate(salts):
        # Scale salt component concentrations
        cation_conc = salt_conc * salt.nu_cation * cation_scale * (salt_ratio**i)
        anion_conc = salt_conc * salt.nu_anion * anion_scale * (salt_ratio**i)
        # Increase solute concentrations
        solute_values[salt.cation] += cation_conc
        solute_values[salt.anion] += anion_conc

    # Only include solutes with non-zero concentrations
    return {k: f"{v} mol/L" for k, v in solute_values.items() if v}


@pytest.fixture(name="volume", params=["1 L"])
def fixture_volume(request: pytest.FixtureRequest) -> str:
    return str(request.param)


@pytest.fixture(name="solution")
def fixture_solution(solutes: dict[str, str], volume: str) -> pyEQL.Solution:
    return pyEQL.Solution(solutes=solutes, volume=volume)


@pytest.fixture(name="cutoff", params=[0.01])
def fixture_cutoff(request: pytest.FixtureRequest) -> float:
    return float(request.param)


@pytest.fixture(name="use_totals", params=[True])
def fixture_use_totals(request: pytest.FixtureRequest) -> bool:
    return bool(request.param)


@pytest.fixture(name="salt_dict")
def fixture_salt_dict(solution: pyEQL.Solution, cutoff: float, use_totals: bool) -> dict[str, dict[str, float | Salt]]:
    salt_dict: dict[str, dict[str, float | Salt]] = solution.get_salt_dict(cutoff=cutoff, use_totals=use_totals)
    return salt_dict


@pytest.mark.xfail(reason="Undecided on ignoring 'HOH' in Solution.get_salt_dict", strict=False)
@pytest.mark.parametrize("salts", [[]])
def test_should_return_empty_dict_for_empty_solution(salt_dict: dict[str, dict[str, float | Salt]]) -> None:
    assert not salt_dict


class TestSaltDictTypes:
    @staticmethod
    @pytest.mark.parametrize("cation_scale", [0.5, 1.0, 1.5])
    def test_should_store_mol_as_floats(salt_dict: dict[str, dict[str, float | Salt]]) -> None:
        mol_values = [d["mol"] for d in salt_dict.values()]
        mol_values_are_floats = [isinstance(mol, float) for mol in mol_values]
        assert all(mol_values_are_floats)

    @staticmethod
    @pytest.mark.parametrize("cation_scale", [0.5, 1.0, 1.5])
    def test_should_store_salt_as_salts(salt_dict: dict[str, dict[str, float | Salt]]) -> None:
        salt_values = [d["salt"] for d in salt_dict.values()]
        salt_values_are_salts = [isinstance(salt, Salt) for salt in salt_values]
        assert all(salt_values_are_salts)


class TestGetSaltDict:
    @staticmethod
    def test_should_match_equimolar_ion_equivalents(
        salts: list[Salt], salt_dict: dict[str, dict[str, float | Salt]]
    ) -> None:
        salts_in_dict = [salt.formula in salt_dict for salt in salts]
        assert all(salts_in_dict)

    @staticmethod
    @pytest.mark.parametrize("cation_scale", [1.1])
    def test_should_match_salts_with_excess_cation(
        salts: list[Salt], salt_dict: dict[str, dict[str, float | Salt]]
    ) -> None:
        salts_in_dict = [salt.formula in salt_dict for salt in salts]
        assert all(salts_in_dict)

    @staticmethod
    @pytest.mark.parametrize("anion_scale", [1.1])
    def test_should_match_salts_with_excess_anion(
        salts: list[Salt], salt_dict: dict[str, dict[str, float | Salt]]
    ) -> None:
        salts_in_dict = [salt.formula in salt_dict for salt in salts]
        assert all(salts_in_dict)

    @staticmethod
    # This parametrization ensures that the solution contains only a single anion other than H+
    @pytest.mark.parametrize(("cation_scale", "salt_ratio"), [(0.0, 0.0)])
    def test_should_match_excess_anions_with_protons(
        salts: list[Salt], salt_dict: dict[str, dict[str, float | Salt]]
    ) -> None:
        acid = Salt("H[+1]", salts[0].anion)
        assert acid.formula in salt_dict

    @staticmethod
    # This parametrization ensures that the solution contains only a single cation other than OH-
    @pytest.mark.parametrize(("anion_scale", "salt_ratio"), [(0.0, 0.0)])
    def test_should_match_excess_cations_with_hydroxide(
        salts: list[Salt], salt_dict: dict[str, dict[str, float | Salt]]
    ) -> None:
        base = Salt(salts[0].cation, "OH[-1]")
        assert base.formula in salt_dict

    @staticmethod
    @pytest.mark.xfail(reason="Undecided on ignoring 'HOH' in Solution.get_salt_dict", strict=False)
    def test_should_not_include_water(salt_dict: dict[str, dict[str, float | Salt]]) -> None:
        assert "HOH" not in salt_dict

    @staticmethod
    # This cutoff should be higher than the concentration of at least one salt (see salt_conc, anion/cation_scale,
    # salt_ratio fixtures)
    # Parametrization over volume is required to ensure that concentrations and not absolute molar quantities are
    # used as a cutoff
    @pytest.mark.parametrize(("cutoff", "volume"), product([0.75, 1.5], ["1 L", "2 L"]))
    def test_should_not_return_salts_with_concentration_below_cutoff(
        salt_dict: dict[str, dict[str, float | Salt]], cutoff: float
    ) -> None:
        salt_concentrations_below_cutoff = [d["mol"] < cutoff for d in salt_dict.values()]
        assert all(salt_concentrations_below_cutoff)

    @staticmethod
    # This parametrization ensures that the concentration of the salt is higher than that of water (~55 M)
    @pytest.mark.parametrize(
        ("anion_scale", "salt_ratio", "salt_conc", "salts"), [(0.0, 0.0, 100.0, [Salt("Na", "Cl")])]
    )
    def test_should_log_warning_for_high_concentrations(solution: pyEQL.Solution) -> None:
        stream = StringIO()
        sh = logging.StreamHandler(stream)
        solution.logger.addHandler(sh)
        solution.logger.setLevel("WARNING")
        _ = solution.get_salt_dict()
        stream.seek(0)
        msg = stream.read()
        assert "H2O(aq) is not the most prominent component in this Solution!" in msg

    @staticmethod
    def test_should_order_salts_by_amount(salt_dict: dict[str, dict[str, float | Salt]]) -> None:
        salt_amounts = [d["mol"] for d in salt_dict.values()]
        assert salt_amounts == sorted(salt_amounts)

    @staticmethod
    def test_should_calculate_correct_concentration_for_salts(
        salt_dict: dict[str, dict[str, float | Salt]], salts: list[Salt], salt_conc: float, salt_ratio: float
    ) -> None:
        salts_have_correct_concentrations = []
        expected_concentrations = dict.fromkeys([s.formula for s in salts], 0.0)

        for i, salt in enumerate(salts):
            expected_concentrations[salt.formula] += salt_conc * (salt_ratio**i)

        for salt, expected in expected_concentrations.items():
            calculated = salt_dict[salt]["mol"]
            salts_have_correct_concentrations.append(calculated == expected)

        assert all(salts_have_correct_concentrations)


class TestGetSaltDictMultipleSalts(TestGetSaltDict):
    @staticmethod
    @pytest.fixture(name="solute_pairs", params=combinations(_SOLUTES, r=2))
    def fixture_solute_pairs(request: pytest.FixtureRequest) -> tuple[tuple[str, str], tuple[str, str]]:
        solute_pairs: tuple[tuple[str, str], tuple[str, str]] = request.param
        return solute_pairs

    @staticmethod
    @pytest.fixture(name="salts")
    def fixture_salts(solute_pairs: tuple[tuple[str, str], tuple[str, str]]) -> list[Salt]:
        major_salt_cation, major_salt_anion = solute_pairs[0]
        minor_salt_cation, minor_salt_anion = solute_pairs[1]
        major_salt = Salt(major_salt_cation, major_salt_anion)
        minor_salt = Salt(minor_salt_cation, minor_salt_anion)

        return [major_salt, minor_salt]

    @staticmethod
    @pytest.mark.parametrize("cation_scale", [2.0])
    def test_should_match_salts_with_excess_cation_if_cation_enough_for_both_anions(
        major_salt: Salt, minor_salt: Salt, salt_dict: dict[str, dict[str, float | Salt]]
    ) -> None:
        mixed_salt = Salt(major_salt.cation, minor_salt.anion)
        assert major_salt.formula in salt_dict
        assert mixed_salt.formula in salt_dict

    @staticmethod
    @pytest.mark.parametrize("anion_scale", [2.0])
    def test_should_match_salts_with_excess_anion_if_anion_enough_for_both_cations(
        major_salt: Salt, minor_salt: Salt, salt_dict: dict[str, dict[str, float | Salt]]
    ) -> None:
        mixed_salt = Salt(minor_salt.cation, major_salt.anion)
        assert major_salt.formula in salt_dict
        assert mixed_salt.formula in salt_dict

    @staticmethod
    @pytest.mark.parametrize("use_totals", [False])
    @pytest.mark.parametrize(
        ("solute_pairs"),
        product(_CATIONS, _CONJUGATE_BASES),
    )
    def test_should_include_salt_for_low_concentration_conjugate_base_when_use_totals_false(
        salt_dict: dict[str, dict[str, float | Salt]], salts: list[Salt]
    ) -> None:
        salts_in_dict = [salt.formula in salt_dict for salt in salts]
        assert all(salts_in_dict)

    @staticmethod
    @pytest.mark.parametrize(
        ("solute_pairs"),
        product(_CATIONS, _CONJUGATE_BASES),
    )
    def test_should_not_include_salt_for_low_concentration_conjugate_base_when_use_totals_true(
        salt_dict: dict[str, dict[str, float | Salt]],
        major_salt: Salt,
        minor_salt: Salt,
    ) -> None:
        assert major_salt.formula in salt_dict
        assert minor_salt.formula not in salt_dict
