from itertools import combinations, product

import pytest

import pyEQL
from pyEQL.salt_ion_match import Salt

_CATIONS = ["Na[+1]", "Ca[+2]", "Fe[+3]", "K[+1]", "Li[+1]", "Cu[+2]", "Ba[+2]"]
_ANIONS = ["Cl[-1]", "SO4[-2]", "PO4[-3]", "ClO[-1]", "NO3[-1]", "CO3[-2]", "MnO4[-2]"]
_SALTS = list(product(_CATIONS[:3], _ANIONS[:3]))
_ACIDS = [("H[+1]", anion) for anion in _ANIONS[:3]]
_BASES = [(cation, "OH[-1]") for cation in _CATIONS[:3]]
_SOLUTES = [*_SALTS, *_ACIDS, *_BASES]
_SOLUTES_LITE = [*product(_CATIONS[:2], _ANIONS[1:3]), _ACIDS[0], _BASES[-1]]


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


# Ratio of the concentration of each Salt in the `salts` fixture to that of the previous Salt in the list
# This must be low enough such that the ordering of cation/anion equivalents coincides with the ordering
# of the corresponding Salt objects in the `salt` fixture. Too high a ratio can make it difficult to
# predict which ions will comprise the major and minor salts.
@pytest.fixture(name="salt_ratio", params=[0.25])
def fixture_salt_ratio(request: pytest.FixtureRequest) -> float:
    return float(request.param)


@pytest.fixture(name="salt_conc_units", params=["mol/L"])
def fixture_salt_conc_units(request: pytest.FixtureRequest) -> str:
    return str(request.param)


# This is the default fixture for parametrizing the solution fixture
@pytest.fixture(name="solutes")
def fixture_solutes(
    salts: list[Salt],
    salt_conc: float,
    cation_scale: float,
    anion_scale: float,
    salt_ratio: float,
    salt_conc_units: str,
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
    return {k: f"{v} {salt_conc_units}" for k, v in solute_values.items() if v}


# This is an alternative way to parametrize the solution fixture
# This fixture is preferred if specific pairs of solutes are required (e.g., salts of a conjugate acid/base pair)
@pytest.fixture(name="solute_pairs", params=combinations(_SOLUTES_LITE, r=2))
def fixture_solute_pairs(request: pytest.FixtureRequest) -> tuple[tuple[str, str], tuple[str, str]]:
    solute_pairs: tuple[tuple[str, str], tuple[str, str]] = request.param
    return solute_pairs


@pytest.fixture(name="volume", params=["1 L"])
def fixture_volume(request: pytest.FixtureRequest) -> str:
    return str(request.param)


@pytest.fixture(name="pH", params=[7.0])
def fixture_pH(request: pytest.FixtureRequest) -> float:
    return float(request.param)


@pytest.fixture(name="solution")
def fixture_solution(solutes: dict[str, str], volume: str, pH: float) -> pyEQL.Solution:
    return pyEQL.Solution(solutes=solutes, volume=volume, pH=pH)


@pytest.fixture(name="cutoff", params=[1e-8])
def fixture_cutoff(request: pytest.FixtureRequest) -> float:
    return float(request.param)


@pytest.fixture(name="use_totals", params=[True])
def fixture_use_totals(request: pytest.FixtureRequest) -> bool:
    return bool(request.param)


@pytest.fixture(name="salt_dict")
def fixture_salt_dict(solution: pyEQL.Solution, cutoff: float, use_totals: bool) -> dict[str, dict[str, float | Salt]]:
    salt_dict: dict[str, dict[str, float | Salt]] = solution.get_salt_dict(cutoff=cutoff, use_totals=use_totals)
    return salt_dict
