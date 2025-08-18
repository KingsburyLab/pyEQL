from itertools import combinations, product
from typing import Literal

import pytest

import pyEQL
from pyEQL.engines import EOS
from pyEQL.salt_ion_match import Salt

_CATIONS = ["Na[+1]", "Ca[+2]", "Fe[+3]", "K[+1]", "Li[+1]", "Cu[+2]", "Ba[+2]", "NH4[+1]"]
_ANIONS = ["Fe(CN)6[-3]", "Cl[-1]", "SO4[-2]", "PO4[-3]", "ClO[-1]", "NO3[-1]", "CO3[-2]", "MnO4[-2]"]
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
@pytest.fixture(
    name="solute_pairs",
    params=[
        (s1, s2)
        for s1, s2 in combinations(_SOLUTES_LITE, r=2)
        if "H[+1]" not in (s1[0], s2[0]) and "OH[-1]" not in (s1[1] and s2[1])
    ],
)
def fixture_solute_pairs(request: pytest.FixtureRequest) -> tuple[tuple[str, str], tuple[str, str]]:
    solute_pairs: tuple[tuple[str, str], tuple[str, str]] = request.param
    return solute_pairs


@pytest.fixture(name="volume", params=["1 L"])
def fixture_volume(request: pytest.FixtureRequest) -> str:
    return str(request.param)


@pytest.fixture(name="pH", params=[7.0])
def fixture_pH(request: pytest.FixtureRequest) -> float:
    return float(request.param)


@pytest.fixture(name="solvent", params=["H2O"])
def fixture_solvent(request: pytest.FixtureRequest) -> str:
    return str(request.param)


@pytest.fixture(name="engine", params=["native"])
def fixture_engine(request: pytest.FixtureRequest) -> str:
    return str(request.param)


@pytest.fixture(name="database", params=[None])
def fixture_database(request: pytest.FixtureRequest) -> str | None:
    return request.param if request.param is None else str(request.param)


@pytest.fixture(name="solution")
def fixture_solution(
    solutes: dict[str, str],
    volume: str,
    pH: float,
    solvent: str | list,
    engine: EOS | Literal["native", "ideal", "phreeqc"],
    database: str | None,
) -> pyEQL.Solution:
    return pyEQL.Solution(solutes=solutes, volume=volume, pH=pH, solvent=solvent, engine=engine, database=database)


# Model Parameters


# Pitzer activity/osmotic parameters
@pytest.fixture(name="alphas")
def fixture_alphas(salt: Salt) -> tuple[float, float]:
    if salt.z_cation >= 2 and salt.z_anion <= -2:
        if salt.z_cation >= 3 or salt.z_anion <= -3:
            alpha1 = 2.0
            alpha2 = 50.0
        else:
            alpha1 = 1.4
            alpha2 = 12
    else:
        alpha1 = 2.0
        alpha2 = 0.0

    return alpha1, alpha2
