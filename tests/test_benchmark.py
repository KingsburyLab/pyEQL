from pathlib import Path

import pytest

from pyEQL.benchmark import BenchmarkEntry, benchmark_engine, calculate_stats, create_entry, load_dataset
from pyEQL.engines import EOS, IdealEOS, NativeEOS
from pyEQL.solution import Solution

datadir = Path(__file__).parent.joinpath(Path(__file__).stem)


@pytest.fixture(
    name="engine",
    params=[
        IdealEOS,
        NativeEOS,
        # PhreeqcEOS - hide for now
    ],
)
def fixture_engine(request: pytest.FixtureRequest) -> EOS:
    engine: type[EOS] = request.param
    return engine()


@pytest.fixture(name="source", params=["molar_conductivity.json", "mean_activity_coefficient.json"])
def fixture_source(request: pytest.FixtureRequest) -> list[Path]:
    return datadir.joinpath(request.param)


class TestLoadDataset:
    @staticmethod
    def test_should_load_dataset_from_file(source: str) -> None:
        dataset = load_dataset(datadir.joinpath(source))
        assert dataset

    @staticmethod
    @pytest.mark.parametrize("source", ["CRC"])
    def test_should_load_dataset_from_internal_source(source: str) -> None:
        dataset = load_dataset(source)
        assert dataset


@pytest.fixture(name="cation", params=["Na[+1]"])
def fixture_cation(request: pytest.FixtureRequest) -> str:
    return request.param


@pytest.fixture(name="anion", params=["Cl[-1]"])
def fixture_anion(request: pytest.FixtureRequest) -> str:
    return request.param


@pytest.fixture(name="conc", params=["0.1 mol/L"])
def fixture_conc(request: pytest.FixtureRequest) -> str:
    return request.param


@pytest.fixture(
    name="solutes",
)
def fixture_solutes(cation: str, anion: str, conc: str) -> dict[str, str]:
    return {cation: conc, anion: conc}


@pytest.fixture(name="solution")
def fixture_solution(solutes: dict[str, str], engine: EOS) -> Solution:
    return Solution(solutes=solutes, engine=engine)


@pytest.fixture(name="solute_property", params=["activity_coefficient", "molar_conductivity"])
def fixture_solute_property(request: pytest.FixtureRequest) -> str:
    return request.param


@pytest.fixture(
    name="solution_property",
    params=[
        "dielectric_constant",
        "debye_length",
        "conductivity",
        "osmotic_coefficient",
        "density",
        "viscosity_dynamic",
        "osmolality",
        "osmolarity",
        "chemical_potential_energy",
        "activity_coefficient_pitzer",
    ],
)
def fixture_solution_property(request: pytest.FixtureRequest) -> str:
    return request.param


@pytest.fixture(
    name="dataset",
)
def fixture_dataset(solution: Solution, solute_property: str, solution_property: str) -> list[BenchmarkEntry]:
    return [create_entry(solution, [solute_property], [solution_property])]


class TestCalculateStats:
    @staticmethod
    def test_should_calculate_zero_rmse_with_engine_dataset(dataset: list[BenchmarkEntry]) -> None:
        results = calculate_stats(dataset)
        assert all(err == 0 for err in results.solute_stats.values())
        assert all(err == 0 for err in results.solution_stats.values())


class TestBenchmarkEngine:
    @staticmethod
    def test_should_benchmark_all_engines(engine: EOS, source: Path) -> None:
        benchmark_results = benchmark_engine(engine, sources=[source])
        assert benchmark_results
