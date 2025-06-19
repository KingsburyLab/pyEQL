from pathlib import Path

import pytest

from pyEQL.benchmark import BenchmarkEntry, benchmark_engine, calculate_stats, create_entry, get_dataset
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


class TestGetDataset:
    @staticmethod
    @pytest.mark.parametrize("source", ["molar_conductivity.json"])
    def test_should_load_dataset_from_file(source: str) -> None:
        dataset = get_dataset(datadir.joinpath(source))
        assert dataset


@pytest.fixture(name="cation", params=["Na[+1]", "NH4[+1]"])
def fixture_cation(request: pytest.FixtureRequest) -> str:
    return request.param


@pytest.fixture(name="anion", params=["Cl[-1]", "SO4[-2]"])
def fixture_anion(request: pytest.FixtureRequest) -> str:
    return request.param


@pytest.fixture(name="conc", params=["0.1 mol/L", "25%"])
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
    params=["dielectric_constant", "debye_length", "conductivity", "osmotic_coefficient", "density"],
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
    def test_should_benchmark_all_engines(engine: EOS) -> None:
        benchmark_results = benchmark_engine(engine, sources=[datadir.joinpath("molar_conductivity.json")])
        assert benchmark_results
