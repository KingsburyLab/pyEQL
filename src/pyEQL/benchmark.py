"""Solution model benchmarking utilities.


Usage:

    python pyeql.benchmark

"""

import json
from collections.abc import Callable
from functools import reduce
from pathlib import Path
from typing import Any, Literal, NamedTuple

import numpy as np

from pyEQL.engines import EOS, IdealEOS, NativeEOS, PhreeqcEOS
from pyEQL.salt_ion_match import Salt
from pyEQL.solution import Solution
from pyEQL.utils import FormulaDict

# TODO: Select and validate data sources
# If all solution reference data are generated from the same solutions, then solutions can be used as an input into
# source creation
SOURCES: list[str] = ["CRC", "IDST", "May2011JCED"]


class _BenchmarkEntry(NamedTuple):
    solution: Solution
    # dict[str, list[tuple[str, float]]]: solute, [(property, value)]
    solute_data: FormulaDict = FormulaDict()
    solution_data: list[tuple[str, float]] = []


# TODO: check tests for missing property checks and values
# TODO: write sub-loaders for different CRC archive types
# TODO: consolidate data files
# TODO: identify/understand theoretical, reference, and package equivalences
# TODO: check tests for other reference databases
# TODO: find other reference databases
# TODO: write loading function for other reference databases
def _load_crc_data(s) -> list[_BenchmarkEntry]:
    # datasets = []

    # solute data
    # load activity coefficient data

    # load molal electrical conductivity

    # load diffusion coefficient data

    # merge solution.solute_data

    # solution data
    # load density data

    # load electrical conductivity

    # ? load osmotic coefficient data

    # merge solution.solution_data

    # previous parametrization of ion pairs
    # cations = [("H+", 1), ("Cs+", 1), ("Li+", 1), ("Rb+", 1), ("K+", 1), ("Na", 1), ("Mg", 2), ("Ba", 2)]
    # anions = [("Cl-", 1), ("I-", 1), ("Br", 1), ("SO4-2", 2)]
    # previous list of concentrations to test, mol/kg
    # conc_list = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5]

    return []


def _get_dataset(source: str) -> list[_BenchmarkEntry]:
    """Load reference dataset.

    Args:
        source: One of "CRC", "IDST", or "May2011JCED" or the path to a file containing reference data. If the latter,
            then the [path must point to a JSON which can be read into a _BenchmarkEntry object.

    Returns:
        A list of _BenchmarkEntry objects one for each data point in the data set.
    """
    match source:
        case "CRC":
            reference = _load_crc_data(source)
        case "IDST":
            pass
        case "May2011JCED":
            pass
        case _:
            with Path(source).open(mode="r", encoding="utf-8") as file:
                data = json.load(file)

            reference: list[_BenchmarkEntry] = []

            for solution, solute_data, solution_data in data:
                reference.append(
                    _BenchmarkEntry(solution=solution, solute_data=solute_data, solution_data=solution_data)
                )

    return reference


def _patch_dataset(
    dataset: list[_BenchmarkEntry], *, engine: EOS | Literal["native", "ideal", "phreeqc"] = "native"
) -> None:
    for data in dataset:
        data.solution.engine = engine


def _rmse(data: list[tuple[float, float]]) -> float:
    return np.std([ref - calc for ref, calc in data])


def _get_solute_property(solution: Solution, solute: str, name: str) -> Any | None:
    return solution.get_property(solute, name)


def _get_mean_activity(solution: Solution) -> float:
    activity_nu_pairs: list[tuple[float, int]] = []

    for salt_dict in solution.get_salt_dict().values():
        salt = Salt.from_dict(salt_dict)
        act_cat = solution.get_activity_coefficient(salt.cation)
        act_an = solution.get_activity_coefficient(salt.anion)
        activity_nu_pairs.extend([(act_an, salt.nu_anion), (act_cat, salt.nu_cation)])

    factor = reduce(lambda x, y: x * y[0] ** y[1], activity_nu_pairs, initial=1.0)
    exponent = 1 / sum(x[1] for x in activity_nu_pairs)
    return factor**exponent


def _get_solution_property(solution: Solution, name: str) -> Any | None:
    if name == "mean_activity":
        return _get_mean_activity(solution)
    if hasattr(solution, name):
        return getattr(solution, name)

    if hasattr(solution, f"get_{name}"):
        return getattr(solution, f"get_{name}")

    msg = f"Property {name} is not supported"
    raise NotImplementedError(msg)


def report_results(
    dataset: list[_BenchmarkEntry], *, metric: Callable[[list[tuple[float, float]]], float] | None = None
) -> None:
    """Report the results of the benchmarking.

    Args:
        dataset: A dictionary mapping 2-tuples (engine, source) to the associated root-means-squared error across
            all data points in the source for the given engine.
        metric: A function that acts on the list of 2-tuples (reference, calculated), which contains reference and
            calculated values. This function should calculate a statistical metric for the list. Defaults to the root-
            mean-squared error.
    """
    metric = metric or _rmse

    # Populate data structure for tracking activity/osmotic coefficient and solvent volume statistics
    # property, (reference, calculated)
    solute_data_pairs: dict[str, list[tuple[float, float]]] = {}
    # property, (reference, calculated)
    solution_data_pairs: dict[str, list[tuple[float, float]]] = {}

    for d in dataset:
        for solute, solute_data in d.solute_data.items():
            for property, reference in solute_data:
                if property not in solute_data_pairs:
                    solute_data_pairs[property] = []

                solute_data_pairs[property].append((reference, _get_solute_property(d.solution, solute, property)))

        for property, reference in d.solution_data:
            if property not in solution_data_pairs:
                solution_data_pairs[property] = []

            solution_data_pairs[property].append((reference, _get_solution_property(d.solution, solute, property)))

    solute_stats = {k: metric(v) for k, v in solute_data_pairs.items()}
    solution_stats = {k: metric(v) for k, v in solution_data_pairs.items()}

    for property, stat in solute_stats.items():
        print(f"{property} statistics: {stat}")

    for property, stat in solution_stats.items():
        print(f"{property} statistics: {stat}")


def main() -> None:
    """Run solution benchmarking logic.

    This function works by reading in reference property values from a list of sources. The reference data is composed
    of an identifying Solution object and a dictionary mapping properties to their reference values. The reference
    values are checked against pyEQL-calculated values for a variety of engines.
    """
    engines: list[EOS] = [IdealEOS, NativeEOS, PhreeqcEOS]
    datasets = [_get_dataset(s) for s in SOURCES]
    results: dict[tuple[str, str], list[_BenchmarkEntry]] = {}

    for engine in engines:
        _patch_dataset(datasets, engine=engine)
        for i, dataset in enumerate(datasets):
            results[(engine.__name__, SOURCES[i])] = report_results(dataset)


if __name__ == "__main__":
    main()
