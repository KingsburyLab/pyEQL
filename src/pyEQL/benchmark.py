"""Solution model benchmarking utilities.

Usage:

>>> from pyEQL.benchmark import benchmark_engine
>>> from pyEQL.engines import IdealEOS

>>> results = benchmark_engine(IdealEOS(), sources=["CRC"])
>>> results["CRC"].solution_data["mean_activity"]
...
"""

import json
from collections.abc import Callable
from functools import reduce
from pathlib import Path
from typing import Any, Literal, NamedTuple

import numpy as np
from pint import Quantity

from pyEQL import ureg
from pyEQL.engines import EOS
from pyEQL.salt_ion_match import Salt
from pyEQL.solution import Solution
from pyEQL.utils import FormulaDict

# TODO: Select and validate data sources
# If all solution reference data are generated from the same solutions, then solutions can be used as an input into
# source creation
INTERNAL_SOURCES: list[str] = ["CRC", "IDST", "JPCRD", "May2011JCED"]


class BenchmarkEntry(NamedTuple):
    """Solution reference data entry.

    Attributes:
        solution: The Solution to which the reference data applies.
        solute_data: A dictionary mapping solutes to a list of solute property-value 2-tuples.
        solution_data: A list of solution property-quantity 2-tuples.

    The property strings in the 2-tuples in `solute_data` and `solution_data` should correspond to properties that can
    be retrieved using the formalisms outlined in `_get_solute_property` and `_get_solution_property`.
    """

    solution: Solution
    solute_data: FormulaDict = FormulaDict()
    solution_data: list[tuple[str, Quantity]] = []


class BenchmarkResults(NamedTuple):
    """Solute and solution stats from :func:`pyEQL.benchmark.benchmark_engine`."""

    solute_stats: dict[str, float]
    solution_stats: dict[str, float]


def get_dataset(source: str | Path) -> list[BenchmarkEntry]:
    """Load reference dataset.

    Args:
        source: One of "CRC", "IDST", "JPCRD", or "May2011JCED" or the path to a file containing reference data. If the latter,
            then the [path must point to a JSON which can be read into a BenchmarkEntry object.

    Returns:
        A list of BenchmarkEntry objects one for each data point in the data set.
    """
    match source:
        case "CRC" | "IDST" | "JPCRD":
            source = Path(__file__).parent.joinpath("database", f"{source}.json")
        case _:
            source = Path(source)

    with source.open(mode="r", encoding="utf-8") as file:
        data: list[tuple[dict, dict[str, list], list[tuple]]] = json.load(file)

    reference: list[BenchmarkEntry] = []

    for solution, solute_data, solution_data in data:
        for k, values in solute_data.items():
            solute_data[k] = [ureg.Quantity(float(x), y) for x, y in values]

        for i, (x, y) in solution_data:
            solution_data[i] = ureg.Quantity(float(x), y)

        reference.append(
            BenchmarkEntry(solution=Solution.from_dict(solution), solute_data=solute_data, solution_data=solution_data)
        )

    return reference


def _patch_dataset(
    dataset: list[BenchmarkEntry], *, engine: EOS | Literal["native", "ideal", "phreeqc"] = "native"
) -> None:
    for data in dataset:
        data.solution.engine = engine


def _rmse(data: list[tuple[float, float]]) -> float:
    return np.std([ref - calc for ref, calc in data])


def _get_solute_property(solution: Solution, solute: str, name: str) -> Any:
    value = solution.get_property(solute, name)

    if value is None:
        msg = f"Solute property: {name} not supported"
        raise ValueError(msg)

    return value


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


def _get_solution_property(solution: Solution, name: str) -> Any:
    if name == "mean_activity":
        return _get_mean_activity(solution)
    if hasattr(solution, name):
        return getattr(solution, name)

    if hasattr(solution, f"get_{name}"):
        return getattr(solution, f"get_{name}")

    msg = f"Property {name} is not supported"
    raise ValueError(msg)


def report_results(
    dataset: list[BenchmarkEntry], *, metric: Callable[[list[tuple[float, float]]], float] | None = None
) -> BenchmarkResults:
    """Report the results of the benchmarking.

    Args:
        dataset: A list of BenchmarkEntry objects
        metric: A function that acts on the list of 2-tuples (reference, calculated), which contains reference and
            calculated values. This function should calculate a statistical metric for the list. Defaults to the root-
            mean-squared error.

    Returns:
        A 2-tuple (`solute_stats`, `solution_stats`) where `solute_stats` and `solution_stats` are dictionaries mapping
        solute and solution properties, respectively, to their benchmark statistics.
    """
    metric = metric or _rmse

    # property: [(reference, calculated)]
    solute_data_pairs: dict[str, list[tuple[float, float]]] = {}
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

            solution_data_pairs[property].append((reference, _get_solution_property(d.solution, property)))

    solute_stats = {k: metric(v) for k, v in solute_data_pairs.items()}
    solution_stats = {k: metric(v) for k, v in solution_data_pairs.items()}

    return BenchmarkResults(solute_stats=solute_stats, solution_stats=solution_stats)


def benchmark_engine(engine: EOS, *, sources: list[str] | None = None) -> BenchmarkResults:
    """Benchmark a modeling engine against reference data.

    Args:
        engine: The modeling engine to benchmark.
        sources: One of INTERNAL_SOURCES or the path to a JSON file that can be read into a list of BenchmarkEntry
            objects. Defaults to INTERNAL_SOURCES.

    Returns:
        A dictionary mapping source names to the corresponding solute and solution statistical metrics.
    """
    sources = sources or INTERNAL_SOURCES
    datasets = [get_dataset(s) for s in sources]
    results: BenchmarkResults = {}

    for i, dataset in enumerate(datasets):
        _patch_dataset(dataset, engine=engine)
        results[sources[i]] = report_results(dataset)

    return results
