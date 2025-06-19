"""Solution model benchmarking utilities.

Example: Benchmark a solution model against reference data

>>> from pyEQL.benchmark import benchmark_engine
>>> from pyEQL.engines import IdealEOS

>>> results = benchmark_engine(IdealEOS(), sources=["CRC"])
>>> results["CRC"].solution_data["mean_activity"]
...

Example: Generate a reference dataset from a solution model

>>> from pyEQL import Solution
>>> from pyEQL.benchmark import calculate_stats
>>> from pyEQL.benchmark import create_entry
>>> from pyEQL.engines import IdealEOS
>>> from pyEQL.engines import NativeEOS

>>> cations = ["H[+1]", "Na[+1]", "Ca[+2]"]
>>> anions = ["OH[-1]", "Cl[-1]", "SO4[-2]"]
>>> concs = ["0.1 mol/L", "25%"]
>>> solutions = []

>>> for ions for product(cations, anions):
...     for conc in concs:
...         solutes = {ion: conc for ion in ions}
...          solutions.append(Solution(solutes=solutes, engine=IdealEOS()))

>>> solute_properties = ["activity_coefficient", "molar_conductivity"]
>>> solution_properties = ["dielectric_constant", "debye_length", "conductivity", "osmotic_coefficient", "density"]
>>> dataset = []

>>> for solution in solutions:
...     entry = create_entry(solution, solute_properties, solution_properties)
...     dataset.append(entry)

>>> for data in dataset:
...     data.solution.engine = NativeEOS()

>>> stats = calculate_stats(dataset)
"""

import json
import math
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
        source: One of "CRC", "IDST", "JPCRD", or "May2011JCED" or the path to a file containing reference data. If the
            latter, then the [path must point to a JSON which can be read into a BenchmarkEntry object.

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
        for solute, values in solute_data.items():
            solute_data[solute] = [(prop, ureg.Quantity(q)) for prop, q in values]

        for i, (prop, q) in enumerate(solution_data):
            solution_data[i] = prop, ureg.Quantity(q)

        reference.append(
            BenchmarkEntry(
                solution=Solution(**solution),
                solute_data=FormulaDict(**solute_data),
                solution_data=solution_data,
            )
        )

    return reference


def _patch_dataset(
    dataset: list[BenchmarkEntry], *, engine: EOS | Literal["native", "ideal", "phreeqc"] = "native"
) -> None:
    for data in dataset:
        data.solution.engine = engine


def _rmse(data: list[tuple[Quantity, Quantity]]) -> float:
    reduced = []

    for ref, calc in data:
        val = (ref - calc) ** 2

        if hasattr(val, "m"):
            val = val.m

        reduced.append(val)
    return math.sqrt(np.mean(reduced))


def _get_solute_property(solution: Solution, solute: str, name: str) -> Any:
    value = solution.get_property(solute, name)

    if value is None:
        if hasattr(solution, name):
            value = getattr(solution, name)
        elif hasattr(solution, f"get_{name}"):
            value = getattr(solution, f"get_{name}")(solute)
        else:
            msg = f"Solute property: {name} not supported"
            raise ValueError(msg)

    return value


def _get_mean_activity_coefficient(solution: Solution) -> float:
    activity_nu_pairs: list[tuple[float, int]] = []

    for salt_dict in solution.get_salt_dict().values():
        _ = salt_dict.pop("mol")
        salt = Salt.from_dict(salt_dict)
        act_cat = solution.get_activity_coefficient(salt.cation)
        act_an = solution.get_activity_coefficient(salt.anion)
        activity_nu_pairs.extend([(act_an, salt.nu_anion), (act_cat, salt.nu_cation)])

    factor = reduce(lambda x, y: x * y[0] ** y[1], activity_nu_pairs, 1.0)
    exponent = 1 / sum(x[1] for x in activity_nu_pairs)
    return factor**exponent


def _get_solution_property(solution: Solution, name: str) -> Any:
    if name == "mean_activity_coefficient":
        return _get_mean_activity_coefficient(solution)
    if hasattr(solution, name):
        return getattr(solution, name)

    if hasattr(solution, f"get_{name}"):
        return getattr(solution, f"get_{name}")()

    msg = f"Property {name} is not supported"
    raise ValueError(msg)


def create_entry(solution: Solution, solute_properties: list[str], solution_properties: list[str]) -> BenchmarkEntry:
    """Create a BenchmarkEntry from a Solution and specified properties.

    Args:
        solution: The Solution from which to create the entry.
        solute_properties: The solute properties to add to the entry.
        solution_properties: The solution properties to add to the entry.
    """
    solute_data = FormulaDict()

    for solute in solution.components:
        solute_data[solute] = []

        for solute_property in solute_properties:
            data = _get_solute_property(solution, solute, solute_property)
            solute_data[solute].append((solute_property, data))

    solution_data = []

    for solution_property in solution_properties:
        data = _get_solution_property(solution, solution_property)
        solution_data.append((solution_property, data))

    return BenchmarkEntry(solution=solution, solute_data=solute_data, solution_data=solution_data)


def calculate_stats(
    dataset: list[BenchmarkEntry], *, metric: Callable[[list[tuple[float, float]]], float] | None = None
) -> BenchmarkResults:
    """Calculate benchmarking statistics.

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
            for prop, reference in solute_data:
                if prop not in solute_data_pairs:
                    solute_data_pairs[prop] = []

                solute_data_pairs[prop].append((reference, _get_solute_property(d.solution, solute, prop)))

        for prop, reference in d.solution_data:
            if prop not in solution_data_pairs:
                solution_data_pairs[prop] = []

            solution_data_pairs[prop].append((reference, _get_solution_property(d.solution, prop)))

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
        key = Path(sources[i]).name
        results[key] = calculate_stats(dataset)

    return results
