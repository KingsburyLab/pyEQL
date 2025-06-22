"""Solution model benchmarking utilities.

Example: Benchmark a solution model against reference data

>>> from pyEQL.benchmark import benchmark_engine
>>> from pyEQL.engines import IdealEOS

>>> results = benchmark_engine(IdealEOS(), sources=["crc"])
>>> results["crc"].solution_stats["mean_activity_coefficient"]
...

Example: Generate a reference dataset from a solution model engine

>>> from pyEQL import Solution
>>> from pyEQL.benchmark import calculate_stats
>>> from pyEQL.benchmark import create_entry
>>> from pyEQL.benchmark import _create_solution_key
>>> from pyEQL.engines import IdealEOS
>>> from pyEQL.engines import NativeEOS

>>> cations = ["H[+1]", "Na[+1]", "Ca[+2]"]
>>> anions = ["OH[-1]", "Cl[-1]", "SO4[-2]"]
>>> concs = ["0.1 mol/L", "25%"]
>>> solute_properties = ["activity_coefficient", "molar_conductivity"]
>>> solution_properties = ["dielectric_constant", "debye_length", "conductivity", "osmotic_coefficient", "density"]
>>> ideal_solutions = []
>>> native_solutions = []

>>> for ions for product(cations, anions):
...     for conc in concs:
...         solutes = {ion: conc for ion in ions}
...         ideal_solutions.append(Solution(solutes=solutes, engine=IdealEOS()))
...         native_solutions.append(Solution(solutes=solutes, engine=NativeEOS()))

>>> ideal_data = {}
>>> native_data = {}

>>> for ideal_solution, native_solution in zip(ideal_solutions, native_solutions, strict=True):
...     ideal_entry = create_entry(ideal_solution, solute_properties, solution_properties)
...     native_entry = create_entry(native_solution, solute_properties, solution_properties)
...     key = _create_solution_key(ideal_solution)
...     ideal_data[key] = ideal_entry
...     native_data[key] = ideal_entry

>>> stats = calculate_stats(ideal_data, native_data)
"""

import json
import math
from collections.abc import Callable
from functools import reduce
from importlib.resources import files
from pathlib import Path
from typing import Any, Literal, NamedTuple

import numpy as np
from pint import Quantity

import pyEQL
from pyEQL import ureg
from pyEQL.engines import EOS
from pyEQL.salt_ion_match import Salt
from pyEQL.solution import Solution
from pyEQL.utils import standardize_formula

INTERNAL_SOURCES: list[str] = ["CRC", "IDST", "JPCRD", "May2011JCED"]


class BenchmarkEntry(NamedTuple):
    """A set of data for a solution.

    Attributes:
        solution: Solution
            The Solution to which the reference data applies.
        solute_data: dict[str, dict[str, Quantity]]
            A dictionary mapping solutes to a dictionary mapping solute properties to a Quantity.
        solution_data: dict[str, Quantity]
            A dictionary mapping solution properties to a Quantity.

    The property strings that serve as keys in ``solute_data`` and ``solution_data`` should correspond to properties
    that can be retrieved using the formalisms outlined in ``_get_solute_property`` and ``_get_solution_property``.
    """

    solution: Solution
    solute_data: dict[str, dict[str, Quantity]] = {}
    solution_data: dict[str, Quantity] = {}


class BenchmarkResults(NamedTuple):
    """Solute and solution stats from :func:`pyEQL.benchmark.calculate_stats`.

    Attributes:
        solute_stats: dict[str, float]
            Benchmarking statistics for solute properties.
        solution_stats: dict[str, float]
            Benchmarking statistics for solution properties.
    """

    solute_stats: dict[str, float]
    solution_stats: dict[str, float]


# TODO: Admittedly, this is an ugly solution to enable Solutions to serve as dictionary keys to speed up the
# calculation of benchmarking stats. For now, we wrap the decision on the final key format in a function
# (_create_solution_key) that produces hashable values from a Solution, and we abstract the key format with a type
# alias (SolutionKey).
SolutionKey = tuple[
    # Composition: (ion, concentration)
    tuple[tuple[str, str], ...],
    # State: temperature, pressure
    # ? should other state variables be checked (e.g., pH, pE)?
    tuple[str, str],
]


def _create_solution_key(solution: Solution) -> SolutionKey:
    vol = solution.volume.magnitude
    composition = []
    for component in sorted(solution.components):
        if component != solution.solvent:
            amount = f"{solution.components[component] / vol} mol/L"
            composition.append((component, amount))
    state = solution.temperature, solution.pressure
    return tuple(composition), state


def load_dataset(
    source: str | Path,
    *,
    solutions: list[Solution] | None = None,
    solute_properties: list[tuple[str, str]] | None = None,
    solution_properties: list[str] | None = None,
) -> dict[SolutionKey, BenchmarkEntry]:
    """Load reference dataset.

    Args:
        source: str | Path
            One of "CRC", "IDST", "JPCRD", or "May2011JCED" or the path to a file containing reference data. If the
            latter, then the [path must point to a JSON which can be read into a BenchmarkEntry object.
        solutions: list[Solution], optional
            The solutions for which data will be loaded from the dataset. If provided, only data corresponding to
            solutions with the same composition, temperature, and pressure will be loaded. If omitted, all
            compositions and conditions in the reference data contained in ``sources`` will be used for the
            benchmarking.
        solute_properties: list[tuple[str, str]], optional
            The solute properties to include in the benchmarking, specified as ``(solute, property)``. The engine will
            only be benchmarked against those solute properties listed here. If omitted, the engine will be benchmarked
            against all solute properties in ``sources``.
        solution_properties: list[str], optional
            The solution properties to include in the benchmarking. The engine will only be benchmarked against those
            solution properties listed here. Defaults to None in which case the engine will be benchmarked against all
            solute properties in ``sources``.

    Returns:
        A dictionary mapping SolutionKey to BenchmarkEntry objects. See the comment over SolutionKey for details about
        its structure.
    """
    match str(source).lower():
        case "crc" | "idst" | "jpcrd":
            source = files("pyEQL").joinpath("database", f"{str(source).lower()}.json")
        case _:
            source = Path(source)

    with source.open(mode="r", encoding="utf-8") as file:
        data: list[tuple[dict, dict[str, list], list[tuple]]] = json.load(file)

    reference: dict[SolutionKey, BenchmarkEntry] = {}
    solution_keys = [] if solutions is None else [_create_solution_key(s) for s in solutions]

    for raw_solution, raw_solute_data, raw_solution_data in data:
        # TODO: handle weight % concentration data
        solution = Solution(**raw_solution)
        soln_key = _create_solution_key(solution)

        if solutions is not None and soln_key not in solution_keys:
            continue

        solute_data: dict[str, dict[str, Quantity]] = {}
        solution_data = {}

        for solute, values in raw_solute_data.items():
            standardized_solute = standardize_formula(solute)
            solute_data[standardized_solute] = {}
            for p, q in values.items():
                if solute_properties is None or (standardized_solute, p) in solute_properties:
                    solute_data[standardized_solute][p] = ureg.Quantity(q)

        for p, q in raw_solution_data.items():
            if solution_properties is None or p in solution_properties:
                solution_data[p] = ureg.Quantity(q)

        reference[soln_key] = BenchmarkEntry(
            solution=solution,
            solute_data=solute_data,
            solution_data=solution_data,
        )

    return reference


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


def _get_activity_coefficient_pitzer(solution: Solution) -> float:
    salt = solution.get_salt()
    param = solution.get_property(salt.formula, "model_parameters.activity_pitzer")
    alpha1 = 2
    alpha2 = 0
    molality = salt.get_effective_molality(solution.ionic_strength)
    temperature = str(solution.temperature)

    return pyEQL.activity_correction.get_activity_coefficient_pitzer(
        solution.ionic_strength,
        molality,
        alpha1,
        alpha2,
        ureg.Quantity(param["Beta0"]["value"]).magnitude,
        ureg.Quantity(param["Beta1"]["value"]).magnitude,
        ureg.Quantity(param["Beta2"]["value"]).magnitude,
        ureg.Quantity(param["Cphi"]["value"]).magnitude,
        salt.z_cation,
        salt.z_anion,
        salt.nu_cation,
        salt.nu_anion,
        temperature,
    )


def _get_solution_property(solution: Solution, name: str) -> Any:
    if name == "mean_activity_coefficient":
        return _get_mean_activity_coefficient(solution)

    if name == "activity_coefficient_pitzer":
        return _get_activity_coefficient_pitzer(solution)

    if hasattr(solution, name):
        return getattr(solution, name)

    if hasattr(solution, f"get_{name}"):
        return getattr(solution, f"get_{name}")()

    msg = f"Property {name} is not supported"
    raise ValueError(msg)


def create_entry(
    solution: Solution,
    solute_properties: list[tuple[str, str]],
    solution_properties: list[str],
) -> BenchmarkEntry:
    """Create a BenchmarkEntry from a Solution and specified properties.

    Args:
        solution: Solution
            The Solution from which to create the entry.
        solute_properties: list[tuple[str, str]]
            The solute properties to add to the entry.
        solution_properties: list[str]
            The solution properties to add to the entry.
    """
    solute_data = {}

    for solute, solute_property in solute_properties:
        standardized_solute = standardize_formula(solute)
        if standardized_solute not in solute_data:
            solute_data[standardized_solute] = {}
        data = _get_solute_property(solution, standardized_solute, solute_property)
        solute_data[standardized_solute][solute_property] = data

    solution_data = {}

    for solution_property in solution_properties:
        data = _get_solution_property(solution, solution_property)
        solution_data[solution_property] = data

    return BenchmarkEntry(solution=solution, solute_data=solute_data, solution_data=solution_data)


def _create_engine_dataset(
    engine: EOS | Literal["native", "ideal", "phreeqc"] = "native",
    datasets: list[dict[SolutionKey, BenchmarkEntry]] | None = None,
) -> dict[SolutionKey, BenchmarkEntry]:
    mapper: dict[
        SolutionKey,
        tuple[
            # solute_properties
            list[tuple[str, str]],
            # solution_properties
            list[str],
        ],
    ] = {}

    for dataset in datasets:
        for key, entry in dataset.items():
            if key not in mapper:
                mapper[key] = [], []

            for solute in entry.solute_data:
                mapper[key][0].extend((solute, prop) for prop in list(entry.solute_data[solute]))
            mapper[key][1].extend(entry.solution_data)

    engine_dataset = {}

    for solution_key, (solute_properties, solution_properties) in mapper.items():
        temperature, pressure = solution_key[1]
        solution = Solution(solutes=dict(solution_key[0]), temperature=temperature, pressure=pressure, engine=engine)
        entry = create_entry(
            solution=solution,
            solute_properties=sorted(set(solute_properties)),
            solution_properties=sorted(set(solution_properties)),
        )
        engine_dataset[key] = entry

    return engine_dataset


def calculate_stats(
    reference: dict[SolutionKey, BenchmarkEntry],
    calculated: dict[SolutionKey, BenchmarkEntry],
    *,
    metric: Callable[[list[tuple[float, float]]], float] | None = None,
) -> BenchmarkResults:
    """Calculate benchmarking statistics.

    Args:
        reference: dict[SolutionKey, BenchmarkEntry]
            The reference data, a dictionary mapping SolutionKeys to BenchmarkEntry object.
        calculated: dict[SolutionKey, BenchmarkEntry]
            The data to be benchmarked against the reference data,  dictionary mapping SolutionKeys to
            BenchmarkEntry object.
        metric: Callable[[list[tuple[float, float]]], float], optional
            A function that acts on the list of 2-tuples (reference, calculated), which contains reference and
            calculated values. This function should calculate a statistical metric for the list. Defaults to the root-
            mean-squared error.

    Returns:
        A 2-tuple (``solute_stats``, ``solution_stats``) where ``solute_stats`` and ``solution_stats`` are dictionaries
        mapping solute and solution properties, respectively, to their benchmark statistics.
    """

    def _rmse(data: list[tuple[Quantity, Quantity]]) -> float:
        return math.sqrt(np.mean([Quantity(ref - calc).m ** 2 for ref, calc in data]))

    metric = metric or _rmse

    # property: [(reference, calculated)]
    solute_data_pairs: dict[str, list[tuple[float, float]]] = {}
    solution_data_pairs: dict[str, list[tuple[float, float]]] = {}

    for key, entry in reference.items():
        for solute in entry.solute_data:
            for prop, value in entry.solute_data[solute].items():
                if prop not in solute_data_pairs:
                    solute_data_pairs[prop] = []
                calculated_value = calculated[key].solute_data[solute][prop]
                solute_data_pairs[prop].append((value, calculated_value))

        for prop, value in entry.solution_data.items():
            if prop not in solution_data_pairs:
                solution_data_pairs[prop] = []
            calculated_value = calculated[key].solution_data[prop]
            solution_data_pairs[prop].append((value, calculated_value))

    solute_stats = {k: metric(v) for k, v in solute_data_pairs.items()}
    solution_stats = {k: metric(v) for k, v in solution_data_pairs.items()}

    return BenchmarkResults(solute_stats=solute_stats, solution_stats=solution_stats)


def benchmark_engine(
    engine: EOS,
    *,
    sources: list[str] | None = None,
    solutions: list[Solution] | None = None,
    solute_properties: list[tuple[str, str]] | None = None,
    solution_properties: list[str] | None = None,
) -> BenchmarkResults:
    """Benchmark a modeling engine against reference data.

    Args:
        engine: EOS
            The modeling engine to benchmark.
        sources: list[str], optional
            One of INTERNAL_SOURCES or the path to a JSON file that can be read into a list of BenchmarkEntry
            objects. Defaults to INTERNAL_SOURCES.
        solutions: list[Solution], optional
            The solutions for which data will be loaded from the dataset. If provided, only data corresponding to
            solutions with the same components, concentrations, and conditions (temperature, pressure) will be loaded.
            If omitted, reference data for all components, concentrations, and conditions (temperature, pressure)
            contained in ``sources`` will be used for the benchmarking.
        solute_properties: list[tuple[str, str]], optional
            The solute properties to include in the benchmarking, specified as ``(solute, property)``. The engine will
            only be benchmarked against those solute properties listed here. If omitted, the engine will be benchmarked
            against all solute properties in ``sources``.
        solution_properties: list[str], optional
            The solution properties to include in the benchmarking. The engine will only be benchmarked against those
            solution properties listed here. Defaults to None in which case the engine will be benchmarked against all
            solute properties in ``sources``.

    Returns:
        A dictionary mapping source names to the corresponding solute and solution statistical metrics.
    """
    sources = sources or INTERNAL_SOURCES
    datasets = [
        load_dataset(
            s, solutions=solutions, solute_properties=solute_properties, solution_properties=solution_properties
        )
        for s in sources
    ]
    engine_dataset = _create_engine_dataset(engine, datasets)
    return {sources[i]: calculate_stats(ref, engine_dataset) for i, ref in enumerate(datasets)}
