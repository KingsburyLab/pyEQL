"""Solution model benchmarking utilities.


Usage:

    python pyeql.benchmark

"""

import json
import math
from collections.abc import Callable
from pathlib import Path

import numpy as np

from pyEQL import ureg
from pyEQL.solution import Solution

# Maps solutes to property values dictionary
ReferenceData = dict[str, dict[str, float]]
# TODO: Select and validate data sources
SOURCES = []
# TODO: revisit tolerance
# relative tolerance between experimental and computed properties for this test file
RTOL = 0.05

s1 = Solution(volume="2 L")
s2 = Solution([["Na+", "4 mol/L"], ["Cl-", "4 mol/L"]], volume="2 L")


# TODO: edit to save Solution objects and add conductivity data to reference data dictionary
def _create_conductivity_data() -> None:
    """This method edits data in place."""
    data: list[tuple[Solution, ReferenceData]] = []
    DEST = Path(__file__).parent.joinpath("database")

    # per CRC handbook - "electrical conductiVity of Water" , conductivity of pure water
    # at 25 and 100 C is 0.0550 and 0.765 uS/cm
    assert np.isclose(s1.conductivity.to("uS/cm").magnitude, 0.055, atol=1e-3)

    # TODO - seems to be a possible bug related to setting temperature here
    # s1.temperature = "100 degC"
    # s2 = Solution(temperature='100 degC')
    # assert np.isclose(s1.conductivity.to('uS/cm').magnitude, 0.765, atol=1e-3)

    # CRC handbook table - "equivalent conductivity of electrolytes in aqueous solution"
    # nacl
    for conc, cond in zip([0.001, 0.05, 0.1], [123.68, 111.01, 106.69], strict=False):
        s1 = Solution({"Na+": f"{conc} mol/L", "Cl-": f"{conc} mol/L"})
        assert np.isclose(
            s1.conductivity.to("S/m").magnitude, conc * cond / 10, atol=0.5
        ), f"Conductivity test failed for NaCl at {conc} mol/L. Result = {s1.conductivity.to('S/m').magnitude}"

    # higher concentration data points from Appelo, 2017 Figure 4.
    s1 = Solution({"Na+": "2 mol/kg", "Cl-": "2 mol/kg"})
    assert np.isclose(s1.conductivity.to("mS/cm").magnitude, 145, atol=10)

    # MgCl2
    for conc, cond in zip([0.001, 0.05, 0.1], [124.15, 114.49, 97.05], strict=False):
        s1 = Solution({"Mg+2": f"{conc} mol/L", "Cl-": f"{2 * conc} mol/L"})
        assert np.isclose(
            s1.conductivity.to("S/m").magnitude, 2 * conc * cond / 10, atol=1
        ), f"Conductivity test failed for MgCl2 at {conc} mol/L. Result = {s1.conductivity.to('S/m').magnitude}"

    # per CRC handbook "standard KCl solutions for calibrating conductivity cells", 0.1m KCl has a conductivity of
    # 12.824 mS/cm at 25 C
    s_kcl = Solution({"K+": "0.1 mol/kg", "Cl-": "0.1 mol/kg"})
    assert np.isclose(s_kcl.conductivity.magnitude, 1.2824, atol=0.25)  # conductivity is in S/m

    # TODO - expected failures due to limited temp adjustment of diffusion coeff
    s_kcl.temperature = "5 degC"
    assert np.isclose(s_kcl.conductivity.magnitude, 0.81837, atol=0.2)

    s_kcl.temperature = "50 degC"
    assert np.isclose(s_kcl.conductivity.magnitude, 1.91809, atol=0.2)

    # TODO - conductivity model not very accurate at high conc.
    s_kcl = Solution({"K+": "1 mol/kg", "Cl-": "1 mol/kg"})
    assert np.isclose(s_kcl.conductivity.magnitude, 10.862, atol=0.45)

    solutions = [s1, s2, s_kcl]

    for solution in solutions:
        reference_data = {"conductivity": solution.conductivity}
        data.append((solution, reference_data))

    filename = DEST.joinpath("CRC_data.json")
    with filename.open(mode="r", encoding="utf-8") as file:
        json.dump(data, file, indent=4)


# TODO: edit to save Solution objects and add diffusion data to reference data dictionary
def _create_diffusion_data() -> None:
    # test ionic strength adjustment
    assert s1.get_diffusion_coefficient("H+") > s2.get_diffusion_coefficient("H+")

    # for Na+, d=122, a1=1.52, a2=3.7, A=1.173802/2.303 at 25 DegC, B = 3.2843078+10
    factor = np.exp(
        -1.52
        * 1.173802
        / 2.303
        * 1
        * np.sqrt(s2.ionic_strength.magnitude)
        / (1 + 3.2843078e10 * np.sqrt(s2.ionic_strength.magnitude) * 3.7 / (1 + s2.ionic_strength.magnitude**0.75))
    )
    assert np.isclose(
        factor * s2.get_diffusion_coefficient("Na+").magnitude,
        s2.get_diffusion_coefficient("Na+").magnitude,
        atol=5e-11,
    )
    s_dilute = Solution({"Na+": "1 mmol/L", "Cl-": "1 mmol/L"})
    assert np.isclose(
        s_dilute.get_diffusion_coefficient("Na+", activity_correction=False).magnitude, 1.334e-9, atol=1e-12
    )
    assert np.isclose(s_dilute.get_transport_number("Na+"), 0.396, atol=1e-3)
    assert np.isclose(s_dilute.get_transport_number("Cl-"), 0.604, atol=1e-3)

    # test setting a default value
    s2.default_diffusion_coeff = 0
    assert s2.get_diffusion_coefficient("Cs+").magnitude == 0
    s2.default_diffusion_coeff = 1e-9
    assert s2.get_diffusion_coefficient("Cs+", activity_correction=False).magnitude == 1e-9
    s2.default_diffusion_coeff = 0
    assert s2.get_diffusion_coefficient("Cs+", activity_correction=True).magnitude < 1e-9
    d25 = s2.get_diffusion_coefficient("Na+", activity_correction=False).magnitude
    nu25 = s2.water_substance.nu
    s2.temperature = "40 degC"
    d40 = s2.get_diffusion_coefficient("Na+", activity_correction=False).magnitude
    nu40 = s2.water_substance.nu
    assert np.isclose(
        d40,
        d25 * np.exp(122 / (273.15 + 40) - 122 / 298.15) * (nu25 / nu40),
        atol=5e-11,
    )

    # test correction factors for concentration, as per Appelo 2017 Fig 5
    D1 = Solution({"Na+": "1 umol/L", "Cl-": "1 umol/L"}).get_diffusion_coefficient("Na+").magnitude
    D2 = Solution({"Na+": "1.7 mol/kg", "Cl-": "1.7 mol/kg"}).get_diffusion_coefficient("Na+").magnitude
    assert np.isclose(D2 / D1, 0.54, atol=1e-2)

    D1 = Solution({"K+": "1 umol/L", "Cl-": "1 umol/L"}).get_diffusion_coefficient("K+").magnitude
    D2 = Solution({"K+": "0.5 mol/kg", "Cl-": "0.5 mol/kg"}).get_diffusion_coefficient("K+").magnitude
    assert np.isclose(D2 / D1, 0.80, atol=1e-2)


def _create_activity_data() -> None:
    # see test_activity.py & test_mixed_electrolyte_activity.py
    pass


def _create_density_data() -> None:
    # see test_density.py
    pass


def _create_osmotic_coefficient_data() -> None:
    # see test_osmotic_coeff.py
    pass


# TODO: edit to save Solution objects and add molar conductivity data to reference data dictionary
def _create_molar_conductivity_data() -> None:
    def test_molar_conductivity_potassium(self):
        # K+ - 73.48 x 10 ** -4 m ** 2 S / mol
        s1 = Solution([["K+", "0.001 mol/L"], ["Cl-", "0.001 mol/L"]], temperature="25 degC")
        result = s1.get_molar_conductivity("K+").to("m**2*S/mol").magnitude
        expected = ureg.Quantity("73.48e-4 m**2 * S / mol").magnitude

        assert np.isclose(result, expected, rtol=RTOL)

    def test_molar_conductivity_sodium(self):
        # Na+ - 50.08 x 10 ** -4 m ** 2 S / mol
        s1 = Solution([["Na+", "0.001 mol/L"], ["Cl-", "0.001 mol/L"]], temperature="25 degC")
        result = s1.get_molar_conductivity("Na+").to("m**2*S/mol").magnitude
        expected = ureg.Quantity("50.08e-4 m**2 * S / mol").magnitude

        assert np.isclose(result, expected, rtol=RTOL)

    def test_molar_conductivity_magnesium(self):
        # Mg+2 - 106 x 10 ** -4 m ** 2 S / mol
        s1 = Solution([["Mg+2", "0.001 mol/L"], ["Cl-", "0.002 mol/L"]], temperature="25 degC")
        result = s1.get_molar_conductivity("Mg+2").to("m**2*S/mol").magnitude
        expected = ureg.Quantity("106e-4 m**2 * S / mol").magnitude

        assert np.isclose(result, expected, atol=0.005)

    def test_molar_conductivity_chloride(self):
        # Cl- - 76.31 x 10 ** -4 m ** 2 S / mol
        s1 = Solution([["Na+", "0.001 mol/L"], ["Cl-", "0.001 mol/L"]], temperature="25 degC")
        result = s1.get_molar_conductivity("Cl-").to("m**2*S/mol").magnitude
        expected = ureg.Quantity("76.31e-4 m**2 * S / mol").magnitude

        assert np.isclose(result, expected, rtol=RTOL)

    def test_molar_conductivity_fluoride(self):
        # F- - 55.4 x 10 ** -4 m ** 2 S / mol
        s1 = Solution([["Na+", "0.001 mol/L"], ["F-", "0.001 mol/L"]], temperature="25 degC")
        result = s1.get_molar_conductivity("F-").to("m**2*S/mol").magnitude
        expected = ureg.Quantity("55.4e-4 m**2 * S / mol").magnitude

        assert np.isclose(result, expected, rtol=RTOL)

    def test_molar_conductivity_sulfate(self):
        # SO4-2 - 160 x 10 ** -4 m ** 2 S / mol
        s1 = Solution([["Na+", "0.002 mol/L"], ["SO4-2", "0.001 mol/L"]], temperature="25 degC")
        result = s1.get_molar_conductivity("SO4-2").to("m**2*S/mol").magnitude
        expected = ureg.Quantity("160.0e-4 m**2 * S / mol").magnitude

        assert np.isclose(result, expected, atol=0.002)

    def test_molar_conductivity_hydroxide(self):
        # OH- - 198 x 10 ** -4 m ** 2 S / mol
        s1 = Solution(temperature="25 degC")
        result = s1.get_molar_conductivity("OH-").to("m**2*S/mol").magnitude
        expected = ureg.Quantity("198e-4 m**2 * S / mol").magnitude

        assert np.isclose(result, expected, rtol=RTOL)

    def test_molar_conductivity_hydrogen(self):
        # H+ - 349.65 x 10 ** -4 m ** 2 S / mol
        s1 = Solution(temperature="25 degC")
        result = s1.get_molar_conductivity("H+").to("m**2*S/mol").magnitude
        expected = ureg.Quantity("349.65e-4 m**2 * S / mol").magnitude

        assert np.isclose(result, expected, rtol=RTOL)

    # molar conductivity of a neutral solute should be zero
    def test_molar_conductivity_neutral(self):
        s1 = Solution([["FeCl3", "0.001 mol/L"]], temperature="25 degC")
        result = s1.get_molar_conductivity("FeCl3").to("m**2*S/mol").magnitude
        expected = ureg.Quantity(0, "m**2 * S / mol").magnitude

        assert round(abs(result - expected), 5) == 0

    # molar conductivity of water should be zero
    def test_molar_conductivity_water(self):
        s1 = Solution(temperature="25 degC")
        result = s1.get_molar_conductivity("H2O").to("m**2*S/mol").magnitude
        expected = ureg.Quantity(0, "m**2 * S / mol").magnitude

        assert round(abs(result - expected), 5) == 0


# TODO: write database files and write loading function
def _load_database(source: str) -> list[tuple[Solution, dict[str, float]]]:
    pass


def _get_reference(source: str) -> list[tuple[Solution, ReferenceData]]:
    """Load reference data.

    Args:
        source: One of "IDST" or "PHREEQC" or the path to a file containing reference data. If the latter, then the
            path must point to a JSON which can be read into a list of Solution, Result pairs. Each Solution can only
            contain the following keys: solutes, volume, temperature, pressure, pH, pE, balance_charge, and solvent.
            The Result item must contain a dictionary mapping property names to their reference values for the
            corresponding Solution.

    Returns:
        A list of Solution, Result 2-tuples.
    """
    match source:
        case "IDST" | "PHREEQC":
            reference = _load_database(source)
        case _:
            with Path(source).open(mode="r", encoding="utf-8") as file:
                data = json.load(file)

            reference: list[tuple[Solution, dict[str, float]]] = []

            for sol, values in data:
                reference.append((Solution(**sol), values))

    return reference


def benchmark_engine(
    references: list[tuple[Solution, ReferenceData]], *, getter: Callable[[Solution, str, str], float] | None = None
) -> dict[str, float]:
    """Calculate the RMSE for a set of solutions with respect to reference data.


    Args:
        references: A dictionary mapping solution properties to their reference values.
        getter: A Callable that should accept three positional arguments: A Solution object, a solute string, and a
            string indicating the name of a solvent property. The Callable should return the property value from the
            Solution. Defaults to the Solution._get_property method.

    Returns:
        A dictionary mapping each key in the reference to the cumulative root-mean-squared error across all solutions.
    """

    def _get_property(s: Solution, solute: str, name: str) -> float:
        return s._get_property(solute, name)

    getter = getter or _get_property

    errors = {}

    for solution, data in references:
        for solute, solute_data in data.items():
            for property, reference_value in solute_data.items():
                calculated = _get_property(solution, solute, property)

                if property not in errors:
                    errors[property] = []

                errors[property].append(math.abs(calculated - reference_value))

    return {k: np.std(values) for k, values in errors.items()}


# TODO: Decide on reporting method and metrics
def report_results(results: dict[tuple[str, str], float]) -> None:
    """Report the results of the benchmarking.

    Args:
        results: A dictionary mapping 2-tuples (engine, source) to the associated root-means-squared error across
            all data points in the source for the given engine.
    """


def run_benchmark() -> None:
    """Solution benchmarking function.

    This function works by reading in reference property values from a list of sources. The reference data is composed
    of an identifying Solution object and a dictionary mapping properties to their reference values. The reference
    values are checked against pyEQL-calculated values for a variety of engines.
    """
    # TODO: create engines (e.g., Pitzer, PHREEQC)
    # see test_phreeqc.py
    engines = []
    references = [(s, _get_reference(s)) for s in SOURCES]
    results: dict[tuple[str, str], float] = {}

    for eng in engines:
        for source, ref in references:
            results[(eng, source)] = benchmark_engine(ref)

    report_results(results)


if __name__ == "__main__":
    run_benchmark()
