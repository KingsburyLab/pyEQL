from importlib import resources

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import pyEQL
import pyEQL.benchmarks as benchmarks


def _load_benchmark_csv(name: str):
    """Load benchmark dataset from pyEQL benchmarks."""
    return pd.read_csv(
        resources.files(benchmarks).joinpath(name),
        skiprows=1,
    )


ac = _load_benchmark_csv("CRC_activity_coefficient.csv")
cond = _load_benchmark_csv("CRC_conductivity.csv")


def _ionic_strength(df):
    """calculate ionic strength of the selected data"""
    z = {"Na[+1]": 1, "K[+1]": 1, "Ca[+2]": 2, "Mg[+2]": 2, "Cl[-1]": -1, "SO4[-2]": -2}
    Is = np.zeros(len(df))
    for k, v in z.items():
        Is += df[k].values * v**2
    return 0.5 * Is


def _salt_label(ions):
    """generate salt label as plot title"""
    a = ions[0].split("[")[0]
    b = ions[1].split("[")[0]
    if a in ["Ca", "Mg"] and b == "Cl":
        return a + "Cl2"
    if a in ["Na", "K"] and b == "SO4":
        return a + "2SO4"
    return a + b


def _pyEQL_sim(ions, conc, prop):
    """pyEQL model result calculation"""
    if prop == "activity_coefficient":
        s = pyEQL.Solution([[i, f"{c} mol/kg"] for i, c in zip(ions, conc, strict=False)])
        nu = [abs(int(i.split("[")[1].split("]")[0])) for i in ions]
        gammas = [s.get_activity_coefficient(i).magnitude for i in ions]
        prod = 1.0
        for g, n in zip(gammas, nu, strict=False):
            prod *= g**n
        return prod ** (1 / sum(nu))

    if prop == "conductivity":
        s = pyEQL.Solution([[i, f"{c} mol/L"] for i, c in zip(ions, conc, strict=False)])
        return s.conductivity.to("mS/cm").magnitude

    raise ValueError(f"Unknown property: {prop}")


def _filter_df(df, ions):
    """filter selected experimental data"""
    m = (df[ions[0]] > 0) & (df[ions[1]] > 0)
    return df[m]


# main benchmark
def benchmark(df, ions, prop, n=200):
    r"""
    Plot experimental benchmark data against pyEQL calculations.

    Parameters:
        df: The experimental benchmark data as a pandas DataFrame.
        ions: The two ions used in the benchmark, specified as a list of
            ion names such as ``["Na[+1]", "Cl[-1]"]``.
        prop: The property to benchmark. Must be either
            ``"activity_coefficient"`` or ``"conductivity"``.
        n: The number of points used to generate the pyEQL simulation curve.

    Returns:
        None.
        Displays a plot comparing the experimental data with the
        corresponding pyEQL calculations.

    Notes:
        The experimental data are filtered to include only measurements
        where both selected ions have non-zero concentrations. Ionic
        strength is calculated from the ion concentrations, and the
        corresponding pyEQL results are generated over a range of ionic
        strengths.

        For ``"activity_coefficient"``, the mean activity coefficient
        calculated by pyEQL is compared with the experimental data.

        For ``"conductivity"``, the conductivity calculated by pyEQL is
        compared with the experimental data. Conductivity is reported in
        mS/cm.

    """

    filtered_df = _filter_df(df, ions)
    Is = _ionic_strength(filtered_df)
    y = filtered_df["mean_activity_coefficient"] if prop == "activity_coefficient" else filtered_df["conductivity"]

    I_sim = np.linspace(Is.min(), Is.max() + 0.1, n)
    base_conc = filtered_df[ions].mean().values
    y_sim = []
    for I_s in I_sim:
        scale = I_s / np.mean(Is) if np.mean(Is) != 0 else 1.0
        conc = base_conc * scale
        y_sim.append(_pyEQL_sim(ions, conc, prop))
    y_sim = np.array(y_sim)

    ylabel_map = {"activity_coefficient": "Mean Activity Coefficient", "conductivity": "Conductivity (mS/cm)"}
    plt.plot(I_sim, y_sim, color="r", lw=2, label="pyEQL")
    plt.scatter(Is, y, s=18, color="k", marker="x", label="Experiment")
    plt.xlabel("Ionic Strength (m)")
    plt.ylabel(ylabel_map[prop])
    plt.title(_salt_label(ions))
    plt.legend()
    plt.tight_layout()
    plt.show()
