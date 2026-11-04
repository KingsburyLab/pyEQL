"""
Regression tests for the industrial wastewater presets.

Each preset is an equilibrated, charge-balanced solution. Re-equilibrating it must not create or destroy any
element other than the one adjusted for charge balance. This guards against species-key collisions in
standardize_formula (e.g. Fe2(OH)2[+4] previously overwriting FeOH[+2]) and against multi-oxidation-state
elements being passed to PHREEQC incorrectly, both of which silently removed mass from the presets.
"""

import pytest
from pymatgen.core.ion import Ion

from pyEQL import Solution

INDUSTRIAL_PRESETS = [
    "ash",
    "batt_mfg",
    "batt_recycling",
    "coal_washing",
    "CRL",
    "drilling",
    "excavation",
    "FGD",
    "flotation",
    "waste_gas",
    "gasification",
    "geothermal",
    "leachate",
    "mine_drainage",
    "mine_tailings",
    "plating",
    "pw_conv",
    "pw_unconv",
    "refining",
    "semiconductor",
    "smelting",
    "tanning",
]


@pytest.mark.parametrize("preset_name", INDUSTRIAL_PRESETS)
def test_preset_reequilibration_conserves_elements(preset_name):
    s = Solution.from_preset(preset_name)
    assert s.balance_charge == "auto"
    assert s._cb_species is not None
    # elements that may legitimately change: the charge-balancing species, plus H and O (water, H+/OH-)
    exempt = {str(el) for el in Ion.from_formula(s._cb_species).elements} | {"H", "O"}
    before = {el: s.get_total_amount(el, "mol").magnitude for el in s.elements if el not in exempt}

    s.equilibrate()

    for el, amt in before.items():
        assert s.get_total_amount(el, "mol").magnitude == pytest.approx(amt, rel=1e-4), (
            f"{preset_name}: total {el} changed on re-equilibration"
        )
