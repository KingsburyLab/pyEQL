import numpy as np

from pyEQL import Solution


def test_unit_conversions():
    sol = Solution({"Ca+2": "1000 ppb"})
    ppb_value = sol.get_amount("Ca+2", "ppb").magnitude
    microgram_per_liter = sol.get_amount("Ca+2", "microgram/L").magnitude
    assert np.isclose(ppb_value, 1000), f"Expected 1000, got {ppb_value}"
    assert np.isclose(microgram_per_liter, 1000), f"Expected 1000, got {microgram_per_liter}"

    sol = Solution({"Ca+2": "1 ppm"})
    ppm_value = sol.get_amount("Ca+2", "ppm").magnitude
    milligram_per_liter = sol.get_amount("Ca+2", "milligram/L").magnitude
    assert np.isclose(ppm_value, 1), f"Expected 1, got {ppm_value}"
    assert np.isclose(milligram_per_liter, 1), f"Expected 1, got {milligram_per_liter}"

    sol = Solution({"Ca+2": "1 ppt"})
    ppt_value = sol.get_amount("Ca+2", "ppt").magnitude
    nanogram_per_liter = sol.get_amount("Ca+2", "nanogram/L").magnitude
    assert np.isclose(ppt_value, 1), f"Expected 1, got {ppt_value}"
    assert np.isclose(nanogram_per_liter, 1), f"Expected 1, got {nanogram_per_liter}"
