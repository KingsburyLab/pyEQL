from pyEQL import Solution


def test_unit_conversions():
    sol = Solution({"Ca+2": "1000 ppb"})
    ppb_value = sol.get_amount("Ca+2", "ppb").magnitude
    microgram_per_liter = sol.get_amount("Ca+2", "microgram/L").magnitude
    print(f"1000 ppb -> {ppb_value} ppb (expected: 1000)")
    print(f"1000 ppb -> {microgram_per_liter} microgram/L (expected: 1000)")

    sol = Solution({"Ca+2": "1 ppm"})
    ppm_value = sol.get_amount("Ca+2", "ppm").magnitude
    milligram_per_liter = sol.get_amount("Ca+2", "milligram/L").magnitude
    print(f"1 ppm -> {ppm_value} ppm (expected: 1)")
    print(f"1 ppm -> {milligram_per_liter} mg/L (expected: 1)")

    sol = Solution({"Ca+2": "1 ppt"})
    ppt_value = sol.get_amount("Ca+2", "ppt").magnitude
    nanogram_per_liter = sol.get_amount("Ca+2", "nanogram/L").magnitude
    print(f"1 ppt -> {ppt_value} ppt (expected: 1)")
    print(f"1 ppt -> {nanogram_per_liter} ng/L (expected: 1)")


if __name__ == "__main__":
    test_unit_conversions()
