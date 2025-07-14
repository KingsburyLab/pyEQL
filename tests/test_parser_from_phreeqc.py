import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

import textwrap
import pytest

from pyEQL import Solution

def test_from_phreeqc_minimal_parser(tmp_path):
    """Smoke-test our new from_phreeqc() on a tiny SOLUTION block."""
    pqi = tmp_path / "sample.pqi"
    pqi.write_text(textwrap.dedent("""
        SOLUTION foo
          pH 7.5
          temp 30
          Ca+2 0.01
          SO4-2 0.01

        END
    """))
    sol = Solution.from_phreeqc(str(pqi))
    assert sol.pH == pytest.approx(7.5, rel=1e-6)
    assert sol.temperature.magnitude == pytest.approx(303.15, rel=1e-6)
    assert sol.get_amount("Ca+2", "mol/L").magnitude == pytest.approx(0.01, rel=1e-6)
    assert sol.get_amount("SO4-2", "mol/L").magnitude == pytest.approx(0.01, rel=1e-6)

def test_from_phreeqc_missing_block_raises(tmp_path):
    """If there is no SOLUTION block, we should get a ValueError."""
    pqi = tmp_path / "empty.pqi"
    pqi.write_text("THIS FILE HAS NO SOLUTION BLOCK\n")
    with pytest.raises(ValueError):
        Solution.from_phreeqc(str(pqi))
