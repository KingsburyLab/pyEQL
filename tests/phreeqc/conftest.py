import pytest

from pyEQL.phreeqc import IS_AVAILABLE

if not IS_AVAILABLE:
    pytest.skip(
        "pyEQL._phreeqc extension not available",
        allow_module_level=True,
    )
