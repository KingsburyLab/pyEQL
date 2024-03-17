"""
Tests of the logging module
"""

import logging

from pyEQL.solution import Solution


def test_logging_kwarg(caplog):
    module_log = logging.getLogger("pyEQL")
    with caplog.at_level(logging.DEBUG, "pyEQL"):
        Solution({"Na+": "4 mol/L", "Cl-": "4 mol/L"}, log_level="debug")
        assert module_log.level == logging.DEBUG
        assert "DEBUG" in caplog.text
    caplog.clear()
    with caplog.at_level(logging.ERROR, "pyEQL"):
        Solution({"Na+": "4 mol/L", "Cl-": "4 mol/L"}, log_level="ERROR")
        assert module_log.level == logging.ERROR
        assert "DEBUG" not in caplog.text
