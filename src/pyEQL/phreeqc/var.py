from typing import Any

from pyEQL._phreeqc import PY_VAR_TYPE, PY_VRESULT, PyVar


class Var:
    def __init__(self, value: Any | None = None):
        self._var = PyVar()
        self._var.var.type = PY_VAR_TYPE.TT_EMPTY
        self.value = value  # will invoke setter

    @property
    def value(self) -> Any:
        match self._var.var.type:
            case PY_VAR_TYPE.TT_EMPTY:
                return None
            case PY_VAR_TYPE.TT_ERROR:
                return self._var.var.vresult
            case PY_VAR_TYPE.TT_LONG:
                return self._var.var.lVal
            case PY_VAR_TYPE.TT_DOUBLE:
                return self._var.var.dVal
            case PY_VAR_TYPE.TT_STRING:
                return self._var.var.sVal
            case _:
                raise RuntimeError("Unknown type")

    @value.setter
    def value(self, value) -> None:
        # If we were previously holding a string, we need to free it by
        # creating a new PyVar
        if self._var.var.type == PY_VAR_TYPE.TT_STRING:
            self._var = PyVar()

        if isinstance(value, PY_VRESULT):
            self._var.var.type = PY_VAR_TYPE.TT_ERROR
            self._var.var.vresult = value
        elif isinstance(value, int):
            self._var.var.type = PY_VAR_TYPE.TT_LONG
            self._var.var.lVal = value
        elif isinstance(value, float):
            self._var.var.type = PY_VAR_TYPE.TT_DOUBLE
            self._var.var.dVal = value
        elif isinstance(value, str):
            self._var.var.type = PY_VAR_TYPE.TT_STRING
            self._var.var.sVal = value
        elif value is None:
            self._var.var.type = PY_VAR_TYPE.TT_EMPTY
        else:
            raise RuntimeError("Unknown type")
