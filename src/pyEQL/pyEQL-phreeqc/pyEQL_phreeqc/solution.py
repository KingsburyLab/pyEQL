from copy import deepcopy
from typing import Any


class Solution:
    def __init__(self, props):
        self._number = -1
        self._input_props = deepcopy(props)
        self._calculated_props = {}

    def __str__(self):
        return "\n".join(f"{k} {v}" for k, v in self._input_props.items())

    def _set_calculated_props(self, props: dict[Any, Any]):
        self._calculated_props = deepcopy(props)

    def _get_calculated_props(self):
        return self._calculated_props

    def _get_calculated_prop(self, which, species: str | None = None):
        if species is not None:
            return self._calculated_props["species"][species][which]
        return self._calculated_props[which]

    def get_activity(self, species):
        return self._get_calculated_prop("ACT", species)

    def get_molality(self, species):
        return self._get_calculated_prop("MOL", species)

    def get_osmotic_coefficient(self):
        return self._get_calculated_prop("OSMOTIC")
