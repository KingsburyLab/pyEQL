from copy import deepcopy
from typing import Any


class PHRQSol:
    def __init__(self, props, phreeqc=None):
        self._phreeqc = phreeqc
        self._number = -1
        self._input_props = props
        self._calculated_props = {}

    def __str__(self):
        return "\n".join(f"{k} {v}" for k, v in self._input_props.items())

    def _set_calculated_props(self, props: dict[Any, Any]):
        self._calculated_props = deepcopy(props)

    def _get_calculated_props(self):
        return self._calculated_props

    def _get_calculated_prop(self, which, species: str | None = None, eq_species: str | None = None):
        if species is not None:
            return self._calculated_props["species"][species][which]
        if eq_species is not None:
            return self._calculated_props["eq_species"][eq_species][which]
        return self._calculated_props[which]

    def get_activity(self, species) -> float:
        return self._get_calculated_prop("ACT", species=species)

    def get_molality(self, species) -> float:
        return self._get_calculated_prop("MOL", species=species)

    def get_osmotic_coefficient(self) -> float:
        return self._get_calculated_prop("OSMOTIC")

    def get_species_list(self) -> list[str]:
        return list(self._calculated_props["species"])

    def get_kgw(self) -> float:
        return self._get_calculated_prop("TOT['water']")

    def get_moles(self, species) -> float:
        return self.get_molality(species) * self.get_kgw()

    def equalize(self, phases: list[str], saturation_indices: list[float], amounts: list[float]) -> None:
        if self._phreeqc is not None:
            self._phreeqc().equalize(self._number, phases, saturation_indices, amounts)

    def si(self, eq_species) -> float:
        return self._get_calculated_prop("SI", eq_species=eq_species)

    """
    The following properties are somewhat redundant, but included in here
    so we can act as a drop-in replacement for PhreeqPython as far as its
    usage in pyEQL.
    """

    @property
    def species(self) -> list[str]:
        return self.get_species_list()

    @property
    def species_moles(self) -> dict[str, float]:
        return {species: self.get_moles(species) for species in self.species}

    def forget(self) -> None:
        if self._phreeqc is not None:
            self._phreeqc().remove_solution(self._number)
