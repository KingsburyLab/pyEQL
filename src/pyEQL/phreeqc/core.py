from importlib import resources
from inspect import cleandoc
from pathlib import Path
from textwrap import dedent, indent
from typing import Any
from weakref import ref

from pyEQL.phreeqc.solution import PHRQSol

SOLUTION_PROPS = (
    "CELL_NO",
    "TOT['water']",
    "OSMOTIC",
)
SPECIES_PROPS = ("MOL", "ACT", "DIFF_C")
EQ_SPECIES_PROPS = ("SI",)


def ext_module():
    """
    Return the compiled extension module is available, else None.
    """
    try:
        from pyEQL._phreeqc import PyIPhreeqc  # noqa: PLC0415
    except ModuleNotFoundError:
        return None
    else:
        return PyIPhreeqc


IS_AVAILABLE = ext_module() is not None


class Phreeqc:
    def __init__(self, database: str = "phreeqc.dat", database_directory: Path | None = None):
        self._ext = ext_module()

        if database_directory is None:
            database_directory = resources.files("pyEQL.phreeqc.database")
        self._ext.load_database(str(database_directory / database))

        self._str = ""
        self._solutions: list[PHRQSol] = []

        from pyEQL.phreeqc.var import Var  # noqa: PLC0415

        # TODO: Is VAR the common denominator for most operations?
        # Here we create one and modify it in operations instead of having
        # the caller create new VARs per operation.
        self._var: Var = Var()

        self.output = PhreeqcOutput(self)

    def __len__(self):
        return len(self._solutions)

    def __getitem__(self, item):
        return self._solutions[item]

    def __getattr__(self, item) -> None:
        """Delegate attribute access to the underlying PyIPhreeqc instance."""
        if hasattr(self._ext, item):
            return getattr(self._ext, item)
        raise AttributeError(f"Phreeqc has no attribute '{item}'")

    def __call__(self, *args, **kwargs):
        self.run_string(self._str)

    def __str__(self):
        return cleandoc(self._str)

    def clear(self):
        self._str = ""

    def accumulate(self, s: str) -> None:
        self._str += dedent(s)

    def _add_solution(self, solution: PHRQSol | list[PHRQSol]) -> PHRQSol | list[PHRQSol]:
        singleton = isinstance(solution, PHRQSol)
        solutions = [solution] if singleton else solution

        for solution in solutions:
            index = len(self)
            solution._phreeqc = ref(self)
            solution._number = index

            # TODO: This should go in the PHRQSol class
            template = (
                "\n"
                + cleandoc(f"""
                SOLUTION {index}
                {{solution}}
                SAVE SOLUTION {index}
                END
            """)
                + "\n"
            )
            _str = template.format(solution=indent(str(solution), "  "))
            self.accumulate(_str)
            self._solutions.append(solution)

        return self._solutions[-1] if singleton else self._solutions

    def remove_solution(self, index: int) -> PHRQSol:
        _str = (
            cleandoc(f"""
            DELETE
              -solution {index}
        """)
            + "\n"
        )
        self.accumulate(_str)
        self()
        return self._solutions.pop(index)

    def add_solution(self, solution: PHRQSol | list[PHRQSol]) -> PHRQSol | list[PHRQSol]:
        solution_punch_line = ", ".join(list(SOLUTION_PROPS))
        species_punch_line = ", ".join([f"{prop}(name$(i))" for prop in SPECIES_PROPS])
        eq_species_punch_line = ", ".join([f"{prop}(name$(j))" for prop in EQ_SPECIES_PROPS])

        self.clear()
        self.accumulate(f"""
            SELECTED_OUTPUT
                -reset false

            USER_PUNCH
            5 PUNCH {solution_punch_line}, EOL_NOTAB$
            10 t = SYS("aq", count, name$, type$, moles)
            20 FOR i = 1 to count
            30 PUNCH name$(i), {species_punch_line}
            40 NEXT i
            50 PUNCH EOL$
            60 p = SYS("phases", count, name$, type$, moles)
            70 FOR j = 1 TO count
            80 PUNCH name$(j), {eq_species_punch_line}
            90 NEXT j
            """)

        return_value = self._add_solution(solution)
        self()
        self._parse_output()

        return return_value

    def _parse_output(self):
        # first line is always the header, but ensure this so we can skip it.
        assert self.output[0][0].startswith("no_heading")

        output = self.output[:][1:]
        for line in output:
            solution_i_props = {"species": {}, "eq_species": {}}
            n_tokens = len(line)
            j = 0

            while j < len(SOLUTION_PROPS):
                solution_i_props[SOLUTION_PROPS[j]] = line[j]
                j += 1

            j += 1  # skip "\n"

            while line[j] != "\n":
                species = line[j]
                solution_i_props["species"][species] = {}
                j += 1
                for prop in SPECIES_PROPS:
                    solution_i_props["species"][species][prop] = line[j]
                    j += 1

            j += 1  # skip "\n"
            while j < n_tokens:
                # We encounter None values here, indicating end of valid
                # entries.
                if line[j] is None:
                    break
                eq_species = line[j]
                solution_i_props["eq_species"][eq_species] = {}
                j += 1
                for prop in EQ_SPECIES_PROPS:
                    solution_i_props["eq_species"][eq_species][prop] = line[j]
                    j += 1

            solution_i = int(solution_i_props["CELL_NO"])
            self[solution_i]._set_calculated_props(solution_i_props)

    def equalize(self, index: int, phases: list[str], saturation_indices: list[float], amounts: list[float]) -> None:
        phases_lines = "\n".join(
            f"{p} {si} {amount}" for (p, si, amount) in zip(phases, saturation_indices, amounts, strict=False)
        )
        self.accumulate(f"""
            USE SOLUTION {index}
            EQUILIBRIUM PHASES 1
            {phases_lines}
            END
            """)
        self()
        self._parse_output()


class PhreeqcOutput:
    def __init__(self, phreeqc: Phreeqc):
        self._phreeqc = phreeqc

    def __getitem__(self, item) -> Any:
        if not isinstance(item, tuple):
            item = (item,)
        while len(item) < 2:
            item += (slice(None),)

        row_idx, col_idx = item

        if isinstance(row_idx, slice):
            row_indices = range(*row_idx.indices(self.shape[0]))
        elif isinstance(row_idx, int):
            row_indices = [row_idx]
        else:
            raise TypeError("Row index must be int or slice")

        if isinstance(col_idx, slice):
            col_indices = range(*col_idx.indices(self.shape[1]))
        elif isinstance(col_idx, int):
            col_indices = [col_idx]
        else:
            raise TypeError("Column index must be int or slice")

        result = []
        for row in row_indices:
            row_values = []
            for col in col_indices:
                self._phreeqc._ext.get_value(row, col, self._phreeqc._var._var.var)
                row_values.append(self._phreeqc._var.value)
            result.append(row_values if len(col_indices) > 1 else row_values[0])

        if len(row_indices) == 1:
            return result[0]
        return result

    @property
    def shape(self) -> tuple[int, int]:
        return self._phreeqc.get_selected_output_row_count(), self._phreeqc.get_selected_output_column_count()
