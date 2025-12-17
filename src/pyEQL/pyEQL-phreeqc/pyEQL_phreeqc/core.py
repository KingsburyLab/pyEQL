from inspect import cleandoc
from pathlib import Path
from textwrap import dedent, indent
from typing import Any

from pyEQL_phreeqc._bindings import PyIPhreeqc
from pyEQL_phreeqc.solution import Solution
from pyEQL_phreeqc.var import Var


class Phreeqc:
    def __init__(self, database: str = "phreeqc.dat", database_directory: Path | None = None):
        self._ext = PyIPhreeqc()

        if database_directory is None:
            database_directory = Path(__file__).parent / "database"
        self._ext.load_database(str(database_directory / database))

        self._str = ""
        self._solutions: list[Solution] = []

        # TODO: Is VAR the common denominator for most operations?
        # Here we create one and modify it in operations instead of having
        # the caller create new VARs per operation.
        self._var: Var = Var()

        self.output = PhreeqcOutput(self)

    def __len__(self):
        return len(self._solutions)

    def __getattr__(self, item) -> None:
        """Delegate attribute access to the underlying PyIPhreeqc instance."""
        if hasattr(self._ext, item):
            return getattr(self._ext, item)
        raise AttributeError(f"Phreeqc has no attribute '{item}'")

    def __call__(self, *args, **kwargs):
        self.run_string(self._str)

    def __str__(self):
        return self._str

    def accumulate(self, s: str) -> None:
        self._str += dedent(s)

    def add_solution(self, solution_dict: dict) -> None:
        solution = Solution(solution_dict)
        index = len(self)
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

    def remove_solution(self, index: int) -> Solution:
        _str = (
            cleandoc(f"""
            DELETE
              -solution {index}
        """)
            + "\n"
        )
        self.accumulate(_str)
        return self._solutions.pop(index)


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
