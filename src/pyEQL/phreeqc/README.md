# `pyEQL.phreeqc`: Python bindings for IPhreeqc

This folder holds pyEQL's own Python interface to
[PHREEQC](https://www.usgs.gov/software/phreeqc-version-3), the USGS geochemical
speciation code. PHREEQC is written in C/C++, so we can't `import` it directly. Instead we
compile a small C++ **extension module** that links against PHREEQC and that Python can
load like any other module.

This document is for contributors who haven't written Python bindings for a C/C++
library before. It covers how the pieces fit together, how the build works, and how to
make the most common changes.

---

## 1. Background: how Python talks to C/C++

A Python "extension module" is a shared library (`.so` on Linux, `.pyd` on Windows,
`.so`/`.dylib` on macOS) that exports one specially named entry point. When you write
`import foo`, Python looks for `foo.py` and also for files like
`foo.cpython-312-x86_64-linux-gnu.so`. If it finds the shared library, it calls that entry
point to get a module object.

Writing that entry point by hand with the raw CPython C API is tedious and easy to get
wrong (reference counting, type conversion, error handling). We use
**[pybind11](https://pybind11.readthedocs.io/)**, a C++ library that generates
all of that glue for us. With pybind11 you write things like:

```cpp
py::class_<MyCppClass>(m, "MyPythonName")
    .def(py::init<>())                       // expose the constructor
    .def("method", &MyCppClass::method);     // expose a method
```

pybind11 then converts arguments and return values between Python and C++ (`str` and
`std::string`, `int` and `int`, `float` and `double`, and so on). It also translates C++
exceptions into Python exceptions; a `std::runtime_error` becomes a `RuntimeError`.

---

## 2. The layers

```
 pyEQL/engines.py  (Phreeqc2026EOS and its subclasses)
        │  imports
        V
 pyEQL.phreeqc            <- pure Python, this folder   (core.py, solution.py)
        │  imports
        V
 pyEQL.phreeqc._iphreeqc  <- compiled extension module  (bindings.cpp + iphreeqc_wrapper.hpp)
        │  statically links
        V
 libIPhreeqc.a            <- USGS IPhreeqc C API         (ext/iphreeqc-<version>/)
        │
        V
 PHREEQC engine           <- USGS C++ code               (ext/iphreeqc-<version>/src/phreeqcpp/)
```

Note the similar-looking names. In the same top-down order as the diagram above:

| Name | What it is                                                                                                                                                                                                                     | Where it's defined |
|---|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---|
| **`pyEQL.phreeqc`** | This folder: the pure-Python, user-facing **package** built on top of `_iphreeqc` (below). Put chemistry logic, convenience methods, and anything "Pythonic" here. This code is much easier to write, debug, and test than C++. | `__init__.py`, `core.py`, `solution.py` |
| **`Phreeqc`** | The main Python **class** in `pyEQL.phreeqc`. It holds one `IPhreeqcWrapper` as `self._wrapper` and forwards unknown attributes to it (see §4).                                                                                         | `core.py` |
| **`pyEQL.phreeqc._iphreeqc`** | The compiled extension **module**, inside the `pyEQL.phreeqc` package (leading underscore = private). It is deliberately thin: it only exposes IPhreeqc functions to Python, with no extra logic.                                                                  | built from `bindings.cpp` |
| **`IPhreeqcWrapper`** | Our class around one IPhreeqc instance. It owns the `id` and turns error codes into exceptions. Defined in C++ and exposed to Python under the same name by `bindings.cpp`.                                                    | `iphreeqc_wrapper.hpp`, `bindings.cpp` |
| **IPhreeqc** | The upstream USGS C API (`CreateIPhreeqc`, `RunString`, ...). Every function takes an integer instance `id`.                                                                                                                   | `ext/iphreeqc-<version>/src/IPhreeqc.h` |

In short: `Phreeqc` is a higher-level Python class that *contains* one
`IPhreeqcWrapper`. `_iphreeqc` and `phreeqc` are modules, not classes.

Not to be confused with **`phreeqpython`**, a separate third-party package that pyEQL also
depends on. It bundles its own copy of IPhreeqc and is unrelated to the code in this
folder.

### Files in this folder

| File | Language | Role                                                                                                                                                                                                                                                                                                                                                    |
|---|---|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `iphreeqc_wrapper.hpp` | C++ | Header-only `IPhreeqcWrapper` class, with no pybind11 dependency. IPhreeqc has a C-style API where every function takes an integer instance `id`. This class owns that `id` by creating it in the constructor and destroying it in the destructor. It also turns IPhreeqc error codes into C++ exceptions and converts C types into standard C++ types. |
| `bindings.cpp` | C++ | The pybind11 module definition (`PYBIND11_MODULE(_iphreeqc, m)`). It only decides which methods of `IPhreeqcWrapper` are visible from Python and under which names. This is the only file CMake compiles; it `#include`s the header.                                                                                                                    |
| `core.py` | Python | `Phreeqc` class. It loads a database, builds PHREEQC input strings, runs them, and parses the selected-output table into `PhreeqcSolution` objects. It also contains `PhreeqcOutput` (a class that allows numpy-style indexing into the selected output).                                                                                               |
| `solution.py` | Python | `PhreeqcSolution`, a solution with input properties and calculated results. Its API mimics the parts of `phreeqpython` that pyEQL uses.                                                                                                                                                                                                                 |
| `__init__.py` | Python | Public exports, plus `__version__`, which tracks the bundled IPhreeqc version.                                                                                                                                                                                                                                                                          |
| `ext/` | C/C++ | Vendored copy of the upstream IPhreeqc source. See `ext/README.md` for where it was downloaded from.                                                                                                                                                                                                                                                    |
| `database/` | data | **Generated at install time**; it exists only in the installed package (e.g. `site-packages/pyEQL/phreeqc/database/`), not in the source tree. CMake copies it from `ext/iphreeqc-*/database/`.                                                                                                                                                         |

### Who uses this?

`Phreeqc2026EOS` in `src/pyEQL/engines.py` imports `Phreeqc` and `PhreeqcSolution` from here. Tests live in
`tests/phreeqc/` (unit tests for this package) and `tests/test_engine_phreeqc2026.py`
(engine-level tests). `tests/phreeqc/conftest.py` skips the tests when the extension isn't
built.

---

## 3. How the build works

The build is driven by **[scikit-build-core](https://scikit-build-core.readthedocs.io/)**,
which runs **CMake** whenever pip or uv builds pyEQL.

1. `pyproject.toml` sets `build-backend = "scikit_build_core.build"` and requires
   `pybind11` at build time.
2. scikit-build-core runs the top-level `CMakeLists.txt` (in the repository root), which:
   * points `IPHREEQC_BASE` at `src/pyEQL/phreeqc/ext/iphreeqc-<version>`;
   * runs `add_subdirectory()` on it, so IPhreeqc's own CMake project builds the static
     library `libIPhreeqc.a`;
   * sets **hidden symbol visibility** on the IPhreeqc target. This matters because
     `phreeqpython` (also a pyEQL dependency) ships its own, different copy of IPhreeqc.
     If both export the same symbol names, loading both into one Python process can
     segfault. Hidden visibility keeps our copy private to `_iphreeqc.so`;
   * runs `pybind11_add_module(_iphreeqc .../bindings.cpp)` and links it against
     `IPhreeqc`;
   * installs `_iphreeqc` to `pyEQL/phreeqc/` in the wheel, next to the Python files
     (CMake component `python`);
   * installs the IPhreeqc `database/` folder to `pyEQL/phreeqc/database/` and adds an
     empty `__init__.py`, so `importlib.resources.files("pyEQL.phreeqc.database")` can find
     it (CMake component `iphreeqc_database`).
3. The result is one self-contained `_iphreeqc.cpython-*.so`. IPhreeqc is linked
   statically, so users don't need a system PHREEQC installation.

### Building without the extension

Set `PYEQL_ENABLE_EXT=OFF` in the environment (wired up in `pyproject.toml` under
`[tool.scikit-build.cmake.define]`) to skip all of the above. Then `import pyEQL.phreeqc._iphreeqc`
fails, `pyEQL.phreeqc.IS_AVAILABLE` is `False`, and the PHREEQC-based engines raise a clear
error instead of crashing.

### Release wheels (cibuildwheel)

Because `_iphreeqc` is compiled, a wheel only works on the OS, CPU architecture, and
Python version it was built for. To spare users from needing a C++ compiler, we publish a
pre-built wheel for every supported combination, using
[cibuildwheel](https://cibuildwheel.pypa.io/). For each target, it runs the normal
`scikit-build-core` build described above in a suitable environment, then installs the
wheel and tests it.

* **When:** `.github/workflows/release.yml` runs when a GitHub release is published. It
  builds the sdist and the wheels, then uploads them all to PyPI.
* **Targets:** the workflow matrix covers Linux (x86_64 and aarch64), Windows (64-bit),
  and macOS (Intel and Apple Silicon). `cibuildwheel.toml` picks the Python versions
  (`build = "cp3*"`, bounded by `requires-python` in `pyproject.toml`) and the targets to
  `skip`.
* **Tests:** each wheel runs only `tests/test_engine_phreeqc2026.py` (`test-command`), as
  a quick check that the extension loads and works. The full suite runs in the regular
  CI workflows.

To reproduce a single release build locally (on Linux this needs Docker):

```bash
uvx cibuildwheel --config-file cibuildwheel.toml --print-build-identifiers   # list targets
uvx cibuildwheel --config-file cibuildwheel.toml --only cp312-manylinux_x86_64
```

Wheels are written to `wheelhouse/`. This is the place to debug failures that occur
only in the release workflow, such as a compiler or platform difference.

---

## 4. A closer look at the C++ code

### Lifetime of an IPhreeqc instance

```cpp
IPhreeqcWrapper()  { id = CreateIPhreeqc(); ... }   // Python: IPhreeqcWrapper()
~IPhreeqcWrapper() { DestroyIPhreeqc(id); }         // Python: object garbage-collected
```

pybind11 holds the C++ object inside the Python object. When Python garbage-collects the
`IPhreeqcWrapper`, the destructor runs and frees the PHREEQC instance.

### Errors

The IPhreeqc C functions return error codes and never throw. Our wrapper checks the code
and throws `std::runtime_error(GetErrorString(id))`, and pybind11 turns that into a Python
`RuntimeError` that carries PHREEQC's own error text. Follow this pattern for any new
function that can fail.

### Selected-output values (`VAR`)

IPhreeqc returns selected-output cells as a C struct called `VAR`. This is a *tagged union*:
a `type` field says which of `lVal` / `dVal` / `sVal` holds the value, and `sVal` is a
heap-allocated C string that must be freed with `VarClear`.

We never expose `VAR` to Python. Instead, `IPhreeqcWrapper::get_value(row, col)` creates
a `VAR`, calls `VarInit`, reads the cell, copies the value out, and calls `VarClear`, all in
one function. It returns a `std::variant<std::monostate, long, double, std::string>`, which
pybind11 (via `<pybind11/stl.h>`) converts to `None` / `int` / `float` / `str`. So from
Python, `get_value(row, col)` just returns a plain value.

`std::variant` needs C++17, which is why the root `CMakeLists.txt` sets
`target_compile_features(_iphreeqc PRIVATE cxx_std_17)`.

### Attribute delegation

`Phreeqc.__getattr__` forwards unknown attributes to the underlying `IPhreeqcWrapper`. So once
a method is bound in C++, it is immediately callable on a `Phreeqc` object
(`Phreeqc().get_component_count()`) without any Python changes.

---

## 5. Development workflow

### First-time setup

From the repository root:

```bash
uv sync --extra testing        # builds _iphreeqc via CMake and installs pyEQL in editable mode
uv run pytest tests/phreeqc tests/phreeqc
```

You need a C++ compiler and CMake ≥ 3.20 on your `PATH`. (scikit-build-core can download
CMake automatically if it's missing.)

### Python-only changes (`core.py`, `solution.py`)

The install is editable, so edits take effect immediately. Just re-run your tests.

### C++ changes (`bindings.cpp`, `iphreeqc_wrapper.hpp`, `CMakeLists.txt`, `ext/`)

The compiled `.so` is **not** rebuilt automatically. It lives in
`.venv/lib/python3.*/site-packages/pyEQL/phreeqc/`, not in this folder. After changing any C++
or CMake file, force a rebuild:

```bash
uv sync --extra testing --reinstall-package pyEQL
```

If your change doesn't seem to take effect, the build is probably stale. Check the
timestamp of `.venv/lib/python3.*/site-packages/pyEQL/phreeqc/_iphreeqc*.so`.

**Tip:** the quickest way to see what the compiled module exposes:

```bash
uv run python -c "import pyEQL.phreeqc._iphreeqc as m; print(dir(m)); print(dir(m.IPhreeqcWrapper))"
```

---

## 6. How-to: expose another IPhreeqc function

IPhreeqc has far more functions than we currently wrap. They are all declared, with
documentation, in `ext/iphreeqc-<version>/src/IPhreeqc.h`. As an example, here is how to
expose PHREEQC's main text output (`SetOutputStringOn` / `GetOutputString`).

**Step 1: add a method to `IPhreeqcWrapper` in `iphreeqc_wrapper.hpp`.** Hide the `id`,
convert `const char*` to `std::string`, and turn error codes into exceptions:

```cpp
    void set_output_string_on(bool on) {
        if (SetOutputStringOn(id, on ? 1 : 0) != IPQ_OK) {
            throw std::runtime_error(get_error_string());
        }
    }

    std::string get_output_string() {
        const char* s = GetOutputString(id);
        return std::string(s ? s : "");   // never construct std::string from NULL
    }
```

**Step 2: bind it in `bindings.cpp`** by adding to the `py::class_<IPhreeqcWrapper>` chain:

```cpp
        .def("set_output_string_on", &IPhreeqcWrapper::set_output_string_on)
        .def("get_output_string", &IPhreeqcWrapper::get_output_string)
```

**Step 3: rebuild** (`uv sync --extra testing --reinstall-package pyEQL`).

**Step 4 (optional): add a Pythonic wrapper in `core.py`** if the raw call isn't
convenient enough. Because of `__getattr__` delegation.

**Step 5: add a test** in `tests/phreeqc/test_phreeqc.py`.

### Rules of thumb for C++/Python interfacing

* **Strings:** take `const std::string&` and return `std::string`, so pybind11 handles
  conversion. Never return a raw `char*` that the caller would have to free.
* **Pointers into IPhreeqc memory** (such as the `const char*` from `GetOutputString`)
  stay valid only until the next call on that instance. Copy them into a `std::string`
  right away, as shown above.
* **Out-parameters** (C functions that write into `int*` or `double*`): return the value
  instead, or a `std::tuple` / `std::vector` / `std::variant`. pybind11 converts these to
  Python objects.
* **Don't expose raw C structs.** Convert them to standard C++ types inside the wrapper,
  as `get_value` does with `VAR`. Otherwise you have to manage their memory from Python.
* **Memory you allocate** (`strdup`, `new`, `VarInit`) must be freed in the same
  function, or by an object that clearly owns it.
* **Keep the header free of pybind11.** Python-specific code belongs in `bindings.cpp`.
* **Keep logic in Python.** The C++ layer should be a thin 1:1 mapping. If you're writing
  loops or parsing in C++, consider moving that work to `core.py`.

---

## 7. How-to: upgrade the bundled IPhreeqc version

IPhreeqc source tarballs can come from two places (see `ext/README.md` for details):

* **GitHub releases** at <https://github.com/phreeqc-dev/iphreeqc/releases>, the preferred
  source, with new releases roughly every 6 months. Download the `iphreeqc-X.Y.Z-NNNNN.tar.gz`
  asset.
* **Official USGS releases** at <https://water.usgs.gov/water-resources/software/PHREEQC>,
  which come out less often. Download the "Linux (any processor)" `.tar.gz`.

Either way, the tarball contains the full, cross-platform C++ source, which is what we
build on every OS.

1. **Download and verify** the tarball, e.g. `iphreeqc-X.Y.Z-NNNNN.tar.gz`, and record its
   md5sum.
2. **Extract** it into `ext/` next to the existing version, as
   `ext/iphreeqc-X.Y.Z-NNNNN/`. Keep the upstream code **unmodified**, so future upgrades
   stay a drop-in replacement. If you need a patch, document it in `ext/README.md`.
   Stage the new folder with
   `git add -f src/pyEQL/phreeqc/ext/iphreeqc-X.Y.Z-NNNNN`. Don't stage its `build/`
   subfolder if you've already built it.
3. **Point CMake at it.** In the root `CMakeLists.txt`, update:
   ```cmake
   set(IPHREEQC_BASE ${CMAKE_CURRENT_SOURCE_DIR}/src/pyEQL/phreeqc/ext/iphreeqc-X.Y.Z-NNNNN)
   ```
   Also check that the include directories (`src/`, `src/phreeqcpp/common`) still exist
   in the new layout, and that upstream's `CMakeLists.txt` still defines a target named
   `IPhreeqc`.
4. **Bump the version** in `__init__.py` (`__version__ = "X.Y.Z-NNNNN"`).
5. **Update `ext/README.md`** with the new filename, download URL, date, and md5sum.
6. **Rebuild and test:**
   ```bash
   uv sync --extra testing --reinstall-package pyEQL
   uv run pytest tests/phreeqc tests/test_engine_phreeqc2026.py
   ```
   If compilation fails, compare the old and new `src/IPhreeqc.h` (e.g. with
   `diff -u`) for changed function signatures, then update `iphreeqc_wrapper.hpp`.
   Numerical changes in test results usually come from database updates in the new
   release (`database/*.dat`), not from the wrapper.
7. **Delete the old `ext/iphreeqc-<old>/` folder** (`git rm -r`) once everything passes,
   so the sdist doesn't ship both. (`sdist.include` in `pyproject.toml` packages all of `ext/`.)

Note: IPhreeqc's build output is written inside the source tree, in
`ext/iphreeqc-<version>/build/`. That folder is git-ignored and safe to delete whenever
you want a clean rebuild.

---

## 8. Troubleshooting

| Symptom | Likely cause |
|---|---|
| `ModuleNotFoundError: No module named 'pyEQL.phreeqc._iphreeqc'`, or `IS_AVAILABLE` is `False` | The extension wasn't built. Check for `PYEQL_ENABLE_EXT=OFF`, or a build failure hidden in pip/uv output (re-run with `-v`). |
| Your C++ change has no effect | You didn't rebuild. See §5. |
| `AttributeError: 'IPhreeqcWrapper' object has no attribute ...` | The method is in `iphreeqc_wrapper.hpp` but not bound in `bindings.cpp`, or the build is stale. |
| Segfault when `phreeqpython` and pyEQL are both imported | Symbol clash between the two IPhreeqc copies. Make sure the visibility settings in the root `CMakeLists.txt` are still applied. |
| `RuntimeError` with PHREEQC text such as "Unknown input" | A PHREEQC input error. Print the generated input with `print(phreeqc_obj)` (`Phreeqc.__str__`) and test it in standalone PHREEQC. |
| Database file not found | `database/` is created at install time. Reinstall, or pass `database_directory=` explicitly. |

## Further reading

* pybind11 documentation: <https://pybind11.readthedocs.io/>
* scikit-build-core documentation: <https://scikit-build-core.readthedocs.io/>
* IPhreeqc API reference: `ext/iphreeqc-<version>/doc/IPhreeqc.pdf` and
  `ext/iphreeqc-<version>/src/IPhreeqc.h`
* PHREEQC input format: `ext/iphreeqc-<version>/doc/Phreeqc_3_2013_manual.pdf`
