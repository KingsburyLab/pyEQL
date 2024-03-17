# pyEQL Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2024-03-17

### Added

- `Solution.__init__`: new keyword argument `log_level` allows user to control the verbosity of log messages by setting
  the level (e.g., ERROR, WARNING, etc.) that will be shown in stdout.
- Docs: added a note about a workaround for Apple M1/M2 Macs proposed by @xiaoxiaozhu123

### Changed

- `Solution.__init__`: The deprecated format for specifying solutes (e.g., `[["Na+", "0.5 mol/L]]`)
  which previously only generated log warning message, now raises a `DeprecationWarning`. Use dict-style input (e.g.,
  `{"Na+":"0.5 mol/L"}`) instead.
- New logo! Updated the `pyEQL` logo (for the first time in 9 years!) to address an obsolete font in the .svg
  and modernize the design.

### Removed

- **BREAKING** All methods and functions (with the exception of `Solution.list_XXX` methods) previously marked with
  deprecation warnings have been removed.

## [0.15.1] - 2024-03-13

### Fixed

- `Solution.get_total_amount`: Fixed an issue in which `ppm` units would fail.

## [0.15.0] - 2024-03-13

### Added

- `utils.interpret_units`: New method to "sanitize" environmental engineering units like ppm to strings that `pint`
  can understand, e.g. ppm -> mg/L. This method is now used in `get_amount` and `get_total_amount` to ensure consistency
  in how they process units.

### Changed

- CI: `pre-commit autoupdate`

### Fixed

- `Solution`: Fixed an issue in which repeated calls to `equilibrate` when using `NativeEOS` or `PHREEQCEOS` would
  change the mass of the `Solution` slightly. This was attributed to the fact that `pyEQL` and `PHREEQC` use slightly
  different molecular weights for water.
- `Solution`: `get_total_amount` and related methods could fail when the oxidation state of an element was
  unknown (e.g., 'Br') (Issue [#116](https://github.com/KingsburyLab/pyEQL/issues/116))

## [0.14.0] - 2024-03-05

### Added

- `NativeEOS` / `PhreeqcEOS`: Added `try`/`catch` so that `pyEQL` can still be used on platforms that PHREEQC does
  not support, such as Apple Silicon. In such cases, functions like `equilibrate` that depend on PHREEQC will
  raise errors, but everything else can still be used.
- CI: Added Apple M1 runner (GitHub: `macos-14`) to the CI tests.

### Fixed

- CI: Addressed several issues in the testing configuration which had resulted in testing
  fewer operating systems x python version combinations than intended. CI tests now
  correctly and comprehensively test every supported version of python on every os
  (macos, windows, ubuntu).
- `utils.FormulaDict`: implemented `__contains__` so that `get()` works correctly in
  python 3.12+. See https://github.com/python/cpython/issues/105524
- Docs: fixed many small problems in documentation causing equations and examples to
  render incorrectly.
- `Solution.from_file`: Add missing `@classmethod` decorator; update documentation.

## [0.13.0] - 2024-03-05

### Fixed

- `equilibrium.alpha()`: Fixed incorrect calculation of acid-base distribution coefficient for multiprotic acids.

## [0.12.2] - 2024-02-25

### Fixed

- `Solution.get_total_amount`: Fix bug that would cause the method to fail if
  units with names not natively understood by `pint` (e.g., 'ppm') were passed.

## [0.12.1] - 2024-02-25

### Fixed

- `Solute.from_formula`: Fix bug in which an uncaught exception could occur when
  if `pymatgen` failed to guess the oxidation state of a solute. (Issue #103 - thanks to @xiaoxiaozhu123 for reporting).
- `Solution.get_total_amount`: Fix bug that would cause the method to fail if
  mass-based units such as mg/L or ppm were requested.

## [0.12.0] - 2024-02-15

### Added

- `Solution`: new kwarg `default_diffusion_coeff` which allows the user to specify a value to use
  when a species diffusion coefficient is missing from the database. By default, the value for NaCl
  salt (1.61e-9 m2/s) is used. This is important for conductivity and transport number calculations,
  which perform weighted summations of diffusion coefficients over every species in the solution.
  Previously, species with missing diffusion coefficients would be excluded from such calculations,
  possibly resulting in inaccuracies, particularly after calling `equilibrate()`, which often
  generates charged complexes such as NaSO4- or MgCl+.

### Fixed

- `Solution.__add__`: Bugfix in the addition operation `+` that could cause problems with
  child classes (i.e., classes that inherit from `Solution`) to work improperly

### Changed

- Removed deprecated `pkg_resources` import in favor of `importlib.resources`

## [0.11.1] - 2023-12-23

### Added

- Add tests for `gibbs_mix` and `entropy_mix` functions. Format docstrings in Google style.

### Fixed

- `Solution.from_preset`: Fixed a packaging error that made this method fail with a `FileNotFoundError`.

### Removed

- `functions.py` is no longer imported into the root namespace. You'll now need to say `from pyEQL.functions import gibbs_mix`
  instead of `from pyEQL import gibbs_mix`

## [0.11.0] - 2023-11-20

### Changed

- `PhreeqcEOS`: performance improvements for the `phreeqc` engine. The `EOS` instance now retains
  the `phreeqpython` solution object in between calls, only re-initializing it if the composition
  of the `Solution` has changed since the previous call.

### Fixed

- `equilibrate`: Fixed several bugs affecting `NativeEOS` and `PhreeqcEOS` in which calling `equilibrate()`
  would mess up the charge balance. This was especially an issue if `balance_charge` was set to something
  other than `pH`.

### Removed

- `equilibrium.equilibrate_phreeqc()` has been removed to reduce redundant code. All its
  was absorbed into `NativeEOS` and `PhreeqcEOS`

## [0.10.1] - 2023-11-12

### Added

- utility function `create_water_substance` with caching to speed up access to IAPWS instances

### Changed

- `Solution.get_diffusion_coefficient`: the default diffusion coefficient (returned when D for a solute is not found in
  the database) is now adjusted for temperature and ionic strength.
- `Solution.water_substance` - use the IAPWS97 model instead of IAPWS95 whenever possible, for a substantial speedup.

## [0.10.0] - 2023-11-12

### Added

- `Solution`: Revamped docstrings for `conductivity`, `get_transport_number`, `get_molar_conductivity`, and
  `get_diffusion_coefficient`.
- `Solution`: new method `get_diffusion_coefficient` for dedicated retrieval of diffusion coefficients. This method
  implements an improved algorithm for temperature adjustment and a new algorithm for adjusting infinite dilution D values
  for ionic strengthe effects. The algorithm is identical to that implemented in PHREEQC >= 3.4.
- Database: empirical parameters for temperature and ionic strength adjustment of diffusion coefficients for 15 solutes
- Added tests for temperature and ionic strength adjustment and conductivity
- Docs: new tutorial notebooks
- Docs: remove duplicate contributing pages (Closes [#68](https://github.com/KingsburyLab/pyEQL/issues/68))
- `Solution`: new method `to_file()` for more convenient saving Solution object to json or yaml files. (@kirill-push)
- `Solution`: new method `from_file()` for more convenient loading Solution object from json or yaml files. (@kirill-push)
- `Solution`: new classmethod `from_preset()` to `replace pyEQL.functions.autogenerate()` and instantiate a solution from a preset composition. (@kirill-push)

### Changed

- `Solution`: method af adjusting diffusion coefficients for temperature was updated (same as used in PHREEQC >= 3.4)
- `Solution.conductvity`: improved equation (same as used in PHREEQC >= 3.4) which is more accurate at higher concentrations

### Fixed

- Database errors with `Cs[+1]` diffusion coefficient and `KBr` Pitzer parameters
- Restored filter that suppresses duplicate log messages

### Deprecated

- `replace pyEQL.functions.autogenerate()` is now deprecated. Use `from_preset` instead.

### Removed

- The `activity_correction` kwarg in `get_transport_number` has been removed, because this now occurs by default and is
  handled in `get_diffusion_coefficient`.

## [0.9.2] - 2023-11-07

### Fixed

- Restored filter that suppresses duplicate log messages

## [0.9.1] - 2023-11-04

### Added

- `format_solutes_dict()` method added into the utils module to help format solutes dictionaries with a unit. (@DhruvDuseja)

### Changed

- Native property database is now instantiated on `pyEQL.__init__` to speed subsequent access by `Solution`
- Numerous performance optimization increased the speed of many `Solution` property and method calls by 3x-10x

## [0.9.0] - 2023-10-17

### Added

- `Solution.print()` added to take the place of the deprecated `list_xxx` methods.

### Changed

- **Breaking** `pyEQL.unit` deprecation machinery has been removed to quiet the warning
  messages on import. The `pyEQL` unit registry was renamed to `pyEQL.ureg` in v0.6.1.
- Significant documentation updates.
- `Solution.components` is now automatically sorted in descending order of amount, for
  consistency with `anions`, `cations`, and `neutrals`.

### Removed

- `Solution.list_solutes()` has been removed. Use `Solution.components.keys()` instead.

### Fixed

- Bugfix in `as_dict` to make serialization via `dumpfn` possible. Previously, `Quantity`
  were not converted to a serializable form. Now, `Quantity` are converted to `str` in
  `as_dict()`.

### Removed

- python 3.8 is no longer supported

## [0.8.1] - 2023-10-01

### Changed

- `from_dict` modified to avoid call to `super()`, making for more robust behavior if `Solution` is inherited.

### Removed

- `copy()` method was removed for consistency with `python` conventions (it returned a deep rather than a
  shallow copy). Use `copy.deepcopy(Solution)` instead.

### Fixed

- Bugfix in `as_dict` in which the `solutes` attribute was saved with `Quantity` rather than `float`
- Simplified `Solution.get_conductivity` to avoid errors in selected cases.
- Required `pymatgen` version was incorrectly set at `2022.8.10` when it should be `2023.8.10`
- Bug in `get_osmotic_coefficient` that caused a `ZeroDivisionError` with an empty solution.

## [0.8.0] - 2023-09-27

### Added

- New electrolyte engine `PhreeqcEOS` provides `phreeqpython` activities within `pyEQL`
- `Solution`: use total element concentrations when performing salt matching (can be disabled via kwarg)
- `Solution`: add speciation support to the native engine via `phreeqpython`
- `Solution`: add keyword argument to enable automatic charge balancing
- `Salt`: class is now MSONable (i.e., serializable via `as_dict` / `from_dict`)
- `Solution`: new properties `anions`, `cations`, `neutrals` provide easy access to subsets `components`.
- `Solution`: improvements to `get_total_amount`.
- `Solution`: new function `get_components_by_element` that lists all species associated with elements in specific
  oxidation states.
- `Solution`: new properties `elements` and `chemical_system`, new function `get_el_amt_dict` to compute the total
  number of moles of each element present in the Solution.

### Changed

- `pH` attribute is now calculated from the H+ concentration rather than its activity. For the old behavior,
  use `Solution.p('H+')` which defaults to applying the activity correction.
- Update `test_salt_ion_match` to `pytest` format and add additional tests
- Update `test_bulk_properties` to `pytest` format

### Removed

- `Solution.list_salts()` has been removed. See `Solution.get_salt_dict()` for equivalent functionality
- `salt_ion_match.generate_salt_list` and `identify_salt` have been removed. See `Solution.get_salt_dict()` and `Solution.get_salt()` for equivalent functionality.

### Fixed

- Bug in `get_transport_number` caused by migration to standardized solute formulas
- Scaling of salt concentrations in `get_salt_dict` was incorrect in some edge cases
- Disable hydrate notation in `standardize_formula`, which caused hydroxides such as 'Ca(OH)3' to be written 'CaO2H.H2O'
- Inconsistent formatting of oxidation states in `get_total_amount` and `Solute`
- Inconsistent return type from `get_property` when `molar_volume` and `diffusion_coefficient` were missing
- Two issues with the formatting of the `H2O(aq)` entry in the database, `pyeql_db.json`

## [0.7.0] - 2023-08-22

### Changed

- `Solution` now more robustly converts any user-supplied formulas into unique values using `pymatgen.core.ion.Ion.reduced_formula`. This means that the `.components` or `solvent` attributes may now differ slightly from whatever is entered during `__init__`. For example, `Solution(solvent='H2O').solvent` gives `H2O(aq)`. This behavior resolved a small bug that could occur when mixing solutions. User supplied formulas passed to `get_amount` or `Solution.components[xxx]` can still be any valid formula. E.g., `Solution.components["Na+"]`, `Solution.components["Na+1"]`, and `Solution.components["Na[+]"]` will all return the same thing.

## [0.6.1] - 2023-08-22

### Added

- `Solution`: enable passing an `EOS` instance to the `engine` kwarg.
- `Solution`: new properties `total_dissolved_solids` and alias `TDS`
- `Solution`: support new units in `get_amount` - ppm, ppb, eq/L, etc.
- `Solution`: implemented arithmetic operations `+` (for mixing two solutions), `*` and `\` for scaling their amounts

### Changed

- `pyEQL.unit` was renamed to `pyEQL.ureg` (short for `UnitRegistry`) for consistency with the `pint` documentation and tutorials.

## [v0.6.0] - 2023-08-15

### Added

- `Solution`: add tests to confirm that solution density changes with temperature and pressure
- `Solution`: add tests for `charge_balance`, `alkalinity`, `hardness`, `osmotic_pressure`, `p()`, and `conductivity`
- `Solution`: `pE` attribute and kwarg
- `Solution`: add support for passing solutes as a `dict`
- Implement extensible system for connecting `Solution` to various activity and speciation
  models. Models can be integrated into pyEQL by implementing an `EOS` class. The desired
  activity model is selected on init. Currently available models are `native` (for pyEQL's
  implementation of Pitzer, which decays gracefully into Debye-Huckel and other models if
  parameters are not available) or `ideal` for a dummy engine that returns unit activity
  coefficients. Support for additional external engines such as [`phreeqpython`](https://github.com/Vitens/phreeqpython) is planned.
- Add `pymatgen`, `monty`, and `maggma` as dependencies
- Add `pre-commit` configuration
- Add pull request template, new GitHub actions, and `tox -e autodocs` environment to serve and update docs in real time
- Add pre-commit configuration and lint with `ruff` using rulesets mostly borrowed from `pymatgen`
- Add more comprehensive platform testing via `tox`

### Changed

- Complete overhaul of the property database. The database is now distributed in a .json file containing serialize `Solute` objects. `Solution` can now be connected to this database (by default) or to any other `maggma` `Store` containing properly formatted data. `database.py`, `parameter.py`, and all the `.tsv` data files have been removed and replaced with `pyeql_db.json`.
- Docs: update, change theme, convert to .md format, and adopt Keep a Changelog format
- Replace `water_properties.py` with [iapws](https://github.com/jjgomera/iapws) package
- Replace `elements.py`` with `pymatgen.core.periodic_table`
- `Solution.charge_balance` now returns in equivalents instead of Coulombs
- Migrate all tests to `pytest`
- Update packaging format to use [pyscaffold](https://pyscaffold.org/en/stable/index.html)

### Deprecated

- `Solution`: new properties `pressure`, `temperature`, `volume`, `pH`, `mass`, `density`, `viscosity_dynamic`, `viscosity_kinematic`, `ionic_strength`, `conductivity`, `debye_length`, `bjerrum_length`, `alkalinity`, `hardness`, `dielectric_constant`, `osmotic_pressure`, `solvent_mass`, `charge_balance` have replaced the corresponding get_XXX and set_XXX (for temperature and pressure) methods, which will be removed in a future release. `get_viscosity_relative` will be removed entirely.
- `Solute`: methods `get_formal_charge()`, `get_name()`, and `get_molecular_weight()` have been
  replaced by direct access to the attributes `charge`, `formula`, and `mw`, respectively.

### Removed

- disable 'verbose' kwarg in `get_activity` and `get_activity_coefficient`

### Fixed

- Fixed various documentation rendering issues
- bug in `alkalinity`

## [0.5.2] - 2020-04-21

- Fix breaking bug introduced by upstream pint change to avogadro_number
- Format project with black
- Misc. linting and docstring changes

## [0.5.0] 2018-09-19

- Implement the effective Pitzer model for improved activity calculations in multicomponent solutions
- Add support for calculation of activity and osmotic coefficients on different scales
- Add support for calculating % by weight to get_amount()
- Added methods for calculating the osmolarity or osmolality of a Solution
- Add the ability to filter list_concentrations() to show only cations or anions
- Add two medical solutions - normal saline and Ringer's lacate -to the autogenerate method
- Add shorthand abbreviations for 'seawater' and 'wastewater' in the autogenerate method
- Enhance automatic test suite to compare results with experimental data based on relative error
- Add test suites for the effective Pitzer model and a multicomponent salt solution
- DEPRECATED get_mole_fraction. Use get_amount() instead
- Fix bug causing get_activity_coefficient to fail if the solute concentration was zero

## [0.4.0] 2016-07-14

- Add ability to calculate dielectric constant based on solution composition for salts
- Add database entries for the viscosity 'B' parameter for 63 more inorganic ions
- Add domestic wastewater and human urine to the autogenerate()method
- Improve entry point for running automated tests (#16, thanks Hernan Grecco)
- Significantly expand documentation of activity correction methods
- Make output of get_osmotic_coefficient more verbose when Pitzer parameters are not found
- Fix bug causing activity corrections for non 1:1 salts to be calculated incorrectly (#15)
- Fix bug causing 'bad operand type' error when calculating osmotic pressure on some systems
- Fix bug causing ValueError exceptions when a solute has zero concentration
- Numerous fixes and corrections in the documentation

## [0.3.1] 2016-02-24

- Fix packaging problems preventing installation from PyPi
- Fix character encoding issue in Erying_viscosity database file

## [0.3.0] 2016-01-15

- Add method to calculate the total concentration of an element in a solution
- Add method to automatically generate certain solutions (like seawater)
- Add method to calculate the hardness of a solution
- Add method to calculate the alkalinity of a solution
- Add method to calculate the charge balance of a solution
- Add method to calculate the Bjerrum length
- Add database entries for hydrated and ionic radii of 23 commonions
- Add database entries for the 'B' parameter in the Jones-Dole viscosity equation for 20 common ions
- Add test suites for solute property methods, hardness, osmotic coefficient, and Debye length
- Improve logging system to work better when using pyEQL interactively
- Improved README with graphics and rich formatting
- Fix bug related to activity and osmotic coefficients for multivalent salts
- Fix bug related to retrieval of water properties
- Documentation enhancements and fixes

## [0.2.2] 2015-08-28

- Fix bug in get_amount() causing no output when mass-based units were specified.

## [0.2.1] 2015-05-06

- Add 93 entries to diffusion coefficient database
- Add 93 entries to Pitzer partial molar volume parameters database
- Add 130 entries to Pitzer activity parameters database
- Change extension for database files from .csv to .tsv
- Corrections and additions to the contributing documentation
- Uploaded to the Python Package Index for easier installation
- Add this changelog

## [0.2.0] 2015-03-26

- First public release
