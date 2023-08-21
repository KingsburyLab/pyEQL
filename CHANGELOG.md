# pyEQL Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

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
- Add pre-commit configuration and lint with `ruff` using  rulesets mostly borrowed from `pymatgen`
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
