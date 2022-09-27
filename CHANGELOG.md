# pyEQL Changelog

## 0.6.0 (in progress)

- Add more comprehensive platform testing via `tox`
- Replace `water_properties.py` with [iapws](https://github.com/jjgomera/iapws) package
- Replace elements.py with `pymatgen.core.periodic_table`
- Add `pymatgen` as a depedency
- Migrate all tests to `pytest`
- Add pre-commit configuration and lint with `black`
- Update packaging to use [pyscaffold](https://pyscaffold.org/en/stable/index.html)

## 0.5.0 (2018-09-19)

- Implement the effective Pitzer model for improved activitycalculations in multicomponent solutions
- Add support for calculation of activity and osmotic coefficientson different scales
- Add support for calculating % by weight to get_amount()
- Added methods for calculating the osmolarity or osmolality of aSolution
- Add the ability to filter list_concentrations() to show onlycations or anions
- Add two medical solutions - normal saline and Ringer's lacate -to the autogenerate method
- Add shorthand abbreviations for 'seawater' and 'wastewater' inthe autogenerate method
- Enhance automatic test suite to compare results withexperimental data based on relative error
- Add test suites for the effective Pitzer model and amulticomponent salt solution
- DEPRECATED get_mole_fraction. Use get_amount() instead
- Fix bug causing get_activity_coefficient to fail if the solute concentration was zero

## 0.4.0 (2016-07-14)

- Add ability to calculate dielectric constant based on solutioncomposition for salts
- Add database entries for the viscosity 'B' parameter for 63 moreinorganic ions
- Add domestic wastewater and human urine to the autogenerate()method
- Improve entry point for running automated tests (#16, thanksHernan Grecco)
- Significantly expand documentation of activity correction methods
- Make output of get_osmotic_coefficient more verbose when Pitzerparameters are not found
- Fix bug causing activity corrections for non 1:1 salts to becalculated incorrectly (#15)
- Fix bug causing 'bad operand type' error when calculatingosmotic pressure on some systems
- Fix bug causing ValueError exceptions when a solute has zeroconcentration
- Numerous fixes and corrections in the documentation

## 0.3.1 (2016-02-24)

- Fix packaging problems preventing installation from PyPi
- Fix character encoding issue in Erying_viscosity database file

## 0.3.0 (2016-01-15)

- Add method to calculate the total concentration of an element ina solution
- Add method to automatically generate certain solutions (likeseawater)
- Add method to calculate the hardness of a solution
- Add method to calculate the alkalinity of a solution
- Add method to calculate the charge balance of a solution
- Add method to calculate the Bjerrum length
- Add database entries for hydrated and ionic radii of 23 commonions
- Add database entries for the 'B' parameter in the Jones-Doleviscosity equation for 20 common ions
- Add test suites for solute property methods, hardness, osmoticcoefficient, and Debye length
- Improve logging system to work better when using pyEQLinteractively
- Improved README with graphics and rich formatting
- Fix bug related to activity and osmotic coefficients formultivalent salts
- Fix bug related to retrieval of water properties
- Documentation enhancements and fixes

## 0.2.2 (2015-08-28)

- Fix bug in get_amount() causing no output when mass-based units were specified.

## 0.2.1 (2015-05-06)

- Add 93 entries to diffusion coefficient database
- Add 93 entries to Pitzer partial molar volume parameters database
- Add 130 entries to Pitzer activity parameters database
- Change extension for database files from .csv to .tsv
- Corrections and additions to the contributing documentation
- Uploaded to the Python Package Index for easier installation
- Add this changelog

## 0.2.0 (2015-03-26)

- First public release
