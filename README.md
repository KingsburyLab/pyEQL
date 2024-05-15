[![Read the Docs](https://img.shields.io/readthedocs/pyeql)](https://pyeql.readthedocs.io/en/latest/)
[![testing](https://github.com/KingsburyLab/pyeql/workflows/testing/badge.svg)](https://github.com/KingsburyLab/pyeql/actions?query=workflow%3Atesting)
[![codecov](https://codecov.io/gh/KingsburyLab/pyeql/branch/main/graph/badge.svg?token=I7RP0QML6S)](https://codecov.io/gh/KingsburyLab/pyeql)
![Supported python versions](https://img.shields.io/badge/python-3.9%20%7C%203.10%20%7C%203.11%20%7C%203.12-blue)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8332915.svg)](https://doi.org/10.5281/zenodo.8332915)
[![PyPI version](https://badge.fury.io/py/pyEQL.svg)](https://badge.fury.io/py/pyEQL)
[![status](https://joss.theoj.org/papers/bdd9e247ea9736a0fdbbd5fe12bef7a6/status.svg)](https://joss.theoj.org/papers/bdd9e247ea9736a0fdbbd5fe12bef7a6)

<img src="pyeql-logo.png" alt="pyEQL logo" style="width:600px;"/>

# A python interface for water chemistry

## Description

**The goal of `pyEQL` is to provide a stable, intuitive, easy to learn python interface
for water chemistry that can be connected to a variety of different modeling engines**

Specifically, `pyEQL` defines a `Solution` class to represent an aqueous
electrolyte solution. The `Solution` class allows the user to manipulate solutions as
Python objects, providing methods to populate them with solutes, calculate
species-specific properties (such as activity and diffusion coefficients),
and retrieve bulk properties (such as density, conductivity, or volume).

```python
>>> from pyEQL import Solution
>>> s1=Solution({"Na+":"1 mol/L", "Cl-": "1 mol/L"})
>>> s1.density
<Quantity(1.03710384, 'kilogram / liter')>
>>> s1.conductivity
<Quantity(8.09523295, 'siemens / meter')>
>>> s1.osmotic_pressure.to('atm')
<Quantity(46.7798197, 'standard_atmosphere')>
>>> s1.get_amount('Na+', 'ug/L')
<Quantity(22989769.3, 'microgram / liter')>
```

`pyEQL` also includes a number of other utilities to support water chemistry analysis,
including a **built-in property database** of diffusion coefficients, activity correction
parameters, and other data on a variety of common electrolytes.

It is designed to be customizable and easy to integrate into projects
that require modeling of chemical thermodyanmics of aqueous solutions.
It aspires to provide a flexible, extensible framework for the user, with a
high level of transparency about data sources and calculation methods.

### Key Features

- Build accurate solution properties using a minimum of inputs. Just specify
  the identity and quantity of a solute and pyEQL will do the rest.

- "Graceful Decay" from more sophisticated, data-intensive modeling approaches
  to simpler, less accurate ones depending on the amount of data supplied.

- Not limited to dilute solutions. pyEQL contains out of the box support for
  the Pitzer Model and other methods for modeling concentrated solutions.

- Built in [database](https://pyeql.readthedocs.io/en/latest/database.html) containing hundreds of model
  parameters and physicochemical properties for different ions.

- Units-aware calculations (by means of the [pint](https://github.com/hgrecco/pint) library)

### Documentation

Detailed documentation is available at [https://pyeql.readthedocs.io/](https://pyeql.readthedocs.io/)

### Dependencies

- Python 3.9+. This project will attempt to adhere to NumPy's
  [NEP 29](https://numpy.org/neps/nep-0029-deprecation_policy.html) deprecation policy
  for older version of Python.
- [pint](https://github.com/hgrecco/pint) - for units-aware calculations
- [pymatgen](https://github.com/materialsproject/pymatgen) - periodic table and chemical formula information
- [phreeqpython](https://github.com/Vitens/phreeqpython) - for PHREEQC-based speciation calculations
- [iapws](https://github.com/jjgomera/iapws/) - equations of state for water
- [monty](https://github.com/materialsvirtuallab/monty) - serialization and deserialization utilities
- [maggma](https://materialsproject.github.io/maggma/) - interface for accessing the property database
- [scipy](https://www.scipy.org/) - for certain nonlinear equation solvers

## <!-- pyscaffold-notes -->

pyEQL is licensed under LGPL.

This project has been set up using PyScaffold 4.5. For details and usage
information on PyScaffold see [https://pyscaffold.org/]().
