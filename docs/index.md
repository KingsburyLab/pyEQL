---
Date: '{{ today }}'
Release: '{{ release }}'
---

![pyeql-logo](../pyeql-logo.png)

# Welcome to pyEQL's documentation!

## Description

pyEQL is a Python library that provides tools for modeling aqueous electrolyte
solutions. It allows the user to manipulate solutions as Python
objects, providing methods to populate them with solutes, calculate
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

pyEQL is designed to be customizable and easy to integrate into projects
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

Contents:

```{toctree}
:maxdepth: 2
quickstart
```

```{toctree}
:caption: User Guide
:maxdepth: 3
installation
creating
chemistry
units
amounts
engines
database
arithmetic
serialization
class_solution
```

```{toctree}
:caption: For Developers
:maxdepth: 1
contributing
utilities
internal
```
