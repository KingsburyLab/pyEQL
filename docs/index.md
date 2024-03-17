---
Date: "{{ today }}"
Release: "{{ release }}"
---

<img src="../pyeql-logo.png" alt="pyEQL logo" style="width:600px;"/>

# `pyEQL`: A python interface for water chemistry

## Description

**The goal of `pyEQL` is to provide a stable, intuitive, easy to learn python interface
for water chemistry that can be connected to a variety of different modeling engines**

Specifically, `pyEQL` defines a `Solution` class to represent an aqueous
electrolyte solution. Virtually all of the user-facing functions in `pyEQL` are accessed
through the `Solution` class. It also includes a number of other utilities to support
water chemistry analysis including a database of diffusion coefficients, activity
correction parameters, and other data on a variety of common electrolytes.

## 1-minute Tutorial

### Install

```bash
pip install pyEQL
```

### Create a Solution

```python
>>> from pyEQL import Solution
>>> s1 = Solution({'Na+':'0.5 mol/kg', 'Cl-': '0.5 mol/kg'},
                         pH=8,
                         temperature = '20 degC',
                         volume='8 L')
```

### Get properties

```python
>>> s1.density
<Quantity(1.03710384, 'kilogram / liter')>
>>> s1.conductivity
<Quantity(8.09523295, 'siemens / meter')>
>>> s1.osmotic_pressure.to('atm')
<Quantity(46.7798197, 'standard_atmosphere')>
>>> s1.get_amount('Na+', 'ug/L')
<Quantity(22989769.3, 'microgram / liter')>
```

## Key Features

`pyEQL` is designed to be customizable and easy to integrate into projects
that require modeling of chemical thermodyanmics of aqueous solutions.
It aspires to provide a flexible, extensible framework for the user, with a
high level of transparency about data sources and calculation methods.

- Build accurate solution properties using a minimum of inputs. Just specify
  the identity and quantity of a solute and pyEQL will do the rest.

- "Graceful Decay" from more sophisticated, data-intensive modeling approaches
  to simpler, less accurate ones depending on the amount of data supplied.

- Not limited to dilute solutions. pyEQL contains out of the box support for
  the Pitzer Model and other methods for modeling concentrated solutions.

- Built in [database](https://pyeql.readthedocs.io/en/latest/database.html) containing hundreds of model
  parameters and physicochemical properties for different ions.

- Customizable [modeling engine system](engines.md) that allows the `Solution` API to
  work with multiple electrolyte models.

- [Units-aware calculations](units.md) (by means of the [pint](https://github.com/hgrecco/pint) library)

## Contents:

```{toctree}
:maxdepth: 2
quickstart
examples/pyeql_demo
tutorials
```

```{toctree}
:caption: User Guide
:maxdepth: 2
installation
creating
chemistry
units
amounts
arithmetic
serialization
```

```{toctree}
:caption: Advanced Topics
:maxdepth: 2
engines
database
mixing
class_solution
internal
```

```{toctree}
:caption: For Developers
:maxdepth: 1
contributing
authors
license
```
