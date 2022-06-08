# pyEQL

![pyeql logo](pyeql-logo.png)

A Python library for solution chemistry

## Description

pyEQL is a Python library that provides tools for modeling aqueous electrolyte
solutions. It allows the user to manipulate solutions as Python
objects, providing methods to populate them with solutes, calculate
species-specific properties (such as activity and diffusion coefficients),
and retreive bulk properties (such as density, conductivity, or volume).

![pyeql demo](pyeql-demo.png)

pyEQL is designed to be customizable and easy to integrate into projects
that require modeling of chemical thermodyanmics of aqueous solutions.
It aspires to provide a flexible, extensible framework for the user, with a
high level of transparency about data sources and calculation methods.

pyEQL runs on Python 3.8+ and is licensed under LGPL.

### Key Features

- Build accurate solution properties using a minimum of inputs. Just specify
  the identity and quantity of a solute and pyEQL will do the rest.

- "Graceful Decay" from more sophisticated, data-intensive modeling approaches
  to simpler, less accurate ones depending on the amount of data supplied.

- Not limited to dilute solutions. pyEQL contains out of the box support for
  the Pitzer Model and other methods for modeling concentrated solutions.

- Extensible database system that allows one to supplement pyEQL's default
  parameters with project-specific data.

- Units-aware calculations (by means of the [pint](https://github.com/hgrecco/pint) library)

### Documentation

Detailed documentation is available at [](https://pyeql.readthedocs.io/)

### Dependencies

- Python 3.8+
- [pint](https://github.com/hgrecco/pint) - for units-awarecalculations
- [scipy](https://www.scipy.org/) - for certain nonlinear equation solvers

<!-- pyscaffold-notes -->

## Note

This project has been set up using PyScaffold 4.0.2. For details and usage
information on PyScaffold see [](https://pyscaffold.org/).
