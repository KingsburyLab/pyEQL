---
title: 'pyEQL: A Python interface for water chemistry'
tags:
  - Python
  - environmental engineering
  - water treatment
  - desalination
  - geochemistry
  - electrolytes
  - electrochemistry
authors:
  - name: Ryan Kingsbury
    orcid: 0000-0002-7168-3967
    equal-contrib: false
    affiliation: 1
affiliations:
 - name: Department of Civil and Environmental Engineering and the Andlinger Center for Energy and the Environment, Princeton University, USA
   index: 1
date: 5 November 2023
bibliography: paper.bib

---

# Summary

The properties and behavior of aqueous solutions -- that is, water containing dissolved
minerals and other solutes -- are vital to understanding natural systems and to developing
new technologies for water purification, wastewater treatment, and sustainable industrial
processes [@Stumm1993)]. `pyEQL` provides object representations for aqueous solutions, creating a
stable, intuitive, and easy to learn interface for calculating properties of solutions
and dissolved solutes. Its purpose is to save researchers time by making a variety of
different models accessible through a single interface and by aggregating hundreds of
properties and model parameters into a built-in database.

# Statement of need

Accurately predicting the thermodynamic and transport properties of complex electrolyte
solutions containing many solutes, especially at moderate to high salt concentrations
commonly encountered in water desalination and resource recovery applications, remains a
major scientific challenge [@rowland_ProgressAqueousSolution_2019]. This challenge is
compounded by the fact that the best available models (such as the Pitzer model [@May2011b])
are difficult to implement on an as-needed basis and require looking up many parameters.
Researchers and practitioners in fields such as water treatment and desalination, electrochemistry,
or environmental engineering need accurate information about electrolyte solutions to perform their work,
but are typically not specialists in solution chemistry or electrolyte thermodynamics.

Available software such as `PHREEQC` [@Charlton2011] (open source), `GeoChemist's Workbench` [@gwb], or
OLI Studio [@oli] implement numerous electrolyte models and contain powerful capabilities for
specialists, but are not highly accessible for routine use by others, due to a
steep learning curve, difficult interoperability with other tools (such as external transport models),
the lack of a freely-available version, and/or limitation to specific operating systems.
Several `python` interfaces to the open-source `PHREEQC` software exist, including IPhreeqC [@Parkhurst2013],
`phreeqpython` [@phreeqpython], and `pyeqion2` [@marcellos2021pyequion; @pyequion2]; however, these interfaces
are either not object-oriented, poorly documented, and/or only offer access to only a limited subset of the `PHREEQC`
parameter databases. There are more subtle limitations as well. For example, `phreeqpython` is unable to calculate
solution conductivity when used in conjunction with the `PHREEQC` `pitzer.dat` database (the most accurate for high
salinity solutions). A researcher seeking quality data on common bulk properties such as density or viscosity or
solute-specific properties such as diffusion coefficient, transport number, or activity coefficient is thus left
to piece together outputs from disparate models and literature -- a time-consuming and error-prone process.

`pyEQL` is designed to free researchers from the tedium of identifying and implementing the relevant models
and compiling the required parameters from literature. It defines a python `Solution` class from which properties
can be easily retrieved. It implements the Pitzer model [@May2011b] for binary salts, with mixing rules [@Mistry2013]
for more complex solutions, and decays gracefully to more approximate models (such as the Debye-Huckel activity
model [@Stumm1993)]) when adequate data is not available. The built-in property database includes Pitzer model
parameters [@May2011b] for more than 100 salts, diffusion coefficients [@CRCdiffusion] for more than 100 solutes,
and an ever-expanding set of additional property data that make the best-available models transparently accessible
to the end user.

![Overview of `pyEQL`'s architecture. Properties such as ionic strength, conductivity, and concentrations are calculated directly by `pyEQL`. Modeling engines are used to calculate non-ideal effects such as activity coefficients, while property database stores necessary parameters. The modular design of the modeling engines and property database facilitate customization.\label{fig:example}](pyEQL_overview.png){ width=80% }

# Example Use Cases

`pyEQL` may be useful to scientists and engineers in various fields broadly related to aqueous
solution chemistry. Specific use cases include, but are not limited to:

- Calculating the osmotic pressure of concentrated salt solutions
- Estimating the speciation of complex electrolyte solutions at different pH values
- Calculating the transport number of a specific ion
- Computing bulk solution characteristics such as ionic strength, alkalinity, or total
  dissolved solids, given the composition of solutes
- Converting concentrations between different unit systems, e.g. moles per L, weight %, parts per million
- Looking up properties of individual ionic species, including molecular weight, diffusion coefficient,
  ionic, hydrated, and van der Waals radii, etc.


# Design Principles

## Return the best answer possible

Recognizing that accurate modeling of complex electrolyte solutions can be difficult or even impossible,
`pyEQL` is designed to return **the best answer possible** given the data and models available. For example,
to calculate the osmotic pressure of a solution, the built-in modeling engine first attempts to use the
Pitzer model, but if parameters are not available, it reverts to a more approximate formula rather than
raising an error. To maintain transpranency, log messages (and where appropriate, warnings) are generated
throughout the codebase to document when assumptions or approximations have to be invoked or when important
model parameters are missing from the database.

## Interoperate with other scientific codes

`pyEQL` is built to be extensible, customizable, and easy to use in conjunction with widely-used scientific `python`
libraries. Specifically, it makes use of `pint` [@pint] to provide automatic unit conversions and leverages codes in
the Materials Project [@Jain2013] ecosystem -- namely, `pymatgen` [@Ong2013] for chemical informatics (e.g., molecular
weight, parsing chemical formulae) and `maggma` [@maggma] for accessing the built-in property database.

# Architecture

## The `Solution` class

The primary user-facing object in `pyEQL` is the `Solution` class. This class contains constituitive relationships
for calculating most solution properties that depend on composition, such as total dissolved solids, ionic strength,
density, conductivity, and many others (\autoref{fig:example}). Calculations that require information about non-idealities
(e.g., activity coefficients) are handled by a "modeling engine" that is stored in `Solution` as an attribute.

## Modeling Engines

Every `Solution` contains a "modeling engine" which inherits from a base class defined in `pyEQL`. Modeling
engines provide methods for calculating non-ideal thermodynamic corrections including solute activity coefficients and
molar volumes as well as performing speciation. The results of these calculations are passed back to `Solution` where they
can be transparently accessed by the user alongside other properties. This modular design facilitates connecting the `Solution`
API to multiple modeling backends or software packages. Currently, the available modeling engines include an ideal solution
approximation, a built-in implementation of the Pitzer model, and the `PHREEQC` modeling engine.

## The Property Database and `Solute` class

`pyEQL` also provides `Solute`, a `dataclass` that defines a structured schema for solute property data.
The database distributed with `pyEQL` is a list of serialized `Solute` objects stored in a `.json` file, which is
accessed via the `maggma` `Store` API [@maggma]. The database used by a particular `Solution` instance can be specified
by keyword argument when the object is created, which makes it possible in principle to use customized databases. Furthermore,
using the `Store` API means that such databases can be stored in any format supported by `maggma` (e.g., Mongo Database,
.json file, etc.).

# Acknowledgements

The author gratefully acknowledges the Persson Research Group at the University
of California, Berkeley, particularly Shyam Dwaraknath, Matthew K. Horton, Donny Winston,
Jason Munro, and Orion Cohen, for their guidance in scientific
software development practices. I also acknowledge Hernan Grecco for helpful discussions
regarding unit conversion and Kirill Pushkarev, Dhruv Duseja and Andrew Rosen for recent contributions.
The author acknowledges partial financial support from Princeton University, Membrion, Inc., and
Bluecell Energy, LLC over the period 2013-2023.

# References

<!-- - JESS review papers? `[@rowland_ProgressAqueousSolution_2019]`
- May Pitzer compilation? `[@May2011b]`
- Effective Pitzer Model `[@Mistry2013]`
- CRC handbook `[@CRCdiffusion]`
- phreeqc `[@Charlton2011]`
- iphreeqc / phreeqpy `[@Parkhurst2013]`
- phreeqpython `[@phreeqpython]`
- pint `[@pint]`
- pymatgen `[@Ong2013]`
- The Materials Project `[@Jain2013]`
- maggma `[@maggma]`
- geochemist's workbench `[@gwb]`
- OLI studio `[@oli]`
- JESS `[@marcellos2021pyequion; @pyequion2]`
- pyequion2 `[@marcellos2021pyequion; @pyequion2]` -->
