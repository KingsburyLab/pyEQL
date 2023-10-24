---
title: 'pyEQL: A Python interface for water chemistry'
tags:
  - Python
  - environmental engineering
  - water treatment
  - desalination
  - geochemistry
  - electrolytes
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
minerals and other solutes -- are vital to understanding natural systems and developing
new technologies for water purification, wastewater treatment, and sustainable industrial
processes. `pyEQL` provides object representations for aqueous solutions, creating a
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

Available software such as `PHREEQC` [@Charlton2011], `GeoChemist's Workbench` [@gwb], or
OLI Studio `[@oli]` implement numerous electrolyte models and contain powerful capabilities for
specialists, but are not highly accessible for routien use by others, due to a
steep learning curve, difficult interoperability with other tools (such as external transport models),
the lack of a freely-available version, and/or limitation to specific operating system platforms.
Although there are `python` interfaces to the open-source `PHREEQC` software, these interfaces
are either not object-oriented (IPhreeqC [@Parkhurst2013]), poorly documented, and/or only offer access
to only a limited subset of the PHREEQC parameter databases (`phreeqpython`, `pyeqion2`). There
are even some limitations within the parent software packages. For example, if `PHREEQC` is
used with the `pitzer.dat` database (the most accurate for high salinity solutions), then
it is unable to calculate the solution's conductivity. A researcher seeking quality
data on common bulk properties such as density or viscosity or solute-specific
properties such as diffusion coefficient, transport number, or activity coefficient is thus left
to piece together outputs from disparate models and literature -- a time-consuming and error-prone process.
`pyEQL` is designed to free researchers from the tedium of identifying appropriate
models, compiling the required parameters from literature, implementing, and validating them.

It defines a python `Solution` class from which properties can be easily retrieved. It
implements the Pitzer model [@May2011b] for binary salts, with mixing rules [@Mistry2013] for
more complex solutions, and decays gracefully to more approximate models (such as the Debye-Huckel activity
model) when adequate data is not available. The built-in property database includes
Pitzer model parameters [@May2011b] for more than 100 salts, diffusion coefficients [CRCdiffusion] for more than
100 solutes, and an ever-expanding set of additional property data that make the best-available models
transparently accessible to the end user.

The package makes us of `pint` [@pint] to provide automatic unit conversions and is designed to
interoperate with widely used scientific codes in the Materials Project [@Jain2013] ecosystem --
namely, `pymatgen` [@Ong2013] for chemical informatics (e.g., molecular weight, parsing chemical
formulae) and the `maggma` `Store` API [@maggma] for accessing the built-in property database.

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

## Facilitate intgegration with other models and databases


# Architecture

## The `Solution` class

## The `Solute` class

## Modeling Engines

only get_activity_coefficient and get_solute_volume and equilibrate

composed of building blocks written in functional style

easy to compose alternative model implementations via inheritance

## Property Database

Serialized `Solute` classes, `maggma`, designed for expansion and customizability


First commit Nov 5, 2013
First release 12/10/14

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.


# Figures

TODO - need a nice figure showing pyEQL connecting to different modeling engines

Block diagram of pyEQL architecture - Solution, Solute, EOS ("modeling engine")
Maybe a ray showing common properties like density, alkalinity, osmotic pressure,
concentrations

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png){ width=20% }
and referenced from text using \autoref{fig:example}.


# Acknowledgements

The author gratefully acknowledges the Persson Research Group at the University
of California, Berkeley, especially Shyam Dwaraknath, Matthew K. Horton, Donny Winston,
Jason Munro, and Orion Cohen, for their guidance in scientific
software development practices. I also acknowledge Hernan Grecco for helpful discussions
regarding unit conversion and XXX and YYY for recent contributions. The author acknowledges
financial support from Princeton University, Membrion, Inc., and
Bluecell Energy, LLC (no longer operating).

# References

- JESS review papers? `[@rowland_ProgressAqueousSolution_2019]`
- May Pitzer compilation? `[@May2011b]`
- Effective Pitzer Model `[@Mistry2013]`
- CRC handbook `[CRCdiffusion]`
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
- pyequion2 `[@marcellos2021pyequion; @pyequion2]`
