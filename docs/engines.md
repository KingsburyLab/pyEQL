# Electrolyte Modeling Engines

## Overview

Every `Solution` is instantiated with an electrolyte modeling "engine", which is a
subclass of [`pyEQL.engine.EOS`](internal.md#modeling-engines-module). The modeling
engine performs three functions:

1. Calculate species activity and osmotic coefficients via `get_activity_coefficient`
   and `get_osmotic_coefficient`, respectively
2. Update the composition of the `Solution` (i.e., speciation) via `equilibrate`
3. Calculate the volume occupied by the solutes via `get_solute_volume`

**All other calculations are performed directly by the `Solution` class and are the
same, regardless of the engine selected**

The purpose of this architecture is to allow `Solution` to provide a **consistent interface**
for working with electrolyte solutions, but allow the underlying models to be customized
as needed to particular use cases.

`pyEQL` currently supports three modeling engines: `ideal`, `native`, and `phreeqc`, which
are selected via the `engine` kwarg to `Solution.__init__()`. Each engine is briefly described below.

```{warning}
If you are using a Mac with an Apple M1, M2, etc. chip (i.e., Arm64 architecture), some features of `pyEQL` will be
unavailable. Specifically, anything which depends on PHREEQC (e.g., the
`equilibrate` method in the native engine and the entire
 [`phreeqc` engine](engines.md#the-phreeqc-engine)) will not work. This is because `phreeqpython` is currently
 not available for this platform. All other functions of `pyEQL` should work as expected.

Feel free to post your experiences or proposed solutions at https://github.com/KingsburyLab/pyEQL/issues/109
```

## The `'native'` engine (Default)

The `native` engine is the default choice and was the only option available prior to
version `0.6.0`.

### Activity and osmotic coefficients

Activity coefficients are calculated using the "effetive Pitzer model" of [Mistry et al.](https://doi.org/10.1016/j.desal.2013.03.015)
when possible. `pyEQL` selects parameters by identifying the predominant salt in the
solution (see [Salt Matching](amounts.md#salt-vs-solute-concentrations)). The ionic
strength is calculated based on all solutes, but only the predominant salt parameters
are used in the Pitzer calculation.

If the required parameters are not available in the [property database](database.md), the `native` engine decays gracefully
through several models more appropriate for dilute solutions, including Davies, Guntelberg,
and Debye-Huckel. See the [module reference](internal.md#modeling-engines-module) for
full details.

### Solute volumes

Solute volumes are also calculated according to the Pitzer model whenever parameters are
available. Specifically, the apparent molar volume of the primary salt is calculated
via Pitzer. The volumes of all other components (except the solvent, water) are added
based on fixed partial molar volumes, if the data are available in the [property database](database.md).
If data are not available, the volume for that solute is not accounted for.

### Speciation

Speciation calculations are provided by PHREEQC via `phreeqpython`. We use the `llnl.dat`
PHREEQC database due to its applicability for moderate salinity water and the large
number of species included (see [Lu et al.](https://doi.org/10.1016/j.earscirev.2021.103888)).
See `pyEQL.equilibrium.eqiulibrate_phreeqc` in the [module reference](internal.md#speciation-functions)
for more details.

```{warning}
Speciation support was added to the `native` engine in `v0.8.0` and should be considered
experimental. Specifically, because the `native` engine uses a non-Pitzer PHREEQC database
for speciation but uses the Pitzer model (when possible) for activity coefficients. As such,
there may be subtle thermodynamic inconsistencies between the activities and the equilibrium
concentrations returned by `equilibrate()`.
```

## The `'phreeqc'` engine

The `phreeqc` engine uses [`phreeqpython`](https://github.com/Vitens/phreeqpython)
for speciation, activity, and volume calculations. The PHREEQC engine
uses the `phreeqc.dat` PHREEQC database by default, although it is possible to instantiate
the engine with other databases such as `llnl.dat`, `pitzer.dat`, etc. See
`pyEQL.equilibrium.eqiulibrate_phreeqc` in the [module reference](internal.md#speciation-functions)
for more details.

### Activity and osmotic coefficients

Activity coefficients are calculated by dividing the PHREEQC activity by the molal
concentration of the solute.

Due to limitations in the `phreeqpython` interface, the osmotic coefficient is always
returned as 1 at present.

```{warning}
The `phreeqc` engine currently returns an osmotic coefficient of 1 and solute volume of
0 for all solutions. There appear to be limitations in the `phreeqpython` interface that
make it difficult to access these properties.
```

### Solute volumes

Due to limitations in the `phreeqpython` interface, solute volumes are ignored (as in
the `ideal` engine). More
research is needed to determine whether this is consistent with intended PHREEQC behavior
(when using the default database) or not.

```{warning}
The `phreeqc` engine currently returns an osmotic coefficient of 1 and solute volume of
0 for all solutions. There appear to be limitations in the `phreeqpython` interface that
make it difficult to access these properties.
```

### Speciation

Speciation calculations are provided by PHREEQC via `phreeqpython`.

## The `'ideal'` engine

The `'ideal'` engine applies ideal solution behavior. Activity and osmotic coefficients
are always equal to 1, solute volumes are always equal to zero, and there is no support
for speciation.

## Custom engines

The modeling engine system is designed to be extensible and customizable. To define a
custom engine, you simply need to inherit from `pyEQL.engines.EOS` (or a pre-existing
engine class) and then populate the abstract methods `get_activity_coefficient)`,
`get_osmotic_coefficient`, `get_solute_volume`, and `equilibrate`.

Equations that implement commonly used models or the above properties (such as the Debye-Huckel
and Pitzer activity models, among others) are available in `pyEQL.activity_correction` and
`pyEQL.equilibrium`, respectively. The idea is that end users can "compose" custom
engine classes by mixing and matching the desired functions from these modules, adding
custom logic as necessary.
