(tutorial)=

# Quickstart

pyEQL creates a new type (`Solution` class) to represent a chemical solution.
It also comes pre-loaded with a database of diffusion coefficients, activity
correction parameters, and other data on a variety of common electrolytes.
Virtually all of the user-facing functions in `pyEQL` are accessed through the
`Solution` class.

## Creating a Solution Object

Create a Solution object by invoking the Solution class:

```{eval-rst}
.. doctest::

   >>> from pyEQL import Solution
   >>> s1 = Solution()
   >>> s1
   <pyEQL.Solution at 0x7f9d188309b0>

```

If no arguments are specified, `pyEQL` creates a 1-L solution of water at
pH 7 and 25 degC.

More usefully, you can specify solutes and bulk properties:

```{eval-rst}
.. doctest::

    >>> s2 = Solution({'Na+':'0.5 mol/kg', 'Cl-': '0.5 mol/kg'}, pH=8, temperature = '20 degC', volume='8 L')
```

See [Creating a Solution](creating.md) for more details.

## Retrieving Solution Properties

### Bulk Solution Properties

pyEQL provides a variety of methods to calculate or look up bulk properties
like temperature, ionic strength, conductivity, and density.

```python
>>> s2.volume
8.071524653929277 liter
>>> s2.density
1.0182802742389558 kilogram/liter
>>> s2.conductivity
4.083570230022633 siemens/meter
>>> s2.ionic_strength
0.500000505903012 mole/kilogram
```

### Individual Solute Properties

You can also retrieve properties for individual solutes (or the solvent, water)

```python
>>> s2.get_amount('Na+','mol/L')
0.4946847550064916 mole/liter
>>> s2.get_activity_coefficient('Na+)
0.6838526233869155
>>> s2.get_activity('Na+')
0.3419263116934578
>>> s2.get_property('Na+','transport.diffusion_coefficient')
1.1206048116287536e-05 centimeter2/second
```

See [Getting Concentrations](amounts.md) for more details.

## Units-Aware Calculations using pint

pyEQL uses [pint](https://github.com/hgrecco/pint) to perform units-aware calculations. The pint library creates
Quantity objects that contain both a magnitude and a unit.

```python
>>> from pyEQL import ureg
>>> test_qty = pyEQL.ureg('1 kg/m**3')
1.0 kilogram/meter3
>>> test_qty.magnitude
1.0
```

Many `pyEQL` methods require physical quantities to be input as strings, then these methods return pint `Quantity` objects. A string quantity must contain both a magnitude and a unit (e.g. '0.5 mol/L').
In general, pint recognizes common abbreviations and SI prefixes. Compound units must follow Python math syntax (e.g. cm\*\*2 not cm2).

See the [Converting Units](units.md) for more details.
