.. _tutorial:


Tutorial
********

pyEQL creates a new type (`Solution` class) to represent a chemical solution. 
It also comes pre-loaded with a database of diffusion coefficients, activity
correction parameters, and other data on a variety of common electrolytes.
Virtually all of the user-facing functions in pyEQL are accessed through the
`Solution` class.

Creating a Solution Object
==========================

Create a Solution object by invoking the Solution class:

.. doctest::

   >>> import pyEQL
   >>> s1 = pyEQL.Solution()
   >>> s1
   <pyEQL.pyEQL.Solution at 0x7f9d188309b0>


If no arguments are specified, pyEQL creates a 1-L solution of water at
pH 7 and 25 degC.

More usefully, you can specify solutes and bulk properties:

.. doctest::

    >>> s2 = pyEQL.Solution([['Na+','0.5 mol/kg'],['Cl-','0.5 mol/kg']],pH=8,temperature = '20 degC', volume='8 L')

Retrieving Solution Properties
==============================
    
Bulk Solution Properties
--------------------------

pyEQL provides a variety of methods to calculate or look up bulk properties
like temperature, ionic strength, conductivity, and density.

    >>> s2.get_volume()
    8.071524653929277 liter
    >>> s2.get_density()
    1.0182802742389558 kilogram/liter 
    >>> s2.get_conductivity()
    4.083570230022633 siemens/meter 
    >>> s2.get_ionic_strength()
    0.500000505903012 mole/kilogram 

Individual Solute Properties
----------------------------

You can also retrieve properties for individual solutes (or the solvent, water)

.. doctest::

    >>> s2.get_amount('Na+','mol/L')
    0.4946847550064916 mole/liter 
    >>> s2.get_activity_coefficient('Na+)
    0.6838526233869155
    >>> s2.get_activity('Na+')
    0.3419263116934578
    >>> s2.get_property('Na+','diffusion_coefficient')
    1.1206048116287536e-05 centimeter2/second

Units-Aware Calculations using pint
===================================

pyEQL uses `pint <https://github.com/hgrecco/pint>`_ to perform units-aware calculations. The pint library creates
Quantity objects that contain both a magnitude and a unit.

    >>> from pyEQL import unit
    >>> test_qty = pyEQL.unit('1 kg/m**3')
    1.0 kilogram/meter3 

Many pyEQL methods require physical quantities to be input as strings, then these methods return pint Quantity objects.
A string quantity must contain both a magnitude and a unit (e.g. '0.5 mol/L').
In general, pint recognizes common abbreviations and SI prefixes. Compound units must follow Python math syntax (e.g. cm**2 not cm2).

Pint Quantity objects have several useful attributes. They can be converted to strings:
    
    >>> str(test_qty)
    '1.0 kg/m**3'

the magnitude, units, or dimensionality can be retrieved via attributes:

    >>> test_qty.magnitude
    1.0
    >>> test_qty.units
    <UnitsContainer({'kilogram': 1.0, 'meter': -3.0})>
    >>> test_qty.dimensionality
    <UnitsContainer({'[length]': -3.0, '[mass]': 1.0})>

See the `pint documentation <http://pint.readthedocs.io/>`_ for more details on creating and manipulating Quantity objects.


Using pyEQL in your projects
============================

To access pyEQL's main features in your project all that is needed is an import statement:

    >>> import pyEQL

In order to directly create Quantity objects, you need to explicitly import the `unit` module:

    >>> from pyEQL import unit
    >>> test_qty = pyEQL.unit('1 kg/m**3')
    1.0 kilogram/meter3 

.. warning:: if you use pyEQL in conjunction with another module that also uses pint for units-aware calculations, you must convert all Quantity objects to strings before passing them to the other module, as pint cannot perform mathematical operations on units that belong to different "registries."  See the `pint documentation <http://pint.readthedocs.io/>`_ for more details.