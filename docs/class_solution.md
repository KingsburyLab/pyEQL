(class-solution)=

# The Solution Class

The `Solution` class defines a pythonic interface for **creating**, **modifying**, and **estimating properties** of
electrolyte solutions. It is the core feature of `pyEQL` and the primary user-facing class.

## Creating a solution

A `Solution` created with no arguments will default to pure water at pH=7, T=25 degrees Celsius, and 1 atm pressure.

```
>>> from pyEQL import Solution
>>> s1 = Solution()
>>> s1.pH
6.998877352386266
```

Alternatively, you can use the `autogenerate()` function to easily create common solutions like seawater:

```
>>> from pyEQL.functions import autogenerate
>>> s2 = autogenerate('seawater')
<pyEQL.solution.Solution object at 0x7f057de6b0a0>
```

You can inspect the solutes present in the solution via the `components` attribute. This comprises a dictionary of solute formula: moles, where 'moles' is the number of moles of that solute in the `Solution`. Note that the solvent (water) is present in `components`, too.

```
>>> s2.components
{'H2O': 55.34455401423017,
 'H+': 7.943282347242822e-09,
 'OH-': 8.207436858780226e-06,
 'Na+': 0.46758273714962967,
 'Mg+2': 0.052661180523467986,
 'Ca+2': 0.010251594148212318,
 'K+': 0.010177468379526856,
 'Sr+2': 9.046483353663286e-05,
 'Cl-': 0.54425785619973,
 'SO4-2': 0.028151873448454025,
 'HCO3-': 0.001712651176926199,
 'Br-': 0.0008395244921424563,
 'CO3-2': 0.00023825904349479546,
 'B(OH)4': 0.0001005389715937341,
 'F-': 6.822478260456777e-05,
 'B(OH)3': 0.0003134669156396757,
 'CO2': 9.515218476861175e-06
 }
```

To get the amount of a specific solute, use `get_amount()` and specify the units you want:

```
>>> s2.get_amount('Na+', 'g/L')
<Quantity(10.6636539, 'gram / liter')>
```

Finally, you can manually create a solution with any list of solutes, temperature, pressure, etc. that you need:

```
>>> from pyEQL import Solution
>>> s1 = Solution(solutes={'Na+':'0.5 mol/kg', 'Cl-': '0.5 mol/kg'},
                  pH=8,
                  temperature = '20 degC',
                  volume='8 L'
                  )
```

## Class reference

The remainder of this page contains detailed information on each of the methods, attributes, and properties in `Solution`. Use the sidebar on the right for easier navigation.

```{eval-rst}
.. autoclass:: pyEQL.Solution
   :members:
   :inherited-members:
   :private-members: _get_property
   :special-members: __init__
   :member-order: bysource
```
