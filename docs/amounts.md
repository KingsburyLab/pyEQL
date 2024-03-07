# Getting Concentrations

## Get the amount of a specific solute

To get the amount of a specific solute, use `get_amount()` and specify the units you want:

```python
>>> from pyEQL import Solution
>>> s = Solution({"Mg+2": "0.5 mol/L", "Cl-": "1.0 mol/L"})
>>> s.get_amount('Mg[+2]', 'mol')
<Quantity(0.5, 'mole')>
```

`get_amount` is highly flexible with respect to the types of units it can interpret. You
can request amounts in moles, mass, or equivalents (i.e., charge-weighted moles) per
unit of mass or volume.

```python
>>> s.get_amount('Mg[+2]', 'M')
<Quantity(0.5, 'molar')>
>>> s.get_amount('Mg[+2]', 'm')
<Quantity(0.506124103, 'mole / kilogram')>
>>> s.get_amount('Mg[+2]', 'eq/L')
<Quantity(1.0, 'mole / liter')>
>>> s.get_amount('Mg[+2]', 'ppm')
<Quantity(12152.5, 'milligram / liter')>
>>> s.get_amount('Mg[+2]', 'ppb')
<Quantity(12152500.0, 'microgram / liter')>
>>> s.get_amount('Mg[+2]', 'ppt')
<Quantity(1.21525e+10, 'nanogram / liter')>
```

```{important}
The unit `'ppt'` is ambiguous in the water community. To most researchers, it means
"parts per trillion" or ng/L, while to many engineers and operators it means "parts
per THOUSAND" or g/L. `pyEQL` interprets `ppt` as **parts per trillion**.
```

You can also request dimensionless concentrations as weight percent (`'%'`),
mole fraction (`'fraction'`) or the total _number_
of particles in the solution (`'count'`, useful for setting up simulation boxes).

```python
>>> s.get_amount('Mg[+2]', '%')
<Quantity(1.17358141, 'dimensionless')>
>>> s.get_amount('Mg[+2]', 'fraction')
<Quantity(0.00887519616, 'dimensionless')>
>>> s.get_amount('Mg[+2]', 'count')
<Quantity(3.01107038e+23, 'dimensionless')>
```

## See all components in the solution

You can inspect the solutes present in the solution via the `components` attribute. This comprises a dictionary of solute formula: moles, where 'moles' is the number of moles of that solute in the `Solution`. Note that the solvent (water) is present in `components`, too.
`components` is reverse sorted, with the most predominant component (i.e., the solvent)
listed first.

```python
>>> from pyEQL import Solution
>>> s = Solution({"Mg+2": "0.5 mol/L", "Cl-": "1.0 mol/L"})
>>> s.components
{'H2O(aq)': 54.83678280993063, 'Cl[-1]': 1.0, 'Mg[+2]': 0.5, 'H[+1]': 1e-07, 'OH[-1]': 1e-07}
```

Similarly, you can use the properties `anions`, `cations`, `neutrals`, and `solvent` to
retrieve subsets of `components`:

```python
>>> s.anions
{'Cl[-1]': 1.0, 'OH[-1]': 1e-07}
>>> s.cations
{'Mg[+2]': 0.5, 'H[+1]': 1e-07}
>>> s.neutrals
{'H2O(aq)': 54.83678280993063}
>>> s.solvent
'H2O(aq)'
```

Like `components`, all of the above dicts are sorted in order of decreasing amount.

## Salt vs. Solute Concentrations

Sometimes the concentration of a dissolved _salt_ (e.g., MgCl2) is of greater interest
than the concentrations of the individual solutes (Mg+2 and Cl-). `pyEQL` has the
ability to interpret a `Solution` composition and represent it as a mixture of salts.

To retrieve only _the predominant salt_ (i.e., the salt with the highest concentration),
use `get_salt`. This returns a `Salt` object with several useful attributes.

```python
>>> from pyEQL import Solution
>>> s = Solution({"Mg+2": "0.4 mol/L", "Na+": "0.1 mol/L", "Cl-": "1.0 mol/L"})
>>> s.get_salt()
<pyEQL.salt_ion_match.Salt object at 0x7f0ded09fd30>
>>> s.get_salt().formula
'MgCl2'
>>> s.get_salt().anion
'Cl[-1]'
>>> s.get_salt().z_cation
2.0
>>> s.get_salt().nu_anion
2
```

To see a `dict` of all the salts in given solution, use `get_salt_dict()`. This method
returns a dict keyed by the salt's formula, where the values are `Salt` objects converted
into dictionaries. All the usual attributes like `anion`, `z_cation` etc. are accessible
in the corresponding keys. Each value also contains a `mol` key giving the moles
of the salt present.

```python
>>> from pyEQL import Solution
>>> s = Solution({"Mg+2": "0.4 mol/L", "Na+": "0.1 mol/L", "Cl-": "1.0 mol/L"})
>>> s.get_salt_dict()
{'MgCl2': {'@module': 'pyEQL.salt_ion_match',
           '@class': 'Salt', '@version': '0.5.2',
           'cation': 'Mg[+2]',
           'anion': 'Cl[-1]',
           'mol': 0.4},
 'NaCl': {'@module': 'pyEQL.salt_ion_match',
          '@class': 'Salt', '@version': '0.5.2',
          'cation': 'Na[+1]',
          'anion': 'Cl[-1]',
          'mol': 0.1},
 'NaOH': {'@module': 'pyEQL.salt_ion_match',
          '@class': 'Salt', '@version': '0.5.2',
          'cation': 'Na[+1]',
          'anion': 'OH[-1]',
          'mol': 1e-07}
}
```

Refer to the [Salt Matching module reference](internal.md#salt-matching-module) for more
details.

## Total Element Concentrations

"Total" concentrations (i.e., concentrations of all species containing a particular
element) are important for certain types of equilibrium calculations. These can
be retrieved via `get_total_amount`. `get_total_amount` takes an element name as
the first argument, and a unit as the second.

```python
>>> from pyEQL import Solution
>>> s = Solution({"Mg+2": "0.5 mol/L", "Cl-": "1.0 mol/L"})
>>> s.equilibrate()
>>> s.components
{'H2O(aq)': 54.85346847938828, 'Cl[-1]': 0.9186683796593457, 'Mg[+2]': 0.41866839204646417, 'MgCl[+1]': 0.08133160795194606, 'OH[-1]': 1.4679440802358093e-07, 'H[+1]': 1.1833989847708719e-07, 'HCl(aq)': 1.2388705241250352e-08, 'MgOH[+1]': 3.9747494391744955e-13, 'O2(aq)': 7.027122927701743e-25, 'HClO(aq)': 1.5544872892067526e-27, 'ClO[-1]': 6.339364938003202e-28, 'H2(aq)': 5.792559717610837e-35, 'ClO2[-1]': 0.0, 'ClO3[-1]': 0.0, 'ClO4[-1]': 0.0, 'HClO2(aq)': 0.0}
>>> s.get_total_amount('Mg', 'mol')
<Quantity(0.5, 'mole')>
```

## Elements present in a `Solution`

If you just want to know the elements present in the `Solution`, use `elements`. This
returns a list of elements, sorted alphabetically.

```python
>>> from pyEQL import Solution
>>> s = Solution({"Mg+2": "0.5 mol/L", "Cl-": "1.0 mol/L"})
>>> s.elements
['Cl', 'H', 'Mg', 'O']
```
