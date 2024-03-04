# Arithmetic Operations

## Addition and Subtraction

You can use the `+` operator to mix (combine) two solutions. The moles of each component
in the two solutions will be added together, and the volume of the mixed solution will
be _approximately_ equal to the sum of the two volumes, depending on the electrolyte
modeling engine used. The pressure and temperature of the mixed solution are computed
as volume-weighted averages.

```python
>>> from pyEQL import Solution
>>> s1 = Solution({"Na+": "0.5 mol/L", "Cl-": "0.5 mol/L"})
>>> s2 = Solution({"Na+": "0.1 mol/L", "Cl-": "0.1 mol/L"})
>>> s1+s2
<pyEQL.solution.Solution object at 0x7f171aee3af0>
>>> (s1+s2).get_amount('Na+', 'mol')
<Quantity(0.6, 'mole')>
>>> (s1+s2).volume
<Quantity(1.99989659, 'liter')>
```

```{note}
Both `Solution` involved in an addition operation must use the same [electrolyte
modeling engine](engines.md).
```

Subtraction is not implemented and will raise a `NotImplementedError`.

```
>>> s1-s2
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/home/ryan/mambaforge/envs/pbx/code/pyEQL/src/pyEQL/solution.py", line 2481, in __sub__
    raise NotImplementedError("Subtraction of solutions is not implemented.")
NotImplementedError: Subtraction of solutions is not implemented.
```

## Multiplication and Division

The `*` and `/` operators scale the volume and all component amounts by a factor.

```python
from pyEQL import Solution
>>> s = Solution({"Na+": "0.2 mol/L", "Cl-": "0.2 mol/L"})
>>> s.volume
<Quantity(1, 'liter')>
>>> s.get_amount('Cl-', 'mol')
<Quantity(0.2, 'mole')>
>>> s*=1.5
<Quantity(1.5, 'liter')>
s.get_amount('Cl-', 'mol')
<Quantity(0.3, 'mole')>
```

The modulo operator `//` is not implemented.
