# Converting Units

`pyEQL` uses [pint](https://pint.readthedocs.io/en/stable/) to automatically interpret and
convert units. For this reason, many quantitative arguments are passed to functions
**as strings** rather than numbers. For example, to specify temperature, you pass
`temperature='298 K'` and NOT `temperature=298`.

## Quantity objects

Most `Solution` class methods return `pint` `Quantity` objects.

```
>>> from pyEQL import Solution
>>> s = Solution()
>>> s.pressure
<Quantity(1, 'standard_atmosphere')>
```

If you want to create a simple `Quantity` not attached to a `Solution`, you can do
so as follows:

```
>>> from pyEQL import ureg
>>> q = ureg.Quantity('1 m')
```

`Quantity` objects have three important attributes: `magnitude`, `units`, and `dimensionality`. To get the numerical value, call magnitude

```
>>> from pyEQL import ureg
>>> q = ureg.Quantity('1 m')
>>> q.magnitude
1
```

Similarly, to get the units, call `units`

```
>>> from pyEQL import ureg
>>> q = ureg.Quantity('1 m')
>>> q.units
<Unit('meter')>
```

To convert from one unit to another, use `to()`:

```
>>> from pyEQL import ureg
>>> q = ureg.Quantity('1 m')
>>> q.to('ft')
<Quantity(3.2808399, 'foot')>
```

If you encounter a `DimensionalityError` when working with pyEQL, it probably
means you are trying to do an operation on two quantities with incompatible units (or
perhaps on a `Quantity` and a regular `float` or `int`). For example, you can't convert
`m` into `m**3`:

```
>>> from pyEQL import ureg
>>> q = ureg.Quantity('1 m')
>>> q.to('m^3')
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/home/ryan/mambaforge/envs/pbx/lib/python3.10/site-packages/pint/facets/plain/quantity.py", line 517, in to
    magnitude = self._convert_magnitude_not_inplace(other, *contexts, **ctx_kwargs)
  File "/home/ryan/mambaforge/envs/pbx/lib/python3.10/site-packages/pint/facets/plain/quantity.py", line 462, in _convert_magnitude_not_inplace
    return self._REGISTRY.convert(self._magnitude, self._units, other)
  File "/home/ryan/mambaforge/envs/pbx/lib/python3.10/site-packages/pint/facets/plain/registry.py", line 961, in convert
    return self._convert(value, src, dst, inplace)
  File "/home/ryan/mambaforge/envs/pbx/lib/python3.10/site-packages/pint/facets/context/registry.py", line 403, in _convert
    return super()._convert(value, src, dst, inplace)
  File "/home/ryan/mambaforge/envs/pbx/lib/python3.10/site-packages/pint/facets/nonmultiplicative/registry.py", line 254, in _convert
    return super()._convert(value, src, dst, inplace)
  File "/home/ryan/mambaforge/envs/pbx/lib/python3.10/site-packages/pint/facets/plain/registry.py", line 1000, in _convert
    raise DimensionalityError(src, dst, src_dim, dst_dim)
pint.errors.DimensionalityError: Cannot convert from 'meter' ([length]) to 'meter ** 3' ([length] ** 3)
```

Refer to the [pint documentation](https://pint.readthedocs.io/en/stable/) for more
etails about working with `Quantity`.

```{note}
Note that the meaning of `ureg` is equivalent in the above `pyEQL` examples and in the [pint documentation](http://pint.readthedocs.io/). `pyEQL` instantiates its own `UnitRegistry` (with custom definitions for solution chemistry) and assigns it to the variable `ureg`. In most `pint` examples, the line `ureq = UnitRegistry()` does the same thing.
```

```{important}
if you use `pyEQL` in conjunction with another module that also uses pint for units-aware calculations, you must convert all `Quantity`  objects to strings before passing them to the other module, as pint cannot perform mathematical operations on units that belong to different "registries."  See the [pint documentation](http://pint.readthedocs.io/) for more details.
```

## Custom Units

`pyEQL` extends the `pint` unit library to include some additional units that are
commonly encountered in solution chemistry.
