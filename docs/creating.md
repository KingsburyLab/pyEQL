# Creating a `Solution`

The `Solution` class defines a pythonic interface for **creating**, **modifying**, and **estimating properties** of electrolyte solutions. It is the core feature of `pyEQL` and the primary user-facing class. There are several ways to create a `Solution`.


## Empty solution

With no input arguments, you get an empty `Solution` at pH 7 and 1 atm pressure.

```python
>>> from pyEQL import Solution
>>> s = Solution()
>>> print(s)
Volume: 1.000 l
Pressure: 1.000 atm
Temperature: 298.150 K
Components: ['H2O(aq)', 'H[+1]', 'OH[-1]']
```

## Manual Creation

Typically, you will create a solution by specifying a list of solutes. Solutes are
passed as a `dict` with amounts given **as strings** that include units (see [units](units.md)). Any unit that can be understood by `get_amount` is valid.

```python
>>> from pyEQL import Solution
>>> s = Solution({"Na+": "0.5 mol/L", "Cl-": "0.5 mol/L"})
```

You can also specify conditions such as temperature, pressure, pH, and pE (redox potential).


Finally, you can manually create a solution with any list of solutes, temperature, pressure, etc. that you need:

```python
>>> from pyEQL import Solution
>>> s1 = Solution(solutes={'Na+':'0.5 mol/kg', 'Cl-': '0.5 mol/kg'},
                  pH=8,
                  temperature = '20 degC',
                  volume='8 L',
                  pE = 4,
                  )
```

## Using a preset

Alternatively, you can use the `Solution.from_preset()` classmethod to easily create common solutions like seawater:

```
>>> from pyEQL import Solution
>>> s2 = Solution.from_preset('seawater')
<pyEQL.solution.Solution object at 0x7f057de6b0a0>
```

## From a dictionary

If you have [converted a `Solution` to a `dict`](serialization.md#serialization-to-dict),
you can re-instantiate it using the `Solution.from_dict()` class method.

## From a file

If you [save a `Solution` to a `.json` file](serialization.md#saving-to-a-json-file),
you can [recreate it](serialization.md#loading-from-a-json-file) using
[monty.serialization.loadfn](https://pythonhosted.org/monty/monty.html#module-monty.serialization)

```python
>>> from monty.serialization import loadfn
>>> s = loadfn('test.json')
print(s)
Volume: 1.000 l
Pressure: 1.000 atm
Temperature: 298.150 K
Components: ['H2O(aq)', 'H[+1]', 'OH[-1]']
```
