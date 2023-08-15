(chemistry)=

# Chemical Formulas

pyEQL interprets the chemical formula of a substance to calculate its molecular
weight and formal charge. The formula is also used as a key to search the
database for parameters (e.g. diffusion coefficient) that are used in subsequent
calculations.

## How to Enter Valid Chemical Formulas

Generally speaking, type the chemical formula of your solute the "normal" way
and `pyEQL` should be able to inerpret it. Internally, `pyEQL` uses the [`pymatgen.core.ion.Ion`](https://pymatgen.org/pymatgen.core.html#pymatgen.core.ion.Ion)
 class to "translate" chemical formulas into a consistent format. Anything that the `Ion` class can understand will
 be processed into a valid formula by `pyEQL`.

Here are some examples:

| Substance | You enter | `pyEQL` understands |
| :--- | :---: | :---: |
| Sodium Chloride | "NaCl", "NaCl(aq)", or "ClNa" | "NaCl(aq)" | 
| Sodium Sulfate | "Na(SO4)2" or "NaS2O8" | "Na(SO4)2(aq)" | 
| Sodium Ion | "Na+", "Na+1", "Na1+", or "Na[+]" | "Na[+1]" | 
| Magnesium Ion | "Mg+2", "Mg++", or "Mg[++]" | "Mg[+2]" | 
| Methanol | "CH3OH", "CH4O" | "'CH3OH(aq)'" | 

Specifically, `pyEQL` uses `Ion.from_formula(<formula>).reduced_formla` (shown in the right hand column of the table) to
identify solutes. Notice that for charged species, the charges are always placed inside square brackets (e.g., `Na[+1]`)
and always include the charge number (even for monovalent ions). Uncharged species are always suffixed by `(aq)` to
disambiguate them from solids.

:::{important}
**When writing multivalent ion formulas, it is strongly recommended that you put the charge number AFTER the + or - 
sign** (e.g., type "Mg+2" NOT "Mg2+"). The latter formula is ambiguous - it could mean $`Mg_2^+`$ or $`Mg^{+2}`$
:::

(manual-testing)=
## Manually testing a formula

If you want to make sure `pyEQL` is understanding your formula correctly, you can manually test it via `pymatgen` as
follows:

```
>>> from pymatgen.core.ion import Ion
>>> Ion.from_formula(<your_formula>).reduced_formula
...
```

## Formulas you will see when using `Solution`

When using the `Solution` class, 

- When creating a `Solution`, you can enter chemical formulas in any format you prefer, as long as `pymatgen` can understand it (see [manual testing](#manually-testing-a-formula)).
- The keys (solute formulas) in `Solution.components` are preserved in the same format the user enters them. So if you entered `Na+` for sodium ion, it will stay that way.
- Arguments to `Solution.get_property` can be entered in any format you prefer. When `pyEQL` queries the database, it will automatically convert the formula to the canonical one from `pymatgen`
- Property data in the database is uniquely identified by the canonical ion formula (output of `Ion.from_formula(<formula>).reduced_formla`, e.g. "Na[+1]" for sodium ion).
