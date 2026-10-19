(chemistry)=

# Writing Formulas

`pyEQL` interprets the chemical formula of a substance to calculate its molecular
weight and formal charge. The formula is also used as a key to search the
[property database](database.md) for parameters (e.g. diffusion coefficient) that are
used in subsequent calculations.

## How to Enter Valid Chemical Formulas

Generally speaking, type the chemical formula of your solute the "normal" way
and `pyEQL` should be able to interpret it. Internally, `pyEQL` uses a utility function `pyEQL.utils.standardize_formula`
to process all formulas into a standard form. The formula is first parsed into element amounts and a charge by the
[`pymatgen.core.ion.Ion`](https://pymatgen.org/pymatgen.core.html#pymatgen.core.ion.Ion) class, so anything that the
`Ion` class can understand will be processed into a valid formula by `pyEQL`.

Here are some examples:

| Substance       |             You enter             | `pyEQL` understands |
| :-------------- | :-------------------------------: | :-----------------: |
| Sodium Chloride |   "NaCl", "NaCl(aq)", or "ClNa"   |     "NaCl(aq)"      |
| Sodium Sulfate  |      "Na2(SO4)" or "Na2SO4"       |    "Na2SO4(aq)"     |
| Sodium Ion      | "Na+", "Na+1", "Na1+", or "Na[+]" |      "Na[+1]"       |
| Magnesium Ion   |    "Mg+2", "Mg++", or "Mg[++]"    |      "Mg[+2]"       |
| Methanol        |          "CH3OH", "CH4O"          |     "CH3OH(aq)"     |
| Phosphate Ion   |         "PO4-3", "PO₄³⁻"          |      "PO4[-3]"      |
| Mercurous Ion   |        "Hg2+2" or "Hg2[2+]"       |      "Hg2[+2]"      |
| Ethane          |         "C2H6" or "CH3CH3"        |     "C2H6(aq)"      |

The standardized formula (shown in the right hand column of the table) is used to identify solutes. It follows these
rules:

- For charged species, the charge is always placed inside square brackets (e.g., `Na[+1]`) and always includes the
  charge number (even for monovalent ions). Uncharged species are always suffixed by `(aq)` to disambiguate them from
  solids.
- The formula is **never reduced** to its simplest ratio of elements, so multimers keep their identity. For example,
  the mercurous ion is `Hg2[+2]` (not `Hg[+1]`) and persulfate is `S2O8[-2]` (not `SO4[-1]`).
- Organic species are written in [Hill order](https://en.wikipedia.org/wiki/Chemical_formula#Hill_system) (carbon,
  then hydrogen, then other elements alphabetically), e.g. `C2H6(aq)` or `C2H3NaO2(aq)` for sodium acetate.
  Carbonates, bicarbonates, cyanides, cyanates, thiocyanates, CO, and CO2 are treated as inorganic.
- Inorganic species are written with elements in order of increasing electronegativity, e.g. `NaCl(aq)` or
  `Al(OH)4[-1]`.
- A few common species use a preferred conventional name instead, e.g. `NH4[+1]`, `H2PO4[-1]`, `SCN[-1]`,
  `CH3COO[-1]`, and `CH3OH(aq)`.

```{important}
**When writing multivalent ion formulas, it is strongly recommended that you put the charge number AFTER the + or -
sign** (e.g., type "Mg+2" NOT "Mg2+"). The latter formula is ambiguous - it could mean $Mg_2^+$ or $Mg^{+2}$.

`pyEQL` interprets ambiguous formulas like this as follows, and logs a warning when it does:

- Monatomic metal ions are read with the final number as the charge if that is a known oxidation state, e.g. "Mg2+"
  becomes `Mg[+2]` and "Fe3+" becomes `Fe[+3]`. "S2-", "Se2-", and "Te2-" become `S[-2]`, `Se[-2]`, and `Te[-2]`.
- Polyatomic ions such as "SO42-" or "Cr2O72-" are read with the final digit as the charge (giving `SO4[-2]` and
  `Cr2O7[-2]`) when the other reading would give a subscript of 10 or more.
- Otherwise, the number is read as a subscript, e.g. "NH4+" is `NH4[+1]` and "I3-" is triiodide, `I3[-1]`. For
  metal oxo-cations such as "VO2+", which could mean `VO2[+1]` or vanadyl (`VO[+2]`), the subscript reading is kept
  and a warning is logged.

You can also use unicode superscript characters (e.g., "Co²⁺"), and `pyEQL` will understand
these regardless of the order of the `+` and the `2`. So you can write "Co²⁺" and it will be
correctly standardized to `Co[+2]`
```

(manual-testing)=

## Manually testing a formula

If you want to make sure `pyEQL` is understanding your formula correctly, you can manually test it as follows:

```
>>> from pyEQL.utils import standardize_formula
>>> standardize_formula(<your_formula>)
...
```

## Formulas you will see when using `Solution`

When using the `Solution` class,

- When creating a `Solution`, you can enter chemical formulas in any format you prefer, as long as `standardize_formula` can understand it (see [manual testing](#manually-testing-a-formula)).
- The keys (solute formulas) in `Solution.components` are standardized. So if you entered `Na+` for sodium ion, it will appear in `components` as `Na[+1]`.
- However, the `components` attribute is a special dictionary that automatically standardizes formulas when accessed. So, you can still enter the formula
  however you want. For example, the following all access or modify the same element in `components`:

  ```python
  >>> Solution.components.get('Na+')
  >>> Solution.components["Na+1"]
  >>> Solution.components.update("Na[+]": 2)
  >>> Solution.components["Na[+1]"]
  ```

- Arguments to `Solution.get_property` can be entered in any format you prefer. When `pyEQL` queries the database, it will automatically standardize the formula.
- Property data in the database is uniquely identified by the standardized ion formula (output of `standardize_formula(<formula>)`, e.g. "Na[+1]" for sodium ion).
