(database)=

# Property Database

pyEQL is distributed with a database of solute properties and model parameters needed to perform
its calculations. The database includes:

- Molecular weight, charge, and other chemical informatics information for any species
- Diffusion coefficients for 104 ions
- Pitzer model activity correction coefficients for 157 salts
- Pitzer model partial molar volume coefficients for 120 salts
- Jones-Dole "B" coefficients for 83 ions
- Hydrated and ionic radii for 23 ions
- Dielectric constant model parameters for 18 ions
- Partial molar volumes for 24 ions

`pyEQL` can automatically infer basic chemical informatics such as molecular weight and charge by passing a solute's formula to `pymatgen.core.ion.Ion` (See [chemical formulas](#chemistry)). For other physicochemical properties, it relies on data compiled into the included database. A list of the data and species covered is available [below](#species-included)

## Format

The database is distributed as a `.json` file containing serialized `Solute` objects that define the schema for aggregated property data (see [below](#the-solute-class)). By default, each instance of `Solution` loads this file as a [`maggma`](https://materialsproject.github.io/maggma/) [`JSONStore`](https://materialsproject.github.io/maggma/reference/stores/#maggma.stores.mongolike.JSONStore) and queries data from it using the [`Store`](https://materialsproject.github.io/maggma/concepts/#store) interface.

If desired, users can point a `Solution` instance to an alternate database by using the `database` keyword argument at creation. The argument should contain either 1) the path to an alternate `.json` file (as a `str`) or 2) a `maggma.Store` instance. The data in the file or `Store` must match the schema defined by `Solute`, with the field `formula` used as the key field (unique identifier).

```
s1 = Solution(database='/path/to/my_database.json')
```

or

```
from maggma.core import JSONStore

db_store = JSONStore('/path/to/my_database.json', key='formula')
s1 = Solution(database=db_store)
```

## The `Solute` class

`pyEQL.Solute` is a [`dataclass`](https://docs.python.org/3/library/dataclasses.html) that defines a schema for organizing solute property data. You can think of the schema as a structured dictionary: `Solute` defines the naming and organization of the keys. You can create a basic `Solute` from just the solute's formula as follows:

```
>>> from pyEQL.solute import Solute
>>> Solute.from_formula('Ti+2')
Solute(formula='Ti[+2]', charge=2, molecular_weight='47.867 g/mol', elements=['Ti'], chemsys='Ti', pmg_ion=Ion: Ti1 +2, formula_html='Ti<sup>+2</sup>', formula_latex='Ti$^{+2}$', formula_hill='Ti', formula_pretty='Ti^+2', oxi_state_guesses=({'Ti': 2.0},), n_atoms=1, n_elements=1, size={'radius_ionic': None, 'radius_hydrated': None, 'radius_vdw': None, 'molar_volume': None}, thermo={'ΔG_hydration': None, 'ΔG_formation': None}, transport={'diffusion_coefficient': None}, model_parameters={'activity_pitzer': {'Beta0': None, 'Beta1': None, 'Beta2': None, 'Cphi': None, 'Max_C': None}, 'molar_volume_pitzer': {'Beta0': None, 'Beta1': None, 'Beta2': None, 'Cphi': None, 'V_o': None, 'Max_C': None}, 'viscosity_jones_dole': {'B': None}})
```

This method uses `pymatgen` to populate the `Solute` with basic chemical information like molecular weight. You can access top-level keys in the schema via attribute, e.g.

```
>>> s.molecular_weight
'47.867 g/mol'
>>> s.charge
2.0
```

Other properties that are present in the schema, but not set, are `None`. For example, here we have not specified a diffusion coefficient. If we inspect the `transport` attribute, we see

```
>>> s.transport
{'diffusion_coefficient': None}
```

You can convert a `Solute` into a regular dictionary using `Solute.as_dict()`

```
>>> s.as_dict()
{'formula': 'Ti[+2]', 'charge': 2, 'molecular_weight': '47.867 g/mol', 'elements': ['Ti'], 'chemsys': 'Ti', 'pmg_ion': Ion: Ti1 +2, 'formula_html': 'Ti<sup>+2</sup>', 'formula_latex': 'Ti$^{+2}$', 'formula_hill': 'Ti', 'formula_pretty': 'Ti^+2', 'oxi_state_guesses': ({'Ti': 2.0},), 'n_atoms': 1, 'n_elements': 1, 'size': {'radius_ionic': None, 'radius_hydrated': None, 'radius_vdw': None, 'molar_volume': None}, 'thermo': {'ΔG_hydration': None, 'ΔG_formation': None}, 'transport': {'diffusion_coefficient': None}, 'model_parameters': {'activity_pitzer': {'Beta0': None, 'Beta1': None, 'Beta2': None, 'Cphi': None, 'Max_C': None}, 'molar_volume_pitzer': {'Beta0': None, 'Beta1': None, 'Beta2': None, 'Cphi': None, 'V_o': None, 'Max_C': None}, 'viscosity_jones_dole': {'B': None}}}
```

## Searching the database

Once you have a created a `Solution`, it will automatically search the database for needed parameters whenever it needs
to perform a calculation. For example, if you call `get_transport_number`, `pyEQL` will search the property database
for diffusion coefficient data to use in the calculation. No user action is needed.

If you want to search the database yourself, or to inspect the values that `pyEQL` uses for a particular parameter, you
can do so via the `get_property` method. First, create a `Solution`

```
>>> from pyEQL import Solution
>>> s1 = pyEQL.Solution
```

Next, call `get_property` with a solute name and the name of the property you need. Valid property names are any key
in the `Solute` schema. Nested keys can be separated by periods, e.g. "model_parameters.activity_pitzer":

```
>>> s1.get_property('Mg+2', 'transport.diffusion_coefficient')
<Quantity(0.00705999997, 'centimeter ** 2 * liter * pascal * second / kilogram / meter ** 2')>
```

If the property exists, it will be returned as a [`pint`](https://pint.readthedocs.io/en/stable/) `Quantity` object, which you can convert to specific units if needed, e.g.

```
>>> s1.get_property('Mg+2', 'transport.diffusion_coefficient').to('m**2/s')
<Quantity(7.05999997e-10, 'meter ** 2 / second')>
```

If the property does not exist in the database, `None` will be returned.

```
>>> s1.get_property('Mg+2', 'transport.randomproperty')
>>>
```

Although the database contains additional context about each and every property value, such as a citation, this information is not currently exposed via the `Solution` interface. Richer methods for exploring and adding to the database may be added in the future.

## Species included

The database currently contains one or more physichochemical properties for each of the solutes listed below. More detailed information about which properties are available for which solutes may be added in the future.

```
 - Ac[+3]
 - Ag(CN)2[-1]
 - AgNO3(aq)
 - Ag[+1]
 - Ag[+2]
 - Ag[+3]
 - Al2(SO4)3(aq)
 - Al[+3]
 - AsO4[-3]
 - Au(CN)2[-1]
 - Au(CN)4[-1]
 - Au[+1]
 - Au[+2]
 - Au[+3]
 - B(H5C6)4[-1]
 - B(OH)3(aq)
 - B(OH)4[-1]
 - BF4[-1]
 - BO2[-1]
 - Ba(ClO4)2(aq)
 - Ba(NO3)2(aq)
 - BaBr2(aq)
 - BaC4O.3H2O(aq)
 - BaCl2(aq)
 - BaI2(aq)
 - Ba[+2]
 - BeSO4(aq)
 - Be[+2]
 - Bi[+3]
 - BrO3[-1]
 - Br[-0.33333333]
 - Br[-1]
 - C2N3[-1]
 - CH3COO[-1]
 - CNO[-1]
 - CN[-1]
 - CO3[-2]
 - CSN[-1]
 - CSeN[-1]
 - Ca(ClO4)2(aq)
 - Ca(NO3)2(aq)
 - CaBr2(aq)
 - CaCl2(aq)
 - CaI2(aq)
 - Ca[+2]
 - Cd(ClO4)2(aq)
 - Cd(NO2)2(aq)
 - Cd(NO3)2(aq)
 - CdSO4(aq)
 - Cd[+2]
 - CeCl3(aq)
 - Ce[+3]
 - Ce[+4]
 - ClO2[-1]
 - ClO3[-1]
 - ClO4[-1]
 - Cl[-1]
 - Co(CN)6[-3]
 - Co(H3N)6[-3]
 - Co(NO3)2(aq)
 - CoBr2(aq)
 - CoCl2(aq)
 - CoI2(aq)
 - Co[+2]
 - Co[+3]
 - Cr(NO3)3(aq)
 - CrCl3(aq)
 - CrO4[-2]
 - Cr[+2]
 - Cr[+3]
 - Cs2SO4(aq)
 - CsBr(aq)
 - CsCl(aq)
 - CsF(aq)
 - CsHC2O.1H2O(aq)
 - CsI(aq)
 - CsNO2(aq)
 - CsNO3(aq)
 - CsOH(aq)
 - Cs[+1]
 - Cu(NO3)2(aq)
 - CuCl2(aq)
 - CuSO4(aq)
 - Cu[+1]
 - Cu[+2]
 - Cu[+3]
 - Dy[+2]
 - Dy[+3]
 - Er[+2]
 - Er[+3]
 - Eu(NO3)3(aq)
 - EuCl3(aq)
 - Eu[+2]
 - Eu[+3]
 - F[-1]
 - Fe(CN)6[-3]
 - Fe(CN)6[-4]
 - FeCl2(aq)
 - FeCl3(aq)
 - Fe[+2]
 - Fe[+3]
 - Ga[+3]
 - GdCl3(aq)
 - Gd[+3]
 - Ge[+2]
 - H2CO3(aq)
 - H2O(aq)
 - H2SNO3[-1]
 - H2SO4(aq)
 - H3O[+1]
 - H4BrN(aq)
 - H4IN(aq)
 - H4N2O3(aq)
 - H4NCl(aq)
 - H4NClO4(aq)
 - H4N[+1]
 - H4SNO4(aq)
 - H5C6O7[-3]
 - H5N2[+1]
 - H8S(NO2)2(aq)
 - HBr(aq)
 - HCO2[-1]
 - HCO3[-1]
 - HCl(aq)
 - HClO4(aq)
 - HF2[-1]
 - HI(aq)
 - HNO3(aq)
 - HO2[-1]
 - HOsO5[-1]
 - HSO3[-1]
 - HSO4[-1]
 - HS[-1]
 - HSeO3[-1]
 - H[+1]
 - Hf[+4]
 - Hg[+2]
 - Ho[+2]
 - Ho[+3]
 - IO3[-1]
 - IO4[-1]
 - I[-1]
 - In[+1]
 - In[+2]
 - In[+3]
 - IrO4[-1]
 - Ir[+3]
 - K2CO3(aq)
 - K2PHO4(aq)
 - K2SO4(aq)
 - K3Fe(CN)6(aq)
 - K3PO4(aq)
 - K4Fe(CN)6(aq)
 - KBr(aq)
 - KBrO3(aq)
 - KCSN(aq)
 - KCl(aq)
 - KClO3(aq)
 - KClO4(aq)
 - KCrO4(aq)
 - KF(aq)
 - KHC2O.1H2O(aq)
 - KHCO3(aq)
 - KI(aq)
 - KNO2(aq)
 - KNO3(aq)
 - KOH(aq)
 - KPO3.1H2O(aq)
 - K[+1]
 - La(NO3)3(aq)
 - LaCl3(aq)
 - La[+3]
 - Li2SO4(aq)
 - LiBr(aq)
 - LiCl(aq)
 - LiClO4(aq)
 - LiHC2O.1H2O(aq)
 - LiI(aq)
 - LiNO2(aq)
 - LiNO3(aq)
 - LiOH(aq)
 - Li[+1]
 - Lu[+3]
 - Mg(ClO4)2(aq)
 - Mg(NO3)2(aq)
 - MgBr2(aq)
 - MgC4O.3H2O(aq)
 - MgCl2(aq)
 - MgI2(aq)
 - MgSO4(aq)
 - Mg[+2]
 - MnCl2(aq)
 - MnO4[-1]
 - MnSO4(aq)
 - Mn[+2]
 - Mn[+3]
 - MoO4[-2]
 - Mo[+3]
 - NO2[-1]
 - NO3[-1]
 - N[-0.33333333]
 - Na2CO3(aq)
 - Na2PHO4(aq)
 - Na2S2O3(aq)
 - Na2SO4(aq)
 - Na3PO4(aq)
 - NaBr(aq)
 - NaBrO3(aq)
 - NaCSN(aq)
 - NaCl(aq)
 - NaClO4(aq)
 - NaCrO4(aq)
 - NaF(aq)
 - NaHC2O.1H2O(aq)
 - NaHC3.2H2O(aq)
 - NaHCO2(aq)
 - NaHCO3(aq)
 - NaI(aq)
 - NaNO2(aq)
 - NaNO3(aq)
 - NaOH(aq)
 - NaPO3.1H2O(aq)
 - Na[+1]
 - Nb[+3]
 - Nd(NO3)3(aq)
 - NdCl3(aq)
 - Nd[+2]
 - Nd[+3]
 - Ni(NO3)2(aq)
 - NiCl2(aq)
 - NiSO4(aq)
 - Ni[+2]
 - Ni[+3]
 - Np[+3]
 - Np[+4]
 - OH[-1]
 - Os[+3]
 - P(HO2)2[-1]
 - P(OH)2[-1]
 - P2O7[-4]
 - P3O10[-5]
 - PF6[-1]
 - PH9(NO2)2(aq)
 - PHO4[-2]
 - PO3F[-2]
 - PO3[-1]
 - PO4[-3]
 - Pa[+3]
 - Pb(ClO4)2(aq)
 - Pb(NO3)2(aq)
 - Pb[+2]
 - Pd[+2]
 - Pm[+2]
 - Pm[+3]
 - Po[+2]
 - PrCl3(aq)
 - Pr[+2]
 - Pr[+3]
 - Pt[+2]
 - Pu[+2]
 - Pu[+4]
 - Ra[+2]
 - Rb2SO4(aq)
 - RbBr(aq)
 - RbCl(aq)
 - RbF(aq)
 - RbHC2O.1H2O(aq)
 - RbI(aq)
 - RbNO2(aq)
 - RbNO3(aq)
 - RbOH(aq)
 - Rb[+1]
 - ReO4[-1]
 - Re[+1]
 - Re[+3]
 - Re[-1]
 - Rh[+3]
 - Ru[+2]
 - Ru[+3]
 - S2O3[-2]
 - SO2[-1]
 - SO3[-1]
 - SO3[-2]
 - SO4[-1]
 - SO4[-2]
 - S[-2]
 - Sb(HO2)2[-1]
 - Sb(OH)6[-1]
 - ScCl3(aq)
 - Sc[+2]
 - Sc[+3]
 - SeO3[-1]
 - SeO4[-1]
 - SeO4[-2]
 - SiF6[-2]
 - SmCl3(aq)
 - Sm[+2]
 - Sm[+3]
 - Sn[+2]
 - Sn[+4]
 - Sr(ClO4)2(aq)
 - Sr(NO3)2(aq)
 - SrBr2(aq)
 - SrCl2(aq)
 - SrI2(aq)
 - Sr[+2]
 - Ta[+3]
 - Tb[+3]
 - TcO4[-1]
 - Tc[+2]
 - Tc[+3]
 - Th(NO3)4(aq)
 - Th[+4]
 - Ti[+2]
 - Ti[+3]
 - Tl(ClO4)3(aq)
 - Tl(NO2)3(aq)
 - Tl(NO3)3(aq)
 - TlH(C3O)2.4H2O(aq)
 - Tl[+1]
 - Tl[+3]
 - Tm[+2]
 - Tm[+3]
 - U(ClO)2(aq)
 - U(ClO5)2(aq)
 - U(NO4)2(aq)
 - UO2[+1]
 - UO2[+2]
 - USO6(aq)
 - U[+3]
 - U[+4]
 - VO2[+1]
 - V[+2]
 - V[+3]
 - WO4[-1]
 - WO4[-2]
 - W[+3]
 - YCl3(aq)
 - YNO3(aq)
 - Y[+3]
 - Yb[+2]
 - Yb[+3]
 - Zn(ClO4)2(aq)
 - Zn(NO3)2(aq)
 - ZnBr2(aq)
 - ZnCl2(aq)
 - ZnI2(aq)
 - ZnSO4(aq)
 - Zn[+2]
 - Zr[+4]
```
