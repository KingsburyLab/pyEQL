import itertools
import re
import warnings
from importlib.resources import files
from pathlib import Path
from typing import Literal

from emmet.core.settings import EmmetSettings
from monty.serialization import loadfn
from mp_api.client.core.settings import MAPIClientSettings
from openpyxl import load_workbook
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.core import Composition, Element
from pymatgen.core.ion import Ion

from pyEQL import Solution
from pyEQL.engines import Phreeqc2026EOS
from pyEQL.pourbaix.pourbaix_diagram import IonEntry
from pyEQL.utils import standardize_formula

_EMMET_SETTINGS = EmmetSettings()
_MAPI_SETTINGS = MAPIClientSettings()

DEFAULT_REFERENCE_SOLIDS = {
    "Li": {"ref_solid": "Li2CO3", "G_ref_solid": -1044.44},
    "Na": {"ref_solid": "Na2CO3", "G_ref_solid": -1044.4},
    "K": {"ref_solid": "KCl", "G_ref_solid": -409.14},
    "Mg": {"ref_solid": "MgCO3", "G_ref_solid": -1012.1},
    "Ca": {"ref_solid": "CaO", "G_ref_solid": -604.03},
    "Cl": {"ref_solid": "KCl", "G_ref_solid": -409.14},
    "S": {"ref_solid": "CaS", "G_ref_solid": -477.4},
    "N": {"ref_solid": "Ca(NO3)2", "G_ref_solid": -743.07},
    "C": {"ref_solid": "Na2CO3", "G_ref_solid": -1044.4},
    "Fe": {"ref_solid": "Fe3O4", "G_ref_solid": -1015.4},
    "Al": {"ref_solid": "AlHO2", "G_ref_solid": -915.85},
    "P": {"ref_solid": "PH3O4", "G_ref_solid": -1119.1},
}


class Pourbaix_api:
    def __init__(
        self,
        mpr,
        ref_solids: dict | None = None,
        ref_db_file: str | Path | None = None,
        ref_xlsx_file: str | Path | None = None,
    ):
        """
        Construct Pourbaix entries from Materials Project database.

        Args:
            mpr: Materials Project API client used to retreieve DFT entries.
            ref_solids: Reference solids used for DFT ion-reference construction. Defaults to DEFAULT_REFERENCE_SOLIDS.
            ref_db_file: Path to the Materials Project ion-reference database. Defaults to the mpr_reference_ion_database.json packaged within pyEQL.
            ref_xlsx_file: Path to the NBS thermodynamic tables. Defaults to the NBS_Tables_Library.xlsx table packaged within pyEQL.
        """
        pbx_dir = files("pyEQL") / "pourbaix"
        self.json_path = str(ref_db_file or pbx_dir / "mpr_reference_ion_database.json")
        self.xlsx_path = str(ref_xlsx_file or pbx_dir / "NBS_Tables_Library.xlsx")
        self.mpr = mpr
        self.ref_solids = {**DEFAULT_REFERENCE_SOLIDS, **(ref_solids or {})}

    @classmethod
    def get_ion_reference_data_for_chemsys(self, chemsys: str | list) -> list[dict]:
        """Download aqueous ion reference data used in the construction of Pourbaix diagrams.

        Use this method to examine the ion reference data and to add additional
        ions if desired. The data returned from this method can be passed to
        get_ion_entries().

        Data are retrieved from the Aqueous Ion Reference Data project
        hosted on MPContribs. Refer to that project and its associated documentation
        for more details about the format and meaning of the data.

        Args:
            chemsys (str or [str]): Chemical system string comprising element
                symbols separated by dashes, e.g., "Li-Fe-O" or List of element
                symbols, e.g., ["Li", "Fe", "O"].

        Returns:
            [dict]: Among other data, each record contains 1) the experimental ion  free energy, 2) the
                formula of the reference solid for the ion, and 3) the experimental free energy of the
                reference solid. All energies are given in kJ/mol. An example is given below.

                {'identifier': 'Li[+]',
                'formula': 'Li[+]',
                'data': {'charge': {'display': '1.0', 'value': 1.0, 'unit': ''},
                'ΔGᶠ': {'display': '-293.71 kJ/mol', 'value': -293.71, 'unit': 'kJ/mol'},
                'MajElements': 'Li',
                'RefSolid': 'Li2O',
                'ΔGᶠRefSolid': {'display': '-561.2 kJ/mol',
                    'value': -561.2,
                    'unit': 'kJ/mol'},
                'reference': 'H. E. Barner and R. V. Scheuerman, Handbook of thermochemical data for
                compounds and aqueous species, Wiley, New York (1978)'}}
        """
        ion_data = loadfn(self.json_path)

        if isinstance(chemsys, str):
            chemsys = chemsys.split("-")
        return [d for d in ion_data if d["data"]["MajElements"] in chemsys]

    @classmethod
    def get_ion_entries(self, pd: PhaseDiagram, ion_ref_data: list[dict] | None = None) -> list[IonEntry]:
        """Retrieve IonEntry objects that can be used in the construction of
        Pourbaix Diagrams. The energies of the IonEntry are calculaterd from
        the solid energies in the provided Phase Diagram to be
        consistent with experimental free energies.

        NOTE! This is an advanced method that assumes detailed understanding
        of how to construct computational Pourbaix Diagrams. If you just want
        to build a Pourbaix Diagram using default settings, use get_pourbaix_entries.

        Args:
            pd: Solid phase diagram on which to construct IonEntry. Note that this
                Phase Diagram MUST include O and H in its chemical system. For example,
                to retrieve IonEntry for Ti, the phase diagram passed here should contain
                materials in the H-O-Ti chemical system. It is also assumed that solid
                energies have already been corrected with MaterialsProjectAqueousCompatibility,
                which is necessary for proper construction of Pourbaix diagrams.
            ion_ref_data: Aqueous ion reference data. If None (default), the data
                are downloaded from the Aqueous Ion Reference Data project hosted
                on MPContribs. To add a custom ionic species, first download
                data using get_ion_reference_data, then add or customize it with
                your additional data, and pass the customized list here.

        Returns:
            [IonEntry]: IonEntry are similar to PDEntry objects. Their energies
                are free energies in eV.
        """
        # determine the chemsys from the phase diagram
        chemsys = "-".join([el.symbol for el in pd.elements])

        # raise ValueError if O and H not in chemsys
        if "O" not in chemsys or "H" not in chemsys:
            raise ValueError(
                f"The phase diagram chemical system must contain O and H! Your diagram chemical system is {chemsys}."
            )

        # ion_data = self.get_ion_reference_data_for_chemsys(chemsys) if not ion_ref_data else ion_ref_data
        ion_data = ion_ref_data if ion_ref_data else self.get_ion_reference_data_for_chemsys(chemsys)

        # position the ion energies relative to most stable reference state
        ion_entries = []
        for _, i_d in enumerate(ion_data):
            formula_cleaned = re.sub(r"\(aq\)|\(l\)|\(g\)", "", i_d["formula"])
            ion = Ion.from_formula(formula_cleaned)
            refs = [e for e in pd.all_entries if e.composition.reduced_formula == i_d["data"]["RefSolid"]]
            if not refs:
                raise ValueError("Reference solid not contained in entry list")
            stable_ref = sorted(refs, key=lambda x: x.energy_per_atom)[0]
            rf = stable_ref.composition.get_reduced_composition_and_factor()[1]

            # TODO - need a more robust way to convert units
            # use pint here?
            if i_d["data"]["ΔGᶠRefSolid"]["unit"] == "kJ/mol":
                # convert to eV/formula unit
                ref_solid_energy = i_d["data"]["ΔGᶠRefSolid"]["value"] / 96.485
            elif i_d["data"]["ΔGᶠRefSolid"]["unit"] == "MJ/mol":
                # convert to eV/formula unit
                ref_solid_energy = i_d["data"]["ΔGᶠRefSolid"]["value"] / 96485
            else:
                raise ValueError(f"Ion reference solid energy has incorrect unit {i_d['data']['ΔGᶠRefSolid']['unit']}")
            solid_diff = pd.get_form_energy(stable_ref) - ref_solid_energy * rf
            elt = i_d["data"]["MajElements"]
            correction_factor = ion.composition[elt] / stable_ref.composition[elt]
            # TODO - need a more robust way to convert units
            # use pint here?
            if i_d["data"]["ΔGᶠ"]["unit"] == "kJ/mol":
                # convert to eV/formula unit
                ion_free_energy = i_d["data"]["ΔGᶠ"]["value"] / 96.485
            elif i_d["data"]["ΔGᶠ"]["unit"] == "MJ/mol":
                # convert to eV/formula unit
                ion_free_energy = i_d["data"]["ΔGᶠ"]["value"] / 96485
            else:
                raise ValueError(f"Ion free energy has incorrect unit {i_d['data']['ΔGᶠ']['unit']}")
            energy = ion_free_energy + solid_diff * correction_factor
            ion_entries.append(IonEntry(ion, energy))

        return ion_entries

    @classmethod
    def get_pourbaix_entries(
        self,
        chemsys: str | list,
        solid_compat="MaterialsProject2020Compatibility",
        use_gibbs: Literal[300] | None = None,
    ):
        """A helper function to get all entries necessary to generate
        a Pourbaix diagram from the rest interface.

        Args:
            chemsys (str or [str]): Chemical system string comprising element
                symbols separated by dashes, e.g., "Li-Fe-O" or List of element
                symbols, e.g., ["Li", "Fe", "O"].
            solid_compat: Compatibility scheme used to pre-process solid DFT energies prior
                to applying aqueous energy adjustments. May be passed as a class (e.g.
                MaterialsProject2020Compatibility) or an instance
                (e.g., MaterialsProject2020Compatibility()). If None, solid DFT energies
                are used as-is. Default: MaterialsProject2020Compatibility
            use_gibbs: Set to 300 (for 300 Kelvin) to use a machine learning model to
                estimate solid free energy from DFT energy (see GibbsComputedStructureEntry).
                This can slightly improve the accuracy of the Pourbaix diagram in some
                cases. Default: None. Note that temperatures other than 300K are not
                permitted here, because MaterialsProjectAqueousCompatibility corrections,
                used in Pourbaix diagram construction, are calculated based on 300 K data.
        """
        # imports are not top-level due to expense
        from pymatgen.entries.compatibility import (  # noqa: PLC0415
            Compatibility,
            MaterialsProject2020Compatibility,
            MaterialsProjectCompatibility,
        )
        from pymatgen.entries.computed_entries import ComputedEntry  # noqa: PLC0415

        from pyEQL.pourbaix.compatibility import MaterialsProjectAqueousCompatibility  # noqa: PLC0415
        from pyEQL.pourbaix.pourbaix_diagram import PourbaixEntry  # noqa: PLC0415

        if solid_compat == "MaterialsProjectCompatibility":
            solid_compat = MaterialsProjectCompatibility()
        elif solid_compat == "MaterialsProject2020Compatibility":
            solid_compat = MaterialsProject2020Compatibility()
        elif isinstance(solid_compat, Compatibility):
            pass
        else:
            raise ValueError(
                "Solid compatibility can only be 'MaterialsProjectCompatibility', "
                "'MaterialsProject2020Compatibility', or an instance of a Compatibility class"
            )

        pbx_entries = []

        if isinstance(chemsys, str):
            chemsys = chemsys.split("-")
        # capitalize and sort the elements
        chemsys = sorted(e.capitalize() for e in chemsys)

        # Get ion entries first, because certain ions have reference
        # solids that aren't necessarily in the chemsys (Na2SO4)

        # download the ion reference data from MPContribs
        ion_data = self.get_ion_reference_data_for_chemsys(chemsys)

        # build the PhaseDiagram for get_ion_entries
        ion_ref_comps = [Ion.from_formula(d["data"]["RefSolid"]).composition for d in ion_data]
        ion_ref_elts = set(itertools.chain.from_iterable(i.elements for i in ion_ref_comps))
        # TODO - would be great if the commented line below would work
        # However for some reason you cannot process GibbsComputedStructureEntry with
        # MaterialsProjectAqueousCompatibility

        # if mpr is None:
        #     # raise ValueError("MPRester object is required")
        #     print("MPRester object is not provided, using default API key")
        #     api_key = os.getenv("MP_API_KEY", "12345678901234567890123456789012")
        #     mpr = MPRester(api_key=api_key)

        ion_ref_entries = self.mpr.get_entries_in_chemsys(list([str(e) for e in ion_ref_elts] + ["O", "H"]))

        # suppress the warning about supplying the required energies; they will be calculated from the
        # entries we get from MPRester
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="You did not provide the required O2 and H2O energies.",
            )
            compat = MaterialsProjectAqueousCompatibility(solid_compat=solid_compat)
        # suppress the warning about missing oxidation states
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Failed to guess oxidation states.*")
            ion_ref_entries = compat.process_entries(ion_ref_entries)  # type: ignore
        # TODO - if the commented line above would work, this conditional block
        # could be removed
        if use_gibbs:
            # replace the entries with GibbsComputedStructureEntry
            from pymatgen.entries.computed_entries import GibbsComputedStructureEntry  # noqa: PLC0415

            ion_ref_entries = GibbsComputedStructureEntry.from_entries(ion_ref_entries, temp=use_gibbs)

        ion_ref_pd = PhaseDiagram(ion_ref_entries)  # type: ignore

        ion_entries = self.get_ion_entries(ion_ref_pd, ion_ref_data=ion_data)
        pbx_entries = [PourbaixEntry(e, f"ion-{n}") for n, e in enumerate(ion_entries)]

        # Construct the solid pourbaix entries from filtered ion_ref entries
        extra_elts = set(ion_ref_elts) - {Element(s) for s in chemsys} - {Element("H"), Element("O")}
        for entry in ion_ref_entries:
            entry_elts = set(entry.composition.elements)
            # Ensure no OH chemsys or extraneous elements from ion references
            if not (entry_elts <= {Element("H"), Element("O")} or extra_elts.intersection(entry_elts)):
                # Create new computed entry
                form_e = ion_ref_pd.get_form_energy(entry)  # type: ignore
                new_entry = ComputedEntry(entry.composition, form_e, entry_id=entry.entry_id)
                pbx_entry = PourbaixEntry(new_entry)
                pbx_entries.append(pbx_entry)

        return pbx_entries

    def ion_pourbaix_entries(self, pbx_entries):
        """
        Newly added PHREEQC speciated ions are defined with their PHREEQC concentrations.
        """

        # imports are not top-level due to expense
        def _normalize_charge(identifier):
            return identifier.replace("[-]", "[-1]").replace("[+]", "[+1]")

        added_ion_conc_map = getattr(self, "added_ion_conc_map", {})

        for entry in pbx_entries:
            phase_type = getattr(entry, "phase_type", None)

            is_ion = (
                ("Ion" in phase_type)
                if isinstance(
                    phase_type,
                    list | tuple | set,
                )
                else phase_type == "Ion"
            )

            if not is_ion:
                continue

            entry_name = _normalize_charge(entry.name)

            if entry_name in added_ion_conc_map:
                old_conc = entry.concentration
                entry.concentration = float(added_ion_conc_map[entry_name])

                print(f"Updated PHREEQC ion concentration: {entry.name}: {old_conc:g} -> {entry.concentration:g} M")

        return pbx_entries

    def generate_solution_objects(self, comp_dict: dict | None = None):
        """
        Args:
            Parsing comp_dict to generate pyEQL solution objects
        Returns:
            List of pyEQL Solution components
        """
        # TODO: Implement the Solution class here to process the comp_dict and do equilibrium calculations
        ion_dict = comp_dict
        default_units = "mol/L"
        custom_eos = Phreeqc2026EOS(phreeqc_db="phreeqc.dat")

        converted_ion_dict = {standardize_formula(ion): f"{val} {default_units}" for ion, val in ion_dict.items()}

        pH_values = [3, 7, 11]  # pH sampling or do we need only one pH?

        excluded_species = {"H[+1]", "OH[-1]", "H2(aq)", "H2O(aq)", "O2(aq)"}

        speciated_ions = {}
        for pH in pH_values:
            sol = Solution(converted_ion_dict, pH=pH, balance_charge="auto", engine=custom_eos)
            try:
                sol.equilibrate()
                print(f"Equilibration succeeded at pH {pH}, here are the proportions:")
            except Exception as e:
                print(f"Equilibration failed at pH {pH} with error: {e}")
                continue
            tds = sol.total_dissolved_solids.magnitude

            for key in sol.components:
                conc_val = sol.get_amount(key, "mg/L").magnitude
                print(f"{key}: {conc_val} / {tds}: {conc_val / tds:.2%}")

                if conc_val / tds < 0.025 or "unk" in key or key in excluded_species:
                    continue

                # try:
                # conc_val_mol_L = sol.get_amount(key, "mol/L").magnitude
                if "[" in key:
                    conc_activity = sol.get_activity(key).magnitude
                else:
                    aq_key = key.removesuffix("(aq)").strip()
                    conc_activity = sol.engine._get_activity(aq_key)
                #     print(f"NaCl: {key}, activity: {conc_activity}")
                # except:
                #     print(f"Key {key}")
                #     print(f"Engine {sol.engine}")

                if key not in speciated_ions:
                    # speciated_ions[key] = conc_val_mol_L
                    speciated_ions[key] = conc_activity
                else:
                    # speciated_ions[key] = max(speciated_ions[key], conc_val_mol_L
                    # )
                    speciated_ions[key] = max(speciated_ions[key], conc_activity)

        speciated_ions = {
            ion: conc
            for ion, conc in speciated_ions.items()
            if ion not in ["H[+1]", "OH[-1]", "H2(aq)", "H2O(aq)", "O2(aq)"]
        }

        print(f"PHREEQC speciated ions: {speciated_ions}")

        return speciated_ions

    @staticmethod
    def _rich_text_formula(value):
        """
        Convert a rich-text chemical formula from Word or Excel documents into pyEQL standardize_formula notation.

        Examples:
        SO4²⁻ -> SO4[-2]
        Mg²⁺ -> Mg[+2]
        """
        # TODO: Move the richtext chemical formula parsing to pyEQL.utils that is attached to the Ion class or standardize_formula function.

        if isinstance(value, str):
            return value.strip()

        formula = ""
        charge = ""

        for run in value:
            text = getattr(run, "text", str(run))
            font = getattr(run, "font", None)

            if getattr(font, "vertAlign", None) == "superscript":
                charge += text
            else:
                formula += text

        formula = formula.strip()
        charge = charge.strip().replace("\N{MINUS SIGN}", "-")

        if not charge:
            return formula

        if formula.count("(") > formula.count(")"):
            formula += ")"

        return f"{formula}[{charge[-1]}{charge[:-1] or '1'}]"

    def NBS_table_ion_data(self):
        """
        Load aqueous ion and ion complexes free energy of formation data from the NBS Table. The data are used to construct the aqueous ion reference database.

        Returns:
            dict: Thermodynamic dictionary with aqueous ion and ion complexes free energy of formation data and entropy data.
        """
        nbs_data = self.xlsx_path
        workbook = load_workbook(nbs_data, data_only=True, rich_text=True)
        worksheet = workbook["NBS Tables"]

        nbs_db = {}
        for row in worksheet.iter_rows(min_row=5, values_only=False):
            if row[4].value not in {"ao", "ai"}:
                continue

            formula = self._rich_text_formula(row[0].value)
            formula = formula.replace("·", "").replace("∙", "")

            if "[" in formula and "]" in formula:
                charge_str = formula[formula.find("[") + 1 : formula.find("]")]
                if charge_str not in {"+", "-"}:
                    charge = float(charge_str)
                    if abs(charge - round(charge)) > 1e-8:
                        continue

            identifier = standardize_formula(formula)

            nbs_db[identifier] = {
                "exp_form_E": {"value": row[8].value, "units": "kJ/mol"},
                "exp_entropy": {"value": row[9].value, "units": "J/(mol*K)"},
            }

        return nbs_db

    def modified_get_ion_reference_data_for_chemsys(self, chemsys: str | list, nbs_db: dict | None = None):
        """
        Modified the Pymatgen's get_ion_reference_data_for_chemsys method to include additional ions from the PHREEQC database, which are not present in mpr_reference_ion_database.json.

        Args:
            chemsys (str or [str]): Chemical system string comprising element
                symbols separated by dashes, e.g., "Li-Fe-O" or List of element
                symbols, e.g., ["Li", "Fe", "O"].

        Returns:
            [dict]: Among other data, each record contains 1) the experimental ion  free energy, 2) the
                formula of the reference solid for the ion, and 3) the experimental free energy of the
                reference solid. All energies are given in kJ/mol. An example of ion complex is given below.

                {'identifier': 'CaSO4(aq)',
                'formula': 'CaSO4(aq)',
                'data': {'charge': {'display': '0.0', 'value': 0.0, 'unit': ''},
                'ΔGᶠ': {'display': '-1298.10 kJ/mol', 'value': -1298.10, 'unit': 'kJ/mol'},
                'MajElements': 'Ca',
                'RefSolid': 'CaO',
                'ΔGᶠRefSolid': {'display': '-604.03 kJ/mol',
                    'value': -604.03,
                    'unit': 'kJ/mol'},
                'reference': 'D. D. Wagman et al., Selected values of chemical thermodynamic properties, NBS Technical note 270, Washington; 1968-1971'}}
        """

        ion_data = loadfn(self.json_path)

        def _normalize_charge(identifier):
            return identifier.replace("[-]", "[-1]").replace("[+]", "[+1]")

        ion_in_sol_init = self.generate_solution_objects()

        if isinstance(ion_in_sol_init, dict):
            ion_conc_map = {_normalize_charge(identifier): float(conc) for identifier, conc in ion_in_sol_init.items()}
            ion_in_sol = list(ion_conc_map.keys())
        else:
            ion_conc_map = {}
            ion_in_sol = [_normalize_charge(identifier) for identifier in ion_in_sol_init]

        if nbs_db is None:
            nbs_db = self.NBS_table_ion_data()

        if isinstance(chemsys, str):
            chemsys = chemsys.split("-")

        existing_identifiers = {
            _normalize_charge(d["formula"]) for d in ion_data if isinstance(d, dict) and "formula" in d
        }

        self.added_ion_conc_map = {}

        for identifier in ion_in_sol:
            # Skip if already present
            if identifier in existing_identifiers:
                continue

            # Skip if not found in NBS db
            if identifier not in nbs_db:
                print(f"Warning: {identifier} not found in NBS database.")
                continue

            # TODO - instead of manually parsing elements, query the pyEQL db or rely on the Solute class
            comp_name = identifier.split("[")[0].split("(")[0]
            comp_name = Composition(comp_name)
            maj_elements = [i.symbol for i in comp_name.elements if i.symbol not in ["H", "O"]]

            maj_element = maj_elements[0]

            if maj_element not in self.ref_solids:
                print(f"Warning: no reference solid mapping for element {maj_element}")
                continue

            ref_solid = self.ref_solids[maj_element]["ref_solid"]
            G_ref_solid = self.ref_solids[maj_element]["G_ref_solid"]

            # TODO - instead of manually parsing elements, query the pyEQL db or rely on the Solute class
            if "[" in identifier and "]" in identifier:
                charge_str = identifier[identifier.find("[") + 1 : identifier.find("]")]
                charge = float(charge_str)
            else:
                charge_str = "0"
                charge = 0.0

            # to remove cases like Cl[-0.3333]
            import math  # noqa: PLC0415

            if not math.isclose(charge, round(charge), abs_tol=1e-8):
                continue

            ion_record = {
                "identifier": identifier,
                "formula": identifier,
                "data": {
                    "charge": {"display": charge_str, "value": charge, "unit": ""},
                    "\u0394G\u1da0": {
                        "display": f"{nbs_db[identifier]['exp_form_E']['value']} {nbs_db[identifier]['exp_form_E']['units']}",
                        "value": float(nbs_db[identifier]["exp_form_E"]["value"]),
                        "unit": nbs_db[identifier]["exp_form_E"]["units"],
                    },
                    "MajElements": maj_element,
                    "RefSolid": ref_solid,
                    "\u0394G\u1da0RefSolid": {
                        "display": f"{G_ref_solid} kJ/mol",
                        "value": G_ref_solid,
                        "unit": "kJ/mol",
                    },
                    "reference": "D. D. Wagman et al., Selected values of chemical thermodynamic properties, NBS Technical note 270, Washington; 1968-1971",
                },
            }

            if identifier in ion_conc_map:
                self.added_ion_conc_map[identifier] = ion_conc_map[identifier]

            ion_data.append(ion_record)

        for d in ion_data:
            maj_element = d["data"]["MajElements"]

            if maj_element in self.ref_solids:
                ref_solid = self.ref_solids[maj_element]

                d["data"]["RefSolid"] = ref_solid["ref_solid"]
                d["data"]["ΔGᶠRefSolid"] = {
                    "display": f"{ref_solid['G_ref_solid']} kJ/mol",
                    "value": ref_solid["G_ref_solid"],
                    "unit": "kJ/mol",
                }

        return [d for d in ion_data if d["data"]["MajElements"] in chemsys]
