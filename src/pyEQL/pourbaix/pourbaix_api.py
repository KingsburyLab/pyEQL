import itertools
import warnings
from importlib.resources import files
from typing import Literal

import pandas as pd
from emmet.core.settings import EmmetSettings
from monty.serialization import loadfn
from mp_api.client.core.settings import MAPIClientSettings
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.analysis.pourbaix_diagram import Ion, IonEntry
from pymatgen.core import Composition, Element

from pyEQL import Solution
from pyEQL.utils import standardize_formula

_EMMET_SETTINGS = EmmetSettings()
_MAPI_SETTINGS = MAPIClientSettings()


class Pourbaix_api:
    def __init__(self, mpr, comp_dict: dict | None = None):
        ref_db_file = files("pyEQL") / "pourbaix" / "mpr_reference_ion_database.json"
        ref_xlsx_file = files("pyEQL") / "pourbaix" / "NBS_Tables_Library.xlsx"
        self.json_path = str(ref_db_file)
        self.xlsx_path = str(ref_xlsx_file)
        self.mpr = mpr
        self.comp_dict = comp_dict

    # @classmethod
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

    # @classmethod
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
            ion = Ion.from_formula(i_d["formula"])
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

    # @classmethod
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
        from pymatgen.analysis.pourbaix_diagram import PourbaixEntry  # noqa: PLC0415
        from pymatgen.entries.compatibility import (  # noqa: PLC0415
            Compatibility,
            MaterialsProject2020Compatibility,
            MaterialsProjectAqueousCompatibility,
            MaterialsProjectCompatibility,
        )
        from pymatgen.entries.computed_entries import ComputedEntry  # noqa: PLC0415

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

    def generate_solution_objects(self):
        """
        Args:
            Parsing comp_dict to generate pyEQL solution objects
        Returns:
            List of pyEQL Solution components
        """
        # TODO: Implement the Solution class here to process the comp_dict and do equilibrium calculations
        ion_dict = self.comp_dict
        default_units = "mol/L"
        converted_ion_dict = {
            standardize_formula(f"{ion}[{Element(ion).common_oxidation_states[0]:+d}]"): f"{val} {default_units}"
            for ion, val in ion_dict.items()
        }

        pH_values = [3, 7, 11]  # pH sampling or do we need only one pH?

        speciated_ions = []
        for pH in pH_values:
            sol = Solution(converted_ion_dict, pH=pH, balance_charge="auto")
            sol.equilibrate()
            ion_names = list(sol.components.keys())
            speciated_ions.append(ion_names)  # get_amount that is higher than a threshold

        speciated_ions = list(set(itertools.chain.from_iterable(speciated_ions)))

        speciated_ions = [
            ion for ion in speciated_ions if ion not in ["H[+1]", "OH[-1]", "H2(aq)", "H2O(aq)", "O2(aq)"]
        ]

        return speciated_ions  # noqa: RET504

    @staticmethod
    def standardize_NBS_formula(x):
        try:
            return standardize_formula(x)
        except Exception as e:
            return f"ERROR: {e}"

    def NBS_table_ion_data(self):
        """
        Docstring for NBS_table_data
        Distinguish between a0 and ai!

        :param self: Description
        """
        nbs_data = self.xlsx_path
        df_dict = pd.read_excel(nbs_data, sheet_name="NBS Tables", skiprows=3)
        df_dict = df_dict[df_dict.iloc[:, 4].isin(["ao", "ai"])]
        df_dict.iloc[:, 0] = df_dict.iloc[:, 0].apply(self.standardize_NBS_formula)  # Unnamed: 8 is delta G_f

        # Dictionary that maps the column names with formation energy
        nbs_db = {}
        for i in range(len(df_dict)):
            identifier = df_dict.iloc[i, 0]
            nbs_db[identifier] = {
                "exp_form_E": {"value": df_dict.iloc[i, 8], "units": "kJ/mol"},
                "exp_entropy": {"value": df_dict.iloc[i, 9], "units": "J/(mol*K)"},
            }

        return nbs_db

    def modified_get_ion_reference_data_for_chemsys(self, chemsys: str | list, nbs_db: dict | None = None):
        """
        Docstring for modified_get_ion_reference_data_for_chemsys

        :param self: Description
        :param chemsys: Description
        :type chemsys: str | list
        """
        ion_data = loadfn(self.json_path)
        ion_in_sol = self.generate_solution_objects()

        if nbs_db is None:
            nbs_db = self.NBS_table_ion_data()

        if isinstance(chemsys, str):
            chemsys = chemsys.split("-")

        # for [i, d] in enumerate(ion_data):
        #     if d['data']['MajElements'] == 'K':
        #         del ion_data[i]

        for entry in ion_in_sol:
            identifier = entry
            if identifier in ion_data:
                continue
            if identifier in nbs_db:
                comp_name = identifier.split("[")[0].split("(")[0]
                comp_name = Composition(comp_name)
                maj_elements = [i.symbol for i in comp_name.elements if i.symbol not in ["H", "O"]]

                matched_ref_solids = [entry for entry in ion_data if entry["data"]["MajElements"] == maj_elements[0]]

                if matched_ref_solids:
                    ref_solid = matched_ref_solids[0]["data"]["RefSolid"]
                    ref_form_E_str = matched_ref_solids[0]["data"]["\u0394G\u1da0RefSolid"]["display"]
                    ref_form_E = float(matched_ref_solids[0]["data"]["\u0394G\u1da0RefSolid"]["value"])

                if "[" in identifier and "]" in identifier:
                    charge_str = identifier[identifier.find("[") + 1 : identifier.find("]")]
                    charge = int(charge_str)
                else:
                    charge_str = "0"
                    charge = 0

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
                        "MajElements": maj_elements[0],
                        "RefSolid": ref_solid,
                        "\u0394G\u1da0RefSolid": {"display": ref_form_E_str, "value": ref_form_E, "unit": "kJ/mol"},
                        "reference": "D. D. Wagman et al., Selected values of chemical thermodynamic properties, NBS Technical note 270, Washington; 1968-1971",
                    },
                }
                ion_data.append(ion_record)
            else:
                print(f"Warning: {identifier} not found in NBS database.")
        # pprint.pprint([d for d in ion_data if d["data"]["MajElements"] in chemsys])
        return [d for d in ion_data if d["data"]["MajElements"] in chemsys]
