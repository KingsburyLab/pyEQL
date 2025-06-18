"""Reference data for benchmarking."""

from pyEQL.benchmark import BenchmarkEntry


# TODO: search for each property in each reference database & sketch sub-loader block
# TODO: identify/understand theoretical, reference, and package equivalences (i.e., how do literature values translate
# TODO:     to Python attributes)
# TODO: verify that _get_solute/ion_property framework will work for reference properties
# TODO: consolidate existing data files
# TODO: remove redundant sources
# TODO: write sub-loaders for different CRC archive types
# TODO: write loading function for other reference databases
def _load_crc_data() -> list[BenchmarkEntry]:
    # Use JSONStore?
    # datasets = []

    # SOLUTE DATA

    # property: molar electrical conductivity
    # Salts: KCl, NaCl, MgCl2, NaF, Na2SO4, OH-, FeCl3, H2O
    # Concentrations: 0.001, 0.002
    # sources: molar_electrical_conductivity/

    # SOLUTION DATA
    # property: mean activity coefficient data
    # Salts: HCl, CsI, BaCl2, LiCl, RbCl, MgCl2, KBr, K2SO4
    # Concentrations: 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5
    # TODO: compile data
    # sources: activity_coefficients/

    # load density data

    # load electrical conductivity

    # ? load osmotic coefficient data

    # merge solution.solution_data

    # previous parametrization of ion pairs
    # cations = [("H+", 1), ("Cs+", 1), ("Li+", 1), ("Rb+", 1), ("K+", 1), ("Na", 1), ("Mg", 2), ("Ba", 2)]
    # anions = [("Cl-", 1), ("I-", 1), ("Br", 1), ("SO4-2", 2)]
    # previous list of concentrations to test, mol/kg
    # conc_list = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5]

    return []


def _load_idst_data() -> list[BenchmarkEntry]:
    pass


def _load_jpcrd_data() -> list[BenchmarkEntry]:
    # SOLUTION DATA
    # property: mean activity coefficient
    # Salts: NaCl
    # Concentrations: 0.1, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6
    # TODO: compile data
    # sources:
    #   - Standard Reference Data Publications 1964-1984, NBS, Sauerwein J. C., Dalton, G. R. "Properties Index" (p. 94)
    #   - Theoretical Mean Activity Coefficients of Strong Electrolytes in Aqueous Solutions from 0 to 100C -
    #     Walter J. Hamer. NSRDS-NBS 24, 271p. (1968).

    # Additional properties from Standard Reference Data Publications 1964-1984, NBS, Sauerwein J. C., Dalton, G. R..
    # "Properties Index" (p. 94):
    # density, diffusion coefficient, electrical resistivity/conductivity, osmotic coefficients, dielectric constants
    # Debye length
    pass


def _load_environ_sci_data() -> list[BenchmarkEntry]:
    # SOLUTION DATA
    # property: debye length
    # Salts: NaCl, Na2SO4
    # Concentrations: 0.1, 10
    # source: doi:10.1021/es400571g
    pass


def _load_fluid_phase_data() -> list[BenchmarkEntry]:
    # SOLUTION DATA
    # property: dielectric constant
    # Salts: NaCl, KBr, LiCl, RbCl
    # Concentrations: 0.5, 1, 2, 2.1, 3.4, 4.4, 5, 6.5, 12
    # source: doi:10.1016/j.fluid.2014.05.037
    pass


def _load_jchem_eng_data1() -> list[BenchmarkEntry]:
    # SOLUTION DATA
    # property: activity coefficient
    # Salts: Na+/K+/NO3-
    # Concentrations: 1:3, 1:1, 3:1
    # source: doi:10.1021/je9004432
    pass


def _load_jchem_eng_data2() -> list[BenchmarkEntry]:
    # SOLUTION DATA
    # property: osmotic coefficient
    # Salts: NH4NO3, Cu(2)SO4
    # Concentrations: 0.25, 0.5, 0.75, 1, 1.5, 2
    # source: doi:10.1021/je2009329
    pass
