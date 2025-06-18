"""Reference data for benchmarking."""

from pyEQL.benchmark import BenchmarkEntry


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
    # TODO: compile data
    # sources: molar_electrical_conductivity/

    # SOLUTION DATA
    # property: mean activity coefficient
    # Salts: HCl, CsI, BaCl2, LiCl, RbCl, MgCl2, KBr, K2SO4
    # Concentrations: 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5
    # TODO: compile data
    # sources: activity_coefficients/

    # property: density
    # Salts: ?
    # Concentrations: ?
    # TODO: compile data
    # sources: density/

    # property: conductivity
    # Salts: NaCl, KCl, MgCl2
    # Concentrations: 0.001, 0.05, 0.1
    # TODO: compile data
    # sources: "Electrical Conductivity of Aqueous Solutions at 20C as a Function of Concentration"

    return []


def _load_idst_data() -> list[BenchmarkEntry]:
    # SOLUTION DATA
    # property: mean activity coefficient
    # Salts: ?
    # Concentrations: ?
    # TODO: compile data
    # sources:
    #   - https://idst.inl.gov

    # property: osmotic coefficients
    # Salts: ?
    # Concentrations: ?
    # TODO: compile data
    # sources:
    #   - https://idst.inl.gov
    pass


def _load_jpcrd_data() -> list[BenchmarkEntry]:
    # SOLUTION DATA
    # property: mean activity coefficient
    # Salts: ?
    # Concentrations: ?
    # TODO: compile data
    # sources:
    #   - Standard Reference Data Publications 1964-1984, NBS, Sauerwein J. C., Dalton, G. R. "Properties Index" (p. 94)
    #   - Theoretical Mean Activity Coefficients of Strong Electrolytes in Aqueous Solutions from 0 to 100C -
    #     Walter J. Hamer. NSRDS-NBS 24, 271p. (1968).
    #   - Thermodynamic Properties of Aqueous Sodium Chloride Solutions — Kenneth S. Pitzer, J. Christopher Peiper, and
    #     R. H. Busey. J Phys Chem Ref Data 13, 1(1984).

    # property: density
    # Salts: ?
    # Concentrations: ?
    # TODO: compile data
    # sources:
    #   - Standard Reference Data Publications 1964-1984, NBS, Sauerwein J. C., Dalton, G. R. "Properties Index" (p. 94)

    # property: dielectric constant
    # Salts: ?
    # Concentrations: ?
    # TODO: compile data
    # sources:
    #   - Standard Reference Data Publications 1964-1984, NBS, Sauerwein J. C., Dalton, G. R. "Properties Index" (p. 94)

    # property: electrical resistivity/conductivity
    # Salts: ?
    # Concentrations: ?
    # TODO: compile data
    # sources:
    #   - Standard Reference Data Publications 1964-1984, NBS, Sauerwein J. C., Dalton, G. R. "Properties Index" (p. 94)

    # property: osmotic coefficients
    # Salts: ?
    # Concentrations: ?
    # TODO: compile data
    # sources:
    #   - Standard Reference Data Publications 1964-1984, NBS, Sauerwein J. C., Dalton, G. R. "Properties Index" (p. 94)

    # property: dielectric constant
    # Salts: ?
    # Concentrations: ?
    # TODO: compile data
    # sources:
    #   - Standard Reference Data Publications 1964-1984, NBS, Sauerwein J. C., Dalton, G. R. "Properties Index" (p. 94)
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
