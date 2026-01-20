from inspect import cleandoc
from pathlib import Path

import numpy as np
from pyEQL_phreeqc import Phreeqc, PHRQSol
from pytest import approx


def test_load_database_internal():
    # `phreeqc.dat` is included with the package, so we don't need to specify
    # `database_directory`.
    Phreeqc(database="phreeqc.dat")


def test_load_database_external():
    # We can always load an external database by specifying
    # `database_directory`.
    Phreeqc(database="phreeqc.dat", database_directory=Path(__file__).parent)


def test_run_sample():
    # Run a generic script string and capture output
    phreeqc = Phreeqc("phreeqc.dat")

    phreeqc.run_string("""
        TITLE Example 11.--Transport and ion exchange.
        SOLUTION 0  CaCl2
            units            mmol/kgw
            temp             25.0
            pH               7.0     charge
            pe               12.5    O2(g)   -0.68
            Ca               0.6
            Cl               1.2
        SOLUTION 1  Initial solution for column
            units            mmol/kgw
            temp             25.0
            pH               7.0     charge
            pe               12.5    O2(g)   -0.68
            Na               1.0
            K                0.2
            N(5)             1.2
            END
        EXCHANGE 1
            equilibrate 1
            X                0.0011
        END
    """)

    assert phreeqc.get_component_count() == 5
    # components are in alphabetic order
    assert phreeqc.get_component(0) == "Ca"
    assert phreeqc.get_component(1) == "Cl"
    assert phreeqc.get_component(2) == "K"
    assert phreeqc.get_component(3) == "N"
    assert phreeqc.get_component(4) == "Na"

    phreeqc.run_string("""
        SELECTED_OUTPUT
            -reset false
        USER_PUNCH
            -headings    cb    H    O    Ca	Cl	K	N	Na
            10 w = TOT("water")
            20 PUNCH CHARGE_BALANCE, TOTMOLE("H"), TOTMOLE("O")
            30 PUNCH w*TOT("Ca")
            40 PUNCH w*TOT("Cl")
            50 PUNCH w*TOT("K")
            60 PUNCH w*TOT("N")
            70 PUNCH w*TOT("Na")
    """)

    phreeqc.run_string("""
        RUN_CELLS; -cells 0-1
    """)

    assert phreeqc.get_selected_output_row_count() == 3
    assert phreeqc.get_selected_output_column_count() == 8

    # We can get values at a specific index, or a slice
    assert phreeqc.output[0, 0] == "cb"
    assert phreeqc.output[0, 1:4] == ["H", "O", "Ca"]
    assert phreeqc.output[0, 5:] == ["K", "N", "Na"]

    assert np.allclose(
        phreeqc.output[1, :],
        np.array(
            [
                2.979292808179192e-18,
                111.01243360409575,
                55.50675186622646,
                0.0006000000000000017,
                0.0012000000000000005,
                0.0,
                0.0,
                0.0,
            ]
        ),
    )

    assert np.allclose(
        phreeqc.output[2],  # single indices are also fine
        np.array(
            [
                -3.395954270633993e-16,
                111.01243360154359,
                55.510351938877264,
                0.0,
                0.0,
                0.00020000000000012212,
                0.0012000000000010212,
                0.0010000000000005584,
            ]
        ),
    )


def test_run_simple():
    phreeqc = Phreeqc()
    # Note: No "SAVE SOLUTION" - we can get the results directly
    phreeqc.run_string("""
        SOLUTION 0
          temp 25.0
          units mol/kgw
          pH 7.0
          pe 8.5
          redox pe
          water 0.9970480319717386
        SELECTED_OUTPUT
          -solution true
        END
    """)
    assert phreeqc.output.shape[0] == 2  # header + 1 solution


def test_run_add_solution():
    phreeqc = Phreeqc()
    solution = PHRQSol(
        {"pH": 7.0, "pe": 8.5, "redox": "pe", "temp": 25.0, "units": "mol/kgw", "water": 0.9970480319717386}
    )
    phreeqc.add_solution(solution)
    assert len(phreeqc) == 1


def test_run_add_delete_solution():
    phreeqc = Phreeqc()
    solution = PHRQSol(
        {"pH": 7.0, "pe": 8.5, "redox": "pe", "temp": 25.0, "units": "mol/kgw", "water": 0.9970480319717386}
    )
    phreeqc.add_solution(solution)
    phreeqc.remove_solution(0)
    assert len(phreeqc) == 0


def test_run_dumpstring():
    phreeqc = Phreeqc()

    phreeqc.run_string("""
        SOLUTION 0
          temp 25.0
        REACTION 1
          CaCl2 1
          Na2CO3 1
        1 mmol
        SAVE SOLUTION 0
        END
    """)

    phreeqc.set_dump_string_on(1)

    phreeqc.run_string("""
        DUMP
          -solution 0
        END
    """)

    dump_string = phreeqc.get_dump_string()
    phreeqc.set_dump_string_on(0)

    # Due to platform specific differences, we only compare the first token
    # from each of the lines below, to the output.
    expected = cleandoc("""
        SOLUTION_RAW                 0 Solution after simulation 1.
          -temp                      25
          -pressure                  1
          -potential                 0
          -total_h                   111.01243359386
          -total_o                   55.509216797548
          -cb                        -1.2200890388064e-09
          -density                   0.99704301397679
          -viscosity                 0.89125921464527
          -viscos_0                  0.89002391825059
          -totals
            C(4)                     0.0010000000027134
            Ca                       0.0010000000010752
            Cl                       0.0019999999999998
            Na                       0.002
            O(0)                     2.8581598475304e-15
          -pH                        10.413844680873
          -pe                        7.3950560340151
          -mu                        0.0045560988049649
          -ah2o                      0.99989811240764
          -mass_water                0.999994868511
          -soln_vol                  1.0029812985132
          -total_alkalinity          0.0020000012222396
          -activities
            C(-4)                    -125.72046751407
            C(4)                     -3.4925949082566
            Ca                       -3.2726648366134
            Cl                       -2.7307519728506
            E                        -7.3950560340151
            H(0)                     -38.767826310689
            Na                       -2.730257199021
            O(0)                     -14.844435883364
          -gammas
        USE mix none
        USE reaction none
        USE reaction_temperature none
        USE reaction_pressure none
    """)

    dump_lines = dump_string.splitlines()
    expected_lines = expected.splitlines()

    assert len(dump_lines) == len(expected_lines)

    for got, exp in zip(dump_lines, expected_lines, strict=False):
        got_first = got.strip().split()[0]
        exp_first = exp.strip().split()[0]
        assert got_first == exp_first


def test_run_logstring():
    phreeqc = Phreeqc()
    phreeqc.set_log_string_on(1)
    phreeqc.run_string("""
        KNOBS
          -logfile true
        SOLUTION 0
          temp 25.0
        REACTION 1
          CaCl2 1
          Na2CO3 1
        1 mmol
        SAVE SOLUTION 0
        END
    """)

    log_string = phreeqc.get_log_string()
    phreeqc.set_log_string_on(0)

    # Due to platform specific differences, we only compare the first token
    # from each of the lines below, to the output.
    expected = cleandoc("""
               -------------------------------------------
               Beginning of initial solution calculations.
               -------------------------------------------

               Initial solution 0.

               Iterations in revise_guesses: 1

               Number of infeasible solutions: 0
               Number of basis changes: 0

               Number of iterations: 0

               -----------------------------------------
               Beginning of batch-reaction calculations.
               -----------------------------------------

               Reaction step 1.

               Overflow: (CO2)2\t1.000000e+03\t3.625196e+00\t-1
               Overflow: CO2\t4.794226e+02\t2.680719e+00\t-1
               Overflow: CaCO3\t1.000000e+03\t3.225283e+00\t-1
               Overflow: CaHCO3+\t1.000000e+03\t3.915297e+00\t-1
               Overflow: HCO3-\t1.000000e+03\t3.329016e+00\t-1
               Overflow: NaHCO3\t1.000000e+03\t3.268854e+00\t-1
               Iterations in revise_guesses: 2

               Number of infeasible solutions: 0
               Number of basis changes: 0

               Number of iterations: 15

               ------------------
               End of simulation.
               ------------------

               ------------------------------------
               Reading input data for simulation 2.
               ------------------------------------

               ---------------------------------
               End of Run after X Seconds.
               ---------------------------------
           """)

    log_lines = log_string.strip("\n").splitlines()
    expected_lines = expected.splitlines()

    assert len(log_lines) == len(expected_lines)

    for got, exp in zip(log_lines, expected_lines, strict=False):
        got = got.strip()
        exp = exp.strip()
        if exp:  # ignore blank lines
            got_first = got.split()[0]
            exp_first = exp.split()[0]

            # The `---` lines seem to have different width depending on
            # platform. In those cases just move on.
            if got_first.startswith("---") and exp_first.startswith("---"):
                continue

            assert got_first == exp_first


def test_speciation_one_solution():
    phreeqc = Phreeqc()

    phreeqc.run_string("""
        SELECTED_OUTPUT
            -reset false

        USER_PUNCH
        10 t = SYS("aq", count, name$, type$, moles)
        20 FOR i = 1 to count
        30 PUNCH name$(i), MOL(name$(i)), ACT(name$(i))
        40 NEXT i

        SOLUTION 0
          temp 25.0
          units mol/kgw
          pH 7.0
          pe 8.5
          redox pe
          water 0.9970480319717386
        END
        """)

    # heading + soln 0 data (regardless of whether we have -headings in USER_PUNCH or not)
    assert phreeqc.get_selected_output_row_count() == 2
    # [<name>, <molality>, <activity>], repeated for each (4) species
    assert phreeqc.get_selected_output_column_count() == 12


def test_speciation_two_solutions():
    phreeqc = Phreeqc()

    phreeqc.run_string("""
        SELECTED_OUTPUT
            -reset false

        USER_PUNCH
        10 t = SYS("aq", count, name$, type$, moles)
        20 FOR i = 1 to count
        30 PUNCH name$(i), MOL(name$(i)), ACT(name$(i))
        40 NEXT i

        SOLUTION 0
          temp 25.0
          units mol/kgw
          pH 7.0
          pe 8.5
          redox pe
          water 0.9970480319717386
        SAVE SOLUTION 0
        END

        SOLUTION 1
          temp 50.0
          units mol/kgw
          pH 10.0
          pe 8.5
          redox pe
          water 0.9970480319717386
        SAVE SOLUTION 1
        END
        """)

    # heading + soln 0 data + soln 1 data (regardless of whether we have -headings in USER_PUNCH or not)
    assert phreeqc.get_selected_output_row_count() == 3
    # [<name>, <molality>, <activity>], repeated for each (4) species
    assert phreeqc.get_selected_output_column_count() == 12


def test_add_solution_input():
    phreeqc = Phreeqc()

    phreeqc.add_solution(
        [
            PHRQSol(
                {"pH": 7.0, "pe": 8.5, "redox": "pe", "temp": 25.0, "units": "mol/kgw", "water": 0.9970480319717386}
            ),
            PHRQSol(
                {"pH": 10.0, "pe": 8.5, "redox": "pe", "temp": 50.0, "units": "mol/kgw", "water": 0.9970480319717386}
            ),
        ]
    )

    expected = cleandoc("""
        SELECTED_OUTPUT
            -reset false

        USER_PUNCH
        5 PUNCH CELL_NO, TOT['water'], OSMOTIC, EOL_NOTAB$
        10 t = SYS("aq", count, name$, type$, moles)
        20 FOR i = 1 to count
        30 PUNCH name$(i), MOL(name$(i)), ACT(name$(i)), DIFF_C(name$(i))
        40 NEXT i
        50 PUNCH EOL$
        60 p = SYS("phases", count, name$, type$, moles)
        70 FOR j = 1 TO count
        80 PUNCH name$(j), SI(name$(j))
        90 NEXT j

        SOLUTION 0
          pH 7.0
          pe 8.5
          redox pe
          temp 25.0
          units mol/kgw
          water 0.9970480319717386
        SAVE SOLUTION 0
        END

        SOLUTION 1
          pH 10.0
          pe 8.5
          redox pe
          temp 50.0
          units mol/kgw
          water 0.9970480319717386
        SAVE SOLUTION 1
        END""")

    assert str(phreeqc) == expected


def test_add_solution_output():
    phreeqc = Phreeqc()

    phreeqc.add_solution(
        [
            PHRQSol(
                {"pH": 7.0, "pe": 8.5, "redox": "pe", "temp": 25.0, "units": "mol/kgw", "water": 0.9970480319717386}
            ),
            PHRQSol(
                {"pH": 10.0, "pe": 8.5, "redox": "pe", "temp": 50.0, "units": "mol/kgw", "water": 0.9970480319717386}
            ),
        ]
    )

    # heading + soln 0 data + soln 1 data (regardless of whether we have -headings in USER_PUNCH or not)
    assert phreeqc.get_selected_output_row_count() == 3
    # <cell_no>, <tot_water>, <osmotic>, "\n",
    # +[<name>, <molality>, <activity>, <diff_c>] repeated for each (4) species
    # +"\n" + [<equilibrium_species>, <si>] repeated for each (3) equilibrium species
    assert phreeqc.get_selected_output_column_count() == 27


def test_kgw():
    phreeqc = Phreeqc()
    phreeqc.add_solution(
        PHRQSol({"pH": 7.0, "pe": 8.5, "redox": "pe", "temp": 25.0, "units": "mol/kgw", "water": 0.9970480319717386})
    )
    assert phreeqc[0].get_kgw() == approx(0.9970480319717386)


def test_speciate():
    phreeqc = Phreeqc()
    phreeqc.add_solution(
        PHRQSol({"pH": 7.0, "pe": 8.5, "redox": "pe", "temp": 25.0, "units": "mol/kgw", "water": 0.9970480319717386})
    )
    assert phreeqc[0].get_species_list() == ["OH-", "H+", "O2", "H2"]


def test_speciate_get_molality():
    solution = PHRQSol(
        {"pH": 7.0, "pe": 8.5, "redox": "pe", "temp": 25.0, "units": "mol/kgw", "water": 0.9970480319717386}
    )
    phreeqc = Phreeqc()
    phreeqc.add_solution(solution)

    activity = phreeqc[0].get_molality("H+")
    assert activity == approx(1.0005246407839175e-07)


def test_speciate_get_moles():
    solution = PHRQSol(
        {"pH": 7.0, "pe": 8.5, "redox": "pe", "temp": 25.0, "units": "mol/kgw", "water": 0.9970480319717386}
    )
    phreeqc = Phreeqc()
    phreeqc.add_solution(solution)

    activity = phreeqc[0].get_moles("H+")
    assert activity == approx(1.0005246407839175e-07 * 0.9970480319717386)


def test_speciate_species_moles():
    solution = PHRQSol(
        {"pH": 7.0, "pe": 8.5, "redox": "pe", "temp": 25.0, "units": "mol/kgw", "water": 0.9970480319717386}
    )
    phreeqc = Phreeqc()
    phreeqc.add_solution(solution)

    species_moles = phreeqc[0].species_moles
    expected = {
        "H+": 9.975711240328358e-08,
        "H2": 7.05855918557026e-35,
        "O2": 8.293083928522462e-25,
        "OH-": 1.009700360204178e-07,
    }
    assert set(species_moles.keys()) == set(expected.keys())
    for k, v in species_moles.items():
        assert v == approx(expected[k])


def test_species_all_props():
    solutions = [
        PHRQSol({"pH": 7.0, "pe": 8.5, "redox": "pe", "temp": 25.0, "units": "mol/kgw", "water": 0.9970480319717386}),
        PHRQSol({"pH": 10.0, "pe": 8.5, "redox": "pe", "temp": 50.0, "units": "mol/kgw", "water": 0.9970480319717386}),
    ]
    phreeqc = Phreeqc()
    phreeqc.add_solution(solutions)

    expected = {
        0: {
            "CELL_NO": 0,
            "TOT['water']": 0.9970480319717386,
            "OSMOTIC": 0.0,
            "species": {
                "H+": {"ACT": 1.0001522689856982e-07, "MOL": 1.0005246407839175e-07},
                "H2": {"ACT": 7.079457681907915e-35, "MOL": 7.079457517820301e-35},
                "O2": {"ACT": 8.317637520771417e-25, "MOL": 8.317637327985348e-25},
                "OH-": {"ACT": 1.0123126727760366e-07, "MOL": 1.0126897880811404e-07},
            },
            "eq_species": {
                "H2(g)": {"SI": -31.048923033678367},
                "H2O(g)": {"SI": -1.5028233204748613},
                "O2(g)": {"SI": -21.18764110890733},
            },
        },
        1: {
            "CELL_NO": 1,
            "TOT['water']": 0.9970480319717386,
            "OSMOTIC": 0.0,
            "species": {
                "H+": {"ACT": 1e-10, "MOL": 1.0197827284798617e-10},
                "H2": {"ACT": 5.626700758118202e-41, "MOL": 0.0},
                "O2": {"ACT": 3.659468532125681e-05, "MOL": 3.6592332322612637e-05},
                "OH-": {"ACT": 0.0005473549298676792, "MOL": 0.0005585111533405827},
            },
            "eq_species": {
                "H2(g)": {"SI": -37.113374354647966},
                "H2O(g)": {"SI": -0.9158658033784839},
                "O2(g)": {"SI": -1.4065718231105961},
            },
        },
    }

    assert len(phreeqc) == len(expected)
    for solution_index, solution_props in expected.items():
        props = phreeqc[solution_index]._get_calculated_props()
        assert set(props.keys()) == set(expected[solution_index].keys())
        for k in expected[solution_index]:
            if k not in ("species", "eq_species"):
                assert props[k] == approx(expected[solution_index][k])

        for k in solution_props["species"]:
            for prop in ("ACT", "MOL"):
                assert np.isclose(props["species"][k][prop], expected[solution_index]["species"][k][prop])

        for k in solution_props["eq_species"]:
            for prop in ("SI",):
                assert np.isclose(props["eq_species"][k][prop], expected[solution_index]["eq_species"][k][prop])


def test_speciate_get_activity():
    solution = PHRQSol(
        {"pH": 7.0, "pe": 8.5, "redox": "pe", "temp": 25.0, "units": "mol/kgw", "water": 0.9970480319717386}
    )
    phreeqc = Phreeqc()
    phreeqc.add_solution(solution)

    activity = phreeqc[0].get_activity("H+")
    assert activity == approx(1.0001522689856982e-07)


def test_get_osmotic_coefficient():
    solution = PHRQSol(
        {"pH": 7.0, "pe": 8.5, "redox": "pe", "temp": 25.0, "units": "mol/kgw", "water": 0.9970480319717386}
    )
    phreeqc = Phreeqc()
    phreeqc.add_solution(solution)

    osmotic_coefficient = phreeqc[0].get_osmotic_coefficient()
    assert osmotic_coefficient == approx(0.0)


def test_get_diffusion_coefficient():
    solution = PHRQSol(
        {"pH": 7.0, "pe": 8.5, "redox": "pe", "temp": 25.0, "units": "mol/kgw", "water": 0.9970480319717386}
    )
    phreeqc = Phreeqc()
    phreeqc.add_solution(solution)

    diffusion_coefficient = phreeqc[0].get_diffusion_coefficient("H+")
    assert diffusion_coefficient == approx(9.31e-09)
