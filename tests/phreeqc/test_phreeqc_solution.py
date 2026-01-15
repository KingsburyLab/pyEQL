from pyEQL_phreeqc import PHRQSol


def test_create_solution():
    PHRQSol(
        {
            "Cl": "4.011842831773806",
            "Na": "4.011842831773806",
            "pH": 7.0,
            "pe": 8.5,
            "redox": "pe",
            "temp": 25.0,
            "units": "mol/kgw",
            "water": 0.9970480319717386,
        }
    )
