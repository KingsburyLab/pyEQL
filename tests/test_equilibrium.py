"""
Tests of pyEQL.equilibrium module

"""
import numpy as np

from pyEQL.equilibrium import alpha


def test_alpha():
    # for each case, test pH << pKa, pH=pKa, pH >> pKa

    # monoprotic acid
    assert np.isclose(alpha(0, 3, [7]), 1, atol=1e-3)
    assert np.isclose(alpha(1, 3, [7]), 0, atol=1e-3)

    assert np.isclose(alpha(0, 7, [7]), 0.5, atol=1e-3)
    assert np.isclose(alpha(1, 7, [7]), 0.5, atol=1e-3)

    assert np.isclose(alpha(0, 10, [7]), 0, atol=1e-3)
    assert np.isclose(alpha(1, 10, [7]), 1, atol=1e-3)

    # carbonic acid, diprotic
    assert np.isclose(alpha(0, 3, [6.35, 10.33]), 1, atol=1e-3)
    assert np.isclose(alpha(1, 3, [6.35, 10.33]), 0, atol=1e-3)
    assert np.isclose(alpha(2, 3, [6.35, 10.33]), 0, atol=1e-3)

    assert np.isclose(alpha(0, 6.35, [6.35, 10.33]), 0.5, atol=1e-3)
    assert np.isclose(alpha(1, 6.35, [6.35, 10.33]), 0.5, atol=1e-3)

    assert np.isclose(alpha(0, 7.35, [6.35, 10.33]), 0.1, atol=1e-2)
    # test sorting
    assert np.isclose(alpha(0, 7.35, [10.33, 6.35]), 0.1, atol=1e-2)
    assert np.isclose(alpha(1, 7.35, [6.35, 10.33]), 0.9, atol=1e-2)
    assert np.isclose(alpha(2, 7.35, [6.35, 10.33]), 0, atol=1e-3)

    assert np.isclose(alpha(1, 10.33, [6.35, 10.33]), 0.5, atol=1e-3)
    assert np.isclose(alpha(2, 10.33, [6.35, 10.33]), 0.5, atol=1e-3)

    assert np.isclose(alpha(0, 12.33, [6.35, 10.33]), 0, atol=1e-3)
    assert np.isclose(alpha(1, 12.33, [6.35, 10.33]), 0.01, atol=1e-2)
    assert np.isclose(alpha(2, 12.33, [6.35, 10.33]), 0.99, atol=1e-2)
