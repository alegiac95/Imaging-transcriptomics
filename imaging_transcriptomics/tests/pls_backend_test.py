import numpy as np

from imaging_transcriptomics.pls_backend import pls_regression


def test_local_pls_backend_runs_and_returns_expected_shapes():
    rs = np.random.RandomState(1)
    X = rs.normal(size=(24, 9))
    y = rs.normal(size=(24, 1))

    result = pls_regression(X, y, n_components=3, n_perm=0, n_boot=0)

    assert result["x_weights"].shape == (9, 3)
    assert result["x_scores"].shape == (24, 3)
    assert result["y_scores"].shape == (24, 3)
    assert result["varexp"].shape == (3,)
    assert np.all(result["varexp"] >= 0)


def test_local_pls_backend_is_deterministic():
    rs = np.random.RandomState(7)
    X = rs.normal(size=(30, 12))
    y = rs.normal(size=(30, 1))

    first = pls_regression(X, y, n_components=4, n_perm=0, n_boot=0)
    second = pls_regression(X, y, n_components=4, n_perm=0, n_boot=0)

    np.testing.assert_allclose(first["x_weights"], second["x_weights"])
    np.testing.assert_allclose(first["x_scores"], second["x_scores"])
    np.testing.assert_allclose(first["y_scores"], second["y_scores"])
    np.testing.assert_allclose(first["varexp"], second["varexp"])
