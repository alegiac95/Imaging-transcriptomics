import numpy as np

from imaging_transcriptomics.pls_backend import fit_prepared_pls1, pls_regression, prepare_pls1


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


def test_prepared_pls_backend_matches_direct_fit():
    rs = np.random.RandomState(11)
    X = rs.normal(size=(28, 10))
    y = rs.normal(size=(28, 1))

    direct = pls_regression(X, y, n_components=3, n_perm=0, n_boot=0)
    prepared = fit_prepared_pls1(prepare_pls1(X), y, n_components=3)

    np.testing.assert_allclose(prepared["x_weights"], direct["x_weights"])
    np.testing.assert_allclose(prepared["x_scores"], direct["x_scores"])
    np.testing.assert_allclose(prepared["y_scores"], direct["y_scores"])
    np.testing.assert_allclose(prepared["varexp"], direct["varexp"])


def test_reduced_prepared_pls_backend_matches_full_weights_and_variance():
    rs = np.random.RandomState(19)
    X = rs.normal(size=(32, 11))
    y = rs.normal(size=(32, 1))

    full = fit_prepared_pls1(prepare_pls1(X), y, n_components=4)
    reduced = fit_prepared_pls1(
        prepare_pls1(X),
        y,
        n_components=4,
        return_full=False,
        return_x_scores=False,
        return_x_weights=True,
    )

    np.testing.assert_allclose(reduced["x_weights"], full["x_weights"])
    np.testing.assert_allclose(reduced["varexp"], full["varexp"])
    assert "x_scores" not in reduced
