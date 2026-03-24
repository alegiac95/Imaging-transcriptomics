from __future__ import annotations

import numpy as np


def resid_yscores(x_scores, y_scores, copy: bool = True):
    """Orthogonalize Y scores with respect to preceding X scores.

    This follows the standard SIMPLS score orthogonalization step.
    """

    x_scores = np.array(x_scores, dtype=float, copy=False)
    y_scores = np.array(y_scores, dtype=float, copy=copy)

    for comp in range(x_scores.shape[1]):
        ui = y_scores[:, [comp]]
        for _ in range(2):
            for j in range(comp):
                tj = x_scores[:, [j]]
                ui = ui - ((tj.T @ ui) * tj)
        y_scores[:, [comp]] = ui

    return y_scores


def _get_mask(X, Y):
    return np.logical_not(np.logical_or(np.all(np.isnan(X), axis=1), np.all(np.isnan(Y), axis=1)))


def _first_svd_component(crosscov):
    """Return the leading left singular vector / value / right singular vector."""

    U, singular_values, Vt = np.linalg.svd(np.asarray(crosscov, dtype=float), full_matrices=False)
    return Vt[:1, :].T, np.diag(singular_values[:1]), U[:, [0]]


def _simpls(X, Y, n_components: int | None = None):
    """Small SIMPLS implementation for the repo's PLS-1 use case."""

    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    if Y.ndim == 1:
        Y = Y.reshape(-1, 1)
    if n_components is None:
        n_components = min(len(X) - 1, X.shape[1])

    X0 = X - X.mean(axis=0, keepdims=True)
    Y0 = Y - Y.mean(axis=0, keepdims=True)
    cov = X0.T @ Y0

    x_loadings = np.zeros((X.shape[1], n_components), dtype=float)
    y_loadings = np.zeros((Y.shape[1], n_components), dtype=float)
    x_scores = np.zeros((X.shape[0], n_components), dtype=float)
    y_scores = np.zeros((X.shape[0], n_components), dtype=float)
    x_weights = np.zeros((X.shape[1], n_components), dtype=float)
    basis = np.zeros((X.shape[1], n_components), dtype=float)

    for comp in range(n_components):
        _, singular_values, right_vector = _first_svd_component(cov)
        ti = X0 @ right_vector
        normti = np.linalg.norm(ti)
        if normti == 0:
            break

        x_weights[:, [comp]] = right_vector / normti
        ti = ti / normti
        x_scores[:, [comp]] = ti

        x_loadings[:, [comp]] = X0.T @ ti
        qi = Y0.T @ ti
        y_loadings[:, [comp]] = qi
        y_scores[:, [comp]] = Y0 @ qi

        vi = x_loadings[:, [comp]]
        for _ in range(2):
            for j in range(comp):
                vj = basis[:, [j]]
                vi = vi - ((vj.T @ vi) * vj)
        vi_norm = np.linalg.norm(vi)
        if vi_norm == 0:
            break
        vi = vi / vi_norm
        basis[:, [comp]] = vi

        cov = cov - (vi @ (vi.T @ cov))
        if comp > 0:
            previous = basis[:, :comp]
            cov = cov - (previous @ (previous.T @ cov))

    y_scores = resid_yscores(x_scores, y_scores)
    beta = x_weights @ y_loadings.T
    beta = np.vstack([Y.mean(axis=0) - (X.mean(axis=0) @ beta), beta])

    pctvar = [
        np.sum(x_loadings ** 2, axis=0) / np.sum(X0 ** 2),
        np.sum(y_loadings ** 2, axis=0) / np.sum(Y0 ** 2),
    ]
    mse = np.zeros((2, n_components + 1), dtype=float)
    mse[0, 0] = np.sum(np.abs(X0) ** 2)
    mse[1, 0] = np.sum(np.abs(Y0) ** 2)
    X0_recon = np.zeros_like(X0)
    Y0_recon = np.zeros_like(Y0)
    for i in range(n_components):
        X0_recon = x_scores[:, : i + 1] @ x_loadings[:, : i + 1].T
        Y0_recon = x_scores[:, : i + 1] @ y_loadings[:, : i + 1].T
        mse[0, i + 1] = np.sum(np.abs(X0 - X0_recon) ** 2)
        mse[1, i + 1] = np.sum(np.abs(Y0 - Y0_recon) ** 2)
    mse /= len(X)

    t2 = np.sum((np.abs(x_scores) ** 2) / np.var(x_scores, axis=0, ddof=1), axis=1)
    x_residuals = X0 - X0_recon
    y_residuals = Y0 - Y0_recon

    return dict(
        x_weights=x_weights,
        x_loadings=x_loadings,
        y_loadings=y_loadings,
        x_scores=x_scores,
        y_scores=y_scores,
        x_residuals=x_residuals,
        y_residuals=y_residuals,
        beta=beta,
        pctvar=pctvar,
        mse=mse,
        t2=t2,
        varexp=np.asarray(pctvar[1], dtype=float),
    )


def pls_regression(
    X,
    Y,
    *,
    n_components: int | None = None,
    n_perm: int = 0,
    n_boot: int = 0,
    seed=None,
    verbose: bool = True,
    **kwargs,
):
    """Small local SIMPLS regression backend for imaging transcriptomics.

    The imaging-transcriptomics codepath only needs the direct SIMPLS fit with
    `n_perm=0` and `n_boot=0`. We intentionally reject unsupported resampling
    modes here so the backend stays small and explicit.
    """

    del seed, verbose, kwargs
    if n_perm != 0 or n_boot != 0:
        raise NotImplementedError(
            "The local PLS backend only implements the direct SIMPLS fit with n_perm=0 and n_boot=0."
        )

    X = np.array(X, dtype=float, copy=True)
    Y = np.array(Y, dtype=float, copy=True)
    if Y.ndim == 1:
        Y = Y.reshape(-1, 1)
    mask = _get_mask(X, Y)
    result = _simpls(X[mask], Y[mask], n_components=n_components)

    # Re-expand score arrays to the original sample axis for compatibility.
    x_scores = np.full((len(Y), result["x_scores"].shape[1]), np.nan, dtype=float)
    y_scores = np.full((len(Y), result["y_scores"].shape[1]), np.nan, dtype=float)
    x_scores[mask] = result["x_scores"]
    y_scores[mask] = result["y_scores"]
    result["x_scores"] = x_scores
    result["y_scores"] = y_scores
    return result
