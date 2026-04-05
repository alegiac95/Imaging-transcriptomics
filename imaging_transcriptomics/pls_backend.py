from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class PreparedPLS1:
    """Prepared X-side data for repeated PLS-1 fits with different Y vectors."""

    X: np.ndarray
    X0: np.ndarray
    x_mean: np.ndarray
    x_ss_total: float


def _get_mask(X, Y):
    return np.logical_not(np.logical_or(np.all(np.isnan(X), axis=1), np.all(np.isnan(Y), axis=1)))


def prepare_pls1(X) -> PreparedPLS1:
    X = np.asarray(X, dtype=float)
    x_mean = X.mean(axis=0, keepdims=True)
    X0 = X - x_mean
    return PreparedPLS1(
        X=X,
        X0=X0,
        x_mean=x_mean,
        x_ss_total=float(np.sum(X0 * X0)),
    )


def _simpls_prepared_pls1_reduced(
    prepared: PreparedPLS1,
    Y,
    n_components: int | None = None,
    *,
    return_x_scores: bool,
    return_x_weights: bool,
):
    """Specialized fast path for repeated PLS-1 fits.

    This reduced path is used during permutation testing, where the code only
    needs explained variance plus a subset of the latent variables.
    """

    X = prepared.X
    X0 = prepared.X0
    y = np.asarray(Y, dtype=float).reshape(-1)
    if n_components is None:
        n_components = min(len(X) - 1, X.shape[1])

    y_mean = float(np.mean(y))
    y0 = y - y_mean
    y_ss_total = float(np.dot(y0, y0))
    cov = X0.T @ y0

    x_scores = np.zeros((X.shape[0], n_components), dtype=float) if return_x_scores else None
    x_weights = np.zeros((X.shape[1], n_components), dtype=float) if return_x_weights else None
    varexp = np.zeros(n_components, dtype=float)
    basis = np.zeros((X.shape[1], n_components), dtype=float)

    for comp in range(n_components):
        cov_norm = float(np.linalg.norm(cov))
        if cov_norm == 0:
            break
        right_vector = cov / cov_norm
        ti = X0 @ right_vector
        normti = float(np.linalg.norm(ti))
        if normti == 0:
            break

        weight = right_vector / normti
        ti = ti / normti
        if return_x_weights:
            x_weights[:, comp] = weight
        if return_x_scores:
            x_scores[:, comp] = ti

        x_loading = X0.T @ ti
        qi = float(np.dot(y0, ti))
        varexp[comp] = (qi * qi) / y_ss_total if y_ss_total else 0.0

        vi = x_loading
        if comp > 0:
            previous = basis[:, :comp]
            vi = vi - (previous @ (previous.T @ vi))
        vi_norm = float(np.linalg.norm(vi))
        if vi_norm == 0:
            break
        basis[:, comp] = vi / vi_norm

        active_basis = basis[:, : comp + 1]
        cov = cov - (active_basis @ (active_basis.T @ cov))

    result = {"varexp": varexp}
    if return_x_weights:
        result["x_weights"] = x_weights
    if return_x_scores:
        result["x_scores"] = x_scores
    return result


def _simpls_prepared_pls1(prepared: PreparedPLS1, Y, n_components: int | None = None, *, return_full: bool = True):
    """Small SIMPLS implementation specialized for the repo's PLS-1 use case."""

    X = prepared.X
    X0 = prepared.X0
    Y = np.asarray(Y, dtype=float)
    if Y.ndim == 1:
        Y = Y.reshape(-1, 1)
    if n_components is None:
        n_components = min(len(X) - 1, X.shape[1])

    y_mean = Y.mean(axis=0, keepdims=True)
    Y0 = Y - y_mean
    y_ss_total = float(np.sum(Y0 * Y0))
    cov = X0.T @ Y0

    x_loadings = np.zeros((X.shape[1], n_components), dtype=float)
    y_loadings = np.zeros((Y.shape[1], n_components), dtype=float)
    x_scores = np.zeros((X.shape[0], n_components), dtype=float)
    y_scores = np.zeros((X.shape[0], n_components), dtype=float)
    x_weights = np.zeros((X.shape[1], n_components), dtype=float)
    basis = np.zeros((X.shape[1], n_components), dtype=float)

    for comp in range(n_components):
        cov_norm = np.linalg.norm(cov)
        if cov_norm == 0:
            break
        right_vector = cov / cov_norm
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
        if comp > 0:
            previous = basis[:, :comp]
            vi = vi - (previous @ (previous.T @ vi))
        vi_norm = np.linalg.norm(vi)
        if vi_norm == 0:
            break
        vi = vi / vi_norm
        basis[:, [comp]] = vi

        active_basis = basis[:, : comp + 1]
        cov = cov - (active_basis @ (active_basis.T @ cov))

    pctvar = [
        np.sum(x_loadings ** 2, axis=0) / prepared.x_ss_total,
        np.sum(y_loadings ** 2, axis=0) / y_ss_total if y_ss_total else np.zeros(n_components, dtype=float),
    ]
    varexp = np.asarray(pctvar[1], dtype=float)
    if not return_full:
        return dict(
            x_weights=x_weights,
            x_scores=x_scores,
            y_scores=y_scores,
            varexp=varexp,
        )

    beta = x_weights @ y_loadings.T
    beta = np.vstack([y_mean - (prepared.x_mean @ beta), beta])

    mse = np.zeros((2, n_components + 1), dtype=float)
    mse[0, 0] = prepared.x_ss_total
    mse[1, 0] = y_ss_total
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
        varexp=varexp,
    )


def fit_prepared_pls1(
    prepared: PreparedPLS1,
    Y,
    *,
    n_components: int | None = None,
    return_full: bool = True,
    return_x_scores: bool = True,
    return_x_weights: bool = True,
):
    """Fit a prepared PLS-1 model.

    ``return_full=True`` preserves the original compatibility-heavy output.
    ``return_full=False`` switches to a reduced permutation-friendly path that
    only materializes the requested arrays.
    """

    if return_full:
        return _simpls_prepared_pls1(prepared, Y, n_components=n_components, return_full=True)
    return _simpls_prepared_pls1_reduced(
        prepared,
        Y,
        n_components=n_components,
        return_x_scores=return_x_scores,
        return_x_weights=return_x_weights,
    )


def fit_prepared_pls1_chunk(
    prepared: PreparedPLS1,
    Y,
    *,
    n_components: int,
    return_x_weights: bool = False,
) -> dict[str, np.ndarray]:
    """Fit repeated reduced PLS-1 models for a block of response vectors.

    This provides a backend seam for permutation-heavy workflows without
    forcing the workflow layer to manage one fit at a time.
    """

    response_block = np.asarray(Y, dtype=float)
    if response_block.ndim == 1:
        response_block = response_block.reshape(-1, 1)
    if response_block.ndim != 2:
        raise ValueError("Y must be a one- or two-dimensional array.")

    n_perm = response_block.shape[1]
    varexp = np.zeros((n_perm, n_components), dtype=float)
    x_weights = (
        np.zeros((n_perm, prepared.X.shape[1], n_components), dtype=float)
        if return_x_weights
        else None
    )

    for index in range(n_perm):
        fit = _simpls_prepared_pls1_reduced(
            prepared,
            response_block[:, index],
            n_components=n_components,
            return_x_scores=False,
            return_x_weights=return_x_weights,
        )
        varexp[index, :] = np.asarray(fit["varexp"], dtype=float)[:n_components]
        if return_x_weights and x_weights is not None:
            x_weights[index, :, :] = np.asarray(fit["x_weights"], dtype=float)[:, :n_components]

    out = {"varexp": varexp}
    if return_x_weights and x_weights is not None:
        out["x_weights"] = x_weights
    return out


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
    prepared = prepare_pls1(X[mask])
    result = fit_prepared_pls1(prepared, Y[mask], n_components=n_components, return_full=True)

    # Re-expand score arrays to the original sample axis for compatibility.
    x_scores = np.full((len(Y), result["x_scores"].shape[1]), np.nan, dtype=float)
    y_scores = np.full((len(Y), result["y_scores"].shape[1]), np.nan, dtype=float)
    x_scores[mask] = result["x_scores"]
    y_scores[mask] = result["y_scores"]
    result["x_scores"] = x_scores
    result["y_scores"] = y_scores
    return result
