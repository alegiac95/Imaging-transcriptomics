import numpy
import numpy as np
from tqdm import tqdm


def pls_regression(X, Y, n_components, **kwargs):
    """Placeholder while i understand how to incorporate the library from Ross Markello."""
    return {}


def bootstrap_pls(x, y, y_perm, dim, iterations=1_000):
    """Return a coefficient of determination (R^2) and a p value from the bootstrapping on the PLS regression.

    The function calculates for each of the N ``iterations`` the value of the coefficient of determination  on the
    permuted data and determines how many times this is at least equal to the coefficient of determination with the
    original data.

    :param x: X data for the PLS regression (in our case the zscore of the average of the ROIs)
    :param y: Y data for the PLS regression (in our case the zscore of the gene expression data)
    :param y_perm: matrix with the permuted Y data. This should have a number of rows equal to the original data and
        number of columns equal to the number of iterations.
    :param int dim: number of PLS components to use for the analysis.
    :param int iterations: number of bootstrap iterations.

    :return R_boot: array of dimension ``dim`` with the coefficient of determination for each PLS component.
    :return p_boot: array of dimension ``dim`` with the p value for each PLS component.
    """
    R_boot = np.zeros(dim)
    p_boot = np.zeros(dim)
    for component in range(1, dim + 1):
        pls_results = pls_regression(x, y,
                                     n_components=component,
                                     n_perm=0,
                                     n_boot=0)
        exp_var = pls_results.get("varexp")
        temp = 100 * exp_var.cumsum(axis=0)
        R_squared = temp[component - 1]
        R_sq = np.zeros(iterations)
        for i in tqdm(range(1000), desc=f"Bootstrapping on PLS component {component}", unit="iteration"):
            y_data = y_perm[:, i].reshape(41, 1)
            _result = pls_regression(x, y_data,
                                     n_components=component,
                                     n_perm=0,
                                     n_boot=0)
            _exp_var = 100 * np.cumsum(_result.get("varexp"))
            R_sq[i] = _exp_var[component - 1]
        R_boot[component - 1] = R_squared
        p_boot[component - 1] = float(len(R_sq[numpy.nonzero(R_sq >= R_squared)])) / iterations
    return R_boot, p_boot


def bootstrap_genes():
    pass
