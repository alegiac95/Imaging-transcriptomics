import numpy
import numpy as np
from tqdm import tqdm


def pls_regression(X, Y, n_components, **kwargs):
    """Placeholder while i understand how to incorporate the library."""
    return {}


def bootstrap_pls(X, Y, Y_perm, dim, iterations=1_000):
    R_boot = np.zeros(dim)
    p_boot = np.zeros(dim)
    for component in range(1, dim + 1):
        pls_results = pls_regression(X, Y,
                                     n_components=component,
                                     n_perm=0,
                                     n_boot=0)
        exp_var = pls_results.get("varexp")
        temp = 100 * exp_var.cumsum(axis=0)
        R_squared = temp[component - 1]
        R_sq = np.zeros(iterations)
        for i in tqdm(range(1000), desc=f"Bootstrapping on PLS component {component}", unit="iteration"):
            y_data = Y_perm[:, i].reshape(41, 1)
            _result = pls_regression(X, y_data,
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
