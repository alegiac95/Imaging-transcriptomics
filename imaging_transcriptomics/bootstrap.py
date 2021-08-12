import numpy
import numpy as np
from tqdm import tqdm
from pyls import pls_regression

from .genes import GeneResults


def correlate(corr1, corr2):
    """Return correlation similar to MATLAB corr function.

    :param corr1: first element to correlate.
    :param corr2: second element to correlate.
    """
    return np.corrcoef(np.hstack(
            corr1, corr2
    ), rowvar=False)[0, 1:]


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


def bootstrap_genes(x, y, n_components, x_norm, genes, n_iterations=1000):
    """Run a bootstrap permuting the genes as well as the data, to assess reliability of the results.



    :param x: data to regress, provided as z scores (e.g., average values from ROIs).
    :param y: data to regress against, provided as z scores (e.g. gene expression data).
    :param n_components: number of components to use for the PLS regression.
    :param x_norm: original x data (before z scoring).
    :param genes: labels of the genes regressed against (e.g., the labels of the y data rows).
    :param n_iterations: number of iterations for the bootstrapping and permutations.
    :return gene_results: GeneResults data structure with all the results of the bootstrapping and the original
    results as well.
    """
    n_genes = 15_633
    gene_index = np.array(list(range(1, n_genes+1)))
    results = pls_regression(x, y, n_components=n_components, n_boot=0, n_perm=0)
    r1 = correlate(results.get("x_scores"), x_norm)
    weights = results.get("x_weights")
    gene_results = GeneResults(n_components, dim1=weights.shape[0], dim2=weights.shape[1])
    scores = results.get("x_scores")
    for i in range(r1.size):
        if r1[i] < 0:
            weights[:, i] *= -1
            scores[:, i] *= -1
    for idx in range(1, n_components+1):
        x = np.argsort(weights[:, idx-1], kind='mergesort')[::-1]
        gene_results.original_results.set_result_values(idx,
                                                        np.sort(weights[:, idx-1], kind='mergesort')[::-1],
                                                        x,
                                                        genes[x],
                                                        gene_index[x])

    # Main genes bootstrap
    for iteration in tqdm(range(n_iterations), desc="Bootstrapping gene list"):
        my_resample = np.random.choice(41, size=41)
        x_perm = x[my_resample, :]
        y_perm = y[my_resample, :]
        results = pls_regression(x_perm, y_perm, n_components=n_components)
        _weights = results.get("x_weights")
        for component in range(1, n_components+1):
            __temp = _weights[:, component-1]
            __new_weights = __temp[gene_results.original_results.x[component-1]]
            __correlation = correlate(
                gene_results.original_results.pls_weights[component-1].reshape(15633, 1),
                __new_weights.reshape(15633, 1)
            )
            if __correlation < 0:
                __new_weights *= -1
            gene_results.boot_results.pls_weights_boot[component-1][:, component-1, iteration] = __new_weights
    return gene_results
