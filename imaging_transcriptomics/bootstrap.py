import logging
import logging.config
import yaml
from pathlib import Path
from itertools import product
from multiprocessing import Pool
from functools import partial

from scipy.stats import spearmanr
import numpy
import numpy as np
from tqdm import tqdm
from pyls import pls_regression

from .genes import GeneResults

cfg_file_path = Path(__file__).parent / "log_config.yaml"
with open(cfg_file_path, "r") as config_file:
    log_cfg = yaml.safe_load(config_file.read())

logging.config.dictConfig(log_cfg)
logger = logging.getLogger("bootstrapping")
logger.setLevel(logging.DEBUG)


def correlate(corr1, corr2):
    """Return correlation similar to MATLAB corr function.

    :param corr1: first element to correlate.
    :param corr2: second element to correlate.
    """
    return np.corrcoef(np.hstack((corr1, corr2)), rowvar=False)[0, 1:]


def bootstrap_pls(x, y, y_perm, dim, iterations=1_000):
    """Return a coefficient of determination (R^2) and a p value from the bootstrapping on the PLS regression.

    The function calculates for each of the N ``iterations`` the value of the coefficient of determination  on the
    permuted data and determines how many times this is at least equal to the coefficient of determination with the
    original data.

    :param x: X data for the PLS regression (in our case the zscore of the gene expression data)
    :param y: Y data for the PLS regression (in our case the zscore of the average of the ROIs)
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
        pls_results = pls_regression(x, y, n_components=component, n_perm=0, n_boot=0)
        exp_var = pls_results.get("varexp")
        temp = 100 * exp_var.cumsum(axis=0)
        R_squared = temp[component - 1]
        R_sq = np.zeros(iterations)
        for i in tqdm(
            range(1000),
            desc=f"Bootstrapping on PLS component {component}",
            unit=" iterations",
        ):
            y_data = y_perm[:, i].reshape(41, 1)
            _result = pls_regression(
                x, y_data, n_components=component, n_perm=0, n_boot=0
            )
            _exp_var = 100 * np.cumsum(_result.get("varexp"))
            R_sq[i] = _exp_var[component - 1]
        R_boot[component - 1] = R_squared
        p_boot[component - 1] = (
            float(len(R_sq[numpy.nonzero(R_sq >= R_squared)])) / iterations
        )
    logger.debug(
        "Computed components p value(s) and coefficient(s) of determination: \n R: %s \n p: %s",
        R_boot,
        p_boot,
    )
    return R_boot, p_boot


def bootstrap_genes(x_data, y, n_components, y_norm, genes, n_iterations=1000):
    """Run a bootstrap permuting the genes as well as the data, to assess reliability of the results.



    :param x_data: data to regress, provided as z scores (e.g., gene expression data).
    :param y: data to regress against, provided as z scores (e.g. average values from ROIs).
    :param n_components: number of components to use for the PLS regression.
    :param y_norm: original y data (before z scoring).
    :param genes: labels of the genes regressed against (e.g., the labels of the y data rows).
    :param n_iterations: number of iterations for the bootstrapping and permutations.
    :return gene_results: GeneResults data structure with all the results of the bootstrapping and the original
    results as well.
    """
    n_genes = 15_633
    gene_index = np.array(list(range(1, n_genes + 1)))
    results = pls_regression(x_data, y, n_components=n_components, n_boot=0, n_perm=0)
    r1 = correlate(
        results.get("x_scores").reshape(41, n_components), y_norm.reshape(41, 1)
    )
    logger.debug("Correlation between original data and regression scores: %s", r1)
    weights = results.get("x_weights")
    gene_results = GeneResults(
        n_components, dim1=weights.shape[0], dim2=weights.shape[1]
    )
    scores = results.get("x_scores")
    for i in range(r1.size):
        if r1[i] < 0:
            weights[:, i] *= -1
            scores[:, i] *= -1
    for idx in range(1, n_components + 1):
        x = np.argsort(weights[:, idx - 1], kind="mergesort")[::-1]
        gene_results.original_results.set_result_values(
            idx,
            np.sort(weights[:, idx - 1], kind="mergesort")[::-1],
            x,
            genes[x],
            gene_index[x],
        )

    # Main genes bootstrap
    for iteration in tqdm(
        range(n_iterations), desc="Bootstrapping gene list", unit=" iterations"
    ):
        my_resample = np.random.choice(41, size=41)
        x_perm = x_data[my_resample, :]
        y_perm = y[my_resample].reshape(41, 1)
        results = pls_regression(
            x_perm, y_perm, n_components=n_components, n_perm=0, n_boot=0
        )
        _weights = results.get("x_weights")
        for component in range(1, n_components + 1):
            __temp = _weights[:, component - 1]
            __new_weights = __temp[gene_results.original_results.index[component - 1]]
            __t_genes = np.array(
                gene_results.original_results.pls_weights[component - 1]
            )
            __correlation = correlate(
                __t_genes.reshape(15633, 1), __new_weights.reshape(15633, 1)
            )
            if __correlation < 0:
                __new_weights *= -1
            gene_results.boot_results.pls_weights_boot[component - 1][
                :, component - 1, iteration
            ] = __new_weights
    return gene_results


def spearman_op(idx, permuted, y_data):
    return spearmanr(permuted[:, idx[0]], y_data[:, idx[1]])[0]


def bootstrap_correlation(x_data, y_data, permuted, labels, n_iterations=1000):
    """
    Bootstrap the results using pearson correlation.

    :param x_data: imaging data
    :param y_data: gene expression data
    :param permuted: permuted matrix of imaging data
    :param labels: labels of the genes (original order)
    :return gene_results: GeneResults class with correlation results.
    """
    gene_results = GeneResults(n_comp=1, dim1=1, dim2=y_data.shape[1])
    n_genes = 15633
    pool = Pool()
    # original correlation
    corr_ = np.zeros(y_data.shape[1])
    for gene in range(n_genes):
        corr_[gene], _ = spearmanr(x_data, y_data[:, gene])
    # order results and set indexes in gene results
    gene_results.original_results.pls_weights = np.sort(corr_, kind="mergesort")[::-1]
    __idx = np.argsort(corr_, kind="mergesort")[::-1]
    gene_results.original_results.pls_gene = labels[__idx]
    gene_results.original_results.gene_id = __idx
    # bootstrap
    __res = np.zeros((15633, 1000))
    _iter = product(range(n_iterations), range(n_genes))
    for ind, result in tqdm(enumerate(pool.imap(partial(
            spearman_op, permuted=permuted, y_data=y_data),
                                                _iter, chunksize=25_000
                                                )),
                            desc="Bootstrapping correlation",
                            unit="iterations"):
        __res.flat[ind] = result
    gene_results.boot_results.pls_weights_boot = __res
    return gene_results
