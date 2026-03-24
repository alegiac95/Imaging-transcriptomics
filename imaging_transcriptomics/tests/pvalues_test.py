import importlib
import sys
import types
from pathlib import Path

import numpy as np
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests


PACKAGE_ROOT = Path(__file__).resolve().parents[1]


def _import_without_package_init(module_name):
    package_name = "imaging_transcriptomics"
    if package_name not in sys.modules:
        package = types.ModuleType(package_name)
        package.__path__ = [str(PACKAGE_ROOT)]
        sys.modules[package_name] = package
    return importlib.import_module(module_name)


genes_module = _import_without_package_init("imaging_transcriptomics.genes")
pls_module = _import_without_package_init("imaging_transcriptomics.pls")


def _make_small_pls_genes(n_iter):
    pls_genes = genes_module.PLSGenes(1, n_iter=n_iter)
    pls_genes.n_genes = 3
    pls_genes.orig = genes_module.OrigPLS(1, 3)
    pls_genes.boot = genes_module.BootPLS(1, 3, n_iter=n_iter)
    return pls_genes


def test_pls_gene_compute_uses_two_sided_pvalues_and_unsorted_bh():
    pls_genes = _make_small_pls_genes(n_iter=3)
    pls_genes.orig.genes[0, :] = np.array(["A", "B", "C"], dtype=object)
    pls_genes.orig.weights[0, :] = np.array([10.0, -10.0, 1.0])
    pls_genes.boot.weights[0, :, :] = np.array(
        [[9.0, 10.0, 11.0],
         [-11.0, -10.0, -9.0],
         [0.0, 1.0, 2.0]]
    )

    pls_genes.compute()

    sorted_z = np.array([10.0, 1.0, -10.0])
    expected_pval = 2 * norm.sf(np.abs(sorted_z))
    expected_fdr = multipletests(expected_pval, method="fdr_bh")[1]

    np.testing.assert_allclose(pls_genes.boot.z_score[0, :], sorted_z)
    np.testing.assert_allclose(pls_genes.boot.pval[0, :], expected_pval)
    np.testing.assert_allclose(pls_genes.boot.pval_corr[0, :], expected_fdr)
    assert pls_genes.boot.genes[0, :].tolist() == ["A", "C", "B"]


def test_corr_gene_sorting_keeps_scores_genes_and_bootstraps_aligned():
    corr_genes = genes_module.CorrGenes(n_iter=3)
    corr_genes.n_genes = 3
    corr_genes.corr = np.array([[0.2, -0.3, 0.1]])
    corr_genes.genes = np.array([["A"], ["B"], ["C"]], dtype=object)
    corr_genes.boot_corr = np.array(
        [[1.0, 1.1, 1.2],
         [2.0, 2.1, 2.2],
         [3.0, 3.1, 3.2]]
    )
    corr_genes.pval = np.zeros((1, 3))
    corr_genes.pval_corr = np.zeros((1, 3))

    corr_genes.sort_genes()

    np.testing.assert_array_equal(corr_genes._index, np.array([0, 2, 1]))
    np.testing.assert_allclose(corr_genes.corr[0, :], np.array([0.2, 0.1, -0.3]))
    assert corr_genes.genes[:, 0].tolist() == ["A", "C", "B"]
    np.testing.assert_allclose(
        corr_genes.boot_corr,
        np.array([[1.0, 1.1, 1.2],
                  [3.0, 3.1, 3.2],
                  [2.0, 2.1, 2.2]])
    )
