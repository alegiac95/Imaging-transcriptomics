import numpy as np
from scipy.stats import spearmanr

from imaging_transcriptomics.corr import _spearman_correlation_bootstrap, _spearman_correlation_matrix, _ranked_gene_expression


def test_fast_spearman_matches_reference():
    rng = np.random.default_rng(1234)
    imaging = rng.normal(size=12)
    genes = rng.normal(size=(12, 5))

    ranked_genes = _ranked_gene_expression(genes)
    fast = _spearman_correlation_matrix(imaging, ranked_genes).reshape(-1)
    ref = np.array([spearmanr(imaging, genes[:, index]).statistic for index in range(genes.shape[1])])
    np.testing.assert_allclose(fast, ref, atol=1e-12)


def test_fast_spearman_bootstrap_matches_reference():
    rng = np.random.default_rng(1234)
    permuted = rng.normal(size=(12, 7))
    genes = rng.normal(size=(12, 4))

    ranked_genes = _ranked_gene_expression(genes)
    fast = _spearman_correlation_bootstrap(permuted, ranked_genes)
    ref = np.column_stack(
        [
            np.array([spearmanr(permuted[:, col], genes[:, row]).statistic for row in range(genes.shape[1])])
            for col in range(permuted.shape[1])
        ]
    )
    np.testing.assert_allclose(fast, ref, atol=1e-12)
