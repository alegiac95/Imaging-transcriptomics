from __future__ import annotations

import logging
from collections import OrderedDict

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from ._logging import get_logger
from .genes import CorrGenes, GeneResults
from .genesets import get_geneset

logger = get_logger("genes")
logger.setLevel(logging.DEBUG)


def _rank_vector(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(values.shape[0], dtype=float)
    ranks[order] = np.arange(values.shape[0], dtype=float)
    return ranks


def _rank_columns(matrix: np.ndarray) -> np.ndarray:
    order = np.argsort(matrix, axis=0, kind="mergesort")
    ranks = np.empty(order.shape, dtype=float)
    template = np.arange(matrix.shape[0], dtype=float)[:, None]
    np.put_along_axis(ranks, order, template, axis=0)
    return ranks


def _standardize_columns(matrix: np.ndarray) -> np.ndarray:
    centered = matrix - matrix.mean(axis=0, keepdims=True)
    scale = matrix.std(axis=0, ddof=1, keepdims=True)
    scale[scale == 0] = 1.0
    return centered / scale


def _standardize_vector(values: np.ndarray) -> np.ndarray:
    centered = values - values.mean()
    scale = values.std(ddof=1)
    if scale == 0:
        return centered
    return centered / scale


def _ranked_gene_expression(gene_exp: np.ndarray) -> np.ndarray:
    return _standardize_columns(_rank_columns(np.asarray(gene_exp, dtype=float)))


def _spearman_correlation_matrix(imaging_data: np.ndarray, ranked_genes: np.ndarray) -> np.ndarray:
    ranked_img = _standardize_vector(_rank_vector(np.asarray(imaging_data, dtype=float).reshape(-1)))
    denom = max(ranked_img.shape[0] - 1, 1)
    return (ranked_genes.T @ ranked_img.reshape(-1, 1)) / denom


def _spearman_correlation_bootstrap(permuted_imaging: np.ndarray, ranked_genes: np.ndarray) -> np.ndarray:
    ranked_perm = _standardize_columns(_rank_columns(np.asarray(permuted_imaging, dtype=float)))
    denom = max(ranked_perm.shape[0] - 1, 1)
    return (ranked_genes.T @ ranked_perm) / denom


class CorrAnalysis:
    """Store correlation analysis results and optional GSEA output."""

    def __init__(self, n_iterations=1000, n_genes=None):
        self.gene_results = GeneResults("corr", n_iter=n_iterations, n_genes=n_genes)

    def bootstrap_correlation(self, imaging_data, permuted_imaging, gene_exp, gene_labels):
        """Run the original and bootstrapped correlation analyses."""

        assert isinstance(self.gene_results.results, CorrGenes)
        logger.info("Calculating correlation on original data.")
        ranked_genes = _ranked_gene_expression(gene_exp)
        self.gene_results.results.corr[:, :] = _spearman_correlation_matrix(imaging_data, ranked_genes).T

        logger.info("Calculating correlation on permuted data.")
        self.gene_results.results.boot_corr[:, :] = _spearman_correlation_bootstrap(permuted_imaging, ranked_genes)
        self.gene_results.results.genes = gene_labels
        self.gene_results.results.sort_genes()
        self.gene_results.results.compute_pval()
        return

    def gsea(self, gene_set="lake", outdir=None, gene_limit=1500, n_perm=1_000):  # pragma: no cover
        """Perform GSEA on the correlation ranking."""

        assert isinstance(self.gene_results.results, CorrGenes)
        logger.info("Performing GSEA.")
        try:
            import gseapy
        except ImportError as exc:  # pragma: no cover - optional dependency
            raise ImportError("gseapy is required to run GSEA analyses.") from exc

        gene_set = get_geneset(gene_set)
        gene_list = list(self.gene_results.results.genes[:, 0].tolist())
        rnk = pd.DataFrame(zip(gene_list, self.gene_results.results.corr[0, :]))
        gsea_results = gseapy.prerank(
            rnk,
            gene_set,
            max_size=gene_limit,
            outdir=None,
            permutation_num=0,
            seed=1234,
        )
        origin_es = gsea_results.res2d.ES.to_numpy()
        boot_es = np.zeros((origin_es.shape[0], n_perm))
        for index in range(n_perm):
            rnk = pd.DataFrame(zip(gene_list, self.gene_results.results.boot_corr[:, index]))
            gsea_res = gseapy.prerank(
                rnk,
                gene_set,
                max_size=gene_limit,
                permutation_num=0,
                no_plot=True,
                outdir=None,
                seed=1234,
            )
            boot_es[:, index] = gsea_res.res2d.ES.values

        boot_nes_mean = np.mean(boot_es, axis=1)
        boot_nes_std = np.std(boot_es, axis=1)
        boot_nes = (boot_nes_mean - origin_es) / boot_nes_std
        p_val = np.zeros(origin_es.shape[0], dtype=float)
        for index in range(origin_es.shape[0]):
            count = (
                np.sum(boot_es[index, :] >= origin_es[index])
                if origin_es[index] >= 0
                else np.sum(boot_es[index, :] <= origin_es[index])
            )
            p_val[index] = (count + 1) / (n_perm + 1)
        _, p_corr, _, _ = multipletests(p_val, method="fdr_bh", is_sorted=False)

        out_df = pd.DataFrame.from_dict(
            OrderedDict(
                Term=gsea_results.res2d.Term.values.tolist(),
                es=gsea_results.res2d.ES.values,
                nes=boot_nes,
                p_val=p_val,
                fdr=p_corr,
            )
        )
        if outdir is not None:
            logger.info("Saving GSEA results.")
            outdir = Path(outdir)
            assert outdir.exists()
            out_df.to_csv(outdir / "gsea_corr_results.tsv", index=False, sep="\t")
