from __future__ import annotations

from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd

from ._logging import get_logger
from .corr_stats import (
    ranked_gene_expression,
    spearman_correlation_bootstrap,
    spearman_correlation_matrix,
)
from .genes import CorrGenes, GeneResults
from .gsea_utils import (
    gsea_style_fdr,
    make_prerank_table,
    nominal_pvalues_from_nulls,
    normalize_enrichment_nulls,
    normalize_enrichment_scores,
    run_prerank,
)
from .genesets import get_geneset
from .ora import ora_from_gene_table

logger = get_logger(__name__)

# Backward-compatible aliases for tests and internal callers that still import
# the older private helper names from this module.
_ranked_gene_expression = ranked_gene_expression
_spearman_correlation_matrix = spearman_correlation_matrix
_spearman_correlation_bootstrap = spearman_correlation_bootstrap


class CorrAnalysis:
    """Store correlation analysis results and optional GSEA output."""

    def __init__(self, n_iterations=1000, n_genes=None):
        self.gene_results = GeneResults("corr", n_iter=n_iterations, n_genes=n_genes)

    def bootstrap_correlation(self, imaging_data, permuted_imaging, gene_exp, gene_labels):
        """Run the original and bootstrapped correlation analyses."""

        assert isinstance(self.gene_results.results, CorrGenes)
        logger.info("Calculating correlation on original data.")
        ranked_genes = ranked_gene_expression(gene_exp)
        self.gene_results.results.corr[:, :] = spearman_correlation_matrix(imaging_data, ranked_genes).T

        logger.info("Calculating correlation on permuted data.")
        self.gene_results.results.boot_corr[:, :] = spearman_correlation_bootstrap(permuted_imaging, ranked_genes)
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
        rnk = make_prerank_table(gene_list, self.gene_results.results.corr[0, :])
        gsea_results = run_prerank(
            gseapy,
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
            rnk = make_prerank_table(gene_list, self.gene_results.results.boot_corr[:, index])
            gsea_res = run_prerank(
                gseapy,
                rnk,
                gene_set,
                max_size=gene_limit,
                permutation_num=0,
                no_plot=True,
                outdir=None,
                seed=1234,
            )
            boot_es[:, index] = gsea_res.res2d.ES.values

        boot_nes = normalize_enrichment_scores(origin_es, boot_es)
        boot_nes_null = normalize_enrichment_nulls(origin_es, boot_es)
        p_val = nominal_pvalues_from_nulls(origin_es, boot_es)
        p_corr = gsea_style_fdr(boot_nes, boot_nes_null)

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

    def ora(self, gene_set="lake", outdir=None, p_threshold=0.05):
        """Perform ORA on positively and negatively associated genes."""

        assert isinstance(self.gene_results.results, CorrGenes)
        logger.info("Performing ORA.")
        gene_table = pd.DataFrame(
            {
                "gene": self.gene_results.results.genes[:, 0],
                "score": self.gene_results.results.corr[0, :],
                "p_value": self.gene_results.results.pval[0, :],
            }
        )
        ora_tables = ora_from_gene_table(
            gene_table,
            gene_set=gene_set,
            score_column="score",
            p_threshold=p_threshold,
        )
        if outdir is not None:
            logger.info("Saving ORA results.")
            outdir = Path(outdir)
            assert outdir.exists()
            for direction, table in ora_tables.items():
                table.to_csv(outdir / f"ora_corr_{direction}.tsv", index=False, sep="\t")
        return ora_tables
