from __future__ import annotations

from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor, as_completed
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
from .genesets import resolve_geneset_resource
from .ora import ora_from_gene_table

logger = get_logger(__name__)

# Backward-compatible aliases for tests and internal callers that still import
# the older private helper names from this module.
_ranked_gene_expression = ranked_gene_expression
_spearman_correlation_matrix = spearman_correlation_matrix
_spearman_correlation_bootstrap = spearman_correlation_bootstrap


def _progress_step(total: int) -> int:
    """Choose a coarse logging interval for long permutation loops."""

    return max(1, total // 10)


def _log_permutation_progress(completed: int, total: int, *, step: int) -> None:
    """Emit periodic progress updates for long correlation permutation loops."""

    if total < 1:
        return
    if completed == total or completed % step == 0:
        logger.info("Processed %d/%d correlation permutations.", completed, total)


def _corr_chunk_bounds(total: int, n_genes: int, n_jobs: int, target_bytes: int = 64 * 1024 * 1024) -> list[tuple[int, int]]:
    """Split correlation permutations into memory-friendly chunks."""

    bytes_per_perm = max(1, int(n_genes) * np.dtype(float).itemsize)
    chunk_size = max(1, target_bytes // bytes_per_perm)
    chunk_size = min(total, chunk_size)
    return [
        (start, min(total, start + chunk_size))
        for start in range(0, total, chunk_size)
    ]


class CorrAnalysis:
    """Run spatial correlation analysis and hold its downstream outputs.

    The object keeps the gene-wise correlation statistics together with any
    optional enrichment results produced from the ranked gene table.
    """

    def __init__(self, n_iterations=1000, n_genes=None, *, store_boot_corr: bool = True, n_jobs: int = 1):
        """Create an empty correlation analysis container."""

        self.n_jobs = max(1, int(n_jobs))
        self.gene_results = GeneResults(
            "corr",
            n_iter=n_iterations,
            n_genes=n_genes,
            store_boot_corr=store_boot_corr,
        )

    def bootstrap_correlation(self, imaging_data, permuted_imaging, gene_exp, gene_labels):
        """Run the original and bootstrapped correlation analyses."""

        assert isinstance(self.gene_results.results, CorrGenes)
        total_permutations = int(permuted_imaging.shape[1])
        progress_step = _progress_step(total_permutations)
        logger.info("Calculating correlation on original data.")
        ranked_genes = ranked_gene_expression(gene_exp)
        self.gene_results.results.corr[:, :] = spearman_correlation_matrix(imaging_data, ranked_genes).T
        self.gene_results.results.genes = gene_labels

        logger.info(
            "Calculating correlation on permuted data (%d permutations, %d job%s).",
            total_permutations,
            self.n_jobs,
            "" if self.n_jobs == 1 else "s",
        )
        bounds = _corr_chunk_bounds(total_permutations, self.gene_results.results.n_genes, self.n_jobs)
        if self.n_jobs == 1 or len(bounds) == 1:
            completed = 0
            for start, end in bounds:
                boot_chunk = spearman_correlation_bootstrap(permuted_imaging[:, start:end], ranked_genes)
                self.gene_results.results.accumulate_boot_corr(boot_chunk, start=start)
                completed += end - start
                _log_permutation_progress(completed, total_permutations, step=progress_step)
        else:
            max_workers = min(self.n_jobs, len(bounds))

            def _corr_chunk(bounds_: tuple[int, int]):
                start, end = bounds_
                return start, end, spearman_correlation_bootstrap(permuted_imaging[:, start:end], ranked_genes)

            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = [executor.submit(_corr_chunk, bound) for bound in bounds]
                completed = 0
                for future in as_completed(futures):
                    start, end, boot_chunk = future.result()
                    self.gene_results.results.accumulate_boot_corr(boot_chunk, start=start)
                    completed += end - start
                    _log_permutation_progress(completed, total_permutations, step=progress_step)
        self.gene_results.results.compute_pval()
        self.gene_results.results.sort_genes()
        return

    def gsea(
        self,
        gene_set="lake",
        outdir=None,
        gene_limit=1500,
        n_perm=1_000,
        geneset_organism: str = "Human",
    ):  # pragma: no cover
        """Run preranked GSEA on the correlation-based gene ranking.

        The observed enrichment score is taken from the ranked correlation
        vector, while NES, nominal p-values, and q-values are recalculated from
        the external permutation nulls stored in ``boot_corr``.
        """

        assert isinstance(self.gene_results.results, CorrGenes)
        if self.gene_results.results.boot_corr is None:
            raise RuntimeError("Correlation GSEA requires stored permutation nulls. Re-run with GSEA enabled.")
        logger.info("Performing GSEA.")
        try:
            import gseapy
        except ImportError as exc:  # pragma: no cover - optional dependency
            raise ImportError("gseapy is required to run GSEA analyses.") from exc

        gene_set = resolve_geneset_resource(gene_set, organism=geneset_organism)
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

    def ora(self, gene_set="lake", outdir=None, p_threshold=0.05, geneset_organism: str = "Human"):
        """Run ORA on positive and negative correlation tails separately."""

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
            geneset_organism=geneset_organism,
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
