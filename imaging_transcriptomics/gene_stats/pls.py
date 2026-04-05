from __future__ import annotations

from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import zscore

from .._logging import get_logger
from ..ensemble import category_scores_many, ensemble_table as build_ensemble_table, prepare_category_sets
from ..gsea_utils import (
    enrichment_scores_many,
    gsea_style_fdr,
    make_prerank_table,
    nominal_pvalues_from_nulls,
    normalize_enrichment_nulls,
    normalize_enrichment_scores,
    prepare_prerank_genesets,
    result_column,
    result_terms,
    run_prerank,
)
from ..genesets import resolve_geneset_resource
from ..pls_backend import pls_regression
from ..stats_utils import bh_fdr, empirical_signed_pvalues, max_t_fwer_abs

logger = get_logger("genes")


def correlate_pls_scores(scores: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Correlate each PLS score vector with the input regional imaging values."""

    stacked = np.hstack((np.asarray(scores, dtype=float), np.asarray(values, dtype=float).reshape(-1, 1)))
    return np.corrcoef(stacked, rowvar=False)[0, 1:]


def rowwise_corrsign(
    reference: np.ndarray,
    candidate: np.ndarray,
    reference_centered: np.ndarray,
    reference_ss: np.ndarray,
) -> np.ndarray:
    """Return sign flips that align candidate rows to reference rows by correlation."""

    centered = candidate - candidate.mean(axis=1, keepdims=True)
    denom = np.sqrt(reference_ss * np.sum(centered * centered, axis=1))
    corr = np.divide(
        np.sum(reference_centered * centered, axis=1),
        denom,
        out=np.zeros(reference.shape[0], dtype=float),
        where=denom != 0,
    )
    return np.where(corr < 0, -1.0, 1.0)


class OrigPLS:
    """Hold original per-component PLS gene rankings and summary statistics."""

    def __init__(self, n_components, n_genes):
        """Allocate arrays for one original PLS gene-ranking result."""

        self.n_components = n_components
        self.weights = np.zeros((n_components, n_genes))
        self.genes = np.zeros((n_components, n_genes), dtype=object)
        self.index = np.zeros((n_components, n_genes), dtype=np.int32)
        self.zscored = np.zeros((n_components, n_genes))


class BootPLS:
    """Hold permutation-derived PLS gene weights and correction outputs."""

    def __init__(self, n_components, n_genes, n_iter=1000, store_weights: bool = True):
        """Allocate arrays for one permutation-derived PLS gene-stat bundle."""

        self.n_components = n_components
        self.n_iter = n_iter
        self.weights = np.zeros((n_components, n_genes, n_iter), dtype=np.float32) if store_weights else None
        self.genes = np.zeros((n_components, n_genes), dtype=object)
        self.weights_sorted = np.zeros((n_components, n_genes))
        self.std = np.zeros((n_components, n_genes))
        self._z_score = np.zeros((n_components, n_genes))
        self.pval = np.zeros((n_components, n_genes))
        self.pval_corr = np.zeros((n_components, n_genes))
        self.pval_fwer = np.zeros((n_components, n_genes))

    @property
    def z_score(self):
        """Return the sorted per-gene z-scores for each PLS component."""

        return self._z_score


class PLSGenes:
    """Store original and permutation-based gene statistics for PLS results."""

    def __init__(self, n_components, n_iter=1000, n_genes=None, store_weights: bool = True):
        """Create the storage backing one PLS gene-analysis workflow."""

        self.n_genes = int(n_genes) if n_genes is not None else 15633
        self.n_components = n_components
        self.n_iter = n_iter
        self.store_weights = bool(store_weights)
        self.orig = OrigPLS(n_components, self.n_genes)
        self.boot = BootPLS(n_components, self.n_genes, n_iter=n_iter, store_weights=self.store_weights)
        self._orig_centered = None
        self._orig_ss = None

    def prepare_from_fit(self, fit_result, scan_data, gene_labels):
        """Align, sort, and z-score original PLS gene weights per component."""

        weights = np.asarray(fit_result.get("x_weights"), dtype=float).copy()
        scores = np.asarray(fit_result.get("x_scores"), dtype=float).copy()
        score_corr = correlate_pls_scores(scores, scan_data)
        for component in range(score_corr.size):
            if score_corr[component] < 0:
                weights[:, component] *= -1
                scores[:, component] *= -1
        for component in range(self.n_components):
            sort_index = np.argsort(weights[:, component], kind="mergesort")[::-1]
            self.orig.index[component, :] = sort_index
            self.orig.genes[component, :] = gene_labels[:, 0][sort_index]
            self.orig.weights[component, :] = weights[:, component][sort_index]
            self.orig.zscored[component, :] = zscore(self.orig.weights[component, :], axis=0, ddof=1)
        self._orig_centered = self.orig.weights - self.orig.weights.mean(axis=1, keepdims=True)
        self._orig_ss = np.sum(self._orig_centered * self._orig_centered, axis=1)
        return

    def store_permuted_weights(self, iteration: int, x_weights):
        """Store one permuted PLS fit after aligning it to the original gene order."""

        if self._orig_centered is None or self._orig_ss is None:
            raise RuntimeError("Call prepare_from_fit() before storing permuted weights.")
        if self.boot.weights is None:
            raise RuntimeError("Permutation weights were not allocated for this PLS run.")
        weights = np.asarray(x_weights, dtype=float).T
        reordered = np.take_along_axis(weights, self.orig.index, axis=1)
        sign = rowwise_corrsign(self.orig.weights, reordered, self._orig_centered, self._orig_ss)
        self.boot.weights[:, :, iteration] = reordered * sign.reshape(-1, 1)
        return

    def boot_genes(self, imaging_data, permuted_imaging, scan_data, gene_exp, gene_labels):
        """Legacy helper to refit PLS for every permutation and store gene weights."""

        logger.info("Performing bootstrapping of the genes.")
        result = pls_regression(
            gene_exp,
            imaging_data.reshape(imaging_data.shape[0], 1),
            n_components=self.n_components,
            n_boot=0,
            n_perm=0,
        )
        self.prepare_from_fit(result, scan_data, gene_labels)
        if permuted_imaging.shape[1] != self.boot.weights.shape[2]:
            raise ValueError("The number of bootstrapped permutations does not match the configured iteration count.")
        for iteration in range(self.boot.weights.shape[2]):
            permuted_values = permuted_imaging[:, iteration]
            iteration_result = pls_regression(
                gene_exp,
                permuted_values.reshape(permuted_values.shape[0], 1),
                n_components=self.n_components,
                n_boot=0,
                n_perm=0,
            )
            self.store_permuted_weights(iteration, iteration_result.get("x_weights"))
        return

    def compute(self, n_jobs: int = 1):
        """Compute sorted PLS gene statistics from original and permuted weights."""

        logger.info("Calculating statistics.")
        if self.boot.weights is None:
            raise RuntimeError("Permutation weights are required to compute PLS gene statistics.")

        max_workers = min(max(1, int(n_jobs)), self.n_components)

        def _compute_component(component: int):
            boot_weights = np.asarray(self.boot.weights[component, :, :], dtype=float)
            orig_weights = np.asarray(self.orig.weights[component, :], dtype=float)
            orig_genes = np.asarray(self.orig.genes[component, :], dtype=object)
            std = boot_weights.std(axis=1, ddof=1)
            safe_std = np.where(std == 0, np.finfo(float).eps, std)
            zscores = orig_weights / safe_std
            indices = np.argsort(zscores, kind="mergesort")[::-1]
            raw_p = empirical_signed_pvalues(orig_weights, boot_weights)
            raw_fwer = max_t_fwer_abs(orig_weights, boot_weights)
            fdr = bh_fdr(raw_p)
            return {
                "component": component,
                "std": std,
                "weights_sorted": orig_weights[indices],
                "zscore_sorted": zscores[indices],
                "genes_sorted": orig_genes[indices],
                "p_sorted": raw_p[indices],
                "fdr_sorted": fdr[indices],
                "fwer_sorted": raw_fwer[indices],
            }

        if max_workers == 1 or self.n_components == 1:
            outputs = [_compute_component(component) for component in range(self.n_components)]
        else:
            outputs = []
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = [executor.submit(_compute_component, component) for component in range(self.n_components)]
                for future in as_completed(futures):
                    outputs.append(future.result())

        for output in outputs:
            component = output["component"]
            self.boot.std[component, :] = output["std"]
            self.boot.weights_sorted[component, :] = output["weights_sorted"]
            self.boot.z_score[component, :] = output["zscore_sorted"]
            self.boot.genes[component, :] = output["genes_sorted"]
            self.boot.pval[component, :] = output["p_sorted"]
            self.boot.pval_corr[component, :] = output["fdr_sorted"]
            self.boot.pval_fwer[component, :] = output["fwer_sorted"]
        return

    def gsea(
        self,
        gene_set="lake",
        outdir=None,
        gene_limit=1500,
        n_iter=1000,
        geneset_organism: str = "Human",
    ):
        """Run preranked GSEA on the PLS gene ranking for each component."""

        assert isinstance(self.orig, OrigPLS)
        assert isinstance(self.boot, BootPLS)
        if self.boot.weights is None:
            raise RuntimeError("PLS GSEA requires stored permutation gene weights. Re-run with GSEA enabled.")
        logger.info("Performing GSEA.")
        try:
            import gseapy
        except ImportError as exc:  # pragma: no cover - optional dependency
            raise ImportError("gseapy is required to run GSEA analyses.") from exc
        gene_set = resolve_geneset_resource(gene_set, organism=geneset_organism)
        for component in range(self.n_components):
            gene_list = list(self.orig.genes[component, :])
            rnk = make_prerank_table(gene_list, self.orig.zscored[component, :])
            gsea_results = run_prerank(
                gseapy,
                rnk,
                gene_set,
                max_size=gene_limit,
                outdir=None,
                seed=1234,
                permutation_num=0,
            )
            term_order = result_terms(gsea_results.res2d)
            prepared_sets = prepare_prerank_genesets(gene_list, gene_set, term_order=term_order)
            origin_es = result_column(gsea_results.res2d, "ES", "es").astype(float)
            boot_scores = zscore(
                np.asarray(self.boot.weights[component, :, :n_iter], dtype=float),
                axis=0,
                ddof=1,
            )
            boot_es = enrichment_scores_many(boot_scores, prepared_sets)
            nes = normalize_enrichment_scores(origin_es, boot_es)
            nes_null = normalize_enrichment_nulls(origin_es, boot_es)
            p_val = nominal_pvalues_from_nulls(origin_es, boot_es)
            p_corr = gsea_style_fdr(nes, nes_null)
            out_data = OrderedDict()
            out_data["Term"] = term_order
            out_data["es"] = origin_es
            out_data["nes"] = nes
            out_data["p_val"] = p_val
            out_data["fdr"] = p_corr
            out_data["genest_size"] = result_column(gsea_results.res2d, "geneset_size", "genest_size")
            out_data["matched_size"] = result_column(gsea_results.res2d, "matched_size")
            out_data["matched_genes"] = result_column(gsea_results.res2d, "matched_genes")
            out_data["ledge_genes"] = result_column(gsea_results.res2d, "ledge_genes")
            out_df = pd.DataFrame.from_dict(out_data)
            if outdir is not None:
                logger.info("Saving GSEA results.")
                output_dir = Path(outdir)
                assert output_dir.exists()
                out_df.to_csv(output_dir / f"gsea_pls{component + 1}_results.tsv", index=False, sep="\t")

    def ora(self, gene_set="lake", outdir=None, p_threshold=0.05, geneset_organism: str = "Human"):
        """Run ORA on positive and negative component gene tails separately."""

        assert isinstance(self.orig, OrigPLS)
        assert isinstance(self.boot, BootPLS)
        logger.info("Performing ORA.")
        from .. import genes as genes_module

        results: list[dict[str, pd.DataFrame]] = []
        for component in range(self.n_components):
            gene_table = pd.DataFrame(
                {
                    "gene": self.boot.genes[component, :],
                    "zscore": self.boot.z_score[component, :],
                    "p_value": self.boot.pval[component, :],
                }
            )
            ora_tables = genes_module.ora_from_gene_table(
                gene_table,
                gene_set=gene_set,
                geneset_organism=geneset_organism,
                score_column="zscore",
                p_threshold=p_threshold,
            )
            if outdir is not None:
                logger.info("Saving ORA results for PLS component %d.", component + 1)
                output_dir = Path(outdir)
                assert output_dir.exists()
                for direction, table in ora_tables.items():
                    table.to_csv(
                        output_dir / f"ora_pls{component + 1}_{direction}.tsv",
                        index=False,
                        sep="\t",
                    )
            results.append(ora_tables)
        return results

    def ensemble(
        self,
        gene_set="lake",
        outdir=None,
        n_iter=1000,
        geneset_organism: str = "Human",
    ) -> list[pd.DataFrame]:
        """Run phenotype-null ensemble enrichment on each retained PLS component."""

        assert isinstance(self.orig, OrigPLS)
        assert isinstance(self.boot, BootPLS)
        if self.boot.weights is None:
            raise RuntimeError("PLS ensemble enrichment requires stored permutation gene weights. Re-run with enrichment enabled.")
        logger.info("Performing ensemble enrichment.")
        gene_set = resolve_geneset_resource(gene_set, organism=geneset_organism)
        outputs: list[pd.DataFrame] = []
        for component in range(self.n_components):
            gene_list = list(self.orig.genes[component, :])
            prepared = prepare_category_sets(gene_list, gene_set)
            observed = category_scores_many(self.orig.zscored[component, :], prepared)[:, 0]
            boot_scores = zscore(
                np.asarray(self.boot.weights[component, :, :n_iter], dtype=float),
                axis=0,
                ddof=1,
            )
            null_scores = category_scores_many(boot_scores, prepared)
            out_df = build_ensemble_table(observed, null_scores, prepared)
            outputs.append(out_df)
            if outdir is not None:
                logger.info("Saving ensemble enrichment results.")
                output_dir = Path(outdir)
                assert output_dir.exists()
                out_df.to_csv(output_dir / f"ensemble_pls{component + 1}_results.tsv", index=False, sep="\t")
        return outputs
