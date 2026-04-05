from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .genesets import as_geneset_mapping
from .stats_utils import bh_fdr, empirical_signed_pvalues


@dataclass(frozen=True)
class PreparedCategorySets:
    """Compact geneset representation for repeated category-score evaluation."""

    terms: tuple[str, ...]
    member_indices: tuple[np.ndarray, ...]
    matched_genes: tuple[tuple[str, ...], ...]


def prepare_category_sets(genes, geneset_resource) -> PreparedCategorySets:
    """Prepare geneset membership indices for category-score enrichment."""

    genes = [str(gene) for gene in np.asarray(genes, dtype=object).reshape(-1).tolist()]
    gene_to_pos = {gene: index for index, gene in enumerate(genes)}
    geneset_mapping = as_geneset_mapping(geneset_resource)

    terms: list[str] = []
    member_indices: list[np.ndarray] = []
    matched_genes: list[tuple[str, ...]] = []
    for term, members in geneset_mapping.items():
        hits = tuple(gene for gene in members if gene in gene_to_pos)
        if not hits:
            continue
        positions = np.asarray([gene_to_pos[gene] for gene in hits], dtype=np.int32)
        terms.append(str(term))
        member_indices.append(positions)
        matched_genes.append(hits)

    if not terms:
        raise ValueError("The selected geneset resource does not overlap the ranked gene universe.")

    return PreparedCategorySets(
        terms=tuple(terms),
        member_indices=tuple(member_indices),
        matched_genes=tuple(matched_genes),
    )


def category_scores_many(scores, prepared: PreparedCategorySets) -> np.ndarray:
    """Return one mean category score per geneset and score vector."""

    score_matrix = np.asarray(scores, dtype=float)
    if score_matrix.ndim == 1:
        score_matrix = score_matrix.reshape(-1, 1)
    if score_matrix.ndim != 2:
        raise ValueError("scores must be a one- or two-dimensional array.")

    out = np.zeros((len(prepared.terms), score_matrix.shape[1]), dtype=float)
    for term_index, indices in enumerate(prepared.member_indices):
        out[term_index, :] = score_matrix[indices, :].mean(axis=0)
    return out


def ensemble_table(
    observed_scores,
    null_scores,
    prepared: PreparedCategorySets,
) -> pd.DataFrame:
    """Build the standard ensemble-enrichment output table."""

    observed = np.asarray(observed_scores, dtype=float).reshape(-1)
    nulls = np.asarray(null_scores, dtype=float)
    if nulls.ndim != 2 or nulls.shape[0] != observed.shape[0]:
        raise ValueError("Null category scores must be a 2D array aligned to the observed scores.")

    null_mean = nulls.mean(axis=1)
    null_sd = nulls.std(axis=1, ddof=1)
    safe_sd = np.where(null_sd == 0, np.finfo(float).eps, null_sd)
    z_score = (observed - null_mean) / safe_sd
    p_val = empirical_signed_pvalues(observed, nulls)
    fdr = bh_fdr(p_val)

    out = pd.DataFrame(
        {
            "Term": list(prepared.terms),
            "category_score": observed,
            "null_mean": null_mean,
            "null_sd": null_sd,
            "z_score": z_score,
            "p_val": p_val,
            "fdr": fdr,
            "matched_size": [int(indices.size) for indices in prepared.member_indices],
            "matched_genes": [";".join(genes) for genes in prepared.matched_genes],
        }
    )
    out["_abs_z"] = np.abs(out["z_score"].to_numpy(dtype=float))
    out = out.sort_values(["fdr", "_abs_z", "Term"], ascending=[True, False, True], kind="mergesort")
    return out.drop(columns="_abs_z").reset_index(drop=True)
