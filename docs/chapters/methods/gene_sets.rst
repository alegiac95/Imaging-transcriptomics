=========================
GSEA and ORA methodology
=========================

The package supports two gene-set analyses on top of the gene-level outputs:

- preranked GSEA
- over-representation analysis (ORA)

GSEA
----

GSEA is run on the ranked gene list from a correlation run or from one PLS
component at a time.

The current implementation:

- builds a deterministic preranked table
- resolves ties by adding tiny stable offsets inside exact score ties
- computes observed enrichment scores with ``gseapy``
- computes external-null enrichment scores by rerunning the same ranking logic
  on the package's bootstrap or permutation nulls
- derives ``NES`` from same-sign null normalization
- derives nominal ``p_val`` from the external null enrichment scores
- derives ``fdr`` using a GSEA-style FDR calculation on normalized null scores

This means the reported ``NES``, ``p_val``, and ``fdr`` are calibrated to the
package's external nulls rather than to GSEApy's internal permutation engine.

ORA
---

ORA works on a selected subset of genes rather than on the full ranked list.
The current selection rule is:

- keep genes with raw ``p_value <= threshold``
- split them into positive and negative sets
- run hypergeometric enrichment separately for ``up`` and ``down``

Each ORA table includes:

- ``overlap_size``
- ``set_size``
- ``selected_size``
- ``universe_size``
- ``enrichment_ratio``
- ``odds_ratio``
- ``odds_ratio_ci_low`` and ``odds_ratio_ci_high``
- ``p_value``
- ``fdr``
- ``overlap_genes``

Confidence intervals
--------------------

ORA reports a 95 percent confidence interval for the odds ratio. The current
implementation uses a log-odds interval with a Haldane-Anscombe correction
when a contingency-table cell is zero.

When to use which
-----------------

- use GSEA when you want a ranking-based test that uses all genes
- use ORA when you want to threshold genes and interpret discrete hit lists
- use both when you want a threshold-free view and a thresholded view of the
  same result
