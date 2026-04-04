====================
Enrichment methods
====================

This page summarizes the three enrichment families exposed by the toolbox:

- ``ensemble``
- ``gsea``
- ``ora``

All three start from a gene-wise signal derived from the imaging analysis, but
they summarize that signal differently and they use different null models.
That distinction matters in imaging transcriptomics, where the inferential
problem is shaped as much by the spatial structure of the phenotype as by the
gene set itself. For a broader discussion of those concerns, see
`Fulcher, Arnatkeviciute and Fornito (2021)
<https://doi.org/10.1038/s41467-021-22862-1>`_.

Overview
--------

``ensemble``
   Tests a category-level score directly against category scores obtained from
   null phenotypes. This is the default enrichment backend in the package.

``gsea``
   Tests whether genes from a category accumulate near one end of a ranked gene
   list by using a running enrichment score.

``ora``
   Tests whether a category is over-represented among a thresholded subset of
   selected genes.

The workflow-facing explanation lives in :doc:`/chapters/workflows/enrichment`.
This page instead focuses on the statistical objects that each method computes.

Ensemble
--------

Category score
~~~~~~~~~~~~~~

Ensemble enrichment starts from one observed gene-wise score vector. In
practice, that score is:

- a gene-map association score for ``corr``
- a component-specific gene weight or aligned component score for ``pls``

For a term :math:`t` with member genes :math:`G_t`, the toolbox computes a
category score as the mean gene score over the genes in that term:

.. math::

   S_t = \frac{1}{|G_t|} \sum_{g \in G_t} s_g

where :math:`s_g` is the observed gene-wise score.

Null calibration
~~~~~~~~~~~~~~~~

The same category score is recomputed across the null phenotype ensemble. If
:math:`s_g^{(b)}` is the gene-wise score from null phenotype :math:`b`, then:

.. math::

   S_t^{(b)} = \frac{1}{|G_t|} \sum_{g \in G_t} s_g^{(b)}

This yields one null distribution per term. The important point is that the
null is defined on the phenotype side, not by randomizing gene membership.

Empirical p-value
~~~~~~~~~~~~~~~~~

The exported empirical p-value is sign-aware:

.. math::

   p_t =
   \begin{cases}
   \frac{1 + \sum_{b=1}^{B} I(S_t^{(b)} \ge S_t)}{B + 1}, & S_t \ge 0 \\
   \frac{1 + \sum_{b=1}^{B} I(S_t^{(b)} \le S_t)}{B + 1}, & S_t < 0
   \end{cases}

The table also reports:

- ``category_score``
- ``null_mean``
- ``null_sd``
- ``z_score``
- ``p_value``
- ``fdr``

where ``fdr`` is the Benjamini-Hochberg correction across terms in the current
run.

Interpretation
~~~~~~~~~~~~~~

This method is often the most natural fit for imaging transcriptomics because
it carries the phenotype null all the way through to the pathway level. That
is why it is the package default.

GSEA
----

What is ranked
~~~~~~~~~~~~~~

GSEA uses the full ranked gene list rather than a hard threshold.

- in ``corr``, genes are ranked by the gene-level association ``score``
- in ``pls``, genes are ranked component by component using the aligned
  component ``zscore``

Observed enrichment score
~~~~~~~~~~~~~~~~~~~~~~~~~

The observed enrichment score (``es``) is computed with ``gseapy`` on a
deterministic preranked table. When exact tied gene scores occur, the toolbox
adds extremely small stable offsets within the tied groups so that the ranking
order is reproducible and duplicate-score warnings are avoided.

External-null calibration
~~~~~~~~~~~~~~~~~~~~~~~~~

The reported significance is not taken from GSEApy's internal permutations.
Instead, the toolbox evaluates enrichment against null gene rankings derived
from the same imaging permutations used by the main workflow. That means the
reported:

- ``nes``
- ``p_val``
- ``fdr``

are all calibrated to the package's own external nulls.

Normalized enrichment score
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Observed enrichment scores are normalized with the same-sign rule used by
classic GSEA:

.. math::

   NES_t =
   \begin{cases}
   ES_t / \mathrm{mean}(ES_{t,\mathrm{null}}^{+}), & ES_t \ge 0 \\
   -ES_t / \mathrm{mean}(ES_{t,\mathrm{null}}^{-}), & ES_t < 0
   \end{cases}

where :math:`ES_t` is the observed enrichment score for term :math:`t`, and
:math:`ES_{t,\mathrm{null}}^{+}` and :math:`ES_{t,\mathrm{null}}^{-}` are the
same-sign null enrichment scores for that term.

Nominal p-value
~~~~~~~~~~~~~~~

The nominal ``p_val`` is computed from the external enrichment-score nulls:

.. math::

   p_t =
   \begin{cases}
   \frac{1 + \sum_{b=1}^{B} I(ES_t^{(b)} \ge ES_t)}{B + 1}, & ES_t \ge 0 \\
   \frac{1 + \sum_{b=1}^{B} I(ES_t^{(b)} \le ES_t)}{B + 1}, & ES_t < 0
   \end{cases}

GSEA-style FDR
~~~~~~~~~~~~~~

The reported GSEA ``fdr`` is not a Benjamini-Hochberg correction on pathway
p-values. It is a GSEA-style q-value computed from the observed NES values and
the pooled null NES values.

This is why the GSEA ``fdr`` column and the ORA or ensemble ``fdr`` columns
should not be interpreted as the same quantity.

Interpretation
~~~~~~~~~~~~~~

GSEA is most useful when the full gene ordering matters and you do not want to
choose a hard selection threshold. It is a ranking-based enrichment test, not
a direct phenotype-ensemble test.

ORA
---

How genes are selected
~~~~~~~~~~~~~~~~~~~~~~

ORA starts from a thresholded gene list rather than from the full ranking.

The current rule is:

1. keep genes with raw gene-level ``p`` below the chosen threshold
2. split them into ``up`` and ``down`` according to the sign of the gene score
3. run enrichment separately for the positive and negative sets

Important: ORA currently thresholds on the raw gene ``p`` value, not on
gene-level ``fdr`` or ``maxT``.

Hypergeometric test
~~~~~~~~~~~~~~~~~~~

For one pathway, ORA defines a standard contingency table with:

- selected genes in the pathway
- selected genes outside the pathway
- unselected genes in the pathway
- unselected genes outside the pathway

The reported ``p_value`` is the upper-tail hypergeometric probability of
observing at least the observed overlap.

Odds ratio and confidence interval
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ORA also reports effect-size terms:

- ``enrichment_ratio``
- ``odds_ratio``
- ``odds_ratio_ci_low``
- ``odds_ratio_ci_high``

The odds ratio is computed from the contingency table, and the 95 percent
confidence interval is derived on the log-odds scale. When a contingency-table
cell is zero, the implementation applies a Haldane-Anscombe ``+0.5`` correction
before computing the interval.

ORA FDR
~~~~~~~

Unlike GSEA, ORA ``fdr`` is a standard Benjamini-Hochberg correction applied to
the ORA pathway p-values within the current ``up`` or ``down`` analysis.

Interpretation
~~~~~~~~~~~~~~

ORA is the easiest enrichment family to explain because it works on a discrete
hit list. It is most useful when you want overlap counts and direction-specific
tables, but it is also the method that depends most strongly on the chosen
threshold and on the chosen gene universe.

Choosing between methods
------------------------

Use ``ensemble`` when:

- you want the enrichment null to inherit the phenotype null
- spatial autocorrelation is the main inferential concern
- you want the most direct category-score interpretation

Use ``gsea`` when:

- you want a threshold-free ranking-based test
- you care about the ordering of all genes, not just the strongest hits
- you want classic GSEA quantities such as ``es`` and ``nes``

Use ``ora`` when:

- you want a discrete hit-list view
- you want overlap counts and odds ratios
- you want separate positive and negative pathway tables
