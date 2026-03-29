=========================
GSEA and ORA methodology
=========================

The toolbox supports two pathway-level analyses built on top of the gene-level
outputs:

- preranked GSEA
- over-representation analysis (ORA)

They answer related but different questions, so it is often useful to run both
on the same workflow.

GSEA
----

What is ranked
~~~~~~~~~~~~~~

GSEA uses the full ranked gene list rather than a hard threshold.

- in ``corr``, genes are ranked by the correlation ``score``
- in ``pls``, genes are ranked component by component using the aligned
  component ``zscore``

Observed enrichment score
~~~~~~~~~~~~~~~~~~~~~~~~~

The observed enrichment score (``es``) is computed with ``gseapy`` on a
deterministic preranked table. When exact tied gene scores occur, the toolbox
adds extremely small stable offsets within the tied groups so that the ranking
order is reproducible and the noisy duplicate-score warning is avoided.

External-null calibration
~~~~~~~~~~~~~~~~~~~~~~~~~

The toolbox does not use GSEApy's internal permutation engine for the final
reported significance. Instead it recomputes enrichment on external null gene
rankings derived from the same imaging permutations used by the main workflow.

This matters because the reported:

- ``nes``
- ``p_val``
- ``fdr``

are all calibrated to the package's own external nulls.

Normalized enrichment score
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Observed enrichment scores are normalized with the same-sign rule used by
classic GSEA and GSEApy:

.. math::

   NES_t =
   \begin{cases}
   ES_t / \mathrm{mean}(ES_{t,\mathrm{null}}^{+}), & ES_t \ge 0 \\
   -ES_t / \mathrm{mean}(ES_{t,\mathrm{null}}^{-}), & ES_t < 0
   \end{cases}

where:

- :math:`ES_t` is the observed enrichment score for term :math:`t`
- :math:`ES_{t,\mathrm{null}}^{+}` are positive null enrichment scores for the same term
- :math:`ES_{t,\mathrm{null}}^{-}` are negative null enrichment scores for the same term

The sign-specific mean prevents positive and negative terms from being
normalized against incompatible null tails.

Nominal GSEA p-value
~~~~~~~~~~~~~~~~~~~~

The nominal ``p_val`` is computed from the external enrichment-score nulls:

.. math::

   p_t =
   \begin{cases}
   \frac{1 + \sum_{b=1}^{B} I(ES_t^{(b)} \ge ES_t)}{B + 1}, & ES_t \ge 0 \\
   \frac{1 + \sum_{b=1}^{B} I(ES_t^{(b)} \le ES_t)}{B + 1}, & ES_t < 0
   \end{cases}

So the tail direction follows the sign of the observed enrichment score.

GSEA-style FDR
~~~~~~~~~~~~~~

The reported GSEA ``fdr`` is not a Benjamini-Hochberg correction on pathway
p-values. It is a GSEA-style q-value computed from the observed NES values and
the pooled null NES values.

Conceptually, for one observed term:

.. math::

   FDR(\mathrm{NES}) =
   \frac{P(\mathrm{null\ NES} \ge \mathrm{NES})}{P(\mathrm{observed\ NES} \ge \mathrm{NES})}

for positive scores, with the symmetric left-tail version for negative scores.

This is why the GSEA ``fdr`` column and the ORA ``fdr`` column should not be
interpreted as the same quantity.

ORA
---

How genes are selected
~~~~~~~~~~~~~~~~~~~~~~

ORA starts from a thresholded gene list rather than from the full ranking.

The current rule is:

1. keep genes with raw gene-level ``p`` below the chosen threshold
2. split them into ``up`` and ``down`` according to the sign of the gene score
3. run enrichment separately for the positive and negative sets

Important: ORA currently thresholds on the raw gene p-value, not on gene-level
``fdr`` or ``maxT``.

Hypergeometric test
~~~~~~~~~~~~~~~~~~~

For one pathway, ORA defines a standard contingency table with:

- selected genes in the pathway
- selected genes outside the pathway
- unselected genes in the pathway
- unselected genes outside the pathway

The reported ``p_value`` is the upper-tail hypergeometric probability of
observing at least the observed overlap:

.. math::

   p = P(X \ge k)

where:

- ``k`` is the observed overlap
- ``X`` follows a hypergeometric distribution with the current universe,
  pathway size, and selected-gene count

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

So:

- GSEA ``fdr`` = GSEA-style NES-based q-value from external nulls
- ORA ``fdr`` = BH-adjusted hypergeometric p-value

Choosing between GSEA and ORA
-----------------------------

Use GSEA when:

- you want a threshold-free ranking-based test
- you do not want to choose a hard gene cutoff
- you care about coordinated weak-to-moderate shifts across many genes

Use ORA when:

- you want a discrete hit list
- you want pathway odds ratios and overlap genes
- you want separate positive and negative subsets with a clear selection rule

Use both when:

- you want a threshold-free and a thresholded view of the same result
- you want to compare broad ranked enrichment with more selective hit-based
  enrichment
