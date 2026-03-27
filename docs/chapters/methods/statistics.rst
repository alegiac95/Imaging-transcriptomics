=======================
Statistics and models
=======================

This page summarizes the current statistical workflows implemented by the
toolbox.

Correlation workflow
--------------------

The correlation analysis compares one regional imaging vector with each gene in
the selected atlas expression matrix.

The current implementation:

- standardizes the imaging vector
- rank-transforms the atlas expression matrix and the imaging vector
- computes a Spearman-like association through vectorized matrix operations
- evaluates significance against spatial or grouped random null maps

The gene table written by ``run_corr()`` contains:

- ``score``: the observed gene-wise correlation
- ``p``: nominal permutation p-value from the gene-specific null
- ``fdr``: Benjamini-Hochberg correction across all genes
- ``maxT``: maxT family-wise correction from the permutation maxima

PLS workflow
------------

The PLS workflow models the relationship between one imaging vector and the
full atlas expression matrix using a local PLS backend.

The current implementation:

- standardizes the imaging vector
- fits PLS on the selected atlas rows
- chooses components either from ``n_components`` or a cumulative variance
  target ``var``
- uses permuted imaging maps to evaluate component-level significance
- derives gene-level statistics for each retained component

Each component table contains:

- ``weight``: gene weight on that component
- ``zscore``: weight divided by the bootstrap standard deviation
- ``p``: two-sided z-based gene p-value
- ``fdr``: Benjamini-Hochberg correction across genes in that component
- ``maxT``: maxT family-wise correction from the component-wise maximum
  absolute null weight per permutation

Gene-level correction
---------------------

The toolbox currently uses two complementary correction schemes at the gene
level:

- ``fdr`` controls the false discovery rate across genes with
  Benjamini-Hochberg
- ``maxT`` controls family-wise error using a permutation-based maxT
  statistic

These columns answer different questions. It is therefore normal for a gene to
have a small nominal ``p`` but a much larger ``maxT``.

Permutation resolution
----------------------

Permutation-based p-values are discrete. For ``B`` permutations, the smallest
possible nominal p-value is:

.. math::

   \frac{1}{B + 1}

With about ``15,677`` genes, gene-wise BH-corrected values can remain high
unless the permutation count is very large. This is a limitation of p-value
resolution rather than a sorting bug.

Gene PCA workflow
-----------------

``run_gene_pca()`` is not an imaging association test. It is a descriptive
workflow that:

- filters the atlas expression matrix to the requested genes
- standardizes each gene across regions
- runs PCA on the resulting ``regions x genes`` matrix
- returns regional component scores, gene loadings, and explained variance

This is useful for turning a gene list into one or more regional expression
patterns that can then be compared with imaging maps outside the core
association workflows.

Interpretation
--------------

- for correlation, the sign of ``score`` shows whether a gene is positively or
  negatively associated with the imaging map
- for PLS, component signs are arbitrary before alignment, so interpretation
  should use the aligned component outputs written by the package
- for gene PCA, component signs are also arbitrary; the regional pattern and
  gene loadings should be interpreted together
