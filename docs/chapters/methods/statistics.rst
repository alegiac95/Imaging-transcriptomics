=======================
Statistics and models
=======================

This page summarizes the statistical machinery used by the main analysis
workflows. It is written to match the current implementation in the toolbox,
so the formulas below describe what the package actually computes rather than
an idealized textbook version.

Notation
--------

Throughout this page:

- ``n`` is the number of selected atlas regions
- ``p`` is the number of genes in the selected atlas expression matrix
- ``B`` is the number of imaging permutations
- ``x`` is the observed regional imaging vector with shape ``n x 1``
- ``X`` is the atlas expression matrix with shape ``n x p``
- ``x^(b)`` is the ``b``\ th permuted imaging vector

The package always works on the selected atlas subset, so ``n`` changes with
the chosen atlas, hemisphere, and region scope.

Correlation workflow
--------------------

The correlation workflow compares one regional imaging vector with every gene
in the atlas expression matrix.

Observed statistic
~~~~~~~~~~~~~~~~~~

The current implementation computes a Spearman-like correlation by ranking the
imaging vector and each gene across regions, standardizing those ranks, and
then taking a dot product:

.. math::

   r_g = \frac{\tilde{x}^{\mathsf T} \tilde{g}}{n - 1}

where:

- ``g`` is one gene-expression column from ``X``
- ``tilde{x}`` is the standardized rank vector of the imaging values
- ``tilde{g}`` is the standardized rank vector of one gene

Two implementation details matter:

- ranks are produced with a stable sort, so exact ties are broken by order
  rather than by average rank
- standardization uses the sample standard deviation with ``ddof=1``

This makes the workflow very close to Spearman correlation while keeping the
computation fully vectorized across genes.

Null distribution
~~~~~~~~~~~~~~~~~

The null model is generated on the imaging side, not on the gene-expression
side. The toolbox creates ``B`` permuted imaging vectors and correlates each
of them with the same ranked gene-expression matrix:

.. math::

   r_g^{(b)} = \frac{{\tilde{x}^{(b)}}^{\mathsf T} \tilde{g}}{n - 1}

This produces one null distribution per gene.

Nominal p-value
~~~~~~~~~~~~~~~

The exported gene-level ``p`` column is an empirical sign-aware p-value. The
tail is chosen from the sign of the observed statistic:

.. math::

   p_g =
   \begin{cases}
   \frac{1 + \sum_{b=1}^{B} I(r_g^{(b)} \ge r_g)}{B + 1}, & r_g \ge 0 \\
   \frac{1 + \sum_{b=1}^{B} I(r_g^{(b)} \le r_g)}{B + 1}, & r_g < 0
   \end{cases}

This is not a two-sided absolute-value test. Instead, it asks whether the
gene is unusually positive or unusually negative in the same direction as the
observed score.

False discovery rate
~~~~~~~~~~~~~~~~~~~~

The exported ``fdr`` column is the Benjamini-Hochberg correction applied to
the nominal gene-level p-values across all genes in the current run.

Conceptually:

1. sort the nominal p-values
2. scale them by the number of genes and their rank
3. enforce monotonicity from the smallest p-value upward

This correction controls the expected false discovery rate across the gene
table, but it is limited by permutation resolution. If ``B`` is small relative
to the number of genes, the smallest attainable adjusted value can remain very
large.

maxT family-wise correction
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The exported ``maxT`` column is a maxT-style family-wise error correction based
on the most extreme absolute null statistic in each permutation:

.. math::

   M^{(b)} = \max_g |r_g^{(b)}|

and then:

.. math::

   p_g^{\mathrm{maxT}} =
   \frac{1 + \sum_{b=1}^{B} I(M^{(b)} \ge |r_g|)}{B + 1}

This is stricter than ``fdr`` because it controls the chance of seeing an
effect at least that large anywhere in the full gene table.

PLS workflow
------------

The PLS workflow is a one-response PLS model where the imaging map is the
response and the atlas expression matrix is the predictor matrix.

Model
~~~~~

Let ``X`` be the ``n x p`` atlas expression matrix and ``y`` the imaging
vector. The toolbox uses a local PLS-1 backend derived from the SIMPLS
algorithm.

For each component ``k`` it finds a weight vector ``w_k`` and component score
vector ``t_k`` such that:

.. math::

   t_k = X w_k

and the component captures covariance between the imaging vector and the gene
expression matrix under the usual SIMPLS orthogonality constraints.

Component count
~~~~~~~~~~~~~~~

The number of retained components is determined in one of two ways:

- ``n_components`` keeps exactly that many components
- ``var`` keeps the smallest number of components whose cumulative explained
  response variance reaches the requested threshold

The stored ``variance_explained`` values are the per-component fractions of
response-side variance explained by the fitted latent variables. The stored
``cumulative_variance`` values are the cumulative sums of those fractions.

Component sign alignment
~~~~~~~~~~~~~~~~~~~~~~~~

PLS component signs are arbitrary. After fitting, the toolbox aligns each
component so that the correlation between the component scores and the imaging
vector is positive. This makes component outputs easier to compare across runs.

Component-level p-values
~~~~~~~~~~~~~~~~~~~~~~~~

PLS significance is evaluated by fitting the same number of components to each
permuted imaging vector. For each permutation the workflow recomputes the
cumulative explained variance curve and compares it with the observed one.

For component ``k``:

.. math::

   p_k =
   \frac{1 + \sum_{b=1}^{B} I(R_{k,\mathrm{perm}}^{(b)} \ge R_{k,\mathrm{obs}})}{B + 1}

where ``R_k`` is the cumulative variance explained through component ``k``.

This means the component p-value is attached to the cumulative model up to that
component, not only to the marginal increment of the single component.

PLS gene-level statistics
~~~~~~~~~~~~~~~~~~~~~~~~~

For each retained component, the package also stores the gene weights from the
observed fit and from every permuted fit. Permuted gene weights are:

- reordered into the original gene ranking
- sign-aligned to the observed component before statistics are computed

The exported columns are:

- ``weight``: the aligned observed gene weight
- ``zscore``: the observed weight divided by the permutation standard deviation
- ``p``: empirical sign-aware p-value from the permutation weight distribution
- ``fdr``: Benjamini-Hochberg correction across genes in that component
- ``maxT``: family-wise correction from the maximum absolute null weight in
  each permutation

The nominal p-value uses the same sign-aware rule as the correlation workflow:

.. math::

   p_g =
   \begin{cases}
   \frac{1 + \sum_{b=1}^{B} I(w_g^{(b)} \ge w_g)}{B + 1}, & w_g \ge 0 \\
   \frac{1 + \sum_{b=1}^{B} I(w_g^{(b)} \le w_g)}{B + 1}, & w_g < 0
   \end{cases}

and the component-wise maxT correction uses:

.. math::

   M_k^{(b)} = \max_g |w_{g,k}^{(b)}|

.. math::

   p_{g,k}^{\mathrm{maxT}} =
   \frac{1 + \sum_{b=1}^{B} I(M_k^{(b)} \ge |w_{g,k}|)}{B + 1}

Important interpretation note: the exported ``zscore`` is descriptive and helps
rank genes, but the inferential quantities are the empirical ``p``, ``fdr``,
and ``maxT`` columns.

Permutation resolution
----------------------

All empirical p-values in the package use the standard ``+1`` correction:

.. math::

   p_{\min} = \frac{1}{B + 1}

This has two practical consequences:

- nominal p-values are discrete
- corrected values can remain large when the number of tests is much larger
  than the number of permutations

For gene-level BH correction over ``m`` genes, the best possible adjusted
value is approximately:

.. math::

   \frac{m}{B + 1}

So with roughly ``15{,}677`` genes, tens of thousands of permutations may still
be insufficient for small adjusted gene-level values.

How to interpret the columns
----------------------------

``score`` in correlation
   Direction and strength of the regional association between the imaging map
   and one gene.

``weight`` in PLS
   Direction and contribution of a gene to one aligned PLS component.

``p``
   Gene-specific empirical p-value from the relevant permutation null.

``fdr``
   Benjamini-Hochberg false discovery rate across genes in the current table.

``maxT``
   Family-wise error correction against the strongest absolute null statistic
   observed in each permutation.

These numbers answer different questions, so it is common to see a small
nominal ``p`` together with a much larger ``maxT``.

Gene PCA
--------

``run_gene_pca()`` is not a hypothesis test. It is a descriptive projection of
selected genes into a lower-dimensional regional expression space.

The workflow:

1. filters the requested genes to those present in the atlas expression matrix
2. optionally removes genes outside the packaged brain-gene filter
3. standardizes each retained gene across regions
4. runs PCA on the resulting ``regions x genes`` matrix

The outputs are:

- regional component scores
- per-gene loadings
- explained and cumulative variance

PCA component signs are arbitrary, so regional scores and gene loadings should
always be interpreted together.
