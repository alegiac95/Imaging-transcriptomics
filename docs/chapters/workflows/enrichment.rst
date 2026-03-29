=====================
Enrichment workflows
=====================

This page unifies user guidance for GSEA and ORA across the correlation and
PLS workflows.

Overview
--------

The toolbox offers two enrichment styles:

- GSEA, which uses the full ranked gene list
- ORA, which uses a thresholded subset of genes

Both are available:

- after ``corr``
- after each retained PLS component

How to enable them
------------------

GSEA
~~~~

GSEA is enabled with ``--gsea`` on the CLI or ``run_gsea=True`` in the API.
It can be used with:

- packaged genesets such as ``lake`` and ``pooled``
- local ``.gmt`` files
- remote Enrichr libraries resolved through ``gseapy``

ORA
~~~

ORA is enabled by giving a threshold with ``--ora-p-threshold``. The threshold
is applied to the raw gene-level p-values produced by the main workflow.

For example:

.. code-block:: bash

   imt corr \
     --input /abs/path/map.nii.gz \
     --space MNI152 \
     --atlas dk \
     --geneset lake \
     --ora-p-threshold 0.01 \
     --output /abs/path/out_dir

and for both ORA and GSEA:

.. code-block:: bash

   imt corr \
     --input /abs/path/map.nii.gz \
     --space MNI152 \
     --atlas dk \
     --geneset GO_Biological_Process_2025 \
     --geneset-organism Human \
     --ora-p-threshold 0.01 \
     --gsea \
     --output /abs/path/out_dir

Choosing a geneset resource
---------------------------

The ``--geneset`` argument accepts:

- packaged entries such as ``lake`` and ``pooled``
- a local GMT file
- a remote Enrichr/GSEApy library name

Useful commands:

.. code-block:: bash

   imt genesets --packaged-only
   imt genesets --organism Human

Interpreting GSEA outputs
-------------------------

The main GSEA columns are:

``Term``
   Pathway or geneset name.

``es``
   Observed enrichment score.

``nes``
   Same-sign normalized enrichment score calibrated to the external nulls used
   by the package.

``p_val``
   Sign-aware nominal enrichment p-value from the external null enrichment
   scores.

``fdr``
   GSEA-style q-value derived from observed and null normalized enrichment
   scores.

The most important practical point is that GSEA ``fdr`` is not a BH correction
on pathway p-values. It is a GSEA-style quantity and should be interpreted as
such.

Interpreting ORA outputs
------------------------

ORA writes separate ``up`` and ``down`` tables.

Important columns:

``overlap_size``
   Number of selected genes that fall in the pathway.

``selected_size``
   Number of genes that passed the raw gene-p threshold in that direction.

``enrichment_ratio``
   Observed overlap divided by expected overlap.

``odds_ratio``
   Contingency-table odds ratio.

``odds_ratio_ci_low`` / ``odds_ratio_ci_high``
   95 percent confidence interval for the odds ratio.

``p_value``
   Hypergeometric enrichment p-value.

``fdr``
   Benjamini-Hochberg correction across ORA terms within that direction.

So ORA ``fdr`` and GSEA ``fdr`` are not the same kind of quantity.

Plots
-----

GSEA plots
~~~~~~~~~~

The package writes dotplots of the top GSEA terms. These are useful for a fast
overview of the strongest enriched terms.

ORA plots
~~~~~~~~~

The package writes a two-row heatmap:

- one row for ``up``
- one row for ``down``

Cells are colored by enrichment significance and annotated with the odds ratio
plus significance stars. The heatmap is intentionally capped to a limited
number of terms so very large libraries such as GO biological process remain
readable.

Common reasons for weak or flat results
---------------------------------------

Flat GSEA ``fdr``
   Often means the library is large, the null is broad, or the enrichment is
   modest rather than that the code failed.

Nearly empty ORA
   Usually means the raw gene-p threshold is too strict or the selected
   pathway library has little overlap with the selected genes.

Huge GSEA runtime
   Usually comes from large remote libraries and repeated external-null
   enrichment calculations, not from the correlation step itself.

Large ORA tables but unreadable plots
   Usually happens with very large libraries. The plot now caps the number of
   displayed terms, but the TSV will still contain the full result table.
