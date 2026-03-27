====================
Methods at a glance
====================

The methods documentation is split into four complementary pages.

Method families
---------------

- gene-level association workflows
- spatial null models
- pathway enrichment analyses
- gene-list PCA

How they fit together
---------------------

``statistics``
   Explains the correlation, PLS, and gene-PCA workflows.

``null_models``
   Explains how imaging permutations are generated before downstream testing.

``gene_sets``
   Explains how GSEA and ORA are layered on top of the gene-level outputs.

In practice, most users move through these pages in that order:

1. understand the main workflow
2. understand the null model it depends on
3. understand how enrichment is computed from the resulting gene table
