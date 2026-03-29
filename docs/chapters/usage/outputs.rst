=======
Outputs
=======

Each workflow writes a small self-contained output bundle designed for both
manual inspection and downstream scripting.

What every run writes
---------------------

Every persisted workflow writes:

- ``README.txt``
- ``metadata.json``
- one or more TSV tables
- one or more PNG plots when plotting dependencies are available

For ``corr`` and ``pls``, the bundle also includes ``regional_values.tsv`` so
the extracted atlas-aligned imaging vector is always inspectable.

Core files
----------

``README.txt``
   Human-readable summary of the run, the main settings, and the generated
   files.

``metadata.json``
   Machine-readable record of atlas choice, hemisphere, region subset, null
   method, permutation count, enrichment settings, and workflow-specific
   bookkeeping.

Workflow-specific outputs
-------------------------

``corr``
   Writes ``corr_genes.tsv`` plus optional ``gsea_corr_results.tsv`` and
   ``ora_corr_up.tsv`` / ``ora_corr_down.tsv``.

``pls``
   Writes ``pls_summary.tsv`` plus one ``pls_component_<n>.tsv`` per retained
   component and optional enrichment tables for each component.

``gene-pca``
   Writes ``gene_pca_scores.tsv``, ``gene_pca_loadings.tsv``,
   ``gene_pca_variance.tsv``, ``matched_genes.txt``,
   ``brain_filtered_genes.txt``, and ``missing_genes.txt``.

``gedar``
   Writes ``gedar_scores.tsv``, ``gedar_genes.tsv``, ``gedar_excluded.tsv``,
   ``matched_genes.txt``, and ``missing_genes.txt``.

How to read the main tables
---------------------------

``corr_genes.tsv``
   Ranked gene-level correlation results with association score, nominal
   p-value, BH FDR, and maxT family-wise correction.

``pls_summary.tsv``
   Component-level explained variance and permutation p-values.

``pls_component_<n>.tsv``
   Gene weights and gene-level statistics for one aligned PLS component.

``gene_pca_scores.tsv``
   Regional PCA component scores for a selected gene list.

``gene_pca_loadings.tsv``
   Gene contributions to each PCA component.

``gedar_scores.tsv``
   Regional GEDAR weighted-expression scores, optionally with separate up/down
   score columns in split mode.

``gedar_genes.tsv``
   Matched genes, retained weights, rank columns, and selection flags used in
   the GEDAR calculation.

Output bundles and interpretation
---------------------------------

The design goal is:

- a small number of plain-text files
- enough metadata to rerun the analysis
- enough tables to inspect the results without opening the code

The files are intended to be complementary:

- ``README.txt`` for a quick human summary
- ``metadata.json`` for scripting and provenance
- TSV tables for detailed analysis
- PNG figures for fast visual checks

See also
--------

For column-level definitions, use the dedicated reference page:
``reference/file_formats``.
