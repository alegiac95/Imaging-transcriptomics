==========
Public API
==========

This page is the compact reference for the supported Python entry points
re-exported from ``imaging_transcriptomics``.

Stable workflow functions
-------------------------

- ``run_corr(data, *, atlas="dk", hemisphere="left", regions="all", source_space=None, input_rh=None, n_permutations=1000, null_method="auto", output_dir=None, run_gsea=False, gene_set="lake", ora_p_threshold=None, seed=1234, n_jobs=1)``
- ``run_pls(data, *, atlas="dk", hemisphere="left", regions="all", source_space=None, input_rh=None, n_components=None, var=None, n_permutations=1000, null_method="auto", output_dir=None, run_gsea=False, gene_set="lake", ora_p_threshold=None, seed=1234, n_jobs=1)``
- ``run_gene_pca(genes, *, atlas="dk", hemisphere="left", regions="all", n_components=3, output_dir=None)``
- ``run_analysis(data, config, *, input_rh=None)``

Configuration and selection helpers
-----------------------------------

- ``build_run_config()``
- ``select_atlas_data()``
- ``extract_scan_data()``
- ``list_atlases()``
- ``describe_atlas()``
- ``get_atlas()``
- ``atlas_table()``

Return types
------------

- ``CorrelationResult``: metadata, regional values, ranked gene table, and
  optional GSEA or ORA outputs
- ``PLSResult``: metadata, regional values, retained components, and
  cumulative variance
- ``PLSComponentResult``: one retained PLS component with its gene table and
  optional enrichment outputs
- ``GenePCAResult``: PCA scores, gene loadings, variance table, and matched
  or missing genes

Quick reference
---------------

``run_corr()``
   Use for gene-wise association between one imaging map and one atlas
   expression matrix.

``run_pls()``
   Use for latent-variable analysis between one imaging map and the atlas
   expression matrix.

``run_gene_pca()``
   Use when you already have a gene list and want a regional expression pattern
   derived from those genes only.

``run_analysis()``
   Use when your code constructs a validated ``RunConfig`` separately and
   wants one dispatcher for ``corr`` or ``pls``.

Notes on stability
------------------

These functions and result objects are the intended supported surface. Internal
modules such as ``corr.py``, ``genes.py``, ``nulls.py``, and plotting helpers
are implementation details and may change more freely.
