==========
Python API
==========

The v2 Python surface is function-first. The stable entry points live in the
top-level package and return typed dataclass-like result objects.

High-level functions
--------------------

- ``run_corr()``
- ``run_pls()``
- ``run_gene()``
- ``run_gene_pca()``
- ``run_gedar()``
- ``run_analysis()``

Typical usage
-------------

Correlation:

.. code-block:: python

   import imaging_transcriptomics as imt

   result = imt.run_corr(
       "/absolute/path/to/map.nii.gz",
       atlas="dk",
       hemisphere="left",
       source_space="MNI152",
       n_permutations=20000,
       null_method="auto",
       run_gsea=True,
       gene_set="lake",
       output_dir="/absolute/path/to/out_dir",
   )

   result.gene_table.head()

PLS:

.. code-block:: python

   import imaging_transcriptomics as imt

   result = imt.run_pls(
       "/absolute/path/to/map.nii.gz",
       atlas="dk",
       hemisphere="left",
       source_space="MNI152",
       n_components=2,
       n_permutations=20000,
       n_jobs=8,
   )

   result.components[0].gene_table.head()

Single-gene query:

.. code-block:: python

   import imaging_transcriptomics as imt

   result = imt.run_gene(
       "RELN",
       atlas="dk",
       hemisphere="left",
       top_n=25,
   )

   result.regional_values.head()
   result.coexpressed_genes.head()

Gene PCA:

.. code-block:: python

   import imaging_transcriptomics as imt

   result = imt.run_gene_pca(
       ["RELN", "GAD1", "SLC1A2", "SV2A"],
       atlas="dk",
       hemisphere="left",
       n_components=3,
   )

   result.variance_table

GEDAR:

.. code-block:: python

   import imaging_transcriptomics as imt

   result = imt.run_gedar(
       "/absolute/path/to/twas.tsv",
       atlas="dk",
       gene_column="gene_name",
       weight_column="z_mean",
       rank_column="pvalue",
       top_percent=5,
       direction="combined",
   )

   result.regional_scores.head()

Supporting helpers
------------------

- ``build_run_config()``
- ``select_atlas_data()``
- ``extract_scan_data()``
- atlas lookup helpers such as ``get_atlas()`` and ``atlas_table()``

Configuration helper
--------------------

``build_run_config()`` validates shared options for correlation and PLS runs.
This is most useful when you want one code path that can dispatch to
``run_analysis()``:

.. code-block:: python

   import imaging_transcriptomics as imt

   config = imt.build_run_config(
       "corr",
       atlas="dk",
       hemisphere="left",
       source_space="MNI152",
       n_permutations=5000,
       run_gsea=True,
   )
   result = imt.run_analysis("/absolute/path/to/map.nii.gz", config)

Result objects
--------------

- ``CorrelationResult``
- ``PLSResult``
- ``PLSComponentResult``
- ``GeneQueryResult``
- ``GenePCAResult``
- ``GEDARResult``

Result structure
----------------

- ``CorrelationResult`` stores ``metadata``, ``regional_values``,
  ``gene_table``, and optional ``gsea_table`` / ``ora_tables``
- ``PLSResult`` stores ``metadata``, ``regional_values``, one
  ``PLSComponentResult`` per retained component, and
  ``cumulative_variance``
- ``GeneQueryResult`` stores one gene's ``regional_values``, the full
  co-expression ``gene_table``, and the filtered ``coexpressed_genes`` table
- ``GenePCAResult`` stores the filtered PCA outputs:
  ``regional_scores``, ``gene_loadings``, ``variance_table``,
  ``matched_genes``, ``brain_filtered_genes``, and ``missing_genes``
- ``GEDARResult`` stores ``regional_scores``, ``gene_table``,
  ``excluded_table``, ``matched_genes``, and ``missing_genes``

Stability guidance
------------------

The functions re-exported from ``imaging_transcriptomics`` are the supported
surface for normal users. Modules such as ``corr.py``, ``genes.py``,
``nulls.py``, and plotting internals are implementation details and may change
more frequently.

Stable API reference
--------------------

Stable workflow functions
~~~~~~~~~~~~~~~~~~~~~~~~~

- ``run_corr(data, *, atlas="dk", hemisphere="left", regions="default", source_space=None, input_rh=None, n_permutations=1000, null_method="auto", output_dir=None, run_gsea=False, gene_set="lake", ora_p_threshold=None, seed=1234, n_jobs=1)``
- ``run_pls(data, *, atlas="dk", hemisphere="left", regions="default", source_space=None, input_rh=None, n_components=None, var=None, n_permutations=1000, null_method="auto", output_dir=None, run_gsea=False, gene_set="lake", ora_p_threshold=None, seed=1234, n_jobs=1)``
- ``run_gene(gene, *, atlas="dk", hemisphere="left", regions="default", zscore_expression=True, top_n=25, fdr_threshold=0.05, output_dir=None)``
- ``run_gene_pca(genes, *, atlas="dk", hemisphere="left", regions="default", n_components=3, output_dir=None)``
- ``run_gedar(weights, *, atlas="dk", hemisphere="both", regions="default", gene_column="gene", weight_column="weight", rank_column=None, rank_mode="ascending", top_percent=None, top_n=None, p_threshold=None, direction="combined", normalize_expression="zscore", normalize_weights="none", output_dir=None)``
- ``run_analysis(data, config, *, input_rh=None)``

Configuration and selection helpers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``build_run_config()``
- ``select_atlas_data()``
- ``extract_scan_data()``
- ``list_atlases()``
- ``describe_atlas()``
- ``get_atlas()``
- ``atlas_table()``

Return types
~~~~~~~~~~~~

- ``CorrelationResult``: metadata, regional values, ranked gene table, and
  optional GSEA or ORA outputs
- ``PLSResult``: metadata, regional values, retained components, and
  cumulative variance
- ``PLSComponentResult``: one retained PLS component with its gene table and
  optional enrichment outputs
- ``GeneQueryResult``: one gene's regional expression vector, the full
  co-expression table, and the selected top co-expressed genes
- ``GenePCAResult``: PCA scores, gene loadings, variance table, and matched
  or missing genes
- ``GEDARResult``: regional weighted-expression scores, matched gene table,
  excluded rows, and gene bookkeeping

Quick reference
~~~~~~~~~~~~~~~

``run_corr()``
   Use for gene-wise association between one imaging map and one atlas
   expression matrix.

``run_pls()``
   Use for latent-variable analysis between one imaging map and the atlas
   expression matrix.

``run_gene()``
   Use when you want one atlas-aligned gene expression profile together with
   the top significantly positively co-expressed genes in that atlas.

``run_gene_pca()``
   Use when you already have a gene list and want a regional expression pattern
   derived from those genes only.

``run_gedar()``
   Use when you have a weighted gene table, such as a TWAS-like result, and
   want a regional weighted-expression score on one atlas.

``run_analysis()``
   Use when your code constructs a validated ``RunConfig`` separately and
   wants one dispatcher for ``corr`` or ``pls``.
