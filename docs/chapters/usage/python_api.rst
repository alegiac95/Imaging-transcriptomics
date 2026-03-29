==========
Python API
==========

The v2 Python surface is function-first. The stable entry points live in the
top-level package and return typed dataclass-like result objects.

High-level functions
--------------------

- ``run_corr()``
- ``run_pls()``
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
- ``GenePCAResult``
- ``GEDARResult``

Result structure
----------------

- ``CorrelationResult`` stores ``metadata``, ``regional_values``,
  ``gene_table``, and optional ``gsea_table`` / ``ora_tables``
- ``PLSResult`` stores ``metadata``, ``regional_values``, one
  ``PLSComponentResult`` per retained component, and
  ``cumulative_variance``
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
