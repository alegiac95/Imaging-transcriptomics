=============
CLI reference
=============

This page summarizes the supported commands and the options users most often
need to reach for.

Commands to document
--------------------

- ``imt atlases``
- ``imt corr``
- ``imt pls``
- ``imt gene-pca``

Synopsis
--------

.. code-block:: text

   imt atlases [--packaged-only]
   imt corr --input PATH [shared options]
   imt pls --input PATH (--ncomp N | --var FRACTION) [shared options]
   imt gene-pca --genes VALUE [gene-pca options]

Shared analysis options
-----------------------

``corr`` and ``pls`` share these main options:

- ``--input``: vector, NIfTI map, or supported surface input
- ``--input-rh``: matching right-hemisphere surface file
- ``--output``: output directory
- ``--atlas``: atlas preset ID
- ``--hemisphere``: ``left`` or ``both``
- ``--regions``: ``all``, ``cort``, or ``cort+sub``
- ``--space``: source space, such as ``MNI152``, ``fsaverage``, ``fsLR``, or ``CIVET``
- ``--permutations``: number of null samples
- ``--seed``: random seed
- ``--jobs``: worker count for parallel work such as PLS permutations
- ``--null-method``: ``auto``, ``vasa``, ``alexander_bloch``, ``moran``, or ``random``
- ``--geneset``: packaged gene set name or local ``.gmt`` file
- ``--ora-p-threshold``: raw gene p-value threshold for ORA
- ``--gsea`` / ``--no-gsea``: explicit GSEA control

PLS-specific options
--------------------

- ``--ncomp``: keep a fixed number of PLS components
- ``--var``: keep enough components to reach a cumulative variance target

Gene-PCA-specific options
-------------------------

- ``--genes``: text file, TSV, CSV, or comma-separated list of gene symbols
- ``--atlas`` / ``--hemisphere`` / ``--regions``: same atlas-selection
  semantics as the other workflows
- ``--ncomp``: maximum number of PCA components to retain

Behavior notes
--------------

- when ``--ora-p-threshold`` is given and neither ``--gsea`` nor ``--no-gsea``
  is provided, the CLI runs ORA only
- ``imagingtranscriptomics`` is still available as a long-form alias for
  ``imt``
