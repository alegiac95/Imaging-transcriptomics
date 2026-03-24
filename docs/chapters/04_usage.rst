.. _Usage:

============
Script usage
============

Once installed, the v2 CLI exposes atlas-aware workflows:

.. code:: bash

    imagingtranscriptomics atlases --packaged-only
    imagingtranscriptomics corr --input /path/to/map.nii.gz --atlas dk --output /path/to/out
    imagingtranscriptomics pls --input /path/to/map.tsv --atlas schaefer-100 --ncomp 2 --output /path/to/out

The shared analysis options are:

- ``--input`` / ``-i``: input regional vector, volumetric NIfTI, or surface file pair.
- ``--input-rh``: right-hemisphere surface file when using surface inputs.
- ``--atlas`` / ``-a``: atlas preset, for example ``dk`` or ``schaefer-100``.
- ``--hemisphere``: ``left`` or ``both``. Packaged expression data preserve the abagen left-only workflow and also expose mirrored right-hemisphere values when available.
- ``--regions`` / ``-r``: ``all``, ``cort``, or ``cort+sub``.
- ``--space``: source space for non-native inputs, for example ``MNI152``, ``fsaverage``, ``fsLR``, or ``CIVET``.
- ``--permutations`` / ``-p``: number of permutations or spatial null samples.
- ``--null-method``: ``auto``, ``vasa``, ``alexander_bloch``, ``moran``, or ``random``.
- ``--seed``: random seed for reproducible permutations and null-model generation.
- ``--geneset`` and ``--no-gsea``: enable or disable optional GSEA output.

``corr`` runs Spearman correlation against all genes. ``pls`` runs the local SIMPLS-based PLS backend and requires either ``--ncomp`` or ``--var``.

Outputs are written as lightweight text and image files:

- ``README.txt``: human-readable run summary.
- ``metadata.json``: machine-readable metadata, including atlas and null-model settings.
- ``regional_values.tsv``: parcellated input values aligned to atlas labels.
- ``corr_genes.tsv`` or ``pls_component_<n>.tsv`` / ``pls_summary.tsv``: analysis tables.
- ``plots/*.png``: overview plots instead of a PDF report.


.. _library:

=======================
Usage as python library
=======================

The v2 API is function-first:

.. code:: python

    import imaging_transcriptomics as imt

    corr_result = imt.run_corr(
        "/path/to/map.nii.gz",
        atlas="dk",
        hemisphere="left",
        null_method="auto",
        output_dir="/path/to/out",
    )

    pls_result = imt.run_pls(
        "/path/to/map.tsv",
        atlas="schaefer-100",
        hemisphere="both",
        n_components=2,
        output_dir="/path/to/out",
    )

You can also build an explicit configuration object:

.. code:: python

    config = imt.build_run_config(
        "corr",
        atlas="dk",
        hemisphere="left",
        regions="all",
        n_permutations=1000,
        null_method="vasa",
        output_dir="/path/to/out",
    )
    result = imt.run_analysis("/path/to/map.nii.gz", config)

Atlas selection and scan extraction are exposed directly:

.. code:: python

    atlas = imt.get_atlas("dk")
    selection = imt.select_atlas_data(atlas="dk", hemisphere="both", regions="all")
    extracted = imt.extract_scan_data("/path/to/map.nii.gz", atlas="dk", hemisphere="left")
