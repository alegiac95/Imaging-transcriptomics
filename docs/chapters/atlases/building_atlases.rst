================
Building atlases
================

Atlas assets are generated from external atlas definitions plus ``abagen``.
The packaged presets in this branch keep the ``abagen``-derived expression
matrices and store them in a format optimized for loading rather than manual
inspection.

What a build needs
------------------

- a volumetric atlas image for ``MNI152`` workflows when available
- a label table with at least region ``id``, ``label``, ``hemisphere``, and
  ``structure``
- surface parcellation files when the atlas supports surface workflows
- optional surface geometry files for surface-based null models

Build policy
------------

The current atlas builder keeps these high-level rules:

- expression is generated with ``abagen``
- left-right mirroring is enabled so the packaged matrix can support
  ``hemisphere="both"``
- atlas labels and expression are versioned together
- per-atlas ``README.txt`` and provenance metadata are written with the build

Output layout
-------------

Each packaged atlas folder contains a subset of these files:

- ``atlas-<name>_labels.csv``
- ``atlas-<name>_gene_expression_data.npz``
- surface parcellation files such as ``.annot`` or ``.label.gii``
- packaged volume images for ``1mm`` and ``2mm`` extraction where available
- ``README.txt`` describing how the atlas was built
- ``provenance.json`` with machine-readable build details

Shared gene labels are stored separately in the atlas shared-data directory so
they can be reused across atlases.

Python entry point
------------------

The public builder entry point is ``build_expression_assets()``, re-exported
from the package root:

.. code-block:: python

   import imaging_transcriptomics as imt

   imt.build_expression_assets(
       "schaefer-200",
       "/absolute/path/to/output_dir",
       atlas_image="/absolute/path/to/atlas.nii.gz",
       atlas_info="/absolute/path/to/labels.csv",
       lr_mirror="leftright",
       n_proc=4,
   )

Validation checks worth keeping
-------------------------------

Before packaging a new atlas build, validate:

- the number of rows in the expression matrix matches the label table
- atlas IDs match between volume or surface files and labels
- hemisphere labels are correct
- cortex versus subcortex structure labels are correct
- left-only and both-hemisphere selection return the expected row counts
- the atlas can be used by ``imt atlases`` and by a small smoke-test run
