================
Included atlases
================

The package ships six atlas presets with packaged expression assets derived
from ``abagen``. All atlas selection happens through the registry, so the same
IDs work in the CLI and Python API.

Current presets
---------------

.. list-table::
   :header-rows: 1
   :widths: 18 10 10 20 20 22

   * - Atlas ID
     - Left regions
     - Both regions
     - Coverage
     - Supported spaces
     - Notes
   * - ``dk``
     - 41
     - 83
     - cortex + subcortex
     - ``MNI152``, ``fsaverage``
     - Legacy-compatible Desikan-Killiany preset with mirrored right-side expression.
   * - ``schaefer-100``
     - 50
     - 100
     - cortex
     - ``MNI152``, ``fsaverage``
     - Packaged Schaefer 100 cortical atlas.
   * - ``schaefer-200``
     - 100
     - 200
     - cortex
     - ``MNI152``, ``fsaverage``, ``fsLR``
     - Higher-resolution Schaefer preset with packaged local build assets.
   * - ``schaefer-400``
     - 200
     - 400
     - cortex
     - ``MNI152``, ``fsaverage``, ``fsLR``
     - Higher-resolution Schaefer preset with packaged local build assets.
   * - ``destrieux``
     - 74
     - 148
     - cortex
     - ``MNI152``, ``fsaverage``
     - Destrieux cortical atlas built from local atlas assets.
   * - ``glasser-360``
     - 180
     - 360
     - cortex
     - ``fsLR``, ``fsaverage``, ``MNI152``
     - Surface-oriented Glasser preset with packaged surface geometry.

Hemisphere modes
----------------

All packaged atlases support:

- ``hemisphere="left"``
- ``hemisphere="both"``

For ``both``, right-hemisphere expression comes from the ``abagen`` left-right
mirror option used at atlas-build time. The package does not estimate a second
independent AHBA matrix for the right hemisphere.

Coverage and region scope
-------------------------

- ``regions="all"`` keeps all atlas rows selected by the hemisphere mode
- ``regions="cort"`` keeps only cortical parcels
- ``regions="cort+sub"`` is accepted for API consistency; it is most relevant
  for ``dk``, which includes subcortex

Expression matrices
-------------------

Packaged expression assets are stored as compressed ``.npz`` matrices plus:

- a label table for atlas rows
- a shared gene-label array
- atlas-specific ``README.txt`` and provenance files

Most atlases use a shared AHBA gene list with about ``15.6k`` genes. The
exact count depends on the packaged build for that atlas.

Choosing an atlas
-----------------

Use ``dk`` when you want:

- a legacy-compatible default
- cortical and subcortical regions
- the lowest-resolution preset

Use a Schaefer atlas when you want:

- cortex-only workflows
- a family of matched resolutions
- better control over coarse versus fine parcellation

Use ``destrieux`` when you want:

- a classical cortical atlas distinct from the Schaefer family

Use ``glasser-360`` when you want:

- a surface-oriented cortical workflow
- a finer modern parcellation for cortical maps
