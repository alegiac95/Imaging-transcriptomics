================
Included atlases
================

The toolbox ships six atlas presets with packaged ``abagen`` expression assets.
The same atlas IDs work in the command-line interface, the Python API, and the
docs examples.

Atlas gallery
-------------

The previews below show the cortical surface coverage of each packaged atlas as
single left-lateral surface renders. They are intended to help you compare
parcel scale and cortical coverage quickly, not to encode biology.

.. grid:: 1 1 2 2
   :gutter: 2

   .. grid-item-card:: ``dk``
      :img-top: ../images/atlases/dk_gallery.png

      Desikan-Killiany with ``83`` bilateral regions (``41`` left-only),
      including cortex and subcortex. Best when you want a low-resolution,
      legacy-compatible default with whole-brain coverage.

      - Coverage: cortex + subcortex
      - Spaces: ``MNI152``, ``fsaverage``

   .. grid-item-card:: ``schaefer-100``
      :img-top: ../images/atlases/schaefer-100_gallery.png

      Coarse Schaefer parcellation with ``100`` bilateral cortical parcels.
      Good for cortex-only workflows when you want a compact, easy-to-interpret
      regional map.

      - Coverage: cortex
      - Spaces: ``MNI152``, ``fsaverage``

   .. grid-item-card:: ``schaefer-200``
      :img-top: ../images/atlases/schaefer-200_gallery.png

      Mid-resolution Schaefer preset with ``200`` bilateral cortical parcels.
      A good default when you want more regional detail without moving to the
      finest cortical scales.

      - Coverage: cortex
      - Spaces: ``MNI152``, ``fsaverage``, ``fsLR``

   .. grid-item-card:: ``schaefer-400``
      :img-top: ../images/atlases/schaefer-400_gallery.png

      Fine Schaefer preset with ``400`` bilateral cortical parcels. Best for
      cortex-only analyses where higher regional granularity matters more than
      simplicity.

      - Coverage: cortex
      - Spaces: ``MNI152``, ``fsaverage``, ``fsLR``

   .. grid-item-card:: ``destrieux``

      Destrieux cortical atlas with ``148`` bilateral regions (``74`` left-only).
      The packaged preset includes expression and label assets, but it does not
      currently ship a surface parcellation file for the gallery preview.

      - Coverage: cortex
      - Spaces: ``MNI152``, ``fsaverage``

   .. grid-item-card:: ``glasser-360``
      :img-top: ../images/atlases/glasser-360_gallery.png

      Surface-oriented multimodal atlas with ``360`` bilateral cortical
      parcels. Use it when you want a modern high-resolution cortical workflow,
      especially for surface-first analyses.

      - Coverage: cortex
      - Spaces: ``fsLR``, ``fsaverage``, ``MNI152``

Atlas comparison
----------------

.. list-table::
   :header-rows: 1
   :widths: 18 10 10 16 20 26

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
