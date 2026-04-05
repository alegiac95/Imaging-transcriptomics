================
Included atlases
================

The toolbox ships six base atlas presets with packaged ``abagen`` expression
assets. The same atlas IDs work in the command-line interface, the Python API,
and the docs examples.

Atlas gallery
-------------

The previews below show the cortical surface coverage of each packaged atlas as
single left-lateral surface renders. They are intended to help you compare
parcel scale and cortical coverage quickly, not to encode biology. All atlas
IDs below now use the same simple interface: ``regions="default"`` or
``regions="cort"`` keeps the cortical atlas, while ``regions="all"`` or
``regions="cort+sub"`` appends the packaged ``aseg`` subcortical add-on.

.. container:: imt-atlas-gallery

   .. grid:: 1 1 2 2
      :gutter: 2

      .. grid-item-card:: ``dk``
         :img-top: ../images/atlases/dk_gallery.png

         Legacy-compatible Desikan-Killiany atlas. ``default`` returns the
         cortical parcels; ``all`` appends the packaged ``aseg`` subcortical
         add-on, for ``83`` bilateral regions total.

         - Coverage: cortex by default, optional ``aseg`` add-on
         - Spaces: ``MNI152``, ``fsaverage``

      .. grid-item-card:: ``schaefer-100``
         :img-top: ../images/atlases/schaefer-100_gallery.png

         Coarse Schaefer parcellation with ``100`` bilateral cortical parcels.
         ``all`` appends the shared ``aseg`` subcortical add-on for ``115``
         bilateral regions total.

         - Coverage: cortex by default, optional ``aseg`` add-on
         - Spaces: ``MNI152``, ``fsaverage``

      .. grid-item-card:: ``schaefer-200``
         :img-top: ../images/atlases/schaefer-200_gallery.png

         Mid-resolution Schaefer preset with ``200`` bilateral cortical parcels.
         ``all`` appends the shared ``aseg`` subcortical add-on for ``215``
         bilateral regions total.

         - Coverage: cortex by default, optional ``aseg`` add-on
         - Spaces: ``MNI152``, ``fsaverage``, ``fsLR``

      .. grid-item-card:: ``schaefer-400``
         :img-top: ../images/atlases/schaefer-400_gallery.png

         Fine Schaefer preset with ``400`` bilateral cortical parcels. ``all``
         appends the shared ``aseg`` subcortical add-on for ``415`` bilateral
         regions total.

         - Coverage: cortex by default, optional ``aseg`` add-on
         - Spaces: ``MNI152``, ``fsaverage``, ``fsLR``

      .. grid-item-card:: ``destrieux``
         :img-top: ../images/atlases/destrieux_gallery.png

         Destrieux cortical atlas with ``148`` bilateral regions (``74`` left-only).
         ``all`` appends the shared ``aseg`` subcortical add-on for ``163``
         bilateral regions total.

         - Coverage: cortex by default, optional ``aseg`` add-on
         - Spaces: ``MNI152``, ``fsaverage``

      .. grid-item-card:: ``glasser-360``
         :img-top: ../images/atlases/glasser-360_gallery.png

         Surface-oriented multimodal atlas with ``360`` bilateral cortical
         parcels. ``all`` appends the shared ``aseg`` subcortical add-on for
         ``375`` bilateral regions total.

         - Coverage: cortex by default, optional ``aseg`` add-on
         - Spaces: ``fsLR``, ``fsaverage``, ``MNI152``

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

- ``regions="default"`` keeps only cortical parcels
- ``regions="cort"`` is an explicit synonym for ``default``
- ``regions="all"`` keeps cortex plus the packaged ``aseg`` subcortical add-on
- ``regions="cort+sub"`` is an explicit synonym for ``all``

Subcortical add-on
------------------

All packaged atlas names can now expose the same shared ``aseg``-derived
subcortical parcels through the region scope, rather than through separate
atlas IDs.

This means:

- the atlas name stays the same
- ``default`` gives you the historical cortical atlas
- ``all`` adds the same subcortical parcel set used in ``dk``
- old ``*-aseg`` names are accepted as compatibility aliases, but they resolve
  to the base atlas name

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
- easy access to both cortical and aseg-derived subcortical regions
- the lowest-resolution preset

Use a Schaefer atlas when you want:

- cortex-first workflows with an optional subcortical extension
- a family of matched resolutions
- better control over coarse versus fine parcellation

Use ``destrieux`` when you want:

- a classical cortical atlas distinct from the Schaefer family

Use ``glasser-360`` when you want:

- a surface-oriented cortical workflow
- a finer modern parcellation for cortical maps
