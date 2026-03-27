===================
Spatial null models
===================

The null-model layer is used by the correlation and PLS workflows to generate
permuted imaging maps that preserve some spatial structure.

Null methods to explain
-----------------------

- ``auto``
- ``vasa``
- ``alexander_bloch``
- ``moran``
- ``random``

How the package uses them
-------------------------

- cortical rows use the requested null model when possible
- non-cortical rows use grouped random shuffling within hemisphere labels
- the same permuted imaging matrix is then reused by the downstream workflow

``auto``
--------

``auto`` is the default because it gives the best available cortical nulls
without forcing the user to know the underlying asset requirements.

The current order is:

1. try ``vasa``
2. if that fails, try ``alexander_bloch``
3. if surface nulls are unavailable, fall back to grouped random shuffling

``vasa`` and ``alexander_bloch``
--------------------------------

These are surface-based cortical nulls provided through ``neuromaps``. They
require:

- the ``maps`` optional dependency
- a compatible surface parcellation for the selected atlas
- surface geometry assets when required by the null implementation

``moran``
---------

``moran`` is also exposed through ``neuromaps`` and is available for
surface-based cortical null generation in the same framework.

``random``
----------

``random`` shuffles values within hemisphere groups. It is the simplest and
least spatially informed option, but it is useful:

- as an explicit baseline
- when surface assets are unavailable
- for debugging or lightweight smoke tests

Cortex versus subcortex
-----------------------

Only cortical rows are candidates for the surface-based null models.
Subcortical or other non-cortical rows are shuffled within hemisphere groups.
This means mixed atlases such as ``dk`` can use a hybrid strategy:

- cortical rows from a spatial null model
- subcortical rows from grouped random shuffling

Practical guidance
------------------

- use ``auto`` unless you need a specific null family
- use ``random`` when you want to avoid surface-null dependencies entirely
- expect a warning when ``auto`` falls back from surface nulls to random
  shuffling
