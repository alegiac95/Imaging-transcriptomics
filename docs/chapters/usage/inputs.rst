======
Inputs
======

The package accepts a small number of input families and validates them
aggressively before analysis starts.

Input families
--------------

- regional vectors
- volumetric ``MNI152`` NIfTI maps
- supported surface inputs
- gene lists for ``gene-pca``

Regional vectors
----------------

- accepted file types: plain text, TSV, and CSV
- expected length: exactly the number of rows in the selected atlas subset
- expected ordering: the same ordering returned by ``select_atlas_data()``

Volumes
-------

- ``MNI152`` NIfTI maps can be parcellated directly
- native-space subject anatomical scans are not registered automatically
- if the volume is not already aligned to the packaged atlas grid, the package
  will only resample when the input is explicitly declared to be in
  ``MNI152``

Surfaces
--------

- surface workflows expect a left and right hemisphere pair
- a single left file can be combined with ``--input-rh`` in the CLI
- cross-space handling depends on ``neuromaps`` and the required atlas assets

Gene lists
----------

``gene-pca`` accepts:

- a text file
- a TSV or CSV file that can be tokenized into gene symbols
- a comma-separated string
- a Python iterable of symbols through the API

Duplicate genes are removed while preserving the first occurrence.
