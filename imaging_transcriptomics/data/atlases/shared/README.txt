Shared atlas gene-label files
=============================

This directory stores shared AHBA gene-name arrays used by the included atlas
expression matrices.

Files
-----
- `genes-ahba-15677.npy`
  Shared gene ordering for DK, Schaefer 100, and Schaefer 200.
- `genes-ahba-15675.npy`
  Shared gene ordering for Schaefer 400, Destrieux, and Glasser 360.

Why this exists
---------------
The atlas `.npz` files only store the numeric matrix. Region metadata stays in
each atlas `labels.csv`, and gene names live here so the same list is not
copied into every atlas file.

Runtime contract
----------------
- Expression matrix rows follow the order of each atlas `labels.csv`
- Expression matrix columns follow the shared gene-label file listed in the
  atlas registry
