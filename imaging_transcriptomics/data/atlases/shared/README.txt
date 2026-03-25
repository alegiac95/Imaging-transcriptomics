Shared atlas gene-label assets
==============================

This directory stores shared AHBA gene label arrays used by the packaged atlas
expression matrices.

Files
-----
- `genes-ahba-15677.npy`
  Shared gene ordering for DK, Schaefer 100, and Schaefer 200.
- `genes-ahba-15675.npy`
  Shared gene ordering for Schaefer 400, Destrieux, and Glasser 360.

Why this exists
---------------
The atlas `.npz` expression archives only store the numeric matrix. Region
metadata lives in each atlas `labels.csv`, and gene names are shared here so
the same list is not duplicated inside every atlas archive.

Runtime contract
----------------
- Expression matrix rows are aligned to the order of each atlas `labels.csv`
- Expression matrix columns are aligned to the shared gene-label file declared
  in the atlas registry
