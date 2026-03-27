.. _Gettingstarted:

===============
Getting started
===============

This page is the fast-entry guide for new users.

It should answer four questions quickly:

1. what the toolbox does
2. which command to run first
3. what kinds of inputs are accepted
4. where to look next in the documentation

Quick start
-----------

Run a simple correlation analysis:

.. code:: bash

    imt corr --input /path/to/your-map.nii.gz --atlas dk --output /path/to/out

Run a simple PLS analysis:

.. code:: bash

    imt pls --input /path/to/your-map.nii.gz --atlas dk --ncomp 1 --no-gsea --output /path/to/out

Run a simple gene-list PCA analysis:

.. code:: bash

    imt gene-pca --genes RELN,GAD1,SLC1A2,SV2A --atlas dk --ncomp 2 --output /path/to/out

Accepted input families
-----------------------

- regional vectors aligned to an included atlas
- volumetric NIfTI maps in ``MNI152``
- supported surface files when ``neuromaps`` support is installed
- gene lists for ``gene-pca``

What to expect from an output folder
------------------------------------

Most workflows write:

- ``README.txt``
- ``metadata.json``
- one or more TSV tables
- one or more PNG plots

Documentation map
-----------------

- :ref:`Installation <Installation>`: install options and extras
- :ref:`Usage <Usage>`: CLI, API, inputs, and outputs
- :ref:`Workflows <workflows>`: end-to-end guides for correlation, PLS, enrichment, and gene PCA
- :ref:`Methods <methods>`: statistical choices and null models
- :ref:`Reference <reference>`: command and output schemas
