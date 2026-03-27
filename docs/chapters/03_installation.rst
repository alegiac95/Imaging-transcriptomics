.. _Installation:

============
Installation
============

This page should become the canonical install guide for all supported setups.

Core requirements
-----------------

- Python ``3.10+``
- optional extras for GSEA and map resampling
- a writable environment for ``neuromaps`` caches if surface nulls are used

Recommended install paths
-------------------------

Minimal install:

.. code:: bash

    pip install imaging-transcriptomics

Repository checkout:

.. code:: bash

    pip install -e .

Optional extras:

.. code:: bash

    pip install -e .[gsea]
    pip install -e .[maps]
    pip install -e .[dev]

Conda environment
-----------------

.. code:: bash

    conda env create -f environment-v2.yml
    conda activate imaging-transcriptomics-v2

What this page should later document in more detail
---------------------------------------------------

- platform-specific notes
- optional dependency tradeoffs
- troubleshooting for ``neuromaps`` caches and surface assets
- editable installs for development
- doc-build and test dependencies
