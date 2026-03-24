.. _Installation:

============
Installation
============

``imaging-transcriptomics`` 2.0 targets Python ``3.10+`` and installs directly from PyPI or from a local checkout.

.. tip::

    Install the package in a dedicated virtual environment or conda environment to keep optional neuroimaging dependencies isolated from other projects.

Minimal installation:

.. code:: bash

    pip install imaging-transcriptomics

Optional extras:

.. code:: bash

    pip install imaging-transcriptomics[gsea]
    pip install imaging-transcriptomics[maps]
    pip install imaging-transcriptomics[dev]

For a full local development environment from the repository root:

.. code:: bash

    pip install -e .[dev,maps,gsea]

Or with conda:

.. code:: bash

    conda env create -f environment-v2.yml
    conda activate imaging-transcriptomics-v2

The package now ships with a local SIMPLS-based PLS backend, so no external ``pyls`` install is required.
