.. _Installation:

============
Installation
============

This page is the canonical install guide for the refactored ``2.0`` toolbox.

Requirements
------------

Core requirements:

- Python ``3.10`` or newer
- a working scientific Python stack
- a writable environment if you plan to use editable installs

Optional extras:

- ``gsea`` for GSEA support
- ``maps`` for ``neuromaps``-based resampling and surface nulls
- ``dev`` for tests and development tools

Recommended repository install
------------------------------

Clone the repository and switch to the active branch:

.. code-block:: bash

   git clone https://github.com/alegiac95/Imaging-transcriptomics.git
   cd Imaging-transcriptomics
   git checkout refactor-v2.0.0

Create the packaged conda environment:

.. code-block:: bash

   conda env create -f environment-v2.yml
   conda activate imaging-transcriptomics-v2

Install the package:

.. code-block:: bash

   pip install -e .

Optional extras
---------------

.. code-block:: bash

   pip install -e .[gsea]
   pip install -e .[maps]
   pip install -e .[dev]

Quick validation
----------------

After installation, these are good first checks:

.. code-block:: bash

   imt --help
   imt atlases --packaged-only
   pytest -q

Platform notes
--------------

- surface-based null models require the ``maps`` extra and a writable
  ``neuromaps`` cache
- cloud-synced folders such as Dropbox or iCloud may trigger macOS permission
  issues when used directly as inputs; copying files to a regular local folder
  is often the simplest fix
- native-space anatomical images are not automatically registered; use derived
  maps in ``MNI152`` or regional vectors
