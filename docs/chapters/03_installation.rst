.. _Installation:

============
Installation
============

This page describes the stable installation routes for the ``2.0`` toolbox.

For most users, the recommended paths are:

- install from PyPI with ``pip``
- install the CLI into an isolated environment with ``uv``
- run the published container images if you want a reproducible runtime with
  no local Python setup
- use extras only for advanced atlas-building or BrainSpace-specific features

Requirements
------------

Core requirements:

- Python ``3.10`` or newer for ``pip`` and ``uv`` installs
- a writable user environment
- Docker, Podman, or Apptainer/Singularity only if you plan to use containers

Advanced extras:

- ``brainspace`` for optional BrainSpace-rendered cortical comparison plots
- ``atlas-build`` for ``abagen``-based atlas asset building
- ``maps`` as a compatibility bundle that installs both advanced extras

Stable install routes
---------------------

.. tab-set::

   .. tab-item:: pip

      This is the simplest route for most users.

      Recommended install:

      .. code-block:: bash

         pip install --upgrade pip
         pip install imaging-transcriptomics

      Upgrade to the latest published release:

      .. code-block:: bash

         pip install --upgrade imaging-transcriptomics

      Add the optional advanced extras only if you need them:

      .. code-block:: bash

         pip install "imaging-transcriptomics[brainspace]"
         pip install "imaging-transcriptomics[atlas-build]"

   .. tab-item:: uv

      ``uv`` is a good choice if you want the CLI in an isolated tool
      environment without managing a project virtual environment by hand.

      Recommended install:

      .. code-block:: bash

         uv tool install imaging-transcriptomics

      Upgrade an existing tool install:

      .. code-block:: bash

         uv tool upgrade imaging-transcriptomics

      If you prefer to work from a project-local ``uv`` environment instead of
      a tool install:

      .. code-block:: bash

         uv venv
         uv pip install imaging-transcriptomics

      Add the optional advanced extras only if you need them:

      .. code-block:: bash

         uv tool install "imaging-transcriptomics[brainspace]"
         uv pip install "imaging-transcriptomics[atlas-build]"

   .. tab-item:: Docker / Podman

      Published container images are useful when you want a reproducible
      runtime with the standard CLI dependencies already installed.

      Standard image:

      .. code-block:: bash

         docker pull ghcr.io/alegiac95/imaging-transcriptomics:latest
         docker run --rm ghcr.io/alegiac95/imaging-transcriptomics:latest --help

      Podman uses the same OCI image:

      .. code-block:: bash

         podman pull ghcr.io/alegiac95/imaging-transcriptomics:latest
         podman run --rm ghcr.io/alegiac95/imaging-transcriptomics:latest --help

      A locked ``uv``-based image can also be published for fully pinned
      environments:

      .. code-block:: bash

         docker pull ghcr.io/alegiac95/imaging-transcriptomics-uv:latest
         docker run --rm ghcr.io/alegiac95/imaging-transcriptomics-uv:latest --help

   .. tab-item:: Apptainer / Singularity

      Apptainer and Singularity use the same GHCR OCI image as Docker and
      Podman. Instead of publishing a separate ``.sif`` artifact, the toolbox
      release flow publishes the OCI image once and lets cluster users derive a
      local ``.sif`` from it.

      Pull a local image from GHCR:

      .. code-block:: bash

         apptainer pull imaging-transcriptomics.sif docker://ghcr.io/alegiac95/imaging-transcriptomics:latest
         apptainer exec imaging-transcriptomics.sif imt --help

      Pull a specific tagged release:

      .. code-block:: bash

         apptainer pull imaging-transcriptomics.sif docker://ghcr.io/alegiac95/imaging-transcriptomics:v2.0.0

      The same approach works with Singularity installations that support
      ``docker://`` sources. The repository still ships
      ``containers/Singularity.def`` for advanced users who want a native
      definition-based rebuild.

What The Default Install Includes
---------------------------------

The standard package install already includes the common runtime stack used by
the main workflows:

- Matplotlib-backed figure generation for cortical, volumetric, enrichment,
  and summary outputs
- ``gseapy`` for preranked GSEA and remote GMT library access
- ``neuromaps`` support for surface inputs, cross-space resampling, and
  cortical spatial nulls

Use the advanced extras only when you specifically need them:

- ``brainspace`` for optional comparison renders such as
  ``*_cortex_brainspace.png``
- ``atlas-build`` for rebuilding atlas expression assets with ``abagen``
- ``maps`` if you want the legacy one-step compatibility bundle for both
  advanced extras

Source install for development
------------------------------

You only need a source install if you want to edit the code, run the tests,
build the docs, or work on new atlas assets.

Editable ``pip`` install:

.. code-block:: bash

   git clone https://github.com/alegiac95/Imaging-transcriptomics.git
   cd Imaging-transcriptomics
   python -m venv .venv
   source .venv/bin/activate
   pip install --upgrade pip
   pip install -e ".[dev]"

If you prefer ``uv`` for development:

.. code-block:: bash

   git clone https://github.com/alegiac95/Imaging-transcriptomics.git
   cd Imaging-transcriptomics
   uv sync --extra dev

If you also want the advanced build/rendering dependencies in a development
environment:

.. code-block:: bash

   pip install -e ".[dev,maps]"
   uv sync --extra dev --extra maps

Quick validation
----------------

After installation, these are good first checks:

.. code-block:: bash

   imt --help
   imt atlases --packaged-only
   imt corr --help

Platform notes
--------------

- surface-based null models and surface inputs are part of the standard
  install through ``neuromaps``
- published GHCR containers already include the standard runtime dependencies
- BrainSpace comparison renders require the ``brainspace`` extra
- atlas-building utilities require the ``atlas-build`` extra
- Apptainer and Singularity users pull from the GHCR OCI image rather than a
  separately published ``.sif`` artifact
- cloud-synced folders such as Dropbox or iCloud may trigger macOS permission
  issues when used directly as inputs; copying files to a regular local folder
  is often the simplest fix
- native-space anatomical images are not automatically registered; use derived
  maps in ``MNI152`` or regional vectors
