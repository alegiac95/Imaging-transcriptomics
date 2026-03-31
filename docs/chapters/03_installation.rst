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

Requirements
------------

Core requirements:

- Python ``3.10`` or newer for ``pip`` and ``uv`` installs
- a writable user environment
- Docker, Podman, or Apptainer/Singularity only if you plan to use containers

Optional extras:

- ``plots`` for saved figures and report-style plot outputs
- ``gsea`` for GSEA and geneset-based enrichment workflows
- ``maps`` for ``neuromaps`` resampling, surface inputs, and spatial nulls

Stable install routes
---------------------

.. tab-set::

   .. tab-item:: pip

      This is the simplest route for most users.

      Recommended full install:

      .. code-block:: bash

         python -m pip install --upgrade pip
         python -m pip install "imaging-transcriptomics[plots,gsea,maps]"

      Minimal core install:

      .. code-block:: bash

         python -m pip install imaging-transcriptomics

      Upgrade to the latest published release:

      .. code-block:: bash

         python -m pip install --upgrade "imaging-transcriptomics[plots,gsea,maps]"

   .. tab-item:: uv

      ``uv`` is a good choice if you want the CLI in an isolated tool
      environment without managing a project virtual environment by hand.

      Recommended full install:

      .. code-block:: bash

         uv tool install "imaging-transcriptomics[plots,gsea,maps]"

      Upgrade an existing tool install:

      .. code-block:: bash

         uv tool upgrade imaging-transcriptomics

      If you prefer to work from a project-local ``uv`` environment instead of
      a tool install:

      .. code-block:: bash

         uv venv
         uv pip install "imaging-transcriptomics[plots,gsea,maps]"

   .. tab-item:: Docker / Podman

      Published container images are useful when you want a reproducible
      runtime with the CLI and runtime extras already installed.

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

      If you use Apptainer or Singularity on shared compute systems, you can
      pull a local ``.sif`` image from the published OCI image:

      .. code-block:: bash

         apptainer pull imaging-transcriptomics.sif docker://ghcr.io/alegiac95/imaging-transcriptomics:latest
         apptainer exec imaging-transcriptomics.sif imt --help

      The same approach works with Singularity installations that support
      ``docker://`` sources.

What The Extras Add
-------------------

The published package can be installed in a minimal core form, but most users
will want the runtime extras:

- ``plots`` adds Matplotlib-backed figure generation for cortical, volumetric,
  enrichment, and summary outputs
- ``gsea`` adds ``gseapy`` for preranked GSEA and remote GMT library access
- ``maps`` adds ``neuromaps``, ``abagen``, and ``brainspace`` support for
  spatial resampling, surface workflows, and atlas-building utilities

If you are unsure, install:

.. code-block:: bash

   python -m pip install "imaging-transcriptomics[plots,gsea,maps]"

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
   python -m pip install --upgrade pip
   python -m pip install -e ".[dev,plots,gsea,maps]"

If you prefer ``uv`` for development:

.. code-block:: bash

   git clone https://github.com/alegiac95/Imaging-transcriptomics.git
   cd Imaging-transcriptomics
   uv sync --extra dev --extra plots --extra gsea --extra maps

Quick validation
----------------

After installation, these are good first checks:

.. code-block:: bash

   imt --help
   imt atlases --packaged-only
   imt corr --help

Platform notes
--------------

- surface-based null models and surface inputs require the ``maps`` extra
- published containers already include the runtime extras used by the CLI
- cloud-synced folders such as Dropbox or iCloud may trigger macOS permission
  issues when used directly as inputs; copying files to a regular local folder
  is often the simplest fix
- native-space anatomical images are not automatically registered; use derived
  maps in ``MNI152`` or regional vectors
