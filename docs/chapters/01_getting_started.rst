.. _Gettingstarted:

===============
Getting started
===============

Once the tool is installed, a simple run looks like:

.. code:: bash

    imagingtranscriptomics pls --input /path/to/your-map.nii.gz --atlas dk --ncomp 1 --no-gsea

This command runs the local PLS backend on the included DK atlas and writes a small output folder next to the input:

- ``README.txt``
- ``metadata.json``
- TSV result tables
- PNG plots

For a quick correlation workflow:

.. code:: bash

    imagingtranscriptomics corr --input /path/to/your-map.nii.gz --atlas dk --null-method auto

If your data are already reduced to atlas regions, you can pass a vector file instead of a NIfTI image. If your data are in a different standard space, install the ``maps`` extra and provide ``--space`` so ``neuromaps`` can resample them first.

For more detail see the :ref:`usage <Usage>` page. The package can also be used directly from Python via :ref:`the library API <library>`.
