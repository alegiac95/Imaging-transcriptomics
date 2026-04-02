.. _Gettingstarted:

==========
Quickstart
==========

This page is the fast-entry guide for new users.

It should answer four questions quickly:

1. what the toolbox does
2. which command to run first
3. what kinds of inputs are accepted
4. where to look next in the documentation

Quick start
-----------

If you only want the shortest route through the docs:

1. install the package
2. choose one workflow
3. run a small example
4. use the workflow guide to interpret the outputs

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

- :doc:`Installation </chapters/03_installation>`: install options and advanced extras
- :doc:`Workflow hub </chapters/05_what_to_do>`: choose between correlation, PLS, GEDAR, gene PCA, and enrichment
- :doc:`CLI guide </chapters/usage/cli>` and :doc:`Python API guide </chapters/usage/python_api>`: how to run analyses
- :doc:`Inputs </chapters/usage/inputs>` and :doc:`Outputs </chapters/usage/outputs>`: accepted data and generated files
- :doc:`Methods </chapters/methods/statistics>`: statistical definitions, null models, and enrichment calculations
- :doc:`FAQ </chapters/08_faq>` and :doc:`Contact Us </chapters/07_contact_us>`: troubleshooting and support
- :doc:`Reference </chapters/reference/citations>` and :doc:`Further reading </chapters/reference/further_reading>`: papers to cite and broader imaging transcriptomics resources
