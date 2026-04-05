.. _Gettingstarted:

==========
Quickstart
==========

Imaging transcriptomics links spatially resolved neuroimaging phenotypes to
regional gene-expression atlases, most often the Allen Human Brain Atlas. In
practice, that usually means asking whether the spatial pattern in a brain map
aligns with the spatial pattern of one gene, many genes, or a broader pathway,
and then checking that alignment against spatially informed null models rather
than naive shuffles. The overall framing used in this toolbox follows the
practical guide by `Arnatkeviciute, Fulcher and Fornito (2019)
<https://doi.org/10.1016/j.neuroimage.2019.01.011>`_ and the broader review by
`Arnatkeviciute et al. (2022) <https://doi.org/10.1016/j.bpsgos.2021.10.002>`_.

This page is the fast-entry guide for new users.

It should answer four questions quickly:

1. what the toolbox does
2. which command to run first
3. what kinds of inputs are accepted
4. where to look next in the documentation



Quick start
-----------

If you only want the shortest route through the docs:

1. :doc:`Install the package </chapters/03_installation>`
2. :doc:`Choose one workflow </chapters/05_what_to_do>`
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

Query one gene and return its regional expression plus the top co-expressed genes:

.. code:: bash

    imt gene --gene RELN --atlas dk --top-n 25 --output /path/to/out

Accepted input families
-----------------------

- regional vectors aligned to an included atlas
- volumetric NIfTI maps in ``MNI152``
- supported surface files when ``neuromaps`` support is installed
- one gene symbol for ``gene``
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
- :doc:`Workflow hub </chapters/05_what_to_do>`: choose between correlation, PLS, single-gene queries, GEDAR, gene PCA, and enrichment
- :doc:`CLI guide </chapters/usage/cli>` and :doc:`Python API guide </chapters/usage/python_api>`: how to run analyses
- :doc:`Inputs </chapters/usage/inputs>` and :doc:`Outputs </chapters/usage/outputs>`: accepted data and generated files
- :doc:`Imaging transcriptomics </chapters/methods/imaging_transcriptomics>`, :doc:`statistics </chapters/methods/statistics>`, :doc:`null models </chapters/methods/null_models>`, and :doc:`enrichment methods </chapters/methods/gene_sets>`: conceptual background and the statistical choices used by the toolbox
- :doc:`FAQ </chapters/08_faq>` and :doc:`Contact Us </chapters/07_contact_us>`: troubleshooting and support
- :doc:`Reference </chapters/reference/citations>` and :doc:`Further reading </chapters/reference/further_reading>`: papers to cite and broader imaging transcriptomics resources
