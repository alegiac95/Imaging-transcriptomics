.. _workflows:

=========
Workflows
=========

In this section you can decide which of the available workflows best suits 
your research question and how to run it.

Once you decide which workflow you want to use, simply click on the respective 
icon on the figure below.

.. only:: html

   .. raw:: html

      <div style="margin: 1.5rem 0 2rem 0; max-width: 1100px;">

   .. raw:: html
      :file: ../_static/workflow_hub.svg

   .. raw:: html

      </div>

.. only:: not html

   .. image:: images/imaging_transcriptomics.png
      :alt: Workflow overview
      :width: 100%

Each workflow page should explain:

- what the workflow is for
- what inputs it needs
- the main CLI and Python entry points
- the key outputs
- how to interpret the results
- the most common pitfalls

.. toctree::
   :hidden:
   :maxdepth: 1

   workflows/correlation
   workflows/pls
   workflows/gedar
   workflows/gene_pca
   workflows/enrichment
