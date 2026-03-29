.. _workflows:

=========
Workflows
=========

This section is the home for end-to-end task guides.

At a glance, the toolbox has two imaging-map workflows, two gene-centered
workflows, and one shared enrichment layer.

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
   :maxdepth: 1

   workflows/correlation
   workflows/pls
   workflows/gene_pca
   workflows/gedar
   workflows/enrichment
