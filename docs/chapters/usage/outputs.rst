=======
Outputs
=======

Each workflow writes a small self-contained output bundle designed for both
manual inspection and downstream scripting.

Core output patterns
--------------------

- ``README.txt``
- ``metadata.json``
- TSV analysis tables
- plot PNGs

Workflow-specific outputs
-------------------------

- correlation gene tables
- PLS summary and per-component tables
- GSEA and ORA tables
- gene PCA score, loading, and variance tables

What is always present
----------------------

- ``README.txt``
- ``metadata.json``
- ``regional_values.tsv`` for ``corr`` and ``pls``
- at least one plot in ``plots/``

What depends on the workflow
----------------------------

``corr``
   Writes ``corr_genes.tsv`` and optional GSEA or ORA tables.

``pls``
   Writes ``pls_summary.tsv``, one ``pls_component_<n>.tsv`` per kept
   component, and optional enrichment outputs for each component.

``gene-pca``
   Writes PCA score, loading, and variance tables plus matched and missing
   gene lists.

See also
--------

For column-level details, use the dedicated reference page:
``reference/file_formats``.
