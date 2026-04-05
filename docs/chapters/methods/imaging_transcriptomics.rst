=======================
Imaging transcriptomics
=======================

Imaging transcriptomics links spatially resolved neuroimaging phenotypes to
regional gene-expression measurements, most commonly from the Allen Human Brain
Atlas. The core question is simple: if a brain map varies across regions, do
any genes, pathways, or weighted gene signatures vary across the same regions
in a similar way?

This general framing has been laid out in a practical guide by
`Arnatkeviciute, Fulcher and Fornito (2019)
<https://doi.org/10.1016/j.neuroimage.2019.01.011>`_ and reviewed in the
context of disease-focused applications by `Arnatkeviciute et al. (2022)
<https://doi.org/10.1016/j.bpsgos.2021.10.002>`_.

The basic idea
--------------

Most imaging-transcriptomic analyses have the same moving parts:

1. one brain phenotype summarized as one value per atlas region
2. one atlas-matched expression matrix with regions by genes
3. one gene-level statistic relating the phenotype to each gene
4. one spatial null model to test whether the association is stronger than
   expected under spatially structured random maps
5. one interpretation layer, often pathway enrichment or gene-signature
   summary

The toolkit in this package follows that same structure. The main workflows
change the gene-level statistic and the final interpretation step, but they all
share the same basic atlas-matching and null-model logic.

What the toolbox does
---------------------

The current workflows cover three common analysis families:

``corr``
   Start from one brain map and rank genes by how closely their regional
   expression follows it.

``pls``
   Start from one brain map and identify multivariate gene components that
   explain structured covariance with that map.

``gene``, ``gene-pca``, and ``gedar``
   Start from one gene, one curated gene list, or one weighted gene signature
   and turn that gene-centric input into regional transcriptomic summaries.

Across those workflows, the package tries to keep three things explicit:

- what is being compared to what
- which null model is actually being used
- which output is descriptive and which output is inferential

Why preprocessing matters
-------------------------

Imaging transcriptomics is sensitive to how the Allen Human Brain Atlas is
processed before any statistical model is fitted. Probe selection, donor
aggregation, missing-data handling, and especially gene normalization can all
change downstream results. The package therefore builds on preprocessed atlas
assets and follows the broader recommendation, articulated by
`Markello et al. (2021) <https://doi.org/10.7554/eLife.72129>`_, that AHBA
processing choices should be made explicit rather than treated as invisible
defaults.

Why spatial nulls matter
------------------------

Brain maps are spatially autocorrelated, so nearby regions often have similar
values. If that structure is ignored, gene-map associations can look more
significant than they really are. This is why the toolbox generates nulls on
the imaging side and reuses them consistently across gene-level, component-
level, and enrichment-level inference. For a broader comparison of the spatial
null families used in the field, see `Markello and Misic (2021)
<https://doi.org/10.1016/j.neuroimage.2021.118052>`_ and the dedicated
:doc:`null-model page <null_models>`.

Why enrichment needs care
-------------------------

Gene-set interpretation is useful, but it is also easy to overstate. In
imaging transcriptomics, enrichment results depend strongly on what is treated
as the gene universe, how the null is defined, and whether the method tests a
ranked signature, a thresholded hit list, or a category score under phenotype
nulls. Those concerns are discussed particularly clearly by
`Fulcher, Arnatkeviciute and Fornito (2021)
<https://doi.org/10.1038/s41467-021-22862-1>`_. The toolbox therefore exposes
multiple enrichment families and documents them separately in
:doc:`enrichment methods <gene_sets>` and in the workflow-facing
:doc:`enrichment guide </chapters/workflows/enrichment>`.

What the toolbox does not claim
-------------------------------

Imaging transcriptomics is best treated as a structured association framework,
not as direct evidence of causality. A strong spatial alignment between a brain
map and a gene, pathway, or cell-type signature does not, by itself, prove
that the gene program causes the imaging phenotype. The outputs are most useful
when interpreted together with prior biology, external genetics, single-cell
data, or independent experimental evidence.

Where to go next
----------------

- :doc:`statistics`: exact definitions of the exported correlation and PLS statistics
- :doc:`null_models`: how cortical and mixed-atlas spatial nulls are generated
- :doc:`gene_sets`: how GSEA, ORA, and ensemble-GCEA differ
