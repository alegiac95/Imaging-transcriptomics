.. _faq:

===
FAQ
===

Which workflow should I use?
----------------------------

Use ``corr`` when you want a ranked gene list for one imaging map.

Use ``pls`` when you want multivariate components rather than one-gene-at-a-time
associations.

Use ``gene-pca`` when you already have a gene list and want its regional
expression pattern.

Use ``gedar`` when you have a weighted gene table, such as a TWAS-style result,
and want a regional weighted-expression score.

How many permutations should I run?
-----------------------------------

Enough for your inferential goal:

- a few thousand are often enough for workflow testing and broad pathway work
- tens of thousands give more stable nominal p-values
- gene-level BH correction across roughly ``15,677`` genes may require hundreds
  of thousands of permutations for small adjusted values to become attainable

Why are many corrected p-values equal to 1?
-------------------------------------------

Usually because permutation resolution is too coarse relative to the number of
tests. This is especially common for gene-level BH correction when the number
of genes is much larger than the number of permutations.

What is the difference between ``fdr`` and ``maxT``?
----------------------------------------------------

``fdr`` controls the expected false discovery rate across many tests.

``maxT`` is a family-wise correction based on the strongest null statistic in
each permutation and is therefore much stricter.

Why did my raw T1w image fail as an input?
------------------------------------------

The toolbox expects a meaningful derived imaging phenotype or a regional vector,
not a native-space anatomical scan used directly as the analysis map.

Use a map already in ``MNI152`` or first transform your data into a valid
regional or voxelwise phenotype in standard space.

Why did a cortical null fall back to random shuffling?
------------------------------------------------------

Usually because the required ``neuromaps`` assets are unavailable locally or
the local ``neuromaps`` installation or cache is unavailable.

Check:

- ``pip install --upgrade imaging-transcriptomics``
- a writable ``NEUROMAPS_DATA`` cache
- that the selected atlas has compatible surface assets for the requested null
  method

Can I use a remote geneset library?
-----------------------------------

Yes. ``--geneset`` accepts:

- packaged sets such as ``lake`` and ``pooled``
- local ``.gmt`` files
- Enrichr/GSEApy library names

Use ``imt genesets --organism Human`` to inspect available remote libraries.
