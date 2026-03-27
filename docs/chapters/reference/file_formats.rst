===================
Output file formats
===================

The toolbox writes plain-text outputs that are easy to inspect and easy to
parse from downstream scripts.

Core files
----------

``README.txt``
   Human-readable run summary describing the method, atlas, settings, and
   main outputs.

``metadata.json``
   Machine-readable summary of run settings such as atlas, hemisphere, region
   scope, null method, permutation count, and enrichment configuration.

``regional_values.tsv``
   Regional values aligned to the selected atlas subset.

Correlation outputs
-------------------

``corr_genes.tsv``
   Ranked gene table with:

   - ``gene``
   - ``score``
   - ``p_value``
   - ``fdr``
   - ``fwer_maxT``

``gsea_corr_results.tsv``
   Correlation GSEA results with enrichment scores, normalized enrichment
   scores, nominal p-values, and GSEA-style FDR.

``ora_corr_up.tsv`` and ``ora_corr_down.tsv``
   ORA results for positive and negative gene subsets.

PLS outputs
-----------

``pls_summary.tsv``
   Per-component summary with explained variance, cumulative variance, and
   component-level permutation p-values.

``pls_component_<n>.tsv``
   Ranked per-component gene table with:

   - ``gene``
   - ``weight``
   - ``zscore``
   - ``p_value``
   - ``fdr``
   - ``fwer_maxT``

``gsea_pls<n>_results.tsv``
   GSEA results for one PLS component.

``ora_pls<n>_up.tsv`` and ``ora_pls<n>_down.tsv``
   ORA results for one PLS component.

Gene PCA outputs
----------------

``gene_pca_scores.tsv``
   Regional PCA scores with one ``PC<n>`` column per retained component.

``gene_pca_loadings.tsv``
   Per-gene loadings with one ``PC<n>`` column per retained component.

``gene_pca_variance.tsv``
   Variance explained and cumulative variance by component.

``matched_genes.txt`` and ``missing_genes.txt``
   The genes from the requested list that were found or not found in the atlas
   expression matrix.
