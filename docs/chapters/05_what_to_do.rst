.. _workflows:

=========
Workflows
=========

In this section you can decide which of the available workflows best suits 
your research question and how to run it.

The workflows are grouped below by the kind of question they answer.

.. raw:: html

   <div class="imt-workflow-groups">
     <section class="imt-workflow-group imt-workflow-group--inferential">
       <div class="imt-workflow-group-head">
         <p class="imt-workflow-group-kicker">Inferential workflows</p>
         <h2>Start from an imaging phenotype and test transcriptomic structure against it</h2>
         <p>Use these when your main question begins with a brain map and you want gene-level, component-level, or pathway-level interpretation.</p>
       </div>
       <div class="imt-workflow-card-grid">
         <a class="imt-workflow-card" href="workflows/correlation.html">
           <h3>Correlation</h3>
           <p>Rank genes by how closely their regional expression follows the imaging map.</p>
         </a>
         <a class="imt-workflow-card" href="workflows/pls.html">
           <h3>PLS</h3>
           <p>Find latent multivariate gene components that align with the imaging map.</p>
         </a>
         <a class="imt-workflow-card" href="workflows/enrichment.html">
           <h3>Enrichment</h3>
           <p>Move from ranked genes to pathway- and cell-type-level interpretation.</p>
         </a>
       </div>
     </section>

     <section class="imt-workflow-group imt-workflow-group--signature">
       <div class="imt-workflow-group-head">
         <p class="imt-workflow-group-kicker">Gene-signature workflows</p>
         <h2>Start from genes or weights and summarize them into regional transcriptomic patterns</h2>
         <p>Use these when your input is already a gene list or a weighted gene signature.</p>
       </div>
       <div class="imt-workflow-card-grid">
         <a class="imt-workflow-card" href="workflows/gene_pca.html">
           <h3>Gene-list PCA</h3>
           <p>Reduce a curated gene set into one or more dominant regional expression components.</p>
         </a>
         <a class="imt-workflow-card" href="workflows/gedar.html">
           <h3>GEDAR</h3>
           <p>Project a weighted gene signature, such as a TWAS table, into a regional transcriptomic score.</p>
         </a>
       </div>
     </section>

     <section class="imt-workflow-group imt-workflow-group--gene">
       <div class="imt-workflow-group-head">
         <p class="imt-workflow-group-kicker">Gene utilities</p>
         <h2>Start from one gene and inspect where it is expressed and who travels with it</h2>
         <p>Use this workflow when your question is gene-centric rather than map-centric.</p>
       </div>
       <div class="imt-workflow-card-grid imt-workflow-card-grid--single">
         <a class="imt-workflow-card" href="workflows/gene.html">
           <h3>Single-gene query</h3>
           <p>Extract one gene's regional expression profile and its strongest positive and negative co-expression partners.</p>
         </a>
       </div>
     </section>
   </div>


.. toctree::
   :hidden:
   :maxdepth: 1

   workflows/correlation
   workflows/pls
   workflows/gene
   workflows/gedar
   workflows/gene_pca
   workflows/enrichment
