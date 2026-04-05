======================
Command-line interface
======================

The command-line interface is the quickest way to run the packaged workflows.
The recommended executable is ``imt``. The longer legacy name
``imagingtranscriptomics`` still works and maps to the same parser.

Commands
--------

- ``imt atlases``: inspect atlas presets
- ``imt genesets``: inspect packaged and remote geneset resources
- ``imt corr``: run correlation analysis
- ``imt pls``: run PLS analysis
- ``imt gene``: query one gene and return its top co-expressed genes
- ``imt gene-pca``: run PCA on atlas expression restricted to a gene list
- ``imt gedar``: run weighted regional GEDAR scoring from a gene-weight table

Common patterns
---------------

The CLI follows a few shared conventions:

- analysis outputs are written to a folder containing ``README.txt``,
  ``metadata.json``, tables, and plots
- ``--atlas``, ``--hemisphere``, and ``--regions`` control which packaged
  expression matrix is selected before analysis
- ``--permutations`` controls the number of null samples used for correlation
  and PLS workflows
- ``--null-method`` controls cortical null generation
- ``--geneset`` is used by both GSEA and ORA

The accepted input families for ``corr`` and ``pls`` are:

- parcel vectors in text, TSV, or CSV format
- volumetric NIfTI maps
- left and right surface files for supported spaces

Raw subject-space anatomical scans are not registered automatically. Native
subject images should be aligned to ``MNI152`` first, then passed with
``--space MNI152``.

Atlas discovery
---------------

List all known atlas presets:

.. code-block:: bash

   imt atlases

List only atlases that are ready to use immediately:

.. code-block:: bash

   imt atlases --packaged-only

Correlation workflow
--------------------

Run a standard correlation analysis on an ``MNI152`` map:

.. code-block:: bash

   imt corr \
     --input /absolute/path/to/map.nii.gz \
     --space MNI152 \
     --atlas dk \
     --hemisphere left \
     --regions all \
     --permutations 20000 \
     --null-method auto \
     --output /absolute/path/to/out_dir

Run ORA only, using a raw gene ``p_value`` threshold:

.. code-block:: bash

   imt corr \
     --input /absolute/path/to/map.nii.gz \
     --space MNI152 \
     --atlas dk \
     --hemisphere left \
     --permutations 20000 \
     --geneset lake \
     --ora-p-threshold 0.05 \
     --output /absolute/path/to/out_dir

Run ORA and GSEA together:

.. code-block:: bash

   imt corr \
     --input /absolute/path/to/map.nii.gz \
     --space MNI152 \
     --atlas dk \
     --hemisphere left \
     --permutations 20000 \
     --geneset lake \
     --ora-p-threshold 0.05 \
     --gsea \
     --output /absolute/path/to/out_dir

PLS workflow
------------

Keep a fixed number of components:

.. code-block:: bash

   imt pls \
     --input /absolute/path/to/map.nii.gz \
     --space MNI152 \
     --atlas dk \
     --ncomp 2 \
     --permutations 20000 \
     --jobs 8 \
     --output /absolute/path/to/out_dir

Choose the number of kept components from a cumulative variance target:

.. code-block:: bash

   imt pls \
     --input /absolute/path/to/map.nii.gz \
     --space MNI152 \
     --atlas dk \
     --var 0.5 \
     --permutations 20000 \
     --output /absolute/path/to/out_dir

Gene PCA workflow
-----------------

Pass a text file with one gene symbol per line:

.. code-block:: bash

   imt gene-pca \
     --genes /absolute/path/to/genes.txt \
     --atlas dk \
     --hemisphere left \
     --ncomp 3 \
     --output /absolute/path/to/out_dir

Pass a comma-separated list directly:

.. code-block:: bash

   imt gene-pca \
     --genes RELN,GAD1,SLC1A2,SV2A \
     --atlas schaefer-100 \
     --hemisphere both \
     --ncomp 3 \
     --output /absolute/path/to/out_dir

Single-gene workflow
--------------------

Query one atlas-aligned gene expression profile and return the top significant
positively co-expressed genes:

.. code-block:: bash

   imt gene \
     --gene RELN \
     --atlas dk \
     --hemisphere left \
     --top-n 25 \
     --output /absolute/path/to/out_dir

Write raw regional expression instead of z-scored regional expression:

.. code-block:: bash

   imt gene \
     --gene MBP \
     --atlas schaefer-200 \
     --hemisphere both \
     --raw-expression \
     --fdr-threshold 0.01 \
     --output /absolute/path/to/out_dir

GEDAR workflow
--------------

Run GEDAR on a TWAS-style weight table:

.. code-block:: bash

   imt gedar \
     --weights /absolute/path/to/twas.tsv \
     --atlas dk \
     --gene-column gene_name \
     --weight-column z_mean \
     --rank-column pvalue \
     --top-percent 5 \
     --direction combined \
     --output /absolute/path/to/out_dir

Return separate up and down GEDAR scores:

.. code-block:: bash

   imt gedar \
     --weights /absolute/path/to/twas.tsv \
     --atlas dk \
     --gene-column gene_name \
     --weight-column z_mean \
     --rank-column pvalue \
     --top-percent 5 \
     --direction split \
     --output /absolute/path/to/out_dir

GSEA and ORA behavior
---------------------

The enrichment defaults are workflow-specific:

- ``corr`` and ``pls`` default to ``ensemble`` unless you request another
  backend
- ``gedar`` defaults to ``none`` and only runs enrichment when you request
  ``--enrichment gsea`` or ``--enrichment ora``
- ``--gsea`` remains as a legacy shortcut for ``--enrichment gsea``
- ``--no-gsea`` disables the legacy GSEA shortcut and leaves the workflow on
  its default or explicit ``--enrichment`` choice

Performance notes
-----------------

- correlation is vectorized and scales well for moderate permutation counts
- PLS permutations are more expensive; use ``--jobs`` to parallelize them
- very large gene-wise correction targets may require hundreds of thousands of
  permutations, which can become expensive even when the workflow is
  optimized
- surface-based null models depend on ``neuromaps`` assets; if they are not
  available, ``auto`` may fall back to random within-hemisphere shuffling

Command reference
-----------------

Supported commands:

- ``imt atlases``
- ``imt genesets``
- ``imt corr``
- ``imt pls``
- ``imt gene``
- ``imt gene-pca``
- ``imt gedar``

Synopsis
~~~~~~~~

.. code-block:: text

   imt atlases [--packaged-only]
   imt genesets [--packaged-only] [--organism NAME]
   imt corr --input PATH [shared options]
   imt pls --input PATH (--ncomp N | --var FRACTION) [shared options]
   imt gene --gene SYMBOL [gene options]
   imt gene-pca --genes VALUE [gene-pca options]
   imt gedar --weights PATH [gedar options]

Shared analysis options
~~~~~~~~~~~~~~~~~~~~~~~

``corr`` and ``pls`` share these main options:

- ``--input``: vector, NIfTI map, or supported surface input
- ``--input-rh``: matching right-hemisphere surface file
- ``--output``: output directory
- ``--atlas``: atlas preset ID
- ``--hemisphere``: ``left`` or ``both``
- ``--regions``: ``default``, ``all``, ``cort``, or ``cort+sub``. ``default``
  and ``cort`` are equivalent, and ``all`` and ``cort+sub`` are equivalent
- ``--space``: source space, such as ``MNI152``, ``fsaverage``, ``fsLR``, or
  ``CIVET``
- ``--permutations``: number of null samples
- ``--seed``: random seed
- ``--jobs``: worker count for parallel work such as PLS permutations
- ``--null-method``: ``auto``, ``vasa``, ``alexander_bloch``, ``moran``, or
  ``random``
- ``--geneset``: packaged gene set name or local ``.gmt`` file
- ``--ora-p-threshold``: raw gene p-value threshold for ORA
- ``--gsea`` / ``--no-gsea``: explicit GSEA control

PLS-specific options
~~~~~~~~~~~~~~~~~~~~

- ``--ncomp``: keep a fixed number of PLS components
- ``--var``: keep enough components to reach a cumulative variance target

Gene-PCA-specific options
~~~~~~~~~~~~~~~~~~~~~~~~~

- ``--genes``: text file, TSV, CSV, or comma-separated list of gene symbols
- ``--atlas`` / ``--hemisphere`` / ``--regions``: same atlas-selection
  semantics as the other workflows
- ``--ncomp``: maximum number of PCA components to retain

Gene-specific options
~~~~~~~~~~~~~~~~~~~~~

- ``--gene``: gene symbol to query in the selected atlas
- ``--atlas`` / ``--hemisphere`` / ``--regions``: same atlas-selection
  semantics as the other workflows
- ``--top-n``: maximum number of significant positively co-expressed genes to
  return
- ``--fdr-threshold``: BH cutoff used to select the top co-expressed genes
- ``--raw-expression``: write raw regional expression instead of z-scored
  regional expression

GEDAR-specific options
~~~~~~~~~~~~~~~~~~~~~~

- ``--weights``: CSV or TSV table containing gene weights
- ``--gene-column``: column containing gene symbols
- ``--weight-column``: column containing the weights used in the GEDAR average
- ``--rank-column``: optional ranking column used for thresholding or top-gene
  selection
- ``--rank-mode``: ``ascending`` or ``descending``
- ``--top-percent`` / ``--top-n`` / ``--p-threshold``: gene-selection rules
- ``--direction``: ``combined``, ``up``, ``down``, or ``split``
- ``--normalize-expression``: ``zscore`` or ``none``
- ``--normalize-weights``: ``none``, ``zscore``, or ``unit``
- ``--enrichment``: ``gsea``, ``ora``, or ``none``. GEDAR does not run
  enrichment unless you request it
- ``--geneset`` / ``--geneset-organism``: enrichment resource for GEDAR GSEA
  or ORA

Reference notes
~~~~~~~~~~~~~~~

- ``imagingtranscriptomics`` is still available as a long-form alias for
  ``imt``
