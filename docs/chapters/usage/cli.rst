=====================
Command-line interface
=====================

The command-line interface is the quickest way to run the packaged workflows.
The recommended executable is ``imt``. The longer legacy name
``imagingtranscriptomics`` still works and maps to the same parser.

Commands
--------

- ``imt atlases``: inspect atlas presets
- ``imt genesets``: inspect packaged and remote geneset resources
- ``imt corr``: run correlation analysis
- ``imt pls``: run PLS analysis
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

The enrichment defaults are intentionally simple:

- if you do not pass ``--ora-p-threshold``, the CLI runs GSEA only when
  ``--gsea`` is on or when the workflow defaults require it
- if you pass ``--ora-p-threshold`` and do nothing else, the CLI runs ORA
  only
- if you want both ORA and GSEA, combine ``--ora-p-threshold`` with
  ``--gsea``
- ``--no-gsea`` forces the workflow to skip GSEA

Performance notes
-----------------

- correlation is vectorized and scales well for moderate permutation counts
- PLS permutations are more expensive; use ``--jobs`` to parallelize them
- very large gene-wise correction targets may require hundreds of thousands of
  permutations, which can become expensive even when the workflow is
  optimized
- surface-based null models depend on ``neuromaps`` assets; if they are not
  available, ``auto`` may fall back to random within-hemisphere shuffling
