=======
Testing
=======

The repository uses ``pytest`` for both targeted and full-suite validation.

Core checks
-----------

Run the full suite:

.. code-block:: bash

   pytest -q

Focused suites
--------------

Useful targeted runs on this branch include:

.. code-block:: bash

   pytest -q imaging_transcriptomics/tests/v2_test.py
   pytest -q imaging_transcriptomics/tests/pvalues_test.py
   pytest -q imaging_transcriptomics/tests/golden_test.py
   pytest -q imaging_transcriptomics/tests/cli_test.py imaging_transcriptomics/tests/plotting_test.py

Coverage
--------

The repository includes a ``.coveragerc`` file configured for the current v2
layout. Typical usage:

.. code-block:: bash

   coverage erase
   coverage run -m pytest
   coverage report -m
   coverage html

Atlas and workflow smoke tests
------------------------------

For changes that affect analysis outputs, useful manual checks include:

- ``imt atlases --packaged-only``
- one short ``imt corr`` run
- one short ``imt pls`` run
- one short ``imt gene-pca`` run
- one short ``imt gedar`` run
