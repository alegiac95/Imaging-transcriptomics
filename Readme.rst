=======================
Imaging Transcriptomics
=======================

Imaging transcriptomics is a methodology that allows to identify patterns of correlation between gene expression and
some property of brain structure or function as measured by neuroimaging (e.g., MRI, fMRI, PET).

The package is written in Python3 and requires the some libraries to run:

* numpy
* pandas
* scipy
* nibabel

You can install the package from ``Pypi`` by running:

.. code-block:: shell
    $ pip install imaging-transcriptomics


This Python package contains a command line interface (cli) and a library to import directly in your Python scripts.
To run the cli the command ``imaging-transcriptomics`` is provided while to import it to a Python script you can use

.. code-block:: python
    import imaging_transcriptomics as imt


