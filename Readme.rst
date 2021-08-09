=======================
Imaging Transcriptomics
=======================

Imaging transcriptomics is a methodology that allows to identify patterns of correlation between gene expression and
some property of brain structure or function as measured by neuroimaging (e.g., MRI, fMRI, PET).

An overview of the methodology is shown in figure below.

.. image:: docs/chapters/images/imaging_transcriptomics.png
    :alt: alternate-text
    :align: center


In brief, average values of the scan are extracted from 41 brain regions as defined by the Desikan-Killiany (DK) atlas.
Regional values are then used to perform partial least squares (PLS) regression with gene expression data from the
Allen Human Brain Atlas (AHBA) mapped to the DK atlas, in the left hemisphere only.

As a result of the PLS regression we obtain the ranked genes list according to the spatial alignment with the
neuroimaging marker of interest.

The package is written in Python3 and requires the some libraries to run:

* numpy
* pandas
* scipy
* nibabel
* netneurotools

You can install the package from ``Pypi`` by running:

.. code-block:: shell
    $ pip install imaging-transcriptomics


This Python package contains a command line interface (cli) and a library to import directly in your Python scripts.
To run the cli the command ``imaging-transcriptomics`` is provided while to import it to a Python script you can use

.. code-block:: python
    import imaging_transcriptomics as imt


