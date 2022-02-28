
.. _Usage:

============
Script usage
============

Once you have installed the package you can run the analysis script as:

.. code:: bash

    imagingtranscriptomics --input /path-to-your-in-file [options] {corr|pls [options]}

The script has some options that allow the user to tune the analysis to their specific application. The options are as follows:

- `--input` (`-i`, **mandatory**): path to the input data. This can be either a neuroimaging scan (*i.e.*, .nii[.gz]) or a text file (*i.e.*, .txt).
  
  .. warning::

      If the input scan is a neuroimaging scan (*i.e.*, .nii, .nii.gz) this is expected to be in the same resolution as the Desikan-Killiany (DK) atlas used which is 1mm isotropic (matrix size 182x218x182). On the other hand if the input is a text file, this must be a text file with one column and no headers, with the rows containing the values of interest in the same order as the DK atlas used.

- `--output` (`-o`, **optional**): path to the output directory. If none is provided the results will be saved in the same folder as the input scan.
- `--regions` (`-r`, **optional**): regions to use for the analysis, can be *cort+sub* (or equivalently *all*) which specifies that all the regions are used, or, alternatively, *cort* for the cortical regions only. The latter is useful with some certain types of data, where the subcortical regions might not be available (*e.g.*, EEG).
- `--no-gsea` (**optional**): specifies whether or not Gene Set Enrihment Analysis should be performed.
- `--geneset` (**optional**): specifies the name of the gene set or the path to the file to use for Gene Set Enrichment Analysis.
  
  .. warning:: 

      The `--geneset` argument will be ignored if you also specify the `--no-gsea` flag. If the GSEA analysis is performed, the name of the gene set, or a path to a custom made gene set, should be given. To lookup the name of the available gene sets or on how to create a custom one refer to the GSEA section. 

- 


.. tip::
    All paths given as input should be given as absolute paths instead of relative paths to avoid any errors in reading the file.

.. _library:

=======================
Usage as python library
=======================

Once installed the library can be used like any other Python package in custom written analysis pipelines.
To the library can be imported by running:

.. code:: python

    import imaging_transcriptomics as imt

Once imported the package will contain the core ``ImagingTranscriptomics``
class, along with other useful functions. To see all the available functions
imported in the library run:

.. code:: python

    dir(imt)

which will display all the functions and modules imported in the library.

ImagingTranscriptomics Class
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``ImagingTranscriptomics`` class is the core class of the entire package and allows you to run the entire analysis on your data.
To start using the class the first step is to initialise it. A way to do this is:

.. code:: python

    my_analysis = imt.ImagingTranscriptomics(my_data, n_components=1)

The ``my_data`` is an array that contains the data you want to analyse (e.g., the average values from the scan). This vector has to have some characteristics, mainly:

* it has to be a ``numpy.array`` with length 41, which corresponds to the number of regions in the left hemisphere of the Desikan-Killiany atlas.
* it has to contain the values you want to analyse but not the ``zscore`` of the values as this is computed automatically during the initialisation.

Alternatively to initialise the class with the number of desired components, you can initialise the class by specifying the amount of variance that you want the components to explain. The software will then select the number of components that explains *at least* the specyfied amount (e.g., you specify a 60% of variance and one component explains 58% while the first two components combined explain 70%, two componets will be selcted).

.. code:: python

   my_analysis = imt.ImagingTranscriptomics(my_data, var=0.6)
   # The amount of varinace can be expressed in different ways and gets converted internally.
   # The following will produce the same results as the above
   my_analysis = imt.ImagingTranscriptomics(my_data, var=60)

Once the class in initialised, you can run the analysis by invoking the ``.run()`` method.

.. code:: python

   my_analysis.run()

There are currently two methods to run the analysis, the first uses PLS regression while the other uses Spearman correlation. The PLS analysis is the default method to analyse the data is PLS, while if yoh want to run the analysis with correlation you can run the command:

.. code:: python 

   my_analysis.run(method="corr")
 
.. note:: Please be aware that running the correlation method is currently much slower than the PLS method. This is due the number of correlation that have to be ran during the permutation analysis. The code running these analysis is leveraging multiprocessing of the processor, by using as many cores of the CPU as possible, but even doing this times of *20min* are not uncommon. 

Once the analysis is completed you can check you results by accessing the attributes of the class.


Other Functions of Imaging_Transcriptomics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``imaging_transcriptomics`` library contains several helpful functions, like:

* ``read_scan``: that allows to read a NIfTI file and returns the data matrix (without any of the header information.
* ``extract_average``: that allows to extract the average value from the left hemisphere of the scan.  
