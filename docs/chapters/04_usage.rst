
.. _Usage:

============
Script usage
============

Once you have installed the package you can run the script as

.. code:: bash

    imagingtranscriptomics -i path_to_file -n 2


The script has some parameters that can be tuned according to the necessity of the analysis, which can be viewed by calling the help function from the script as
:code:`imagingtranscriptomics -h` or :code:`imagingtranscriptomics --help`.


Here we will describe in more details all the parameters of the script.
The parameters of the script are:

``-i`` (``--input``)     Scan on which you want to perform the analysis. It is *recommended* that you provide an absolute path to your scan (e.g., :code:`~/Desktop/myscan.nii.gz`) instead of a relative one (e.g., :code:`myscan.nii.gz`) to avoid errors. The input file *must* be an imaging file in NIfTI format (both :code:`.nii` and :code:`.nii.gz` formats are supported), be in MNI152 space and have a matrix dimension of 182x218x182. If your image has a different matrix size you can reslice it to match this dimension with your preferred method. (A quick method is to use *fslview* and reslice to match the dimension of the included *MNI152_1mm* brain atlas).

.. warning:: The input scan must have a predefined dimension **(182x218x182)** and be in **MNI152** space. If the input scan is not in the required dimension the script will throw an error. You should always check the dimension before running the script and eventually reslice or spatially normalise your image to the matching dimensions with your preferred method (e.g., SPM, FSL, ANTS).

``-n`` (``--ncomp``)     Number of PLS components to use for the analysis. The parameter *must* be an **integer** between 1 and 15, otherwise an error will occur. *Please note* that in PLS regression the first component is not necessarily the component explaining the most amount of variance, as in PCA. Example: running :code:`imaging-transcriptomics -i path_to_file -n 2` will run the script on your imaging file selecting the two first components for the analysis.


``-v`` (``--variance``)  Total amount of variance you want your components to explain. The code will automatically select the number of components that explain at least the variance you specify. The parameter *must* be an **integer** between 10 and 100, which represents the percentage of explained variance. Example: if you run :code:`imaging-transcriptomics -i path_to_file -v 30` and the first 3 components explain 10%, 25% and 3%, respectively, of the total variance the script will use 2 components, even if they explain 35% (which is a bit more than specified) of the total variance.

.. warning:: Please note that the **-v** and **-n** parameters are mutually exclusive and only one has to be provided as argument for the script, otherwise an error will occur.

Optional additional parameters that can be provided are:

``-o`` (``--out``)   Path where you want to save the results, if no path is provided the results will be saved in the same path as the input scan. When the code is finished running the results will be saved in a folder named *Imt_myscanname* and will contain all the results (.csv files, .pdf report and images in .png format). If you run the script multiple times you will have more than one results folder with a trailing number for the run (e.g., *Imt_myscanname* for the first run and *Imt_myscanname_1* for the second run).

``--verbose`` Sets the output logging level to debug mode and shows all debug values and steps

``--suppress`` Sets the logging level to warning and will display only eventual warning messages.


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

..code:: python

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


