
.. _Usage:

============
Script usage
============

Once you have installed the package you can run the script as

.. code:: bash

    imaging-transcriptomics -i path_to_file -n 2


The script has some parameters that can be tuned according to the necessity of the analysis, which can be viewed by calling the help function from the script as
:code:`imaging-transcriptomics -h` or :code:`imaging-transcriptomics --help`.


Here we will describe in more details all the parameters of the script.
The parameters of the script are:

``-i`` (``--input``)     Scan on which you want to perform the analysis. It is *recommended* that you provide an absolute path to your scan (e.g., :code:`~/Desktop/myscan.nii.gz`) instead of a relative one (e.g., :code:`myscan.nii.gz`) to avoid errors. The input file *must* be an imaging file in NIfTI format (both :code:`.nii` and :code:`.nii.gz` formats are supported), be in MNI152 space and have a matrix dimension of 182x218x182. If your image has a different matrix size you can reslice it to match this dimension with your preferred method. (A quick method is to use *fslview* and reslice to match the dimension of the included *MNI152_1mm* brain atlas).

.. warning:: The input scan must have a predefined dimension **(182x218x182)** and be in **MNI152** space. If the input scan is not in the required dimension the script will throw an error. You should always check the dimension before running the script and eventually reslice or spatially normalise your image to the matching dimensions with your preferred method (e.g., SPM, FSL, ANTS).

``-n`` (``--ncomp``)     Number of PLS components to use for the analysis. The parameter *must* be an **integer** between 1 and 15, otherwise an error will occur. *Please note* that in PLS regression the first component is not necessarily the component explaining the most amount of variance, as in PCA. Example: running :code:`imaging-transcriptomics -i path_to_file -n 2` will run the script on your imaging file selecting the two first components for the analysis.


``-v`` (``--variance``)  Total amount of variance you want your components to explain. The code will automatically select the number of components that explain at least the variance you specify. The parameter *must* be an **integer** between 10 and 100, which represents the percentage of explained variance. Example: if you run :code:`imaging-transcriptomics -i path_to_file -v 30` and the first 3 components explain 10%, 25% and 3%, respectively, of the total variance the script will use 2 components, even if they explain 35% (which is a bit more than specified) of the total variance.

.. warning:: Please note that the **-v** and **-n** parameters are mutually exclusive and only one has to be provided as argument for the script, otherwise an error will occur.

Optional additional parameters that can be provided are:

``-o`` (``--output``)   Path where you want to save the results, if no path is provided the results will be saved in the same path as the input scan. When the code is finished running the results will be saved in a folder named *vh_myscanname* and will contain all the results (.csv files, .pdf report and images in .png format). If you run the script multiple times you will have more than one results folder with a trailing number for the run (e.g., *vh_myscanname* for the first run and *vh_myscanname_1* for the second run).

``--verbose`` Sets the output logging level to debug mode and shows all debug values and steps

``--suppress`` Sets the logging level to warning and will display only eventual warning messages.
