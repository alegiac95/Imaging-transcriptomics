.. _Gettingstarted:

===============
Getting started
===============

Once the tool is installed you can run the analysis by calling the script from the terminal as:

.. code:: bash

    imaging-transcriptomics -i path-to-your-file.nii -n 1

This is the most simple way to run the script and will permorm the analysis with 1 PLS component on your file and save
the results in a folder named *vh_file_name* in the same path as the original scan file.
It might be that running this will not hold much of the total variance of the scan, however this can be used as a
"first quick estimation". In the resulting path there will be a plot with the variance explained by the first 15
components independently and cumulatively, that can be used to tune consequent analyses, if needed.

For more information on the use have a look at the :ref:`usage <Usage>` page. You can also have a  deeper look at the
:ref:`methods <imgtrans>` and on :ref:`what to do with the results from the script <whatdo>`.

For more advanced use, or to integrate it in your python workflow, you can use the :ref:`python module <library>`.
