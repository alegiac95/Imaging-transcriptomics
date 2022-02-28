.. _faq:

===
FAQ
===



#. **How can I install the imaging transcriptomics package?**
    The short answer is: you can install it via ``pip``. For more details on how to install refer to the :ref:`installation section. <Installation>`

#. **Why does the analysis use only the left hemisphere?**
    The analysis relies on the left hemisphere only due to the genetic data used. The Allen Human Brain Atlas (AHBA) has a discrepancy in data acquisition between left and right hemisphere resulting in a lot of missing data in the right hemisphere. Given that the brain is not symmetrical, we decided to not mirror data from one hemisphere to the other and constrain the analysis to this hemisphere only.

#. **Why did you use the pypls library instead of some more maintained PLS library, e.g., sklearn?**
    We used pypls instead of sklearn because the latter one, and most of the other available, are implemented using the NIPALS algorithm, while pypls uses the SIMPLS.
    One of the main advantages of the SIMPLS algorithm in respect to the NIPALS is that is is less time consuming.

#. **Can I run the ImaginTranscriptomics analysis on just the cortical areas without the subcortical areas?**
    Yes, check out the main page on the use to get an idea on how to do this.
