.. raw:: html

   <div class="imt-index-brand">
     <img src="_static/imaging_transcriptomics_toolbox_logo.svg" alt="Imaging Transcriptomics Toolbox logo" class="only-light">
     <img src="_static/imaging_transcriptomics_toolbox_logo_dark.png" alt="Imaging Transcriptomics Toolbox logo" class="only-dark">
   </div>

.. raw:: html

   <div class="imt-index-badges">
     <a href="https://github.com/alegiac95/Imaging-transcriptomics"><img src="https://img.shields.io/badge/GitHub-repository-172842?style=flat-square&amp;logo=github" alt="GitHub repository"></a>
     <a href="https://imaging-transcriptomics.readthedocs.io/"><img src="https://readthedocs.org/projects/imaging-transcriptomics/badge/?version=latest" alt="Docs status"></a>
     <a href="https://github.com/alegiac95/Imaging-transcriptomics/releases"><img src="https://img.shields.io/github/v/release/alegiac95/Imaging-transcriptomics?style=flat-square" alt="Latest release"></a>
     <a href="https://codecov.io/gh/alegiac95/Imaging-transcriptomics" > <img src="https://codecov.io/gh/alegiac95/Imaging-transcriptomics/branch/refactor-v2.0.0/graph/badge.svg?token=5VA10XQKAY"/></a>
     <a href="https://pypi.org/project/imaging-transcriptomics/"><img src="https://img.shields.io/pypi/v/imaging-transcriptomics?style=flat-square&amp;logo=pypi&amp;logoColor=white" alt="PyPI version"></a>
     <a href="https://github.com/alegiac95/Imaging-transcriptomics/pkgs/container/imaging-transcriptomics"><img src="https://img.shields.io/badge/Docker-GHCR-2496ED?style=flat-square&amp;logo=docker&amp;logoColor=white" alt="Docker image on GHCR"></a>
     <a href="chapters/03_installation.html"><img src="https://img.shields.io/badge/Apptainer-from_GHCR-0E2F5A?style=flat-square&amp;color=0E2F5A&amp;labelColor=EEF4FA" alt="Apptainer and Singularity via GHCR"></a>
     <a href="https://pypi.org/project/imaging-transcriptomics/"><img src="https://img.shields.io/badge/python-%3E%3D3.10-1565C0?style=flat-square&amp;logo=python&amp;logoColor=white" alt="Python version"></a>
     <a href="https://doi.org/10.5281/zenodo.5507505"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.5507505.svg" alt="Zenodo DOI"></a>
     <a href="https://github.com/alegiac95/Imaging-transcriptomics/blob/main/LICENSE"><img src="https://img.shields.io/badge/license-MIT-F28C28?style=flat-square" alt="License"></a>
   </div>


The Imaging Transcriptomics Toolbox helps researchers integrate neuroimaging data with transcriptomic information from the Allen Human Brain Atlas. 
Built to be practical, accessible, and extensible, it provides a streamlined framework for exploring the molecular basis of brain imaging findings 
and developing reproducible imaging transcriptomics workflows.



.. raw:: html

   <section class="imt-section-intro">
     <img src="_static/undraw_outer-space_qey5.svg" alt="Illustration for getting started with the toolbox">
     <div class="imt-section-copy">
       <p class="imt-section-kicker">Start here</p>
       <h2>Take your first steps without getting lost in the details</h2>
       <p>
         If you are new to the toolbox, begin with the quickstart, move through
         installation, and then use the workflow hub to choose the analysis that
         matches your data and question.
       </p>
       <p class="imt-section-links">
         <a href="chapters/01_getting_started.html">Quickstart</a>
         <span>•</span>
         <a href="chapters/03_installation.html">Installation</a>
         <span>•</span>
         <a href="chapters/05_what_to_do.html">Workflow hub</a>
       </p>
     </div>
   </section>


.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Start here

   chapters/01_getting_started
   chapters/03_installation
   chapters/05_what_to_do


.. raw:: html

   <section class="imt-section-intro imt-section-intro--reverse">
     <img src="_static/undraw_data-at-work_3tbf.svg" alt="Illustration for using the toolbox in practice">
     <div class="imt-section-copy">
       <p class="imt-section-kicker">Using the toolbox</p>
       <h2>Work with real inputs, inspect outputs, and move smoothly between the CLI and Python API</h2>
       <p>
         Once you know which workflow you need, this section helps you run the
         toolbox in practice: choose the right interface, prepare the inputs,
         understand the atlases, and interpret the output bundle.
       </p>
       <p class="imt-section-links">
         <a href="chapters/usage/cli.html">CLI</a>
         <span>•</span>
         <a href="chapters/usage/python_api.html">Python API</a>
         <span>•</span>
         <a href="chapters/usage/inputs.html">Inputs</a>
         <span>•</span>
         <a href="chapters/usage/outputs.html">Outputs</a>
       </p>
     </div>
   </section>

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Using the toolbox

   chapters/usage/cli
   chapters/usage/python_api
   chapters/usage/inputs
   chapters/usage/outputs
   chapters/atlases/included_atlases
   chapters/atlases/building_atlases

.. raw:: html

   <section class="imt-section-intro">
     <img src="_static/undraw_problem-solving_1kpx_flipped.svg" alt="Illustration for understanding the methods behind the toolbox">
     <div class="imt-section-copy">
       <p class="imt-section-kicker">Methods</p>
       <h2>Understand the models, nulls, and enrichment logic before you interpret the results</h2>
       <p>
         This section explains how the toolbox computes association scores,
         generates null maps, and evaluates enrichment so you can understand
         what each workflow is actually testing.
       </p>
       <p class="imt-section-links">
         <a href="chapters/methods/statistics.html">Statistics</a>
         <span>•</span>
         <a href="chapters/methods/null_models.html">Null models</a>
         <span>•</span>
         <a href="chapters/methods/gene_sets.html">Gene enrichment</a>
       </p>
     </div>
   </section>

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Methods

   chapters/methods/statistics
   chapters/methods/null_models
   chapters/methods/gene_sets

.. raw:: html

   <section class="imt-section-intro imt-section-intro--reverse">
     <img src="_static/undraw_book-lover_m9n3.svg" alt="Illustration for citations and references" class="imt-section-image--flip">
     <div class="imt-section-copy">
       <p class="imt-section-kicker">Reference</p>
       <h2>Keep the key papers close when you write up or compare results</h2>
       <p>
         Use the reference section to find the core toolbox citations, the
         dedicated GEDAR paper, and broader reading material on imaging
         transcriptomics and companion resources.
       </p>
       <p class="imt-section-links">
         <a href="chapters/reference/citations.html">Reference</a>
         <span>•</span>
         <a href="chapters/reference/further_reading.html">Further reading</a>
       </p>
     </div>
   </section>

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Reference

   chapters/reference/citations
   chapters/reference/further_reading

.. raw:: html

   <section class="imt-section-intro">
     <img src="_static/undraw_version-control_e4yu.svg" alt="Illustration for developing and contributing to the toolbox">
     <div class="imt-section-copy">
       <p class="imt-section-kicker">Development</p>
       <h2>Test, extend, and contribute without guessing how the project fits together</h2>
       <p>
         If you want to work on the package itself, start with the testing
         guide and then move to the contributing notes for coding conventions,
         local workflows, and release-facing project structure.
       </p>
       <p class="imt-section-links">
         <a href="chapters/development/testing.html">Testing</a>
         <span>•</span>
         <a href="chapters/development/contributing.html">Contributing</a>
       </p>
     </div>
   </section>

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Development

   chapters/development/testing
   chapters/development/contributing

.. raw:: html

   <section class="imt-section-intro imt-section-intro--reverse">
     <img src="_static/undraw_questions_52ic.svg" alt="Illustration for frequently asked questions and support">
     <div class="imt-section-copy">
       <p class="imt-section-kicker">FAQ and Contact</p>
       <h2>Get unstuck quickly when something is unclear, unexpected, or hard to reproduce</h2>
       <p>
         This section collects the most common questions about inputs,
         permutations, and workflow interpretation, and it also points to the
         best place to ask for help when something still is not behaving as
         expected.
       </p>
       <p class="imt-section-links">
         <a href="chapters/08_faq.html">FAQ</a>
         <span>•</span>
         <a href="chapters/07_contact_us.html">Contact Us</a>
       </p>
     </div>
   </section>

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: FAQ and Contact

   chapters/08_faq
   chapters/07_contact_us
