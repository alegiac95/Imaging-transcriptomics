.. raw:: html

   <div class="imt-index-brand">
     <img src="_static/imaging_transcriptomics_toolbox_logo.svg" alt="Imaging Transcriptomics Toolbox logo" class="only-light">
     <img src="_static/imaging_transcriptomics_toolbox_logo_dark.png" alt="Imaging Transcriptomics Toolbox logo" class="only-dark">
   </div>


The Imaging Transcriptomics Toolbox helps researchers integrate neuroimaging data with transcriptomic information from the Allen Human Brain Atlas. 
Built to be practical, accessible, and extensible, it provides a streamlined framework for exploring the molecular basis of brain imaging findings and developing reproducible imaging transcriptomics workflows.

This documentation is organized around what most users actually need to do:
install the toolbox, choose a workflow, understand the inputs and outputs, and
interpret the resulting statistics.

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
   :caption: Development

   chapters/development/testing
   chapters/development/contributing



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

.. toctree::
   :maxdepth: 1
   :caption: Methods

   chapters/methods/statistics
   chapters/methods/null_models
   chapters/methods/gene_sets

.. toctree::
   :maxdepth: 1
   :caption: Reference and help

   chapters/reference/public_api
   chapters/reference/cli_reference
   chapters/reference/file_formats
   chapters/08_faq
   chapters/07_contact_us

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
