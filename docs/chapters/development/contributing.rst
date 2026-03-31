============
Contributing
============

Contributions should keep the repository easier to use, easier to test, and
easier to understand.

Branch and commit hygiene
-------------------------

- keep commits grouped by concern
- prefer small refactors over large mixed commits
- avoid bundling feature work, docs, and unrelated cleanup together unless the
  change is genuinely inseparable

Code style
----------

- prefer the function-first public API exposed from the top-level package
- keep implementation details inside the internal modules
- add docstrings for user-facing functions and short descriptive docstrings for
  internal helpers when they clarify intent
- prefer clear array-oriented implementations in the statistical code paths

Tests and docs
--------------

- add or update tests for any behavior change
- update the docs when a workflow, output schema, or statistical definition
  changes
- keep CLI help, README examples, and docs examples consistent

Proposing a new CLI workflow
----------------------------

If you want to propose a new top-level workflow for ``imt``, it helps to start
from a short design note before opening a large implementation PR.

Good workflow proposals usually answer four questions up front:

- what scientific question the workflow answers
- what the expected inputs and outputs are
- how the workflow differs from the existing commands
- which statistical assumptions or null models it depends on

Proposal skeleton
~~~~~~~~~~~~~~~~~

Contributors can use the following template as a starting point:

.. code-block:: text

   Workflow name
   -------------
   Proposed CLI command:
   Example: imt <new-command>

   1. Scientific goal
   - What does this workflow estimate or test?
   - What kind of user question should it answer?

   2. Inputs
   - Required inputs:
   - Optional inputs:
   - Supported spaces / atlas assumptions:
   - Whether it operates on imaging maps, regional vectors, gene tables, or both:

   3. Outputs
   - Main result tables:
   - Main plots:
   - Metadata or bookkeeping files:

   4. Statistics
   - Core score / model:
   - Null model or resampling strategy:
   - How p values are computed:
   - How multiple-comparison correction is handled:

   5. CLI shape
   - Required flags:
   - Optional flags:
   - Which existing shared flags should be reused:
   - Expected default output directory:

   6. Relationship to existing workflows
   - Why is this not just corr / pls / gene-pca / gedar?
   - Does it share code with an existing workflow?

   7. Validation
   - Small reproducible example:
   - Expected edge cases:
   - Expected tests:

Implementation checklist
~~~~~~~~~~~~~~~~~~~~~~~~

For a new CLI workflow, contributors will usually need to touch most of the
following areas:

- add the parser entry in [parser.py](/Users/alessiogiacomel/Imaging-transcriptomics/imaging_transcriptomics/cli_support/parser.py)
- dispatch the parsed arguments in [runners.py](/Users/alessiogiacomel/Imaging-transcriptomics/imaging_transcriptomics/cli_support/runners.py)
- implement the workflow entry point in [workflows](/Users/alessiogiacomel/Imaging-transcriptomics/imaging_transcriptomics/workflows)
- expose the public function from the top-level package or API layer if it is
  meant to be used from Python
- add tests for parser behavior, workflow behavior, and at least one smoke
  path
- document the command in the CLI guide, reference pages, outputs page, and
  any method pages needed to explain the statistics

Contributor note
~~~~~~~~~~~~~~~~

As a rule of thumb, a new top-level command should only be added when it
represents a genuinely distinct workflow. If the behavior is just a small
variation of an existing method, it is often better to add a new option or
mode to the existing command rather than expanding the CLI surface.
