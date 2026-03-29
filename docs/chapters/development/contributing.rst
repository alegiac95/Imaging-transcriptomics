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
