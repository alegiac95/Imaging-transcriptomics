# Imaging Transcriptomics 2.0

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6364963.svg)](https://doi.org/10.5281/zenodo.6364963)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/doc/versions/)

`imaging-transcriptomics` links brain maps to Allen Human Brain Atlas gene-expression data. Version `2.0.0` keeps the package lighter and easier to use: simpler function-based entry points, built-in atlas handling, optional `neuromaps` resampling, and text-plus-plot outputs instead of PDF reports.

## What Changed In V2

- Functional entry points: `run_corr()` and `run_pls()`
- Included atlases with `abagen`-derived expression data
- Left-only or both-hemisphere expression data
- Direct support for vectors, text tables, NIfTI files, and surface files
- `README.txt` plus TSV tables and plot PNGs as the default outputs
- Local SIMPLS-based PLS backend
- A smaller v2-only public API centered on `RunConfig`, `run_corr()`, and `run_pls()`

## Included Atlases

Included in this branch and ready to use:

- `dk`
- `schaefer-100`
- `schaefer-200`
- `schaefer-400`
- `destrieux`
- `glasser-360`

The expression matrices are still generated from `abagen`. For two-hemisphere analyses, the right side comes from the `abagen` mirror option (`lr_mirror="leftright"`).

## Installation

Create an environment first if you want an isolated install:

```bash
conda env create -f environment-v2.yml
conda activate imaging-transcriptomics-v2
```

Install the package:

```bash
pip install -e .
```

Optional extras:

- `pip install -e .[gsea]` for GSEA support
- `pip install -e .[maps]` for `neuromaps` and `abagen`
- `pip install -e .[dev]` for tests and tooling

## Python API

Correlation:

```python
import numpy as np
import imaging_transcriptomics as imt

scan = np.linspace(-1.0, 1.0, 41)
result = imt.run_corr(
    scan,
    atlas="dk",
    hemisphere="left",
    regions="all",
    n_permutations=1000,
    output_dir="out_corr",
)
```

PLS:

```python
import numpy as np
import imaging_transcriptomics as imt

scan = np.linspace(-1.0, 1.0, 83)
result = imt.run_pls(
    scan,
    atlas="dk",
    hemisphere="both",
    regions="all",
    n_components=2,
    n_permutations=1000,
    output_dir="out_pls",
)
```

Inspect available atlases:

```python
import imaging_transcriptomics as imt

print(imt.atlas_table(packaged_only=True))
print(imt.describe_atlas("dk"))
```

## CLI

List atlases:

```bash
imagingtranscriptomics atlases --packaged-only
```

Run correlation:

```bash
imagingtranscriptomics corr \
  --input /abs/path/scan.nii.gz \
  --atlas dk \
  --hemisphere left \
  --null-method vasa \
  --regions all \
  --output /abs/path/out_dir
```

Run PLS:

```bash
imagingtranscriptomics pls \
  --input /abs/path/scan.nii.gz \
  --atlas schaefer-100 \
  --hemisphere both \
  --null-method auto \
  --ncomp 2 \
  --output /abs/path/out_dir
```

Available null methods are `auto`, `vasa`, `alexander_bloch`, `moran`, and `random`. The default `auto` mode tries `vasa` first for cortical data and falls back to random shuffling within each hemisphere if a surface-based method is not available locally.

## Inputs And Resampling

Version 2 accepts:

- regional vectors as NumPy arrays or text tables
- volumetric NIfTI data in `MNI152`
- non-MNI or surface data through `neuromaps` when the `maps` extra is installed

For maps already in `MNI152`, the package can extract region values directly from the included atlas image. For other standard spaces or surface inputs, `neuromaps` can resample the input before analysis.

## Outputs

Each run writes:

- `README.txt`
- `metadata.json`
- `regional_values.tsv`
- analysis tables such as `corr_genes.tsv`, `pls_summary.tsv`, `pls_component_<n>.tsv`
- plot PNGs in `plots/`

If GSEA is enabled, matching `gsea_*.tsv` tables are also written.

PDF reporting was intentionally removed in this refactor.

## Development

Run the targeted test suite used for this refactor:

```bash
pytest -q imaging_transcriptomics/tests/v2_test.py
pytest -q imaging_transcriptomics/tests/pvalues_test.py
pytest -q imaging_transcriptomics/tests/auto_test.py imaging_transcriptomics/tests/plotting_test.py
```

## Citing

If you publish work using this toolbox, please cite:

- Martins D, Giacomel A, Williams SCR, Turkheimer F, Dipasquale O, Veronese M. *Imaging transcriptomics: Convergent cellular, transcriptomic, and molecular neuroimaging signatures in the healthy adult human brain.* Cell Reports. [https://doi.org/10.1016/j.celrep.2021.110173](https://doi.org/10.1016/j.celrep.2021.110173)
- Giacomel A, Martins D. *Imaging-transcriptomics: Second release update (v1.0.2).* Zenodo. [https://doi.org/10.5281/zenodo.5726839](https://doi.org/10.5281/zenodo.5726839)
- Giacomel A, Martins D, Frigo M, Turkheimer F, Williams SCR, Dipasquale O, Veronese M. *Integrating neuroimaging and gene expression data using the imaging transcriptomics toolbox.* STAR Protocols. [https://doi.org/10.1016/j.xpro.2022.101315](https://doi.org/10.1016/j.xpro.2022.101315)
