# Imaging Transcriptomics 2.0

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6364963.svg)](https://doi.org/10.5281/zenodo.6364963)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/doc/versions/)

`imaging-transcriptomics` links brain maps to Allen Human Brain Atlas gene-expression data. This branch is a `2.0.0` refactor with a lighter functional API, atlas-aware scan extraction, optional neuromaps-based resampling, and text-plus-plot outputs instead of PDF reports.

## What Changed In V2

- Functional entry points: `run_corr()` and `run_pls()`
- Packaged atlas presets with abagen-derived expression data
- Left-only or mirrored left+right hemisphere expression matrices
- Direct vector, text-table, NIfTI, and surface-input handling
- `README.txt` plus TSV tables and plot PNGs as the default outputs
- Local SIMPLS-based PLS backend
- A smaller v2-only public API centered on `RunConfig`, `run_corr()`, and `run_pls()`

## Packaged Atlases

Ready to run in this branch:

- `dk`
- `schaefer-100`

Preset definitions included for local abagen builds:

- `schaefer-200`
- `schaefer-400`
- `destrieux`
- `glasser-360`

The packaged expression matrices are still generated from `abagen`. For bilateral analyses, the right hemisphere is represented through the mirrored expression strategy supported by `abagen` (`lr_mirror="leftright"`).

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

- `pip install -e .[gsea]` for gene-set enrichment analysis
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

Inspect atlas presets:

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

Available null methods are `auto`, `vasa`, `alexander_bloch`, `moran`, and `random`. The default `auto` mode prefers `vasa` for cortical parcellated data and falls back to within-hemisphere random shuffles when a surface null model is unavailable locally.

## Inputs And Resampling

The v2 scan layer accepts:

- regional vectors as NumPy arrays or text tables
- volumetric NIfTI data in `MNI152`
- non-MNI or surface data through `neuromaps` when the `maps` extra is installed

For volumetric maps already in `MNI152`, the package can extract regional values directly with the packaged atlas image. For cross-space or surface workflows, `neuromaps` is used to resample/parcellate the input before analysis.

## Outputs

Each run writes:

- `README.txt`
- `metadata.json`
- `regional_values.tsv`
- analysis tables such as `corr_genes.tsv`, `pls_summary.tsv`, `pls_component_<n>.tsv`
- plot PNGs in `plots/`

If GSEA is enabled, the corresponding `gsea_*.tsv` tables are also written.

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
