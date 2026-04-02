# Imaging Transcriptomics 2.0

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6364963.svg)](https://doi.org/10.5281/zenodo.6364963)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Maintainer](https://img.shields.io/badge/maintainer-alegiac95-blue)](https://github.com/alegiac95)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/doc/versions/)
[![Documentation Status](https://readthedocs.org/projects/imaging-transcriptomics/badge/?version=latest)](https://imaging-transcriptomics.readthedocs.io/en/latest/?badge=latest)

<p align="center">
  <img src="docs/_static/imaging_transcriptomics_toolbox_logo.svg" alt="Imaging Transcriptomics Toolbox logo" width="720">
</p>

Imaging transcriptomics links regional neuroimaging maps to Allen Human Brain Atlas gene-expression data. Version `2.0` provides a function-based Python API, a simpler command-line interface, packaged atlas assets, and lightweight text, table, and figure outputs.

The toolbox supports:

- spatial correlation between imaging maps and regional gene expression
- partial least squares (PLS) workflows with atlas-aware null models
- gene set enrichment analysis (GSEA) and over-representation analysis (ORA)
- gene-list PCA
- GEDAR weighted regional expression scoring from TWAS-style gene tables

The GEDAR workflow in this release is inspired by:

- Giacomel A, Powell TR, Duarte RRR, et al. *Transcriptome-informed brain cartography of polygenic risk and association with brain structure in major psychiatric disorders.* Molecular Psychiatry (2026). [https://doi.org/10.1038/s41380-026-03497-4](https://doi.org/10.1038/s41380-026-03497-4)

## Installation

The recommended install routes are the same as in the docs:

- `pip` for the standard Python package
- `uv` for an isolated CLI install
- `Docker` / `Podman` via GHCR for a reproducible container runtime
- `Apptainer` / `Singularity` by pulling from the GHCR OCI image

### pip

```bash
pip install --upgrade pip
pip install imaging-transcriptomics
```

### uv

```bash
uv tool install imaging-transcriptomics
```

### Docker

```bash
docker pull ghcr.io/alegiac95/imaging-transcriptomics:latest
docker run --rm ghcr.io/alegiac95/imaging-transcriptomics:latest --help
```

### Apptainer / Singularity

```bash
apptainer pull imaging-transcriptomics.sif docker://ghcr.io/alegiac95/imaging-transcriptomics:latest
apptainer exec imaging-transcriptomics.sif imt --help
```

### What the default install includes

The standard install already includes the common runtime stack used by the main workflows:

- `matplotlib` for figure generation
- `gseapy` for GSEA
- `neuromaps` for surface inputs, resampling, and cortical null models

Advanced extras are only needed for specific use cases:

- `imaging-transcriptomics[brainspace]` for optional BrainSpace cortical comparison renders
- `imaging-transcriptomics[atlas-build]` for rebuilding atlas assets with `abagen`
- `imaging-transcriptomics[maps]` as a compatibility bundle for both advanced extras
- `imaging-transcriptomics[dev]` for tests and development tools

## Quick start

The CLI entry point is `imt`. The longer `imagingtranscriptomics` command still works, but `imt` is the preferred interface in the docs.

### Choose a workflow

- `imt corr` or `run_corr()` for map-to-gene correlation and optional enrichment
- `imt pls` or `run_pls()` for multivariate gene components
- `imt gene-pca` or `run_gene_pca()` for PCA on a selected gene list
- `imt gedar` or `run_gedar()` for weighted regional transcriptomic scoring
- `imt atlases` and `imt genesets` to inspect packaged resources

### CLI example

```bash
imt corr \
  --input /abs/path/map.nii.gz \
  --space MNI152 \
  --atlas dk \
  --hemisphere left \
  --regions default \
  --permutations 1000 \
  --null-method auto \
  --output /abs/path/out_corr
```

### Python example

```python
import numpy as np
import imaging_transcriptomics as imt

scan = np.linspace(-1.0, 1.0, 41)
result = imt.run_corr(
    scan,
    atlas="dk",
    hemisphere="left",
    regions="default",
    n_permutations=1000,
    output_dir="out_corr",
)
```

## Packaged atlases

The repository ships with packaged atlas assets and expression matrices for:

- `dk`
- `schaefer-100`
- `schaefer-200`
- `schaefer-400`
- `destrieux`
- `glasser-360`

All atlas names use the same region-scope interface:

- `regions="default"` or `regions="cort"` keeps the cortical atlas
- `regions="all"` or `regions="cort+sub"` adds the packaged `aseg` subcortical extension

Expression matrices are derived with `abagen`. For `hemisphere="both"`, the right hemisphere values are obtained through `abagen` mirroring (`lr_mirror="leftright"`).

## Inputs and outputs

Accepted inputs include:

- regional vectors as arrays or text/tabular files
- volumetric NIfTI maps in `MNI152`
- supported surface files through the default `neuromaps` runtime support

Each workflow writes a simple result bundle with:

- `README.txt`
- `metadata.json`
- workflow-specific TSV tables
- matched and missing gene lists when relevant
- plot PNGs in `plots/`

Common examples include:

- `corr_genes.tsv`
- `pls_summary.tsv`
- `gene_pca_scores.tsv`
- `gedar_scores.tsv`

## Documentation

The documentation is organized by user task and is the best place for the full workflow and reference material:

- [Docs home](https://imaging-transcriptomics.readthedocs.io/)
- [Quickstart](https://imaging-transcriptomics.readthedocs.io/en/latest/chapters/01_getting_started.html)
- [Installation](https://imaging-transcriptomics.readthedocs.io/en/latest/chapters/03_installation.html)
- [Workflow hub](https://imaging-transcriptomics.readthedocs.io/en/latest/chapters/05_what_to_do.html)
- [CLI](https://imaging-transcriptomics.readthedocs.io/en/latest/chapters/usage/cli.html)
- [Python API](https://imaging-transcriptomics.readthedocs.io/en/latest/chapters/usage/python_api.html)
- [Methods](https://imaging-transcriptomics.readthedocs.io/en/latest/chapters/methods/statistics.html)
- [Included atlases](https://imaging-transcriptomics.readthedocs.io/en/latest/chapters/atlases/included_atlases.html)

The CLI help is also intentionally descriptive:

```bash
imt --help
imt corr --help
imt pls --help
imt gene-pca --help
imt gedar --help
```

## Troubleshooting

- raw subject T1w scans are not valid direct analysis inputs; use a derived map in `MNI152`, a regional vector, or a supported surface map
- if cortical null generation falls back unexpectedly, set a writable `NEUROMAPS_DATA` cache and confirm `neuromaps` is installed
- if a cloud-synced input under Dropbox or iCloud raises a macOS permission error, copy it to a regular local folder first
- for gene-level FDR with correlation, very low permutation counts can make adjusted p-values collapse to `1`

For problems with the software, please [open an issue on GitHub](https://github.com/alegiac95/Imaging-transcriptomics/issues).

## Development

For local development:

```bash
pip install -e ".[dev]"
pytest -q
```

Useful focused suites:

```bash
pytest -q imaging_transcriptomics/tests/v2_test.py
pytest -q imaging_transcriptomics/tests/golden_test.py
pytest -q imaging_transcriptomics/tests/cli_test.py imaging_transcriptomics/tests/plotting_test.py
```

## Citing

If you publish work using this toolbox, please cite:

- Martins D, Giacomel A, Williams SCR, Turkheimer F, Dipasquale O, Veronese M. *Imaging transcriptomics: Convergent cellular, transcriptomic, and molecular neuroimaging signatures in the healthy adult human brain.* Cell Reports. [https://doi.org/10.1016/j.celrep.2021.110173](https://doi.org/10.1016/j.celrep.2021.110173)
- Giacomel A, Martins D. *Imaging-transcriptomics: Second release update (v1.0.2).* Zenodo. [https://doi.org/10.5281/zenodo.5726839](https://doi.org/10.5281/zenodo.5726839)
- Giacomel A, Martins D, Frigo M, Turkheimer F, Williams SCR, Dipasquale O, Veronese M. *Integrating neuroimaging and gene expression data using the imaging transcriptomics toolbox.* STAR Protocols. [https://doi.org/10.1016/j.xpro.2022.101315](https://doi.org/10.1016/j.xpro.2022.101315)
- Giacomel A, Powell TR, Duarte RRR, et al. *Transcriptome-informed brain cartography of polygenic risk and association with brain structure in major psychiatric disorders.* Molecular Psychiatry (2026). [https://doi.org/10.1038/s41380-026-03497-4](https://doi.org/10.1038/s41380-026-03497-4)
