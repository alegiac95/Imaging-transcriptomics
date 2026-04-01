# Containers

This directory contains reproducible runtime definitions for:

- `Dockerfile`: standard Docker runtime built with `pip`
- `Containerfile`: Podman / Buildah runtime built with `pip`
- `Dockerfile.uv`: Docker runtime built from the committed `uv.lock`
- `Singularity.def`: Apptainer / Singularity definition file

Each environment installs the toolbox with the runtime extras used by the CLI:

- `plots`
- `gsea`
- `maps`

Typical local commands:

```bash
docker build -f containers/Dockerfile -t imaging-transcriptomics:docker .
docker run --rm imaging-transcriptomics:docker --help
```

```bash
docker build -f containers/Dockerfile.uv -t imaging-transcriptomics:uv .
docker run --rm imaging-transcriptomics:uv --help
```

```bash
podman build -f containers/Containerfile -t imaging-transcriptomics:podman .
podman run --rm imaging-transcriptomics:podman --help
```

```bash
apptainer build imaging-transcriptomics.sif containers/Singularity.def
apptainer exec imaging-transcriptomics.sif imt --help
```
