FROM python:3.11-slim

ENV MPLBACKEND=Agg \
    PIP_NO_CACHE_DIR=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

ARG APP_DIR=/opt/imaging-transcriptomics

WORKDIR ${APP_DIR}

RUN apt-get update \
    && apt-get install --yes --no-install-recommends ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY pyproject.toml README.md ./
COPY imaging_transcriptomics imaging_transcriptomics

RUN python -m pip install --upgrade pip \
    && python -m pip install .[plots,gsea,maps]

RUN useradd --create-home --uid 1000 imt \
    && chown -R imt:imt ${APP_DIR}

USER imt

ENTRYPOINT ["imt"]
CMD ["--help"]

