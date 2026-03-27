from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

import imaging_transcriptomics.corr as corr_module
import imaging_transcriptomics.genes as genes_module
from imaging_transcriptomics._logging import get_logger
from imaging_transcriptomics.corr import CorrAnalysis
from imaging_transcriptomics.genes import PLSGenes


def test_package_logger_uses_namespaced_library_safe_defaults():
    logger = get_logger("genes")
    package_logger = logging.getLogger("imaging_transcriptomics")

    assert logger.name == "imaging_transcriptomics.genes"
    assert any(isinstance(handler, logging.NullHandler) for handler in package_logger.handlers)


def test_corr_ora_logs_progress_and_saving(tmp_path: Path, monkeypatch, caplog):
    analysis = CorrAnalysis(n_iterations=2, n_genes=2)
    analysis.gene_results.results.genes[:, 0] = ["GENE1", "GENE2"]
    analysis.gene_results.results.corr[0, :] = [0.4, -0.3]
    analysis.gene_results.results.pval[0, :] = [0.01, 0.02]

    monkeypatch.setattr(
        corr_module,
        "ora_from_gene_table",
        lambda *args, **kwargs: {"up": pd.DataFrame({"Term": ["Path"]}), "down": pd.DataFrame({"Term": ["Path"]})},
    )

    with caplog.at_level(logging.INFO, logger="imaging_transcriptomics.corr"):
        analysis.ora(gene_set="lake", outdir=tmp_path, p_threshold=0.05)

    assert "Performing ORA." in caplog.text
    assert "Saving ORA results." in caplog.text


def test_pls_ora_logs_progress_and_component_saves(tmp_path: Path, monkeypatch, caplog):
    genes = PLSGenes(n_components=1, n_iter=2, n_genes=2)
    genes.orig.genes[0, :] = ["GENE1", "GENE2"]
    genes.orig.zscored[0, :] = [1.2, -1.1]
    genes.boot.pval[0, :] = [0.01, 0.02]

    monkeypatch.setattr(
        genes_module,
        "ora_from_gene_table",
        lambda *args, **kwargs: {"up": pd.DataFrame({"Term": ["Path"]}), "down": pd.DataFrame({"Term": ["Path"]})},
    )

    with caplog.at_level(logging.INFO, logger="imaging_transcriptomics.genes"):
        genes.ora(gene_set="lake", outdir=tmp_path, p_threshold=0.05)

    assert "Performing ORA." in caplog.text
    assert "Saving ORA results for PLS component 1." in caplog.text
