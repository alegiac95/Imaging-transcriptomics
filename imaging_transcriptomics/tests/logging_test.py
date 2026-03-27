from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

import imaging_transcriptomics.corr as corr_module
import imaging_transcriptomics.genes as genes_module
import imaging_transcriptomics.nulls as nulls_module
from imaging_transcriptomics._logging import get_logger
from imaging_transcriptomics.corr import CorrAnalysis
from imaging_transcriptomics.genes import PLSGenes
from imaging_transcriptomics.models import AtlasSelection, AtlasSpec, ExtractedScan
from imaging_transcriptomics.pls import PLSAnalysis


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


def test_permute_scan_values_logs_start_and_finish(caplog):
    spec = AtlasSpec(
        id="dummy",
        label="Dummy",
        description="Dummy atlas",
        family="test",
        volumetric_space="MNI152",
        supported_spaces=("MNI152",),
        packaged=True,
        buildable=False,
        default_hemisphere="left",
        has_subcortex=True,
        n_regions_left=2,
        n_regions_both=2,
    )
    labels = pd.DataFrame(
        {
            "id": [1, 2],
            "label": ["A", "B"],
            "hemisphere": ["L", "L"],
            "structure": ["subcortex", "subcortex"],
        }
    )
    expression = pd.DataFrame({"id": [1, 2], "Region": ["A", "B"], "GENE1": [0.1, 0.2]})
    selection = AtlasSelection(
        atlas=spec,
        hemisphere="left",
        regions="all",
        labels=labels,
        expression=expression,
        gene_labels=np.array([["GENE1"]], dtype=object),
    )
    extracted = ExtractedScan(values=np.array([1.0, 2.0]), selection=selection, source="array")

    with caplog.at_level(logging.INFO, logger="imaging_transcriptomics.nulls"):
        nulls_module.permute_scan_values(extracted, 3, null_method="random", seed=1234)

    assert "Generating 3 permuted maps with null method 'random'." in caplog.text
    assert "Generating within-hemisphere permutations for 2 non-cortical regions." in caplog.text
    assert "Finished generating permuted maps with resolved method 'random'." in caplog.text


def test_pls_boot_logs_permutation_progress(monkeypatch, caplog):
    def _fake_fit(_gene_exp, imaging_data, n_components, **kwargs):
        del kwargs
        n_regions = np.asarray(imaging_data, dtype=float).shape[0]
        return {
            "x_scores": np.ones((n_regions, n_components), dtype=float),
            "x_weights": np.ones((3, n_components), dtype=float),
            "varexp": np.linspace(0.3, 0.1, n_components, dtype=float),
        }

    monkeypatch.setattr(PLSAnalysis, "_fit_pls", staticmethod(_fake_fit))

    imaging = np.array([0.1, 0.2, 0.3, 0.4], dtype=float)
    gene_exp = np.arange(12, dtype=float).reshape(4, 3)
    analysis = PLSAnalysis(imaging, gene_exp, n_components=2, var=None, n_iter=5, n_jobs=1)
    permuted = np.tile(imaging.reshape(-1, 1), (1, 5))

    with caplog.at_level(logging.INFO, logger="imaging_transcriptomics.pls"):
        analysis.boot_pls(imaging, permuted, analysis._prepared_gene_exp)

    assert "Calculating PLS with permuted data (5 permutations, 1 job)." in caplog.text
    assert "Processed 5/5 PLS permutations." in caplog.text
    assert "Finished PLS permutation fits." in caplog.text
