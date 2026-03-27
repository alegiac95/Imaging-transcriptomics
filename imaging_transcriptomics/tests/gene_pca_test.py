from __future__ import annotations

from pathlib import Path

import pandas as pd

from imaging_transcriptomics import load_gene_labels, run_gene_pca
from imaging_transcriptomics.exceptions import PlottingUnavailableError
from imaging_transcriptomics.serialization import write_result_bundle


def test_run_gene_pca_writes_expected_outputs(tmp_path: Path):
    genes = load_gene_labels("dk").reshape(-1).tolist()[:5]
    result = run_gene_pca(
        genes,
        atlas="dk",
        hemisphere="left",
        regions="all",
        n_components=3,
        output_dir=tmp_path,
    )

    assert result.regional_scores.shape[0] == 41
    assert result.variance_table.shape[0] == 3
    assert list(result.variance_table.columns) == [
        "component",
        "variance_explained",
        "cumulative_variance",
    ]
    assert set(result.matched_genes) == set(genes)
    assert (tmp_path / "gene_pca_scores.tsv").exists()
    assert (tmp_path / "gene_pca_loadings.tsv").exists()
    assert (tmp_path / "gene_pca_variance.tsv").exists()
    assert (tmp_path / "metadata.json").exists()
    assert (tmp_path / "README.txt").exists()
    assert (tmp_path / "matched_genes.txt").exists()
    assert (tmp_path / "missing_genes.txt").exists()
    assert (tmp_path / "plots" / "gene_pca_variance.png").exists()
    assert (tmp_path / "plots" / "gene_pca_pc1_regions.png").exists()
    assert (tmp_path / "plots" / "gene_pca_pc1_loadings.png").exists()
    assert (tmp_path / "plots" / "gene_pca_pc3_regions.png").exists()
    assert (tmp_path / "plots" / "gene_pca_pc3_loadings.png").exists()


def test_run_gene_pca_reports_missing_genes(tmp_path: Path):
    genes = load_gene_labels("dk").reshape(-1).tolist()[:3] + ["NOT_A_REAL_GENE"]
    result = run_gene_pca(
        genes,
        atlas="dk",
        hemisphere="left",
        regions="all",
        n_components=2,
        output_dir=tmp_path,
    )

    assert result.missing_genes == ("NOT_A_REAL_GENE",)
    metadata = pd.read_json(tmp_path / "metadata.json", typ="series")
    assert int(metadata["n_genes_missing"]) == 1


def test_run_gene_pca_accepts_file_input(tmp_path: Path):
    genes = load_gene_labels("dk").reshape(-1).tolist()[:4]
    gene_file = tmp_path / "genes.txt"
    gene_file.write_text("\n".join(genes) + "\n")

    result = run_gene_pca(
        gene_file,
        atlas="dk",
        hemisphere="left",
        regions="all",
        n_components=2,
    )

    assert len(result.matched_genes) == 4
    assert result.variance_table.shape[0] == 2


def test_gene_pca_bundle_writes_tables_when_plotting_is_unavailable(tmp_path: Path, monkeypatch):
    import imaging_transcriptomics.plotting as plotting_module

    genes = load_gene_labels("dk").reshape(-1).tolist()[:4]
    result = run_gene_pca(
        genes,
        atlas="dk",
        hemisphere="left",
        regions="all",
        n_components=2,
    )

    monkeypatch.setattr(
        plotting_module,
        "save_result_plots",
        lambda *args, **kwargs: (_ for _ in ()).throw(PlottingUnavailableError("plots unavailable")),
    )

    write_result_bundle(result, tmp_path)

    assert (tmp_path / "gene_pca_scores.tsv").exists()
    assert (tmp_path / "gene_pca_loadings.tsv").exists()
    assert (tmp_path / "gene_pca_variance.tsv").exists()
    assert (tmp_path / "metadata.json").exists()
    assert (tmp_path / "README.txt").exists()
    assert not (tmp_path / "plots" / "gene_pca_variance.png").exists()
    assert "Plot PNGs were skipped" in (tmp_path / "README.txt").read_text()
