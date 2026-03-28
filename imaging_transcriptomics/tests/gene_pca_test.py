from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from imaging_transcriptomics import load_gene_labels, run_gene_pca, select_atlas_data
from imaging_transcriptomics.exceptions import PlottingUnavailableError
from imaging_transcriptomics.gene_expression import load_brain_gene_symbols
from imaging_transcriptomics.serialization import write_result_bundle


def _genes_in_brain_filter(n: int) -> list[str]:
    brain_genes = set(load_brain_gene_symbols())
    genes = load_gene_labels("dk").reshape(-1).tolist()
    selected = [gene for gene in genes if gene.upper() in brain_genes]
    if len(selected) < n:
        raise AssertionError("Not enough DK genes were found in the packaged brain gene filter.")
    return selected[:n]


def test_run_gene_pca_writes_expected_outputs(tmp_path: Path):
    genes = _genes_in_brain_filter(5)
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
    assert result.brain_filtered_genes == ()
    assert (tmp_path / "gene_pca_scores.tsv").exists()
    assert (tmp_path / "gene_pca_loadings.tsv").exists()
    assert (tmp_path / "gene_pca_variance.tsv").exists()
    assert (tmp_path / "metadata.json").exists()
    assert (tmp_path / "README.txt").exists()
    assert (tmp_path / "matched_genes.txt").exists()
    assert (tmp_path / "brain_filtered_genes.txt").exists()
    assert (tmp_path / "missing_genes.txt").exists()
    assert (tmp_path / "plots" / "gene_pca_variance.png").exists()
    assert (tmp_path / "plots" / "gene_pca_pc1_brain.png").exists()
    assert (tmp_path / "plots" / "gene_pca_pc1_regions.png").exists()
    assert (tmp_path / "plots" / "gene_pca_pc1_loadings.png").exists()
    assert (tmp_path / "plots" / "gene_pca_pc3_brain.png").exists()
    assert (tmp_path / "plots" / "gene_pca_pc3_regions.png").exists()
    assert (tmp_path / "plots" / "gene_pca_pc3_loadings.png").exists()


def test_run_gene_pca_reports_missing_genes(tmp_path: Path):
    genes = _genes_in_brain_filter(3) + ["NOT_A_REAL_GENE"]
    result = run_gene_pca(
        genes,
        atlas="dk",
        hemisphere="left",
        regions="all",
        n_components=2,
        output_dir=tmp_path,
    )

    assert result.missing_genes == ()
    assert result.brain_filtered_genes == ("NOT_A_REAL_GENE",)
    metadata = pd.read_json(tmp_path / "metadata.json", typ="series")
    assert int(metadata["n_genes_missing"]) == 0
    assert int(metadata["n_genes_filtered"]) == 1


def test_run_gene_pca_accepts_file_input(tmp_path: Path):
    genes = _genes_in_brain_filter(4)
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

    genes = _genes_in_brain_filter(4)
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
    assert (tmp_path / "brain_filtered_genes.txt").exists()
    assert not (tmp_path / "plots" / "gene_pca_variance.png").exists()
    assert "Plot PNGs were skipped" in (tmp_path / "README.txt").read_text()


def test_run_gene_pca_handles_nonfinite_gene_columns(monkeypatch):
    selection = select_atlas_data(atlas="dk", hemisphere="left", regions="all")
    expression = selection.expression.copy()
    brain_genes = set(load_brain_gene_symbols())
    genes = [gene for gene in expression.columns[2:].tolist() if gene.upper() in brain_genes][:3]
    expression.loc[:, genes[0]] = np.nan
    expression.loc[:, genes[1]] = np.inf
    expression.loc[:, genes[2]] = 1.0

    monkeypatch.setattr(
        "imaging_transcriptomics.gene_pca.select_atlas_data",
        lambda **kwargs: selection.__class__(
            atlas=selection.atlas,
            hemisphere=selection.hemisphere,
            regions=selection.regions,
            labels=selection.labels,
            expression=expression,
            gene_labels=selection.gene_labels,
        ),
    )

    result = run_gene_pca(genes, atlas="dk", hemisphere="left", regions="all", n_components=2)

    assert result.variance_table.shape[0] == 2
    assert np.isfinite(result.regional_scores.filter(like="PC").to_numpy(dtype=float)).all()
    assert np.isfinite(result.gene_loadings.filter(like="PC").to_numpy(dtype=float)).all()


def test_run_gene_pca_tracks_brain_filtered_genes(tmp_path: Path):
    brain_genes = set(load_brain_gene_symbols())
    all_genes = load_gene_labels("dk").reshape(-1).tolist()
    inside = next(gene for gene in all_genes if gene.upper() in brain_genes)
    outside = next(gene for gene in all_genes if gene.upper() not in brain_genes)

    result = run_gene_pca(
        [inside, outside],
        atlas="dk",
        hemisphere="left",
        regions="all",
        n_components=1,
        output_dir=tmp_path,
    )

    assert result.matched_genes == (inside,)
    assert result.brain_filtered_genes == (outside,)
    metadata = pd.read_json(tmp_path / "metadata.json", typ="series")
    assert int(metadata["n_genes_filtered"]) == 1
    assert (tmp_path / "brain_filtered_genes.txt").read_text().strip() == outside
