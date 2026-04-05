from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from imaging_transcriptomics import load_gene_labels, run_gedar, select_atlas_data
from imaging_transcriptomics.gene_expression import load_brain_gene_symbols
from imaging_transcriptomics.workflows import gedar as gedar_workflow


def _genes_in_brain_filter(n: int) -> list[str]:
    brain_genes = set(load_brain_gene_symbols())
    genes = load_gene_labels("dk").reshape(-1).tolist()
    selected = [gene for gene in genes if gene.upper() in brain_genes]
    if len(selected) < n:
        raise AssertionError("Not enough DK genes were found in the packaged brain gene filter.")
    return selected[:n]


def _write_test_gmt(path: Path, rows: dict[str, list[str]]) -> Path:
    lines = [
        "\t".join([term, "na", *genes])
        for term, genes in rows.items()
    ]
    path.write_text("\n".join(lines) + "\n")
    return path


def test_run_gedar_matches_manual_weighted_average():
    genes = _genes_in_brain_filter(4)
    weights = pd.DataFrame(
        {
            "gene": genes,
            "weight": [1.0, 2.0, -1.5, 0.5],
            "fdr": [0.20, 0.01, 0.10, 0.03],
        }
    )

    result = run_gedar(
        weights,
        atlas="dk",
        hemisphere="left",
        regions="all",
        gene_column="gene",
        weight_column="weight",
        rank_column="fdr",
        top_n=2,
        direction="combined",
        normalize_expression="zscore",
        normalize_weights="none",
    )

    selection = select_atlas_data(atlas="dk", hemisphere="left", regions="all", zscore_expression=True)
    kept_genes = [genes[1], genes[3]]
    expected_weights = np.array([2.0, 0.5], dtype=float)
    expected = (
        selection.expression.loc[:, kept_genes].to_numpy(dtype=float) @ expected_weights
    ) / expected_weights.sum()

    np.testing.assert_allclose(result.regional_scores["score"].to_numpy(dtype=float), expected)
    assert result.matched_genes == tuple(genes)
    assert result.missing_genes == ()
    assert result.excluded_table.empty
    assert result.gene_table["selected"].astype(bool).sum() == 2
    assert result.gene_table.loc[result.gene_table["selected"], "gene"].tolist() == kept_genes


def test_run_gedar_respects_direction_filter():
    genes = _genes_in_brain_filter(4)
    weights = pd.DataFrame(
        {
            "gene": genes,
            "weight": [1.0, -2.0, -0.5, 0.2],
        }
    )

    result = run_gedar(
        weights,
        atlas="dk",
        hemisphere="left",
        regions="all",
        direction="down",
        normalize_expression="zscore",
        normalize_weights="none",
    )

    selected = result.gene_table.loc[result.gene_table["selected"]].copy()
    assert result.excluded_table.empty
    assert selected["input_weight"].lt(0).all()
    assert selected["weight"].gt(0).all()
    assert set(selected["gene"]) == {genes[1], genes[2]}

    selection = select_atlas_data(atlas="dk", hemisphere="left", regions="all", zscore_expression=True)
    kept_genes = [genes[1], genes[2]]
    expected_weights = np.abs(np.array([-2.0, -0.5], dtype=float))
    expected = (
        selection.expression.loc[:, kept_genes].to_numpy(dtype=float) @ expected_weights
    ) / expected_weights.sum()
    np.testing.assert_allclose(result.regional_scores["score"].to_numpy(dtype=float), expected)


def test_run_gedar_split_returns_separate_up_and_down_scores():
    genes = _genes_in_brain_filter(4)
    weights = pd.DataFrame(
        {
            "gene": genes,
            "weight": [1.0, -2.0, -0.5, 0.2],
            "p": [0.04, 0.01, 0.03, 0.02],
        }
    )

    result = run_gedar(
        weights,
        atlas="dk",
        hemisphere="left",
        regions="all",
        rank_column="p",
        rank_mode="ascending",
        top_percent=100,
        direction="split",
        normalize_expression="zscore",
        normalize_weights="none",
    )

    assert {"score_up", "score_up_z", "score_down", "score_down_z"}.issubset(result.regional_scores.columns)
    assert {"weight_up", "selected_up", "weight_down", "selected_down"}.issubset(result.gene_table.columns)
    assert result.gene_table["selected_up"].astype(bool).sum() == 2
    assert result.gene_table["selected_down"].astype(bool).sum() == 2
    assert result.gene_table.loc[result.gene_table["selected_up"], "input_weight"].gt(0).all()
    assert result.gene_table.loc[result.gene_table["selected_down"], "input_weight"].lt(0).all()


def test_run_gedar_drops_duplicate_and_invalid_rows_automatically():
    genes = _genes_in_brain_filter(4)
    weights = pd.DataFrame(
        {
            "gene_name": [genes[0], genes[0], genes[1], genes[2], None, genes[3]],
            "z_mean": [0.5, 2.0, np.nan, -1.0, 0.1, 1.5],
            "pvalue": [0.20, 0.01, 0.03, np.inf, 0.05, 0.02],
        }
    )

    result = run_gedar(
        weights,
        atlas="dk",
        hemisphere="left",
        regions="all",
        gene_column="gene_name",
        weight_column="z_mean",
        rank_column="pvalue",
        rank_mode="ascending",
        normalize_expression="zscore",
        normalize_weights="none",
    )

    assert result.gene_table["gene"].tolist() == [genes[0], genes[3]]
    assert result.gene_table["input_weight"].tolist() == [2.0, 1.5]
    assert set(result.excluded_table["exclusion_reason"]) == {"duplicate_symbol", "invalid_weight", "invalid_rank", "missing_gene"}
    assert int((result.excluded_table["exclusion_reason"] == "duplicate_symbol").sum()) == 1
    duplicate_row = result.excluded_table.loc[result.excluded_table["exclusion_reason"] == "duplicate_symbol"].iloc[0]
    assert int(duplicate_row["kept_input_row"]) == 2
    assert tuple(result.requested_genes) == (genes[0], genes[3])


def test_run_gedar_excludes_genes_outside_packaged_brain_filter():
    brain_genes = set(load_brain_gene_symbols())
    all_genes = load_gene_labels("dk").reshape(-1).tolist()
    inside = next(gene for gene in all_genes if gene.upper() in brain_genes)
    outside = next(gene for gene in all_genes if gene.upper() not in brain_genes)

    weights = pd.DataFrame(
        {
            "gene": [inside, outside],
            "weight": [1.0, 5.0],
        }
    )

    result = run_gedar(
        weights,
        atlas="dk",
        hemisphere="left",
        regions="all",
        normalize_expression="zscore",
        normalize_weights="none",
    )

    assert result.gene_table["gene"].tolist() == [inside]
    assert result.requested_genes == (inside,)
    assert set(result.excluded_table["exclusion_reason"]) == {"outside_brain_filter"}
    assert result.excluded_table.loc[0, "input_gene"] == outside
    assert pd.isna(result.excluded_table.loc[0, "kept_input_row"])


def test_run_gedar_writes_expected_outputs(tmp_path: Path):
    genes = _genes_in_brain_filter(5)
    weights = pd.DataFrame(
        {
            "gene": genes,
            "weight": [0.8, -0.3, 0.5, -0.2, 0.1],
            "rank": [0.05, 0.01, 0.02, 0.2, 0.5],
        }
    )

    result = run_gedar(
        weights,
        atlas="dk",
        hemisphere="both",
        regions="all",
        rank_column="rank",
        top_percent=60,
        output_dir=tmp_path,
    )

    assert result.regional_scores.shape[0] == 83
    assert (tmp_path / "gedar_scores.tsv").exists()
    assert (tmp_path / "gedar_genes.tsv").exists()
    assert (tmp_path / "gedar_excluded.tsv").exists()
    assert (tmp_path / "metadata.json").exists()
    assert (tmp_path / "README.txt").exists()
    assert (tmp_path / "matched_genes.txt").exists()
    assert (tmp_path / "missing_genes.txt").exists()
    assert (tmp_path / "plots" / "gedar_scores.png").exists()
    assert (tmp_path / "plots" / "gedar_brain.png").exists()


def test_run_gedar_defaults_to_no_enrichment():
    genes = _genes_in_brain_filter(4)
    weights = pd.DataFrame(
        {
            "gene": genes,
            "weight": [1.0, -2.0, 0.5, 0.2],
        }
    )

    result = run_gedar(
        weights,
        atlas="dk",
        hemisphere="left",
        regions="all",
    )

    assert result.enrichment_method == "none"
    assert result.gsea_table is None
    assert result.ora_tables is None


def test_run_gedar_gsea_uses_full_matched_gene_ranking(monkeypatch):
    genes = _genes_in_brain_filter(4)
    weights = pd.DataFrame(
        {
            "gene": genes,
            "weight": [1.0, -2.0, 0.5, 0.2],
            "rank": [0.01, 0.02, 0.03, 0.50],
        }
    )
    captured: dict[str, object] = {}

    def _fake_gsea(gene_table, *, gene_set, geneset_organism, gene_limit=1500, n_perm=1000):
        captured["genes"] = gene_table["gene"].tolist()
        captured["selected"] = int(gene_table["selected"].astype(bool).sum())
        captured["gene_set"] = gene_set
        captured["geneset_organism"] = geneset_organism
        return pd.DataFrame(
            {
                "Term": ["TEST_TERM"],
                "es": [0.5],
                "nes": [1.2],
                "p_val": [0.03],
                "fdr": [0.05],
            }
        )

    monkeypatch.setattr(gedar_workflow, "_run_gedar_gsea", _fake_gsea)

    result = run_gedar(
        weights,
        atlas="dk",
        hemisphere="left",
        regions="all",
        rank_column="rank",
        top_n=2,
        enrichment_method="gsea",
        gene_set="pooled",
        geneset_organism="Human",
    )

    assert captured["genes"] == genes
    assert captured["selected"] == 2
    assert captured["gene_set"] == "pooled"
    assert result.enrichment_method == "gsea"
    assert result.gsea_table is not None
    assert result.gsea_table.loc[0, "Term"] == "TEST_TERM"


def test_run_gedar_ora_splits_selected_up_and_down_genes(tmp_path: Path):
    genes = _genes_in_brain_filter(5)
    geneset_path = _write_test_gmt(
        tmp_path / "gedar_test.gmt",
        {
            "UP_TERM": [genes[0], genes[2]],
            "DOWN_TERM": [genes[1]],
        },
    )
    weights = pd.DataFrame(
        {
            "gene": genes,
            "weight": [1.0, -2.0, 0.5, -0.2, 0.1],
            "rank": [0.01, 0.02, 0.03, 0.5, 0.6],
        }
    )

    result = run_gedar(
        weights,
        atlas="dk",
        hemisphere="left",
        regions="all",
        rank_column="rank",
        top_n=3,
        enrichment_method="ora",
        gene_set=str(geneset_path),
        output_dir=tmp_path,
    )

    assert result.enrichment_method == "ora"
    assert result.ora_tables is not None
    assert "UP_TERM" in result.ora_tables["up"]["Term"].tolist()
    assert "DOWN_TERM" in result.ora_tables["down"]["Term"].tolist()
    assert (tmp_path / "ora_gedar_up.tsv").exists()
    assert (tmp_path / "ora_gedar_down.tsv").exists()
    assert (tmp_path / "plots" / "ora_gedar_heatmap.png").exists()
