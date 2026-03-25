from pathlib import Path

import numpy as np
import pytest

from imaging_transcriptomics import (
    RunConfig,
    atlas_table,
    build_run_config,
    extract_scan_data,
    run_analysis,
    run_corr,
    select_atlas_data,
)
import imaging_transcriptomics.api as api
import imaging_transcriptomics.nulls as spatial_nulls


def test_packaged_atlas_table_contains_ready_to_run_presets():
    table = atlas_table(packaged_only=True)
    assert set(table["id"]) == {
        "dk",
        "schaefer-100",
        "schaefer-200",
        "schaefer-400",
        "destrieux",
        "glasser-360",
    }


def test_select_atlas_data_supports_left_and_both_hemispheres():
    left = select_atlas_data(atlas="dk", hemisphere="left", regions="all")
    both = select_atlas_data(atlas="dk", hemisphere="both", regions="all")

    assert left.n_regions == 41
    assert both.n_regions == 83
    assert set(left.labels["hemisphere"]) == {"L"}
    assert {"L", "R", "B"}.issubset(set(both.labels["hemisphere"]))


def test_extract_scan_data_accepts_vectors_for_both_hemispheres():
    extracted = extract_scan_data(np.linspace(0.0, 1.0, 83), atlas="dk", hemisphere="both", regions="all")
    assert extracted.values.shape == (83,)
    assert extracted.selection.n_regions == 83
    assert extracted.source_kind == "vector"


def test_run_corr_writes_readme_tables_and_plots(tmp_path, monkeypatch):
    def fake_permutations(extracted, n_permutations, *, null_method="auto", seed=1234):
        del seed
        values = extracted.values - np.mean(extracted.values)
        std = np.std(values, ddof=1)
        zvalues = values / std if std else values
        return np.tile(zvalues.reshape(-1, 1), (1, n_permutations)), null_method

    monkeypatch.setattr(api, "_permute_scan_values", fake_permutations)
    result = run_corr(
        np.linspace(-1.0, 1.0, 41),
        atlas="dk",
        hemisphere="left",
        regions="all",
        n_permutations=8,
        null_method="moran",
        output_dir=tmp_path,
        run_gsea=False,
    )

    assert result.gene_table.shape[0] == select_atlas_data(atlas="dk", hemisphere="left", regions="all").gene_labels.shape[0]
    assert result.metadata.null_method == "moran"
    assert (tmp_path / "README.txt").exists()
    assert (tmp_path / "metadata.json").exists()
    assert (tmp_path / "regional_values.tsv").exists()
    assert (tmp_path / "corr_genes.tsv").exists()
    assert (tmp_path / "plots" / "regional_values.png").exists()
    assert (tmp_path / "plots" / "corr_top_genes.png").exists()
    assert (tmp_path / "plots" / "corr_distribution.png").exists()


def test_auto_null_method_falls_back_to_random(monkeypatch):
    extracted = extract_scan_data(
        np.linspace(-1.0, 1.0, 41),
        atlas="dk",
        hemisphere="left",
        regions="all",
    )

    def boom(*args, **kwargs):
        raise RuntimeError("surface nulls unavailable")

    monkeypatch.setattr(spatial_nulls, "generate_surface_nulls", boom)
    with pytest.warns(RuntimeWarning, match="Falling back"):
        permuted, used_method = api._permute_scan_values(
            extracted,
            n_permutations=6,
            null_method="auto",
        )

    assert permuted.shape == (41, 6)
    assert used_method == "random"


def test_explicit_null_method_raises_when_surface_nulls_fail(monkeypatch):
    extracted = extract_scan_data(
        np.linspace(-1.0, 1.0, 41),
        atlas="dk",
        hemisphere="left",
        regions="all",
    )

    def boom(*args, **kwargs):
        raise RuntimeError("surface nulls unavailable")

    monkeypatch.setattr(spatial_nulls, "generate_surface_nulls", boom)
    with pytest.raises(RuntimeError, match="Unable to generate cortical nulls"):
        api._permute_scan_values(
            extracted,
            n_permutations=6,
            null_method="vasa",
        )


def test_build_run_config_normalizes_inputs(tmp_path):
    config = build_run_config(
        "corr",
        atlas="DK",
        hemisphere="left",
        regions="all",
        output_dir=tmp_path,
        n_permutations=8,
    )

    assert isinstance(config, RunConfig)
    assert config.atlas == "dk"
    assert config.output_dir == tmp_path


def test_run_analysis_accepts_explicit_config(monkeypatch):
    def fake_permutations(extracted, n_permutations, *, null_method="auto", seed=1234):
        del seed
        values = extracted.values - np.mean(extracted.values)
        std = np.std(values, ddof=1)
        zvalues = values / std if std else values
        return np.tile(zvalues.reshape(-1, 1), (1, n_permutations)), null_method

    monkeypatch.setattr(api, "_permute_scan_values", fake_permutations)
    config = build_run_config("corr", atlas="dk", hemisphere="left", regions="all", n_permutations=4)
    result = run_analysis(np.linspace(-1.0, 1.0, 41), config)

    assert result.metadata.method == "corr"
    assert result.metadata.n_permutations == 4
