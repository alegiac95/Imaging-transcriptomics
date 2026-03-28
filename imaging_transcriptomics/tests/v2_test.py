from pathlib import Path
import sys
import types

import nibabel as nib
import numpy as np
import pandas as pd
import pytest
import warnings

from imaging_transcriptomics import (
    RunConfig,
    atlas_table,
    build_run_config,
    extract_scan_data,
    run_analysis,
    run_corr,
    run_pls,
    select_atlas_data,
)
import imaging_transcriptomics.api as api
from imaging_transcriptomics._compat import suppress_pkg_resources_deprecation
from imaging_transcriptomics.atlas_registry import get_atlas
from imaging_transcriptomics.corr import CorrAnalysis
from imaging_transcriptomics.genes import PLSGenes
from imaging_transcriptomics.gsea_utils import (
    gsea_style_fdr,
    make_prerank_table,
    normalize_enrichment_nulls,
    normalize_enrichment_scores,
    run_prerank,
)
from imaging_transcriptomics.genesets import get_geneset
import imaging_transcriptomics.nulls as spatial_nulls
import imaging_transcriptomics.pls as pls_module
from imaging_transcriptomics.surfaces import infer_surface_density


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


def test_select_atlas_data_keeps_glasser_expression_finite():
    selection = select_atlas_data(atlas="glasser-360", hemisphere="both", regions="all")
    matrix = selection.expression.iloc[:, 2:].to_numpy(dtype=float)

    assert matrix.shape[0] == 360
    assert np.isfinite(matrix).all()


def test_surface_density_is_inferred_from_packaged_parcellations():
    assert infer_surface_density(get_atlas("dk"), "left") == "10k"
    assert infer_surface_density(get_atlas("schaefer-100"), "left") == "10k"
    assert infer_surface_density(get_atlas("schaefer-200"), "left") == "164k"
    assert infer_surface_density(get_atlas("schaefer-400"), "left") == "164k"
    assert infer_surface_density(get_atlas("glasser-360"), "left") == "32k"


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
        ora_p_threshold=1.0,
    )

    assert result.gene_table.shape[0] == select_atlas_data(atlas="dk", hemisphere="left", regions="all").gene_labels.shape[0]
    assert result.metadata.null_method == "moran"
    assert list(result.gene_table.columns) == ["gene", "score", "p", "fdr", "maxT"]
    assert (tmp_path / "README.txt").exists()
    assert (tmp_path / "metadata.json").exists()
    assert (tmp_path / "regional_values.tsv").exists()
    assert (tmp_path / "corr_genes.tsv").exists()
    assert (tmp_path / "ora_corr_up.tsv").exists()
    assert (tmp_path / "ora_corr_down.tsv").exists()
    assert (tmp_path / "plots" / "regional_values_brain.png").exists()
    assert (tmp_path / "plots" / "regional_values_cortex.png").exists()
    assert (tmp_path / "plots" / "corr_top_genes.png").exists()
    assert (tmp_path / "plots" / "corr_distribution.png").exists()
    assert (tmp_path / "plots" / "ora_corr_heatmap.png").exists()


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


def test_run_pls_reuses_original_and_permuted_fits_once(monkeypatch):
    fit_calls = []

    def fake_permutations(extracted, n_permutations, *, null_method="auto", seed=1234):
        del seed
        values = extracted.values - np.mean(extracted.values)
        std = np.std(values, ddof=1)
        zvalues = values / std if std else values
        return np.tile(zvalues.reshape(-1, 1), (1, n_permutations)), null_method

    def fake_fit(gene_exp, imaging_data, n_components: int, **kwargs):
        del kwargs
        fit_calls.append(int(n_components))
        n_samples = np.asarray(imaging_data, dtype=float).shape[0]
        n_genes = gene_exp.X.shape[1] if hasattr(gene_exp, "X") else np.asarray(gene_exp, dtype=float).shape[1]
        x_scores = np.tile(np.linspace(-1.0, 1.0, n_samples, dtype=float).reshape(-1, 1), (1, n_components))
        x_weights = np.tile(np.linspace(-1.0, 1.0, n_genes, dtype=float).reshape(-1, 1), (1, n_components))
        varexp = np.linspace(0.2, 0.05, n_components, dtype=float)
        return {
            "x_scores": x_scores,
            "x_weights": x_weights,
            "varexp": varexp,
        }

    monkeypatch.setattr(api, "_permute_scan_values", fake_permutations)
    monkeypatch.setattr(pls_module.PLSAnalysis, "_fit_pls", staticmethod(fake_fit))

    result = run_pls(
        np.linspace(-1.0, 1.0, 41),
        atlas="dk",
        hemisphere="left",
        regions="all",
        n_components=1,
        n_permutations=4,
        run_gsea=False,
        n_jobs=2,
    )

    assert result.metadata.method == "pls"
    assert len(result.components) == 1
    assert list(result.components[0].gene_table.columns) == ["gene", "weight", "zscore", "p", "fdr", "maxT"]
    assert fit_calls.count(1) == 4
    assert len(fit_calls) == 5


def test_extract_scan_data_rejects_native_like_nifti_with_clear_error(tmp_path: Path):
    image = nib.Nifti1Image(np.zeros((32, 40, 36), dtype=np.float32), np.eye(4))
    path = tmp_path / "native_like_T1w.nii.gz"
    nib.save(image, path)

    with pytest.raises(ValueError, match="Register the image to MNI152 first"):
        extract_scan_data(path, atlas="dk", hemisphere="left", regions="all")


def test_extract_scan_data_resamples_explicit_mni_volume(tmp_path: Path):
    affine = np.array(
        [
            [-3.0, 0.0, 0.0, 90.0],
            [0.0, 3.0, 0.0, -126.0],
            [0.0, 0.0, 3.0, -72.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        dtype=float,
    )
    image = nib.Nifti1Image(np.ones((61, 73, 61), dtype=np.float32), affine)
    path = tmp_path / "mni_like_3mm.nii.gz"
    nib.save(image, path)

    extracted = extract_scan_data(
        path,
        atlas="dk",
        hemisphere="left",
        regions="all",
        source_space="MNI152",
    )

    assert extracted.values.shape == (41,)
    assert np.isfinite(extracted.values).all()
    assert np.allclose(extracted.values, 1.0, atol=1e-3)


def test_extract_scan_data_accepts_multidot_nifti_filename(tmp_path: Path):
    affine = np.array(
        [
            [-3.0, 0.0, 0.0, 90.0],
            [0.0, 3.0, 0.0, -126.0],
            [0.0, 0.0, 3.0, -72.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        dtype=float,
    )
    image = nib.Nifti1Image(np.ones((61, 73, 61), dtype=np.float32), affine)
    path = tmp_path / "5-HTT.mean.bmax.mrtm2.nopvc.MNI152.sm5.nii.gz"
    nib.save(image, path)

    extracted = extract_scan_data(
        path,
        atlas="dk",
        hemisphere="left",
        regions="all",
        source_space="MNI152",
    )

    assert extracted.values.shape == (41,)
    assert np.isfinite(extracted.values).all()


def test_corr_gsea_writes_results_with_outdir(tmp_path: Path, monkeypatch):
    seen_duplicate_flags = []
    es_sequence = [
        [0.4, -0.6],
        [0.2, -0.3],
        [0.4, -0.9],
    ]

    class _FakeResult:
        def __init__(self, es):
            self.res2d = pd.DataFrame(
                {
                    "Term": ["TermA", "TermB"],
                    "ES": np.asarray(es, dtype=float),
                }
            )

    def _fake_prerank(rnk, *args, **kwargs):
        frame = rnk if isinstance(rnk, pd.DataFrame) else pd.DataFrame(rnk)
        seen_duplicate_flags.append(frame.iloc[:, 1].duplicated().any())
        return _FakeResult(es_sequence[len(seen_duplicate_flags) - 1])

    fake_gseapy = types.SimpleNamespace(prerank=_fake_prerank)
    monkeypatch.setitem(sys.modules, "gseapy", fake_gseapy)

    analysis = CorrAnalysis(n_iterations=2, n_genes=3)
    results = analysis.gene_results.results
    results.genes = np.array([["G1"], ["G2"], ["G3"]], dtype=object)
    results.corr = np.array([[0.4, -0.1, 0.2]], dtype=float)
    results.boot_corr = np.array(
        [
            [0.5, 0.3],
            [-0.2, -0.1],
            [0.1, 0.4],
        ],
        dtype=float,
    )

    analysis.gsea(gene_set="lake", outdir=tmp_path, n_perm=2)

    output = tmp_path / "gsea_corr_results.tsv"
    assert output.exists()
    table = pd.read_csv(output, sep="\t")
    assert list(table.columns) == ["Term", "es", "nes", "p_val", "fdr"]
    assert table.shape == (2, 5)
    expected_nes = normalize_enrichment_scores(
        np.array(es_sequence[0], dtype=float),
        np.array(es_sequence[1:], dtype=float).T,
    )
    expected_fdr = gsea_style_fdr(
        expected_nes,
        normalize_enrichment_nulls(
            np.array(es_sequence[0], dtype=float),
            np.array(es_sequence[1:], dtype=float).T,
        ),
    )
    np.testing.assert_allclose(table["nes"].to_numpy(dtype=float), expected_nes)
    np.testing.assert_allclose(table["p_val"].to_numpy(dtype=float), np.array([2 / 3, 2 / 3]))
    np.testing.assert_allclose(table["fdr"].to_numpy(dtype=float), expected_fdr)
    assert seen_duplicate_flags == [False, False, False]


def test_make_prerank_table_breaks_ties_deterministically():
    table = make_prerank_table(
        ["G1", "G2", "G3", "G4", "G5"],
        [0.5, 0.5, -0.2, -0.2, -0.2],
    )

    assert list(table["gene"]) == ["G1", "G2", "G3", "G4", "G5"]
    assert not table["score"].duplicated().any()
    assert table.loc[0, "score"] < table.loc[1, "score"]
    assert table.loc[2, "score"] < table.loc[3, "score"] < table.loc[4, "score"]


def test_normalize_enrichment_scores_matches_gseapy_same_sign_rule():
    es = np.array([0.6, -0.8], dtype=float)
    esnull = np.array(
        [
            [0.3, 0.6, -0.2],
            [-0.4, -0.8, 0.5],
        ],
        dtype=float,
    )

    nes = normalize_enrichment_scores(es, esnull)

    np.testing.assert_allclose(nes, np.array([4 / 3, -4 / 3], dtype=float))


def test_gsea_style_fdr_matches_expected_tail_ratio():
    nes = np.array([2.0, 1.5, -1.0, -2.0], dtype=float)
    nesnull = np.array(
        [
            [2.5, 1.2],
            [1.6, -0.4],
            [-1.5, 0.2],
            [-2.5, -1.2],
        ],
        dtype=float,
    )

    fdr = gsea_style_fdr(nes, nesnull)

    np.testing.assert_allclose(fdr, np.array([0.5, 0.5, 0.75, 0.5], dtype=float))


def test_pls_gsea_uses_external_nulls_for_nes(tmp_path: Path, monkeypatch):
    es_sequence = [
        [0.6, -0.5],
        [0.3, -0.25],
        [0.9, -0.75],
    ]
    permutation_nums = []

    class _FakeResult:
        def __init__(self, es):
            self.res2d = pd.DataFrame(
                {
                    "es": np.asarray(es, dtype=float),
                    "nes": np.array([99.0, -99.0], dtype=float),
                    "pval": np.array([0.01, 0.02], dtype=float),
                    "fdr": np.array([0.05, 0.06], dtype=float),
                    "geneset_size": np.array([30, 40], dtype=int),
                    "matched_size": np.array([5, 6], dtype=int),
                    "matched_genes": np.array(["G1;G2", "G2;G3"], dtype=object),
                    "ledge_genes": np.array(["G1", "G2"], dtype=object),
                },
                index=["TermA", "TermB"],
            )

    def _fake_prerank(rnk, *args, **kwargs):
        permutation_nums.append(kwargs.get("permutation_num"))
        return _FakeResult(es_sequence[len(permutation_nums) - 1])

    fake_gseapy = types.SimpleNamespace(prerank=_fake_prerank)
    monkeypatch.setitem(sys.modules, "gseapy", fake_gseapy)

    pls = PLSGenes(1, n_iter=2, n_genes=3)
    pls.orig.genes[0, :] = np.array(["G1", "G2", "G3"], dtype=object)
    pls.orig.zscored[0, :] = np.array([1.0, 0.2, -0.5], dtype=float)
    pls.boot.weights[0, :, 0] = np.array([1.2, 0.1, -0.6], dtype=float)
    pls.boot.weights[0, :, 1] = np.array([0.8, 0.3, -0.7], dtype=float)

    pls.gsea(gene_set="lake", outdir=tmp_path, n_iter=2)

    output = tmp_path / "gsea_pls1_results.tsv"
    assert output.exists()
    table = pd.read_csv(output, sep="\t")
    expected_nes = normalize_enrichment_scores(
        np.array(es_sequence[0], dtype=float),
        np.array(es_sequence[1:], dtype=float).T,
    )
    expected_fdr = gsea_style_fdr(
        expected_nes,
        normalize_enrichment_nulls(
            np.array(es_sequence[0], dtype=float),
            np.array(es_sequence[1:], dtype=float).T,
        ),
    )
    np.testing.assert_allclose(table["nes"].to_numpy(dtype=float), expected_nes)
    np.testing.assert_allclose(table["fdr"].to_numpy(dtype=float), expected_fdr)
    assert permutation_nums == [0, 0, 0]


def test_run_prerank_suppresses_duplicate_score_warning(tmp_path: Path, monkeypatch, capsys):
    monkeypatch.setenv("MPLCONFIGDIR", str(tmp_path / "mplconfig"))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    gseapy = pytest.importorskip("gseapy")

    raw_rnk = pd.DataFrame(
        {
            "gene": ["SLC1A2", "GPC5", "ADGRV1", "RNF219-AS1", "ZNF98"],
            "score": [0.5, 0.5, -0.2, -0.2, -0.2],
        }
    )

    run_prerank(
        gseapy,
        raw_rnk,
        get_geneset("lake"),
        outdir=None,
        no_plot=True,
        permutation_num=0,
        min_size=1,
        max_size=5000,
        seed=1234,
    )

    captured = capsys.readouterr()
    assert "Duplicated values found in preranked stats" not in captured.err


def test_suppress_pkg_resources_deprecation_filters_only_target_warning():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with suppress_pkg_resources_deprecation():
            warnings.warn(
                "pkg_resources is deprecated as an API. See https://setuptools.pypa.io/en/latest/pkg_resources.html.",
                UserWarning,
            )
            warnings.warn("some other warning", UserWarning)

    messages = [str(item.message) for item in caught]
    assert "some other warning" in messages
    assert not any("pkg_resources is deprecated as an API" in message for message in messages)
