from pathlib import Path

import numpy as np
import pandas as pd

from imaging_transcriptomics.models import (
    AnalysisMetadata,
    CorrelationResult,
    GenePCAResult,
    PLSComponentResult,
    PLSResult,
)
from imaging_transcriptomics.plotting import (
    _ora_heatmap_frame,
    _load_surface_parcellation,
    _surface_value_frames,
    _vertex_values_for_hemisphere,
    plot_brain_volume_map,
    plot_cortical_surface_map,
    save_result_plots,
)
from imaging_transcriptomics import select_atlas_data
from imaging_transcriptomics.atlas_registry import get_atlas
from imaging_transcriptomics.outputs import bundle as bundle_output
from imaging_transcriptomics.outputs import brain as brain_output


def _gsea_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Term": ["Astrocytes", "Microglia", "Neurons"],
            "es": [0.3, -0.2, 0.25],
            "nes": [1.8, -1.4, 1.1],
            "p_val": [0.001, 0.01, 0.03],
            "fdr": [0.005, 0.02, 0.04],
        }
    )


def _ora_tables() -> dict[str, pd.DataFrame]:
    return {
        "up": pd.DataFrame(
            {
                "Term": ["Astrocytes", "Neurons"],
                "overlap_size": [8, 5],
                "set_size": [100, 80],
                "selected_size": [30, 30],
                "universe_size": [15000, 15000],
                "enrichment_ratio": [4.0, 3.2],
                "odds_ratio": [4.5, 3.4],
                "odds_ratio_ci_low": [2.1, 1.7],
                "odds_ratio_ci_high": [9.6, 6.8],
                "p_value": [0.001, 0.01],
                "fdr": [0.01, 0.04],
                "overlap_genes": ["A;B", "C;D"],
            }
        ),
        "down": pd.DataFrame(
            {
                "Term": ["Microglia"],
                "overlap_size": [6],
                "set_size": [120],
                "selected_size": [25],
                "universe_size": [15000],
                "enrichment_ratio": [3.0],
                "odds_ratio": [3.2],
                "odds_ratio_ci_low": [1.4],
                "odds_ratio_ci_high": [7.1],
                "p_value": [0.02],
                "fdr": [0.05],
                "overlap_genes": ["E;F"],
            }
        ),
    }


def test_save_result_plots_writes_corr_gsea_dotplot(tmp_path: Path):
    result = CorrelationResult(
        metadata=AnalysisMetadata(
            method="corr",
            atlas_id="dk",
            atlas_label="Desikan-Killiany (83 regions)",
            hemisphere="left",
            regions="all",
            source="array",
            source_kind="vector",
            source_space=None,
            n_permutations=8,
        ),
        regional_values=pd.DataFrame({"id": [1, 2, 3], "label": ["a", "b", "c"], "value": [0.1, -0.2, 0.3]}),
        gene_table=pd.DataFrame(
            {
                "gene": ["A", "B", "C"],
                "score": [0.5, -0.4, 0.2],
                "p_value": [0.01, 0.02, 0.03],
                "fdr": [0.02, 0.03, 0.04],
                "fwer_maxT": [0.05, 0.1, 0.2],
            }
        ),
        gsea_table=_gsea_table(),
        ora_tables=_ora_tables(),
    )

    paths = save_result_plots(result, tmp_path)

    assert tmp_path.joinpath("plots", "regional_values_brain.png").exists()
    assert tmp_path.joinpath("plots", "gsea_corr_dotplot.png").exists()
    assert tmp_path.joinpath("plots", "ora_corr_heatmap.png").exists()
    assert any(path.name == "regional_values_brain.png" for path in paths)
    assert any(path.name == "gsea_corr_dotplot.png" for path in paths)
    assert any(path.name == "ora_corr_heatmap.png" for path in paths)


def test_save_result_plots_writes_pls_gsea_dotplot(tmp_path: Path):
    component = PLSComponentResult(
        index=1,
        explained_variance=0.2,
        p_value=0.01,
        gene_table=pd.DataFrame(
            {
                "gene": ["A", "B", "C"],
                "weight": [0.6, -0.3, 0.1],
                "zscore": [2.0, -1.8, 0.5],
                "p_value": [0.01, 0.02, 0.1],
                "fdr": [0.02, 0.04, 0.2],
            }
        ),
        gsea_table=_gsea_table(),
        ora_tables=_ora_tables(),
    )
    result = PLSResult(
        metadata=AnalysisMetadata(
            method="pls",
            atlas_id="dk",
            atlas_label="Desikan-Killiany (83 regions)",
            hemisphere="left",
            regions="all",
            source="array",
            source_kind="vector",
            source_space=None,
            n_permutations=8,
            n_components=1,
        ),
        regional_values=pd.DataFrame({"id": [1, 2, 3], "label": ["a", "b", "c"], "value": [0.1, -0.2, 0.3]}),
        components=(component,),
        cumulative_variance=np.array([0.2]),
    )

    paths = save_result_plots(result, tmp_path)

    assert tmp_path.joinpath("plots", "regional_values_brain.png").exists()
    assert tmp_path.joinpath("plots", "gsea_pls1_dotplot.png").exists()
    assert tmp_path.joinpath("plots", "ora_pls1_heatmap.png").exists()
    assert any(path.name == "regional_values_brain.png" for path in paths)
    assert any(path.name == "gsea_pls1_dotplot.png" for path in paths)
    assert any(path.name == "ora_pls1_heatmap.png" for path in paths)


def test_save_result_plots_writes_gene_pca_plots(tmp_path: Path):
    result = GenePCAResult(
        atlas_id="dk",
        atlas_label="Desikan-Killiany (83 regions)",
        hemisphere="left",
        regions="all",
        requested_genes=("A", "B", "C"),
        regional_scores=pd.DataFrame(
            {
                "id": [1, 2, 3],
                "label": ["a", "b", "c"],
                "PC1": [0.3, -0.2, 0.1],
                "PC2": [0.1, 0.0, -0.1],
            }
        ),
        gene_loadings=pd.DataFrame(
            {
                "gene": ["A", "B", "C"],
                "PC1": [0.8, -0.4, 0.2],
                "PC2": [0.1, 0.5, -0.3],
            }
        ),
        variance_table=pd.DataFrame(
            {
                "component": [1, 2],
                "variance_explained": [0.6, 0.2],
                "cumulative_variance": [0.6, 0.8],
            }
        ),
        matched_genes=("A", "B", "C"),
        brain_filtered_genes=(),
        missing_genes=(),
    )

    paths = save_result_plots(result, tmp_path)

    assert tmp_path.joinpath("plots", "gene_pca_variance.png").exists()
    assert tmp_path.joinpath("plots", "gene_pca_pc1_brain.png").exists()
    assert tmp_path.joinpath("plots", "gene_pca_pc1_regions.png").exists()
    assert tmp_path.joinpath("plots", "gene_pca_pc1_loadings.png").exists()
    assert tmp_path.joinpath("plots", "gene_pca_pc2_brain.png").exists()
    assert tmp_path.joinpath("plots", "gene_pca_pc2_regions.png").exists()
    assert tmp_path.joinpath("plots", "gene_pca_pc2_loadings.png").exists()
    assert any(path.name == "gene_pca_variance.png" for path in paths)
    assert any(path.name == "gene_pca_pc1_brain.png" for path in paths)


def test_plot_cortical_surface_map_writes_glasser_surface_plot(tmp_path: Path):
    selection = select_atlas_data(atlas="glasser-360", hemisphere="left", regions="default")
    regional = selection.labels.assign(value=np.linspace(-1.0, 1.0, selection.n_regions))

    path = plot_cortical_surface_map(
        regional,
        atlas_id="glasser-360",
        value_column="value",
        title="Glasser cortical map",
        output_path=tmp_path / "glasser_cortex.png",
    )

    assert path is not None
    assert path.exists()


def test_plot_cortical_surface_map_writes_bilateral_surface_plot(tmp_path: Path):
    selection = select_atlas_data(atlas="glasser-360", hemisphere="both", regions="default")
    regional = selection.labels.assign(value=np.linspace(-1.0, 1.0, selection.n_regions))

    path = plot_cortical_surface_map(
        regional,
        atlas_id="glasser-360",
        value_column="value",
        title="Glasser cortical map",
        output_path=tmp_path / "glasser_cortex.png",
    )

    assert path is not None
    assert path.exists()


def test_plot_brain_volume_map_keeps_missing_subcortex_visible(tmp_path: Path):
    selection = select_atlas_data(atlas="dk", hemisphere="both", regions="all")
    cortical = selection.labels.loc[selection.labels["structure"] == "cortex"].copy()
    cortical["score_z"] = np.linspace(-1.0, 1.0, cortical.shape[0])

    path = plot_brain_volume_map(
        cortical,
        atlas_id="dk",
        value_column="score_z",
        title="GEDAR brain map",
        output_path=tmp_path / "gedar_brain.png",
    )

    assert path is not None
    assert path.exists()


def test_surface_mesh_paths_passes_requested_mesh_kind(monkeypatch):
    class FakeAtlas:
        surface_geometry = None
        surface_space = "fsaverage"
        surface_density = "10k"

    seen: list[str] = []

    def fake_fetch(space: str, density: str, mesh_kind: str = "pial"):
        seen.append(mesh_kind)
        return ("left.surf.gii", "right.surf.gii")

    monkeypatch.setattr(brain_output, "get_atlas", lambda atlas_id: FakeAtlas())
    monkeypatch.setattr(brain_output, "fetch_standard_surface_meshes", fake_fetch)

    paths = brain_output.surface_mesh_paths("dk", mesh_kind="inflated")

    assert paths == ("left.surf.gii", "right.surf.gii")
    assert seen == ["inflated"]


def test_save_result_plots_writes_brainspace_comparison_when_available(tmp_path: Path, monkeypatch):
    selection = select_atlas_data(atlas="glasser-360", hemisphere="both", regions="default")
    regional_scores = selection.labels.assign(PC1=np.linspace(-1.0, 1.0, selection.n_regions))
    result = GenePCAResult(
        atlas_id="glasser-360",
        atlas_label="Glasser 360",
        hemisphere="both",
        regions="all",
        requested_genes=("A", "B", "C"),
        regional_scores=regional_scores,
        gene_loadings=pd.DataFrame({"gene": ["A", "B", "C"], "PC1": [0.5, -0.3, 0.2]}),
        variance_table=pd.DataFrame(
            {
                "component": [1],
                "variance_explained": [0.6],
                "cumulative_variance": [0.6],
            }
        ),
        matched_genes=("A", "B", "C"),
        brain_filtered_genes=(),
        missing_genes=(),
    )

    def fake_brainspace(*args, output_path: Path, **kwargs):
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_bytes(b"brainspace")
        return output_path

    monkeypatch.setattr(bundle_output, "plot_cortical_surface_map_brainspace", fake_brainspace)

    paths = save_result_plots(result, tmp_path)

    assert tmp_path.joinpath("plots", "gene_pca_pc1_cortex.png").exists()
    assert tmp_path.joinpath("plots", "gene_pca_pc1_cortex_brainspace.png").exists()
    assert any(path.name == "gene_pca_pc1_cortex_brainspace.png" for path in paths)


def test_ora_heatmap_frame_limits_large_term_sets_to_top_25():
    terms = [f"Term {index}" for index in range(40)]
    frame = pd.DataFrame(
        {
            "Term": terms,
            "odds_ratio": np.linspace(2.0, 5.0, 40),
            "fdr": np.linspace(0.2, 0.9, 40),
        }
    )

    heatmap = _ora_heatmap_frame({"up": frame, "down": pd.DataFrame()})

    assert heatmap is not None
    term_order, matrix, annotations = heatmap
    assert len(term_order) == 25
    assert matrix.shape == (2, 25)
    assert annotations.shape == (2, 25)


def test_ora_heatmap_frame_prefers_significant_terms_only_when_present():
    frame = pd.DataFrame(
        {
            "Term": ["sig_a", "sig_b", "nonsig_a", "nonsig_b"],
            "odds_ratio": [4.0, 3.5, 2.0, 1.8],
            "fdr": [0.001, 0.04, 0.2, 0.7],
        }
    )

    heatmap = _ora_heatmap_frame({"up": frame, "down": pd.DataFrame()})

    assert heatmap is not None
    term_order, _, _ = heatmap
    assert term_order == ["sig_a", "sig_b"]


def test_surface_value_mapping_handles_bilateral_dk_surface_ids():
    selection = select_atlas_data(atlas="dk", hemisphere="both", regions="all")
    regional = selection.labels.assign(value=np.linspace(-1.0, 1.0, selection.n_regions))
    frames = _surface_value_frames(regional)
    atlas = get_atlas("dk")

    label_array, code_to_name = _load_surface_parcellation(str(atlas.surface_paths[1]))
    vertex_values = _vertex_values_for_hemisphere(
        frames["right"],
        value_column="value",
        label_array=label_array,
        code_to_name=code_to_name,
    )

    assert np.isfinite(vertex_values).any()
