from pathlib import Path

import numpy as np
import pandas as pd

from imaging_transcriptomics.models import (
    AnalysisMetadata,
    CorrelationResult,
    PLSComponentResult,
    PLSResult,
)
from imaging_transcriptomics.plotting import save_result_plots


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
        regional_values=pd.DataFrame({"label": ["a", "b", "c"], "value": [0.1, -0.2, 0.3]}),
        gene_table=pd.DataFrame(
            {
                "gene": ["A", "B", "C"],
                "score": [0.5, -0.4, 0.2],
                "p_value": [0.01, 0.02, 0.03],
                "fdr": [0.02, 0.03, 0.04],
            }
        ),
        gsea_table=_gsea_table(),
    )

    paths = save_result_plots(result, tmp_path)

    assert tmp_path.joinpath("plots", "gsea_corr_dotplot.png").exists()
    assert any(path.name == "gsea_corr_dotplot.png" for path in paths)


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
        regional_values=pd.DataFrame({"label": ["a", "b", "c"], "value": [0.1, -0.2, 0.3]}),
        components=(component,),
        cumulative_variance=np.array([0.2]),
    )

    paths = save_result_plots(result, tmp_path)

    assert tmp_path.joinpath("plots", "gsea_pls1_dotplot.png").exists()
    assert any(path.name == "gsea_pls1_dotplot.png" for path in paths)
