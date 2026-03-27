from pathlib import Path

import numpy as np
import pandas as pd

from imaging_transcriptomics.ora import ora_from_gene_table


def test_ora_from_gene_table_splits_up_and_down_and_uses_p_threshold(tmp_path: Path):
    gmt = tmp_path / "toy.gmt"
    gmt.write_text(
        "\n".join(
            [
                "PathUp\tna\tA\tB\tF",
                "PathDown\tna\tC\tD\tE",
                "Mixed\tna\tA\tC\tE",
            ]
        )
        + "\n"
    )

    gene_table = pd.DataFrame(
        {
            "gene": ["A", "B", "C", "D", "E", "F"],
            "score": [1.2, 0.8, -1.1, -0.7, 0.1, -0.2],
            "p_value": [0.01, 0.02, 0.01, 0.03, 0.04, 0.04],
        }
    )

    tables = ora_from_gene_table(
        gene_table,
        gene_set=str(gmt),
        score_column="score",
        p_threshold=0.05,
    )

    up = tables["up"]
    down = tables["down"]

    assert up.iloc[0]["Term"] == "PathUp"
    assert int(up.iloc[0]["overlap_size"]) == 2
    assert int(up.iloc[0]["selected_size"]) == 3
    assert int(up.iloc[0]["universe_size"]) == 6
    assert set(up.iloc[0]["overlap_genes"].split(";")) == {"A", "B"}
    assert float(up.iloc[0]["odds_ratio"]) == 4.0
    assert float(up.iloc[0]["odds_ratio_ci_low"]) < float(up.iloc[0]["odds_ratio"])
    assert float(up.iloc[0]["odds_ratio_ci_high"]) > float(up.iloc[0]["odds_ratio"])

    assert down.iloc[0]["Term"] == "PathDown"
    assert int(down.iloc[0]["overlap_size"]) == 2
    assert int(down.iloc[0]["selected_size"]) == 3
    assert set(down.iloc[0]["overlap_genes"].split(";")) == {"C", "D"}
    assert float(down.iloc[0]["odds_ratio"]) == 4.0
    assert float(down.iloc[0]["odds_ratio_ci_low"]) < float(down.iloc[0]["odds_ratio"])
    assert float(down.iloc[0]["odds_ratio_ci_high"]) > float(down.iloc[0]["odds_ratio"])
