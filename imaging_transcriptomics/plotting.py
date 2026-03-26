from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .models import CorrelationResult, PLSComponentResult, PLSResult


def _save(fig, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return path


def _ranked_gene_frame(gene_table: pd.DataFrame, score_column: str, top_n: int) -> pd.DataFrame:
    top = gene_table.head(top_n).copy()
    bottom = gene_table.tail(top_n).copy().iloc[::-1]
    ranked = pd.concat([top, bottom], ignore_index=True)
    ranked["_color"] = ["#0f766e" if value >= 0 else "#b91c1c" for value in ranked[score_column]]
    return ranked


def _barh_plot(
    labels: pd.Series,
    values: pd.Series,
    colors: list[str],
    *,
    title: str,
    xlabel: str,
    output_path: Path,
) -> Path:
    fig, ax = plt.subplots(figsize=(9, 8))
    ax.barh(labels, values, color=colors)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.grid(axis="x", alpha=0.2)
    return _save(fig, output_path)


def _gsea_dot_frame(gsea_table: pd.DataFrame, top_n: int) -> pd.DataFrame:
    required = {"Term", "nes", "fdr"}
    missing = sorted(required - set(gsea_table.columns))
    if missing:
        raise ValueError(f"GSEA table is missing required columns: {', '.join(missing)}")

    ranked = gsea_table.copy()
    ranked["fdr"] = pd.to_numeric(ranked["fdr"], errors="coerce")
    ranked["nes"] = pd.to_numeric(ranked["nes"], errors="coerce")
    ranked = ranked.dropna(subset=["Term", "nes", "fdr"])
    if ranked.empty:
        return ranked

    ranked = ranked.assign(
        _score=-np.log10(np.clip(ranked["fdr"].to_numpy(dtype=float), 1e-300, 1.0)),
        _abs_nes=np.abs(ranked["nes"].to_numpy(dtype=float)),
    )
    ranked = ranked.sort_values(["fdr", "_abs_nes"], ascending=[True, False], kind="mergesort").head(top_n)
    return ranked.sort_values("nes", ascending=True, kind="mergesort").reset_index(drop=True)


def _ora_heatmap_frame(ora_tables: dict[str, pd.DataFrame]) -> tuple[list[str], np.ndarray] | None:
    frames: dict[str, pd.DataFrame] = {}
    score_maps: dict[str, dict[str, float]] = {}
    ranking_rows: list[pd.DataFrame] = []
    for direction in ("up", "down"):
        table = ora_tables.get(direction)
        if table is None or table.empty:
            frames[direction] = pd.DataFrame(columns=["Term", "fdr", "enrichment_ratio"])
            score_maps[direction] = {}
            continue
        frame = table.copy()
        frame["fdr"] = pd.to_numeric(frame["fdr"], errors="coerce")
        frame["enrichment_ratio"] = pd.to_numeric(frame["enrichment_ratio"], errors="coerce")
        frame = frame.dropna(subset=["Term", "fdr", "enrichment_ratio"])
        if frame.empty:
            frames[direction] = frame
            score_maps[direction] = {}
            continue
        frame["_score"] = -np.log10(np.clip(frame["fdr"].to_numpy(dtype=float), 1e-300, 1.0))
        frames[direction] = frame
        score_maps[direction] = dict(zip(frame["Term"].astype(str), frame["_score"], strict=False))
        ranking_rows.append(frame[["Term", "fdr", "enrichment_ratio"]])

    if not ranking_rows:
        return None

    combined = pd.concat(ranking_rows, ignore_index=True)
    term_order = (
        combined.groupby("Term", sort=False)
        .agg(min_fdr=("fdr", "min"), max_ratio=("enrichment_ratio", "max"))
        .sort_values(["min_fdr", "max_ratio"], ascending=[True, False], kind="mergesort")
        .index.astype(str)
        .tolist()
    )
    if not term_order:
        return None

    matrix = np.full((2, len(term_order)), np.nan, dtype=float)
    for row_index, direction in enumerate(("up", "down")):
        mapping = score_maps[direction]
        for col_index, term in enumerate(term_order):
            if term in mapping:
                matrix[row_index, col_index] = mapping[term]
    return term_order, matrix


def plot_gsea_dotplot(
    gsea_table: pd.DataFrame,
    output_path: Path,
    *,
    title: str,
    top_n: int = 20,
) -> Path | None:
    ranked = _gsea_dot_frame(gsea_table, top_n=top_n)
    if ranked.empty:
        return None

    sizes = 80 + 70 * ranked["_score"].to_numpy(dtype=float)
    colors = ranked["nes"].to_numpy(dtype=float)

    fig_height = max(4.5, 0.35 * ranked.shape[0] + 1.5)
    fig, ax = plt.subplots(figsize=(9, fig_height))
    scatter = ax.scatter(
        ranked["nes"],
        ranked["Term"],
        s=sizes,
        c=colors,
        cmap="RdBu_r",
        edgecolor="black",
        linewidth=0.35,
        alpha=0.9,
    )
    ax.axvline(0, color="#64748b", lw=1, alpha=0.6)
    ax.set_title(title)
    ax.set_xlabel("Normalized enrichment score")
    ax.set_ylabel("Gene set")
    ax.grid(axis="x", alpha=0.2)
    colorbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    colorbar.set_label("NES")
    return _save(fig, output_path)


def plot_ora_heatmap(
    ora_tables: dict[str, pd.DataFrame],
    output_path: Path,
    *,
    title: str,
) -> Path | None:
    heatmap = _ora_heatmap_frame(ora_tables)
    if heatmap is None:
        return None
    terms, matrix = heatmap
    fig_width = max(8.0, 0.45 * len(terms))
    fig, ax = plt.subplots(figsize=(fig_width, 2.8))
    cmap = plt.get_cmap("YlOrRd").copy()
    cmap.set_bad("#f8fafc")
    image = ax.imshow(matrix, aspect="auto", cmap=cmap, interpolation="nearest")
    ax.set_title(title)
    ax.set_xlabel("Gene set")
    ax.set_ylabel("Direction")
    ax.set_yticks([0, 1], labels=["up", "down"])
    ax.set_xticks(np.arange(len(terms)), labels=terms, rotation=60, ha="right")
    colorbar = fig.colorbar(image, ax=ax, pad=0.02)
    colorbar.set_label("-log10(FDR)")
    return _save(fig, output_path)


def plot_regional_values(regional_values: pd.DataFrame, output_dir: Path) -> Path:
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(
        np.arange(regional_values.shape[0]),
        regional_values["value"].to_numpy(),
        color="#1f4e79",
        lw=1.5,
    )
    ax.set_title("Regional values")
    ax.set_xlabel("Region index")
    ax.set_ylabel("Value")
    ax.grid(alpha=0.2)
    return _save(fig, output_dir / "plots" / "regional_values.png")


def plot_correlation_ranking(gene_table: pd.DataFrame, output_dir: Path, top_n: int = 20) -> Path:
    ranked = _ranked_gene_frame(gene_table, "score", top_n)
    return _barh_plot(
        ranked["gene"],
        ranked["score"],
        ranked["_color"].tolist(),
        title="Top correlation genes",
        xlabel="Spearman correlation",
        output_path=output_dir / "plots" / "corr_top_genes.png",
    )


def plot_correlation_distribution(gene_table: pd.DataFrame, output_dir: Path) -> Path:
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(gene_table["score"], bins=60, color="#1f4e79", alpha=0.8)
    ax.set_title("Distribution of gene correlations")
    ax.set_xlabel("Spearman correlation")
    ax.set_ylabel("Count")
    return _save(fig, output_dir / "plots" / "corr_distribution.png")


def plot_pls_variance(result: PLSResult, output_dir: Path) -> list[Path]:
    explained = np.array([component.explained_variance for component in result.components], dtype=float)
    cumulative = np.array(result.cumulative_variance[: explained.shape[0]], dtype=float)

    fig_bar, ax_bar = plt.subplots(figsize=(8, 4))
    ax_bar.bar(np.arange(1, explained.shape[0] + 1), 100 * explained, color="#d97706")
    ax_bar.set_title("PLS component variance")
    ax_bar.set_xlabel("Component")
    ax_bar.set_ylabel("Variance explained (%)")

    fig_line, ax_line = plt.subplots(figsize=(8, 4))
    ax_line.plot(
        np.arange(1, cumulative.shape[0] + 1),
        100 * cumulative,
        marker="o",
        color="#b45309",
    )
    ax_line.set_title("Cumulative PLS variance")
    ax_line.set_xlabel("Component")
    ax_line.set_ylabel("Cumulative variance (%)")
    ax_line.set_ylim(bottom=0)
    ax_line.grid(alpha=0.2)

    return [
        _save(fig_bar, output_dir / "plots" / "pls_variance.png"),
        _save(fig_line, output_dir / "plots" / "pls_cumulative_variance.png"),
    ]


def plot_pls_component(component: PLSComponentResult, output_dir: Path, top_n: int = 20) -> Path:
    ranked = _ranked_gene_frame(component.gene_table, "zscore", top_n)
    return _barh_plot(
        ranked["gene"],
        ranked["zscore"],
        ranked["_color"].tolist(),
        title=f"PLS component {component.index} gene weights",
        xlabel="Z-score",
        output_path=output_dir / "plots" / f"pls_component_{component.index}_genes.png",
    )


def save_result_plots(result: CorrelationResult | PLSResult, output_dir: Path) -> list[Path]:
    paths: list[Path] = [plot_regional_values(result.regional_values, output_dir)]
    if isinstance(result, CorrelationResult):
        paths.append(plot_correlation_ranking(result.gene_table, output_dir))
        paths.append(plot_correlation_distribution(result.gene_table, output_dir))
        if result.gsea_table is not None:
            gsea_path = plot_gsea_dotplot(
                result.gsea_table,
                output_dir / "plots" / "gsea_corr_dotplot.png",
                title="Top GSEA terms",
            )
            if gsea_path is not None:
                paths.append(gsea_path)
        if result.ora_tables is not None:
            ora_path = plot_ora_heatmap(
                result.ora_tables,
                output_dir / "plots" / "ora_corr_heatmap.png",
                title="ORA up/down heatmap",
            )
            if ora_path is not None:
                paths.append(ora_path)
        return paths
    paths.extend(plot_pls_variance(result, output_dir))
    for component in result.components:
        paths.append(plot_pls_component(component, output_dir))
        if component.gsea_table is not None:
            gsea_path = plot_gsea_dotplot(
                component.gsea_table,
                output_dir / "plots" / f"gsea_pls{component.index}_dotplot.png",
                title=f"PLS component {component.index} GSEA",
            )
            if gsea_path is not None:
                paths.append(gsea_path)
        if component.ora_tables is not None:
            ora_path = plot_ora_heatmap(
                component.ora_tables,
                output_dir / "plots" / f"ora_pls{component.index}_heatmap.png",
                title=f"PLS component {component.index} ORA up/down heatmap",
            )
            if ora_path is not None:
                paths.append(ora_path)
    return paths
