from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..models import GEDARResult, GenePCAResult, GeneQueryResult, PLSComponentResult, PLSResult
from .common import (
    ACCENT,
    NEGATIVE,
    NEUTRAL,
    POSITIVE,
    matplotlib_backend,
    regional_colors,
    save_figure,
    shorten_labels,
    style_axes,
)


def ranked_gene_frame(gene_table: pd.DataFrame, score_column: str, top_n: int) -> pd.DataFrame:
    """Return the top and bottom tails of one ranked gene table."""

    top = gene_table.head(top_n).copy()
    bottom = gene_table.tail(top_n).copy().iloc[::-1]
    ranked = pd.concat([top, bottom], ignore_index=True)
    ranked["_color"] = [POSITIVE if value >= 0 else NEGATIVE for value in ranked[score_column]]
    return ranked


def barh_plot(
    labels: pd.Series,
    values: pd.Series,
    colors: list[str],
    *,
    title: str,
    xlabel: str,
    output_path: Path,
) -> Path:
    """Write a lollipop-style horizontal ranking plot."""

    _, plt = matplotlib_backend()
    values_array = pd.Series(values).to_numpy(dtype=float)
    labels_array = shorten_labels(pd.Series(labels))
    y_positions = np.arange(len(labels_array))
    fig_height = max(5.5, 0.32 * len(labels_array) + 1.5)
    fig, ax = plt.subplots(figsize=(9.5, fig_height))
    ax.hlines(y_positions, 0, values_array, color=colors, lw=2.4, alpha=0.85)
    ax.scatter(values_array, y_positions, s=55, color=colors, edgecolor="white", linewidth=0.8, zorder=3)
    ax.axvline(0, color="#64748b", lw=1.1, alpha=0.8)
    ax.set_yticks(y_positions, labels=labels_array)
    style_axes(ax, title=title, xlabel=xlabel, ylabel="Gene", grid_axis="x")
    return save_figure(fig, output_path)


def plot_regional_values(regional_values: pd.DataFrame, output_dir: Path) -> Path:
    """Write the legacy region-index line plot for one regional table."""

    _, plt = matplotlib_backend()
    x = np.arange(regional_values.shape[0])
    y = regional_values["value"].to_numpy(dtype=float)
    colors, _ = regional_colors(regional_values)
    fig, ax = plt.subplots(figsize=(12, 5.2))
    ax.plot(x, y, color="#94a3b8", lw=1.2, alpha=0.9, zorder=1)
    ax.scatter(x, y, c=colors, s=28, edgecolor="white", linewidth=0.5, zorder=2)
    ax.axhline(0, color="#64748b", lw=1, alpha=0.7)
    style_axes(ax, title="Regional values", xlabel="Region index", ylabel="Value", grid_axis="y")
    return save_figure(fig, output_dir / "plots" / "regional_values.png")


def plot_region_profile(
    regional_values: pd.DataFrame,
    *,
    value_column: str,
    title: str,
    ylabel: str,
    output_path: Path,
) -> Path:
    """Write a generic region-index profile plot for one regional value column."""

    _, plt = matplotlib_backend()
    x = np.arange(regional_values.shape[0])
    y = regional_values[value_column].to_numpy(dtype=float)
    colors, _ = regional_colors(regional_values)
    fig, ax = plt.subplots(figsize=(12, 5.2))
    ax.plot(x, y, color="#94a3b8", lw=1.2, alpha=0.9, zorder=1)
    ax.scatter(x, y, c=colors, s=28, edgecolor="white", linewidth=0.5, zorder=2)
    ax.axhline(0, color="#64748b", lw=1, alpha=0.7)
    style_axes(ax, title=title, xlabel="Region index", ylabel=ylabel, grid_axis="y")
    return save_figure(fig, output_path)


def plot_correlation_ranking(gene_table: pd.DataFrame, output_dir: Path, top_n: int = 20) -> Path:
    """Write the top positive and negative correlation genes plot."""

    ranked = ranked_gene_frame(gene_table, "score", top_n)
    return barh_plot(
        ranked["gene"],
        ranked["score"],
        ranked["_color"].tolist(),
        title="Top correlation genes",
        xlabel="Spearman correlation",
        output_path=output_dir / "plots" / "corr_top_genes.png",
    )


def plot_correlation_distribution(gene_table: pd.DataFrame, output_dir: Path) -> Path:
    """Write the histogram of gene-wise correlation scores."""

    _, plt = matplotlib_backend()
    fig, ax = plt.subplots(figsize=(8, 5))
    scores = gene_table["score"].to_numpy(dtype=float)
    ax.hist(scores, bins=60, color=ACCENT, alpha=0.8, edgecolor="white", linewidth=0.4)
    ax.axvline(0, color="#64748b", lw=1, alpha=0.7)
    ax.axvline(np.median(scores), color="#f59e0b", lw=1.3, alpha=0.9, linestyle="--")
    style_axes(ax, title="Distribution of gene correlations", xlabel="Spearman correlation", ylabel="Count", grid_axis="y")
    return save_figure(fig, output_dir / "plots" / "corr_distribution.png")


def plot_gene_query_ranking(result: GeneQueryResult, output_dir: Path) -> Path | None:
    """Write the top positive and negative co-expression tails for one seed gene."""

    if result.coexpressed_genes.empty and result.anticorrelated_genes.empty:
        return None
    positive = result.coexpressed_genes.copy()
    positive["_color"] = [POSITIVE] * positive.shape[0]
    negative = result.anticorrelated_genes.copy()
    negative["_color"] = [NEGATIVE] * negative.shape[0]
    ranked = pd.concat([negative, positive], ignore_index=True)
    return barh_plot(
        ranked["gene"],
        ranked["score"],
        ranked["_color"].tolist(),
        title=f"Top co-expression tails with {result.gene}",
        xlabel="Spearman correlation",
        output_path=output_dir / "plots" / "gene_coexpression_top_genes.png",
    )


def plot_gene_query_distribution(result: GeneQueryResult, output_dir: Path) -> Path:
    """Write the distribution of co-expression scores for one seed gene."""

    _, plt = matplotlib_backend()
    fig, ax = plt.subplots(figsize=(8, 5))
    scores = result.gene_table["score"].to_numpy(dtype=float)
    ax.hist(scores, bins=60, color=ACCENT, alpha=0.8, edgecolor="white", linewidth=0.4)
    ax.axvline(0, color="#64748b", lw=1, alpha=0.7)
    ax.axvline(np.median(scores), color="#f59e0b", lw=1.3, alpha=0.9, linestyle="--")
    style_axes(
        ax,
        title=f"Distribution of co-expression with {result.gene}",
        xlabel="Spearman correlation",
        ylabel="Count",
        grid_axis="y",
    )
    return save_figure(fig, output_dir / "plots" / "gene_coexpression_distribution.png")


def plot_gene_query_matrix(result: GeneQueryResult, output_dir: Path) -> Path:
    """Write a seed-plus-top-genes co-expression heatmap."""

    _, plt = matplotlib_backend()
    matrix = result.coexpression_matrix.copy()
    labels = shorten_labels(matrix.index.astype(str))
    n_genes = matrix.shape[0]
    fig_size = max(4.8, min(10.5, 0.42 * n_genes + 2.8))
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    image = ax.imshow(
        matrix.to_numpy(dtype=float),
        cmap="RdBu_r",
        vmin=-1.0,
        vmax=1.0,
        interpolation="nearest",
    )
    ax.set_xticks(np.arange(n_genes), labels=labels, rotation=45, ha="right")
    ax.set_yticks(np.arange(n_genes), labels=labels)
    ax.set_xticks(np.arange(-0.5, n_genes, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_genes, 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=0.9)
    ax.tick_params(which="minor", bottom=False, left=False)
    style_axes(
        ax,
        title=f"Co-expression matrix around {result.gene}",
        xlabel="Gene",
        ylabel="Gene",
        grid_axis="",
    )
    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    colorbar.set_label("Spearman correlation")
    return save_figure(fig, output_dir / "plots" / "gene_coexpression_matrix.png")


def plot_pls_variance(result: PLSResult, output_dir: Path) -> list[Path]:
    """Write the per-component and cumulative PLS variance summaries."""

    _, plt = matplotlib_backend()
    explained = np.array([component.explained_variance for component in result.components], dtype=float)
    cumulative = np.array(result.cumulative_variance[: explained.shape[0]], dtype=float)
    component_p = np.array([component.p_value for component in result.components], dtype=float)

    positions = np.arange(1, explained.shape[0] + 1)
    fig_bar, ax_bar = plt.subplots(figsize=(8.5, 4.4))
    bars = ax_bar.bar(positions, 100 * explained, color="#d97706", edgecolor="white", linewidth=0.7)
    for bar, p_value in zip(bars, component_p, strict=False):
        ax_bar.text(
            bar.get_x() + bar.get_width() / 2.0,
            bar.get_height() + 0.6,
            f"p={p_value:.3g}",
            ha="center",
            va="bottom",
            fontsize=8,
            color=NEUTRAL,
        )
    style_axes(ax_bar, title="PLS component variance", xlabel="Component", ylabel="Variance explained (%)", grid_axis="y")

    fig_line, ax_line = plt.subplots(figsize=(8.5, 4.4))
    ax_line.plot(
        positions,
        100 * cumulative,
        marker="o",
        color="#b45309",
        lw=2,
    )
    style_axes(ax_line, title="Cumulative PLS variance", xlabel="Component", ylabel="Cumulative variance (%)", grid_axis="y")
    ax_line.set_ylim(bottom=0)

    return [
        save_figure(fig_bar, output_dir / "plots" / "pls_variance.png"),
        save_figure(fig_line, output_dir / "plots" / "pls_cumulative_variance.png"),
    ]


def plot_pls_component(component: PLSComponentResult, output_dir: Path, top_n: int = 20) -> Path:
    """Write the top positive and negative genes for one PLS component."""

    ranked = ranked_gene_frame(component.gene_table, "zscore", top_n)
    return barh_plot(
        ranked["gene"],
        ranked["zscore"],
        ranked["_color"].tolist(),
        title=f"PLS component {component.index} gene weights",
        xlabel="Z-score",
        output_path=output_dir / "plots" / f"pls_component_{component.index}_genes.png",
    )


def plot_gene_pca_variance(result: GenePCAResult, output_dir: Path) -> Path:
    """Write the two-panel gene-PCA variance summary."""

    _, plt = matplotlib_backend()
    components = result.variance_table["component"].to_numpy(dtype=int)
    explained = 100 * result.variance_table["variance_explained"].to_numpy(dtype=float)
    cumulative = 100 * result.variance_table["cumulative_variance"].to_numpy(dtype=float)
    fig, (ax_bar, ax_line) = plt.subplots(1, 2, figsize=(10.2, 4.6), gridspec_kw={"wspace": 0.18})
    fig.patch.set_facecolor("white")

    bar_color = POSITIVE
    line_color = "#475569"
    fill_color = "#cbd5e1"

    bars = ax_bar.bar(
        components,
        explained,
        color=bar_color,
        edgecolor="white",
        linewidth=0.8,
        width=0.62,
        zorder=2,
    )
    style_axes(
        ax_bar,
        title="Explained By Component",
        xlabel="Component",
        ylabel="Variance (%)",
        grid_axis="y",
    )
    ax_bar.set_xticks(components)
    ax_bar.set_xlim(0.5, components.max() + 0.5 if components.size else 1.5)
    ax_bar.set_ylim(0, max(float(np.max(explained)) * 1.2, 10.0) if explained.size else 10.0)

    for bar, value in zip(bars, explained, strict=False):
        ax_bar.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + ax_bar.get_ylim()[1] * 0.025,
            f"{value:.1f}%",
            ha="center",
            va="bottom",
            fontsize=9,
            color=NEUTRAL,
        )

    ax_line.fill_between(components, cumulative, 0, color=fill_color, alpha=0.45, zorder=1)
    ax_line.plot(
        components,
        cumulative,
        color=line_color,
        marker="o",
        markersize=7,
        lw=2.2,
        markerfacecolor="white",
        markeredgewidth=2.0,
        zorder=3,
    )
    style_axes(
        ax_line,
        title="Cumulative Variance",
        xlabel="Component",
        ylabel="Variance (%)",
        grid_axis="y",
    )
    ax_line.set_xticks(components)
    ax_line.set_xlim(0.5, components.max() + 0.5 if components.size else 1.5)
    ax_line.set_ylim(0, 100)

    for component, value in zip(components, cumulative, strict=False):
        ax_line.text(
            component,
            min(value + 3.0, 98.0),
            f"{value:.1f}%",
            ha="center",
            va="bottom",
            fontsize=9,
            color=line_color,
        )

    fig.suptitle("Gene PCA Variance Summary", x=0.06, y=0.98, ha="left", fontweight="bold", fontsize=12)
    return save_figure(fig, output_dir / "plots" / "gene_pca_variance.png")


def plot_gene_pca_regional_component(result: GenePCAResult, output_dir: Path, component_index: int) -> Path:
    """Write the region-index plot for one gene-PCA component."""

    _, plt = matplotlib_backend()
    column = f"PC{component_index}"
    x = np.arange(result.regional_scores.shape[0])
    y = result.regional_scores[column].to_numpy(dtype=float)
    colors, _ = regional_colors(result.regional_scores)
    fig, ax = plt.subplots(figsize=(10.5, 4.2))
    ax.plot(x, y, color="#94a3b8", lw=1.2, alpha=0.9)
    ax.scatter(x, y, c=colors, s=26, edgecolor="white", linewidth=0.5, zorder=2)
    ax.axhline(0, color="#64748b", lw=1, alpha=0.6)
    style_axes(ax, title=f"Regional expression pattern for {column}", xlabel="Region index", ylabel=f"{column} score", grid_axis="y")
    return save_figure(fig, output_dir / "plots" / f"gene_pca_{column.lower()}_regions.png")


def plot_gene_pca_loadings(result: GenePCAResult, output_dir: Path, component_index: int, top_n: int = 15) -> Path:
    """Write the top positive and negative loadings for one gene-PCA component."""

    column = f"PC{component_index}"
    frame = result.gene_loadings.sort_values(column, ascending=False, kind="mergesort")
    ranked = pd.concat([frame.head(top_n), frame.tail(top_n)], ignore_index=True)
    colors = [POSITIVE if value >= 0 else NEGATIVE for value in ranked[column]]
    return barh_plot(
        ranked["gene"],
        ranked[column],
        colors,
        title=f"Top gene loadings for {column}",
        xlabel="Loading",
        output_path=output_dir / "plots" / f"gene_pca_{column.lower()}_loadings.png",
    )


def plot_gedar_regional_scores(
    result: GEDARResult,
    output_dir: Path,
    *,
    value_column: str = "score_z",
    title: str = "GEDAR regional score",
    filename: str = "gedar_scores.png",
) -> Path:
    """Write the region-index GEDAR score plot."""

    _, plt = matplotlib_backend()
    x = np.arange(result.regional_scores.shape[0])
    y = result.regional_scores[value_column].to_numpy(dtype=float)
    colors, _ = regional_colors(result.regional_scores)
    fig, ax = plt.subplots(figsize=(10.5, 4.2))
    ax.plot(x, y, color="#94a3b8", lw=1.2, alpha=0.9)
    ax.scatter(x, y, c=colors, s=26, edgecolor="white", linewidth=0.5, zorder=2)
    ax.axhline(0, color="#64748b", lw=1, alpha=0.6)
    style_axes(ax, title=title, xlabel="Region index", ylabel="Z-scored score", grid_axis="y")
    return save_figure(fig, output_dir / "plots" / filename)


def plot_gedar_weights(
    result: GEDARResult,
    output_dir: Path,
    *,
    weight_column: str = "weight",
    selected_column: str = "selected",
    title: str = "GEDAR gene weights",
    filename: str = "gedar_weights.png",
    top_n: int = 20,
) -> Path | None:
    """Write the GEDAR gene-weight tail plot for the selected genes."""

    selected = result.gene_table.loc[result.gene_table[selected_column].astype(bool)].copy()
    if selected.empty:
        return None
    ordered = selected.sort_values(weight_column, ascending=False, kind="mergesort")
    ranked = pd.concat([ordered.head(top_n), ordered.tail(top_n)], ignore_index=True)
    colors = [POSITIVE if value >= 0 else NEGATIVE for value in ranked[weight_column].to_numpy(dtype=float)]
    return barh_plot(
        ranked["gene"],
        ranked[weight_column],
        colors,
        title=title,
        xlabel="Weight used in projection",
        output_path=output_dir / "plots" / filename,
    )
