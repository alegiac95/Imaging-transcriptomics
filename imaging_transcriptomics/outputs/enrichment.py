from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .common import BACKGROUND, matplotlib_backend, safe_neglog10, save_figure, shorten_labels, style_axes


def gsea_dot_frame(gsea_table: pd.DataFrame, top_n: int) -> pd.DataFrame:
    """Prepare the ranked subset used by the GSEA dot plot."""

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
        _score=safe_neglog10(ranked["fdr"].to_numpy(dtype=float)),
        _abs_nes=np.abs(ranked["nes"].to_numpy(dtype=float)),
    )
    ranked = ranked.sort_values(["fdr", "_abs_nes"], ascending=[True, False], kind="mergesort").head(top_n)
    return ranked.sort_values("nes", ascending=True, kind="mergesort").reset_index(drop=True)


def ensemble_dot_frame(ensemble_table: pd.DataFrame, top_n: int) -> pd.DataFrame:
    """Prepare the ranked subset used by the ensemble-enrichment dot plot."""

    required = {"Term", "category_score", "z_score", "fdr"}
    missing = sorted(required - set(ensemble_table.columns))
    if missing:
        raise ValueError(f"Ensemble table is missing required columns: {', '.join(missing)}")

    ranked = ensemble_table.copy()
    ranked["fdr"] = pd.to_numeric(ranked["fdr"], errors="coerce")
    ranked["z_score"] = pd.to_numeric(ranked["z_score"], errors="coerce")
    ranked["category_score"] = pd.to_numeric(ranked["category_score"], errors="coerce")
    ranked = ranked.dropna(subset=["Term", "category_score", "z_score", "fdr"])
    if ranked.empty:
        return ranked

    ranked = ranked.assign(
        _score=safe_neglog10(ranked["fdr"].to_numpy(dtype=float)),
        _abs_z=np.abs(ranked["z_score"].to_numpy(dtype=float)),
    )
    ranked = ranked.sort_values(["fdr", "_abs_z"], ascending=[True, False], kind="mergesort").head(top_n)
    return ranked.sort_values("z_score", ascending=True, kind="mergesort").reset_index(drop=True)


def ora_stars(fdr: float) -> str:
    """Return significance stars for one ORA q-value."""

    if not np.isfinite(fdr):
        return ""
    if fdr <= 0.001:
        return "***"
    if fdr <= 0.01:
        return "**"
    if fdr <= 0.05:
        return "*"
    return ""


def format_odds_ratio(value: float) -> str:
    """Format an odds ratio annotation for the ORA heatmap."""

    if np.isnan(value):
        return ""
    if np.isinf(value):
        return "inf"
    return f"{value:.2f}"


def ora_heatmap_frame(
    ora_tables: dict[str, pd.DataFrame],
    *,
    top_n: int = 25,
    significant_fdr: float = 0.05,
) -> tuple[list[str], np.ndarray, np.ndarray] | None:
    """Return the ORA heatmap matrices used by the two-row up/down plot."""

    value_maps: dict[str, dict[str, float]] = {}
    annotation_maps: dict[str, dict[str, str]] = {}
    ranking_rows: list[pd.DataFrame] = []
    for direction in ("up", "down"):
        table = ora_tables.get(direction)
        if table is None or table.empty:
            value_maps[direction] = {}
            annotation_maps[direction] = {}
            continue
        frame = table.copy()
        frame["fdr"] = pd.to_numeric(frame["fdr"], errors="coerce")
        frame["odds_ratio"] = pd.to_numeric(frame["odds_ratio"], errors="coerce")
        frame = frame.dropna(subset=["Term", "fdr", "odds_ratio"])
        if frame.empty:
            value_maps[direction] = {}
            annotation_maps[direction] = {}
            continue
        frame["_plot_value"] = safe_neglog10(frame["fdr"].to_numpy(dtype=float))
        frame["_annotation"] = [
            f"{format_odds_ratio(odds)}{ora_stars(fdr)}"
            for odds, fdr in zip(
                frame["odds_ratio"].to_numpy(dtype=float),
                frame["fdr"].to_numpy(dtype=float),
                strict=False,
            )
        ]
        value_maps[direction] = dict(zip(frame["Term"].astype(str), frame["_plot_value"], strict=False))
        annotation_maps[direction] = dict(zip(frame["Term"].astype(str), frame["_annotation"], strict=False))
        ranking_rows.append(frame[["Term", "fdr", "odds_ratio"]])

    if not ranking_rows:
        return None

    combined = pd.concat(ranking_rows, ignore_index=True)
    ranked_terms = (
        combined.groupby("Term", sort=False)
        .agg(min_fdr=("fdr", "min"), max_ratio=("odds_ratio", "max"))
        .sort_values(["min_fdr", "max_ratio"], ascending=[True, False], kind="mergesort")
    )
    significant_terms = ranked_terms.loc[ranked_terms["min_fdr"] <= float(significant_fdr)]
    selected_terms = significant_terms.head(int(top_n)) if not significant_terms.empty else ranked_terms.head(int(top_n))
    term_order = selected_terms.index.astype(str).tolist()
    if not term_order:
        return None

    matrix = np.full((2, len(term_order)), np.nan, dtype=float)
    annotations = np.full((2, len(term_order)), "", dtype=object)
    for row_index, direction in enumerate(("up", "down")):
        value_mapping = value_maps[direction]
        annotation_mapping = annotation_maps[direction]
        for col_index, term in enumerate(term_order):
            if term in value_mapping:
                matrix[row_index, col_index] = value_mapping[term]
                annotations[row_index, col_index] = annotation_mapping[term]
    return term_order, matrix, annotations


def plot_gsea_dotplot(
    gsea_table: pd.DataFrame,
    output_path: Path,
    *,
    title: str,
    top_n: int = 20,
) -> Path | None:
    """Write a GSEA dot plot for the top-ranked enrichment terms."""

    _, plt = matplotlib_backend()
    ranked = gsea_dot_frame(gsea_table, top_n=top_n)
    if ranked.empty:
        return None

    sizes = 90 + 75 * ranked["_score"].to_numpy(dtype=float)
    colors = ranked["nes"].to_numpy(dtype=float)

    fig_height = max(4.5, 0.35 * ranked.shape[0] + 1.5)
    fig, ax = plt.subplots(figsize=(9, fig_height))
    scatter = ax.scatter(
        ranked["nes"],
        shorten_labels(ranked["Term"]),
        s=sizes,
        c=colors,
        cmap="RdBu_r",
        edgecolor="black",
        linewidth=0.35,
        alpha=0.9,
    )
    ax.axvline(0, color="#64748b", lw=1, alpha=0.6)
    style_axes(ax, title=title, xlabel="Normalized enrichment score", ylabel="Gene set", grid_axis="x")
    colorbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    colorbar.set_label("NES")
    return save_figure(fig, output_path)


def plot_ensemble_dotplot(
    ensemble_table: pd.DataFrame,
    output_path: Path,
    *,
    title: str,
    top_n: int = 20,
) -> Path | None:
    """Write an ensemble-enrichment dot plot for the top-ranked terms."""

    _, plt = matplotlib_backend()
    ranked = ensemble_dot_frame(ensemble_table, top_n=top_n)
    if ranked.empty:
        return None

    sizes = 90 + 75 * ranked["_score"].to_numpy(dtype=float)
    colors = ranked["category_score"].to_numpy(dtype=float)

    fig_height = max(4.5, 0.35 * ranked.shape[0] + 1.5)
    fig, ax = plt.subplots(figsize=(9, fig_height))
    scatter = ax.scatter(
        ranked["z_score"],
        shorten_labels(ranked["Term"]),
        s=sizes,
        c=colors,
        cmap="RdBu_r",
        edgecolor="black",
        linewidth=0.35,
        alpha=0.9,
    )
    ax.axvline(0, color="#64748b", lw=1, alpha=0.6)
    style_axes(ax, title=title, xlabel="Category z-score", ylabel="Gene set", grid_axis="x")
    colorbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    colorbar.set_label("Category score")
    return save_figure(fig, output_path)


def plot_ora_heatmap(
    ora_tables: dict[str, pd.DataFrame],
    output_path: Path,
    *,
    title: str,
) -> Path | None:
    """Write the two-row ORA heatmap used by correlation and PLS outputs."""

    _, plt = matplotlib_backend()
    heatmap = ora_heatmap_frame(ora_tables)
    if heatmap is None:
        return None
    terms, matrix, annotations = heatmap
    matrix = np.asarray(matrix, dtype=float)
    fig_width = max(8.0, 0.45 * len(terms))
    fig, ax = plt.subplots(figsize=(fig_width, 3.0))
    cmap = plt.get_cmap("YlGnBu").copy()
    cmap.set_bad(BACKGROUND)
    vmax = np.nanmax(matrix) if np.any(np.isfinite(matrix)) else 1.0
    image = ax.imshow(matrix, aspect="auto", cmap=cmap, interpolation="nearest", vmin=0, vmax=max(vmax, 1.0))
    style_axes(ax, title=title, xlabel="Gene set", ylabel="Direction", grid_axis="")
    ax.set_yticks([0, 1], labels=["up", "down"])
    ax.set_xticks(np.arange(len(terms)), labels=shorten_labels(terms, max_length=34), rotation=55, ha="right")
    for row_index in range(matrix.shape[0]):
        for col_index in range(matrix.shape[1]):
            label = str(annotations[row_index, col_index])
            if not label:
                continue
            value = matrix[row_index, col_index]
            color = "white" if np.isfinite(value) and image.norm(value) > 0.52 else "#0f172a"
            ax.text(
                col_index,
                row_index,
                label,
                ha="center",
                va="center",
                fontsize=7,
                color=color,
            )
    colorbar = fig.colorbar(image, ax=ax, pad=0.02)
    colorbar.set_label("-log10(FDR)")
    return save_figure(fig, output_path)
