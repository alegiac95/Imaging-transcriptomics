from __future__ import annotations

import os
import subprocess
import tempfile
from functools import lru_cache
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colormaps
from matplotlib.colors import Normalize
from matplotlib.patches import FancyBboxPatch
from PIL import Image

import pySimpleBrainPlot as psbp


REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "docs" / "chapters" / "images"
STATIC_DIR = REPO_ROOT / "docs" / "_static"
MPLCONFIGDIR = Path(os.environ.get("MPLCONFIGDIR", "/tmp/mplconfig"))
ARROW_SVG = STATIC_DIR / "undraw_arrow.svg"

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": [
            "Avenir Next",
            "Avenir",
            "Trebuchet MS",
            "Gill Sans",
            "Helvetica Neue",
            "Arial",
            "DejaVu Sans",
        ],
    }
)


def _brain_values(seed: int, *, scale: float = 1.9) -> np.ndarray:
    rng = np.random.default_rng(seed)
    x = np.linspace(0, 2 * np.pi, 68, endpoint=False)
    signal = (
        0.9 * np.sin(1.7 * x + rng.uniform(0, np.pi))
        + 0.45 * np.cos(3.2 * x + rng.uniform(0, np.pi))
        + 0.18 * rng.normal(size=x.shape[0])
    )
    signal /= np.max(np.abs(signal))
    return signal * scale


def _smooth_profile(values: np.ndarray, window: int = 3) -> np.ndarray:
    if values.size <= 1:
        return values
    window = min(window, values.size)
    if window <= 1:
        return values
    pad_left = window // 2
    pad_right = window - 1 - pad_left
    padded = np.pad(values, (pad_left, pad_right), mode="edge")
    kernel = np.ones(window) / window
    return np.convolve(padded, kernel, mode="valid")


def _blocky_matrix(rows: int, cols: int, seed: int, *, scale: float = 1.0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    row_profile = rng.normal(0, 0.9, rows)
    col_profile = rng.normal(0, 0.9, cols)
    row_profile = _smooth_profile(row_profile)
    col_profile = _smooth_profile(col_profile)
    data = np.outer(row_profile, col_profile) + 0.28 * rng.normal(size=(rows, cols))
    data /= np.max(np.abs(data))
    return data * scale


def _trim_image(img: Image.Image, *, pad: int = 6) -> Image.Image:
    bbox = img.getchannel("A").getbbox()
    if bbox is None:
        return img
    left, top, right, bottom = bbox
    return img.crop(
        (
            max(left - pad, 0),
            max(top - pad, 0),
            min(right + pad, img.width),
            min(bottom + pad, img.height),
        )
    )


def _compose_lateral_views(img: Image.Image) -> np.ndarray:
    width, height = img.size
    half_w = width // 2
    half_h = height // 2

    left_lateral = _trim_image(img.crop((0, 0, half_w, half_h)), pad=4)
    right_lateral = _trim_image(img.crop((0, half_h, half_w, height)), pad=4)

    panel_height = max(left_lateral.height, right_lateral.height)
    gap = 18
    canvas = Image.new(
        "RGBA",
        (left_lateral.width + right_lateral.width + gap, panel_height),
        (255, 255, 255, 0),
    )
    canvas.paste(left_lateral, (0, (panel_height - left_lateral.height) // 2), left_lateral)
    canvas.paste(
        right_lateral,
        (left_lateral.width + gap, (panel_height - right_lateral.height) // 2),
        right_lateral,
    )
    return np.asarray(canvas)


@lru_cache(maxsize=1)
def _load_arrow_image() -> np.ndarray:
    with tempfile.TemporaryDirectory(prefix="imt-workflow-arrow-") as tmp:
        tmpdir = Path(tmp)
        png_path = tmpdir / "arrow.png"
        subprocess.run(["rsvg-convert", str(ARROW_SVG), "-o", str(png_path)], check=True)
        img = Image.open(png_path).convert("RGBA")
        img = _trim_image(img, pad=0)
        return np.asarray(img)


def _render_brain_png(stem: str, values: np.ndarray, *, cmap: str = "RdBu_r") -> np.ndarray:
    MPLCONFIGDIR.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="imt-workflow-brain-") as tmp:
        tmpdir = Path(tmp)
        psbp.plotBrain(
            "aparc",
            values.copy(),
            vmin=-np.max(np.abs(values)),
            vmax=np.max(np.abs(values)),
            save_path=str(tmpdir),
            save_file=stem,
            cm=cmap,
            scaling=0.03,
            viewer=False,
        )
        svg_path = tmpdir / f"{stem}_aparc.svg"
        png_path = tmpdir / f"{stem}_aparc.png"
        subprocess.run(["rsvg-convert", str(svg_path), "-o", str(png_path)], check=True)
        img = Image.open(png_path).convert("RGBA")

        # Remove the built-in colorbar, trim the canvas, and keep only the lateral views.
        img = img.crop((0, 0, int(img.width * 0.67), img.height))
        img = _trim_image(img)
        return _compose_lateral_views(img)


def _setup_canvas() -> plt.Figure:
    fig = plt.figure(figsize=(13.6, 3.7), dpi=220)
    fig.patch.set_alpha(0)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    return fig


def _add_stage_title(fig: plt.Figure, x: float, title: str, subtitle: str) -> None:
    fig.text(
        x,
        0.905,
        title,
        ha="center",
        va="bottom",
        fontsize=13.5,
        fontweight="bold",
        color="#0E2F5A",
    )
    fig.text(
        x,
        0.852,
        subtitle,
        ha="center",
        va="top",
        fontsize=8.9,
        color="#5B6472",
    )


def _add_arrow(fig: plt.Figure, x0: float, x1: float, *, y: float = 0.54, width: float = 0.064) -> None:
    arrow = _load_arrow_image()
    height = width * (arrow.shape[0] / arrow.shape[1])
    left = ((x0 + x1) / 2) - (width / 2)
    ax = fig.add_axes([left, y - height / 2, width, height])
    ax.imshow(arrow)
    ax.axis("off")


def _add_image(fig: plt.Figure, image: np.ndarray, rect: list[float]) -> None:
    ax = fig.add_axes(rect)
    ax.imshow(image)
    ax.axis("off")


def _style_tile_panel(ax: plt.Axes) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_facecolor("none")


def _add_heatmap(
    fig: plt.Figure,
    rect: list[float],
    data: np.ndarray,
    *,
    cmap: str,
    title: str | None = None,
    left_label: str | None = None,
    left_label_x: float = -0.16,
    title_pad: float = 4.0,
    rounded: bool = True,
) -> None:
    ax = fig.add_axes(rect)
    if rounded:
        bg = FancyBboxPatch(
            (0.0, 0.0),
            1.0,
            1.0,
            boxstyle="round,pad=0.02,rounding_size=0.05",
            transform=ax.transAxes,
            facecolor="#F4F8FC",
            edgecolor="#D7E3F0",
            linewidth=0.8,
            zorder=-20,
        )
        ax.add_patch(bg)

    ax.imshow(data, cmap=cmap, aspect="equal", interpolation="nearest", zorder=1)
    ax.set_anchor("C")
    n_rows, n_cols = data.shape
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color=(1, 1, 1, 0.92), linewidth=1.4)
    ax.tick_params(which="minor", bottom=False, left=False)
    _style_tile_panel(ax)
    if title:
        ax.set_title(title, fontsize=8.5, color="#5B6472", pad=title_pad)
    if left_label:
        ax.text(
            left_label_x,
            0.5,
            left_label,
            transform=ax.transAxes,
            rotation=90,
            ha="center",
            va="center",
            fontsize=8.0,
            color="#5B6472",
        )


def _add_gene_column(
    fig: plt.Figure,
    rect: list[float],
    genes: list[str],
    values: list[float] | None = None,
    *,
    neutral: bool = False,
    heading: str | None = None,
) -> None:
    ax = fig.add_axes(rect)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    cmap = colormaps["RdBu_r"]
    norm = Normalize(vmin=-2, vmax=2)
    if heading:
        ax.text(0.5, 0.98, heading, ha="center", va="top", fontsize=8.4, color="#5B6472")
    y_positions = np.linspace(0.84, 0.14, len(genes))
    for idx, (gene, y) in enumerate(zip(genes, y_positions)):
        if gene == "...":
            ax.text(0.5, y, gene, ha="center", va="center", fontsize=13, color="#7E8EA7")
            continue
        if neutral:
            color = "#1565C0"
        else:
            score = values[idx] if values is not None else 0.0
            color = cmap(norm(score))
        ax.text(0.5, y, gene, ha="center", va="center", fontsize=10.7, fontweight="bold", color=color)


def _add_dotplot(fig: plt.Figure, rect: list[float], terms: list[str], scores: list[float], sizes: list[float]) -> None:
    ax = fig.add_axes(rect)
    y = np.arange(len(terms))
    ax.scatter(scores, y, s=np.asarray(sizes) * 18, c=scores, cmap="RdYlBu_r", edgecolors="none")
    ax.set_yticks(y)
    ax.set_yticklabels(terms, fontsize=8.3, color="#3F4B5C")
    ax.set_xticks([])
    ax.invert_yaxis()
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_facecolor("none")


def _add_component_summary(fig: plt.Figure, rect: list[float]) -> None:
    ax = fig.add_axes(rect)
    vals = np.array([0.31, 0.18, 0.09])
    colors = ["#1565C0", "#19B5E8", "#90BE6D"]
    ax.bar(np.arange(3), vals, color=colors, width=0.62)
    ax.set_xticks(np.arange(3))
    ax.set_xticklabels(["PLS1", "PLS2", "PLS3"], fontsize=8.2, color="#3F4B5C")
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_facecolor("none")


def _add_running_sum_panel(fig: plt.Figure, rect: list[float]) -> None:
    ax = fig.add_axes(rect)
    x = np.linspace(0.06, 0.94, 360)
    y = np.interp(
        x,
        [0.06, 0.18, 0.33, 0.46, 0.58, 0.76, 0.94],
        [0.02, 0.05, 0.24, 0.54, 0.16, -0.18, -0.06],
    )
    y = _smooth_profile(y, window=17)
    ax.fill_between(x, 0, y, where=y >= 0, color="#DCE8F4", alpha=0.95, lw=0)
    ax.fill_between(x, 0, y, where=y < 0, color="#EEF4FA", alpha=0.85, lw=0)
    ax.plot(x, y, color="#1565C0", lw=2.25)
    hit_positions = np.array([0.11, 0.18, 0.29, 0.35, 0.42, 0.51, 0.73, 0.86])
    rank_strip = np.linspace(1, -1, 300).reshape(1, -1)
    ax.imshow(
        rank_strip,
        extent=(0.08, 0.92, -0.70, -0.57),
        cmap="RdBu_r",
        aspect="auto",
        interpolation="nearest",
        alpha=0.9,
        zorder=0,
    )
    for xpos in hit_positions:
        ax.vlines(xpos, -0.79, -0.57, color="#5B6472", lw=1.15, alpha=0.88)
    ax.axhline(0, color="#D7E3F0", lw=1.0)
    ax.text(0.03, 0.95, "running ES", transform=ax.transAxes, ha="left", va="top", fontsize=8.2, color="#5B6472")
    ax.text(0.08, -0.52, "pathway hits", ha="left", va="center", fontsize=7.5, color="#5B6472")
    ax.text(0.92, -0.74, "ranked genes", ha="right", va="center", fontsize=7.5, color="#5B6472")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_facecolor("none")


def _add_ranked_hits_panel(fig: plt.Figure, rect: list[float]) -> None:
    ax = fig.add_axes(rect)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.22, 0.96, "up", ha="center", va="top", fontsize=8.3, color="#5B6472")
    ax.text(0.78, 0.96, "down", ha="center", va="top", fontsize=8.3, color="#5B6472")
    up_genes = ["RELN", "GAD1", "PVALB", "..."]
    down_genes = ["MBP", "GFAP", "MOG", "..."]
    up_vals = [1.8, 1.2, 0.8, 0.0]
    down_vals = [-1.7, -1.1, -0.8, 0.0]
    cmap = colormaps["RdBu_r"]
    norm = Normalize(vmin=-2, vmax=2)
    y_positions = [0.78, 0.58, 0.38, 0.19]
    for gene, score, y in zip(up_genes, up_vals, y_positions):
        color = "#7E8EA7" if gene == "..." else cmap(norm(score))
        ax.text(0.22, y, gene, ha="center", va="center", fontsize=10.1 if gene != "..." else 13, fontweight="bold" if gene != "..." else None, color=color)
    for gene, score, y in zip(down_genes, down_vals, y_positions):
        color = "#7E8EA7" if gene == "..." else cmap(norm(score))
        ax.text(0.78, y, gene, ha="center", va="center", fontsize=10.1 if gene != "..." else 13, fontweight="bold" if gene != "..." else None, color=color)


def _add_term_score_bars(fig: plt.Figure, rect: list[float], *, title: str = "category score") -> None:
    ax = fig.add_axes(rect)
    terms = ["Synapse", "Interneuron", "Myelin", "Stress"]
    scores = np.array([1.9, 1.15, -0.85, -1.35])
    colors = ["#1565C0" if value >= 0 else "#E67E5F" for value in scores]
    y = np.arange(len(terms))
    ax.barh(y, scores, color=colors, height=0.56)
    ax.axvline(0, color="#D7E3F0", lw=1.0)
    ax.set_yticks(y)
    ax.set_yticklabels(terms, fontsize=8.1, color="#3F4B5C")
    ax.set_xticks([])
    ax.invert_yaxis()
    ax.set_title(title, fontsize=8.5, color="#5B6472", pad=4.0)
    ax.tick_params(axis="y", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_facecolor("none")


def _add_weighted_strip(
    fig: plt.Figure,
    rect: list[float],
    *,
    title: str | None = None,
    left_label: str | None = None,
    left_label_x: float = -0.16,
    title_pad: float = 4.0,
) -> None:
    vals = np.linspace(-1.6, 1.6, 12).reshape(-1, 1)
    _add_heatmap(
        fig,
        rect,
        vals,
        cmap="RdBu_r",
        title=title,
        left_label=left_label,
        left_label_x=left_label_x,
        title_pad=title_pad,
    )


def _add_null_panel(fig: plt.Figure, rect: list[float]) -> None:
    ax = fig.add_axes(rect)
    x = np.linspace(-2.6, 2.6, 300)
    specs = [
        {"sigma": 0.7, "obs": 1.25, "color": "#C7254E", "offset": 2.1},
        {"sigma": 0.9, "obs": 0.55, "color": "#E67E5F", "offset": 1.1},
        {"sigma": 0.8, "obs": -1.15, "color": "#1565C0", "offset": 0.1},
    ]
    for spec in specs:
        density = np.exp(-0.5 * (x / spec["sigma"]) ** 2)
        density /= density.max()
        baseline = spec["offset"]
        curve = baseline + density * 0.58
        ax.fill_between(x, baseline, curve, color="#DCE8F4", lw=0)
        ax.plot(x, curve, color="#8BAED1", lw=1.6)
        ax.vlines(spec["obs"], baseline, baseline + 0.56, color=spec["color"], lw=2.1)

    ax.text(-2.45, 2.86, "null scores", fontsize=8.4, color="#5B6472", ha="left", va="top")
    ax.text(2.4, 2.76, "obs", fontsize=7.6, color="#5B6472", ha="right", va="center")
    ax.set_xlim(-2.6, 2.6)
    ax.set_ylim(-0.05, 3.05)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_facecolor("none")


def _save(fig: plt.Figure, name: str) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / name, bbox_inches="tight", pad_inches=0.12, transparent=True)
    plt.close(fig)


def render_correlation() -> None:
    fig = _setup_canvas()
    brain = _render_brain_png("corr_input", _brain_values(2))
    _add_stage_title(fig, 0.135, "Imaging map", "Atlas-aligned cortical pattern")
    _add_stage_title(fig, 0.395, "Regional values", "Extract one atlas vector")
    _add_stage_title(fig, 0.655, "Spatial nulls", "Permute the imaging map")
    _add_stage_title(fig, 0.885, "Gene ranking", "Observed scores and empirical p")

    _add_image(fig, brain, [0.035, 0.22, 0.215, 0.48])
    _add_heatmap(
        fig,
        [0.315, 0.24, 0.165, 0.44],
        _blocky_matrix(12, 2, 21, scale=1.2),
        cmap="RdBu_r",
        title="regional vector",
        left_label="regions",
        left_label_x=-0.22,
    )
    _add_null_panel(fig, [0.565, 0.23, 0.145, 0.43])
    _add_gene_column(fig, [0.805, 0.17, 0.15, 0.58], ["RELN", "SLC1A2", "GAD1", "...", "MBP", "SNAP25", "GFAP"], [1.7, 1.1, 0.8, 0.0, -0.6, -1.1, -1.6], heading="top genes")

    _add_arrow(fig, 0.25, 0.315)
    _add_arrow(fig, 0.48, 0.565)
    _add_arrow(fig, 0.715, 0.805)
    _save(fig, "correlation_pipeline_story.png")


def render_pls() -> None:
    fig = _setup_canvas()
    brain = _render_brain_png("pls_input", _brain_values(7))
    _add_stage_title(fig, 0.135, "Imaging map", "Atlas-aligned cortical pattern")
    _add_stage_title(fig, 0.395, "Expression matrix", "Regions by genes")
    _add_stage_title(fig, 0.655, "Latent components", "Find multivariate gene axes")
    _add_stage_title(fig, 0.885, "Component genes", "Inspect aligned weights")

    _add_image(fig, brain, [0.035, 0.22, 0.215, 0.48])
    _add_heatmap(fig, [0.31, 0.21, 0.175, 0.49], _blocky_matrix(12, 10, 7), cmap="viridis", title="expression matrix", left_label="regions")
    _add_heatmap(
        fig,
        [0.575, 0.285, 0.09, 0.33],
        _blocky_matrix(6, 3, 31, scale=1.2),
        cmap="RdBu_r",
        title="gene weights",
        left_label="genes",
    )
    _add_component_summary(fig, [0.655, 0.30, 0.11, 0.30])
    _add_gene_column(fig, [0.81, 0.17, 0.145, 0.58], ["CAMK2A", "RELN", "GRIN2B", "...", "MBP", "GFAP", "PDYN"], [1.6, 1.1, 0.7, 0.0, -0.5, -1.0, -1.5], heading="PLS1")

    _add_arrow(fig, 0.25, 0.31)
    _add_arrow(fig, 0.49, 0.575)
    _add_arrow(fig, 0.765, 0.81)
    _save(fig, "pls_pipeline_story.png")


def render_gedar() -> None:
    fig = _setup_canvas()
    brain = _render_brain_png("gedar_output", _brain_values(11))
    _add_stage_title(fig, 0.12, "Weighted genes", "TWAS-style effect table")
    _add_stage_title(fig, 0.37, "Atlas matching", "Keep brain-expressed hits")
    _add_stage_title(fig, 0.62, "Weighted score", "Project the signature regionally")
    _add_stage_title(fig, 0.87, "Regional map", "Visualize the GEDAR pattern")

    _add_gene_column(fig, [0.04, 0.17, 0.145, 0.58], ["CACNA1C", "GRIA1", "RELN", "...", "MBP", "GFAP", "SST"], [1.7, 1.1, 0.7, 0.0, -0.7, -1.1, -1.5], heading="z or beta")
    _add_heatmap(fig, [0.285, 0.21, 0.175, 0.49], _blocky_matrix(12, 6, 41), cmap="viridis", title="matched expression", left_label="regions")
    _add_weighted_strip(
        fig,
        [0.565, 0.26, 0.072, 0.42],
        title="weighted score",
        left_label="regions",
        title_pad=14.0,
    )
    _add_gene_column(fig, [0.635, 0.26, 0.085, 0.42], ["R1", "R2", "R3", "...", "R10", "R11", "R12"], neutral=True, heading="regions")
    _add_image(fig, brain, [0.765, 0.22, 0.205, 0.48])

    _add_arrow(fig, 0.195, 0.285)
    _add_arrow(fig, 0.46, 0.565)
    _add_arrow(fig, 0.72, 0.765)
    _save(fig, "gedar_pipeline_story.png")


def render_gene_pca() -> None:
    fig = _setup_canvas()
    brain = _render_brain_png("gene_pca_output", _brain_values(15))
    _add_stage_title(fig, 0.12, "Gene list", "Curated markers or pathways")
    _add_stage_title(fig, 0.37, "Atlas matrix", "Regional expression of retained genes")
    _add_stage_title(fig, 0.62, "PCA summary", "Components and variance")
    _add_stage_title(fig, 0.87, "Regional pattern", "Project component scores back to cortex")

    _add_gene_column(fig, [0.045, 0.17, 0.145, 0.58], ["RELN", "GAD1", "SLC1A2", "...", "VIP", "PVALB", "SST"], neutral=True, heading="selected genes")
    _add_heatmap(fig, [0.285, 0.21, 0.175, 0.49], _blocky_matrix(12, 7, 51), cmap="viridis", title="expression matrix", left_label="regions")
    _add_heatmap(fig, [0.555, 0.30, 0.078, 0.28], _blocky_matrix(5, 3, 61, scale=1.2), cmap="RdBu_r", title="loadings", left_label="genes")
    _add_component_summary(fig, [0.632, 0.30, 0.102, 0.30])
    _add_image(fig, brain, [0.758, 0.22, 0.212, 0.48])

    _add_arrow(fig, 0.195, 0.285)
    _add_arrow(fig, 0.46, 0.555)
    _add_arrow(fig, 0.724, 0.758)
    _save(fig, "gene_pca_pipeline_story.png")


def render_gene() -> None:
    fig = _setup_canvas()
    brain = _render_brain_png("gene_output", _brain_values(23))
    _add_stage_title(fig, 0.12, "Gene query", "Select one atlas gene")
    _add_stage_title(fig, 0.37, "Atlas profile", "Regional expression vector")
    _add_stage_title(fig, 0.64, "Co-expression tails", "Positive and negative partners")
    _add_stage_title(fig, 0.88, "Regional map", "Visualize the queried gene")

    _add_gene_column(fig, [0.05, 0.22, 0.12, 0.48], ["RELN"], neutral=True, heading="seed gene")
    _add_heatmap(
        fig,
        [0.28, 0.24, 0.16, 0.42],
        _blocky_matrix(12, 1, 81, scale=1.5),
        cmap="RdBu_r",
        title="expression vector",
        left_label="regions",
    )
    _add_gene_column(
        fig,
        [0.565, 0.17, 0.155, 0.58],
        ["SLC1A2", "GAD1", "PVALB", "...", "MBP", "GFAP", "PDYN"],
        [1.7, 1.1, 0.7, 0.0, -0.7, -1.1, -1.6],
        heading="top hits",
    )
    _add_heatmap(
        fig,
        [0.69, 0.25, 0.09, 0.38],
        _blocky_matrix(8, 8, 82, scale=1.2),
        cmap="RdBu_r",
        title="matrix",
        left_label="genes",
    )
    _add_image(fig, brain, [0.80, 0.22, 0.19, 0.48])

    _add_arrow(fig, 0.18, 0.28)
    _add_arrow(fig, 0.44, 0.565)
    _add_arrow(fig, 0.78, 0.80)
    _save(fig, "gene_pipeline_story.png")


def render_enrichment() -> None:
    fig = _setup_canvas()
    brain = _render_brain_png("enrichment_input", _brain_values(19))
    _add_stage_title(fig, 0.105, "Regional map", "Output from corr or PLS")
    _add_stage_title(fig, 0.31, "Gene signature", "Continuous scores or hit lists")
    _add_stage_title(fig, 0.54, "GSEA", "Rank-based pathway shifts")
    _add_stage_title(fig, 0.73, "ORA", "Thresholded overlap tests")
    _add_stage_title(fig, 0.90, "Ensemble", "Category scores under null phenotypes")

    _add_image(fig, brain, [0.025, 0.22, 0.16, 0.48])
    _add_gene_column(fig, [0.215, 0.18, 0.12, 0.56], ["RELN", "CAMK2A", "GAD1", "...", "GFAP", "MBP", "PDYN"], [1.8, 1.2, 0.7, 0.0, -0.6, -1.1, -1.7], heading="ranked genes")
    _add_running_sum_panel(fig, [0.45, 0.25, 0.12, 0.36])
    _add_heatmap(fig, [0.66, 0.30, 0.10, 0.28], _blocky_matrix(2, 4, 71, scale=1.2), cmap="RdBu_r", title="up / down")
    _add_term_score_bars(fig, [0.835, 0.24, 0.12, 0.40], title="null-aware terms")

    _add_arrow(fig, 0.185, 0.215, width=0.06)
    _add_arrow(fig, 0.335, 0.45, width=0.06)
    _add_arrow(fig, 0.57, 0.66, width=0.058)
    _add_arrow(fig, 0.765, 0.835, width=0.058)
    _save(fig, "enrichment_pipeline_story.png")


def render_enrichment_gsea() -> None:
    fig = _setup_canvas()
    brain = _render_brain_png("enrichment_gsea_input", _brain_values(25))
    _add_stage_title(fig, 0.13, "Regional map", "Observed corr or PLS output")
    _add_stage_title(fig, 0.39, "Ranked genes", "Use the whole continuous signature")
    _add_stage_title(fig, 0.66, "Running-sum test", "Track where pathway genes fall")
    _add_stage_title(fig, 0.88, "Pathway table", "ES, NES, p, and q")

    _add_image(fig, brain, [0.035, 0.22, 0.195, 0.48])
    _add_gene_column(fig, [0.285, 0.17, 0.165, 0.58], ["RELN", "GAD1", "SLC1A2", "...", "GFAP", "MBP", "PDYN"], [1.8, 1.2, 0.8, 0.0, -0.6, -1.1, -1.7], heading="full ranking")
    _add_running_sum_panel(fig, [0.57, 0.24, 0.13, 0.38])
    _add_dotplot(fig, [0.82, 0.22, 0.14, 0.48], ["Synapse", "Interneuron", "Glutamate", "Myelin"], [1.7, 1.0, -0.9, -1.4], [5, 4, 3, 2])

    _add_arrow(fig, 0.23, 0.285)
    _add_arrow(fig, 0.45, 0.57)
    _add_arrow(fig, 0.70, 0.82)
    _save(fig, "enrichment_gsea_story.png")


def render_enrichment_ora() -> None:
    fig = _setup_canvas()
    brain = _render_brain_png("enrichment_ora_input", _brain_values(27))
    _add_stage_title(fig, 0.13, "Regional map", "Observed corr or PLS output")
    _add_stage_title(fig, 0.39, "Significant genes", "Threshold positive and negative tails")
    _add_stage_title(fig, 0.66, "Set overlap", "Count hits against each pathway")
    _add_stage_title(fig, 0.88, "ORA summary", "Odds ratios and BH-FDR")

    _add_image(fig, brain, [0.035, 0.22, 0.195, 0.48])
    _add_ranked_hits_panel(fig, [0.28, 0.18, 0.17, 0.56])
    _add_heatmap(fig, [0.56, 0.24, 0.13, 0.38], _blocky_matrix(3, 4, 73, scale=1.2), cmap="RdBu_r", title="overlap counts", left_label="terms")
    _add_heatmap(fig, [0.81, 0.28, 0.13, 0.30], _blocky_matrix(2, 4, 74, scale=1.2), cmap="RdBu_r", title="up / down")

    _add_arrow(fig, 0.23, 0.28)
    _add_arrow(fig, 0.45, 0.56)
    _add_arrow(fig, 0.69, 0.81)
    _save(fig, "enrichment_ora_story.png")


def render_enrichment_ensemble() -> None:
    fig = _setup_canvas()
    brain = _render_brain_png("enrichment_ensemble_input", _brain_values(29))
    _add_stage_title(fig, 0.13, "Regional map", "Observed phenotype of interest")
    _add_stage_title(fig, 0.39, "Gene scores", "Correlate or weight genes once")
    _add_stage_title(fig, 0.66, "Null ensemble", "Recompute category scores on null phenotypes")
    _add_stage_title(fig, 0.88, "Category table", "Empirical p and FDR")

    _add_image(fig, brain, [0.035, 0.22, 0.195, 0.48])
    _add_gene_column(fig, [0.285, 0.17, 0.155, 0.58], ["RELN", "CAMK2A", "GAD1", "...", "GFAP", "MBP", "PDYN"], [1.8, 1.2, 0.7, 0.0, -0.6, -1.1, -1.7], heading="gene scores")
    _add_null_panel(fig, [0.545, 0.23, 0.155, 0.43])
    _add_term_score_bars(fig, [0.81, 0.22, 0.15, 0.42], title="term score")

    _add_arrow(fig, 0.23, 0.285)
    _add_arrow(fig, 0.44, 0.545)
    _add_arrow(fig, 0.70, 0.81)
    _save(fig, "enrichment_ensemble_story.png")


def main() -> None:
    render_correlation()
    render_pls()
    render_gene()
    render_gedar()
    render_gene_pca()
    render_enrichment()
    render_enrichment_gsea()
    render_enrichment_ora()
    render_enrichment_ensemble()


if __name__ == "__main__":
    main()
