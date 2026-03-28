from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..exceptions import PlottingUnavailableError


POSITIVE = "#0f766e"
NEGATIVE = "#b91c1c"
ACCENT = "#1d4ed8"
NEUTRAL = "#334155"
GRID = "#cbd5e1"
BACKGROUND = "#f8fafc"
MISSING = "#e2e8f0"


def matplotlib_backend():
    """Return the lazy matplotlib imports used by the plotting layer."""

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise PlottingUnavailableError(
            "Plot writing requires the optional plotting dependencies. "
            "Install imaging-transcriptomics[plots]."
        ) from exc
    return matplotlib, plt


def save_figure(fig, path: Path) -> Path:
    """Write one figure to disk and close it consistently."""

    _, plt = matplotlib_backend()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def style_axes(ax, *, title: str, xlabel: str, ylabel: str, grid_axis: str = "y") -> None:
    """Apply the shared plotting style used across result figures."""

    ax.set_title(title, loc="left", fontweight="bold", fontsize=12)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_facecolor("white")
    if grid_axis:
        ax.grid(axis=grid_axis, color=GRID, alpha=0.45, linewidth=0.8)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("#94a3b8")
    ax.spines["bottom"].set_color("#94a3b8")


def safe_neglog10(values: np.ndarray) -> np.ndarray:
    """Convert p-like values to ``-log10`` while guarding against zeros."""

    return -np.log10(np.clip(np.asarray(values, dtype=float), 1e-300, 1.0))


def shorten_labels(labels: pd.Series | list[str], max_length: int = 42) -> list[str]:
    """Shorten long labels without losing their left-most identifying text."""

    out: list[str] = []
    for label in pd.Series(labels).astype(str).tolist():
        if len(label) <= max_length:
            out.append(label)
        else:
            out.append(f"{label[: max_length - 1]}…")
    return out


def regional_colors(regional_values: pd.DataFrame) -> tuple[list[str], list[str]]:
    """Return display colors for region-wise plots based on hemisphere/structure."""

    if {"hemisphere", "structure"}.issubset(regional_values.columns):
        colors = []
        for hemisphere, structure in zip(
            regional_values["hemisphere"].astype(str),
            regional_values["structure"].astype(str),
            strict=False,
        ):
            if structure == "subcortex":
                colors.append("#7c3aed" if hemisphere == "L" else "#a855f7")
            else:
                colors.append(POSITIVE if hemisphere == "L" else ACCENT)
        legend = ["L cortex", "R cortex", "L/R subcortex"]
        return colors, legend
    return [ACCENT] * regional_values.shape[0], []


def zscore_for_plot(values: np.ndarray) -> np.ndarray:
    """Return a display-friendly z-scored copy of one value vector."""

    vector = np.asarray(values, dtype=float).reshape(-1)
    centered = vector - np.nanmean(vector)
    scale = np.nanstd(vector, ddof=1)
    if not np.isfinite(scale) or scale == 0:
        return centered
    return centered / scale
