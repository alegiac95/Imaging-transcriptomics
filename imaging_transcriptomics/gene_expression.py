from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from .atlas_registry import get_atlas
from .exceptions import AtlasAssetError, ConfigurationError
from .models import AtlasSelection, AtlasSpec, HemisphereMode, RegionScope


VALID_HEMISPHERES = {"left", "both"}
VALID_REGION_SCOPES = {"default", "all", "cort", "cort+sub"}
DATA_DIR = Path(__file__).resolve().parent / "data"
AHPA_BRAIN_GENES_PATH = DATA_DIR / "filters" / "AHPA_mrna_brain.tsv"


def _packaged_atlas_spec(atlas: str) -> AtlasSpec:
    """Return a packaged atlas specification or raise if assets are missing."""

    spec = get_atlas(atlas)
    if not spec.packaged or spec.expression_path is None or spec.labels_path is None:
        raise AtlasAssetError(
            f"Atlas '{spec.id}' is defined but does not ship packaged expression assets in this branch."
        )
    return spec


def _read_labels(labels_path: Path) -> pd.DataFrame:
    """Load atlas labels into a writable DataFrame copy."""

    labels = _read_labels_cached(str(labels_path)).copy()
    return labels


@lru_cache(maxsize=None)
def _read_labels_cached(labels_path: str) -> pd.DataFrame:
    """Read and normalize a labels table once per path."""

    labels = pd.read_csv(labels_path)
    if "Unnamed: 0" in labels.columns:
        labels = labels.drop(columns=["Unnamed: 0"])
    if labels.columns[0] == "":
        labels = labels.iloc[:, 1:]
    return labels


def _read_expression_values(expression_path: Path) -> np.ndarray:
    """Load the packaged regional-by-gene matrix into a writable array copy."""

    return _read_expression_values_cached(str(expression_path)).copy()


@lru_cache(maxsize=None)
def _read_expression_values_cached(expression_path: str) -> np.ndarray:
    """Read atlas expression values from NPZ or legacy CSV assets."""

    path = Path(expression_path)
    if path.suffix == ".npz":
        with np.load(path, allow_pickle=False) as archive:
            return archive["values"].astype(np.float32, copy=False)

    expression = pd.read_csv(path)
    first_col = str(expression.columns[0]).lstrip("\ufeff")
    if first_col != expression.columns[0]:
        expression = expression.rename(columns={expression.columns[0]: first_col})
    return expression.iloc[:, 2:].to_numpy(dtype=np.float32, copy=False)


def _read_gene_labels(path: Path) -> np.ndarray:
    """Load the shared atlas gene labels into a writable array copy."""

    return _read_gene_labels_cached(str(path)).copy()


@lru_cache(maxsize=None)
def _read_gene_labels_cached(path: str) -> np.ndarray:
    """Read gene labels once per packaged file path."""

    return np.load(path, allow_pickle=False).astype(str, copy=False)


def load_brain_gene_symbols() -> tuple[str, ...]:
    """Load the packaged AHPA brain gene filter as uppercase gene symbols."""

    return _load_brain_gene_symbols_cached(str(AHPA_BRAIN_GENES_PATH))


@lru_cache(maxsize=1)
def _load_brain_gene_symbols_cached(path: str) -> tuple[str, ...]:
    """Read the packaged AHPA brain gene list once per process."""

    table = pd.read_csv(path, sep="\t")
    if "Gene" not in table.columns:
        raise AtlasAssetError("The packaged AHPA brain gene file is missing the required 'Gene' column.")
    genes = (
        table["Gene"]
        .astype("string")
        .dropna()
        .str.strip()
        .str.upper()
    )
    genes = genes.loc[genes != ""].drop_duplicates()
    return tuple(genes.tolist())


def _filter_labels(
    labels: pd.DataFrame,
    hemisphere: HemisphereMode,
    regions: RegionScope,
) -> pd.DataFrame:
    """Restrict atlas labels to the requested hemisphere and region scope."""

    filtered = labels.copy()
    if hemisphere == "left":
        filtered = filtered.loc[filtered["hemisphere"] == "L"]
    elif hemisphere == "both":
        filtered = filtered.loc[filtered["hemisphere"].isin(["L", "R", "B"])]
    else:
        raise ConfigurationError("hemisphere must be either 'left' or 'both'.")

    if regions in {"default", "cort"}:
        filtered = filtered.loc[filtered["structure"] == "cortex"]
    elif regions not in VALID_REGION_SCOPES:
        raise ConfigurationError("regions must be 'default', 'all', 'cort', or 'cort+sub'.")
    return filtered


def _region_lookup_names(labels: pd.DataFrame) -> np.ndarray:
    """Build display names that include hemisphere prefixes when needed."""

    names = labels["label"].astype(str).to_numpy(dtype=str, copy=True)
    hemispheres = labels["hemisphere"].astype(str).to_numpy(dtype=str, copy=False)
    mask = np.isin(hemispheres, ["L", "R"])
    names[mask] = np.char.add(np.char.add(hemispheres[mask], "_"), names[mask])
    return names


def _zscore_expression_columns(values: np.ndarray) -> np.ndarray:
    """Z-score atlas expression per gene while tolerating missing regions.

    Missing values are ignored when estimating the mean and standard deviation.
    Any entries that remain non-finite after normalization are filled with 0,
    which corresponds to the mean after z-scoring.
    """

    numeric = np.asarray(values, dtype=np.float32)
    finite_mask = np.isfinite(numeric)
    safe = np.where(finite_mask, numeric, 0.0)
    counts = finite_mask.sum(axis=0, keepdims=True)
    means = np.divide(
        safe.sum(axis=0, keepdims=True),
        counts,
        out=np.zeros((1, numeric.shape[1]), dtype=np.float32),
        where=counts > 0,
    )
    centered = np.where(finite_mask, numeric - means, 0.0)
    centered_ss = np.sum(centered * centered, axis=0, keepdims=True)
    scale = np.sqrt(
        np.divide(
            centered_ss,
            np.maximum(counts - 1, 1),
            out=np.zeros((1, numeric.shape[1]), dtype=np.float32),
            where=counts > 1,
        )
    )
    scale[scale == 0] = 1.0
    standardized = centered / scale
    standardized[~np.isfinite(standardized)] = 0.0
    return standardized.astype(np.float32, copy=False)


def _select_expression_rows(
    spec: AtlasSpec,
    hemisphere: HemisphereMode,
    regions: RegionScope,
    zscore_expression: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return filtered atlas labels and the matching expression frame."""

    labels = _read_labels(spec.labels_path)
    filtered_labels = _filter_labels(labels, hemisphere=hemisphere, regions=regions)
    row_positions = filtered_labels.index.to_numpy(dtype=np.int32, copy=False)
    values = _read_expression_values(spec.expression_path)
    gene_labels = _read_gene_labels(spec.gene_labels_path)
    selected_values = np.asarray(values[row_positions, :], dtype=np.float32)
    if zscore_expression:
        selected_values = _zscore_expression_columns(selected_values)
    filtered_labels = filtered_labels.reset_index(drop=True)
    selected = pd.DataFrame(selected_values, columns=gene_labels)
    selected.insert(0, "Region", _region_lookup_names(filtered_labels))
    selected.insert(0, "id", filtered_labels["id"].to_numpy(dtype=np.int32, copy=True))
    return filtered_labels, selected


def load_expression_frame(
    atlas: str = "dk",
    hemisphere: HemisphereMode = "left",
    regions: RegionScope = "default",
    zscore_expression: bool = True,
) -> pd.DataFrame:
    """Load a packaged atlas expression table as a DataFrame.

    The returned frame contains the atlas region ``id`` and display ``Region``
    columns followed by one column per gene. By default, each gene is z-scored
    across the selected regions before returning the table.
    """

    spec = _packaged_atlas_spec(atlas)
    _, expression = _select_expression_rows(
        spec,
        hemisphere=hemisphere,
        regions=regions,
        zscore_expression=zscore_expression,
    )
    return expression


def load_gene_labels(atlas: str = "dk") -> np.ndarray:
    """Load the packaged gene labels for an atlas as a column vector."""

    spec = _packaged_atlas_spec(atlas)
    if spec.gene_labels_path is None:
        raise AtlasAssetError(f"Atlas '{spec.id}' does not define a packaged shared gene labels file.")
    return _read_gene_labels(spec.gene_labels_path).astype(object).reshape(-1, 1)


def select_atlas_data(
    atlas: str = "dk",
    hemisphere: HemisphereMode = "left",
    regions: RegionScope = "default",
    zscore_expression: bool = True,
) -> AtlasSelection:
    """Load the packaged labels and expression data for one atlas subset.

    This is the main atlas-loading helper used by the analysis workflows. It
    returns both the filtered region metadata and the selected expression matrix
    wrapped in an :class:`~imaging_transcriptomics.models.AtlasSelection`.
    """

    if hemisphere not in VALID_HEMISPHERES:
        raise ConfigurationError("hemisphere must be either 'left' or 'both'.")
    if regions not in VALID_REGION_SCOPES:
        raise ConfigurationError("regions must be one of 'default', 'all', 'cort', or 'cort+sub'.")

    spec = _packaged_atlas_spec(atlas)
    labels, expression = _select_expression_rows(
        spec,
        hemisphere=hemisphere,
        regions=regions,
        zscore_expression=zscore_expression,
    )
    return AtlasSelection(
        atlas=spec,
        hemisphere=hemisphere,
        regions=regions,
        labels=labels,
        expression=expression,
        gene_labels=expression.columns[2:].to_numpy(dtype=object).reshape(-1, 1),
    )


def expression_matrix(selection: AtlasSelection) -> np.ndarray:
    """Extract only the numeric expression matrix from an atlas selection."""

    return selection.expression.iloc[:, 2:].to_numpy(dtype=float).copy()
