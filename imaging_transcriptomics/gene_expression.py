from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import zscore

from .atlas_registry import get_atlas
from .exceptions import AtlasAssetError, ConfigurationError
from .models import AtlasSelection, AtlasSpec, HemisphereMode, RegionScope


VALID_HEMISPHERES = {"left", "both"}
VALID_REGION_SCOPES = {"all", "cort", "cort+sub"}


def _packaged_atlas_spec(atlas: str) -> AtlasSpec:
    spec = get_atlas(atlas)
    if not spec.packaged or spec.expression_path is None or spec.labels_path is None:
        raise AtlasAssetError(
            f"Atlas '{spec.id}' is defined but does not ship packaged expression assets in this branch."
        )
    return spec


def _read_labels(labels_path: Path) -> pd.DataFrame:
    labels = _read_labels_cached(str(labels_path)).copy()
    return labels


@lru_cache(maxsize=None)
def _read_labels_cached(labels_path: str) -> pd.DataFrame:
    labels = pd.read_csv(labels_path)
    if "Unnamed: 0" in labels.columns:
        labels = labels.drop(columns=["Unnamed: 0"])
    if labels.columns[0] == "":
        labels = labels.iloc[:, 1:]
    return labels


def _read_expression_values(expression_path: Path) -> np.ndarray:
    return _read_expression_values_cached(str(expression_path)).copy()


@lru_cache(maxsize=None)
def _read_expression_values_cached(expression_path: str) -> np.ndarray:
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
    return _read_gene_labels_cached(str(path)).copy()


@lru_cache(maxsize=None)
def _read_gene_labels_cached(path: str) -> np.ndarray:
    return np.load(path, allow_pickle=False).astype(str, copy=False)


def _filter_labels(
    labels: pd.DataFrame,
    hemisphere: HemisphereMode,
    regions: RegionScope,
) -> pd.DataFrame:
    filtered = labels.copy()
    if hemisphere == "left":
        filtered = filtered.loc[filtered["hemisphere"] == "L"]
    elif hemisphere == "both":
        filtered = filtered.loc[filtered["hemisphere"].isin(["L", "R", "B"])]
    else:
        raise ConfigurationError("hemisphere must be either 'left' or 'both'.")

    if regions == "cort":
        filtered = filtered.loc[filtered["structure"] == "cortex"]
    elif regions not in VALID_REGION_SCOPES:
        raise ConfigurationError("regions must be 'all', 'cort', or 'cort+sub'.")
    return filtered


def _region_lookup_names(labels: pd.DataFrame) -> np.ndarray:
    names = labels["label"].astype(str).to_numpy(dtype=str, copy=True)
    hemispheres = labels["hemisphere"].astype(str).to_numpy(dtype=str, copy=False)
    mask = np.isin(hemispheres, ["L", "R"])
    names[mask] = np.char.add(np.char.add(hemispheres[mask], "_"), names[mask])
    return names


def _select_expression_rows(
    spec: AtlasSpec,
    hemisphere: HemisphereMode,
    regions: RegionScope,
    zscore_expression: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    labels = _read_labels(spec.labels_path)
    filtered_labels = _filter_labels(labels, hemisphere=hemisphere, regions=regions)
    row_positions = filtered_labels.index.to_numpy(dtype=np.int32, copy=False)
    values = _read_expression_values(spec.expression_path)
    gene_labels = _read_gene_labels(spec.gene_labels_path)
    selected_values = np.asarray(values[row_positions, :], dtype=np.float32)
    if zscore_expression:
        selected_values = zscore(selected_values, axis=0, ddof=1).astype(np.float32, copy=False)
    filtered_labels = filtered_labels.reset_index(drop=True)
    selected = pd.DataFrame(selected_values, columns=gene_labels)
    selected.insert(0, "Region", _region_lookup_names(filtered_labels))
    selected.insert(0, "id", filtered_labels["id"].to_numpy(dtype=np.int32, copy=True))
    return filtered_labels, selected


def load_expression_frame(
    atlas: str = "dk",
    hemisphere: HemisphereMode = "left",
    regions: RegionScope = "all",
    zscore_expression: bool = True,
) -> pd.DataFrame:
    spec = _packaged_atlas_spec(atlas)
    _, expression = _select_expression_rows(
        spec,
        hemisphere=hemisphere,
        regions=regions,
        zscore_expression=zscore_expression,
    )
    return expression


def load_gene_labels(atlas: str = "dk") -> np.ndarray:
    spec = _packaged_atlas_spec(atlas)
    if spec.gene_labels_path is None:
        raise AtlasAssetError(f"Atlas '{spec.id}' does not define a packaged shared gene labels file.")
    return _read_gene_labels(spec.gene_labels_path).astype(object).reshape(-1, 1)


def select_atlas_data(
    atlas: str = "dk",
    hemisphere: HemisphereMode = "left",
    regions: RegionScope = "all",
) -> AtlasSelection:
    if hemisphere not in VALID_HEMISPHERES:
        raise ConfigurationError("hemisphere must be either 'left' or 'both'.")
    if regions not in VALID_REGION_SCOPES:
        raise ConfigurationError("regions must be one of 'all', 'cort', or 'cort+sub'.")

    spec = _packaged_atlas_spec(atlas)
    labels, expression = _select_expression_rows(
        spec,
        hemisphere=hemisphere,
        regions=regions,
        zscore_expression=True,
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
    return selection.expression.iloc[:, 2:].to_numpy(dtype=float).copy()
