from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import zscore

from .atlas_registry import get_atlas
from .models import AtlasSelection, AtlasSpec, HemisphereMode, RegionScope


VALID_HEMISPHERES = {"left", "both"}
VALID_REGION_SCOPES = {"all", "cort", "cort+sub"}


def _packaged_atlas_spec(atlas: str) -> AtlasSpec:
    spec = get_atlas(atlas)
    if not spec.packaged or spec.expression_path is None or spec.labels_path is None:
        raise FileNotFoundError(
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


def _read_expression(expression_path: Path) -> pd.DataFrame:
    return _read_expression_cached(str(expression_path)).copy()


@lru_cache(maxsize=None)
def _read_expression_cached(expression_path: str) -> pd.DataFrame:
    path = Path(expression_path)
    if path.suffix == ".npz":
        with np.load(path, allow_pickle=False) as archive:
            ids = archive["ids"].astype(np.int32, copy=False)
            regions = archive["regions"].astype(str, copy=False)
            genes = archive["genes"].astype(str, copy=False)
            values = archive["values"].astype(np.float32, copy=False)
        frame = pd.DataFrame(values, columns=genes)
        frame.insert(0, "Region", regions)
        frame.insert(0, "id", ids)
        return frame

    expression = pd.read_csv(path)
    first_col = str(expression.columns[0]).lstrip("\ufeff")
    if first_col != expression.columns[0]:
        expression = expression.rename(columns={expression.columns[0]: first_col})
    return expression


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
        raise ValueError("hemisphere must be either 'left' or 'both'.")

    if regions == "cort":
        filtered = filtered.loc[filtered["structure"] == "cortex"]
    elif regions not in VALID_REGION_SCOPES:
        raise ValueError("regions must be 'all', 'cort', or 'cort+sub'.")
    return filtered.reset_index(drop=True)


def _resolve_expression_lookup(
    labels: pd.DataFrame,
    expression: pd.DataFrame,
) -> tuple[list[str], str]:
    if "Region" not in expression.columns:
        return labels["label"].astype(str).tolist(), "label"

    region_names = expression["Region"].astype(str)
    region_set = set(region_names)
    raw_names = labels["label"].astype(str).tolist()
    if all(name in region_set for name in raw_names):
        return raw_names, "Region"

    lookup: list[str] = []
    for row in labels.itertuples(index=False):
        label = str(row.label)
        hemisphere = str(row.hemisphere)
        if label in region_set:
            lookup.append(label)
        elif hemisphere in {"L", "R"} and f"{hemisphere}_{label}" in region_set:
            lookup.append(f"{hemisphere}_{label}")
        elif hemisphere == "B" and label in region_set:
            lookup.append(label)
        else:
            lookup.append(label)
    return lookup, "Region"


def _select_expression_rows(
    spec: AtlasSpec,
    hemisphere: HemisphereMode,
    regions: RegionScope,
    zscore_expression: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    labels = _read_labels(spec.labels_path)
    expression = _read_expression(spec.expression_path)
    filtered_labels = _filter_labels(labels, hemisphere=hemisphere, regions=regions)

    lookup_names, index_col = _resolve_expression_lookup(filtered_labels, expression)
    expression_index = expression.set_index(index_col)
    missing = sorted(set(lookup_names) - set(expression_index.index.astype(str)))
    if missing:
        raise ValueError(
            f"Atlas '{spec.id}' is missing expression rows for labels: {', '.join(missing[:5])}"
        )
    selected = expression_index.loc[lookup_names].copy().reset_index()
    if zscore_expression:
        zscored = pd.DataFrame(
            zscore(selected.iloc[:, 2:].to_numpy(dtype=float), axis=0, ddof=1),
            columns=selected.columns[2:],
            index=selected.index,
        )
        selected = pd.concat([selected.iloc[:, :2].copy(), zscored], axis=1)
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
    expression = _read_expression(spec.expression_path)
    return expression.columns[2:].to_numpy(dtype=object).reshape(-1, 1)


def select_atlas_data(
    atlas: str = "dk",
    hemisphere: HemisphereMode = "left",
    regions: RegionScope = "all",
) -> AtlasSelection:
    if hemisphere not in VALID_HEMISPHERES:
        raise ValueError("hemisphere must be either 'left' or 'both'.")
    if regions not in VALID_REGION_SCOPES:
        raise ValueError("regions must be one of 'all', 'cort', or 'cort+sub'.")

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
