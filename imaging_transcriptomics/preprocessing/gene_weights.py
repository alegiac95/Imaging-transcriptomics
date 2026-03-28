"""Preprocessing helpers for weighted gene tables used by GEDAR-like workflows."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ..exceptions import ConfigurationError
from ..models import RankMode


@dataclass(frozen=True)
class CleanedWeightTable:
    """One cleaned weight table plus its excluded audit rows."""

    table: pd.DataFrame
    excluded_table: pd.DataFrame
    requested_genes: tuple[str, ...]


@dataclass(frozen=True)
class MatchedWeightTable:
    """One cleaned weight table after matching to atlas gene labels."""

    table: pd.DataFrame
    requested_genes: tuple[str, ...]
    matched_genes: tuple[str, ...]
    missing_genes: tuple[str, ...]


def load_weights_table(weights) -> tuple[pd.DataFrame, str]:
    """Load a gene-weight table from a DataFrame or a delimited text file."""

    if isinstance(weights, pd.DataFrame):
        return weights.copy(), "<dataframe>"

    path = Path(str(weights)).expanduser()
    if not path.exists() or not path.is_file():
        raise ConfigurationError(
            "weights must be a pandas DataFrame or a path to a CSV/TSV/text table."
        )
    return pd.read_csv(path, sep=None, engine="python"), str(path)


def validate_rank_selection_args(
    *,
    rank_column: str | None,
    top_percent: float | None,
    top_n: int | None,
    p_threshold: float | None,
) -> None:
    """Ensure optional rank-based selection settings are internally consistent."""

    specified = [value is not None for value in (top_percent, top_n, p_threshold)]
    if sum(specified) > 1:
        raise ConfigurationError("Use at most one of top_percent, top_n, or p_threshold.")
    if any(specified) and rank_column is None:
        raise ConfigurationError("rank_column is required when applying top_percent, top_n, or p_threshold.")
    if top_percent is not None and not 0 < float(top_percent) <= 100:
        raise ConfigurationError("top_percent must be in the interval (0, 100].")
    if top_n is not None and int(top_n) < 1:
        raise ConfigurationError("top_n must be at least 1.")
    if p_threshold is not None and not np.isfinite(float(p_threshold)):
        raise ConfigurationError("p_threshold must be a finite number.")


def clean_weight_table(
    table: pd.DataFrame,
    *,
    gene_column: str,
    weight_column: str,
    rank_column: str | None,
    rank_mode: RankMode,
) -> CleanedWeightTable:
    """Normalize a weight table and audit excluded rows before atlas matching."""

    cleaned = table.copy().reset_index(drop=True)
    cleaned.insert(0, "input_row", np.arange(1, cleaned.shape[0] + 1, dtype=int))

    normalized_gene = cleaned[gene_column].astype("string").fillna("").str.strip()
    cleaned["input_gene"] = normalized_gene.astype(str)
    cleaned["_normalized_gene"] = normalized_gene.astype(str).str.upper()
    cleaned["_weight_value"] = pd.to_numeric(cleaned[weight_column], errors="coerce")
    if rank_column is not None:
        cleaned["_rank_value"] = pd.to_numeric(cleaned[rank_column], errors="coerce")

    excluded_columns = list(cleaned.columns) + ["exclusion_reason", "kept_input_row"]
    excluded_frames: list[pd.DataFrame] = []

    def _exclude(mask: pd.Series | np.ndarray, reason: str) -> None:
        mask_array = np.asarray(mask, dtype=bool)
        if not np.any(mask_array):
            return
        frame = cleaned.loc[mask_array].copy()
        frame["exclusion_reason"] = reason
        frame["kept_input_row"] = pd.NA
        excluded_frames.append(frame)

    invalid_gene_mask = cleaned["input_gene"].eq("") | cleaned["_normalized_gene"].isin({"", "NAN", "<NA>"})
    _exclude(invalid_gene_mask, "missing_gene")

    invalid_weight_mask = ~np.isfinite(cleaned["_weight_value"].to_numpy(dtype=float))
    invalid_weight_mask &= ~np.asarray(invalid_gene_mask, dtype=bool)
    _exclude(invalid_weight_mask, "invalid_weight")

    invalid_rank_mask = np.zeros(cleaned.shape[0], dtype=bool)
    if rank_column is not None:
        invalid_rank_mask = ~np.isfinite(cleaned["_rank_value"].to_numpy(dtype=float))
        invalid_rank_mask &= ~np.asarray(invalid_gene_mask, dtype=bool)
        invalid_rank_mask &= ~np.asarray(invalid_weight_mask, dtype=bool)
        _exclude(invalid_rank_mask, "invalid_rank")

    valid_mask = ~(
        np.asarray(invalid_gene_mask, dtype=bool)
        | np.asarray(invalid_weight_mask, dtype=bool)
        | np.asarray(invalid_rank_mask, dtype=bool)
    )
    valid = cleaned.loc[valid_mask].copy()
    if valid.empty:
        excluded = pd.concat(excluded_frames, ignore_index=True) if excluded_frames else pd.DataFrame(columns=excluded_columns)
        if not excluded.empty:
            excluded = excluded.loc[:, excluded_columns]
        return CleanedWeightTable(table=valid, excluded_table=excluded, requested_genes=())

    if rank_column is not None:
        valid = valid.sort_values(
            by=["_rank_value", "input_row"],
            ascending=[rank_mode == "ascending", True],
            kind="mergesort",
        )
    else:
        valid = valid.sort_values(by="input_row", kind="mergesort")

    duplicate_mask = valid["_normalized_gene"].duplicated(keep="first")
    if duplicate_mask.any():
        duplicate_rows = valid.loc[duplicate_mask].copy()
        kept_rows = valid.loc[~duplicate_mask, ["_normalized_gene", "input_row"]].rename(
            columns={"input_row": "kept_input_row"}
        )
        duplicate_rows = duplicate_rows.merge(kept_rows, on="_normalized_gene", how="left")
        duplicate_rows["exclusion_reason"] = "duplicate_symbol"
        excluded_frames.append(duplicate_rows)
        valid = valid.loc[~duplicate_mask].copy()

    valid = valid.sort_values(by="input_row", kind="mergesort").reset_index(drop=True)
    valid[gene_column] = valid["input_gene"]
    valid[weight_column] = valid["_weight_value"].to_numpy(dtype=float)
    if rank_column is not None:
        valid[rank_column] = valid["_rank_value"].to_numpy(dtype=float)

    excluded = pd.concat(excluded_frames, ignore_index=True) if excluded_frames else pd.DataFrame(columns=excluded_columns)
    if not excluded.empty:
        excluded = excluded.sort_values(by="input_row", kind="mergesort").loc[:, excluded_columns].reset_index(drop=True)
    else:
        excluded = pd.DataFrame(columns=excluded_columns)

    requested_genes = tuple(valid[gene_column].astype(str).tolist())
    valid = valid.drop(
        columns=["_normalized_gene", "_weight_value"] + (["_rank_value"] if rank_column is not None else [])
    )
    return CleanedWeightTable(table=valid, excluded_table=excluded, requested_genes=requested_genes)


def apply_brain_gene_filter_to_weights(
    cleaned: CleanedWeightTable,
    *,
    gene_column: str,
    brain_gene_symbols: set[str],
) -> CleanedWeightTable:
    """Restrict one cleaned weight table to genes present in the brain filter."""

    mask = cleaned.table[gene_column].astype(str).str.upper().isin(brain_gene_symbols)
    if bool(np.all(mask.to_numpy(dtype=bool))):
        return cleaned

    kept = cleaned.table.loc[mask].reset_index(drop=True)
    excluded = cleaned.table.loc[~mask].copy()
    excluded["exclusion_reason"] = "outside_brain_filter"
    excluded["kept_input_row"] = pd.NA

    if cleaned.excluded_table.empty:
        combined = excluded
    else:
        combined = pd.concat([cleaned.excluded_table, excluded], ignore_index=True, sort=False)

    combined = combined.sort_values(by="input_row", kind="mergesort").reset_index(drop=True)
    return CleanedWeightTable(
        table=kept,
        excluded_table=combined,
        requested_genes=tuple(kept[gene_column].astype(str).tolist()),
    )


def match_weight_table_to_genes(
    table: pd.DataFrame,
    *,
    gene_column: str,
    weight_column: str,
    rank_column: str | None,
    available_genes: np.ndarray,
) -> MatchedWeightTable:
    """Match one cleaned weight table to the canonical atlas gene symbols."""

    lookup = {gene.upper(): gene for gene in np.asarray(available_genes, dtype=str).tolist()}
    requested = tuple(table[gene_column].astype(str).tolist())
    matched_rows: list[dict[str, object]] = []
    matched_genes: list[str] = []
    missing_genes: list[str] = []

    for _, row in table.iterrows():
        gene = str(row[gene_column])
        canonical = lookup.get(gene.upper())
        if canonical is None:
            missing_genes.append(gene)
            continue
        matched_genes.append(canonical)
        payload: dict[str, object] = {
            "gene": canonical,
            "input_gene": gene,
            "input_weight": float(row[weight_column]),
        }
        if rank_column is not None:
            payload["rank_value"] = float(row[rank_column])
        matched_rows.append(payload)

    return MatchedWeightTable(
        table=pd.DataFrame(matched_rows),
        requested_genes=requested,
        matched_genes=tuple(matched_genes),
        missing_genes=tuple(missing_genes),
    )


def apply_rank_selection_mask(
    frame: pd.DataFrame,
    *,
    rank_mode: RankMode,
    top_percent: float | None,
    top_n: int | None,
    p_threshold: float | None,
) -> np.ndarray:
    """Select rows from one matched table using rank-based filters."""

    if rank_mode not in {"ascending", "descending"}:
        raise ConfigurationError("rank_mode must be either 'ascending' or 'descending'.")

    selected = np.ones(frame.shape[0], dtype=bool)
    if "rank_value" not in frame.columns:
        return selected

    values = frame["rank_value"].to_numpy(dtype=float)
    order = np.argsort(values, kind="mergesort")
    if rank_mode == "descending":
        order = order[::-1]

    selected[:] = False
    if top_percent is not None:
        keep = max(1, int(np.ceil(frame.shape[0] * float(top_percent) / 100.0)))
        selected[order[:keep]] = True
        return selected
    if top_n is not None:
        selected[order[: int(top_n)]] = True
        return selected
    if p_threshold is not None:
        if rank_mode == "ascending":
            return values <= float(p_threshold)
        return values >= float(p_threshold)

    selected[:] = True
    return selected
