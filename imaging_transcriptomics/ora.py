from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom, norm
from statsmodels.stats.multitest import multipletests

from .genesets import get_geneset


ORA_COLUMNS = [
    "Term",
    "overlap_size",
    "set_size",
    "selected_size",
    "universe_size",
    "enrichment_ratio",
    "odds_ratio",
    "odds_ratio_ci_low",
    "odds_ratio_ci_high",
    "p_value",
    "fdr",
    "overlap_genes",
]


def _empty_ora_frame() -> pd.DataFrame:
    """Return an empty ORA result frame with the standard output columns."""

    return pd.DataFrame({column: pd.Series(dtype=object) for column in ORA_COLUMNS})


def _parse_gmt(path: Path) -> dict[str, set[str]]:
    """Parse a GMT file into uppercase gene-set membership tables."""

    genesets: dict[str, set[str]] = {}
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        entries = line.split("\t")
        if len(entries) < 3:
            continue
        term = entries[0]
        genes = {gene.strip().upper() for gene in entries[2:] if gene.strip()}
        if genes:
            genesets[term] = genes
    if not genesets:
        raise ValueError(f"No valid gene sets were found in {path}.")
    return genesets


def load_ora_genesets(gene_set: str) -> dict[str, set[str]]:
    """Resolve and load a GMT geneset resource for ORA."""

    resolved = get_geneset(gene_set)
    path = Path(resolved)
    if not path.exists() or not path.is_file() or path.suffix != ".gmt":
        raise ValueError(
            "ORA currently supports packaged gene sets or local .gmt files. "
            f"Received {gene_set!r}, which resolved to {resolved!r}."
        )
    return _parse_gmt(path)


def _selected_genes(
    gene_table: pd.DataFrame,
    *,
    score_column: str,
    p_value_column: str,
    p_threshold: float,
    direction: str,
) -> tuple[dict[str, str], set[str]]:
    """Select significant genes for one ORA direction and preserve display names."""

    selected = gene_table.loc[gene_table[p_value_column] <= p_threshold].copy()
    if direction == "up":
        selected = selected.loc[selected[score_column] > 0]
    elif direction == "down":
        selected = selected.loc[selected[score_column] < 0]
    else:
        raise ValueError(f"Unsupported ORA direction: {direction}")
    selected["gene_upper"] = selected["gene"].astype(str).str.upper()
    mapping = dict(zip(selected["gene_upper"], selected["gene"].astype(str), strict=False))
    return mapping, set(mapping)


def _odds_ratio(
    overlap_size: int,
    set_size: int,
    selected_size: int,
    universe_size: int,
) -> float:
    """Compute the contingency-table odds ratio for one ORA term."""

    a = float(overlap_size)
    b = float(selected_size - overlap_size)
    c = float(set_size - overlap_size)
    d = float(universe_size - selected_size - set_size + overlap_size)
    denominator = b * c
    numerator = a * d
    if denominator == 0:
        if numerator > 0:
            return float("inf")
        return float("nan")
    return numerator / denominator


def _odds_ratio_confidence_interval(
    overlap_size: int,
    set_size: int,
    selected_size: int,
    universe_size: int,
    *,
    confidence: float = 0.95,
) -> tuple[float, float]:
    """Compute a log-odds confidence interval with sparse-table correction."""

    a = float(overlap_size)
    b = float(selected_size - overlap_size)
    c = float(set_size - overlap_size)
    d = float(universe_size - selected_size - set_size + overlap_size)
    if min(a, b, c, d) < 0:
        return float("nan"), float("nan")

    # Haldane-Anscombe continuity correction for sparse tables.
    if min(a, b, c, d) == 0:
        a += 0.5
        b += 0.5
        c += 0.5
        d += 0.5

    z_value = float(norm.ppf(0.5 + confidence / 2.0))
    log_or = np.log((a * d) / (b * c))
    standard_error = np.sqrt((1.0 / a) + (1.0 / b) + (1.0 / c) + (1.0 / d))
    interval = z_value * standard_error
    return float(np.exp(log_or - interval)), float(np.exp(log_or + interval))


def ora_from_gene_table(
    gene_table: pd.DataFrame,
    *,
    gene_set: str,
    score_column: str,
    p_value_column: str = "p_value",
    p_threshold: float = 0.05,
) -> dict[str, pd.DataFrame]:
    """Run ORA on significant positive and negative gene sets from one table."""

    if not 0 < float(p_threshold) <= 1:
        raise ValueError("ORA p-threshold must be in the interval (0, 1].")
    required = {"gene", score_column, p_value_column}
    missing = sorted(required - set(gene_table.columns))
    if missing:
        raise ValueError(f"Gene table is missing required columns: {', '.join(missing)}")

    genesets = load_ora_genesets(gene_set)
    working = gene_table.copy()
    working["gene"] = working["gene"].astype(str)
    working["gene_upper"] = working["gene"].str.upper()
    universe_map = dict(zip(working["gene_upper"], working["gene"], strict=False))
    universe = set(universe_map)
    universe_size = len(universe)

    results: dict[str, pd.DataFrame] = {}
    for direction in ("up", "down"):
        selected_map, selected = _selected_genes(
            working,
            score_column=score_column,
            p_value_column=p_value_column,
            p_threshold=float(p_threshold),
            direction=direction,
        )
        selected_size = len(selected)
        if selected_size == 0:
            results[direction] = _empty_ora_frame()
            continue

        rows: list[dict[str, object]] = []
        for term, genes in genesets.items():
            genes_in_universe = genes & universe
            set_size = len(genes_in_universe)
            if set_size == 0:
                continue
            overlap = selected & genes_in_universe
            overlap_size = len(overlap)
            if overlap_size == 0:
                continue
            expected = (set_size / universe_size) * selected_size
            enrichment_ratio = overlap_size / expected if expected else np.nan
            odds_ratio = _odds_ratio(
                overlap_size,
                set_size,
                selected_size,
                universe_size,
            )
            odds_ratio_ci_low, odds_ratio_ci_high = _odds_ratio_confidence_interval(
                overlap_size,
                set_size,
                selected_size,
                universe_size,
            )
            p_value = float(
                hypergeom.sf(overlap_size - 1, universe_size, set_size, selected_size)
            )
            overlap_genes = ";".join(sorted(selected_map[gene] for gene in overlap))
            rows.append(
                {
                    "Term": term,
                    "overlap_size": overlap_size,
                    "set_size": set_size,
                    "selected_size": selected_size,
                    "universe_size": universe_size,
                    "enrichment_ratio": enrichment_ratio,
                    "odds_ratio": odds_ratio,
                    "odds_ratio_ci_low": odds_ratio_ci_low,
                    "odds_ratio_ci_high": odds_ratio_ci_high,
                    "p_value": p_value,
                    "overlap_genes": overlap_genes,
                }
            )

        if not rows:
            results[direction] = _empty_ora_frame()
            continue

        frame = pd.DataFrame(rows)
        _, fdr, _, _ = multipletests(frame["p_value"].to_numpy(dtype=float), method="fdr_bh")
        frame["fdr"] = fdr
        frame = frame[ORA_COLUMNS].sort_values(
            ["fdr", "p_value", "enrichment_ratio"],
            ascending=[True, True, False],
            kind="mergesort",
        ).reset_index(drop=True)
        results[direction] = frame
    return results
