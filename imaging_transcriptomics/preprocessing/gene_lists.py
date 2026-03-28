"""Reusable helpers for parsing and resolving user-provided gene lists."""

from __future__ import annotations

from pathlib import Path
from typing import Collection, Iterable

import numpy as np
import pandas as pd


def dedupe_preserve_order(values: Iterable[str]) -> list[str]:
    """Drop empty and duplicate gene symbols while preserving their order."""

    seen: set[str] = set()
    deduped: list[str] = []
    for value in values:
        gene = str(value).strip()
        if not gene:
            continue
        key = gene.upper()
        if key in seen:
            continue
        seen.add(key)
        deduped.append(gene)
    return deduped


def parse_gene_tokens(text: str) -> list[str]:
    """Parse gene symbols from free text, one-per-line files, or CSV-like text."""

    tokens: list[str] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        normalized = line.replace(",", " ").replace("\t", " ")
        tokens.extend(part for part in normalized.split() if part)
    return dedupe_preserve_order(tokens)


def load_gene_list(genes: str | Path | Iterable[str]) -> list[str]:
    """Load a gene list from a path, iterable, or inline comma-separated string."""

    if isinstance(genes, (list, tuple, set, np.ndarray, pd.Series)):
        return dedupe_preserve_order(str(item) for item in genes)
    path = Path(str(genes)).expanduser()
    if path.exists() and path.is_file():
        return parse_gene_tokens(path.read_text())
    return parse_gene_tokens(str(genes))


def partition_brain_filtered_genes(
    requested_genes: Iterable[str],
    brain_gene_symbols: Collection[str],
) -> tuple[list[str], list[str]]:
    """Split genes into retained brain genes and filtered-out non-brain genes."""

    brain_genes = {str(gene).upper() for gene in brain_gene_symbols}
    retained: list[str] = []
    filtered: list[str] = []
    for gene in requested_genes:
        if str(gene).upper() in brain_genes:
            retained.append(str(gene))
        else:
            filtered.append(str(gene))
    return retained, filtered


def resolve_selected_genes(
    requested_genes: Iterable[str],
    available_genes: np.ndarray,
) -> tuple[list[str], list[str]]:
    """Split requested genes into canonical atlas matches and missing symbols."""

    lookup: dict[str, str] = {}
    for gene in np.asarray(available_genes, dtype=str).tolist():
        lookup.setdefault(gene.upper(), gene)

    matched: list[str] = []
    missing: list[str] = []
    for gene in requested_genes:
        canonical = lookup.get(str(gene).upper())
        if canonical is None:
            missing.append(str(gene))
        else:
            matched.append(canonical)
    return matched, missing
