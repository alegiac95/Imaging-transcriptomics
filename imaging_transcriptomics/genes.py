"""Public facade for gene-statistics containers used by analysis workflows.

The heavy lifting for correlation- and PLS-specific gene statistics lives in
the internal :mod:`imaging_transcriptomics.gene_stats` package. This module
keeps the historical import surface stable for the rest of the package and for
downstream users importing :mod:`imaging_transcriptomics.genes`.
"""

from .ora import ora_from_gene_table
from .gene_stats.correlation import CorrGenes
from .gene_stats.factory import GeneResults
from .gene_stats.pls import BootPLS, OrigPLS, PLSGenes

__all__ = [
    "BootPLS",
    "CorrGenes",
    "GeneResults",
    "OrigPLS",
    "PLSGenes",
    "ora_from_gene_table",
]
