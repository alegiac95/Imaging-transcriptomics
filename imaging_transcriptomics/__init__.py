__version__ = "0.1.0"  # Do not change line for the version or setup.py won't work

from . import errors
from . import bootstrap
from . import inputs
from . import permutations
from . import transcriptomics
from . import genes
from . import reporting

from .transcriptomics import ImagingTranscriptomics
from .bootstrap import (bootstrap_pls,
                        bootstrap_genes)
from .inputs import (load_gene_expression,
                     load_gene_labels,
                     get_components,
                     read_scan,
                     extract_average)
from .genes import GeneResults
