from pathlib import Path

import pandas as pd
import numpy as np
from scipy.stats import zscore

from .inputs import (load_gene_expression,
                     load_gene_labels)


class GeneResults(dict):
    """Class to save the gene results."""
    def __init__(self, n_comp):
        super().__init__()
        self.pls_weights = [None] * n_comp
        self.pls_gene = [None] * n_comp
        self.gene_id = [None] * n_comp


class ImagingTranscriptomics:
    def __init__(self, scan_data, **kwargs):
        """Initialise the imaging transcriptomics class with the input scan's data and number of components or variance
        explained.

        :param array-like scan_data: average values in the ROI defined by the Desikan-Killiany atlas.
        :param int n_components: number of components to use for the PLS regression.
        :param int variance: total explained variance by the PLS components.
        """
        self.scan_data = zscore(scan_data, ddof=1, axis=0)
        self.n_components = kwargs.get("n_components")
        self.var = kwargs.get("variance")
        self.__cortical = zscore(scan_data[0:34], ddof=1, axis=0)
        self.__subcortical = zscore(scan_data[34:], ddof=1, axis=0)
        self.__permuted = None
        self.__gene_expression = load_gene_expression()
        self.__gene_labels = load_gene_labels()

    def __permute_data(self, iterations=1_000):
        """Permute the scan data for the analysis.

        :param int iterations: number of iterations to perform in the permutations.
        """
        self.__permuted = np.zeros((self.scan_data.shape[0], iterations))
        # subcortical
        sub_permuted = np.array(
            [np.random.permutation(self.__subcortical) for _ in range(iterations)]
        ).reshape(7, iterations)
        self.__permuted[34:, :] = sub_permuted
        # cortical

    def save_permutations(self, path):
        """Save the permutations to a csv file at a specified path.

        :param path: Path used to save the permutations, this *should* also include the name of the file, e.g.,
        "~/Documents/my_permuted.csv"
        """
        if self.__permuted is not None:
            pd.DataFrame(self.__permuted).to_csv(Path(path), header=None, index=False)
        else:
            raise AttributeError("There are no permutations of the scan available to save. Before saving the "
                                 "permutations you need to compute them.")
        return

    def run(self, n_iter=1_000):
        """Run the analysis of the imaging scan.

        :param int n_iter: number of permutations to make.
        """
        self.__permute_data(iterations=n_iter)
        pass
