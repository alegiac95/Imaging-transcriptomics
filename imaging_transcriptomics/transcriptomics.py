import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import zscore
from pyls import pls_regression

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from netneurotools import freesurfer, stats

from .inputs import (load_gene_expression,
                     load_gene_labels,
                     get_components)
from .bootstrap import bootstrap_pls, bootstrap_genes


class ImagingTranscriptomics:
    def __init__(self, scan_data, **kwargs):
        """Initialise the imaging transcriptomics class with the input scan's data and number of components or variance
        explained.

        :param array-like scan_data: average values in the ROI defined by the Desikan-Killiany atlas.
        :param int n_components: number of components to use for the PLS regression.
        :param int variance: total explained variance by the PLS components.
        """
        self.scan_data = scan_data
        self.zscore_data = zscore(scan_data, ddof=1, axis=0)
        self.n_components = kwargs.get("n_components")
        self.var = kwargs.get("variance")
        self.__cortical = self.zscore_data[0:34].reshape(34, 1)
        self.__subcortical = self.zscore_data[34:].reshape(7, 1)
        self.__gene_expression = load_gene_expression()
        self.__gene_labels = load_gene_labels()
        # Initialise with defaults for later
        self.__permuted = None
        self.r_boot = None
        self.p_boot = None
        self.gene_results = None
        self.var_components = None

    def __permute_data(self, iterations=1_000):
        """Permute the scan data for the analysis.

        The permutations are computed into cortical and subcortical regions separately and then merged. This is done
        to maintain the spatial autocorrelation in the cortical regions for more accuracy.
        To compute the cortical permutations the library python package ``netneurotools`` developed by R. Markello is
        used. For more information about the methods used you can refer to the official `documentation of the
        package. <https://netneurotools.readthedocs.io/en/latest/>_`

        :param int iterations: number of iterations to perform in the permutations.
        """
        self.__permuted = np.zeros((self.zscore_data.shape[0], iterations))
        # subcortical
        sub_permuted = np.array(
            [np.random.permutation(self.__subcortical) for _ in range(iterations)]
        ).reshape(7, iterations)
        self.__permuted[34:, :] = sub_permuted
        # Cortical
        # Annotation file for the Desikan-Killiany atlas in fs5
        annot_lh = Path(__file__).resolve().parent.parent / "data/fsa5_lh_aparc.annot"
        annot_rh = Path(__file__).resolve().parent.parent / "data/fsa5_rh_aparc.annot"
        # Get the parcel centroids of the Desikan-Killiany atlas
        parcel_centroids, parcel_hemi = freesurfer.find_parcel_centroids(
            lhannot=annot_lh,
            rhannot=annot_rh,
            version='fsaverage5',
            surf='sphere',
            method="surface")
        # Mask the results to have only the left hemisphere
        left_hemi_mask = parcel_hemi == 0
        parcel_centroids, parcel_hemi = parcel_centroids[left_hemi_mask], parcel_hemi[left_hemi_mask]
        # Get the spin samples
        spins = stats.gen_spinsamples(parcel_centroids, parcel_hemi, n_rotate=iterations, method='vasa', seed=1234)
        cort_permuted = np.array(self.__cortical[spins]).reshape(34, iterations)
        self.__permuted[0:34, :] = cort_permuted

    def save_permutations(self, path):
        """Save the permutations to a csv file at a specified path.

        :param path: Path used to save the permutations, this *should* also include the name of the file, e.g.,
        "~/Documents/my_permuted.csv"
        """
        if self.__permuted is None:
            raise AttributeError("There are no permutations of the scan available to save. Before saving the "
                                 "permutations you need to compute them.")
        pd.DataFrame(self.__permuted).to_csv(Path(path), header=None, index=False)

    def pls_all_components(self):
        """Compute a PLS regression with all components.

        After the regression is estimated, either the number of components or the estimated percentage of variance
        given by the components is estimated, depending on what is set by the user in the __init__() method.
        """
        results = pls_regression(self.__gene_expression, self.zscore_data.reshape(41, 1),
                                 n_components=15, n_perm=0, n_boot=0)
        var_exp = results.get("varexp")
        if self.n_components is None and self.var != 0.0:
            self.n_components = get_components((self.var / 100), var_exp)
        elif self.var is None and self.n_components != 0:
            self.var = np.cumsum(var_exp)[self.n_components-1]
        self.var_components = var_exp

    def run(self, n_iter=1_000):
        """Run the analysis of the imaging scan.

        :param int n_iter: number of permutations to make.
        """
        self.pls_all_components()
        self.__permute_data(iterations=n_iter)
        self.r_boot, self.p_boot = bootstrap_pls(self.__gene_expression,
                                                 self.zscore_data.reshape(41, 1),
                                                 self.__permuted,
                                                 self.n_components,
                                                 iterations=n_iter)
        self.gene_results = bootstrap_genes(self.__gene_expression,
                                            self.zscore_data.reshape(41, 1),
                                            self.n_components,
                                            self.scan_data,
                                            self.__gene_labels,
                                            n_iter)
        self.gene_results.boot_results.compute_values(self.n_components,
                                                      self.gene_results.original_results.pls_weights,
                                                      self.gene_results.original_results.pls_gene)
