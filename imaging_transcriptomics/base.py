import numpy as np
from pathlib import Path
from scipy.stats import zscore, spearmanr
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from netneurotools import freesurfer, stats
from pyls import pls_regression
import pandas as pd
from .inputs import load_gene_expression, extract_average, read_scan


def load_gene_labels():
    """Return an array with the gene labels.
    The gene labels are available in the ``data`` sub-folder.

    :return: numpy array with the labels of the genes.
    """
    genes_labels_path = (
        Path(__file__).resolve().parent / "data" / "gene_expression_labels.txt"
    )
    return pd.read_fwf(genes_labels_path, header=None)[0].tolist()


class Imt:
    # --------- INITIALIZATION --------- #
    def __init__(self,
                 scan_data,
                 regions="cort+sub",
                 method="pls",
                 **kwargs):
        """ImagingTranscriptomics class for imaging transcriptomics analysis.

        :param np.array scan_data: imaging scan data.
        :param str regions: regions to be used for analysis.
        :param str method: method to run the analysis, can be either "pls"
        for pls regression or "corr" cor simple correlation analysis.
        :param kwargs: additional arguments for the method. This include:
            * "n_components": number of components for pls regression.
            * "var": variance explained for pls regression.
            * "n_permutations": number of permutations for permutation test.

        """
        if regions == "cort+sub" or regions == "all":
            assert scan_data.shape == (41,)
            self.zscore_data = zscore(scan_data, axis=0, ddof=1)
            self._cortical = self.zscore_data[:34]
            self._subcortical = self.zscore_data[34:]
        elif regions == "cort":
            assert scan_data.shape == (34,) or scan_data.shape == (41,)
            self.zscore_data = zscore(scan_data, axis=0, ddof=1)
            self._cortical = scan_data
            self._subcortical = None
        self._regions = regions
        self.scan_data = scan_data
        if method not in ["pls", "corr"]:
            raise ValueError(
                "The method must be either pls or corr."
                "Please choose either pls or corr to run the analysis."
            )
        else:
            self._method = method
        if self._method == "pls":
            if "n_components" not in kwargs and "var" not in kwargs:
                raise ValueError("You must specify either the variance or "
                                 "the number of components for pls regression")
            else:
                self._analysis = PLSAnalysis(kwargs.get("n_components"),
                                             kwargs.get("var"))
        elif self._method == "corr":
            self._analysis = CorrAnalysis()
        self.gene_expression = load_gene_expression(self._regions)
        self.gene_labels = load_gene_labels()
        self._permutations = None

    @classmethod
    def from_scan(cls, scan_path, method="pls", regions="cort+sub", **kwargs):
        """Initialise an ImagingTranscriptomics object from a NIfTI scan.
        The extracted data corresponds to the average of the ROIs in the DK
        atlas.
        """
        scan_data = extract_average(read_scan(scan_path))
        return cls(scan_data, method=method, regions=regions, **kwargs)

    @classmethod
    def from_file(cls, file_path, method="pls", regions="cort+sub", **kwargs):
        """Initialise an ImagingTranscriptomics object from a text file.
        The file should contain a column with the data you want to use for
        the analysis.
        """
        scan_data = np.loadtxt(file_path)
        return cls(scan_data, method=method, regions=regions, **kwargs)

    # --------- PROPERTIES --------- #
    @property
    def method(self):
        return self._method

    # --------- METHODS --------- #
    def permute_data(self, n_permutations=1000):
        """Permute the imaging data maintaining spatial autocorrelation for
        the cortical regions. The permutation is done using the 
        netneurotools Python package.

        :param int n_permutations: number of permutations.
        """
        _permuted = np.zeros((self.zscore_data.shape[0], n_permutations))
        if self._subcortical is not None:
            sub_permuted = np.array(
                [np.random.permutation(self._subcortical) for _ in
                 range(n_permutations)]
            ).reshape(7, n_permutations)
            _permuted[34:, :] = sub_permuted
        # Cortical
        # Annotation file for the Desikan-Killiany atlas in fs5
        annot_lh = Path(__file__).resolve().parent / "data/fsa5_lh_aparc.annot"
        annot_rh = Path(__file__).resolve().parent / "data/fsa5_rh_aparc.annot"
        # Get the parcel centroids of the Desikan-Killiany atlas
        parcel_centroids, parcel_hemi = freesurfer.find_parcel_centroids(
            lhannot=annot_lh,
            rhannot=annot_rh,
            version="fsaverage5",
            surf="sphere",
            method="surface",
        )
        # Mask the results to have only the left hemisphere
        left_hemi_mask = parcel_hemi == 0
        parcel_centroids, parcel_hemi = (
            parcel_centroids[left_hemi_mask],
            parcel_hemi[left_hemi_mask],
        )
        # Get the spin samples
        spins = stats.gen_spinsamples(
            parcel_centroids, parcel_hemi,
            n_rotate=n_permutations,
            method="vasa",
            seed=1234
        )
        cort_permuted = np.array(self._cortical[spins]).reshape(34,
                                                                n_permutations)
        _permuted[0:34, :] = cort_permuted
        self._permutations = _permuted

    # --------- RUN ANALYSIS --------- #
    def run(self):
        pass


class PLSAnalysis:
    def __init__(self, n_components: int, var: float):
        self.n_components = n_components
        self.var = self.check_var(var)
        if self.var is None and self.n_components is None:
            raise ValueError("You must specify either the variance or "
                             "the number of components for pls regression")
        self.components_var = np.zeros(15)
        self._p_val = np.zeros(self.n_components)
        self._r2 = np.zeros(self.n_components)
        self.gene_results = GeneResults("pls", n_components=self.n_components)

    @staticmethod
    def check_var(var: float):
        """ Check if the variance is between 0 and 1.

        :param float var: variance.
        """
        if var < 0:
            raise ValueError("The variance must be a positive number.")
        elif var > 1:
            raise ValueError("The variance must be a number between 0 and 1.")
        return var

    def set_coef(self, data, gene_exp):
        """Set the coefficients for the PLS regression. The function will
        estimate the variance or the number of components depending on the
        non missing parameter through a PLS regression with 15 components.

        """
        res = pls_regression(data, gene_exp,
                             n_components=15,
                             n_perm=0,
                             n_boot=0)
        explained_var = res.get("varexp")
        if self.var is None:
            self.var = np.cumsum(explained_var)[self.n_components-1]
        if self.n_components is None:
            dim = 1
            cumulative_var = np.cumsum(explained_var)
            while cumulative_var[dim - 1] < self.var:
                dim += 1
            self.n_components = dim
        self.components_var = explained_var

    @property
    def p_val(self):
        return self._p_val

    @property
    def r2(self):
        return self._r2

    @p_val.setter
    def p_val(self, p_val):
        self._p_val = p_val

    @r2.setter
    def r2(self, r2):
        self._r2 = r2


class CorrAnalysis:
    def __init__(self):
        self.gene_results = GeneResults("corr")


class GeneResults:
    def __init__(self, method, **kwargs):
        """Initialize the results of the analysis. Depending on the method
        used, the results will have underlying result classes, which account
        for the different analysis methods.

        :param str method: the method used for the analysis.
        :param kwargs: Additional parameters, for the initialisation. If the
        method is "pls" among the kwargs you *MUST* specify the number of
        components used, for the initialisation of the pls class.
        """
        self.method = method
        if self.method == "pls":
            self.results = PLSGenes(kwargs.get("n_components"))
        elif self.method == "corr":
            self.results = CorrGenes()
        else:
            raise ValueError(f"The method {method} is not supported.")


class PLSGenes:
    def __init__(self, n_components):
        self.n_components = n_components


class CorrGenes:
    pass
