import logging.config
import logging
import warnings
from pathlib import Path
import pickle

import numpy as np
import yaml
from scipy.stats import zscore

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from netneurotools import freesurfer, stats

from .inputs import load_gene_expression, load_gene_labels, \
    extract_average, read_scan
from .pls import PLSAnalysis
from .corr import CorrAnalysis

cfg_file_path = Path(__file__).parent / "log_config.yaml"
with open(cfg_file_path, "r") as config_file:
    log_cfg = yaml.safe_load(config_file.read())

logging.config.dictConfig(log_cfg)
logger = logging.getLogger("transcriptomics")
logger.setLevel(logging.DEBUG)


class ImagingTranscriptomics:
    # --------- INITIALIZATION --------- #
    def __init__(self,
                 scan_data,
                 regions="cort+sub",
                 method="pls",
                 **kwargs):
        """ImagingTranscriptomics class for imaging transcriptomics analysis.

        :param np.array scan_data: imaging scan data.
        :param str regions: regions to be used for analysis. These can be
        "cort+sub" (or "all") which will perform the analysis on the
        cortical and subcortical regions, or "cort" which will only perform
        the analysis on the cortical regions.
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
            self._cortical = self.zscore_data if scan_data.shape == (34,) else\
                self.zscore_data[:34]
            self._subcortical = None
        self._regions = regions
        self.scan_data = scan_data
        self.gene_expression = load_gene_expression(self._regions)
        self.gene_labels = load_gene_labels()
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
                if self._regions == "all" or self._regions == "cort+sub":
                    self.analysis = PLSAnalysis(self.zscore_data,
                                                self.gene_expression,
                                                kwargs.get("n_components"),
                                                kwargs.get("var"))
                else:
                    self.analysis = PLSAnalysis(self._cortical,
                                                self.gene_expression,
                                                kwargs.get("n_components"),
                                                kwargs.get("var"))
        elif self._method == "corr":
            self.analysis = CorrAnalysis()
        self._permutations = None
        self._permutation_ind = None

    @classmethod
    def from_scan(cls, scan_path, method="pls", regions="cort+sub", **kwargs):
        """Initialise an ImagingTranscriptomics object from a NIfTI scan.
        The extracted data corresponds to the average of the ROIs in the DK
        atlas.

        :param str scan_path: path to the NIfTI scan.
        :param str method: method to run the analysis, can be either "pls"
        for pls regression or "corr" cor simple correlation analysis.
        :param str regions: regions to be used for analysis. These can be
        "cort+sub" (or "all") which will perform the analysis on the
        cortical and subcortical regions, or "cort" which will only perform
        the analysis on the cortical regions.
        :param kwargs: additional arguments for the method. This include:
            * "n_components": number of components for pls regression.
            * "var": variance explained for pls regression.
            * "n_permutations": number of permutations for permutation test.
        :return: ImagingTranscriptomics object.
        """
        if not Path(scan_path).exists():
            raise FileNotFoundError("The specified scan file does not exist.")
        if method not in ["pls", "corr"]:
            raise ValueError(
                "The method must be either pls or corr."
                "Please choose either pls or corr to run the analysis."
            )
        if regions not in ["cort+sub", "cort", "all"]:
            raise ValueError(
                "The regions must be either cort+sub, cort or all."
                "Please choose one of these to run the analysis."
            )
        scan_data = extract_average(read_scan(scan_path))
        return cls(scan_data, method=method, regions=regions, **kwargs)

    @classmethod
    def from_file(cls, file_path, method="pls", regions="cort+sub", **kwargs):
        """Initialise an ImagingTranscriptomics object from a text file.
        The file should contain a column with the data you want to use for
        the analysis.

        :param str file_path: path to the text file. The text file should
        contain a column with the data you want to use for the analysis.
        :param str method: method to run the analysis, can be either "pls"
        for pls regression or "corr" cor simple correlation analysis.
        :param str regions: regions to be used for analysis. These can be
        "cort+sub" (or "all") which will perform the analysis on the
        cortical and subcortical regions, or "cort" which will only perform
        the analysis on the cortical regions.
        :param kwargs: additional arguments for the method. This include:
            * "n_components": number of components for pls regression.
            * "var": variance explained for pls regression.
            * "n_permutations": number of permutations for permutation test.
        :return: ImagingTranscriptomics object.
        """
        if not Path(file_path).exists():
            raise FileNotFoundError("The specified file does not exist.")
        if not Path(file_path).is_file():
            raise ValueError(f"{file_path} is not a file.")
        scan_data = np.loadtxt(file_path)
        if scan_data.shape[0] > 41:
            scan_data = scan_data[:41]
        return cls(scan_data, method=method, regions=regions, **kwargs)

    # --------- PROPERTIES --------- #
    @property
    def method(self):  # pragma: no cover, simply returns stuff
        return self._method

    @property
    def gene_results(self):  # pragma: no cover, simply returns stuff
        return self.analysis.gene_results

    # --------- METHODS --------- #
    def permute_data(self, n_permutations=1000):
        """Permute the imaging data maintaining spatial autocorrelation for
        the cortical regions. The permutation is done using the
        netneurotools Python package.

        :param int n_permutations: number of permutations.

        """
        logger.info("Permuting data.")
        _permuted = np.zeros((self.zscore_data.shape[0], n_permutations))
        _perm_indexes = np.zeros((self.zscore_data.shape[0],
                                  n_permutations), dtype=np.int32)
        # Calculate the permutations on the subcortical regions.
        if self._subcortical is not None:
            sub_permuted = np.zeros((self._subcortical.shape[0],
                                     n_permutations))
            for i in range(n_permutations):
                sub_resample = np.random.choice(7, size=7)
                _perm_indexes[34:, i] = sub_resample + 34  # keep into
                # account the shift of the subcortical given by the cortical
                # regions.
                sub_permuted[:, i] = self._subcortical[sub_resample]
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
        _perm_indexes[:34, :] = spins
        _permuted[0:34, :] = cort_permuted
        self._permutations = _permuted
        self._permutation_ind = _perm_indexes
        return

    def _make_output_dir(self, output_dir, name=""):
        """Create the output directory if it does not exist.

        :param str output_dir: path to the output directory.
        :param str name: name of the output directory, if not provided
        disregard.
        """
        outdir = Path(output_dir) / f"Imt_{name}_{self.method}"
        outdir.mkdir(exist_ok=True)
        return outdir

    def _save_object(self, outdir, name):
        """Save the object as a pickle file.

        :param str outdir: path to the output directory.
        :param str name: name of the output file.
        """
        outfile = Path(outdir) / f"{name}.pkl"
        with open(outfile, "wb") as f:
            pickle.dump(self, f)
        return

    # --------- RUN ANALYSIS --------- #
    def gsea(self, outdir=None, gene_set="lake", gene_limit=500):

        if self.method == "corr":
            self.analysis.gsea(gene_set=gene_set, outdir=outdir,
                               gene_limit=gene_limit)
        elif self.method == "pls":
            self.gene_results.results.gsea(gene_set=gene_set,
                                           outdir=outdir,
                                           gene_limit=gene_limit)

    def run(self, outdir=None, scan_name="", gsea=True,
            gene_set="lake", save_res=True, n_cpu=4,
            gene_limit=500):  # pragma: no cover
        """Method to run the imaging transcriptomics analysis.

        :param str outdir: path to the output directory, if not provided the
        results will be saved in the current directory.
        :param str scan_name: name of the scan, if not provided the name will
        be ignored. Is used only to create the output directory.
        :param bool gsea: if True, run the GSEA analysis, if False the GSEA
        analysis is skipped.
        :param str gene_set: gene set to use for the GSEA analysis.
        :param bool save_res: if True, save the results in a directory,
        if False the results are not saved.
        :param int n_cpu: number of CPUs to use for the analysis (only for
        correlation analysis).
        :param int gene_limit: number of genes to use for the GSEA analysis.
        """
        logger.info(f"Running the {self.method} analysis.")
        # Create the permuted data matrix
        self.permute_data()
        # Check if the ouput directory is provided and create the output folder
        if save_res:
            if outdir is not None:
                outdir = self._make_output_dir(outdir, name=scan_name)
            else:
                outdir = self._make_output_dir(Path.cwd(), name=scan_name)
        # Run the analysis
        # CORRELATION
        if self._method == "corr":
            # Select the data or slice of data
            if self._regions == "cort":
                _d = self._cortical
                if self._permutations.shape[0] == 41:
                    _d_perm = self._permutations[0:34, :]
                else:
                    _d_perm = self._permutations
            elif self._regions == "cort+sub" or self._regions == "all":
                _d = self.zscore_data
                _d_perm = self._permutations
            self.analysis.bootstrap_correlation(_d, _d_perm,
                                                self.gene_expression,
                                                self.gene_labels,
                                                n_cpu=n_cpu)
            if save_res:
                self._save_object(outdir, f"{self.method}_analysis")
                self.analysis.save_results(outdir=outdir)
            if gsea:
                self.gsea(gene_set=gene_set, outdir=outdir)
        # PLS
        elif self._method == "pls":
            # Select the data or slice of data
            if self._regions == "cort":
                _d = self._cortical
                if self._permutations.shape[0] == 41:
                    _d_perm = self._permutations[0:34, :]
                else:
                    _d_perm = self._permutations
            elif self._regions == "cort+sub" or self._regions == "all":
                _d = self.zscore_data
                _d_perm = self._permutations
            assert isinstance(self.analysis, PLSAnalysis)
            self.analysis.boot_pls(_d, _d_perm, self.gene_expression)
            if self._regions == "cort":
                _orig = self.scan_data if self.scan_data.shape[0] == 34 else \
                    self.scan_data[0:34, :]
            elif self._regions == "cort+sub" or self._regions == "all":
                _orig = self.scan_data
            self.gene_results.results.boot_genes(_d,
                                                 _d_perm,
                                                 _orig,
                                                 self.gene_expression,
                                                 self.gene_labels)
            self.gene_results.results.compute()
            if save_res:
                self._save_object(outdir, f"{self.method}_analysis")
                self.analysis.save_results(outdir=outdir)
            if gsea:
                self.gsea(gene_set=gene_set, outdir=outdir,
                          gene_limit=gene_limit)
