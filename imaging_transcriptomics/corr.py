import numpy as np
from pathlib import Path
from scipy.stats import spearmanr
import pandas as pd
import logging
from statsmodels.stats.multitest import multipletests
from multiprocessing import Pool
from functools import partial
from itertools import product
from collections import OrderedDict
import gseapy
from gseapy.plot import gseaplot
from .genes import GeneResults, CorrGenes
from .inputs import get_geneset
import yaml
import logging
import logging.config

np.random.seed(1234)

cfg_file_path = Path(__file__).parent / "log_config.yaml"
with open(cfg_file_path, "r") as config_file:
    log_cfg = yaml.safe_load(config_file.read())

logging.config.dictConfig(log_cfg)
logger = logging.getLogger("genes")
logger.setLevel(logging.DEBUG)


# --------- CORR ANALYSIS --------- #
def _spearman_op(idx, permuted, y):  # pragma: no cover, used for multiprocess
    """Wrapper for the spearman correlation function, to be used for
    parallel computation."""
    return spearmanr(permuted[:, idx[0]], y[:, idx[1]])[0]


class CorrAnalysis:
    """Class to store all the results of the correlation analysis.
    The only field contained is the a GeneResults object, storing all
    information about the correlation analysis.
    """
    def __init__(self):
        self.gene_results = GeneResults("corr")
        self.gsea_res = None

    # --------- COMPUTE FUNCTIONS --------- #
    def bootstrap_correlation(self, imaging_data, permuted_imaging,
                              gene_exp, gene_labels, n_cpu=4):
        """Perform bootstrapping on the correlation.

        The function first calculates the correlation between the imaging
        vector and each of the genes. Then, it performs 1000 bootstrapping
        iterations of the same correaltion only using the permuted imaging
        data. The function then calls the methods of the GeneResults class
        to order the genes and calculate the p-values.

        :param np.ndarray imaging_data: the imaging data. This vector
        represents the original data to be correlated, and can have either
        length 41 or 34, depending on whether the subcortical regions want
        to be included in the analysis.
        :param np.ndarray permuted_imaging: the imaging data permuted. This
        is a matrix with shape (n_imaging, 1000), where `n_imaging`
        represents the length of the imaging vector.
        :param np.ndarray gene_exp: the gene expression data.
        :param list gene_labels: the labels of the genes.
        :param int n_cpu: the number of CPUs to use for the parallelization
        """
        assert isinstance(self.gene_results.results, CorrGenes)
        logger.info("Calculating correlation on original data.")
        for i in range(self.gene_results.n_genes):
            self.gene_results.results.corr[0, i], _ = spearmanr(
                imaging_data, gene_exp[:, i])
        pool = Pool(n_cpu)
        _ind = product(range(1000), range(self.gene_results.n_genes))

        logger.info("Calculating correlation on permuted data.")
        for ind, res in enumerate(pool.imap(partial(_spearman_op,
                                                    permuted=permuted_imaging,
                                                    y=gene_exp),
                                            _ind, chunksize=15633)):
            self.gene_results.results.boot_corr.flat[ind] = res
        self.gene_results.results.genes = gene_labels
        self.gene_results.results.sort_genes()
        self.gene_results.results.compute_pval()
        return

    def gsea(self, gene_set="lake", outdir=None):  # pragma: no cover,
        # long to process (t > 1 h) - TESTED on run
        """Perform GSEA on the correlation.

        The function runs a first gsea with the data and then runs the same
        analysis using the permuted data. The analysis of the permuted data
        is used ot calculate the p-values for the Enrichment scores (ES).

        :param str gene_set: the gene set to use for the analysis.
        :param str outdir: the directory where to save the results.
        """
        assert isinstance(self.gene_results.results, CorrGenes)
        logger.info("Performing GSEA.")
        gene_set = get_geneset(gene_set)
        # prepare the gene_list as a list of strings
        gene_list = [
            gene for gene in self.gene_results.results.genes[:, 0].tolist()
        ]
        # perform the GSEA on original results
        rnk = pd.DataFrame(
                zip(gene_list, self.gene_results.results.corr[0, :]))
        gsea_results = gseapy.prerank(rnk, gene_set,
                                      outdir=None,
                                      permutation_num=1000,
                                      seed=1234)
        _origin_es = gsea_results.res2d.es.to_numpy()
        _boot_es = np.zeros((_origin_es.shape[0], 1000))
        # perform the GSEA on the permutations
        for i in range(1000):
            rnk = pd.DataFrame(
                zip(gene_list, self.gene_results.results.boot_corr[:, i])
            )
            _gsea_res = gseapy.prerank(rnk, gene_set,
                                       permutation_num=1,
                                       no_plot=True,
                                       outdir=None,
                                       seed=1234)
            _boot_es[:, i] = _gsea_res.res2d.es.to_numpy()
        # calculate the p-value
        _p_val = np.zeros((_origin_es.shape[0],))
        for i in range(_origin_es.shape[0]):
            _p_val[i] = np.sum(_boot_es[i, :] >= _origin_es[i]) / 1000
        # calculate the p-value corrected
        _, _p_corr, _, _ = multipletests(_p_val, method='fdr_bh')
        # Prepare data to save
        _out_data = OrderedDict()
        _out_data["Term"] = gsea_results.res2d.axes[0].to_list()
        _out_data["es"] = gsea_results.res2d.values[:, 0]
        _out_data["nes"] = gsea_results.res2d.values[:, 1]
        _out_data["p_val"] = _p_val
        _out_data["fdr"] = _p_corr
        _out_data["genest_size"] = gsea_results.res2d.values[:, 4]
        _out_data["matched_size"] = gsea_results.res2d.values[:, 5]
        _out_data["matched_genes"] = gsea_results.res2d.values[:, 6]
        _out_data["ledge_genes"] = gsea_results.res2d.values[:, 7]
        out_df = pd.DataFrame.from_dict(_out_data)
        if outdir is not None:
            logger.info("Saving GSEA results.")
            outdir = Path(outdir)
            assert outdir.exists()
            out_df.to_csv(outdir / "gsea_corr_results.tsv", index=False,
                          sep="\t")
            for i in range(len(gsea_results.res2d.index)):
                term = gsea_results.res2d.index[i]
                gsea_results.results[term]["pval"] = _p_val[i]
                gsea_results.results[term]["fdr"] = _p_corr[i]
                gseaplot(rank_metric=gsea_results.ranking,
                         term=term,
                         **gsea_results.results[term],
                         ofname=f"{outdir}/{term}_corr_prerank.pdf")

    # --------- SAVE FUNCTIONS --------- #
    def save_results(self, outdir):  # pragma: no cover, simply saves stuff
        """Save the results to a file."""
        outdir = Path(outdir)
        assert outdir.exists()
        # Create the data to save
        logger.info("Saving results.")
        data_to_save = zip(self.gene_results.results.genes[:, 0],
                           self.gene_results.results.corr[0, :],
                           self.gene_results.results.pval[0, :],
                           self.gene_results.results.pval_corr[0, :])
        # Save the data
        pd.DataFrame(
            data_to_save, columns=["Gene", "Corr", "Pval", "Pval_corr"]
        ).to_csv(outdir / "corr_genes.tsv", index=False, sep='\t')
        return
