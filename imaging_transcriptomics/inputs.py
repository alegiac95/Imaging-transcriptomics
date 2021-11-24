from pathlib import Path
import logging
import logging.config
import yaml

import nibabel as nib
import numpy as np
import pandas as pd
from scipy.stats import zscore

from .errors import CheckShape, CheckPath, CheckExtension, CheckVariance


cfg_file_path = Path(__file__).parent / "log_config.yaml"
with open(cfg_file_path, "r") as config_file:
    log_cfg = yaml.safe_load(config_file.read())

logging.config.dictConfig(log_cfg)
logger = logging.getLogger("inputs")
logger.setLevel(logging.DEBUG)


# Imaging data
@CheckPath
@CheckExtension
def read_scan(path):
    """Return the imaging file associated to the input scan.

    Uses the `Nibabel <https://nipy.org/nibabel/>`_ to read the input imaging file and get the voxel data.

    :param str path: path of the imaging file to analyze.

    :return: Numpy matrix with the voxel of the input scan.
    """
    logger.debug("Reading scan: %s", path)
    data = nib.load(Path(path)).get_fdata()
    return data


@CheckShape
def extract_average(imaging_matrix):
    """Extract the average value of the ROIs from the imaging scan.

    The values are extracted from the left hemisphere only, since the data of the Allen Human Brain Atlas are
    available for that hemisphere only for all donors.

    :param imaging_matrix: matrix with the voxels of the image.

    :return: numpy array with the average value from 41 brain regions.
    """
    n_regions = 41
    logger.debug("Extracting average from scan.")
    atlas_data = nib.load(
        Path(__file__).resolve().parent
        / "data"
        / "atlas-desikankilliany_1mm_MNI152.nii.gz"
    ).get_fdata()
    data = np.zeros(n_regions)
    for i in range(1, n_regions + 1):
        data[i - 1] = np.mean(imaging_matrix[np.where(atlas_data == i)])
    logger.debug("Extracted values are: \n%s", data)
    return np.array(data)


# Other parameters
@CheckVariance
def get_components(target_variance, explained_var):
    """Return the optimal number of PLS components to explain a desired amount of the total variance.

    :param float target_variance: total explained variance
    :param np.array explained_var: array with explained variance by each individual pls component.

    :return int dim: optimal number of components to extract the desired explained variance.
    """
    dim = 1
    cumulative_var = np.cumsum(explained_var)
    while cumulative_var[dim - 1] < target_variance:
        dim += 1
    logger.debug(
        "Extracted variance is %s with %s components", cumulative_var[dim - 1], dim
    )
    return dim


# Gene expression data
def load_gene_expression():
    """Return matrix with gene expression data.

    The data have been previously normalised and are available in the  ``data`` sub-folder.

    :return: numpy array with the gene expression data.
    """
    logger.debug("Loading gene_expression data.")
    expression_file_path = (
        Path(__file__).resolve().parent / "data" / "gene_expression_data.csv"
    )
    expression_data = pd.read_csv(expression_file_path, sep=",")
    my_data_x = expression_data.iloc[0:41, 2:].to_numpy()
    return zscore(my_data_x, ddof=1)


def load_gene_labels():
    """Return an array with the gene labels.
    The gene labels are available in the ``data`` sub-folder.

    :return: numpy array with the labels of the genes.
    """
    logger.debug("Loading gene labels.")
    genes_labels_path = (
        Path(__file__).resolve().parent / "data" / "gene_expression_labels.txt"
    )
    return pd.read_fwf(genes_labels_path, header=None).to_numpy()
