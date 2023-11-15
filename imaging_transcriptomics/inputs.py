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

    Uses the `Nibabel <https://nipy.org/nibabel/>`_
    to read the input imaging file and get the voxel data.

    :param str path: path of the imaging file to analyze.
    :return: Numpy matrix with the voxel of the input scan.
    """
    logger.debug("Reading scan: %s", path)
    return nib.load(Path(path)).get_fdata()


@CheckShape
def get_vox_size(imaging_matrix):
    """
    Get the voxel dimensions in mm from the matrix dimensions (i.e., (91,
    109, 91) = 2mm
    :param imaging_matrix:
    :return str mm:
    """
    if imaging_matrix.shape == (91, 109, 91):
        mm = "2mm"
    elif imaging_matrix.shape == (182, 218, 182):
        mm = "1mm"
    else:
        raise ValueError("This voxel dimension is not supported!")
    return mm


# Atlas data
def load_atlas_imaging(atlas_name="DK", vox_dim="1mm"):
    """
    Load the name and path of the atlas used for the analysis.

    :param vox_dim: voxel dimension of an image to analyse.
    :param atlas_name: Name of the atlas to use for analysis. Available
    atlases, so far, are: Desikan-Killiany ("DK_1mm"), Shaefer 100 (
    "Shaefer_100"),
    Shaefer 200 ("Shaefer_200"), Shaefer 400 ("Shaefer_400").
    :return:
    """
    # Check if the atlas specified is one of the available atlases
    atlas_list = ["DK", "Schaefer_100", "Schaefer_200", "Schaefer_400"]
    atlas_base_path = Path(__file__).resolve().parent / "data" / "atlases"
    if atlas_name not in atlas_list:
        raise FileExistsError(f"The atlas {atlas_name} is not a valid "
                              f"atlas!")
    elif atlas_name == "DK":
        # Desikan-Killiany atlas
        n_regions = 41
        atlas_path = atlas_base_path / "DK"
        atlas_path = atlas_path / "atlas-DK_1mm.nii.gz" if vox_dim == "1mm" \
            else atlas_path / "atlas-DK_2mm.nii.gz"
    elif atlas_name == "Schaefer_100":
        n_regions = 50
        atlas_path = atlas_base_path / "Schaefer_100"
        atlas_path = atlas_path / \
                     "atlas-Schaefer_100_1mm.nii.gz" if vox_dim == "1mm" \
            else atlas_path / "atlas-Schaefer_100_2mm.nii.gz"
    elif atlas_name == "Schaefer_200":
        n_regions = 100
        atlas_path = atlas_base_path / "Schaefer_200"
    elif atlas_name == "Schaefer_400":
        n_regions = 200
        atlas_path = atlas_base_path / "Schaefer_400"

    return n_regions, atlas_path


@CheckShape
def extract_average(imaging_matrix, atlas="DK"):
    """Extract the average value of the ROIs from the imaging scan.

    The values are extracted from the left hemisphere only, since the data of
    the Allen Human Brain Atlas are available for that hemisphere only for
    all donors.

    :param atlas: atlas to use for the parcellation of the scan
    :param imaging_matrix: matrix with the voxels of the image.

    :return: numpy array with the average value from 41 brain regions.
    """
    vox_dim = get_vox_size(imaging_matrix)
    # Load atlas
    n_regions, atlas_path = load_atlas_imaging(atlas, vox_dim=vox_dim)
    logger.debug(f"Extracting average from {atlas} of scan")
    atlas_data = nib.load(atlas_path).get_fdata()
    # create empty array to store the regional average values
    data = np.zeros(n_regions)
    for i in range(1, n_regions + 1):
        data[i - 1] = np.nanmean(imaging_matrix[np.where(atlas_data == i)])
    logger.debug(f"Extracted values are: {data}")
    return np.array(data)


def get_annot_files(atlas="DK"):
    atlas_base_path = Path(__file__).resolve().parent / "data" / "atlases" / \
                      atlas
    lh_annot = list(atlas_base_path.glob(f"atlas-{atlas}*_lh_aparc.annot"))[0]
    rh_annot = list(atlas_base_path.glob(f"atlas-{atlas}*_rh_aparc.annot"))[0]
    return str(lh_annot), str(rh_annot)


# Gene expression data
def load_gene_expression(regions="cort+sub", atlas="DK"):
    """Return matrix with gene expression data.

    The data have been previously normalised and are available in the
    ``data/atlases`` sub-folder.

    :param regions: String indicating whether the regions of only cortical
    or cortical and subcortical are used.
    :param atlas:
    :return: numpy array with the gene expression data.
    """
    if regions not in ["cort+sub", "cort", "all"]:
        raise ValueError("The regions must be either 'cort+sub', 'cort' or "
                         "'all'.")
    logger.debug("Loading gene_expression data.")
    # TODO: change file loadings based on atlas, and retrieve automatically
    #  the gene labels from header.
    expression_file_path = (
            Path(__file__).resolve().parent / f"data/atlases/{atlas}/"
                                              f"atlas-{atlas}_gene_expression_data.csv")
    expression_data = pd.read_csv(expression_file_path, sep=",")
    if atlas == "DK":
        my_data_x = expression_data.iloc[0:41, 2:].to_numpy()
        my_data = zscore(my_data_x, ddof=1)
        if regions in ["cort+sub", "all"]:
            my_data = my_data
        elif regions == "cort":
            my_data = my_data[:34, :]
    elif atlas == "Schaefer_100":
        my_data_x = expression_data.iloc[0:50, 2:].to_numpy()
        my_data = zscore(my_data_x, ddof=1)
        if regions in ["cort+sub", "all"]:
            my_data = my_data
        elif regions == "cort":
            my_data = my_data[:34, :]
    return my_data


def load_gene_labels(atlas="DK"):
    """Return an array with the gene labels.
    The gene labels are available in the ``data`` sub-folder.

    :return: numpy array with the labels of the genes.
    """
    logger.debug("Loading gene labels.")
    expression_file_path = (
            Path(__file__).resolve().parent /
            f"data/atlases/{atlas}/atlas-{atlas}_gene_expression_data.csv")
    expression_labels = pd.read_csv(expression_file_path, sep=",")
    return np.array(expression_labels.columns[2:]).reshape(
        (expression_labels.columns[2:].shape[0], 1))


def get_geneset(gene_set: str):
    """Returns the path to the geneset file, if it is not one of those
    defined in gseapy.

    :param str gene_set: the name of the geneset. If the geneset is "lake"
    or "pooled", the genesets are provided with the package data, otherwise
    the gene sets are from the gseapy package.
    """
    if gene_set.lower() == "lake":
        return str(Path(__file__).parent / "data" / "geneset_LAKE.gmt")
    elif gene_set.lower() == "pooled":
        return str(Path(__file__).parent / "data" /
                   "geneset_Pooled.gmt")
    else:
        return gene_set
