from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import nibabel as nib
from fpdf import FPDF
from enigmatoolbox.utils.parcellation import parcel_to_surface
from enigmatoolbox.plotting.surface_plotting import plot_cortical, \
    plot_subcortical

from .errors import CheckPath
from .genes import GeneResults
from .transcriptomics import ImagingTranscriptomics


@CheckPath
def make_folder(path, folder_name: str):
    """Create a folder based on path and name to use.

    The function checks if in the directory there are already folders with the same name. If this is the case the new
    folder will have a trialing underscore with a progressive number (e.g., _2)

    :param path: path where to save the folder.
    :param folder_name: name of the folder to save

    :return folder_made: path of the folder created.
    """
    n_folds = len(list(Path(path).absolute().glob(f"{folder_name}*")))
    if n_folds != 0:
        folder_name = f"{folder_name}_{n_folds}"
    folder_made = Path(path).absolute() / folder_name
    folder_made.mkdir()
    return folder_made


@CheckPath
def make_plots(path, limit_x, data_y):
    """Generate the plots for the explained variance by each component,
    one with the cumulative sum and one with the variance explained by each
    individual component.

    :param path: path where the plots will be saved.
    :param limit_x: number of PLS components.
    :param data_y: data to plot on the y axis.
    :return:
    """
    if not isinstance(path, Path):
        path = Path(path)
    varexp = 100 * np.cumsum(data_y)

    # Plot cumulative percentage variance
    plt.plot(range(1, 16), varexp, marker="o", color="sandybrown")
    plt.plot(limit_x, varexp[limit_x - 1], "o", color="red")
    plt.vlines(
        limit_x,
        varexp[0] - 10,
        varexp[limit_x - 1],
        colors="lightgrey",
        linestyles="dashed",
    )
    plt.hlines(varexp[limit_x - 1], 0, limit_x, colors="lightgrey",
               linestyles="dashed")
    plt.title("Cumulative variance explained by PLS components")
    plt.ylabel("Total explained variance (%)")
    plt.xlabel("Number of PLS components")
    plt.xlim(0, 15)
    plt.ylim(varexp[0] - 10, 105)
    plt.grid(True)
    plt.savefig(path / "cumulative_variance.png", dpi=1200)
    plt.close()

    # Plot individual explained variance by each component
    plt.bar(range(1, 16), 100 * data_y, color="sandybrown")
    for index, value in enumerate(data_y):
        plt.text(index + 0.5, 100 * value, "{:.1f}".format(100 * value))
    plt.title("Individual variance explained by PLS components")
    plt.xlabel("PLS component")
    plt.ylabel("Variance (%)")
    plt.savefig(path / "individual_variance.png", dpi=1200)
    plt.close()
    return


def extract_from_atlas(image):
    """Extract the average signal from the image using the DK
    atlas.
    """
    # Atlas to use for parcellation
    atlas_image = Path(__file__).parent / "data" / \
                  "atlas-desikankilliany_1mm_MNI152.nii.gz"
    atlas_image = nib.load(atlas_image).get_fdata()
    max_region = int(atlas_image.max())

    new_order = pd.read_csv("alessio_order_enigma.csv")["Enigma_order"]
    new_order = new_order.to_numpy() - 1
    image = nib.load(Path(image)).get_fdata()
    e_values = np.zeros(max_region)
    for i in range(1, max_region):
        e_values[i] = np.mean(image[np.where(atlas_image==i)])
    return e_values[new_order]


def plot_brain(image):
    average_val = extract_from_atlas(image)
    avg_fs5 = parcel_to_surface(average_val, 'aparc_fsa5')
    plot_cortical(array_name=avg_fs5,
                  surface_name='fsa5',
                  size=(800, 400),
                  color_bar=True,
                  color_range=(-0.5, 0.5),
                  screenshot=True,
                  filename=f"cortical.png",
                  transparent_bg=True)
    plot_subcortical(array_name=average_val[67:],
                     size=(800, 400),
                     color_bar=True,
                     color_range=(-0.5, 0.5),
                     screenshot=True,
                     filename=f"subcortical.png",
                     transparent_bg=True)


def pls_components(data):
    """Return a string with the number of components, its cumulative
    variance and its p value from the analysis."""
    if not isinstance(data, ImagingTranscriptomics):
        raise TypeError("The data must be an ImagingTranscriptomics object.")
    n_components = data.analysis.n_components
    res_str = "PLS Component: R2 pval "
    for i in range(n_components):
        res_str += f"PLS {i+1}: " \
                   f"{data.analysis.r2[i]:.3f}" \
                   f" {data.analysis.p_val[i]:.4f}"
    return res_str


def make_pdf(transcriptomics_data, save_dir, name="report", scanname="", ):
    if not isinstance(transcriptomics_data, ImagingTranscriptomics):
        raise TypeError("The data must be an ImagingTranscriptomics object.")
    save_dir = Path(save_dir)
    WIDTH = 210  # mm
    HEIGHT = 297  # mm
    MARGIN = 5  # mm
    report = FPDF(orientation="P", unit="mm", format="A4")
    report.set_font("Arial", size=11)
    report.add_page()
    # Header
    report.image(str(Path(__file__).parent / "resources" / "header.png"),
                 x=0, y=0, w=WIDTH)
    # Info on the analysis
    analysis_date = datetime.now().strftime("%d-%m-%Y")
    report.ln(25)
    report.cell(WIDTH/2, txt=f"Scan name: {scanname}")
    report.cell(WIDTH/2, txt=f"Analysis date:{analysis_date}")
    # PLS plots
    if transcriptomics_data.method == "pls":
        report.ln(10)
        make_plots(save_dir,
                   transcriptomics_data.analysis.n_components,
                   transcriptomics_data.analysis.components_var)
        report.cell(WIDTH, txt="PLS Analysis")
        report.image(str(save_dir / "individual_variance.png"),
                     x=MARGIN, y=53, w=WIDTH/2)
        report.image(str(save_dir / "cumulative_variance.png"),
                     x=WIDTH/2 + MARGIN, y=53,  w=WIDTH/2)
        report.ln(90)
        report.cell(WIDTH, txt=pls_components(transcriptomics_data))
    # MAKE GENE RESULTS TABLE
    report.ln(80)
    report.cell(WIDTH, txt="For the gene results refer to the spreadsheet in "
                           "the report folder.")
    report.output(str(save_dir / f"{name}.pdf"), "F")
