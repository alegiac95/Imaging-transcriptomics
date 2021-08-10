from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from fpdf import FPDF

from .errors import CheckPath


class PDF(FPDF):
    """Class to generate a PDF report for the imaging-transcriptomics script.
    """
    def header(self):
        """The header will contain always the title."""
        self.rect(10, 10, 190, 280)
        self.line(10, 50, 200, 50)
        self.set_font("Helvetica", "B", 14)
        self.cell(w=0, h=15, align="C", txt="Imaging Transcriptomics Analysis Report", ln=True)

    def analysis_info(self, filename, date, filepath):
        """Info on the analysis performed. Information included are the name of
        the scan of the input, date of the analysis and the original path of the scan."""
        self.set_font("Courier", "", 10)
        self.cell(w=100, h=8, align="L", txt=f"  Scan Name: {filename}")
        self.cell(w=100, h=8, align="L", txt=f"  Date: {date}", ln=True)
        self.cell(w=100, h=10, align="L", txt=f"  File Path: {filepath}")

    def pls_regression(self, path_plots):
        """Include the plots of the pls components."""
        self.ln(20)
        self.set_font("Helvetica", "BU", 12)
        self.cell(w=0, h=10, align="L", txt="-PLS Regression")
        self.ln(10)
        self.image(Path(path_plots) / "individual_variance.png", x=12, w=120)
        self.image(Path(path_plots) / "cumulative_variance.png", x=12, w=120)


@CheckPath
def make_folder(path, folder_name: str):
    """Create a folder based on path and name to use.

    The function checks if in the directory there are already folders with the same name. If this is the case the new
    folder will have a trialing underscore with a progressive number (e.g., _2)

    :param path: path where to save the folder.
    :param folder_name: name of the folder to save

    :return folder_made: path of the folder created.
    """
    n_folds = len(list(Path(path).absolute().glob(f"{folder_name}_*")))
    if n_folds != 0:
        folder_name = f"{folder_name}_{n_folds+1}"
    folder_made = Path(path).absolute() / folder_name
    folder_made.mkdir()
    return folder_made


@CheckPath
def make_plots(path, limit_x, data_y):
    """Generate the plots for the explained variance by each component, one with the cumulative sum and one with the
    variance explained by each individual component.

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
    plt.plot(limit_x, varexp[limit_x - 1], 'o', color="red")
    plt.vlines(limit_x, varexp[0] - 10, varexp[limit_x - 1], colors="lightgrey", linestyles="dashed")
    plt.hlines(varexp[limit_x - 1], 0, limit_x, colors="lightgrey", linestyles="dashed")
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
    plt.xlabel("PLS components")
    plt.ylabel("Variance (%)")
    plt.savefig(path / "individual_variance.png", dpi=1200)
    plt.close()
    return


def create_pdf():
    pass
