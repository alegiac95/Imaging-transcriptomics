from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from fpdf import FPDF

from .errors import CheckPath
from .genes import GeneResults


class PDF(FPDF):
    """Class to generate a PDF report for the imaging-transcriptomics script."""

    def header(self):
        """The header will contain always the title."""
        self.rect(10, 10, 190, 280)
        self.line(10, 50, 200, 50)
        self.set_font("Helvetica", "B", 14)
        self.cell(
            w=0, h=15, align="C", txt="Imaging Transcriptomics Analysis Report", ln=True
        )

    def analysis_info(self, filename, date, filepath):
        """Info on the analysis performed. Information included are the name of
        the scan of the input, date of the analysis and the original path of the scan.

        :param filename: name of the scan used for analysis.
        :param date: date when the analysis was performed.
        :param filepath: absolute path of the imaging scan used for analysis.
        """
        self.set_font("Courier", "", 10)
        self.cell(w=100, h=8, align="L", txt=f"  Scan Name: {filename}")
        self.cell(w=100, h=8, align="L", txt=f"  Date: {date}", ln=True)
        self.cell(w=100, h=10, align="L", txt=f"  File Path: {filepath}")

    def pls_regression(self, path_plots):
        """Include the plots of the pls components.

        :param path_plots: path where the .png plots are located.
        """
        self.ln(20)
        self.set_font("Helvetica", "BU", 12)
        self.cell(w=0, h=10, align="L", txt="-PLS Regression")
        self.ln(10)
        self.image(Path(path_plots) / "individual_variance.png", x=15, w=120)
        self.image(Path(path_plots) / "cumulative_variance.png", x=15, w=120)

    def reproducibility_line(self, cli_commands, version):
        """Create a string with the command used from the command line to run the analysis.

        :param cli_commands: commands given in the cli to run the analysis.
        :param version: version of the software used to run the analysis.
        """
        pass


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
    plt.plot(limit_x, varexp[limit_x - 1], "o", color="red")
    plt.vlines(
        limit_x,
        varexp[0] - 10,
        varexp[limit_x - 1],
        colors="lightgrey",
        linestyles="dashed",
    )
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


def create_pdf(filepath, save_dir):
    """Create a PDF report.

    Creates a report for the imaging transcriptomics analysis with some plots and details of what was run by the user.

    :param filepath: path of the image used for analysis.
    :param save_dir: path where to save the report (same as where the plots are located).
    """
    if not isinstance(filepath, Path):
        filepath = Path(filepath)
    analysis_date = datetime.now().strftime("%d-%m-%Y")

    report = PDF(orientation="P", unit="mm", format="A4")
    report.add_page()
    report.analysis_info(filename=filepath.name, date=analysis_date, filepath=filepath)
    report.pls_regression(path_plots=save_dir)
    report.output(save_dir / "Report.pdf", "F")


def create_csv(analysis_results, n_comp, save_dir):
    """Create .csv files for the results of each component.

    The function creates a different csv file for each of the components used for PLS regression (if 2 components are
    used the files 'PLS1.csv' and 'PLS2.csv' will be created.

    :param analysis_results: GeneResults data structure with the results of bootstrapping.
    :param int n_comp: number of components used in the regression.
    :param save_dir: path where the output will be saved.
    :return:
    """
    if not isinstance(analysis_results, GeneResults):
        raise TypeError("The data are not of the GeneResults class.")
    for i in range(n_comp):
        data = np.vstack(
            (
                np.array(analysis_results.boot_results.pls_genes[i].reshape(1, 15633)),
                np.array(analysis_results.boot_results.z_scores[i]),
                np.array(analysis_results.boot_results.pval[i]),
                np.array(analysis_results.boot_results.pval_corrected[i]),
            )
        ).T
        data = pd.DataFrame(data, columns=["Gene ID", "Z", "p", "p corrected"])
        data.to_csv(save_dir / f"PLS{i+1}.csv", index=False)


def create_corr_csv(analysis_results, save_dir):
    """Create a csv file for the correlation coefficients.

    :param analysis_results: GeneResults data structure with the results of bootstrapping.
    :param save_dir: path where the output will be saved.
    :return:
    """
    if not isinstance(analysis_results, GeneResults):
        raise TypeError("The data are not of the GeneResults class.")
    if not isinstance(save_dir, Path):
        save_dir = Path(save_dir)
    data = np.vstack(
        (
            np.array(analysis_results.boot_results.pls_genes.reshape(1, 15633)),
            np.array(analysis_results.original_results.pls_weights),
            np.array(analysis_results.boot_results.pval),
            np.array(analysis_results.boot_results.pval_corrected),
        )
    ).T
    data = pd.DataFrame(data, columns=["Gene ID", "Correlation coefficient", "p", "p corrected"])
    data.to_csv(save_dir / "Correlation_coefficients.csv", index=False)
