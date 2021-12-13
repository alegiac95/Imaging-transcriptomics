import numpy as np
from pyls import pls_regression

from tqdm import tqdm
from .genes import GeneResults

np.random.seed(1234)


# --------- PLS ANALYSIS --------- #
class PLSAnalysis:
    """Class for performing PLS regression on the imaging data.
    This class contains fields for
    """
    def __init__(self, imaging_data, gene_exp,  n_components: int, var: float):
        self.var, self.n_components, self.components_var = self.set_coef(
            imaging_data, gene_exp, n_components=n_components, var=var)
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

    @staticmethod
    def set_coef(data, gene_exp, var=None, n_components=None):
        """Set the coefficients for the PLS regression. The function will
        estimate the variance or the number of components depending on the
        non missing parameter through a PLS regression with 15 components.

        :param data: imaging data.
        :param gene_exp: gene expression data.
        :param var: variance.
        :param n_components: number of components.
        """
        res = pls_regression(gene_exp, data.reshape(data.shape[0], 1),
                             n_components=15,
                             n_perm=0,
                             n_boot=0)
        explained_var = res.get("varexp")
        if var is None:
            var = np.cumsum(explained_var)[n_components-1]
        if n_components is None:
            dim = 1
            cumulative_var = np.cumsum(explained_var)
            while cumulative_var[dim - 1] < var:
                dim += 1
            n_components = dim
        return var, n_components, explained_var

    @property
    def p_val(self):
        """P-value of the pls regression given by the bootstrapping on the
        PLS components."""
        return self._p_val

    @property
    def r2(self):
        """R2 of the pls regression given by the bootstrapping on the PLS
        components"""
        return self._r2

    @p_val.setter
    def p_val(self, p_val):
        self._p_val = p_val

    @r2.setter
    def r2(self, r2):
        self._r2 = r2

    def boot_pls(self, imaging_data, permuted_imaging, gene_exp):
        """Bootstrapping on the PLS components.

        :param imaging_data: imaging data. Allows the user to specify the
        data to use (e.g., with only cortical regions this can be only the
        cortical vector, other wise the whole data).
        :param permuted_imaging: imaging data permuted. Allows the user to
        specify the data to use (e.g., with only cortical regions this can be
        only the cortical vector, other wise the whole data).
        :param gene_exp: gene expression data.
        """
        # Iterate over the components
        for component in range(1, self.n_components + 1):
            _res = pls_regression(gene_exp, imaging_data.reshape(
                imaging_data.shape[0], 1),
                                  n_components=component,
                                  n_perm=0, n_boot=0)
            _exp_var = _res.get("varexp")
            _temp = _exp_var.cumsum(axis=0)
            _R = _temp[component - 1]
            _R_sq = np.zeros(1000)
            # Iterate over the permutations
            for i in tqdm(range(1000),
                          desc=f"Bootstrapping on PLS component {component}",
                          unit=" iterations"):
                y_data = permuted_imaging[:, i].reshape(
                    imaging_data.shape[0], 1)
                _res = pls_regression(gene_exp, y_data, n_components=component,
                                      n_perm=0, n_boot=0)
                _exp_var = _res.get("varexp")
                _R_sq[i] = _exp_var.cumsum(axis=0)[component - 1]
            # Set the results
            self.r2[component - 1] = _R
            self.p_val[component - 1] = np.sum(_R_sq >= _R) / 1000
        return


