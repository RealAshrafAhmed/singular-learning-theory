import numpy as np
from scipy import stats
import pandas as pd
from sklearn_extensions.confidence_regions import CRProvider


class CovarianceProvider:
    def __init__(self, fisher_mat, reg_level=1e-12):
        self.fisher_mat = fisher_mat
        self.reg_level = reg_level

    def reg_fisher_mat(self):
        # Add minimal regularization for numerical stability
        reg = self.reg_level * np.trace(self.fisher_mat) / len(self.fisher_mat)
        result = self.fisher_mat + reg * np.eye(len(self.fisher_mat))
        # print(f"unregularized fisher is")
        # print(self.fisher_mat)
        # print(f"regularized fisher is")
        # print(result)
        return result

    def cov_mat(self, params_index=None):
        # Invert Fisher matrix to get full covariance
        reg_fmatrix = self.reg_fisher_mat()

        try:
            full_cov = np.linalg.inv(reg_fmatrix)    
        except np.linalg.LinAlgError as e:
            raise ValueError(f"Fisher matrix is singular - cannot invert, cause{e}")

        # print(f"full cov is {full_cov}")
        if not params_index:
            params_index = range(0, len(reg_fmatrix))

        # Extract marginal covariance submatrix
        marginal_cov = full_cov[np.ix_(params_index, params_index)]
        # print(f"marginal cov is {marginal_cov}")
        return marginal_cov

    def standard_errors(self, params_index=None):
        # Extract standard errors (square root of diagonal)
        sigmas = np.sqrt(np.diag(self.cov_mat(params_index)))
        return sigmas


class GaussianCRProvider(CRProvider):

    def __init__(self, mle, cov_provider: CovarianceProvider):
        self.mle = mle
        self.cov_provider=cov_provider

    @property
    def type(self) -> str:
        return "gaussian"

    def intervals(self, params_index=None, confidence_level=0.95):
        # Extract standard errors (square root of diagonal)
        sigmas = self.cov_provider.standard_errors(params_index)

        # Get the z-score for the desired confidence level
        # For a two-tailed test, we use (1 + confidence_level)/2
        alpha = self.alpha_from_level(confidence_level)
        z_score = stats.norm.ppf(1-alpha)
        
        # Compute confidence intervals
        n_params = len(self.mle)
        intervals = np.zeros((n_params+1, 2)) # extra param for the constrained weight

        for i in range(len(sigmas)):
            intervals[i, 0] = self.mle[i] - z_score * sigmas[i]  # lower bound
            intervals[i, 1] = self.mle[i] + z_score * sigmas[i]  # upper bound
        
        return intervals

                         
        