import numpy as np
from scipy import stats
from scipy.stats import chi2

import pandas as pd
from sklearn_extensions.confidence_regions import CRProvider


def profile_ci_regions(grid, profile_ll, threshold):
    # Get indices where above threshold
    idx = np.where(profile_ll > threshold)[0]

    # Group consecutive indices
    intervals = np.split(grid[idx], np.where(np.diff(idx) != 1)[0] + 1)
    ci_intervals = [(region[0], region[-1]) for region in intervals if len(region) > 0]
    return ci_intervals


class profile_ll_factory:
    def __init__(self, model, data):
        self.model = model
        self.data = data

    def new_profile_ll(self):
        def result(params_index, params_value):
            pass

        return result


class ProfileLikelihoodCRProvider(CRProvider):
    def __init__(self, mle, profile_ll, param_linespaces):
        """
        num: the number of points controlling the discretization of a given parameter
        """
        self.mle = mle
        self.profile_ll = profile_ll
        self.num = num
        self.param_linespaces = param_linespaces

    @property
    def type(self) -> str:
        return "profile_likelihood"

    def intervals(self, params_index=None, confidence_level=0.95):
        # Extract standard errors (square root of diagonal)
        sigmas = self.cov_provider.standard_errors(params_index)

        # Get the z-score for the desired confidence level
        # For a two-tailed test, we use (1 + confidence_level)/2
        alpha = self.alpha_from_level(confidence_level)
        z_score = stats.norm.ppf(1 - alpha)

        # Compute confidence intervals
        n_params = len(self.mle)
        intervals = np.zeros(
            (n_params + 1, 2)
        )  # extra param for the constrained weight

        for pindex in range(len(sigmas)):
            param_linespace = self.param_linespaces[pindex]
            for pvalue in param_linespace:
                profile_ll_values = self.profile_ll(
                    params_index=[pindex], params_value=[self.mle[pindex]]
                )
                threshold = np.max(profile_ll) - chi2.ppf(confidence_level, 1) / 2
                pintervals = profile_ci_regions(
                    grid=param_linespace,
                    profile_ll=profile_ll_values,
                    threshold=threshold,
                )

            if len(pintervals) > 1:
                raise Exception(
                    f"Profile likelihood intervals are {len(pintervals)} not 1"
                )
            intervals[pindex, 0] = pintervals[0][0]  # lower bound
            intervals[pindex, 1] = pintervals[0][1]  # upper bound

        return intervals
