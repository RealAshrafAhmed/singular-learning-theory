import numpy as np
from scipy import stats
import pandas as pd
from tqdm.notebook import tqdm
from sklearn.utils import resample
from sklearn_extensions.confidence_regions import CRProvider


class BootstrapCRProvider(CRProvider):
    def __init__(self, n_params, estimates):
        self.estimates = estimates
        self.n_params = n_params

    @property
    def type(self) -> str:
        return "bootstrap"

    def intervals(self, params_index=None, confidence_level=0.95):
        alpha = self.alpha_from_level(confidence_level)
        print(f"alpha_lower={alpha}, alpha_upper={1-alpha}")

        # Compute confidence intervals
        intervals = np.zeros((self.n_params+1, 2)) # extra param for the constrained weight
        
        for i in range(self.n_params):
            intervals[i, 0] = np.quantile(self.estimates[:,i], alpha)
            intervals[i, 1] = np.quantile(self.estimates[:,i], 1-alpha)
        
        return intervals


class BootstrapCRProviderFactory:

    def __init__(self, data, point_estimator, n_params, n_iterations=10**4, random_state=None):
        self.data=data
        self.point_estimator=point_estimator
        self.n_params = n_params
        self.samples_per_iteration=len(data)
        self.n_iterations=n_iterations
        self.random_state=random_state

    def new_provider(self) -> BootstrapCRProvider:
        estimates = np.zeros((self.n_iterations, self.n_params))
        for i in tqdm(range(self.n_iterations), desc=f"bootstrap "):
            X_boot = resample(self.data, random_state=self.random_state, n_samples=self.samples_per_iteration)
            estimate, estimate_detail = self.point_estimator(X_boot)
            estimates[i, :]= estimate

        return BootstrapCRProvider(n_params=self.n_params, estimates=estimates)

    


