from dataclasses import dataclass
from typing import Optional
import numpy as np
import pymc as pm
import arviz as az
from scipy.special import logsumexp

import models 
import utils

@dataclass
class BayesianModel:
    """
    This class bundles together a PyMC model, the observations used to fit it,
    and the resulting inference data from sampling.
    
    Attributes:
    -----------
    model : pymc.Model
        The PyMC model object containing priors, likelihood, and model structure
    observations : numpy.ndarray
        The observed data used to fit the model
    inference_data : arviz.InferenceData, optional
        The ArviZ InferenceData object containing posterior samples and diagnostics.
        Initially None until sampling is performed.
        
    Methods:
    --------
    sample(draws=1000, tune=1000, **kwargs)
        Sample from the posterior using PyMC's default sampler
    WBIC()
        Compute and return the WBIC
    """
    model: pm.Model
    observations: np.ndarray
    inference_data: Optional[az.InferenceData] = None
    
    def sample(self, draws: int = 10000, tune: int = 10000, **kwargs) -> az.InferenceData:
        """
        Sample from the posterior distribution.
        
        Parameters:
        -----------
        draws : int, default=10000
            Number of posterior samples to draw
        tune : int, default=10000
            Number of tuning/warmup samples
        **kwargs
            Additional arguments passed to pm.sample()
            
        Returns:
        --------
        arviz.InferenceData
            The inference data object containing samples and diagnostics
        """
        with self.model:
            self.inference_data = pm.sample(draws=draws, tune=tune, **kwargs)
        return self.inference_data

    def WBIC(self):
        if self.inference_data is None:
           print("Warning: Posterior has not yet been sampled; sampling with default parameters now.")
           self.sample(
            draws=2500,
            tune=2500,
            chains=4,
            max_treedepth=50,
            target_accept=.995
           )
        # log_likelihood has shape (n_chain, n_draw, n_data)
        # so summing over axis=2 is computing log p(w|X) = sum(log p(w|xi))
        ll = np.sum(self.inference_data.posterior['log_likelihood'].values, axis=2)
        n_chain, n_draw = ll.shape[0], ll.shape[1]
        # sometimes ll has an axis of size 1 and sometimes it doesn't, not sure why.
        ll = ll.reshape(n_chain, n_draw) #get rid of axis of size 1 
        random_mask = np.random.randint(0, n_chain, size=n_draw)
        ll_samples = [ll[random_mask[i], i] for i in range(n_draw)]
        WBIC = np.mean(ll_samples)

        return WBIC
        
    def empirical_loss(self):
        if self.inference_data is None:
            print("Sample the posterior first!")
            return 1
        ll = self.inference_data.posterior['log_likelihood'].values
        n_chain, n_draw, n_data = ll.shape
        chain_indices = np.random.randint(0, n_chain, size=n_draw)
        draw_indices = np.arange(n_draw)
        selected_ll = ll[chain_indices, draw_indices, :]

        ## This log_mean_exp is computing the average risk per observation
        ## then we average that over all samples.
        log_mean_exp = logsumexp(selected_ll, axis=0) - np.log(n_draw)
        empirical_loss = -np.mean(log_mean_exp)

        return empirical_loss
 
    def learning_coefficient(self, method='imai'):
        n_data = self.observations.size
        ## Takio & Suzuki
        if method == 'empirical_loss':
            WBIC = self.WBIC()
            Tn = self.empirical_loss()
            LC = (-WBIC - n_data*Tn)/np.log(n_data)
            return LC
        if method == 'imai':
            ll = np.sum(self.inference_data.posterior['log_likelihood'].values, axis=2)
            n_chain, n_draw = ll.shape
            random_mask = np.random.randint(0, n_chain, size=n_draw)
            ll_samples = [ll[random_mask[i], i] for i in range(n_draw)]
            LC = np.var(ll_samples)/(np.log(n_data)**2)
            return LC 
            

if __name__ == "__main__":
    np.random.seed(42)
    X = np.random.normal(loc=0, scale=1, size=500)
    model = models.tempered_gaussian_mixture(
        X,
        n_components=2,
        beta = 1/np.log(X.size),
    )
    gmm = BayesianModel(
        model=model,
        observations=X,
    )
    gmm.sample(
        draws=2000,
        tune=1000,
        chains=4,
        max_treedepth=50,
        target_accept=.995
    )
    gmm.WBIC()
    LC = gmm.learning_coefficient()
    print(f'estimated WBIC: {gmm.WBIC()}')
    print(f'estimated LC: {LC}')