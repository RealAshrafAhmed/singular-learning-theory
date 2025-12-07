import numpy as np
import pymc as pm
from pymc import logp
import pytensor.tensor as pt


class TemperedPoissonMixture():
  def __init__(self, 
               X,
               beta=None,
               n_components=2,
               weights_prior_params=None,
               mus_log_prior_mean=0,
               mus_log_prior_sigma=2):
    self.X = X
    if not beta:
      self.likelihood_power = 1/np.log(len(X))
    else:
      self.likelihood_power = beta
      
    self.n_components = n_components
    if weights_prior_params:
      assert n_components == len(weights_prior_params), "n_components must match weights prior params length!"
      self.weights_prior_params=weights_prior_params
    else:
      self.weights_prior_params=np.ones(n_components)

    self.mus_log_prior_mean=mus_log_prior_mean
    self.mus_log_prior_sigma=mus_log_prior_sigma

  def __enter__(self):
    self.model = pm.Model()
    self.model.__enter__()  # This puts it on PyMC's context stack
    # raw_weights = pm.Normal("raw_weights", mu=0, sigma=4, shape=self.n_components)
    
    # # Softmax to get simplex
    # weights = pm.Deterministic("weights", pt.special.softmax(raw_weights))
    
    weights = pm.Dirichlet("weights", a=self.weights_prior_params)  # 2 mixture weights    
    # mus = [pm.Gamma(f"m{i}", alpha=self.mus_prior_alpha, beta=self.mus_prior_beta) for i in range(self.n_components)]

    # Rates via exp (log-normal parameterization)
    log_mus = pm.Normal("log_mus", 
                        mu=self.mus_log_prior_mean, 
                        sigma=self.mus_log_prior_sigma, 
                        shape=self.n_components)
    
    mus = pm.Deterministic("mus", pt.exp(log_mus))
    
    components = [pm.Poisson.dist(mu=mu) for mu in mus]
    like_dist = pm.Mixture.dist(w=weights, comp_dists=components)

    # Compute tempered log-likelihood and sum over observations
    tempered_log_likelihood = self.likelihood_power * pm.logp(like_dist, self.X).sum()
    pm.Potential("tempered_log_likelihood", tempered_log_likelihood)
    return self
    
  def __exit__(self, exc_type, exc_value, traceback):
      # Exit the PyMC model context
      self.model.__exit__(exc_type, exc_value, traceback)
      return False  # Don't suppress exceptions

  def sample(self, **kwargs):
    defaults = {
        "nuts_sampler": "numpyro",
        "chains": 4,
        "max_treedepth": 50,
        "tune": 8000,
        "target_accept": 0.999
    }
    merged = {**defaults, **kwargs}
    return pm.sample(**merged)