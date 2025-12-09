import numpy as np
import pymc as pm
from pymc import logp


class TemperedBinomialMixture():
  def __init__(self, 
               X,
               n_trials,
               beta=None,
               n_components=2,
               weights_prior_params=None):
    self.X = X
    self.n_trials=n_trials
    
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

  def __enter__(self):
    self.model = pm.Model()
    self.model.__enter__()  # This puts it on PyMC's context stack
    # raw_weights = pm.Normal("raw_weights", mu=0, sigma=4, shape=self.n_components)
    
    # # Softmax to get simplex
    # weights = pm.Deterministic("weights", pt.special.softmax(raw_weights))
    
    weights = pm.Dirichlet("weights", a=self.weights_prior_params)

    probs = [pm.Uniform(f"p{i}", 0, 1) for i in range(self.n_components)]
    components = [pm.Binomial.dist(n=self.n_trials, p=p) for p in probs]
    
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

  def wbic(self, **kwargs):
    samples = self.sample(**kwargs)
    

