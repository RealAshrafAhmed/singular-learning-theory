import numpy as np
import pymc as pm
from pymc import logp
import pytensor as pt


class TemperedPoissonMixture():
  def __init__(self, 
               n_components=2, 
               weights_prior_params=np.ones(2),
               mus_prior_alpha=0.001,
               mus_prior_beta=0.001
              ):
    self.n_components = n_components
    if weights_prior_params:
      assert n_components == len(weights_prior_params), "n_components must match weights prior params length!"
      self.weights_prior_params=weights_prior_params
    else:
      self.weights_prior_params=np.ones(n_components)
      
    self.mus_prior_alpha=mus_prior_alpha
    self.mus_prior_beta=mus_prior_beta

  def __enter__(self, beta=None, X):
    if not beta:
        likelihood_power = 1/np.log(len(X))
    else:
        likelihood_power = beta

    model = pm.Model()

    weights = pm.Dirichlet("weights", a=self.weights_prior_params)  # 2 mixture weights    
    mus = [pm.Gamma(f"m{i}", alpha=self.mus_prior_alpha, beta=self.mus_prior_beta) for i in range(n_components)]
    components = [pm.Poisson.dist(mu=mu) for mu in mus]
    like = pm.Mixture("like", w=weights, comp_dists=components)

    tempered_log_likelihood = likelihood_power * logp(like, X)
    pm.Potential("tempered_log_likelihood", tempered_log_likelihood)
    self.model = model
    return self.model
    
  def __exit__(self, type, value, traceback):
      self.model = None


  def sample(**kwargs):
    return self.model.sample(**kwargs | {"nuts_sampler": "numpyro",
                                         "chains"=4,
                                         "max_treedepth": 50,
                                         "tune": 8000,
                                         "target_accept": 0.999})