import numpy as np
import pymc as pm
from pymc import logp
from scipy_extensions import mixbinom
from pymc_extensions import pmx
from joblib import Parallel, delayed
import warnings


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

  def wbic(self, trace):
    prob_names = [f"p{i}" for i in range(self.n_components)]
    weights = pmx.column_stack_vars(trace, ["weights"])
    probs = pmx.column_stack_vars(trace, prob_names)
    log_likelihood = mixbinom.log_likelihood(weights, probs, n=self.n_trials, x=self.X)
    return -log_likelihood.mean()


def __approx_tempered_nll(X, n_trials, beta, nuts_sampler):
  with warnings.catch_warnings():
    warnings.filterwarnings("ignore", message="Some donated buffers were not usable")

    with TemperedBinomialMixture(X=X, n_trials=n_trials, beta=beta) as model:
      idata = model.sample(draws=1000,
                           chains=1, 
                           max_tree_depth=50,
                           tune=2000,
                           target_accept=0.99,
                           cores=1,
                           nuts_sampler=nuts_sampler,
                           progressbar=False,
                           compute_convergence_checks=False
                          )
        
      return model.wbic(idata)


# def free_energy(n_trials, X, betas, nuts_sampler="nutpie"):
#   """Compute free energy using thermodynamic integral"""
#   if len(betas)==0:
#     betas = np.linspace(0, 1, 30) **2

#   wbic_betas=[]
#   for beta in betas:
#     wbic_betas.append(__approx_tempered_nll(X, n_trials, beta, nuts_sampler))

#   # Compute integration numerically using the trapezoidal rule
#   return np.trapz(wbic_betas, betas)

def free_energy(n_trials, X, betas, parallel_n_jobs=1, parallel_verbose=10, nuts_sampler="nutpie"):
  """Compute free energy using thermodynamic integral"""
  if len(betas)==0:
    betas = np.linspace(0, 1, 30) **2

  results = Parallel(n_jobs=parallel_n_jobs, verbose=parallel_verbose)(
    delayed(__approx_tempered_nll)(
      X, n_trials, beta, nuts_sampler
    )
    for beta in betas
  )

  # Compute integration numerically using the trapezoidal rule
  return np.trapz(results, betas)

    

