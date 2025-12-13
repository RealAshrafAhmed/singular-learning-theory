from scipy.stats import poisson
from scipy.special import logsumexp
from pymc_extensions import pmx
import numpy as np

def logpmf(weights, mus, x):
  """
  Compute a probability mass function of binomial mixture
  x: samples
  n_trials: number of binomial trials
  params_weights: the mixing weights of each component, must add to 1
  params_probs: the probability of each component in the mixture
  log: return log probability, default is True
  flat: returna the probilities not just each component
  """
  negative_weights = [w >= 0 and w <= 1 for w in weights]
  assert np.all(negative_weights), (f"Weights {weights} must be between 0 and 1 inclusive")

  x = np.atleast_1d(x)
  result = poisson.logpmf(k=x[:, np.newaxis], mu=mus)
  return logsumexp(result + np.log(weights), axis=1)


def log_likelihood(weights, mus, x):
  """
  Compute the log likelihood mass for a binomial mixture
  using a batch of weights and probs

  This function calculates the log-likelihood for ALL observations (x) 
  across ALL psterior draws (weights, probs) at once.

  Args:
      weights (np.ndarray): Array of mixing weights, shape (N, 2).
      probs (np.ndarray): Array of component probabilities, shape (N, 2).
      n_trials (int): Number of binomial trials (n).
      x (np.ndarray): Observed data samples, shape (N,).

  Returns:
      np.ndarray: Array of summed log-likelihoods, shape (N,).
      Each element is the log-likelihood of the entire observed dataset (x)
      under one set of posterior parameters (a single draw).
  """
  x = np.atleast_1d(x)  # Shape (N_obs,)

  # 1. Expand dimensions for broadcasting:
  # x: (N_obs, 1, 1) -> Observations axis (0), Draw axis (1), Component axis (2)
  # weights_draws: (1, N_draws, 2)
  # probs_draws: (1, N_draws, 2)

  # 2. Calculate the log PMF for each observation, draw, and component:
  # binom.logpmf(k=x, n=n_trials, p=probs_draws)
  # k=x[..., np.newaxis, np.newaxis] broadcasts x against the posterior draws
  # p=probs_draws[np.newaxis, ...] broadcasts draws against observations
  
  # result_per_component will have shape (N_obs, N_draws, 2)
  # N_obs: axis 0 (Observations)
  # N_draws: axis 1 (Posterior Draws)
  # 2: axis 2 (Mixture Components)
  result_per_component = poisson.logpmf(
      k=x[:, np.newaxis, np.newaxis],  # Shape (N_obs, 1, 1)
      mu=mus[np.newaxis, ...]   # Shape (1, N_draws, 2)
  )
  
  # 3. Add log-weights to the result:
  # np.log(weights_draws) has shape (N_draws, 2)
  # When added to result_per_component (N_obs, N_draws, 2), broadcasting works:
  # (N_obs, N_draws, 2) + (1, N_draws, 2) -> (N_obs, N_draws, 2)
  log_weighted_result = result_per_component + np.log(weights[np.newaxis, ...])

  # 4. Perform logsumexp over the components (axis 2):
  # This sums up the component probabilities for each observation *and* draw.
  # log_likelihood_per_obs_and_draw shape: (N_obs, N_draws)
  log_likelihood_per_obs_and_draw = logsumexp(log_weighted_result, axis=2)

  # 5. Sum the log-likelihoods over all observations (axis 0):
  # This gives the total log-likelihood for the *entire dataset* for each draw.
  # total_log_likelihood_per_draw shape: (N_draws,)
  total_log_likelihood_per_draw = np.sum(log_likelihood_per_obs_and_draw, axis=0)
  
  return total_log_likelihood_per_draw


def n_components(weights, mus):
  """
  Computes the number of components from a given list of parameter values
  """
  # if any of the weights is 1, we just have a single component
  if len([w for w in weights if w == 1]) > 0:
    return 1

  indices_of_positive_weights = [i for i, w in enumerate(weights) if w > 0]
  equal_mus = 0
  for i in range(len(indices_of_positive_weights)):
    for j in range(i, len(indices_of_positive_weights)):
      if i == j: # doesn't count, we need different indices
        continue

      if mus[i] == mus[j]:
        equal_mus += 1

  return len(indices_of_positive_weights)-equal_mus


def rlct(d:int, r: int, k: int):
  assert d == 1, "only supports when the dimension of data is 1"
  return (3*r + k-2)/4


def afe(X, mus, weights, rlct):
  n=len(X)
  data_dim = X.ndim # the dimension of data
  
  average_log_likelihood = None
  log_p = logpmf(weights=weights, mus=mus, x=X) # sample loglikelihood under the parameter
  average_log_likelihood = log_p.mean()

  return -n*average_log_likelihood+rlct*np.log(n)
# def create_mixbinom_profile_kl(
#     n_components, n_trials, truth, x_param, y_param, profile_grid_size=20
# ):
#     """
#     Compute profile KL divergence

#     Automatically infers which parameters to profile over based on:
#     - Parameters that exist in the model (probs and weights)
#     - Parameters NOT being varied (x_param, y_param)

#     Parameters:
#     -----------
#     x_param, y_param: dict
#         Parameters to vary on grid, format: {"name": "weights"/"probs", "index": int}
#     profile_grid_size: int
#         Number of grid points for each profiled parameter
#     """
#     k = np.arange(0, n_trials + 1, dtype=int)

#     truth_probs = truth[:n_components]
#     truth_weights = truth[n_components:]

#     log_q = mixbinom_logpmf(n_trials, weights=truth_weights, probs=truth_probs)(k)

#     # Infer which parameters to profile over
#     profile_params = _infer_profile_params(n_components, x_param, y_param)

#     # Create grid for each profile parameter
#     profile_grids = {}
#     for param_key, param_info in profile_params.items():
#         profile_grids[param_key] = np.linspace(0.01, 0.99, profile_grid_size)

#     def fun(x, y):
#         if x < 0 or x > 1 or y < 0 or y > 1:
#             return np.nan

#         min_kl = np.inf

#         # Iterate over all combinations of profile parameters
#         for profile_values in _iterate_profile_combinations(
#             profile_grids, profile_params
#         ):
#             m_probs = np.zeros(n_components)
#             m_weights = np.zeros(n_components)

#             # Set x and y parameters
#             m_probs, m_weights = __set_param(
#                 m_probs, m_weights, x_param["name"], x_param["index"], x
#             )
#             m_probs, m_weights = __set_param(
#                 m_probs, m_weights, y_param["name"], y_param["index"], y
#             )

#             # Set all profile parameters
#             for param_key, value in profile_values.items():
#                 param_info = profile_params[param_key]
#                 m_probs, m_weights = __set_param(
#                     m_probs, m_weights, param_info["name"], param_info["index"], value
#                 )

#             # Handle weight constraint (last weight is determined)
#             m_weights[n_components - 1] = 1 - np.sum(m_weights[:-1])

#             # Validation
#             if not np.isclose(np.sum(m_weights), 1.0, atol=1e-6):
#                 continue
#             if np.any(m_weights < 0) or np.any(m_weights > 1):
#                 continue
#             if np.any(m_probs < 0) or np.any(m_probs > 1):
#                 continue

#             # Compute KL
#             log_p = mixbinom_logpmf(n_trials, weights=m_weights, probs=m_probs)(k)

#             mask = (
#                 np.isfinite(log_q)
#                 & np.isfinite(log_p)
#                 & (log_q > -100)
#                 & (log_p > -100)
#             )

#             if np.any(mask):
#                 kl = np.sum(np.exp(log_q[mask]) * (log_q[mask] - log_p[mask]))
#                 min_kl = min(min_kl, kl)

#                 # Early stopping if we found a perfect match
#                 if min_kl < 1e-8:
#                     break

#         # return max(0, min_kl) if np.isfinite(min_kl) else np.nan
#         return max(0, min_kl)

#     return np.vectorize(fun, otypes=[np.float64]) 