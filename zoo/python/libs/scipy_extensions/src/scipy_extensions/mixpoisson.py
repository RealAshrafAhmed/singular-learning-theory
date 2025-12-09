from scipy.stats import poisson
from scipy.special import logsumexp
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