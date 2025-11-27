from typing import Callable, List
from enum import Enum
import numpy as np
from numpy.typing import NDArray
from itertools import permutations
from scipy.stats import binom, chi2
from scipy.special import logsumexp
from sklearn.mixture._base import BaseMixture
from sklearn.utils.validation import check_X_y, check_array
from sklearn.utils import check_random_state
import warnings
from joblib import Parallel, delayed
import multiprocessing as mp
from sklearn_extensions.random import dirichlet
from sklearn_extensions.mixture._fun import __set_params, __set_param
from sklearn.utils._array_api import (
    get_namespace,
    get_namespace_and_device
)

from scipy.optimize import minimize_scalar


def create_mixbinom_profile_kl(n_components, n_trials, truth, 
                                       x_param, y_param,
                                       profile_grid_size=20):
    """
    Compute profile KL divergence
    
    Automatically infers which parameters to profile over based on:
    - Parameters that exist in the model (probs and weights)
    - Parameters NOT being varied (x_param, y_param)
    
    Parameters:
    -----------
    x_param, y_param: dict
        Parameters to vary on grid, format: {"name": "weights"/"probs", "index": int}
    profile_grid_size: int
        Number of grid points for each profiled parameter
    """
    k = np.arange(0, n_trials + 1, dtype=int)
    
    truth_probs = truth[:n_components]
    truth_weights = truth[n_components:]
    
    log_q = mixbinom_logpmf(n_trials, weights=truth_weights, probs=truth_probs)(k)
    
    # Infer which parameters to profile over
    profile_params = _infer_profile_params(n_components, x_param, y_param)
    
    # Create grid for each profile parameter
    profile_grids = {}
    for param_key, param_info in profile_params.items():
        profile_grids[param_key] = np.linspace(0.01, 0.99, profile_grid_size)
    
    def fun(x, y):
        if x < 0 or x > 1 or y < 0 or y > 1:
            return np.nan
        
        min_kl = np.inf
        
        # Iterate over all combinations of profile parameters
        for profile_values in _iterate_profile_combinations(profile_grids, profile_params):
            m_probs = np.zeros(n_components)
            m_weights = np.zeros(n_components)
            
            # Set x and y parameters
            m_probs, m_weights = __set_param(m_probs, m_weights,
                                             x_param["name"], x_param["index"], x)
            m_probs, m_weights = __set_param(m_probs, m_weights,
                                             y_param["name"], y_param["index"], y)
            
            # Set all profile parameters
            for param_key, value in profile_values.items():
                param_info = profile_params[param_key]
                m_probs, m_weights = __set_param(m_probs, m_weights,
                                                 param_info["name"],
                                                 param_info["index"],
                                                 value)
            
            # Handle weight constraint (last weight is determined)
            m_weights[n_components-1] = 1 - np.sum(m_weights[:-1])
            
            # Validation
            if not np.isclose(np.sum(m_weights), 1.0, atol=1e-6):
                continue
            if np.any(m_weights < 0) or np.any(m_weights > 1):
                continue
            if np.any(m_probs < 0) or np.any(m_probs > 1):
                continue
            
            # Compute KL
            log_p = mixbinom_logpmf(n_trials, weights=m_weights, probs=m_probs)(k)
            
            mask = (np.isfinite(log_q) & np.isfinite(log_p) & 
                   (log_q > -100) & (log_p > -100))
            
            if np.any(mask):
                kl = np.sum(np.exp(log_q[mask]) * (log_q[mask] - log_p[mask]))
                min_kl = min(min_kl, kl)
                
                # Early stopping if we found a perfect match
                if min_kl < 1e-8:
                    break
        
        # return max(0, min_kl) if np.isfinite(min_kl) else np.nan
        return max(0, min_kl)
    
    return np.vectorize(fun, otypes=[np.float64])


def _check_weights(weights, atol=1e-6):
    if not np.isclose(np.sum(weights), 1.0, atol=atol):
        return False
    if np.any(weights < 0) or np.any(weights > 1):
        return False
    if np.any(weights < 0) or np.any(weights > 1):
        return False

    return True


def mixbinom_kl_divergence(n_trials: int) -> Callable[[List[NDArray[np.float64]], NDArray[np.float64]], NDArray[np.float64]]:
    k = np.arange(0, n_trials + 1, dtype=int)  # precompute
    
    def _parse_params(params: NDArray) -> tuple[NDArray, NDArray]:
        """Parse parameter array into probs and weights."""
        n_components = (params.shape[-1] + 1) // 2
        probs = params[..., :n_components]
        weights_partial = params[..., n_components:]
        weights_last = 1 - weights_partial.sum(axis=-1, keepdims=True)
        weights = np.concatenate([weights_partial, weights_last], axis=-1)
        return probs, weights
    
    def _compute_log_pmf(probs: NDArray, weights: NDArray) -> NDArray:
        """
        Compute log PMF for mixture binomial.
        probs: (..., n_components)
        weights: (..., n_components)
        returns: (..., n_trials + 1)
        """
        # binom.logpmf with broadcasting
        # probs[..., np.newaxis]: (..., n_components, 1)
        # k: (n_trials + 1,)
        # result: (..., n_components, n_trials + 1)
        log_pmf_components = binom.logpmf(k, n_trials, probs[..., np.newaxis])
        log_weights = np.log(weights)[..., np.newaxis]  # (..., n_components, 1)
        
        # logsumexp over components
        log_pmf = logsumexp(log_weights + log_pmf_components, axis=-2)
        return log_pmf
    
    def fun(references: List[NDArray[np.float64]], queries: NDArray[np.float64]) -> NDArray[np.float64]:
        """
        Compute KL divergence from each query to each reference.
        
        Args:
            references: list of 1D parameter vectors [p0, p1, ..., w0, w1, ...]
            queries: (n_queries, dims) array of parameter vectors
        
        Returns:
            (n_queries, n_refs) array of KL divergences
        """
        refs = np.stack(references)  # (n_refs, dims)
        
        ref_probs, ref_weights = _parse_params(refs)
        query_probs, query_weights = _parse_params(queries)
        
        # Compute log PMFs
        log_p = _compute_log_pmf(ref_probs, ref_weights)  # (n_refs, n_trials + 1)
        log_q = _compute_log_pmf(query_probs, query_weights)  # (n_queries, n_trials + 1)
        
        # KL(q || p) = sum_k q(k) * (log q(k) - log p(k))
        # Broadcasting for all pairs
        q = np.exp(log_q)  # (n_queries, n_trials + 1)
        
        log_q_expanded = log_q[:, np.newaxis, :]  # (n_queries, 1, n_trials + 1)
        log_p_expanded = log_p[np.newaxis, :, :]  # (1, n_refs, n_trials + 1)
        q_expanded = q[:, np.newaxis, :]          # (n_queries, 1, n_trials + 1)
        
        # Handle numerical issues
        diff = log_q_expanded - log_p_expanded
        mask = np.isfinite(diff) & (log_q_expanded > -100)
        
        kl_terms = np.where(mask, q_expanded * diff, 0.0)
        kl = np.sum(kl_terms, axis=-1)  # (n_queries, n_refs)
        
        return kl
    
    return fun


def _infer_profile_params(n_components, x_param, y_param):
    """
    Infer which parameters to profile over
    
    Returns:
    --------
    dict: {param_key: {"name": str, "index": int}}
        param_key format: "probs_0", "weights_1", etc.
    """
    profile_params = {}
    
    # Track which parameters are being varied
    varied_params = {
        f"{x_param['name']}_{x_param['index']}",
        f"{y_param['name']}_{y_param['index']}"
    }
    
    # All probs except the ones being varied
    for i in range(n_components):
        param_key = f"probs_{i}"
        if param_key not in varied_params:
            profile_params[param_key] = {"name": "probs", "index": i}
    
    # All weights except the ones being varied and the last one (constrained)
    for i in range(n_components - 1):  # Last weight is 1 - sum(others)
        param_key = f"weights_{i}"
        if param_key not in varied_params:
            profile_params[param_key] = {"name": "weights", "index": i}
    
    return profile_params


def _iterate_profile_combinations(profile_grids, profile_params):
    """
    Generate all combinations of profile parameter values
    
    Yields:
    -------
    dict: {param_key: value}
    """
    from itertools import product
    
    if not profile_params:
        # No parameters to profile - just yield empty dict
        yield {}
        return
    
    # Get parameter keys in consistent order
    param_keys = sorted(profile_params.keys())
    
    # Get grid values for each parameter
    grid_values = [profile_grids[key] for key in param_keys]
    
    # Generate all combinations
    for values in product(*grid_values):
        yield {key: val for key, val in zip(param_keys, values)}


class SetType(Enum):
    Point = 1
    AlgebraicVariety = 2


class Divergence(Enum):
    Euclidean = "Euclidean"
    KL = "KL"


def mixbinom_divergence(n_trials: int, 
                        kind: Divergence=None,
                        allow_marginals=False):
    match kind:
        case Divergence.KL:
            if allow_marginals:
                raise Exception("marginalization not allowed for KL")
            return mixbinom_kl_divergence(n_trials)
        
        case _: # default to euclidean distance
            def euclidean_div(references: List[np.float64], queries: NDArray[np.float64]) -> NDArray[np.float64]:
                refs = np.stack(references)  # list of 1D → (n_refs, dims)
                if allow_marginals:
                    to_shape = queries.shape[1]
                    refs = refs[:, :to_shape]

                diff = queries[:, np.newaxis, :] - refs[np.newaxis, :, :]
                return np.linalg.norm(diff, axis=-1)

            return euclidean_div


def mixbinom_logpmf(n_trials, weights, probs, flat=True):
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
    assert np.all(negative_weights), f"Weights {weights} must be between 0 and 1 inclusive"
    n_trials=n_trials

    def logpmf(x):
        x = np.atleast_1d(x)
        result = binom.logpmf(k=x[:, np.newaxis], n=n_trials, p=probs)
        return logsumexp(result + np.log(weights), axis=1)

    return logpmf


def bad_mixbinom_logpmf(n_trials, weights, probs, flat=True, numerical_stability=None):
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
    assert np.all(negative_weights), f"Weights {weights} must be between 0 and 1 inclusive"
    n_components = len(weights)
    n_trials=n_trials

    def pmf(x):
        x = np.atleast_1d(x)
        n_samples=len(x)
        result = np.zeros((n_samples, n_components)) # +1 to capture the mixture
        for k in range(len(weights)):
            result[:, k] = binom.logpmf(x, n=n_trials, p=probs[k])
            # if log:
            #     # xx = binom.logpmf(k=x, n=n_trials, p=params_probs[k])
            #     # print(xx.shape)
            # else:
            #     result[:, k] = binom.pmf(k=x, n=n_trials, p=params_probs[k])
            # total = np.zeros(n_samples)
            # for i in range(n_components):
                    # print(result[:, i], params_weights[i])
                    # print(result[:, i]*params_weights[i])
            # print(f"result={result}, rz(result)={rz(result)}, np.log(params_weights)={np.log(params_weights)}, rz(np.log(params_weights))={rz(np.log(params_weights))}")
        total = logsumexp(result + np.log(weights), axis=1)
        return np.sum(total)

    return np.vectorize(pmf)


def _estimate_binomial_parameters(n_components, X, resp, reg_prob, random_state, xp):
    if xp is None:
        import numpy as xp

    if resp is None:
    # Initialize mixing weights uniformly
        weights_ = dirichlet(k=n_components, n_samples=1, alpha=1.0, xp=xp)
    else:
        weights_ = resp

    random_state = check_random_state(random_state)
    # Init probs randomly
    probs_ = random_state.uniform(0.1, 0.9, n_components)
    
    # Ensure probabilities are in the regularized range
    probs_ = np.clip(probs_, reg_prob, 1 - reg_prob)
    return weights_, probs_

        
class BinomialMixture(BaseMixture):
    """
    Binomial Mixture Model using Expectation-Maximization algorithm.
    
    This class implements a mixture of binomial distributions that follows
    scikit-learn's estimator interface.
    
    Parameters
    ----------
    n_components : int, default=2
        The number of mixture components.

    n_trials : int or array-like, default=1
        Number of trials for binomial distribution. If int, same for all samples.
        If array-like, must have same length as X.

    max_iter : int, default=100
        Maximum number of EM iterations.
        
    tol : float, default=1e-6
        Convergence tolerance.
        
    random_state : int, RandomState instance or None, default=None
        Random seed for initialization.
        
    init_params : str, default='random'
        Method for initialization ('random' or 'kmeans').
        
    reg_covar : float, default=1e-6
        Regularization parameter to avoid numerical issues.
    
    Attributes
    ----------
    weights_ : array-like of shape (n_components,)
        Mixing weights for each component.
    probs_ : array-like of shape (n_components,)
        Success probabilities for each binomial component.
    converged_ : bool
        True if the EM algorithm converged.
    n_iter_ : int
        Number of iterations performed.
    lower_bound_ : float
        Log-likelihood of the best fit.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn_ext.mixture import BinomialMixture
    >>> X = np.array([[8, 12], [18, 7], [15, 16]])
    >>> gm = BinomialMixture(n_components=2, random_state=0).fit(X)
    >>> gm.probs_
    array([[.3,  .1],
           [.2,  .6]])
    >>> gm.predict([[0, 0], [12, 3]])
    array([1, 0])
    """
    def __init__(
        self,
        n_trials,
        probs_init=None,
        reg_prob=1e-9,
        n_components=1,
        *,
        tol=1e-8,
        reg_covar=1e-9,
        max_iter=10**5,
        n_init=30,
        weights_init=None,
        random_state=None,
        warm_start=False,
        verbose=0,
        verbose_interval=10,
        enforce_ordering=True,
    ):
        super().__init__(
            n_components=n_components,
            tol=tol,
            reg_covar=reg_covar,
            max_iter=max_iter,
            n_init=n_init,
            init_params="random",
            random_state=random_state,
            warm_start=warm_start,
            verbose=verbose,
            verbose_interval=verbose_interval,
        )
        self.reg_prob=reg_prob
        self.probs_init=probs_init
        self.weights_init=weights_init

        if (
            weights_init is not None
            and probs_init is not None
        ):
            self.probs_ = self.probs_init
            self.weights_ =self.weights_init

        self.n_trials=n_trials
        self.enforce_ordering = enforce_ordering
    
    def _check_parameters(self, X):
        """Check input parameters."""
        if self.n_components is not None and self.n_components < 1:
            raise ValueError("n_components must be >= 1")

        if self.max_iter is not None and self.max_iter < 1:
            raise ValueError("max_iter must be >= 1")

        if self.tol is not None and self.tol <= 0:
            raise ValueError("tol must be > 0")
    
    def _initialize(self, X, resp, xp=None):
        """Initialization of the Binomial mixture parameters.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)

        resp : array-like of shape (n_samples, n_components)
        """

        self.weights_, self.probs_ = _estimate_binomial_parameters(n_components=self.n_components, 
                                                                   X=X, 
                                                                   resp=resp, 
                                                                   reg_prob=self.reg_prob, 
                                                                   random_state=self.random_state, 
                                                                   xp=xp)

    def _initialize_parameters(self, X, random_state, xp=None):
        # If all the initial parameters are all provided, then there is no need to run
        # the initialization.
        compute_resp = (
            self.weights_init is None
            or self.probs_init is None
        )
        if compute_resp:
            self.weights_, self.probs_ = _estimate_binomial_parameters(n_components=self.n_components, 
                                                                       X=X, 
                                                                       resp=None, 
                                                                       reg_prob=self.reg_prob, 
                                                                       random_state=random_state, 
                                                                       xp=xp)
        else:
            self._initialize(X, None, xp=xp)


    def _get_parameters(self):
        return (
            self.weights_,
            self.probs_
        )

    def _set_parameters(self, params, xp=None):
        xp, _, device_ = get_namespace_and_device(params, xp=xp)
        (
            self.weights_,
            self.probs_
        ) = params

    def _estimate_log_prob(self, X, xp=None):
        result = self._compute_log_likelihood(X=X)
        return result

    def _estimate_log_weights(self, xp=None):
        xp, _ = get_namespace(self.weights_, xp=xp)
        return xp.log(self.weights_)
            
    def _compute_log_likelihood(self, X):
        """Compute log-likelihood for each sample and component."""
        n_samples = X.shape[0]
        log_prob = np.zeros((n_samples, self.n_components))
        obs = X[:,0]
        # n_trials = X[:,1]
        
        for k in range(self.n_components):            
            # Compute log probability for binomial distribution
            log_prob[:, k] = binom.logpmf(obs, self.n_trials, self.probs_[k])
        
        return log_prob

    def _compute_lower_bound(self, _, log_prob_norm):
        return log_prob_norm

    def _enforce_ordering(self):
        if not self.enforce_ordering:
            return
        """Simple function to enforce p1 >= p2 >= ... >= pK ordering."""
        # Get sorted indices
        sort_indices = np.argsort(self.probs_)
        self.probs_ = self.probs_[sort_indices]
        
        if(not isinstance(self.weights_, np.ndarray)):
            self.weights_ = np.array(self.weights_)
            
        self.weights_ = self.weights_[sort_indices]
        
    def _m_step(self, X, log_resp):
        """Maximization step."""
        n_samples = X.shape[0]

        resp = np.exp(log_resp)
        # Update mixing weights
        resp_sum = np.sum(resp, axis=0) + 10 * np.finfo(resp.dtype).eps # ensure some numerical stability
        self.weights_ = resp_sum / n_samples

        # adjust the weights percision
        self.weights_ = [np.clip(w, self.reg_covar, 1 - self.reg_covar) for w in self.weights_]
        
        # Update success probabilities
        for k in range(self.n_components):
            # Weighted maximum likelihood estimation
            # Different number of trials for each sample
            weighted_successes = np.sum(resp[:, k] * X[:,0])
            weighted_trials = np.sum(resp[:, k] * self.n_trials)
            
            if weighted_trials > 0:
                self.probs_[k] = weighted_successes / weighted_trials
            else:
                self.probs_[k] = 0.5  # Default value
            
            # Ensure valid probability range
            self.probs_[k] = np.clip(self.probs_[k], self.reg_covar, 1 - self.reg_covar)

        self._enforce_ordering()

    
    # def predict(self, X):
    #     """
    #     Predict the labels for the data samples in X.
        
    #     Parameters
    #     ----------
    #     X : array-like of shape (n_samples,)
    #         Number of successes for each sample.
            
    #     Returns
    #     -------
    #     labels : array-like of shape (n_samples,)
    #         Component labels for each sample.
    #     """
    #     X = check_array(X, ensure_2d=False)
    #     if X.ndim == 2 and X.shape[1] == 1:
    #         X = X.ravel()
        
    #     resp, _ = self._e_step(X)
    #     return np.argmax(resp, axis=1)
    
    # def predict_proba(self, X):
    #     """
    #     Predict posterior probabilities of each component for each observation.
        
    #     Parameters
    #     ----------
    #     X : array-like of shape (n_samples,)
    #         Number of successes for each sample.
            
    #     Returns
    #     -------
    #     resp : array-like of shape (n_samples, n_components)
    #         Posterior probabilities of each mixture component for each observation.
    #     """
    #     X = check_array(X, ensure_2d=False)
    #     if X.ndim == 2 and X.shape[1] == 1:
    #         X = X.ravel()
        
    #     resp, _ = self._e_step(X)
    #     return resp
    
    # def score(self, X, y=None):
    #     """
    #     Compute the average log-likelihood of the samples.
        
    #     Parameters
    #     ----------
    #     X : array-like of shape (n_samples,)
    #         Number of successes for each sample.
    #     y : ignored
    #         Not used, present for API consistency.
            
    #     Returns
    #     -------
    #     log_likelihood : float
    #         Average log-likelihood of the samples under the model.
    #     """
    #     X = check_array(X, ensure_2d=False)
    #     if X.ndim == 2 and X.shape[1] == 1:
    #         X = X.ravel()
        
    #     _, log_likelihood = self._e_step(X)
    #     return log_likelihood / len(X)
    
    def sample(self, n_samples=1, random_state=None):
        """
        Generate random samples from the fitted mixture model.
        
        Parameters
        ----------
        n_samples : int, default=1
            Number of samples to generate.
        random_state : int, RandomState instance or None, default=None
            Random seed.
            
        Returns
        -------
        X : array-like of shape (n_samples,)
            Generated samples.
        y : array-like of shape (n_samples,)
            Component labels for each sample.
        """
        random_state = check_random_state(random_state)
        
        # Sample component assignments
        component_samples = random_state.choice(self.n_components, 
                                                size=n_samples, 
                                                p=self.weights_)
        
        # Generate samples from each component
        samples = np.zeros(n_samples, dtype=int)
        
        for k in range(self.n_components):
            mask = component_samples == k
            n_component_samples = np.sum(mask)
            if n_component_samples > 0:
                if isinstance(self.n_trials, (int, float)):
                    n_trials = self.n_trials
                else:
                    # For varying n_trials, use the mean
                    n_trials = int(np.mean(self.n_trials))
                
                samples[mask] = random_state.binomial(n_trials, 
                                                      self.probs_[k], 
                                                      size=n_component_samples
                )
        
        return samples, component_samples

    def point_estimate(self):
        self._enforce_ordering()
        params_detail = {}
        for i in range(self.n_components):
            params_detail[f"p{i}"] = {
                "pindex": i,
                "val": self.probs_[i]
            }

        for i in range(self.n_components-1):
            params_detail[f"w{i}"] = {
                "pindex": self.n_components+i,
                "val": self.weights_[i]
            }

        params_detail[f"w{self.n_components-1}"] = {
            "pindex": self.n_components*2-1,
            "val": 1-np.sum(self.weights_)
        }

        estimate = np.zeros(self.n_components*2)
        estimate[:self.n_components] = self.probs_
        estimate[self.n_components:] = self.weights_
        return estimate, params_detail
    
    def fisher_matrix(self, X):
        d = self.n_components*2-1 # the size of the parameter space. -1 since one mixing component is a function of the others
        m = np.zeros((d, d))
        for i in range(d):
            for j in range(d):
                m[i,j]=self._fisher_tensor(i=i, j=j)(x=X).sum()

        return m

    def _fisher_tensor(self, i, j):
        """Fixed version of fisher_tensor with correct parameter indexing."""
        weights = self.weights_
        probs = self.probs_
        n_components = self.n_components
        n_trials = self.n_trials
        
        def __partial_logp_partial_prob(x, prob_index):
            """Derivative w.r.t. probability parameter prob_index."""
            binomial_term = weights[prob_index] * binom.pmf(k=x, n=n_trials, p=probs[prob_index])
            score_factor = (x/probs[prob_index] - (n_trials-x)/(1-probs[prob_index]))
            result= (1/mixbinom_logpmf(n_trials=n_trials, 
                                      weights=self.weights_, 
                                      probs=self.probs_, 
                                      flat=True)(x)) * binomial_term * score_factor
            return result
        
        def __partial_logp_partial_weight(x):
            """Derivative w.r.t. the FREE weight parameter (w_0)."""
            # Since w_1 = 1 - w_0, derivative w.r.t. w_0 is:
            # d/dw_0 log[w_0*p_0 + (1-w_0)*p_1] = (p_0 - p_1) / mixture_prob
            prob_diff = (binom.pmf(k=x, n=n_trials, p=probs[0]) - 
                         binom.pmf(k=x, n=n_trials, p=probs[1]))
            result = prob_diff / mixbinom_logpmf(n_trials=n_trials, 
                                               weights=self.weights_, 
                                               probs=self.probs_, 
                                               flat=True)(x)
            return result
        
        def __tensor(x):
            # Parameter ordering: [p_0, p_1, w_0]
            
            # Determine parameter types for i and j
            is_i_prob = i < n_components  # i=0,1 are probabilities
            is_j_prob = j < n_components  # j=0,1 are probabilities
            
            if is_i_prob and is_j_prob:
                # Both are probability parameters
                if i == j:
                    # Diagonal: same probability parameter
                    return __partial_logp_partial_prob(x, i)**2
                else:
                    # Off-diagonal: different probability parameters
                    return (__partial_logp_partial_prob(x, i) * 
                            __partial_logp_partial_prob(x, j))
            
            elif is_i_prob and not is_j_prob:
                # i is probability, j is weight
                return (__partial_logp_partial_prob(x, i) * 
                        __partial_logp_partial_weight(x))
            
            elif not is_i_prob and is_j_prob:
                # i is weight, j is probability  
                return (__partial_logp_partial_weight(x) * 
                        __partial_logp_partial_prob(x, j))
            
            else:
                # Both are weight parameters (only one free weight for 2-component)
                return __partial_logp_partial_weight(x)**2
        
        return np.vectorize(__tensor)
        
    # def aic(self, X):
    #     """Akaike Information Criterion for the current model on the input X."""
    #     return -2 * self.score(X) * len(X) + 2 * (2 * self.n_components - 1)
    
    # def bic(self, X):
    #     """Bayesian Information Criterion for the current model on the input X."""
    #     return (-2 * self.score(X) * len(X) + 
    #             (2 * self.n_components - 1) * np.log(len(X)))
