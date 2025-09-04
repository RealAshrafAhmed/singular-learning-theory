import numpy as np
from scipy.stats import binom
from scipy.special import logsumexp
from sklearn.mixture._base import BaseMixture
from sklearn.utils.validation import check_X_y, check_array
from sklearn.utils import check_random_state
import warnings
from joblib import Parallel, delayed
import multiprocessing as mp
from sklearn_ext.random import dirichlet
from sklearn.utils._array_api import (
    get_namespace,
    get_namespace_and_device
)


        
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
        probs_init=None,
        reg_prob=1e-6,
        n_components=1,
        n_trials=None,
        *,
        tol=1e-3,
        reg_covar=1e-6,
        max_iter=100,
        n_init=1,
        weights_init=None,
        random_state=None,
        warm_start=False,
        verbose=0,
        verbose_interval=10,
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
    
    def _check_parameters(self, X):
        """Check input parameters."""
        if self.n_components < 1:
            raise ValueError("n_components must be >= 1")

        if self.max_iter < 1:
            raise ValueError("max_iter must be >= 1")

        if self.tol <= 0:
            raise ValueError("tol must be > 0")
    
    def _initialize(self, X, resp, xp=None):
        """Initialization of the Binomial mixture parameters.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)

        resp : array-like of shape (n_samples, n_components)
        """
        n_samples = X.shape[0]
        
        # Initialize mixing weights uniformly
        self.weights_ = dirichlet(k=self.n_components, n_samples=1, alpha=1.0, xp=xp)

        random_state = check_random_state(self.random_state)
        # Init probs randomly
        self.probs_ = random_state.uniform(0.1, 0.9, self.n_components)
        
        # Ensure probabilities are in the regularized range
        self.probs_ = np.clip(self.probs_, self.reg_prob, 1 - self.reg_prob)

    def _initialize_parameters(self, X, random_state, xp=None):
        # If all the initial parameters are all provided, then there is no need to run
        # the initialization.
        # compute_resp = (
        #     self.weights_init is None
        #     or self.probs_init is None
        # )
        # if compute_resp:
        #     super()._initialize_parameters(X, random_state)
        # else:
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
        n_trials = X[:,1]
        
        for k in range(self.n_components):            
            # Compute log probability for binomial distribution
            log_prob[:, k] = binom.logpmf(obs, n_trials, self.probs_[k])
        
        return log_prob

    def _compute_lower_bound(self, _, log_prob_norm):
        return log_prob_norm

    def _enforce_ordering(self):
        """Simple function to enforce p1 >= p2 >= ... >= pK ordering."""
        # Get sorted indices
        sort_indices = np.argsort(self.probs_)[::-1]
        self.probs_ = self.probs_[sort_indices]
        self.weights_ = self.weights_[sort_indices]

    def _m_step(self, X, log_resp):
        """Maximization step."""
        n_samples = X.shape[0]

        resp = np.exp(log_resp)
        # Update mixing weights
        resp_sum = np.sum(resp, axis=0) + 10 * np.finfo(resp.dtype).eps # ensure some numerical stability
        self.weights_ = resp_sum / n_samples
        
        # Update success probabilities
        for k in range(self.n_components):
            # Weighted maximum likelihood estimation
            # Different number of trials for each sample
            weighted_successes = np.sum(resp[:, k] * X[:,0])
            weighted_trials = np.sum(resp[:, k] * X[:,1])
            
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

    # def aic(self, X):
    #     """Akaike Information Criterion for the current model on the input X."""
    #     return -2 * self.score(X) * len(X) + 2 * (2 * self.n_components - 1)
    
    # def bic(self, X):
    #     """Bayesian Information Criterion for the current model on the input X."""
    #     return (-2 * self.score(X) * len(X) + 
    #             (2 * self.n_components - 1) * np.log(len(X)))

# def fit_new(X, **kwargs):
#     # Fit main GMM
#     bmm = SingleBinomialMixture(**kwargs)
#     bmm.fit(X)
    
#     return  {
#         "probs": bmm.probs_.flatten(),
#         "weights": bmm.weights_.flatten(),
#         "lower_bound": bmm.lower_bound_,
#         "converged": bmm.converged_
#     }

# class BinomialMixture():
#     def __init__(self, n_init=1, parallelism_n_jobs=-1, parallelism_verbose=1, **kwargs):
#         self.n_init=n_init
#         self.parallelism_verbose=parallelism_verbose
#         self.parallelism_n_jobs=parallelism_n_jobs
#         self.kwargs=kwargs

#     def fit(self, X):
#         # results = Parallel(n_jobs=self.parallelism_n_jobs, verbose=self.parallelism_verbose)(
#         #     delayed(fit_new)(X, **self.kwargs) 
#         #     for i in range(self.n_init)
#         # )
#         results = []
#         for i in range(self.n_init):
#             # Fit main GMM
#             print(self.kwargs)
#             bmm = SingleBinomialMixture(**self.kwargs)
#             bmm.fit(X)
            
#             results.append({
#                 "probs": bmm.probs_.flatten(),
#                 "weights": bmm.weights_.flatten(),
#                 "lower_bound": bmm.lower_bound_,
#                 "converged": bmm.converged_
#             })

#         # choose the best result
#         best_lower_bound=None
#         best_probs=None
#         best_weights=None
#         for i in range(len(results)):
#             if best_lower_bound is None or lower_bound > best_lower_bound:
#                 self.probs_ = results[i]["probs"]
#                 self.weights_ = results[i]["weights"]
#                 self.lower_bound_ = results[i]["lower_bound"]
#                 self.converged_ = results[i]["converged"]



                
        