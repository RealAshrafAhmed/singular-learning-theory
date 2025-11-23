import numpy as np
from scipy.stats import poisson, chi2
from scipy.special import logsumexp
from sklearn.mixture._base import BaseMixture
from sklearn.utils.validation import check_X_y, check_array
from sklearn.utils import check_random_state
import warnings
from joblib import Parallel, delayed
import multiprocessing as mp
from sklearn_extensions.random import dirichlet
from sklearn_extensions.mixture._fun import __set_param, __set_params
from sklearn.utils._array_api import (
    get_namespace,
    get_namespace_and_device
)


def create_kl_on_grid(n_components, n_trials, truth, fixed_params, x_param, y_param):
    # n_components=n_components
    n_trials = n_trials
    # truth = truth
    # fixed_params = fixed_params
    x_param = x_param
    y_param = y_param
    k_values = np.arange(0, n_trials)
    truth_probs = truth[:n_components]
    truth_weights = truth[n_components:]
    model_probs = np.zeros(n_components)
    model_weights = np.zeros(n_components)
    model_probs, model_weights = __set_params(model_probs, model_weights, fixed_params)

    def fun(X, Y):
        # Z = np.zeros(len(X))
        Z = np.zeros(X.shape)
        for index in range(Z.shape[0]):
            for sindex in range(Z.shape[1]):
                # h[index, sindex]=ell(xindex=i, yindex=j, x=, y=)
            # print(X[index])
                m_probs, m_weights = __set_param(model_probs, model_weights,  x_param["name"], x_param["index"], X[index, sindex])
                m_probs, m_weights = __set_param(m_probs, m_weights,  y_param["name"], y_param["index"], Y[index, sindex])
                m_weights[n_components-1] = 1-np.sum(m_weights[:-1])
                # print(model_weights)
                k_data = np.zeros(len(k_values))
                for i, k in enumerate(k_values):
                    log_q = binommix_logpmf(k, n_trials, truth_weights, truth_probs)
                    log_p = binommix_logpmf(k, n_trials, m_weights, m_probs)
                    k_data[i] = np.exp(log_q)*(log_q-log_p)
        
                Z[index, sindex] = np.sum(k_data)

        return Z

    return fun


def mixpoisson_logpmf(x, params_mus, params_weights, flat=True):
    """
    Compute a probability mass function of binomial mixture
    x: samples
    n_trials: number of binomial trials
    params_weights: the mixing weights of each component, must add to 1
    params_probs: the probability of each component in the mixture
    log: return log probability, default is True
    flat: returna the probilities not just each component
    """
    x = np.atleast_1d(x)
    n_components = len(params_mus)
    size = len(x)
    result = np.zeros((size, n_components)) # +1 to capture the mixture
    
    for pindex in range(len(params_weights)):
        result[:, pindex] = poisson.logpmf(k=x, mu=params_mus[pindex])
        negative_weights = np.any(params_weights < 0)
        if(negative_weights):
            print(f"weightts={params_weights} has ne gative values, x={s}")

        total = logsumexp(result + np.log(params_weights), axis=1)
        return np.sum(total)


def mixpoisson_sample(size=1, mus=np.ones(2), weights=np.ones(2)*0.5, random_state=None):
    """
    Generate random samples from the fitted mixture model.
    
    Parameters
    ----------
    size : int, default=1
        Number of samples to generate.

    mus : int, default=1
        rates parameters

    weights : int, default=1
        mixing weights

    random_state : int, RandomState instance or None, default=None
        Random seed.
        
    Returns
    -------
    X : array-like of shape (size,)
        Generated samples.
    y : array-like of shape (size,)
        Component labels for each sample.
    """
    random_state = check_random_state(random_state)
    n_components = len(mus)
    assert len(mus) == len(weights), "Number of components in mus and weights must be the same"
    
    # Sample component assignments
    component_samples = random_state.choice(n_components, 
                                            size=size, 
                                            p=weights)
    
    # Generate samples from each component
    samples = np.zeros(size, dtype=int)
    
    for k in range(n_components):
        mask = component_samples == k
        n_component_samples = np.sum(mask)
        if n_component_samples > 0:
            samples[mask] = random_state.poisson(mus[k], size=n_component_samples)
    
    return samples, component_samples


def _estimate_mixpoisson_parameters(n_components, X, resp, random_state, xp):
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
    mus_ = probs_/(1-probs_)
    
    # Ensure probabilities are in the regularized range
    mus_ = np.clip(mus_, np.finfo(float).tiny, np.finfo(np.float32).max)
    return weights_, mus_

        
class PoissonMixture(BaseMixture):
    """
    Poisson Mixture Model using Expectation-Maximization algorithm.
    
    This class implements a mixture of poisson distributions that follows
    scikit-learn's estimator interface.
    
    Parameters
    ----------
    n_components : int, default=2
        The number of mixture components.

    max_iter : int, default=100
        Maximum number of EM iterations.
        
    tol : float, default=1e-6
        Convergence tolerance.
        
    random_state : int, RandomState instance or None, default=None
        Random seed for initialization.
        
    init_params : str, default='random'
        Method for initialization ('random' or 'kmeans').
    
    Attributes
    ----------
    weights_ : array-like of shape (n_components,)
        Mixing weights for each component.
    mus_ : array-like of shape (n_components,)
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
    >>> from sklearn_ext.mixture import PoissonMixture
    >>> X = np.array([[8], [18], [15]])
    >>> gm = PoissonMixture(n_components=2, random_state=0).fit(X)
    >>> gm.mus_
    array([[.3,  .1],
           [.2,  .6]])
    >>> gm.predict([[0, 0], [12, 3]])
    array([1, 0])
    """
    def __init__(
        self,
        mus_init=None,
        n_components=1,
        *,
        tol=1e-5,
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
            max_iter=max_iter,
            reg_covar=None,
            n_init=n_init,
            init_params="random",
            random_state=random_state,
            warm_start=warm_start,
            verbose=verbose,
            verbose_interval=verbose_interval,
        )
        self.mus_init=mus_init
        self.weights_init=weights_init

        if (
            weights_init is not None
            and probs_init is not None
        ):
            self.mus_ = self.mus_init
            self.weights_ =self.weights_init

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

        self.weights_, self.probs_ = _estimate_poisson_parameters(n_components=self.n_components, 
                                                                  X=X,
                                                                  resp=resp, 
                                                                  random_state=self.random_state, 
                                                                  xp=xp)

    def _initialize_parameters(self, X, random_state, xp=None):
        # If all the initial parameters are all provided, then there is no need to run
        # the initialization.
        compute_resp = (
            self.weights_init is None
            or self.mus_init is None
        )
        if compute_resp:
            self.weights_, self.mus_ = _estimate_mixpoisson_parameters(n_components=self.n_components, 
                                                                     X=X, 
                                                                     resp=None,
                                                                     random_state=random_state, 
                                                                     xp=xp)
        else:
            self._initialize(X, None, xp=xp)

    def _get_parameters(self):
        return (
            self.weights_,
            self.mus_
        )

    def _set_parameters(self, params, xp=None):
        xp, _, device_ = get_namespace_and_device(params, xp=xp)
        (
            self.weights_,
            self.mus_
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
        
        for k in range(self.n_components):            
            log_prob[:, k] = poisson.logpmf(obs, mu=self.mus_[k])
        
        return log_prob

    def _compute_lower_bound(self, _, log_mu_norm):
        return log_mu_norm

    def _enforce_ordering(self):
        if not self.enforce_ordering:
            return
        """Simple function to enforce p1 >= p2 >= ... >= pK ordering."""
        # Get sorted indices
        sort_indices = np.argsort(self.mus_)[::-1]
        self.mus_ = self.mus_[sort_indices]
        self.weights_ = self.weights_[sort_indices]
        
    def _m_step(self, X, log_resp):
        """Maximization step."""
        old_mus = self.mus_
        old_weights = self.weights_
        n_samples = X.shape[0]
        # Update mixing weights
        resp = np.exp(log_resp)
        # resp_sum = np.sum(resp, axis=0) + 10 * np.finfo(resp.dtype).eps # ensure some numerical stability
        resp_sum = np.sum(resp, axis=0)
        self.weights_ = resp_sum / n_samples
        # Update rates
        for k in range(self.n_components):
            # Weighted maximum likelihood estimation
            weighted_means = np.sum(resp[:, k] * X[:,0])
            
            # if means > 0:
            self.mus_[k] = weighted_means / np.sum(resp[:,k])
            # else:
                # self.mus_[k] = 0.5  # Default value

            self.mus_[k] = np.clip(self.mus_[k], np.finfo(float).tiny, np.finfo(np.float32).max)

        # print(f"changed weights from {old_weights} to {self.weights_}")
        # print(f"changed mus from {old_mus} to {self.mus_}")
        self._enforce_ordering()
    
    def sample(self, n_samples=1, random_state=None):
        return mixpoisson_sample(mus=self.mus_, weights=self.weights_, size=n_samples, random_state=random_state)

    def point_estimate(self):
        self._enforce_ordering()
        params_detail = {}
        for i in range(self.n_components):
            params_detail[f"m{i}"] = {
                "pindex": i,
                "val": self.mus_[i]
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
        estimate[:self.n_components] = self.mus_
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
            result= (1/binommix_logpmf(x=x, 
                                      n_trials=n_trials, 
                                      params_weights=self.weights_, 
                                      params_probs=self.probs_, 
                                      flat=True)) * binomial_term * score_factor
            return result
        
        def __partial_logp_partial_weight(x):
            """Derivative w.r.t. the FREE weight parameter (w_0)."""
            # Since w_1 = 1 - w_0, derivative w.r.t. w_0 is:
            # d/dw_0 log[w_0*p_0 + (1-w_0)*p_1] = (p_0 - p_1) / mixture_prob
            prob_diff = (binom.pmf(k=x, n=n_trials, p=probs[0]) - 
                         binom.pmf(k=x, n=n_trials, p=probs[1]))
            result = prob_diff / binommix_logpmf(x=x, 
                                               n_trials=n_trials, 
                                               params_weights=self.weights_, 
                                               params_probs=self.probs_, 
                                               flat=True)
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
