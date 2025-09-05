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

    def diff_prob(self, x, probs, weights, pindex):
        return weights[pindex]*(x[:, 0]/probs[index]-(n_trials-x[:, 0])/(1-probs[pindex]))*binom.pmf(x[:, 0], n=n_trials, p=probs[pindex])

    def diff2_prob(self, x, probs, weights, pindex, n_trials):
        return weights[pindex]*(-x/(probs[index])^2-(n_trials-x)/(1-probs[pindex])^2+
                                (x/probs[index]-(n_trials-x)/(1-probs[pindex]))^2
                               )*binom.pmf(x, n=n_trials, p=probs[pindex])

    def diff2_weight(self, x, probs, weights, pindex, n_trials):
        return 0

    def fisher_metric(self, n_components, n_trials, probs, weights):
        def tensor(x):
            d = n_components*2-1 # the size of the parameter space. -1 since one mixing component is a function of the others
            m = np.zeros((dim, dim))
            for i in range(n_components):
                for j in range(n_components):
                    if i==j: # diagonal elements
                        if i>d/2: # weights second derivative
                            m[i,j] = self.diff2_weight(x=x, n_trials=n_trials, probs=probs, weights=weights)
                        else: # probs second derivative
                            m[i,j] = self.diff2_prob(x=x, n_trials=n_trials, probs=probs, weights=weights)
                    elif j == i+math.ceiling(d/2): # non-diagonal elements, when the weight index match the probs index
                        m[i,j]=(x/probs[j]-(n_trials-x)/(1-probs[j]))*binom.pmf(x, n=n_trials, p=probs[j])
                    else:
                        m[i,j]=0
            return m
        return tensor
        
    def fisher_matrix(self, X, probs, weights):
        if not probs:
            probs = self.probs_

        if not weights:
            weights = self.weights_
        
        sort_indices = np.argsort(probs)[::-1]
        probs = probs[sorted_indices]
        weights = weights[sorted_indices]
        fisher_tensor = self.fisher_metric(n_components=self.n_components, n_trials=self.n_trials, probs=probs, weights=weights)
        result=[]
        for k in x[:, 0]: # compute matrix over the sample and sum
            result.append(fisher_tensor(k))

        print(result)
        return np.sum(result)

    def compute_fisher_confidence_region(self, X, confidence_level=0.95):
        """
        Compute 95% confidence region using observed Fisher information.
        
        Parameters:
        -----------
        X : numpy array
            observations
        confidence_level : float
            Confidence level (default 0.95)
        
        Returns:
        --------
        dict with covariance matrix, ellipse parameters, etc.
        """
        observed_fisher_mat = self.fisher_matrix(X=X, probs=self.probs_, weights=self.weights_)
        
        # CORRECTION: The asymptotic covariance is J(θ̂)^(-1), not J(θ̂)^(-1)/√n
        # The Fisher information already scales with sample size
        try:
            # Invert Fisher information to get covariance matrix
            cov_matrix = np.linalg.inv(observed_fisher_mat)
            
            # Chi-square critical value for d parameters
            d = len(self.n_components*2-1)
            chi2_critical = chi2.ppf(confidence_level, df=d)
            
            return {
                'cov_matrix': cov_matrix,
                'chi2_critical': chi2_critical,
                'confidence_level': confidence_level,
                'standard_errors': np.sqrt(np.diag(cov_matrix))
            }
            
        except np.linalg.LinAlgError as e:
            raise IllegalException(f"Warning: Fisher information matrix is singular!, cause {e}")
        
    # def aic(self, X):
    #     """Akaike Information Criterion for the current model on the input X."""
    #     return -2 * self.score(X) * len(X) + 2 * (2 * self.n_components - 1)
    
    # def bic(self, X):
    #     """Bayesian Information Criterion for the current model on the input X."""
    #     return (-2 * self.score(X) * len(X) + 
    #             (2 * self.n_components - 1) * np.log(len(X)))
