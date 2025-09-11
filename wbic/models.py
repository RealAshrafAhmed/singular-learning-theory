import numpy as np
import pymc as pm
from pymc import logp
import pytensor as pt

def tempered_normal_mixture_known_variance(data, beta=1.0, 
                            n_components=3,
                            weights_prior_alpha=np.full(3, 0.1),
                            mean_prior_mu=pt.tensor.zeros((3,)),
                            mean_prior_cov=pt.tensor.eye(3)*2):
    """
    Create a tempered Gaussian mixture model with known variance components.
    
    This function creates a Bayesian mixture model where each component has unit variance
    (sigma=1) and applies tempering to the likelihood via a power parameter beta.
    
    Parameters:
    -----------
    data : array-like
        Observed data (1D array)
    beta : float, default=1.0
        Inverse temperature parameter for tempering. If beta=0, uses 1/log(n_data).
        Lower values flatten the likelihood, higher values sharpen it.
    n_components : int, default=3
        Number of mixture components
    weights_prior_alpha : array-like, default=np.full(3, 0.1)
        Dirichlet prior parameters for mixture weights
    mean_prior_mu : tensor, default=zeros(3)
        Prior mean for component means (multivariate normal)
    mean_prior_cov : tensor, default=2*I_3
        Prior covariance matrix for component means
        
    Returns:
    --------
    model : pymc.Model
        PyMC model object ready for sampling
        
    Notes:
    ------
    - All components have fixed unit variance (sigma=1)
    - Uses tempered likelihood: L(data|params)^beta
    - Component means follow multivariate normal prior
    """
    if not beta:
        likelihood_power = 1/np.log(len(data))
    else:
        likelihood_power = beta

    model = pm.Model()
    with model:
        # Priors for mixture weights
        weights = pm.Dirichlet("weights", a=weights_prior_alpha, shape=n_components)
    
        # Priors for component means
        mus = pm.MvNormal("mus", mu=mean_prior_mu, cov=mean_prior_cov)

        likelihood = pm.NormalMixture("likelihood", w=weights, mu=mus, sigma=1, observed=data)
        tempered_log_likelihood = likelihood_power * logp(likelihood, data)
        pm.Potential("tempered_likelihood", tempered_log_likelihood)
    return model

def tempered_normal_mixture_model(data, beta=1.0, n_components=2, estimate_sigma=True):
    """
    Create a tempered Gaussian mixture model with flexible variance estimation.
    
    This function creates a Bayesian mixture model where component variances can either
    be estimated individually for each component or shared across components. The 
    likelihood is tempered using the beta parameter.
    
    Parameters:
    -----------
    data : array-like
        Observed data (1D array)
    beta : float, default=1.0
        Inverse temperature parameter for tempering the likelihood.
        Lower values flatten the likelihood, higher values sharpen it.
    n_components : int, default=2
        Number of mixture components
    estimate_sigma : bool, default=True
        If True, estimate individual sigmas for each component.
        If False, use shared sigma across all components.
        
    Returns:
    --------
    model : pymc.Model
        PyMC model object ready for sampling
        
    Notes:
    ------
    - Component means are given broad normal priors centered on data range
    - Component variances use half-normal priors
    - Uses tempered likelihood: L(data|params)^beta
    """
    n_data = len(data)
    
    # Data statistics for priors
    data_min, data_max = np.min(data), np.max(data)
    data_range = data_max - data_min
    data_center = (data_min + data_max) / 2
    
    with pm.Model() as model:
        # Mixture weights (must sum to 1)
        if n_components == 1:
            weights = pt.tensor.constant([1.0])
        else:
            weights = pm.Dirichlet("weights", a=np.ones(n_components), shape=n_components)
        
        # Component means
        if n_components == 1:
            mus = pm.Normal("mus", mu=data_center, sigma=data_range, shape=1)
        else:
            mus = pm.Normal("mus", 
                           mu=data_center, 
                           sigma=data_range, 
                           shape=n_components)
        
        # Component standard deviations
        if estimate_sigma:
            sigmas = pm.HalfNormal("sigmas", sigma=data_range/4, shape=n_components)
        else:
            sigma_shared = pm.HalfNormal("sigma", sigma=data_range/4)
            sigmas = pt.tensor.full(n_components, sigma_shared)
            
        # Mixture likelihood
        if n_components == 1:
            likelihood = pm.Normal("likelihood", mu=mus[0], sigma=sigmas[0], observed=data)
        else:
            likelihood = pm.NormalMixture("likelihood", w=weights, mu=mus, sigma=sigmas, observed=data)
        
        # Apply tempering
        tempered_log_likelihood = beta * logp(likelihood, data)
        pm.Potential("tempered_likelihood", tempered_log_likelihood)
        
    return model

def tempered_binomial_mixture_model(data, n_trials=10, beta=1.0, n_components=2):
    """
    Create a tempered binomial mixture model.
    
    This function creates a Bayesian mixture model for binomial data where each
    component has its own success probability parameter. The likelihood is 
    tempered using the beta parameter.
    
    Parameters:
    -----------
    data : array-like
        Observed success counts (1D array of integers)
    n_trials : int or array-like
        Number of trials for each observation. If int, assumes same number
        of trials for all observations.
    beta : float, default=1.0
        Inverse temperature parameter for tempering the likelihood.
        Lower values flatten the likelihood, higher values sharpen it.
    n_components : int, default=2
        Number of mixture components
        
    Returns:
    --------
    model : pymc.Model
        PyMC model object ready for sampling
        
    Notes:
    ------
    - Each component has its own success probability with Beta(1,1) prior
    - Mixture weights follow uniform Dirichlet prior
    - Uses tempered likelihood: L(data|params)^beta
    """
    n_data = len(data)
    
    with pm.Model() as model:
        # Mixture weights (must sum to 1)
        if n_components == 1:
            weights = pt.tensor.constant([1.0])
        else:
            weights = pm.Dirichlet("weights", a=np.ones(n_components), shape=n_components)
        
        # Component success probabilities
        ps = pm.Beta("ps", alpha=1, beta=1, shape=n_components)
        
        # Mixture likelihood
        likelihood = pm.BinomialMixture("likelihood", w=weights, n=n_trials, p=ps)
        
        # Apply tempering
        tempered_log_likelihood = beta * logp(like, data)
        pm.Potential("tempered_likelihood", tempered_log_likelihood)
        
    return model

def tempered_poisson_mixture_model(data, beta=1.0, n_components=2):
    """
    Create a tempered Poisson mixture model.
    
    This function creates a Bayesian mixture model for count data where each
    component follows a Poisson distribution with its own rate parameter.
    The likelihood is tempered using the beta parameter.
    
    Parameters:
    -----------
    data : array-like
        Observed count data (1D array of non-negative integers)
    beta : float, default=1.0
        Inverse temperature parameter for tempering the likelihood.
        Lower values flatten the likelihood, higher values sharpen it.
    n_components : int, default=2
        Number of mixture components
        
    Returns:
    --------
    model : pymc.Model
        PyMC model object ready for sampling
        
    Notes:
    ------
    - Each component has its own rate parameter with Gamma prior
    - Gamma prior parameters chosen based on data statistics
    - Mixture weights follow uniform Dirichlet prior
    - Uses tempered likelihood: L(data|params)^beta
    """
    n_data = len(data)
    
    # Data statistics for priors
    data_mean = np.mean(data)
    data_std = np.std(data)
    
    # Gamma prior parameters (method of moments)
    # For Gamma(alpha, beta): mean = alpha/beta, var = alpha/beta^2
    if data_std > 0:
        gamma_beta = data_mean / (data_std**2)
        gamma_alpha = data_mean * gamma_beta
    else:
        # If no variance, use weakly informative prior
        gamma_alpha = data_mean + 0.1
        gamma_beta = 1.0
    
    with pm.Model() as model:
        # Mixture weights (must sum to 1)
        if n_components == 1:
            weights = pt.tensor.constant([1.0])
        else:
            weights = pm.Dirichlet("weights", a=np.ones(n_components), shape=n_components)
        
        # Component rate parameters
        lambdas = pm.Gamma("lambdas", alpha=gamma_alpha, beta=gamma_beta, shape=n_components)
        
        # Mixture likelihood
        like = pm.PoissonMixture("likelihood", w=weights, mu=lambdas)
        
        # Apply tempering
        tempered_log_likelihood = beta * logp(like, data)
        pm.Potential("tempered_likelihood", tempered_log_likelihood)
        
    return model