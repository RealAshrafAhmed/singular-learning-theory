import numpy as np

def dirichlet(k, n_samples=1, alpha=1.0, xp=None):
    """
    Generate probability vectors using Dirichlet distribution
    
    Parameters:
    -----------
    k : int
        Dimension of simplex (number of components)
    n_samples : int
        Number of probability vectors to generate
    alpha : float or array-like
        Concentration parameter(s) for Dirichlet
    random_state : int or None
        Random seed
    
    Returns:
    --------
    prob_vectors : array of shape (n_samples, k)
        Each row sums to 1
    """
    if not xp:
        xp = np
    
    # If alpha is scalar, make it a vector
    if xp.isscalar(alpha):
        alpha = xp.full(k, alpha)
    
    prob_vectors = xp.random.dirichlet(alpha, size=n_samples)
    
    if n_samples == 1:
        return prob_vectors[0]  # Return 1D array for single sample
    
    return prob_vectors