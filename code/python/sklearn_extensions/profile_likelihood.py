import numpy as np
from sklearn.base import clone
from typing import Callable, Tuple, Dict, Any, List
from tqdm.notebook import tqdm


def create_profile_likelihood(model, X, random_state) -> Callable:
    """
    Factory function that creates a profile likelihood evaluator
    
    Parameters:
    -----------
    model : sklearn estimator
        Base model (e.g., GaussianMixture)
    X : array-like
        Training data
    max_iter : int
        Maximum EM iterations
    tol : float
        Convergence tolerance
    n_init : int
        Number of random initializations
        
    Returns:
    --------
    profile_ll : callable
        Function that takes fixed parameters and returns (log_likelihood, optimized_params)
        
    Examples:
    ---------
    >>> from sklearn.mixture import GaussianMixture
    >>> model = GaussianMixture(n_components=2)
    >>> profile_ll = create_profile_likelihood(model, X)
    >>> 
    >>> # Fix weights_[0] = 0.3
    >>> ll, params = profile_ll({'weights_': (0, 0.3)})
    >>> 
    >>> # Fix means_[0, 0] = 2.5
    >>> ll, params = profile_ll({'means_': ((0, 0), 2.5)})
    >>> 
    >>> # Fix entire weights array
    >>> ll, params = profile_ll({'weights_': np.array([0.3, 0.7])})
    """
    
    def profile_ll(fixed_params: List[Tuple[str, int]], fixed_value: np.ndarray):
        """
        Compute profile likelihood with fixed parameters
        
        Parameters:
        -----------
        fixed_params : dict
            Parameters to fix, formats:
            - {'param_name': (index, value)} - fix param_name[index] = value
            - {'param_name': ((i, j), value)} - fix param_name[i, j] = value
            - {'param_name': array} - fix entire parameter array
            
        Returns:
        --------
        log_likelihood : float
            Profile log-likelihood at fixed parameters
        optimized_params : dict
            All optimized model parameters (including fixed ones)
        """
        assert len(fixed_params) == len(fixed_value.shape), "fixed_params must match the column size of fixed_values."
        best_ll = -np.inf
        best_params = None
        
        for _ in range(model.n_init):
            # Clone and initialize
            m = clone(model)
            m._initialize_parameters(X, random_state)

            for i, fixed_param in enumerate(fixed_params):
                param_name, index = fixed_param
                param_array = getattr(m, param_name)
                param_array[index] = fixed_value[i]
            
            # Constrained EM
            prev_ll = -np.inf
            for iteration in range(model.max_iter):
                # E-step
                log_prob_norm, log_resp = m._e_step(X)
                
                # M-step
                m._m_step(X, log_resp)
                
                # Restore fixed parameters
                old_p, _  = m.point_estimate()
                for i, fixed_param in enumerate(fixed_params):
                    param_name, index = fixed_param
                    param_array = getattr(m, param_name)
                    param_array[index] = fixed_value[i]

                # print(f"old_p={np.round(old_p,3)}, new_p={np.round(m.point_estimate()[0],3)}")
                
                # Check convergence
                current_ll = np.sum(log_prob_norm)
                if abs(current_ll - prev_ll) < model.tol:
                    break
                prev_ll = current_ll
            
            # Keep best
            if current_ll > best_ll:
                best_ll = current_ll
                best_params, _ = m.point_estimate()

        return best_ll, best_params
    

    def v_profile_ll(fixed_params: List[Tuple[str, int]], fixed_values: np.ndarray):
        if fixed_values.ndim < 2:
            return profile_ll(fixed_params, fixed_values)

        n_rows, _ = fixed_values.shape
        plls = np.zeros(n_rows)
        fits = []
        for i, fixed_value in tqdm(enumerate(fixed_values), desc=f"profile_ll "):
            pll, fit = profile_ll(fixed_params=fixed_params, fixed_value=fixed_value)
            plls[i] = pll
            fits.append(fit)

        return plls, np.array(fits)

    return v_profile_ll
            
    
    