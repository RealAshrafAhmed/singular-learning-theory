import numpy as np
from numpy.typing import NDArray
from itertools import permutations


def generate_mixture_permutations(params: NDArray[np.float64]) -> NDArray[np.float64]:
    """
    Generate all permutations of mixture model parameters due to label switching.
    
    Args:
        params: (2k-1,) array [p0, p1, ..., p_{k-1}, w0, w1, ..., w_{k-2}]
                where w_{k-1} = 1 - sum(w_i) is implicit (constrained)
    
    Returns:
        (k!, 2k-1) array of all equivalent parameterizations
    """
    n = len(params)
    k = (n + 1) // 2  # number of components
    
    # Split into component params and weights
    p = params[:k]
    w_partial = params[k:]
    w_full = np.append(w_partial, 1 - w_partial.sum())  # recover constrained weight
    
    result = []
    for perm in permutations(range(k)):
        perm = list(perm)
        new_p = p[perm]
        new_w = w_full[perm][:-1]  # permute then drop last (constrained)
        result.append(np.concatenate([new_p, new_w]))
    
    return result