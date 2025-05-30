import numpy as np
import pymc as pm
from pymc import logp

# inverse temperature factor
def tempered_normal_mixture(beta, data, n_components=3):
    if not beta:
        likelihood_power = 1/np.log(len(data))
    else:
        likelihood_power = beta

    basic_model = pm.Model()
    with basic_model:
        # Priors for mixture weights
        # Dirichlet distribution ensures weights sum to 1
        weights = pm.Dirichlet("rho", np.ones(n_components))
    
        # Priors for component means
        mus = pm.Normal("mu",
                        mu=0,  # Centered around the data mean
                        sigma=10,        # Broad prior
                        shape=n_components)
    
        # X_obs = pm.NormalMixture("X_obs", 
                                 # w=weights, 
                                 # mu=mus,
                                 # sigma=np.ones(3),
                                 # observed=x_data,
                                 # weight=c/np.log(len(x_data))
    
        like = pm.NormalMixture("like", w=weights, mu=mus, sigma=1)
        # print(f"Type of mixture_dist: {type(mixture_dist)}")
        # print(f"Does mixture_dist have 'logp' attribute? {'logp' in dir(mixture_dist)}")
        tempered_log_likelihood = likelihood_power * logp(like, data)
        pm.Potential("tempered_likelihood", tempered_log_likelihood)
    return basic_model
