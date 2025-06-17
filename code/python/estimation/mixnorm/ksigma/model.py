import numpy as np
import pymc as pm
from pymc import logp
import pytensor as pt

# inverse temperature factor
def tempered_normal_mixture(beta, data, 
                            n_components=3,
                            weights_prior_alpha=np.full(3, 0.1),
                            mean_prior_mu=pt.tensor.zeros((3,)),
                            mean_prior_cov=pt.tensor.eye(3)*2):
    if not beta:
        likelihood_power = 1/np.log(len(data))
    else:
        likelihood_power = beta

    model = pm.Model()
    with model:
        # Priors for mixture weights
        # scale = pm.HalfCauchy("unconstrained_weights_scale", beta=2.5)
        # uncons_w_mean = pm.MvNormal("ucons_w_mean", mu=pt.tensor.zeros((n_components,)),
        #                             cov=weights_raw_prior_scale)
        
        # unconstrained_weights = pm.MvNormal("unconstrained_weights",
        #                                     mu=pt.tensor.zeros((n_components,)),
        #                                     # mu=uncons_w_mean,
        #                                     cov=weights_raw_prior_scale)
        
        # weights = pm.Deterministic("weights", pm.math.softmax(unconstrained_weights))
        weights = pm.Dirichlet("weights", a=weights_prior_alpha, shape=n_components)
    
        # Priors for component means
        mus = pm.MvNormal("mus", mu=mean_prior_mu, cov=mean_prior_cov)
        like = pm.NormalMixture("like", w=weights, mu=mus, sigma=1)
        # print(f"Type of mixture_dist: {type(mixture_dist)}")
        # print(f"Does mixture_dist have 'logp' attribute? {'logp' in dir(mixture_dist)}")
        tempered_log_likelihood = likelihood_power * logp(like, data)
        pm.Potential("tempered_likelihood", tempered_log_likelihood)
        # metropolis_step_mu = pm.Metropolis(vars=[weights])
    return model
