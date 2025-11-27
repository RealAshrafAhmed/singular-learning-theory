import numpy as np
import pymc as pm
from pymc import logp
<<<<<<< HEAD
<<<<<<< HEAD:code/python/estimation/mixnorm/ksigma/model.py
=======
>>>>>>> 276c671 (typer command)
import pytensor as pt

# inverse temperature factor
def tempered_normal_mixture(beta, data, 
                            n_components=3,
                            weights_prior_alpha=np.full(3, 0.1),
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 63a42fa (added pymc logging per sample)
                            mean_prior_mu=pt.tensor.zeros((3,)),
                            mean_prior_cov=pt.tensor.eye(3)*2):
=======

# inverse temperature factor
def tempered_normal_mixture(beta, data, n_components=3):
>>>>>>> 73837af (tempered normal mixture):code/python/normal-mixtures-1d/known-sigma/model.py
=======
                            mean_prior_cov=pt.tensor.eye(3)*2):
>>>>>>> 276c671 (typer command)
    if not beta:
        likelihood_power = 1/np.log(len(data))
    else:
        likelihood_power = beta

<<<<<<< HEAD
<<<<<<< HEAD:code/python/estimation/mixnorm/ksigma/model.py
=======
>>>>>>> 276c671 (typer command)
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
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
        mus = pm.MvNormal("mus", mu=mean_prior_mu, cov=mean_prior_cov)
=======
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
    
>>>>>>> 73837af (tempered normal mixture):code/python/normal-mixtures-1d/known-sigma/model.py
=======
        mus = pm.MvNormal("mu", 
=======
        mus = pm.MvNormal("mus", 
>>>>>>> 30724ab (fixed output)
                          mu=pt.tensor.zeros((n_components,)),
                          cov=mean_prior_cov)
    
>>>>>>> 276c671 (typer command)
=======
        mus = pm.MvNormal("mus", mu=mean_prior_mu, cov=mean_prior_cov)
>>>>>>> 63a42fa (added pymc logging per sample)
        like = pm.NormalMixture("like", w=weights, mu=mus, sigma=1)
        # print(f"Type of mixture_dist: {type(mixture_dist)}")
        # print(f"Does mixture_dist have 'logp' attribute? {'logp' in dir(mixture_dist)}")
        tempered_log_likelihood = likelihood_power * logp(like, data)
        pm.Potential("tempered_likelihood", tempered_log_likelihood)
<<<<<<< HEAD
<<<<<<< HEAD:code/python/estimation/mixnorm/ksigma/model.py
        # metropolis_step_mu = pm.Metropolis(vars=[weights])
    return model
=======
    return basic_model
>>>>>>> 73837af (tempered normal mixture):code/python/normal-mixtures-1d/known-sigma/model.py
=======
        # metropolis_step_mu = pm.Metropolis(vars=[weights])
    return model
>>>>>>> 276c671 (typer command)
