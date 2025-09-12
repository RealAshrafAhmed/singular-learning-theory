import pymc as pm
import numpy as np
import pytensor.tensor as pt

def tempered_gaussian(
        observations,
        beta=1.0,
    ):
    with pm.Model() as model:
        # Priors
        mu = pm.Normal("mu", mu=0, sigma=10)
        sigma = pm.HalfNormal("sigma", sigma=10)

        # This has to be wrapped in pm.Deterministic
        # Otherwise it breaks the tempering for some reason?
        likelihood = pm.Deterministic(
            "log_likelihood", pm.logp(pm.Normal.dist(mu=mu,sigma=sigma), observations)
        )

        # Tempering
        obs_dist = pm.Normal.dist(mu=mu, sigma=sigma)
        tempered_ll = pm.Potential('tempered_ll', beta * pm.logp(obs_dist, observations))
                                  
    return model

def tempered_gaussian_mixture(
        data,
        sigma=None,
        beta=1.0,
        n_components = 2
    ):
    n_data = len(data)
    data_min, data_max = np.min(data), np.max(data)
    data_range = data_max - data_min
    with pm.Model() as model:
        ## Priors
        if n_components == 1:
            weights = pt.constant([1.0])
        else:
            weights = pm.Dirichlet("weights", a=np.ones(n_components), shape=n_components)

        mus = pm.Normal(
            "mus", 
            mu=np.mean(data), 
            sigma=data_range, 
            shape=n_components
        )
        if sigma == None:
            sigmas = pm.HalfNormal(
                "sigmas",
                sigma=data_range/4,
                shape=n_components
            )
        else:
            sigmas = pt.constant([sigma]*n_components)
        
        likelihood = pm.Deterministic(
            "log_likelihood",
            pm.logp(pm.NormalMixture.dist(
                w=weights,
                mu=mus,
                sigma=sigmas,
                ), data
            )
        )

        ## Tempered potential
        obs_pdf  = pm.NormalMixture.dist(
            w=weights,
            mu=mus,
            sigma=sigmas,
        )
        tempered_ll = pm.Potential('tempered_ll', beta * pm.logp(obs_pdf, data))
                                  
    return model