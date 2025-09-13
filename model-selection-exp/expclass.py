from dataclasses import dataclass
from typing import Optional
import numpy as np
import pymc as pm 
import arviz as az

import models 
import utils
import wbic 

from sklearn.mixture import GaussianMixture
def gaussian_mixture_loglikelihood(observations, n_components):
    gmm = GaussianMixture(n_components=n_components)
    X = observations.reshape(-1,1) #sklearn expects 2d arrays
    gmm.fit(X)
    log_likelihood = np.sum(gmm.score_samples(X))
    return log_likelihood

def marginal_loglikelihood(X, mle, lc):
    m = 1
    return mle - lc*np.log(X.size) + (m-1)*np.log(np.log(X.size)) 

def sbic(X, i, mles, lc, store={}, p=None):
    if p is None:
        # Uniform priors - all p[j] are equal, so they cancel out in ratios
        p_i = 1.0
        p_j = 1.0
    else:
        p_i = p[i]
        p_j = p  # assuming same prior for all j < i for simplicity
    
    if i == 1:
        # Base case: minimal model
        Lii = marginal_loglikelihood(X, mles[i-1], lc[i-1])
        return Lii
    
    # Get Lii for model i
    Lii = marginal_loglikelihood(X, mles[i-1], lc[i-1])
    
    # Compute sum of L[j] * p[j] for all j < i
    sum_Lj_pj = 0
    sum_Lij_Lj_pj = 0
    
    for j in range(1, i):
        if j in store:
            Lj = store[j]
        else:
            Lj = sBIC(X, j, store, p)
            store[j] = Lj
        
        # Get Lij (marginal likelihood of data under model i, constrained to submodel j)
        Lij = marginal_loglikelihood(X, mles[i-1], lc[j-1])
        
        sum_Lj_pj += np.exp(Lj) * p_j
        sum_Lij_Lj_pj += np.exp(Lij + Lj) * p_j
    
    # Quadratic formula coefficients
    a = p_i
    b = -np.exp(Lii) * p_i + sum_Lj_pj
    c = -sum_Lij_Lj_pj
    
    # Solve: a * L[i]^2 + b * L[i] + c = 0
    # Take positive root: L[i] = (-b + sqrt(b^2 - 4ac)) / (2a)
    discriminant = b**2 - 4*a*c
    
    if discriminant < 0:
        raise ValueError(f"Negative discriminant: {discriminant}")
    
    Li_value = (-b + np.sqrt(discriminant)) / (2*a)
    
    if Li_value <= 0:
        # This can happen when the value is like -1e-200...
        Li_value = 1e-200
        #raise ValueError(f"Non-positive probability: {Li_value}")
    
    # Convert back to log space
    Li = np.log(Li_value)
    
    return Li   


class MixtureModelExperiment:
    def __init__(self, observations, n_components, dist='Gaussian'):
        self.submodels = []
        self.observations = observations
        self.idata = []
        self.n_components = n_components
        for k in range(n_components):
            match dist:
                case 'Gaussian':
                    model = models.tempered_gaussian_mixture(
                        observations,
                        n_components = k+1,
                        beta = 1/np.log(observations.size)
                    )
                case _:
                    raise ValueError
            gmm = wbic.BayesianModel(
                model = model,
                observations = observations
            )
            self.submodels.append(gmm)
        self.assign_learning_coefficients()

    def assign_learning_coefficients(self):
        # This formula only valid for Gaussian mixtures.
        k = 1
        for model in self.submodels:
            lc = (self.n_components + 2*k - 1)/2
            model.relative_learning_coefficient = lc
            k += 1
    
    def compute_sbic(self):
        k = 0
        mles = np.zeros(self.n_components)
        lc = np.zeros(self.n_components)
        store = {}
        for model in self.submodels:
            # computing mles in this loop is ok because ith sbic
            # only needs mles for j<=i
            mles[k] = gaussian_mixture_loglikelihood(self.observations, k+1)
            lc[k] = model.relative_learning_coefficient
            store[k+1] = sbic(
                self.observations,
                k+1,
                lc=lc,
                mles=mles,
                store=store,
            )
            model.sbic = store[k+1]
            model.mle = mles[k]
    
    def compute_wbic(self, sampling_params={'draws':10000, 'tune':5000}):
        k=1
        for model in self.submodels:
            print(f"Sampling submodel {k}/{self.n_components}.")
            model.sample(**sampling_params)
            model.wbic = model.WBIC()
            k += 1

    def compute_bic(self):
        k = 1
        for model in self.submodels:
            if model.mle == None:
                model.mle = gaussian_mixture_loglikelihood(self.observations, k)

            d = 2*k + k-1
            model.bic = model.mle - np.log(self.observations.size)*d/2
            k+=1

