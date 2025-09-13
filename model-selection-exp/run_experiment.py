import numpy as np
import pymc as pm 
import arviz as az
import pandas as pd

import os
import time
from pathlib import Path
import logging

from expclass import MixtureModelExperiment
from generate_data import generate_gmm_samples
"""
EXPERIMENT PARAMETERS:
"""

experiment_folder = "test_experiment"

n_data = 500
n_repetitions = 10
true_components = 3
max_components = 5

sampling_params = {
    'draws':5000,
    'tune':2500,
}

"""
RUN EXPERIMENT:
"""

# Set up folder structure
def setup_experiment(namespace):
    """Set up directory structure and logging for experiment."""
    base_path = Path(namespace)
    
    (base_path / "data").mkdir(parents=True, exist_ok=True)
    (base_path / "figs").mkdir(parents=True, exist_ok=True)
    (base_path / "logs").mkdir(parents=True, exist_ok=True)
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(base_path / "logs" / "experiment.log"),
            logging.StreamHandler()
        ]
    )
    
    return logging.getLogger(__name__)

print(f"Setting up experiment in folder {experiment_folder}")
setup_experiment(experiment_folder)

## initialize everything
for i in range(n_repetitions):
    observations, true_params = generate_gmm_samples(n_data, true_components)

    experiment = MixtureModelExperiment(
        observations,
        max_components
    )

    experiment.compute_sbic()
    experiment.compute_wbic(sampling_params)
    experiment.compute_bic()

    results = {
        'Components':[],
        'BIC':[],
        'sBIC':[],
        'WBIC':[],
        'LOO':[],
    }

    k=1
    for model in experiment.submodels:
        results['Components'].append(k)
        results['BIC'].append(model.bic)
        results['sBIC'].append(model.sbic)
        results['WBIC'].append(model.wbic)
        idata = model.inference_data
        idata.add_groups(
            log_likelihood={"log_likelihood":idata.posterior['log_likelihood'].values}
        )
        elpd = az.loo(idata)
        results['LOO'].append(elpd.elpd_loo)
        k += 1

    df = pd.DataFrame(results)
    df.to_csv(experiment_folder+f"/data/trial{i}_information_criteria.csv")
    print(df)

