
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from collections import namedtuple
ModelInits = namedtuple('ModelInits', ['max_iter', 'tol', 'random_state'])


class bootstrapper:
    def __init__(self, model_class, model_inits: ModelInits)
        assert model_class, "model_class must be a type of                  
                 : int =10**6, tol: float =10e-4, , type):
        self.max_iter = max_iter
        self.tol = tol
        self.parallelism_n_jobs = parallelism_n_jobs

    def boot(X,  n_bootstrap=10000, , parallelism_n_jobs: int =-1)
        # compute_mle(x_obs, max_iter=10**6, tol=10e-4):
        # Bootstrap
        # Parallel execution with progress bar
        results = Parallel(n_jobs=parallelism_n_jobs, verbose=0)(
            delayed(resample_and_fit)(X, i, len(x_obs), max_iter, tol) 
            for i in range(n_bootstrap)
        )
    
        # Filter out None results (non-converged)
        converged_results = [r for r in results if r is not None]
    
        # print(f"Converged bootstrap samples: {len(converged_results)}/{n_bootstrap}")
        # print(f"Successful bootstrap samples: {len(p1_boot)}/{n_bootstrap}")
        
        if len(converged_results) < 50:  # Too few successful samples
            print("Warning: Very few successful bootstrap samples!")
            return {"error": "Bootstrap failed - too few successful samples"}
        
        estimates = {"ss": len(x_obs)}
    
        for pname in ["p1", "p2", "rho1", "rho2"]:
            estimates[pname] = {
                "mle": mle[pname],
                "bootstraps": [r[pname] for r in converged_results],
                "ci_lower": np.percentile([r[pname] for r in converged_results], 2.5),
                "ci_upper": np.percentile([r[pname] for r in converged_results], 97.5),
            }
        
        return estimates