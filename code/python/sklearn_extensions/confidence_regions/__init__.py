import numpy as np
import pandas as pd


class ConfidenceInterval:
    def __init__(self, intervals):
        # print(f"intervals={intervals}")
        self.lower_bounds = intervals[:, 0]
        self.upper_bounds = intervals[:, 1]

    def contains(self, values):
        # print(f"values={values}")
        within_intervals = np.all((values >= self.lower_bounds) & (values <= self.upper_bounds))
        return within_intervals


class CRProvider:
    @property
    def type(self) -> str:
        pass

    def intervals(self, params_index=None):
        pass

    def alpha_from_level(self, confidence_level):
        return round((1-confidence_level)/2, 3)


def confidence_estimates(
    truth, 
    data, 
    mle,
    cov_provider,
    cr_providers, 
    confidence_levels=[0.90, 0.95, 0.99]
):
    print(f"truth={truth}")
    fisher_mat = cov_provider.reg_fisher_mat()
    print(f"fisher_mat={fisher_mat}")
    cov_mat = cov_provider.cov_mat()
    sigmas = cov_provider.standard_errors()

    estimate_quality = []
    for k,v in mle.items():
        pindex = v["pindex"]
        pval = v["val"]
        if pindex < fisher_mat.shape[0]: # only include the unconstrained parameters
            estimate_quality.append({
                "pname": k,
                "index": pindex,
                "truth": truth[pindex],
                "estimate": pval,
                "fisher_val": fisher_mat[pindex, pindex],
                "variance": cov_mat[pindex, pindex],
                "stderr": sigmas[pindex]
            })

    cr_estimate = []
    for confidence_level in confidence_levels:
        for cr_provider in cr_providers:
            intervals = cr_provider.intervals(confidence_level=confidence_level)
            for k,v in mle.items():
                pindex = v["pindex"]
                pval = v["val"]
                if pindex < fisher_mat.shape[0]: # only include the unconstrained parameters
                    cr_estimate.append({
                        "pname": k,
                        "index": pindex,
                        "truth": truth[pindex],
                        "estimate": pval,
                        "method": cr_provider.type,
                        "level": confidence_level,
                        "lb": intervals[pindex][0],
                        "ub": intervals[pindex][1]
                    })
        

    return pd.DataFrame(estimate_quality), pd.DataFrame(cr_estimate)

