import numpy as np
import arviz as az


def column_stack_vars(idata, var_names, group="posterior"):
  draws_xarray = az.extract(
      idata,
      var_names=var_names,
      group=group,
      combined=True
  )
  
  if len(var_names) == 1:
    return np.column_stack(draws_xarray.to_numpy())
  else:
    v_draws = [draws_xarray[var_name].to_numpy() for var_name in var_names]
    return np.column_stack(v_draws)


