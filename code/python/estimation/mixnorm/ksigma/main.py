import typer
from typing import List
from typing_extensions import Annotated

app = typer.Typer(help="Estimates RLCT using toru's estimator, see ...", rich_help_panel="User Actions")
    
@app.command()
def rlctx(
    n: Annotated[int, 
        typer.Argument(
            help="Number of observations to use, must be pre-generated"
        )],

    targetdir: Annotated[str, 
        typer.Argument(
            help="Output directory for writting output data"
        )],

    # trials: Annotated[int, 
    #     typer.Option(
    #         help="Number of trials used, i.e. 10 means we are running 10 trials at the same temperature"
    #     )]=10,

    invtemp_scaling_factor: Annotated[int, 
        typer.Option(
            help="Inverse temperature scaling factor"
        )] = 1,
    
    mixture_components: Annotated[int, 
        typer.Option(
            help="Number of normal mixture components to use"
        )] = 3,
    
    mean_prior_cov_scaling: Annotated[float, 
        typer.Option(
            help="The variance of the prior of the mean components"
        )] = 4,

    pymc_chains: Annotated[int, 
        typer.Option(
            help="Number of parallel chains to generate"
        )] = 4,

    pymc_draws: Annotated[int, 
        typer.Option(
            help="Number of draws per chain to generate"
        )] = 10000,

    pymc_tune: Annotated[int, 
        typer.Option(
            help="Number of chain for HMC tunning"
        )] = 4000,

    pymc_cores: Annotated[int, 
        typer.Option(
            help="Number of cores for pymc to use"
        )] = 4,

    pymc_progressbar: Annotated[bool, 
        typer.Option(
            help="Whether or not to show pymc progress bar, for cluster run disable this otherwise the logs will be spammed"
        )] = True
):
    import numpy as np
    from pathlib import Path
    import pymc as pm
    import pytensor as pt
    import arviz as az

    from estimation.mixnorm.ksigma.model import tempered_normal_mixture

    pathstr = "estimation/mixnorm/ksigma/data/observations/n{n}.csv"
    x_data = np.loadtxt(Path(pathstr.format(n=n)), delimiter=',')
    print(f"Loaded observations for n={n}")

    model = tempered_normal_mixture(
        beta=invtemp_scaling_factor/np.log(len(x_data)), 
        data=x_data,
        n_components=mixture_components,
        weights_prior_alpha=np.full(mixture_components, 0.1),
        mean_prior_cov=pt.tensor.eye(mixture_components)*mean_prior_cov_scaling
        )
    print(f"Running on PyMC v{pm.__version__}")
    print(f"Created a tempered {mixture_components} normal mixture with {invtemp_scaling_factor} inverse temperature scale")
    print(model.str_repr())
    idata = None
    with model:
        idata = pm.sample(draws=pymc_draws,
                          tune=pymc_tune, 
                          chains=pymc_chains,
                          cores=pymc_cores,
                          max_treedepth=50,
                          target_accept=.995,
                          progressbar=pymc_progressbar)

    print(az.summary(idata, var_names=["weights", "mus"], round_to=2))
    # because of how az.extract does not extract what we want, 
    # let's extract every variable, flatten the column index and then join on draw and chain
    like_df = idata.posterior["like"].to_dataframe().reset_index()
    weights_df = idata.posterior["weights"].to_dataframe().unstack(level="weights_dim_0").reset_index()
    # flatten the column index
    weights_df.columns = weights_df.columns.map(lambda x: f"{x[0]}_{x[1]}" if isinstance(x, tuple) else x)
    weights_df.rename(columns={'chain_': 'chain', 'draw_': 'draw'}, inplace=True)
    
    mus_df = idata.posterior["mus"].to_dataframe().unstack(level="mus_dim_0").reset_index()
    # flatten the column index
    mus_df.columns = mus_df.columns.map(lambda x: f"{x[0]}_{x[1]}" if isinstance(x, tuple) else x)
    mus_df.rename(columns={'chain_': 'chain', 'draw_': 'draw'}, inplace=True)

    results = like_df.copy()
    results = results.merge(weights_df, on=['chain', 'draw'], how='inner')
    flat_results = results.merge(mus_df, on=['chain', 'draw'], how='inner')

    outputfile = f"{targetdir}/posterior_samples_n{n}.csv"
    print(f"Saving {len(flat_results)} samples with shape {flat_results.shape} in {outputfile}.")
    flat_results.to_csv(outputfile, index=False)
    print("We are done here!")


