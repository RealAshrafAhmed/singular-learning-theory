import typer
from typing import List
from typing_extensions import Annotated

app = typer.Typer(help="Estimates RLCT using toru's estimator, see ...", rich_help_panel="User Actions")
    
@app.command()
def rlctx(
    n: Annotated[int, 
        typer.Argument(
            help="A comma separated list of random sample sizes."
        )],

    targetdir: Annotated[str, 
        typer.Argument(
            help="Output directory for where we should write the output data"
        )],

    # trials: Annotated[int, 
    #     typer.Option(
    #         help="Number of trials used, i.e. 10 means we are running 10 trials at the same temperature"
    #     )]=10,

    invtemp_scaling_factor: Annotated[int, 
        typer.Option(
            help="A comma separated list of scaling factors. Default is 1"
        )] = 1,
    
    mixture_components: Annotated[int, 
        typer.Option(
            help="Number of normal mixture components."
        )] = 3,
    
    mean_prior_cov_scaling: Annotated[float, 
        typer.Option(
            help="the variance of the mvn for the weights raw prior"
        )] = 4,

    parallel_chains: Annotated[int, 
        typer.Option(
            help="the variance of the mvn for the weights raw prior"
        )] = 4,

    draws_per_chain: Annotated[int, 
        typer.Option(
            help="the variance of the mvn for the weights raw prior"
        )] = 10000,

    pmi_cores: Annotated[int, 
        typer.Option(
            help="the variance of the mvn for the weights raw prior"
        )] = 4,

    pymc_progressbar: Annotated[bool, 
        typer.Option(
            help="Show pymc progress bar"
        )] = False
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
        idata = pm.sample(draws=draws_per_chain,
                          tune=4000, 
                          chains=parallel_chains,
                          cores=pmi_cores,
                          max_treedepth=50,
                          target_accept=.995,
                          progressbar=pymc_progressbar)

    print(az.summary(idata, var_names=["weights", "mu"], round_to=2))
    results = az.extract(idata, group="posterior").to_dataframe()
    outputfile = f"{targetdir}/posterior_samples_n{n}.csv"
    print(f"Saving {len(results)} samples in {outputfile}.")
    results.to_csv(outputfile, index=False)
    print("We are done here!")


