#!/bin/bash
#SBATCH --mail-user=kdruscit@uwaterloo.ca
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --job-name="LCT_estimation"
#SBATCH --partition=cpu_pr3
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err


pushd ~/singular-learning-theory/code/python
N=100
DRAWS=1000
CHAINS=2
CORES=10
TUNE=1000
echo "N=$N, Draws=$DRAWS, Chains=$CHAINS, Cores=$CORES, Tune=$TUNE."
echo "Running on: $(hostname)"
date --iso-8601=seconds
srun ./.venv/bin/python -m estimation.main mixnorm_ksigma rlctx $N .\
                --pymc-draws $DRAWS\
               	--pymc-chains $CHAINS\
                --pymc-cores $CORES\
                --pymc-tune $TUNE

echo "Complete: $(date --iso-8601=seconds)"
