#!/bin/bash
#SBATCH --mail-user=kdruscit@uwaterloo.ca
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --job-name="LCT_estimation"
#SBATCH --partition=cpu_pr3
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --output=./slurm_output/%x-%j.out
#SBATCH --error=./slurm_output/%x-%j.err

srun ./.venv/bin/python -m estimation.main mixnorm_ksigma rlctx 1000 .\
                --pymc-draws 5000\
               	--pymc-chains 16\
                --pymc-cores 256\
                --pymc-tune 1000
