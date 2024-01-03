#!/bin/bash
#SBATCH --mem=10G     # memory; default unit is megabytes
#SBATCH --time=00-05:00           # time (DD-HH:MM)
#SBATCH --output=out-%x.out
#SBATCH --mail-user=melina.ribaud@free.fr # Send email updates to you or someone else
#SBATCH --mail-type=ALL
#SBATCH --ntasks-per-node=1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
module load gcc/9.3.0 r/4.0.2              # Adjust version and add the gcc module used for installing packages.

Rscript simu_V_fast.R
