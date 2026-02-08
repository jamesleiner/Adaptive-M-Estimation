#!/bin/bash
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 36:00:00
#SBATCH --mail-type=ALL
#SBATCH --array=1-1000
VAR=$((SLURM_ARRAY_TASK_ID))
echo $VAR

module load anaconda3
Rscript run_simulation.R $VAR