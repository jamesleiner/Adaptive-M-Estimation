# Code to replicate experiments for "Adaptive Off-Policy Inference for M-Estimators Under Model Misspecification"


This repository contains code to replicate all figures and experiments found in this paper

## Instructions
* First, run the script [`load_data.R`](load_data.R)  in order to compile the R dataset needed for semi-synthetic experiments. Note that this script compiles SAS datasets available for download from the Osteoarthritis Initiative
* The experiments are designed to be run on a cluster and then compiled. 
  * The script [`run_simulation.R`](run_simulation.R) takes a single argument from 1-3000. Each number corresponds to 30 different variations of hyperparameters replicated over 100 trials each.
  * We used a SLURM-managed batch computing environemnt to submit these jobs. See [`run_sims.sbatch`](run_sims.sbatch) for the exact command used to submit batch jobs.
  * After compilation, the script [`make_figs.R`](make_figs.R) compiles the figures available in the paper. 
