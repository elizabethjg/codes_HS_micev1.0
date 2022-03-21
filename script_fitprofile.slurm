#!/bin/bash
#SBATCH --job-name=fnonrex
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --time=3-00:00:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

. /etc/profile
source $HOME/.bashrc

##cd /home/elizabeth/MICE/HS-lensing/codes_HS_micev1.0/

conda activate myenv

srun ./run_fit_2h_2q
srun ./run_fit_forz_400.sh Hz_relaxed
srun ./run_fit_forz_400.sh Lz_relaxed
srun ./run_fit_forz_400.sh Lz
srun ./run_fit_forz_250.sh Lz
srun ./run_fit_forz_400.sh Hz
srun ./run_fit_forz_250.sh Hz
