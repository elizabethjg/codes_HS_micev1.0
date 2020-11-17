#!/bin/bash
#SBATCH --job-name='bin_130'
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --time=6:00:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar
##SBATCH --lm_min=13.0

. /etc/profile
source $HOME/.bashrc

cd /mnt/clemente/lensing/HALO_SHAPE/MICE_v1.0/codes_HS_micev1.0/

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

conda activate py3env

srun python -u forGroup_qprofile.py -sample $SLURM_JOB_NAME -lM_min 13.0 -lM_max 13.2 -ncores $SLURM_CPUS_PER_TASK
