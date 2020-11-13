#!/bin/bash
#SBATCH --job-name=pMice
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --time=6:00:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar

. /etc/profile
source $HOME/.bashrc

cd /mnt/clemente/lensing/HALO_SHAPE/MICE_v1.0/codes_HS_micev1.0/

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

conda activate py3env

srun python -u forGroup_profile.py -sample 'Mbin6_zH_cM'  -sample 'bin14' -ncores 56
