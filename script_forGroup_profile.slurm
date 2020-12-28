#!/bin/bash
#SBATCH --job-name='bin_138'
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --time=2-00:00
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

##srun python -u forGroup_qprofile.py -sample 'bin_130' -lM_min 13.0 -lM_max 13.2 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample 'bin_132' -lM_min 13.2 -lM_max 13.4 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample 'bin_134' -lM_min 13.4 -lM_max 13.6 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample 'bin_136' -lM_min 13.6 -lM_max 13.8 -ncores $SLURM_CPUS_PER_TASK
srun python -u forGroup_qprofile.py -sample 'bin_138' -lM_min 13.8 -lM_max 14.0 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample 'bin_140' -lM_min 14.0 -lM_max 14.2 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample 'bin_144' -lM_min 14.4 -lM_max 15.0 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample 'bin_146' -lM_min 14.6 -lM_max 15.0 -ncores $SLURM_CPUS_PER_TASK

