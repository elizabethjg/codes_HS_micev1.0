#!/bin/bash
#SBATCH --job-name='all2'
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --time=3-00:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar
##SBATCH --lm_min=13.0

. /etc/profile
source $HOME/.bashrc

cd /mnt/clemente/lensing/HALO_SHAPE/MICE_v1.0/codes_HS_micev1.0/

export OMP_NUM_THREADS=56
export MKL_NUM_THREADS=56

conda activate py3env

##srun python -u forGroup_qprofile.py -sample 'ebin_130' -rmin 0. -lM_min 13.0 -lM_max 13.2 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample 'ebin_132' -rmin 0. -lM_min 13.2 -lM_max 13.4 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample 'ebin_134' -rmin 0. -lM_min 13.4 -lM_max 13.6 -ncores $SLURM_CPUS_PER_TASK
srun python -u forGroup_qprofile.py -sample 'ebin_136' -rmin 0. -lM_min 13.6 -lM_max 13.8 -ncores 56
srun python -u forGroup_qprofile.py -sample 'ebin_138' -rmin 0. -lM_min 13.8 -lM_max 14.0 -ncores 56
srun python -u forGroup_qprofile.py -sample 'ebin_140' -rmin 0. -lM_min 14.0 -lM_max 14.2 -ncores 56
srun python -u forGroup_qprofile.py -sample 'ebin_142' -rmin 0. -lM_min 14.2 -lM_max 14.4 -ncores 56
srun python -u forGroup_qprofile.py -sample 'ebin_144' -rmin 0. -lM_min 14.4 -lM_max 16.0 -ncores 56

srun python -u forGroup_qprofile.py -sample 'ebin_138_dens' -rmax 30. -lM_min 13.8 -lM_max 14.0 -ncores 56
srun python -u forGroup_qprofile.py -sample 'ebin_140_dens' -rmax 30. -lM_min 14.0 -lM_max 14.2 -ncores 56
srun python -u forGroup_qprofile.py -sample 'ebin_142_dens' -rmax 30. -lM_min 14.2 -lM_max 14.4 -ncores 56
srun python -u forGroup_qprofile.py -sample 'ebin_144_dens' -rmax 30. -lM_min 14.4 -lM_max 16.0 -ncores 56


##srun python -u forGroup_qprofile.py -sample 'bin_130_iso' -rmin 30. -lM_min 13.0 -lM_max 13.2 -ncores $SLURM_CPUS_PER_TASK

##srun python -u forGroup_qprofile.py -sample 'bin_132' -lM_min 13.2 -lM_max 13.4 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample 'bin_134' -lM_min 13.4 -lM_max 13.6 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample 'bin_136' -lM_min 13.6 -lM_max 13.8 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample 'bin_138' -lM_min 13.8 -lM_max 14.0 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample 'bin_140' -lM_min 14.0 -lM_max 14.2 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample 'bin_142' -lM_min 14.2 -lM_max 14.4 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample 'bin_144' -lM_min 14.4 -lM_max 15.0 -ncores $SLURM_CPUS_PER_TASK


