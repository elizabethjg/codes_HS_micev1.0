#!/bin/bash
#SBATCH --job-name=LM_Hz10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --time=3-00:00:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar


. /etc/profile
source $HOME/.bashrc

##cd /mnt/projects/lensing/HALO_SHAPE/MICEv1.0/codes_HS_micev1.0/

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

conda activate myenv
##srun python -u forGroup_qprofile.py -sample 'LM_Hz' -lM_min 13.5 -lM_max 14.0 -z_min 0.3 -z_max 0.4 -ncores 56  -miscen True -soff 0.10 
##srun python -u forGroup_qprofile.py -sample 'HM_Hz' -lM_min 14.0 -lM_max 14.5 -z_min 0.3 -z_max 0.4 -ncores 56  -misalign 10
##srun python -u forGroup_qprofile.py -sample 'LM_Hz' -lM_min 14.0 -lM_max 14.5 -z_min 0.3 -z_max 0.4 -ncores 56  -misalign 30
srun python -u forGroup_qprofile.py -sample 'LM_Hz' -lM_min 13.5 -lM_max 14.0 -z_min 0.3 -z_max 0.4 -ncores 56  -misalign 20



##srun python -u forGroup_qprofile.py -sample 'HM_Lz_Nmis' -lM_min 14.0 -lM_max 14.5 -z_min 0.1 -z_max 0.2  -ncores 56 -miscen True -soff 0.14 
##srun python -u forGroup_qprofile.py -sample 'LM_Lz' -lM_min 13.5 -lM_max 14.0 -z_min 0.1 -z_max 0.2  -ncores 56 -miscen True -soff 0.10 -misalign 20
##srun python -u forGroup_qprofile.py -sample 'HM_Mz' -lM_min 14.0 -lM_max 14.5 -z_min 0.2 -z_max 0.3  -ncores 56 -miscen True -soff 0.14 -misalign 20
##srun python -u forGroup_qprofile.py -sample 'LM_Mz' -lM_min 13.5 -lM_max 14.0 -z_min 0.2 -z_max 0.3  -ncores 56 -miscen True -soff 0.10 -misalign 20
##srun python -u forGroup_qprofile.py -sample 'HM_Hz' -lM_min 14.0 -lM_max 14.5 -z_min 0.3 -z_max 0.4 -ncores 56  -miscen True -soff 0.14 -misalign 20
##srun python -u forGroup_qprofile.py -sample 'LM_Hz' -lM_min 13.5 -lM_max 14.0 -z_min 0.3 -z_max 0.4 -ncores 56  -miscen True -soff 0.10 -misalign 20
