#!/bin/bash
#SBATCH --job-name=relax
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --time=2-00:00:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar


. /etc/profile
source $HOME/.bashrc

cd /mnt/projects/lensing/HALO_SHAPE/MICEv1.0/codes_HS_micev1.0/

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

conda activate py3env

srun python -u forGroup_qprofile.py -sample 'HM' -lM_min 14.4 -lM_max 15.0 -z_min 0.1 -z_max 0.2 -ncores 56 -relax 'True' -q_min 0.4
srun python -u forGroup_qprofile.py -sample 'MM' -lM_min 14.0 -lM_max 14.2 -z_min 0.1 -z_max 0.2 -ncores 56 -relax 'True' -q_min 0.4
srun python -u forGroup_qprofile.py -sample 'LM' -lM_min 13.5 -lM_max 13.6 -z_min 0.1 -z_max 0.2 -ncores 56 -relax 'True' -q_min 0.4

##srun python -u forGroup_qprofile.py -sample 'HM_Hz' -lM_min 14.0 -lM_max 14.5 -z_min 0.3 -z_max 0.35 -ncores 56 -domap 'True' -nbins 10000 -ROUT 5000
##srun python -u forGroup_qprofile.py -sample 'HM_Lz' -lM_min 14.0 -lM_max 14.5 -z_min 0.1 -z_max 0.2  -ncores 56 -domap 'True' -nbins 10000 -ROUT 5000
##srun python -u forGroup_qprofile.py -sample 'HM_Hz' -lM_min 14.0 -lM_max 14.5 -z_min 0.3 -z_max 0.35 -ncores 56 -relax 'True' -domap 'True' -nbins 10000 -ROUT 5000
##srun python -u forGroup_qprofile.py -sample 'HM_Lz' -lM_min 14.0 -lM_max 14.5 -z_min 0.1 -z_max 0.2  -ncores 56 -relax 'True' -domap 'True' -nbins 10000 -ROUT 5000
##
##srun python -u forGroup_qprofile.py -sample 'LM_Hz' -lM_min 13.5 -lM_max 14. -z_min 0.3 -z_max 0.35 -ncores 56 -domap 'True' -nbins 10000 -ROUT 5000
##srun python -u forGroup_qprofile.py -sample 'LM_Lz' -lM_min 13.5 -lM_max 14. -z_min 0.1 -z_max 0.2  -ncores 56 -domap 'True' -nbins 10000 -ROUT 5000
##srun python -u forGroup_qprofile.py -sample 'LM_Hz' -lM_min 13.5 -lM_max 14. -z_min 0.3 -z_max 0.35 -ncores 56 -relax 'True' -domap 'True' -nbins 10000 -ROUT 5000
##srun python -u forGroup_qprofile.py -sample 'LM_Lz' -lM_min 13.5 -lM_max 14. -z_min 0.1 -z_max 0.2  -ncores 56 -relax 'True' -domap 'True' -nbins 10000 -ROUT 5000
