#!/bin/bash
#SBATCH --job-name=allq
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
##srun python -u forGroup_qprofile.py -sample 'LM_Hz' -lM_min 13.5 -lM_max 14.0 -z_min 0.3 -z_max 0.4 -ncores 56  -misalign 30

##srun python -u forGroup_qprofile.py -sample 'all' -lM_min 13.5 -lM_max 14.5 -z_min 0.1 -z_max 0.3 -ncores 56 
srun python -u forGroup_qprofile.py -sample 'all_q1' -lM_min 13.5 -lM_max 14.5 -z_min 0.1 -z_max 0.3 -ncores 56 -q_min 0.0 -q_max 0.3 
srun python -u forGroup_qprofile.py -sample 'all_q2' -lM_min 13.5 -lM_max 14.5 -z_min 0.1 -z_max 0.3 -ncores 56 -q_min 0.3 -q_max 0.5 
srun python -u forGroup_qprofile.py -sample 'all_q3' -lM_min 13.5 -lM_max 14.5 -z_min 0.1 -z_max 0.3 -ncores 56 -q_min 0.5 -q_max 0.6 
srun python -u forGroup_qprofile.py -sample 'all_q4' -lM_min 13.5 -lM_max 14.5 -z_min 0.1 -z_max 0.3 -ncores 56 -q_min 0.6 -q_max 0.7 
srun python -u forGroup_qprofile.py -sample 'all_q5' -lM_min 13.5 -lM_max 14.5 -z_min 0.1 -z_max 0.3 -ncores 56 -q_min 0.7 -q_max 0.8 
srun python -u forGroup_qprofile.py -sample 'all_q6' -lM_min 13.5 -lM_max 14.5 -z_min 0.1 -z_max 0.3 -ncores 56 -q_min 0.8 -q_max 1.0 

##srun python -u forGroup_qprofile.py -sample 'HM_Lz' -lM_min 14.0 -lM_max 14.5 -z_min 0.1 -z_max 0.2 -ncores 56 -relax True -miscen True -soff 0.14 -misalign 20
##srun python -u forGroup_qprofile.py -sample 'LM_Lz' -lM_min 13.5 -lM_max 14.0 -z_min 0.1 -z_max 0.2 -ncores 56 -relax True -miscen True -soff 0.10 -misalign 20
##srun python -u forGroup_qprofile.py -sample 'LM_Mz' -lM_min 13.5 -lM_max 14.0 -z_min 0.2 -z_max 0.3 -ncores 56 -relax True -miscen True -soff 0.10 -misalign 20
##srun python -u forGroup_qprofile.py -sample 'HM_Mz' -lM_min 14.0 -lM_max 14.5 -z_min 0.2 -z_max 0.3 -ncores 56 -relax True -miscen True -soff 0.14 -misalign 20
##srun python -u forGroup_qprofile.py -sample 'HM_Hz' -lM_min 14.0 -lM_max 14.5 -z_min 0.3 -z_max 0.4 -ncores 56 -relax True -miscen True -soff 0.14 -misalign 20
##srun python -u forGroup_qprofile.py -sample 'LM_Hz' -lM_min 13.5 -lM_max 14.0 -z_min 0.3 -z_max 0.4 -ncores 56 -relax True -miscen True -soff 0.10 -misalign 20
