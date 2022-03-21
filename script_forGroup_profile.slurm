#!/bin/bash
#SBATCH --job-name=addsamps
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

srun python -u forGroup_qprofile.py -sample 'HM_Mz' -lM_min 14.0 -lM_max 14.5 -z_min 0.2 -z_max 0.3  -ncores 56 -relax 'True' 
srun python -u forGroup_qprofile.py -sample 'LM_Mz' -lM_min 13.5 -lM_max 14.0 -z_min 0.2 -z_max 0.3  -ncores 56 -relax 'True' 
srun python -u forGroup_qprofile.py -sample 'LM_HHz' -lM_min 13.5 -lM_max 14.0 -z_min 0.35 -z_max 0.4 -ncores 56 -relax 'True'
srun python -u forGroup_qprofile.py -sample 'HM_HHz' -lM_min 14.0 -lM_max 14.5 -z_min 0.35 -z_max 0.4 -ncores 56 -relax 'True'
srun ./run_fit_forz_400.sh HHz_relaxed
srun ./run_fit_forz_400.sh Mz_relaxed
srun ./run_fit_forz_250.sh HHz_relaxed
srun ./run_fit_forz_250.sh Mz_relaxed

##srun python -u forGroup_qprofile.py -sample 'HM_Lz_soff' -lM_min 14.0 -lM_max 14.5 -z_min 0.1 -z_max 0.2 -ncores 56 -relax 'True' -miscen 'True' -soff 0.3  -misaling 20
##srun python -u forGroup_qprofile.py -sample 'HM_Lz_round_obl' -lM_min 14.0 -lM_max 14.5 -z_min 0.1 -z_max 0.2  -ncores 56 -relax 'True' -q_min 0.6 -T_max 0.8
##srun python -u forGroup_qprofile.py -sample 'LM_Lz_round_obl' -lM_min 13.5 -lM_max 14.0 -z_min 0.1 -z_max 0.2  -ncores 56 -relax 'True' -q_min 0.6 -T_max 0.8
##srun python -u forGroup_qprofile.py -sample 'LM_Hz_round_obl' -lM_min 13.5 -lM_max 14.0 -z_min 0.3 -z_max 0.35 -ncores 56 -relax 'True' -q_min 0.6 -T_max 0.8
##srun python -u forGroup_qprofile.py -sample 'HM_Hz_round_obl' -lM_min 14.0 -lM_max 14.5 -z_min 0.3 -z_max 0.35 -ncores 56 -relax 'True' -q_min 0.6 -T_max 0.8

##srun python -u forGroup_qprofile.py -sample 'HM_Hz' -lM_min 14.0 -lM_max 14.5 -z_min 0.3 -z_max 0.35 -ncores 56 -domap 'True' -nbins 10000 -ROUT 5000
##srun python -u forGroup_qprofile.py -sample 'HM_Lz' -lM_min 14.0 -lM_max 14.5 -z_min 0.1 -z_max 0.2  -ncores 56 -domap 'True' -nbins 10000 -ROUT 5000
##srun python -u forGroup_qprofile.py -sample 'HM_Hz' -lM_min 14.0 -lM_max 14.5 -z_min 0.3 -z_max 0.35 -ncores 56 -relax 'True' -domap 'True' -nbins 10000 -ROUT 5000
##srun python -u forGroup_qprofile.py -sample 'HM_Lz' -lM_min 14.0 -lM_max 14.5 -z_min 0.1 -z_max 0.2  -ncores 56 -relax 'True' -domap 'True' -nbins 10000 -ROUT 5000

##srun python -u forGroup_qprofile.py -sample 'LM_Hz' -lM_min 13.5 -lM_max 14. -z_min 0.3 -z_max 0.35 -ncores 56 -domap 'True' -nbins 10000 -ROUT 5000
##srun python -u forGroup_qprofile.py -sample 'LM_Lz' -lM_min 13.5 -lM_max 14. -z_min 0.1 -z_max 0.2  -ncores 56 -domap 'True' -nbins 10000 -ROUT 5000
##srun python -u forGroup_qprofile.py -sample 'LM_Hz' -lM_min 13.5 -lM_max 14. -z_min 0.3 -z_max 0.35 -ncores 56 -relax 'True' -domap 'True' -nbins 10000 -ROUT 5000
##srun python -u forGroup_qprofile.py -sample 'LM_Lz' -lM_min 13.5 -lM_max 14. -z_min 0.1 -z_max 0.2  -ncores 56 -relax 'True' -domap 'True' -nbins 10000 -ROUT 5000
