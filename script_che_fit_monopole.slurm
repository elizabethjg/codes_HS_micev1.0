#!/bin/bash
#SBATCH --job-name=monomis134
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=3-00:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar


export OMP_NUM_THREADS=10

module load softwares/anaconda3/5.0-intel-2018.0
source activate py3env

srun python -u fit_profile_monopole.py -folder '/home/elizabeth/MICEv2.0/profiles/' -file 'profile_ebin_134.fits' -ncores $SLURM_CPUS_PER_TASK

##srun -u python fit_profile_monopole.py -folder '/home/elizabeth/MICEv2.0/profiles/' -file 'profile_ebin_142.fits' -ROUT 2500. -ncores 10 
##srun -u python compute_miss.py


##srun -u python fit_profile_monopole_misscentred_pcc.py -folder '/home/elizabeth/profiles_indcat/' \
##-file 'profile_'$SLURM_JOB_NAME'.fits' -ncores 15                                    
##srun -u python fit_profile_monopole_misscentred_pcc.py -folder '/home/elizabeth/profiles_indcat/' \
##-file 'profile_'$SLURM_JOB_NAME'.fits' -Yanmiss 'True' -ncores 15                                    
##srun -u python fit_profile_monopole_misscentred_pcc.py -folder '/home/elizabeth/profiles_new/' \
##-file 'profile_'$SLURM_JOB_NAME'_cM.fits' -ncores 15                                    
##srun -u python fit_profile_monopole_misscentred_pcc.py -folder '/home/elizabeth/profiles_newanalysis/' \
##-file 'profile_'$SLURM_JOB_NAME'.fits' -Yanmiss 'True' -ncores 15                                    
##srun -u python fit_profile_monopole_misscentred_pcc.py -folder '/home/elizabeth/profiles_new/' \
##-file 'profile_'$SLURM_JOB_NAME'_cM.fits' -Yanmiss 'True' -ncores 15                                    
