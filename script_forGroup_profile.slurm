#!/bin/bash
#SBATCH --job-name=Mbin_zcM
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --time=6:00:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar

. /etc/profile
source $HOME/.bashrc

cd /mnt/clemente/lensing/RodriguezGroups/codes_RGroups/

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

conda activate py2env

##srun python -u forGroup_profile.py -sample 'Mbin1_cM'  -lMH_min 12.5 -lMH_max 13.0 -C_BG_min 3. -z_max 0.2 -ncores 56
##srun python -u forGroup_profile.py -sample 'Mbin2_cM'  -lMH_min 13.0 -lMH_max 13.3 -C_BG_min 3. -z_max 0.2 -ncores 56
##srun python -u forGroup_profile.py -sample 'Mbin3_cM'  -lMH_min 13.3 -lMH_max 13.6 -C_BG_min 3. -z_max 0.2 -ncores 56
##srun python -u forGroup_profile.py -sample 'Mbin4_cM'  -lMH_min 13.6 -lMH_max 13.9 -C_BG_min 3. -z_max 0.2 -ncores 56
##srun python -u forGroup_profile.py -sample 'Mbin5_cM'  -lMH_min 13.9 -lMH_max 14.2 -C_BG_min 3. -z_max 0.2 -ncores 56
##srun python -u forGroup_profile.py -sample 'Mbin6_cM'  -lMH_min 14.2 -lMH_max 15.5 -C_BG_min 3. -z_max 0.2 -ncores 56

srun python -u forGroup_profile.py -sample 'Mbin1_zL_cM'  -lMH_min 12.5 -lMH_max 13.0 -C_BG_min 3. -z_max 0.13 -ncores 45
srun python -u forGroup_profile.py -sample 'Mbin2_zL_cM'  -lMH_min 13.0 -lMH_max 13.3 -C_BG_min 3. -z_max 0.13 -ncores 45
srun python -u forGroup_profile.py -sample 'Mbin3_zL_cM'  -lMH_min 13.3 -lMH_max 13.6 -C_BG_min 3. -z_max 0.13 -ncores 45
srun python -u forGroup_profile.py -sample 'Mbin4_zL_cM'  -lMH_min 13.6 -lMH_max 13.9 -C_BG_min 3. -z_max 0.13 -ncores 45
srun python -u forGroup_profile.py -sample 'Mbin5_zL_cM'  -lMH_min 13.9 -lMH_max 14.2 -C_BG_min 3. -z_max 0.13 -ncores 45
srun python -u forGroup_profile.py -sample 'Mbin6_zL_cM'  -lMH_min 14.2 -lMH_max 15.5 -C_BG_min 3. -z_max 0.13 -ncores 45
                                                                                                
srun python -u forGroup_profile.py -sample 'Mbin1_zH_cM'  -lMH_min 12.5 -lMH_max 13.0 -C_BG_min 3. -z_min 0.13 -z_max 0.2 -ncores 45
srun python -u forGroup_profile.py -sample 'Mbin2_zH_cM'  -lMH_min 13.0 -lMH_max 13.3 -C_BG_min 3. -z_min 0.13 -z_max 0.2 -ncores 45
srun python -u forGroup_profile.py -sample 'Mbin3_zH_cM'  -lMH_min 13.3 -lMH_max 13.6 -C_BG_min 3. -z_min 0.13 -z_max 0.2 -ncores 45
srun python -u forGroup_profile.py -sample 'Mbin4_zH_cM'  -lMH_min 13.6 -lMH_max 13.9 -C_BG_min 3. -z_min 0.13 -z_max 0.2 -ncores 45
srun python -u forGroup_profile.py -sample 'Mbin5_zH_cM'  -lMH_min 13.9 -lMH_max 14.2 -C_BG_min 3. -z_min 0.13 -z_max 0.2 -ncores 45
srun python -u forGroup_profile.py -sample 'Mbin6_zH_cM'  -lMH_min 14.2 -lMH_max 15.5 -C_BG_min 3. -z_min 0.13 -z_max 0.2 -ncores 45

 python -u forGroup_profile.py -sample 'Mbin6b_cM'     -lMH_min 14.2 -lMH_max 14.5 -C_BG_min 3. -z_max 0.2 -ncores 48
 python -u forGroup_profile.py -sample 'Mbin6b_zH'     -lMH_min 14.2 -lMH_max 14.5 -C_BG_min 0. -z_min 0.13 -z_max 0.2 -ncores 48
 python -u forGroup_profile.py -sample 'Mbin6b_zH_cM'  -lMH_min 14.2 -lMH_max 14.5 -C_BG_min 3. -z_min 0.13 -z_max 0.2 -ncores 48
 python -u forGroup_profile.py -sample 'Mbin6b_zL'     -lMH_min 14.2 -lMH_max 14.5 -C_BG_min 0. -z_max 0.13 -ncores 48
 python -u forGroup_profile.py -sample 'Mbin6b_zL_cM'  -lMH_min 14.2 -lMH_max 14.5 -C_BG_min 3. -z_max 0.13 -ncores 48
 python -u forGroup_profile.py -sample 'Mbin7b_zL'     -lMH_min 14.5 -lMH_max 15. -C_BG_min 0. -z_max 0.13 -ncores 48
 python -u forGroup_profile.py -sample 'Mbin7b_zL_cM'  -lMH_min 14.5 -lMH_max 15. -C_BG_min 3. -z_max 0.13 -ncores 48
 python -u forGroup_profile.py -sample 'Mbin7b_zH'     -lMH_min 14.5 -lMH_max 15. -C_BG_min 0. -z_min 0.13 -z_max 0.2 -ncores 48
 python -u forGroup_profile.py -sample 'Mbin7b_zH_cM'  -lMH_min 14.5 -lMH_max 15. -C_BG_min 3. -z_min 0.13 -z_max 0.2 -ncores 48
 
