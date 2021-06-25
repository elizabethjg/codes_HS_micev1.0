#!/bin/bash
#SBATCH --job-name=new_rel
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=55
#SBATCH --time=6:00:00
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

srun python -u forGroup_qprofile.py -sample 'test_v1_z15' -lM_min 13.5 -z_max 0.5 -ncores 55  -vmice 1
srun python -u forGroup_qprofile.py -sample 'test_v1_z56' -lM_min 13.5 -z_min 0.5 -z_max 0.6  -ncores 55  -vmice 1
srun python -u forGroup_qprofile.py -sample 'test_v2_z15' -lM_min 13.5 -z_max 0.5 -ncores 55  -vmice 2
srun python -u forGroup_qprofile.py -sample 'test_v2_z56' -lM_min 13.5 -z_min 0.5 -z_max 0.6  -ncores 55  -vmice 2

##srun python -u forGroup_qprofile.py -sample 'test_newparam'    -lM_min 13.5 -z_max 2 -ncores 56 -new_version 0

##srun python -u compute_miss.py -sample pru_iso -ncores $SLURM_CPUS_PER_TASK
##srun python -u compute_miss.py -sample pru_dens -ncores $SLURM_CPUS_PER_TASK
##srun python -u compute_miss.py -sample pru_elip -ncores $SLURM_CPUS_PER_TASK
##srun python -u compute_miss.py -sample pru_round -ncores $SLURM_CPUS_PER_TASK

##srun python -u forGroup_qprofile_individual.py -sample '136_145_zc' -lM_min 13.9 -lM_max 14.5 -z_max 0.12 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample '140_145_sel4' -lM_min 14.0 -lM_max 14.5 -z_max 0.25 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_map.py -sample '140_145_sel4' -lM_min 14.0 -lM_max 14.5 -z_max 0.25 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_map.py -sample '136_145_nc' -lM_min 13.6 -lM_max 14.5 -z_max 0.12 -ncores $SLURM_CPUS_PER_TASK

##srun python -u forGroup_map.py -sample 'bin_134' -lM_min 13.4 -lM_max 13.8 -ncores $SLURM_CPUS_PER_TASK

##srun python -u fit_profile_fullmodel_allmis.py -file 'profile_ebin_142.fits' -ncores $SLURM_CPUS_PER_TASK

##srun python -u compute_miss.py -sample $SLURM_JOB_NAME -ncores $SLURM_CPUS_PER_TASK
##srun python -u compute_miss.py -sample pru_round -ncores $SLURM_CPUS_PER_TASK
##srun python -u compute_miss.py -sample pru_elip -ncores $SLURM_CPUS_PER_TASK
##srun python -u compute_miss.py -sample pru_dens -ncores $SLURM_CPUS_PER_TASK
##srun python -u compute_miss.py -sample pru_iso -ncores $SLURM_CPUS_PER_TASK

##srun python -u forGroup_qprofile.py -sample 'pru_hz' -lM_min 13.4 -lM_max 14. -z_min 0.30 -z_max 0.31 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample 'pru_round_hz' -q_min 0.7 -lM_min 13.4 -lM_max 14. -z_min 0.30 -z_max 0.31 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_map.py -sample 'pru_round_hz' -q_min 0.7 -lM_min 13.4 -lM_max 14. -z_min 0.30 -z_max 0.31 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample 'pru_elip_hz' -q_max 0.5 -lM_min 13.4 -lM_max 14. -z_min 0.30 -z_max 0.31 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_map.py -sample 'pru_elip_hz' -q_max 0.5 -lM_min 13.4 -lM_max 14. -z_min 0.30 -z_max 0.31 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample 'pru_dens_hz' -rmax 30. -lM_min 13.4 -lM_max 14. -z_min 0.30 -z_max 0.31 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_map.py -sample 'pru_dens_hz' -rmax 30. -lM_min 13.4 -lM_max 14. -z_min 0.30 -z_max 0.31 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample 'pru_iso_hz' -rmin 30. -lM_min 13.4 -lM_max 14. -z_min 0.30 -z_max 0.31 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_map.py -sample 'pru_iso_hz'-rmin 30. -lM_min 13.4 -lM_max 14. -z_min 0.30 -z_max 0.31 -ncores $SLURM_CPUS_PER_TASK

##srun python -u forGroup_map.py -sample 'bin_144' -lM_min 14.4 -lM_max 14.6 -ncores 56
##srun python -u forGroup_map.py -sample 'bin_142' -lM_min 14.2 -lM_max 14.4 -ncores 56
##srun python -u forGroup_map.py -sample 'bin_140' -lM_min 14.0 -lM_max 14.2 -ncores 56
##srun python -u forGroup_map.py -sample 'bin_138' -lM_min 13.8 -lM_max 14.0 -ncores 56
##srun python -u forGroup_map.py -sample 'bin_136' -lM_min 13.6 -lM_max 13.8 -ncores 56

##srun python -u forGroup_qprofile.py -sample 'ebin_130' -rmin 0. -lM_min 13.0 -lM_max 13.2 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample 'ebin_132' -rmin 0. -lM_min 13.2 -lM_max 13.4 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample 'ebin_134' -rmin 0. -lM_min 13.4 -lM_max 13.6 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample 'ebin_136' -rmin 0. -lM_min 13.6 -lM_max 13.8 -ncores 56
##srun python -u forGroup_qprofile.py -sample 'ebin_138' -rmin 0. -lM_min 13.8 -lM_max 14.0 -ncores 56
##srun python -u forGroup_qprofile.py -sample 'ebin_140' -rmin 0. -lM_min 14.0 -lM_max 14.2 -ncores 56
##srun python -u forGroup_qprofile.py -sample 'ebin_142' -rmin 0. -lM_min 14.2 -lM_max 14.4 -ncores 56
##srun python -u forGroup_qprofile.py -sample 'ebin_144' -rmin 0. -lM_min 14.4 -lM_max 16.0 -ncores 56
##
##srun python -u forGroup_qprofile.py -sample 'ebin_138_dens' -rmax 30. -lM_min 13.8 -lM_max 14.0 -ncores 56
##srun python -u forGroup_qprofile.py -sample 'ebin_140_dens' -rmax 30. -lM_min 14.0 -lM_max 14.2 -ncores 56
##srun python -u forGroup_qprofile.py -sample 'ebin_142_dens' -rmax 30. -lM_min 14.2 -lM_max 14.4 -ncores 56
##srun python -u forGroup_qprofile.py -sample 'ebin_144_dens' -rmax 30. -lM_min 14.4 -lM_max 16.0 -ncores 56


##srun python -u forGroup_qprofile.py -sample 'bin_130_iso' -rmin 30. -lM_min 13.0 -lM_max 13.2 -ncores $SLURM_CPUS_PER_TASK

##srun python -u forGroup_qprofile.py -sample 'bin_132' -lM_min 13.2 -lM_max 13.4 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample 'bin_134' -lM_min 13.4 -lM_max 13.6 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample 'bin_136' -lM_min 13.6 -lM_max 13.8 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample 'bin_138' -lM_min 13.8 -lM_max 14.0 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample 'bin_140' -lM_min 14.0 -lM_max 14.2 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample 'bin_142' -lM_min 14.2 -lM_max 14.4 -ncores $SLURM_CPUS_PER_TASK
##srun python -u forGroup_qprofile.py -sample 'bin_144' -lM_min 14.4 -lM_max 15.0 -ncores $SLURM_CPUS_PER_TASK


