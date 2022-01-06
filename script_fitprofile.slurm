#!/bin/bash
#SBATCH --job-name=HM_Lz
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

python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$SLURM_JOB_NAME.fits -RIN 250 -ROUT 2000 -ncores 7
python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$SLURM_JOB_NAME_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7
python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$SLURM_JOB_NAME.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced'
python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$SLURM_JOB_NAME_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced'
python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$SLURM_JOB_NAME.fits -RIN 250 -ROUT 2000 -ncores 7
python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$SLURM_JOB_NAME_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7
python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$SLURM_JOB_NAME.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced'
python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$SLURM_JOB_NAME_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced'

##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_HM_Lz.fits -RIN 250 -ROUT 2000 -ncores 56
##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_HM_Hz.fits -RIN 250 -ROUT 2000 -ncores 56
##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_LM_Lz.fits -RIN 250 -ROUT 2000 -ncores 56
##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_LM_Hz.fits -RIN 250 -ROUT 2000 -ncores 56
####
##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_HM_Lz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 56
##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_HM_Hz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 56
##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_LM_Lz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 56
##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_LM_Hz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 56
####
##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_HM_Lz.fits -RIN 250 -ROUT 2000 -ncores 56 -ang 'reduced'
##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_HM_Hz.fits -RIN 250 -ROUT 2000 -ncores 56 -ang 'reduced'
##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_LM_Lz.fits -RIN 250 -ROUT 2000 -ncores 56 -ang 'reduced'
##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_LM_Hz.fits -RIN 250 -ROUT 2000 -ncores 56 -ang 'reduced'
####
##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_HM_Lz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 56 -ang 'reduced'
##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_HM_Hz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 56 -ang 'reduced'
##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_LM_Lz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 56 -ang 'reduced'
##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_LM_Hz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 56 -ang 'reduced'
####
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_HM_Lz.fits -RIN 250 -ROUT 2000 -ncores 56
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_HM_Hz.fits -RIN 250 -ROUT 2000 -ncores 56
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_LM_Lz.fits -RIN 250 -ROUT 2000 -ncores 56
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_LM_Hz.fits -RIN 250 -ROUT 2000 -ncores 56
####
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_HM_Lz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 56
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_HM_Hz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 56
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_LM_Lz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 56
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_LM_Hz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 56
####
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_HM_Lz.fits -RIN 250 -ROUT 2000 -ncores 56 -ang 'reduced'
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_HM_Hz.fits -RIN 250 -ROUT 2000 -ncores 56 -ang 'reduced'
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_LM_Lz.fits -RIN 250 -ROUT 2000 -ncores 56 -ang 'reduced'
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_LM_Hz.fits -RIN 250 -ROUT 2000 -ncores 56 -ang 'reduced'
####
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_HM_Lz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 56 -ang 'reduced'
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_HM_Hz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 56 -ang 'reduced'
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_LM_Lz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 56 -ang 'reduced'
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_LM_Hz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 56 -ang 'reduced'
