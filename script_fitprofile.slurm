#!/bin/bash
#SBATCH --job-name=onlyt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --time=2-00:00:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

. /etc/profile
source $HOME/.bashrc

##cd /home/elizabeth/MICE/HS-lensing/codes_HS_micev1.0/

conda activate py3env


python fit_allprofiles_centred_sequential_2h.py -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Lz_relaxed.fits -RIN 250 -ROUT 5000  -ncores 15 -nit 100 -components tangential
##srun ./run_fit.sh $SLURM_JOB_NAME


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
