#!/bin/bash
#SBATCH --job-name=fit_simul
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --time=3-00:00:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

. /etc/profile
source $HOME/.bashrc

##cd /home/elizabeth/MICE/HS-lensing/codes_HS_micev1.0/

conda activate myenv


##python fit_allprofiles_centred_sequential_2h.py -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Lz_relaxed.fits -RIN 250 -ROUT 5000  -ncores 15 -nit 100 
##srun ./run_fit_2h.sh $SLURM_JOB_NAME
##srun ./run_fit_2h_2q.sh $SLURM_JOB_NAME
srun ./run_fit_misal.sh $SLURM_JOB_NAME
##srun ./run_fit_micen.sh $SLURM_JOB_NAME


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
