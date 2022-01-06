#!/bin/bash

conda activate py3env

python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$SLURM_JOB_NAME.fits -RIN 250 -ROUT 2000 -ncores 7 &
python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$SLURM_JOB_NAME_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7 &
python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$SLURM_JOB_NAME.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced' &
python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$SLURM_JOB_NAME_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced' &
python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$SLURM_JOB_NAME.fits -RIN 250 -ROUT 2000 -ncores 7 &
python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$SLURM_JOB_NAME_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7 &
python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$SLURM_JOB_NAME.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced' &
python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$SLURM_JOB_NAME_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced' &
