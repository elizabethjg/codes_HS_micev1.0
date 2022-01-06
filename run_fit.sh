#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate py3env

python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz.fits -RIN 250 -ROUT 2000 -ncores 7 &
python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7 &
python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced' &
python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced' &
python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz.fits -RIN 250 -ROUT 2000 -ncores 7 &
python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7 &
python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced' &
python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced'

python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Lz.fits -RIN 250 -ROUT 2000 -ncores 7 &
python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Lz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7 &
python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Lz.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced' &
python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Lz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced' &
python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Lz.fits -RIN 250 -ROUT 2000 -ncores 7 &
python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Lz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7 &
python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Lz.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced' &
python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Lz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced'
