#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv

python fit_allprofiles_centred_sequential_2h_2q.py -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz_relaxed_miscen.fits -RIN 250 -ROUT 5000  -ncores 10 -nit 250 &
python fit_allprofiles_centred_sequential_2h_2q.py -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Hz_relaxed_miscen.fits -RIN 250 -ROUT 5000  -ncores 10 -nit 250 &
python fit_allprofiles_centred_sequential_2h_2q.py -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz_relaxed_miscen.fits -RIN 250 -ROUT 5000  -ncores 10 -nit 250 -ang reduced &
python fit_allprofiles_centred_sequential_2h_2q.py -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Hz_relaxed_miscen.fits -RIN 250 -ROUT 5000  -ncores 10 -nit 250 -ang reduced &
python fit_allprofiles_centred_sequential_2h_2q.py -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Lz_relaxed_miscen.fits -RIN 250 -ROUT 5000  -ncores 10 -nit 250 &
python fit_allprofiles_centred_sequential_2h_2q.py -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Lz_relaxed_miscen.fits -RIN 250 -ROUT 5000  -ncores 10 -nit 250 &
python fit_allprofiles_centred_sequential_2h_2q.py -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Lz_relaxed_miscen.fits -RIN 250 -ROUT 5000  -ncores 10 -nit 250 -ang reduced &
python fit_allprofiles_centred_sequential_2h_2q.py -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Lz_relaxed_miscen.fits -RIN 250 -ROUT 5000  -ncores 10 -nit 250 -ang reduced &
wait

##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz.fits -RIN 250 -ROUT 2000 -ncores 7 &
##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7 &
##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced' &
##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced' &
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz.fits -RIN 250 -ROUT 2000 -ncores 7 &
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7 &
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced' &
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced'
