#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv

python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_all_q4.fits -RIN 350 -ROUT 1500 -ncores 5 -nit 250 -ang reduced &
python fit_allprofiles_centred_sequential_2h_2q_Ein.py -folder ~/MICE/HS-lensing/profiles/ -file profile_all_q4.fits -RIN 350 -ROUT 5000 -ncores 5 -nit 250 -ang reduced &
python fit_allprofiles_centred_sequential_2h_2q_woc.py -folder ~/MICE/HS-lensing/profiles/ -file profile_all_q4.fits -RIN 350 -ROUT 5000 -ncores 5 -nit 250 -ang reduced &
python fit_allprofiles_centred_sequential_2h_2q.py     -folder ~/MICE/HS-lensing/profiles/ -file profile_all_q4.fits -RIN 350 -ROUT 5000 -ncores 5 -nit 250 -ang reduced &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_all_q5.fits -RIN 350 -ROUT 1500 -ncores 5 -nit 250 -ang reduced &
python fit_allprofiles_centred_sequential_2h_2q_Ein.py -folder ~/MICE/HS-lensing/profiles/ -file profile_all_q5.fits -RIN 350 -ROUT 5000 -ncores 5 -nit 250 -ang reduced &
python fit_allprofiles_centred_sequential_2h_2q_woc.py -folder ~/MICE/HS-lensing/profiles/ -file profile_all_q5.fits -RIN 350 -ROUT 5000 -ncores 5 -nit 250 -ang reduced &
python fit_allprofiles_centred_sequential_2h_2q.py     -folder ~/MICE/HS-lensing/profiles/ -file profile_all_q5.fits -RIN 350 -ROUT 5000 -ncores 5 -nit 250 -ang reduced &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_all_q6.fits -RIN 350 -ROUT 1500 -ncores 5 -nit 250 -ang reduced &
python fit_allprofiles_centred_sequential_2h_2q_Ein.py -folder ~/MICE/HS-lensing/profiles/ -file profile_all_q6.fits -RIN 350 -ROUT 5000 -ncores 5 -nit 250 -ang reduced &
python fit_allprofiles_centred_sequential_2h_2q_woc.py -folder ~/MICE/HS-lensing/profiles/ -file profile_all_q6.fits -RIN 350 -ROUT 5000 -ncores 5 -nit 250 -ang reduced &
python fit_allprofiles_centred_sequential_2h_2q.py     -folder ~/MICE/HS-lensing/profiles/ -file profile_all_q6.fits -RIN 350 -ROUT 5000 -ncores 5 -nit 250 -ang reduced &
wait

##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz.fits -RIN 250 -ROUT 2000 -ncores 7 &
##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7 &
##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced' &
##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced' &
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz.fits -RIN 250 -ROUT 2000 -ncores 7 &
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7 &
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced' &
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced'
