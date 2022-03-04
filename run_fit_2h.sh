#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv

python fit_allprofiles_centred_sequential_2h.py -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz_relaxed.fits -RIN 250 -ROUT 5000  -ncores 10 -nit 250 &
python fit_allprofiles_centred_sequential_2h.py -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Hz_relaxed.fits -RIN 250 -ROUT 5000  -ncores 10 -nit 250 &
python fit_allprofiles_centred_sequential_2h.py -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz_relaxed.fits -RIN 250 -ROUT 5000  -ncores 10 -nit 250 -ang reduced &
python fit_allprofiles_centred_sequential_2h.py -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Hz_relaxed.fits -RIN 250 -ROUT 5000  -ncores 10 -nit 250 -ang reduced &
python fit_allprofiles_centred_sequential_2h.py -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz.fits -RIN 250 -ROUT 5000  -ncores 10 -nit 250 &
python fit_allprofiles_centred_sequential_2h.py -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Hz.fits -RIN 250 -ROUT 5000  -ncores 10 -nit 250 &
python fit_allprofiles_centred_sequential_2h.py -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz.fits -RIN 250 -ROUT 5000  -ncores 10 -nit 250 -ang reduced &
python fit_allprofiles_centred_sequential_2h.py -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Hz.fits -RIN 250 -ROUT 5000  -ncores 10 -nit 250 -ang reduced &
wait
