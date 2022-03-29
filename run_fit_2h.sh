#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv

python fit_allprofiles_centred_sequential_2h_2q.py     -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_HHz.fits -RIN 350 -ROUT 5000  -ncores 6 -nit 250 -ang reduced &
python fit_allprofiles_centred_sequential_2h_2q.py     -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_HHz.fits -RIN 350 -ROUT 5000  -ncores 6 -nit 250 -ang reduced &
python fit_allprofiles_centred_sequential_2h_2q.py     -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Mz.fits -RIN 350 -ROUT 5000  -ncores 6 -nit 250 -ang reduced &
python fit_allprofiles_centred_sequential_2h_2q.py     -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Lz.fits -RIN 350 -ROUT 5000  -ncores 6 -nit 250 -ang reduced &
python fit_allprofiles_centred_sequential_2h_2q.py     -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Lz.fits -RIN 350 -ROUT 5000  -ncores 6 -nit 250 -ang reduced &
python fit_allprofiles_centred_sequential_2h_2q.py     -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Mz.fits -RIN 350 -ROUT 5000  -ncores 6 -nit 250 -ang reduced &
wait

