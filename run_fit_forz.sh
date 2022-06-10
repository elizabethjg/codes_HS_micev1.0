#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv

python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_all.fits      -RIN 350 -ROUT 1500 -ncores 6 -nit 250 &
python fit_allprofiles_centred_sequential_2h_2q_Ein.py -folder ~/MICE/HS-lensing/profiles/ -file profile_all.fits      -RIN 350 -ROUT 5000 -ncores 6 -nit 250 &
python fit_allprofiles_centred_sequential_2h_2q_woc.py -folder ~/MICE/HS-lensing/profiles/ -file profile_all.fits      -RIN 350 -ROUT 5000 -ncores 6 -nit 250 &
python fit_allprofiles_centred_sequential_2h_2q.py     -folder ~/MICE/HS-lensing/profiles/ -file profile_all.fits      -RIN 350 -ROUT 5000 -ncores 6 -nit 250 &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_all_qmed.fits -RIN 350 -ROUT 1500 -ncores 6 -nit 250 &
python fit_allprofiles_centred_sequential_2h_2q_Ein.py -folder ~/MICE/HS-lensing/profiles/ -file profile_all_qmed.fits -RIN 350 -ROUT 5000 -ncores 6 -nit 250 &
python fit_allprofiles_centred_sequential_2h_2q_woc.py -folder ~/MICE/HS-lensing/profiles/ -file profile_all_qmed.fits -RIN 350 -ROUT 5000 -ncores 6 -nit 250 &
python fit_allprofiles_centred_sequential_2h_2q.py     -folder ~/MICE/HS-lensing/profiles/ -file profile_all_qmed.fits -RIN 350 -ROUT 5000 -ncores 6 -nit 250 &
wait
