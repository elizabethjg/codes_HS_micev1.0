#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv

python fit_allprofiles_centred_sequential_2h_2q.py -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Lz_qcut.fits -RIN 350 -ROUT 5000 -ncores 5 -nit 250 -qext '_des' &
python fit_allprofiles_centred_sequential_2h_2q.py -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Mz_qcut.fits -RIN 350 -ROUT 5000 -ncores 5 -nit 250 -qext '_des' &
python fit_allprofiles_centred_sequential_2h_2q.py -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Hz_qcut.fits -RIN 450 -ROUT 5000 -ncores 5 -nit 250 -qext '_des' &
python fit_allprofiles_centred_sequential_2h_2q.py -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Lz_qcut.fits -RIN 350 -ROUT 5000 -ncores 5 -nit 250 -qext '_des' &
python fit_allprofiles_centred_sequential_2h_2q.py -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Mz_qcut.fits -RIN 350 -ROUT 5000 -ncores 5 -nit 250 -qext '_des' &
python fit_allprofiles_centred_sequential_2h_2q.py -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz_qcut.fits -RIN 400 -ROUT 5000 -ncores 5 -nit 250 -qext '_des' &
wait
