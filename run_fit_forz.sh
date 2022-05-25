#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv

python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz.fits -RIN 400 -ROUT 1000 -ncores 7 -nit 250 -qext _mis10 &
python fit_allprofiles_centred_sequential_2h_2q_Ein.py -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz.fits -RIN 400 -ROUT 5000 -ncores 7 -nit 250 -qext _mis10 &
python fit_allprofiles_centred_sequential_2h_2q_woc.py -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz.fits -RIN 400 -ROUT 5000 -ncores 7 -nit 250 -qext _mis10 &
python fit_allprofiles_centred_sequential_2h_2q.py     -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz.fits -RIN 400 -ROUT 5000 -ncores 7 -nit 250 -qext _mis10 &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz.fits -RIN 400 -ROUT 1000 -ncores 7 -nit 250 -qext _mis20 &
python fit_allprofiles_centred_sequential_2h_2q_Ein.py -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz.fits -RIN 400 -ROUT 5000 -ncores 7 -nit 250 -qext _mis20 &
python fit_allprofiles_centred_sequential_2h_2q_woc.py -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz.fits -RIN 400 -ROUT 5000 -ncores 7 -nit 250 -qext _mis20 &
python fit_allprofiles_centred_sequential_2h_2q.py     -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz.fits -RIN 400 -ROUT 5000 -ncores 7 -nit 250 -qext _mis20 &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz.fits -RIN 400 -ROUT 1000 -ncores 7 -nit 250 -qext _mis30 &
python fit_allprofiles_centred_sequential_2h_2q_Ein.py -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz.fits -RIN 400 -ROUT 5000 -ncores 7 -nit 250 -qext _mis30 &
python fit_allprofiles_centred_sequential_2h_2q_woc.py -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz.fits -RIN 400 -ROUT 5000 -ncores 7 -nit 250 -qext _mis30 &
python fit_allprofiles_centred_sequential_2h_2q.py     -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz.fits -RIN 400 -ROUT 5000 -ncores 7 -nit 250 -qext _mis30 &
wait
