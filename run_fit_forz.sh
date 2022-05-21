#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv

python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz.fits -RIN 400 -ROUT 1000 -ncores 7 -nit 250 -qext _mis20 &
python fit_allprofiles_centred_sequential_2h_2q_Ein.py -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz.fits -RIN 400 -ROUT 5000 -ncores 7 -nit 250 -qext _mis20 &
python fit_allprofiles_centred_sequential_2h_2q_woc.py -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz.fits -RIN 400 -ROUT 5000 -ncores 7 -nit 250 -qext _mis20 &
python fit_allprofiles_centred_sequential_2h_2q.py     -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz.fits -RIN 400 -ROUT 5000 -ncores 7 -nit 250 -qext _mis20 &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Hz.fits -RIN 450 -ROUT 2000 -ncores 7 -nit 250 -qext _mis20 &
python fit_allprofiles_centred_sequential_2h_2q_Ein.py -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Hz.fits -RIN 450 -ROUT 5000 -ncores 7 -nit 250 -qext _mis20 &
python fit_allprofiles_centred_sequential_2h_2q_woc.py -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Hz.fits -RIN 450 -ROUT 5000 -ncores 7 -nit 250 -qext _mis20 &
python fit_allprofiles_centred_sequential_2h_2q.py     -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Hz.fits -RIN 450 -ROUT 5000 -ncores 7 -nit 250 -qext _mis20 &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Mz.fits -RIN 350 -ROUT 1000 -ncores 7 -nit 250 -qext _mis20 &
python fit_allprofiles_centred_sequential_2h_2q_Ein.py -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Mz.fits -RIN 350 -ROUT 5000 -ncores 7 -nit 250 -qext _mis20 &
python fit_allprofiles_centred_sequential_2h_2q_woc.py -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Mz.fits -RIN 350 -ROUT 5000 -ncores 7 -nit 250 -qext _mis20 &
python fit_allprofiles_centred_sequential_2h_2q.py     -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Mz.fits -RIN 350 -ROUT 5000 -ncores 7 -nit 250 -qext _mis20 &
wait
