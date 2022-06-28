#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv

python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz_qcut.fits -RIN 400 -ROUT 1000 -ncores 7 -nit 250  &
python fit_allprofiles_centred_sequential_2h_2q_Ein.py -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz_qcut.fits -RIN 400 -ROUT 5000 -ncores 7 -nit 250  &
python fit_allprofiles_centred_sequential_2h_2q_woc.py -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz_qcut.fits -RIN 400 -ROUT 5000 -ncores 7 -nit 250  &
python fit_allprofiles_centred_sequential_2h_2q.py     -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz_qcut.fits -RIN 400 -ROUT 5000 -ncores 7 -nit 250  &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Hz_qcut.fits -RIN 450 -ROUT 2000 -ncores 7 -nit 250  &
python fit_allprofiles_centred_sequential_2h_2q_Ein.py -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Hz_qcut.fits -RIN 450 -ROUT 5000 -ncores 7 -nit 250  &
python fit_allprofiles_centred_sequential_2h_2q_woc.py -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Hz_qcut.fits -RIN 450 -ROUT 5000 -ncores 7 -nit 250  &
python fit_allprofiles_centred_sequential_2h_2q.py     -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Hz_qcut.fits -RIN 450 -ROUT 5000 -ncores 7 -nit 250  &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Mz_qcut.fits -RIN 350 -ROUT 2000 -ncores 7 -nit 250  &
python fit_allprofiles_centred_sequential_2h_2q_Ein.py -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Mz_qcut.fits -RIN 350 -ROUT 5000 -ncores 7 -nit 250  &
python fit_allprofiles_centred_sequential_2h_2q_woc.py -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Mz_qcut.fits -RIN 350 -ROUT 5000 -ncores 7 -nit 250  &
python fit_allprofiles_centred_sequential_2h_2q.py     -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Mz_qcut.fits -RIN 350 -ROUT 5000 -ncores 7 -nit 250  &

##python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Lz_relaxed.fits -RIN 350 -ROUT 1000 -ncores 7 -nit 250 -qext _mis20_miscen &
##python fit_allprofiles_centred_sequential_2h_2q_Ein.py -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Lz_relaxed.fits -RIN 350 -ROUT 5000 -ncores 7 -nit 250 -qext _mis20_miscen &
##python fit_allprofiles_centred_sequential_2h_2q_woc.py -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Lz_relaxed.fits -RIN 350 -ROUT 5000 -ncores 7 -nit 250 -qext _mis20_miscen &
##python fit_allprofiles_centred_sequential_2h_2q.py     -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Lz_relaxed.fits -RIN 350 -ROUT 5000 -ncores 7 -nit 250 -qext _mis20_miscen &
##python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Lz_relaxed.fits -RIN 350 -ROUT 2000 -ncores 7 -nit 250 -qext _mis20_miscen &
##python fit_allprofiles_centred_sequential_2h_2q_Ein.py -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Lz_relaxed.fits -RIN 350 -ROUT 5000 -ncores 7 -nit 250 -qext _mis20_miscen &
##python fit_allprofiles_centred_sequential_2h_2q_woc.py -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Lz_relaxed.fits -RIN 350 -ROUT 5000 -ncores 7 -nit 250 -qext _mis20_miscen &
##python fit_allprofiles_centred_sequential_2h_2q.py     -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Lz_relaxed.fits -RIN 350 -ROUT 5000 -ncores 7 -nit 250 -qext _mis20_miscen &
##python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Mz_relaxed.fits -RIN 350 -ROUT 1000 -ncores 7 -nit 250 -qext _mis20_miscen &
##python fit_allprofiles_centred_sequential_2h_2q_Ein.py -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Mz_relaxed.fits -RIN 350 -ROUT 5000 -ncores 7 -nit 250 -qext _mis20_miscen &
##python fit_allprofiles_centred_sequential_2h_2q_woc.py -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Mz_relaxed.fits -RIN 350 -ROUT 5000 -ncores 7 -nit 250 -qext _mis20_miscen &
##python fit_allprofiles_centred_sequential_2h_2q.py     -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Mz_relaxed.fits -RIN 350 -ROUT 5000 -ncores 7 -nit 250 -qext _mis20_miscen &

##python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_all_q1_relaxed.fits -RIN 350 -ROUT 1500 -ncores 5 -nit 250 -ang reduced &
##python fit_allprofiles_centred_sequential_2h_2q_Ein.py -folder ~/MICE/HS-lensing/profiles/ -file profile_all_q1_relaxed.fits -RIN 350 -ROUT 5000 -ncores 5 -nit 250 -ang reduced &
##python fit_allprofiles_centred_sequential_2h_2q_woc.py -folder ~/MICE/HS-lensing/profiles/ -file profile_all_q1_relaxed.fits -RIN 350 -ROUT 5000 -ncores 5 -nit 250 -ang reduced &
##python fit_allprofiles_centred_sequential_2h_2q.py     -folder ~/MICE/HS-lensing/profiles/ -file profile_all_q1_relaxed.fits -RIN 350 -ROUT 5000 -ncores 5 -nit 250 -ang reduced &
##python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_all_q2_relaxed.fits -RIN 350 -ROUT 1500 -ncores 5 -nit 250 -ang reduced &
##python fit_allprofiles_centred_sequential_2h_2q_Ein.py -folder ~/MICE/HS-lensing/profiles/ -file profile_all_q2_relaxed.fits -RIN 350 -ROUT 5000 -ncores 5 -nit 250 -ang reduced &
##python fit_allprofiles_centred_sequential_2h_2q_woc.py -folder ~/MICE/HS-lensing/profiles/ -file profile_all_q2_relaxed.fits -RIN 350 -ROUT 5000 -ncores 5 -nit 250 -ang reduced &
##python fit_allprofiles_centred_sequential_2h_2q.py     -folder ~/MICE/HS-lensing/profiles/ -file profile_all_q2_relaxed.fits -RIN 350 -ROUT 5000 -ncores 5 -nit 250 -ang reduced &
##python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_all_q3_relaxed.fits -RIN 350 -ROUT 1500 -ncores 5 -nit 250 -ang reduced &
##python fit_allprofiles_centred_sequential_2h_2q_Ein.py -folder ~/MICE/HS-lensing/profiles/ -file profile_all_q3_relaxed.fits -RIN 350 -ROUT 5000 -ncores 5 -nit 250 -ang reduced &
##python fit_allprofiles_centred_sequential_2h_2q_woc.py -folder ~/MICE/HS-lensing/profiles/ -file profile_all_q3_relaxed.fits -RIN 350 -ROUT 5000 -ncores 5 -nit 250 -ang reduced &
##python fit_allprofiles_centred_sequential_2h_2q.py     -folder ~/MICE/HS-lensing/profiles/ -file profile_all_q3_relaxed.fits -RIN 350 -ROUT 5000 -ncores 5 -nit 250 -ang reduced &
wait
