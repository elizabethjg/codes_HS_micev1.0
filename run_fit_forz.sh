#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv

python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz.fits -RIN 400 -ROUT 1000 -ncores 5 -nit 250 &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Hz.fits -RIN 450 -ROUT 2000 -ncores 5 -nit 250 &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Mz.fits -RIN 350 -ROUT 1000 -ncores 10 -nit 250 &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Lz.fits -RIN 350 -ROUT 1000  -ncores 7 -nit 250 &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Lz.fits -RIN 350 -ROUT 2000  -ncores 7 -nit 250 &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Mz.fits -RIN 350 -ROUT 2000  -ncores 7 -nit 250 &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz.fits -RIN 400 -ROUT 1000 -ncores 5 -nit 250 -qext _miscen &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Hz.fits -RIN 450 -ROUT 2000 -ncores 5 -nit 250 -qext _miscen &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Mz.fits -RIN 350 -ROUT 1000 -ncores 10 -nit 250 -qext _miscen &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Lz.fits -RIN 350 -ROUT 1000  -ncores 7 -nit 250 -qext _miscen &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Lz.fits -RIN 350 -ROUT 2000  -ncores 7 -nit 250 -qext _miscen &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Mz.fits -RIN 350 -ROUT 2000  -ncores 7 -nit 250 -qext _miscen &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz.fits -RIN 400 -ROUT 1000 -ncores 5 -nit 250 -qext _mis20 &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Hz.fits -RIN 450 -ROUT 2000 -ncores 5 -nit 250 -qext _mis20 &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Mz.fits -RIN 350 -ROUT 1000 -ncores 10 -nit 250 -qext _mis20 &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Lz.fits -RIN 350 -ROUT 1000  -ncores 7 -nit 250 -qext _mis20 &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Lz.fits -RIN 350 -ROUT 2000  -ncores 7 -nit 250 -qext _mis20 &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Mz.fits -RIN 350 -ROUT 2000  -ncores 7 -nit 250 -qext _mis20 &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Hz.fits -RIN 400 -ROUT 1000 -ncores 5 -nit 250 -qext _mis20_miscen &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Hz.fits -RIN 450 -ROUT 2000 -ncores 5 -nit 250 -qext _mis20_miscen &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Mz.fits -RIN 350 -ROUT 1000 -ncores 10 -nit 250 -qext _mis20_miscen &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_LM_Lz.fits -RIN 350 -ROUT 1000  -ncores 7 -nit 250 -qext _mis20_miscen &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Lz.fits -RIN 350 -ROUT 2000  -ncores 7 -nit 250 -qext _mis20_miscen &
python fit_allprofiles_centred_sequential.py           -folder ~/MICE/HS-lensing/profiles/ -file profile_HM_Mz.fits -RIN 350 -ROUT 2000  -ncores 7 -nit 250 -qext _mis20_miscen &
wait

##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz.fits -RIN 250 -ROUT 2000 -ncores 7 &
##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7 &
##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced' &
##python fit_allprofiles_centred_withoutc200.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced' &
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz.fits -RIN 250 -ROUT 2000 -ncores 7 &
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7 &
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced' &
##python fit_allprofiles_centred.py -folder '/home/elizabeth/MICE/HS-lensing/profiles/' -file profile_$1_Hz_relaxed.fits -RIN 250 -ROUT 2000 -ncores 7 -ang 'reduced'
