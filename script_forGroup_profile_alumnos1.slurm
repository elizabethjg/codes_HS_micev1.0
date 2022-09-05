#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv

##python -u forGroup_qprofile.py -sample 'HM_Lz_qcut_n6'  -lM_min 14.0 -lM_max 14.5 -z_min 0.1 -z_max 0.2 -ncores 32 -q_min 0.5 -nbins 30 -addnoise True -nback 6.3
##python -u forGroup_qprofile.py -sample 'LM_Lz_qcut_n6'  -lM_min 13.5 -lM_max 14.0 -z_min 0.1 -z_max 0.2 -ncores 32 -q_min 0.5 -nbins 30 -addnoise True -nback 6.3
python -u forGroup_qprofile.py -sample 'LM_Mz_qcut_n6'  -lM_min 13.5 -lM_max 14.0 -z_min 0.2 -z_max 0.3 -ncores 60 -q_min 0.5 -nbins 30 -addnoise True -nback 6.3
##python -u forGroup_qprofile.py -sample 'HM_Lz_qcut_n30' -lM_min 14.0 -lM_max 14.5 -z_min 0.1 -z_max 0.2 -ncores 32 -q_min 0.5 -nbins 40 -addnoise True -nback 30.
python -u forGroup_qprofile.py -sample 'LM_Lz_qcut_n30' -lM_min 13.5 -lM_max 14.0 -z_min 0.1 -z_max 0.2 -ncores 60 -q_min 0.5 -nbins 40 -addnoise True -nback 30.
python -u forGroup_qprofile.py -sample 'LM_Mz_qcut_n30' -lM_min 13.5 -lM_max 14.0 -z_min 0.2 -z_max 0.3 -ncores 60 -q_min 0.5 -nbins 40 -addnoise True -nback 30.
##srun python -u forGroup_qprofile.py -sample 'HM_Mz_qcut' -lM_min 14.0 -lM_max 14.5 -z_min 0.2 -z_max 0.3 -ncores 32 -q_min 0.5
##srun python -u forGroup_qprofile.py -sample 'HM_Hz_qcut' -lM_min 14.0 -lM_max 14.5 -z_min 0.3 -z_max 0.4 -ncores 32 -q_min 0.5
##srun python -u forGroup_qprofile.py -sample 'LM_Hz_qcut' -lM_min 13.5 -lM_max 14.0 -z_min 0.3 -z_max 0.4 -ncores 32 -q_min 0.5                                                                                                                          

