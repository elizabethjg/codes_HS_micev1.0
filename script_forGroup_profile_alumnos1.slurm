#!/bin/bash
. /etc/profile
source /home/elizabeth/.bashrc

conda activate myenv

python -u forGroup_qprofile.py -sample 'all'      -lM_min 13.5 -lM_max 14.5 -z_min 0.1 -z_max 0.3 -ncores 32 -relax True
python -u forGroup_qprofile.py -sample 'all_qmin' -lM_min 13.5 -lM_max 14.5 -z_min 0.1 -z_max 0.3 -ncores 32 -q_min 0.0 -q_max 0.5 -relax True
