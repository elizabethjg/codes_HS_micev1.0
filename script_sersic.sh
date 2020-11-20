#!/bin/bash
srun python -u forGroup_qprofile.py -sample 'bin_146' -lM_min 14.6 -lM_max 15.0 -h_cosmo 0.7 -ncores 35
srun python -u forGroup_qprofile.py -sample 'bin_144' -lM_min 14.4 -lM_max 14.6 -h_cosmo 0.7 -ncores 35
srun python -u forGroup_qprofile.py -sample 'bin_140' -lM_min 14.0 -lM_max 14.2 -h_cosmo 0.7 -ncores 35
srun python -u forGroup_qprofile.py -sample 'bin_136' -lM_min 13.6 -lM_max 13.8 -h_cosmo 0.7 -ncores 35
srun python -u forGroup_qprofile.py -sample 'bin_132' -lM_min 13.2 -lM_max 13.4 -h_cosmo 0.7 -ncores 35
