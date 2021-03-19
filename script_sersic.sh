#!/bin/bash
python -u forGroup_qprofile.py -sample '136_145' -lM_min 13.6 -lM_max 14.5 -z_max 0.12 -ncores 30
python -u forGroup_qprofile_individual.py -sample '136_145' -lM_min 13.6 -lM_max 14.5 -z_max 0.12 -ncores 30
python -u forGroup_map.py -sample '136_145' -lM_min 13.9 -lM_max 13.6 -z_max 0.12 -ncores 30
