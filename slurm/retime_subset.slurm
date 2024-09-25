#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=carbon-kp
#SBATCH --partition=carbon-kp
#SBATCH -o ./logs/csl_retime_%A.out

source activate noaa_csl
cd /uufs/chpc.utah.edu/common/home/u0890904/NOAA_CSL/noaa_csl

python -u retime_subset.py

